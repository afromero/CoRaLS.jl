#!/usr/bin/env julia

"""Record matched direct and reflected signals for deterministic facet tilts.

One valid event geometry (surface point, payload position, CR axis, and Xmax)
is drawn per event ID.  Every requested facet tilt then evaluates that exact
same candidate with the same facet azimuth.  This makes event-by-event changes
relative to the zero-tilt facet meaningful: an event can be retained, lost, or
gained at a larger tilt.

By default the facet azimuth is sampled isotropically once per event and is
held fixed across all tilt magnitudes.  ``--facet-aspect`` instead fixes that
azimuth relative to the payload-projected tangent direction: 0 degrees points
toward the payload, 180 degrees points away, and +90 degrees is sideways.
"""

import Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "..")))
include(joinpath(@__DIR__, "facet_diagnostic_common.jl"))
using Unitful: uconvert, ustrip

const DEFAULT_TRIALS = 20_000
const DEFAULT_ENERGY = 10.0
const DEFAULT_REGION = "polar:south,-87,1.0"
const DEFAULT_TILTS = (0.0, 5.0, 10.0, 20.0, 30.0, 45.0)
const DEFAULT_THETA0S = (-90.0, -48.0)
const DEFAULT_SEED = 20261109
const DEFAULT_OUTPUT = joinpath(@__DIR__, "..", "data", "terrain_diagnostics", "fixed_facet_slope_step1.csv")
const DEFAULT_WRITE_BATCH_EVENTS = 5_000

function usage()
    println("""
Usage: julia --project=. scripts/record_fixed_facet_slope_step1.jl [options]

  --ntrials N          shared valid event geometries (default: $DEFAULT_TRIALS)
  --energy-eev E       fixed CR energy (default: $DEFAULT_ENERGY)
  --region CONFIG      CoRaLS region (default: $DEFAULT_REGION)
  --facet-tilts A,B    deterministic slope magnitudes in degrees (default: $(join(DEFAULT_TILTS, ',')))
  --theta0s A,B        LPDA boresights in degrees (default: $(join(DEFAULT_THETA0S, ',')))
  --facet-aspect D     optional fixed aspect relative to the payload; omit for one
                        shared isotropic azimuth draw per event
  --seed N             random seed (default: $DEFAULT_SEED)
  --write-batch-events N
                       number of paired candidate events buffered before appending
                       CSV rows (default: $DEFAULT_WRITE_BATCH_EVENTS)
  --output PATH        output CSV path
""")
end

function options(arguments)
    cfg = Dict{Symbol, Any}(
        :ntrials=>DEFAULT_TRIALS, :energy=>DEFAULT_ENERGY, :region=>DEFAULT_REGION,
        :tilts=>collect(DEFAULT_TILTS), :theta0s=>collect(DEFAULT_THETA0S),
        :facet_aspect=>nothing, :seed=>DEFAULT_SEED, :output=>DEFAULT_OUTPUT,
        :write_batch_events=>DEFAULT_WRITE_BATCH_EVENTS,
    )
    i = 1
    while i <= length(arguments)
        arg = arguments[i]
        arg == "--help" && (usage(); exit(0))
        arg in ("--ntrials", "--energy-eev", "--region", "--facet-tilts", "--theta0s",
            "--facet-aspect", "--seed", "--write-batch-events", "--output") || error("unknown option: $arg")
        i == length(arguments) && error("$arg requires a value")
        value = arguments[i + 1]
        if arg == "--ntrials"; cfg[:ntrials] = parse(Int, value)
        elseif arg == "--energy-eev"; cfg[:energy] = parse(Float64, value)
        elseif arg == "--region"; cfg[:region] = value
        elseif arg == "--facet-tilts"; cfg[:tilts] = parse.(Float64, split(value, ','))
        elseif arg == "--theta0s"; cfg[:theta0s] = parse.(Float64, split(value, ','))
        elseif arg == "--facet-aspect"; cfg[:facet_aspect] = parse(Float64, value)
        elseif arg == "--seed"; cfg[:seed] = parse(Int, value)
        elseif arg == "--write-batch-events"; cfg[:write_batch_events] = parse(Int, value)
        else; cfg[:output] = value
        end
        i += 2
    end
    cfg[:ntrials] > 0 || error("--ntrials must be positive")
    cfg[:write_batch_events] > 0 || error("--write-batch-events must be positive")
    cfg[:energy] > 0 || error("--energy-eev must be positive")
    !isempty(cfg[:tilts]) && all(0.0 .<= cfg[:tilts] .<= 90.0) ||
        error("--facet-tilts must lie in [0,90]")
    !isempty(cfg[:theta0s]) && all(-90.0 .<= cfg[:theta0s] .<= 0.0) ||
        error("--theta0s must lie in [-90,0]")
    create_region(cfg[:region])
    return cfg
end

"""Append buffered rows while retaining event-major CSV ordering."""
function flush_rows!(rows, output, has_written)
    isempty(rows) && return has_written
    CSV.write(output, DataFrame(rows); append=has_written, writeheader=!has_written)
    empty!(rows)
    return true
end

payload_field_uV_per_m(event) = ustrip(norm(event.pol) * event.dν * sum(event.Ef) / (μV / m))

longitude_deg(phi) = mod(rad2deg(phi) + 180.0, 360.0) - 180.0

function payload_geometry(candidate)
    surface_theta, surface_phi, _ = CoRaLS.cartesian_to_spherical(candidate.surface)
    payload_theta, payload_phi, _ = CoRaLS.cartesian_to_spherical(candidate.antenna)
    _, _, payload_altitude = CoRaLS.cartesian_to_latlonalt(candidate.antenna)
    # A local north/east basis defines a payload bearing in the event sky.
    north_reference = typeof(candidate.flat_normal)(0.0, 0.0, 1.0)
    north = north_reference - (north_reference ⋅ candidate.flat_normal) * candidate.flat_normal
    if norm(north) < 1e-12
        north_reference = typeof(candidate.flat_normal)(1.0, 0.0, 0.0)
        north = north_reference - (north_reference ⋅ candidate.flat_normal) * candidate.flat_normal
    end
    north /= norm(north)
    east = north × candidate.flat_normal
    east /= norm(east)
    payload_zenith = rad2deg(acos(clamp(candidate.flat_normal ⋅ candidate.view, -1.0, 1.0)))
    payload_bearing = mod(rad2deg(atan(candidate.view ⋅ east, candidate.view ⋅ north)), 360.0)
    return (
        surface_lat_deg=90.0 - rad2deg(surface_theta),
        surface_lon_deg=longitude_deg(surface_phi),
        payload_lat_deg=90.0 - rad2deg(payload_theta),
        payload_lon_deg=longitude_deg(payload_phi),
        payload_altitude_km=ustrip(uconvert(km, payload_altitude)),
        payload_x_km=ustrip(uconvert(km, candidate.antenna[1])),
        payload_y_km=ustrip(uconvert(km, candidate.antenna[2])),
        payload_z_km=ustrip(uconvert(km, candidate.antenna[3])),
        payload_zenith_at_surface_deg=payload_zenith,
        payload_elevation_at_surface_deg=90.0 - payload_zenith,
        payload_bearing_at_surface_deg=payload_bearing,
    )
end

function route_row(event_id, route, tilt, aspect, theta0, candidate, normal, signal, response, cfg, geometry)
    status = failure_name(signal)
    common = (
        event_id=event_id, route=route, facet_tilt_deg=tilt, facet_aspect_deg=aspect,
        lpda_theta0_deg=theta0, region=cfg[:region], energy_EeV=cfg[:energy],
        raw_exit_angle_theta_r_deg=exit_angle_deg(candidate, normal),
        backfacing=CoRaLS.is_backfacing_facet(normal, candidate.view),
        status=status, geometry...,
    )
    if signal isa CoRaLS.TrialFailed
        return (; common..., triggered=false, theta_i_deg=missing, theta_emit_deg=missing,
            theta_ice_deg=missing, psi_deg=missing, payload_elevation_signal_deg=missing,
            field_uV_per_m=missing, lpda_voltage_uV=missing, lpda_threshold_uV=missing,
            lpda_snr_proxy=missing)
    end
    voltage = ustrip(response.voltage / μV)
    threshold = ustrip(response.threshold / μV)
    return (; common..., triggered=response.voltage > response.threshold,
        theta_i_deg=signal.θ_i, theta_emit_deg=signal.θ_emit,
        theta_ice_deg=route == "reflected" ? signal.θ_ice : missing,
        psi_deg=signal.ψ, payload_elevation_signal_deg=signal.θ_el,
        field_uV_per_m=payload_field_uV_per_m(signal), lpda_voltage_uV=voltage,
        lpda_threshold_uV=threshold, lpda_snr_proxy=voltage / threshold)
end

function main()
    cfg = options(ARGS)
    Random.seed!(cfg[:seed])
    region = create_region(cfg[:region])
    spacecraft = CircularOrbit(FACET_DEFAULT_ALTITUDE)
    energy = cfg[:energy] * EeV
    material = diagnostic_material_options()
    triggers = Dict(theta0 => make_lpda(theta0; return_details=true) for theta0 in cfg[:theta0s])
    println("Fixed-facet paired diagnostic: $(cfg[:ntrials]) shared candidates; tilts=$(join(cfg[:tilts], ','))°; routes=direct,R")
    output = abspath(cfg[:output])
    ispath(output) && error("output already exists: $output")
    mkpath(dirname(output))

    # Use an independent RNG for facet aspects.  Candidate draws therefore
    # retain their old seed sequence while the recorder can process and flush
    # one candidate at a time instead of retaining every geometry in memory.
    aspect_rng = MersenneTwister(cfg[:seed] + 1)
    rows_per_event = 2 * length(cfg[:tilts]) * length(cfg[:theta0s])
    rows = NamedTuple[]
    sizehint!(rows, cfg[:write_batch_events] * rows_per_event)
    has_written = false

    for event_id in 1:cfg[:ntrials]
        candidate = draw_candidate(region, spacecraft, energy;
            ice_depth=material.ice_depth, densitymodel=material.densitymodel)
        aspect = isnothing(cfg[:facet_aspect]) ? rand(aspect_rng) * 360.0 : cfg[:facet_aspect]
        geometry = payload_geometry(candidate)
        for tilt in cfg[:tilts]
            normal = facet_normal(candidate.flat_normal, candidate.view, tilt, aspect)
            direct, reflected = evaluate_facet(candidate, energy, normal, true, material)
            for theta0 in cfg[:theta0s]
                direct_response = direct isa CoRaLS.TrialFailed ? nothing : triggers[theta0](direct)
                reflected_response = reflected isa CoRaLS.TrialFailed ? nothing : triggers[theta0](reflected)
                push!(rows, route_row(event_id, "direct", tilt, aspect, theta0,
                    candidate, normal, direct, direct_response, cfg, geometry))
                push!(rows, route_row(event_id, "reflected", tilt, aspect, theta0,
                    candidate, normal, reflected, reflected_response, cfg, geometry))
            end
        end
        if event_id % cfg[:write_batch_events] == 0
            has_written = flush_rows!(rows, output, has_written)
            println("  recorded $event_id / $(cfg[:ntrials]) candidates")
        end
    end
    flush_rows!(rows, output, has_written)
    println("Wrote $output")
end

main()
