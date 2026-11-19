"""Shared helpers for controlled local-facet diagnostics.

These helpers intentionally condition on events that already have a visible
surface point and Xmax above the ice.  They therefore measure trigger and
escape-response fractions, *not* an absolute acceptance or event rate.
"""

using CoRaLS
using CSV
using DataFrames
using Interpolations
using LinearAlgebra
using Random
using Statistics
using Unitful: EeV, MHz, km, m, μV

const FACET_DEFAULT_ALTITUDE = 50km
const FACET_DEFAULT_ICE_DEPTH = 5m

"""Return a unit normal tilted from `flat_normal` toward the payload by `tilt_deg`.

`aspect_deg=0` is toward the payload projection; 180 degrees is away.  Other
aspects rotate around the spherical-surface normal.  This is deliberately a
controlled local boundary normal, not a self-consistent topographic surface.
"""
function facet_normal(flat_normal, view, tilt_deg, aspect_deg)
    toward = view - (view ⋅ flat_normal) * flat_normal
    if norm(toward) < 1e-12
        # The sub-payload point has no preferred tangent direction.  It is not
        # encountered for the near-horizon diagnostic events, but retain a
        # deterministic basis for completeness.
        reference = abs(flat_normal[3]) < 0.9 ?
            typeof(flat_normal)(0.0, 0.0, 1.0) : typeof(flat_normal)(1.0, 0.0, 0.0)
        toward = reference - (reference ⋅ flat_normal) * flat_normal
    end
    toward /= norm(toward)
    sideways = flat_normal × toward
    tangent = cosd(aspect_deg) * toward + sind(aspect_deg) * sideways
    normal = cosd(tilt_deg) * flat_normal + sind(tilt_deg) * tangent
    return normal / norm(normal)
end

"""Draw one event geometry usable by both direct and reflected scalar paths."""
function draw_candidate(region, spacecraft, energy; ice_depth=FACET_DEFAULT_ICE_DEPTH,
        densitymodel=LSB_Diviner_Density())
    while true
        # Sampling inside an AOI avoids spending almost all diagnostic draws
        # rejecting whole-Moon points for a 3-degree polar cap.
        surface = region isa CoRaLS.AOIRegion ?
            random_point_in_aoi(region.aoi) : random_point_on_sphere(Rmoon)
        CoRaLS.is_in_region(surface, region) || continue
        antenna = get_position(spacecraft)
        CoRaLS.is_visible(surface, antenna) || continue

        flat_normal = surface / norm(surface)
        axis = random_direction(flat_normal)
        xmax = CoRaLS.propagate_to_Xmax(surface, axis, energy, densitymodel)
        norm(xmax) <= Rmoon || continue
        depth = norm(surface) - norm(xmax)
        depth <= ice_depth || continue
        view = (antenna - surface) / norm(antenna - surface)
        return (; surface, antenna, axis, xmax, flat_normal, view)
    end
end

function draw_candidates(ntrials, region, spacecraft, energy; kwargs...)
    return [draw_candidate(region, spacecraft, energy; kwargs...) for _ in 1:ntrials]
end

function diagnostic_material_options(; ice_depth=FACET_DEFAULT_ICE_DEPTH)
    attenuation_table = CSV.read(joinpath(pkgdir(CoRaLS), "data", "allOnes.csv"), DataFrame)
    satten = linear_interpolation(attenuation_table[:, 1], attenuation_table[:, 2])
    indexmodel = LSB_DivinerIndex()
    return (
        ice_depth=ice_depth,
        ice_thickness=1m,
        Nice=1.265,
        Nbed=regolith_index(indexmodel, ice_depth),
        ν_min=150MHz,
        ν_max=800MHz,
        dν=10MHz,
        tand_mag=0.00025,
        tanδnorm=0.0005,
        slopemodel=NoSlope(),
        roughnessmodel=NoRoughness(),
        iceroughness=NoIceRoughness(),
        low_temp_corr_factor=0.9,
        satten=satten,
        indexmodel=indexmodel,
        densitymodel=LSB_Diviner_Density(),
        fieldmodel=ARW(),
        divergencemodel=MixedFieldDivergence(),
    )
end

"""Evaluate both paths for one pre-drawn event and one prescribed facet."""
function evaluate_facet(candidate, energy, normal, reject_backfacing, material)
    direct = CoRaLS.compute_direct(
        ScalarGeometry(), energy, candidate.surface, candidate.axis,
        candidate.antenna, candidate.xmax;
        material...,
        surface_normal_override=normal,
        reject_backfacing=reject_backfacing,
    )
    reflected = CoRaLS.compute_reflected(
        ScalarGeometry(), energy, candidate.surface, candidate.axis,
        candidate.antenna, candidate.xmax, material.ice_depth;
        material...,
        surface_normal_override=normal,
        reject_backfacing=reject_backfacing,
    )
    return direct, reflected
end

failure_name(result) = result isa CoRaLS.TrialFailed ? string(result) : "signal"
did_trigger(result, trigger) = !(result isa CoRaLS.TrialFailed) && trigger(result)

function exit_angle_deg(candidate, normal)
    return rad2deg(acos(clamp(candidate.view ⋅ normal, -1.0, 1.0)))
end

function make_lpda(theta0_deg; altitude=FACET_DEFAULT_ALTITUDE, kwargs...)
    return LPDA(Nant=8, Ntrig=3, θ0=theta0_deg, altitude=altitude,
        skyfrac=0.15; kwargs...)
end

"""Reconstruct the scalar reflected ray at its approximate in-medium source.

`compute_reflected` uses a mirror-path approximation and stores its selected
emission angle, but not the vector itself.  This reproduces that construction
so backward diagnostics can sample axes around the reflected radio direction.
"""
function reflected_radio_direction(candidate, reflected, ice_depth)
    depth = reflected.depth
    theta_emit = deg2rad(reflected.θ_emit)
    tangent_surface = candidate.view - (candidate.view ⋅ candidate.flat_normal) * candidate.flat_normal
    tangent_surface /= norm(tangent_surface)
    x1 = (ice_depth - depth) * tan(theta_emit)
    x2 = ice_depth * tan(theta_emit)
    source = candidate.surface - (2.0 * ice_depth - depth) * candidate.flat_normal -
        (x1 + x2) * tangent_surface
    source_normal = source / norm(source)
    tangent_source = candidate.view - (candidate.view ⋅ source_normal) * source_normal
    tangent_source /= norm(tangent_source)
    direction = -cos(theta_emit) * source_normal + sin(theta_emit) * tangent_source
    return direction / norm(direction), source_normal
end

function parse_bool_modes(value)
    modes = Bool[]
    for token in split(value, ',')
        normalized = lowercase(strip(token))
        normalized == "reject" && push!(modes, true)
        normalized == "allow" && push!(modes, false)
        normalized in ("reject", "allow") || error("--backfacing-modes must use reject and/or allow")
    end
    isempty(modes) && error("--backfacing-modes cannot be empty")
    return unique(modes)
end

backfacing_mode_name(reject) = reject ? "physical_reject" : "counterfactual_allow"
