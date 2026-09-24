#!/usr/bin/env julia

import Pkg
# Activate the project relative to this script so this also works when Slurm
# starts the job from a different directory.
Pkg.activate(normpath(joinpath(@__DIR__, "..")))

using CoRaLS
using Unitful: km, m, sr, EeV, MHz, cm, μV, °
using Unitful
using Dates
using CSV
using DataFrames
using Interpolations

#–– Parse command-line arguments ––#
if length(ARGS) ∉ (8, 9)
    println(stderr, "Usage: julia generic_acceptance.jl <altitude_index> <energy_multiplier> <ice_depth_m> <antennas> <triggers> <pointing_angle_deg> <frequency_min_MHz> <trials_exponent> [slope_model]")
    println(stderr, "Slope models: ", join(available_slope_models(), ", "))
    exit(1)
end

# convert to Unitful quantities
altitude = 5*parse(Float64, ARGS[1])km
energyMult = parse(Float64, ARGS[2])
ice_depth = parse(Float64, ARGS[3])m
antNum = parse(Int, ARGS[4])
trigNum = parse(Int, ARGS[5])
angle = parse(Float64, ARGS[6])
freqMin = parse(Float64, ARGS[7])MHz
ntrials = 10^parse(Int, ARGS[8])
slope_model_name = length(ARGS) == 9 ? ARGS[9] : "gaussian_0"
slope_model = slope_model_from_name(slope_model_name)
ENERGY1 = 0.1 * energyMult * EeV
ENERGY2 = 100 * energyMult * EeV

# Array rate sweeps should not collect or overwrite large event files. To save
# triggered events deliberately, set CORALS_SAVE_EVENTS=true and provide a
# unique CORALS_SAVEFILE path in the Slurm submission environment.
save_events = lowercase(get(ENV, "CORALS_SAVE_EVENTS", "false")) in ("1", "true", "yes")
savefile = get(ENV, "CORALS_SAVEFILE", "")
if save_events && isempty(savefile)
    @warn "CORALS_SAVE_EVENTS is true but CORALS_SAVEFILE is empty; events will not be written."
end

#df = CSV.read(joinpath(@__DIR__, "../data/Mare_Cuboid_Efield_linear_perm.csv"), DataFrame)
df = CSV.read(joinpath(@__DIR__, "../data/allOnes.csv"), DataFrame)
satten = linear_interpolation(df[:, 1], df[:, 2])

#–– Set up your run ––#
#region   = create_region("polar:south,-80,0.0557")
#region   = create_region("polar:south,-85,0.114")
#region   = create_region("polar:south,-80,0.1600")
#region   = create_region("polar:south,-80,1.0")
region = PSRRegion()
#region = NonPSRRegion()

#region = PSR_NRegion()
#region = Non_PSR_NRegion()

#region = MareRegion()
#region = HighlandRegion()
#function always_true((lat, lon))
#    return true
#end
#WholeMoonMare = CustomRegion(always_true, 0.162, 6.1512e6km^2)
#region   = WholeMoonHighlands
#region = WholeMoonMare
#sc       = CircularOrbit(altitude)
sc       = FixedPlatform(-80, 0, 50km)
trigger  = LPDA(Nant=antNum, Ntrig=trigNum, θ0=angle, altitude=altitude, skyfrac=0.15)
#trigger  = ANITA(Nant=antNum, Ntrig=trigNum, θ0=angle, altitude=altitude, skyfrac=0.15)
#trigger  = magnitude_trigger(100μV / m)
## Potentially pass phi0 (with actual phi character) for aligning antennas
kws      = Dict(
    :min_energy => ENERGY1, #30.0EeV,
    :max_energy => ENERGY2,#31.0EeV,
    :ν_min      => freqMin,
    :ν_max      => 800MHz,
    :dν         => 10MHz,
    :min_count  => 10, # 50 -> 10
    :max_tries  => 10, # 100 -> 10
    :simple_area=> false,
    :tand_mag=> 0.000,
    :tanδnorm=> 0.000,
    :slopemodel=> slope_model,
    :roughnessmodel=> GaussianRoughness(0),
    :iceroughness=> GaussianIceRoughness(0.0cm),
    :low_temp_corr_factor=> 1.0,
    :satten=>satten,
    #:ϕ0=>pi/2
)

if altitude <= 10km
				ntrials *= 30
else
				ntrials *= 3
end

nbins   = 16
println("ntrials: ", ntrials)

#–– Run acceptance ––#
A = acceptance(ntrials, nbins;
    region=region,
    spacecraft=sc,
    trigger=trigger,
    ice_depth=ice_depth,
    ice_thickness=1.0m,
    Nice=1.265, ##regolith_index(StrangwayIndex(), ice_depth), ## Peter on call (Sept. 11, 2025) said ice permittivity should be 2.2
    Nbed=regolith_index(LSB_DivinerIndex(), ice_depth), # For bedrock: use Nbed=2.56, this is based on anorthosite P. Linton uses for rocks in volume scattering sims
    kws...,
    indexmodel=LSB_DivinerIndex(), # MACHTAY try this for changing index of refraction,
    densitymodel=LSB_Diviner_Density(),
    #indexmodel=ConstantIndex(), # MACHTAY try this for changing index of refraction
    save_events=save_events,
    savetriggered=save_events,
    savefile=savefile,
)

#–– Define MCSE ––#
function mcse(count::Int, n::Int)
    p = count / n
    return sqrt(p * (1 - p) / n)
end

#–– Compute spectra + errors ––#
r_spectra = differential_spectrum(A.energies, A.rAΩ, 1yr)
r_error   = differential_spectrum(A.energies, mcse.(A.rcount, A.ntrials) .* A.gAΩ, 1yr)

d_spectra = differential_spectrum(A.energies, A.dAΩ, 1yr)
d_error   = differential_spectrum(A.energies, mcse.(A.dcount, A.ntrials) .* A.gAΩ, 1yr)

test_spectra = get_spectrum(A.energies, A.rAΩ, 1yr)
test_spectra_error = get_spectrum(A.energies, mcse.(A.rcount, A.ntrials) .* A.gAΩ, 1yr)
#–– Print results ––#
timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
println("# Run on $timestamp")
#println("# altitude = $(altitude), ice_depth = $(ice_depth)")
#println("r_spectra = ", r_spectra)
#println("r_error   = ", r_error)
#println("d_spectra = ", d_spectra)
#println("d_error   = ", d_error)
#println("alt (km), ice_depth (m), r_count, r_err, d_count, d_err")
println("Energy (EeV), Altitude (km), Ice Depth (m), Slope Model, Reflected Rate, Reflected Error, Direct Rate, Direct Error, ARW Reflected Rate, ARW Reflected Error")
for i in 1:length(r_spectra)
    println(ustrip(A.energies[i]), ", ", ustrip(altitude), ", ", ustrip(ice_depth), ", ", slope_model_name, ", ", r_spectra[i], ", ", r_error[i], ", ", d_spectra[i], ", ", d_error[i], ", ", test_spectra[i], ", ", test_spectra_error[i])
end
