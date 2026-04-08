#!/usr/bin/env julia

import Pkg
# Activate the CoRaLS.jl project
Pkg.activate("/users/PAS0654/machtay1/CoRaLS/main/CoRaLS.jl")    # adjust if your script lives elsewhere

using CoRaLS
using Unitful: km, m, sr, EeV, MHz
using Unitful
using Dates

#–– Parse command-line arguments ––#
if length(ARGS) != 8
    println(stderr, "Usage: julia acceptance.jl <altitude_km> <ice_depth_m>")
    exit(1)
end

# convert to Unitful quantities
altitude = 5*parse(Float64, ARGS[1])km    ## Altitude of detector above surface
energyMult = parse(Float64, ARGS[2])      ## Scales up energy range (usually set to 1)
ice_depth = parse(Float64, ARGS[3])m      ## Ice depth in meters
antNum = parse(Int, ARGS[4])              ## Numer of Antennas
trigNum = parse(Int, ARGS[5])             ## Number of triggering antennas required 
angle = parse(Float64, ARGS[6])           ## Pointing angle from detector of antennas (-90 = straight down)
freqMin = parse(Float64, ARGS[7])MHz      ## Lowest sensitive frequency (change based on antenna size)
ntrials = 10^parse(Int, ARGS[8])          ## Number of trials for each energy bin to try
ENERGY1 = 0.01 * energyMult * EeV         ## Lowest energy bin
ENERGY2 = 100 * energyMult * EeV          ## Highest energy bin

#–– Set up your run ––#
region   = create_region("polar:south,-80,0.0557")
sc       = CircularOrbit(altitude)
trigger  = LPDA(Nant=antNum, Ntrig=trigNum, θ0=angle, altitude=altitude, skyfrac=0)
## Potentially pass phi0 (with actual phi character) for aligning antennas
kws      = Dict(
    :min_energy => ENERGY1, #30.0EeV,
    :max_energy => ENERGY2,#31.0EeV,
    :ν_min      => freqMin,
    :ν_max      => 1000MHz,
    :dν         => 30MHz,
    :min_count  => 10, # 50 -> 10
    :max_tries  => 10, # 100 -> 10
    :simple_area=> true,
)

if altitude <= 10km
				ntrials *= 30
else
				ntrials *= 3
end

nbins   = 8
println("ntrials: ", ntrials)

#–– Run acceptance ––#
A = acceptance(ntrials, nbins;
    region=region,
    spacecraft=sc,
    trigger=trigger,
    ice_depth=ice_depth,
    ice_thickness=1.0m,
    Nice=1.350, 
    #Nice=1.483, 
    #Nice=regolith_index(StrangwayIndex(), ice_depth), ## Peter on call (Sept. 11, 2025) said ice permittivity should be 2.2
    #Nbed=regolith_index(StrangwayIndex(), ice_depth), # For bedrock: use Nbed=2.56, this is based on anorthosite P. Linton uses for rocks in volume scattering sims
    kws...,
    indexmodel=StrangwayIndex(), # MACHTAY try this for changing index of refraction,
    #densitymodel=StrangwayDensity(),
    #indexmodel=ConstantIndex(), # MACHTAY try this for changing index of refraction
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

#–– Print results ––#
timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
println("# Run on $timestamp")
#println("# altitude = $(altitude), ice_depth = $(ice_depth)")
#println("r_spectra = ", r_spectra)
#println("r_error   = ", r_error)
#println("d_spectra = ", d_spectra)
#println("d_error   = ", d_error)
#println("alt (km), ice_depth (m), r_count, r_err, d_count, d_err")
println("Energy (EeV), Altitude (km), Ice Depth (m), Angle (deg), Reflected Count, Reflected Error, Direct Count, Direct Error")
for i in 1:length(r_spectra)
    println(ustrip(A.energies[i]), ", ", ustrip(altitude), ", ", ustrip(ice_depth), ", ", ustrip(angle), ", ", r_spectra[i], ", ", r_error[i], ", ", d_spectra[i], ", ", d_error[i])
end
