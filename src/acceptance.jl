using JLD2
using ProgressMeter
"""
    acceptance.jl

This file is part of the `CoRaLS` module and contains functions and structures for calculating the acceptance of the system for sub-surface Ultra-High-Energy Cosmic Ray (UHECR) reflections from ice. It also includes utilities for analyzing trials and computing differential spectra.

## Main Components

- `struct Acceptance`: A structure to hold results of acceptance calculations, including trials, altitude, energies, and different acceptance types (direct and reflected).

- `function acceptance(...)`: Calculates the acceptance of CoRaLS for sub-surface UHECR reflections. It involves simulations of cosmic ray interactions, triggering conditions, and aggregating results across multiple trials and energy bins. 

- `function trials_passed(...)`: Computes the number of successful trials for both direct and reflected acceptances. 

- `function differential_spectrum(...)`: Calculates the differential spectrum of detected UHECR events based on acceptance calculations and observation time.
"""

"""
A structure to hold results of acceptance calculations, including trials, altitude, energies, and different acceptance types (direct and reflected).
"""
mutable struct Acceptance
    ntrials
    region::AbstractRegion
    spacecraft::Spacecraft
    energies
    gAΩ
    dAΩ
    rAΩ
    dcount
    rcount
    failtypes
    dfailed
    rfailed
    direct
    reflected
end

"""
acceptance(ntrials, nbins; region, spacecraft, trigger, ...)

Calculates the acceptance of CoRaLS for sub-surface UHECR reflections. It involves simulations of cosmic ray interactions, triggering conditions, and aggregating results across multiple trials and energy bins.
"""
function acceptance(ntrials::Int, nbins::Int; 
    region::AbstractRegion=AllPSR,
    spacecraft::Spacecraft=CircularOrbit(50.0km),
    trigger=magnitude_trigger(67μV / m),
    min_energy=0.1EeV,
    max_energy=600.0EeV,
    min_count=10,
    max_tries=1,
    save_events=false,
    simple_area=false,
    savefile="",
    kwargs...
    )
    # Init arrays
    dcount = zeros(Int, nbins)
    rcount = zeros(Int, nbins)
    dfailed = zeros(Int, nbins, length(instances(TrialFailed)))
    rfailed = zeros(Int, nbins, length(instances(TrialFailed)))
    tries_left = ones(Int, nbins) .* max_tries  # Num times to retry remaining
    devents = StructArray{Direct}(undef, 0)
    revents = StructArray{Reflected}(undef, 0)
    energies = (10.0 .^ range(log10(min_energy / 1.0EeV), log10(max_energy / 1.0EeV), length=nbins + 1))EeV

    @info "Calculating acceptance using $(ntrials) trials across $(nbins) bins with up to $(max_tries-1) retries"

    @showprogress 1 "Simulating..." for bin = 1:nbins
        # Run acceptance loop at least once, retry if min_count not met up to max_tries times
        while tries_left[bin] > 0 && (tries_left[bin] == max_tries || rcount[bin] < min_count)
            # loop over the number of trials in this bin
            for i = 1:ntrials
                # throw a random cosmic ray trial and get the signal at the payload
                direct, reflected = throw_cosmicray(
                    #sample_auger(energies[bin], energies[bin+1]),
                    energies[bin], 
                    trigger, region, spacecraft; simple_area=simple_area, kwargs...)
                if direct isa TrialFailed
                    dfailed[bin, Int(direct)] += 1
                else
                    if trigger(direct)
                        direct.triggered = true
                        dcount[bin] += 1
                    end
                    # Save all events that get to payload regardless of trigger
                    save_events && push!(devents, direct)
                end
                
                if reflected isa TrialFailed
                    rfailed[bin, Int(reflected)] += 1
                else
                    if trigger(reflected)
                        reflected.triggered = true
                        rcount[bin] += 1
                    end
                    save_events && push!(revents, reflected)
                end
            end # end ntrials loop
            tries_left[bin] -= 1
        end # end retry loop
    end # end nbins loop

    # Compute acceptance (pi * A_collected) sr [km^2 sr]
    #  factor of pi from integral(cos(theta)) possible angles of incoming CRs
    #  A_collected : area over which CRs are sampled (usually whole Moon)
    #  Other acceptance factors are computed by MC rejection (visibility, triggering, etc)
    if simple_area && region isa AOIRegion
        # We threw CRs only in the AOIRegion, so A_collected is that AOI surface area
        # Note: Visibility to spacecraft is still computed via rejection sampling
        gAΩ = pi * sr * region.aoi_frac * aoi_area(region.aoi)
        #gAΩ = pi * sr * aoi_area(region.aoi)
    else
        # Otherwise we threw cosmic rays for whole moon and computed area with MC
        gAΩ = pi * sr * region_area(WholeMoonRegion())  # [km^2 sr]
    end

    # Count fraction of passed trials per E bin
    ntrials_per_bin = ntrials .* (max_tries .- tries_left)
    dAΩ = gAΩ .* (dcount ./ ntrials_per_bin)  # [km^2 sr]
    rAΩ = gAΩ .* (rcount ./ ntrials_per_bin)  # [km^2 sr]

    AΩ = Acceptance(ntrials_per_bin, region, spacecraft, energies, gAΩ, dAΩ, rAΩ, dcount, rcount, instances(TrialFailed), dfailed, rfailed, devents, revents)
    
    if !isempty(savefile)
        save_acceptance(AΩ, savefile)
    end

    return AΩ
end


struct OldAcceptance
    ntrials::Int64 # per bin
    altitude
    gAΩ
    energies
    dAΩ
    direct
    rAΩ
    reflected
end

"""
old_acceptance(ntrials, nbins; ...)

Calculates the acceptance of CoRaLS for sub-surface UHECR reflections. It involves simulations of cosmic ray interactions, triggering conditions, and aggregating results across multiple trials and energy bins.
"""
function old_acceptance(ntrials, nbins; min_energy=0.1EeV,
    max_energy=600.0EeV,
    altitude=20.0km,
    trigger=magnitude_trigger(100μV / m),
    save_events=true,
    kwargs...)

    # create a StructArray for out Detection type
    # this is a vector with the length of the number of threads
    # we concat them together at the end to make a single collection of events
    devents = Vector{StructArray{Direct}}()
    revents = Vector{StructArray{Reflected}}()

    # setup the struct arrays for each thread
    for i in 1:Threads.nthreads()
        push!(devents, StructArray{Direct}(undef, 0))
        push!(revents, StructArray{Reflected}(undef, 0))
    end

    @info "Calculating acceptance using $(ntrials) trials across $(nbins) bins..."

    # we define the lunar poles as everything within 10 degrees
    θpole = deg2rad(10.0)

    # the distance to the horizon from the SC's altitude
    θmax = -horizon_angle(altitude)

    # we sample events only in the lunar poles
    A = spherical_cap_area(θpole, Rmoon)

    # the maximum total acceptance is πA km^2 steradians
    # however the spacecraft is in the *opposite* hemisphere from the
    # cosmic ray trial 50% of the time, so we pick up another factor of 0.5
    # since we are only simulating one "combined" pole (to be more efficient).
    AΩ = 0.5 * π * sr * A * ones(nbins)

    # fold in the fraction of the polar caps that have PSRs
    # assuming 30_000 km^2 of total PSR regions
    AΩ *= 30_000km^2 / spherical_cap_area(θpole, Rmoon)  # TODO: 30000 is both poles, divinding by only 1 sph cap 

    # and account for the fraction of a single hemisphere of the CRs
    # orbit that is visible from the poles, since we only simulate these
    # "visible" SC angles
    AΩ *= (θpole + θmax) / (π / 2.0)  # TODO:rederive this - max area visible by payload

    # construct the bins we sample from in log-space
    energies = (10.0 .^ range(log10(min_energy / 1.0EeV), log10(max_energy / 1.0EeV), length=nbins + 1))EeV

    # construct arrays containing the effective area for direct
    # and reflected events at each energy
    dAΩ = AΩ .* ones(nbins)
    rAΩ = AΩ .* ones(nbins)

    # try the each number of event that was requested
    @showprogress 1 "Simulating..." for bin = 1:nbins

        # the number of trials that have passed for direct and reflected
        dpassed = Threads.Atomic{Int}(0.0)
        rpassed = Threads.Atomic{Int}(0.0)

        # sample the Auger spectrum for UHECR energies within this bin
        Ecr = [sample_auger(energies[bin], energies[bin+1]) for i = 1:ntrials]

        # loop over the number of trials in this bin
        Threads.@threads for i = 1:ntrials

            # throw a random cosmic ray trial and get the signal at the payload
            direct, reflected = throw_cosmicray(Ecr[i], trigger; altitude=altitude, kwargs...)

            # check for a direct trigger
            if !(direct isa TrialFailed)
                if trigger(direct)
                    direct.triggered = true
                    Threads.atomic_add!(dpassed, 1)
                end
                save_events && push!(devents[Threads.threadid()], direct)
            end

            # check for a reflected trigger
            if !(reflected isa TrialFailed)
                if trigger(reflected)
                    reflected.triggered = true
                    Threads.atomic_add!(rpassed, 1)
                end
                save_events && push!(revents[Threads.threadid()], reflected)
            end
        end # end loop over trials

        # calculate the acceptance at this energy
        dAΩ[bin] *= dpassed[] / ntrials
        rAΩ[bin] *= rpassed[] / ntrials

    end # end loop over events

    # and combine the event structs across all threads
    devents = vcat(devents...)
    revents = vcat(revents...)

    return OldAcceptance(ntrials, altitude, AΩ, energies, dAΩ, devents, rAΩ, revents)

end


"""
    trials_passed

Calculate the number of successful trials in the direct and reflected acceptances.
"""
trials_passed(AΩ) = (AΩ.ntrials * (AΩ.dAΩ ./ AΩ.gAΩ), AΩ.ntrials * (AΩ.rAΩ ./ AΩ.gAΩ))


"""
    differential_spectrum(energies, AΩ, T)

Calculate the differential spectrum of detected UHECR events based on acceptance calculations and T, the total observation time of the experiment.

Calculates the differential spectrum of detected UHECR events .
"""
function differential_spectrum(energies, AΩ, T)

    # array to store the spectrum in each bin
    spectrum = zeros(length(AΩ))

    # loop over each bin
    for bin = 1:(length(energies)-1)

        # calculate the average of the auger spectrum in this bin
        avg = mean(auger_spectrum.(range(energies[bin], stop=energies[bin+1], length=100)))

        # calculate dN/DE using this average value of the spectrum
        # spectrum[bin] = (T * AΩ[bin] * avg) |> (EeV^(-1))
        spectrum[bin] = (T * AΩ[bin] * avg * (energies[bin+1] - energies[bin])) |> NoUnits

    end

    # and we are done
    return spectrum
end

# I/O helpers
"""
    save_acceptance(A::Acceptance, filename::String)

Save the acceptance data to a file (only .jld2 supported).
"""
function save_acceptance(A::Acceptance, filename::String)
    if filename[end-4:end] != ".jld2"
        throw("Only .jld2 supported.")
    end
    @save filename A
end

"""
    load_acceptance(filename::String)::Acceptance

Load the acceptance data from a .jld2 file.
"""
function load_acceptance(filename::String)::Acceptance
    @load filename A
    return A
end

"""
    merge_acceptance(A::Acceptance, B::Acceptance)::Acceptance

Merge two Acceptance structures into one, combining their trial counts and results.
"""
function merge_acceptance(A::Acceptance, B::Acceptance)

    # Check if Acceptance structs are compatible
    # TODO: need to include more metadata in struct to say this for certain
    if !((A.region == B.region) &  
        (A.spacecraft == B.spacecraft) &
        (A.energies == B.energies) &
        (A.gAΩ == B.gAΩ) & 
        (A.failtypes == B.failtypes))
        throw("Acceptances not compatible.")
    end

    C = deepcopy(A)
    C.ntrials = A.ntrials .+ B.ntrials
    C.dcount = A.dcount .+ B.dcount
    C.rcount = A.rcount .+ B.rcount
    C.dfailed = A.dfailed .+ B.dfailed
    C.rfailed = A.rfailed .+ B.rfailed
    C.direct = vcat(A.direct, B.direct)
    C.reflected = vcat(A.reflected, B.reflected)

    # Compute new acceptance from combined ntrials and counts
    C.dAΩ = C.gAΩ .* (C.dcount ./ C.ntrials)
    C.rAΩ = C.gAΩ .* (C.rcount ./ C.ntrials)
    return C
end