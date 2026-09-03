using DelimitedFiles
using Unitful: Ω, °, MHz

"""
    detector.jl

This file is part of the `CoRaLS` module, focusing on modeling and simulating antenna-based detector systems for cosmic ray detection. It includes functions for constructing antenna simulations and specific models for different types of antennas, such as LPDAs and ANITA horns, based on their antenna factor data and operational parameters.

## Main Components
- `create_antenna`: Constructs an antenna simulation based on antenna factor data. This function sets up the parameters for an antenna system, including its frequency range, noise characteristics, and response factors.

- `LPDA`: Creates a model for the Low-Profile Dipole Array (LPDA) used in cosmic ray detection simulations. Configures the LPDA based on given parameters or defaults.

- `ANITA`: Creates a model for the ANtarctic Impulsive Transient Antenna (ANITA). Configures the ANITA system based on given parameters or defaults.
"""

"""
    create_antenna(fname, skyfrac; kwargs...)

Construct an antenna-factor based antenna simulation. This function sets up the parameters for an antenna system, including its frequency range, noise characteristics, and response factors, based on antenna factor data from a file.

# Arguments
- `fname`: Filename containing antenna factor data.
- `skyfrac`: Fraction of the antenna's view that is sky.
- `kwargs`: Additional keyword arguments including altitude, frequency range, angular resolutions, and trigger settings.

# Returns
- A function that evaluates the antenna trigger condition for a given event.
"""
function create_antenna(fname;
    skyfrac=0.03, altitude=20.0km,
    ν_min=150MHz, ν_max=800MHz,
    σθ=29.09, σϕ=31.28, θ0=-40.0, ϕ0=-1.0,
    Nant=8, Ntrig=3, Tlna=100K, Tmoon=85K, SNR=4.0)

    # current LPDA design is 150-600 MHz
    # ANITA horns are 100-800 MHz
    ν_ant = ν_min:10MHz:ν_max

    # the antenna-receiver matching
    Γ = 0.9

    Tant = antenna_temp(
    ν_ant,
    50km,
    θ0 * u"°",
    0°,
    29.09°,
    31.28°,
    sky_temperature(ν_ant);
    T_moon = 85K,
    )

    # construct the system temperature
    #Tsys = Tlna .+ Γ * (Tsky .* skyfrac .+ Tmoon * moonfrac)
    Tsys = Tlna .+ Γ*Tant

    # the current LPDA design is nominally a 120 ohm design
    # but the antenna factors are referenced to 50 ohm so
    # calculate the noise at 50 ohm directly
    Z = 50Ω

    # load the antenna factor from files
    AFs = readdlm("$(@__DIR__)/../data/$(fname).txt", skipstart=1)

    # create an interpolator and interpolate it onto our frequencies
    AF = LinearInterpolation(AFs[:, 1] ./ 1e6, AFs[:, 2])(ν_ant ./ MHz) ./ m

    # construct the noise voltage spectrum - this is units of Volts
    Vn = sqrt.(2* Ntrig * k_b .* Tsys .* Z * (ν_ant[2] - ν_ant[1])) .|> μV

    # we also need the total summed noise voltage
    Vntot::typeof(1.0μV) = sqrt(sum(Vn .^ 2))
    
    #Vntot::typeof(1.0μV) = sqrt(2*Ntrig*k_b*(60 + 0.9*85)K*Z*132.687MHz) |> μV

    # construct the azimuthal angles given the number of antennas
    ϕ = range(0, 2π, length=Nant + 1)[1:Nant]
    if ϕ0 > 0
        # If ϕ0 is given, align all antennas along the same azimuthal direction
        ϕ = ϕ0 * ones(Nant)
    end

    #θ0 = [-55°, -50°, -55°, -50°, -55°, -50°, -55°, -50°]
    
    # construct the boresight vectors for each antenna
    boresight = spherical_to_cartesian.(π / 2.0 .+ deg2rad(-θ0) * ones(Nant), ϕ, 1)
    #boresight = spherical_to_cartesian.(π / 2.0 .+ deg2rad.(-θ0), ϕ, 1)

    # construct the perpendicular (i.e. H-pol) antenna axis
    Haxis = [boresight[i] × SA[0, 0, 1] for i = 1:Nant]
    Haxis ./= norm.(Haxis) # make sure it is normalized

    # and now cross the Haxis with the boresight to get the Vaxis
    Vaxis = Haxis .× boresight
    Vaxis ./= norm.(Vaxis) # make sure it is normalized

    # we now construct the function that evaluates the trigger
    # use a let-block to improve performance on earlier Julia versions
    trigger = let ν_ant = ν_ant, Vntot = Vntot, AF = AF, Nant = Nant,
        Ntrig = Ntrig, Haxis = Haxis, Vaxis = Vaxis, boresight = boresight

        function (event::AbstractSignal)

            # make sure our simulated LPDA frequencies are
            # the same as the electric field frequencies
            @assert (ν_ant[2] - ν_ant[1]) == event.dν

            # interpolate the electric field onto the antennas frequencies
            # VmMHz = (1V/m/MHz) .* LinearInterpolation((event.ν_min:event.dν:event.ν_max) ./ MHz,
            #                                           event.Ef ./ (V/m/MHz))(ν_ant ./ MHz)
            VmMHz = (1μV / m / MHz) .* extrapolate(interpolate(((event.ν_min:event.dν:event.ν_max) ./ MHz,),
                    event.Ef ./ (μV / m / MHz),
                    Gridded(Linear())), 0.0)(ν_ant ./ MHz)

            # and calculate the electric field spectrum in V/m
            #Ef = norm(event.pol) * VmMHz .* (ν_ant[2] - ν_ant[1]) # now in Volts/m
            ## VmMHz is a spectral density, so multiply by frequency bin size to get V/m
            Ef = VmMHz .* (ν_ant[2] - ν_ant[1]) # now in Volts/m

            # construct the vector to the event from SC
            ## ant_view takes in the lunar centric coordinates for the event and the spacecraft and 
            ## returns the view vector in spacecraft coordinates.
            antenna_view = ant_view((event.sc_r, event.sc_θ, event.sc_ϕ), (Rmoon, deg2rad(event.θ), deg2rad(event.ϕ)))
            view = antenna_view.view
            R = antenna_view.rotation
            
            ## This is the old way of calculating the view vector. However this uses the events ϕ
            ## which is a lunar centric value, but θ_el is a spacecraft centric value.
            
            #view = spherical_to_cartesian.(π / 2.0 + deg2rad(-event.θ_el),
            #    deg2rad(event.ϕ), 1)

            # this array stores the relative antenna response factors due
            # to the off-axis beam response
            scale = ones(Nant)

            # Scalar off-axis response factors for each antenna
            scale_V = ones(Float64, Nant)
            scale_H = ones(Float64, Nant)

            # Each antenna receives a full frequency spectrum
            Nfreq = length(Ef)

            V_H = Matrix{typeof(Ef[1] / AF[1])}(undef, Nfreq, Nant)
            V_V = Matrix{typeof(Ef[1] / AF[1])}(undef, Nfreq, Nant)
            V_CP = Matrix{typeof(Ef[1] / AF[1])}(undef, Nfreq, Nant)

            # loop over every antenna to calculate θ, ϕ and then the off-axis weighting
            for i = 1:Nant

                # project the view vector onto the plane defined by Haxis
                # i.e. the vertical plane centered at the boresight
                Vproj = view - (view ⋅ Haxis[i]) * Haxis[i]

                # the angle between `Vproj` and boresight is the elevation angle
                θev = acos(clamp((Vproj / norm(Vproj)) ⋅ boresight[i], -1, 1)) |> rad2deg

                # project the view vector onto the plane defined by Vaxis
                # i.e. the horizontal plane centered at the boresight
                Hproj = view - (view ⋅ Vaxis[i]) * Vaxis[i]

                # the angle between `Hproj` and boresight is the azimuth angle
                ϕev = acos(clamp((Hproj / norm(Hproj)) ⋅ boresight[i], -1, 1)) |> rad2deg


                EH, EV, _, _, _ = project_event_pol(Ef, event.pol, R, Haxis[i], Vaxis[i])
                # apply the off-axis antenna response at these angles
                #scale[i] *= 0.074 + (1.0 - 0.074)exp(-θev * θev / (2σθ * σθ))
                #scale[i] *= 0.062 + (1.0 - 0.062)exp(-ϕev * ϕev / (2σϕ * σϕ))

                ## PL NOTE: My XF sims of the CoRaLS fullscale antenna seem to imply theta_ev and phi_ev were flipped

                scale[i] *= 0.0125 + (1.0 - 0.0125)exp(-θev * θev / (2σθ * σθ))
                scale[i] *= 0.08395 + (1.0 - 0.08395)exp(-ϕev * ϕev / (2σϕ * σϕ))

                scale_V[i] *= 0.0125 + (1.0 - 0.0125)exp(-θev * θev / (2σθ * σθ))
                scale_V[i] *= 0.08395 + (1.0 - 0.08395)exp(-ϕev * ϕev / (2σϕ * σϕ))

                #Rotating the polarization by 90 deg gives flips the off axis response
                scale_H[i] *= 0.0125 + (1.0 - 0.0125)exp(-ϕev * ϕev / (2σθ * σθ))
                scale_H[i] *= 0.08395 + (1.0 - 0.08395)exp(-θev * θev / (2σϕ * σϕ))

                ## Determine the voltage on vertical and horizontal RF lines
                V_V[:, i] .= EV.*scale_V[i] ./ AF
                V_H[:, i] .= EH.*scale_H[i] ./ AF

                ## When you form the LCP-RCP basis on a purely linear signal
                ## you get equal power in LCP RCP. Subtracting them will
                ## give an additional factor of 2, as it is analagous
                ## to summing twice the number of antennas
                V_CP[:, i] .= 2*hypot.(V_H[i], V_V[i]) / sqrt(2)

            end

            #####################################################
            #####################################################

            ## Old way of calculating the signal voltage peak on the antenna

            # we now pick the best Ntrig antennas to calculate the trigger with
            #Nsum = sum(sort!(scale, rev=true)[1:Ntrig])

            # apply the antenna factors to the electric field
            # to convert it to Volts, and scale it by the number of antennas
            # double the field because we have two polarizations per antenna
            #Vs = norm(event.pol)*(2.0 * Nsum) * (Ef ./ AF)

            # and lastly sum up the electric field
            #Vf = sum(Vs)

            #######################################################
            #######################################################

            #Vf = best_adjacent_triplet_value(V_CP).magnitude
            Vf = best_adjacent_antennas(V_CP, Ntrig).magnitude .|> μV

            # we now pick the best Ntrig antennas to calculate the trigger with
            #Vf = sum(sort!(V_CP, rev=true)[1:Ntrig])

            # apply the antenna factors to the electric field
            # to convert it to Volts, and scale it by the number of antennas
            # double the field because we have two polarizations per antenna
            #Vs = (2.0 * Nsum) * (Ef ./ AF)

            # and lastly sum up the electric field
            #Vf = sum(Vs)

            # do the same process above for the reflection from the bottom of the layer
            if typeof(event) == Reflected
                scale = ones(Nant)

                # Scalar off-axis response factors for each antenna
                scale_V = ones(Float64, Nant)
                scale_H = ones(Float64, Nant)

                # Each antenna receives a full frequency spectrum
                Nfreq = length(Ef)

                V_H = Matrix{typeof(Ef[1] / AF[1])}(undef, Nfreq, Nant)
                V_V = Matrix{typeof(Ef[1] / AF[1])}(undef, Nfreq, Nant)
                V_CP = Matrix{typeof(Ef[1] / AF[1])}(undef, Nfreq, Nant)

                for i = 1:Nant

                    # project the view vector onto the plane defined by Haxis
                    # i.e. the vertical plane centered at the boresight
                    Vproj = view - (view ⋅ Haxis[i]) * Haxis[i]

                    # the angle between `Vproj` and boresight is the elevation angle
                    θev = acos(clamp((Vproj / norm(Vproj)) ⋅ boresight[i], -1, 1)) |> rad2deg

                    # project the view vector onto the plane defined by Vaxis
                    # i.e. the horizontal plane centered at the boresight
                    Hproj = view - (view ⋅ Vaxis[i]) * Vaxis[i]

                    # the angle between `Hproj` and boresight is the azimuth angle
                    ϕev = acos(clamp((Hproj / norm(Hproj)) ⋅ boresight[i], -1, 1)) |> rad2deg


                    EH, EV, _, _, _ = project_event_pol(Ef, event.polsub, R, Haxis[i], Vaxis[i])
                    # apply the off-axis antenna response at these angles
                    #scale[i] *= 0.074 + (1.0 - 0.074)exp(-θev * θev / (2σθ * σθ))
                    #scale[i] *= 0.062 + (1.0 - 0.062)exp(-ϕev * ϕev / (2σϕ * σϕ))

                    ## PL NOTE: My XF sims of the CoRaLS fullscale antenna seem to imply theta_ev and phi_ev were flipped

                    scale[i] *= 0.0125 + (1.0 - 0.0125)exp(-θev * θev / (2σθ * σθ))
                    scale[i] *= 0.08395 + (1.0 - 0.08395)exp(-ϕev * ϕev / (2σϕ * σϕ))

                    scale_V[i] *= 0.0125 + (1.0 - 0.0125)exp(-θev * θev / (2σθ * σθ))
                    scale_V[i] *= 0.08395 + (1.0 - 0.08395)exp(-ϕev * ϕev / (2σϕ * σϕ))

                    #Rotating the polarization by 90 deg gives flips the off axis response
                    scale_H[i] *= 0.0125 + (1.0 - 0.0125)exp(-ϕev * ϕev / (2σθ * σθ))
                    scale_H[i] *= 0.08395 + (1.0 - 0.08395)exp(-θev * θev / (2σϕ * σϕ))

                    ## Determine the voltage on vertical and horizontal RF lines
                    V_V[:, i] .= EV.*scale_V[i] ./ AF
                    V_H[:, i] .= EH.*scale_H[i] ./ AF

                    V_CP[:, i] .= 2*hypot.(V_H[i], V_V[i]) / sqrt(2)

                end

                #Vf = best_adjacent_triplet_value(V_CP).magnitude
                Vfsub = best_adjacent_antennas(V_CP, Ntrig).magnitude .|> μV






                # and calculate the electric field spectrum for the bottom reflection in V/m
                #Efsub = norm(event.polsub) * VmMHz .* (ν_ant[2] - ν_ant[1]) # now in Volts/m

                # and calculate the corresponding voltage at the receiver
                #Vssub = (2.0 * Nsum) * (Efsub ./ AF)

                # and sum it up for the trigger
                #Vfsub = sum(Vssub)
                #Vfsub = 0.0μV
            else
                Vfsub = 0.0μV
            end

            # realize a random signal+noise spectrum
            # this is in volts already, not V/MHz
            # throw for two different events - surface and subsurface
            # Vf = rician(Vs ./ μV, Vn ./ μV) .* μV
            # Vfsub = rician(Vssub ./ μV, Vn ./ μV) .* μV

            # throw for a Rician with signal and noise
            # println(sum(Vs))
            # println(Vntot)
            # Vf = rician( sum(Vs) / μV,  Vntot / μV)μV
            # println(Vf)
            # Vfsub = 0.0μV

            # println(sum(Vs) / Vntot)

            # println(sum(Vs))
            # println(sum(Vssub))
            # println(Vntot)

            # check that the total integrated voltage is greater
            # than the SNR times the summed noise voltage
            # we only need to check the greater of the top
            # and bottom surface reflections
            # return max(sqrt(sum(Vf.^2)), sqrt(sum(Vfsub.^2))) > (SNR*Vntot)
            return max(Vf, Vfsub) > (SNR * Vntot)
            # return sum(Vs) > (SNR * Vntot)
            #return max(sum(Vs), sum(Vssub)) > (SNR*Vntot)
            # return max(Vf, Vfsub) > (SNR*Vntot)

        end # end trigger function
    end # end let-block

    return trigger
end # end create_antenna


# """
# Construct the boresight vectors for each antenna spaced evenly in az.
# """
# function create_boresight(orientation::EvenAzimuth, Nant, θ0)
#     ϕ = range(0, 2π, length=Nant + 1)[1:Nant]  # Evenly spaced in azimuth
#     return spherical_to_cartesian.(π / 2.0 .+ deg2rad(-θ0) * ones(Nant), ϕ, 1)
# end

# """
# Construct the boresight vectors for each antenna spaced evenly in az.
# """
# function create_boresight(orientation::EvenAzimuth, Nant, θ0)
#     ϕ = range(0, 2π, length=Nant + 1)[1:Nant]  # Evenly spaced in azimuth
#     return spherical_to_cartesian.(π / 2.0 .+ deg2rad(-θ0) * ones(Nant), ϕ, 1)
# end    

# create a model for LPDA's and for ANITA horns
"""
    LPDA(;kwargs...)

Create a model for the Low-Profile Dipole Array (LPDA) used in cosmic ray detection simulations. Configures the LPDA based on given parameters or defaults.

# Keyword Arguments
- `kwargs`: Parameters for the LPDA model, such as altitude, frequency range, and noise settings.

# Returns
- A configured LPDA antenna model.
"""
LPDA(; kwargs...) = create_antenna("dualLPDA_antennaFacsV3"; kwargs...)

"""
    ANITA(;kwargs...)

Create a model for the ANtarctic Impulsive Transient Antenna (ANITA) used in cosmic ray detection simulations. Configures the ANITA system based on given parameters or defaults.

# Keyword Arguments
- `kwargs`: Parameters for the ANITA model, including frequency range and other antenna-specific settings.

# Returns
- A configured ANITA antenna model.
"""
ANITA(; kwargs...) = create_antenna("ANITA_antennaFacs"; ν_min=150MHz, ν_max=800MHz, kwargs...)
