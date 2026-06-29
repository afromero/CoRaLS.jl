using Distributions
using Unitful: EeV, g, cm, C
"""
    cosmic_ray.jl

This file is part of the `CoRaLS` module and deals with modeling and analysis of cosmic rays. It includes functions to estimate key parameters of cosmic ray interactions, such as the depth of shower maximum (Xmax), and to propagate cosmic rays to their point of maximum interaction. The file also provides tools to model the charge excess profile of cosmic ray-induced air showers.

# Contents
- `estimate_Xmax`: Estimate mean and variation of Xmax.
- `propagate_to_Xmax`: Calculate the trajectory of a cosmic ray to Xmax.
- `charge_excess_profile`: Model the charge excess profile of cosmic ray showers.
"""

"""
    estimate_Xmax(Ecr::Unitful.Quantity)

Estimate the mean and variation in Xmax for a given cosmic ray energy `Ecr`.

Xmax is the maximum depth in the atmosphere where a cosmic ray-induced air shower reaches its peak intensity. This function uses the Auger fits to Xmax and σ(Xmax) from:
    https://www.icrc2019.org/uploads/1/1/9/0/
        119067782/yushkov_auger_icrc2019_mass1_talk_final.pdf

This is only valid up to 1e20.5 (VERIFY THIS).

# Arguments
- `Ecr`: Energy of the cosmic ray in EeV.

# Returns
- Tuple of (mean, sigma) of Xmax in g/cm^2.
"""
function estimate_Xmax(Ecr)

    # calculate the cosmic ray energy in log-space
    lEcr = 18.0 + log10(Ecr / 1EeV)

    # check that the range is sensible
    if lEcr > 20.8
        error("1e20.8 is above the range of `estimate_Xmax`")
    end

    # estimate the mean using a two-element linear fit
    ## PL NOTE: If we are using the Yushkov ICRC2019, then these are not the right equations
    ##          Yushkov states in their proceeding that the observed rate of change for Xmax is 
    ##          77g/cm^2 below lEcr = 18.32, and 26g/cm^2 above

    ## Old Remy fit of Xmax
    #if lEcr <= 19.0
    #    mean = 730.0g / cm^2 + ((760.0g / cm^2 - 730.0g / cm^2) / (19.0 - 18.0)) * (lEcr - 18.0)
    #else
    #    mean = 760.0g / cm^2 + ((780.0g / cm^2 - 760.0g / cm^2) / (19.6 - 19.0)) * (lEcr - 19.0)
    #end

    ## This is using P. Linton's fit of Yushkov ICRC(2019) proceeding data. 
    ## See CoRaLS elog entry 27 for presentation (https://elog.phys.hawaii.edu/elog/corals/)
    if lEcr <= 18.32
        mean = 730.0g / cm^2 + ((77g / cm^2)) * (lEcr - 18.0)
    else
        mean = 760.0g / cm^2 + ((26g / cm^2)) * (lEcr - 19.0)
    end

    ## Old Remy fit of Xmax
    # and estimate the sigma using a similar two-piece linear fit
    #if lEcr <= 18.51
    #    sigma = 59.0g / cm^2 + ((51.0g / cm^2 - 59.0g / cm^2) / (18.51 - 18.0)) * (lEcr - 18.0)
    #else
    #    sigma = 51.0g / cm^2 + ((25.0g / cm^2 - 51.0g / cm^2) / (19.71 - 18.51)) * (lEcr - 18.51)
    #end

    ## and estimate the sigma using a similar two-piece linear fit
    ## This is using P. Linton's fit of Yushkov ICRC(2019) proceeding data. 
    ## See CoRaLS elog entry 27 for presentation (https://elog.phys.hawaii.edu/elog/corals/)
    if lEcr <= 18.32
        sigma = 106.262g / cm^2 - 2.68g / cm^2 * (lEcr)
    else
        sigma = 459.423g / cm^2 - 21.958g / cm^2 * (lEcr)
    end

    # return the mean and a positive standard deviation
    return mean, max(sigma, 0.0g / cm^2)
end

"""
    propagate_to_Xmax(origin, direction, Ecr, densitymodel)

Propagate a cosmic ray from the surface to its shower maximum (Xmax).

This function calculates the point of shower maximum based on the cosmic ray's energy and direction, considering the density model of the medium through which it propagates.

# Arguments
- `origin`: Starting point of the cosmic ray.
- `direction`: Direction vector of the cosmic ray.
- `Ecr`: Energy of the cosmic ray.
- `densitymodel`: Model of the density of the medium.

# Returns
- The position vector where the cosmic ray reaches Xmax.
"""
function propagate_to_Xmax(origin, direction, Ecr, densitymodel)

    # get the mean Xmax and variation at this cosmic ray energy
    Xmax, σ = estimate_Xmax(Ecr)

    # throw for the actual grammage at shower max - assume
    # normally distributed around `Xmax` with mean `σ`
    Xint = rand(Normal(Xmax / (g / cm^2), σ / (g / cm^2)))g / cm^2

    # to do this integration in 1D, we need the polar angle
    # of the cosmic ray w.r.t the local surface
    # this assumes that the cosmic ray interacts almost immediately
    # near the surface since it is uses the local tangent plane
    # at the interaction location for all depth calculations
    # This is easy to make in 3D but is >3x as computationally expensive
    # We also assume `direction` is normalized
    ctheta = (origin / norm(origin)) ⋅ (-direction)

    # we integrate grammage here
    X = 0g / cm^2

    # our step size for this integration
    dl = 2cm

    # now propagate along `l` until we get to Xint
    for l in 0m:dl:20m

        # calculate the depth in this 1D system
        depth = l * ctheta

        # get the density at this depth
        density = regolith_density(densitymodel, depth)

        # add the integrated grammage to our total grammage
        X += dl * density

        # and check for an interaction
        if X >= Xint
            # calculate Xmax which is `l` along direction
            return origin + l * direction
        end

    end

    # if we get here, we didn't reach Xmax, this should never happen in
    # the current simulation.
    error("Xmax was not reached during propagation of the cosmic ray.")

end

"""
    charge_excess_profile(Ecr::Unitful.Quantity)

Return the profile of excess shower charge (in Coulombs) for a cosmic ray of energy `Ecr`.

This function calculates the charge profile based on parameters derived from ZHAireS simulations of cosmic ray air showers.

# Arguments
- `Ecr`: Energy of the cosmic ray in EeV.

# Returns
- Function `Q(z)` that evaluates the shower charge profile at a given distance `z`.
"""
function charge_excess_profile(Ecr::typeof(1.0EeV))

    # TODO: These need to be run at a wider variety of energies
    # these are taken from a ZHAireS shower at 1e20
    Xmax = 1733.16g / cm^2
    X0 = 1165.35g / cm^2
    λ = 87.8409g / cm^2
    ρ = 1.8g / cm^3

    # calculate the total charge in the shower
    Qmax = 0.2 * 0.62 * ((Ecr / 1GeV) |> NoUnits) * (1.602176634e-19C)

    # and return a  function that evaluates the shower profile at a given distance
    return Q(z) = Qmax * (((z .- (X0 / ρ)) ./ ((Xmax / ρ) - (X0 / ρ))) .^ (((Xmax / ρ) - (X0 / ρ)) / (λ / ρ))) .* exp.(((Xmax / ρ) .- z) ./ (λ / ρ))

end
