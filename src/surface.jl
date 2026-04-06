using Rotations
using Distributions

"""
Abstract type for regolith surface roughness model.
"""
abstract type RoughnessModel end

"""
No surface roughness (smooth)
"""
struct NoRoughness <: RoughnessModel end

"""
Gaussian surface roughness model.

This is characterized by Gaussian "sigma" of the
angular distribution of the surface.
"""
struct GaussianRoughness <: RoughnessModel
    σ::Float64 # in radians
    GaussianRoughness(σ) = new(deg2rad(σ))
end

"""
Abstract type for regolith surface slope model.
"""
abstract type SlopeModel end

"""
No surface slope (smooth)
"""
struct NoSlope <: SlopeModel end

"""
Gaussian surface slope model.

This is characterized by Gaussian "sigma" of the
angular distribution of the surface in degrees.
"""
struct GaussianSlope <: SlopeModel
    σ::Float64 # in radians
    GaussianSlope(σ) = new(deg2rad(σ))
end


"""
Rayleigh surface slope model.

"""
struct RayleighSlope <: SlopeModel
    σ::Float64 # in radians
    RayleighSlope(σ) = new(deg2rad(σ))
end

"""
    random_surface(slope::NoSlope, normal)

Generate a random normal vector for a surface slope that would
originally be pointing along `normal` in the absence of roughness.

This just returns `normal` since this implements a "No Slope" model.
"""
function random_surface_normal(slope::NoSlope, normal)
    return normal
end

"""
    random_surface(slope::GaussianSlope, normal)

Generate a random normal vector for a surface slope that would
originally be pointing along `normal` in the absence of roughness.
"""
function random_surface_normal(slope::GaussianSlope, normal)

    # check that normal is accurately normalized
    @toggled_assert norm(normal) ≈ 1.0

    # draw a random polar angle consistent with the slope
    theta = abs(rand(Normal(0.0, slope.σ))) # in radians

    # draw an azimuthal angle uniformaly in the desired range
    phi = rand(Uniform(0.0, 2π))

    # this is the vector in the z-hat space
    direction = spherical_to_cartesian(theta, phi, 1.0)

    # we must now rotate this so that z-hat is aligned with normal
    # we do this with an axis angle representation

    # this is our z-hat vector in the original coordinate system
    zhat = SA[0.0, 0.0, 1.0]

    # construct the angle between z-hat and `normal` - already normalized
    θ = acos(zhat ⋅ normal)

    # and construct the *right-handed* axis
    axis = normal × zhat

    # have to check that axis can be properly normalized
    if norm(axis) < 1e-6
        # if this is true, we are in the same coordinate system
        return direction
    end

    # otherwide, normalize, do the rotation
    rotated = AngleAxis(-θ, (axis / norm(axis))...) * direction
    # rotated /= norm(rotated)
    # fix some precision issues

    # check that this vector is always less than 6-sigma
    # @toggled_assert acos(rotated ⋅ normal) < (6.0*slope.σ)
    @toggled_assert norm(rotated) ≈ 1.0

    return rotated
end



"""
    random_surface(slope::RayleighSlope, normal)

Generate a random normal vector for a surface slope that would
originally be pointing along `normal` in the absence of roughness.
"""
function random_surface_normal(slope::RayleighSlope, normal)

    # check that normal is accurately normalized
    @toggled_assert norm(normal) ≈ 1.0

    # draw a random polar angle consistent with the slope
    theta = abs(rand(Rayleigh(slope.σ))) # in radians

    # draw an azimuthal angle uniformaly in the desired range
    phi = rand(Uniform(0.0, 2π))

    # this is the vector in the z-hat space
    direction = spherical_to_cartesian(theta, phi, 1.0)

    # we must now rotate this so that z-hat is aligned with normal
    # we do this with an axis angle representation

    # this is our z-hat vector in the original coordinate system
    zhat = SA[0.0, 0.0, 1.0]

    # construct the angle between z-hat and `normal` - already normalized
    θ = acos(zhat ⋅ normal)

    # and construct the *right-handed* axis
    axis = normal × zhat

    # have to check that axis can be properly normalized
    if norm(axis) < 1e-6
        # if this is true, we are in the same coordinate system
        return direction
    end

    # otherwide, normalize, do the rotation
    rotated = AngleAxis(-θ, (axis / norm(axis))...) * direction
    # rotated /= norm(rotated)
    # fix some precision issues

    # check that this vector is always less than 6-sigma
    # @toggled_assert acos(rotated ⋅ normal) < (6.0*slope.σ)
    @toggled_assert norm(rotated) ≈ 1.0

    return rotated
end

"""
    surface_transmission(::NoRoughness, ν, E, θ_i)

Apply roughness to a simulated electric field given the frequencies,
the electric field, and the angle of incidence at the surface. This
returns the frequencies and electric field back to the caller.

Without roughness, we just return ν and E as given.
"""
function surface_transmission(roughness::NoRoughness, divergencemodel, θ_i, n, args...)

    # return zero if we are beyond TIR
    θ_i > asin(1.0 / n) && return 0, 0

    # get the fresnel coefficients
    tpar = divergence_tpar(divergencemodel, θ_i, n, args...)
    tperp = divergence_tperp(divergencemodel, θ_i, n, args...)

    return tpar, tperp

end

"""
    surface_transmission(::GaussianRoughness, ν, E, θ_i)

Apply a simple model for surface transmission through a rough
surface by taking the average of the transmission coefficient
over the Gaussian

This is based on the "Diffuse reflection by rough surfaces: an introduction"
by Sylvain, Pg. 671
"""
function surface_transmission(roughness::GaussianRoughness, divergencemodel, θ_i, n, args...)

    # we want to ignore any trials that are outside of TIR
    θtir = asin(1.0 / n) # we are going into vacuum

    # only do this near the horizon
    # if abs(θ_i - θtir) > 2.0 * roughness.σ
    if θ_i < 0.999 * θtir
        return (divergence_tpar(divergencemodel, θ_i, n, args...),
            divergence_tperp(divergencemodel, θ_i, n, args...))
    end

    # the number of samples that we throw
    N = 50

    # generate the N random samples from this Gaussian around θ_i
    # σ, in roughness, is *already* in radians.
    Θ = rand(Normal(θ_i, roughness.σ), N)

    # this is the average transmission coefficient that we build
    # one for each polarization
    Tpar = 0.0
    Tperp = 0.0

    # loop over each incident angle
    for θ in Θ
        # if we haven't TIR'd
        if θ < θtir

            # and calculate the coefficients
            Tpar += divergence_tpar(divergencemodel, θ, n, args...)
            Tperp += divergence_tperp(divergencemodel, θ, n, args...)
        end
    end

    # and convert it to an average
    Tpar /= N
    Tperp /= N

    # and return the pair of coefficients
    return Tpar, Tperp

end
