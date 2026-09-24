using Rotations
using Distributions
using CSV

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

const DEFAULT_5M_SLOPE_DISTRIBUTION = normpath(joinpath(
    @__DIR__, "..", "data", "south_polar_5m_slope_distribution.csv",
))

"""
    DataDrivenSlope([path])

Empirical surface-slope model sampled from the binned 5 m/pixel south-polar
DEM analysis stored in `data/south_polar_5m_slope_distribution.csv`.

The distribution is area-weighted over valid terrain in the 87--90°S coverage
used in that analysis.  It is not restricted to PSR pixels and does not retain
spatial correlations: each trial independently samples a polar slope angle and
uses an isotropic slope azimuth.  Within a 0.5° histogram bin, the slope angle
is sampled uniformly.
"""
struct DataDrivenSlope <: SlopeModel
    bin_edges::Vector{Float64} # radians, one more element than cdf
    cdf::Vector{Float64}
end

function DataDrivenSlope(path::AbstractString=DEFAULT_5M_SLOPE_DISTRIBUTION)
    table = CSV.File(path)
    left = Float64.(table.bin_left_deg)
    right = Float64.(table.bin_right_deg)
    weights = Float64.(table.count)

    isempty(weights) && throw(ArgumentError("Slope distribution at $path has no bins."))
    length(left) == length(right) == length(weights) ||
        throw(ArgumentError("Slope distribution at $path has inconsistent column lengths."))
    all(isfinite, left) && all(isfinite, right) && all(isfinite, weights) ||
        throw(ArgumentError("Slope distribution at $path contains a non-finite value."))
    all(right .> left) ||
        throw(ArgumentError("Slope distribution at $path has an invalid bin."))
    all(right[1:end-1] .== left[2:end]) ||
        throw(ArgumentError("Slope distribution at $path must have contiguous bins."))
    all(weight -> weight > 0.0, weights) ||
        throw(ArgumentError("Slope distribution at $path has a non-positive bin count."))

    cdf = cumsum(weights)
    cdf ./= cdf[end]
    bin_edges = deg2rad.(vcat(left, right[end]))

    return DataDrivenSlope(bin_edges, cdf)
end

const _AVAILABLE_SLOPE_MODELS = (
    "no_slope",
    "gaussian_0",
    "gaussian_7p6",
    "rayleigh_5p37",
    "data_5m",
)

"""
    available_slope_models()

Return the canonical names accepted by [`slope_model_from_name`](@ref).
"""
available_slope_models() = _AVAILABLE_SLOPE_MODELS

"""
    slope_model_from_name(name)

Construct a slope model from a stable command-line-friendly name.  In addition
to the canonical names returned by [`available_slope_models`](@ref), `none`,
`data_driven_5m`, and `ldsm_5m` are accepted aliases.
"""
function slope_model_from_name(name::AbstractString)
    normalized_name = lowercase(strip(name))

    if normalized_name in ("no_slope", "none")
        return NoSlope()
    elseif normalized_name == "gaussian_0"
        return GaussianSlope(0.0)
    elseif normalized_name in ("gaussian_7p6", "gaussian_7.6")
        return GaussianSlope(7.6)
    elseif normalized_name in ("rayleigh_5p37", "rayleigh_5.37")
        return RayleighSlope(5.37)
    elseif normalized_name in ("data_5m", "data_driven_5m", "ldsm_5m")
        return DataDrivenSlope()
    end

    options = join(available_slope_models(), ", ")
    throw(ArgumentError("Unknown slope model '$name'. Choose one of: $options."))
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
    random_surface_normal(slope::DataDrivenSlope, normal)

Generate a random normal by drawing a polar slope angle from a binned empirical
slope distribution and an isotropic azimuth.
"""
function random_surface_normal(slope::DataDrivenSlope, normal)

    # check that normal is accurately normalized
    @toggled_assert norm(normal) ≈ 1.0

    # Choose a histogram bin by its area-weighted count, then interpolate
    # uniformly within that 0.5 degree bin to avoid discretizing the angles.
    bin = searchsortedfirst(slope.cdf, rand())
    theta = slope.bin_edges[bin] + rand() * (
        slope.bin_edges[bin + 1] - slope.bin_edges[bin]
    )
    phi = rand(Uniform(0.0, 2π))
    direction = spherical_to_cartesian(theta, phi, 1.0)

    # Rotate the z-hat-frame draw so that z-hat is aligned with `normal`.
    zhat = SA[0.0, 0.0, 1.0]
    θ = acos(zhat ⋅ normal)
    axis = normal × zhat

    if norm(axis) < 1e-6
        return direction
    end

    rotated = AngleAxis(-θ, (axis / norm(axis))...) * direction
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
    #if abs(θ_i) < 0.999 * θtir
    #    return (divergence_tpar(divergencemodel, θ_i, n, args...),
    #        divergence_tperp(divergencemodel, θ_i, n, args...))
    #end

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
        if abs(θ) < θtir

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
