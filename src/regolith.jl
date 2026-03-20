using Interpolations
using Unitful
import Unitful: g, cm

"""
Abstract type for regolith index models.
"""
abstract type RegolithIndex end

"""
Constant regolith refractive index.
"""
struct ConstantIndex <: RegolithIndex end

"""
A "surface" and "deep" regolith model.

This uses a small refractive index at the surface
and a larger one for below the surface.
"""
struct SurfaceDeepIndex <: RegolithIndex end

"""
A refractive index model from Olhoeft & Strangway.

This uses Eq. 1 from Olhoeft and Strangway. (1975)
"""
struct StrangwayIndex <: RegolithIndex end


"""
Olhoeft and Strangway Curve B regolith index (1975)
"""
struct StrangwayIndexCB <: RegolithIndex end

"""
Chang'E-3 regolith index 
"""
struct CE3Index <: RegolithIndex end

"""
Diviner regolith index 
"""
struct DivinerIndex <: RegolithIndex end

"""
Combined Lunar Source Book + Diviner regolith index 
"""
struct LSB_DivinerIndex <: RegolithIndex end

"""
Lunar Source Book regolith index 
"""
struct LunarSourceBookIndex <: RegolithIndex end

"""
Chang'E-4 regolith index 
"""
struct CE4Index <: RegolithIndex end

"""
Abstract type for regolith density models.
"""
abstract type RegolithDensity end

"""
Constant regolith density
"""
struct ConstantDensity <: RegolithDensity end

"""
A regolith density model from Olhoeft & Strangway.

This is "Curve A" from O&S.
"""
struct StrangwayDensity <: RegolithDensity end

"""
A regolith density model from Olhoeft & Strangway.

This is "Curve B" from O&S.
"""
struct StrangwayDensityCB <: RegolithDensity end


struct CE3LPRDensity <: RegolithDensity end

struct DivinerRadiDensity <: RegolithDensity end

struct LunarSourceBookDensity <: RegolithDensity end

struct LSB_Diviner_Density <: RegolithDensity end

struct CE4LPRDensity_Dong2020 <: RegolithDensity end


"""
A logarithmic regolith density model from PG's ARIA MC.

This is *WRONG* - do *NOT* use this except when
trying to reproduce or evaluate results from early
August 2021.
"""
struct OldIncorrectDensity <: RegolithDensity end

"""
    regolith_density(::ConstantDensity, depth)

Calculate the density of the regolith in g/cm^3.
"""
function regolith_density(::ConstantDensity, depth)
    # a sanity check for when this called outside the regolith
    depth < 0.0m && return 0.0g / cm^3

    return 1.8g / cm^3
end

function regolith_density(::DivinerRadiDensity, depth)
    # a sanity check for when this called outside the regolith
    depth < 0.0m && return 0.0g / cm^3

    rhod = 1.8g / cm^3
    rhos = 1.1g / cm^3
    H = 5cm
    return rhod - (rhod - rhos) * exp(-depth/H)
end


function regolith_density(::LSB_Diviner_Density, depth)
    # a sanity check for when this called outside the regolith
    depth < 0.0m && return 0.0g / cm^3

    rhod = 1.92g / cm^3
    rhos = 1.1g / cm^3
    H = 6cm
    return rhod - (rhod - rhos) * exp(-depth/H)
end


function regolith_density(::LunarSourceBookDensity, depth)
    # a sanity check for when this called outside the regolith
    depth < 0.0m && return 0.0g / cm^3
    num = depth + 12.2cm
    den = depth + 18cm
    rho0 = 1.92g / cm^3
    rho = rho0 * num / den
    return rho
end

function regolith_density(::CE4LPRDensity_Dong2020, depth)
    # a sanity check for when this called outside the regolith
    depth < 0.0m && return 0.0g / cm^3
    num = 630*depth + 3250cm
    den = depth + 1480cm
    rho0 = log10(1.919)
    rho = (1/rho0) * log10(num / den)g / cm^3
    return rho
end

"""
    regolith_density(::OldIncorrectDensity, depth)

Calculate the density of the regolith in g/cm^3 for depth.
"""
function regolith_density(::OldIncorrectDensity, depth)
    # a sanity check for when this called outside the regolith
    depth < 0.0m && return 0.0g / cm^3

    # these constants stolen from Peter's CoRaLS MC
    k = 0.121g / cm^3
    rho0 = 1.27g / cm^3
    return rho0 + k * log((depth + 1cm) / cm |> Unitful.NoUnits)
end


"""
    regolith_density(::CE3LPRDensity, depth)

Calculate the density of the regolith in g/cm^3 for depth.
"""
function regolith_density(::CE3LPRDensity, depth)
    # a sanity check for when this called outside the regolith
    depth < 0.0m && return 0.0g / cm^3
    rhod = 2.25g / cm^3
    rhos = 0.85g / cm^3
    H = 103cm
    return rhod - (rhod - rhos) * exp(-depth/H)
end

"""
    regolith_density(::StrangwayDensity, depth)

Calculate the density of the regolith in g/cm^3 for depth
using the Olhoeft & Strangway curve.
"""
function regolith_density(::StrangwayDensity, depth)

    # a sanity check for when this called outside the regolith
    depth < 0.0cm && return 0.0g / cm^3

    # we define a surface density for the minimum value of the LUT
    depth <= 2.3e-16cm && return 0.80015g / cm^3

    # a sanity check for when this called outside the regolith
    depth > 23.80m && return StrangwayDensityLUT(23.80m)

    # and now just evaluate the interpolator at this depth
    return StrangwayDensityLUT(depth)

end


"""
    regolith_density(::StrangwayDensity, depth)

Calculate the density of the regolith in g/cm^3 for depth
using the Olhoeft & Strangway curve.
"""
function regolith_density(::StrangwayDensityCB, depth)

    # a sanity check for when this called outside the regolith
    depth < 0.0cm && return 0.0g / cm^3

    # we define a surface density for the minimum value of the LUT
    depth <= 2.3e-16cm && return 1.4g / cm^3

    # a sanity check for when this called outside the regolith
    depth > 23.80m && return StrangwayDensityCB_LUT(23.80m)

    # and now just evaluate the interpolator at this depth
    return StrangwayDensityCB_LUT(depth)

end

"""
    create_density_lut(N)

Create a LUT of the density as a function of depth
using "Curve A" from Olhoeft & Strangway. `N` is the
number of points in the LUT.
"""
function create_density_lut(::StrangwayDensity, N=100)

    # first we create a function to evaluate O&S parameterization
    A1 = 0.0323cm
    b1 = 4.29cm^3 / g
    A2 = 1.9e-39cm
    b2 = 60.2cm^3 / g

    # this is our O&S model for z(ρ)
    z(ρ) = -1.0cm + A1 * exp(b1 * ρ) + A2 * exp(b2 * ρ)

    # calculate the surface density for the O&S model
    # solution for z(ρ) = 0
    ρ_min = -log(A1 / cm) / b1

    # we stop at 1.61g/cm^3 which occurs at 24m below the surface
    ρ_max = 1.61g / cm^3

    # create an array of densities between the surface
    # density, and
    ρ = range(ρ_min, stop=ρ_max, length=N)

    # and evaluate the depths for each of these densities
    depth = z.(ρ)

    return LinearInterpolation(depth, ρ)
end


"""
    create_density_lut(N)

Create a LUT of the density as a function of depth
using "Curve B" from Olhoeft & Strangway. `N` is the
number of points in the LUT.
"""
function create_density_lut(::StrangwayDensityCB, N=100)

    # first we create a function to evaluate O&S parameterization
    A1 = 1.63e-5cm
    b1 = 7.87cm^3 / g
    A2 = 2.46e-28cm
    b2 = 35.7cm^3 / g

    # this is our O&S model for z(ρ)
    z(ρ) = -1.0cm + A1 * exp(b1 * ρ) + A2 * exp(b2 * ρ)

    # calculate the surface density for the O&S model
    # solution for z(ρ) = 0
    ρ_min = -log(A1 / cm) / b1

    # we stop at 1.61g/cm^3 which occurs at 24m below the surface
    ρ_max = 1.997g / cm^3

    # create an array of densities between the surface
    # density, and
    ρ = range(ρ_min, stop=ρ_max, length=N)

    # and evaluate the depths for each of these densities
    depth = z.(ρ)

    return LinearInterpolation(depth, ρ)
end




"""
    regolith_index(::ConstantIndex, depth)

A constant refractive index as a function of depth.
"""
function regolith_index(::ConstantIndex, depth)
    # a sanity check for when this called outside the regolith
    depth < 0.0m && return 1.0

    # otherwise, return a constant index
    return 1.7
end

"""
    regolith_index(::SurfaceDeepIndex, depth)

A two-part "surface" and "deep" refractive index model.
"""
function regolith_index(::SurfaceDeepIndex, depth)
    # a sanity check for when this called outside the regolith
    depth < 0.0m && return 1.0

    if depth <= 0.1cm
        return 1.27
    else
        return 1.7
    end
end


"""
    regolith_index(::IndexModel, depth, low_temp_corr_factor=0.9)

A refractive index model from Olhoeft & Strangway.

Dielectric constant fit is done by Olhoeft and Strangway
Peter estimated a 10% reduction for lunar PSR's due to the
~80 K temperatures compared to the typical lunar temperatures
of O&S, hence the low_temp_corr_factor.
"""
function regolith_index(::StrangwayIndex, depth, low_temp_corr_factor=0.9)
    # get the density at this depth - we need this in g/cm^3
    ρ = regolith_density(StrangwayDensity(), depth) / (g / cm^3)
    K = low_temp_corr_factor*(1.93^ρ)

    return sqrt(K)
end

function regolith_index(::StrangwayIndexCB, depth, low_temp_corr_factor=0.9)
    # get the density at this depth - we need this in g/cm^3
    ρ = regolith_density(StrangwayDensityCB(), depth) / (g / cm^3)
    K = low_temp_corr_factor*(1.93^ρ)

    return sqrt(K)
end

function regolith_index(::CE3Index, depth, low_temp_corr_factor=0.9)
    # get the density at this depth - we need this in g/cm^3
    ρ = regolith_density(CE3LPRDensity(), depth) / (g / cm^3)
    K = low_temp_corr_factor*(1.93^ρ)

    return sqrt(K)
end

function regolith_index(::CE4Index, depth, low_temp_corr_factor=0.9)
    # get the density at this depth - we need this in g/cm^3
    ρ = regolith_density(CE4LPRDensity_Dong2020(), depth) / (g / cm^3)
    K = low_temp_corr_factor*(1.93^ρ)

    return sqrt(K)
end

function regolith_index(::DivinerIndex, depth, low_temp_corr_factor=0.9)
    # get the density at this depth - we need this in g/cm^3
    ρ = regolith_density(DivinerRadiDensity(), depth) / (g / cm^3)
    K = low_temp_corr_factor*(1.93^ρ)

    return sqrt(K)
end

function regolith_index(::LSB_DivinerIndex, depth, low_temp_corr_factor=0.9)
    # get the density at this depth - we need this in g/cm^3
    ρ = regolith_density(LSB_Diviner_Density(), depth) / (g / cm^3)
    K = low_temp_corr_factor*(1.93^ρ)

    return sqrt(K)
end

function regolith_index(::LunarSourceBookIndex, depth, low_temp_corr_factor=0.9)
    # get the density at this depth - we need this in g/cm^3
    ρ = regolith_density(LunarSourceBookDensity(), depth) / (g / cm^3)
    K = low_temp_corr_factor*(1.93^ρ)

    return sqrt(K)
end



# create the LUT for the density as a function of depth
const StrangwayDensityLUT = create_density_lut(StrangwayDensity())
const StrangwayDensityCB_LUT = create_density_lut(StrangwayDensityCB())
