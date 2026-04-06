using NPZ
using Interpolations
using CSV
using DataFrames
using Unitful: g, cm, m, C, MeV, GeV, V, MHz, NoUnits

"""
    efield.jl

This file focuses on electric field models and calculations pertaining to cosmic ray interactions. It includes functions for loading field data, defining field and profile models, and calculating electric fields based on various models.

## Main Components
- `load_ARW`: Loads ARW E-field data and constructs an interpolator for the electric field.

- `FieldModel` and `ProfileModel`: Abstract types for electric field and profile models.

- `GaussianProfile`, `GaisserHillasProfile`: Specific profile models for charge excess in cosmic ray showers.

- `FORTE`, `JAM`, `ARW`: Specific field models, each with their own characteristics and assumptions.

- `regolith_field`: Calculates the integrated electric field at the payload for different field models.

- `attenuation_length`: Computes the attenuation length for radio waves in the regolith.

"""

"""
    load_ARW()

Load the ARW (Askaryan Radio Wave) electric field data and construct an interpolator for the electric field. This function is responsible for reading the ARW data from a file and setting up a linear interpolation model for later use in electric field calculations.

# Returns
- A function that takes an angle (θ) and frequency (ν) and returns the interpolated electric field value.
"""
function load_ARW()

    # load the data file
    data = npzread("$(@__DIR__)/../data/Efield_LUT_1e20_onesided_20210724.npz")

    # construct the interpolator
    itp = LinearInterpolation((data["view_angle_array"], data["frq_list"] ./ 1e6),
        transpose(data["RE_clean"]))

    # construct the function that does the interpolation
    return (θ, ν) -> itp(θ, (ν ./ MHz) .|> NoUnits)
end

"""
The data file for the practical and accurate calculation.
"""
const ARW_RE = load_ARW()

const DEFAULT_SATTEN = linear_interpolation(CSV.read(joinpath(@__DIR__, "../data/CE4_Cuboid_Efield_linear_perm.csv"), DataFrame)[:, 1], CSV.read(joinpath(@__DIR__, "../data/CE4_Cuboid_Efield_linear_perm.csv"), DataFrame)[:, 2])

"""
This type is the base for all instances of electric field models
for the cosmic ray emission.
"""
abstract type FieldModel end

"""
This type is the base for all longitudinal profile models.
"""
abstract type ProfileModel end

"""
Gaussian charge excess profile.
"""
struct GaussianProfile <: ProfileModel end

"""
Gaisser-Hillas charge excess profile.
"""
struct GaisserHillasProfile <: ProfileModel end

"""
The FORTE (N. Lehtinen, PRD 2004) field model.

This requires an argument to determine what charge profile
to use - currently Gaussian or GaisserHillas
"""
struct FORTE{T<:ProfileModel} <: FieldModel end

"""
The JAM field model.

This is model is *wrong* and should never be used for production
results. The angular width is *way* too wide.

Parameterization from Alvarez-muniz et al arxiv:astro-ph/0512337 2005.
"""
struct JAM <: FieldModel end

"""
The "ARW" model from Practical & Accurate Calculations...
"""
struct ARW <: FieldModel end

"""
    regolith_field(::ARW, Ecr, θ, Drego, Dvacuum; kwargs...)

Calculate the integrated electric field at the payload using the ARW model. This function takes into account the cosmic ray energy, viewing angle, and distances through regolith and vacuum to compute the electric field.

# Arguments
- `Ecr`: Energy of the cosmic ray.
- `θ`: Viewing angle.
- `Drego`: Distance through regolith.
- `Dvacuum`: Distance through vacuum.
- `kwargs`: Additional keyword arguments for fine-tuning the model.

# Returns
- Tuple of frequencies and corresponding electric field values.
"""
function regolith_field(::ARW, Ecr, θ,
    Drego, Dvacuum;
    n=sqrt(3.0), ν_min=30.0MHz, ν_max=300.0MHz,
    density=1.8g / cm^3, satten=DEFAULT_SATTEN, kwargs...)

    # TODO: Currently dν is fixed to 10 MHz by other parts of the code
    dν = 10MHz

    # calculate the array of frequencies that we calculate at
    ν = ν_min:dν:ν_max

    # calculate the actual Cherenkov angle for the desired efield
    θtrue = acos(1.0 / n)

    # calculate the Cherenkov angle for the simulation interpolants
    θsim = acos(1.0 / 2.02)

    # and now map (θ - θtrue) into the θsim space
    θactual = clamp(θsim + (θ - θtrue), 0.0, π)

    # # evaluate R|E| at the θ and across a frequency band
    # RE = ARW_RE(θactual, ν) * (1V / 1MHz)

    # # scale by the energy - the interpolator was done at 1e20 eV
    # RE *= ((Ecr / 100.0EeV) |> NoUnits)

    # # correct this simulation using the new shower profiles
    # RE *= 1.1222

    # # scale by the distance to the shower and explicitly force to
    # # normal electric-field spectral units
    # E = (RE ./ (Drego + Dvacuum)) .|> V / m / MHz

    # # get the attenuation length at this density
    # Latten = attenuation_length(ν, n, density)

    # # and apply the attenuation length due to proapgation in regolith
    # E .*= exp.(-Drego ./ Latten)

    factor = 1.1222 * Ecr / 100.0EeV / (Drego + Dvacuum)
    E = factor * ARW_RE(θactual, ν) * (1V / 1MHz)
    Latten = attenuation_length(ν, n, density; kwargs...)

    @. E *= exp(-Drego / Latten)

    if (satten !== nothing)
        @. E *= (satten(ustrip(ν)))^(Drego / m)
    end
    #@. E *= (satten.(ustrip.(ν)))^(Drego / m)



    return ν, E

end

"""
    regolith_field(::FORTE{T}, Ecr, θ, Drego, Dvacuum; kwargs...) where T <: ProfileModel

Calculate the integrated electric field at the payload using the FORTE model. This function integrates a specified longitudinal charge excess profile (Gaussian or Gaisser-Hillas) to determine the electric field.

# Arguments
- `Ecr`: Energy of the cosmic ray.
- `θ`: Viewing angle.
- `Drego`: Distance through regolith.
- `Dvacuum`: Distance through vacuum.
- `kwargs`: Additional keyword arguments for the model.

# Returns
- Tuple of frequencies and corresponding electric field values, if implemented.
"""

function regolith_field(::FORTE{GaisserHillasProfile}, Ecr, θ,
    Drego, Dvacuum;
    n=sqrt(3.0), ν_min=30.0MHz, ν_max=300.0MHz, dν=10.0MHz,
    Latten=35.0m, kwargs...)

    error("Gaisser-Hillas not implemented!")
end

function regolith_field(::FORTE{GaussianProfile}, Ecr, θ,
    Drego, Dvacuum;
    n=sqrt(3.0), ν_min=30.0MHz, ν_max=300.0MHz, dν=10.0MHz,
    density=1.8g / cm^3, kwargs...)

    # calculate the total charge in the shower for a particle with this energy
    Q = 0.2 * 0.62 * ((Ecr / 1GeV) |> NoUnits) * (1.602176634e-19C)

    # we assume the regolith has a relative permittivity of 1
    μ_r = 1.0

    # Assume a fixed shower length for now (1.32)
    # TODO: Incorporate the acutal shower length for this energy
    L = 1.32m

    # calculate the array of frequencies that we calculate at
    ν = ν_min:dν:ν_max

    # and the associated wave vectors
    k = 2π .* n .* ν / c_0

    # construct the functional form of the E-field
    RE = sqrt(2π) * (μ_0 .* μ_r .* Q .* L .* ν) .*
         sin(θ) .* exp.(-0.5(k .* L .* (cos(θ) .- (1.0 / n))) .^ 2.0)

    # scale by the distance to the shower and explicitly force to
    # normal electric-field spectral units
    E = (RE ./ (Drego + Dvacuum)) .|> V / m / MHz

    # get the attenuation length at this density
    Latten = attenuation_length(ν, n, density; kwargs...)

    # and apply the attenuation length due to proapgation in regolith
    E .*= exp.(-Drego ./ Latten)

    return ν, E

end

"""
    regolith_field(::JAM, Ecr, θ,
                        Drego, Dvacuum;
                        n = sqrt(2.8), ν_min=30.0MHz, ν_max=300.0MHz, dν=10.0MHz,
                        Latten = 35.0m)

Calculate the integrated electric field at the payload using the JAM model.

*Do not use this* - this vastly overestimates the width of the "beam",
and allows significantly more events to pass.
"""
function regolith_field(::JAM, Ecr, θ,
    Drego, Dvacuum;
    n=sqrt(3.0), ν_min=30.0MHz, ν_max=300.0MHz, dν=10.0MHz,
    density=1.8g / cm^3, kwargs...)

    # constants given in Alvarez-muniz
    X0 = 22.59g / cm^2
    RM = 11.70g / cm^2
    ρ = 1.80g / cm^3
    Ec = 40.0MeV
    kL = 23.74
    kR = 1.41
    alpha = 1.23
    beta = 2.70
    kDelta = 17.78
    kE = 3.30e-16V / cm / MHz

    # the cherenkov angle in this medium
    θc = acos(1.0 / n)

    # calculate the array of frequencies that we calculate at
    ν = ν_min:dν:ν_max

    # the width of the angular distribution of the electric field (radians)
    Δθ = ((c_0 ./ ν) * (ρ / (kDelta * X0)) .|> Unitful.NoUnits) * (1.0 / sqrt(n * n - 1.0))

    # the turnover frequency due to the longitudinal development
    # this is not currently used since we the longitudinal shower
    # development is so compact in the regolith (I think?)
    νL = (ρ / (kL * X0)) * (c_0 / abs(1.0 - n * cos(θ)))

    # this is the turnover frequency due to the radial development
    # which we consider to be dominant for regolith
    νR = (ρ / (kR * RM)) * (c_0 / sqrt(n * n - 1.0))

    # calculate the R*(electric field) at the Cherenkov angle
    REc = kE * (Ecr / Ec |> NoUnits) .* (X0 / ρ) .* (ν ./ 1MHz) .* sin(θc) .* (1.0 ./ (1.0 .+ (ν ./ νR) .^ alpha))

    # modulate this by the angular distribution
    RE = REc .* (sin(θ) ./ sin(θc)) .* exp.(-((θ - θc) ./ Δθ) .^ 2.0)

    # and divide out the distance to get the efield at the payload
    E = RE ./ (Drego + Dvacuum)

    # get the attenuation length at this density
    Latten = attenuation_length(ν, n, density; kwargs...)

    # and apply the attenuation length due to proapgation in regolith
    E .*= exp.(-Drego ./ Latten)

    # Peter's code includes the ZHSfac - let's apply it here
    # TODO: Check if this actually needed in this formalism.
    E *= sqrt(2.0) / 2.0

    # and integrate this over frequency and we are done
    return ν, E

end

"""
    attenuation_length(ν, n, density; tanδnorm=6.5e-4)

Calculate the attenuation length for radio waves in regolith at given frequencies. This function considers the refractive index and density of the regolith, as well as the loss tangent, to estimate how far radio waves can travel before being significantly attenuated.

# Arguments
- `ν`: Array of frequencies.
- `n`: Refractive index of the medium.
- `density`: Density of the medium.
- `tanδnorm`: Normalized loss tangent (default: 6.5e-4).

# Returns
- Array of attenuation lengths corresponding to each frequency.
"""
function attenuation_length(ν, n, density; tanδnorm=0.001, tand_mag=0.0005, kwargs...)

    ## Note from PL: 
    ## Lab lunar sample data at room temperature has submature highland 
    ## electric loss tangent at 1e-3 when density normalized. The Submature Highland
    ## magnetic loss tangent seems to be around 3e-3 or a bit higher at the 300-800 MHz range.

    # 1e-3 is Peter's estimate for the density-normalized loss tangent
    # for polar region temperatures. Closer to 7e-4 for PSRs.
    # here we scale it up by the density
    # assuming a density of 1.27g/cm^3 at the surface.

    ## It is unclear how losses should should scale with temperature. Using
    ## Barmatz data on Lunar simulants, highland permittivity loss goes down ~40%
    ## Magnetic highland loss didnt see any noticable change.

    ## PL NOTE: From Lunar sample data, it does not seem like magnetic losses scale much with density
    ## In low loss limit total loss of each is ~ their sum
    tanδ = tand_mag + tanδnorm * (density / (1.0g / cm^3)) |> NoUnits

    # we need the wavelength for the attenuation length calculation
    λ = c_0 ./ (n .* ν)

    # since tanδ is small, we use an approximation to calculate
    # the attenuation length. I got this from PG on Slack.
    Latten = λ ./ (π * tanδ)

    return Latten

end
