"""
    ice.jl

Models the surface roughness of ice, which is a critical factor in simulating the reflection and transmission of radio waves on ice surfaces. It includes abstract and concrete types representing different models of ice roughness and functions to apply these models to simulated electric fields.

## Main Components

- `IceRoughnessModel`: An abstract type for general ice roughness models.

- `NoIceRoughness`: A model representing no surface roughness (smooth ice).

- `GaussianIceRoughness`: A model representing Gaussian ice surface roughness, characterized by a standard deviation of surface angles.

- `ice_roughness`: Functions to apply roughness models to electric fields, modifying them based on the specified roughness model.

"""
abstract type IceRoughnessModel end

"""
    NoIceRoughness <: IceRoughnessModel

Represents a model of ice surface with no roughness (smooth surface). This model assumes a perfectly smooth ice surface, affecting the reflection and transmission of radio waves accordingly.
"""
struct NoIceRoughness <: IceRoughnessModel end

"""
    GaussianIceRoughness(σ::typeof(0.0cm)) <: IceRoughnessModel

A model representing Gaussian ice surface roughness, characterized by a standard deviation (σ) of the angular distribution of the surface.

# Arguments
- `σ`: Standard deviation of the surface angles, representing the degree of roughness.
"""
struct GaussianIceRoughness <: IceRoughnessModel
    σ::typeof(0.0cm)
end

"""
    ice_roughness(model::IceRoughnessModel, ν, E, θ_i[, n])

Apply the specified ice roughness model to a simulated electric field. The function modifies the electric field based on the roughness characteristics of the model.

# Arguments
- `model`: An instance of `IceRoughnessModel` to apply.
- `ν`: Frequencies of the electric field.
- `E`: Electric field amplitudes.
- `θ_i`: Angle of incidence at the ice surface.
- `n`: Optional local refractive index.  It is used by roughness models that
  depend on the wavelength in the medium, and ignored by `NoIceRoughness`.

# Returns
- Modified frequencies and electric field based on the roughness model. If the model is `NoIceRoughness`, the function returns the given frequencies and electric field. If the model is `GaussianIceRoughness`, the function modifies the electric field based on the "Diffuse reflection by rough surfaces: an introduction" by Sylvain, Pg. 671.
"""
function ice_roughness(::NoIceRoughness, ν, E, θ_i)
    return ν, E
end

# Reflected rays supply the local refractive index so that roughness models can
# evaluate their in-medium wavelength.  The smooth-surface model is a no-op,
# but it must accept the same call shape.
function ice_roughness(::NoIceRoughness, ν, E, θ_i, n)
    return ν, E
end

function ice_roughness(roughness::GaussianIceRoughness, ν, E, θ_i, n=sqrt(3.0))
    # calculate the wavelength for each frequency
    # λ = c_0 ./ (ν .* n)

    # calculate `g` - which is the parameter to our exponential model
    # g = (4π .* (roughness.σ ./ λ) .* cos(θ_i)) .^ 2.0
    factor = 4π * roughness.σ * cos(θ_i) * (n / c_0)

    # apply the roughness and return (@. does all vector operators in 1 step for efficiency)
    E_out = @. E * exp(-0.5 * (factor * ν)^2)
    return ν, E_out 
end
