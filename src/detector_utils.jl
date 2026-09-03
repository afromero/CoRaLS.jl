"""
    detector_utils.jl

This file is created have various functions that are helpful in the detector module

"""

using Unitful: MHz, K, °, km
using StaticArrays
using LinearAlgebra


"""
    antenna_temp(
        ν,
        alt,
        θ0,
        ϕ0,
        σθ,
        σϕ,
        Tsky;
        T_moon = 85u"K",
        lunar_radius = 1737.4u"km",
        nθ = 181,
        nϕ = 360,
        θ_floor = 0.0125,
        ϕ_floor = 0.08395,
        response_is_field = false,
    )

Calculate antenna temperature at every frequency in `ν`.

Inputs
------
- `ν`: array of frequencies
- `alt`: spacecraft altitude above the lunar surface
- `θ0`: antenna elevation pointing direction
- `ϕ0`: antenna azimuth pointing direction
- `σθ`: beam width in the elevation plane
- `σϕ`: beam width in the azimuth plane
- `Tsky`: array of sky temperatures, one value per frequency

Coordinate convention
---------------------
- `θ = 0`: spacecraft horizontal
- `θ = -90°`: nadir, toward the Moon
- `θ = +90°`: zenith, away from the Moon
- `ϕ`: azimuth around the spacecraft vertical

Returns
-------
An array of antenna temperatures with the same length as `ν`.
"""
function antenna_temp(
    ν::AbstractVector,
    alt,
    θ0,
    ϕ0,
    σθ,
    σϕ,
    Tsky::AbstractVector;
    T_moon = 85K,
    lunar_radius = 1737.4km,
    nθ::Integer = 181,
    nϕ::Integer = 360,
    θ_floor::Real = 0.0125,
    ϕ_floor::Real = 0.08395,
    response_is_field::Bool = true,
    nsurf = sqrt(3.0),
)

    # ------------------------------------------------------------
    # Input checks
    # ------------------------------------------------------------

    isempty(ν) &&
        throw(ArgumentError("ν cannot be empty"))

    length(ν) == length(Tsky) ||
        throw(
            ArgumentError(
                "ν and Tsky must have the same length. " *
                "Got length(ν)=$(length(ν)) and length(Tsky)=$(length(Tsky))."
            )
        )

    nθ > 0 ||
        throw(ArgumentError("nθ must be positive"))

    nϕ > 0 ||
        throw(ArgumentError("nϕ must be positive"))

    0.0 <= θ_floor <= 1.0 ||
        throw(ArgumentError("θ_floor must lie between 0 and 1"))

    0.0 <= ϕ_floor <= 1.0 ||
        throw(ArgumentError("ϕ_floor must lie between 0 and 1"))

    # ------------------------------------------------------------
    # Convert angles to dimensionless radians
    # ------------------------------------------------------------

    θ0_rad = _ant_angle_to_rad(θ0)
    ϕ0_rad = _ant_angle_to_rad(ϕ0)
    σθ_rad = _ant_angle_to_rad(σθ)
    σϕ_rad = _ant_angle_to_rad(σϕ)

    σθ_rad > 0 ||
        throw(ArgumentError("σθ must be positive"))

    σϕ_rad > 0 ||
        throw(ArgumentError("σϕ must be positive"))

    # ------------------------------------------------------------
    # Lunar geometry
    # ------------------------------------------------------------

    spacecraft_radius = lunar_radius + alt

    spacecraft_radius > lunar_radius ||
        throw(
            ArgumentError(
                "The spacecraft must be above the lunar surface."
            )
        )

    # Apparent angular radius of the Moon as seen by the spacecraft.
    α_moon = asin(
        ustrip(
            uconvert(
                Unitful.NoUnits,
                lunar_radius / spacecraft_radius,
            )
        )
    )

    cos_α_moon = cos(α_moon)

    # ------------------------------------------------------------
    # Define the antenna coordinate basis
    # ------------------------------------------------------------

    # Local spacecraft coordinates:
    #
    # x-y plane: spacecraft horizontal
    # +z: away from the Moon
    # -z: toward the Moon
    #
    # θ0 is elevation from horizontal.
    boresight = @SVector [
        cos(θ0_rad) * cos(ϕ0_rad),
        cos(θ0_rad) * sin(ϕ0_rad),
        sin(θ0_rad),
    ]

    # Unit vector corresponding to increasing elevation at boresight.
    elevation_axis = @SVector [
        -sin(θ0_rad) * cos(ϕ0_rad),
        -sin(θ0_rad) * sin(ϕ0_rad),
        cos(θ0_rad),
    ]

    # Unit vector corresponding to increasing azimuth at boresight.
    azimuth_axis = @SVector [
        -sin(ϕ0_rad),
        cos(ϕ0_rad),
        0.0,
    ]

    # Nadir points toward the center of the Moon.
    nadir = @SVector [
        0.0,
        0.0,
        -1.0,
    ]

    # ------------------------------------------------------------
    # Integrate the antenna response over the sky
    # ------------------------------------------------------------

    Δθ = π / nθ
    Δϕ = 2π / nϕ

    beam_on_moon = 0.0
    beam_on_sky = 0.0

    for iθ in 1:nθ

        # Midpoint of this elevation bin.
        θ = -π / 2 + (iθ - 0.5) * Δθ

        cosθ = cos(θ)
        sinθ = sin(θ)

        for iϕ in 1:nϕ

            # Midpoint of this azimuth bin.
            ϕ = -π + (iϕ - 0.5) * Δϕ

            direction = @SVector [
                cosθ * cos(ϕ),
                cosθ * sin(ϕ),
                sinθ,
            ]

            # Components of the sampled direction in the antenna frame.
            forward_component =
                dot(direction, boresight)

            elevation_component =
                dot(direction, elevation_axis)

            azimuth_component =
                dot(direction, azimuth_axis)

            # Signed off-boresight angles in the two principal planes.
            θev = atan(
                elevation_component,
                forward_component,
            )

            ϕev = atan(
                azimuth_component,
                forward_component,
            )

            # Antenna response in the elevation plane.
            Bθ =
                θ_floor +
                (1.0 - θ_floor) *
                exp(
                    -(θev^2) /
                    (2.0 * σθ_rad^2)
                )

            # Antenna response in the azimuth plane.
            Bϕ =
                ϕ_floor +
                (1.0 - ϕ_floor) *
                exp(
                    -(ϕev^2) /
                    (2.0 * σϕ_rad^2)
                )

            beam_response = Bθ * Bϕ

            # Antenna temperature must be weighted by power response.
            # Set response_is_field=true only if Bθ and Bϕ describe
            # electric-field amplitudes instead of power.
            if response_is_field
                beam_response = beam_response^2
            end

            # Solid-angle element for elevation/azimuth coordinates:
            #
            # dΩ = cos(θ) dθ dϕ
            # It is cosine because we are using -pi/2 ro pi/2 for elevation instead of 0 to pi
            dΩ = cosθ * Δθ * Δϕ

            weight = beam_response * dΩ

            # A direction sees the Moon when its angular separation
            # from nadir is less than the apparent lunar radius.
            sees_moon =
                dot(direction, nadir) >= cos_α_moon

            if sees_moon
                beam_on_moon += weight
            else
                beam_on_sky += weight
            end
        end
    end

    total_beam = beam_on_moon + beam_on_sky

    total_beam > 0 ||
        throw(
            ErrorException(
                "The integrated antenna response is zero."
            )
        )

    # Normalize the beam fractions.
    f_moon = beam_on_moon / total_beam
    f_sky = beam_on_sky / total_beam

    # ------------------------------------------------------------
    # Calculate antenna temperature at every frequency
    # ------------------------------------------------------------

    Tant = similar(Tsky)

    ## Calculate power fraction that will be reflected off lunar surface
    R = abs2((nsurf - 1)/(nsurf + 1))

    for i in eachindex(ν, Tsky)
        Tant[i] =
            f_moon * ((1-R)*T_moon + R*Tsky[i]) +
            f_sky * Tsky[i]
    end

    return Tant
end


"""
    ant_view(spacecraft_sph, event_sph)

Calculate the unit vector pointing from the spacecraft to an event and
express it in spacecraft-local coordinates.

Both inputs must have the form

    (r, θ, ϕ)

using the lunar-centric spherical convention:

- `r`: distance from the lunar center
- `θ`: colatitude
    - 0 at the north pole
    - π/2 at the equator
    - π at the south pole
- `ϕ`: azimuth
    - 0 on the lunar-centric +x axis
    - increasing from +x toward +y

The spacecraft-local coordinate system is:

- local +x: increasing colatitude, locally southward
- local +y: increasing azimuth, locally eastward
- local +z: radially outward from the lunar center

Returns a named tuple:

    (
        view = view_local,
        rotation = R,
        view_global = view_global_hat,
    )

The rotation matrix transforms any lunar-centric Cartesian vector into
spacecraft-local coordinates:

    vector_local = R * vector_global
"""
function ant_view(
    spacecraft_sph,
    event_sph,
)
    r_sc, θ_sc, ϕ_sc = spacecraft_sph
    r_ev, θ_ev, ϕ_ev = event_sph

    θs = _ant_angle_to_rad(θ_sc)
    ϕs = _ant_angle_to_rad(ϕ_sc)

    θe = _ant_angle_to_rad(θ_ev)
    ϕe = _ant_angle_to_rad(ϕ_ev)

    # ------------------------------------------------------------
    # Lunar-centric Cartesian positions
    # ------------------------------------------------------------

    spacecraft_cart = @SVector [
        r_sc * sin(θs) * cos(ϕs),
        r_sc * sin(θs) * sin(ϕs),
        r_sc * cos(θs),
    ]

    event_cart = @SVector [
        r_ev * sin(θe) * cos(ϕe),
        r_ev * sin(θe) * sin(ϕe),
        r_ev * cos(θe),
    ]

    # Vector pointing from spacecraft to event in lunar-centric
    # Cartesian coordinates.
    view_global = event_cart - spacecraft_cart

    view_distance = norm(view_global)

    iszero(view_distance) &&
        throw(
            ArgumentError(
                "The spacecraft and event coordinates describe the same point."
            )
        )

    view_global_hat = view_global / view_distance

    # ------------------------------------------------------------
    # Spacecraft-local orthonormal basis
    # ------------------------------------------------------------

    # Local +x: increasing colatitude, locally southward.
    x_local = @SVector [
        cos(θs) * cos(ϕs),
        cos(θs) * sin(ϕs),
        -sin(θs),
    ]

    # Local +y: increasing azimuth, locally eastward.
    y_local = @SVector [
        -sin(ϕs),
        cos(ϕs),
        0.0,
    ]

    # Local +z: radially outward from the lunar center.
    z_local = @SVector [
        sin(θs) * cos(ϕs),
        sin(θs) * sin(ϕs),
        cos(θs),
    ]

    # ------------------------------------------------------------
    # Rotation matrix
    #
    # Each row is one spacecraft-local basis vector expressed in
    # lunar-centric Cartesian coordinates.
    #
    # Therefore:
    #
    #     vector_local = R * vector_global
    # ------------------------------------------------------------

    R = @SMatrix [
        x_local[1]  x_local[2]  x_local[3]
        y_local[1]  y_local[2]  y_local[3]
        z_local[1]  z_local[2]  z_local[3]
    ]

    view_local = R * view_global_hat

    # Remove small accumulated floating-point normalization errors.
    view_local = view_local / norm(view_local)

    return (
        view = view_local,
        rotation = R,
        view_global = view_global_hat,
    )
end


"""
    project_event_pol(Ef, pol, rotation, Haxis, Vaxis)

Rotate an event polarization vector into spacecraft-local coordinates and
project the electric-field spectral density onto an antenna's H and V axes.

Returns a named tuple containing `EH`, `EV`, `Hproj`, `Vproj`, and
`pol_local`.
"""
function project_event_pol(
    Ef::AbstractVector,
    pol::AbstractVector,
    rotation::AbstractMatrix,
    Haxis::AbstractVector,
    Vaxis::AbstractVector,
)
    # Rotate the polarization vector from lunar-centric Cartesian
    # coordinates into spacecraft-local coordinates.
    pol_local = rotation * SVector{3}(pol)

    # Normalize the antenna polarization axes.
    Ĥ = normalize(SVector{3}(Haxis))
    V̂ = normalize(SVector{3}(Vaxis))

    # Signed polarization projection coefficients.
    Hproj = dot(pol_local, Ĥ)
    Vproj = dot(pol_local, V̂)

    # Project the scalar electric-field spectrum onto H and V.
    EH = Ef .* Hproj
    EV = Ef .* Vproj

    return (
        EH = EH,
        EV = EV,
        Hproj = Hproj,
        Vproj = Vproj,
        pol_local = pol_local,
    )
end



"""
    best_adjacent_triplet(V)

For every adjacent triplet in a circular antenna array, coherently sum
the three antenna spectra and then sum over frequency to obtain one
scalar voltage.

`V` must have dimensions:

    Nfreq × Nant

Returns the best scalar sum, its triplet indices, and the corresponding
frequency spectrum.
"""
function best_adjacent_triplet_value(V::AbstractMatrix)
    _, Nant = size(V)

    Nant >= 3 ||
        throw(ArgumentError("At least three antennas are required."))

    coherent_sums = Vector{eltype(V)}(undef, Nant)

    for i in 1:Nant
        i_left  = mod1(i - 1, Nant)
        i_right = mod1(i + 1, Nant)

        coherent_sums[i] =
            sum(V[:, i_left]) +
            sum(V[:, i]) +
            sum(V[:, i_right])
    end

    best_center = argmax(abs.(coherent_sums))

    return (
        value = coherent_sums[best_center],
        magnitude = abs(coherent_sums[best_center]),
        indices = (
            mod1(best_center - 1, Nant),
            best_center,
            mod1(best_center + 1, Nant),
        ),
        all_values = coherent_sums,
    )
end


"""
    best_adjacent_antennas(V, Nadj, dν)

Find the strongest coherently summed group of `Nadj` adjacent antennas
in a circular antenna array.

Inputs
------
- `V`: voltage spectral-density matrix with dimensions `Nfreq × Nant`
- `Nadj`: number of adjacent antennas to coherently sum
- `dν`: frequency-bin width

For every contiguous circular group:

1. Sum the antenna spectra frequency-by-frequency.
2. Integrate the resulting spectrum over frequency using `sum(...) * dν`.
3. Choose the group with the largest absolute integrated voltage.

For odd `Nadj`, each group has one center antenna.

For even `Nadj`, each group has two middle antennas, returned as
`middle_left` and `middle_right`.

Returns
-------
A named tuple containing:

- `value`: signed integrated voltage of the best group
- `magnitude`: absolute integrated voltage of the best group
- `indices`: antenna indices in the best group
- `spectrum`: coherently summed spectral density of the best group
- `center`: center antenna for odd `Nadj`, otherwise `nothing`
- `middle_left`: left-middle antenna for even `Nadj`, otherwise `nothing`
- `middle_right`: right-middle antenna for even `Nadj`, otherwise `nothing`
- `all_values`: integrated values for all possible groups
"""
function best_adjacent_antennas(
    V::AbstractMatrix,
    Nadj::Integer,
)
    Nfreq, Nant = size(V)

    1 <= Nadj <= Nant ||
        throw(ArgumentError(
            "Nadj must satisfy 1 ≤ Nadj ≤ Nant. " *
            "Received Nadj=$Nadj and Nant=$Nant."
        ))

    # Each column stores the summed spectrum for one candidate group.
    group_spectra = similar(V, Nfreq, Nant)

    # Integration changes voltage spectral density into voltage.
    integrated_type = typeof(zero(eltype(V)))
    group_values = Vector{integrated_type}(undef, Nant)

    # Candidate group `start` consists of:
    #
    # start, start+1, ..., start+Nadj-1
    #
    # with circular wrapping.
    for start in 1:Nant
        indices = mod1.(start .+ (0:Nadj-1), Nant)

        # Coherent frequency-by-frequency sum.
        group_spectra[:, start] .= vec(
            sum(V[:, indices], dims=2)
        )

        # All frequency bins are assumed to be in phase.
        group_values[start] =
            sum(group_spectra[:, start])
    end

    # A large negative pulse is as significant as a large positive pulse.
    best_start = argmax(abs.(group_values))

    best_indices = collect(
        mod1.(best_start .+ (0:Nadj-1), Nant)
    )

    if isodd(Nadj)
        center_position = (Nadj + 1) ÷ 2

        center = best_indices[center_position]
        middle_left = nothing
        middle_right = nothing
    else
        middle_left_position = Nadj ÷ 2
        middle_right_position = middle_left_position + 1

        center = nothing
        middle_left = best_indices[middle_left_position]
        middle_right = best_indices[middle_right_position]
    end

    return (
        value = group_values[best_start],
        magnitude = abs(group_values[best_start]),
        indices = best_indices,
        spectrum = copy(group_spectra[:, best_start]),
        center = center,
        middle_left = middle_left,
        middle_right = middle_right,
        all_values = group_values,
    )
end


"""
Convert an angle to a dimensionless value in radians.

Plain numerical values are assumed to already be radians.
"""
function _ant_angle_to_rad(
    angle::Unitful.AbstractQuantity,
)
    return ustrip(uconvert(u"rad", angle))
end

function _ant_angle_to_rad(
    angle::Real,
)
    return Float64(angle)
end