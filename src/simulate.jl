using Optim
using Distributions
using StructArrays
using LinearAlgebra
using NumericalIntegration
using Unitful
using Unitful: m, km, μV, sr
"""
This script implements various functions and types for simulating cosmic ray events and their corresponding radio frequency (RF) signals. Below is a summary of the main components of the script:

1. Enums:
   - `TrialFailed`: Enumerates possible reasons why a cosmic ray trial failed.
   - `EventGeometry`: Enumerates possible event geometries (direct or reflected).

2. Abstract Types:
   - `GeometryImplementation`: Abstract type for event simulation implementation.
   - `AbstractSignal`: Abstract type for simulated signals, both direct and reflected.

3. Concrete Types:
   - `VectorGeometry`: 3D vector geometry implementation without approximations.
   - `ScalarGeometry`: Approximate scalar geometry based on Andres' MC.

4. Signal Structs:
   - `Direct`: Information stored for a direct RF detection.
   - `Reflected`: Information stored for a reflected RF detection.

5. Functions:
   - `throw_cosmicray`: Simulates a single cosmic ray trial with a given energy, calculating the field at the payload for direct and reflected emissions.
   - `compute_direct`: Computes the 'direct' RF solution using the scalar geometry.
   - `compute_reflected`: Computes the 'reflected' RF solution using the scalar geometry.

6. Miscellaneous:
   - Overloaded `show` function for custom printing of `Direct` and `Reflected` structs.
"""

"""
An enum for possible reasons why the a cosmic ray trial failed.
"""
@enum TrialFailed TIR = 1 XmaxAfterIce = 2 NoXmax = 3 Upgoing = 4 NotVisible = 5 NotInRegion = 6

"""
An enum for possible event geometries (i.e. direct or reflected)
"""
@enum EventGeometry DirectEvent = 0 ReflectedEvent = 1

"""
Abstract type for event simulation implementation.
"""
abstract type GeometryImplementation end

"""
3D vector geometry implementation.

This is the original method of the UH MC - there are no
approximations used in this solution.
"""
struct VectorGeometry <: GeometryImplementation end

"""
*Approximate* scalar geometry based on Andres' MC.

This makes the approximation that the CR impact point is
the same as the refraction point - this is an error of
order <10m out of the 30km distance to the payload.
"""
struct ScalarGeometry <: GeometryImplementation end

"""
Abstract type for simulated signals, both direct and reflected.
"""
abstract type AbstractSignal end

"""
The information stored for a direct RF detection.
"""
mutable struct Direct <: AbstractSignal
    Ecr::typeof(0.0EeV) # the energy of the cosmic ray
    θ::Float64 # the Lunar-centric angle of the cosmic ray (deg)
    ϕ::Float64 # the Lunar-centric azimuthal angle (deg)
    θ_z::Float64 # the CR zenith angle (deg) at θ
    pol::typeof(SA[0.0, 0.0, 0.0]) # the magnitude & polarization vector at the payload
    ν_min::typeof(1.0MHz) # the minimum frequency of this e-field
    dν::typeof(1.0MHz) # the frequency bin spacing of this e-field
    ν_max::typeof(1.0MHz) # the maximum frequency of this e-field
    Ef # the total electric field spectral density at the payload
    θpol::Float64 # the polarization angle (deg) at the payload, 0 for H-Pol, 90 for V-Pol
    θ_el::Float64 # the observed elevation angle below the payload (deg)
    ψ::Float64 # the off-axis angle (deg)
    depth::typeof(0.0km) # the depth below the surface of Xmax
    Drego::typeof(0.0km) # the distance in the regolith
    Dvacuum::typeof(0.0km) # the distance in vacum
    θ_i::Float64 # the incident angle at the surface
    θ_emit::Float64 # the zenith angle of the emitted ray (deg)
    tpar::Float64 # the parallel Fresnel transmission coefficient
    tperp::Float64 # the perpendicular Fresnel transmission coefficient
    triggered::Bool # did this event trigger
    sc_altitude::typeof(0.0km) # Altitude of spacecraft for this event in Moon centric coordinates
    sc_θ::Float64 # Lunar-centric theta angle of spacecraft
    sc_ϕ::Float64 # Lunar-centric phi angle of spacecraft
end

"""
The information stored for a reflected RF detection.
"""
mutable struct Reflected <: AbstractSignal
    Ecr::typeof(0.0EeV) # the energy of the cosmic ray
    θ::Float64 # the Lunar-centric angle of the cosmic ray (deg)
    ϕ::Float64 # the Lunar-centric azimuthal angle (deg)
    θ_z::Float64 # the zenith angle (deg) at θ
    pol::typeof(SA[0.0, 0.0, 0.0]) # the magnitude & polarization vector
    ν_min::typeof(1.0MHz) # the minimum frequency of this e-field
    dν::typeof(1.0MHz) # the frequency bin spacing of this e-field
    ν_max::typeof(1.0MHz) # the maximum frequency of this e-field
    Ef # the total electric field spectral density at the payload
    polsub::typeof(SA[0.0, 0.0, 0.0]) # the polarization vector from the subsurface
    θpol::Float64 # the polarization angle (deg) at the payload, 0 for H-Pol, 90 for V-Pol
    θpolsub::Float64 # the polarization angle (deg) at the payload, 0 for H-Pol, 90 for V-Pol
    θ_el::Float64 # the observed elevation angle below the payload (deg)
    ψ::Float64 # the off-axis angle (deg)
    depth::typeof(0.0km) # the depth below the surface of Xmax
    Drego::typeof(0.0km) # the distance in the regolith
    Dvacuum::typeof(0.0km) # the distance in vacum
    θ_i::Float64 # the incident angle at the surface (deg)
    θ_emit::Float64 # the zenith angle of the emitted ray (deg)
    θ_ice::Float64 # the incident angle at the ice (deg)
    tpar::Float64 # the parallel Fresnel transmission coefficient
    tperp::Float64 # the perpendicular Fresnel transmission coefficient
    rpar::Float64 # the parallel Fresnel reflection coefficient
    rperp::Float64 # the perpendicular Fresnel reflection coefficient
    subpar::Float64 # the parallel Fresnel coefficient from the subsurface refl.
    subperp::Float64 # the perpendicular Fresnel coefficient from the subsurface refl.
    triggered::Bool # did this event trigger
    sc_altitude::typeof(0.0km) # Altitude of spacecraft for this event in Moon centric coordinates
    sc_θ::Float64 # Lunar-centric theta angle of spacecraft
    sc_ϕ::Float64 # Lunar-centric phi angle of spacecraft
end

"""
Simplifiy the printing of Signal structs.
"""
Base.show(io::IO, event::Direct) = print(io, "DirectEvent($(event.Ecr), $(event.θ_el))")
Base.show(io::IO, event::Reflected) = print(io, "ReflectedEvent($(event.Ecr), $(event.θ_el))")


"""
Simulate a single cosmic ray trial with a given energy.

This calculates the field at the payload for the direct and reflected emission,
(if they exist) and returns whether the event triggered.
"""
function throw_cosmicray(Ecr, trigger, region, spacecraft; simple_area=false, kwargs...)
    SC = get_position(spacecraft)
    if simple_area && region isa AOIRegion
        surface = random_point_in_aoi(region.aoi)
        if !is_visible(surface, SC)
            return (NotVisible, NotVisible)
        end
    elseif !simple_area && region isa AOIRegion
        surface = random_point_in_aoi(region.aoi)
        if !is_visible(surface, SC)
            return (NotVisible, NotVisible)
        elseif !is_in_region(surface, region)
            return (NotInRegion, NotInRegion)
        end
    else
        surface = random_point_on_sphere(Rmoon)
        if !is_visible(surface, SC)
            return (NotVisible, NotVisible)
        elseif !is_in_region(surface, region)
            return (NotInRegion, NotInRegion)
        end
    end
    return propagate_cosmicray(Ecr, surface, SC, trigger; kwargs...)
end


"""
Simulate a single cosmic ray trial with a given energy.

This calculates the field at the payload for the direct and reflected emission,
(if they exist) and returns whether the event triggered.
"""
function throw_cosmicray(Ecr, trigger; altitude=20.0km, kwargs...)

    # get the maximum angle that we sample CR impact points from
    θmax = -horizon_angle(altitude)

    # and assume that the polar regions extend out to 10 degrees from the poles
    θpole = deg2rad(10.0)

    # sample a point on the surface of North or South Poles
    # surface = rand() > 0.5 ? random_north_pole_point() : random_south_pole_point()

    # check if this point lies within a PSR
    # within_psr = point_impacts_psr(surface)

    # if we didn't hit a PSR, then we throw away this trial
    # within_psr == false && return NotInRegion, NotInRegion

    # throw for a random origin point on a single hemisphere
    # surface = random_point_on_cap(θmax; r=Rmoon)
    surface = random_point_on_cap(θpole; r=Rmoon)
    # surface = random_point_on_cap(π; r=Rmoon)

    # calculate the sampled value of θ for the surface point
    θ, _, _ = cartesian_to_spherical(surface)

    # check if the point lies outside the polar cap
    # ((θ > θpole) && θ < (π - θpole)) && return (NotInRegion, NotInRegion)
    θ > θpole && return (NotInRegion, NotInRegion)

    # draw a random polar angle for the spacecraft location
    # the edges of the polar cap are visible from twice its
    # opening angle away from the pole.
    # This is drawing an angle from a *circular orbit* so its
    # uniform in θ *not* uniform in cos(θ) like an area would be.
    # we assume this is symmetric around the poles so we still
    # only simulate one of the poles.
    θsc = rand(Uniform(0, θpole + θmax))

    # if the SC is not in view of the lunar pols
    # ((θsc > 2.0*θpole) && θsc < (π - 2.0*θpole)) && (NotVisible, NotVisible)
    # θsc > 2.0*θpole && (NotVisible, NotVisible)

    # and finally calculate the 3D position of the spacecraft
    # SC = spherical_to_cartesian(θsc, 0., Rmoon+altitude)
    SC = spherical_to_cartesian(θsc, rand(Uniform(0, 2π)), Rmoon + altitude)

    # calculate the total central angle between the SC and the event
    Δσ = atan(norm(SC × surface), (SC ⋅ surface))

    # if this point is not within the horizon of the SC, then we can't see it.
    # abs(θ - θsc) > θmax && return (NotVisible, NotVisible)
    abs(Δσ) > θmax && return (NotVisible, NotVisible)
    return propagate_cosmicray(Ecr, surface, SC, trigger; kwargs...)
end

"""
Propagate cosmic ray into surface and compute direct and reflected RF.
"""
function propagate_cosmicray(Ecr, surface, SC, trigger;
    ice_depth=5m,
    geometrymodel=ScalarGeometry(),
    indexmodel=StrangwayIndex(), fieldmodel=ARW(),
    divergencemodel=MixedFieldDivergence(),
    densitymodel=StrangwayDensity(),
    slopemodel=GaussianSlope(7.6),
    roughnessmodel=GaussianRoughness(2.0),
    iceroughness=GaussianIceRoughness(2.0cm),
    low_temp_corr_factor=1.0,
    kwargs...)
    # throw for a random incident cosmic ray direction in cos^2
    direction = random_direction(surface / norm(surface))

    # if it's an upgoing cosmic ray, reject the trial
    if surface ⋅ direction > 0km
        @debug "Got an upgoing cosmic ray trial! This should never happen."
        return Upgoing, Upgoing
    end

    # we now step the cosmic ray into the regolith until Xmax
    Xmax = propagate_to_Xmax(surface, direction, Ecr, densitymodel)

    # check that we didn't "leave" the regolith due to a highly
    # inclined shower being "too" long.
    if norm(Xmax) > Rmoon
        return NoXmax, NoXmax
    end

    # find the depth of the cosmic ray
    depth = norm(surface) - norm(Xmax)

    # check that Xmax occured before the ice layer
    if depth > ice_depth
        return XmaxAfterIce, XmaxAfterIce
    end

    # compute the solution for the "direct" RF
    direct = compute_direct(geometrymodel, Ecr, surface, direction, SC, Xmax;
        indexmodel=indexmodel,
        fieldmodel=fieldmodel,
        slopemodel=slopemodel,
        roughnessmodel=roughnessmodel,
        densitymodel=densitymodel,
        divergencemodel=divergencemodel,
        low_temp_corr_factor=low_temp_corr_factor,
        kwargs...)

    # compute the solution for the "reflected" RF
    reflected = compute_reflected(geometrymodel, Ecr, surface, direction, SC, Xmax, ice_depth;
        indexmodel=indexmodel,
        fieldmodel=fieldmodel,
        slopemodel=slopemodel,
        roughnessmodel=roughnessmodel,
        iceroughness=iceroughness,
        densitymodel=densitymodel,
        divergencemodel=divergencemodel,
        low_temp_corr_factor=low_temp_corr_factor,
        kwargs...)

    return direct, reflected
end

"""
Compute the 'direct' RF solution using the scalar geometry.
"""
function compute_direct(::ScalarGeometry,
    Ecr, origin, axis, antenna, Xmax;
    indexmodel=SurfaceDeepIndex(),
    densitymodel=StrangwayDensity(),
    fieldmodel=ARW(),
    divergencemodel=MixedFieldDivergence(),
    slopemodel=NoSlope(),
    roughnessmodel=NoRoughness(),
    ν_min=150MHz, ν_max=600MHz,
    low_temp_corr_factor=1.0,
    kwargs...)
    # Normal vector to the origin point on the surface, relative to flat spherical Moon
    flat_normal = origin / norm(origin)

    # throw for a random surface normal at this location relative to surface slope distribution
    surface_normal = random_surface_normal(slopemodel, flat_normal)

    # the observation vector is from the surface to the spacecraft
    obs = antenna - origin

    # the normalized view vector
    view = obs / norm(obs)

    # the exit angle at the surface to hit the space-craft
    # this is the "refracted" angle
    θ_r = acos(view ⋅ surface_normal)

    # since we have surface roughness, we can have TIR or shadowed geometries
    θ_r > π / 2.0 && return TIR

    # calculate the depth of Xmax w.r.t to the surface
    depth = norm(origin) - norm(Xmax)

    # get some needed refractive indices
    #Nsurf = regolith_index(indexmodel, 0.0m)
    #NXmax = regolith_index(indexmodel, depth)

    Nsurf = regolith_index(indexmodel, 0.0m; low_temp_corr_factor = low_temp_corr_factor)
    NXmax = regolith_index(indexmodel, depth; low_temp_corr_factor = low_temp_corr_factor)

    # we then calculate the incident angle at the surface consistent
    # with the refracted wave leaving the surface at θ_r
    θ_i = asin(sin(θ_r) / Nsurf)

    # but the distance is strongly driven by the angle below the surface
    θ_reg = asin(Nsurf * sin(θ_i) / NXmax)

    # the horizontal distance from the exit point to Xmax
    x = depth * tan(θ_reg)

    # use this to calculate the total distance in the regolith
    Drego = sqrt(x * x + depth * depth)

    # and the total distance in vacuum
    Dvacuum = norm(obs)

    # project this  onto the local horizontal - this is used to
    # shift the location of Xmax back along the local horizontal
    proj = view - (view ⋅ flat_normal) * flat_normal
    proj /= norm(proj) # make sure it is normalized

    # use Snell's law for the emission angle at Xmax (same as θ_reg)
    θ_emit = θ_reg

    # Compute the normal and proj vectors at xmax, use to calculate the emission vector
    normal_xmax = Xmax / norm(Xmax)
    proj_xmax = view - (view ⋅ normal_xmax) * normal_xmax
    proj_xmax /= norm(proj_xmax) # make sure it is normalized

    # calculate the emission vector from our two component vectors
    emit = cos(θ_emit) * normal_xmax + sin(θ_emit) * proj_xmax
    emit /= norm(emit) # make sure it is normalized

    # calculate the off-axis angle - angle between emission and axis
    ψ = acos(emit ⋅ axis)

    # calculate the electric field at the surface
    # we use the divergence formula below to account for the distance to the spacecraft
    ν, E = regolith_field(fieldmodel, Ecr, ψ, Drego, Dvacuum;
        n=NXmax,
        density=regolith_density(densitymodel, depth),
        ν_min=ν_min, ν_max=ν_max, kwargs...)

    # calculate the incident polarization vector for the direct emission
    # pol = (axis × emit) × emit
    pol = -emit × (emit × axis)
    pol /= norm(pol) # fixes some tiny rounding errors

    # the vector incident at the surface
    incident = cos(θ_i) * surface_normal + sin(θ_i) * proj
    incident /= norm(incident) # make sure it is normalized

    # incident at the surface
    niperp = incident × surface_normal
    niperp /= norm(niperp)

    # this vector defines "parallel" to the plane of incident i.e. "pokey"
    # again assuming that "origin" === "surface"
    nipar = niperp × incident
    nipar /= norm(nipar)

    # the perp. for the refracted emission
    ntperp = view × surface_normal
    ntperp /= norm(ntperp)

    # this vector defines "parallel" to the plane of incident i.e. "pokey"
    # again assuming that "origin" === "surface"
    ntpar = ntperp × view
    ntpar /= norm(ntpar)

    # @assert θ_r ≈ asin(Nsurf * sin(θ_i)) "$(θ_r) ≈ $(asin(Nsurf * sin(θ_i)))"

    # calculate the modified Fresnel coefficients for transmission
    # tpar = divergence_tpar(divergencemodel, θ_i, Nsurf, Drego, Dvacuum)
    # tperp = divergence_tperp(divergencemodel, θ_i, Nsurf, Drego, Dvacuum)
    tpar, tperp = surface_transmission(roughnessmodel, divergencemodel,
        θ_i, Nsurf, Drego, Dvacuum)

    # project `pol` onto `npar` and `nperp`, apply Fresnel,
    # and then recombine into the transmitted polarization vector
    # poltr = tperp*(pol ⋅ nperp)*nperp + tpar*(pol ⋅ npar)*npar
    poltr = tperp * (pol ⋅ niperp) * ntperp + tpar * (pol ⋅ nipar) * ntpar

    θpol = atan(poltr ⋅ ntpar, poltr ⋅ ntperp)

    # calculate the zenith angle of the cosmic ray
    # flip the axis so we get the complement of the angle
    zenith = acos(flat_normal ⋅ (-axis))

    # get the Lunar-centric angle of the cosmic ray impact point
    θ, ϕ, _ = cartesian_to_spherical(origin)

    # calculate the below horizon angle at the payload
    el = -(acos(-view ⋅ (antenna / norm(antenna))) - pi / 2.0)
    sc_θ, sc_ϕ, sc_altitude = cartesian_to_spherical(antenna)

    # and construct and return the signal
    return Direct(Ecr, rad2deg(θ), rad2deg(ϕ), rad2deg(zenith),
        poltr, ν_min, 10MHz, ν_max, E .|> (μV / m / MHz), rad2deg(θpol), rad2deg(el), rad2deg(ψ),
        depth, Drego, Dvacuum, rad2deg(θ_i), rad2deg(θ_emit), tpar, tperp, false, , sc_altitude, sc_θ, sc_ϕ)

end


"""
Compute the 'reflected' RF solution using the scalar solution.
"""
function compute_reflected(::ScalarGeometry,
    Ecr, origin, axis, antenna, Xmax, ice_depth;
    Nice=1.305,
    Nbed=1.6, # 1.6 for regolith, 2.8 for bedrock
    ice_thickness=1.0m,
    indexmodel=SurfaceDeepIndex(),
    densitymodel=StrangwayDensity(),
    fieldmodel=ARW(),
    divergencemodel=MixedFieldDivergence(),
    slopemodel=NoSlope(),
    roughnessmodel=NoRoughness(),
    iceroughness=NoIceRoughness(),
    low_temp_corr_factor=1.0,
    ν_min=150MHz, ν_max=600MHz,
    kwargs...)
    # Normal vector to the origin point on the surface, relative to flat spherical Moon
    flat_normal = origin / norm(origin)

    # throw for a random surface normal at this location relative to surface slope distribution
    surface_normal = random_surface_normal(slopemodel, flat_normal)

    # the observation vector is from the surface to the spacecraft
    obs = antenna - origin

    # the normalized view vector
    view = obs / norm(obs)

    # the exit angle at the surface to hit the space-craft
    # this is the "refracted" angle
    θ_r = acos(view ⋅ surface_normal)

    # check that with a random slope, we aren't beyond TIR
    θ_r > π / 2.0 && return TIR

    # calculate the depth of Xmax w.r.t to the surface
    depth = norm(origin) - norm(Xmax)

    # we need the refractive index at the surface and at the layer several times
    #Nsurf = regolith_index(indexmodel, 0.0m)
    #NXmax = regolith_index(indexmodel, depth)
    #Nrego_at_ice = regolith_index(indexmodel, ice_depth)

    Nsurf = regolith_index(indexmodel, 0.0m; low_temp_corr_factor=low_temp_corr_factor)
    NXmax = regolith_index(indexmodel, depth; low_temp_corr_factor=low_temp_corr_factor)
    Nrego_at_ice = regolith_index(indexmodel, ice_depth; low_temp_corr_factor=low_temp_corr_factor)
    
    # we then calculate the incident angle at the surface consistent
    # with the refracted wave leaving the surface at θ_r
    θ_i = asin(sin(θ_r) / Nsurf)

    # but the distance is strongly driven by the angle below the surface
    θ_reg = asin(Nsurf * sin(θ_i) / NXmax)

    # we have to split the calculation of the distance into two parts
    # pre- and post- reflection
    x1 = (ice_depth - depth) * tan(θ_reg) # pre-reflected
    x2 = ice_depth * tan(θ_reg) # post-reflection

    # check thate we aren't violating TIR
    (x1 + x2) > (2.0 * ice_depth - depth) / sqrt(regolith_index(indexmodel, 0.0m)^2.0 - 1) && return TIR

    # the total distance is calculated over each step
    Drego = sqrt(x1 * x1 + (ice_depth - depth) * (ice_depth - depth)) +
            sqrt(x2 * x2 + ice_depth * ice_depth)

    # and the total distance in vacuum
    Dvacuum = norm(obs)

    # project the viw  onto the local horizontal - this is used to
    # shift the location of Xmax back along the local horizontal
    # this defines the horizontal unit-vector in our local coordinate system
    proj = view - (view ⋅ flat_normal) * flat_normal
    proj /= norm(proj) # make sure it is normalized

    # calculate the approximate location of Xmax in our system
    # start at the surface - go "down" by `depth` and then "back"
    # along `flat`
    source = origin - (2.0 * ice_depth - depth) * flat_normal - (x1 + x2) * proj

    # use Snell's law for the emission angle at Xmax (same as θ_reg)
    θ_emit = θ_reg

    # Compute the normal and proj vectors at xmax, use to calculate the emission vector
    normal_xmax = source / norm(source)
    proj_xmax = view - (view ⋅ normal_xmax) * normal_xmax
    proj_xmax /= norm(proj_xmax) # make sure it is normalized

    # calculate the emission vector from our two component vectors
    emit = -cos(θ_emit) * normal_xmax + sin(θ_emit) * proj_xmax
    emit /= norm(emit) # make sure it is normalized

    # calculate the off-axis angle - angle between emission and axis
    ψ = acos(emit ⋅ axis)

    # calculate the electric field at the spacecraft
    ν, E = regolith_field(fieldmodel, Ecr, ψ, Drego, Dvacuum;
        n=NXmax,
        density=regolith_density(densitymodel, depth),
        ν_min=ν_min, ν_max=ν_max, kwargs...)


    # calculate the incident polarization vector for the emission
    pol = -emit × (emit × axis)
    pol /= norm(pol) # fixes some tiny rounding errors

    # calculate the angle at the ice layer using Snell's law
    θ_ice = asin(Nsurf * sin(θ_i) / Nrego_at_ice)

    # the vector incident at the surface
    incident = cos(θ_i) * surface_normal + sin(θ_i) * proj
    incident /= norm(incident) # make sure it is normalized

    # incident at the surface
    # TODO: which of these nipar and niperp is correct?  same for direct or no?
    #niperp = emit × surface_normal
    niperp = incident × surface_normal
    niperp /= norm(niperp)

    # this vector defines "parallel" to the plane of incident i.e. "pokey"
    # again assuming that "origin" === "surface"
    #nipar = niperp × emit
    nipar = niperp × incident
    nipar /= norm(nipar)

    # the perp. for the refracted emission
    ntperp = view × surface_normal
    ntperp /= norm(ntperp)

    # this vector defines "parallel" to the plane of incident i.e. "pokey"
    # again assuming that "origin" === "surface"
    ntpar = ntperp × view
    ntpar /= norm(ntpar)

    # apply roughness to the simulated electric field
    ν, E = ice_roughness(iceroughness, ν, E, θ_ice, NXmax)

    # calculate the Fresnel reflection coefficients at the ice
    rpar, rperp = fresnel_coeffs(θ_ice, Nrego_at_ice, Nice)[1:2]

    # calculate the modified Fresnel coefficients for transmission
    # use the refractive index at the surface - this handles the random
    # properties of surface roughness as well
    tpar, tperp = surface_transmission(roughnessmodel, divergencemodel,
        θ_i, Nsurf, Drego, Dvacuum)

    # project `pol` onto `npar` and `nperp`, apply Fresnel,
    # and then recombine into the transmitted polarization vector
    poltr = rperp * tperp * (pol ⋅ niperp) * ntperp + rpar * tpar * (pol ⋅ nipar) * ntpar

    # and lastly calculate the polarization angle above the surface
    θpol = atan(poltr ⋅ ntpar, poltr ⋅ ntperp)

    # construct the parallel and perpendicular coefficients for the
    # reflection of the *bottom* of the subsurface layer
    if ice_thickness > 0.0m

        # get the transmission from the regolith into the ice
        sub_tpar, sub_tperp = fresnel_coeffs(θ_ice, Nrego_at_ice, Nice)[3:4]

        # get the refracted angle in the ice layer using Snell's law
        θ_bed = asin(Nrego_at_ice * sin(θ_ice) / Nice)

        # get the reflection coefficient from at the ice->regolith
        # or ice->bedrock interface
        bed_rpar, bed_rperp = fresnel_coeffs(θ_bed, Nice, Nbed)[1:2]

        # and then the transmission coefficient from the ice back
        # into the top of the regolith
        ice_tpar, ice_tperp = fresnel_coeffs(θ_bed, Nice, Nrego_at_ice)[3:4]

        # and finally stack them all together
        # # `sub` for "sub"surface layer.
        subpar = sub_tpar * bed_rpar * ice_tpar
        subperp = sub_tperp * bed_rperp * ice_tperp

        # construct the new polarization vector above the surface
        polsub = subperp * tperp * (pol ⋅ niperp) * ntperp + subpar * tpar * (pol ⋅ nipar) * ntpar

        # and lastly calculate the polarization angle above the surface for this reflection
        θpolsub = atan(polsub ⋅ ntpar, polsub ⋅ ntperp)

    else
        # otherwise, we don't treat the subsurface reflection
        Efsub = 0.0μV / m
        polsub = SA[0.0, 0.0, 0.0]
        θpolsub = 0.0
        subpar = subperp = 0.0
    end

    # calculate the zenith angle of the cosmic ray
    # flip the axis so we get the complement of the angle
    zenith = acos(flat_normal ⋅ (-axis))

    # get the Lunar-centric angle of the cosmic ray impact point
    θ, ϕ, _ = cartesian_to_spherical(origin)

    # calculate the below horizon angle at the payload
    el = -(acos(-view ⋅ (antenna / norm(antenna))) - pi / 2.0)
    sc_θ, sc_ϕ, sc_altitude = cartesian_to_spherical(antenna)

    # and construct and return the signal
    return Reflected(Ecr, rad2deg(θ), rad2deg(ϕ), rad2deg(zenith),
        poltr,
        ν_min, 10MHz, ν_max,
        E .|> (μV / m / MHz), polsub, rad2deg(θpol), rad2deg(θpolsub), rad2deg(el), rad2deg(ψ),
        depth, Drego, Dvacuum, rad2deg(θ_i), rad2deg(θ_emit), rad2deg(θ_ice),
        tpar, tperp, rpar, rperp, subpar, subperp, false, sc_altitude, sc_θ, sc_ϕ)

end

