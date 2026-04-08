using Rotations
using StaticArrays
using Distributions
using ToggleableAsserts
using Parameters
using Unitful
import Unitful: km

const Rmoon = 1737.4km

# AOI: Area of Interest defines the location of a Regions. Used for determining is_in_region()
abstract type AreaOfInterest end

struct Quadrangle <: AreaOfInterest
    min_lat::Float64
    max_lat::Float64
    min_lon::Float64
    max_lon::Float64
end

abstract type Pole <: AreaOfInterest end
struct SouthPoleAOI <: Pole
    max_lat::Float64
end

struct NorthPoleAOI <: Pole
    min_lat::Float64
end

struct Circle <: AreaOfInterest
    center_lat::Float64
    center_lon::Float64
    radius::typeof(1.0km)
end

struct SphericalCap <: AreaOfInterest
    center_lat::Float64
    center_lon::Float64
    dlat::Float64  # Cap radius (max_lat = center_lat + dlat; min_lat = center_lat - dlat)
end

# Define region types and functions
#  Regions define contiguous locations on the Moon defined by some solid angle in lat / lon space
#  Assigning an Area of Interest within the Region will limit sampling to that fractional area
#  Note: AOIs are assumed to be isotropic throughout the region and supplied as the fraction: Area(AOI) / Area(Region)
#  This approximation is used for Permanantly Shadowed Craters (AOIs) within PolarRegions
#  If the given orbit spends similar amounts of time over the the polar region (e.g. 80S - 90S), an isotropic AOI is ok
#  Do not use AOIs where strong biases exist in spatial distribution within the region and orbit. 
#    E.g., using an AOI for PSRs relative to the WholeMoonRegion would under-count for polar orbits that spend more time at the poles where PSRs are actually found
abstract type AbstractRegion end
abstract type AOIRegion <: AbstractRegion end

@with_kw struct WholeMoonRegion <: AbstractRegion 
    aoi_frac::Float64 = 1.
end

@with_kw struct PolarRegion <: AOIRegion 
    aoi::Pole = SouthPoleAOI(-80)
    aoi_frac::Float64 = 1.
end

@with_kw struct Region <: AOIRegion
    aoi::AreaOfInterest
    aoi_frac::Float64 = 1.
end

struct CustomRegion <: AbstractRegion
    criteria::Function
    aoi_frac::Float64
    surface_area::typeof(1.0km^2)  # [km^2]
end

# PSR Definitions from the following paper giving PSR area in south and north.
# Ref: Mazarico, E., G. A. Neumann, D. E. Smith, M. T. Zuber, and M. H. Torrence (2011) “Illumination conditions of the lunar polar regions using LOLA topography.” Icarus, 211, no. 2: 1066-1081. https://doi.org/10.1016/j.icarus.2010.10.030 
# Chance of hitting PSR = area of PSR / area 10 of spherical cap 10 degrees from pole (equiprobable for C.R. to hit anywhere on pole)
#  TODO: to improve accuracy, could update this to lookup table of lat/lons in PSRs

# South pole: 5.57% of area < -80 degrees is PSR (16055 km^2)
# SouthPolePSR = (prob=0.0557, area=1.6055e4km^2) -> CustomRegion((lat, lon)->(lat <= -80) && rand() < prob, area)
SouthPolePSR = PolarRegion(SouthPoleAOI(-80), 0.0557)

# North pole: 4.46% of area > 80 degrees is PSR (12866 km^2)
# NorthPolePSR = (prob=0.0446, area=1.2866e4km^2) -> CustomRegion((lat, lon)->(lat >= 80) && rand() < prob, area)
NorthPolePSR = PolarRegion(NorthPoleAOI(80), 0.0446)

function always_true(lat, lon)
    return true
end

## 
WholeMoonMare = CustomRegion(always_true, 0.162, 6.1512e6km^2)
WholeMoonHighlands = CustomRegion(always_true, 0.8370751, 31.7522035e6km^2)

# Both poles: 5.02% of area within 10 degrees of either pole is PSR (28921 km^2)
AllPSR = CustomRegion((lat, lon)->(abs(lat) >= 80), 0.0502, 2.8921e4km^2)
# AllPSR = [SouthPolePSR(), NorthPolePSR()]  # TODO: implement list of regions

function create_region(config::String)
    config = lowercase(config)
    if config == "whole_moon"
        return WholeMoonRegion()
    elseif startswith(config, "quad:")
        # Format: "quad:min_lat,max_lat,min_lon,max_lon,aoi_frac(optional)"
        vals = parse.(Float64, split(config[6:end], ","))
        aoi = Quadrangle(vals[1], vals[2], vals[3], vals[4])
        aoi_frac = length(vals) >= 5 ? vals[5] : 1.0
        return Region(aoi, aoi_frac)
    elseif startswith(config, "cap:")
        # Format: "cap:center_lat,center_lon,dlat,aoi_frac(optional)"
        vals = parse.(Float64, split(config[5:end], ","))
        aoi = SphericalCap(vals[1], vals[2], vals[3])
        aoi_frac = length(vals) >= 4 ? vals[4] : 1.0
        return Region(aoi, aoi_frac)
    elseif startswith(config, "circle:")
        # Format: "circle:center_lat,center_lon,radius_km,aoi_frac(optional)"
        vals = parse.(Float64, split(config[8:end], ","))
        aoi = Circle(vals[1], vals[2], vals[3]*km)
        aoi_frac = length(vals) >= 4 ? vals[4] : 1.0
        return Region(aoi, aoi_frac)
    elseif startswith(config, "polar:")
        # Format: "polar:south,max_lat,aoi_frac(optional)"
        vals = split(config[7:end], ",")
        if vals[1] == "south"
            pole = SouthPoleAOI(parse(Float64,vals[2]))
        elseif vals[1] == "north"
            pole = NorthPoleAOI(parse(Float64,vals[2]))
        end
        aoi_frac = length(vals) >= 3 ? parse(Float64,vals[3]) : 1.0
        return PolarRegion(pole, aoi_frac)
    elseif config == "psr:south"
        return SouthPolePSR
    elseif config == "psr:north"
        return NorthPolePSR
    elseif config == "psr:all"
        return AllPSR
    else
        throw(ArgumentError("Invalid region configuration: $config"))
    end
end

# Region location filtering
function is_in_region(surface, region::AbstractRegion)
    if region isa WholeMoonRegion
        in_region = true
    elseif region isa CustomRegion
        lat, lon = cartesian_to_latlon(surface)
        in_region = region.criteria(lat, lon)
    else
        in_region = is_in_aoi(surface, region.aoi)
    end

    # Probabalistically account for fractional area of interest
    if region.aoi_frac < 1.0 && in_region && rand() > region.aoi_frac
        in_region = false
    end
    return in_region
end

function is_in_aoi(surface, aoi::Union{Circle, SphericalCap})
    center = latlon_to_cartesian(aoi.center_lat, aoi.center_lon)
    gc = atan(norm((center / norm(center)) × surface), 
                ((center / norm(center)) ⋅ surface))
    radius =  isa(aoi, Circle) ? aoi.radius : deg2rad(aoi.dlat) * Rmoon
    return gc * Rmoon <= radius
end

function is_in_aoi(surface, aoi::Union{NorthPoleAOI,SouthPoleAOI,Quadrangle})
    lat, lon = cartesian_to_latlon(surface)
    if aoi isa NorthPoleAOI
        return lat >= aoi.min_lat
    elseif aoi isa SouthPoleAOI
        return lat <= aoi.max_lat
    elseif aoi isa Quadrangle
        return (aoi.min_lat <= lat <= aoi.max_lat) && (aoi.min_lon <= lon <= aoi.max_lon)
    end
end

# AOIRegion sampling
function random_point_in_aoi(aoi::Pole)
    if aoi isa NorthPoleAOI
        theta_max = deg2rad(90 - aoi.min_lat)  # co-latitude
        theta_min = 0.0
    elseif aoi isa SouthPoleAOI
        theta_max = π
        theta_min = deg2rad(90 - aoi.max_lat)  # co-latitude
    end
    theta, phi = random_angles_on_cap(theta_max; theta_min=theta_min)
    return spherical_to_cartesian(theta, phi, Rmoon)
end

function random_point_in_aoi(aoi::Union{Circle,SphericalCap})
    # Calculate angle subtended by region
    dtheta = isa(aoi, Circle) ? aoi.radius / Rmoon : deg2rad(aoi.dlat)
    # Sample (dtheta) on cap at north pole
    theta, phi = random_angles_on_cap(dtheta)
    p = spherical_to_cartesian(theta, phi, Rmoon)
    # Rotate to center at (center_lat, center_lon)
    center_theta = deg2rad(90 - aoi.center_lat)
    center_phi = deg2rad(aoi.center_lon)
    center_vec = spherical_to_cartesian(center_theta, center_phi, 1.0)
    return rotate_vector_to_align(p, SA[0.0, 0.0, 1.0], center_vec)
end

function random_point_in_aoi(aoi::Quadrangle)
    # Uniform sampling in lon, and in sin(lat)
    lat1 = deg2rad(aoi.min_lat)
    lat2 = deg2rad(aoi.max_lat)
    lon1 = deg2rad(aoi.min_lon)
    lon2 = deg2rad(aoi.max_lon)
    # Uniform in sin(lat)
    s = rand()
    lat = asin(s * (sin(lat2) - sin(lat1)) + sin(lat1))
    lon = rand(Uniform(lon1, lon2))
    return spherical_to_cartesian(π/2 - lat, lon, Rmoon)
end

# Return area of region
function region_area(region::WholeMoonRegion)
    return 4 * π * Rmoon^2
end


function region_area(region::AOIRegion)
    return aoi_area(region.aoi)
end


function region_area(region::CustomRegion)
    return region.area
end

# Define spacecraft types and location sampling
abstract type Spacecraft end

struct FixedPlatform <: Spacecraft
    lat::Float64
    lon::Float64
    altitude::typeof(1.0km)
end

struct CircularOrbit <: Spacecraft 
    altitude::typeof(1.0km)
end

struct EllipticalOrbit <: Spacecraft
    periapse::typeof(1.0km)
    apoapse::typeof(1.0km)
    inclination::Float64
end

struct SampledOrbit <: Spacecraft
    latlonalt
end

function create_spacecraft(config::String)
    if startswith(config, "fixed:")
        lat, lon, alt = parse.(Float64, split(config[7:end], ","))
        return FixedPlatform(lat, lon, alt * km)
    elseif startswith(config, "orbit:")
        alt = uparse(split(config, ":")[2])
        altitude = (alt isa Unitful.AbstractQuantity) ? alt : alt * km
        return CircularOrbit(altitude)
    elseif startswith(config, "elliptical:")
        pai = split(split(config, ":")[2], ",")
        periapse = uparse(pai[1])
        apoapse = uparse(pai[2])
        inclination = parse(Float64, pai[3])
        return EllipticalOrbit(periapse, apoapse, inclination)
    elseif startswith(config, "file:")
        filename = config[6:end]
        return parse_orbit(filename)
    else
        throw(ArgumentError("Invalid spacecraft configuration"))
    end
end

"""
    parse_orbit(fname)

Parse orbital info from a CSV file with columns: time, longitude, latitude, altitude.

Returns SampledOrbit matrix of lon,lat,alt positions to randomly sample.
"""
function parse_orbit(fname="lro_orbit_1yr_2010.csv")
    data = readdlm("$(@__DIR__)/../data/$(fname)", ',', skipstart=1)
    return SampledOrbit(data[:, 2:end])
end

# Not yet implemented - no use for orbital datetimes yet
# function parse_orbit(fname="lro_orbit_1yr_2010.csv", fmt=dateformat"yyyy-mm-dd HH:MM:SS.ssssss \UTC")
#     data = readdlm("$(@__DIR__)/../data/$(fname)", ',', skipstart=1)
#     datetime = [Dates.DateTime(dt, fmt) for dt in data[:, 1]]
#     return SampledOrbit(data[:, 2:end])
# end

function get_position(spacecraft::FixedPlatform)
    return latlon_to_cartesian(spacecraft.lat, spacecraft.lon, Rmoon + spacecraft.altitude)
end

function get_position(spacecraft::CircularOrbit)
    return random_point_on_cicular_orbit(Rmoon + spacecraft.altitude)
end

function get_position(spacecraft::EllipticalOrbit)
    return random_point_on_elliptical_orbit(
        spacecraft.periapse, spacecraft.apoapse, spacecraft.inclination; rbody=Rmoon)
end

function get_position(spacecraft::SampledOrbit)
    idx = rand(1:size(spacecraft.latlonalt)[1])
    lon, lat, alt = spacecraft.latlonalt[idx, :]  # [deg, deg, km]
    return latlon_to_cartesian(lat, lon, Rmoon + alt*km)
end

function is_visible(surface, spacecraft)
    surface_normal = normalize(surface)  
    to_spacecraft = normalize(spacecraft - surface)
    return dot(surface_normal, to_spacecraft)  > 0  # central angle < 90 deg
end

"""
    random_point_on_sphere()

Return cartesian point (x, y, z) randomly distributed on sphere of radius r.

This uses the triple gaussian method which is slightly faster than the 
spherical coords method. https://mathworld.wolfram.com/SpherePointPicking.html
"""
function random_point_on_sphere(r=1)
    xyz = SVector{3, Float64}(randn(), randn(), randn())
    return r * xyz / norm(xyz)
end

"""
Return cartesian vector to spacecraft at altitude.

Note: this samples uniformly in latitude which does not give uniformly
distributed points on the sphere (overdensity of points on the poles which 
is expected for circular polar orbits).
"""
function random_point_on_cicular_orbit(r=Rmoon)
    φ = π * rand()  # Co-latitute uniform on [0, π]
    λ = 2π * rand() # Longitude uniform on [0, 2π]
    return spherical_to_cartesian(φ, λ, r)
end

"""
    random_point_on_elliptical_orbit(periapse, apoapse, inclination; rbody=Rmoon)

Return cartesian vector to spacecraft on an elliptical orbit.

Parameters:
- periapse: minimum altitude above reference body surface
- apoapse: maximum altitude above reference body surface
- inclination: orbital inclination in degrees
- rbody: radius of the central body (default: Rmoon)

Note: This samples uniformly in true anomaly which gives a non-uniform
distribution of points along the orbit (points are more densely distributed
near periapse than apoapse, which is physically correct for constant
angular momentum).
"""
function random_point_on_elliptical_orbit(periapse, apoapse, inclination; rbody=Rmoon)
    # Convert altitudes to radii from center
    rp = periapse + rbody
    ra = apoapse + rbody
    
    # Calculate orbital elements
    a = (rp + ra) / 2  # semi-major axis
    e = (ra - rp) / (ra + rp)  # eccentricity
    
    # Generate random true anomaly (θ)
    θ = 2π * rand()
    
    # Calculate radius at this true anomaly
    r = a * (1 - e^2) / (1 + e * cos(θ))
    
    # Generate random position in orbital plane
    x_orbit = r * cos(θ)
    y_orbit = r * sin(θ)
    z_orbit = 0.0km
    
    # Rotate by inclination around x-axis
    inc = deg2rad(inclination)
    x = x_orbit
    y = y_orbit * cos(inc) - z_orbit * sin(inc)
    z = y_orbit * sin(inc) + z_orbit * cos(inc)
    
    # Add random rotation around z-axis for the longitude of ascending node (Ω)
    Ω = 2π * rand()
    x_final = x * cos(Ω) - y * sin(Ω)
    y_final = x * sin(Ω) + y * cos(Ω)
    z_final = z
    
    
    return SVector{3}([x_final, y_final, z_final])
end

"""
Return cartesian vector to spacecraft at altitude.

Note: this samples uniformly in latitude which does not give uniformly
distributed points on the sphere (overdensity of points on the poles which 
is expected for circular polar orbits).
"""
function random_point_on_elliptical_orbit(r=Rmoon)
    φ = π * rand()  # Co-latitute uniform on [0, π]
    λ = 2π * rand() # Longitude uniform on [0, 2π]
    return spherical_to_cartesian(φ, λ, r)
end

"""
    random_north_pole_point()

Draw a random point on the North lunar pole
"""
random_north_pole_point() = random_point_on_cap(deg2rad(10.0); r=Rmoon)

"""
    random_south_pole_point()

Draw a random point on the South lunar pole
"""
function random_south_pole_point()
    # sample the (theta, phi) angles on the North pole
    theta, phi = random_angles_on_cap(deg2rad(10.0))

    # rotate this to be at the south pole
    theta = π - theta

    # and return them in Cartesian coordinates
    return spherical_to_cartesian(theta, phi, Rmoon)
end

"""
Generate a (theta, phi) pair on a spherical cap uniformly in solid angle.

All angles are in radians.
"""
function random_angles_on_cap(theta_max; theta_min=0.0, phi_min=0.0, phi_max=2 * pi)

    # draw a random polar angle uniformly in solid angle
    theta = acos(rand(Uniform(cos(theta_max), cos(theta_min))))

    # draw an azimuthal angle uniformaly in the desired range
    phi = rand(Uniform(phi_min, phi_max))

    return theta, phi

end

"""
Sample a random point on a spherical cap uniformly in solid angle.

We assume that the spherical cap is centered on the +z-pole of the sphere.
"""
function random_point_on_cap(theta_max; r=Rmoon, kwargs...)

    # sample the (theta, phi) angles on the cap
    theta, phi = random_angles_on_cap(theta_max; kwargs...)

    # and return them in Cartesian coordinates
    return spherical_to_cartesian(theta, phi, r)

end


"""
Convert a spherical point (theta, phi, r) into Cartesian coordinates.
"""
function spherical_to_cartesian(theta, phi, r)
    return SVector{3}(r * [sin(theta) * cos(phi),
        sin(theta) * sin(phi),
        cos(theta)])
end

"""
Convert a cartesian point (x, y, z) into a spherical point (theta, phi, r).
"""
function cartesian_to_spherical(point::SVector{3})
    x, y, z = point
    r = sqrt(x * x + y * y + z * z)
    return SVector{3}([acos(z / r), atan(y, x), r])
end

"""
Convert a cartesian point (x, y, z) into (lat, lon).
"""
function cartesian_to_latlon(point::SVector{3})
    x, y, z = point
    lat = rad2deg(asin(z / norm(point)))
    lon = rad2deg(atan(y, x))
    return lat, lon
end

"""
Helper function to convert orbital position (x, y, z) to (lat, lon, altitude).
"""
function cartesian_to_latlonalt(point::SVector{3}; rbody=Rmoon)
    lat, lon = cartesian_to_latlon(point)
    r = norm(point)
    altitude = r - rbody
    return lat, lon, altitude
end

"""
Convert a (lat, lon) to cartesian (x, y, z) vector at radius r.
"""
function latlon_to_cartesian(lat::Float64, lon::Float64)
    lat_rad = deg2rad(lat)
    lon_rad = deg2rad(lon)
    x = cos(lat_rad) * cos(lon_rad)
    y = cos(lat_rad) * sin(lon_rad)
    z = sin(lat_rad)
    return SVector{3}(x, y, z)
end

function latlon_to_cartesian(lat::Float64, lon::Float64, r::Union{Float64, typeof(1.0km)})
    return r * latlon_to_cartesian(lat, lon)
end

"""
    horizon_angle(altitude; R=Rmoon)

Calculate the angle below the horizontal of the horizon from a given altitude.

This returns "negative" angles, in radians, for "below" the horizontal.

# Arguments
- `height`: The altitude above the surface (with units).`

"""
function horizon_angle(altitude; R=Rmoon)
    return -acos(R / (R + altitude))
end

"""
    random_vector()

Generate a random Cartesian direction distributed uniformly on the unit-sphere
with the pole of the sphere aligned with `normal`
"""
function random_vector(normal=SA[0.0, 0.0, 1.0])
    @toggled_assert norm(normal) ≈ 1.0
    theta = acos(rand(Uniform(0.0, 1.0)))
    phi = rand(Uniform(0.0, 2π))
    direction = spherical_to_cartesian(theta + π / 2.0, phi, 1.0)
    return rotate_vector_to_align(direction, SA[0.0, 0.0, 1.0], normal)
end

"""
    random_direction()

Generate a random Cartesian direction distributed uniformly on the unit-sphere
with the pole of the sphere aligned with `normal` weighted by cos^2(\theta)
"""
function random_direction(normal=SA[0.0, 0.0, 1.0])
    @toggled_assert norm(normal) ≈ 1.0
    theta = acos(sqrt(rand(Uniform(0.0, 1.0))))
    phi = rand(Uniform(0.0, 2π))
    direction = spherical_to_cartesian(theta + π / 2.0, phi, 1.0)
    return rotate_vector_to_align(direction, SA[0.0, 0.0, 1.0], normal)
end

"""
    spherical_cap_area(theta, r)
Area of a spherical cap with central angle `theta` and radius `r`.
"""
function spherical_cap_area(theta, r=Rmoon)
    return 2π * r * r * (1.0 - cos(theta))
end


"""
    intersect_with_sphere(start, direction, radius)

Propagate a vector from `start` to `radius` along `direction`.

Assumes `direction` is normalized.

See the below link for the formalism behind this implementation:
    https://www.scratchapixel.com/lessons/3d-basic-rendering/
        minimal-ray-tracer-rendering-simple-shapes/ray-sphere-intersection

"""
function intersect_with_sphere(start, direction, radius)

    # the distance from the origin to the perpendicular radial line
    LD = start ⋅ direction

    # the squared distance from the origin that we intersect that line
    d2 = (start ⋅ start) - LD * LD

    # if d^2 is greater than radius^2, then we never intersect
    if d2 > radius * radius
        return nothing
    end

    # compute the offset from the central intersection
    thc = sqrt(radius * radius - d2)

    # there are two solutions. -LD +/- thc

    #  -LD + thc is the first solution
    #  -LD - thc is the second solution

    # when we are INSIDE the sphere, we want the far/second solution
    # as this is the solution along the ray path.
    # when we are OUTSIDE the sphere, we want the closest/first solution
    # as that is the location of first intersection

    # we use the radius of the starting location to determine which solution
    # we pick so this method works from inside and outside
    t = -LD - sign(norm(start) - radius) * thc

    # `t` is now the length alongth `direction` that we must step to reach the intersection
    intersection = start + t * direction

    # a sanity check that we are actually at the correct radius
    @toggled_assert (norm(intersection) - radius) / radius < 1e-6 "Radius of intersection is invalid"

    return intersection


end

"""
    propagate_and_refract(start, antenna, direction, scale, indexmodel)
Propagate from `start` along `direction` to `antenna`, refracting at the surface.
"""
function propagate_and_refract(start, direction, antenna, scale, indexmodel)

    # scale the z-coordinate to change the incident angle at the surface
    emit = direction + (scale[1] .* (start / norm(start)))
    emit /= norm(emit)

    # propagate it to the surface
    surface = intersect_with_sphere(start, emit, Rmoon)

    # we need the surface normal, not the actual surface location
    normal = surface / norm(surface)

    # calculate the incident angle at the surface - they are already normalized
    θ_i = acos(normal ⋅ emit)

    # check for TIR - send wave along the surface
    if θ_i > asin(1.0 / regolith_index(indexmodel, 0.0m))
        # return emit, surface, nothing, nothing, 0.
        θ_r = pi / 2.0
    else
        # use Snell's law to calculate the actual refracted angle
        θ_r = asin(regolith_index(indexmodel, 0.0m) * sin(θ_i))
    end

    # and calculate the actual refracted angle to the payload
    θ_true = acos(((antenna - surface) / norm(antenna - surface)) ⋅ normal)

    # and return the residuals between the two
    return emit, surface, θ_i, θ_r, abs(θ_r - θ_true)

end


function aoi_area(aoi::Quadrangle)
    # Area of a lat-lon rectangle on a sphere
    lat1 = deg2rad(aoi.min_lat)
    lat2 = deg2rad(aoi.max_lat)
    lon1 = deg2rad(aoi.min_lon)
    lon2 = deg2rad(aoi.max_lon)
    return abs(Rmoon^2 * (sin(lat2) - sin(lat1)) * (lon2 - lon1))
end

function aoi_area(aoi::Circle)
    # Area of a spherical cap with central angle theta = radius / Rmoon
    theta = aoi.radius / Rmoon
    return spherical_cap_area(theta, Rmoon)
end

function aoi_area(aoi::SphericalCap)
    # dlat is the cap radius in degrees
    theta = deg2rad(aoi.dlat)
    return spherical_cap_area(theta, Rmoon)
end

function aoi_area(aoi::NorthPoleAOI)
    # Cap from min_lat to pole
    theta = deg2rad(90 - aoi.min_lat)
    return spherical_cap_area(theta, Rmoon)
end

function aoi_area(aoi::SouthPoleAOI)
    # Cap from max_lat to south pole
    theta = deg2rad(90 + aoi.max_lat)
    return spherical_cap_area(theta, Rmoon)
end

"""
    rotate_vector_to_align(vec, from, to)

Rotate `vec` so that the direction `from` is aligned with `to`.
If `from` and `to` are already aligned, returns `vec` unchanged.
Numerically stable for nearly parallel vectors.
"""
function rotate_vector_to_align(vec::SVector{3}, from::SVector{3}, to::SVector{3})
    axis = cross(from, to)
    axis_norm = norm(axis)
    if axis_norm < 1e-8
        return SVector{3}(vec)
    end
    angle = acos(clamp(dot(from, to), -1.0, 1.0))
    rot = AngleAxis(angle, (axis / axis_norm)...)
    return SVector{3}(rot * vec)
end
