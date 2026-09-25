using Rotations
using Distributions
using CSV
using NPZ

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

const DEFAULT_80M_SOUTH_POLAR_SLOPE_MAP = normpath(joinpath(
    @__DIR__, "..", "data", "ldsm_87s_slope_80mpp.npz",
))

const DEFAULT_40M_SOUTH_POLAR_SLOPE_MAP = normpath(joinpath(
    @__DIR__, "..", "data", "ldsm_87s_slope_40mpp.npz",
))

const DEFAULT_5M_SOUTH87_SLOPE_SHARDS = normpath(joinpath(
    @__DIR__, "..", "data", "ldsm_87s_slope_5mpp_1deg_south87_shards",
))

const DEFAULT_5M_SOUTH87_DEM_NORMAL_SHARDS = normpath(joinpath(
    @__DIR__, "..", "data", "ldem_87s_5mpp_normal_1deg_2deg_south87_shards",
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

"""
    SouthPolarSlopeMap([path]; min_latitude_deg=-87.0,
                       fallback=GaussianSlope(7.6))

Location-aware slope model backed by the native 80 m/pixel south-polar LDSM
slope-angle map. The map stores only slope magnitude, so a finite map value
sets the polar tilt and the tilt azimuth is drawn isotropically. The map is
used from `min_latitude_deg` through the south pole; invalid map pixels and
locations outside that range use `fallback`.

The NPZ contains a `Float32` 2500 x 2500 slope-angle array (about 25 MB when
loaded) plus the polar-stereographic raster transform. It is loaded once when
this model is constructed, not on every trial.
"""
struct SouthPolarSlopeMap{T <: AbstractMatrix{<:Real}, F <: SlopeModel} <: SlopeModel
    slope_deg::T
    pixel_size_m::Float64
    origin_x_m::Float64
    origin_y_m::Float64
    radius_m::Float64
    min_latitude_deg::Float64
    fallback::F
end

function SouthPolarSlopeMap(path::AbstractString=DEFAULT_80M_SOUTH_POLAR_SLOPE_MAP;
    min_latitude_deg=-87.0,
    fallback::SlopeModel=GaussianSlope(7.6),
)
    isfile(path) || throw(ArgumentError(
        "South-polar slope map $path is missing. Build local terrain products using " *
        "scripts/terrain/README.md."
    ))
    table = npzread(path)
    required = ("slope_deg", "pixel_size_m", "origin_x_m", "origin_y_m", "radius_m")
    missing = filter(key -> !haskey(table, key), required)
    isempty(missing) || throw(ArgumentError(
        "South-polar slope map at $path is missing: $(join(missing, ", "))"
    ))

    slope_deg = table["slope_deg"]
    ndims(slope_deg) == 2 || throw(ArgumentError("slope_deg at $path must be two-dimensional."))
    pixel_size_m = Float64(only(table["pixel_size_m"]))
    origin_x_m = Float64(only(table["origin_x_m"]))
    origin_y_m = Float64(only(table["origin_y_m"]))
    radius_m = Float64(only(table["radius_m"]))
    pixel_size_m > 0.0 || throw(ArgumentError("pixel_size_m at $path must be positive."))
    radius_m > 0.0 || throw(ArgumentError("radius_m at $path must be positive."))
    -90.0 <= min_latitude_deg <= 90.0 ||
        throw(ArgumentError("min_latitude_deg must be between -90 and 90 degrees."))

    return SouthPolarSlopeMap(
        slope_deg, pixel_size_m, origin_x_m, origin_y_m, radius_m,
        Float64(min_latitude_deg), fallback,
    )
end

"""
    SouthPolarSlopeMap40m([path]; kwargs...)

Load the native 40 m/pixel south-polar LDSM slope-angle map. This has the
same location lookup and randomized-azimuth treatment as
[`SouthPolarSlopeMap`](@ref), whose default map is the 80 m/pixel product.
"""
function SouthPolarSlopeMap40m(
    path::AbstractString=DEFAULT_40M_SOUTH_POLAR_SLOPE_MAP; kwargs...
)
    return SouthPolarSlopeMap(path; kwargs...)
end

# Location-aware slope model backed by 1-degree-quantized UInt8 LDSM shards.
# Only a bounded least-recently-used cache of shards is held in memory.
mutable struct ShardedSouthPolarSlopeMap{F <: SlopeModel} <: SlopeModel
    directory::String
    shape::Tuple{Int, Int}
    shard_size_px::Int
    pixel_size_m::Float64
    origin_x_m::Float64
    origin_y_m::Float64
    radius_m::Float64
    min_latitude_deg::Float64
    slope_step_deg::Float64
    nodata_code::UInt8
    fallback::F
    cache_limit::Int
    cache::Dict{Tuple{Int, Int}, Matrix{UInt8}}
    cache_recency::Vector{Tuple{Int, Int}}
    cache_lock::ReentrantLock
end

function ShardedSouthPolarSlopeMap(
    directory::AbstractString=DEFAULT_5M_SOUTH87_SLOPE_SHARDS;
    fallback::SlopeModel=GaussianSlope(7.6),
    cache_shards::Integer=16,
)
    cache_shards > 0 || throw(ArgumentError("cache_shards must be positive."))
    manifest_path = joinpath(directory, "manifest.npz")
    isfile(manifest_path) || throw(ArgumentError(
        "South-polar slope shards are missing at $directory. Build local terrain products using " *
        "scripts/terrain/README.md."
    ))
    table = npzread(manifest_path)
    required = (
        "shape", "shard_size_px", "pixel_size_m", "origin_x_m", "origin_y_m",
        "radius_m", "min_latitude_deg", "slope_step_deg", "nodata_code",
    )
    missing = filter(key -> !haskey(table, key), required)
    isempty(missing) || throw(ArgumentError(
        "Slope-shard manifest at $manifest_path is missing: $(join(missing, ", "))"
    ))
    shape_values = Int.(table["shape"])
    length(shape_values) == 2 || throw(ArgumentError("manifest shape must have two entries."))
    shape = (shape_values[1], shape_values[2])
    all(>(0), shape) || throw(ArgumentError("manifest shape entries must be positive."))

    return ShardedSouthPolarSlopeMap(
        String(directory), shape, Int(only(table["shard_size_px"])),
        Float64(only(table["pixel_size_m"])), Float64(only(table["origin_x_m"])),
        Float64(only(table["origin_y_m"])), Float64(only(table["radius_m"])),
        Float64(only(table["min_latitude_deg"])), Float64(only(table["slope_step_deg"])),
        UInt8(only(table["nodata_code"])), fallback, Int(cache_shards),
        Dict{Tuple{Int, Int}, Matrix{UInt8}}(), Tuple{Int, Int}[], ReentrantLock(),
    )
end

function _shard_filename(map::ShardedSouthPolarSlopeMap, shard_row::Int, shard_col::Int)
    row_name = lpad(string(shard_row), 3, '0')
    col_name = lpad(string(shard_col), 3, '0')
    return joinpath(map.directory, "shard_r$(row_name)_c$(col_name).npz")
end

function _cached_slope_shard(
    map::ShardedSouthPolarSlopeMap, shard_row::Int, shard_col::Int
)
    key = (shard_row, shard_col)
    lock(map.cache_lock)
    try
        if haskey(map.cache, key)
            index = findfirst(==(key), map.cache_recency)
            deleteat!(map.cache_recency, index)
            push!(map.cache_recency, key)
            return map.cache[key]
        end

        filename = _shard_filename(map, shard_row, shard_col)
        isfile(filename) || throw(ArgumentError("Missing slope shard $filename"))
        codes = npzread(filename)["slope_code"]
        codes isa Matrix{UInt8} || throw(ArgumentError(
            "Slope shard $filename must contain a UInt8 slope_code matrix."
        ))
        if length(map.cache_recency) == map.cache_limit
            evicted = popfirst!(map.cache_recency)
            delete!(map.cache, evicted)
        end
        map.cache[key] = codes
        push!(map.cache_recency, key)
        return codes
    finally
        unlock(map.cache_lock)
    end
end

function slope_degrees_at(map::ShardedSouthPolarSlopeMap, normal)
    @toggled_assert norm(normal) ≈ 1.0
    lat, lon = cartesian_to_latlon(normal)
    lat <= map.min_latitude_deg || return nothing
    lat_rad = deg2rad(lat)
    lon_rad = deg2rad(mod(lon, 360.0))
    ρ = 2.0 * map.radius_m * tan((π / 2.0 + lat_rad) / 2.0)
    isfinite(ρ) || return nothing
    x = ρ * sin(lon_rad)
    y = -ρ * cos(lon_rad)
    isfinite(x) && isfinite(y) || return nothing

    col = floor(Int, (x - map.origin_x_m) / map.pixel_size_m) + 1
    row = floor(Int, (map.origin_y_m - y) / map.pixel_size_m) + 1
    1 <= row <= map.shape[1] && 1 <= col <= map.shape[2] || return nothing

    shard_row = (row - 1) ÷ map.shard_size_px
    shard_col = (col - 1) ÷ map.shard_size_px
    codes = _cached_slope_shard(map, shard_row, shard_col)
    local_row = (row - 1) % map.shard_size_px + 1
    local_col = (col - 1) % map.shard_size_px + 1
    checkbounds(Bool, codes, local_row, local_col) || return nothing
    code = codes[local_row, local_col]
    return code == map.nodata_code ? nothing : Float64(code) * map.slope_step_deg
end

function random_surface_normal(slope::ShardedSouthPolarSlopeMap, normal)
    @toggled_assert norm(normal) ≈ 1.0
    slope_deg = slope_degrees_at(slope, normal)
    slope_deg === nothing && return random_surface_normal(slope.fallback, normal)
    theta = deg2rad(slope_deg)
    phi = rand(Uniform(0.0, 2π))
    direction = spherical_to_cartesian(theta, phi, 1.0)
    zhat = SA[0.0, 0.0, 1.0]
    θ = acos(clamp(zhat ⋅ normal, -1.0, 1.0))
    axis = normal × zhat
    norm(axis) < 1e-6 && return direction
    rotated = AngleAxis(-θ, (axis / norm(axis))...) * direction
    @toggled_assert norm(rotated) ≈ 1.0
    return rotated
end

"""
    ShardedSouthPolarDEMNormalMap([directory]; fallback=GaussianSlope(7.6),
                                  cache_shards=16)

Location-aware 5 m/pixel south-polar surface model built from the LDEM
elevation raster.  Unlike the LDSM slope-angle maps, each map pixel includes
both a slope magnitude and an azimuth for the outward-normal tilt, derived
from central differences of the DEM.  Thus valid pixels use a deterministic
map-derived surface normal rather than drawing an isotropic slope azimuth.

The default map is restricted to 87--90 degrees south and stores 1 degree
slope magnitude plus 2 degree normal-tilt azimuth in compact UInt16 NPZ
shards.  Only `cache_shards` decoded shards are retained in memory. Missing
or out-of-range pixels use `fallback`.
"""
mutable struct ShardedSouthPolarDEMNormalMap{F <: SlopeModel} <: SlopeModel
    directory::String
    shape::Tuple{Int, Int}
    shard_size_px::Int
    pixel_size_m::Float64
    origin_x_m::Float64
    origin_y_m::Float64
    radius_m::Float64
    min_latitude_deg::Float64
    slope_step_deg::Float64
    aspect_step_deg::Float64
    aspect_bins::Int
    nodata_code::UInt16
    fallback::F
    cache_limit::Int
    cache::Dict{Tuple{Int, Int}, Matrix{UInt16}}
    cache_recency::Vector{Tuple{Int, Int}}
    cache_lock::ReentrantLock
end

function ShardedSouthPolarDEMNormalMap(
    directory::AbstractString=DEFAULT_5M_SOUTH87_DEM_NORMAL_SHARDS;
    fallback::SlopeModel=GaussianSlope(7.6),
    cache_shards::Integer=16,
)
    cache_shards > 0 || throw(ArgumentError("cache_shards must be positive."))
    manifest_path = joinpath(directory, "manifest.npz")
    isfile(manifest_path) || throw(ArgumentError(
        "South-polar DEM-normal shards are missing at $directory. Build local terrain products using " *
        "scripts/terrain/README.md."
    ))
    table = npzread(manifest_path)
    required = (
        "shape", "shard_size_px", "pixel_size_m", "origin_x_m", "origin_y_m",
        "radius_m", "min_latitude_deg", "slope_step_deg", "aspect_step_deg",
        "aspect_bins", "nodata_code",
    )
    missing = filter(key -> !haskey(table, key), required)
    isempty(missing) || throw(ArgumentError(
        "DEM-normal shard manifest at $manifest_path is missing: $(join(missing, ", "))"
    ))
    shape_values = Int.(table["shape"])
    length(shape_values) == 2 || throw(ArgumentError("manifest shape must have two entries."))
    shape = (shape_values[1], shape_values[2])
    all(>(0), shape) || throw(ArgumentError("manifest shape entries must be positive."))
    aspect_bins = Int(only(table["aspect_bins"]))
    aspect_bins > 0 || throw(ArgumentError("manifest aspect_bins must be positive."))

    return ShardedSouthPolarDEMNormalMap(
        String(directory), shape, Int(only(table["shard_size_px"])),
        Float64(only(table["pixel_size_m"])), Float64(only(table["origin_x_m"])),
        Float64(only(table["origin_y_m"])), Float64(only(table["radius_m"])),
        Float64(only(table["min_latitude_deg"])), Float64(only(table["slope_step_deg"])),
        Float64(only(table["aspect_step_deg"])), aspect_bins,
        UInt16(only(table["nodata_code"])), fallback, Int(cache_shards),
        Dict{Tuple{Int, Int}, Matrix{UInt16}}(), Tuple{Int, Int}[], ReentrantLock(),
    )
end

function _dem_normal_shard_filename(
    map::ShardedSouthPolarDEMNormalMap, shard_row::Int, shard_col::Int
)
    row_name = lpad(string(shard_row), 3, '0')
    col_name = lpad(string(shard_col), 3, '0')
    return joinpath(map.directory, "shard_r$(row_name)_c$(col_name).npz")
end

function _cached_dem_normal_shard(
    map::ShardedSouthPolarDEMNormalMap, shard_row::Int, shard_col::Int
)
    key = (shard_row, shard_col)
    lock(map.cache_lock)
    try
        if haskey(map.cache, key)
            index = findfirst(==(key), map.cache_recency)
            deleteat!(map.cache_recency, index)
            push!(map.cache_recency, key)
            return map.cache[key]
        end

        filename = _dem_normal_shard_filename(map, shard_row, shard_col)
        isfile(filename) || throw(ArgumentError("Missing DEM-normal shard $filename"))
        codes = npzread(filename)["normal_code"]
        codes isa Matrix{UInt16} || throw(ArgumentError(
            "DEM-normal shard $filename must contain a UInt16 normal_code matrix."
        ))
        if length(map.cache_recency) == map.cache_limit
            evicted = popfirst!(map.cache_recency)
            delete!(map.cache, evicted)
        end
        map.cache[key] = codes
        push!(map.cache_recency, key)
        return codes
    finally
        unlock(map.cache_lock)
    end
end

function _dem_normal_code_at(map::ShardedSouthPolarDEMNormalMap, normal)
    @toggled_assert norm(normal) ≈ 1.0
    lat, lon = cartesian_to_latlon(normal)
    lat <= map.min_latitude_deg || return nothing
    lat_rad = deg2rad(lat)
    lon_rad = deg2rad(mod(lon, 360.0))
    ρ = 2.0 * map.radius_m * tan((π / 2.0 + lat_rad) / 2.0)
    isfinite(ρ) || return nothing
    x = ρ * sin(lon_rad)
    y = -ρ * cos(lon_rad)
    isfinite(x) && isfinite(y) || return nothing

    col = floor(Int, (x - map.origin_x_m) / map.pixel_size_m) + 1
    row = floor(Int, (map.origin_y_m - y) / map.pixel_size_m) + 1
    1 <= row <= map.shape[1] && 1 <= col <= map.shape[2] || return nothing
    shard_row = (row - 1) ÷ map.shard_size_px
    shard_col = (col - 1) ÷ map.shard_size_px
    codes = _cached_dem_normal_shard(map, shard_row, shard_col)
    local_row = (row - 1) % map.shard_size_px + 1
    local_col = (col - 1) % map.shard_size_px + 1
    checkbounds(Bool, codes, local_row, local_col) || return nothing
    code = codes[local_row, local_col]
    return code == map.nodata_code ? nothing : Int(code)
end

function slope_degrees_at(map::ShardedSouthPolarDEMNormalMap, normal)
    code = _dem_normal_code_at(map, normal)
    code === nothing && return nothing
    return (code ÷ map.aspect_bins) * map.slope_step_deg
end

function random_surface_normal(slope::ShardedSouthPolarDEMNormalMap, normal)
    @toggled_assert norm(normal) ≈ 1.0
    code = _dem_normal_code_at(slope, normal)
    code === nothing && return random_surface_normal(slope.fallback, normal)
    slope_code, aspect_code = divrem(code, slope.aspect_bins)
    theta = deg2rad(slope_code * slope.slope_step_deg)
    alpha = deg2rad(aspect_code * slope.aspect_step_deg)

    # The stored azimuth is measured in the south-polar stereographic grid.
    # At a location with longitude λ, projected +x and +y map to the local
    # tangent vectors below.  This turns the DEM's normal direction into the
    # global Cartesian basis used throughout CoRaLS.
    lat, lon = cartesian_to_latlon(normal)
    lat_rad = deg2rad(lat)
    lon_rad = deg2rad(mod(lon, 360.0))
    east = SA[-sin(lon_rad), cos(lon_rad), 0.0]
    north = SA[-sin(lat_rad) * cos(lon_rad), -sin(lat_rad) * sin(lon_rad), cos(lat_rad)]
    map_x = cos(lon_rad) * east + sin(lon_rad) * north
    map_y = sin(lon_rad) * east - cos(lon_rad) * north
    tilt = cos(alpha) * map_x + sin(alpha) * map_y
    surface_normal = cos(theta) * normal + sin(theta) * tilt
    @toggled_assert norm(surface_normal) ≈ 1.0
    return surface_normal
end

"""
    slope_degrees_at(map, normal)

Return the stored LDSM slope magnitude, in degrees, at a spherical surface
normal. Return `nothing` outside the configured south-polar range, outside the
raster, or at a missing raster pixel.
"""
function slope_degrees_at(map::SouthPolarSlopeMap, normal)
    @toggled_assert norm(normal) ≈ 1.0

    lat, lon = cartesian_to_latlon(normal)
    lat <= map.min_latitude_deg || return nothing

    lat_rad = deg2rad(lat)
    lon_rad = deg2rad(mod(lon, 360.0))
    # South-pole stereographic coordinates, matching the PSR-LUT convention.
    ρ = 2.0 * map.radius_m * tan((π / 2.0 + lat_rad) / 2.0)
    isfinite(ρ) || return nothing
    x = ρ * sin(lon_rad)
    y = -ρ * cos(lon_rad)
    isfinite(x) && isfinite(y) || return nothing

    col = floor(Int, (x - map.origin_x_m) / map.pixel_size_m) + 1
    row = floor(Int, (map.origin_y_m - y) / map.pixel_size_m) + 1
    checkbounds(Bool, map.slope_deg, row, col) || return nothing

    slope_deg = map.slope_deg[row, col]
    return isfinite(slope_deg) ? Float64(slope_deg) : nothing
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
    random_surface_normal(slope::SouthPolarSlopeMap, normal)

Use the stored local slope magnitude when `normal` lies within the configured
map region. Since this LDSM product has no aspect layer, draw the tilt azimuth
isotropically. Missing and out-of-range locations defer to the model's
configured fallback slope model.
"""
function random_surface_normal(slope::SouthPolarSlopeMap, normal)
    @toggled_assert norm(normal) ≈ 1.0

    slope_deg = slope_degrees_at(slope, normal)
    slope_deg === nothing && return random_surface_normal(slope.fallback, normal)

    theta = deg2rad(slope_deg)
    phi = rand(Uniform(0.0, 2π))
    direction = spherical_to_cartesian(theta, phi, 1.0)
    zhat = SA[0.0, 0.0, 1.0]
    θ = acos(clamp(zhat ⋅ normal, -1.0, 1.0))
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
