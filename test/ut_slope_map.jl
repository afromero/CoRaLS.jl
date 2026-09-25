using CoRaLS: SouthPolarSlopeMap, SouthPolarSlopeMap40m, ShardedSouthPolarSlopeMap,
    ShardedSouthPolarDEMNormalMap, GaussianSlope, slope_degrees_at,
    random_surface_normal, latlon_to_cartesian, dot, norm
using NPZ: npzwrite

const _TEST_MAP_PIXEL_SIZE_M = 2.0e6
const _TEST_MAP_ORIGIN_M = 1.0e6

function _write_test_slope_map(path)
    npzwrite(path, Dict(
        "slope_deg" => reshape(Float32[12.0], 1, 1),
        "pixel_size_m" => Float64[_TEST_MAP_PIXEL_SIZE_M],
        "origin_x_m" => Float64[-_TEST_MAP_ORIGIN_M],
        "origin_y_m" => Float64[_TEST_MAP_ORIGIN_M],
        "radius_m" => Float64[1_737_400.0],
    ))
end

function _write_test_slope_shards(directory)
    mkpath(directory)
    npzwrite(joinpath(directory, "manifest.npz"), Dict(
        "shape" => Int64[1, 1], "shard_size_px" => Int64[1],
        "pixel_size_m" => Float64[_TEST_MAP_PIXEL_SIZE_M],
        "origin_x_m" => Float64[-_TEST_MAP_ORIGIN_M],
        "origin_y_m" => Float64[_TEST_MAP_ORIGIN_M], "radius_m" => Float64[1_737_400.0],
        "min_latitude_deg" => Float64[-87.0], "slope_step_deg" => Float64[1.0],
        "nodata_code" => UInt8[255],
    ))
    npzwrite(joinpath(directory, "shard_r000_c000.npz"), Dict(
        "slope_code" => reshape(UInt8[12], 1, 1),
    ))
end

function _write_test_dem_normal_shards(directory)
    mkpath(directory)
    npzwrite(joinpath(directory, "manifest.npz"), Dict(
        "shape" => Int64[1, 1], "shard_size_px" => Int64[1],
        "pixel_size_m" => Float64[_TEST_MAP_PIXEL_SIZE_M],
        "origin_x_m" => Float64[-_TEST_MAP_ORIGIN_M],
        "origin_y_m" => Float64[_TEST_MAP_ORIGIN_M], "radius_m" => Float64[1_737_400.0],
        "min_latitude_deg" => Float64[-87.0], "slope_step_deg" => Float64[1.0],
        "aspect_step_deg" => Float64[2.0], "aspect_bins" => Int64[180],
        "nodata_code" => UInt16[65535],
    ))
    # 10 degree tilt and 0 degree projected-grid azimuth.
    npzwrite(joinpath(directory, "shard_r000_c000.npz"), Dict(
        "normal_code" => reshape(UInt16[1800], 1, 1),
    ))
end

@testset "south-polar slope map" begin
    mktempdir() do directory
        path = joinpath(directory, "slope.npz")
        _write_test_slope_map(path)
        map = SouthPolarSlopeMap(path, fallback=GaussianSlope(0.0))
        inside = latlon_to_cartesian(-89.0, 0.0)
        outside = latlon_to_cartesian(-86.0, 0.0)
        @test slope_degrees_at(map, inside) == 12.0
        @test slope_degrees_at(map, outside) === nothing
        sampled_normal = random_surface_normal(map, inside)
        sampled_tilt_deg = rad2deg(acos(clamp(dot(sampled_normal, inside), -1.0, 1.0)))
        @test sampled_tilt_deg ≈ 12.0 atol=1e-10
        @test random_surface_normal(map, outside) ≈ outside
    end
end

@testset "40 m south-polar slope map constructor" begin
    mktempdir() do directory
        path = joinpath(directory, "slope.npz")
        _write_test_slope_map(path)
        map = SouthPolarSlopeMap40m(path, fallback=GaussianSlope(0.0))
        @test slope_degrees_at(map, latlon_to_cartesian(-89.0, 0.0)) == 12.0
    end
end

@testset "sharded south-polar slope map" begin
    mktempdir() do directory
        _write_test_slope_shards(directory)
        map = ShardedSouthPolarSlopeMap(directory, fallback=GaussianSlope(0.0), cache_shards=1)
        @test slope_degrees_at(map, latlon_to_cartesian(-89.0, 0.0)) == 12.0
        @test random_surface_normal(map, latlon_to_cartesian(-86.0, 0.0)) ≈ latlon_to_cartesian(-86.0, 0.0)
    end
end

@testset "south-polar DEM normal map" begin
    mktempdir() do directory
        _write_test_dem_normal_shards(directory)
        map = ShardedSouthPolarDEMNormalMap(directory, fallback=GaussianSlope(0.0), cache_shards=1)
        inside = latlon_to_cartesian(-89.0, 0.0)
        outside = latlon_to_cartesian(-86.0, 0.0)
        @test slope_degrees_at(map, inside) == 10.0
        sampled_normal = random_surface_normal(map, inside)
        @test norm(sampled_normal) ≈ 1.0
        sampled_tilt_deg = rad2deg(acos(clamp(dot(sampled_normal, inside), -1.0, 1.0)))
        @test sampled_tilt_deg ≈ 10.0 atol=1e-10
        @test random_surface_normal(map, outside) ≈ outside
    end
end
