using CoRaLS: DataDrivenSlope, NoSlope, available_slope_models, random_surface_normal,
    slope_model_from_name, spherical_to_cartesian
using Random

@testset verbose = true "surface.jl" begin
    @testset "Data-driven 5 m slope model" begin
        Random.seed!(20260917)

        slope = DataDrivenSlope()
        normal = spherical_to_cartesian(0.0, 0.0, 1.0)
        samples = [random_surface_normal(slope, normal) for _ in 1:20_000]
        angles_deg = rad2deg.([acos(clamp(sample[3], -1.0, 1.0)) for sample in samples])

        # These 5 m analysis values allow for Monte Carlo variation and the
        # uniform-in-bin sampling used by DataDrivenSlope.
        @test all(0.0 .<= angles_deg .< 70.5)
        @test sum(angles_deg) / length(angles_deg) ≈ 10.186 atol = 0.25
        @test sqrt(sum(abs2, angles_deg) / length(angles_deg)) ≈ 12.012 atol = 0.25
    end

    @testset "Slope model names" begin
        @test "data_5m" in available_slope_models()
        @test slope_model_from_name("data_5m") isa DataDrivenSlope
        @test slope_model_from_name("none") isa NoSlope
        @test_throws ArgumentError slope_model_from_name("not-a-model")
    end
end
