using CoRaLS
using CoRaLS: Unitful, @u_str, s, yr, H, eV, V, K, MHz, cm, m, km, sr, ustrip, c_0
using Test
# using Random
# Set a random seed for reproducible tests
# SEED = 10  # trunc(Int, time())
# println("# Random seed ($SEED), on VERSION == $VERSION")
# Random.seed!(SEED)

@testset "CoRaLS.jl Tests" begin
    include("ut_acceptance.jl");
    include("ut_geometry.jl");
    include("ut_fresnel.jl");
    include("ut_ice.jl");    
    include("ut_spectrum.jl");
    include("ut_plots.jl");
end

# Cleanup plot outputs so only visually changed PNG files remain in git diff.
include("plots/visual_diff_test_plots.jl")