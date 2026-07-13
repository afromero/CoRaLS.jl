using Test
# Pulls testsets from test/plots/*.jl and runs them
plots_dir = joinpath(@__DIR__, "plots")
plot_test_files = sort(filter(f -> startswith(f, "plot_") && endswith(f, ".jl"), readdir(plots_dir)))

@testset "Plotting Unittests (saved to: test/plots/test_plots)" begin
    for plot_test_file in plot_test_files
        include(joinpath("plots", plot_test_file))
    end
end
