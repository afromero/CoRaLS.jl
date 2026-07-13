using CoRaLS
using PyPlot
using Unitful
using Unitful: EeV, eV, km, sr, yr, ustrip
using Random
using Test

function plot_auger_2021_comparison()
    Random.seed!(54120001)
    # Raw counts from Table 10 of Auger 2021 paper
    J_expected = [6.431e-14, 3.191e-14, 1.577e-14, 7.643e-15, 3.650e-15, 1.739e-15, 8.32e-16, 3.90e-16, 1.85e-16, 8.87e-17, 4.14e-17, 1.9e-17, 8.47e-18, 4.17e-18, 1.929e-18, 9.041e-19, 4.294e-19, 2.167e-19, 1.226e-19, 6.82e-20, 3.79e-20, 2.07e-20, 1.04e-20, 0.53e-20, 2.49e-21, 1.25e-21, 5.99e-22, 1.95e-22, 8.1e-23, 1.8e-23, 5.5e-24, 2.9e-24]
    bin_centers = range(17.05, 20.15+0.05, length=length(J_expected))
    bins = 10 .^ bin_centers

    # Auger spectrum parameterization
    J_auger = auger_spectrum_2021.(bins .* eV)  # [1 / (km^2 sr yr eV)]
    J_auger = round.(ustrip.(J_auger), sigdigits=3)  # cancel units for plot

    # Show rejection sampled spectrum
    ntrials = 1e4
    bin_edges = 10 .^ range(17, 20.25, length(J_expected)+1)
    rsamp = [sample_auger(u"EeV"(10^17.0 * eV), u"EeV"(10^20.2 * eV), auger_func=auger_spectrum_2021) / eV |> NoUnits for _ = 1:ntrials]


    # Plot
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["axes.axisbelow"] = true
    rcParams["axes.grid"] = true
    rcParams["grid.linestyle"] = "dashed"
    rcParams["xtick.direction"] = "in"
    rcParams["ytick.direction"] = "in"

    fig, axs = plt.subplots(1, 1, figsize=(8, 6))
    axs.scatter(bins, J_expected, color="r",  label="J (Table 10, Auger 2020)")
    axs.loglog(bins, J_auger, color="k", linestyle="--", label="Auger 2021 parameterization")
    axs.hist(rsamp, bins=bin_edges, density=true, histtype="step", label="Rejection sampled", zorder=0)

    [axs.annotate(eJ, (E, eJ), rotation=30, color="r") for (E, eJ) in zip(bins, J_expected)]
    [axs.annotate(J, (E, J), ha="right", va="top", rotation=30) for (E, J) in zip(bins, J_auger)]
    axs.set_xlabel("Energy [eV]")
    axs.set_ylabel("J [1 / (km^2 yr sr eV)]")
    axs.set_xlim(5e16, 5e20)
    axs.set_ylim(1e-26, 1e-12)
    axs.legend()
    fig.savefig("$(@__DIR__)/test_plots/auger21_comparison.png")

    @test J_auger ≈ J_expected rtol = 0.025  # 13% offset

end

function plot_auger_2020_comparison()
    Random.seed!(54120002)
    # Raw counts from Fig. 7 of Auger 2020 paper
    raw_count = [83143, 47500, 28657, 17843, 12435, 8715, 6050, 4111, 2620, 1691, 991, 624, 372, 156, 83, 24, 9, 6, 0, 0]

    # Init arrays
    log_E_min = 0.4
    E_bins = 20
    E_bin_size = 0.1
    log_E_max = log_E_min + E_bins * E_bin_size
    log_bins = range(log_E_min, log_E_max, length = E_bins + 1)
    log_bin_centers = log_bins[1:end-1] .+ 0.05
    bins = 10 .^ log_bins .* 1e18
    bin_energy = 10 .^ log_bin_centers .* 1e18
    bin_width = diff(bins)

    # Auger spectrum parameterization
    J_auger = auger_spectrum_2020.(bin_energy .* eV)  # [1 / (km^2 sr yr eV)]
    ntrials = 1e4
    rsamp = [sample_auger(u"EeV".(bin_energy[1].* eV), u"EeV".(bin_energy[end].* eV), auger_func=auger_spectrum_2020) / eV |> NoUnits for _ = 1:ntrials]

    # Trapezoidal integration over energy to get counts
    # Note: the Auger data is Jraw, but the spectral fit computes J
    J_bins = auger_spectrum_2020.(bins .* eV)  # [1 / (km^2 sr yr eV)]
    count_rate = bin_width * eV .* (J_bins[1:end-1] .+ J_bins[2:end]) ./ 2 # [1 / (km^2 sr yr)]
    count = floor.(Int, count_rate * 60400 * km^2 * sr * yr)  # effective experiment sampling area and time

    # Plot
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["axes.axisbelow"] = true
    rcParams["axes.grid"] = true
    rcParams["grid.linestyle"] = "dashed"
    rcParams["xtick.direction"] = "in"
    rcParams["ytick.direction"] = "in"

    fig, axs = plt.subplots(2, 1, figsize=(6, 8))
    axs[1].loglog(bin_energy, ustrip.(J_auger), label="Auger 2020 parameterization")
    axs[1].hist(rsamp, bins=bins, density=true, histtype="step", label="Auger 2020 sampled")
    
    axs[1].set_ylabel("J (E) [1 / (km^2 yr sr eV)]")
    # axs[1].set_ylim(1e-25, 1e-17)
    axs[1].legend()
    axs[2].scatter(bin_energy, raw_count, color="r",  label="Raw counts, Fig. 7, Auger 2020")
    axs[2].loglog(bin_energy, count, color="k", linestyle="--", label="Counts integrated from J")
    [axs[2].annotate(rc, (E, rc), rotation=30, color="r") for (E, rc) in zip(bin_energy, raw_count)]
    [axs[2].annotate(c, (E, c), ha="right", va="top", rotation=30) for (E, c) in zip(bin_energy, count)]
    axs[2].set_xlabel("Energy [eV]")
    axs[2].set_ylabel("Events per bin")
    axs[2].set_xlim(1.5e18, 2e20)
    axs[2].set_ylim(5e-1, 1e6)
    axs[2].legend()
    fig.savefig("$(@__DIR__)/test_plots/auger20_comparison.png")
end

function plot_sample_power_law()
    Random.seed!(43142432)

    nbins = 100
    samples = 10^6
    bin_edges = range(1.0, 1001.0, length = nbins+1)
    bin_centers = (bin_edges[1:end-1] .+ bin_edges[2:end]) ./ 2

    ## Define the functions
    y_inv    = ((bin_centers .^ (-1)) ./ log(1001/1)) .* samples .* (1001 - 1) ./ nbins
    y_inv_sq = ((bin_centers .^ (-2)) .* 1001/1000)   .* samples .* (1001 - 1) ./ nbins

    ## Sample from the functions and then histogam them
    inv_sample    = sample_power_law(-1.0, samples, min_value = 1.0EeV, max_value = 1001.0EeV)
    inv_sq_sample = sample_power_law(-2.0, samples, min_value = 1.0EeV, max_value = 1001.0EeV)
    
    ## Make some plots
    fig1, ax1 = plt.subplots(1, 1, figsize=(8, 6))
    ax1.hist(ustrip.(inv_sample), bins=collect(bin_edges), label="samples")
    ax1.set_yscale("log")
    ax1.set_title("gamma = -1")
    ax1.set_xlabel("x")
    ax1.set_ylabel("counts")
    ax1.legend()
    fig1.savefig("$(@__DIR__)/test_plots/sample_E_-1_power_law_flux.png")
    close(fig1)

    fig2, ax2 = plt.subplots(1, 1, figsize=(8, 6))
    ax2.hist(ustrip.(inv_sq_sample), bins=collect(bin_edges), label="samples")
    ax2.set_yscale("log")
    ax2.set_title("gamma = -2")
    ax2.set_xlabel("x")
    ax2.set_ylabel("counts")
    ax2.legend()
    fig2.savefig("$(@__DIR__)/test_plots/sample_E_-2_power_law_flux.png")
    close(fig2)

    ## For now, this test just makes the plots and we'll implement a more sophisticated on later
    @test 1 == 1
end

@testset "Test spectrum.jl plots" begin
    plot_auger_2020_comparison()
    plot_auger_2021_comparison()
    plot_sample_power_law()

    @test 1 == 1  # Finished with no errors
end
