using CoRaLS: sample_power_law 
using CoRaLS: auger_spectrum_2020 
using CoRaLS: auger_spectrum_2021 
using Unitful: g, cm, m, km, sr, μV, V, eV, GeV, EeV, Hz, MHz, W, K, NoUnits, yr, ustrip

@testset verbose = true "spectrum.jl" begin
    @testset "Test auger_spectrum_2021 flux" begin
        # Raw counts from Table 10 of Auger 2021 paper
        J_expected = [6.431e-14, 3.191e-14, 1.577e-14, 7.643e-15, 3.650e-15, 1.739e-15, 8.32e-16, 3.90e-16, 1.85e-16, 8.87e-17, 4.14e-17, 1.9e-17, 8.47e-18, 4.17e-18, 1.929e-18, 9.041e-19, 4.294e-19, 2.167e-19, 1.226e-19, 6.82e-20, 3.79e-20, 2.07e-20, 1.04e-20, 0.53e-20, 2.49e-21, 1.25e-21, 5.99e-22, 1.95e-22, 8.1e-23, 1.8e-23, 5.5e-24, 2.9e-24]
        log_bin_centers = range(17.05, 20.15+0.05, length=length(J_expected))
        bins = 10 .^ log_bin_centers

        # Compute flux from the Auger spectrum
        J_auger = auger_spectrum_2021.(bins .* eV)  # [1 / (km^2 sr yr eV)]
        J_auger = round.(ustrip.(J_auger), sigdigits=3)  # cancel units for comparison

        # Agrees to within 2.5% - not sure what causes the offset, maybe the binning? 
        @test J_auger ≈ J_expected rtol = 0.025  # 2.5% agreement
    end

    @testset "Test auger_spectrum_2020 flux" begin
        # Raw counts from Fig. 7 of Auger 2020 paper
        raw_count = [83143, 47500, 28657, 17843, 12435, 8715, 6050, 4111, 2620, 1691, 991, 624, 372, 156, 83, 24, 9, 6, 0, 0]
        
        # Based on code from "The energy spectrum" Jupyter tutorial on the Auger open 
        # data website https://www.auger.org/opendata/ (Accessed 07/2024). Uses data
        # from Auger Open Data summary.zip (https://doi.org/10.5281/zenodo.10488964)k
        log_E_min = 0.4
        E_bins = 20
        E_bin_size = 0.1
        log_E_max = log_E_min + E_bins * E_bin_size
        log_bins = range(log_E_min, log_E_max, length = E_bins + 1)
        log_bin_centers = log_bins[1:end-1] .+ 0.05
        bins = 10 .^ log_bins .* 1e18
        bin_energy = 10 .^ log_bin_centers .* 1e18
        bin_width = diff(bins)
        # Compute flux from the Auger spectrum
        # Note: the Auger data is Jraw, but the spectral fit computes J (which we want)
        # Technically we're comparing two different things but works out to same
        # shape and < 5% differnce so parameterization looks good
        J_bins = auger_spectrum_2020.(bins .* eV)  # [1 / (km^2 sr yr eV)]

        # Apply effective sampling area-time and trapezoidal integration over energy to get counts
        count_rate = bin_width * eV .* (J_bins[1:end-1] .+ J_bins[2:end]) ./ 2 # [1 / (km^2 sr yr)]
        count = count_rate * 60400 * km^2 * sr * yr  # Auger 2020 paper

        # 5% tolerance, some offset probably due to Jraw -> J conversion
        @test count ≈ raw_count rtol = 0.05  
    end

end
;
