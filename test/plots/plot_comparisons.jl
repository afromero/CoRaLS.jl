using CoRaLS
using Unitful: km, sr, EeV, MHz, m, μV, cm, ustrip
using Random
using PyPlot

function plot_ice_depth(ntrials=1000, nbins=7, ice_depths=[1, 5, 10, 100]m)
    Random.seed!(54110001)
    figname = "$(@__DIR__)/test_plots/ice_depth_experiment.png"
    region = create_region("polar:south,-80,0.0557")
    altitude = 10km
    title="Spole PSR detections\naltitude=$(altitude) (4 nadir 0.3-1GHz)"

    kws = Dict(:min_energy=>0.5EeV,:max_energy=>500.0EeV,
            :ν_min=>300MHz,:ν_max=>1000MHz,:dν=>30MHz,
            :min_count=>50,:max_tries=>100,:simple_area=>true)

    A_arr = []
    for ice_depth in ice_depths
        sc = CircularOrbit(altitude)  
        trigger = LPDA(Nant=4, Ntrig=4, θ0=-90.0, altitude=altitude, skyfrac=0)
        A = acceptance(ntrials, nbins, region=region, spacecraft=sc; trigger=trigger, ice_depth=ice_depth, kws...)
        push!(A_arr, A)
    end

    fig, ax = plot_rate_experiment(A_arr, ice_depths; xlabel="Altitude [km]")
    ax.set_title(title)
    fig.savefig(figname)
end

function plot_altitudes(ntrials=1000, nbins=7, altitudes=[1, 10, 50, 100, 500, 1000]km)
    Random.seed!(54110002)
    figname = "$(@__DIR__)/test_plots/altitude_experiment.png"
    region = create_region("polar:south,-80,0.0557")
    ice_depth = 10.0m
    title="Spole PSR detections\nice_depth=$(ice_depth) (4 nadir 0.3-1GHz)"

    kws = Dict(:ice_depth=>ice_depth,:min_energy=>0.5EeV,:max_energy=>500.0EeV,
            :ν_min=>300MHz,:ν_max=>1000MHz,:dν=>30MHz,
            :min_count=>50,:max_tries=>100,:simple_area=>true)

    A_arr = []
    for altitude in altitudes
        sc = CircularOrbit(altitude)  
        trigger = LPDA(Nant=4, Ntrig=4, θ0=-90.0, altitude=altitude, skyfrac=0)
        A = acceptance(ntrials, nbins, region=region, spacecraft=sc; trigger=trigger, kws...)
        push!(A_arr, A)
    end

    fig, ax = plot_rate_experiment(A_arr, altitudes; xlabel="Altitude [km]")
    ax.set_title(title)
    fig.savefig(figname)
end

# Compare acceptance for different detector geometry
# Next: look as fuction of angle with co-aligned antennas
#  SNR: trigger rate at kHz with nano sec window
function plot_detector_geoms(; verbose=false, from_dir="")
    Random.seed!(54110003)
    # Setup params
    figname = "$(@__DIR__)/test_plots/detector_geom_comparison.png"
    savefile = "$(@__DIR__)/../../tmp/compare_geoms.jld2"
    mkpath(dirname(savefile))
    nbins = 7
    ntrials = 50000
    max_tries = 10000
    min_count = 20
    altitude = 50km
    
    # Make triggers
    trigger_mag = magnitude_trigger(67μV/m)
    trigger_octo = LPDA(Nant=8, Ntrig=3, θ0=-48.0, altitude=altitude)
    trigger_nadir = LPDA(Nant=4, Ntrig=4, θ0=-90.0, altitude=altitude, skyfrac=0)
    trigger_solo = LPDA(Nant=1, Ntrig=1, θ0=-90.0, altitude=altitude, skyfrac=0)

    if !isempty(from_dir)
        # Read acceptance from multiple files and merge event lists into one Acceptance struct
        files = readdir(from_dir)
        geom_files = filter(f -> startswith(f, "compare_geoms_"), files)
        acc_solo = load_acceptance(from_dir * geom_files[1])
        for f in geom_files[2:end]
            acc_solo = merge_acceptance(acc_solo, load_acceptance("./tmp/" * f))
        end
    else
        # Run acceptance and collect events
        region = AllPSR
        sc = CircularOrbit(altitude)
        kws = Dict(:ice_depth=>6.0m, :min_energy=>3.0EeV, :max_energy=>200.0EeV,
                :ν_min=>300MHz, :ν_max=>1000MHz, :region=>region, :spacecraft=>sc,
                :min_count=>min_count, :max_tries=>max_tries, :savefile=>savefile)   
        acc_solo = acceptance(ntrials, nbins; trigger=trigger_solo, save_events=true, kws...)
    end

    # retrigger saved events on each other detector
    acc_mag = retrigger(acc_solo, trigger_mag)
    acc_octo = retrigger(acc_solo, trigger_octo)
    acc_nadir = retrigger(acc_solo, trigger_nadir)

    if verbose
        @info "Saved events:", length(acc_solo.reflected), length(acc_solo.direct)
        @info "magnitude counts (r, d):", acc_mag.rcount, acc_mag.dcount
        @info "octo (r,d):", acc_octo.rcount, acc_octo.dcount
        @info "nadir x4 counts (r,d):", acc_nadir.rcount, acc_nadir.dcount
        @info "nadir x1 counts (r,d):", acc_solo.rcount, acc_solo.dcount
    end

    # Make plots
    blues = get_cmap(:Blues_r)
    oranges = get_cmap(:Oranges_r)

    fig, ax = plt.subplots(figsize=(8, 8))
    f, ax = plot_acceptance(acc_mag; min_count=min_count, ax=ax, name="thresh=67μV/m", colorrefl=blues(0.1), colordirect=oranges(0.1))
    f, ax = plot_acceptance(acc_octo; min_count=min_count, ax=ax, name="N=8,trig=3,dip=48°", linestyle="--", colorrefl=blues(0.2), colordirect=oranges(0.2))
    f, ax = plot_acceptance(acc_nadir; min_count=min_count, ax=ax, name="N=4,trig=4,dip=90°", linestyle="-.", colorrefl=blues(0.4), colordirect=oranges(0.4))
    f, ax = plot_acceptance(acc_solo; min_count=min_count, ax=ax, name="N=1,trig=1,dip=90°", linestyle=":", colorrefl=blues(0.6), colordirect=oranges(0.6))
    ax.set_title("50km polar orbit, AllPSR, ice_depth=6m")
    fig.savefig(figname)
end

function plot_detector_coaligned(;verbose=false)
    Random.seed!(54110004)
    # Setup params
    figname = "$(@__DIR__)/test_plots/detector_coaligned.png"
    savefile = "$(@__DIR__)/../../tmp/compare_coaligned.jld2"
    mkpath(dirname(savefile))
    min_count = 20
    altitude = 50km
    region = AllPSR
    sc = CircularOrbit(altitude)
    kws = Dict(:ice_depth=>6.0m, :min_energy=>3.0EeV, :max_energy=>200.0EeV,
              :ν_min=>300MHz, :ν_max=>1000MHz, :region=>region, :spacecraft=>sc,
              :min_count=>min_count, :max_tries=>1000, :savefile=>savefile)
    trigger_n = LPDA(Nant=4, Ntrig=4, θ0=-45.0, ϕ0=0, altitude=altitude)
    trigger_ne = LPDA(Nant=4, Ntrig=4, θ0=-45.0, ϕ0=45, altitude=altitude)
    trigger_e = LPDA(Nant=4, Ntrig=4, θ0=-45.0, ϕ0=90, altitude=altitude)
    trigger_se = LPDA(Nant=4, Ntrig=4, θ0=-45.0, ϕ0=135, altitude=altitude)
    trigger_s = LPDA(Nant=4, Ntrig=4, θ0=-45.0, ϕ0=180, altitude=altitude)
    trigger_nadir = LPDA(Nant=4, Ntrig=4, θ0=-90.0, altitude=altitude)

    # Run acceptance and collect events
    ntrials = 1000
    nbins = 7
    a_nadir = acceptance(ntrials, nbins; trigger=trigger_nadir, save_events=true, kws...)
    @info "Saved events:", length(a_nadir.reflected), length(a_nadir.direct)

    # retrigger saved events on each other detector
    a_n = retrigger(a_nadir, trigger_n)
    a_ne = retrigger(a_nadir, trigger_ne)
    a_e = retrigger(a_nadir, trigger_e)
    a_se = retrigger(a_nadir, trigger_se)
    a_s = retrigger(a_nadir, trigger_s)

    if verbose
        @info "nadir x4 counts (r,d):", a_nadir.rcount, a_nadir.dcount
        @info "θ=45, N x4 counts (r,d):", a_n.rcount, a_n.dcount
        @info "θ=45, NE x4 counts (r,d):", a_ne.rcount, a_ne.dcount
        @info "θ=45, E x4 counts (r,d):", a_e.rcount, a_e.dcount
        @info "θ=45, SE x4 counts (r,d):", a_se.rcount, a_se.dcount
        @info "θ=45, S x4 counts (r,d):", a_s.rcount, a_s.dcount
    end

    # Make plots
    blues = get_cmap(:Blues_r)
    oranges = get_cmap(:Oranges_r)

    fig, ax = plt.subplots(figsize=(11, 9))
    f, ax = plot_acceptance(a_nadir; min_count=min_count, ax=ax, name="Nadir", colorrefl=blues(1), colordirect=oranges(1))
    f, ax = plot_acceptance(a_n; min_count=min_count, ax=ax, name="dip=45°, az=0°", linestyle="--", colorrefl=blues(0.85), colordirect=oranges(0.85))
    f, ax = plot_acceptance(a_ne; min_count=min_count, ax=ax, name="dip=45°, az=45°", linestyle="-.", colorrefl=blues(0.6), colordirect=oranges(0.6))
    f, ax = plot_acceptance(a_e; min_count=min_count, ax=ax, name="dip=45°, az=90°", linestyle=":", colorrefl=blues(0.45), colordirect=oranges(0.45))
    f, ax = plot_acceptance(a_se; min_count=min_count, ax=ax, name="dip=45°, az=135°", linestyle="-.", colorrefl=blues(0.3), colordirect=oranges(0.3))
    f, ax = plot_acceptance(a_s; min_count=min_count, ax=ax, name="dip=45°, az=180", linestyle="--", colorrefl=blues(0.1), colordirect=oranges(0.1))

    ax.set_title("4 coaligned antennas, 50km polar orbit, AllPSR, ice_depth=6m")
    fig.savefig(figname)
end

@testset "Test comparison plots" begin
    plot_ice_depth()
    plot_altitudes();
    # plot_detector_geoms(; from_dir="/mnt/d/data/corals/")
    plot_detector_coaligned()
    @test 1 == 1  # Finished with no errors
end