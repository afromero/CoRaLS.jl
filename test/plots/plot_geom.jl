using CoRaLS
using PyPlot
using Unitful: eV, km, sr, yr, ustrip
using Test
using LinearAlgebra
using Random

function plot_orbit_sampling()
    Random.seed!(54265324)

    # Plot
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["axes.axisbelow"] = true
    rcParams["axes.grid"] = true
    rcParams["grid.linestyle"] = "dashed"
    rcParams["xtick.direction"] = "in"
    rcParams["ytick.direction"] = "in"
    
    fig = plt.figure(figsize=(8, 8))

    ntrials = 1000
    ax = fig.add_subplot(111, projection="3d")     
    orbit = create_spacecraft("file:lro_orbit_1yr_2010.csv")
    SCs = [get_position(orbit) for _ in 1:ntrials] ./ km
    
    for sc in SCs
        sc = normalize(sc)
        ax.scatter(sc[1], sc[2], sc[3], color="k", s=20, alpha=0.5)
    end
    ax.set_xlim([-1,1])
    ax.set_ylim([-1,1])
    ax.set_zlim([-1,1])
    ax.set_aspect("equal")
    plt.tight_layout()
    fig.savefig("$(@__DIR__)/test_plots/orbit_sampled.png")
end


function plot_orbit_circular()
    Random.seed!(54265325)

    # Plot
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["axes.axisbelow"] = true
    rcParams["axes.grid"] = true
    rcParams["grid.linestyle"] = "dashed"
    rcParams["xtick.direction"] = "in"
    rcParams["ytick.direction"] = "in"
    
    fig = plt.figure(figsize=(8, 8))

    ntrials = 1000
    orbit = create_spacecraft("orbit:50km")
    ax = fig.add_subplot(111, projection="3d")            
    SCs = [get_position(orbit) for _ in 1:ntrials] ./ km
    _ = [ax.scatter(sc..., color="k", s=20, alpha=0.5) for sc in SCs]
    ax.set_aspect("equal")
    plt.tight_layout()
    fig.savefig("$(@__DIR__)/test_plots/orbit_circular.png")
end

function plot_orbit_elliptical()
    Random.seed!(54265326)

    # Plot
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["axes.axisbelow"] = true
    rcParams["axes.grid"] = true
    rcParams["grid.linestyle"] = "dashed"
    rcParams["xtick.direction"] = "in"
    rcParams["ytick.direction"] = "in"
    
    fig = plt.figure(figsize=(8, 14))
    fig.suptitle("Periapse:3000km, Apoapse:8000km, inc:30°", fontsize=16)

    ntrials = 1000
    orbit = create_spacecraft("elliptical:3000km,8000km,30.0")
    SCs = [get_position(orbit) for _ in 1:ntrials] ./ km

    # Side profile
    ax = fig.add_subplot(311, projection="3d") 
    ax.view_init(elev=0, azim=25)           
    _ = [ax.scatter(sc..., color="k", s=20, alpha=0.5) for sc in SCs]
    ax.set_aspect("equal")

    # Oblique view
    ax = fig.add_subplot(312, projection="3d") 
    ax.view_init(elev=30, azim=25)           
    _ = [ax.scatter(sc..., color="k", s=20, alpha=0.5) for sc in SCs]
    ax.set_aspect("equal")

    # Oblique view
    ax = fig.add_subplot(313, projection="3d") 
    ax.view_init(elev=90, azim=0)           
    _ = [ax.scatter(sc..., color="k", s=20, alpha=0.5) for sc in SCs]
    ax.set_aspect("equal")

    plt.tight_layout()
    fig.savefig("$(@__DIR__)/test_plots/orbit_elliptical.png")
end


function plot_spherical_sampling()
    Random.seed!(54265327)

    # Plot
    rcParams = PyPlot.PyDict(PyPlot.matplotlib."rcParams")
    rcParams["axes.axisbelow"] = true
    rcParams["axes.grid"] = true
    rcParams["grid.linestyle"] = "dashed"
    rcParams["xtick.direction"] = "in"
    rcParams["ytick.direction"] = "in"
    
    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111, projection="3d")
    for i = 1:1000
        surface = random_point_on_sphere()
    
        ax.scatter(surface[1], surface[2], surface[3], color="k", s=20, alpha=0.5)
    end
    ax.set_xlim([-1,1])
    ax.set_ylim([-1,1])
    ax.set_zlim([-1,1])
    ax.set_aspect("equal")
    plt.tight_layout()
    fig.savefig("$(@__DIR__)/test_plots/spherical_sampling.png")

end

@testset "Test geometry plots" begin

    plot_orbit_sampling()
    plot_orbit_circular()
    plot_orbit_elliptical()
    plot_spherical_sampling()
    @test 1 == 1  # Finished with no errors
end