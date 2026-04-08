using PyPlot
using Unitful: m, km

const colordirect = "#F88D07"
const colorrefl = "#0772F8"

"""
`plots.jl` - Visualization Tools for CoRaLS Simulation Data

Key Functionalities:
- Plotting incident angles at the ice and the surface, along with the corresponding Fresnel coefficients.
- Visualizing polarization angles for direct and reflected events, using both standard and polar coordinates.
- Generating histograms for off-axis angles of both direct and reflected events.
- Displaying direct and reflected acceptance of the CoRaLS project through line plots.
- Creating differential event spectra plots, which are crucial for analyzing the frequency and distribution of cosmic ray events over different energy levels.
"""


"""
    plot_incident_angles(AΩ::Acceptance)

Plot the incident angles at the ice and the surface for both direct and reflected events, along with the corresponding Fresnel coefficients.

Arguments:
- `AΩ::Acceptance`: An object containing the acceptance calculation results.

Generates a two-part plot. The left plot shows histograms of incident angles at the surface for direct and reflected events. The right plot shows incident angles at the ice layer for reflected events. Behind these histograms, the function plots the transmission and reflection coefficients as a function of the incident angle.
"""
function plot_incident_angles(AΩ::Union{Acceptance, OldAcceptance})

    # create the figure
    fig, axes = plt.subplots(1, 2, figsize=(8, 3))

    # the bins for the tx histogram
    bins = 0:1:90

    # create the right-hand axes for the Fresnel coeffs
    raxes = [axes[1].twinx(), axes[2].twinx()]

    # histogram the observed incident angles
    axes[1].hist(AΩ.direct.θ_i, histtype="step",
        label=L"Direct $\theta_i$",
        color=colordirect, density=true, bins=bins, zorder=3)
    axes[1].hist(AΩ.reflected.θ_i, histtype="step",
        label=L"Refl. $\theta_i$", ls="dashed",
        color=colorrefl, density=true, bins=bins, zorder=3)

    # calculate the transmission coefficients at the surface
    Nsurf = regolith_index(StrangwayIndex(), 0.0m)
    Θ = 0:0.1:90.0 # incident angles at the surface
    T = [surface_transmission(GaussianRoughness(2.0),
        MixedFieldDivergence(),
        θ, Nsurf, 6.0m, 100.0km) for θ in deg2rad.(Θ)]

    # plot the transmission coefficients behind
    raxes[1].plot(Θ, [t[1] for t in T], alpha=0.6, zorder=0,
        color="grey", ls="dashed", label=L"Direct $t_{\parallel}$")
    raxes[1].plot(Θ, [t[2] for t in T], alpha=0.6, zorder=0,
        color="darkgrey", ls="dotted", label=L"Direct $t_{\perp}$")

    # histogram the observed incident angles at the ice
    axes[2].hist(AΩ.reflected.θ_ice, histtype="step",
        label=L"$\theta_{\mathrm{ice}}$",
        color=colorrefl, density=true, bins=bins, zorder=3)

    # calculate the reflection coefficients behind
    Nreg = regolith_index(StrangwayIndex(), 6.0m)
    Nice = 1.305
    Θ = 0.0:0.01:rad2deg(asin(Nice / Nreg)) # incident angles at the ice
    rpar = [fresnel_coeffs(θ, Nreg, Nice)[1] for θ in deg2rad.(Θ)]
    rperp = [fresnel_coeffs(θ, Nreg, Nice)[2] for θ in deg2rad.(Θ)]

    # plot the reflection coefficients behind
    raxes[2].plot(Θ, abs.(rpar), alpha=0.6, zorder=0,
        color="grey", ls="dashed", label=L"$|r_{\parallel}|$")
    raxes[2].plot(Θ, abs.(rperp), alpha=0.6, zorder=0,
        color="darkgrey", ls="dotted", label=L"$|r_{\perp}|$")

    # set the axis ticks and ranges
    axes[1].set(xlim=[0, 60.0], ylim=[1e-5, 0.9])
    axes[1].set_xticks([0, 15.0, 30.0, 45.0, 60.0])
    axes[2].set(xlim=[0, 60.0], ylim=[1e-5, 0.9])
    axes[2].set_xticks([0, 15.0, 30.0, 45.0, 60.0])
    [rax.set_ylim([0.0, 1.1]) for rax in raxes]
    [ax.set_yscale("log") for ax in axes]

    # labels
    axes[1].set_xlabel("Incident Angle at Surface [deg]")
    axes[1].set_ylabel("Event Density")
    raxes[1].set_ylabel("Fresnel Coeff.")
    axes[2].set_xlabel("Incident Angle at Ice [deg]")
    axes[2].set_ylabel("Event Density")
    raxes[2].set_ylabel("Fresnel Coeff.")

    # legends
    axes[1].legend(loc="lower left")
    raxes[2].legend(loc="upper left")

    # grids
    [ax.set_axisbelow(true) for ax in axes]
    [ax.grid(which="both", ls="dashed") for ax in axes]
    fig.tight_layout()

end

"""
    plot_polarization_angle(AΩ::Acceptance)

Plot histograms of polarization angles for direct and reflected events.

Arguments:
- `AΩ::Acceptance`: An object containing the acceptance calculation results.
"""
function plot_polarization_angle(AΩ::Union{Acceptance, OldAcceptance})

    # create the figure - plot this in polar coordinates
    fig = plt.figure(figsize=(6, 3))
    axes = [fig.add_subplot(1, 2, 1, polar=true),
        fig.add_subplot(1, 2, 2, polar=true)]

    # an array of colors
    colors = [colordirect, colorrefl]

    # axes[1].hist(AΩ.direct.θpol .|> deg2rad, density=true, color=colordirect, histtype="step")
    # axes[2].hist(AΩ.reflected.θpol .|> deg2rad, density=true, color=colorrefl, histtype="step")

    axes[1].hist(AΩ.reflected.θpol .|> deg2rad, density=true, color=colorrefl, histtype="step")
    axes[2].hist(AΩ.reflected.θpolsub .|> deg2rad, density=true, color="navy", histtype="step")

    # for (i, events) = enumerate([AΩ.direct, AΩ.reflected])

    #     # compute the polarization angle - 0 == H-Pol, 90 = V-Pol
    #     χ = deg2rad.(polarization_angle(events))

    #     # histogram the polarization angle
    #     # axes[i].hist(χ[trig.(events) .== true], bins=45, density=true, color=colors[i], histtype="step")
    #     axes[i].hist(χ, density=true, color=colors[i], histtype="step")

    # end
    # end loop over events

    # and make it pretty
    [ax.set(xlabel="Polarization Angle [deg]") for ax in axes]
    [ax.set_xticklabels(["+H", "", "+V", "", "-H", "", "-V", ""]) for ax in axes]
    axes[1].set_title("Top of Ice")
    axes[2].set_title("Bottom of Ice")
    # fig.suptitle("UH MC")
    fig.tight_layout()

    return fig, axes

end

"""
    plot_offaxis_angle(AΩ::Acceptance)

Generate a plot showing the off-axis angle distribution of direct and reflected events in the CoRaLS simulations.

Arguments:
- `AΩ::Acceptance`: An object containing the acceptance calculation results.
"""
function plot_offaxis_angle(AΩ::Union{Acceptance, OldAcceptance})

    # create the figure
    fig, ax = plt.subplots(figsize=(4, 4))

    # plot the histograms
    ax.hist(AΩ.reflected.ψ, histtype="step", color=colorrefl, alpha=0.3, label="Reflected (all)")
    ax.hist(AΩ.direct.ψ, histtype="step", color=colordirect, alpha=0.3, label="Direct (all)")
    ax.hist(AΩ.reflected.ψ[AΩ.reflected.triggered.==true],
        histtype="step", color=colorrefl, label="Reflected (trig)")
    ax.hist(AΩ.direct.ψ[AΩ.direct.triggered.==true],
        histtype="step", color=colordirect, label="Direct (trig)")

    # and make it pretty
    ax.set(xlabel="Off-Axis Angle [deg]", ylabel="Counts",
        xticks=[0.0, 30.0, 60.0, 90.0, 120.0, 150.0, 180.0],
        yscale="log", xlim=[0.0, 180.0])

    ax.legend(bbox_to_anchor=(0, 1.02, 1, 0.2), loc="lower left",
              mode="expand", borderaxespad=0, ncol=2)
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")

    return fig, ax
end

"""
    plot_acceptance(AΩ::Acceptance)

Plot the direct and reflected acceptance of the CoRaLS project.

Arguments:
- `AΩ::Acceptance`: An object containing the acceptance calculation results.

This function generates a line plot illustrating the acceptance rates for both direct and reflected events across different energy levels.
"""
function plot_acceptance(AΩ::Union{Acceptance, OldAcceptance}; min_count=0, ax=nothing, name="", 
                         colorrefl=colorrefl, colordirect=nothing, kwargs...)

    # create the figure
    if ax === nothing
        fig, ax = plt.subplots(figsize=(4, 4))
    else
        fig = gcf()
    end
    colordirect = isnothing(colordirect) ? colorrefl : colordirect

    if AΩ isa Acceptance
        rcount = AΩ.rcount
        dcount = AΩ.dcount 
    elseif AΩ isa OldAcceptance
        rcount = Int.(round.(AΩ.rAΩ ./ AΩ.gAΩ .* AΩ.ntrials))
        dcount = Int.(round.(AΩ.dAΩ ./ AΩ.gAΩ .* AΩ.ntrials))
    end
    ridx = rcount .>= min_count
    didx = dcount .>= min_count

    # Compute errors
    rerr = mcse.(rcount, AΩ.ntrials) .* AΩ.gAΩ ./ km^2 ./ sr
    derr = mcse.(dcount, AΩ.ntrials) .* AΩ.gAΩ ./ km^2 ./ sr

    # plot the lines
    logE = log10.(AΩ.energies * 1e18 / EeV)
    E = (logE[1:end-1] + logE[2:end]) / 2  # Bin centers in logspace
    ax.plot(E[ridx], AΩ.rAΩ[ridx] ./ km^2 ./ sr, 
        color=colorrefl, label=name*" Reflected", drawstyle="steps-mid"; kwargs...)
    ax.errorbar(
            E[ridx], AΩ.rAΩ[ridx] ./ km^2 ./ sr, yerr=rerr[ridx], 
            fmt="none", ecolor=colorrefl, alpha=0.5, capsize=2, elinewidth=1,
        )

    if haskey(kwargs, :linestyle) || haskey(kwargs, :ls)
        ax.plot(E[didx], AΩ.dAΩ[didx] ./ km^2 ./ sr,
            color=colordirect, label=name*" Direct", drawstyle="steps-mid"; kwargs...)
    else
        ax.plot(E[didx], AΩ.dAΩ[didx] ./ km^2 ./ sr,
            color=colordirect, label=name*" Direct", drawstyle="steps-mid", linestyle="--"; kwargs...)
    end
    ax.errorbar(
            E[didx], AΩ.dAΩ[didx] ./ km^2 ./ sr, yerr=derr[didx],
            fmt="none", ecolor=colordirect, alpha=0.5, capsize=2, elinewidth=1,
        )
    
    # and finally pretty it up
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")
    ax.set(xlabel=L"Energy [$\log_{10}$(eV)]",
        ylabel=L"Acceptance [km$^2$ sr]",
        yscale="log")
    ax.legend()
    ax.set_xlim(logE[1], logE[end])
    ax[:xaxis][:set_major_formatter](PyPlot.matplotlib.ticker.StrMethodFormatter("{x:.4g}"))
    return fig, ax

end

"""
    plot_differential_spectrum(AΩ::Acceptance, T::Unitful.Time)

Create a plot of the differential event spectrum of CoRaLS events given an acceptance object and a mission duration.

Arguments:
- `AΩ::Acceptance`: An object containing the acceptance calculation results.
- `T::Unitful.Time`: The total observation time of the experiment.

This function plots the estimated number of events per energy bin over the duration of the mission.
"""
function plot_differential_spectrum(AΩ::Union{Acceptance, OldAcceptance}, T)

    # create the figure
    fig, ax = plt.subplots(figsize=(4, 4))

    # the fraction of events that have PSRs
    psr_frac = psr_fraction(AΩ.altitude)

    # calculate the differential event spectrums
    # devents = differential_spectrum(AΩ.energies, psr_frac*AΩ.dAΩ, 2yr)
    revents = differential_spectrum(AΩ.energies, (30_000 / 20_000) * psr_frac * AΩ.rAΩ, 2yr)
    PIRevents = differential_spectrum(AΩ.energies, 6.3 * psr_frac * AΩ.rAΩ, 2yr)

    # plot the direct events
    # ax.hist(18.0 .+ log10.(AΩ.energies[1:end-1] ./ 1.0EeV),
    #          18.0 .+ log10.(AΩ.energies ./ 1.0EeV),
    #          weights=devents,
    #          histtype="step", color="#F88D07",
    #         label="Direct, N=$(Int(round(sum(devents))))")
    ax.hist(18.0 .+ log10.(AΩ.energies[1:end-1] ./ 1.0EeV),
        18.0 .+ log10.(AΩ.energies ./ 1.0EeV),
        weights=PIRevents,
        histtype="step", color="#F88D07",
        label="PIRs, N=$(Int(round(sum(PIRevents))))")

    # plot the reflected events
    ax.hist(18.0 .+ log10.(AΩ.energies[1:end-1] ./ 1.0EeV),
        18.0 .+ log10.(AΩ.energies ./ 1.0EeV),
        weights=revents,
        histtype="step", color="#0772F8",
        label="PSRs, N=$(Int(round(sum(revents))))")

    # and finally pretty it up
    ax.set_axisbelow(true)
    ax.grid(which="both", linestyle="dashed")
    ax.set(xlabel=L"Energy [$\log_{10}$(eV)]",
        xlim=[18.0 + log10(AΩ.energies[1] / EeV), 18.0 + log10(AΩ.energies[end] / EeV)],
        ylabel="Events per bin / $(T)",
        xticks=[17.0, 17.5, 18.0, 18.5, 19.0, 19.5, 20.0, 20.5, 21.0],
        yscale="log")
    ax.legend()

    return fig, ax

end

"""
    mcse(count, ntrials, [gAΩ])

Compute the Monte Carlo Standard Error as count or acceptance if gAΩ is given.

Uses Monte Carlo Standard Error for binomial proportion: sqrt(p*(1-p)/n)
"""
function mcse(count::Real, ntrials::Integer)
    p = count ./ ntrials
    return sqrt.(p .* (1 .- p) ./ ntrials)
end

"""
    plot_rate_experiment(AΩ::Acceptance, xs; xlabel)

Plot the counts in events / yr for an array of acceptance runs computed at
some variable x (e.g., ice depths, altitudes).

Arguments:
- AΩs : Array of acceptance objects.
- xs : Array of paramter values for the x axis corresponding to each AΩ.

This function plots the estimated number of events per energy bin over the duration of the mission.
"""
function plot_rate_experiment(AΩs, xs; xlabel="")
    xs = ustrip.(xs)

    r_spectra = [differential_spectrum(AΩ.energies, AΩ.rAΩ, 1yr) for AΩ in AΩs]
    r_spectra_mat = hcat(r_spectra...)'  # matrix [xs, energies]
    d_spectra = [differential_spectrum(AΩ.energies, AΩ.dAΩ, 1yr) for AΩ in AΩs]
    d_spectra_mat = hcat(d_spectra...)'

    r_mcse = [differential_spectrum(AΩ.energies, mcse.(AΩ.rcount, AΩ.ntrials)*AΩ.gAΩ, 1yr) for AΩ in AΩs]
    r_mcse_mat = hcat(r_mcse...)'  # matrix [xs, energies]
    # d_mcse = [differential_spectrum(AΩ.energies, mcse(AΩ.dcount, AΩ.ntrials, AΩ.gAΩ), 1yr) for AΩ in AΩs]
    # d_mcse_mat = hcat(d_mcse...)'

    cmap = PyPlot.get_cmap("plasma")

    # Energy array (take bin centers in log space)
    logE = log10.(AΩs[1].energies * 1e18 / EeV)
    E = (logE[1:end-1] + logE[2:end]) / 2
    nE = length(E)

    fig, ax = plt.subplots(figsize=(5, 4))

    # Show total and label it at its peak
    total_r = sum(r_spectra_mat; dims=2)
    ax.plot(xs, total_r, color="k", linewidth=2, alpha=0.6)
    y_peak_r = maximum(total_r)
    idx_peak_r = argmax(total_r)
    x_peak_r = xs[idx_peak_r]
    ha = x_peak_r == xs[1] ? "left" : "center"
    ax.text(x_peak_r, y_peak_r, "Total", color="k", fontsize=10, va="bottom", ha=ha)

    # Show each energy curve
    for i in 1:nE
        norm_val = (i-1) / nE
        color = cmap(norm_val)
        ax.plot(xs, r_spectra_mat[:, i], color=color, alpha=0.8)
        ax.errorbar(
            xs, r_spectra_mat[:, i], yerr=r_mcse_mat[:, i],
            fmt="none", ecolor=color, alpha=0.5, capsize=2, elinewidth=1,
        )
        ax.plot(xs, d_spectra_mat[:, i], color=color, alpha=0.8, linestyle="--")

        # Find peak for reflected spectrum
        y_peak_r = maximum(r_spectra_mat[:, i])
        idx_peak_r = argmax(r_spectra_mat[:, i])
        x_peak_r = xs[idx_peak_r]
        ha = x_peak_r == xs[1] ? "left" : "center"
        ha = x_peak_r == xs[end] ? "right" : ha
        ax.text(
            x_peak_r, y_peak_r * 1.01,  # Slightly offset in y
            "10^"*string(round(E[i], sigdigits=4))*" eV", 
            color=color, fontsize=8, va="bottom", ha=ha
        )
    end
    # Show colorbar
    sm = PyPlot.matplotlib.cm.ScalarMappable(
        cmap=cmap, 
        norm=PyPlot.matplotlib.colors.Normalize(vmin=logE[1], vmax=logE[end])
    )
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, pad=0.02)
    cbar.set_label("log₁₀(Energy [eV])")
    cbar.set_ticks(logE)
    cbar.set_ticklabels(string.(round.(logE, sigdigits=3)))

    # Plot params
    ax.set(
        xlabel=xlabel,
        ylabel="Events / yr",
        axisbelow=true,
        ylim=0,
        xlim=[xs[1], xs[end]]
    )
    ax.grid(true, which="major", linestyle="dashed")
    ax.tick_params(
        axis="both",
        which="major",
        direction="in",
        top=true,
        right=true
    )
    plt.xscale("log")

    return fig, ax
end

function plot_event_outcomes(A::Acceptance; title="")
    # Get CR energy bin centers in log space
    logE = log10.(A.energies * 1e18 / EeV)
    E = (logE[1:end-1] + logE[2:end]) / 2

    # Outcomes are accepted counts, nottriggered, and failtypes (vs E)
    # Not triggered = (Ntrials - accepted - failed) by process of elimination
    outcomes = vcat("Accepted", "NotTriggered", String.(Symbol.(A.failtypes))...)
    d_not_triggered = (A.ntrials .- A.dcount .- sum(A.dfailed;dims=2))
    r_not_triggered = (A.ntrials .- A.rcount .- sum(A.rfailed;dims=2))
    d_outcomes = hcat(A.dcount, d_not_triggered, A.dfailed)
    r_outcomes = hcat(A.rcount, r_not_triggered, A.rfailed)
    
    colors = PyPlot.cm.get_cmap("tab10")(0:length(outcomes))
    markers = [".", "x", "v", "^", "<", ">", "+", "D", "p", "h"]

    fig, ax = subplots(1, 1, figsize=(5, 4))
    for (j, outcome) in enumerate(outcomes)
        color = colors[j,:]
        marker = markers[mod1(j, length(markers))]
        y_d = d_outcomes[:, j] ./ A.ntrials
        y_r = r_outcomes[:, j] ./ A.ntrials
        y_rmcse = mcse.(r_outcomes[:, j],  A.ntrials)
        
        # Plot each line only if events exist
        if !all(iszero, y_r) 
            ax.plot(E, y_r, color=color, marker=marker, label="$(outcome)", alpha=0.8, drawstyle="steps-mid")
            ax.errorbar(
                    E, y_r, yerr=y_rmcse,fmt="none", ecolor=color, 
                    alpha=0.5, capsize=2, elinewidth=1
                )
        end
        !all(iszero, y_d) && ax[:plot](E, y_d, color=color, marker=marker, linestyle="--",  alpha=0.8, drawstyle="steps-mid")
    end

    # Axes setup
    ax.legend(loc="best", title="(Refl: —, Direct: ---)", title_fontsize=8, fontsize=8, framealpha=0.3)
    ax.set(
        title=title,
        xlabel="log₁₀(Energy [eV])",
        ylabel="P (Event Outcome / Ntrials)",
        axisbelow=true,
        xlim=(logE[1], logE[end]),
        ylim=(5e-7,2),
        yscale="log"
    )
    ax[:xaxis][:set_major_formatter](PyPlot.matplotlib.ticker.StrMethodFormatter("{x:.4g}"))
    ax.grid(true, which="major", linestyle="dashed")
    ax.tick_params(
        axis="both",
        which="major",
        direction="in",
        top=true,
        right=true
    )
    return fig, ax
end