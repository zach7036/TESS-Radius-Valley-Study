def save_figures(
    nominal: EvalResult,
    overlap: EvalResult,
    robust: pd.DataFrame,
    me_mc: pd.DataFrame,
    scan: pd.DataFrame,
) -> None:
    # Figure 1: which planets are admitted by interval overlap.
    pool = nominal.pool[[
        "planet_key", "stype", "Planet Radius Value", "Planet Radius Uncertainty",
        "analysis_center", "valley"
    ]].copy().rename(columns={"valley": "nominal"})
    ov = overlap.pool[["planet_key", "valley"]].rename(columns={"valley": "overlap"})
    pool = pool.merge(ov, on="planet_key", how="left")
    pool["distance"] = np.abs(pool["Planet Radius Value"] - pool["analysis_center"])
    pool["overlap_only"] = pool["overlap"] & ~pool["nominal"]
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    for t, marker in zip(TYPE_ORDER, ["o", "s", "^"]):
        g = pool[pool["stype"] == t]
        ax.scatter(
            g["Planet Radius Uncertainty"],
            g["distance"],
            marker=marker,
            alpha=0.65,
            label=f"{t} dwarfs",
        )
    lim = max(pool["Planet Radius Uncertainty"].max(), pool["distance"].max()) * 1.03
    xx = np.linspace(0, lim, 200)
    ax.axhline(DEFAULT_WIDTH, linestyle="--", linewidth=1.2, label="nominal window")
    ax.plot(xx, DEFAULT_WIDTH + xx, linestyle=":", linewidth=1.5, label="overlap boundary")
    oo = pool[pool["overlap_only"]]
    ax.scatter(
        oo["Planet Radius Uncertainty"],
        oo["distance"],
        facecolors="none",
        edgecolors="black",
        s=85,
        linewidths=1.1,
        label="overlap-only label",
    )
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)
    ax.set_xlabel("Reported 1-sigma radius uncertainty (Earth radii)")
    ax.set_ylabel("Absolute distance from spectral-type valley center (Earth radii)")
    ax.set_title("Interval overlap expands the effective valley window")
    ax.legend(fontsize=8, ncol=2)
    fig.tight_layout()
    fig.savefig(FIG / "figure1_overlap_geometry.png", dpi=220)
    plt.close(fig)

    # Figure 2: observed and null rates for both classification methods.
    labels, obs, low, high, null = [], [], [], [], []
    for r, name in [(nominal, "Nominal"), (overlap, "1-sigma overlap")]:
        s = r.summary
        labels += [f"{name}\nvalley", f"{name}\nnon-valley"]
        obs += [s["valley_peas_rate"], s["nonvalley_peas_rate"]]
        low += [s["valley_wilson_low"], s["nonvalley_wilson_low"]]
        high += [s["valley_wilson_high"], s["nonvalley_wilson_high"]]
        null += [s["valley_null_rate"], s["nonvalley_null_rate"]]
    x = np.arange(len(labels))
    yerr = np.vstack([np.array(obs) - np.array(low), np.array(high) - np.array(obs)])
    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    ax.errorbar(x, obs, yerr=yerr, fmt="o", capsize=4, label="observed (95% Wilson CI)")
    ax.scatter(x, null, marker="x", s=65, label="type-stratified random-pair null")
    ax.set_xticks(x, labels)
    ax.set_ylim(0, max(high) * 1.18)
    ax.set_ylabel("Fraction of adjacent pairs within 10% in radius")
    ax.set_title("The apparent valley deficit depends on membership rule")
    ax.legend(fontsize=9)
    fig.tight_layout()
    fig.savefig(FIG / "figure2_rates.png", dpi=220)
    plt.close(fig)

    # Figure 3: robustness at 10% pea tolerance, all dispositions, rmax=4, no SNR threshold.
    sub = robust[
        (robust["pea_tolerance"] == 0.10)
        & (robust["disposition"] == "all")
        & (robust["rmax"] == 4.0)
        & (robust["snrmin"].isna())
        & (robust["fmax"].isin([0.10, 0.15, 0.20, 0.25]))
    ].copy()
    fig, ax = plt.subplots(figsize=(8.0, 5.2))
    for method, linestyle in [("nominal", "-"), ("overlap", "--")]:
        for fmax, marker in zip([0.10, 0.15, 0.20, 0.25], ["o", "s", "^", "D"]):
            g = sub[(sub["method"] == method) & (sub["fmax"] == fmax)].sort_values("width_rearth")
            ax.plot(
                g["width_rearth"],
                g["enhancement_contrast_non_minus_valley"],
                marker=marker,
                linestyle=linestyle,
                label=f"{method}, frac. error <= {int(fmax*100)}%",
            )
    ax.axhline(0, linewidth=1)
    ax.set_xlabel("Valley half-width (Earth radii)")
    ax.set_ylabel("Null-adjusted similarity contrast")
    ax.set_title("Membership-rule sensitivity across radius-quality cuts")
    ax.legend(fontsize=7.5, ncol=2)
    fig.tight_layout()
    fig.savefig(FIG / "figure3_robustness.png", dpi=220)
    plt.close(fig)

    # Figure 4: measurement-error propagated contrast.
    valid = me_mc[me_mc["valid"]]
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    bins = np.linspace(
        np.nanquantile(valid["contrast"], 0.005),
        np.nanquantile(valid["contrast"], 0.995),
        55,
    )
    for dist in ["normal", "lognormal", "normal+center"]:
        g = valid[valid["distribution"] == dist]
        ax.hist(g["contrast"], bins=bins, density=True, histtype="step", linewidth=1.6, label=dist)
    ax.axvline(0, linewidth=1)
    ax.axvline(nominal.summary["enhancement_contrast_non_minus_valley"], linestyle="--", linewidth=1.2, label="nominal point estimate")
    ax.axvline(overlap.summary["enhancement_contrast_non_minus_valley"], linestyle=":", linewidth=1.4, label="overlap point estimate")
    ax.set_xlabel("Null-adjusted similarity contrast")
    ax.set_ylabel("Monte Carlo density")
    ax.set_title("Measurement-error propagation does not yield a decisive TESS contrast")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(FIG / "figure4_measurement_error_mc.png", dpi=220)
    plt.close(fig)

    # Figure 5: placebo center scan.
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    for method, linestyle in [("nominal", "-"), ("overlap", "--")]:
        g = scan[scan["method"] == method].sort_values("center_shift")
        ax.plot(g["center_shift"], g["enhancement_contrast_non_minus_valley"], linestyle=linestyle, label=method)
    ax.axvline(0, linewidth=1, linestyle=":")
    ax.axhline(0, linewidth=1)
    ax.set_xlabel("Common shift applied to M/K/G valley centers (Earth radii)")
    ax.set_ylabel("Null-adjusted similarity contrast")
    ax.set_title("The overlap-rule contrast is not sharply localized to the published valley centers")
    ax.legend()
    fig.tight_layout()
    fig.savefig(FIG / "figure5_center_scan.png", dpi=220)
    plt.close(fig)

