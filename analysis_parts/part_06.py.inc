def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", type=Path, default=ROOT / "tois_mast.csv")
    ap.add_argument("--bootstrap", type=int, default=100_000)
    ap.add_argument("--permutations", type=int, default=100_000)
    ap.add_argument("--measurement-mc", type=int, default=20_000)
    ap.add_argument("--center-mc", type=int, default=50_000)
    ap.add_argument("--batch-size", type=int, default=500)
    args = ap.parse_args()

    raw, planets, pairs, meta = load_data(args.data)
    data_hash = sha256(args.data)
    meta.update(
        {
            "sha256": data_hash,
            "raw_rows": int(len(raw)),
            "eligible_multiplanet_rows": int(len(planets)),
            "eligible_systems": int(planets["TIC"].nunique()),
            "all_adjacent_pairs_before_radius_quality": int(len(pairs)),
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
            "scipy": scipy.__version__,
        }
    )

    base = select_pairs(pairs)
    nominal = evaluate_selected(base, planets, "nominal", DEFAULT_WIDTH, DEFAULT_PEA_TOL)
    overlap = evaluate_selected(base, planets, "overlap", DEFAULT_WIDTH, DEFAULT_PEA_TOL)
    certain = evaluate_selected(base, planets, "certain", DEFAULT_WIDTH, DEFAULT_PEA_TOL)
    for r in [nominal, overlap, certain]:
        r.summary.update(
            {
                "fmax": DEFAULT_FERR,
                "rmin": DEFAULT_RMIN,
                "rmax": DEFAULT_RMAX,
                "snrmin": None,
                "disposition": "all",
            }
        )

    # Decompose overlap into nominal and overlap-only pair labels.
    comp = nominal.pairs[["TIC", "inner_key", "outer_key", "peas", "valley_pair"]].rename(columns={"valley_pair": "nominal_valley"})
    comp = comp.merge(
        overlap.pairs[["TIC", "inner_key", "outer_key", "valley_pair"]].rename(columns={"valley_pair": "overlap_valley"}),
        on=["TIC", "inner_key", "outer_key"],
        how="left",
    )
    comp["classification_group"] = np.select(
        [comp["nominal_valley"], comp["overlap_valley"]],
        ["nominal valley", "overlap-only"],
        default="neither",
    )
    decomposition = (
        comp.groupby("classification_group", observed=True)
        .agg(n_pairs=("peas", "size"), n_peas=("peas", "sum"))
        .reset_index()
    )
    decomposition["peas_rate"] = decomposition["n_peas"] / decomposition["n_pairs"]

    poolcomp = nominal.pool[["planet_key", "stype", "Planet Radius Value", "Planet Radius Uncertainty", "analysis_center", "valley"]].rename(columns={"valley": "nominal_valley"})
    poolcomp = poolcomp.merge(overlap.pool[["planet_key", "valley"]].rename(columns={"valley": "overlap_valley"}), on="planet_key")
    poolcomp["distance_from_center"] = np.abs(poolcomp["Planet Radius Value"] - poolcomp["analysis_center"])
    poolcomp["overlap_only"] = poolcomp["overlap_valley"] & ~poolcomp["nominal_valley"]
    # Probability that the latent radius lies in the fixed valley window under
    # a Gaussian approximation to the reported 1-sigma uncertainty.
    z_hi = (poolcomp["analysis_center"] + DEFAULT_WIDTH - poolcomp["Planet Radius Value"]) / poolcomp["Planet Radius Uncertainty"]
    z_lo = (poolcomp["analysis_center"] - DEFAULT_WIDTH - poolcomp["Planet Radius Value"]) / poolcomp["Planet Radius Uncertainty"]
    poolcomp["latent_valley_probability"] = norm.cdf(z_hi) - norm.cdf(z_lo)
    pmap = poolcomp.set_index("planet_key")["latent_valley_probability"]
    comp["inner_latent_probability"] = comp["inner_key"].map(pmap)
    comp["outer_latent_probability"] = comp["outer_key"].map(pmap)
    comp["latent_valley_pair_probability"] = 1 - (1-comp["inner_latent_probability"])*(1-comp["outer_latent_probability"])

    # Main resampling and falsification analyses.
    boot_nom = cluster_bootstrap(nominal, args.bootstrap, SEED + 1)
    boot_ov = cluster_bootstrap(overlap, args.bootstrap, SEED + 2)
    loo_nom = leave_one_system_out(nominal, planets)
    loo_ov = leave_one_system_out(overlap, planets)
    perm_nom = permutation_test(nominal, args.permutations, 2500, SEED + 3)
    perm_ov = permutation_test(overlap, args.permutations, 2500, SEED + 4)
    mc_normal = measurement_error_mc(nominal, args.measurement_mc, args.batch_size, SEED + 5, "normal")
    mc_lognormal = measurement_error_mc(nominal, args.measurement_mc, args.batch_size, SEED + 6, "lognormal")
    mc_normal_center = measurement_error_mc(nominal, args.measurement_mc, args.batch_size, SEED + 8, "normal", sample_centers=True)
    me_mc = pd.concat([mc_normal, mc_lognormal, mc_normal_center], ignore_index=True)
    center_mc = center_uncertainty_mc(base, planets, args.center_mc, args.batch_size, SEED + 7)
    robust = robustness_grid(pairs, planets)
    scan = center_scan(base, planets, np.round(np.arange(-0.60, 0.6001, 0.02), 2))
    widgets = pd.concat([widget_tests(nominal), widget_tests(overlap)], ignore_index=True)

    # Non-adjacent placebo pairs. These are deliberately analyzed with the
    # same system-cluster bootstrap and within-type radius-label permutation
    # used for adjacent pairs, because pairs from the same system are not
    # independent.
    nonadj = nonadjacent_pairs(planets)
    if len(nonadj):
        nonadj_nom = evaluate_selected(nonadj, planets, "nominal", DEFAULT_WIDTH, DEFAULT_PEA_TOL)
        nonadj_ov = evaluate_selected(nonadj, planets, "overlap", DEFAULT_WIDTH, DEFAULT_PEA_TOL)
        nonadj_boot_nom = cluster_bootstrap(nonadj_nom, args.bootstrap, SEED + 11)
        nonadj_boot_ov = cluster_bootstrap(nonadj_ov, args.bootstrap, SEED + 12)
        nonadj_perm_nom = permutation_test(nonadj_nom, args.permutations, 2500, SEED + 13)
        nonadj_perm_ov = permutation_test(nonadj_ov, args.permutations, 2500, SEED + 14)
    else:
        nonadj_nom = nonadj_ov = None
        nonadj_boot_nom = nonadj_boot_ov = None
        nonadj_perm_nom = nonadj_perm_ov = None

    # Exact permutation p values.
    for label, res, perm in [("nominal", nominal, perm_nom), ("overlap", overlap, perm_ov)]:
        valid = perm[perm["valid"]]
        obs_c = res.summary["enhancement_contrast_non_minus_valley"]
        obs_d = res.summary["raw_rate_difference_non_minus_valley"]
        res.summary["permutation_p_contrast"] = float((1 + (valid["contrast"] >= obs_c).sum()) / (1 + len(valid)))
        res.summary["permutation_p_raw_difference"] = float((1 + (valid["raw_difference"] >= obs_d).sum()) / (1 + len(valid)))

    # Add bootstrap summaries and leave-one-out ranges.
    for res, boot, loo in [(nominal, boot_nom, loo_nom), (overlap, boot_ov, loo_ov)]:
        vb = boot[boot["valid"]]
        res.summary["cluster_bootstrap_contrast"] = qsummary(vb["contrast"])
        res.summary["cluster_bootstrap_raw_difference"] = qsummary(vb["raw_difference"])
        res.summary["leave_one_system_out_contrast_min"] = float(loo["contrast"].min())
        res.summary["leave_one_system_out_contrast_max"] = float(loo["contrast"].max())
        res.summary["leave_one_system_out_fisher_p_min"] = float(loo["fisher_p"].min())
        res.summary["leave_one_system_out_fisher_p_max"] = float(loo["fisher_p"].max())

    # Correlation-aware summaries for the non-adjacent placebo.
    if nonadj_nom is not None:
        for res, boot, perm in [
            (nonadj_nom, nonadj_boot_nom, nonadj_perm_nom),
            (nonadj_ov, nonadj_boot_ov, nonadj_perm_ov),
        ]:
            vb = boot[boot["valid"]]
            vp = perm[perm["valid"]]
            obs_c = res.summary["enhancement_contrast_non_minus_valley"]
            obs_d = res.summary["raw_rate_difference_non_minus_valley"]
            res.summary["cluster_bootstrap_contrast"] = qsummary(vb["contrast"])
            res.summary["cluster_bootstrap_raw_difference"] = qsummary(vb["raw_difference"])
            res.summary["permutation_p_contrast"] = float(
                (1 + (vp["contrast"] >= obs_c).sum()) / (1 + len(vp))
            )
            res.summary["permutation_p_raw_difference"] = float(
                (1 + (vp["raw_difference"] >= obs_d).sum()) / (1 + len(vp))
            )

    # Measurement- and center-uncertainty summaries.
    mc_summary = {}
    for dist in ["normal", "lognormal", "normal+center"]:
        g = me_mc[(me_mc["distribution"] == dist) & me_mc["valid"]]
        mc_summary[dist] = {
            "contrast": qsummary(g["contrast"]),
            "raw_difference": qsummary(g["raw_difference"]),
            "valley_peas_rate": qsummary(g["valley_peas_rate"]),
            "nonvalley_peas_rate": qsummary(g["nonvalley_peas_rate"]),
            "n_valley": qsummary(g["n_valley"]),
        }
    center_valid = center_mc[center_mc["valid"]]
    center_summary = {
        "contrast": qsummary(center_valid["contrast"]),
        "raw_difference": qsummary(center_valid["raw_difference"]),
        "n_valley": qsummary(center_valid["n_valley"]),
    }

    overlap_only = poolcomp[poolcomp["overlap_only"]]
    inflation_summary = {
        "nominal_valley_pairs": nominal.summary["n_valley_pairs"],
        "overlap_valley_pairs": overlap.summary["n_valley_pairs"],
        "pair_inflation_factor": overlap.summary["n_valley_pairs"] / max(nominal.summary["n_valley_pairs"], 1),
        "nominal_valley_planets": nominal.summary["n_pool_valley_planets"],
        "overlap_valley_planets": overlap.summary["n_pool_valley_planets"],
        "planet_inflation_factor": overlap.summary["n_pool_valley_planets"] / max(nominal.summary["n_pool_valley_planets"], 1),
        "overlap_only_planets": int(len(overlap_only)),
        "overlap_only_distance_median": float(overlap_only["distance_from_center"].median()),
        "overlap_only_distance_min": float(overlap_only["distance_from_center"].min()),
        "overlap_only_distance_max": float(overlap_only["distance_from_center"].max()),
        "overlap_only_frac_error_median": float((overlap_only["Planet Radius Uncertainty"] / overlap_only["Planet Radius Value"]).median()),
        "expected_latent_valley_planets": float(poolcomp["latent_valley_probability"].sum()),
        "expected_latent_valley_pairs": float(comp["latent_valley_pair_probability"].sum()),
        "overlap_label_precision_expected": float(comp.loc[comp["overlap_valley"], "latent_valley_pair_probability"].sum() / comp["overlap_valley"].sum()),
        "nominal_label_precision_expected": float(comp.loc[comp["nominal_valley"], "latent_valley_pair_probability"].sum() / comp["nominal_valley"].sum()),
        "overlap_label_expected_sensitivity": float(comp.loc[comp["overlap_valley"], "latent_valley_pair_probability"].sum() / comp["latent_valley_pair_probability"].sum()),
        "nominal_label_expected_sensitivity": float(comp.loc[comp["nominal_valley"], "latent_valley_pair_probability"].sum() / comp["latent_valley_pair_probability"].sum()),
    }

    # Center-scan specificity diagnostics.
    scan_summary = {}
    for method in ["nominal", "overlap"]:
        g = scan[scan["method"] == method].dropna(subset=["enhancement_contrast_non_minus_valley"])
        at0 = float(g.loc[np.isclose(g["center_shift"], 0), "enhancement_contrast_non_minus_valley"].iloc[0])
        scan_summary[method] = {
            "contrast_at_published_centers": at0,
            "max_contrast": float(g["enhancement_contrast_non_minus_valley"].max()),
            "shift_at_max": float(g.loc[g["enhancement_contrast_non_minus_valley"].idxmax(), "center_shift"]),
            "percentile_of_published_center": float((g["enhancement_contrast_non_minus_valley"] <= at0).mean()),
            "fraction_positive": float((g["enhancement_contrast_non_minus_valley"] > 0).mean()),
        }

    power_nominal = power_curve(nominal.summary["n_valley_pairs"], nominal.summary["n_nonvalley_pairs"], nominal.summary["nonvalley_peas_rate"])
    power_overlap = power_curve(overlap.summary["n_valley_pairs"], overlap.summary["n_nonvalley_pairs"], overlap.summary["nonvalley_peas_rate"])
    power_summary = {
        "nominal_max_power_if_valley_rate_zero": float(power_nominal.loc[power_nominal["p_valley"] == 0, "power"].iloc[0]),
        "nominal_power_at_observed_rates": exact_fisher_power(nominal.summary["n_valley_pairs"], nominal.summary["n_nonvalley_pairs"], nominal.summary["valley_peas_rate"], nominal.summary["nonvalley_peas_rate"]),
        "overlap_max_power_if_valley_rate_zero": float(power_overlap.loc[power_overlap["p_valley"] == 0, "power"].iloc[0]),
        "overlap_power_at_observed_rates": exact_fisher_power(overlap.summary["n_valley_pairs"], overlap.summary["n_nonvalley_pairs"], overlap.summary["valley_peas_rate"], overlap.summary["nonvalley_peas_rate"]),
    }

    all_summary = {
        "metadata": meta,
        "sample": {
            "selected_pairs": int(len(base)),
            "selected_systems": int(base["TIC"].nunique()),
            "selected_unique_planets": int(len(make_pool(base, planets))),
            "pairs_by_type": {str(k): int(v) for k, v in base.groupby("stype", observed=True).size().to_dict().items()},
            "fractional_radius_error_quantiles": {str(k): float(v) for k, v in planets["frac_radius_error"].replace([np.inf, -np.inf], np.nan).quantile([0, .1, .25, .5, .75, .9, .95, 1]).to_dict().items()},
        },
        "nominal": nominal.summary,
        "overlap": overlap.summary,
        "certain": certain.summary,
        "classification_inflation": inflation_summary,
        "measurement_error_monte_carlo": mc_summary,
        "center_uncertainty_monte_carlo": center_summary,
        "center_scan": scan_summary,
        "binomial_power_diagnostic": power_summary,
        "nonadjacent_nominal": nonadj_nom.summary if nonadj_nom else None,
        "nonadjacent_overlap": nonadj_ov.summary if nonadj_ov else None,
    }

    # Save tables and raw simulation summaries.
    clean_cols = [
        "TIC", "Full TOI ID", "TOI Disposition", "Orbital Period (days) Value",
        "Planet Radius Value", "Planet Radius Uncertainty", "frac_radius_error",
        "Effective Temperature Value", "Surface Gravity Value", "Star Radius Value",
        "Signal-to-noise", "stype", "valley_center", "planet_key"
    ]
    planets[clean_cols].to_csv(TAB / "eligible_multiplanet_planets.csv", index=False)
    base.to_csv(TAB / "primary_adjacent_pairs.csv", index=False)
    nominal.pairs.to_csv(TAB / "primary_pairs_nominal_classification.csv", index=False)
    overlap.pairs.to_csv(TAB / "primary_pairs_overlap_classification.csv", index=False)
    poolcomp.to_csv(TAB / "planet_classification_comparison.csv", index=False)
    decomposition.to_csv(TAB / "overlap_decomposition.csv", index=False)
    comp.to_csv(TAB / "pair_classification_probability_comparison.csv", index=False)
    power_nominal.to_csv(TAB / "power_curve_nominal.csv", index=False)
    power_overlap.to_csv(TAB / "power_curve_overlap.csv", index=False)
    robust.to_csv(TAB / "robustness_grid.csv", index=False)
    scan.to_csv(TAB / "center_placebo_scan.csv", index=False)
    widgets.to_csv(TAB / "widget_tests.csv", index=False)
    loo_nom.to_csv(TAB / "leave_one_system_out_nominal.csv", index=False)
    loo_ov.to_csv(TAB / "leave_one_system_out_overlap.csv", index=False)
    boot_nom.to_csv(TAB / "cluster_bootstrap_nominal.csv.gz", index=False, compression="gzip")
    boot_ov.to_csv(TAB / "cluster_bootstrap_overlap.csv.gz", index=False, compression="gzip")
    perm_nom.to_csv(TAB / "permutation_nominal.csv.gz", index=False, compression="gzip")
    perm_ov.to_csv(TAB / "permutation_overlap.csv.gz", index=False, compression="gzip")
    me_mc.to_csv(TAB / "measurement_error_mc.csv.gz", index=False, compression="gzip")
    center_mc.to_csv(TAB / "center_uncertainty_mc.csv.gz", index=False, compression="gzip")
    if len(nonadj):
        nonadj.to_csv(TAB / "nonadjacent_placebo_pairs.csv", index=False)
        nonadj_boot_nom.to_csv(TAB / "nonadjacent_cluster_bootstrap_nominal.csv.gz", index=False, compression="gzip")
        nonadj_boot_ov.to_csv(TAB / "nonadjacent_cluster_bootstrap_overlap.csv.gz", index=False, compression="gzip")
        nonadj_perm_nom.to_csv(TAB / "nonadjacent_permutation_nominal.csv.gz", index=False, compression="gzip")
        nonadj_perm_ov.to_csv(TAB / "nonadjacent_permutation_overlap.csv.gz", index=False, compression="gzip")

    with (ROOT / "summary.json").open("w", encoding="utf-8") as f:
        json.dump(all_summary, f, indent=2, allow_nan=True)

    save_figures(nominal, overlap, robust, me_mc, scan)

    # Concise machine-readable text summary.
    print(json.dumps(all_summary, indent=2, allow_nan=True))


if __name__ == "__main__":
    main()
