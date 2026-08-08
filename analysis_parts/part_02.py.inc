def select_pairs(
    pairs: pd.DataFrame,
    fmax: float = DEFAULT_FERR,
    rmin: float = DEFAULT_RMIN,
    rmax: float = DEFAULT_RMAX,
    snrmin: Optional[float] = None,
    disposition: Literal["all", "confirmed", "pc_only"] = "all",
) -> pd.DataFrame:
    x = pairs[
        pairs["inner_radius"].between(rmin, rmax)
        & pairs["outer_radius"].between(rmin, rmax)
    ].copy()
    if np.isfinite(fmax):
        x = x[
            np.isfinite(x["inner_frac_error"])
            & np.isfinite(x["outer_frac_error"])
            & (x["inner_frac_error"] <= fmax)
            & (x["outer_frac_error"] <= fmax)
        ].copy()
    if snrmin is not None:
        x = x[
            np.isfinite(x["inner_snr"])
            & np.isfinite(x["outer_snr"])
            & (x["inner_snr"] >= snrmin)
            & (x["outer_snr"] >= snrmin)
        ].copy()
    if disposition == "confirmed":
        x = x[
            x["inner_disposition"].isin(["CP", "KP"])
            & x["outer_disposition"].isin(["CP", "KP"])
        ].copy()
    elif disposition == "pc_only":
        x = x[
            (x["inner_disposition"] == "PC") & (x["outer_disposition"] == "PC")
        ].copy()
    return x


def classify_planets(
    radii: np.ndarray,
    errors: np.ndarray,
    centers: np.ndarray,
    width: float,
    method: Literal["nominal", "overlap", "certain"],
) -> np.ndarray:
    if method == "nominal":
        return np.abs(radii - centers) <= width
    if method == "overlap":
        return (
            np.isfinite(errors)
            & (radii + errors >= centers - width)
            & (radii - errors <= centers + width)
        )
    if method == "certain":
        return (
            np.isfinite(errors)
            & (radii - errors >= centers - width)
            & (radii + errors <= centers + width)
        )
    raise ValueError(method)


def make_pool(x: pd.DataFrame, planets: pd.DataFrame) -> pd.DataFrame:
    keys = pd.unique(
        pd.concat([x["inner_key"], x["outer_key"]], ignore_index=True)
    )
    pool = planets[planets["planet_key"].isin(keys)].copy()
    # Deterministic ordering is useful for simulations.
    return pool.sort_values(["stype", "TIC", "Orbital Period (days) Value"]).reset_index(
        drop=True
    )


def distinct_match_probability(
    radii_a: np.ndarray,
    ids_a: np.ndarray,
    radii_b: np.ndarray,
    ids_b: np.ndarray,
    pea_tol: float,
) -> float:
    if len(radii_a) == 0 or len(radii_b) == 0:
        return float("nan")
    valid = ids_a[:, None] != ids_b[None, :]
    if not valid.any():
        return float("nan")
    similar = np.abs(np.log(radii_a[:, None] / radii_b[None, :])) < np.log(
        1 + pea_tol
    )
    return float(similar[valid].mean())


@dataclass
class EvalResult:
    summary: dict
    pairs: pd.DataFrame
    pool: pd.DataFrame
    p_non_by_type: dict
    p_val_by_type: dict


def evaluate_selected(
    x: pd.DataFrame,
    planets: pd.DataFrame,
    method: Literal["nominal", "overlap", "certain"] = "nominal",
    width: float = DEFAULT_WIDTH,
    pea_tol: float = DEFAULT_PEA_TOL,
    center_override: Optional[dict] = None,
) -> EvalResult:
    x = x.copy()
    centers_pair = x["stype"].map(center_override or CENTERS).to_numpy(float)
    x["analysis_center"] = centers_pair
    va = classify_planets(
        x["inner_radius"].to_numpy(float),
        x["inner_error"].to_numpy(float),
        centers_pair,
        width,
        method,
    )
    vb = classify_planets(
        x["outer_radius"].to_numpy(float),
        x["outer_error"].to_numpy(float),
        centers_pair,
        width,
        method,
    )
    x["inner_valley"] = va
    x["outer_valley"] = vb
    x["valley_pair"] = va | vb
    x["peas"] = x["abs_log_radius_ratio"] < np.log(1 + pea_tol)

    pool = make_pool(x, planets)
    pool_centers = pool["stype"].map(center_override or CENTERS).to_numpy(float)
    pool["analysis_center"] = pool_centers
    pool["valley"] = classify_planets(
        pool["Planet Radius Value"].to_numpy(float),
        pool["Planet Radius Uncertainty"].to_numpy(float),
        pool_centers,
        width,
        method,
    )

    p_non, p_val = {}, {}
    for t, g in pool.groupby("stype", observed=True):
        non = g[~g["valley"]]
        val = g[g["valley"]]
        allg = g
        p_non[t] = distinct_match_probability(
            non["Planet Radius Value"].to_numpy(float),
            non["planet_key"].to_numpy(str),
            non["Planet Radius Value"].to_numpy(float),
            non["planet_key"].to_numpy(str),
            pea_tol,
        )
        p_val[t] = distinct_match_probability(
            val["Planet Radius Value"].to_numpy(float),
            val["planet_key"].to_numpy(str),
            allg["Planet Radius Value"].to_numpy(float),
            allg["planet_key"].to_numpy(str),
            pea_tol,
        )

    x["null_probability"] = [
        p_val[t] if is_v else p_non[t]
        for t, is_v in zip(x["stype"], x["valley_pair"])
    ]
    v = x[x["valley_pair"]]
    n = x[~x["valley_pair"]]
    kv, kn = int(v["peas"].sum()), int(n["peas"].sum())
    nv, nn = len(v), len(n)
    obs_v = kv / nv if nv else float("nan")
    obs_n = kn / nn if nn else float("nan")
    null_v = v["null_probability"].mean() if nv else float("nan")
    null_n = n["null_probability"].mean() if nn else float("nan")
    enhancement_v = obs_v - null_v if nv else float("nan")
    enhancement_n = obs_n - null_n if nn else float("nan")
    contrast = enhancement_n - enhancement_v if nv and nn else float("nan")
    raw_diff = obs_n - obs_v if nv and nn else float("nan")
    fisher_p = (
        fisher_exact([[kv, nv - kv], [kn, nn - kn]], alternative="less").pvalue
        if nv and nn
        else float("nan")
    )
    mw_p = (
        mannwhitneyu(
            v["abs_log_radius_ratio"],
            n["abs_log_radius_ratio"],
            alternative="greater",
        ).pvalue
        if nv and nn
        else float("nan")
    )
    ci_v = wilson(kv, nv)
    ci_n = wilson(kn, nn)
    summary = {
        "method": method,
        "width_rearth": width,
        "pea_tolerance": pea_tol,
        "n_pairs": int(len(x)),
        "n_systems": int(x["TIC"].nunique()),
        "n_pool_planets": int(len(pool)),
        "n_pool_valley_planets": int(pool["valley"].sum()),
        "n_valley_pairs": nv,
        "n_valley_peas": kv,
        "valley_peas_rate": obs_v,
        "valley_wilson_low": ci_v[0],
        "valley_wilson_high": ci_v[1],
        "valley_null_rate": null_v,
        "valley_enhancement": enhancement_v,
        "n_nonvalley_pairs": nn,
        "n_nonvalley_peas": kn,
        "nonvalley_peas_rate": obs_n,
        "nonvalley_wilson_low": ci_n[0],
        "nonvalley_wilson_high": ci_n[1],
        "nonvalley_null_rate": null_n,
        "nonvalley_enhancement": enhancement_n,
        "raw_rate_difference_non_minus_valley": raw_diff,
        "enhancement_contrast_non_minus_valley": contrast,
        "fisher_one_sided_p": float(fisher_p),
        "mannwhitney_one_sided_p": float(mw_p),
        "median_abs_log_ratio_valley": float(v["abs_log_radius_ratio"].median())
        if nv
        else float("nan"),
        "median_abs_log_ratio_nonvalley": float(n["abs_log_radius_ratio"].median())
        if nn
        else float("nan"),
        "n_valley_pairs_by_type": {
            str(k): int(vv)
            for k, vv in v.groupby("stype", observed=True).size().to_dict().items()
        },
        "n_nonvalley_pairs_by_type": {
            str(k): int(vv)
            for k, vv in n.groupby("stype", observed=True).size().to_dict().items()
        },
        "null_non_by_type": {str(k): float(vv) for k, vv in p_non.items()},
        "null_valley_by_type": {str(k): float(vv) for k, vv in p_val.items()},
    }
    return EvalResult(summary, x, pool, p_non, p_val)

