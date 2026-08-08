def evaluate(
    pairs: pd.DataFrame,
    planets: pd.DataFrame,
    fmax: float = DEFAULT_FERR,
    rmin: float = DEFAULT_RMIN,
    rmax: float = DEFAULT_RMAX,
    snrmin: Optional[float] = None,
    disposition: Literal["all", "confirmed", "pc_only"] = "all",
    method: Literal["nominal", "overlap", "certain"] = "nominal",
    width: float = DEFAULT_WIDTH,
    pea_tol: float = DEFAULT_PEA_TOL,
    center_override: Optional[dict] = None,
) -> EvalResult:
    x = select_pairs(pairs, fmax, rmin, rmax, snrmin, disposition)
    out = evaluate_selected(x, planets, method, width, pea_tol, center_override)
    out.summary.update(
        {
            "fmax": fmax,
            "rmin": rmin,
            "rmax": rmax,
            "snrmin": snrmin,
            "disposition": disposition,
        }
    )
    return out


def cluster_bootstrap(result: EvalResult, b: int, seed: int) -> pd.DataFrame:
    x = result.pairs.copy()
    systems = np.array(sorted(x["TIC"].unique()))
    agg = []
    for tic in systems:
        g = x[x["TIC"] == tic]
        v = g[g["valley_pair"]]
        n = g[~g["valley_pair"]]
        agg.append(
            [
                len(v),
                v["peas"].sum(),
                v["null_probability"].sum(),
                len(n),
                n["peas"].sum(),
                n["null_probability"].sum(),
            ]
        )
    agg = np.asarray(agg, dtype=float)
    rng = np.random.default_rng(seed)
    # Resample systems with replacement by multinomial weights.
    w = rng.multinomial(len(systems), np.full(len(systems), 1 / len(systems)), size=b)
    totals = w @ agg
    nv, kv, qv, nn, kn, qn = totals.T
    valid = (nv > 0) & (nn > 0)
    obs_v = np.divide(kv, nv, out=np.full_like(kv, np.nan), where=nv > 0)
    obs_n = np.divide(kn, nn, out=np.full_like(kn, np.nan), where=nn > 0)
    null_v = np.divide(qv, nv, out=np.full_like(qv, np.nan), where=nv > 0)
    null_n = np.divide(qn, nn, out=np.full_like(qn, np.nan), where=nn > 0)
    raw = obs_n - obs_v
    contrast = (obs_n - null_n) - (obs_v - null_v)
    return pd.DataFrame(
        {
            "valid": valid,
            "raw_difference": raw,
            "contrast": contrast,
            "n_valley": nv,
            "n_nonvalley": nn,
        }
    )


def leave_one_system_out(result: EvalResult, planets: pd.DataFrame) -> pd.DataFrame:
    rows = []
    x = result.pairs
    method = result.summary["method"]
    width = result.summary["width_rearth"]
    pea_tol = result.summary["pea_tolerance"]
    for tic in sorted(x["TIC"].unique()):
        y = x[x["TIC"] != tic].copy()
        r = evaluate_selected(y, planets, method=method, width=width, pea_tol=pea_tol)
        rows.append(
            {
                "left_out_TIC": tic,
                "n_pairs": r.summary["n_pairs"],
                "n_valley_pairs": r.summary["n_valley_pairs"],
                "contrast": r.summary["enhancement_contrast_non_minus_valley"],
                "raw_difference": r.summary["raw_rate_difference_non_minus_valley"],
                "fisher_p": r.summary["fisher_one_sided_p"],
            }
        )
    return pd.DataFrame(rows)


def permutation_test(
    result: EvalResult,
    b: int,
    batch_size: int,
    seed: int,
) -> pd.DataFrame:
    """Shuffle radius/error records among planet nodes within spectral type.

    This preserves the global type-specific radius distribution, uncertainty
    distribution, pair graph, and multiplicity while destroying within-system
    radius arrangement.
    """
    x = result.pairs.reset_index(drop=True)
    pool = result.pool.reset_index(drop=True)
    method = result.summary["method"]
    width = result.summary["width_rearth"]
    pea_tol = result.summary["pea_tolerance"]
    key_to_idx = {k: i for i, k in enumerate(pool["planet_key"])}
    edge_a = x["inner_key"].map(key_to_idx).to_numpy(int)
    edge_b = x["outer_key"].map(key_to_idx).to_numpy(int)
    edge_type = x["stype"].to_numpy(str)
    type_idx = {t: np.where(pool["stype"].to_numpy(str) == t)[0] for t in TYPE_ORDER}
    r0 = pool["Planet Radius Value"].to_numpy(float)
    e0 = pool["Planet Radius Uncertainty"].to_numpy(float)
    c0 = pool["analysis_center"].to_numpy(float)
    v0 = classify_planets(r0, e0, c0, width, method)
    # Null match probabilities are invariant to permutation because each type's
    # multiset of radius records and labels is unchanged.
    pN = result.p_non_by_type
    pV = result.p_val_by_type
    qN_edge = np.array([pN[t] for t in edge_type], dtype=float)
    qV_edge = np.array([pV[t] for t in edge_type], dtype=float)
    rng = np.random.default_rng(seed)
    rows = []
    done = 0
    logthr = np.log(1 + pea_tol)
    while done < b:
        m = min(batch_size, b - done)
        rmat = np.empty((m, len(pool)), dtype=np.float64)
        vmat = np.empty((m, len(pool)), dtype=bool)
        for t, idx in type_idx.items():
            if len(idx) == 0:
                continue
            order = np.argsort(rng.random((m, len(idx))), axis=1)
            rmat[:, idx] = r0[idx][order]
            vmat[:, idx] = v0[idx][order]
        pea = np.abs(np.log(rmat[:, edge_b] / rmat[:, edge_a])) < logthr
        val = vmat[:, edge_a] | vmat[:, edge_b]
        nv = val.sum(axis=1)
        nn = (~val).sum(axis=1)
        kv = (pea & val).sum(axis=1)
        kn = (pea & ~val).sum(axis=1)
        qv = (val * qV_edge[None, :]).sum(axis=1)
        qn = ((~val) * qN_edge[None, :]).sum(axis=1)
        valid = (nv > 0) & (nn > 0)
        obs_v = np.divide(kv, nv, out=np.full(m, np.nan), where=nv > 0)
        obs_n = np.divide(kn, nn, out=np.full(m, np.nan), where=nn > 0)
        null_v = np.divide(qv, nv, out=np.full(m, np.nan), where=nv > 0)
        null_n = np.divide(qn, nn, out=np.full(m, np.nan), where=nn > 0)
        rows.append(
            pd.DataFrame(
                {
                    "valid": valid,
                    "raw_difference": obs_n - obs_v,
                    "contrast": (obs_n - null_n) - (obs_v - null_v),
                    "n_valley": nv,
                }
            )
        )
        done += m
    return pd.concat(rows, ignore_index=True)


def _dynamic_null_probabilities(
    rmat: np.ndarray,
    vmat: np.ndarray,
    type_indices: dict[str, np.ndarray],
    logthr: float,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray]]:
    """Per-draw distinct-object null probabilities for each spectral type."""
    pN, pV = {}, {}
    for t, idx in type_indices.items():
        if len(idx) == 0:
            continue
        rt = rmat[:, idx]
        vt = vmat[:, idx]
        m, n = rt.shape
        sim = np.abs(np.log(rt[:, :, None] / rt[:, None, :])) < logthr
        diag = np.eye(n, dtype=bool)[None, :, :]
        sim &= ~diag
        non = ~vt
        mask_nn = non[:, :, None] & non[:, None, :] & ~diag
        den_n = mask_nn.sum(axis=(1, 2))
        num_n = (sim & mask_nn).sum(axis=(1, 2))
        pN[t] = np.divide(
            num_n,
            den_n,
            out=np.full(m, np.nan, dtype=float),
            where=den_n > 0,
        )
        mask_va = vt[:, :, None] & ~diag
        den_v = mask_va.sum(axis=(1, 2))
        num_v = (sim & mask_va).sum(axis=(1, 2))
        pV[t] = np.divide(
            num_v,
            den_v,
            out=np.full(m, np.nan, dtype=float),
            where=den_v > 0,
        )
    return pN, pV


def measurement_error_mc(
    nominal_result: EvalResult,
    b: int,
    batch_size: int,
    seed: int,
    distribution: Literal["normal", "lognormal"] = "normal",
    sample_centers: bool = False,
) -> pd.DataFrame:
    """Propagate reported radius errors and reclassify latent radii each draw."""
    x = nominal_result.pairs.reset_index(drop=True)
    pool = nominal_result.pool.reset_index(drop=True)
    width = nominal_result.summary["width_rearth"]
    pea_tol = nominal_result.summary["pea_tolerance"]
    key_to_idx = {k: i for i, k in enumerate(pool["planet_key"])}
    edge_a = x["inner_key"].map(key_to_idx).to_numpy(int)
    edge_b = x["outer_key"].map(key_to_idx).to_numpy(int)
    edge_types = x["stype"].to_numpy(str)
    r0 = pool["Planet Radius Value"].to_numpy(float)
    s0 = pool["Planet Radius Uncertainty"].to_numpy(float)
    c0 = pool["analysis_center"].to_numpy(float)
    pool_types = pool["stype"].to_numpy(str)
    type_indices = {t: np.where(pool_types == t)[0] for t in TYPE_ORDER}
    rng = np.random.default_rng(seed)
    logthr = np.log(1 + pea_tol)
    rows = []
    done = 0
    while done < b:
        m = min(batch_size, b - done)
        if distribution == "normal":
            rmat = rng.normal(r0, s0, size=(m, len(pool)))
            # With f<=0.2, negatives are virtually absent; resample any that occur.
            bad = rmat <= 0
            while bad.any():
                rmat[bad] = rng.normal(
                    np.broadcast_to(r0, rmat.shape)[bad],
                    np.broadcast_to(s0, rmat.shape)[bad],
                )
                bad = rmat <= 0
        elif distribution == "lognormal":
            frac = s0 / r0
            sigma_ln = np.sqrt(np.log1p(frac * frac))
            mu_ln = np.log(r0) - 0.5 * sigma_ln * sigma_ln
            rmat = rng.lognormal(mu_ln, sigma_ln, size=(m, len(pool)))
        else:
            raise ValueError(distribution)
        if sample_centers:
            centers_mat = np.empty((m, len(pool)), dtype=float)
            for t, idx in type_indices.items():
                draw_t = rng.normal(CENTERS[t], CENTER_SIGMA[t], size=m)
                centers_mat[:, idx] = draw_t[:, None]
        else:
            centers_mat = np.broadcast_to(c0, (m, len(pool)))
        vmat = np.abs(rmat - centers_mat) <= width
        pea = np.abs(np.log(rmat[:, edge_b] / rmat[:, edge_a])) < logthr
        val = vmat[:, edge_a] | vmat[:, edge_b]
        pN, pV = _dynamic_null_probabilities(rmat, vmat, type_indices, logthr)
        qN = np.empty((m, len(x)), dtype=float)
        qV = np.empty((m, len(x)), dtype=float)
        for j, t in enumerate(edge_types):
            qN[:, j] = pN[t]
            qV[:, j] = pV[t]
        nv = val.sum(axis=1)
        nn = (~val).sum(axis=1)
        kv = (pea & val).sum(axis=1)
        kn = (pea & ~val).sum(axis=1)
        qv = np.where(val, qV, 0.0).sum(axis=1)
        qn = np.where(~val, qN, 0.0).sum(axis=1)
        valid = (nv > 0) & (nn > 0) & np.isfinite(qv) & np.isfinite(qn)
        obs_v = np.divide(kv, nv, out=np.full(m, np.nan), where=nv > 0)
        obs_n = np.divide(kn, nn, out=np.full(m, np.nan), where=nn > 0)
        null_v = np.divide(qv, nv, out=np.full(m, np.nan), where=nv > 0)
        null_n = np.divide(qn, nn, out=np.full(m, np.nan), where=nn > 0)
        rows.append(
            pd.DataFrame(
                {
                    "valid": valid,
                    "distribution": distribution + ("+center" if sample_centers else ""),
                    "n_valley": nv,
                    "n_valley_peas": kv,
                    "valley_peas_rate": obs_v,
                    "n_nonvalley": nn,
                    "n_nonvalley_peas": kn,
                    "nonvalley_peas_rate": obs_n,
                    "raw_difference": obs_n - obs_v,
                    "contrast": (obs_n - null_n) - (obs_v - null_v),
                }
            )
        )
        done += m
    return pd.concat(rows, ignore_index=True)

