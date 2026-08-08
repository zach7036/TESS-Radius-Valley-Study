def center_uncertainty_mc(
    base_pairs: pd.DataFrame,
    planets: pd.DataFrame,
    b: int,
    batch_size: int,
    seed: int,
    width: float = DEFAULT_WIDTH,
    pea_tol: float = DEFAULT_PEA_TOL,
) -> pd.DataFrame:
    """Propagate published uncertainties in the M/K/G valley centers."""
    x = base_pairs.reset_index(drop=True)
    pool = make_pool(x, planets).reset_index(drop=True)
    key_to_idx = {k: i for i, k in enumerate(pool["planet_key"])}
    edge_a = x["inner_key"].map(key_to_idx).to_numpy(int)
    edge_b = x["outer_key"].map(key_to_idx).to_numpy(int)
    edge_types = x["stype"].to_numpy(str)
    r0 = pool["Planet Radius Value"].to_numpy(float)
    pool_types = pool["stype"].to_numpy(str)
    type_indices = {t: np.where(pool_types == t)[0] for t in TYPE_ORDER}
    pea = np.abs(
        np.log(
            x["outer_radius"].to_numpy(float)
            / x["inner_radius"].to_numpy(float)
        )
    ) < np.log(1 + pea_tol)
    rng = np.random.default_rng(seed)
    logthr = np.log(1 + pea_tol)
    rows = []
    done = 0
    while done < b:
        m = min(batch_size, b - done)
        sampled = {
            t: rng.normal(CENTERS[t], CENTER_SIGMA[t], size=m) for t in TYPE_ORDER
        }
        centers_mat = np.empty((m, len(pool)), dtype=float)
        for t, idx in type_indices.items():
            centers_mat[:, idx] = sampled[t][:, None]
        rmat = np.broadcast_to(r0, (m, len(pool)))
        vmat = np.abs(rmat - centers_mat) <= width
        pN, pV = _dynamic_null_probabilities(rmat, vmat, type_indices, logthr)
        val = vmat[:, edge_a] | vmat[:, edge_b]
        qN = np.empty((m, len(x)), dtype=float)
        qV = np.empty((m, len(x)), dtype=float)
        for j, t in enumerate(edge_types):
            qN[:, j] = pN[t]
            qV[:, j] = pV[t]
        nv = val.sum(axis=1)
        nn = (~val).sum(axis=1)
        kv = (val & pea[None, :]).sum(axis=1)
        kn = ((~val) & pea[None, :]).sum(axis=1)
        qv = np.where(val, qV, 0).sum(axis=1)
        qn = np.where(~val, qN, 0).sum(axis=1)
        valid = (nv > 0) & (nn > 0) & np.isfinite(qv) & np.isfinite(qn)
        obs_v = np.divide(kv, nv, out=np.full(m, np.nan), where=nv > 0)
        obs_n = np.divide(kn, nn, out=np.full(m, np.nan), where=nn > 0)
        null_v = np.divide(qv, nv, out=np.full(m, np.nan), where=nv > 0)
        null_n = np.divide(qn, nn, out=np.full(m, np.nan), where=nn > 0)
        frame = {
            "valid": valid,
            "n_valley": nv,
            "raw_difference": obs_n - obs_v,
            "contrast": (obs_n - null_n) - (obs_v - null_v),
        }
        for t in TYPE_ORDER:
            frame[f"center_{t}"] = sampled[t]
        rows.append(pd.DataFrame(frame))
        done += m
    return pd.concat(rows, ignore_index=True)


def robustness_grid(pairs: pd.DataFrame, planets: pd.DataFrame) -> pd.DataFrame:
    """Targeted one-factor-at-a-time and membership-grid robustness checks."""
    specs = []
    # Core two-dimensional membership grid used in Figure 3.
    for method in ["nominal", "overlap"]:
        for fmax in [0.10, 0.15, 0.20, 0.25, np.inf]:
            for width in [0.08, 0.10, 0.12, 0.15, 0.20]:
                specs.append(dict(method=method, fmax=fmax, width=width, pea_tol=0.10, rmax=4.0, snrmin=None, disposition="all"))
        # One-factor perturbations around the prespecified primary analysis.
        for pea_tol in [0.05, 0.10, 0.15, 0.20]:
            specs.append(dict(method=method, fmax=0.20, width=0.10, pea_tol=pea_tol, rmax=4.0, snrmin=None, disposition="all"))
        for rmax in [3.0, 4.0, 6.0]:
            specs.append(dict(method=method, fmax=0.20, width=0.10, pea_tol=0.10, rmax=rmax, snrmin=None, disposition="all"))
        for snrmin in [7.1, 10.0, 15.0]:
            specs.append(dict(method=method, fmax=0.20, width=0.10, pea_tol=0.10, rmax=4.0, snrmin=snrmin, disposition="all"))
        for disp in ["all", "confirmed", "pc_only"]:
            specs.append(dict(method=method, fmax=0.20, width=0.10, pea_tol=0.10, rmax=4.0, snrmin=None, disposition=disp))
    # Remove duplicate primary rows while preserving order.
    seen, unique_specs = set(), []
    for d in specs:
        key = tuple((k, d[k]) for k in sorted(d))
        if key not in seen:
            seen.add(key); unique_specs.append(d)
    rows = []
    for d in unique_specs:
        r = evaluate(pairs, planets, **d)
        rows.append(r.summary)
    return pd.DataFrame(rows)


def center_scan(
    base_pairs: pd.DataFrame,
    planets: pd.DataFrame,
    shifts: np.ndarray,
    width: float = DEFAULT_WIDTH,
    pea_tol: float = DEFAULT_PEA_TOL,
) -> pd.DataFrame:
    rows = []
    for method in ["nominal", "overlap"]:
        for d in shifts:
            centers = {t: CENTERS[t] + float(d) for t in TYPE_ORDER}
            r = evaluate_selected(
                base_pairs,
                planets,
                method=method,
                width=width,
                pea_tol=pea_tol,
                center_override=centers,
            )
            s = r.summary.copy()
            s["center_shift"] = float(d)
            rows.append(s)
    return pd.DataFrame(rows)


def widget_tests(result: EvalResult) -> pd.DataFrame:
    """Additive, multiplicative, and logit transfer of non-valley enhancement."""
    s = result.summary
    v = result.pairs[result.pairs["valley_pair"]]
    k = int(v["peas"].sum())
    q = v["null_probability"].to_numpy(float)
    obs_n = s["nonvalley_peas_rate"]
    null_n = s["nonvalley_null_rate"]
    delta = obs_n - null_n
    factor = obs_n / null_n if null_n > 0 else float("nan")
    logit_shift = logit(obs_n) - logit(null_n)
    models = {
        "additive": np.clip(q + delta, 0, 1),
        "multiplicative": np.clip(q * factor, 0, 1),
        "logit": expit(logit(q) + logit_shift),
    }
    rows = []
    for name, probs in models.items():
        rows.append(
            {
                "method": s["method"],
                "model": name,
                "n_valley_pairs": len(v),
                "observed_valley_peas": k,
                "expected_valley_peas": float(probs.sum()),
                "expected_rate": float(probs.mean()),
                "one_sided_poisson_binomial_p": poisson_binomial_cdf(k, probs),
            }
        )
    return pd.DataFrame(rows)


def nonadjacent_pairs(
    planets: pd.DataFrame,
    fmax: float = DEFAULT_FERR,
    rmin: float = DEFAULT_RMIN,
    rmax: float = DEFAULT_RMAX,
) -> pd.DataFrame:
    rows = []
    for tic, g in planets.groupby("TIC", sort=True):
        g = g.sort_values("Orbital Period (days) Value").reset_index(drop=True)
        for i in range(len(g)):
            for j in range(i + 2, len(g)):
                a, b = g.iloc[i], g.iloc[j]
                if not (rmin <= a["Planet Radius Value"] <= rmax and rmin <= b["Planet Radius Value"] <= rmax):
                    continue
                if not (
                    np.isfinite(a["frac_radius_error"])
                    and np.isfinite(b["frac_radius_error"])
                    and a["frac_radius_error"] <= fmax
                    and b["frac_radius_error"] <= fmax
                ):
                    continue
                rows.append(
                    {
                        "TIC": int(tic),
                        "stype": a["stype"],
                        "inner_key": a["planet_key"],
                        "outer_key": b["planet_key"],
                        "inner_toi": a["Full TOI ID"],
                        "outer_toi": b["Full TOI ID"],
                        "inner_period": a["Orbital Period (days) Value"],
                        "outer_period": b["Orbital Period (days) Value"],
                        "period_ratio": b["Orbital Period (days) Value"] / a["Orbital Period (days) Value"],
                        "inner_radius": a["Planet Radius Value"],
                        "outer_radius": b["Planet Radius Value"],
                        "inner_error": a["Planet Radius Uncertainty"],
                        "outer_error": b["Planet Radius Uncertainty"],
                        "inner_frac_error": a["frac_radius_error"],
                        "outer_frac_error": b["frac_radius_error"],
                        "inner_snr": a["Signal-to-noise"],
                        "outer_snr": b["Signal-to-noise"],
                        "inner_disposition": a["TOI Disposition"],
                        "outer_disposition": b["TOI Disposition"],
                        "valley_center": a["valley_center"],
                        "abs_log_radius_ratio": abs(math.log(b["Planet Radius Value"] / a["Planet Radius Value"])),
                        "radius_ratio": b["Planet Radius Value"] / a["Planet Radius Value"],
                        "peas_10pct": abs(math.log(b["Planet Radius Value"] / a["Planet Radius Value"])) < math.log(1.1),
                    }
                )
    return pd.DataFrame(rows)


def exact_fisher_power(nv: int, nn: int, pv: float, pn: float, alpha: float = 0.05) -> float:
    """Exact independent-binomial power for the one-sided Fisher test."""
    power = 0.0
    for kv in range(nv + 1):
        pkv = binom.pmf(kv, nv, pv)
        for kn in range(nn + 1):
            pkn = binom.pmf(kn, nn, pn)
            p = fisher_exact([[kv, nv-kv], [kn, nn-kn]], alternative="less").pvalue
            if p < alpha:
                power += pkv * pkn
    return float(power)


def power_curve(nv: int, nn: int, pn: float) -> pd.DataFrame:
    vals = np.linspace(0, pn, 101)
    return pd.DataFrame({
        "p_valley": vals,
        "p_nonvalley": pn,
        "power": [exact_fisher_power(nv, nn, float(pv), pn) for pv in vals],
    })


def qsummary(values: pd.Series | np.ndarray) -> dict:
    a = np.asarray(values, dtype=float)
    a = a[np.isfinite(a)]
    if len(a) == 0:
        return {"n": 0}
    return {
        "n": int(len(a)),
        "mean": float(np.mean(a)),
        "median": float(np.median(a)),
        "q025": float(np.quantile(a, 0.025)),
        "q975": float(np.quantile(a, 0.975)),
        "prob_gt_zero": float(np.mean(a > 0)),
    }

