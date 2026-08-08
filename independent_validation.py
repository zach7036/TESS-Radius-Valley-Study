#!/usr/bin/env python3
"""Independent, deliberately simple validation of the primary result.

This script does not import tess_valley_analysis.py. It reconstructs the main
hard labels, pair counts, pea rates, Fisher tests, and analytic expected latent
valley count from the exported primary pair table.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, norm

ROOT = Path(__file__).resolve().parent
P = ROOT / "tables" / "primary_adjacent_pairs.csv"
CENTERS = {"M": 1.64, "K": 1.75, "G": 1.86}
WIDTH = 0.10

x = pd.read_csv(P)
center = x["stype"].map(CENTERS).to_numpy(float)
ri = x["inner_radius"].to_numpy(float)
ro = x["outer_radius"].to_numpy(float)
si = x["inner_error"].to_numpy(float)
so = x["outer_error"].to_numpy(float)
peas = np.abs(np.log(ro / ri)) < math.log(1.10)

nom_i = np.abs(ri-center) <= WIDTH
nom_o = np.abs(ro-center) <= WIDTH
ov_i = (ri+si >= center-WIDTH) & (ri-si <= center+WIDTH)
ov_o = (ro+so >= center-WIDTH) & (ro-so <= center+WIDTH)

p_i = norm.cdf((center+WIDTH-ri)/si)-norm.cdf((center-WIDTH-ri)/si)
p_o = norm.cdf((center+WIDTH-ro)/so)-norm.cdf((center-WIDTH-ro)/so)
p_pair = 1-(1-p_i)*(1-p_o)

out = {"source": str(P), "n_pairs": int(len(x)), "duplicate_edges": int(x.duplicated(["TIC","inner_key","outer_key"]).sum())}
for name, v in {"nominal": nom_i|nom_o, "overlap": ov_i|ov_o}.items():
    nv = int(v.sum()); nn = int((~v).sum())
    kv = int((peas&v).sum()); kn = int((peas&~v).sum())
    out[name] = {
        "n_valley_pairs": nv,
        "n_valley_peas": kv,
        "valley_rate": kv/nv,
        "n_nonvalley_pairs": nn,
        "n_nonvalley_peas": kn,
        "nonvalley_rate": kn/nn,
        "raw_difference": kn/nn-kv/nv,
        "fisher_one_sided_p": float(fisher_exact([[kv,nv-kv],[kn,nn-kn]], alternative="less").pvalue),
    }
out["analytic_expected_latent_valley_pairs"] = float(p_pair.sum())
out["all_period_ratios_gt_one"] = bool((x["period_ratio"] > 1).all())
out["all_frac_errors_within_primary_cut"] = bool(((x["inner_frac_error"] <= .20) & (x["outer_frac_error"] <= .20)).all())
out["all_radii_within_primary_range"] = bool((x[["inner_radius","outer_radius"]].min(axis=1) >= .5).all() and (x[["inner_radius","outer_radius"]].max(axis=1) <= 4).all())

summary = json.load(open(ROOT / "summary.json"))
checks = {
    "nominal_counts_match": out["nominal"]["n_valley_pairs"] == summary["nominal"]["n_valley_pairs"],
    "overlap_counts_match": out["overlap"]["n_valley_pairs"] == summary["overlap"]["n_valley_pairs"],
    "nominal_fisher_matches": math.isclose(out["nominal"]["fisher_one_sided_p"], summary["nominal"]["fisher_one_sided_p"], rel_tol=0, abs_tol=1e-15),
    "overlap_fisher_matches": math.isclose(out["overlap"]["fisher_one_sided_p"], summary["overlap"]["fisher_one_sided_p"], rel_tol=0, abs_tol=1e-15),
    "analytic_expected_matches": math.isclose(out["analytic_expected_latent_valley_pairs"], summary["classification_inflation"]["expected_latent_valley_pairs"], rel_tol=0, abs_tol=1e-12),
}
out["checks"] = checks
out["all_checks_pass"] = bool(all(checks.values()))

path = ROOT / "independent_validation.json"
path.write_text(json.dumps(out, indent=2))
print(json.dumps(out, indent=2))
