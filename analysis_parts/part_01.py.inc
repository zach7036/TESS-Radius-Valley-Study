#!/usr/bin/env python3
"""Independent TESS test of radius-valley disruption of peas-in-a-pod size similarity.

This script uses the fixed MAST TOI CSV snapshot shipped alongside it. It creates
all cleaned tables, statistical summaries, robustness checks, simulations, and
figures used in the accompanying report.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, Literal, Optional, Sequence, Tuple

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy
from scipy.stats import fisher_exact, mannwhitneyu, norm, binom

ROOT = Path(__file__).resolve().parent
FIG = ROOT / "figures"
TAB = ROOT / "tables"
FIG.mkdir(exist_ok=True)
TAB.mkdir(exist_ok=True)

CENTERS = {"M": 1.64, "K": 1.75, "G": 1.86}
CENTER_SIGMA = {"M": 0.03, "K": 0.11, "G": 0.06}
TYPE_ORDER = ["M", "K", "G"]
DEFAULT_WIDTH = 0.10
DEFAULT_FERR = 0.20
DEFAULT_RMIN = 0.50
DEFAULT_RMAX = 4.00
DEFAULT_PEA_TOL = 0.10
SEED = 20260808


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def parse_header(path: Path) -> dict:
    lines = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("#"):
                break
            lines.append(line.strip("# \n"))
    out = {"header_lines": lines}
    for line in lines:
        if line.startswith("File Created:"):
            out["file_created"] = line.split(":", 1)[1].strip()
        elif line.startswith("Sectors:"):
            out["sectors"] = line.split(":", 1)[1].strip()
        elif line.startswith("Date of previous TOI list release"):
            out["previous_release"] = line.rsplit(" ", 1)[-1].strip()
    return out


def wilson(k: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n <= 0:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1 + z * z / n
    center = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return center - half, center + half


def poisson_binomial_cdf(k: int, probs: Sequence[float]) -> float:
    """Exact P(X <= k) for independent nonidentical Bernoulli trials."""
    probs = np.asarray(probs, dtype=float)
    if len(probs) == 0:
        return float("nan")
    pmf = np.array([1.0])
    for p in probs:
        pmf = np.convolve(pmf, np.array([1 - p, p]))
    return float(pmf[: k + 1].sum())


def logit(x: np.ndarray | float) -> np.ndarray | float:
    x = np.clip(x, 1e-9, 1 - 1e-9)
    return np.log(x / (1 - x))


def expit(x: np.ndarray | float) -> np.ndarray | float:
    return 1 / (1 + np.exp(-x))


def _numeric(df: pd.DataFrame, columns: Iterable[str]) -> None:
    for c in columns:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")


def load_data(path: Path) -> tuple[pd.DataFrame, pd.DataFrame, dict]:
    meta = parse_header(path)
    raw = pd.read_csv(path, comment="#")
    nums = [
        "TIC",
        "Orbital Period (days) Value",
        "Orbital Period Uncertainty",
        "Planet Radius Value",
        "Planet Radius Uncertainty",
        "Star Radius Value",
        "Star Radius Uncertainty",
        "Surface Gravity Value",
        "Surface Gravity Uncertainty",
        "Effective Temperature Value",
        "Effective Temperature Uncertainty",
        "Signal-to-noise",
        "Centroid Offset",
        "Planet Number",
    ]
    _numeric(raw, nums)

    planets = raw[raw["TOI Disposition"].isin(["PC", "CP", "KP"])].copy()
    planets = planets[
        np.isfinite(planets["Orbital Period (days) Value"])
        & (planets["Orbital Period (days) Value"] > 0)
        & np.isfinite(planets["Planet Radius Value"])
        & (planets["Planet Radius Value"] > 0)
    ].copy()
    # Main-sequence G/K/M host cuts. The temperature boundaries match the
    # TESS radius-valley study used for the adopted centers.
    planets = planets[
        planets["Effective Temperature Value"].between(2600, 6040, inclusive="left")
        & planets["Surface Gravity Value"].between(4.0, 5.5)
        & planets["Star Radius Value"].between(0.1, 1.5)
    ].copy()
    planets["stype"] = pd.cut(
        planets["Effective Temperature Value"],
        [-np.inf, 3880, 5340, 6040],
        labels=["M", "K", "G"],
        right=False,
    ).astype(str)
    planets["valley_center"] = planets["stype"].map(CENTERS).astype(float)
    planets["center_sigma"] = planets["stype"].map(CENTER_SIGMA).astype(float)
    planets["frac_radius_error"] = (
        planets["Planet Radius Uncertainty"] / planets["Planet Radius Value"]
    )
    planets["planet_key"] = planets["Full TOI ID"].astype(str)
    planets["TIC"] = planets["TIC"].astype("int64")

    # Retain systems with at least two eligible planet candidates/planets.
    multiplicity = planets.groupby("TIC").size()
    planets = planets[planets["TIC"].isin(multiplicity[multiplicity >= 2].index)].copy()

    pairs = []
    for tic, g in planets.groupby("TIC", sort=True):
        g = g.sort_values("Orbital Period (days) Value")
        for i in range(len(g) - 1):
            a, b = g.iloc[i], g.iloc[i + 1]
            pairs.append(
                {
                    "TIC": int(tic),
                    "stype": a["stype"],
                    "inner_key": a["planet_key"],
                    "outer_key": b["planet_key"],
                    "inner_toi": a["Full TOI ID"],
                    "outer_toi": b["Full TOI ID"],
                    "inner_period": a["Orbital Period (days) Value"],
                    "outer_period": b["Orbital Period (days) Value"],
                    "period_ratio": b["Orbital Period (days) Value"]
                    / a["Orbital Period (days) Value"],
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
                }
            )
    pairs = pd.DataFrame(pairs)
    pairs["abs_log_radius_ratio"] = np.abs(
        np.log(pairs["outer_radius"] / pairs["inner_radius"])
    )
    pairs["radius_ratio"] = pairs["outer_radius"] / pairs["inner_radius"]
    pairs["peas_10pct"] = pairs["abs_log_radius_ratio"] < np.log(1.10)
    return raw, planets, pairs, meta

