# TESS Radius-Valley Size-Similarity Study

**Completed computational observational study — August 8, 2026**  
**Status:** Original analysis; not peer reviewed

> **Headline finding:** In a quality-controlled public TESS multi-planet sample, a hard 1σ interval-overlap rule expands the number of adjacent pairs labeled as radius-valley-inclusive from **10 to 24**, even though direct uncertainty propagation implies only **10.33 expected latent valley-inclusive pairs**. That hard classification produces a borderline-significant apparent loss of neighboring-planet size similarity, but the effect is not supported when measurement uncertainty is propagated probabilistically and is even stronger in a nonadjacent placebo.

This repository documents the completed study:

> **Hard uncertainty-overlap labeling can overstate evidence for a radius-valley break in TESS multi-planet size similarity**

The astrophysical conclusion is deliberately conservative: **current TESS data do not independently confirm or refute** the previously reported Kepler radius-valley neighbor anomaly. The stronger contribution is methodological—hard uncertainty-overlap labeling can make cross-survey evidence look stronger than the underlying membership probabilities justify.

## Main result

The primary sample contains **74 period-adjacent pairs in 58 systems** (132 unique planets/candidates).

| Analysis | Valley-inclusive result | Null-adjusted contrast | Evidence |
|---|---:|---:|---|
| Nominal radii | 1/10 size-matched | `C = 0.056` | permutation `p = 0.419`; cluster-bootstrap 95% interval `[-0.217, 0.249]` |
| Hard 1σ overlap | 2/24 size-matched | `C = 0.170` | Fisher `p = 0.0327`; permutation `p = 0.0458`; cluster interval `[-0.022, 0.349]` |
| Probabilistic error propagation | 10.33 expected latent valley pairs | median `C = 0.037` | 95% simulation interval `[-0.231, 0.229]` |
| Nonadjacent overlap placebo | 0/13 size-matched | `C = 0.362` | Fisher `p = 0.0054`; permutation `p ≈ 2×10⁻⁵` |

`C` is the difference between the null-adjusted size-similarity enhancement of non-valley and valley-inclusive pairs. Positive values indicate less excess similarity among valley-inclusive pairs.

## Why the hard-overlap rule matters

The overlap rule obeys:

```text
[r - σ, r + σ] intersects [c - w, c + w]
    iff
|r - c| <= w + σ
```

Larger measurement uncertainties therefore make a planet *more likely* to receive the binary “in-valley” label. In this sample it labels **24 adjacent pairs**, while the reported uncertainty model implies only **10.33 expected latent valley-inclusive pairs**. The estimated expected positive predictive value of those overlap labels is **36.2%**.

When radius uncertainty is propagated directly rather than collapsed into binary overlap labels, the evidence for a similarity break disappears.

## Scientific contribution

The study contributes a **qualified negative replication plus a methodological result**:

> In the current public TESS multi-planet sample, the adjacent radius-valley size-similarity anomaly is not independently established under nominal or probabilistically propagated radii. A borderline signal appears under hard 1σ interval-overlap labeling, but that rule substantially inflates valley membership and the resulting signal fails measurement-error and adjacency-specific falsification tests.

A result-specific literature search through August 8, 2026 did not identify a prior empirical report of this exact TESS adjacent-pair analysis, the **2.4× overlap-label inflation**, its comparison with probabilistic latent membership, or the nonadjacent falsification result. The novelty claim is therefore **provisional pending independent review**.

## Read the study

- [`report.md`](report.md) — complete paper index
- [`paper/`](paper/) — seven ordered sections containing the full scientific report
- [`summary.json`](summary.json) — machine-readable principal statistics
- [`search_log.md`](search_log.md) — literature and novelty-search record
- [`figures/`](figures/) — five vector figures

## Reproduce the analysis

The exact analysis source is included. Because the connected publishing interface has per-file transport constraints, the original script is stored byte-for-byte in six ordered chunks under [`analysis_parts/`](analysis_parts/); [`tess_valley_analysis.py`](tess_valley_analysis.py) executes those chunks in the original order. Recombining the six chunks gives SHA-256:

```text
3e91c9eff0d0099f3607391f36155c8248c5223ecdb23056329b62cf4b04d6a3
```

That is identical to the original completed analysis script used for this study.

Other reproducibility files:

- [`independent_validation.py`](independent_validation.py) — separately written validation implementation
- [`independent_validation.json`](independent_validation.json) — validation output; all audited checks pass
- [`tables/primary_adjacent_pairs.csv`](tables/primary_adjacent_pairs.csv) — the 74-pair primary analysis dataset needed by the independent validator
- [`requirements.txt`](requirements.txt) — pinned dependencies
- [`live_endpoint_check.json`](live_endpoint_check.json) — provenance and checksum verification
- [`download_data.py`](download_data.py) — downloads the public MAST source catalog and checks whether it matches the original study snapshot

The source data are the public cumulative **TESS Objects of Interest** catalog from MAST. Download it to `tois_mast.csv`:

```bash
python download_data.py
```

The exact catalog used in the study has SHA-256:

```text
6797d7ea49dedd95cf5a5711fe4c93553f218268d3ce68f9fbde569d8ff06a7e
```

The MAST file served on August 8, 2026 matched this hash and identified itself as created March 14, 2026, covering TESS Sectors 1–97. If the live catalog changes, `download_data.py` warns that the run is an **updated replication**, not an exact rerun of the frozen snapshot.

Then run:

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -r requirements.txt

python tess_valley_analysis.py \
  --data tois_mast.csv \
  --bootstrap 50000 \
  --permutations 50000 \
  --measurement-mc 10000 \
  --center-mc 20000 \
  --batch-size 500

python independent_validation.py
```

Fixed random seed: `20260808`.

The main pipeline deterministically regenerates the complete cleaned tables, robustness grid, center scan, bootstrap/permutation samples, measurement-error simulations, summary, and publication figures from the public input and fixed seed. Large Monte Carlo draw files are intentionally regenerated rather than committed.

## Validation and falsification tests

The study includes graph-preserving permutations, system-cluster bootstrapping, leave-one-system-out analysis, Gaussian and lognormal measurement-error propagation, radius-valley-center uncertainty propagation, quality-threshold grids, confirmed-planets-only sensitivity checks, alternative valley widths, nonadjacent-pair placebos, and an independent implementation of the primary counts and latent-membership calculation.

## Limitations

The nominal TESS valley sample is small and underpowered, the catalog includes planet candidates as well as confirmed planets, and published stellar/planetary measurements are heterogeneous. This study therefore **does not claim that the original Kepler result is false**. It establishes that the present TESS sample is inconclusive under uncertainty-aware classification and documents a cross-survey classification failure mode worth avoiding in future demographic studies.

## Data source

Mikulski Archive for Space Telescopes (MAST), cumulative TESS Objects of Interest catalog:

`https://archive.stsci.edu/missions/tess/catalogs/toi/tois.csv`

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff).
