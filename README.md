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

The apparent hard-overlap replication is fragile because the rule changes the effective radius window for each object:

```text
[r - σ, r + σ] intersects [c - w, c + w]
    iff
|r - c| <= w + σ
```

Larger measurement uncertainties therefore make a planet *more likely* to receive the binary “in-valley” label. In this sample the overlap rule labels **24 adjacent pairs**, while the reported uncertainty model implies only **10.33 expected latent valley-inclusive pairs**. The estimated positive predictive value of the overlap labels is **36.2%**.

When radius uncertainty is propagated directly rather than collapsed into binary overlap labels, the evidence for a similarity break disappears.

## Scientific contribution

The study contributes a **qualified negative replication plus a methodological result**:

> In the current public TESS multi-planet sample, the adjacent radius-valley size-similarity anomaly is not independently established under nominal or probabilistically propagated radii. A borderline signal appears under hard 1σ interval-overlap labeling, but that rule substantially inflates valley membership and the resulting signal fails measurement-error and adjacency-specific falsification tests.

A result-specific literature search through August 8, 2026 did not identify a prior empirical report of this exact TESS adjacent-pair analysis, the **2.4× overlap-label inflation**, its comparison with probabilistic latent membership, or the nonadjacent falsification result. The novelty claim is therefore **provisional pending independent review**.

## Read the study

- [`report.md`](report.md) — complete scientific paper in Markdown
- [`summary.json`](summary.json) — machine-readable statistical results
- [`search_log.md`](search_log.md) — literature/novelty search record

## Reproduce the analysis

The analysis code is fully included:

- [`tess_valley_analysis.py`](tess_valley_analysis.py) — primary pipeline
- [`independent_validation.py`](independent_validation.py) — independently written validation implementation
- [`requirements.txt`](requirements.txt) — pinned dependencies
- [`live_endpoint_check.json`](live_endpoint_check.json) — provenance and snapshot checksum

The source data are the public cumulative **TESS Objects of Interest** catalog from MAST. Download it to `tois_mast.csv`:

```bash
python download_data.py
```

The study used a byte-for-byte catalog with SHA-256:

```text
6797d7ea49dedd95cf5a5711fe4c93553f218268d3ce68f9fbde569d8ff06a7e
```

The MAST file served on August 8, 2026 matched this hash and identified itself as created March 14, 2026, covering TESS Sectors 1–97. If the live catalog changes, `download_data.py` will warn that the result is an **updated replication**, not an exact rerun of the frozen snapshot.

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

## Included derived data

The `tables/` directory contains the core pair-level samples, classification comparisons, robustness grid, center placebo scan, leave-one-system-out results, and power curves. Large Monte Carlo/permutation draw files are intentionally not required for version control because the included scripts deterministically regenerate them from the stated seed and inputs.

## Validation and falsification tests

The study includes:

- graph-preserving permutations
- system-cluster bootstrapping
- leave-one-system-out analysis
- Gaussian measurement-error propagation
- lognormal error-model sensitivity testing
- radius-valley-center uncertainty propagation
- quality-threshold grids
- confirmed-planets-only sensitivity checks
- alternative valley widths
- nonadjacent-pair placebos
- independent implementation of the core result

## Limitations

The nominal TESS valley sample is small and underpowered, the catalog includes planet candidates as well as confirmed planets, and published stellar/planetary measurements are heterogeneous. This study therefore **does not claim that the original Kepler result is false**. It establishes that the present TESS sample is inconclusive under uncertainty-aware classification and documents a cross-survey classification failure mode worth avoiding in future demographic studies.

## Data source

Mikulski Archive for Space Telescopes (MAST), cumulative TESS Objects of Interest catalog:

`https://archive.stsci.edu/missions/tess/catalogs/toi/tois.csv`

## Citation

Citation metadata are provided in [`CITATION.cff`](CITATION.cff).