# TESS Radius-Valley Size-Similarity Study

**Completed computational observational study — August 8, 2026**  
**Status:** Original analysis; not peer reviewed

> **Headline finding:** In a quality-controlled public TESS multi-planet sample, a hard 1σ interval-overlap rule expands the number of adjacent pairs labeled as radius-valley-inclusive from **10 to 24**, even though direct uncertainty propagation implies only **10.33 expected latent valley-inclusive pairs**. That hard classification produces a borderline-significant apparent loss of neighboring-planet size similarity, but the effect is not supported when measurement uncertainty is propagated probabilistically and is even stronger in a nonadjacent placebo.

## Study

**Hard uncertainty-overlap labeling can overstate evidence for a radius-valley break in TESS multi-planet size similarity**

The primary sample contains **74 period-adjacent pairs in 58 systems** (132 unique planets/candidates).

| Analysis | Valley-inclusive result | Null-adjusted contrast | Evidence |
|---|---:|---:|---|
| Nominal radii | 1/10 size-matched | `C = 0.056` | permutation `p = 0.419`; cluster-bootstrap 95% interval `[-0.217, 0.249]` |
| Hard 1σ overlap | 2/24 size-matched | `C = 0.170` | Fisher `p = 0.0327`; permutation `p = 0.0458`; cluster interval `[-0.022, 0.349]` |
| Probabilistic error propagation | 10.33 expected latent valley pairs | median `C = 0.037` | 95% simulation interval `[-0.231, 0.229]` |
| Nonadjacent overlap placebo | 0/13 size-matched | `C = 0.362` | Fisher `p = 0.0054`; permutation `p ≈ 2×10⁻⁵` |

The astrophysical conclusion is deliberately conservative: **current TESS data do not independently confirm or refute** the previously reported Kepler radius-valley neighbor anomaly. The stronger contribution is methodological—hard uncertainty-overlap labeling can make cross-survey evidence look stronger than the underlying membership probabilities justify.

The overlap rule obeys

```text
[r - σ, r + σ] intersects [c - w, c + w]
    iff
|r - c| <= w + σ
```

so larger measurement uncertainties make an object more likely to receive the binary “in-valley” label. Here the overlap rule labels 24 adjacent pairs, while the uncertainty model implies only 10.33 expected latent valley-inclusive pairs. Direct uncertainty propagation does not support the hard-label result.

## Read the complete paper

- [`report.md`](report.md) — paper index and concise result
- [`paper/`](paper/) — seven ordered sections containing the complete scientific report
- [`summary.json`](summary.json) — machine-readable principal statistics
- [`search_log.md`](search_log.md) — candidate screen, literature search, and result-specific novelty review
- [`figures/`](figures/) — five vector figures

## Reproducibility

The complete analysis logic is published. [`tess_valley_analysis.py`](tess_valley_analysis.py) is a launcher that executes the original analysis source stored in six ordered chunks under [`analysis_parts/`](analysis_parts/).

The original monolithic analysis script has SHA-256:

```text
3e91c9eff0d0099f3607391f36155c8248c5223ecdb23056329b62cf4b04d6a3
```

GitHub's text-file transport omitted the terminal newline from chunks 01–05 while leaving all substantive source text unchanged. To reconstruct the original byte-for-byte file, append one newline to each of `part_01.py.inc` through `part_05.py.inc`, then concatenate parts 01–06 in order. The result has the SHA-256 above. The launcher itself is unaffected because it compiles each chunk separately; this execution path was tested and reproduced the 74-pair sample, nominal counts, overlap counts, and expected latent pair count.

Additional reproducibility files:

- [`independent_validation.py`](independent_validation.py) — separately written primary-result validator
- [`independent_validation.json`](independent_validation.json) — validation output; all audited checks pass
- [`tables/primary_adjacent_pairs.csv`](tables/primary_adjacent_pairs.csv) — primary 74-pair dataset used by the validator
- [`requirements.txt`](requirements.txt) — pinned Python dependencies
- [`live_endpoint_check.json`](live_endpoint_check.json) — catalog provenance and checksum record
- [`download_data.py`](download_data.py) — downloads the public MAST source catalog and checks whether it matches the study snapshot

### Exact source-data snapshot

The source is the public cumulative **TESS Objects of Interest** catalog from MAST. The exact study snapshot has SHA-256:

```text
6797d7ea49dedd95cf5a5711fe4c93553f218268d3ce68f9fbde569d8ff06a7e
```

On August 8, 2026, the live MAST endpoint served the same byte-for-byte catalog; its embedded header identifies the file as created March 14, 2026 and covering Sectors 1–97.

Download/check it with:

```bash
python download_data.py
```

If MAST has since changed the catalog, the downloader warns that the run is an **updated replication** rather than an exact reproduction.

### Run the study

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

The main pipeline regenerates the cleaned tables, hard classifications, robustness grid, center scan, cluster-bootstrap/permutation samples, measurement-error simulations, statistical summary, and figures. Large Monte Carlo draw files are regenerated rather than committed.

## Validation and falsification

The completed study includes graph-preserving permutations, system-cluster bootstrapping, leave-one-system-out analysis, Gaussian and lognormal measurement-error propagation, radius-valley-center uncertainty propagation, quality-threshold grids, confirmed-planets-only sensitivity checks, alternative valley widths, nonadjacent-pair placebos, and an independent implementation of the primary counts and latent-membership calculation.

## Novelty and limitations

A result-specific literature search through August 8, 2026 did not locate a prior empirical TESS report of this exact adjacent-pair analysis, the **2.4× overlap-label inflation**, its comparison with probabilistic latent membership, or the nonadjacent falsification result. The novelty claim is **provisional pending independent review**.

The nominal TESS valley sample is small and underpowered, the catalog includes planet candidates as well as confirmed planets, and the underlying stellar/planetary measurements are heterogeneous. This repository **does not claim the original Kepler result is false**.

## Data source

MAST cumulative TOI catalog:

`https://archive.stsci.edu/missions/tess/catalogs/toi/tois.csv`

## Citation

See [`CITATION.cff`](CITATION.cff).
