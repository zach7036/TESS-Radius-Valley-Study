# Full scientific report

## Hard uncertainty-overlap labeling can overstate evidence for a radius-valley break in TESS multi-planet size similarity

**An independent cross-survey reanalysis with measurement-error propagation and adjacency placebos**  
**Study date:** 8 August 2026  
**Status:** Original completed computational study; not peer reviewed

The full paper is stored in ordered sections so it remains easy to read and audit on GitHub:

1. [Abstract, introduction, literature/novelty review, and hypotheses](paper/01_abstract_introduction_literature.md)
2. [Data and materials](paper/02_data_materials.md)
3. [Methods](paper/03_methods.md)
4. [Analysis and primary results](paper/04_analysis_results.md)
5. [Robustness/falsification tests, discussion, and limitations](paper/05_robustness_discussion_limitations.md)
6. [Novel contribution, scientific significance, and reproducibility](paper/06_contribution_significance_reproducibility.md)
7. [Conclusion, references, and appendices](paper/07_conclusion_references_appendices.md)

## Core result

Using 74 quality-controlled period-adjacent TESS pairs in 58 systems, nominal radius-valley membership produced a small, statistically inconclusive contrast (`C = 0.056`, permutation `p = 0.419`). A hard 1-sigma interval-overlap rule expanded the number of valley-inclusive pairs from 10 to 24 and produced a borderline contrast (`C = 0.170`, permutation `p = 0.0458`). Direct uncertainty propagation instead implied only 10.33 expected latent valley-inclusive pairs and a median contrast of `C = 0.037` with a 95% simulation interval spanning zero (`-0.231` to `0.229`). The hard-overlap effect was also stronger among nonadjacent pairs, contrary to adjacency localization.

**Conclusion:** current TESS evidence is insufficient to independently confirm or refute the reported Kepler astrophysical effect, while the analysis does identify a quantifiable methodological failure mode: hard interval-overlap classification can overstate evidence when measurement uncertainties are appreciable.

See [`README.md`](README.md) for the repository overview and [`summary.json`](summary.json) for machine-readable results.
