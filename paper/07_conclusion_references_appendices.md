## 14. Conclusion

Previous research established that Kepler planets in compact systems generally exhibit physical size similarity and reported that this enhancement disappears for adjacent pairs containing a radius-valley planet. Whether the same effect appears independently in TESS, and whether the hard interval-overlap membership rule remains trustworthy with noisier radii, were unresolved.

This study analyzed the current public MAST TOI catalog and obtained 74 quality-controlled adjacent pairs in 58 G/K/M-dwarf systems. Nominal radii produced a small positive but nonsignificant null-adjusted contrast (`C = 0.056`, permutation `p = 0.419`, cluster-bootstrap 95% interval `-0.217 to 0.249`). The test is underpowered and therefore does not falsify the Kepler result.

A hard 1-sigma overlap rule produced an apparently stronger contrast (`C = 0.170`, permutation `p = 0.0458`), but it expanded the valley window by each object's uncertainty, labeled 24 pairs where the latent-radius model expected 10.33, failed stricter quality controls, disappeared under direct error propagation, and became stronger among nonadjacent pairs. The adversarial checks therefore do not support interpreting that p value as an independent astrophysical replication.

What has now been learned is specific: **with the current TESS TOI sample, evidence for an adjacency-specific radius-valley break in planetary size similarity is inconclusive, and hard 1-sigma interval-overlap classification can overstate the evidence relative to probabilistic treatment of the same measurements.**

## References

1. Weiss, L. M., et al. (2018). Peas in a Pod: Planets in a Kepler Multi-Planet System Are Similar in Size and Regularly Spaced. *The Astronomical Journal*, 155, 48. DOI: 10.3847/1538-3881/aa9ff6.
2. Weiss, L. M., & Petigura, E. A. (2020). The Kepler Peas in a Pod Pattern Is Astrophysical. *The Astrophysical Journal Letters*, 893, L1. DOI: 10.3847/2041-8213/ab7c69.
3. Otegi, J. F., Helled, R., & Bouchy, F. (2022). The similarity of multi-planet systems. *Astronomy & Astrophysics*, 658, A107. DOI: 10.1051/0004-6361/202142110.
4. Chance, Q., & Ballard, S. (2026). Evidence that Planets in the Radius Gap Do Not Resemble Their Neighbors. *The Astronomical Journal*. DOI: 10.3847/1538-3881/ae77ec; arXiv:2410.02150v2.
5. Goyal, A. V., & Wang, S. (2024). Peas-in-a-Pod Across the Radius Valley: Rocky Systems are Less Uniform in Mass but More Uniform in Size and Spacing. *The Astrophysical Journal Letters*, 968, L4. DOI: 10.3847/2041-8213/ad4f6e.
6. Parashivamurthy, H. M., & Mulders, G. D. (2025). Radius valley scaling among low-mass stars with TESS. *Astronomy & Astrophysics*, 703, A8. DOI: 10.1051/0004-6361/202554006.
7. Mikulski Archive for Space Telescopes. TESS Objects of Interest cumulative catalog, `tois.csv`, accessed 8 August 2026. https://archive.stsci.edu/missions/tess/catalogs/toi/tois.csv.
8. Ricker, G. R., et al. (2015). Transiting Exoplanet Survey Satellite (TESS). *Journal of Astronomical Telescopes, Instruments, and Systems*, 1, 014003. DOI: 10.1117/1.JATIS.1.1.014003.
9. Fulton, B. J., et al. (2017). The California-Kepler Survey. III. A Gap in the Radius Distribution of Small Planets. *The Astronomical Journal*, 154, 109. DOI: 10.3847/1538-3881/aa80eb.

## Appendix A. Primary definitions at a glance

| Element | Definition |
|---|---|
| Adjacent pair | Consecutive eligible TOIs after ordering by orbital period within TIC |
| Primary radius range | 0.5-4.0 Earth radii for both endpoints |
| Primary maximum fractional radius uncertainty | 20% for both endpoints |
| Pea / size match | Absolute log radius ratio below ln(1.10) |
| Valley centers | M 1.64, K 1.75, G 1.86 Earth radii |
| Valley half-width | 0.10 Earth radii |
| Nominal membership | Best-fit radius lies inside fixed window |
| Overlap membership | Reported 1-sigma interval intersects fixed window |
| Certain membership | Reported 1-sigma interval lies wholly inside fixed window |
| Valley-inclusive pair | At least one endpoint is classified as valley |
| Primary contrast | Non-valley physical enhancement minus valley physical enhancement |
| Positive contrast | Less physical size similarity around valley-labeled planets |

## Appendix B. Interpretation guardrails

- A nonsignificant nominal result is not proof of no effect because power is low.
- A nominally significant overlap result is not automatically a discovery because classification is uncertainty-dependent and robustness is weak.
- The nonadjacent placebo is a falsification test, not a claim of a new nonadjacent astrophysical population.
- Expected label precision and expected latent counts depend on the stated Gaussian measurement model; they are not external truth labels.
- The novelty assessment is based on a thorough but necessarily fallible literature search and should be independently checked before journal submission.
