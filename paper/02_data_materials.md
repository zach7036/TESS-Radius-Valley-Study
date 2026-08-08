## 4. Data and materials

### 4.1 Catalog

The analysis used the public cumulative TESS Objects of Interest CSV distributed by the Mikulski Archive for Space Telescopes (MAST) [7]. The file contains 7,913 rows. Its embedded header states:

- File created: 2026-03-14
- Covered sectors: S0001-S0097
- Previous TOI-list release: 2026-03-05

The endpoint was downloaded again on 8 August 2026. The live download and the frozen analysis file had the identical SHA-256 checksum:

`6797d7ea49dedd95cf5a5711fe4c93553f218268d3ce68f9fbde569d8ff06a7e`

Thus the study used the exact file then served by the live MAST endpoint, although the archive file's own creation date was March 2026. Both the frozen input and the endpoint-verification record are included.

### 4.2 Host and planet eligibility

Rows were retained when the TOI disposition was PC, CP, or KP; orbital period and planet radius were positive and finite; and the host satisfied broad main-sequence G/K/M cuts:

- 2,600 K <= effective temperature < 6,040 K
- 4.0 <= log g <= 5.5
- 0.1 <= stellar radius <= 1.5 solar radii

Spectral types followed the temperature boundaries used in the TESS valley study [6]: M below 3,880 K, K from 3,880 to below 5,340 K, and G from 5,340 to below 6,040 K. Systems required at least two eligible TOIs. This produced 393 planet rows in 172 systems and 221 period-adjacent edges before the primary radius and precision restrictions.

### 4.3 Primary pair sample

For both members of an adjacent pair, the primary analysis required:

- 0.5 <= planet radius <= 4.0 Earth radii
- finite fractional radius uncertainty <= 20%

The resulting sample contained:

| Quantity | Value |
|---|---:|
| Period-adjacent pairs | 74 |
| Host systems | 58 |
| Unique planet candidates/planets | 132 |
| G-dwarf pairs | 21 |
| K-dwarf pairs | 42 |
| M-dwarf pairs | 11 |
| Confirmed/known endpoints (CP or KP) | 46 |
| Planet-candidate endpoints (PC) | 86 |

The 20% threshold was chosen before the final comparison as a compromise: stricter thresholds sharply reduce the already small TESS multi-planet sample, whereas no threshold admits broad error bars that make a hard overlap rule nearly synonymous with uncertainty size. Results at 10%, 15%, 25%, and no finite cutoff are reported rather than hidden.
