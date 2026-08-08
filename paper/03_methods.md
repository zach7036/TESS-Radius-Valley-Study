## 5. Methods

### 5.1 Pair construction

Within each TIC system, eligible planets were ordered by reported orbital period. Consecutive planets formed adjacent edges. Importantly, adjacency was established before the primary radius-precision cut, so a low-quality intervening TOI was not silently skipped to create a false adjacent edge. Each primary edge therefore joins consecutive members of the eligible candidate list and has period ratio greater than one.

A pair was called size-matched, or a "pea," when

`abs[ln(R_outer / R_inner)] < ln(1.10)`,

that is, when the radii differed by less than 10% multiplicatively. Sensitivity checks used 5%, 15%, and 20% thresholds.

### 5.2 Radius-valley centers

The adopted fixed centers were the three-component TESS mixture minima reported by Parashivamurthy and Mulders [6]:

| Host type | Center (Earth radii) | Published 1-sigma uncertainty |
|---|---:|---:|
| M | 1.64 | 0.03 |
| K | 1.75 | 0.11 |
| G | 1.86 | 0.06 |

The primary half-width was `w = 0.10` Earth radii, matching the narrow window in the Kepler analysis [4]. Center uncertainty was propagated in a dedicated simulation, and all centers were shifted together from -0.60 to +0.60 Earth radii as a specificity placebo.

### 5.3 Three membership definitions

Let `r` be the reported radius, `sigma` its reported 1-sigma uncertainty, `c` the spectral-type center, and `w = 0.10`.

1. **Nominal:** `abs(r - c) <= w`.
2. **1-sigma overlap:** `[r - sigma, r + sigma]` intersects `[c - w, c + w]`.
3. **Certain:** `[r - sigma, r + sigma]` lies wholly inside `[c - w, c + w]`.

The overlap rule has a crucial algebraic property:

`[r - sigma, r + sigma] intersects [c - w, c + w]` if and only if `abs(r - c) <= w + sigma`.

It is therefore not a fixed valley definition. Every object's accepted distance from the valley grows linearly with its uncertainty. A noisier radius can be farther from the valley center and still receive the same binary label as a precise radius near the center.

A pair was valley-inclusive when either endpoint satisfied the selected membership rule.

![Figure 1. The hard overlap rule expands the accepted distance from the valley center from a fixed 0.10 Earth radii to 0.10 plus each object's reported uncertainty. Open circles are overlap-only labels in the eligible multi-candidate pool before the primary uncertainty cut.](../figures/figure1_overlap_geometry.svg)

### 5.4 Type-stratified random-pair baseline

Raw size-match rates depend on where radii lie in the population distribution. To separate that intrinsic contribution from within-system organization, the analysis followed the logic of the Kepler difference-in-differences test [4]. For each spectral type:

- The non-valley baseline was the probability that two distinct randomly drawn non-valley planets were within 10% in radius.
- The valley-inclusive baseline was the probability that a randomly drawn valley planet and a distinct randomly drawn planet from the full type-specific pool were within 10%.

For each pair class, the **physical enhancement** was observed size-match rate minus its type-weighted random-pair baseline. The primary contrast was

`C = (observed_non - null_non) - (observed_valley - null_valley)`.

Positive `C` means the ordinary physical similarity enhancement is weaker for valley-inclusive pairs. This is the directional replication statistic.

### 5.5 Classical and correlation-aware tests

The analysis reported:

- Wilson 95% intervals for observed binary size-match rates.
- One-sided Fisher exact tests of the raw rate deficit.
- One-sided Mann-Whitney tests on absolute log radius ratios.
- 50,000 within-spectral-type permutations that shuffled paired radius/error records among planet nodes while preserving the system graph, multiplicities, spectral-type radius distribution, and uncertainty distribution.
- 50,000 system-cluster bootstrap draws, resampling TIC systems rather than treating all pair edges as independent.
- Leave-one-system-out estimates.

The Monte Carlo standard error of the primary overlap permutation p value, 0.0458 from 50,000 draws, is approximately 0.00094; a rough 95% Monte Carlo error band is therefore +/-0.0018.

### 5.6 Probabilistic measurement-error analysis

Hard labels discard how likely a radius is to lie inside the valley. For each planet, the Gaussian latent-membership probability was

`P(c - w <= R_true <= c + w)`

under `R_true ~ Normal(r, sigma)`. Pair membership probability was `1 - (1 - p_inner)(1 - p_outer)`. Summing those probabilities gave the expected latent number of valley-inclusive pairs.

The full effect was also recomputed over 10,000 draws under each of three models:

1. positive-truncated normal radii with the reported standard deviation;
2. lognormal radii matched to the reported mean-scale fractional uncertainty;
3. normal radii while also sampling the published G/K/M center uncertainties.

Each draw reclassified planets, recomputed which edges were size-matched, and recalculated the type-specific random-pair baseline. A separate 20,000-draw simulation varied only the published center uncertainties.

### 5.7 Robustness and falsification tests

The potential result was challenged with:

- radius-error cutoffs of 10%, 15%, 20%, 25%, and no finite cutoff;
- valley half-widths of 0.08, 0.10, 0.12, 0.15, and 0.20 Earth radii;
- size-match thresholds of 5%, 10%, 15%, and 20%;
- maximum planet radii of 3, 4, and 6 Earth radii;
- signal-to-noise thresholds of 7.1, 10, and 15;
- all accepted dispositions, confirmed/known planets only, and PC-only pairs;
- additive, multiplicative, and logit models for transferring the non-valley enhancement to valley pairs;
- common center shifts across +/-0.60 Earth radii;
- nonadjacent-pair placebo tests;
- an independent implementation that reconstructed hard labels and Fisher tests without importing the main analysis code.
