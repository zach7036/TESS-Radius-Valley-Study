# Hard uncertainty-overlap labeling can overstate evidence for a radius-valley break in TESS multi-planet size similarity

**An independent cross-survey reanalysis with measurement-error propagation and adjacency placebos**

**Study date:** 8 August 2026  
**Study type:** Completed computational observational study  
**Status:** Original analysis; not peer reviewed  
**Data:** Public TESS Objects of Interest catalog from MAST  
**Code and complete outputs:** Included with this report

## Abstract

Planets in compact multi-planet systems are often similar in size, a pattern called "peas in a pod." A 2026 Kepler analysis reported that this physical size-similarity enhancement disappears specifically for adjacent pairs containing a planet in the exoplanet radius valley, with a difference-in-differences test of p = 0.001. Whether that result appears independently in TESS systems had not been reported. This study tested it using the current public MAST TESS Objects of Interest catalog available on 8 August 2026. After excluding nonplanet dispositions, restricting hosts to main-sequence G/K/M dwarfs, requiring both radii to lie between 0.5 and 4 Earth radii, and requiring fractional radius uncertainty no greater than 20%, the primary sample contained 74 period-adjacent pairs in 58 systems.

With nominal radii, 1 of 10 radius-valley-inclusive pairs (10.0%, 95% Wilson interval 1.8%-40.4%) and 16 of 64 non-valley pairs (25.0%, 16.0%-36.8%) were size-matched to within 10%. The type-stratified, null-adjusted similarity contrast was 0.056; a within-spectral-type network permutation test gave p = 0.419 and a system-cluster bootstrap gave a 95% interval of -0.217 to 0.249. Thus the nominal TESS data did not independently confirm the Kepler effect, but the test was underpowered: even an extreme true valley size-match rate of zero would have produced only 32.6% power in a simple Fisher-test diagnostic.

A hard rule copied from the Kepler study - labeling a planet as a valley member whenever its reported 1-sigma radius interval overlaps the valley window - changed the result. It produced 24 labeled pairs and a contrast of 0.170 (Fisher p = 0.0327; permutation p = 0.0458). However, the rule is algebraically equivalent to widening the valley window from w to w + sigma for each planet. It therefore admitted 14 additional pairs while a Gaussian measurement model implied only 10.33 latent valley-inclusive pairs in the entire sample. The hard overlap labels had an expected positive predictive value of 36.2%, their cluster-bootstrap interval crossed zero (-0.022 to 0.349), and the signal disappeared under stricter uncertainty, signal-to-noise, and confirmed-only cuts. Direct propagation of radius errors produced a median contrast of 0.037 with a 95% interval of -0.231 to 0.229. Most importantly, the overlap rule produced an even larger deficit among nonadjacent pairs, contrary to the claimed localization to immediate neighbors.

The completed result is therefore a qualified negative replication and a methodological finding: current TESS TOI data do not decisively establish an adjacency-specific radius-valley break in size similarity, and hard 1-sigma interval-overlap labeling can create stronger-looking evidence than probabilistic error propagation supports. This does not refute the Kepler result; it identifies a cross-survey classification failure mode that matters when radius uncertainties are larger.

## 1. Introduction

Two population-level regularities organize the observed architecture of small exoplanets. First, the radius distribution is bimodal, with a relative scarcity of planets between the super-Earth and sub-Neptune populations - the radius valley. Second, planets within the same compact multi-planet system tend to resemble one another in radius and to have comparatively regular spacing, the "peas-in-a-pod" pattern [1-3]. Both patterns constrain theories of planet formation, atmospheric loss, migration, and late dynamical evolution.

Chance and Ballard recently connected the two phenomena. Using 1,719 Kepler planets orbiting 690 stars, they identified valley membership with stellar-mass-dependent centers and a narrow radius window, then compared observed size similarity with the similarity expected from the underlying radius distribution [4]. Their revised analysis reported that non-valley adjacent pairs retain a physical excess of size matches, whereas valley-inclusive pairs do not; the equality of the two enhancements was rejected at p = 0.001. They further reported that the effect weakened when classification criteria were relaxed and that the ordinary near-unity peak returned for nonadjacent pairs, suggesting localization to the immediate neighbors of valley planets [4].

That is a scientifically consequential claim. If robust across surveys, it would imply that planets occupying the valley have a systematically different architectural context, potentially reflecting stochastic evolution such as giant impacts. But an external replication is nontrivial. TESS and Kepler observe different stellar populations, TESS multi-candidate samples are smaller, and many public TESS radius estimates are less precise. These differences make uncertainty handling central rather than incidental.

This study therefore asked a strategically narrow question that could be answered completely with public data and local computation: **Does the reported loss of adjacent-planet size similarity around radius-valley planets appear independently in the current TESS TOI sample, and how does the conclusion depend on treating radius uncertainty as a hard label versus a probability?**

The main result is not a clean positive replication. Nominal and measurement-error-propagated TESS analyses are inconclusive. A borderline positive result appears only under a hard 1-sigma interval-overlap definition. That rule expands the accepted radius region in direct proportion to measurement uncertainty, inflates the hard-labeled pair count 2.4-fold relative to nominal classification, and produces a stronger effect among nonadjacent pairs. The defensible conclusion is therefore methodological and negative: the present TESS evidence is not adjacency-specific or robust enough to confirm the Kepler result.

## 2. Literature and novelty review

### 2.1 What was already known

Weiss et al. established that planets within Kepler multi-planet systems tend to be similar in size and regularly spaced [1]. Subsequent work showed that the pattern is astrophysical rather than solely a detection artifact and that it can be studied in additional planet properties [2,3]. Otegi et al. examined radii, masses, densities, and period ratios in a sample enriched by TESS follow-up and found broad intra-system similarity, but did not test whether planets *inside* the radius valley lose similarity specifically with adjacent neighbors [3].

Goyal and Wang compared systems lying entirely below the valley with systems lying entirely above it [5]. That addresses whether rocky and volatile-rich systems differ in overall uniformity, not whether a pair containing a planet in the low-occurrence valley itself behaves anomalously.

Chance and Ballard directly tested the in-valley question with Kepler [4]. Their work therefore establishes the closest precedent and the target for independent evaluation. Their default rule required a planet's 1-sigma radius interval to overlap a +/-0.10 Earth-radius window around a fitted valley center, together with fractional radius uncertainty no greater than 5%. The revised paper reported 39 mixed or valley-inclusive adjacent pairs and 769 non-valley pairs, a null-adjusted difference-in-differences p value of 0.001, dilution under looser membership criteria, and recovery of the near-unity peak in nonadjacent pairs [4].

Separately, Parashivamurthy and Mulders measured the TESS radius valley for low-mass stars and reported three-component mixture minima of 1.64 +/- 0.03, 1.75 +/- 0.11, and 1.86 +/- 0.06 Earth radii for M, K, and G dwarfs, respectively [6]. Those survey-specific values make a TESS replication possible without importing a Kepler-only valley location.

### 2.2 Candidate-question screen

The study question was selected only after rejecting alternatives that were already answered or not resolvable with the available sample.

| Candidate | Status after literature search | Decision |
|---|---|---|
| Do TESS multi-planet systems exhibit general size similarity? | Substantially addressed by Otegi et al. and related work. | Rejected as non-novel. |
| Is size uniformity different in systems wholly below versus wholly above the radius valley? | Addressed by Goyal and Wang. | Rejected as non-novel. |
| Do Kepler pairs containing an in-valley planet lack the ordinary similarity enhancement? | Directly addressed by Chance and Ballard. | Rejected as a new discovery target; retained as the replication target. |
| Is the effect independently present within TESS M-dwarf pairs alone? | Scientifically attractive, but the primary quality-controlled sample contains only 11 M-dwarf adjacent pairs. | Rejected as too underpowered for a definitive standalone study. |
| Does the Kepler radius-valley neighbor result replicate in TESS, and does hard interval-overlap classification remain valid with noisier radii? | No prior empirical TESS replication or quantification of this hard-label inflation was found. Public data and complete computation were available. | Selected. |

### 2.3 Pre-result and result-specific novelty searches

Before analysis, searches combined synonyms for TESS, TOIs, radius valley/radius gap, adjacent planets/neighbors, size similarity/uniformity, and peas in a pod. After the result emerged, additional searches used the specific concepts "1-sigma overlap," "uncertainty-overlap classification," "radius-valley membership uncertainty," "TESS radius-valley neighbor similarity," and the Chance and Ballard arXiv identifier. The searches were repeated on 8 August 2026 across arXiv, journal pages, and web-indexed scholarly results.

The closest results remained: (i) general TESS multi-planet similarity [3], (ii) comparisons of systems on opposite sides of the valley [5], (iii) the Kepler in-valley neighbor analysis [4], and (iv) TESS measurement of the valley location itself [6]. No located work reported an empirical TESS test of the specific adjacent-pair anomaly, quantified the 2.4-fold pair-label inflation caused by the hard overlap rule in a TESS sample, compared it with latent-membership probabilities, or tested adjacency localization with TESS nonadjacent pairs. The novelty claim must remain provisional because no literature search can prove nonexistence, but the result appears previously unreported as of the search date.

## 3. Research question and hypotheses

### 3.1 Primary question

Among period-adjacent TESS planet candidates and confirmed planets orbiting main-sequence G, K, and M dwarfs, is the physical excess of size-matched pairs lower when at least one member lies in the spectral-type-specific radius valley?

### 3.2 Methodological question

Does a hard label based on overlap between a reported 1-sigma radius interval and the valley window yield evidence consistent with direct probabilistic propagation of the same radius uncertainties?

### 3.3 Falsifiable hypotheses

**H1, astrophysical replication hypothesis.** Valley-inclusive adjacent pairs have a smaller null-adjusted size-similarity enhancement than non-valley adjacent pairs. The planned directional statistic is positive when this is true.

**H2, adjacency-localization hypothesis.** If the effect reflects the immediate architectural neighborhood of a valley planet, it should weaken or disappear for nonadjacent pairs from the same quality-controlled systems.

**H3, uncertainty-consistency hypothesis.** If a hard 1-sigma overlap result reflects latent valley membership rather than classification geometry, its effect estimate should remain compatible with Monte Carlo propagation of reported radius uncertainties and should not strengthen merely as lower-precision objects are admitted.
