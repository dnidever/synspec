# APOGEE-grid specialization and full-band checks — November 13, 2026

**Accepted:** evaluate the standard molecular broadening multiply directly in MOLOP, while preserving the EXOMOL/general path. **Rejected:** the fully split/vectorized core-and-wing loop, which was slower. No large cache is enabled and no production compiler flags are changed.

## Target supplied by the user

2800 < Teff < 5500 K; −0.5 < log g < 3.0; −1 < [Fe/H] < +0.5; wavelength range 15000–17000 Å. Boundary models are intentionally included as numerical stress tests. The target is recorded in apogee_grid.json. Production wavelength sampling was not specified, so the existing 0.04 Å internal sampling control was retained for consistent comparisons; this is not a statement of the production grid resolution.

All timings below cover the complete 2000 Å interval, not the earlier 100 Å or 1000 Å windows. They are SYNSPEC process wall times, including its initialization/read/selection/transfer/output work, excluding atmosphere preprocessing, compilation, external convolution, case preparation, and grid-file writing.

## Small code change

For IVDWLI=0, the original GVDW routine returns GWM(line,list)*VDWC(depth). MOLOP now evaluates that same expression directly instead of calling the more general routine for every contributing line. The DR17 seven-parameter molecular/water lists take this standard path. For IVDWLI != 0, the original GVDW call, including the EXOMOL H2/He temperature law, remains unchanged.

The coefficient is read from the existing PRFQUA common block; no extra per-line storage or population cache was introduced. Populations, stimulated emission, Doppler widths, selection thresholds and line-profile formulas are unchanged. The large experimental factor cache remains opt-in and was disabled for all timing comparisons here.

## Representative full-band measurements

| Model | Compiler | Previous code | New code | Less time | Measured pairs |
|---|---|---:|---:|---:|---:|
| cool2800 | O2 | 22.466 s | 19.671 s | 12.4% | 1 |
| cool2800 | native | 21.148 s | 18.853 s | 10.8% | 1 |
| giant3500 | O2 | 13.258 s | 12.033 s | 9.2% | 3 |
| giant3500 | native | 12.525 s | 11.316 s | 9.7% | 3 |
| arcturus4286 | O2 | 9.500 s | 9.156 s | 3.6% | 1 |
| arcturus4286 | native | 8.654 s | 8.484 s | 2.0% | 1 |

Native means O3 -march=native; both variants use -fno-automatic and strict math. The 3500 K medians have three paired measured repeats per flag set, following excluded warm-ups. The other two models have single screening pairs, not statistically precise medians. Native code is CPU-specific. Wall and process CPU timings are saved in analysis.json and runs.json.

### Cumulative baseline, on the same full-band 3500 K workload

- Corrected original algorithm, no optimization: **34.893 s**, one measured control run.
- Current improved algorithm with O2: **12.033 s**, three-pair median.
- Current improved algorithm with native tuning: **11.316 s**, three-pair median.

The corresponding cumulative ratios are about **2.90×** and **3.08×**. These combine compiler and algorithm improvements, not just the latest multiply specialization. The unoptimized denominator is a single run, so the ratios provide context rather than a robust full-grid speedup estimate. Its printed spectrum also matched. This explicitly avoids comparing different wavelength spans or treating O2 as an unoptimized baseline.

## Actual MARCS grid checks

The user supplied a local MARCS model tree. Eight boundary/near-boundary models were selected from unfilled .mod files at [Fe/H]=−1 and +0.5. Five exact named corners were available. For three unavailable low-gravity corners, the nearest available gravity at the same temperature and metallicity was used, and the actual header parameters were recorded. No model structure was synthesized or silently relabeled.

| Actual Teff | Actual log g | Header [Fe/H] | Requested low/high log g | O2 previous | O2 new | Pairs |
|---:|---:|---:|---:|---:|---:|---:|
| 2800 | -0.50 | -1.00 | -0.5 | 20.202 s | 17.853 s | 1 |
| 2800 | 3.00 | -1.00 | +3.0 | 16.691 s | 14.750 s | 1 |
| 5500 | 0.00 | -1.00 | -0.5 | 6.284 s | 6.302 s | 3 |
| 5500 | 3.00 | -1.00 | +3.0 | 6.370 s | 6.441 s | 3 |
| 2800 | 0.00 | +0.50 | -0.5 | 29.822 s | 23.222 s | 1 |
| 2800 | 3.00 | +0.50 | +3.0 | 24.510 s | 21.580 s | 1 |
| 5500 | 1.00 | +0.50 | -0.5 | 7.233 s | 7.126 s | 1 |
| 5500 | 3.00 | +0.50 | +3.0 | 7.322 s | 7.262 s | 1 |

All **12 corner/near-corner comparison pairs** produced finite, exactly equal printed flux and continuum arrays. At 2800 K the single-pair time savings were roughly 12–22%; those values are screening results, not grid-wide guarantees. At 5500 K the effect is approximately neutral. Two initially slow 5500 K / [Fe/H]=−1 comparisons were repeated; the large apparent regressions did not repeat, and the resulting three-pair medians differ by about 0–1%. All raw observations, including the slow first runs, are retained.

Several .filled files in the model tree have header parameters different from their filenames (including substitutions in gravity, composition, or metallicity). They were not used as if they were independent exact-corner atmospheres. This is a provenance observation, not a modification or scientific judgment about the production filling policy. Missing unfilled hot/very-low-gravity combinations are still not validated as exact corners.

## Numerical checks

- The accepted-change series has 24 physical executions, including warm-ups; every flux/continuum array is finite and identical to its saved current-O2 reference at printed precision.
- All 12 additional actual-MARCS paired comparisons pass identically.
- Both fixed-structure metallicity-like stress pairs (all metals rescaled to the −1/+0.5 endpoints relative to the 3500 K model) pass. These are opacity/equilibrium stresses, not new self-consistent atmosphere models.
- The general-purpose molecular/water list case and a DR17 case with the optional cache enabled also match their saved reference outputs.
- Bounds/FP-checked molecular regressions pass in both cache modes. They now use the actual immutable GVDW implementation with initialized standard and EXOMOL parameters instead of a broadening stub, and exercise temperature/density coefficient changes. The existing 16,551 hydrogen interval checks still pass.

The model/case builder records model and line-list hashes, actual header identifiers and MARCS header metallicities. Its produced input deck was byte-compared with the measured 3500 K case. Header/source provenance is preferred over filename labels.

## Rejected vectorization experiment

The core/far-wing partition candidate did get all three molecular wavelength ranges vectorized with 32-byte vectors under native tuning. However, its 3500 K full-band native median was 14.997 s versus 12.583 s for the current loop, across three measured pairs. It was also slower in the other screening models. Thus generated SIMD did not translate into a synthesis gain; the partition overhead and generated code costs were not separately isolated. This candidate is retained in the reproduction archive only and is not in the accepted source.

## Scope and reproduction

GNU Fortran 11.1.0 on this Intel Mac, same runtime/SDK setup as earlier. No fast-math. Large lists are linked externally. The provisional DR17 atomic binary still omits the previously documented malformed damping row; this data issue is unchanged. These comparisons do not certify NLTE, spherical radiative transfer, every atmosphere-grid combination, all compilers, or internal bitwise equivalence beyond printed output precision. Global bounds-checking limitations documented earlier remain.

Prepare a full-band external input deck with:

```sh
python3 benchmarks/prepare_apogee_case.py --model /path/to/model.mod --linelist-dir /path/to/converted-dr17-lists --output /path/to/new-case
```

The builder preserves atmosphere abundances, enforces target temperature/gravity bounds and MARCS header metallicity bounds, and does not interpolate missing models. The production sampling control can be supplied explicitly with --sampling. The existing comparison runners accept externally prepared decks/manifests. The reproduction archive contains the actual scripts and raw timing/validation metadata used for this experiment, plus source snapshots; model and line-list data remain external.

The optimization branch is the review destination. main, the original project checkout, production compiler defaults, and the user’s MARCS/line-list files are not changed.
