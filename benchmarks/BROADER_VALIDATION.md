# Broader validation and improved-build profile — November 13, 2026

**All 24 baseline/improved comparison pairs passed with identical finite printed flux and continuum arrays.** This broadens the numerical regression coverage to four atmospheres and two wavelength windows, but does not establish physical correctness for all stellar parameters.

## Atmospheres and wavelength coverage

Existing local atmosphere files were used: an Arcturus model (4286 K, log g=1.66), a solar Kurucz dwarf (5777 K, log g=4.437), a cool MARCS giant (2800 K, log g=0), and the repository MARCS giant test model (6000 K, log g=1). The two windows were H-band 16000–16100 Å and H-beta 4840–4880 Å. Microturbulence and abundances came from each model; the synthesis setup remained LTE/plane-parallel.

The general-purpose atomic, solar-isotope molecular, TiO and H2O lists were used for both builds because these files cover optical and infrared wavelengths. This is a comparative code regression, not a claim that solar isotope ratios are scientifically appropriate for each giant or that the line data are calibrated to every model. The DR17-specific lists do not cover the optical window and were not substituted there.

Each case had one excluded warm-up and two measured runs per build, in shuffled sequential order. Baseline is the repaired O2 build before the two profile-loop changes; improved is the low-memory O2 version with the optional large cache disabled. These short-window medians are illustrative, not statistically precise speedup estimates.

| Case | Baseline median | Improved median | Less time |
|---|---:|---:|---:|
| arcturus-Hband | 4.121 s | 3.984 s | +3.3% |
| arcturus-Hbeta | 3.731 s | 3.578 s | +4.1% |
| solar_dwarf-Hband | 4.550 s | 4.346 s | +4.5% |
| solar_dwarf-Hbeta | 3.549 s | 3.529 s | +0.6% |
| cool_giant-Hband | 4.815 s | 4.681 s | +2.8% |
| cool_giant-Hbeta | 7.841 s | 7.530 s | +4.0% |
| warm_giant-Hband | 3.827 s | 3.659 s | +4.4% |
| warm_giant-Hbeta | 2.861 s | 2.939 s | -2.7% |

Most of these setup-heavy short cases show only about 3–4% lower wall time. The 6000 K H-beta case was about 3% slower in two measured repeats; the change is not universally faster. This does not contradict the earlier five-repeat, 1,000 Å DR17 result of about 13% less time, where profile work occupies more of the run. System load, CPU frequency and cache state were not controlled.

## Re-profile of the improved O2 executable

A separate -O2 -g -fno-omit-frame-pointer build of the accepted source was sampled on the original 3500 K / 15500–16500 Å DR17 reference. The large factor cache was disabled, and explicit stage timers were not present. macOS sample requested 1 ms observations for up to 30 s and ended when synthesis exited. It collected 6,978 main-thread stack observations. The sampled flux matched the saved baseline exactly at printed precision: True.

| Largest exclusive sampled leaves | Observations |
|---|---:|
| MOLOP | 2,292 |
| exp | 1,498 |
| HYDLIN | 675 |
| OPACON | 432 |
| RUSSEL | 322 |
| GVDW | 194 |
| LINEQS | 139 |

These are sampled observations, not exact per-routine seconds, and counts from a different run cannot be subtracted as a controlled speedup measurement. The call graph still identifies molecular line-strength exponentials and profile accumulation as major costs. The hydrogen search is much less prominent after the binary lookup. VOIGTK itself has only 24 exclusive leaves now because the common branch is evaluated directly in MOLOP; that work has not disappeared, it is attributed to MOLOP. Raw symbols and source locations are in the supplied sample report.

## Hydrogen interpolation-width experiment

A separate, uncommitted candidate precomputed reciprocal table-interval widths per hydrogen line/depth, with an original-division fallback for a zero-width interval. On three measured wide-reference pairs it changed the median from 8.008 s to 7.862 s (1.8% less time). One comparison for each of the eight broader cases showed mixed timing changes, including a slower solar-dwarf H-band run. All 12 candidate/current comparison pairs, including the warm-up, retained identical finite printed outputs.

This small result is not strong enough to add another interpolation path to the accepted source. The candidate remains in the reproduction archive only. No reciprocal-width calculation is enabled in the code being prepared for review.

## Accepted changes and checks

The optimization branch contains the three array-bound corrections, the molecular small-damping Voigt path, the hydrogen binary lookup, and the optional factor cache disabled by default. The focused molecular regression now uses immutable pre-optimization MOLOP/VOIGTK routines stored in benchmarks/reference_molop_voigt.f, rather than a moving git HEAD. Both cache modes passed the bounds/FP checked regression, and 16,551 hydrogen interval tests passed.

The main branch and production compiler defaults are unchanged. No large line-list data are committed. Broader atmosphere-grid, NLTE, spherical transfer and alternative compiler validation remain outside this experiment; the previously documented general bounds-checking and DR17 atomic-data limitations are not declared resolved.

## Next optimization direction

Prioritize the remaining molecular opacity/exponential work. Re-testing vectorization/native compiler options on the new inlined wavelength loop is a useful next experiment: removing scalar calls may change what the compiler can vectorize. Any resulting gain must be measured again with spectrum checks; no additional speedup from that idea is asserted here.
