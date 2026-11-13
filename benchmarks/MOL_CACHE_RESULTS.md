# Molecular-opacity cache experiment — November 13, 2026

The selected implementation reduces the measured 1,000 Å synthesis median from **9.587 s to 8.430 s**: **12.1% less wall time**, or **1.137× throughput**, on top of the repaired O2 build. Five measured reference pairs were run sequentially in shuffled order, after one excluded warm-up per build. Median paired CPU-time speedup is 1.121×.

## Implementation

MOLOP now caches only EXP(GFM − EXCLM/T), per molecular line and atmospheric depth. Tables are filled lazily when a depth is first used at a temperature. RRMOL, DOPMOL, STIM, GVDW, damping, line cutoffs and profile accumulation remain on the ordinary calculation path, preserving their abundance and wavelength dependence. Cached factors retain double precision; the multiplication order is preserved.

Each molecular-list load increments a generation counter. A new generation, changed line count, or changed depth count releases the old table and attempts fresh allocation. Temperature changes refill the corresponding depth. Allocated factor and temperature arrays share a 256 MiB total cap. If the budget is insufficient or allocation fails, that list uses the original exponential calculation. This bounds additional table storage even with the source maximum of nine million molecular lines. All production line-parameter writes to GFM/EXCLM are in INMOLI, which invokes invalidation.

The source change is in the separate local perf/compiler-benchmarks checkout. MOLCACHE.FOR supplies the declarations, and the makefile dependency includes it. The production optimization flags were not changed. This remains single-process/serial state, consistent with the surrounding saved/common-block implementation; it is not a new thread-safe or persistent-synthesis interface.

## Measured memory tradeoff

A fresh-process measurement of this same broad case gave peak resident memory **76.6 MiB baseline versus 219.6 MiB cached**, an increase of **143.0 MiB**. These are measured peaks from RUSAGE_CHILDREN in a fresh Python helper for each executable, not the table-allocation limit. The 256 MiB cap is an upper bound; this workload uses less.

## Validation

All ten physical comparison pairs for the selected 256 MiB implementation passed: one warm-up, five reference repeats, and Fe −0.2 dex, Fe +0.2 dex, C +0.2 dex, N +0.2 dex. Every flux and continuum array was finite and exactly equal to its uncached counterpart at printed precision, on identical wavelength grids. The original general-purpose line-list case also matched its saved baseline. These checks do not establish internal bitwise identity or scientific validity across all stellar parameters.

A focused regression compares the real cached MOLOP routine to the original MOLOP using controlled broadening/profile stubs. It passes with bounds checking and invalid/zero/overflow traps enabled. It exercises repeated calls, list reload with changed gf at unchanged temperature, temperature changes, population/Doppler/stimulated-emission changes, separate lists and depths, changed dimensions, and a list that exceeds the allocation cap and must fall back. It is an isolated cache test, not a replacement for the physical comparisons.

Run it with `python3 benchmarks/test_molop_cache.py`, supplying --compiler or --reference-source as needed. The default reference comes from git HEAD and therefore assumes HEAD predates the cache; an explicit pre-cache source also works. The documented full-application legacy interpolation bounds declaration remains outside this patch.

## Experiments and timings

| Layout | Broad baseline median | Broad cached median | Result |
|---|---:|---:|---|
| 8,192-slot direct-mapped cache | 9.575 s | 9.672 s | -1.0% less time |
| Tables capped at 64 MiB | 9.635 s | 9.581 s | 0.6% less time |
| Tables capped at 256 MiB (selected) | 9.587 s | 8.430 s | 12.1% less time |

The sparse layout had no repeatable benefit and was discarded. The 64 MiB layout improved the short 100 Å case by about 3.8%, but the broad water list exceeded its cap and fell back. Raising the bound allowed caching that dominant list. Earlier short-window sparse timings include severe system-load outliers; raw values are retained rather than presented as a reliable gain. Each layout used freshly paired baselines; comparisons across separate experiments are not a controlled speedup measurement.

Abundance variations with the selected implementation (one measured pair each):

| Variation | Baseline | Cached |
|---|---:|---:|
| Fe_minus | 9.175 s | 8.301 s |
| Fe_plus | 9.752 s | 8.812 s |
| C_plus | 10.230 s | 8.570 s |
| N_plus | 9.244 s | 8.542 s |

## Scope and next work

GNU Fortran 11.1.0 with -fno-automatic -O2 on this Mac, warm local data, same repaired source except this cache, fixed 3500 K MARCS structure and 15500–16500 Å DR17 window. Compilation and input preparation are excluded from timings. The existing malformed atomic row omission is unchanged. Timings are measured, not extrapolated full-grid estimates. No additional compiler flags or physics shortcuts were introduced.

The next kernel experiment should examine molecular Voigt/profile accumulation or hydrogen profiles. This cache saves only one part of molecular-opacity work, and its moderate gain has a material memory cost. Broader atmosphere and wavelength validation is needed before production adoption.

main, the original project checkout, and remote PR #1 remain unchanged. No commit or push was made, and no line-list data are packaged.
