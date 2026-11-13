# Repaired SYNSPEC compiler benchmark — November 13, 2026

Measured on this Mac with GNU Fortran (Homebrew GCC 11.1.0_1) 11.1.0. All builds use the same repaired source, `-fno-automatic`, one OpenMP thread, and the same local runtime libraries. No profiling, debug checks, or floating-point traps were enabled in these timed builds. The source SHA-256 is `a885cbb91a24f88d635827803ff93655b95b1693d3e61ff36587d6d73bcbd079`. The tested source is saved alongside this report.

The three array corrections are included. These results compare compiler options on repaired code; they do not compare corrected code against a corrupt optimized executable.

| Build | 100 Å median | Speedup | 1,000 Å median | Speedup |
|---|---:|---:|---:|---:|
| No optimization | 6.533 s | 1.00× | 25.965 s | 1.00× |
| `-O2` | 3.024 s | 2.16× | 9.710 s | 2.67× |
| `-O3` | 3.181 s | 2.05× | 11.139 s | 2.33× |
| `-O3 -march=native` | 2.891 s | 2.26× | 11.032 s | 2.35× |
| Native + `-flto` | 2.904 s | 2.25× | 10.161 s | 2.56× |
| Native + `-ffast-math` | 3.908 s | 1.67× | 13.240 s | 1.96× |

## Method

100 Å: 16000–16100 Å, five measured repeats per build. 1,000 Å: 15500–16500 Å, three measured repeats per build. Every workload/build has one untimed warm-up excluded from medians. Build order within each measurement round was shuffled using seed 74219. Executables ran sequentially; compilation completed before measurement. Times include process startup, initialization, line-list reading, synthesis, and ordinary file output; they exclude compilation, case preparation, spectrum analysis, and warm-ups. These are warm-cache throughput measurements, not controlled cold-start tests. System load and CPU frequency were not controlled. Small differences between nearby configurations should not be treated as statistically established.

The workload is the existing LTE plane-parallel synthesis of a 3500 K MARCS giant, log g=0, [M/H]=−0.5, microturbulence=2 km/s, with SDSS DR17 atomic, giant molecular, and BC8.5V water lists. Sampling control is 0.04 Å; the actual adaptive output grids contain 2,839 and 27,842 samples. The atomic binary still omits the one malformed source record documented earlier. These are two windows of one atmosphere, not full-grid or full-H-band validation.

## Spectrum validation

Every measured run produced finite flux and continuum, identical wavelength grids, and passed a numerical screen of 1e-5 relative flux/continuum difference and 1e-5 absolute continuum-normalized flux difference.
In fact, every printed flux value matched the baseline exactly in both windows. See spectrum-checks.json for per-run continuum comparisons, normalized-flux metrics, and SHA-256 hashes. Equality at text-output precision does not establish equality of internal floating-point calculations or validate all stellar parameters. Fast-math remains an experimental choice; passing these cases is not general scientific approval.

## Interpretation

Optimization provides a substantial measured gain. Native tuning and link-time optimization must be compared against ordinary O2/O3 using the timing spread, rather than assuming that more flags are faster. Native binaries target this CPU and should not be distributed indiscriminately. Fast-math was slower than ordinary O2/O3 in these measurements and offers no demonstrated reason for adoption here.

All raw wall and process-CPU times are saved in runs.json, including excluded warm-ups; summary.csv includes min/max ranges. No extrapolated full-grid timings are asserted. main, production makefile, and the remote PR remain unchanged. No large line-list files were committed or packaged.
