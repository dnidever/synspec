# Local SYNSPEC benchmark — November 12, 2026

Source: dnidever/synspec, branch `perf/compiler-benchmarks`, commit `dceaedda1871cdf92859081d5ee607fdf17f552e`. Draft PR #1 remains unchanged. The separate checkout is in `work/synspec`; the original checkout and main were not modified. No data were committed or uploaded.

## Available data

The requested local directory is readable. It contains general-purpose lists, not the APOGEE-specific `synspec.20200921nlte.19`, `synspec.20200921nltegiant.20`, or `synspec.h2o-BC8.5V.20` files.

| File | Bytes | Validated binary records | Record payload |
|---|---:|---:|---|
| gfATO.19.11 | 262,270,664 | 2,260,954 | 10 float64 + 7 int32 |
| gfMOLsun.20.11 | 241,955,480 | 6,048,887 | 4 float64 |
| gfTiO.20.11 | 333,292,880 | 8,332,322 | 4 float64 |
| H2O-8.20.11 | 786,739,440 | 9,834,243 | 9 float64 |

All are little-endian Fortran sequential unformatted files, with matching 4-byte leading/trailing record lengths validated throughout. Total: 1,624,258,464 bytes. The `.00` files are small ASCII summaries with source counts and endpoint records, not full ASCII line lists. The atomic `.00` reports 2,309,497 lines, 48,543 more than the actual binary; cause unknown. The other counts agree. Checksums and wavelength endpoints in nm are in `linelist-inventory.json`.

## Reproducible case

Local MARCS atmosphere: `s3500_g+0.0_m1.0_t02_x3_z-0.50_a+0.00_c+0.00_n+0.00_o+0.00_r+0.00_s+0.00.mod`, from `/Users/nidever/projects/synple/models`. Teff 3500 K, log g 0.0, [M/H] −0.5, microturbulence 2 km/s, 56 depths. LTE plane-parallel transfer using the MARCS structure; this is not a validation of spherical transfer. Interval 16000–16100 Å, requested sampling control 0.04 Å, cutoff 250 Å, selection strength 1e-4. All four lists were linked read-only from their original locations. This is a cool-giant H-band workload using the available solar molecular list, not a scientifically validated APOGEE giant line-list configuration.

The branch's Python input writers generated fort.5, fort.8, fort.55 and tas. The local runner preserves symlinks to avoid repeatedly copying 1.6 GB of data; the committed comparison harness normally dereferences them. No Fortran source changes were made.

## Measured wall times

Mac x86_64, macOS 13.7.2; Homebrew GNU Fortran 11.1.0. Every build uses `-fno-automatic`, with optimization added as below, one OpenMP thread. Three sequential runs per build, baseline first, then O2, then O3. Cache state and competing machine load were not controlled. Times exclude compilation and input preparation.

| Build | Measured seconds | Median | Flux validation |
|---|---|---:|---|
| Baseline (no -O flag) | 9.973, 9.381, 10.056 | 9.973 | 2838 finite values, identical hashes across runs |
| -O2 | 4.955, 4.847, 4.916 | 4.881* | all 2838 flux values NaN, every run |
| -O3 | 5.035, 4.475, 4.611 | 4.611 | all 2838 flux values NaN, every run |

*O2 median excludes the first run, during which a sampling attempt was rejected. All nine runs returned exit code 0. All continuum outputs are finite and have identical hashes across builds. The optimized flux outputs are invalid: these timings **do not establish a scientifically usable speedup**. Root cause is not yet diagnosed. Finite baseline output alone does not establish physical accuracy.

Initial build attempts failed to locate libSystem. Adding `-L/Library/Developer/CommandLineTools/SDKs/MacOSX.sdk/usr/lib` resolved linking. Execution also requires `DYLD_LIBRARY_PATH=/usr/local/Cellar/gcc/11.1.0_1/lib/gcc/11` because the compiler's embedded runtime path is stale. These adjustments are local to this workspace.

## Routine-level sampling

GNU gprof is absent. macOS `/usr/bin/sample` succeeded with expanded process access. A separate baseline run was sampled at requested 1 ms intervals for up to 12 seconds, ending when the process exited. 9,417 main-thread stack observations were collected. The largest exclusive leaf counts were:

| Routine | Leaf observations | Share of observed stacks |
|---|---:|---:|
| RUSSEL | 953 | 10.1% |
| LINEQS | 914 | 9.7% |
| exp (math library) | 789 | 8.4% |
| OPACON | 739 | 7.8% |
| HYDLIN | 595 | 6.3% |
| MOLOP | 422 | 4.5% |
| SIGK | 367 | 3.9% |
| VOIGTK | 199 | 2.1% |

These are sampled stack shares, **not measured per-routine seconds or exact phase timings**. Some stacks are unresolved or mis-unwound, so inclusive phase totals are not reliable. Do not sum nested calls. The raw full report is `sample-baseline.txt`; the short O2 sample is retained only for diagnosis of the invalid optimized build. No extrapolated full-H-band timings are asserted.

## Next experiment

Diagnose the O2/O3 NaNs with floating-point traps and a debug build before accepting compiler optimizations. Profile repeated finite-output runs with debug symbols and frame pointers, then add explicit CPU/wall phase timers if robust inclusive timings are required. Obtain the three APOGEE-specific line lists named above (or locate an existing copy) before claiming an APOGEE-production benchmark. Preserve the current general-list case as a regression reproducer.

## Saved artifacts

`timings.json` and `validation.json` hold all runs and output hashes; `comparison.json` summarizes the first run of each build. The reproduction archive contains the small input deck, scripts, metadata, build logs, and representative spectra. Large line lists and auxiliary opacity data remain external; the original local paths must remain available to reproduce the run.
