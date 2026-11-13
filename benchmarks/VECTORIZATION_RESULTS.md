# Strict-math vectorization/compiler comparison — November 13, 2026

## Setup

Source is the repaired, low-memory profile-loop implementation at commit 1f623b7. The large exponential cache was disabled (SYNSPEC_MOL_CACHE=0), and no fast-math options were used. All builds retained -fno-automatic and used GNU Fortran 11.1.0 on this Intel Mac. Compiler diagnostics were enabled with -fopt-info-vec-all. They do not instrument the running program.

Variants: O2; O3; O3 -march=native; and O3 -fno-tree-vectorize as a diagnostic control. The control disables loop and basic-block/SLP vectorization together. Generic O3 and native tuning differ in additional CPU-specific choices, not vector width alone. Production compiler defaults were not modified.

The 3500 K / 15500–16500 Å DR17 reference used one excluded warm-up and five measured repeats per variant. The four-atmosphere, two-window matrix used one excluded warm-up and two measured repeats of O2/O3/native, with the general-purpose local lists. Order was shuffled within rounds, runs were sequential, and builds completed before measurement. Wall time includes process startup and ordinary synthesis I/O but excludes compilation, case preparation and validation. Process CPU time was also recorded.

## Wide-reference measured medians

| Variant | Median wall time | Min–max | Speedup over O2 | Paired CPU speedup over O2 |
|---|---:|---:|---:|---:|
| O2 | 10.118 s | 8.927–12.115 s | 1.000× | 1.000× |
| O3 | 8.972 s | 8.266–12.196 s | 1.128× | 1.085× |
| native | 9.461 s | 8.241–10.683 s | 1.069× | 1.080× |
| O3-no-vector | 11.094 s | 8.948–12.210 s | 0.912× | 0.903× |

Elapsed times varied substantially. A late load-average observation was about 7.7 on eight logical CPUs, consistent with contention; it is not a CPU-utilization measurement and does not explain every difference. Ratios of medians and medians of paired ratios are both in summary.json; neither is a confidence interval. Small differences should not be promoted to a universal compiler recommendation. The diagnostic control tests aggregate vectorization effects throughout the application, not just the molecular wavelength loop.

## Decision

Keep O2 as the current portable baseline. O3 looks promising on the wide reference, but the broad-case ranking is mixed and this run has substantial timing noise. Native tuning has no consistent advantage. No production flag change is justified by these results alone. A quieter repeat would be needed before treating the additional wide-reference difference as a reliable gain.

The compiler reports identify a concrete next source-level experiment: split the molecular small-damping core and far-wing paths into branch-free wavelength ranges, then check whether they vectorize without changing formulas or boundary decisions. That experiment has not been implemented in this step.

## Broader matrix

| Case | O2 median | O3 median | Native median |
|---|---:|---:|---:|
| arcturus-Hband | 5.907 s | 7.139 s | 7.260 s |
| arcturus-Hbeta | 4.945 s | 4.956 s | 5.600 s |
| solar_dwarf-Hband | 7.039 s | 5.586 s | 6.385 s |
| solar_dwarf-Hbeta | 4.459 s | 3.775 s | 3.771 s |
| cool_giant-Hband | 5.250 s | 4.781 s | 4.676 s |
| cool_giant-Hbeta | 8.384 s | 8.960 s | 8.239 s |
| warm_giant-Hband | 3.886 s | 3.940 s | 3.802 s |
| warm_giant-Hbeta | 3.092 s | 3.062 s | 3.039 s |

These are two-repeat medians for short, setup-heavy windows. They expand correctness coverage and show workload dependence; they are not precise rankings of close configurations. Atmosphere/window definitions and line-data limitations are the same as BROADER_VALIDATION.md.

## What actually vectorized

GCC 11 emits no vectorization successes for the O2 or O3-no-vector builds in these reports. O3 reports 16-byte vectors in selected loops; native tuning additionally reports 32-byte vectors. These are compiler claims, not measurements that every generated loop version executed.

- MOLOP: the final emissivity-accumulation loop vectorizes. Cache-table initialization/bookkeeping also has vectorized code, but the optional cache was disabled during these runs.
- The expensive small-damping wavelength loop still does **not** vectorize: GCC reports control flow in the loop. The larger-damping loop remains blocked as well.
- HYDLIN: population setup and profile copying have vectorized paths. The main lookup/profile loops still contain control flow and scalar transcendental work.

Thus the earlier hypothesis is only partly realized: removing a scalar Voigt call helped scalar execution, but O3 alone has not vectorized the principal molecular profile-accumulation loop. A gain from O3 cannot be attributed specifically to vectorizing that hot loop. See compiler-analysis.json for exact source lines and selected missed-vectorization messages, and the compressed reproduction archive for complete compiler logs.

## Spectrum checks

Across 96 successful physical executions, 66 non-O2 outputs were compared with the same-case/same-round O2 reference, including warm-ups. **66/66** have finite flux and continuum arrays, identical wavelength grids, and exactly identical printed values. Full SHA-256 values and numerical metrics are in spectrum-checks.json. Finite/relative-error screening also exists, but exact printed equality is reported separately.

Equality at output precision is not internal bitwise equivalence or universal physical validation. No claim is made for NLTE, spherical transfer, a full atmosphere grid, other compilers, or other CPUs. Native binaries are specific to the build CPU.

## Reproduction and limitations

Use benchmarks/compare_vectorization.py with a JSON case manifest, an external prepared wide-case directory, and a fresh output directory. The manifest is an array of objects containing name and case (directory) keys; other provenance fields may be included. The runner preserves symlinks to large local input files.

A local compiler wrapper supplied the macOS SDK library path, and --runtime-library-dir selected the installed Fortran runtime. The first launch omitted that runtime path and failed in the dynamic loader before synthesis; it is excluded from all measurements. Existing matching binaries were reused after correcting the runtime environment, without recompilation. No failed launch is counted as a synthesis speed measurement.

No source kernels, production build flags, or main-branch files were changed for this experiment. Large line lists and model/opacity directories remain external.
