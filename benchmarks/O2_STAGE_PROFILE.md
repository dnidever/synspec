# Repaired O2 stage profile and abundance batch — November 13, 2026

**The next optimization target is molecular-line opacity (MOLOP), followed by hydrogen-line profiles (HYDLIN).** Initial line-list preparation and radiative transfer are relatively small in this workload.

## Measured workload

Same repaired GNU Fortran 11.1.0 O2 source as the compiler comparison, with -fno-automatic. 3500 K MARCS giant, log g=0, [M/H]=−0.5, microturbulence=2 km/s; LTE plane-parallel synthesis over 15500–16500 Å with DR17 atomic, giant molecular, and BC8.5V water lists. The previously documented single malformed atomic record remains omitted. The atmosphere structure stays fixed during the abundance variations.

Seven input cases were each run with both ordinary O2 and an instrumented O2 build: three reference repeats, Fe −0.2 dex, Fe +0.2 dex, C +0.2 dex, and N +0.2 dex relative to the corresponding abundance in the input deck. Variations use ICHEMC=1 / fort.56, exercising the existing abundance-change path. All start as independent fresh processes. This is not an implemented persistent-synthesis benchmark.

The uninstrumented reference median was **10.961 s**; instrumented reference median **10.404 s**. The five-spectrum batch (one reference plus four abundance variations) took **52.176 s** in total, summing individual uninstrumented process wall times. Case creation and Python-side validation are excluded. The separate measured runs, not the sampled run, supply stage timings.

## Measured CPU breakdown

Mean CPU seconds across three instrumented reference runs. Percentages use the mean timed-region CPU total. The rows below are disjoint: the molecular and hydrogen entries include their untimed descendants; separately timed descendants are subtracted from their parents.

| Work | CPU seconds | CPU share |
|---|---:|---:|
| Molecular-line opacity (MOLOP) | 4.375 | 43.8% |
| Hydrogen-line opacity (HYDLIN) | 2.572 | 25.8% |
| Tabulated continuum opacity (OPACON) | 0.914 | 9.1% |
| Molecular/ionization equilibrium (MOLEQ) | 0.615 | 6.2% |
| Other atomic-line opacity (LINOP) | 0.141 | 1.4% |
| Radiative-transfer driver (RTECD) | 0.051 | 0.5% |
| Spectrum output (OUTPRI) | 0.035 | 0.4% |
| Initial atomic-list preparation (INILIN) | 0.023 | 0.2% |
| Initial molecular-list preparation (INMOLI) | 0.123 | 1.2% |
| Remaining setup, continuum work, bookkeeping and unclassified work | 1.140 | 11.4% |

Total mean instrumented CPU time: **9.989 s**. Full per-routine inclusive and exclusive CPU/wall timings and call counts are in summary.json and the individual stage-times.csv files. OPACON does not include every continuum operation; some remain in OPAC. INILIN/INMOLI include initial reading/selection, not all subsequent per-chunk line selection. Do not add inclusive RESOLV, OPAC, MOLOP, and HYDLIN times together.

## Sampling evidence

macOS sample collected 8,175 main-thread stack observations from a separate instrumented O2 build with debug symbols and frame pointers. The requested interval was 1 ms, for up to 30 s, ending when synthesis exited. Largest exclusive leaf counts: MOLOP 2,112; exp 1,688; HYDLIN 1,246; VOIGTK 612; OPACON 430; RUSSEL 292; GVDW 233; LINEQS 142. These are stack observations, not exact CPU seconds. exp and VOIGTK are descendants of the opacity stages, so their counts must not be added to the stage percentages.

The call tree attributes 716 exp leaf observations to the MOLOP line-strength calculation, and 584 to a hydrogen-profile calculation. Molecular work also spends visible time in GVDW and VOIGTK. This supports examining repeated exponential/broadening calculations and profile accumulation before changing the transfer solver. Raw symbols and generated-source line locations are in sample-O2.txt and instrumented.f.

## Reuse opportunities and limits

Main-program setup before RESOLV averages 2.339 CPU s, or 23.4% of the timed CPU total. An idealized removal of *all* that setup would permit at most **1.31×** speedup for this reference case, ignoring implementation overhead. This is a CPU-time upper bound, not a measured persistent-process speedup. Much of the setup includes equilibrium and continuum calculations that depend on abundance and cannot simply be reused.

Initial INILIN + INMOLI accounts for about 1.5% of CPU time. Reusing only these entire initial preparation phases would therefore have a small ceiling in this already warm-cache DR17 workload; caching raw files alone removes only a portion. Later per-chunk line selection has not been isolated as a separate timer. Different line lists, shorter wavelength intervals, or many atmosphere models can change this balance.

A hypothetical doubling of MOLOP alone would produce about 1.28× end-to-end CPU speedup under this measured breakdown; doubling both MOLOP and HYDLIN would yield about 1.53×. These are Amdahl-law estimates, not measured optimization results.

## Validation and timing limitations

All seven instrumented/uninstrumented pairs had identical finite printed flux and continuum arrays. Every abundance-variation spectrum differed from the reference, confirming the inputs took effect. The sampled reference spectrum also matched the ordinary reference exactly. These checks establish noninterference on these cases at text-output precision, not internal bitwise equality or universal scientific accuracy.

Timers use SYSTEM_CLOCK and CPU_TIME around selected CALL sites with a nesting stack. Inclusive time covers descendants; exclusive time subtracts other instrumented descendants only. Hot scalar routines such as exp and VOIGTK were sampled instead of timed on every call. Timer calls have nonzero overhead, particularly system CPU accounting. The timed median happened to be lower than the ordinary median; machine-load variation prevents using that difference as an overhead measurement or claiming zero overhead. Frame-pointer settings differ between the two builds. CPU shares are diagnostic approximations, and process wall time includes startup outside the Fortran timed region. No instrumented performance gain is claimed.

## Next implementation experiment

Optimize MOLOP first: quantify repeated per-line/depth exponential and damping computations across frequency chunks; test a bounded cache or hoisting of demonstrably invariant factors, preserving abundance/equilibrium invalidation. Then inspect the expensive HYDLIN profile loop. Each change should run against both reference and abundance-varied cases, with flux/continuum comparisons and uninstrumented wall-time measurements. No such kernel optimization or persistent-state cache has been implemented in this step.

Production Fortran was not modified by instrumentation: a generated source copy lives under work/stage-profile. main and the remote PR remain unchanged; no large line lists are included in the reproduction archive.
