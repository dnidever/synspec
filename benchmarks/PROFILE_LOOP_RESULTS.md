# Low-memory profile-loop improvements — November 13, 2026

**Keep the molecular Voigt fast path and hydrogen binary lookup; leave the large molecular-factor cache off by default.**

## Measured result

| Case | Repaired O2 baseline | New default O2 | Less time | Repeats |
|---|---:|---:|---:|---:|
| 100A | 2.961 s | 2.712 s | 8.4% | 3 |
| 1000A | 9.485 s | 8.282 s | 12.7% | 5 |

Peak resident memory for the 1,000 Å reference was **76.6 MiB baseline and 76.5 MiB new default**, measured in separate fresh processes. There is no material measured memory increase. This replaces the previous default cache experiment, which increased peak memory to about 220 MiB.

The initial independent three-repeat experiment measured these contributions (each against the same baseline series):

| Change | Median | Less time |
|---|---:|---:|
| voigt | 9.141 s | 4.0% |
| hydrogen | 8.640 s | 9.2% |
| combined | 8.299 s | 12.8% |

These contributions are measured alternatives, not additive percentages. Final combined wide-window statistics incorporate two additional paired runs for five repeats total. Warm-ups were excluded, build order was shuffled, and all executions were sequential. Raw process CPU times also support the improvement: final median paired CPU speedup is approximately 1.15×. Compilation and case preparation are excluded. Warm-cache measurements on one Mac are not controlled cold-start or full-grid measurements.

## Changes

1. **MOLOP small-damping path:** AGAM is fixed for a line/depth. When AGAM < 0.2, select that Voigt branch once outside the frequency loop, then evaluate the existing lookup polynomial or far-wing expression directly. Larger damping still calls VOIGTK. Constants, arithmetic order, frequency limits, and accumulation order are preserved. The duplicate fast path is deliberately checked against the original VOIGTK implementation in regression tests.

2. **HYDLIN lookup:** replace the per-wavelength linear scan of the ordered hydrogen profile table with a binary lower-bound search. The first upper endpoint >= AL is selected, with the same low/high tail interval clamping and equality convention as the old scan. The interpolation expression and profile exponentials remain unchanged. This relies on ordered profile grids, as supplied by the tested tables; broader table/version validation remains prudent.

3. **Optional cache:** the previously tested factor cache allocates no tables by default. Enable it explicitly with `SYNSPEC_MOL_CACHE=1`. It retains its 256 MiB table cap and original-path fallback, plus temperature/list/depth/line-count invalidation. Values other than exactly 1 leave it disabled. The environment setting is read once per process. Both enabled and disabled modes were checked; this step does not claim an additional measured speedup for combining the cache with these loops.

## Validation

All physical comparison pairs produced finite flux and continuum arrays exactly equal at printed precision on identical wavelength grids. Coverage includes the independent four-build 1,000 Å comparison (including warm-ups), two further wide pairs, the 100 Å case, Fe ±0.2 dex, C +0.2 dex, and N +0.2 dex. An additional general-purpose-list case matched its saved baseline. A full 1,000 Å run with the cache explicitly enabled also matched the baseline. See analysis.json and extra-validation.json.

The hydrogen search regression checks **16,551** intervals, exact endpoints, neighboring values and tails over table lengths 2–55 against the original linear rule. The molecular regression uses the real MOLOP and original VOIGTK with synthetic lookup tables and a controlled broadening stub, exercises both damping branches and far wings, and repeats cache-invalidation/fallback tests. Both tests pass under bounds checking; the molecular test additionally enables floating-point invalid/zero/overflow traps. These focused checks do not claim that the full application is globally bounds-clean.

Reproduction commands from the checkout:

```sh
python3 benchmarks/test_profile_loops.py --compiler /path/to/gfortran
python3 benchmarks/test_hydrogen_search.py --source src/synspec54.f --compiler /path/to/gfortran
```

The molecular test defaults to git HEAD as the pre-change reference. After committing these changes, provide `--reference /path/to/pre-change-source.f` to keep an independent old reference. The source archive supplies the saved pre-cache, repaired baseline. The legacy `test_molop_cache.py` entry point delegates to the updated molecular regression.

## Scope

Same 3500 K MARCS giant, log g=0, [M/H]=−0.5, microturbulence=2 km/s, LTE plane-parallel setup and DR17 lists used earlier. Wide interval 15500–16500 Å; short interval 16000–16100 Å. GNU Fortran 11.1.0, -fno-automatic -O2. Abundance variations change opacity/equilibrium inputs with a fixed atmosphere structure. The known malformed atomic-row omission and separate global bounds-checking limitations remain unchanged. Printed-array equality does not establish internal bitwise equality or universal physical validation.

The local optimization branch contains these changes and the earlier three array corrections. main, the original project checkout, and remote PR #1 remain unchanged; no commit, push, or large-data upload was performed. Production compiler defaults were not changed.
