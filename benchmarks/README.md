# SYNSPEC performance experiments

This directory contains **non-production**, reproducible optimization experiments. The `main` branch and the existing `src/makefile` are intentionally unchanged.

## APOGEE optimization target

Prioritize full 15000–17000 Å syntheses for Teff 2800–5500 K, log g −0.5 to 3.0, and [Fe/H] −1 to +0.5. The user-specified bounds and current benchmark sampling control are recorded in [apogee_grid.json](apogee_grid.json); boundary models are included as stress tests. The production sampling was not specified, so 0.04 Å remains a benchmark control rather than a production requirement.

[APOGEE_GRID_RESULTS.md](APOGEE_GRID_RESULTS.md) documents the latest full-band comparisons, actual MARCS boundary/near-boundary models, and the accepted standard-broadening shortcut. It avoids a general GVDW call for the standard width expression, preserves the EXOMOL path, and adds no large per-line storage. The fully vectorized core/wing partition was slower and was not adopted.

Prepare a deck outside the checkout with `prepare_apogee_case.py --model /path/to/model.mod --linelist-dir /path/to/converted-dr17-lists --output /path/to/new-case`. The builder hashes inputs and records actual header parameters, including whether an input is a .filled substitute. It preserves the supplied atmosphere and abundances; it does not interpolate missing grid points.

## Current results and code changes

The optimization branch now includes three bounds corrections (HYDLIN population storage, STATE0 ionization initialization, and explicit-atom abundance factors), a molecular Voigt fast path, binary lookup of hydrogen profile intervals, and a standard molecular-broadening fast path. The large molecular exponential cache is **disabled by default**; `SYNSPEC_MOL_CACHE=1` enables the experimental 256 MiB capped cache.

On one 3500 K giant / 1,000 Å DR17 case, O2 reduced synthesis time by roughly 2–2.7× relative to no optimization. The low-memory profile changes then reduced O2 wall time by about 13%, with no material peak-memory increase. These are measured case-specific results, not full-grid guarantees. Four additional atmospheres and H-band/H-beta windows produced identical finite printed spectra in 24 baseline/improved comparison pairs; short-window timing changes ranged from about 4% faster to 3% slower.

See [BROADER_VALIDATION.md](BROADER_VALIDATION.md) for the latest matrix and re-profile, and [PROFILE_LOOP_RESULTS.md](PROFILE_LOOP_RESULTS.md) for kernel details. Earlier experiment reports describe the state at the time of each run. Main and production compiler defaults are not changed by this branch.

Run focused regression checks with:

```sh
python3 benchmarks/test_profile_loops.py
python3 benchmarks/test_hydrogen_search.py
```

Use `--compiler` when a compiler wrapper is needed (for example, to find the installed macOS SDK). The molecular test uses the committed, immutable pre-optimization routines in `reference_molop_voigt.f`, so its numerical oracle stays independent after commits. The tests use synthetic profile tables and initialized broadening parameters; full physical comparisons are documented separately. An all-application bounds-clean claim is not made: the legacy INTERP length-one dummy-array declarations remain a limitation of global bounds checks.

## Initial experiment: compiler flags

Use `python3 benchmarks/compare_builds.py --help` for all options. The script compiles `src/synspec54.f` into isolated temporary directories, avoiding the package installation targets or modification of `bin/`.

Example (compile all variants, then run each on the **same** input deck):

```bash
python3 benchmarks/compare_builds.py --build-only --output results/compiler
python3 benchmarks/compare_builds.py --input-dir /path/to/reproducible-case --output results/compiler --repeats 3
```

The case directory should contain all inputs expected by an ordinary SYNSPEC run (including any relative-path files). It is copied separately for every measurement. You may pass `--command` to run an explicit shell command instead of executing the compiled binary directly, e.g. to redirect units and standard input: `--command './synspec54 < fort.5'`. Run the command from the copied case directory. The script sets `OMP_NUM_THREADS=1` by default.

The first run always compiles a baseline version. The other configurations are `O2`, `O3`, `native`, and `fast-math`. The `native` build is CPU-specific and should not be redistributed to machines with different instruction sets. `fast-math` may change numerical results and must not be accepted without stringent physical validation.

For each successful run the script records wall-clock time, exit status, stdout/stderr, and SHA-256 hashes of output files. A changed hash is **not** a quantitative spectrum comparison: it indicates that a dedicated wavelength-aligned flux/continuum comparison is needed. The benchmark harness never claims results for a variant that failed to compile or run.

### Scope and safety

- The current `src/makefile` has `OPT=-fno-automatic` and adds `-mcmodel=medium` on Linux. These flags are retained for all initial gfortran variants.
- Compare complete, realistic H-band syntheses before estimating end-to-end speedup.
- Measure output differences at spectrum/continuum level, including relative error, changes in line cores and equivalent widths. Set scientific tolerances before adopting `-ffast-math`.
- Later phases: instrument `RESOLV`, `LINOP`, `MOLOP`, `VOIGTK`, `RTE` and `RTEDFE`; investigate line-opacity and SIMD changes as separate commits.
- Record processor, compiler version, options, code commit, atmosphere, linelists, wavelength interval, runtime, and spectra checksums.
- Current measured results are linked above. Earlier Voigt microbenchmark numbers from a separate machine are not full-SYNSPEC benchmarks.


## Vectorization comparison on the improved code

[Vectorization results](VECTORIZATION_RESULTS.md) compare strict-math O2, O3, native tuning, and an O3/no-vectorization control. All 66 non-reference output comparisons match O2 at printed precision. The timings are noisy and workload-dependent, so production compiler defaults remain unchanged. The principal molecular profile loop still does not vectorize because of control flow.

To reproduce with externally prepared input decks:

```sh
python3 benchmarks/compare_vectorization.py --case-manifest /path/to/cases.json --wide-case /path/to/wide-case --output results/vectorization
python3 benchmarks/analyze_vectorization.py --results results/vectorization
```

The manifest is a JSON array containing `name` and `case` directory keys. NumPy is needed for analysis. On this Mac, `--compiler` selects the local SDK wrapper and `--runtime-library-dir` selects the matching Fortran runtime directory. Use a new output directory for a fresh run. `--reuse-builds` only accepts binaries whose recorded source, compiler version, platform and flags match; it does not overwrite existing case directories.
