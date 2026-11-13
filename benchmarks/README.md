# SYNSPEC performance experiments

This directory contains **non-production**, reproducible optimization experiments. The `main` branch and the existing `src/makefile` are intentionally unchanged.

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
- No benchmark results are asserted in this repository yet. Earlier Voigt microbenchmark numbers from a separate machine are *not* full-SYNSPEC benchmarks.
