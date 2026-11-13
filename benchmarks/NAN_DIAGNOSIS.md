# SYNSPEC NaN diagnosis — November 13, 2026

The NaNs are caused by an array overrun in HYDLIN, not by the DR17 line-list substitution.

## Reproducer and causal evidence

On the existing 3500 K MARCS giant / 16000–16100 Å DR17 case, the original O2 executable returns exit code 0 but writes NaN for all 2,839 emergent-flux samples. The continuum remains finite.

An O2 build with debug symbols and `-ffpe-trap=invalid,zero,overflow` stops in `exp`, called by HYDLIN at original source line 5574. The call chain is MAIN → RESOLV → OPAC → HYDLIN → exp. The failing expression is:

```
EMTRA=PJ(J)*WNHINT(I,ID)*II*XJJ*EXP(CPJ*(XII-XJJ)*T1)
```

HYDLIN declares `PJ(40)` at original line 5434, then initializes `PJ(IL)` for `IL=1,50` at original lines 5488–5492. Ten writes therefore exceed the array. Local diagnostics show temperature variables already corrupted when the exponential fails: at depth 36, I=3, J=5, T is approximately 0.0686 and T1 approximately 0.0700, yielding an exponential argument far beyond the representable range. The benchmark atmosphere temperatures are thousands of kelvin. Compiler optimization changes local storage layout, exposing the corruption differently.

The minimal correction increases only HYDLIN's population array from 40 to 50 elements, preserving its existing initialization and all physics formulas. The other PJ(40) declaration belongs to a separate routine whose initialization stops at 40 and was left unchanged. With this one correction, the O2 run completes with floating-point traps enabled, produces finite values, and matches the original unoptimized baseline exactly at every printed wavelength and flux sample. The continuum matches as well. All six verification combinations passed: O0, O2, and O3 on both DR17 and general-list cases, with floating-point traps enabled and exact equality of printed flux and continuum values. See nan-fix-verification.json for the measurements.

This is equality at the precision of the text output; it is not a claim of bitwise equality of internal calculations or full physical validation. All measurements use the same GNU Fortran 11.1.0 runtime and SDK adjustments documented in the earlier report. Verification run timings include debug/trap options and are single runs, not revised performance medians.

## Changes and regression protection

The separate local checkout on `perf/compiler-benchmarks` has the one-array correction in `src/synspec54.f`, plus `benchmarks/check_spectrum.py`. The checker rejects NaN/infinity, nonmatching or nonmonotonic grids, and flux disagreements beyond explicit tolerances. It was exercised against the repaired spectrum (passes, zero printed difference) and the original O2 NaN spectrum (fails). It requires NumPy.

`main`, the original checkout, and the remote PR are unchanged. No commit or push was made. No line-list data are included in the patch.

## Other defects exposed by bounds checking

Full `-fcheck=all` stops earlier in initialization, independently of the HYDLIN reproducer:

- STATE0 original line 1716 reads `ENEV(I,9)` even though its second dimension is MI1=8. The loop fills INPOT through nine ionization entries. A temporary diagnostic guard bypassed this read only to continue checking; no speculative ninth-ion handling was put into the final patch.
- After that diagnostic guard, STATE0 original line 1778 writes `RELAB(JJ,ID)` through NATOMS=99, although RELAB's first dimension is MATEX=30. This prevents claiming that the entire application is bounds-clean.

These are separate follow-up defects. The minimal HYDLIN correction resolves the observed optimized-build NaNs without relying on those diagnostic changes. The remaining defects deserve their own fixes and scientific review. The DR17 malformed atomic damping entry and its single-record omission also remain as previously documented.
