# Initialization bounds corrections — November 13, 2026

Two changes were made in the separate local checkout on `perf/compiler-benchmarks`, retaining the earlier HYDLIN population-array correction.

## ENEV / INPOT

STATE0 now loads and classifies the `MI1` available ionization energies instead of looping to nine and reading beyond ENEV's eight columns. INPOT's final entry is explicitly assigned the category of the last available threshold: `INPOT(I,MION0)=INPOT(I,MI1)`. No extra ionization energy or array storage is invented.

STATE uses the final INPOT entry when computing the highest ion stage, while its corresponding energy access is ENEV(I,J-1), the eighth energy for stage nine. Reusing the last available threshold category makes that final entry deterministic and avoids reading neighboring COMMON-block data. Categories for the first eight entries are unchanged. High-temperature/non-LTE validation of this previously undefined final-entry behavior remains outside the tested cool-giant LTE cases.

## RELAB

STATE0 initializes RELAB over MATEX, the allocated explicit-atom count, instead of NATOMS, the chemical-element count. In this workload those counts are 30 and 99 respectively. Existing abundance-change writes use IATEX, the explicit-atom mapping, and the population reader uses IATM. The correction therefore preserves the existing explicit-atom indexing without expanding COMMON-block storage.

## Verification

The full bounds-checked executable passes both repaired initialization sections. It subsequently stops in INTERP, whose dummy arguments are declared X(1), Y(1), XX(1), YY(1) despite callers supplying longer arrays. This is a separate legacy declaration issue, not evidence that the actual caller allocation has length one. No interpolation changes were included.

To check initialization without suppressing its checks or changing interpolation, STATE0 was extracted unchanged into a separate compilation unit and compiled with `-O2 -g -fcheck=all -ffpe-trap=invalid,zero,overflow`. The remainder was compiled with the same options except `-fcheck=all`. This executable completed the physical synthesis cases. Thus bounds-checking coverage is specifically STATE0; the application is not claimed to be globally bounds-clean.

O0, O2 (STATE0 checked), and O3 builds were compared against the saved original finite baseline for both DR17 and general-purpose lists. All six runs completed with floating-point traps enabled; flux and continuum outputs were finite and identical at printed precision. Raw run measurements and hashes are in initialization-fix-verification.json. These single-run diagnostic timings are not revised benchmark medians.

An additional ICHEMC=1 case changes the Fe abundance through fort.56, exercising the RELAB abundance-change path. Both the previous HYDLIN-only repair and the new repair produced identical finite spectra, and the modified-abundance spectrum differs from the ordinary case. Its detailed result is recorded in initialization-abundance-check.json.

## Scope

The patch contains all three demonstrated array corrections: ENEV/INPOT initialization, RELAB initialization, and the earlier HYDLIN PJ(50) change. The file initialization-only-fixes.patch isolates just the two changes from this turn. No production source outside the separate workspace checkout was modified. main, the remote branch, and draft PR #1 remain unchanged; no commit, push, or large-data upload was performed. The DR17 malformed atomic damping entry remains a separate data issue.
