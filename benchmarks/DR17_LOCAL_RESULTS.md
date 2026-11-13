# APOGEE DR17 benchmark update

Use the official SDSS DR17 lists for the APOGEE benchmark: https://dr17.sdss.org/sas/dr17/apogee/spectro/speclib/linelists/synspec/

Downloaded the atomic list (5,331,110 bytes), giant molecular list (3,536,568 bytes), and BC8.5V water list (34,360,024 bytes), a total of 43,227,702 bytes. All three match the archive's published SHA-1 checksums. SHA-256 manifests are supplied separately. The solar molecular list is an alternative to the giant list, not an additional component. The BC9.0V/BC9.5V water lists were not downloaded or evaluated.

The SDSS README documents isotope ratios embedded in molecular log(gf), including 12C/13C = 7 for the giant list. Thus these are a more appropriate APOGEE benchmark input than the general-purpose solar molecular list previously available. Atmospheric abundance and isotope assumptions still need review for scientific production.

## Conversion issues found

The repository's converter uses a 102-character first-record buffer, but the atomic records have 103 characters. This truncates the final integer and causes incorrect format detection (10 reals + 1 integer instead of 10 reals + 7 integers). A separate local converter copy expands the buffer to 4096 characters; production source is unchanged. This correction is also saved as an uncommitted experimental file under benchmarks/list2bin_dr17.f on the local optimization branch.

The atomic source additionally contains a nonnumeric `s` radiative-damping value on row 33407, at 1630.3785 nm, species 60.01. The converter skips this record through its existing error handler. Therefore the current binary has 51,260 of the documented 51,261 atomic records. **The current benchmark is provisional; resolve this field before production use.** No numeric replacement was invented. Every remaining converted record was compared byte-for-byte to values parsed from the original text. The molecular and water binaries contain 49,119 and 483,944 records, respectively; every record was verified.

## Measured benchmark

Same source commit, compiler, atmosphere, and 16000–16100 Å settings as the earlier report, replacing four general lists with these three DR17 lists. Separate run directories preserve symlinks to local data. The baseline took 6.784, 6.619, and 6.525 seconds, median **6.619 seconds**, and produced 2,839 finite flux values each time. This is a different scientific workload, so the difference from the earlier 9.973-second median is not a compiler speedup. The O2 runs took 2.933, 4.281, and 3.021 seconds but produced NaN at every flux sample in all three runs. The optimized-build failure persists with the DR17 lists; these times do not establish a usable speedup. Raw timing/validation results are in dr17-timings.json. No new routine-level profile was collected for DR17.

Data remain in work/dr17-linelists, outside the Git checkout. main and PR #1 remain unchanged. No data were committed or uploaded.
