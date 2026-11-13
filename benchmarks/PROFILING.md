# Profiling the actual SYNSPEC code

Run on a representative synthesis input case, ideally a MARCS giant and a 100 Å H-band window. For full-grid throughput, repeat with a broad H-band interval and a line-rich, molecularly complex atmosphere.

## GNU gprof: first pass

```bash
bash benchmarks/profile_gprof.sh /path/to/case './synspec54 < fort.5'
```

The case directory is copied into `results/gprof-<timestamp>/run`. If the input is supplied through another filename, shell redirection, or a wrapper, change the second argument. The command must run the local `./synspec54` executable produced by the script. Output includes `gprof.txt`, `metadata.txt`, compile logs and run logs. The script compiles with `-O2 -g -pg` and preserves the existing static storage and Linux memory-model flags.

**Interpretation:** The gprof flat profile gives approximate self time; the call graph estimates inclusive time and call counts. A call graph may attribute time imperfectly, particularly with recursion, compiler inlining, or low sample counts. Profiling with `-pg` changes runtime; compare relative costs, not its absolute wall time against uninstrumented builds. Run several times, and target runs lasting at least tens of seconds. On systems without `gprof`, Linux `perf record -g` on a `-O2 -g -fno-omit-frame-pointer` build is an alternative (subject to system permissions).

**Nested routines:** `RESOLV` includes `LINOP` and `MOLOP`. `RTECD` and `RTEDFE` may also be nested. Do not add inclusive timings together.

**Next phase:** add optional wall/CPU timers at explicit phase boundaries (atmosphere load, line selection, opacity, atomic and molecular line opacity, transfer, output), and measure initial model setup versus per-abundance variation. Confirm that diagnostic output and counters do not change the synthesized spectrum. This script does not alter the Fortran source.

The GNU gprof script has not yet produced a full-program profile in this macOS environment. Actual macOS sampling and explicit CPU/wall stage-timer results are documented in O2_STAGE_PROFILE.md and BROADER_VALIDATION.md.
