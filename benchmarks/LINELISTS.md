# Line-list preparation for SYNSPEC benchmarks

The existing upstream `linelists/makefile` and this project's `python/synspec/linelists/makefile` use an authenticated FTP archive hosted at `ftp.ll.iac.es` and convert gzipped ASCII files into `.11` binary format via the project's `list2bin` executable. The files can be very large and should **not** be committed.

For HGE, prioritize the APOGEE-specific atomic list with the red-giant molecular list:

- `synspec.20200921nlte.19.gz`
- `synspec.20200921nltegiant.20.gz`
- `synspec.h2o-BC8.5V.20.gz`

The solar molecular list is an alternative, not an extra list to stack alongside the giant list for the same synthesis.

## Prepare

```bash
make -C src list2bin
python3 benchmarks/prepare_linelists.py --set hband-giant
```

Alternative presets: `--set hband-solar`, `--set standard-giant`, `--set standard-solar`. Files go into `benchmark-data/linelists/` by default. The helper preserves source checksums and output SHA-256 values in `manifest.json`, creates the legacy `.00` file, and avoids loading compressed/decompressed contents entirely into RAM. Use `--base-url` for an authorized alternate mirror, or transfer the archives manually into the output directory first. The source FTP service may be unreachable from some environments; in that case the script fails, and does not invent a successful download.

```bash
python3 benchmarks/prepare_linelists.py --set hband-giant \
    --base-url 'https://your-authorized-mirror.example/linelists'
```

After preparation, generate a real synthesis input deck via `python/synspec/synthesis.py`'s `do_synspec` with `save=True` or a persistent workdir, and profile it with `benchmarks/profile_gprof.sh`. The `synthesize` wrapper needs a valid atmosphere file, line-list paths, and may depend on separately installed atmospheric data and atomic-opacity files. The existing test atmosphere `python/synspec/test/s6000_g+1.0_m0.5.mod` offers a small initial functional test but is not representative of a cool red-giant H-band benchmark.

Do not publish the source FTP credentials in generated logs. The upstream Makefiles contain a historical FTP login already available publicly; this script stores only the source URL in the manifest, so review/redact any private alternative mirror URLs before sharing.

Do not use downloaded files as a scientifically authoritative HGE line list without validating data versions, line-list formatting and license/access permissions.


## Official DR17 source used in local validation

The public [SDSS DR17 SYNSPEC directory](https://dr17.sdss.org/sas/dr17/apogee/spectro/speclib/linelists/synspec/) supplies `synspec.20200921nlte.atoms`, `synspec.20200921nltegiant.molec`, and `synspec.h2o-BC8.5V.molec`. These are plain text, not the gzipped filenames expected by the FTP preparation helper. Local downloads were verified against SDSS's SHA-1 manifest and converted separately; no large files are in Git.

The atomic records are 103 characters, while the original converter reads only 102 characters of the first record and therefore detects the wrong format. The experimental `benchmarks/list2bin_dr17.f` expands that buffer. The atomic source also has a nonnumeric `s` damping field on row 33407, which the existing conversion error handler skips; the provisional binary therefore has 51,260 rather than 51,261 records. This data issue is unresolved and the benchmark must not be mistaken for an authoritative production line-list preparation. General-purpose lists were used for the broader optical/H-band matrix. See DR17_LOCAL_RESULTS.md and BROADER_VALIDATION.md.
