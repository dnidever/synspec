#!/usr/bin/env python3
"""Prepare APOGEE/H-band SYNSPEC binary line lists from Carlos's distribution.

The files are large; the downloadable archives are streamed to disk and never
stored in Git. This helper invokes the existing list2bin executable.
"""
import argparse
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
import urllib.request

BASE = "ftp://carlos:allende@ftp.ll.iac.es/linelists/"
SETS = {
    "hband-giant": ["synspec.20200921nlte.19",
                    "synspec.20200921nltegiant.20",
                    "synspec.h2o-BC8.5V.20"],
    "hband-solar": ["synspec.20200921nlte.19",
                    "synspec.20200921nltesun.20",
                    "synspec.h2o-BC8.5V.20"],
    "standard-giant": ["gfATO.19", "gfMOLgiant.20",
                       "gfTiO.20", "H2O-8.20"],
    "standard-solar": ["gfATO.19", "gfMOLsun.20",
                       "gfTiO.20", "H2O-8.20"],
}

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def prepare_one(name, target, converter, base, keep_gz=False):
    binary = target / (name + ".11")
    if binary.is_file() and binary.stat().st_size:
        print("Already present:", binary, flush=True)
        return {"name": name, "binary": str(binary),
                "bytes": binary.stat().st_size, "sha256": digest(binary),
                "skipped_existing": True}
    archive = target / (name + ".gz")
    source = base.rstrip("/") + "/" + name + ".gz"
    if not archive.exists():
        print("Downloading", source, flush=True)
        with urllib.request.urlopen(source, timeout=120) as response, archive.open("wb") as fout:
            shutil.copyfileobj(response, fout, 1024 * 1024)
    archive_hash = digest(archive)
    # list2bin emits fort.12 in its current directory; isolate the invocation.
    with tempfile.TemporaryDirectory(prefix="synspec-linelist-") as temp:
        tempdir = Path(temp)
        text_path = tempdir / name
        with gzip.open(archive, "rb") as compressed, text_path.open("wb") as outfile:
            shutil.copyfileobj(compressed, outfile, 1024 * 1024)
        # Legacy .00 file records number of records + first/last line.
        with text_path.open("rb") as stream:
            count = 0
            first = last = b""
            for line in stream:
                count += 1
                if count == 1:
                    first = line
                last = line
        (target / (name + ".00")).write_bytes(
            (str(count) + " " + name + "\n").encode() + first + last
        )
        with text_path.open("rb") as instream, (tempdir / "convert.log").open("wb") as log:
            proc = subprocess.run([str(converter)], cwd=tempdir,
                                  stdin=instream, stdout=log, stderr=subprocess.STDOUT)
        if proc.returncode != 0 or not (tempdir / "fort.12").is_file():
            shutil.copy2(tempdir / "convert.log", target / (name + ".convert-failed.log"))
            raise RuntimeError(f"list2bin conversion failed for {name}: {proc.returncode}")
        shutil.move(str(tempdir / "fort.12"), binary)
    if not keep_gz:
        archive.unlink()
    return {"name": name, "binary": str(binary), "bytes": binary.stat().st_size,
            "sha256": digest(binary), "archive_sha256": archive_hash,
            "source": source}

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--set", choices=SETS, default="hband-giant")
    p.add_argument("--output", type=Path, default=Path("benchmark-data/linelists"))
    p.add_argument("--list2bin", type=Path, default=Path("bin/list2bin"))
    p.add_argument("--base-url", default=BASE, help="Original FTP mirror or user-provided HTTPS mirror")
    p.add_argument("--keep-archives", action="store_true")
    args = p.parse_args()
    target = args.output.resolve()
    converter = args.list2bin.resolve()
    if not converter.is_file():
        p.error(f"Converter missing: {converter}. Build using 'make -C src list2bin'.")
    target.mkdir(parents=True, exist_ok=True)
    manifest = {"selection": args.set, "files": []}
    for name in SETS[args.set]:
        row = prepare_one(name, target, converter, args.base_url, args.keep_archives)
        manifest["files"].append(row)
        (target / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print("Prepared", len(manifest["files"]), "lists at", target)

if __name__ == "__main__":
    main()
