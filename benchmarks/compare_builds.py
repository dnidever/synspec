#!/usr/bin/env python3
"""Build and optionally benchmark SYNSPEC 54 with reproducible gfortran flags.

Does not modify src/makefile, src/*.f, or tracked production executables.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "src" / "synspec54.f"
FLAGS = {
    "baseline": [],
    "O2": ["-O2"],
    "O3": ["-O3"],
    "native": ["-O3", "-march=native"],
    "fast-math": ["-O3", "-march=native", "-ffast-math"],
}
INPUT_EXCLUSIONS = {"synspec54", "benchmark-metadata.json"}

def checksum(path):
    digest = hashlib.sha256()
    with path.open("rb") as file:
        for block in iter(lambda: file.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def record_outputs(case_dir):
    outputs = {}
    for file in sorted(case_dir.rglob("*")):
        if file.is_file() and file.name not in INPUT_EXCLUSIONS:
            outputs[str(file.relative_to(case_dir))] = checksum(file)
    return outputs

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default="gfortran")
    parser.add_argument("--output", type=Path, default=Path("results/compiler"))
    parser.add_argument("--input-dir", type=Path, help="Complete reproducible SYNSPEC run directory")
    parser.add_argument("--command", default="./synspec54", help="Shell command to run within each copied input directory")
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--build-only", action="store_true")
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument("--variants", nargs="+", choices=list(FLAGS), default=list(FLAGS))
    args = parser.parse_args()

    if args.repeats < 1:
        parser.error("--repeats must be positive")
    if not args.build_only and args.input_dir is None:
        parser.error("--input-dir is required unless --build-only is set")
    out = args.output.resolve()
    inputs = args.input_dir.resolve() if args.input_dir else None
    if inputs and not inputs.is_dir():
        parser.error("--input-dir must exist")
    if inputs and (out == inputs or inputs in out.parents):
        parser.error("--output must not be inside --input-dir")
    if inputs and (inputs == ROOT or ROOT in inputs.parents):
        parser.error("--input-dir must not contain the repository")
    out.mkdir(parents=True, exist_ok=True)

    version = subprocess.run([args.compiler, "--version"], capture_output=True, text=True)
    metadata = {
        "source": str(SOURCE.relative_to(ROOT)),
        "source_sha256": checksum(SOURCE),
        "compiler": args.compiler,
        "compiler_version": version.stdout.splitlines()[0] if version.stdout else "",
        "machine": platform.platform() + " / " + platform.machine(),
        "variants": FLAGS,
        "input_dir": str(inputs) if inputs else None,
        "command": args.command,
        "repeats": args.repeats,
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")
    rows = []
    buildroot = out / "build"
    buildroot.mkdir(exist_ok=True)
    for variant in args.variants:
        build_dir = buildroot / variant
        build_dir.mkdir(exist_ok=True)
        exe = build_dir / "synspec54"
        flags = ["-fno-automatic"]
        if sys.platform.startswith("linux"):
            flags.append("-mcmodel=medium")
        flags += FLAGS[variant]
        cmd = [args.compiler, *flags, str(SOURCE), "-o", str(exe)]
        print("BUILD", variant, shlex.join(cmd), flush=True)
        built = subprocess.run(cmd, cwd=SOURCE.parent, capture_output=True, text=True)
        (build_dir / "compile.stdout").write_text(built.stdout)
        (build_dir / "compile.stderr").write_text(built.stderr)
        if built.returncode != 0:
            print("  compile failed; see compile.stderr", file=sys.stderr)
            rows.append({"variant": variant, "repeat": "", "build_returncode": built.returncode,
                         "run_returncode": "", "seconds": "", "output_hashes": ""})
            continue
        if args.build_only:
            rows.append({"variant": variant, "repeat": "", "build_returncode": 0,
                         "run_returncode": "", "seconds": "", "output_hashes": ""})
            continue
        for run_id in range(1, args.repeats + 1):
            case = out / "runs" / variant / str(run_id)
            if case.exists():
                shutil.rmtree(case)
            shutil.copytree(inputs, case)
            shutil.copy2(exe, case / "synspec54")
            before = record_outputs(case)
            env = os.environ.copy()
            env.setdefault("OMP_NUM_THREADS", "1")
            start = time.perf_counter()
            try:
                result = subprocess.run(args.command, cwd=case, shell=True, text=True,
                                        capture_output=True, timeout=args.timeout, env=env)
                code = result.returncode
                stdout, stderr = result.stdout, result.stderr
            except subprocess.TimeoutExpired as exc:
                code = 124
                stdout = (exc.stdout or b"").decode(errors="replace") if isinstance(exc.stdout, bytes) else (exc.stdout or "")
                stderr = (exc.stderr or b"").decode(errors="replace") if isinstance(exc.stderr, bytes) else (exc.stderr or "")
            elapsed = time.perf_counter() - start
            (case / "run.stdout").write_text(stdout)
            (case / "run.stderr").write_text(stderr)
            after = record_outputs(case)
            generated = {k: v for k, v in after.items() if before.get(k) != v and k not in ("run.stdout", "run.stderr")}
            (case / "output_hashes.json").write_text(json.dumps(generated, indent=2) + "\n")
            rows.append({"variant": variant, "repeat": run_id, "build_returncode": 0,
                         "run_returncode": code, "seconds": f"{elapsed:.6f}",
                         "output_hashes": str(case / "output_hashes.json")})
            print(f"  run {run_id}: exit={code}, wall={elapsed:.3f} s", flush=True)
    with (out / "results.csv").open("w", newline="") as file:
        writer = csv.DictWriter(file, fieldnames=["variant", "repeat", "build_returncode",
                                                  "run_returncode", "seconds", "output_hashes"])
        writer.writeheader()
        writer.writerows(rows)
    print("Results:", out / "results.csv")

if __name__ == "__main__":
    main()
