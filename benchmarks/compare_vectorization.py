#!/usr/bin/env python3
"""Compare strict-math compiler variants on prepared SYNSPEC input decks.

Manifest: JSON array of objects with 'name' and 'case' (directory) keys.
Large input files are preserved as symlinks; no line data are downloaded.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import resource
import shutil
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
FLAGS = {'O2': ['-O2'], 'O3': ['-O3'],
         'native': ['-O3', '-march=native'],
         'O3-no-vector': ['-O3', '-fno-tree-vectorize']}

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--compiler', default='gfortran')
    p.add_argument('--case-manifest', type=Path, required=True)
    p.add_argument('--wide-case', type=Path)
    p.add_argument('--runtime-library-dir', type=Path,
                   help='Optional macOS directory for this compiler runtime')
    p.add_argument('--reuse-builds', action='store_true',
                   help='Reuse matching binaries in the output directory')
    p.add_argument('--wide-repeats', type=int, default=5)
    p.add_argument('--matrix-repeats', type=int, default=2)
    p.add_argument('--timeout', type=int, default=300)
    p.add_argument('--output', type=Path, default=ROOT/'results/vectorization')
    args = p.parse_args()
    if min(args.wide_repeats, args.matrix_repeats, args.timeout) < 1:
        p.error('Repeats and timeout must be positive')
    out = args.output.resolve()
    # New output directories prevent accidental replacement of prior runs.
    source = ROOT/'src/synspec54.f'
    cases = json.loads(args.case_manifest.read_text())
    schedule = []
    if args.wide_case:
        schedule.append(('DR17-wide', args.wide_case.resolve(),
                         args.wide_repeats, list(FLAGS)))
    for c in cases:
        name = c['name']
        if not name or name in ('.', '..') or '/' in name or '\\' in name:
            p.error('Case names must be single directory names')
        schedule.append((name, Path(c['case']).resolve(),
                         args.matrix_repeats, ['O2', 'O3', 'native']))
    for _, case, _, _ in schedule:
        if not case.is_dir() or case == out or case in out.parents or out in case.parents:
            p.error('Each input case must exist and be separate from output')
    previous = None
    if args.reuse_builds:
        previous = json.loads((out/'metadata.json').read_text())
    out.mkdir(parents=True, exist_ok=args.reuse_builds)
    metadata = dict(source_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),
                    compiler=subprocess.check_output([args.compiler, '--version'], text=True).splitlines()[0],
                    platform=platform.platform(), flags=FLAGS,
                    common_flags=['-fno-automatic'], seed=5137,
                    runtime_library_dir=str(args.runtime_library_dir.resolve()) if args.runtime_library_dir else None,
                    warmups=1, schedule=[dict(name=n,case=str(c),repeats=r,variants=v)
                                         for n,c,r,v in schedule])
    if previous and any(previous.get(k) != metadata[k] for k in
                        ['source_sha256', 'compiler', 'platform', 'flags', 'common_flags']):
        raise RuntimeError('Existing build metadata does not match current source/compiler')
    (out/'metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    exe = {}
    for name, flags in FLAGS.items():
        b = out/'build'/name
        b.mkdir(parents=True, exist_ok=args.reuse_builds)
        exe[name] = b/'synspec54'
        if args.reuse_builds:
            if not exe[name].is_file():
                raise RuntimeError('Missing prebuilt executable: '+str(exe[name]))
            print('REUSING', name, flush=True)
            continue
        common = ['-fno-automatic']
        if sys.platform.startswith('linux'):
            common.append('-mcmodel=medium')
        command = [args.compiler, *common, *flags,
                   '-fopt-info-vec-all='+str(b/'vectorization.txt'),
                   str(source), '-o', str(exe[name])]
        (b/'command.json').write_text(json.dumps(command, indent=2)+'\n')
        with (b/'compile.log').open('wb') as log:
            result = subprocess.run(command, cwd=source.parent, stdout=log, stderr=subprocess.STDOUT)
        if result.returncode:
            raise RuntimeError('Compilation failed: '+name)
        print('BUILT', name, flush=True)
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', SYNSPEC_MOL_CACHE='0')
    if args.runtime_library_dir:
        if sys.platform != 'darwin':
            raise RuntimeError('--runtime-library-dir currently targets macOS')
        env['DYLD_LIBRARY_PATH'] = str(args.runtime_library_dir.resolve())
    rows = []
    rng = random.Random(5137)
    for name, case, repeats, variants in schedule:
        for repeat in range(repeats+1):
            order = variants.copy()
            rng.shuffle(order)
            for variant in order:
                d = out/'runs'/name/str(repeat)/variant
                shutil.copytree(case, d, symlinks=True)
                before = resource.getrusage(resource.RUSAGE_CHILDREN)
                with (d/'fort.5').open('rb') as inp, (d/'stdout').open('wb') as log:
                    t = time.perf_counter()
                    try:
                        result = subprocess.run([str(exe[variant])], cwd=d, env=env,
                                                stdin=inp, stdout=log, stderr=subprocess.STDOUT,
                                                timeout=args.timeout)
                        code = result.returncode
                    except subprocess.TimeoutExpired:
                        code = 124
                    elapsed = time.perf_counter()-t
                after = resource.getrusage(resource.RUSAGE_CHILDREN)
                row = dict(case=name, repeat=repeat, warmup=repeat==0, variant=variant,
                           seconds=elapsed, cpu_seconds=after.ru_utime+after.ru_stime-before.ru_utime-before.ru_stime,
                           returncode=code, run_dir=str(d))
                rows.append(row)
                (out/'runs.json').write_text(json.dumps(rows, indent=2)+'\n')
                print(name, repeat, variant, round(elapsed,3), 'exit', code, flush=True)
                if code:
                    raise RuntimeError('Synthesis failed; see '+str(d/'stdout'))
    print('COMPLETE', flush=True)

if __name__ == '__main__':
    main()
