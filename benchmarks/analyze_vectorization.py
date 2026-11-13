#!/usr/bin/env python3
"""Validate compiler-comparison spectra and summarize vectorization reports.

Requires NumPy. Run after compare_vectorization.py completes.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import statistics

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
ROUTINES = ['MOLOP', 'HYDLIN', 'OPAC', 'OPACON', 'RTECD', 'INILIN', 'INMOLI']


def run_dir(results, row):
    return results/'runs'/row['case']/str(row['repeat'])/row['variant']


def compare(results, row, reference):
    directory = run_dir(results, row)
    refdir = run_dir(results, reference)
    got = [np.loadtxt(directory/f, ndmin=2) for f in ['fort.7', 'fort.17']]
    ref = [np.loadtxt(refdir/f, ndmin=2) for f in ['fort.7', 'fort.17']]
    finite = all(bool(np.isfinite(a).all()) for a in got+ref)
    shape_ok = all(a.shape == b.shape and a.shape[1] == 2 and len(a) > 0
                   for a, b in zip(got, ref))
    same_grid = shape_ok and all(np.array_equal(a[:, 0], b[:, 0])
                                 for a, b in zip(got, ref))
    check = {k: row[k] for k in ['case', 'repeat', 'warmup', 'variant']}
    check.update(finite=finite, same_grid=bool(same_grid),
                 exact_printed_arrays=all(np.array_equal(a, b) for a, b in zip(got, ref)))
    for name, filename in [('flux', 'fort.7'), ('continuum', 'fort.17')]:
        check[name+'_sha256'] = hashlib.sha256((directory/filename).read_bytes()).hexdigest()
    if finite and same_grid:
        a, c = got
        b, e = ref
        norm_a = a[:, 1]/np.interp(a[:, 0], c[:, 0], c[:, 1])
        norm_b = b[:, 1]/np.interp(b[:, 0], e[:, 0], e[:, 1])
        check.update(
            max_relative_flux_difference=float(np.max(abs(a[:, 1]-b[:, 1])/np.maximum(abs(b[:, 1]), 1e-99))),
            max_relative_continuum_difference=float(np.max(abs(c[:, 1]-e[:, 1])/np.maximum(abs(e[:, 1]), 1e-99))),
            max_normalized_flux_difference=float(np.max(abs(norm_a-norm_b))))
    check['passes_numerical_screen'] = bool(finite and same_grid and all(
        check[k] <= 1e-5 for k in ['max_relative_flux_difference',
                                   'max_relative_continuum_difference',
                                   'max_normalized_flux_difference']))
    return check


def summarize(rows, checks):
    summary = []
    for case in dict.fromkeys(x['case'] for x in rows):
        refs = [x for x in rows if x['case'] == case and x['variant'] == 'O2' and not x['warmup']]
        base = statistics.median(x['seconds'] for x in refs)
        for variant in ['O2', 'O3', 'native', 'O3-no-vector']:
            group = [x for x in rows if x['case'] == case and x['variant'] == variant and not x['warmup']]
            if not group:
                continue
            med = statistics.median(x['seconds'] for x in group)
            paired, paired_cpu = [], []
            for row in group:
                ref = next(a for a in refs if a['repeat'] == row['repeat'])
                paired.append(ref['seconds']/row['seconds'])
                paired_cpu.append(ref['cpu_seconds']/row['cpu_seconds'])
            selected = [x for x in checks if x['case'] == case and x['variant'] == variant and not x['warmup']]
            summary.append(dict(
                case=case, variant=variant, repeats=len(group), median_seconds=med,
                min_seconds=min(x['seconds'] for x in group),
                max_seconds=max(x['seconds'] for x in group),
                speedup_vs_O2=base/med, time_saved_percent=100*(1-med/base),
                median_paired_speedup=statistics.median(paired),
                median_paired_cpu_speedup=statistics.median(paired_cpu),
                all_finite_exact=all(x['finite'] and x['exact_printed_arrays'] for x in selected),
                all_pass_numerical_screen=all(x['passes_numerical_screen'] for x in selected)))
    return summary


def compiler_analysis(results, source):
    lines = source.read_text().splitlines()
    regions = {}
    for name in ROUTINES:
        start = next(i+1 for i, line in enumerate(lines)
                     if re.match(r'^\s*subroutine\s+'+name+r'\b', line, re.I))
        end = next(i+1 for i in range(start, len(lines))
                   if re.match(r'^\s*end\s*$', lines[i], re.I))
        regions[name] = (start, end)
    variants = {}
    for variant in ['O2', 'O3', 'native', 'O3-no-vector']:
        successes, misses = [], {}
        report = results/'build'/variant/'vectorization.txt'
        if not report.exists():
            continue
        for line in report.read_text(errors='replace').splitlines():
            match = re.match(r'(.+?):(\d+):(\d+): (optimized|missed): (.*)', line)
            if not match:
                continue
            filename, number, _, kind, message = match.groups()
            if not filename.endswith('synspec54.f'):
                continue
            number = int(number)
            routine = next((n for n, (a, b) in regions.items() if a <= number <= b), None)
            if kind == 'optimized':
                successes.append(dict(line=number, routine=routine, message=message))
            elif routine and ('not vectorized:' in message or "couldn't vectorize loop" in message):
                items = misses.setdefault(routine, [])
                if len(items) < 25:
                    items.append(dict(line=number, message=message))
        variants[variant] = dict(optimized=successes, missed_hot_loops=misses)
    return dict(regions=regions, variants=variants)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--results', type=Path, required=True)
    p.add_argument('--output', type=Path)
    p.add_argument('--source', type=Path, default=ROOT/'src/synspec54.f')
    args = p.parse_args()
    results = args.results.resolve()
    output = (args.output or results).resolve()
    output.mkdir(parents=True, exist_ok=True)
    metadata = json.loads((results/'metadata.json').read_text())
    if hashlib.sha256(args.source.read_bytes()).hexdigest() != metadata['source_sha256']:
        raise RuntimeError('Source does not match the measured build')
    rows = json.loads((results/'runs.json').read_text())
    checks = []
    for row in rows:
        if row['returncode']:
            raise RuntimeError('A synthesis failed; inspect raw results before analysis')
        reference = next(x for x in rows if x['case'] == row['case'] and
                         x['repeat'] == row['repeat'] and x['variant'] == 'O2')
        checks.append(compare(results, row, reference))
    data = {'spectrum-checks.json': checks, 'summary.json': summarize(rows, checks),
            'compiler-analysis.json': compiler_analysis(results, args.source)}
    for name, value in data.items():
        (output/name).write_text(json.dumps(value, indent=2)+'\n')
    print('Checked', len(checks), 'executions;',
          sum(x['finite'] and x['exact_printed_arrays'] for x in checks),
          'finite and exactly equal at printed precision')
    if not all(x['passes_numerical_screen'] for x in checks):
        raise SystemExit('Numerical screening failed; inspect spectrum-checks.json')


if __name__ == '__main__':
    main()
