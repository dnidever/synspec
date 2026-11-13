#!/usr/bin/env python3
"""Reject nonfinite SYNSPEC output and compare spectra on identical grids."""
import argparse
import json
from pathlib import Path
import numpy as np

def check(reference, candidate, rtol=1e-5, atol=0.0):
    ref = np.loadtxt(reference)
    got = np.loadtxt(candidate)
    if ref.ndim != 2 or ref.shape[1] != 2 or got.shape != ref.shape:
        raise ValueError('Expected matching two-column spectra')
    if not np.isfinite(ref).all() or not np.isfinite(got).all():
        raise ValueError('Spectrum contains NaN or infinity')
    if not np.array_equal(ref[:, 0], got[:, 0]):
        raise ValueError('Wavelength grids differ; align them before comparison')
    if not np.all(np.diff(got[:, 0]) > 0):
        raise ValueError('Wavelength grid is not strictly increasing')
    difference = np.abs(ref[:, 1] - got[:, 1])
    passed = bool(np.allclose(ref[:, 1], got[:, 1], rtol=rtol, atol=atol))
    return dict(passed=passed, samples=len(got), max_absolute_difference=float(difference.max()),
                rtol=rtol, atol=atol)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('reference', type=Path)
    parser.add_argument('candidate', type=Path)
    parser.add_argument('--rtol', type=float, default=1e-5)
    parser.add_argument('--atol', type=float, default=0.0)
    args = parser.parse_args()
    if args.rtol < 0 or args.atol < 0:
        parser.error('Tolerances must be nonnegative')
    try:
        result = check(args.reference, args.candidate, args.rtol, args.atol)
    except ValueError as exc:
        parser.exit(1, str(exc) + '\n')
    print(json.dumps(result, indent=2))
    raise SystemExit(0 if result['passed'] else 1)
