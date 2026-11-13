#!/usr/bin/env python3
"""Prepare an external LTE/plane-parallel full-band APOGEE benchmark deck.

Requires the repository Python dependencies and already converted DR17 lists.
Does not generate/interpolate atmospheres or certify line-list correctness.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--model', type=Path, required=True)
    p.add_argument('--linelist-dir', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--grid-config', type=Path,
                   default=ROOT/'benchmarks/apogee_grid.json')
    p.add_argument('--sampling', type=float,
                   help='Override benchmark sampling control in Angstroms')
    args = p.parse_args()
    config = json.loads(args.grid_config.read_text())
    sampling = (args.sampling if args.sampling is not None else
                config['benchmark_sampling_control_angstrom'])
    if sampling <= 0:
        p.error('Sampling control must be positive')
    model = args.model.resolve()
    directory = args.linelist_dir.resolve()
    files = [directory/n for n in ['synspec.20200921nlte.19.11',
                                  'synspec.20200921nltegiant.20.11',
                                  'synspec.h2o-BC8.5V.20.11']]
    for file in [model, *files]:
        if not file.is_file():
            p.error('Required input missing: '+str(file))
    output = args.output.resolve()
    if output.exists():
        p.error('Use a new output directory')
    if output == ROOT or ROOT in output.parents:
        p.error('Keep prepared decks and large input links outside the checkout')
    sys.path.insert(0, str(ROOT/'python'))
    from synspec import models, synthesis
    kind, teff, logg, vmicro, abundances, nd, atmosphere = models.read_model(str(model))
    if not config['teff_bounds_k'][0] <= teff <= config['teff_bounds_k'][1]:
        p.error('Atmosphere temperature is outside the target benchmark bounds')
    if not config['logg_bounds_cgs'][0] <= logg <= config['logg_bounds_cgs'][1]:
        p.error('Atmosphere gravity is outside the target benchmark bounds')
    header = model.read_text().splitlines()
    composition = {}
    if kind == 'marcs':
        fields = header[6].split()
        composition = dict(model_fe_h_header=float(fields[0]),
                           model_alpha_fe_header=float(fields[1]))
        if not config['fe_h_bounds_dex'][0] <= composition['model_fe_h_header'] <= config['fe_h_bounds_dex'][1]:
            p.error('MARCS header metallicity is outside the target bounds')
    metadata = dict(model=str(model), model_sha256=digest(model),
                    model_identifier_header=header[0],
                    filled_model_input=model.name.endswith('.filled'),
                    **composition,
                    teff=float(teff), logg=float(logg), vmicro=float(vmicro),
                    depths=int(nd), atmosphere_type=kind,
                    iron_abundance_nh=float(abundances[25]),
                    wavelength_angstrom=config['wavelength_angstrom'],
                    sampling_control_angstrom=sampling,
                    target_grid=config,
                    linelists=[dict(path=str(f), bytes=f.stat().st_size,
                                    sha256=digest(f)) for f in files],
                    geometry='LTE plane-parallel transfer; structure as supplied')
    output.mkdir(parents=True)
    cwd = Path.cwd()
    try:
        os.chdir(output)
        (output/'data').symlink_to(ROOT/'python/synspec/data')
        synthesis.write5(teff, logg, abundances)
        synthesis.write8(teff, logg, nd, atmosphere, kind)
        synthesis.write55(config['wavelength_angstrom'], dw=sampling,
                          imode=0, cutoff0=250., strength=1e-4,
                          vmicro=vmicro, linelist=[str(f) for f in files],
                          atmostype=kind)
        synthesis.writetas('tas', nd, [str(f) for f in files])
        synthesis.create_links([str(f) for f in files])
        (output/'benchmark-metadata.json').write_text(json.dumps(metadata, indent=2)+'\n')
    finally:
        os.chdir(cwd)
    print(output)


if __name__ == '__main__':
    main()
