#!/usr/bin/env python3
"""Exercise real MOLOP/VOIGTK with synthetic tables and bounds/FP checks."""
import argparse
from pathlib import Path
import subprocess, tempfile, os, sys

def routine(source, name):
    start=source.index('      SUBROUTINE '+name+'(')
    end=source.index('\n      END\n',start)+len('\n      END\n')
    return source[start:end]

ROOT=Path(__file__).resolve().parents[1]

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--compiler',default='gfortran')
    p.add_argument('--root',type=Path,default=ROOT)
    p.add_argument('--source',type=Path)
    p.add_argument('--reference','--reference-source',dest='reference',type=Path)
    p.add_argument('--driver',type=Path)
    a=p.parse_args()
    a.source=a.source or a.root/'src/synspec54.f'
    a.driver=a.driver or a.root/'benchmarks/molop_cache_regression.f'
    cur=a.source.read_text()
    ref=(a.reference or a.root/'benchmarks/reference_molop_voigt.f').read_text()
    start=ref.index('      FUNCTION VOIGTK(')
    end=ref.index('\n      END\n',start)+len('\n      END\n')
    with tempfile.TemporaryDirectory(prefix='profile-regression-') as t:
        t=Path(t)
        (t/'reference.f').write_text(routine(ref,'MOLOP').replace('SUBROUTINE MOLOP(', 'SUBROUTINE MOLOP_REF(',1)+'\n'+ref[start:end])
        (t/'candidate.f').write_text(routine(cur,'MOLOP')+'\n'+cur[cur.index('      SUBROUTINE MOLBCLR('):])
        flags=['-fno-automatic','-O2','-fcheck=all','-ffpe-trap=invalid,zero,overflow']
        if sys.platform.startswith('linux'):flags.append('-mcmodel=medium')
        exe=t/'check'
        subprocess.run([a.compiler,*flags,'-I',str(a.root/'src'),str(a.driver),str(t/'reference.f'),str(t/'candidate.f'),'-o',str(exe)],check=True)
        for enabled in [False,True]:
            env=os.environ.copy();env.pop('SYNSPEC_MOL_CACHE',None)
            if enabled:env['SYNSPEC_MOL_CACHE']='1'
            print('Cache enabled:',enabled,flush=True)
            subprocess.run([str(exe)],env=env,check=True)
if __name__=='__main__':main()
