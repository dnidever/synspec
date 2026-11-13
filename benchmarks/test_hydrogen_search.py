#!/usr/bin/env python3
"""Check the production hydrogen lookup against the old linear search."""
import argparse, pathlib, tempfile, subprocess, sys
p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=pathlib.Path,default=pathlib.Path(__file__).resolve().parents[1]/'src/synspec54.f');p.add_argument('--compiler',default='gfortran');a=p.parse_args()
s=a.source.read_text();start=s.index('C              First tabulated upper endpoint >= AL; clamp tails.');end=s.index('               IW1=IW0+1',start)+len('               IW1=IW0+1');body=s[start:end]
program='''      PROGRAM TESTSEARCH
      IMPLICIT NONE
      DOUBLE PRECISION WLHYD(1,55),AL,EPS
      INTEGER N,I,J,K,EXPECTED,GOT,COUNT
      COUNT=0
      EPS=1.D-10
      DO N=2,55
         DO I=1,N
            WLHYD(1,I)=DBLE(I*I)/7.D0
         END DO
         DO J=-10,210
            AL=DBLE(J)*WLHYD(1,N)/200.D0
            CALL CHECKONE(AL,WLHYD,N,COUNT)
         END DO
         DO I=1,N
            DO J=-1,1
               AL=WLHYD(1,I)+DBLE(J)*EPS
               CALL CHECKONE(AL,WLHYD,N,COUNT)
            END DO
         END DO
      END DO
      WRITE(*,*) 'PASS: hydrogen interval cases',COUNT
      END
      SUBROUTINE CHECKONE(AL,WLHYD,N,COUNT)
      IMPLICIT NONE
      DOUBLE PRECISION AL,WLHYD(1,55)
      INTEGER N,COUNT,K,EXPECTED,GOT
      EXPECTED=N-1
      DO K=1,N-1
         IF(AL.LE.WLHYD(1,K+1)) THEN
            EXPECTED=K
            EXIT
         END IF
      END DO
      CALL LOOKUP(AL,WLHYD,N,GOT)
      IF(GOT.NE.EXPECTED) STOP 'LOOKUP MISMATCH'
      COUNT=COUNT+1
      END
      SUBROUTINE LOOKUP(AL,WLHYD,NWL,IW0)
      IMPLICIT NONE
      DOUBLE PRECISION AL,WLHYD(1,55)
      INTEGER NWL,IW0,IW1,ILO,IHI,IMID,ILINE
      ILINE=1
'''+body+'''
      END
'''
with tempfile.TemporaryDirectory(prefix='hydrogen-lookup-') as d:
 d=pathlib.Path(d);(d/'test.f').write_text(program);exe=d/'test';subprocess.run([a.compiler,'-O2','-fcheck=all',str(d/'test.f'),'-o',str(exe)],check=True);subprocess.run([str(exe)],check=True)
