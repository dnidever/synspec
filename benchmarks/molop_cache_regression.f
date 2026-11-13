      PROGRAM TESTCACHE
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
      DOUBLE PRECISION A(MFREQ),E(MFREQ),B(MFREQ),F(MFREQ)
      PARAMETER (MVOI=2001)
      COMMON/VOITAB/H0TAB(MVOI),H1TAB(MVOI),H2TAB(MVOI)
      DO IV=1,MVOI
         H0TAB(IV)=EXP(-((DBLE(IV)-1.D0)/200.D0)**2)
         H1TAB(IV)=0.01D0
         H2TAB(IV)=0.002D0
      END DO
      ND=2
      NFREQ=8
      NFREQS=8
      DFRCON=1.D0
      TMOLIM=8000.D0
      DO J=1,2
         CALL MOLBCLR(J)
         NLINM0(J)=8193
         NLINML(J)=2
         INACTM(J)=0
         DO I=1,2
            IL=I
            IF(I.EQ.2) IL=8193
            INMLIN(I,J)=IL
            IJCMTR(I,J)=5
            INDATM(IL,J)=1
            GFM(IL,J)=-2.D0
            EXCLM(IL,J)=100.D0
            GRM(IL,J)=0.02D0
            GSM(IL,J)=0.D0
            FREQM(IL,J)=5.D0
         END DO
      END DO
      DO ID=1,2
         TEMP(ID)=3500.D0
         ELEC(ID)=1.D0
         DOPMOL(1,ID)=1.D0
         RRMOL(1,ID)=1.D0
         STIM(ID)=0.5D0
         PLAN(ID)=2.D0
      END DO
      DO I=1,NFREQ
         FREQ(I)=DBLE(I)
      END DO
      DO ITEST=1,10
         IF(ITEST.EQ.2) TEMP(1)=4200.D0
         IF(ITEST.EQ.3) THEN
            RRMOL(1,1)=2.D0
            STIM(1)=0.7D0
            DOPMOL(1,1)=0.8D0
         END IF
         IF(ITEST.EQ.4) THEN
            GFM(1,1)=-1.D0
            CALL MOLBCLR(1)
         END IF
         IF(ITEST.EQ.5) THEN
            NLINML(1)=1
            INMLIN(1,1)=8193
         END IF
         IF(ITEST.EQ.7) THEN
            ND=10
            NLINM0(1)=MLINM0
         END IF
         IF(ITEST.EQ.8) THEN
            NLINM0(1)=8193
            ND=3
         END IF
         IF(ITEST.EQ.9) THEN
            GRM(1,1)=0.8D0
            GRM(8193,2)=0.4D0
         END IF
         IF(ITEST.EQ.10) FREQ(NFREQ)=100.D0
         DO J=1,2
            DO ID=1,2
               CALL MOLOP(ID,A,E,0.001D0,J)
               CALL MOLOP_REF(ID,B,F,0.001D0,J)
               DO I=1,NFREQ
                  IF(A(I).NE.B(I).OR.E(I).NE.F(I)) THEN
                     WRITE(*,*) 'MISMATCH',ITEST,J,ID,I
                     STOP 1
                  END IF
               END DO
            END DO
         END DO
      END DO
      WRITE(*,*) 'PASS: profile branches and cache invalidation'
      END
      DOUBLE PRECISION FUNCTION GVDW(IL,ILIST,ID)
      INTEGER IL,ILIST,ID
      GVDW=0.01D0
      END
