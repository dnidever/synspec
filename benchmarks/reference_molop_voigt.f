C     Immutable pre-optimization routines from commit dceaedda.
C     Kept as a numerical regression oracle, not linked in production.
      SUBROUTINE MOLOP(ID,ABLIN,EMLIN,AVAB,ILIST)
C     ===========================================
C
C     Total molecular line opacity (ABLIN) and emissivity (EMLIN)
C
      INCLUDE 'PARAMS.FOR'
      INCLUDE 'MODELP.FOR'
      INCLUDE 'SYNTHP.FOR'
      INCLUDE 'LINDAT.FOR'
      PARAMETER (UN     = 1.,
     *           EXT0   = 3.17,
     *           TEN    = 10.)
      DIMENSION ABLIN(MFREQ),EMLIN(MFREQ)
C
      DO IJ=1,NFREQ
         ABLIN(IJ)=0.
         EMLIN(IJ)=0.
      END DO
C
      if(temp(id).gt.tmolim) return
      IF(NLINML(ILIST).EQ.0) RETURN
      if(inactm(ilist).ne.0) return
C
C     overall loop over contributing lines
C
      TEM1=UN/TEMP(ID)
      ANE=ELEC(ID)
      DO I=1,NLINML(ILIST)
         IL=INMLIN(I,ILIST)
         IMOL=INDATM(IL,ILIST)
         DOP1=DOPMOL(IMOL,ID)
         AGAM=(GRM(IL,ILIST)+GSM(IL,ILIST)*ANE+
     *         GVDW(IL,ILIST,ID))*DOP1
         FR0=FREQM(IL,ILIST)
         AB0=EXP(GFM(IL,ILIST)-EXCLM(IL,ILIST)*TEM1)*RRMOL(IMOL,ID)*
     *           DOP1*STIM(ID)
C
C        set up limiting frequencies where the line I is supposed to
C        contribute to the opacity
C
         EX0=AB0/AVAB*AGAM
         EXT=EXT0
         IF(EX0.GT.TEN) EXT=SQRT(EX0)
         EXT=EXT/DOP1
         XIJEXT=DFRCON*EXT+1.5
         IJ1=int(MAX(float(IJCMTR(I,ILIST))-XIJEXT,3.))
         IJ2=int(MIN(float(IJCMTR(I,ILIST))+XIJEXT,float(NFREQS)))
         IF(IJ1.LT.NFREQ.AND.IJ2.GT.2) THEN
            DO IJ=IJ1,IJ2
               XF=ABS(FREQ(IJ)-FR0)*DOP1
               ABLIN(IJ)=ABLIN(IJ)+AB0*VOIGTK(AGAM,XF)
            END DO
         END IF
      END DO
C
      DO IJ=3,NFREQ
         EMLIN(IJ)=EMLIN(IJ)+ABLIN(IJ)*PLAN(ID)
      END DO
C
      RETURN
      END

      FUNCTION VOIGTK(A,V)
C     ====================
C
C     Voigt function after Kurucz (in Computational Astrophysics)
C
      INCLUDE 'PARAMS.FOR'
      PARAMETER (MVOI=2001)
      PARAMETER (ONE=1., THREE=3., TEN=10., FIFTN=15., TWOH=200.,
     *           C14142=1.4142, C11283=1.12838, C15=1.5,C32=3.2,
     *           C05642=0.5642,C79788=0.79788,C02=0.2,C14=1.4,
     *           C37613=0.37613,C23=2./3.,
     *           CV1=-.122727278,CV2=.532770573,CV3=-.96284325,
     *           CV4=.979895032)
      COMMON/VOITAB/H0TAB(MVOI),H1TAB(MVOI),H2TAB(MVOI)
      IV=int(V*TWOH+C15)
      IF(A.LT.C02) THEN
         IF(V.LE.TEN) THEN
            VOIGTK=(H2TAB(IV)*A+H1TAB(IV))*A+H0TAB(IV)
          ELSE
            VOIGTK=C05642*A/(V*V)
         END IF
         RETURN
      END IF
c     IF(A.GT.C14) GO TO 10
c     IF(A+V.GT.C32) GO TO 10
      IF(A.LE.C14.AND.A+V.LE.C32) THEN
         VV=V*V
         HH1=H1TAB(IV)+H0TAB(IV)*C11283
         HH2=H2TAB(IV)+HH1*C11283-H0TAB(IV)
         HH3=(ONE-H2TAB(IV))*C37613-HH1*C23*VV+HH2*C11283
         HH4=(THREE*HH3-HH1)*C37613+H0TAB(IV)*C23*VV*VV
         VOIGTK=((((HH4*A+HH3)*A+HH2)*A+HH1)*A+H0TAB(IV))*
     *          (((CV1*A+CV2)*A+CV3)*A+CV4)
       ELSE
         AA=A*A
         VV=V*V
         U=(AA+VV)*C14142
         UU=U*U
         VOIGTK=((((AA-TEN*VV)*AA*THREE+FIFTN*VV*VV)/UU+
     *          THREE*VV-AA)/UU+ONE)*A*C79788/U
      END IF
      RETURN
      END
