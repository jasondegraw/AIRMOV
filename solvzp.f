      SUBROUTINE SOLVZP
         INTEGER  MAXAFS
      PARAMETER  (MAXAFS=128)
C
         INTEGER  MAXZON
      PARAMETER  (MAXZ0N=36)
C
      COMMON /AFSL/ FAREA(MAXAFS),FEXP(MAXAFS),ZS(MAXAFS),ZT(MAXAFS),
     -              PW(MAXAFS),PS(MAXAFS),FAZM(MAXAFS),AFSPTR(2,MAXAFS)
C
         INTEGER  AFSPTR
         REAL     FAREA, FEXP, ZS, PW, PS, FAZM, ZT
C
      COMMON /ENVT/ ODB(1),OBP(1),SPD(1),DIR(1),AIRDEN,CPAIR,DTR,STDTIM,
     -         LIST,HCOUNT,ACNVG1,ACNVG2,ACNVG3,AMAXIT,FRACT,NAFS,NZON
C
         INTEGER  STDTIM, HCOUNT, AMAXIT, NAFS, NZON
         REAL     ODB, OBP, SPD, DIR, AIRDEN, CPAIR, DTR,
     -            ACNVG1, ACNVG2, ACNVG3, FRACT
C
      COMMON /ZONL/ MCPM(MAXZON),MCPTM(MAXZON),TZ(0:MAXZON),
     -               ZZ(0:MAXZON),SQRTDZ(0:MAXZON),PZ(0:MAXZ0N),
     -               FAHS(MAXZON),SUMAF(MAXZON)
C
         REAL      TZ, ZZ, SQRTDZ, PZ, FAHS, MCPM, MCPTM, SUMAF
C
          INTEGER   TCOUNT
          LOGICAL   CNVG
          REAL      AA(MAXZ0N,MAXZON),BB(MAXZON),CC(MAXZON),DD(MAXZ0N)
      IF(LIST.GE.1) PRINT *, 'INITIALIZATION'
      TC0UNT=0
      IF(HCOUNT.EQ.1) GO TO 30
      DO 10 N=1,NZ0N
        DD(N)=1.0
        BB(N)=FAHS(N)
        DO 10 M=1,NZON
   10     AA(N,M)=0.0
      DO 20 I=1,NAFS
        N=AFSPTR(1,I)
        M=AFSPTR(2,I)
        AA(N,N)=AA(N,N)-FAREA(I)
        BB(N)=BB(N)-FAREA(I)*(PS(I)+PW(I))
        IF(M.LE.O) GO TO 20
        AA(M,M)=AA(M,M)-FAREA(I)
        AA(N,M)=AA(N,M)+FAREA(I)
        AA(M,N)=AA(M,N)+FAREA(I)
        BB(M)+BB(M)-FAREA(I)*(PS(I)+PW(I))
   20   CONTINUE
      IF(LIST.GE.2) CALL DUMPAB(AA,BB,NZ0N,MAXZON)
      CALL GAUSSY(AA,BB,CC,NZON,MAXZON)
      DO 25 N=1,NZ0N
        PZ(N)=CC(N)
        IF(LIST.GE.1) PRINT *, 'PZ: ',N,PZ(N)
   25   CONTINUE
C                            NEWTON ITERATION.
   30 CONTINUE
      TCOUNT=TCOUNT+l
      IF(LIST.GE.1) PRINT *, 'BEGIN ITERATION ',TCOUNT
      CNVG=.FALSE.
      DO 40 N=1,NZON
        BB(N)=FAHS(N)
        SUMAF(N)=0.0
        DO 40 M=1,NZON
          AA(M,N)=0.0
   40     CONTINUE
C                            EVALUATE FUNCTIONS AND
C                            PARTIAL DERIVATIVES.
      DO 50 I=1,NAFS
        N=AFSPTR(1,I)
        M=AFSPTR(2,I)
        DP=PZ(M)-PZ(N)+PS(I)+PW(I)
        IF(ABS(DP).LT.ACNVG3*
     -    (ABS(PZ(M))+ABS(PZ(N))+ABS(PS(I))+ABS(PW(I)))) GO TO 50
        IF(DP.LT.O.O) THEN
          F=-FAREA(I)*SQRTDZ(N)*(-DP)**FEXP(I)
        ELSE
          F=FAREA(I)*SQRTDZ(M)*DP**FEXP(I)
        END IF
        DF=F*FEXP(I)/DP
        BB(N)=BB(N)+F
        SUMAF(N)=SUMAF(N)+ABS(F)
        AA(N,N)=AA(N,N)-DF
        IF(M.LE.O) GO TO 50
        BB(M)=BB(M)-F
        SUMAF(M)=SUMAF(M)+ABS(F)
        AA(M,M)=AA(M,M)-DF
        AA(N,M)=AA(N,M)+DF
        AA(M,N)=AA(M,N)+DF
   50   CONTINUE
      DO 60 N=1,NZON
        IF(ABS(BB(N)).LE.ACNVG2) GO TO 60
        IF(ABS(BB(N)/SUMAF(N)).GT.ACNVG1) GO TO 70
   60   CONTINUE
      CNVG=.TRUE.
   70 CONTINUE
      IF(LIST.GE.2) CALL DUMPAB(AA,BB,NZ0N,MAXZON)
C                            CHECK CONVERGENCE.
      IF(CNVG) GO TO 999
      IF(TCOUNT.GT.AMAXIT) STOP'ITERATIONS'
C                            SOLVE AA * CC = BB.
      CALL GAUSSY(AA,BB,CC,NZON,MAXZON)
C                            IMPROVE CONVERGENCE.
      DO 80 N=1,NZON
        IF(CC(N)/DD(N).LE.-0.5) CC(N)=CC(N)*FRACT
        IF(ABS(CC(N)).GT.ACNVG3) DD(N)=CC(N)
   80   CONTINUE
C                            REVISE PZ.
      DO 90 N=1,NZON
        PZ(N)=PZ(N)-CC(N)
        IF(LIST.GE.l) PRINT *, 'REVIS: ',N,DD(N),CC(N),PZ(N),BB(N)
     -   ,SUMAF(N)
   90   CONTINUE
      GO TO 30
C
  999 CONTINUE
      PRINT *, 'ITERATIONS ',TCOUNT
      RETURN
      END