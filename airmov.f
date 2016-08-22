      SUBROUTINE AIRMOV
C
         REAL MCP
C
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
     -              ZZ(0:MAXZON),SQRTDZ(0:MAXZON),PZ(0:MAXZON),
     -              FAHS(MAXZON),SUMAF(MAXZON)
C
         REAL     TZ, ZZ, SORTDZ, PZ, FAHS, MCPM, MCPTM, SUMAF
C
C                            COMPUTE ZONE AIR DENSITIES.
      DO 10 N=0,NZON
        SQRTDZ(N)=SQRT(0.0034838*OBP(STDTIM)/(TZ(N)+273.15))
   10   CONTINUE
C                            COMPUTE CONSTANT DELTA-P.
      DO 20 1=1,NAFS
        N=AFSPTR(1,1)
        M=AFSPTR(2,1)
        IF(M.E0.0 .AND. HCOUNT.EQ.l) THEN
C                              WIND PRESSURE.
          V=SPD(STDTIM)*(0.1*ZT(I))**0.143
          W=0.5*AIRDEN*V*V
          Y=AMAX1(FAZM(I),DIR(STDTIM))-AMIN1(FAZM(I),DIR(STDTIM))
          IF(Y.GT.180.) Y=360.-Y
          IF(Y.LE.90.) X=0.75-1.05/90.*Y
          IF(Y.GT.90.) X=0.15/90.*Y-0.45
          PW(I)=X*W
        END IF
C                              STACK PRESSURE.
        PS(I)=9.80*((ZZ(M)-ZS(I))*SQRTDZ(M)*SQRTDZ(M)-
     -             (ZZ(N)-ZS(I))*SQRTDZ(N)*SQRTDZ(N))
      IF(LIST.GE.1) PRINT *, 'CONSTS: ',I,PS(I),PW(I)
   20   CONTINUE
C                            COMPUTE ZONE PRESSURES TO CONVERGENCE.
      CALL SOLVZP
C                            COMPUTE ZONE FLOWS.
      DO 30 N=1,NZON
        MCPM(N)=0.0
        MCPTM(N)=0.0
   30   CONTINUE
      DO 40 1=1 ,NAFS
        N=AFSPTR(1,I)
        M=AFSPTR(2,I)
        DP=PS(I)+PW(I)+PZ(M)-PZ(N)
        IF(DP.GT.O.O) THEN
          MCP=CPAIR*FAREA(I)*SQRTDZ(M)*DP**FEXP(I)
          MCPM(N)=MCPM(N)+MCP
          MCPTM(N)=MCPTM(N)+MCP*TZ(M)
        ELSE IF(M.GT.O) THEN
          MCP=CPAIR*FAREA(I)*SQRTDZ(N)*(-DP)**FEXP(I)
          MCPM(M)=MCPM(M)+MCP
          MCPTM(M)=MCPTM(M)+MCP*TZ(N)
        END IF
   40   CONTINUE
      DO 50 N=1,NZON
        IF(LIST.GE.l) PRINT *, 'MCPM: ',N,PZ(N),MCPM(N)
   50   CONTINUE
C                            IF (LIST.GE.1) PRINT FLOWS.
      DO 60 1=1,NAFS
        N=AFSPTR(1,I)
        M=AFSPTR(2,I)
        F=0.0
        DP=PS(I)+PW(I)+PZ(M)-PZ(N)
        IF(DP.LT.O.O) THEN
          F=-FAREA(I)*SQRTDZ(N)*(-DP)**FEXP(I)
          IF(LIST.GE.l) PRINT *, 'FLOW: ',I,M,FAREA(I),SQRTDZ(N),DP,F
        ELSE IF(DP.GT.O.O) THEN
          F=FAREA(I)*SORTDZ(M)*DP**FEXP(I)
          IF(LIST.GE.l) PRINT *, 'FLOW: ',I,M,FAREA(I),SQRTDZ(M),DP,F
        END IF
   60   CONTINUE
      RETURN
      END