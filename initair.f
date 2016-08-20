      SUBROUTINE INITAIR
C
         REAL     Z, T, F
C
         INTEGER  MAXAFS
      PARAMETER  (MAXAFS=128)
C
         INTEGER  MAXZ0N
      PARAMETER  (MAXZON=36)
C
      COMMON /AFSL/ FAREA(MAXAFS),FEXP(MAXAFS),ZS(MAXAFS),ZT(MAXAFS),
    -               PW(MAXAFS),PS(MAXAFS),FAZM(MAXAFS),AFSPTR(2,MAXAFS)
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
     -              ZZ(0:MAXZON),SORTDZ(0:MAXZON),PZ(0:MAXZON),
     -              FAHS(MAXZON),SUMAF(MAXZON)
C
         REAL     TZ, ZZ, SQRTDZ, PZ, FAHS, MCPM, MCPTM, SUMAF
C
C                                      INITIALIZE VARIABLES.
      HC0UNT=1
      STDTIM=1
      AMAXIT=20
      FRACT=.55
      LIST=1
      NAFS=0
      NZON=0
      DTR=0.0174533
      ACNVG1=.0001
      ACNVG2=.000005
      ACNVG3=.00001
      CPAIR=1004.
      ODB(1)=0.
      0BP(1)=100000.
      SPD(1)=5.
      DIR(1)=270.
      SQRT2=SQRT(2.0)
      WRITE(*,101)
  101 F0RMAT('1',8X,'N',6X,'ZZ',6X,'TZ',4X,'FAHS')
   10   CONTINUE
C                            READ ZONE DATA:
C                              HEIGHT, TEMP, SYSTEM FLOW
        READ *, Z,T,F
        IF (Z.LT.0.0) GO TO 20
        NZON=NZON+l
        WRITE(*,102) NZON,Z,T,F
  102   FORMAT(' ZON:',I4,3F8.2)
        ZZ(NZ0N)=Z
        TZ(NZ0N)=T
        FAHS(NZON)=F
        GO TO 10
   20 CONTINUE
      WRITE(*,103)
  103 FORMAT('O',8X,'I   N   M',7X,'A',7X,'X',7X,'C',5X,'AZM',
     -  6X,'ZS',6X,'ZT')
   30 CONTINUE
C                            READ OPENING DATA:
C                              NEAR SIDE ZONE, FAR SIDE ZONE,
C                              AREA, EXPONENT, FLOW COEFFICIENT,
C                              AZIMUTH, HEIGHT, SURFACE HEIGHT
        READ *, N,M,A,X,C,AZ,Z1,Z2
        IF (N.LE.O) GO TO 40
        NAFS=NAFS+1
        WRITE(*, 104) NAFS,N,M,A,X,C,AZ,Z1,Z2
  104   FORMAT(' AFS:',3I4,3F8.4,3F8.2)
        AFSPTR(1,NAFS)=N
        AFSPTR(2,NAFS)=M
        FAREA(NAFS)=SQRT2*C*A
        FEXP(NAFS)=X
        FAZM(NAFS)=AZ
        ZS(NAFS)=Z1
        ZT(NAFS)=Z2
        GO TO 30
   40 CONTINUE
C                             COMPUTE VARIABLES.
      TZ(0)=ODB(1)
      AIRDEN=0.0034838*OBP(1)/( 0DB(1)+273.15)
      ZZ(0)=0.0
      PZ(0)=0.0
      DO 50 1=1,NAFS
        PS(I)=0.0
        PW(I)=0.0
   50   CONTINUE
      RETURN
      END
