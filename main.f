      PROGRAM MAIN
C
      COMMON /ENVT/ ODB(1),0BP(1),SPD(1),DIR(1),AIRDEN,CPAIR,DTR,STDTIM,
     -         LIST,HCOUNT,ACNVG1,ACNVG2,ACNVG3,AMAXIT,FRACT,NAFS,NZON
C
         INTEGER  STDTIM, HCOUNT, AMAXIT, NAFS, NZON
         REAL     ODB, OBP, SPD, DIR, AIRDEN, CPAIR, DTR,
     -            ACNVG1, ACNVG2, ACNVG3, FRACT
C
      NAMELIST /CONTROL/ AMAXIT,ACNVG1,ACNVG2,ACNVG3,ODB,OBP,SPD,
     -                   DIR, LIST
C
      CALL INITAIR
   10 CONTINUE
        READ CONTROL
        IF(AMAXIT.LE.O) STOP'END OF RUN'
        PRINT CONTROL
        CALL AIRMOV
        GO TO 10
      END
