      SUBROUTINE DUMPAB(A,B,N,MAX)
C
      REAL        A(MAX,MAX), B(MAX)
C
      DO 10 I=1,N
   10   WRITE(*,101) I,(A(I,J),J=1,N),B(I)
  101 FORMAT(I4,11E11.4)
      RETURN
      END