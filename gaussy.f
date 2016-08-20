      SUBROUTINE GAUSSY(A,B,X,N,MAX)
C
         REAL     A(MAX,MAX), B(MAX), X(MAX)
C
C                            GAUSS ELIMINATION.
C                            NO PIVOTING.
      DO 10 K=1,N-1
        DO 10 I=K+1,N
          D=A(I,K)/A(K,K)
          B(I)=B(I)-B(K)*D
          DO 10 J=K+1,N
   10       A(I,J)=A(I,J)-A(K,J)*D
C                            BACK SUBSTITUTION.
      IF(ABS(A(N,N)).LT.1.E-12) THEN
        X(N)=0.0
      ELSE
        X(N)=B(N)/A(N,N)
      END IF
      DO 30 I=N-1,1,-1
        D=0.0
        DO 20 J=I+1,N
   20     D=D+A(I,J)*X(J)
   30   X(I)=(B(I)-D)/A(I,I)
C
      RETURN
      END