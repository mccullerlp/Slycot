      SUBROUTINE CONVF( LDA, A, E, F )
      
      DOUBLE PRECISION   A(LDA)
      REAL*16            E(LDA)
      REAL*10            F(LDA)
      
      DO 10 I = 1, LDA
          E(I)=A(I)
          F(I)=A(I)
10    CONTINUE
C
      RETURN
      END
