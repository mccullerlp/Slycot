      SUBROUTINE CONVFF( LDA, A, E, F )
      
      DOUBLE PRECISION   A(LDA)
      REAL*16            E(LDA)
      REAL*10            F(LDA)
      
      DO 10 I = 1, LDA
          A(I)=F(I)
          E(I)=F(I)
10    CONTINUE
C
      RETURN
      END