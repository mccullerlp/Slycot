      SUBROUTINE CONVF( LDA, A, E, F )
      
      DOUBLE PRECISION   A(LDA, LDA)
      REAL*16            E(LDA, LDA)
      REAL*10            F(LDA, LDA)
      
      DO 20 I = 1, LDA
         DO 10 J = 1, LDA
          E(I, J)=A(I, J)
          F(I, J)=A(I, J)
10       CONTINUE
20    CONTINUE
C
      RETURN
      END
