      SUBROUTINE CONVFF( LDA, A, E, F )
      
      DOUBLE PRECISION   A(LDA, *)
      REAL*16            E(LDA, *)
      REAL*10            F(LDA, *)
      
      DO 20 I = 1, LDA
         DO 10 J = 1, LDA
            A(I, J)=F(I, J)
            E(I, J)=F(I, J)
 10      CONTINUE
 20   CONTINUE
C
      RETURN
      END
