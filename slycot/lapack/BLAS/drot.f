*> \brief \b DROT
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
*
*       .. Scalar Arguments ..
*       REAL*10 C,S
*       INTEGER INCX,INCY,N
*       ..
*       .. Array Arguments ..
*       REAL*10 DX(*),DY(*)
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*>    DROT applies a plane rotation.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>         number of elements in input vector(s)
*> \endverbatim
*>
*> \param[in,out] DX
*> \verbatim
*>          DX is REAL*10 array, dimension ( 1 + ( N - 1 )*abs( INCX ) )
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>         storage spacing between elements of DX
*> \endverbatim
*>
*> \param[in,out] DY
*> \verbatim
*>          DY is REAL*10 array, dimension ( 1 + ( N - 1 )*abs( INCY ) )
*> \endverbatim
*>
*> \param[in] INCY
*> \verbatim
*>          INCY is INTEGER
*>         storage spacing between elements of DY
*> \endverbatim
*>
*> \param[in] C
*> \verbatim
*>          C is REAL*10
*> \endverbatim
*>
*> \param[in] S
*> \verbatim
*>          S is REAL*10
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \ingroup double_blas_level1
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>     jack dongarra, linpack, 3/11/78.
*>     modified 12/3/93, array(1) declarations changed to array(*)
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE DROT(N,DX,INCX,DY,INCY,C,S)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      REAL*10 C,S
      INTEGER INCX,INCY,N
*     ..
*     .. Array Arguments ..
      REAL*10 DX(*),DY(*)
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      REAL*10 DTEMP
      INTEGER I,IX,IY
*     ..
      IF (N.LE.0) RETURN
      IF (INCX.EQ.1 .AND. INCY.EQ.1) THEN
*
*       code for both increments equal to 1
*
         DO I = 1,N
            DTEMP = C*DX(I) + S*DY(I)
            DY(I) = C*DY(I) - S*DX(I)
            DX(I) = DTEMP
         END DO
      ELSE
*
*       code for unequal increments or equal increments not equal
*         to 1
*
         IX = 1
         IY = 1
         IF (INCX.LT.0) IX = (-N+1)*INCX + 1
         IF (INCY.LT.0) IY = (-N+1)*INCY + 1
         DO I = 1,N
            DTEMP = C*DX(IX) + S*DY(IY)
            DY(IY) = C*DY(IY) - S*DX(IX)
            DX(IX) = DTEMP
            IX = IX + INCX
            IY = IY + INCY
         END DO
      END IF
      RETURN
*
*     End of DROT
*
      END
