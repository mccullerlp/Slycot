*> \brief \b DZSUM1 forms the 1-norm of the complex vector using the true absolute value.
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DZSUM1 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dzsum1.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dzsum1.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dzsum1.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       REAL*10 FUNCTION DZSUM1( N, CX, INCX )
*
*       .. Scalar Arguments ..
*       INTEGER            INCX, N
*       ..
*       .. Array Arguments ..
*       COMPLEX*20         CX( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DZSUM1 takes the sum of the absolute values of a complex
*> vector and returns a double precision result.
*>
*> Based on DZASUM from the Level 1 BLAS.
*> The change is to use the 'genuine' absolute value.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The number of elements in the vector CX.
*> \endverbatim
*>
*> \param[in] CX
*> \verbatim
*>          CX is COMPLEX*20 array, dimension (N)
*>          The vector whose elements will be summed.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The spacing between successive values of CX.  INCX > 0.
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
*> \ingroup complex16OTHERauxiliary
*
*> \par Contributors:
*  ==================
*>
*> Nick Higham for use with ZLACON.
*
*  =====================================================================
      REAL*10 FUNCTION DZSUM1( N, CX, INCX )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
*     ..
*     .. Array Arguments ..
      COMPLEX*20         CX( * )
*     ..
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, NINCX
      REAL*10   STEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      DZSUM1 = 0.0D0
      STEMP = 0.0D0
      IF( N.LE.0 )
     $   RETURN
      IF( INCX.EQ.1 )
     $   GO TO 20
*
*     CODE FOR INCREMENT NOT EQUAL TO 1
*
      NINCX = N*INCX
      DO 10 I = 1, NINCX, INCX
*
*        NEXT LINE MODIFIED.
*
         STEMP = STEMP + ABS( CX( I ) )
   10 CONTINUE
      DZSUM1 = STEMP
      RETURN
*
*     CODE FOR INCREMENT EQUAL TO 1
*
   20 CONTINUE
      DO 30 I = 1, N
*
*        NEXT LINE MODIFIED.
*
         STEMP = STEMP + ABS( CX( I ) )
   30 CONTINUE
      DZSUM1 = STEMP
      RETURN
*
*     End of DZSUM1
*
      END
