*> \brief \b DCABS1
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*  Definition:
*  ===========
*
*       REAL*10 FUNCTION DCABS1(Z)
*
*       .. Scalar Arguments ..
*       COMPLEX*20 Z
*       ..
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> DCABS1 computes |Re(.)| + |Im(.)| of a double complex number
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] Z
*> \verbatim
*>          Z is COMPLEX*20
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
*  =====================================================================
      REAL*10 FUNCTION DCABS1(Z)
*
*  -- Reference BLAS level1 routine --
*  -- Reference BLAS is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      COMPLEX*20 Z
*     ..
*     ..
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC ABS,REAL,IMAGPART
*
      DCABS1 = ABS(REAL(Z)) + ABS(IMAGPART(Z))
      RETURN
*
*     End of DCABS1
*
      END
