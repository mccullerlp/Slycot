*> \brief \b ZLARFG generates an elementary reflector (Householder matrix).
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZLARFG + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/zlarfg.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/zlarfg.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/zlarfg.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZLARFG( N, ALPHA, X, INCX, TAU )
*
*       .. Scalar Arguments ..
*       INTEGER            INCX, N
*       COMPLEX*16         ALPHA, TAU
*       ..
*       .. Array Arguments ..
*       COMPLEX*16         X( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZLARFG generates a complex elementary reflector H of order n, such
*> that
*>
*>       H**H * ( alpha ) = ( beta ),   H**H * H = I.
*>              (   x   )   (   0  )
*>
*> where alpha and beta are scalars, with beta real, and x is an
*> (n-1)-element complex vector. H is represented in the form
*>
*>       H = I - tau * ( 1 ) * ( 1 v**H ) ,
*>                     ( v )
*>
*> where tau is a complex scalar and v is a complex (n-1)-element
*> vector. Note that H is not hermitian.
*>
*> If the elements of x are all zero and alpha is real, then tau = 0
*> and H is taken to be the unit matrix.
*>
*> Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the elementary reflector.
*> \endverbatim
*>
*> \param[in,out] ALPHA
*> \verbatim
*>          ALPHA is COMPLEX*16
*>          On entry, the value alpha.
*>          On exit, it is overwritten with the value beta.
*> \endverbatim
*>
*> \param[in,out] X
*> \verbatim
*>          X is COMPLEX*16 array, dimension
*>                         (1+(N-2)*abs(INCX))
*>          On entry, the vector x.
*>          On exit, it is overwritten with the vector v.
*> \endverbatim
*>
*> \param[in] INCX
*> \verbatim
*>          INCX is INTEGER
*>          The increment between elements of X. INCX > 0.
*> \endverbatim
*>
*> \param[out] TAU
*> \verbatim
*>          TAU is COMPLEX*16
*>          The value tau.
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
*> \ingroup larfg
*
*  =====================================================================
      SUBROUTINE zlarfg( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      COMPLEX*16         ALPHA, TAU
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      parameter( one = 1.0d+0, zero = 0.0d+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY3, DZNRM2
      COMPLEX*16         ZLADIV
      EXTERNAL           dlamch, dlapy3, dznrm2, zladiv
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          abs, dble, dcmplx, imagpart, sign
*     ..
*     .. External Subroutines ..
      EXTERNAL           zdscal, zscal
*     ..
*     .. Executable Statements ..
*
      IF( n.LE.0 ) THEN
         tau = zero
         RETURN
      END IF
*
      xnorm = dznrm2( n-1, x, incx )
      alphr = dble( alpha )
      alphi = imagpart( alpha )
*
      IF( xnorm.EQ.zero .AND. alphi.EQ.zero ) THEN
*
*        H  =  I
*
         tau = zero
      ELSE
*
*        general case
*
         beta = -sign( dlapy3( alphr, alphi, xnorm ), alphr )
         safmin = dlamch( 'S' ) / dlamch( 'E' )
         rsafmn = one / safmin
*
         knt = 0
         IF( abs( beta ).LT.safmin ) THEN
*
*           XNORM, BETA may be inaccurate; scale X and recompute them
*
   10       CONTINUE
            knt = knt + 1
            CALL zdscal( n-1, rsafmn, x, incx )
            beta = beta*rsafmn
            alphi = alphi*rsafmn
            alphr = alphr*rsafmn
            IF( (abs( beta ).LT.safmin) .AND. (knt .LT. 20) )
     $         GO TO 10
*
*           New BETA is at most 1, at least SAFMIN
*
            xnorm = dznrm2( n-1, x, incx )
            alpha = dcmplx( alphr, alphi )
            beta = -sign( dlapy3( alphr, alphi, xnorm ), alphr )
         END IF
         tau = dcmplx( ( beta-alphr ) / beta, -alphi / beta )
         alpha = zladiv( dcmplx( one ), alpha-beta )
         CALL zscal( n-1, alpha, x, incx )
*
*        If ALPHA is subnormal, it may lose relative accuracy
*
         DO 20 j = 1, knt
            beta = beta*safmin
 20      CONTINUE
         alpha = beta
      END IF
*
      RETURN
*
*     End of ZLARFG
*
      END
