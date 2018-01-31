C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
      SUBROUTINE DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, A, LDA, B, LDB,
     $                   VL, VU, IL, IU, ABSTOL, M, W, Z, LDZ, WORK,
     $                   LWORK, IWORK, IFAIL, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, ITYPE, IU, LDA, LDB, LDZ, LWORK, M, N
      DOUBLE PRECISION   ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            IFAIL( * ), IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * ), W( * ), WORK( * ),
     $                   Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  DSYGVX computes selected eigenvalues, and optionally, eigenvectors
*  of a real generalized symmetric-definite eigenproblem, of the form
*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A
*  and B are assumed to be symmetric and B is also positive definite.
*  Eigenvalues and eigenvectors can be selected by specifying either a
*  range of values or a range of indices for the desired eigenvalues.
*
*  Arguments
*  =========
*
*  ITYPE   (input) INTEGER
*          Specifies the problem type to be solved:
*          = 1:  A*x = (lambda)*B*x
*          = 2:  A*B*x = (lambda)*x
*          = 3:  B*A*x = (lambda)*x
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  RANGE   (input) CHARACTER*1
*          = 'A': all eigenvalues will be found.
*          = 'V': all eigenvalues in the half-open interval (VL,VU]
*                 will be found.
*          = 'I': the IL-th through IU-th eigenvalues will be found.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A and B are stored;
*          = 'L':  Lower triangle of A and B are stored.
*
*  N       (input) INTEGER
*          The order of the matrix pencil (A,B).  N >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA, N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*
*          On exit, the lower triangle (if UPLO='L') or the upper
*          triangle (if UPLO='U') of A, including the diagonal, is
*          destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  B       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA, N)
*          On entry, the symmetric matrix B.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of B contains the
*          upper triangular part of the matrix B.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of B contains
*          the lower triangular part of the matrix B.
*
*          On exit, if INFO <= N, the part of B containing the matrix is
*          overwritten by the triangular factor U or L from the Cholesky
*          factorization B = U**T*U or B = L*L**T.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  VL      (INPUT) DOUBLE PRECISION
*  VU      (INPUT) DOUBLE PRECISION
*          If RANGE='V', the lower and upper bounds of the interval to
*          be searched for eigenvalues. VL < VU.
*          Not referenced if RANGE = 'A' or 'I'.
*
*  IL      (input) INTEGER
*  IU      (input) INTEGER
*          If RANGE='I', the indices (in ascending order) of the
*          smallest and largest eigenvalues to be returned.
*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
*          Not referenced if RANGE = 'A' or 'V'.
*
*  ABSTOL  (INPUT) DOUBLE PRECISION
*          The absolute error tolerance for the eigenvalues.
*          An approximate eigenvalue is accepted as converged
*          when it is determined to lie in an interval [a,b]
*          of width less than or equal to
*
*                  ABSTOL + EPS *   max( |a|,|b| ) ,
*
*          where EPS is the machine precision.  If ABSTOL is less than
*          or equal to zero, then  EPS*|T|  will be used in its place,
*          where |T| is the 1-norm of the tridiagonal matrix obtained
*          by reducing A to tridiagonal form.
*
*          Eigenvalues will be computed most accurately when ABSTOL is
*          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
*          If this routine returns with INFO>0, indicating that some
*          eigenvectors did not converge, try setting ABSTOL to
*          2*DLAMCH('S').
*
*  M       (output) INTEGER
*          The total number of eigenvalues found.  0 <= M <= N.
*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
*
*  W       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          On normal exit, the first M elements contain the selected
*          eigenvalues in ascending order.
*
*  Z       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDZ, MAX(1,M))
*          If JOBZ = 'N', then Z is not referenced.
*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
*          contain the orthonormal eigenvectors of the matrix A
*          corresponding to the selected eigenvalues, with the i-th
*          column of Z holding the eigenvector associated with W(i).
*          The eigenvectors are normalized as follows:
*          if ITYPE = 1 or 2, Z**T*B*Z = I;
*          if ITYPE = 3, Z**T*inv(B)*Z = I.
*
*          If an eigenvector fails to converge, then that column of Z
*          contains the latest approximation to the eigenvector, and the
*          index of the eigenvector is returned in IFAIL.
*          Note: the user must ensure that at least max(1,M) columns are
*          supplied in the array Z; if RANGE = 'V', the exact value of M
*          is not known in advance and an upper bound must be used.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  WORK    (WORKSPACE/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= max(1,8*N).
*          For optimal efficiency, LWORK >= (NB+3)*N,
*          where NB is the blocksize for DSYTRD returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace) INTEGER array, dimension (5*N)
*
*  IFAIL   (output) INTEGER array, dimension (N)
*          If JOBZ = 'V', then if INFO = 0, the first M elements of
*          IFAIL are zero.  If INFO > 0, then IFAIL contains the
*          indices of the eigenvectors that failed to converge.
*          If JOBZ = 'N', then IFAIL is not referenced.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  DPOTRF or DSYEVX returned an error code:
*             <= N:  if INFO = i, DSYEVX failed to converge;
*                    i eigenvectors failed to converge.  Their indices
*                    are stored in array IFAIL.
*             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
*                    minor of order i of B is not positive definite.
*                    The factorization of B could not be completed and
*                    no eigenvalues or eigenvectors were computed.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Mark Fahey, Department of Mathematics, Univ. of Kentucky, USA
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, LQUERY, UPPER, VALEIG, WANTZ
      CHARACTER          TRANS
      INTEGER            LOPT, LWKOPT, NB
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DPOTRF, DSYEVX, DSYGST, DTRMM, DTRSM, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      UPPER = LSAME( UPLO, 'U' )
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
      LQUERY = ( LWORK.EQ.-1 )
*
      INFO = 0
      IF( ITYPE.LT.0 .OR. ITYPE.GT.3 ) THEN
         INFO = -1
      ELSE IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -3
      ELSE IF( .NOT.( UPPER .OR. LSAME( UPLO, 'L' ) ) ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( VALEIG .AND. N.GT.0 ) THEN
         IF( VU.LE.VL )
     $      INFO = -11
      ELSE IF( INDEIG .AND. IL.LT.1 ) THEN
         INFO = -12
      ELSE IF( INDEIG .AND. ( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) ) THEN
         INFO = -13
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -18
      ELSE IF( LWORK.LT.MAX( 1, 8*N ) .AND. .NOT.LQUERY ) THEN
         INFO = -20
      END IF
*
      IF( INFO.EQ.0 ) THEN
         NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
         LWKOPT = ( NB+3 )*N
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYGVX', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      M = 0
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
*     Form a Cholesky factorization of B.
*
      CALL DPOTRF( UPLO, N, B, LDB, INFO )
      IF( INFO.NE.0 ) THEN
         INFO = N + INFO
         RETURN
      END IF
*
*     Transform problem to standard eigenvalue problem and solve.
*
      CALL DSYGST( ITYPE, UPLO, N, A, LDA, B, LDB, INFO )
      CALL DSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU, ABSTOL,
     $             M, W, Z, LDZ, WORK, LWORK, IWORK, IFAIL, INFO )
      LOPT = WORK( 1 )
*
      IF( WANTZ ) THEN
*
*        Backtransform eigenvectors to the original problem.
*
         IF( INFO.GT.0 )
     $      M = INFO - 1
         IF( ITYPE.EQ.1 .OR. ITYPE.EQ.2 ) THEN
*
*           For A*x=(lambda)*B*x and A*B*x=(lambda)*x;
*           backtransform eigenvectors: x = inv(L)'*y or inv(U)*y
*
            IF( UPPER ) THEN
               TRANS = 'N'
            ELSE
               TRANS = 'T'
            END IF
*
            CALL DTRSM( 'Left', UPLO, TRANS, 'Non-unit', N, M, ONE, B,
     $                  LDB, Z, LDZ )
*
         ELSE IF( ITYPE.EQ.3 ) THEN
*
*           For B*A*x=(lambda)*x;
*           backtransform eigenvectors: x = L*y or U'*y
*
            IF( UPPER ) THEN
               TRANS = 'T'
            ELSE
               TRANS = 'N'
            END IF
*
            CALL DTRMM( 'Left', UPLO, TRANS, 'Non-unit', N, M, ONE, B,
     $                  LDB, Z, LDZ )
         END IF
      END IF
*
*     Set WORK(1) to optimal workspace size.
*
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of DSYGVX
*
      END
      SUBROUTINE DSYEVX( JOBZ, RANGE, UPLO, N, A, LDA, VL, VU, IL, IU,
     $                   ABSTOL, M, W, Z, LDZ, WORK, LWORK, IWORK,
     $                   IFAIL, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, RANGE, UPLO
      INTEGER            IL, INFO, IU, LDA, LDZ, LWORK, M, N
      DOUBLE PRECISION   ABSTOL, VL, VU
*     ..
*     .. Array Arguments ..
      INTEGER            IFAIL( * ), IWORK( * )
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  DSYEVX computes selected eigenvalues and, optionally, eigenvectors
*  of a real symmetric matrix A.  Eigenvalues and eigenvectors can be
*  selected by specifying either a range of values or a range of indices
*  for the desired eigenvalues.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  RANGE   (input) CHARACTER*1
*          = 'A': all eigenvalues will be found.
*          = 'V': all eigenvalues in the half-open interval (VL,VU]
*                 will be found.
*          = 'I': the IL-th through IU-th eigenvalues will be found.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (INPUT/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA, N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*          On exit, the lower triangle (if UPLO='L') or the upper
*          triangle (if UPLO='U') of A, including the diagonal, is
*          destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  VL      (INPUT) DOUBLE PRECISION
*  VU      (INPUT) DOUBLE PRECISION
*          If RANGE='V', the lower and upper bounds of the interval to
*          be searched for eigenvalues. VL < VU.
*          Not referenced if RANGE = 'A' or 'I'.
*
*  IL      (input) INTEGER
*  IU      (input) INTEGER
*          If RANGE='I', the indices (in ascending order) of the
*          smallest and largest eigenvalues to be returned.
*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
*          Not referenced if RANGE = 'A' or 'V'.
*
*  ABSTOL  (INPUT) DOUBLE PRECISION
*          The absolute error tolerance for the eigenvalues.
*          An approximate eigenvalue is accepted as converged
*          when it is determined to lie in an interval [a,b]
*          of width less than or equal to
*
*                  ABSTOL + EPS *   max( |a|,|b| ) ,
*
*          where EPS is the machine precision.  If ABSTOL is less than
*          or equal to zero, then  EPS*|T|  will be used in its place,
*          where |T| is the 1-norm of the tridiagonal matrix obtained
*          by reducing A to tridiagonal form.
*
*          Eigenvalues will be computed most accurately when ABSTOL is
*          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
*          If this routine returns with INFO>0, indicating that some
*          eigenvectors did not converge, try setting ABSTOL to
*          2*DLAMCH('S').
*
*          See "Computing Small Singular Values of Bidiagonal Matrices
*          with Guaranteed High Relative Accuracy," by Demmel and
*          Kahan, LAPACK Working Note #3.
*
*  M       (output) INTEGER
*          The total number of eigenvalues found.  0 <= M <= N.
*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
*
*  W       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (N)
*          On normal exit, the first M elements contain the selected
*          eigenvalues in ascending order.
*
*  Z       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDZ, MAX(1,M))
*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
*          contain the orthonormal eigenvectors of the matrix A
*          corresponding to the selected eigenvalues, with the i-th
*          column of Z holding the eigenvector associated with W(i).
*          If an eigenvector fails to converge, then that column of Z
*          contains the latest approximation to the eigenvector, and the
*          index of the eigenvector is returned in IFAIL.
*          If JOBZ = 'N', then Z is not referenced.
*          Note: the user must ensure that at least max(1,M) columns are
*          supplied in the array Z; if RANGE = 'V', the exact value of M
*          is not known in advance and an upper bound must be used.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  WORK    (WORKSPACE/OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= max(1,8*N).
*          For optimal efficiency, LWORK >= (NB+3)*N,
*          where NB is the max of the blocksize for DSYTRD and DORMTR
*          returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace) INTEGER array, dimension (5*N)
*
*  IFAIL   (output) INTEGER array, dimension (N)
*          If JOBZ = 'V', then if INFO = 0, the first M elements of
*          IFAIL are zero.  If INFO > 0, then IFAIL contains the
*          indices of the eigenvectors that failed to converge.
*          If JOBZ = 'N', then IFAIL is not referenced.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, then i eigenvectors failed to converge.
*                Their indices are stored in array IFAIL.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ALLEIG, INDEIG, LOWER, LQUERY, VALEIG, WANTZ
      CHARACTER          ORDER
      INTEGER            I, IINFO, IMAX, INDD, INDE, INDEE, INDIBL,
     $                   INDISP, INDIWO, INDTAU, INDWKN, INDWRK, ISCALE,
     $                   ITMP1, J, JJ, LLWORK, LLWRKN, LOPT, LWKOPT, NB,
     $                   NSPLIT
      DOUBLE PRECISION   ABSTLL, ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN,
     $                   SIGMA, SMLNUM, TMP1, VLL, VUU
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANSY
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DCOPY, DLACPY, DORGTR, DORMTR, DSCAL, DSTEBZ,
     $                   DSTEIN, DSTEQR, DSTERF, DSWAP, DSYTRD, XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      LOWER = LSAME( UPLO, 'L' )
      WANTZ = LSAME( JOBZ, 'V' )
      ALLEIG = LSAME( RANGE, 'A' )
      VALEIG = LSAME( RANGE, 'V' )
      INDEIG = LSAME( RANGE, 'I' )
      LQUERY = ( LWORK.EQ.-1 )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( ALLEIG .OR. VALEIG .OR. INDEIG ) ) THEN
         INFO = -2
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE
         IF( VALEIG ) THEN
            IF( N.GT.0 .AND. VU.LE.VL )
     $         INFO = -8
         ELSE IF( INDEIG ) THEN
            IF( IL.LT.1 .OR. IL.GT.MAX( 1, N ) ) THEN
               INFO = -9
            ELSE IF( IU.LT.MIN( N, IL ) .OR. IU.GT.N ) THEN
               INFO = -10
            END IF
         END IF
      END IF
      IF( INFO.EQ.0 ) THEN
         IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
            INFO = -15
         ELSE IF( LWORK.LT.MAX( 1, 8*N ) .AND. .NOT.LQUERY ) THEN
            INFO = -17
         END IF
      END IF
*
      IF( INFO.EQ.0 ) THEN
         NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
         NB = MAX( NB, ILAENV( 1, 'DORMTR', UPLO, N, -1, -1, -1 ) )
         LWKOPT = ( NB+3 )*N
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYEVX', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      M = 0
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      IF( N.EQ.1 ) THEN
         WORK( 1 ) = 7
         IF( ALLEIG .OR. INDEIG ) THEN
            M = 1
            W( 1 ) = A( 1, 1 )
         ELSE
            IF( VL.LT.A( 1, 1 ) .AND. VU.GE.A( 1, 1 ) ) THEN
               M = 1
               W( 1 ) = A( 1, 1 )
            END IF
         END IF
         IF( WANTZ )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = MIN( SQRT( BIGNUM ), ONE / SQRT( SQRT( SAFMIN ) ) )
*
*     Scale matrix to allowable range, if necessary.
*
      ISCALE = 0
      ABSTLL = ABSTOL
      VLL = VL
      VUU = VU
      ANRM = DLANSY( 'M', UPLO, N, A, LDA, WORK )
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         IF( LOWER ) THEN
            DO 10 J = 1, N
               CALL DSCAL( N-J+1, SIGMA, A( J, J ), 1 )
   10       CONTINUE
         ELSE
            DO 20 J = 1, N
               CALL DSCAL( J, SIGMA, A( 1, J ), 1 )
   20       CONTINUE
         END IF
         IF( ABSTOL.GT.0 )
     $      ABSTLL = ABSTOL*SIGMA
         IF( VALEIG ) THEN
            VLL = VL*SIGMA
            VUU = VU*SIGMA
         END IF
      END IF
*
*     Call DSYTRD to reduce symmetric matrix to tridiagonal form.
*
      INDTAU = 1
      INDE = INDTAU + N
      INDD = INDE + N
      INDWRK = INDD + N
      LLWORK = LWORK - INDWRK + 1
      CALL DSYTRD( UPLO, N, A, LDA, WORK( INDD ), WORK( INDE ),
     $             WORK( INDTAU ), WORK( INDWRK ), LLWORK, IINFO )
      LOPT = 3*N + WORK( INDWRK )
*
*     If all eigenvalues are desired and ABSTOL is less than or equal to
*     zero, then call DSTERF or DORGTR and SSTEQR.  If this fails for
*     some eigenvalue, then try DSTEBZ.
*
      IF( ( ALLEIG .OR. ( INDEIG .AND. IL.EQ.1 .AND. IU.EQ.N ) ) .AND.
     $    ( ABSTOL.LE.ZERO ) ) THEN
         CALL DCOPY( N, WORK( INDD ), 1, W, 1 )
         INDEE = INDWRK + 2*N
         IF( .NOT.WANTZ ) THEN
            CALL DCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
            CALL DSTERF( N, W, WORK( INDEE ), INFO )
         ELSE
            CALL DLACPY( 'A', N, N, A, LDA, Z, LDZ )
            CALL DORGTR( UPLO, N, Z, LDZ, WORK( INDTAU ),
     $                   WORK( INDWRK ), LLWORK, IINFO )
            CALL DCOPY( N-1, WORK( INDE ), 1, WORK( INDEE ), 1 )
            CALL DSTEQR( JOBZ, N, W, WORK( INDEE ), Z, LDZ,
     $                   WORK( INDWRK ), INFO )
            IF( INFO.EQ.0 ) THEN
               DO 30 I = 1, N
                  IFAIL( I ) = 0
   30          CONTINUE
            END IF
         END IF
         IF( INFO.EQ.0 ) THEN
            M = N
            GO TO 40
         END IF
         INFO = 0
      END IF
*
*     Otherwise, call DSTEBZ and, if eigenvectors are desired, SSTEIN.
*
      IF( WANTZ ) THEN
         ORDER = 'B'
      ELSE
         ORDER = 'E'
      END IF
      INDIBL = 1
      INDISP = INDIBL + N
      INDIWO = INDISP + N
      CALL DSTEBZ( RANGE, ORDER, N, VLL, VUU, IL, IU, ABSTLL,
     $             WORK( INDD ), WORK( INDE ), M, NSPLIT, W,
     $             IWORK( INDIBL ), IWORK( INDISP ), WORK( INDWRK ),
     $             IWORK( INDIWO ), INFO )
*
      IF( WANTZ ) THEN
         CALL DSTEIN( N, WORK( INDD ), WORK( INDE ), M, W,
     $                IWORK( INDIBL ), IWORK( INDISP ), Z, LDZ,
     $                WORK( INDWRK ), IWORK( INDIWO ), IFAIL, INFO )
*
*        Apply orthogonal matrix used in reduction to tridiagonal
*        form to eigenvectors returned by DSTEIN.
*
         INDWKN = INDE
         LLWRKN = LWORK - INDWKN + 1
         CALL DORMTR( 'L', UPLO, 'N', N, M, A, LDA, WORK( INDTAU ), Z,
     $                LDZ, WORK( INDWKN ), LLWRKN, IINFO )
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
   40 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = M
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
*     If eigenvalues are not in order, then sort them, along with
*     eigenvectors.
*
      IF( WANTZ ) THEN
         DO 60 J = 1, M - 1
            I = 0
            TMP1 = W( J )
            DO 50 JJ = J + 1, M
               IF( W( JJ ).LT.TMP1 ) THEN
                  I = JJ
                  TMP1 = W( JJ )
               END IF
   50       CONTINUE
*
            IF( I.NE.0 ) THEN
               ITMP1 = IWORK( INDIBL+I-1 )
               W( I ) = W( J )
               IWORK( INDIBL+I-1 ) = IWORK( INDIBL+J-1 )
               W( J ) = TMP1
               IWORK( INDIBL+J-1 ) = ITMP1
               CALL DSWAP( N, Z( 1, I ), 1, Z( 1, J ), 1 )
               IF( INFO.NE.0 ) THEN
                  ITMP1 = IFAIL( I )
                  IFAIL( I ) = IFAIL( J )
                  IFAIL( J ) = ITMP1
               END IF
            END IF
   60    CONTINUE
      END IF
*
*     Set WORK(1) to optimal workspace size.
*
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of DSYEVX
*
      END
      SUBROUTINE DLACPY( UPLO, M, N, A, LDA, B, LDB )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  DLACPY copies all or part of a two-dimensional matrix A to another
*  matrix B.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the part of the matrix A to be copied to B.
*          = 'U':      Upper triangular part
*          = 'L':      Lower triangular part
*          Otherwise:  All of the matrix A
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (INPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDA,N)
*          The m by n matrix A.  If UPLO = 'U', only the upper triangle
*          or trapezoid is accessed; if UPLO = 'L', only the lower
*          triangle or trapezoid is accessed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (OUTPUT) DOUBLE PRECISION ARRAY, DIMENSION (LDB,N)
*          On exit, B = A in the locations specified by UPLO.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,M).
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
      RETURN
*
*     End of DLACPY
*
      END
