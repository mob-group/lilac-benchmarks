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
C   This subroutine makes a 3*NATOMS square matrix to project translations,
C   rotations and the direction of the gradient vector out of the Hessian.
C
      SUBROUTINE PROJH(Q,NATOMS,ATMASS,VNEW,PROJGRAD)
      USE MODHESS
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4, NATOMS, IPVT(6)
      DOUBLE PRECISION DUMMY, TMASS, DET(2), ZWORK(6), RCOND,
     &                 Q(3*NATOMS), ATMASS(NATOMS), VNEW(3*NATOMS), RMS,
     &                 BMAT(3*NATOMS,6), SMAT(6,6), PMAT(3*NATOMS,3*NATOMS), TEMP(3*NATOMS,6), TMAT(3*NATOMS,3*NATOMS)
      LOGICAL PROJGRAD

      BMAT(1:3*NATOMS,1:6)=0.0D0
      DO J1=1,NATOMS
         J3=3*(J1-1)
         DUMMY=SQRT(ATMASS(J1))
         VNEW(J3+1)=VNEW(J3+1)/SQRT(ATMASS(J1))
         VNEW(J3+2)=VNEW(J3+2)/SQRT(ATMASS(J1))
         VNEW(J3+3)=VNEW(J3+3)/SQRT(ATMASS(J1))
         BMAT(J3+1,1)=DUMMY
         BMAT(J3+2,2)=DUMMY
         BMAT(J3+3,3)=DUMMY
         BMAT(J3+2,4)= Q(J3+3)*DUMMY
         BMAT(J3+3,4)=-Q(J3+2)*DUMMY
         BMAT(J3+1,5)=-Q(J3+3)*DUMMY
         BMAT(J3+3,5)= Q(J3+1)*DUMMY
         BMAT(J3+1,6)= Q(J3+2)*DUMMY
         BMAT(J3+2,6)=-Q(J3+1)*DUMMY
      ENDDO
C
C  Normalising BMAT appears to be unnecessary.
C
      DO J1=1,6
         DUMMY=0.0D0
         DO J3=1,3*NATOMS
            DUMMY=DUMMY+BMAT(J3,J1)**2
         ENDDO
         DUMMY=SQRT(DUMMY)
         DO J3=1,3*NATOMS
            BMAT(J3,J1)=BMAT(J3,J1)/DUMMY
         ENDDO
         DUMMY=0.0D0
         DO J3=1,3*NATOMS
         ENDDO
!        PRINT '(A,I6,A,G20.10)',' projh> overlap of gradient vector with basis vector ',J1,' is: ',DUMMY
      ENDDO
      DO J1=1,6
         DO J2=1,6
            DUMMY=0.0D0
            DO J3=1,3*NATOMS
               DUMMY=DUMMY+BMAT(J3,J1)*BMAT(J3,J2)
            ENDDO
            SMAT(J1,J2)=DUMMY
         ENDDO
      ENDDO

!     PRINT '(A)',' projh> SMAT:'
!     PRINT '(6G20.10)',SMAT(1:6,1:6)

      CALL DGECO(SMAT,6,6,IPVT,RCOND,ZWORK)
      CALL DGEDI(SMAT,6,6,IPVT,DET,ZWORK,11)

!     PRINT '(A,G20.10)',' projh> SINV: DET=',DET(1)
!     PRINT '(6G20.10)',SMAT(1:6,1:6)
C
C  SMAT should now contain the inverse of the original SMAT
C
      DO J1=1,3*NATOMS
         DO J2=1,6
            TEMP(J1,J2)=BMAT(J1,1)*SMAT(J2,1)+BMAT(J1,2)*SMAT(J2,2)+BMAT(J1,3)*SMAT(J2,3)
     &                 +BMAT(J1,4)*SMAT(J2,4)+BMAT(J1,5)*SMAT(J2,5)+BMAT(J1,6)*SMAT(J2,6)
         ENDDO
      ENDDO
 
      DO J1=1,3*NATOMS
         DO J2=1,3*NATOMS
C
C  The next four lines should be unnecessary, but pgi miscompiles and gives SEGV
C  without them !!!!!!!
C
            DUMMY=0.0D0
            DO J3=1,6
               DUMMY=DUMMY+TEMP(J1,J3)*BMAT(J2,J3)
            ENDDO

            PMAT(J1,J2)=TEMP(J1,1)*BMAT(J2,1)+TEMP(J1,2)*BMAT(J2,2)+TEMP(J1,3)*BMAT(J2,3)
     &                 +TEMP(J1,4)*BMAT(J2,4)+TEMP(J1,5)*BMAT(J2,5)+TEMP(J1,6)*BMAT(J2,6)
         ENDDO
      ENDDO
    
      IF (PROJGRAD) THEN
         PRINT '(A)',' projh> projecting gradient direction from Hessian'
         DUMMY=0.0D0
         DO J1=1,3*NATOMS
            DUMMY=DUMMY+VNEW(J1)**2
         ENDDO
         DO J1=1,3*NATOMS
            DO J2=1,3*NATOMS
               PMAT(J2,J1)=PMAT(J2,J1)+VNEW(J2)*VNEW(J1)/DUMMY
            ENDDO
         ENDDO
      ENDIF
C
C  Make the matrix I - P
C
      PMAT(1:3*NATOMS,1:3*NATOMS)=-PMAT(1:3*NATOMS,1:3*NATOMS)
      DO J1=1,3*NATOMS
         PMAT(J1,J1)=1.0D0+PMAT(J1,J1)
      ENDDO
C
C  Project the Hessian. All the matrices should be symmetric.
C
      DO J1=1,3*NATOMS
         DO J2=1,3*NATOMS
            DUMMY=0.0D0
            DO J3=1,3*NATOMS
               DUMMY=DUMMY+HESS(J1,J3)*PMAT(J3,J2)
            ENDDO
            TMAT(J1,J2)=DUMMY
         ENDDO
      ENDDO
      DO J1=1,3*NATOMS
         DO J2=1,3*NATOMS
            DUMMY=0.0D0
            DO J3=1,3*NATOMS
               DUMMY=DUMMY+PMAT(J1,J3)*TMAT(J3,J2)
            ENDDO
            HESS(J1,J2)=DUMMY
         ENDDO
      ENDDO

      RETURN
      END

      subroutine dgedi(a,lda,n,ipvt,det,work,job)
      integer lda,n,ipvt(*),job
      double precision a(lda,*),det(2),work(*)
c
c     dgedi computes the determinant and inverse of a matrix
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        work    double precision(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     double precision(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. dabs(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if dgeco has set rcond .gt. 0.0 or dgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,dswap
c     fortran dabs,mod
c
c     internal variables
c
      double precision t
      double precision ten
      integer i,j,k,kb,kp1,l,nm1
c
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = 1.0d0
         det(2) = 0.0d0
         ten = 10.0d0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (det(1) .eq. 0.0d0) go to 60
   10       if (dabs(det(1)) .ge. 1.0d0) go to 20
               det(1) = ten*det(1)
               det(2) = det(2) - 1.0d0
            go to 10
   20       continue
   30       if (dabs(det(1)) .lt. ten) go to 40
               det(1) = det(1)/ten
               det(2) = det(2) + 1.0d0
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            a(k,k) = 1.0d0/a(k,k)
            t = -a(k,k)
            call dscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = 0.0d0
               call daxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = 0.0d0
  110       continue
            do 120 j = kp1, n
               t = work(j)
               call daxpy(n,t,a(1,j),1,a(1,k),1)
  120       continue
            l = ipvt(k)
            if (l .ne. k) call dswap(n,a(1,k),1,a(1,l),1)
  130    continue
  140    continue
  150 continue
      return
      end

      subroutine dgeco(a,lda,n,ipvt,rcond,z)
      integer lda,n,ipvt(*)
      double precision a(lda,*),z(*)
      double precision rcond
c
c     dgeco factors a double precision matrix by gaussian elimination
c     and estimates the condition of the matrix.
c
c     if  rcond  is not needed, dgefa is slightly faster.
c     to solve  a*x = b , follow dgeco by dgesl.
c     to compute  inverse(a)*c , follow dgeco by dgesl.
c     to compute  determinant(a) , follow dgeco by dgedi.
c     to compute  inverse(a) , follow dgeco by dgedi.
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        rcond   double precision
c                an estimate of the reciprocal condition of  a .
c                for the system  a*x = b , relative perturbations
c                in  a  and  b  of size  epsilon  may cause
c                relative perturbations in  x  of size  epsilon/rcond .
c                if  rcond  is so small that the logical expression
c                           1.0 + rcond .eq. 1.0
c                is true, then  a  may be singular to working
c                precision.  in particular,  rcond  is zero  if
c                exact singularity is detected or the estimate
c                underflows.
c
c        z       double precision(n)
c                a work vector whose contents are usually unimportant.
c                if  a  is close to a singular matrix, then  z  is
c                an approximate null vector in the sense that
c                norm(a*z) = rcond*norm(a)*norm(z) .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     linpack dgefa
c     blas daxpy,ddot,dscal,dasum
c     fortran dabs,dmax1,dsign
c
c     internal variables
c
      double precision ddot,ek,t,wk,wkm
      double precision anorm,s,dasum,sm,ynorm
      integer info,j,k,kb,kp1,l
c
c
c     compute 1-norm of a
c
      anorm = 0.0d0
      do 10 j = 1, n
         anorm = dmax1(anorm,dasum(n,a(1,j),1))
   10 continue
c
c     factor
c
      call dgefa(a,lda,n,ipvt,info)
c
c     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) .
c     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e .
c     trans(a)  is the transpose of a .  the components of  e  are
c     chosen to cause maximum local growth in the elements of w  where
c     trans(u)*w = e .  the vectors are frequently rescaled to avoid
c     overflow.
c
c     solve trans(u)*w = e
c
      ek = 1.0d0
      do 20 j = 1, n
         z(j) = 0.0d0
   20 continue
      do 100 k = 1, n
         if (z(k) .ne. 0.0d0) ek = dsign(ek,-z(k))
         if (dabs(ek-z(k)) .le. dabs(a(k,k))) go to 30
            s = dabs(a(k,k))/dabs(ek-z(k))
            call dscal(n,s,z,1)
            ek = s*ek
   30    continue
         wk = ek - z(k)
         wkm = -ek - z(k)
         s = dabs(wk)
         sm = dabs(wkm)
         if (a(k,k) .eq. 0.0d0) go to 40
            wk = wk/a(k,k)
            wkm = wkm/a(k,k)
         go to 50
   40    continue
            wk = 1.0d0
            wkm = 1.0d0
   50    continue
         kp1 = k + 1
         if (kp1 .gt. n) go to 90
            do 60 j = kp1, n
               sm = sm + dabs(z(j)+wkm*a(k,j))
               z(j) = z(j) + wk*a(k,j)
               s = s + dabs(z(j))
   60       continue
            if (s .ge. sm) go to 80
               t = wkm - wk
               wk = wkm
               do 70 j = kp1, n
                  z(j) = z(j) + t*a(k,j)
   70          continue
   80       continue
   90    continue
         z(k) = wk
  100 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
c     solve trans(l)*y = w
c
      do 120 kb = 1, n
         k = n + 1 - kb
         if (k .lt. n) z(k) = z(k) + ddot(n-k,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 110
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
  110    continue
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
  120 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
c
      ynorm = 1.0d0
c
c     solve l*v = y
c
      do 140 k = 1, n
         l = ipvt(k)
         t = z(l)
         z(l) = z(k)
         z(k) = t
         if (k .lt. n) call daxpy(n-k,t,a(k+1,k),1,z(k+1),1)
         if (dabs(z(k)) .le. 1.0d0) go to 130
            s = 1.0d0/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  130    continue
  140 continue
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
c     solve  u*z = v
c
      do 160 kb = 1, n
         k = n + 1 - kb
         if (dabs(z(k)) .le. dabs(a(k,k))) go to 150
            s = dabs(a(k,k))/dabs(z(k))
            call dscal(n,s,z,1)
            ynorm = s*ynorm
  150    continue
         if (a(k,k) .ne. 0.0d0) z(k) = z(k)/a(k,k)
         if (a(k,k) .eq. 0.0d0) z(k) = 1.0d0
         t = -z(k)
         call daxpy(k-1,t,a(1,k),1,z(1),1)
  160 continue
c     make znorm = 1.0
      s = 1.0d0/dasum(n,z,1)
      call dscal(n,s,z,1)
      ynorm = s*ynorm
c
      if (anorm .ne. 0.0d0) rcond = ynorm/anorm
      if (anorm .eq. 0.0d0) rcond = 0.0d0
      return
      end

      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(*),info
      double precision a(lda,*)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end

