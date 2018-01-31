!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of OPTIM.
!
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!   js850> This computes the curvature of the potential energy surface along the
!   direction vec at point coords.  This function is minimized by in xmylbfgs in
!   order to compute the lowest eigenvector.  This is known as Rayleigh-Ritz 
!   optimization because the curvature is defined via the Hessian matrix using
!   the Rayleigh-Ritz ratio
!
!           lambda(x, v) = v^T * H(x) * v / abs(v)**2
!
!   expanding the Hessian (H = d^2 E(x) / dx_i dx_j) in the derivatives gives
!   the approximation for the curvature used below.  In order to use LBFGS to
!   minimize lambda we need the derivative of lambda
!   with respect to v
!
!           d lambda(x,v) / dv = 2* v * H - 2 * lambda * v
!
!   Again expanding the hessian in powers alows us to compute this without
!   knowing the actual hessian
!
!
      SUBROUTINE SECDIAG(VEC,COORDS,ENERGY,grad,GL,DIAG,GTEST,XRMS)
      USE COMMONS
      USE KEY
      USE MODCHARMM
      USE PORFUNCS
      IMPLICIT NONE
      INTEGER J1
      LOGICAL GTEST, FPLUS, FMINUS
      DOUBLE PRECISION, INTENT(INOUT) :: COORDS(3*NATOMS) ! the point at which to compute the curvature (this should be intent(IN), but I couldn't easily guarantee that orthogopt didn't modify it)
      DOUBLE PRECISION, INTENT(IN) :: VEC(3*NATOMS) ! the direction along which to compute the curvature (the current best estimate for the lowest eigenvector)
      DOUBLE PRECISION, INTENT(IN) :: ENERGY ! the energy at point coords
      DOUBLE PRECISION, INTENT(IN) :: GRAD(3*NATOMS) ! the gradient at coords (will be nonsense unless nsecdiag > 2)
      DOUBLE PRECISION, INTENT(OUT) :: DIAG ! the curvature
      DOUBLE PRECISION, INTENT(OUT) :: GL(3*NATOMS) ! the gradient of the curvature
      DOUBLE PRECISION, INTENT(OUT) :: XRMS ! the rms of GL
      DOUBLE PRECISION DUMMY3(3*NATOMS), DIFF, DIAG2, DIAG3, &
                       EPLUS,EMINUS,GRAD1(3*NATOMS),LOCALV(3*NATOMS), &
                       RMS,GRAD2(3*NATOMS), VECL, ZETA, PROJ
      double precision diag4 ! the curvature computed with the first order forward finite differences method

      DIAG = 0.0D0
      DIAG2 = 0.0D0
      DIAG3 = 0.0D0
!      EPLUS = 0.0D0
      EMINUS = 0.0D0
!      XRMS = 0.0D0
      DIFF=1.0D-3
      IF (NIH2LEPST .or. NIMET .or. NIHEAM7T .or. NIHLEPST .or. NIHPAIRONLYT) DIFF=1.0D-5
      IF (PYADD2T) DIFF=1.0D-5
      IF (CHARMMDFTBT) DIFF=1.0D-2
      IF (AMHT) DIFF=1.0D-2
      IF (ZSYM(NATOMS).EQ.'GO') DIFF=2.0D-3
      IF (GAMESSUK.OR.GAMESSUS) DIFF=1.0D-2
      IF (CADPAC) DIFF=1.0D-2
      IF (CASTEP) DIFF=0.01D0
      IF (GAUSSIAN03.OR.GAUSSIAN09.OR.GAUSSIAN16) DIFF=1.0D-2
      IF (QCHEM) DIFF=0.01D0
      IF (VASP) DIFF=0.01D0
      IF (ONETEP) DIFF=0.01D0
      IF (CP2K) DIFF=0.001D0 
      IF (CPMD) DIFF=0.04D0
!     IF (CHRMMT) DIFF=5.0D-2
      IF (CHRMMT) DIFF=0.01D0
!     IF (DFTBT) DIFF=1.0D-3
!
!  Must read VEC into LOCALV because we are going to play with the vector in
!  question. This would mess up cases where we need to retry with a smaller
!  step size, because we cannot undo the previous step cleanly if the magnitude
!  of VEC is changed! 1/7/04 DJW
!

      LOCALV(1:NOPT)=VEC(1:NOPT)
!     PRINT '(A)','secdiag> vec before orthogopt'
!     PRINT '(6G20.10)',LOCALV(1:NOPT)
      IF (NFREEZE.LT.3.AND..NOT.MIEFT .AND. (.NOT. NOTRANSROTT)) THEN
        CALL ORTHOGOPT(LOCALV,COORDS,.TRUE.)
      ELSE
         CALL VECNORM(LOCALV,NOPT)
      ENDIF
!     PRINT '(A)','secdiag> vec after orthogopt'
!     PRINT '(6G20.10)',LOCALV(1:NOPT)

      IF (FREEZE) THEN
         DO J1=1,NATOMS
            IF (FROZEN(J1)) THEN
               IF (VARIABLES) THEN
                  LOCALV(J1)=0.0D0
               ELSE
                  LOCALV(3*(J1-1)+1)=0.0D0
                  LOCALV(3*(J1-1)+2)=0.0D0
                  LOCALV(3*(J1-1)+3)=0.0D0
               ENDIF
            ENDIF
         ENDDO
      ENDIF

      VECL=1.0D0
      ZETA=DIFF

      ! compute the energy and gradient at coords + zeta*vec
      DO J1=1,NOPT
         DUMMY3(J1)=COORDS(J1)+ZETA*LOCALV(J1)
      ENDDO

      IF (CPMD) THEN
         INQUIRE(FILE='RESTART.1.plus',EXIST=FPLUS)
         IF (FPLUS) THEN
            CALL SYSTEM(' cp RESTART.1.plus RESTART.1 ')
         ELSE
            CALL SYSTEM(' cp RESTART.1.save RESTART.1 ')
         ENDIF
      ELSE IF (CASTEP) THEN
!        INQUIRE(FILE=SYS(1:LSYS) // '.wvfn.plus',EXIST=FPLUS)
!        IF (FPLUS) THEN
!           CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.wvfn.plus ' // SYS(1:LSYS) // '.wvfn.1 ')
!        ELSE
!           CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.wvfn.save ' // SYS(1:LSYS) // '.wvfn.1 ')
!        ENDIF
      ENDIF
!     PRINT*,'DUMMY3:'
      CALL POTENTIAL(DUMMY3,EPLUS,GRAD1,GTEST,.FALSE.,RMS,.FALSE.,.FALSE.)

      IF (CPMD) CALL SYSTEM(' cp RESTART.1 RESTART.1.plus ')
!     IF (CASTEP) CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.wvfn.1 ' // SYS(1:LSYS) // '.wvfn.plus ')

      IF (NSECDIAG .LE. 2) THEN
         ! use the second order central central differences expasion
         ! this requires an extra potential call GMINUS
         DO J1=1,NOPT
            DUMMY3(J1)=COORDS(J1)-ZETA*LOCALV(J1)
         ENDDO
         IF (CPMD) THEN
            INQUIRE(FILE='RESTART.1.minus',EXIST=FMINUS)
            IF (FMINUS) THEN
               CALL SYSTEM(' cp RESTART.1.minus RESTART.1 ')
            ELSE
               CALL SYSTEM(' cp RESTART.1.save RESTART.1 ')
            ENDIF
         ELSE IF (CASTEP) THEN
!           INQUIRE(FILE=SYS(1:LSYS) // '.wvfn.minus',EXIST=FMINUS)
!           IF (FMINUS) THEN
!              CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.wvfn.minus ' // SYS(1:LSYS) // '.wvfn.1 ')
!           ELSE
!              CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.wvfn.save ' // SYS(1:LSYS) // '.wvfn.1 ')
!           ENDIF
         ENDIF
!        PRINT*,'DUMMY3:'
         CALL POTENTIAL(DUMMY3,EMINUS,GRAD2,GTEST,.FALSE.,RMS,.FALSE.,.FALSE.)
         IF (CPMD) CALL SYSTEM(' cp RESTART.1 RESTART.1.minus ')
!        IF (CASTEP) CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.wvfn.1 ' // SYS(1:LSYS) // '.wvfn.minus ')
         ! js850> diag is a second order central differences expansion of the curvature using energies
         DIAG=(EPLUS+EMINUS-2.0D0*ENERGY)/(ZETA**2*VECL)
!        IF (DEBUG) WRITE(*,'(A,G15.5)') 'secdiag> DIFF',ZETA
         ! js850> diag2 is a second order central differences expansion of the curvature using gradients
         ! diag and diag2 are of the same order, but diag2 probably has fewer numerical precision issues
         DIAG2=0.0D0
         DO J1=1,NOPT
            DIAG2=DIAG2+(GRAD1(J1)-GRAD2(J1))*LOCALV(J1)
!           WRITE(*,'(A,I4,4F20.10)') 'J1,GRAD1,GRAD2,LOCALV,DIAG2=',J1,GRAD1(J1),GRAD2(J1),LOCALV(J1),DIAG2
         ENDDO
         DIAG2=DIAG2/(2.0D0*ZETA)
         DIAG3=2*(DIAG-DIAG2/2)
!     IF (.NOT.GTEST) WRITE(*,'(A,6F20.10)') 'D,D2,D3,E+,E-,E=',DIAG,DIAG2,DIAG3,EPLUS,EMINUS,ENERGY
!
!  Although DIAG3 is a more accurate estimate of the diagonal second derivative, it
!  cannot be differentiated analytically.
!
      ELSE ! (nsecdiag .gt. 2)
         ! here we use a lower order forward finite differences expansion where
         ! we don't need GMINUS.  We do need GRAD though
         DIAG4 = sum((grad1(:) - grad(:)) * localv(:)) / zeta
      ENDIF
      IF (GTEST) THEN
         DO J1=1,NOPT
            IF (NSECDIAG.gt.2) THEN
               ! first order forward finite differences
               GL(J1)=2.d0 * (GRAD1(J1)-GRAD(J1))/(ZETA*VECL**2)-2.0D0*DIAG4*LOCALV(J1)/VECL**2
            ELSEIF (NSECDIAG.EQ.2) THEN
               ! second order central finite differences
               GL(J1)=(GRAD1(J1)-GRAD2(J1))/(ZETA*VECL**2)-2.0D0*DIAG2*LOCALV(J1)/VECL**2
            ELSE
               ! second order central finite differences
               GL(J1)=(GRAD1(J1)-GRAD2(J1))/(ZETA*VECL**2)-2.0D0*DIAG*LOCALV(J1)/VECL**2
            ENDIF
!           WRITE(*,'(A,I4,4G16.7)') 'secdiag> J1,GRAD1,GRAD2,LOCALV,GL=',J1,GRAD1(J1),GRAD2(J1),LOCALV(J1),GL(J1)
         ENDDO
         IF (NFREEZE.LT.3.AND..NOT.MIEFT .AND. (.NOT. NOTRANSROTT)) CALL ORTHOGOPT(GL,COORDS,.FALSE.)
!        PRINT *,'secdiag> before proj stuff GL:'
!        PRINT '(3F20.10)',GL(1:NOPT)
!        CALL ORTHOGOPT(GL,COORDS,.FALSE.) ! seems to do some good for MSEVB
!
!  Project out any component of the gradient along LOCALV (which is a unit vector)
!  This is a big improvement for DFTB.  js850> The parallel components have already
!  been projected out, this simply improves the numerical precision
!
         PROJ=0.0D0
         DO J1=1,NOPT
            PROJ=PROJ+GL(J1)*LOCALV(J1)
         ENDDO
         DO J1=1,NOPT
            GL(J1)=GL(J1)-PROJ*LOCALV(J1)
         ENDDO
         XRMS=0.0D0
         DO J1=1,NOPT
            XRMS=XRMS+GL(J1)**2
         ENDDO
         XRMS=DSQRT(XRMS/NOPT)
!        PRINT *,'secdiag> after proj stuff GL:'
!        PRINT '(3F20.10)',GL(1:NOPT)
         IF (DEBUG) THEN
            IF(NSECDIAG.NE.3) THEN
               WRITE(*,'(A,3G15.5,3G20.12,G10.3)') 'D,D2,D3,E+,E-,E,RMS=',DIAG,DIAG2,DIAG3,EPLUS,EMINUS,ENERGY,XRMS
               WRITE(*,'(A,G20.10)') 'predicted gradient component=',(EPLUS-EMINUS)/(2*ZETA)
            ELSE  ! forward difference method: don't need to print so much
               WRITE(*,'(A,3G15.5,3G20.12,G10.3)') 'D4,E+,E,RMS=',DIAG4,EPLUS,ENERGY,XRMS
            ENDIF
         ENDIF
      ENDIF
!     PRINT '(A)','LOCALV:'
!     PRINT '(3G20.10)',LOCALV(1:NOPT)
      IF (NSECDIAG.EQ.2) DIAG=DIAG2
      IF (NSECDIAG.GT.2) DIAG=DIAG4

      RETURN
      END
