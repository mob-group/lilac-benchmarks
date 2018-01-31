!
!
!
 SUBROUTINE Poly2(nAtoms, vCoords, vGrad, energy, gTest, sTest)
! **********************************************************************
! *                                                                    *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       08/05/2014                                             *
! * Version:    1.0                                                    *
! *                                                                    *
! **********************************************************************
!
   USE KEY, ONLY: HARMPOLYT, HARMPOLY_BONLEN, HARMPOLY_K
   USE MODHESS, ONLY: Hess
   IMPLICIT NONE
!
! -------------------------------------------------------------------
!
!
!   number of polymer beads
!
   INTEGER, INTENT(IN) ::  nAtoms
!
!   coordinate array (x1, y1, x2, y2, ... znAtoms)
!
   DOUBLE PRECISION, INTENT(IN), DIMENSION(3*nAtoms) :: vCoords 
!
!   switches - gradient and hessian
!
   LOGICAL, INTENT(IN) :: gTest, sTest
!
! -------------------------------------------------------------------
!
   DOUBLE PRECISION :: energy
!
!   energy gradient
!
   DOUBLE PRECISION, DIMENSION(3*nAtoms) :: vGrad ! IN and OUT
!
! -------------------------------------------------------------------
!
   INTEGER :: i, j, i1, i2
   DOUBLE PRECISION :: k, r0, r1, alpha, beta, x12, y12, z12
   DOUBLE PRECISION :: eCon, rp, prefac, r2
!
! -------------------------------------------------------------------
!
   WRITE(*,*) "poly2 potential called, lj energy", energy
   IF (.NOT.HARMPOLYT) THEN 
      WRITE(*,*) "error> HARMPOLYTnot set up"
      STOP
   ENDIF
!
!  bonded contributions
!
!
!   parameters of potential
!
   k = HARMPOLY_K
   r0 = HARMPOLY_BONLEN
!   r0 = 2.0**(1.0/6.0)
!
!  0.5*k*(|r|-r0)^2
!
   cycBond: DO i = 2,nAtoms
!
!   for each bond (there are nAtoms-1 of them)
!
      r1 = MAX(SQRT(SUM((vCoords(3*i-5:3*i-3)-vCoords(3*i-2:3*i))**2)),1.0d-32)
      alpha = k*(r1-r0)/r1
      beta = k*r0/(r1**3)
!
!   energy
!
      energy = energy + 0.5*k*(r1-r0)**2
      x12 = vCoords(3*i-2)-vCoords(3*i-5)
      y12 = vCoords(3*i-1)-vCoords(3*i-4)
      z12 = vCoords(3*i)  -vCoords(3*i-3)
      vGrad(3*i-5) = vGrad(3*i-5) - alpha*x12
      vGrad(3*i-4) = vGrad(3*i-4) - alpha*y12
      vGrad(3*i-3) = vGrad(3*i-3) - alpha*z12
      vGrad(3*i-2) = vGrad(3*i-2) + alpha*x12
      vGrad(3*i-1) = vGrad(3*i-1) + alpha*y12
      vGrad(3*i)   = vGrad(3*i)   + alpha*z12
!
!   contribution to Hessian
!
      IF (sTest) THEN
         Hess(3*i-5,3*i-5) = Hess(3*i-5,3*i-5) + alpha + beta*x12*x12
         Hess(3*i-5,3*i-4) = Hess(3*i-5,3*i-4)         + beta*x12*y12
         Hess(3*i-5,3*i-3) = Hess(3*i-5,3*i-3)         + beta*x12*z12
         Hess(3*i-5,3*i-2) = Hess(3*i-5,3*i-2) - alpha - beta*x12*x12
         Hess(3*i-5,3*i-1) = Hess(3*i-5,3*i-1)         - beta*x12*y12
         Hess(3*i-5,3*i)   = Hess(3*i-5,3*i)           - beta*x12*z12
         
         Hess(3*i-4,3*i-5) = Hess(3*i-4,3*i-5)         + beta*y12*x12
         Hess(3*i-4,3*i-4) = Hess(3*i-4,3*i-4) + alpha + beta*y12*y12
         Hess(3*i-4,3*i-3) = Hess(3*i-4,3*i-3)         + beta*y12*z12
         Hess(3*i-4,3*i-2) = Hess(3*i-4,3*i-2)         - beta*y12*x12
         Hess(3*i-4,3*i-1) = Hess(3*i-4,3*i-1) - alpha - beta*y12*y12
         Hess(3*i-4,3*i)   = Hess(3*i-4,3*i)           - beta*y12*z12
         
         Hess(3*i-3,3*i-5) = Hess(3*i-3,3*i-5)         + beta*z12*x12
         Hess(3*i-3,3*i-4) = Hess(3*i-3,3*i-4)         + beta*z12*y12
         Hess(3*i-3,3*i-3) = Hess(3*i-3,3*i-3) + alpha + beta*z12*z12
         Hess(3*i-3,3*i-2) = Hess(3*i-3,3*i-2)         - beta*z12*x12
         Hess(3*i-3,3*i-1) = Hess(3*i-3,3*i-1)         - beta*z12*y12
         Hess(3*i-3,3*i)   = Hess(3*i-3,3*i)   - alpha - beta*z12*z12
         
         Hess(3*i-2,3*i-5) = Hess(3*i-2,3*i-5) - alpha - beta*x12*x12
         Hess(3*i-2,3*i-4) = Hess(3*i-2,3*i-4)         - beta*x12*y12
         Hess(3*i-2,3*i-3) = Hess(3*i-2,3*i-3)         - beta*x12*z12
         Hess(3*i-2,3*i-2) = Hess(3*i-2,3*i-2) + alpha + beta*x12*x12
         Hess(3*i-2,3*i-1) = Hess(3*i-2,3*i-1)         + beta*x12*y12
         Hess(3*i-2,3*i)   = Hess(3*i-2,3*i)           + beta*x12*z12
         
         Hess(3*i-1,3*i-5) = Hess(3*i-1,3*i-5)         - beta*y12*x12
         Hess(3*i-1,3*i-4) = Hess(3*i-1,3*i-4) - alpha - beta*y12*y12
         Hess(3*i-1,3*i-3) = Hess(3*i-1,3*i-3)         - beta*y12*z12
         Hess(3*i-1,3*i-2) = Hess(3*i-1,3*i-2)         + beta*y12*x12
         Hess(3*i-1,3*i-1) = Hess(3*i-1,3*i-1) + alpha + beta*y12*y12
         Hess(3*i-1,3*i)   = Hess(3*i-1,3*i)           + beta*y12*z12
         
         Hess(3*i,  3*i-5) = Hess(3*i,  3*i-5)         - beta*z12*x12
         Hess(3*i,  3*i-4) = Hess(3*i,  3*i-4)         - beta*z12*y12
         Hess(3*i,  3*i-3) = Hess(3*i,  3*i-3) - alpha - beta*z12*z12
         Hess(3*i,  3*i-2) = Hess(3*i,  3*i-2)         + beta*z12*x12
         Hess(3*i,  3*i-1) = Hess(3*i,  3*i-1)         + beta*z12*y12
         Hess(3*i,  3*i)   = Hess(3*i,  3*i)   + alpha + beta*z12*z12
      ENDIF
   ENDDO cycBond
!
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE Poly2
!
!
!
 SUBROUTINE PORE8(nAtoms, vCoords, vGrad, energy, gTest, sTest)
! **********************************************************************
! *                                                                    *
! *                                                                    *
! **********************************************************************
! *                                                                    *
! * Author:     Boris Fackovec                                         *
! * Date:       08/05/2015                                             *
! * Version:    1.0                                                    *
! *                                                                    *
! **********************************************************************
!
   USE KEY, ONLY: PORE8T, PORE8_AXIS, PORE8_ENERGY
   USE MODHESS, ONLY: Hess
   IMPLICIT NONE
!
! -------------------------------------------------------------------
!
!
!   number of polymer beads
!
   INTEGER, INTENT(IN) ::  nAtoms
!
!   coordinate array (x1, y1, x2, y2, ... znAtoms)
!
   DOUBLE PRECISION, INTENT(IN), DIMENSION(3*nAtoms) :: vCoords 
!
!   switches - gradient and hessian
!
   LOGICAL, INTENT(IN) :: gTest, sTest
!
! -------------------------------------------------------------------
!
   DOUBLE PRECISION :: energy
!
!   energy gradient
!
   DOUBLE PRECISION, DIMENSION(3*nAtoms) :: vGrad 
!
! -------------------------------------------------------------------
!
   INTEGER :: i, j, i1, i2
   DOUBLE PRECISION :: prefac, r2
!
! -------------------------------------------------------------------
!
   WRITE(*,*) "pore8 potential called, lj energy", energy
   prefac = PORE8_ENERGY ! prefactor
   IF (PORE8_AXIS.EQ.1) THEN
      i1 = 2
      i2 = 3
   ELSEIF (PORE8_AXIS.EQ.2) THEN
      i1 = 1
      i2 = 2
   ELSEIF (PORE8_AXIS.EQ.3) THEN
      i1 = 1
      i2 = 3
   ELSE
      WRITE(*,*) "error> PORE8_AXIS other than {1,2,3}"
   ENDIF
   cycPore8: DO i=1,nAtoms
      r2 = (vCoords(3*i-2)**2 + vCoords(3*i-1)**2)
      energy = energy + prefac*(r2**4)
      vGrad(3*i-3+i1) = vGrad(3*i-3+i1) + 8*prefac*(r2**3)*vCoords(3*i-3+i1)
      vGrad(3*i-3+i2) = vGrad(3*i-3+i2) + 8*prefac*(r2**3)*vCoords(3*i-3+i2)
!      vGrad(3*i-2) += vGrad(3*i-2) + 8*prefac*(r2**3)*vCoords(3*i-2)
!      vGrad(3*i-1) += vGrad(3*i-1) + 8*prefac*(r2**3)*vCoords(3*i-1)
      IF(sTest) THEN
         Hess(3*i-3+i1,3*i-3+i1) = Hess(3*i-3+i1,3*i-3+i1) + &
              8*prefac*(r2**3) + 48*prefac*(r2**2)*(vCoords(3*i-3+i1)**2)
         Hess(3*i-3+i1,3*i-3+i2) = Hess(3*i-3+i1,3*i-3+i2) + &
              48*prefac*(r2**2)*vCoords(3*i-3+i1)*vCoords(3*i-3+i2)
         Hess(3*i-3+i2,3*i-3+i1) = Hess(3*i-3+i2,3*i-3+i1) + &
              48*prefac*(r2**2)*vCoords(3*i-2)*vCoords(3*i-1)
         Hess(3*i-3+i2,3*i-3+i2) = Hess(3*i-3+i2,3*i-3+i2) + &
              8*prefac*(r2**3) + 48*prefac*(r2**2)*(vCoords(3*i-3+i2)**2)
!         Hess(3*i-2,3*i-2) = Hess(3*i-2,3*i-2) + 8*prefac*(r2**3) + &
!              48*prefac*(r2**2)*(vCoords(3*i-2)**2)
!         Hess(3*i-2,3*i-1) = Hess(3*i-2,3*i-1) + &
!              48*prefac*(r2**2)*vCoords(3*i-2)*vCoords(3*i-1)
!         Hess(3*i-1,3*i-2) = Hess(3*i-1,3*i-2) + &
!              48*prefac*(r2**2)*vCoords(3*i-2)*vCoords(3*i-1)
!         Hess(3*i-1,3*i-1) = Hess(3*i-1,3*i-1) + 8*prefac*(r2**3) + &
!              48*prefac*(r2**2)*(vCoords(3*i-1)**2)
      ENDIF
   ENDDO cycPore8
!   
!   WRITE(*,*) vCoords
!   WRITE(*,*) vGrad
   RETURN
!
! -------------------------------------------------------------------
!
 END SUBROUTINE PORE8
!
