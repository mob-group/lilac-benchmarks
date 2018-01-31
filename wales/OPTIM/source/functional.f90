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
SUBROUTINE FUNCTIONAL(Z,VNEW,POTEL,GTEST,STEST)
USE MODHESS
USE COMMONS
IMPLICIT NONE
LOGICAL GTEST, STEST
DOUBLE PRECISION Z(3*NATOMS), POTEL
DOUBLE PRECISION VNEW(3*NATOMS), G12, G11, DELTA, T

! G11=Z(1)
G12=Z(1)
DELTA=PARAM1
T=PARAM2
G11=PARAM3

! PRINT '(A,2G20.10)','G11,G12=',G11,G12
! PRINT '(A,2G20.10)','delta,t=',delta,t
! PRINT '(A,G20.10)','sqrt arg=',-((-2.0D0 + g11)*g11) - g12**2

IF (-((-2.0D0 + g11)*g11) - g12**2.LT.0.0D0) THEN
   PRINT '(A)','functional> ERROR *** negative square root argument'
   VNEW(1)=0.0D0
   HESS(1,1)=1.0D0
   RETURN
ENDIF

POTEL= (delta*(-2.0D0 + g11) + delta*g11 + (2.0D0 + 2.0D0*(-2 + g11)*g11 - g12**2*(-1.0D0 + Sqrt(-((-2.0D0 + g11)*g11) - g12**2)))/ &
  &            ((-1.0D0 + g11)**2 + g12**2) - 4.0D0*g12*t)/2.0D0

IF (GTEST) THEN
!  VNEW(1)= delta - ((-1.0D0 + g11)*g12**2*(-1.0D0 + (-2.0D0 + g11)*g11 + g12**2 - 2.0D0*Sqrt(-((-2.0D0 + g11)*g11) - g12**2)))/ &
!    &   (2.0D0*Sqrt(-((-2.0D0 + g11)*g11) - g12**2)*((-1.0D0 + g11)**2 + g12**2)**2)
!  VNEW(2)=((g12*(2.0D0 + (2.0D0*(-2.0D0 + g11)*g11 + 3.0D0*g12**2)/Sqrt(-((-2.0D0 + g11)*g11) - g12**2)))/ &
!    & ((-1.0D0 + g11)**2 + g12**2) + &
!    &    (-4.0D0*(-1.0D0 + g11)**2*g12 + 2.0D0*g12**3.0D0*(-1.0D0 + Sqrt(-((-2.0D0 + g11)*g11) - g12**2)))/ &
!    &     ((-1.0D0 + g11)**2 + g12**2)**2 - 4.0D0*t)/2.0D0
   VNEW(1)=((g12*(2.0D0 + (2.0D0*(-2.0D0 + g11)*g11 + 3.0D0*g12**2)/Sqrt(-((-2.0D0 + g11)*g11) - g12**2)))/ &
     & ((-1.0D0 + g11)**2 + g12**2) + &
     &    (-4.0D0*(-1.0D0 + g11)**2*g12 + 2.0D0*g12**3.0D0*(-1.0D0 + Sqrt(-((-2.0D0 + g11)*g11) - g12**2)))/ &
     &     ((-1.0D0 + g11)**2 + g12**2)**2 - 4.0D0*t)/2.0D0
ENDIF

IF (STEST) THEN
    HESS(1,1)=(-((g12**2*(-8*g11 + 4*g11**2 + 3*g12**2))/((2*g11 - g11**2 - g12**2)**1.5*(1 - 2*g11 + g11**2 + g12**2))) - &
     &    (2*g12**2*(2 + (2*(-2 + g11)*g11 + 3*g12**2)/Sqrt(-((-2 + g11)*g11) - g12**2)))/((-1 + g11)**2 + g12**2)**2 + &
     &    (2 + (2*(-2 + g11)*g11 + 3*g12**2)/Sqrt(-((-2 + g11)*g11) - g12**2))/((-1 + g11)**2 + g12**2) + &
     &    (-4*(-1 + g11)**2 - (2*g12**4)/Sqrt(-((-2 + g11)*g11) - g12**2) + 6*g12**2*(-1 + Sqrt(-((-2 + g11)*g11) - g12**2)))/((-1 + g11)**2 + g12**2)**2 - &
     &    (4*g12*(-4*(-1 + g11)**2*g12 + 2*g12**3*(-1 + Sqrt(-((-2 + g11)*g11) - g12**2))))/((-1 + g11)**2 + g12**2)**3)/2.
ENDIF

! PRINT'(A,4F20.10)','POTEL,grad,hess=',POTEL,VNEW(1),HESS(1,1)

RETURN

END SUBROUTINE FUNCTIONAL
