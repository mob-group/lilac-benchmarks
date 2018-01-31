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
      SUBROUTINE CTEST(N, VAR, V, ENERGY, GTEST, STEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION VAR(3*N), V(3*N), ENERGY, X, Y, R, PHI, E
      PARAMETER (E=2.718281828459045D0)
      LOGICAL CART, METRIC, EFFHESS, GTEST, STEST
      PARAMETER (CART=.TRUE.,METRIC=.FALSE.,EFFHESS=.FALSE.)

      IF (CART) THEN
         X=VAR(1)
         Y=VAR(2)
         ENERGY=(-2*x*y)/(E**(-E + Sqrt(x**2 + y**2))**2*(x**2 + y**2))
         IF (GTEST) THEN
            V(1)=(2*y*((-2*E*x**2)/(x**2 + y**2)**1.5 + (2*x**4 - y**2 + x**2*(1 + 2*y**2))/(x**2 + y**2)**2))/
     1  E**(E - Sqrt(x**2 + y**2))**2
            V(2)=(2*x*((-2*E*y**2)/(x**2 + y**2)**1.5 + (y**2 + 2*y**4 + x**2*(-1 + 2*y**2))/(x**2 + y**2)**2))/
     1  E**(E - Sqrt(x**2 + y**2))**2
         ENDIF
         IF (STEST) THEN
            HESS(1,1)=(4*x*y*((-4*x**2 + (3 - 2*(2 + E**2)*x**2)*(x**2 + y**2) + (3 - 2*x**2)*(x**2 + y**2)**2)/
     1       (x**2 + y**2)**3 + (E*(4*x**4 - 3*y**2 + x**2*(2 + 4*y**2)))/(x**2 + y**2)**2.5))/
     1       E**(E - Sqrt(x**2 + y**2))**2
            HESS(1,2)=(2*((2*E*(-y**4 + x**4*(-1 + 4*y**2) + x**2*y**2*(3 + 4*y**2)))/(x**2 + y**2)**2.5 + 
     1      (y**4 + 2*y**6 + x**6*(2 - 4*y**2) - 2*x**2*y**2*(3 + (1 + 2*E**2)*y**2 + 2*y**4) - 
     1         x**4*(-1 + (2 + 4*E**2)*y**2 + 8*y**4))/(x**2 + y**2)**3))/E**(E - Sqrt(x**2 + y**2))**2
            HESS(2,1)=HESS(1,2)
            HESS(2,2)=(4*x*y*(-((-3*(x**2 + x**4) + (1 + 2*x**2*(-1 + E**2 + x**2))*y**2 + 
     1           (1 + 2*E**2 + 4*x**2)*y**4 + 2*y**6)/(x**2 + y**2)**3) + 
     1      (E*(2*y**2 + 4*y**4 + x**2*(-3 + 4*y**2)))/(x**2 + y**2)**2.5))/
     1  E**(E - Sqrt(x**2 + y**2))**2
         ENDIF
      ELSE
         R=VAR(1)
         PHI=VAR(2)
         ENERGY=-(Sin(2*phi)/E**(-E + r)**2)
         IF (GTEST) THEN
            V(1)=(2*(-E + r)*Sin(2*phi))/EXP((E - r)**2)
            V(2)=(-2*Cos(2*phi))/EXP((E - r)**2)
         ENDIF
         IF (STEST) THEN
            IF (METRIC) THEN
               IF (EFFHESS) THEN
                  HESS(1,1)=(-2*(-1 + 2*E**2 + 2*r*(-2*E + r))*Sin(2*phi))/EXP((E - r)**2)
                  HESS(1,2)=(2*(1 + 2*r*(-E + r))*Cos(2*phi))/(EXP((E - r)**2)*r)
                  HESS(2,1)=(2*(1 + 2*r*(-E + r))*Cos(2*phi))/(EXP((E - r)**2)*r**3)
                  HESS(2,2)=(2*(2 - E*r + r**2)*Sin(2*phi))/(EXP((E - r)**2)*r**2)
C                 WRITE(*,'(A,F20.10)') 'A11=',HESS(1,1)
C                 WRITE(*,'(A,F20.10)') 'A12=',HESS(1,2)
C                 WRITE(*,'(A,F20.10)') 'A21=',HESS(2,1)
C                 WRITE(*,'(A,F20.10)') 'A21=',HESS(2,2)
               ELSE
                  HESS(1,1)=(-2*(-1 + 2*E**2 + 2*r*(-2*E + r))*Sin(2*phi))/E**(E - r)**2
                  HESS(1,2)=(4*(-E + r)*Cos(2*phi))/E**(E - r)**2
                  HESS(2,1)=(4*(-E + r)*Cos(2*phi))/(E**(E - r)**2*r**2)
                  HESS(2,2)=(4*Sin(2*phi))/(E**(E - r)**2*r**2)
               ENDIF
            ELSE
               IF (EFFHESS) THEN
                  HESS(1,1)=(-2*(-1 + 2*E**2 + 2*r*(-2*E + r))*Sin(2*phi))/E**(E - r)**2
                  HESS(1,2)=(4*(-E + r)*Cos(2*phi))/E**(E - r)**2-V(2)/R
                  HESS(2,1)=HESS(1,2)
                  HESS(2,2)=(4*Sin(2*phi))/E**(E - r)**2+R*V(1)
               ELSE
                  HESS(1,1)=(-2*(-1 + 2*E**2 + 2*r*(-2*E + r))*Sin(2*phi))/E**(E - r)**2
                  HESS(1,2)=(4*(-E + r)*Cos(2*phi))/E**(E - r)**2
                  HESS(2,1)=HESS(1,2)
                  HESS(2,2)=(4*Sin(2*phi))/E**(E - r)**2
               ENDIF
            ENDIF
         ENDIF
      ENDIF

      RETURN
      END
