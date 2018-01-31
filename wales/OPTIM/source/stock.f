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
!
!  Energy and gradient for the Stockmayer potential using two polar
!  coordinates for the dipole direction
!
      SUBROUTINE STOCK(NATOMS,X,V,ESTOCK,GTEST,STEST)
      USE KEY,ONLY : STOCKMU, STOCKLAMBDA, STOCKSPIN, STOCKZTOL, STOCKMAXSPIN
      USE MODHESS
      IMPLICIT NONE
      LOGICAL GTEST,STEST
      LOGICAL ZALIGNTEST
      INTEGER J1, J2, K1, K2, NAT2, J3, J4, REALNATOMS, NATOMS, OFFSET, SPINITER
      DOUBLE PRECISION X(3*NATOMS), V(3*NATOMS), ESTOCK
      DOUBLE PRECISION DPRAND, S, SIGMA1, SIGMA2, THETA1, THETA2, Q(4), PI, M(3,3)
      DOUBLE PRECISION MU2, DUMMY
      DOUBLE PRECISION X1, X2, Y1, Y2, Z1, Z2, T1, P1, T2, P2
      DOUBLE PRECISION R12, R122, R125, R126, R127, R128, R129, R1211, R1212, R1214, R1216
      DOUBLE PRECISION N1DOTN2, N1DOTR12, N2DOTR12, CP1P2, SP1PP2, SP1MP2
      DOUBLE PRECISION CT1, CT2, ST1, ST2, CP1, CP2, SP1, SP2, DX, DY, DZ, DX2, DY2, DZ2

      REALNATOMS=NATOMS/2
      OFFSET = 3*REALNATOMS
      MU2=STOCKMU**2
      ESTOCK=0.0D0
      V(1:6*REALNATOMS)=0.0D0

      IF (STOCKSPIN) THEN
         PI = ATAN2(1.0D0, 1.0D0) * 4.0D0
!        Check that no spins are too closely aligned with z, since that would
!        introduce redundant phi angles.
         SPINITER = 1
         DO
            IF (.NOT.ZALIGNTEST(NATOMS, X)) THEN
               EXIT
            END IF
            IF (SPINITER.GT.STOCKMAXSPIN) THEN
               PRINT*,'WARNING: Randomisation of orientation failed to remove z-aligned dipoles'
               EXIT
            END IF

!           Random quaternion
            S = DPRAND();
            SIGMA1 = SQRT(1.0D0-S);
            SIGMA2 = SQRT(S);
            THETA1 = 2.0D0*PI*DPRAND();
            THETA2 = 2.0*PI*DPRAND();
            Q(1) = COS(THETA2) * SIGMA2;
            Q(2) = SIN(THETA1) * SIGMA1;
            Q(3) = COS(THETA1) * SIGMA1;
            Q(4) = SIN(THETA2) * SIGMA2;

!           Rotation matrix corresponding to the quaternion
            M(1,1) = Q(1)*Q(1) + Q(2)*Q(2) - Q(3)*Q(3) - Q(4)*Q(4)
            M(2,1) = 2.0*(Q(2)*Q(3) + Q(1)*Q(4))
            M(3,1) = 2.0*(Q(2)*Q(4) - Q(1)*Q(3))
            M(1,2) = 2.0*(Q(2)*Q(3) - Q(1)*Q(4))
            M(2,2) = Q(1)*Q(1) - Q(2)*Q(2) + Q(3)*Q(3) - Q(4)*Q(4)
            M(3,2) = 2.0*(Q(3)*Q(4) + Q(1)*Q(2))
            M(1,3) = 2.0*(Q(2)*Q(4) + Q(1)*Q(3))
            M(2,3) = 2.0*(Q(3)*Q(4) - Q(1)*Q(2))
            M(3,3) = Q(1)*Q(1) - Q(2)*Q(2) - Q(3)*Q(3) + Q(4)*Q(4)

            CALL NEWROTGEOMSTOCK(NATOMS, X, M, 0.0D0, 0.0D0, 0.0D0)

            SPINITER = SPINITER + 1
         END DO

         IF ( (SPINITER.GT.1).AND.(SPINITER.LE.STOCKMAXSPIN) ) THEN
            PRINT*,'WARNING: Orientation has been randomised to remove z-aligned dipoles'
         END IF
      END IF
      
      DO J1=1,REALNATOMS
         J3=3*J1
         X1=X(J3-2)
         Y1=X(J3-1)
         Z1=X(J3)
         T1=X(OFFSET+J3-2)
         CT1=cos(T1)
         ST1=sin(T1)
         P1=X(OFFSET+J3-1)
         CP1=cos(P1)
         SP1=sin(P1)
         DO J2=J1+1,REALNATOMS
            J4=3*J2
            X2=X(J4-2)
            Y2=X(J4-1)
            Z2=X(J4)
            T2=X(OFFSET+J4-2)
            CT2=cos(T2)
            ST2=sin(T2)
            P2=X(OFFSET+J4-1)
            CP2=cos(P2)
            SP2=sin(P2)
            n1dotr12=st1*cp1*(x1-x2)+st1*sp1*(y1-y2)+ct1*(z1-z2)
            n2dotr12=st2*cp2*(x1-x2)+st2*sp2*(y1-y2)+ct2*(z1-z2)
            n1dotn2=st1*cp1*st2*cp2 + st1*sp1*st2*sp2 + ct1*ct2
            R122=(X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2
            R12=SQRT(R122)
            R126=R122**3
            R127=R126*R12
            R129=R127*R122
            R1212=R126**2
            R125=R126/R12
            R1214=R1212*R122
           
            DUMMY= (MU2*R127*(-3*n1dotr12*n2dotr12 + n1dotn2*R122) +
     &             (4 - 4*stocklambda*R126))/R1212
         
            ESTOCK=ESTOCK+DUMMY

! derivatives for positions

      DUMMY = (-3*(-((-16 + 8*stocklambda*R126 + 
     &   mu2*(5*n1dotr12*n2dotr12*R127 - n1dotn2*R129))*(x1 - x2))
     &   + mu2*n2dotr12*R129*CP1*ST1 + mu2*n1dotr12*R129*CP2*ST2))/R1214
      V(J3-2)=V(J3-2)+DUMMY
      V(J4-2)=V(J4-2)-DUMMY

      DUMMY = (-3*(-((-16 + 8*stocklambda*R126 + 
     &   mu2*(5*n1dotr12*n2dotr12*R127 - n1dotn2*R129))*(y1 - y2))
     &   + mu2*n2dotr12*R129*SP1*ST1 + mu2*n1dotr12*R129*SP2*ST2))/R1214
      V(J3-1)=V(J3-1)+DUMMY
      V(J4-1)=V(J4-1)-DUMMY

      DUMMY = (-3*(-((-16 + 8*stocklambda*R126 + 
     &   mu2*(5*n1dotr12*n2dotr12*R127 - n1dotn2*R129))*(z1 - z2))
     &   + mu2*n2dotr12*R129*CT1 + mu2*n1dotr12*R129*CT2))/R1214
      V(J3)=V(J3)+DUMMY
      V(J4)=V(J4)-DUMMY

! derivatives for angular variables of atom J1

      V(OFFSET+J3-2) = V(OFFSET+J3-2) +
     &   (mu2*((3*n2dotr12*z1 - 3*n2dotr12*z2 - R122*CT2)*ST1 + 
     &   CP1*CT1*(3*n2dotr12*(-x1 + x2) + R122*CP2*ST2) + 
     &   CT1*SP1*(-3*n2dotr12*y1 + 3*n2dotr12*y2 + R122*SP2*ST2)))/R125

      V(OFFSET+J3-1) = V(OFFSET+J3-1) +
     &   (mu2*ST1*(SP1*(3*n2dotr12*(x1 - x2) - R122*CP2*ST2) + 
     &   CP1*(3*n2dotr12*(-y1 + y2) + R122*SP2*ST2)))/R125

! derivatives for angular variables of atom J2

      V(OFFSET+J4-2) = V(OFFSET+J4-2) +
     &   (mu2*(CP2*CT2*(3*n1dotr12*(-x1 + x2) + R122*CP1*ST1) + 
     &   CT2*SP2*(-3*n1dotr12*y1 + 3*n1dotr12*y2 + R122*SP1*ST1) + 
     &   (3*n1dotr12*z1 - 3*n1dotr12*z2 - R122*CT1)*ST2))/R125

      V(OFFSET+J4-1) = V(OFFSET+J4-1) +
     &   (mu2*(SP2*(3*n1dotr12*(x1 - x2) - R122*CP1*ST1) + 
     &   CP2*(3*n1dotr12*(-y1 + y2) + R122*SP1*ST1))*ST2)/R125

         ENDDO
      ENDDO

!     ***********************
!     Analytic Hessian matrix
!     ***********************
      IF (STEST) THEN
         HESS(1:3*NATOMS, 1:3*NATOMS) = 0.0D0
         DO J1=1, REALNATOMS
            J3 = J1*3
!           Coordinates and functions of first atom
            x1 = X(J3 - 2)
            y1 = X(J3 - 1)
            z1 = X(J3)
            t1 = X(OFFSET + J3 - 2)
            p1 = X(OFFSET + J3 - 1)
            st1 = SIN(t1)
            ct1 = COS(t1)
            sp1 = SIN(p1)
            cp1 = COS(p1)

            DO J2=1, REALNATOMS
               IF (J1 == J2) CYCLE
               J4 = J2*3

!              Coordinates and functions of second atom
               x2 = X(J4 - 2)
               y2 = X(J4 - 1)
               z2 = X(J4)
               t2 = X(OFFSET + J4 - 2)
               p2 = X(OFFSET + J4 - 1)
               st2 = SIN(t2)
               ct2 = COS(t2)
               sp2 = SIN(p2)
               cp2 = COS(p2)

!              Coordinates and functions related to both atoms
               dx = x1 - x2
               dy = y1 - y2
               dz = z1 - z2
               dx2 = dx * dx
               dy2 = dy * dy
               dz2 = dz * dz
               n1dotr12 = st1*cp1*dx + st1*sp1*dy + ct1*dz
               n2dotr12 = st2*cp2*dx + st2*sp2*dy + ct2*dz
               n1dotn2 = st1*cp1*st2*cp2 + st1*sp1*st2*sp2 + ct1*ct2
               r122 = dx2 + dy2 + dz2
               r12 = SQRT(r122)
               r125 = r122*r122*r12
               r126 = r125*r12
               r127 = r126*r12
               r128 = r126*r122
               r129 = r127*r122
               r1211 = r126*r125
               r1212 = r126*r126
               r1214 = r1212*r122
               r1216 = r128*r128
               sp1pp2 = SIN(p1 + p2)
               sp1mp2 = SIN(p1 - p2)
               cp1p2 = COS(p1 - p2)

!              [1] The five completely diagonal terms: same atom, same coordinate
!              x1,x1
               HESS(J3-2, J3-2) = HESS(J3-2, J3-2) +
     &            (-3*(MU2*n1dotn2*R1211 + 16*R122 - 8*stocklambda*R128 +
     &            dx2*(-224 + 64*stocklambda*R126 + 35*MU2*n1dotr12*n2dotr12*R127 -
     &            5*MU2*n1dotn2*R129) -
     &            MU2*R129*(5*n1dotr12*(n2dotr12 + 2*cp2*dx*st2) +
     &            2*cp1*st1*(5*dx*n2dotr12 - cp2*R122*st2))))/R1216
!              y1,y1
               HESS(J3-1, J3-1) = HESS(J3-1, J3-1) +
     &            (3*(-(MU2*n1dotn2*R1211) - 16*R122 + 8*stocklambda*R128 +
     &            dy2*(224 - 64*stocklambda*R126 +
     &            5*MU2*(-7*n1dotr12*n2dotr12*R127 + n1dotn2*R129)) +
     &            MU2*R129*(5*n1dotr12*(n2dotr12 + 2*dy*sp2*st2) +
     &            2*sp1*st1*(5*dy*n2dotr12 - R122*sp2*st2))))/R1216
!              z1,z1
               HESS(J3, J3) = HESS(J3, J3) +
     &            (-3*(MU2*n1dotn2*R1211 + 16*R122 +
     &            dz2*(-224 + 64*stocklambda*R126 + 35*MU2*n1dotr12*n2dotr12*R127) -
     &            8*stocklambda*R128) + 3*MU2*(5*
     &            (dz2*n1dotn2 + 2*ct2*dz*n1dotr12 + (2*ct1*dz + n1dotr12)*n2dotr12)
     &            - 2*ct1*ct2*R122)*R129)/R1216
!              t1,t1
               HESS(OFFSET+J3-2, OFFSET+J3-2) = HESS(OFFSET+J3-2, OFFSET+J3-2) +
     &            (MU2*(3*ct1*dz*n2dotr12 - ct1*ct2*R122 + 3*cp1*dx*n2dotr12*st1 +
     &            3*dy*n2dotr12*sp1*st1 - cp1p2*R122*st1*st2))/R125
!              p1,p1
               HESS(OFFSET+J3-1, OFFSET+J3-1) = HESS(OFFSET+J3-1, OFFSET+J3-1) +
     &            (MU2*st1*(3*cp1*dx*n2dotr12 + 3*dy*n2dotr12*sp1 - cp1p2*R122*st2))/R125

!              [2] Off-diagonal terms on the diagonal blocks: same atom, different coordinates
!              x1,y1
               DUMMY =
     &            (3*(MU2*R129*(5*cp1*dy*n2dotr12*st1 + 5*cp2*dy*n1dotr12*st2 -
     &            R122*sp1pp2*st1*st2) +
     &            dx*(dy*(224 - 64*stocklambda*R126 + 5*MU2*(-7*n1dotr12*n2dotr12*R127 +
     &            n1dotn2*R129)) + 5*MU2*R129*(n2dotr12*sp1*st1 + n1dotr12*sp2*st2))))/
     &            R1216
               HESS(J3-2, J3-1) = HESS(J3-2, J3-1) + DUMMY
               HESS(J3-1, J3-2) = HESS(J3-1, J3-2) + DUMMY
!              x1,z1
               DUMMY =
     &            (-3*(dx*dz*(-224 + 64*stocklambda*R126 + 35*MU2*n1dotr12*n2dotr12*R127) -
     &            5*dx*MU2*(dz*n1dotn2 + ct2*n1dotr12 + ct1*n2dotr12)*R129 +
     &            MU2*R129*(cp1*(-5*dz*n2dotr12 + ct2*R122)*st1 +
     &            cp2*(-5*dz*n1dotr12 + ct1*R122)*st2)))/R1216
               HESS(J3-2, J3) = HESS(J3-2, J3) + DUMMY
               HESS(J3, J3-2) = HESS(J3, J3-2) + DUMMY
!              y1,z1
               DUMMY =
     &            (3*(dy*dz*(224 - 64*stocklambda*R126 - 35*MU2*n1dotr12*n2dotr12*R127) +
     &            5*dy*MU2*(dz*n1dotn2 + ct2*n1dotr12 + ct1*n2dotr12)*R129 +
     &            MU2*R129*(-(R122*(ct2*sp1*st1 + ct1*sp2*st2)) +
     &            5*dz*(n2dotr12*sp1*st1 + n1dotr12*sp2*st2))))/R1216
               HESS(J3-1, J3) = HESS(J3-1, J3) + DUMMY
               HESS(J3, J3-1) = HESS(J3, J3-1) + DUMMY
!              x1,t1
               DUMMY =
     &            (-3*MU2*(st1*(5*dx*dz*n2dotr12 - ct2*dx*R122 - cp2*dz*R122*st2) +
     &            ct1*sp1*(-5*dx*dy*n2dotr12 + cp2*dy*R122*st2 + dx*R122*sp2*st2) +
     &            cp1*ct1*(-5*dx2*n2dotr12 + R122*(n2dotr12 + 2*cp2*dx*st2))))/R127
               HESS(J3-2, OFFSET+J3-2) = HESS(J3-2, OFFSET+J3-2) + DUMMY
               HESS(OFFSET+J3-2, J3-2) = HESS(OFFSET+J3-2, J3-2) + DUMMY
!              x1,p1
               DUMMY =
     &            (3*MU2*st1*(5*cp1*dx*dy*n2dotr12 - cp1*R122*(cp2*dy + dx*sp2)*st2 +
     &            sp1*(-5*dx2*n2dotr12 + R122*(n2dotr12 + 2*cp2*dx*st2))))/R127
               HESS(J3-2, OFFSET+J3-1) = HESS(J3-2, OFFSET+J3-1) + DUMMY
               HESS(OFFSET+J3-1, J3-2) = HESS(OFFSET+J3-1, J3-2) + DUMMY
!              y1,t1
               DUMMY =
     &            (-3*MU2*(cp1*ct1*(-5*dx*dy*n2dotr12 + cp2*dy*R122*st2 + dx*R122*sp2*st2) +
     &            st1*(5*dy*dz*n2dotr12 - ct2*dy*R122 - dz*R122*sp2*st2) +
     &            ct1*sp1*(-5*dy2*n2dotr12 + R122*(n2dotr12 + 2*dy*sp2*st2))))/R127
               HESS(J3-1, OFFSET+J3-2) = HESS(J3-1, OFFSET+J3-2) + DUMMY
               HESS(OFFSET+J3-2, J3-1) = HESS(OFFSET+J3-2, J3-1) + DUMMY
!              y1,p1
               DUMMY =
     &            (3*MU2*st1*(n2dotr12*(5*cp1*dy2 - cp1*R122 - 5*dx*dy*sp1) +
     &            R122*(cp2*dy*sp1 - 2*cp1*dy*sp2 + dx*sp1*sp2)*st2))/R127
               HESS(J3-1, OFFSET+J3-1) = HESS(J3-1, OFFSET+J3-1) + DUMMY
               HESS(OFFSET+J3-1, J3-1) = HESS(OFFSET+J3-1, J3-1) + DUMMY
!              z1,t1
               DUMMY =
     &            (-3*MU2*((5*dz2*n2dotr12 - (2*ct2*dz + n2dotr12)*R122)*st1 +
     &            cp1*ct1*(-5*dx*dz*n2dotr12 + ct2*dx*R122 + cp2*dz*R122*st2) +
     &            ct1*sp1*(-5*dy*dz*n2dotr12 + ct2*dy*R122 + dz*R122*sp2*st2)))/R127
               HESS(J3, OFFSET+J3-2) = HESS(J3, OFFSET+J3-2) + DUMMY
               HESS(OFFSET+J3-2, J3) = HESS(OFFSET+J3-2, J3) + DUMMY
!              z1,p1
               DUMMY =
     &            (3*MU2*st1*((5*dz*n2dotr12 - ct2*R122)*(cp1*dy - dx*sp1) +
     &            dz*R122*sp1mp2*st2))/R127
               HESS(J3, OFFSET+J3-1) = HESS(J3, OFFSET+J3-1) + DUMMY
               HESS(OFFSET+J3-1, J3) = HESS(OFFSET+J3-1, J3) + DUMMY
!              t1,p1
               DUMMY =
     &            -((ct1*MU2*(3*cp1*dy*n2dotr12 - 3*dx*n2dotr12*sp1 + R122*sp1mp2*st2))/R125)
               HESS(OFFSET+J3-2, OFFSET+J3-1) = HESS(OFFSET+J3-2, OFFSET+J3-1) + DUMMY
               HESS(OFFSET+J3-1, OFFSET+J3-2) = HESS(OFFSET+J3-1, OFFSET+J3-2) + DUMMY

!              [3] Diagonal elements on off-diagonal blocks: different particle, same coordinate
!              x1,x2
               HESS(J3-2, J4-2) =
     &            (3*(MU2*n1dotn2*R1211 + 16*R122 - 8*stocklambda*R128 +
     &            dx2*(-224 + 64*stocklambda*R126 + 35*MU2*n1dotr12*n2dotr12*R127 -
     &            5*MU2*n1dotn2*R129) - MU2*R129*(5*n1dotr12*(n2dotr12 + 2*cp2*dx*st2) +
     &            2*cp1*st1*(5*dx*n2dotr12 - cp2*R122*st2))))/R1216
!              y1,y2
               HESS(J3-1, J4-1) =
     &            (3*(-224*dy2 + MU2*n1dotn2*R1211 + 16*R122 + 64*dy2*stocklambda*R126 +
     &            35*dy2*MU2*n1dotr12*n2dotr12*R127 - 8*stocklambda*R128 -
     &            5*MU2*(dy2*n1dotn2 + n1dotr12*n2dotr12)*R129 +
     &            2*MU2*R129*(R122*sp1*sp2*st1*st2 - 5*dy*(n2dotr12*sp1*st1 +
     &            n1dotr12*sp2*st2))))/R1216
!              z1,z2
               HESS(J3, J4) =
     &            (3*(-224*dz2 + MU2*n1dotn2*R1211 + 16*R122 + 64*dz2*stocklambda*R126 +
     &            35*dz2*MU2*n1dotr12*n2dotr12*R127 - 8*stocklambda*R128 -
     &            MU2*(5*(dz2*n1dotn2 + 2*ct2*dz*n1dotr12 + (2*ct1*dz + n1dotr12)*
     &            n2dotr12) - 2*ct1*ct2*R122)*R129))/R1216
!              t1,t2
               HESS(OFFSET+J3-2, OFFSET+J4-2) =
     &            (MU2*(-3*(cp1*ct1*dx + ct1*dy*sp1 - dz*st1)*(cp2*ct2*dx + ct2*dy*sp2 -
     &            dz*st2) + R122*(cp1p2*ct1*ct2 + st1*st2)))/R125
!              p1,p2
               HESS(OFFSET+J3-1, OFFSET+J4-1) =
     &            (MU2*(cp1p2*R122 - 3*(cp1*dy - dx*sp1)*(cp2*dy - dx*sp2))*st1*st2)/R125

!              [4] Completely off-diagonal terms: different particle, different coordinate
!              x1,y2 and y1,x2
               HESS(J3-2, J4-1) =
     &            (3*(dx*dy*(-224 + 64*stocklambda*R126 + 35*MU2*n1dotr12*n2dotr12*R127 -
     &            5*MU2*n1dotn2*R129) - MU2*R129*(5*n2dotr12*(cp1*dy + dx*sp1)*st1 +
     &            (5*cp2*dy*n1dotr12 + 5*dx*n1dotr12*sp2 - R122*sp1pp2*st1)*st2)))/R1216
               HESS(J3-1, J4-2) = HESS(J3-2, J4-1)
!              x1,z2 and z1,x2
               HESS(J3-2, J4) =
     &            (3*(dx*dz*(-224 + 64*stocklambda*R126 + 35*MU2*n1dotr12*n2dotr12*R127) -
     &            5*dx*MU2*(dz*n1dotn2 + ct2*n1dotr12 + ct1*n2dotr12)*R129 +
     &            MU2*R129*(cp1*(-5*dz*n2dotr12 + ct2*R122)*st1 +
     &            cp2*(-5*dz*n1dotr12 + ct1*R122)*st2)))/R1216
               HESS(J3, J4-2) = HESS(J3-2, J4)
!              y1,z2 and z1,y2
               HESS(J3-1, J4) =
     &            (3*(dy*dz*(-224 + 64*stocklambda*R126 + 35*MU2*n1dotr12*n2dotr12*R127) -
     &            5*dy*MU2*(dz*n1dotn2 + ct2*n1dotr12 + ct1*n2dotr12)*R129 +
     &            MU2*R129*(R122*(ct2*sp1*st1 + ct1*sp2*st2) -
     &            5*dz*(n2dotr12*sp1*st1 + n1dotr12*sp2*st2))))/R1216
               HESS(J3, J4-1) = HESS(J3-1, J4)
!              x1,t2
               HESS(J3-2, OFFSET+J4-2) =
     &            (-3*MU2*(ct2*sp2*(-5*dx*dy*n1dotr12 + cp1*dy*R122*st1 + dx*R122*sp1*st1) +
     &            cp2*ct2*(-5*dx2*n1dotr12 + R122*(n1dotr12 + 2*cp1*dx*st1)) +
     &            (5*dx*dz*n1dotr12 - ct1*dx*R122 - cp1*dz*R122*st1)*st2))/R127
!              x1,p2
               HESS(J3-2, OFFSET+J4-1) =
     &            (-3*MU2*(-(n1dotr12*(5*cp2*dx*dy + (-5*dx2 + R122)*sp2)) +
     &            R122*(cp1*cp2*dy + cp2*dx*sp1 - 2*cp1*dx*sp2)*st1)*st2)/R127
!              y1,t2
               HESS(J3-1, OFFSET+J4-2) =
     &            (-3*MU2*(cp2*ct2*(-5*dx*dy*n1dotr12 + cp1*dy*R122*st1 + dx*R122*sp1*st1) +
     &            ct2*sp2*(-5*dy2*n1dotr12 + R122*(n1dotr12 + 2*dy*sp1*st1)) +
     &            (5*dy*dz*n1dotr12 - ct1*dy*R122 - dz*R122*sp1*st1)*st2))/R127
!              y1,p2
               HESS(J3-1, OFFSET+J4-1) =
     &            (3*MU2*(n1dotr12*(5*cp2*dy2 - cp2*R122 - 5*dx*dy*sp2) +
     &            R122*(-2*cp2*dy*sp1 + cp1*dy*sp2 + dx*sp1*sp2)*st1)*st2)/R127
!              z1,t2
               HESS(J3, OFFSET+J4-2) =
     &            (-3*MU2*(cp2*ct2*(-5*dx*dz*n1dotr12 + ct1*dx*R122 + cp1*dz*R122*st1) +
     &            ct2*sp2*(-5*dy*dz*n1dotr12 + ct1*dy*R122 + dz*R122*sp1*st1) +
     &            (5*dz2*n1dotr12 - (2*ct1*dz + n1dotr12)*R122)*st2))/R127
!              z1,p2
               HESS(J3, OFFSET+J4-1) =
     &            (3*MU2*((5*dz*n1dotr12 - ct1*R122)*(cp2*dy - dx*sp2) - dz*R122*sp1mp2*st1)*
     &            st2)/R127
!              t1,x2
               HESS(OFFSET+J3-2, J4-2) =
     &            (3*MU2*(-(st1*(-5*dx*dz*n2dotr12 + ct2*dx*R122 + cp2*dz*R122*st2)) +
     &            ct1*sp1*(-5*dx*dy*n2dotr12 + cp2*dy*R122*st2 + dx*R122*sp2*st2) +
     &            cp1*ct1*(-5*dx2*n2dotr12 + R122*(n2dotr12 + 2*cp2*dx*st2))))/R127
!              p1,x2
               HESS(OFFSET+J3-1, J4-2) =
     &            (3*MU2*st1*(-(n2dotr12*(5*cp1*dx*dy + (-5*dx2 + R122)*sp1)) +
     &            R122*(cp1*cp2*dy - 2*cp2*dx*sp1 + cp1*dx*sp2)*st2))/R127
!              t1,y2
               HESS(OFFSET+J3-2, J4-1) =
     &            (3*MU2*(cp1*ct1*(-5*dx*dy*n2dotr12 + cp2*dy*R122*st2 + dx*R122*sp2*st2) +
     &            st1*(5*dy*dz*n2dotr12 - ct2*dy*R122 - dz*R122*sp2*st2) +
     &            ct1*sp1*(-5*dy2*n2dotr12 + R122*(n2dotr12 + 2*dy*sp2*st2))))/R127
!              p1,y2
               HESS(OFFSET+J3-1, J4-1) =
     &            (3*MU2*st1*(5*dx*dy*n2dotr12*sp1 - R122*sp1*(cp2*dy + dx*sp2)*st2 +
     &            cp1*(-5*dy2*n2dotr12 + R122*(n2dotr12 + 2*dy*sp2*st2))))/R127
!              t1,z2
               HESS(OFFSET+J3-2, J4) =
     &            (3*MU2*((5*dz2*n2dotr12 - (2*ct2*dz + n2dotr12)*R122)*st1 +
     &            cp1*ct1*(-5*dx*dz*n2dotr12 + ct2*dx*R122 + cp2*dz*R122*st2) +
     &            ct1*sp1*(-5*dy*dz*n2dotr12 + ct2*dy*R122 + dz*R122*sp2*st2)))/R127
!              p1,z2
               HESS(OFFSET+J3-1, J4) =
     &            (3*MU2*st1*(-((5*dz*n2dotr12 - ct2*R122)*(cp1*dy - dx*sp1)) -
     &            dz*R122*sp1mp2*st2))/R127
!              t1,p2
               HESS(OFFSET+J3-2, OFFSET+J4-1) =
     &            (MU2*(ct1*R122*sp1mp2 + 3*(-(cp2*dy) + dx*sp2)*(cp1*ct1*dx + ct1*dy*sp1 -
     &            dz*st1))*st2)/R125
!              p1,t2
               HESS(OFFSET+J3-1, OFFSET+J4-2) =
     &            (MU2*st1*(-(ct2*R122*sp1mp2) + 3*(-(cp1*dy) + dx*sp1)*
     &            (cp2*ct2*dx + ct2*dy*sp2 - dz*st2)))/R125

            END DO

         END DO
      ENDIF

      RETURN
      END
C
C  Orthogonalise VEC1 to overall translations and rotations.
C  Still not working properly for no apparent reason!
C
      SUBROUTINE ORTHOGSTOCK(VEC1,COORDS,OTEST) 
      USE COMMONS
      USE KEY
      USE VECCK
      IMPLICIT NONE
      INTEGER J2, J3, NCHECK
      DOUBLE PRECISION COORDS(*), VEC1(*), DUMMY1, DUMMY2, ROOTN, VDOT, DUMMYX, DUMMYY, DUMMYZ, DDOT,
     1                 CMX, CMY, CMZ, AMASS(NATOMS), TMASS, RMASS(NATOMS), THETA, PHI, DUMMY,
     &                 ETX(3*NATOMS), ETY(3*NATOMS), ETZ(3*NATOMS), ERX(3*NATOMS), ERY(3*NATOMS), ERZ(3*NATOMS)
      LOGICAL OTEST

      IF (.NOT.ALLOCATED(VECCHK)) ALLOCATE(VECCHK(3*NATOMS,30))
      ROOTN=SQRT(1.0D0*(NATOMS/2)) ! factor of 2 for Stockmayer
      NCHECK=0
      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      TMASS=0.0D0
      DO J2=1,(NATOMS/2) ! factor of 2 for Stockmayer
         AMASS(J2)=1.0D0
         RMASS(J2)=1.0D0
         IF (MASST) AMASS(J2)=ATMASS(J2)
         IF (MASST) RMASS(J2)=SQRT(ATMASS(J2))
         TMASS=TMASS+AMASS(J2)
      ENDDO
C
C  If MASST then the coordinates already have a square root of the mass in them
C
      DO J2=1,(NATOMS/2) ! factor of 2 for Stockmayer
         CMX=CMX+COORDS(3*(J2-1)+1)*RMASS(J2)
         CMY=CMY+COORDS(3*(J2-1)+2)*RMASS(J2)
         CMZ=CMZ+COORDS(3*(J2-1)+3)*RMASS(J2)
      ENDDO
      CMX=CMX/TMASS
      CMY=CMY/TMASS
      CMZ=CMZ/TMASS
C
C  Orthogonalize to known eigenvectors corresponding to negative eigenvalues
C  for CHECKINDEX run with NOHESS.
C
      IF (CHECKINDEX) THEN
         DO J2=1,NMDONE
            DUMMY1=0.0D0
            DO J3=1,NOPT
               DUMMY1=DUMMY1+VECCHK(J3,J2)*VEC1(J3)
            ENDDO
            VDOT=MAX(VDOT,ABS(DUMMY1))
            DUMMY2=0.0D0
            DO J3=1,NOPT
               VEC1(J3)=VEC1(J3)-DUMMY1*VECCHK(J3,J2)
               DUMMY2=DUMMY2+VEC1(J3)**2
            ENDDO
            DUMMY2=1.0D0/DSQRT(DUMMY2)
            DO J3=1,NOPT
               VEC1(J3)=VEC1(J3)*DUMMY2
            ENDDO
         ENDDO
         IF (OTEST) CALL VECNORM(VEC1,NOPT)
      ENDIF
!
!  Gram-Scmidt orthogonalisation of overall trans/rot eigenvectors.
!
      ETX(1:3*NATOMS)=0.0D0; ETY(1:3*NATOMS)=0.0D0; ETZ(1:3*NATOMS)=0.0D0
      ERX(1:3*NATOMS)=0.0D0; ERY(1:3*NATOMS)=0.0D0; ERZ(1:3*NATOMS)=0.0D0

      DO J2=1,(NATOMS/2) 
         ETX(3*(J2-1)+1)=1.0D0/ROOTN
         ETY(3*(J2-1)+2)=1.0D0/ROOTN
         ETZ(3*(J2-1)+3)=1.0D0/ROOTN
      ENDDO

      DUMMYX=0.0D0
      DUMMYY=0.0D0
      DUMMYZ=0.0D0
      DO J2=1,(NATOMS/2) ! rotation 
         J3=3*J2
         ERX(3*(J2-1)+2)=  COORDS(J3)  -CMZ
         ERX(3*(J2-1)+3)=-(COORDS(J3-1)-CMY)
         DUMMYX=DUMMYX+ERX(3*(J2-1)+2)**2+ERX(3*(J2-1)+3)**2
         ERY(3*(J2-1)+1)=  COORDS(J3)  -CMZ
         ERY(3*(J2-1)+3)=-(COORDS(J3-2)-CMX)
         DUMMYY=DUMMYY+ERY(3*(J2-1)+1)**2+ERY(3*(J2-1)+3)**2
         ERZ(3*(J2-1)+1)=  COORDS(J3-1)-CMY
         ERZ(3*(J2-1)+2)=-(COORDS(J3-2)-CMX)
         DUMMYZ=DUMMYZ+ERZ(3*(J2-1)+1)**2+ERZ(3*(J2-1)+2)**2
      ENDDO
      DO J2=(NATOMS/2)+1,NATOMS ! xrotation for theta, phi coordinates
         THETA=COORDS(3*(J2-1)+1)
         PHI=COORDS(3*(J2-1)+2)
         DUMMY=TAN(THETA)
         ERX(3*(J2-1)+1)=SIN(PHI)
         IF (DUMMY.NE.0.0D0) ERX(3*(J2-1)+2)=COS(PHI)/DUMMY
         DUMMYX=DUMMYX+ERX(3*(J2-1)+1)**2+ERX(3*(J2-1)+2)**2
         ERY(3*(J2-1)+1)=COS(PHI)
         IF (DUMMY.NE.0.0D0) ERY(3*(J2-1)+2)=-SIN(PHI)/DUMMY
         DUMMYY=DUMMYY+ERY(3*(J2-1)+1)**2+ERY(3*(J2-1)+2)**2
         ERZ(3*(J2-1)+2)=-1.0D0
         DUMMYZ=DUMMYZ+ERZ(3*(J2-1)+2)**2
       ENDDO 
       DUMMYX=SQRT(DUMMYX); DUMMYY=SQRT(DUMMYY); DUMMYZ=SQRT(DUMMYZ)
       DO J2=1,3*NATOMS
          ERX(J2)=ERX(J2)/DUMMYX
          ERY(J2)=ERY(J2)/DUMMYY
          ERZ(J2)=ERZ(J2)/DUMMYZ
       ENDDO
!
! Orthogonalise ERY to ERX
!
       DUMMY=0.0D0
       DO J2=1,3*NATOMS
          DUMMY=DUMMY+ERX(J2)*ERY(J2)
       ENDDO
       DUMMY2=0.0D0
       DO J2=1,3*NATOMS
          ERY(J2)=ERY(J2)-DUMMY*ERX(J2)
          DUMMY2=DUMMY2+ERY(J2)**2
       ENDDO
       DUMMY2=SQRT(DUMMY2)
       DO J2=1,3*NATOMS
          ERY(J2)=ERY(J2)/DUMMY2
       ENDDO
!
! Orthogonalise ERZ to ERX
!
       DUMMY=0.0D0
       DO J2=1,3*NATOMS
          DUMMY=DUMMY+ERX(J2)*ERZ(J2)
       ENDDO
       DUMMY2=0.0D0
       DO J2=1,3*NATOMS
          ERZ(J2)=ERZ(J2)-DUMMY*ERX(J2)
          DUMMY2=DUMMY2+ERZ(J2)**2
       ENDDO
       DUMMY2=SQRT(DUMMY2)
       DO J2=1,3*NATOMS
          ERZ(J2)=ERZ(J2)/DUMMY2
       ENDDO
!
! Orthogonalise ERZ to ERY
!
       DUMMY=0.0D0
       DO J2=1,3*NATOMS
          DUMMY=DUMMY+ERY(J2)*ERZ(J2)
       ENDDO
       DUMMY2=0.0D0
       DO J2=1,3*NATOMS
          ERZ(J2)=ERZ(J2)-DUMMY*ERY(J2)
          DUMMY2=DUMMY2+ERZ(J2)**2
       ENDDO
       DUMMY2=SQRT(DUMMY2)
       DO J2=1,3*NATOMS
          ERZ(J2)=ERZ(J2)/DUMMY2
       ENDDO
!
!  Check orthogonalisation of trans/rot
!
!       DUMMY=DDOT(3*NATOMS,ETX,1,ETX,1)
!       PRINT *,'etx,etx=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETY,1,ETY,1)
!       PRINT *,'ety,ety=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETZ,1,ETZ,1)
!       PRINT *,'etz,etz=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ERX,1,ERX,1)
!       PRINT *,'erx,erx=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ERY,1,ERY,1)
!       PRINT *,'ery,ery=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ERZ,1,ERZ,1)
!       PRINT *,'erz,erz=',DUMMY
!
!       DUMMY=DDOT(3*NATOMS,ETX,1,ETY,1)
!       PRINT *,'etx,ety=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETX,1,ETZ,1)
!       PRINT *,'etx,etz=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETX,1,ERX,1)
!       PRINT *,'etx,erx=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETX,1,ERY,1)
!       PRINT *,'etx,ery=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETX,1,ERZ,1)
!       PRINT *,'etx,erz=',DUMMY
!
!       DUMMY=DDOT(3*NATOMS,ETY,1,ETZ,1)
!       PRINT *,'ety,etz=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETY,1,ERX,1)
!       PRINT *,'ety,erx=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETY,1,ERY,1)
!       PRINT *,'ety,ery=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETY,1,ERZ,1)
!       PRINT *,'ety,erz=',DUMMY
!    
!       DUMMY=DDOT(3*NATOMS,ETZ,1,ERX,1)
!       PRINT *,'etz,erx=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETZ,1,ERY,1)
!       PRINT *,'etz,ery=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ETZ,1,ERZ,1)
!       PRINT *,'etz,erz=',DUMMY
!
!       DUMMY=DDOT(3*NATOMS,ERX,1,ERY,1)
!       PRINT *,'erx,ery=',DUMMY
!       DUMMY=DDOT(3*NATOMS,ERX,1,ERZ,1)
!       PRINT *,'erx,erz=',DUMMY
!
!       DUMMY=DDOT(3*NATOMS,ERY,1,ERZ,1)
!       PRINT *,'ery,erz=',DUMMY
!
!  Project trans/rot set out of VEC1
!
       DUMMY=DDOT(3*NATOMS,ETX,1,VEC1,1)
       DO J2=1,3*NATOMS
          VEC1(J2)=VEC1(J2)-DUMMY*ETX(J2)
       ENDDO
       DUMMY=DDOT(3*NATOMS,ETY,1,VEC1,1)
       DO J2=1,3*NATOMS
          VEC1(J2)=VEC1(J2)-DUMMY*ETY(J2)
       ENDDO
       DUMMY=DDOT(3*NATOMS,ETZ,1,VEC1,1)
       DO J2=1,3*NATOMS
          VEC1(J2)=VEC1(J2)-DUMMY*ETZ(J2)
       ENDDO
       DUMMY=DDOT(3*NATOMS,ERX,1,VEC1,1)
       DO J2=1,3*NATOMS
          VEC1(J2)=VEC1(J2)-DUMMY*ERX(J2)
       ENDDO
       DUMMY=DDOT(3*NATOMS,ERY,1,VEC1,1)
       DO J2=1,3*NATOMS
          VEC1(J2)=VEC1(J2)-DUMMY*ERY(J2)
       ENDDO
       DUMMY=DDOT(3*NATOMS,ERZ,1,VEC1,1)
       DO J2=1,3*NATOMS
          VEC1(J2)=VEC1(J2)-DUMMY*ERZ(J2)
       ENDDO
!
! The sixth coordinate for each particle is a dummy.
!
      DO J2=(NATOMS/2)+1,NATOMS
         VEC1(3*(J2-1)+3)=0.0D0
      ENDDO
      IF (OTEST) CALL VECNORM(VEC1,NOPT)

      RETURN
      END

C
      SUBROUTINE SHIFTSTOCK(COORDS,NATOMS)
      USE KEY
      USE MODHESS
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4, NOPT, NATOMS, NIADD
      DOUBLE PRECISION CMX, CMY, CMZ, COORDS(3*NATOMS), DUMMYX, DUMMYY, DUMMYZ, DDOT, THETA, PHI, DUMMY, VEC6(3*NATOMS),
     &                 ETX(3*NATOMS), ETY(3*NATOMS), ETZ(3*NATOMS), ERX(3*NATOMS), ERY(3*NATOMS), ERZ(3*NATOMS), ROOTN

      CMX=0.0D0
      CMY=0.0D0
      CMZ=0.0D0
      DO J2=1,(NATOMS/2)
         CMX=CMX+COORDS(3*(J2-1)+1)
         CMY=CMY+COORDS(3*(J2-1)+2)
         CMZ=CMZ+COORDS(3*(J2-1)+3)
      ENDDO
      CMX=CMX/(NATOMS/2)
      CMY=CMY/(NATOMS/2)
      CMZ=CMZ/(NATOMS/2)
      ETX(1:3*NATOMS)=0.0D0; ETY(1:3*NATOMS)=0.0D0; ETZ(1:3*NATOMS)=0.0D0
      ERX(1:3*NATOMS)=0.0D0; ERY(1:3*NATOMS)=0.0D0; ERZ(1:3*NATOMS)=0.0D0

      ROOTN=SQRT(1.0D0*(NATOMS/2)) ! factor of 2 for Stockmayer
      DO J2=1,(NATOMS/2) 
         ETX(3*(J2-1)+1)=1.0D0/ROOTN
         ETY(3*(J2-1)+2)=1.0D0/ROOTN
         ETZ(3*(J2-1)+3)=1.0D0/ROOTN
      ENDDO

      DUMMYX=0.0D0
      DUMMYY=0.0D0
      DUMMYZ=0.0D0
      DO J2=1,(NATOMS/2) ! rotation 
         J3=3*J2
         ERX(3*(J2-1)+2)=  COORDS(J3)  -CMZ
         ERX(3*(J2-1)+3)=-(COORDS(J3-1)-CMY)
         DUMMYX=DUMMYX+ERX(3*(J2-1)+2)**2+ERX(3*(J2-1)+3)**2
         ERY(3*(J2-1)+1)=  COORDS(J3)  -CMZ
         ERY(3*(J2-1)+3)=-(COORDS(J3-2)-CMX)
         DUMMYY=DUMMYY+ERY(3*(J2-1)+1)**2+ERY(3*(J2-1)+3)**2
         ERZ(3*(J2-1)+1)=  COORDS(J3-1)-CMY
         ERZ(3*(J2-1)+2)=-(COORDS(J3-2)-CMX)
         DUMMYZ=DUMMYZ+ERZ(3*(J2-1)+1)**2+ERZ(3*(J2-1)+2)**2
      ENDDO
      DO J2=(NATOMS/2)+1,NATOMS ! xrotation for theta, phi coordinates
         THETA=COORDS(3*(J2-1)+1)
         PHI=COORDS(3*(J2-1)+2)
         DUMMY=TAN(THETA)
         ERX(3*(J2-1)+1)=SIN(PHI)
         IF (DUMMY.NE.0.0D0) ERX(3*(J2-1)+2)=COS(PHI)/DUMMY
         DUMMYX=DUMMYX+ERX(3*(J2-1)+1)**2+ERX(3*(J2-1)+2)**2
         ERY(3*(J2-1)+1)=COS(PHI)
         IF (DUMMY.NE.0.0D0) ERY(3*(J2-1)+2)=-SIN(PHI)/DUMMY
         DUMMYY=DUMMYY+ERY(3*(J2-1)+1)**2+ERY(3*(J2-1)+2)**2
         ERZ(3*(J2-1)+2)=-1.0D0
         DUMMYZ=DUMMYZ+ERZ(3*(J2-1)+2)**2
      ENDDO 
      DUMMYX=SQRT(DUMMYX); DUMMYY=SQRT(DUMMYY); DUMMYZ=SQRT(DUMMYZ)
      DO J2=1,3*NATOMS
         ERX(J2)=ERX(J2)/DUMMYX
         ERY(J2)=ERY(J2)/DUMMYY
         ERZ(J2)=ERZ(J2)/DUMMYZ
      ENDDO
!
! The sixth coordinate for each particle is a dummy.
!
      DO J2=(NATOMS/2)+1,NATOMS
         HESS(3*(J2-1)+3,3*(J2-1)+3)=HESS(3*(J2-1)+3,3*(J2-1)+3)+SHIFTL(1)
      ENDDO

C     PRINT*,'unshifted hess:'
C     WRITE(*,'(6F15.5)') ((HESS(J1,J2),J1=1,NOPT),J2=1,NOPT)

      DO J1=1,3*NATOMS
         DO J2=1,3*NATOMS
            HESS(J2,J1)=HESS(J2,J1)
     &                             +SHIFTL(1)*ETX(J2)*ETX(J1)
     &                             +SHIFTL(2)*ETY(J2)*ETY(J1)
     &                             +SHIFTL(3)*ETZ(J2)*ETZ(J1)
     &                             +SHIFTL(4)*ERX(J2)*ERX(J1)
     &                             +SHIFTL(5)*ERY(J2)*ERY(J1)
     &                             +SHIFTL(6)*ERZ(J2)*ERZ(J1)
         ENDDO
      ENDDO

C     PRINT*,'shifted hess:'
C     WRITE(*,'(6F15.5)') ((HESS(J1,J2),J1=1,NOPT),J2=1,NOPT)
C     PRINT*,'coordinates:'
C     WRITE(*,'(6F15.5)') (Q(J1),J1=1,NOPT)

      RETURN
      END


C     Test for z-aligned dipoles, returning TRUE if found and FALSE otherwise.
      LOGICAL FUNCTION ZALIGNTEST(NATOMS, X)
         USE KEY,ONLY : STOCKZTOL, STOCKMAXSPIN

         IMPLICIT NONE
         INTEGER, INTENT(IN) :: NATOMS
         DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS)

         INTEGER J1

         DO J1 = NATOMS/2+1, NATOMS
            IF ( (1.0D0 - ABS(COS(X(J1*3-2)))) .LT. STOCKZTOL) THEN
               ZALIGNTEST = .TRUE.
               RETURN
            END IF
         END DO

         ZALIGNTEST = .FALSE.

      END FUNCTION ZALIGNTEST
