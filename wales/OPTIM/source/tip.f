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

C
C  Energy and gradient for rigid body TIPnP using new rigid body
C  derivative functions etc.
C
      SUBROUTINE TIPnP(TIPID,X,V,ETIP,GTEST,SECT)

      USE commons
      use key, only : efield
      USE MODHESS

      IMPLICIT NONE

! Subroutine arguments

      integer, intent(IN) :: TIPID                   ! TIPnP, n=1-5
      DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS)    ! COORDINATES
      DOUBLE PRECISION, INTENT(OUT) :: V(3*NATOMS)   ! GRADIENTS
      DOUBLE PRECISION, INTENT(OUT) :: ETIP          ! ENERGY
      LOGICAL, intent(IN) :: GTEST                   ! Calculate gradients
      LOGICAL, intent(IN) :: SECT                    ! Calculate second derivatives

! Local variables

      INTEGER J1, J2, K1, K2, NAT2, NSITE, J3, J4, J3MIN, J4MIN
      DOUBLE PRECISION X1, X2, Y1, Y2, Z1, Z2, FLJ, DFLJ, D2FLJ, FC, DFC, D2FC,
     1                 pairEnergy, EPS2, RALPHA12, RALPHA22, RDIST, 
     2                 M1, L1, N1, M2, L2, N2, ALPHA1, CA1, CA2, ALPHA2, S1, S2, C3A1, C3A2,
     3                 L12, M12, N12, L22, M22, N22, SIGMA, C6, C12,
     4                 GX1, GY1, GZ1, GL1, GM1, GN1, GL2, GM2, GN2, C2A2, C2A1,
     5                 DIST,
     6                 P1X, P1Y, P1Z, P2X, P2Y, P2Z, XDUMM, Q1, Q2,
     7                 CHARGE(5), HTOKJ, SITE(5,3)
      PARAMETER (HTOKJ=1389.354848D0) ! conversion factor for coulomb energy and gradients

      integer :: RBcoord1, RBcoord2, index1, index2
      double precision :: firstDerivatives(12)
      double precision :: pairGrad(12), pairHess(12,12)
      double precision :: dEdRij, d2EdRij2
      double precision :: singleGrad(6), singleHess(6,6)  ! For use with external fields

      double precision :: zReal, efieldInteraction, pDOTxzero, RALPHA1, eVtoKJ
      parameter (eVtoKJ=96.4847D0) ! eV to kJ/mol conversion, input field in V/A

! External functions

      double precision :: dRijdXrb, d2RijdXrb1dXrb2
      double precision :: dzdXrb, d2zdXrb1dXrb2
      
C
C  Statement functions.
C  Site-site energy terms - Coulombic and LJ. XDUMM is the site-site distance.
C
      FLJ(XDUMM)=(C12/XDUMM**6-C6)/XDUMM**6
      DFLJ(XDUMM)=-6.0D0*(-C6+2.0D0*C12/XDUMM**6)/XDUMM**7
      D2FLJ(XDUMM)=6.0D0*(26.0d0*C12/XDUMM**6 - 7.0d0*C6)/XDUMM**8

      FC(Q1,Q2,XDUMM)=Q1*Q2/XDUMM
      DFC(Q1,Q2,XDUMM)=-Q1*Q2/XDUMM**2
      D2FC(Q1,Q2,XDUMM)=2.0d0*Q1*Q2/XDUMM**3

C     PRINT*,'coords in TIP:'
C     WRITE(*,'(3F20.10)') (X(J1),J1=1,3*NATOMS)
C
C Distinguish TIP1P, TIP2P, TIP3P, TIP4P
C
      SITE(1,1)=0.0D0
      SITE(1,2)=0.0D0
      SITE(1,3)=0.0D0
C      SITE(1,3)=0.065098D0
      SITE(2,1)=0.756950327D0
      SITE(2,2)=0.0D0
      SITE(2,3)=-0.5858822760D0
C      SITE(2,3)=-0.5207842760D0
      SITE(3,1)=-0.756950327D0
      SITE(3,2)=0.0D0
      SITE(3,3)=-0.5858822760D0
C      SITE(3,3)=-0.5207842760D0
      SITE(4,1)=0.0D0
      SITE(4,2)=0.0D0
      CHARGE(1)=0.0D0
      NSITE=4

      IF (TIPID.EQ.1) THEN
         C6=2510.4D0 ! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
         C12=2426720.0D0
         SITE(4,3)=0.0D0
C        SITE(4,3)=0.065098D0
         CHARGE(2)=0.4D0
         CHARGE(3)=0.4D0
         CHARGE(4)=-0.8D0
      ELSE IF (TIPID.EQ.2) THEN
         C6=2510.4D0 ! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
         C12=2907880.0D0
         SITE(4,3)=-0.15D0
C        SITE(4,3)=-0.084902D0
         CHARGE(2)=0.535D0
         CHARGE(3)=0.535D0
         CHARGE(4)=-1.07D0
      ELSE IF (TIPID.EQ.3) THEN
         C6=2489.48D0 ! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
         C12=2435088.0D0
         SITE(4,3)=0.0D0
C        SITE(4,3)=0.065098D0
         CHARGE(2)=0.417D0
         CHARGE(3)=0.417D0
         CHARGE(4)=-0.834D0
      ELSE IF (TIPID.EQ.4) THEN
         C6=2552.24D0 ! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
         C12=2510.4D3
         SITE(4,3)=-0.15D0
C        SITE(4,3)=-0.084902D0
         CHARGE(2)=0.52D0
         CHARGE(3)=0.52D0
         CHARGE(4)=-1.04D0
      ELSE IF (TIPID.EQ.5) THEN
! C6  = 590.3472412D0 kcal/mol Angstrom**6
! C12 = 544546.6644D0 kcal/mol Angstrom**12 
         C6=2470.012857D0 ! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
         C12=2278383.244D0
         NSITE=5
         SITE(4,1)=0.0D0; SITE(4,2)= 0.571543301D0; SITE(4,3)=0.404151276D0
         SITE(5,1)=0.0D0; SITE(5,2)=-0.571543301D0; SITE(5,3)=0.404151276D0
C         SITE(4,1)=0.0D0; SITE(4,2)= 0.571543301D0; SITE(4,3)=0.469249307D0
C         SITE(5,1)=0.0D0; SITE(5,2)=-0.571543301D0; SITE(5,3)=0.469249307D0
         CHARGE(2)=0.241D0
         CHARGE(3)=0.241D0
         CHARGE(4)=-0.241D0
         CHARGE(5)=-0.241D0
      ENDIF

      NAT2=NATOMS/2
      ETIP=0.0D0
      V = 0.0d0
      HESS = 0.0d0

C     PRINT*,'site 1', site(1,1),site(1,2),site(1,3)
C     PRINT*,'site 2', site(2,1),site(2,2),site(2,3)
C     PRINT*,'site 3', site(3,1),site(3,2),site(3,3)
C     PRINT*,'site 4', site(4,1),site(4,2),site(4,3)

C
C  Potential energy first.
C
      DO J1=1,NAT2
         X1=X(3*(J1-1)+1)
         Y1=X(3*(J1-1)+2)
         Z1=X(3*(J1-1)+3)
         L1=X(3*(NAT2+J1-1)+1)
         M1=X(3*(NAT2+J1-1)+2)
         N1=X(3*(NAT2+J1-1)+3)
         L12=L1**2
         M12=M1**2
         N12=N1**2
         ALPHA1=SQRT(L12+M12+N12)
         RALPHA1=1.0d0/MAX(ALPHA1,1.0D-10)
         RALPHA12=1.0D0/MAX(ALPHA1**2,1.0D-20)
         CA1=COS(ALPHA1)
         C2A1=CA1
         IF (ALPHA1.LT.0.0010D0) THEN
C           C3A1=-ALPHA1/2+ALPHA1**3/24 ! bug
            C3A1=-0.5D0+ALPHA1**2/24.0D0
            S1=1.0D0-ALPHA1**2/6
         ELSE
            C3A1=(CA1-1.0D0)/ALPHA1**2
            S1=SIN(ALPHA1)/ALPHA1
         ENDIF
C        WRITE(*,'(A,6F15.5)') 'ALPHA1,RALPHA12,CA1,C2A1,C3A1,S1=',ALPHA1,RALPHA12,CA1,C2A1,C3A1,S1
         DO J2=J1+1,NAT2
            X2=X(3*(J2-1)+1)
            Y2=X(3*(J2-1)+2)
            Z2=X(3*(J2-1)+3)
            L2=X(3*(NAT2+J2-1)+1)
            M2=X(3*(NAT2+J2-1)+2)
            N2=X(3*(NAT2+J2-1)+3)
            L22=L2**2
            M22=M2**2
            N22=N2**2
            ALPHA2=SQRT(L22+M22+N22)
            RALPHA22=1.0D0/MAX(ALPHA2**2,1.0D-20)
            CA2=COS(ALPHA2)
            C2A2=CA2
            IF (ALPHA2.LT.0.00001D0) THEN
C              C3A2=-ALPHA2/2+ALPHA2**3/24
               C3A2=-0.5D0+ALPHA2**2/24.0D0
               S2=1.0D0-ALPHA2**2/6
            ELSE
               C3A2=(CA2-1.0D0)/ALPHA2**2
               S2=SIN(ALPHA2)/ALPHA2
            ENDIF
C
C  O-O LJ contribution.
C
            P1X=SITE(1,1)
            P1Y=SITE(1,2)
            P1Z=SITE(1,3)
            P2X=SITE(1,1)
            P2Y=SITE(1,2)
            P2Z=SITE(1,3)
            DIST= 
     1  Sqrt((c2a1*p1x-c3a1*l1*(l1*p1x+m1*p1y+n1*p1z)-c2a2*p2x+c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)+ 
     1      (n1*p1y - m1*p1z)*s1 - (n2*p2y - m2*p2z)*s2 + x1 - x2)**2 + 
     1   (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y+c3a2*m2*(l2*p2x+m2*p2y+n2*p2z)+ 
     1      (-(n1*p1x) + l1*p1z)*s1 - (-(n2*p2x) + l2*p2z)*s2 + y1-y2)**2 + 
     1   (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z+c3a2*n2*(l2*p2x+m2*p2y+n2*p2z)+ 
     1      (m1*p1x - l1*p1y)*s1 - (m2*p2x - l2*p2y)*s2 + z1 - z2)**2)
C           PRINT*,'coordinates of molecule pair:'
C           WRITE(*,'(3F20.10)') X1,Y1,Z1
C           WRITE(*,'(3F20.10)') L1,M1,N1
C           WRITE(*,'(3F20.10)') X2,Y2,Z2
C           WRITE(*,'(3F20.10)') L2,M2,N2
C           PRINT*,'distance=',DIST
            pairEnergy = FLJ(DIST)

            IF (GTEST.OR.SECT) THEN
               dEdRij = DFLJ(DIST)
               RDIST = 1.0D0/DIST

               firstDerivatives = 0.0d0

               do RBcoord1 = 1, 6
                  firstDerivatives(RBcoord1) = dRijdXrb(RBcoord1,
     1                 p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)                  
               enddo

! Know dRij/dx2 = -dRij/dx1, Rij/dy2 = -dRij/dy1, Rij/dz2 = -dRij/dz1

               do RBcoord1 = 10, 12
                  firstDerivatives(RBcoord1) = dRijdXrb(RBcoord1,
     1                 p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
               enddo

               pairGrad = dEdRij*firstDerivatives

! Second derivatives
               
               if (SECT) then 
                  firstDerivatives(7) = -firstDerivatives(1)
                  firstDerivatives(8) = -firstDerivatives(2)
                  firstDerivatives(9) = -firstDerivatives(3)
                
                  d2EdRij2 = D2FLJ(DIST)

                  pairHess = 0.0d0

                  do RBcoord1 = 1, 12
                     do RBcoord2 = RBcoord1, 12

                        pairHess(RBcoord1,RBcoord2) = dEdRij*
     1                  d2RijdXrb1dXrb2(RBcoord1,RBcoord2,
     1                                  p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22) +
     1                  d2EdRij2*firstDerivatives(RBcoord1)*firstDerivatives(RBcoord2)
                     enddo
                  enddo
                        
               endif
            ENDIF
C
C  Sum over charged sites.
C
            DO K1=2,NSITE
               P1X=SITE(K1,1)
               P1Y=SITE(K1,2)
               P1Z=SITE(K1,3)
               DO K2=2,NSITE
                  P2X=SITE(K2,1)
                  P2Y=SITE(K2,2)
                  P2Z=SITE(K2,3)
                  DIST=
     1  Sqrt((c2a1*p1x-c3a1*l1*(l1*p1x+m1*p1y+n1*p1z)-c2a2*p2x+c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)+ 
     1      (n1*p1y - m1*p1z)*s1 - (n2*p2y - m2*p2z)*s2 + x1 - x2)**2 + 
     1   (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y+c3a2*m2*(l2*p2x+m2*p2y+n2*p2z)+ 
     1      (-(n1*p1x) + l1*p1z)*s1 - (-(n2*p2x) + l2*p2z)*s2 + y1-y2)**2 + 
     1   (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z+c3a2*n2*(l2*p2x+m2*p2y+n2*p2z)+ 
     1      (m1*p1x - l1*p1y)*s1 - (m2*p2x - l2*p2y)*s2 + z1 - z2)**2)

                  pairEnergy = pairEnergy + HTOKJ*FC(CHARGE(K1),CHARGE(K2),DIST)

                  IF (GTEST .or. SECT) THEN
                     dEdRij=HTOKJ*DFC(CHARGE(K1),CHARGE(K2),DIST)
                     RDIST=1.0D0/DIST

                     firstDerivatives = 0.0d0

                     do RBcoord1 = 1, 6
                        firstDerivatives(RBcoord1) = dRijdXrb(RBcoord1,
     1                       p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                       C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     enddo

                     do RBcoord1 = 10, 12
                        firstDerivatives(RBcoord1) = dRijdXrb(RBcoord1,
     1                       p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                       C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     enddo                

                     pairGrad = pairGrad + dEdRij*firstDerivatives

! Second derivatives
               
                     if (SECT) then                     
                        firstDerivatives(7) = -firstDerivatives(1)
                        firstDerivatives(8) = -firstDerivatives(2)
                        firstDerivatives(9) = -firstDerivatives(3)

                        d2EdRij2 = HTOKJ*D2FC(CHARGE(K1),CHARGE(K2),DIST)
                        
                        do RBcoord1 = 1, 12
                           do RBcoord2 = RBcoord1, 12
                              pairHess(RBcoord1,RBcoord2) = pairHess(RBcoord1,RBcoord2) +
     1                             dEdRij*
     1                             d2RijdXrb1dXrb2(RBcoord1,RBcoord2,
     1                                  p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22) +
     1                             d2EdRij2*firstDerivatives(RBcoord1)*firstDerivatives(RBcoord2)
                           enddo
                        enddo                        

                     endif
                  ENDIF
               ENDDO
            ENDDO
   
            ETIP = ETIP + pairEnergy

            IF (GTEST) THEN
               do RBcoord1 = 1, 3
                  V(3*(J1-1)+RBcoord1) = V(3*(J1-1)+RBcoord1) + pairGrad(RBcoord1)              ! COM 1
                  V(3*(NAT2+J1-1)+RBcoord1) = V(3*(NAT2+J1-1)+RBcoord1) + pairGrad(RBcoord1+3)  ! alpha 1
                  V(3*(J2-1)+RBcoord1) = V(3*(J2-1)+RBcoord1) - pairGrad(RBcoord1)              ! COM 2
                  V(3*(NAT2+J2-1)+RBcoord1) = V(3*(NAT2+J2-1)+RBcoord1) + pairGrad(RBcoord1+9)  ! alpha 2
               enddo
            ENDIF

            IF (SECT) THEN
! It is necessary to convert to the slightly strange Hessian numbering scheme

               do RBcoord1 = 1, 12

                  if (RBcoord1.lt.4) then
                     index1 = 3*(J1-1)+RBcoord1
                  else if (RBcoord1.lt.7) then
                     index1 = 3*(NAT2+J1-1)+RBcoord1-3
                  else if (RBcoord1.lt.10) then
                     index1 = 3*(J2-1)+RBcoord1-6
                  else
                     index1 = 3*(NAT2+J2-1)+RBcoord1-9
                  endif

                  do RBcoord2 = RBcoord1, 12

                     if (RBcoord2.lt.4) then
                        index2 = 3*(J1-1)+RBcoord2
                     else if (RBcoord2.lt.7) then
                        index2 = 3*(NAT2+J1-1)+RBcoord2-3
                     else if (RBcoord2.lt.10) then
                        index2 = 3*(J2-1)+RBcoord2-6
                     else
                        index2 = 3*(NAT2+J2-1)+RBcoord2-9
                     endif

                     HESS(index1,index2) = HESS(index1,index2) + pairHess(RBcoord1,RBcoord2)
                  enddo
               enddo
            END IF
         ENDDO

! Applied uniform electric field along the z-direction

         if (efield.ne.0.0d0) then

            efieldInteraction = 0.0d0           
            singleGrad = 0.0d0
            singleHess = 0.0d0

            DO K1=2,NSITE
               P1X=SITE(K1,1)
               P1Y=SITE(K1,2)
               P1Z=SITE(K1,3)     ! z-coordinate of this site in the reference geometry       
               
               pDOTxzero = L1*P1X+M1*P1Y+N1*P1Z
               zReal = Z1 + P1Z*CA1 + N1*RALPHA12*pDOTxzero*(1.0d0-CA1) + (P1X*M2-P1Y*L2)*S1
               efieldInteraction = efieldInteraction + zReal*CHARGE(K1)

! Gradients if necessary

               if (GTEST) then                  
                  do RBcoord1 = 3, 6
                     singleGrad(RBcoord1) = singleGrad(RBcoord1) +
     &                 CHARGE(K1)*dzdXrb(RBcoord1,L1,M1,N1,P1X,P1Y,P1Z,CA1,S1,RALPHA12,pDOTxzero)
                  enddo
               end if

! Likewise second derivatives
! We already know most of these are zero

               if (SECT) then      
                  do RBcoord1 = 4, 6
                     do RBcoord2 = RBcoord1, 6
                        singleHess(RBcoord1,RBcoord2) = singleHess(RBcoord1,RBcoord2) + CHARGE(K1)*
     &                   d2zdXrb1dXrb2(RBcoord1,RBcoord2,L1,M1,N1,P1X,P1Y,P1Z,CA1,S1,RALPHA12,pDOTxzero)
                     enddo
                  enddo
               endif
            ENDDO

            ETIP = ETIP - eVtoKJ*efield*efieldInteraction

            if (GTEST) then
               V(3*(J1-1)+3) = V(3*(J1-1)+3) - eVtoKJ*efield*singleGrad(3)
               V(3*(NAT2+J1-1)+1) = V(3*(NAT2+J1-1)+1) - eVtoKJ*efield*singleGrad(4)
               V(3*(NAT2+J1-1)+2) = V(3*(NAT2+J1-1)+2) - eVtoKJ*efield*singleGrad(5)
               V(3*(NAT2+J1-1)+3) = V(3*(NAT2+J1-1)+3) - eVtoKJ*efield*singleGrad(6)
            endif

            if (SECT) then
               do RBcoord1 = 4, 6
                  index1 = 3*(NAT2+J1-1)+RBcoord1-3

                  do RBcoord2 = RBcoord1, 6
                     index2 = 3*(NAT2+J1-1)+RBcoord2-3                     
                     HESS(index1,index2) = HESS(index1,index2) - eVtoKJ*efield*singleHess(RBcoord1,RBcoord2)
                  enddo
               enddo

            endif
         endif

      ENDDO

      IF (SECT) THEN
! Fill in blank entries in the upper triangle of Hessian
         DO J1=1,NAT2
            DO J2=J1+1,NAT2
               HESS(3*(J2-1)+1,3*(NAT2+J1-1)+1)=HESS(3*(NAT2+J1-1)+1,3*(J2-1)+1)
               HESS(3*(J2-1)+2,3*(NAT2+J1-1)+1)=HESS(3*(NAT2+J1-1)+1,3*(J2-1)+2)
               HESS(3*(J2-1)+3,3*(NAT2+J1-1)+1)=HESS(3*(NAT2+J1-1)+1,3*(J2-1)+3)
               HESS(3*(J2-1)+1,3*(NAT2+J1-1)+2)=HESS(3*(NAT2+J1-1)+2,3*(J2-1)+1)
               HESS(3*(J2-1)+2,3*(NAT2+J1-1)+2)=HESS(3*(NAT2+J1-1)+2,3*(J2-1)+2)
               HESS(3*(J2-1)+3,3*(NAT2+J1-1)+2)=HESS(3*(NAT2+J1-1)+2,3*(J2-1)+3)
               HESS(3*(J2-1)+1,3*(NAT2+J1-1)+3)=HESS(3*(NAT2+J1-1)+3,3*(J2-1)+1)
               HESS(3*(J2-1)+2,3*(NAT2+J1-1)+3)=HESS(3*(NAT2+J1-1)+3,3*(J2-1)+2)
               HESS(3*(J2-1)+3,3*(NAT2+J1-1)+3)=HESS(3*(NAT2+J1-1)+3,3*(J2-1)+3)
            ENDDO
         ENDDO

! Symmetrise Hessian
         DO J1=1,3*NATOMS
            DO J2=J1+1,3*NATOMS
               HESS(J2,J1)=HESS(J1,J2)
            ENDDO
         ENDDO
      ENDIF

!      WRITE(*,'(A,G20.10)') 'energy=',ETIP
!      PRINT*,'coords:'
!      WRITE(*,'(I6,G20.10)') (J1,X(J1),J1=1,3*NATOMS)
C     PRINT*,'gradient:'
C     WRITE(*,'(I6,G20.10)') (J1,V(J1),J1=1,3*NATOMS)

!      DO J1=1,3*NATOMS
!         print *, j1, V(J1)
!      ENDDO

      RETURN
      END

! ****************************************************************************************************

!  SUBROUTINE to convert TIP oxygen and DV coordinates to Cartesians.

      SUBROUTINE TIPIO(X1, Y1, Z1, L1, M1, N1, COORDS)
      IMPLICIT NONE
      DOUBLE PRECISION X1, Y1, Z1, COORDS(*), C2A1, M1, L1, N1, ALPHA1, CA1, S1, C3A1, L12, M12, N12

      L12=L1**2
      M12=M1**2
      N12=N1**2
      ALPHA1=SQRT(L12+M12+N12)
      CA1=COS(ALPHA1)
      C2A1=CA1
      IF (ALPHA1.LT.0.0001D0) THEN
C        C3A1=(-ALPHA1/2+ALPHA1**3/24) ! bug
         C3A1=(-0.5D0+ALPHA1**2/24.0D0)
         S1=(1.0D0-ALPHA1**2/6)
      ELSE
         C3A1=(CA1-1.0D0)/ALPHA1**2
         S1=SIN(ALPHA1)/ALPHA1
      ENDIF

      COORDS(1) = X1
      COORDS(2) = Y1
      COORDS(3) = Z1    
      COORDS(4) = 0.756950327*c2a1 - c3a1*l1*(0.756950327*l1 - 0.585882276*n1) + 0.585882276*m1*s1 + X1
      COORDS(5) = -(c3a1*m1*(0.756950327*l1 - 0.585882276*n1)) + (-0.585882276*l1 - 0.756950327*n1)*s1 + Y1
      COORDS(6) = -0.585882276*c2a1 - c3a1*(0.756950327*l1 - 0.585882276*n1)*n1 + 0.756950327*m1*s1 + Z1
      COORDS(7) = -0.756950327*c2a1 + c3a1*l1*(0.756950327*l1 + 0.585882276*n1) + 0.585882276*m1*s1 + X1
      COORDS(8) = c3a1*m1*(0.756950327*l1 + 0.585882276*n1) + (-0.585882276*l1 + 0.756950327*n1)*s1 + Y1
      COORDS(9) = -0.585882276*c2a1 + c3a1*(0.756950327*l1 + 0.585882276*n1)*n1 - 0.756950327*m1*s1 + Z1

      RETURN
      END

! ****************************************************************************************************
