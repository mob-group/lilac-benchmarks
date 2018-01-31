C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
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

C
C  Energy and gradient for rigid body molecule - 6-site virus pentamer.
C  Variables have been reordered to the OPTIM convention.
C
      SUBROUTINE FCAPSID(NATOMS,X,V,ECAPSID,GTEST,SECT)
      USE MODHESS
      IMPLICIT NONE
      LOGICAL GTEST,SECT
      INTEGER J1, J2, K1, K2, NATOMS, NAT2
      DOUBLE PRECISION X(3*NATOMS), V(3*NATOMS), X1, X2, Y1, Y2, Z1, Z2, FATT, DFATT, FREP, DFREP,
     1                 ECAPSID, DUMMY, EPS2, RALPHA12, RALPHA22, RDIST, SITE(6,3), DDFATT, RHO,
     2                 M1, L1, N1, M2, L2, N2, ALPHA1, RAD, CA1, CA2, ALPHA2, S1, S2, C3A1, C3A2,
     3                 L12, M12, N12, L22, M22, N22, DDFREP, DUMMY3, 
     4                 GX1, GY1, GZ1, GL1, GM1, GN1, GL2, GM2, GN2, EPSEFF, C2A2, C2A1,
     5                 D1, D2, D3, D4, D5, D6, D7, D8, D9, DA, DB, DC, DIST, DUMMY2,
     6                 P1X, P1Y, P1Z, P2X, P2Y, P2Z, XDUMM, SIGMA, HEIGHT,
     1           D1S, D2S, D3S, D4S, D5S, D6S, D7S, D8S, D9S, DAS, DBS, DCS,
     7           DD11, DD12, DD13, DD14, DD15, DD16, DD17, DD18, DD19, DD1A, DD1B, DD1C,
     7           DD22, DD23, DD24, DD25, DD26, DD27, DD28, DD29, DD2A, DD2B, DD2C,
     7           DD33, DD34, DD35, DD36, DD37, DD38, DD39, DD3A, DD3B, DD3C,
     7           DD44, DD45, DD46, DD47, DD48, DD49, DD4A, DD4B, DD4C,
     7           DD55, DD56, DD57, DD58, DD59, DD5A, DD5B, DD5C,
     7           DD66, DD67, DD68, DD69, DD6A, DD6B, DD6C,
     7           DD77, DD78, DD79, DD7A, DD7B, DD7C,
     7           DD88, DD89, DD8A, DD8B, DD8C,
     7           DD99, DD9A, DD9B, DD9C,
     7           DDAA, DDAB, DDAC,
     7           DDBB, DDBC,
     7           DDCC,
     7           S11, S12, S13, S14, S15, S16, S17, S18, S19, S1A, S1B, S1C,
     7           S22, S23, S24, S25, S26, S27, S28, S29, S2A, S2B, S2C,
     7           S33, S34, S35, S36, S37, S38, S39, S3A, S3B, S3C,
     7           S44, S45, S46, S47, S48, S49, S4A, S4B, S4C,
     7           S55, S56, S57, S58, S59, S5A, S5B, S5C,
     7           S66, S67, S68, S69, S6A, S6B, S6C,
     7           S77, S78, S79, S7A, S7B, S7C,
     7           S88, S89, S8A, S8B, S8C,
     7           S99, S9A, S9B, S9C,
     7           SAA, SAB, SAC,
     7           SBB, SBC,
     7           SCC
      LOGICAL EVAP
      COMMON /EV/ EVAP
      COMMON /CAPS/ RHO, EPS2, RAD, HEIGHT
C
C  Attractive and repulsive function definitions and their derivatives.
C
      FATT(RHO,XDUMM)=-1.0D0 + (1.0D0 - EXP(RHO*(1.0D0 - XDUMM)))**2
      DFATT(RHO,XDUMM)=2.0D0*(-EXP(2.0D0*RHO*(1.0D0-XDUMM)) + EXP(RHO*(1.0D0-XDUMM)))*RHO
      DDFATT(RHO,XDUMM)=-2.0D0*(-2.0D0*EXP(2.0D0*RHO*(1.0D0-XDUMM)) + EXP(RHO*(1.0D0-XDUMM)))*RHO**2
      FREP(SIGMA,XDUMM)=(SIGMA/XDUMM)**12
      DFREP(SIGMA,XDUMM)=-12.0D0*(SIGMA/XDUMM)**12/XDUMM
      DDFREP(SIGMA,XDUMM)=156.0D0*(SIGMA/XDUMM)**12/XDUMM**2
C
C Derivatives of R(i,j) with respect to rigid body coordinates.
C P(i) = (p1x,p1y,p1z), P(j)=(p2x,p2y,p2z) are the site cordinates in
C the reference geometry for the molecule at the origin.
C
C      D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1    rdist*(c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
C     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - n2*p2y*s2 + 
C     1    m2*p2z*s2 + x1 - x2)

C      D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1    rdist*(c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + n2*p2x*s2 - 
C     1    l2*p2z*s2 + y1 - y2)

C      D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1    rdist*(c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - m2*p2x*s2 + 
C     1    l2*p2y*s2 + z1 - z2)

C      D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1      rdist*((c3a1*(2*l1*p1x*(-1 + l12*ralpha12) + 
C     1            (m1*p1y + n1*p1z)*(-1 + 2*l12*ralpha12)) + 
C     1         (l1*l12*p1x + l12*m1*p1y - l1*n1*p1y + l1*m1*p1z + l12*n1*p1z)*
C     1          ralpha12*s1 + l1*(c2a1*(n1*p1y - m1*p1z)*ralpha12 - p1x*s1))*
C     1       (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
C     1         c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
C     1         n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + 
C     1      (c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1         c3a1*m1*(-p1x + 2*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
C     1         (-(l1*p1y) + p1z + l1* (l1*m1*p1x + n1*p1x + m1**2*p1y - l1*p1z + m1*n1*p1z)*ralpha12)*
C     1          s1)*(c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1         c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
C     1         n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + 
C     1      (-(c2a1*l1*(-(m1*p1x) + l1*p1y)*ralpha12) + 
C     1         c3a1*n1*(-p1x + 2*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
C     1         (-p1y - l1*p1z + l1*(-(m1*p1x) + l1*n1*p1x + l1*p1y + m1*n1*p1y + n1**2*p1z)*ralpha12)*s1)*
C     1       (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1         c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - 
C     1         m2*p2x*s2 + l2*p2y*s2 + z1 - z2))

C      D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1 rdist*((-(c2a1*m1*(-(n1*p1y) + m1*p1z)*ralpha12) + 
C     1         c3a1*l1*(-p1y + 2*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
C     1         (-(m1*p1x) - p1z + m1* (l1**2*p1x + l1*m1*p1y - n1*p1y + m1*p1z + l1*n1*p1z)*ralpha12)*
C     1          s1)*(c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
C     1         c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
C     1         n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + (-(c3a1*(l1*p1x + 2*m1*p1y + n1*p1z)) + 
C     1         c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
C     1         (-(m1*p1y) + (l1*m12*p1x + m1*n1*p1x + m1*m12*p1y - l1*m1*p1z + 
C     1               m12*n1*p1z)*ralpha12)*s1)*
C     1       (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1         c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
C     1         n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + (c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1         c3a1*n1*(-p1y + 2*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
C     1         (p1x - m1*(m1 - l1*n1)*p1x*ralpha12 + 
C     1            m1*((l1 + m1*n1)*p1y*ralpha12 + p1z*(-1 + n1**2*ralpha12)))*s1)*
C     1       (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1         c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - 
C     1         m2*p2x*s2 + l2*p2y*s2 + z1 - z2))

C      D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1 rdist*((c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1         c3a1*l1*(-p1z + 2*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + (-(n1*p1x) + p1y + n1*
C     1             (l1**2*p1x + l1*m1*p1y - n1*p1y + m1*p1z + l1*n1*p1z)*ralpha12)*
C     1          s1)*(c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
C     1         c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
C     1         n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + 
C     1      (-(c2a1*n1*(n1*p1x - l1*p1z)*ralpha12) + 
C     1         c3a1*m1*(-p1z + 2*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
C     1         (-p1x - n1*p1y + n1*(l1*m1*p1x + n1*p1x + m1**2*p1y - l1*p1z + 
C     1               m1*n1*p1z)*ralpha12)*s1)*
C     1       (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1         c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
C     1         n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + 
C     1      (-(c3a1*(l1*p1x + m1*p1y + 2*n1*p1z)) + 
C     1         c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1         2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
C     1         (-(n1*p1z) + (-(m1*n1*p1x) + l1*n12*p1x + l1*n1*p1y + m1*n12*p1y + 
C     1               n1*n12*p1z)*ralpha12)*s1)*
C     1       (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1         c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - 
C     1         m2*p2x*s2 + l2*p2y*s2 + z1 - z2))

C      D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1   -rdist*(c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
C     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - n2*p2y*s2 + 
C     1    m2*p2z*s2 + x1 - x2)

C      D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1   -rdist*(c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + n2*p2x*s2 - 
C     1    l2*p2z*s2 + y1 - y2)

C      D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1   -rdist*(c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - m2*p2x*s2 + 
C     1    l2*p2y*s2 + z1 - z2)

C      DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1         rdist*((c3a2*(2*l2*p2x + m2*p2y + n2*p2z - 
C     1            2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) - 
C     1         (l2*l22*p2x + l22*m2*p2y - l2*n2*p2y + l2*m2*p2z + l22*n2*p2z)*
C     1          ralpha22*s2 + l2*(-(c2a2*n2*p2y*ralpha22) + c2a2*m2*p2z*ralpha22 + 
C     1            p2x*s2))*(c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - 
C     1         c2a2*p2x + c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + 
C     1         (n1*p1y - m1*p1z)*s1 - n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + 
C     1      (-(c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22) + 
C     1         c3a2*m2*(p2x - 2*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) - (-(l2*p2y) + p2z + l2*
C     1             (l2*m2*p2x + n2*p2x + m2**2*p2y - l2*p2z + m2*n2*p2z)*ralpha22)*
C     1          s2)*(c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1         c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
C     1         n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + 
C     1      (c2a2*l2*(-(m2*p2x) + l2*p2y)*ralpha22 + 
C     1         c3a2*n2*(p2x - 2*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
C     1    (p2y + l2*p2z - l2*(-(m2*p2x) + l2*n2*p2x + l2*p2y + m2*n2*p2y + n2**2*p2z)*ralpha22)*s2)*
C     1       (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - m2*p2x*s2 + l2*p2y*s2 + z1 - z2))

C      DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1 rdist*((c2a2*m2*(-(n2*p2y) + m2*p2z)*ralpha22 + 
C     1         c3a2*l2*(p2y - 2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
C     1         (m2*p2x + p2z - m2*(l2**2*p2x + l2*m2*p2y - n2*p2y + m2*p2z + l2*n2*p2z)*ralpha22)*s2)*
C     1       (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
C     1         c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
C     1         n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + (c3a2*(l2*p2x + 2*m2*p2y + n2*p2z - 
C     1            2*m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) - 
C     1         (l2*m22*p2x + m2*n2*p2x + m2*m22*p2y - l2*m2*p2z + m22*n2*p2z)*
C     1          ralpha22*s2 + m2*(c2a2*(n2*p2x - l2*p2z)*ralpha22 + p2y*s2))*
C     1       (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1         c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
C     1         n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + (-(c2a2*m2*(m2*p2x - l2*p2y)*ralpha22) + 
C     1         c3a2*n2*(p2y - 2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
C     1         (-p2x + m2*p2z - m2*(-(m2*p2x) + l2*n2*p2x + l2*p2y + m2*n2*p2y + 
C     1               n2**2*p2z)*ralpha22)*s2)*
C     1       (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1         c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - 
C     1         m2*p2x*s2 + l2*p2y*s2 + z1 - z2))

C      DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1 rdist*((-(c2a2*n2*(n2*p2y - m2*p2z)*ralpha22) + 
C     1         c3a2*l2*(p2z - 2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
C     1         (n2*p2x - p2y - n2*(l2**2*p2x + l2*m2*p2y - n2*p2y + m2*p2z + l2*n2*p2z)*ralpha22)*s2)*
C     1       (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
C     1         c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
C     1         n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + (c2a2*n2*(n2*p2x - l2*p2z)*ralpha22 + 
C     1         c3a2*m2*(p2z - 2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
C     1         (p2x + n2*p2y - n2*(l2*m2*p2x + n2*p2x + m2**2*p2y - l2*p2z + m2*n2*p2z)*ralpha22)*s2)*
C     1       (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1         c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
C     1         n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + (c3a2*(l2*p2x + m2*p2y + 2*n2*p2z - 
C     1            2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
C     1         (n2*(m2*p2x - l2*p2y) - n22*(l2*p2x + m2*p2y + n2*p2z))*ralpha22*
C     1          s2 + n2*(-(c2a2*m2*p2x*ralpha22) + c2a2*l2*p2y*ralpha22 + p2z*s2))*
C     1       (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1         c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - 
C     1         m2*p2x*s2 + l2*p2y*s2 + z1 - z2))
C 
C  Second derivatives of R(i,j) with respect to rigid body coordinates.
C 

C      DD11(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = rdist
C     1 -D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD12(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1 -D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD13(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1 -D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD14(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1  rdist*(-(c3a1*l1*p1x) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1 c2a1*l1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1 2*c3a1*l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1x*s1 - 
C     1 l1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1 l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)
C     1 -D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD15(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1 rdist*(-(c3a1*l1*p1y) + c2a1*m1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1 2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1x*s1 - p1z*s1 - 
C     1 m1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1 l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)
C     1 -D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD16(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1 rdist*(-(c3a1*l1*p1z) + c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1 2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1x*s1 + p1y*s1 - 
C     1 n1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1 l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)
C     1 -D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD17(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = -rdist
C     1 -D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD18(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1 -D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD19(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1 -D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD1A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = rdist*(c3a2*l2*p2x + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1 c2a2*l2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1 2*c3a2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2x*s2 + 
C     1 l2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1 l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)
C     1 -D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD1B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = rdist*(c3a2*l2*p2y - c2a2*m2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1 2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2x*s2 + p2z*s2 + 
C     1 m2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1 l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)
C     1 -D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD1C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = rdist*(c3a2*l2*p2z - c2a2*n2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1 2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2x*s2 - p2y*s2 + 
C     1 n2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1 l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)
C     1 -D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD22(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = rdist
C     1 -D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD23(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1 -D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD24(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = rdist*(-(c3a1*m1*p1x) + c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1 2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1y*s1 + p1z*s1 - 
C     1 l1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1 l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)
C     1 -D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD25(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = rdist*(-(c3a1*m1*p1y) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1 c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1 2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1y*s1 - 
C     1 m1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1 m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)
C     1 -D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD26(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = rdist*(-(c3a1*m1*p1z) + c2a1*n1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1 2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 - n1*p1y*s1 - 
C     1 n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1 m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)
C     1 -D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD27(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1 -D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD28(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = -rdist
C     1 -D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD29(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 0
C     1 -D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD2A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = rdist*(c3a2*m2*p2x - c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1 2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2y*s2 - p2z*s2 + 
C     1 l2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1 l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)
C     1 -D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD2B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = rdist*(c3a2*m2*p2y + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1 c2a2*m2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1 2*c3a2*m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2y*s2 + 
C     1 m2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1 m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)
C     1 -D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD2C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = rdist*(c3a2*m2*p2z - c2a2*n2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1 2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 + n2*p2y*s2 + 
C     1 n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1 m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)
C     1 -D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD33(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = rdist
C     1 -D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD34(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = rdist*(-(c3a1*n1*p1x) + c2a1*l1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1 2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 - l1*p1z*s1 - 
C     1 l1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1 l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)
C     1 -D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD35(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = rdist*(-(c3a1*n1*p1y) + c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1 2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + p1x*s1 - m1*p1z*s1 - 
C     1 m1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1 m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)
C     1 -D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD36(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = rdist*(-(c3a1*n1*p1z) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1 c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1 2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1z*s1 - 
C     1 n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1 n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)
C     1 -D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD37(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1 -D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD38(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1 -D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD39(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = -rdist
C     1 -D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD3A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = rdist*(c3a2*n2*p2x - c2a2*l2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1 2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 + l2*p2z*s2 + 
C     1 l2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1 l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)
C     1 -D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD3B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = rdist*(c3a2*n2*p2y - c2a2*m2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1 2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - p2x*s2 + m2*p2z*s2 + 
C     1 m2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1 m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)
C     1 -D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD3C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = rdist*(c3a2*n2*p2z + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1 c2a2*n2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1 2*c3a2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2z*s2 + 
C     1 n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1 n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)
C     1 -D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD44(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(-(c3a1*l1*p1x) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1     c2a1*l1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1     2*c3a1*l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1x*s1 - 
C     1     l1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1     l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)**2 + 
C     1 2*(-(c3a1*m1*p1x) + c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1     2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1y*s1 + 
C     1     p1z*s1 - l1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1     l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)**2 + 
C     1 2*(-(c3a1*n1*p1x) + c2a1*l1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1     2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 - 
C     1     l1*p1z*s1 - l1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1     l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)**2 + 
C     1 2*(-2*c3a1*p1x - (3*c2a1*l12*(n1*p1y - m1*p1z))/alpha1**4 + 
C     1    (c2a1*l1**3*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    (8*c3a1*l1**3*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    c2a1*l12*p1x*ralpha12 + 4*c3a1*l12*p1x*ralpha12 + 
C     1    c2a1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1    6*c3a1*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 + 
C     1    (3*l12*(n1*p1y - m1*p1z)*s1)/alpha1**4 - 
C     1    (5*l1**3*(l1*p1x + m1*p1y + n1*p1z)*s1)/alpha1**4 + 
C     1    3*l12*p1x*ralpha12*s1 - (n1*p1y - m1*p1z)*ralpha12*s1 - 
C     1    l12*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    3*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
C     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
C     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
C     1 2*((-3*c2a1*l12*(-(n1*p1x) + l1*p1z))/alpha1**4 + 
C     1    (c2a1*l12*m1*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    (8*c3a1*l12*m1*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 + 
C     1    4*c3a1*l1*m1*p1x*ralpha12 - c2a1*l12*p1y*ralpha12 + 
C     1    2*c2a1*l1*p1z*ralpha12 + c2a1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1    2*c3a1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 + 
C     1    (3*l12*(-(n1*p1x) + l1*p1z)*s1)/alpha1**4 - 
C     1    (5*l12*m1*(l1*p1x + m1*p1y + n1*p1z)*s1)/alpha1**4 + 
C     1    2*l1*m1*p1x*ralpha12*s1 + l12*p1y*ralpha12*s1 - 
C     1    2*l1*p1z*ralpha12*s1 - (-(n1*p1x) + l1*p1z)*ralpha12*s1 - 
C     1    l12*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
C     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
C     1 2*((-3*c2a1*l12*(m1*p1x - l1*p1y))/alpha1**4 + 
C     1    (c2a1*l12*n1*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    (8*c3a1*l12*n1*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 + 
C     1    4*c3a1*l1*n1*p1x*ralpha12 - 2*c2a1*l1*p1y*ralpha12 + 
C     1    c2a1*(m1*p1x - l1*p1y)*ralpha12 - c2a1*l12*p1z*ralpha12 + 
C     1    2*c3a1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
C     1    (3*l12*(m1*p1x - l1*p1y)*s1)/alpha1**4 - p1z*s1 - 
C     1    (5*l12*n1*(l1*p1x + m1*p1y + n1*p1z)*s1)/alpha1**4 + 
C     1    2*l1*n1*p1x*ralpha12*s1 + 2*l1*p1y*ralpha12*s1 - 
C     1    (m1*p1x - l1*p1y)*ralpha12*s1 - 
C     1    l12*(m1*p1x - l1*p1y)*ralpha12*s1 + l12*p1z*ralpha12*s1 + 
C     1    n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
C     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
C     1 -D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD45(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(-(c3a1*l1*p1x) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1    c2a1*l1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1    2*c3a1*l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1x*s1 - 
C     1    l1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (-(c3a1*l1*p1y) + c2a1*m1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1x*s1 - 
C     1    p1z*s1 - m1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) + 
C     1 2*(-(c3a1*m1*p1x) + c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1y*s1 + 
C     1    p1z*s1 - l1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (-(c3a1*m1*p1y) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1    c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1    2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1y*s1 - 
C     1    m1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) + 
C     1 2*(-(c3a1*n1*p1x) + c2a1*l1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 - 
C     1    l1*p1z*s1 - l1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (-(c3a1*n1*p1y) + c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + p1x*s1 - 
C     1    m1*p1z*s1 - m1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) + 
C     1 2*(-(c3a1*p1y) - (3*c2a1*l1*m1*(n1*p1y - m1*p1z))/alpha1**4 + 
C     1    (c2a1*l12*m1*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    (8*c3a1*l12*m1*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    c2a1*l1*m1*p1x*ralpha12 + 2*c3a1*l1*m1*p1x*ralpha12 + 
C     1    2*c3a1*l12*p1y*ralpha12 - c2a1*l1*p1z*ralpha12 + 
C     1    2*c3a1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
C     1    (3*l1*m1*(n1*p1y - m1*p1z)*s1)/alpha1**4 - 
C     1    (5*l12*m1*(l1*p1x + m1*p1y + n1*p1z)*s1)/alpha1**4 + 
C     1    2*l1*m1*p1x*ralpha12*s1 + l12*p1y*ralpha12*s1 + 
C     1    l1*p1z*ralpha12*s1 - l1*m1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
C     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
C     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
C     1 2*(-(c3a1*p1x) - (3*c2a1*l1*m1*(-(n1*p1x) + l1*p1z))/alpha1**4 + 
C     1    (c2a1*l1*m12*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    (8*c3a1*l1*m12*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 + 
C     1    2*c3a1*m12*p1x*ralpha12 - c2a1*l1*m1*p1y*ralpha12 + 
C     1    2*c3a1*l1*m1*p1y*ralpha12 + c2a1*m1*p1z*ralpha12 + 
C     1    2*c3a1*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
C     1    (3*l1*m1*(-(n1*p1x) + l1*p1z)*s1)/alpha1**4 - 
C     1    (5*l1*m12*(l1*p1x + m1*p1y + n1*p1z)*s1)/alpha1**4 + 
C     1    m12*p1x*ralpha12*s1 + 2*l1*m1*p1y*ralpha12*s1 - 
C     1    m1*p1z*ralpha12*s1 - l1*m1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
C     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
C     1 2*((-3*c2a1*l1*m1*(m1*p1x - l1*p1y))/alpha1**4 + 
C     1    (c2a1*l1*m1*n1*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    (8*c3a1*l1*m1*n1*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 + 
C     1    c2a1*l1*p1x*ralpha12 + 2*c3a1*m1*n1*p1x*ralpha12 - 
C     1    c2a1*m1*p1y*ralpha12 + 2*c3a1*l1*n1*p1y*ralpha12 - 
C     1    c2a1*l1*m1*p1z*ralpha12 + 
C     1    (3*l1*m1*(m1*p1x - l1*p1y)*s1)/alpha1**4 - 
C     1    (5*l1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*s1)/alpha1**4 - 
C     1    l1*p1x*ralpha12*s1 + m1*n1*p1x*ralpha12*s1 + m1*p1y*ralpha12*s1 + 
C     1    l1*n1*p1y*ralpha12*s1 - l1*m1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1    l1*m1*p1z*ralpha12*s1)*
C     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
C     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
C     1 -D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*rdist

C      DD46(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(-(c3a1*l1*p1x) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1    c2a1*l1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1    2*c3a1*l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1x*s1 - 
C     1    l1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (-(c3a1*l1*p1z) + c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1x*s1 + 
C     1    p1y*s1 - n1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) + 
C     1 2*(-(c3a1*m1*p1x) + c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1y*s1 + 
C     1    p1z*s1 - l1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (-(c3a1*m1*p1z) + c2a1*n1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 - 
C     1    n1*p1y*s1 - n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) + 
C     1 2*(-(c3a1*n1*p1x) + c2a1*l1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 - 
C     1    l1*p1z*s1 - l1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (-(c3a1*n1*p1z) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1    c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1    2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1z*s1 - 
C     1    n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1    n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) + 
C     1 2*(-(c3a1*p1z) - (3*c2a1*l1*n1*(n1*p1y - m1*p1z))/alpha1**4 + 
C     1    (c2a1*l12*n1*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    (8*c3a1*l12*n1*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    c2a1*l1*n1*p1x*ralpha12 + 2*c3a1*l1*n1*p1x*ralpha12 + 
C     1    c2a1*l1*p1y*ralpha12 + 2*c3a1*l12*p1z*ralpha12 + 
C     1    2*c3a1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
C     1    (3*l1*n1*(n1*p1y - m1*p1z)*s1)/alpha1**4 - 
C     1    (5*l12*n1*(l1*p1x + m1*p1y + n1*p1z)*s1)/alpha1**4 + 
C     1    2*l1*n1*p1x*ralpha12*s1 - l1*p1y*ralpha12*s1 + 
C     1    l12*p1z*ralpha12*s1 - l1*n1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
C     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
C     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
C     1 2*((-3*c2a1*l1*n1*(-(n1*p1x) + l1*p1z))/alpha1**4 + 
C     1    (c2a1*l1*m1*n1*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    (8*c3a1*l1*m1*n1*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    c2a1*l1*p1x*ralpha12 + 2*c3a1*m1*n1*p1x*ralpha12 - 
C     1    c2a1*l1*n1*p1y*ralpha12 + 2*c3a1*l1*m1*p1z*ralpha12 + 
C     1    c2a1*n1*p1z*ralpha12 + 
C     1    (3*l1*n1*(-(n1*p1x) + l1*p1z)*s1)/alpha1**4 - 
C     1    (5*l1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*s1)/alpha1**4 + 
C     1    l1*p1x*ralpha12*s1 + m1*n1*p1x*ralpha12*s1 + 
C     1    l1*n1*p1y*ralpha12*s1 + l1*m1*p1z*ralpha12*s1 - 
C     1    n1*p1z*ralpha12*s1 - l1*n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1)*
C     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
C     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
C     1 2*(-(c3a1*p1x) - (3*c2a1*l1*n1*(m1*p1x - l1*p1y))/alpha1**4 + 
C     1    (c2a1*l1*n12*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    (8*c3a1*l1*n12*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 + 
C     1    2*c3a1*n12*p1x*ralpha12 - c2a1*n1*p1y*ralpha12 - 
C     1    c2a1*l1*n1*p1z*ralpha12 + 2*c3a1*l1*n1*p1z*ralpha12 + 
C     1    2*c3a1*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
C     1    (3*l1*n1*(m1*p1x - l1*p1y)*s1)/alpha1**4 - 
C     1    (5*l1*n12*(l1*p1x + m1*p1y + n1*p1z)*s1)/alpha1**4 + 
C     1    n12*p1x*ralpha12*s1 + n1*p1y*ralpha12*s1 - 
C     1    l1*n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 2*l1*n1*p1z*ralpha12*s1 + 
C     1    l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
C     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
C     1 -D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD47(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = -(rdist*(-(c3a1*l1*p1x) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1 c2a1*l1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1 2*c3a1*l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1x*s1 - 
C     1 l1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1 l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1))
C     1 -D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD48(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = -(rdist*(-(c3a1*m1*p1x) + c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1 2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1y*s1 + 
C     1 p1z*s1 - l1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1 l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1))
C     1 -D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD49(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = -(rdist*(-(c3a1*n1*p1x) + c2a1*l1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1 2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 - 
C     1 l1*p1z*s1 - l1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1 l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1))
C     1 -D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD4A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(-(c3a1*l1*p1x) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1    c2a1*l1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1    2*c3a1*l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1x*s1 - 
C     1    l1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*l2*p2x + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1    c2a2*l2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1    2*c3a2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2x*s2 + 
C     1    l2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(-(c3a1*m1*p1x) + c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1y*s1 + 
C     1    p1z*s1 - l1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*m2*p2x - c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2y*s2 - 
C     1    p2z*s2 + l2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(-(c3a1*n1*p1x) + c2a1*l1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 - 
C     1    l1*p1z*s1 - l1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*n2*p2x - c2a2*l2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 + 
C     1    l2*p2z*s2 + l2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)))/2.
C     1 -D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD4B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(-(c3a1*l1*p1x) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1    c2a1*l1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1    2*c3a1*l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1x*s1 - 
C     1    l1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*l2*p2y - c2a2*m2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2x*s2 + 
C     1    p2z*s2 + m2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(-(c3a1*m1*p1x) + c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1y*s1 + 
C     1    p1z*s1 - l1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*m2*p2y + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1    c2a2*m2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1    2*c3a2*m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2y*s2 + 
C     1    m2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(-(c3a1*n1*p1x) + c2a1*l1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 - 
C     1    l1*p1z*s1 - l1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*n2*p2y - c2a2*m2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - p2x*s2 + 
C     1    m2*p2z*s2 + m2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)))/2.
C     1 -D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD4C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(-(c3a1*l1*p1x) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1    c2a1*l1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1    2*c3a1*l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1x*s1 - 
C     1    l1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*l2*p2z - c2a2*n2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2x*s2 - 
C     1    p2y*s2 + n2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(-(c3a1*m1*p1x) + c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1y*s1 + 
C     1    p1z*s1 - l1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*m2*p2z - c2a2*n2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 + 
C     1    n2*p2y*s2 + n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(-(c3a1*n1*p1x) + c2a1*l1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 - 
C     1    l1*p1z*s1 - l1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*n2*p2z + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1    c2a2*n2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1    2*c3a2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2z*s2 + 
C     1    n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1    n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)))/2.
C     1 -D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD55(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(-(c3a1*l1*p1y) + 
C     1     c2a1*m1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1     2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1x*s1 - 
C     1     p1z*s1 - m1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1     l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)**2 + 
C     1 2*(-(c3a1*m1*p1y) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1     c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1     2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1y*s1 - 
C     1     m1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1     m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)**2 + 
C     1 2*(-(c3a1*n1*p1y) + c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1     2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + p1x*s1 - 
C     1     m1*p1z*s1 - m1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1     m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)**2 + 
C     1 2*((-3*c2a1*m12*(n1*p1y - m1*p1z))/alpha1**4 + 
C     1    (c2a1*l1*m12*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    (8*c3a1*l1*m12*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    c2a1*m12*p1x*ralpha12 + 4*c3a1*l1*m1*p1y*ralpha12 - 
C     1    2*c2a1*m1*p1z*ralpha12 + c2a1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1    2*c3a1*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 + 
C     1    (3*m12*(n1*p1y - m1*p1z)*s1)/alpha1**4 - 
C     1    (5*l1*m12*(l1*p1x + m1*p1y + n1*p1z)*s1)/alpha1**4 + 
C     1    m12*p1x*ralpha12*s1 + 2*l1*m1*p1y*ralpha12*s1 + 
C     1    2*m1*p1z*ralpha12*s1 - (n1*p1y - m1*p1z)*ralpha12*s1 - 
C     1    m12*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
C     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
C     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
C     1 2*(-2*c3a1*p1y - (3*c2a1*m12*(-(n1*p1x) + l1*p1z))/alpha1**4 + 
C     1    (c2a1*m1**3*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    (8*c3a1*m1**3*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    c2a1*m12*p1y*ralpha12 + 4*c3a1*m12*p1y*ralpha12 + 
C     1    c2a1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1    6*c3a1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 + 
C     1    (3*m12*(-(n1*p1x) + l1*p1z)*s1)/alpha1**4 - 
C     1    (5*m1**3*(l1*p1x + m1*p1y + n1*p1z)*s1)/alpha1**4 + 
C     1    3*m12*p1y*ralpha12*s1 - (-(n1*p1x) + l1*p1z)*ralpha12*s1 - 
C     1    m12*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    3*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
C     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
C     1 2*((-3*c2a1*m12*(m1*p1x - l1*p1y))/alpha1**4 + 
C     1    (c2a1*m12*n1*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    (8*c3a1*m12*n1*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 + 
C     1    2*c2a1*m1*p1x*ralpha12 + 4*c3a1*m1*n1*p1y*ralpha12 + 
C     1    c2a1*(m1*p1x - l1*p1y)*ralpha12 - c2a1*m12*p1z*ralpha12 + 
C     1    2*c3a1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
C     1    (3*m12*(m1*p1x - l1*p1y)*s1)/alpha1**4 - p1z*s1 - 
C     1    (5*m12*n1*(l1*p1x + m1*p1y + n1*p1z)*s1)/alpha1**4 - 
C     1    2*m1*p1x*ralpha12*s1 + 2*m1*n1*p1y*ralpha12*s1 - 
C     1    (m1*p1x - l1*p1y)*ralpha12*s1 - 
C     1    m12*(m1*p1x - l1*p1y)*ralpha12*s1 + m12*p1z*ralpha12*s1 + 
C     1    n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
C     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
C     1 -D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD56(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(-(c3a1*l1*p1y) + 
C     1    c2a1*m1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1x*s1 - 
C     1    p1z*s1 - m1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (-(c3a1*l1*p1z) + c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1x*s1 + 
C     1    p1y*s1 - n1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) + 
C     1 2*(-(c3a1*m1*p1y) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1    c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1    2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1y*s1 - 
C     1    m1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (-(c3a1*m1*p1z) + c2a1*n1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 - 
C     1    n1*p1y*s1 - n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) + 
C     1 2*(-(c3a1*n1*p1y) + c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + p1x*s1 - 
C     1    m1*p1z*s1 - m1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (-(c3a1*n1*p1z) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1    c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1    2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1z*s1 - 
C     1    n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1    n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) + 
C     1 2*((-3*c2a1*m1*n1*(n1*p1y - m1*p1z))/alpha1**4 + 
C     1    (c2a1*l1*m1*n1*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    (8*c3a1*l1*m1*n1*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    c2a1*m1*n1*p1x*ralpha12 + c2a1*m1*p1y*ralpha12 + 
C     1    2*c3a1*l1*n1*p1y*ralpha12 + 2*c3a1*l1*m1*p1z*ralpha12 - 
C     1    c2a1*n1*p1z*ralpha12 + 
C     1    (3*m1*n1*(n1*p1y - m1*p1z)*s1)/alpha1**4 - 
C     1    (5*l1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*s1)/alpha1**4 + 
C     1    m1*n1*p1x*ralpha12*s1 - m1*p1y*ralpha12*s1 + 
C     1    l1*n1*p1y*ralpha12*s1 + l1*m1*p1z*ralpha12*s1 + 
C     1    n1*p1z*ralpha12*s1 - m1*n1*(n1*p1y - m1*p1z)*ralpha12*s1)*
C     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
C     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
C     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
C     1 2*(-(c3a1*p1z) - (3*c2a1*m1*n1*(-(n1*p1x) + l1*p1z))/alpha1**4 + 
C     1    (c2a1*m12*n1*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    (8*c3a1*m12*n1*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    c2a1*m1*p1x*ralpha12 - c2a1*m1*n1*p1y*ralpha12 + 
C     1    2*c3a1*m1*n1*p1y*ralpha12 + 2*c3a1*m12*p1z*ralpha12 + 
C     1    2*c3a1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
C     1    (3*m1*n1*(-(n1*p1x) + l1*p1z)*s1)/alpha1**4 - 
C     1    (5*m12*n1*(l1*p1x + m1*p1y + n1*p1z)*s1)/alpha1**4 + 
C     1    m1*p1x*ralpha12*s1 + 2*m1*n1*p1y*ralpha12*s1 + 
C     1    m12*p1z*ralpha12*s1 - m1*n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
C     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
C     1 2*(-(c3a1*p1y) - (3*c2a1*m1*n1*(m1*p1x - l1*p1y))/alpha1**4 + 
C     1    (c2a1*m1*n12*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    (8*c3a1*m1*n12*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 + 
C     1    c2a1*n1*p1x*ralpha12 + 2*c3a1*n12*p1y*ralpha12 - 
C     1    c2a1*m1*n1*p1z*ralpha12 + 2*c3a1*m1*n1*p1z*ralpha12 + 
C     1    2*c3a1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
C     1    (3*m1*n1*(m1*p1x - l1*p1y)*s1)/alpha1**4 - 
C     1    (5*m1*n12*(l1*p1x + m1*p1y + n1*p1z)*s1)/alpha1**4 - 
C     1    n1*p1x*ralpha12*s1 + n12*p1y*ralpha12*s1 - 
C     1    m1*n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 2*m1*n1*p1z*ralpha12*s1 + 
C     1    m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
C     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
C     1 -D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD57(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = -(rdist*(-(c3a1*l1*p1y) + c2a1*m1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1 2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1x*s1 - 
C     1 p1z*s1 - m1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1 l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1))
C     1 -D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD58(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = -(rdist*(-(c3a1*m1*p1y) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1 c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1 2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1y*s1 - 
C     1 m1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1 m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1))
C     1 -D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD59(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = -(rdist*(-(c3a1*n1*p1y) + c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1 2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + p1x*s1 - 
C     1 m1*p1z*s1 - m1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1 m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1))
C     1 -D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD5A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(-(c3a1*l1*p1y) + 
C     1    c2a1*m1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1x*s1 - 
C     1    p1z*s1 - m1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*l2*p2x + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1    c2a2*l2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1    2*c3a2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2x*s2 + 
C     1    l2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(-(c3a1*m1*p1y) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1    c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1    2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1y*s1 - 
C     1    m1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*m2*p2x - c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2y*s2 - 
C     1    p2z*s2 + l2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(-(c3a1*n1*p1y) + c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + p1x*s1 - 
C     1    m1*p1z*s1 - m1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*n2*p2x - c2a2*l2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 + 
C     1    l2*p2z*s2 + l2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)))/2.
C     1 -D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD5B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(-(c3a1*l1*p1y) + 
C     1    c2a1*m1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1x*s1 - 
C     1    p1z*s1 - m1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*l2*p2y - c2a2*m2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2x*s2 + 
C     1    p2z*s2 + m2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(-(c3a1*m1*p1y) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1    c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1    2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1y*s1 - 
C     1    m1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*m2*p2y + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1    c2a2*m2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1    2*c3a2*m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2y*s2 + 
C     1    m2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(-(c3a1*n1*p1y) + c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + p1x*s1 - 
C     1    m1*p1z*s1 - m1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*n2*p2y - c2a2*m2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - p2x*s2 + 
C     1    m2*p2z*s2 + m2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)))/2.
C     1 -D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD5C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(-(c3a1*l1*p1y) + 
C     1    c2a1*m1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1x*s1 - 
C     1    p1z*s1 - m1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*l2*p2z - c2a2*n2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2x*s2 - 
C     1    p2y*s2 + n2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(-(c3a1*m1*p1y) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1    c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1    2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1y*s1 - 
C     1    m1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*m2*p2z - c2a2*n2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 + 
C     1    n2*p2y*s2 + n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(-(c3a1*n1*p1y) + c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + p1x*s1 - 
C     1    m1*p1z*s1 - m1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*n2*p2z + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1    c2a2*n2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1    2*c3a2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2z*s2 + 
C     1    n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1    n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)))/2.
C     1 -D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD66(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(-(c3a1*l1*p1z) + 
C     1     c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1     2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1x*s1 + 
C     1     p1y*s1 - n1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1     l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)**2 + 
C     1 2*(-(c3a1*m1*p1z) + c2a1*n1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1     2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 - 
C     1     n1*p1y*s1 - n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1     m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)**2 + 
C     1 2*(-(c3a1*n1*p1z) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1     c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1     2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1z*s1 - 
C     1     n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1     n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)**2 + 
C     1 2*((-3*c2a1*n12*(n1*p1y - m1*p1z))/alpha1**4 + 
C     1    (c2a1*l1*n12*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    (8*c3a1*l1*n12*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    c2a1*n12*p1x*ralpha12 + 2*c2a1*n1*p1y*ralpha12 + 
C     1    4*c3a1*l1*n1*p1z*ralpha12 + c2a1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1    2*c3a1*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 + 
C     1    (3*n12*(n1*p1y - m1*p1z)*s1)/alpha1**4 - 
C     1    (5*l1*n12*(l1*p1x + m1*p1y + n1*p1z)*s1)/alpha1**4 + 
C     1    n12*p1x*ralpha12*s1 - 2*n1*p1y*ralpha12*s1 + 
C     1    2*l1*n1*p1z*ralpha12*s1 - (n1*p1y - m1*p1z)*ralpha12*s1 - 
C     1    n12*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
C     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
C     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
C     1 2*((-3*c2a1*n12*(-(n1*p1x) + l1*p1z))/alpha1**4 + 
C     1    (c2a1*m1*n12*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    (8*c3a1*m1*n12*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    2*c2a1*n1*p1x*ralpha12 - c2a1*n12*p1y*ralpha12 + 
C     1    4*c3a1*m1*n1*p1z*ralpha12 + c2a1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1    2*c3a1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 + 
C     1    (3*n12*(-(n1*p1x) + l1*p1z)*s1)/alpha1**4 - 
C     1    (5*m1*n12*(l1*p1x + m1*p1y + n1*p1z)*s1)/alpha1**4 + 
C     1    2*n1*p1x*ralpha12*s1 + n12*p1y*ralpha12*s1 + 
C     1    2*m1*n1*p1z*ralpha12*s1 - (-(n1*p1x) + l1*p1z)*ralpha12*s1 - 
C     1    n12*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
C     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
C     1 2*((-3*c2a1*n12*(m1*p1x - l1*p1y))/alpha1**4 - 2*c3a1*p1z + 
C     1    (c2a1*n1**3*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 - 
C     1    (8*c3a1*n1**3*(l1*p1x + m1*p1y + n1*p1z))/alpha1**4 + 
C     1    c2a1*(m1*p1x - l1*p1y)*ralpha12 - c2a1*n12*p1z*ralpha12 + 
C     1    4*c3a1*n12*p1z*ralpha12 + 
C     1    6*c3a1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
C     1    (3*n12*(m1*p1x - l1*p1y)*s1)/alpha1**4 - p1z*s1 - 
C     1    (5*n1**3*(l1*p1x + m1*p1y + n1*p1z)*s1)/alpha1**4 - 
C     1    (m1*p1x - l1*p1y)*ralpha12*s1 - 
C     1    n12*(m1*p1x - l1*p1y)*ralpha12*s1 + 3*n12*p1z*ralpha12*s1 + 
C     1    3*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
C     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
C     1 -D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD67(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = -(rdist*(-(c3a1*l1*p1z) + c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1 2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1x*s1 + 
C     1 p1y*s1 - n1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1 l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1))
C     1 -D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD68(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = -(rdist*(-(c3a1*m1*p1z) + c2a1*n1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1 2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 - 
C     1 n1*p1y*s1 - n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1 m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1))
C     1 -D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD69(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = -(rdist*(-(c3a1*n1*p1z) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1 c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1 2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1z*s1 - 
C     1 n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1 n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1))
C     1 -D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD6A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(-(c3a1*l1*p1z) + 
C     1    c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1x*s1 + 
C     1    p1y*s1 - n1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*l2*p2x + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1    c2a2*l2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1    2*c3a2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2x*s2 + 
C     1    l2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(-(c3a1*m1*p1z) + c2a1*n1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 - 
C     1    n1*p1y*s1 - n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*m2*p2x - c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2y*s2 - 
C     1    p2z*s2 + l2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(-(c3a1*n1*p1z) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1    c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1    2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1z*s1 - 
C     1    n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1    n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*n2*p2x - c2a2*l2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 + 
C     1    l2*p2z*s2 + l2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)))/2.
C     1 -D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD6B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(-(c3a1*l1*p1z) + 
C     1    c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1x*s1 + 
C     1    p1y*s1 - n1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*l2*p2y - c2a2*m2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2x*s2 + 
C     1    p2z*s2 + m2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(-(c3a1*m1*p1z) + c2a1*n1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 - 
C     1    n1*p1y*s1 - n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*m2*p2y + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1    c2a2*m2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1    2*c3a2*m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2y*s2 + 
C     1    m2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(-(c3a1*n1*p1z) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1    c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1    2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1z*s1 - 
C     1    n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1    n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*n2*p2y - c2a2*m2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - p2x*s2 + 
C     1    m2*p2z*s2 + m2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)))/2.
C     1 -D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD6C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(-(c3a1*l1*p1z) + 
C     1    c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
C     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1x*s1 + 
C     1    p1y*s1 - n1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
C     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*l2*p2z - c2a2*n2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2x*s2 - 
C     1    p2y*s2 + n2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(-(c3a1*m1*p1z) + c2a1*n1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
C     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 - 
C     1    n1*p1y*s1 - n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
C     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*m2*p2z - c2a2*n2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 + 
C     1    n2*p2y*s2 + n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(-(c3a1*n1*p1z) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
C     1    c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
C     1    2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1z*s1 - 
C     1    n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
C     1    n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
C     1  (c3a2*n2*p2z + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1    c2a2*n2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1    2*c3a2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2z*s2 + 
C     1    n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1    n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)))/2.
C     1 -D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD77(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = rdist
C     1 -D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD78(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = 
C     1 -D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD79(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) =
C     1 -D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD7A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = -(rdist*(c3a2*l2*p2x + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1 c2a2*l2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1 2*c3a2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2x*s2 + 
C     1 l2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1 l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2))
C     1 -D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD7B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = -(rdist*(c3a2*l2*p2y - c2a2*m2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1 2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2x*s2 + 
C     1 p2z*s2 + m2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1 l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2))
C     1 -D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD7C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = -(rdist*(c3a2*l2*p2z - c2a2*n2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1 2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2x*s2 - 
C     1 p2y*s2 + n2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1 l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2))
C     1 -D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD88(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = rdist
C     1 -D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD89(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) =
C     1 -D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD8A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = -(rdist*(c3a2*m2*p2x - c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1 2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2y*s2 - 
C     1 p2z*s2 + l2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1 l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2))
C     1 -D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD8B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = -(rdist*(c3a2*m2*p2y + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1 c2a2*m2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1 2*c3a2*m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2y*s2 + 
C     1 m2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1 m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2))
C     1 -D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD8C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = -(rdist*(c3a2*m2*p2z - c2a2*n2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1 2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 + 
C     1 n2*p2y*s2 + n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1 m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2))
C     1 -D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD99(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) = rdist
C     1 -D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD9A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = -(rdist*(c3a2*n2*p2x - c2a2*l2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1 2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 + 
C     1 l2*p2z*s2 + l2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1 l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2))
C     1 -D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD9B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = -(rdist*(c3a2*n2*p2y - c2a2*m2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1 2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - p2x*s2 + 
C     1 m2*p2z*s2 + m2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1 m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2))
C     1 -D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DD9C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = -(rdist*(c3a2*n2*p2z + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1 c2a2*n2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1 2*c3a2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2z*s2 + 
C     1 n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1 n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2))
C     1 -D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DDAA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(c3a2*l2*p2x + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1     c2a2*l2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1     2*c3a2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2x*s2 + 
C     1     l2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1     l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)**2 + 
C     1 2*(c3a2*m2*p2x - c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1     2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2y*s2 - 
C     1     p2z*s2 + l2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1     l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)**2 + 
C     1 2*(c3a2*n2*p2x - c2a2*l2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1     2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 + 
C     1     l2*p2z*s2 + l2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1     l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)**2 + 
C     1 2*(2*c3a2*p2x + (3*c2a2*l22*(n2*p2y - m2*p2z))/alpha2**4 - 
C     1    (c2a2*l2**3*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    (8*c3a2*l2**3*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    c2a2*l22*p2x*ralpha22 - 4*c3a2*l22*p2x*ralpha22 - 
C     1    c2a2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1    6*c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 - 
C     1    (3*l22*(n2*p2y - m2*p2z)*s2)/alpha2**4 + 
C     1    (5*l2**3*(l2*p2x + m2*p2y + n2*p2z)*s2)/alpha2**4 - 
C     1    3*l22*p2x*ralpha22*s2 + (n2*p2y - m2*p2z)*ralpha22*s2 + 
C     1    l22*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    3*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
C     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
C     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
C     1 2*((3*c2a2*l22*(-(n2*p2x) + l2*p2z))/alpha2**4 - 
C     1    (c2a2*l22*m2*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    (8*c3a2*l22*m2*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 - 
C     1    4*c3a2*l2*m2*p2x*ralpha22 + c2a2*l22*p2y*ralpha22 - 
C     1    2*c2a2*l2*p2z*ralpha22 - c2a2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1    2*c3a2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 - 
C     1    (3*l22*(-(n2*p2x) + l2*p2z)*s2)/alpha2**4 + 
C     1    (5*l22*m2*(l2*p2x + m2*p2y + n2*p2z)*s2)/alpha2**4 - 
C     1    2*l2*m2*p2x*ralpha22*s2 - l22*p2y*ralpha22*s2 + 
C     1    2*l2*p2z*ralpha22*s2 + (-(n2*p2x) + l2*p2z)*ralpha22*s2 + 
C     1    l22*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
C     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
C     1 2*((3*c2a2*l22*(m2*p2x - l2*p2y))/alpha2**4 - 
C     1    (c2a2*l22*n2*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    (8*c3a2*l22*n2*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 - 
C     1    4*c3a2*l2*n2*p2x*ralpha22 + 2*c2a2*l2*p2y*ralpha22 - 
C     1    c2a2*(m2*p2x - l2*p2y)*ralpha22 + c2a2*l22*p2z*ralpha22 - 
C     1    2*c3a2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - 
C     1    (3*l22*(m2*p2x - l2*p2y)*s2)/alpha2**4 + p2z*s2 + 
C     1    (5*l22*n2*(l2*p2x + m2*p2y + n2*p2z)*s2)/alpha2**4 - 
C     1    2*l2*n2*p2x*ralpha22*s2 - 2*l2*p2y*ralpha22*s2 + 
C     1    (m2*p2x - l2*p2y)*ralpha22*s2 + 
C     1    l22*(m2*p2x - l2*p2y)*ralpha22*s2 - l22*p2z*ralpha22*s2 - 
C     1    n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
C     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
C     1 -DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DDAB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(c3a2*l2*p2x + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1    c2a2*l2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1    2*c3a2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2x*s2 + 
C     1    l2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c3a2*l2*p2y - c2a2*m2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2x*s2 + 
C     1    p2z*s2 + m2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(c3a2*m2*p2x - c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2y*s2 - 
C     1    p2z*s2 + l2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c3a2*m2*p2y + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1    c2a2*m2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1    2*c3a2*m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2y*s2 + 
C     1    m2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(c3a2*n2*p2x - c2a2*l2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 + 
C     1    l2*p2z*s2 + l2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c3a2*n2*p2y - c2a2*m2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - p2x*s2 + 
C     1    m2*p2z*s2 + m2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(c3a2*p2y + (3*c2a2*l2*m2*(n2*p2y - m2*p2z))/alpha2**4 - 
C     1    (c2a2*l22*m2*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    (8*c3a2*l22*m2*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    c2a2*l2*m2*p2x*ralpha22 - 2*c3a2*l2*m2*p2x*ralpha22 - 
C     1    2*c3a2*l22*p2y*ralpha22 + c2a2*l2*p2z*ralpha22 - 
C     1    2*c3a2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - 
C     1    (3*l2*m2*(n2*p2y - m2*p2z)*s2)/alpha2**4 + 
C     1    (5*l22*m2*(l2*p2x + m2*p2y + n2*p2z)*s2)/alpha2**4 - 
C     1    2*l2*m2*p2x*ralpha22*s2 - l22*p2y*ralpha22*s2 - 
C     1    l2*p2z*ralpha22*s2 + l2*m2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
C     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
C     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
C     1 2*(c3a2*p2x + (3*c2a2*l2*m2*(-(n2*p2x) + l2*p2z))/alpha2**4 - 
C     1    (c2a2*l2*m22*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    (8*c3a2*l2*m22*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 - 
C     1    2*c3a2*m22*p2x*ralpha22 + c2a2*l2*m2*p2y*ralpha22 - 
C     1    2*c3a2*l2*m2*p2y*ralpha22 - c2a2*m2*p2z*ralpha22 - 
C     1    2*c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - 
C     1    (3*l2*m2*(-(n2*p2x) + l2*p2z)*s2)/alpha2**4 + 
C     1    (5*l2*m22*(l2*p2x + m2*p2y + n2*p2z)*s2)/alpha2**4 - 
C     1    m22*p2x*ralpha22*s2 - 2*l2*m2*p2y*ralpha22*s2 + 
C     1    m2*p2z*ralpha22*s2 + l2*m2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
C     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
C     1 2*((3*c2a2*l2*m2*(m2*p2x - l2*p2y))/alpha2**4 - 
C     1    (c2a2*l2*m2*n2*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    (8*c3a2*l2*m2*n2*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 - 
C     1    c2a2*l2*p2x*ralpha22 - 2*c3a2*m2*n2*p2x*ralpha22 + 
C     1    c2a2*m2*p2y*ralpha22 - 2*c3a2*l2*n2*p2y*ralpha22 + 
C     1    c2a2*l2*m2*p2z*ralpha22 - 
C     1    (3*l2*m2*(m2*p2x - l2*p2y)*s2)/alpha2**4 + 
C     1    (5*l2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*s2)/alpha2**4 + 
C     1    l2*p2x*ralpha22*s2 - m2*n2*p2x*ralpha22*s2 - m2*p2y*ralpha22*s2 - 
C     1    l2*n2*p2y*ralpha22*s2 + l2*m2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1    l2*m2*p2z*ralpha22*s2)*
C     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
C     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
C     1 -DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DDAC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(c3a2*l2*p2x + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1    c2a2*l2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1    2*c3a2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2x*s2 + 
C     1    l2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c3a2*l2*p2z - c2a2*n2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2x*s2 - 
C     1    p2y*s2 + n2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(c3a2*m2*p2x - c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2y*s2 - 
C     1    p2z*s2 + l2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c3a2*m2*p2z - c2a2*n2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 + 
C     1    n2*p2y*s2 + n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(c3a2*n2*p2x - c2a2*l2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 + 
C     1    l2*p2z*s2 + l2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c3a2*n2*p2z + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1    c2a2*n2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1    2*c3a2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2z*s2 + 
C     1    n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1    n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(c3a2*p2z + (3*c2a2*l2*n2*(n2*p2y - m2*p2z))/alpha2**4 - 
C     1    (c2a2*l22*n2*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    (8*c3a2*l22*n2*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    c2a2*l2*n2*p2x*ralpha22 - 2*c3a2*l2*n2*p2x*ralpha22 - 
C     1    c2a2*l2*p2y*ralpha22 - 2*c3a2*l22*p2z*ralpha22 - 
C     1    2*c3a2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - 
C     1    (3*l2*n2*(n2*p2y - m2*p2z)*s2)/alpha2**4 + 
C     1    (5*l22*n2*(l2*p2x + m2*p2y + n2*p2z)*s2)/alpha2**4 - 
C     1    2*l2*n2*p2x*ralpha22*s2 + l2*p2y*ralpha22*s2 - 
C     1    l22*p2z*ralpha22*s2 + l2*n2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
C     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
C     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
C     1 2*((3*c2a2*l2*n2*(-(n2*p2x) + l2*p2z))/alpha2**4 - 
C     1    (c2a2*l2*m2*n2*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    (8*c3a2*l2*m2*n2*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    c2a2*l2*p2x*ralpha22 - 2*c3a2*m2*n2*p2x*ralpha22 + 
C     1    c2a2*l2*n2*p2y*ralpha22 - 2*c3a2*l2*m2*p2z*ralpha22 - 
C     1    c2a2*n2*p2z*ralpha22 - 
C     1    (3*l2*n2*(-(n2*p2x) + l2*p2z)*s2)/alpha2**4 + 
C     1    (5*l2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*s2)/alpha2**4 - 
C     1    l2*p2x*ralpha22*s2 - m2*n2*p2x*ralpha22*s2 - 
C     1    l2*n2*p2y*ralpha22*s2 - l2*m2*p2z*ralpha22*s2 + 
C     1    n2*p2z*ralpha22*s2 + l2*n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2)*
C     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
C     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
C     1 2*(c3a2*p2x + (3*c2a2*l2*n2*(m2*p2x - l2*p2y))/alpha2**4 - 
C     1    (c2a2*l2*n22*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    (8*c3a2*l2*n22*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 - 
C     1    2*c3a2*n22*p2x*ralpha22 + c2a2*n2*p2y*ralpha22 + 
C     1    c2a2*l2*n2*p2z*ralpha22 - 2*c3a2*l2*n2*p2z*ralpha22 - 
C     1    2*c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - 
C     1    (3*l2*n2*(m2*p2x - l2*p2y)*s2)/alpha2**4 + 
C     1    (5*l2*n22*(l2*p2x + m2*p2y + n2*p2z)*s2)/alpha2**4 - 
C     1    n22*p2x*ralpha22*s2 - n2*p2y*ralpha22*s2 + 
C     1    l2*n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 2*l2*n2*p2z*ralpha22*s2 - 
C     1    l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
C     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
C     1 -DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DDBB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(c3a2*l2*p2y - c2a2*m2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1     2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2x*s2 + 
C     1     p2z*s2 + m2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1     l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)**2 + 
C     1 2*(c3a2*m2*p2y + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1     c2a2*m2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1     2*c3a2*m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2y*s2 + 
C     1     m2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1     m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)**2 + 
C     1 2*(c3a2*n2*p2y - c2a2*m2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1     2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - p2x*s2 + 
C     1     m2*p2z*s2 + m2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1     m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)**2 + 
C     1 2*((3*c2a2*m22*(n2*p2y - m2*p2z))/alpha2**4 - 
C     1    (c2a2*l2*m22*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    (8*c3a2*l2*m22*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    c2a2*m22*p2x*ralpha22 - 4*c3a2*l2*m2*p2y*ralpha22 + 
C     1    2*c2a2*m2*p2z*ralpha22 - c2a2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1    2*c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 - 
C     1    (3*m22*(n2*p2y - m2*p2z)*s2)/alpha2**4 + 
C     1    (5*l2*m22*(l2*p2x + m2*p2y + n2*p2z)*s2)/alpha2**4 - 
C     1    m22*p2x*ralpha22*s2 - 2*l2*m2*p2y*ralpha22*s2 - 
C     1    2*m2*p2z*ralpha22*s2 + (n2*p2y - m2*p2z)*ralpha22*s2 + 
C     1    m22*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
C     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
C     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
C     1 2*(2*c3a2*p2y + (3*c2a2*m22*(-(n2*p2x) + l2*p2z))/alpha2**4 - 
C     1    (c2a2*m2**3*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    (8*c3a2*m2**3*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    c2a2*m22*p2y*ralpha22 - 4*c3a2*m22*p2y*ralpha22 - 
C     1    c2a2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1    6*c3a2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 - 
C     1    (3*m22*(-(n2*p2x) + l2*p2z)*s2)/alpha2**4 + 
C     1    (5*m2**3*(l2*p2x + m2*p2y + n2*p2z)*s2)/alpha2**4 - 
C     1    3*m22*p2y*ralpha22*s2 + (-(n2*p2x) + l2*p2z)*ralpha22*s2 + 
C     1    m22*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    3*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
C     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
C     1 2*((3*c2a2*m22*(m2*p2x - l2*p2y))/alpha2**4 - 
C     1    (c2a2*m22*n2*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    (8*c3a2*m22*n2*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 - 
C     1    2*c2a2*m2*p2x*ralpha22 - 4*c3a2*m2*n2*p2y*ralpha22 - 
C     1    c2a2*(m2*p2x - l2*p2y)*ralpha22 + c2a2*m22*p2z*ralpha22 - 
C     1    2*c3a2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - 
C     1    (3*m22*(m2*p2x - l2*p2y)*s2)/alpha2**4 + p2z*s2 + 
C     1    (5*m22*n2*(l2*p2x + m2*p2y + n2*p2z)*s2)/alpha2**4 + 
C     1    2*m2*p2x*ralpha22*s2 - 2*m2*n2*p2y*ralpha22*s2 + 
C     1    (m2*p2x - l2*p2y)*ralpha22*s2 + 
C     1    m22*(m2*p2x - l2*p2y)*ralpha22*s2 - m22*p2z*ralpha22*s2 - 
C     1    n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
C     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
C     1 -DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DDBC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(c3a2*l2*p2y - c2a2*m2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2x*s2 + 
C     1    p2z*s2 + m2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c3a2*l2*p2z - c2a2*n2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2x*s2 - 
C     1    p2y*s2 + n2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(c3a2*m2*p2y + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1    c2a2*m2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1    2*c3a2*m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2y*s2 + 
C     1    m2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c3a2*m2*p2z - c2a2*n2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 + 
C     1    n2*p2y*s2 + n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*(c3a2*n2*p2y - c2a2*m2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - p2x*s2 + 
C     1    m2*p2z*s2 + m2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c3a2*n2*p2z + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1    c2a2*n2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1    2*c3a2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2z*s2 + 
C     1    n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1    n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
C     1 2*((3*c2a2*m2*n2*(n2*p2y - m2*p2z))/alpha2**4 - 
C     1    (c2a2*l2*m2*n2*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    (8*c3a2*l2*m2*n2*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    c2a2*m2*n2*p2x*ralpha22 - c2a2*m2*p2y*ralpha22 - 
C     1    2*c3a2*l2*n2*p2y*ralpha22 - 2*c3a2*l2*m2*p2z*ralpha22 + 
C     1    c2a2*n2*p2z*ralpha22 - 
C     1    (3*m2*n2*(n2*p2y - m2*p2z)*s2)/alpha2**4 + 
C     1    (5*l2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*s2)/alpha2**4 - 
C     1    m2*n2*p2x*ralpha22*s2 + m2*p2y*ralpha22*s2 - 
C     1    l2*n2*p2y*ralpha22*s2 - l2*m2*p2z*ralpha22*s2 - 
C     1    n2*p2z*ralpha22*s2 + m2*n2*(n2*p2y - m2*p2z)*ralpha22*s2)*
C     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
C     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
C     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
C     1 2*(c3a2*p2z + (3*c2a2*m2*n2*(-(n2*p2x) + l2*p2z))/alpha2**4 - 
C     1    (c2a2*m22*n2*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    (8*c3a2*m22*n2*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    c2a2*m2*p2x*ralpha22 + c2a2*m2*n2*p2y*ralpha22 - 
C     1    2*c3a2*m2*n2*p2y*ralpha22 - 2*c3a2*m22*p2z*ralpha22 - 
C     1    2*c3a2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - 
C     1    (3*m2*n2*(-(n2*p2x) + l2*p2z)*s2)/alpha2**4 + 
C     1    (5*m22*n2*(l2*p2x + m2*p2y + n2*p2z)*s2)/alpha2**4 - 
C     1    m2*p2x*ralpha22*s2 - 2*m2*n2*p2y*ralpha22*s2 - 
C     1    m22*p2z*ralpha22*s2 + m2*n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
C     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
C     1 2*(c3a2*p2y + (3*c2a2*m2*n2*(m2*p2x - l2*p2y))/alpha2**4 - 
C     1    (c2a2*m2*n22*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    (8*c3a2*m2*n22*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 - 
C     1    c2a2*n2*p2x*ralpha22 - 2*c3a2*n22*p2y*ralpha22 + 
C     1    c2a2*m2*n2*p2z*ralpha22 - 2*c3a2*m2*n2*p2z*ralpha22 - 
C     1    2*c3a2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - 
C     1    (3*m2*n2*(m2*p2x - l2*p2y)*s2)/alpha2**4 + 
C     1    (5*m2*n22*(l2*p2x + m2*p2y + n2*p2z)*s2)/alpha2**4 + 
C     1    n2*p2x*ralpha22*s2 - n22*p2y*ralpha22*s2 + 
C     1    m2*n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 2*m2*n2*p2z*ralpha22*s2 - 
C     1    m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
C     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
C     1 -DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1 DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist

C      DDCC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist) 
C     1 = (rdist*(2*(c3a2*l2*p2z - c2a2*n2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1     2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2x*s2 - 
C     1     p2y*s2 + n2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1     l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)**2 + 
C     1 2*(c3a2*m2*p2z - c2a2*n2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1     2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 + 
C     1     n2*p2y*s2 + n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1     m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)**2 + 
C     1 2*(c3a2*n2*p2z + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
C     1     c2a2*n2*(m2*p2x - l2*p2y)*ralpha22 - 
C     1     2*c3a2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2z*s2 + 
C     1     n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
C     1     n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)**2 + 
C     1 2*((3*c2a2*n22*(n2*p2y - m2*p2z))/alpha2**4 - 
C     1    (c2a2*l2*n22*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    (8*c3a2*l2*n22*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    c2a2*n22*p2x*ralpha22 - 2*c2a2*n2*p2y*ralpha22 - 
C     1    4*c3a2*l2*n2*p2z*ralpha22 - c2a2*(n2*p2y - m2*p2z)*ralpha22 - 
C     1    2*c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 - 
C     1    (3*n22*(n2*p2y - m2*p2z)*s2)/alpha2**4 + 
C     1    (5*l2*n22*(l2*p2x + m2*p2y + n2*p2z)*s2)/alpha2**4 - 
C     1    n22*p2x*ralpha22*s2 + 2*n2*p2y*ralpha22*s2 - 
C     1    2*l2*n2*p2z*ralpha22*s2 + (n2*p2y - m2*p2z)*ralpha22*s2 + 
C     1    n22*(n2*p2y - m2*p2z)*ralpha22*s2 - 
C     1    l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
C     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
C     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
C     1 2*((3*c2a2*n22*(-(n2*p2x) + l2*p2z))/alpha2**4 - 
C     1    (c2a2*m2*n22*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    (8*c3a2*m2*n22*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    2*c2a2*n2*p2x*ralpha22 + c2a2*n22*p2y*ralpha22 - 
C     1    4*c3a2*m2*n2*p2z*ralpha22 - c2a2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
C     1    2*c3a2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 - 
C     1    (3*n22*(-(n2*p2x) + l2*p2z)*s2)/alpha2**4 + 
C     1    (5*m2*n22*(l2*p2x + m2*p2y + n2*p2z)*s2)/alpha2**4 - 
C     1    2*n2*p2x*ralpha22*s2 - n22*p2y*ralpha22*s2 - 
C     1    2*m2*n2*p2z*ralpha22*s2 + (-(n2*p2x) + l2*p2z)*ralpha22*s2 + 
C     1    n22*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
C     1    m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
C     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
C     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
C     1 2*((3*c2a2*n22*(m2*p2x - l2*p2y))/alpha2**4 + 2*c3a2*p2z - 
C     1    (c2a2*n2**3*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 + 
C     1    (8*c3a2*n2**3*(l2*p2x + m2*p2y + n2*p2z))/alpha2**4 - 
C     1    c2a2*(m2*p2x - l2*p2y)*ralpha22 + c2a2*n22*p2z*ralpha22 - 
C     1    4*c3a2*n22*p2z*ralpha22 - 
C     1    6*c3a2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - 
C     1    (3*n22*(m2*p2x - l2*p2y)*s2)/alpha2**4 + p2z*s2 + 
C     1    (5*n2**3*(l2*p2x + m2*p2y + n2*p2z)*s2)/alpha2**4 + 
C     1    (m2*p2x - l2*p2y)*ralpha22*s2 + 
C     1    n22*(m2*p2x - l2*p2y)*ralpha22*s2 - 3*n22*p2z*ralpha22*s2 - 
C     1    3*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
C     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
C     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
C     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
C     1 -DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)*
C     1  DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist)* rdist
C
C  The six reference site positions per capped pentagon. These need to
C  be multiplied by RAD.
C
      SITE(1,1)=1.0D0
      SITE(1,2)=0.0D0
      SITE(1,3)=0.0D0
      SITE(2,1)=((-1.0D0 + Sqrt(5.0D0)))/4.0D0
      SITE(2,2)=(Sqrt((5.0D0 + Sqrt(5.0D0))/2.0D0))/2.0D0
      SITE(2,3)=0.0D0
      SITE(3,1)=((-1.0D0 - Sqrt(5.0D0)))/4.0D0
      SITE(3,2)=(Sqrt((5.0D0 - Sqrt(5.0D0))/2.0D0))/2.0D0
      SITE(3,3)=0.0D0
      SITE(4,1)=((-1 - Sqrt(5.0D0)))/4.0D0
      SITE(4,2)=-(Sqrt((5.0D0 - Sqrt(5.0D0))/2.))/2.0D0
      SITE(4,3)=0.0D0
      SITE(5,1)=((-1 + Sqrt(5.0D0)))/4.0D0
      SITE(5,2)=-(Sqrt((5.0D0 + Sqrt(5.0D0))/2.))/2.0D0
      SITE(5,3)=0.0D0
      SITE(6,1)=0.0D0
      SITE(6,2)=0.0D0
      SITE(6,3)=HEIGHT

      NAT2=NATOMS/2
      EVAP=.FALSE.
      ECAPSID=0.0D0
      DO J1=1,3*NATOMS
         V(J1)=0.0D0
         IF (SECT) THEN
            DO J2=1,3*NATOMS
               HESS(J2,J1)=0.0D0
            ENDDO
         ENDIF
      ENDDO
      
C     EPSEFF=EPS2*(1.0D0+RAD*SQRT((5.0D0+SQRT(5.0D0))/2.0D0))**12
      SIGMA=(1.0D0+RAD*SQRT((5.0D0+SQRT(5.0D0))/2.0D0))
      EPSEFF=EPS2
C     WRITE(*,'(A,4G15.5)') 'RHO,RAD,EPS2,EPSEFF=',RHO,RAD,EPS2,EPSEFF
C
C  Potential energy first.
C
      DO J1=1,NAT2-1
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
         RALPHA12=1.0D0/MAX(ALPHA1**2,1.0D-10)
         CA1=COS(ALPHA1)
         C2A1=CA1
         IF (ALPHA1.LT.0.0001D0) THEN
C           C3A1=-ALPHA1/2+ALPHA1**3/24 ! bug spotted by Tim!
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
            RALPHA22=1.0D0/MAX(ALPHA2**2,1.0D-10)
            CA2=COS(ALPHA2)
            C2A2=CA2
            IF (ALPHA2.LT.0.0001D0) THEN
C              C3A2=-ALPHA2/2+ALPHA2**3/24 ! bug spotted by Tim!
               C3A2=-0.5D0+ALPHA2**2/24.0D0 ! bug spotted by Tim!
               S2=1.0D0-ALPHA2**2/6
            ELSE
               C3A2=(CA2-1.0D0)/ALPHA2**2
               S2=SIN(ALPHA2)/ALPHA2
            ENDIF
C
C  Repulsive site contribution.
C
            P1X=SITE(6,1)*RAD
            P1Y=SITE(6,2)*RAD
            P1Z=SITE(6,3)*RAD
            P2X=SITE(6,1)*RAD
            P2Y=SITE(6,2)*RAD
            P2Z=SITE(6,3)*RAD
            DIST= 
     1  Sqrt((c2a1*p1x-c3a1*l1*(l1*p1x+m1*p1y+n1*p1z)-c2a2*p2x+c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)+ 
     1      (n1*p1y - m1*p1z)*s1 - (n2*p2y - m2*p2z)*s2 + x1 - x2)**2 + 
     1   (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y+c3a2*m2*(l2*p2x+m2*p2y+n2*p2z)+ 
     1      (-(n1*p1x) + l1*p1z)*s1 - (-(n2*p2x) + l2*p2z)*s2 + y1-y2)**2 + 
     1   (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z+c3a2*n2*(l2*p2x+m2*p2y+n2*p2z)+ 
     1      (m1*p1x - l1*p1y)*s1 - (m2*p2x - l2*p2y)*s2 + z1 - z2)**2)
C           PRINT*,'coordinates of pentamer pair:'
C           WRITE(*,'(3F20.10)') X1,Y1,Z1
C           WRITE(*,'(3F20.10)') L1,M1,N1
C           WRITE(*,'(3F20.10)') X2,Y2,Z2
C           WRITE(*,'(3F20.10)') L2,M2,N2
C           PRINT*,'epseff,distance=',EPSEFF,DIST
            DUMMY=EPSEFF*FREP(SIGMA,DIST)
            IF (GTEST.OR.SECT) THEN
               DUMMY2=EPSEFF*DFREP(SIGMA,DIST)
               RDIST=1.0D0/DIST
               D1S=D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
               D2S=D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
               D3S=D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
               D4S=D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
               D5S=D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
               D6S=D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
               D7S=-D1S
               D8S=-D2S
               D9S=-D3S
               DAS=DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
               DBS=DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
               DCS=DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
               GX1=DUMMY2*D1S
               GY1=DUMMY2*D2S
               GZ1=DUMMY2*D3S
               GL1=DUMMY2*D4S
               GM1=DUMMY2*D5S
               GN1=DUMMY2*D6S
               GL2=DUMMY2*DAS
               GM2=DUMMY2*DBS
               GN2=DUMMY2*DCS
            ENDIF
            IF (SECT) THEN
               DUMMY2=EPSEFF*DFREP(SIGMA,DIST)
               DUMMY3=EPSEFF*DDFREP(SIGMA,DIST)
               RDIST=1.0D0/DIST

               S11=DUMMY2*DD11(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*D1S
               S12=DUMMY2*DD12(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*D2S
               S13=DUMMY2*DD13(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*D3S
               S14=DUMMY2*DD14(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*D4S
               S15=DUMMY2*DD15(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*D5S
               S16=DUMMY2*DD16(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*D6S
               S17=DUMMY2*DD17(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*D7S
               S18=DUMMY2*DD18(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*D8S
               S19=DUMMY2*DD19(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*D9S
               S1A=DUMMY2*DD1A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*DAS
               S1B=DUMMY2*DD1B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*DBS
               S1C=DUMMY2*DD1C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*DCS
               S22=DUMMY2*DD22(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*D2S
               S23=DUMMY2*DD23(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*D3S
               S24=DUMMY2*DD24(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*D4S
               S25=DUMMY2*DD25(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*D5S
               S26=DUMMY2*DD26(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*D6S
               S27=DUMMY2*DD27(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*D7S
               S28=DUMMY2*DD28(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*D8S
               S29=DUMMY2*DD29(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*D9S
               S2A=DUMMY2*DD2A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*DAS
               S2B=DUMMY2*DD2B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*DBS
               S2C=DUMMY2*DD2C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*DCS
               S33=DUMMY2*DD33(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*D3S
               S34=DUMMY2*DD34(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*D4S
               S35=DUMMY2*DD35(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*D5S
               S36=DUMMY2*DD36(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*D6S
               S37=DUMMY2*DD37(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*D7S
               S38=DUMMY2*DD38(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*D8S
               S39=DUMMY2*DD39(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*D9S
               S3A=DUMMY2*DD3A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*DAS
               S3B=DUMMY2*DD3B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*DBS
               S3C=DUMMY2*DD3C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*DCS
               S44=DUMMY2*DD44(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D4S*D4S
               S45=DUMMY2*DD45(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D4S*D5S
               S46=DUMMY2*DD46(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D4S*D6S
               S47=DUMMY2*DD47(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D4S*D7S
               S48=DUMMY2*DD48(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D4S*D8S
               S49=DUMMY2*DD49(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D4S*D9S
               S4A=DUMMY2*DD4A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D4S*DAS
               S4B=DUMMY2*DD4B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D4S*DBS
               S4C=DUMMY2*DD4C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D4S*DCS
               S55=DUMMY2*DD55(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D5S*D5S
               S56=DUMMY2*DD56(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D5S*D6S
               S57=DUMMY2*DD57(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D5S*D7S
               S58=DUMMY2*DD58(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D5S*D8S
               S59=DUMMY2*DD59(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D5S*D9S
               S5A=DUMMY2*DD5A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D5S*DAS
               S5B=DUMMY2*DD5B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D5S*DBS
               S5C=DUMMY2*DD5C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D5S*DCS
               S66=DUMMY2*DD66(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D6S*D6S
               S67=DUMMY2*DD67(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D6S*D7S
               S68=DUMMY2*DD68(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D6S*D8S
               S69=DUMMY2*DD69(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D6S*D9S
               S6A=DUMMY2*DD6A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D6S*DAS
               S6B=DUMMY2*DD6B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D6S*DBS
               S6C=DUMMY2*DD6C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D6S*DCS
               S77=DUMMY2*DD77(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D7S*D7S
               S78=DUMMY2*DD78(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D7S*D8S
               S79=DUMMY2*DD79(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D7S*D9S
               S7A=DUMMY2*DD7A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D7S*DAS
               S7B=DUMMY2*DD7B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D7S*DBS
               S7C=DUMMY2*DD7C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D7S*DCS
               S88=DUMMY2*DD88(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D8S*D8S
               S89=DUMMY2*DD89(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D8S*D9S
               S8A=DUMMY2*DD8A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D8S*DAS
               S8B=DUMMY2*DD8B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D8S*DBS
               S8C=DUMMY2*DD8C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D8S*DCS
               S99=DUMMY2*DD99(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D9S*D9S
               S9A=DUMMY2*DD9A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D9S*DAS
               S9B=DUMMY2*DD9B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D9S*DBS
               S9C=DUMMY2*DD9C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D9S*DCS
               SAA=DUMMY2*DDAA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*DAS*DAS
               SAB=DUMMY2*DDAB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*DAS*DBS
               SAC=DUMMY2*DDAC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*DAS*DCS
               SBB=DUMMY2*DDBB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*DBS*DBS
               SBC=DUMMY2*DDBC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*DBS*DCS
               SCC=DUMMY2*DDCC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*DCS*DCS
            ENDIF
C
C  Sum over Morse sites.
C
            DO K1=1,5
               P1X=SITE(K1,1)*RAD
               P1Y=SITE(K1,2)*RAD
               P1Z=SITE(K1,3)*RAD
               DO K2=1,5
                  P2X=SITE(K2,1)*RAD
                  P2Y=SITE(K2,2)*RAD
                  P2Z=SITE(K2,3)*RAD
                  DIST=
     1  Sqrt((c2a1*p1x-c3a1*l1*(l1*p1x+m1*p1y+n1*p1z)-c2a2*p2x+c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)+ 
     1      (n1*p1y - m1*p1z)*s1 - (n2*p2y - m2*p2z)*s2 + x1 - x2)**2 + 
     1   (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y+c3a2*m2*(l2*p2x+m2*p2y+n2*p2z)+ 
     1      (-(n1*p1x) + l1*p1z)*s1 - (-(n2*p2x) + l2*p2z)*s2 + y1-y2)**2 + 
     1   (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z+c3a2*n2*(l2*p2x+m2*p2y+n2*p2z)+ 
     1      (m1*p1x - l1*p1y)*s1 - (m2*p2x - l2*p2y)*s2 + z1 - z2)**2)
                  DUMMY=DUMMY+FATT(RHO,DIST)
                  DUMMY2=DFATT(RHO,DIST)
                  IF (GTEST) THEN
                     RDIST=1.0D0/DIST
                     GX1=GX1+DUMMY2*D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GY1=GY1+DUMMY2*D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GZ1=GZ1+DUMMY2*D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GL1=GL1+DUMMY2*D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GM1=GM1+DUMMY2*D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GN1=GN1+DUMMY2*D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GL2=GL2+DUMMY2*DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GM2=GM2+DUMMY2*DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     GN2=GN2+DUMMY2*DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                  ENDIF
                  IF (SECT) THEN
                     DUMMY2=DFATT(RHO,DIST)
                     DUMMY3=DDFATT(RHO,DIST)
                     RDIST=1.0D0/DIST
                     D1S=D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     D2S=D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     D3S=D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     D4S=D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     D5S=D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     D6S=D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     D7S=D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     D8S=D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     D9S=D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     DAS=DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     DBS=DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
                     DCS=DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)

                     S11=S11+DUMMY2*DD11(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*D1S
                     S12=S12+DUMMY2*DD12(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*D2S
                     S13=S13+DUMMY2*DD13(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*D3S
                     S14=S14+DUMMY2*DD14(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*D4S
                     S15=S15+DUMMY2*DD15(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*D5S
                     S16=S16+DUMMY2*DD16(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*D6S
                     S17=S17+DUMMY2*DD17(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*D7S
                     S18=S18+DUMMY2*DD18(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*D8S
                     S19=S19+DUMMY2*DD19(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*D9S
                     S1A=S1A+DUMMY2*DD1A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*DAS
                     S1B=S1B+DUMMY2*DD1B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*DBS
                     S1C=S1C+DUMMY2*DD1C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D1S*DCS
                     S22=S22+DUMMY2*DD22(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*D2S
                     S23=S23+DUMMY2*DD23(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*D3S
                     S24=S24+DUMMY2*DD24(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*D4S
                     S25=S25+DUMMY2*DD25(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*D5S
                     S26=S26+DUMMY2*DD26(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*D6S
                     S27=S27+DUMMY2*DD27(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*D7S
                     S28=S28+DUMMY2*DD28(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*D8S
                     S29=S29+DUMMY2*DD29(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*D9S
                     S2A=S2A+DUMMY2*DD2A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*DAS
                     S2B=S2B+DUMMY2*DD2B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*DBS
                     S2C=S2C+DUMMY2*DD2C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D2S*DCS
                     S33=S33+DUMMY2*DD33(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*D3S
                     S34=S34+DUMMY2*DD34(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*D4S
                     S35=S35+DUMMY2*DD35(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*D5S
                     S36=S36+DUMMY2*DD36(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*D6S
                     S37=S37+DUMMY2*DD37(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*D7S
                     S38=S38+DUMMY2*DD38(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*D8S
                     S39=S39+DUMMY2*DD39(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*D9S
                     S3A=S3A+DUMMY2*DD3A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*DAS
                     S3B=S3B+DUMMY2*DD3B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*DBS
                     S3C=S3C+DUMMY2*DD3C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D3S*DCS
                     S44=S44+DUMMY2*DD44(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D4S*D4S
                     S45=S45+DUMMY2*DD45(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D4S*D5S
                     S46=S46+DUMMY2*DD46(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D4S*D6S
                     S47=S47+DUMMY2*DD47(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D4S*D7S
                     S48=S48+DUMMY2*DD48(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D4S*D8S
                     S49=S49+DUMMY2*DD49(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D4S*D9S
                     S4A=S4A+DUMMY2*DD4A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D4S*DAS
                     S4B=S4B+DUMMY2*DD4B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D4S*DBS
                     S4C=S4C+DUMMY2*DD4C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D4S*DCS
                     S55=S55+DUMMY2*DD55(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D5S*D5S
                     S56=S56+DUMMY2*DD56(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D5S*D6S
                     S57=S57+DUMMY2*DD57(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D5S*D7S
                     S58=S58+DUMMY2*DD58(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D5S*D8S
                     S59=S59+DUMMY2*DD59(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D5S*D9S
                     S5A=S5A+DUMMY2*DD5A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D5S*DAS
                     S5B=S5B+DUMMY2*DD5B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D5S*DBS
                     S5C=S5C+DUMMY2*DD5C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D5S*DCS
                     S66=S66+DUMMY2*DD66(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D6S*D6S
                     S67=S67+DUMMY2*DD67(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D6S*D7S
                     S68=S68+DUMMY2*DD68(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D6S*D8S
                     S69=S69+DUMMY2*DD69(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D6S*D9S
                     S6A=S6A+DUMMY2*DD6A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D6S*DAS
                     S6B=S6B+DUMMY2*DD6B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D6S*DBS
                     S6C=S6C+DUMMY2*DD6C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D6S*DCS
                     S77=S77+DUMMY2*DD77(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D7S*D7S
                     S78=S78+DUMMY2*DD78(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D7S*D8S
                     S79=S79+DUMMY2*DD79(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D7S*D9S
                     S7A=S7A+DUMMY2*DD7A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D7S*DAS
                     S7B=S7B+DUMMY2*DD7B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D7S*DBS
                     S7C=S7C+DUMMY2*DD7C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D7S*DCS
                     S88=S88+DUMMY2*DD88(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D8S*D8S
                     S89=S89+DUMMY2*DD89(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D8S*D9S
                     S8A=S8A+DUMMY2*DD8A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D8S*DAS
                     S8B=S8B+DUMMY2*DD8B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D8S*DBS
                     S8C=S8C+DUMMY2*DD8C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D8S*DCS
                     S99=S99+DUMMY2*DD99(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D9S*D9S
                     S9A=S9A+DUMMY2*DD9A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D9S*DAS
                     S9B=S9B+DUMMY2*DD9B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D9S*DBS
                     S9C=S9C+DUMMY2*DD9C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*D9S*DCS
                     SAA=SAA+DUMMY2*DDAA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*DAS*DAS
                     SAB=SAB+DUMMY2*DDAB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*DAS*DBS
                     SAC=SAC+DUMMY2*DDAC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*DAS*DCS
                     SBB=SBB+DUMMY2*DDBB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*DBS*DBS
                     SBC=SBC+DUMMY2*DDBC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*DBS*DCS
                     SCC=SCC+DUMMY2*DDCC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)+DUMMY3*DCS*DCS
                  ENDIF 
               ENDDO
            ENDDO
   
            ECAPSID=ECAPSID+DUMMY

            IF (GTEST) THEN
               V(3*(J1-1)+1)=V(3*(J1-1)+1)+GX1
               V(3*(J1-1)+2)=V(3*(J1-1)+2)+GY1
               V(3*(J1-1)+3)=V(3*(J1-1)+3)+GZ1
               V(3*(NAT2+J1-1)+1)=V(3*(NAT2+J1-1)+1)+GL1
               V(3*(NAT2+J1-1)+2)=V(3*(NAT2+J1-1)+2)+GM1
               V(3*(NAT2+J1-1)+3)=V(3*(NAT2+J1-1)+3)+GN1
               V(3*(J2-1)+1)=V(3*(J2-1)+1)-GX1
               V(3*(J2-1)+2)=V(3*(J2-1)+2)-GY1
               V(3*(J2-1)+3)=V(3*(J2-1)+3)-GZ1
               V(3*(NAT2+J2-1)+1)=V(3*(NAT2+J2-1)+1)+GL2
               V(3*(NAT2+J2-1)+2)=V(3*(NAT2+J2-1)+2)+GM2
               V(3*(NAT2+J2-1)+3)=V(3*(NAT2+J2-1)+3)+GN2
            ENDIF
            IF (SECT) THEN
               HESS(3*(J1-1)+1,3*(J1-1)+1)=HESS(3*(J1-1)+1,3*(J1-1)+1)+S11
               HESS(3*(J1-1)+1,3*(J1-1)+2)=HESS(3*(J1-1)+1,3*(J1-1)+2)+S12
               HESS(3*(J1-1)+1,3*(J1-1)+3)=HESS(3*(J1-1)+1,3*(J1-1)+3)+S13
               HESS(3*(J1-1)+1,3*(NAT2+J1-1)+1)=HESS(3*(J1-1)+1,3*(NAT2+J1-1)+1)+S14
               HESS(3*(J1-1)+1,3*(NAT2+J1-1)+2)=HESS(3*(J1-1)+1,3*(NAT2+J1-1)+2)+S15
               HESS(3*(J1-1)+1,3*(NAT2+J1-1)+3)=HESS(3*(J1-1)+1,3*(NAT2+J1-1)+3)+S16
               HESS(3*(J1-1)+1,3*(J2-1)+1)=HESS(3*(J1-1)+1,3*(J2-1)+1)+S17
               HESS(3*(J1-1)+1,3*(J2-1)+2)=HESS(3*(J1-1)+1,3*(J2-1)+2)+S18
               HESS(3*(J1-1)+1,3*(J2-1)+3)=HESS(3*(J1-1)+1,3*(J2-1)+3)+S19
               HESS(3*(J1-1)+1,3*(NAT2+J2-1)+1)=HESS(3*(J1-1)+1,3*(NAT2+J2-1)+1)+S1A
               HESS(3*(J1-1)+1,3*(NAT2+J2-1)+2)=HESS(3*(J1-1)+1,3*(NAT2+J2-1)+2)+S1B
               HESS(3*(J1-1)+1,3*(NAT2+J2-1)+3)=HESS(3*(J1-1)+1,3*(NAT2+J2-1)+3)+S1C

               HESS(3*(J1-1)+2,3*(J1-1)+2)=HESS(3*(J1-1)+2,3*(J1-1)+2)+S22
               HESS(3*(J1-1)+2,3*(J1-1)+3)=HESS(3*(J1-1)+2,3*(J1-1)+3)+S23
               HESS(3*(J1-1)+2,3*(NAT2+J1-1)+1)=HESS(3*(J1-1)+2,3*(NAT2+J1-1)+1)+S24
               HESS(3*(J1-1)+2,3*(NAT2+J1-1)+2)=HESS(3*(J1-1)+2,3*(NAT2+J1-1)+2)+S25
               HESS(3*(J1-1)+2,3*(NAT2+J1-1)+3)=HESS(3*(J1-1)+2,3*(NAT2+J1-1)+3)+S26
               HESS(3*(J1-1)+2,3*(J2-1)+1)=HESS(3*(J1-1)+2,3*(J2-1)+1)+S27
               HESS(3*(J1-1)+2,3*(J2-1)+2)=HESS(3*(J1-1)+2,3*(J2-1)+2)+S28
               HESS(3*(J1-1)+2,3*(J2-1)+3)=HESS(3*(J1-1)+2,3*(J2-1)+3)+S29
               HESS(3*(J1-1)+2,3*(NAT2+J2-1)+1)=HESS(3*(J1-1)+2,3*(NAT2+J2-1)+1)+S2A
               HESS(3*(J1-1)+2,3*(NAT2+J2-1)+2)=HESS(3*(J1-1)+2,3*(NAT2+J2-1)+2)+S2B
               HESS(3*(J1-1)+2,3*(NAT2+J2-1)+3)=HESS(3*(J1-1)+2,3*(NAT2+J2-1)+3)+S2C

               HESS(3*(J1-1)+3,3*(J1-1)+3)=HESS(3*(J1-1)+3,3*(J1-1)+3)+S33
               HESS(3*(J1-1)+3,3*(NAT2+J1-1)+1)=HESS(3*(J1-1)+3,3*(NAT2+J1-1)+1)+S34
               HESS(3*(J1-1)+3,3*(NAT2+J1-1)+2)=HESS(3*(J1-1)+3,3*(NAT2+J1-1)+2)+S35
               HESS(3*(J1-1)+3,3*(NAT2+J1-1)+3)=HESS(3*(J1-1)+3,3*(NAT2+J1-1)+3)+S36
               HESS(3*(J1-1)+3,3*(J2-1)+1)=HESS(3*(J1-1)+3,3*(J2-1)+1)+S37
               HESS(3*(J1-1)+3,3*(J2-1)+2)=HESS(3*(J1-1)+3,3*(J2-1)+2)+S38
               HESS(3*(J1-1)+3,3*(J2-1)+3)=HESS(3*(J1-1)+3,3*(J2-1)+3)+S39
               HESS(3*(J1-1)+3,3*(NAT2+J2-1)+1)=HESS(3*(J1-1)+3,3*(NAT2+J2-1)+1)+S3A
               HESS(3*(J1-1)+3,3*(NAT2+J2-1)+2)=HESS(3*(J1-1)+3,3*(NAT2+J2-1)+2)+S3B
               HESS(3*(J1-1)+3,3*(NAT2+J2-1)+3)=HESS(3*(J1-1)+3,3*(NAT2+J2-1)+3)+S3C

               HESS(3*(NAT2+J1-1)+1,3*(NAT2+J1-1)+1)=HESS(3*(NAT2+J1-1)+1,3*(NAT2+J1-1)+1)+S44
               HESS(3*(NAT2+J1-1)+1,3*(NAT2+J1-1)+2)=HESS(3*(NAT2+J1-1)+1,3*(NAT2+J1-1)+2)+S45
               HESS(3*(NAT2+J1-1)+1,3*(NAT2+J1-1)+3)=HESS(3*(NAT2+J1-1)+1,3*(NAT2+J1-1)+3)+S46
               HESS(3*(NAT2+J1-1)+1,3*(J2-1)+1)=HESS(3*(NAT2+J1-1)+1,3*(J2-1)+1)+S47
               HESS(3*(NAT2+J1-1)+1,3*(J2-1)+2)=HESS(3*(NAT2+J1-1)+1,3*(J2-1)+2)+S48
               HESS(3*(NAT2+J1-1)+1,3*(J2-1)+3)=HESS(3*(NAT2+J1-1)+1,3*(J2-1)+3)+S49
               HESS(3*(NAT2+J1-1)+1,3*(NAT2+J2-1)+1)=HESS(3*(NAT2+J1-1)+1,3*(NAT2+J2-1)+1)+S4A
               HESS(3*(NAT2+J1-1)+1,3*(NAT2+J2-1)+2)=HESS(3*(NAT2+J1-1)+1,3*(NAT2+J2-1)+2)+S4B
               HESS(3*(NAT2+J1-1)+1,3*(NAT2+J2-1)+3)=HESS(3*(NAT2+J1-1)+1,3*(NAT2+J2-1)+3)+S4C

               HESS(3*(NAT2+J1-1)+2,3*(NAT2+J1-1)+2)=HESS(3*(NAT2+J1-1)+2,3*(NAT2+J1-1)+2)+S55
               HESS(3*(NAT2+J1-1)+2,3*(NAT2+J1-1)+3)=HESS(3*(NAT2+J1-1)+2,3*(NAT2+J1-1)+3)+S56
               HESS(3*(NAT2+J1-1)+2,3*(J2-1)+1)=HESS(3*(NAT2+J1-1)+2,3*(J2-1)+1)+S57
               HESS(3*(NAT2+J1-1)+2,3*(J2-1)+2)=HESS(3*(NAT2+J1-1)+2,3*(J2-1)+2)+S58
               HESS(3*(NAT2+J1-1)+2,3*(J2-1)+3)=HESS(3*(NAT2+J1-1)+2,3*(J2-1)+3)+S59
               HESS(3*(NAT2+J1-1)+2,3*(NAT2+J2-1)+1)=HESS(3*(NAT2+J1-1)+2,3*(NAT2+J2-1)+1)+S5A
               HESS(3*(NAT2+J1-1)+2,3*(NAT2+J2-1)+2)=HESS(3*(NAT2+J1-1)+2,3*(NAT2+J2-1)+2)+S5B
               HESS(3*(NAT2+J1-1)+2,3*(NAT2+J2-1)+3)=HESS(3*(NAT2+J1-1)+2,3*(NAT2+J2-1)+3)+S5C

               HESS(3*(NAT2+J1-1)+3,3*(NAT2+J1-1)+3)=HESS(3*(NAT2+J1-1)+3,3*(NAT2+J1-1)+3)+S66
               HESS(3*(NAT2+J1-1)+3,3*(J2-1)+1)=HESS(3*(NAT2+J1-1)+3,3*(J2-1)+1)+S67
               HESS(3*(NAT2+J1-1)+3,3*(J2-1)+2)=HESS(3*(NAT2+J1-1)+3,3*(J2-1)+2)+S68
               HESS(3*(NAT2+J1-1)+3,3*(J2-1)+3)=HESS(3*(NAT2+J1-1)+3,3*(J2-1)+3)+S69
               HESS(3*(NAT2+J1-1)+3,3*(NAT2+J2-1)+1)=HESS(3*(NAT2+J1-1)+3,3*(NAT2+J2-1)+1)+S6A
               HESS(3*(NAT2+J1-1)+3,3*(NAT2+J2-1)+2)=HESS(3*(NAT2+J1-1)+3,3*(NAT2+J2-1)+2)+S6B
               HESS(3*(NAT2+J1-1)+3,3*(NAT2+J2-1)+3)=HESS(3*(NAT2+J1-1)+3,3*(NAT2+J2-1)+3)+S6C

               HESS(3*(J2-1)+1,3*(J2-1)+1)=HESS(3*(J2-1)+1,3*(J2-1)+1)+S77
               HESS(3*(J2-1)+1,3*(J2-1)+2)=HESS(3*(J2-1)+1,3*(J2-1)+2)+S78
               HESS(3*(J2-1)+1,3*(J2-1)+3)=HESS(3*(J2-1)+1,3*(J2-1)+3)+S79
               HESS(3*(J2-1)+1,3*(NAT2+J2-1)+1)=HESS(3*(J2-1)+1,3*(NAT2+J2-1)+1)+S7A
               HESS(3*(J2-1)+1,3*(NAT2+J2-1)+2)=HESS(3*(J2-1)+1,3*(NAT2+J2-1)+2)+S7B
               HESS(3*(J2-1)+1,3*(NAT2+J2-1)+3)=HESS(3*(J2-1)+1,3*(NAT2+J2-1)+3)+S7C

               HESS(3*(J2-1)+2,3*(J2-1)+2)=HESS(3*(J2-1)+2,3*(J2-1)+2)+S88
               HESS(3*(J2-1)+2,3*(J2-1)+3)=HESS(3*(J2-1)+2,3*(J2-1)+3)+S89
               HESS(3*(J2-1)+2,3*(NAT2+J2-1)+1)=HESS(3*(J2-1)+2,3*(NAT2+J2-1)+1)+S8A
               HESS(3*(J2-1)+2,3*(NAT2+J2-1)+2)=HESS(3*(J2-1)+2,3*(NAT2+J2-1)+2)+S8B
               HESS(3*(J2-1)+2,3*(NAT2+J2-1)+3)=HESS(3*(J2-1)+2,3*(NAT2+J2-1)+3)+S8C

               HESS(3*(J2-1)+3,3*(J2-1)+3)=HESS(3*(J2-1)+3,3*(J2-1)+3)+S99
               HESS(3*(J2-1)+3,3*(NAT2+J2-1)+1)=HESS(3*(J2-1)+3,3*(NAT2+J2-1)+1)+S9A
               HESS(3*(J2-1)+3,3*(NAT2+J2-1)+2)=HESS(3*(J2-1)+3,3*(NAT2+J2-1)+2)+S9B
               HESS(3*(J2-1)+3,3*(NAT2+J2-1)+3)=HESS(3*(J2-1)+3,3*(NAT2+J2-1)+3)+S9C

               HESS(3*(NAT2+J2-1)+1,3*(NAT2+J2-1)+1)=HESS(3*(NAT2+J2-1)+1,3*(NAT2+J2-1)+1)+SAA
               HESS(3*(NAT2+J2-1)+1,3*(NAT2+J2-1)+2)=HESS(3*(NAT2+J2-1)+1,3*(NAT2+J2-1)+2)+SAB
               HESS(3*(NAT2+J2-1)+1,3*(NAT2+J2-1)+3)=HESS(3*(NAT2+J2-1)+1,3*(NAT2+J2-1)+3)+SAC

               HESS(3*(NAT2+J2-1)+2,3*(NAT2+J2-1)+2)=HESS(3*(NAT2+J2-1)+2,3*(NAT2+J2-1)+2)+SBB
               HESS(3*(NAT2+J2-1)+2,3*(NAT2+J2-1)+3)=HESS(3*(NAT2+J2-1)+2,3*(NAT2+J2-1)+3)+SBC

               HESS(3*(NAT2+J2-1)+3,3*(NAT2+J2-1)+3)=HESS(3*(NAT2+J2-1)+3,3*(NAT2+J2-1)+3)+SCC
            ENDIF
         ENDDO
      ENDDO

      IF (SECT) THEN
         DO J1=1,NATOMS/2
            DO J2=J1+1,NATOMS/2
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
         DO J1=1,3*NATOMS
            DO J2=J1+1,3*NATOMS
               HESS(J2,J1)=HESS(J1,J2)
            ENDDO
         ENDDO
      ENDIF

C     WRITE(*,'(A,G20.10)') 'energy=',ECAPSID
C     PRINT*,'coords:'
C     WRITE(*,'(I6,G20.10)') (J1,X(J1),J1=1,3*NATOMS)
C     PRINT*,'gradient:'
C     WRITE(*,'(I6,G20.10)') (J1,V(J1),J1=1,3*NATOMS)
C     WRITE(*,'(A,2G20.10)') 'ALPHA1,ALPHA2=',ALPHA1,ALPHA2

      RETURN
      END
