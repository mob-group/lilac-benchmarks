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
C Derivatives of R(i,j) with respect to rigid body coordinates of two
C molecules x1,y1,z1,n1,l1,m1 and x2,y2,z2,l2,m2,n2.
C P(i) = (p1x,p1y,p1z), P(j)=(p2x,p2y,p2z) are the site cordinates in
C the reference geometry for the molecule at the origin and rdist is the inverse
C actual distance between the two sites in question.
C
C This should enable rigid body systems to be coded in a
C straightforward, systematic way.
C

! #######################################################################################################

! Wrapper function for the rigid body first derivatives
! coordNo runs from 1 to 12, with 1 specifying dRij/dx1, 2 specifying dRij/dy1 etc

      function dRijdXrb(coordNo,
     1                  p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)

       implicit none

       DOUBLE PRECISION :: DRIJDXRB
       integer :: coordNo
       DOUBLE PRECISION :: P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,INV_R 
       DOUBLE PRECISION :: C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22

! External functions
       
       DOUBLE PRECISION :: D1,D2,D3,D4,D5,D6,D7,D8,D9,DA,DB,DC

       select case (coordNo)
          case (1)
             dRijdXrb = D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (2)
             dRijdXrb = D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)      
          case (3)
             dRijdXrb = D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (4)
             dRijdXrb = D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (5)
             dRijdXrb = D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (6)
             dRijdXrb = D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (7)
             dRijdXrb = D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (8)
             dRijdXrb = D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (9)
             dRijdXrb = D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (10)
             dRijdXrb = DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (11)
             dRijdXrb = DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (12)
             dRijdXrb = DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case default
             print *, 'Out of range coordinate index passed to dRijdXrb'
             stop

       end select

       return
            
      end function dRijdXrb

! #######################################################################################################

! Wrapper function for the rigid body second derivatives
! coordNo`s run from 1 to 12, with 1 specifying dRij/dx1, 2 specifying dRij/dy1 etc

      function d2RijdXrb1dXrb2(coordNo1,coordNo2,
     1                    p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)

       implicit none

       DOUBLE PRECISION :: D2RIJDXRB1DXRB2
       integer :: coordNo1,coordNo2
       DOUBLE PRECISION :: P1X,P1Y,P1Z,P2X,P2Y,P2Z,X1,Y1,Z1,N1,L1,M1,X2,Y2,Z2,L2,M2,N2,INV_R 
       DOUBLE PRECISION :: C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22

! Local variables

       integer :: firstCoord, secondCoord

! External functions

       DOUBLE PRECISION :: DD11,DD12,DD13,DD14,DD15,DD16,DD17,DD18,DD19,DD1A,DD1B,DD1C
       DOUBLE PRECISION :: DD22,DD23,DD24,DD25,DD26,DD27,DD28,DD29,DD2A,DD2B,DD2C
       DOUBLE PRECISION :: DD33,DD34,DD35,DD36,DD37,DD38,DD39,DD3A,DD3B,DD3C
       DOUBLE PRECISION :: DD44,DD45,DD46,DD47,DD48,DD49,DD4A,DD4B,DD4C      
       DOUBLE PRECISION :: DD55,DD56,DD57,DD58,DD59,DD5A,DD5B,DD5C 
       DOUBLE PRECISION :: DD66,DD67,DD68,DD69,DD6A,DD6B,DD6C
       DOUBLE PRECISION :: DD77,DD78,DD79,DD7A,DD7B,DD7C
       DOUBLE PRECISION :: DD88,DD89,DD8A,DD8B,DD8C
       DOUBLE PRECISION :: DD99,DD9A,DD9B,DD9C
       DOUBLE PRECISION :: DDAA,DDAB,DDAC
       DOUBLE PRECISION :: DDBB,DDBC
       DOUBLE PRECISION :: DDCC

       if (coordNo1.gt.coordNo2) then
          firstCoord = coordNo2
          secondCoord = coordNo1
       else
          firstCoord = coordNo1
          secondCoord = coordNo2          
       end if

       select case (firstCoord)
       case (1)
          select case (secondCoord)
          case (1)   
             d2RijdXrb1dXrb2 = DD11(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (2)
             d2RijdXrb1dXrb2 = DD12(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (3)
             d2RijdXrb1dXrb2 = DD13(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (4)
             d2RijdXrb1dXrb2 = DD14(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (5)
             d2RijdXrb1dXrb2 = DD15(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (6)
             d2RijdXrb1dXrb2 = DD16(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (7)
             d2RijdXrb1dXrb2 = DD17(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (8)
             d2RijdXrb1dXrb2 = DD18(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (9)
             d2RijdXrb1dXrb2 = DD19(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (10)
             d2RijdXrb1dXrb2 = DD1A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (11)
             d2RijdXrb1dXrb2 = DD1B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (12)
             d2RijdXrb1dXrb2 = DD1C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case default
             print *, 'Out of range second coordinate index passed to d2RijdXrb1dXrb2'
             stop
          end select

       case (2)
          select case (secondCoord)
          case (2)
             d2RijdXrb1dXrb2 = DD22(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (3)
             d2RijdXrb1dXrb2 = DD23(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (4)
             d2RijdXrb1dXrb2 = DD24(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (5)
             d2RijdXrb1dXrb2 = DD25(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (6)
             d2RijdXrb1dXrb2 = DD26(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (7)
             d2RijdXrb1dXrb2 = DD27(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (8)
             d2RijdXrb1dXrb2 = DD28(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (9)
             d2RijdXrb1dXrb2 = DD29(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (10)
             d2RijdXrb1dXrb2 = DD2A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (11)
             d2RijdXrb1dXrb2 = DD2B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (12)
             d2RijdXrb1dXrb2 = DD2C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case default
             print *, 'Out of range second coordinate index passed to d2RijdXrb1dXrb2'
             stop
          end select
      
       case (3)
          select case (secondCoord)
          case (3)
             d2RijdXrb1dXrb2 = DD33(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (4)
             d2RijdXrb1dXrb2 = DD34(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (5)
             d2RijdXrb1dXrb2 = DD35(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (6)
             d2RijdXrb1dXrb2 = DD36(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (7)
             d2RijdXrb1dXrb2 = DD37(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (8)
             d2RijdXrb1dXrb2 = DD38(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (9)
             d2RijdXrb1dXrb2 = DD39(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (10)
             d2RijdXrb1dXrb2 = DD3A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (11)
             d2RijdXrb1dXrb2 = DD3B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (12)
             d2RijdXrb1dXrb2 = DD3C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case default
             print *, 'Out of range second coordinate index passed to d2RijdXrb1dXrb2'
             stop
          end select

       case (4)
          select case (secondCoord)
          case (4)
             d2RijdXrb1dXrb2 = DD44(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (5)
             d2RijdXrb1dXrb2 = DD45(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (6)
             d2RijdXrb1dXrb2 = DD46(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (7)
             d2RijdXrb1dXrb2 = DD47(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (8)
             d2RijdXrb1dXrb2 = DD48(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (9)
             d2RijdXrb1dXrb2 = DD49(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (10)
             d2RijdXrb1dXrb2 = DD4A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (11)
             d2RijdXrb1dXrb2 = DD4B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (12)
             d2RijdXrb1dXrb2 = DD4C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case default
             print *, 'Out of range second coordinate index passed to d2RijdXrb1dXrb2'
             stop
          end select

       case (5)
          select case (secondCoord)
          case (5)
             d2RijdXrb1dXrb2 = DD55(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (6)
             d2RijdXrb1dXrb2 = DD56(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (7)
             d2RijdXrb1dXrb2 = DD57(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (8)
             d2RijdXrb1dXrb2 = DD58(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (9)
             d2RijdXrb1dXrb2 = DD59(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (10)
             d2RijdXrb1dXrb2 = DD5A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (11)
             d2RijdXrb1dXrb2 = DD5B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (12)
             d2RijdXrb1dXrb2 = DD5C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case default
             print *, 'Out of range second coordinate index passed to d2RijdXrb1dXrb2'
             stop
          end select

       case (6)
          select case (secondCoord)
          case (6)
             d2RijdXrb1dXrb2 = DD66(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (7)
             d2RijdXrb1dXrb2 = DD67(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (8)
             d2RijdXrb1dXrb2 = DD68(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (9)
             d2RijdXrb1dXrb2 = DD69(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (10)
             d2RijdXrb1dXrb2 = DD6A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (11)
             d2RijdXrb1dXrb2 = DD6B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (12)
             d2RijdXrb1dXrb2 = DD6C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case default
             print *, 'Out of range second coordinate index passed to d2RijdXrb1dXrb2'
             stop
          end select

       case (7)
          select case (secondCoord)
          case (7)
             d2RijdXrb1dXrb2 = DD77(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (8)
             d2RijdXrb1dXrb2 = DD78(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (9)
             d2RijdXrb1dXrb2 = DD79(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (10)
             d2RijdXrb1dXrb2 = DD7A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (11)
             d2RijdXrb1dXrb2 = DD7B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (12)
             d2RijdXrb1dXrb2 = DD7C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case default
             print *, 'Out of range second coordinate index passed to d2RijdXrb1dXrb2'
             stop
          end select

       case (8)
          select case (secondCoord)
          case (8)
             d2RijdXrb1dXrb2 = DD88(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (9)
             d2RijdXrb1dXrb2 = DD89(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (10)
             d2RijdXrb1dXrb2 = DD8A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (11)
             d2RijdXrb1dXrb2 = DD8B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (12)
             d2RijdXrb1dXrb2 = DD8C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case default
             print *, 'Out of range second coordinate index passed to d2RijdXrb1dXrb2'
             stop
          end select

       case (9)
          select case (secondCoord)
          case (9)
             d2RijdXrb1dXrb2 = DD99(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (10)
             d2RijdXrb1dXrb2 = DD9A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (11)
             d2RijdXrb1dXrb2 = DD9B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (12)
             d2RijdXrb1dXrb2 = DD9C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case default
             print *, 'Out of range second coordinate index passed to d2RijdXrb1dXrb2'
             stop
          end select

       case (10)
          select case (secondCoord)
          case (10)
             d2RijdXrb1dXrb2 = DDAA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (11)
             d2RijdXrb1dXrb2 = DDAB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (12)
             d2RijdXrb1dXrb2 = DDAC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case default
             print *, 'Out of range second coordinate index passed to d2RijdXrb1dXrb2'
             stop
          end select

       case (11)
          select case (secondCoord)
          case (11)
             d2RijdXrb1dXrb2 = DDBB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case (12)
             d2RijdXrb1dXrb2 = DDBC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case default
             print *, 'Out of range second coordinate index passed to d2RijdXrb1dXrb2'
             stop
          end select

       case (12)
          select case (secondCoord)
          case (12)
             d2RijdXrb1dXrb2 = DDCC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,inv_R,
     1                              C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
          case default
             print *, 'Out of range second coordinate index passed to d2RijdXrb1dXrb2'
             stop
          end select

       case default
          print *, 'Out of range first coordinate index passed to d2RijdXrb1dXrb2'
          stop
          
       end select
       
       return
            
      end function d2RijdXrb1dXrb2

! #######################################################################################################

! drij/dx1

      DOUBLE PRECISION FUNCTION D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist

      D1=
     @    rdist*(c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     @    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - n2*p2y*s2 + 
     @    m2*p2z*s2 + x1 - x2)
      
      RETURN
      END

! drij/dy1

      DOUBLE PRECISION FUNCTION D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist

      D2 =
     @    rdist*(c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     @    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + n2*p2x*s2 - 
     @    l2*p2z*s2 + y1 - y2)

      RETURN
      END

! drij/dz1

      DOUBLE PRECISION FUNCTION D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist

      D3=
     @    rdist*(c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - m2*p2x*s2 + 
     @    l2*p2y*s2 + z1 - z2)

      RETURN
      END

! drij/dl1

      DOUBLE PRECISION FUNCTION D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist

      D4 =
     @      rdist*((c3a1*(2*l1*p1x*(-1 + l12*ralpha12) + 
     @            (m1*p1y + n1*p1z)*(-1 + 2*l12*ralpha12)) + 
     @         (l1*l12*p1x + l12*m1*p1y - l1*n1*p1y + l1*m1*p1z + l12*n1*p1z)*
     @          ralpha12*s1 + l1*(c2a1*(n1*p1y - m1*p1z)*ralpha12 - p1x*s1))*
     @       (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     @         c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     @         n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + 
     @      (c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     @         c3a1*m1*(-p1x + 2*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @         (-(l1*p1y) + p1z + l1* (l1*m1*p1x + n1*p1x + m1**2*p1y - l1*p1z + m1*n1*p1z)*ralpha12)*
     @          s1)*(c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     @         c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
     @         n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + 
     @      (-(c2a1*l1*(-(m1*p1x) + l1*p1y)*ralpha12) + 
     @         c3a1*n1*(-p1x + 2*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @         (-p1y - l1*p1z + l1*(-(m1*p1x) + l1*n1*p1x + l1*p1y + m1*n1*p1y + n1**2*p1z)*ralpha12)*s1)*
     @       (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @         c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - 
     @         m2*p2x*s2 + l2*p2y*s2 + z1 - z2))

      RETURN
      END

! drij/dm1

      DOUBLE PRECISION FUNCTION D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist

      D5=
     @ rdist*((-(c2a1*m1*(-(n1*p1y) + m1*p1z)*ralpha12) + 
     @         c3a1*l1*(-p1y + 2*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @         (-(m1*p1x) - p1z + m1* (l1**2*p1x + l1*m1*p1y - n1*p1y + m1*p1z + l1*n1*p1z)*ralpha12)*
     @          s1)*(c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     @         c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     @         n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + (-(c3a1*(l1*p1x + 2*m1*p1y + n1*p1z)) + 
     @         c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
     @         (-(m1*p1y) + (l1*m12*p1x + m1*n1*p1x + m1*m12*p1y - l1*m1*p1z + 
     @               m12*n1*p1z)*ralpha12)*s1)*
     @       (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     @         c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
     @         n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + (c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
     @         c3a1*n1*(-p1y + 2*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @         (p1x - m1*(m1 - l1*n1)*p1x*ralpha12 + 
     @            m1*((l1 + m1*n1)*p1y*ralpha12 + p1z*(-1 + n1**2*ralpha12)))*s1)*
     @       (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @         c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - 
     @         m2*p2x*s2 + l2*p2y*s2 + z1 - z2))

      RETURN
      END

! drij/dn1

      DOUBLE PRECISION FUNCTION D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist

      D6=
     @ rdist*((c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
     @         c3a1*l1*(-p1z + 2*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + (-(n1*p1x) + p1y + n1*
     @             (l1**2*p1x + l1*m1*p1y - n1*p1y + m1*p1z + l1*n1*p1z)*ralpha12)*
     @          s1)*(c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     @         c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     @         n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + 
     @      (-(c2a1*n1*(n1*p1x - l1*p1z)*ralpha12) + 
     @         c3a1*m1*(-p1z + 2*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12) + 
     @         (-p1x - n1*p1y + n1*(l1*m1*p1x + n1*p1x + m1**2*p1y - l1*p1z + 
     @               m1*n1*p1z)*ralpha12)*s1)*
     @       (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     @         c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
     @         n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + 
     @      (-(c3a1*(l1*p1x + m1*p1y + 2*n1*p1z)) + 
     @         c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
     @         2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
     @         (-(n1*p1z) + (-(m1*n1*p1x) + l1*n12*p1x + l1*n1*p1y + m1*n12*p1y + 
     @               n1*n12*p1z)*ralpha12)*s1)*
     @       (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @         c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - 
     @         m2*p2x*s2 + l2*p2y*s2 + z1 - z2))

      RETURN
      END

! drij/dx2

      DOUBLE PRECISION FUNCTION D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist

      D7=
     @   -rdist*(c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     @    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - n2*p2y*s2 + 
     @    m2*p2z*s2 + x1 - x2)

      RETURN
      END

! drij/dy2

      DOUBLE PRECISION FUNCTION D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist

      D8=
     @   -rdist*(c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     @    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + n2*p2x*s2 - 
     @    l2*p2z*s2 + y1 - y2)

      RETURN
      END

! drij/dz2

      DOUBLE PRECISION FUNCTION D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist

      D9=
     @   -rdist*(c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - m2*p2x*s2 + 
     @    l2*p2y*s2 + z1 - z2)

      RETURN
      END

! drij/dl2

      DOUBLE PRECISION FUNCTION DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist

      DA=
     @         rdist*((c3a2*(2*l2*p2x + m2*p2y + n2*p2z - 
     @            2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) - 
     @         (l2*l22*p2x + l22*m2*p2y - l2*n2*p2y + l2*m2*p2z + l22*n2*p2z)*
     @          ralpha22*s2 + l2*(-(c2a2*n2*p2y*ralpha22) + c2a2*m2*p2z*ralpha22 + 
     @            p2x*s2))*(c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - 
     @         c2a2*p2x + c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + 
     @         (n1*p1y - m1*p1z)*s1 - n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + 
     @      (-(c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22) + 
     @         c3a2*m2*(p2x - 2*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) - (-(l2*p2y) + p2z + l2*
     @             (l2*m2*p2x + n2*p2x + m2**2*p2y - l2*p2z + m2*n2*p2z)*ralpha22)*
     @          s2)*(c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     @         c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
     @         n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + 
     @      (c2a2*l2*(-(m2*p2x) + l2*p2y)*ralpha22 + 
     @         c3a2*n2*(p2x - 2*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
     @    (p2y + l2*p2z - l2*(-(m2*p2x) + l2*n2*p2x + l2*p2y + m2*n2*p2y + n2**2*p2z)*ralpha22)*s2)*
     @       (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - m2*p2x*s2 + l2*p2y*s2 + z1 - z2))

      RETURN
      END

! drij/dm2

      DOUBLE PRECISION FUNCTION DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist

      DB=
     @ rdist*((c2a2*m2*(-(n2*p2y) + m2*p2z)*ralpha22 + 
     @         c3a2*l2*(p2y - 2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
     @         (m2*p2x + p2z - m2*(l2**2*p2x + l2*m2*p2y - n2*p2y + m2*p2z + l2*n2*p2z)*ralpha22)*s2)*
     @       (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     @         c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     @         n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + (c3a2*(l2*p2x + 2*m2*p2y + n2*p2z - 
     @            2*m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) - 
     @         (l2*m22*p2x + m2*n2*p2x + m2*m22*p2y - l2*m2*p2z + m22*n2*p2z)*
     @          ralpha22*s2 + m2*(c2a2*(n2*p2x - l2*p2z)*ralpha22 + p2y*s2))*
     @       (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     @         c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
     @         n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + (-(c2a2*m2*(m2*p2x - l2*p2y)*ralpha22) + 
     @         c3a2*n2*(p2y - 2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
     @         (-p2x + m2*p2z - m2*(-(m2*p2x) + l2*n2*p2x + l2*p2y + m2*n2*p2y + 
     @               n2**2*p2z)*ralpha22)*s2)*
     @       (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @         c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - 
     @         m2*p2x*s2 + l2*p2y*s2 + z1 - z2))

      RETURN
      END

! drij/dn2

      DOUBLE PRECISION FUNCTION DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist

      DC=
     @ rdist*((-(c2a2*n2*(n2*p2y - m2*p2z)*ralpha22) + 
     @         c3a2*l2*(p2z - 2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
     @         (n2*p2x - p2y - n2*(l2**2*p2x + l2*m2*p2y - n2*p2y + m2*p2z + l2*n2*p2z)*ralpha22)*s2)*
     @       (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     @         c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     @         n2*p2y*s2 + m2*p2z*s2 + x1 - x2) + (c2a2*n2*(n2*p2x - l2*p2z)*ralpha22 + 
     @         c3a2*m2*(p2z - 2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
     @         (p2x + n2*p2y - n2*(l2*m2*p2x + n2*p2x + m2**2*p2y - l2*p2z + m2*n2*p2z)*ralpha22)*s2)*
     @       (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     @         c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) - n1*p1x*s1 + l1*p1z*s1 + 
     @         n2*p2x*s2 - l2*p2z*s2 + y1 - y2) + (c3a2*(l2*p2x + m2*p2y + 2*n2*p2z - 
     @            2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22) + 
     @         (n2*(m2*p2x - l2*p2y) - n22*(l2*p2x + m2*p2y + n2*p2z))*ralpha22*
     @          s2 + n2*(-(c2a2*m2*p2x*ralpha22) + c2a2*l2*p2y*ralpha22 + p2z*s2))*
     @       (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     @         c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + m1*p1x*s1 - l1*p1y*s1 - 
     @         m2*p2x*s2 + l2*p2y*s2 + z1 - z2))

       RETURN
       END

! Second derivatives of R(i,j) with respect to rigid body coordinates.

       DOUBLE PRECISION FUNCTION DD11(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                                C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D1

       DD11 = rdist*(1.0d0 -
     1               D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)**2)

       return
      end FUNCTION DD11

      DOUBLE PRECISION FUNCTION DD12(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D1, D2

       DD12 = -rdist *
     1        D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1        D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)

       return
      end FUNCTION DD12

      DOUBLE PRECISION FUNCTION DD13(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D1, D3

       DD13 = -rdist *
     1        D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1        D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)

       return
      end FUNCTION DD13

      DOUBLE PRECISION FUNCTION DD14(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D1, D4

       DD14 = rdist*(-(c3a1*l1*p1x) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1        c2a1*l1*(n1*p1y - m1*p1z)*ralpha12 + 
     1        2*c3a1*l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1x*s1 - 
     1        l1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1        l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) -
     1        D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1        D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD14

      DOUBLE PRECISION FUNCTION DD15(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist      
       DOUBLE PRECISION D1, D5
       
       DD15 = rdist*(-(c3a1*l1*p1y) + c2a1*m1*(n1*p1y - m1*p1z)*ralpha12 + 
     1        2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1x*s1 - p1z*s1 - 
     1        m1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1        l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) -
     1        D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1        D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD15

      DOUBLE PRECISION FUNCTION DD16(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D1, D6

       DD16 = rdist*(-(c3a1*l1*p1z) + c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
     1        2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1x*s1 + p1y*s1 - 
     1        n1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1        l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) -
     1        D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1        D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD16

      DOUBLE PRECISION FUNCTION DD17(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D1, D7

       DD17 = -rdist*(1.0d0 + 
     1                D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                   C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1                D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                   C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22))
      
       return
      end FUNCTION DD17

      DOUBLE PRECISION FUNCTION DD18(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D1, D8

       DD18 = -D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1            C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1         D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1            C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD18

      DOUBLE PRECISION FUNCTION DD19(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D1, D9

       DD19 = -D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1            C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1         D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1            C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD19

      DOUBLE PRECISION FUNCTION DD1A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D1, DA

       DD1A = rdist*(c3a2*l2*p2x + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1        c2a2*l2*(n2*p2y - m2*p2z)*ralpha22 - 
     1        2*c3a2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2x*s2 + 
     1        l2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1        l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) -
     1        D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1        DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD1A

      DOUBLE PRECISION FUNCTION DD1B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D1, DB

       DD1B = rdist*(c3a2*l2*p2y - c2a2*m2*(n2*p2y - m2*p2z)*ralpha22 - 
     1        2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2x*s2 + p2z*s2 + 
     1        m2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1        l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) -
     1        D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1        DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD1B

      DOUBLE PRECISION FUNCTION DD1C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D1, DC

       DD1C = rdist*(c3a2*l2*p2z - c2a2*n2*(n2*p2y - m2*p2z)*ralpha22 - 
     1        2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2x*s2 - p2y*s2 + 
     1        n2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1        l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) -
     1        D1(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1        DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD1C

      DOUBLE PRECISION FUNCTION DD22(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
      IMPLICIT NONE
      DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
      DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
      DOUBLE PRECISION D2

      DD22 = rdist*(1.0d0 -
     1              D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                 C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)**2)

      return
      end FUNCTION DD22

      DOUBLE PRECISION FUNCTION DD23(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D2, D3

       DD23 = -rdist *
     1        D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1        D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)

       return
      end FUNCTION DD23

      DOUBLE PRECISION FUNCTION DD24(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D2, D4

      DD24 
     1 = rdist*(-(c3a1*m1*p1x) + c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1 2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1y*s1 + p1z*s1 - 
     1 l1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1 l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)
     1 -D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1  D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD24

      DOUBLE PRECISION FUNCTION DD25(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D2, D5

      DD25  
     1 = rdist*(-(c3a1*m1*p1y) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1 c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1 2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1y*s1 - 
     1 m1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1 m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)
     1 -D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1  D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD25

      DOUBLE PRECISION FUNCTION DD26(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D2, D6

      DD26  
     1 = rdist*(-(c3a1*m1*p1z) + c2a1*n1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1 2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 - n1*p1y*s1 - 
     1 n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1 m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)
     1 -D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1  D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD26

      DOUBLE PRECISION FUNCTION DD27(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D2, D7

       DD27 = -D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1            C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1         D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1            C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD27

      DOUBLE PRECISION FUNCTION DD28(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D2, D8

       DD28 = -rdist*(1.0d0 + 
     1                D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                   C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1                D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                   C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22))
      
       return
      end FUNCTION DD28

      DOUBLE PRECISION FUNCTION DD29(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D2, D9

       DD29 = -D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1            C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1         D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1            C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD29

      DOUBLE PRECISION FUNCTION DD2A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D2, DA

      DD2A  
     1 = rdist*(c3a2*m2*p2x - c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1 2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2y*s2 - p2z*s2 + 
     1 l2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1 l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)
     1 -D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1  DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD2A

      DOUBLE PRECISION FUNCTION DD2B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D2, DB

      DD2B  
     1 = rdist*(c3a2*m2*p2y + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1 c2a2*m2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1 2*c3a2*m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2y*s2 + 
     1 m2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1 m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)
     1 -D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1  DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD2B

      DOUBLE PRECISION FUNCTION DD2C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D2, DC

      DD2C  
     1 = rdist*(c3a2*m2*p2z - c2a2*n2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1 2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 + n2*p2y*s2 + 
     1 n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1 m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)
     1 -D2(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1  DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD2C

      DOUBLE PRECISION FUNCTION DD33(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D3

       DD33 = rdist*(1.0d0 -
     1               D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)**2)

       return
      end FUNCTION DD33

      DOUBLE PRECISION FUNCTION DD34(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D3, D4

      DD34  
     1 = rdist*(-(c3a1*n1*p1x) + c2a1*l1*(m1*p1x - l1*p1y)*ralpha12 + 
     1 2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 - l1*p1z*s1 - 
     1 l1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1 l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)
     1 -D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1  D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD34

      DOUBLE PRECISION FUNCTION DD35 (p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D3, D5

      DD35  
     1 = rdist*(-(c3a1*n1*p1y) + c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
     1 2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + p1x*s1 - m1*p1z*s1 - 
     1 m1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1 m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)
     1 -D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1  D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD35

      DOUBLE PRECISION FUNCTION DD36(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D3, D6

      DD36  
     1 = rdist*(-(c3a1*n1*p1z) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1 c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
     1 2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1z*s1 - 
     1 n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1 n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)
     1 -D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1  D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD36

      DOUBLE PRECISION FUNCTION DD37(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D3, D7

       DD37 = -D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1            C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1         D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1            C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD37

      DOUBLE PRECISION FUNCTION DD38(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D3, D8

       DD38 = -D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1            C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1         D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1            C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD38

      DOUBLE PRECISION FUNCTION DD39(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D3, D9

       DD39 = -rdist*(1.0d0 + 
     1                D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                   C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1                D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                   C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22))
      
       return
      end FUNCTION DD39

      DOUBLE PRECISION FUNCTION DD3A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D3, DA

      DD3A  
     1 = rdist*(c3a2*n2*p2x - c2a2*l2*(m2*p2x - l2*p2y)*ralpha22 - 
     1 2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 + l2*p2z*s2 + 
     1 l2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1 l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)
     1 -D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1  DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD3A

      DOUBLE PRECISION FUNCTION DD3B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D3, DB

      DD3B  
     1 = rdist*(c3a2*n2*p2y - c2a2*m2*(m2*p2x - l2*p2y)*ralpha22 - 
     1 2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - p2x*s2 + m2*p2z*s2 + 
     1 m2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1 m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)
     1 -D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1  DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD3B

      DOUBLE PRECISION FUNCTION DD3C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D3, DC

      DD3C  
     1 = rdist*(c3a2*n2*p2z + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1 c2a2*n2*(m2*p2x - l2*p2y)*ralpha22 - 
     1 2*c3a2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2z*s2 + 
     1 n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1 n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)
     1 -D3(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1  DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD3C

      DOUBLE PRECISION FUNCTION DD44(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D4

      DD44  
     1 = (rdist*(2*(-(c3a1*l1*p1x) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1     c2a1*l1*(n1*p1y - m1*p1z)*ralpha12 + 
     1     2*c3a1*l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1x*s1 - 
     1     l1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1     l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)**2 + 
     1 2*(-(c3a1*m1*p1x) + c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1     2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1y*s1 + 
     1     p1z*s1 - l1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1     l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)**2 + 
     1 2*(-(c3a1*n1*p1x) + c2a1*l1*(m1*p1x - l1*p1y)*ralpha12 + 
     1     2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 - 
     1     l1*p1z*s1 - l1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1     l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)**2 + 
     1 2*(-2*c3a1*p1x - (3*c2a1*l12*(n1*p1y - m1*p1z))*(ralpha12**2) + 
     1    (c2a1*l1**3*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    (8*c3a1*l1**3*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    c2a1*l12*p1x*ralpha12 + 4*c3a1*l12*p1x*ralpha12 + 
     1    c2a1*(n1*p1y - m1*p1z)*ralpha12 + 
     1    6*c3a1*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 + 
     1    (3*l12*(n1*p1y - m1*p1z)*s1)*(ralpha12**2) - 
     1    (5*l1**3*(l1*p1x + m1*p1y + n1*p1z)*s1)*(ralpha12**2) + 
     1    3*l12*p1x*ralpha12*s1 - (n1*p1y - m1*p1z)*ralpha12*s1 - 
     1    l12*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    3*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
     1 2*((-3*c2a1*l12*(-(n1*p1x) + l1*p1z))*(ralpha12**2) + 
     1    (c2a1*l12*m1*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    (8*c3a1*l12*m1*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) + 
     1    4*c3a1*l1*m1*p1x*ralpha12 - c2a1*l12*p1y*ralpha12 + 
     1    2*c2a1*l1*p1z*ralpha12 + c2a1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1    2*c3a1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 + 
     1    (3*l12*(-(n1*p1x) + l1*p1z)*s1)*(ralpha12**2) - 
     1    (5*l12*m1*(l1*p1x + m1*p1y + n1*p1z)*s1)*(ralpha12**2) + 
     1    2*l1*m1*p1x*ralpha12*s1 + l12*p1y*ralpha12*s1 - 
     1    2*l1*p1z*ralpha12*s1 - (-(n1*p1x) + l1*p1z)*ralpha12*s1 - 
     1    l12*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
     1 2*((-3*c2a1*l12*(m1*p1x - l1*p1y))*(ralpha12**2) + 
     1    (c2a1*l12*n1*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    (8*c3a1*l12*n1*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) + 
     1    4*c3a1*l1*n1*p1x*ralpha12 - 2*c2a1*l1*p1y*ralpha12 + 
     1    c2a1*(m1*p1x - l1*p1y)*ralpha12 - c2a1*l12*p1z*ralpha12 + 
     1    2*c3a1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
     1    (3*l12*(m1*p1x - l1*p1y)*s1)*(ralpha12**2) - p1z*s1 - 
     1    (5*l12*n1*(l1*p1x + m1*p1y + n1*p1z)*s1)*(ralpha12**2) + 
     1    2*l1*n1*p1x*ralpha12*s1 + 2*l1*p1y*ralpha12*s1 - 
     1    (m1*p1x - l1*p1y)*ralpha12*s1 - 
     1    l12*(m1*p1x - l1*p1y)*ralpha12*s1 + l12*p1z*ralpha12*s1 + 
     1    n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
     1 -D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1  D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD44

      DOUBLE PRECISION FUNCTION DD45(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D4, D5

      DD45  
     1 = (rdist*(2*(-(c3a1*l1*p1x) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1    c2a1*l1*(n1*p1y - m1*p1z)*ralpha12 + 
     1    2*c3a1*l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1x*s1 - 
     1    l1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (-(c3a1*l1*p1y) + c2a1*m1*(n1*p1y - m1*p1z)*ralpha12 + 
     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1x*s1 - 
     1    p1z*s1 - m1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) + 
     1 2*(-(c3a1*m1*p1x) + c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1y*s1 + 
     1    p1z*s1 - l1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (-(c3a1*m1*p1y) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1    c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1    2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1y*s1 - 
     1    m1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) + 
     1 2*(-(c3a1*n1*p1x) + c2a1*l1*(m1*p1x - l1*p1y)*ralpha12 + 
     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 - 
     1    l1*p1z*s1 - l1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (-(c3a1*n1*p1y) + c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + p1x*s1 - 
     1    m1*p1z*s1 - m1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) + 
     1 2*(-(c3a1*p1y) - (3*c2a1*l1*m1*(n1*p1y - m1*p1z))*(ralpha12**2) + 
     1    (c2a1*l12*m1*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    (8*c3a1*l12*m1*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    c2a1*l1*m1*p1x*ralpha12 + 2*c3a1*l1*m1*p1x*ralpha12 + 
     1    2*c3a1*l12*p1y*ralpha12 - c2a1*l1*p1z*ralpha12 + 
     1    2*c3a1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
     1    (3*l1*m1*(n1*p1y - m1*p1z)*s1)*(ralpha12**2) - 
     1    (5*l12*m1*(l1*p1x + m1*p1y + n1*p1z)*s1)*(ralpha12**2) + 
     1    2*l1*m1*p1x*ralpha12*s1 + l12*p1y*ralpha12*s1 + 
     1    l1*p1z*ralpha12*s1 - l1*m1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
     1 2*(-(c3a1*p1x) - (3*c2a1*l1*m1*(-(n1*p1x) + l1*p1z))*(ralpha12**2) + 
     1    (c2a1*l1*m12*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    (8*c3a1*l1*m12*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) + 
     1    2*c3a1*m12*p1x*ralpha12 - c2a1*l1*m1*p1y*ralpha12 + 
     1    2*c3a1*l1*m1*p1y*ralpha12 + c2a1*m1*p1z*ralpha12 + 
     1    2*c3a1*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
     1    (3*l1*m1*(-(n1*p1x) + l1*p1z)*s1)*(ralpha12**2) - 
     1    (5*l1*m12*(l1*p1x + m1*p1y + n1*p1z)*s1)*(ralpha12**2) + 
     1    m12*p1x*ralpha12*s1 + 2*l1*m1*p1y*ralpha12*s1 - 
     1    m1*p1z*ralpha12*s1 - l1*m1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
     1 2*((-3*c2a1*l1*m1*(m1*p1x - l1*p1y))*(ralpha12**2) + 
     1    (c2a1*l1*m1*n1*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    (8*c3a1*l1*m1*n1*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) + 
     1    c2a1*l1*p1x*ralpha12 + 2*c3a1*m1*n1*p1x*ralpha12 - 
     1    c2a1*m1*p1y*ralpha12 + 2*c3a1*l1*n1*p1y*ralpha12 - 
     1    c2a1*l1*m1*p1z*ralpha12 + 
     1    (3*l1*m1*(m1*p1x - l1*p1y)*s1)*(ralpha12**2) - 
     1    (5*l1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*s1)*(ralpha12**2) - 
     1    l1*p1x*ralpha12*s1 + m1*n1*p1x*ralpha12*s1 + m1*p1y*ralpha12*s1 + 
     1    l1*n1*p1y*ralpha12*s1 - l1*m1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1    l1*m1*p1z*ralpha12*s1)*
     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
     1 -D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1  D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*rdist

       return
      end FUNCTION DD45

      DOUBLE PRECISION FUNCTION DD46(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D4, D6

      DD46  
     1 = (rdist*(2*(-(c3a1*l1*p1x) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1    c2a1*l1*(n1*p1y - m1*p1z)*ralpha12 + 
     1    2*c3a1*l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1x*s1 - 
     1    l1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (-(c3a1*l1*p1z) + c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1x*s1 + 
     1    p1y*s1 - n1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) + 
     1 2*(-(c3a1*m1*p1x) + c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1y*s1 + 
     1    p1z*s1 - l1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (-(c3a1*m1*p1z) + c2a1*n1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 - 
     1    n1*p1y*s1 - n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) + 
     1 2*(-(c3a1*n1*p1x) + c2a1*l1*(m1*p1x - l1*p1y)*ralpha12 + 
     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 - 
     1    l1*p1z*s1 - l1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (-(c3a1*n1*p1z) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1    c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
     1    2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1z*s1 - 
     1    n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1    n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) + 
     1 2*(-(c3a1*p1z) - (3*c2a1*l1*n1*(n1*p1y - m1*p1z))*(ralpha12**2) + 
     1    (c2a1*l12*n1*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    (8*c3a1*l12*n1*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    c2a1*l1*n1*p1x*ralpha12 + 2*c3a1*l1*n1*p1x*ralpha12 + 
     1    c2a1*l1*p1y*ralpha12 + 2*c3a1*l12*p1z*ralpha12 + 
     1    2*c3a1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
     1    (3*l1*n1*(n1*p1y - m1*p1z)*s1)*(ralpha12**2) - 
     1    (5*l12*n1*(l1*p1x + m1*p1y + n1*p1z)*s1)*(ralpha12**2) + 
     1    2*l1*n1*p1x*ralpha12*s1 - l1*p1y*ralpha12*s1 + 
     1    l12*p1z*ralpha12*s1 - l1*n1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
     1 2*((-3*c2a1*l1*n1*(-(n1*p1x) + l1*p1z))*(ralpha12**2) + 
     1    (c2a1*l1*m1*n1*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    (8*c3a1*l1*m1*n1*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    c2a1*l1*p1x*ralpha12 + 2*c3a1*m1*n1*p1x*ralpha12 - 
     1    c2a1*l1*n1*p1y*ralpha12 + 2*c3a1*l1*m1*p1z*ralpha12 + 
     1    c2a1*n1*p1z*ralpha12 + 
     1    (3*l1*n1*(-(n1*p1x) + l1*p1z)*s1)*(ralpha12**2) - 
     1    (5*l1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*s1)*(ralpha12**2) + 
     1    l1*p1x*ralpha12*s1 + m1*n1*p1x*ralpha12*s1 + 
     1    l1*n1*p1y*ralpha12*s1 + l1*m1*p1z*ralpha12*s1 - 
     1    n1*p1z*ralpha12*s1 - l1*n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1)*
     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
     1 2*(-(c3a1*p1x) - (3*c2a1*l1*n1*(m1*p1x - l1*p1y))*(ralpha12**2) + 
     1    (c2a1*l1*n12*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    (8*c3a1*l1*n12*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) + 
     1    2*c3a1*n12*p1x*ralpha12 - c2a1*n1*p1y*ralpha12 - 
     1    c2a1*l1*n1*p1z*ralpha12 + 2*c3a1*l1*n1*p1z*ralpha12 + 
     1    2*c3a1*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
     1    (3*l1*n1*(m1*p1x - l1*p1y)*s1)*(ralpha12**2) - 
     1    (5*l1*n12*(l1*p1x + m1*p1y + n1*p1z)*s1)*(ralpha12**2) + 
     1    n12*p1x*ralpha12*s1 + n1*p1y*ralpha12*s1 - 
     1    l1*n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 2*l1*n1*p1z*ralpha12*s1 + 
     1    l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
     1 -D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD46

      DOUBLE PRECISION FUNCTION DD47(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D4, D7

      DD47  
     1 = -(rdist*(-(c3a1*l1*p1x) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1 c2a1*l1*(n1*p1y - m1*p1z)*ralpha12 + 
     1 2*c3a1*l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1x*s1 - 
     1 l1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1 l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1))
     1 -D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD47

      DOUBLE PRECISION FUNCTION DD48(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D4, D8

      DD48  
     1 = -(rdist*(-(c3a1*m1*p1x) + c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1 2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1y*s1 + 
     1 p1z*s1 - l1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1 l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1))
     1 -D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD48

      DOUBLE PRECISION FUNCTION DD49(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D4, D9

      DD49  
     1 = -(rdist*(-(c3a1*n1*p1x) + c2a1*l1*(m1*p1x - l1*p1y)*ralpha12 + 
     1 2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 - 
     1 l1*p1z*s1 - l1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1 l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1))
     1 -D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD49

      DOUBLE PRECISION FUNCTION DD4A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D4, DA

      DD4A  
     1 = (rdist*(2*(-(c3a1*l1*p1x) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1    c2a1*l1*(n1*p1y - m1*p1z)*ralpha12 + 
     1    2*c3a1*l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1x*s1 - 
     1    l1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*l2*p2x + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1    c2a2*l2*(n2*p2y - m2*p2z)*ralpha22 - 
     1    2*c3a2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2x*s2 + 
     1    l2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(-(c3a1*m1*p1x) + c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1y*s1 + 
     1    p1z*s1 - l1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*m2*p2x - c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2y*s2 - 
     1    p2z*s2 + l2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(-(c3a1*n1*p1x) + c2a1*l1*(m1*p1x - l1*p1y)*ralpha12 + 
     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 - 
     1    l1*p1z*s1 - l1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*n2*p2x - c2a2*l2*(m2*p2x - l2*p2y)*ralpha22 - 
     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 + 
     1    l2*p2z*s2 + l2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)))/2.
     1 -D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD4A

      DOUBLE PRECISION FUNCTION DD4B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D4, DB

      DD4B  
     1 = (rdist*(2*(-(c3a1*l1*p1x) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1    c2a1*l1*(n1*p1y - m1*p1z)*ralpha12 + 
     1    2*c3a1*l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1x*s1 - 
     1    l1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*l2*p2y - c2a2*m2*(n2*p2y - m2*p2z)*ralpha22 - 
     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2x*s2 + 
     1    p2z*s2 + m2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(-(c3a1*m1*p1x) + c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1y*s1 + 
     1    p1z*s1 - l1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*m2*p2y + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1    c2a2*m2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1    2*c3a2*m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2y*s2 + 
     1    m2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(-(c3a1*n1*p1x) + c2a1*l1*(m1*p1x - l1*p1y)*ralpha12 + 
     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 - 
     1    l1*p1z*s1 - l1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*n2*p2y - c2a2*m2*(m2*p2x - l2*p2y)*ralpha22 - 
     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - p2x*s2 + 
     1    m2*p2z*s2 + m2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)))/2.
     1 -D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD4B

      DOUBLE PRECISION FUNCTION DD4C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D4, DC

      DD4C  
     1 = (rdist*(2*(-(c3a1*l1*p1x) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1    c2a1*l1*(n1*p1y - m1*p1z)*ralpha12 + 
     1    2*c3a1*l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1x*s1 - 
     1    l1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    l12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*l2*p2z - c2a2*n2*(n2*p2y - m2*p2z)*ralpha22 - 
     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2x*s2 - 
     1    p2y*s2 + n2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(-(c3a1*m1*p1x) + c2a1*l1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - l1*p1y*s1 + 
     1    p1z*s1 - l1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*m2*p2z - c2a2*n2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 + 
     1    n2*p2y*s2 + n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(-(c3a1*n1*p1x) + c2a1*l1*(m1*p1x - l1*p1y)*ralpha12 + 
     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 - 
     1    l1*p1z*s1 - l1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*n2*p2z + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1    c2a2*n2*(m2*p2x - l2*p2y)*ralpha22 - 
     1    2*c3a2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2z*s2 + 
     1    n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1    n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)))/2.
     1 -D4(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD4C

      DOUBLE PRECISION FUNCTION DD55(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D5

      DD55  
     1 = (rdist*(2*(-(c3a1*l1*p1y) + 
     1     c2a1*m1*(n1*p1y - m1*p1z)*ralpha12 + 
     1     2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1x*s1 - 
     1     p1z*s1 - m1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1     l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)**2 + 
     1 2*(-(c3a1*m1*p1y) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1     c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1     2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1y*s1 - 
     1     m1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1     m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)**2 + 
     1 2*(-(c3a1*n1*p1y) + c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
     1     2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + p1x*s1 - 
     1     m1*p1z*s1 - m1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1     m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)**2 + 
     1 2*((-3*c2a1*m12*(n1*p1y - m1*p1z))*(ralpha12**2) + 
     1    (c2a1*l1*m12*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    (8*c3a1*l1*m12*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    c2a1*m12*p1x*ralpha12 + 4*c3a1*l1*m1*p1y*ralpha12 - 
     1    2*c2a1*m1*p1z*ralpha12 + c2a1*(n1*p1y - m1*p1z)*ralpha12 + 
     1    2*c3a1*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 + 
     1    (3*m12*(n1*p1y - m1*p1z)*s1)*(ralpha12**2) - 
     1    (5*l1*m12*(l1*p1x + m1*p1y + n1*p1z)*s1)*(ralpha12**2) + 
     1    m12*p1x*ralpha12*s1 + 2*l1*m1*p1y*ralpha12*s1 + 
     1    2*m1*p1z*ralpha12*s1 - (n1*p1y - m1*p1z)*ralpha12*s1 - 
     1    m12*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
     1 2*(-2*c3a1*p1y - (3*c2a1*m12*(-(n1*p1x) + l1*p1z))*(ralpha12**2) + 
     1    (c2a1*m1**3*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    (8*c3a1*m1**3*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    c2a1*m12*p1y*ralpha12 + 4*c3a1*m12*p1y*ralpha12 + 
     1    c2a1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1    6*c3a1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 + 
     1    (3*m12*(-(n1*p1x) + l1*p1z)*s1)*(ralpha12**2) - 
     1    (5*m1**3*(l1*p1x + m1*p1y + n1*p1z)*s1)*(ralpha12**2) + 
     1    3*m12*p1y*ralpha12*s1 - (-(n1*p1x) + l1*p1z)*ralpha12*s1 - 
     1    m12*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    3*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
     1 2*((-3*c2a1*m12*(m1*p1x - l1*p1y))*(ralpha12**2) + 
     1    (c2a1*m12*n1*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    (8*c3a1*m12*n1*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) + 
     1    2*c2a1*m1*p1x*ralpha12 + 4*c3a1*m1*n1*p1y*ralpha12 + 
     1    c2a1*(m1*p1x - l1*p1y)*ralpha12 - c2a1*m12*p1z*ralpha12 + 
     1    2*c3a1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
     1    (3*m12*(m1*p1x - l1*p1y)*s1)*(ralpha12**2) - p1z*s1 - 
     1    (5*m12*n1*(l1*p1x + m1*p1y + n1*p1z)*s1)*(ralpha12**2) - 
     1    2*m1*p1x*ralpha12*s1 + 2*m1*n1*p1y*ralpha12*s1 - 
     1    (m1*p1x - l1*p1y)*ralpha12*s1 - 
     1    m12*(m1*p1x - l1*p1y)*ralpha12*s1 + m12*p1z*ralpha12*s1 + 
     1    n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
     1 -D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD55

      DOUBLE PRECISION FUNCTION DD56(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D5, D6

      DD56  
     1 = (rdist*(2*(-(c3a1*l1*p1y) + 
     1    c2a1*m1*(n1*p1y - m1*p1z)*ralpha12 + 
     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1x*s1 - 
     1    p1z*s1 - m1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (-(c3a1*l1*p1z) + c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1x*s1 + 
     1    p1y*s1 - n1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) + 
     1 2*(-(c3a1*m1*p1y) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1    c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1    2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1y*s1 - 
     1    m1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (-(c3a1*m1*p1z) + c2a1*n1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 - 
     1    n1*p1y*s1 - n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) + 
     1 2*(-(c3a1*n1*p1y) + c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + p1x*s1 - 
     1    m1*p1z*s1 - m1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (-(c3a1*n1*p1z) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1    c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
     1    2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1z*s1 - 
     1    n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1    n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1) + 
     1 2*((-3*c2a1*m1*n1*(n1*p1y - m1*p1z))*(ralpha12**2) + 
     1    (c2a1*l1*m1*n1*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    (8*c3a1*l1*m1*n1*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    c2a1*m1*n1*p1x*ralpha12 + c2a1*m1*p1y*ralpha12 + 
     1    2*c3a1*l1*n1*p1y*ralpha12 + 2*c3a1*l1*m1*p1z*ralpha12 - 
     1    c2a1*n1*p1z*ralpha12 + 
     1    (3*m1*n1*(n1*p1y - m1*p1z)*s1)*(ralpha12**2) - 
     1    (5*l1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*s1)*(ralpha12**2) + 
     1    m1*n1*p1x*ralpha12*s1 - m1*p1y*ralpha12*s1 + 
     1    l1*n1*p1y*ralpha12*s1 + l1*m1*p1z*ralpha12*s1 + 
     1    n1*p1z*ralpha12*s1 - m1*n1*(n1*p1y - m1*p1z)*ralpha12*s1)*
     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
     1 2*(-(c3a1*p1z) - (3*c2a1*m1*n1*(-(n1*p1x) + l1*p1z))*(ralpha12**2) + 
     1    (c2a1*m12*n1*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    (8*c3a1*m12*n1*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    c2a1*m1*p1x*ralpha12 - c2a1*m1*n1*p1y*ralpha12 + 
     1    2*c3a1*m1*n1*p1y*ralpha12 + 2*c3a1*m12*p1z*ralpha12 + 
     1    2*c3a1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
     1    (3*m1*n1*(-(n1*p1x) + l1*p1z)*s1)*(ralpha12**2) - 
     1    (5*m12*n1*(l1*p1x + m1*p1y + n1*p1z)*s1)*(ralpha12**2) + 
     1    m1*p1x*ralpha12*s1 + 2*m1*n1*p1y*ralpha12*s1 + 
     1    m12*p1z*ralpha12*s1 - m1*n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
     1 2*(-(c3a1*p1y) - (3*c2a1*m1*n1*(m1*p1x - l1*p1y))*(ralpha12**2) + 
     1    (c2a1*m1*n12*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    (8*c3a1*m1*n12*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) + 
     1    c2a1*n1*p1x*ralpha12 + 2*c3a1*n12*p1y*ralpha12 - 
     1    c2a1*m1*n1*p1z*ralpha12 + 2*c3a1*m1*n1*p1z*ralpha12 + 
     1    2*c3a1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
     1    (3*m1*n1*(m1*p1x - l1*p1y)*s1)*(ralpha12**2) - 
     1    (5*m1*n12*(l1*p1x + m1*p1y + n1*p1z)*s1)*(ralpha12**2) - 
     1    n1*p1x*ralpha12*s1 + n12*p1y*ralpha12*s1 - 
     1    m1*n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 2*m1*n1*p1z*ralpha12*s1 + 
     1    m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
     1 -D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD56

      DOUBLE PRECISION FUNCTION DD57(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D5, D7

      DD57  
     1 = -(rdist*(-(c3a1*l1*p1y) + c2a1*m1*(n1*p1y - m1*p1z)*ralpha12 + 
     1 2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1x*s1 - 
     1 p1z*s1 - m1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1 l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1))
     1 -D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD57

      DOUBLE PRECISION FUNCTION DD58(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D5, D8

      DD58  
     1 = -(rdist*(-(c3a1*m1*p1y) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1 c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1 2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1y*s1 - 
     1 m1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1 m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1))
     1 -D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD58

      DOUBLE PRECISION FUNCTION DD59(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D5, D9

      DD59  
     1 = -(rdist*(-(c3a1*n1*p1y) + c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
     1 2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + p1x*s1 - 
     1 m1*p1z*s1 - m1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1 m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1))
     1 -D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD59

      DOUBLE PRECISION FUNCTION DD5A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D5, DA

      DD5A  
     1 = (rdist*(2*(-(c3a1*l1*p1y) + 
     1    c2a1*m1*(n1*p1y - m1*p1z)*ralpha12 + 
     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1x*s1 - 
     1    p1z*s1 - m1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*l2*p2x + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1    c2a2*l2*(n2*p2y - m2*p2z)*ralpha22 - 
     1    2*c3a2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2x*s2 + 
     1    l2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(-(c3a1*m1*p1y) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1    c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1    2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1y*s1 - 
     1    m1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*m2*p2x - c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2y*s2 - 
     1    p2z*s2 + l2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(-(c3a1*n1*p1y) + c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + p1x*s1 - 
     1    m1*p1z*s1 - m1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*n2*p2x - c2a2*l2*(m2*p2x - l2*p2y)*ralpha22 - 
     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 + 
     1    l2*p2z*s2 + l2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)))/2.
     1 -D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD5A

      DOUBLE PRECISION FUNCTION DD5B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D5, DB

      DD5B  
     1 = (rdist*(2*(-(c3a1*l1*p1y) + 
     1    c2a1*m1*(n1*p1y - m1*p1z)*ralpha12 + 
     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1x*s1 - 
     1    p1z*s1 - m1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*l2*p2y - c2a2*m2*(n2*p2y - m2*p2z)*ralpha22 - 
     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2x*s2 + 
     1    p2z*s2 + m2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(-(c3a1*m1*p1y) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1    c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1    2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1y*s1 - 
     1    m1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*m2*p2y + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1    c2a2*m2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1    2*c3a2*m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2y*s2 + 
     1    m2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(-(c3a1*n1*p1y) + c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + p1x*s1 - 
     1    m1*p1z*s1 - m1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*n2*p2y - c2a2*m2*(m2*p2x - l2*p2y)*ralpha22 - 
     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - p2x*s2 + 
     1    m2*p2z*s2 + m2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)))/2.
     1 -D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD5B

      DOUBLE PRECISION FUNCTION DD5C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D5, DC

      DD5C  
     1 = (rdist*(2*(-(c3a1*l1*p1y) + 
     1    c2a1*m1*(n1*p1y - m1*p1z)*ralpha12 + 
     1    2*c3a1*l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1x*s1 - 
     1    p1z*s1 - m1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    l1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*l2*p2z - c2a2*n2*(n2*p2y - m2*p2z)*ralpha22 - 
     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2x*s2 - 
     1    p2y*s2 + n2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(-(c3a1*m1*p1y) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1    c2a1*m1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1    2*c3a1*m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - m1*p1y*s1 - 
     1    m1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    m12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*m2*p2z - c2a2*n2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 + 
     1    n2*p2y*s2 + n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(-(c3a1*n1*p1y) + c2a1*m1*(m1*p1x - l1*p1y)*ralpha12 + 
     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + p1x*s1 - 
     1    m1*p1z*s1 - m1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*n2*p2z + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1    c2a2*n2*(m2*p2x - l2*p2y)*ralpha22 - 
     1    2*c3a2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2z*s2 + 
     1    n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1    n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)))/2.
     1 -D5(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD5C

      DOUBLE PRECISION FUNCTION DD66(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D6

      DD66  
     1 = (rdist*(2*(-(c3a1*l1*p1z) + 
     1     c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
     1     2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1x*s1 + 
     1     p1y*s1 - n1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1     l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)**2 + 
     1 2*(-(c3a1*m1*p1z) + c2a1*n1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1     2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 - 
     1     n1*p1y*s1 - n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1     m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)**2 + 
     1 2*(-(c3a1*n1*p1z) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1     c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
     1     2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1z*s1 - 
     1     n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1     n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)**2 + 
     1 2*((-3*c2a1*n12*(n1*p1y - m1*p1z))*(ralpha12**2) + 
     1    (c2a1*l1*n12*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    (8*c3a1*l1*n12*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    c2a1*n12*p1x*ralpha12 + 2*c2a1*n1*p1y*ralpha12 + 
     1    4*c3a1*l1*n1*p1z*ralpha12 + c2a1*(n1*p1y - m1*p1z)*ralpha12 + 
     1    2*c3a1*l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 + 
     1    (3*n12*(n1*p1y - m1*p1z)*s1)*(ralpha12**2) - 
     1    (5*l1*n12*(l1*p1x + m1*p1y + n1*p1z)*s1)*(ralpha12**2) + 
     1    n12*p1x*ralpha12*s1 - 2*n1*p1y*ralpha12*s1 + 
     1    2*l1*n1*p1z*ralpha12*s1 - (n1*p1y - m1*p1z)*ralpha12*s1 - 
     1    n12*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    l1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
     1 2*((-3*c2a1*n12*(-(n1*p1x) + l1*p1z))*(ralpha12**2) + 
     1    (c2a1*m1*n12*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    (8*c3a1*m1*n12*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    2*c2a1*n1*p1x*ralpha12 - c2a1*n12*p1y*ralpha12 + 
     1    4*c3a1*m1*n1*p1z*ralpha12 + c2a1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1    2*c3a1*m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1y*s1 + 
     1    (3*n12*(-(n1*p1x) + l1*p1z)*s1)*(ralpha12**2) - 
     1    (5*m1*n12*(l1*p1x + m1*p1y + n1*p1z)*s1)*(ralpha12**2) + 
     1    2*n1*p1x*ralpha12*s1 + n12*p1y*ralpha12*s1 + 
     1    2*m1*n1*p1z*ralpha12*s1 - (-(n1*p1x) + l1*p1z)*ralpha12*s1 - 
     1    n12*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    m1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
     1 2*((-3*c2a1*n12*(m1*p1x - l1*p1y))*(ralpha12**2) - 2*c3a1*p1z + 
     1    (c2a1*n1**3*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) - 
     1    (8*c3a1*n1**3*(l1*p1x + m1*p1y + n1*p1z))*(ralpha12**2) + 
     1    c2a1*(m1*p1x - l1*p1y)*ralpha12 - c2a1*n12*p1z*ralpha12 + 
     1    4*c3a1*n12*p1z*ralpha12 + 
     1    6*c3a1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 + 
     1    (3*n12*(m1*p1x - l1*p1y)*s1)*(ralpha12**2) - p1z*s1 - 
     1    (5*n1**3*(l1*p1x + m1*p1y + n1*p1z)*s1)*(ralpha12**2) - 
     1    (m1*p1x - l1*p1y)*ralpha12*s1 - 
     1    n12*(m1*p1x - l1*p1y)*ralpha12*s1 + 3*n12*p1z*ralpha12*s1 + 
     1    3*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
     1 -D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD66

      DOUBLE PRECISION FUNCTION DD67(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D6, D7

      DD67  
     1 = -(rdist*(-(c3a1*l1*p1z) + c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
     1 2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1x*s1 + 
     1 p1y*s1 - n1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1 l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1))
     1 -D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD67

      DOUBLE PRECISION FUNCTION DD68(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D6, D8

      DD68  
     1 = -(rdist*(-(c3a1*m1*p1z) + c2a1*n1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1 2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 - 
     1 n1*p1y*s1 - n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1 m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1))
     1 -D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD68

      DOUBLE PRECISION FUNCTION DD69(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D6, D9

      DD69  
     1 = -(rdist*(-(c3a1*n1*p1z) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1 c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
     1 2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1z*s1 - 
     1 n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1 n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1))
     1 -D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD69

      DOUBLE PRECISION FUNCTION DD6A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D6, DA

      DD6A  
     1 = (rdist*(2*(-(c3a1*l1*p1z) + 
     1    c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1x*s1 + 
     1    p1y*s1 - n1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*l2*p2x + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1    c2a2*l2*(n2*p2y - m2*p2z)*ralpha22 - 
     1    2*c3a2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2x*s2 + 
     1    l2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(-(c3a1*m1*p1z) + c2a1*n1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 - 
     1    n1*p1y*s1 - n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*m2*p2x - c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2y*s2 - 
     1    p2z*s2 + l2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(-(c3a1*n1*p1z) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1    c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
     1    2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1z*s1 - 
     1    n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1    n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*n2*p2x - c2a2*l2*(m2*p2x - l2*p2y)*ralpha22 - 
     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 + 
     1    l2*p2z*s2 + l2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)))/2.
     1 -D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD6A

      DOUBLE PRECISION FUNCTION DD6B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D6, DB

      DD6B  
     1 = (rdist*(2*(-(c3a1*l1*p1z) + 
     1    c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1x*s1 + 
     1    p1y*s1 - n1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*l2*p2y - c2a2*m2*(n2*p2y - m2*p2z)*ralpha22 - 
     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2x*s2 + 
     1    p2z*s2 + m2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(-(c3a1*m1*p1z) + c2a1*n1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 - 
     1    n1*p1y*s1 - n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*m2*p2y + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1    c2a2*m2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1    2*c3a2*m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2y*s2 + 
     1    m2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(-(c3a1*n1*p1z) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1    c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
     1    2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1z*s1 - 
     1    n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1    n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*n2*p2y - c2a2*m2*(m2*p2x - l2*p2y)*ralpha22 - 
     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - p2x*s2 + 
     1    m2*p2z*s2 + m2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)))/2.
     1 -D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD6B

      DOUBLE PRECISION FUNCTION DD6C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D6, DC

      DD6C  
     1 = (rdist*(2*(-(c3a1*l1*p1z) + 
     1    c2a1*n1*(n1*p1y - m1*p1z)*ralpha12 + 
     1    2*c3a1*l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1x*s1 + 
     1    p1y*s1 - n1*(n1*p1y - m1*p1z)*ralpha12*s1 + 
     1    l1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*l2*p2z - c2a2*n2*(n2*p2y - m2*p2z)*ralpha22 - 
     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2x*s2 - 
     1    p2y*s2 + n2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(-(c3a1*m1*p1z) + c2a1*n1*(-(n1*p1x) + l1*p1z)*ralpha12 + 
     1    2*c3a1*m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - p1x*s1 - 
     1    n1*p1y*s1 - n1*(-(n1*p1x) + l1*p1z)*ralpha12*s1 + 
     1    m1*n1*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*m2*p2z - c2a2*n2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 + 
     1    n2*p2y*s2 + n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(-(c3a1*n1*p1z) - c3a1*(l1*p1x + m1*p1y + n1*p1z) + 
     1    c2a1*n1*(m1*p1x - l1*p1y)*ralpha12 + 
     1    2*c3a1*n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12 - n1*p1z*s1 - 
     1    n1*(m1*p1x - l1*p1y)*ralpha12*s1 + 
     1    n12*(l1*p1x + m1*p1y + n1*p1z)*ralpha12*s1)*
     1  (c3a2*n2*p2z + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1    c2a2*n2*(m2*p2x - l2*p2y)*ralpha22 - 
     1    2*c3a2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2z*s2 + 
     1    n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1    n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)))/2.
     1 -D6(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD6C

      DOUBLE PRECISION FUNCTION DD77(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D7

       DD77 = rdist*(1.0d0 -
     1               D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)**2)

       return
      end FUNCTION DD77

      DOUBLE PRECISION FUNCTION DD78(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D7, D8

       DD78 = -rdist *
     1        D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1        D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)

       return
      end FUNCTION DD78

      DOUBLE PRECISION FUNCTION DD79(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D7, D9

       DD79 = -rdist *
     1        D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1        D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)

       return
      end FUNCTION DD79

      DOUBLE PRECISION FUNCTION DD7A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D7, DA

      DD7A  
     1 = -(rdist*(c3a2*l2*p2x + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1 c2a2*l2*(n2*p2y - m2*p2z)*ralpha22 - 
     1 2*c3a2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2x*s2 + 
     1 l2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1 l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2))
     1 -D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD7A

      DOUBLE PRECISION FUNCTION DD7B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D7, DB

      DD7B  
     1 = -(rdist*(c3a2*l2*p2y - c2a2*m2*(n2*p2y - m2*p2z)*ralpha22 - 
     1 2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2x*s2 + 
     1 p2z*s2 + m2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1 l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2))
     1 -D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD7B

      DOUBLE PRECISION FUNCTION DD7C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D7, DC

      DD7C  
     1 = -(rdist*(c3a2*l2*p2z - c2a2*n2*(n2*p2y - m2*p2z)*ralpha22 - 
     1 2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2x*s2 - 
     1 p2y*s2 + n2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1 l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2))
     1 -D7(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD7C

      DOUBLE PRECISION FUNCTION DD88(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D8

       DD88 = rdist*(1.0d0 -
     1               D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)**2)

       return
      end FUNCTION DD88

      DOUBLE PRECISION FUNCTION DD89(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D8, D9

       DD89 = -rdist *
     1        D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1        D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1           C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)

       return
      end FUNCTION DD89

      DOUBLE PRECISION FUNCTION DD8A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D8, DA

      DD8A  
     1 = -(rdist*(c3a2*m2*p2x - c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1 2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2y*s2 - 
     1 p2z*s2 + l2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1 l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2))
     1 -D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD8A

      DOUBLE PRECISION FUNCTION DD8B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D8, DB

      DD8B  
     1 = -(rdist*(c3a2*m2*p2y + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1 c2a2*m2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1 2*c3a2*m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2y*s2 + 
     1 m2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1 m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2))
     1 -D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD8B
      
      DOUBLE PRECISION FUNCTION DD8C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D8, DC

      DD8C  
     1 = -(rdist*(c3a2*m2*p2z - c2a2*n2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1 2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 + 
     1 n2*p2y*s2 + n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1 m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2))
     1 -D8(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD8C

      DOUBLE PRECISION FUNCTION DD99(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D9

       DD99 = rdist*(1.0d0 -
     1               D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                  C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)**2)

       return
      end FUNCTION DD99

      DOUBLE PRECISION FUNCTION DD9A(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D9, DA

      DD9A  
     1 = -(rdist*(c3a2*n2*p2x - c2a2*l2*(m2*p2x - l2*p2y)*ralpha22 - 
     1 2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 + 
     1 l2*p2z*s2 + l2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1 l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2))
     1 -D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD9A

      DOUBLE PRECISION FUNCTION DD9B(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D9, DB

      DD9B  
     1 = -(rdist*(c3a2*n2*p2y - c2a2*m2*(m2*p2x - l2*p2y)*ralpha22 - 
     1 2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - p2x*s2 + 
     1 m2*p2z*s2 + m2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1 m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2))
     1 -D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD9B

      DOUBLE PRECISION FUNCTION DD9C(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION D9, DC

      DD9C  
     1 = -(rdist*(c3a2*n2*p2z + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1 c2a2*n2*(m2*p2x - l2*p2y)*ralpha22 - 
     1 2*c3a2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2z*s2 + 
     1 n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1 n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2))
     1 -D9(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DD9C

      DOUBLE PRECISION FUNCTION DDAA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION DA

      DDAA  
     1 = (rdist*(2*(c3a2*l2*p2x + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1     c2a2*l2*(n2*p2y - m2*p2z)*ralpha22 - 
     1     2*c3a2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2x*s2 + 
     1     l2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1     l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)**2 + 
     1 2*(c3a2*m2*p2x - c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1     2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2y*s2 - 
     1     p2z*s2 + l2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1     l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)**2 + 
     1 2*(c3a2*n2*p2x - c2a2*l2*(m2*p2x - l2*p2y)*ralpha22 - 
     1     2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 + 
     1     l2*p2z*s2 + l2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1     l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)**2 + 
     1 2*(2*c3a2*p2x + (3*c2a2*l22*(n2*p2y - m2*p2z))*(ralpha22**2) - 
     1    (c2a2*l2**3*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    (8*c3a2*l2**3*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    c2a2*l22*p2x*ralpha22 - 4*c3a2*l22*p2x*ralpha22 - 
     1    c2a2*(n2*p2y - m2*p2z)*ralpha22 - 
     1    6*c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 - 
     1    (3*l22*(n2*p2y - m2*p2z)*s2)*(ralpha22**2) + 
     1    (5*l2**3*(l2*p2x + m2*p2y + n2*p2z)*s2)*(ralpha22**2) - 
     1    3*l22*p2x*ralpha22*s2 + (n2*p2y - m2*p2z)*ralpha22*s2 + 
     1    l22*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    3*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
     1 2*((3*c2a2*l22*(-(n2*p2x) + l2*p2z))*(ralpha22**2) - 
     1    (c2a2*l22*m2*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    (8*c3a2*l22*m2*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) - 
     1    4*c3a2*l2*m2*p2x*ralpha22 + c2a2*l22*p2y*ralpha22 - 
     1    2*c2a2*l2*p2z*ralpha22 - c2a2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1    2*c3a2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 - 
     1    (3*l22*(-(n2*p2x) + l2*p2z)*s2)*(ralpha22**2) + 
     1    (5*l22*m2*(l2*p2x + m2*p2y + n2*p2z)*s2)*(ralpha22**2) - 
     1    2*l2*m2*p2x*ralpha22*s2 - l22*p2y*ralpha22*s2 + 
     1    2*l2*p2z*ralpha22*s2 + (-(n2*p2x) + l2*p2z)*ralpha22*s2 + 
     1    l22*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
     1 2*((3*c2a2*l22*(m2*p2x - l2*p2y))*(ralpha22**2) - 
     1    (c2a2*l22*n2*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    (8*c3a2*l22*n2*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) - 
     1    4*c3a2*l2*n2*p2x*ralpha22 + 2*c2a2*l2*p2y*ralpha22 - 
     1    c2a2*(m2*p2x - l2*p2y)*ralpha22 + c2a2*l22*p2z*ralpha22 - 
     1    2*c3a2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - 
     1    (3*l22*(m2*p2x - l2*p2y)*s2)*(ralpha22**2) + p2z*s2 + 
     1    (5*l22*n2*(l2*p2x + m2*p2y + n2*p2z)*s2)*(ralpha22**2) - 
     1    2*l2*n2*p2x*ralpha22*s2 - 2*l2*p2y*ralpha22*s2 + 
     1    (m2*p2x - l2*p2y)*ralpha22*s2 + 
     1    l22*(m2*p2x - l2*p2y)*ralpha22*s2 - l22*p2z*ralpha22*s2 - 
     1    n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
     1 -DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DDAA

      DOUBLE PRECISION FUNCTION DDAB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION DA, DB

      DDAB  
     1 = (rdist*(2*(c3a2*l2*p2x + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1    c2a2*l2*(n2*p2y - m2*p2z)*ralpha22 - 
     1    2*c3a2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2x*s2 + 
     1    l2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c3a2*l2*p2y - c2a2*m2*(n2*p2y - m2*p2z)*ralpha22 - 
     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2x*s2 + 
     1    p2z*s2 + m2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(c3a2*m2*p2x - c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2y*s2 - 
     1    p2z*s2 + l2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c3a2*m2*p2y + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1    c2a2*m2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1    2*c3a2*m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2y*s2 + 
     1    m2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(c3a2*n2*p2x - c2a2*l2*(m2*p2x - l2*p2y)*ralpha22 - 
     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 + 
     1    l2*p2z*s2 + l2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c3a2*n2*p2y - c2a2*m2*(m2*p2x - l2*p2y)*ralpha22 - 
     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - p2x*s2 + 
     1    m2*p2z*s2 + m2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(c3a2*p2y + (3*c2a2*l2*m2*(n2*p2y - m2*p2z))*(ralpha22**2) - 
     1    (c2a2*l22*m2*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    (8*c3a2*l22*m2*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    c2a2*l2*m2*p2x*ralpha22 - 2*c3a2*l2*m2*p2x*ralpha22 - 
     1    2*c3a2*l22*p2y*ralpha22 + c2a2*l2*p2z*ralpha22 - 
     1    2*c3a2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - 
     1    (3*l2*m2*(n2*p2y - m2*p2z)*s2)*(ralpha22**2) + 
     1    (5*l22*m2*(l2*p2x + m2*p2y + n2*p2z)*s2)*(ralpha22**2) - 
     1    2*l2*m2*p2x*ralpha22*s2 - l22*p2y*ralpha22*s2 - 
     1    l2*p2z*ralpha22*s2 + l2*m2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
     1 2*(c3a2*p2x + (3*c2a2*l2*m2*(-(n2*p2x) + l2*p2z))*(ralpha22**2) - 
     1    (c2a2*l2*m22*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    (8*c3a2*l2*m22*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) - 
     1    2*c3a2*m22*p2x*ralpha22 + c2a2*l2*m2*p2y*ralpha22 - 
     1    2*c3a2*l2*m2*p2y*ralpha22 - c2a2*m2*p2z*ralpha22 - 
     1    2*c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - 
     1    (3*l2*m2*(-(n2*p2x) + l2*p2z)*s2)*(ralpha22**2) + 
     1    (5*l2*m22*(l2*p2x + m2*p2y + n2*p2z)*s2)*(ralpha22**2) - 
     1    m22*p2x*ralpha22*s2 - 2*l2*m2*p2y*ralpha22*s2 + 
     1    m2*p2z*ralpha22*s2 + l2*m2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
     1 2*((3*c2a2*l2*m2*(m2*p2x - l2*p2y))*(ralpha22**2) - 
     1    (c2a2*l2*m2*n2*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    (8*c3a2*l2*m2*n2*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) - 
     1    c2a2*l2*p2x*ralpha22 - 2*c3a2*m2*n2*p2x*ralpha22 + 
     1    c2a2*m2*p2y*ralpha22 - 2*c3a2*l2*n2*p2y*ralpha22 + 
     1    c2a2*l2*m2*p2z*ralpha22 - 
     1    (3*l2*m2*(m2*p2x - l2*p2y)*s2)*(ralpha22**2) + 
     1    (5*l2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*s2)*(ralpha22**2) + 
     1    l2*p2x*ralpha22*s2 - m2*n2*p2x*ralpha22*s2 - m2*p2y*ralpha22*s2 - 
     1    l2*n2*p2y*ralpha22*s2 + l2*m2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1    l2*m2*p2z*ralpha22*s2)*
     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
     1 -DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DDAB

      DOUBLE PRECISION FUNCTION DDAC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION DA, DC

      DDAC  
     1 = (rdist*(2*(c3a2*l2*p2x + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1    c2a2*l2*(n2*p2y - m2*p2z)*ralpha22 - 
     1    2*c3a2*l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2x*s2 + 
     1    l2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    l22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c3a2*l2*p2z - c2a2*n2*(n2*p2y - m2*p2z)*ralpha22 - 
     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2x*s2 - 
     1    p2y*s2 + n2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(c3a2*m2*p2x - c2a2*l2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + l2*p2y*s2 - 
     1    p2z*s2 + l2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c3a2*m2*p2z - c2a2*n2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 + 
     1    n2*p2y*s2 + n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(c3a2*n2*p2x - c2a2*l2*(m2*p2x - l2*p2y)*ralpha22 - 
     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 + 
     1    l2*p2z*s2 + l2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c3a2*n2*p2z + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1    c2a2*n2*(m2*p2x - l2*p2y)*ralpha22 - 
     1    2*c3a2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2z*s2 + 
     1    n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1    n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(c3a2*p2z + (3*c2a2*l2*n2*(n2*p2y - m2*p2z))*(ralpha22**2) - 
     1    (c2a2*l22*n2*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    (8*c3a2*l22*n2*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    c2a2*l2*n2*p2x*ralpha22 - 2*c3a2*l2*n2*p2x*ralpha22 - 
     1    c2a2*l2*p2y*ralpha22 - 2*c3a2*l22*p2z*ralpha22 - 
     1    2*c3a2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - 
     1    (3*l2*n2*(n2*p2y - m2*p2z)*s2)*(ralpha22**2) + 
     1    (5*l22*n2*(l2*p2x + m2*p2y + n2*p2z)*s2)*(ralpha22**2) - 
     1    2*l2*n2*p2x*ralpha22*s2 + l2*p2y*ralpha22*s2 - 
     1    l22*p2z*ralpha22*s2 + l2*n2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
     1 2*((3*c2a2*l2*n2*(-(n2*p2x) + l2*p2z))*(ralpha22**2) - 
     1    (c2a2*l2*m2*n2*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    (8*c3a2*l2*m2*n2*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    c2a2*l2*p2x*ralpha22 - 2*c3a2*m2*n2*p2x*ralpha22 + 
     1    c2a2*l2*n2*p2y*ralpha22 - 2*c3a2*l2*m2*p2z*ralpha22 - 
     1    c2a2*n2*p2z*ralpha22 - 
     1    (3*l2*n2*(-(n2*p2x) + l2*p2z)*s2)*(ralpha22**2) + 
     1    (5*l2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*s2)*(ralpha22**2) - 
     1    l2*p2x*ralpha22*s2 - m2*n2*p2x*ralpha22*s2 - 
     1    l2*n2*p2y*ralpha22*s2 - l2*m2*p2z*ralpha22*s2 + 
     1    n2*p2z*ralpha22*s2 + l2*n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2)*
     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
     1 2*(c3a2*p2x + (3*c2a2*l2*n2*(m2*p2x - l2*p2y))*(ralpha22**2) - 
     1    (c2a2*l2*n22*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    (8*c3a2*l2*n22*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) - 
     1    2*c3a2*n22*p2x*ralpha22 + c2a2*n2*p2y*ralpha22 + 
     1    c2a2*l2*n2*p2z*ralpha22 - 2*c3a2*l2*n2*p2z*ralpha22 - 
     1    2*c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - 
     1    (3*l2*n2*(m2*p2x - l2*p2y)*s2)*(ralpha22**2) + 
     1    (5*l2*n22*(l2*p2x + m2*p2y + n2*p2z)*s2)*(ralpha22**2) - 
     1    n22*p2x*ralpha22*s2 - n2*p2y*ralpha22*s2 + 
     1    l2*n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 2*l2*n2*p2z*ralpha22*s2 - 
     1    l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
     1 -DA(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DDAC

      DOUBLE PRECISION FUNCTION DDBB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION DB

      DDBB  
     1 = (rdist*(2*(c3a2*l2*p2y - c2a2*m2*(n2*p2y - m2*p2z)*ralpha22 - 
     1     2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2x*s2 + 
     1     p2z*s2 + m2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1     l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)**2 + 
     1 2*(c3a2*m2*p2y + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1     c2a2*m2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1     2*c3a2*m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2y*s2 + 
     1     m2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1     m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)**2 + 
     1 2*(c3a2*n2*p2y - c2a2*m2*(m2*p2x - l2*p2y)*ralpha22 - 
     1     2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - p2x*s2 + 
     1     m2*p2z*s2 + m2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1     m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)**2 + 
     1 2*((3*c2a2*m22*(n2*p2y - m2*p2z))*(ralpha22**2) - 
     1    (c2a2*l2*m22*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    (8*c3a2*l2*m22*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    c2a2*m22*p2x*ralpha22 - 4*c3a2*l2*m2*p2y*ralpha22 + 
     1    2*c2a2*m2*p2z*ralpha22 - c2a2*(n2*p2y - m2*p2z)*ralpha22 - 
     1    2*c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 - 
     1    (3*m22*(n2*p2y - m2*p2z)*s2)*(ralpha22**2) + 
     1    (5*l2*m22*(l2*p2x + m2*p2y + n2*p2z)*s2)*(ralpha22**2) - 
     1    m22*p2x*ralpha22*s2 - 2*l2*m2*p2y*ralpha22*s2 - 
     1    2*m2*p2z*ralpha22*s2 + (n2*p2y - m2*p2z)*ralpha22*s2 + 
     1    m22*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
     1 2*(2*c3a2*p2y + (3*c2a2*m22*(-(n2*p2x) + l2*p2z))*(ralpha22**2) - 
     1    (c2a2*m2**3*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    (8*c3a2*m2**3*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    c2a2*m22*p2y*ralpha22 - 4*c3a2*m22*p2y*ralpha22 - 
     1    c2a2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1    6*c3a2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 - 
     1    (3*m22*(-(n2*p2x) + l2*p2z)*s2)*(ralpha22**2) + 
     1    (5*m2**3*(l2*p2x + m2*p2y + n2*p2z)*s2)*(ralpha22**2) - 
     1    3*m22*p2y*ralpha22*s2 + (-(n2*p2x) + l2*p2z)*ralpha22*s2 + 
     1    m22*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    3*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
     1 2*((3*c2a2*m22*(m2*p2x - l2*p2y))*(ralpha22**2) - 
     1    (c2a2*m22*n2*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    (8*c3a2*m22*n2*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) - 
     1    2*c2a2*m2*p2x*ralpha22 - 4*c3a2*m2*n2*p2y*ralpha22 - 
     1    c2a2*(m2*p2x - l2*p2y)*ralpha22 + c2a2*m22*p2z*ralpha22 - 
     1    2*c3a2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - 
     1    (3*m22*(m2*p2x - l2*p2y)*s2)*(ralpha22**2) + p2z*s2 + 
     1    (5*m22*n2*(l2*p2x + m2*p2y + n2*p2z)*s2)*(ralpha22**2) + 
     1    2*m2*p2x*ralpha22*s2 - 2*m2*n2*p2y*ralpha22*s2 + 
     1    (m2*p2x - l2*p2y)*ralpha22*s2 + 
     1    m22*(m2*p2x - l2*p2y)*ralpha22*s2 - m22*p2z*ralpha22*s2 - 
     1    n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
     1 -DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DDBB

      DOUBLE PRECISION FUNCTION DDBC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION DB, DC

      DDBC  
     1 = (rdist*(2*(c3a2*l2*p2y - c2a2*m2*(n2*p2y - m2*p2z)*ralpha22 - 
     1    2*c3a2*l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2x*s2 + 
     1    p2z*s2 + m2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    l2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c3a2*l2*p2z - c2a2*n2*(n2*p2y - m2*p2z)*ralpha22 - 
     1    2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2x*s2 - 
     1    p2y*s2 + n2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(c3a2*m2*p2y + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1    c2a2*m2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1    2*c3a2*m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + m2*p2y*s2 + 
     1    m2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    m22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c3a2*m2*p2z - c2a2*n2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 + 
     1    n2*p2y*s2 + n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*(c3a2*n2*p2y - c2a2*m2*(m2*p2x - l2*p2y)*ralpha22 - 
     1    2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - p2x*s2 + 
     1    m2*p2z*s2 + m2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1    m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c3a2*n2*p2z + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1    c2a2*n2*(m2*p2x - l2*p2y)*ralpha22 - 
     1    2*c3a2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2z*s2 + 
     1    n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1    n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2) + 
     1 2*((3*c2a2*m2*n2*(n2*p2y - m2*p2z))*(ralpha22**2) - 
     1    (c2a2*l2*m2*n2*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    (8*c3a2*l2*m2*n2*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    c2a2*m2*n2*p2x*ralpha22 - c2a2*m2*p2y*ralpha22 - 
     1    2*c3a2*l2*n2*p2y*ralpha22 - 2*c3a2*l2*m2*p2z*ralpha22 + 
     1    c2a2*n2*p2z*ralpha22 - 
     1    (3*m2*n2*(n2*p2y - m2*p2z)*s2)*(ralpha22**2) + 
     1    (5*l2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*s2)*(ralpha22**2) - 
     1    m2*n2*p2x*ralpha22*s2 + m2*p2y*ralpha22*s2 - 
     1    l2*n2*p2y*ralpha22*s2 - l2*m2*p2z*ralpha22*s2 - 
     1    n2*p2z*ralpha22*s2 + m2*n2*(n2*p2y - m2*p2z)*ralpha22*s2)*
     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
     1 2*(c3a2*p2z + (3*c2a2*m2*n2*(-(n2*p2x) + l2*p2z))*(ralpha22**2) - 
     1    (c2a2*m22*n2*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    (8*c3a2*m22*n2*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    c2a2*m2*p2x*ralpha22 + c2a2*m2*n2*p2y*ralpha22 - 
     1    2*c3a2*m2*n2*p2y*ralpha22 - 2*c3a2*m22*p2z*ralpha22 - 
     1    2*c3a2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - 
     1    (3*m2*n2*(-(n2*p2x) + l2*p2z)*s2)*(ralpha22**2) + 
     1    (5*m22*n2*(l2*p2x + m2*p2y + n2*p2z)*s2)*(ralpha22**2) - 
     1    m2*p2x*ralpha22*s2 - 2*m2*n2*p2y*ralpha22*s2 - 
     1    m22*p2z*ralpha22*s2 + m2*n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
     1 2*(c3a2*p2y + (3*c2a2*m2*n2*(m2*p2x - l2*p2y))*(ralpha22**2) - 
     1    (c2a2*m2*n22*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    (8*c3a2*m2*n22*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) - 
     1    c2a2*n2*p2x*ralpha22 - 2*c3a2*n22*p2y*ralpha22 + 
     1    c2a2*m2*n2*p2z*ralpha22 - 2*c3a2*m2*n2*p2z*ralpha22 - 
     1    2*c3a2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - 
     1    (3*m2*n2*(m2*p2x - l2*p2y)*s2)*(ralpha22**2) + 
     1    (5*m2*n22*(l2*p2x + m2*p2y + n2*p2z)*s2)*(ralpha22**2) + 
     1    n2*p2x*ralpha22*s2 - n22*p2y*ralpha22*s2 + 
     1    m2*n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 2*m2*n2*p2z*ralpha22*s2 - 
     1    m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
     1 -DB(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1 DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1    C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DDBC

      DOUBLE PRECISION FUNCTION DDCC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1                               C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)
       IMPLICIT NONE
       DOUBLE PRECISION C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22
       DOUBLE PRECISION p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist
       DOUBLE PRECISION DC

      DDCC 
     1 = (rdist*(2*(c3a2*l2*p2z - c2a2*n2*(n2*p2y - m2*p2z)*ralpha22 - 
     1     2*c3a2*l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2x*s2 - 
     1     p2y*s2 + n2*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1     l2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)**2 + 
     1 2*(c3a2*m2*p2z - c2a2*n2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1     2*c3a2*m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 + 
     1     n2*p2y*s2 + n2*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1     m2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)**2 + 
     1 2*(c3a2*n2*p2z + c3a2*(l2*p2x + m2*p2y + n2*p2z) - 
     1     c2a2*n2*(m2*p2x - l2*p2y)*ralpha22 - 
     1     2*c3a2*n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + n2*p2z*s2 + 
     1     n2*(m2*p2x - l2*p2y)*ralpha22*s2 - 
     1     n22*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)**2 + 
     1 2*((3*c2a2*n22*(n2*p2y - m2*p2z))*(ralpha22**2) - 
     1    (c2a2*l2*n22*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    (8*c3a2*l2*n22*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    c2a2*n22*p2x*ralpha22 - 2*c2a2*n2*p2y*ralpha22 - 
     1    4*c3a2*l2*n2*p2z*ralpha22 - c2a2*(n2*p2y - m2*p2z)*ralpha22 - 
     1    2*c3a2*l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2x*s2 - 
     1    (3*n22*(n2*p2y - m2*p2z)*s2)*(ralpha22**2) + 
     1    (5*l2*n22*(l2*p2x + m2*p2y + n2*p2z)*s2)*(ralpha22**2) - 
     1    n22*p2x*ralpha22*s2 + 2*n2*p2y*ralpha22*s2 - 
     1    2*l2*n2*p2z*ralpha22*s2 + (n2*p2y - m2*p2z)*ralpha22*s2 + 
     1    n22*(n2*p2y - m2*p2z)*ralpha22*s2 - 
     1    l2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c2a1*p1x - c3a1*l1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2x + 
     1    c3a2*l2*(l2*p2x + m2*p2y + n2*p2z) + (n1*p1y - m1*p1z)*s1 - 
     1    (n2*p2y - m2*p2z)*s2 + x1 - x2) + 
     1 2*((3*c2a2*n22*(-(n2*p2x) + l2*p2z))*(ralpha22**2) - 
     1    (c2a2*m2*n22*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    (8*c3a2*m2*n22*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    2*c2a2*n2*p2x*ralpha22 + c2a2*n22*p2y*ralpha22 - 
     1    4*c3a2*m2*n2*p2z*ralpha22 - c2a2*(-(n2*p2x) + l2*p2z)*ralpha22 - 
     1    2*c3a2*m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 + p2y*s2 - 
     1    (3*n22*(-(n2*p2x) + l2*p2z)*s2)*(ralpha22**2) + 
     1    (5*m2*n22*(l2*p2x + m2*p2y + n2*p2z)*s2)*(ralpha22**2) - 
     1    2*n2*p2x*ralpha22*s2 - n22*p2y*ralpha22*s2 - 
     1    2*m2*n2*p2z*ralpha22*s2 + (-(n2*p2x) + l2*p2z)*ralpha22*s2 + 
     1    n22*(-(n2*p2x) + l2*p2z)*ralpha22*s2 - 
     1    m2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c2a1*p1y - c3a1*m1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2y + 
     1    c3a2*m2*(l2*p2x + m2*p2y + n2*p2z) + (-(n1*p1x) + l1*p1z)*s1 - 
     1    (-(n2*p2x) + l2*p2z)*s2 + y1 - y2) + 
     1 2*((3*c2a2*n22*(m2*p2x - l2*p2y))*(ralpha22**2) + 2*c3a2*p2z - 
     1    (c2a2*n2**3*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) + 
     1    (8*c3a2*n2**3*(l2*p2x + m2*p2y + n2*p2z))*(ralpha22**2) - 
     1    c2a2*(m2*p2x - l2*p2y)*ralpha22 + c2a2*n22*p2z*ralpha22 - 
     1    4*c3a2*n22*p2z*ralpha22 - 
     1    6*c3a2*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22 - 
     1    (3*n22*(m2*p2x - l2*p2y)*s2)*(ralpha22**2) + p2z*s2 + 
     1    (5*n2**3*(l2*p2x + m2*p2y + n2*p2z)*s2)*(ralpha22**2) + 
     1    (m2*p2x - l2*p2y)*ralpha22*s2 + 
     1    n22*(m2*p2x - l2*p2y)*ralpha22*s2 - 3*n22*p2z*ralpha22*s2 - 
     1    3*n2*(l2*p2x + m2*p2y + n2*p2z)*ralpha22*s2)*
     1  (c2a1*p1z - c3a1*n1*(l1*p1x + m1*p1y + n1*p1z) - c2a2*p2z + 
     1    c3a2*n2*(l2*p2x + m2*p2y + n2*p2z) + (m1*p1x - l1*p1y)*s1 - 
     1    (m2*p2x - l2*p2y)*s2 + z1 - z2)))/2.
     1 -DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)*
     1  DC(p1x,p1y,p1z,p2x,p2y,p2z,x1,y1,z1,n1,l1,m1,x2,y2,z2,l2,m2,n2,rdist,
     1     C2A2,RALPHA12,RALPHA22,C3A2,S2,C2A1,C3A1,S1,N12,L12,M12,L22,M22,N22)* rdist

       return
      end FUNCTION DDCC

! ****************************************************************************************************

! trj25 annotated derivatives of individual Cartesian coordinates for use with electric fields
! See notebook entries circa 18/3/05

! ****************************************************************************************************

      function dxdpx (px, py, pz, x0, y0, z0, 
     &                cos_alpha, sin_alpha, inv_alpha, inv_alpha2, 
     &                pdotX0)

        implicit none

        DOUBLE PRECISION :: DXDPX
        DOUBLE PRECISION, INTENT(IN) :: PX, PY, PZ, X0, Y0, Z0
        DOUBLE PRECISION, INTENT(IN) :: COS_ALPHA, SIN_ALPHA, INV_ALPHA, INV_ALPHA2
        DOUBLE PRECISION, INTENT(IN) :: PDOTX0

        dxdpx = pdotX0*(1.0d0-cos_alpha)*inv_alpha2 - px*x0*sin_alpha*inv_alpha +
     &          px*inv_alpha2*(x0*(1.0d0-cos_alpha) + 
     &          pdotX0*px*inv_alpha*(sin_alpha-2.0d0*(1.0d0-cos_alpha)*inv_alpha)) +
     &          (pz*y0-py*z0)*px*inv_alpha2*(cos_alpha-sin_alpha*inv_alpha)

        return

      end function dxdpx

! ****************************************************************************************************

      function dxdpy (px, py, pz, x0, y0, z0, 
     &                cos_alpha, sin_alpha, inv_alpha, inv_alpha2, 
     &                pdotX0)

        implicit none

        DOUBLE PRECISION :: DXDPY
        DOUBLE PRECISION, INTENT(IN) :: PX, PY, PZ, X0, Y0, Z0
        DOUBLE PRECISION, INTENT(IN) :: COS_ALPHA, SIN_ALPHA, INV_ALPHA, INV_ALPHA2
        DOUBLE PRECISION, INTENT(IN) :: PDOTX0

        dxdpy = px*inv_alpha2*(y0*(1.0d0-cos_alpha) + 
     &          pdotX0*py*inv_alpha*(sin_alpha-2.0d0*(1.0d0-cos_alpha)*inv_alpha)) -
     &          (py*x0+z0)*sin_alpha*inv_alpha +
     &          (pz*y0-py*z0)*py*inv_alpha2*(cos_alpha-sin_alpha*inv_alpha)
        
        return

      end function dxdpy

! ****************************************************************************************************

      function dxdpz (px, py, pz, x0, y0, z0, 
     &                cos_alpha, sin_alpha, inv_alpha, inv_alpha2, 
     &                pdotX0)

        implicit none

        DOUBLE PRECISION :: DXDPZ
        DOUBLE PRECISION, INTENT(IN) :: PX, PY, PZ, X0, Y0, Z0
        DOUBLE PRECISION, INTENT(IN) :: COS_ALPHA, SIN_ALPHA, INV_ALPHA, INV_ALPHA2
        DOUBLE PRECISION, INTENT(IN) :: PDOTX0

        dxdpz = (y0-pz*x0)*sin_alpha*inv_alpha +
     &          px*inv_alpha2*(z0*(1.0d0-cos_alpha) + 
     &          pdotX0*pz*inv_alpha*(sin_alpha-2.0d0*(1.0d0-cos_alpha)*inv_alpha)) +
     &          (pz*y0-py*z0)*pz*inv_alpha2*(cos_alpha-sin_alpha*inv_alpha)

        return

      end function dxdpz

! ****************************************************************************************************

      function dydpx (px, py, pz, x0, y0, z0, 
     &                cos_alpha, sin_alpha, inv_alpha, inv_alpha2, 
     &                pdotX0)

        implicit none

        DOUBLE PRECISION :: DYDPX
        DOUBLE PRECISION, INTENT(IN) :: PX, PY, PZ, X0, Y0, Z0
        DOUBLE PRECISION, INTENT(IN) :: COS_ALPHA, SIN_ALPHA, INV_ALPHA, INV_ALPHA2
        DOUBLE PRECISION, INTENT(IN) :: PDOTX0

        dydpx = (z0-px*y0)*sin_alpha*inv_alpha +
     &          py*inv_alpha2*(x0*(1.0d0-cos_alpha) + 
     &          pdotX0*px*inv_alpha*(sin_alpha-2.0d0*(1.0d0-cos_alpha)*inv_alpha)) +
     &          (px*z0-pz*x0)*px*inv_alpha2*(cos_alpha-sin_alpha*inv_alpha)

        return

      end function dydpx

! ****************************************************************************************************

      function dydpy (px, py, pz, x0, y0, z0, 
     &                cos_alpha, sin_alpha, inv_alpha, inv_alpha2, 
     &                pdotX0)

        implicit none

        DOUBLE PRECISION :: DYDPY
        DOUBLE PRECISION, INTENT(IN) :: PX, PY, PZ, X0, Y0, Z0
        DOUBLE PRECISION, INTENT(IN) :: COS_ALPHA, SIN_ALPHA, INV_ALPHA, INV_ALPHA2
        DOUBLE PRECISION, INTENT(IN) :: PDOTX0

        dydpy = pdotX0*(1.0d0-cos_alpha)*inv_alpha2 - py*y0*sin_alpha*inv_alpha +
     &          py*inv_alpha2*(y0*(1.0d0-cos_alpha) + 
     &          pdotX0*py*inv_alpha*(sin_alpha-2.0d0*(1.0d0-cos_alpha)*inv_alpha)) +
     &          (px*z0-pz*x0)*py*inv_alpha2*(cos_alpha-sin_alpha*inv_alpha)

        return

      end function dydpy

! ****************************************************************************************************

      function dydpz (px, py, pz, x0, y0, z0, 
     &                cos_alpha, sin_alpha, inv_alpha, inv_alpha2, 
     &                pdotX0)

        implicit none

        DOUBLE PRECISION :: DYDPZ
        DOUBLE PRECISION, INTENT(IN) :: PX, PY, PZ, X0, Y0, Z0
        DOUBLE PRECISION, INTENT(IN) :: COS_ALPHA, SIN_ALPHA, INV_ALPHA, INV_ALPHA2
        DOUBLE PRECISION, INTENT(IN) :: PDOTX0

        dydpz = py*inv_alpha2*(z0*(1.0d0-cos_alpha) + 
     &          pdotX0*pz*inv_alpha*(sin_alpha-2.0d0*(1.0d0-cos_alpha)*inv_alpha)) -
     &          (pz*y0+x0)*sin_alpha*inv_alpha +
     &          (px*z0-pz*x0)*pz*inv_alpha2*(cos_alpha-sin_alpha*inv_alpha)

        return

      end function dydpz

! ****************************************************************************************************

! Wrapper function for dzdXrb

      function dzdXrb (coordNo,px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)

        implicit none

        DOUBLE PRECISION :: DZDXRB
        integer, intent(IN) :: coordNo
        DOUBLE PRECISION, INTENT(IN) :: PX,PY,PZ,X0,Y0,Z0,COS_ALPHA,SINA_A,INV_ALPHA2,PDOTX0
        
! External functions

        DOUBLE PRECISION :: DZDPX,DZDPY,DZDPZ
        
        select case (coordNo)
           case (1:2)
              dzdXrb = 0.0d0
           case (3)
              dzdXrb = 1.0d0
           case (4)
              dzdXrb = dzdpx(px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)
           case (5)
              dzdXrb = dzdpy(px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)              
           case (6)
              dzdXrb = dzdpz(px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)
           case default
              print *, 'Out of range coordinate index passed to dzdXrb'
              stop
        end select

        return

      end function dzdXrb


! ****************************************************************************************************

      function dzdpx (px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)

        implicit none

        DOUBLE PRECISION :: DZDPX
        DOUBLE PRECISION, INTENT(IN) :: PX, PY, PZ, X0, Y0, Z0
        DOUBLE PRECISION, INTENT(IN) :: COS_ALPHA, INV_ALPHA2
        DOUBLE PRECISION, INTENT(IN) :: SINA_A    ! SIN(ALPHA)/ALPHA
        DOUBLE PRECISION, INTENT(IN) :: PDOTX0

        dzdpx = pz*inv_alpha2*(x0*(1.0d0-cos_alpha) + 
     &          pdotX0*px*(sina_a-2.0d0*(1.0d0-cos_alpha)*inv_alpha2)) -
     &          sina_a*(px*z0+y0) + 
     &          (py*x0-px*y0)*px*inv_alpha2*(cos_alpha-sina_a)
        
        return
 
      end function dzdpx

! ****************************************************************************************************

      function dzdpy (px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)

        implicit none

        DOUBLE PRECISION :: DZDPY
        DOUBLE PRECISION, INTENT(IN) :: PX, PY, PZ, X0, Y0, Z0
        DOUBLE PRECISION, INTENT(IN) :: COS_ALPHA, INV_ALPHA2
        DOUBLE PRECISION, INTENT(IN) :: SINA_A    ! SIN(ALPHA)/ALPHA
        DOUBLE PRECISION, INTENT(IN) :: PDOTX0

        dzdpy = sina_a*(x0-z0*py) + pz*inv_alpha2*(y0*(1.0d0-cos_alpha) + 
     &          pdotX0*py*(sina_a-2.0d0*(1.0d0-cos_alpha)*inv_alpha2)) +
     &          (py*x0-px*y0)*py*inv_alpha2*(cos_alpha-sina_a)
        
        return

      end function dzdpy

! ****************************************************************************************************

      function dzdpz (px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)

        implicit none

        DOUBLE PRECISION :: DZDPZ
        DOUBLE PRECISION, INTENT(IN) :: PX, PY, PZ, X0, Y0, Z0
        DOUBLE PRECISION, INTENT(IN) :: COS_ALPHA, INV_ALPHA2
        DOUBLE PRECISION, INTENT(IN) :: SINA_A    ! SIN(ALPHA)/ALPHA
        DOUBLE PRECISION, INTENT(IN) :: PDOTX0

        dzdpz = pdotX0*(1.0d0-cos_alpha)*inv_alpha2 - pz*z0*sina_a +
     &          pz*inv_alpha2*(z0*(1.0d0-cos_alpha) + 
     &          pdotX0*pz*(sina_a-2.0d0*(1.0d0-cos_alpha)*inv_alpha2)) +
     &          (py*x0-px*y0)*pz*inv_alpha2*(cos_alpha-sina_a)

        return

      end function dzdpz

! ****************************************************************************************************

! And selected second derivatives

! ****************************************************************************************************

! Wrapper function for the rigid body second derivatives
! coordNo`s run from 1 to 6, with 1 specifying dz/dxcom, 2 specifying dz/dycom etc

      function d2zdXrb1dXrb2(coordNo1,coordNo2,
     1                       px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)

       implicit none

       DOUBLE PRECISION :: D2ZDXRB1DXRB2
       integer :: coordNo1,coordNo2
       DOUBLE PRECISION :: PX,PY,PZ,X0,Y0,Z0,COS_ALPHA,SINA_A,INV_ALPHA2,PDOTX0

! Local variables

       integer :: firstCoord, secondCoord

! External functions

       DOUBLE PRECISION :: D2ZDPX2,D2ZDPXDPY,D2ZDPXDPZ,D2ZDPY2,D2ZDPYDPZ,D2ZDPZ2

! Check for range

       if (coordNo1.lt.1 .or. coordNo1.gt.6) then
          print *, 'Out of range coordNo1 passed to d2zdXrb1dXrb2:', coordNo1
          stop
       endif
       if (coordNo2.lt.1 .or. coordNo2.gt.6) then
          print *, 'Out of range coordNo1 passed to d2zdXrb1dXrb2:', coordNo2
          stop
       endif

       if (coordNo1.gt.coordNo2) then
          firstCoord = coordNo2
          secondCoord = coordNo1
       else
          firstCoord = coordNo1
          secondCoord = coordNo2          
       end if

       select case (firstCoord)
       case (4)
          select case (secondCoord)
          case (4)   
             d2zdXrb1dXrb2 = d2zdpx2(px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)
          case (5)
             d2zdXrb1dXrb2 = d2zdpxdpy(px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)
          case (6)
             d2zdXrb1dXrb2 = d2zdpxdpz(px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)
          case default
             print *, 'Should not get here'
             stop
          end select

       case (5)
          select case (secondCoord)
          case (5)
             d2zdXrb1dXrb2 = d2zdpy2(px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)
          case (6)
             d2zdXrb1dXrb2 = d2zdpydpz(px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)
          case default
             print *, 'Should not get here'
             stop
          end select

       case (6)
          select case (secondCoord)
          case (6)
             d2zdXrb1dXrb2 = d2zdpz2(px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)
          case default
             print *, 'Should not get here'
             stop
          end select

       case default
          d2zdXrb1dXrb2 = 0.0d0

       end select

       return

      end function d2zdXrb1dXrb2

! ****************************************************************************************************

      function d2zdpx2 (px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)

        implicit none

        DOUBLE PRECISION :: D2ZDPX2
        DOUBLE PRECISION, INTENT(IN) :: PX, PY, PZ, X0, Y0, Z0
        DOUBLE PRECISION, INTENT(IN) :: COS_ALPHA, INV_ALPHA2
        DOUBLE PRECISION, INTENT(IN) :: SINA_A    ! SIN(ALPHA)/ALPHA
        DOUBLE PRECISION, INTENT(IN) :: PDOTX0

        d2zdpx2 = -2.0d0*px*(y0+px*(py*x0-px*y0)*inv_alpha2)*cos_alpha*inv_alpha2 +
     &            (8.0d0*(px**2)*inv_alpha2-2.0d0)*pz*(1.0d0-cos_alpha)*pdotX0*(inv_alpha2**2) +
     &            (2.0d0*px*y0 + (3.0d0*(px**2)*inv_alpha2-1.0d0)*(py*x0-px*y0))*sina_a*inv_alpha2 +
     &            (py*x0-px*y0)*(cos_alpha*(1.0d0-(px**2)*inv_alpha2)-(px**2)*sina_a)*inv_alpha2 +
     &            z0*((px**2)*inv_alpha2*(sina_a-cos_alpha) - sina_a) -
     &            4.0d0*px*pz*(x0*(1.0d0-cos_alpha) + px*pdotX0*sina_a)*(inv_alpha2**2) +
     &            pz*(2.0d0*px*x0*sina_a+pdotX0*((px**2)*inv_alpha2*(cos_alpha-sina_a)+sina_a))*inv_alpha2
        
        return
 
      end function d2zdpx2

! ****************************************************************************************************

      function d2zdpxdpy (px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)

        implicit none

        DOUBLE PRECISION :: D2ZDPXDPY
        DOUBLE PRECISION, INTENT(IN) :: PX, PY, PZ, X0, Y0, Z0
        DOUBLE PRECISION, INTENT(IN) :: COS_ALPHA, INV_ALPHA2
        DOUBLE PRECISION, INTENT(IN) :: SINA_A    ! SIN(ALPHA)/ALPHA
        DOUBLE PRECISION, INTENT(IN) :: PDOTX0

        d2zdpxdpy = inv_alpha2*(cos_alpha*(px*x0 - py*y0 - px*py*z0) + 
     &                          sina_a*(py*pz*x0+py*y0+px*pz*y0+px*py*z0-px*py*(py*x0-px*y0)-px*x0)) +
     &              (inv_alpha2**2)*(cos_alpha*px*py*(pz*pdotX0 - 3.0d0*(py*x0 - px*y0)) +
     &                               sina_a*px*py*(3.0d0*(py*x0 - px*y0) - 5.0d0*pz*pdotX0) -
     &                               2.0d0*(1.0d0 - cos_alpha)*pz*(py*x0 + px*y0)) +
     &              (inv_alpha2**3)*8.0d0*px*py*pz*(1.0d0 - cos_alpha)*pdotX0 
 
        return
 
      end function d2zdpxdpy

! ****************************************************************************************************

      function d2zdpxdpz (px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)

        implicit none

        DOUBLE PRECISION :: D2ZDPXDPZ
        DOUBLE PRECISION, INTENT(IN) :: PX, PY, PZ, X0, Y0, Z0
        DOUBLE PRECISION, INTENT(IN) :: COS_ALPHA, INV_ALPHA2
        DOUBLE PRECISION, INTENT(IN) :: SINA_A    ! SIN(ALPHA)/ALPHA
        DOUBLE PRECISION, INTENT(IN) :: PDOTX0

        d2zdpxdpz = inv_alpha2*(x0*(1.0d0-cos_alpha) - cos_alpha*pz*(y0+px*z0) + 
     &                          sina_a*(pz*(pz*x0+y0-px*(py*x0-px*y0)+2.0d0*px*z0)+px*pdotX0)) +
     &              (inv_alpha2**2)*(cos_alpha*px*pz*(pz*pdotX0-3.0d0*(py*x0-px*y0)) + 
     &                               sina_a*px*pz*(3.0d0*(py*x0-px*y0)-5.0d0*pz*pdotX0) -
     &                               2.0d0*(1.0d0-cos_alpha)*((pz**2)*x0+px*pz*z0+px*pdotX0)) +
     &              (inv_alpha2**3)*8.0d0*px*(pz**2)*(1.0d0 - cos_alpha)*pdotX0 

        return
 
      end function d2zdpxdpz

! ****************************************************************************************************

      function d2zdpy2 (px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)

        implicit none

        DOUBLE PRECISION :: D2ZDPY2
        DOUBLE PRECISION, INTENT(IN) :: PX, PY, PZ, X0, Y0, Z0
        DOUBLE PRECISION, INTENT(IN) :: COS_ALPHA, INV_ALPHA2
        DOUBLE PRECISION, INTENT(IN) :: SINA_A    ! SIN(ALPHA)/ALPHA
        DOUBLE PRECISION, INTENT(IN) :: PDOTX0

        d2zdpy2 = 2.0d0*py*(x0-py*(py*x0-px*y0)*inv_alpha2)*cos_alpha*inv_alpha2 +
     &            (8.0d0*(py**2)*inv_alpha2-2.0d0)*pz*(1.0d0-cos_alpha)*pdotX0*(inv_alpha2**2) +
     &            (-2.0d0*py*x0 + (3.0d0*(py**2)*inv_alpha2-1.0d0)*(py*x0-px*y0))*sina_a*inv_alpha2 +
     &            (py*x0-px*y0)*(cos_alpha*(1.0d0-(py**2)*inv_alpha2)-(py**2)*sina_a)*inv_alpha2 +
     &            z0*((py**2)*inv_alpha2*(sina_a-cos_alpha) - sina_a) -
     &            4.0d0*py*pz*(y0*(1.0d0-cos_alpha) + py*pdotX0*sina_a)*(inv_alpha2**2) +
     &            pz*(2.0d0*py*y0*sina_a+pdotX0*((py**2)*inv_alpha2*(cos_alpha-sina_a)+sina_a))*inv_alpha2
        
        return
 
      end function d2zdpy2

! ****************************************************************************************************

      function d2zdpydpz (px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)

        implicit none

        DOUBLE PRECISION :: D2ZDPYDPZ
        DOUBLE PRECISION, INTENT(IN) :: PX, PY, PZ, X0, Y0, Z0
        DOUBLE PRECISION, INTENT(IN) :: COS_ALPHA, INV_ALPHA2
        DOUBLE PRECISION, INTENT(IN) :: SINA_A    ! SIN(ALPHA)/ALPHA
        DOUBLE PRECISION, INTENT(IN) :: PDOTX0

        d2zdpydpz = inv_alpha2*(y0*(1.0d0-cos_alpha) + cos_alpha*pz*(x0-py*z0) +
     &                          sina_a*(pz*(pz*y0+2.0d0*py*z0-x0-py*(py*x0 - px*y0)) + py*pdotX0)) +
     &              (inv_alpha2**2)*(cos_alpha*py*pz*(pz*pdotX0-3.0d0*(py*x0-px*y0)) +
     &                               sina_a*py*pz*(3.0d0*(py*x0 - px*y0) - 5.0d0*pz*pdotX0) -
     &                               2.0d0*(1.0d0 - cos_alpha)*(pz*(pz*y0+py*z0) + py*pdotX0)) +
     &              (inv_alpha2**3)*8.0d0*py*(pz**2)*(1.0d0 - cos_alpha)*pdotX0 

        return
 
      end function d2zdpydpz 

! ****************************************************************************************************

      function d2zdpz2 (px,py,pz,x0,y0,z0,cos_alpha,sina_a,inv_alpha2,pdotX0)

        implicit none

        DOUBLE PRECISION :: D2ZDPZ2
        DOUBLE PRECISION, INTENT(IN) :: PX, PY, PZ, X0, Y0, Z0
        DOUBLE PRECISION, INTENT(IN) :: COS_ALPHA, INV_ALPHA2
        DOUBLE PRECISION, INTENT(IN) :: SINA_A    ! SIN(ALPHA)/ALPHA
        DOUBLE PRECISION, INTENT(IN) :: PDOTX0

        d2zdpz2 = -2.0d0*(pz**2)*(py*x0-px*y0)*cos_alpha*(inv_alpha2**2) +
     &            pz*(8.0d0*(pz**2)*inv_alpha2 - 6.0d0)*(1.0d0-cos_alpha)*pdotX0*(inv_alpha2**2) +
     &            (3.0d0*(pz**2)*inv_alpha2 - 1.0d0)*(py*x0 - px*y0)*sina_a*inv_alpha2 +
     &            (py*x0-px*y0)*(cos_alpha-(pz**2)*(cos_alpha*inv_alpha2-sina_a))*inv_alpha2 +
     &            z0*((pz**2)*inv_alpha2*(sina_a-cos_alpha) - sina_a) +
     &            2.0d0*inv_alpha2*(1.0d0-2.0d0*(pz**2)*inv_alpha2)*(z0*(1.0d0-cos_alpha)+pz*pdotX0*sina_a) +
     &            pz*(2.0d0*pz*z0*sina_a+pdotX0*((pz**2)*inv_alpha2*(cos_alpha-sina_a)+sina_a))*inv_alpha2
        
        return
 
      end function d2zdpz2       

! ****************************************************************************************************






