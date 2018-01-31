!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999-2009 David J. Wales
!   This file is part of PATHSAMPLE.
!
!   PATHSAMPLE is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   PATHSAMPLE is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!

!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! Reads a database of minima and calculates the backbone dihedral angle omega for
! PHE(167)CA-PHE(167)C-Pro(168)N-Pro(168)CA.  
! 
      SUBROUTINE CALCORDER(NATOMS,NMIN,NTS,UMIN,UTS,DEBUG)
      IMPLICIT NONE
      INTEGER J1,J2,NMIN,NATOMS,UMIN,UTS,NTS
      DOUBLE PRECISION LOCALPOINTS(3*NATOMS),MAG,F(3),G(3),H(3),A(3),B(3),C(3),CST,PHI 
      LOGICAL DEBUG

      WRITE(*,*) 'calcorder> Calculating order paramters for CypA' 
      OPEN(UNIT=91,FILE='cis.dbase',STATUS='UNKNOWN')
      OPEN(UNIT=92,FILE='trans.dbase',STATUS='UNKNOWN')     
      DO J1=1,NMIN
        READ(UMIN,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
      	F(1)=LOCALPOINTS(3*(1556-1)+1)-LOCALPOINTS(3*(1564-1)+1)
       	F(2)=LOCALPOINTS(3*(1556-1)+2)-LOCALPOINTS(3*(1564-1)+2)
        F(3)=LOCALPOINTS(3*(1556-1)+3)-LOCALPOINTS(3*(1564-1)+3)
        MAG=F(1)*F(1)+F(2)*F(2)+F(3)*F(3)
        MAG=SQRT(MAG)
        F(1)=F(1)/MAG
        F(2)=F(2)/MAG
        F(3)=F(3)/MAG
        G(1)=LOCALPOINTS(3*(1564-1)+1)-LOCALPOINTS(3*(1566-1)+1)
        G(2)=LOCALPOINTS(3*(1564-1)+2)-LOCALPOINTS(3*(1566-1)+2)
        G(3)=LOCALPOINTS(3*(1564-1)+3)-LOCALPOINTS(3*(1566-1)+3)
        MAG=G(1)*G(1)+G(2)*G(2)+G(3)*G(3)
        MAG=SQRT(MAG)
        G(1)=G(1)/MAG
        G(2)=G(2)/MAG
        G(3)=G(3)/MAG
        H(1)=LOCALPOINTS(3*(1568-1)+1)-LOCALPOINTS(3*(1566-1)+1)
        H(2)=LOCALPOINTS(3*(1568-1)+2)-LOCALPOINTS(3*(1566-1)+2)
        H(3)=LOCALPOINTS(3*(1568-1)+3)-LOCALPOINTS(3*(1566-1)+3)
        MAG=H(1)*H(1)+H(2)*H(2)+H(3)*H(3)
        MAG=SQRT(MAG)
        H(1)=H(1)/MAG
        H(2)=H(2)/MAG
        H(3)=H(3)/MAG
        A(1)=(F(2)*G(3)-F(3)*G(2))
        A(2)=(F(3)*G(1)-F(1)*G(3))
        A(3)=(F(1)*G(2)-F(2)*G(1))
        MAG=A(1)*A(1)+A(2)*A(2)+A(3)*A(3)
        MAG=SQRT(MAG)
        A(1)=A(1)/MAG
        A(2)=A(2)/MAG
        A(3)=A(3)/MAG
        B(1)=H(2)*G(3)-H(3)*G(2)
        B(2)=H(3)*G(1)-H(1)*G(3)
        B(3)=H(1)*G(2)-H(2)*G(1)
        MAG=B(1)*B(1)+B(2)*B(2)+B(3)*B(3)
        MAG=SQRT(MAG)
        B(1)=B(1)/MAG
        B(2)=B(2)/MAG
        B(3)=B(3)/MAG
        CST=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)
        IF (ABS(CST).GE.1.0D0) CST=CST/ABS(CST)
        PHI=ACOS(CST)
        C(1)=A(2)*B(3)-A(3)*B(2)
        C(2)=A(3)*B(1)-A(1)*B(3)
        C(3)=A(1)*B(2)-A(2)*B(1)
        IF(G(1)*C(1)+G(2)*C(2)+G(3)*C(3).GT.0.0D0) PHI=-PHI
        PHI=(PHI*180.0D0)/3.14159265
        IF(PHI.GE.-25.0D0.AND.PHI.LE.25.0D0) THEN
           WRITE(91,*) J1,PHI
        ELSEIF(PHI.GE.150.0D0.OR.PHI.LE.-150.0D0) THEN
           WRITE(92,*) J1,PHI
        ENDIF	
      ENDDO
      CLOSE(1)
      CLOSE(2)

      RETURN 
      END
