!   NEB module is an implementation of the nudged elastic band method for performing double-ended pathway searches.
!   Copyright (C) 2003-2006 Semen A. Trygubenko and David J. Wales
!   This file is part of NEB module. NEB module is part of OPTIM.
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
SUBROUTINE OLDNEBGRADIENT
USE KEYGRAD
USE NEBDATA
USE KEYNEB,ONLY: NIMAGE
USE TANGENT
USE NEBUTILS
USE KEY, ONLY: UNRST, FROZEN, FREEZE, NEBK, BULKT, TWOD, VARIABLES
USE COMMONS, ONLY : PARAM1, PARAM2, PARAM3
USE GRADIENTS
IMPLICIT NONE
INTEGER J1, J2
DOUBLE PRECISION DUMMY1, DUMMY2

!
! Image 1 is START
! Movable images are numbers 2 to NIMAGE+1
! Image NIMAGE+2 is FINISH
! The following assignments are made at the start of NEWNEB
!
! X         => XYZ(NOPT+1:NOPT*(NIMAGE+1))
! EIMAGE    => EEE(2:NIMAGE+1)
! G         => GGG(NOPT+1:NOPT*(NIMAGE+1))
! GSPR      => SSS(NOPT+1:NOPT*(NIMAGE+1))
! RMSFIMAGE => RRR(2:NIMAGE+1)
! TANPTR => TANVEC
! EEE = 0.0D0
! EEE(1)=EINITIAL
! EEE(NIMAGE+2)=EFINAL
! 
! TRUEPOTEG calculates energy and gradients for images
! coordinates are in XYZ(1:NOPT*(NIMAGE+2))
! gradient is in GGG(1:NOPT*(NIMAGE+2)) and also saved in TRUEGRAD(1:NOPT*(NIMAGE+2))
! image potential energies are in EEE(1:NIMAGE+2)
!
CALL TRUEPOTEG
!
! MEPTANGENT calculates tangent vector according to Henkelmann and Jonsson, JCP, 113, 9978, 2000
! and stores it in TANVEC(1:NOPT,1:NIMAGE)
! The tangent vector for image J1 is stored in TANVEC(:,J1-1)
!
CALL MEPTANGENT        

!
!  Gradient of the potential perpendicular to the tangent vector.
!
ETOTAL=0.0D0
DO J1=2,NIMAGE+1
   DUMMY1=0.0D0
   DO J2=1,NOPT
      DUMMY1=DUMMY1+GGG(J2+NOPT*(J1-1))*TANVEC(J2,J1-1) 
   ENDDO
   DO J2=1,NOPT
      GGG(J2+NOPT*(J1-1))=GGG(J2+NOPT*(J1-1))-DUMMY1*TANVEC(J2,J1-1)
      RMS=RMS+GGG(J2+NOPT*(J1-1))**2
   ENDDO
   ETOTAL=ETOTAL+EEE(J1)
!  PRINT '(A,I6)','oldnebgrad> GGG perp for image ',J1
!  PRINT '(6G20.10)',GGG(1+NOPT*(J1-1):J1*NOPT)
!  PRINT '(A,I6)','oldnebgrad> tanvec for image ',J1
!  PRINT '(6G20.10)',TANVEC(1:NOPT,J1-1)
ENDDO
RMS=SQRT(RMS/(NIMAGE*NOPT*1.0D0))
!
!  Spring energy derivatives.
!
DO J1=2,NIMAGE+1
   DUMMY1=0.0D0
   DUMMY2=0.0D0
   IF (BULKT) THEN
      DO J2=1,NATOMS
         DUMMY1=DUMMY1+ &
  &                 (XYZ(NOPT*(J1-1)+3*(J2-1)+1)-XYZ(NOPT*J1+3*(J2-1)+1) &
  &    -PARAM1*NINT((XYZ(NOPT*(J1-1)+3*(J2-1)+1)-XYZ(NOPT*J1+3*(J2-1)+1)/PARAM1)))**2 &
  &                +(XYZ(NOPT*(J1-1)+3*(J2-1)+2)-XYZ(NOPT*J1+3*(J2-1)+2) &
  &    -PARAM2*NINT((XYZ(NOPT*(J1-1)+3*(J2-1)+2)-XYZ(NOPT*J1+3*(J2-1)+2)/PARAM2)))**2 
         IF (.NOT.TWOD) DUMMY1=DUMMY1+ &
  &                 (XYZ(NOPT*(J1-1)+3*(J2-1)+3)-XYZ(NOPT*J1+3*(J2-1)+3) &
  &    -PARAM3*NINT((XYZ(NOPT*(J1-1)+3*(J2-1)+3)-XYZ(NOPT*J1+3*(J2-1)+3)/PARAM3)))**2

         DUMMY2=DUMMY2+ &
  &                 (XYZ(NOPT*(J1-1)+3*(J2-1)+1)-XYZ(NOPT*(J1-2)+3*(J2-1)+1) &
  &    -PARAM1*NINT((XYZ(NOPT*(J1-1)+3*(J2-1)+1)-XYZ(NOPT*(J1-2)+3*(J2-1)+1)/PARAM1)))**2 &
  &                +(XYZ(NOPT*(J1-1)+3*(J2-1)+2)-XYZ(NOPT*(J1-2)+3*(J2-1)+2) &
  &    -PARAM2*NINT((XYZ(NOPT*(J1-1)+3*(J2-1)+2)-XYZ(NOPT*(J1-2)+3*(J2-1)+2)/PARAM2)))**2 
         IF (.NOT.TWOD) DUMMY1=DUMMY1+ &
  &                 (XYZ(NOPT*(J1-1)+3*(J2-1)+3)-XYZ(NOPT*(J1-2)+3*(J2-1)+3) &
  &    -PARAM3*NINT((XYZ(NOPT*(J1-1)+3*(J2-1)+3)-XYZ(NOPT*(J1-2)+3*(J2-1)+3)/PARAM3)))**2
      ENDDO
   ELSE
      DO J2=1,NOPT
         DUMMY1=DUMMY1+(XYZ(J2+NOPT*(J1-1))-XYZ(J2+NOPT*J1))**2
         DUMMY2=DUMMY2+(XYZ(J2+NOPT*(J1-1))-XYZ(J2+NOPT*(J1-2)))**2
      ENDDO
   ENDIF
   DUMMY1=SQRT(DUMMY1)*NEWNEBK(J1-1)
   DUMMY2=SQRT(DUMMY2)*NEWNEBK(J1)
   DO J2=1,NOPT
      GGG(J2+NOPT*(J1-1))=GGG(J2+NOPT*(J1-1))-(DUMMY1-DUMMY2)*TANVEC(J2,J1-1)
   ENDDO
!  PRINT '(A,I6)','oldnebgrad> GGG g perp + with parallel spring derivatives for image ',J1
!  PRINT '(6G20.10)',GGG(1+NOPT*(J1-1):J1*NOPT)
ENDDO
!   
! Set gradients on frozen atoms to zero.
!   
IF (FREEZE) THEN
   DO J1=1,NIMAGE+2
      DO J2=1,NATOMS
         IF (FROZEN(J2)) THEN
            IF (VARIABLES) THEN
               GGG(NOPT*(J1-1)+J2)=0.0D0
            ELSE
               GGG(NOPT*(J1-1)+3*(J2-1)+1)=0.0D0
               GGG(NOPT*(J1-1)+3*(J2-1)+2)=0.0D0
               GGG(NOPT*(J1-1)+3*(J2-1)+3)=0.0D0
            ENDIF
         ENDIF
      ENDDO
   ENDDO
ENDIF

IF (ORT) THEN ! REMOVE OVERALL ROTATION/TRANSLATION BY FREEZING 6 DEGREES OF FREEDOM
   G(1::NOPT)=0.0D0
   G(2::NOPT)=0.0D0
   G(3::NOPT)=0.0D0
   G(4::NOPT)=0.0D0
   G(5::NOPT)=0.0D0
   G(7::NOPT)=0.0D0
!  DO J1=2,NIMAGE+1  
!     GGG(NOPT*(J1-1)+1)=0.0D0
!     GGG(NOPT*(J1-1)+2)=0.0D0
!     GGG(NOPT*(J1-1)+3)=0.0D0
!     GGG(NOPT*(J1-1)+4)=0.0D0
!     GGG(NOPT*(J1-1)+5)=0.0D0
!     GGG(NOPT*(J1-1)+7)=0.0D0
!  ENDDO
ENDIF

RMS=0.0D0
DO J1=NOPT+1,(NIMAGE+1)*NOPT
   RMS=RMS+GGG(J1)**2
ENDDO
RMS=SQRT(RMS/(NIMAGE*NOPT))

END SUBROUTINE OLDNEBGRADIENT
