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
SUBROUTINE INTGRADLJ(ETOTAL,XYZ,GGG,IMGFREEZE,RMS,PEDGE)
USE KEY, ONLY: FROZEN, FREEZE, FREEZENODEST, INTIMAGE, INTLJDEL, INTLJEPS, KINT, IMSEPMAX
USE COMMONS, ONLY: NATOMS, NOPT, DEBUG
IMPLICIT NONE
           
INTEGER :: J1,J2,NI2,NI1,NJ2,NJ1,NMAXINT,NMININT,NLJINT(INTIMAGE+2),J3,NSHIFT
DOUBLE PRECISION :: ETOTAL, RMS, R6, COS2, SIN2, DCUT
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,D2,D1,DI,D2SQ,G2(3),SPGRAD(3),DPLUS,ESPRING
DOUBLE PRECISION G1(3),G1INT(3),G2INT(3),G1MIN(3)
DOUBLE PRECISION DUMMY, GRAD(3), D1SQ, DISQ
DOUBLE PRECISION MAXINT,MININT,LJEINT(INTIMAGE+2),RMSIMAGE(INTIMAGE+2)
LOGICAL NOINT, LPRINT ! , EDGEINT(INTIMAGE+1,NATOMS,NATOMS)
DOUBLE PRECISION XYZ(NOPT*(INTIMAGE+2)), GGG(NOPT*(INTIMAGE+2)), XLOCAL(6*NATOMS), GLOCAL(6*NATOMS)
DOUBLE PRECISION VEC1(3), VEC2(3), VEC1M2(3), E1, E2, EINT, VDUM(3)
LOGICAL IMGFREEZE(INTIMAGE), PEDGE
DOUBLE PRECISION R1APR2BMR2AMR1BSQ, R1AMR1BSQ, R2AMR2BSQ, R1AMR1BDR2AMR2B, R1AMR1BDR2AMR2BSQ

NLJINT(1:INTIMAGE+2)=0     ! number of energy contributions from internal minima
! LJEINT(1:INTIMAGE+2)=0.0D0 ! energy contributions from internal minima
GGG(1:NOPT*(INTIMAGE+2))=0.0D0
ETOTAL=0.0D0 ! normal LJ energy summed over images
LPRINT=.FALSE.
DCUT=0.1D0
!
! Energy and forces.
!
! For J1 we consider the line segment between image J1-1 and J1.
! There are INTIMAGE+1 line segments in total, with an energy contribution
! and corresponding gradient terms for each. 
! A and B refer to atoms, 1 and 2 to images J1-1 and J1 (or J1-2 and J1-1 below).
!
! IMGFREEZE(1:INTIMAGE) refers to the images excluding end points!
!
DO J1=2,INTIMAGE+2
   GLOCAL(1:2*NOPT)=0.0D0
   IF (FREEZENODEST) THEN ! IMGFREEZE is not allocated otherwise!
      IF (J1.EQ.2) THEN
         IF (IMGFREEZE(1)) CYCLE
      ELSE IF (J1.EQ.INTIMAGE+2) THEN
         IF (IMGFREEZE(INTIMAGE)) CYCLE
      ELSE
         IF (IMGFREEZE(J1-2).AND.IMGFREEZE(J1-1)) CYCLE
      ENDIF
   ENDIF
   NSHIFT=NOPT*(J1-2)
   XLOCAL(1:2*NOPT)=XYZ(NSHIFT+1:NSHIFT+2*NOPT)
   DO J2=1,NATOMS
      NI1=3*(J2-1)
      NI2=NOPT+NI1
      R1AX=XLOCAL(NI1+1); R1AY=XLOCAL(NI1+2); R1AZ=XLOCAL(NI1+3)
      R2AX=XLOCAL(NI2+1); R2AY=XLOCAL(NI2+2); R2AZ=XLOCAL(NI2+3)
      DO J3=J2+1,NATOMS
         NJ1=3*(J3-1)
         NJ2=NOPT+NJ1
         R1BX=XLOCAL(NJ1+1); R1BY=XLOCAL(NJ1+2); R1BZ=XLOCAL(NJ1+3)
         R2BX=XLOCAL(NJ2+1); R2BY=XLOCAL(NJ2+2); R2BZ=XLOCAL(NJ2+3)
         VEC1(1)=R1AX-R1BX; VEC1(2)=R1AY-R1BY; VEC1(3)=R1AZ-R1BZ
         VEC2(1)=R2AX-R2BX; VEC2(2)=R2AY-R2BY; VEC2(3)=R2AZ-R2BZ
         VEC1M2(1:3)=VEC1(1:3)-VEC2(1:3)
         R1APR2BMR2AMR1BSQ=VEC1M2(1)**2+VEC1M2(2)**2+VEC1M2(3)**2
         NOINT=.TRUE.
         D1SQ=R1AX**2+R1AY**2+R1AZ**2+R1BX**2+R1BY**2+R1BZ**2-2*(R1AX*R1BX+R1AY*R1BY+R1AZ*R1BZ)
         G1(1:3)=VEC1(1:3)/D1SQ ! this is actually derivative of D1 divided by an extra D1
         IF (J1.GT.2) THEN ! images 1 and 2 correspond to J1-1 and J1
            R6=1.0D0/D1SQ**3
            GRAD(1:3)=-24.0D0*R6*(2.0D0*R6-1.0D0)*G1(1:3)
            ETOTAL=ETOTAL+4*R6*(R6-1.0D0)
            GLOCAL(NI1+1:NI1+3)=GLOCAL(NI1+1:NI1+3)+GRAD(1:3)
            GLOCAL(NJ1+1:NJ1+3)=GLOCAL(NJ1+1:NJ1+3)-GRAD(1:3)
         ENDIF
         IF (R1APR2BMR2AMR1BSQ.NE.0.0D0) THEN
            COS2=(VEC1(1)*VEC1M2(1)+VEC1(2)*VEC1M2(2)+VEC1(3)*VEC1M2(3))/R1APR2BMR2AMR1BSQ
            IF ((COS2.GT.0.0D0).AND.(COS2.LT.1.0D0)) THEN ! internal minimum
               NOINT=.FALSE.
               R1AMR1BDR2AMR2B=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
               R1AMR1BDR2AMR2BSQ=R1AMR1BDR2AMR2B**2
               R1AMR1BSQ=VEC1(1)**2+VEC1(2)**2+VEC1(3)**2
               R2AMR2BSQ=VEC2(1)**2+VEC2(2)**2+VEC2(3)**2
               DISQ=MAX((-R1AMR1BDR2AMR2BSQ + R1AMR1BSQ*R2AMR2BSQ)/R1APR2BMR2AMR1BSQ,0.0D0)
               DI=SQRT(DISQ); D1=SQRT(D1SQ)
               IF (D1-DI.LE.INTLJDEL) THEN
                  NOINT=.TRUE.
               ELSE
                  D2SQ=R2AX**2+R2AY**2+R2AZ**2+R2BX**2+R2BY**2+R2BZ**2-2*(R2AX*R2BX+R2AY*R2BY+R2AZ*R2BZ)
                  D2=SQRT(D2SQ)
                  IF (D2-DI.LE.INTLJDEL) THEN
                     NOINT=.TRUE.
                  ELSE
                     DUMMY=DISQ*R1APR2BMR2AMR1BSQ**2
                     G2(1:3)=VEC2(1:3)/D2SQ
                     VDUM(1:3)=(R1AMR1BDR2AMR2BSQ - R1AMR1BSQ*R2AMR2BSQ)*VEC1M2(1:3)
                     G1INT(1:3)=( VDUM(1:3)+R1APR2BMR2AMR1BSQ*(VEC1(1:3)*R2AMR2BSQ &
  &                                 -R1AMR1BDR2AMR2B*VEC2(1:3)))/DUMMY
                     G2INT(1:3)=(-VDUM(1:3)+R1APR2BMR2AMR1BSQ*(VEC2(1:3)*R1AMR1BSQ &
  &                                 -R1AMR1BDR2AMR2B*VEC1(1:3)))/DUMMY
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
!
! terms for image J1
! D1 is the distance for image 1 (R1..) J1-1
! D2 is the distance for image 2 (R2..) J1
! There are two terms for each image 2 <= J1 <= INTIMAGE+1
!
         IF (.NOT.NOINT) THEN
            DUMMY=INTLJEPS*((-D1+DI+INTLJDEL)*(-D2+DI+INTLJDEL)/DISQ**2)**2
            GRAD(1:3)=DUMMY*(-4*G1INT(1:3)+(2*(-D1*G1(1:3)+DI*G1INT(1:3)))/(-D1+DI+INTLJDEL) &
     &               +(2*DI*G1INT(1:3))/(-D2+DI+INTLJDEL))

            GLOCAL(NI1+1:NI1+3)=GLOCAL(NI1+1:NI1+3)+GRAD(1:3)
            GLOCAL(NJ1+1:NJ1+3)=GLOCAL(NJ1+1:NJ1+3)-GRAD(1:3)

            GRAD(1:3)=DUMMY*(-4*G2INT(1:3)+(2*DI*G2INT(1:3))/(-D1+DI+INTLJDEL)+(2*(-D2*G2(1:3) &
     &               +DI*G2INT(1:3)))/(-D2+DI+INTLJDEL))

            GLOCAL(NI2+1:NI2+3)=GLOCAL(NI2+1:NI2+3)+GRAD(1:3)
            GLOCAL(NJ2+1:NJ2+3)=GLOCAL(NJ2+1:NJ2+3)-GRAD(1:3)

            ETOTAL=ETOTAL+DUMMY

!           NLJINT(J1-1)=NLJINT(J1-1)+1
!           NLJINT(J1)=NLJINT(J1)+1
!           IF (PEDGE.AND.(.NOT.EDGEINT(J1-1,J2,J3))) THEN
!              PRINT '(A,I6,A,2I6,A)','edge ',J1-1,' atoms ',J2,J3,' changed from false to true'
!              PRINT '(A,4G20.10)','D1,D2,DI,c2=',D1,D2,DI,COS2
!           ENDIF
!           EDGEINT(J1-1,J2,J3)=.TRUE.
!        ELSE
!           IF (PEDGE.AND.EDGEINT(J1-1,J2,J3)) THEN
!              PRINT '(A,I6,A,2I6,A)','edge ',J1-1,' atoms ',J2,J3,' changed from true to false'
!              PRINT '(A,3G20.10)','D1SQ,DISQ,c2=',D1,DI,COS2
!           ENDIF
!           EDGEINT(J1-1,J2,J3)=.FALSE.
         ENDIF
      ENDDO
   ENDDO
   GGG(NSHIFT+1:NSHIFT+2*NOPT)=GGG(NSHIFT+1:NSHIFT+2*NOPT)+GLOCAL(1:2*NOPT)
ENDDO
!
! Spring energy. Set EEE(J1) and ESPRING dividing up the pairwise
! energy terms between images except for the end points.
!
ESPRING=0.0D0
IF (KINT.NE.0.0D0) THEN
   DO J1=1,INTIMAGE+1 ! sum over edges from J1 to J1+1
      NI1=NOPT*(J1-1)
      NI2=NOPT*J1
!
!  Edge between J1 and J1+1
!
      DPLUS=0.0D0
      DO J2=1,NATOMS
         DPLUS=DPLUS+(XYZ(NI1+3*(J2-1)+1)-XYZ(NI2+3*(J2-1)+1))**2 &
  &                 +(XYZ(NI1+3*(J2-1)+2)-XYZ(NI2+3*(J2-1)+2))**2 &
  &                 +(XYZ(NI1+3*(J2-1)+3)-XYZ(NI2+3*(J2-1)+3))**2
      ENDDO
      DPLUS=SQRT(DPLUS)
      IF (DPLUS.GT.IMSEPMAX) THEN
         DUMMY=KINT*0.5D0*(DPLUS-IMSEPMAX)**2
         ESPRING=ESPRING+DUMMY
         DUMMY=KINT*(DPLUS-IMSEPMAX)/DPLUS
         DO J2=1,NATOMS
            SPGRAD(1:3)=DUMMY*(XYZ(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)-XYZ(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3))
            GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)=GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)+SPGRAD(1:3)
            GGG(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3)=GGG(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3)-SPGRAD(1:3)
         ENDDO
      ENDIF
   ENDDO
ENDIF
ETOTAL=ETOTAL+ESPRING
!
! Set gradients on frozen atoms to zero.
!
IF (FREEZE) THEN
   DO J1=2,INTIMAGE+1  
      DO J2=1,NATOMS
         IF (FROZEN(J2)) THEN
            GGG(NOPT*(J1-1)+3*(J2-1)+1)=0.0D0
            GGG(NOPT*(J1-1)+3*(J2-1)+2)=0.0D0
            GGG(NOPT*(J1-1)+3*(J2-1)+3)=0.0D0
         ENDIF
      ENDDO
   ENDDO
ENDIF
!
! Set gradients to zero for start and finish images.
!
RMS=0.0D0
RMSIMAGE(1:INTIMAGE+2)=0.0D0
DO J1=2,INTIMAGE+1
   DO J2=1,NOPT
      RMSIMAGE(J1)=RMSIMAGE(J1)+GGG(NOPT*(J1-1)+J2)**2
   ENDDO
   RMS=RMS+RMSIMAGE(J1)
!  PRINT '(A,I6,G20.10,I6)',' intgradlj> J1,square grad,int terms=',J1,RMSIMAGE(J1),NLJINT(J1)
ENDDO
IF (INTIMAGE.GT.0) THEN
   GGG(1:NOPT)=0.0D0
   GGG((INTIMAGE+1)*NOPT+1:(INTIMAGE+2)*NOPT)=0.0D0
ENDIF
IF (INTIMAGE.NE.0) RMS=SQRT(RMS/(NOPT*INTIMAGE))

END SUBROUTINE INTGRADLJ
