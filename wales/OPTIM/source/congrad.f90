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
SUBROUTINE CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
USE KEY, ONLY: FROZEN, FREEZE, NREPI, NREPJ, NNREPULSIVE, &
  &            NCONSTRAINT, CONI, CONJ, INTCONSTRAINTDEL, CONDISTREF, INTCONSTRAINTREP, CONDISTREFLOCAL, &
  &            CONACTIVE, INTCONSTRAINREPCUT, NREPCUT,INTIMAGE, KINT, IMSEPMAX, ATOMACTIVE, JMAXCON, &
  &            INTFREEZET, INTFROZEN, CONCUTLOCAL, CONCUT, CONCUTABST, CONCUTABS, CONCUTFRACT, CONCUTFRAC, &
  &  FREEZENODEST, INTSPRINGACTIVET, INTMINFAC, NCONOFF, CONOFFLIST, CONOFFTRIED, INTCONMAX, ECON, EREP, ESPRING, &
  &  CONVERGECONTEST, CONVERGEREPTEST, FCONTEST, FREPTEST, QCIINTREPMINSEP, QCIAVDEV
USE COMMONS, ONLY: NATOMS, NOPT, DEBUG
USE PORFUNCS
IMPLICIT NONE
           
INTEGER :: J1,J2,NI2,NI1,NJ2,NJ1,ISTAT,MYUNIT,JMAX,IMAX,JMAXNOFF,IMAXNOFF,NMAXINT,NMININT
DOUBLE PRECISION :: ETOTAL, RMS, EMAX, EMAXNOFF, FMAX, FMIN, SEPARATION
INTEGER OFFSET1, OFFSET2
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,D2,D1
DOUBLE PRECISION G1(3),G2(3),DINT,G1INT(3),G2INT(3)
DOUBLE PRECISION DUMMY, REPGRAD(3), D12, DSQ2, DSQ1, DSQI
LOGICAL NOINT
DOUBLE PRECISION XYZ((3*NATOMS)*(INTIMAGE+2)), GGG((3*NATOMS)*(INTIMAGE+2)), EEE(INTIMAGE+2), CCLOCAL
DOUBLE PRECISION GGGR((3*NATOMS)*(INTIMAGE+2)), EREP1,EREP2
LOGICAL IMGFREEZE(INTIMAGE)
DOUBLE PRECISION DPLUS, SPGRAD(3), DCUT, r1amr1bdr2amr2b,r1apr2bmr2amr1bsq,CUTMAX,DISTMAX
DOUBLE PRECISION CONDMAX, CONREFMAX, CONCUTMAX, DUMMY2, CCLOCAL2, DVEC(INTIMAGE+1), DEVIATION(INTIMAGE+1)
DOUBLE PRECISION GLOCAL1(3*NATOMS), GLOCAL2(3*NATOMS), XYZ1(3*NATOMS), XYZ2(3*NATOMS)

EEE(1:INTIMAGE+2)=0.0D0
GGG(1:(3*NATOMS)*(INTIMAGE+2))=0.0D0
ECON=0.0D0; EREP=0.0D0
NMAXINT=0; NMININT=0 ! not actually used
MYUNIT=6

EMAX=-1.0D200
FMAX=-1.0D200
FMIN=1.0D200
! EEMAX(1:INTIMAGE+2)=-1.0D200
! JJMAX(1:INTIMAGE+2)=-1
JMAX=-1
IMAX=-1
EMAXNOFF=-1.0D200
JMAXNOFF=-1
IMAXNOFF=-1
IF (INTCONSTRAINTDEL.EQ.0.0D0) GOTO 531
!
!  Constraint energy and forces.
!
DO J2=1,NCONSTRAINT
   IF (.NOT.CONACTIVE(J2)) CYCLE
      CCLOCAL=CONCUTLOCAL(J2)
      IF (CONCUTABST) CCLOCAL=CCLOCAL+CONCUTABS
      IF (CONCUTFRACT) CCLOCAL=CCLOCAL+CONCUTFRAC*CONDISTREFLOCAL(J2)
!
! For J1 we consider the line segment between image J1-1 and J1.
! There are INTIMAGE+1 line segments in total, with an energy contribution
! and corresponding gradient terms for each. 
! A and B refer to atoms, 2 refers to image J1.
!
   DO J1=2,INTIMAGE+1
!  DO J1=1,INTIMAGE+2  ! checking for zero!
      IF (FREEZENODEST) THEN ! IMGFREEZE is not allocated otherwise!
         IF (J1.EQ.2) THEN
            IF (IMGFREEZE(1)) THEN
!              IF (J2.EQ.1) PRINT '(A)','J1=2 and IMGFREEZE(1)=T cycle'
               CYCLE
            ENDIF
         ELSE IF (J1.EQ.INTIMAGE+2) THEN
            IF (IMGFREEZE(INTIMAGE)) THEN
!              IF (J2.EQ.1) PRINT '(A)','J1=INTIMAGE+2 and IMGFREEZE(INTIMAGE)=T cycle'
               CYCLE
            ENDIF
         ELSE
            IF (IMGFREEZE(J1-2).AND.IMGFREEZE(J1-1)) THEN
!              IF (J2.EQ.1) PRINT '(A,I6,A)','J1=',J1,' IMGFREEZE(J1-2)=T and IMGFREEZE(J1-1)=T cycle'
               CYCLE
            ENDIF
         ENDIF
      ENDIF
      NI1=(3*NATOMS)*(J1-1)+3*(CONI(J2)-1)
      NJ1=(3*NATOMS)*(J1-1)+3*(CONJ(J2)-1)
      R2AX=XYZ(NI1+1); R2AY=XYZ(NI1+2); R2AZ=XYZ(NI1+3)
      R2BX=XYZ(NJ1+1); R2BY=XYZ(NJ1+2); R2BZ=XYZ(NJ1+3)
      D2=SQRT((R2AX-R2BX)**2+(R2AY-R2BY)**2+(R2AZ-R2BZ)**2)
      DUMMY=D2-CONDISTREFLOCAL(J2)  
      DUMMY2=DUMMY**2
      CCLOCAL2=CCLOCAL**2
!
! Reduced form for penalty function and gradient - multiply through by 1/(CCLOCAL**2*INTCONSTRAINTDEL)
! Should save CCLOCAL**2
! Could save DUMMY**2-CCLOCAL**2 outside IF block and base the test on difference of squares
!
      IF (DUMMY2.GT.CCLOCAL2) THEN 
         G2(1)=(R2AX-R2BX)/D2
         G2(2)=(R2AY-R2BY)/D2
         G2(3)=(R2AZ-R2BZ)/D2
!
!        REPGRAD(1:3)=2*INTCONSTRAINTDEL*((DUMMY/CCLOCAL)**2-1.0D0)*DUMMY*G2(1:3)
!        DUMMY=INTCONSTRAINTDEL*(DUMMY**2-CCLOCAL**2)**2/(2.0D0*CCLOCAL**2)
!
         REPGRAD(1:3)=2*(DUMMY2-CCLOCAL2)*DUMMY*G2(1:3)/CCLOCAL2**2
         DUMMY=(DUMMY2-CCLOCAL2)**2/(2.0D0*CCLOCAL2**2)
         EEE(J1)=EEE(J1)  +DUMMY
         ECON=ECON        +DUMMY
         IF (DUMMY.GT.EMAX) THEN
            IMAX=J1
            JMAX=J2
            EMAX=DUMMY
            CONDMAX=D2
            CONREFMAX=CONDISTREFLOCAL(J2)
            CONCUTMAX=CCLOCAL
         ENDIF
         IF (DUMMY.GT.EMAXNOFF) THEN
            IF (.NOT.CONOFFTRIED(J2)) THEN
               IMAXNOFF=J1
               JMAXNOFF=J2
               EMAXNOFF=DUMMY
            ENDIF
         ENDIF
         GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
         GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
         DUMMY2=MINVAL(REPGRAD)
         IF (DUMMY2.LT.FMIN) FMIN=DUMMY
         DUMMY2=MAXVAL(REPGRAD)
         IF (DUMMY2.GT.FMAX) FMAX=DUMMY
      ENDIF
   ENDDO
ENDDO
IF (-FMIN.GT.FMAX) FMAX=-FMIN
FCONTEST=FMAX
IF (JMAX.GT.0) THEN
   WRITE(*,'(A,I6,A,I6,A,2I6,A,G15.5,A,3G15.5,A,G15.5)') ' congrad> Highest constraint for any image in ',IMAX, &
 & ' con ',JMAX, ' atoms ',CONI(JMAX),CONJ(JMAX),' value=',EMAX,' d,ref,cutoff=',CONDMAX,CONREFMAX,CONCUTMAX, &
 & ' max grad=',FMAX

ENDIF
CONVERGECONTEST=EMAX
! IF (JMAXNOFF.GT.0) THEN
!    WRITE(*,'(A,I6,A,I6,A,2I8,A,G20.10,A,I8)') ' congrad> Highest constraint contribution never turned off for any image in image ',IMAXNOFF, &
!  & ' constraint ',JMAXNOFF, &
!  &                              ' atoms ',CONI(JMAX),CONJ(JMAX),' value=',EMAX,' off=',NCONOFF
! ELSEIF (JMAX.GT.0) THEN
!    JMAXNOFF=JMAX
!    WRITE(*,'(A,I6,A,I6,A,2I8,A,I8)') ' congrad> *** Using highest constraint contribution for any image, setting NCONOFF to 0'
!    CONOFFTRIED(1:INTCONMAX)=.FALSE.
!    NCONOFF=0
! ENDIF
! JMAXCON=JMAXNOFF
531 CONTINUE

! GGG(1:(3*NATOMS))=0.0D0                            ! can delete when loop range above changes
! GGG((3*NATOMS)*(INTIMAGE+1)+1:(3*NATOMS)*(INTIMAGE+2))=0.0D0 ! can delete when loop range above changes

! INTCONST=INTCONSTRAINREPCUT**13

EMAX=-1.0D200
FMAX=-1.0D200
FMIN=1.0D200
JMAX=-1
IMAX=-1

GGGR(1:3*NATOMS*(INTIMAGE+2))=0.0D0
IF (INTCONSTRAINTREP.EQ.0.0D0) GOTO 654
DO J1=2,INTIMAGE+1
! DO J1=2,INTIMAGE+2 ! we don't do anything for INTIMAGE+2 any more
!  DO J1=1,INTIMAGE+2 ! can change when zero energies are confirmed for end images
   IF (FREEZENODEST) THEN
      IF (J1.EQ.2) THEN
         IF (IMGFREEZE(1)) CYCLE
!     ELSE IF (J1.EQ.INTIMAGE+2) THEN
!        IF (IMGFREEZE(INTIMAGE)) CYCLE
      ELSE
         IF (IMGFREEZE(J1-2).AND.IMGFREEZE(J1-1)) CYCLE
      ENDIF
   ENDIF
   OFFSET2=(3*NATOMS)*(J1-1)
   OFFSET1=(3*NATOMS)*(J1-2)
   XYZ1(1:3*NATOMS)=XYZ(OFFSET1+1:OFFSET1+3*NATOMS)
   XYZ2(1:3*NATOMS)=XYZ(OFFSET2+1:OFFSET2+3*NATOMS)
   GLOCAL1(1:3*NATOMS)=0.0D0
   GLOCAL2(1:3*NATOMS)=0.0D0
   EREP1=0.0D0
   EREP2=0.0D0
   DO J2=1,NNREPULSIVE
! !  INTCONST=NREPCUT(J2)**13
!    INTCONST=NREPCUT(J2)**3
!      IF (INTFROZEN(NREPI(J2)).AND.INTFROZEN(NREPJ(J2))) THEN
!!        WRITE(*, '(A,I6,A,2I6)') ' congrad> ERROR *** repulsion ',J2,' between frozen atoms ',NREPI(J2),NREPJ(J2)
!         STOP
!      ENDIF

!     NI1=OFFSET1+3*(NREPI(J2)-1)
!     NI2=OFFSET2+3*(NREPI(J2)-1)
!     NJ1=OFFSET1+3*(NREPJ(J2)-1)
!     NJ2=OFFSET2+3*(NREPJ(J2)-1)

      NI1=3*(NREPI(J2)-1)
      NI2=3*(NREPI(J2)-1)
      NJ1=3*(NREPJ(J2)-1)
      NJ2=3*(NREPJ(J2)-1)

!     G1(1:3)=XYZ(NI1+1:NI1+3)-XYZ(NJ1+1:NJ1+3)
!     G2(1:3)=XYZ(NI2+1:NI2+3)-XYZ(NJ2+1:NJ2+3)

      G1(1:3)=XYZ1(NI1+1:NI1+3)-XYZ1(NJ1+1:NJ1+3)
      G2(1:3)=XYZ2(NI2+1:NI2+3)-XYZ2(NJ2+1:NJ2+3)
!
! Squared distance between atoms A and B for theta=0 - distance in image 2
!
      DSQ2=G2(1)**2 + G2(2)**2 + G2(3)**2
!
! Squared distance between atoms A and B for theta=Pi/2 - distance in image 1
!
      DSQ1=G1(1)**2 + G1(2)**2 + G1(3)**2
      DCUT=NREPCUT(J2)**2
      IF ((DSQ1.GT.DCUT).AND.(DSQ2.GT.DCUT)) CYCLE ! don't look for an internal minimum if both repulsions outside cutoff
      IF (ABS(NREPI(J2)-NREPJ(J2)).GT.QCIINTREPMINSEP) THEN ! don't check for internal minimum in distance - atoms too close for chain crossing.
         r1apr2bmr2amr1bsq=0.0D0
      ELSE
!        WRITE(*,'(A,I6,A,I6,A,2G20.10,A,G20.10)') 'distances in ',J1-1,' and ',J1,' are ',SQRT(DSQ1),SQRT(DSQ2),' cutoff=',NREPCUT(J2)
!        WRITE(*,'(A,I6,A,2I6,A,6G15.7)') 'image ',J1-1,' atoms ',NREPI(J2),NREPJ(J2),' coords ',XYZ(NI1+1:NI1+3),XYZ(NJ1+1:NJ1+3)
!        WRITE(*,'(A,I6,A,2I6,A,6G15.7)') 'image ',J1  ,' atoms ',NREPI(J2),NREPJ(J2),' coords ',XYZ(NI2+1:NI2+3),XYZ(NJ2+1:NJ2+3)
         r1amr1bdr2amr2b=G1(1)*G2(1)+G1(2)*G2(2)+G1(3)*G2(3)
         r1apr2bmr2amr1bsq=DSQ1+DSQ2-2.0D0*r1amr1bdr2amr2b
      ENDIF
!
! Is the denominator of the d^2 internal minimum term close to zero?
!
      IF (r1apr2bmr2amr1bsq.LT.1.0D-50) THEN
         NOINT=.TRUE.
         D1=SQRT(DSQ1)
         D2=SQRT(DSQ2)
         G2(1:3)=G2(1:3)/D2
         G1(1:3)=G1(1:3)/D1
      ELSE
         CALL MINMAXD2R(D2,D1,DINT,DSQ2,DSQ1,DSQI,G1,G2,G1INT,G2INT,NOINT,.FALSE.,NREPCUT(J2),r1amr1bdr2amr2b,r1apr2bmr2amr1bsq)
      ENDIF
!
! Skip image INTIMAGE+2 - no non-zero gradients on other images and no energy contributions.
!
! Multiply energy and gradient by NREPCUT**2/INTCONSTRAINTREP to put values on a common scale.
!
!     IF ((D2.LT.NREPCUT(J2)).AND.(J1.LT.INTIMAGE+2)) THEN ! terms for image J1 - non-zero derivatives only for J1
      IF (D2.LT.NREPCUT(J2)) THEN ! terms for image J1 - non-zero derivatives only for J1
!        D12=DSQ2**6
         D12=DSQ2
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*D2-13.0D0*INTCONSTRAINREPCUT)/INTCONST)
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*D2-13.0D0*NREPCUT(J2))/INTCONST)
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(2.0D0*D2-3.0D0*NREPCUT(J2))/INTCONST)
!        DUMMY=(1.0D0/D12+(2.0D0*D2-3.0D0*NREPCUT(J2))/INTCONST)*NREPCUT(J2)**2
         DUMMY=NREPCUT(J2)**2/D12+2.0D0*D2/NREPCUT(J2)-3.0D0  
!        EEE(J1)=EEE(J1)+DUMMY
         EREP1=EREP1+DUMMY
         IF (DUMMY.GT.EMAX) THEN
            IMAX=J1
            JMAX=J2
            CUTMAX=NREPCUT(J2)
            DISTMAX=D2
            EMAX=DUMMY
         ENDIF
         EREP=EREP+DUMMY
!        DUMMY=-12.0D0*INTCONSTRAINTREP*(1.0D0/(D2*D12)-1.0D0/INTCONST)
!        DUMMY=-2.0D0*INTCONSTRAINTREP*(1.0D0/(D2*D12)-1.0D0/INTCONST)
!        DUMMY=-2.0D0*(1.0D0/(D2*D12)-1.0D0/INTCONST)*NREPCUT(J2)**2

         DUMMY=-2.0D0*(NREPCUT(J2)**2/(D2*D12)-1.0D0/NREPCUT(J2))
         REPGRAD(1:3)=DUMMY*G2(1:3)
!        GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
!        GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
         GLOCAL2(NI2+1:NI2+3)=GLOCAL2(NI2+1:NI2+3)+REPGRAD(1:3)
         GLOCAL2(NJ2+1:NJ2+3)=GLOCAL2(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
!
! For internal minima we are counting edges. 
! Edge J1 is between images J1-1 and J1, starting from J1=2.
! Energy contributions are shared evenly, except for
! edge 1, which was assigned to image 2, and edge INTIMAGE+1, which
! was assigned to image INTIMAGE+1. Gradients are set to zero for the end images.
! 20/11/17 - changed to turn off internal minimum checks involving the fixed endpoints. DJW
!
      DUMMY=0.0D0
      IF ((.NOT.NOINT).AND.(DINT.LT.NREPCUT(J2)).AND.(J1.NE.2)) THEN
!        D12=DSQI**6
         D12=DSQI
! !        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*DINT-13.0D0*NREPCUT(J2))/INTCONST)
!        DUMMY=INTMINFAC*INTCONSTRAINTREP*(1.0D0/D12+(2.0D0*DINT-3.0D0*NREPCUT(J2))/INTCONST)
!        DUMMY=INTMINFAC*(1.0D0/D12+(2.0D0*DINT-3.0D0*NREPCUT(J2))/INTCONST)*NREPCUT(J2)**2
         DUMMY=INTMINFAC*(NREPCUT(J2)**2/D12+2.0D0*DINT/NREPCUT(J2)-3.0D0)
!        IF (J1.EQ.2) THEN
!           EEE(J1)=EEE(J1)+DUMMY
!        ELSE ! IF (J1.LT.INTIMAGE+2) THEN
!           EEE(J1)=EEE(J1)+DUMMY/2.0D0
!           EEE(J1-1)=EEE(J1-1)+DUMMY/2.0D0
!        ENDIF
         EREP2=EREP2+DUMMY
         EREP=EREP+DUMMY
         IF (DUMMY.GT.EMAX) THEN
            IMAX=J1
            JMAX=J2
            EMAX=DUMMY
            DISTMAX=DINT
            CUTMAX=NREPCUT(J2)
         ENDIF
! !        DUMMY=-12.0D0*INTCONSTRAINTREP*(1.0D0/(DINT*D12)-1.0D0/INTCONST)
!        DUMMY=-2.0D0*INTCONSTRAINTREP*(1.0D0/(DINT*D12)-1.0D0/INTCONST)
!        DUMMY=-2.0D0*(1.0D0/(DINT*D12)-1.0D0/INTCONST)*NREPCUT(J2)**2
         DUMMY=-2.0D0*(NREPCUT(J2)**2/(DINT*D12)-1.0D0/NREPCUT(J2))
         REPGRAD(1:3)=INTMINFAC*DUMMY*G1INT(1:3)
!        PRINT '(A,4I6,2G15.5)','in1 J1,J2,REPI,REPJ,REPGRAD,NREPCUT=',J1,J2,NREPI(J2),NREPJ(J2), &
! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2),NREPCUT(J2)
!
! Gradient contributions for image J1-1
!

!        GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
!        GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
         GLOCAL1(NI1+1:NI1+3)=GLOCAL1(NI1+1:NI1+3)+REPGRAD(1:3)
         GLOCAL1(NJ1+1:NJ1+3)=GLOCAL1(NJ1+1:NJ1+3)-REPGRAD(1:3)
         DUMMY2=MINVAL(REPGRAD)
!
! Gradient contributions for image J1
!
         REPGRAD(1:3)=INTMINFAC*DUMMY*G2INT(1:3)
!        GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
!        GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
         GLOCAL2(NI2+1:NI2+3)=GLOCAL2(NI2+1:NI2+3)+REPGRAD(1:3)
         GLOCAL2(NJ2+1:NJ2+3)=GLOCAL2(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
   ENDDO
   GGGR(OFFSET1+1:OFFSET1+3*NATOMS)=GGGR(OFFSET1+1:OFFSET1+3*NATOMS)+GLOCAL1(1:3*NATOMS)
   GGGR(OFFSET2+1:OFFSET2+3*NATOMS)=GGGR(OFFSET2+1:OFFSET2+3*NATOMS)+GLOCAL2(1:3*NATOMS)
   EEE(J1)=EEE(J1)+EREP1
   IF (J1.EQ.2) THEN
      EEE(J1)=EEE(J1)+EREP2
   ELSE ! IF (J1.LT.INTIMAGE+2) THEN
      EEE(J1)=EEE(J1)+EREP2/2.0D0
      EEE(J1-1)=EEE(J1-1)+EREP2/2.0D0
   ENDIF
ENDDO
FMIN=MINVAL(GGGR(3*NATOMS+1:3*NATOMS*(INTIMAGE+1)))
FMAX=MAXVAL(GGGR(3*NATOMS+1:3*NATOMS*(INTIMAGE+1)))
IF (-FMIN.GT.FMAX) FMAX=-FMIN
FREPTEST=FMAX

GGG(1:3*NATOMS*(INTIMAGE+2))=GGG(1:3*NATOMS*(INTIMAGE+2))+GGGR(1:3*NATOMS*(INTIMAGE+2))

IF (JMAX.GT.0) THEN
   WRITE(*,'(A,I6,A,I6,A,2I6,A,G15.5,A,2G15.5,A,G15.5)') ' congrad> Highest repulsion  for any image in ',IMAX, &  
 &  ' ind ',JMAX,' atoms ',NREPI(JMAX),NREPJ(JMAX),' value=',EMAX,' d,cutoff=',DISTMAX,CUTMAX, &
 &  ' max grad=',FMAX
ENDIF
CONVERGEREPTEST=EMAX
654 CONTINUE

!
! Spring energy. Set EEE(J1) and ESPRING dividing up the pairwise
! energy terms between images except for the end points.
!
ESPRING=0.0D0
IF (KINT.NE.0.0D0) THEN
   DO J1=1,INTIMAGE+1 ! sum over edges from J1 to J1+1
      NI1=(3*NATOMS)*(J1-1)
      NI2=(3*NATOMS)*J1
!
!  Edge between J1 and J1+1
!
      DPLUS=0.0D0
      DO J2=1,NATOMS
         IF ((.NOT.INTSPRINGACTIVET).OR.ATOMACTIVE(J2)) THEN 
            DPLUS=DPLUS+(XYZ(NI1+3*(J2-1)+1)-XYZ(NI2+3*(J2-1)+1))**2 &
  &                    +(XYZ(NI1+3*(J2-1)+2)-XYZ(NI2+3*(J2-1)+2))**2 &
  &                    +(XYZ(NI1+3*(J2-1)+3)-XYZ(NI2+3*(J2-1)+3))**2
         ENDIF
      ENDDO
      DPLUS=SQRT(DPLUS)
      DVEC(J1)=DPLUS
!     IF (DPLUS.GT.IMSEPMAX) THEN
!        DUMMY=KINT*0.5D0*(DPLUS-IMSEPMAX)**2
         DUMMY=KINT*0.5D0*DPLUS**2
         IF (DUMMY.GT.EMAX) THEN
            IMAX=J1
            EMAX=DUMMY
         ENDIF
         ESPRING=ESPRING+DUMMY
!        IF (J1.EQ.1) THEN
!           EEE(2)=EEE(2)+DUMMY
!        ELSEIF (J1.EQ.INTIMAGE+1) THEN
!           EEE(INTIMAGE+1)=EEE(INTIMAGE+1)+DUMMY
!        ELSE
!           EEE(J1)=EEE(J1)+DUMMY/2.0D0
!           EEE(J1+1)=EEE(J1+1)+DUMMY/2.0D0
!        ENDIF
!        DUMMY=KINT*(DPLUS-IMSEPMAX)/DPLUS
         DUMMY=KINT
         DO J2=1,NATOMS
            IF ((.NOT.INTSPRINGACTIVET).OR.ATOMACTIVE(J2)) THEN 
               SPGRAD(1:3)=DUMMY*(XYZ(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)-XYZ(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3))
               GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)=GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)+SPGRAD(1:3)
               GGG(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3)=GGG(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3)-SPGRAD(1:3)
            ENDIF
         ENDDO
!     ENDIF
   ENDDO
ENDIF
SEPARATION=SUM(DVEC(1:INTIMAGE+1))
DEVIATION(1:INTIMAGE+1)=ABS(100*((INTIMAGE+1)*DVEC(1:INTIMAGE+1)/SEPARATION-1.0D0))
QCIAVDEV=SUM(DEVIATION)/(INTIMAGE+1)
WRITE(*,'(A,I6,A,I6,A,2I6)') ' congrad> Highest spring  contribution for any image in image ',IMAX
IF (DEBUG) THEN
   WRITE(*, '(A,3G20.10)') ' congrad> ECON,EREP,ESPRING=',ECON,EREP,ESPRING
!  WRITE(*,'(A)') '   edge         gap                deviation      '
!  WRITE(*,'(I6,3X,G20.10,G20.10)') (J1,DVEC(J1),DEVIATION(J1),J1=1,INTIMAGE+1)
   WRITE(*, '(A,2G20.10)') ' congrad> mean gap and mean deviation=',SEPARATION/(INTIMAGE+1),QCIAVDEV
ENDIF
!
! Set gradients on frozen atoms to zero.
!
IF (FREEZE) THEN
   DO J1=2,INTIMAGE+1  
      DO J2=1,NATOMS
         IF (FROZEN(J2)) THEN
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+1)=0.0D0
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+2)=0.0D0
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+3)=0.0D0
         ENDIF
      ENDDO
   ENDDO
ENDIF
!
! Set gradients on locally frozen atoms to zero.
!
IF (INTFREEZET) THEN
   DO J1=2,INTIMAGE+1  
      DO J2=1,NATOMS
         IF (INTFROZEN(J2)) THEN
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+1)=0.0D0
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+2)=0.0D0
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+3)=0.0D0
         ENDIF
      ENDDO
   ENDDO
ENDIF
!
! Set gradients to zero for start and finish images.
!
GGG(1:(3*NATOMS))=0.0D0
GGG((INTIMAGE+1)*(3*NATOMS)+1:(INTIMAGE+2)*(3*NATOMS))=0.0D0
RMS=0.0D0
DO J1=2,INTIMAGE+1
   DO J2=1,(3*NATOMS)
      RMS=RMS+GGG((3*NATOMS)*(J1-1)+J2)**2
   ENDDO
ENDDO
IF (INTIMAGE.NE.0) THEN
   RMS=SQRT(RMS/((3*NATOMS)*INTIMAGE))
ENDIF
!
! For INTIMAGE images there are INTIMAGE+2 replicas including the end points,
! and INTIMAGE+1 line segements, with associated energies stored in EEE(2:INTIMAGE+2)
!
ETOTAL=0.0D0
DO J1=2,INTIMAGE+1
   ETOTAL=ETOTAL+EEE(J1)
ENDDO

END SUBROUTINE CONGRAD

SUBROUTINE MINMAXD2(D2,D1,DINT,DSQ2,DSQ1,G1,G2,G1INT,G2INT,NOINT,DEBUG,r1amr1bdr2amr2b,r1apr2bmr2amr1bsq)
USE KEY, ONLY : CHECKCONINT
IMPLICIT NONE
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,D2,D1,DINT
DOUBLE PRECISION G1(3),G2(3),G1INT(3),G2INT(3)
DOUBLE PRECISION DSQ2, DSQ1, DSQI, r1apr2bmr2amr1bsq, r1amr1bsq, r2amr2bsq
DOUBLE PRECISION r1amr1bdr2amr2b, r1amr1bdr2amr2bsq, DUMMY, DUMMY2
LOGICAL NOINT, DEBUG

!
! Is there an internal extremum?
! THIS TEST IS NOT NEEDED IF CHECKCONINT IS FALSE. SHOULD SKIP.
!
! IF (r1apr2bmr2amr1bsq.EQ.0.0D0) THEN ! now done in calling routine
!
DUMMY=(DSQ1-r1amr1bdr2amr2b)/r1apr2bmr2amr1bsq
NOINT=.TRUE.
IF ((DUMMY.GT.0.0D0).AND.(DUMMY.LT.1.0D0)) NOINT=.FALSE.
G1INT(1:3)=0.0D0
G2INT(1:3)=0.0D0
D2=SQRT(DSQ2)
D1=SQRT(DSQ1)
DSQI=1.0D10
DINT=1.0D10
IF (.NOT.NOINT) THEN
   r1amr1bdr2amr2bsq=r1amr1bdr2amr2b**2
   DUMMY2=r1amr1bdr2amr2bsq - DSQ1*DSQ2
   DSQI=MAX(-DUMMY2/r1apr2bmr2amr1bsq,0.0D0)
   DUMMY=r1apr2bmr2amr1bsq**2
   DINT=SQRT(DSQI)
   IF (DINT.LE.0.0D0) THEN
      NOINT=.TRUE.
   ELSE
     DUMMY2=r1amr1bdr2amr2bsq - DSQ1*DSQ2
     DUMMY=DUMMY*DINT ! Convert derivatives of distance^2 to derivative of distance.
     G1INT(1:3)= (DUMMY2*(G1(1:3) - G2(1:3)) + r1apr2bmr2amr1bsq*(G1(1:3)*DSQ2 -G2(1:3)*r1amr1bdr2amr2b))/DUMMY
     G2INT(1:3)= (DUMMY2*(G2(1:3) - G1(1:3)) + r1apr2bmr2amr1bsq*(G2(1:3)*DSQ1 -G1(1:3)*r1amr1bdr2amr2b))/DUMMY
   ENDIF
ENDIF
!
! Convert derivatives of distance^2 to derivative of distance.
! We have cancelled a factor of two above and below!
!
G2(1:3)=G2(1:3)/D2
G1(1:3)=G1(1:3)/D1
! IF (.NOT.NOINT) THEN
!    G1INT(1:3)=G1INT(1:3)/DINT
!    G2INT(1:3)=G2INT(1:3)/DINT
! ENDIF

END SUBROUTINE MINMAXD2

SUBROUTINE MINMAXD2R(D2,D1,DINT,DSQ2,DSQ1,DSQI,G1,G2,G1INT,G2INT,NOINT,DEBUG,INTCONSTRAINREPCUT,r1amr1bdr2amr2b,r1apr2bmr2amr1bsq)   
IMPLICIT NONE
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,D2,D1,DINT, DUMMY2
DOUBLE PRECISION G1(3),G2(3),G1INT(3),G2INT(3),INTCONSTRAINREPCUT
DOUBLE PRECISION DSQ2, DSQ1, DSQI, r1apr2bmr2amr1bsq, r1amr1bsq, r2amr2bsq
DOUBLE PRECISION r1amr1bdr2amr2b, r1amr1bdr2amr2bsq, DUMMY
LOGICAL NOINT, DEBUG

DUMMY=(DSQ1-r1amr1bdr2amr2b)/r1apr2bmr2amr1bsq
NOINT=.TRUE.
IF ((DUMMY.GT.0.0D0).AND.(DUMMY.LT.1.0D0)) NOINT=.FALSE.
G1INT(1:3)=0.0D0
G2INT(1:3)=0.0D0
D2=SQRT(DSQ2)
D1=SQRT(DSQ1)
DSQI=1.0D10
DINT=1.0D10
IF (.NOT.NOINT) THEN
   r1amr1bdr2amr2bsq=r1amr1bdr2amr2b**2
   DUMMY2=r1amr1bdr2amr2bsq - DSQ1*DSQ2
   DSQI=MAX(-DUMMY2/r1apr2bmr2amr1bsq,0.0D0)
   DUMMY=r1apr2bmr2amr1bsq**2
   DINT=SQRT(DSQI)
   IF (DINT.LE.0.0D0) THEN
      NOINT=.TRUE.
   ELSEIF (DINT.LE.INTCONSTRAINREPCUT) THEN ! skip otherwise
      DUMMY2=r1amr1bdr2amr2bsq - DSQ1*DSQ2
      DUMMY=DUMMY*DINT ! to convert derivatives of distance^2 to derivative of distance.
      G1INT(1:3)= (DUMMY2*(G1(1:3) - G2(1:3)) + r1apr2bmr2amr1bsq*(G1(1:3)*DSQ2 -G2(1:3)*r1amr1bdr2amr2b))/DUMMY
      G2INT(1:3)= (DUMMY2*(G2(1:3) - G1(1:3)) + r1apr2bmr2amr1bsq*(G2(1:3)*DSQ1 -G1(1:3)*r1amr1bdr2amr2b))/DUMMY
   ENDIF
ENDIF
!
! Convert derivatives of distance^2 to derivative of distance.
! We have cancelled a factor of two above and below!
!
G2(1:3)=G2(1:3)/D2
G1(1:3)=G1(1:3)/D1

END SUBROUTINE MINMAXD2R

!
! This version of congrad tests for an internal minimum in the
! constraint distances as well as the repulsions.
!
SUBROUTINE CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
USE KEY, ONLY: FROZEN, FREEZE, NREPI, NREPJ, NNREPULSIVE, &
  &            NCONSTRAINT, CONI, CONJ, INTCONSTRAINTDEL, CONDISTREF, INTCONSTRAINTREP, CONDISTREFLOCAL, &
  &            CONACTIVE, INTCONSTRAINREPCUT, NREPCUT,FREEZENODEST, INTIMAGE, ATOMACTIVE, KINT, IMSEPMAX, &
  &            INTFREEZET, INTFROZEN, REPI, REPJ, CONCUT, CONCUTLOCAL, &
  &            CONCUTABS, CONCUTABST, CONCUTFRAC, CONCUTFRACT, INTMINFAC, INTSPRINGACTIVET, CHECKCONINT, JMAXCON, &
  &            NCONOFF, CONOFFTRIED, INTCONMAX, KINTENDS, ECON, EREP, ESPRING, &
  &  CONVERGECONTEST, CONVERGEREPTEST
USE COMMONS, ONLY: NATOMS, NOPT, DEBUG
IMPLICIT NONE
           
INTEGER :: J1,J2,NI2,NI1,NJ2,NJ1,NMAXINT,NMININT,NCONINT(INTIMAGE+2),NREPINT(INTIMAGE+2),JMAX,IMAX,JMAXNOFF,IMAXNOFF
DOUBLE PRECISION :: ETOTAL, RMS, EMAX, EMAXNOFF
INTEGER JJMAX(INTIMAGE+2)
DOUBLE PRECISION  EEMAX(INTIMAGE+2)
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,D2,D1
DOUBLE PRECISION G1(3),G2(3),DINT,G1INT(3),G2INT(3)
DOUBLE PRECISION DUMMY, REPGRAD(3), INTCONST, D12, DSQ2, DSQ1, DSQI
DOUBLE PRECISION CONE(INTIMAGE+2), REPE(INTIMAGE+2),MAXINT,MININT,REPEINT(INTIMAGE+2),CONEINT(INTIMAGE+2),RMSIMAGE(INTIMAGE+2)
LOGICAL NOINT, LPRINT
DOUBLE PRECISION XYZ((3*NATOMS)*(INTIMAGE+2)), GGG((3*NATOMS)*(INTIMAGE+2)), EEE(INTIMAGE+2)
LOGICAL IMGFREEZE(INTIMAGE), PRINTE
DOUBLE PRECISION DPLUS, SPGRAD(3), CCLOCAL, DCUT, r1amr1bdr2amr2b,r1apr2bmr2amr1bsq
DOUBLE PRECISION CONDMAX, CONREFMAX, CONCUTMAX

PRINTE=.FALSE.
111 CONTINUE

EEE(1:INTIMAGE+2)=0.0D0
CONE(1:INTIMAGE+2)=0.0D0
REPE(1:INTIMAGE+2)=0.0D0
NCONINT(1:INTIMAGE+2)=0
NREPINT(1:INTIMAGE+2)=0
REPEINT(1:INTIMAGE+2)=0.0D0
CONEINT(1:INTIMAGE+2)=0.0D0
GGG(1:(3*NATOMS)*(INTIMAGE+2))=0.0D0
ECON=0.0D0; EREP=0.0D0
LPRINT=.TRUE.
LPRINT=.FALSE.
!
!  Constraint energy and forces.
!
! For J1 we consider the line segment between image J1-1 and J1.
! There are INTIMAGE+1 line segments in total, with an energy contribution
! and corresponding gradient terms for each. 
! A and B refer to atoms, 1 and 2 to images J1-1 and J1 corresponding to J1-2 and J1-1 below.
!
! IMGFREEZE(1:INTIMAGE) refers to the images excluding end points!
!
EMAX=-1.0D200
EMAXNOFF=-1.0D200
EEMAX(1:INTIMAGE+2)=-1.0D200
JJMAX(1:INTIMAGE+2)=-1
JMAX=-1
IMAX=-1
JMAXNOFF=-1
IMAXNOFF=-1
DO J2=1,NCONSTRAINT
   IF (.NOT.CONACTIVE(J2)) CYCLE
   CCLOCAL=CONCUTLOCAL(J2)
   IF (CONCUTABST) CCLOCAL=CCLOCAL+CONCUTABS
   IF (CONCUTFRACT) CCLOCAL=CCLOCAL+CONCUTFRAC*CONDISTREFLOCAL(J2)
!!!!!!!!!!!!!!!!!!!!!!!!!! DEBUG
!  IF (INTFROZEN(CONI(J2)).AND.INTFROZEN(CONJ(J2))) THEN
!     WRITE(*, '(A,I6,A,2I6)') ' congrad2> ERROR *** constraint ',J2,' between frozen atoms ',CONI(J2),CONJ(J2)
!     STOP
!  ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!! DEBUG
   DO J1=2,INTIMAGE+2
      IF (FREEZENODEST) THEN ! IMGFREEZE is not allocated otherwise!
         IF (J1.EQ.2) THEN
            IF (IMGFREEZE(1)) THEN
!              IF (J2.EQ.1) PRINT '(A)','J1=2 and IMGFREEZE(1)=T cycle'
               CYCLE
            ENDIF
         ELSE IF (J1.EQ.INTIMAGE+2) THEN
            IF (IMGFREEZE(INTIMAGE)) THEN
!              IF (J2.EQ.1) PRINT '(A)','J1=INTIMAGE+2 and IMGFREEZE(INTIMAGE)=T cycle'
               CYCLE
            ENDIF
         ELSE
            IF (IMGFREEZE(J1-2).AND.IMGFREEZE(J1-1)) THEN
!              IF (J2.EQ.1) PRINT '(A,I6,A)','J1=',J1,' IMGFREEZE(J1-2)=T and IMGFREEZE(J1-1)=T cycle'
               CYCLE
            ENDIF
         ENDIF
      ENDIF
      NI1=(3*NATOMS)*(J1-2)+3*(CONI(J2)-1)
      NI2=(3*NATOMS)*(J1-1)+3*(CONI(J2)-1)
      NJ1=(3*NATOMS)*(J1-2)+3*(CONJ(J2)-1)
      NJ2=(3*NATOMS)*(J1-1)+3*(CONJ(J2)-1)

      G1(1:3)=XYZ(NI1+1:NI1+3)-XYZ(NJ1+1:NJ1+3)
      G2(1:3)=XYZ(NI2+1:NI2+3)-XYZ(NJ2+1:NJ2+3)
!
! Squared distance between atoms A and B for theta=0 - distance in image 2
!
      DSQ2=G2(1)**2 + G2(2)**2 + G2(3)**2
!
! Squared distance between atoms A and B for theta=Pi/2 - distance in image 1
!
      DSQ1=G1(1)**2 + G1(2)**2 + G1(3)**2
      r1amr1bdr2amr2b=G1(1)*G2(1)+G1(2)*G2(2)+G1(3)*G2(3)
!
! Is there an internal extremum?
!
      r1apr2bmr2amr1bsq=DSQ1+DSQ2-2.0D0*r1amr1bdr2amr2b

      IF ((.NOT.CHECKCONINT).OR.(r1apr2bmr2amr1bsq.LT.1.0D-50)) THEN
         NOINT=.TRUE.
         D1=SQRT(DSQ1)
         D2=SQRT(DSQ2)
         G2(1:3)=G2(1:3)/D2
         G1(1:3)=G1(1:3)/D1
      ELSE
         CALL MINMAXD2(D2,D1,DINT,DSQ2,DSQ1,G1,G2,G1INT,G2INT,NOINT,DEBUG,r1amr1bdr2amr2b,r1apr2bmr2amr1bsq)
      ENDIF
!
! Need to include both D2 and D1 contributions if they are both outside tolerance.
! Otherwise we get discontinuities if they are very close and swap over.
!
!     CONCUT=CONCUTFRAC*CONDISTREF(J2)
!
! terms for image J1 - non-zero derivatives only for J1. D2 is the distance for image J1.
!
!     IF (LPRINT) WRITE(*, '(A,I6,5G15.5)') &
! &       'J1,D2,D1,DINT,MIN diff,CONCUT=',J1,D2,D1,DINT,ABS(D2-CONDISTREFLOCAL(J2)),CCLOCAL
      IF ((ABS(D2-CONDISTREFLOCAL(J2)).GT.CCLOCAL).AND.(J1.LT.INTIMAGE+2)) THEN 
         DUMMY=D2-CONDISTREFLOCAL(J2)  
         REPGRAD(1:3)=2*INTCONSTRAINTDEL*((DUMMY/CCLOCAL)**2-1.0D0)*DUMMY*G2(1:3)
         DUMMY=INTCONSTRAINTDEL*(DUMMY**2-CCLOCAL**2)**2/(2.0D0*CCLOCAL**2)
         IF (DUMMY.GT.EMAX) THEN
            IMAX=J1
            JMAX=J2
            EMAX=DUMMY
            CONDMAX=D2
            CONREFMAX=CONDISTREFLOCAL(J2)
            CONCUTMAX=CCLOCAL
         ENDIF
         IF (DUMMY.GT.EMAXNOFF) THEN
            IF (.NOT.CONOFFTRIED(J2)) THEN
               IMAXNOFF=J1
               JMAXNOFF=J2
               EMAXNOFF=DUMMY
            ENDIF
         ENDIF
         IF (DUMMY.GT.EEMAX(J1)) THEN
            JJMAX(J1)=J2
            EEMAX(J1)=DUMMY
         ENDIF
         EEE(J1)=EEE(J1)+DUMMY
         CONE(J1)=CONE(J1)+DUMMY
         ECON=ECON      +DUMMY
         IF (LPRINT) WRITE(*, '(A,4I6,G15.5)') 'min J1,J2,CONI,CONJ,REPGRAD=',J1,J2,CONI(J2),CONJ(J2), &
  &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2)
         GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
         GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
! !
! ! Don't add energy contributions to EEE(2) from D1, since the gradients are non-zero only for image 1.
! !
! ! terms for image J1-1 - non-zero derivatives only for J1-1. D1 is the distance for image J1-1.
! !
! !     IF (LPRINT) WRITE(*, '(A,I6,5G15.5)') &
! ! &       'J1,D2,D1,DINT,MAX diff,CCLOCAL=',J1,D2,D1,DINT,ABS(D1-CONDISTREFLOCAL(J2)),CCLOCAL
!       IF ((ABS(D1-CONDISTREFLOCAL(J2)).GT.CCLOCAL).AND.(J1.GT.2)) THEN  
!          DUMMY=D1-CONDISTREFLOCAL(J2)  
!          REPGRAD(1:3)=2*INTCONSTRAINTDEL*((DUMMY/CCLOCAL)**2-1.0D0)*DUMMY*G1(1:3)
!          DUMMY=INTCONSTRAINTDEL*(DUMMY**2-CCLOCAL**2)**2/(2.0D0*CCLOCAL**2)
! !        IF (PRINTE.AND.(DUMMY.GT.1.0D-4)) THEN
! !           WRITE(*, '(A,2I6,2L5,G20.10)') 'A CONI,CONJ,INTFROZEN(CONI),INTFROZEN(CONJ),DUMMY=', &
! ! &                                       CONI(J2),CONJ(J2),INTFROZEN(CONI(J2)),INTFROZEN(CONJ(J2)),DUMMY
! !        ENDIF
!          IF (DUMMY.GT.EMAX) THEN
!             IMAX=J1
!             JMAX=J2
!             EMAX=DUMMY
!          ENDIF
!          IF (DUMMY.GT.EMAXNOFF) THEN
!             IF (.NOT.CONOFFTRIED(J2)) THEN
!                IMAXNOFF=J1
!                JMAXNOFF=J2
!                EMAXNOFF=DUMMY
!             ENDIF
!          ENDIF
!          IF (DUMMY.GT.EEMAX(J1-1)) THEN
!             JJMAX(J1-1)=J2
!             EEMAX(J1-1)=DUMMY
!          ENDIF
!          EEE(J1-1)=EEE(J1-1)+DUMMY
!          CONE(J1-1)=CONE(J1-1)+DUMMY
!          ECON=ECON      +DUMMY
!          IF (LPRINT) WRITE(*, '(A,4I6,G15.5)') 'max J1,J2,CONI,CONJ,REPGRAD=',J1,J2,CONI(J2),CONJ(J2), &
!   &         SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2)
!          GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
!          GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
!       ENDIF
!
      IF (CHECKCONINT.AND.(.NOT.NOINT).AND.(ABS(DINT-CONDISTREFLOCAL(J2)).GT.CCLOCAL)) THEN
         DUMMY=DINT-CONDISTREFLOCAL(J2)  
         REPGRAD(1:3)=2*INTMINFAC*INTCONSTRAINTDEL*((DUMMY/CCLOCAL)**2-1.0D0)*DUMMY*G1INT(1:3)
         GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
         GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
         REPGRAD(1:3)=2*INTMINFAC*INTCONSTRAINTDEL*((DUMMY/CCLOCAL)**2-1.0D0)*DUMMY*G2INT(1:3)
         DUMMY=INTMINFAC*INTCONSTRAINTDEL*(DUMMY**2-CCLOCAL**2)**2/(2.0D0*CCLOCAL**2)
!        WRITE(*,'(A,3I7,9F13.5)') 'J1,NI1,NJ1,INTMINFAC,INTCONSTRAINTDEL,DUMMY,GGG=',J1,NI1,NJ1,INTMINFAC,INTCONSTRAINTDEL, &
! &            DUMMY, GGG(NI1+1:NI1+3),GGG(NJ1+1:NJ1+3)   
!        IF (PRINTE.AND.(DUMMY.GT.1.0D-4)) THEN
!           WRITE(*, '(A,2I6,2L5,G20.10)') 'CONI,CONJ,INTFROZEN(CONI),INTFROZEN(CONJ),DUMMY=', &
! &                                       CONI(J2),CONJ(J2),INTFROZEN(CONI(J2)),INTFROZEN(CONJ(J2)),DUMMY
!        ENDIF
         ECON=ECON+DUMMY
         IF (DUMMY.GT.EMAX) THEN
            IMAX=J1
            JMAX=J2
            EMAX=DUMMY
            CONDMAX=DINT
            CONREFMAX=CONDISTREFLOCAL(J2)
            CONCUTMAX=CCLOCAL
         ENDIF
         IF (DUMMY.GT.EMAXNOFF) THEN
            IF (.NOT.CONOFFTRIED(J2)) THEN
               IMAXNOFF=J1
               JMAXNOFF=J2
               EMAXNOFF=DUMMY
            ENDIF
         ENDIF
         IF (DUMMY.GT.EEMAX(J1-1)) THEN
            JJMAX(J1-1)=J2
            EEMAX(J1-1)=DUMMY
         ENDIF
         IF (DUMMY.GT.EEMAX(J1)) THEN
            JJMAX(J1)=J2
            EEMAX(J1)=DUMMY
         ENDIF
         IF (J1.EQ.2) THEN
            EEE(J1)=EEE(J1)+DUMMY
            CONEINT(J1)=CONEINT(J1)+DUMMY
            NCONINT(J1)=NCONINT(J1)+1
         ELSE IF (J1.LT.INTIMAGE+2) THEN
            EEE(J1)=EEE(J1)+DUMMY/2.0D0
            EEE(J1-1)=EEE(J1-1)+DUMMY/2.0D0
            CONEINT(J1)=CONEINT(J1)+DUMMY/2.0D0
            CONEINT(J1-1)=CONEINT(J1-1)+DUMMY/2.0D0
            NCONINT(J1)=NCONINT(J1)+1
            NCONINT(J1-1)=NCONINT(J1-1)+1
         ELSE IF (J1.EQ.INTIMAGE+2) THEN
            EEE(J1-1)=EEE(J1-1)+DUMMY
            CONEINT(J1-1)=CONEINT(J1-1)+DUMMY
            NCONINT(J1-1)=NCONINT(J1-1)+1
         ENDIF
!        WRITE(*, '(A,4I6,G15.5)') 'in2 J1,J2,CONI,CONJ,REPGRAD=',J1,J2,CONI(J2),CONJ(J2), &
! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2)
         GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
         GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
!        WRITE(*, '(A,3I7,9G13.5)') 'J1,NI2,NJ2,INTMINFAC,INTCONSTRAINTDEL,DUMMY,GGG=',J1,NI2,NJ2,INTMINFAC,INTCONSTRAINTDEL, & 
! &            DUMMY, GGG(NI2+1:NI2+3),GGG(NJ2+1:NJ2+3)   
!        WRITE(*,'(A,2G20.10)') 'in intmin block EEE(J1),EEE(J1-1)=',EEE(J1),EEE(J1-1)
      ENDIF
   ENDDO
ENDDO
IF (JMAX.GT.0) THEN
   WRITE(*,'(A,I6,A,I6,A,2I6,A,G20.10,A,3G20.10)') ' congrad> Highest constraint contribution for any image in image ',IMAX, &
 & ' constraint ',JMAX, &
 &                              ' atoms ',CONI(JMAX),CONJ(JMAX),' value=',EMAX,' d,ref,cutoff=',CONDMAX,CONREFMAX,CONCUTMAX

ENDIF
CONVERGECONTEST=EMAX/INTCONSTRAINTDEL
! IF (JMAXNOFF.GT.0) THEN
!    WRITE(*,'(A,I6,A,I6,A,2I8,A,I8)') ' congrad2> Highest constraint contribution never turned off for any image in image ',IMAXNOFF, &
!  & ' constraint ',JMAXNOFF, &
!  &                              ' atoms ',CONI(JMAX),CONJ(JMAX),' off=',NCONOFF
! ELSEIF (JMAX.GT.0) THEN
!    JMAXNOFF=JMAX
!    WRITE(*,'(A,I6,A,I6,A,2I8,A,I8)') ' congrad2> *** Using highest constraint contribution for any image, setting NCONOFF to 0'
!    CONOFFTRIED(1:INTCONMAX)=.FALSE.
!    NCONOFF=0
! ENDIF
JMAXCON=JMAXNOFF

! INTCONST=INTCONSTRAINREPCUT**13

EMAX=-1.0D200
EEMAX(1:INTIMAGE+2)=-1.0D200
JJMAX(1:INTIMAGE+2)=-1
JMAX=-1
IMAX=-1
DO J2=1,NNREPULSIVE
!  INTCONST=NREPCUT(J2)**13
   INTCONST=NREPCUT(J2)**3
   DO J1=2,INTIMAGE+2
      IF (FREEZENODEST) THEN
         IF (J1.EQ.2) THEN
            IF (IMGFREEZE(1)) CYCLE
         ELSE IF (J1.EQ.INTIMAGE+2) THEN
            IF (IMGFREEZE(INTIMAGE)) CYCLE
         ELSE
            IF (IMGFREEZE(J1-2).AND.IMGFREEZE(J1-1)) CYCLE
         ENDIF
      ENDIF
!     IF (INTFROZEN(NREPI(J2)).AND.INTFROZEN(NREPJ(J2))) THEN
!        WRITE(*, '(A,I6,A,2I6)') ' congrad2> ERROR *** repulsion ',J2,' between frozen atoms ',NREPI(J2),NREPJ(J2)
!        STOP
!     ENDIF
      NI1=(3*NATOMS)*(J1-2)+3*(NREPI(J2)-1)
      NI2=(3*NATOMS)*(J1-1)+3*(NREPI(J2)-1)
      NJ1=(3*NATOMS)*(J1-2)+3*(NREPJ(J2)-1)
      NJ2=(3*NATOMS)*(J1-1)+3*(NREPJ(J2)-1)

      G1(1:3)=XYZ(NI1+1:NI1+3)-XYZ(NJ1+1:NJ1+3) 
      G2(1:3)=XYZ(NI2+1:NI2+3)-XYZ(NJ2+1:NJ2+3) 
!
! Squared distance between atoms A and B for theta=0 - distance in image 2
!
      DSQ2=G2(1)**2 + G2(2)**2 + G2(3)**2
!
! Squared distance between atoms A and B for theta=Pi/2 - distance in image 1
!
      DSQ1=G1(1)**2 + G1(2)**2 + G1(3)**2
      DCUT=NREPCUT(J2)**2
      IF ((DSQ1.GT.DCUT).AND.(DSQ2.GT.DCUT)) CYCLE ! don't look for an internal minimum if both repulsions outside cutoff
      r1amr1bdr2amr2b=G1(1)*G2(1)+G1(2)*G2(2)+G1(3)*G2(3)
!
! Is there an internal extremum?
!
      r1apr2bmr2amr1bsq=DSQ1+DSQ2-2.0D0*r1amr1bdr2amr2b

      IF (r1apr2bmr2amr1bsq.LT.1.0D-50) THEN
!        D1=1.0D100; D2=1.0D100; 
         NOINT=.TRUE.  
         D1=SQRT(DSQ1)
         D2=SQRT(DSQ2)
         G2(1:3)=G2(1:3)/D2
         G1(1:3)=G1(1:3)/D1
      ELSE
         CALL MINMAXD2R(D2,D1,DINT,DSQ2,DSQ1,DSQI,G1,G2,G1INT,G2INT,NOINT,.FALSE.,NREPCUT(J2),r1amr1bdr2amr2b,r1apr2bmr2amr1bsq)
      ENDIF
!??????????????????????????????????????????????????????????????????????????????
!
!  WHY ARE WE DOING both D2 and D1? Isn't this double counting? See CONGRAD routine above
!
!
! Skip image INTIMAGE+2 - no non-zero gradients on other images and no energy contributions.
!
      IF ((D2.LT.NREPCUT(J2)).AND.(J1.LT.INTIMAGE+2)) THEN ! terms for image J1 - non-zero derivatives only for J1
!        D12=DSQ2**6
         D12=DSQ2
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*D2-13.0D0*INTCONSTRAINREPCUT)/INTCONST)
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*D2-13.0D0*NREPCUT(J2))/INTCONST)
         DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(2.0D0*D2-3.0D0*NREPCUT(J2))/INTCONST)
         EEE(J1)=EEE(J1)+DUMMY
!        IF (PRINTE.AND.(DUMMY.GT.1.0D-4)) THEN
!           WRITE(*, '(A,2I6,2L5,G20.10)') 'R1 NREPI,NREPJ,INTFROZEN(NREPI),INTFROZEN(NREPJ),DUMMY=', &
! &                                     NREPI(J2),NREPJ(J2),INTFROZEN(NREPI(J2)),INTFROZEN(NREPJ(J2)),DUMMY
!        ENDIF
         IF (DUMMY.GT.EMAX) THEN
            IMAX=J1
            JMAX=J2
            EMAX=DUMMY
         ENDIF
         IF (DUMMY.GT.EEMAX(J1)) THEN
            JJMAX(J1)=J2
            EEMAX(J1)=DUMMY
         ENDIF
         REPE(J1)=REPE(J1)+DUMMY
         EREP=EREP+DUMMY
!        DUMMY=-12.0D0*INTCONSTRAINTREP*(1.0D0/(D2*D12)-1.0D0/INTCONST)
         DUMMY=-2.0D0*INTCONSTRAINTREP*(1.0D0/(D2*D12)-1.0D0/INTCONST)
         REPGRAD(1:3)=DUMMY*G2(1:3)
!        WRITE(*, '(A,4I6,G15.5)') 'min J1,J2,REPI,REPJ,REPGRAD=',J1,J2,NREPI(J2),NREPJ(J2), &
! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2)
         GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
         GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
      DUMMY=0.0D0
!
! SURELY THIS BLOCK JUST DUPLICATES THE TERMS FOR THE D2 BLOCK?
!
! !
! ! Don't add energy contributions to EEE(2) from D1, since the gradients are non-zero only for image 1.
! !
! !     IF ((D1.LT.INTCONSTRAINREPCUT).AND.(J1.GT.2)) THEN ! terms for image J1-1 - non-zero derivatives only for J1-1
!       IF ((D1.LT.NREPCUT(J2)).AND.(J1.GT.2)) THEN ! terms for image J1-1 - non-zero derivatives only for J1-1
! !        D12=DSQ1**6
!          D12=DSQ1
! !        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*D1-13.0D0*INTCONSTRAINREPCUT)/INTCONST)
! !        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*D1-13.0D0*NREPCUT(J2))/INTCONST)
!          DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(2.0D0*D1-3.0D0*NREPCUT(J2))/INTCONST)
!          EEE(J1-1)=EEE(J1-1)+DUMMY
! !        IF (PRINTE.AND.(DUMMY.GT.1.0D-4)) THEN
! !           WRITE(*, '(A,2I6,2L5,G20.10)') 'R2 NREPI,NREPJ,INTFROZEN(NREPI),INTFROZEN(NREPJ),DUMMY=', &
! ! &                                     NREPI(J2),NREPJ(J2),INTFROZEN(NREPI(J2)),INTFROZEN(NREPJ(J2)),DUMMY
! !        ENDIF
!          IF (DUMMY.GT.EMAX) THEN
!             IMAX=J1
!             JMAX=J2
!             EMAX=DUMMY
!          ENDIF
!          IF (DUMMY.GT.EEMAX(J1-1)) THEN
!             JJMAX(J1-1)=J2
!             EEMAX(J1-1)=DUMMY
!          ENDIF
!          REPE(J1-1)=REPE(J1-1)+DUMMY
!          EREP=EREP+DUMMY
! !        DUMMY=-12.0D0*INTCONSTRAINTREP*(1.0D0/(D1*D12)-1.0D0/INTCONST)
!          DUMMY=-2.0D0*INTCONSTRAINTREP*(1.0D0/(D1*D12)-1.0D0/INTCONST)
!          REPGRAD(1:3)=DUMMY*G1(1:3)
! !        WRITE(*, '(A,4I6,G15.5)') 'max J1,J2,REPI,REPJ,REPGRAD=',J1,J2,NREPI(J2),NREPJ(J2), &
! ! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2)
!          GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
!          GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
!       ENDIF
      DUMMY=0.0D0
!     IF ((.NOT.NOINT).AND.(DINT.LT.INTCONSTRAINREPCUT)) THEN
      IF ((.NOT.NOINT).AND.(DINT.LT.NREPCUT(J2))) THEN
!        D12=DSQI**6
         D12=DSQI
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*DINT-13.0D0*INTCONSTRAINREPCUT)/INTCONST)
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*DINT-13.0D0*NREPCUT(J2))/INTCONST)
         DUMMY=INTMINFAC*INTCONSTRAINTREP*(1.0D0/D12+(2.0D0*DINT-3.0D0*NREPCUT(J2))/INTCONST)
         EREP=EREP+DUMMY
!        IF (PRINTE.AND.(DUMMY.GT.1.0D-4)) THEN
!           WRITE(*, '(A,2I6,2L5,G20.10)') 'R3 NREPI,NREPJ,INTFROZEN(NREPI),INTFROZEN(NREPJ),DUMMY=', &
! &                                     NREPI(J2),NREPJ(J2),INTFROZEN(NREPI(J2)),INTFROZEN(NREPJ(J2)),DUMMY
!        ENDIF
         IF (DUMMY.GT.EMAX) THEN
            IMAX=J1
            JMAX=J2
            EMAX=DUMMY
         ENDIF
         IF (DUMMY.GT.EEMAX(J1)) THEN
            JJMAX(J1)=J2
            EEMAX(J1)=DUMMY
         ENDIF
         IF (DUMMY.GT.EEMAX(J1-1)) THEN
            JJMAX(J1-1)=J2
            EEMAX(J1-1)=DUMMY
         ENDIF
         IF (J1.EQ.2) THEN
            EEE(J1)=EEE(J1)+DUMMY
            REPEINT(J1)=REPEINT(J1)+DUMMY
            NREPINT(J1)=NREPINT(J1)+1
         ELSE IF (J1.LT.INTIMAGE+2) THEN
            EEE(J1)=EEE(J1)+DUMMY/2.0D0
            EEE(J1-1)=EEE(J1-1)+DUMMY/2.0D0
            REPEINT(J1)=REPEINT(J1)+DUMMY/2.0D0
            REPEINT(J1-1)=REPEINT(J1-1)+DUMMY/2.0D0
            NREPINT(J1)=NREPINT(J1)+1
            NREPINT(J1-1)=NREPINT(J1-1)+1
         ELSE IF (J1.EQ.INTIMAGE+2) THEN
            EEE(J1-1)=EEE(J1-1)+DUMMY
            REPEINT(J1-1)=REPEINT(J1-1)+DUMMY
            NREPINT(J1-1)=NREPINT(J1-1)+1
         ENDIF
!        DUMMY=-12.0D0*INTCONSTRAINTREP*(1.0D0/(DINT*D12)-1.0D0/INTCONST)
         DUMMY=-2.0D0*INTCONSTRAINTREP*(1.0D0/(DINT*D12)-1.0D0/INTCONST)
         REPGRAD(1:3)=INTMINFAC*DUMMY*G1INT(1:3)
         GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
         GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
         REPGRAD(1:3)=INTMINFAC*DUMMY*G2INT(1:3)
         GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
         GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
   ENDDO
ENDDO
IF (JMAX.GT.0) THEN
   WRITE(*,'(A,I6,A,I6,A,2I6)') ' congrad2> Highest repulsive  contribution for any image in image ',IMAX, &
 &  ' pair index ', &
 &                                JMAX,' atoms ',NREPI(JMAX),NREPJ(JMAX)
ENDIF
CONVERGEREPTEST=EMAX/INTCONSTRAINTREP
!
! Spring energy. Set EEE(J1) and ESPRING dividing up the pairwise
! energy terms between images except for the end points.
!
ESPRING=0.0D0
EMAX=0.0D0
IMAX=0
IF (KINT.NE.0.0D0) THEN
   DO J1=1,INTIMAGE+1 ! sum over edges from J1 to J1+1
      NI1=(3*NATOMS)*(J1-1)
      NI2=(3*NATOMS)*J1
!
!  Edge between J1 and J1+1
!
      DPLUS=0.0D0
      DO J2=1,NATOMS
         IF ((.NOT.INTSPRINGACTIVET).OR.ATOMACTIVE(J2)) THEN 
            DPLUS=DPLUS+(XYZ(NI1+3*(J2-1)+1)-XYZ(NI2+3*(J2-1)+1))**2 &
  &                    +(XYZ(NI1+3*(J2-1)+2)-XYZ(NI2+3*(J2-1)+2))**2 &
  &                    +(XYZ(NI1+3*(J2-1)+3)-XYZ(NI2+3*(J2-1)+3))**2
         ENDIF
!        WRITE(*,'(A,2I8,G20.10)') 'J1,J2,DPLUS: ',J1,J2,DPLUS
      ENDDO
      DPLUS=SQRT(DPLUS)
!     IF (DPLUS.GT.IMSEPMAX) THEN
!        DUMMY=KINT*0.5D0*(DPLUS-IMSEPMAX)**2
         DUMMY=KINT*0.5D0*DPLUS**2
         IF (DUMMY.GT.EMAX) THEN
            IMAX=J1
            EMAX=DUMMY
         ENDIF
!        DUMMY=KINT*0.5D0*DPLUS**2
         ESPRING=ESPRING+DUMMY
         IF (DUMMY.GT.EEMAX(J1+1)) THEN
            EEMAX(J1+1)=DUMMY
         ENDIF

!        WRITE(*,'(A,4G20.10)') 'DPLUS,IMSEPMAX,DUMMY,ESPRING=',DPLUS,IMSEPMAX,DUMMY,ESPRING
!        DUMMY=KINT*(DPLUS-IMSEPMAX)/DPLUS
         DUMMY=KINT
         DO J2=1,NATOMS
            IF ((.NOT.INTSPRINGACTIVET).OR.ATOMACTIVE(J2)) THEN 
               SPGRAD(1:3)=DUMMY*(XYZ(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)-XYZ(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3))
               GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)=GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)+SPGRAD(1:3)
               GGG(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3)=GGG(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3)-SPGRAD(1:3)
            ENDIF
         ENDDO
!     ENDIF
   ENDDO
ENDIF
! IF (KINTENDS.GT.0.0D0) THEN
! !
! ! Extra terms for the two fixed end points.
! !
!    DO J1=2,INTIMAGE+1 ! sum over images
!       NI1=(3*NATOMS)*(J1-1)
!       NI2=(3*NATOMS)*(INTIMAGE+1)
! !
! !  Spring between 1 and J1
! !
!       DPLUS=0.0D0
!       DO J2=1,NATOMS
!          IF ((.NOT.INTSPRINGACTIVET).OR.ATOMACTIVE(J2)) THEN
!             DPLUS=DPLUS+(XYZ(NI1+3*(J2-1)+1)-XYZ(3*(J2-1)+1))**2 &
!   &                    +(XYZ(NI1+3*(J2-1)+2)-XYZ(3*(J2-1)+2))**2 &
!   &                    +(XYZ(NI1+3*(J2-1)+3)-XYZ(3*(J2-1)+3))**2
!          ENDIF
! !        WRITE(*,'(A,2I8,G20.10)') 'J1,J2,DPLUS: ',J1,J2,DPLUS
!       ENDDO
!       DPLUS=SQRT(DPLUS)
! !     IF (DPLUS.GT.IMSEPMAX) THEN
! !        DUMMY=KINTENDS*0.5D0*(DPLUS-IMSEPMAX)**2
!          DUMMY=KINTENDS*0.5D0*DPLUS**2
!          DUMMY=DUMMY*(1.0D0/(J1-1))**2 ! (INTIMAGE/(J1-1))**2
!          IF (DUMMY.GT.EMAX) THEN
!             IMAX=J1
!             EMAX=DUMMY
!          ENDIF
! !        DUMMY=KINTENDS*0.5D0*DPLUS**2
!          ESPRING=ESPRING+DUMMY
!          IF (DUMMY.GT.EEMAX(J1)) THEN
!             EEMAX(J1)=DUMMY
!          ENDIF
! 
! !        WRITE(*,'(A,4G20.10)') 'DPLUS,IMSEPMAX,DUMMY,ESPRING=',DPLUS,IMSEPMAX,DUMMY,ESPRING
! !        DUMMY=KINTENDS*(DPLUS-IMSEPMAX)/DPLUS
!          DUMMY=KINTENDS
!          DO J2=1,NATOMS
!             IF ((.NOT.INTSPRINGACTIVET).OR.ATOMACTIVE(J2)) THEN
!                SPGRAD(1:3)=DUMMY*(XYZ(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)-XYZ(3*(J2-1)+1:3*(J2-1)+3))
!                GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)=GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)+SPGRAD(1:3)
!             ENDIF
!          ENDDO
! !     ENDIF
! !
! !  Spring between INTIMAGE+2 and J1
! !
!       DPLUS=0.0D0
!       DO J2=1,NATOMS
!          IF ((.NOT.INTSPRINGACTIVET).OR.ATOMACTIVE(J2)) THEN
!             DPLUS=DPLUS+(XYZ(NI1+3*(J2-1)+1)-XYZ(NI2+3*(J2-1)+1))**2 &
!   &                    +(XYZ(NI1+3*(J2-1)+2)-XYZ(NI2+3*(J2-1)+2))**2 &
!   &                    +(XYZ(NI1+3*(J2-1)+3)-XYZ(NI2+3*(J2-1)+3))**2
!          ENDIF
! !        WRITE(*,'(A,2I8,G20.10)') 'J1,J2,DPLUS: ',J1,J2,DPLUS
!       ENDDO
!       DPLUS=SQRT(DPLUS)
! !     IF (DPLUS.GT.IMSEPMAX) THEN
! !        DUMMY=KINTENDS*0.5D0*(DPLUS-IMSEPMAX)**2
!          DUMMY=KINTENDS*0.5D0*DPLUS**2
!          DUMMY=DUMMY*(INTIMAGE/(INTIMAGE+2-J1))**2
!          IF (DUMMY.GT.EMAX) THEN
!             IMAX=J1
!             EMAX=DUMMY
!          ENDIF
! !        DUMMY=KINTENDS*0.5D0*DPLUS**2
!          ESPRING=ESPRING+DUMMY
!          IF (DUMMY.GT.EEMAX(J1)) THEN
!             EEMAX(J1)=DUMMY
!          ENDIF
! 
! !        WRITE(*,'(A,4G20.10)') 'DPLUS,IMSEPMAX,DUMMY,ESPRING=',DPLUS,IMSEPMAX,DUMMY,ESPRING
! !        DUMMY=KINTENDS*(DPLUS-IMSEPMAX)/DPLUS
!          DUMMY=KINTENDS
!          DO J2=1,NATOMS
!             IF ((.NOT.INTSPRINGACTIVET).OR.ATOMACTIVE(J2)) THEN
!                SPGRAD(1:3)=DUMMY*(XYZ(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)-XYZ(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3))
!                GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)=GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)+SPGRAD(1:3)
!             ENDIF
!          ENDDO
! !     ENDIF
!    ENDDO
! 
! ENDIF
WRITE(*,'(A,I6,A,I6,A,2I6)') ' congrad2> Highest spring  contribution for any image in image ',IMAX
IF (DEBUG) WRITE(*, '(A,3G20.10)') 'congrad2> ECON,EREP,ESPRING=',ECON,EREP,ESPRING
!
! Set gradients on frozen atoms to zero.
!
IF (FREEZE) THEN
   DO J1=2,INTIMAGE+1  
      DO J2=1,NATOMS
         IF (FROZEN(J2)) THEN
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+1)=0.0D0
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+2)=0.0D0
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+3)=0.0D0
         ENDIF
      ENDDO
   ENDDO
ENDIF
!
! Set gradients on locally frozen atoms to zero.
!
IF (INTFREEZET) THEN
   DO J1=2,INTIMAGE+1
      DO J2=1,NATOMS
         IF (INTFROZEN(J2)) THEN
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+1)=0.0D0
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+2)=0.0D0
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+3)=0.0D0
         ENDIF
      ENDDO
   ENDDO
ENDIF
!
! Set gradients to zero for start and finish images.
!
IF (INTIMAGE.GT.0) THEN
   GGG(1:(3*NATOMS))=0.0D0
   GGG((INTIMAGE+1)*(3*NATOMS)+1:(INTIMAGE+2)*(3*NATOMS))=0.0D0
ENDIF
RMS=0.0D0
RMSIMAGE(1:INTIMAGE+2)=0.0D0
DO J1=2,INTIMAGE+1
   DO J2=1,(3*NATOMS)
      RMSIMAGE(J1)=RMSIMAGE(J1)+GGG((3*NATOMS)*(J1-1)+J2)**2
   ENDDO
   RMS=RMS+RMSIMAGE(J1)
   IF (LPRINT) WRITE(*, '(A,I6,2G20.10)') ' congrad2> J1,EEE,RMSIMAGE=',J1,EEE(J1),RMSIMAGE(J1)
ENDDO
IF (INTIMAGE.NE.0) THEN
   RMS=SQRT(RMS/((3*NATOMS)*INTIMAGE))
ENDIF
!
! For INTIMAGE images there are INTIMAGE+2 replicas including the end points,
! and INTIMAGE+1 line segements, with associated energies stored in EEE(2:INTIMAGE+2)
!

ETOTAL=0.0D0
MAXINT=-1.0D100
MININT=1.0D100
DO J1=2,INTIMAGE+1
   ETOTAL=ETOTAL+EEE(J1)
!  IF (DEBUG) PRINT '(A,I6,A,4G15.5)',' congrad2> con/rep and con/rep int image ', &
! &      J1,' ',CONE(J1),REPE(J1),CONEINT(J1),REPEINT(J1)
   IF (CONEINT(J1)+REPEINT(J1).LT.MININT) THEN
      MININT=CONEINT(J1)+REPEINT(J1)
      NMININT=J1
   ENDIF
   IF (CONEINT(J1)+REPEINT(J1).GT.MAXINT) THEN
      MAXINT=CONEINT(J1)+REPEINT(J1)
      NMAXINT=J1
   ENDIF
ENDDO
! IF (DEBUG) PRINT '(A,G20.10,A,2I6)',' congrad2> largest  internal energy=',MAXINT,' for image ',NMAXINT
! IF (DEBUG) PRINT '(A,G20.10,A,2I6)',' congrad2> smallest internal energy=',MININT,' for image ',NMININT
IF (INTIMAGE.EQ.0) ETOTAL=EEE(1)+EEE(2)

! IF ((RMS.LT.1.0D-50).AND.(ETOTAL.GT.0.1D0).AND.(INTIMAGE.GT.0.0D0)) THEN
!    WRITE(*, '(A,2G20.10)') 'ETOTAL,RMS=',ETOTAL,RMS
!    WRITE(*, '(A,G20.10)') 'ECON=',ECON
!    WRITE(*, '(A,G20.10)') 'EREP=',EREP
!    IF (PRINTE) STOP
!    PRINTE=.TRUE.
!    GOTO 111
! ENDIF

END SUBROUTINE CONGRAD2

SUBROUTINE INTMINONLY(R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,DINT,NOINT)
IMPLICIT NONE
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,DINT,DUMMY
DOUBLE PRECISION DSQI, r1apr2bmr2amr1bsq, r1amr1bsq, r2amr2bsq, r1amr1bdr2amr2b, r1amr1bdr2amr2bsq
LOGICAL NOINT
!
! Is there an internal extremum?
!
! PRINT '(A,4G20.10)','r1ax,r1bx,r2ax,r2bx=',r1ax,r1bx,r2ax,r2bx
! PRINT '(A,G20.10)','(r1ax-r1bx-r2ax+r2bx)**2=',(r1ax-r1bx-r2ax+r2bx)**2
! PRINT '(A,4G20.10)','r1ay,r1by,r2ay,r2by=',r1ay,r1by,r2ay,r2by
! PRINT '(A,G20.10)','(r1ay-r1by-r2ay+r2by)**2=',(r1ay-r1by-r2ay+r2by)**2
! PRINT '(A,4G20.10)','r1az,r1bz,r2az,r2bz=',r1az,r1bz,r2az,r2bz
! PRINT '(A,G20.10)','(r1az-r1bz-r2az+r2bz)**2=',(r1az-r1bz-r2az+r2bz)**2
r1apr2bmr2amr1bsq=(r1ax-r1bx-r2ax+r2bx)**2+(r1ay-r1by-r2ay+r2by)**2+(r1az-r1bz-r2az+r2bz)**2
NOINT=.TRUE.
DINT=1.0D100
IF (r1apr2bmr2amr1bsq.EQ.0.0D0) THEN
   RETURN ! just to skip the internal solution
ELSE
   DUMMY=((r1ax-r1bx)*(r1ax-r1bx-r2ax+r2bx)+ &
 &      (r1ay-r1by)*(r1ay-r1by-r2ay+r2by)+(r1az-r1bz)*(r1az-r1bz-r2az+r2bz))/r1apr2bmr2amr1bsq
ENDIF
IF ((DUMMY.GT.0.0D0).AND.(DUMMY.LT.1.0D0)) NOINT=.FALSE.
IF (.NOT.NOINT) THEN
   r1amr1bdr2amr2b=(r1ax-r1bx)*(r2ax-r2bx)+(r1ay-r1by)*(r2ay-r2by)+(r1az-r1bz)*(r2az-r2bz)
   r1amr1bdr2amr2bsq=r1amr1bdr2amr2b**2
   r1amr1bsq=(r1ax - r1bx)**2 + (r1ay - r1by)**2 + (r1az - r1bz)**2
   r2amr2bsq=(r2ax - r2bx)**2 + (r2ay - r2by)**2 + (r2az - r2bz)**2
   DSQI=MAX((-r1amr1bdr2amr2bsq + r1amr1bsq*r2amr2bsq)/r1apr2bmr2amr1bsq,0.0D0)
   DINT=SQRT(DSQI)
ENDIF

RETURN

END SUBROUTINE INTMINONLY

!
! Call this version for additional repulsive terms between constraints.
!
SUBROUTINE CONGRAD3(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
USE KEY, ONLY: FROZEN, FREEZE, NREPI, NREPJ, NNREPULSIVE, &
  &            NCONSTRAINT, CONI, CONJ, INTCONSTRAINTDEL, CONDISTREF, INTCONSTRAINTREP, CONDISTREFLOCAL, &
  &            CONACTIVE, INTCONSTRAINREPCUT, NREPCUT,INTIMAGE, KINT, IMSEPMAX, ATOMACTIVE, QCINOREPINT, &
  &            INTFREEZET, INTFROZEN, CONCUTLOCAL, CONCUT, CONCUTABST, CONCUTABS, CONCUTFRACT, CONCUTFRAC, &
  &  INTMINFAC, FREEZENODEST, INTSPRINGACTIVET, QCIADDREP, QCIADDREPCUT, QCIADDREPEPS, QCIBONDS
USE COMMONS, ONLY: NATOMS, NOPT, DEBUG
USE PORFUNCS
IMPLICIT NONE
           
INTEGER :: J1,J2,NI2,NI1,NJ2,NJ1,NMAXINT,NMININT,NREPINT(INTIMAGE+2),ISTAT,J3,J4,J5,NINTMIN,NINTMIN2,MYUNIT
DOUBLE PRECISION :: ECON, EREP, ETOTAL, RMS
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,D2,D1
DOUBLE PRECISION G1(3),G2(3),DINT,G1INT(3),G2INT(3)
DOUBLE PRECISION DUMMY, REPGRAD(3), INTCONST, D12, DSQ2, DSQ1, DSQI
DOUBLE PRECISION CONE(INTIMAGE+2), REPE(INTIMAGE+2),MAXINT,MININT,REPEINT(INTIMAGE+2),RMSIM(INTIMAGE+2)
LOGICAL NOINT
DOUBLE PRECISION XYZ((3*NATOMS)*(INTIMAGE+2)), GGG((3*NATOMS)*(INTIMAGE+2)), EEE(INTIMAGE+2), CCLOCAL
LOGICAL IMGFREEZE(INTIMAGE)
DOUBLE PRECISION DPLUS, ESPRING, SPGRAD(3), X1, X2, Y1, Y2, Z1, Z2, DCUT, r1amr1bdr2amr2b,r1apr2bmr2amr1bsq

EEE(1:INTIMAGE+2)=0.0D0
CONE(1:INTIMAGE+2)=0.0D0
REPE(1:INTIMAGE+2)=0.0D0
REPEINT(1:INTIMAGE+2)=0.0D0
NREPINT(1:INTIMAGE+2)=0
GGG(1:(3*NATOMS)*(INTIMAGE+2))=0.0D0
ECON=0.0D0; EREP=0.0D0
MYUNIT=6

IF (QCIADDREP.LT.1) THEN
   WRITE(*,'(A,I6)') 'congrad3> ERROR congrad3 called with no QCIADDREP=',QCIADDREP
   STOP
ENDIF
!
!  Constraint energy and forces.
!
OPEN(UNIT=852,FILE='test.xyz',STATUS='UNKNOWN')
INTCONST=QCIADDREPCUT**3
DO J2=1,NCONSTRAINT
   IF (.NOT.CONACTIVE(J2)) CYCLE
   CCLOCAL=CONCUTLOCAL(J2)
   IF (CONCUTABST) CCLOCAL=CCLOCAL+CONCUTABS
   IF (CONCUTFRACT) CCLOCAL=CCLOCAL+CONCUTFRAC*CONDISTREFLOCAL(J2)
   DO J1=2,INTIMAGE+1
      NI1=(3*NATOMS)*(J1-1)+3*(CONI(J2)-1)
      NJ1=(3*NATOMS)*(J1-1)+3*(CONJ(J2)-1)
      R2AX=XYZ(NI1+1); R2AY=XYZ(NI1+2); R2AZ=XYZ(NI1+3)
      R2BX=XYZ(NJ1+1); R2BY=XYZ(NJ1+2); R2BZ=XYZ(NJ1+3)
      D2=SQRT((R2AX-R2BX)**2+(R2AY-R2BY)**2+(R2AZ-R2BZ)**2)
      IF (ABS(D2-CONDISTREFLOCAL(J2)).GT.CCLOCAL) THEN 
         DUMMY=D2-CONDISTREFLOCAL(J2)  
         G2(1)=(R2AX-R2BX)/D2
         G2(2)=(R2AY-R2BY)/D2
         G2(3)=(R2AZ-R2BZ)/D2
         REPGRAD(1:3)=2*INTCONSTRAINTDEL*((DUMMY/CCLOCAL)**2-1.0D0)*DUMMY*G2(1:3)
         DUMMY=INTCONSTRAINTDEL*(DUMMY**2-CCLOCAL**2)**2/(2.0D0*CCLOCAL**2)
         EEE(J1)=EEE(J1)  +DUMMY
         ECON=ECON        +DUMMY
         CONE(J1)=CONE(J1)+DUMMY
         GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
         GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
      ENDIF
!     WRITE(MYUNIT,'(A,2I6,5G20.10)') 'J1,J2,D2,CONDISTREFLOCAL,CCLOCAL,EEE,CONE=',J1,J2,D2,CONDISTREFLOCAL(J2),CCLOCAL,EEE(J1),CONE(J1)
      IF (J2.GT.QCIBONDS) CYCLE
      DO J3=J2+1,QCIBONDS
         IF (.NOT.CONACTIVE(J3)) CYCLE
         IF (CONI(J3).EQ.CONI(J2)) CYCLE ! no extra terms for bonds with a common atom
         IF (CONI(J3).EQ.CONJ(J2)) CYCLE ! no extra terms for bonds with a common atom
         IF (CONJ(J3).EQ.CONI(J2)) CYCLE ! no extra terms for bonds with a common atom
         IF (CONJ(J3).EQ.CONJ(J2)) CYCLE ! no extra terms for bonds with a common atom
         NI2=(3*NATOMS)*(J1-1)+3*(CONI(J3)-1)
         NJ2=(3*NATOMS)*(J1-1)+3*(CONJ(J3)-1)
         DO J4=1,QCIADDREP
            X1=(J4*XYZ(NI1+1)+(QCIADDREP+1-J4)*XYZ(NJ1+1))/(QCIADDREP+1.0D0)
            Y1=(J4*XYZ(NI1+2)+(QCIADDREP+1-J4)*XYZ(NJ1+2))/(QCIADDREP+1.0D0)
            Z1=(J4*XYZ(NI1+3)+(QCIADDREP+1-J4)*XYZ(NJ1+3))/(QCIADDREP+1.0D0)
            DO J5=1,QCIADDREP
               X2=(J5*XYZ(NI2+1)+(QCIADDREP+1-J5)*XYZ(NJ2+1))/(QCIADDREP+1.0D0)
               Y2=(J5*XYZ(NI2+2)+(QCIADDREP+1-J5)*XYZ(NJ2+2))/(QCIADDREP+1.0D0)
               Z2=(J5*XYZ(NI2+3)+(QCIADDREP+1-J5)*XYZ(NJ2+3))/(QCIADDREP+1.0D0)
               D2=SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
               IF (D2.LT.QCIADDREPCUT) THEN 
                  D12=D2**2
                  DUMMY=QCIADDREPEPS*(1.0D0/D12+(2.0D0*D2-3.0D0*QCIADDREPCUT)/INTCONST)
                  EEE(J1)=EEE(J1)+DUMMY
                  REPE(J1)=REPE(J1)+DUMMY
                  EREP=EREP+DUMMY
!                 WRITE(*,'(A,4I6,A,2I6,A,2G20.10)') 'congrad3> atoms: ',CONI(J2),CONJ(J2),CONI(J3),CONJ(J3), &
! &                     ' sites ',J4,J5,' dist,erep=',D2,DUMMY   
                  DUMMY=-2.0D0*QCIADDREPEPS*(1.0D0/(D2*D12)-1.0D0/INTCONST)
                  G2(1)=(X1-X2)/D2
                  G2(2)=(Y1-Y2)/D2
                  G2(3)=(Z1-Z2)/D2
                  REPGRAD(1:3)=DUMMY*G2(1:3)
                  GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+J4*REPGRAD(1:3)/(QCIADDREP+1.0D0) ! forces on the four atoms involved in the two constraints
                  GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)+(QCIADDREP+1-J4)*REPGRAD(1:3)/(QCIADDREP+1.0D0)
                  GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)-J5*REPGRAD(1:3)/(QCIADDREP+1.0D0)
                  GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-(QCIADDREP+1-J5)*REPGRAD(1:3)/(QCIADDREP+1.0D0)
               ENDIF
            ENDDO
         ENDDO
      ENDDO
   ENDDO
ENDDO
CLOSE(852)

GGG(1:(3*NATOMS))=0.0D0                                      ! can delete when loop range above changes
GGG((3*NATOMS)*(INTIMAGE+1)+1:(3*NATOMS)*(INTIMAGE+2))=0.0D0 ! can delete when loop range above changes

! INTCONST=INTCONSTRAINREPCUT**13

DO J2=1,NNREPULSIVE
!  INTCONST=NREPCUT(J2)**13
   INTCONST=NREPCUT(J2)**3
   DO J1=2,INTIMAGE+2
!  DO J1=1,INTIMAGE+2 ! can change when zero energies are confirmed for end images
      IF (FREEZENODEST) THEN
         IF (J1.EQ.2) THEN
            IF (IMGFREEZE(1)) CYCLE
         ELSE IF (J1.EQ.INTIMAGE+2) THEN
            IF (IMGFREEZE(INTIMAGE)) CYCLE
         ELSE
            IF (IMGFREEZE(J1-2).AND.IMGFREEZE(J1-1)) CYCLE
         ENDIF
      ENDIF
      IF (INTFROZEN(NREPI(J2)).AND.INTFROZEN(NREPJ(J2))) THEN
         WRITE(*, '(A,I6,A,2I6)') ' congrad3> ERROR *** repulsion ',J2,' between frozen atoms ',NREPI(J2),NREPJ(J2)
         STOP
      ENDIF
!     WRITE(*,'(A,2I8,6G20.10)') 'congrad3> B J1,J2,GGG(1:6)=',J1,J2,GGG(1:6)
      NI2=(3*NATOMS)*(J1-1)+3*(NREPI(J2)-1)
      NJ2=(3*NATOMS)*(J1-1)+3*(NREPJ(J2)-1)
      R2AX=XYZ(NI2+1); R2AY=XYZ(NI2+2); R2AZ=XYZ(NI2+3)
      R2BX=XYZ(NJ2+1); R2BY=XYZ(NJ2+2); R2BZ=XYZ(NJ2+3)
      D2=SQRT((R2AX-R2BX)**2+(R2AY-R2BY)**2+(R2AZ-R2BZ)**2)
      IF (D2.LT.NREPCUT(J2)) THEN ! term for image J1
!        D12=D2**12
         D12=D2**2
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*D2-13.0D0*NREPCUT(J2))/INTCONST)
         DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(2.0D0*D2-3.0D0*NREPCUT(J2))/INTCONST)
         EEE(J1)=EEE(J1)+DUMMY
         REPE(J1)=REPE(J1)+DUMMY
         EREP=EREP+DUMMY
!        DUMMY=-12.0D0*INTCONSTRAINTREP*(1.0D0/(D2*D12)-1.0D0/INTCONST)
         DUMMY=-2.0D0*INTCONSTRAINTREP*(1.0D0/(D2*D12)-1.0D0/INTCONST)
         G2(1)=(R2AX-R2BX)/D2
         G2(2)=(R2AY-R2BY)/D2
         G2(3)=(R2AZ-R2BZ)/D2
         REPGRAD(1:3)=DUMMY*G2(1:3)
         GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
         GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
!     WRITE(MYUNIT,'(A,2I6,4G20.10)') 'J1,J2,D2,NREPCUT,EEE,REPE=',J1,J2,D2,NREPCUT(J2),EEE(J1),REPE(J1)
!
! For internal minima we are counting edges. 
! Edge J1 is between images J1-1 and J1, starting from J1=2.
! Energy contributions are shared evenly, except for
! edge 1, which is assigned to image 2, and edge INTIMAGE+1, which
! is assigned to image INTIMAGE+1. Gradients are set to zero for
! the end images.
!
      IF (J1.EQ.1) CYCLE
      IF (QCINOREPINT) CYCLE
      NI1=(3*NATOMS)*(J1-2)+3*(NREPI(J2)-1)
      NJ1=(3*NATOMS)*(J1-2)+3*(NREPJ(J2)-1)

      G1(1:3)=XYZ(NI1+1:NI1+3)-XYZ(NJ1+1:NJ1+3)
      G2(1:3)=XYZ(NI2+1:NI2+3)-XYZ(NJ2+1:NJ2+3)
!
! Squared distance between atoms A and B for theta=0 - distance in image 2
!
      DSQ2=G2(1)**2 + G2(2)**2 + G2(3)**2
!
! Squared distance between atoms A and B for theta=Pi/2 - distance in image 1
!
      DSQ1=G1(1)**2 + G1(2)**2 + G1(3)**2
      DCUT=NREPCUT(J2)**2
      IF ((DSQ1.GT.DCUT).AND.(DSQ2.GT.DCUT)) CYCLE ! don't look for an internal minimum if both repulsions outside cutoff
      r1amr1bdr2amr2b=G1(1)*G2(1)+G1(2)*G2(2)+G1(3)*G2(3)
!
! Is there an internal extremum?
!
      r1apr2bmr2amr1bsq=DSQ1+DSQ2-2.0D0*r1amr1bdr2amr2b

!     IF (r2ax**2+r2ay**2+r2az**2+r2bx**2+r2by**2+r2bz**2-2*(r2ax*r2bx+r2ay*r2by+r2az*r2bz).EQ.0.0D0) THEN
      IF (r1apr2bmr2amr1bsq.LT.1.0D-10) THEN
!        D1=1.0D100; D2=1.0D100;
         NOINT=.TRUE.
         D1=SQRT(DSQ1)
         D2=SQRT(DSQ2)
         G2(1:3)=G2(1:3)/D2
         G1(1:3)=G1(1:3)/D1
      ELSE
         CALL MINMAXD2R(D2,D1,DINT,DSQ2,DSQ1,DSQI,G1,G2,G1INT,G2INT,NOINT,.FALSE.,NREPCUT(J2),r1amr1bdr2amr2b,r1apr2bmr2amr1bsq) 
         IF (.NOT.NOINT) THEN
!           WRITE(*,'(A,I6,A,I6,A,2I6,A,2G20.10)') 'congrad3> internal minimum images ',J1-1,' and ',J1,' atoms: ',NREPI(J2),NREPJ(J2), &
! &                        ' distance,cutoff=',DINT,NREPCUT(J2)
            NINTMIN=NINTMIN+1
         ENDIF
      ENDIF
      IF ((.NOT.NOINT).AND.(DINT.LT.NREPCUT(J2))) THEN
         NINTMIN2=NINTMIN2+1
!        D12=DSQI**6
         D12=DSQI
!        DUMMY=INTCONSTRAINTREP*(1.0D0/D12+(12.0D0*DINT-13.0D0*NREPCUT(J2))/INTCONST)
         DUMMY=INTMINFAC*INTCONSTRAINTREP*(1.0D0/D12+(2.0D0*DINT-3.0D0*NREPCUT(J2))/INTCONST)
         IF (J1.EQ.2) THEN
            EEE(J1)=EEE(J1)+DUMMY
            REPEINT(J1)=REPEINT(J1)+DUMMY
            NREPINT(J1)=NREPINT(J1)+1
         ELSE IF (J1.LT.INTIMAGE+2) THEN
            EEE(J1)=EEE(J1)+DUMMY/2.0D0
            EEE(J1-1)=EEE(J1-1)+DUMMY/2.0D0
            REPEINT(J1)=REPEINT(J1)+DUMMY/2.0D0
            REPEINT(J1-1)=REPEINT(J1-1)+DUMMY/2.0D0
            NREPINT(J1)=NREPINT(J1)+1
            NREPINT(J1-1)=NREPINT(J1-1)+1
         ELSE IF (J1.EQ.INTIMAGE+2) THEN
            EEE(J1-1)=EEE(J1-1)+DUMMY
            REPEINT(J1-1)=REPEINT(J1-1)+DUMMY
            NREPINT(J1-1)=NREPINT(J1-1)+1
         ENDIF
         EREP=EREP+DUMMY
!        DUMMY=-12.0D0*INTCONSTRAINTREP*(1.0D0/(DINT*D12)-1.0D0/INTCONST)
         DUMMY=-2.0D0*INTCONSTRAINTREP*(1.0D0/(DINT*D12)-1.0D0/INTCONST)
         REPGRAD(1:3)=INTMINFAC*DUMMY*G1INT(1:3)
!        PRINT '(A,4I6,2G15.5)','in1 J1,J2,REPI,REPJ,REPGRAD,NREPCUT=',J1,J2,NREPI(J2),NREPJ(J2), &
! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2),NREPCUT(J2)
!
! Gradient contributions for image J1-1
!
         GGG(NI1+1:NI1+3)=GGG(NI1+1:NI1+3)+REPGRAD(1:3)
         GGG(NJ1+1:NJ1+3)=GGG(NJ1+1:NJ1+3)-REPGRAD(1:3)
         REPGRAD(1:3)=INTMINFAC*DUMMY*G2INT(1:3)
!        PRINT '(A,4I6,2G15.5)','in1 J1,J2,REPI,REPJ,REPGRAD,NREPCUT=',J1,J2,NREPI(J2),NREPJ(J2), &
! &                              SQRT(REPGRAD(1)**2+REPGRAD(2)**2+REPGRAD(3)**2),NREPCUT(J2)
!
! Gradient contributions for image J1
!
         GGG(NI2+1:NI2+3)=GGG(NI2+1:NI2+3)+REPGRAD(1:3)
         GGG(NJ2+1:NJ2+3)=GGG(NJ2+1:NJ2+3)-REPGRAD(1:3)
      ENDIF
   ENDDO
ENDDO
!
! Spring energy. Set EEE(J1) and ESPRING dividing up the pairwise
! energy terms between images except for the end points.
!
ESPRING=0.0D0
IF (KINT.NE.0.0D0) THEN
   DO J1=1,INTIMAGE+1 ! sum over edges from J1 to J1+1
      NI1=(3*NATOMS)*(J1-1)
      NI2=(3*NATOMS)*J1
!
!  Edge between J1 and J1+1
!
      DPLUS=0.0D0
      DO J2=1,NATOMS
         IF ((.NOT.INTSPRINGACTIVET).OR.ATOMACTIVE(J2)) THEN 
            DPLUS=DPLUS+(XYZ(NI1+3*(J2-1)+1)-XYZ(NI2+3*(J2-1)+1))**2 &
  &                    +(XYZ(NI1+3*(J2-1)+2)-XYZ(NI2+3*(J2-1)+2))**2 &
  &                    +(XYZ(NI1+3*(J2-1)+3)-XYZ(NI2+3*(J2-1)+3))**2
         ENDIF
      ENDDO
      DPLUS=SQRT(DPLUS)
      IF (DPLUS.GT.IMSEPMAX) THEN
!        DUMMY=KINT*0.5D0*(DPLUS-IMSEPMAX)**2
         DUMMY=KINT*0.5D0*DPLUS**2
         ESPRING=ESPRING+DUMMY
!        DUMMY=KINT*(DPLUS-IMSEPMAX)/DPLUS
         DUMMY=KINT
         DO J2=1,NATOMS
            IF ((.NOT.INTSPRINGACTIVET).OR.ATOMACTIVE(J2)) THEN 
               SPGRAD(1:3)=DUMMY*(XYZ(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)-XYZ(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3))
               GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)=GGG(NI1+3*(J2-1)+1:NI1+3*(J2-1)+3)+SPGRAD(1:3)
               GGG(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3)=GGG(NI2+3*(J2-1)+1:NI2+3*(J2-1)+3)-SPGRAD(1:3)
            ENDIF
         ENDDO
      ENDIF
   ENDDO
ENDIF
!
! Set gradients on frozen atoms to zero.
!
IF (FREEZE) THEN
   DO J1=2,INTIMAGE+1  
      DO J2=1,NATOMS
         IF (FROZEN(J2)) THEN
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+1)=0.0D0
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+2)=0.0D0
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+3)=0.0D0
         ENDIF
      ENDDO
   ENDDO
ENDIF
!
! Set gradients on locally frozen atoms to zero.
!
IF (INTFREEZET) THEN
   DO J1=2,INTIMAGE+1  
      DO J2=1,NATOMS
         IF (INTFROZEN(J2)) THEN
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+1)=0.0D0
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+2)=0.0D0
            GGG((3*NATOMS)*(J1-1)+3*(J2-1)+3)=0.0D0
         ENDIF
      ENDDO
   ENDDO
ENDIF
!
! Set gradients to zero for start and finish images.
!
GGG(1:(3*NATOMS))=0.0D0
GGG((INTIMAGE+1)*(3*NATOMS)+1:(INTIMAGE+2)*(3*NATOMS))=0.0D0
RMS=0.0D0
DO J1=2,INTIMAGE+1
   RMSIM(J1)=0.0D0
   DO J2=1,(3*NATOMS)
      RMS=RMS+GGG((3*NATOMS)*(J1-1)+J2)**2
      RMSIM(J1)=RMSIM(J1)+GGG((3*NATOMS)*(J1-1)+J2)**2
   ENDDO
   RMSIM(J1)=SQRT(RMSIM(J1)/(3*NATOMS))
ENDDO
IF (INTIMAGE.NE.0) THEN
   RMS=SQRT(RMS/((3*NATOMS)*INTIMAGE))
ENDIF
!
! For INTIMAGE images there are INTIMAGE+2 replicas including the end points,
! and INTIMAGE+1 line segements, with associated energies stored in EEE(2:INTIMAGE+2)
!
ETOTAL=0.0D0
MAXINT=-1.0D100
MININT=1.0D100
DO J1=2,INTIMAGE+1
   ETOTAL=ETOTAL+EEE(J1)
!  WRITE(*, '(A,I6,A,3G20.10)') ' congrad3> con/rep/RMS image ',J1,' ',CONE(J1),REPE(J1),RMSIM(J1)
   IF (REPEINT(J1).LT.MININT) THEN
      MININT=REPEINT(J1)
      NMININT=J1
   ENDIF
   IF (REPE(J1).GT.MAXINT) THEN
      MAXINT=REPE(J1)
      NMAXINT=J1
   ENDIF
ENDDO

END SUBROUTINE CONGRAD3
