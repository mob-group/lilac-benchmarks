!   Copyright (C) 2003-2010 David J. Wales
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
SUBROUTINE INTLBFGSLJ(QSTART,QFINISH,NMINFOUND,NTSFOUND,MINFOUND,TSFOUND)
USE PORFUNCS
USE KEY, ONLY : FREEZENODEST, FREEZETOL, MAXINTBFGS, &
     & INTLJTOL, INTIMAGE, INTMUPDATE, INTDGUESS, BFGSTST, NSTEPS, &
     & ATOMACTIVE, INTLJSTEPS, DUMPINTXYZ, DUMPINTXYZFREQ, &
     & MUPDATE, BFGSSTEPS, INTTST, IMSEPMIN, IMSEPMAX
USE COMMONS, ONLY: NATOMS, NOPT, ZSYM, DEBUG
USE MODEFOL

IMPLICIT NONE 

DOUBLE PRECISION, INTENT(IN) :: QSTART(NOPT), QFINISH(NOPT)  ! The two end points
DOUBLE PRECISION EDUMMY,EVALMIN,EVALMAX
INTEGER D, U
DOUBLE PRECISION DF, DMAX, DMIN
INTEGER JMAX, JMIN, INTIMAGESAVE
LOGICAL KNOWE, KNOWG, KNOWH, PTEST, MFLAG, PRINTOPTIMIZETS
COMMON /KNOWN/ KNOWE, KNOWG, KNOWH

INTEGER J2,POINT,BOUND,NPT,CP,I,ISTAT,J1,J3,NIMAGEFREEZE,NACTIVE,NBEST,NMINFOUND,NTSFOUND
INTEGER NSIDE,ITDONE
DOUBLE PRECISION, DIMENSION(3*NATOMS) :: LGDUMMY, VECS, XDIAG
INTEGER NDUMMY
INTEGER NITERDONE, EXITSTATUS, J4
DOUBLE PRECISION :: YS,YY,SQ,YR,BETA,GNORM,DDOT,DUMMY,STPMIN,PREVGRAD,EPLUS,EMINUS, & 
  &                 DINCREMENT, STARTTIME, TIME0, STIME, FTIME, &
  &                 ETOTAL, RMS, STEPTOT, &
  &                 VNEW(NOPT), ENERGY, RMS2, EREAL, LOCALCOORDS(3*NATOMS)
DOUBLE PRECISION,DIMENSION(INTMUPDATE)     :: RHO1,ALPHA
DOUBLE PRECISION :: EOLD, DIFF, DIST, DTOTAL
LOGICAL TSCONVERGED
DOUBLE PRECISION, POINTER :: X(:), G(:)
DOUBLE PRECISION, ALLOCATABLE :: EINT(:)
!
! These declarations have to match those in NEB/ntc.f90
!
TYPE MINFOUNDTYPE
   DOUBLE PRECISION,POINTER :: E
   DOUBLE PRECISION,POINTER :: COORD(:)
END TYPE MINFOUNDTYPE
INTEGER,PARAMETER :: NMINMAX = 3000 ! Maximal number of min to be checked in one intlbfgslj run
TYPE (MINFOUNDTYPE) :: MINFOUND(NMINMAX)

INTEGER,PARAMETER :: NTSMAX = 3000 ! Maximal number of ts to be checked in one intlbfgslj run
TYPE TSFOUNDTYPE
     DOUBLE PRECISION,POINTER :: E
     DOUBLE PRECISION,POINTER :: EVALMIN
     DOUBLE PRECISION,POINTER :: COORD(:)
     DOUBLE PRECISION,POINTER :: VECS(:)
END TYPE TSFOUNDTYPE

TYPE (TSFOUNDTYPE) :: TSFOUND(NTSMAX)
!
! efk: for freezenodes
!
DOUBLE PRECISION :: TESTG, TOTGNORM
INTEGER :: IM
!
! Dimensions involving INTIMAGE
!
INTEGER, PARAMETER :: MAXINTIMAGE=200
DOUBLE PRECISION, ALLOCATABLE :: STEPIMAGE(:), &
  &              GTMP(:), DIAG(:), STP(:), SEARCHSTEP(:,:), GDIF(:,:), GLAST(:)
DOUBLE PRECISION, ALLOCATABLE, TARGET :: XYZ(:), GGG(:), DPTMP(:), D2TMP(:,:)
LOGICAL, ALLOCATABLE :: IMGFREEZE(:)

! LOGICAL EDGEINT(INTIMAGE+1,NATOMS,NATOMS)

ALLOCATE(GTMP(3*NATOMS*INTIMAGE), &
  &      DIAG(3*NATOMS*INTIMAGE), STP(3*NATOMS*INTIMAGE), SEARCHSTEP(0:INTMUPDATE,NOPT*INTIMAGE), &
  &      GDIF(0:INTMUPDATE,NOPT*INTIMAGE), &
  &      XYZ(NOPT*(INTIMAGE+2)), GGG(NOPT*(INTIMAGE+2)), IMGFREEZE(INTIMAGE), &
  &      STEPIMAGE(INTIMAGE))

INTIMAGESAVE=INTIMAGE
CALL MYCPU_TIME(STIME,.FALSE.)
PRINT '(A,I6)',' intlbfgslj> Maximum number of steps for LJ interp phase is ',INTLJSTEPS
PREVGRAD=1.0D100
IF (FREEZENODEST) IMGFREEZE(1:INTIMAGE)=.FALSE.
D=NOPT*INTIMAGE
U=INTMUPDATE
NITERDONE=1

PRINT '(A)',' intlbfgslj> why are we in this routine?'

IF ( D<=0 ) THEN
   PRINT *, 'd is not positive, d=',d
   CALL TSUMMARY
   STOP
ENDIF
IF ( U<=0 ) THEN
   PRINT *, 'u is not positive, u=',u
   CALL TSUMMARY
   STOP
ENDIF
IF (INTLJSTEPS < 0) THEN
   PRINT '(1x,a)', 'Maximal number of iterations is less than zero! Stop.'
   CALL TSUMMARY
   STOP
ENDIF
!
! XYZ, GGG include the end point images
! X, G do not.
!
X=>XYZ(NOPT+1:NOPT*(INTIMAGE+1))
G=>GGG(NOPT+1:NOPT*(INTIMAGE+1))
!
! Initialise XYZ
!
XYZ(1:NOPT)=QSTART(1:NOPT)
XYZ(NOPT*(INTIMAGE+1)+1:NOPT*(INTIMAGE+2))=QFINISH(1:NOPT)
DO J1=1,INTIMAGE+2
   XYZ((J1-1)*NOPT+1:J1*NOPT)=((INTIMAGE+2-J1)*QSTART(1:NOPT)+(J1-1)*QFINISH(1:NOPT))/(INTIMAGE+1)
ENDDO
!
! Calculate initial constraints.
!
IF (.NOT.ALLOCATED(ATOMACTIVE)) ALLOCATE(ATOMACTIVE(NATOMS))
DUMMY=1.0D100

DO J1=1,NATOMS
   DF=SQRT((XYZ(3*(J1-1)+1)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+1))**2 &
  &       +(XYZ(3*(J1-1)+2)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+2))**2 &
  &       +(XYZ(3*(J1-1)+3)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+3))**2)
   IF (DF.LT.DUMMY) THEN
      NBEST=J1
      DUMMY=DF
   ENDIF
ENDDO
IF (DEBUG) PRINT '(A,I6,A,F15.5)',' intlbfgslj> Smallest overall motion for atom ',NBEST,' dist=',DUMMY
! EDGEINT(1:INTIMAGE+1,1:NATOMS,1:NATOMS)=.FALSE.
!!!!!!!!!!!!!!!!!!! all in one from linear !!!!!!!!!!!!!!!!!
NACTIVE=NATOMS
ATOMACTIVE(1:NATOMS)=.TRUE.
CALL INTGRADLJ(ETOTAL,XYZ,GGG,IMGFREEZE,RMS,.FALSE.)
EOLD=ETOTAL
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

IF (DEBUG) WRITE(*,'(A26,A20,A15,A13,A9)') 'Iter','Energy per image','RMS Force','Step'

DO ! Main do loop with counter NITERDONE, initially set to one

   MAIN: IF (NITERDONE==1) THEN
        POINT = 0
        DIAG(1:D)=INTDGUESS
        SEARCHSTEP(0,:)= -G*DIAG                           ! NR STEP FOR DIAGONAL INVERSE HESSIAN
        GTMP           = SEARCHSTEP(0,:)
        GNORM          = MAX(SQRT(DOT_PRODUCT(G,G)),1.0D-100)
        STP            = MIN(1.0D0/GNORM, GNORM)           ! MAKE THE FIRST GUESS FOR THE STEP LENGTH CAUTIOUS
   ELSE MAIN
        BOUND=NITERDONE-1
        IF (NITERDONE.GT.INTMUPDATE) BOUND=INTMUPDATE
        YS=DOT_PRODUCT( GDIF(NPT/D,:), SEARCHSTEP(NPT/D,:)  )
        IF (YS==0.0D0) YS=1.0D0
    
        ! Update estimate of diagonal inverse Hessian elements.
        ! We divide by both YS and YY at different points, so they had better not be zero!

        YY=DOT_PRODUCT( GDIF(NPT/D,:) , GDIF(NPT/D,:) )
        IF (YY==0.0D0) YY=1.0D0
!       DIAG = ABS(YS/YY)
        DIAG = YS/YY
      
        ! COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980, "Updating quasi-Newton matrices with limited storage",
        ! Mathematics of Computation, Vol.35, No.151, pp. 773-782
        CP= POINT; IF (POINT==0) CP = INTMUPDATE
        RHO1(CP)=1.0D0/YS
        GTMP = -G
        CP= POINT 
                   
        DO I= 1,BOUND 
!            CP = CP - 1; IF (CP == -1) CP = M - 1
             CP = CP - 1; IF (CP == -1) CP = INTMUPDATE - 1
             SQ= DOT_PRODUCT( SEARCHSTEP(CP,:),GTMP )
             ALPHA(CP+1) = RHO1(CP+1) * SQ
             GTMP        = -ALPHA(CP+1)*GDIF(CP,:) + GTMP
        ENDDO
              
        GTMP=DIAG*GTMP

        DO I=1,BOUND
             YR= DOT_PRODUCT( GDIF(CP,:) , GTMP )
             BETA= RHO1(CP+1)*YR
             BETA= ALPHA(CP+1)-BETA
             GTMP = BETA*SEARCHSTEP(CP,:) + GTMP
             CP=CP+1
!            IF (CP==M) CP=0
             IF (CP==INTMUPDATE) CP=0
        ENDDO
              
        STP(1:D) = 1.0D0
   ENDIF MAIN

   !  Store the new search direction
   IF (NITERDONE.GT.1) SEARCHSTEP(POINT,:)=GTMP
      
!
! If the number of images has changed since G was declared then G is not the same
! size as Gtmp and Dot_Product cannot be used.
!
!  IF (Dot_Product(G,Gtmp)/SQRT( Dot_Product(G,G)*Dot_Product(Gtmp,Gtmp) ) > 0.0D0) THEN
!
!  Separate sqrt;s to avoid overflow.
!
   IF (DDOT(D,G,1,GTMP,1)/MAX(1.0D-100,SQRT( DDOT(D,G,1,G,1))*SQRT(DDOT(D,GTMP,1,GTMP,1)) ) > 0.0D0) THEN
        IF (DEBUG) PRINT*,'Search direction has positive projection onto gradient - reversing step'
        GTMP=-GTMP
        SEARCHSTEP(POINT,:)=GTMP
   ENDIF
   GTMP=G

!  We should apply the maximum LBFGS step to each image separately.
!  However, using different scale factors for different images leads to huge
!  discontinuities! Now take the minimum scale factor for all images. DJW 26/11/07

   STPMIN=1.0D0
   DO J2=1,INTIMAGE
      STEPIMAGE(J2) = SQRT(DOT_PRODUCT(SEARCHSTEP(POINT,NOPT*(J2-1)+1:NOPT*J2),SEARCHSTEP(POINT,NOPT*(J2-1)+1:NOPT*J2)))
      DUMMY=STEPIMAGE(J2)
      IF (STEPIMAGE(J2) > MAXINTBFGS) THEN
           STP(NOPT*(J2-1)+1:NOPT*J2) = MAXINTBFGS/STEPIMAGE(J2)
           STPMIN=MIN(STPMIN,STP(NOPT*(J2-1)+1))
      ENDIF
!     PRINT '(A,I8,3G20.10)','image,initial step size,STP,prod=',J2,DUMMY,STP(NOPT*(J2-1)+1),STEPIMAGE(J2)*STP(NOPT*(J2-1)+1)
   ENDDO
   STP(1:D)=STPMIN

! EFK: decide whether to freeze some nodes
   IF (FREEZENODEST) THEN
      TOTGNORM=SQRT(DOT_PRODUCT(G(1:NOPT*INTIMAGE),G(1:NOPT*INTIMAGE))/INTIMAGE)
      NIMAGEFREEZE=0
      DO IM=1,INTIMAGE
         TESTG=SQRT(DOT_PRODUCT(G(NOPT*(IM-1)+1:NOPT*IM),G(NOPT*(IM-1)+1:NOPT*IM)))
         IMGFREEZE(IM)=.FALSE.
         IF (TOTGNORM.NE.0.0D0) THEN
            IF (TESTG/TOTGNORM.LT.FREEZETOL) THEN
!              IF (DEBUG) PRINT '(A,I6,2G20.10)', ' intlbfgslj> Freezing image: ', IM, TESTG, TOTGNORM
               IMGFREEZE(IM)=.TRUE.
               STEPIMAGE(IM)=0.0D0
               NIMAGEFREEZE=NIMAGEFREEZE+1
               STP(NOPT*(IM-1)+1:NOPT*IM)=0.0D0
            ENDIF
         ENDIF
      ENDDO
      IF (DEBUG) PRINT '(2(A,I6))', ' intlbfgslj> Number of frozen images=',NIMAGEFREEZE,' / ',INTIMAGE
   ENDIF
   !  We now have the proposed step - update geometry and calculate new gradient
20 X(1:D) = X(1:D) + STP(1:D)*SEARCHSTEP(POINT,1:D)

!  IF ((RMS.LT.INTRMSTOL*1.0D10).AND.(MOD(NITERDONE,10).EQ.0).AND.(NSTEPSMAX-NITERDONE.GT.100)) &
! &            CALL CHECKSEP(NMAXINT,NMININT,INTIMAGE,XYZ,NOPT,NATOMS)
   IF (MOD(NITERDONE,50).EQ.0) THEN
      DMAX=0.0D0
      DMIN=HUGE(1.0D0)
      DO J1=1,INTIMAGE+1
         DUMMY=0.0D0
         DO J2=1,3*NATOMS
            IF (ATOMACTIVE((J2-1)/3+1)) THEN
                  DUMMY=DUMMY+( XYZ((J1-1)*3*NATOMS+J2) - XYZ(J1*3*NATOMS+J2) )**2
         ENDIF
         ENDDO
         DUMMY=SQRT(DUMMY)
         IF (DUMMY.GT.DMAX) THEN
            DMAX=DUMMY
            JMAX=J1
         ENDIF
         IF (DUMMY.LT.DMIN) THEN
            DMIN=DUMMY
            JMIN=J1
         ENDIF
         IF (DEBUG) PRINT '(A,I6,A,I6,A,G20.10)',' intlbfgs> distance between images ', &
  &                                               J1,' and ',J1+1,' is ',DUMMY
      ENDDO
      IF ((DMAX.GT.IMSEPMAX).AND.(INTIMAGE.LT.MAXINTIMAGE)) THEN
         PRINT '(A,I6,A,I6)',' intlbfgs> Add an image between ',JMAX,' and ',JMAX+1
         ALLOCATE(DPTMP(3*NATOMS*(INTIMAGE+2)))
         DPTMP(1:3*NATOMS*(INTIMAGE+2))=XYZ(1:3*NATOMS*(INTIMAGE+2))
         DEALLOCATE(XYZ)
         ALLOCATE(XYZ(3*NATOMS*(INTIMAGE+3)))
         XYZ(1:3*NATOMS*JMAX)=DPTMP(1:3*NATOMS*JMAX)
         XYZ(3*NATOMS*JMAX+1:3*NATOMS*(JMAX+1))=(DPTMP(3*NATOMS*(JMAX-1)+1:3*NATOMS*JMAX) &
  &                                            + DPTMP(3*NATOMS*JMAX+1:3*NATOMS*(JMAX+1)))/2.0D0
         XYZ(3*NATOMS*(JMAX+1)+1:3*NATOMS*(INTIMAGE+3))=DPTMP(3*NATOMS*JMAX+1:3*NATOMS*(INTIMAGE+2))
!
! Save step-taking memories in SEARCHSTEP and GDIF.
! These arrays run from 0 to INTMUPDATE over memories and
! 1:NOPT*INTIMAGE over only the variable images.
!
         DEALLOCATE(DPTMP)
         ALLOCATE(D2TMP(0:INTMUPDATE,1:NOPT*INTIMAGE))
         D2TMP(0:INTMUPDATE,1:NOPT*INTIMAGE)=SEARCHSTEP(0:INTMUPDATE,1:NOPT*INTIMAGE)
         DEALLOCATE(SEARCHSTEP)
         ALLOCATE(SEARCHSTEP(0:INTMUPDATE,1:NOPT*(INTIMAGE+1)))
         DO J1=1,INTMUPDATE
            IF (JMAX.GT.1) SEARCHSTEP(J1,1:3*NATOMS*(JMAX-1))=D2TMP(J1,1:3*NATOMS*(JMAX-1))
            IF (JMAX.LT.INTIMAGE+1) SEARCHSTEP(J1,3*NATOMS*JMAX+1:3*NATOMS*(INTIMAGE+1))= &
  &              D2TMP(J1,3*NATOMS*(JMAX-1)+1:3*NATOMS*INTIMAGE)
            SEARCHSTEP(J1,3*NATOMS*(JMAX-1)+1:3*NATOMS*JMAX)= &
  &                          D2TMP(J1,3*NATOMS*(MIN(JMAX,INTIMAGE)-1)+1:3*NATOMS*MIN(JMAX,INTIMAGE))
         ENDDO
         D2TMP(0:INTMUPDATE,1:NOPT*INTIMAGE)=GDIF(0:INTMUPDATE,1:NOPT*INTIMAGE)
         DEALLOCATE(GDIF)
         ALLOCATE(GDIF(0:INTMUPDATE,1:NOPT*(INTIMAGE+1)))
         DO J1=1,INTMUPDATE
            IF (JMAX.GT.1) GDIF(J1,1:3*NATOMS*(JMAX-1))=D2TMP(J1,1:3*NATOMS*(JMAX-1))
            IF (JMAX.LT.INTIMAGE+1) GDIF(J1,3*NATOMS*JMAX+1:3*NATOMS*(INTIMAGE+1))= &
  &              D2TMP(J1,3*NATOMS*(JMAX-1)+1:3*NATOMS*INTIMAGE)
            GDIF(J1,3*NATOMS*(JMAX-1)+1:3*NATOMS*JMAX)= &
  &                    D2TMP(J1,3*NATOMS*(MIN(JMAX,INTIMAGE)-1)+1:3*NATOMS*MIN(JMAX,INTIMAGE))
         ENDDO
         DEALLOCATE(D2TMP)

         DEALLOCATE(GTMP,GGG,DIAG,STP,STEPIMAGE,IMGFREEZE)
         ALLOCATE(GTMP(3*NATOMS*(INTIMAGE+1)), &
  &               DIAG(3*NATOMS*(INTIMAGE+1)), STP(3*NATOMS*(INTIMAGE+1)), &
  &               IMGFREEZE(INTIMAGE+1), &
  &               STEPIMAGE(INTIMAGE+1), GGG(3*NATOMS*(INTIMAGE+3)))
         GGG(1:3*NATOMS*(INTIMAGE+3))=0.0D0
         GTMP(1:3*NATOMS*(INTIMAGE+1))=0.0D0
         DIAG(1:3*NATOMS*(INTIMAGE+1))=0.0D0
         STP(1:3*NATOMS*(INTIMAGE+1))=0.0D0
         IMGFREEZE(1:INTIMAGE+1)=.FALSE.
         STEPIMAGE(1:INTIMAGE+1)=0.0D0

         X=>XYZ(NOPT+1:NOPT*(INTIMAGE+2))
         G=>GGG(NOPT+1:NOPT*(INTIMAGE+2))
         INTIMAGE=INTIMAGE+1
         D=NOPT*INTIMAGE
      ELSEIF ((DMIN.LT.IMSEPMIN).AND.(INTIMAGE.GT.1)) THEN
         IF (JMIN.EQ.1) JMIN=2
         PRINT '(A,I6,A,I6)',' intlbfgs> Remove image ',JMIN
         ALLOCATE(DPTMP(3*NATOMS*(INTIMAGE+2)))
         DPTMP(1:3*NATOMS*(INTIMAGE+2))=XYZ(1:3*NATOMS*(INTIMAGE+2))
         DEALLOCATE(XYZ)
         ALLOCATE(XYZ(3*NATOMS*(INTIMAGE+1)))
         XYZ(1:3*NATOMS*(JMIN-1))=DPTMP(1:3*NATOMS*(JMIN-1))
         XYZ(3*NATOMS*(JMIN-1)+1:3*NATOMS*(INTIMAGE+1))=DPTMP(3*NATOMS*JMIN+1:3*NATOMS*(INTIMAGE+2))

         DEALLOCATE(DPTMP)
!
! Save step-taking memories in SEARCHSTEP and GDIF.
! These arrays run from 0 to INTMUPDATE over memories and
! 1:NOPT*INTIMAGE over only the variable images.
!
         ALLOCATE(D2TMP(0:INTMUPDATE,1:NOPT*INTIMAGE))
         D2TMP(0:INTMUPDATE,1:NOPT*INTIMAGE)=SEARCHSTEP(0:INTMUPDATE,1:NOPT*INTIMAGE)
         DEALLOCATE(SEARCHSTEP)
         ALLOCATE(SEARCHSTEP(0:INTMUPDATE,1:NOPT*(INTIMAGE-1)))
         DO J1=1,INTMUPDATE
            SEARCHSTEP(J1,1:3*NATOMS*(JMIN-2))=D2TMP(J1,1:3*NATOMS*(JMIN-2))
            SEARCHSTEP(J1,3*NATOMS*(JMIN-2)+1:3*NATOMS*(INTIMAGE-1))= &
  &                  D2TMP(J1,3*NATOMS*(JMIN-1)+1:3*NATOMS*INTIMAGE)
         ENDDO
         D2TMP(0:INTMUPDATE,1:NOPT*INTIMAGE)=GDIF(0:INTMUPDATE,1:NOPT*INTIMAGE)
         DEALLOCATE(GDIF)
         ALLOCATE(GDIF(0:INTMUPDATE,1:NOPT*(INTIMAGE-1)))
         DO J1=1,INTMUPDATE
            GDIF(J1,1:3*NATOMS*(JMIN-2))=D2TMP(J1,1:3*NATOMS*(JMIN-2))
            GDIF(J1,3*NATOMS*(JMIN-2)+1:3*NATOMS*(INTIMAGE-1))= &
  &                  D2TMP(J1,3*NATOMS*(JMIN-1)+1:3*NATOMS*INTIMAGE)
         ENDDO
         DEALLOCATE(D2TMP)

         DEALLOCATE(GTMP,GGG,DIAG,STP,STEPIMAGE,IMGFREEZE)
         ALLOCATE(GTMP(3*NATOMS*(INTIMAGE-1)), &
  &               DIAG(3*NATOMS*(INTIMAGE-1)), STP(3*NATOMS*(INTIMAGE-1)), &
  &               IMGFREEZE(INTIMAGE-1), &
  &               STEPIMAGE(INTIMAGE-1), GGG(3*NATOMS*(INTIMAGE+1)))
         GGG(1:3*NATOMS*(INTIMAGE+1))=0.0D0
         GTMP(1:3*NATOMS*(INTIMAGE-1))=0.0D0
         DIAG(1:3*NATOMS*(INTIMAGE-1))=0.0D0
         STP(1:3*NATOMS*(INTIMAGE-1))=0.0D0
         IMGFREEZE(1:INTIMAGE-1)=.FALSE.
         STEPIMAGE(1:INTIMAGE-1)=0.0D0

         X=>XYZ(NOPT+1:NOPT*(INTIMAGE))
         G=>GGG(NOPT+1:NOPT*(INTIMAGE))
         INTIMAGE=INTIMAGE-1
         D=NOPT*INTIMAGE
      ENDIF
   ENDIF
   CALL INTGRADLJ(ETOTAL,XYZ,GGG,IMGFREEZE,RMS,.TRUE.)
   EOLD=ETOTAL
   STEPTOT = SUM(STEPIMAGE)/INTIMAGE
   IF (DEBUG) THEN
      WRITE(*,'(A,I6,2G20.10,F9.3)') ' intlbfgslj> steps: ',NITERDONE,ETOTAL/INTIMAGE,RMS,STEPTOT
      CALL FLUSH(6)
   ENDIF

   EXITSTATUS=0
   INTDGUESS=DIAG(1) ! should be ok for subsequent runs of the same system DJW
   IF (RMS<=INTLJTOL.AND.NITERDONE>1) EXITSTATUS=1
   IF (NITERDONE==INTLJSTEPS) EXITSTATUS=2

   IF (EXITSTATUS > 0) THEN  
      IF (EXITSTATUS.EQ.1) THEN ! add active atom 
         IF (NACTIVE.LT.NATOMS) THEN 
            GOTO 777
         ENDIF
!        CALL MYCPU_TIME(FTIME,.FALSE.)
!        PRINT '(A,I6,A,F10.1)',' intlbfgslj> converged at step ',NITERDONE,' time=',FTIME-STIME
      ELSEIF (EXITSTATUS.EQ.2) THEN 
         CALL MYCPU_TIME(FTIME,.FALSE.)
         PRINT '(A,F10.1)',' intlbfgslj> Failed to achieve requested RMS convergence, time=',FTIME-STIME
      ENDIF
      EXIT
   ENDIF
   777 CONTINUE
!
! Compute the new step and gradient change
!
   NPT=POINT*D
   SEARCHSTEP(POINT,:) = STP*SEARCHSTEP(POINT,:)
   GDIF(POINT,:)=G-GTMP
   POINT=POINT+1; IF (POINT==INTMUPDATE) POINT=0

   IF (DUMPINTXYZ.AND.MOD(NITERDONE,DUMPINTXYZFREQ)==0) CALL INTRWG(NACTIVE,NITERDONE,INTIMAGE,XYZ)
   PREVGRAD=RMS

   NITERDONE=NITERDONE+1
   IF (NITERDONE.GT.INTLJSTEPS) EXIT

ENDDO ! end of main do loop over counter NITERDONE

IF (EXITSTATUS.EQ.1) THEN
   CALL MYCPU_TIME(FTIME,.FALSE.)
   WRITE(*,'(A,I6,A,G20.10,A,G15.5,A,F10.1)') ' intlbfgslj> Converged in ',NITERDONE,' steps, energy/image=',ETOTAL/INTIMAGE, &
  &                               ' RMS=',RMS,' time=',FTIME-STIME
ELSEIF (EXITSTATUS.EQ.2) THEN
   WRITE(*,'(A,I6,A,G20.10,A,G20.10)') ' intlbfgslj> After ',NITERDONE,' steps, energy per image=',ETOTAL/INTIMAGE, &
  &                               ' RMS gradient=',RMS
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!        DIFF=1.0D-5
!        CALL INTGRADLJ(ETOTAL,XYZ,GGG,IMGFREEZE,RMS,.TRUE.)
!        EOLD=ETOTAL
!        ALLOCATE(GLAST(D))
!        GLAST(1:D)=G(1:D)
!        PRINT '(A,I6)',' intlbfgslj> analytic and numerical gradients: D=',D
!        DO J2=1,D
!           X(J2)=X(J2)+DIFF
!           CALL INTGRADLJ(ETOTAL,XYZ,GGG,IMGFREEZE,RMS,.TRUE.)
!           EPLUS=ETOTAL
!           X(J2)=X(J2)-2.0D0*DIFF
!           CALL INTGRADLJ(ETOTAL,XYZ,GGG,IMGFREEZE,RMS,.TRUE.)
!           EMINUS=ETOTAL
!           X(J2)=X(J2)+DIFF
!           IF (ABS(GLAST(J2)).NE.0.0D0) THEN
!              IF (100.0D0*ABS((GLAST(J2)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GLAST(J2)).GT.10.0D0) THEN
!                 WRITE(*,'(A,3I8,3G20.10)') 'error ',(J2-1)/NOPT+1,(J2-NOPT*((J2-1)/NOPT)-1)/3+1,J2, &
!    &                                 GLAST(J2),(EPLUS-EMINUS)/(2.0D0*DIFF), &
!    &                                 (EPLUS-EMINUS)/(2.0D0*DIFF*GLAST(J2))
!              ELSE
!                 WRITE(*,'(A,3I8,3G20.10)') 'OK    ',(J2-1)/NOPT+1,(J2-NOPT*((J2-1)/NOPT)-1)/3+1,J2, &
!    &                                       GLAST(J2),(EPLUS-EMINUS)/(2.0D0*DIFF), &
!    &                                       (EPLUS-EMINUS)/(2.0D0*DIFF*GLAST(J2))
!              ENDIF
!           ENDIF
!        ENDDO
!        DEALLOCATE(GLAST)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Linear interpolation for real potential
!
DINCREMENT=0.01D0
DTOTAL=0.0D0
OPEN(UNIT=753,FILE='intenergy',STATUS='UNKNOWN')
!
! local maxima must have NSIDE higher energies on each side
! This has the desirable side-effect that we don't bother with 
! images that are essentially collapsed on each other - their
! spacing will probably be < DINCREMENT, or 5*DINCREMENT.
!
NSIDE=5      
INTTST=.TRUE. ! try passing local maxima back as ts guesses
NTSFOUND=0
NMINFOUND=0
PRINTOPTIMIZETS=DEBUG
DO J1=1,INTIMAGE+1
   DUMMY=0.0D0
   DO J2=1,3*NATOMS
      DUMMY=DUMMY+( XYZ((J1-1)*3*NATOMS+J2) - XYZ(J1*3*NATOMS+J2) )**2
   ENDDO
   DUMMY=SQRT(DUMMY)
   DIST=0.0D0
   PRINT '(A,I6,A,I6,A,G20.10)',' intlbfgslj> distance between images ',J1,' and ',J1+1,' is ',DUMMY
   NDUMMY=DUMMY/DINCREMENT+1
   ALLOCATE(EINT(NDUMMY))
   J3=1
 
   intloop: DO 
      LOCALCOORDS(1:3*NATOMS)=((DUMMY-DIST)*XYZ((J1-1)*3*NATOMS+1:J1*3*NATOMS)+ &
  &                                    DIST*XYZ(J1*3*NATOMS+1:(J1+1)*3*NATOMS))/DUMMY
      CALL POTENTIAL(LOCALCOORDS,EREAL,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
      If (DEBUG) PRINT '(A,3G20.10)',' intlbfgslj> ',DTOTAL+DIST,EREAL
      WRITE(753,'(3G20.10)') DTOTAL+DIST,EREAL
      DIST=DIST+DINCREMENT
      EINT(J3)=EREAL
      IF (INTTST) THEN
         IF (J3-NSIDE.GT.0) THEN
            DO J4=MAX(J3-2*NSIDE,1),J3
               IF (J4.EQ.J3-NSIDE) CYCLE
               IF (EINT(J3-NSIDE).LT.EINT(J4)) GOTO 432
            ENDDO
!
! We have a ts candidate. Try optimising it!
!
            CALL MYCPU_TIME(STARTTIME,.FALSE.)
            KNOWG=.FALSE.
            KNOWE=.FALSE. ! to be safe!
            LOCALCOORDS(1:NOPT)= &
  &                 ((DUMMY-(J3-NSIDE-1)*DINCREMENT)*XYZ((J1-1)*NOPT+1:J1*NOPT)+ &
  &                         (J3-NSIDE-1)*DINCREMENT *XYZ(J1*NOPT+1:(J1+1)*NOPT))/DUMMY
            IF (BFGSTST) THEN
               VECS(1:NOPT)=(XYZ((J1-1)*NOPT+1:J1*NOPT)-XYZ(J1*NOPT+1:(J1+1)*NOPT))/DUMMY
               CALL BFGSTS(NSTEPS,LOCALCOORDS,  &
  &               EDUMMY,LGDUMMY,TSCONVERGED,RMS,EVALMIN,EVALMAX,VECS,ITDONE,.TRUE.,PRINTOPTIMIZETS)
            ELSE
               CALL EFOL(LOCALCOORDS,TSCONVERGED,NSTEPS,EDUMMY,ITDONE,EVALMIN,DEBUG,XDIAG,2)
            ENDIF
            CALL MYCPU_TIME(TIME0,.FALSE.)
            IF (TSCONVERGED) THEN
               NTSFOUND=NTSFOUND+1
!
! Save coordinates and direction vector between images to use as starting guess
! for the eigenvector.
!
               ALLOCATE(TSFOUND(NTSFOUND)%E,TSFOUND(NTSFOUND)%COORD(NOPT), &
  &                     TSFOUND(NTSFOUND)%EVALMIN,TSFOUND(NTSFOUND)%VECS(NOPT))
               TSFOUND(NTSFOUND)%VECS(1:NOPT)=VECS(1:NOPT)
               TSFOUND(NTSFOUND)%COORD(1:NOPT)=LOCALCOORDS(1:NOPT)
               TSFOUND(NTSFOUND)%E=EDUMMY
               TSFOUND(NTSFOUND)%EVALMIN=EVALMIN
               PRINT '(A,I6,A,G20.10,A,F10.1)',' intlbfgslj> transition state found, iterations=',ITDONE, &
  &                                  ' energy=',EDUMMY,' time=',TIME0-STARTTIME
            ENDIF
432         CONTINUE
         ENDIF
      ENDIF
      J3=J3+1
      IF (DIST.GT.DUMMY) EXIT INTLOOP
      IF (J3.GT.NDUMMY) THEN
         PRINT '(A,I6)',' intlbfgslj> ERROR *** number of interpolated energies should not be ',J3
      ENDIF
   ENDDO intloop
   DTOTAL=DTOTAL+DUMMY
   DEALLOCATE(EINT)
ENDDO

LOCALCOORDS(1:3*NATOMS)=XYZ((INTIMAGE+1)*3*NATOMS+1:(INTIMAGE+2)*3*NATOMS)
CALL POTENTIAL(LOCALCOORDS,EREAL,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
IF (DEBUG) PRINT '(A,3G20.10)',' intlbfgslj> ',DTOTAL,EREAL
WRITE(753,'(3G20.10)') DTOTAL,EREAL
CLOSE(753)

IF (.NOT.INTTST) THEN
   PTEST=.FALSE.
   PRINT '(A,I8)',' intlbfgslj> minimising all the images - results written to images.min'
   OPEN(987,FILE='images.min',STATUS='UNKNOWN')
   WRITE(987,'(I6)') NATOMS
   WRITE(987,'(A)') 'start - image 1'
   WRITE(987,'(A,3G20.10)') (ZSYM(J2),XYZ(3*(J2-1)+1:3*(J2-1)+3),J2=1,NATOMS)
   DO J1=2,INTIMAGE+1
      KNOWG=.FALSE.
      KNOWE=.FALSE. 
!     PRINT '(A,I8,A,F20.10)',' intlbfgslj> minimising image ',J1

!     BSMIN=.TRUE.
!     DEBUG=.TRUE.
!     CALL ODESD(100,XYZ(NOPT*(J1-1)+1:NOPT*J1),MFLAG,ITDONE,.TRUE.)
!     DEBUG=.FALSE.
!     BSMIN=.FALSE.

!     KNOWG=.FALSE.
!     KNOWE=.FALSE. 
!     PTEST=.TRUE.
      CALL MYLBFGS(NOPT,MUPDATE,XYZ(NOPT*(J1-1)+1:NOPT*J1),.FALSE., &
   &               MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,.TRUE.,ITDONE,PTEST,VNEW,.TRUE.,.FALSE.)
!     PTEST=.FALSE.

      IF (MFLAG) THEN
         NMINFOUND=NMINFOUND+1
!
!  We have to communicate the minima found back to tryconnect using the data structure
!  set up for new transition states.
!  Added new variable MINFOUND to allow for this check in tryconnect.
!  It seems impossible to make intlbfgs see isnewmin and addnewmin for some reason.
!
         ALLOCATE(MINFOUND(NMINFOUND)%E,MINFOUND(NMINFOUND)%COORD(NOPT))
         MINFOUND(NMINFOUND)%COORD(1:NOPT)=XYZ(NOPT*(J1-1)+1:NOPT*J1)
         MINFOUND(NMINFOUND)%E=EREAL
         WRITE(987,'(I6)') NATOMS
         WRITE(987,'(A,I5)') 'image ',J1
         WRITE(987,'(A,3G20.10)') (ZSYM(J2), MINFOUND(NMINFOUND)%COORD(3*(J2-1)+1:3*(J2-1)+3),J2=1,NATOMS)
      ENDIF
   ENDDO
   WRITE(987,'(I6)') NATOMS
   WRITE(987,'(A)') 'finish - image INTIMAGE+2'
   WRITE(987,'(A,3G20.10)') (ZSYM(J2),XYZ(NOPT*(INTIMAGE+1)+3*(J2-1)+1:NOPT*(INTIMAGE+1)+3*(J2-1)+3),J2=1,NATOMS)
   CLOSE(987)
ENDIF

DEALLOCATE(GTMP, DIAG, STP, SEARCHSTEP, GDIF, XYZ, GGG, IMGFREEZE, STEPIMAGE)
INTIMAGE=INTIMAGESAVE

END SUBROUTINE INTLBFGSLJ
