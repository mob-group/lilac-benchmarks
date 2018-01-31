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
MODULE MINIMISER3
     IMPLICIT NONE
     CONTAINS

SUBROUTINE NEBBFGSINT(D,U)
     USE KEYMINIMIZER
     USE KEYLBFGS
     USE KEYNEB,ONLY: NIMAGE,MOREPRINTING,DUMPNEBXYZ,DUMPNEBXYZFREQ,DUMPNEBEOS,DUMPNEBEOSFREQ,DEBUG,&
     &DUMPNEBPTS,DUMPNEBPTSFREQ
     USE NEBDATA
     USE NEBUTILS
     USE GRADIENTS
     USE NEBOUTPUT
     USE MODUNRES
     USE PORFUNCS
     USE KEY,ONLY : NEBMUPDATE, MAXNEBBFGS, NEBDGUESS
     IMPLICIT NONE 

     INTEGER,INTENT(IN) :: D,U ! DIMENSIONALITY OF THE PROBLEM (NINTS*NIMAGE) AND NUMBER OF UPDATES
     INTEGER :: J1,J2,POINT,BOUND,NPT,CP,I,ISTAT
     DOUBLE PRECISION :: YS,YY,SQ,YR,BETA,GNORM
     DOUBLE PRECISION,DIMENSION(D)     :: GTMP,DIAG,STP,GINT
     DOUBLE PRECISION,DIMENSION(U)     :: RHO1,ALPHA
     DOUBLE PRECISION,DIMENSION(0:U,D) :: SEARCHSTEP,GDIF
     DOUBLE PRECISION :: XINTN(D), OLDCARTN(NOPT*NIMAGE), OLDQN(D+2*NINTS)

!  Store base set of Cartesian coordinates

     OLDCARTN = X

! Transform coordinates to internals Xint.  (Note G is always in internals for unres)

     DO J1=1,NIMAGE
         DO J2=1,NRES
            C(1,J2)=X(6*(J2-1)+1+NOPT*(J1-1))
            C(2,J2)=X(6*(J2-1)+2+NOPT*(J1-1))
            C(3,J2)=X(6*(J2-1)+3+NOPT*(J1-1))
            C(1,J2+NRES)=X(6*(J2-1)+4+NOPT*(J1-1))
            C(2,J2+NRES)=X(6*(J2-1)+5+NOPT*(J1-1))
            C(3,J2+NRES)=X(6*(J2-1)+6+NOPT*(J1-1))
         END DO
         CALL UPDATEDC
         CALL INT_FROM_CART(.TRUE.,.FALSE.)
         CALL GEOM_TO_VAR(NINTS,XINTN(NINTS*(J1-1)+1:NINTS*J1))
     ENDDO

! Q
     DO J2=1,NRES
        C(1,J2)=XYZ(6*(J2-1)+1)
        C(2,J2)=XYZ(6*(J2-1)+2)
        C(3,J2)=XYZ(6*(J2-1)+3)
        C(1,J2+NRES)=XYZ(6*(J2-1)+4)
        C(2,J2+NRES)=XYZ(6*(J2-1)+5)
        C(3,J2+NRES)=XYZ(6*(J2-1)+6)
     END DO
     CALL UPDATEDC
     CALL INT_FROM_CART(.TRUE.,.FALSE.)
     CALL GEOM_TO_VAR(NINTS,OLDQN(1:NINTS))

! FIN
     DO J2=1,NRES
        C(1,J2)=XYZ(6*(J2-1)+1+NOPT*(NIMAGE+1))
        C(2,J2)=XYZ(6*(J2-1)+2+NOPT*(NIMAGE+1))
        C(3,J2)=XYZ(6*(J2-1)+3+NOPT*(NIMAGE+1))
        C(1,J2+NRES)=XYZ(6*(J2-1)+4+NOPT*(NIMAGE+1))
        C(2,J2+NRES)=XYZ(6*(J2-1)+5+NOPT*(NIMAGE+1))
        C(3,J2+NRES)=XYZ(6*(J2-1)+6+NOPT*(NIMAGE+1))
     END DO
     CALL UPDATEDC
     CALL INT_FROM_CART(.TRUE.,.FALSE.)
     CALL GEOM_TO_VAR(NINTS,OLDQN(NINTS*(NIMAGE+1)+1:NINTS*(NIMAGE+2)))

     OLDQN(NINTS+1:NINTS*(NIMAGE+1)) = XINTN(1:NINTS*NIMAGE)  !  STORE OLD INTERNALS

     CALL INITIALISE2(OLDQN)
     CALL NEBGRADIENT

! G (which points to either whole or the appropriate part of ggg) is of dimension 3*natoms*(function of nimages)
! so need to make array Gint, which matches Diag etc.
     DO J1=1,NIMAGE
        GINT(NINTS*(J1-1)+1:NINTS*J1)=G(NOPT*(J1-1)+1:NOPT*(J1-1)+NINTS)
     ENDDO

     DO NITERDONE=1,NITERMAX
     IF (BADTAU) EXIT

     MAIN: IF (NITERDONE==1) THEN
          POINT = 0
          IF (.NOT.DIAGCO) DIAG(1:D)=NEBDGUESS
          SEARCHSTEP(0,:)= -GINT*DIAG                        ! NR STEP FOR DIAGONAL INVERSE HESSIAN
          GTMP           = SEARCHSTEP(0,:)
          GNORM          = SQRT(DOT_PRODUCT(GINT,GINT))
          STP            = MIN(1.0D0/GNORM, GNORM)           ! MAKE THE FIRST GUESS FOR THE STEP LENGTH CAUTIOUS
     ELSE MAIN
          BOUND=NITERDONE-1
          IF (NITERDONE.GT.NEBMUPDATE) BOUND=NEBMUPDATE
          YS=DOT_PRODUCT( GDIF(NPT/D,:), SEARCHSTEP(NPT/D,:)  )
          IF (YS==0.0D0) YS=1.0D0
      
          ! Update estimate of diagonal inverse Hessian elements.
          ! We divide by both YS and YY at different points, so they had better not be zero!
          IF (.NOT.DIAGCO) THEN
               YY=DOT_PRODUCT( GDIF(NPT/D,:) , GDIF(NPT/D,:) )
               IF (YY==0.0D0) YY=1.0D0
!              DIAG = ABS(YS/YY)
               DIAG = YS/YY
          ELSE
               CALL CHECKINPUT2
          ENDIF
      
          ! COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980, "Updating quasi-Newton matrices with limited storage",
          ! Mathematics of Computation, Vol.35, No.151, pp. 773-782
          CP= POINT; IF (POINT==0) CP = NEBMUPDATE
          RHO1(CP)=1.0D0/YS
          GTMP = -GINT
          CP= POINT 
                   
          DO I= 1,BOUND 
               CP = CP - 1; IF (CP == -1) CP = NEBMUPDATE - 1
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
               IF (CP==NEBMUPDATE) CP=0
          ENDDO
              
          STP = 1.0D0
     ENDIF MAIN

     !  Store the new search direction
     IF (NITERDONE.GT.1) SEARCHSTEP(POINT,:)=GTMP
      
     IF (DOT_PRODUCT(GINT,GTMP)/SQRT( DOT_PRODUCT(GINT,GINT)*DOT_PRODUCT(GTMP,GTMP) ) > 0.0D0) THEN
          IF (MOREPRINTING) PRINT*,'Search direction has positive projection onto gradient - reversing step'
          GTMP=-GTMP
          SEARCHSTEP(POINT,:)=GTMP
     ENDIF
     GTMP=GINT

     !  We should apply the maximum DNEB LBFGS step to each image separately.
     DO J2=1,NIMAGE
          STEPIMAGE(J2) = SQRT(DOT_PRODUCT(SEARCHSTEP(POINT,NINTS*(J2-1)+1:NINTS*J2),SEARCHSTEP(POINT,NINTS*(J2-1)+1:NINTS*J2)))
          IF (STEPIMAGE(J2) > MAXNEBBFGS) THEN
               STP(NINTS*(J2-1)+1:NINTS*J2) = MAXNEBBFGS/STEPIMAGE(J2)
               STEPIMAGE(J2) = STP(NINTS*(J2-1)+1)*STEPIMAGE(J2)
          ENDIF
     ENDDO

     !  We now have the proposed step - update geometry and calculate new gradient
     XINTN = XINTN + STP*SEARCHSTEP(POINT,:)

     DO J1=1,NIMAGE
        CALL VAR_TO_GEOM(NINTS,XINTN(NINTS*(J1-1)+1:NINTS*J1))
! get cartesians
        CALL CHAINBUILD
        DO J2=1,NRES
           X(6*(J2-1)+1+NOPT*(J1-1))=C(1,J2)
           X(6*(J2-1)+2+NOPT*(J1-1))=C(2,J2)
           X(6*(J2-1)+3+NOPT*(J1-1))=C(3,J2)
           X(6*(J2-1)+4+NOPT*(J1-1))=C(1,J2+NRES)
           X(6*(J2-1)+5+NOPT*(J1-1))=C(2,J2+NRES)
           X(6*(J2-1)+6+NOPT*(J1-1))=C(3,J2+NRES)
        ENDDO
     ENDDO

     CALL NEBGRADIENT

! G (which points to either whole or the appropriate part of ggg) is of dim 3*natoms*(function of nimages)
! so need to make Gint, which matches Diag etc.
     DO J1=1,NIMAGE
        GINT(NINTS*(J1-1)+1:NINTS*J1)=G(NOPT*(J1-1)+1:NOPT*(J1-1)+NINTS)
     ENDDO

! DAE step finished, reset oldcart,oldq
     OLDCARTN = X
     OLDQN(NINTS+1:NINTS*(NIMAGE+1)) = XINTN

! jmc     call ImageDistribution
     CALL INTERNALIMAGEDISTRIBUTION(OLDQN)
     STEPTOT = SUM(STEPIMAGE)/NIMAGE
     IF (MOREPRINTING) THEN
        WRITE(*,'(i6,2f20.10,f8.2,a,f8.3,f9.3)') NIterDone,ETOTAL/Nimage,RMS,&
     &    AVDEV,'%',Separation, StepTot
        CALL FLUSH(6)
     ENDIF
     IF (DEBUG) CALL DUMPFILES("m")

     CALL TERMINATERUN2
     IF (EXITSTATUS > 0) EXIT


     ! Compute the new step and gradient change
     NPT=POINT*D
     SEARCHSTEP(POINT,:) = STP*SEARCHSTEP(POINT,:)
     GDIF(POINT,:)=GINT-GTMP
     POINT=POINT+1; IF (POINT==NEBMUPDATE) POINT=0
     IF (DUMPNEBXYZ.AND.MOD(NITERDONE,DUMPNEBXYZFREQ)==0) CALL RWG("w",.False.,NITERDONE)
     IF (DUMPNEBPTS.AND.MOD(NITERDONE,DUMPNEBPTSFREQ)==0) CALL SAVEBANDCOORD
     IF (DUMPNEBEOS.AND.MOD(NITERDONE,DUMPNEBEOSFREQ)==0) CALL WRITEPROFILE(NITERDONE)
     END DO

     IF (DEBUG) CALL DUMPFILES("e")

contains
     SUBROUTINE INITIALISE2(OLDQN) ! IF OLDQN IS NOT MADE AN ARGUMENT THEN IFC MISCOMPILES
          DOUBLE PRECISION :: OLDQN(D+2*NINTS)

          CALL CHECKINPUT2
          IF (DEBUG) CALL DUMPFILES("b")

          IF (MOREPRINTING) THEN
               PRINT '(/1x,a)', repeat(">",21)//" ENTERING L-BFGS MINIMISER "//repeat(">",21)
               WRITE(*,'(1x,a)') repeat("_",91)
               WRITE(*,'(a6,a20,a20,a9,2a8,a9)') 'Iter','Energy','RMS Force','Av.Dev','Path','Step'
               WRITE(*,'(1x,a)') repeat("_",91)
          ENDIF
     END SUBROUTINE INITIALISE2

     SUBROUTINE CHECKINPUT2
          
          INTEGER :: I

          IF ( D<=0 .OR. NEBMUPDATE<=0 ) THEN ! D AVAILABLE BY HOST ASSOCIATION
               WRITE(*,'(" IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)")')
               CALL TSUMMARY
               STOP
          ENDIF

          IF (NITERMAX < 0) THEN
               PRINT '(1x,a)', 'Maximal number of iterations is less than zero! Stop.'
               CALL TSUMMARY
               STOP
          ENDIF
          
          IF (DIAGCO) THEN
               PRINT*,'using estimate of the inverse diagonal elements'
               DO I=1,D
                    IF (DIAG(I)<=0.0D0) THEN
                         WRITE(*,'(" THE",I5,"-TH DIAGONAL ELEMENT OF THE INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE")') I
                         CALL TSUMMARY
                         STOP
                    ENDIF
               ENDDO
          ENDIF
     END SUBROUTINE CHECKINPUT2


     SUBROUTINE TERMINATERUN2
          EXITSTATUS=0

          NEBDGUESS=DIAG(1) ! SAVED FOR SUBSEQUENT CALLS - SHOULD BE OK FOR THE SAME SYSTEM?
          IF ( RMS <= RMSTOL.AND. NITERDONE>1.AND. .NOT.NITERDONE < NITERMIN ) THEN
               EXITSTATUS=1
               RETURN
          ENDIF
          
          IF (NITERDONE == NITERMAX) THEN
               EXITSTATUS=2
               RETURN
          ENDIF
     END SUBROUTINE TERMINATERUN2
END SUBROUTINE NEBBFGSINT

END MODULE MINIMISER3
