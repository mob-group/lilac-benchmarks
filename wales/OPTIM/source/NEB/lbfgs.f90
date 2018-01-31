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
MODULE MINIMISER1
     IMPLICIT NONE
     CONTAINS

SUBROUTINE NEBBFGS(D,U,NPERSIST,PERSISTENT)
     USE KEYMINIMIZER
     USE KEYLBFGS
     USE NEBDATA
     USE KEYNEB,ONLY: NIMAGE,MOREPRINTING,DUMPNEBXYZ,DUMPNEBXYZFREQ,DUMPNEBEOS,DUMPNEBEOSFREQ,DEBUG,&
     & DUMPNEBPTS,DUMPNEBPTSFREQ
     USE NEBUTILS
     USE GRADIENTS
     USE NEBOUTPUT
     USE PORFUNCS
     USE MODCHARMM, ONLY : CHRMMT,CHECKOMEGAT
     USE KEY, ONLY : FREEZENODEST, FREEZETOL, DESMDEBUG, MAXNEBBFGS, NEBMUPDATE, NEBDGUESS, &
          & DESMAXEJUMP,INTEPSILON, DESMAXAVGE, KADJUSTFRAC, KADJUSTFRQ, DNEBSWITCH, KADJUSTTOL, NEBRESEEDT, &
          & NEBRESEEDINT, NEBRESEEDEMAX, NEBRESEEDBMAX, NEBKFINAL, NEBFACTOR, &
          & NREPMAX, ORDERI, ORDERJ, EPSALPHA, NREPULSIVE, DISTREF, ADDREPT, NEBKINITIAL, &
          & NEBRESEEDDEL1, NEBRESEEDDEL2, NEBRESEEDPOW1, NEBRESEEDPOW2, REPPOW, &
          & BULKT, REDOTSIM, NREPMAX, COLDFUSIONLIMIT, DNEBEFRAC,VARIABLES, NEBMAXERISE
     USE INTCOMMONS, ONLY : DESMINT, NNZ, KD, NINTC, PREVDIH, DIHINFO
     USE COMMONS, ONLY: REDOPATHNEB
! hk286
     USE GENRIGID

     IMPLICIT NONE 

     INTEGER,INTENT(IN) :: D,U  ! DIMENSIONALITY OF THE PROBLEM AND NUMBER OF UPDATES
     INTEGER NPERSIST           ! number of persistent minima
     INTEGER PERSISTTHRESH      ! persistence threshold
     LOGICAL PERSISTENT(NIMAGE+2) ! logical array to identify persistent minima
     INTEGER NTIMESMIN(NIMAGE+2)  ! number of consecutive steps this image is a local minimum
     LOGICAL NOTNEW 
     DOUBLE PRECISION DIJ, DMAX, DS, DF, EWORST, EMAX, BESTE, ENOLD
     DOUBLE PRECISION, ALLOCATABLE :: REPTEMP(:)
     INTEGER, ALLOCATABLE :: IREPTEMP(:)
     INTEGER MAXIM, NDECREASE, NFAIL, NRED
     LOGICAL KNOWE, KNOWG, KNOWH, ADDATOM
     COMMON /KNOWN/ KNOWE, KNOWG, KNOWH

     INTEGER :: J2,POINT,BOUND,NPT,CP,I,ISTAT,J1,JMINUS,JPLUS,J3,JDO,J4,NPEPFAIL,NBADTOTAL
     INTEGER AT1(NATOMS),AT2(NATOMS),AT3(NATOMS),AT4(NATOMS)
     DOUBLE PRECISION :: YS,YY,SQ,YR,BETA,GNORM,DDOT,DUMMY,STPMIN,PREVGRAD,MEANSEP,EMINUS,EPLUS, &
  &                      LCOORDS(NOPT), ENERGY, INVDTOACTIVE(NATOMS), STIME
     DOUBLE PRECISION,DIMENSION(D)     :: GTMP,DIAG,STP
     DOUBLE PRECISION,DIMENSION(U)     :: RHO1,ALPHA
     DOUBLE PRECISION,DIMENSION(0:U,D) :: SEARCHSTEP,GDIF

     ! efk: for freezenodes
     DOUBLE PRECISION :: TESTG, TOTGNORM
     INTEGER :: IM

     ! efk: for internals
     LOGICAL :: FAILED, INTPTEST, SKIPPED, AMIDEFAIL
     DOUBLE PRECISION :: STEPCART(NOPT), AVGE, COORDS(NOPT)
     DOUBLE PRECISION :: TMPRMS, TESTE, LGDUMMY(NOPT)

     CALL MYCPU_TIME(STIME,.FALSE.)

     NTIMESMIN(1:NIMAGE)=0 ! number of consecutive steps the image is identified as a local minimum in the profile
!    PERSISTTHRESH=50      ! persistence identification threshold
     PERSISTTHRESH=HUGE(1)      ! persistence identification threshold
     NPERSIST=0
     PREVGRAD=1.0D100
     INTPTEST = .FALSE.
     IF (DESMDEBUG) MOREPRINTING = .TRUE.
     ADDATOM=.FALSE.
     NFAIL=0
     IF (FREEZENODEST) IMGFREEZE(1:NIMAGE)=.FALSE.
     BESTE=HUGE(1.0D0)

     CALL CHECKINPUT
     IF (DEBUG) CALL DUMPFILES("b")
     IF (DUMPNEBXYZ) CALL RWG("w",.False.,0)
!     IF (DUMPNEBEOS) CALL WRITEPROFILE(-1)  ! sn402: moved this because we
!     currently don't know the energies of the images
     IF (MOREPRINTING) THEN
          WRITE(*,'(A)') '                  Iter   Energy per image      RMS Force       Av.Dev    Path     Step length/images'
     ENDIF
     IF (PREVGRAD.LT.DNEBSWITCH) THEN
        CALL OLDNEBGRADIENT
     ELSE
        CALL NEBGRADIENT
     ENDIF
    
     IF (DUMPNEBEOS) THEN
        CALL DISTANCES  ! sn402: We need to know the distances between the images to write the profile
        CALL WRITEPROFILE(-1)
     ENDIF

     ENOLD=ETOTAL
     IF (ETOTAL/NIMAGE.LT.COLDFUSIONLIMIT) THEN
        WRITE(*,'(A,2G20.10)') ' lbfgs> Cold fusion diagnosed - step discarded, energy, limit=',ETOTAL/NIMAGE,COLDFUSIONLIMIT
        IF (DEBUG) CALL DUMPFILES("e")
        IF (KADJUSTFRQ.GT.0) WRITE(*,'(A,G20.10)') ' lbfgs> Final DNEB force constant ',NEWNEBK(1)
        RETURN
     ENDIF

     NITERDONE=1
!    DO NITERDONE=1,MAX(NITERMAX,NITERMIN)
     DO ! main do loop with counter NITERDONE
!    IF (BADTAU) EXIT
!
     MAIN: IF (NITERDONE==1) THEN
          POINT = 0
          IF (.NOT.DIAGCO) DIAG(1:D)=NEBDGUESS
          SEARCHSTEP(0,:)= -G*DIAG                           ! NR STEP FOR DIAGONAL INVERSE HESSIAN
          GTMP           = SEARCHSTEP(0,:)
          GNORM          = MAX(SQRT(DOT_PRODUCT(G,G)),1.0D-100)
          STP            = MIN(1.0D0/GNORM, GNORM)           ! MAKE THE FIRST GUESS FOR THE STEP LENGTH CAUTIOUS
     ELSE MAIN
          BOUND=NITERDONE-1
          IF (NITERDONE.GT.NEBMUPDATE) BOUND=NEBMUPDATE
          YS=DOT_PRODUCT( GDIF(NPT/D,:), SEARCHSTEP(NPT/D,:)  )
          IF (ABS(YS).LT.1.0D-20) THEN
              YS=SIGN(1.0D-20,YS)
          ENDIF 
 
          ! Update estimate of diagonal inverse Hessian elements.
          ! We divide by both YS and YY at different points, so they had better not be zero!
          IF (.NOT.DIAGCO) THEN
               YY=DOT_PRODUCT( GDIF(NPT/D,:) , GDIF(NPT/D,:) )
                IF (ABS(YY).LT.1.0D-20) THEN
                    YY=SIGN(1.0D-20,YY)
                ENDIF
               DIAG = YS/YY
          ELSE
               CALL CHECKINPUT
          ENDIF
 
          ! COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980, "Updating quasi-Newton matrices with limited storage",
          ! Mathematics of Computation, Vol.35, No.151, pp. 773-782
          CP= POINT; IF (POINT==0) CP = NEBMUPDATE
          RHO1(CP)=1.0D0/YS
          GTMP = -G
          CP= POINT 
 
          DO I= 1,BOUND 
!              CP = CP - 1; IF (CP == -1) CP = M - 1
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
!              IF (CP==M) CP=0
               IF (CP==NEBMUPDATE) CP=0
          ENDDO

          STP(1:D) = 1.0D0
     ENDIF MAIN

     !  Store the new search direction
     IF (NITERDONE.GT.1) SEARCHSTEP(POINT,:)=GTMP
      
!
! If the number of images has changed since G was declared then G is not the same
! size as Gtmp and Dot_Product cannot be used.
!
!    IF (Dot_Product(G,Gtmp)/SQRT( Dot_Product(G,G)*Dot_Product(Gtmp,Gtmp) ) > 0.0D0) THEN
!
!  Separate sqrt;s to avoid overflow.
!
     IF (DDOT(D,G,1,GTMP,1)/MAX(1.0D-100,SQRT( DDOT(D,G,1,G,1))*SQRT(DDOT(D,GTMP,1,GTMP,1)) ) > 0.0D0) THEN
          IF (MOREPRINTING) PRINT*,'Search direction has positive projection onto gradient - reversing step'
          GTMP=-GTMP
          SEARCHSTEP(POINT,:)=GTMP
     ENDIF
     GTMP=G

!  We should apply the maximum DNEB LBFGS step to each image separately.
!  However, using different scale factors for different images leads to huge
!  discontinuities! Now take the minimum scale factor for all images. DJW 26/11/07
!
!  SEARCHSTEP(POINT,1:D) the original unscaled step - vector length atomic degrees of freedom times number of interpolating images
!  STEPIMAGE(J2)         length of step for image J2 before any scaling
!  STP(1:D)              a scaling factor with equal values for all entries in each interpolating image. Scales the step
!                        for each image separately. Dimension could be number of images?!
!

     STPMIN=1.0D0
     DO J2=1,NIMAGE
          STEPIMAGE(J2) = SQRT(DOT_PRODUCT(SEARCHSTEP(POINT,NOPT*(J2-1)+1:NOPT*J2),SEARCHSTEP(POINT,NOPT*(J2-1)+1:NOPT*J2)))
          DUMMY=STEPIMAGE(J2)
          IF (STEPIMAGE(J2) > MAXNEBBFGS) THEN
               STP(NOPT*(J2-1)+1:NOPT*J2) = MAXNEBBFGS/STEPIMAGE(J2)
               STPMIN=MIN(STPMIN,STP(NOPT*(J2-1)+1))
          ENDIF
!         PRINT '(A,I8,3G20.10)','image,initial step size,STP,prod=',J2,DUMMY,STP(NOPT*(J2-1)+1),STEPIMAGE(J2)*STP(NOPT*(J2-1)+1)
     ENDDO
     STP(1:D)=STPMIN
! EFK: decide whether to freeze some nodes
     IF (FREEZENODEST.AND.(.NOT.NEBRESEEDT)) THEN
        TOTGNORM = SQRT(DOT_PRODUCT(G(1:NOPT*NIMAGE),G(1:NOPT*NIMAGE))/NIMAGE)
        DO IM = 1,NIMAGE
           TESTG = SQRT(DOT_PRODUCT(G(NOPT*(IM-1)+1:NOPT*IM),G(NOPT*(IM-1)+1:NOPT*IM)))
           IF(TESTG/TOTGNORM.LT.FREEZETOL) THEN
              IF (MOREPRINTING) PRINT '(A,I6,2G20.10)', ' lbfgs> Freezing image: ', IM, TESTG, TOTGNORM
              IMGFREEZE(IM)=.TRUE.
              STEPIMAGE(IM)=0.0D0
              STP(NOPT*(IM-1)+1:NOPT*IM)=0.0D0
           ELSE
              IMGFREEZE(IM) = .FALSE.
           ENDIF
        ENDDO
        IF (REDOPATHNEB) THEN
           IMGFREEZE(REDOTSIM)=.TRUE.
           STEPIMAGE(REDOTSIM)=0.0D0
           STP(NOPT*(REDOTSIM-1)+1:NOPT*REDOTSIM)=0.0D0
           IF (MOREPRINTING) PRINT '(A,I6)',' lbfgs> Freezing variable image ',REDOTSIM ! REDOTSIM+1 if we count starting image!
        ENDIF
     ELSEIF (FREEZENODEST.AND.NEBRESEEDT) THEN
        DO IM=1,NIMAGE
           IF (IMGFREEZE(IM)) THEN
               STEPIMAGE(IM)=0.0D0
               STP(NOPT*(IM-1)+1:NOPT*IM) = 0.0D0
           ENDIF
        ENDDO
     ENDIF
     !  We now have the proposed step - update geometry and calculate new gradient
     IF (DESMINT) THEN
        DO IM = 1,NIMAGE
           FAILED = .TRUE.
           DO WHILE (FAILED)
              PREVDIH => DIHINFO(IM+1,:)
              CALL TRANSBACKDELTA(STP(NOPT*(IM-1)+1:NOPT*IM)*SEARCHSTEP(POINT,NOPT*(IM-1)+1:NOPT*IM),&
                   & STEPCART,XCART(NOPT*(IM-1)+1:NOPT*IM),NINTC, &
                   & NOPT,NNZ,KD,FAILED,INTPTEST,INTEPSILON)
              CALL POTENTIAL(XCART(NOPT*(IM-1)+1:NOPT*IM) +STEPCART,TESTE,LGDUMMY,.FALSE.,.FALSE.,TMPRMS,.FALSE.,.FALSE.)

              IF (TESTE-EEE(IM+1).GT.DESMAXEJUMP.AND.EEE(IM+1).LT.0) THEN
                 IF (MOREPRINTING) print*, 'Too great an energy increase, skip image.', IM, TESTE, EEE(IM+1)
                 STP(NOPT*(IM-1)+1:NOPT*IM) = 0.0D0     
                 STEPIMAGE(IM) = 0.0D0
                 SKIPPED = .TRUE.
              ENDIF
              IF (FAILED) THEN
                 STP(NOPT*(IM-1)+1:NOPT*IM) = STP(NOPT*(IM-1)+1:NOPT*IM)*0.1D0        
                 STEPIMAGE(IM) = STEPIMAGE(IM)*0.1D0
                 IF (MOREPRINTING) print*, 'neblbfgs>> Decreasing step to ', IM, STEPIMAGE(IM)                 
              ENDIF
           ENDDO
           IF (.NOT.SKIPPED) XCART(NOPT*(IM-1)+1:NOPT*IM) = XCART(NOPT*(IM-1)+1:NOPT*IM) +STEPCART
           SKIPPED = .FALSE.
        ENDDO
     ENDIF
     NDECREASE=0
20   X(1:D) = X(1:D) + STP(1:D)*SEARCHSTEP(POINT,1:D)

     IF (PREVGRAD.LT.DNEBSWITCH) THEN
        CALL OLDNEBGRADIENT
     ELSE
        CALL NEBGRADIENT
     ENDIF

     NRED = 0
     DO WHILE ((ETOTAL-ENOLD).GT.NEBMAXERISE*NIMAGE)
        X(1:D) = X(1:D) - STP(1:D)*SEARCHSTEP(POINT,1:D)
        STP(1:D) = STP(1:D)/10
        X(1:D) = X(1:D) + STP(1:D)*SEARCHSTEP(POINT,1:D)
        NRED = NRED + 1
        IF (PREVGRAD.LT.DNEBSWITCH) THEN
           CALL OLDNEBGRADIENT
        ELSE
           CALL NEBGRADIENT
        END IF
        IF (NRED .EQ. 10) THEN
            EXIT
        END IF
     END DO

     ENOLD=ETOTAL

     IF (ETOTAL/NIMAGE.LT.COLDFUSIONLIMIT) THEN
        WRITE(*,'(A,2G20.10)') ' lbfgs> Cold fusion diagnosed - step discarded, energy, limit=',ETOTAL/NIMAGE,COLDFUSIONLIMIT
        IF (DEBUG) CALL DUMPFILES("e")
        IF (KADJUSTFRQ.GT.0) PRINT '(A,G20.10)',' lbfgs> Final DNEB force constant ',NEWNEBK(1)
        RETURN
     ENDIF

     CALL IMAGEDISTRIBUTION
!
!  Dynamic adjustment of NEWNEBK values. Local values are changed by KADJUSTFRAC if the
!  corresponding separation is outside a fraction KADJUSTTOL of the average value.
!  The adjustment is done every KADJUSTFRQ cycles.
!
     IF (KADJUSTFRQ.GT.0) THEN
        IF (MOD(NITERDONE,KADJUSTFRQ).EQ.0) THEN ! dynamic adjustment of K
           IF (AVDEV.GT.KADJUSTTOL) THEN
              NEWNEBK(1:NIMAGE+1)=NEWNEBK(1:NIMAGE+1)*KADJUSTFRAC
              IF (DEBUG) PRINT '(A,G20.10)',' lbfgs> Increasing DNEB force constant to ',NEWNEBK(1)
           ELSEIF (AVDEV.LT.KADJUSTTOL) THEN
              NEWNEBK(1:NIMAGE+1)=NEWNEBK(1:NIMAGE+1)/KADJUSTFRAC
              IF (DEBUG) PRINT '(A,G20.10)',' lbfgs> Decreasing DNEB force constant to ',NEWNEBK(1)
           ENDIF

!            MEANSEP=SEPARATION/(NIMAGE+1)
!            DO J1=1,NIMAGE+1
!               IF (ABS((DVEC(J1)-MEANSEP)/MAX(MEANSEP,1.0D-100)).GT.KADJUSTTOL) THEN 
!                  IF (DVEC(J1).GT.MEANSEP) THEN
!                     NEWNEBK(J1)=MIN(1.0D2 ,NEWNEBK(J1)*(1.0D0+KADJUSTFRAC))
!                     PRINT '(A,F12.4,A,I8,A,F12.4,A,F12.4)',' lbfgs> adjusting NEB force constant to ',NEWNEBK(J1), &
!    &                                                    ' for gap ',J1,' value=',DVEC(J1),' mean value=',MEANSEP
!                  ENDIF
! !                IF (DVEC(J1).LT.MEANSEP) NEWNEBK(J1)=MAX(1.0D-2,NEWNEBK(J1)*(1.0D0-KADJUSTFRAC))
!               ENDIF
!            ENDDO
        ENDIF
     ELSE
!
!  Alternative turning on of NEWNEBK.
!
        DO J1=1,NIMAGE+1
           NEWNEBK(J1)=MIN(NEBKFINAL,NEWNEBK(J1)*NEBFACTOR)
        ENDDO
!       IF (DEBUG) PRINT '(A,G20.10)','lbfgs> DNEB force constant for image 1 is ',NEWNEBK(1)
     ENDIF

!  SEARCHSTEP(POINT,1:D) the original unscaled step - vector length atomic degrees of freedom times number of interpolating images
!  STEPIMAGE(J2)         length of step for image J2 before any scaling
!  STP(1:D)              a scaling factor with equal values for all entries in each interpolating image. Scales the step

     STEPTOT=0.0D0
     DO J1=1,D
        STEPTOT=STEPTOT+(SEARCHSTEP(POINT,J1)*STP(J1))**2
     ENDDO
     STEPTOT=SQRT(STEPTOT)/NIMAGE

     IF (MOREPRINTING) THEN
        IF (PREVGRAD.LT.DNEBSWITCH) THEN
           WRITE(*,'(A,I6,2G20.10,F8.2,A,F8.3,2X,G18.5,A)') ' lbfgs> steps: ',NITERDONE,ETOTAL/NIMAGE,RMS,AVDEV,'%',SEPARATION, &
  &                                                     STEPTOT, ' singly-nudged'
        ELSE
           WRITE(*,'(A,I6,2G20.10,F8.2,A,F8.3,2X,G18.5)') ' lbfgs> steps: ',NITERDONE,ETOTAL/NIMAGE,RMS,AVDEV,'%',SEPARATION, STEPTOT
        ENDIF
        CALL FLUSH(6)
     ENDIF
     IF (ETOTAL/NIMAGE.LT.BESTE) BESTE=ETOTAL/NIMAGE
     IF (DEBUG) CALL DUMPFILES("m")


     CALL TERMINATERUN

!
!  Check to see if the profile actually has a maximum in it. If not try
!  reducing the force constants by 10 per cent.
!
!    EMAX=-1.0D100
!    DO J1=1,NIMAGE+2
!       IF (EEE(J1).GT.EMAX) THEN
!          EMAX=EEE(J1)
!          MAXIM=J1
!       ENDIF
!    ENDDO
!    PRINT '(A,I6,A,G20.10)',' lbfgs> MAXIM=',MAXIM,' energy=',EMAX
!    IF ((MAXIM.EQ.1).OR.(MAXIM.EQ.NIMAGE+2)) THEN
!       DO J1=1,NIMAGE+1
!          NEWNEBK(J1)=NEWNEBK(J1)*0.9D0
!       ENDDO
!       PRINT '(A,G20.10)',' lbfgs> No maximum in DNEB profile - reducing force constant to ',NEWNEBK(1)
!       EXITSTATUS=0
!    ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     IF (EXITSTATUS > 0.AND.((ETOTAL/NIMAGE.LT.BESTE*DNEBEFRAC).OR.(NITERDONE.GT.10*NITERMAX).OR.DNEBEFRAC.EQ.0)) THEN  
        IF (KADJUSTFRQ.GT.0) PRINT '(A,G20.10)',' lbfgs> Final DNEB force constant ',NEWNEBK(1)
        EXIT
     ENDIF
     777 CONTINUE
!
! Compute the new step and gradient change
!
     NPT=POINT*D
     SEARCHSTEP(POINT,:) = STP*SEARCHSTEP(POINT,:)
     GDIF(POINT,:)=G-GTMP
     POINT=POINT+1; IF (POINT==NEBMUPDATE) POINT=0
!
! Optional reseeding of bad images.
! Note that the numbering from 1 to NIMAGE in X and EIMAGE is different
! from XYZ and EEE, which go from 1 to NIMAGE+2.
!
! Should delete all this and tidy up if it simply doesn't work! DJW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
     IF (NEBRESEEDT.AND.MOD(NITERDONE,NEBRESEEDINT)==0) THEN
        IF (BULKT) THEN
           PRINT '(A)','lbfgs> ERROR - minimum image convention needs to be coded in NEB/lbfgs.f90'
           STOP
        ENDIF
!! DJW 
!       IF (.NOT.ADDREPT) THEN
           ADDREPT=.TRUE.
           PRINT '(A)','lbfgs> Turning on repulsive/constraint terms for bad images/cis peptide bonds'
!       ELSE
!
! Turn off repulsive terms for a thaw phase.
!
!! DJW try keeping them ??
!
!          ADDREPT=.FALSE.
!          PRINT '(A)','lbfgs> Turning off repulsive/constraint terms for bad images/cis peptide bonds'
!          FREEZENODEST=.FALSE.
!          IMGFREEZE(1:NIMAGE)=.FALSE.
!          PRINT '(A,I6,A)','lbfgs> Unfreezing all images'
!          GOTO 555
!       ENDIF
        BADIMAGE(1:NIMAGE+2)=.FALSE.
        BADPEPTIDE(1:NIMAGE+2)=.FALSE.
        NPEPFAIL=0
        NBADTOTAL=0
        DO J1=1,NIMAGE+2
           IF (CHRMMT.AND.CHECKOMEGAT) THEN 
               AMIDEFAIL=.FALSE.
! hk286 
               IF (RIGIDINIT) THEN
                  CALL TRANSFORMRIGIDTOC (1, NRIGIDBODY, COORDS(1:3*NATOMS), XYZ((J1-1)*3*NATOMS+1:(J1-1)*3*NATOMS+DEGFREEDOMS))
               ELSE
! hk286
                  COORDS(1:NOPT)=XYZ((J1-1)*NOPT+1:J1*NOPT)
               ENDIF
!              CALL CHECKOMEGA(COORDS,AMIDEFAIL)
               CALL DJWCHECKOMEGA(COORDS,AMIDEFAIL,NPEPFAIL,AT1,AT2,AT3,AT4)
               IF (AMIDEFAIL) THEN
                  PRINT '(2(A,I6))',' potential> WARNING *** cis peptide bond(s) detected for image ',J1,' total=',NPEPFAIL
                  BADPEPTIDE(J1)=.TRUE.
              ENDIF
           ENDIF
           IF (EEE(J1).GT.NEBRESEEDEMAX) BADIMAGE(J1)=.TRUE.
           IF ((EEE(J1)-EEE(1).GT.NEBRESEEDBMAX).AND.(EEE(J1)-EEE(NIMAGE+2).GT.NEBRESEEDBMAX)) BADIMAGE(J1)=.TRUE.
           IF (BADIMAGE(J1).OR.BADPEPTIDE(J1)) NBADTOTAL=NBADTOTAL+1
        ENDDO
        IF (BADIMAGE(1).OR.BADPEPTIDE(1)) PRINT '(A)','lbfgs> WARNING - bad starting image will not be replaced'
        IF (BADIMAGE(NIMAGE+2).OR.BADPEPTIDE(NIMAGE+2)) PRINT '(A)','lbfgs> WARNING - bad final image will not be replaced'
        J1=2
        imageloop: DO 
           IF (BADIMAGE(J1).AND.(.NOT.VARIABLES)) THEN
              DO J2=J1-1,1,-1
                 IF ((.NOT.BADIMAGE(J2)).OR.(J2.EQ.1)) THEN
                    JMINUS=J2
                    EMINUS=EEE(J2)
                    EXIT
                 ENDIF
              ENDDO
              DO J2=J1+1,NIMAGE+2
                 IF ((.NOT.BADIMAGE(J2)).OR.(J2.EQ.NIMAGE+2)) THEN
                    JPLUS=J2
                    EPLUS=EEE(J2)
                    EXIT
                 ENDIF
              ENDDO
              PRINT '(A,I6,A,G20.10,2(A,I6),A,2F20.10)','lbfgs> DNEB bad image ',J1,' energy=',EEE(J1),' bracketed by ', &
  &                                    JMINUS,',',JPLUS,' energies ',EMINUS,EPLUS
!
!  Consider only the highest energy bad image to avoid duplication.
!
!             JDO=J1
              EWORST=-1.0D100
              DO J2=JMINUS+1,JPLUS-1
                 IF (EEE(J2).GT.EWORST) THEN
                    EWORST=EEE(J2)
                    JDO=J2
                 ENDIF
              ENDDO
              PRINT '(A,I6,A,G20.10)','lbfgs> highest bad image in this range is number ',JDO,' energy=',EEE(JDO)
              DMAX=-1.0D0
              DO J2=1,NATOMS
                 J3loop: DO J3=J2+1,NATOMS
                    DO J4=1,NREPULSIVE
                       IF ((ORDERI(J4).EQ.J2).AND.(ORDERJ(J4).EQ.J3)) CYCLE J3loop
                    ENDDO
                    DIJ=SQRT((XYZ((JDO-1)*NOPT+3*(J2-1)+1)-XYZ((JDO-1)*NOPT+3*(J3-1)+1))**2 &
  &                         +(XYZ((JDO-1)*NOPT+3*(J2-1)+2)-XYZ((JDO-1)*NOPT+3*(J3-1)+2))**2 &
  &                         +(XYZ((JDO-1)*NOPT+3*(J2-1)+3)-XYZ((JDO-1)*NOPT+3*(J3-1)+3))**2) 
                    DS =SQRT((XYZ(3*(J2-1)+1)-XYZ(3*(J3-1)+1))**2 &
  &                         +(XYZ(3*(J2-1)+2)-XYZ(3*(J3-1)+2))**2 &
  &                         +(XYZ(3*(J2-1)+3)-XYZ(3*(J3-1)+3))**2) 
                    DF =SQRT((XYZ((NIMAGE+1)*NOPT+3*(J2-1)+1)-XYZ((NIMAGE+1)*NOPT+3*(J3-1)+1))**2 &
  &                         +(XYZ((NIMAGE+1)*NOPT+3*(J2-1)+2)-XYZ((NIMAGE+1)*NOPT+3*(J3-1)+2))**2 &
  &                         +(XYZ((NIMAGE+1)*NOPT+3*(J2-1)+3)-XYZ((NIMAGE+1)*NOPT+3*(J3-1)+3))**2) 
                    DUMMY=MIN(ABS(DIJ-DF),ABS(DIJ-DS))/DIJ
!
!  We need DIJ to be less than both or greater than both distances in the end points.
!
                    IF (((DIJ-DF)*(DIJ-DS).GT.0.0D0).AND.(DUMMY.GT.DMAX)) THEN
                       DMAX=DUMMY
!                      PRINT '(A,2I6,5G20.10)','lbfgs> I,J,DIJ,DS,DF,DUMMY=',J2,J3,DIJ,DS,DF,DUMMY
                    ENDIF
                 ENDDO J3loop
              ENDDO
              DO J2=1,NATOMS
                 DO J3=J2+1,NATOMS
                    DIJ=SQRT((XYZ((JDO-1)*NOPT+3*(J2-1)+1)-XYZ((JDO-1)*NOPT+3*(J3-1)+1))**2 &
  &                         +(XYZ((JDO-1)*NOPT+3*(J2-1)+2)-XYZ((JDO-1)*NOPT+3*(J3-1)+2))**2 &
  &                         +(XYZ((JDO-1)*NOPT+3*(J2-1)+3)-XYZ((JDO-1)*NOPT+3*(J3-1)+3))**2) 
                    DS =SQRT((XYZ(3*(J2-1)+1)-XYZ(3*(J3-1)+1))**2 &
  &                         +(XYZ(3*(J2-1)+2)-XYZ(3*(J3-1)+2))**2 &
  &                         +(XYZ(3*(J2-1)+3)-XYZ(3*(J3-1)+3))**2) 
                    DF =SQRT((XYZ((NIMAGE+1)*NOPT+3*(J2-1)+1)-XYZ((NIMAGE+1)*NOPT+3*(J3-1)+1))**2 &
  &                         +(XYZ((NIMAGE+1)*NOPT+3*(J2-1)+2)-XYZ((NIMAGE+1)*NOPT+3*(J3-1)+2))**2 &
  &                         +(XYZ((NIMAGE+1)*NOPT+3*(J2-1)+3)-XYZ((NIMAGE+1)*NOPT+3*(J3-1)+3))**2) 
                    DUMMY=MIN(ABS(DIJ-DF),ABS(DIJ-DS))/DIJ
                    IF ((DIJ-DF)*(DIJ-DS).GT.0.0D0) THEN
                       IF (DUMMY.GT.0.9D0*DMAX) THEN
                          NOTNEW=.FALSE.
                          reploop: DO J4=1,NREPULSIVE
                             IF ((ORDERI(J4).EQ.J2).AND.(ORDERJ(J4).EQ.J3)) THEN
                                NOTNEW=.TRUE.
                                EXIT reploop 
                             ENDIF
                          ENDDO reploop
                          IF (.NOT.NOTNEW) THEN
!
!  Add order parameter for bad image to current list.
!
!! DJW debug maximum two repulsions
!
                             IF (NREPULSIVE.EQ.4) GOTO 888
                             NREPULSIVE=NREPULSIVE+1
                             PRINT '(A,2I6,A,I6)','lbfgs> Adding repulsive order parameter for atoms ',J2,J3,'  total=',NREPULSIVE
                             IF (NREPULSIVE.GT.NREPMAX) THEN
                                ALLOCATE(IREPTEMP(NREPMAX))
               
                                IREPTEMP(1:NREPMAX)=ORDERI(1:NREPMAX)
                                DEALLOCATE(ORDERI)
                                ALLOCATE(ORDERI(2*NREPMAX))
                                ORDERI(1:NREPMAX)=IREPTEMP(1:NREPMAX)
               
                                IREPTEMP(1:NREPMAX)=ORDERJ(1:NREPMAX)
                                DEALLOCATE(ORDERJ)
                                ALLOCATE(ORDERJ(2*NREPMAX))
                                ORDERJ(1:NREPMAX)=IREPTEMP(1:NREPMAX)
               
                                IREPTEMP(1:NREPMAX)=REPPOW(1:NREPMAX)
                                DEALLOCATE(REPPOW)
                                ALLOCATE(REPPOW(2*NREPMAX))
                                REPPOW(1:NREPMAX)=IREPTEMP(1:NREPMAX)
               
                                DEALLOCATE(IREPTEMP)
                                ALLOCATE(REPTEMP(1:NREPMAX))
               
                                REPTEMP(1:NREPMAX)=EPSALPHA(1:NREPMAX)
                                DEALLOCATE(EPSALPHA)
                                ALLOCATE(EPSALPHA(2*NREPMAX))
                                EPSALPHA(1:NREPMAX)=REPTEMP(1:NREPMAX)
               
                                REPTEMP(1:NREPMAX)=DISTREF(1:NREPMAX)
                                DEALLOCATE(DISTREF)
                                ALLOCATE(DISTREF(2*NREPMAX))
                                DISTREF(1:NREPMAX)=REPTEMP(1:NREPMAX)

                                NREPMAX=2*NREPMAX
                                DEALLOCATE(REPTEMP)
                             ENDIF
                             ORDERI(NREPULSIVE)=J2
                             ORDERJ(NREPULSIVE)=J3
                             EPSALPHA(NREPULSIVE)=NEBRESEEDDEL1
                             REPPOW(NREPULSIVE)=NEBRESEEDPOW1
                             DISTREF(NREPULSIVE)=DIJ
                             PRINT '(A,2G20.10)','lbfgs> bond length and factor for bad image are   ', &
  &                                                DIJ,EPSALPHA(NREPULSIVE)
888                          CONTINUE
                          ENDIF
                       ENDIF
                    ENDIF
                 ENDDO
              ENDDO

              J2=JMINUS
              IF (EPLUS.GT.EMINUS) J2=JPLUS

              X((JDO-2)*NOPT+1:(JDO-1)*NOPT)=XYZ((J2-1)*NOPT+1:J2*NOPT)
              EIMAGE(JDO-1)=EIMAGE(J2-1)
              PRINT '(A,I6,A,G20.10)','lbfgs> image ',JDO,' reinterpolated energy=',EIMAGE(JDO-1)
!
! Change all the other bad images in this range.
!
              DO J3=JMINUS+1,JDO-1
!
!  Reset to MINUS end point
!
                 LCOORDS(1:NOPT)=XYZ((JMINUS-1)*NOPT+1:JMINUS*NOPT)
!                XYZ((J3-1)*NOPT+1:J3*NOPT)=LCOORDS(1:NOPT)
                 X((J3-2)*NOPT+1:(J3-1)*NOPT)=LCOORDS(1:NOPT)
                 CALL POTENTIAL(LCOORDS,ENERGY,LGDUMMY,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!                EEE(J3)=ENERGY
                 EIMAGE(J3-1)=ENERGY
                 PRINT '(A,I6,A,G20.10)','lbfgs> image ',J3,' reinterpolated energy=',ENERGY
              ENDDO
              DO J3=JDO+1,JPLUS-1
!
!  Reset to PLUS end point
!
                 LCOORDS(1:NOPT)=XYZ((JPLUS-1)*NOPT+1:JPLUS*NOPT)
!                XYZ((J3-1)*NOPT+1:J3*NOPT)=LCOORDS(1:NOPT)
                 X((J3-2)*NOPT+1:(J3-1)*NOPT)=LCOORDS(1:NOPT)
                 CALL POTENTIAL(LCOORDS,ENERGY,LGDUMMY,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!                EEE(J3)=ENERGY
                 EIMAGE(J3-1)=ENERGY
                 PRINT '(A,I6,A,G20.10)','lbfgs> image ',J3,' reinterpolated energy=',ENERGY
              ENDDO
!             NEWNEBK(JMINUS:JPLUS-1)=NEBKINITIAL
!             PRINT '(A,G20.10)','lbfgs> for bad images DNEB force constants reset to ',NEBKINITIAL
              J1=JPLUS+1
           ELSE
              J1=J1+1
           ENDIF
           IF (J1.GE.NIMAGE+2) EXIT imageloop
        ENDDO imageloop
!
!  Now add harmonic constraints for any cis peptide bonds in this range if we haven't covered them already.
!
        DO J2=1,NPEPFAIL
           PRINT '(A,I6,A,4I6)','lbfgs> checking for new constraint terms for cis peptide number ',J2,' atoms ', &
           &     AT1(J2),AT2(J2),AT3(J2),AT4(J2)
           NOTNEW=.FALSE.
           reploop2: DO J4=1,NREPULSIVE
              IF ((ORDERI(J4).EQ.AT1(J2)).AND.(ORDERJ(J4).EQ.AT4(J2))) THEN
                 NOTNEW=.TRUE.
                 EXIT reploop2
              ENDIF
              IF ((ORDERI(J4).EQ.AT4(J2)).AND.(ORDERJ(J4).EQ.AT1(J2))) THEN
                 NOTNEW=.TRUE.
                 EXIT reploop2
              ENDIF
           ENDDO reploop2
           IF (.NOT.NOTNEW) THEN
!
!  Add order parameter for bad image to current list.
!
              NREPULSIVE=NREPULSIVE+3
              IF (NREPULSIVE.GT.NREPMAX) THEN
                 ALLOCATE(IREPTEMP(NREPMAX))
               
                 IREPTEMP(1:NREPMAX)=ORDERI(1:NREPMAX)
                 DEALLOCATE(ORDERI)
                 ALLOCATE(ORDERI(2*NREPMAX))
                 ORDERI(1:NREPMAX)=IREPTEMP(1:NREPMAX)
               
                 IREPTEMP(1:NREPMAX)=ORDERJ(1:NREPMAX)
                 DEALLOCATE(ORDERJ)
                 ALLOCATE(ORDERJ(2*NREPMAX))
                 ORDERJ(1:NREPMAX)=IREPTEMP(1:NREPMAX)
               
                 IREPTEMP(1:NREPMAX)=REPPOW(1:NREPMAX)
                 DEALLOCATE(REPPOW)
                 ALLOCATE(REPPOW(2*NREPMAX))
                 REPPOW(1:NREPMAX)=IREPTEMP(1:NREPMAX)
               
                 DEALLOCATE(IREPTEMP)
                 ALLOCATE(REPTEMP(1:NREPMAX))
               
                 REPTEMP(1:NREPMAX)=EPSALPHA(1:NREPMAX)
                 DEALLOCATE(EPSALPHA)
                 ALLOCATE(EPSALPHA(2*NREPMAX))
                 EPSALPHA(1:NREPMAX)=REPTEMP(1:NREPMAX)
               
                 REPTEMP(1:NREPMAX)=DISTREF(1:NREPMAX)
                 DEALLOCATE(DISTREF)
                 ALLOCATE(DISTREF(2*NREPMAX))
                 DISTREF(1:NREPMAX)=REPTEMP(1:NREPMAX)

                 NREPMAX=2*NREPMAX
                 DEALLOCATE(REPTEMP)
              ENDIF
              PRINT '(A,2I6,A,I6)','lbfgs> cis Adding harmonic constraint for atoms ',AT1(J2),AT4(J2),'  total=',NREPULSIVE-2
              ORDERI(NREPULSIVE-2)=AT1(J2)
              ORDERJ(NREPULSIVE-2)=AT4(J2)
              EPSALPHA(NREPULSIVE-2)=NEBRESEEDDEL2
              REPPOW(NREPULSIVE-2)=NEBRESEEDPOW2
              DIJ=3.8D0
              DISTREF(NREPULSIVE-2)=DIJ 
              PRINT '(A,2G20.10)','lbfgs> bond length and factor for bad image are   ', &
  &                                       DIJ,EPSALPHA(NREPULSIVE-2)

              PRINT '(A,2I6,A,I6)','lbfgs> cis Adding harmonic constraint for atoms ',AT1(J2),AT3(J2),'  total=',NREPULSIVE-1
              ORDERI(NREPULSIVE-1)=AT1(J2)
              ORDERJ(NREPULSIVE-1)=AT3(J2)
              EPSALPHA(NREPULSIVE-1)=NEBRESEEDDEL2
              REPPOW(NREPULSIVE-1)=NEBRESEEDPOW2
              DIJ=2.5D0
              DISTREF(NREPULSIVE-1)=DIJ 
              PRINT '(A,2G20.10)','lbfgs> bond length and factor for bad image are   ', &
  &                                       DIJ,EPSALPHA(NREPULSIVE-1)

              PRINT '(A,2I6,A,I6)','lbfgs> cis Adding harmonic constraint for atoms ',AT2(J2),AT4(J2),'  total=',NREPULSIVE
              ORDERI(NREPULSIVE)=AT2(J2)
              ORDERJ(NREPULSIVE)=AT4(J2)
              EPSALPHA(NREPULSIVE)=NEBRESEEDDEL2
              REPPOW(NREPULSIVE)=NEBRESEEDPOW2
              DIJ=2.4D0
              DISTREF(NREPULSIVE)=DIJ 
              PRINT '(A,2G20.10)','lbfgs> bond length and factor for bad image are   ', &
  &                                       DIJ,EPSALPHA(NREPULSIVE)
           ENDIF
        ENDDO

!
! debug comment DJW
!
!
          FREEZENODEST=.FALSE.
!
!         IF (NBADTOTAL.EQ.-1) THEN
!            FREEZENODEST=.FALSE.
!            IMGFREEZE(1:NIMAGE)=.FALSE.
!            PRINT '(A,I6,A)','lbfgs> No bad images - unfreezing everything and removing repulsive and constraint terms'
!  
! !  Could turn off the repulsive and constraint terms as well?
! !
!            ADDREPT=.FALSE.
!         ELSE
!            FREEZENODEST=.TRUE.
!            NIMAGEFREEZE=0
!            DO J2=2,NIMAGE+1
!               IF (BADIMAGE(J2).OR.BADPEPTIDE(J2)) THEN
!                  IMGFREEZE(J2-1)=.FALSE.
!               ELSE
!                  IF ((.NOT.BADIMAGE(J2-1)).OR.(.NOT.BADIMAGE(J2+1))) THEN
!                     NIMAGEFREEZE=NIMAGEFREEZE+1
!                     IMGFREEZE(J2-1)=.TRUE.
! !                   STEPIMAGE(J2-1)=0.0D0
! !                   STP(NOPT*(J2-1-1)+1:NOPT*J2-1)=0.0D0
!                  ENDIF
!               ENDIF
!             ENDDO
!            PRINT '(A,I6,A)','lbfgs> Freezing ',NIMAGE-NBADTOTAL,' good images'
!         ENDIF

555     CONTINUE
     ENDIF

     IF (DUMPNEBXYZ.AND.MOD(NITERDONE,DUMPNEBXYZFREQ)==0) CALL RWG("w",.False.,NITERDONE)
     IF (DUMPNEBPTS.AND.MOD(NITERDONE,DUMPNEBPTSFREQ)==0) CALL SAVEBANDCOORD
     IF (DUMPNEBEOS.AND.MOD(NITERDONE,DUMPNEBEOSFREQ)==0) CALL WRITEPROFILE(NITERDONE)
!
!  Check for persistent minima
!
     DO J1=2,NIMAGE+1
        IF ((EEE(J1).LT.EEE(J1-1)).AND.(EEE(J1).LT.EEE(J1+1))) THEN
           NTIMESMIN(J1)=NTIMESMIN(J1)+1
        ELSE
           NTIMESMIN(J1)=0
        ENDIF
     ENDDO
     NPERSIST=0
     PERSISTENT(1:NIMAGE+2)=.FALSE.
!    IF (MOD(NITERDONE,PERSISTTHRESH)/2.EQ.0) THEN ! to give more than one image a chance to become persistent
!       DO J1=2,NIMAGE+1
!          IF (NTIMESMIN(J1).GE.PERSISTTHRESH) THEN
!             NPERSIST=NPERSIST+1
!             PERSISTENT(J1)=.TRUE.
!             PRINT '(3(A,I8))',' nebbfgs> identified persistent minimum, image number ',J1,' minimum for ',NTIMESMIN(J1),' steps'
!          ENDIF
!       ENDDO
!    ENDIF
     IF (NPERSIST.GT.0) THEN
        IF (DEBUG) CALL DUMPFILES("e")
        IF (KADJUSTFRQ.GT.0) PRINT '(A,G20.10)',' lbfgs> Final DNEB force constant ',NEWNEBK(1)
        RETURN
     ENDIF
!
!  End of check for persistent minima
!
     PREVGRAD=RMS

     NITERDONE=NITERDONE+1
     IF (NITERDONE.GT.MAX(NITERMAX,NITERMIN)) THEN
        IF ((ETOTAL/NIMAGE.LT.BESTE*DNEBEFRAC).OR.(NITERDONE.GT.10*NITERMAX).OR.DNEBEFRAC.EQ.0) EXIT
     ENDIF

  ENDDO ! end of main do loop over variable NITERDONE

     IF (DEBUG) CALL DUMPFILES("e")

     CONTAINS

     SUBROUTINE CHECKINPUT
          INTEGER :: I

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
     END SUBROUTINE CHECKINPUT

     SUBROUTINE TERMINATERUN
          EXITSTATUS=0
          NEBDGUESS=DIAG(1) ! should be ok for subsequent runs of the same system DJW
          IF ( RMS <= RMSTOL.AND. NITERDONE>1.AND. .NOT.NITERDONE < NITERMIN ) THEN
               EXITSTATUS=1
               RETURN
          ENDIF
          
          AVGE = SUM(EEE(2:NIMAGE+1))/NIMAGE
          IF (NITERDONE == NITERMAX.AND.NITERDONE.GE.NITERMIN.AND.AVGE.LT.DESMAXAVGE) THEN
               EXITSTATUS=2
               RETURN
          ENDIF
     END SUBROUTINE TERMINATERUN

END SUBROUTINE NEBBFGS

END MODULE MINIMISER1
