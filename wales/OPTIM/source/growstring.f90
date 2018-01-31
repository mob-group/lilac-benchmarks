MODULE GROWSTRINGUTILS
  USE COMMONS
  USE KEY, ONLY : DESMDEBUG, PERMDIST, TWOD, RIGIDBODY, BULKT, &
       & FREEZENODEST, FREEZETOL,INTEPSILON, DESMAXEJUMP,      &
       & DESMAXAVGE,CLSTRINGT,CLSTRINGTST
  USE CUBSPLSTRING
  USE GSDATA
  USE SPFUNCTS, ONLY : SHIFTZERO, ROTATEZERO
  USE INTCOMMONS, ONLY : NINTC, NINTIM, INTINTERPT, NATINT, DESMINT, NNZ, KD, NDIH, ALIGNDIR, &
       & DIHINFO, PREVDIH, BACKTCUTOFF, INTERPBACKTCUT, MINBACKTCUT
  USE PORFUNCS
!
! Moved this line to prevent a segmentation fault in gfortran
!
! USE INTCUTILS, ONLY : CART2INT, INTINTERPOLATE
  IMPLICIT NONE
  
CONTAINS
    
  SUBROUTINE GROWSTRING(RS, RF, INPUTIM, FINALX, FINALE, FINALTAN,RMS,MFLAG)
    ! RS and RF are start and end coordinates, respectively
    ! Implement growing string method using EVOLVESTRING subroutine
    ! FINALTAN is a 2d array pointer with coordinates as first dim, image as second

    USE INTCUTILS, ONLY : CART2INT, INTINTERPOLATE
    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(INOUT) :: RS(3*NATOMS), RF(3*NATOMS)
    DOUBLE PRECISION, INTENT(OUT) :: RMS
    INTEGER, INTENT(IN) :: INPUTIM
    DOUBLE PRECISION, POINTER :: FINALX(:), FINALE(:), FINALTAN(:,:)
    LOGICAL, INTENT(OUT) :: MFLAG

    TYPE(IMGNODE), POINTER :: DUMMYP    

    DOUBLE PRECISION :: PREVBACKTCUTOFF

    LOGICAL :: GROWLEFT
    ! make sure diff is big enough to hold either cartesians or internals
    DOUBLE PRECISION :: TMPRMS, PARAM, DIFF(3*NATOMS+NINTC), TOL, DIST
    DOUBLE PRECISION :: TMPXYZ(1:3*NATOMS*INPUTIM)
    INTEGER :: IM, CRD

    DOUBLE PRECISION :: STPINT(NINTC), STPCART(3*NATOMS)
    LOGICAL :: FAILED

    INTEGER :: ISTAT

    ! alignment stuff
    DOUBLE PRECISION :: DISTF, DIST2, RMAT(3,3)
    CHARACTER(LEN=5) :: ZSYMSAVE
    COMMON /SYS/ ZSYMSAVE

    ! testing only
    DOUBLE PRECISION :: TESTINT(NINTC)

    LOGICAL::REGENE

    ! --------------------------------------------    

    IF (.NOT.EVOLVESTRINGT.AND.INTINTERPT) THEN
       print*, 'ERROR! Growing strings not implemented with INTINTERP'
       STOP
    ENDIF

    PTEST = DESMDEBUG.OR.DEBUG
    INTPTEST = .FALSE.

    IF (PTEST) THEN
       CALL DUMPCOORDS(RS, 'start.xyz')
       CALL DUMPCOORDS(RF, 'fin.xyz')
    ENDIF

    ! Initialize some stuff; TIM & NC are declared in GSDATA
    TIM = INPUTIM ! total number of images
    IF (DESMINT) THEN ! number of coordinates
       NC = NINTC 
    ELSE
       NC = 3*NATOMS
    ENDIF
    NLIM = 0; NRIM = 0; NIM = 0
    M = GSUPDATE
    MAXLEN = MAXLENPERIM*TIM
    IF (GSMAXTOTITD.GE.0) MAXTOTSTEPS = GSMAXTOTITD*TIM
    TOTSTEPS = 0

    ! zero coordinates 1:4 and 7:8 of both start and end structures      
    IF (FIXATMS.AND.PREROTATE.AND..NOT.DESMINT) THEN
       ! shift molecule so 1st atm is at origin
       CALL SHIFTZERO(RS, NATOMS, 1)
       CALL SHIFTZERO(RF, NATOMS, 1)

       ! rotate arount y axis to zero x component of 3rd atom       
       CALL ROTATEZERO(RS, NATOMS, 3, 2, 1)
       CALL ROTATEZERO(RF, NATOMS, 3, 2, 1)

       ! rotate arount x axis to zero y component of 3rd atom
       CALL ROTATEZERO(RS, NATOMS, 3, 1, 2)
       CALL ROTATEZERO(RF, NATOMS, 3, 1, 2)

       ! now 3rd atom is on z axis
       ! rotate around z axis to zero x component of 2nd atom
       CALL ROTATEZERO(RS, NATOMS, 2, 3, 1)
       CALL ROTATEZERO(RF, NATOMS, 2, 3, 1)
    ENDIF

    ! make the endpoints
    NULLIFY(FIRST,LAST, LEFTFRONT, RIGHTFRONT)
    CALL NEWNODE(FIRST); CALL NEWNODE(LAST)        
    LEFTFRONT => FIRST; RIGHTFRONT => LAST

    IF (DESMINT) THEN
       DIHINFO(:,:) = 0.0D0
       
       FIRST%XYZCART(1:3*NATOMS) = RS(:)
       LAST%XYZCART(1:3*NATOMS) = RF(:)

       PREVDIH => DIHINFO(1,:)
       CALL CART2INT(FIRST%XYZCART(1:3*NATOMS),FIRST%XYZ(1:NC))
       
       !align to first image
       ALIGNDIR = .TRUE.
       DIHINFO(TIM+2,:) = DIHINFO(1,:)
       PREVDIH => DIHINFO(TIM+2,:)
       
       CALL CART2INT(LAST%XYZCART(1:3*NATOMS),LAST%XYZ(1:NC))                                          

       ALIGNDIR = .FALSE.
       
       FIRST%PREVXYZ(:) = FIRST%XYZ(:); FIRST%PREVXCART(:) = FIRST%XYZCART(:)
       LAST%PREVXYZ(:) = LAST%XYZ(:); LAST%PREVXCART(:) = LAST%XYZCART(:)
    ELSE
       FIRST%XYZ(1:NC) = RS(1:NC); FIRST%PREVXYZ(1:NC) = RS(1:NC)    
       LAST%XYZ(1:NC) = RF(1:NC); LAST%PREVXYZ(1:NC) = RF(1:NC)
    ENDIF

    DIFF(1:NC) = LAST%XYZ(1:NC) - FIRST%XYZ(1:NC)

    STRINGLEN = SQRT(DOT_PRODUCT(DIFF(1:NC),DIFF(1:NC)))
    IF(PTEST) print*, 'Starting string len: ', STRINGLEN

    FIRST%ARC = 0.0D0; FIRST%CHORD = 0.0D0; FIRST%IND = 0          

    ! get initial energies
    IF (DESMINT) THEN
       CALL POTENTIAL(RS(:), FIRST%E, FIRST%GCART, .FALSE., .FALSE., &
            & TMPRMS, .FALSE., .FALSE.)    
       CALL POTENTIAL(RF(:), LAST%E, LAST%GCART, .FALSE., .FALSE., &
            & TMPRMS, .FALSE., .FALSE.)            
    ELSE
       CALL POTENTIAL(RS(:), FIRST%E, FIRST%GRAD, .FALSE., .FALSE., &
            & TMPRMS, .FALSE., .FALSE.)    
       CALL POTENTIAL(RF(:), LAST%E, LAST%GRAD, .FALSE., .FALSE., &
            & TMPRMS, .FALSE., .FALSE.)
    ENDIF

    IF (CUBSPLT.OR.TANTYPE.EQ.3) THEN
       IF (DESMINT) THEN
          print*, 'ERROR! DESMINT not implemented with cubic splines yet!'
          STOP
       ENDIF

       ALLOCATE(COEFF(4,TIM+2,NC), ABSC(TIM+2))
       ABSC(1) = 0.0D0; COEFF(1,1,:) = FIRST%XYZ(:)
    ENDIF

    IF (EVOLVESTRINGT) THEN       

       ! for just evolving a linear interpolation       
       NIM = TIM
       NLIM = TIM/2; NRIM = TIM-TIM/2
       IF(CLSTRINGT) THEN
         NIM=TIM-1
         NLIM=NIM
         NRIM=0
       ENDIF
       JOINED = .TRUE.
       DO IM = 1,NLIM
          CALL EXTENDLEFT
       ENDDO
       DO IM = 1,NRIM
          CALL EXTENDRIGHT
       ENDDO

       ! move first and last to images rather than endpoints
       FIRST => FIRST%NEXT
       IF(.NOT.CLSTRINGT) THEN
          LAST => LAST%PREV 
       ENDIF
       REGENE=.TRUE.
       DO WHILE (REGENE)
          REGENE=.FALSE.
          IF (INTINTERPT) THEN
             CALL INTINTERPOLATE(RS,RF,NINTIM,NIM,TMPXYZ(:),PTEST, FAILED)
             DUMMYP => FIRST
             DO IM = 1,NIM
                DUMMYP%XYZ(:) = TMPXYZ(NC*(IM-1)+1:NC*IM)
                DUMMYP => DUMMYP%NEXT
             ENDDO      
          ELSE
             BACKTCUTOFF = INTERPBACKTCUT
             DUMMYP => FIRST
             DO IM = 1,NIM
                DUMMYP%XYZ = DUMMYP%PREV%XYZ + DIFF(1:NC)*1.0D0/(NIM+1)
                IF (DESMINT) THEN
                   DIHINFO(IM+1,:) = DIHINFO(IM,:)
                   PREVDIH => DIHINFO(IM+1,:)

                   STPINT(:) = DIFF(1:NC)*1.0D0/(NIM+1)
                   CALL TRANSBACKDELTA(STPINT,STPCART,DUMMYP%PREV%XYZCART,NINTC,3*NATOMS,NNZ,KD,FAILED,INTPTEST,INTEPSILON)
                   DUMMYP%XYZCART(:) = DUMMYP%PREV%XYZCART(:) + STPCART(:)
                   
                   ! align to start struct
                   CALL NEWMINDIST(RF(:),DUMMYP%XYZCART,NATOMS,DIST,.FALSE.,.FALSE.,ZSYM(1),.FALSE.,.FALSE.,.FALSE.,RMAT)
                ENDIF             

                DUMMYP => DUMMYP%NEXT

             ENDDO
             BACKTCUTOFF = MINBACKTCUT
          ENDIF
          CALL CHECKENERGIES(REGENE) ! check energies aren't infinite
          IF(REGENE) THEN
             RF=LAST%XYZ(1:NC)
             IF(DESMINT) RF=LAST%XYZCART(1:NC)
             DIFF(1:NC) = LAST%XYZ(1:NC) - FIRST%PREV%XYZ(1:NC)
             STRINGLEN = SQRT(DOT_PRODUCT(DIFF(1:NC),DIFF(1:NC)))
             IF(PTEST) print*, 'Recalc starting string len: ', STRINGLEN
          ENDIF
       ENDDO
       IF (PTEST) CALL TSUMMARY

       CALL EVOLVESTRING(GROWLEFT, RMS, MFLAG)      

    ELSE

       JOINED = .FALSE.

       DO WHILE(.NOT.JOINED)
          IF (NIM.EQ.0) THEN

             ! make first images          
             CALL EXTENDLEFT; CALL EXTENDRIGHT

             FIRST => LEFTFRONT; LAST => RIGHTFRONT 
             
             FIRST%XYZ(1:NC) = FIRST%PREV%XYZ(:) + 1.0D0/(TIM+1)*DIFF(1:NC)
             LAST%XYZ(1:NC) = LAST%NEXT%XYZ(:) - 1.0D0/(TIM+1)*DIFF(1:NC)

             IF (DESMINT) THEN
                STPINT(:) = DIFF(1:NC)*1.0D0/(TIM+1)

                DIHINFO(2,:) = DIHINFO(1,:)
                PREVDIH => DIHINFO(2,:)
                
                CALL TRANSBACKDELTA(STPINT,STPCART,FIRST%PREV%XYZCART,NINTC,3*NATOMS,NNZ,KD,FAILED,INTPTEST,INTEPSILON)
                FIRST%XYZCART(:) = FIRST%PREV%XYZCART(:) + STPCART(:)

                DIHINFO(TIM+1,:) = DIHINFO(TIM+2,:)
                PREVDIH => DIHINFO(2,:)
                
                CALL TRANSBACKDELTA(-STPINT,STPCART,FIRST%PREV%XYZCART,NINTC,3*NATOMS,NNZ,KD,FAILED,INTPTEST,INTEPSILON)
                LAST%XYZCART(:) = LAST%NEXT%XYZCART(:) + STPCART(:)
             ENDIF             
             NLIM = 1; NRIM = 1; NIM = 2

             LEFTFRONT%LF = .TRUE.; RIGHTFRONT%RF = .TRUE.             
             
             DUMMYP => FIRST
             DO IM = 1,2
                CALL GETIMGPOT(DUMMYP, TMPRMS)
                DUMMYP => DUMMYP%NEXT
             ENDDO
          ELSE
             ! Assume string interpolated, so diffs, chords, arcs are set
             IF (GROWLEFT) THEN ! grow from the left

                IF (PTEST) print*, 'GROWING LEFT'             
                CALL EXTENDLEFT

                IF(CUBSPLT) THEN
                   PARAM = FOLLOWARCNEWT((RIGHTFRONT%ARC-LEFTFRONT%PREV%ARC)*1.0/(TIM-NIM+1), NLIM, 0)
                   CALL GETSPLVAL(PARAM, NLIM, 0, LEFTFRONT%XYZ)                
                ELSE
                   LEFTFRONT%XYZ = LEFTFRONT%PREV%XYZ + RIGHTFRONT%DIFF * 1.0/(TIM-NIM+1)
                   IF (DESMINT) THEN
                      DIHINFO(NLIM+2,:) = DIHINFO(NLIM+1,:)
                      PREVDIH => DIHINFO(NLIM+2,:)

                      STPINT(:) = RIGHTFRONT%DIFF * 1.0/(TIM-NIM+1)
                      CALL TRANSBACKDELTA(STPINT,STPCART,LEFTFRONT%PREV%XYZCART,NINTC, &
                           & 3*NATOMS,NNZ,KD,FAILED,INTPTEST,INTEPSILON)
                      LEFTFRONT%XYZCART(:) = LEFTFRONT%PREV%XYZCART(:) + STPCART(:)
                   ENDIF
                ENDIF
                NLIM = NLIM + 1; NIM = NIM+1
                
                CALL GETIMGPOT(LEFTFRONT, TMPRMS)

             ELSE
                IF (PTEST) print*, 'GROWING RIGHT'

                CALL EXTENDRIGHT

                IF(CUBSPLT) THEN
                   PARAM = FOLLOWARCNEWT((RIGHTFRONT%NEXT%ARC-LEFTFRONT%ARC)*DBLE(TIM-NIM)/(TIM-NIM+1), NLIM+1, 1)
                   CALL GETSPLVAL(PARAM, NLIM, 0, RIGHTFRONT%XYZ)
                ELSE                   
                   RIGHTFRONT%XYZ = RIGHTFRONT%NEXT%XYZ - RIGHTFRONT%NEXT%DIFF*1.0/(TIM-NIM+1)
                   IF (DESMINT) THEN
                      DIHINFO(TIM-NRIM+1,:) = DIHINFO(TIM-NRIM+2,:)
                      PREVDIH => DIHINFO(TIM-NRIM+1,:)
                      
                      STPINT(:) = - RIGHTFRONT%NEXT%DIFF*1.0/(TIM-NIM+1)
                      CALL TRANSBACKDELTA(STPINT,STPCART,RIGHTFRONT%NEXT%XYZCART,NINTC, &
                           & 3*NATOMS,NNZ,KD,FAILED,INTPTEST,INTEPSILON)
                      RIGHTFRONT%XYZCART(:) = RIGHTFRONT%NEXT%XYZCART(:) + STPCART(:)
                   ENDIF
                ENDIF
                NRIM = NRIM +1; NIM = NIM+1
                
                CALL GETIMGPOT(RIGHTFRONT, TMPRMS)
                
             ENDIF
          ENDIF

          IF (NIM.EQ.TIM) THEN
             IF (PTEST) print*, 'GS> String has joined!'
             JOINED = .TRUE.
          ENDIF

          CALL EVOLVESTRING(GROWLEFT, RMS, MFLAG)
          
          IF (GSMAXTOTITD.GE.0.AND.TOTSTEPS.GT.MAXTOTSTEPS) EXIT
       ENDDO
    ENDIF

    IF(PTEST) print '(A,I6,2G20.10)', 'Final string length, len per image: ', TIM, STRINGLEN, STRINGLEN/(TIM+1)

    ! output results

    CALL OUTPUT(FINALX, FINALE, FINALTAN)

    ! if for some reason outputting before full string grown:
    IF (NIM.LT.TIM.AND..NOT.CLSTRINGT) THEN 
       DO IM = NIM+1, TIM
          FINALE(IM) = LAST%NEXT%E
          IF (DESMINT) THEN
             FINALX(NC*(IM-1)+1:NC*IM) = LAST%NEXT%XYZCART
          ELSE
             FINALX(NC*(IM-1)+1:NC*IM) = LAST%NEXT%XYZ
          ENDIF
       ENDDO
    ENDIF

    CALL DELETESTRING
    
    IF (CUBSPLT.OR.TANTYPE.EQ.3) DEALLOCATE(COEFF, ABSC)
    IF (CLSTRINGTST) CALL INEXACT_NEWTON(FINALX(NC*(TIM-1)+1:NC*TIM),FINALTAN(1:NC,TIM))
  END SUBROUTINE GROWSTRING

  SUBROUTINE EVOLVESTRING(GROWLEFT, FNRMS, MFLAG)
    ! Evolve string until either one of the frontier rms or the entire rms fperp is below a given tolerance
    ! lots of the variables are defined in gsdata
    ! returns growleft true if left side converged before right
    ! returns mflag true if successfully converged

    IMPLICIT NONE

    TYPE(IMGNODE), POINTER :: START, DUMMYP
    DOUBLE PRECISION, INTENT(OUT) :: FNRMS
    LOGICAL, INTENT(OUT) :: MFLAG, GROWLEFT
    INTEGER :: I, IM, IMCOUNT, IMT,IMC,NIMC
    LOGICAL :: REPARAM

    INTEGER :: BOUND, EVOLIT, P, K, POINT
    DOUBLE PRECISION :: SLENGTH, GAMMA, YY, SY, BETA, FNORM2, RMS
    DOUBLE PRECISION, DIMENSION(0:GSUPDATE-1,NC*TIM):: Y, S
    DOUBLE PRECISION, DIMENSION(0:GSUPDATE-1) :: RO, ALPHA
    DOUBLE PRECISION, DIMENSION(NC*TIM) :: FPERP, NRDIR, HF, DIAG, Q, STP, PREVFPERP

    DOUBLE PRECISION :: FNRMSL, FNRMSR
    DOUBLE PRECISION, DIMENSION(NC*TIM) :: FULLGRAD, PREVFULLGRAD

    DOUBLE PRECISION :: AVGE, MAXE, PREVE
    INTEGER :: ISTAT
    
    ! benchmarks only
    DOUBLE PRECISION :: TOTE

    ! for internals
    LOGICAL :: FAILED,REGENE
    DOUBLE PRECISION :: STPINT(NC), STPCART(3*NATOMS)

    ! testing only
    INTEGER :: CRD

    ! -------------------------------------------    

    ! do some initializing stuff
    MFLAG = .FALSE.

    !---------
    ! Evolution step, using LBFGS method as described in Nocedal & Wright, 2000, Ch.9
    !---------
    EVOLIT = 0 ! current string evolution iteration (the index k in Peters, et al)
    K = 0
    NIMC=NIM
    IF(CLSTRINGT) NIMC=TIM
    IF (HESSGRAD) THEN
       DUMMYP => FIRST
       DO IM = 1,NIMC
          FULLGRAD(NC*(IM-1)+1:NC*IM) = DUMMYP%GRAD(:)
          DUMMYP=>DUMMYP%NEXT
       ENDDO
    ENDIF
    REGENE=.FALSE.
    DO         
       EVOLIT = EVOLIT + 1
       TOTSTEPS = TOTSTEPS + 1

       CALL INTERPOLATE          ! Interpolate string, reset DIFFS, ARCS, CHORDS       
       IF(.NOT.REGENE) THEN
          CALL DECIDEREPARAM(REPARAM) ! decide whether to reparameterize       
       ELSE
          REPARAM=.TRUE.
          REGENE=.FALSE.
          print*, "String truncated"
       ENDIF


       IF (REPARAM) CALL REPARAMETRIZE          

       IMCOUNT = NIM ! imcount is the number of movable images
       START => FIRST ! start of movable images
       
       ! get energies and grad vectors
       ! if not reparametrizing, can just use energies and grads from end of last cycle
       ! assumes energies were precalculated before entering EVOLVESTRING routine
       IF (REPARAM) THEN
          DUMMYP => START
          DO IM = 1,IMCOUNT
             CALL GETIMGPOT(DUMMYP,RMS)
             IF (HESSGRAD) FULLGRAD(NC*(IM-1)+1:NC*IM) = DUMMYP%GRAD(:)
             DUMMYP => DUMMYP%NEXT
          ENDDO
          ! reset Hessian approximation
          K = 0
       ENDIF

       CALL GETTANGENTS
       FNRMSL = 0.0D0; FNRMSR = 0.0D0
       FNORM2 = 0.0D0       

       DUMMYP => START
       DO IM = 1,NIMC
          ! recalculate fperp and fnorms
          FPERP(NC*(IM-1)+1:NC*IM) = - DUMMYP%GRAD(:) + &
               & DOT_PRODUCT(DUMMYP%GRAD, DUMMYP%TGT)*DUMMYP%TGT(:)       
          IF(IM==TIM.AND.CLSTRINGT) FPERP(NC*(IM-1)+1:NC*IM) = FPERP(NC*(IM-1)+1:NC*IM) + &
                      DOT_PRODUCT(DUMMYP%GRAD, DUMMYP%TGT)*DUMMYP%TGT(:)

          IF (FIXATMS) THEN
             FPERP(NC*(IM-1)+1:NC*(IM-1)+4) = 0.0D0
             FPERP(NC*(IM-1)+7:NC*(IM-1)+8) = 0.0D0
          ENDIF

          DUMMYP%FNORM2 = DOT_PRODUCT(FPERP(NC*(IM-1)+1:NC*IM),FPERP(NC*(IM-1)+1:NC*IM))

          FNORM2 = FNORM2+DUMMYP%FNORM2
          IF (JOINED.OR.DUMMYP%LF) THEN 
             FNRMSL = FNRMSL+DUMMYP%FNORM2                         
          ELSE IF (DUMMYP%RF) THEN
             FNRMSR = FNRMSR + DUMMYP%FNORM2
          ENDIF
          DUMMYP => DUMMYP%NEXT
       ENDDO

       IF(CLSTRINGT) THEN
          FNRMS = SQRT(FNORM2/IMCOUNT+1)
       ELSE
          FNRMS = SQRT(FNORM2/IMCOUNT)
       ENDIF

       ! decide whether to freeze some nodes b/c not contributing significantly to FNRMS
       IF (FREEZENODEST) THEN
          DUMMYP => START
          DO IM = 1,IMCOUNT
             DUMMYP%FREEZE = (SQRT(DUMMYP%FNORM2)/FNRMS.LT.FREEZETOL)
             IF (DUMMYP%FREEZE.AND.PTEST) &
                  & print '(A,I4,L2,2F15.5)', 'Freezing image: ', IM, DUMMYP%FREEZE,SQRT(DUMMYP%FNORM2),FNRMS
             DUMMYP => DUMMYP%NEXT
          ENDDO
       ENDIF

       FNRMS = FNRMS / SQRT(DBLE(NC))

       IF (JOINED) THEN
          FNRMSL = SQRT(FNRMSL/NIMC/NC)
       ELSE
          FNRMSL = SQRT(FNRMSL/NC)
          FNRMSR = SQRT(FNRMSR/NC)
       ENDIF
       IF (PTEST) THEN
          ! get average and maximum energy
          AVGE = 0.0D0
          MAXE = -1.0D10
          DUMMYP => FIRST
          DO IM = 1,NIMC
             AVGE = AVGE + DUMMYP%E
             IF (DUMMYP%E.GT.MAXE) MAXE = DUMMYP%E
             DUMMYP => DUMMYP%NEXT
          ENDDO
          AVGE = AVGE/NIMC

          IF (DUMPGSALL) THEN
             CALL DUMPGSPATH(EVOLIT)
          ELSE IF (PTEST) THEN
             CALL DUMPGSPATH(-1)
          ENDIF
          IF (JOINED) THEN
             print '(A,1x,3G20.10,1x,2I6)', '>>> FNRMS, MAXE, AVGE, NIM, EVOLIT: ', FNRMS, MAXE, AVGE, NIM, EVOLIT
          ELSE
             print '(A,1x,3G20.10,1x,2I6)', '>>> FNRMSL, FNRMSR, FNRMS, NIM, EVOLIT: ', FNRMSL, FNRMSR, FNRMS, NIM, EVOLIT
          ENDIF

          CALL FLUSH(6)
       ENDIF

       ! check for convergence
       IF (((.NOT.CLSTRINGT).OR.CLSTRINGTST).AND.JOINED.AND.FNRMSL.LT.GSCONV) THEN
          IF (PTEST) print*, 'SUCCESS! - whole string converged!'
          MFLAG = .TRUE.
          EXIT
       ELSEIF (CLSTRINGT.AND.FNRMSL.LT.GSCONV*10d0.AND.SQRT(LAST%FNORM2).LT.GSCONV*SQRT(DBLE(NC))) THEN
          IF (PTEST) print*, 'SUCCESS! - the last image converged!'
          MFLAG = .TRUE.
          EXIT
       ELSE IF (.NOT.JOINED.AND.FNRMSL.LT.GSGROWTOL) THEN
          IF (PTEST) print*, 'SUCCESS! - left end converged!'
          GROWLEFT = .TRUE.
          MFLAG = .TRUE.
          EXIT 
       ELSE IF (.NOT.JOINED.AND.FNRMSR.LT.GSGROWTOL) THEN
          IF (PTEST) print*, 'SUCCESS! - right end converged!'
          MFLAG = .TRUE.
          GROWLEFT = .FALSE.
          EXIT  
       ENDIF

       IF (GSMAXTOTITD.GE.0.AND.TOTSTEPS.GT.MAXTOTSTEPS) THEN
          print*, 'Total number of steps for entire string exceeded maximum.', NIM, TOTSTEPS, MAXTOTSTEPS
          IF(PTEST) print*, 'String length: ', STRINGLEN
          EXIT
       ELSE IF(JOINED.AND.EVOLIT.GT.GSITERDENSITY*NIMC.AND.AVGE.LT.DESMAXAVGE) THEN
          IF (PTEST) THEN 
             print*, 'Number of final evolution steps exceeded maximum; average energy below MAXAVGE.'
             IF (PTEST) print*, 'String length: ', STRINGLEN
          ENDIF
          GROWLEFT = (FNRMSL <= FNRMSR)
          EXIT
       ELSE IF (.NOT.JOINED.AND.EVOLIT.GT.MAXGROWSTEPS) THEN
          print*, 'Number of growing evolution steps exceeded maximum.'
          IF (PTEST) print*, 'String length: ', STRINGLEN
          GROWLEFT = (FNRMSL <= FNRMSR)
          EXIT
       ENDIF


       ! get the NR steps and move images       
       lbfgs: IF (.NOT.NOLBFGS) THEN
       main: IF(K.EQ.0) THEN
          ! Use DGUESS as the initial guess for the inverse Hessian diagonal
          HF(:) = GSDGUESS*FPERP(:)
             
          ! Make the first guess for the step length cautious
          IF (FNORM2.EQ.0.0D0) THEN
             FNORM2=1.0D0 ! exact zero is presumably wrong!
             PRINT '(A)','WARNING - FNORM was zero in lbfgs,&
                  & resetting to one'
          ENDIF
          STP(:) = MIN(1.0D0/SQRT(FNORM2), SQRT(FNORM2))
          
          POINT = 0
       
       ELSE main
          P = MODULO(POINT-1,M) ! position in lbfgs lists
          BOUND = MIN(K,M) ! how far backwards to look for Hessian  
          
          IF (HESSGRAD) THEN
             Y(P,:) = FULLGRAD(:) - PREVFULLGRAD(:)
          ELSE
             Y(P,:) = PREVFPERP(:) - FPERP(:)
          ENDIF
        
          ! Get initial diagonal Hessian estimate; eq. 9.6 in Nocedal, et al.
          SY = DOT_PRODUCT(S(P,:),Y(P,:))               
          YY = DOT_PRODUCT(Y(P,:), Y(P,:))               
          IF (SY.EQ.0.0D0) SY = 1.0D0
          IF (YY.EQ.0.0D0) YY = 1.0D0
          GAMMA = SY / YY
          DIAG(:) = GAMMA

          ! update RO
          RO(P) = 1.0D0/SY
             
          ! Get HF = Hk*fperp (two-loop algorithm 9.1 in Nocedal et al)
          Q(:) = FPERP(:)             
                    
          DO I = 1,BOUND                             
             ALPHA(P) = RO(P)*DOT_PRODUCT(S(P,:),Q(:))
             Q(:) = Q(:) - ALPHA(P)*Y(P,:)
             P = MODULO(P-1,M)
          ENDDO          

          HF(:) = DIAG(:)*Q(:)

          DO I = 1, BOUND
             P = MODULO(P+1,M)
             BETA = RO(P)*DOT_PRODUCT(Y(P,:), HF(:))
             HF(:) = HF(:) + (ALPHA(P) - BETA)*S(P,:)                
          ENDDO
          STP(:) = 1.0
       ENDIF main
       ENDIF lbfgs

       IF (HESSGRAD) THEN
          DUMMYP => START
          DO IM = 1,NIM
             NRDIR(NC*(IM-1)+1:NC*IM) = HF(NC*(IM-1)+1:NC*IM) - &
                  & DOT_PRODUCT(HF(NC*(IM-1)+1:NC*IM), DUMMYP%TGT)*DUMMYP%TGT(:)
             DUMMYP => DUMMYP%NEXT
          ENDDO
          IF(CLSTRINGT) THEN
             IM=TIM
             NRDIR(NC*(IM-1)+1:NC*IM) = HF(NC*(IM-1)+1:NC*IM) - 2D0*&
                      & DOT_PRODUCT(HF(NC*(IM-1)+1:NC*IM),DUMMYP%TGT)*DUMMYP%TGT(:)         
          ENDIF
       ELSE IF (NOLBFGS) THEN
          NRDIR(:) = FPERP(:)
          STP(:) = 1.0D0
       ELSE
          NRDIR(:) = HF(:)
       ENDIF

       IF (FIXATMS) THEN
          DO IM = 1,NIMC
             NRDIR(NC*(IM-1)+1:NC*(IM-1)+4) = 0.0D0
             NRDIR(NC*(IM-1)+7:NC*(IM-1)+8) = 0.0D0
          ENDDO
       ENDIF
       
       IF (DOT_PRODUCT(NRDIR,FPERP).LT.0) THEN 
          IF(PTEST) print*, 'NRDIR has negative projection onto fperp. Reversing step.'
          NRDIR(:) = -NRDIR(:)
       ENDIF

       ! take no more than max step for each individual image
       IF (HESSGRAD) PREVFULLGRAD(:) = FULLGRAD(:)

       DUMMYP => START
       DO IM = 1,NIMC
          IF (DUMMYP%FREEZE) THEN
             STP(NC*(IM-1)+1:NC*IM) = 0.0D0
             S(POINT,NC*(IM-1)+1:NC*IM) = 0.0D0
          ELSE
             DUMMYP%PREVXYZ(1:NC) = DUMMYP%XYZ(1:NC)
             IF (DESMINT) DUMMYP%PREVXCART(1:3*NATOMS) = DUMMYP%XYZCART(1:3*NATOMS)

             SLENGTH = SQRT(DOT_PRODUCT(NRDIR(NC*(IM-1)+1:NC*IM), NRDIR(NC*(IM-1)+1:NC*IM)))
             STP(NC*(IM-1)+1:NC*IM) = MIN(GSMXSTP/SLENGTH, STP(NC*(IM-1)+1))
             
             S(POINT,NC*(IM-1)+1:NC*IM) = STP(NC*(IM-1)+1:NC*IM)*NRDIR(NC*(IM-1)+1:NC*IM)
             DUMMYP%XYZ(:) = DUMMYP%XYZ(:) + S(POINT,NC*(IM-1)+1:NC*IM)   

             IF (DESMINT) THEN
                PREVDIH => DIHINFO(IM+1,:)
                CALL TRANSBACKDELTA(S(POINT,NC*(IM-1)+1:NC*IM),STPCART,DUMMYP%XYZCART,NINTC, &
                     & 3*NATOMS,NNZ,KD,FAILED,INTPTEST,INTEPSILON)                
                DUMMYP%XYZCART = DUMMYP%XYZCART + STPCART
             ENDIF

             PREVE = DUMMYP%E
             CALL GETIMGPOT(DUMMYP,RMS)

             IF ((PREVE.LT.0.AND.DUMMYP%E-PREVE.GT.DESMAXEJUMP).OR.(DESMINT.AND.FAILED)) THEN
                IF (PTEST) print*, 'Energy jump too big for image or failed transbackdelta. Skipping this image.', &
                     & IM, PREVE, DUMMYP%E
                DUMMYP%XYZ(:) = DUMMYP%PREVXYZ(:)
                IF (DESMINT) DUMMYP%XYZCART(:) = DUMMYP%PREVXCART(:)
                STP(NC*(IM-1)+1:NC*IM) = 0.0D0
                S(POINT,NC*(IM-1)+1:NC*IM) = 0.0D0
                CALL GETIMGPOT(DUMMYP,RMS)
             ENDIF
             IF (CLSTRINGT) THEN
                IF (DUMMYP%E<DUMMYP%PREV%E) THEN
                   IMC=IM
                   IF(IM==1) THEN
                      IMC=2
                      DUMMYP => DUMMYP%NEXT
                   ENDIF
                   DO IMT=1,TIM-IMC+1
                      DUMMYP%XYZ=DUMMYP%PREV%XYZ
                      DUMMYP%E=DUMMYP%PREV%E
                      DUMMYP%GRAD=DUMMYP%PREV%GRAD
                      DUMMYP%PREVXYZ=DUMMYP%PREV%PREVXYZ
                      IF(DESMINT) THEN
                         DUMMYP%XYZCART = DUMMYP%PREV%XYZCART
                         DUMMYP%GCART = DUMMYP%PREV%GCART
                         DUMMYP%PREVXCART = DUMMYP%PREV%PREVXCART
                      ENDIF
                      DUMMYP => DUMMYP%NEXT
                   ENDDO
                   REGENE=.TRUE.
                   EXIT
                ENDIF
             ENDIF

             IF (HESSGRAD) FULLGRAD(NC*(IM-1)+1:NC*IM) = DUMMYP%GRAD
          ENDIF
          DUMMYP => DUMMYP%NEXT
       ENDDO

       IF (.NOT.HESSGRAD) PREVFPERP(:) = FPERP(:)       

       ! update iteration number and point in lbfgs lists
       K = K + 1
       POINT = MODULO(POINT+1,M)          

!       DUMMYP => DUMMYP%NEXT
       
    ENDDO

    IF (PTEST) CALL DUMPGSPATH(-1)

  END SUBROUTINE EVOLVESTRING

  SUBROUTINE DECIDEREPARAM(REPARAM)
    ! decide whether or not to reparametrize string
    
    IMPLICIT NONE
    
    LOGICAL, INTENT(OUT) :: REPARAM
    INTEGER :: IM,NIMC
    DOUBLE PRECISION :: OFFSET, GOODLEN, ACTLEN, NEWOFFSET, MINOFFSET, MAXOFFSET
    TYPE(IMGNODE), POINTER :: DUMMYP
    
    REPARAM = .FALSE.
    
    MINOFFSET = DBLE(NIM)
    MAXOFFSET = 0.0D0
    NIMC=TIM
    IF(CLSTRINGT) NIMC=NIM
    DUMMYP => FIRST%PREV
    DO IM = 1,NIM+1
       DUMMYP => DUMMYP%NEXT

       ! get desired length
       IF (.NOT.JOINED.AND.DUMMYP%RF) THEN 
          GOODLEN = STRINGLEN*(1.0-DBLE(NIM)/(NIMC+1))
       ELSE
          GOODLEN = STRINGLEN*1.0D0/(NIMC+1)
       ENDIF
       
       ! get actual length
       IF (CUBSPLT) THEN 
          ACTLEN = DUMMYP%ARC-DUMMYP%PREV%ARC
       ELSE
          ACTLEN = DUMMYP%CHORD
       ENDIF
       
       NEWOFFSET  = ACTLEN/GOODLEN
       IF (NEWOFFSET.LT.MINOFFSET) THEN
          MINOFFSET = NEWOFFSET
       ELSE IF (NEWOFFSET.GT.MAXOFFSET) THEN
          MAXOFFSET = NEWOFFSET
       ENDIF

    ENDDO       
    
    IF (MINOFFSET/MAXOFFSET.LT.REPARAMTOL) THEN
       IF(PTEST) print '(A,5F20.10)', ' minoffset, maxoffset, ratio, reparamtol, total length: ', &
  &              MINOFFSET, MAXOFFSET, MINOFFSET/MAXOFFSET, REPARAMTOL, STRINGLEN
       REPARAM = .TRUE.
    ENDIF    
  END SUBROUTINE DECIDEREPARAM

  SUBROUTINE REPARAMETRIZE
    ! reparametrize string
    ! assume chords, arcs, diffs, stringlen have all been set already

    IMPLICIT NONE

    INTEGER :: IM,NIMC
    TYPE(IMGNODE), POINTER :: CURN1, CURN2
    DOUBLE PRECISION :: PARAM, STPINT(NINTC), STPCART(3*NATOMS)
    LOGICAL :: FAILED

    ! alignment stuff
    DOUBLE PRECISION :: DISTF, DIST, DIST2, RMAT(3,3)
    CHARACTER(LEN=5) :: ZSYMSAVE
    COMMON /SYS/ ZSYMSAVE

    DOUBLE PRECISION :: STPDIST1, STPDIST2, STEP(NC)

    IF(PTEST) print*, 'reparametrizing string...'
    NIMC=TIM
    IF(CLSTRINGT) NIMC=NIM
    ! save current coordinates in prevxyz
    CURN1 => FIRST
    DO IM = 1, NIM       
       CURN1%PREVXYZ(1:NC) = CURN1%XYZ(1:NC)
       IF (DESMINT) CURN1%PREVXCART(1:3*NATOMS) = CURN1%XYZCART(1:3*NATOMS)
       CURN1%IND = IM
       CURN1 => CURN1%NEXT
    ENDDO
    IF(CLSTRINGT) THEN
       LAST%IND = NIM+1
    ELSE
       LAST%NEXT%IND = NIM+1
    ENDIF
    CURN1 => FIRST; CURN2 => FIRST

    DO IM = 1,NLIM
       CURN2 => CURN2%PREV

       DO WHILE (CURN2%ARC.LE.STRINGLEN*DBLE(IM)/(NIMC+1))
          CURN2 => CURN2%NEXT          
       ENDDO

       IF (CUBSPLT) THEN
          PARAM = FOLLOWARCNEWT(STRINGLEN*DBLE(IM)/(NIMC+1)-CURN2%PREV%ARC,CURN2%PREV%IND,0)
          CALL GETSPLVAL(PARAM,CURN2%PREV%IND,0,CURN1%XYZ)
       ELSE
          STEP(:) = CURN2%DIFF(:)/CURN2%CHORD *(STRINGLEN*DBLE(IM)/(NIMC+1)-CURN2%PREV%ARC)
          CURN1%XYZ(:) = CURN2%PREV%PREVXYZ(:) + STEP(:)
          IF (DESMINT) THEN
             STPDIST1 = SQRT(DOT_PRODUCT(STEP,STEP))
             STPINT(:) = CURN1%XYZ - CURN1%PREVXYZ
             STPDIST2 = SQRT(DOT_PRODUCT(STPINT,STPINT))

             IF (NATINT) THEN
                PREVDIH => DIHINFO(IM+1,:)
                PREVDIH(:) = 0.0D0
             ENDIF

             IF (STPDIST1.LT.STPDIST2) THEN 
                ! take step from previous image position
                CALL TRANSBACKDELTA(STEP,STPCART,CURN2%PREV%PREVXCART,NINTC,3*NATOMS,NNZ,KD,FAILED,INTPTEST,INTEPSILON)
                IF (FAILED) THEN
                   print*, 'TRANSBACKDELTA failed in REPARAMETRIZE w/ STPDIST1', STPDIST1, STPDIST2
                   STOP
                ENDIF
                CURN1%XYZCART = CURN2%PREV%PREVXCART + STPCART
             ELSE
                ! take step from this image position
                CALL TRANSBACKDELTA(STPINT,STPCART,CURN1%PREVXCART,NINTC,3*NATOMS,NNZ,KD,FAILED,INTPTEST,INTEPSILON)
                IF (FAILED) THEN
                   print*, 'TRANSBACKDELTA failed in REPARAMETRIZE w/ STPDIST2', STPDIST1, STPDIST2
                   STOP
                ENDIF
                CURN1%XYZCART = CURN1%PREVXCART + STPCART
             ENDIF
          ENDIF
       ENDIF       
       CURN1 => CURN1%NEXT
    ENDDO
    CURN1 => LAST
    IF(.NOT.CLSTRINGT) CURN2 => LAST%NEXT
    DO IM = 1,NRIM
       DO WHILE (CURN2%ARC.GE.STRINGLEN*(1.0-DBLE(IM)/(NIMC+1)))
          CURN2 => CURN2%PREV
       ENDDO
       CURN2 => CURN2%NEXT       
       IF (CUBSPLT) THEN
          PARAM = FOLLOWARCNEWT(STRINGLEN*(1.0D0-DBLE(IM)/(TIM+1))-CURN2%PREV%ARC,CURN2%IND,1)
          CALL GETSPLVAL(PARAM,CURN2%PREV%IND,0,CURN1%XYZ)
       ELSE          
          STEP(:) = - CURN2%DIFF(:)/CURN2%CHORD*(CURN2%ARC - STRINGLEN*(1.0-DBLE(IM)/(TIM+1)))
          CURN1%XYZ(:) = CURN2%PREVXYZ(:) + STEP(:)
          IF (DESMINT) THEN
             STPDIST1 = SQRT(DOT_PRODUCT(STEP,STEP))
             STPINT(:) = CURN1%XYZ - CURN1%PREVXYZ
             STPDIST2 = SQRT(DOT_PRODUCT(STPINT,STPINT))
                          
             IF (NATINT) THEN
                PREVDIH => DIHINFO(TIM-IM+2,:)
                PREVDIH(:) = 0.0D0
             ENDIF

             IF (STPDIST1.LT.STPDIST2) THEN
                CALL TRANSBACKDELTA(STEP,STPCART,CURN2%PREVXCART,NINTC,3*NATOMS,NNZ,KD,FAILED,INTPTEST,INTEPSILON)
                IF (FAILED) THEN
                   print*, 'TRANSBACKDELTA failed in REPARAMETRIZE w/ STPDIST1', STPDIST1, STPDIST2
                   STOP
                ENDIF
                CURN1%XYZCART = CURN2%PREVXCART + STPCART
             ELSE                
                CALL TRANSBACKDELTA(STPINT,STPCART,CURN1%PREVXCART,NINTC,3*NATOMS,NNZ,KD,FAILED,INTPTEST,INTEPSILON)
                IF (FAILED) THEN
                   print*, 'TRANSBACKDELTA failed in REPARAMETRIZE w/ STPDIST2', STPDIST1, STPDIST2
                   STOP
                ENDIF
                CURN1%XYZCART = CURN1%PREVXCART + STPCART
             ENDIF
          ENDIF
       ENDIF

       CURN1 => CURN1%PREV         
    ENDDO

    ! align cartesian structures along string
    IF (DESMINT) THEN
       CURN1 => FIRST%PREV
       DO IM = 1,NIM+1
          CURN1 => CURN1%NEXT
          CALL NEWMINDIST(FIRST%PREV%XYZCART,CURN1%XYZCART,NATOMS,DIST,.FALSE.,.FALSE.,ZSYM(1),.FALSE.,.FALSE.,.FALSE.,RMAT)
       ENDDO
    ENDIF

    CALL INTERPOLATE

    IF (PTEST) print*, 'String length, NIM after reparam: ', STRINGLEN, NIM
    
    RETURN
  END SUBROUTINE REPARAMETRIZE

  SUBROUTINE INTERPOLATE
    ! interpolate string
    !recalculate diffs, chords, arcs, stringlen
    ! if using internal coords, these are all in internals

    IMPLICIT NONE

    TYPE(IMGNODE), POINTER :: DUMMYP
    INTEGER :: IM  

    IF (CUBSPLT.OR.TANTYPE.EQ.3) THEN
       CALL INTERPCUBSPL
       RETURN
    ENDIF

    DUMMYP => FIRST%PREV
    DO IM = 1, NIM+1    
       DUMMYP => DUMMYP%NEXT

       DUMMYP%DIFF(:) = DUMMYP%XYZ(:) - DUMMYP%PREV%XYZ(:)       
       DUMMYP%CHORD = SQRT(DOT_PRODUCT(DUMMYP%DIFF, DUMMYP%DIFF))
       DUMMYP%ARC = DUMMYP%PREV%ARC + DUMMYP%CHORD

    ENDDO

    STRINGLEN = DUMMYP%ARC

    IF (STRINGLEN.GT.MAXLEN) THEN
       print*, 'STRINGLEN exceeds MAXLEN', STRINGLEN, MAXLEN
       STOP
    ENDIF

  END SUBROUTINE INTERPOLATE

  SUBROUTINE GETTANGENTS
    IMPLICIT NONE
    ! get tangent vectors; entire tangent vector for all images is normalized
    ! Assumes images, diffs, chords, arcs are all properly allocated and calculated

    INTEGER :: IM
    DOUBLE PRECISION :: WP, WM
    TYPE(IMGNODE), POINTER :: CURN

    CURN => FIRST
    DO IM = 1,NIM
       IF (CUBSPLT.OR.TANTYPE.EQ.3) THEN

          CALL GETSPLVAL(0.0D0, IM, 1, CURN%TGT)
          
       ELSE IF (TANTYPE.EQ.1) THEN
          IF (CURN%E > CURN%PREV%E .AND. CURN%NEXT%E > CURN%E) THEN
             WP = 1.0; WM = 0.0
          ELSE IF (CURN%E < CURN%PREV%E .AND. CURN%NEXT%E < CURN%E) THEN
             WP = 0.0; WM = 1.0
          ELSE
             WP = MAX(ABS(CURN%NEXT%E-CURN%E), ABS(CURN%PREV%E-CURN%E))
             WM = MIN(ABS(CURN%NEXT%E-CURN%E), ABS(CURN%PREV%E-CURN%E))
          ENDIF
          IF (WP.EQ.0.0D0.AND.WM.EQ.0.0D0) WP = 1.0D0
          CURN%TGT = WP*CURN%NEXT%DIFF(:) + WM*CURN%DIFF(:)
          
       ELSE IF (TANTYPE.EQ.2) THEN
          CURN%TGT = CURN%DIFF(:)/SQRT(DOT_PRODUCT(CURN%DIFF(:),CURN%DIFF(:))) + &
               & CURN%NEXT%DIFF(:)/SQRT(DOT_PRODUCT(CURN%NEXT%DIFF(:),CURN%NEXT%DIFF(:)))
       ELSE IF (TANTYPE.EQ.4) THEN
          CURN%TGT = CURN%NEXT%DIFF + CURN%DIFF
       ELSE
          print*, 'TANTYPE must be 1,2,3,or 4 if using linear interpolation: ', TANTYPE
          STOP
       ENDIF

       CURN%TNORM = SQRT(DOT_PRODUCT(CURN%TGT,CURN%TGT))
       
       IF (CURN%TNORM.EQ.0.0D0) THEN
          print*, 'ERROR: TNORM is zero. Images too close together. This is image: ', IM
          STOP
       ENDIF

       CURN%TGT = CURN%TGT/CURN%TNORM

       CURN => CURN%NEXT
    ENDDO
    IF(CLSTRINGT) THEN
       CURN%TGT = CURN%DIFF(:)
       CURN%TNORM = SQRT(DOT_PRODUCT(CURN%TGT,CURN%TGT))
       CURN%TGT = CURN%TGT/CURN%TNORM
    ENDIF
    
  END SUBROUTINE GETTANGENTS

  SUBROUTINE DUMPGSPATH(N)

    IMPLICIT NONE
    INTEGER, INTENT(IN) :: N
    INTEGER :: I,J
    TYPE(IMGNODE), POINTER :: DUMMYP
    CHARACTER*30 :: FNAME

    IF (N.LT.0) THEN
       FNAME = 'gspathway.xyz'
    ELSE
       WRITE(FNAME,'(A,I0.3,A)') 'gspathway.', N, '.xyz'
    ENDIF

    ! output pathway
    OPEN(UNIT=45,FILE=FNAME,STATUS='UNKNOWN')

    DUMMYP => FIRST%PREV
    DO I = 0,NIM+1
       WRITE(45,'(I6)') NATOMS
       WRITE(45,'(A,G25.15)') ' Energy= ', DUMMYP%E
       IF (DESMINT) THEN
          WRITE(45,'(A3,3G20.10)') ('LA ',DUMMYP%XYZCART(3*(J-1)+1:3*(J-1)+3),J=1,NATOMS)
       ELSE
          WRITE(45,'(A3,3G20.10)') ('LA ',DUMMYP%XYZ(3*(J-1)+1:3*(J-1)+3),J=1,NATOMS)
       ENDIF
       IF (I.EQ.NIM+1) EXIT
       DUMMYP => DUMMYP%NEXT
    ENDDO
    CLOSE(45)

    OPEN(UNIT=45,FILE='imgenergies.out', STATUS='UNKNOWN')
    DUMMYP => FIRST%PREV
    DO I=0,NIM+1
       WRITE(45,'(I6,G20.10)') I+1, DUMMYP%E
       IF (I.EQ.NIM+1) EXIT
       DUMMYP => DUMMYP%NEXT
    ENDDO
    CLOSE(45)

  END SUBROUTINE DUMPGSPATH    

  SUBROUTINE OUTPUT(FINALX, FINALE, FINALTAN)
    ! output the final coordinates and energy

    IMPLICIT NONE

    DOUBLE PRECISION, POINTER :: FINALX(:), FINALE(:), FINALTAN(:,:)
    INTEGER :: IM,NIMC
    TYPE(IMGNODE), POINTER :: DUMMYP
    DOUBLE PRECISION :: TNORM

    DUMMYP => FIRST
    NIMC=NIM
    IF(CLSTRINGT) NIMC=TIM
    DO IM = 1,NIMC
       FINALE(IM) = DUMMYP%E
       IF (DESMINT) THEN
          FINALX(3*NATOMS*(IM-1)+1:3*NATOMS*IM) = DUMMYP%XYZCART(:)
       ELSE
          FINALX(3*NATOMS*(IM-1)+1:3*NATOMS*IM) = DUMMYP%XYZ(:)
       ENDIF
       FINALTAN(1:NC, IM) = DUMMYP%TGT(:)
       DUMMYP => DUMMYP%NEXT
    ENDDO           
    RETURN
  END SUBROUTINE OUTPUT

  SUBROUTINE DUMPCOORDS(X, FNAME)
    IMPLICIT NONE
    CHARACTER (*) :: FNAME
    LOGICAL :: APPEND
    DOUBLE PRECISION :: X(3*NATOMS)
    INTEGER :: A

    ! given coordinate array X for a single molecule, dump into file FNAME
    OPEN (UNIT = 55, FILE = FNAME, STATUS = 'UNKNOWN')
    WRITE(55,'(I6)') NATOMS
    WRITE(55,'(A)') ' '
    WRITE(55,'(A3,3F24.15)') ('AX ',X(3*(A-1)+1:3*(A-1)+3),A=1,NATOMS)
    CLOSE(55)
  END SUBROUTINE DUMPCOORDS

  SUBROUTINE CHECKENERGIES(REGENE)
    ! calculate energies for points on the current string
    ! if any are too large, randomly perturb the images

    TYPE (IMGNODE), POINTER :: DUMMYP
    INTEGER :: IM, J
    DOUBLE PRECISION :: RMS, HARVEST, DPRAND
    LOGICAL::REGENE

    DUMMYP => FIRST
    DO IM = 1,NIM
       IF (DESMINT) PREVDIH => DIHINFO(IM+1,:)
       CALL GETIMGPOT(DUMMYP,RMS)
       IF (-DUMMYP%E.LT.-HUGE(DUMMYP%E)) THEN
          IF (DESMINT) THEN
             print*, 'ERROR: shouldnt have infinite energies when working in internals'
             STOP
          ENDIF

          ! energy is NaN or infinite
          PRINT *, "IMAGE",IM," IS BAD! - trying to lower it's energy..."
          DO J=1,3*NATOMS ! CHANGING GEOMETRY RANDOMLY
             HARVEST=DPRAND()
             DUMMYP%XYZ(J) = DUMMYP%XYZ(J) + HARVEST*0.01
          ENDDO

          CALL GETIMGPOT(DUMMYP, RMS)

          IF (-DUMMYP%E.LT.-HUGE(DUMMYP%E)) THEN
             print*, 'FAILED!'
             CALL TSUMMARY
             STOP
          ENDIF
       ENDIF
       IF (CLSTRINGT) THEN
          IF (DUMMYP%E<DUMMYP%PREV%E) THEN
             REGENE=.TRUE.
             IF(IM==1) THEN
                DUMMYP => DUMMYP%NEXT
             ENDIF
             LAST%XYZ  = DUMMYP%PREV%XYZ
             LAST%GRAD = DUMMYP%PREV%GRAD
             LAST%E    = DUMMYP%PREV%E
             LAST%PREVXYZ = DUMMYP%PREV%PREVXYZ
             IF(DESMINT) THEN
                LAST%XYZCART = DUMMYP%PREV%XYZCART
                LAST%GCART = DUMMYP%PREV%GCART
                LAST%PREVXCART = DUMMYP%PREV%PREVXCART
             ENDIF
             RETURN
          ENDIF
       ENDIF

       DUMMYP => DUMMYP%NEXT
    ENDDO    
    IF (CLSTRINGT) THEN
       IF (LAST%E<LAST%PREV%E) THEN
          REGENE=.TRUE.
          LAST%XYZ  = DUMMYP%PREV%XYZ
          LAST%GRAD = DUMMYP%PREV%GRAD
          LAST%E    = DUMMYP%PREV%E
          LAST%PREVXYZ = DUMMYP%PREV%PREVXYZ
          IF(DESMINT) THEN
             LAST%XYZCART = DUMMYP%PREV%XYZCART
             LAST%GCART = DUMMYP%PREV%GCART
             LAST%PREVXCART = DUMMYP%PREV%PREVXCART
          ENDIF
          RETURN
       ENDIF
    ENDIF
  END SUBROUTINE CHECKENERGIES

  SUBROUTINE GETIMGPOT(IMGP, RMS)
    ! run potential for the given string image, given XYZCART
    ! set IMGP%GRAD, IMGP%GCART, IMGP%E
    ! returns RMS

    IMPLICIT NONE

    TYPE(IMGNODE), POINTER :: IMGP
    DOUBLE PRECISION, INTENT(OUT) :: RMS
    DOUBLE PRECISION :: DUMINT(NINTC)
    LOGICAL :: NOCOOR, NODERV

    IF (DESMINT) THEN
       CALL POTENTIAL(IMGP%XYZCART, IMGP%E, IMGP%GCART, .TRUE., .FALSE., RMS, .FALSE., .FALSE.)

       NOCOOR = .TRUE.; NODERV = .FALSE.
       CALL TRANSFORM(IMGP%XYZCART,IMGP%GCART,DUMINT,IMGP%GRAD, &
                        & NINTC,3*NATOMS,NNZ,NOCOOR,NODERV,KD,INTEPSILON)
    ELSE
       CALL POTENTIAL(IMGP%XYZ, IMGP%E, IMGP%GRAD, .TRUE., .FALSE., RMS, .FALSE., .FALSE.)
    ENDIF

    RETURN
  END SUBROUTINE GETIMGPOT

!==================================================================================================|
! The inexact newton method. Used by climbing string method to acceleate the convergence to the TS.| 
! W. Ren and E. Vanden-Eijnden J. Chem. Phys. 138, 134105(2013)                                    |
!==================================================================================================|
  SUBROUTINE INEXACT_NEWTON(x0,v0)
    USE KEY ,ONLY :STTSRMSCONV,ST_TSSTEP,LAN_DIST,LANSTEP,LANCONV,LANFACTOR
    IMPLICIT NONE
    INTEGER::STTSIT,LANIT
    DOUBLE PRECISION::RMS, be, e, s_old, s, c_old, c, ze_old, ze, al, be_old, sp, tmp, ga, de, ep
    DOUBLE PRECISION::x0(NC),v0(NC),ETS
    DOUBLE PRECISION,ALLOCATABLE,SAVE,DIMENSION(:)::p,v_old,v,w,g0,x,Hv,g,u
    
    IF(.NOT.ALLOCATED(p)) THEN
       ALLOCATE(p(NC),v(NC),v_old(NC),w(NC),g0(NC),x(NC),Hv(NC),g(NC),u(NC))
    ENDIF
    
    STTSIT=0
    RMS=STTSRMSCONV+1D0
    DO 
       p(:)=0d0                                 
       v_old(:)=0d0                             
       w(:)=0d0                                 
       be=1d0                                   
       e=1d0                                    
       s_old=0d0                                
       s=0d0                                    
       c_old=-1d0                               
       c=-1d0                                   
       ze_old=1d0                               
       ze=1d0
       CALL POTENTIAL(x0,ETS,g0,.TRUE., .FALSE., RMS, .FALSE., .FALSE.) 
       v=-g0/sqrt(sum(g0*g0)) 
       print '(A,1x,2G20.10,1x,I6)', '>>> INEXACT NEWTON RMS,Energy,',RMS,ETS,STTSIT
       IF(RMS<STTSRMSCONV ) THEN
          IF (PTEST) print*, 'Converged to TS'
          EXIT
       ENDIF
       LANIT=0   
       DO
          x=x0+v*LAN_DIST
          CALL POTENTIAL(x,ETS,g,.TRUE., .FALSE., RMS, .FALSE., .FALSE.)   
          Hv=(g-g0)/LAN_DIST*LANFACTOR ! The factor "LANFACTOR" is needed to control the stepsize.
          al=sum(v*Hv)                   
          u=Hv-al*v-be*v_old                    
          be_old=be                                
          be=sqrt(sum(u**2d0))                     
          v_old=v                                  
          v=u/be                                   
          ep=s_old*be_old                          
          tmp=-c_old*be_old                        
          de=c*tmp+s*al                            
          tmp=s*tmp-c*al                           
          ga=sqrt(tmp**2d0+be**2d0)                
          c_old=c                                  
          c=tmp/ga                                 
          s_old=s                                  
          s=be/ga                                  
          tmp=-(ep*ze_old+de*ze)/ga                
          ze_old=ze                                
          ze=tmp                                   
          w=s_old*w-c_old*v_old                    
          e=e*dabs(c_old*s/c)
          IF(e<=LANCONV) THEN                    
             EXIT                                  
          ENDIF                                  
          p=p+ze*(c*w+s*v)                               
          LANIT=LANIT+1
          IF(LANIT==LANSTEP) THEN
             IF (PTEST) print*, 'LANCSOZ FAIL TO CONVERGE'
             EXIT      
          ENDIF
       ENDDO
       STTSIT=STTSIT+1
       IF(STTSIT>ST_TSSTEP) THEN
          IF (PTEST) print*, 'FAIL TO CONVERGE'
          EXIT
       ENDIF
       x0=x0+p
    ENDDO
  END SUBROUTINE INEXACT_NEWTON


END MODULE GROWSTRINGUTILS

