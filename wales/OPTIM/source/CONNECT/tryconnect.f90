!   CONNECT module is an implementation of a connection algorithm for finding rearrangement pathways.
!   Copyright (C) 2003-2006 Semen A. Trygubenko and David J. Wales
!   This file is part of CONNECT module. CONNECT module is part of OPTIM.
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
MODULE TRYCONNECTMODULE
     IMPLICIT NONE
     CONTAINS

     SUBROUTINE TRYCONNECT(JS,JF,TSREDO,REDOPATH,REDOPATHXYZ,USEINT,USEINTLJ)
          USE NEWNEBMODULE
          USE PORFUNCS
          USE CONNECTDATA
          USE KEYCONNECT
          USE CONNECTUTILS
          USE NEBTOCONNECT
          USE KEYNEB, ONLY : NIMAGE, NITERMAX, READGUESS
          USE KEY, ONLY : UNRST, FILTH, FILTHSTR, DUMPALLPATHS, TWOD, MAXTSENERGY, RIGIDBODY, &
  &                   MAXMAXBARRIER, DIJKSTRALOCAL, INTNTRIESMAX, &
  &                   PERMDIST, MAXBARRIER, GROWSTRINGT, BULKT, AMBERT, NABT, AMBER12T, FREEZE, FROZEN, NFREEZE, &
  &                   PUSHOFF, MAXBFGS, MIN1REDO, MIN2REDO, REDOKADD, D1INIT, D2INIT, REDOE1, REDOE2, &
  &                   AMHT, SEQ, INTIMAGE, NEBKFINAL, INTIMAGEINCR, MAXINTIMAGE, SETCHIRAL, REDOTS, VARIABLES, &
  &                   FLATPATHT, QCIIMAGE, PUSHOPTT
          USE MODGUESS
          USE MODUNRES
          USE MODCHARMM, ONLY : CHRMMT,NCHENCALLS,CHECKOMEGAT,CHECKCHIRALT,NOCISTRANS,MINOMEGA
          USE MODAMBER9, ONLY : NOCISTRANSRNA, NOCISTRANSDNA, GOODSTRUCTURE1, GOODSTRUCTURE2, CISARRAY1, CISARRAY2,CHIARRAY1,CHIARRAY2
          USE MODMEC
          USE KEYUTILS
          USE COMMONS, ONLY : NINTS, PARAM1, PARAM2, PARAM3, ZSYM, EVDISTTHRESH, REDOPATHNEB, DEBUG
          USE PORFUNCS
          USE AMHGLOBALS, ONLY : NMRES
          USE GSDATA, ONLY : GSITERDENSITY, GSCURITERD=>ITERD
! hk286
          USE GENRIGID
          USE CHIRALITY, ONLY: CIS_TRANS_CHECK, CHIRALITY_CHECK
          IMPLICIT NONE
          DOUBLE PRECISION RMAT(3,3), DISTFAC

          INTEGER,INTENT(IN) :: JS,JF
          LOGICAL,INTENT(IN) :: USEINT, USEINTLJ
          LOGICAL, ALLOCATABLE :: FOUNDBEFORE(:), DOAGAIN(:)
          
          INTEGER         :: I,I1,UNIQUE=0,MINPLUSPOS,MINMINUSPOS,J1,NC1,NC2,GLY_COUNT,LINTIMAGE,SAMEAS,ISTAT
          DOUBLE PRECISION         :: EDUMMY,EDUMMY2,TMPTS(NOPT),SAVEPUSHOFF,SAVEMAXBFGS
          DOUBLE PRECISION,POINTER :: QPLUS(:),QMINUS(:),EPLUS,EMINUS
          LOGICAL         :: PLUSNEW,MINUSNEW,PATHFAILT,AMIDEFAIL,CHIRALFAIL,RERUN
          CHARACTER       :: ITSTRING*80, EOFSSTRING*80
          DOUBLE PRECISION :: STARTINT(NINTS), FINISHINT(NINTS),DUM, LDUMMY
          LOGICAL REDOPATH, REDOPATHXYZ, PERMUTE, EXISTS, NORERUN, PERMTEST
          DOUBLE PRECISION TSREDO(NOPT), ETSPREV, ETSDUM, DIST2, QP(NOPT), QM(NOPT), LGDUMMY(NOPT)
          DOUBLE PRECISION NEBKFINALSAVE
          INTEGER INVERT, INDEX(NATOMS), J2, IMATCH
          CHARACTER(LEN=2) ZDUM
          CHARACTER(LEN=10) DUMMYS          
          CHARACTER(LEN=80) TRYFNAME          
          DOUBLE PRECISION, POINTER :: PINTERPCOORDS(:), PENERGY
          LOGICAL MINNEW
          INTEGER POSITION, NMINSAVE, NMINORIG
          DOUBLE PRECISION DIST1P, DIST1M, DIST2P, DIST2M
          LOGICAL PATHFAILED
          INTEGER :: NREDOPATHTRIES1=1
          INTEGER :: NREDOPATHTRIES2=1
          DOUBLE PRECISION :: REDOSTRETCH=5.0D0
          TYPE (MINFOUNDTYPE) :: MYMINFOUND(NMINMAX)
          TYPE (TSFOUNDTYPE) :: MYTSFOUND(NTSMAX)
! hk286
          DOUBLE PRECISION :: XRIGIDCOORDS(DEGFREEDOMS), XCOORDS(NOPT)
          DOUBLE PRECISION :: DISTFLAT, PATHLENGTH(3), EOFS(3), FLATSUM2, FLATSUM4, MINIM

!
!  We want to return here to rerun the path for transition states read in with REDOPATH
!  in case the connection fails.
!
          NC1=0 ! counter for changes of PUSHOFF value
          NC2=0 ! counter for changes of BFGSSTEP value
          SAVEPUSHOFF=PUSHOFF
          SAVEMAXBFGS=MAXBFGS
          NMINORIG=NMIN
          RERUN=.FALSE.
10        CONTINUE

          IF (GROWSTRINGT) THEN
             IF (NCONDONE.EQ.1.AND.FCD) THEN
                GSCURITERD = GSITERDENSITY
             ELSE
                GSCURITERD = ITERDENSITY
             ENDIF
          ENDIF

          IF (CHRMMT) NCHENCALLS = 999 ! update non-bonded list on next call to potential.
!
!  Subroutine CHECKPAIR puts the endpoints into optimal alignment.
!
          IF (.NOT.REDOPATH) THEN
             CALL CHECKPAIR(JS,JF,PERMTEST)
          ELSE
            PERMTEST=.FALSE.
            IF (ABS(MI(JS)%DATA%E-MI(JF)%DATA%E) < EDIFFTOL) PERMTEST=.TRUE. ! must initialise this logical
            CALL ALIGN_DECIDE(TSREDO,MIN1REDO,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
            D1INIT=D
            CALL ALIGN_DECIDE(TSREDO,MIN2REDO,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
            D2INIT=D

             IF (DEBUG) PRINT '(A,2F20.10)',' tryconnect> Initial distances of transition state to minima are :',D1INIT,D2INIT
          ENDIF

          IF (GUESSPATHT) THEN
            IF (UNRST) THEN
               DO J1=1,NRES
                  C(1,J1)=MI(JS)%DATA%X(6*(J1-1)+1)
                  C(2,J1)=MI(JS)%DATA%X(6*(J1-1)+2)
                  C(3,J1)=MI(JS)%DATA%X(6*(J1-1)+3)
                  C(1,J1+NRES)=MI(JS)%DATA%X(6*(J1-1)+4)
                  C(2,J1+NRES)=MI(JS)%DATA%X(6*(J1-1)+5)
                  C(3,J1+NRES)=MI(JS)%DATA%X(6*(J1-1)+6)
               ENDDO
               CALL UPDATEDC
               CALL INT_FROM_CART(.TRUE.,.FALSE.)
               CALL GEOM_TO_VAR(NINTS,STARTINT(1:NINTS)) 
               DO J1=1,NRES
                  C(1,J1)=MI(JF)%DATA%X(6*(J1-1)+1)
                  C(2,J1)=MI(JF)%DATA%X(6*(J1-1)+2)
                  C(3,J1)=MI(JF)%DATA%X(6*(J1-1)+3)
                  C(1,J1+NRES)=MI(JF)%DATA%X(6*(J1-1)+4)
                  C(2,J1+NRES)=MI(JF)%DATA%X(6*(J1-1)+5)
                  C(3,J1+NRES)=MI(JF)%DATA%X(6*(J1-1)+6)
               ENDDO
               CALL UPDATEDC
               CALL INT_FROM_CART(.TRUE.,.FALSE.)
               CALL GEOM_TO_VAR(NINTS,FINISHINT(1:NINTS))
               CALL GUESSPATH(STARTINT,FINISHINT,NINTS,EDIFFTOL,NATOMS)
            ELSE
               CALL GUESSPATH(MI(JS)%DATA%X,MI(JF)%DATA%X,NOPT,EDIFFTOL,NATOMS)
            ENDIF
          ! how many images to use?
            IF (NIMAGE > IMAGEMAX) PRINT*,'WARNING - Nimage is greater than ImageMax'
            IF (NIMAGE < 2       ) PRINT*,'WARNING - Nimage is < 2'
            NIMAGE=NINTERP
!           NIterMax = Nimage*IterDensity ! try zero neb iterations if we have a GUESSPATH path
            NITERMAX = 0
            IF (NINTERP.LT.2) THEN ! NO IMAGES FROM GUESSPATH - REVERT TO USUAL SCHEME
               IF (.NOT.(NCONDONE==1 .AND. FCD)) THEN
!                   NIMAGE=IMAGEDENSITY*MI(JF)%DATA%D(JS) &
!                       +IMAGEINCR*IMAGEDENSITY*MI(JF)%DATA%D(JS)*(MI(JF)%DATA%NTRIES(JS)-1)
                  DISTFAC=MI(JF)%DATA%D(JS)
                  IF (DIJKSTRALOCAL.NE.1.0D0) THEN
                     CALL ALIGN_DECIDE(MI(JS)%DATA%X,MI(JF)%DATA%X,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD, &
     &                                   DISTFAC,DIST2,RIGIDBODY,RMAT)
                  ENDIF
                  IF (.NOT.(USEINT.OR.USEINTLJ)) THEN
                     NIMAGE=DISTFAC*(IMAGEDENSITY+IMAGEINCR*(MI(JF)%DATA%NTRIES(JS)-1))
                     IF (NIMAGE >= IMAGEMAX) THEN
                        NIMAGE = IMAGEMAX
                        MI(JF)%DATA%NTRIES(JS)=NTRIESMAX ! no point trying again with the same number of images
                     ENDIF
                     IF (NIMAGE < 2       ) NIMAGE = 2
                     NITERMAX = NIMAGE*ITERDENSITY
                  ENDIF
               ENDIF

            ENDIF
          ELSEIF (MECCANOT) THEN
          ! how many images to use?
            NIMAGE=NINT(MIN(MECIMDENS*MI(JF)%DATA%D(JS),MECMAXIMAGES*1.0D0)) ! IMAGE DENSITY TIMES DISTANCE
!           if (Nimage > ImageMax) PRINT*,'WARNING - Nimage is greater than ImageMax'
            IF (NIMAGE < 1       ) NIMAGE=1
            NITERMAX=NINT(MIN(NIMAGE*MECITDENS,MECMAXIT*1.0D0)) ! NUMBER OF IMAGES TIMES ITERATION DENSITY
            IF (UNRST) THEN
               DO J1=1,NRES
                  C(1,J1)=MI(JS)%DATA%X(6*(J1-1)+1)
                  C(2,J1)=MI(JS)%DATA%X(6*(J1-1)+2)
                  C(3,J1)=MI(JS)%DATA%X(6*(J1-1)+3)
                  C(1,J1+NRES)=MI(JS)%DATA%X(6*(J1-1)+4)
                  C(2,J1+NRES)=MI(JS)%DATA%X(6*(J1-1)+5)
                  C(3,J1+NRES)=MI(JS)%DATA%X(6*(J1-1)+6)
               ENDDO
               CALL UPDATEDC
               CALL INT_FROM_CART(.TRUE.,.FALSE.)
               CALL GEOM_TO_VAR(NINTS,STARTINT(1:NINTS))
               DO J1=1,NRES
                  C(1,J1)=MI(JF)%DATA%X(6*(J1-1)+1)
                  C(2,J1)=MI(JF)%DATA%X(6*(J1-1)+2)
                  C(3,J1)=MI(JF)%DATA%X(6*(J1-1)+3)
                  C(1,J1+NRES)=MI(JF)%DATA%X(6*(J1-1)+4)
                  C(2,J1+NRES)=MI(JF)%DATA%X(6*(J1-1)+5)
                  C(3,J1+NRES)=MI(JF)%DATA%X(6*(J1-1)+6)
               ENDDO
               CALL UPDATEDC
               CALL INT_FROM_CART(.TRUE.,.FALSE.)
               CALL GEOM_TO_VAR(NINTS,FINISHINT(1:NINTS))
               CALL UNMECCANO(.TRUE.,DEBUG,DUM,.FALSE.,STARTINT,FINISHINT,DUM,DUM,DUM,DUM)
            ELSE
               CALL MECCANO(.TRUE.,DEBUG,DUM,.FALSE.,MI(JS)%DATA%X,MI(JF)%DATA%X,DUM,DUM,DUM,DUM)
            ENDIF
            NITERMAX = 0 ! TRY ZERO NEB ITERATIONS IF WE HAVE A MECCANO PATH

          ELSEIF (REDOPATH.AND.(.NOT.REDOPATHNEB)) THEN
            NIMAGE=1
            NITERMAX = 0 ! TRY ZERO NEB ITERATIONS IF WE HAVE POINTS IN TSREDO
          ELSEIF (NCONDONE==1 .AND. FCD) THEN ! FIRST CYCLE DIFFERENT - PARAMETERS SUPPLIED USING NEWNEB
                                          ! keyword or newneb defaults will be used instead
               WRITE(*,'(A)') " tryconnect> First DNEB calculation will use parameters from the NEWNEB line in odata"
          ELSE
            IF (REDOPATHNEB) THEN
               NIMAGE=(D1INIT+D2INIT)*IMAGEDENSITY
            ELSE
               NIMAGE=MI(JF)%DATA%D(JS)*(IMAGEDENSITY+IMAGEINCR*(MI(JF)%DATA%NTRIES(JS)-1))
            ENDIF
            IF (.NOT.(USEINT.OR.USEINTLJ)) THEN
!              PRINT '(A,F10.2,G20.10,F10.2,2I8)',' tryconnect> IMAGEDENSITY,dist,IMAGEINCR,tries,NIMAGE=', &
! &                       IMAGEDENSITY,MI(JF)%DATA%D(JS), &
! &                       IMAGEINCR,MI(JF)%DATA%NTRIES(JS),NIMAGE
                  IF (NIMAGE >= IMAGEMAX) THEN
                     NIMAGE = IMAGEMAX
!
! No point trying again with the same number of images
! Actually, there probably is if we're using constraints because we could obtain
! a different initial alignment for a different number of interpolation images.
!
!                    MI(JF)%DATA%NTRIES(JS)=NTRIESMAX 
                  ENDIF
                  IF (NIMAGE < 2       ) NIMAGE = 2
                  NITERMAX = NIMAGE*ITERDENSITY
             ENDIF
          ENDIF
         
          ! book-keeping :-)
          IF (.NOT.(MECCANOT.OR.(REDOPATH.AND.(.NOT.REDOPATHNEB)))) THEN             
             IF (GROWSTRINGT) THEN
                WRITE(CHR,'(i7)') INT(NIMAGE*GSCURITERD)
                WRITE(*,'(/a)',advance='no') ' tryconnect> '//trim(adjustl(chr))//'-iteration GS run for minima '
             ELSE IF (USEINT) THEN
                WRITE(*,'(/a)',advance='no') ' tryconnect> Interpolation for minima '
             ELSE IF (USEINTLJ) THEN
                WRITE(*,'(/a)',advance='no') ' tryconnect> Interpolation with interpLJ potential for minima '
             ELSE
                WRITE(CHR,'(I7)') NITERMAX
                WRITE(*,'(/a)',advance='no') ' tryconnect> '//trim(adjustl(chr))//'-iteration DNEB run for minima '
             ENDIF
             WRITE(CHR,'(i5)') js
             WRITE(*,'(a)',advance='no') trim(adjustl(chr))
          
             IF (MI(JS)%DATA%S) THEN
                  WRITE(*,'(a)',advance='no') '_S'
             ELSEIF (MI(JS)%DATA%F) THEN
                  WRITE(*,'(a)',advance='no') '_F'
             ELSE
                  WRITE(*,'(a)',advance='no') '_U'
             ENDIF
             WRITE(CHR,'(i5)') jf
             WRITE(*,'(a)',advance='no') ' and '//trim(adjustl(chr))
             IF (MI(JF)%DATA%S) THEN
                  WRITE(*,'(a)',advance='no') '_S'
             ELSEIF (MI(JF)%DATA%F) THEN
                  WRITE(*,'(a)',advance='no') '_F'
             ELSE
                  WRITE(*,'(a)',advance='no') '_U'
             ENDIF
          
             ! getting ts candidates from NEB

             IF (USEINT.OR.USEINTLJ) THEN
                LINTIMAGE=INTIMAGE
!               INTIMAGE=MIN(LINTIMAGE+INTIMAGEINCR*(MI(JF)%DATA%INTNTRIES(JS)-1),MAXINTIMAGE)
!               IF (INTIMAGE.GE.MAXINTIMAGE) MI(JF)%DATA%INTNTRIES(JS)=INTNTRIESMAX ! no point trying again 
                WRITE(CHR,'(i5)') INTIMAGE
             ELSE
                WRITE(CHR,'(i5)') NIMAGE
             ENDIF
             IF (.NOT.(USEINT.OR.USEINTLJ)) THEN
                WRITE(*,'(a)',advance='no') ' using '//TRIM(ADJUSTL(CHR))//' images '
!               IF (MI(JF)%DATA%NTRIES(JS) > 1) THEN
                     WRITE(CHR,'(i5)') mi(jf)%data%ntries(js)
                     WRITE(*,'(a)',advance='no') ' (DNEB attempt #'//TRIM(ADJUSTL(CHR))//') '
!               ENDIF
             ELSE
                WRITE(*,'(a)',advance='no') ' using '//TRIM(ADJUSTL(CHR))//' initial images '
                WRITE(CHR,'(i5)') MI(JF)%DATA%INTNTRIES(JS)
                WRITE(*,'(a)',advance='no') ' (interp attempt #'//TRIM(ADJUSTL(CHR))//') '
             ENDIF
             WRITE(*,'(a)') ' ...'
          ENDIF

          EVDISTTHRESH=-0.2D0 ! DJW ! not being used at the moment

          IF (EVDISTTHRESH.GT.0.0D0) THEN
             NFREEZE=0
             DO J1=1,NATOMS
                IF (VARIABLES) THEN
                   LDUMMY=(MI(JS)%DATA%X(J1)-MI(JF)%DATA%X(J1))**2 
                ELSE
                   LDUMMY=(MI(JS)%DATA%X(3*(J1-1)+1)-MI(JF)%DATA%X(3*(J1-1)+1))**2 &
     &             +(MI(JS)%DATA%X(3*(J1-1)+2)-MI(JF)%DATA%X(3*(J1-1)+2))**2 &
     &             +(MI(JS)%DATA%X(3*(J1-1)+3)-MI(JF)%DATA%X(3*(J1-1)+3))**2
                ENDIF
                IF (LDUMMY.GE.EVDISTTHRESH) THEN
                   PRINT '(A,I8,A,F20.10)',' tryconnect> displacement of atom ',J1,' between end points is ',SQRT(LDUMMY)
                ELSE
                   FROZEN(J1)=.TRUE.
                   NFREEZE=NFREEZE+1
                ENDIF
             ENDDO
             IF (NFREEZE.GT.0) THEN
                FREEZE=.TRUE.
                PRINT '(A,I8,A)',' tryconnect> initially freezing ',NFREEZE,' atoms'
             ENDIF
          ENDIF

!
! No point trying again with the same number of images if we're not using constraints.
!
          IF ((.NOT.(USEINT.OR.USEINTLJ)).AND.(NIMAGE>=IMAGEMAX)) MI(JF)%DATA%NTRIES(JS)=NTRIESMAX 
          IF (REDOPATHXYZ) THEN
             NORERUN=.FALSE.
!
! NConDone should be NTS+1, otherwise one previous run has failed to get a new ts from the path.<n>.xyz file,
! which should never happen?!
!
             IF (NCONDONE.NE.NTS+1) THEN
                PRINT '(2(A,I8))',' tryconnect> ERROR - NCONDONE=',NConDone,' NTS=',NTS
                STOP
             ENDIF
             CALL MKFNAMES(NTS+1,FILTH,FILTHSTR,ITSTRING,EOFSSTRING)
             INQUIRE(FILE=TRIM(ADJUSTL(ITSTRING)),EXIST=EXISTS)
             IF (EXISTS) THEN ! allows for rerun with different energy difference criterion for
                              ! consecutive frames in path without redoing original path
                PRINT '(2A)',' tryconnect> Reading data for ts from existing file ',TRIM(ADJUSTL(itstring))
                NORERUN=.TRUE.
             ELSE
                PRINT '(3A)',' tryconnect> No file ',TRIM(ADJUSTL(itstring)),' found'
                REDOPATH=.FALSE.
                REDOPATHXYZ=.FALSE.
                STOP ! REDOPATHXYZ is only going to work after a run from redopoints, in which case
                     ! we should finish after reading the last path.
             ENDIF
             IF (NORERUN) THEN
                NTSFOUND=1
                PRINT *,'tryconnect Z NOPT=',NOPT
                ALLOCATE(TSFOUND(NTSFOUND)%E,TSFOUND(NTSFOUND)%COORD(NOPT),&
   &                     TSFOUND(NTSFOUND)%EVALMIN,TSFOUND(NTSFOUND)%VECS(NOPT))
                OPEN(UNIT=89,FILE=ITSTRING,STATUS='OLD')
                ETSPREV=1.0D100
                TSFOUND(1)%E=-1.0D100
                DO
                   READ(89,*)
                   READ(89,*) DUMMYS, ETSDUM
                   IF ((TSFOUND(1)%E.GT.ETSDUM).AND.(TSFOUND(1)%E.GT.ETSPREV)) GOTO 111
                   DO J1=1,NATOMS
                      READ(89,*) ZDUM,TSFOUND(1)%COORD(3*(J1-1)+1),TSFOUND(1)%COORD(3*(J1-1)+2),TSFOUND(1)%COORD(3*(J1-1)+3)
                   ENDDO
                   ETSPREV=TSFOUND(1)%E
                   TSFOUND(1)%E=ETSDUM
                ENDDO
111             CONTINUE
                IF (DEBUG) PRINT *,'TSfound(1)%E=',TSfound(1)%E
                PRINT '(A,G20.10)',' tryconnect> Ets=',TSfound(1)%E
                CLOSE(89)
             ELSE
                IF (PERMTEST.AND.(MOD(NIMAGE,2).NE.0).AND.(.NOT.REDOPATH).AND.(.NOT.READGUESS)) THEN
                   PRINT '(A)',' tryconnect> Changing to an even number of images for possible permutational isomer'
                   NIMAGE=NIMAGE+1
                ENDIF
                IF (UNRST) THEN 
!                  CALL NEWNEB(REDOPATH,TSREDO,-1.0D100,MI(JS)%DATA%X,-1.0D100,MI(JF)%DATA%X,.TRUE.,NINTSIN=NINTS)
                   PRINT *,'tryconnect> unres currently not allowed DJW'
                   STOP
                ELSEIF (USEINT) THEN
!                  NEBKFINALSAVE=NEBKFINAL
!                  IF (MOD(MI(JF)%DATA%NTRIES(JS),3).EQ.1) THEN
!                     NEBKFINAL=NEBKFINAL*10.0D0
!                     PRINT '(A,F20.10)',' tryconnect> Increasing NEB force constant to ',NEBKFINAL
!                  ELSEIF (MOD(MI(JF)%DATA%NTRIES(JS),3).EQ.2) THEN
!                     NEBKFINAL=NEBKFINAL/10.0D0
!                     PRINT '(A,F20.10)',' tryconnect> Decreasing NEB force constant to ',NEBKFINAL
!                  ENDIF
                   CALL INTLBFGS(MI(JS)%DATA%X,MI(JF)%DATA%X)
                   NIMAGE=QCIIMAGE
                   IF (.NOT.(NCONDONE==1 .AND. FCD)) NITERMAX=NIMAGE*ITERDENSITY
               WRITE(*,'(A)') " tryconnect> First DNEB calculation will use parameters from the NEWNEB line in odata"

                   CALL NEWNEB(REDOPATH,TSREDO,-1.0D100,MI(JS)%DATA%X,-1.0D100,MI(JF)%DATA%X,.TRUE.)
!                  NEBKFINAL=NEBKFINALSAVE
!                  INTIMAGE=LINTIMAGE
                ELSEIF (USEINTLJ) THEN
                   CALL INTLBFGSLJ(MI(JS)%DATA%X,MI(JF)%DATA%X,NMINFOUND,NTSFOUND,MINFOUND,TSFOUND)
!                  INTIMAGE=LINTIMAGE
                ELSE
                   CALL NEWNEB(REDOPATH,TSREDO,-1.0D100,MI(JS)%DATA%X,-1.0D100,MI(JF)%DATA%X,.TRUE.)
                ENDIF
             ENDIF
          ELSE
             IF (PERMTEST.AND.(MOD(NIMAGE,2).NE.0).AND.(.NOT. REDOPATH).AND.(.NOT.READGUESS)) THEN
                PRINT '(A)',' tryconnect> Changing to an even number of images for possible permutational isomer'
                NIMAGE=NIMAGE+1
             ENDIF
             IF (UNRST) THEN 
!                CALL NEWNEB(REDOPATH,TSREDO,MI(JS)%DATA%E,MI(JS)%DATA%X,MI(JF)%DATA%E,MI(JF)%DATA%X,.TRUE.,NINTSIN=NINTS)
                 PRINT *,'tryconnect> unres currently not allowed DJW'
                 STOP
             ELSEIF (USEINT) THEN
!               NEBKFINALSAVE=NEBKFINAL
!               IF (MOD(MI(JF)%DATA%NTRIES(JS),3).EQ.1) THEN
!                  NEBKFINAL=NEBKFINAL*10.0D0
!                  PRINT '(A,F20.10)',' tryconnect> Increasing NEB force constant to ',NEBKFINAL
!               ELSEIF (MOD(MI(JF)%DATA%NTRIES(JS),3).EQ.2) THEN
!                  NEBKFINAL=NEBKFINAL/10.0D0
!                  PRINT '(A,F20.10)',' tryconnect> Decreasing NEB force constant to ',NEBKFINAL
!               ENDIF
                CALL INTLBFGS(MI(JS)%DATA%X,MI(JF)%DATA%X)
!               NEBKFINAL=NEBKFINALSAVE
                NIMAGE=QCIIMAGE
                IF (.NOT.(NCONDONE==1 .AND. FCD)) NITERMAX=NIMAGE*ITERDENSITY
                CALL NEWNEB(REDOPATH,TSREDO,-1.0D100,MI(JS)%DATA%X,-1.0D100,MI(JF)%DATA%X,.TRUE.)
             ELSEIF (USEINTLJ) THEN
                CALL INTLBFGSLJ(MI(JS)%DATA%X,MI(JF)%DATA%X,NMINFOUND,NTSFOUND,MINFOUND,TSFOUND)
!               INTIMAGE=LINTIMAGE
             ELSE
                IF (REDOPATHNEB) THEN
                   IF (DEBUG) THEN
                      PRINT '(A,2F20.10)',' tryconnect> calling newneb with min-sad-min distances:',D1INIT,D2INIT
                   ENDIF
                   CALL NEWNEB(REDOPATH,TSREDO,REDOE1,MIN1REDO,REDOE2,MIN2REDO,.TRUE.)
                ELSE
                   CALL NEWNEB(REDOPATH,TSREDO,MI(JS)%DATA%E,MI(JS)%DATA%X,MI(JF)%DATA%E,MI(JF)%DATA%X,.TRUE.)
                ENDIF
             ENDIF
          ENDIF
!
!  We may have new minima rather than new ts. 
!  Deal with this first.
!
          IF (NMINFOUND.GT.0) THEN
!
! Check if it;s a new min - follow what we do in newconnect.f90
! Save data for new min using horrible write/read trick to save in static pointer variables
! Create distances - fudging NEB neighbour values 
! 
             NMINSAVE=NMIN
             IF(DEBUG) WRITE(*,*) "tryconnect> ", NMINFOUND, "new minima found" ! sn402
             DO I=1,NMINFOUND
                NULLIFY(PINTERPCOORDS,PENERGY)
                ALLOCATE(PINTERPCOORDS(NOPT),PENERGY)
                OPEN(UNIT=781,FILE='minscratch',STATUS='UNKNOWN')
                WRITE(781,*) MINFOUND(I)%E,MINFOUND(I)%COORD(1:NOPT)
                REWIND(781)
                READ(781,*) PENERGY,PINTERPCOORDS
                CLOSE(781)
                CALL ISNEWMIN(PENERGY,PINTERPCOORDS,POSITION,MINNEW,REDOPATH,PERMUTE,INVERT,INDEX,IMATCH)
                IF (MINNEW) THEN
                   CALL ADDNEWMIN(PENERGY,PINTERPCOORDS)
!                  PRINT*, PINTERPCOORDS(:)
                   WRITE(*,'(A,I7,A,G20.10)') ' tryconnect> added new minimum ',NMIN,' energy=',PENERGY
                   IF (USEINT.OR.USEINTLJ) THEN ! write minimum to file min<n> to enable debugging
                      WRITE(TRYFNAME,'(I6)') NMIN
                      TRYFNAME='min' // TRIM(ADJUSTL(TRYFNAME))
                      OPEN(UNIT=89,FILE=TRYFNAME,STATUS='UNKNOWN')
                      IF (AMHT) THEN
                         GLY_COUNT=0
                         DO J2=1,NMRES
                            IF (SEQ(J2).EQ.8) THEN
                               WRITE(89,'(3G25.15)') PINTERPCOORDS(9*(J2-1)+1-GLY_COUNT*3), &
     &                                               PINTERPCOORDS(9*(J2-1)+2-GLY_COUNT*3), &
     &                                               PINTERPCOORDS(9*(J2-1)+3-GLY_COUNT*3)
                               WRITE(89,'(3G25.15)') PINTERPCOORDS(9*(J2-1)+1-GLY_COUNT*3), &
     &                                               PINTERPCOORDS(9*(J2-1)+2-GLY_COUNT*3), &
     &                                               PINTERPCOORDS(9*(J2-1)+3-GLY_COUNT*3)
                               WRITE(89,'(3G25.15)') PINTERPCOORDS(9*(J2-1)+4-GLY_COUNT*3), &
     &                                               PINTERPCOORDS(9*(J2-1)+5-GLY_COUNT*3), &
     &                                               PINTERPCOORDS(9*(J2-1)+6-GLY_COUNT*3)
                               GLY_COUNT=GLY_COUNT +1
                            ELSE
                               WRITE(89,'(3G25.15)') PINTERPCOORDS(9*(J2-1)+1-GLY_COUNT*3), &
     &                                               PINTERPCOORDS(9*(J2-1)+2-GLY_COUNT*3), &
     &                                               PINTERPCOORDS(9*(J2-1)+3-GLY_COUNT*3)
                               WRITE(89,'(3G25.15)') PINTERPCOORDS(9*(J2-1)+4-GLY_COUNT*3), &
     &                                               PINTERPCOORDS(9*(J2-1)+5-GLY_COUNT*3), &
     &                                               PINTERPCOORDS(9*(J2-1)+6-GLY_COUNT*3)
                               WRITE(89,'(3G25.15)') PINTERPCOORDS(9*(J2-1)+7-GLY_COUNT*3), &
     &                                               PINTERPCOORDS(9*(J2-1)+8-GLY_COUNT*3), &
     &                                               PINTERPCOORDS(9*(J2-1)+9-GLY_COUNT*3)
                            ENDIF
                         ENDDO
                      ELSE
                         WRITE(89,'(3G25.15)') PINTERPCOORDS(1:NOPT)
                      ENDIF
                      CLOSE(89)
                   ENDIF
!                  DO J2=1,NMIN-1
!                     PRINT '(A,2I8,G20.10)','min1 min2 d=',NMIN,J2,MI(NMIN)%DATA%D(J2)
!                  ENDDO
                ELSE
                   WRITE(*,'(A,I7)') ' tryconnect> found old minimum ',POSITION
                ENDIF
                NULLIFY(PINTERPCOORDS,PENERGY)
                DEALLOCATE(MINFOUND(I)%E,MINFOUND(I)%COORD)
             ENDDO
             if(allocated(FOUNDBEFORE)) DEALLOCATE(FOUNDBEFORE,DOAGAIN)
             RETURN ! assumes that we have no TS if we have new minima. Probably OK.
          ENDIF

! saving new ts into ts rack; otherwise - free memory immediately
   
          NTSOLD=NTS
          UNIQUE=0
          IF (ALLOCATED(FOUNDBEFORE)) DEALLOCATE(FOUNDBEFORE)
          IF (ALLOCATED(DOAGAIN)) DEALLOCATE(DOAGAIN)
          ALLOCATE(FOUNDBEFORE(NTSFOUND),DOAGAIN(NTSFOUND))
          FOUNDBEFORE(1:NTSFOUND)=.FALSE.
          DOAGAIN(1:NTSFOUND)=.FALSE.
          DO I=1,NTSFOUND
!              PRINT '(A,2G20.10)',' tryconnect> TSFOUND(i)%E, MAXTSENERGY=',TSfound(i)%E, MAXTSENERGY
! hk286
               AMIDEFAIL=.FALSE.
               CHIRALFAIL=.FALSE.
               IF (CHRMMT) THEN
                  IF (RIGIDINIT) THEN
                     XRIGIDCOORDS(1:DEGFREEDOMS) = TSFOUND(I)%COORD(1:DEGFREEDOMS)
                     CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, XCOORDS, XRIGIDCOORDS)
                     IF (CHECKOMEGAT) &
                          CALL CHECKOMEGA(XCOORDS,AMIDEFAIL)
                     IF (CHECKCHIRALT) &
                          CALL CHECKCHIRAL(XCOORDS,CHIRALFAIL)
                  ELSE
                     IF (CHECKOMEGAT) &
                          CALL CHECKOMEGA(TSFOUND(I)%COORD,AMIDEFAIL)
                     IF (CHECKCHIRALT) &
                          CALL CHECKCHIRAL(TSFOUND(I)%COORD,CHIRALFAIL)
                  ENDIF
               END IF
               IF (CHRMMT .AND. AMIDEFAIL) THEN
                  PRINT '(A,G20.10,A)',' tryconnect> Transition state with energy ',TSfound(i)%E, &
  &                                   ' ignored, amide cis-trans isomerisation detected'
                  DEALLOCATE(TSFOUND(I)%E,TSFOUND(I)%COORD,TSFOUND(I)%EVALMIN,TSFOUND(I)%VECS)
               ELSEIF (CHRMMT .AND. CHIRALFAIL) THEN
                  PRINT '(A,G20.10,A)',' tryconnect> Transition state with energy ',TSfound(i)%E, &
  &                                    ' ignored, inversion of a chiral CA center detected.'
                  DEALLOCATE(TSFOUND(I)%E,TSFOUND(I)%COORD,TSFOUND(I)%EVALMIN,TSFOUND(I)%VECS)
               ELSEIF (TSFOUND(I)%E.GT.MAXTSENERGY) THEN
                  PRINT '(A,G20.10,A)',' tryconnect> Transition state with energy ',TSfound(i)%E,' ignored'
                  DEALLOCATE(TSFOUND(I)%E,TSFOUND(I)%COORD,TSFOUND(I)%EVALMIN,TSFOUND(I)%VECS)
!
! Allow redopath to add the same transition state more than once.
! For use with different PUSHOFF and BFGSSTEP values in case the connection fails
! to give the minima pair that we actually want.
!
!              ELSEIF ( NCISNEWTS(TSFOUND(I),SAMEAS).OR.REDOPATH ) THEN
               ELSE
                  IF ((.NOT.NCISNEWTS(TSFOUND(I),SAMEAS)).AND.(.NOT.REDOPATH)) FOUNDBEFORE(I)=.TRUE.
                  IF (FOUNDBEFORE(I)) THEN
                     TS(SAMEAS)%DATA%HITS=TS(SAMEAS)%DATA%HITS+1
                     IF ((.NOT.PUSHOPTT).AND.(TS(SAMEAS)%DATA%HITS.LT.REDOTS)) THEN
                        PRINT '(A,I6,A,G20.10)',' tryconnect> Try the path again for ts ',I,' with pushoff=',PUSHOFF/1.0D1
                        IF (NTS==TSRACKSIZE) CALL REALLOCATETSRACK
                        NTS=NTS+1; UNIQUE=UNIQUE+1
                        DOAGAIN(UNIQUE)=.TRUE.
                        TS(NTS)%DATA%E => TSFOUND(I)%E
                        TS(NTS)%DATA%X => TSFOUND(I)%COORD
                        TS(NTS)%DATA%EVALMIN => TSFOUND(I)%EVALMIN
                        TS(NTS)%DATA%VECS => TSFOUND(I)%VECS
                        TS(NTS)%DATA%BAD=.FALSE.
                        TS(NTS)%DATA%HITS=0
                        NULLIFY(TSFOUND(I)%E,TSFOUND(I)%COORD,TSFOUND(I)%EVALMIN,TSFOUND(I)%VECS)
                     ELSE
                        PRINT '(2(A,I6))',' tryconnect> Will not repeat search for TS ',I, &
  &                                       ' same as TS ',SAMEAS
                        DEALLOCATE(TSFOUND(I)%E,TSFOUND(I)%COORD,TSFOUND(I)%EVALMIN,TSFOUND(I)%VECS)
                     ENDIF
                  ELSE
                     IF (NTS==TSRACKSIZE) CALL REALLOCATETSRACK
                     NTS=NTS+1; UNIQUE=UNIQUE+1; NTSGOOD=NTSGOOD+1
                     TS(NTS)%DATA%E => TSFOUND(I)%E
                     TS(NTS)%DATA%X => TSFOUND(I)%COORD
                     TS(NTS)%DATA%EVALMIN => TSFOUND(I)%EVALMIN
                     TS(NTS)%DATA%VECS => TSFOUND(I)%VECS
                     TS(NTS)%DATA%BAD=.FALSE.
                     TS(NTS)%DATA%HITS=0
                     NULLIFY(TSFOUND(I)%E,TSFOUND(I)%COORD,TSFOUND(I)%EVALMIN,TSFOUND(I)%VECS)
                  ENDIF
!              ELSE
!                 IF (NCONDONE==1) PRINT *, 'Discarded TS #',i
!                 DEALLOCATE(TSFOUND(I)%E,TSFOUND(I)%COORD,TSFOUND(I)%EVALMIN,TSFOUND(I)%VECS)
               ENDIF
          ENDDO

!         CALL DUMPTS ! ts.xyz file is not used for anything?

          ! print info as to how many TS are actually useful
!          IF (UNIQUE==NTSFOUND.AND..NOT.UNIQUE==0) THEN
!               IF (NTSFOUND==1) THEN
!                  IF (RERUN) THEN
!                    WRITE(*,'(A)') ' tryconnect> rerunning path for this TS'
!                  ELSE
!                    WRITE(*,'(A)') ' tryconnect> TS appears to be new'
!                  ENDIF
!               ELSE
!                    WRITE(*,'(A)') ' tryconnect> All of TS found appear to be new'
!               ENDIF
!          ELSEIF (UNIQUE < NTSFOUND) THEN
!               WRITE(CHR,'(i7)') unique 
!               WRITE(*,'(1X,A)') trim(adjustl(chr))//' of TS found appear to be new.'
!          ELSEIF (UNIQUE ==0 .AND..NOT.NTSFOUND==0) THEN
!               WRITE(*,'(1X,A)') ' tryconnect> All of TS found are already known'
!          ENDIF

          ! path run for all unique ts
!
!sy349: if the path got from dneb is regarded as flat, the highest-energy image
!will be directly regarded as the transition state, and the two points will be
!regarded connected. The judgement of whether the path is flat is in NEB/output.f90
!
          IF (FLATPATHT) THEN
             WRITE(*,'(A45,I3,A45,I3)') 'tryconnect > flat path found between minimum ', JS, ' and ', JF
             CALL NEWCONNECTION(JS,JF,NTS)
             CALL SETDISTANCE(JS,JF,0.0D0)
             IF (INTERPCOSTFUNCTION) CALL SETINTERP(JS,JF,0.0D0)
             CALL MKFNAMES(NTS,FILTH,FILTHSTR,ITSTRING,EOFSSTRING)

             EOFS(1)=MI(JS)%DATA%E
             EOFS(2)=TS(NTS)%DATA%E
             EOFS(3)=MI(JF)%DATA%E

             PATHLENGTH(1)=0.0D0
             DISTFLAT=0.0D0
             IF (BULKT) THEN
                DO I1=1,NATOMS
                   DISTFLAT=DISTFLAT &
   &               +(MI(JS)%DATA%X(3*(I1-1)+1)-TS(NTS)%DATA%X(3*(I1-1)+1) &
   &            -PARAM1*NINT((MI(JS)%DATA%X(3*(I1-1)+1)-TS(NTS)%DATA%X(3*(I1-1)+1))/PARAM1))**2 &
   &               +(MI(JS)%DATA%X(3*(I1-1)+2)-TS(NTS)%DATA%X(3*(I1-1)+2) &
   &            -PARAM2*NINT((MI(JS)%DATA%X(3*(I1-1)+2)-TS(NTS)%DATA%X(3*(I1-1)+2))/PARAM2))**2 &
   &               +(MI(JS)%DATA%X(3*(I1-1)+3)-TS(NTS)%DATA%X(3*(I1-1)+3) &
   &            -PARAM3*NINT((MI(JS)%DATA%X(3*(I1-1)+3)-TS(NTS)%DATA%X(3*(I1-1)+3))/PARAM3))**2
                ENDDO
                PATHLENGTH(2)=SQRT(DISTFLAT)

                DISTFLAT=0.0D0
                DO I1=1,NATOMS
                   DISTFLAT=DISTFLAT &
   &               +(MI(JF)%DATA%X(3*(I1-1)+1)-TS(NTS)%DATA%X(3*(I1-1)+1) &
   &            -PARAM1*NINT((MI(JF)%DATA%X(3*(I1-1)+1)-TS(NTS)%DATA%X(3*(I1-1)+1))/PARAM1))**2 &
   &               +(MI(JF)%DATA%X(3*(I1-1)+2)-TS(NTS)%DATA%X(3*(I1-1)+2) &
   &            -PARAM2*NINT((MI(JF)%DATA%X(3*(I1-1)+2)-TS(NTS)%DATA%X(3*(I1-1)+2))/PARAM2))**2 &
   &               +(MI(JF)%DATA%X(3*(I1-1)+3)-TS(NTS)%DATA%X(3*(I1-1)+3) &
   &            -PARAM3*NINT((MI(JF)%DATA%X(3*(I1-1)+3)-TS(NTS)%DATA%X(3*(I1-1)+3))/PARAM3))**2
                ENDDO
                PATHLENGTH(3)=SQRT(DISTFLAT)+PATHLENGTH(2)

                DISTFLAT=0.0D0
                DO I1=1,NATOMS
                   DISTFLAT=DISTFLAT &
   &               +(MI(JS)%DATA%X(3*(I1-1)+1)-MI(JF)%DATA%X(3*(I1-1)+1) &
   &            -PARAM1*NINT((MI(JS)%DATA%X(3*(I1-1)+1)-MI(JF)%DATA%X(3*(I1-1)+1))/PARAM1))**2 &
   &               +(MI(JS)%DATA%X(3*(I1-1)+2)-MI(JF)%DATA%X(3*(I1-1)+2) &
   &            -PARAM2*NINT((MI(JS)%DATA%X(3*(I1-1)+2)-MI(JF)%DATA%X(3*(I1-1)+2))/PARAM2))**2 &
   &               +(MI(JS)%DATA%X(3*(I1-1)+3)-MI(JF)%DATA%X(3*(I1-1)+3) &
   &            -PARAM3*NINT((MI(JS)%DATA%X(3*(I1-1)+3)-MI(JF)%DATA%X(3*(I1-1)+3))/PARAM3))**2
                ENDDO
                FLATSUM2=DISTFLAT

                DISTFLAT=0.0D0
                DO I1=1,NATOMS
                   DISTFLAT=DISTFLAT &
   &               +((MI(JS)%DATA%X(3*(I1-1)+1)-MI(JF)%DATA%X(3*(I1-1)+1) &
   &            -PARAM1*NINT((MI(JS)%DATA%X(3*(I1-1)+1)-MI(JF)%DATA%X(3*(I1-1)+1))/PARAM1))**2 &
   &               +(MI(JS)%DATA%X(3*(I1-1)+2)-MI(JF)%DATA%X(3*(I1-1)+2) &
   &            -PARAM2*NINT((MI(JS)%DATA%X(3*(I1-1)+2)-MI(JF)%DATA%X(3*(I1-1)+2))/PARAM2))**2 &
   &               +(MI(JS)%DATA%X(3*(I1-1)+3)-MI(JF)%DATA%X(3*(I1-1)+3) &
   &            -PARAM3*NINT((MI(JS)%DATA%X(3*(I1-1)+3)-MI(JF)%DATA%X(3*(I1-1)+3))/PARAM3))**2)**2
                ENDDO
                FLATSUM4=DISTFLAT
             ELSE
                DO I1=1,NOPT
                   DISTFLAT=DISTFLAT+(MI(JS)%DATA%X(I1)-TS(NTS)%DATA%X(I1))**2
                ENDDO
                PATHLENGTH(2)=SQRT(DISTFLAT)

                DISTFLAT=0.0D0
                DO I1=1,NOPT
                   DISTFLAT=DISTFLAT+(MI(JF)%DATA%X(I1)-TS(NTS)%DATA%X(I1))**2
                ENDDO
                PATHLENGTH(3)=SQRT(DISTFLAT)+PATHLENGTH(2)
                
                DISTFLAT=0.0D0
                DO I1=1,NOPT
                   DISTFLAT=DISTFLAT+(MI(JS)%DATA%X(I1)-MI(JF)%DATA%X(I1))**2
                ENDDO
                FLATSUM2=DISTFLAT

                DISTFLAT=0.0D0
                DO I1=1,NATOMS
                   DISTFLAT=DISTFLAT &
   &                  +((MI(JS)%DATA%X(3*(I1-1)+1)-MI(JF)%DATA%X(3*(I1-1)+1))**2 &
   &                  +(MI(JS)%DATA%X(3*(I1-1)+1)-MI(JF)%DATA%X(3*(I1-1)+1))**2 &
   &                  +(MI(JS)%DATA%X(3*(I1-1)+1)-MI(JF)%DATA%X(3*(I1-1)+1))**2)**2
                ENDDO
                FLATSUM4=DISTFLAT
             ENDIF
             
             OPEN(UNIT=250,FILE=EOFSSTRING,STATUS='UNKNOWN')
             DO I1=1,3
                WRITE(250,'(2G20.10,I6)') PATHLENGTH(I1), EOFS(I1), I1
             ENDDO
             CLOSE(250)

             OPEN (UNIT=251,FILE=ITSTRING,STATUS='UNKNOWN')

             WRITE(251,'(I6)') NATOMS
             WRITE(251,'(A,G25.15)') 'Energy=', EOFS(1)
             DO I1=1,NATOMS
                WRITE(251,'(A2,4X,3G20.10)') &
   &               ZSYM(I1),MI(JS)%DATA%X(3*(I1-1)+1),MI(JS)%DATA%X(3*(I1-1)+2), &
   &               MI(JS)%DATA%X(3*(I1-1)+3)
             ENDDO

             WRITE(251,'(I6)') NATOMS
             WRITE(251,'(A,G25.15)') 'Energy=', EOFS(2)
             DO I1=1,NATOMS
                WRITE(251,'(A2,4X,3G20.10)') &
   &                 ZSYM(I1),TS(NTS)%DATA%X(3*(I1-1)+1),TS(NTS)%DATA%X(3*(I1-1)+2), &
   &                 TS(NTS)%DATA%X(3*(I1-1)+3)
             ENDDO

             WRITE(251,'(I6)') NATOMS
             WRITE(251,'(A,G25.15)') 'Energy=', EOFS(3)
             DO I1=1,NATOMS
                WRITE(251,'(A2,4X,3G20.10)') &
   &                 ZSYM(I1),MI(JF)%DATA%X(3*(I1-1)+1),MI(JF)%DATA%X(3*(I1-1)+2), &
   &                 MI(JF)%DATA%X(3*(I1-1)+3)
             ENDDO
             CLOSE(251)
             
             TS(NTS)%DATA%SLENGTH=PATHLENGTH(3)
             TS(NTS)%DATA%DISP=SQRT(FLATSUM2)
             TS(NTS)%DATA%GAMMA=FLATSUM4*NATOMS/FLATSUM2**2
             TS(NTS)%DATA%NTILDE=FLATSUM2**2/FLATSUM4

             IF (DUMPALLPATHS) CALL MAKEALLPATHINFO(TS(NTS)%DATA%X,MI(JS)%DATA%X,MI(JF)%DATA%X,EOFS(2),EOFS(1),EOFS(3),FRQSTS,FRQSPLUS,FRQSMINUS)

         ELSE 
          DO I=NTS-UNIQUE+1,NTS
               WRITE(CHR,'(i5)') I
               WRITE(*,'(/1x,a)') '>>>>>  Path run for ts '//trim(adjustl(chr))//' ...'
               IF (DOAGAIN(I-NTS+UNIQUE).AND.(.NOT.REDOPATH)) THEN
                  PUSHOFF=PUSHOFF/10.0D0
                  WRITE(*,'(A,G20.10)') ' tryconnect> Trying this transition state path again with pushoff=',PUSHOFF
               ENDIF
               ALLOCATE( QPLUS(NOPT),QMINUS(NOPT),EPLUS,EMINUS )
               CALL MKFNAMES(I,FILTH,FILTHSTR,ITSTRING,EOFSSTRING)
               EDUMMY=TS(I)%DATA%E
               TMPTS(1:NOPT)=TS(I)%DATA%X(1:NOPT)
               ! structure in ts(i)%data%X is a stationary point which is why we don't need to store G and rms for it
               LGDUMMY(1:NOPT)=0.0D0; RMS=0.0D0 ! we must initialize them here, however 
               NORERUN=.FALSE.
               IF (REDOPATH) THEN
!                 CALL MKFNAMES(NCONDONE,FILTH,FILTHSTR,ITSTRING,EOFSSTRING)
                  IF (REDOPATHXYZ) THEN
                     INQUIRE(FILE=TRIM(ADJUSTL(ITSTRING)),EXIST=EXISTS)
                     IF (EXISTS) THEN ! ALLOWS FOR RERUN WITH DIFFERENT ENERGY DIFFERENCE CRITERION FOR
                                      ! consecutive frames in path without redoing original path
                        PRINT '(2A)',' tryconnect> Reading data for minima from existing file ',TRIM(ADJUSTL(itstring))
                        NORERUN=.TRUE.
                     ENDIF
                  ENDIF
               ENDIF
               IF (.NOT.NORERUN) THEN
!                REDOKADD=.TRUE.
! hk286
                 IF (RIGIDINIT) THEN
                 ! Need this to be correctly set in PATH because we call POTENTIAL and do minimisations.
                    ATOMRIGIDCOORDT = .FALSE.
                 ENDIF

!                PRINT '(A)','tryconnect> calling PATH'
                 CALL PATH(TMPTS,EDUMMY,LGDUMMY,RMS,TS(I)%DATA%EVALMIN,TS(I)%DATA%VECS,  &
                    & .FALSE.,QPLUS,QMINUS,DEBUG,EDUMMY2,EPLUS,EMINUS, &
                    & TS(I)%DATA%SLENGTH,TS(I)%DATA%DISP,TS(I)%DATA%GAMMA,TS(I)%DATA%NTILDE,FRQSTS,FRQSPLUS, &
                    & FRQSMINUS,ITSTRING,EOFSSTRING,PATHFAILT)
! hk286
                 IF (RIGIDINIT) THEN
                    ATOMRIGIDCOORDT = .TRUE.
                    XRIGIDCOORDS(1:DEGFREEDOMS) = QPLUS(1:DEGFREEDOMS)
                    CALL TRANSFORMRIGIDTOC (1, NRIGIDBODY, QPLUS, XRIGIDCOORDS)
                    XRIGIDCOORDS(1:DEGFREEDOMS) = QMINUS(1:DEGFREEDOMS)
                    CALL TRANSFORMRIGIDTOC (1, NRIGIDBODY, QMINUS, XRIGIDCOORDS)
                 ENDIF
! hk286
                 QP(1:NOPT)=QPLUS(1:NOPT)
                 QM(1:NOPT)=QMINUS(1:NOPT)
!                REDOKADD=.FALSE.
               ELSE
                  OPEN(UNIT=89,FILE=ITSTRING,STATUS='OLD')
                  READ(89,*)
                  READ(89,*) DUMMYS, EPLUS
                  DO J1=1,NATOMS
                     READ(89,*) ZDUM,QPLUS(3*(J1-1)+1),QPLUS(3*(J1-1)+2),QPLUS(3*(J1-1)+3)
                  ENDDO
                  DO
                     READ(89,*,END=99)
                     READ(89,*) DUMMYS, EMINUS
                     DO J1=1,NATOMS
                        READ(89,*) ZDUM,QMINUS(3*(J1-1)+1),QMINUS(3*(J1-1)+2),QMINUS(3*(J1-1)+3)
                     ENDDO
                  ENDDO
99                CONTINUE
                  IF (DEBUG) PRINT '(A,G20.10)','Eplus=',EPLUS
                  IF (DEBUG) PRINT '(A,G20.10)','Eminus=',Eminus
                  CLOSE(89)
               ENDIF
               DEALLOCATE(TS(I)%DATA%VECS)
! AMBER 12
               IF (AMBER12T.AND.NOCISTRANS) THEN
                   GOODSTRUCTURE1=.TRUE.
                   GOODSTRUCTURE2=.TRUE.
! Check for chiral inversions in DNA/RNA, since cis/trans refers to the substituents on the ribose rings.
                 IF(NOCISTRANSRNA.OR.NOCISTRANSDNA) THEN
                   CALL CHIRALITY_CHECK(QPLUS, GOODSTRUCTURE1)
                   CALL CHIRALITY_CHECK(QMINUS, GOODSTRUCTURE2)
                   IF((.NOT.GOODSTRUCTURE1).OR.(.NOT.GOODSTRUCTURE2)) THEN
                    PRINT '(A,G20.10,A)',' tryconnect> Transition state with energy ',TS(I)%DATA%E, &
  &                                     ' ignored, cis-trans isomerisation in (deoxy)ribose ring or chirality inversion detected.'
                    DEALLOCATE(TS(I)%DATA%EVALMIN)
                    DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                    TS(I)%DATA%BAD=.TRUE.
                    NTSGOOD=NTSGOOD-1
                    CYCLE
                   END IF
                 ELSE
! Otherwise, check for cis/trans flips of peptide bonds
                   CALL CIS_TRANS_CHECK(QMINUS, GOODSTRUCTURE1)
                   CALL CIS_TRANS_CHECK(QPLUS, GOODSTRUCTURE2)

                   IF((.NOT.GOODSTRUCTURE1).OR.(.NOT.GOODSTRUCTURE2)) THEN
                      WRITE(*,'(A)') ' tryconnect> Cis-trans isomerisation of a peptide bond detected '//&
  &                                                '(wrt. the original structure), rejecting'
                    PRINT '(A,G20.10,A)',' tryconnect> Transition state with energy ',TS(I)%DATA%E, &
  &                                     ' ignored, cis-trans isomerisation detected in one or more peptide bonds.'
                    DEALLOCATE(TS(I)%DATA%EVALMIN)
                    DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                    TS(I)%DATA%BAD=.TRUE.
                    NTSGOOD=NTSGOOD-1
                    CYCLE
                   END IF
                 END IF
               END IF
               IF ((CHECKCHIRALT.OR.SETCHIRAL).AND.AMBER12T) THEN
                  GOODSTRUCTURE1=.TRUE.
                  GOODSTRUCTURE2=.TRUE.
                  CALL CHIRALITY_CHECK(QMINUS, GOODSTRUCTURE1)
                  CALL CHIRALITY_CHECK(QPLUS, GOODSTRUCTURE2)
                  IF ((.NOT.GOODSTRUCTURE1).OR.(.NOT.GOODSTRUCTURE2)) THEN
                    WRITE(*,'(A)') ' tryconnect> Chirality inversion detected in at least one of the carbon centres, rejecting'
                    DEALLOCATE(TS(I)%DATA%EVALMIN)
                    DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                    TS(I)%DATA%BAD=.TRUE.
                    NTSGOOD=NTSGOOD-1
                    CYCLE
                  ENDIF
               ENDIF
! AMBER 9
               IF ((AMBERT.OR.NABT).AND.NOCISTRANS) THEN
                 IF(NOCISTRANSRNA) THEN
                   CALL CHECK_CISTRANS_RNA(QPLUS,NATOMS,ZSYM,GOODSTRUCTURE1)
                   CALL CHECK_CISTRANS_RNA(QMINUS,NATOMS,ZSYM,GOODSTRUCTURE2)
                   IF(.NOT.GOODSTRUCTURE1.OR..NOT.GOODSTRUCTURE2) THEN
                    PRINT '(A,G20.10,A)',' tryconnect> Transition state with energy ',TS(I)%DATA%E, &
  &                                     ' ignored, cis-trans isomerisation detected in the RNA ribose ring.'
                    DEALLOCATE(TS(I)%DATA%EVALMIN)
                    DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                    TS(I)%DATA%BAD=.TRUE.
                    NTSGOOD=NTSGOOD-1
                    CYCLE
                   END IF
                 ELSE IF(NOCISTRANSDNA) THEN
                   CALL CHECK_CISTRANS_DNA(QPLUS,NATOMS,ZSYM,GOODSTRUCTURE1)
                   CALL CHECK_CISTRANS_DNA(QMINUS,NATOMS,ZSYM,GOODSTRUCTURE2)
                   IF(.NOT.GOODSTRUCTURE1.OR..NOT.GOODSTRUCTURE2) THEN
                    PRINT '(A,G20.10,A)',' tryconnect> Transition state with energy ',TS(I)%DATA%E, &
  &                                     ' ignored, cis-trans isomerisation detected in the DNA deoxyribose ring.'
                    DEALLOCATE(TS(I)%DATA%EVALMIN)
                    DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                    TS(I)%DATA%BAD=.TRUE.
                    NTSGOOD=NTSGOOD-1
                    CYCLE
                   END IF
                 ELSE
                   CALL CHECK_CISTRANS_PROTEIN(MI(1)%DATA%X(1:NOPT),NATOMS,GOODSTRUCTURE1,MINOMEGA,CISARRAY2)
                   CALL CHECK_CISTRANS_PROTEIN(QMINUS,NATOMS,GOODSTRUCTURE2,MINOMEGA,CISARRAY1)
                   CISARRAY1=CISARRAY1-CISARRAY2
                   GOODSTRUCTURE1=.TRUE.
                   DO J1=1,NATOMS
                    IF(CISARRAY1(J1)/=0) THEN
                      GOODSTRUCTURE1=.FALSE.
                      WRITE(*,'(A,I6)') ' tryconnect> MINUS minimum: cis-trans isomerisation '// &
  &                                                   'of a peptide bond detected involving atom ', J1
                    END IF
                   END DO

                   CALL CHECK_CISTRANS_PROTEIN(QPLUS,NATOMS,GOODSTRUCTURE2,MINOMEGA,CISARRAY1)
                   CISARRAY1=CISARRAY1-CISARRAY2
                   GOODSTRUCTURE2=.TRUE.
                   DO J1=1,NATOMS
                    IF(CISARRAY1(J1)/=0) THEN
                      GOODSTRUCTURE2=.FALSE.
                      WRITE(*,'(A,I6)') ' tryconnect> PLUS minimum: cis-trans isomerisation '// &
  &                                                   'of a peptide bond detected involving atom ', J1
                    END IF
                   END DO
                   IF(.NOT.GOODSTRUCTURE1.OR..NOT.GOODSTRUCTURE2) THEN
                      WRITE(*,'(A)') ' tryconnect> Cis-trans isomerisation of a peptide bond detected '//&
  &                                                '(wrt. the original structure), rejecting'
                    PRINT '(A,G20.10,A)',' tryconnect> Transition state with energy ',TS(I)%DATA%E, &
  &                                     ' ignored, cis-trans isomerisation detected in one or more peptide bonds.'
                    DEALLOCATE(TS(I)%DATA%EVALMIN)
                    DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                    TS(I)%DATA%BAD=.TRUE.
                    NTSGOOD=NTSGOOD-1
                    CYCLE
                   END IF
                 END IF
               END IF
               IF (CHECKCHIRALT.AND.(AMBERT.OR.NABT)) THEN
                  CALL CHECK_CHIRALITY(QMINUS,NATOMS,GOODSTRUCTURE1)
                  CALL CHECK_CHIRALITY(QPLUS,NATOMS,GOODSTRUCTURE2)
                  IF (.NOT.GOODSTRUCTURE1.OR..NOT.GOODSTRUCTURE2) THEN
                    WRITE(*,'(A)') ' tryconnect> Chirality inversion detected in at least one of the carbon centres, rejecting'
                    DEALLOCATE(TS(I)%DATA%EVALMIN)
                    DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                    TS(I)%DATA%BAD=.TRUE.
                    NTSGOOD=NTSGOOD-1
                    CYCLE
                  ENDIF
               ENDIF
!mo361> Discard transition state if any chiral atoms are inverted compared to finish structure. 
               IF (SETCHIRAL.AND.(AMBERT .OR. NABT)) THEN
                  CALL SET_CHECK_CHIRAL(QPLUS,NATOMS,GOODSTRUCTURE1,CHIARRAY2)
                  CHIARRAY2=CHIARRAY1-CHIARRAY2
                  GOODSTRUCTURE1=.TRUE.
                  DO J1=1,NATOMS
                     IF (CHIARRAY2(J1)/=0) THEN
                        GOODSTRUCTURE1=.FALSE.
                        WRITE(*,'(A,I6)') ' tryconnect> chirality inversion at ', J1
                     ENDIF
                  ENDDO

                  CALL SET_CHECK_CHIRAL(QMINUS,NATOMS,GOODSTRUCTURE2,CHIARRAY2)
                  CHIARRAY2=CHIARRAY1-CHIARRAY2
                  GOODSTRUCTURE1=.TRUE.
                  DO J1=1,NATOMS
                     IF (CHIARRAY2(J1)/=0) THEN
                        GOODSTRUCTURE1=.FALSE.
                        WRITE(*,'(A,I6)') ' tryconnect> chirality inversion at ', J1
                     ENDIF
                  ENDDO
                  IF (.NOT.GOODSTRUCTURE1.OR..NOT.GOODSTRUCTURE2) THEN
                    WRITE(*,'(A)') ' tryconnect> Chirality inversion detected in at least one of the carbon centres, rejecting'
                    DEALLOCATE(TS(I)%DATA%EVALMIN)
                    DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                    TS(I)%DATA%BAD=.TRUE.
                    NTSGOOD=NTSGOOD-1
                    CYCLE
                  ENDIF
               ENDIF

! hk286 - QPLUS already in atom coordinates
! hk286 - QMINUS already in atom coordinates
               AMIDEFAIL=.FALSE.
               IF (CHRMMT .AND. CHECKOMEGAT) &
                  CALL CHECKOMEGA(QPLUS,AMIDEFAIL)
               IF (CHRMMT .AND. CHECKOMEGAT .AND. .NOT.AMIDEFAIL) &
                  CALL CHECKOMEGA(QMINUS,AMIDEFAIL)
               CHIRALFAIL=.FALSE.
               IF (CHRMMT .AND. CHECKCHIRALT) &
                  CALL CHECKCHIRAL(QPLUS,CHIRALFAIL)
               IF (CHRMMT .AND. CHECKCHIRALT .AND. .NOT.CHIRALFAIL) &
                  CALL CHECKCHIRAL(QMINUS,CHIRALFAIL)
               IF (CHRMMT .AND. AMIDEFAIL) THEN
                  PRINT '(A,G20.10,A)',' tryconnect> Transition state with energy ',TS(I)%DATA%E, &
  &                                     ' ignored, cis-trans isomerisation of an amide-bond detected.'
                  DEALLOCATE(TS(I)%DATA%EVALMIN)
                  DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                  TS(I)%DATA%BAD=.TRUE.
                  NTSGOOD=NTSGOOD-1
                  CYCLE
               ELSEIF (CHRMMT .AND. CHIRALFAIL) THEN
                  PRINT '(A,G20.10,A)',' tryconnect> Transition state with energy ',TS(I)%DATA%E, &
  &                                        ' ignored, inversion of a chiral CA center detected.'
                  DEALLOCATE(TS(I)%DATA%EVALMIN)
                  DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                  TS(I)%DATA%BAD=.TRUE.
                  NTSGOOD=NTSGOOD-1
                  CYCLE
               ELSEIF (PATHFAILT) THEN
                    DEALLOCATE(TS(I)%DATA%EVALMIN)
                    DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                    WRITE(*,'(A,I10,G20.10,A)') ' tryconnect> path failed for ts ', I, TS(I)%DATA%E, ', ts ignored' 
                    TS(I)%DATA%BAD=.TRUE.
                    NTSGOOD=NTSGOOD-1
                    CYCLE
               ELSEIF (TS(I)%DATA%E-MIN(EPLUS,EMINUS).GT.MAXMAXBARRIER) THEN
                  PRINT '(2(A,G20.10))',' tryconnect> Transition state with energy ',TS(I)%DATA%E,' ignored, maximum barrier=', &
  &                                      TS(I)%DATA%E-MIN(EPLUS,EMINUS)
                  DEALLOCATE(TS(I)%DATA%EVALMIN)
                  DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                  TS(I)%DATA%BAD=.TRUE.
                  NTSGOOD=NTSGOOD-1
!                 STOP !!! DJW
                  CYCLE
               ELSEIF (TS(I)%DATA%E-MAX(EPLUS,EMINUS).GT.MAXBARRIER) THEN
                  PRINT '(2(A,G20.10))',' tryconnect> Transition state with energy ',TS(I)%DATA%E,' ignored, minimum barrier=', &
  &                                      TS(I)%DATA%E-MAX(EPLUS,EMINUS)
                  DEALLOCATE(TS(I)%DATA%EVALMIN)
                  DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                  TS(I)%DATA%BAD=.TRUE.
                  NTSGOOD=NTSGOOD-1
                  CYCLE
               ELSE
                  IF (UNRST) CALL TESTSAMEMIN(EPLUS,QPLUS,EMINUS,QMINUS,PATHFAILT)
                  IF (PATHFAILT) THEN
                      DEALLOCATE(TS(I)%DATA%EVALMIN)
                      DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                    WRITE(*,*) ' tryconnect1> path failed for ts ', I, TS(I)%DATA%E, ', ts ignored'
                      TS(I)%DATA%BAD=.TRUE.
                      NTSGOOD=NTSGOOD-1
                      CYCLE 
                  ENDIF
               ENDIF

333            CONTINUE
               CALL ISNEWMIN(EPLUS,QPLUS,MINPLUSPOS,PLUSNEW,REDOPATH,PERMUTE,INVERT,INDEX,IMATCH)
               CALL ISNEWMIN(EMINUS,QMINUS,MINMINUSPOS,MINUSNEW,REDOPATH,PERMUTE,INVERT,INDEX,IMATCH)
!
! The above check will not discover the case when the plus minimum is new, and is the same
! as the minus minimum.
!
               IF (PLUSNEW.AND.MINUSNEW) THEN
                  IF (ABS(EMINUS-EPLUS) < EDIFFTOL) THEN
                     CALL ALIGN_DECIDE(QMINUS,QPLUS,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
                     IF (D<GEOMDIFFTOL) THEN
                        MINUSNEW=.FALSE.
                        MINMINUSPOS=MINPLUSPOS
                     ENDIF
                  ENDIF
               ENDIF

               EDUMMY=TS(I)%DATA%E
               TMPTS(1:NOPT)=TS(I)%DATA%X(1:NOPT)
               IF (DUMPALLPATHS) CALL MAKEALLPATHINFO(TMPTS,QPLUS,QMINUS,EDUMMY,EPLUS,EMINUS,FRQSTS,FRQSPLUS,FRQSMINUS)
               
               WRITE(CHR,'(i7)') MinPlusPos
               WRITE(CHR2,'(i7)') MinMinusPos
               100 FORMAT (8X,A,T65,A)
               IF ( .NOT.PLUSNEW .AND. .NOT.MINUSNEW ) THEN
                    WRITE(*,100) 'Known (#'//trim(adjustl(chr))//')','Known (#'//trim(adjustl(chr2))//')'
                    CALL NEWCONNECTION(MINPLUSPOS,MINMINUSPOS,I)
!                   IF (ALLOCATED(QPLUS)) DEALLOCATE(QPLUS)
!                   IF (ALLOCATED(QMINUS)) DEALLOCATE(QMINUS)
!                   DEALLOCATE(QPLUS,QMINUS)
                    CALL SETDISTANCE(MINPLUSPOS,MINMINUSPOS,0.0D0)
                    IF (INTERPCOSTFUNCTION) CALL SETINTERP(MINPLUSPOS,MINMINUSPOS,0.0D0)
               ELSE IF ( PLUSNEW .AND. MINUSNEW ) THEN
                    WRITE(CHR2,'(i7)') MinPlusPos+1
                    WRITE(*,100) '*NEW* (Placed in '//trim(adjustl(chr))//')','*NEW* (Placed in '//trim(adjustl(chr2))//')'

                    CALL ADDNEWMIN(EPLUS,QPLUS)
                    CALL ADDNEWMIN(EMINUS,QMINUS)
                    CALL NEWCONNECTION(MINPLUSPOS,MINPLUSPOS+1,I)
                    MI(MINPLUSPOS+1)%DATA%D(MINPLUSPOS)=0.0D0
                    IF (INTERPCOSTFUNCTION) MI(MINPLUSPOS+1)%DATA%INTERP(MINPLUSPOS)=0.0D0
               ELSE IF ( PLUSNEW .OR. MINUSNEW ) THEN
                    IF ( PLUSNEW ) THEN
                         WRITE(*,100) '*NEW* (Placed in '//trim(adjustl(chr))//')','Known (#'//trim(adjustl(chr2))//')'
                         CALL ADDNEWMIN(EPLUS,QPLUS)
                         DEALLOCATE(EMINUS,QMINUS)
                    ELSE
                         WRITE(*,100) 'Known (#'//trim(adjustl(chr))//')','*NEW* (Placed in '//trim(adjustl(chr2))//')'
                         DEALLOCATE(EPLUS,QPLUS)
                         CALL ADDNEWMIN(EMINUS,QMINUS)
                    ENDIF
                    CALL NEWCONNECTION(MINPLUSPOS,MINMINUSPOS,I)
                    CALL SETDISTANCE(MINPLUSPOS,MINMINUSPOS,0.0D0)
                    IF (INTERPCOSTFUNCTION) CALL SETINTERP(MINPLUSPOS,MINMINUSPOS,0.0D0)
               ENDIF
               IF (DOAGAIN(I-NTS+UNIQUE).AND.(.NOT.REDOPATH)) THEN
                  PUSHOFF=SAVEPUSHOFF
               ENDIF
          ENDDO
        ENDIF  

          !IF(ALLOCATED(FLATPATHT)) DEALLOCATE(FLATPATHT)
!
! Allow for new pathway calculation with different PUSHOFF and MAXBFGS
!
          IF (REDOPATH.AND.(.NOT.PATHFAILT)) THEN
             CALL ALIGN_DECIDE(QP,MIN1REDO,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
             DIST1P=D

             CALL ALIGN_DECIDE(QM,MIN1REDO,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
             DIST1M=D

             CALL ALIGN_DECIDE(QM,MIN2REDO,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
             DIST2M=D

             CALL ALIGN_DECIDE(QP,MIN2REDO,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
             DIST2P=D

             PATHFAILED=.FALSE.
             IF ((DIST1P.GT.GEOMDIFFTOL).AND.(DIST1M.GT.GEOMDIFFTOL)) THEN
                PRINT '(A)',' tryconnect> path failed to match first minimum'
                PATHFAILED=.TRUE.
             ENDIF
             IF ((DIST2P.GT.GEOMDIFFTOL).AND.(DIST2M.GT.GEOMDIFFTOL)) THEN
                PRINT '(A)',' tryconnect> path failed to match second minimum'
                PATHFAILED=.TRUE.
             ENDIF
             IF (PATHFAILED) THEN
                NC1=NC1+1
                IF (PUSHOPTT.OR.(NC1.GT.2*NREDOPATHTRIES1)) THEN
                   NC2=NC2+1
                   NC1=0
                ENDIF
                IF (NREDOPATHTRIES1.EQ.0) THEN
                   PUSHOFF=SAVEPUSHOFF
                ELSEIF (NC1.GT.NREDOPATHTRIES1) THEN
                   PUSHOFF=SAVEPUSHOFF+(NC1-NREDOPATHTRIES1)*(REDOSTRETCH*SAVEPUSHOFF-SAVEPUSHOFF)/NREDOPATHTRIES1
                ELSE
                   PUSHOFF=SAVEPUSHOFF-NC1*(SAVEPUSHOFF-SAVEPUSHOFF/REDOSTRETCH)/NREDOPATHTRIES1
                ENDIF
                IF (NREDOPATHTRIES2.EQ.0) THEN
                   MAXBFGS=SAVEMAXBFGS
                ELSEIF (NC2.GT.NREDOPATHTRIES2) THEN
                   MAXBFGS=SAVEMAXBFGS+(NC2-NREDOPATHTRIES2)*(REDOSTRETCH*SAVEMAXBFGS-SAVEMAXBFGS)/NREDOPATHTRIES2
                ELSE
                   MAXBFGS=SAVEMAXBFGS-NC2*(SAVEMAXBFGS-SAVEMAXBFGS/REDOSTRETCH)/NREDOPATHTRIES2
                ENDIF
                IF (NC2.LE.2*NREDOPATHTRIES2) THEN
                   PRINT '(2(A,F15.5))',' tryconnect> Redo path with pushoff=',PUSHOFF,' maxbfgs=',MAXBFGS
                   RERUN=.TRUE.
                   GOTO 10
                ENDIF
             ENDIF
           ENDIF

           RERUN=.FALSE.
           PUSHOFF=SAVEPUSHOFF
           MAXBFGS=SAVEMAXBFGS
           DEALLOCATE(FOUNDBEFORE,DOAGAIN)
 
     END SUBROUTINE TRYCONNECT

END MODULE TRYCONNECTMODULE
