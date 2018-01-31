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

     SUBROUTINE TRYCONNECT(JS,JF,TSREDO,REDOPATH,REDOPATHXYZ)
          USE NEWNEBMODULE
          USE CONNECTDATA
          USE KEYCONNECT
          USE CONNECTUTILS
          USE NEBTOCONNECT
          USE KEYNEB, ONLY: NIMAGE, NITERMAX
          USE KEY, ONLY: UNRST, DEBUG, FILTH, FILTHSTR, DUMPALLPATHS, TWOD, MAXTSENERGY, RIGIDBODY, &
  &                   PERMDIST, MAXBARRIER, GROWSTRINGT, DESMDEBUG, BULKT
          USE MODGUESS
          USE MODUNRES
          USE MODCHARMM, ONLY: CHRMMT,NCHENCALLS,CHECKOMEGAT,CHECKCHIRALT
          USE MODMEC
          USE KEYUTILS
          USE COMMONS, ONLY: NINTS, PARAM1, PARAM2, PARAM3
          USE PORFUNCS
          USE GSDATA, ONLY: GSITERDENSITY, GSCURITERD=>ITERD
          IMPLICIT NONE
          DOUBLE PRECISION RMAT(3,3)

          INTEGER,INTENT(IN) :: JS,JF
          
          INTEGER         :: I,UNIQUE=0,MINPLUSPOS,MINMINUSPOS,J1
          DOUBLE PRECISION         :: EDUMMY,EDUMMY2,TMPTS(3*NATOMS)
          DOUBLE PRECISION,POINTER :: QPLUS(:),QMINUS(:),EPLUS,EMINUS
          LOGICAL         :: PLUSNEW,MINUSNEW,PATHFAILT,AMIDEFAIL,CHIRALFAIL
          CHARACTER       :: ITSTRING*80, EOFSSTRING*80
          CHARACTER       :: ITSTRINGP*80, EOFSSTRINGP*80
          CHARACTER       :: ITSTRINGM*80, EOFSSTRINGM*80
          DOUBLE PRECISION :: STARTINT(NINTS), FINISHINT(NINTS),DUM
          LOGICAL REDOPATH, REDOPATHXYZ, PERMUTE, EXISTS, NORERUN, PERMONCE
          DOUBLE PRECISION TSREDO(3*NATOMS), QTEMP(3*NATOMS), QLOCAL(3*NATOMS), ETSPREV, ETSDUM, DIST2
          CHARACTER(LEN=5) ZTEMP(NATOMS)
          CHARACTER(LEN=80) DSTRING, TMPSTRING
          INTEGER INVERT, INDEX(NATOMS), J2, IMATCH
          CHARACTER(LEN=2) ZDUM
          CHARACTER(LEN=10) DUMMYS          

          IF (GROWSTRINGT) THEN
             IF (NCONDONE.EQ.1.AND.FCD) THEN
                GSCURITERD = GSITERDENSITY
             ELSE
                GSCURITERD = ITERDENSITY
             ENDIF
          ENDIF

          IF (CHRMMT) NCHENCALLS = 999 ! update non-bonded list on next call to potential.

          IF (.NOT.REDOPATH) CALL CHECKPAIR(JS,JF)

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
               CALL GUESSPATH(MI(JS)%DATA%X,MI(JF)%DATA%X,3*NATOMS,EDIFFTOL,NATOMS)
            ENDIF
          ! how many images to use?
            IF (NIMAGE > IMAGEMAX) PRINT*,'WARNING - Nimage is greater than ImageMax'
            IF (NIMAGE < 2       ) PRINT*,'WARNING - Nimage is < 2'
            NIMAGE=NINTERP
!           NIterMax = Nimage*IterDensity ! try zero neb iterations if we have a GUESSPATH path
            NITERMAX = 0
            IF (NINTERP.LT.2) THEN ! NO IMAGES FROM GUESSPATH - REVERT TO USUAL SCHEME
               IF (.NOT.(NCONDONE==1 .AND. FCD)) THEN
                    NIMAGE=IMAGEDENSITY*MI(JF)%DATA%D(JS) &
                        +IMAGEINCR*IMAGEDENSITY*MI(JF)%DATA%D(JS)*(MI(JF)%DATA%NTRIES(JS)-1)
                    IF (NIMAGE >= IMAGEMAX) THEN
                       NIMAGE = IMAGEMAX
!                      mi(jf)%data%ntries(js)=NTriesMax ! no point trying again with the same number of images
                    ENDIF
                    IF (NIMAGE < 2       ) NIMAGE = 2
                    NITERMAX = NIMAGE*ITERDENSITY
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

          ELSEIF (REDOPATH) THEN
            NIMAGE=1
            NITERMAX = 0 ! TRY ZERO NEB ITERATIONS IF WE HAVE POINTS IN TSREDO
          ELSEIF (NCONDONE==1 .AND. FCD) THEN ! FIRST CYCLE DIFFERENT - PARAMETERS SUPPLIED USING NEWNEB
                                          ! keyword or newneb defaults will be used instead
               PRINT *, "First cycle will be done using externally supplied parameters"
          ELSE
            NIMAGE=IMAGEDENSITY*MI(JF)%DATA%D(JS) &
             +IMAGEINCR*IMAGEDENSITY*MI(JF)%DATA%D(JS)*(MI(JF)%DATA%NTRIES(JS)-1)
               IF (NIMAGE >= IMAGEMAX) THEN
                  NIMAGE = IMAGEMAX
!                 mi(jf)%data%ntries(js)=NTriesMax ! no point trying again with the same number of images
               ENDIF
               IF (NIMAGE < 2       ) NIMAGE = 2
               NITERMAX = NIMAGE*ITERDENSITY
          ENDIF
          
          ! book-keeping :-)
          IF (.NOT.(MECCANOT.OR.REDOPATH)) THEN             
             IF (GROWSTRINGT) THEN
                WRITE(CHR,'(i7)') INT(NIMAGE*GSCURITERD)
                WRITE(*,'(/1x,a)',advance='no') '>>>>>  '//trim(adjustl(chr))//'-iteration GS run for minima '
             ELSE
                WRITE(CHR,'(i7)') NIterMax
                WRITE(*,'(/1x,a)',advance='no') '>>>>>  '//trim(adjustl(chr))//'-iteration DNEB run for minima '
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
             WRITE(CHR,'(i5)') Nimage
             WRITE(*,'(a)',advance='no') ' using '//trim(adjustl(chr))//' images '
             IF (MI(JF)%DATA%NTRIES(JS) > 1) THEN
                  WRITE(CHR,'(i5)') mi(jf)%data%ntries(js)
                  WRITE(*,'(a)',advance='no') '(attempt #'//trim(adjustl(chr))//') '
             ENDIF
             WRITE(*,'(a)') ' ...'
          ENDIF

          IF (NIMAGE >= IMAGEMAX) MI(JF)%DATA%NTRIES(JS)=NTRIESMAX ! NO POINT TRYING AGAIN WITH THE SAME NUMBER OF IMAGES
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
             IF (EXISTS) THEN ! ALLOWS FOR RERUN WITH DIFFERENT ENERGY DIFFERENCE CRITERION FOR
                              ! consecutive frames in path without redoing original path
                PRINT '(2A)',' tryconnect> Reading data for ts from existing file ',TRIM(ADJUSTL(itstring))
                NORERUN=.TRUE.
             ELSE
                PRINT '(3A)',' tryconnect> No file ',TRIM(ADJUSTL(itstring)),' found'
                REDOPATH=.FALSE.
                REDOPATHXYZ=.FALSE.
                STOP ! REDOPATHXYZ IS ONLY GOING TO WORK AFTER A RUN FROM REDOPOINTS, IN WHICH CASE
                     ! we should finish after reading the last path.
             ENDIF
             IF (NORERUN) THEN
                NTSFOUND=1
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
                IF (UNRST) THEN 
                   CALL NEWNEB(REDOPATH,TSREDO,-1.0D100,MI(JS)%DATA%X,-1.0D100,MI(JF)%DATA%X,NINTSIN=NINTS)
                ELSE
                   CALL NEWNEB(REDOPATH,TSREDO,-1.0D100,MI(JS)%DATA%X,-1.0D100,MI(JF)%DATA%X)
                ENDIF
             ENDIF
          ELSE
             IF (UNRST) THEN 
                 CALL NEWNEB(REDOPATH,TSREDO,MI(JS)%DATA%E,MI(JS)%DATA%X,MI(JF)%DATA%E,MI(JF)%DATA%X,NINTSIN=NINTS)
             ELSE
                 CALL NEWNEB(REDOPATH,TSREDO,MI(JS)%DATA%E,MI(JS)%DATA%X,MI(JF)%DATA%E,MI(JF)%DATA%X)
             ENDIF
          ENDIF

          ! saving new ts into ts rack; otherwise - free memory immediately
          NTSOLD=NTS
          UNIQUE=0
          DO I=1,NTSFOUND
!              PRINT '(A,2G20.10)',' tryconnect> TSfound(i)%E, MAXTSENERGY=',TSfound(i)%E, MAXTSENERGY
               AMIDEFAIL=.FALSE.
               IF (CHRMMT .AND. CHECKOMEGAT) &
                  CALL CHECKOMEGA(TSFOUND(I)%COORD,AMIDEFAIL)
               CHIRALFAIL=.FALSE.
               IF (CHRMMT .AND. CHECKCHIRALT) &
                  CALL CHECKCHIRAL(TSFOUND(I)%COORD,CHIRALFAIL)
               IF (CHRMMT .AND. AMIDEFAIL) THEN
                 PRINT '(A,G20.10,A)',' tryconnect> Transition state with energy ',TSfound(i)%E, &
  &                                   ' ignored, cis-trans isomerisation of an amide bond detected.'
                 DEALLOCATE(TSFOUND(I)%E,TSFOUND(I)%COORD,TSFOUND(I)%EVALMIN,TSFOUND(I)%VECS)
               ELSEIF (CHRMMT .AND. CHIRALFAIL) THEN
                  PRINT '(A,G20.10,A)',' tryconnect> Transition state with energy ',TSfound(i)%E, &
  &                                    ' ignored, inversion of a chiral CA center detected.'
                  DEALLOCATE(TSFOUND(I)%E,TSFOUND(I)%COORD,TSFOUND(I)%EVALMIN,TSFOUND(I)%VECS)
               ELSEIF (TSFOUND(I)%E.GT.MAXTSENERGY) THEN
                  PRINT '(A,G20.10,A)',' tryconnect> Transition state with energy ',TSfound(i)%E,' ignored'
                  DEALLOCATE(TSFOUND(I)%E,TSFOUND(I)%COORD,TSFOUND(I)%EVALMIN,TSFOUND(I)%VECS)
               ELSEIF ( ISNEWTS(TSFOUND(I)) ) THEN
                    IF (NTS==TSRACKSIZE) CALL REALLOCATETSRACK
                    NTS=NTS+1; UNIQUE=UNIQUE+1
                    TS(NTS)%DATA%E => TSFOUND(I)%E
                    TS(NTS)%DATA%X => TSFOUND(I)%COORD
                    TS(NTS)%DATA%EVALMIN => TSFOUND(I)%EVALMIN
                    TS(NTS)%DATA%VECS => TSFOUND(I)%VECS
                    TS(NTS)%DATA%BAD=.FALSE.
                    NULLIFY(TSFOUND(I)%E,TSFOUND(I)%COORD,TSFOUND(I)%EVALMIN,TSFOUND(I)%VECS)
               ELSE
                    IF (NCONDONE==1) PRINT *, 'Discarded TS #',i
                    DEALLOCATE(TSFOUND(I)%E,TSFOUND(I)%COORD,TSFOUND(I)%EVALMIN,TSFOUND(I)%VECS)
               ENDIF
          ENDDO

          CALL DUMPTS

          ! print info as to how many TS are actually useful
          IF (UNIQUE==NTSFOUND.AND..NOT.UNIQUE==0) THEN
               IF (NTSFOUND==1) THEN
                    WRITE(*,'(1x,a)') 'TS appears to be new'
               ELSE
                    WRITE(*,'(1x,a)') 'All of TS found appear to be new'
               ENDIF
          ELSEIF (UNIQUE < NTSFOUND) THEN
               WRITE(CHR,'(i7)') unique 
               WRITE(*,'(1x,a)') trim(adjustl(chr))//' of TS found appear to be new.'
          ELSEIF (UNIQUE ==0 .AND..NOT.NTSFOUND==0) THEN
               WRITE(*,'(1x,a)') 'All of TS found are already known'
          ENDIF

          ! path run for all unique ts
          DO I=NTS-UNIQUE+1,NTS
               WRITE(CHR,'(i5)') i
               PRINT '(/1x,a)', '>>>>>  Path run for ts '//trim(adjustl(chr))//' ...'
               ALLOCATE( QPLUS(NOPT),QMINUS(NOPT),EPLUS,EMINUS )
               CALL MKFNAMES(I,FILTH,FILTHSTR,ITSTRING,EOFSSTRING)
               EDUMMY=TS(I)%DATA%E
               TMPTS=TS(I)%DATA%X
               ! structure in ts(i)%data%X is a stationary point which is why we don't need to store G and rms for it
               GDUMMY(1:3*NATOMS)=0.0D0; RMS=0.0D0 ! WE MUST INITIALIZE THEM HERE, HOWEVER 
               NORERUN=.FALSE.
               IF (REDOPATH) THEN
                  CALL MKFNAMES(NCONDONE,FILTH,FILTHSTR,ITSTRING,EOFSSTRING)
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
                  CALL PATH(TMPTS,EDUMMY,GDUMMY,RMS,TS(I)%DATA%EVALMIN,TS(I)%DATA%VECS,  &
                    & .FALSE.,QPLUS,QMINUS,DEBUG,EDUMMY2,EPLUS,EMINUS, &
                    & TS(I)%DATA%SLENGTH,TS(I)%DATA%DISP,TS(I)%DATA%GAMMA,TS(I)%DATA%NTILDE,FRQSTS,FRQSPLUS, &
                    & FRQSMINUS,ITSTRING,EOFSSTRING,PATHFAILT)
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
                  PRINT '(A,G20.10,A,G20.10)',' tryconnect> E+=',Eplus,'                      E-=',Eminus
                  CLOSE(89)
               ENDIF
               DEALLOCATE(TS(I)%DATA%VECS)
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
                  CYCLE
               ELSEIF (CHRMMT .AND. CHIRALFAIL) THEN
                  PRINT '(A,G20.10,A)',' tryconnect> Transition state with energy ',TS(I)%DATA%E, &
  &                                        ' ignored, inversion of a chiral CA center detected.'
                  DEALLOCATE(TS(I)%DATA%EVALMIN)
                  DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                  TS(I)%DATA%BAD=.TRUE.
                  CYCLE
               ELSEIF (PATHFAILT) THEN
                    DEALLOCATE(TS(I)%DATA%EVALMIN)
                    DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                    TS(I)%DATA%BAD=.TRUE.
                    CYCLE
               ELSEIF (TS(I)%DATA%E-MAX(EPLUS,EMINUS).GT.MAXBARRIER) THEN
                  PRINT '(2(A,G20.10))',' tryconnect> Transition state with energy ',TS(I)%DATA%E,' ignored, minimum barrier=', &
  &                                      TS(I)%DATA%E-MAX(EPLUS,EMINUS)
                  DEALLOCATE(TS(I)%DATA%EVALMIN)
                  DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                  TS(I)%DATA%BAD=.TRUE.
                  CYCLE
               ELSE
                  IF (UNRST) CALL TESTSAMEMIN(EPLUS,QPLUS,EMINUS,QMINUS,PATHFAILT)
                  IF (PATHFAILT) THEN
                      DEALLOCATE(TS(I)%DATA%EVALMIN)
                      DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
                      TS(I)%DATA%BAD=.TRUE.
                      CYCLE 
                  ENDIF
               ENDIF

!              PERMONCE=.FALSE. ! A degenerate rearrangement of a permuational isomer could otherwise
!                               ! cause an infinite loop.
333            CALL ISNEWMIN(EPLUS,QPLUS,MINPLUSPOS,PLUSNEW,REDOPATH,PERMUTE,INVERT,INDEX,IMATCH)
!              IF ((.NOT.(PERMUTE.AND.REDOPATH)).OR.PERMONCE.OR.REDOPATHXYZ) THEN
                  CALL ISNEWMIN(EMINUS,QMINUS,MINMINUSPOS,MINUSNEW,REDOPATH,PERMUTE,INVERT,INDEX,IMATCH)
!              ENDIF

!                IF (PERMUTE.AND.REDOPATH.AND.(.NOT.PERMONCE).AND.(.NOT.REDOPATHXYZ)) THEN
!                   PERMONCE=.TRUE.
!                   IF (IMATCH.EQ.2) THEN ! permute the final minimum instead!
!                      PRINT '(A)',' tryconnect> Permuting the finish minimum to align it'
!                      DO J2=1,NATOMS
!                         QTEMP(3*(J2-1)+1)=INVERT*MI(IMATCH)%DATA%X(3*(INDEX(J2)-1)+1)
!                         IF (TWOD) THEN
!                            QTEMP(3*(J2-1)+2)=MI(IMATCH)%DATA%X(3*(INDEX(J2)-1)+2)
!                            QTEMP(3*(J2-1)+3)=MI(IMATCH)%DATA%X(3*(INDEX(J2)-1)+3)
!                         ELSE
!                            QTEMP(3*(J2-1)+2)=INVERT*MI(IMATCH)%DATA%X(3*(INDEX(J2)-1)+2)
!                            QTEMP(3*(J2-1)+3)=INVERT*MI(IMATCH)%DATA%X(3*(INDEX(J2)-1)+3)
!                         ENDIF
!                      ENDDO
!                      MI(IMATCH)%DATA%X(1:NOPT)=QTEMP(1:NOPT)
! !
! ! recalculate all distances 
! !
!                      IF (PERMDIST) THEN
!                         CALL MINPERMDIST(MI(1)%DATA%X,MI(2)%DATA%X, NATOMS, &
!   &                                DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
!                         D=SQRT(D)
!                      ELSE
!                         CALL NEWMINDIST(MI(1)%DATA%X,MI(2)%DATA%X,NATOMS,D,.FALSE.,.FALSE.,'AX   ',.True.,RIGIDBODY,DEBUG,RMAT)
!                      ENDIF
!                      MI(2)%DATA%D(1)=D
!                      IF (INTERPCOSTFUNCTION) MI(2)%DATA%INTERP(1)=INTERPVALUE(MI(1)%DATA%X(:),MI(2)%DATA%X(:),D)
!                      DO J2=3,NMIN
!                         CALL NEWMINDIST(MI(J2)%DATA%X,MI(2)%DATA%X,NATOMS,D,.FALSE.,.FALSE.,'AX   ',.True.,RIGIDBODY,DEBUG,RMAT)
!                         MI(J2)%DATA%D(2)=D
!                         IF (INTERPCOSTFUNCTION) MI(J2)%DATA%INTERP(2)=INTERPVALUE(MI(J2)%DATA%X(:),MI(2)%DATA%X(:),D)
!                      ENDDO
!                      GOTO 333
!                   ELSE
!                      PRINT '(A)',' tryconnect> Permuting the pathway to align it'
!                      DO J2=1,NATOMS
!                         QTEMP(3*(J2-1)+1)=INVERT*TS(I)%DATA%X(3*(INDEX(J2)-1)+1)
!                         IF (TWOD) THEN
!                            QTEMP(3*(J2-1)+2)=TS(I)%DATA%X(3*(INDEX(J2)-1)+2)
!                            QTEMP(3*(J2-1)+3)=TS(I)%DATA%X(3*(INDEX(J2)-1)+3)
!                         ELSE
!                            QTEMP(3*(J2-1)+2)=INVERT*TS(I)%DATA%X(3*(INDEX(J2)-1)+2)
!                            QTEMP(3*(J2-1)+3)=INVERT*TS(I)%DATA%X(3*(INDEX(J2)-1)+3)
!                         ENDIF
!                      ENDDO
!                      TS(I)%DATA%X(1:NOPT)=QTEMP(1:NOPT)
!                      DO J2=1,NATOMS
!                         QTEMP(3*(J2-1)+1)=INVERT*QPLUS(3*(INDEX(J2)-1)+1)
!                         IF (TWOD) THEN
!                            QTEMP(3*(J2-1)+2)=QPLUS(3*(INDEX(J2)-1)+2)
!                            QTEMP(3*(J2-1)+3)=QPLUS(3*(INDEX(J2)-1)+3)
!                         ELSE
!                            QTEMP(3*(J2-1)+2)=INVERT*QPLUS(3*(INDEX(J2)-1)+2)
!                            QTEMP(3*(J2-1)+3)=INVERT*QPLUS(3*(INDEX(J2)-1)+3)
!                         ENDIF
!                      ENDDO
!                      QPLUS(1:NOPT)=QTEMP(1:NOPT)
!                      DO J2=1,NATOMS
!                         QTEMP(3*(J2-1)+1)=INVERT*QMINUS(3*(INDEX(J2)-1)+1)
!                         IF (TWOD) THEN
!                            QTEMP(3*(J2-1)+2)=QMINUS(3*(INDEX(J2)-1)+2)
!                            QTEMP(3*(J2-1)+3)=QMINUS(3*(INDEX(J2)-1)+3)
!                         ELSE
!                            QTEMP(3*(J2-1)+2)=INVERT*QMINUS(3*(INDEX(J2)-1)+2)
!                            QTEMP(3*(J2-1)+3)=INVERT*QMINUS(3*(INDEX(J2)-1)+3)
!                         ENDIF
!                      ENDDO
!                      QMINUS(1:NOPT)=QTEMP(1:NOPT)
! !
! !  We need to fix all the points in the path.n.xyz file as well, or we will have
! !  problems when we try to reconstruct the full path.
! !
!                      OPEN(UNIT=81,FILE=ITSTRING,STATUS='OLD')
!                      TMPSTRING=TRIM(ADJUSTL("temp." // ITSTRING))
!                      OPEN(UNIT=82,FILE=TMPSTRING,STATUS='UNKNOWN')
! 222                  READ(81,'(A80)',END=444) DSTRING
!                      WRITE(82,'(A80)') DSTRING
!                      READ(81,'(A80)') DSTRING
!                      WRITE(82,'(A80)') DSTRING
!                      DO J2=1,NATOMS
!                         READ(81,*) ZTEMP(J2),QTEMP(3*(J2-1)+1),QTEMP(3*(J2-1)+2),QTEMP(3*(J2-1)+3)
!                      ENDDO
!                      DO J2=1,NATOMS
!                         QLOCAL(1)=INVERT*QTEMP(3*(INDEX(J2)-1)+1)
!                         IF (TWOD) THEN
!                            QLOCAL(2)=QTEMP(3*(INDEX(J2)-1)+2)
!                            QLOCAL(3)=QTEMP(3*(INDEX(J2)-1)+3)
!                         ELSE
!                            QLOCAL(2)=INVERT*QTEMP(3*(INDEX(J2)-1)+2)
!                            QLOCAL(3)=INVERT*QTEMP(3*(INDEX(J2)-1)+3)
!                         ENDIF
!                         WRITE(82,'(A2,4X,3F20.10)') ZTEMP(J2),QLOCAL(1),QLOCAL(2),QLOCAL(3)
!                      ENDDO
!                      GOTO 222
! 444                  CLOSE(81)
!                      CLOSE(82)
!                      CALL SYSTEM('mv ' // TMPSTRING // ' ' // ITSTRING)
!                      GOTO 333
!                   ENDIF
!                ENDIF

               EDUMMY=TS(I)%DATA%E
               TMPTS=TS(I)%DATA%X
               IF (DUMPALLPATHS) CALL MAKEALLPATHINFO(TMPTS,QPLUS,QMINUS,EDUMMY,EPLUS,EMINUS,FRQSTS,FRQSPLUS,FRQSMINUS)
               
               WRITE(CHR,'(i7)') MinPlusPos
               WRITE(CHR2,'(i7)') MinMinusPos
               100 FORMAT (8X,A,T65,A)
               IF ( .NOT.PLUSNEW .AND. .NOT.MINUSNEW ) THEN
                    WRITE(*,100) 'Known (#'//trim(adjustl(chr))//')','Known (#'//trim(adjustl(chr2))//')'
                    CALL NEWCONNECTION(MINPLUSPOS,MINMINUSPOS,I)
                    DEALLOCATE(QPLUS,QMINUS,EPLUS,EMINUS)
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
          ENDDO

     END SUBROUTINE TRYCONNECT

END MODULE TRYCONNECTMODULE
