!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999- David J. Wales
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
! MC trial move is all atom random cartesian perturbation.
! We read geometries from file path.xyz. 
! Some of this sampling is based on GMIN bspt routine.
!
! IMPORTANT VARIABLES
!
! XO, YO, ZO, VOLD = the old coordinates and energy. i.e. the saved markov state
! at the end of the previous step
!
! X, Y, Z, VNEW    = the new perturbed coordinates and energy.  If the step is
! rejected, these are set back to XO, YO, ZO, VOLD, so at the end of the
! iteration these are the state of the markov chain
!
! RECOUNT  =  the boolean which records whether the step is accepted or
! rejected.  If the step is rejected, RECOUNT==.TRUE., and the previous markov
! state is recounted.
!
! This is version 2 of mcpath, intended for use with short overlapping blocks
! corresponding to configurations in the path.xyz file. No bias potential this
! time. Or perhaps just use the pe of the configuration that is closest since we
! are tracking this, and rejecting moves outside the current set of frames.
!

SUBROUTINE MCPATH2
USE KEY
USE COMMONS
USE PORFUNCS
IMPLICIT NONE
INTEGER J1, J2, NSTRUCTREF, HIGHESTREF, LUNIT, K1, K2, J3, NNODES, I, ISTAT
CHARACTER(LEN=80) :: ARG
CHARACTER(LEN=80),ALLOCATABLE,DIMENSION(:) :: NODENAME
CHARACTER(LEN=80) :: USERNAME
CHARACTER(LEN=100) :: WORKINGDIRECTORY
INTEGER GETUNIT, IBININDEX, NDUMMY, NPATH, NTRIES
INTEGER STARTSAMPLE, ENDSAMPLE, NAPPEND
INTEGER DOBLOCK, MCPATHCPU
INTEGER, ALLOCATABLE :: NEAREST(:), NEARESTF(:), PID(:)
DOUBLE PRECISION, ALLOCATABLE :: BINLABELQORDER(:), NEARESTFW(:), SLENGTH(:), QFRAME(:), QFRAMEAV(:), &
  &                              QFRAMESD(:), QFBEST(:)
DOUBLE PRECISION, ALLOCATABLE :: EOFS(:), PATHFRAMES(:,:), LEOFS(:)
DOUBLE PRECISION, ALLOCATABLE :: QORDERHIST(:), QORDERVISITS(:), NORM1(:)
DOUBLE PRECISION, ALLOCATABLE :: SMOOTH2D(:,:), SMOOTH2D2(:,:)
DOUBLE PRECISION, ALLOCATABLE :: COORDSREF(:,:), VREF(:), LCOORDSREF(:,:), LVREF(:)
DOUBLE PRECISION, ALLOCATABLE :: PROBS(:), PROBQ(:), PROBSQ(:,:), PROBSFAKE(:), ZNORM(:), PROBSPAIR(:), ZNORMPAIR(:)
DOUBLE PRECISION PROBQFROMS(MCPATHBINS), SMOOTHQFROMS(MCPATHBINS), PROBQFROMS2(MCPATHBINS), SMOOTHQFROMS2(MCPATHBINS)
DOUBLE PRECISION PROBQFROMSPAIR(MCPATHBINS), SMOOTHQFROMSPAIR(MCPATHBINS), PROBQFROMSPAIR2(MCPATHBINS), SMOOTHQFROMSPAIR2(MCPATHBINS)
DOUBLE PRECISION, ALLOCATABLE :: SMOOTHS(:), SMOOTHQ(:), SMOOTHSQ(:,:)
DOUBLE PRECISION, ALLOCATABLE :: SMOOTHSPUP(:), SMOOTHSPDOWN(:),SMOOTHSPFIXUP(:),SMOOTHSPFIXDOWN(:)
DOUBLE PRECISION GRAD(3*NATOMS), RMS, DUMMY, DIST2, RMAT(3,3), CDUMMY(3*NATOMS), &
  &              DUMMY3, DUMMY4, DUMMY2, QORDERHISTINT, SLEEPTIME1, ZSAVE, DUMMY6, DUMMY7, &
  &              NORM2(MCPATHBINS), DUMMY1, DUMMY5, GWIDTHS, GWIDTHQ, DUMMYS, DUMMYQ, DPRAND
DOUBLE PRECISION, ALLOCATABLE :: DSTRUCT(:)
LOGICAL, ALLOCATABLE :: ALLZEROS(:), ALLZEROQ(:), ALLZEROSTOT(:)
LOGICAL YESNO, LASTMIN, LASTTS, PERMDISTSAVE, LPERMDISTSAVE, REDOT, KILLED, LDEBUG, LINVART, OFAILS, OFAILQ, FIXZ
CHARACTER (LEN=8) SDUMMY
CHARACTER(LEN=80) REPSTRING, FSTRING
CHARACTER (LEN=80) FNAME
CHARACTER (LEN=5) SSYM
CHARACTER(LEN=10) PATHSTR, APPSTRING
DOUBLE PRECISION QTEMP(3,NATOMS),QORDER
INTEGER NRUNS, PIDDONE, STATUS, NRUNNING, NEWJOB
CHARACTER(LEN=80) SLEEPSTRING1
CHARACTER(LEN=10) CONNSTR

GRAD(1:3*NATOMS)=0.0D0
SLEEPTIME1=1.0D0
LDEBUG=.TRUE.
WRITE(SLEEPSTRING1,'(A,F20.10)') 'sleep ',SLEEPTIME1
RMS=0.0D0
PERMDISTSAVE=PERMDIST 
LPERMDISTSAVE=LPERMDIST
LINVART=.FALSE. ! linear variables can give negative probabilities!

PRINT '(A,G20.10)','mcpath2> Monte Carlo sampling around path.xyz file, canonical temperature=',MCPATHTEMP
PRINT '(A,G20.10)','mcpath2> Distance limit on sampling from configurations on the pathway= ',MCPATHDMAX
PRINT '(A,G20.10)','mcpath2> Maximum initial MC Cartesian step size=                          ',MCPATHSTEP
PRINT '(A,G20.10)','mcpath2> Target acceptance ratio=                                         ',MCPATHACCRATIO
PRINT '(A,I10)','mcpath2> Number of bins for Q histograms=                          ',MCPATHBINS
PRINT '(A,I10,A,F15.5,A)','mcpath2> Number of equilibration steps=                                   ',MCPATHEQUIL, &
  &                     ' = ',MCPATHEQUIL/1.0D6,' Ms'
PRINT '(A,I10,A,F15.5,A)','mcpath2> Number of MC production sampling steps=                          ',MCPATHSTEPS, &
  &                     ' = ',MCPATHSTEPS/1.0D6,' Ms'
PRINT '(A,I10)','mcpath2> Print frequency interval=                                        ',MCPATHPRTFRQ
PRINT '(A,G20.10)','mcpath2> Bottom of order parameter range=                                 ',MCPATHQMIN
PRINT '(A,G20.10)','mcpath2> Top of order parameter range=                                    ',MCPATHQMAX
PRINT '(A,G20.10)','mcpath2> Number of consecutive frames in each sample                      ',MCPATHBLOCK
PRINT '(A,G20.10)','mcpath2> Number of path.xyz frames in block overlaps                      ',MCPATHOVER
IF (MCPATHDOBLOCK.GT.0) PRINT '(A,G20.10)','mcpath2> Single run for block                                             ', &
  & MCPATHDOBLOCK 
PRINT '(A,G20.10)','mcpath2> Threshold fraction for neglect of bins based upon visits         ',MCPATHNEGLECT
PRINT '(A,G20.10)','mcpath2> Converge condition for WHAM LBFGS chi^2 fit                      ',MCPATHTOL
IF (NCPU.GT.0) PRINT '(A,G20.10)','mcpath2> Running on local node, number of cores=                          ',NCPU
IF (PBST) PRINT '(A,G20.10)','mcpath2> Distributing OPTIM jobs to blocks using pbs'
IF (MCBIAST) THEN
   PRINT '(A)', 'mcpath2> Bias function will be constructed and applied'
ENDIF
IF (MCPATHOVER.GE.MCPATHBLOCK) THEN
   PRINT '(A)','mcath2> ERROR *** overlap must be less than block size'
   STOP
ENDIF
!
! Parse path.xyz to discover the number of frames saved.
! Identify the number of maxima and minima and set NSTRUCTREF.
! Read the reference coordinates and energies into
! COORDSREF and VREF arrays after declaring them.
! Set integer highest to the highest energy reference structure.
!
INQUIRE(FILE='path.xyz',EXIST=YESNO)
IF (.NOT.YESNO) THEN
   PRINT '(A)','mcpath2> ERROR *** no path.xyz file'
   STOP
ENDIF
LUNIT=GETUNIT()
OPEN(LUNIT,FILE='path.xyz',STATUS='OLD')
NPATH=0
DO 
   READ(LUNIT,*,END=20) SDUMMY
   READ(LUNIT,*) SDUMMY
   READ(LUNIT,*) (SSYM,CDUMMY(1:3),J1=1,NATOMS)
   NPATH=NPATH+1
ENDDO
20 PRINT '(A,I10)','mcpath2> Number of frames in path.xyz=',NPATH
NRUNS=(NPATH-MCPATHOVER)/(MCPATHBLOCK-MCPATHOVER)
IF (NRUNS.EQ.1) THEN
   PRINT '(A)','mcpath2> Running a single replica - P(Q) from direct q bin visits is valid'
ELSE
   PRINT '(A)','mcpath2> Running multiple replicas/windows - P(Q) from direct q bin visits will contain systematic errors'
ENDIF
IF (MOD(NPATH-MCPATHOVER,MCPATHBLOCK-MCPATHOVER).NE.0) NRUNS=NRUNS+1
PRINT '(A,G20.10)','mcpath2> Number of MC runs required                                       ',NRUNS
CLOSE(LUNIT)
ALLOCATE(EOFS(NPATH),PATHFRAMES(NPATH,3*NATOMS),NORM1(NPATH))
ALLOCATE(PROBS(NPATH),PROBQ(MCPATHBINS),PROBSQ(NPATH,MCPATHBINS),PROBSPAIR(NPATH))
ALLOCATE(ALLZEROS(NPATH),ALLZEROQ(MCPATHBINS),ALLZEROSTOT(NPATH))
ALLOCATE(SMOOTHS(NPATH),SMOOTHSPUP(NPATH),SMOOTHSPDOWN(NPATH),SMOOTHSPFIXUP(NPATH),SMOOTHSPFIXDOWN(NPATH), &
  &      SMOOTHQ(MCPATHBINS),SMOOTHSQ(NPATH,MCPATHBINS))
ALLOCATE(SLENGTH(NPATH),QFRAME(NPATH),QFRAMEAV(NPATH),QFRAMESD(NPATH))
LUNIT=GETUNIT()
OPEN(LUNIT,FILE='path.xyz',STATUS='OLD')
DO J1=1,NPATH
   READ(LUNIT,*) SDUMMY
   READ(LUNIT,*) SDUMMY,EOFS(J1)
   DO J2=1,NATOMS
      READ(LUNIT,'(A5,1X,3F20.10)') SSYM,PATHFRAMES(J1,3*(J2-1)+1:3*(J2-1)+3)
   ENDDO
ENDDO
CLOSE(LUNIT)
!
! Just in case things aren't permutationally aligned and should be, let's check this here.
!
DUMMY2=-1.0D0
SLENGTH(1)=0.0D0
DO J1=2,NPATH
   NTRIES=0
863 CONTINUE
   NTRIES=NTRIES+1
   CALL ALIGN_DECIDE(PATHFRAMES(J1-1,1:3*NATOMS),PATHFRAMES(J1,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
   IF (DEBUG) PRINT '(A,I6,A,I6,A,G20.10)','mcpath2> Distance between frames ',J1,' and ',J1-1,' after alignment is ',DUMMY
   IF (DEBUG) PRINT '(3(A,I6),A,G20.10)','mcpath2> Distance between frames ',J1,' and ',J1-1, &
  &                 ' after alignment attempt ',NTRIES,' is ',DUMMY
   IF ((DUMMY.GT.2.0D0*MAXBFGS).AND.(NTRIES.LT.10)) GOTO 863
   PRINT '(3(A,I6),A,G20.10)','mcpath2> Distance between frames ',J1,' and ',J1-1,' after alignment attempt ',NTRIES,' is ',DUMMY
   IF (DUMMY.GT.DUMMY2) DUMMY2=DUMMY
   IF (DUMMY.GT.1.0D0) THEN
      PRINT '(A,I6)','frame ',J1-1
      PRINT '(3G20.10)',PATHFRAMES(J1-1,1:3*NATOMS)
      PRINT '(A,I6)','frame ',J1
      PRINT '(3G20.10)',PATHFRAMES(J1,1:3*NATOMS)
      STOP
   ENDIF
   SLENGTH(J1)=SLENGTH(J1-1)+DUMMY
ENDDO
PRINT '(A,G20.10)','mcpath2> largest distance between frames in path.xyz is now',DUMMY2
!
! Should now be safe to turn permutations off!
!
PERMDISTSAVE=.FALSE. 
LPERMDISTSAVE=.FALSE.
PERMDIST=.FALSE. 
LPERMDIST=.FALSE.
DUMMY2=-1.0D0
DO J1=2,NPATH
   CALL ALIGN_DECIDE(PATHFRAMES(J1-1,1:3*NATOMS),PATHFRAMES(J1,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
   IF (DEBUG) PRINT '(A,I6,A,I6,A,G20.10)','mcpath2> Distance between frames ',J1,' and ',J1-1,' is now ',DUMMY
   IF (DUMMY.GT.DUMMY2) DUMMY2=DUMMY
ENDDO
PRINT '(A,G20.10)','mcpath2> check: largest distance between frames in path.xyz is ',DUMMY2

IF (DUMMY2.GT.1.0D0) STOP !!! debug DJW

NSTRUCTREF=NPATH

ALLOCATE(COORDSREF(NSTRUCTREF,3*NATOMS),VREF(NSTRUCTREF),NEAREST(NSTRUCTREF))
ALLOCATE(DSTRUCT(NSTRUCTREF))
ALLOCATE(QORDERHIST(MCPATHBINS),BINLABELQORDER(MCPATHBINS),QORDERVISITS(MCPATHBINS))
ALLOCATE(NEARESTF(NPATH),NEARESTFW(NPATH))

DO J1=1,NPATH
   VREF(J1)=EOFS(J1) 
   COORDSREF(J1,1:3*NATOMS)=PATHFRAMES(J1,1:3*NATOMS)
ENDDO

HIGHESTREF=1
LASTMIN=.TRUE.
LASTTS=.FALSE.
PRINT '(A,I10,A,G20.10)','mcpath2> minimum assumed             for frame ',1,' energy=',VREF(1)
!
! Assume we will always have an order parameter. Need to replace the call below appropriately.
!
IF (.TRUE.) THEN
   DO J1=1,NPATH
      DO K1=1,3
         DO K2=1,NATOMS
            QTEMP(K1,K2)=PATHFRAMES(J1,3*(K2-1)+K1)
         ENDDO
      ENDDO
      CALL GETORDER(QTEMP,QORDER,NATOMS)
      QFRAME(J1)=QORDER
   ENDDO

   DO K1=1,3
      DO K2=1,NATOMS
         QTEMP(K1,K2)=COORDSREF(1,3*(K2-1)+K1)
      ENDDO
   ENDDO
   CALL GETORDER(QTEMP,QORDER,NATOMS)
   PRINT '(A,G20.10)','mcpath2> Q=',QORDER
ENDIF

DO J1=2,NPATH-1
   IF (EOFS(J1-1)+EDIFFTOL*10.0D0<EOFS(J1).AND.EOFS(J1)>EOFS(J1+1)+EDIFFTOL*10.0D0) THEN
      IF (VREF(J1).GT.VREF(HIGHESTREF)) HIGHESTREF=J1
      PRINT '(A,I10,2(A,G20.10))','mcpath2> transition state identified for frame ',J1,' energy=', &
  &         VREF(J1)
      IF (.NOT.LASTMIN) THEN
         PRINT '(A)','mcpath2> ERROR *** previous stationary point identified was not a minimum'
         STOP
      ENDIF
      LASTMIN=.FALSE.
      LASTTS=.TRUE.
      IF (.TRUE.) THEN
         DO K1=1,3
            DO K2=1,NATOMS
               QTEMP(K1,K2)=COORDSREF(J1,3*(K2-1)+K1)
            ENDDO
         ENDDO
         CALL GETORDER(QTEMP,QORDER,NATOMS)
         PRINT '(A,G20.10)','mcpath2> Q=',QORDER
      ENDIF
   ENDIF
   IF (EOFS(J1-1)>EOFS(J1).AND.EOFS(J1)<=EOFS(J1+1)) THEN
      IF (VREF(J1).GT.VREF(HIGHESTREF)) HIGHESTREF=J1
      CALL ALIGN_DECIDE(COORDSREF(J1-1,1:3*NATOMS),COORDSREF(J1,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
      PRINT '(A,I10,2(A,G20.10))','mcpath2> minimum identified for frame ',J1,' energy=', &
  &         VREF(J1)
      IF (.NOT.LASTTS) THEN
         PRINT '(A)','mcpath2> ERROR *** previous stationary point identified was not a transition state'
         STOP
      ENDIF
      LASTTS=.FALSE.
      LASTMIN=.TRUE.
      IF (.TRUE.) THEN
         DO K1=1,3
            DO K2=1,NATOMS
               QTEMP(K1,K2)=COORDSREF(J1,3*(K2-1)+K1)
            ENDDO
         ENDDO
         CALL GETORDER(QTEMP,QORDER,NATOMS)
         PRINT '(A,G20.10)','mcpath2> Q=',QORDER
      ENDIF
   ENDIF
ENDDO
CALL ALIGN_DECIDE(COORDSREF(NPATH-1,1:3*NATOMS),COORDSREF(NPATH,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
      PRINT '(A,I10,2(A,G20.10))','mcpath2> minimum assumed for frame ',NPATH,' energy=', &
  &         VREF(NPATH)
!
! This interval will depend upon the order parameter. 
!
QORDERHISTINT=(MAX(MCPATHQMAX,MCPATHQMIN)-MIN(MCPATHQMAX,MCPATHQMIN))/MCPATHBINS 
DO J1=1,MCPATHBINS
! these Q values correspond to the bottom of the Q bins
   BINLABELQORDER(J1)=MIN(MCPATHQMAX,MCPATHQMIN)+QORDERHISTINT*(J1-1.0D0) 
ENDDO

REDOT=.FALSE.
INQUIRE(FILE='redo',EXIST=YESNO)

IF (MCPATHEQUIL.EQ.0) THEN
   PRINT '(A)','mcpath2> zero equilibration steps specified - reading previous results from disc'
   GOTO 963
ENDIF
IF (YESNO) THEN
   PRINT '(A)','mcpath2> File redo detected in working directory - reading previous results from disc'
   REDOT=.TRUE.
   GOTO 963
ENDIF

PRINT '(4(A,I10))','mcpath2> Total number of reference structures=',NSTRUCTREF
CALL FLUSH(6)
!
! Align path frames if necessary. Should not be executed, since PERMDISTSAVE and LPERMDISTSAVE
! were turned off above!
!
IF (PERMDISTSAVE.OR.LPERMDISTSAVE) THEN
   CALL ALIGN_DECIDE(COORDSREF(1,1:3*NATOMS),PATHFRAMES(1,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
   PRINT '(A,F20.10)','mcpath2> first distance=',DUMMY
   DO J1=2,NPATH
      CALL ALIGN_DECIDE(PATHFRAMES(J1-1,1:3*NATOMS),PATHFRAMES(J1,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
      PRINT '(A,I10,F20.10)','mcpath2> frame and distance=',J1,DUMMY
   ENDDO
!
! Need to resave COORDSREF to get consistent permutational alignment.
!
   DO J1=1,NSTRUCTREF 
      COORDSREF(J1,1:3*NATOMS)=PATHFRAMES(J1,1:3*NATOMS)
   ENDDO
   PERMDIST=.FALSE.
   LPERMDIST=.FALSE. ! freeze permutations of references - should already be optimal now
ENDIF
IF (.NOT.MCBIAST) GOTO 987

PRINT '(A,I10,A)','mcpath2> ',NSTRUCTREF,' reference structures:'
DO J1=1,NSTRUCTREF
   PRINT '(A,I10,A,G20.10)','frame ',J1,' energy=',VREF(J1)
ENDDO
!
! The number of extra reference structures is now fixed.
!
987 CONTINUE ! jump here if no bias
DOBLOCK=1
IF (MCPATHDOBLOCK.GT.0) DOBLOCK=MCPATHDOBLOCK
IF (DOBLOCK.GT.NRUNS) THEN
   PRINT '(A)','mcpath2> ERROR *** specified block > number of runs'
   PRINT '(A)','mcpath2> DOBLOCK,MCPATHDOBLOCK,NRUNS=',DOBLOCK,MCPATHDOBLOCK,NRUNS
ENDIF

IF (PBST.OR.(NCPU.GT.0)) THEN
   PRINT '(A)','mcpath2> Copying original odata file to odata.save'
   CALL MYSYSTEM(STATUS,LDEBUG,'cp odata odata.save')
   IF (PBST) THEN
      INQUIRE(FILE='nodes.info',EXIST=YESNO)
      IF (.NOT.YESNO) THEN
         PRINT '(A)','mcpath2> No nodes.info file - stop'
         STOP
      ENDIF
      LUNIT=GETUNIT()
      OPEN(UNIT=LUNIT,FILE='nodes.info',STATUS='OLD')
      READ(LUNIT,*) NNODES
      IF (ALLOCATED(NODENAME)) DEALLOCATE(NODENAME) 
      ALLOCATE(NODENAME(NNODES),PID(NNODES)) 
      WRITE(*,'(A,I2,a)') 'mcpath2> Following ',NNODES,' nodes are available:'
      DO I=1,NNODES
         READ(LUNIT,'(A)') ARG
         PRINT '(A)', TRIM(ADJUSTL(ARG))
         NODENAME(I)=TRIM(ADJUSTL(ARG))
      ENDDO
      READ(LUNIT,'(A)') USERNAME
      READ(LUNIT,'(A)') WORKINGDIRECTORY
      CLOSE(LUNIT)
      WRITE(*,'(2A)') 'mcpath2> Working in directory ',TRIM(ADJUSTL(WORKINGDIRECTORY))
      CALL FLUSH(6)
   ELSE
      PRINT '(A)','mcpath2> Interactive run - not checking for nodes.info file'
      ALLOCATE(NODENAME(NCPU),PID(NCPU)) 
      NNODES=NCPU
   ENDIF
   IF (NNODES.GT.NRUNS) NNODES=NRUNS
!
! Initial job submission - one OPTIM job per core.
!
   NRUNNING=0
   DO J3=1,NNODES
      IF (SLEEPTIME1.GT.0.0D0) CALL SYSTEM(TRIM(ADJUSTL(SLEEPSTRING1)))
      CALL FLUSH(6)
      CALL FORK_SUBR(PID(J3)) ! PID is zero in the child, non-zero in the parent
      IF ((PID(J3).NE.0)) WRITE(*,'(A,I8)') 'mcpath2> forked connect run process id=',PID(J3)
      CALL FLUSH(6)
!     PRINT *,'J3,DOBLOCK,PID=',J3,DOBLOCK,PID(J3)
      IF (PID(J3).EQ.0) CALL SUBMITOPTIMJOB(J3,DEBUG,DOBLOCK,'OPTIM.',NNODES,NODENAME,USERNAME,WORKINGDIRECTORY)
      NRUNNING=NRUNNING+1
      ENDSAMPLE=DOBLOCK*MCPATHBLOCK-(DOBLOCK-1)*MCPATHOVER
      ENDSAMPLE=MIN(ENDSAMPLE,NPATH)
!     PRINT *,'NRUNNING,ENDSAMPLE=',NRUNNING,ENDSAMPLE
!     PRINT *,'initial loop DOBLOCK=',DOBLOCK
      IF (DOBLOCK*MCPATHBLOCK-(DOBLOCK-1)*MCPATHOVER.LE.NPATH) THEN
         DOBLOCK=DOBLOCK+1
      ELSE
         IF (DOBLOCK.GT.NRUNS) PID(J3:NNODES)=-1
         DOBLOCK=DOBLOCK+1
         EXIT
      ENDIF
   ENDDO
   DO WHILE (NRUNNING.GT.0) 
      IF (SLEEPTIME1.GT.0.0D0) CALL SYSTEM(TRIM(ADJUSTL(SLEEPSTRING1))) ! to allow jobs to catch up for demos
      KILLED=.FALSE.
      CALL FLUSH(6) ! the child process may duplicate output without this line!
  113 CONTINUE
      CALL FLUSH(6) ! the child process may duplicate output without this line!
      CALL WAIT_SUBR(PIDDONE,STATUS)
      NRUNNING=NRUNNING-1
   11 CONTINUE
      IF (DEBUG) THEN
         PRINT '(A,2I8)','mcpath2> PIDDONE,STATUS,PID=',PIDDONE,STATUS
         PRINT '(I10)', PID(1:NNODES)
      ENDIF
      CALL FLUSH(6) ! the child process may duplicate output without this line!
      IF (PIDDONE.GT.0) THEN
         IF (DEBUG) PRINT '(A,I8,A,I6)','cycle2> PID ',PIDDONE,' has finished with exit status ',STATUS
         DO J2=1,NNODES
!           PRINT *,'J2,PID(J2),PIDDONE=',J2,PID(J2),PIDDONE
            IF (PIDDONE.EQ.PID(J2)) THEN
               IF (STATUS.NE.0) KILLED=.TRUE. ! Incomplete OPTIM jobs would ideally return a non-zero exit code
               NEWJOB=J2
               IF (DEBUG) PRINT '(2(A,I8))','cycle2> PID ',PIDDONE,' has finished on cpu ',J2
               GOTO 10
            ENDIF
         ENDDO
         PRINT*,'ERROR - PID of completed child process not recognised: ',PIDDONE
         STOP
      ELSE
112      CALL FLUSH(6) ! the child process may duplicate output without this line!
         PRINT '(A,I20)','mcpath2> WARNING - WAIT returned system error code ',-PIDDONE!
!
! Try calling wait again to see if this fixes things.
! For very short OPTIM jobs WAIT may have trouble keeping up!
! 
         CALL MYSYSTEM(STATUS,DEBUG,' sleep 1')
114      CONTINUE
         CALL FLUSH(6) ! the child process may duplicate output without this line!
         CALL WAIT_SUBR(PIDDONE,STATUS)
         PRINT '(2(A,I8))','cycle2> on calling wait again pid=',PIDDONE,' status=',STATUS
         IF (PIDDONE.GT.0) GOTO 11
PRINT *,'mcpath2> here D'
!        STOP 
      ENDIF
10    CONTINUE
!  
!  Identify OPTIM jobs that did not terminate with exit code 0
!  In such cases KILLED should be .TRUE.
!  
      CALL FLUSH(6) ! the child process may duplicate output without this line!
      WRITE(*,'(3(A,I8))') 'mcpath2> analysing result of OPTIM run on CPU ',NEWJOB,' for process id ',PID(NEWJOB)
      WRITE(CONNSTR,'(I10)') PID(NEWJOB)
      IF (KILLED) THEN
         WRITE(*,'(3(A,I8))') 'mcpath2> OPTIM run on CPU ',NEWJOB,' was unsuccessful'
         STOP
      ENDIF

      ENDSAMPLE=DOBLOCK*MCPATHBLOCK-(DOBLOCK-1)*MCPATHOVER
      ENDSAMPLE=MIN(ENDSAMPLE,NPATH)
!     PRINT *,'DOBLOCK,NRUNNING,ENDSAMPLE=',DOBLOCK,NRUNNING,ENDSAMPLE
      IF (DOBLOCK.LE.NRUNS) THEN
         CALL FLUSH(6)
         CALL FORK_SUBR(PID(NEWJOB))
         IF (PID(NEWJOB).EQ.0) CALL SUBMITOPTIMJOB(NEWJOB,DEBUG,DOBLOCK,'OPTIM.',NNODES,NODENAME,USERNAME,WORKINGDIRECTORY)
         NRUNNING=NRUNNING+1
         DOBLOCK=DOBLOCK+1
!        PRINT *,'new job submitted on CPU ',NEWJOB,' NRUNNING,ENDSAMPLE=',NRUNNING,ENDSAMPLE
      ELSE
         DOBLOCK=DOBLOCK+1
         PID(NEWJOB)=-1
      ENDIF
   ENDDO
   DEALLOCATE(NODENAME,PID)
   PRINT '(A)','mcpath2> Restoring original odata file from odata.save'
   CALL MYSYSTEM(STATUS,LDEBUG,'cp odata.save odata')
ELSE

   753 CONTINUE ! jump back here if we are doing multiple blocks

   STARTSAMPLE=(DOBLOCK-1)*(MCPATHBLOCK-MCPATHOVER)+1
   ENDSAMPLE=DOBLOCK*MCPATHBLOCK-(DOBLOCK-1)*MCPATHOVER
   ENDSAMPLE=MIN(ENDSAMPLE,NPATH)
   PRINT '(A,I6)','mcpath2> Sampling block ',DOBLOCK
   PRINT '(A,I6,A,I6)','mcpath2> Sampling block defined by frames ',STARTSAMPLE,' and ',ENDSAMPLE

   CALL SAMPLEBLOCK(ENDSAMPLE,STARTSAMPLE,COORDSREF,VREF,NSTRUCTREF,DOBLOCK,NPATH,QFRAMEAV,QFRAMESD,QORDERHIST, &
  &              EOFS,QORDERHISTINT,BINLABELQORDER,SLENGTH,QFRAME,PROBSFAKE)

   PRINT *,'DOBLOCK=',DOBLOCK
   IF ((ENDSAMPLE.LT.NPATH).AND.(MCPATHDOBLOCK.EQ.0)) THEN
      DOBLOCK=DOBLOCK+1
!     DEALLOCATE(LCOORDSREF,LVREF,LEOFS)
      GOTO 753
   ENDIF
ENDIF

DEALLOCATE(COORDSREF,VREF,NEAREST,DSTRUCT)
DEALLOCATE(NEARESTF,NEARESTFW)
DEALLOCATE(QORDERHIST,QORDERVISITS)
DEALLOCATE(PATHFRAMES)

963 CONTINUE

NAPPEND=1
DO WHILE (NAPPEND.LE.2048) ! start of do while over intermediate results

   PATHSTR='f.block'
   IF (NAPPEND.EQ.2048) THEN
      APPSTRING=''
   ELSE
      WRITE(APPSTRING,'(I4)') NAPPEND
      APPSTRING='.' // TRIM(ADJUSTL(APPSTRING))
   ENDIF
! 
!  Assume that data will exist for all the blocks if it is there for block1.
!
   FNAME=TRIM(ADJUSTL('mcpath.f.block1'//TRIM(ADJUSTL(APPSTRING))))
   INQUIRE(FILE=FNAME,EXIST=YESNO)
   IF (YESNO) THEN
      PRINT '(A)','mcpath2> Calculating distributions for mcpath files with extension '//TRIM(ADJUSTL(APPSTRING))
   ELSE
      NAPPEND=NAPPEND*2
      CYCLE
   ENDIF
   IF (.TRUE.) THEN ! recover QFRAMEAV and QFRAMESD from the mcpath.f files. Blocks may overlap. 
                    ! Use the entry with the most visits.
                    ! These values are <Q> and sqrt(<Q^2>-<Q)^2) for frame s in block=window=replica r.
      ALLOCATE(QFBEST(NPATH))
      QFBEST(1:NPATH)=-1.0D0
      QFRAMEAV(1:NPATH)=0.0D0
      QFRAMESD(1:NPATH)=0.0D0

      DO J1=1,NRUNS
         WRITE(REPSTRING,'(I6)') J1
         WRITE(FSTRING,'(A)') 'mcpath.' // TRIM(ADJUSTL(PATHSTR)) // TRIM(ADJUSTL(REPSTRING)) // TRIM(ADJUSTL(APPSTRING))
         INQUIRE(FILE=TRIM(ADJUSTL(FSTRING)),EXIST=YESNO)
         IF (.NOT.YESNO) THEN
            PRINT '(3A)','File ',TRIM(ADJUSTL(FSTRING)),' not found'
            STOP
         ENDIF
         OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FSTRING)),STATUS='OLD')
         STARTSAMPLE=(J1-1)*(MCPATHBLOCK-MCPATHOVER)+1
         ENDSAMPLE=J1*MCPATHBLOCK-(J1-1)*MCPATHOVER
         ENDSAMPLE=MIN(ENDSAMPLE,NPATH)
!        PRINT *,'J1,STARTSAMPLE,ENDSAMPLE=',STARTSAMPLE,ENDSAMPLE
         DO J2=1,STARTSAMPLE-1
            READ(1,*) NDUMMY
         ENDDO
         DO J2=STARTSAMPLE,ENDSAMPLE
            READ(1,*) NDUMMY,DUMMY1,DUMMY2,DUMMY3,DUMMY4,DUMMY5,DUMMY,DUMMY6,DUMMY7
!
! DUMMY is the % visits to this s bin in the run. Use the block with the
! most visits for Q mean and SD.
!
            IF(DUMMY.GT.QFBEST(J2)) QFRAMEAV(J2)=DUMMY6
            IF(DUMMY.GT.QFBEST(J2)) QFRAMESD(J2)=DUMMY7 
         ENDDO
         CLOSE(1)
      ENDDO
      PRINT '(A)','mcpath2> Q frame mean and standard deviation read from mcpath.f files.'
      PRINT '(A)','mcpath2> Q frame mean, frame value, ratio, and standard deviation:'
      PRINT '(I6,4G20.10)',(J1,QFRAMEAV(J1),QFRAME(J1),QFRAMEAV(J1)/MAX(1.0D-100,QFRAME(J1)),QFRAMESD(J1),J1=1,NPATH)
      DEALLOCATE(QFBEST)
   ENDIF
!
! Call WHAM for successive windows.
! (1) Starting from the bottom: pair up
! (2) Simultaneous optimisation of bin weights for normalisations fixed from pair up.
! (3) Starting from the top: pair down
! (4) Simultaneous optimisation of bin weights for normalisations fixed from pair down.
! (5) Unconstrained optimisation starting from (4).
!
! Items (1) to (4) have been commented and moved to the
! bottom of this file.
!
   GWIDTHS=MCPATHGWS
   PRINT '(A,G20.10)','mcpath2> standard deviation sigma for P(s) bins smoothing is ',SQRT(GWIDTHS)
   IF (ALLOCATED(ZNORM)) DEALLOCATE(ZNORM)
   IF (ALLOCATED(ZNORMPAIR)) DEALLOCATE(ZNORMPAIR)
   ALLOCATE(ZNORM(NRUNS),ZNORMPAIR(NRUNS))
!
! First pair up. Not any more - all pair code has been commented and moved to the
! bottom of this file.
!
   DO J1=1,NPATH
      IF (LINVART) THEN ! positive probabilities!
         PROBS(J1)=(DPRAND()+1.0D0)/2.0D0
      ELSE
         PROBS(J1)=DPRAND()
!        PROBS(J1)=-100.0D0
      ENDIF
   ENDDO
   DO J1=1,NRUNS
      ZNORM(J1)=DPRAND()
   ENDDO
   IF (LINVART) ZNORM(1)=1.0D0
   IF (LINVART) PROBS(1)=1.0D0
   IF (.NOT.LINVART) ZNORM(1)=0.0D0
   IF (.NOT.LINVART) PROBS(1)=0.0D0

   FIXZ=.FALSE.
   ALLZEROSTOT(1:NPATH)=.TRUE.
   CALL MYWHAM(NRUNS,NPATH,PATHSTR,PROBS,APPSTRING,LINVART,OFAILS,1,NRUNS,ZNORM,FIXZ,ALLZEROSTOT)
   IF (OFAILS) THEN
      PRINT '(A)','mcpath> Pofs and PofsQ will not be produced'
      GOTO 864
   ENDIF
   IF (DEBUG) PRINT '(A,2I6,A,2G20.10)','after all windows mywham PROBS for unconstrained:'
   IF (DEBUG) PRINT '(G20.10)',PROBS(1:NPATH)
   IF (DEBUG) PRINT '(A,2I6,A,2G20.10)','after all windows mywham znorm vector for unconstrained:'
   IF (DEBUG) PRINT '(G20.10)',ZNORM(1:NRUNS)

   IF (.NOT.LINVART) THEN
      DUMMY=PROBS(1)
      DO J1=1,NPATH
         PROBS(J1)=EXP(PROBS(J1)-DUMMY)
         IF (ALLZEROSTOT(J1)) PRINT '(A,I6,A,G20.10)','mcpath2> s bin ',J1,' was not visited, will zero probability of ',PROBS(J1)
         IF (ALLZEROSTOT(J1)) PROBS(J1)=0.0D0
      ENDDO
   ENDIF

   DUMMY=-1.0D100
   IF (DEBUG) PRINT '(A)','mcpath2> frames, s, ln(P(s)),P(s) from WHAM unconstrained:'
   DO J1=1,NPATH
      IF (DEBUG) PRINT '(I6,3G20.10)',J1,EOFS(J1),LOG(MAX(1.0D-200,PROBS(J1))),PROBS(J1)
      IF (PROBS(J1).GT.DUMMY) DUMMY=PROBS(J1)
   ENDDO
!
   CALL MAKESMOOTH(NPATH,PROBS,GWIDTHS,SMOOTHS,SLENGTH)

   LUNIT=GETUNIT()
   IF (LINVART) THEN
      FNAME='Pofs.pair.lin'//TRIM(ADJUSTL(APPSTRING))
   ELSE
      FNAME='Pofs.pair.ln'//TRIM(ADJUSTL(APPSTRING))
   ENDIF
   OPEN(LUNIT,FILE=FNAME,STATUS='UNKNOWN')
   DUMMY=0.0D0
   DO J1=1,NPATH
      WRITE(LUNIT,'(I6,9G20.10)') J1,SLENGTH(J1),SMOOTHSPUP(J1),SMOOTHSPDOWN(J1),SMOOTHSPFIXUP(J1),SMOOTHSPFIXDOWN(J1), &
   &                              SMOOTHS(J1),PROBS(J1),QFRAME(J1),EOFS(J1)
      DUMMY=DUMMY+SMOOTHS(J1)
   ENDDO
   CLOSE(LUNIT)
   PRINT '(A,G20.10)','mcpath2> P(s) normalisation=',DUMMY

   LUNIT=GETUNIT()
   IF (LINVART) THEN
      FNAME='Pofs.lin'//TRIM(ADJUSTL(APPSTRING))
   ELSE
      FNAME='Pofs.ln'//TRIM(ADJUSTL(APPSTRING))
   ENDIF
   OPEN(LUNIT,FILE=FNAME,STATUS='UNKNOWN')
   DO J1=1,NPATH
      IF (ABS(SMOOTHS(J1)).LT.1.0D-90) SMOOTHS(J1)=0.0D0 ! to make sure format is OK for gnuplot
      WRITE(LUNIT,'(I6,5G20.10)') J1,SLENGTH(J1),SMOOTHS(J1),PROBS(J1),QFRAME(J1),EOFS(J1)
   ENDDO
   CLOSE(LUNIT)

   PROBQFROMS(1:MCPATHBINS)=0.0D0
   DO J1=1,NPATH
      IBININDEX=INT((QFRAME(J1)-MIN(MCPATHQMAX,MCPATHQMIN))/QORDERHISTINT)+1
      PROBQFROMS(IBININDEX)=PROBQFROMS(IBININDEX)+QFRAME(J1)*PROBS(J1)
   ENDDO
   PROBQFROMS2(1:MCPATHBINS)=0.0D0
 
   DO J1=1,NPATH ! commented DJW
      IBININDEX=INT((QFRAMEAV(J1)-MIN(MCPATHQMAX,MCPATHQMIN))/QORDERHISTINT)+1 ! uncommented DJW
      PROBQFROMS2(IBININDEX)=PROBQFROMS2(IBININDEX)+QFRAMEAV(J1)*PROBS(J1) ! uncommented DJW
   ENDDO ! uncommented DJW
!
! We have a standard deviation for Q in each s bin, so why not use this to
! smooth the corresponding Q(s)?
!
!    DO J1=1,NPATH
!!       
! ! Need to normalise over accessible Q bins ! commented DJW
! !
!       QFRAMESD(J1)=MAX(QFRAMESD(J1),1.0D-100) ! otherwise one day we will divide by zero. Trust me. DJW
! !     QFRAMESD(J1)=SQRT(0.001D0) ! debug DJW
!       DUMMY=0.0D0
!       DO J2=1,MCPATHBINS
!          DUMMY=DUMMY+EXP(-(BINLABELQORDER(J2)-QFRAMEAV(J1))**2/(2.0D0*QFRAMESD(J1)**2))
!       ENDDO
!       DO J2=1,MCPATHBINS
!          PROBQFROMS2(J2)=PROBQFROMS2(J2)+DUMMY*EXP(-(BINLABELQORDER(J2)-QFRAMEAV(J1))**2/(2.0D0*QFRAMESD(J1)**2))*PROBS(J1)
!       ENDDO
!    ENDDO

   PROBQFROMSPAIR(1:MCPATHBINS)=PROBQFROMS(1:MCPATHBINS)
   PROBQFROMSPAIR2(1:MCPATHBINS)=PROBQFROMS2(1:MCPATHBINS)

   864 CONTINUE ! jump here if P(s) overlap failed and OFAILS is true

   PATHSTR='Q.block'
   IF (ALLOCATED(ZNORM)) DEALLOCATE(ZNORM)
   IF (ALLOCATED(ZNORMPAIR)) DEALLOCATE(ZNORMPAIR)
   ALLOCATE(ZNORM(NRUNS))

   DO J1=1,MCPATHBINS
      IF (LINVART) THEN ! positive probabilities!
         PROBQ(J1)=(DPRAND()+1.0D0)/2.0D0
      ELSE 
         PROBQ(J1)=DPRAND()
      ENDIF
   ENDDO
   DO J1=1,NRUNS
      ZNORM(J1)=DPRAND()
   ENDDO
!  IF (LINVART) ZNORM(1)=1.0D0
!  IF (.NOT.LINVART) ZNORM(1)=0.0D0

   FIXZ=.FALSE.
   ALLZEROQ(1:MCPATHBINS)=.TRUE.
   CALL MYWHAM(NRUNS,MCPATHBINS,PATHSTR,PROBQ,APPSTRING,LINVART,OFAILQ,1,NRUNS,ZNORM,FIXZ,ALLZEROQ)
   IF (OFAILQ) THEN
      PRINT '(A)','mcpath> PofQ and PofsQ will not be produced'
      GOTO 862
   ENDIF
   IF (.NOT.LINVART) THEN
      DUMMY=PROBQ(1)
      DO J1=1,MCPATHBINS
         PROBQ(J1)=EXP(PROBQ(J1)-DUMMY)
         IF (ALLZEROQ(J1)) PRINT '(A,I6,A,G20.10)','mcpath2> Q bin ',J1,' was not visited, will zero probability of ',PROBQ(J1)
         IF (ALLZEROQ(J1)) PROBQ(J1)=0.0D0
      ENDDO
   ENDIF

   DUMMY=-1.0D100
   DO J1=1,MCPATHBINS
      IF (PROBQ(J1).GT.DUMMY) DUMMY=PROBQ(J1)
   ENDDO
   DO J1=1,MCPATHBINS
      PROBQ(J1)=PROBQ(J1)/DUMMY
   ENDDO
   GWIDTHQ=MCPATHGWQ
   SMOOTHQ(1:MCPATHBINS)=0.0D0
   SMOOTHQFROMS(1:MCPATHBINS)=0.0D0
   SMOOTHQFROMSPAIR(1:MCPATHBINS)=0.0D0
   SMOOTHQFROMS2(1:MCPATHBINS)=0.0D0
   SMOOTHQFROMSPAIR2(1:MCPATHBINS)=0.0D0
   DUMMY=0.0D0
   DO J1=1,MCPATHBINS
      DUMMY=DUMMY+PROBQ(J1)
   ENDDO
   PROBQ(1:MCPATHBINS)=PROBQ(1:MCPATHBINS)/DUMMY
   DUMMYQ=1.0D0/(SQRT(GWIDTHQ)*2.50662827463100D0)
!
! Normalise each contribution over accessible bins to conserve probabilty at the end points
!
   DO J2=1,MCPATHBINS
      NORM2(J2)=0.0D0
      DO J1=1,MCPATHBINS
         NORM2(J2)=NORM2(J2)+EXP(-(BINLABELQORDER(J1)-BINLABELQORDER(J2))**2/(2.0D0*GWIDTHQ))
      ENDDO
   ENDDO
   DO J1=1,MCPATHBINS
      DO J2=1,MCPATHBINS
         SMOOTHQ(J1)=SMOOTHQ(J1)+PROBQ(J2)*EXP(-((BINLABELQORDER(J1)-BINLABELQORDER(J2))**2/(2.0D0*GWIDTHQ)))/NORM2(J2)
      ENDDO
   ENDDO

   IF (OFAILS) THEN
      PROBQFROMS(1:MCPATHBINS)=0.0D0
      SMOOTHQFROMS(1:MCPATHBINS)=0.0D0
      PROBQFROMS2(1:MCPATHBINS)=0.0D0
      SMOOTHQFROMS2(1:MCPATHBINS)=0.0D0
   ELSE
      DUMMY=0.0D0
      DO J1=1,MCPATHBINS
         DUMMY=DUMMY+PROBQFROMS(J1)
      ENDDO
      PROBQFROMS(1:MCPATHBINS)=PROBQFROMS(1:MCPATHBINS)/DUMMY
      DUMMY=0.0D0
      DO J1=1,MCPATHBINS
         DUMMY=DUMMY+PROBQFROMSPAIR(J1)
      ENDDO
      PROBQFROMSPAIR(1:MCPATHBINS)=PROBQFROMSPAIR(1:MCPATHBINS)/DUMMY
      DO J1=1,MCPATHBINS
         DO J2=1,MCPATHBINS
            SMOOTHQFROMS(J1)=SMOOTHQFROMS(J1)+ &
  &                           PROBQFROMS(J2)*EXP(-((BINLABELQORDER(J1)-BINLABELQORDER(J2))**2/(2.0D0*GWIDTHQ)))/NORM2(J2)
            SMOOTHQFROMSPAIR(J1)=SMOOTHQFROMSPAIR(J1)+ &
  &                           PROBQFROMSPAIR(J2)*EXP(-((BINLABELQORDER(J1)-BINLABELQORDER(J2))**2/(2.0D0*GWIDTHQ)))/NORM2(J2)
         ENDDO
      ENDDO

      DUMMY=0.0D0
      DO J1=1,MCPATHBINS
         DUMMY=DUMMY+PROBQFROMS2(J1)
      ENDDO
      PROBQFROMS2(1:MCPATHBINS)=PROBQFROMS2(1:MCPATHBINS)/DUMMY
      DUMMY=0.0D0
      DO J1=1,MCPATHBINS
         DUMMY=DUMMY+PROBQFROMSPAIR2(J1)
      ENDDO
      PROBQFROMSPAIR2(1:MCPATHBINS)=PROBQFROMSPAIR2(1:MCPATHBINS)/DUMMY ! uncommented DJW
      DO J1=1,MCPATHBINS
!        SMOOTHQFROMS2(J1)=PROBQFROMS2(J1) ! because hopefully we have already used the SD for the s bin to smooth. commented DJW
         DO J2=1,MCPATHBINS
            SMOOTHQFROMS2(J1)=SMOOTHQFROMS2(J1)+ & ! uncommented DJW
  &                           PROBQFROMS2(J2)*EXP(-((BINLABELQORDER(J1)-BINLABELQORDER(J2))**2/(2.0D0*GWIDTHQ)))/NORM2(J2) ! uncommented DJW
            SMOOTHQFROMSPAIR2(J1)=SMOOTHQFROMSPAIR2(J1)+ &
  &                           PROBQFROMSPAIR2(J2)*EXP(-((BINLABELQORDER(J1)-BINLABELQORDER(J2))**2/(2.0D0*GWIDTHQ)))/NORM2(J2)
         ENDDO
      ENDDO
   ENDIF

   LUNIT=GETUNIT()
   IF (LINVART) THEN
      FNAME='PofQ.lin'//TRIM(ADJUSTL(APPSTRING))
   ELSE
      FNAME='PofQ.ln'//TRIM(ADJUSTL(APPSTRING))
   ENDIF
   OPEN(LUNIT,FILE=FNAME,STATUS='UNKNOWN')
   DUMMY=0.0D0
   DUMMY2=0.0D0
   DUMMY5=0.0D0
   DO J1=1,MCPATHBINS
      IF (ABS(SMOOTHQ(J1)).LT.1.0D-90) SMOOTHQ(J1)=0.0D0 ! to make gnuplot happy
      IF (ABS(SMOOTHQFROMS(J1)).LT.1.0D-90) SMOOTHQFROMS(J1)=0.0D0
      IF (ABS(SMOOTHQFROMS2(J1)).LT.1.0D-90) SMOOTHQFROMS2(J1)=0.0D0
      WRITE(LUNIT,'(I6,7G20.10)') J1,BINLABELQORDER(J1),SMOOTHQ(J1), &
  &                            SMOOTHQFROMS(J1),SMOOTHQFROMS2(J1),PROBQ(J1),PROBQFROMS(J1),PROBQFROMS2(J1)
      DUMMY=DUMMY+SMOOTHQ(J1)
      DUMMY2=DUMMY2+SMOOTHQFROMS(J1)
      DUMMY5=DUMMY5+SMOOTHQFROMS2(J1)
   ENDDO
   PRINT '(A,3G20.10)','mcpath2> P(Q) and P(Q) from P(s) (twice) normalisation=',DUMMY,DUMMY2,DUMMY5
   IF (NRUNS.GT.1) PRINT '(A,3G20.10)','mcpath2> WARNING *** for multiple replicas P(Q) from direct visits may contain systematic errors!'
   CLOSE(LUNIT)

   DO J1=1,NPATH
      NORM1(J1)=0.0D0
      DO J2=1,NPATH
         DO J3=1,MCPATHBINS
            IBININDEX=INT((QFRAME(J1)-MIN(MCPATHQMAX,MCPATHQMIN))/QORDERHISTINT)+1
            NORM1(J1)=NORM1(J1)+EXP(-(SLENGTH(J1)-SLENGTH(J2))**2/(2.0D0*GWIDTHS)) &
  &                                 *EXP(-((BINLABELQORDER(IBININDEX)-BINLABELQORDER(J3))**2/(2.0D0*GWIDTHQ)))
         ENDDO
      ENDDO
   ENDDO

   ALLOCATE(SMOOTH2D2(NPATH,MCPATHBINS),SMOOTH2D(NPATH,MCPATHBINS))
   SMOOTH2D(1:NPATH,1:MCPATHBINS)=0.0D0
   SMOOTH2D2(1:NPATH,1:MCPATHBINS)=0.0D0
   IF (.NOT.OFAILS) THEN
      DO J2=1,NPATH
         DO J3=1,MCPATHBINS
            DO J1=1,NPATH
               IBININDEX=INT((QFRAME(J1)-MIN(MCPATHQMAX,MCPATHQMIN))/QORDERHISTINT)+1
               SMOOTH2D(J2,J3)=SMOOTH2D(J2,J3)+PROBS(J1)*EXP(-(SLENGTH(J1)-SLENGTH(J2))**2/(2.0D0*GWIDTHS)) &
  &                                 *EXP(-((BINLABELQORDER(IBININDEX)-BINLABELQORDER(J3))**2/(2.0D0*GWIDTHQ)))/NORM1(J1)
               IBININDEX=INT((QFRAMEAV(J1)-MIN(MCPATHQMAX,MCPATHQMIN))/QORDERHISTINT)+1
               SMOOTH2D2(J2,J3)=SMOOTH2D2(J2,J3)+PROBS(J1)*EXP(-(SLENGTH(J1)-SLENGTH(J2))**2/(2.0D0*GWIDTHS)) &
  &                                 *EXP(-((BINLABELQORDER(IBININDEX)-BINLABELQORDER(J3))**2/(2.0D0*GWIDTHQ)))/NORM1(J1)
            ENDDO
         ENDDO
      ENDDO

      LUNIT=GETUNIT()
      IF (LINVART) THEN
         FNAME='PofsQ.lin'//TRIM(ADJUSTL(APPSTRING))
      ELSE
         FNAME='PofsQ.ln'//TRIM(ADJUSTL(APPSTRING))
      ENDIF
      OPEN(LUNIT,FILE=FNAME,STATUS='UNKNOWN')
      DUMMY=0.0D0
      DUMMY5=0.0D0
      DO J1=1,NPATH
         DO J2=1,MCPATHBINS
            WRITE(LUNIT,'(2I6,4G20.10)') J1,J2,SLENGTH(J1),BINLABELQORDER(J2),SMOOTH2D(J1,J2),SMOOTH2D2(J1,J2)
            DUMMY=DUMMY+SMOOTH2D(J1,J2)
            DUMMY5=DUMMY5+SMOOTH2D2(J1,J2)
         ENDDO
      ENDDO
      PRINT '(A,2G20.10)','mcpath2> P(s,Q) normalisations=',DUMMY,DUMMY5
      CLOSE(LUNIT)
   ENDIF

862 CONTINUE ! jump here if P(Q) overlap failed and OFAILQ is true.

   IF (ALLOCATED(SMOOTH2D)) DEALLOCATE(SMOOTH2D)
   IF (ALLOCATED(SMOOTH2D2)) DEALLOCATE(SMOOTH2D2)

   NAPPEND=2*NAPPEND

ENDDO ! end of do while over intermediate results

IF (LINVART) THEN
   LINVART=.FALSE.
   GOTO 963
ENDIF

DEALLOCATE(BINLABELQORDER,NORM1,QFRAME,QFRAMEAV,QFRAMESD)
IF (ALLOCATED(ZNORM)) DEALLOCATE(ZNORM)
IF (ALLOCATED(ZNORMPAIR)) DEALLOCATE(ZNORMPAIR)
IF (ALLOCATED(PROBSPAIR)) DEALLOCATE(PROBSPAIR)
DEALLOCATE(ALLZEROS,ALLZEROQ)

RETURN

END SUBROUTINE MCPATH2


SUBROUTINE SAMPLEBLOCK(ENDSAMPLE,STARTSAMPLE,COORDSREF,VREF,NSTRUCTREF,DOBLOCK,NPATH,QFRAMEAV,QFRAMESD,QORDERHIST, &
  &                    EOFS,QORDERHISTINT,BINLABELQORDER,SLENGTH,QFRAME,PROBSFAKE)
USE KEY
USE COMMONS
USE PORFUNCS
IMPLICIT NONE
INTEGER ENDSAMPLE,STARTSAMPLE,NSTRUCTREF,DOBLOCK,NPATH
DOUBLE PRECISION COORDSREF(NSTRUCTREF,3*NATOMS),VREF(NSTRUCTREF),PROBSFAKE(NPATH)
INTEGER LSTRUCTREF
INTEGER GETUNIT, IACCEPT
DOUBLE PRECISION, ALLOCATABLE :: LCOORDSREF(:,:), LVREF(:)
DOUBLE PRECISION, ALLOCATABLE :: LEOFS(:)
DOUBLE PRECISION DUMMY, GRAD(3*NATOMS), RMS
INTEGER K1, J1, NDUMMY, K2, LUNIT, J2
DOUBLE PRECISION QTEMP(3,NATOMS),QORDER,EOFS(NPATH)
DOUBLE PRECISION MCCOORDS(3*NATOMS), MCCOORDSO(3*NATOMS)
DOUBLE PRECISION SLENGTH(NPATH),QFRAME(NPATH),QFRAMEAV(NPATH),QFRAMESD(NPATH)
INTEGER NEAREST(NPATH), NEARESTF(NPATH), NTOT, IMAGEMINFO, K, NFAILS, IMAGEMIN, IMAGEMAX
INTEGER IREJDIST, IREJFRAME, IREJCHECK, NUPDATE, LNPATH, NCHECKINT, IMAGEMIN2, IMAGEMINF, IBININDEX
INTEGER IMAGEMINO, IMCSTEP, STATUS
DOUBLE PRECISION DIST2, RMAT(3,3), DPRAND, OPOTEL, IMAGEDIFF
DOUBLE PRECISION NEARESTFW(NPATH)
DOUBLE PRECISION QORDERHIST(MCPATHBINS),BINLABELQORDER(MCPATHBINS),QORDERVISITS(MCPATHBINS)
DOUBLE PRECISION DUMMY3, WOLD, VOLD, VNEW, WNEW, WCOMP
CHARACTER (LEN=8) SDUMMY
CHARACTER (LEN=80) FNAME, FNAME2
DOUBLE PRECISION WAC, DISTMIN, DISTMAX, VNEWSAVE, WNEWSAVE, XDISTMIN, XDISTMAX 
DOUBLE PRECISION QORDERHISTINT, DISTMINF, MEANBIAS, MEANBIASEQ, QORDERMEAN, QORDERSQMEAN
DOUBLE PRECISION QORDERSQMEAN2, QORDERMEAN2, QORDERO, RANDOM
LOGICAL RECOUNT
DOUBLE PRECISION X(NATOMS), Y(NATOMS), Z(NATOMS), XO(NATOMS), YO(NATOMS), ZO(NATOMS)

!
! NUPDATE specifies the interval for dynamically altering the maximum step size.
!
NUPDATE=100
PRINT '(A,I10,A)','mcpath2> Step size will be adjusted every ',NUPDATE ,' MC steps during equilibration phase'
LSTRUCTREF=ENDSAMPLE-STARTSAMPLE+1

ALLOCATE(LCOORDSREF(LSTRUCTREF,3*NATOMS),LVREF(LSTRUCTREF))
DO J1=1,LSTRUCTREF
   LCOORDSREF(J1,1:3*NATOMS)=COORDSREF(J1+STARTSAMPLE-1,1:3*NATOMS)
   LVREF(J1)=VREF(J1+STARTSAMPLE-1)
   CALL POTENTIAL(LCOORDSREF(J1,1:3*NATOMS),DUMMY,GRAD,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
   DO K1=1,3
      DO K2=1,NATOMS
         QTEMP(K1,K2)=LCOORDSREF(J1,3*(K2-1)+K1)
      ENDDO
   ENDDO
   CALL GETORDER(QTEMP,QORDER,NATOMS)
   PRINT '(A,I6,A,G20.10,A,G20.10,A,F15.5,A,I6)','mcpath2> Local reference ',J1,' energy=',DUMMY,' should be ', &
  &                         VREF(J1+STARTSAMPLE-1),' Q=',QORDER, &
  &                        ' frame=',J1+STARTSAMPLE-1
ENDDO

LNPATH=ENDSAMPLE-STARTSAMPLE+1
PRINT '(A,I6)','mcpath2> Setting local number of path frames to ',LNPATH
ALLOCATE(LEOFS(LNPATH))
LEOFS(1:LNPATH)=EOFS(STARTSAMPLE:ENDSAMPLE)
PRINT '(A,I6)','DOBLOCK=',DOBLOCK
PRINT '(A,2I6)','STARTSAMPLE,ENDSAMPLE=',STARTSAMPLE,ENDSAMPLE
!
! Start from the middle or ends of a block.
!
IF (MCPATHSTART.EQ.0) THEN
   MCCOORDS(1:3*NATOMS)=COORDSREF((ENDSAMPLE+STARTSAMPLE)/2,1:3*NATOMS)
   PRINT '(A,I6)','mcpath2> Starting from middle of block at frame ',(ENDSAMPLE+STARTSAMPLE)/2
   IF (MCBIAST) WOLD=-VREF((ENDSAMPLE+STARTSAMPLE)/2)
ELSEIF (MCPATHSTART.LT.0) THEN
   MCCOORDS(1:3*NATOMS)=COORDSREF(STARTSAMPLE,1:3*NATOMS)
   PRINT '(A,I6)','mcpath2> Starting from beginning of block at frame ',STARTSAMPLE
   IF (MCBIAST) WOLD=-VREF(STARTSAMPLE)
ELSE
   MCCOORDS(1:3*NATOMS)=COORDSREF(ENDSAMPLE,1:3*NATOMS)
   PRINT '(A,I6)','mcpath2> Starting from final min at frame ',ENDSAMPLE
   IF (MCBIAST) WOLD=-VREF(ENDSAMPLE)
ENDIF

CALL POTENTIAL(MCCOORDS,VOLD,GRAD,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
PRINT '(A,G20.10)','mcpath2> Initial configuration energy             = ',VOLD

VNEW=VOLD
IF (MCBIAST) WNEW=WOLD

IF (MCBIAST) THEN
   VOLD=VNEW+WNEW
   PRINT '(A,G20.10)','mcpath2> Initial bias energy                      = ',WOLD
ELSE
   WOLD=0.0D0
ENDIF

DO J1=1,NATOMS
   X(J1)=MCCOORDS(3*(J1-1)+1)
   Y(J1)=MCCOORDS(3*(J1-1)+2)
   Z(J1)=MCCOORDS(3*(J1-1)+3)
ENDDO
WRITE(*,'(A)') 'mcpath2> Using hardcoded value (1) as random number seed'
CALL SDPRND(1)

NTOT=0
IACCEPT=0
IMAGEDIFF=0
IREJDIST=0
IREJFRAME=0
IREJCHECK=0
XDISTMIN=1.0D100
XDISTMAX=-1.0D100
MEANBIAS=0.0D0
MEANBIASEQ=0.0D0
QORDERHIST(1:MCPATHBINS)=0.0D0
QORDERVISITS(1:MCPATHBINS)=0.0D0
QORDERMEAN=0.0D0
QORDERSQMEAN=0.0D0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Main loop over steps.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IMCSTEP=0
NEAREST(1:NSTRUCTREF)=0
NEARESTF(1:NPATH)=0
NEARESTFW(1:NPATH)=0.0D0
QFRAMEAV(1:NPATH)=0.0D0
QFRAMESD(1:NPATH)=0.0D0

DO 
   IMCSTEP=IMCSTEP+1
   IF (IMCSTEP.GT.MCPATHEQUIL+MCPATHSTEPS) EXIT
   RECOUNT=.FALSE.
   NFAILS=0
961  CONTINUE
   DO K=1, NATOMS
      XO(K)=X(K)
      YO(K)=Y(K)
      ZO(K)=Z(K)
      MCCOORDSO(3*(K-1)+1)=XO(K)
      MCCOORDSO(3*(K-1)+2)=YO(K)
      MCCOORDSO(3*(K-1)+3)=ZO(K)
   ENDDO
   IF (.FALSE.) THEN ! debug DJW
      CALL POTENTIAL(MCCOORDSO,OPOTEL,GRAD,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
      IF (ABS(OPOTEL-VOLD+WOLD).GT.EDIFFTOL) THEN
         PRINT '(A,G20.10,A,G20.10,A,I10)','mcpath2> ERROR *** energy for coordinates in MCCOORDSO=', &
     &                     OPOTEL,' but Markov energy=',VOLD-WOLD,'IMCSTEP=',IMCSTEP
         STOP
      ENDIF
   ENDIF
!
! Take a Monte Carlo trial step
!
   DO K=1,NATOMS
      RANDOM=DPRAND()
      X(K)=X(K)+(2.0D0*RANDOM-1.0D0)*MCPATHSTEP
      RANDOM=DPRAND()
      Y(K)=Y(K)+(2.0D0*RANDOM-1.0D0)*MCPATHSTEP
      RANDOM=DPRAND()
      Z(K)=Z(K)+(2.0D0*RANDOM-1.0D0)*MCPATHSTEP
   ENDDO
!
! Make sure no frozen atoms have moved
!
   IF (FREEZE) THEN
      DO J1=1,NATOMS
         IF (FROZEN(J1)) THEN
            X(J1)=XO(J1)
            Y(J1)=YO(J1)
            Z(J1)=ZO(J1)
         ENDIF
      ENDDO
   ENDIF
!
! copy post-move coordinates
!
   DO K=1,NATOMS
      MCCOORDS(3*(K-1)+1)=X(K)
      MCCOORDS(3*(K-1)+2)=Y(K)
      MCCOORDS(3*(K-1)+3)=Z(K)
   ENDDO
!
! Calculate distances to references.
! Check that we are within MCPATHDMAX of at least one reference geometry.
! If not, reject and recount. We also reject steps outside the current block.
!
  RECOUNT=.TRUE.
  DISTMIN=1.0D100
  IMAGEMIN=1
  DISTMAX=-1.0D100
  IMAGEMAX=1

! DO J1=1,LSTRUCTREF
!    CALL ALIGN_DECIDE(MCCOORDS,LCOORDSREF(J1,1:3*NATOMS),NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT, &
! &       TWOD,DUMMY3,DIST2,RIGIDBODY,RMAT)
!    IF (DEBUG) PRINT '(A,I10,A,G20.10)','mcpath2> distance for local reference structure ',J1,' is ',DUMMY3
!    IF (DUMMY3.LT.DISTMIN) THEN
!       DISTMIN=DUMMY3
!       IMAGEMIN=J1
!    ENDIF
!    IF (DUMMY3.GT.DISTMAX) THEN
!       DISTMAX=DUMMY3
!       IMAGEMAX=J1
!    ENDIF
! ENDDO

!
! First try previous closest frame and neighbours.
!
! IF (IMCSTEP.GT.1) THEN
!    DO J1=MAX(MAX(1,STARTSAMPLE-MCPATHBLOCK+MCPATHOVER),IMAGEMINO+(STARTSAMPLE-1)-1), &
! &        MIN(MIN(NPATH,ENDSAMPLE+MCPATHBLOCK-MCPATHOVER),IMAGEMINO+(STARTSAMPLE-1)+1)
!       CALL ALIGN_DECIDE(MCCOORDS,COORDSREF(J1,1:3*NATOMS),NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT, &
! &       TWOD,DUMMY3,DIST2,RIGIDBODY,RMAT)
!       IF (DEBUG) PRINT '(A,I10,A,G20.10)','mcpath2> distance for local reference structure ',J1,' is ',DUMMY3
!       PRINT '(A,I10,A,G20.10)','mcpath2> distance for local reference structure ',J1,' is ',DUMMY3
!       IF (DUMMY3.LT.DISTMIN) THEN
!          DISTMIN=DUMMY3
!          IMAGEMIN=J1
!       ENDIF
!       IF (DUMMY3.GT.DISTMAX) THEN
!          DISTMAX=DUMMY3
!          IMAGEMAX=J1
!       ENDIF
!    ENDDO
! ENDIF
!
! Check all frames every NCHECKINT steps.
!
  NCHECKINT=1
  IF ((MOD(IMCSTEP,NCHECKINT).EQ.0).OR.(IMCSTEP.EQ.1)) THEN
     DISTMIN=1.0D100
     DISTMAX=-1.0D100
     DO J1=MAX(1,STARTSAMPLE-MCPATHBLOCK+MCPATHOVER),MIN(NPATH,ENDSAMPLE+MCPATHBLOCK-MCPATHOVER)
        CALL ALIGN_DECIDE(MCCOORDS,COORDSREF(J1,1:3*NATOMS),NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT, &
  &       TWOD,DUMMY3,DIST2,RIGIDBODY,RMAT)
        IF (DEBUG) PRINT '(A,I10,A,G20.10)','mcpath2> distance for local reference structure ',J1,' is ',DUMMY3
!       PRINT '(A,I10,A,G20.10)','mcpath2> distance for local reference structure ',J1,' is ',DUMMY3
        IF (DUMMY3.LT.DISTMIN) THEN
           DISTMIN=DUMMY3
           IMAGEMIN2=J1
        ENDIF
        IF (DUMMY3.GT.DISTMAX) THEN
           DISTMAX=DUMMY3
           IMAGEMAX=J1
        ENDIF
     ENDDO
!    IF (IMAGEMIN2.NE.IMAGEMIN) THEN
!       IF (IMCSTEP.GT.1) THEN
!          PRINT '(A,I10,A,I10)','mcpath2> Closest image assignment corrected from ',IMAGEMIN,' to ',IMAGEMIN2
!             IMAGEDIFF=IMAGEDIFF+1
!       ENDIF
        IMAGEMIN=IMAGEMIN2
!    ENDIF
  ENDIF

  IMAGEMIN=IMAGEMIN-(STARTSAMPLE-1)
  IMAGEMAX=IMAGEMAX-(STARTSAMPLE-1)

  IMAGEMINF=IMAGEMIN+STARTSAMPLE-1
  DISTMINF=DISTMIN

  IF (IMCSTEP.EQ.1) THEN
     IMAGEMINO=IMAGEMIN ! initialise this information at first step
     PRINT '(A,I6,A,G20.10)','mcpath2> Initial configuration is closest to path frame ',IMAGEMINF,' distance=',DISTMINF
     IF ((MCPATHTS.GE.0).AND.((IMAGEMINF.LT.STARTSAMPLE).OR.(IMAGEMINF.GT.ENDSAMPLE))) THEN
        PRINT '(A)','mcath> ERROR *** Initial configuration lies outside the frames defined by the two bracketing stationary points'
        DO K=1, NATOMS
           X(K)=XO(K)
           Y(K)=YO(K)
           Z(K)=ZO(K)
        ENDDO
        DO K=1,NATOMS
           MCCOORDS(3*(K-1)+1)=X(K)
           MCCOORDS(3*(K-1)+2)=Y(K)
           MCCOORDS(3*(K-1)+3)=Z(K)
        ENDDO
        VNEW=VOLD
        WNEW=WOLD
        IMAGEMIN=IMAGEMINO
        IMAGEMINF=IMAGEMINFO
        DO K1=1,3
           DO K2=1,NATOMS
              QTEMP(K1,K2)=MCCOORDS(3*(K2-1)+K1)
           ENDDO
        ENDDO
        CALL GETORDER(QTEMP,QORDERO,NATOMS)
        QORDER=QORDERO
        NFAILS=NFAILS+1
        IF (NFAILS.GT.1000) THEN ! prevent infinite loop!
           PRINT '(A)','mcpath2> Too many failures - give up!'
           STOP
        ENDIF

        GOTO 961
     ENDIF
     IMAGEMINFO=IMAGEMINF
  ELSE
     IMAGEMINF=IMAGEMIN+STARTSAMPLE-1
     DISTMINF=DISTMIN
  ENDIF

  IF (DISTMIN.LT.MCPATHDMAX) RECOUNT=.FALSE.
  IF (DEBUG.AND.(DISTMIN.GT.MCPATHDMAX)) PRINT '(A,2G20.10,L5)','DISTMIN,MCPATHDMAX,RECOUNT=',DISTMIN,MCPATHDMAX,RECOUNT
  IF (IMCSTEP.GT.MCPATHEQUIL) THEN
     IF (RECOUNT) IREJDIST=IREJDIST+1
  ENDIF
  IF ((MCPATHTS.GE.0).AND.((IMAGEMINF.LT.STARTSAMPLE).OR.(IMAGEMINF.GT.ENDSAMPLE))) THEN
     IF (DEBUG) PRINT '(A)','mcpath2> reject: closest frame is outside the range defined by the two end frames'
     RECOUNT=.TRUE.
     IF (IMCSTEP.GT.MCPATHEQUIL) IREJFRAME=IREJFRAME+1
  ENDIF

  IF (IMCSTEP.GT.MCPATHEQUIL) THEN
     IF (DISTMIN.LT.XDISTMIN) XDISTMIN=DISTMIN
     IF (DISTMAX.GT.XDISTMAX) XDISTMAX=DISTMAX
  ENDIF
  IF (DEBUG.AND.(MOD(IMCSTEP-1,MCPATHPRTFRQ).EQ.0)) THEN
     PRINT '(A,2G20.10,A,2I6)','mcpath2> Minimum and maximum reference distances are ',DISTMIN,DISTMAX, &
  &                            ' for local structures ',IMAGEMIN,IMAGEMAX
     PRINT '(A,G20.10,A,I6)','mcpath2> Minimum path frame distance is ',DISTMINF, &
  &                            ' for frame structure ',IMAGEMINF
  ENDIF
!
! At this point all we have done is take a step, energies have not been computed yet and
! also Metropolis check has not been done.
!
! The perturbed coordinates are in both MCCOORDS and X, Y, Z.
! The old coordinates are in XO, YO, ZO.
!
! New and old instantaneous energies will be in    VNEW    VOLD
!
! Compute energy of perturbed coordinates and apply Metropolis check.
! If the step is to be rejected, RECOUNT is set to TRUE
!
  IF (RECOUNT) THEN
     VNEW=VOLD
     WNEW=WOLD
  ELSE
     CALL POTENTIAL(MCCOORDS,VNEW,GRAD,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!    PRINT '(A,G20.10)','mcpath2> after potential call VNEW=',VNEW
     IF (MCBIAST) THEN
        WNEW=-LVREF(IMAGEMIN)
        VNEW=VNEW+WNEW
     ENDIF
     WCOMP=(VNEW-VOLD)/MCPATHTEMP
     DUMMY=MIN(1.0D0,EXP(-WCOMP))
     RANDOM=DPRAND()
     IF (RANDOM.GT.DUMMY) RECOUNT=.TRUE. ! RECOUNT is initialised to .FALSE. at the top of the loop
  ENDIF
  IF ((IMCSTEP.EQ.1).AND.RECOUNT) THEN
     PRINT '(A)','mcpath2> WARNING *** first step was initially rejected - accepting it - check initial step size '
     PRINT '(A,G20.10,A)','mcpath2> Sqrt(3)*maximum step size is ',1.73D0*MCPATHSTEP,' versus 0.1 used in bias potential'
     RECOUNT=.FALSE.
  ENDIF
  VNEWSAVE=VNEW  ! this value is saved so it can be printed if the step is rejected
  WNEWSAVE=WNEW  ! this value is saved so it can be printed if the step is rejected
  IF (RECOUNT) THEN ! reject move
!
! overwriting post-move coordinates X,Y,Z why over writing post-move coordinates X,Y,Z?
! needed because XO, YO and ZO are initialized to X,Y,Z at the top of the loop.
! Post-move coordinates are lost at this point. If needed, add another variable for post-move coords.
!
     DO K=1, NATOMS
        X(K)=XO(K)
        Y(K)=YO(K)
        Z(K)=ZO(K)
     ENDDO
     DO K=1,NATOMS
        MCCOORDS(3*(K-1)+1)=X(K)
        MCCOORDS(3*(K-1)+2)=Y(K)
        MCCOORDS(3*(K-1)+3)=Z(K)
     ENDDO
     VNEW=VOLD
     WNEW=WOLD
     IMAGEMIN=IMAGEMINO
     IMAGEMINF=IMAGEMINFO
     IF (IMCSTEP.EQ.1) THEN ! initialise QORDERO, the saved value
        DO K1=1,3
           DO K2=1,NATOMS
              QTEMP(K1,K2)=MCCOORDS(3*(K2-1)+K1)
           ENDDO
        ENDDO
        CALL GETORDER(QTEMP,QORDERO,NATOMS)
     ENDIF
     QORDER=QORDERO
  ELSE ! accept move
! since move is accepted, no change to MCCOORDS
     IACCEPT=IACCEPT+1
     IMAGEMINO=IMAGEMIN
     IMAGEMINFO=IMAGEMINF
     DO K1=1,3
        DO K2=1,NATOMS
           QTEMP(K1,K2)=MCCOORDS(3*(K2-1)+K1)
        ENDDO
     ENDDO
     QORDERO=QORDER
     CALL GETORDER(QTEMP,QORDER,NATOMS)
  ENDIF ! closes IF (RECOUNT)

!
!  Must not accumulate statistics until we have equilibrated for MCPATHEQUIL steps.
!
  IF (IMCSTEP.GT.MCPATHEQUIL) THEN
     NEAREST(IMAGEMIN)=NEAREST(IMAGEMIN)+1
     NEARESTF(IMAGEMINF)=NEARESTF(IMAGEMINF)+1
     QFRAMEAV(IMAGEMINF)=QFRAMEAV(IMAGEMINF)+QORDER
     QFRAMESD(IMAGEMINF)=QFRAMESD(IMAGEMINF)+QORDER**2
!    PRINT '(A,I10,2G20.10,I10,G20.10)','IMAGEMINF,QORDER,QFRAMEAV,NEARESTF,QFRAMEAV/NEARESTF=', &
! &              IMAGEMINF,QORDER,QFRAMEAV(IMAGEMINF),NEARESTF(IMAGEMINF),QFRAMEAV(IMAGEMINF)/NEARESTF(IMAGEMINF)

     IF (MCBIAST) MEANBIAS=MEANBIAS+EXP((WNEW-MEANBIASEQ)/MCPATHTEMP) 
     IF (.TRUE.) THEN
        IBININDEX=INT((QORDER-MIN(MCPATHQMAX,MCPATHQMIN))/QORDERHISTINT)+1
        IF (MCBIAST) THEN
           IF (IBININDEX.LE.MCPATHBINS) QORDERHIST(IBININDEX)=QORDERHIST(IBININDEX)+EXP((WNEW-MEANBIASEQ)/MCPATHTEMP)
           QORDERMEAN=QORDERMEAN+QORDER*EXP((WNEW-MEANBIASEQ)/MCPATHTEMP)
           QORDERSQMEAN=QORDERSQMEAN+QORDER**2*EXP((WNEW-MEANBIASEQ)/MCPATHTEMP)
           NEARESTFW(IMAGEMINF)=NEARESTFW(IMAGEMINF)+EXP((WNEW-MEANBIASEQ)/MCPATHTEMP)
        ELSE  
           IF (IBININDEX.LE.MCPATHBINS) QORDERHIST(IBININDEX)=QORDERHIST(IBININDEX)+1.0D0
           QORDERMEAN=QORDERMEAN+QORDER
           QORDERSQMEAN=QORDERSQMEAN+QORDER**2
           NEARESTFW(IMAGEMINF)=NEARESTFW(IMAGEMINF)+1.0D0
        ENDIF
        IF (IBININDEX.LE.MCPATHBINS) QORDERVISITS(IBININDEX)=QORDERVISITS(IBININDEX)+1.0D0
     ENDIF
  ELSE
     IF (MCBIAST) MEANBIASEQ=MEANBIASEQ+WNEW
  ENDIF

  SDUMMY='ACC'
  IF (RECOUNT) SDUMMY='REJ'
  IF (MOD(IMCSTEP-1,MCPATHPRTFRQ).EQ.0) THEN
     IF (MCBIAST) THEN
        IF (IMCSTEP.GT.MCPATHEQUIL) THEN
           WRITE(*, '(I10,A,F20.10,A,F20.10,A,3(A,G20.10))') IMCSTEP,' Vn=', VNEWSAVE,' Vo=',VOLD, &
  &            ' ' // SDUMMY,' Wn=',WNEWSAVE,' Wo=',WOLD,' <exp(W)>=',MEANBIAS/IMCSTEP
        ELSE
           WRITE(*, '(I10,A,F20.10,A,F20.10,A,3(A,F20.10))') IMCSTEP,' Vn=', VNEWSAVE,' Vo=',VOLD, &
  &            ' ' // SDUMMY,' Wn=',WNEWSAVE,' Wo=',WOLD
        ENDIF
     ELSE
        WRITE(*, '(I10,A,F20.10,A,F20.10,A)') IMCSTEP,' Vn=', VNEWSAVE,' Vo=',VOLD,' ' // SDUMMY
     ENDIF
  ENDIF
  IF (DEBUG) CALL FLUSH(6)
! IF (VNEWSAVE.GT.170.0D0) THEN
!    PRINT '(A)',' check VNEWSAVE value !'
!    STOP ! DJW debug
! ENDIF
!
!     --- pre-move
!            coord = XO,YO,ZO   energy = VOLD
!     --- post-move and post-Metropolis
!            coord = X, Y, Z    energy = VNEW (if accepted)
!            coord and energy are lost if rejected
!     --- Outcome of Metropolis and other checks
!            RECOUNT = FALSE if accepted
!                    = TRUE  if rejected
!     --- Markov state post Metropolis
!            coord = MCCOORDS      energy = VNEW
!
!
! This is where the step size is adjusted. All adjustments
! are make during equilibration (IMCSTEP<MCPATHEQUIL). Updates are
! made every NUPDATE steps based on acceptance ratio during previous
! NUPDATE steps.
!

! PRINT '(A,4G20.10)','VNEW,VOLD,WNEW,WOLD=',VNEW,VOLD,WNEW,WOLD
  VOLD=VNEW         ! saving current Markov state
  WOLD=WNEW

  IF ((IMCSTEP.LE.MCPATHEQUIL).AND.(MOD(IMCSTEP,NUPDATE).EQ.0)) THEN ! update MC step size if not fixed
     WAC=1.0*IACCEPT/NUPDATE
     IF (WAC.LT.MCPATHACCRATIO-0.1D0) THEN
        MCPATHSTEP=MCPATHSTEP*0.9D0
     ENDIF
     IF (WAC.GT.MCPATHACCRATIO+0.1D0) THEN
        MCPATHSTEP=MCPATHSTEP*1.1D0
     ENDIF
     IACCEPT=0
     IF (MOD(IMCSTEP-1,MCPATHPRTFRQ).EQ.0) WRITE(*,'(A,G20.10,A,G20.10,A,I6,A)')  &
 &                   'mcpath2> adjusting step-size> current step size = ',&
 &                   MCPATHSTEP,' acceptance ratio = ', WAC ,' over ', NUPDATE, ' steps'
  ENDIF ! step size update
  IF (IMCSTEP.EQ.MCPATHEQUIL) THEN
     WRITE(*, '(A)') 'mcpath2> ---------- Equilibration done '
     WRITE(*, '(A,F20.10,A,I10,A,G20.10)') 'mcpath2> Temperature = ', MCPATHTEMP, &
  &                 ' MCSteps = ', IMCSTEP,' MarkovEner = ', VNEW
     WRITE(*, '(A,G20.10,A,G20.10,A,I6,A)') 'mcpath2> Final step size = ',MCPATHSTEP, &
  &'  corresponding to acceptance ratio = ', WAC ,' over previous ', NUPDATE, ' steps'
     WRITE(*, '(A,G20.10)') 'mcpath2>   compare with target acceptance ratio = ', MCPATHACCRATIO
     MEANBIASEQ=MEANBIASEQ/MCPATHEQUIL
     IF (MCBIAST) THEN
        WRITE(*,'(A,G20.10)') 'mcpath2> <W> over equilibration phase=',MEANBIASEQ
     ENDIF
     WRITE(*, '(A)') 'mcpath2> ---------- Starting production run '
  ENDIF

!
! Collect any further averages here.
!
   IF (IMCSTEP.GT.MCPATHEQUIL) THEN
   ENDIF

IF ((IMCSTEP-MCPATHEQUIL.EQ.500000).OR. &
  & (IMCSTEP-MCPATHEQUIL.EQ.1000000).OR. &
  & (IMCSTEP-MCPATHEQUIL.EQ.2000000).OR. &
  & (IMCSTEP-MCPATHEQUIL.EQ.4000000).OR. &
  & (IMCSTEP-MCPATHEQUIL.EQ.8000000).OR. & 
  & (IMCSTEP-MCPATHEQUIL.EQ.16000000).OR. &
  & (IMCSTEP-MCPATHEQUIL.EQ.32000000).OR. &
  & (IMCSTEP-MCPATHEQUIL.EQ.64000000).OR.  &
  & (IMCSTEP-MCPATHEQUIL.EQ.128000000).OR. &
  & (IMCSTEP-MCPATHEQUIL.EQ.256000000).OR. &
  & (IMCSTEP-MCPATHEQUIL.EQ.512000000).OR. &
  & (IMCSTEP-MCPATHEQUIL.EQ.1024000000)) THEN
   LUNIT=GETUNIT()
   WRITE(FNAME,'(I6)') DOBLOCK
   WRITE(FNAME2,'(I6)') INT((IMCSTEP-MCPATHEQUIL)*1.0D0/1.0D6)
   FNAME='mcpath.Q.block' // TRIM(ADJUSTL(FNAME)) // '.' // TRIM(ADJUSTL(FNAME2))
   OPEN(LUNIT,FILE=FNAME,STATUS='UNKNOWN')
   QORDERMEAN2=0.0D0
   QORDERSQMEAN2=0.0D0
   DO J1=1,MCPATHBINS
      IF (MCBIAST) THEN
         WRITE(LUNIT,'(4G20.10)') BINLABELQORDER(J1),QORDERHIST(J1)/MEANBIAS, &
  &                               -LOG(MAX(1.0D-100,QORDERHIST(J1)/MEANBIAS)), &
  &                               QORDERVISITS(J1)/(IMCSTEP-MCPATHEQUIL)
         QORDERMEAN2=QORDERMEAN2+BINLABELQORDER(J1)*QORDERHIST(J1)/MEANBIAS
         QORDERSQMEAN2 =QORDERSQMEAN2+BINLABELQORDER(J1)**2*QORDERHIST(J1)/MEANBIAS
      ELSE
         WRITE(LUNIT,'(4G20.10)') BINLABELQORDER(J1),QORDERHIST(J1)/(IMCSTEP-MCPATHEQUIL), &
  &                               -LOG(MAX(1.0D-100,QORDERHIST(J1)/(IMCSTEP-MCPATHEQUIL))), &
  &                               QORDERVISITS(J1)/(IMCSTEP-MCPATHEQUIL)
         QORDERMEAN2=QORDERMEAN2+BINLABELQORDER(J1)*QORDERHIST(J1)/(IMCSTEP-MCPATHEQUIL)
         QORDERSQMEAN2=QORDERSQMEAN2+BINLABELQORDER(J1)**2*QORDERHIST(J1)/(IMCSTEP-MCPATHEQUIL)
      ENDIF
   ENDDO
   CLOSE(LUNIT)
   PRINT '(A)','mcpath2> Statistics for closest reference structure:'
   PRINT '(2I10,F20.5,A3)',(J1,NEAREST(J1),1.0D2*NEAREST(J1)/(IMCSTEP-MCPATHEQUIL),' % ',J1=1,LSTRUCTREF)

   LUNIT=GETUNIT()
   WRITE(FNAME,'(I6)') DOBLOCK
   WRITE(FNAME2,'(I6)') INT((IMCSTEP-MCPATHEQUIL)*1.0D0/1.0D6)
   FNAME='mcpath.f.block' // TRIM(ADJUSTL(FNAME)) // '.' // TRIM(ADJUSTL(FNAME2))
   OPEN(LUNIT,FILE=FNAME,STATUS='UNKNOWN')
   DO J1=1,NPATH
      IF (MCBIAST) THEN
         WRITE(LUNIT,'(I6,8G20.10)') J1,SLENGTH(J1),EOFS(J1),QFRAME(J1),NEARESTFW(J1)/MEANBIAS, &
  &          -LOG(MAX(1.0D-100,NEARESTFW(J1)/MEANBIAS)),1.0D2*NEARESTF(J1)/(IMCSTEP-MCPATHEQUIL), &
  &  QFRAMEAV(J1)/MAX(1,NEARESTF(J1)),SQRT(MAX(1.0D-100,QFRAMESD(J1)/MAX(1,NEARESTF(J1))-(QFRAMEAV(J1)/MAX(1,NEARESTF(J1)))**2))  
      ELSE
         WRITE(LUNIT,'(I6,8G20.10)') J1,SLENGTH(J1),EOFS(J1),QFRAME(J1),NEARESTFW(J1)/(IMCSTEP-MCPATHEQUIL), &
  &          -LOG(MAX(1.0D-100,NEARESTFW(J1)/(IMCSTEP-MCPATHEQUIL))),1.0D2*NEARESTF(J1)/(IMCSTEP-MCPATHEQUIL), &
  &  QFRAMEAV(J1)/MAX(1,NEARESTF(J1)),SQRT(MAX(1.0D-100,QFRAMESD(J1)/MAX(1,NEARESTF(J1))-(QFRAMEAV(J1)/MAX(1,NEARESTF(J1)))**2))   
      ENDIF
   ENDDO
   CLOSE(LUNIT)
ENDIF

ENDDO ! end of MC loop

WRITE (*,'(A)') 'mcpath2> Exited main MC loop. '

IF (.TRUE.) THEN
   LUNIT=GETUNIT()
   WRITE(FNAME,'(I6)') DOBLOCK
   FNAME='mcpath.Q.block' // TRIM(ADJUSTL(FNAME))
   OPEN(LUNIT,FILE=FNAME,STATUS='UNKNOWN')
   IF (MCBIAST) THEN
      QORDERMEAN=QORDERMEAN/MEANBIAS
      QORDERSQMEAN=QORDERSQMEAN/MEANBIAS
   ELSE
      QORDERMEAN=QORDERMEAN/MCPATHSTEPS
      QORDERSQMEAN=QORDERSQMEAN/MCPATHSTEPS
   ENDIF
   QORDERMEAN2=0.0D0
   QORDERSQMEAN2=0.0D0
   DO J1=1,MCPATHBINS
      IF (MCBIAST) THEN
         WRITE(LUNIT,'(4G20.10)') BINLABELQORDER(J1),QORDERHIST(J1)/MEANBIAS, &
  &                               -LOG(MAX(1.0D-100,QORDERHIST(J1)/MEANBIAS)), &
  &                               QORDERVISITS(J1)/MCPATHSTEPS
         QORDERMEAN2=QORDERMEAN2+BINLABELQORDER(J1)*QORDERHIST(J1)/MEANBIAS
         QORDERSQMEAN2 =QORDERSQMEAN2+BINLABELQORDER(J1)**2*QORDERHIST(J1)/MEANBIAS
      ELSE
         WRITE(LUNIT,'(4G20.10)') BINLABELQORDER(J1),QORDERHIST(J1)/MCPATHSTEPS, &
  &                               -LOG(MAX(1.0D-100,QORDERHIST(J1)/MCPATHSTEPS)), &
  &                               QORDERVISITS(J1)/MCPATHSTEPS
         QORDERMEAN2=QORDERMEAN2+BINLABELQORDER(J1)*QORDERHIST(J1)/MCPATHSTEPS
         QORDERSQMEAN2=QORDERSQMEAN2+BINLABELQORDER(J1)**2*QORDERHIST(J1)/MCPATHSTEPS
      ENDIF
   ENDDO
   CLOSE(LUNIT)

   LUNIT=GETUNIT()
   WRITE(FNAME,'(I6)') DOBLOCK
   FNAME='mcpath.f.block' // TRIM(ADJUSTL(FNAME)) 
   PRINT *,'DOBLOCK,FNAME=',DOBLOCK,FNAME
   OPEN(LUNIT,FILE=FNAME,STATUS='UNKNOWN')
   DO J1=1,NPATH
      IF (MCBIAST) THEN
         WRITE(LUNIT,'(I6,8G20.10)') J1,SLENGTH(J1),EOFS(J1),QFRAME(J1),NEARESTFW(J1)/MEANBIAS, &
  &          -LOG(MAX(1.0D-100,NEARESTFW(J1)/MEANBIAS)),1.0D2*NEARESTF(J1)/MCPATHSTEPS, &
  &      QFRAMEAV(J1)/MAX(1,NEARESTF(J1)),SQRT(MAX(1.0D-100,QFRAMESD(J1)/MAX(1,NEARESTF(J1))-(QFRAMEAV(J1)/MAX(1,NEARESTF(J1)))**2)) 
      ELSE
         WRITE(LUNIT,'(I6,8G20.10)') J1,SLENGTH(J1),EOFS(J1),QFRAME(J1),NEARESTFW(J1)/MCPATHSTEPS, &
  &          -LOG(MAX(1.0D-100,NEARESTFW(J1)/MCPATHSTEPS)),1.0D2*NEARESTF(J1)/MCPATHSTEPS, &
  &      QFRAMEAV(J1)/MAX(1,NEARESTF(J1)),SQRT(MAX(1.0D-100,QFRAMESD(J1)/MAX(1,NEARESTF(J1))-(QFRAMEAV(J1)/MAX(1,NEARESTF(J1)))**2)) 

      ENDIF
   ENDDO
   CLOSE(LUNIT)

ENDIF
!
! Printing summary
!
WRITE(*,'(A,I10,A,I10,A,F15.5,A)') 'mcpath2> ',IACCEPT, ' steps accepted out of ', &
   &      MCPATHSTEPS+MCPATHEQUIL, ' i.e. ',IACCEPT*100.0D0/(MCPATHSTEPS+MCPATHEQUIL),'%'
! WRITE(*,'(A,I10,A,I10)') 'mcpath2> ',IMAGEDIFF, ' changes in minimum distance assigment for check interval ',NCHECKINT
WRITE(*,'(A,G20.10)') 'bspt> Final stepsize ',MCPATHSTEP
PRINT '(A)','mcpath2> Number of production run steps for which configuration was closest to a given reference structure:'
PRINT '(2I10,F20.5,A3)',(J1,NEAREST(J1),1.0D2*NEAREST(J1)/MCPATHSTEPS,' % ',J1=1,LSTRUCTREF)
PRINT '(A)','mcpath2> Number of production run steps for which configuration was closest to a given path frame:'
PRINT *,'STARTSAMPLE=',STARTSAMPLE
PRINT *,'ENDSAMPLE=',ENDSAMPLE
PRINT '(2I10,F20.5,A3)',(J1,NEARESTF(J1),1.0D2*NEARESTF(J1)/MCPATHSTEPS,' % ',J1=STARTSAMPLE,ENDSAMPLE)
NDUMMY=0
DO J1=STARTSAMPLE,ENDSAMPLE
   IF (NEARESTF(J1).EQ.0) NDUMMY=NDUMMY+1
ENDDO
IF (NDUMMY.GT.0) PRINT '(A,I10,A)','mcpath2> WARNING *** ',NDUMMY,' frames were not visited'
PRINT '(A,2G20.10)','mcpath2> Smallest and largest distances from references seen were ',XDISTMIN,XDISTMAX
PRINT '(A,I10,A,F10.2,A)','mcpath2> Production steps rejected on distance criterion=',IREJDIST,' i.e. ', &
  &                        1.0D2*IREJDIST/MCPATHSTEPS,'% '
PRINT '(A,I10,A,F10.2,A)','mcpath2> Production steps rejected on closest path frame criterion=',IREJFRAME,' i.e. ', &
  &                        1.0D2*IREJFRAME/MCPATHSTEPS,'% '
PRINT '(A,3G20.10)','mcpath2> <Q>, <Q^2> and sigma=',QORDERMEAN,QORDERSQMEAN,SQRT(QORDERSQMEAN-QORDERMEAN**2)
PRINT '(A,3G20.10,A)','mcpath2> <Q>, <Q^2> and sigma=',QORDERMEAN2,QORDERSQMEAN2,SQRT(QORDERSQMEAN2-QORDERMEAN2**2),' from Q bins'

DEALLOCATE(LCOORDSREF,LVREF,LEOFS)

IF (PBST.OR.(MCPATHDOBLOCK.GT.0)) STOP ! exit from forked process in this case

END SUBROUTINE SAMPLEBLOCK

SUBROUTINE SUBMITOPTIMJOB(CONNID,DEBUG,DOBLOCK,JOBSTRING,NNODES,NODENAME,USERNAME,WORKINGDIRECTORY)
USE PORFUNCS
USE KEY,ONLY : PBST
IMPLICIT NONE
INTEGER CONNID,DOBLOCK,LUNIT,GETUNIT,NNODES
LOGICAL DEBUG, DUMMYRUNT, LDEBUG
CHARACTER(LEN=10) CONNSTR1, CONNSTR2
CHARACTER(LEN=*) JOBSTRING
CHARACTER(LEN=256) MYJOBSTRING
INTEGER CHILDPID, STATUS
CHARACTER(LEN=80) EXEC
CHARACTER(LEN=80) FPOO
CHARACTER(LEN=80) NODENAME(NNODES)
CHARACTER(LEN=80) :: USERNAME
CHARACTER(LEN=100) :: WORKINGDIRECTORY

DUMMYRUNT=.FALSE.
LDEBUG=.FALSE.
EXEC='OPTIM'

CALL GETPID_SUBR(CHILDPID)
! PRINT '(A,2I8)','in SUBMITOPTIMJOB, CHILDPID,CONNID=',CHILDPID,CONNID
WRITE(CONNSTR1,'(I10)') CHILDPID
WRITE(CONNSTR2,'(I10)') CONNID

FPOO='odata.'//TRIM(ADJUSTL(CONNSTR1)) 
LUNIT=GETUNIT()
OPEN(LUNIT,FILE=TRIM(ADJUSTL(FPOO)),STATUS='UNKNOWN')
WRITE(LUNIT,'(A,I10)') 'MCPATHDOBLOCK ',DOBLOCK
CLOSE(LUNIT)
CALL MYSYSTEM(STATUS,LDEBUG,'sed -e "/PBS/d" -e "/CPUS/d" odata.save >> odata.' // TRIM(ADJUSTL(CONNSTR1)))
IF (STATUS.NE.0) PRINT '(A,I8)','submitoptimjob> WARNING - exit status=',STATUS

MYJOBSTRING=TRIM(ADJUSTL(EXEC))//' '//TRIM(ADJUSTL(CONNSTR1))//' > '//TRIM(ADJUSTL(JOBSTRING))//TRIM(ADJUSTL(CONNSTR1))

IF (LDEBUG) PRINT '(2A)','submitoptimjob> real myjobstring=',TRIM(ADJUSTL(MYJOBSTRING))
IF (DUMMYRUNT) MYJOBSTRING='sleep 5'
IF (LDEBUG) PRINT '(2A)','submitoptimjob> Job string=',TRIM(ADJUSTL(MYJOBSTRING))

IF (PBST) THEN
   CALL SSHSUBMIT(CONNID,STATUS,TRIM(ADJUSTL(MYJOBSTRING)),TRIM(ADJUSTL(CONNSTR1)),LDEBUG,NNODES,NODENAME,USERNAME,WORKINGDIRECTORY)
ELSE
   CALL MYSYSTEM(STATUS,LDEBUG,TRIM(ADJUSTL(MYJOBSTRING)))
ENDIF

IF (STATUS.NE.0) PRINT '(A,I8)','submitoptimjob> WARNING - '//TRIM(ADJUSTL(MYJOBSTRING))//' exit status=',STATUS
CALL EXIT(STATUS)

END SUBROUTINE SUBMITOPTIMJOB

SUBROUTINE MYSYSTEM(STATUS,DEBUG,JOBSTRING)
USE PORFUNCS
IMPLICIT NONE
LOGICAL DEBUG
INTEGER STATUS
CHARACTER(LEN=*) JOBSTRING

IF (DEBUG) WRITE(*,'(A)') 'mysystem> '//TRIM(ADJUSTL(JOBSTRING)) 
! WRITE(*,'(A)') 'mysystem> '//TRIM(ADJUSTL(JOBSTRING)) 
CALL SYSTEM(JOBSTRING,STATUS)

! IF (DEBUG) PRINT '(A,I6)','command '//JOBSTRING//' exit status=',STATUS
! IF (STATUS.NE.0) PRINT '(A,I8)','mysystem> WARNING - '//JOBSTRING//' exit status=',STATUS

RETURN
END SUBROUTINE MYSYSTEM

SUBROUTINE SSHSUBMIT(ICPU,STAT,JOBSTRING,CONNSTR1,LDEBUG,NNODES,NODENAME,USERNAME,WORKINGDIRECTORY)
USE PORFUNCS
USE KEY, ONLY: SSHT
IMPLICIT NONE
CHARACTER(LEN=1000) :: TOTALJOBSTRING
CHARACTER(LEN=200) :: PATHSTRING
INTEGER,INTENT(IN) :: ICPU, NNODES
LOGICAL,INTENT(IN) :: LDEBUG
INTEGER,INTENT(OUT) :: STAT
CHARACTER(LEN=*) :: JOBSTRING,CONNSTR1
CHARACTER(LEN=150) COPYFILES
CHARACTER(LEN=80) NODENAME(NNODES)
CHARACTER(LEN=80) :: USERNAME
CHARACTER(LEN=100) :: WORKINGDIRECTORY

INTEGER :: MYSTATUS
CHARACTER(LEN=80) :: NODE

NODE=TRIM(ADJUSTL(NODENAME(ICPU)))
!  
!  Job submission now changed to use a single system call and one large job string.
!  Putting all the data transfer etc. in the original rsh, so that it runs on
!  the compute node, should avoid any other rcp, rsync or rsh !
!

! PRINT *,'USERNAME=',USERNAME
! PRINT *,'WORKINGDIRECTORY=',WORKINGDIRECTORY

PATHSTRING='/scratch/' // TRIM(ADJUSTL(USERNAME)) // '/' // CONNSTR1

!
! Build up the complete rsh command step by step:
! (1) make the scratch directory on the node. -p flag means no error is generated if the directory already exists.
          IF (SSHT) THEN
             TOTALJOBSTRING= 'ssh ' // TRIM(node) // ' " mkdir -p ' // TRIM(ADJUSTL(PATHSTRING))
          ELSE
             TOTALJOBSTRING= 'rsh ' // TRIM(node) // ' " mkdir -p ' // TRIM(ADJUSTL(PATHSTRING))
          ENDIF
! (2) move to the WORKINGDIRECTORY (saves unpicking the COPYFILES list!)
          TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; cd ' // TRIM(ADJUSTL(WORKINGDIRECTORY))
! (3) copy data from WORKINGDIRECTORY to the scratch directory on the node.
!     Note that if any file is missing an error condition will result, and subsequent commands will fail.
      COPYFILES=' path.xyz '
          TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) &
  &          // ' ; cp -r  *.' // connstr1 // ' ' // TRIM(ADJUSTL(COPYFILES)) // ' ' // TRIM(ADJUSTL(PATHSTRING))
! (4) move to the scratch directory on the node
          TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; cd ' // TRIM(ADJUSTL(PATHSTRING))
! (5) run the OPTIM job
          TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; ' // JOBSTRING
! (6) copy results back
          IF (LDEBUG) THEN ! copy everything back
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; cp  mcpath.* *.' // connstr1 &
   &                      // ' ' // TRIM(ADJUSTL(WORKINGDIRECTORY))
          ELSE 
             TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; cp  mcpath.* *.' // connstr1 &
   &                      // ' ' // TRIM(ADJUSTL(WORKINGDIRECTORY))
          ENDIF
! (7) remove the scratch directory
          TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; rm -r ' // TRIM(ADJUSTL(PATHSTRING)) // ' " '
!     or don;t rm it for debugging
!         TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; ' // ' " '
!         TOTALJOBSTRING=TRIM(ADJUSTL(TOTALJOBSTRING)) // ' ; ls ' // TRIM(ADJUSTL(PATHSTRING)) // ' " '
          IF (LDEBUG) PRINT '(2A)', 'nodes> complete job string: ',TRIM(ADJUSTL(TOTALJOBSTRING))
! (8) submit the job for real
!         PRINT *,' totaljobstring in sshsubmit:' 
!         PRINT *,TOTALJOBSTRING
          CALL SYSTEM(TRIM(ADJUSTL(TOTALJOBSTRING)),MYSTATUS)
          STAT=MYSTATUS
END SUBROUTINE SSHSUBMIT


SUBROUTINE MAKESMOOTH(NPATH,PROBS,GWIDTHS,SMOOTHS,SLENGTH)
IMPLICIT NONE
INTEGER NPATH, J1, J2
DOUBLE PRECISION PROBS(NPATH),GWIDTHS,SMOOTHS(NPATH),SLENGTH(NPATH),NORM1(NPATH)
DOUBLE PRECISION DUMMY,DUMMYS

DUMMY=0.0D0
DO J1=1,NPATH
   DUMMY=DUMMY+PROBS(J1)
ENDDO
PROBS(1:NPATH)=PROBS(1:NPATH)/DUMMY
!
! Gaussian normalisation factor
!
DUMMYS=1.0D0/(SQRT(GWIDTHS)*2.50662827463100D0)
SMOOTHS(1:NPATH)=0.0D0
!
! Normalise each contribution over accessible bins to conserve probabilty at the end points
!
DO J2=1,NPATH
   NORM1(J2)=0.0D0
   DO J1=1,NPATH
      NORM1(J2)=NORM1(J2)+EXP(-(SLENGTH(J1)-SLENGTH(J2))**2/(2.0D0*GWIDTHS))
   ENDDO
ENDDO
DO J1=1,NPATH
   DO J2=1,NPATH
      SMOOTHS(J1)=SMOOTHS(J1)+PROBS(J2)*EXP(-(SLENGTH(J1)-SLENGTH(J2))**2/(2.0D0*GWIDTHS))/NORM1(J2)
   ENDDO
ENDDO

RETURN

END SUBROUTINE MAKESMOOTH


!    GOTO 965 ! skip to the usual fully uncontrained wham from random initial guess - no Pofs.pair output
! 
!    FIXZ=.FALSE.
!    ALLZEROSTOT(1:NPATH)=.TRUE.
!    DO J1=1,NRUNS-1
!       IF (DEBUG) PRINT '(A,2I6,A,2G20.10)','calling mywham with initial normalisations for ',J1,J1+1,' of ',ZNORM(J1:J1+1)
!       ZSAVE=ZNORM(J1)
!       CALL MYWHAM(NRUNS,NPATH,PATHSTR,PROBS,APPSTRING,LINVART,OFAILS,J1,J1+1,ZNORM,FIXZ,ALLZEROS)
!       DO J2=1,NPATH
!          IF (.NOT.ALLZEROS(J2)) ALLZEROSTOT(J2)=.FALSE.
!       ENDDO
!       IF (.NOT.LINVART) THEN
!          ZNORM(J1+1)=ZNORM(J1+1)-ZNORM(J1)+ZSAVE
!          ZNORM(J1)=ZSAVE
!       ELSE
!          STOP
!       ENDIF
!       IF (DEBUG) PRINT '(A,2I6,A,2G20.10)','after mywham shifting znorm vector:'
!       IF (DEBUG) PRINT '(G20.10)',ZNORM(J1:J1+1)
! !     PRINT '(A,2I6,A,2G20.10)','after mywham PROBS:'
! !     PRINT '(G20.10)',PROBS(1:NPATH)
!    ENDDO
!    IF (OFAILS) THEN
!       PRINT '(A)','mcpath> Pofs and PofsQ will not be produced'
!       GOTO 864
!    ENDIF
!    IF (DEBUG) PRINT '(A,2I6,A,2G20.10)','after final pairwise mywham pair up PROBS:'
!    IF (DEBUG) PRINT '(G20.10)',PROBS(1:NPATH)
!    IF (DEBUG) PRINT '(A,2I6,A,2G20.10)','after final pairwise mywham pair up znorm vector:'
!    IF (DEBUG) PRINT '(G20.10)',ZNORM(1:NRUNS)
!    PRINT '(A,2I6,A,2G20.10)','after final pairwise mywham pair up znorm vector:'
!    PRINT '(G20.10)',ZNORM(1:NRUNS)
!    PROBSPAIR(1:NPATH)=PROBS(1:NPATH) ! will serve as initial guess for all window WHAM below.
!    ZNORMPAIR(1:NRUNS)=ZNORM(1:NRUNS) ! will serve as initial guess for all window WHAM below.
! 
!    IF (.NOT.LINVART) THEN
!       DUMMY=PROBS(1)
!       DO J1=1,NPATH
!          PROBS(J1)=EXP(PROBS(J1)-DUMMY)
!          IF (ALLZEROSTOT(J1)) PRINT '(A,I6,A,G20.10)','mcpath2> s bin ',J1,' was not visited, will zero probability of ',PROBS(J1)
!          IF (ALLZEROSTOT(J1)) PROBS(J1)=0.0D0
!          IF (ALLZEROSTOT(J1)) PROBSPAIR(J1)=-50.0D0 ! this is a ln
!       ENDDO
!    ENDIF
! 
!    DUMMY=-1.0D100
!    IF (DEBUG) PRINT '(A)','mcpath2> frames, s, ln(P(s)),P(s) from WHAM:'
!    DO J1=1,NPATH
!       IF (DEBUG) PRINT '(I6,3G20.10)',J1,EOFS(J1),LOG(MAX(1.0D-200,PROBS(J1))),PROBS(J1)
!       IF (PROBS(J1).GT.DUMMY) DUMMY=PROBS(J1)
!    ENDDO
! !
! ! GWIDTH is sigma^2, the variance
! !
!    CALL MAKESMOOTH(NPATH,PROBS,GWIDTHS,SMOOTHSPUP,SLENGTH)
! !
! ! Now the WHAM for all windows using the pair up initial guesses.
! !  
!    PROBS(1:NPATH)=PROBSPAIR(1:NPATH)
!    ZNORM(1:NRUNS)=ZNORMPAIR(1:NRUNS)
! 
!    IF (DEBUG) PRINT '(A,2I6,A,2G20.10)','initial guess for ln PROBS variables from pair up:'
!    IF (DEBUG) PRINT '(G20.10)',PROBS(1:NPATH)
!    IF (DEBUG) PRINT '(A,2I6,A,2G20.10)','initial guess for ln normalisation variables from pair up:'
!    IF (DEBUG) PRINT '(G20.10)',ZNORM(1:NRUNS)
!    FIXZ=.TRUE.
!    ALLZEROSTOT(1:NPATH)=.TRUE.
!    CALL MYWHAM(NRUNS,NPATH,PATHSTR,PROBS,APPSTRING,LINVART,OFAILS,1,NRUNS,ZNORM,FIXZ,ALLZEROSTOT)
!    IF (OFAILS) THEN
!       PRINT '(A)','mcpath> Pofs and PofsQ will not be produced'
!       GOTO 864
!    ENDIF
!    IF (DEBUG) PRINT '(A,2I6,A,2G20.10)','after all windows mywham PROBS for pair up fix:'
!    IF (DEBUG) PRINT '(G20.10)',PROBS(1:NPATH)
!    IF (DEBUG) PRINT '(A,2I6,A,2G20.10)','after all windows mywham znorm vector for pair up fix:'
!    IF (DEBUG) PRINT '(G20.10)',ZNORM(1:NRUNS)
! 
!    IF (.NOT.LINVART) THEN
!       DUMMY=PROBS(1)
!       DO J1=1,NPATH
!          PROBS(J1)=EXP(PROBS(J1)-DUMMY)
!          IF (ALLZEROSTOT(J1)) PRINT '(A,I6,A,G20.10)','mcpath2> s bin ',J1,' was not visited, will zero probability of ',PROBS(J1)
!          IF (ALLZEROSTOT(J1)) PROBS(J1)=0.0D0
!       ENDDO
!    ENDIF
! 
!    DUMMY=-1.0D100
!    IF (DEBUG) PRINT '(A)','mcpath2> frames, s, ln(P(s)),P(s) from WHAM pair up fix:'
!    DO J1=1,NPATH
!       IF (DEBUG) PRINT '(I6,3G20.10)',J1,EOFS(J1),LOG(MAX(1.0D-200,PROBS(J1))),PROBS(J1)
!       IF (PROBS(J1).GT.DUMMY) DUMMY=PROBS(J1)
!    ENDDO
! !
!    CALL MAKESMOOTH(NPATH,PROBS,GWIDTHS,SMOOTHSPFIXUP,SLENGTH)
! !
! ! Now for pair down.
! !
!    DO J1=1,NPATH
!       IF (LINVART) THEN ! positive probabilities!
!          PROBS(J1)=(DPRAND()+1.0D0)/2.0D0
!       ELSE
!          PROBS(J1)=DPRAND()
!       ENDIF
!    ENDDO
!    DO J1=1,NRUNS
!       ZNORM(J1)=DPRAND()
!    ENDDO
!    IF (LINVART) ZNORM(1)=1.0D0
!    IF (LINVART) PROBS(1)=1.0D0
!    IF (.NOT.LINVART) ZNORM(NRUNS)=0.0D0 ! count down
!    IF (.NOT.LINVART) PROBS(NPATH)=0.0D0
!    ALLZEROSTOT(1:NPATH)=.TRUE.
! 
!    FIXZ=.FALSE.
!    DO J1=NRUNS-1,1,-1 ! count down
!       IF (DEBUG) PRINT '(A,2I6,A,2G20.10)','calling mywham with initial normalisations for ',J1,J1+1,' of ',ZNORM(J1:J1+1)
!       ZSAVE=ZNORM(J1+1) ! count down
!       CALL MYWHAM(NRUNS,NPATH,PATHSTR,PROBS,APPSTRING,LINVART,OFAILS,J1,J1+1,ZNORM,FIXZ,ALLZEROS)
!       DO J2=1,NPATH
!          IF (.NOT.ALLZEROS(J2)) ALLZEROSTOT(J2)=.FALSE.
!       ENDDO
!       IF (.NOT.LINVART) THEN
!          ZNORM(J1)=ZNORM(J1)-ZNORM(J1+1)+ZSAVE ! count down
!          ZNORM(J1+1)=ZSAVE ! count down
!       ELSE
!          STOP
!       ENDIF
!       PRINT '(A,2I6,A,2G20.10)','after mywham and shifting znorm vector:'
!       PRINT '(G20.10)',ZNORM(J1:J1+1)
! !     PRINT '(A,2I6,A,2G20.10)','after mywham PROBS:'
! !     PRINT '(G20.10)',PROBS(1:NPATH)
!    ENDDO
!    IF (OFAILS) THEN
!       PRINT '(A)','mcpath> Pofs and PofsQ will not be produced'
!       GOTO 864
!    ENDIF
!    IF (DEBUG) PRINT '(A,2I6,A,2G20.10)','after final pairwise mywham PROBS for pair down:'
!    IF (DEBUG) PRINT '(G20.10)',PROBS(1:NPATH)
!    IF (DEBUG) PRINT '(A,2I6,A,2G20.10)','after final pairwise mywham znorm vector for pair down:'
!    IF (DEBUG) PRINT '(G20.10)',ZNORM(1:NRUNS)
!    PRINT '(A,2I6,A,2G20.10)','after final pairwise mywham znorm vector for pair down:'
!    PRINT '(G20.10)',ZNORM(1:NRUNS)
!    PROBSPAIR(1:NPATH)=PROBS(1:NPATH) ! will serve as initial guess for all window WHAM below.
!    ZNORMPAIR(1:NRUNS)=ZNORM(1:NRUNS) ! will serve as initial guess for all window WHAM below.
! 
!    IF (.NOT.LINVART) THEN
!       DUMMY=PROBS(1)
!       DO J1=1,NPATH
!          PROBS(J1)=EXP(PROBS(J1)-DUMMY)
!          IF (ALLZEROSTOT(J1)) PRINT '(A,I6,A,G20.10)','mcpath2> s bin ',J1,' was not visited, will zero probability of ',PROBS(J1)
!          IF (ALLZEROSTOT(J1)) PROBS(J1)=0.0D0
!          IF (ALLZEROSTOT(J1)) PROBSPAIR(J1)=-50.0D0 ! this is a ln
!       ENDDO
!    ENDIF
! 
!    DUMMY=-1.0D100
!    IF (DEBUG) PRINT '(A)','mcpath2> frames, s, ln(P(s)),P(s) from WHAM for pair down:'
!    DO J1=1,NPATH
!       IF (DEBUG) PRINT '(I6,3G20.10)',J1,EOFS(J1),LOG(MAX(1.0D-200,PROBS(J1))),PROBS(J1)
!       IF (PROBS(J1).GT.DUMMY) DUMMY=PROBS(J1)
!    ENDDO
! !
! ! GWIDTH is sigma^2, the variance
! !
!    GWIDTHS=MCPATHGWS
!    CALL MAKESMOOTH(NPATH,PROBS,GWIDTHS,SMOOTHSPDOWN,SLENGTH)
! 
!    PROBQFROMS(1:MCPATHBINS)=0.0D0
!    DO J1=1,NPATH
!       IBININDEX=INT((QFRAME(J1)-MIN(MCPATHQMAX,MCPATHQMIN))/QORDERHISTINT)+1
!       PROBQFROMS(IBININDEX)=PROBQFROMS(IBININDEX)+QFRAME(J1)*PROBS(J1)
!    ENDDO
!    PROBQFROMS2(1:MCPATHBINS)=0.0D0
!    DO J1=1,NPATH
!       IBININDEX=INT((QFRAMEAV(J1)-MIN(MCPATHQMAX,MCPATHQMIN))/QORDERHISTINT)+1
!       PROBQFROMS2(IBININDEX)=PROBQFROMS2(IBININDEX)+QFRAMEAV(J1)*PROBS(J1)
!    ENDDO
!    PROBQFROMSPAIR(1:MCPATHBINS)=PROBQFROMS(1:MCPATHBINS)
!    PROBQFROMSPAIR2(1:MCPATHBINS)=PROBQFROMS2(1:MCPATHBINS)
! !
! ! Now the WHAM for all windows using the pairwise initial guesses.
! !  
!    PROBS(1:NPATH)=PROBSPAIR(1:NPATH)
!    ZNORM(1:NRUNS)=ZNORMPAIR(1:NRUNS)
! 
!    IF (DEBUG) PRINT '(A,2I6,A,2G20.10)','initial guess for ln PROBS variables from pair down:'
!    IF (DEBUG) PRINT '(G20.10)',PROBS(1:NPATH)
!    IF (DEBUG) PRINT '(A,2I6,A,2G20.10)','initial guess for ln normalisation variables from pair down:'
!    IF (DEBUG) PRINT '(G20.10)',ZNORM(1:NRUNS)
!    FIXZ=.TRUE.
!    ALLZEROSTOT(1:NPATH)=.TRUE.
!    CALL MYWHAM(NRUNS,NPATH,PATHSTR,PROBS,APPSTRING,LINVART,OFAILS,1,NRUNS,ZNORM,FIXZ,ALLZEROSTOT)
!    IF (OFAILS) THEN
!       PRINT '(A)','mcpath> Pofs and PofsQ will not be produced'
!       GOTO 864
!    ENDIF
!    IF (DEBUG) PRINT '(A,2I6,A,2G20.10)','after all windows mywham PROBS:'
!    IF (DEBUG) PRINT '(G20.10)',PROBS(1:NPATH)
!    PRINT '(A,2I6,A,2G20.10)','after all windows mywham initialised by pair down znorm vector:'
!    PRINT '(G20.10)',ZNORM(1:NRUNS)
! 
!    IF (.NOT.LINVART) THEN
!       DUMMY=PROBS(1)
!       DO J1=1,NPATH
!          PROBS(J1)=EXP(PROBS(J1)-DUMMY)
!          IF (ALLZEROSTOT(J1)) PRINT '(A,I6,A,G20.10)','mcpath2> s bin ',J1,' was not visited, will zero probability of ',PROBS(J1)
!          IF (ALLZEROSTOT(J1)) PROBS(J1)=0.0D0
!       ENDDO
!    ENDIF
! 
!    DUMMY=-1.0D100
!    IF (DEBUG) PRINT '(A)','mcpath2> frames, s, ln(P(s)),P(s) from WHAM pair down fix:'
!    DO J1=1,NPATH
!       IF (DEBUG) PRINT '(I6,3G20.10)',J1,EOFS(J1),LOG(MAX(1.0D-200,PROBS(J1))),PROBS(J1)
!       IF (PROBS(J1).GT.DUMMY) DUMMY=PROBS(J1)
!    ENDDO
! !
!    CALL MAKESMOOTH(NPATH,PROBS,GWIDTHS,SMOOTHSPFIXDOWN,SLENGTH)
! 
! !
! ! Now unconstrained starting from previous variables for pair down fix.
! !
!    PROBS(1:NPATH)=PROBSPAIR(1:NPATH) ! reset to initialise fully unconstrained wham
!    ZNORM(1:NRUNS)=ZNORMPAIR(1:NRUNS)
! 
!    IF (DEBUG) PRINT '(A,2I6,A,2G20.10)','initial guess for ln PROBS variables from pair down fix:'
!    IF (DEBUG) PRINT '(G20.10)',PROBS(1:NPATH)
!    IF (DEBUG) PRINT '(A,2I6,A,2G20.10)','initial guess for ln normalisation variables from pair down fix:'
!    IF (DEBUG) PRINT '(G20.10)',ZNORM(1:NRUNS)
! 
! 965 CONTINUE
