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

SUBROUTINE MCPATH
USE KEY
USE COMMONS
USE PORFUNCS
IMPLICIT NONE
INTEGER J1, J2, NSTRUCTREF, HIGHESTREF, NTOT, NUPDATE, IMCSTEP, K, LUNIT, K1, K2, IMAX, NTS, NMIN, J3, NFAILS
INTEGER GETUNIT, IACCEPT, IBININDEX, NDUMMY, IMAGEMIN, IMAGEMAX, IMAGEMINO, IREJDIST, NPATH, NTRIES, IREJFRAME
INTEGER IMAGEMINF, IMAGEMINFO, IMAGEMINSTART, NDIRECTION, STARTSAMPLE, ENDSAMPLE
INTEGER LSTRUCTREF, IREJCHECK, NCHECKSTRUCT, NCHECKPLUS, NCHECKMINUS, DOTS, STARTTS, ENDTS, LNPATH, NCLOSEP, NCLOSEM
INTEGER, ALLOCATABLE :: NEAREST(:), REFSTRUCTFRAME(:), TEMPN(:), NEARESTF(:), LREFSTRUCTFRAME(:)
DOUBLE PRECISION, ALLOCATABLE :: BINLABELQORDER(:), TEMPS(:), TEMPS2(:,:), NEARESTFW(:), SLENGTH(:), QFRAME(:)
DOUBLE PRECISION, ALLOCATABLE :: EOFS(:), PATHFRAMES(:,:), LCOORDSCHECK(:,:), LEOFS(:)
DOUBLE PRECISION, ALLOCATABLE :: QORDERHIST(:), QORDERVISITS(:), NORM1(:)
DOUBLE PRECISION, ALLOCATABLE :: DISTFRAME(:,:), LDISTFRAME(:,:), SMOOTH2D(:,:)
DOUBLE PRECISION, ALLOCATABLE :: COORDSREF(:,:), VREF(:), LCOORDSREF(:,:), LVREF(:)
DOUBLE PRECISION, ALLOCATABLE :: PROBS(:), PROBQ(:), PROBSQ(:,:), ZNORM(:)
DOUBLE PRECISION PROBQFROMS(MCPATHBINS), SMOOTHQFROMS(MCPATHBINS)
DOUBLE PRECISION, ALLOCATABLE :: SMOOTHS(:), SMOOTHQ(:), SMOOTHSQ(:,:)
LOGICAL, ALLOCATABLE :: ALLZEROS(:), ALLZEROQ(:)
DOUBLE PRECISION MCCOORDS(3*NATOMS), MCCOORDSO(3*NATOMS), X(NATOMS), Y(NATOMS), Z(NATOMS), &
  &                                                   XO(NATOMS), YO(NATOMS), ZO(NATOMS), &
  &              GRAD(3*NATOMS), RMS, VNEW, VOLD, OPOTEL, RANDOM, DPRAND, WNEW, WOLD, &
  &              WCOMP, DUMMY, VNEWSAVE, WNEWSAVE, WAC, DIST2, RMAT(3,3), CDUMMY(3*NATOMS), &
  &              DISTMINF, DISTMIN, DISTMAX, XDISTMIN, XDISTMAX, DUMMY3, DPLUS, DMINUS, DISTFMIN, DUMMY4, &
  &              TEMP(MCPATHBINS), DUMMY2, DMIN, MEANBIAS, MEANBIASEQ, QORDERHISTINT, QORDERO, SUM1, SUM2, &
  &              BIASBEST, BIASSTEP, SUMBEST, MAXDEV, DISTCHECK, QORDERMEAN, QORDERSQMEAN, NORM2(MCPATHBINS), &
  &              GWIDTHS, GWIDTHQ, QORDERMEAN2, QORDERSQMEAN2, DUMMYS, DUMMYQ
DOUBLE PRECISION, ALLOCATABLE :: DSTRUCT(:)
LOGICAL RECOUNT, YESNO, BIAST, LASTMIN, LASTTS, PERMDISTSAVE, LPERMDISTSAVE, LINVART, OFAILS, OFAILQ, FIXZ
CHARACTER (LEN=8) SDUMMY
CHARACTER (LEN=80) FNAME, FNAME2
CHARACTER (LEN=5) SSYM
CHARACTER(LEN=10) PATHSTR
DOUBLE PRECISION QTEMP(3,NATOMS),QORDER

GRAD(1:3*NATOMS)=0.0D0
RMS=0.0D0
PERMDISTSAVE=PERMDIST 
LPERMDISTSAVE=LPERMDIST

PRINT '(A,G20.10)','mcpath> Monte Carlo sampling around path.xyz file, canonical temperature=',MCPATHTEMP
PRINT '(A,G20.10)','mcpath> Distance limit on sampling from reference points on the pathway= ',MCPATHDMAX
PRINT '(A,G20.10)','mcpath> Maximum initial MC Cartesian step size=                          ',MCPATHSTEP
PRINT '(A,G20.10)','mcpath> Target acceptance ratio=                                         ',MCPATHACCRATIO
PRINT '(A,I10)','mcpath> Number of bins for distance histograms=                          ',MCPATHBINS
PRINT '(A,I10,A,F15.5,A)','mcpath> Number of equilibration steps=                                   ',MCPATHEQUIL, &
  &                     ' = ',MCPATHEQUIL/1.0D6,' Ms'
PRINT '(A,I10,A,F15.5,A)','mcpath> Number of MC production sampling steps=                          ',MCPATHSTEPS, &
  &                     ' = ',MCPATHSTEPS/1.0D6,' Ms'
PRINT '(A,I10)','mcpath> Print frequency interval=                                        ',MCPATHPRTFRQ
PRINT '(A,G20.10)','mcpath> Bottom of order parameter range=                                 ',MCPATHQMIN
PRINT '(A,G20.10)','mcpath> Top of order parameter range=                                    ',MCPATHQMAX
IF (MCBIAST) THEN
   PRINT '(A)', 'mcpath> Bias function will be costructed and applied'
   PRINT '(A,G20.10)','mcpath> Starting exponential bias factor=                                         ',BIASFAC
   PRINT '(A,G20.10)','mcpath> Additional reference structures will be added for bias deviation          ',MCADDDEV
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
   PRINT '(A)','mcpath> ERROR *** no path.xyz file'
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
20 PRINT '(A,I10)','mcpath> Number of frames in path.xyz=',NPATH
CLOSE(LUNIT)
ALLOCATE(EOFS(NPATH),PATHFRAMES(NPATH,3*NATOMS),NORM1(NPATH))
ALLOCATE(PROBS(NPATH),PROBQ(MCPATHBINS),PROBSQ(NPATH,MCPATHBINS))
ALLOCATE(ALLZEROS(NPATH),ALLZEROQ(MCPATHBINS))
ALLOCATE(SMOOTHS(NPATH),SMOOTHQ(MCPATHBINS),SMOOTHSQ(NPATH,MCPATHBINS))
ALLOCATE(SLENGTH(NPATH),QFRAME(NPATH))
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
163 CONTINUE
   NTRIES=NTRIES+1
   CALL ALIGN_DECIDE(PATHFRAMES(J1-1,1:3*NATOMS),PATHFRAMES(J1,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
   IF (DEBUG) PRINT '(A,I6,A,I6,A,G20.10)','mcpath> Distance between frames ',J1,' and ',J1-1,' after alignment is ',DUMMY
   IF (DEBUG) PRINT '(3(A,I6),A,G20.10)','mcpath> Distance between frames ',J1,' and ',J1-1, &
  &                 ' after alignment attempt ',NTRIES,' is ',DUMMY
   IF ((DUMMY.GT.2.0D0*MAXBFGS).AND.(NTRIES.LT.10000)) GOTO 163
   PRINT '(3(A,I6),A,G20.10)','mcpath> Distance between frames ',J1,' and ',J1-1,' after alignment attempt ',NTRIES,' is ',DUMMY
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
PRINT '(A,G20.10)','mcpath> largest distance between frames in path.xyz is now',DUMMY2
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
   IF (DEBUG) PRINT '(A,I6,A,I6,A,G20.10)','mcpath> Distance between frames ',J1,' and ',J1-1,' is now ',DUMMY
   IF (DUMMY.GT.DUMMY2) DUMMY2=DUMMY
ENDDO
PRINT '(A,G20.10)','mcpath> largest distance between frames in path.xyz is now',DUMMY2

IF (DUMMY2.GT.1.0D0) STOP !!! debug DJW

NSTRUCTREF=2
DO J1=2,NPATH-1
   IF (EOFS(J1-1)+EDIFFTOL*10.0D0<EOFS(J1).AND.EOFS(J1)>EOFS(J1+1)+EDIFFTOL*10.0D0) NSTRUCTREF=NSTRUCTREF+1
   IF (EOFS(J1-1)+EDIFFTOL*10.0D0>EOFS(J1).AND.EOFS(J1)<EOFS(J1+1)+EDIFFTOL*10.0D0) NSTRUCTREF=NSTRUCTREF+1
ENDDO

ALLOCATE(COORDSREF(NSTRUCTREF,3*NATOMS),VREF(NSTRUCTREF),NEAREST(NSTRUCTREF),REFSTRUCTFRAME(NSTRUCTREF))
ALLOCATE(DSTRUCT(NSTRUCTREF))
ALLOCATE(QORDERHIST(MCPATHBINS),BINLABELQORDER(MCPATHBINS),QORDERVISITS(MCPATHBINS))
ALLOCATE(NEARESTF(NPATH),NEARESTFW(NPATH))

NSTRUCTREF=1
NMIN=1
NTS=0
VREF(1)=EOFS(1) ! first minimum
REFSTRUCTFRAME(1)=1
COORDSREF(1,1:3*NATOMS)=PATHFRAMES(1,1:3*NATOMS)
HIGHESTREF=1
LASTMIN=.TRUE.
LASTTS=.FALSE.
PRINT '(A,I10,A,G20.10)','mcpath> minimum assumed             for frame ',1,' energy=',VREF(NSTRUCTREF)
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
   PRINT '(A,G20.10)','mcpath> Q=',QORDER
ENDIF

DO J1=2,NPATH-1
   IF (EOFS(J1-1)+EDIFFTOL*10.0D0<EOFS(J1).AND.EOFS(J1)>EOFS(J1+1)+EDIFFTOL*10.0D0) THEN
      NSTRUCTREF=NSTRUCTREF+1
      REFSTRUCTFRAME(NSTRUCTREF)=J1
      VREF(NSTRUCTREF)=EOFS(J1) 
      IF (VREF(NSTRUCTREF).GT.VREF(HIGHESTREF)) HIGHESTREF=NSTRUCTREF
      COORDSREF(NSTRUCTREF,1:3*NATOMS)=PATHFRAMES(J1,1:3*NATOMS)
      CALL ALIGN_DECIDE(COORDSREF(NSTRUCTREF-1,1:3*NATOMS),COORDSREF(NSTRUCTREF,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
      PRINT '(A,I10,2(A,G20.10))','mcpath> transition state identified for frame ',J1,' energy=', &
  &         VREF(NSTRUCTREF),' distance to previous minimum=',DUMMY
      NTS=NTS+1
      IF (.NOT.LASTMIN) THEN
         PRINT '(A)','mcpath> ERROR *** previous stationary point identified was not a minimum'
         STOP
      ENDIF
      LASTMIN=.FALSE.
      LASTTS=.TRUE.
      IF (.TRUE.) THEN
         DO K1=1,3
            DO K2=1,NATOMS
               QTEMP(K1,K2)=COORDSREF(NSTRUCTREF,3*(K2-1)+K1)
            ENDDO
         ENDDO
         CALL GETORDER(QTEMP,QORDER,NATOMS)
         PRINT '(A,G20.10)','mcpath> Q=',QORDER
      ENDIF
   ENDIF
   IF ((EOFS(J1-1)>EOFS(J1).AND.EOFS(J1)<=EOFS(J1+1)).AND.(ABS(EOFS(J1)-VREF(NSTRUCTREF)).GT.EDIFFTOL)) THEN
      NSTRUCTREF=NSTRUCTREF+1
      REFSTRUCTFRAME(NSTRUCTREF)=J1
      VREF(NSTRUCTREF)=EOFS(J1) 
      IF (VREF(NSTRUCTREF).GT.VREF(HIGHESTREF)) HIGHESTREF=NSTRUCTREF
      COORDSREF(NSTRUCTREF,1:3*NATOMS)=PATHFRAMES(J1,1:3*NATOMS)
      CALL ALIGN_DECIDE(COORDSREF(NSTRUCTREF-1,1:3*NATOMS),COORDSREF(NSTRUCTREF,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
      PRINT '(A,I10,2(A,G20.10))','mcpath> minimum identified for frame ',J1,' energy=', &
  &         VREF(NSTRUCTREF),' distance to previous transition state=',DUMMY
      NMIN=NMIN+1
      IF (.NOT.LASTTS) THEN
         PRINT '(A)','mcpath> ERROR *** previous stationary point identified was not a transition state'
         STOP
      ENDIF
      LASTTS=.FALSE.
      LASTMIN=.TRUE.
      IF (.TRUE.) THEN
         DO K1=1,3
            DO K2=1,NATOMS
               QTEMP(K1,K2)=COORDSREF(NSTRUCTREF,3*(K2-1)+K1)
            ENDDO
         ENDDO
         CALL GETORDER(QTEMP,QORDER,NATOMS)
         PRINT '(A,G20.10)','mcpath> Q=',QORDER
      ENDIF
   ENDIF
ENDDO
!
! This interval will depend upon the order parameter. 
!
QORDERHISTINT=(MAX(MCPATHQMAX,MCPATHQMIN)-MIN(MCPATHQMAX,MCPATHQMIN))/MCPATHBINS 
DO J1=1,MCPATHBINS
! these Q values correspond to the bottom of the Q bins
   BINLABELQORDER(J1)=MIN(MCPATHQMAX,MCPATHQMIN)+QORDERHISTINT*(J1-1.0D0) 
ENDDO

INQUIRE(FILE='redo',EXIST=YESNO)

IF (MCPATHEQUIL.EQ.0) THEN
   PRINT '(A)','mcpath> zero equilibration steps specified - reading previous results from disc'
   GOTO 963
ENDIF
IF (YESNO) THEN
   PRINT '(A)','mcpath> File redo detected in working directory - reading previous results from disc'
   GOTO 963
ENDIF

! IF (NTS.EQ.1) THEN ! just one path - permutations should be OK already
!    PERMDISTSAVE=.FALSE. 
!    LPERMDISTSAVE=.FALSE.
!    PERMDIST=.FALSE. 
!    LPERMDIST=.FALSE.
! ENDIF
IF (ABS(EOFS(NPATH)-VREF(NSTRUCTREF)).GT.EDIFFTOL) THEN
   NSTRUCTREF=NSTRUCTREF+1
   REFSTRUCTFRAME(NSTRUCTREF)=NPATH
   VREF(NSTRUCTREF)=EOFS(NPATH) ! last minimum
   COORDSREF(NSTRUCTREF,1:3*NATOMS)=PATHFRAMES(NPATH,1:3*NATOMS)
      CALL ALIGN_DECIDE(COORDSREF(NSTRUCTREF-1,1:3*NATOMS),COORDSREF(NSTRUCTREF,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
   PRINT '(A,I10,2(A,G20.10))','mcpath> minimum assumed    for frame ',NPATH,' energy=', &
  &         VREF(NSTRUCTREF),' distance to previous transition state=',DUMMY
   NMIN=NMIN+1
   IF (.NOT.LASTTS) THEN
      PRINT '(A)','mcpath> ERROR *** previous stationary point identified was not a transition state'
      STOP
   ENDIF
   IF (.TRUE.) THEN
      DO K1=1,3
         DO K2=1,NATOMS
            QTEMP(K1,K2)=COORDSREF(NSTRUCTREF,3*(K2-1)+K1)
         ENDDO
      ENDDO
      CALL GETORDER(QTEMP,QORDER,NATOMS)
      PRINT '(A,G20.10)','mcpath> Q=',QORDER
   ENDIF
ENDIF
PRINT '(4(A,I10))','mcpath> Total number of reference structures=',NSTRUCTREF,' with ',NTS,' ts and ',NMIN,' minima'
CALL FLUSH(6)
IF (MCPATHTS.GT.NTS) THEN
   PRINT '(A,I6,A)','mcpath> ERROR *** odata file asks for ts ',MCPATHTS,' which does not exist'
   STOP
ENDIF
!
! Check whether any of the reference structures are separated by distances greater than MCPATHDMAX
! If so, add more references! We might want to tie the mcpath simulation closer to the path than
! the separation between two stationary points.
!
654 CONTINUE
DO J1=1,NSTRUCTREF
!
! Check distance to original frame is zero!
!
   CALL ALIGN_DECIDE(COORDSREF(J1,1:3*NATOMS),PATHFRAMES(REFSTRUCTFRAME(J1),1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
   PRINT '(A,I6,A,I6,A,F20.10)','mcpath> distance between reference structure ',J1,' and corresponding path frame ', &
  &                 REFSTRUCTFRAME(J1),' is ',DUMMY
!
! Find closest references forwards and backwards.
!
   DPLUS=1.0D100
   DMINUS=1.0D100
   DO J2=1,NSTRUCTREF
      IF (J2.EQ.J1) CYCLE

      CALL ALIGN_DECIDE(COORDSREF(J1,1:3*NATOMS),COORDSREF(J2,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
      IF (DEBUG) PRINT '(2(A,I10),A,G20.10,A,I6,A,I6)','mcpath> distance between reference structures ',J1,' and ',J2,' is ',DUMMY, &
  &   ' frames: ',REFSTRUCTFRAME(J1),' and ',REFSTRUCTFRAME(J2)
      IF (REFSTRUCTFRAME(J1).LT.REFSTRUCTFRAME(J2)) THEN
         IF (DUMMY.LT.DPLUS) THEN
            NCLOSEP=J2
            DPLUS=DUMMY
         ENDIF
      ELSE
         IF (DUMMY.LT.DMINUS) THEN
            NCLOSEM=J2
            DMINUS=DUMMY
         ENDIF
      ENDIF
   ENDDO
   IF (J1.EQ.1) DMINUS=-1.0D0
   IF (J1.EQ.NMIN+NTS) DPLUS=-1.0D0
   IF (J1.NE.NMIN+NTS) PRINT '(A,I6,A,I6,A,G20.10,A,I6)','mcpath> Closest plus reference for ',J1,' is ', &
  &                           NCLOSEP,' at distance ',DPLUS, &
  &                                  ' and frame ',REFSTRUCTFRAME(NCLOSEP)
   IF (J1.NE.1) PRINT '(A,I6,A,I6,A,G20.10,A,I6)','mcpath> Closest minus reference for ',J1,' is ',NCLOSEM,' at distance ',DMINUS, &
  &                                  ' and frame ',REFSTRUCTFRAME(NCLOSEM)

   IF (DPLUS.GT.MCPATHADDREF) THEN ! add another reference structure between J1 and NCLOSEP
      DUMMY2=1.0D100
      NDUMMY=-1
      rloopp: DO J2=MIN(REFSTRUCTFRAME(J1),REFSTRUCTFRAME(NCLOSEP))+1,MAX(REFSTRUCTFRAME(J1),REFSTRUCTFRAME(NCLOSEP))-1 
!
! Check whether this frame is already a reference.
!
         DO J3=1,NSTRUCTREF
            IF (REFSTRUCTFRAME(J3).EQ.J2) CYCLE rloopp
         ENDDO
         CALL ALIGN_DECIDE(COORDSREF(J1,1:3*NATOMS),PATHFRAMES(J2,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
         DUMMY3=DUMMY
         CALL ALIGN_DECIDE(COORDSREF(NCLOSEP,1:3*NATOMS),PATHFRAMES(J2,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
         DUMMY3=ABS(DUMMY3-DUMMY)
!        PRINT '(A,G20.10,A,I6)','mcpath> Difference of distances is ',DUMMY3,' for frame ',J2
         IF (DUMMY3.LT.DUMMY2) THEN
            DUMMY2=DUMMY3
            NDUMMY=J2
            IF (DEBUG) PRINT '(A,G20.10,A,I6)','mcpath> Smallest difference of distances so far is ',DUMMY2,' for frame ',NDUMMY
         ENDIF
      ENDDO rloopp

      IF (NDUMMY.EQ.-1) THEN
         PRINT '(A,I6)','mcpath> Cannot satisfy required conditions - increase MCPATHDMAX'
         STOP
      ENDIF

      PRINT '(A,I6,A,G20.10)','mcpath> Adding additional reference from frame ',NDUMMY,' to satisfy maximum separation of ',MCPATHADDREF

      ALLOCATE(TEMPS(NSTRUCTREF))
      TEMPS(1:NSTRUCTREF)=VREF(1:NSTRUCTREF)
      DEALLOCATE(VREF)
      ALLOCATE(VREF(NSTRUCTREF+1))
      VREF(1:NSTRUCTREF)=TEMPS(1:NSTRUCTREF)
      DEALLOCATE(TEMPS)

      DEALLOCATE(DSTRUCT)
      ALLOCATE(DSTRUCT(NSTRUCTREF+1))
            
      DEALLOCATE(NEAREST)
      ALLOCATE(NEAREST(NSTRUCTREF+1))

      ALLOCATE(TEMPN(NSTRUCTREF))
      TEMPN(1:NSTRUCTREF)=REFSTRUCTFRAME(1:NSTRUCTREF)
      DEALLOCATE(REFSTRUCTFRAME)
      ALLOCATE(REFSTRUCTFRAME(NSTRUCTREF+1))
      REFSTRUCTFRAME(1:NSTRUCTREF)=TEMPN(1:NSTRUCTREF)
      DEALLOCATE(TEMPN)

      ALLOCATE(TEMPS2(NSTRUCTREF,3*NATOMS))
      TEMPS2(1:NSTRUCTREF,1:3*NATOMS)=COORDSREF(1:NSTRUCTREF,1:3*NATOMS)
      DEALLOCATE(COORDSREF)
      ALLOCATE(COORDSREF(NSTRUCTREF+1,3*NATOMS))
      COORDSREF(1:NSTRUCTREF,1:3*NATOMS)=TEMPS2(1:NSTRUCTREF,1:3*NATOMS)
      DEALLOCATE(TEMPS2)

      NSTRUCTREF=NSTRUCTREF+1
!
! Must not reorder minima and transition states from reference frames. Procedures below
! assume additional frames are added at the end!
!
      VREF(NSTRUCTREF)=EOFS(NDUMMY)
      REFSTRUCTFRAME(NSTRUCTREF)=NDUMMY
      COORDSREF(NSTRUCTREF,1:3*NATOMS)=PATHFRAMES(NDUMMY,1:3*NATOMS)
      GOTO 654
   ENDIF
   IF (DMINUS.GT.MCPATHADDREF) THEN ! add another reference structure between J1 and NCLOSEM
      DUMMY2=1.0D100
      NDUMMY=-1
      rloopm: DO J2=MIN(REFSTRUCTFRAME(J1),REFSTRUCTFRAME(NCLOSEM))+1,MAX(REFSTRUCTFRAME(J1),REFSTRUCTFRAME(NCLOSEM))-1 
!
! Check whether this frame is already a reference.
!
         DO J3=1,NSTRUCTREF
            IF (REFSTRUCTFRAME(J3).EQ.J2) CYCLE rloopm
         ENDDO
         CALL ALIGN_DECIDE(COORDSREF(J1,1:3*NATOMS),PATHFRAMES(J2,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
         DUMMY3=DUMMY
         CALL ALIGN_DECIDE(COORDSREF(NCLOSEM,1:3*NATOMS),PATHFRAMES(J2,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
         DUMMY3=ABS(DUMMY3-DUMMY)
!        PRINT '(A,G20.10,A,I6)','mcpath> Difference of distances is ',DUMMY3,' for frame ',J2
         IF (DUMMY3.LT.DUMMY2) THEN
            DUMMY2=DUMMY3
            NDUMMY=J2
            IF (DEBUG) PRINT '(A,G20.10,A,I6)','mcpath> Smallest difference of distances so far is ',DUMMY2,' for frame ',NDUMMY
         ENDIF
      ENDDO rloopm

      IF (NDUMMY.EQ.-1) THEN
         PRINT '(A,I6)','mcpath> Cannot satisfy required conditions - increase MCPATHDMAX'
         STOP
      ENDIF

      PRINT '(A,I6,A,G20.10)','mcpath> Adding additional reference from frame ',NDUMMY,' to satisfy maximum separation of ',MCPATHADDREF

      ALLOCATE(TEMPS(NSTRUCTREF))
      TEMPS(1:NSTRUCTREF)=VREF(1:NSTRUCTREF)
      DEALLOCATE(VREF)
      ALLOCATE(VREF(NSTRUCTREF+1))
      VREF(1:NSTRUCTREF)=TEMPS(1:NSTRUCTREF)
      DEALLOCATE(TEMPS)

      DEALLOCATE(DSTRUCT)
      ALLOCATE(DSTRUCT(NSTRUCTREF+1))
            
      DEALLOCATE(NEAREST)
      ALLOCATE(NEAREST(NSTRUCTREF+1))

      ALLOCATE(TEMPN(NSTRUCTREF))
      TEMPN(1:NSTRUCTREF)=REFSTRUCTFRAME(1:NSTRUCTREF)
      DEALLOCATE(REFSTRUCTFRAME)
      ALLOCATE(REFSTRUCTFRAME(NSTRUCTREF+1))
      REFSTRUCTFRAME(1:NSTRUCTREF)=TEMPN(1:NSTRUCTREF)
      DEALLOCATE(TEMPN)

      ALLOCATE(TEMPS2(NSTRUCTREF,3*NATOMS))
      TEMPS2(1:NSTRUCTREF,1:3*NATOMS)=COORDSREF(1:NSTRUCTREF,1:3*NATOMS)
      DEALLOCATE(COORDSREF)
      ALLOCATE(COORDSREF(NSTRUCTREF+1,3*NATOMS))
      COORDSREF(1:NSTRUCTREF,1:3*NATOMS)=TEMPS2(1:NSTRUCTREF,1:3*NATOMS)
      DEALLOCATE(TEMPS2)

      NSTRUCTREF=NSTRUCTREF+1
!
! Must not reorder minima and transition states from reference frames. Procedures below
! assume additional frames are added at the end!
!
!     DO J2=NSTRUCTREF,J1+1,-1
!        VREF(J2)=VREF(J2-1)
!        REFSTRUCTFRAME(J2)=REFSTRUCTFRAME(J2-1)
!        COORDSREF(J2,1:3*NATOMS)=COORDSREF(J2-1,1:3*NATOMS)
!     ENDDO
      VREF(NSTRUCTREF)=EOFS(NDUMMY)
      REFSTRUCTFRAME(NSTRUCTREF)=NDUMMY
      COORDSREF(NSTRUCTREF,1:3*NATOMS)=PATHFRAMES(NDUMMY,1:3*NATOMS)
      GOTO 654
   ENDIF
ENDDO
PRINT '(A,I6)','mcpath> Total number of reference structures is now: ',NSTRUCTREF
!
! Align path frames if necessary. Should not be executed, since PERMDISTSAVE and LPERMDISTSAVE
! were turned off above!
!
IF (PERMDISTSAVE.OR.LPERMDISTSAVE) THEN
   CALL ALIGN_DECIDE(COORDSREF(1,1:3*NATOMS),PATHFRAMES(1,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
   PRINT '(A,F20.10)','mcpath> first distance=',DUMMY
   DO J1=2,NPATH
      CALL ALIGN_DECIDE(PATHFRAMES(J1-1,1:3*NATOMS),PATHFRAMES(J1,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
      PRINT '(A,I10,F20.10)','mcpath> frame and distance=',J1,DUMMY
   ENDDO
!
! Need to resave COORDSREF to get consistent permutational alignment.
!
   DO J1=1,NSTRUCTREF 
      COORDSREF(J1,1:3*NATOMS)=PATHFRAMES(REFSTRUCTFRAME(J1),1:3*NATOMS)
   ENDDO
   PERMDIST=.FALSE.
   LPERMDIST=.FALSE. ! freeze permutations of references - should already be optimal now
ENDIF
IF (.NOT.MCBIAST) GOTO 987
765 CONTINUE ! jump back here after adding extra reference structures
!
! Should only need to calculate distances for last reference if we have just added one.
! This was not working!! DISTFRAME was never set.
!
IF (.NOT.ALLOCATED(DISTFRAME))   ALLOCATE(DISTFRAME(NPATH,NSTRUCTREF))
!
! Shepard interpolation using just the pe at
! each reference structure and the distance.
!
! IF (NSTRUCTREF.EQ.NMIN+NTS) THEN
   DO J2=1,NSTRUCTREF
      DO J1=1,NPATH
         CALL ALIGN_DECIDE(COORDSREF(J2,1:3*NATOMS),PATHFRAMES(J1,1:3*NATOMS),NATOMS,DEBUG, &
  &                       PARAM1,PARAM2,PARAM3,BULKT,TWOD,DISTFRAME(J1,J2),DIST2,RIGIDBODY,RMAT)
      ENDDO
   ENDDO
! ELSE
!    DO J1=1,NPATH
!       CALL ALIGN_DECIDE(COORDSREF(NSTRUCTREF,1:3*NATOMS),PATHFRAMES(J1,1:3*NATOMS),NATOMS,DEBUG, &
!   &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DISTFRAME(J1,NSTRUCTREF),DIST2,RIGIDBODY,RMAT)
!    ENDDO
! ENDIF

BIASSTEP=1.0D0
CALL SUMSQ(BIASFAC,DISTFRAME,NPATH,NSTRUCTREF,SUM1,VREF,EOFS,MAXDEV,IMAX)
BIASBEST=BIASFAC
SUMBEST=SUM1
NDUMMY=0
 DO WHILE (NDUMMY<51)
   NDUMMY=NDUMMY+1
   BIASFAC=MAX(BIASFAC*(1.0D0+BIASSTEP),0.1D0)
   CALL SUMSQ(BIASFAC,DISTFRAME,NPATH,NSTRUCTREF,SUM2,VREF,EOFS,MAXDEV,IMAX)
   IF (DEBUG) PRINT '(4(A,G20.10))','mcpath> For bias exponent ',BIASFAC,' sum=',SUM2,' best=',SUMBEST,' step=',BIASSTEP
   IF (SUM2.LT.SUMBEST) THEN
      BIASBEST=BIASFAC
      SUMBEST=SUM2
   ELSE IF (SUM2.LT.SUM1) THEN
      BIASSTEP=BIASSTEP*1.01D0
   ELSE
      BIASSTEP=-BIASSTEP/1.01D0
   ENDIF
   IF (ABS(BIASSTEP).GT.0.5D0) BIASSTEP=SIGN(0.5D0,BIASSTEP)
   SUM1=SUM2
   IF (ABS(BIASSTEP).LT.1.0D-10) EXIT
ENDDO
PRINT '(2(A,G20.10))','mcpath> Smallest deviation=',SUMBEST,' for bias exponent=',BIASBEST
BIASFAC=BIASBEST
CALL SUMSQ(BIASFAC,DISTFRAME,NPATH,NSTRUCTREF,SUM2,VREF,EOFS,MAXDEV,IMAX)
PRINT '(A,G20.10,A,I10)','mcpath> Largest squared deviation=',MAXDEV,' for frame ',IMAX

IF (MAXDEV.GT.MCADDDEV) THEN

   ALLOCATE(TEMPS(NSTRUCTREF))
   TEMPS(1:NSTRUCTREF)=VREF(1:NSTRUCTREF)
   DEALLOCATE(VREF)
   ALLOCATE(VREF(NSTRUCTREF+1))
   VREF(1:NSTRUCTREF)=TEMPS(1:NSTRUCTREF)
   VREF(1+NSTRUCTREF)=EOFS(IMAX)
   DEALLOCATE(TEMPS)

   DEALLOCATE(DSTRUCT)
   ALLOCATE(DSTRUCT(NSTRUCTREF+1))

   DEALLOCATE(NEAREST)
   ALLOCATE(NEAREST(NSTRUCTREF+1))

   ALLOCATE(TEMPN(NSTRUCTREF))
   TEMPN(1:NSTRUCTREF)=REFSTRUCTFRAME(1:NSTRUCTREF)
   DEALLOCATE(REFSTRUCTFRAME)
   ALLOCATE(REFSTRUCTFRAME(NSTRUCTREF+1))
   REFSTRUCTFRAME(1:NSTRUCTREF)=TEMPN(1:NSTRUCTREF)
   REFSTRUCTFRAME(1+NSTRUCTREF)=IMAX
   PRINT '(A,I6,A,I6,A)','mcpath> Adding reference structure number ',NSTRUCTREF+1,' from frame ',IMAX,' to improve bias potential'
   DEALLOCATE(TEMPN)

   ALLOCATE(TEMPS2(NSTRUCTREF,3*NATOMS))
   TEMPS2(1:NSTRUCTREF,1:3*NATOMS)=COORDSREF(1:NSTRUCTREF,1:3*NATOMS)
   DEALLOCATE(COORDSREF)
   ALLOCATE(COORDSREF(NSTRUCTREF+1,3*NATOMS))
   COORDSREF(1:NSTRUCTREF,1:3*NATOMS)=TEMPS2(1:NSTRUCTREF,1:3*NATOMS)
   COORDSREF(1+NSTRUCTREF,1:3*NATOMS)=PATHFRAMES(IMAX,1:3*NATOMS)
   DEALLOCATE(TEMPS2)
   
   ALLOCATE(TEMPS2(NPATH,NSTRUCTREF))
   TEMPS2(1:NPATH,1:NSTRUCTREF)=DISTFRAME(1:NPATH,1:NSTRUCTREF)
   DEALLOCATE(DISTFRAME)
   ALLOCATE(DISTFRAME(NPATH,NSTRUCTREF+1))
   DISTFRAME(1:NPATH,1:NSTRUCTREF)=TEMPS2(1:NPATH,1:NSTRUCTREF)
   DEALLOCATE(TEMPS2)

   NSTRUCTREF=NSTRUCTREF+1

   GOTO 765

ENDIF

PRINT '(A,I10,A)','mcpath> ',NSTRUCTREF,' reference structures:'
DO J1=1,NSTRUCTREF
   IF (J1.LE.NMIN+NTS) THEN
      IF (MOD(J1,2).EQ.1) THEN
         PRINT '(A,I10,A,G20.10)','minimum           at frame ',REFSTRUCTFRAME(J1),' energy=',VREF(J1)
      ELSE
         PRINT '(A,I10,A,G20.10)','transition  state at frame ',REFSTRUCTFRAME(J1),' energy=',VREF(J1)
      ENDIF
   ELSE
      PRINT '(A,I10,A,G20.10)','intervening state at frame ',REFSTRUCTFRAME(J1),' energy=',VREF(J1)
   ENDIF
ENDDO
!
! The number of extra reference structures is now fixed.
! However, we could reoptimise the bias exponent for different
! ts regions.
! Visualise new bias potential along the reaction pathway.
!
OPEN(LUNIT,FILE='Bias.new',STATUS='UNKNOWN')
DO J1=1,NPATH 
   DUMMY=0.0D0
   DUMMY2=0.0D0
   DMIN=1.0D100
   DO J2=1,NSTRUCTREF
      CALL ALIGN_DECIDE(COORDSREF(J2,1:3*NATOMS),PATHFRAMES(J1,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DSTRUCT(J2),DIST2,RIGIDBODY,RMAT)
      IF (DSTRUCT(J2).LT.DMIN) DMIN=DSTRUCT(J2)
     
   ENDDO
   DO J2=1,NSTRUCTREF
      DUMMY=DUMMY+VREF(J2)*EXP(BIASFAC/MAX(1.0D-2,DSTRUCT(J2))-BIASFAC/MAX(1.0D-2,DMIN))
      DUMMY2=DUMMY2+       EXP(BIASFAC/MAX(1.0D-2,DSTRUCT(J2))-BIASFAC/MAX(1.0D-2,DMIN))
!     PRINT '(A,2I6,5G15.5)','J1,J2, VREF,DSTRUCT,DUMMY,DUMMY2,DMIN=',J1,J2,VREF(J2),DSTRUCT(J2),DUMMY,DUMMY2,DMIN
   ENDDO

   DO K1=1,3
      DO K2=1,NATOMS
         QTEMP(K1,K2)=PATHFRAMES(J1,3*(K2-1)+K1)
      ENDDO
   ENDDO
   CALL GETORDER(QTEMP,QORDER,NATOMS)

   WRITE(LUNIT,'(I10,4G20.10)') J1,EOFS(J1),DUMMY/MAX(1.0D-10,DUMMY2),EOFS(J1)-DUMMY/MAX(1.0D-10,DUMMY2),QORDER

ENDDO
CLOSE(LUNIT)

! STOP !!! debug DJW
987 CONTINUE ! jump here if no bias
!
! NUPDATE specifies the interval for dynamically altering the maximum step size.
!
NUPDATE=100

STARTTS=MCPATHTS
ENDTS=MCPATHTS
IF (MCPATHTS.LE.0) STARTTS=1
IF (MCPATHTS.LE.0) ENDTS=NTS
IF (MCPATHTS.EQ.0) PRINT '(A)','mcpath> Sampling all transition states sequentially'
IF (MCPATHTS.LT.0) PRINT '(A)','mcpath> Sampling the whole path in one go'
! ENDTS=MIN(ENDTS,2) ! DJW debug
! NTS=MIN(NTS,2) ! DJW debug
! PRINT '(A)',' debug - sampling first two ts only'
DOTS=STARTTS
PRINT *,'STARTTS,ENDTS,NTS,DOTS=',STARTTS,ENDTS,NTS,DOTS
753 CONTINUE ! jump back here if we are doing multiple TS for MCPATHTS=0
STARTSAMPLE=1
ENDSAMPLE=2*NTS+1 ! last minimum
IF (MCPATHTS.GE.0) THEN
   PRINT '(A,I6)','mcpath> Sampling pathway for transition state ',DOTS
   IF (DOTS.EQ.1) THEN
      STARTSAMPLE=1
      IF (NTS.EQ.1) THEN
         ENDSAMPLE=3
         PRINT '(A,I6,A,I6)','mcpath> Sampling region defined by minimum ',1,' and minimum 2'
      ELSE
         ENDSAMPLE=4
         PRINT '(A,I6,A,I6)','mcpath> Sampling region defined by minimum ',1,' and transition state ',2
      ENDIF
   ELSEIF (DOTS.EQ.NTS) THEN
      STARTSAMPLE=2*NTS-2
      ENDSAMPLE=2*NTS+1
      PRINT '(A,I6,A,I6)','mcpath> Sampling region defined by transition state ',DOTS-1,' and minimum',2*DOTS+1
   ELSE
      STARTSAMPLE=2*DOTS-2
      ENDSAMPLE=2*DOTS+2
      PRINT '(A,I6,A,I6)','mcpath> Sampling region defined by transition states ',DOTS-1,' and ',DOTS+1
   ENDIF
   PRINT '(A,I6,A,I6,A,2I6)','        which are reference structures    ',STARTSAMPLE,' and ',ENDSAMPLE,' frames ', &
  &                     REFSTRUCTFRAME(STARTSAMPLE),REFSTRUCTFRAME(ENDSAMPLE)
ELSE
   PRINT '(A,I6,A,I6)','mcpath> Sampling region defined by minima ',1,' and ',NMIN
   PRINT '(A,I6,A,I6)','        which are reference structures    ',1,' and ',NPATH
ENDIF
PRINT '(A,I10,A)','mcpath> Step size will be adjusted every ',NUPDATE ,' MC steps during equilibration phase'
LSTRUCTREF=3
IF (NSTRUCTREF.GT.NMIN+NTS) THEN
   IF (MCPATHTS.GE.0) THEN
      DO J1=1,NSTRUCTREF
         IF ((REFSTRUCTFRAME(J1).LT.REFSTRUCTFRAME(ENDSAMPLE)).AND.(REFSTRUCTFRAME(J1).GT.REFSTRUCTFRAME(STARTSAMPLE))) THEN
            LSTRUCTREF=LSTRUCTREF+1
         ENDIF
      ENDDO
   ELSE
      LSTRUCTREF=NSTRUCTREF
   ENDIF
ENDIF

ALLOCATE(LCOORDSREF(LSTRUCTREF,3*NATOMS),LVREF(LSTRUCTREF),LREFSTRUCTFRAME(LSTRUCTREF))
IF (MCPATHTS.GE.0) THEN
   LCOORDSREF(1,1:3*NATOMS)=COORDSREF(STARTSAMPLE,1:3*NATOMS)
   CALL POTENTIAL(LCOORDSREF(1,1:3*NATOMS),DUMMY,GRAD,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
   DO K1=1,3
      DO K2=1,NATOMS
         QTEMP(K1,K2)=LCOORDSREF(1,3*(K2-1)+K1)
      ENDDO
   ENDDO
   CALL GETORDER(QTEMP,QORDER,NATOMS)
   PRINT '(3(A,G20.10),A,I6)','mcpath> Local reference 1 energy=',DUMMY,' should be ',VREF(STARTSAMPLE),' Q=',QORDER, &
  &                        ' frame=',REFSTRUCTFRAME(STARTSAMPLE)
   LVREF(1)=VREF(STARTSAMPLE)
   LREFSTRUCTFRAME(1)=REFSTRUCTFRAME(STARTSAMPLE)
   LCOORDSREF(2,1:3*NATOMS)=COORDSREF(2*DOTS,1:3*NATOMS)
   CALL POTENTIAL(LCOORDSREF(2,1:3*NATOMS),DUMMY,GRAD,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
   LVREF(2)=VREF(2*DOTS)
   LREFSTRUCTFRAME(2)=REFSTRUCTFRAME(2*DOTS)
   DO K1=1,3
      DO K2=1,NATOMS
         QTEMP(K1,K2)=LCOORDSREF(2,3*(K2-1)+K1)
      ENDDO
   ENDDO
   CALL GETORDER(QTEMP,QORDER,NATOMS)
   PRINT '(3(A,G20.10),A,I6)','mcpath> Local reference 2 energy=',DUMMY,' should be ',VREF(2*DOTS),' Q=',QORDER, &
  &                        ' frame=',REFSTRUCTFRAME(2*DOTS)
   LCOORDSREF(3,1:3*NATOMS)=COORDSREF(ENDSAMPLE,1:3*NATOMS)
   CALL POTENTIAL(LCOORDSREF(3,1:3*NATOMS),DUMMY,GRAD,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
   LVREF(3)=VREF(ENDSAMPLE)
   LREFSTRUCTFRAME(3)=REFSTRUCTFRAME(ENDSAMPLE)
   DO K1=1,3
      DO K2=1,NATOMS
         QTEMP(K1,K2)=LCOORDSREF(3,3*(K2-1)+K1)
      ENDDO
   ENDDO
   CALL GETORDER(QTEMP,QORDER,NATOMS)
   PRINT '(3(A,G20.10),A,I6)','mcpath> Local reference 3 energy=',DUMMY,' should be ',VREF(ENDSAMPLE),' Q=',QORDER, &
  &                        ' frame=',REFSTRUCTFRAME(ENDSAMPLE)
         CALL ALIGN_DECIDE(LCOORDSREF(1,1:3*NATOMS),LCOORDSREF(2,1:3*NATOMS),NATOMS,DEBUG, &
  &                  PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
   PRINT '(A,G20.10)','minimum distance between reference structures 1 and 2=',DUMMY
      CALL ALIGN_DECIDE(LCOORDSREF(3,1:3*NATOMS),LCOORDSREF(2,1:3*NATOMS),NATOMS,DEBUG, &
  &                  PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
   PRINT '(A,G20.10)','minimum distance between reference structures 3 and 2=',DUMMY
   LSTRUCTREF=3
   LNPATH=REFSTRUCTFRAME(ENDSAMPLE)-REFSTRUCTFRAME(STARTSAMPLE)+1
   PRINT '(A,I6)','mcpath> Setting local number of path frames to ',LNPATH
   IF (NSTRUCTREF.GT.NMIN+NTS) THEN
      DO J1=1,NSTRUCTREF
         IF (J1.EQ.2*DOTS) CYCLE ! don't add the transition state again!
         IF ((REFSTRUCTFRAME(J1).LT.REFSTRUCTFRAME(ENDSAMPLE)).AND.(REFSTRUCTFRAME(J1).GT.REFSTRUCTFRAME(STARTSAMPLE))) THEN
            LSTRUCTREF=LSTRUCTREF+1
            LCOORDSREF(LSTRUCTREF,1:3*NATOMS)=COORDSREF(J1,1:3*NATOMS)
            LVREF(LSTRUCTREF)=VREF(J1)
            LREFSTRUCTFRAME(LSTRUCTREF)=REFSTRUCTFRAME(J1)
            DO K1=1,3
               DO K2=1,NATOMS
                  QTEMP(K1,K2)=LCOORDSREF(LSTRUCTREF,3*(K2-1)+K1)
               ENDDO
            ENDDO
            CALL GETORDER(QTEMP,QORDER,NATOMS)
            PRINT '(A,I6,A,I6,A,I6,A,G20.10)','mcpath> Added intermediate reference structure ',J1,' to list for ts ',DOTS, &
  &                          ' frame number=',REFSTRUCTFRAME(J1),' Q=',QORDER
         ENDIF
      ENDDO
   ENDIF
ELSE
   LCOORDSREF(1:LSTRUCTREF,1:3*NATOMS)=COORDSREF(1:LSTRUCTREF,1:3*NATOMS)
   LVREF(1:LSTRUCTREF)=VREF(1:LSTRUCTREF)
   LREFSTRUCTFRAME(1:LSTRUCTREF)=REFSTRUCTFRAME(1:LSTRUCTREF)
   LNPATH=NPATH
   PRINT '(A,I6)','mcpath> Setting local number of path frames to ',LNPATH
   DO J1=1,NSTRUCTREF
      DO K1=1,3
         DO K2=1,NATOMS
            QTEMP(K1,K2)=LCOORDSREF(J1,3*(K2-1)+K1)
         ENDDO
      ENDDO
      CALL GETORDER(QTEMP,QORDER,NATOMS)
      PRINT '(A,I6,A,I6,A,G20.10)','mcpath> Added intermediate reference structure ',J1, &
  &                    ' frame number=',REFSTRUCTFRAME(J1),' Q=',QORDER
   ENDDO
ENDIF
!
! Identify neighbouring reference structures in the global
! list if required by MCPATHSCHECK>0. MC steps will be rejected
! if these neighbouring reference structures have shorter distances.
! They could be transition states or added non-stationary points.
!
NCHECKSTRUCT=0
IF (MCPATHSCHECK.GT.0) THEN
   PRINT '(A)','mcpath> MPATHCHECK > 0 has not been adapted to overlappgin regions yet'
   STOP
!
! Add transition states and minima up to MCPATHSCHECK steps away from the two end
! stationary points bracketing the current ts.
! Unless MCPATHTS is < 0 we are
! Sampling region defined by reference structures    ',STARTSAMPLE,' and ',ENDSAMPLE
!
   NCHECKMINUS=0
   NCHECKPLUS=0
! would need to change this...
   NCHECKMINUS=MIN(MCPATHSCHECK,2*DOTS-2)
   NCHECKPLUS=MIN(MCPATHSCHECK,NTS+NMIN-2*DOTS-1)
   PRINT '(A,I6,A,I6,A)','mcpath> Using ',NCHECKMINUS,' and ',NCHECKPLUS, &
  &           ' stationary point reference structures before and after the region in question'
   NCHECKSTRUCT=NCHECKMINUS+NCHECKPLUS
   PRINT '(A,I6,A,I6,A)','mcpath> Will reject configurations closer to these ',NCHECKSTRUCT,' neighbouring structures'
   ALLOCATE(LCOORDSCHECK(NCHECKSTRUCT,3*NATOMS))
   NDUMMY=0
   DO J1=MAX(2*DOTS-1-NCHECKMINUS,1),2*DOTS-2
      NDUMMY=NDUMMY+1
      PRINT '(A,I6,A)','mcpath> Adding reference structure ',J1,' to check list on minus side'
      LCOORDSCHECK(NDUMMY,1:3*NATOMS)=COORDSREF(J1,1:3*NATOMS)
   ENDDO
   PRINT '(A,4I6)','DOTS,NCHECKSTRUCT,NTS+NMIN,min=',DOTS,NCHECKSTRUCT,NTS+NMIN,MIN(2*DOTS+1+NCHECKSTRUCT,NTS+NMIN)
   DO J1=2*DOTS+2,MIN(2*DOTS+1+NCHECKPLUS,NTS+NMIN)
      NDUMMY=NDUMMY+1
      PRINT '(A,I6,A)','mcpath> Adding reference structure ',J1,' to check list on plus side'
      LCOORDSCHECK(NDUMMY,1:3*NATOMS)=COORDSREF(J1,1:3*NATOMS)
   ENDDO
ENDIF
!
! Reoptimise BIASFAC for this particular ts and pathway
!
! IF (MCPATHTS.EQ.0) THEN
IF ((NTS.GT.1).AND.(MCPATHTS.GE.0)) THEN
!
! We need to create LDISTFRAME for use in BIASFAC reoptimisation
!
   ALLOCATE(LDISTFRAME(LNPATH,LSTRUCTREF))
   ALLOCATE(LEOFS(LNPATH))
   LEOFS(1:LNPATH)=EOFS(REFSTRUCTFRAME(STARTSAMPLE):REFSTRUCTFRAME(ENDSAMPLE))
   DO J1=1,LNPATH
      DO J2=1,LSTRUCTREF
         CALL ALIGN_DECIDE(LCOORDSREF(J2,1:3*NATOMS), &
  &                    PATHFRAMES(J1+REFSTRUCTFRAME(STARTSAMPLE)-1,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,LDISTFRAME(J1,J2),DIST2,RIGIDBODY,RMAT)
      ENDDO
   ENDDO
   BIASSTEP=1.0D0
   BIASSTEP=1.0D0
   CALL SUMSQ(BIASFAC,LDISTFRAME,LNPATH,LSTRUCTREF,SUM1,LVREF,LEOFS,MAXDEV,IMAX)
   BIASBEST=BIASFAC
   SUMBEST=SUM1
   NDUMMY=0
   DO WHILE (NDUMMY<51)
      NDUMMY=NDUMMY+1
      BIASFAC=BIASFAC+BIASSTEP
      CALL SUMSQ(BIASFAC,LDISTFRAME,LNPATH,LSTRUCTREF,SUM2,LVREF,LEOFS,MAXDEV,IMAX)
      IF (DEBUG) PRINT '(4(A,G20.10))','mcpath> For bias exponent ',BIASFAC,' sum=',SUM2,' best=',SUMBEST,' step=',BIASSTEP
      IF (SUM2.LT.SUMBEST) THEN
         BIASBEST=BIASFAC
         SUMBEST=SUM2
      ELSE IF (SUM2.LT.SUM1) THEN
         BIASSTEP=BIASSTEP*1.05D0
      ELSE
         BIASSTEP=-BIASSTEP/2.0D0
      ENDIF
      SUM1=SUM2
      IF (ABS(BIASSTEP).LT.1.0D-10) EXIT
   ENDDO
   PRINT '(2(A,G20.10))','mcpath> Smallest deviation=',SUMBEST,' for bias exponent=',BIASBEST
   BIASFAC=BIASBEST
   CALL SUMSQ(BIASFAC,LDISTFRAME,LNPATH,LSTRUCTREF,SUM2,LVREF,LEOFS,MAXDEV,IMAX)
   PRINT '(A,G20.10,A,I10)','mcpath> Largest squared deviation=',MAXDEV,' for frame ',IMAX
!
! Plot local part of the bias function with reoptimised exponent
!
   WRITE(FNAME,'(I6)') DOTS
   FNAME='Bias.new.ts' // TRIM(ADJUSTL(FNAME))
   OPEN(LUNIT,FILE=FNAME,STATUS='UNKNOWN')
   PRINT '(A,I6)','DOTS=',DOTS
   PRINT '(A,I6)','REFSTRUCTFRAME(STARTSAMPLE)=',REFSTRUCTFRAME(STARTSAMPLE)
   DO J1=1,LNPATH
      DUMMY=0.0D0
      DUMMY2=0.0D0
      DMIN=1.0D100
      DO J2=1,LSTRUCTREF
         CALL ALIGN_DECIDE(LCOORDSREF(J2,1:3*NATOMS), &
  &                    PATHFRAMES(J1+REFSTRUCTFRAME(STARTSAMPLE)-1,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DSTRUCT(J2),DIST2,RIGIDBODY,RMAT)
         IF (DSTRUCT(J2).LT.DMIN) DMIN=DSTRUCT(J2)
      ENDDO
      DO J2=1,LSTRUCTREF
         DUMMY=DUMMY+LVREF(J2)*EXP(BIASFAC/MAX(1.0D-2,DSTRUCT(J2))-BIASFAC/MAX(1.0D-2,DMIN))
         DUMMY2=DUMMY2+        EXP(BIASFAC/MAX(1.0D-2,DSTRUCT(J2))-BIASFAC/MAX(1.0D-2,DMIN))
!     PRINT '(A,2I6,5G15.5)','J1,J2,LVREF,DSTRUCT,DUMMY,DUMMY2,DMIN=',J1,J2,LVREF(J2),DSTRUCT(J2),DUMMY,DUMMY2,DMIN
      ENDDO

      DO K1=1,3
         DO K2=1,NATOMS
            QTEMP(K1,K2)=PATHFRAMES(J1+REFSTRUCTFRAME(STARTSAMPLE)-1,3*(K2-1)+K1)
         ENDDO
      ENDDO
      CALL GETORDER(QTEMP,QORDER,NATOMS)

      WRITE(LUNIT,'(I10,4G20.10)') J1,LEOFS(J1),DUMMY/MAX(1.0D-10,DUMMY2),LEOFS(J1)-DUMMY/MAX(1.0D-10,DUMMY2),QORDER
   
   ENDDO
   CLOSE(LUNIT)
ENDIF ! end of condition for replotting with multiple ts

!
! Start from the specified transition state.
!
IF (MCPATHTS.GE.0) THEN
   IF (MCPATHSTART.EQ.0) THEN
      MCCOORDS(1:3*NATOMS)=COORDSREF(2*DOTS,1:3*NATOMS)
      PRINT '(A,I6)','mcpath> Starting from transition state at frame ',REFSTRUCTFRAME(2*DOTS)
   ELSEIF (MCPATHSTART.LT.0) THEN
      MCCOORDS(1:3*NATOMS)=COORDSREF(2*DOTS-1,1:3*NATOMS)
      PRINT '(A,I6)','mcpath> Starting from previous minimum at frame ',REFSTRUCTFRAME(2*DOTS-1)
   ELSE
      MCCOORDS(1:3*NATOMS)=COORDSREF(2*DOTS+1,1:3*NATOMS)
      PRINT '(A,I6)','mcpath> Starting from next minimum at frame ',REFSTRUCTFRAME(2*DOTS+1)
   ENDIF
ELSE
   IF (MCPATHSTART.EQ.0) THEN
      MCCOORDS(1:3*NATOMS)=COORDSREF(2,1:3*NATOMS)
      PRINT '(A,I6)','mcpath> Starting from transition state at frame ',REFSTRUCTFRAME(2)
   ELSEIF (MCPATHSTART.LT.0) THEN
      MCCOORDS(1:3*NATOMS)=COORDSREF(1,1:3*NATOMS)
      PRINT '(A,I6)','mcpath> Starting from min1 at frame ',REFSTRUCTFRAME(1)
   ELSE
      MCCOORDS(1:3*NATOMS)=COORDSREF(NSTRUCTREF,1:3*NATOMS)
      PRINT '(A,I6)','mcpath> Starting from final min at frame ',REFSTRUCTFRAME(NSTRUCTREF)
   ENDIF
ENDIF

CALL POTENTIAL(MCCOORDS,VOLD,GRAD,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
PRINT '(A,G20.10)','mcpath> Initial configuration energy             = ',VOLD

IF (MCBIAST) THEN
   DUMMY=0.0D0
   DUMMY2=0.0D0
   DMIN=1.0D100
   DO J2=1,LSTRUCTREF
      CALL ALIGN_DECIDE(LCOORDSREF(J2,1:3*NATOMS),MCCOORDS,NATOMS,DEBUG, &
  &                  PARAM1,PARAM2,PARAM3,BULKT,TWOD,DSTRUCT(J2),DIST2,RIGIDBODY,RMAT)
      IF (DSTRUCT(J2).LT.DMIN) DMIN=DSTRUCT(J2)
   ENDDO
   DO J2=1,LSTRUCTREF
      DUMMY=DUMMY+LVREF(J2)*EXP(BIASFAC/MAX(1.0D-2,DSTRUCT(J2))-BIASFAC/MAX(1.0D-2,DMIN))
      DUMMY2=DUMMY2+        EXP(BIASFAC/MAX(1.0D-2,DSTRUCT(J2))-BIASFAC/MAX(1.0D-2,DMIN))
!     PRINT '(A,I6,3G20.10)','J2,DSTRUCT,DUMMY,DUMMY2=',J2,DSTRUCT(J2),DUMMY,DUMMY2
!     PRINT '(A,2G20.10)','MAX(1.0D-2,DSTRUCT(J2)),MAX(1.0D-2,DMIN)=',MAX(1.0D-2,DSTRUCT(J2)),MAX(1.0D-2,DMIN)
   ENDDO
   WOLD=-DUMMY/MAX(1.0D-10,DUMMY2)
   PRINT '(A,G20.10)','mcpath> Initial bias energy                      = ',WOLD
   VOLD=VOLD+WOLD
ELSE
   WOLD=0.0D0
   WNEW=0.0D0
ENDIF

VNEW=VOLD
WNEW=WOLD
DO J1=1,NATOMS
   X(J1)=MCCOORDS(3*(J1-1)+1)
   Y(J1)=MCCOORDS(3*(J1-1)+2)
   Z(J1)=MCCOORDS(3*(J1-1)+3)
ENDDO
WRITE(*,'(A)') 'mcpath> Using hardcoded value (1) as random number seed'
CALL SDPRND(1)

NTOT=0
IACCEPT=0
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
         PRINT '(A,G20.10,A,G20.10,A,I10)','mcpath> ERROR *** energy for coordinates in MCCOORDSO=', &
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
! If not, reject and recount.
!
  RECOUNT=.TRUE.
  DISTMIN=1.0D100
  IMAGEMIN=1
  IMAGEMINF=1
  DISTMAX=-1.0D100
  IMAGEMAX=1
  DO J1=1,LSTRUCTREF
     CALL ALIGN_DECIDE(MCCOORDS,LCOORDSREF(J1,1:3*NATOMS),NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT, &
  &       TWOD,DUMMY3,DIST2,RIGIDBODY,RMAT)
     IF (DEBUG) PRINT '(A,I10,A,G20.10)','mcpath> distance for local reference structure ',J1,' is ',DUMMY3
     IF (DUMMY3.LT.DISTMIN) THEN
        DISTMIN=DUMMY3
        IMAGEMIN=J1
     ENDIF
     IF (DUMMY3.GT.DISTMAX) THEN
        DISTMAX=DUMMY3
        IMAGEMAX=J1
     ENDIF
  ENDDO

  IF (IMCSTEP.EQ.1) THEN
     IMAGEMINO=IMAGEMIN ! initialise this information at first step
     DISTMINF=1.0D100
     IMAGEMINF=1
     DO J1=1,NPATH
        CALL ALIGN_DECIDE(MCCOORDS,PATHFRAMES(J1,1:3*NATOMS),NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT, &
     &       TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
        IF (DEBUG) PRINT '(A,I10,A,G20.10)','mcpath> distance for path frame ',J1,' is ',DUMMY
        IF (DUMMY.LT.DISTMINF) THEN
           DISTMINF=DUMMY
           IMAGEMINF=J1
        ENDIF
     ENDDO
     PRINT '(A,I6,A,G20.10)','mcpath> Initial configuration is closest to path frame ',IMAGEMINF,' distance=',DISTMINF
     IF ((MCPATHTS.GE.0).AND.((IMAGEMINF.LT.REFSTRUCTFRAME(STARTSAMPLE)).OR.(IMAGEMINF.GT.REFSTRUCTFRAME(ENDSAMPLE)))) THEN
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
           PRINT '(A)','mcpath> Too many failures - give up!'
           STOP
        ENDIF

        GOTO 961
     ENDIF
     IMAGEMINFO=IMAGEMINF
  ELSE
     IMAGEMINSTART=IMAGEMINFO
     NDIRECTION=1
! 
! Get distance for previous best frame, initialised at IMAGEMINFO
!
     CALL ALIGN_DECIDE(MCCOORDS,PATHFRAMES(IMAGEMINSTART,1:3*NATOMS),NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT, &
  &    TWOD,DISTMINF,DIST2,RIGIDBODY,RMAT)
!    IF (DEBUG) PRINT '(A,I10,A,G20.10)','mcpath> distance for path frame ',IMAGEMINSTART,' is ',DISTMINF
263 CONTINUE
     IF ((IMAGEMINSTART+NDIRECTION.GT.0).AND.(IMAGEMINSTART+NDIRECTION.LE.NPATH)) THEN
        CALL ALIGN_DECIDE(MCCOORDS,PATHFRAMES(IMAGEMINSTART+NDIRECTION,1:3*NATOMS),NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT, &
  &                      TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
     ELSE
        DUMMY=1.0D100
     ENDIF
     IF (DEBUG) PRINT '(A,I10,A,G20.10)','mcpath> distance for path frame ',IMAGEMINSTART+NDIRECTION,' is ',DUMMY
     IF (DUMMY.LT.DISTMINF) THEN
        DISTMINF=DUMMY
        IMAGEMINSTART=IMAGEMINSTART+NDIRECTION
!
! Closest frame has changed. Need to go back and check again so that we identify the
! closest frame consistently and avoid getting stuck! It needs to be a local minimum at least!
!
        GOTO 263
     ENDIF
     IF ((IMAGEMINSTART-NDIRECTION.GT.0).AND.(IMAGEMINSTART+NDIRECTION.LE.NPATH)) THEN
        CALL ALIGN_DECIDE(MCCOORDS,PATHFRAMES(IMAGEMINSTART-NDIRECTION,1:3*NATOMS),NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT, &
  &                      TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
          ELSE
        DUMMY=1.0D100
     ENDIF
     IF (DEBUG) PRINT '(A,I10,A,G20.10)','mcpath> distance for path frame ',IMAGEMINSTART-NDIRECTION,' is ',DUMMY
     IF (DUMMY.LT.DISTMINF) THEN
        DISTMINF=DUMMY
        IMAGEMINSTART=IMAGEMINSTART-NDIRECTION
        NDIRECTION=-NDIRECTION
!
! Closest frame has changed. Need to go back and check again so that we identify the
! closest frame consistently and avoid getting stuck!
!
        GOTO 263
     ENDIF
     IMAGEMINF=IMAGEMINSTART
  ENDIF

  IF ((DISTMINF-DISTMIN.GT.GEOMDIFFTOL).AND.(LREFSTRUCTFRAME(IMAGEMIN).NE.IMAGEMINF)) THEN
     PRINT '(A,I6,A,I6,A,G15.5,A,G15.5,A,I6)','mcpath> WARNING *** Distance to local reference structure ',IMAGEMIN,' frame ', &
  &              LREFSTRUCTFRAME(IMAGEMIN),' of ',DISTMIN,' is less than ',DISTMINF,' for current frame ',IMAGEMINF
!
! Try for a local minimum from this reference structure.
!
     IMAGEMINSTART=LREFSTRUCTFRAME(IMAGEMIN)
     NDIRECTION=1
     DISTMINF=DISTMIN
764 CONTINUE
     IF ((IMAGEMINSTART+NDIRECTION.GT.0).AND.(IMAGEMINSTART+NDIRECTION.LE.NPATH)) THEN
        CALL ALIGN_DECIDE(MCCOORDS,PATHFRAMES(IMAGEMINSTART+NDIRECTION,1:3*NATOMS),NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT, &
  &                      TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
     ELSE
        DUMMY=1.0D100
     ENDIF
     IF (DEBUG) PRINT '(A,I10,A,G20.10)','mcpath> distance for path frame ',IMAGEMINSTART+NDIRECTION,' is ',DUMMY
     IF (DUMMY.LT.DISTMINF) THEN
        DISTMINF=DUMMY
        IMAGEMINSTART=IMAGEMINSTART+NDIRECTION
!
! Closest frame has changed. Need to go back and check again so that we identify the
! closest frame consistently and avoid getting stuck! It needs to be a local minimum at least!
!
        GOTO 764
     ENDIF
     IF ((IMAGEMINSTART-NDIRECTION.GT.0).AND.(IMAGEMINSTART+NDIRECTION.LE.NPATH)) THEN
        CALL ALIGN_DECIDE(MCCOORDS,PATHFRAMES(IMAGEMINSTART-NDIRECTION,1:3*NATOMS),NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT, &
  &                      TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
          ELSE
        DUMMY=1.0D100
     ENDIF
     IF (DEBUG) PRINT '(A,I10,A,G20.10)','mcpath> distance for path frame ',IMAGEMINSTART-NDIRECTION,' is ',DUMMY
     IF (DUMMY.LT.DISTMINF) THEN
        DISTMINF=DUMMY
        IMAGEMINSTART=IMAGEMINSTART-NDIRECTION
        NDIRECTION=-NDIRECTION
!
! Closest frame has changed. Need to go back and check again so that we identify the
! closest frame consistently and avoid getting stuck!
!
        GOTO 764
     ENDIF
     IMAGEMINF=IMAGEMINSTART
     PRINT '(A,I6,A,G15.5)','mcpath> WARNING *** Resetting closest frame to ',IMAGEMINF,' distance ',DISTMINF

!    CALL ALIGN_DECIDE(MCCOORDS,PATHFRAMES(LREFSTRUCTFRAME(IMAGEMIN),1:3*NATOMS),NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT, &
! &    TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
!    IF (LREFSTRUCTFRAME(IMAGEMIN)-1.GT.0) THEN
!       CALL ALIGN_DECIDE(MCCOORDS,PATHFRAMES(LREFSTRUCTFRAME(IMAGEMIN)-1,1:3*NATOMS),NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT, &
! &                      TWOD,DUMMY2,DIST2,RIGIDBODY,RMAT)
!    ELSE
!       DUMMY2=1.0D100
!    ENDIF
!    IF (LREFSTRUCTFRAME(IMAGEMIN)+1.LE.NPATH) THEN
!       CALL ALIGN_DECIDE(MCCOORDS,PATHFRAMES(LREFSTRUCTFRAME(IMAGEMIN)+1,1:3*NATOMS),NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT, &
! &                      TWOD,DUMMY3,DIST2,RIGIDBODY,RMAT)
!    ELSE
!       DUMMY3=1.0D100
!    ENDIF
!    IF ((DUMMY.LT.DUMMY2).AND.(DUMMY.LT.DUMMY3)) THEN
!       PRINT '(A,I6,A)','mcpath> WARNING *** Resetting to frame ',LREFSTRUCTFRAME(IMAGEMIN),' for proposed step'
!       DISTMINF=DISTMIN
!       IMAGEMINF=LREFSTRUCTFRAME(IMAGEMIN)
!    ELSE
!       PRINT '(A)','mcpath> WARNING *** Cannot reset - proposed frame is not a local minimum'
!       PRINT '(A,3G20.10)','mcpath> distances: ',DUMMY2,DUMMY,DUMMY3
!    ENDIF
  ENDIF

  IF (DISTMIN.LT.MCPATHDMAX) RECOUNT=.FALSE.
  IF (DEBUG.AND.(DISTMIN.GT.MCPATHDMAX)) PRINT '(A,2G20.10,L5)','DISTMIN,MCPATHDMAX,RECOUNT=',DISTMIN,MCPATHDMAX,RECOUNT
  IF (IMCSTEP.GT.MCPATHEQUIL) THEN
     IF (RECOUNT) IREJDIST=IREJDIST+1
  ENDIF
  IF ((MCPATHTS.GE.0).AND.((IMAGEMINF.LT.REFSTRUCTFRAME(STARTSAMPLE)).OR.(IMAGEMINF.GT.REFSTRUCTFRAME(ENDSAMPLE)))) THEN
     IF (DEBUG) PRINT '(A)','mcpath> reject: closest frame is outside the range defined by the two neighbouring minima'
     RECOUNT=.TRUE.
     IF (IMCSTEP.GT.MCPATHEQUIL) IREJFRAME=IREJFRAME+1
  ENDIF
!
! Check rejection condition on neighbouring reference strucures.
!
  DO J1=1,NCHECKSTRUCT
     CALL ALIGN_DECIDE(MCCOORDS,LCOORDSCHECK(J1,1:3*NATOMS),NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT, &
  &       TWOD,DISTCHECK,DIST2,RIGIDBODY,RMAT)
     IF (DEBUG) PRINT '(A,I10,A,G20.10)','mcpath> distance for check reference structure ',J1,' is ',DISTCHECK
     IF (DISTCHECK.LT.DISTMIN) THEN ! reject step
        IREJCHECK=IREJCHECK+1
        RECOUNT=.TRUE.
        PRINT '(A)','mcpath> Rejecting step because configuration is nearer to check structure ',J1
        PRINT '(A,F15.5,A,F15.5)','mcpath> Distance is ',DISTCHECK,' compared to ',DISTMIN
        EXIT
     ENDIF
  ENDDO

  IF (IMCSTEP.GT.MCPATHEQUIL) THEN
     IF (DISTMIN.LT.XDISTMIN) XDISTMIN=DISTMIN
     IF (DISTMAX.GT.XDISTMAX) XDISTMAX=DISTMAX
  ENDIF
  IF (DEBUG.AND.(MOD(IMCSTEP-1,MCPATHPRTFRQ).EQ.0)) THEN
     PRINT '(A,2G20.10,A,2I6)','mcpath> Minimum and maximum reference distances are ',DISTMIN,DISTMAX, &
  &                            ' for local structures ',IMAGEMIN,IMAGEMAX
     PRINT '(A,G20.10,A,I6)','mcpath> Minimum path frame distance is ',DISTMINF, &
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
!    PRINT '(A,G20.10)','mcpath> after potential call VNEW=',VNEW
     IF (MCBIAST) THEN
        DUMMY=0.0D0
        DUMMY2=0.0D0
        DMIN=1.0D100
        DO J2=1,LSTRUCTREF
           CALL ALIGN_DECIDE(LCOORDSREF(J2,1:3*NATOMS),MCCOORDS,NATOMS,DEBUG, &
  &                       PARAM1,PARAM2,PARAM3,BULKT,TWOD,DSTRUCT(J2),DIST2,RIGIDBODY,RMAT)
           IF (DSTRUCT(J2).LT.DMIN) DMIN=DSTRUCT(J2)
        ENDDO
        DO J2=1,LSTRUCTREF
           DUMMY=DUMMY+LVREF(J2)*EXP(BIASFAC/MAX(1.0D-2,DSTRUCT(J2))-BIASFAC/MAX(1.0D-2,DMIN))
           DUMMY2=DUMMY2+       EXP(BIASFAC/MAX(1.0D-2,DSTRUCT(J2))-BIASFAC/MAX(1.0D-2,DMIN))
!          PRINT '(A,I6,3G20.10)','J2,DSTRUCT,DUMMY,DUMMY2=',J2,DSTRUCT(J2),DUMMY,DUMMY2
!          PRINT '(A,2G20.10)','MAX(1.0D-2,DSTRUCT(J2)),MAX(1.0D-2,DMIN)=',MAX(1.0D-2,DSTRUCT(J2)),MAX(1.0D-2,DMIN)
        ENDDO
        WNEW=-DUMMY/MAX(1.0D-10,DUMMY2)
!       PRINT '(A,G20.10)','mcpath> after potential call WNEW=',WNEW
        VNEW=VNEW+WNEW
!       PRINT '(A,G20.10)','mcpath> after potential call VNEW changed to ',VNEW
     ENDIF
     WCOMP=(VNEW-VOLD)/MCPATHTEMP
     DUMMY=MIN(1.0D0,EXP(-WCOMP))
     RANDOM=DPRAND()
!    PRINT '(A,2G20.10)','RANDOM,DUMMY=',RANDOM,DUMMY
     IF (RANDOM.GT.DUMMY) RECOUNT=.TRUE. ! RECOUNT is initialised to .FALSE. at the top of the loop
  ENDIF
  IF ((IMCSTEP.EQ.1).AND.RECOUNT) THEN
     PRINT '(A)','mcpath> WARNING *** first step was initially rejected - accepting it - check initial step size '
     PRINT '(A,G20.10,A)','mcpath> Sqrt(3)*maximum step size is ',1.73D0*MCPATHSTEP,' versus 0.1 used in bias potential'
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
  IF (VNEWSAVE.GT.170.0D0) THEN
     PRINT '(A)',' check VNEWSAVE value !'
     STOP ! DJW debug
  ENDIF
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
     IF (MOD(IMCSTEP-1,MCPATHPRTFRQ).EQ.0) WRITE(*,'(A,G20.10,A,G20.10,A,I6,A)') 'mcpath> adjusting step-size> current step size = ',&
 &                   MCPATHSTEP,' acceptance ratio = ', WAC ,' over ', NUPDATE, ' steps'
  ENDIF ! step size update
  IF (IMCSTEP.EQ.MCPATHEQUIL) THEN
     WRITE(*, '(A)') 'mcpath> ---------- Equilibration done '
     WRITE(*, '(A,F20.10,A,I10,A,G20.10)') 'mcpath> Temperature = ', MCPATHTEMP, &
  &                 ' MCSteps = ', IMCSTEP,' MarkovEner = ', VNEW
     WRITE(*, '(A,G20.10,A,G20.10,A,I6,A)') 'mcpath> Final step size = ',MCPATHSTEP, &
  &'  corresponding to acceptance ratio = ', WAC ,' over previous ', NUPDATE, ' steps'
     WRITE(*, '(A,G20.10)') 'mcpath>   compare with target acceptance ratio = ', MCPATHACCRATIO
     MEANBIASEQ=MEANBIASEQ/MCPATHEQUIL
     IF (MCBIAST) THEN
        WRITE(*,'(A,G20.10)') 'mcpath> <W> over equilibration phase=',MEANBIASEQ
     ENDIF
     WRITE(*, '(A)') 'mcpath> ---------- Starting production run '
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
   WRITE(FNAME,'(I6)') DOTS
   WRITE(FNAME2,'(I6)') INT((IMCSTEP-MCPATHEQUIL)*1.0D0/1.0D6)
   FNAME='mcpath.Q.ts' // TRIM(ADJUSTL(FNAME)) // '.' // TRIM(ADJUSTL(FNAME2))
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
   PRINT '(A)','mcpath> Statistics for closest reference structure:'
   PRINT '(2I10,F20.5,A3)',(J1,NEAREST(J1),1.0D2*NEAREST(J1)/(IMCSTEP-MCPATHEQUIL),' % ',J1=1,LSTRUCTREF)

   LUNIT=GETUNIT()
   WRITE(FNAME,'(I6)') DOTS
   WRITE(FNAME2,'(I6)') INT((IMCSTEP-MCPATHEQUIL)*1.0D0/1.0D6)
   FNAME='mcpath.f.ts' // TRIM(ADJUSTL(FNAME)) // '.' // TRIM(ADJUSTL(FNAME2))
   OPEN(LUNIT,FILE=FNAME,STATUS='UNKNOWN')
   DO J1=1,NPATH
      IF (MCBIAST) THEN
         WRITE(LUNIT,'(I6,6G20.10)') J1,SLENGTH(J1),EOFS(J1),QFRAME(J1),NEARESTFW(J1)/MEANBIAS, &
  &          -LOG(MAX(1.0D-100,NEARESTFW(J1)/MEANBIAS)),1.0D2*NEARESTF(J1)/(IMCSTEP-MCPATHEQUIL)
      ELSE
         WRITE(LUNIT,'(I6,6G20.10)') J1,SLENGTH(J1),EOFS(J1),QFRAME(J1),NEARESTFW(J1)/(IMCSTEP-MCPATHEQUIL), &
  &          -LOG(MAX(1.0D-100,NEARESTFW(J1)/(IMCSTEP-MCPATHEQUIL))),1.0D2*NEARESTF(J1)/(IMCSTEP-MCPATHEQUIL)
      ENDIF
   ENDDO
   CLOSE(LUNIT)
ENDIF

ENDDO ! end of MC loop

WRITE (*,'(A)') 'mcpath> Exited main MC loop. '

IF (.TRUE.) THEN
   LUNIT=GETUNIT()
   WRITE(FNAME,'(I6)') DOTS
   FNAME='mcpath.Q.ts' // TRIM(ADJUSTL(FNAME))
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
   WRITE(FNAME,'(I6)') DOTS
   FNAME='mcpath.f.ts' // TRIM(ADJUSTL(FNAME)) 
   OPEN(LUNIT,FILE=FNAME,STATUS='UNKNOWN')
   DO J1=1,NPATH
      IF (MCBIAST) THEN
         WRITE(LUNIT,'(I6,6G20.10)') J1,SLENGTH(J1),EOFS(J1),QFRAME(J1),NEARESTFW(J1)/MEANBIAS, &
  &          -LOG(MAX(1.0D-100,NEARESTFW(J1)/MEANBIAS)),1.0D2*NEARESTF(J1)/MCPATHSTEPS
      ELSE
         WRITE(LUNIT,'(I6,6G20.10)') J1,SLENGTH(J1),EOFS(J1),QFRAME(J1),NEARESTFW(J1)/MCPATHSTEPS, &
  &          -LOG(MAX(1.0D-100,NEARESTFW(J1)/MCPATHSTEPS)),1.0D2*NEARESTF(J1)/MCPATHSTEPS

      ENDIF
   ENDDO
   CLOSE(LUNIT)

ENDIF
!
! Printing summary
!
WRITE(*,'(A,I10,A,I10,A,F15.5,A)') 'mcpath> ',IACCEPT, ' steps accepted out of ', &
   &      MCPATHSTEPS+MCPATHEQUIL, ' i.e. ',IACCEPT*100.0D0/(MCPATHSTEPS+MCPATHEQUIL),'%'
WRITE(*,'(A,G20.10)') 'bspt> Final stepsize ',MCPATHSTEP
PRINT '(A)','mcpath> Number of production run steps for which configuration was closest to a given reference structure:'
PRINT '(2I10,F20.5,A3)',(J1,NEAREST(J1),1.0D2*NEAREST(J1)/MCPATHSTEPS,' % ',J1=1,LSTRUCTREF)
PRINT '(A)','mcpath> Number of production run steps for which configuration was closest to a given path frame:'
PRINT '(2I10,F20.5,A3)',(J1,NEARESTF(J1),1.0D2*NEARESTF(J1)/MCPATHSTEPS,' % ', &
  &             J1=1,REFSTRUCTFRAME(STARTSAMPLE),REFSTRUCTFRAME(ENDSAMPLE))
NDUMMY=0
DO J1=REFSTRUCTFRAME(STARTSAMPLE),REFSTRUCTFRAME(ENDSAMPLE)
   IF (NEARESTF(J1).EQ.0) NDUMMY=NDUMMY+1
ENDDO
IF (NDUMMY.GT.0) PRINT '(A,I10,A)','mcpath> WARNING *** ',NDUMMY,' frames were not visited'
PRINT '(A,2G20.10)','mcpath> Smallest and largest distances from references seen were ',XDISTMIN,XDISTMAX
PRINT '(A,I10,A,F10.2,A)','mcpath> Production steps rejected on distance criterion=',IREJDIST,' i.e. ', &
  &                        1.0D2*IREJDIST/MCPATHSTEPS,'% '
PRINT '(A,I10,A,F10.2,A)','mcpath> Production steps rejected on closest path frame criterion=',IREJFRAME,' i.e. ', &
  &                        1.0D2*IREJFRAME/MCPATHSTEPS,'% '
IF (NCHECKSTRUCT.GT.0) THEN
   PRINT '(A,I10,A,F10.2,A)','mcpath> Production steps rejected on distance criterion for neighbouring reference structures=', &
  &      IREJCHECK,' i.e. ',1.0D2*IREJCHECK/MCPATHSTEPS,'% '
ENDIF
PRINT '(A,3G20.10)','mcpath> <Q>, <Q^2> and sigma=',QORDERMEAN,QORDERSQMEAN,SQRT(QORDERSQMEAN-QORDERMEAN**2)
PRINT '(A,3G20.10,A)','mcpath> <Q>, <Q^2> and sigma=',QORDERMEAN2,QORDERSQMEAN2,SQRT(QORDERSQMEAN2-QORDERMEAN2**2),' from Q bins'

PRINT *,'DOTS,ENDTS=',DOTS,ENDTS
IF ((MCPATHTS.GE.0).AND.(DOTS.LT.ENDTS)) THEN
   DOTS=DOTS+1
   DEALLOCATE(LCOORDSREF,LVREF,LDISTFRAME,LEOFS,LREFSTRUCTFRAME)
   GOTO 753
ENDIF

DEALLOCATE(COORDSREF,VREF,NEAREST,DSTRUCT)
DEALLOCATE(NEARESTF,NEARESTFW)
DEALLOCATE(QORDERHIST,QORDERVISITS)
IF (ALLOCATED(DISTFRAME)) DEALLOCATE(DISTFRAME)
DEALLOCATE(PATHFRAMES)

963 CONTINUE

! NTS=MIN(2,NTS) ! DJW debug
! PRINT '(A)',' debug - setting two ts only for wham'
PATHSTR='f.ts'
LINVART=.FALSE.

CALL SDPRND(0) ! initialise
FIXZ=.FALSE.
IF (MCPATHTS.LT.0) THEN
   ALLOCATE(ZNORM(1))
   DO J1=1,NPATH
      IF (LINVART) THEN ! positive probabilities!
         PROBS(J1)=(DPRAND()+1.0D0)/2.0D0
      ELSE
         PROBS(J1)=DPRAND()
      ENDIF
   ENDDO
   IF (LINVART) ZNORM(1)=1.0D0
   IF (.NOT.LINVART) ZNORM(1)=0.0D0
   CALL MYWHAM(1,NPATH,PATHSTR,PROBS,'          ',LINVART,OFAILS,1,1,ZNORM,FIXZ,ALLZEROS)
ELSE
   ALLOCATE(ZNORM(NTS))
   DO J1=1,NPATH
      IF (LINVART) THEN ! positive probabilities!
         PROBS(J1)=(DPRAND()+1.0D0)/2.0D0
      ELSE
         PROBS(J1)=DPRAND()
      ENDIF
   ENDDO
   DO J1=1,NTS
      ZNORM(J1)=DPRAND()
   ENDDO
   IF (LINVART) ZNORM(1)=1.0D0
   IF (.NOT.LINVART) ZNORM(1)=0.0D0
   CALL MYWHAM(NTS,NPATH,PATHSTR,PROBS,'          ',LINVART,OFAILS,1,NTS,ZNORM,FIXZ,ALLZEROS)
ENDIF

IF (.NOT.LINVART) THEN
   DO J1=1,NPATH
      PROBS(J1)=EXP(PROBS(J1)-PROBS(1))
   ENDDO
ENDIF

IF (OFAILS) THEN
   PRINT '(A)','mcpath> Pofs and PofsQ will not be produced'
   GOTO 864
ENDIF
DUMMY=-1.0D100
PRINT '(A)','mcpath> frames, s, ln(P(s)),P(s) from WHAM:'
DO J1=1,NPATH
   PRINT '(I6,3G20.10)',J1,EOFS(J1),LOG(PROBS(J1)),PROBS(J1)
   IF (PROBS(J1).GT.DUMMY) DUMMY=PROBS(J1)
ENDDO
PRINT '(A)','mcpath> frames, s, P(s) rescaled:'
DO J1=1,NPATH
   PROBS(J1)=PROBS(J1)/DUMMY
   PRINT '(I6,2G20.10)',J1,EOFS(J1),PROBS(J1)
ENDDO
!
! GWIDTH is sigma^2, the variance
!
GWIDTHS=SLENGTH(NPATH)*1.0D0/NPATH
PRINT '(A,G20.10)','mcpath> standard deviation sigma for P(s) bins smoothing is ',SQRT(GWIDTHS)
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
LUNIT=GETUNIT()
OPEN(LUNIT,FILE='Pofs',STATUS='UNKNOWN')
DUMMY=0.0D0
DO J1=1,NPATH
   WRITE(LUNIT,'(I6,5G20.10)') J1,SLENGTH(J1),MAX(1.0D-50,SMOOTHS(J1)),PROBS(J1),QFRAME(J1),EOFS(J1)
   DUMMY=DUMMY+SMOOTHS(J1)
ENDDO
CLOSE(LUNIT)
PRINT '(A,G20.10)','mcpath2> P(s) normalisation=',DUMMY

PROBQFROMS(1:MCPATHBINS)=0.0D0
DO J1=1,NPATH
   IBININDEX=INT((QFRAME(J1)-MIN(MCPATHQMAX,MCPATHQMIN))/QORDERHISTINT)+1
   PROBQFROMS(IBININDEX)=PROBQFROMS(IBININDEX)+QFRAME(J1)*PROBS(J1)
ENDDO

864 CONTINUE ! jump here if P(s) overlap failed and OFAILS is true
PATHSTR='Q.ts'
IF (ALLOCATED(ZNORM)) DEALLOCATE(ZNORM)
ALLOCATE(ZNORM(NTS))
CALL MYWHAM(NTS,MCPATHBINS,PATHSTR,PROBQ,'          ',LINVART,OFAILQ,1,NTS,ZNORM,FIXZ,ALLZEROQ)
IF (OFAILQ) THEN
   PRINT '(A)','mcpath> PofQ and PofsQ will not be produced'
   GOTO 862
ENDIF
IF (.NOT.LINVART) THEN
   DO J1=1,MCPATHBINS
      PROBQ(J1)=EXP(PROBQ(J1)-PROBQ(1))
   ENDDO
ENDIF

DUMMY=-1.0D100
DO J1=1,MCPATHBINS
   IF (PROBQ(J1).GT.DUMMY) DUMMY=PROBQ(J1)
ENDDO
DO J1=1,MCPATHBINS
   PROBQ(J1)=PROBQ(J1)/DUMMY
ENDDO
GWIDTHQ=QORDERHISTINT*0.05D0
SMOOTHQ(1:MCPATHBINS)=0.0D0
SMOOTHQFROMS(1:MCPATHBINS)=0.0D0
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
ELSE
   DUMMY=0.0D0
   DO J1=1,MCPATHBINS
      DUMMY=DUMMY+PROBQFROMS(J1)
   ENDDO
   PROBQFROMS(1:MCPATHBINS)=PROBQFROMS(1:MCPATHBINS)/DUMMY
   DO J1=1,MCPATHBINS
      DO J2=1,MCPATHBINS
         SMOOTHQFROMS(J1)=SMOOTHQFROMS(J1)+ &
  &                           PROBQFROMS(J2)*EXP(-((BINLABELQORDER(J1)-BINLABELQORDER(J2))**2/(2.0D0*GWIDTHQ)))/NORM2(J2)
      ENDDO
   ENDDO
ENDIF
LUNIT=GETUNIT()
OPEN(LUNIT,FILE='PofQ',STATUS='UNKNOWN')
DUMMY=0.0D0
DUMMY2=0.0D0
DO J1=1,MCPATHBINS
   WRITE(LUNIT,'(I6,5G20.10)') J1,BINLABELQORDER(J1),MAX(1.0D-50,SMOOTHQ(J1)), &
  &                            MAX(1.0D-50,SMOOTHQFROMS(J1)),PROBQ(J1),PROBQFROMS(J1)
   DUMMY=DUMMY+SMOOTHQ(J1)
   DUMMY2=DUMMY2+SMOOTHQFROMS(J1)
ENDDO
PRINT '(A,2G20.10)','mcpath2> P(Q) and P(Q) from P(s) normalisation=',DUMMY,DUMMY2
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

ALLOCATE(SMOOTH2D(NPATH,MCPATHBINS))
SMOOTH2D(1:NPATH,1:MCPATHBINS)=0.0D0
IF (.NOT.OFAILS) THEN
   DO J2=1,NPATH
      DO J3=1,MCPATHBINS
         DO J1=1,NPATH
            IBININDEX=INT((QFRAME(J1)-MIN(MCPATHQMAX,MCPATHQMIN))/QORDERHISTINT)+1
            SMOOTH2D(J2,J3)=SMOOTH2D(J2,J3)+PROBS(J1)*EXP(-(SLENGTH(J1)-SLENGTH(J2))**2/(2.0D0*GWIDTHS)) &
  &                                 *EXP(-((BINLABELQORDER(IBININDEX)-BINLABELQORDER(J3))**2/(2.0D0*GWIDTHQ)))/NORM1(J1)
         ENDDO
      ENDDO
   ENDDO

   LUNIT=GETUNIT()
   OPEN(LUNIT,FILE='PofsQ',STATUS='UNKNOWN')
   DUMMY=0.0D0
   DO J1=1,NPATH
      DO J2=1,MCPATHBINS
         WRITE(LUNIT,'(2I6,3G20.10)') J1,J2,SLENGTH(J1),BINLABELQORDER(J2),MAX(1.0D-50,SMOOTH2D(J1,J2))
         DUMMY=DUMMY+SMOOTH2D(J1,J2)
      ENDDO
   ENDDO
   PRINT '(A,G20.10)','mcpath2> P(s,Q) normalisation=',DUMMY
   CLOSE(LUNIT)
ENDIF

862 CONTINUE ! jump here if P(Q) overlap failed and OFAILQ is true.

DEALLOCATE(BINLABELQORDER,SMOOTH2D,NORM1)
IF (ALLOCATED(ZNORM)) DEALLOCATE(ZNORM)

RETURN

END SUBROUTINE MCPATH

SUBROUTINE SUMSQ(BIASFAC,LDISTFRAME,LNPATH,LSTRUCTREF,SUM,LVREF,LEOFS,MAXDEV,IMAX)
USE KEY, ONLY : GEOMDIFFTOL, EDIFFTOL
IMPLICIT NONE
INTEGER J1,J2,LNPATH,LSTRUCTREF,IMAX
DOUBLE PRECISION BIASFAC, LDISTFRAME(LNPATH,LSTRUCTREF), DUMMY, DUMMY2, SUM, LVREF(LSTRUCTREF), &
  &              LEOFS(LNPATH), DMIN, MAXDEV, DUMMY3
LOGICAL SPTEST

SUM=0.0D0
MAXDEV=-1.0D0
IMAX=1
DO J1=1,LNPATH
   SPTEST=.FALSE.
   DUMMY=0.0D0
   DUMMY2=0.0D0
   DMIN=1.0D100
   DO J2=1,LSTRUCTREF
      IF (LDISTFRAME(J1,J2).LT.DMIN) DMIN=LDISTFRAME(J1,J2)
      IF ((LDISTFRAME(J1,J2).LT.GEOMDIFFTOL).AND.(ABS(LEOFS(J1)-LVREF(J2)).LT.EDIFFTOL)) SPTEST=.TRUE.
   ENDDO
   DO J2=1,LSTRUCTREF
      DUMMY=DUMMY+LVREF(J2)*EXP(BIASFAC/MAX(1.0D-2,LDISTFRAME(J1,J2))-BIASFAC/MAX(1.0D-2,DMIN))
      DUMMY2=DUMMY2+        EXP(BIASFAC/MAX(1.0D-2,LDISTFRAME(J1,J2))-BIASFAC/MAX(1.0D-2,DMIN))
   ENDDO
   DUMMY3=(LEOFS(J1)-DUMMY/MAX(1.0D-10,DUMMY2))**2
   SUM=SUM+DUMMY3
   IF (SPTEST) SUM=SUM+10.0D0*DUMMY3 ! increase the weight for stationary points
   IF (DUMMY3.GT.MAXDEV) THEN
      MAXDEV=DUMMY3
      IMAX=J1
   ENDIF
ENDDO

! STOP

RETURN

END SUBROUTINE SUMSQ

SUBROUTINE GETORDER(QTEMP,QORDER,NATOMS)
USE KEY
IMPLICIT NONE
INTEGER NATOMS
DOUBLE PRECISION QTEMP(3,NATOMS),Q4,Q6,QORDER

IF (COLLAGENOP) THEN
    CALL COLLAGENOP_CALC(QTEMP,QORDER,NATOMS)
ELSE
    CALL QORDER_LJ(QTEMP,Q4,Q6)
    QORDER=Q6
ENDIF

END SUBROUTINE GETORDER

!  WHAM routine using lbfgs minimisation for chisq.
!  Minimise chi^2 statistic to extract best probabilities.
!
SUBROUTINE MYWHAM(NREP,NBINS,PATHSTR,EVAR,APPSTRING,LINVART,OFAIL,NRMIN,NRMAX,ZNORM,FIXZ,ALLZERO)
USE KEY
IMPLICIT NONE
DOUBLE PRECISION DIFF, CHI2, RMS, VTOTAL2, NEGLECT, WTOT, VTOT
DOUBLE PRECISION, ALLOCATABLE :: VTOTAL(:), VAR(:), WIJ(:,:), BESTVAR(:)
DOUBLE PRECISION, ALLOCATABLE :: WEIGHT(:,:)
DOUBLE PRECISION, ALLOCATABLE :: GRAD(:), DGRAD(:)
DOUBLE PRECISION, ALLOCATABLE :: SUMVISITS(:), VISITS(:,:)
!
DOUBLE PRECISION, ALLOCATABLE :: RMIN(:)
INTEGER, ALLOCATABLE :: VISITSALL(:)
INTEGER, ALLOCATABLE :: VISITST(:)
DOUBLE PRECISION, ALLOCATABLE :: MAXVISITS(:)
DOUBLE PRECISION CHI2PLUS, CHI2MINUS, TOL, DPRAND, DUMMY, &
  &              DUMMY1, DUMMY2, DUMMY3, DUMMY4, DUMMY5
INTEGER NBINS, J1, J2, NREP, NVAR, NEVAR, NCOUNT, MUPDATES, ITMAX, ITDONE, &
  &     J3, MERGEI, NBINSSAVE, NDUMMY, NBINSAVE, NRMIN, NRMAX
DOUBLE PRECISION EVAR(NBINS), ZNORM(NREP)
CHARACTER(LEN=80) FSTRING, FNAME, REPSTRING, OFILE
CHARACTER(LEN=10) PATHSTR, APPSTRING
CHARACTER(LEN=1) DUMMYSTRING
LOGICAL, ALLOCATABLE :: OVERLAP(:)
LOGICAL ALLZERO(NBINS)
LOGICAL CONVERGED, YESNO, LINVART, OFAIL, FIXZ
DOUBLE PRECISION, ALLOCATABLE :: EMIN(:)

!
! These are the conditions for a single ts = replica. However,
! the best fit should still work!
!
! IF (NREP.EQ.1) RETURN
! IF (MCPATHTS.NE.0) RETURN
NEGLECT=MCPATHNEGLECT

543 CONTINUE

PRINT '(3(A,I4),2I6)','replicas=',NREP,' # bins=',NBINS,' replica range: ',NRMIN,NRMAX
NBINSAVE=NBINS
MERGEI=1
NBINSSAVE=NBINS
NBINS=NBINS/MERGEI
IF (MERGEI.NE.1) PRINT '(3(A,I4))','bins will be merged in groups of ',MERGEI,' to give ',NBINS
PRINT '(A,G20.10,A)','Bins will be neglected if the number of visits is less than ',NEGLECT,' times the largest value'

ALLOCATE(VTOTAL(NREP))
ALLOCATE(VISITST(NBINS))
ALLOCATE(VISITSALL(NBINSSAVE))
ALLOCATE(WEIGHT(NBINS,NREP),VISITS(NBINS,NREP))
ALLOCATE(SUMVISITS(NBINS))
ALLOCATE(OVERLAP(NREP))
!
!!!!!!!!!!!!!!! Read data from mcpath.<PATHSTR><n> files !!!!!!!!!!!!!!!!
! mcpath.f:
! with bias
!        WRITE(LUNIT,'(I6,8G20.10)') J1,SLENGTH(J1),EOFS(J1),QFRAME(J1),NEARESTFW(J1)/MEANBIAS, &
! &          -LOG(MAX(1.0D-100,NEARESTFW(J1)/MEANBIAS)),1.0D2*NEARESTF(J1)/MCPATHSTEPS, &
! &      QFRAMEAV(J1)/MAX(1,NEARESTF(J1)),SQRT(MAX(1.0D-100,QFRAMESD(J1)/MAX(1,NEARESTF(J1))-(QFRAMEAV(J1)/MAX(1,NEARESTF(J1)))**2))
!     ELSE
!        WRITE(LUNIT,'(I6,8G20.10)') J1,SLENGTH(J1),EOFS(J1),QFRAME(J1),NEARESTFW(J1)/MCPATHSTEPS, &
! &          -LOG(MAX(1.0D-100,NEARESTFW(J1)/MCPATHSTEPS)),1.0D2*NEARESTF(J1)/MCPATHSTEPS, &
! &      QFRAMEAV(J1)/MAX(1,NEARESTF(J1)),SQRT(MAX(1.0D-100,QFRAMESD(J1)/MAX(1,NEARESTF(J1))-(QFRAMEAV(J1)/MAX(1,NEARESTF(J1)))**2))
!
!     ENDIF
!      
! mcpath.Q:
! with bias
! BINLABELQORDER(J1),QORDERHIST(J1)/MEANBIAS,-LOG(MAX(1.0D-100,QORDERHIST(J1)/MEANBIAS)),QORDERVISITS(J1)/MCPATHSTEPS
! without bias
! BINLABELQORDER(J1),QORDERHIST(J1)/MCPATHSTEPS,-LOG(MAX(1.0D-100,QORDERHIST(J1)/MCPATHSTEPS)),QORDERVISITS(J1)/MCPATHSTEPS

!
VISITS(1:NBINS,1:NREP)=0.0D0
WEIGHT(1:NBINS,1:NREP)=0.0D0
VTOTAL2=0.0D0
DO J1=NRMIN,NRMAX
   WRITE(REPSTRING,'(I6)') J1
   WRITE(FSTRING,'(A)') 'mcpath.' // TRIM(ADJUSTL(PATHSTR)) // TRIM(ADJUSTL(REPSTRING)) // TRIM(ADJUSTL(APPSTRING))
   INQUIRE(FILE=TRIM(ADJUSTL(FSTRING)),EXIST=YESNO)
   IF (.NOT.YESNO) THEN
      PRINT '(3A)','File ',TRIM(ADJUSTL(FSTRING)),' not found'
      STOP
   ENDIF
   OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FSTRING)),STATUS='OLD')
   VTOTAL(J1)=0.0D0
   DO J2=1,NBINS
      WTOT=0.0D0; VTOT=0.0D0
      DO J3=1,MERGEI
         IF (PATHSTR(1:1).EQ.'f') THEN
            READ(1,*) NDUMMY,DUMMY1,DUMMY2,DUMMY3,DUMMY4,DUMMY5,DUMMY
         ELSEIF (PATHSTR(1:1).EQ.'Q') THEN
            READ(1,*) DUMMY3,DUMMY4,DUMMY5,DUMMY
         ELSE
            PRINT '(A,A)','mcpath> Unrecognised path string: ' // TRIM(ADJUSTL(PATHSTR))
            STOP
         ENDIF
         IF (PATHSTR(1:1).EQ.'f') THEN
            VTOT=VTOT+DUMMY/1.0D2 ! last column is a percent
         ELSE
            VTOT=VTOT+DUMMY
         ENDIF
         WTOT=WTOT+DUMMY4
      ENDDO
      VISITS(J2,J1)=VTOT
!     IF (VTOT.GT.0.0D0) VISITS(J2,J1)=1.0D0 ! debug DJW
      WEIGHT(J2,J1)=WTOT
      VTOTAL(J1)=VTOTAL(J1)+VISITS(J2,J1) ! summed over replicas
   ENDDO
   DO J2=MERGEI*NBINS,NBINSAVE-1 ! in case MERGEI is not an integer divisor of NBINSAVE
      READ(1,*) DUMMYSTRING
   ENDDO
ENDDO
!
! Check whether the distributions actually overlap!
!
OFAIL=.FALSE.
IF (NRMAX-NRMIN.GT.1) THEN
   OVERLAP(1:NREP)=.FALSE.
   DO J1=NRMIN,NRMAX
      DO J2=J1+1,NRMAX
         dotp: DO J3=1,NBINS
            IF ((VISITS(J3,J1).GT.0.0D0).AND.(VISITS(J3,J2).GT.0.0D0)) THEN
               OVERLAP(J1)=.TRUE.
               OVERLAP(J2)=.TRUE.
               EXIT dotp
            ENDIF
         ENDDO dotp
      ENDDO
      IF (.NOT.OVERLAP(J1)) THEN
         PRINT '(A,I6,A)','mywham> WARNING *** replica ',J1,' has no overlap with the others'
         OFAIL=.TRUE.
      ENDIF
   ENDDO
ENDIF
IF (OFAIL) RETURN

IF (.NOT.ALLOCATED(MAXVISITS)) THEN
   ALLOCATE(MAXVISITS(NREP))
ENDIF
MAXVISITS(1:NREP)=-1.0D0
DO J1=NRMIN,NRMAX
   DO J2=1,NBINS
      IF (VISITS(J2,J1).GT.MAXVISITS(J1)) MAXVISITS(J1)=VISITS(J2,J1)
   ENDDO
   WRITE(*,'(A,I5,A,F20.1,A,F20.10,A,F15.5)') 'for replica ',J1,' sum % visits=',VTOTAL(J1), &
  &        ' neglect thresh ',MAX(0.0D0,MAXVISITS(J1)*NEGLECT),' max % visits=',MAXVISITS(J1)

ENDDO
DO J1=NRMIN,NRMAX
   DO J2=1,NBINS
      IF ((1.0D0*VISITS(J2,J1))/(1.0D0*MAXVISITS(J1)).LT.NEGLECT) THEN
         VISITS(J2,J1)=0.0D0
      ENDIF
   ENDDO
ENDDO
DO J1=NRMIN,NRMAX
! sum histogram
   VTOT=0
   DO J2=1,NBINS
      VTOT = VTOT + VISITS(J2,J1)
   ENDDO
!  DO J2=1,NBINS
!     VISITS(J2,J1) = VISITS(J2,J1)/VTOT
!  ENDDO
ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Solve the simplest lsq fit problem for direct visits to bins
! 
IF (ALLOCATED(WIJ)) DEALLOCATE(WIJ)
ALLOCATE(WIJ(NREP,NBINS))
DO J1=1,NBINS
   DO J2=NRMIN,NRMAX
      IF (VISITS(J1,J2).GT.0.0D0) THEN
         IF (LINVART) THEN
            WIJ(J2,J1)=WEIGHT(J1,J2)
         ELSE
            WIJ(J2,J1)=LOG(WEIGHT(J1,J2))
         ENDIF
      ELSE
         WIJ(J2,J1)=0.0D0
      ENDIF
   ENDDO
ENDDO 
ALLZERO(1:NBINS)=.TRUE.
NEVAR=0
iloop: DO J1=1,NBINS
   DO J2=NRMIN,NRMAX
     IF (VISITS(J1,J2).GT.0) THEN
        ALLZERO(J1)=.FALSE.
        NEVAR=NEVAR+1
        CYCLE iloop
     ENDIF
   ENDDO
ENDDO iloop

NVAR=NEVAR+NREP
PRINT '(A,I10)','mywham> number of bin weight variables=',NEVAR
PRINT '(A,I8)','mywham> number of unknown constants of proportionality for all replicas=',NREP
IF (LINVART) THEN
   PRINT '(A,I8)','mywham> Optimising probabilities, not ln P'
ELSE
   PRINT '(A,I8)','mywham> Optimising ln P'
ENDIF
IF (NEVAR.EQ.0) STOP
!
! Initialise variables. 
!
IF (ALLOCATED(BESTVAR)) DEALLOCATE(BESTVAR)
IF (ALLOCATED(VAR)) DEALLOCATE(VAR)
ALLOCATE(VAR(NVAR),BESTVAR(NVAR))
!
! Initialise variables. Guesses assumed set in calling routine.
!
NCOUNT=0
DO J1=1,NBINS
   IF (ALLZERO(J1)) CYCLE 
   NCOUNT=NCOUNT+1
   VAR(NCOUNT)=EVAR(J1)
ENDDO

VAR(NEVAR+1:NEVAR+NREP)=ZNORM(1:NREP)
!
! Test analytic derivatives
!
IF (.FALSE.) THEN
   ALLOCATE(GRAD(NVAR),DGRAD(NVAR))
   CALL GETCHI(NBINS,NEVAR,NREP,WIJ,VAR,CHI2,GRAD,ALLZERO,VISITS,RMS,LINVART,NRMIN,NRMAX,FIXZ)
   DIFF=0.0001D0
   DO J1=1,NVAR
      VAR(J1)=VAR(J1)+DIFF
      CALL GETCHI(NBINS,NEVAR,NREP,WIJ,VAR,CHI2PLUS,DGRAD,ALLZERO,VISITS,RMS,LINVART,NRMIN,NRMAX,FIXZ)
      VAR(J1)=VAR(J1)-2.0D0*DIFF
      CALL GETCHI(NBINS,NEVAR,NREP,WIJ,VAR,CHI2MINUS,DGRAD,ALLZERO,VISITS,RMS,LINVART,NRMIN,NRMAX,FIXZ)
      VAR(J1)=VAR(J1)+DIFF
      IF (GRAD(J1).NE.0.0D0) THEN
         PRINT '(A,I5,3G20.10)','J1,num,anal,rat=',J1,(CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF),GRAD(J1), &
  &                               (CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF*GRAD(J1))
      ELSE
         PRINT '(A,I5,3G20.10)','J1,num,anal=',J1,(CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF),GRAD(J1)
      ENDIF
   ENDDO
!  STOP
ENDIF

MUPDATES=400
TOL=MCPATHTOL

ITMAX=10000

CALL WHAMLBFGS(NVAR,MUPDATES,VAR,TOL,ITDONE,ITMAX,CHI2,CONVERGED,NBINS,NEVAR, &
  &               NREP,WIJ,ALLZERO,VISITS,RMS,LINVART,NRMIN,NRMAX,FIXZ)

IF (LINVART) THEN
   DO J1=1,NEVAR
      PRINT '(A,I8,A,G20.10)','mcpath> after WHAMLBFGS variable ',J1,' is ',VAR(J1)
      IF (VAR(J1).LT.0.0D0) PRINT '(A,I8,A,G20.10)','mcpath> after WHAMLBFGS variable ',J1,' is negative, value=',VAR(J1)
   ENDDO
ENDIF

BESTVAR(1:NVAR)=VAR(1:NVAR)

NCOUNT=0

DO J1=1,NBINS
   IF (.NOT.ALLZERO(J1)) THEN
      NCOUNT=NCOUNT+1
      IF (LINVART) THEN
         IF (BESTVAR(NCOUNT).LT.-1.0D-10) THEN
            PRINT '(A,I8,A,G20.10)','mcpath> WARNING *** Bin ',J1,' probability is negative: value=',BESTVAR(NCOUNT)
         ENDIF
         IF (BESTVAR(1).NE.0.0D0) VAR(NCOUNT)=BESTVAR(NCOUNT)/BESTVAR(1)
      ELSE
         VAR(NCOUNT)=EXP(BESTVAR(NCOUNT)-BESTVAR(1))
      ENDIF
      EVAR(J1)=BESTVAR(NCOUNT)
   ENDIF
ENDDO
PRINT '(A)','Final bin weights - dumping relative values to weights.' // TRIM(ADJUSTL(PATHSTR))
NCOUNT=0
WRITE(FNAME,'(A)') 'weights.'//TRIM(ADJUSTL(PATHSTR)) // TRIM(ADJUSTL(APPSTRING))
OPEN(UNIT=1,FILE=FNAME,STATUS='UNKNOWN')
DO J1=1,NBINS
   IF (.NOT.ALLZERO(J1)) THEN
      WRITE(1,'(I8,2G20.10)') J1,EVAR(J1)
   ELSE
      WRITE(1,'(I8,2G20.10)') J1,0.0D0
   ENDIF
ENDDO
CLOSE(1)
IF (LINVART) THEN
!  PRINT '(A)','Final bin weights:'
ELSE
!  PRINT '(A)','Final ln(bin weights):'
ENDIF
! PRINT '(I8,G20.10)',(J1,EVAR(J1),J1=1,NBINS)
IF (LINVART) THEN
!  PRINT '(A)','Final normalisation terms:'
ELSE
   PRINT '(A)','Final ln(normalisation) terms:'
ENDIF
  PRINT '(I8,G20.10)',(J1,VAR(NEVAR+J1),J1=NRMIN,NRMAX)
ZNORM(NRMIN:NRMAX)=BESTVAR(NEVAR+NRMIN:NEVAR+NRMAX)

END SUBROUTINE MYWHAM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GETCHI(NBINS,NEVAR,NREP,WIJ,VAR,CHI2,GRAD,ALLZERO,VISITS,RMS,LINVART,NRMIN,NRMAX,FIXZ)
IMPLICIT NONE
INTEGER NEVAR, NREP, NBINS, NCOUNT, J1, J2, NRMIN,NRMAX
DOUBLE PRECISION GRAD(NEVAR+NREP), VAR(NEVAR+NREP), CHI2, RIJ, RMS
DOUBLE PRECISION WIJ(NREP,NBINS)
DOUBLE PRECISION VISITS(NBINS,NREP)
LOGICAL ALLZERO(NBINS),LINVART,FIXZ

CHI2=0.0D0
GRAD(1:NEVAR+NREP)=0.0D0
NCOUNT=0
IF (LINVART) THEN
   elooplin: DO J1=1,NBINS
      IF (ALLZERO(J1)) CYCLE elooplin
      NCOUNT=NCOUNT+1
      DO J2=NRMIN,NRMAX
         RIJ=WIJ(J2,J1)*VAR(NEVAR+J2)-VAR(NCOUNT)
         CHI2=CHI2+VISITS(J1,J2)*RIJ**2
         GRAD(NCOUNT)=GRAD(NCOUNT)    -VISITS(J1,J2)*RIJ ! factor of 2 included outside loop
         GRAD(NEVAR+J2)=GRAD(NEVAR+J2)+VISITS(J1,J2)*RIJ*WIJ(J2,J1)
      ENDDO
   ENDDO elooplin
ELSE
   eloop: DO J1=1,NBINS
      IF (ALLZERO(J1)) CYCLE eloop
      NCOUNT=NCOUNT+1
      DO J2=NRMIN,NRMAX
         RIJ=WIJ(J2,J1)-(VAR(NCOUNT)-VAR(NEVAR+J2))
         CHI2=CHI2+VISITS(J1,J2)*RIJ**2
         GRAD(NCOUNT)=GRAD(NCOUNT)    -VISITS(J1,J2)*RIJ
         GRAD(NEVAR+J2)=GRAD(NEVAR+J2)+VISITS(J1,J2)*RIJ
      ENDDO
   ENDDO eloop
ENDIF

! ss2029 > normalize CHI2 by dividing by number of variables 
CHI2=CHI2/(NEVAR+NRMAX-NRMIN) 

GRAD(1:NEVAR+NREP)=2*GRAD(1:NEVAR+NREP)/(NEVAR+NRMAX-NRMIN)
IF (FIXZ) GRAD(NEVAR+1:NEVAR+NREP)=0.0D0 ! freeze all normalisations
! GRAD(NEVAR+NRMIN)=0.0D0 ! freeze first normalisation
RMS=0.0D0
DO J1=1,NEVAR+NREP
   RMS=RMS+GRAD(J1)**2
ENDDO
RMS=SQRT(RMS/(NEVAR+NREP))

END SUBROUTINE GETCHI

SUBROUTINE WHAMLBFGS(N,M,X,EPS,ITDONE,ITMAX,ENERGY,CONVERGED,NBINS,NEVAR,NREP,WIJ,ALLZERO,VISITS,RMS,LINVART,NRMIN,NRMAX,FIXZ)
IMPLICIT NONE
INTEGER N,M,J1,ITMAX,ITDONE,NFAIL,NRMIN,NRMAX
DOUBLE PRECISION X(N),G(N),DIAG(N),W(N*(2*M+1)+2*M),SLENGTH,DDOT
DOUBLE PRECISION EPS,ENERGY,ENEW,GNEW(N),RMS
LOGICAL CONVERGED, FIXZ
DOUBLE PRECISION GNORM,STP,YS,YY,SQ,YR,BETA,OVERLAP,DOT1,DOT2,DGUESS,MAXBFGS
INTEGER ITER,POINT,ISPT,IYPT,BOUND,NPT,CP,I,INMC,IYCN,ISCN,NDECREASE
LOGICAL MFLAG, FAILED
INTEGER NBINS, NEVAR, NREP
DOUBLE PRECISION WIJ(NREP,NBINS)
DOUBLE PRECISION VISITS(NBINS,NREP)
LOGICAL ALLZERO(NBINS), LINVART
DOUBLE PRECISION DIFF, CHI2, CHI2PLUS, CHI2MINUS, GRAD(N), DGRAD(N)

DGUESS=1.0D0
MAXBFGS=1000.0D0
ITER=0
ITDONE=0
NFAIL=0
FAILED=.FALSE.
CALL GETCHI(NBINS,NEVAR,NREP,WIJ,X,ENERGY,G,ALLZERO,VISITS,RMS,LINVART,NRMIN,NRMAX,FIXZ)

! WRITE(*,'(A,2G20.10,A,I6,A)') 'mylbfgs> chi^2 and RMS force=',ENERGY,RMS,' after ',ITDONE,' steps'

10    MFLAG=.FALSE.
IF (RMS.LE.EPS) THEN
   MFLAG=.TRUE.
   IF (MFLAG) THEN
      WRITE(*,'(A,I8,A,G15.5,A,G15.5)') 'mylbfgs> chi^2 converged in ',ITDONE, &
  &                                                ' steps. Value=',ENERGY,' RMS force=',RMS
      IF (ITER.GT.0) WRITE(*,'(A,F20.10)') 'mylbfgs> Diagonal inverse Hessian elements are now ',DIAG(1)
      CONVERGED=.TRUE.
      RETURN
   ENDIF
ENDIF

IF ((ITDONE.EQ.ITMAX).OR.FAILED) THEN
   WRITE(*,'(A,G15.7,A,G15.7)') 'mylbfgs> **WARNING - chi^2 did not converge, value=',ENERGY,' RMS force=',RMS
   WRITE(*,'(A,F20.10)') 'mylbfgs> Diagonal inverse Hessian elements are now ',DIAG(1)
   CONVERGED=.FALSE.
   RETURN
ENDIF

IF (ITER.EQ.0) THEN
   IF (N.LE.0.OR.M.LE.0) THEN
      WRITE(*,240)
240   FORMAT('xmylbfgs> IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
      STOP
   ENDIF
   POINT=0
   MFLAG=.FALSE.
   DO I=1,N
      DIAG(I)=DGUESS
   ENDDO
   ISPT= N+2*M
   IYPT= ISPT+N*M
!
!  NR step for diagonal inverse Hessian
!
   DO I=1,N
      W(ISPT+I)= -G(I)*DIAG(I)
      W(I)= -G(I)*DIAG(I)
   ENDDO
   GNORM= DSQRT(DDOT(N,G,1,G,1))
!
!  Make the first guess for the step length cautious.
!
   IF (GNORM.EQ.0.0D0) THEN
      GNORM=1.0D0 ! exact zero is presumably wrong!
      PRINT '(A)','WARNING - GNORM was zero in xmylbfgs, resetting to one'
   ENDIF
   STP=MIN(GNORM,1.0D0/GNORM)
ELSE 
   BOUND=ITER
   IF (ITER.GT.M) BOUND=M
   YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
   IF (YS.EQ.0.0D0) YS=1.0D0
!
!  Update estimate of diagonal inverse Hessian elements
!  We divide by both YS and YY at different points, so
!  they had better not be zero!
!
   YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
   IF (YY.EQ.0.0D0) YY=1.0D0
   DO I=1,N
      DIAG(I)= YS/YY
!     DIAG(I)= ABS(YS/YY)
   ENDDO
!
!     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
!     "Updating quasi-Newton matrices with limited storage",
!     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
!     ---------------------------------------------------------
!
   CP= POINT
   IF (POINT.EQ.0) CP=M
   W(N+CP)= 1.0D0/YS
   DO I=1,N
      W(I)= -G(I)
   ENDDO
   CP= POINT
   DO I= 1,BOUND
      CP=CP-1
      IF (CP.EQ. -1)CP=M-1
      SQ= DDOT(N,W(ISPT+CP*N+1),1,W,1)
      INMC=N+M+CP+1
      IYCN=IYPT+CP*N
      W(INMC)= W(N+CP+1)*SQ
      CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
   ENDDO
  
   DO I=1,N
      W(I)=DIAG(I)*W(I)
   ENDDO

   DO I=1,BOUND
      YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
      BETA= W(N+CP+1)*YR
      INMC=N+M+CP+1
      BETA= W(INMC)-BETA
      ISCN=ISPT+CP*N
      CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
      CP=CP+1
      IF (CP.EQ.M) CP=0
   ENDDO
   STP=1.0D0
ENDIF
!
!  Store the new search direction
!
! W(1)=0.0D0 ! freeze first weight
DO I=1,N
   W(ISPT+POINT*N+I)= W(I)
ENDDO

DOT1=SQRT(DDOT(N,G,1,G,1))
DOT2=SQRT(DDOT(N,W,1,W,1))
OVERLAP=0.0D0
IF (DOT1*DOT2.NE.0.0D0) OVERLAP=DDOT(N,G,1,W,1)/(DOT1*DOT2)

IF (OVERLAP.GT.0.0D0) THEN
   PRINT*,'Search direction has positive projection onto gradient - reversing step'
   DO I=1,N
      W(ISPT+POINT*N+I)= -W(I)
   ENDDO
ENDIF
      
DO I=1,N
   W(I)=G(I)
ENDDO
SLENGTH=0.0D0
DO J1=1,N
   SLENGTH=SLENGTH+W(ISPT+POINT*N+J1)**2
ENDDO
SLENGTH=SQRT(SLENGTH)
IF (STP*SLENGTH.GT.MAXBFGS) STP=MAXBFGS/SLENGTH
!
!  We now have the proposed step.
!
DO J1=1,N
   X(J1)=X(J1)+STP*W(ISPT+POINT*N+J1)
ENDDO 
NDECREASE=0
20    CONTINUE
CALL GETCHI(NBINS,NEVAR,NREP,WIJ,X,ENEW,GNEW,ALLZERO,VISITS,RMS,LINVART,NRMIN,NRMAX,FIXZ)

IF (ENEW.EQ.0.0D0) ENEW=1.0D-100 ! to prevent divide by zero
IF (((ENEW-ENERGY)/ABS(ENEW).LE.1.0D-6)) THEN
   ITER=ITER+1
   ITDONE=ITDONE+1
   ENERGY=ENEW
   DO J1=1,N
      G(J1)=GNEW(J1)
   ENDDO
!  WRITE(*,'(A,2G20.10,A,I6,A,F13.6)') &
! &                'mylbfgs> chi^2 and RMS force=',ENERGY,RMS,' after ',ITDONE,' steps, step:',STP*SLENGTH
ELSE
!
!  chi^2 increased - try again with a smaller step size?
!
   IF (NDECREASE.GT.2) THEN  
      NFAIL=NFAIL+1
      PRINT*,' LBFGS step cannot find a lower value, NFAIL=',NFAIL

      CALL GETCHI(NBINS,NEVAR,NREP,WIJ,X,CHI2,GRAD,ALLZERO,VISITS,RMS,LINVART,NRMIN,NRMAX,FIXZ)
      DIFF=0.0001D0
      DO J1=1,N
         X(J1)=X(J1)+DIFF
         CALL GETCHI(NBINS,NEVAR,NREP,WIJ,X,CHI2PLUS,DGRAD,ALLZERO,VISITS,RMS,LINVART,NRMIN,NRMAX,FIXZ)
         X(J1)=X(J1)-2.0D0*DIFF
         CALL GETCHI(NBINS,NEVAR,NREP,WIJ,X,CHI2MINUS,DGRAD,ALLZERO,VISITS,RMS,LINVART,NRMIN,NRMAX,FIXZ)
         X(J1)=X(J1)+DIFF
         IF (GRAD(J1).NE.0.0D0) THEN
            PRINT '(A,I5,3G20.10)','J1,num,anal,rat=',J1,(CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF),GRAD(J1), &
  &                               (CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF*GRAD(J1))
         ELSE
            PRINT '(A,I5,3G20.10)','J1,num,anal=',J1,(CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF),GRAD(J1)
         ENDIF
      ENDDO
      STOP

      ITER=0  !  try resetting
      DO J1=1,N
         X(J1)=X(J1)-STP*W(ISPT+POINT*N+J1)
      ENDDO 
      IF (NFAIL.GT.1) FAILED=.TRUE.
      GOTO 10
   ENDIF
   DO J1=1,N
      X(J1)=X(J1)-0.9*STP*W(ISPT+POINT*N+J1)
   ENDDO 
   NDECREASE=NDECREASE+1
   STP=STP/10.0D0
!  WRITE(*,'(A,G20.10,A,G20.10,A,F15.8)') &
! &                         'mylbfgs> chi^2 increased from ',ENERGY,' to ',ENEW,' decreasing step to ',STP*SLENGTH
   GOTO 20
ENDIF
!
!  Compute the new step and gradient change
!
NPT=POINT*N
DO I=1,N
   W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
   W(IYPT+NPT+I)= G(I)-W(I)
ENDDO

POINT=POINT+1
IF (POINT.EQ.M) POINT=0
GOTO 10

RETURN

END SUBROUTINE WHAMLBFGS
