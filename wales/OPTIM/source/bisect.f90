!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999-2006 David J. Wales
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
!  BISECT finds a local minima between two starting points in a systematic way.
!
SUBROUTINE BISECT_OPT(NATOMS,ESTART,CSTART,EFINISH,CFINISH,DENDPT)
USE KEY,ONLY : BISECTSTEPS, MUPDATE, BULKT, TWOD, RIGIDBODY, BISECTMINDIST, BISECTMAXATTEMPTS, DUMPDATAT, &
  &            BFGSSTEPS, BISECTDEBUG, BISECTMAXENERGY, GEOMDIFFTOL, EDIFFTOL, PERMDIST, GMAX, BHSFRAC, &
  &            REOPTIMISEENDPOINTS
USE MODCHARMM,ONLY : CHRMMT, ICINTERPT, CHECKOMEGAT,CHECKCHIRALT,INTMINT,NOCISTRANS
USE MODAMBER9, ONLY : NOCISTRANSRNA, GOODSTRUCTURE1
USE PORFUNCS
USE COMMONS,ONLY: PARAM1,PARAM2,PARAM3,ZSYM,DEBUG
! USE MODEFOL
IMPLICIT NONE
INTEGER ITDONE, NATOMS, J1, J2, ISTAT, K, NSTEPS, NBISECT, MAXBISECT, NCOUNT, NMIN, NDO
DOUBLE PRECISION CSTART(3*NATOMS), CFINISH(3*NATOMS)
DOUBLE PRECISION COORDS(3*NATOMS), VNEW(3*NATOMS), DSTART, DFINISH, RMAT(3,3), &
  &              EREAL, RMS2, EBEST, EPREV, RMS, ENERGY, DELTAX(3*NATOMS), &
  &              TEMPSTART(3*NATOMS), TEMPFINISH(3*NATOMS), CTEMP(3*NATOMS)
DOUBLE PRECISION SFRAC, DSTARTSAVE, DFINISHSAVE, ESTART, EFINISH, SLIMIT, FLIMIT, MAXDIST
DOUBLE PRECISION DIST, DIST2, GMAXSAVE, DENDPT
DOUBLE PRECISION XVECS(3*NATOMS),XENERGY,XEVALMIN,XVNEW(3*NATOMS)
LOGICAL XMFLAG

LOGICAL MFLAG, ATEST, PTEST, HITS, HITF, CHIRALFAIL, AMIDEFAIL
INTEGER :: MAXMIN=10
DOUBLE PRECISION, ALLOCATABLE :: EMIN(:), MINCOORDS(:,:), DISTANCES(:)
LOGICAL, ALLOCATABLE :: TRIED(:)
LOGICAL KNOWE, KNOWG, KNOWH
COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
DOUBLE PRECISION, ALLOCATABLE :: VDP(:), VDIST(:,:)
LOGICAL, ALLOCATABLE :: VLOG(:)
DOUBLE PRECISION BISECTENERGY
COMMON /BISECTE/ BISECTENERGY
CHARACTER(LEN=80) FNAMEF
CHARACTER(LEN=20) EFNAME

IF (CHRMMT.AND.INTMINT) THEN
   PRINT '(A)',' bisect> CHARMM internal coodinate minimisation not allowed with BISECT'
   STOP
ENDIF
IF (.NOT.ALLOCATED(EMIN)) ALLOCATE(EMIN(MAXMIN))
IF (.NOT.ALLOCATED(MINCOORDS)) ALLOCATE(MINCOORDS(3*NATOMS,MAXMIN))
IF (.NOT.ALLOCATED(DISTANCES)) ALLOCATE(DISTANCES(MAXMIN))
IF (.NOT.ALLOCATED(TRIED)) ALLOCATE(TRIED(MAXMIN))
EMIN(1)=ESTART; EMIN(2)=EFINISH
MINCOORDS(1:3*NATOMS,1)=CSTART(1:3*NATOMS)
MINCOORDS(1:3*NATOMS,2)=CFINISH(1:3*NATOMS)
NMIN=2
DISTANCES(1)=DENDPT
TRIED(1:MAXMIN)=.FALSE.
!
! Move through list of minima to find the largest distance for which we have not yet
! tried to quench to an intervening minimum. Stop if this gap is larger than BISECTMINDIST
! Must distinguish jumps that simply change SFRAC and where we want to try a new pair of
! minima.
!
NSTEPS=0

11 CONTINUE ! jump here to try new endpoints
IF (BISECTDEBUG) THEN
   PRINT '(A,I8,A)',' bisect> After ',NSTEPS,' steps '
   PRINT '(A)','   minimum       energy       distance to next min       tried'
   DO J1=1,NMIN-1
      PRINT '(I8,F20.10,F20.10,7X,L5)',J1,EMIN(J1),DISTANCES(J1),TRIED(J1)
   ENDDO
   PRINT '(I8,F20.10)',NMIN,EMIN(NMIN)
ENDIF
NSTEPS=NSTEPS+1
SLIMIT=1.0D0
FLIMIT=0.0D0
IF (NSTEPS.GT.BISECTSTEPS) GOTO 13
NCOUNT=1
NBISECT=0
SFRAC=0.5D0
NDO=-1
MAXDIST=-10.0D0
DO J1=1,NMIN-1
   IF (TRIED(J1)) CYCLE
   IF (DISTANCES(J1).LT.BISECTMINDIST) CYCLE
   IF (DISTANCES(J1).GT.MAXDIST) THEN
      MAXDIST=DISTANCES(J1)
      NDO=J1
   ENDIF
ENDDO
IF (NDO.EQ.-1) THEN
   PRINT '(A)',' bisect> no eligible gaps left in bisect'
   STOP !!! need to dump minima first !!! DJW
ENDIF
12 CONTINUE ! jump here for same endpoints, different SFRAC
NBISECT=NBISECT+1
PRINT '(A,I8,A,I8,A,F15.5,2(A,I8))',' bisect> trying minima in positions ',NDO,' and ',NDO+1,' fraction=',SFRAC,' step ', &
  &                        NSTEPS,' attempt ',NBISECT
HITS=.FALSE.
HITF=.FALSE.
TEMPSTART(1:3*NATOMS)=MINCOORDS(1:3*NATOMS,NDO)
TEMPFINISH(1:3*NATOMS)=MINCOORDS(1:3*NATOMS,NDO+1)
IF (BULKT) THEN
   DO K=1,NATOMS
      DELTAX(3*(K-1)+1)=MINCOORDS(3*(K-1)+1,NDO+1) - MINCOORDS(3*(K-1)+1,NDO) &
  &       -PARAM1*NINT((MINCOORDS(3*(K-1)+1,NDO+1) - MINCOORDS(3*(K-1)+1,NDO))/PARAM1)
      DELTAX(3*(K-1)+2)=MINCOORDS(3*(K-1)+2,NDO+1) - MINCOORDS(3*(K-1)+2,NDO) &
  &       -PARAM2*NINT((MINCOORDS(3*(K-1)+2,NDO+1) - MINCOORDS(3*(K-1)+2,NDO))/PARAM2)
      DELTAX(3*(K-1)+3)=MINCOORDS(3*(K-1)+3,NDO+1) - MINCOORDS(3*(K-1)+3,NDO) &
  &       -PARAM3*NINT((MINCOORDS(3*(K-1)+3,NDO+1) - MINCOORDS(3*(K-1)+3,NDO))/PARAM3)
   ENDDO
   COORDS(1:3*NATOMS)=SFRAC*MINCOORDS(1:3*NATOMS,NDO)+(1.0D0-SFRAC)*DELTAX(1:3*NATOMS)
ELSE
   COORDS(1:3*NATOMS)=SFRAC*CSTART(1:3*NATOMS)+(1.0D0-SFRAC)*CFINISH(1:3*NATOMS)
   IF(CHRMMT.AND.ICINTERPT) CALL ICINTERPOL(COORDS,TEMPSTART,TEMPFINISH,SFRAC)
ENDIF
!
! Minimsation
!
KNOWE=.FALSE.
PTEST=DEBUG
CALL MYLBFGS(3*NATOMS,MUPDATE,COORDS,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,.TRUE.,ITDONE,PTEST,VNEW,.TRUE.,.FALSE.)
!
! Action for AMIDEFAIL or CHIRALFAIL in interpolated minimum.
!
IF (CHRMMT) THEN
   AMIDEFAIL=.FALSE.
   IF (CHECKOMEGAT) THEN
      CALL CHECKOMEGA(COORDS,AMIDEFAIL)
      IF(BISECTDEBUG) PRINT '(A,L5)',' bisect> interpolated geometry: AMIDEFAIL  = ',AMIDEFAIL
   ENDIF
   CHIRALFAIL=.FALSE.
   IF (CHECKCHIRALT) THEN
      CALL CHECKCHIRAL(COORDS,CHIRALFAIL)
      IF(BISECTDEBUG) PRINT '(A,L5)',' bisect> interpolated geometry: CHIRALFAIL = ',CHIRALFAIL
   ENDIF
   IF (CHIRALFAIL.OR.AMIDEFAIL) THEN
      SFRAC=0.5D0+(-1.D0)**NBISECT*NCOUNT*0.1D0
      IF ((SFRAC.LE.0).OR.(SFRAC.GE.1.D0)) THEN
         PRINT *,' bisect> bisection attempt failed for this pair of minima'
         GOTO 11
      ENDIF
      IF(MOD(NBISECT,2).EQ.0) NCOUNT=NCOUNT+1
      PRINT '(A,F8.3)',' bisect> retry interpolation, SFRAC= ',SFRAC
      GOTO 12
   ENDIF
ENDIF
IF (.NOT.MFLAG) THEN
   SFRAC=0.5D0+(-1.D0)**NBISECT*NCOUNT*0.1D0
   IF ((SFRAC.LE.0).OR.(SFRAC.GE.1.D0)) THEN
      PRINT *,' bisect> bisection attempt failed for this pair of minima'
      GOTO 11
   ENDIF
   IF(MOD(NBISECT,2).EQ.0) NCOUNT=NCOUNT+1
   PRINT '(A,F8.3)',' bisect> retry interpolation, SFRAC= ',SFRAC
   GOTO 12
ENDIF
IF ((.NOT.MFLAG).OR.(EREAL.GT.BISECTMAXENERGY)) THEN
   TRIED(NDO)=.TRUE.
   GOTO 11
ENDIF

WRITE(*,*) "sn402: Changed BISECT so that it calls ALIGN_DECIDE instead of MINPERMDIST"
WRITE(*,*) "I haven't tested this change and am not certain whether it's sensible." 
WRITE(*,*) "Please set FASTOVERLAP and check carefully that this part of the code is working as you expect, then remove these messages!"
CALL ALIGN_DECIDE(COORDS,TEMPSTART,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
DSTART=DIST
CALL ALIGN_DECIDE(COORDS,TEMPFINISH,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
DFINISH=DIST
! PRINT '(A,4G18.8,I8)','EREAL,DS,DF,SFRAC,ITDONE=',EREAL,DSTART,DFINISH,SFRAC,ITDONE
IF (.NOT.MFLAG) THEN
   PRINT '(A,G20.10,A,I8)',' bisect> WARNING - initial quench failed to converge, energy=',EREAL, &
  &  ' lbfgs steps=',ITDONE
   PRINT '(A,2G15.5)',' bisect> distances from bracketing minima: ',DSTART,DFINISH
ENDIF
IF (BISECTDEBUG) PRINT '(A,G20.10,A,I8)',' bisect> quench energy=',EREAL,' lbfgs steps=',ITDONE
IF (BISECTDEBUG) PRINT '(A,2G15.5)',' bisect> distances from bracketing minima: ',DSTART,DFINISH

IF ((DSTART.LT.GEOMDIFFTOL).AND.(ABS(EREAL-EMIN(NDO)).LT.EDIFFTOL)) THEN
   IF (NBISECT.EQ.BISECTMAXATTEMPTS) THEN
      IF (BISECTDEBUG) PRINT '(A)',' bisect> Bisection limit reached - abandon interpolation for this pair'
      CALL FLUSH(6)
      TRIED(NDO)=.TRUE.
      GOTO 11
   ENDIF
   IF (HITF) THEN
      IF (BISECTDEBUG) PRINT '(A)',' bisect> Both end points hit - abandon interpolation for this pair'
      CALL FLUSH(6)
      TRIED(NDO)=.TRUE.
      GOTO 11
   ENDIF
   SLIMIT=MIN(SFRAC,SLIMIT)
   SFRAC=(SFRAC+FLIMIT)/2.0D0
   IF (BISECTDEBUG) PRINT '(A,G20.10)',' bisect> Initial guess minimised to first end point - change fraction to ',SFRAC
   HITS=.TRUE.
   GOTO 12
ENDIF
IF ((DFINISH.LT.GEOMDIFFTOL).AND.(ABS(EREAL-EMIN(NDO+1)).LT.EDIFFTOL)) THEN
   IF (NBISECT.EQ.BISECTMAXATTEMPTS) THEN
      IF (BISECTDEBUG) PRINT '(A)',' bisect> Bisection limit reached - abandon interpolation for this pair'
      CALL FLUSH(6)
      TRIED(NDO)=.TRUE.
      GOTO 11
   ENDIF
   IF (HITS) THEN
      IF (BISECTDEBUG) PRINT '(A)',' bisect> Both end points hit - abandon interpolation for this pair'
      CALL FLUSH(6)
      TRIED(NDO)=.TRUE.
      GOTO 11
   ENDIF
   FLIMIT=MAX(SFRAC,FLIMIT)
   SFRAC=(SFRAC+SLIMIT)/2.0D0
   IF (BISECTDEBUG) PRINT '(A,G20.10)',' bisect> Initial guess minimised to second end point - change fraction to ',SFRAC
   HITF=.TRUE.
   GOTO 12
ENDIF
!
! Now check that the minimum isn;t identical to any of the others. Don't add it to the list if so.
!
DO J1=1,NMIN
   IF (ABS(EREAL-EMIN(J1)).GT.EDIFFTOL) CYCLE
   CTEMP(1:3*NATOMS)=MINCOORDS(1:3*NATOMS,J1)
   CALL ALIGN_DECIDE(COORDS,CTEMP,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
   IF (DIST.LT.GEOMDIFFTOL) THEN
      PRINT '(A,I8,A)',' bisect> quench minimum is the same as number ',J1,' not adding it to the list'
      TRIED(NDO)=.TRUE.
      GOTO 11
   ENDIF
ENDDO

!
! If we have reached this point then we need to add the new minimum to the list and 
! then go back to the beginning. Must set the TRIED values appropriately. Connections to 
! the new minimum will be untried.
!
NMIN=NMIN+1
IF (NMIN.GT.MAXMIN) THEN

   ALLOCATE(VDP(MAXMIN),VLOG(MAXMIN),VDIST(3*NATOMS,MAXMIN))

   VDP(1:MAXMIN)=EMIN(1:MAXMIN)
   DEALLOCATE(EMIN)
   ALLOCATE(EMIN(2*MAXMIN))
   EMIN(1:MAXMIN)=VDP(1:MAXMIN)

   VDP(1:MAXMIN)=DISTANCES(1:MAXMIN)
   DEALLOCATE(DISTANCES)
   ALLOCATE(DISTANCES(2*MAXMIN))
   DISTANCES(1:MAXMIN)=VDP(1:MAXMIN)

   VDIST(1:3*NATOMS,1:MAXMIN)=MINCOORDS(1:3*NATOMS,1:MAXMIN)
   DEALLOCATE(MINCOORDS)
   ALLOCATE(MINCOORDS(3*NATOMS,2*MAXMIN))
   MINCOORDS(1:3*NATOMS,1:MAXMIN)=VDIST(1:3*NATOMS,1:MAXMIN)

   VLOG(1:MAXMIN)=TRIED(1:MAXMIN)
   DEALLOCATE(TRIED)
   ALLOCATE(TRIED(2*MAXMIN))
   TRIED(1:MAXMIN)=VLOG(1:MAXMIN)

   MAXMIN=2*MAXMIN
   DEALLOCATE(VDP,VLOG,VDIST)
ENDIF
!
!  Move existing entries down after NDO.
!
DO J1=NMIN,NDO+2,-1
   EMIN(J1)=EMIN(J1-1)
   DISTANCES(J1)=DISTANCES(J1-1)
   TRIED(J1)=TRIED(J1-1)
   MINCOORDS(1:3*NATOMS,J1)=MINCOORDS(1:3*NATOMS,J1-1)
ENDDO
EMIN(NDO+1)=EREAL
DISTANCES(NDO)=DSTART
DISTANCES(NDO+1)=DFINISH
TRIED(NDO)=.FALSE.
TRIED(NDO+1)=.FALSE.
MINCOORDS(1:3*NATOMS,NDO+1)=COORDS(1:3*NATOMS)

GOTO 11

13 CONTINUE

IF (DUMPDATAT) THEN 
   FNAMEF='points.final.bh'
   EFNAME='  '
   DO J1=2,NMIN-1 ! don't include the end point minima
      REOPTIMISEENDPOINTS=.FALSE.
      BISECTENERGY=EMIN(J1)
      CTEMP(1:3*NATOMS)=MINCOORDS(1:3*NATOMS,J1)
      CALL GEOPT(FNAMEF,EFNAME,CTEMP,XVECS,XMFLAG,XENERGY,XEVALMIN,XVNEW) ! file names are dummies here
      CALL FLUSH(881)    ! flush min.data.info unit again
   ENDDO
ENDIF

STOP

CALL FLUSH(6)

END SUBROUTINE BISECT_OPT
