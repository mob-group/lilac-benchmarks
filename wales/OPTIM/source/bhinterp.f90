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
!  BHINTERP finds a local minimum between two end point minima using
!  basin-hopping global optimisation based on the actual PE plus the
!  energy of two srings connected to CSTART and CFINISH
!  Assume we only use the LBFGS minimiser.
!  Allow a different force constant, BHK, from NEB phase.
!
SUBROUTINE BHINTERP(CSTART,CFINISH,NCOORDS,NATOMS,INTERPCOORDS,SUCCESS,DSTART,DFINISH,FINALENERGY,ESTART,EFINISH,DENDPT)
USE KEY,ONLY : BHSTEPS, BHCONV, BHTEMP, BHACCREJ, MUPDATE, BULKT, TWOD, RIGIDBODY, &
  &            BHK, BFGSSTEPS, BHMAXENERGY, GEOMDIFFTOL, EDIFFTOL, PERMDIST, GMAX, BHSFRAC, AMBERT, NABT, &
  &  BHINTERPUSELOWEST, BHCHECKENERGYT, BHSTEPSMIN, CLOSESTALIGNMENT, BHDEBUG, BHSTEPSIZE
USE MODCHARMM,ONLY : CHRMMT, ICINTERPT, CHECKOMEGAT,CHECKCHIRALT,INTMINT,NOCISTRANS,MINOMEGA
USE MODAMBER9, ONLY : NOCISTRANSRNA, NOCISTRANSDNA, GOODSTRUCTURE1,GOODSTRUCTURE2,CISARRAY1,&
CISARRAY2, lenic, NICTOT, AMBERICT, AMBSTEPT, nphia, nphih, AMBPERTT
USE PORFUNCS
USE COMMONS,ONLY: PARAM1,PARAM2,PARAM3,ZSYM,DEBUG
USE MODEFOL
USE SPFUNCTS, ONLY : DUMPCOORDS
IMPLICIT NONE
INTEGER, INTENT(IN) :: NCOORDS
DOUBLE PRECISION CSTART(NCOORDS), CFINISH(NCOORDS)
DOUBLE PRECISION COORDS(NCOORDS), VNEW(NCOORDS), COORDSO(NCOORDS), DSTART, DFINISH, RMAT(3,3), ESTRING, ESPREV, RANDOM, DPRAND, &
  &              BESTCOORDS(NCOORDS), EREAL, RMS2, P0, ETOTAL, EBEST, EPREV, RMS, ENERGY, DELTAX(NCOORDS)
DOUBLE PRECISION INTERPCOORDS(NCOORDS), FINALENERGY, SFRAC, DSTARTSAVE, DFINISHSAVE, ESTART, EFINISH, SLIMIT, FLIMIT
DOUBLE PRECISION DIST, DIST2, GMAXSAVE, DENDPT
DOUBLE PRECISION,DIMENSION(:), ALLOCATABLE:: PHI_START, PHI_FINISH
INTEGER ITDONE, NSUCCESS, NFAIL, NSUCCESST, NFAILT, NATOMS, J1, J2, ISTAT, K, NSTEPS, NBISECT, MAXBISECT, NTRY, NCOUNT
LOGICAL MFLAG, ATEST, SUCCESS, PTEST, HITS, HITF, CHIRALFAIL, AMIDEFAIL, EXITBHRUN
LOGICAL KNOWE, KNOWG, KNOWH
COMMON /KNOWN/ KNOWE, KNOWG, KNOWH

! 
! Calculate initial energy
!
IF (CHRMMT.AND.INTMINT) THEN
   PRINT '(A)','bhinterp> CHARMM internal coodinate minimisation not allowed with BHINTERP'
   STOP
ENDIF
MAXBISECT=5
NBISECT=0
NTRY=0
NCOUNT=1
SLIMIT=1.0D0
FLIMIT=0.0D0
SFRAC=BHSFRAC ! BHSFRAC is the fixed value specified on the BHINTERP line
HITS=.FALSE.
HITF=.FALSE.
1 CONTINUE
! NSTEPS=100
! DO J1=1,NSTEPS-1

! SFRAC=(J1*1.0D0)/NSTEPS
IF(BHDEBUG) PRINT '(A,F8.3)',' bhinterp> start interpolation, SFRAC= ',SFRAC
10 CONTINUE
IF (BULKT) THEN
   DO K=1,NATOMS
      DELTAX(3*(K-1)+1)=CFINISH(3*(K-1)+1) - CSTART(3*(K-1)+1) &
  &       -PARAM1*NINT((CFINISH(3*(K-1)+1) - CSTART(3*(K-1)+1))/PARAM1)
      DELTAX(3*(K-1)+2)=CFINISH(3*(K-1)+2) - CSTART(3*(K-1)+2) &
  &       -PARAM2*NINT((CFINISH(3*(K-1)+2) - CSTART(3*(K-1)+2))/PARAM2)
      DELTAX(3*(K-1)+3)=CFINISH(3*(K-1)+3) - CSTART(3*(K-1)+3) &
  &       -PARAM3*NINT((CFINISH(3*(K-1)+3) - CSTART(3*(K-1)+3))/PARAM3)
   ENDDO
   COORDS(1:NCOORDS)=SFRAC*CSTART(1:NCOORDS)+(1.0D0-SFRAC)*DELTAX(1:NCOORDS)
ELSE
   COORDS(1:NCOORDS)=SFRAC*CSTART(1:NCOORDS)+(1.0D0-SFRAC)*CFINISH(1:NCOORDS)
   IF(CHRMMT.AND.ICINTERPT) CALL ICINTERPOL(COORDS,CSTART,CFINISH,SFRAC)
   IF ((AMBERT.OR.NABT).AND.AMBERICT) THEN
      PRINT*, "AMBER internal coordinate interpolation"
      IF (.NOT.ALLOCATED(NICTOT)) THEN
          CALL SETDIHEAM() 
      ENDIF
      CALL TAKESTEPAMDIHED(COORDS, CSTART, CFINISH,SFRAC) !msb50
   ENDIF
ENDIF
IF(BHDEBUG) CALL DUMPCOORDS(COORDS, 'midpoint.xyz', .FALSE.)
!bs360
IF (CHRMMT) THEN
   AMIDEFAIL=.FALSE.
   IF (CHECKOMEGAT) THEN
      CALL CHECKOMEGA(COORDS,AMIDEFAIL)
      IF(BHDEBUG) PRINT '(A,L5)','           interpolated geometry: AMIDEFAIL  = ',AMIDEFAIL
   ENDIF
   CHIRALFAIL=.FALSE.
   IF (CHECKCHIRALT) THEN
      CALL CHECKCHIRAL(COORDS,CHIRALFAIL)
      IF(BHDEBUG) PRINT '(A,L5)','           interpolated geometry: CHIRALFAIL = ',CHIRALFAIL
   ENDIF
   IF (CHIRALFAIL.OR.AMIDEFAIL) THEN
       NTRY=NTRY+1
      SFRAC=BHSFRAC+(-1.D0)**NTRY*NCOUNT*0.05D0
      IF ((SFRAC.LE.0).OR.(SFRAC.GE.1.D0)) THEN
         PRINT *,' bhinterp> could not find acceptable interpolated minimum'
         RETURN
      ENDIF
      IF(MOD(NTRY,2).EQ.0) NCOUNT=NCOUNT+1
      PRINT '(A,F8.3)','bhinterp> retry interpolation, SFRAC= ',SFRAC
      GOTO 10
   ENDIF
ENDIF

KNOWE=.FALSE.
PTEST=DEBUG

! WRITE(99,*) NATOMS
! WRITE(99,'(A)') 'start:'
! DO K=1,NATOMS
!    WRITE(99,'(A,3G20.10)') 'LA ',CSTART(3*(K-1)+1),CSTART(3*(K-1)+2),CSTART(3*(K-1)+3)
! ENDDO
! WRITE(99,*) NATOMS
! WRITE(99,'(A)') 'finish:'
! DO K=1,NATOMS
!    WRITE(99,'(A,3G20.10)') 'LA ',CFINISH(3*(K-1)+1),CFINISH(3*(K-1)+2),CFINISH(3*(K-1)+3)
! ENDDO
! WRITE(99,*) NATOMS
! WRITE(99,'(A)') 'interp:'
! DO K=1,NATOMS
!    WRITE(99,'(A,3G20.10)') 'LA ',COORDS(3*(K-1)+1),COORDS(3*(K-1)+2),COORDS(3*(K-1)+3)
! ENDDO

GMAXSAVE=GMAX; GMAX=BHCONV ! mylbfgs now uses GMAX so that we can change this parameter via changep
CALL MYLBFGS(NCOORDS,MUPDATE,COORDS,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS, &
     &                   .TRUE.,ITDONE,PTEST,VNEW,.TRUE.,.FALSE.)
GMAX=GMAXSAVE
!bs360
IF (CHRMMT) THEN
   AMIDEFAIL=.FALSE.
   IF (CHECKOMEGAT) THEN
      CALL CHECKOMEGA(COORDS,AMIDEFAIL)
      IF(BHDEBUG) PRINT '(A,L5)','           minimised interpolated geometry: AMIDEFAIL  = ',AMIDEFAIL
   ENDIF
   CHIRALFAIL=.FALSE.
   IF (CHECKCHIRALT) THEN
      CALL CHECKCHIRAL(COORDS,CHIRALFAIL)
      IF(BHDEBUG) PRINT '(A,L5)','           minimised interpolated geometry: CHIRALFAIL = ',CHIRALFAIL
   ENDIF
   IF (CHIRALFAIL.OR.AMIDEFAIL) THEN
      NTRY=NTRY+1
      SFRAC=BHSFRAC+(-1.D0)**NTRY*NCOUNT*0.05D0
      IF ((SFRAC.LE.0).OR.(SFRAC.GE.1.D0)) THEN
         PRINT *,'bhinterp> could not find acceptable interpolated minimum'
         RETURN
      ENDIF
      IF(MOD(NTRY,2).EQ.0) NCOUNT=NCOUNT+1
      PRINT '(A,F8.3)','bhinterp> retry interpolation, SFRAC= ',SFRAC
      GOTO 10
   ENDIF
ENDIF

WRITE(*,*) "sn402: Changed BHINTERP so that it calls ALIGN_DECIDE instead of MINPERMDIST"
WRITE(*,*) "I haven't tested this change and am not certain whether it's sensible." 
WRITE(*,*) "Please check carefully with FASTOVERLAP set that this part of the code is working as you expect, then remove these messages!"
!msb50 changed order CSTART, COORDS because this is only for distance so shouldn't matter
CALL ALIGN_DECIDE(CSTART,COORDS,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
DSTART=DIST
CALL ALIGN_DECIDE(CFINISH,COORDS,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
DFINISH=DIST

ESPREV=0.5D0*BHK*(DSTART**2+DFINISH**2)
! PRINT '(A,5G18.8,I8)','EREAL,ESPREV,DS,DF,SFRAC,ITDONE=',EREAL,ESPREV,DSTART,DFINISH,SFRAC,ITDONE

IF (.NOT.MFLAG) THEN
   PRINT '(2(A,G20.10),A,I8)','bhinterp> WARNING - initial quench failed to converge, energy=',EREAL, &
  &  ' spring energy=',ESPREV,' lbfgs steps=',ITDONE
   PRINT '(A,2G15.5)',' bhinterp> Initial distances: ',DSTART,DFINISH
ENDIF
IF (BHDEBUG) PRINT '(2(A,G20.10),A,I8)',' bhinterp> Initial energy=',EREAL,' spring energy=',ESPREV,' lbfgs steps=',ITDONE
IF (BHDEBUG) PRINT '(A,2G15.5)',' bhinterp> Initial distances: ',DSTART,DFINISH
!
! Check for suspicously small distances, which might converge to an end point if further minimised
!
IF (MIN(DSTART,DFINISH).LT.1.1D0) THEN
   IF (BHDEBUG) PRINT '(A,G20.10)',' bhinterp> one of the initial distances is suspiciously small: reoptimise'
   DSTARTSAVE=DSTART; DFINISHSAVE=DFINISH
   GMAXSAVE=GMAX; GMAX=BHCONV
   NCOUNT=0
   DO WHILE (GMAX.GT.GMAXSAVE) 
      GMAX=GMAX/10.0D0
      CALL MYLBFGS(NCOORDS,MUPDATE,COORDS,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS, &
     &                      .TRUE.,ITDONE,PTEST,VNEW,.TRUE.,.FALSE.)
      WRITE(*,*) "sn402: Changed BHINTERP so that it calls ALIGN_DECIDE instead of MINPERMDIST"
      WRITE(*,*) "I haven't tested this change and am not certain whether it's sensible." 
      WRITE(*,*) "Please check carefully with FASTOVERLAP set that this part of the code is working as you expect, then remove these messages!"
      CALL ALIGN_DECIDE(COORDS,CSTART,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
      DSTART=DIST
      CALL ALIGN_DECIDE(COORDS,CFINISH,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
      DFINISH=DIST
      IF (BHDEBUG) PRINT '(A,G20.10,A,2G15.5)',' bhinterp> for RMS condition ',GMAX,' distances: ',DSTART,DFINISH
      IF ((ABS(DSTARTSAVE-DSTART)*100/DSTARTSAVE.LT.10.0D0).AND. &
          (ABS(DFINISHSAVE-DFINISH)*100/DFINISHSAVE.LT.10.0D0)) EXIT ! guess that we really have converged
      IF ((DSTART.LT.GEOMDIFFTOL).OR.(DFINISH.LT.GEOMDIFFTOL)) EXIT ! we have hit an end point
    
      DSTARTSAVE=DSTART; DFINISHSAVE=DFINISH
   END DO
   GMAX=GMAXSAVE
ENDIF
IF (DSTART.LT.GEOMDIFFTOL) THEN
   IF (NBISECT.EQ.MAXBISECT) THEN
      IF (BHDEBUG) PRINT '(A)',' bhinterp> Bisection limit reached - abandon interpolation'
      CALL FLUSH(6)
      SUCCESS=.FALSE.
      RETURN
   ENDIF
   IF (HITF) THEN
      IF (BHDEBUG) PRINT '(A)',' bhinterp> Both end points hit - abandon interpolation'
      CALL FLUSH(6)
      SUCCESS=.FALSE.
      RETURN
   ENDIF
   NBISECT=NBISECT+1
   SLIMIT=MIN(SFRAC,SLIMIT)
   SFRAC=(SFRAC+FLIMIT)/2.0D0
   IF (BHDEBUG) PRINT '(A,G20.10)',' bhinterp> Initial guess minimised to starting end point - change SFRAC to ',SFRAC
   HITS=.TRUE.
   GOTO 1
ENDIF
IF (DFINISH.LT.GEOMDIFFTOL) THEN
   IF (NBISECT.EQ.MAXBISECT) THEN
      IF (BHDEBUG) PRINT '(A)',' bhinterp> Bisection limit reached - abandon interpolation'
      CALL FLUSH(6)
      SUCCESS=.FALSE.
      RETURN
   ENDIF
   IF (HITS) THEN
      IF (BHDEBUG) PRINT '(A)',' bhinterp> Both end points hit - abandon interpolation'
      CALL FLUSH(6)
      SUCCESS=.FALSE.
      RETURN
   ENDIF
   NBISECT=NBISECT+1
   FLIMIT=MAX(SFRAC,FLIMIT)
   SFRAC=(SFRAC+SLIMIT)/2.0D0
   IF (BHDEBUG) PRINT '(A,G20.10)',' bhinterp> Initial guess minimised to finish end point - change SFRAC to ',SFRAC
   HITF=.TRUE.
   GOTO 1
ENDIF

EPREV=EREAL
COORDSO(1:NCOORDS)=COORDS(1:NCOORDS)
IF (BHINTERPUSELOWEST) THEN
   EBEST=EPREV ! ignore the string energy in determining the "best" intermediate minimum
ELSE
   EBEST=EPREV+ESPREV
ENDIF
BESTCOORDS(1:NCOORDS)=COORDS(1:NCOORDS)
NSUCCESS=0; NFAIL=0; NSUCCESST=0; NFAILT=0
EXITBHRUN=.FALSE.

DO J1=1,BHSTEPS
!
! Perturb current geometry.
!
!  IF (AMBERT .AND. ICINTERPT) THEN
   IF ((AMBERT.OR.NABT) .AND. AMBSTEPT) THEN
      IF (.NOT.ALLOCATED(NICTOT)) THEN
          CALL SETDIHEAM() 
      ENDIF
      IF (AMBPERTT .AND. .NOT.AMBERICT) THEN
   !msb50 - note - Hydrogens may have been swapped around so angles don't 
   ! conform to the ones from internal interpolation in all cases (swapped)

          IF (.NOT.ALLOCATED(PHI_START)) ALLOCATE(PHI_START(nphia +nphih))
          IF (.NOT.ALLOCATED(PHI_FINISH)) ALLOCATE(PHI_FINISH(nphia+nphih))
          CALL CHGETICVAL(CSTART, PHI_START, PHI_START,PHI_START,PHI_START,PHI_START, .TRUE.)
          CALL CHGETICVAL(CFINISH, PHI_FINISH, PHI_FINISH,PHI_FINISH,PHI_FINISH,PHI_FINISH, .TRUE.)
          !DO K = 1, 104
          !   PRINT*,k,"start", PHI_START(k),"FINISH",PHI_FINISH(k)
          !ENDDO
          CALL GET_TWISTABLE(PHI_START, PHI_FINISH)
      ENDIF
     PRINT*, "AMBER TAKESTEP"
     CALL TAKESTEPAMM(COORDS,BHDEBUG, BHSTEPSIZE)
   ELSE IF (CHRMMT) THEN
      CALL TAKESTEPCH(COORDS) ! changed back due to problems with CYPA!
!     DO J2=1,NCOORDS
!        COORDS(J2)=COORDS(J2)+BHSTEPSIZE*(DPRAND()-0.5D0)*2.0D0
!     ENDDO
   ELSE
      DO J2=1,NCOORDS
         COORDS(J2)=COORDS(J2)+BHSTEPSIZE*(DPRAND()-0.5D0)*2.0D0
      ENDDO
   ENDIF
   KNOWE=.FALSE.
   GMAXSAVE=GMAX; GMAX=BHCONV ! mylbfgs now uses GMAX so that we can change this parameter via changep
   CALL MYLBFGS(NCOORDS,MUPDATE,COORDS,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS, &
     &                     .TRUE.,ITDONE,PTEST,VNEW,.TRUE.,.FALSE.)
   GMAX=GMAXSAVE
      WRITE(*,*) "sn402: Changed BHINTERP so that it calls ALIGN_DECIDE instead of MINPERMDIST"
      WRITE(*,*) "I haven't tested this change and am not certain whether it's sensible." 
      WRITE(*,*) "Please check carefully with FASTOVERLAP set that this part of the code is working as you expect, then remove these messages!"
   CALL ALIGN_DECIDE(COORDS,CSTART,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
   DSTART=DIST
   CALL ALIGN_DECIDE(COORDS,CFINISH,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
   DFINISH=DIST
   ESTRING=0.5D0*BHK*(DSTART**2+DFINISH**2)
   ETOTAL=EREAL+ESTRING
!  PRINT '(A,4G20.10)',' bhinterp> E, E+ES, distances=',EREAL,ETOTAL,DSTART,DFINISH
! WRITE(99,*) 864
! WRITE(99,'(A,4G20.10)') 'initial interp, E, E+ES, dists=',EREAL,ESPREV,DSTART,DFINISH
! DO K=1,NATOMS
!    WRITE(99,'(A,3G20.10)') 'LA ',COORDS(3*(K-1)+1),COORDS(3*(K-1)+2),COORDS(3*(K-1)+3)
! ENDDO
!
! Accept/reject step. 
! Don;t allow the interpolated minimum to be either of the endpoints!
!
   IF ((DSTART.LT.GEOMDIFFTOL).OR.(DFINISH.LT.GEOMDIFFTOL)) THEN
      ATEST=.FALSE.
!  ELSEIF ((DSTART.GT.DSTARTSAVE).OR.(DFINISH.GT.DFINISHSAVE)) THEN
!     ATEST=.FALSE.
   ELSEIF (EREAL+ESTRING.LT.EPREV+ESPREV) THEN
      ATEST=.TRUE.
   ELSE
       RANDOM=DPRAND()
       IF (DEXP(-(EREAL+ESTRING-EPREV-ESPREV)/BHTEMP).GT.RANDOM) THEN
          ATEST=.TRUE.
       ELSE
          ATEST=.FALSE.
       ENDIF
   ENDIF
   IF (ATEST) THEN
      NSUCCESS=NSUCCESS+1
      EPREV=EREAL
      ESPREV=ESTRING
      COORDSO(1:NCOORDS)=COORDS(1:NCOORDS)
!
! Record best solution found so far.
!
      IF (BHINTERPUSELOWEST) THEN
         IF (EPREV.LT.EBEST) THEN
            EBEST=EPREV
            BESTCOORDS(1:NCOORDS)=COORDS(1:NCOORDS)
         ENDIF
      ELSE
         IF (EPREV+ESPREV.LT.EBEST) THEN
            EBEST=EPREV+ESPREV
            BESTCOORDS(1:NCOORDS)=COORDS(1:NCOORDS)
         ENDIF
      ENDIF
      IF (BHDEBUG) PRINT '(I6,4(A,G15.5),A,I6,A)',J1,' En= ', EREAL, ' Eo= ',EPREV,  &
  &            ' ETn=', ETOTAL,' ETo=',EPREV+ESPREV,' iter ',ITDONE,' AC'
      IF (BHDEBUG) PRINT '(A,2G15.5)','       distances from start and finish: ',DSTART,DFINISH
   ELSE 
      NFAIL=NFAIL+1
      IF (BHDEBUG) PRINT '(I6,4(A,G15.5),A,I6,A)',J1,' En= ', EREAL, ' Eo= ',EPREV,  &
  &            ' ETn=', ETOTAL,' ETo=',EPREV+ESPREV,' iter ',ITDONE,' RJ'
      IF (BHDEBUG) PRINT '(A,2G15.5)','       distances from start and finish: ',DSTART,DFINISH
      COORDS(1:NCOORDS)=COORDSO(1:NCOORDS)
   ENDIF
!
   IF (BHINTERPUSELOWEST.AND.BHCHECKENERGYT) THEN
      IF (EBEST.LT.BHMAXENERGY .AND. J1.GE.BHSTEPSMIN) EXITBHRUN=.TRUE.
   ENDIF
   IF (EXITBHRUN) EXIT
!
! Adjust step size.
!
   IF (MOD(J1,50).EQ.0) THEN
      P0=(1.D0*NSUCCESS)/(1.D0*(NSUCCESS+NFAIL))
      IF (P0.GT.BHACCREJ) BHSTEPSIZE=BHSTEPSIZE*1.05D0
      IF (P0.LT.BHACCREJ) BHSTEPSIZE=BHSTEPSIZE/1.05D0
      IF (BHDEBUG) PRINT '(A,G20.10,A,G20.10)',' bhinterp> acceptance ratio=',P0,' maximum step size=',BHSTEPSIZE
      NSUCCESST=NSUCCESST+NSUCCESS
      NFAILT=NFAILT+NFAIL
      NSUCCESS=0; NFAIL=0
   ENDIF

!  IF ((DSTART.LT.DENDPT).AND.(DFINISH.LT.DENDPT)) EXIT ! stop if both distances are < end point distance
   
ENDDO ! main loop over BH steps
!
!  Reoptimise best minimum to global RMS tolerance. Here we actually use GMAX not BHCONV!
!
KNOWE=.FALSE.
CALL MYLBFGS(NCOORDS,MUPDATE,BESTCOORDS,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS, &
  &          .TRUE.,ITDONE,PTEST,VNEW,.TRUE.,.FALSE.)
      WRITE(*,*) "sn402: Changed BHINTERP so that it calls ALIGN_DECIDE instead of MINPERMDIST"
      WRITE(*,*) "I haven't tested this change and am not certain whether it's sensible." 
      WRITE(*,*) "Please check carefully with FASTOVERLAP set that this part of the code is working as you expect, then remove these messages!"
CALL ALIGN_DECIDE(BESTCOORDS,CSTART,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
DSTART=DIST
CALL ALIGN_DECIDE(BESTCOORDS,CFINISH,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
DFINISH=DIST
!
! Don't allow the interpolated minimum to be either of the endpoints!
!
IF (BHDEBUG) PRINT '(A,2G20.10)',' bhinterp> DSTART,DFINISH=',DSTART,DFINISH
IF ((DSTART.LT.GEOMDIFFTOL).OR.(DFINISH.LT.GEOMDIFFTOL)) THEN
   SUCCESS=.FALSE.
   IF (BHDEBUG) PRINT '(A)',' bhinterp> Tight quench optimised to an end point'
ELSE
! sf344> AMBER stuff
   IF (NOCISTRANS.AND.(AMBERT.OR.NABT)) THEN
      GOODSTRUCTURE1=.TRUE.
      IF (NOCISTRANSRNA) THEN
         CALL CHECK_CISTRANS_RNA(BESTCOORDS,NATOMS,ZSYM,GOODSTRUCTURE1)
         IF (.NOT.GOODSTRUCTURE1) THEN
            SUCCESS=.FALSE.
            IF (BHDEBUG) PRINT '(A)', ' bhinterp> cis-trans isomerisation in RNA ribose ring detected, abandoning structure'
         ELSE
            SUCCESS=.TRUE.
         ENDIF
      ELSE IF (NOCISTRANSDNA) THEN
         CALL CHECK_CISTRANS_DNA(BESTCOORDS,NATOMS,ZSYM,GOODSTRUCTURE1)
         IF (.NOT.GOODSTRUCTURE1) THEN
            SUCCESS=.FALSE.
            IF (BHDEBUG) PRINT '(A)', ' bhinterp> cis-trans isomerisation in DNA deoxyribose ring detected, abandoning structure'
         ELSE
            SUCCESS=.TRUE.
         ENDIF
      ELSE
          CALL CHECK_CISTRANS_PROTEIN(CSTART,NATOMS,GOODSTRUCTURE1,MINOMEGA,CISARRAY1)
          CALL CHECK_CISTRANS_PROTEIN(BESTCOORDS,NATOMS,GOODSTRUCTURE2,MINOMEGA,CISARRAY2)
          CISARRAY1=CISARRAY1-CISARRAY2
          GOODSTRUCTURE1=.TRUE.
          DO J1=1,NATOMS
                IF(CISARRAY1(J1)/=0) THEN
                  GOODSTRUCTURE1=.FALSE.
                  WRITE(*,'(A,I6)') ' bhinterp> cis-trans isomerisation of a peptide bond detected involving atom ', J1
                END IF
          END DO
          IF(.NOT.GOODSTRUCTURE1) THEN
            SUCCESS=.FALSE.
            WRITE(*,'(A)') ' bhinterp> Cis-trans isomerisation of a peptide bond detected, rejecting'
          ELSE
            SUCCESS=.TRUE.
          END IF
      ENDIF
   ELSE
      SUCCESS=.TRUE.
   ENDIF
ENDIF
IF (.NOT.MFLAG) THEN
   PRINT '(A)',' bhinterp> WARNING - final tight minimisation failed to converge'
   SUCCESS=.FALSE.
ELSE
   IF (BHDEBUG) PRINT '(A,G20.10,A,I8)',' bhinterp> Tight quench: potential energy= ',EREAL,' iterations ',ITDONE
   IF(BHDEBUG) CALL DUMPCOORDS(COORDS, 'bhpoint.xyz', .FALSE.)
ENDIF
IF (EREAL.GT.BHMAXENERGY) THEN
   SUCCESS=.FALSE.
   IF (BHDEBUG) PRINT '(A,G20.10)',' bhinterp> Minimum exceeds energy threshold: ',BHMAXENERGY
ENDIF
INTERPCOORDS(1:NCOORDS)=BESTCOORDS(1:NCOORDS)
FINALENERGY=EREAL
CALL FLUSH(6)
END SUBROUTINE BHINTERP
