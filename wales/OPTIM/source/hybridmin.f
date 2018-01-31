C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C***********************************************************************
C
C  This subroutine is designed to perform a hybrid eigenvector-following
C  optimization for large systems where we do not want to diagonalize the whole Hessian.
C  In contrast to bfgsts.f, in this routine we go DOWNHILL.
C
C***********************************************************************
C
      SUBROUTINE HYBRIDMIN(ITMAX,COORDS,ENERGY,GRAD,MFLAG,RMS,EVALMIN,EVALMAX,VECTS,ITER,POTCALL,PTEST)
      USE COMMONS
      USE KEY
      USE VECCK
      USE ZWK
      USE MODCHARMM
      USE MODHESS
      use porfuncs
      IMPLICIT NONE

      INTEGER J1, J2, INEG, J, ITER, NS, I, K1, NBFGS, ITMAX, ITDONE, J3
      DOUBLE PRECISION GRAD(3*NATOMS),ENERGY,COORDS(3*NATOMS),SUM,AVG,FOBNEW,DUMMY,
     1                 RMS,PSTEP,DPRAND, EVALMAX,EVALMIN,STPMAG,VECSP(3*NATOMS),ESAVE,
     3                 SCALE,RMS2,VECS(3*NATOMS),SOVER,EREAL,XRAT(3*NATOMS),VECTS(3*NATOMS),
     4                 XFOB(3*NATOMS),STEP(3*NATOMS),SAVECOORDS(3*NATOMS),CSTEP(3*NATOMS),
     5                 RAT1,RAT2,TEMP,LP1,LP2,LP,VECL(3*NATOMS), DOTOPT
      DOUBLE PRECISION DELTAS, DELTAT, DIAGEXP(3*NATOMS), DELTASP, TPAR
      CHARACTER(LEN=80) :: FNAME
      LOGICAL MFLAG, AWAY, PVFLAG, RESET, STEST, POTCALL, PTEST, CONVERGED
      COMMON /PVF/ PVFLAG
C
C  Assign enough memory to WORK for a blocksize of 32 to be possible.
C  This is for DSYEVR.
C
      INTEGER IWORK(33*3*NATOMS)
      INTEGER ILWORK, LWORK, NFOUND, ISUPPZ(2*3*NATOMS)
      INTEGER INFO, ISTAT, NEVSSAVE, NEVLSAVE
      DOUBLE PRECISION WORK(33*3*NATOMS), ABSTOL, DIAG(3*NATOMS), DLAMCH, CEIGSAVE
      LOGICAL KNOWE, KNOWG, KNOWH
      COMMON /KNOWN/ KNOWE, KNOWG, KNOWH

      LWORK=33*3*NATOMS
      ILWORK=33*3*NATOMS
      IF (ZSYM(1)(1:1).EQ.'W') THEN
         PRINT*,'HYBRIDMIN procedures have not been programmed for TIP potentials'
         CALL FLUSH(6)
         STOP
      ENDIF
      IF ((NUSEEV.GT.1).AND.(.NOT.NOIT.OR.NOHESS)) THEN
         WRITE(*,'(A)') 'For NUSEEV > 1 you must use a Hessian and NOIT with HYBRIDMIN'
         CALL FLUSH(6)
         STOP
      ENDIF
!
! Check for consistent convergence criteria.
!
      IF (GMAX.NE.CONVR) THEN
         PRINT '(2(A,G20.10),A)','bfgsts> WARNING - GMAX ',GMAX,' is different from CONVR ',CONVR,' - resetting'
         GMAX=MIN(CONVR,GMAX)
         CONVR=MIN(CONVR,GMAX)
      ENDIF
!
!  Reset maximum step sizes in case this isn't the first call to EFOL.
!
      IF (DEBUG) PRINT '(A,G20.10)',' hybridmin> resetting maximum step sizes to ',HMMXSTP
      STPMAX(1:NOPT)=HMMXSTP

      ITER=1
      CONVERGED=.TRUE.
      VECS(1:NOPT)=VECTS(1:NOPT) ! so that VECTS is not overwritten
C
C  If DUMPV is .TRUE. then vector.dump is already open and attached to unit 44.
C  SGI compiler won;t allow us to attach it to another unit.
C  vector.dump could contain multiple dumps for more than the last step, so we 
C  have to make sure we get the results for the last step.
C
      IF (READV.AND.(ITER.EQ.1)) THEN
         IF (DUMPV) THEN
            REWIND(44)
211         READ(44,*,END=111) EVALMIN
            READ(44,*) (VECS(J1),J1=1,NOPT)
            GOTO 211
111         CONTINUE
         ELSE
            IF (FILTH.EQ.0) THEN
               FNAME='vector.dump'
            ELSE
               WRITE(FNAME,'(A)') 'vector.dump.'//TRIM(ADJUSTL(FILTHSTR))
            ENDIF
   
            OPEN(UNIT=45,FILE=FNAME,STATUS='OLD')
20          READ(45,*,END=10) EVALMIN
            READ(45,*) (VECS(J1),J1=1,NOPT)
            GOTO 20
10          CLOSE(45)
            DIAG(1)=EVALMIN
         ENDIF
         WRITE(*,'(A,F20.10)') ' hybridmin> Reaction vector read successfully. Eigenvalue=   ',EVALMIN
         NS=100
         READV=.FALSE. ! needed for runs that do bfgsts then path
      ENDIF
      IF (FREEZE) THEN
         DO J1=1,NATOMS
            IF (FROZEN(J1)) THEN
               IF (VARIABLES) THEN
                  VECS(J1)=0.0D0
               ELSE
                  VECS(3*(J1-1)+1)=0.0D0
                  VECS(3*(J1-1)+2)=0.0D0
                  VECS(3*(J1-1)+3)=0.0D0
               ENDIF
            ENDIF
         ENDDO
      ENDIF

      STEST=.TRUE.
      KNOWG=.FALSE.
      IF (NOHESS) STEST=.FALSE.
      IF (READV.AND.(ITER.EQ.1).AND.(NUSEEV.LE.1)) STEST=.FALSE.
90    IF (PTEST) PRINT*
      NUP=NUSEEV
      IF (PTEST) WRITE(*,11) ITER
11          FORMAT (' hybridmin> Beginning of optimization cycle ', I4,'.',/
     1              ' ---------------------------------------')
      FIXIMAGE=.FALSE.
      IF ((FIXAFTER.GT.0).AND.(ITER.GT.FIXAFTER)) FIXIMAGE=.TRUE.
      IF (CHRMMT) NCHENCALLS = 999 ! make sure non-bond list is updated at the start of each TS search.
      IF (POTCALL) THEN
         IF (PV) THEN
            IF (.NOT.KNOWE) CALL POTENTIAL(COORDS,ENERGY,GRAD,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            PVFLAG=.FALSE.
            CALL PVOPT(COORDS,ENERGY,GRAD)
         ENDIF
         IF ((.NOT.KNOWG).OR.(STEST.AND.(.NOT.KNOWH))) CALL POTENTIAL(COORDS,ENERGY,GRAD,.TRUE.,STEST,RMS,PTEST,.FALSE.)
      ENDIF
      IF ((.NOT.VARIABLES).AND.NOIT.AND.STEST) THEN
         IF (ZSYM(NATOMS).EQ.'SY') THEN
            CALL SHIFTSTOCK(COORDS,NATOMS)
         ELSEIF (ZSYM(NATOMS).EQ.'TH') THEN
            CALL SHIFTHTH(COORDS,NATOMS)
! hk286
         ELSEIF (GTHOMSONT) THEN
            CALL SHIFTHGTH(COORDS,NATOMS)
         ELSE
            CALL SHIFTH(COORDS,.TRUE.,NOPT,NATOMS,ATMASS)
            IF (TWOD) THEN
               DO J1=1,NATOMS
                  J2=3*J1
                  DO J3=1,NATOMS
                     HESS(J2,3*(J3-1)+1)=0.0D0
                     HESS(J2,3*(J3-1)+2)=0.0D0
                     HESS(J2,3*(J3-1)+3)=0.0D0
                     HESS(3*(J3-1)+1,J2)=0.0D0
                     HESS(3*(J3-1)+2,J2)=0.0D0
                     HESS(3*(J3-1)+3,J2)=0.0D0
                  ENDDO
                  HESS(J2,J2)=SHIFTV
               ENDDO
            ENDIF
         ENDIF
      ENDIF
      DO I=1,NOPT
         SAVECOORDS(I)=COORDS(I)
      ENDDO
      IF (POTCALL) THEN
         IF (.NOT.NOHESS) THEN 
            IF (.NOT.NOIT) THEN
               CEIGSAVE=CEIG; CEIG=HMCEIG; NEVSSAVE=NEVS; NEVS=HMNEVS; NEVLSAVE=NEVL; NEVL=HMNEVL ! change common variables
               CALL ITEIG(ITER,COORDS,VECS,EVALMIN,EVALMAX,NS,SOVER,PTEST,VECL,CONVERGED)
               CEIG=CEIGSAVE; NEVS=NEVSSAVE; NEVL=NEVLSAVE ! reset common variables
               DIAG(1)=EVALMIN
            ELSE
               ABSTOL=DLAMCH('Safe  minimum')
               IF (ITER.GT.1) THEN
                  DO J1=1,NOPT
                     VECSP(J1)=VECS(J1)
                  ENDDO
               ENDIF
               CALL DSYEVR('V','I','U',NOPT,HESS,SIZE(HESS,1),0.0D0,1.0D0,1,NUSEEV,ABSTOL,NFOUND,DIAG,
     1                        ZWORK,3*NATOMS,ISUPPZ,WORK,
     2                        LWORK, IWORK, ILWORK, INFO )
               IF (INFO.NE.0) PRINT*,'WARNING - INFO=',INFO,' in DSYEVR'
C              PRINT*,'Optimal and actual values of LWORK=',WORK(1),LWORK
C              PRINT*,'Optimal and actual values of ILWORK=',IWORK(1),ILWORK
               EVALMIN=DIAG(1)
               SOVER=0.0D0
               DO J1=1,NOPT
                  VECS(J1)=ZWORK(J1,1)
                  IF (ITER.GT.1) SOVER=SOVER+VECS(J1)*VECSP(J1)
               ENDDO
               IF (PTEST) THEN
                  IF (ITER.GT.1) THEN
                     WRITE(*,'(A,F15.7,A,F15.7)') ' Smallest eigenvalue=',EVALMIN,' overlap with previous vector=',SOVER
                  ELSE
                     WRITE(*,'(A,F15.7,A,F15.7)') ' Smallest eigenvalue=',EVALMIN
                  ENDIF
               ENDIF
            ENDIF
         ELSE
            CEIGSAVE=CEIG; CEIG=HMCEIG; NEVSSAVE=NEVS; NEVS=HMNEVS ! change common variables
            CALL BEIG(ITER,COORDS,ENERGY,VECS,EVALMIN,NS,SOVER,PTEST,CONVERGED)
            CEIG=CEIGSAVE; NEVS=NEVSSAVE  ! reset common variables
            DIAG(1)=EVALMIN
         ENDIF
      ENDIF
!
!  Return if the eigenvalue goes above a certain threshold.
!
      IF (EVALMIN.GT.HMEVMAX) THEN
         FIXIMAGE=.FALSE.
         MFLAG=.FALSE.
         RETURN
      ENDIF

      CALL VECNORM(VECS,NOPT)
      DUMMY=0.0D0
      DO J1=1,NOPT
         ZWORK(J1,1)=VECS(J1)
         DUMMY=DUMMY+VECS(J1)**2
      ENDDO
!
!  Save XFOB for scaling the maximum allowed step along each eigenvector.
!
      DO I=1,NUSEEV
         XFOB(I)=0.0D0
         DO J = 1,NOPT
            XFOB(I)=XFOB(I)+GRAD(J)*ZWORK(J,I)
         ENDDO
      ENDDO
C
C  Dump the smallest non-zero eigenvalue and eigenvector 
C  in file vector.dump, if required.
C
      IF (DUMPV) THEN
         IF (.NOT.ALLSTEPS) REWIND(44)
         WRITE(44,'(E20.10)') EVALMIN
         WRITE(44,'(3F20.10)') (ZWORK(J1,1),J1=1,NOPT)
         CALL FLUSH(44)
      ENDIF
C
C  Take an eigenvector-following or Page-McIver step downhill along all known eigenvectors.
C  Then do LBFGS minimization in the tangent space if required. 
C
      INEG=0
      DO J1=1,NUSEEV
         IF (DIAG(J1).LT.0.0D0) INEG=INEG+1
      ENDDO
      IF (DEBUG.AND.(INEG.GT.0)) WRITE(*,'(A,I5,A)') ' hybridmin> There are at least ',INEG,' negative Hessian eigenvalues'
C
C  Take a step away from a stationary point along the appropriate
C  Hessian eigenvector. This enables us to start from converged transition states.
C
      AWAY=.FALSE.
      IF (RMS.LT.PUSHCUT) THEN
         IF ((INEG.NE.0).AND.CONVERGED) THEN ! don;t try pushoff if CONVERGED is .FALSE.
            IF ((ITER.EQ.1).AND.(INEG.EQ.1)) THEN
               IF (PTEST) THEN
                  IF (IVEC.GE.0) PRINT '(A)',' hybridmin> Stepping away from minimum along softest mode + direction'
                  IF (IVEC.LT.0) PRINT '(A)',' hybridmin> Stepping away from minimum along softest mode - direction'
               ENDIF
               AWAY=.TRUE.
            ELSE IF (MOD(ITER-1,4).EQ.0) THEN
               DO J1=1,NUSEEV
                  IF (DIAG(J1).LT.0.0D0) THEN
                     WRITE(*,'(A,I6)') ' hybridmin> Stepping away from saddle point along eigenvector ',J1
                     IF (PUSHOFF.EQ.0.0D0) THEN
                        XFOB(J1)=STPMAX(J1)
                     ELSE
                        XFOB(J1)=PUSHOFF
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
         ENDIF
      ENDIF
      IF (HMMETHOD.EQ.'EF') THEN
C
C  EF determination of steps
C
         DO I=NUSEEV,1,-1
            STEP(I)=0.0D0
            LP1=DABS(DIAG(I))/2.0D0
            IF (MASST) THEN
               SUM=0
               DO J2=1,NATOMS
                  SUM=SUM+ATMASS(J2)
               ENDDO
               AVG=SUM/NATOMS
               LP2=1.0D0 + 4.0D0*((XFOB(I)/DIAG(I))**2)/AVG
            ELSE
               LP2=1.0D0 + 4.0D0*(XFOB(I)/DIAG(I))**2
            ENDIF
            LP=LP1*(1.0D0+DSQRT(LP2))
            IF (DEBUG) WRITE(*,'(A,I4,A,4X,F19.10)') ' hybridmin> Mode ',I,' will be searched downhill. Eigenvalue=',DIAG(I)
            STEP(I)=-XFOB(I)/LP
            IF (AWAY.AND.(I.EQ.1)) THEN
               IF (PUSHOFF.NE.0.0D0) THEN
                  STEP(I)=PUSHOFF
               ELSE
                  STEP(I)=STPMAX(1)/1.0D1
               ENDIF
               IF (IVEC.LT.0) STEP(I)=-STEP(I)
            ENDIF
         ENDDO
      ELSE IF (HMMETHOD.EQ.'PM') THEN
         IF (AWAY) THEN
            DO I=1,NUSEEV
               IF (PUSHOFF.NE.0.0D0) THEN
                  STEP(I)=PUSHOFF
               ELSE
                  STEP(I)=STPMAX(1)/1.0D1
               ENDIF
               IF (IVEC.LT.0) STEP(I)=-STEP(I)
            ENDDO
         ELSE
C
C  Step and scaling determined by Page-McIver scheme.
C  STPMAX(1) is dynamically adjusted via a trust radius scheme.
C     
            DELTAT=STPMAX(1)/(100.0D0*RMS*SQRT(1.0D0*NOPT))
            DELTAS=RMS*SQRT(1.0D0*NOPT)*DELTAT/2.0D0
            DO J1=1,NUSEEV
               DIAGEXP(J1)=DEXP(-2*DIAG(J1)*DELTAT)
            ENDDO
            J1=0

666         J1=J1+1
            DELTASP=DELTAS
            TEMP=0.0D0
            DO J2=1,NUSEEV
               TEMP=TEMP+XFOB(J2)**2*DIAGEXP(J2)**J1
            ENDDO 
            DELTAS=DELTAS+DSQRT(TEMP)*DELTAT
!
!  We need to escape from the loop if the integral has effectively converged
!  or if the arc length exceeds the allowed value. 
!  
            IF ((DELTAS.LT.STPMAX(1)).AND.(J1.LT.100000).AND.((DELTAS-DELTASP)/DELTASP.GT.1.0D-10)) GOTO 666

            TPAR=J1*DELTAT
            IF (PTEST) WRITE(*,'(A,G16.6,A,F15.3,A,I6)') ' efol> Estimated arc length=',DELTAS,' for t=',TPAR,
     &                               ' integration steps=',J1

            STEP(1:NOPT)=0.0D0
            DO J1=1,NUSEEV
               STEP(J1)=XFOB(J1)*(DEXP(-DIAG(J1)*TPAR)-1)/DIAG(J1)
            ENDDO
         ENDIF
      ELSE 
         PRINT '(A,A)',' hybridmin> unrecognised hybrid minimisation step method: ',HMMETHOD
         STOP
      ENDIF
C
C  Scale according to step size in ev basis:
C
      STPMAG=MAX(DSQRT(DOTOPT(STEP(1),STEP(1),NUSEEV)),1.0D-10)
      IF (PTEST) WRITE(*,'(A,2F12.6)') ' hybridmin> % of step and gradient along softest mode=',  
     &                              ABS(STEP(1))*100.0D0/STPMAG,ABS(XFOB(1))*100.0D0/(RMS*SQRT(1.0D0*NOPT))
      IF (.NOT.AWAY) THEN 
         DO J1=1,NUSEEV
            SCALE=MIN(STPMAX(J1)/MAX(DABS(STEP(J1)),1D-10),1.0D0)
            IF (DEBUG) PRINT '(A,I8,A,2G20.10)',' hybridmin> scaled and unscaled steps for mode ',J1,' are: ',
     &                                            STEP(J1),SCALE*STEP(J1)
            STEP(J1)=SCALE*STEP(J1)
         ENDDO
      ENDIF
C
C  Convert the steps to the Cartesian rather than the EV basis and put in CSTEP(NOPT+1:2*NOPT)
C
      DO J=1,NOPT
         CSTEP(J)=0.0D0
         DO I=1,NUSEEV
            CSTEP(J)=CSTEP(J)+STEP(I)*ZWORK(J,I)
         ENDDO
      ENDDO
C
C  Put new coordinates into COORDS.
C
      DO J=1,NOPT
         COORDS(J)=SAVECOORDS(J)+CSTEP(J)
      ENDDO
      KNOWE=.FALSE.
      KNOWG=.FALSE.
      KNOWH=.FALSE.

      ESAVE=ENERGY
      CALL POTENTIAL(COORDS,ENERGY,GRAD,.TRUE.,.FALSE.,RMS,PTEST,.FALSE.)
C
C  Adjust maximum step size.
C
      IF (ENERGY.GT.ESAVE+MAXERISE) THEN
         PRINT '(A,G20.10,A,G20.10,A)',' hybridmin> energy increased from ',ESAVE,' to ',ENERGY,' - reversing step'
         DO J1=1,NOPT
            COORDS(J1)=COORDS(J1)-CSTEP(J1)
            STPMAX(J1)=MAX(STPMAX(J1)/2.0D0,MINMAX)
         ENDDO
         ENERGY=ESAVE
      ELSE
         DO J1=1,NUSEEV
            FOBNEW=0.0D0
            DO J2=1,NOPT
               FOBNEW=FOBNEW+GRAD(J2)*ZWORK(J2,J1)
            ENDDO
            XRAT(J1)=MIN(DABS(1.0D0-(FOBNEW-XFOB(J1))/(STEP(J1)*DIAG(J1))),DABS(1.0D0-(-FOBNEW-XFOB(J1))/(STEP(J1)*DIAG(J1))))
!           PRINT '(A,I8,5G15.5)','J1,XRAT, FOBNEW,XFOB,STEP,DIAG=',J1,XRAT(J1),FOBNEW,XFOB(J1),STEP(J1),DIAG(J1)
!           PRINT '(A,I8,5G15.5)','J1,XRAT,-FOBNEW,XFOB,STEP,DIAG=',J1,XRAT(J1),-FOBNEW,XFOB(J1),STEP(J1),DIAG(J1)
            IF (XRAT(J1).GT.TRAD) THEN
               STPMAX(J1)=MAX(STPMAX(J1)/1.1D0,MINMAX)
            ELSE
               STPMAX(J1)=MIN(STPMAX(J1)*1.1D0,MAXMAX)
            ENDIF
         ENDDO
      ENDIF
      CALL DUMPP(COORDS,ENERGY)
      CALL FLUSH(6)
C
C Summarize
C
      IF (PTEST) THEN
         WRITE(*,30)
30       FORMAT(1X,79('-'))
         WRITE(*,40)
40       FORMAT(' Vector      Gradient        Secder       Step          Max step    Trust ratio')
         WRITE(*,30)
         DO I=NUSEEV,1,-1
            WRITE(*,50) I,XFOB(I),DIAG(I),STEP(I),STPMAX(I),XRAT(I)
         ENDDO
50       FORMAT(1X,I4,1X,E15.6,1X,E15.6,1X,E13.6,1X,E12.6,1X,E15.6)
         WRITE(*,30)
      ENDIF
C
C  Tangent space minimization next.
C  Uphill direction is projected out of the step in mylbfgs
C  The next IF block allows for zero tangent space steps in the initial phase
C
      IF ((HMNBFGSMAX1.EQ.0).AND.((1.0D0-DABS(SOVER).GT.0.0001D0))) THEN
         FIXIMAGE=.FALSE.
         ITER=ITER+1
         MFLAG=.FALSE.
         IF (ITER.GT.ITMAX) RETURN
         GOTO 90
      ENDIF

      IF (((1.0D0-DABS(SOVER).GT.0.0001D0).OR.(INEG.NE.0)).OR.(ITER.EQ.1)) THEN
          NBFGS=HMNBFGSMAX1
      ELSE
          NBFGS=HMNBFGSMAX2
      ENDIF

      RESET=.FALSE.
      IF (ITER.EQ.1) RESET=.TRUE.
      IF (CHRMMT.AND.INTMINT) THEN
         CALL MYLBFGS(NINTS,MUPDATE,COORDS,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,NBFGS,
     1                RESET,ITDONE,PTEST,GRAD,.TRUE.,.TRUE.)
      ELSE
         CALL MYLBFGS(NOPT,MUPDATE,COORDS,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,NBFGS,
     1                RESET,ITDONE,PTEST,GRAD,.TRUE.,.TRUE.)
      ENDIF
      RMS2=RMS ! save subspace RMS
      RMS=DSQRT(DOTOPT(GRAD(1),GRAD(1),NOPT))/SQRT(1.0D0*NOPT) ! true RMS
      IF (PTEST) WRITE(*,'(A,F15.7,A,F15.7,A,F15.7)') ' hybridmin> Total RMS gradient=',RMS,
     &         ' subspace gradient=',RMS2,' unscaled second order step length=',STPMAG
      IF (MFLAG) THEN
         IF (((RMS.GT.CONVR).OR.(INEG.GT.0)).OR.(STPMAG.GT.CONVU)) MFLAG=.FALSE.
         IF (MFLAG) RETURN
      ENDIF 

      ITER=ITER+1
      IF (ITER.GT.ITMAX) RETURN
      GOTO 90

      RETURN
      END
