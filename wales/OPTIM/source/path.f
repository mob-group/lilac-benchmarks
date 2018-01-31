C
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
C jmc Note nothing has been done here to fix the unres pathlength coordinate resetting problem...
C
      SUBROUTINE PATH(Q,ENERGY,VNEW,RMS,EVTS,VECS,POTCALL,QPLUS,QMINUS,PTEST,ETS,EPLUS,EMINUS,SLENGTH,DISP,GAMMA,NTILDE,
     1                FRQSTS,FRQSPLUS,FRQSMINUS,ITSTRING,EOFSSTRING,PATHFAILT)
      USE COMMONS
      USE KEY
      USE SYMINF
      USE modcharmm
      USE MODUNRES
      USE MODHESS
      USE MODEFOL
      USE PORFUNCS
      USE AMHGLOBALS, ONLY : NMRES
! hk286
      USE GENRIGID
      
      USE MODCUDALBFGS, ONLY : CUDA_LBFGS_WRAPPER

      IMPLICIT NONE

      DOUBLE PRECISION ENERGY, VNEW(NOPT), RMS, EVTS, VECS(NOPT)
      DOUBLE PRECISION Q(NOPT), QPLUS(NOPT), QMINUS(NOPT)
      DOUBLE PRECISION ETS, EPLUS, EMINUS, SLENGTH, DISP, GAMMA, NTILDE
      DOUBLE PRECISION FRQSTS(NOPT), FRQSPLUS(NOPT), FRQSMINUS(NOPT)
      DOUBLE PRECISION DUMMYA(NOPT), DUMMYB(NOPT)
       
      LOGICAL POTCALL, PTEST
      LOGICAL PATHFAILT
      CHARACTER(LEN=*) ITSTRING,EOFSSTRING

      INTEGER NSTEPPLUS, ITDONE, NSTEPMINUS, J1, J2, NPATHFRAME, NATOMSSAVE, NEWINR, GETUNIT, NEG,
     1        INEG, HORDER, INFO, IVECSAVE, IVEC2SAVE, IPOT, J3, NFPLUS, NFMINUS, RECLEN, ISTAT, NUSE
      INTEGER, PARAMETER :: NFMAX=20000

      DOUBLE PRECISION, ALLOCATABLE :: EOFS(:), PATHLENGTH(:), EOFSFRAMEP(:), EOFSFRAMEM(:)
      DOUBLE PRECISION EVALMAX, RANDOM, QE, PINIT, POPT, PENERGY,
     1                 DIAG(NOPT), EREAL, RMS2, STEP(NOPT), QINIT(NOPT), 
     2                 TEMP, MINIM, SUM2, SUM4, EVPLUS, EVMINUS, 
     4                 SPLUS, SMINUS, STS, STEMP, ETEMP, DUMMY, TIME, TIME0,
     5                 OVEC(3), H1VEC(3), H2VEC(3),
     6                 TEMPA(9*NATOMS), CAPSCOORDS1(18), CAPSCOORDS2(18), 
     7                 DPRAND, PPLUS, PMINUS, lambdats, lambdap, lambdam, distp, distm, RMAT(3,3) !, P(3)
      LOGICAL CONNECTT, DUMPPATH, READPATH, CALCRATES, STOPFIRST, ETEST, NOSHIFTSAVE, NOHESSSAVE
      DOUBLE PRECISION TEMPERATURE, HRED, DIHE, ALLANG, LASTE
      INTEGER NCONNECT
      COMMON /CONN/ STOPFIRST, CONNECTT, NCONNECT, DUMPPATH, READPATH, CALCRATES, TEMPERATURE, HRED
      DOUBLE PRECISION CAPSRHO, EPS2, RAD, HEIGHT
      COMMON /CAPS/ CAPSRHO, EPS2, RAD, HEIGHT

C jmc
      CHARACTER(LEN=80) ITSTRING2
      DOUBLE PRECISION NEWINT(NINTS),TSINT(NINTS)
      DOUBLE PRECISION PEPCOORDS(3*NATOMS), INERTIA(3,3)
      INTEGER K1,K2,KD,NNZ,NINTB
C jmc

C    LOCAL AMH VARIABLES
      INTEGER GLY_COUNT
C      CHARACTER(LEN=5) TARFL

      LOGICAL MFLAG, BFGSTSTSAVE
      LOGICAL PATHT, DRAGT
      COMMON /RUNTYPE/ DRAGT, PATHT, NPATHFRAME
      LOGICAL KNOWE, KNOWG, KNOWH
      COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
      INTEGER NATOMSIMUL, LUNIT
      CHARACTER(LEN=5) ZSYMSAVE
      COMMON /SYS/ ZSYMSAVE
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: Q1, Q2, QW
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: QFRAMEP, QFRAMEM
! hk286
      DOUBLE PRECISION :: XCOORDS(NOPT), XRIGIDCOORDS(DEGFREEDOMS), XCOORDSA(9*NATOMS/2), TMPCOORDS(NOPT)

      PATHFAILT=.FALSE.

      IF (BULKT.AND..NOT.VASP) THEN
         IF (TWOD) THEN
            IF ((PARAM1*PARAM2.EQ.0.0D0).AND.(BULK_BOXVEC(1)*BULK_BOXVEC(2).EQ.0.0D0)) THEN
               PRINT '(A)',' path> ERROR - BULKT is true but a box parameter is zero' 
               STOP
            ENDIF
         ENDIF
      ENDIF
      
      IF (VASP) THEN
         CALL SYSTEM(' cp WAVECAR WAVECAR.trans')
      ENDIF
C  Calls to dumpp dump energies and coordinates on units 1 and 2. Must rewind for path to work.
C  If PRINTPTS is .FALSE. then reduce the I/O to a bare minimum. Use
C  QPLUS, QINIT and QMINUS for the stationary points rather than reading through
C  points saved on disk.
C
      IF (PRINTPTS) THEN 
         REWIND(1)  ! points file
         REWIND(2)  ! energies file
      ENDIF
      BFGSTSTSAVE=BFGSTST
      IVECSAVE=IVEC
      IVEC2SAVE=IVEC2
C
C  Plus side first. 
C
      IVEC=1
      DO J1=1,NOPT
         QINIT(J1)=Q(J1)
      ENDDO
!     IF (DEBUG) PRINT*,'ts points in path:'
!     IF (DEBUG) WRITE(*,'(3G20.10)') (Q(J1),J1=1,NOPT)

      CALL MYCPU_TIME(TIME0,.FALSE.)

!      Do the '+ side' stepoff and minimisation first.

      ! Select which method will be used to compute the path
      PLUSSIDET=.TRUE.
      IF (HYBRIDMINT) THEN
!
! We must have the ts energy somewhere already. Should;t really need the call to potential.
!
         CALL POTENTIAL(Q,ETS,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
         ENERGY=ETS
         POTCALL=.TRUE.
         CALL DUMPP(Q,ENERGY)
         CALL HYBRIDMIN(HMNSTEPS,Q,ENERGY,VNEW,MFLAG,RMS,EVTS,EVALMAX,VECS,ITDONE,POTCALL,PTEST)
         NSTEPPLUS=ITDONE
         IF (.NOT.MFLAG) THEN
            IF (PTEST) PRINT '(A,I8,A)','path> switching to LBFGS minimisation after ',NSTEPPLUS,' hybrid minimisation steps'
            KNOWE=.FALSE.
            KNOWG=.FALSE.
            GMAX=MIN(CONVR,GMAX)
            WRITE(*,'(A,G20.10,A)') ' path> RMS convergence reset to ',GMAX,' for LBFGS and EF phases for consistency'
            IF (CHRMMT.AND.INTMINT) THEN
               CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ELSE
               CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            END IF
            NSTEPPLUS=NSTEPPLUS+ITDONE
         ENDIF
      ELSE IF (BFGSTST) THEN  ! Another pushoff step finder
         BFGSSTEP=.TRUE.
         IF (UNRST) THEN
            CALL INTBFGSTS(NSTEPS,Q,ENERGY,VNEW,MFLAG,RMS,EVTS,EVALMAX,VECS,ITDONE,POTCALL,PTEST)
         ELSE
!           sn402: to compute the pushoff step from the transition state?
            CALL BFGSTS(NSTEPS,Q,ENERGY,VNEW,MFLAG,RMS,EVTS,EVALMAX,VECS,ITDONE,POTCALL,PTEST)
         ENDIF
         ETS=ENERGY
!        IF (DEBUG) PRINT*,'ts step off plus points in path:'
!        IF (DEBUG) WRITE(*,'(3G20.10)') (Q(J1),J1=1,NOPT)
         DO J1=1,NOPT
            STEP(J1)=Q(J1)-QINIT(J1)
         ENDDO
         IF (PUSHOPTT) THEN  ! set REDOTS to zero? maybe not
            PINIT=PUSHOFF
            CALL GOLDEN(VECS,QINIT,ETS,PINIT,POPT,PENERGY)
            Q(1:NOPT)=QINIT(1:NOPT)+POPT*VECS(1:NOPT)
            STEP(1:NOPT)=POPT*VECS(1:NOPT)
            WRITE(6,'(A,3G25.15)') ' path> golden section + pushoff, energy, delta, step: ',
     &                                 PENERGY,PENERGY-ETS,POPT
         ENDIF
         
         IF (RKMIN) RMS=1.0D0
         MFLAG=.FALSE.
         BFGSSTEP=.FALSE.
         BFGSTST=.FALSE.
         IF (BFGSMINT) THEN
            IF (UNRST.OR.(CHRMMT.AND.INTMINT)) THEN
                CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                       .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ELSE
!               sn402: to perform the actual minimisation.
                IF (CUDAT) THEN
                    CALL CUDA_LBFGS_WRAPPER(NOPT,MUPDATE,Q,MFLAG,ENERGY,RMS,BFGSSTEPS,ITDONE,.FALSE.)
                    WRITE(*,'(A,G25.17)') ' path> Final energy is ',ENERGY
                ELSE
!                   PRINT '(A)','path> calling MYBFGSTS + side'
                    CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                       .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
                END IF
            ENDIF
            NSTEPPLUS=ITDONE+1
         ELSE IF (BBRSDMT) THEN
            CALL BBRSDM(Q,MFLAG,ITDONE,ENERGY,RMS,.FALSE.,VNEW,PTEST)
            NSTEPPLUS=ITDONE+1
C
C DAE to switch to BFGS after NSTEPS sd
C
            IF (.NOT.MFLAG) THEN
               KNOWE=.FALSE.
               KNOWG=.FALSE.
               IF (CHRMMT.AND.INTMINT) THEN
                  CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ELSE
                  CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ENDIF
               NSTEPPLUS=NSTEPPLUS+ITDONE
            ENDIF
         ELSE IF (BSMIN.OR.RKMIN) THEN
            NUSE=NSTEPS
            IF (PATHSDSTEPS.GT.0) NUSE=PATHSDSTEPS
            CALL ODESD(NUSE,Q,MFLAG,ITDONE,PTEST)
            NSTEPPLUS=ITDONE+1
C
C DAE to switch to BFGS after NSTEPS sd
C
            IF (.NOT.MFLAG) THEN
               KNOWE=.FALSE.
               KNOWG=.FALSE.
               IF (CHRMMT.AND.INTMINT) THEN
                  CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ELSE
                  CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ENDIF
               NSTEPPLUS=NSTEPPLUS+ITDONE
            ELSE
!
!  ODESD does not return the energy!
!
               CALL POTENTIAL(Q,ENERGY,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            ENDIF
C end DAE
         ELSE
!           IF (INR.LT.6) THEN
!              NEWINR=0 ! so we can use Page-McIver
!           ELSE
!              NEWINR=INR
!           ENDIF
            IF (PMPATHT) THEN
               NEWINR=PMPATHINR
            ELSE
               NEWINR=0 
            ENDIF
C bs360 (29/07/08): ACE does not converge with SD. I will investigate it later.
            IF (ACESOLV) THEN
               NUSE=20
            ELSE
               NUSE=NSTEPS
            ENDIF
            IF ((PATHSDSTEPS.GT.0).AND.(INR.LE.6)) NUSE=PATHSDSTEPS
            NOSHIFTSAVE=NOSHIFT; NOHESSSAVE=NOHESS
            NOSHIFT=.FALSE.; NOHESS=.FALSE.
            CALL EFOL(Q,MFLAG,NUSE,ENERGY,ITDONE,EVPLUS,PTEST,FRQSPLUS,NEWINR)
            NOSHIFT=NOSHIFTSAVE; NOHESS=NOHESSSAVE
            NSTEPPLUS=ITDONE
!
!  Switch to LBFGS if SD did not finish.
!
            IF (.NOT.MFLAG) THEN
               IF (PTEST) PRINT '(A,I8,A)','path> switching to LBFGS minimisation after ',NUSE,' steepest-descent steps'
               KNOWE=.FALSE.
               KNOWG=.FALSE.
               GMAX=MIN(CONVR,GMAX)
               WRITE(*,'(A,G20.10,A)') ' path> RMS convergence reset to ',GMAX,' for LBFGS and EF phases for consistency'
               IF (CHRMMT.AND.INTMINT) THEN
                  CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ELSE
                  CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ENDIF
               NSTEPPLUS=NSTEPPLUS+ITDONE
            ENDIF
         ENDIF
      ELSE  ! Another method for finding the pushoff step
         EFOLSTEP=.TRUE.
         CALL EFOL(Q,MFLAG,1,ENERGY,ITDONE,EVTS,PTEST,FRQSTS,0)
         EFOLSTEP=.FALSE.
         ETS=ENERGY
         VECS(1:NOPT)=0.0D0
         DUMMY=0.0D0
         DO J1=1,NOPT
            STEP(J1)=Q(J1)-QINIT(J1)
            VECS(J1)=STEP(J1)
            DUMMY=DUMMY+VECS(J1)**2
         ENDDO
         DUMMY=SQRT(DUMMY)
         DO J1=1,NOPT
            VECS(J1)=VECS(J1)/DUMMY
         ENDDO
         IF (PUSHOPTT) THEN
            PINIT=PUSHOFF
            CALL GOLDEN(VECS,QINIT,ETS,PINIT,POPT,PENERGY)  
            Q(1:NOPT)=QINIT(1:NOPT)+POPT*VECS(1:NOPT)
            STEP(1:NOPT)=POPT*VECS(1:NOPT)
            WRITE(6,'(A,3G25.15)') ' path> golden section + pushoff, energy, delta, step: ',
     &                                 PENERGY,PENERGY-ETS,POPT
         ENDIF
         KNOWE=.FALSE.
         KNOWG=.FALSE. ! we don`t know the gradient at the point in Q because VNEW isn`t passed from efol?
         IF (BFGSMINT) THEN
C           NOSHIFT=.FALSE.
            IF (CHRMMT.AND.INTMINT) THEN
               CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ELSE
               CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            END IF

            NSTEPPLUS=ITDONE+1
         ELSE IF (BBRSDMT) THEN
            CALL BBRSDM(Q,MFLAG,ITDONE,ENERGY,RMS,.FALSE.,VNEW,PTEST)
            NSTEPPLUS=ITDONE+1
C
C DAE to switch to BFGS after NSTEPS sd
C
            IF (.NOT.MFLAG) THEN
               KNOWE=.FALSE.
               KNOWG=.FALSE.
               IF (CHRMMT.AND.INTMINT) THEN
                  CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ELSE
                  CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ENDIF
               NSTEPPLUS=NSTEPPLUS+ITDONE
            ENDIF
         ELSE IF (BSMIN.OR.RKMIN) THEN
C           NOSHIFT=.FALSE.
            NUSE=NSTEPS
            IF (PATHSDSTEPS.GT.0) NUSE=PATHSDSTEPS
            CALL ODESD(NUSE,Q,MFLAG,ITDONE,PTEST)
            NSTEPPLUS=ITDONE+1
C DAE to switch to BFGS after NSTEPS sd
            IF (.NOT.MFLAG) THEN
               KNOWE=.FALSE.
               KNOWG=.FALSE.
               IF (CHRMMT.AND.INTMINT) THEN
                  CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ELSE
                  CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ENDIF
               NSTEPPLUS=NSTEPPLUS+ITDONE
            ELSE
!
!  ODESD does not return the energy!
!
               CALL POTENTIAL(Q,ENERGY,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            ENDIF
C end DAE
         ELSE
!           IF (INR.LT.6) THEN
!              NEWINR=0 ! so we can use Page-McIver
!           ELSE
!              NEWINR=INR
!           ENDIF
            IF (PMPATHT) THEN
               NEWINR=PMPATHINR
            ELSE
               NEWINR=0
            ENDIF
            NUSE=NSTEPS
            IF ((PATHSDSTEPS.GT.0).AND.(INR.LE.6)) NUSE=PATHSDSTEPS
            NOSHIFTSAVE=NOSHIFT; NOHESSSAVE=NOHESS
            NOSHIFT=.FALSE.; NOHESS=.FALSE.
            CALL EFOL(Q,MFLAG,NUSE,ENERGY,ITDONE,EVPLUS,PTEST,FRQSPLUS,NEWINR)
            NOSHIFT=NOSHIFTSAVE; NOHESS=NOHESSSAVE
            NSTEPPLUS=ITDONE ! bug fix DJW 7/10/08
!
!  Switch to LBFGS if SD did not finish.
!
            IF (.NOT.MFLAG) THEN
               IF (PTEST) PRINT '(A,I8,A)','path> switching to LBFGS minimisation after ',NUSE,' steepest-descent steps'
               KNOWE=.FALSE.
               KNOWG=.FALSE.
               GMAX=MIN(CONVR,GMAX)
               WRITE(*,'(A,G20.10,A)') ' path> RMS convergence reset to ',GMAX,' for LBFGS and EF phases for consistency'
               IF (CHRMMT.AND.INTMINT) THEN
                  CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ELSE
                  CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                         .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
               ENDIF
               NSTEPPLUS=NSTEPPLUS+ITDONE
            ENDIF
         ENDIF
      ENDIF
      PLUSSIDET=.FALSE.

      CALL MYCPU_TIME(TIME,.FALSE.)
      IF (MFLAG) THEN
         WRITE(*,'(A,I20,A,G20.10,2X,A,F11.2)') ' Plus  side of path:    ',NSTEPPLUS,' steps. Energy=',ENERGY,' time=',TIME-TIME0
         CALL FLUSH(6)
      ELSE
         WRITE(*,'(A,I20,A)') ' Plus  side of path failed to converge in ',NSTEPPLUS,' steps'
c        STOP
         PATHFAILT=.TRUE.
         BFGSTST=BFGSTSTSAVE
         IVEC=IVECSAVE
         IVEC2=IVEC2SAVE
C
C jmc Note that before the plus-side minimization above, BFGSTST is set to false, so if we`re doing a connect run (and BFGSTS
C was initially true), the next ts search will mess up, calling efol not intbfgsts.  Reset to BFGSTSTSAVE and also reset IVEC
C and IVEC2 here (as at the end of this subroutine).
C
         IF (ALLOCATED(Q1)) DEALLOCATE(Q1)
         IF (ALLOCATED(Q2)) DEALLOCATE(Q2)
         IF (ALLOCATED(QW)) DEALLOCATE(QW)
         IF (ALLOCATED(QFRAMEP)) DEALLOCATE(QFRAMEP)
         IF (ALLOCATED(QFRAMEM)) DEALLOCATE(QFRAMEM)
         IF (ALLOCATED(EOFS)) DEALLOCATE(EOFS, PATHLENGTH, EOFSFRAMEP, EOFSFRAMEM)
         RETURN
      ENDIF

C Check Hessian index  

      IF ((MFLAG.AND.CHECKINDEX).AND.(BFGSMINT.OR.BFGSTST.OR.BSMIN.OR.RKMIN)) THEN
         IF (NOHESS) THEN
            CALL CHECKIND2(Q,MFLAG,INEG,ENERGY)
         ELSE
C
C  We need the Hessian in CHECKIND. 
C
            IF (BFGSMINT.OR.BSMIN.OR.RKMIN) CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
            CALL CHECKIND(Q,MFLAG,INEG,ENERGY,EVPLUS,EVALMAX,.FALSE.)
         ENDIF
C
C sf344> discard the pathway if one of the minima has a non-zero Hessian index.
C
            IF (INEG>0) THEN
            ! sn402: debug statement
            IF(DEBUG) WRITE(*,*) "path> + direction minimum has non-zero Hessian index. Discarding path."
               PATHFAILT=.TRUE.
               BFGSTST=BFGSTSTSAVE
               IVEC=IVECSAVE
               IVEC2=IVEC2SAVE
               IF (ALLOCATED(Q1)) DEALLOCATE(Q1)
               IF (ALLOCATED(Q2)) DEALLOCATE(Q2)
               IF (ALLOCATED(QW)) DEALLOCATE(QW)
               IF (ALLOCATED(QFRAMEP)) DEALLOCATE(QFRAMEP)
               IF (ALLOCATED(QFRAMEM)) DEALLOCATE(QFRAMEM)
               IF (ALLOCATED(EOFS)) DEALLOCATE(EOFS, PATHLENGTH, EOFSFRAMEP, EOFSFRAMEM)
               RETURN
            END IF
      ENDIF
C
C  Minus side. 
C
!     This time we skip the compute-step stage, and just use the negative of the step that
!     was computed for the + side of the path. (Except for UNRES which is different)
      IF (VASP) THEN
         CALL SYSTEM(' cp WAVECAR.trans WAVECAR')
      ENDIF 

      MINUSSIDET=.TRUE.

      IVEC=-1
      IF (REDOPATH1) THEN
         REDOPATH1=.FALSE.
         REDOPATH2=.TRUE.
      ELSEIF (REDOPATH2) THEN
         REDOPATH1=.TRUE.
         REDOPATH2=.FALSE.
      ENDIF
      DO J1=1,NOPT
         QPLUS(J1)=Q(J1)
         IF (.NOT.UNRST) Q(J1)=QINIT(J1)-STEP(J1)
         IF (HYBRIDMINT) Q(J1)=QINIT(J1)
      ENDDO
      EPLUS=ENERGY
!     IF (DEBUG) PRINT*,'ts step off minus points in path:'
!     IF (DEBUG) WRITE(*,'(3G20.10)') (Q(J1),J1=1,NOPT)

      IF (PUSHOPTT) THEN 
         PINIT=-PUSHOFF
         CALL GOLDEN(VECS,QINIT,ETS,PINIT,POPT,PENERGY)
         Q(1:NOPT)=QINIT(1:NOPT)+POPT*VECS(1:NOPT)
         STEP(1:NOPT)=POPT*VECS(1:NOPT)
         WRITE(6,'(A,3G25.15)') ' path> golden section - pushoff, energy, delta, step: ',
     &                              PENERGY,PENERGY-ETS,POPT
      ENDIF

      IF (UNRST) THEN ! jmc new intstep stuff 
         DO J1=1,nres
            c(1,J1)=QINIT(6*(J1-1)+1)
            c(2,J1)=QINIT(6*(J1-1)+2)
            c(3,J1)=QINIT(6*(J1-1)+3)
            c(1,J1+nres)=QINIT(6*(J1-1)+4)
            c(2,J1+nres)=QINIT(6*(J1-1)+5)
            c(3,J1+nres)=QINIT(6*(J1-1)+6)
         END DO
         CALL UPDATEDC
         CALL int_from_cart(.true.,.false.)
         CALL geom_to_var(NINTS,TSINT)
         NEWINT=TSINT-INTSTEP
         CALL var_to_geom(NINTS,NEWINT)
         CALL chainbuild
         DO J1=1,nres
            Q(6*(J1-1)+1)=c(1,J1)
            Q(6*(J1-1)+2)=c(2,J1)
            Q(6*(J1-1)+3)=c(3,J1)
            Q(6*(J1-1)+4)=c(1,J1+nres)
            Q(6*(J1-1)+5)=c(2,J1+nres)
            Q(6*(J1-1)+6)=c(3,J1+nres)
         END DO
      ENDIF ! jmc end new stuff

      KNOWE=.FALSE.
      KNOWG=.FALSE.

      CALL MYCPU_TIME(TIME0,.FALSE.)
      
      ! Do the minimisation exactly as before.
      ! HYBRIDMINT BFGSMIT BBRSDMT BSMIN RKMIN;  ELSE
      IF (HYBRIDMINT) THEN
         POTCALL=.TRUE.
         ENERGY=ETS
         CALL HYBRIDMIN(HMNSTEPS,Q,ENERGY,VNEW,MFLAG,RMS,EVTS,EVALMAX,VECS,ITDONE,POTCALL,PTEST)
         NSTEPMINUS=ITDONE
         IF (.NOT.MFLAG) THEN
            IF (PTEST) PRINT '(A,I8,A)','path> switching to LBFGS minimisation after ',NSTEPMINUS,' hybrid minimisation steps'
            KNOWE=.FALSE.
            KNOWG=.FALSE.
            GMAX=MIN(CONVR,GMAX)
            WRITE(*,'(A,G20.10,A)') ' path> RMS convergence reset to ',GMAX,' for LBFGS and EF phases for consistency'
            IF (CHRMMT.AND.INTMINT) THEN
               CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ELSE
               CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ENDIF
            NSTEPMINUS=NSTEPMINUS+ITDONE
         ENDIF
      ELSE IF (BFGSMINT) THEN
         IF (UNRST.OR.(CHRMMT.AND.INTMINT)) THEN
            CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                   .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
         ELSE
            IF (CUDAT) THEN
                CALL CUDA_LBFGS_WRAPPER(NOPT,MUPDATE,Q,MFLAG,ENERGY,RMS,BFGSSTEPS,ITDONE,.FALSE.)
                WRITE(*,'(A,G25.17)') ' path> Final energy is ',ENERGY
            ELSE
!                   PRINT '(A)','path> calling MYBFGSTS - side'
                CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                   .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            END IF
         ENDIF
         NSTEPMINUS=ITDONE+1
      ELSE IF (BBRSDMT) THEN
         CALL BBRSDM(Q,MFLAG,ITDONE,ENERGY,RMS,.FALSE.,VNEW,PTEST)
         NSTEPMINUS=ITDONE+1
C
C DAE to switch to BFGS after NSTEPS sd
C
         IF (.NOT.MFLAG) THEN
            KNOWE=.FALSE.
            KNOWG=.FALSE.
            IF (CHRMMT.AND.INTMINT) THEN
               CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ELSE
               CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ENDIF
            NSTEPMINUS=NSTEPMINUS+ITDONE
         ENDIF
      ELSE IF (BSMIN.OR.RKMIN) THEN
         NUSE=NSTEPS
         IF (PATHSDSTEPS.GT.0) NUSE=PATHSDSTEPS
         CALL ODESD(NUSE,Q,MFLAG,ITDONE,PTEST)
         NSTEPMINUS=ITDONE+1
C DAE to switch to BFGS after NSTEPS sd
         IF (.NOT.MFLAG) THEN
            KNOWE=.FALSE.
            KNOWG=.FALSE.
            IF (CHRMMT.AND.INTMINT) THEN
               CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ELSE
               CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ENDIF
            NSTEPMINUS=NSTEPMINUS+ITDONE
         ELSE
!
!  ODESD does not return the energy!
!
            CALL POTENTIAL(Q,ENERGY,VNEW,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
         ENDIF
C end DAE
      ELSE
!        IF (INR.LT.6) THEN
!           NEWINR=0 ! so we can use Page-McIver
!        ELSE
!           NEWINR=INR
!        ENDIF
         IF (PMPATHT) THEN
            NEWINR=PMPATHINR
         ELSE
            NEWINR=0
         ENDIF

C bs360 (29/07/08): ACE does not converge with SD. I will investigate it later.
         IF (ACESOLV) THEN
            NUSE=20
         ELSE
            NUSE=NSTEPS
         ENDIF
         IF ((PATHSDSTEPS.GT.0).AND.(INR.LE.6)) NUSE=PATHSDSTEPS
         NOSHIFTSAVE=NOSHIFT; NOHESSSAVE=NOHESS
         NOSHIFT=.FALSE.; NOHESS=.FALSE.
         CALL EFOL(Q,MFLAG,NUSE,ENERGY,ITDONE,EVMINUS,PTEST,FRQSMINUS,NEWINR)
         NOSHIFT=NOSHIFTSAVE; NOHESS=NOHESSSAVE
         NSTEPMINUS=ITDONE
!
!  Switch to LBFGS if SD did not finish.
!
         IF (.NOT.MFLAG) THEN
            IF (PTEST) PRINT '(A,I8,A)','path> switching to LBFGS minimisation after ',NUSE,' steepest-descent steps'
            KNOWE=.FALSE.
            KNOWG=.FALSE.
            GMAX=MIN(CONVR,GMAX)
            WRITE(*,'(A,G20.10,A)') ' path> RMS convergence set to ',GMAX,' for consistency'
            IF (CHRMMT.AND.INTMINT) THEN
               CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ELSE
               CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,PTEST,VNEW,.FALSE.,.FALSE.)
            ENDIF
            NSTEPMINUS=NSTEPMINUS+ITDONE
         ENDIF
      ENDIF
      MINUSSIDET=.FALSE.

      CALL MYCPU_TIME(TIME,.FALSE.)

      IF (MFLAG) THEN
         WRITE(*,'(A,I20,A,G20.10,2X,A,F11.2)') ' Minus side of path:    ',NSTEPMINUS,' steps. Energy=',ENERGY,' time=',TIME-TIME0
         CALL FLUSH(6)
      ELSE
         WRITE(*,'(A,I20,A)') ' Minus side of path failed to converge in ',NSTEPMINUS,' steps'
c        STOP
         PATHFAILT=.TRUE.
         BFGSTST=BFGSTSTSAVE
         IVEC=IVECSAVE
         IVEC2=IVEC2SAVE
C jmc Note that before the plus-side minimization above, BFGSTST is set to false, so if we`re doing a connect run (and BFGSTS 
C was initially true), the next ts search will mess up, calling efol not intbfgsts.  Reset to BFGSTSTSAVE and also reset IVEC
C and IVEC2 here (as at the end of this subroutine).
         IF (ALLOCATED(Q1)) DEALLOCATE(Q1)
         IF (ALLOCATED(Q2)) DEALLOCATE(Q2)
         IF (ALLOCATED(QW)) DEALLOCATE(QW)
         IF (ALLOCATED(QFRAMEP)) DEALLOCATE(QFRAMEP)
         IF (ALLOCATED(QFRAMEM)) DEALLOCATE(QFRAMEM)
         IF (ALLOCATED(EOFS)) DEALLOCATE(EOFS, PATHLENGTH, EOFSFRAMEP, EOFSFRAMEM)
         RETURN
      ENDIF


C
C  Check Hessian index 
C
      IF ((MFLAG.AND.CHECKINDEX).AND.(BFGSMINT.OR.BFGSTST.OR.BSMIN.OR.RKMIN)) THEN
         IF (NOHESS) THEN
            CALL CHECKIND2(Q,MFLAG,INEG,ENERGY)
         ELSE
C
C  We need the Hessian in CHECKIND. 
C
            IF (BFGSMINT.OR.BSMIN.OR.RKMIN) CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
            CALL CHECKIND(Q,MFLAG,INEG,ENERGY,EVMINUS,EVALMAX,.FALSE.)
         ENDIF
C
C sf344> discard the pathway if one of the minima has a non-zero Hessian index.
C
            IF (INEG>0) THEN
            ! sn402: debug statement
            IF(DEBUG) WRITE(*,*) "path> The minus-side minimum has non-zero Hessian index. Discarding path."
               PATHFAILT=.TRUE.
               BFGSTST=BFGSTSTSAVE
               IVEC=IVECSAVE
               IVEC2=IVEC2SAVE
               IF (ALLOCATED(Q1)) DEALLOCATE(Q1)
               IF (ALLOCATED(Q2)) DEALLOCATE(Q2)
               IF (ALLOCATED(QW)) DEALLOCATE(QW)
               IF (ALLOCATED(QFRAMEP)) DEALLOCATE(QFRAMEP)
               IF (ALLOCATED(QFRAMEM)) DEALLOCATE(QFRAMEM)
               IF (ALLOCATED(EOFS)) DEALLOCATE(EOFS, PATHLENGTH, EOFSFRAMEP, EOFSFRAMEM)
               RETURN
            END IF
      ENDIF

      DO J1=1,NOPT
         QMINUS(J1)=Q(J1)
      ENDDO
      EMINUS=ENERGY
C
C  The total number of energies and coordinates is NSTEPPLUS + NSTEPMINUS + 1 for the transition state.
C  The rest of this subroutine is post-processing 
C  Minimise IO if PRINTPTS is .FALSE.
C
      IF (ZSYMSAVE(1:1).EQ.'W') THEN  !  WCOMMENT
         ALLOCATE(QFRAMEP(9*(NATOMS/2),NFMAX),QFRAMEM(9*(NATOMS/2),NFMAX))
      ELSE
         ALLOCATE(QFRAMEP(NOPT,NFMAX),QFRAMEM(NOPT,NFMAX))
!        ALLOCATE(QFRAMEP(3*NATOMS,NFMAX),QFRAMEM(3*NATOMS,NFMAX))
      ENDIF
      IF (.NOT.PRINTPTS) THEN 
         ALLOCATE(EOFS(3), PATHLENGTH(3), EOFSFRAMEP(3), EOFSFRAMEM(3))
         NSTEPPLUS=1
         NSTEPMINUS=1
         NFPLUS=1
         NFMINUS=1
         EOFS(1)=EPLUS
         EOFS(2)=ETS
         EOFS(3)=EMINUS
         EOFSFRAMEP(1)=EPLUS
         EOFSFRAMEM(1)=EMINUS
         PATHLENGTH(1)=0.0D0
         DUMMY=0.0D0
C        PRINT*,'WARNING - S,N,gamma are not calculated from Cartesian coordinates here'
         DO J2=1,NOPT
            DUMMY=DUMMY+(QINIT(J2)-QPLUS(J2))**2
         ENDDO
         PATHLENGTH(2)=SQRT(DUMMY)
         DUMMY=0.0D0
         DO J2=1,NOPT
            DUMMY=DUMMY+(QINIT(J2)-QMINUS(J2))**2
         ENDDO
         PATHLENGTH(3)=SQRT(DUMMY)+PATHLENGTH(2)
         DO J2=1,NOPT
            QFRAMEP(J2,1)=QPLUS(J2)
            QFRAMEM(J2,1)=QMINUS(J2)
         ENDDO
         GOTO 555
      ENDIF

        ! sn402: presumably the energies and coordinates have been written to these units
        ! in the minimiser routine - whatever that was.
      REWIND(1)
      REWIND(2)
      J1=0
      DO 
         READ(2,*,END=975) DUMMY
         J1=J1+1
      ENDDO
975   IF (DEBUG) PRINT '(A,I6)',' path> number of entries in EofS file=',J1
      REWIND(2)
      ALLOCATE(EOFS(J1), PATHLENGTH(J1), EOFSFRAMEP(J1), EOFSFRAMEM(J1))
      READ(2,*) EOFS(NSTEPPLUS+1)
      IF (ZSYMSAVE(1:1).EQ.'W') THEN  !  WCOMMENT
         ALLOCATE(Q1(9*(NATOMS/2)),Q2(9*(NATOMS/2)),QW(9*(NATOMS/2)))
C        READ(1,*) (Q1(J1),J1=1,9*NATOMS)
         READ(1,*) (Q1(J1),J1=1,9*(NATOMS/2))
      ELSE
         ALLOCATE(Q1(NOPT),Q2(NOPT))
         READ(1,*) (Q1(J1),J1=1,NOPT)
      ENDIF
      PATHLENGTH(NSTEPPLUS+1)=0.0D0
C     
C  The number of frames for each side of the path, NPATHFRAME, is now treated 
C  in an average way. We must dump the three stationary points at the very least,
C  and other frames on the two paths are then dumped with a probability proportional
C  to NPATHFRAME. If NPATHFRAME is > NFMAX then use NFMAX in its place, otherwise
C  we don;t have enough storage declared for the frames.
C
C  The new variable FRAMEDIST can be used to exclude frames with configurations that
C  are separated by less than FRAMEDIST distance units.
C
!     PPLUS=MIN(MIN(MIN(NPATHFRAME,NSTEPPLUS),NFMAX)*1.0D0/(1.0D0*(NSTEPPLUS-1)),1.0D0)
!     PMINUS=MIN(MIN(MIN(NPATHFRAME,NSTEPMINUS),NFMAX)*1.0D0/(1.0D0*(NSTEPMINUS-1)),1.0D0)

      if (NSTEPPLUS.lt.NPATHFRAME .or. NSTEPPLUS.eq.1) then
         PPLUS = 1.0d0
      else
         PPLUS = MIN(MIN(NPATHFRAME,NFMAX)*1.0D0/(1.0D0*(NSTEPPLUS-1)),1.0D0)
      endif

      if (NSTEPMINUS.lt.NPATHFRAME .or. NSTEPMINUS.eq.1) then
         PMINUS = 1.0d0
      else
         PMINUS=MIN(MIN(NPATHFRAME,NFMAX)*1.0D0/(1.0D0*(NSTEPMINUS-1)),1.0D0)
      endif

      IF (PTEST) WRITE(*,'(A,F10.4,A,F10.4,A)') ' Frames will be dumped to points.path.xyz with probability ',
     1                          PPLUS,'/',PMINUS,' steps on the plus/minus sides'
      NFPLUS=0
      DO J1=1,NSTEPPLUS
         READ(2,*) EOFS(NSTEPPLUS+1-J1)
         IF (ZSYMSAVE(1:1).EQ.'W') THEN  
C           READ(1,*) (Q2(J2),J2=1,9*NATOMS)  !  WCOMMENT
            READ(1,*) (Q2(J2),J2=1,9*(NATOMS/2))
         ELSE
            READ(1,*) (Q2(J2),J2=1,NOPT)
         ENDIF

         ! Only print the frame if ETEST is TRUE.
         ETEST=.FALSE.
         IF ((J1.EQ.1).OR.(J1.EQ.NSTEPPLUS)) THEN ! always take the end points
            ETEST=.TRUE.
         ELSE
            IF (ABS(EOFS(NSTEPPLUS+1-J1)-LASTE).GE.FRAMEEDIFF) ETEST=.TRUE. ! energy must have changed enough
         ENDIF
         IF (ETEST) LASTE=EOFS(NSTEPPLUS+1-J1)
         ! This is where we test whether to print the frame:
         IF (((DPRAND().LE.PPLUS).OR.(J1.EQ.NSTEPPLUS)).AND.(NFPLUS.LT.NFMAX).AND.ETEST) THEN
            IF ((J1.EQ.NSTEPPLUS).OR.(NFPLUS.LT.NFMAX-1)) THEN ! save room for the last endpoint
               NFPLUS=NFPLUS+1

               ! Save the energy of the frame
               EOFSFRAMEP(NFPLUS)=EOFS(NSTEPPLUS+1-J1)
               LASTE=EOFS(NSTEPPLUS+1-J1)

               ! Save the coordinates of the frame
               IF (ZSYMSAVE(1:1).EQ.'W') THEN
C                 DO J2=1,9*NATOMS  !  WCOMMENT
                  DO J2=1,9*(NATOMS/2)
                     QFRAMEP(J2,NFPLUS)=Q2(J2)
                  ENDDO
               ELSE IF (ZSYMSAVE(1:2).EQ.'CD') THEN
                   DO J2=1,NATOMS/2
                      CALL CAPSIDIO(Q2(3*(J2-1)+1),Q2(3*(J2-1)+2),Q2(3*(J2-1)+3),
     1                              Q2(3*(NATOMS/2+J2-1)+1),Q2(3*(NATOMS/2+J2-1)+2),Q2(3*(NATOMS/2+J2-1)+3),
     2                              CAPSCOORDS2,RAD,HEIGHT)
                      DO J3=1,18
                         QFRAMEP(18*(J2-1)+J3,NFPLUS)=CAPSCOORDS2(J3)
                      ENDDO
                  ENDDO
               ELSE
                  DO J2=1,NOPT
                     QFRAMEP(J2,NFPLUS)=Q2(J2)
                  ENDDO
               ENDIF
!              PRINT '(A,I6,A,I6,A,G20.10)','dumping plus frame ',J1,' NFPLUS=',NFPLUS,' energy=',EOFSFRAMEP(NFPLUS)
            ENDIF
         ENDIF
         TEMP=0.0D0
         ! sn402: Calculate the distance to the previously entered frame (as TEMP)
         IF (RIGIDINIT .AND. .NOT. ATOMRIGIDCOORDT) THEN
            CALL RB_DISTANCE(TEMP,Q1(1:DEGFREEDOMS),Q2(1:DEGFREEDOMS),DUMMYA,DUMMYB,.FALSE.)
            TEMP = TEMP**2
         ELSE IF (BULKT.AND..NOT.VASP) THEN
            DO J2=1,NATOMS
               TEMP=TEMP+MINIM(Q2(3*(J2-1)+1),Q1(3*(J2-1)+1),PARAM1)**2
     1                  +MINIM(Q2(3*(J2-1)+2),Q1(3*(J2-1)+2),PARAM2)**2
               IF (.NOT.TWOD) TEMP=TEMP+MINIM(Q2(3*(J2-1)+3),Q1(3*(J2-1)+3),PARAM3)**2
            ENDDO
         ELSE
            IF (ZSYMSAVE(1:1).EQ.'W') THEN
C              DO J2=1,9*NATOMS  !  WCOMMENT
               DO J2=1,9*(NATOMS/2)
                  TEMP=TEMP+(Q2(J2)-Q1(J2))**2
               ENDDO
            ELSE IF (ZSYMSAVE.EQ.'CD') THEN
               DO J2=1,NATOMS/2
                  CALL CAPSIDIO(Q2(3*(J2-1)+1),Q2(3*(J2-1)+2),Q2(3*(J2-1)+3),
     1                          Q2(3*(NATOMS/2+J2-1)+1),Q2(3*(NATOMS/2+J2-1)+2),Q2(3*(NATOMS/2+J2-1)+3),
     2                          CAPSCOORDS2,RAD,HEIGHT)
                  CALL CAPSIDIO(Q1(3*(J2-1)+1),Q1(3*(J2-1)+2),Q1(3*(J2-1)+3),
     1                          Q1(3*(NATOMS/2+J2-1)+1),Q1(3*(NATOMS/2+J2-1)+2),Q1(3*(NATOMS/2+J2-1)+3),
     2                          CAPSCOORDS1,RAD,HEIGHT)
                  DO J3=1,18
                      TEMP=TEMP+(CAPSCOORDS1(J3)-CAPSCOORDS2(J3))**2
                  ENDDO
               ENDDO
            ELSE IF((PYT.OR.PYGPERIODICT.OR.PYBINARYT.OR.MULTISITEPYT).AND.UNIAXT) THEN
C           calculate path lengths wrt. Cartesian coordinate of the centre and two coordinates along the symmetry axis only
               CALL UNIAXGETPATHLENGTH(Q1,Q2,TEMP)
            ELSE
               DO J2=1,NOPT
                  TEMP=TEMP+(Q2(J2)-Q1(J2))**2
               ENDDO
            ENDIF
         ENDIF
         ! Compute the total pathlength from the TS to this frame
         PATHLENGTH(NSTEPPLUS+1-J1)=PATHLENGTH(NSTEPPLUS+2-J1)-SQRT(TEMP)
         IF (ZSYMSAVE(1:1).EQ.'W') THEN
C           DO J2=1,9*NATOMS  !  WCOMMENT
            DO J2=1,9*(NATOMS/2)
               Q1(J2)=Q2(J2)
            ENDDO
         ELSE
            ! Copy the current frame to the "previous frame" array
            DO J2=1,NOPT
               Q1(J2)=Q2(J2)
            ENDDO
         ENDIF
      ENDDO  ! Loop over + direction frames
      IF (PTEST) WRITE(*,'(A,I6,A,I6,A)') ' Transition state will be frame number ',NFPLUS+1

      ! Now set up for the - side of the path
      IF (ZSYMSAVE(1:1).EQ.'W') THEN
C        DO J1=1,NATOMS  !  WCOMMENT
         DO J1=1,NATOMS/2
            CALL CONVERT(QINIT(3*(J1-1)+1),QINIT(3*(J1-1)+2),QINIT(3*(J1-1)+3),
C    1                   QINIT(3*(NATOMS+J1-1)+1),QINIT(3*(NATOMS+J1-1)+2),QINIT(3*(NATOMS+J1-1)+3),
     1                   QINIT(3*(NATOMS/2+J1-1)+1),QINIT(3*(NATOMS/2+J1-1)+2),QINIT(3*(NATOMS/2+J1-1)+3),
     2                   OVEC,H1VEC,H2VEC)
            Q1(9*(J1-1)+1)=OVEC(1)
            Q1(9*(J1-1)+2)=OVEC(2)
            Q1(9*(J1-1)+3)=OVEC(3)
            Q1(9*(J1-1)+4)=H1VEC(1)
            Q1(9*(J1-1)+5)=H1VEC(2)
            Q1(9*(J1-1)+6)=H1VEC(3)
            Q1(9*(J1-1)+7)=H2VEC(1)
            Q1(9*(J1-1)+8)=H2VEC(2)
            Q1(9*(J1-1)+9)=H2VEC(3)
         ENDDO
      ELSE
         DO J1=1,NOPT
            Q1(J1)=QINIT(J1)
         ENDDO
      ENDIF

      NFMINUS=0
      DO J1=1,NSTEPMINUS
         READ(2,*) EOFS(NSTEPPLUS+1+J1)
         IF (ZSYMSAVE(1:1).EQ.'W') THEN
C           READ(1,*) (Q2(J2),J2=1,9*NATOMS)  !  WCOMMENT
            READ(1,*) (Q2(J2),J2=1,9*(NATOMS/2))
         ELSE
            READ(1,*) (Q2(J2),J2=1,NOPT)
         ENDIF

         ETEST=.FALSE.
         IF ((J1.EQ.1).OR.(J1.EQ.NSTEPMINUS)) THEN ! always take the end points
            ETEST=.TRUE.
         ELSE
            IF (ABS(EOFS(NSTEPPLUS+1+J1)-LASTE).GE.FRAMEEDIFF) ETEST=.TRUE. ! energy must have changed enough
         ENDIF
         IF (ETEST) LASTE=EOFS(NSTEPPLUS+1+J1)
         RANDOM=DPRAND()
         IF (((RANDOM.LE.PMINUS).OR.(J1.EQ.NSTEPMINUS)).AND.(NFMINUS.LE.NFMAX).AND.ETEST) THEN
            IF ((J1.EQ.NSTEPMINUS).OR.(NFMINUS.LT.NFMAX-1)) THEN ! save a space for the stationary point at the end
               NFMINUS=NFMINUS+1
               EOFSFRAMEM(NFMINUS)=EOFS(NSTEPPLUS+1+J1)
               LASTE=EOFS(NSTEPPLUS+1+J1)
               IF (ZSYMSAVE(1:1).EQ.'W') THEN
C                 DO J2=1,9*NATOMS  !  WCOMMENT
                     DO J2=1,9*(NATOMS/2)
                     QFRAMEM(J2,NFMINUS)=Q2(J2)
                  ENDDO
               ELSE IF (ZSYMSAVE.EQ.'CD') THEN
                   DO J2=1,NATOMS/2
                      CALL CAPSIDIO(Q2(3*(J2-1)+1),Q2(3*(J2-1)+2),Q2(3*(J2-1)+3),
     1                              Q2(3*(NATOMS/2+J2-1)+1),Q2(3*(NATOMS/2+J2-1)+2),Q2(3*(NATOMS/2+J2-1)+3),
     2                              CAPSCOORDS2,RAD,HEIGHT)
                      DO J3=1,18
                         QFRAMEM(18*(J2-1)+J3,NFMINUS)=CAPSCOORDS2(J3)
                      ENDDO
                  ENDDO
               ELSE
                  DO J2=1,NOPT
                     QFRAMEM(J2,NFMINUS)=Q2(J2)
                  ENDDO
               ENDIF
!              PRINT '(A,I6,A,I6,A,G20.10)','dumping minus frame ',J1,' NFMINUS=',NFMINUS,' energy=',EOFSFRAMEM(NFMINUS)

            ENDIF
         ENDIF
         TEMP=0.0D0
         IF (RIGIDINIT .AND. .NOT. ATOMRIGIDCOORDT) THEN
            CALL RB_DISTANCE(TEMP,Q1(1:DEGFREEDOMS),Q2(1:DEGFREEDOMS),DUMMYA,DUMMYB,.FALSE.)
            TEMP = TEMP**2
         ELSE IF (BULKT.AND..NOT.VASP) THEN
            DO J2=1,NATOMS
               TEMP=TEMP+MINIM(Q2(3*(J2-1)+1),Q1(3*(J2-1)+1),PARAM1)**2
     1                  +MINIM(Q2(3*(J2-1)+2),Q1(3*(J2-1)+2),PARAM2)**2
               IF (.NOT.TWOD) TEMP=TEMP+MINIM(Q2(3*(J2-1)+3),Q1(3*(J2-1)+3),PARAM3)**2
            ENDDO
         ELSE
            IF (ZSYMSAVE(1:1).EQ.'W') THEN
C              DO J2=1,9*NATOMS  !  WCOMMENT
               DO J2=1,9*(NATOMS/2)
                  TEMP=TEMP+(Q2(J2)-Q1(J2))**2
               ENDDO
            ELSE IF (ZSYMSAVE.EQ.'CD') THEN
               DO J2=1,NATOMS/2
                  CALL CAPSIDIO(Q2(3*(J2-1)+1),Q2(3*(J2-1)+2),Q2(3*(J2-1)+3),
     1                          Q2(3*(NATOMS/2+J2-1)+1),Q2(3*(NATOMS/2+J2-1)+2),Q2(3*(NATOMS/2+J2-1)+3),
     2                          CAPSCOORDS2,RAD,HEIGHT)
                  CALL CAPSIDIO(Q1(3*(J2-1)+1),Q1(3*(J2-1)+2),Q1(3*(J2-1)+3),
     1                          Q1(3*(NATOMS/2+J2-1)+1),Q1(3*(NATOMS/2+J2-1)+2),Q1(3*(NATOMS/2+J2-1)+3),
     2                          CAPSCOORDS1,RAD,HEIGHT)
                  DO J3=1,18
                      TEMP=TEMP+(CAPSCOORDS1(J3)-CAPSCOORDS2(J3))**2
                  ENDDO
               ENDDO
            ELSE IF((PYT.OR.PYGPERIODICT.OR.PYBINARYT.OR.MULTISITEPYT).AND.UNIAXT) THEN
C           calculate path lengths wrt. Cartesian coordinate of the centre and two coordinates along the symmetry axis only
               CALL UNIAXGETPATHLENGTH(Q1,Q2,TEMP)
            ELSE
               DO J2=1,NOPT
                  TEMP=TEMP+(Q2(J2)-Q1(J2))**2
               ENDDO
            ENDIF
         ENDIF
         PATHLENGTH(NSTEPPLUS+1+J1)=PATHLENGTH(NSTEPPLUS+J1)+SQRT(TEMP)
         IF (ZSYMSAVE(1:1).EQ.'W') THEN
C           DO J2=1,9*NATOMS  !  WCOMMENT
            DO J2=1,9*(NATOMS/2)
               Q1(J2)=Q2(J2)
            ENDDO
         ELSE
            DO J2=1,NOPT
               Q1(J2)=Q2(J2)
            ENDDO
         ENDIF
      ENDDO

555   CONTINUE
      ! Write the data we have just obtained to the EofS file for this path
      OPEN(UNIT=3,FILE=EOFSSTRING,STATUS='UNKNOWN')
      WRITE(3,'(2G25.15,I6)') (PATHLENGTH(J1),EOFS(J1),J1,J1=1,NSTEPPLUS+NSTEPMINUS+1)
      CLOSE(3)
      
      SUM2=0.0D0
      SUM4=0.0D0
C     NDUMMY=1 !  WCOMMENT
C     IF (ZSYMSAVE(1:1).EQ.'W') NDUMMY=9  !  WCOMMENT 9 should have been 3?
      NATOMSIMUL=NATOMS
      IF (ZSYMSAVE(1:1).EQ.'W') NATOMSIMUL=3*(NATOMS/2)
      IF (ZSYMSAVE(1:1).EQ.'W') THEN
         DO J1=1,NATOMSIMUL  !  WCOMMENT
            SUM2=SUM2+
     1           (QFRAMEP(3*(J1-1)+1,NFPLUS)-QFRAMEM(3*(J1-1)+1,NFMINUS))**2
     2          +(QFRAMEP(3*(J1-1)+2,NFPLUS)-QFRAMEM(3*(J1-1)+2,NFMINUS))**2
     3          +(QFRAMEP(3*(J1-1)+3,NFPLUS)-QFRAMEM(3*(J1-1)+3,NFMINUS))**2
            SUM4=SUM4+
     1          ((QFRAMEP(3*(J1-1)+1,NFPLUS)-QFRAMEM(3*(J1-1)+1,NFMINUS))**2
     2          +(QFRAMEP(3*(J1-1)+2,NFPLUS)-QFRAMEM(3*(J1-1)+2,NFMINUS))**2
     3          +(QFRAMEP(3*(J1-1)+3,NFPLUS)-QFRAMEM(3*(J1-1)+3,NFMINUS))**2)**2
         ENDDO
!        sn402: to test
      ELSE IF (RIGIDINIT .AND. .NOT. ATOMRIGIDCOORDT) THEN
          CALL RB_DISTANCE(SUM2,QPLUS(1:DEGFREEDOMS),Q(1:DEGFREEDOMS),DUMMYA,DUMMYB,.FALSE.)
          SUM2 = SUM2**2
          SUM4 = SUM2**2
      ELSE IF (BULKT.AND..NOT.VASP) THEN
         IF (TWOD) THEN
            DO J1=1,NATOMSIMUL  !  WCOMMENT
               SUM2=SUM2+
     1              MINIM(QPLUS(3*(J1-1)+1),Q(3*(J1-1)+1),PARAM1)**2
     2             +MINIM(QPLUS(3*(J1-1)+2),Q(3*(J1-1)+2),PARAM2)**2
               SUM4=SUM4+
     1             (MINIM(QPLUS(3*(J1-1)+1),Q(3*(J1-1)+1),PARAM1)**2
     2             +MINIM(QPLUS(3*(J1-1)+2),Q(3*(J1-1)+2),PARAM2)**2)**2
            ENDDO
         ELSE
C           DO J1=1,NATOMS*NDUMMY
            DO J1=1,NATOMSIMUL  !  WCOMMENT
               SUM2=SUM2+
     1              MINIM(QPLUS(3*(J1-1)+1),Q(3*(J1-1)+1),PARAM1)**2
     2             +MINIM(QPLUS(3*(J1-1)+2),Q(3*(J1-1)+2),PARAM2)**2
     3             +MINIM(QPLUS(3*(J1-1)+3),Q(3*(J1-1)+3),PARAM3)**2
               SUM4=SUM4+
     1             (MINIM(QPLUS(3*(J1-1)+1),Q(3*(J1-1)+1),PARAM1)**2
     2             +MINIM(QPLUS(3*(J1-1)+2),Q(3*(J1-1)+2),PARAM2)**2
     3             +MINIM(QPLUS(3*(J1-1)+3),Q(3*(J1-1)+3),PARAM3)**2)**2
            ENDDO
         ENDIF
      ELSE IF (ZSYMSAVE.EQ.'CD') THEN
         DO J2=1,NATOMS/2
            CALL CAPSIDIO(QPLUS(3*(J2-1)+1),QPLUS(3*(J2-1)+2),QPLUS(3*(J2-1)+3),
     1                    QPLUS(3*(NATOMS/2+J2-1)+1),QPLUS(3*(NATOMS/2+J2-1)+2),QPLUS(3*(NATOMS/2+J2-1)+3),
     2                    CAPSCOORDS2,RAD,HEIGHT)
            CALL CAPSIDIO(Q(3*(J2-1)+1),Q(3*(J2-1)+2),Q(3*(J2-1)+3),
     1                    Q(3*(NATOMS/2+J2-1)+1),Q(3*(NATOMS/2+J2-1)+2),Q(3*(NATOMS/2+J2-1)+3),
     2                    CAPSCOORDS1,RAD,HEIGHT)
! khs26> This previously used to go up to 18, which is out-of-bounds. There are 18 coordinates
! but only 6 atoms.
            DO J1=1,6
               SUM2=SUM2+
     1              (CAPSCOORDS1(3*(J1-1)+1)-CAPSCOORDS2(3*(J1-1)+1))**2
     2             +(CAPSCOORDS1(3*(J1-1)+2)-CAPSCOORDS2(3*(J1-1)+2))**2
     3             +(CAPSCOORDS1(3*(J1-1)+3)-CAPSCOORDS2(3*(J1-1)+3))**2
               SUM4=SUM4+
     1             ((CAPSCOORDS1(3*(J1-1)+1)-CAPSCOORDS2(3*(J1-1)+1))**2
     2             +(CAPSCOORDS1(3*(J1-1)+2)-CAPSCOORDS2(3*(J1-1)+2))**2
     3             +(CAPSCOORDS1(3*(J1-1)+3)-CAPSCOORDS2(3*(J1-1)+3))**2)**2
            ENDDO
         ENDDO
      ELSEIF (RINGPOLYMERT.OR.VARIABLES) THEN
         SUM2=1.0D0
         SUM4=1.0D0
      ELSE
C        DO J1=1,NATOMS*NDUMMY  !  WCOMMENT
         DO J1=1,NATOMSIMUL
            SUM2=SUM2+
     1           (QPLUS(3*(J1-1)+1)-Q(3*(J1-1)+1))**2
     2          +(QPLUS(3*(J1-1)+2)-Q(3*(J1-1)+2))**2
     3          +(QPLUS(3*(J1-1)+3)-Q(3*(J1-1)+3))**2
            SUM4=SUM4+
     1          ((QPLUS(3*(J1-1)+1)-Q(3*(J1-1)+1))**2
     2          +(QPLUS(3*(J1-1)+2)-Q(3*(J1-1)+2))**2
     3          +(QPLUS(3*(J1-1)+3)-Q(3*(J1-1)+3))**2)**2
         ENDDO
      ENDIF

      ETS=EOFS(NSTEPPLUS+1)
      EPLUS=EOFS(1)
      EMINUS=EOFS(NSTEPPLUS+NSTEPMINUS+1)
      PRINT*

      IF (RINGPOLYMERT) THEN
         WRITE(*,'(A)')
     1'         E+        Ets - E+           Ets       Ets - E-           E-'
         WRITE(*,60) EPLUS,ETS-EPLUS,ETS,ETS-EMINUS,EMINUS
      ELSE
         WRITE(*,'(A)')
     1'         E+        Ets - E+           Ets       Ets - E-           E-          S       D' //
     2 '      gamma   ~N'
         WRITE(*,60) EPLUS,ETS-EPLUS,ETS,ETS-EMINUS,EMINUS,
     1         PATHLENGTH(NSTEPPLUS+NSTEPMINUS+1)-PATHLENGTH(1),SQRT(SUM2),SUM4*NATOMSIMUL/SUM2**2,SUM2**2/SUM4
60       FORMAT(F17.7,G12.5,F17.7,G12.5,F17.7,F8.3,F8.3,F8.3,F8.3)
      ENDIF
C
C tvb Calculation of catastrophe ratios
      IF (RATIOS) THEN
            NEG=1
            IF (EOFS(1).GT.EOFS(NSTEPPLUS+NSTEPMINUS+1)) THEN ! swap sides of path
               NEG=-1
            ENDIF
            DO J1=1,NSTEPPLUS+NSTEPMINUS+1
               PATHLENGTH(J1)=NEG*PATHLENGTH(J1)
            ENDDO
            LUNIT=GETUNIT()
            OPEN(UNIT=LUNIT,FILE='EofS.fold',STATUS='UNKNOWN')
            DO J1=1,NSTEPPLUS+NSTEPMINUS+1
               WRITE(LUNIT,'(3G20.10)') PATHLENGTH(J1),EOFS(J1),EOFS(J1)-MAX(EOFS(1),EOFS(NSTEPPLUS+NSTEPMINUS+1))
            ENDDO
            CLOSE(LUNIT)

            ! sn402: Added IF blocks to deal with rigid body systems here
            IF(RIGIDINIT) THEN
                WRITE(*,*) "path> Warning: RIGIDINIT not tested with RATIOS"
                ! sn402: Not sure whether ATMASS is always set here, and if it isn't then this call won't work. Needs review.
                CALL GENRIGID_EIGENVALUES(QPLUS, ATMASS, DIAG, INFO)
                ! We need EPLUS (and possibly VNEW?) as well as the eigenvalues, so we need an extra potential call.
                CALL POTENTIAL(QPLUS,EPLUS,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            ELSEIF(RBAAT) THEN
                WRITE(*,*) "path> Warning: RATIOS not tested with rigid body potentials"
                CALL NRMLMD(QPLUS, DIAG, .FALSE.)
                ! We need EPLUS (and possibly VNEW?) as well as the eigenvalues, so we need an extra potential call.
                CALL POTENTIAL(QPLUS,EPLUS,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            ELSE
                CALL POTENTIAL(QPLUS,EPLUS,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
                CALL DSYEV('V','U',NOPT,HESS,SIZE(HESS,1),DIAG,TEMPA,3*NOPT,INFO)
            ENDIF
            CALL EIGENSORT_VAL_ASC(DIAG,HESS,NOPT,NOPT)

            IF (DEBUG) THEN
               PRINT '(A)','+ min energy'
               PRINT '(G20.10)',EPLUS
               PRINT '(A)','+ min eigenvalues:'
               PRINT '(3G20.10)',DIAG(1:NOPT)
            ENDIF

            LAMBDAP=DIAG(NOPT-NZERO)

            ! sn402: Added IF blocks to deal with rigid body systems here
            IF(RIGIDINIT) THEN
                WRITE(*,*) "path> Warning: RIGIDINIT not tested with RATIOS"
                ! sn402: Not sure whether ATMASS is always set here, and if it isn't then this call won't work. Needs review.
                CALL GENRIGID_EIGENVALUES(QMINUS, ATMASS, DIAG, INFO)
                CALL POTENTIAL(QMINUS,EMINUS,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            ELSEIF(RBAAT) THEN
                WRITE(*,*) "path> Warning: RATIOS not tested with rigid body potentials"
                CALL NRMLMD(QMINUS, DIAG, .FALSE.)
                CALL POTENTIAL(QMINUS,EMINUS,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            ELSE
                CALL POTENTIAL(QMINUS,EMINUS,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
                CALL DSYEV('N','U',NOPT,HESS,SIZE(HESS,1),DIAG,TEMPA,3*NOPT,INFO)
            ENDIF
            CALL EIGENSORT_VAL_ASC(DIAG,HESS,NOPT,NOPT)

            IF (DEBUG) THEN
               PRINT '(A)','- min energy'
               PRINT '(G20.10)',EMINUS
               PRINT '(A)','- min eigenvalues:'
               PRINT '(3G20.10)',DIAG(1:NOPT)
            ENDIF

            LAMBDAM=DIAG(NOPT-NZERO)

            ! sn402: Added IF blocks to deal with rigid body systems here
            IF(RIGIDINIT) THEN
                WRITE(*,*) "path> Warning: RIGIDINIT not tested with RATIOS"
                ! sn402: Not sure whether ATMASS is always set here, and if it isn't then this call won't work. Needs review.
                CALL GENRIGID_EIGENVALUES(QINIT, ATMASS, DIAG, INFO)
                CALL POTENTIAL(QINIT,EOFS(NSTEPPLUS+1),VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            ELSEIF(RBAAT) THEN
                WRITE(*,*) "path> Warning: RATIOS not tested with rigid body potentials"
                CALL NRMLMD(QINIT, DIAG, .FALSE.)
                CALL POTENTIAL(QINIT,EOFS(NSTEPPLUS+1),VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            ELSE
                CALL POTENTIAL(QINIT,EOFS(NSTEPPLUS+1),VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
                CALL DSYEV('N','U',NOPT,HESS,SIZE(HESS,1),DIAG,TEMPA,3*NOPT,INFO)
            ENDIF
            CALL EIGENSORT_VAL_ASC(DIAG,HESS,NOPT,NOPT)

            IF (DEBUG) THEN
               PRINT '(A)','ts energy'
               PRINT '(G20.10)',EOFS(NSTEPPLUS+1)
               PRINT '(A)','ts eigenvalues:'
               PRINT '(3G20.10)',DIAG(1:NOPT)
            ENDIF
            PRINT '(A,I8,G20.10)','for ts: NOPT,DIAG=',NOPT,DIAG(NOPT)
            LAMBDATS=DIAG(NOPT)
            CALL NEWMINDIST(QPLUS,QINIT,NATOMS,DISTP,.FALSE.,.FALSE.,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
            CALL NEWMINDIST(QMINUS,QINIT,NATOMS,DISTM,.FALSE.,.FALSE.,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
            LUNIT=GETUNIT()
            IF (MULTIJOBT) THEN
               OPEN(UNIT=LUNIT,FILE='folddata',STATUS='UNKNOWN',POSITION='APPEND')
            ELSE
               OPEN(UNIT=LUNIT,FILE='folddata',STATUS='UNKNOWN')
            ENDIF
C           Folddata: Em1, Em2, Ets, Sm1, Sm2, Evts, Ev1, Ev2, Sm1min, Sm2min
            WRITE(LUNIT,'(10F13.7)') EOFS(1), EOFS(NSTEPPLUS+NSTEPMINUS+1), EOFS(NSTEPPLUS+1), PATHLENGTH(1),
     &                               PATHLENGTH(NSTEPPLUS+NSTEPMINUS+1), LAMBDATS, LAMBDAP, LAMBDAM, DISTP, DISTM
            CLOSE(LUNIT)
            LUNIT=GETUNIT()
            IF (MULTIJOBT) THEN
               OPEN(UNIT=LUNIT,FILE='foldratios',STATUS='UNKNOWN',POSITION='APPEND')
            ELSE
               OPEN(UNIT=LUNIT,FILE='foldratios',STATUS='UNKNOWN')
            ENDIF
            IF (NEG.EQ.1) THEN
               WRITE(LUNIT,'(12G15.5)') 
     &            6*(Eofs(nstepplus+1)-Eofs(1))/(ABS(lambdats)*PATHLENGTH(1)**2),
     &           -6*(Eofs(nstepplus+nstepminus+1)-Eofs(nstepplus+1))/(ABS(lambdats)*PATHLENGTH(NSTEPPLUS+NSTEPMINUS+1)**2),
     &            6*(Eofs(nstepplus+1)-Eofs(1))/(ABS(lambdats)*DISTP**2),
     &           -6*(Eofs(nstepplus+nstepminus+1)-Eofs(nstepplus+1))/(ABS(lambdats)*DISTM**2),
     &            6*(Eofs(nstepplus+1)-Eofs(1))/(lambdap*PATHLENGTH(1)**2),
     &           -6*(Eofs(nstepplus+nstepminus+1)-Eofs(nstepplus+1))/(lambdam*PATHLENGTH(NSTEPPLUS+NSTEPMINUS+1)**2),
     &            6*(Eofs(nstepplus+1)-Eofs(1))/(lambdap*DISTP**2),
     &           -6*(Eofs(nstepplus+nstepminus+1)-Eofs(nstepplus+1))/(lambdam*DISTM**2),
     &            PATHLENGTH(1),PATHLENGTH(NSTEPPLUS+NSTEPMINUS+1),DISTP, DISTM
            ELSE
               WRITE(LUNIT,'(12G15.5)') 
     &           -6*(Eofs(nstepplus+nstepminus+1)-Eofs(nstepplus+1))/(ABS(lambdats)*PATHLENGTH(NSTEPPLUS+NSTEPMINUS+1)**2),
     &            6*(Eofs(nstepplus+1)-Eofs(1))/(ABS(lambdats)*PATHLENGTH(1)**2),
     &           -6*(Eofs(nstepplus+nstepminus+1)-Eofs(nstepplus+1))/(ABS(lambdats)*DISTM**2),
     &            6*(Eofs(nstepplus+1)-Eofs(1))/(ABS(lambdats)*DISTP**2),
     &           -6*(Eofs(nstepplus+nstepminus+1)-Eofs(nstepplus+1))/(lambdam*PATHLENGTH(NSTEPPLUS+NSTEPMINUS+1)**2),
     &            6*(Eofs(nstepplus+1)-Eofs(1))/(lambdap*PATHLENGTH(1)**2),
     &           -6*(Eofs(nstepplus+nstepminus+1)-Eofs(nstepplus+1))/(lambdam*DISTM**2),
     &            6*(Eofs(nstepplus+1)-Eofs(1))/(lambdap*DISTP**2),
     &            PATHLENGTH(NSTEPPLUS+NSTEPMINUS+1),PATHLENGTH(1),DISTM, DISTP
            ENDIF
            CLOSE(LUNIT)
      ENDIF
C end tvb
C (i.e. end of IF(RATIOS)

      SLENGTH=PATHLENGTH(NSTEPPLUS+NSTEPMINUS+1)-PATHLENGTH(1)
      DISP=SQRT(SUM2)
C     GAMMA=SUM4*NDUMMY*NATOMS/SUM2**2  !  WCOMMENT
      GAMMA=SUM4*NATOMSIMUL/SUM2**2
      NTILDE=SUM2**2/SUM4
      IF (CHECKINDEX.AND.RATIOS.AND..NOT.CONNECTT) THEN
         SMINUS=PATHLENGTH(NSTEPPLUS+NSTEPMINUS+1)
         SPLUS=PATHLENGTH(1)
         STS=PATHLENGTH(NSTEPPLUS+1)
         IF (EPLUS.LT.EMINUS) THEN
            ETEMP=EPLUS
            STEMP=SPLUS
            DUMMY=EVPLUS
            EPLUS=EMINUS
            SPLUS=SMINUS
            EVPLUS=EVMINUS
            EMINUS=ETEMP
            SMINUS=STEMP
            EVMINUS=DUMMY
         ENDIF
         PRINT*
         WRITE(*,'(A)')
     1      '         evts         evplus/ts    evminus/ts  del e plus      del e minus      s plus         s minus          frat +'
     2        // '    frat -  '

         WRITE(*,'(A5,3F13.7,4G16.8,2F10.4)') 'fold ',EVTS,EVPLUS/ABS(EVTS),EVMINUS/ABS(EVTS),ETS-EPLUS,ETS-EMINUS,ABS(STS-SPLUS),
     1              ABS(SMINUS-STS),6*(ETS-EPLUS)/(ABS(EVTS)*(STS-SPLUS)**2),6*(ETS-EMINUS)/(ABS(EVTS)*(STS-SMINUS)**2)

         WRITE(*,'(A)') '        del e plus      del e minus     s plus           s minus       '
     1        // '   evrat     dev1      dev2      dev3      rat+      rat-'

         WRITE(*,'(A5,4G16.8,6F10.4)') 'cusp ',ETS-EPLUS,ETS-EMINUS,ABS(STS-SPLUS),ABS(SMINUS-STS),
     1                           EVPLUS*EVMINUS/((EVPLUS+EVMINUS)*EVTS),
     2                           (ETS-EPLUS)*EVMINUS**3*(EVMINUS+2.0D0*EVPLUS)/((ETS-EMINUS)*EVPLUS**3*(EVPLUS+2.0D0*EVMINUS)),
     3                           (ETS-EPLUS)*EVTS**3*(EVTS+2.0D0*EVPLUS)/((ETS-EMINUS)*(EVTS-EVPLUS)*(EVPLUS+EVTS)**3),
     4                           (ETS-EPLUS)*(EVTS-EVMINUS)*(EVMINUS+EVTS)**3/((ETS-EMINUS)*EVTS**3*(2.0D0*EVMINUS+EVTS)),
     5                           12*(ETS-EPLUS) /((EVPLUS -EVTS)*(STS-SPLUS )**2),
     6                           12*(ETS-EMINUS)/((EVMINUS-EVTS)*(STS-SMINUS)**2)
!       sf344> now swap back EPLUS and EMINUS to their original value. Just to be sure that the plus and minus minima are never
!              saved with swapped energies in the path.info file!
         IF (EPLUS.LT.EMINUS) THEN
            ETEMP=EPLUS
            STEMP=SPLUS
            DUMMY=EVPLUS
            EPLUS=EMINUS
            SPLUS=SMINUS
            EVPLUS=EVMINUS
            EMINUS=ETEMP
            SMINUS=STEMP
            EVMINUS=DUMMY
         ENDIF
      ENDIF

      ! Open .xyz file for this path.
      OPEN(UNIT=3,FILE=ITSTRING,STATUS='UNKNOWN')

      ! ZSYMSAVE == W CD; ELSE  
C     IF (ZSYMSAVE(1:1).EQ.'W') THEN
C        WRITE(3,'(I6)') 3*NATOMS
C        WRITE(3,'(g20.10)') EOFS(1)
C        DO J1=1,NATOMS
C           CALL CONVERT(QPLUS(3*(J1-1)+1),QPLUS(3*(J1-1)+2),QPLUS(3*(J1-1)+3),
C    1                   QPLUS(3*(NATOMS+J1-1)+1),QPLUS(3*(NATOMS+J1-1)+2),QPLUS(3*(NATOMS+J1-1)+3),OVEC,H1VEC,H2VEC)
C           WRITE(3,'(A2,4X,3g20.10)') 'O  ',OVEC(1),OVEC(2),OVEC(3)
C           WRITE(3,'(A2,4X,3g20.10)') 'H  ',H1VEC(1),H1VEC(2),H1VEC(3)
C           WRITE(3,'(A2,4X,3g20.10)') 'H  ',H2VEC(1),H2VEC(2),H2VEC(3)
C        ENDDO
C     ELSE IF (ZSYMSAVE.EQ.'CD') THEN
C        WRITE(3,'(I6)') NATOMS*6/2
C        WRITE(3,'(g20.10)') EOFS(1)
C        DO J2=1,NATOMS/2
C           CALL CAPSIDIO(QPLUS(3*(J2-1)+1),QPLUS(3*(J2-1)+2),QPLUS(3*(J2-1)+3),
C    1                    QPLUS(3*(NATOMS/2+J2-1)+1),QPLUS(3*(NATOMS/2+J2-1)+2),QPLUS(3*(NATOMS/2+J2-1)+3),CAPSCOORDS2,RAD,HEIGHT)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(1),CAPSCOORDS2(2),CAPSCOORDS2(3)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(4),CAPSCOORDS2(5),CAPSCOORDS2(6)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(7),CAPSCOORDS2(8),CAPSCOORDS2(9)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(10),CAPSCOORDS2(11),CAPSCOORDS2(12)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(13),CAPSCOORDS2(14),CAPSCOORDS2(15)
C           WRITE(3,'(A2,4X,3g20.10)') 'C4  ',CAPSCOORDS2(16),CAPSCOORDS2(17),CAPSCOORDS2(18)
C        ENDDO
C     ELSE
C        WRITE(3,'(I6)') NATOMS
C        WRITE(3,'(g20.10)') EOFS(1)
C        WRITE(3,'(A2,4X,3g20.10)') (ZSYM(J1),QPLUS(3*(J1-1)+1),QPLUS(3*(J1-1)+2),QPLUS(3*(J1-1)+3),J1=1,NATOMS)
C     ENDIF

      ! Write + frame coordinates for this path
      DO J1=NFPLUS,1,-1
         IF (ZSYMSAVE(1:1).EQ.'W') THEN
C           WRITE(3,'(I6)') 3*NATOMS  !  WCOMMENT
            WRITE(3,'(I6)') 3*(NATOMS/2)
            WRITE(3,'(A,g25.15)') 'Energy=',EOFSFRAMEP(J1)
C           IF (J1.EQ.NFPLUS) THEN 
C              WRITE(3,'(g20.10)') EOFS(1)
C           ELSE
C              WRITE(3,'(A)') ' '
C           ENDIF
            WRITE(3,'(A2,4X,3g20.10)') 
     1   ('O  ',QFRAMEP(9*(J2-1)+1,J1),QFRAMEP(9*(J2-1)+2,J1),QFRAMEP(9*(J2-1)+3,J1),
     1    'H  ',QFRAMEP(9*(J2-1)+4,J1),QFRAMEP(9*(J2-1)+5,J1),QFRAMEP(9*(J2-1)+6,J1),
C    2    'H  ',QFRAMEP(9*(J2-1)+7,J1),QFRAMEP(9*(J2-1)+8,J1),QFRAMEP(9*(J2-1)+9,J1),J2=1,NATOMS)
     2    'H  ',QFRAMEP(9*(J2-1)+7,J1),QFRAMEP(9*(J2-1)+8,J1),QFRAMEP(9*(J2-1)+9,J1),J2=1,NATOMS/2)
         ELSE IF (ZSYMSAVE.EQ.'CD') THEN
            WRITE(3,'(I6)') NATOMS*6/2
            WRITE(3,'(A,g25.15)') 'Energy=',EOFSFRAMEP(J1)
C           IF (J1.EQ.NFPLUS) THEN
C              WRITE(3,'(g20.10)') EOFS(1)
C           ELSE
C              WRITE(3,'(A)') ' '
C           ENDIF
            WRITE(3,'(A)') ' '
            DO J2=1,NATOMS/2
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEP(18*(J2-1)+1,J1),QFRAMEP(18*(J2-1)+2,J1),QFRAMEP(18*(J2-1)+3,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEP(18*(J2-1)+4,J1),QFRAMEP(18*(J2-1)+5,J1),QFRAMEP(18*(J2-1)+6,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEP(18*(J2-1)+7,J1),QFRAMEP(18*(J2-1)+8,J1),QFRAMEP(18*(J2-1)+9,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEP(18*(J2-1)+10,J1),QFRAMEP(18*(J2-1)+11,J1),QFRAMEP(18*(J2-1)+12,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEP(18*(J2-1)+13,J1),QFRAMEP(18*(J2-1)+14,J1),QFRAMEP(18*(J2-1)+15,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C4  ',QFRAMEP(18*(J2-1)+16,J1),QFRAMEP(18*(J2-1)+17,J1),QFRAMEP(18*(J2-1)+18,J1)
            ENDDO
         ELSEIF (STOCKT) THEN
            WRITE(3,'(I6)') (NATOMS/2)
            WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEP(J1)
            DO J2=1,(NATOMS/2)
               WRITE(3,'(A2,4X,3G20.10,A13,3G20.10)')
     &         ZSYM(J2),QFRAMEP(3*(J2-1)+1,J1),QFRAMEP(3*(J2-1)+2,J1),QFRAMEP(3*(J2-1)+3,J1),
     &                        ' atom_vector ',
     &                        SIN(QFRAMEP(3*((NATOMS/2)+J2-1)+1,J1))*COS(QFRAMEP(3*((NATOMS/2)+J2-1)+2,J1)),
     &                        SIN(QFRAMEP(3*((NATOMS/2)+J2-1)+1,J1))*SIN(QFRAMEP(3*((NATOMS/2)+J2-1)+2,J1)),
     &                        COS(QFRAMEP(3*((NATOMS/2)+J2-1)+1,J1))
            ENDDO

!         ELSEIF (STOCKAAT) THEN
!            WRITE(3,'(I6)') (NATOMS/2)
!            WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEP(J1)
!            DO J2=1,(NATOMS/2)
               
!               P(:) = QFRAMEP(3*((NATOMS/2)+J2-1)+1:3*((NATOMS/2)+J2-1)+3,J1)
!               CALL ROTMAT(P(:), RMAT(:,:))
!               P(:) = RMAT(:,3)     ! THE DIPOLE-VECTOR HAS TO BE ALONG THE Z-AXIS IN THE BODY-FRAME

!               WRITE(3,'(A1,4X,3G20.10,A13,3G20.10)')
!     &         'O', QFRAMEP(3*(J2-1)+1,J1), QFRAMEP(3*(J2-1)+2,J1), QFRAMEP(3*(J2-1)+3,J1),
!     &         ' atom_vector ', P(1), P(2), P(3)
!            ENDDO
!         ELSEIF (RBAAT) THEN
!            Q(:) = QFRAMEP(:,J1)
!            QE   = EOFSFRAMEP(J1)
!            CALL RBPATHFRAME(Q,QE)
!            WRITE(3,'(I6)') NATOMS
!            WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEP(J1)
!            DO J2=1,(NATOMS/2)
!               WRITE(3,'(A1,4X,6G20.10)')
!     &         'O', QFRAMEP(3*(J2-1)+1,J1), QFRAMEP(3*(J2-1)+2,J1), QFRAMEP(3*(J2-1)+3,J1),
!     &         QFRAMEP(3*((NATOMS/2)+J2-1)+1,J1), QFRAMEP(3*((NATOMS/2)+J2-1)+2,J1),
!     &         QFRAMEP(3*((NATOMS/2)+J2-1)+3,J1) 
!            ENDDO
!            DO J2=1,NATOMS
!               WRITE(3,'(A1,4X,3G20.10)')
!     &         'O', QFRAMEP(3*(J2-1)+1,J1), QFRAMEP(3*(J2-1)+2,J1), QFRAMEP(3*(J2-1)+3,J1)
!            ENDDO
         ELSEIF (AMHT) THEN
            GLY_COUNT = 0
            DO J2=1,NMRES
               IF (SEQ(J2).EQ.8) GLY_COUNT = GLY_COUNT +1
            ENDDO
            WRITE(3,'(I6)') NATOMS+GLY_COUNT
            WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEP(J1)

            GLY_COUNT = 0
            DO J2=1,NMRES
              IF (SEQ(J2).EQ.8) THEN
          WRITE(3,'(a5,1x,3f20.10)') 'C1   ',QFRAMEP(9*(J2-1)+1-GLY_COUNT*3,J1),QFRAMEP(9*(J2-1)+2-GLY_COUNT*3,J1),
     &                                  QFRAMEP(9*(J2-1)+3-GLY_COUNT*3,J1)
          WRITE(3,'(a5,1x,3f20.10)') 'C1   ',QFRAMEP(9*(J2-1)+1-GLY_COUNT*3,J1),QFRAMEP(9*(J2-1)+2-GLY_COUNT*3,J1),
     &                                  QFRAMEP(9*(J2-1)+3-GLY_COUNT*3,J1)
          WRITE(3,'(a5,1x,3f20.10)') 'O    ',QFRAMEP(9*(J2-1)+4-GLY_COUNT*3,J1),QFRAMEP(9*(J2-1)+5-GLY_COUNT*3,J1),
     &                                  QFRAMEP(9*(J2-1)+6-GLY_COUNT*3,J1)
                GLY_COUNT = GLY_COUNT +1
              ELSE
          WRITE(3,'(a5,1x,3f20.10)') 'C1   ',QFRAMEP(9*(J2-1)+1-GLY_COUNT*3,J1),QFRAMEP(9*(J2-1)+2-GLY_COUNT*3,J1),
     &                                  QFRAMEP(9*(J2-1)+3-GLY_COUNT*3,J1)
          WRITE(3,'(a5,1x,3f20.10)') 'C2   ',QFRAMEP(9*(J2-1)+4-GLY_COUNT*3,J1),QFRAMEP(9*(J2-1)+5-GLY_COUNT*3,J1),
     &                                  QFRAMEP(9*(J2-1)+6-GLY_COUNT*3,J1)
          WRITE(3,'(a5,1x,3f20.10)') 'O    ',QFRAMEP(9*(J2-1)+7-GLY_COUNT*3,J1),QFRAMEP(9*(J2-1)+8-GLY_COUNT*3,J1),
     &                                  QFRAMEP(9*(J2-1)+9-GLY_COUNT*3,J1)
              ENDIF
          ENDDO
         ELSEIF (RINGPOLYMERT.AND.(RPSYSTEM(1:4).EQ.'AECK')) THEN
            WRITE(3,'(I6)') NOPT
            WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEP(J1)
            WRITE(3,'(A2,4X,3G20.10)') (ZSYM(J2),QFRAMEP(J2,J1),0.0D0,0.0D0,J2=1,NOPT)
! hk286
         ELSEIF (RIGIDINIT) THEN
            WRITE(3,'(I6)') NATOMS
            WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEP(J1)
            CALL TRANSFORMRIGIDTOC (1, NRIGIDBODY, XCOORDS, QFRAMEP(1:DEGFREEDOMS,J1))
            WRITE(3,'(A2,4X,3G20.10)') 
     &           (ZSYM(J2),XCOORDS(3*(J2-1)+1),XCOORDS(3*(J2-1)+2),XCOORDS(3*(J2-1)+3),J2=1,NATOMS)

         ELSEIF (GTHOMSONT) THEN
            WRITE(3,'(I6)') NATOMS
            WRITE(3,'(A,G25.15)') '1Energy=',EOFSFRAMEP(J1)
            CALL GTHOMSONANGTOC(XCOORDSA,QFRAMEP(1:3*NATOMS,J1),NATOMS)
            WRITE(3,'(A2,4X,3G20.10)') 
     &           (ZSYM(J2),XCOORDSA(3*(J2-1)+1),XCOORDSA(3*(J2-1)+2),XCOORDSA(3*(J2-1)+3),J2=1,NATOMS)

         ELSE
            WRITE(3,'(I6)') NATOMS
            WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEP(J1)
C           IF (J1.EQ.NFPLUS) THEN
C              WRITE(3,'(g20.10)') EOFS(1)
C           ELSE
C              WRITE(3,'(A)') ' '
C           ENDIF
            if (ZSYM(NATOMS).eq.'SV') then
               do j2=1, NATOMS
                  if (MOD(j2,3).eq.0) then
                     WRITE(3,'(A,4X,3g20.10)') 'O',QFRAMEP(3*(J2-1)+1,J1),QFRAMEP(3*(J2-1)+2,J1),QFRAMEP(3*(J2-1)+3,J1)
                  else
                     WRITE(3,'(A,4X,3g20.10)') 'H',QFRAMEP(3*(J2-1)+1,J1),QFRAMEP(3*(J2-1)+2,J1),QFRAMEP(3*(J2-1)+3,J1)
                  endif
               enddo
            ELSEIF (VARIABLES) THEN
               WRITE(3,'(G20.10)') QFRAMEP(1:NOPT,J1)
            ELSE
               WRITE(3,'(A2,4X,3G20.10)') 
     &         (ZSYM(J2),QFRAMEP(3*(J2-1)+1,J1),QFRAMEP(3*(J2-1)+2,J1),QFRAMEP(3*(J2-1)+3,J1),J2=1,NATOMS)
            ENDIF
         ENDIF
      ENDDO


      ! Now write coords of the transition state to file.
!
!  ZSYMSAVE=W,CD  
!  STOCKT RBAAT AMHT RINGPOLYMERT ELSE
      IF (ZSYMSAVE(1:1).EQ.'W') THEN
C        WRITE(3,'(I6)') 3*NATOMS ! WCOMMENT
         WRITE(3,'(I6)') 3*(NATOMS/2)
         WRITE(3,'(A,g25.15)') 'Energy=',EOFS(NSTEPPLUS+1)
C        DO J1=1,NATOMS ! WCOMMENT
         DO J1=1,NATOMS/2
            CALL CONVERT(QINIT(3*(J1-1)+1),QINIT(3*(J1-1)+2),QINIT(3*(J1-1)+3),
C    1                   QINIT(3*(NATOMS+J1-1)+1),QINIT(3*(NATOMS+J1-1)+2),QINIT(3*(NATOMS+J1-1)+3),
     1                   QINIT(3*(NATOMS/2+J1-1)+1),QINIT(3*(NATOMS/2+J1-1)+2),QINIT(3*(NATOMS/2+J1-1)+3),
     2                   OVEC,H1VEC,H2VEC)
            WRITE(3,'(A2,4X,3g20.10)') 'O  ',OVEC(1),OVEC(2),OVEC(3)
            WRITE(3,'(A2,4X,3g20.10)') 'H  ',H1VEC(1),H1VEC(2),H1VEC(3)
            WRITE(3,'(A2,4X,3g20.10)') 'H  ',H2VEC(1),H2VEC(2),H2VEC(3)
         ENDDO
      ELSE IF (ZSYMSAVE.EQ.'CD') THEN
         WRITE(3,'(I6)') NATOMS*6/2
         WRITE(3,'(A,g25.15)') 'Energy=',EOFS(NSTEPPLUS+1)
         DO J2=1,NATOMS/2
            CALL CAPSIDIO(QINIT(3*(J2-1)+1),QINIT(3*(J2-1)+2),QINIT(3*(J2-1)+3),
     1                    QINIT(3*(NATOMS/2+J2-1)+1),QINIT(3*(NATOMS/2+J2-1)+2),QINIT(3*(NATOMS/2+J2-1)+3),
     2                    CAPSCOORDS2,RAD,HEIGHT)
            WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(1),CAPSCOORDS2(2),CAPSCOORDS2(3)
            WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(4),CAPSCOORDS2(5),CAPSCOORDS2(6)
            WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(7),CAPSCOORDS2(8),CAPSCOORDS2(9)
            WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(10),CAPSCOORDS2(11),CAPSCOORDS2(12)
            WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(13),CAPSCOORDS2(14),CAPSCOORDS2(15)
            WRITE(3,'(A2,4X,3g20.10)') 'C4  ',CAPSCOORDS2(16),CAPSCOORDS2(17),CAPSCOORDS2(18)
         ENDDO
      ELSEIF (STOCKT) THEN
         WRITE(3,'(I6)') (NATOMS/2)
         WRITE(3,'(A,g25.15)') 'Energy=',EOFS(NSTEPPLUS+1)
         DO J2=1,(NATOMS/2)
            WRITE(3,'(A2,4X,3G20.10,A13,3G20.10)')
     &         ZSYM(J2),QINIT(3*(J2-1)+1),QINIT(3*(J2-1)+2),QINIT(3*(J2-1)+3),
     &                                  ' atom_vector ',
     &                                  SIN(QINIT(3*((NATOMS/2)+J2-1)+1))*COS(QINIT(3*((NATOMS/2)+J2-1)+2)),
     &                                  SIN(QINIT(3*((NATOMS/2)+J2-1)+1))*SIN(QINIT(3*((NATOMS/2)+J2-1)+2)),
     &                                  COS(QINIT(3*((NATOMS/2)+J2-1)+1))
         ENDDO
!      ELSEIF (RBAAT) THEN
!         Q(:) = QINIT(:)
!         QE   = EOFSFRAMEP(NSTEPPLUS+1)
!         CALL RBPATHFRAME(Q,QE)
!         WRITE(3,'(I6)') NATOMS
!         WRITE(3,'(A,G25.15)') 'Energy=',EOFS(NSTEPPLUS+1)
!         DO J2=1,(NATOMS/2)
!            WRITE(3,'(A1,4X,6G20.10)')
!     &      'O', QINIT(3*(J2-1)+1), QINIT(3*(J2-1)+2), QINIT(3*(J2-1)+3),
!     &       QINIT(3*((NATOMS/2)+J2-1)+1), QINIT(3*((NATOMS/2)+J2-1)+2),
!     &       QINIT(3*((NATOMS/2)+J2-1)+3)
!         ENDDO
!         DO J2=1,NATOMS
!            WRITE(3,'(A1,4X,3G20.10)')
!     &      'O', QINIT(3*(J2-1)+1), QINIT(3*(J2-1)+2), QINIT(3*(J2-1)+3)
!         ENDDO

      ELSEIF (AMHT) THEN
            GLY_COUNT = 0
            DO J2=1,NMRES
               IF (SEQ(J2).EQ.8) GLY_COUNT = GLY_COUNT +1
            ENDDO
            WRITE(3,'(I6)') NATOMS+GLY_COUNT
            WRITE(3,'(A,G25.15)') 'Energy=',EOFS(NSTEPPLUS+1)
            GLY_COUNT = 0
            DO J2=1,NMRES
              IF (SEQ(J2).EQ.8) THEN
        WRITE(3,'(a5,1x,3f20.10)') 'C1   ',QINIT(9*(J2-1)+1-GLY_COUNT*3),QINIT(9*(J2-1)+2-GLY_COUNT*3),QINIT(9*(J2-1)+3-GLY_COUNT*3)
        WRITE(3,'(a5,1x,3f20.10)') 'C1   ',QINIT(9*(J2-1)+1-GLY_COUNT*3),QINIT(9*(J2-1)+2-GLY_COUNT*3),QINIT(9*(J2-1)+3-GLY_COUNT*3)
        WRITE(3,'(a5,1x,3f20.10)') 'O    ',QINIT(9*(J2-1)+4-GLY_COUNT*3),QINIT(9*(J2-1)+5-GLY_COUNT*3),QINIT(9*(J2-1)+6-GLY_COUNT*3)
                GLY_COUNT = GLY_COUNT +1
              ELSE
        WRITE(3,'(a5,1x,3f20.10)') 'C1   ',QINIT(9*(J2-1)+1-GLY_COUNT*3),QINIT(9*(J2-1)+2-GLY_COUNT*3),QINIT(9*(J2-1)+3-GLY_COUNT*3)
        WRITE(3,'(a5,1x,3f20.10)') 'C2   ',QINIT(9*(J2-1)+4-GLY_COUNT*3),QINIT(9*(J2-1)+5-GLY_COUNT*3),QINIT(9*(J2-1)+6-GLY_COUNT*3)
        WRITE(3,'(a5,1x,3f20.10)') 'O    ',QINIT(9*(J2-1)+7-GLY_COUNT*3),QINIT(9*(J2-1)+8-GLY_COUNT*3),QINIT(9*(J2-1)+9-GLY_COUNT*3)
              ENDIF
          ENDDO

      ELSEIF (RINGPOLYMERT.AND.(RPSYSTEM(1:4).EQ.'AECK')) THEN
         WRITE(3,'(I6)') NOPT
         WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEP(NSTEPPLUS+1)
         WRITE(3,'(A2,4X,3G20.10)') (ZSYM(J2),QINIT(J2),0.0D0,0.0D0,J2=1,NOPT)

! hk286
      ELSEIF (RIGIDINIT) THEN
         WRITE(3,'(I6)') NATOMS
         WRITE(3,'(A,g25.15)') 'Energy=',EOFS(NSTEPPLUS+1)
         CALL TRANSFORMRIGIDTOC (1, NRIGIDBODY, XCOORDS, QINIT(1:DEGFREEDOMS))
         WRITE(3,'(A2,4X,3G20.10)') 
     &        (ZSYM(J2),XCOORDS(3*(J2-1)+1),XCOORDS(3*(J2-1)+2),XCOORDS(3*(J2-1)+3),J2=1,NATOMS)

         ELSEIF (GTHOMSONT) THEN
            WRITE(3,'(I6)') NATOMS
            WRITE(3,'(A,G25.15)') '2Energy=',EOFS(NSTEPPLUS+1)  
            CALL GTHOMSONANGTOC(XCOORDSA,QINIT(1:3*NATOMS),NATOMS)
            WRITE(3,'(A2,4X,3G20.10)') 
     &           (ZSYM(J2),XCOORDSA(3*(J2-1)+1),XCOORDSA(3*(J2-1)+2),XCOORDSA(3*(J2-1)+3),J2=1,NATOMS)

      ELSE
         WRITE(3,'(I6)') NATOMS
         WRITE(3,'(A,g25.15)') 'Energy=',EOFS(NSTEPPLUS+1)

         if (ZSYM(NATOMS).eq.'SV') then
            do j1=1, NATOMS
               if (MOD(j1,3).eq.0) then
                  WRITE(3,'(A,4X,3g20.10)') 'O',QINIT(3*(J1-1)+1),QINIT(3*(J1-1)+2),QINIT(3*(J1-1)+3)
               else
                  WRITE(3,'(A,4X,3g20.10)') 'H',QINIT(3*(J1-1)+1),QINIT(3*(J1-1)+2),QINIT(3*(J1-1)+3)
               endif
            enddo
         ELSEIF (VARIABLES) THEN
            WRITE(3,'(G20.10)') QINIT(1:NOPT)
         ELSE
            WRITE(3,'(A2,4X,3G20.10)') (ZSYM(J1),QINIT(3*(J1-1)+1),QINIT(3*(J1-1)+2),QINIT(3*(J1-1)+3),J1=1,NATOMS)
         endif
      ENDIF

      ! Now - path frames
      DO J1=1,NFMINUS
         IF (ZSYMSAVE(1:1).EQ.'W') THEN
C           WRITE(3,'(I6)') 3*NATOMS  !  WCOMMENT
            WRITE(3,'(I6)') 3*(NATOMS/2)
            WRITE(3,'(A,g25.15)') 'Energy=',EOFSFRAMEM(J1)
C           IF (J1.EQ.NFMINUS) THEN
C              WRITE(3,'(g20.10)') EOFS(NSTEPPLUS+NSTEPMINUS+1)
C           ELSE
C              WRITE(3,'(A)') ' '
C           ENDIF
            WRITE(3,'(A2,4X,3g20.10)') 
     1            ('O  ',QFRAMEM(9*(J2-1)+1,J1),QFRAMEM(9*(J2-1)+2,J1),QFRAMEM(9*(J2-1)+3,J1),
     1             'H  ',QFRAMEM(9*(J2-1)+4,J1),QFRAMEM(9*(J2-1)+5,J1),QFRAMEM(9*(J2-1)+6,J1),
C    2             'H  ',QFRAMEM(9*(J2-1)+7,J1),QFRAMEM(9*(J2-1)+8,J1),QFRAMEM(9*(J2-1)+9,J1),J2=1,NATOMS)
     2             'H  ',QFRAMEM(9*(J2-1)+7,J1),QFRAMEM(9*(J2-1)+8,J1),QFRAMEM(9*(J2-1)+9,J1),J2=1,NATOMS/2)
         ELSE IF (ZSYMSAVE.EQ.'CD') THEN
            WRITE(3,'(I6)') NATOMS*6/2
            WRITE(3,'(A,g25.15)') 'Energy=',EOFSFRAMEM(J1)
C           IF (J1.EQ.NFMINUS) THEN
C              WRITE(3,'(g20.10)') EOFS(NSTEPPLUS+NSTEPMINUS+1)
C           ELSE
C              WRITE(3,'(A)') ' '
C           ENDIF
            DO J2=1,NATOMS/2
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEM(18*(J2-1)+1,J1),QFRAMEM(18*(J2-1)+2,J1),QFRAMEM(18*(J2-1)+3,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEM(18*(J2-1)+4,J1),QFRAMEM(18*(J2-1)+5,J1),QFRAMEM(18*(J2-1)+6,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEM(18*(J2-1)+7,J1),QFRAMEM(18*(J2-1)+8,J1),QFRAMEM(18*(J2-1)+9,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEM(18*(J2-1)+10,J1),QFRAMEM(18*(J2-1)+11,J1),QFRAMEM(18*(J2-1)+12,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C1  ',QFRAMEM(18*(J2-1)+13,J1),QFRAMEM(18*(J2-1)+14,J1),QFRAMEM(18*(J2-1)+15,J1)
               WRITE(3,'(A2,4X,3g20.10)') 'C4  ',QFRAMEM(18*(J2-1)+16,J1),QFRAMEM(18*(J2-1)+17,J1),QFRAMEM(18*(J2-1)+18,J1)
            ENDDO
         ELSEIF (STOCKT) THEN
            WRITE(3,'(I6)') (NATOMS/2)
            WRITE(3,'(A,g25.15)') 'Energy=',EOFSFRAMEM(J1)
            DO J2=1,(NATOMS/2)
               WRITE(3,'(A2,4X,3G20.10,A13,3G20.10)')
     &         ZSYM(J2),QFRAMEM(3*(J2-1)+1,J1),QFRAMEM(3*(J2-1)+2,J1),QFRAMEM(3*(J2-1)+3,J1),
     &                  ' atom_vector ',
     &                  SIN(QFRAMEM(3*((NATOMS/2)+J2-1)+1,J1))*COS(QFRAMEM(3*((NATOMS/2)+J2-1)+2,J1)),
     &                  SIN(QFRAMEM(3*((NATOMS/2)+J2-1)+1,J1))*SIN(QFRAMEM(3*((NATOMS/2)+J2-1)+2,J1)),
     &                  COS(QFRAMEM(3*((NATOMS/2)+J2-1)+1,J1))
            ENDDO
!         ELSEIF (RBAAT) THEN
!            Q(:) = QFRAMEM(:,J1)
!            QE   = EOFSFRAMEM(J1)
!            CALL RBPATHFRAME(Q,QE)
!            WRITE(3,'(I6)') NATOMS
!            WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEM(J1)
!            DO J2=1,(NATOMS/2)
!               WRITE(3,'(A1,4X,6G20.10)')
!     &         'O', QFRAMEM(3*(J2-1)+1,J1), QFRAMEM(3*(J2-1)+2,J1), QFRAMEM(3*(J2-1)+3,J1),
!     &         QFRAMEM(3*((NATOMS/2)+J2-1)+1,J1), QFRAMEM(3*((NATOMS/2)+J2-1)+2,J1),
!     &         QFRAMEM(3*((NATOMS/2)+J2-1)+3,J1)
!            ENDDO
!            DO J2=1,NATOMS
!               WRITE(3,'(A1,4X,3G20.10)')
!     &         'O', QFRAMEM(3*(J2-1)+1,J1), QFRAMEM(3*(J2-1)+2,J1), QFRAMEM(3*(J2-1)+3,J1)
!            ENDDO

         ELSEIF (AMHT) THEN
            GLY_COUNT = 0
            DO J2=1,NMRES
               IF (SEQ(J2).EQ.8) GLY_COUNT = GLY_COUNT +1
            ENDDO
            WRITE(3,'(I6)') NATOMS+GLY_COUNT
            WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEM(J1)

            GLY_COUNT = 0
            DO J2=1,NMRES
              IF (SEQ(J2).EQ.8) THEN
          WRITE(3,'(a5,1x,3f20.10)') 'C1   ',QFRAMEM(9*(J2-1)+1-GLY_COUNT*3,J1),QFRAMEM(9*(J2-1)+2-GLY_COUNT*3,J1),
     &                                  QFRAMEM(9*(J2-1)+3-GLY_COUNT*3,J1)
          WRITE(3,'(a5,1x,3f20.10)') 'C1   ',QFRAMEM(9*(J2-1)+1-GLY_COUNT*3,J1),QFRAMEM(9*(J2-1)+2-GLY_COUNT*3,J1),
     &                                  QFRAMEM(9*(J2-1)+3-GLY_COUNT*3,J1)
          WRITE(3,'(a5,1x,3f20.10)') 'O    ',QFRAMEM(9*(J2-1)+4-GLY_COUNT*3,J1),QFRAMEM(9*(J2-1)+5-GLY_COUNT*3,J1),
     &                                  QFRAMEM(9*(J2-1)+6-GLY_COUNT*3,J1)
                GLY_COUNT = GLY_COUNT +1
              ELSE
          WRITE(3,'(a5,1x,3f20.10)') 'C1   ',QFRAMEM(9*(J2-1)+1-GLY_COUNT*3,J1),QFRAMEM(9*(J2-1)+2-GLY_COUNT*3,J1),
     &                                  QFRAMEM(9*(J2-1)+3-GLY_COUNT*3,J1)
          WRITE(3,'(a5,1x,3f20.10)') 'C2   ',QFRAMEM(9*(J2-1)+4-GLY_COUNT*3,J1),QFRAMEM(9*(J2-1)+5-GLY_COUNT*3,J1),
     &                                  QFRAMEM(9*(J2-1)+6-GLY_COUNT*3,J1)
          WRITE(3,'(a5,1x,3f20.10)') 'O    ',QFRAMEM(9*(J2-1)+7-GLY_COUNT*3,J1),QFRAMEM(9*(J2-1)+8-GLY_COUNT*3,J1),
     &                                  QFRAMEM(9*(J2-1)+9-GLY_COUNT*3,J1)
              ENDIF
          ENDDO
         ELSEIF (RINGPOLYMERT.AND.(RPSYSTEM(1:4).EQ.'AECK')) THEN
            WRITE(3,'(I6)') NOPT
            WRITE(3,'(A,G25.15)') 'Energy=',EOFSFRAMEM(J1)
            WRITE(3,'(A2,4X,3G20.10)') (ZSYM(J2),QFRAMEM(J2,J1),0.0D0,0.0D0,J2=1,NOPT)
! hk286
         ELSEIF (RIGIDINIT) THEN
            WRITE(3,'(I6)') NATOMS
            WRITE(3,'(A,g25.15)') 'Energy=',EOFSFRAMEM(J1)
            CALL TRANSFORMRIGIDTOC (1, NRIGIDBODY, XCOORDS, QFRAMEM(1:DEGFREEDOMS,J1))
            WRITE(3,'(A2,4X,3g20.10)') 
     &           (ZSYM(J2),XCOORDS(3*(J2-1)+1),XCOORDS(3*(J2-1)+2),XCOORDS(3*(J2-1)+3),J2=1,NATOMS)

         ELSEIF (GTHOMSONT) THEN
            WRITE(3,'(I6)') NATOMS
            WRITE(3,'(A,G25.15)') '3Energy=',EOFSFRAMEM(J1)
            CALL GTHOMSONANGTOC(XCOORDSA,QFRAMEM(1:3*NATOMS,J1),NATOMS)
            WRITE(3,'(A2,4X,3G20.10)') 
     &           (ZSYM(J2),XCOORDSA(3*(J2-1)+1),XCOORDSA(3*(J2-1)+2),XCOORDSA(3*(J2-1)+3),J2=1,NATOMS)

         ELSE
            WRITE(3,'(I6)') NATOMS
            WRITE(3,'(A,g25.15)') 'Energy=',EOFSFRAMEM(J1)
C           IF (J1.EQ.NFMINUS) THEN
C              WRITE(3,'(g20.10)') EOFS(NSTEPPLUS+NSTEPMINUS+1)
C           ELSE
C              WRITE(3,'(A)') ' '
C           ENDIF

            if (ZSYM(NATOMS).eq.'SV') then
               do j2=1, NATOMS
                  if (MOD(j2,3).eq.0) then
                     WRITE(3,'(A,4X,3g20.10)') 'O',QFRAMEM(3*(J2-1)+1,J1),QFRAMEM(3*(J2-1)+2,J1),QFRAMEM(3*(J2-1)+3,J1)
                  else
                     WRITE(3,'(A,4X,3g20.10)') 'H',QFRAMEM(3*(J2-1)+1,J1),QFRAMEM(3*(J2-1)+2,J1),QFRAMEM(3*(J2-1)+3,J1)
                  endif
               enddo
            ELSEIF (VARIABLES) THEN
               WRITE(3,'(G20.10)') QFRAMEM(1:NOPT,J1)
            else
               WRITE(3,'(A2,4X,3g20.10)') 
     &         (ZSYM(J2),QFRAMEM(3*(J2-1)+1,J1),QFRAMEM(3*(J2-1)+2,J1),QFRAMEM(3*(J2-1)+3,J1),J2=1,NATOMS)
            endif
         ENDIF
      ENDDO
C     IF (ZSYMSAVE(1:1).EQ.'W') THEN
C        WRITE(3,'(I6)') 3*NATOMS
C        WRITE(3,'(g20.10)') EOFS(NSTEPPLUS+NSTEPMINUS+1)
C        DO J1=1,NATOMS
C           CALL CONVERT(Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3),
C    1                   Q(3*(NATOMS+J1-1)+1),Q(3*(NATOMS+J1-1)+2),Q(3*(NATOMS+J1-1)+3),OVEC,H1VEC,H2VEC)
C           WRITE(3,'(A2,4X,3g20.10)') 'O  ',OVEC(1),OVEC(2),OVEC(3)
C           WRITE(3,'(A2,4X,3g20.10)') 'H  ',H1VEC(1),H1VEC(2),H1VEC(3)
C           WRITE(3,'(A2,4X,3g20.10)') 'H  ',H2VEC(1),H2VEC(2),H2VEC(3)
C        ENDDO
C     ELSE IF (ZSYMSAVE.EQ.'CD') THEN
C        WRITE(3,'(I6)') NATOMS*6/2
C        WRITE(3,'(g20.10)') EOFS(NSTEPPLUS+NSTEPMINUS+1)
C        DO J2=1,NATOMS/2
C           CALL CAPSIDIO(Q(3*(J2-1)+1),Q(3*(J2-1)+2),Q(3*(J2-1)+3),
C    1                    Q(3*(NATOMS/2+J2-1)+1),Q(3*(NATOMS/2+J2-1)+2),Q(3*(NATOMS/2+J2-1)+3),CAPSCOORDS2,RAD,HEIGHT)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(1),CAPSCOORDS2(2),CAPSCOORDS2(3)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(4),CAPSCOORDS2(5),CAPSCOORDS2(6)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(7),CAPSCOORDS2(8),CAPSCOORDS2(9)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(10),CAPSCOORDS2(11),CAPSCOORDS2(12)
C           WRITE(3,'(A2,4X,3g20.10)') 'C1  ',CAPSCOORDS2(13),CAPSCOORDS2(14),CAPSCOORDS2(15)
C           WRITE(3,'(A2,4X,3g20.10)') 'C4  ',CAPSCOORDS2(16),CAPSCOORDS2(17),CAPSCOORDS2(18)
C        ENDDO
C     ELSE
C        WRITE(3,'(I6)') NATOMS
C        WRITE(3,'(g20.10)') EOFS(NSTEPPLUS+NSTEPMINUS+1)
C        WRITE(3,'(A2,4X,3g20.10)') (ZSYM(J1),Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3),J1=1,NATOMS)
C     ENDIF

      CLOSE(3)  ! Finished writing to xyz file.

C jmc for unres, to put in the dummy peptide groups
      IF (UNRST) THEN
         WRITE(ITSTRING2,'(A)') 'unr.'//TRIM(ADJUSTL(ITSTRING))
         OPEN(UNIT=3,FILE=ITSTRING2,STATUS='UNKNOWN')
         DO J1=NFPLUS,1,-1
            DO K1=1,(NATOMS/2)-1
               DO K2=1,3
                  PEPCOORDS(6*(K1-1)+K2)=(2.0D0*QFRAMEP(6*(K1-1)+K2,J1)+QFRAMEP(6*K1+K2,J1))/3.0D0
                  PEPCOORDS(6*(K1-1)+K2+3)=(QFRAMEP(6*(K1-1)+K2,J1)+2.0D0*QFRAMEP(6*K1+K2,J1))/3.0D0
               END DO
            END DO
            WRITE(3,'(I6)') 2*NATOMS-2
            WRITE(3,'(A,g25.15)') 'Energy=',EOFSFRAMEP(J1)
C           IF (J1.EQ.NFPLUS) THEN
C              WRITE(3,'(g20.10)') EOFS(1)
C           ELSE
C              WRITE(3,'(A)') ' '
C           ENDIF
            WRITE(3,'(A2,4X,3g20.10)') (ZSYM(J2),QFRAMEP(3*(J2-1)+1,J1),QFRAMEP(3*(J2-1)+2,J1),QFRAMEP(3*(J2-1)+3,J1),J2=1,NATOMS)
            WRITE(3,'(A2,4X,3g20.10)') ('O ',PEPCOORDS(6*(K1-1)+1),PEPCOORDS(6*(K1-1)+2),PEPCOORDS(6*(K1-1)+3),K1=1,(NATOMS/2)-1)
            WRITE(3,'(A2,4X,3g20.10)') ('N ',PEPCOORDS(6*(K1-1)+4),PEPCOORDS(6*(K1-1)+5),PEPCOORDS(6*(K1-1)+6),K1=1,(NATOMS/2)-1)
         ENDDO
         DO K1=1,(NATOMS/2)-1
            DO K2=1,3
               PEPCOORDS(6*(K1-1)+K2)=(2.0D0*QINIT(6*(K1-1)+K2)+QINIT(6*K1+K2))/3.0D0
               PEPCOORDS(6*(K1-1)+K2+3)=(QINIT(6*(K1-1)+K2)+2.0D0*QINIT(6*K1+K2))/3.0D0
            END DO
         END DO
         WRITE(3,'(I6)') 2*NATOMS-2
         WRITE(3,'(A,g25.15)') 'Energy=',EOFS(NSTEPPLUS+1)
         WRITE(3,'(A2,4X,3g20.10)') (ZSYM(J1),QINIT(3*(J1-1)+1),QINIT(3*(J1-1)+2),QINIT(3*(J1-1)+3),J1=1,NATOMS)
         WRITE(3,'(A2,4X,3g20.10)') ('O ',PEPCOORDS(6*(K1-1)+1),PEPCOORDS(6*(K1-1)+2),PEPCOORDS(6*(K1-1)+3),K1=1,(NATOMS/2)-1)
         WRITE(3,'(A2,4X,3g20.10)') ('N ',PEPCOORDS(6*(K1-1)+4),PEPCOORDS(6*(K1-1)+5),PEPCOORDS(6*(K1-1)+6),K1=1,(NATOMS/2)-1)
         DO J1=1,NFMINUS
            DO K1=1,(NATOMS/2)-1
               DO K2=1,3
                  PEPCOORDS(6*(K1-1)+K2)=(2.0D0*QFRAMEM(6*(K1-1)+K2,J1)+QFRAMEM(6*K1+K2,J1))/3.0D0
                  PEPCOORDS(6*(K1-1)+K2+3)=(QFRAMEM(6*(K1-1)+K2,J1)+2.0D0*QFRAMEM(6*K1+K2,J1))/3.0D0
               END DO
            END DO
            WRITE(3,'(I6)') 2*NATOMS-2
            WRITE(3,'(A,g25.15)') 'Energy=',EOFSFRAMEM(J1)
C           IF (J1.EQ.NFMINUS) THEN
C              WRITE(3,'(g20.10)') EOFS(NSTEPPLUS+NSTEPMINUS+1)
C           ELSE
C              WRITE(3,'(A)') ' '
C           ENDIF
            WRITE(3,'(A2,4X,3g20.10)') (ZSYM(J2),QFRAMEM(3*(J2-1)+1,J1),QFRAMEM(3*(J2-1)+2,J1),QFRAMEM(3*(J2-1)+3,J1),J2=1,NATOMS)
            WRITE(3,'(A2,4X,3g20.10)') ('O ',PEPCOORDS(6*(K1-1)+1),PEPCOORDS(6*(K1-1)+2),PEPCOORDS(6*(K1-1)+3),K1=1,(NATOMS/2)-1)
            WRITE(3,'(A2,4X,3g20.10)') ('N ',PEPCOORDS(6*(K1-1)+4),PEPCOORDS(6*(K1-1)+5),PEPCOORDS(6*(K1-1)+6),K1=1,(NATOMS/2)-1)
         ENDDO
         CLOSE(3)
      END IF ! unrst
!
! This is where the path.info file is dumped for PATH runs without a
! CONNECT or NEWCONNECT keyword. It is a triple for both DUMPPATH and
! DUMPALLPATHS.
!
      IF ((DUMPPATH.OR.DUMPALLPATHS).AND.(.NOT.CONNECTT)) THEN
         IF (UNRST) WRITE(*,'(A)') '*** NOTE - pathlengths calculated from saved Cartesian coords will be rubbish
     & as they have been placed in the standard unres orientation.'
         IF (ZSYMSAVE.EQ.'CD') WRITE(*,'(A)') 'WARNING, symmetry and normal modes not implemented properly for CAPSID'


         CALL MAKE_PATHINFO_POINT(QPLUS, EPLUS, 'minhess.plus', .FALSE.)

         CALL MAKE_PATHINFO_POINT(QINIT, ETS, 'tshess', .TRUE.)

         CALL MAKE_PATHINFO_POINT(QMINUS, EMINUS, 'minhess.minus', .TRUE.)

         CLOSE(88)

      else if (machine.and..not.connectt) then
C SAT this is for the case when we need points for minima to be output in binary format, but do not want expensive Hessian
C diagonalization, which is required to produce "path.info" file
         inquire(iolength=reclen) (diag(J1),J1=1,3*Natoms)
         open(unit=38,file="points1.out",status='unknown',form='unformatted',access='direct',recl=reclen)
         write(38,rec=1) (QPLUS(J2),J2=1,NOPT)
         close(38)
         open(unit=38,file="points2.out",status='unknown',form='unformatted',access='direct',recl=reclen)
         write(38,rec=1) (QMINUS(J2),J2=1,NOPT)
         close(38)
      endif

      BFGSTST=BFGSTSTSAVE
      IVEC=IVECSAVE
      IVEC2=IVEC2SAVE

      IF (ALLOCATED(Q1)) DEALLOCATE(Q1)
      IF (ALLOCATED(Q2)) DEALLOCATE(Q2)
      IF (ALLOCATED(QW)) DEALLOCATE(QW)
      IF (ALLOCATED(QFRAMEP)) DEALLOCATE(QFRAMEP)
      IF (ALLOCATED(QFRAMEM)) DEALLOCATE(QFRAMEM)
      IF (ALLOCATED(EOFS)) DEALLOCATE(EOFS, PATHLENGTH, EOFSFRAMEP, EOFSFRAMEM)

      RETURN
      END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE MAKE_PATHINFO_POINT(Q, E, HESSDUMP_FNAME, ISTS)
      USE COMMONS
      USE KEY
      USE SYMINF
      USE modcharmm
      USE MODUNRES
      USE MODHESS
      USE GENRIGID

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)   :: Q(NOPT)          ! Coordinate vector for the point that we want to write out
      DOUBLE PRECISION, INTENT(IN)   :: E                ! Energy for the point that we want to write out
      CHARACTER(*), INTENT(IN) :: HESSDUMP_FNAME   ! Name of the file to which the hessian will be dumped
                                                         ! (if DUMPHESS is set)
      LOGICAL, INTENT(IN)            :: ISTS             ! TRUE if the point is a transition state, false otherwise

      ! Element labels for writing to path.info files
      CHARACTER(LEN=5) :: ZSYMSAVE
      COMMON /SYS/ ZSYMSAVE

      ! Variables for calls to POTENTIAL and SYMMETRY
      DOUBLE PRECISION :: VNEW(NOPT), DIAG(NOPT), RMS
      DOUBLE PRECISION :: INERTIA(3,3)
      INTEGER          :: HORDER

      ! Variables for dumping Hessian
      INTEGER :: LUNIT, GETUNIT

      ! Variables for old water potential
      INTEGER                                     :: NATOMSSAVE, IPOT
      DOUBLE PRECISION                            :: OVEC(3), H1VEC(3), H2VEC(3)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: QW

      ! Variables for UNRES, AMH
      INTEGER          :: KD, NNZ, NINTB, GLY_COUNT
      DOUBLE PRECISION :: DIHE, ALLANG

      ! Dummy variables
      DOUBLE PRECISION :: XCOORDS(NOPT), TMPCOORDS(NOPT), TEMPA(9*NATOMS)
      INTEGER          :: J1, J2, INFO

      ! We have identified a stationary point (minimum or ts) and now wish to write its information to a path.info file.
      ! We know the energy and coordinates of the point, we need to determine the point group and the normal mode frequencies.
      ! This is done here.
      ! (Note, sometimes we may have already calculated frequencies but if so they have probably been shifted, so we don't
      ! use them. We calculate new ones instead)


      ! First, we deal with a couple of the more special cases. Both of these should be largely obsolete now: the GENRIGID
      ! versions of rigid body potentials (including the TIP potentials) are more likely to be maintained.
      IF (ZSYMSAVE(1:1).EQ.'W') THEN
         IF (ZSYMSAVE.EQ.'W4') IPOT=4
         IF (ZSYMSAVE.EQ.'W3') IPOT=3
         IF (ZSYMSAVE.EQ.'W2') IPOT=2
         IF (ZSYMSAVE.EQ.'W1') IPOT=1
C        DO J2=1,NATOMS
         DO J2=1,NATOMS/2 ! WCOMMENT
            CALL CONVERT(Q(3*(J2-1)+1),Q(3*(J2-1)+2),Q(3*(J2-1)+3),
C    1                Q(3*(NATOMS+J2-1)+1),Q(3*(NATOMS+J2-1)+2),Q(3*(NATOMS+J2-1)+3),
     1                Q(3*(NATOMS/2+J2-1)+1),Q(3*(NATOMS/2+J2-1)+2),Q(3*(NATOMS/2+J2-1)+3),
     2                OVEC,H1VEC,H2VEC)
            QW(9*(J2-1)+1)=OVEC(1) ! WCOMMENT
            QW(9*(J2-1)+2)=OVEC(2)
            QW(9*(J2-1)+3)=OVEC(3)
            QW(9*(J2-1)+4)=H1VEC(1)
            QW(9*(J2-1)+5)=H1VEC(2)
            QW(9*(J2-1)+6)=H1VEC(3)
            QW(9*(J2-1)+7)=H2VEC(1)
            QW(9*(J2-1)+8)=H2VEC(2)
            QW(9*(J2-1)+9)=H2VEC(3)
         ENDDO
C        NATOMS=NATOMS*3 ! WCOMMENT
         NATOMSSAVE=NATOMS
         NATOMS=(NATOMS/2)*3
         CALL SYMMETRY(HORDER,.FALSE.,QW,INERTIA) ! WCOMMENT
C        NATOMS=NATOMS/3 ! WCOMMENT
         NATOMS=2*(NATOMS/3)
         NATOMS=NATOMSSAVE
         WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
C        CALL H2OMODES(NATOMS,IPOT,Q,DIAG) ! WCOMMENT
         CALL H2OMODES(NATOMS/2,IPOT,Q,DIAG)
C        WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,6*NATOMS) ! WCOMMENT
         IF (.NOT.NOFRQS) WRITE(88,'(3G20.10)') (DIAG(J2),J2=1,3*NATOMS)

      ! sn402: For the small number of potentials coded in the fully-rigid form (RBAAT), we can use the NRMLMD subroutine
      ! directly to calculate the Hessian and calculate its eigenvalues. We write squared angular frequencies in whatever
      ! time unit is specified by FRQCONV (see documentation and comments in keywords.f)
! hk286 - compute potential for normal modes, notice the toggle
      ELSE IF (RBAAT.AND..NOT.PYGPERIODICT.AND..NOT.SANDBOXT) THEN
         CALL POTENTIAL(Q,E,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
         WRITE(88, '(G20.10)') E
         ! The following line should definitely not be hard-coded in!
         WRITE(88, '(A,A)') "1  ", "C1" ! TEMP Solution

         IF (.NOT.NOFRQS) THEN
            RBAANORMALMODET = .TRUE.           
            CALL NRMLMD (Q, DIAG, .FALSE.)
            RBAANORMALMODET = .FALSE.
            WRITE(88,'(3G20.10)') (DIAG(J2)*FRQCONV2,J2=1,3*NATOMS)  ! sn402: Added the frequency conversion
         ENDIF
! hk286

      !!!!! This is the start of the main block to calculate normal mode frequencies. !!!!!
      !!!!! A similar procedure is followed for all remaining potentials. !!!!!

      ELSE

         IF (CHRMMT .OR. AMBERT .OR. AMBER12T .OR. NABT .OR. AMHT) THEN
            IF((AMBERT .OR. AMBER12T) .AND. MACROCYCLET) THEN
               CALL SYMMETRY(HORDER,.FALSE.,Q,INERTIA)
            ELSE
               HORDER=1
               FPGRP='C1'
            ENDIF
            ! We will shortly want to calculate the Hessian eigenvalues. For most systems, we start by computing the hessian
            ! using either MAKENUMHESS or POTENTIAL. With rigid bodies, the hessian is calculated within GENRIGID_EIGENVALUES,
            ! so we can save on a call to POTENTIAL if we do the complete diagonalisation here instead.
            IF (RIGIDINIT.AND.(.NOT.NOFRQS)) THEN
               ! This returns the square frequencies in internal units.
               CALL GENRIGID_EIGENVALUES(Q, ATMASS, DIAG, INFO)
               IF (DIAG(1).LT.DIAG(DEGFREEDOMS)) THEN
                  CALL EIGENSORT_VAL_ASC(DIAG(1:DEGFREEDOMS),HESS(1:DEGFREEDOMS,1:DEGFREEDOMS),DEGFREEDOMS,DEGFREEDOMS)
               ENDIF
            ELSE
               IF (ENDNUMHESS) THEN
                  CALL MAKENUMHESS(Q,NATOMS)
               ELSEIF (.NOT.NOFRQS .AND. (.NOT.(AMBERT .OR. AMBER12T))) THEN
                  ! sn402: I'm not sure why we don't do this for AMBER, but I'm retaining the original behaviour here.
                  CALL POTENTIAL(Q,E,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
               ENDIF
               IF (.NOT.NOFRQS) CALL MASSWT2(NATOMS,ATMASS,Q,VNEW,.TRUE.)
            ENDIF
         ELSE IF (UNRST) THEN
            IF(RIGIDINIT) THEN
               WRITE(*,*) "path> ERROR: RIGIDINIT option not coded for UNRES in path.f"
               STOP
            ENDIF
            DO J2=1,nres
               c(1,J2)=Q(6*(J2-1)+1)
               c(2,J2)=Q(6*(J2-1)+2)
               c(3,J2)=Q(6*(J2-1)+3)
               c(1,J2+nres)=Q(6*(J2-1)+4)
               c(2,J2+nres)=Q(6*(J2-1)+5)
               c(3,J2+nres)=Q(6*(J2-1)+6)
            ENDDO
            CALL UPDATEDC
            CALL int_from_cart(.true.,.false.)
            CALL chainbuild
            HORDER=1
            FPGRP='C1'
            IF (ENDNUMHESS) THEN
               CALL MAKENUMINTHESS(NINTS,NATOMS)
               CALL GETSTUFF(KD,NNZ,NINTB)
               CALL INTSECDET(Q,3*NATOMS,KD,NNZ,NINTB,DIAG)
            ELSEIF (.NOT.NOFRQS) THEN
               CALL POTENTIAL(Q,E,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
            ENDIF
         ELSEIF (GTHOMSONT) THEN
            CALL GTHOMSONANGTOC(TMPCOORDS,Q,NATOMS)
            CALL SYMMETRY(HORDER,.FALSE.,TMPCOORDS,INERTIA)
            IF (.NOT.NOFRQS) THEN
               IF (ENDNUMHESS) THEN
                  CALL MAKENUMHESS(Q,NATOMS)
               ELSE
                  CALL POTENTIAL(Q,E,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
               ENDIF
            ENDIF
         ELSE
            IF (RIGIDINIT.AND.(.NOT.NOFRQS)) THEN
              ! sn402: see comment above (in the IF(AMBER) block)
              CALL GENRIGID_EIGENVALUES(Q, ATMASS, DIAG, INFO)
              IF (DIAG(1).LT.DIAG(DEGFREEDOMS)) THEN
                 CALL EIGENSORT_VAL_ASC(DIAG(1:DEGFREEDOMS),HESS(1:DEGFREEDOMS,1:DEGFREEDOMS),DEGFREEDOMS,DEGFREEDOMS)
              ENDIF
              CALL SYMMETRY(HORDER,.FALSE.,Q,INERTIA) ! sn402: I think this needs to be here?
            ELSE
               IF (ENDNUMHESS) THEN
                   CALL MAKENUMHESS(Q,NATOMS)
               ELSEIF (.NOT.NOFRQS) THEN
                   CALL POTENTIAL(Q,E,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
               ENDIF
               IF (VARIABLES) THEN
                   HORDER=1
                   FPGRP='C1'
               ELSE
                   CALL SYMMETRY(HORDER,.FALSE.,Q,INERTIA)
               ENDIF
               IF (.NOT.NOFRQS) CALL MASSWT2(NATOMS,ATMASS,Q,VNEW,.TRUE.)
            ENDIF
         ENDIF

         ! Dump the Hessian to a file, if we're doing that.
         IF (HESSDUMPT) THEN
            LUNIT=GETUNIT()
            OPEN(LUNIT,FILE=TRIM(ADJUSTL(HESSDUMP_FNAME)),STATUS='UNKNOWN',POSITION ='APPEND')
            WRITE(LUNIT,'(6G20.10)') HESS(1:NOPT,1:NOPT)
            CLOSE(LUNIT)
         ENDIF

         ! Diagonalise the Hessian to obtain the squared frequencies.
         ! In the case of NOFRQS, we obviously don't need to do the diagonalisation!
         ! In the case of UNRST and RIGIDINIT, we've already done it.
         IF (.NOT.(UNRST.OR.NOFRQS.OR.RIGIDINIT)) THEN
            CALL DSYEV('N','U',NOPT,HESS,SIZE(HESS,1),DIAG,TEMPA,3*NOPT,INFO)
            IF (DIAG(1).LT.DIAG(NOPT)) CALL EIGENSORT_VAL_ASC(DIAG,HESS,NOPT,NOPT)
         ENDIF


         ! We now have everything we need to write this point to the path.info file:
         !    (i) The energy of the point
         !    (ii) The point group
         !    (iii) The Hessian eigenvalues
         !    (iv) The coordinates
         ! We now proceed to write these to the file.

         ! (i), (ii) Write the first two header lines for this point: The energy of the point,
         ! followed by the point group order and symbol.

         ! (i)
         IF (UNRST .AND. (.NOT. ISTS)) THEN
            IF (CALCDIHE) THEN
                CALL UNRESCALCDIHEREF(DIHE,ALLANG,Q)
            ELSE
                DIHE=0.5D0 ! dummy order param for pathsample related purposes
            ENDIF
            IF(MACHINE) THEN
               WRITE(88) E, DIHE
            ELSE
               WRITE(88,'(2G25.15)') E, DIHE
            ENDIF
         ELSE
            IF(MACHINE) THEN
                WRITE(88) E
            ELSE
                WRITE(88,'(F25.15)') E
            ENDIF
         ENDIF
         ! (ii)
         if (machine) then
              WRITE(88) HORDER,FPGRP
         else
              WRITE(88,'(I6,1X,A4)') HORDER,FPGRP
         endif

         ! (iii) Now write the frequencies. In line with what PATHSAMPLE is expecting, these are squared angular frequencies.
         ! The default units are (rad/s)^2 for AMBER, CHARMM and a few others (see keywords.f), and (rad/[time])^2 for everything
         ! else, where [time] denotes the internal time unit for the potential in question.
         ! FRQCONV can be set manually if you want your square frequencies in a different unit.
! hk286
         IF (.NOT. NOFRQS) THEN
             IF (GTHOMSONT) THEN
                if (machine) then
                   WRITE(88) (1.0D10, J2=1, NATOMS)
                   WRITE(88) (DIAG(J2)*FRQCONV2,J2=1,2*NATOMS)
                else
                   WRITE(88,'(3G20.10)') (1.0D10, J2=1, NATOMS)
                   WRITE(88,'(3G20.10)') (DIAG(J2)*FRQCONV2,J2=1,2*NATOMS)
                endif
             ELSE
                 if (machine) then
                     IF (RIGIDINIT) THEN
                         WRITE(88) (DIAG(J2)*FRQCONV2,J2=1,DEGFREEDOMS)
                     ELSE
                         WRITE(88) (DIAG(J2)*FRQCONV2,J2=1,NOPT)
                     ENDIF
                 else
                     IF (RIGIDINIT) THEN
                         WRITE(88,'(3G20.10)') (DIAG(J2)*FRQCONV2,J2=1,DEGFREEDOMS)
                     ELSE
                         WRITE(88,'(3G20.10)') (DIAG(J2)*FRQCONV2,J2=1,NOPT)
                     ENDIF
                 endif
             ENDIF
         ENDIF

      ENDIF

      ! (iv) Finally, write the coordinates of the point.

      IF (MACHINE) then
          IF (GTHOMSONT) THEN
            CALL GTHOMSONANGTOC(TMPCOORDS, Q, NATOMS)
            WRITE(88) (TMPCOORDS(J2), J2=1, 3*NATOMS)
          ELSE IF (RIGIDINIT) THEN
             CALL TRANSFORMRIGIDTOC (1, NRIGIDBODY, XCOORDS, Q)
             WRITE(88) (XCOORDS(J2),J2=1,3*NATOMS)
          ELSE
             WRITE(88) (Q(J2),J2=1,NOPT)
          ENDIF
      ELSEIF (AMHT) THEN

!  THIS IS FOR PLACE HOLDING C-BETAS FOR GLYCINE IN AMH
         GLY_COUNT = 0

         DO J2=1, NRES_AMH_TEMP
            IF (SEQ(J2).EQ.8) THEN
               WRITE(88,*)Q(9*(J2-1)+1-GLY_COUNT*3),
     &           Q(9*(J2-1)+2-GLY_COUNT*3),Q(9*(J2-1)+3-GLY_COUNT*3)
               WRITE(88,*)Q(9*(J2-1)+1-GLY_COUNT*3),
     &           Q(9*(J2-1)+2-GLY_COUNT*3),Q(9*(J2-1)+3-GLY_COUNT*3)
               WRITE(88,*)Q(9*(J2-1)+4-GLY_COUNT*3),
     &               Q(9*(J2-1)+5-GLY_COUNT*3),Q(9*(J2-1)+6-GLY_COUNT*3)
               GLY_COUNT = GLY_COUNT + 1
            ELSE
              WRITE(88,*)Q(9*(J2-1)+1-GLY_COUNT*3),
     &           Q(9*(J2-1)+2-GLY_COUNT*3),Q(9*(J2-1)+3-GLY_COUNT*3)
              WRITE(88,*)Q(9*(J2-1)+4-GLY_COUNT*3),
     &           Q(9*(J2-1)+5-GLY_COUNT*3),Q(9*(J2-1)+6-GLY_COUNT*3)
             WRITE(88,*)Q(9*(J2-1)+7-GLY_COUNT*3),
     &           Q(9*(J2-1)+8-GLY_COUNT*3),Q(9*(J2-1)+9-GLY_COUNT*3)
            ENDIF
         ENDDO
      ELSE
          IF (GTHOMSONT) THEN
            CALL GTHOMSONANGTOC(TMPCOORDS, Q, NATOMS)
            WRITE(88,'(3F25.15)') (TMPCOORDS(J2),J2=1, 3*NATOMS)
           ELSE IF (RIGIDINIT) THEN
              CALL TRANSFORMRIGIDTOC (1, NRIGIDBODY, XCOORDS, Q)
              WRITE(88,'(3F25.15)') (XCOORDS(J2),J2=1,3*NATOMS)
           ELSE
              WRITE(88,'(3F25.15)') (Q(J2),J2=1,NOPT)
           ENDIF
      ENDIF

      END SUBROUTINE
