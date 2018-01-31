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
      SUBROUTINE GEOPT(FNAMEF,EFNAMEF,Q,VECS,MFLAG,ENERGY,EVALMIN,VNEW)
      USE COMMONS
      USE KEY
      USE MODTWOEND
      USE MODHESS
      USE modcharmm
      use porfuncs
      use binaryio
      USE MODEFOL
      USE MODNEB,ONLY : NEWCONNECTT
      USE AMHGLOBALS, ONLY : NMRES
! hk286
      USE GENRIGID
      use internals_wrapper
      USE MODCUDALBFGS, ONLY : CUDA_LBFGS_WRAPPER
      USE MODCUDABFGSTS, ONLY : CUDA_BFGSTS_WRAPPER
      USE CHIRALITY, ONLY: CIS_TRANS_CHECK, CHIRALITY_CHECK

      IMPLICIT NONE
      INTEGER FCALL, INEG, ECALL, NPCALL, SCALL, ITDONE, INFO, NEXMODES, II1, I1, NUMHB, J1, HORDER, J3, NSTEP, J2, K1, K2,
     &        IMAX, IDONE, ITS 
      DOUBLE PRECISION  DIAG(3*NATOMS), ENERGY, EREAL, EVALMAX, EVALMIN, ETIME, FTIME, STIME, Q(3*NATOMS), OMAX,
     1                  RMS, RMS2, VECS(3*NATOMS), VNEW(3*NATOMS), IT(3,3), IV(3,3), AMASS, ATMASSSAVE(NATOMS), DUMMY,
     2                  EVALUES(3*NATOMS),TEMPA(9*NATOMS), DUMQ(3*NATOMS), ITX, ITY, ITZ, RMSD, RGYR, PROD, DPRAND, TSDISP,
     &                  EVSAVE(3*NATOMS), ESAVE, TSDISPSAVE, MINFRQ2, MINCURVE, GRAD(3*NATOMS), RPBN, QFAC, TMPC(3*NATOMS)
      DOUBLE PRECISION DIFF, ORDERPLUS, ORDERMINUS, DIST, DIST2, RMAT(3,3), QPATH(3*NATOMS)
      DOUBLE PRECISION, ALLOCATABLE :: ORDERDERIV(:),ORDER(:),ORDERSAVE(:)
      DOUBLE PRECISION QSAVE(3*NATOMS)
      CHARACTER(LEN=80) FNAMEF
      CHARACTER(LEN=20) EFNAMEF
      CHARACTER(LEN=22) ITSTRING
      CHARACTER(LEN=5)  ADUMMY
      LOGICAL MFLAG, ZT(3*NATOMS), PERMDISTLOCAL, EWARN
      LOGICAL :: BTEST, LSELECT, LNATIVE, INERTIAT=.TRUE., GOODSTRUC1, GOODSTRUC2
!     AMH LOCAL VARIABLES
      INTEGER :: NRES,I_RES, GLY_COUNT, IPOT
!      CHARACTER(LEN=5) :: TARFL
      CHARACTER(LEN=2) :: SDUMMY 
      REAL(8) QTEMP(3,NATOMS), Q4, Q6

      COMMON /PCALL/ NPCALL, ECALL, FCALL, SCALL, ETIME, FTIME, STIME
      LOGICAL KNOWE, KNOWG, KNOWH
      COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
      LOGICAL CONNECTT, DUMPPATH, READPATH, CALCRATES, STOPFIRST, PROJGRAD
      DOUBLE PRECISION TEMPERATURE, HRED, XIMAGE(GCIMAGE,3*NATOMS), XINITIAL(3*NATOMS+1)
      INTEGER NCONNECT, NEWINR, ISTAT
      COMMON /CONN/ STOPFIRST, CONNECTT, NCONNECT, DUMPPATH, READPATH, CALCRATES, TEMPERATURE, HRED
      DOUBLE PRECISION BHENERGY
      COMMON /BHINTE/ BHENERGY
      DOUBLE PRECISION BISECTENERGY
      COMMON /BISECTE/ BISECTENERGY
      LOGICAL PATHT, DRAGT, LZT(NOPT)
      INTEGER NPATHFRAME
      COMMON /RUNTYPE/ DRAGT, PATHT, NPATHFRAME
C
C  Storage for DSYEVR for LOWESTFRQ calculation
C
      INTEGER IWORK(33*3*NATOMS)
      INTEGER ILWORK, LWORK, NFOUND, ISUPPZ(2*3*NATOMS)
      DOUBLE PRECISION WORK(33*3*NATOMS), ABSTOL, DLAMCH
      DOUBLE PRECISION, ALLOCATABLE :: ZSAVE(:,:), ZWK(:,:)

      DOUBLE PRECISION :: msb50AR(3*natoms)

      INTEGER KD,NNZ,NINTB ! jmc
      DOUBLE PRECISION DIHE,ALLANG ! JMC

      DOUBLE PRECISION ETS,EPLUS,EMINUS
      COMMON /OEPATH/ ETS,EPLUS,EMINUS

! hk286
      DOUBLE PRECISION :: XCOORDS (3*NATOMS)

! cs778
      DOUBLE PRECISION :: SOVER

      INTEGER :: NSTARTHESS ! sn402

      PROD = 0.0D0  ! Initialise here, so it can be used by every other IF block

      INFO=0
      LWORK=33*3*NATOMS
      ILWORK=33*3*NATOMS
      IF (NENDHESS.LE.0) NENDHESS=NOPT
C
C  *************** two-ended pathways ********************
C
      TTDONE=.FALSE.
      IF (TWOENDS) THEN
         CALL TWOEND(ENERGY,VNEW,VECS,Q)
         TTDONE=.TRUE.
      ENDIF
C
C  *************** BFGS minimization ********************
C
      IF ((BHINTERPT.OR.BISECTT).AND.(.NOT.NEWCONNECTT).AND.(.NOT.REOPTIMISEENDPOINTS)) THEN ! do nothing
      ELSE IF ((ENDHESS.OR.ENDNUMHESS).AND.(NSTEPS.EQ.0)) THEN ! no geometry optimisation
         MFLAG=.TRUE.
      ELSE IF (HYBRIDMINT) THEN
         NSTEP=0
         CALL HYBRIDMIN(HMNSTEPS,Q,ENERGY,VNEW,MFLAG,RMS,EVALMIN,EVALMAX,VECS,ITDONE,.TRUE.,.TRUE.)
         NSTEP=ITDONE
         IF (.NOT.MFLAG) THEN
            PRINT '(A,I8,A)',' geopt> switching to LBFGS minimisation after ',NSTEP,' hybrid minimisation steps'
            KNOWE=.FALSE.
            KNOWG=.FALSE.
            IF (INTWRAP_USEINTERNALS()) THEN
               CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,.TRUE.,VNEW,.FALSE.,.FALSE.)
            ELSE
               CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                      .TRUE.,ITDONE,.TRUE.,VNEW,.FALSE.,.FALSE.)
            ENDIF
            NSTEP=NSTEP+ITDONE
         ENDIF
      ELSE IF ((BFGSMINT.AND.(.NOT.BFGSTST).AND.(.NOT.BFGSSTEP).AND. 
     &                    (.NOT.(INR.EQ.1)).AND.(.NOT.(INR.EQ.2))).OR.(BFGSTST.AND.(HINDEX.EQ.0))) THEN
         IF (UNRST.OR.INTWRAP_USEINTERNALS()) THEN
            CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                   .TRUE.,ITDONE,.TRUE.,VNEW,.FALSE.,.FALSE.)
         ELSE
            IF (CUDAT) THEN
               CALL CUDA_LBFGS_WRAPPER(NOPT,MUPDATE,Q,MFLAG,ENERGY,RMS,BFGSSTEPS,ITDONE,.FALSE.)
               WRITE(*,'(A,G25.17)') ' geopt> Final energy is ',ENERGY
            ELSE
               CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                   .TRUE.,ITDONE,.TRUE.,VNEW,.FALSE.,.FALSE.)
            END IF
         ENDIF
C
C  ******** great circle interpolation between end points ********************
C
      ELSE IF (GREATCIRCLET.AND.(.NOT.CONNECTT)) THEN
          IF (UNRST) THEN
             PRINT *,'not available'
             CALL FLUSH(6)
             STOP
          ELSE
             XINITIAL(1:3*NATOMS+1)= 0.0D0
             XINITIAL(3*NATOMS)= 1.0D0
             XINITIAL(3*NATOMS+1)= 1.0D0
             CALL GCLBFGS(Q,FIN,XIMAGE,3*NATOMS+1,GCUPDATE,XINITIAL,.FALSE.,GCCONV,MFLAG,ENERGY,RMS,
     1                   GCSTEPS,.TRUE.,ITDONE,.TRUE.)
          ENDIF
C
C  *************** morphing via coarse-grained EF ****************************
C
      ELSE IF (MORPHT.AND.(.NOT.CONNECTT)) THEN
          IF (UNRST) THEN
             PRINT *,'not available'
             CALL FLUSH(6)
             STOP
C            CALL INTMORPH(MSTEPS,Q,ENERGY,VNEW,MFLAG,RMS,EVALMIN,EVALMAX,VECS,ITDONE,.TRUE.,.TRUE.)
          ELSE
             CALL MORPH(MSTEPS,Q,FIN,ENERGY,VNEW,MFLAG,RMS,ITDONE,.TRUE.)
          ENDIF
C
C  *************** hybrid EF ****************************
C
      ELSE IF (BFGSTST) THEN
          IF (UNRST) THEN
             PRINT '(A)', ' geopt> setting random initial vector for eigenvector'
             DO J1=1,NINTS
                VECS(J1)=DPRAND()*2-1.0D0
             ENDDO
             CALL VECNORM(VECS,NINTS)
             CALL INTBFGSTS(NSTEPS,Q,ENERGY,VNEW,MFLAG,RMS,EVALMIN,EVALMAX,VECS,ITDONE,.TRUE.,.TRUE.) 
          ELSE
             PRINT '(A)', ' geopt> setting random initial vector for eigenvector'
             DO J1=1,NOPT
                VECS(J1)=DPRAND()*2-1.0D0
             ENDDO 
             CALL VECNORM(VECS,NOPT)

             IF (CUDAT) THEN
                 CALL CUDA_BFGSTS_WRAPPER(NSTEPS,Q,ENERGY,MFLAG,RMS,EVALMIN,VECS,ITDONE)
             ELSE
                 CALL BFGSTS(NSTEPS,Q,ENERGY,VNEW,MFLAG,RMS,EVALMIN,EVALMAX,VECS,ITDONE,.TRUE.,.TRUE.)
             END IF
             IF (MFLAG) THEN
                 WRITE(*,'(1X,A,I6)') 'Converged to TS (number of iterations):     ',ITDONE
                 write(*,*) "Energy:", ENERGY
             ELSE
                 WRITE(*,'(1X,A,I6)') 'DID NOT converge to TS (number of iterations):     ',ITDONE
             END IF
             WRITE(*,'(A,F20.10)') ' geopt> Energy of TS is:                    ',ENERGY
          ENDIF

C
C  We may now have the energy and gradient, so the first energy and gradient call in
C  mylbfgs for BFGSSTEP may be unnecessary
C
          IF (BFGSSTEP) THEN
C
C  Switch to appropriate minimization after pushoff.
C
             IF (.NOT.(RKMIN.OR.BSMIN)) BFGSMINT=.TRUE.
             IF (RKMIN) RMS=1.0D0
             MFLAG=.FALSE.
C            NOSHIFT=.TRUE.
             BFGSSTEP=.FALSE.
             BFGSTST=.FALSE.
             IF (.NOT.BFGSMINT) WRITE(*,'(A,26X,F20.10)') ' Energy after pushoff=  ',ENERGY
             IF (BFGSMINT) THEN
                IF (UNRST.OR.INTWRAP_USEINTERNALS()) THEN
                   CALL MYLBFGS(NINTS,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,
     1                       BFGSSTEPS,.TRUE.,ITDONE,.TRUE.,VNEW,.FALSE.,.FALSE.)
                ELSE
                   CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                          .TRUE.,ITDONE,.TRUE.,VNEW,.FALSE.,.FALSE.)
                ENDIF
             ELSE IF (BSMIN.OR.RKMIN) THEN
                CALL ODESD(NSTEPS,Q,MFLAG,ITDONE,.TRUE.)
             ENDIF
          ENDIF

C
C Using xmylbfgs to locate the lowest eigenvalue and the corresponding eigenvector
C
      ELSE IF (EIGENONLY) THEN
         TSTART=TSTART-1.0d0
          PRINT '(A)', ' geopt> locate the lowest eigenvector with lbfgs'
          DO J1=1,NOPT
             VECS(J1)=DPRAND()*2-1.0D0
          ENDDO
          IF (READV) THEN
             PRINT '(A)', ' geopt> read vector.dump as initial vector'
             IF(.NOT.DUMPV) OPEN(44,FILE="vector.dump",STATUS='UNKNOWN')
             REWIND(44)
             READ(44,*,END=111) EVALMIN
             READ(44,*,END=111) (VECS(J1),J1=1,NOPT)
             IF(.NOT.DUMPV) CLOSE(44)
111          CONTINUE
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
          CALL VECNORM(VECS,NOPT)
          SOVER=0.0d0
          CALL BEIG(ITDONE,Q,ENERGY,VECS,EVALMIN,J1,SOVER,.TRUE.,MFLAG)
          IF (DUMPV) THEN
            REWIND(44)
            WRITE(44,'(E20.10)') EVALMIN
            WRITE(44,'(3F20.10)') (VECS(J1),J1=1,NOPT)
            CALL FLUSH(44)
          ENDIF
          
      ELSE IF ((BSMIN.OR.RKMIN).AND.(.NOT.BFGSTST).AND.(.NOT.BFGSSTEP).AND.
     &          (.NOT.(INR.EQ.1)).AND.(.NOT.(INR.EQ.2))) THEN
         CALL ODESD(NSTEPS,Q,MFLAG,ITDONE,.TRUE.)
         IF (.NOT.MFLAG) THEN
            BSMIN=.FALSE.
            RKMIN=.FALSE.
            BFGSMINT=.TRUE.
            KNOWE=.FALSE.
            KNOWG=.FALSE.
            CALL MYLBFGS(NOPT,MUPDATE,Q,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,BFGSSTEPS,
     1                   .TRUE.,ITDONE,.TRUE.,VNEW,.FALSE.,.FALSE.)
         ENDIF
      ELSE
C
C  *************** EF optimization (or Page-McIver second order SD) ***************
C
         NEWINR=INR
         CALL EFOL(Q,MFLAG,NSTEPS,ENERGY,ITDONE,EVALMIN,.TRUE.,DIAG,NEWINR)
      ENDIF
C
C  *************** End of optimisation possibilities. Checkindex if required. *****
C
      IF ((MFLAG.AND.CHECKINDEX).AND.(BFGSMINT.OR.BFGSTST.OR.BSMIN.OR.RKMIN)) THEN
         INEG=0
         IF (BFGSTST) INEG=1
         IF (NOHESS) THEN
            CALL CHECKIND2(Q,MFLAG,INEG,ENERGY)
         ELSE
C
C  We need the Hessian in CHECKIND.
C
            IF (BFGSMINT.OR.BSMIN.OR.RKMIN) CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
            CALL CHECKIND(Q,MFLAG,INEG,ENERGY,EVALMIN,EVALMAX,.FALSE.)
         ENDIF
      ENDIF
!
! If we have just finished the BHINTERPT or BISECTT procedure and geopt
! is being called to make min.data.info MFLAG is not set, and could be false.
!
      IF ((BHINTERPT.OR.BISECTT).AND.(.NOT.NEWCONNECTT).AND.(.NOT.REOPTIMISEENDPOINTS)) MFLAG=.TRUE.
      IF (.NOT.MFLAG) GOTO 11 ! skip the calculations for min.data.info construction

! Set up the required information (i.e. frequencies) for min.data.info to be printed
C
C DAE
C may want to calculate rate from system which has no second derivatives in the potential
C There are two possibilites for doing this
C One is constructing a numerical hessian and diagonalise - remember mass weighting
C Or could make approximation that only the highest (lowest?!) frequency mode is important and use
C the XMYLBFGS routine to give estimate of this for stationary points - again need mass weighting
C
C Alternatively may want to do BFGS minimisation on a system where second derivatives
C are available but only use the sec.der. at the end
C
C In all cases aim to to mimic output of efol of 'Log product ...', so that
C grep in Filthy or pathsample can cope without modification.
C should remove need for paul;s gmfreq and crap (sic!) programs which he introduced to do
C mass weighting
C
C May want to rewrite above CHECKINDEX code to avoid duplication of effort
C

! If we're not attempting to calculate frequencies, set the log product to 1 as a dummy, and do no more.
      IF (NOFRQS) THEN  
         PROD=1.0D0
!
! REOPTIMISEENDPOINTS is set to false after a bhinterp run is finished
! so that second derivatives are only calculated at the end.
!

! ENDHESS should be set if you want to print frequencies to min.data.info
!     ELSEIF (ENDHESS .AND. ((.NOT. (REOPTIMISEENDPOINTS.AND.(BHINTERPT.OR.BISECTT))).AND.DUMPDATAT)) THEN
      ELSEIF (ENDHESS .AND. (.NOT. (REOPTIMISEENDPOINTS.AND.(BHINTERPT.OR.BISECTT)))) THEN

! This first block calculates the Hessian (either analytical or numerical) to extract the frequencies.
! jmc49 skipping the calculation of the Hessian here if RIGIDINIT is true (unless LOWESTFRQT is also true)
! to avoid duplication of effort because GENRIGID_EIGENVALUES will calculate everything later 
       IF ((.NOT. RIGIDINIT) .OR. (RIGIDINIT .AND. LOWESTFRQT)) THEN
         IF (ENDNUMHESS) THEN
            IF (UNRST) THEN
               CALL MAKENUMINTHESS(NINTS,NATOMS)
               WRITE (*,'(A)') ' geopt> Value below calculated from numerical hessian' 
               CALL GETSTUFF(KD,NNZ,NINTB)
               CALL INTSECDET(Q,3*NATOMS,KD,NNZ,NINTB,EVALUES)
               NEXMODES=0
C jmc how many eigenvalues? Use NEXMODES to get right number for internals...
               IF (BFGSMINT) THEN
                  NEXMODES=3*NATOMS-NINTS
               ELSEIF (BFGSTST) THEN
                  NEXMODES=3*NATOMS-NINTS+1
               END IF
            ELSEIF (RINGPOLYMERT) THEN
               CALL MAKENUMHESSRP(Q,NOPT)
               WRITE (*,'(A)') ' geopt> Value below calculated from numerical hessian'
            ELSE
               CALL MAKENUMHESS(Q,NATOMS)
               WRITE (*,'(A)') ' geopt> Value below calculated from numerical hessian'
            ENDIF
         ELSE
C
C might want to do this if analytical hessian is available e.g. after a BFGSMIN run
C
! hk286 - tells potential that we switch to moving frame - needed for normal modes calculations
            RBAANORMALMODET = .TRUE.
            CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)  
            RBAANORMALMODET = .FALSE.
! hk286
            WRITE (*,'(A)') ' geopt> Value below calculated from true hessian'
         ENDIF
       ENDIF

! Next, we diagonalise the Hessian. We need to know how many zero eigenvalues there will be (NEXMODES), 
! so we can ignore these when we come to print output.
         IF (.NOT.UNRST) THEN
C
C dae DSEYV sorts eigenvalues so 6 zero ones should be at bottom
C therefore do product up to 3*NATOMS-6
C this will only work for isolated non-linear systems - fine for charmm
C For general case need code involving ZT from efol
C

! sn402: possible todo: add an IF block here for RIGIDINIT?
            NEXMODES=6
            IF (BFGSMINT.AND.(.NOT.BFGSTST)) THEN
               NEXMODES=6
            ELSEIF (BFGSTST) THEN
               ! If we have a transition state, there will be a negative mode as well as the 6 zeros. Exclude this.
               NEXMODES=7
            ELSEIF (INR.EQ.2) THEN ! this will detect higher order saddles
               NEXMODES=7
            ENDIF
            ! No rotational modes in bulk systems
            IF (BULKT) NEXMODES=NEXMODES-3
            IF (BULKT.AND.TWOD) NEXMODES=NATOMS+2
            IF (PULLT.OR.EFIELDT) NEXMODES=4
! hk286
            IF (GTHOMSONT .AND. (GTHOMMET .EQ. 5) ) NEXMODES = 3
            IF (GTHOMSONT .AND. (GTHOMMET < 5) ) NEXMODES = 1
            IF (TWOD) NEXMODES=NEXMODES+NATOMS
            IF (FREEZE) THEN
               NEXMODES=3*NFREEZE
            ENDIF
            IF (RBAAT) THEN
             IF(UNIAXT.AND..NOT.SANDBOXT) THEN
               IF (EFIELDT) THEN
                  NEXMODES = NATOMS/2 + 4
               ELSE
                  NEXMODES = NATOMS/2 + 6
               ENDIF
             ELSE IF(SANDBOXT) THEN
               NEXMODES=6
               DO J1=1,NATOMS/2
                 IF (UNIAXARRAY(J1)) THEN
                  NEXMODES = NEXMODES + 1
                 END IF
               END DO
             ELSE
               IF (EFIELDT) THEN
                  NEXMODES = 4
               ELSE
                  NEXMODES = 6
               ENDIF
             END IF
               IF (STOCKAAT) NEXMODES = NEXMODES + NATOMS/2
            ENDIF
            IF (RINGPOLYMERT) THEN
               IF (RPSYSTEM(1:4).EQ.'AECK') THEN
                  NEXMODES=0
               ELSE
                  NEXMODES=6
               ENDIF
               IF (BFGSTST.OR.(INR.EQ.2)) NEXMODES=NEXMODES+1
               IF (BFGSMINT) NEXMODES=NEXMODES+2
            ENDIF
            IF (STEALTHYT) THEN
               IF (ENERGY.LT.1.0D-8) THEN
                  NEXMODES=NATOMS*3 - STM*2
               ELSE
                  NEXMODES=3
               ENDIF
            ENDIF
            IF (MIEFT .OR. NOTRANSROTT) NEXMODES=0
            IF (VARIABLES) NEXMODES=NZERO
            IF (BFGSTST) NEXMODES=NEXMODES+1 !sf344> so that the log product for transition state geometries can be calculated and dumped correctly
            WRITE(*,'(A,I6)') ' geopt> Number of zero/imaginary eigenvalues assumed to be ',NEXMODES

! Keyword-specific block
            IF (LOWESTFRQT) THEN
C
C  Calculate lowest non-zero eigenvalue and dump to min.data.info file
C  No mass-weighting here!
C  'U' specifies that the upper triangle contains the Hessian.
C  Need to save HESS before this call and restore afterwards.
C
    
               ABSTOL=DLAMCH('Safe  minimum')
               IF (ALLOCATED(ZWK)) DEALLOCATE(ZWK)
!              ALLOCATE(ZWK(NOPT,NEXMODES+1))
               ALLOCATE(ZWK(1,1))
               ALLOCATE(ZSAVE(NOPT,NOPT))
               ZSAVE(1:NOPT,1:NOPT)=HESS(1:NOPT,1:NOPT)

               IF (RIGIDINIT) THEN
                  RBAANORMALMODET = .TRUE.
                  CALL GENRIGID_EIGENVALUES(Q, ATMASS, EVALUES, INFO)
                  RBAANORMALMODET = .FALSE.
                  WRITE(*,*) "geopt> This combination of keywords (LOWESTFRQ) has not been checked - sn402"
               ELSEIF (RBAAT.AND..NOT.PYGPERIODICT.AND..NOT.SANDBOXT) THEN  ! hk286
                  RBAANORMALMODET = .TRUE.
                  CALL NRMLMD (Q, EVALUES, .TRUE.)
                  RBAANORMALMODET = .FALSE.
                  PRINT *, "geopt> This combination of keywords has not been checked carefully - hk286"
               ELSE
                  CALL DSYEVR('N','I','U',NOPT,HESS,NOPT,0.0D0,1.0D0,1,NEXMODES+1,ABSTOL,
     &                        NFOUND,EVALUES,
     &                        ZWK,NOPT,ISUPPZ,WORK,
     &                        LWORK, IWORK, ILWORK, INFO )
               ENDIF
               MINCURVE=EVALUES(NEXMODES+1)
               DEALLOCATE(ZWK)
               HESS(1:NOPT,1:NOPT)=ZSAVE(1:NOPT,1:NOPT)
               DEALLOCATE(ZSAVE)
               PRINT '(A,G20.10)',' geopt> lowest non-zero positive eigenvalue=',MINCURVE
               DO J1=1,NEXMODES+1
                  PRINT '(I8,G20.10)',J1,EVALUES(J1)
               ENDDO
            ENDIF  ! End of IF(LOWESTFRQT)

!           PRINT '(A,I6,A)',' geopt> Ignoring the ',NEXMODES,' lowest Hessian eigenvalues'
!           WRITE(*,*) 'ATMASS=', ATMASS(:)

! hk286
! Mass-weight the Hessian
! We don't mass-weight rigid body systems here. It's done in the GENRIGID_EIGENVALUES routine (NRMLMD for the old RBAA scheme)
            IF (.NOT. RIGIDINIT) THEN
               IF (CHRMMT) THEN
                  CALL MASSWT2(NATOMS,ATMASS,Q,VNEW,.TRUE.) ! just does the Hessian
               ELSEIF (RINGPOLYMERT) THEN
                  CALL MASSWTRP(NOPT,RPMASSES,RPDOF) ! just does the Hessian
               ELSEIF (RBAAT.OR.VARIABLES) THEN ! hk286
               ELSE
                  PRINT '(A)',' geopt> calling MASSWT' !! DJW
                  CALL MASSWT(NATOMS,ATMASS,Q,VNEW,.TRUE.) ! does the Hessian, coordinates and gradient vector
               ENDIF
            ENDIF

! Keyword-specific block
C
C calling DSEYV with 'N' as we;re not interested in the eigenvectors
C Hard-wired calculation of eigenvalues and eigenvectors for diala system
C in order to get phi and psi numerical derivatives.
C
C In the new scheme we only include the derivative that is largest in magnitude to the
C corresponding order parameter.
C
            IF (ORDERPARAMT) THEN

               IF (RBAAT .OR. RIGIDINIT) THEN
                  WRITE(*,*) "geopt> Warning- ORDERPARAM probably gives the wrong answers when used with rigid body systems"
               ENDIF

C If AMBER or NAB, need to undo the mass weighting  
               IF(.NOT.CHRMMT) THEN
                 DO J1=1,NATOMS
                    AMASS=1/SQRT(ATMASS(J1))
                    J3=3*J1
                    Q(J3-2)=AMASS*Q(J3-2)
                    Q(J3-1)=AMASS*Q(J3-1)
                    Q(J3)=AMASS*Q(J3)
                 ENDDO
               ENDIF
               IDONE=0
               QSAVE(1:3*NATOMS)=Q(1:3*NATOMS)  !save the ts coordinates
               TSDISP=0.0D0
               OPEN(UNIT=9123,FILE='order.info',STATUS='UNKNOWN')
               IF ((NEXMODES.EQ.7).AND.(PATHT)) THEN
                  IF (FILTH.EQ.0) THEN
                     ITSTRING='points.path.xyz'
                  ELSE
                     WRITE(ITSTRING,'(A)') 'points.path.xyz.'//TRIM(ADJUSTL(FILTHSTR))
                  ENDIF
                  OPEN(UNIT=9124,FILE=ITSTRING,STATUS='UNKNOWN')
               ENDIF
               IF (PATHT) ALLOCATE(ORDERSAVE(NORDER))
333            ALLOCATE(ORDERDERIV(NORDER),ORDER(NORDER))
               CALL DSYEV('V','U',3*NATOMS,HESS,3*NATOMS,EVALUES,TEMPA,9*NATOMS,INFO)
               IF (EVALUES(1).LT.EVALUES(3*NATOMS)) CALL EIGENSORT_VAL_ASC(EVALUES,HESS,3*NATOMS,3*NATOMS)
               IF (IDONE.EQ.0) THEN
                  ESAVE=ENERGY
                  EVSAVE(1:3*NATOMS)=EVALUES(1:3*NATOMS)
                  VECS(1:3*NATOMS)=HESS(1:3*NATOMS,3*NATOMS)
               ENDIF
               PRINT '(A)',' geopt> projected Hessian eigenvalues:'
               PRINT '(6G20.10)',EVALUES(1:3*NATOMS)
               DIFF=1.0D-2
               IF (PATHT) THEN
                  DO I1=1,NORDER
                     IF ((CHRMMT).AND.(WHICHORDER(I1).EQ.'DIHE')) THEN
                        CALL GETDIHE(Q,ORDER(I1),ORDERNUM(I1))
                     ELSEIF ((AMBERT.OR.NABT).AND.(WHICHORDER(I1).EQ.'DIHE')) THEN
                        IF (I1.EQ.1) CALL AMBERDIHEDR(Q,NATOMS,5,7,9,15,ORDER(I1))
                        IF (I1.EQ.2) CALL AMBERDIHEDR(Q,NATOMS,7,9,15,17,ORDER(I1))
                     ENDIF
                  ENDDO
                  QPATH(1:3*NATOMS)=Q(1:3*NATOMS) ! save the sd path configuration
                  IF (IDONE.EQ.0) THEN
                     ORDERSAVE(1:NORDER)=ORDER(1:NORDER)
                     ITS=0
                  ELSEIF (IDONE.GT.0) THEN
                     DUMMY=0.D0
                     DO I1=1,NORDER
                        DUMMY=DUMMY+(ORDERSAVE(I1)-ORDER(I1))**2
                     ENDDO
                     IF (SQRT(DUMMY).LT.0.1D0) ITS=1 
                  ENDIF
                  PRINT '(A,6G20.10)','order parameters with path configuration:        ',ORDER(1:NORDER) 
C TSDISP calculation for PATHT only as a check for how much the Cartesian steepest descent path
C diverges from the reaction path 
                  TSDISPSAVE=TSDISP
                  DO J1=1,2
                     TSDISP=TSDISPSAVE*(-1.D0)**(J1+ITS)
                     DO J2=1,NATOMS
                        Q(3*(J2-1)+1)=QSAVE(3*(J2-1)+1)+(TSDISP*VECS(3*(J2-1)+1))/SQRT(ATMASS(J2))
                        Q(3*(J2-1)+2)=QSAVE(3*(J2-1)+2)+(TSDISP*VECS(3*(J2-1)+2))/SQRT(ATMASS(J2))
                        Q(3*(J2-1)+3)=QSAVE(3*(J2-1)+3)+(TSDISP*VECS(3*(J2-1)+3))/SQRT(ATMASS(J2))
                     ENDDO
                     DO I1=1,NORDER
                        IF ((CHRMMT).AND.(WHICHORDER(I1).EQ.'DIHE')) THEN
                           CALL GETDIHE(Q,ORDER(I1),ORDERNUM(I1))
                        ELSEIF ((AMBERT.OR.NABT).AND.(WHICHORDER(I1).EQ.'DIHE')) THEN
                           IF (I1.EQ.1) CALL AMBERDIHEDR(Q,NATOMS,5,7,9,15,ORDER(I1))
                           IF (I1.EQ.2) CALL AMBERDIHEDR(Q,NATOMS,7,9,15,17,ORDER(I1))
                        ENDIF
                     ENDDO
                     PRINT '(A,7G20.10)','TSDISP and order parameters: ',TSDISP,ORDER(1:NORDER) 
                  ENDDO
C use steepest descent configurations for derivative calculations
                  TSDISP=0.D0
                  Q(1:3*NATOMS)=QPATH(1:3*NATOMS)
                  DO I1=1,NORDER
                     IF ((CHRMMT).AND.(WHICHORDER(I1).EQ.'DIHE')) THEN
                        CALL GETDIHE(Q,ORDER(I1),ORDERNUM(I1))
                     ELSEIF ((AMBERT.OR.NABT).AND.(WHICHORDER(I1).EQ.'DIHE')) THEN
                        IF (I1.EQ.1) CALL AMBERDIHEDR(Q,NATOMS,5,7,9,15,ORDER(I1))
                        IF (I1.EQ.2) CALL AMBERDIHEDR(Q,NATOMS,7,9,15,17,ORDER(I1))
                     ENDIF
                  ENDDO
               ELSE
                  DO J2=1,NATOMS
                     Q(3*(J2-1)+1)=QSAVE(3*(J2-1)+1)+(TSDISP*VECS(3*(J2-1)+1))/SQRT(ATMASS(J2))
                     Q(3*(J2-1)+2)=QSAVE(3*(J2-1)+2)+(TSDISP*VECS(3*(J2-1)+2))/SQRT(ATMASS(J2))
                     Q(3*(J2-1)+3)=QSAVE(3*(J2-1)+3)+(TSDISP*VECS(3*(J2-1)+3))/SQRT(ATMASS(J2))
                  ENDDO
                  DO I1=1,NORDER
                     IF ((CHRMMT).AND.(WHICHORDER(I1).EQ.'DIHE')) THEN
                        CALL GETDIHE(Q,ORDER(I1),ORDERNUM(I1))
                     ELSEIF ((AMBERT.OR.NABT).AND.(WHICHORDER(I1).EQ.'DIHE')) THEN
                        IF (I1.EQ.1) CALL AMBERDIHEDR(Q,NATOMS,5,7,9,15,ORDER(I1))
                        IF (I1.EQ.2) CALL AMBERDIHEDR(Q,NATOMS,7,9,15,17,ORDER(I1))
                     ELSE 
!
! align_decide selects a subroutine that will align the second coordinate array with the first. We don't want to change Q
! so choose the order carefully.
!
                        IF (I1.EQ.1) THEN
                           CALL ALIGN_DECIDE(Q,POINTSDECA,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,ORDER(1),DIST2,
     &                                                RIGIDBODY,RMAT)
                        ELSE IF (I1.EQ.2) THEN
                           CALL ALIGN_DECIDE(Q,POINTSICOS,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,ORDER(2),DIST2,
     &                                                RIGIDBODY,RMAT)
                        ELSE IF (I1.GT.2) THEN
                           DO K1=1,3
                              DO K2=1,NATOMS
                                 QTEMP(K1,K2)=Q(3*(K2-1)+K1)
                              ENDDO
                           ENDDO
                           CALL QORDER_LJ(QTEMP,Q4,Q6)
                           ORDER(3)=Q6
!                          ORDER(4)=Q6
                        ENDIF
                     ENDIF
                  ENDDO
                  IF (DEBUG) THEN
                     PRINT '(A)','         E                 dist deca           dist icos             Q6 unperturbed'
                     PRINT '(5G20.10)',ENERGY,ORDER(1:3)
                  ENDIF
               ENDIF 
               ORDERDERIV(1:NORDER)=0.D0
               DO J1=1,3*NATOMS ! cycle over normal modes
                  OMAX=0.0D0
                  IMAX=1
                  DO I1=1,NORDER ! cycle over order parameters
                     DO J2=1,NATOMS
                        Q(3*(J2-1)+1)=QSAVE(3*(J2-1)+1)+(TSDISP*VECS(3*(J2-1)+1)+DIFF*HESS(3*(J2-1)+1,J1))/SQRT(ATMASS(J2))
                        Q(3*(J2-1)+2)=QSAVE(3*(J2-1)+2)+(TSDISP*VECS(3*(J2-1)+2)+DIFF*HESS(3*(J2-1)+2,J1))/SQRT(ATMASS(J2))
                        Q(3*(J2-1)+3)=QSAVE(3*(J2-1)+3)+(TSDISP*VECS(3*(J2-1)+3)+DIFF*HESS(3*(J2-1)+3,J1))/SQRT(ATMASS(J2))
                     ENDDO
                     IF ((CHRMMT).AND.(WHICHORDER(I1).EQ.'DIHE')) THEN
                        CALL GETDIHE(Q,ORDERPLUS,ORDERNUM(I1))
                     ELSEIF ((AMBERT.OR.NABT).AND.(WHICHORDER(I1).EQ.'DIHE')) THEN
                        IF (I1.EQ.1) CALL AMBERDIHEDR(Q,NATOMS,5,7,9,15,ORDERPLUS)
                        IF (I1.EQ.2) CALL AMBERDIHEDR(Q,NATOMS,7,9,15,17,ORDERPLUS)
                     ELSE 
!
! Fix the permutation to save time and avoid discontinuities.
!
                        PERMDISTLOCAL=PERMDIST
                        PERMDIST=.FALSE.
                        IF (I1.EQ.1) THEN
                           CALL ALIGN_DECIDE(POINTSDECA,Q,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,ORDERPLUS,DIST2,
     &                                                RIGIDBODY,RMAT)
                        ELSE IF (I1.EQ.2) THEN
                           CALL ALIGN_DECIDE(POINTSICOS,Q,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,ORDERPLUS,DIST2,
     &                                                RIGIDBODY,RMAT)
                        ELSE IF (I1.GT.2) THEN
                           DO K1=1,3
                              DO K2=1,NATOMS
                                 QTEMP(K1,K2)=Q(3*(K2-1)+K1)
                              ENDDO
                           ENDDO
                           CALL QORDER_LJ(QTEMP,Q4,Q6)
                           IF (I1.EQ.3) ORDERPLUS=Q6
!                          IF (I1.EQ.4) ORDERPLUS=Q6
                        ENDIF
                        PERMDIST=PERMDISTLOCAL
                     ENDIF
                     DO J2=1,NATOMS
                        Q(3*(J2-1)+1)=QSAVE(3*(J2-1)+1)+(TSDISP*VECS(3*(J2-1)+1)-DIFF*HESS(3*(J2-1)+1,J1))/SQRT(ATMASS(J2))
                        Q(3*(J2-1)+2)=QSAVE(3*(J2-1)+2)+(TSDISP*VECS(3*(J2-1)+2)-DIFF*HESS(3*(J2-1)+2,J1))/SQRT(ATMASS(J2))
                        Q(3*(J2-1)+3)=QSAVE(3*(J2-1)+3)+(TSDISP*VECS(3*(J2-1)+3)-DIFF*HESS(3*(J2-1)+3,J1))/SQRT(ATMASS(J2))
                     ENDDO
                     IF ((CHRMMT).AND.(WHICHORDER(I1).EQ.'DIHE')) THEN
                        CALL GETDIHE(Q,ORDERMINUS,ORDERNUM(I1))
                     ELSEIF ((AMBERT.OR.NABT).AND.(WHICHORDER(I1).EQ.'DIHE')) THEN
                        IF (I1.EQ.1) CALL AMBERDIHEDR(Q,NATOMS,5,7,9,15,ORDERMINUS)
                        IF (I1.EQ.2) CALL AMBERDIHEDR(Q,NATOMS,7,9,15,17,ORDERMINUS)
                     ELSE
!
! Fix the permutation to save time and avoid discontinuities.
!
                        PERMDISTLOCAL=PERMDIST
                        PERMDIST=.FALSE.
                        IF (I1.EQ.1) THEN
                           CALL ALIGN_DECIDE(POINTSDECA,Q,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,ORDERMINUS,DIST2,
     &                                                RIGIDBODY,RMAT)
                        ELSE IF (I1.EQ.2) THEN
                           CALL ALIGN_DECIDE(POINTSICOS,Q,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,ORDERMINUS,DIST2,
     &                                                RIGIDBODY,RMAT)
                        ELSE IF (I1.GT.2) THEN
                           DO K1=1,3
                              DO K2=1,NATOMS
                                 QTEMP(K1,K2)=Q(3*(K2-1)+K1)
                              ENDDO
                           ENDDO
                           CALL QORDER_LJ(QTEMP,Q4,Q6)
                           IF (I1.EQ.3) ORDERMINUS=Q6
!                          IF (I1.EQ.4) ORDERMINUS=Q6
                        ENDIF
                        PERMDIST=PERMDISTLOCAL
                     ENDIF
                     IF (J1.LE.3*NATOMS-NEXMODES) THEN
                        DUMMY=((ORDERPLUS-ORDERMINUS)/(2*DIFF))**2/EVALUES(J1)
                        IF (ABS(DUMMY).GT.ABS(OMAX)) THEN
                           OMAX=DUMMY
                           IMAX=I1
                        ENDIF
!
! This is what we need to calculate the full A_\alpha factor for all modes in LJ75.
!
                         ORDERDERIV(I1)=ORDERDERIV(I1)+DUMMY
                        IF (DEBUG.AND.(I1.EQ.1)) PRINT '(A,2I6,3G20.10)',
     &                              ' geopt> mode,oparam,eigenvalue,deriv contribution,total: ',J1,I1,EVALUES(J1),
     &                              ((ORDERPLUS-ORDERMINUS)/(2*DIFF))**2/EVALUES(J1),ORDERDERIV(I1)
                     ENDIF
!                    PRINT '(A,2I8,2G20.10,I8,G20.10)','J1,I1,ORDERPLUS,ORDERMINUS,IMAX,ORDERDERIV(I1)=',
!    &                                                  J1,I1,ORDERPLUS,ORDERMINUS,IMAX,ORDERDERIV(I1)
                  ENDDO
!
! This is what we have been using for projection onto two order parameters.
!
!                 ORDERDERIV(IMAX)=ORDERDERIV(IMAX)+OMAX

!                 IF ((J1.LE.3*NATOMS-NEXMODES).AND.(DEBUG)) PRINT '(A,I6,A,G20.10,A,I6,A,G20.10)',' geopt> for mode ',J1,
!    &                                     ' maximum contribution is ',OMAX,
!    &                                     ' for mode ',IMAX,' total=',ORDERDERIV(IMAX)
!                 PRINT '(A,I8,A,G20.10)','after mode ',J1,' ORDERDERIV(1)=',ORDERDERIV(1)
               ENDDO
!              PRINT '(A,G20.10)','after mode loop ORDERDERIV(1)=',ORDERDERIV(1)
C
C  Print information for FES calculation here.
C
               PROD=0.0D0
               EWARN=.FALSE.
               DO I1=1,NENDHESS-NEXMODES
                  IF (I1.GT.1) THEN
                     IF (EVALUES(I1-1).NE.0.0D0) THEN
                        IF (ABS(EVALUES(I1)/EVALUES(I1-1)).LT.1.0D-2) THEN
                           PRINT '(A,G20.10,A,G20.10)',' geopt> WARNING - decrease in magnitude of eigenvalues from ',EVALUES(I1-1),
     &                                    ' to ',EVALUES(I1)
                           PRINT '(A)',' geopt> WARNING - this could indicate a stationary point of the wrong index'
                           EWARN=.TRUE.
                        ENDIF
                     ENDIF
                  ENDIF
                  IF (EVALUES(I1).GT.0.0D0) THEN
                     PROD=PROD+DLOG(EVALUES(I1))
                  ELSE
                     IF (I1.LT.(NENDHESS-NEXMODES)) PRINT *,'Higher order saddle detected: eigenvalue ',EVALUES(I1)
                     ! jmc put in this test mainly for pathsample purposes...
                  ENDIF
               ENDDO
               IF (CHRMMT.OR.AMBERT.OR.NABT.OR.AMBER12T.OR.SDT.OR.TTM3T) THEN
C
C if charmm need to convert this to (radian/s)^2, rather than charmm unit
C conversion factor for this is 4.184 x 10^26
C same for AMBER and for Stillinger-David and TTM3.
C
                  PROD=PROD+(NENDHESS-NEXMODES)*DLOG(4.184D26)
                  WRITE (*,'(A,G20.10)') ' geopt> Scaling product of eigenvalues to SI units (radian/s)^2 by ',
     &                                  (3*NATOMS-NEXMODES)*DLOG(4.184D26)
               ENDIF
               IF (BOWMANT) THEN
C
C If Bowman need to convert this to (radian/s)^2
C conversion factor for this is 2625.47 x 10^26
C
                  PROD=PROD+(NENDHESS-NEXMODES)*DLOG(2625.47D26)
                  WRITE (*,'(A,G20.10)') ' geopt> Scaling product of eigenvalues to SI units (radian/s)^2 by ',
     &                                  (3*NATOMS-NEXMODES)*DLOG(2625.47D26)
               ENDIF
               IF (PATHT) TSDISP=TSDISPSAVE
               WRITE(9123,'(3G20.10)') ENERGY,PROD,ABS(TSDISP)
               WRITE(*,'(A,9G20.10)') 'summary: ',ENERGY,ORDER(1:NORDER),ORDERDERIV(1:NORDER)
               DO I1=1,NORDER
                  PRINT '(A,I6,A,2G20.10)',' geopt> order parameter ',I1,' and derivative term: ',ORDER(I1),ORDERDERIV(I1)
                  WRITE(9123,'(2G20.10)') ORDER(I1),ORDERDERIV(I1)
               ENDDO
               DEALLOCATE(ORDER,ORDERDERIV)
               IDONE=IDONE+1
               IF ((NEXMODES.EQ.7).AND.(IDONE.LT.1000)) THEN
                  IF (PATHT) THEN
                     READ(9124,*,END=444) 
                     READ(9124,*,END=444) 
                     DO J2=1,NATOMS
                        READ(9124,*,END=444) ADUMMY,Q(3*(J2-1)+1),Q(3*(J2-1)+2),Q(3*(J2-1)+3)
                     ENDDO
                     CALL ALIGN_DECIDE(Q,QSAVE,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
                     DUMMY=0.D0
                     DO J2=1,NATOMS
                        DUMMY=DUMMY+ATMASS(J2)*( (Q(3*(J2-1)+1)-QSAVE(3*(J2-1)+1))**2 
     &                                          +(Q(3*(J2-1)+2)-QSAVE(3*(J2-1)+2))**2
     &                                          +(Q(3*(J2-1)+3)-QSAVE(3*(J2-1)+3))**2 )
C                        DUMMY=DUMMY+( (Q(3*(J2-1)+1)-QSAVE(3*(J2-1)+1))**2
C     &                                          +(Q(3*(J2-1)+2)-QSAVE(3*(J2-1)+2))**2
C     &                                          +(Q(3*(J2-1)+3)-QSAVE(3*(J2-1)+3))**2 )
                     ENDDO
                     TSDISP=SQRT(DUMMY)
                     WRITE(9125,'(A,I5,F10.2)') 'IDONE, DIST: ',IDONE,DIST
                     DO J2=1,NATOMS
                        WRITE(9125,'(3F10.2,A,3F10.2)') Q(3*(J2-1)+1),Q(3*(J2-1)+2),Q(3*(J2-1)+3),'   ',
     &                  QSAVE(3*(J2-1)+1),QSAVE(3*(J2-1)+2),QSAVE(3*(J2-1)+3)
                     ENDDO
                     WRITE(9125,*)'TSDISP= ',TSDISP
                     WRITE(9125,*)
C                     PRINT *,'TSDISP= ',TSDISP
                  ELSE
                     TSDISP=1.0D-1*((IDONE+1)/2)
                     IF (MOD(IDONE,2).EQ.0) TSDISP=-TSDISP
                     DO J2=1,NATOMS
                        Q(3*(J2-1)+1)=QSAVE(3*(J2-1)+1)+TSDISP*VECS(3*(J2-1)+1)/SQRT(ATMASS(J2))
                        Q(3*(J2-1)+2)=QSAVE(3*(J2-1)+2)+TSDISP*VECS(3*(J2-1)+2)/SQRT(ATMASS(J2))
                        Q(3*(J2-1)+3)=QSAVE(3*(J2-1)+3)+TSDISP*VECS(3*(J2-1)+3)/SQRT(ATMASS(J2))
                     ENDDO
                  ENDIF
                  PRINT '(A,G20.10)',' geopt> displacing TS geometry along reaction coordinate by ',TSDISP
                  IF (ENDNUMHESS) THEN 
C call to potential to get the energy, then numerical Hessian
                     CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
                     WRITE (*,'(A)') ' geopt> Calculating from numerical hessian'
                     IF (UNRST) THEN 
                        CALL MAKENUMINTHESS(NINTS,NATOMS)
                        CALL GETSTUFF(KD,NNZ,NINTB)
                        CALL INTSECDET(Q,3*NATOMS,KD,NNZ,NINTB,EVALUES)
                     ELSEIF (RINGPOLYMERT) THEN
                        CALL MAKENUMHESSRP(Q,NOPT)
                        WRITE (*,'(A)') ' geopt> Value below calculated from numerical hessian'
                     ELSE
                        CALL MAKENUMHESS(Q,NATOMS)
                     ENDIF
                  ELSE
                     CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
                     WRITE (*,'(A,2G20.10,A,G20.10)') ' geopt> Calculating analytical hessian, energy and RMS=',ENERGY,RMS,
     &                                       ' energy change=',ENERGY-ESAVE
                  ENDIF
                  IF (CHRMMT) THEN
                     CALL MASSWT2(NATOMS,ATMASS,Q,VNEW,.TRUE.)  ! just does the Hessian
                  ELSEIF (RINGPOLYMERT) THEN
                     PRINT '(A)','geopt> ERROR *** ring polymer incomaptible with order parameters'
                     STOP
                  ELSE
                     CALL MASSWT(NATOMS,ATMASS,Q,VNEW,.TRUE.)   ! Hessian, coordinates and gradient vector
                  ENDIF
                  PROJGRAD=.TRUE.
                  IF (RMS.LT.CONVR) PROJGRAD=.FALSE.
                  IF (VARIABLES) THEN
                     PRINT '(A)','geopt> Not projecting the Hessian'
                  ELSE
                     CALL PROJH(Q,NATOMS,ATMASS,VNEW,PROJGRAD)
                  ENDIF
                  IF (PATHT) THEN
                     GOTO 333
                  ELSE
                     IF (ENERGY-ESAVE.LT.0.0D0) THEN
                        GOTO 333
                     ELSE
                        GOTO 444
                     ENDIF
                  ENDIF
               ENDIF
444            CONTINUE
               CLOSE(UNIT=9123)
               CLOSE(UNIT=9124)
               ENERGY=ESAVE
               EVALUES(1:3*NATOMS)=EVSAVE(1:3*NATOMS)
!!!!!!!!!!!!!!!!!!!!!!!
!
! End of order parameter block
!
!!!!!!!!!!!!!!!!!!!!!!!

! Back into the main branch of the subroutine.
! This is where we actually calculate the eigenvectors of the (mass-weighted) Hessian.
            ELSE
               IF (DUMPV) THEN ! use job type 'V' to get eigenvectors
! hk286
                  IF (RIGIDINIT) THEN
                     RBAANORMALMODET = .TRUE.
                     ! Eigenvectors will be calculated because DUMPV = .TRUE.
                     CALL GENRIGID_EIGENVALUES(Q, ATMASS, EVALUES, INFO)
                     RBAANORMALMODET = .FALSE.
                  ELSEIF (RBAAT.AND..NOT.PYGPERIODICT.AND..NOT.SANDBOXT) THEN
                     RBAANORMALMODET = .TRUE.
                     ! Last argument to NRMLMD is EIGENVECTORT. Set this to .TRUE. so eigenvectors are calculated
                     CALL NRMLMD (Q, EVALUES, .TRUE.)
                     RBAANORMALMODET = .FALSE.
                  ELSE
                     CALL DSYEV('V','U',NOPT,HESS,NOPT,EVALUES,TEMPA,9*NATOMS,INFO)
                  ENDIF
               ELSE   ! Then we don't need to calculate the eigenvectors. Eigenvalues only.
! hk286 - optional number of ENDHESS not yet implemented
                  IF (RIGIDINIT) THEN
                     RBAANORMALMODET = .TRUE.
                     CALL GENRIGID_EIGENVALUES(Q, ATMASS, EVALUES, INFO)
                     RBAANORMALMODET = .FALSE.
                  ELSE
                     IF ((NENDHESS.GE.NOPT).OR.(VARIABLES.AND.(NENDHESS.GE.NATOMS))) THEN
C                       CALL DSYEV('N','U',NOPT,HESS,NOPT,EVALUES,TEMPA,9*NATOMS,INFO)
C  csw34> changed the call used here to be the same as if NENDHESS < NOPT below
C         as this leads to a 20% speed increase!
                        ABSTOL=DLAMCH('Safe  minimum')
                        ALLOCATE(ZWK(1,1)) ! not referenced for job type 'N'
                        IF (RBAAT.AND..NOT.PYGPERIODICT.AND..NOT.SANDBOXT) THEN ! hk286
                           RBAANORMALMODET = .TRUE.
                           CALL NRMLMD (Q, EVALUES, .TRUE.)
                           RBAANORMALMODET = .FALSE.
                        ELSE
                           CALL DSYEVR('N','I','U',NOPT,HESS,NOPT,0.0D0,1.0D0,1,NOPT,ABSTOL,
     &                          NFOUND,EVALUES,ZWK,NOPT,ISUPPZ,WORK,LWORK,IWORK,ILWORK,INFO )
                        ENDIF
                        DEALLOCATE(ZWK)
                     ELSE
                        ABSTOL=DLAMCH('Safe  minimum')
                        ALLOCATE(ZWK(1,1)) ! not referenced for job type 'N'
                        IF (RBAAT.AND..NOT.PYGPERIODICT.AND..NOT.SANDBOXT) THEN ! hk286
                           RBAANORMALMODET = .TRUE.
                           CALL NRMLMD (Q, EVALUES, .TRUE.)
                           RBAANORMALMODET = .FALSE.
                        ELSE
                           CALL DSYEVR('N','I','U',NOPT,HESS,NOPT,0.0D0,1.0D0,1,NENDHESS,ABSTOL,
     &                          NFOUND,EVALUES,ZWK,NOPT,ISUPPZ,WORK,LWORK,IWORK,ILWORK,INFO )
                        ENDIF
                        DEALLOCATE(ZWK)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF

! By this point we have an array EVALUES containing the eigenvalues (probably ordered from smallest to largest) and
! if VDUMP is set, the matrix HESS contains the corresponding eigenvectors as its columns.

C
C  MASSWT2 and MASSWTRP do not mass weight Q and VNEW, but MASSWT does. Need to undo this
C  if DUMPDATAT is .TRUE. for comparison with pathsample (which uses unit masses).
C  Probably best to always undo it, in case we need non-mass-weighted Q somewhere
C  further down.
C
C           IF (DUMPDATAT.AND.(.NOT.(CHRMMT.OR.RINGPOLYMERT))) THEN
            IF (.NOT. (CHRMMT .OR. RBAAT)) THEN ! hk286
               IF (.NOT.(RIGIDINIT.OR.VARIABLES)) THEN ! jmc49 bugfix, because the mass weighting is not done if rigidinit 
                  DO J1=1,NATOMS
                     AMASS=1/SQRT(ATMASS(J1))
                     J3=3*J1
                     Q(J3-2)=AMASS*Q(J3-2)
                     Q(J3-1)=AMASS*Q(J3-1)
                     Q(J3)=AMASS*Q(J3)
                  ENDDO
               ENDIF
            ENDIF

            ! This orders the eigenvalues from largest to smallest (despite the name of the subroutine...)
            if (evalues(1).lt.evalues(NENDHESS)) call eigensort_val_asc(evalues,hess,NENDHESS,3*natoms)
            IF (INFO.NE.0) PRINT*,'WARNING - INFO=',INFO,' in DSYEV'

         ENDIF   ! End of IF(.NOT. UNRST) (That was a very long block!)
!
! Shift zero eigenvalues for Thomson problem. DJW
!
         IF (GTHOMSONT) THEN
            EVALUES(2*NATOMS+1:3*NATOMS)=1.0D10
            CALL EIGENSORT_VAL_ASC(EVALUES,HESS,NENDHESS,3*NATOMS)
         ENDIF


         IF (CASTEP.OR.ONETEP) THEN
C
C  Added transformation back to Cartesian basis for Hessian eigenvectors,
C  as in "Energy Landscapes" equation (2.51). Otherwise the eigenvector
C  components refer to mass-weighted coordinates, not Cartesians. DJW 7/11/09
C  The eigenvectors of the mass-weighted Hessian correspond to the A matrix
C  components A_{alpha gamma} for eigenvector gamma, and these vectors
C  are orthonormal.
C  Second index of HESS labels the eigenvector, first index runs over components.
C  The transformed eigenvectors in the Cartesian basis are not orthogonal.
C
C  The relative Cartesian displacements for mass-weighted Hessian eigenvector
C  gamma are A_(alpha gamma}/sqrt(m_alpha) where m_alpha is the mass of the
C  atom with component alpha. These are also the relative displacements for
C  atoms corresponding to motion in mode gamma.
C  To put ke of k_gamma into mode gamma choose the Cartesian velocity
C  components as +/- sqrt(2k_gamma) A_{alpha gamma}/sqrt(m_alpha).
C
            DO J1=1,NATOMS ! sum over components
               AMASS=1/SQRT(ATMASS(J1))
               J3=3*J1
               DO J2=1,3*NATOMS ! sum over eigenvectors
                  HESS(J3-2,J2)=HESS(J3-2,J2)*AMASS
                  HESS(J3-1,J2)=HESS(J3-1,J2)*AMASS
                  HESS(J3  ,J2)=HESS(J3  ,J2)*AMASS
               ENDDO
            ENDDO
            PRINT '(A)','geopt> Normalised eigenvectors of mass-weighted Hessian have been transformed to Cartesian components'

         ELSEIF (RINGPOLYMERT.AND.PATHT) THEN
C
C  For ring polymer TS calculate the quantum instanton Im F rate constants for the forward and backward
C  processes. The RP Hessian should have been mass weighted appropriately.
C  Assume mass, length and energy in atomic units, otherwise we need a conversion factor
C  for hbar. PROD already contains the ln product of positive Hessian eigenvalues.
C  Each eigenvalue is an angular frequency squared.
C  We need the energies of the + amd - minima as well, so a path calculation is required.
C
            QFAC=LOG(RPIMAGES/RPBETA)
            IF (RPIMAGES.GT.1) THEN
               DUMMY=0.0D0
               DO J1=1,RPDOF
                  DUMMY=DUMMY+LOG(RPMASSES(J1))
               ENDDO
               DUMMY=DUMMY/RPDOF ! ln of geometric mean mass for RPDOF degrees of freedom
               RPBN=0.0D0
               DO J2=1,RPDOF ! images 1 and RPIMAGES
                  RPBN=RPBN+(Q(J2)-Q(RPDOF*(RPIMAGES-1)+J2))**2
               ENDDO
               DO J1=1,RPIMAGES-1 ! images J1 and J1+1
                  DO J2=1,RPDOF
                     RPBN=RPBN+(Q(RPDOF*(J1-1)+J2)-Q(RPDOF*J1+J2))**2
                  ENDDO
               ENDDO
!
!  This isn't right - the formula in section V has reciprocal factors of g
!  in the frequencies as well. Need to check further.
!
               QFAC=QFAC+(0.5D0)*DUMMY+0.5D0*LOG(RPBN*RPIMAGES/(6.283185307D0*RPBETA))
     &                           -0.5D0*PROD-(RPIMAGES-2)*LOG(RPBETA/RPIMAGES)
            ELSE
            ENDIF
            PRINT '(2(A,G20.10))',' geopt> ln(k_instanton^+ * Q^+)=',QFAC-(ETS-EPLUS)*RPBETA/RPIMAGES,' E+=',EPLUS
            PRINT '(2(A,G20.10))',' geopt> ln(k_instanton^- * Q^-)=',QFAC-(ETS-EMINUS)*RPBETA/RPIMAGES,' E-=',EMINUS
         ENDIF

         ! To avoid lots of messy IF blocks, we define the following variable to help us skip over dummy values with RIGIDINIT
         IF(RIGIDINIT) THEN
             ! This is the first entry in EVALUES which does not correspond to a dummy value
             NSTARTHESS = NENDHESS-DEGFREEDOMS+1
         ELSE
             NSTARTHESS = 1  ! No dummy values, so we start reading from the beginning
         ENDIF

         ! Do some intermediate output printing to show the computed eigenvalues
         ! sn402: I decided not to completely remove the hard-coded unit conversions here, to avoid breaking backwards compatibility
         IF(DEBUG.OR.AMHT.OR.CASTEP.OR.VASP.OR.QCHEM.OR.RINGPOLYMERT.OR.ONETEP.OR.TTM3T.OR.SDT.OR.BOWMANT.OR.CHRMMT) THEN
! hk286
            IF (ONETEP) THEN
               PRINT '(A,I6,A)',' geopt> Frequencies in wavenumbers, assuming hartree and bohr units for input:'
               DO J1=NSTARTHESS,NENDHESS
                  IF (EVALUES(J1).GT.0.0D0) THEN
                     PRINT '(I6,2G20.10)',J1,SQRT(EVALUES(J1)*9.3757D29)/(2*3.141592654D0),
     &                                      SQRT(EVALUES(J1)*9.3757D29)/(2*3.141592654D0*2.998D10)
                  ELSE
                     PRINT '(I6,2(G20.10,A2))',J1,SQRT(-EVALUES(J1)*9.3757D29)/(2*3.141592654D0),' i',
     &                                      SQRT(-EVALUES(J1)*9.3757D29)/(2*3.141592654D0*2.998D10),' i'
                  ENDIF
               ENDDO
            ELSEIF (CASTEP) THEN
               PRINT '(A,I6,A)',' geopt> ',NENDHESS,
     &          ' normal mode frequencies in Hz and wavenumbers, assuming eV and Angstrom units for input:'
               DO J1=NSTARTHESS,NENDHESS
                  IF (EVALUES(J1).GT.0.0D0) THEN
                     PRINT '(I6,2G20.10)',J1,SQRT(EVALUES(J1)*9.75586D27)/(2*3.141592654D0),
     &                                      SQRT(EVALUES(J1)*9.75586D27)/(2*3.141592654D0*2.998D10)
                  ELSE
                     PRINT '(I6,2(G20.10,A2))',J1,SQRT(-EVALUES(J1)*9.75586D27)/(2*3.141592654D0),' i',
     &                                   SQRT(-EVALUES(J1)*9.75586D27)/(2*3.141592654D0*2.998D10),' i'
                  ENDIF
               ENDDO
            ELSEIF (CHRMMT) THEN
               PRINT '(A,I6,A)',' geopt> ',3*NATOMS,' normal mode frequencies in Hz and wavenumbers'
               DO J1=NSTARTHESS,3*NATOMS
                  IF (EVALUES(J1).GT.0.0D0) THEN
                     PRINT '(I6,2G20.10E3)',J1,SQRT(EVALUES(J1)*4.184D26)/(2*3.141592654D0),
     &                 SQRT(EVALUES(J1)*4.184D26)/(2*3.141592654D0*2.998D10)
                  ELSE
                     PRINT '(I6,2(G20.10,A2))',J1,SQRT(-EVALUES(J1)*4.184D26)/(2*3.141592654D0),' i',
     &                 SQRT(-EVALUES(J1)*4.184D26)/(2*3.141592654D0*2.998D10),' i'
                  ENDIF
               ENDDO
            ELSEIF (TTM3T.OR.SDT) THEN
! hk286
               PRINT '(A,I6,A)',' geopt> ',NENDHESS-NSTARTHESS + 1, ' normal mode frequencies in Hz and wavenumbers'
               DO J1=NSTARTHESS,NENDHESS
                  IF (EVALUES(J1).GT.0.0D0) THEN
                     PRINT '(I6,2G20.10E3)',J1,SQRT(EVALUES(J1)*4.184D26)/(2*3.141592654D0),
     &                       SQRT(EVALUES(J1)*4.184D26)/(2*3.141592654D0*2.998D10)
                  ELSE
                     PRINT '(I6,2(G20.10,A2))',J1,SQRT(-EVALUES(J1)*4.184D26)/(2*3.141592654D0),' i',
     &                       SQRT(-EVALUES(J1)*4.184D26)/(2*3.141592654D0*2.998D10),' i'
                  ENDIF
               ENDDO
            ELSEIF (BOWMANT) THEN
               ! sn402: not sure whether the Bowman potential is supposed to be compatible with genrigid
               ! I've changed the array indices so that this bit of code is compatible.
               PRINT '(A,I6,A)',' geopt> ',NENDHESS-NSTARTHESS+1, ' normal mode frequencies in Hz and wavenumbers'
               DO J1=NSTARTHESS,NENDHESS
                  IF (EVALUES(J1).GT.0.0D0) THEN
                     PRINT '(I6,2G20.10)',J1,SQRT(EVALUES(J1)*2625.47D26)/(2*3.141592654D0),
     &                    SQRT(EVALUES(J1)*2625.47D26)/(2*3.141592654D0*2.998D10)
                  ELSE
                     PRINT '(I6,2(G20.10,A2))',J1,SQRT(-EVALUES(J1)*2625.47D26)/(2*3.141592654D0),' i',
     &                    SQRT(-EVALUES(J1)*2625.47D26)/(2*3.141592654D0*2.998D10),' i'
                  ENDIF
               ENDDO
            ELSE
               PRINT '(A,I6,A)',' geopt> ',(NENDHESS-NSTARTHESS+1),' Hessian eigenvalues:'
               PRINT '(6G20.10)',EVALUES(NSTARTHESS:NENDHESS)
               PRINT '(A,I6,A)',' geopt> Frequencies in specified units:'
               PRINT '(6G20.10)',DSQRT(ABS(EVALUES(NSTARTHESS:NENDHESS)))*FRQCONV
            ENDIF
         ENDIF  ! DEBUG block


         ! ZT and LZT determines which normal mode frequencies will be printed out. Generally, this means identifying the
         ! zero eigenvalues, which will be at the front (or perhaps the end) of the eigenvalue vector.
         ! In this case, all the normal modes are set to print - but this conflicts with the documentation which suggests
         ! that only non-zero eigenvalues should be printed.
         LZT(1:NOPT)=.TRUE.
         ! sn402: Adding this bit to make the behaviour match that described in the documentation.
         DO J1=1,NOPT
            IF (ABS(EVALUES(J1)).LT.EVCUT) LZT(J1)=.FALSE. ! EVCUT determines the minimum magnitude for a non-zero eigenvalue
         ENDDO
         ! If we've sorted eigenvalues in descending order, the first set will be dummies for RIGIDINIT.
         IF(RIGIDINIT) LZT(:NENDHESS-DEGFREEDOMS) = .FALSE.

! sn402: This is a problem, because it has the effect of undoing the sorting from the earlier subroutine call. Since we've already
! got the data we would obtain from this call, I'm going to comment it out. This has required editing some earlier code, but that's
! better than doing another hessian diagonalisation unnecessarily.
! hk286 - compute normal modes if rigid body angle-axis framework is used, in preparation for VDUMP subroutine
!         IF (RBAAT) THEN
!            RBAANORMALMODET = .TRUE.
!            CALL NRMLMD (Q, EVALUES, .TRUE.)
!            RBAANORMALMODET = .FALSE.
!         ENDIF
! hk286
         IF (DUMPV) THEN
           IF (RIGIDINIT) THEN
              ! sn402: should this be LZT here? Looks to me as though ZT is undefined.
              ! GENRIGID_EIGENVALUES is set up so that all dummy entries in DIAG (from rigidified degrees of freedom) are at
              ! the end of the vector,so will never be printed.
! sn402: changing the next line on a trial basis
!               CALL GENRIGID_VDUMP(DIAG(1:DEGFREEDOMS),ZT(1:DEGFREEDOMS),NOPT,DEGFREEDOMS)
               CALL GENRIGID_VDUMP(EVALUES,LZT,NOPT,NOPT)
            ELSE
               CALL VDUMP(EVALUES,LZT,NOPT,3*NATOMS)
            ENDIF
         ENDIF
         PROD=0.0D0

! Now we calculate the log product of positive frequencies, which will be written to min.data.info. This must be done regardless of whether we are dumping eigenvectors.
! Unit conversions are added to the product in a separate block - see below. At the moment, the eigenvalues should be in internal units.
! hk286
         ! MINFRQ2 is the (log of the) smallest positive hessian eigenvalue (i.e. smallest square frequency)
         IF (NENDHESS-NEXMODES.GT.0) THEN
            MINFRQ2=LOG(EVALUES(NENDHESS-NEXMODES))
         ELSE
            ! Dummy value, if there are no positive eigenvalues.
            MINFRQ2=1.0D0
         ENDIF
         EWARN=.FALSE.

         DO I1=NSTARTHESS,NENDHESS-NEXMODES
            IF (I1.GT.NSTARTHESS) THEN
               IF (EVALUES(I1-1).NE.0.0D0) THEN
                  IF (ABS(EVALUES(I1)/EVALUES(I1-1)).LT.1.0D-3) THEN
                     PRINT '(A,G20.10,A,G20.10)',' geopt> WARNING - decrease in magnitude of eigenvalues from ',EVALUES(I1-1),
     &                    ' to ',EVALUES(I1)
                     PRINT '(A)',' geopt> WARNING - this could indicate a stationary point of the wrong index'
                     EWARN=.TRUE.
                  ENDIF
               ENDIF
            ENDIF
            ! Calculate the actual log product here!
            IF (EVALUES(I1).GT.0.0D0) THEN
               PROD=PROD+DLOG(EVALUES(I1))
            ELSE
               IF (I1.LT.(NENDHESS-NEXMODES)) PRINT *,'Higher order saddle detected: eigenvalue ',EVALUES(I1)
               ! jmc put in this test mainly for pathsample purposes...
            ENDIF
         ENDDO

! Applying a unit conversion to PROD for the potentials which need it. Remember that PROD is a sum of logs, so multiplying
! the frequencies by a conversion factor is the same as adding the conversion factor to the sum (which is what we do here).
! NB: FRQCONV is defined as the conversion factor to go from internal units to SI (or another unit system) for the frequencies,
! not the square frequencies. However, PROD contains a log of squared frequencies. So we need to square FRQCONV (i.e. double its log)
         IF (FRQCONV.NE.1.0D0) THEN

            MINFRQ2=MINFRQ2+LOG(FRQCONV)*2
! hk286
            IF (RIGIDINIT) THEN
               PROD=PROD+(DEGFREEDOMS-NEXMODES)*DLOG(FRQCONV)*2
               WRITE (*,'(A,G20.10)') ' geopt> ln product scaled to desired frequency units by ',
     &                                (DEGFREEDOMS-NEXMODES)*DLOG(FRQCONV)*2
            ELSE
               PROD=PROD+(3*NATOMS-NEXMODES)*DLOG(FRQCONV)*2
               WRITE (*,'(A,G20.10)') ' geopt> ln product scaled to desired frequency units by ',(3*NATOMS-NEXMODES)*DLOG(FRQCONV)*2
            ENDIF

         ENDIF

! hk286
         IF (RIGIDINIT) THEN
            IF (NENDHESS-NEXMODES.GT.0) WRITE(*,'(A,I8,A,F20.10)') ' geopt> Log product of ',NENDHESS-NEXMODES-3*NATOMS+DEGFREEDOMS,
     &           ' positive Hessian eigenvalues=',PROD
         ELSE
            IF (NENDHESS-NEXMODES.GT.0) WRITE(*,'(A,I8,A,F20.10)') ' geopt> Log product of ',NENDHESS-NEXMODES,
     &           ' positive Hessian eigenvalues=',PROD
         ENDIF

         !kr366> dumps frqs into dump.frqs, used with GETMINFRQS keyword from PATHSAMPLE
         IF (DUMPFRQST) THEN
             IF (FILTH.EQ.0) THEN
                OPEN(91220,FILE='frqs.dump',POSITION='APPEND',ACTION='WRITE',STATUS='UNKNOWN')
             ELSE
                OPEN(91220,FILE='frqs.dump.'//TRIM(ADJUSTL(FILTHSTR)),POSITION='APPEND',ACTION='WRITE',STATUS='UNKNOWN')
             END IF
             WRITE(91220,*)  PROD
             CLOSE(91220)
         ENDIF

      ENDIF  ! End of IF(ENDHESS .AND. ...) (This block contains most of the subroutine, at the time of writing)

      IF (CHRMMT.AND.CALCDIHE) THEN
         STOP 'Necessary CHARMM routines not implemented yet for NSEG>1'
C         LSELECT=.FALSE.
C         CALL CHCALCRGYR(RGYR,Q,LSELECT)
C         LNATIVE=.FALSE.
C         CALL CHCALCNUMHB(NUMHB,Q,LNATIVE)
C         CALL CHCALCRMSD(RMSD,Q)
C         WRITE(*,'(A,4X,F20.10)') 'Final Rmsd from reference structure=',RMSD
C         WRITE(*,'(A,4X,F20.10)') 'Final Radius of gyration=',RGYR
C         IF (LNATIVE) THEN
C            WRITE(*,'(A,4X,I6)') 'Final Number of native hydrogen bonds=',NUMHB
C         ELSE
C            WRITE(*,'(A,4X,I6)') 'Final Number of cross-chain hydrogen bonds=',NUMHB
C         ENDIF
      ELSEIF (UNRST.AND.CALCDIHE) THEN
         CALL UNRESCALCDIHEREF(DIHE,ALLANG,Q)
         CALL UNRESCALCRGYR(RGYR,Q)
         WRITE(*,'(A,4X,F20.10)') ' Dihedral angle order parameter=',DIHE
         WRITE(*,'(A,4X,F20.10)') ' All angle order parameter=',ALLANG
         WRITE(*,'(A,4X,F20.10)') ' Radius of gyration=',RGYR
      ENDIF

! jmc49 Unconditionally executing the lines below, to make sure the Cartesian coords 
! are written correctly to points.final etc.
! hk286
      IF (RIGIDINIT .AND. (ATOMRIGIDCOORDT .EQV. .FALSE.)) THEN
         ! DUMQ will hold the cartesian coordinates. We don't change ATOMRIGIDCOORDT because we're only changing
         ! coordinates of a dummy array, not Q itself.
         CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, DUMQ, Q(1:DEGFREEDOMS))
      ELSE IF (GTHOMSONT) THEN
!jwrm2> With GTHOMSON, want Cartesian coordinates for moment of inertia calculation
         CALL GTHOMSONANGTOC(DUMQ(1:3*NATOMS), Q(1:3*NATOMS), NATOMS)
      ELSE
         DO II1=1,NOPT
            DUMQ(II1)=Q(II1)
         ENDDO
      ENDIF

      INERTIAT=.TRUE.
C     IF ((UNRST.OR.CHRMMT).AND.INERTIAT) THEN

      ! The next block is used for dumping to a min.data.info file. The frequencies calculated
      ! previously will be used, unless NOFRQS is set.
      IF ((DUMPDATAT).AND.INERTIAT) THEN
C
C  Pathsample uses unit masses, so we need to change to unit mass temporarily here,
C  and then put the right values back, otherwise the frequencies will be wrong
C  when we call geopt again from newconnect in a bhinterp run!!!!
C
         DO J1=1,NATOMS
            ATMASSSAVE(J1)=ATMASS(J1)
            ATMASS(J1)=1.0D0 
         ENDDO

         ! Get the inertia tensor
         IF (VARIABLES) THEN
            ITX=1.0D0; ITY=1.0D0; ITZ=1.0D0
         ELSE
            CALL INERTIA2(DUMQ,ITX,ITY,ITZ)
         ENDIF
         WRITE(*,'(A,4X,F20.10)') ' geopt> X component of inertia tensor=',ITX
         WRITE(*,'(A,4X,F20.10)') ' geopt> Y component of inertia tensor=',ITY
         WRITE(*,'(A,4X,F20.10)') ' geopt> Z component of inertia tensor=',ITZ
         HORDER=1
         DO J1=1,NATOMS
            ATMASS(J1)=ATMASSSAVE(J1)
         ENDDO
      ENDIF

! jwrm2> With GTHOMSON, coordinates are expected to be polar for symmetry
      IF (GTHOMSONT) THEN
         DO I1 = 1, 3*NATOMS
            DUMQ(I1) = Q(I1)
         END DO
      END IF

      IF (CHRMMT) THEN
         CALL CHARMMDUMP(DUMQ,FNAMEF,MACHINE)
         IF(INTMINT) THEN
           CALL BONDLENGTHDUMP(DUMQ)
         ENDIF
      ELSE IF (UNRST) THEN
C
C jmc myunresdump just outputs Calpha and side chain coords in plain xyz format, so they can be
C read in again in a filthy run, whereas unresdump3 adds dummy O and N atoms for visualization purposes.
C
         CALL MYUNRESDUMP(DUMQ,FNAMEF)
         CALL UNRESDUMP3(DUMQ,'unr.'//TRIM(ADJUSTL(FNAMEF)))
      ELSE
!
! Call symmetry before dumping the file, otherwise the wrong point group is reported in odata.new and picked up by Filthy_Phyllis
!
         IF (.NOT.(VARIABLES.OR.RINGPOLYMERT).AND.(.NOT.(ZSYM(NATOMS).EQ.'TH'))) CALL SYMMETRY(HORDER,.TRUE.,DUMQ,IT)
         if (MACHINE) then
             call WriteOutFile(DUMQ,FNAMEF)
         else
             CALL DUMPIT(DUMQ,FNAMEF)
         endif
      ENDIF
11    CONTINUE
      IF ((BHINTERPT.OR.BISECTT).AND.(.NOT.NEWCONNECTT).AND.(.NOT.REOPTIMISEENDPOINTS)) THEN ! do nothing
      ! sn402: I think MFLAG is used when we have converged to a minimum and want to write it to min.data.info.
      ELSE IF (MFLAG) THEN
        ! csw34> Mass weighted normal mode dumping 
         IF (ENDHESS .AND. (.NOT. (REOPTIMISEENDPOINTS.AND.(BHINTERPT.OR.BISECTT)))) THEN
! in this case we have already called vdump. We don;t want to do it again!
! This IF statement (ENDHESS AND...) is equivalent to the one near the start of the subroutine.
! So quite a lot of the code here duplicates the earlier behaviour.
         ELSEIF (DUMPV.AND.ALLVECTORS) THEN
                DO J1=1,3*NATOMS
                   ZT(J1)=.TRUE.
                ENDDO
! Set the first six normal modes to not be printed. For a non-linear
! molecule, these correspond to the pure rotations and translations
                DO J1=1,MIN(6,NOPT)
                   ZT(J1)=.FALSE.
                ENDDO

! Here we compute the frequencies/eigenvalues themselves.                            
! Can't mass weight an already diagonalised Hessian so call potential again
! hk286 - tells potential that we switch to moving frame - needed for normal modes calculations
                RBAANORMALMODET = .TRUE.
                CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
                RBAANORMALMODET = .FALSE.
! hk286 - do better than putting this here?
                IF (RIGIDINIT) THEN
                   CALL GENRIGID_EIGENVALUES(Q, ATMASS, DIAG, INFO)
                ELSE
! Mass weight hessian - using MASSWT would also change the coordinates
                   IF (MWVECTORS) THEN
                      IF (RINGPOLYMERT) THEN
                         CALL MASSWTRP(NOPT,RPMASSES,RPDOF)
                      ELSE
                         ! sn402: we're wasting a bit of time here when RBAAT is TRUE - we call this function
                         ! but we then overwrite the mass-weighted hessian in NRMLMD
                         CALL MASSWT2(NATOMS,ATMASS,Q,VNEW,.TRUE.)
                      ENDIF
                   ENDIF
! Diagonalise
! hk286 - normal modes calculation, for rigid body angle-axis the diagonalization is rather complex
                   IF (RBAAT.AND..NOT.PYGPERIODICT.AND..NOT.SANDBOXT) THEN
                      RBAANORMALMODET = .TRUE.
                      CALL NRMLMD (Q, DIAG, .TRUE.)
                      RBAANORMALMODET = .FALSE.
                   ELSE
                      CALL DSYEV('V','U',NOPT,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
                   ENDIF
! hk286
! If we're freezing atoms, some zero eiganvalue modes can creep in that
! are NOT real i.e. one atoms moving and all others stationary. Here, we
! check and remove these using ZT
                   IF (FREEZE) THEN
                      DO J1=1,3*NATOMS
                         IF (ABS(DIAG(J1)).LT.0.000001) ZT(J1)=.FALSE.
                      ENDDO
                   ENDIF
                ENDIF
! Dump vectors - we will only get here if ENDHESS was set, and only get to the earlier block if it wasn't.
!hk286
                IF (RIGIDINIT) THEN
                   CALL GENRIGID_VDUMP(DIAG(1:DEGFREEDOMS),ZT(1:DEGFREEDOMS),NOPT,DEGFREEDOMS)
                ELSE
                   CALL VDUMP(DIAG,ZT,NOPT,3*NATOMS)
                ENDIF
! hk286 - do I need to change the following for CHARMM and AMBER?

! If using CHARMM, call the CHARMMDUMPMODES subroutine to output a PDB
! containing one frame per mode and scaled atomic and residue
! displacements - coords are in Q here!
                IF (CHRMMT.AND.MWVECTORS) CALL CHARMMDUMPMODES(Q,DIAG,ZT,NOPT,3*NATOMS)
! AMBER call to routine in amberinterface.f
                IF ((AMBERT.OR.NABT).AND.MWVECTORS) THEN
                        CALL A9DUMPMODES(DIAG,ZT,NOPT,3*NATOMS)
                ENDIF
              IF(NORMALMODET) CALL VISNORMODES(Q)
         ENDIF ! End of IF(ENDHESS AND...)

         GOODSTRUC1 = .TRUE.
         GOODSTRUC2 = .TRUE.
         IF (AMBER12T) THEN
            IF (RIGIDINIT .AND. (ATOMRIGIDCOORDT .EQV. .FALSE.)) THEN
               CALL TRANSFORMRIGIDTOC(1,NRIGIDBODY,DUMQ,Q(1:DEGFREEDOMS))
               CALL CIS_TRANS_CHECK(DUMQ, GOODSTRUC1)
               CALL CHIRALITY_CHECK(DUMQ, GOODSTRUC2)
            ELSE
               CALL CIS_TRANS_CHECK(Q, GOODSTRUC1)
               CALL CHIRALITY_CHECK(Q, GOODSTRUC2)
            ENDIF
         END IF
         IF (.NOT. GOODSTRUC1) THEN
            WRITE(*,'(A)') ' geopt> cistrans problem' 
         END IF
         IF (.NOT. GOODSTRUC2) THEN
            WRITE(*,'(A)') ' geopt> chirality problem'
         END IF

         PRINT*
         WRITE(*,'(A)') ' geopt>                          **** CONVERGED ****'
         PRINT*
         CALL FLUSH(6)

!
! Print TIP frequencies in cm-1 for CoM/Euler angle representation for debugging angle-axis.
!
         IF (DEBUG.AND.(ZSYM(1)(1:1).EQ.'W')) THEN
            IF (ZSYM(1)(1:2).EQ.'W4') IPOT=4
            IF (ZSYM(1)(1:2).EQ.'W3') IPOT=3
            IF (ZSYM(1)(1:2).EQ.'W2') IPOT=2
            IF (ZSYM(1)(1:2).EQ.'W1') IPOT=1
            CALL H2OMODES(NATOMS/2,IPOT,Q,DIAG)
            PRINT '(A,I6,A)',' geopt> TIP normal mode frequencies in wavenumbers'
            DO J1=1,3*NATOMS
               IF (DIAG(J1).GT.0.0D0) THEN
                  PRINT '(I6,2G20.10)',J1,DIAG(J1)
               ELSE
                  PRINT '(I6,2(G20.10,A2))',J1,-DIAG(J1),' i'
               ENDIF
            ENDDO
         ENDIF
      ELSE  ! End of IF(MFLAG)
         ! If we haven't converged to a minimum, print a bit and exit.
         IF (GRADSQ) WRITE(*,'(A,4F20.10)') ' g^2, RMS force and real energy and RMS=',ENERGY,RMS2,EREAL,RMS
         PRINT*
         IF (BFGSMINT.OR.(BFGSTST.AND.(HINDEX.EQ.0))) THEN
              PRINT*,ITDONE,' steps completed without convergence to required tolerance'
         ELSE
              PRINT*,NSTEPS,' steps completed without convergence to required tolerance'
         ENDIF
         PRINT*
         CALL FLUSH(6)
!        STOP
         IF (.NOT.MULTIJOBT) CALL EXIT(0)
      ENDIF
      CALL FLUSH(6)

! If we are making a min.data.info file, do it here.
      IF (DUMPDATAT) THEN 
!        IF (EWARN.AND.MULTIJOBT) MFLAG=.FALSE.
         IF (MFLAG) THEN
!
! min.data.info file is now opened on unit 881 in keyword.f
!
!           OPEN(UNIT=100,FILE='min.data.info',STATUS='UNKNOWN')
            ! sn402: need Cartesian coords to calculate the inertia tensor
            IF (RIGIDINIT .AND. (ATOMRIGIDCOORDT .EQV. .FALSE.)) THEN
              CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, DUMQ, Q(1:DEGFREEDOMS))
            ENDIF
            IF (BHINTERPT.AND.(.NOT.REOPTIMISEENDPOINTS)) ENERGY=BHENERGY
            IF (BISECTT.AND.(.NOT.REOPTIMISEENDPOINTS)) ENERGY=BISECTENERGY
            IF (LOWESTFRQT) THEN
               IF (CHRMMT.OR.UNRST.OR.AMBERT.OR.NABT.OR.AMBER12T) THEN
                  ! NOTE: PROD is a log product of square frequencies in the square of the units defined by FRQCONV.
                  WRITE(881,'(2F20.10,I6,5F20.10)') ENERGY,PROD,HORDER,ITX,ITY,ITZ,MINCURVE,MINFRQ2
               ELSE
                  CALL INERTIA2(Q,ITX,ITY,ITZ)
                  ! NOTE: PROD is a log product of square frequencies in the square of the units defined by FRQCONV.
                  WRITE(881,'(2F20.10,I6,5F20.10)') ENERGY,PROD,HORDER,ITX,ITY,ITZ,MINCURVE,MINFRQ2
               ENDIF
            ELSE   
               IF (CHRMMT.OR.UNRST.OR.AMBERT.OR.NABT.OR.AMBER12T) THEN
                  ! NOTE: PROD is a log product of square frequencies in the square of the unit defined by FRQCONV.
                  ! The default is SI units (s^-2) in this case.
                  WRITE(881,'(2F20.10,I6,4F20.10)') ENERGY,PROD,HORDER,ITX,ITY,ITZ
               ELSE
                  IF(RIGIDINIT .AND. (ATOMRIGIDCOORDT .EQV. .FALSE.)) THEN
                      CALL INERTIA2(DUMQ,ITX,ITY,ITZ)
                  ELSE IF (GTHOMSONT) THEN
                     CALL GTHOMSONANGTOC(TMPC(1:3*NATOMS), Q(1:3*NATOMS), NATOMS)
                     CALL INERTIA2(TMPC,ITX,ITY,ITZ)
                  ENDIF

                  ! Write the header line to min.data.info. PROD should be 1 if NOFRQS is set, or a positive real value
                  ! if ENDHESS is set. If neither of these is set, PROD is probably undefined (for some reason we don't
                  ! calculate it in the MFLAG section - I'll look into adding it) -sn402
                  ! If we do have ENDHESS and not NOFRQS, PROD is a log product of square frequencies in the square of
                  ! the unit defined by FRQCONV.
                  WRITE(881,'(2F20.10,I6,4F20.10)') ENERGY,PROD,HORDER,ITX,ITY,ITZ
               ENDIF
            ENDIF
            NRES=NMRES

            IF (AMHT) THEN
               GLY_COUNT = 0
               SDUMMY='AM'
               DO J2=1,NRES
                 IF (SEQ(J2).EQ.8) THEN
76                format(A2,3(F20.10))
                     WRITE(881,*)Q(9*(J2-1)+1-GLY_COUNT*3),Q(9*(J2-1)+2-GLY_COUNT*3),Q(9*(J2-1)+3-GLY_COUNT*3)
                     WRITE(881,*)Q(9*(J2-1)+1-GLY_COUNT*3),Q(9*(J2-1)+2-GLY_COUNT*3),Q(9*(J2-1)+3-GLY_COUNT*3)
                     WRITE(881,*)Q(9*(J2-1)+4-GLY_COUNT*3),Q(9*(J2-1)+5-GLY_COUNT*3),Q(9*(J2-1)+6-GLY_COUNT*3)
                   GLY_COUNT = GLY_COUNT +1
                 ELSE
                     WRITE(881,*)Q(9*(J2-1)+1-GLY_COUNT*3),Q(9*(J2-1)+2-GLY_COUNT*3),Q(9*(J2-1)+3-GLY_COUNT*3)   
                     WRITE(881,*)Q(9*(J2-1)+4-GLY_COUNT*3),Q(9*(J2-1)+5-GLY_COUNT*3),Q(9*(J2-1)+6-GLY_COUNT*3)
                     WRITE(881,*)Q(9*(J2-1)+7-GLY_COUNT*3),Q(9*(J2-1)+8-GLY_COUNT*3),Q(9*(J2-1)+9-GLY_COUNT*3)
                 ENDIF

               ENDDO
               CALL FLUSH(881)
            ELSE

               ! Now we write the coordinates, which is the main part of the min.data.info file. 

               IF (DUMPDATA_MACHINET) THEN
                   RECCOUNT = RECCOUNT + 1
               ENDIF

! hk286
               IF (RIGIDINIT) THEN
                  CALL TRANSFORMRIGIDTOC (1, NRIGIDBODY, XCOORDS, Q(1:DEGFREEDOMS))
                  IF (DUMPDATA_MACHINET) THEN
!                      WRITE(882, REC=RECCOUNT) XCOORDS(1:3*NATOMS)
                      WRITE(882, REC=RECCOUNT) XCOORDS(1:3*NATOMS)
                  ELSE
                      WRITE(881,'(3F25.15)') XCOORDS(1:3*NATOMS)
                  ENDIF

               ELSE IF (GTHOMSONT) THEN
                  CALL GTHOMSONANGTOC(TMPC(1:3*NATOMS), Q(1:3*NATOMS), NATOMS)
                  IF (DUMPDATA_MACHINET) THEN
                      WRITE(882, REC=RECCOUNT) TMPC(1:3*NATOMS)
                  ELSE
                      WRITE(881,'(3F25.15)') TMPC(1:3*NATOMS)
                  ENDIF

               ELSE IF (VARIABLES) THEN
                  IF (DUMPDATA_MACHINET) THEN
                      WRITE(882, REC=RECCOUNT) Q(1:NATOMS)
                  ELSE
                      WRITE(881,'(F25.15)') Q(1:NATOMS)
                  ENDIF

               ELSE
                  IF (DUMPDATA_MACHINET) THEN
                      WRITE(882, REC=RECCOUNT) Q(1:3*NATOMS)
                  ELSE
                      WRITE(881,'(3F25.15)') Q(1:3*NATOMS)
                  ENDIF
               ENDIF
               CALL FLUSH(881)
               IF (DUMPDATA_MACHINET) CALL FLUSH(882)
            ENDIF
         ELSE
            PRINT '(A)',' geopt> WARNING - DUMPDATA is set, but MFLAG is false; not creating entry in min.data.info'
         ENDIF
      ENDIF

      RETURN

      END  ! END GEOPT!

      SUBROUTINE MASSWT(NATOMS,ATMASS,Q,VNEW,STEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER NATOMS, J1, J2, J3, J4
      LOGICAL STEST
      DOUBLE PRECISION ATMASS(NATOMS),Q(3*NATOMS),VNEW(3*NATOMS),AMASS,BMASS,PMASS

      IF (SIZE(HESS,2).LT.3*NATOMS) RETURN ! for example, for 2D XY model

      IF (.NOT.ALLOCATED(HESS)) THEN
         PRINT '(A)','masswt2> ERROR - HESS has not been allocated. Do you need NOFRQS, ENDHESS, or ENDNUMHESS in odata?'
         STOP
      ENDIF
      DO J1=1,NATOMS
!        PRINT *,'J1,ATMASS=',J1,ATMASS(J1)
         AMASS=1.0D0/SQRT(ATMASS(J1))
         BMASS=SQRT(ATMASS(J1))
         J3=3*J1
         Q(J3-2)=BMASS*Q(J3-2)
         Q(J3-1)=BMASS*Q(J3-1)
         Q(J3)=BMASS*Q(J3)
         VNEW(J3-2)=VNEW(J3-2)*AMASS
         VNEW(J3-1)=VNEW(J3-1)*AMASS
         VNEW(J3)=VNEW(J3)*AMASS
         IF (STEST) THEN
            DO J2=J1,NATOMS
               BMASS=1.0D0/SQRT(ATMASS(J2))
               PMASS=AMASS*BMASS
               J4=3*J2
               IF (J1.EQ.J2) THEN
                  HESS(J3-2,J4-2)=PMASS*HESS(J3-2,J4-2)
                  HESS(J3-2,J4-1)=PMASS*HESS(J3-2,J4-1)
                  HESS(J3-2,J4)  =PMASS*HESS(J3-2,J4)
                  HESS(J3-1,J4-2)=PMASS*HESS(J3-1,J4-2)
                  HESS(J3-1,J4-1)=PMASS*HESS(J3-1,J4-1)
                  HESS(J3-1,J4)  =PMASS*HESS(J3-1,J4)
                  HESS(J3,  J4-2)=PMASS*HESS(J3,  J4-2)
                  HESS(J3,  J4-1)=PMASS*HESS(J3,  J4-1)
                  HESS(J3,  J4)  =PMASS*HESS(J3,  J4)
               ELSE
                  HESS(J3-2,J4-2)=PMASS*HESS(J3-2,J4-2)
                  HESS(J4-2,J3-2)=HESS(J3-2,J4-2)
                  HESS(J3-2,J4-1)=PMASS*HESS(J3-2,J4-1)
                  HESS(J4-1,J3-2)=HESS(J3-2,J4-1)
                  HESS(J3-2,J4)=PMASS*HESS(J3-2,J4)
                  HESS(J4,J3-2)=HESS(J3-2,J4)
                  HESS(J3-1,J4-2)=PMASS*HESS(J3-1,J4-2)
                  HESS(J4-2,J3-1)=HESS(J3-1,J4-2)
                  HESS(J3-1,J4-1)=PMASS*HESS(J3-1,J4-1)
                  HESS(J4-1,J3-1)=HESS(J3-1,J4-1)
                  HESS(J3-1,J4)=PMASS*HESS(J3-1,J4)
                  HESS(J4,J3-1)=HESS(J3-1,J4)
                  HESS(J3,J4-2)=PMASS*HESS(J3,J4-2)
                  HESS(J4-2,J3)=HESS(J3,  J4-2)
                  HESS(J3,J4-1)=PMASS*HESS(J3,J4-1)
                  HESS(J4-1,J3)=HESS(J3,  J4-1)
                  HESS(J3,J4)=PMASS*HESS(J3,J4)
                  HESS(J4,J3)=HESS(J3,J4)
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      RETURN
      END

      SUBROUTINE MASSWT2(NATOMS,ATMASS,Q,VNEW,STEST)
      USE MODHESS
      IMPLICIT NONE
      INTEGER NATOMS, J1, J2, J3, J4
      LOGICAL STEST
      DOUBLE PRECISION ATMASS(NATOMS),Q(3*NATOMS),VNEW(3*NATOMS),AMASS,BMASS,PMASS

      IF (SIZE(HESS,2).LT.3*NATOMS) RETURN ! for example, for 2D XY model or VARIABLES
      IF (.NOT.ALLOCATED(HESS)) THEN
         PRINT '(A)','masswt2> ERROR - HESS has not been allocated. Do you need NOFRQS, ENDHESS, or ENDNUMHESS in odata?'
         STOP
      ENDIF
      DO J1=1,NATOMS
         AMASS=1/SQRT(ATMASS(J1))
         BMASS=SQRT(ATMASS(J1))
         J3=3*J1
C DAE don't mass weight the coordinates or gradients
         IF (STEST) THEN
            DO J2=J1,NATOMS
               BMASS=1/SQRT(ATMASS(J2))
               PMASS=AMASS*BMASS
               J4=3*J2
               IF (J1.EQ.J2) THEN
                  HESS(J3-2,J4-2)=PMASS*HESS(J3-2,J4-2)
                  HESS(J3-2,J4-1)=PMASS*HESS(J3-2,J4-1)
                  HESS(J3-2,J4)  =PMASS*HESS(J3-2,J4)
                  HESS(J3-1,J4-2)=PMASS*HESS(J3-1,J4-2)
                  HESS(J3-1,J4-1)=PMASS*HESS(J3-1,J4-1)
                  HESS(J3-1,J4)  =PMASS*HESS(J3-1,J4)
                  HESS(J3,  J4-2)=PMASS*HESS(J3,  J4-2)
                  HESS(J3,  J4-1)=PMASS*HESS(J3,  J4-1)
                  HESS(J3,  J4)  =PMASS*HESS(J3,  J4)
               ELSE
                  HESS(J3-2,J4-2)=PMASS*HESS(J3-2,J4-2)
                  HESS(J4-2,J3-2)=HESS(J3-2,J4-2)
                  HESS(J3-2,J4-1)=PMASS*HESS(J3-2,J4-1)
                  HESS(J4-1,J3-2)=HESS(J3-2,J4-1)
                  HESS(J3-2,J4)=PMASS*HESS(J3-2,J4)
                  HESS(J4,J3-2)=HESS(J3-2,J4)
                  HESS(J3-1,J4-2)=PMASS*HESS(J3-1,J4-2)
                  HESS(J4-2,J3-1)=HESS(J3-1,J4-2)
                  HESS(J3-1,J4-1)=PMASS*HESS(J3-1,J4-1)
                  HESS(J4-1,J3-1)=HESS(J3-1,J4-1)
                  HESS(J3-1,J4)=PMASS*HESS(J3-1,J4)
                  HESS(J4,J3-1)=HESS(J3-1,J4)
                  HESS(J3,J4-2)=PMASS*HESS(J3,J4-2)
                  HESS(J4-2,J3)=HESS(J3,  J4-2)
                  HESS(J3,J4-1)=PMASS*HESS(J3,J4-1)
                  HESS(J4-1,J3)=HESS(J3,  J4-1)
                  HESS(J3,J4)=PMASS*HESS(J3,J4)
                  HESS(J4,J3)=HESS(J3,J4)
               ENDIF
            ENDDO
         ENDIF
      ENDDO

      RETURN
      END

      SUBROUTINE MASSWTRP(NOPT,RPMASSES,RPDOF) ! mass-weight upper-triangular part of HESS
      USE MODHESS
      IMPLICIT NONE
      INTEGER J1, J2, RPDOF, NOPT, MINDEX
      DOUBLE PRECISION RPMASSES(RPDOF),AMASS,BMASS,PMASS

      IF (.NOT.ALLOCATED(HESS)) THEN
         PRINT '(A)','masswt2> ERROR - HESS has not been allocated. Do you need NOFRQS, ENDHESS, or ENDNUMHESS in odata?'
         STOP
      ENDIF
      DO J1=1,NOPT
         MINDEX=MOD(J1-1,RPDOF)+1 ! masses are the same for the RPDOF degrees of freedom in each bead
         AMASS=1.0D0/SQRT(RPMASSES(MINDEX))
         DO J2=1,J1
            MINDEX=MOD(J2-1,RPDOF)+1
            BMASS=1/SQRT(RPMASSES(MINDEX))
            PMASS=AMASS*BMASS
            HESS(J2,J1)=PMASS*HESS(J2,J1)
         ENDDO
      ENDDO

      RETURN
      END
!
! Somewhat better eigenvalues for twice the cpu time. Uses Richardson extrapolation.
!
      SUBROUTINE MAKENUMHESS2(X,NATOMS)
      USE MODHESS
      USE MODCHARMM
      USE KEY,ONLY : FROZEN, CASTEP, AMHT, ONETEP, ENDNUMHESSDELTA, VARIABLES
      USE COMMONS,ONLY : DEBUG, NOPT
      use porfuncs
      IMPLICIT NONE
      LOGICAL KNOWE, KNOWG, KNOWH
      COMMON /KNOWN/ KNOWE, KNOWG, KNOWH

      INTEGER I1,J1,NATOMS,ISTAT
      DOUBLE PRECISION X(NOPT)
      DOUBLE PRECISION DUM(NOPT),GRAD1(NOPT),GRAD2(NOPT),DELTA,RMS,ENERGY
      DOUBLE PRECISION GRAD1B(NOPT),GRAD2B(NOPT)

      IF (DEBUG) WRITE(*,'(A,I10)') ' makenumhess2> Making numerical hessian using Richardson extrapolation dimension=',NOPT
      DO I1=1,NOPT
         DUM(I1)=X(I1)
      ENDDO

      DELTA=ENDNUMHESSDELTA
      PRINT '(A,G20.10)',' makenumhess2> displacement=',DELTA
C
      IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1

      IF (.NOT.ALLOCATED(HESS)) ALLOCATE(HESS(NOPT,NOPT))
      DO I1=1,NOPT
         CALL FLUSH(6)
         IF (VARIABLES.AND.FROZEN(I1)) THEN
            DO J1=I1,NOPT
               HESS(I1,J1)=0.0D0
               HESS(J1,I1)=0.0D0
            ENDDO
         ELSE IF ((.NOT.VARIABLES).AND.FROZEN((I1-1)/3+1)) THEN
            DO J1=I1,NOPT
               HESS(I1,J1)=0.0D0
               HESS(J1,I1)=0.0D0
            ENDDO
         ELSE
            DUM(I1)=X(I1)-DELTA
            CALL POTENTIAL(DUM,ENERGY,GRAD1,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            DUM(I1)=X(I1)-2.0D0*DELTA
            CALL POTENTIAL(DUM,ENERGY,GRAD1B,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            DUM(I1)=X(I1)+DELTA
            CALL POTENTIAL(DUM,ENERGY,GRAD2,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            DUM(I1)=X(I1)+2.0D0*DELTA
            CALL POTENTIAL(DUM,ENERGY,GRAD2B,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            DUM(I1)=X(I1)
            DO J1=I1,NOPT
               HESS(I1,J1)=(4.0D0*(GRAD2(J1)-GRAD1(J1))-(GRAD2B(J1)-GRAD1B(J1))/2.0D0)/(2.0D0*DELTA*3.0D0)
               HESS(J1,I1)=HESS(I1,J1)
            ENDDO
         ENDIF
      ENDDO

      IF (DEBUG) WRITE(*,'(A)') ' makenumhess2> Hessian made'
      KNOWH=.TRUE.

      RETURN
      END

      SUBROUTINE MAKENUMHESS(X,NATOMS)
C
C dae
C
      USE MODHESS
      USE MODCHARMM
      USE KEY,ONLY : FROZEN, CASTEP, AMHT, ONETEP, ENDNUMHESS2, QCHEM, VASP, VARIABLES
      USE SBMDATA,ONLY : SBMHESSATOM,SBMFIRSTDD
      USE COMMONS,ONLY : DEBUG, NOPT
      use porfuncs
      IMPLICIT NONE
      LOGICAL KNOWE, KNOWG, KNOWH
      COMMON /KNOWN/ KNOWE, KNOWG, KNOWH

      INTEGER I1,J1,NATOMS,ISTAT
      DOUBLE PRECISION X(NOPT)
      DOUBLE PRECISION DUM(NOPT),GRAD1(NOPT),GRAD2(NOPT),DELTA,RMS,ENERGY

      IF (ENDNUMHESS2) THEN
         CALL MAKENUMHESS2(X,NATOMS)
         RETURN
      END IF

      IF (DEBUG) WRITE(*,'(A,I6)') ' makenumhess> Making numerical hessian dimension=',NOPT
      DO I1=1,NOPT
         DUM(I1)=X(I1)
      ENDDO

      DELTA=1.0D-6
      IF (CASTEP) DELTA=0.001D0
      IF (QCHEM) DELTA=0.001D0
      IF (VASP) DELTA=0.001D0
      IF (ONETEP) DELTA=0.0001D0
      IF (AMHT) DELTA=1.0D-6
C
      IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
 
C
C DAE having two potential calls for each hessian evaluation increases the accuracy
C compared to the analytical solution significantly (around 4sf), relative to the alternative
C of evaluating the gradient once at the beginning then once for each element. But obviously it also
C doubles the number of potential evaluations ((3N)^2/2 for N atoms), so may be not worth it for
C a long run on a big system.
C
      IF (.NOT.ALLOCATED(HESS)) ALLOCATE(HESS(NOPT,NOPT))
!      SBMFIRSTDD = .TRUE.
      DO I1=1,NOPT
!        write(*,*) 'working on atom', I1
         CALL FLUSH(6)
         IF (VARIABLES) THEN
            IF (FROZEN(I1)) THEN
               DO J1=I1,NOPT
                  HESS(I1,J1)=0.0D0
                  HESS(J1,I1)=0.0D0
               ENDDO
            ENDIF
         ELSE IF ((.NOT.VARIABLES).AND.FROZEN((I1-1)/3+1)) THEN
            DO J1=I1,NOPT
               HESS(I1,J1)=0.0D0
               HESS(J1,I1)=0.0D0
            ENDDO
         ELSE
            DUM(I1)=X(I1)-DELTA
!            SBMHESSATOM=int(I1/3)
            CALL POTENTIAL(DUM,ENERGY,GRAD1,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            DUM(I1)=X(I1)+DELTA
            CALL POTENTIAL(DUM,ENERGY,GRAD2,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            DUM(I1)=X(I1)
            DO J1=I1,NOPT
               HESS(I1,J1)=(GRAD2(J1)-GRAD1(J1))/(2.0D0*DELTA)
               HESS(J1,I1)=HESS(I1,J1)
            ENDDO
         ENDIF
      ENDDO
!      SBMHESSATOM=-1
      IF (DEBUG) WRITE(*,'(A)') ' makenumhess> Hessian made'
      KNOWH=.TRUE.
 
      RETURN
      END

      SUBROUTINE MAKENUMHESSRP(X,NOPT)
C
C dae
C
      USE MODHESS
      USE MODCHARMM
      USE KEY,ONLY : FROZEN, CASTEP, AMHT
      USE COMMONS,ONLY : DEBUG
      use porfuncs
      IMPLICIT NONE
      LOGICAL KNOWE, KNOWG, KNOWH
      COMMON /KNOWN/ KNOWE, KNOWG, KNOWH

      INTEGER I1,J1,NOPT,ISTAT
      DOUBLE PRECISION X(NOPT)
      DOUBLE PRECISION DUM(NOPT),GRAD1(NOPT),GRAD2(NOPT),DELTA,RMS,ENERGY

      IF (DEBUG) WRITE(*,'(A)') ' makenumhess> Making numerical hessian'
      DO I1=1,NOPT
         DUM(I1)=X(I1)
      ENDDO

      DELTA=1.0D-5
      IF (.NOT.ALLOCATED(HESS)) ALLOCATE(HESS(NOPT,NOPT))
      DO I1=1,NOPT
         CALL FLUSH(6)
         IF (FROZEN((I1-1)/3+1)) THEN
            DO J1=I1,NOPT
               HESS(I1,J1)=0.0D0
               HESS(J1,I1)=0.0D0
            ENDDO
         ELSE
            DUM(I1)=X(I1)-DELTA
            CALL POTENTIAL(DUM,ENERGY,GRAD1,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            DUM(I1)=X(I1)+DELTA
            CALL POTENTIAL(DUM,ENERGY,GRAD2,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            DUM(I1)=X(I1)
            DO J1=I1,NOPT
               HESS(I1,J1)=(GRAD2(J1)-GRAD1(J1))/(2.0D0*DELTA)
               HESS(J1,I1)=HESS(I1,J1)
            ENDDO
         ENDIF
      ENDDO

      IF (DEBUG) WRITE(*,'(A)') ' makenumhessrp> Hessian made'
      KNOWH=.TRUE.
 
      RETURN
      END

      SUBROUTINE MAKENUMINTHESS(NOPT,NATOMS)
      USE MODHESS
      USE COMMONS,ONLY : DEBUG
C
C jmc assuming that unres internal coordinates have already been updated.
C This subroutine will not affect the stored internal or Cartesian coordinates.
C
C dae
C
      IMPLICIT NONE

      INTEGER I1,J1,NATOMS,NOPT
      DOUBLE PRECISION X(NOPT)
      DOUBLE PRECISION DUM2(3*NATOMS),GRAD1(3*NATOMS),GRAD2(3*NATOMS),DELTA,RMS,ENERGY
      DOUBLE PRECISION DUM(NOPT)

      IF (DEBUG) WRITE(*,*) 'Making numerical hessian .=.'
      CALL geom_to_var(NOPT,X)

      DUM=X

      DUM2=1.0D0

      DELTA=1.0D-6
C
C dae having two potential calls for each hessian evaluation increases the accuracy
C compared to the analytical solution significantly (around 4sf), relative to the alternative
C of evaluating the gradient once at the beginning then once for each element. But obviously it also
C doubles the number of potential evaluations ((3N)^2/2 for N atoms), so may be not worth it for
C a long run on a big system.
C
      IF (.NOT.ALLOCATED(HESS)) ALLOCATE(HESS(3*NATOMS,3*NATOMS))
      DO I1=1,NOPT
         DUM(I1)=X(I1)-DELTA
         CALL var_to_geom(NOPT,DUM)
         CALL chainbuild
         CALL POTENTIAL(DUM2,ENERGY,GRAD1,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
         DUM(I1)=X(I1)+DELTA
         CALL var_to_geom(NOPT,DUM)
         CALL chainbuild
         CALL POTENTIAL(DUM2,ENERGY,GRAD2,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
         DUM(I1)=X(I1)
         DO J1=I1,NOPT
            HESS(I1,J1)=(GRAD2(J1)-GRAD1(J1))/(2.0D0*DELTA)
            HESS(J1,I1)=HESS(I1,J1)
         ENDDO
      ENDDO

      CALL var_to_geom(NOPT,X)
      CALL chainbuild

      IF (DEBUG) WRITE(*,*) 'Hessian made'

      RETURN
      END

C ********************************************************************************************
C ********************************************************************************************
C ********************************************************************************************

      subroutine BONDLENGTHDUMP(X)
         use intcutils, only : CART2INT
         use intcommons, only : NBDS
         use commons, only : NINTS, NATOMS

         implicit none

         double precision X(3*NATOMS)
         double precision XINT(NINTS)
         integer i1
         call CART2INT(X,XINT)        
         open (UNIT=973,FILE='bonds.final')
         print *, "Dumping",NBDS,"bonds"
         do i1=1,NBDS
            write (973,*) XINT(i1)
         enddo
         close(973)           
      end subroutine
