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
!   Eigenvector-following optimization algorithm.
!
MODULE MODEFOL
      IMPLICIT NONE
      CONTAINS
      SUBROUTINE EFOL(QTS,MFLAG,ITMAX,ENERGY,ITER,SMALL,PTEST,DIAG,INRIN)
      USE COMMONS
      USE KEY
      USE SYMINF
      USE MODHESS
      USE ZWK
      use porfuncs
      USE modcharmm, ONLY: CHRMMT
      IMPLICIT NONE

!     integer,intent(in),optional :: INRIN
      INTEGER,INTENT(IN) :: INRIN
      LOGICAL DONE, TSTEST, MINTEST, NRTEST, ZT(3*NATOMS),& 
     &        PZT(3*NATOMS), SDTEST, AWAY, MFLAG, TEST1, TEST2, PTEST
      INTEGER OMODE, K1, IMODE, I, J, INEG, IASSIGN, IM, ICOUNT, NATOMSSAVE, ISTAT,&
     &        J1, J2, INFO, ITMAX, ITER, FRAME, HORDER, INRSAVE, NEV
      DOUBLE PRECISION TEMPA(9*NATOMS), SSTPMAG, STPMAG, TPAR, &
     &   FOB(3*NATOMS), CSTEP(3*NATOMS), STEP(3*NATOMS), TEMP,&
     &   RMS, PROD, SUM, RAT1, RAT2, DELTASP, EOLD, EPER, ENERGY, AMASS, R1, R2, &
     &   VEC(3*NATOMS), AV(6), SMALL, Z0, QTS(3*NATOMS), OVEC(3), H1VEC(3), H2VEC(3),&
     &   DELE, SVH, DELTAT, RAT(3*NATOMS), AVG, DUMMY, INERTIA(3,3), &
     &   DIAGEXP(3*NATOMS), VNEW(3*NATOMS), E1, E2, SCALE, DELTAS, DOTOPT, PENERGY, POPT, PINIT
!    &   VECX(3*NATOMS), VECY(3*NATOMS), VECZ(3*NATOMS), DUMMYX, DUMMYY, DUMMYZ
      DOUBLE PRECISION DIAG(3*NATOMS), PSTEP(3*NATOMS), PFOB(3*NATOMS)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: QW
      LOGICAL KNOWE, KNOWG, KNOWH
      COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
      CHARACTER(LEN=5) ZSYMSAVE
      COMMON /SYS/ ZSYMSAVE
!
!  Assign enough memory to WORK for a blocksize of 32 to be possible. 
!  This is for DSYEVR.
!
      INTEGER ILWORK, LWORK, NFOUND, ISUPPZ(2*3*NATOMS)
      INTEGER IWORK(33*3*NATOMS)
      DOUBLE PRECISION WORK(33*3*NATOMS), ABSTOL

!     COMMON /VN/ VNEW(3*NATOMS)
      DOUBLE PRECISION LP,LP1,LP2
      CHARACTER ESTRING*87, GPSTRING*80, NSTRING*80, FSTRING*80
      COMMON /STRINGS/ ESTRING, GPSTRING, NSTRING, FSTRING
      LOGICAL PVFLAG
      COMMON /PVF/ PVFLAG
!
!  IVEC tells which eigenmode to follow (0,1- lowest eigenvalue,
!          2- next lowest eigenvalue NOPT- highest eigenvalue)
!  INR=0 - minimise
!  INR=1 - Newton-Raphson
!  INR=2 - transition state
!  INR=3 - minimise but don t do any reorientations (reaction paths)
!          and use pseudo-third derivative correction.
!  INR=4 - ts search but and use pseudo-third derivative correction.
!  INR=5 - as for INR=0 but move to principal axes first
!  INR=6 - Steepest descent - minimisation. 
!  INR=7 - Steepest descent - can stop at a saddle point.
!  INR=8 - Steepest descent/ascent search for a transition state
!  STPMAX is the maximum step size 
!  ISTCRT = 0 uses dot product step length (in Cartesian basis)
!           1 uses maximum Cartesian displacement as criterion
!           2 uses maximum e/vector displacement as criterion
!           3 (SD) - total arc length is STPMAX(1)
!           The above all use a single trust radius for dynamic scaling.  
!     10 scales in the EV basis using a trust radius for each direction. 
!
!  MODEDOWNT specifies that the mode in question should be followed downhill
!
      LWORK=33*3*NATOMS
      ILWORK=33*3*NATOMS
      INRSAVE=INR
      INR=INRIN
!
! for CHARMM: update nonbonded list at the start of each optimization
!
      IF(CHRMMT) CALL UPDATENBONDS(QTS)
!
!  Reset maximum step sizes in case this isn't the first call to EFOL.
!
      IF (DEBUG) PRINT '(A,G20.10)',' efol> resetting maximum step sizes to ',MXSTP
      DO J1=1,NOPT
         STPMAX(J1)=MXSTP
      ENDDO
!
!
!  This call to gmetry is needed for REDO runs with TIP potentials, otherwise if
!  an Euler angle goes bad the step is messed up.
!
!     CALL GMETRY(0,VEC,QTS)
      VEC(1:NOPT)=0.0D0
      MFLAG=.FALSE.
      MINTEST=.FALSE.
      TSTEST=.FALSE.
      SDTEST=.FALSE.
      NRTEST=.FALSE.
      FRAME=1

      IF (ZSYMSAVE(1:1).EQ.'W') ALLOCATE(QW(9*(NATOMS/2)))
      IF (INR.EQ.0) THEN
         IF (PTEST) WRITE(*,10) 
         MINTEST=.TRUE.
      ELSE IF (INR.EQ.1) THEN
         IF (PTEST) WRITE(*,20) 
         NRTEST=.TRUE.
      ELSE IF (INR.EQ.2) THEN
         IF (PTEST) WRITE(*,10) 
         TSTEST=.TRUE.
      ELSE IF (INR.EQ.3) THEN
         IF (PTEST) WRITE(*,10)
         MINTEST=.TRUE.
      ELSE IF (INR.EQ.4) THEN
         IF (PTEST) WRITE(*,10)
         TSTEST=.TRUE.
      ELSE IF (INR.EQ.6) THEN
         IF (PTEST) WRITE(*,30)
         SDTEST=.TRUE.
         MINTEST=.TRUE.
         CONVU=0.0D0
      ELSE IF (INR.EQ.7) THEN
         IF (PTEST) WRITE(*,30)
         SDTEST=.TRUE.
         NRTEST=.TRUE.
         CONVU=0.0D0
      ELSE IF (INR.EQ.8) THEN
         IF (PTEST) WRITE(*,30)
         SDTEST=.TRUE.
         TSTEST=.TRUE.
         CONVU=0.0D0
      ELSE
         PRINT '(A,I8)','efol> INR value not recognised,  set to 0 from ',INR
         INR=0
         MINTEST=.TRUE.
         IF (PTEST) WRITE(*,10) 
      ENDIF
10    FORMAT(' Updating structure with eigenvector-following steps.')
20    FORMAT(' Updating structure with pseudo-Newton-Raphson steps.')
30    FORMAT(' Updating structure with quadratic steepest-descent steps.')
      IF (MASST.AND.PTEST) PRINT*,'Mass weighting is active'
      ITER=1

40    IF (PTEST) PRINT*
      IF (PTEST) WRITE(*,60) ITER
60    FORMAT (' Beginning of optimization cycle ', I4,'.',/&
     &' -------------------------------------')
      IF (PV.AND.(.NOT.BFGSSTEP)) THEN
         CALL POTENTIAL(QTS,ENERGY,VNEW,.FALSE.,.FALSE.,RMS,PTEST,.FALSE.)
         PVFLAG=.FALSE.
         CALL PVOPT(QTS,ENERGY,VNEW)
      ENDIF
      CALL POTENTIAL(QTS,ENERGY,VNEW,.TRUE.,.TRUE.,RMS,PTEST,.FALSE.)
      CALL DUMPP(QTS,ENERGY)
!
!  GDIIS step if required. This really just doesn't work!
!
!     IF (DTEST) THEN
!        CALL  DIIS(NOPT,SVEC,ITER-1,QSAVE,VNEW,QTS,RMS,ENERGY,DONE)
!        IF (PTEST) WRITE(*,'(A,F20.10)') ' Energy before DIIS=',ENERGY
!        CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.TRUE.,RMS,PTEST,.FALSE.)
!        IF (PTEST) WRITE(*,'(A,F20.10)') ' Energy after DIIS=',ENERGY
!     ENDIF
!
!  Transformation to mass weighted coordinates if required.
!
      IF (MASST) CALL MASSWT(NATOMS,ATMASS,QTS,VNEW,.TRUE.)

      IF ((.NOT.VARIABLES).AND.(.NOT.MIEFT).AND.(.NOT.NOTRANSROTT).AND.(.NOT.PHI4MODT).AND. &
     &            (.NOT.(RINGPOLYMERT.AND.(RPSYSTEM(1:4).EQ.'AECK')))) THEN
         IF (ZSYM(NATOMS).EQ.'SY') THEN
            CALL SHIFTSTOCK(QTS,NATOMS)
         ELSEIF (RBAAT) THEN
            CALL SHIFTRIGID(QTS,NATOMS)
         ELSEIF (ZSYM(NATOMS).EQ.'TH') THEN
            CALL SHIFTHTH(QTS,NATOMS)
! hk286
         ELSEIF (GTHOMSONT) THEN
            CALL SHIFTHGTH(QTS,NATOMS)
         ELSE
            CALL SHIFTH(QTS,.TRUE.,NOPT,NATOMS,ATMASS)
         ENDIF
      ENDIF
      IF (FREEZE.AND.VARIABLES) THEN
         DO J1=1,NATOMS
            IF (FROZEN(J1)) HESS(J1,J1)=SHIFTL(1)
         ENDDO
      ENDIF

      IF (VARIABLES.AND.FOURDPBCT) CALL SHIFTFOURDPBCT
      IF (VARIABLES.AND.THREEDPBCT) CALL SHIFTTHREEDPBCT
      IF (VARIABLES.AND.TWODPBCT) CALL SHIFTTWODPBCT
      IF (VARIABLES.AND.ONEDPBCT) CALL SHIFTONEDPBCT
      IF (VARIABLES.AND.INVTONEDPBCT) CALL SHIFTONEDPBCT
      IF (VARIABLES.AND.INVTTWODPBCT) CALL SHIFTTWODPBCT
!      IF (VARIABLES.AND.EX1DT) CALL SHIFTEX1D
      IF (ISTCRT.NE.10) THEN
!        IF (PTEST) WRITE(*,'(A,4F20.10)') 'ENERGY,EOLD,ENERGY-EOLD,DELE=',ENERGY,EOLD,ENERGY-EOLD,DELE
         IF ((ENERGY-EOLD.NE.0.0D0).AND.(EOLD.NE.0.0D0)) THEN
            EPER=DABS((ENERGY-EOLD-DELE)/(ENERGY-EOLD))
            IF (EPER.GT.TRAD) THEN
               STPMAX(1)=MAX(STPMAX(1)/1.1D0,MINMAX)
            ELSE
               STPMAX(1)=MIN(STPMAX(1)*1.1D0,MAXMAX)
            ENDIF
            IF (PTEST) WRITE(*,'(3(A,F12.6))') ' Maximum step size=',STPMAX(1),' trust radius=',TRAD,' calculated ratio=', EPER
         ENDIF
      ENDIF
!
!  Diagonalize the relevant hessian and determine the number of
!  negative eigenvalues.
!
!  DSYEV and DSYEVR are supposed to order the eigenvalues and corresponding 
!  eigenvectors in ascending order, but we need the eigenvalues
!  in descending order. 
!

!     IF (LANCZOST) THEN
!        IF (PTEST) WRITE(*,55) ACCLAN, SHIFTLAN
!55      FORMAT(' Matrix diagonalisation by the Lanczos method, accuracy=',E12.5,' shift=',E12.5)
!C       CALL LANCZOS(NOPT,HESS,DIAG,1.0D0,-1.0D0,NFOUND)
!        IF (NFOUND.NE.NOPT) THEN
!           IF (PTEST) WRITE(*,56) NFOUND
!56         FORMAT(' *** WARNING *** Only ',I4,' eigenvalues found by Lanczos')
!           IF (ISTCRT.EQ.10) PRINT*,' SCALE 2 is recommended!'
!        ENDIF
!     ELSE
!        NFOUND=NOPT

      IF ((NUSEEV.GT.0).AND.(NUSEEV.LT.NOPT)) THEN
         CALL DSYEVR('V','I','U',NOPT,HESS,SIZE(HESS,1),0.0D0,1.0D0,1,NUSEEV,ABSTOL,NFOUND,DIAG, &
     &                        ZWORK,3*NATOMS,ISUPPZ,WORK, &
     &                        LWORK, IWORK, ILWORK, INFO )
              IF (INFO.NE.0) PRINT*,'WARNING - INFO=',INFO,' in DSYEVR'
!             PRINT '(A,F15.5,I10)','Optimal and actual value of LWORK=',WORK(1),LWORK
!             PRINT '(A,I10,I10)','Optimal and actual value of ILWORK=',IWORK(1),ILWORK
!
!  Put the NUSEEV eigenvalues and eigenvectors into the same places as they are expected to
!  be for full digagonalisation.
!  Note the different ordering! We have values for the lowest NUSEEV eigenvalues
!  in elements 1:NUSEEV, but we want to put them in NOPT-NUSEEV+1:NOPT
!  Do this using dummy DIAG values and sorting.
!
         HESS=0.0D0 
         DO I=1,NUSEEV
            DO J=1,NOPT
               HESS(J,I)=ZWORK(J,I)
            ENDDO
         ENDDO
         NEV=NUSEEV
         DIAG(NEV+1:NOPT)=1.0D100 ! put non-existent dummy eigenvalues in dummy elements, then sort
!        CALL EIGENSORT_VAL_ASC(DIAG,HESS,NOPT,3*NATOMS)
      ELSE
         CALL DSYEV('V','U',NOPT,HESS,SIZE(HESS,1),DIAG,TEMPA,9*NATOMS,INFO)
         IF (INFO.NE.0) PRINT*,'WARNING - INFO=',INFO,' in DSYEV'
!
!  The sort order given by DSYEV seems to vary with platform!
! 
         if (diag(1).lt.diag(nopt)) call eigensort_val_asc(diag,hess,nopt,3*natoms)
         NEV=NOPT
      ENDIF
!
!  Find eigenvalue demanded by value of ivec (first pass only)
!  and count number of negative eigenvalues in hessian (all passes).
!
      DO J=1,NOPT
         ZT(J)=.TRUE.
      ENDDO
!
!  Determine the zero eigenvalues. Assume the zeros have been
!  shifted to the top of the range.
!  However, some systems should have three zeros, some none
!  at all, and linear systems one less then usual.
!  Rotating systems have one!
!
      IF (ZSYM(1).EQ.'TH') THEN
         NZERO=3
         ZT(1)=.FALSE.
         ZT(2)=.FALSE.
         ZT(3)=.FALSE.
! hk286
      ELSEIF (GTHOMSONT .AND. (GTHOMMET < 5)) THEN
         NZERO=1
         ZT(1)=.FALSE.
      ELSEIF (GTHOMSONT .AND. (GTHOMMET .EQ. 5)) THEN
         NZERO=3+NATOMS
         ZT(1)=.FALSE.
         ZT(2)=.FALSE.
         ZT(3)=.FALSE.
      ELSEIF (ZSYM(1).EQ.'SY') THEN
         NZERO=(NATOMS/2)+6
         DO J1=1,NZERO
            ZT(J1)=.FALSE.
         ENDDO
      ELSE IF (PULLT.OR.EFIELDT) THEN
         NZERO=4
         ZT(1)=.FALSE.
         ZT(2)=.FALSE.
         ZT(3)=.FALSE.
         ZT(4)=.FALSE.
      ELSE IF (MKTRAPT.OR.MIEFT.OR.NOTRANSROTT) THEN
         NZERO=0
      ELSE IF ((EYTRAPT.OR.(ZSYM(1).EQ.'BE')).AND.(.NOT.TWOD)) THEN
         NZERO=3
         ZT(1)=.FALSE.
         ZT(2)=.FALSE.
         ZT(3)=.FALSE.
      ELSE IF (ZSYM(1).EQ.'CK') THEN
         DO J1=1,NATOMS+1
            ZT(J1)=.FALSE.
         ENDDO
         NZERO=NATOMS+1
      ELSE IF (RTEST) THEN
         ZT(1)=.FALSE.
         ZT(2)=.FALSE.
         NZERO=2
         IF (JZ.NE.0.0D0) THEN
            ZT(3)=.FALSE.
            ZT(4)=.FALSE.
            NZERO=4
         ENDIF
      ELSE IF (BULKT) THEN
         ZT(1)=.FALSE.
         ZT(2)=.FALSE.
         ZT(3)=.FALSE.
         NZERO=3
         IF (TWOD) THEN
            DO J1=1,NATOMS+2
               ZT(J1)=.FALSE.
            ENDDO
            NZERO=NATOMS+2
         ENDIF
      ELSE IF ((FPGRP.EQ.'DXh'.OR.FPGRP.EQ.'CXv').AND.(ZSYM(NATOMS)(1:1).NE.'W')) THEN
         ZT(1)=.FALSE.
         ZT(2)=.FALSE.
         ZT(3)=.FALSE.
         ZT(4)=.FALSE.
         ZT(5)=.FALSE.
         NZERO=5
      ELSE IF (VARIABLES.OR.PHI4MODT) THEN
         DO J1=1,NZERO
            ZT(J1)=.FALSE.
         ENDDO
      ELSE IF (RINGPOLYMERT) THEN
         DO J1=1,NZERO
            ZT(J1)=.FALSE.
         ENDDO
      ELSE IF (FIELDT) THEN
         ZT(1)=.FALSE.
         ZT(2)=.FALSE.
         ZT(3)=.FALSE.
         NZERO=3
      ELSE IF (TWOD) THEN
         IF (EYTRAPT.OR.(ZSYM(1).EQ.'BE')) THEN
            DO J1=1,NATOMS+1
               ZT(J1)=.FALSE.
            ENDDO
            NZERO=NATOMS+1
         ELSE
            DO J1=1,NATOMS+3
               ZT(J1)=.FALSE.
            ENDDO
            NZERO=NATOMS+3
         ENDIF
      ELSE IF (FREEZE) THEN
         NZERO=3*NFREEZE
         DO J1=1,NZERO
            ZT(J1)=.FALSE.
         ENDDO
      ELSEIF(GBT.OR.GBDT.OR.(PYGT.AND.UNIAXT).OR.((PYT.OR.PYGPERIODICT.OR.PYBINARYT.OR.MULTISITEPYT).AND.UNIAXT).OR.STOCKAAT) THEN

         IF (EFIELDT) THEN
            NZERO = NATOMS/2 + 4
         ELSE
            NZERO = NATOMS/2 + 6 
         ENDIF
         DO J1=1,NZERO
            ZT(J1)=.FALSE.
         ENDDO
      
      ELSE
         ZT(1)=.FALSE.
         ZT(2)=.FALSE.
         ZT(3)=.FALSE.
         ZT(4)=.FALSE.
         ZT(5)=.FALSE.
         ZT(6)=.FALSE.
         NZERO=6
      ENDIF
      DO J1=1,NOPT-NEV
         ZT(J1)=.FALSE.
      ENDDO

      IF (VALUEST.AND.(MOD(ITER-1,NVALUES).EQ.0)) THEN
         IF (PTEST) WRITE(*,'(A)') ' Eigenvalues of the Hessian matrix:'
         IF (PTEST) WRITE (*,'(6(G12.5,1X))') (DIAG(I),I=MAX(NZERO+1,NOPT-NEV+1),NOPT)
      ENDIF
      IF (VECTORST.AND.(MOD(ITER-1,NVECTORS).EQ.0)) THEN
         IF (PTEST) WRITE(*,'(A)') ' Eigenvectors of the Hessian matrix:'
         CALL HESSOUT(1)
      ENDIF
!     PRINT '(A)','Predicted zero eigenvectors for Thomson'
!     VECX(1:3*NATOMS)=0.0D0; VECY(1:3*NATOMS)=0.0D0; VECZ(1:3*NATOMS)=0.0D0
!     DO J1=1,3*NATOMS,2
!        VECX(J1)=SIN(QTS(J1+1))
!        VECY(J1)=COS(QTS(J1+1))
!     ENDDO
!     DO J1=2,3*NATOMS,2
!        IF (SIN(QTS(J1)).NE.0.0D0) THEN
!           VECX(J1)= COS(QTS(J1-1))*COS(QTS(J1))/SIN(QTS(J1-1))
!           VECY(J1)=-COS(QTS(J1-1))*SIN(QTS(J1))/SIN(QTS(J1-1))
!        ENDIF
!        VECZ(J1)=1.0D0   
!     ENDDO
!     CALL VECNORM(VECX,3*NATOMS)
!     CALL VECNORM(VECY,3*NATOMS)
!     CALL VECNORM(VECZ,3*NATOMS)
!     PRINT '(3F15.6)',(VECX(J1),VECY(J1),VECZ(J1),J1=1,3*NATOMS)
!     DO J1=1,3*NATOMS
!        DUMMYX=0.0D0
!        DUMMYY=0.0D0
!        DUMMYZ=0.0D0
!        DO J2=1,3*NATOMS
!           DUMMYX=DUMMYX+HESS(J2,J1)*VECX(J2)
!           DUMMYY=DUMMYY+HESS(J2,J1)*VECY(J2)
!           DUMMYZ=DUMMYZ+HESS(J2,J1)*VECZ(J2)
!        ENDDO
!        PRINT '(A,I6,A,3F20.10)','x,y,z dot products for eigenvector ',J1,' are ',DUMMYX,DUMMYY,DUMMYZ
!     ENDDO

      IF (EVCUT.NE.0.0D0) THEN
         DO J1=1,NOPT
            IF (DABS(DIAG(J1)).LT.EVCUT) THEN
               IF (PTEST) PRINT*,'Eigenvalue cutoff - no step will be taken for mode ',J1
               ZT(J1)=.FALSE.
            ENDIF
         ENDDO
      ENDIF
!
!  Count negative eigenvalues and find the smallest.
!
      INEG=0
      SMALL=1.0D20
      DO J=NOPT,NOPT-NEV+1,-1
         IF (ZT(J)) THEN
            IF (DIAG(J).LT.0.0D0) INEG=INEG+1
            IF (DIAG(J).LT.SMALL) THEN
               SMALL=DIAG(J)
            ELSE
               IF (DIAG(J).GT.0.0D0) GOTO 70
            ENDIF
         ENDIF
      ENDDO
70    CONTINUE

      IF (KEEPINDEX.AND.(ITER.EQ.1)) THEN
         HINDEX=INEG
         IF (PTEST) WRITE(*,'(A,I4)') ' Searching for a saddle with Hessian index=',HINDEX
      ENDIF
!
!  For transition state searches set IMODE to the appropriate mode to be
!  followed. If IVEC=0 this is the softest mode for each step. 
!  IMODE needs to be set here for:
!  (a) TSTEST and ITER=1,
!  (b) TSTEST and ITER>1 if IVEC=0.
!  (c) MINTEST or SDTEST and IVEC not 0.
!
      IMODE=0
      IASSIGN=0
      DONE=.FALSE.
      DO J=NOPT,1,-1
         IF (ZT(J)) THEN
            IASSIGN=IASSIGN+1
            IF ((IASSIGN.GE.ABS(IVEC)).AND.(.NOT.(DONE))) THEN
               IMODE=J
               DONE=.TRUE.
            ENDIF
         ENDIF
      ENDDO
!
!  On later passes, determine overlap between Hessian eigenvectors
!  and VEC (saved from previous step) for transition state searches.
!
      IF (ITER.GT.1.AND.TSTEST.AND.(IVEC.NE.0)) THEN
         Z0=0.0D0
         DO I=1,NOPT
            SVH=DOTOPT(VEC,HESS(1,I),NOPT)
            IF ((DABS(SVH).GT.Z0).AND.(ZT(I))) THEN
               Z0=DABS(SVH)
               IM=I
            ENDIF
         ENDDO

         IF (PTEST) WRITE(*,100) Z0,IM,DIAG(IM),SMALL
100      FORMAT(' Largest overlap=',F8.5,' for vector ',I4,' eigenvalue=',F14.7,' Smallest eigenvalue=',F14.7)
         IF (Z0.LT.0.8D0) THEN
            IF (IM.EQ.OMODE) THEN
               IF (PTEST) PRINT*,'Small overlap, but with same eigenvector'
            ELSE IF (IM.GT.OMODE) THEN
               IF (PTEST) WRITE(*,'(A)') ' Small overlap with softer eigenvector accepted'
            ELSE
               IF (ZT(OMODE)) THEN
                  IF (PTEST) PRINT*,'Small overlap, follow previous mode'
                  IM=OMODE

!                 IF (PTEST) PRINT*,'Small overlap - backtracking'
!                 DO J1=1,NOPT
!                    CSTEP(J1)=CSTEP(J1)/2.0D0
!                    PSTEP(J1)=PSTEP(J1)/2.0D0
!                    STPMAX(J1)=MAX(MINMAX,STPMAX(J1)/2.0D0)
!                    FOB(J1)=PFOB(J1)
!                    ZT(J1)=PZT(J1)
!                    QTS(J1)=QTS(J1)-CSTEP(J1)
!                 ENDDO
!                 IF (ALLOCATED(QW)) DEALLOCATE(QW)
!                 RETURN

               ENDIF
            ENDIF
!CC          IF (PTEST) PRINT*,'Small overlap - switch to softset mode'
!CC          IVEC=0
!CC          IM=IMODE
         ENDIF

         IF (IM.LT.OMODE-8) THEN
            IF (PTEST) PRINT*,'resetting mode followed to ',OMODE
            IM=OMODE
         ENDIF
!
!  This should save us from following a vector with a zero eigenvalue.
!  IMODE is set to a mode with non-zero eigenvalue before we enter this block.
!
         IF (ZT(IM)) THEN
            IMODE=IM
         ELSE IF (ZT(OMODE)) THEN
            IMODE=OMODE
         ENDIF
         OMODE=IMODE
         IF (PTEST) THEN
            IF (MODEDOWNT) THEN
               PRINT*,'Mode to be followed downhill=',IMODE
            ELSE
               PRINT*,'Mode to be followed uphill=',IMODE
            ENDIF
         ENDIF
      ELSE
         OMODE=IMODE
      ENDIF

      DO I=1,NOPT
         FOB(I)=0.0D0
         DO J=1,NOPT
            FOB(I)=FOB(I)+VNEW(J)*HESS(J,I)
         ENDDO
      ENDDO
!
!  Calculate the localisation index for each pair of normal modes and get the mean and
!  standard deviation. The sum is only over eigenvectors with negative (nonzero) eigenvalues.
!
!     DUMMY1=0.0D0
!     DUMMY2=0.0D0
!     DUMMY5=0.0D0
!     DUMMY6=0.0D0
!     NCONTRIB=0
!     NPCONTRIB=0
!     PDUMMY1=0.0D0
!     NPLUS=0
!     NMINUS=0
!     DO J1=1,NOPT
!        DUMMY3=0.0D0
!        DUMMY4=0.0D0
!        DO J3=1,NOPT
!           DUMMY3=DUMMY3+HESS(J3,J1)**2
!           DUMMY4=DUMMY4+HESS(J3,J1)**4
!        ENDDO
!        IF ((DIAG(J1).LT.0.0D0).AND.ZT(J1)) THEN
!           DUMMY5=DUMMY5+DUMMY3**2/DUMMY4
!           NMINUS=NMINUS+1
!           DO J2=J1+1,NOPT
!              IF ((DIAG(J2).LT.0.0D0).AND.ZT(J2)) THEN
!                 DUMMY3=0.0D0
!                 DO J3=1,NOPT
!                    DUMMY3=DUMMY3+ABS(HESS(J3,J1)*HESS(J3,J2))
!                 ENDDO
!                 DUMMY1=DUMMY1+DUMMY3
!                 DUMMY2=DUMMY2+DUMMY3**2
!                 NCONTRIB=NCONTRIB+1
!              ENDIF
!           ENDDO
!        ELSE IF (ZT(J1)) THEN
!           NPLUS=NPLUS+1
!           DUMMY6=DUMMY6+DUMMY3**2/DUMMY4
!           DO J2=J1+1,NOPT
!              IF ((DIAG(J2).GT.0.0D0).AND.ZT(J2)) THEN
!                 DUMMY3=0.0D0
!                 DO J3=1,NOPT
!                    DUMMY3=DUMMY3+ABS(HESS(J3,J1)*HESS(J3,J2))
!                 ENDDO
!                 PDUMMY1=PDUMMY1+DUMMY3
!                 NPCONTRIB=NPCONTRIB+1
!              ENDIF
!           ENDDO
!        ENDIF
!     ENDDO
!     PRINT*,'INEG*(INEG-1)/2,NCONTRIB=',INEG*(INEG-1)/2,NCONTRIB
!     IF (NCONTRIB.GT.1) THEN
!        WRITE(*,'(A,I9,5G20.10)') 'localisation ',
!    1              INEG,DUMMY1/NCONTRIB,SQRT((DUMMY2-DUMMY1**2/NCONTRIB)/(NCONTRIB-1)),
!    2              DUMMY5/MAX(1,NMINUS),DUMMY6/NPLUS,PDUMMY1/NPCONTRIB
!     ELSE IF (NCONTRIB.GT.0) THEN
!        WRITE(*,'(A,I9,5G20.10)') 'localisation ',
!    1              INEG,DUMMY1/NCONTRIB,0.0D0,DUMMY5/MAX(1,NMINUS),DUMMY6/NPLUS,PDUMMY1/NPCONTRIB
!     ELSE
!        WRITE(*,'(A,I9,5G20.10)') 'localisation ',
!    1              INEG,0.0D0,0.0D0,DUMMY5/MAX(1,NMINUS),DUMMY6/NPLUS,PDUMMY1/NPCONTRIB
!     ENDIF
!
!  Find the vector of STPMAX values by comparing predicted and
!  actual second derivatives for each eigenvector.
!
      IF (ITER.EQ.1) THEN
         DO J1=1,NOPT
            RAT(J1)=0.0D0  !  Initialise RAT - otherwise we can;t print it.
         ENDDO
      ENDIF
      IF ((ITER.GT.1).AND.(ISTCRT.EQ.10)) THEN
!         K1=-1
!         DO J1=1,10
!90          K1=K1+1
!            IF (.NOT.ZT(NOPT-K1)) GOTO 90
!            K2=-1
!            DO J2=1,10
!91             K2=K2+1
!               IF (.NOT.PZT(NOPT-K2)) GOTO 91
!C              TMAT(J2,J1)=ABS(( FOB(NOPT-K1)-PFOB(NOPT-K2))/(PSTEP(NOPT-K2)*DIAG(NOPT-K1))-1.0D0)
!               TMAT(J2,J1)=MIN( ABS(( FOB(NOPT-K1)-PFOB(NOPT-K2))/(PSTEP(NOPT-K2)*DIAG(NOPT-K1))-1.0D0),
!     1                          ABS((-FOB(NOPT-K1)-PFOB(NOPT-K2))/(PSTEP(NOPT-K2)*DIAG(NOPT-K1))-1.0D0) )
!            ENDDO
!         ENDDO
!         IF (PTEST) WRITE(*,'(A)') 'Matrix of trust ratios for five softest non-zero modes:'
!         IF (PTEST) WRITE(*,'(10G12.4)') ((TMAT(K2,K1),K2=1,10),K1=1,10)
            
         DO J1=1,NOPT
            TEMPA(J1)=STPMAX(J1)
         ENDDO
         K1=0
         DO J1=1,NOPT
            RAT(J1)=0.0D0
            IF (ZT(J1)) THEN
150            K1=K1+1
!              IF (PTEST) PRINT*,'K1,PZT,PFOB,PSTEP=',K1,PZT(K1),PFOB(K1),PSTEP(K1)
               IF ((.NOT.PZT(K1)).AND.(K1.LT.NOPT)) GOTO 150
               IF (DABS(PSTEP(K1)).GT.1.0D-40) THEN
!
!  Allow for possible phase change in the eigenvector. Just take the smaller value.
!
!  If the ts has a point group with degeneracies then these modes can mix 
!  differently between steps, which may cause issues with the trust ratios.
!
                  R1=( FOB(J1)-PFOB(K1))/(DIAG(J1)*PSTEP(K1))
                  R2=(-FOB(J1)-PFOB(K1))/(DIAG(J1)*PSTEP(K1))
                  RAT1=DABS(DABS(( FOB(J1)-PFOB(K1))/(DIAG(J1)*PSTEP(K1)))-1.0D0)
                  RAT2=DABS(DABS((-FOB(J1)-PFOB(K1))/(DIAG(J1)*PSTEP(K1)))-1.0D0)
                  RAT(J1)=MIN(RAT1,RAT2)
!                 IF (PTEST) WRITE(*,'(A,2I4,5E15.7)') 'J1,K1,FOB,PFOB,PSTEP,RAT1,DIAG=', 
!     &                                                 J1,K1,FOB(J1),PFOB(K1),PSTEP(K1),RAT1,DIAG(J1)
!                 IF (RAT(J1).GT.1.0D0) WRITE(*,'(A,2I4,5E15.7)') 'J1,K1,FOB,PFOB,PSTEP,RAT1,DIAG=',
!     &                                                            J1,K1,FOB(J1),PFOB(K1),PSTEP(K1),RAT1,DIAG(J1)
                  IF (RAT(J1).GT.TRAD) THEN
                     STPMAX(J1)=MAX(TEMPA(K1)/1.11D0,MINMAX)
                  ELSE
                     IF (MASST) THEN
                        SUM=0
                        DO J2=1,NATOMS
                           SUM=SUM+ATMASS(J2)
                        ENDDO
                        AVG=SQRT(SUM/NATOMS)
!                       IF (PTEST) PRINT *,'the average is',AVG
                        STPMAX(J1)=MIN(MAX(TEMPA(K1)*1.09D0,MINMAX),AVG*MAXMAX)
                     ELSE
                        STPMAX(J1)=MIN(MAX(TEMPA(K1)*1.09D0,MINMAX),MAXMAX)
                     ENDIF
                  ENDIF
               ELSE
                  STPMAX(J1)=MIN(MAX(TEMPA(K1),MINMAX),MAXMAX)
               ENDIF
            ELSE
            ENDIF
         ENDDO
      ENDIF

      IF (PGRAD.AND.(MOD(ITER-1,NGRADIENTS).EQ.0)) THEN
         IF (PTEST) WRITE(*,160)
160      FORMAT(' Gradients along Hessian eigenvectors: ')
         IF (PTEST) WRITE(*,'(6(F12.5,1X))')(FOB(I),I=MAX(NZERO+1,NOPT-NEV+1),NOPT)
      ENDIF
      IF (PTEST) WRITE(*,170)INEG
170   FORMAT(' Number of negative eigenvalues=',I3)
      IF (PTEST) WRITE(NSTRING,170)INEG
      IF (DUMPV) CALL VDUMP(DIAG,ZT,NOPT,3*NATOMS)
!
! Calculate step:
!
      SUM=0.0D0
      PROD=0.0D0
      DO J1=1,NOPT
         IF (ZT(J1)) THEN
            SUM=SUM+DABS(DIAG(J1))
            IF (DIAG(J1).GT.0.0D0) PROD=PROD+DLOG(DIAG(J1))
         ENDIF
      ENDDO
      IF ((FPGRP.EQ.'DXh'.OR.FPGRP.EQ.'CXv').AND.(ZSYM(NATOMS)(1:1).NE.'W')) THEN
         SUM=SUM/MAX(NOPT-5,1)
      ELSE
         SUM=SUM/MAX(NOPT-6,1)
      ENDIF
      IF (PTEST) WRITE(*,235) SUM,PROD
235   FORMAT(' Mean modulus of positive Hessian eigenvalues=',4X,F20.10,/,&
     &       ' Log product of positive Hessian eigenvalues =',4X,F20.10)
      IF (TSTEST) THEN
         IF (HINDEX.GE.1) THEN
            IF (PTEST) PRINT '(A,I8,A)',' efol> ',HINDEX,' modes will be searched uphill'
         ELSE
            IF (PTEST) WRITE(*,210) IMODE,DIAG(IMODE)
210         FORMAT(' Mode ',I4,' will be searched uphill. Eigenvalue=',4X,F19.10)
         ENDIF
!
! Save eigenvector being followed for use on next step
!
         DO I=1,NOPT
            VEC(I)=HESS(I,IMODE)
         ENDDO
!
!  If the convergence criteria are met then return now.
!
      ENDIF
      TEMP=-1.0D0
      RMS=0.0D0
      DO J1=1,NOPT
         RMS=RMS+VNEW(J1)**2
      ENDDO
      RMS=DSQRT(RMS/NOPT)
!
!  Take a step away from a stationary point along the appropriate
!  Hessian eigenvector. This enables us to start from converged minima.
!  Distinguish the case where we want to take a very small step away
!  from a transition state from others where we want a big displacement
!  to get unstuck. 
!
      AWAY=.FALSE.
      IF (RMS.LT.PUSHCUT) THEN
         IF (TSTEST.AND.(INEG.NE.HINDEX)) THEN
            IF (MOD(ITER-1,4).EQ.0) AWAY=.TRUE.
         ELSE IF (MINTEST.AND.(INEG.NE.0)) THEN
!
! This causes the step along the eigenvector with negative eigenvalue not to be scaled.
! Can lead to large energy rises in pathway calculations. Turn off 14/7/17.
! Paths look much smoother!
! However, we need the pushoff if it is the first step from a ts.
!           IF (MOD(ITER-1,4).EQ.0) AWAY=.TRUE. 
            IF (EFOLSTEP) AWAY=.TRUE. 
         ENDIF
      ENDIF
      IF ((.NOT.SDTEST).OR.AWAY) THEN
!
!  EF determination of steps
!
         ICOUNT=0
         DO I=NOPT,1,-1
            STEP(I)=0.0D0
            IF (ZT(I)) THEN
               IF (DABS(DIAG(NZERO+1)/MAX(DABS(DIAG(I)),1.0D-10)).GT.TEMP) TEMP=DABS(DIAG(NZERO+1)/MAX(DABS(DIAG(I)),1.0D-10)) 
               LP1=DABS(DIAG(I))/2.0D0
               IF (MASST) THEN
                  SUM=0
                  DO J2=1,NATOMS
                     SUM=SUM+ATMASS(J2)
                  ENDDO
                  AVG=SUM/NATOMS
                  LP2=1.0D0 + 4.0D0*((FOB(I)/DIAG(I))**2)/AVG
               ELSE
                  LP2=1.0D0 + 4.0D0*(FOB(I)/DIAG(I))**2
               ENDIF
               LP=LP1*(1.0D0+DSQRT(LP2))
               IF ((I.EQ.IMODE).AND.TSTEST) THEN
                  IF (MODEDOWNT) THEN
                     IF (PTEST) WRITE(*,'(A,I4,A,4X,F19.10)') ' Mode ',I,' will be searched downhill. Eigenvalue=',DIAG(I)
                  ELSE
                     IF (PTEST) WRITE(*,'(A,I4,A,4X,F19.10)') ' Mode ',I,' will be searched uphill. Eigenvalue=',DIAG(I)
                     LP=-LP
                     ICOUNT=ICOUNT+1
                  ENDIF
               ELSE IF (TSTEST.AND.(ICOUNT.LT.HINDEX).AND.(HINDEX.GE.1)) THEN
                  IF (PTEST) WRITE(*,'(A,I4,A,4X,F19.10)') ' Mode ',I,' will be searched uphill. Eigenvalue=',DIAG(I)
                  LP=-LP
                  ICOUNT=ICOUNT+1
               ENDIF
!
!  Pseudo-Newton-Raphson
!
               IF (NRTEST) LP=LP*DIAG(I)/DABS(DIAG(I))
               STEP(I)=-FOB(I)/LP
!
!  Minimise to remove zero eigenvalues ! DJW
!
!              IF (ABS(DIAG(I)).LT.1.0D-1) THEN
!                 STEP(I)=-FOB(I)*STPMAX(I)/ABS(FOB(I))
!                 PRINT '(A,I8,A,G20.10)',' efol> minimising for mode ',I,' step=',STEP(I)
!              ENDIF
            ENDIF
         ENDDO
!
!  Take action if we are heading for a stationary point of the wrong index.
!
         IF (AWAY) THEN
            IF (TSTEST) THEN
               IF (INEG.EQ.0) THEN
                  IF ((IVEC.GE.0).AND.PTEST) PRINT*,'Stepping away from minimum along mode ',IMODE,' + direction'
                  IF ((IVEC.LT.0).AND.PTEST) PRINT*,'Stepping away from minimum along mode ',IMODE,' - direction'
                  IF (PUSHOFF.NE.0.0D0) THEN
                     STEP(IMODE)=PUSHOFF
                  ELSE
                     STEP(IMODE)=STPMAX(IMODE)/10.0D0
                  ENDIF
                  IF (IVEC.LT.0) STEP(IMODE)=-STEP(IMODE)
               ELSE 
!
!  Step off along all the modes with negative eigenvalue except the smallest.
!
                  IF (IVEC.EQ.0) THEN
                     DO J1=1,NOPT-1
                        IF (ZT(J1).AND.(DIAG(J1).LT.0.0D0)) THEN
                           IF (PTEST) PRINT*,'Stepping away from higher order saddle along mode ',J1
                           IF (PUSHOFF.NE.0.0D0) THEN
                              STEP(J1)=PUSHOFF
                           ELSE
                              STEP(J1)=STPMAX(J1)/10.0D0
                           ENDIF
                        ENDIF
                     ENDDO
                  ELSE
!
!  Step off only along the mode specified by IMODE.
!
                     IF ((IVEC.GE.0).AND.PTEST) PRINT '(A,I6,A)',' efol> Stepping away from saddle along mode ',IMODE,' + direction'
                     IF ((IVEC.LT.0).AND.PTEST) PRINT '(A,I6,A)',' efol> Stepping away from saddle along mode ',IMODE,' - direction'
                     IF (PUSHOFF.NE.0.0D0) THEN
                        STEP(IMODE)=PUSHOFF
                     ELSE
                        STEP(IMODE)=STPMAX(IMODE)/10.0D0
                     ENDIF
                     IF (IVEC.LT.0) STEP(IMODE)=-STEP(IMODE)
                  ENDIF
               ENDIF
!           ELSE IF ((MINTEST.OR.SDTEST).AND.(.NOT.PUSHOPTT)) THEN
            ELSE IF (MINTEST.OR.SDTEST) THEN
!
!  Step off along all the modes with negative eigenvalues.
!
               IF (IVEC.EQ.0) THEN
                  DO J1=1,NOPT
                     IF (ZT(J1).AND.(DIAG(J1).LT.0.0D0)) THEN
                        IF (PTEST) PRINT '(A,I6)',' efol> Stepping away from saddle along mode ',J1
                        IF (PUSHOFF.NE.0.0D0) THEN
                           STEP(J1)=PUSHOFF
                        ELSE
                           STEP(J1)=STPMAX(J1)/10.0D0
                        ENDIF
                     ENDIF
                  ENDDO
               ELSE
!
!  Step off only along the mode specified.
!
                  IF ((IVEC.GE.0).AND.PTEST) PRINT '(A,I6,A)',' efol> Stepping away from saddle along mode ',IMODE,' + direction'
                  IF ((IVEC.LT.0).AND.PTEST) PRINT '(A,I6,A)',' efol> Stepping away from saddle along mode ',IMODE,' - direction'
!                   IF (PUSHOPTT) THEN
!                      PINIT=PUSHOFF
!                      IF (IVEC.LT.0) PINIT=-PUSHOFF
! !            ZWORK(J1,1) ??
!                      CALL GOLDEN(ZWORK(1:NOPT,IMODE),QTS,ENERGY,PINIT,POPT,PENERGY)
!                      IF (IVEC.GE.0) THEN
!                         WRITE(6,'(A,3G25.15)') 'path> golden section + pushoff, energy and delta energy: ', &
!      &                                 POPT,PENERGY,PENERGY-ENERGY
!                      ELSE
!                         WRITE(6,'(A,3G25.15)') 'path> golden section - pushoff, energy and delta energy: ', &
!      &                                 POPT,PENERGY,PENERGY-ENERGY
!                      ENDIF
!                      STEP(IMODE)=POPT
!                   ELSE 
                     IF (PUSHOFF.NE.0.0D0) THEN
                        STEP(IMODE)=PUSHOFF
                    ELSE
                        STEP(IMODE)=STPMAX(IMODE)/10.0D0
                    ENDIF
                    IF (IVEC.LT.0) STEP(IMODE)=-STEP(IMODE)
!                 ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF (EFSTEPST.AND.(MOD(ITER-1,EFSTEPS).EQ.0).AND.(.NOT.SDTEST)) THEN
            DO I=NZERO+1,NOPT
                IF (PTEST) WRITE(*,360) I, STEP(I)
360             FORMAT(' Unscaled step for mode ',I4,'=',F20.10)
            ENDDO
         ENDIF
         IF ((TEMP.GT.1.0D0).AND.PTEST) WRITE(*,366) TEMP
366      FORMAT(' Largest modulus ratio of non-zero e/values= ',5X,G20.10)   

      ELSE
!
!  Steepest descent/ascent step - Page-McIver method
!
!  Outline of method:
!  Value of parameter t is determined by the arc length via
!  a differential equation (33) from the PM paper.
!  If all eigenvalues are +ve then the integral is bounded.
!  To converge to a saddle point the gradient must have no component
!  in eigendirections with negative eigenvalues. Such cases can
!  obviously cause numerical problems!
!
!  Zero gradient components should be conserved, so try setting steps
!  to zero if the gradient component is less than CONVR and INR=7.
!
!  STPMAX(1) is dynamically adjusted via a trust radius scheme.
!     
         DELTAT=STPMAX(1)/(100.0D0*RMS*SQRT(1.0D0*NOPT))
         DELTAS=RMS*SQRT(1.0D0*NOPT)*DELTAT/2.0D0
         DO J1=1,NOPT
            DIAGEXP(J1)=DEXP(-2*DIAG(J1)*DELTAT)
         ENDDO
         IF (INR.EQ.8) DIAGEXP(IMODE)=DEXP(2*DIAG(IMODE)*DELTAT)
         J1=0

666      J1=J1+1
         DELTASP=DELTAS
         TEMP=0.0D0
         DO J2=1,NOPT
            IF (ZT(J2)) TEMP=TEMP+FOB(J2)**2*DIAGEXP(J2)**J1
         ENDDO
         DELTAS=DELTAS+DSQRT(TEMP)*DELTAT
!
!  We need to escape from the loop if the integral has effectively converged
!  or if the arc length exceeds the allowed value.
!
         IF ((DELTAS.LT.STPMAX(1)).AND.(J1.LT.100000).AND.((DELTAS-DELTASP)/DELTASP.GT.1.0D-10)) GOTO 666

         TPAR=J1*DELTAT
         IF (PTEST) WRITE(*,'(A,G15.6,A,G15.6,A,I6)') ' efol> Estimated arc length=',DELTAS,' for t=',TPAR,' integration steps=',J1

         STEP(1:NOPT)=0.0D0
         DO J1=NOPT-NEV+1,NOPT 
            IF (.NOT.((DABS(FOB(J1)).LT.CONVR).AND.(INR.EQ.7)).AND.ZT(J1)) &
  &            STEP(J1)=FOB(J1)*(DEXP(-DIAG(J1)*TPAR)-1)/DIAG(J1)
         ENDDO
         IF (INR.EQ.8) STEP(IMODE)=FOB(IMODE)*(DEXP(DIAG(IMODE)*TPAR)-1)/DIAG(IMODE)
         STPMAG=DSQRT(DOTOPT(STEP(1),STEP(1),NOPT))
         SSTPMAG=STPMAG
         SCALE=1.0D0
         IF (PTEST) WRITE(*,'(A,2F12.6)') ' efol> % of step and gradient along softest mode=',  &
     &                                       ABS(STEP(NOPT))*100.0D0/STPMAG,ABS(FOB(NOPT))*100.0D0/(RMS*SQRT(1.0D0*NOPT))
      ENDIF
!
!  Calculate parts of the numerically predicted energy change before
!  scaling for possible dynamic adjustment of the step size. The
!  predicted energy change is calculated below once we have scaled
!  everything.
!
      E1=0.0D0
      DO J1=1,NOPT
         E1=E1+STEP(J1)*FOB(J1)
      ENDDO
      E2=0.0D0
      DO J1=1,NOPT
         E2=E2+DIAG(J1)*STEP(J1)**2
      ENDDO
      E2=E2/2.0D0

!
!  Calculate scaling factor from the steps in the normal mode
!  basis. Seems to work best for scaling of step in ev basis
!  for Cartesians and coordinate basis for internals.
!
!  Scale according to step size in ev basis:
!
      IF (ISTCRT.EQ.10) THEN
         CALL VSTAT(STEP(1),AV,NOPT,3*NATOMS)
         STPMAG=MAX(AV(1),1D-10)
         DO J1=1,NOPT
            SCALE=MIN(STPMAX(J1)/MAX(DABS(STEP(J1)),1D-10),1.0D0)
            IF ((J1.NE.IMODE).OR.(.NOT.AWAY)) STEP(J1)=SCALE*STEP(J1)
            IF (ABS(DIAG(J1)).LT.EVCUT) THEN
               STEP(J1)=0.0D0
!              STEP(J1)=-FOB(J1)*MAXMAX/ABS(FOB(J1))
!              PRINT '(A,I8,A,G20.10)',' efol> eigenvalue ',J1,' is too small; push off downhill ',STEP(J1)
!              PRINT '(A,I8,4G20.10)','J1,FOB,DIAG,MAXMAX,STEP=',J1,FOB(J1),DIAG(J1),MAXMAX,STEP(J1)
            ENDIF
         ENDDO
         CALL VSTAT(STEP(1),AV,NOPT,3*NATOMS)
         SSTPMAG=MAX(AV(1),1D-10)
         SCALE=1.0D0
      ELSE IF ((ISTCRT.EQ.3).AND.AWAY) THEN
         STPMAG=PUSHOFF
         SSTPMAG=PUSHOFF
         SCALE=1.0D0
      ENDIF

      IF (ISTCRT.EQ.10) THEN
         IF (SUMMARYT.AND.(MOD(ITER-1,NSUMMARY).EQ.0)) THEN
            IF (PTEST) WRITE(*,290)
290         FORMAT(79('-'))
            IF (PTEST) WRITE(*,'(A)') 'Vector      Gradient        Secder       Step          Max step    Trust ratio'
            IF (PTEST) WRITE(*,290)
            DO I=NOPT-NEV+1,NOPT
               !              IF (ZT(I)) THEN
               IF (PTEST) WRITE(*,310) I,FOB(I),DIAG(I),STEP(I),STPMAX(I),RAT(I)
310            FORMAT(I4,1X,E15.6,1X,E15.6,1X,E13.6,1X,E12.6,1X,E15.6)
!              ENDIF
            ENDDO
            IF (PTEST) WRITE(*,290)
         ENDIF
      ELSE IF (ISTCRT.EQ.3) THEN
         IF (SUMMARYT.AND.(MOD(ITER-1,NSUMMARY).EQ.0)) THEN
            IF (PTEST) WRITE(*,291)
291         FORMAT(51('-'))
            IF (PTEST) WRITE(*,'(A)') 'Vector      Gradient        Secder       Step'
            IF (PTEST) WRITE(*,291)
            DO I=NOPT-NEV+1,NOPT
!              IF (ZT(I)) THEN
                  IF (PTEST) WRITE(*,'(I4,1X,E15.6,1X,E15.6,1X,E13.6)') I,FOB(I),DIAG(I),STEP(I)
!              ENDIF
            ENDDO
            IF (PTEST) WRITE(*,291)
         ENDIF
      ELSE IF (ISTCRT.EQ.2) THEN
         CALL VSTAT(STEP(1),AV,NOPT,3*NATOMS)
         SCALE=MIN(STPMAX(1)/MAX(AV(1),1D-10),1.D0)
         STPMAG=MAX(AV(1),1D-10)
         CALL SCDOT(SCALE,STEP(1),3*NATOMS)
         SSTPMAG=SCALE*STPMAG
      ENDIF
!
!  Find the step in the Cartesian basis and store in CSTEP
!
      DO J=1,NOPT
         DUMMY=0.0D0
         DO I=1,NOPT
            DUMMY=DUMMY+STEP(I)*HESS(J,I)
         ENDDO
         CSTEP(J)=DUMMY
      ENDDO
!
!  Scaling if required in the Cartesian basis - according to the
!  maximum component
!
      IF (ISTCRT.EQ.1) THEN
         CALL VSTAT(CSTEP(1),AV(1),NOPT,6*NATOMS)
         SCALE=MIN(STPMAX(1)/MAX(AV(1),1D-10),1.D0)
         STPMAG=MAX(AV(1),1D-10)
         CALL SCDOT(SCALE,CSTEP(1),6*NATOMS)
         SSTPMAG=SCALE*STPMAG
!
!   Scaling in the Cartesian basis in terms of the total step.
!
      ELSE IF (ISTCRT.EQ.0) THEN
         CALL VSTAT(CSTEP(1),AV(1),NOPT,6*NATOMS)
         STPMAG=DSQRT(DOTOPT(CSTEP(1),CSTEP(1),NOPT))
         SCALE=MIN(STPMAX(1)/MAX(STPMAG,1D-10),1.D0)
         CALL SCDOT(SCALE,CSTEP(1),6*NATOMS)
         SSTPMAG=SCALE*STPMAG
      ENDIF
      IF ((.NOT.SDTEST).AND.PTEST) WRITE(*,'(A,12X,2G20.10)') ' The maximum scaled/unscaled step is: ',SSTPMAG,STPMAG
!
!  We can now check for convergence since we know the RMS force, the
!  Hessian index and the magnitude of the proposed step. If converged,
!  we don't take the step.
!
      IF ((CONVU.EQ.0.0D0).OR.SDTEST) THEN
         TEST1=RMS.LT.CONVR
      ELSE IF (CONVR.EQ.0.0D0) THEN
         TEST1=STPMAG.LT.CONVU
      ELSE
         TEST1=(STPMAG.LT.CONVU).AND.(RMS.LT.CONVR)
      ENDIF
      IF (INDEXT) THEN
         TEST2=NRTEST.OR.(MINTEST.AND.(INEG.EQ.0)).OR.(TSTEST.AND.(INEG.EQ.HINDEX))
      ELSE
         TEST2=.TRUE.
      ENDIF
      IF (TEST1.AND.TEST2) MFLAG=.TRUE.
      IF (ITER.LT.NSTEPMIN) MFLAG=.FALSE.
      IF (PV.AND.(.NOT.PVFLAG)) MFLAG=.FALSE.
!
!  Don't call symmetry if we're doing Fenske-Hall, or if the RMS force is too high.
!
      IF ((ZSYM(NATOMS).NE.'FH').AND.(.NOT.VARIABLES).AND.(.NOT.FIELDT).AND.(.NOT.RINGPOLYMERT) &
     &        .AND.(.NOT.PHI4MODT) &
     &        .AND.(PTEST).AND.(((RMS.LT.SYMCUT).OR.(ITER.EQ.NSTEPS)).OR.MFLAG)) THEN
         IF (ZSYMSAVE(1:1).EQ.'W') THEN
!           DO J2=1,NOPT
!              QSAVE(J2)=QTS(J2)
!           ENDDO
!           DO J2=1,NATOMS  !  WCOMMENT
            DO J2=1,NATOMS/2
               CALL CONVERT(QTS(3*(J2-1)+1),QTS(3*(J2-1)+2),QTS(3*(J2-1)+3),&
!    &                      QTS(3*(NATOMS+J2-1)+1),QTS(3*(NATOMS+J2-1)+2),QTS(3*(NATOMS+J2-1)+3),&
     &                      QTS(3*(NATOMS/2+J2-1)+1),QTS(3*(NATOMS/2+J2-1)+2),QTS(3*(NATOMS/2+J2-1)+3),&
     &                      OVEC,H1VEC,H2VEC)
               QW(9*(J2-1)+1)=OVEC(1)
               QW(9*(J2-1)+2)=OVEC(2)
               QW(9*(J2-1)+3)=OVEC(3)
               QW(9*(J2-1)+4)=H1VEC(1)
               QW(9*(J2-1)+5)=H1VEC(2)
               QW(9*(J2-1)+6)=H1VEC(3)
               QW(9*(J2-1)+7)=H2VEC(1)
               QW(9*(J2-1)+8)=H2VEC(2)
               QW(9*(J2-1)+9)=H2VEC(3)
            ENDDO
!           NATOMS=NATOMS*3  ! WCOMMENT
            NATOMSSAVE=NATOMS
            NATOMS=(NATOMS/2)*3
            CALL SYMMETRY(HORDER,.TRUE.,QW,INERTIA)
!           NATOMS=NATOMS/3  ! WCOMMENT
            NATOMS=2*(NATOMS/3)
            NATOMS=NATOMSSAVE
!           DO J2=1,NOPT
!              QTS(J2)=QSAVE(J2)
!           ENDDO
         ELSE
            CALL SYMMETRY(HORDER,.TRUE.,QTS,INERTIA) 
         ENDIF
      ENDIF
      IF (MFLAG) THEN
         IF (ALLOCATED(QW)) THEN
            DEALLOCATE(QW)
         ENDIF
         GOTO 1111
      ENDIF
!
!  Calculate predicted change in energy. Take the step using VADD.
!
      DELE=E1*SCALE+E2*SCALE**2
      CALL VADD(QTS,QTS,CSTEP,NOPT,1)
      KNOWE=.FALSE.
      KNOWG=.FALSE.
      KNOWH=.FALSE.
!
! Summarize
!
      IF ((ISTCRT.NE.10).AND.(ISTCRT.NE.3)) THEN
         IF (PTEST) WRITE(*,'(A,F8.5)') ' Scale factor set to: ',SCALE
         IF (PTEST) WRITE(*,'(A,F15.6)') ' Maximum step length allowed is ',STPMAX(1)
         IF (SUMMARYT.AND.(MOD(ITER-1,NSUMMARY).EQ.0)) THEN
            IF (PTEST) WRITE(*,490)
            IF (PTEST) WRITE(*,480)
480         FORMAT('Parameter',T20,'dV/dR',T32,'Step',T46,'Rold',T56,'Rnew')
            IF (PTEST) WRITE(*,490)
490         FORMAT(64('-'))
            DO I=1,NOPT
               IF (ZT(I)) THEN
                  IF (PTEST) WRITE(*,500)' ',VNEW(I),CSTEP(I),(QTS(I)-CSTEP(I)),QTS(I)
500              FORMAT(T7,A,T17,F10.6,T29,F10.5,T41,F10.5,T53,F10.5)
               ENDIF
            ENDDO
            IF (PTEST) WRITE(*,490)
         ENDIF
      ENDIF

      IF ((FIXAFTER.GT.0).AND.(ITER.GT.FIXAFTER)) FIXIMAGE=.TRUE.
      CALL FLUSH(6)
!
!  Undo mass weighting if necessary.
!
      IF (MASST) THEN
         DO J1=1,NATOMS
            AMASS=1/SQRT(ATMASS(J1))
            QTS(3*J1-2)=AMASS*QTS(3*J1-2)
            QTS(3*J1-1)=AMASS*QTS(3*J1-1)
            QTS(3*J1)=AMASS*QTS(3*J1)
         ENDDO
      ENDIF

      CALL FLUSH(6)
! What was gmetry doing?? DJW
!
!     CALL GMETRY(ITER,VEC,QTS)
!
!     Print out the coordinates and distance matrix
!
      IF (ADMT.AND.(MOD(ITER-1,NADM).EQ.0)) CALL ADM(QTS)
      EOLD=ENERGY
      DO J1=1,NOPT
         PSTEP(J1)=STEP(J1)
         PFOB(J1)=FOB(J1)
         PZT(J1)=ZT(J1)
      ENDDO
      ITER=ITER+1
      IVEC=IVEC2
      IF (ITER.GE.ITMAX+1) THEN
         IF (ALLOCATED(QW)) THEN
            DEALLOCATE(QW)
         ENDIF
         GOTO 1111
      ENDIF

      GOTO 40

! 1111  if (present(INRin)) then ! SAT

1111  INR=INRSAVE
!     endif

      END SUBROUTINE EFOL
END MODULE MODEFOL
