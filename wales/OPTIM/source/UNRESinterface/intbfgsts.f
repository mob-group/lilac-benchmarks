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
C  This subroutine is designed to perform a hybrid eigenvector-following/gradient
C  minimization optimization for large systems where we do not want to 
C  diagonalize the whole Hessian.
C  This subroutine (as opposed to bfgsts) works in internal coordinates with the 
C  unres potential. jmc March 2003
C
C***********************************************************************
C
      SUBROUTINE INTBFGSTS(ITMAX,COORDS,ENERGY,GRAD,MFLAG,RMS,EVALMIN,EVALMAX,VECS,ITER,POTCALL,PTEST)
      USE COMMONS
      USE KEY
      USE VECCK
      USE ZWK
      USE MODUNRES
      USE MODHESS
      IMPLICIT NONE

C     INCLUDE 'lparams.h'
      INTEGER J1, J2, J4, INEG, J, ITER, NS, I, K1, NBFGS, ITMAX, FRAME, ITDONE
      DOUBLE PRECISION SCRATCH(6*NATOMS),GRAD(3*NATOMS),ENERGY,COORDS(3*NATOMS),SUM,AVG,FOBNEW,DUMMY,
     1                 RMS,FOB,PSTEP,AX,BX,TOL,PSTEPNEW,DPRAND,FIXDIR(3*NATOMS),TEMPA(9*NATOMS),
     2                 EVALMAX, EVALMIN, STPMAG,COORDSN(3*NATOMS),VECSP(3*NATOMS),
     3                 XP1,XP2,E1,E2,DELE,SCALE,EPER,EOLD,RMS2,VECS(3*NATOMS),SOVER,EREAL,XRAT(3*NATOMS),
     4                 XFOB(3*NATOMS),STEP(3*NATOMS),XPSTEP(3*NATOMS),CSTEP(3*NATOMS),PFOB(3*NATOMS),PCSTEP(3*NATOMS),
     5                 AV(6),SSTPMAG,RAT1,RAT2,TEMP,LP1,LP2,LP,VECL(3*NATOMS)
C     EXTERNAL OP, IOVECT
      CHARACTER(LEN=80) :: FNAME
      LOGICAL MFLAG, AWAY, LINE, PVFLAG, RESET, STEST, POTCALL, PTEST, FIXDSAVE
      COMMON /PVF/ PVFLAG
      LOGICAL CONNECTT, DUMPPATH, READPATH, CALCRATES, STOPFIRST, CONVERGED
      INTEGER NCONNECT
      DOUBLE PRECISION TEMPERATURE, HRED
      COMMON /CONN/ STOPFIRST, CONNECTT, NCONNECT, DUMPPATH, READPATH, CALCRATES, TEMPERATURE, HRED
      DOUBLE PRECISION TMPINT(NINTS) 
C     COMMON /FORPATH/ INTSTEP ! jmc for use in path, now in modunres

C
C  Assign enough memory to WORK for a blocksize of 32 to be possible.
C  This is for DSYEVR.
C
      INTEGER ILWORK, LWORK, NFOUND, ISUPPZ(2*3*NATOMS)
      INTEGER IWORK(33*3*NATOMS), INFO, ISTAT
      DOUBLE PRECISION WORK(33*3*NATOMS), ABSTOL, DIAG(3*NATOMS), DLAMCH
      LOGICAL KNOWE, KNOWG, KNOWH
      COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
      SAVE

      LWORK=33*3*NATOMS
      ILWORK=33*3*NATOMS

      IF (.NOT.NOHESS) THEN
         WRITE(*,'(A)') '** WARNING - INTMIN only works with NOHESS at present!!'
         WRITE(*,'(A)') '** Stopping...'
         STOP
      ENDIF

      IF (FIXD) THEN
          WRITE(*,'(A)') '** WARNING - INTMIN and FIXD set: incompatible, therefore stopping.'
          STOP
      ENDIF

      IF ((ZSYM(1)(1:1).EQ.'W').AND.(.NOT.BFGSSTEP)) THEN
         PRINT*,'BFGSTS procedures have not been programmed for TIP potentials'
         STOP
      ENDIF
      FIXDSAVE=FIXD
      IF ((HINDEX.GT.1).AND.(.NOT.NOIT)) THEN
         WRITE(*,'(A)') 'To search for higher index saddles you must use NOIT or SEARCH 2'
         STOP
      ENDIF

      DO J2=1,nres
         c(1,J2)=COORDS(6*(J2-1)+1)
         c(2,J2)=COORDS(6*(J2-1)+2)
         c(3,J2)=COORDS(6*(J2-1)+3)
         c(1,J2+nres)=COORDS(6*(J2-1)+4)
         c(2,J2+nres)=COORDS(6*(J2-1)+5)
         c(3,J2+nres)=COORDS(6*(J2-1)+6)
      END DO
      CALL UPDATEDC
      CALL int_from_cart(.true.,.false.)
      CALL chainbuild

      ITER=1
      frame=1
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
            READ(44,*) (VECS(J1),J1=1,NINTS)
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
C jmc remember vector.dump contains eigenvector in internal coordinates...
            READ(45,*) (VECS(J1),J1=1,NINTS)
            GOTO 20
10          CLOSE(45)
         ENDIF
         WRITE(*,'(A,F20.10)') ' Reaction vector read successfully. Eigenvalue=   ',EVALMIN
         NS=100
      ENDIF

      STEST=.TRUE.
      IF (NOHESS) STEST=.FALSE.
      IF (READV.AND.(ITER.EQ.1)) STEST=.FALSE.
90    IF (PTEST) PRINT*
      NUP=HINDEX
      IF ((.NOT.BFGSSTEP).AND.PTEST) WRITE(*,11) ITER
11          FORMAT (' intBFGSTS> Beginning of optimization cycle ', I4,'.',/
     1              ' ------------------------------------------')
      FIXIMAGE=.FALSE.
      IF ((FIXAFTER.GT.0).AND.(ITER.GT.FIXAFTER)) FIXIMAGE=.TRUE.
      IF (POTCALL) THEN
         IF (PV.AND.(.NOT.BFGSSTEP)) THEN
            IF (.NOT.KNOWE) CALL POTENTIAL(COORDS,ENERGY,GRAD,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            PVFLAG=.FALSE.
            CALL PVOPT(COORDS,ENERGY,GRAD)
         ENDIF
         IF ((.NOT.KNOWG).OR.(STEST.AND.(.NOT.KNOWH))) CALL POTENTIAL(COORDS,ENERGY,GRAD,.TRUE.,STEST,RMS,PTEST,.FALSE.)
      ENDIF
      CALL DUMPP(COORDS,ENERGY)
      IF ((.NOT.VARIABLES).AND.NOIT.AND.STEST) CALL SHIFTH(COORDS,.TRUE.,NOPT,NATOMS,ATMASS)

C     DO J1=1,NOPT
C        VECSP(J1)=0.0D0  ! otherwise VECSP is not initialised
C     ENDDO

C jmc Put internals into COORDS array here if we've come here from path, as we don't go through 
C the do loop below where this would be done normally.
      IF (BFGSSTEP) CALL geom_to_var(NINTS,COORDS(1:NINTS))

      IF (.NOT.((READV.AND.BFGSSTEP).OR.(.NOT.POTCALL))) THEN
         IF (FIXD.AND.(ITER.EQ.1)) THEN
            IF (CONNECTT) THEN
               DO J1=1,NOPT
                  FIXDIR(J1)=VECS(J1)
               ENDDO
            ELSE
               DO J1=1,NOPT
                  DUMMY=2*(DPRAND()-0.5D0)
                  FIXDIR(J1)=DUMMY
               ENDDO
            ENDIF
            CALL ORTHOGOPT(FIXDIR,COORDS,.TRUE.)
            IF (.TRUE.) CALL HSMOVE(COORDS,COORDSN,FIXDIR,PSTEP,PTEST)
            EVALMIN=1.0D0
         ELSE IF (FIXD) THEN
            CALL ORTHOGOPT(FIXDIR,COORDS,.TRUE.)
            IF (.TRUE.) CALL HSMOVE(COORDS,COORDSN,FIXDIR,PSTEP,PTEST)
            EVALMIN=1.0D0
         ENDIF
         IF (.NOT.NOHESS) THEN 
            IF (.NOT.NOIT) THEN
!              CALL ITEIG(ITER,COORDS,VECS,EVALMIN,EVALMAX,NS,SOVER,PTEST,VECL,CONVERGED)
               PRINT'(A)',' Iterative scheme for eigenvalue not available with Hessian in internals'
               STOP
            ELSE
               ABSTOL=DLAMCH('Safe  minimum')
               IF (ITER.GT.1) THEN
                  DO J1=1,NOPT
                     VECSP(J1)=VECS(J1)
                  ENDDO
               ENDIF
               CALL DSYEVR('V','I','U',NOPT,HESS,SIZE(HESS,1),0.0D0,1.0D0,1,HINDEX,ABSTOL,NFOUND,DIAG,ZWORK,3*NATOMS,ISUPPZ,WORK,
     1                        LWORK, IWORK, ILWORK, INFO )
               IF (INFO.NE.0) PRINT*,'WARNING - INFO=',INFO,' in DSYEVR'
C              PRINT*,'Optimal and actual values of LWORK=',WORK(1),LWORK
C              PRINT*,'Optimal and actual values of ILWORK=',IWORK(1),ILWORK
               EVALMIN=DIAG(1)
               SOVER=0.0D0
               DO J1=1,NOPT
                  VECS(J1)=ZWORK(J1,1)
                  IF (ITER.GT.1) SOVER=SOVER+VECS(J1)*VECSP(J1)
               ENDDO
               IF (PTEST) WRITE(*,'(A,F15.7,A,F15.7)') ' Smallest eigenvalue=',EVALMIN,' overlap with previous vector=',SOVER

               DO I=1,NOPT
                  SCRATCH(I) = COORDS(I)
               ENDDO
               DO I=1,HINDEX
                  XFOB(I)=0.0D0
                  DO J = 1,NOPT
                     XFOB(I)=XFOB(I)+GRAD(J)*ZWORK(J,I)
                  ENDDO
               ENDDO
C
C  Find the vector of STPMAX values by comparing predicted and
C  actual second derivatives for each eigenvector.
C
               IF ((ITER.GT.1).AND.(ISTCRT.EQ.10)) THEN
                  DO J1=1,NOPT
                     TEMPA(J1)=STPMAX(J1)
                  ENDDO
                  K1=0
                  DO J1=1,HINDEX
                     XRAT(J1)=0.0D0
                     K1=K1+1
                     IF (DABS(XPSTEP(K1)).GT.1.0D-40) THEN
C
C  Allow for possible phase change in the eigenvector. Just take the smaller value.
C
                        RAT1=DABS(( XFOB(J1)-PFOB(K1))/(DIAG(J1)*XPSTEP(K1))-1.0D0)
                        RAT2=DABS((-XFOB(J1)-PFOB(K1))/(DIAG(J1)*XPSTEP(K1))-1.0D0)
                        XRAT(J1)=MIN(RAT1,RAT2)
C                       WRITE(*,'(A,2I4,5E15.7)') 'J1,K1,FOB,PFOB,PSTEP,RAT1,DIAG=', J1,K1,XFOB(J1),PFOB(K1),PSTEP(K1),RAT1,DIAG(J1)
                        IF (XRAT(J1).GT.TRAD) THEN
                           STPMAX(J1)=MAX(TEMPA(K1)/1.11D0,MINMAX)
                        ELSE
                           IF (MASST) THEN
                              SUM=0
                              DO J2=1,NATOMS
                                 SUM=SUM+ATMASS(J2)
                              ENDDO
                              AVG=SQRT(SUM/NATOMS)
C                             PRINT *,'the average is',AVG
                              STPMAX(J1)=MIN(MAX(TEMPA(K1)*1.09D0,MINMAX),AVG*MAXMAX)
                           ELSE
                              STPMAX(J1)=MIN(MAX(TEMPA(K1)*1.09D0,MINMAX),MAXMAX)
                           ENDIF
                        ENDIF
                     ELSE
                        STPMAX(J1)=MAX(TEMPA(K1),MINMAX)
                     ENDIF
                  ENDDO
               ENDIF
               TEMP=-1.0D0
            ENDIF
C
C  Do not use.
C
C           CALL ITEIGN(ITER,COORDS,VECS,EVALMIN,EVALMAX,PTEST)
C
C  Lanczos routine - does work but seems to be slower? For NFIG=3 it is only
C  a bit slower than ITEIG with CEIG=0.01.
C
C           NVAL=-10
C           ANV=ABS(NVAL)
C           NFIG=3
C           NPERM=0
C           IF (ITER.GT.1) NPERM=1
C           MAXOP=NEVS
C           NBLOCK=1
C           IF ((ANV.GT.MANV).OR.(NBLOCK.GT.MAXBLOCK)) THEN
C              WRITE(*,'(A)') ' Too many eigenvectors or too large a Lanczos blocksize requested'
c              STOP
C           ENDIF
C
C           CALL DNLASO(HESS, Q, NATOMS, OP, IOVECT, 3*MXATMS, NVAL, NFIG, NPERM,
C    *                  NMVAL, VAL, NMVEC, LVEC, NBLOCK, MAXOP, MAXJ, LWORK,
C    *                  IND, IERR)
C           WRITE(*,'(I3,A,100F20.10)') NPERM,' eigenvectors determined: ',(VAL(J1,1),J1=1,NPERM)
C           IF (IERR.NE.0) THEN
C              WRITE(*,'(A,I4)') 'Lanczos call completed with error code ',IERR
C           ELSE
C              WRITE(*,'(A,F20.10)') ' Smallest eigenvalue=',VAL(1,1)
C           ENDIF

C           EVALMIN=VAL(1,1)
C           DO J1=1,NOPT
C              VECS(J1)=LVEC(J1,1)
C           ENDDO
         ELSE
C jmc
C Put internals into COORDS here.  COORDS array should be unchanged on being passed through BEIG.
C Since we're passing the energy corresponding to these coords, I'm assuming that the unres c and 
C geometry arrays have already been updated correctly (in mylbfgs at the end of each 
C iteration and before this routine is called - 18/10/03 in fact at the start of this subroutine...)
            CALL geom_to_var(NINTS,COORDS(1:NINTS))
            CALL INTBEIG(ITER,COORDS,ENERGY,VECS,EVALMIN,NS,SOVER,PTEST,CONVERGED)
C
C  The following two routines are also obsolete.
C
C           CALL POWEIG(ITER,COORDS,ENERGY,VECS,EVALMIN)
C           CALL MCEIG(ITER,COORDS,ENERGY,VECS,EVALMIN)
         ENDIF
      ENDIF

      IF ((EVALMIN.LT.0.0D0).AND.FIXD) THEN
         IF (PTEST) WRITE(*,'(A)') ' Negative eigenvalue, changing to hybrid EF'
         FIXD=.FALSE.
         CALL VECNORM(VECS,NOPT)
         DO J1=1,NOPT
            ZWORK(J1,1)=VECS(J1)
         ENDDO
      ELSE IF (FIXD) THEN
C     IF (FIXD) THEN
         DO J1=1,NOPT
            ZWORK(J1,1)=FIXDIR(J1)
            VECS(J1)=FIXDIR(J1)
         ENDDO
      ELSE ! jmc this is the only option that works here...
         CALL VECNORM(VECS(1:NINTS),NINTS) ! jmc
         DUMMY=0.0D0
         DO J1=1,NINTS
            ZWORK(J1,1)=VECS(J1)
            DUMMY=DUMMY+VECS(J1)**2
         ENDDO
      ENDIF
C
C  Dump the smallest non-zero eigenvalue and eigenvector 
C  in file vector.dump, if required.
C
      IF (DUMPV) THEN
         IF (.NOT.ALLSTEPS) REWIND(44)
         WRITE(44,'(E20.10)') EVALMIN
         WRITE(44,'(3F20.10)') (ZWORK(J1,1),J1=1,NINTS)
      ENDIF
C
C  Take an eigenvector-following step uphill along the direction corresponding
C  to the smallest non-zero eigenvalue. Then do a line minimization along the
C  gradient vector with component along the uphill direction projected out. 
C  Do we need a new gradient vector or will the old one from the point before
C  stepping uphill do?
C

      FOB=0.0D0
      DO J1=1,NINTS
         FOB=FOB+ZWORK(J1,1)*GRAD(J1)
      ENDDO
      XFOB(1)=FOB
      IF (FIXD) GOTO 666
      IF (HINDEX.LE.1) THEN
         IF (EVALMIN.LT.0.0D0) THEN
C           PRINT*,'There is at least one negative eigenvalue'
            INEG=1
         ELSE
C           PRINT*,'There are no negative eigenvalues'
            INEG=0
         ENDIF
      ELSE
         INEG=0
         DO J1=1,HINDEX
            IF (DIAG(J1).LT.0.0D0) INEG=INEG+1
         ENDDO
         WRITE(*,'(A,I5,A)') ' There are at least ',INEG,' negative Hessian eigenvalues'
      ENDIF
C
C  Take a step away from a stationary point along the appropriate
C  Hessian eigenvector. This enables us to start from converged minima.
C  Distinguish the case where we want to take a very small step away
C  from a transition state from others where we want a big displacement
C  to get unstuck. 
C
      AWAY=.FALSE.
      IF ((RMS.LT.PUSHCUT).AND.(.NOT.FIXD)) THEN
         IF ((INEG.NE.HINDEX).AND.CONVERGED) THEN
            IF (HINDEX.LE.1) THEN
               IF ((ITER.EQ.1).AND.(INEG.EQ.0)) THEN
                  IF (PTEST) THEN
                     IF (IVEC.GE.0) PRINT*,'Stepping away from minimum along softest mode + direction'
                     IF (IVEC.LT.0) PRINT*,'Stepping away from minimum along softest mode - direction'
                  ENDIF
                  AWAY=.TRUE.
               ELSE IF (MOD(ITER-1,4).EQ.0) THEN
                  PRINT*,'Stepping away from solution of wrong index'
                  IF (PUSHOFF.EQ.0.0D0) THEN
                     FOB=STPMAX(1)
                  ELSE 
                     FOB=PUSHOFF
                  ENDIF
               ENDIF
            ELSE
               DO J1=1,HINDEX
                  IF (DIAG(J1).GT.0.0D0) THEN
                     WRITE(*,'(A,I6)') 'Stepping away from solution of wrong index along eigenvector ',J1
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
      IF (BFGSSTEP) THEN
         IF (INEG.EQ.0) THEN
            PRINT*,'****WARNING - BFGSSTEP set for a minimum'
         ELSE
            IF (PTEST) THEN
               IF (IVEC.GE.0) PRINT*,'Stepping away from saddle along softest mode + direction'
               IF (IVEC.LT.0) PRINT*,'Stepping away from saddle along softest mode - direction'
            ENDIF
            AWAY=.TRUE.
         ENDIF
      ENDIF
C
C  EF determination of steps
C
      IF (HINDEX.LE.1) THEN
         XP1=DABS(EVALMIN)/2.0D0
         IF (MASST) THEN
            SUM=0
            DO J2=1,NATOMS
               SUM=SUM+ATMASS(J2)
            ENDDO
            AVG=SUM/NATOMS
            XP2=1.0D0 + 4.0D0*(FOB/EVALMIN)**2/AVG
         ELSE
            XP2=1.0D0 + 4.0D0*(FOB/EVALMIN)**2
         ENDIF
         XP1=-XP1*(1.0D0+DSQRT(XP2))
         PSTEP=-FOB/XP1
   
         IF (AWAY) THEN
            IF (PUSHOFF.NE.0.0D0) THEN
               PSTEP=PUSHOFF
            ELSE
               PSTEP=STPMAX(1)/1.0D1
            ENDIF
            IF (IVEC.LT.0) PSTEP=-PSTEP
         ENDIF
c        WRITE(*,'(A,F20.10)') ' Unscaled step=',PSTEP
C
C  Scale according to step size in ev basis:
C
         STPMAG=ABS(PSTEP)
         IF (.NOT.AWAY) THEN
            SCALE=MIN(STPMAX(1)/MAX(DABS(PSTEP),1D-10),1.0D0)
            PSTEP=SCALE*PSTEP
         ENDIF
      ELSE
         DO I=HINDEX,1,-1
            STEP(I)=0.0D0
            NZERO=0
            IF (DABS(DIAG(NZERO+1)/MAX(DABS(DIAG(I)),1.0D-10)).GT.TEMP) TEMP=DABS(DIAG(NZERO+1)/MAX(DABS(DIAG(I)),1.0D-10))
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
            WRITE(*,'(A,I4,A,4X,F19.10)') ' Mode ',I,' will be searched uphill. Eigenvalue=',DIAG(I)
            LP=-LP
            STEP(I)=-XFOB(I)/LP
         ENDDO
         DO J1=1,NOPT
            SCRATCH(NOPT+J1)=STEP(J1)
         ENDDO
C
C  Convert the steps to the Cartesian rather than the EV basis
C
         DO J=1,NOPT
            SCRATCH(J+NOPT)=0.0D0
            DO I=1,HINDEX
               SCRATCH(J+NOPT)=SCRATCH(J+NOPT)+STEP(I)*ZWORK(J,I)
            ENDDO
         ENDDO
C
C  Scale according to step size in ev basis:
C
C jmc         CALL VSTAT(STEP(1),AV,NOPT,3*NATOMS)
C will almost never get to here in my runs, so probably not worth worrying about, but think should have the dimension of 
C the arrays as 1...
         CALL VSTAT(STEP(1),AV,1,1)
         STPMAG=MAX(AV(1),1D-10)
         DO J1=1,NOPT
            SCALE=MIN(STPMAX(J1)/MAX(DABS(STEP(J1)),1D-10),1.0D0)
            STEP(J1)=SCALE*STEP(J1)
         ENDDO
C jmc         CALL VSTAT(STEP(1),AV,NOPT,3*NATOMS)
C will almost never get to here in my runs, so probably not worth worrying about, but think should have the dimension of 
C the arrays as 1...
         CALL VSTAT(STEP(1),AV,1,1)
         SSTPMAG=MAX(AV(1),1D-10)
         SCALE=1.0D0
         E1=0.0D0
         DO J1=1,NOPT
            E1=E1+STEP(J1)*XFOB(J1)
         ENDDO
         E2=0.0D0
         DO J1=1,NOPT
            E2=E2+DIAG(J1)*STEP(J1)**2
         ENDDO
         E2=E2/2.0D0
         DELE=E1*SCALE+E2*SCALE**2
C jmc         CALL VADD(SCRATCH(1),SCRATCH(1),SCRATCH(NOPT+1),NOPT,1)
C will almost never get to here in my runs, so probably not worth worrying about, but think should have the dimension of 
C the arrays as 1...
         CALL VADD(SCRATCH(1),SCRATCH(1),SCRATCH(NINTS+1),1,1)

         IF (EFSTEPST.AND.(MOD(ITER-1,EFSTEPS).EQ.0)) THEN
            DO I=NZERO+1,HINDEX
                WRITE(*,360) I, STEP(I)
360             FORMAT(' Unscaled step for mode ',I3,'=',F20.10)
            ENDDO
         ENDIF
      ENDIF
C
C  Use MAXMAX until we have a negative eigenvalue. Be sure to
C  step uphill!
C
      LINE=.FALSE.
      IF (CONVERGED) THEN
         IF (((EVALMIN.GT.0.0D0).OR.FIXD).AND.(HINDEX.LE.1)) THEN
            PSTEP=FOB*MAXMAX/ABS(FOB)
         ELSE IF (HINDEX.GT.1) THEN
            DO J1=1,HINDEX
               IF (DIAG(J1).GT.0.0D0) STEP(J1)=XFOB(J1)*MAXMAX/ABS(XFOB(J1))
            ENDDO
         ENDIF
      ENDIF

      IF (HINDEX.LE.1) THEN
         E1=PSTEP*FOB
         E2=EVALMIN*PSTEP**2/2.0D0
         DELE=E1+E2
      ENDIF
C
C  Regenerate full Q vector. 
C
666   IF (.NOT.FIXD) THEN
         IF (HINDEX.LE.1) THEN
            DO J=1,NINTS
               COORDS(J)=COORDS(J)+PSTEP*ZWORK(J,1)
               INTSTEP(J)=PSTEP*ZWORK(J,1) ! jmc for use in path
            ENDDO
         ELSE
C
C  Unpack SCRATCH(1:NOPT) to regenerate full Q vector.
C  CSTEP contains the step in the Cartesian basis.
C
            DO J=1,NOPT
               COORDS(J)=SCRATCH(J)
               CSTEP(J)=SCRATCH(NOPT+J)
            ENDDO
         ENDIF
      ELSE
         DO J=1,NOPT
            COORDS(J)=COORDSN(J)
         ENDDO
      ENDIF
      KNOWE=.FALSE.
      KNOWG=.FALSE.
      KNOWH=.FALSE.

      IF (.NOT.LINE) THEN
         EOLD=ENERGY
         IF (.NOT.BFGSSTEP) THEN
C
C  Optimising the box lengths here would change the critical eigenvalue and
C  eigenvector.
C
C           IF (PV) THEN
C              CALL POTENTIAL(COORDS,ENERGY,GRAD,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
C              PVFLAG=.FALSE.
C              CALL PVOPT(COORDS,ENERGY,GRAD)
C           ENDIF
C
C jmc update Cartesians now
            CALL var_to_geom(NINTS,COORDS(1:NINTS))
!           TMPINT=COORDS(1:NINTS)
!           CALL var_to_geom(NINTS,TMPINT)
            CALL chainbuild
C jmc added this do loop to put cartesians into coords for dumpp
            DO J4=1,nres
               COORDS(6*(J4-1)+1)=c(1,J4)
               COORDS(6*(J4-1)+2)=c(2,J4)
               COORDS(6*(J4-1)+3)=c(3,J4)
               COORDS(6*(J4-1)+4)=c(1,J4+nres)
               COORDS(6*(J4-1)+5)=c(2,J4+nres)
               COORDS(6*(J4-1)+6)=c(3,J4+nres)
            END DO

            CALL POTENTIAL(COORDS,ENERGY,GRAD,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            CALL DUMPP(COORDS,ENERGY)
!           COORDS(1:NINTS)=TMPINT ! jmc put internals back into coords
         ENDIF
      ENDIF
      CALL FLUSH(6)

      FOBNEW=0.0D0
      DO J1=1,NINTS
         FOBNEW=FOBNEW+GRAD(J1)*ZWORK(J1,1)
      ENDDO
C
C  Only scale if we have a -ve eigenvalue
C
      IF ((EOLD.NE.0.0D0).AND.(EVALMIN.LT.0.0D0).AND.(.NOT.FIXD)) THEN
         EPER=MIN(DABS(1.0D0-(FOBNEW-FOB)/(PSTEP*EVALMIN)),DABS(1.0D0-(-FOBNEW-FOB)/(PSTEP*EVALMIN)))
C        WRITE(*,'(A,4F20.10)') 'FOB,FOBNEW,PSTEP,EVALMIN=',FOB,FOBNEW,PSTEP,EVALMIN
C        WRITE(*,'(A,3F20.10)') 'EPER,EPER1,EPER2=',
C    1          EPER,DABS(1.0D0-(FOBNEW-FOB)/(PSTEP*EVALMIN)),DABS(1.0D0-(-FOBNEW-FOB)/(PSTEP*EVALMIN))
         IF (EPER.GT.TRAD) THEN
            STPMAX(1)=MAX(STPMAX(1)/1.1D0,MINMAX)
C           WRITE(*,'(A,E12.4,A,E12.4)') ' Decreasing allowed EF step; trust radius=',TRAD,' calculated ratio=',EPER
         ELSE
            STPMAX(1)=MIN(STPMAX(1)*1.1D0,MAXMAX)
C           WRITE(*,'(A,E12.4,A,E12.4)') ' Increasing allowed EF step; trust radius=',TRAD,' calculated ratio=',EPER
         ENDIF
      ENDIF
C
C Summarize
C
C     IF (SUMMARYT.AND.(MOD(ITER-1,NSUMMARY).EQ.0)) THEN
         IF (PTEST) THEN
            IF (HINDEX.LE.1) THEN
               WRITE(*,30)
30             FORMAT(1X,79('-'))
               WRITE(*,40)
40             FORMAT(' Vector      Gradient        Secder       Step          Max step    Trust ratio')
               WRITE(*,30)
               WRITE(*,50) 1,FOB,EVALMIN,PSTEP,STPMAX(1),EPER
50             FORMAT(1X,I4,1X,E15.6,1X,E15.6,1X,E13.6,1X,E12.6,1X,E15.6)
               WRITE(*,30)
            ELSE
               WRITE(*,30)
               WRITE(*,40)
               WRITE(*,30)
               DO I=HINDEX,1,-1
                  WRITE(*,50) I,XFOB(I),DIAG(I),STEP(I),STPMAX(I),XRAT(I)
               ENDDO
               WRITE(*,30)
            ENDIF
         ELSE IF (.NOT.BFGSSTEP) THEN
C           WRITE(*,'(A,I6,A,F20.10,A,F15.7,A,F12.4,A,F12.4)') 
C    1          ' Cycle ',ITER,' E=',ENERGY,' RMS=',RMS,' eigenvalue=',EVALMIN,' overlap=',SOVER
         ENDIF
C     ENDIF
C
      IF (BFGSSTEP) THEN
C jmc put Cartesians back into COORDS
!        TMPINT=COORDS(1:NINTS)
!        CALL var_to_geom(NINTS,TMPINT)
         CALL var_to_geom(NINTS,COORDS(1:NINTS))
         CALL chainbuild
         DO J4=1,nres
            COORDS(6*(J4-1)+1)=c(1,J4)
            COORDS(6*(J4-1)+2)=c(2,J4)
            COORDS(6*(J4-1)+3)=c(3,J4)
            COORDS(6*(J4-1)+4)=c(1,J4+nres)
            COORDS(6*(J4-1)+5)=c(2,J4+nres)
            COORDS(6*(J4-1)+6)=c(3,J4+nres)
         END DO
         FIXIMAGE=.FALSE.
         FIXD=FIXDSAVE
         RETURN
      ENDIF
C
C  Tangent space minimization next.
C  Uphill direction is projected out of the step in mylbfgs
C  The next IF block allows for zero tangent space steps in the initial phase
C
      IF ((NBFGSMAX1.EQ.0).AND.((1.0D0-DABS(SOVER).GT.0.0001D0).OR.(INEG.EQ.0))) THEN
         FIXIMAGE=.FALSE.
         ITER=ITER+1
         FIXD=FIXDSAVE
         IF (ITER.GT.ITMAX) RETURN
         GOTO 90
      ENDIF

      IF (((1.0D0-DABS(SOVER).GT.0.0001D0).OR.(HINDEX.NE.INEG)).OR.(INEG.EQ.0).OR.(TWOENDS.AND.(ITER.EQ.1)).OR.(ITER.EQ.1)) THEN
          NBFGS=NBFGSMAX1
      ELSE
          NBFGS=NBFGSMAX2
      ENDIF

      RESET=.FALSE.
      IF (ITER.EQ.1) RESET=.TRUE.
      RMS2=RMS
C jmc put Cartesians back into COORDS
C jmc 13/7/04 now no need, as not putting internals back into COORDS above.
!     TMPINT=COORDS(1:NINTS)
!     CALL var_to_geom(NINTS,TMPINT)
c     CALL var_to_geom(NINTS,COORDS(1:NINTS))
c     CALL chainbuild
c     DO J4=1,nres
c        COORDS(6*(J4-1)+1)=c(1,J4)
c        COORDS(6*(J4-1)+2)=c(2,J4)
c        COORDS(6*(J4-1)+3)=c(3,J4)
c        COORDS(6*(J4-1)+4)=c(1,J4+nres)
c        COORDS(6*(J4-1)+5)=c(2,J4+nres)
c        COORDS(6*(J4-1)+6)=c(3,J4+nres)
c     END DO
      CALL MYLBFGS(NINTS,MUPDATE,COORDS,.FALSE.,GMAX,MFLAG,ENERGY,RMS2,EREAL,RMS,NBFGS,
     1             RESET,ITDONE,PTEST,GRAD,.FALSE.,.FALSE.)

      IF (REPELTST) CALL REPELSP(COORDS,RMS,INEG,MFLAG)

      IF (MFLAG) THEN
         IF (((RMS.GT.CONVR).OR.(INEG.EQ.0)).OR.(STPMAG.GT.CONVU)) MFLAG=.FALSE.
         IF (PTEST) WRITE(*,'(A,F15.7,A,F15.7,A,F15.7)') ' Total gradient=',RMS,' subspace gradient=',RMS2,' EF step length=',STPMAG
         IF (MFLAG) THEN
C jmc testing!!
c     CALL MAKENUMINTHESS(MYHESS,NINTS,NATOMS)
c     ABSTOL=DLAMCH('Safe  minimum')
c     CALL DSYEVR('V','I','U',NINTS,HESS,3*NATOMS,0.0D0,1.0D0,1,HINDEX,ABSTOL,NFOUND,DIAG
c     CALL DSYEVR('V','I','U',NINTS,MYHESS,3*NATOMS,0.0D0,1.0D0,1,50,ABSTOL,NFOUND,DIAG
c    1,ZMAT,3*NATOMS,ISUPPZ,WORK,
c    2                        LWORK, IWORK, ILWORK, INFO )
c     IF (INFO.NE.0) PRINT*,'WARNING - INFO=',INFO,' in DSYEVR'
c     DO J4=1,50
c     PRINT *,'EVALMIN ',DIAG(J4)
c     END DO
            FIXD=FIXDSAVE
            RETURN
         ENDIF
      ENDIF
      ITER=ITER+1
      IF (ITER.GT.ITMAX) THEN
         FIXD=FIXDSAVE
         RETURN
      ENDIF
C     IF ((AMBER).AND.(MOVIE)) CALL AMOVIEDUMP(FRAME)
      EOLD=ENERGY
      DO J1=1,NINTS
         PCSTEP(J1)=CSTEP(J1)
         XPSTEP(J1)=STEP(J1)
         PFOB(J1)=XFOB(J1)
      ENDDO
      GOTO 90

      RETURN
      END
