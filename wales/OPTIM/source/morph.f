C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2005 David J. Wales
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
C  This subroutine is designed to morph an initial geometry into a final 
C  endpoint avoiding clashes. 
C
C***********************************************************************
C
      SUBROUTINE MORPH(ITMAX,COORDS,QFINISH,ENERGY,GRAD,MFLAG,RMS,ITER,POTCALL)
      USE COMMONS
      USE KEY
      USE ZWK
      USE MODTWOEND
      USE MODCHARMM
      USE PORFUNCS
      IMPLICIT NONE

      INTEGER J1, J, ITER, NBFGS, ITMAX, ITDONE, ISTAT, NREPEL, NPU
      INTEGER ::  REPMAX=10
      DOUBLE PRECISION GRAD(3*NATOMS),ENERGY,COORDS(3*NATOMS),RMS,PSTEP,DISTS,DISTF,
     &                 PROJF,PROJS,DUM,QSTART(3*NATOMS),QFINISH(1:3*NATOMS),EOLD,RMS2,EREAL,SVEC(3*NATOMS),
     &                 DSAVE(ITMAX),MNBFGSMAX1SAVE,MNBFGSMAX2SAVE,STARTTIME,FINISHTIME,RMAT(3,3),DIST2, GMAXSAVE
      LOGICAL MFLAG, PVFLAG, RESET, POTCALL, PTEST, PUSH, PULL, FINISHED
      COMMON /PVF/ PVFLAG
      CHARACTER(LEN=5) ZSYMSAVE
      COMMON /SYS/ ZSYMSAVE
      LOGICAL KNOWE, KNOWG, KNOWH
      COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
      COMMON /MORPHDATA/ PUSH, PULL

      IF (DEBUG) OPEN(UNIT=995,FILE='morph.xyz',STATUS='UNKNOWN')
      CALL MYCPU_TIME(STARTTIME,.TRUE.)
      IF (ZSYM(1)(1:1).EQ.'W') THEN
         PRINT*,'MORPH procedures have not been programmed for TIP potentials'
         STOP
      ENDIF
      PRINT '(A,F8.2,A,F12.2,A,F10.2)',' morph> Interpolation using maximum step ',MORPHMXSTP,' energy maximum ',
     &                                   MORPHEMAX,' maximum rise ',MORPHERISE
      PULL=.TRUE.
      PUSH=.FALSE.
      NPU=0
      FINISHED=.FALSE.
      MNBFGSMAX1SAVE=MNBFGSMAX1
      MNBFGSMAX2SAVE=MNBFGSMAX2

      OPEN(UNIT=44,FILE='morph.points',STATUS='UNKNOWN')

      NUP=1 ! mylbfgs uses this for projection of the search direction
      QSTART(1:3*NATOMS)=COORDS(1:3*NATOMS)
      WRITE(*,*) "sn402: Changed MORPH so that it calls ALIGN_DECIDE instead of MINPERMDIST"
      WRITE(*,*) "I haven't tested this change and am not certain whether it's sensible." 
      WRITE(*,*) "Please set FASTOVERLAP and check carefully that this part of the code is working as you expect."
      WRITE(*,*) "Then remove these messages!"
      CALL ALIGN_DECIDE(QSTART,QFINISH,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DISTF,DIST2,RIGIDBODY,RMAT)
      IF (DEBUG) PRINT '(A,F20.10)',' morph> initial distance between start and finish=',DISTF
!     QSTART(1:3*NATOMS)=COORDS(1:3*NATOMS)
      NREPEL=0
      PRINTPTS=.TRUE.
      ITER=1
      PTEST=.FALSE.
      IF (DEBUG) PTEST=.TRUE.
90    CONTINUE
      IF (POTCALL) THEN
         IF (PV) THEN
            IF (.NOT.KNOWE) CALL POTENTIAL(COORDS,ENERGY,GRAD,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
            PVFLAG=.FALSE.
            CALL PVOPT(COORDS,ENERGY,GRAD)
         ENDIF
         IF (.NOT.KNOWG) CALL POTENTIAL(COORDS,ENERGY,GRAD,.TRUE.,.FALSE.,RMS,PTEST,.FALSE.)
      ENDIF
!     CALL DUMPP(COORDS,ENERGY)
      IF (DEBUG) WRITE(995,'(I6)') NATOMS
      IF (DEBUG) WRITE(995,'(G20.10)') ENERGY
      IF (DEBUG) WRITE(995,'(A3,3G20.10)') ('LA ',COORDS(3*(J1-1)+1),COORDS(3*(J1-1)+2),COORDS(3*(J1-1)+3),J1=1,NATOMS)

      CALL NEWMINDIST(COORDS,QFINISH,NATOMS,DISTF,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
      CALL NEWMINDIST(COORDS,QSTART,NATOMS,DISTS,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)

      PROJF=0.0D0
      PROJS=0.0D0
      DO J1=1,NOPT
         IF (DISTF.GT.0.0D0) PROJF=PROJF+GRAD(J1)*(QFINISH(J1)-COORDS(J1))/DISTF
         IF (DISTS.GT.0.0D0) PROJS=PROJS+GRAD(J1)*(COORDS(J1)-QSTART(J1))/DISTS
      ENDDO
      IF (PULL.AND.(PROJF.GT.0.0D0).AND.(ENERGY.GT.MORPHEMAX)) THEN
         PUSH=.TRUE.
         PULL=.FALSE.
         NPU=0
      ELSEIF (PUSH.AND.(PROJS.GT.0.0D0).AND.(ENERGY.GT.MORPHEMAX)) THEN
         PUSH=.FALSE.
         PULL=.TRUE.
         NPU=0
      ENDIF
      IF (NPU.GT.100) THEN
         IF (ABS(DSAVE(NPU-100)-DISTF).LT.0.1D0) THEN
            PRINT '(A,2G20.10)',' morph> Attempting to unstick: distances=',DSAVE(NPU-100),DISTF
            NPU=0
            MNBFGSMAX1=MAX(1,MNBFGSMAX1-1)
            MNBFGSMAX2=MAX(1,MNBFGSMAX2-1)
         ENDIF
      ENDIF
      IF (PULL) THEN
         PSTEP=MIN(DISTF,MORPHMXSTP)
         DO J=1,NOPT
           IF (DISTF.NE.0.0D0) SVEC(J)=(QFINISH(J)-COORDS(J))/DISTF
         ENDDO
      ELSEIF (PUSH) THEN
         PSTEP=MIN(DISTF,MORPHMXSTP)
         DO J=1,NOPT
           IF (DISTS.NE.0.0D0) SVEC(J)=(COORDS(J)-QSTART(J))/DISTS
         ENDDO
      ENDIF
C
C  Regenerate full Q vector. 
C
      DO J=1,NOPT
         COORDS(J)=COORDS(J)+PSTEP*SVEC(J)
      ENDDO

      KNOWE=.FALSE.
      KNOWG=.FALSE.
      KNOWH=.FALSE.
      EOLD=ENERGY
77    CALL POTENTIAL(COORDS,ENERGY,GRAD,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
      IF (ENERGY-EOLD.GT.MORPHERISE) THEN
         IF (PTEST) PRINT '(2(A,G20.10))',' morph> Energy rose by ',ENERGY-EOLD,' reducing step size to ',PSTEP/2.0D0
         DO J=1,NOPT
           COORDS(J)=COORDS(J)-PSTEP*SVEC(J)/2.0D0
         ENDDO
         PSTEP=PSTEP/2.0D0
         KNOWE=.FALSE.
         KNOWG=.FALSE.
         KNOWH=.FALSE.
         GOTO 77
      ENDIF
      CALL FLUSH(6)
C
C Summarize
C
      IF (PUSH.AND.DEBUG) PRINT '(A,I8,A,F15.6,6(A,F8.2))', 'PUSH ',ITER,' E=',ENERGY,' Step=',PSTEP,
     &     ' DS=',DISTS,' DF=',DISTF,' RMS=',RMS,' PF=',PROJF,' PS=',PROJS
      IF (PULL.AND.DEBUG) PRINT '(A,I8,A,F15.6,6(A,F8.2))', 'PULL ',ITER,' E=',ENERGY,' Step=',PSTEP,
     &     ' DS=',DISTS,' DF=',DISTF,' RMS=',RMS,' PF=',PROJF,' PS=',PROJS
      NPU=NPU+1
      DSAVE(NPU)=DISTF
C
C  Tangent space minimization next.
C  Step direction is projected out of the step in mylbfgs
C  The next IF block allows for zero tangent space steps in the initial phase
C
      IF (DISTF.GT.CONVU) THEN
          NBFGS=MNBFGSMAX1
      ELSE
          NBFGS=MNBFGSMAX2
      ENDIF

      RESET=.FALSE.
      RESET=.TRUE.
      RMS2=RMS ! needs to be set because MYLBFGS tests it for convergence !
      IF (ITER.EQ.1) RESET=.TRUE.
      ZWORK(1:3*NATOMS,1)=SVEC(1:3*NATOMS) ! used for the projection direction in mylbfgs
      GMAXSAVE=GMAX; GMAX=CONVR ! mylbfgs now uses GMAX so that we can change this parameter via changep
      IF (CHRMMT.AND.INTMINT) THEN
         CALL MYLBFGS(NINTS,MUPDATE,COORDS,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,NBFGS,
     1                RESET,ITDONE,PTEST,GRAD,.TRUE.,.TRUE.)
      ELSE
         CALL MYLBFGS(NOPT,MUPDATE,COORDS,.FALSE.,MFLAG,ENERGY,RMS2,EREAL,RMS,NBFGS,
     1                RESET,ITDONE,PTEST,GRAD,.TRUE.,.TRUE.)
      ENDIF
      GMAX=GMAXSAVE
      WRITE(44,'(3G25.15)') COORDS(1:NOPT)
      IF (MFLAG) THEN
         IF (((RMS.GT.CONVR).OR.(PSTEP.GT.CONVU))) MFLAG=.FALSE.
         IF (MFLAG) THEN
            CLOSE(44)
            IF (DEBUG) CLOSE(995)
            MNBFGSMAX1=MNBFGSMAX1SAVE
            MNBFGSMAX2=MNBFGSMAX2SAVE
            FINISHED=.TRUE.
            CALL MYCPU_TIME(FINISHTIME,.FALSE.)
            PRINT '(A,I8,A,F12.2)',' morph> Interpolation succeeded in ',ITER,' steps, time=',FINISHTIME-STARTTIME
            RETURN
         ENDIF
      ENDIF
      ITER=ITER+1
      IF (ITER.GT.ITMAX) THEN
         CLOSE(44)
         IF (DEBUG) CLOSE(995)
         MNBFGSMAX1=MNBFGSMAX1SAVE
         MNBFGSMAX2=MNBFGSMAX2SAVE
         FINISHED=.FALSE.
         CALL MYCPU_TIME(FINISHTIME,.FALSE.)
         PRINT '(A,I8,A,F12.2)',' morph> Interpolation failed to converge in ',ITER,' steps, time taken=',FINISHTIME-STARTTIME
!        STOP
         RETURN
      ENDIF
      EOLD=ENERGY
      GOTO 90

      RETURN
      END
