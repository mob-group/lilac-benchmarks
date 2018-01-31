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
      SUBROUTINE TWOEND(ENERGY,GRAD,VECS,Q)
      use porfuncs
      USE COMMONS
      USE KEY
      USE MODTWOEND
      USE ZWK
      USE MODCHARMM
      USE MODHESS
      IMPLICIT NONE
      INTEGER J1, ITER, NS, J3, ITDONE
      DOUBLE PRECISION DIST, ENERGY, RMS, VECS(3*NATOMS), EVALMIN, EVALMAX, 
     1                 GRAD(3*NATOMS), SOVER, RMS2, EREAL, PROJP,
     2                 EP, STARTP(3*NATOMS), DELTAP, DELTA, RVEC(3*NATOMS), DUMMY1, PROJ, Q(3*NATOMS), VECL(3*NATOMS)
C
C  Assign enough memory to WORK for a blocksize of 32 to be possible.
C  This is for DSYEVR.
C
      INTEGER ILWORK, LWORK, NFOUND, ISUPPZ(2*3*NATOMS)
      INTEGER IWORK(33*3*NATOMS), INFO
      DOUBLE PRECISION WORK(33*3*NATOMS),  ABSTOL, DIAG(3*NATOMS), ZMAT(3*NATOMS,3*NATOMS), DLAMCH
      LOGICAL PVFLAG, MFLAG, RESET, SWITCH,CONVERGED           !  is a full ZMAT really needed?
      COMMON /PVF/ PVFLAG
C     COMMON /WORK/ ZMAT

      LWORK=33*3*NATOMS
      ILWORK=33*3*NATOMS
C
C  Bulk principal image code assumes that PARAM1, PARAM2 and PARAM3 are the box lengths.
C
      OPEN(UNIT=87,FILE='candidate',STATUS='UNKNOWN')

      DO J1=1,NATOMS
         J3=3*(J1-1)
         VECS(J3+1)=FIN(J3+1)-Q(J3+1)
         VECS(J3+2)=FIN(J3+2)-Q(J3+2)
         VECS(J3+3)=FIN(J3+3)-Q(J3+3)
         LANV(J3+1)=0
         LANV(J3+2)=0
         LANV(J3+3)=0
         IF (BULKT) THEN
            LANV(J3+1)=NINT(VECS(J3+1)/PARAM1)
            LANV(J3+2)=NINT(VECS(J3+2)/PARAM2)
            LANV(J3+3)=NINT(VECS(J3+3)/PARAM3)
         ENDIF
         VECS(J3+1)=VECS(J3+1)-PARAM1*LANV(J3+1)
         VECS(J3+2)=VECS(J3+2)-PARAM2*LANV(J3+2)
         VECS(J3+3)=VECS(J3+3)-PARAM3*LANV(J3+3)
         START(J3+1)=Q(J3+1)
         START(J3+2)=Q(J3+2)
         START(J3+3)=Q(J3+3)
      ENDDO

C     DO J1=1,NOPT
C        LANV(J1)=0
C        IF (BULKT) LANV(J1)=NINT((FIN(J1)-Q(J1))/PARAM3)
C        VECS(J1)=FIN(J1)-Q(J1)-PARAM3*LANV(J1)
C        START(J1)=Q(J1)
C     ENDDO

      CALL ORTHOGOPT(VECS,Q,.TRUE.)
      ITER=1
C
C  Initial distance
C
      DIST=0.0D0
      DO J1=1,NOPT
         DIST=DIST+(FIN(J1)-START(J1)-PARAM3*LANV(J1))**2
C        WRITE(*,'(A,I6,2F20.10,I4,F20.10)') 'J1,FIN,START,LANV,DIST=',J1,FIN(J1),START(J1),LANV(J1),DIST
      ENDDO
      DIST=SQRT(DIST)
      EP=0.0D0
      PROJP=0.0D0
      DELTAP=0.0D0
      DO J1=1,NOPT
         STARTP(J1)=START(J1)
      ENDDO
C
C  Minimisation subject to force constant FORCE
C
10    SWITCH=.FALSE.
      DIST=0.0D0
      DO J1=1,NOPT
         DIST=DIST+(FIN(J1)-START(J1)-PARAM3*LANV(J1))**2
      ENDDO
      DIST=SQRT(DIST)

      NUP=1
      DO J1=1,NOPT
         ZWORK(J1,1)=ZMAT(J1,1)
      ENDDO
C
C  Calculate required force constant
C
      DUMMY1=0.0D0
      DO J1=1,NOPT
         RVEC(J1)=FIN(J1)-START(J1)
         DUMMY1=DUMMY1+RVEC(J1)**2
      ENDDO
      DUMMY1=1.0D0/SQRT(DUMMY1)
C     WRITE(*,'(A,F20.10)') 'DIST=',1.0D0/DUMMY1
      DO J1=1,NOPT
         RVEC(J1)=RVEC(J1)*DUMMY1
      ENDDO
      DUMMY1=0.0D0
      DO J1=1,NOPT
        DUMMY1=DUMMY1+RVEC(J1)*GRAD(J1)
      ENDDO
      PROJ=DUMMY1
C     PRINT*,'Projection of gradient=',DUMMY1
C     FORCE=-DUMMY1-MAX(FINC*ABS(DUMMY1),FINC)
      FORCE=-DUMMY1-FINC-ABS(DUMMY1)*FINC
      IF (DUMMY1.LT.-FINC-ABS(DUMMY1)*FINC) FORCE=0.0D0
      RESET=.FALSE.
      IF (ITER.EQ.1) RESET=.TRUE.
      IF (UNRST.OR.CHRMMT) THEN
         PRINT*,'ERROR - twoend not compatible with UNRES or CHARMM'
         STOP
      ENDIF
      CALL MYLBFGS(NOPT,MUPDATE,START,.FALSE.,RMSTWO,MFLAG,ENERGY,RMS2,EREAL,RMS,NTWO,RESET,
     1             ITDONE,.FALSE.,GRAD,.FALSE.,.FALSE.)
      WRITE(23,'(3F20.10)') (START(J1),J1=1,NOPT)
      IF (NOHESS) THEN
         CALL BEIG(ITER,START,ENERGY,VECS,EVALMIN,NS,SOVER,.TRUE.,CONVERGED)
      ELSE 
         CALL POTENTIAL(START,ENERGY,GRAD,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
         IF (NOIT) THEN
            IF (.NOT.VARIABLES) CALL SHIFTH(START,.TRUE.,NOPT,NATOMS,ATMASS)

            ABSTOL=DLAMCH('Safe  minimum')
            CALL DSYEVR('V','I','U',NOPT,HESS,SIZE(HESS,1),0.0D0,1.0D0,1,1,ABSTOL,NFOUND,DIAG,ZMAT,3*NATOMS,ISUPPZ,WORK,
     1                     LWORK, IWORK, ILWORK, INFO )
            IF (INFO.NE.0) PRINT*,'WARNING - INFO=',INFO,' in DSYEVR'
C           PRINT*,'Optimal and actual values of LWORK=',WORK(1),LWORK
C           PRINT*,'Optimal and actual values of ILWORK=',IWORK(1),ILWORK
            EVALMIN=DIAG(1)

C           WRITE(*,'(A,F20.10)') ' Smallest eigenvalue=',EVALMIN
            DO J1=1,NOPT
               VECS(J1)=ZMAT(J1,1)
            ENDDO
C           PRINT*,'Eigenvalues:'
C           WRITE(*,'(6F15.5)') (DIAG(J1),J1=1,6)
         ELSE
            CALL ITEIG(ITER,START,VECS,EVALMIN,EVALMAX,NS,SOVER,.FALSE.,VECL,CONVERGED)
         ENDIF
      ENDIF
      DELTA=ENERGY-EP
      DELTA=PROJ-PROJP
      WRITE(*,'(I4,A,F15.7,A,F12.5,A,F15.5,A,F15.5,A,F12.5,A,F12.5,A,F12.5)') 
     1  ITER,' k=',FORCE,' D=',DIST,' E=',ENERGY,' proj=',PROJ,' delta=',DELTA,' evalue=',EVALMIN
      IF ((DELTA.LT.TWOEVAL).AND.(ITER.GT.2)) THEN
         SWITCH=.TRUE.
      ELSE
         DO J1=1,NOPT
            STARTP(J1)=START(J1)
         ENDDO
      ENDIF
      EP=ENERGY
      PROJP=PROJ
      DELTAP=DELTA
      WRITE(4,'(3F20.10)') (START(J1),J1=1,NOPT)
      IF (DIST.LT.0.05) THEN
         PRINT*,'final geometry reached in twoends'
         CLOSE(87)
         STOP
         RETURN
      ENDIF
      IF (DUMPV) THEN
         IF (.NOT.ALLSTEPS) REWIND(44)
         WRITE(44,'(E20.10)') EVALMIN
         WRITE(44,'(3F20.10)') (VECS(J1),J1=1,NOPT)
      ENDIF
      IF (SWITCH) THEN
         WRITE(*,'(A)') ' Candidate for transition state search located'
C        DO J1=1,NOPT
C           Q(J1)=START(J1)
C        ENDDO
C        TWOENDS=.FALSE.
         WRITE(87,'(3F20.10)') (START(J1),J1=1,NOPT)
C        RETURN
      ENDIF
      ITER=ITER+1
      IF (ITER.GT.NTWOITER) THEN
         PRINT*,'Maximum steps exceeded in double ended search'
         STOP
      ENDIF

      IF (BULKT) THEN
         DO J1=1,NATOMS
            J3=3*(J1-1)
            LANV(J3+1)=NINT((FIN(J3+1)-Q(J3+1))/PARAM1)
            LANV(J3+2)=NINT((FIN(J3+2)-Q(J3+2))/PARAM2)
            LANV(J3+3)=NINT((FIN(J3+3)-Q(J3+3))/PARAM3)
         ENDDO
      ENDIF
      GOTO 10

      CLOSE(87)
      RETURN
      END
