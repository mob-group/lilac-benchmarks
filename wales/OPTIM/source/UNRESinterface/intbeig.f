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
C  This subroutine performs LBFGS minimization to estimate the
C  smallest eigenvalue and eigenvector of the Hessian without second
C  derivatives...
C  ... in internal coordinates (jmc March 03).
C
C***********************************************************************
C
      SUBROUTINE INTBEIG(ITER,COORDS,ENERGY,VECS,EVALMIN,NS,SOVER,PTEST,CONVERGED)

      USE COMMONS
      USE KEY
      USE MODNEB
      USE MODTWOEND
      IMPLICIT NONE

      INTEGER J1, ISEED, NS, NITS
      DOUBLE PRECISION ENERGY,COORDS(3*NATOMS),VEC(3*NATOMS),VECS(3*NATOMS),DPRAND,SOVER,EVALMIN,EPS,FRET
      DOUBLE PRECISION DIAG(3*NATOMS),WORK(3*NATOMS*(2*XMUPDATE+1)+2*XMUPDATE)
      PARAMETER (EPS=1.D-6)
      LOGICAL PTEST, CONVERGED
      INTEGER ITER
      COMMON /IS/ ISEED
      SAVE

!     VEC=0.0D0 ! jmc testing

      IF ((ITER.EQ.1).AND.(.NOT.READV).AND.(.NOT.TWOENDS).AND.(.NOT.TTDONE).AND.(.NOT.(NEWCONNECTT.OR.NEWNEBT))) THEN
         DO J1=1,NINTS
C jmc VEC(J1) contains a random number between 1 and -1
            VEC(J1)=2*(DPRAND()-0.5D0)
C           VEC(J1)=0.0D0
         ENDDO
C        VEC(1)=1.0D0
      ELSE
         DO J1=1,NINTS
C jmc VECS contains previous vector.
            VEC(J1)=VECS(J1)
         ENDDO
      ENDIF

C jmc can't use freeze...
      IF (FREEZE) THEN
         PRINT *,"** WARNING: Trying to use FREEZE with INTMIN; please don't!" ! jmc
         STOP
         DO J1=1,NATOMS
            IF (FROZEN(J1)) THEN
               VEC(3*(J1-1)+1)=0.0D0
               VEC(3*(J1-1)+2)=0.0D0
               VEC(3*(J1-1)+3)=0.0D0
            ENDIF
         ENDDO
      ENDIF

C     CALL VECNORM(VEC,NINTS)

      CALL INTXMYLBFGS(NINTS,XMUPDATE,VEC,.FALSE.,DIAG,CEIG,WORK,ENERGY,COORDS,NITS,NEVS,FRET,PTEST,CONVERGED)
      CALL VECNORM(VEC,NINTS)

      EVALMIN=FRET
      NS=NITS
      SOVER=0.0D0
      DO J1=1,NINTS
         SOVER=SOVER+VECS(J1)*VEC(J1)
         VECS(J1)=VEC(J1)
      ENDDO
      IF (PTEST) WRITE(*,'(A,F15.7)') ' Overlap with previous vector=',SOVER

      FIXIMAGE=.FALSE.

      RETURN
      END
