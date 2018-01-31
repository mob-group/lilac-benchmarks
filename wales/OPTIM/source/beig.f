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
***********************************************************************
C
C  This subroutine performs LBFGS minimization to estimate the
C  smallest eigenvalue and eigenvector of the Hessian without second
C  derivatives.  It calls xmylbfgs to perform the minimization.  xmylbfgs
C  minimizes secdiag as the objective function
C
C  NOTE: This subroutine uses the 'SAVE' statement.  All variables
C  will save their value between calls
C
C***********************************************************************
C
      SUBROUTINE BEIG(ITER,COORDS,ENERGY,VECS,EVALMIN,NS,SOVER,PTEST,CONVERGED)
      USE COMMONS
      USE KEY
      USE MODNEB
      USE MODTWOEND
      use porfuncs
      USE GENRIGID ! hk286
      IMPLICIT NONE

      INTEGER, INTENT(IN) ::             ITER ! the number of iterations of Hybrid eigenvector following we have already done
      INTEGER, INTENT(OUT) ::            NS ! the number of iterations spent in LBFGS
      DOUBLE PRECISION, INTENT(INOUT) :: VECS(3*NATOMS) ! The eigenvector
      DOUBLE PRECISION, INTENT(OUT) ::   EVALMIN ! the eigenvalue
      DOUBLE PRECISION, INTENT(OUT) ::   SOVER ! the overlap with the previous eigenvector
!     DOUBLE PRECISION, INTENT(IN) ::    COORDS(3*NATOMS) ! the coordinates at which the eigenvector is being calculated
      DOUBLE PRECISION COORDS(3*NATOMS) ! the coordinates at which the eigenvector is being calculated
      DOUBLE PRECISION, INTENT(INOUT) :: ENERGY ! the energy at COORDS
      LOGICAL, INTENT(IN) ::             PTEST ! Flag for printing status info
      INTEGER J1, ISEED, NITS
      DOUBLE PRECISION VEC(3*NATOMS),DPRAND,EPS,FRET
      DOUBLE PRECISION DIAG(3*NATOMS),WORK(3*NATOMS*(2*XMUPDATE+1)+2*XMUPDATE)
      PARAMETER (EPS=1.D-6)
      LOGICAL CONVERGED
      COMMON /IS/ ISEED
      SAVE

      IF (RIGIDINIT) THEN ! hk286
         CALL BEIGLOCALRIGID(ITER, COORDS,ENERGY,VECS,EVALMIN,NS,SOVER,PTEST,CONVERGED)
         RETURN
      ENDIF

      ! initialize the eigenvector
      IF ((ITER.EQ.1).AND.(.NOT.READV).AND.(.NOT.TWOENDS).AND.(.NOT.TTDONE).AND.(.NOT.(NEWCONNECTT.OR.NEWNEBT.OR.GROWSTRINGT))) THEN
         DO J1=1,NOPT
            VEC(J1)=2*(DPRAND()-0.5D0)
C           VEC(J1)=0.0D0
         ENDDO
C        VEC(1)=1.0D0
      ELSE
         DO J1=1,NOPT
            VEC(J1)=VECS(J1)
         ENDDO
      ENDIF
      IF (FREEZE) THEN
         DO J1=1,NATOMS
            IF (FROZEN(J1)) THEN
               VEC(3*(J1-1)+1)=0.0D0
               VEC(3*(J1-1)+2)=0.0D0
               VEC(3*(J1-1)+3)=0.0D0
            ENDIF
         ENDDO
      ENDIF

C     CALL VECNORM(VEC,NOPT)
      IF (CPMD) CALL SYSTEM(' mv RESTART.1 RESTART.1.save ')
      IF (QCHEM) CALL SYSTEM(' cp ' // SYS(1:LSYS) // '  ' // SYS(1:LSYS) // '.old ')
      IF (VASP) THEN
         CALL SYSTEM(' cp POSCAR POSCAR.save ')
         CALL SYSTEM(' cp POTCAR POTCAR.save ')
         CALL SYSTEM(' cp INCAR INCAR.save ')
         CALL SYSTEM(' cp KPOINTS KPOINTS.save ')
      ENDIF
      IF (CASTEP) CALL SYSTEM(' mv ' // SYS(1:LSYS) // '.check ' // SYS(1:LSYS) // '.check.save ')
      IF (CASTEP) CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.cell ' // SYS(1:LSYS) // '.cell.save ')
      IF (ONETEP) CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.dat ' // SYS(1:LSYS) // '.dat.save ')
      IF (CP2K) CALL SYSTEM(' cp ' // SYS(1:LSYS) // '.inp ' // SYS(1:LSYS) // '.inp.save ') 
      CALL XMYLBFGS(NOPT,XMUPDATE,VEC,.FALSE.,DIAG,CEIG,WORK,ENERGY,COORDS,NITS,NEVS,FRET,PTEST,CONVERGED)
      CALL VECNORM(VEC,NOPT) 

      EVALMIN=FRET
      NS=NITS
      SOVER=0.0D0
      DO J1=1,NOPT
         SOVER=SOVER+VECS(J1)*VEC(J1)
         VECS(J1)=VEC(J1)
      ENDDO
      IF (PTEST) WRITE(*,'(A,F15.7)') 'beig> Overlap with previous vector=',SOVER

      IF (CPMD) CALL SYSTEM(' rm RESTART.1.plus ; rm RESTART.1.minus ; mv RESTART.1.save RESTART.1 ')
C     IF (CASTEP) CALL SYSTEM(' rm ' // SYS(1:LSYS) // '.check.plus  ; rm ' 
C    1                               // SYS(1:LSYS) // '.check.minus ; mv ' 
C    2                               // SYS(1:LSYS) // '.check.save ' // SYS(1:LSYS) // '.check ')

      FIXIMAGE=.FALSE.

      RETURN
      END
