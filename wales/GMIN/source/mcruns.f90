!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
SUBROUTINE MCRUNS
   USE COMMONS, ONLY: NATOMS, PTMC, BSPT, ONE_ATOM_TAKESTEP, CHECKDT, COORDS, &
                      NRUNS, MCSTEPS, TFAC, RELAXFQ, NSAVE, GBHT, AMBERMUTATIONT, &
                      NMUTATION
   USE GENRIGID, ONLY: RIGIDINIT
   USE QMODULE, ONLY: QMINP
   USE PREC, ONLY: INT32, REAL64
   USE AMBER12_MUTATIONS, ONLY: AMBERMUT_CURR_LOWEST,FINISH_AMBERMUT
   USE AMBER12_INTERFACE_MOD, ONLY: AMBER12_FINISH
   IMPLICIT NONE
! Variables
   LOGICAL                    :: LOPEN
   REAL(REAL64), ALLOCATABLE  :: QMINPTMP(:, :)
   INTEGER(INT32)             :: I

   IF (PTMC .OR. BSPT) THEN
      IF (ONE_ATOM_TAKESTEP) THEN
         CALL PTMC_ONE_ATOM()
      ELSE
         CALL PTBASINSAMPLING()
      END IF
      RETURN
   END IF

   IF (GBHT) THEN ! ds656> Generalised basin-hopping run.
      CALL MC_GBH(MCSTEPS(1))
      CALL FINALQ()
      CALL FINALIO()
      RETURN
   ENDIF

!   IF (MFETT) THEN
!      CALL MFET()
!      CALL FINALIO()
!      RETURN
!   ENDIF

   INQUIRE(UNIT=1, OPENED=LOPEN)
   IF (LOPEN) THEN
      WRITE(*,'(A,I2,A)') 'mcruns> ERROR *** Unit ', 1, ' is not free '
      STOP
   ENDIF

! jwrm2> If checking derivatives, call CHECKD, which then exits
   IF (CHECKDT) THEN
      CALL CHECKD(COORDS(:, :))
   END IF

!
!  NRUNS > 1 is an obsolete option! DJW
!

   DO I = 1, NRUNS
      CALL MC(MCSTEPS(I), TFAC(I))
   END DO

!     DO J1=1,NPAR
!        CLOSE(DUMPVUNIT(J1))
!        CLOSE(DUMPXYZUNIT(J1))
!     ENDDO
!     DUMPT=.FALSE.

   IF (AMBERMUTATIONT) THEN
      CALL FINISH_AMBERMUT()
      CALL AMBER12_FINISH()
      RETURN
   ENDIF

   IF (RELAXFQ) THEN
      ALLOCATE(QMINPTMP(NSAVE, 3*NATOMS))
      QMINPTMP(:,:) = QMINP(:,:)
      RELAXFQ = .FALSE.
      RIGIDINIT = .TRUE.
      CALL FINALQ()
      CALL FINALIO()
      QMINP(:,:) = QMINPTMP(:,:)
      DEALLOCATE(QMINPTMP)
      RELAXFQ = .TRUE.
      RIGIDINIT = .FALSE.
   END IF
   CALL FINALQ()
   CALL FINALIO()
   RETURN

END SUBROUTINE MCRUNS
