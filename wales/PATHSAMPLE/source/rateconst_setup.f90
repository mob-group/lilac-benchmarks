!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999-2009 David J. Wales
!   This file is part of PATHSAMPLE.
!
!   PATHSAMPLE is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   PATHSAMPLE is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! Taking the common code to deal with dodgy transition states and consequently set up the sum of rates out of each minimum 
! out of the individual routines that calculate rate constants and putting it here.
!
! *** NOTE THAT variations on this code still live in a few individual files: Pfold.f, getusepair.f90 and KMC.a2b.f.
! Changes here may need to be propagated to those files too.

      SUBROUTINE RATECONST_SETUP(LKSUM,DEADTS,NDEAD,COUNTNULL,CUT_UNDERFLOW)
      USE PORFUNCS
      USE COMMONS
      IMPLICIT NONE

      INTEGER :: J1, NDEAD
      DOUBLE PRECISION :: LKSUM(NMIN)
      DOUBLE PRECISION, INTENT(IN) :: CUT_UNDERFLOW
      LOGICAL :: DEADTS(NTS)
      LOGICAL, INTENT(IN) :: COUNTNULL

!!!!!!!!!!!!!!!!!!!   initial setup       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
      NDEAD=0
      DO J1=1,NMIN
         IF (NCONN(J1).LE.NCONNMIN) THEN
            NDEAD=NDEAD+1
            IF (DEBUG) PRINT*,'discarding minimum ',J1,' with ',NCONN(J1),' connections'
         ENDIF
      ENDDO
      IF (.NOT.SHANNONT) PRINT '(3(A,I8))','RATECONST_SETUP>',NDEAD,' minima with ',NCONNMIN, &
  &             ' connections or fewer will not be considered'
!
!  Flag transition states to dead-end minima as DEAD and calculate KSUM values.
!  NCONN only counted non-degenerate rearrangements as connections.
!
      DO J1=1,NMIN
         LKSUM(J1)=0.0D0
      ENDDO
      DO J1=1,NTS
         CALL CHECKTS(ETS(J1),EMIN(PLUS(J1)),EMIN(MINUS(J1)),KPLUS(J1),KMINUS(J1),NCONN(PLUS(J1)),NCONN(MINUS(J1)), &
                      PLUS(J1),MINUS(J1),COUNTNULL,CUT_UNDERFLOW,DEADTS(J1))
!
! This block is to allow checks for removal of ts in Dijkstra when downhill barrier
! is over MAXDOWNBARRIER
!
         IF (ALLOCATED(SHIFTABLE)) THEN
            IF (SHIFTABLE(J1)) DEADTS(J1)=.TRUE.
            IF (NTRYING.NE.0) DEADTS(NTRYING)=.TRUE.
         ENDIF

         IF ((.NOT.DEADTS(J1)) .AND. (PLUS(J1).NE.MINUS(J1))) THEN
!           LKSUM(PLUS(J1))=LKSUM(PLUS(J1))+EXP(KPLUS(J1)-KMEAN)
!           LKSUM(MINUS(J1))=LKSUM(MINUS(J1))+EXP(KMINUS(J1)-KMEAN)
            LKSUM(PLUS(J1))=LKSUM(PLUS(J1))+EXP(KPLUS(J1))
            LKSUM(MINUS(J1))=LKSUM(MINUS(J1))+EXP(KMINUS(J1))
         ENDIF
      ENDDO
      DO J1=1,NMIN
         IF (LKSUM(J1).GT.0.0D0) THEN
!           LKSUM(J1)=LOG(LKSUM(J1))+KMEAN
            LKSUM(J1)=LOG(LKSUM(J1))
         ENDIF
      ENDDO
!
!!!!!!!!!!!!!!!!!!!   end of initial setup   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      RETURN
      END SUBROUTINE RATECONST_SETUP
