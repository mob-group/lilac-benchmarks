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

SUBROUTINE PERC(P,NATOMS,PERCCUT,PERCT,DEBUG,MYUNIT,RIGID)
  USE COMMONS, ONLY : MODULARCURRENTN, MODULART, OPPT
  USE OPP_MOD, ONLY: OPP_MODE
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: NATOMS, MYUNIT
  LOGICAL, INTENT(IN)  :: DEBUG, RIGID
  LOGICAL, INTENT(OUT) :: PERCT
  DOUBLE PRECISION, INTENT(IN)  :: P(3*NATOMS), PERCCUT
  INTEGER NSITES, J1, J2, NCYCLE, DMIN1, DMAX1, NUNCON1
  INTEGER, ALLOCATABLE :: NDIST1(:)
  DOUBLE PRECISION DUMMY
  LOGICAL CHANGED
  LOGICAL, ALLOCATABLE :: CON(:,:)

  IF(MODULART) THEN
    NSITES = MODULARCURRENTN 
  ELSE IF (OPPT) THEN
! jwrm2> Need to remove alchemical and unit cell degrees of freedom.
    SELECT CASE (OPP_MODE)
    CASE(1)
      NSITES = NATOMS - 2
    CASE(2)
      NSITES = NATOMS - 3
    CASE(3)
      NSITES = NATOMS - 1
    END SELECT
  ELSE
    NSITES = NATOMS
  END IF
  IF (RIGID) THEN
    NSITES = NSITES/2
  ENDIF
! WRITE(MYUNIT,'(A,2L5,3I6)') 'MODULART,RIGID,NSITES,NATOMS,MODULARCURRENTN=',MODULART,RIGID,NSITES,NATOMS,MODULARCURRENTN
  IF(ALLOCATED(NDIST1)) DEALLOCATE(NDIST1)
  ALLOCATE(NDIST1(NSITES), CON(NSITES,NSITES))

  CON(1:NSITES,1:NSITES)=.FALSE.
  DO J1=1,NSITES
    DO J2=J1+1,NSITES
      DUMMY=(P(3*(J2-1)+1)-P(3*(J1-1)+1))**2+(P(3*(J2-1)+2)-P(3*(J1-1)+2))**2+(P(3*(J2-1)+3)-P(3*(J1-1)+3))**2
      IF (DUMMY.LT.PERCCUT) THEN
        CON(J2,J1)=.TRUE.
        CON(J1,J2)=.TRUE.
!       IF (DEBUG) WRITE(MYUNIT,'(A,2I8)') 'perc> connecting atoms ',J1,J2
      ENDIF
    ENDDO
  ENDDO

! 
! Check that we have a percolating constraint network.
!
  NDIST1(1:NSITES)=1000000
  NDIST1(1)=0
  NCYCLE=0
5 CHANGED=.FALSE.
  NCYCLE=NCYCLE+1
  DMIN1=100000
  DMAX1=0
  NUNCON1=0
  DO J1=1, NSITES
    IF (NDIST1(J1).EQ.0) CYCLE ! minimum 1
    DO J2=1, NSITES
      IF (CON(J2,J1)) THEN
        IF (NDIST1(J2)+1.LT.NDIST1(J1)) THEN
          CHANGED=.TRUE.
          NDIST1(J1)=NDIST1(J2)+1
        ENDIF
        IF (NDIST1(J1)+1.LT.NDIST1(J2)) THEN
          CHANGED=.TRUE.
          NDIST1(J2)=NDIST1(J1)+1
        ENDIF
      ENDIF
    ENDDO
    IF ((NDIST1(J1).GT.DMAX1).AND.(NDIST1(J1).NE.1000000)) DMAX1=NDIST1(J1)
    IF (NDIST1(J1).LT.DMIN1) DMIN1=NDIST1(J1)
    IF (NDIST1(J1).EQ.1000000) NUNCON1=NUNCON1+1
  ENDDO
!  PRINT *,'DMIN1,DMAX1,NUNCON1,NCYCLE,CHANGED=',DMIN1,DMAX1,NUNCON1,NCYCLE,CHANGED
  IF (CHANGED) GOTO 5
  IF (DEBUG) WRITE(MYUNIT,'(3(A,I8))') 'perc> steps to atom 1 converged in ',NCYCLE-1, &
     &                    ' cycles; maximum=',DMAX1,' disconnected=',NUNCON1
  PERCT=.TRUE.
  IF (NUNCON1.GT.0) THEN
     PERCT=.FALSE.
!    IF (DEBUG) WRITE(MYUNIT,'(3G20.10)') P(1:3*NATOMS)
     IF (DEBUG) THEN
        DO J1=1,NSITES
           WRITE(MYUNIT,'(2I6,3G20.10)') J1,NDIST1(J1),P(3*(J1-1)+1),P(3*(J1-1)+2),P(3*(J1-1)+3)
        ENDDO
     ENDIF
  ENDIF

  DEALLOCATE(NDIST1,CON)

END SUBROUTINE PERC

SUBROUTINE PERCGROUP(P,NATOMS,PERCCUT,DEBUG,MYUNIT,RIGID)
! sf344> group the structures into connected parts, in order to move such parts together during MC step taking
! (should be useful for hierarchical global optimization where parts of the cluster should be moved 
! around as a rigid body)

  USE COMMONS, ONLY : atomingroup, modulart, modularcurrentn
  IMPLICIT NONE

  INTEGER, INTENT(IN)  :: NATOMS, MYUNIT
  LOGICAL, INTENT(IN)  :: DEBUG, RIGID
  DOUBLE PRECISION, INTENT(IN)  :: PERCCUT
  INTEGER NSITES, J1, J2, NCYCLE, DMIN1, DMAX1, NUNCON1, groupindex, nsitesskip,nunconnsave,k1,k2
  INTEGER, ALLOCATABLE :: NDIST1(:)
  DOUBLE PRECISION DUMMY, COORDSNEW(3*NATOMS), P(3*NATOMS) !P is now intent inout
  LOGICAL CHANGED
  LOGICAL, ALLOCATABLE :: CON(:,:)
 

  IF(MODULART) THEN
    NSITES = MODULARCURRENTN 
  ELSE
    NSITES = NATOMS
  END IF
  IF (RIGID) THEN
    NSITES = NSITES/2
  ENDIF

COORDSNEW(:)=P(:)
P(:)=0.0D0  !reset P, we will populate it later in an ordered fashion
NUNCONNSAVE=NSITES
groupindex = 0
! WRITE(MYUNIT,'(A,2L5,3I6)') 'MODULART,RIGID,NSITES,NATOMS,MODULARCURRENTN=',MODULART,RIGID,NSITES,NATOMS,MODULARCURRENTN
  IF(ALLOCATED(NDIST1)) DEALLOCATE(NDIST1)
  IF(ALLOCATED(CON)) DEALLOCATE(CON)
  ALLOCATE(NDIST1(NSITES), CON(NSITES,NSITES))
  IF(ALLOCATED(atomingroup)) DEALLOCATE(atomingroup)
  ALLOCATE(atomingroup(nsites)) ! this will list the group index for each particle
atomingroup(:)=0
DO
! sf344>  after the first loop we have NUNNCONSAVE particles disconnected from particle 1, we need 
! to check percolation only for these. Array P is being reordered so that each 
! detected disconnected group is populating it sequentially. We take always AMERI...khm..CONNECTED particles first.
 groupindex = groupindex + 1
 NSITESSKIP=NSITES-NUNCONNSAVE !this is 0 for the first loop
  IF(DEBUG) WRITE(MYUNIT,*) 'percgroup> Building group ', groupindex, nunconnsave, nsitesskip
! NSITES=NUNCONNSAVE

  CON(1:NSITES,1:NSITES)=.FALSE.
  DO J1=NSITESSKIP+1,NSITES
    DO J2=J1+1,NSITES
      DUMMY=(COORDSNEW(3*(J2-1)+1)-COORDSNEW(3*(J1-1)+1))**2+&
            (COORDSNEW(3*(J2-1)+2)-COORDSNEW(3*(J1-1)+2))**2+&
            (COORDSNEW(3*(J2-1)+3)-COORDSNEW(3*(J1-1)+3))**2
      IF (DUMMY.LT.PERCCUT) THEN
        CON(J2,J1)=.TRUE.
        CON(J1,J2)=.TRUE.
      ENDIF
    ENDDO
  ENDDO

! 
! Check that we have a percolating constraint network.
!
  NDIST1(NSITESSKIP+1:NSITES)=1000000
  NDIST1(NSITESSKIP+1)=0
  NCYCLE=0
55 CHANGED=.FALSE.
  NCYCLE=NCYCLE+1
  DMIN1=100000
  DMAX1=0
  NUNCON1=0
  DO J1=NSITESSKIP+1, NSITES
    IF (NDIST1(J1).EQ.0) CYCLE ! first particle in the group
    DO J2=NSITESSKIP+1, NSITES
      IF (CON(J2,J1)) THEN
        IF (NDIST1(J2)+1.LT.NDIST1(J1)) THEN
          CHANGED=.TRUE.
          NDIST1(J1)=NDIST1(J2)+1
        ENDIF
        IF (NDIST1(J1)+1.LT.NDIST1(J2)) THEN
          CHANGED=.TRUE.
          NDIST1(J2)=NDIST1(J1)+1
        ENDIF
      ENDIF
    ENDDO
    IF ((NDIST1(J1).GT.DMAX1).AND.(NDIST1(J1).NE.1000000)) DMAX1=NDIST1(J1)
    IF (NDIST1(J1).LT.DMIN1) DMIN1=NDIST1(J1)
    IF (NDIST1(J1).EQ.1000000) NUNCON1=NUNCON1+1
  ENDDO
!  PRINT *,'DMIN1,DMAX1,NUNCON1,NCYCLE,CHANGED=',DMIN1,DMAX1,NUNCON1,NCYCLE,CHANGED
  IF (CHANGED) GOTO 55
  IF (DEBUG) WRITE(MYUNIT,'(3(A,I8))') ' percgroup> steps to first atom in the group converged in ',NCYCLE-1, &
     &                    ' cycles; maximum=',DMAX1,' disconnected=',NUNCON1
  NUNCONNSAVE=NUNCON1
  IF (NUNCON1.GT.0) THEN
! reorder the coordinates
     K1=NSITESSKIP
     K2=0
!        WRITE(MYUNIT,*) 'before reordering', k1, k2, nsitesskip,nuncon1
!        WRITE(MYUNIT,*) NDIST1(:)
!        WRITE(MYUNIT,*) 'original coordinates: '
!        WRITE(MYUNIT,*) COORDSNEW(:)
     DO J1=NSITESSKIP+1,NSITES
!        WRITE(MYUNIT,*) 'sf344 debug> J1=', J1, NDIST1(J1)
        IF(NDIST1(J1)<1000000) THEN
          K1=K1+1
          atomingroup(k1)=groupindex
          P(3*K1-2:3*K1)=COORDSNEW(3*J1-2:3*J1) !copy over coordinates for connected group
          IF(RIGID) THEN
            !copy over angle-axis coordinates as well
            P(3*K1-2+3*NATOMS/2:3*K1+3*NATOMS/2)=COORDSNEW(3*J1-2+3*NATOMS/2:3*J1+3*NATOMS/2)
          END IF
        ELSE !disconnected particle, coordinates should go to the end of the array
          K2=K2+1
!          WRITE(MYUNIT,*) 'sf344 debug> ', 3*NSITES-3*NUNCON1+3*K2-2
!          WRITE(MYUNIT,*) COORDSNEW(3*J1-2:3*J1)
          P(3*NSITES-3*NUNCON1+3*K2-2:3*NSITES-3*NUNCON1+3*K2)=COORDSNEW(3*J1-2:3*J1)
          IF(RIGID) THEN
            P(3*NSITES-3*NUNCON1+3*K2+3*NATOMS/2-2:3*NSITES-3*NUNCON1+3*K2+3*NATOMS/2)=COORDSNEW(3*J1-2+3*NATOMS/2:3*J1+3*NATOMS/2)
          END IF
        END IF
     END DO
!        WRITE(MYUNIT,*) 'after reordering', k1, k2, nsitesskip
!        WRITE(MYUNIT,*) 'reordered coordinates: '
!        WRITE(MYUNIT,*) P(:)
! save the reordered coordinate array
     COORDSNEW(:)=P(:)
!    IF (DEBUG) WRITE(MYUNIT,'(3G20.10)') P(1:3*NATOMS)
!     IF (DEBUG) THEN !this debug printing does not make sense anymore, we already reordered the coordinates
!        DO J1=1,NSITES
!           WRITE(MYUNIT,'(2I6,3G20.10)') J1,NDIST1(J1),P(3*(J1-1)+1),P(3*(J1-1)+2),P(3*(J1-1)+3)
!        ENDDO
!     ENDIF
  ELSE !NUNCON1==0, so every atom is connected in the last group, or the very last atom is disconnected
   atomingroup(NSITESSKIP+1:NSITES)=groupindex
     IF(DEBUG) THEN
        WRITE(MYUNIT,*) 'percgroup> finished grouping, particles in the following groups: '
        WRITE(MYUNIT,*) atomingroup(:)
     END IF
     EXIT !exit the outer loop
  ENDIF
END DO

!        WRITE(MYUNIT,*) 'percgroup> final reordered coordinates: '
!        WRITE(MYUNIT,*) P(:)

  DEALLOCATE(NDIST1,CON)
END SUBROUTINE PERCGROUP
