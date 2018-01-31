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
! linear interpolation for unres with iterative refinement of
! longer or shorter route around each dihedral angle
SUBROUTINE UNGUESSPATH(START,FINISH,NCOORDS,EDIFFTOLLOCAL,NATOMS)
USE KEY
USE MODCHARMM
USE MODGUESS
USE MODUNRES
use porfuncs
IMPLICIT NONE
INTEGER, INTENT(IN) :: NCOORDS, NATOMS
DOUBLE PRECISION, INTENT(IN) :: START(NCOORDS), FINISH(NCOORDS), EDIFFTOLLOCAL
DOUBLE PRECISION, ALLOCATABLE :: INTERPENERGY(:), INTERPENERGYSAVE(:), PEPCOORDS(:), SMALLSTEP(:), LARGESTEP(:), DELTA(:)
DOUBLE PRECISION :: MAXE, SUME, NEWMAXE
INTEGER :: JMAXE(1), NEWJMAXE(1)
INTEGER :: ISTAT, NLONG
LOGICAL, ALLOCATABLE :: SHORTER(:)
DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793D0
CHARACTER(LEN=80) LFNAME

INTEGER :: J1

IF (.NOT.UNRST) THEN
   WRITE(*,'(A)') 'ERROR - unguesspath is only for UNRES potential'
   STOP
ENDIF
IF (UNRST) ALLOCATE(PEPCOORDS(3*NATOMS),SMALLSTEP(NCOORDS),LARGESTEP(NCOORDS),DELTA(NCOORDS))
ALLOCATE(SHORTER(NCOORDS))
! define the angular step for the shortest and longest angular routes between start and
! finish for UNRES.
DO J1=1,NCOORDS
   IF (ABS(FINISH(J1)-START(J1)).LE.PI) THEN
      SMALLSTEP(J1)=FINISH(J1)-START(J1)
      IF (FINISH(J1).GT.START(J1)) THEN
         LARGESTEP(J1)=FINISH(J1)-START(J1)-2*PI
      ELSE
         LARGESTEP(J1)=FINISH(J1)-START(J1)+2*PI
      ENDIF
   ELSE
      LARGESTEP(J1)=FINISH(J1)-START(J1)
      IF (FINISH(J1).GT.START(J1)) THEN
         SMALLSTEP(J1)=FINISH(J1)-START(J1)-2*PI
      ELSE
         SMALLSTEP(J1)=FINISH(J1)-START(J1)+2*PI
      ENDIF
   ENDIF
   DELTA(J1)=SMALLSTEP(J1)
   SHORTER(J1)=.TRUE.
ENDDO
NINTERP=MAXGCYCLES
WRITE(*,'(2(A,I6))') ' Number of interior points for interpolated path in UNGUESSPATH=',NINTERP
ALLOCATE(INTERPENERGY(NINTERP+2),INTERPENERGYSAVE(NINTERP+2))
! to optimise the order of the steps we need to recalculate INTERPENERGY for all configurations
! The energy at steps 1 and NINTERP+2 are just the energies of the START and FINISH minima

CALL GENERGIES(1,NINTERP+2,.FALSE.) ! fill the initial INTERPENERGY vector with interpolated energies
INTERPENERGYSAVE(1:NINTERP+2)=INTERPENERGY(1:NINTERP+2)

MAXE=MAXVAL(INTERPENERGY(2:NINTERP+1))
JMAXE=MAXLOC(INTERPENERGY(2:NINTERP+1)) ! note that MAXLOC returns a vector of dimension one, not a scalar
SUME=SUM(INTERPENERGY(2:NINTERP+1))

WRITE(*,'(A,2G20.10)') 'maximum energy and mean of intervening points=',MAXE,SUME/(NINTERP)
WRITE(*,'(A,I6)') 'maximum is at location ',JMAXE(1)
PRINT*,'energies:'
WRITE(*,'(I6,G20.10)') (J1,INTERPENERGY(J1),J1=1,NINTERP+2)
! try changing the direction (SHORTER=true or false) for each degree of freedom in turn
! until we get through them all with no improvement

main: DO 
   DO J1=1,NCOORDS ! try changing the interpolation direction for each coordinate in turn
   
! first try using LARGESTEP instead of SMALLSTEP for the degree of freedom
! corresponding to JMAXE(1). Or try SMALLSTEP instead of LARGESTEP if we've changed this
! an odd number of times already!

      IF (SHORTER(J1)) THEN
!        WRITE(*,'(A,I5)') ' trying longer interpolation for coordinate ',J1
         DELTA(J1)=LARGESTEP(J1)
         SHORTER(J1)=.FALSE.
      ELSE
!        WRITE(*,'(A,I5)') ' trying shorter interpolation for coordinate ',J1
         DELTA(J1)=SMALLSTEP(J1)
         SHORTER(J1)=.TRUE.
      ENDIF
      CALL GENERGIES(2,NINTERP+1,.FALSE.) 
      NEWMAXE=MAXVAL(INTERPENERGY(2:NINTERP+1))
      NEWJMAXE=MAXLOC(INTERPENERGY(2:NINTERP+1))
      IF (MAXE-NEWMAXE.GT.EDIFFTOLLOCAL) THEN ! accept change
         WRITE(*,'(A,G20.10,A,I6)') ' maximum interior energy is now ',NEWMAXE,' at position ',NEWJMAXE(1)
         MAXE=NEWMAXE
         INTERPENERGYSAVE(2:NINTERP+1)=INTERPENERGY(2:NINTERP+1)
         CYCLE main ! go back and start from the first coordinate again
      ELSE
!        WRITE(*,'(A,G20.10,A,I6))') ' maximum interior energy is now ',NEWMAXE,' at position ',NEWJMAXE(1)
         INTERPENERGY(2:NINTERP+1)=INTERPENERGYSAVE(2:NINTERP+1)
         IF (SHORTER(J1)) THEN
            DELTA(J1)=LARGESTEP(J1)
            SHORTER(J1)=.FALSE.
         ELSE
            DELTA(J1)=SMALLSTEP(J1)
            SHORTER(J1)=.TRUE.
         ENDIF
      ENDIF
   ENDDO
   EXIT main
ENDDO main

SUME=SUM(INTERPENERGY(2:NINTERP+1))
WRITE(*,'(A,2G20.10)') 'maximum energy and mean=',MAXE,SUME/(NINTERP)
WRITE(*,'(A,I6)') 'maximum is at location ',JMAXE(1)
PRINT*,' energies:'
WRITE(*,'(I6,G20.10)') (J1,INTERPENERGY(J1),J1=1,NINTERP+2)
NLONG=0
DO J1=1,NCOORDS
   IF (.NOT.SHORTER(J1)) NLONG=NLONG+1
ENDDO
PRINT '(A,I5)',' number of longer interpolations=',NLONG
! dump xyz path
CALL GENERGIES(1,NINTERP+2,.TRUE.) 

DEALLOCATE(INTERPENERGY,INTERPENERGYSAVE,SHORTER)
IF (ALLOCATED(PEPCOORDS)) DEALLOCATE(PEPCOORDS,DELTA,SMALLSTEP,LARGESTEP)
CALL FLUSH(6)

CONTAINS
   SUBROUTINE GENERGIES(BOTTOM,TOP,DUMPPATH)
   INTEGER, INTENT(IN) :: BOTTOM, TOP
   INTEGER :: K1, K2, LBOTTOM, LTOP, K3
   LOGICAL, INTENT(IN) :: DUMPPATH
   DOUBLE PRECISION :: X(3*NATOMS), RMS, GRAD(3*NATOMS)
               ! there should be a way of not 
               ! declaring the Hessian unless it is actually needed? Use (*,*) in potential
               ! and routines called from it?
   DOUBLE PRECISION :: ENERGY

   LBOTTOM=BOTTOM
   LTOP=TOP
   IF (LBOTTOM.GT.LTOP) THEN
      WRITE(*,'(A)') ' ERROR - LBOTTOM,LTOP=',LBOTTOM,LTOP
      STOP
   ENDIF
   IF (DUMPPATH) THEN
      IF (FILTH.EQ.0) THEN
         OPEN(UNIT=7,FILE='guess.xyz',STATUS='UNKNOWN')
      ELSE
         LFNAME='guess.xyz.'//TRIM(ADJUSTL(FILTHSTR))
         OPEN(UNIT=7,FILE=TRIM(ADJUSTL(LFNAME)),STATUS='UNKNOWN')
      ENDIF
      IF (FILTH.EQ.0) THEN
         OPEN(UNIT=8,FILE='guess.unres.xyz',STATUS='UNKNOWN')
      ELSE
         LFNAME='guess.unres.xyz.'//TRIM(ADJUSTL(FILTHSTR))
         OPEN(UNIT=8,FILE=TRIM(ADJUSTL(LFNAME)),STATUS='UNKNOWN')
      ENDIF
      IF (LBOTTOM.NE.1) THEN
         PRINT*,'here'
         CALL var_to_geom(NCOORDS,START)
         CALL chainbuild
         DO K2=1,NRES
            WRITE(7,'(3G20.10)') C(1,K2),C(2,K2),C(3,K2) ! backbone
            WRITE(7,'(3G20.10)') C(1,K2+NRES),C(2,K2+NRES),C(3,K2+NRES) ! side chains
         ENDDO
         DO K2=1,(NATOMS/2)-1 ! jmc add peptide atoms...
            DO K3=1,3
               PEPCOORDS(6*(K2-1)+K3)=(2.0D0*C(K3,K2)+C(K3,K2+1))/3.0D0
               PEPCOORDS(6*(K2-1)+K3+3)=(C(K3,K2)+2.0D0*C(K3,K2+1))/3.0D0
            ENDDO
         ENDDO
         WRITE(8,'(I6)') 2*NATOMS-2
         WRITE(8,'(A)') ' '
         WRITE(8,'(A2,3G20.10)') ('C  ',C(1,K2),C(2,K2),C(3,K2),K2=1,NATOMS)
         WRITE(8,'(A2,4X,3F20.10)') ('O ',PEPCOORDS(6*(K2-1)+1),PEPCOORDS(6*(K2-1)+2),PEPCOORDS(6*(K2-1)+3) &
                                   ,K2=1,(NATOMS/2)-1)
         WRITE(8,'(A2,4X,3F20.10)') ('N ',PEPCOORDS(6*(K2-1)+4),PEPCOORDS(6*(K2-1)+5),PEPCOORDS(6*(K2-1)+6) &
                                   ,K2=1,(NATOMS/2)-1)
      ENDIF
   ENDIF
   X(1:NCOORDS)=START(1:NCOORDS)
! Remove exact atom superpositions using some random numerical noise with DPRAND
   DO K1=LBOTTOM,LTOP
      DO K2=1,NCOORDS ! linear interpolation
         X(K2)=START(K2)+DELTA(K2)*(K1-1)/(NINTERP+1)
      ENDDO
      CALL var_to_geom(NCOORDS,X)
      CALL chainbuild
      IF (DUMPPATH) THEN
         DO K2=1,NRES
            WRITE(7,'(3G20.10)') C(1,K2),C(2,K2),C(3,K2) ! backbone
            WRITE(7,'(3G20.10)') C(1,K2+NRES),C(2,K2+NRES),C(3,K2+NRES) ! side chains
         ENDDO
         DO K2=1,(NATOMS/2)-1 ! jmc add peptide atoms...
            DO K3=1,3
               PEPCOORDS(6*(K2-1)+K3)=(2.0D0*C(K3,K2)+C(K3,K2+1))/3.0D0
               PEPCOORDS(6*(K2-1)+K3+3)=(C(K3,K2)+2.0D0*C(K3,K2+1))/3.0D0
            ENDDO
         ENDDO
         WRITE(8,'(I6)') 2*NATOMS-2
         WRITE(8,'(A)') ' '
         WRITE(8,'(A2,3G20.10)') ('C  ',C(1,K2),C(2,K2),C(3,K2),K2=1,NATOMS)
         WRITE(8,'(A2,4X,3F20.10)') ('O ',PEPCOORDS(6*(K2-1)+1),PEPCOORDS(6*(K2-1)+2),PEPCOORDS(6*(K2-1)+3) &
                                   ,K2=1,(NATOMS/2)-1)
         WRITE(8,'(A2,4X,3F20.10)') ('N ',PEPCOORDS(6*(K2-1)+4),PEPCOORDS(6*(K2-1)+5),PEPCOORDS(6*(K2-1)+6) &
                                   ,K2=1,(NATOMS/2)-1)
      ELSE
         CALL POTENTIAL(X,ENERGY,GRAD,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
         INTERPENERGY(K1)=ENERGY
      ENDIF
   ENDDO
   IF (DUMPPATH) CLOSE(7)
   IF (DUMPPATH.AND.UNRST) CLOSE(8)

   END SUBROUTINE GENERGIES

END SUBROUTINE UNGUESSPATH
