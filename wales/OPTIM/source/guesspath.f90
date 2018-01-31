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
SUBROUTINE GUESSPATH(START,FINISH,NCOORDS,EDIFFTOLLOCAL,NATOMS)
USE KEY
USE MODCHARMM
USE MODGUESS
USE MODUNRES
use porfuncs
IMPLICIT NONE
INTEGER, INTENT(IN) :: NCOORDS, NATOMS
DOUBLE PRECISION, INTENT(IN) :: START(NCOORDS), FINISH(NCOORDS), EDIFFTOLLOCAL
INTEGER, ALLOCATABLE :: STEPORDER(:)
DOUBLE PRECISION, ALLOCATABLE :: INTERPENERGY(:), INTERPENERGYSAVE(:), PEPCOORDS(:), SMALLSTEP(:), LARGESTEP(:), DELTA(:)
DOUBLE PRECISION :: MAXE, SUME, NEWMAXE
INTEGER :: JMAXE(1), NEWJMAXE(1)
INTEGER :: NDUMMY, ISTAT
LOGICAL :: IMPROVED
LOGICAL, ALLOCATABLE :: SHORTER(:)
DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793D0
CHARACTER(LEN=80) LFNAME

LOGICAL :: SEQUENTIAL(NCOORDS), DOSEQ(NCOORDS)
INTEGER :: NSEQ, J1, NCOUNT

IF (CHRMMT) THEN
   WRITE(*,'(A)') 'For CHARMM there may be more internal coordinates than Cartesians and the dimensions need attention'
   STOP
ENDIF
IF (UNRST) ALLOCATE(PEPCOORDS(3*NATOMS),SMALLSTEP(NCOORDS),LARGESTEP(NCOORDS),DELTA(NCOORDS))
ALLOCATE(SHORTER(NCOORDS))
! first determine the number of coordinates to be treated sequentially according to GTHRESHOLD
NSEQ=0
! define the angular step for the shortest and longest angular routes between start and
! finish for UNRES.
IF (UNRST) THEN
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
      SEQUENTIAL(J1)=.FALSE.
      DOSEQ(J1)=.FALSE.
      IF (ABS(SMALLSTEP(J1)).GT.GTHRESHOLD) THEN
         SEQUENTIAL(J1)=.TRUE.
         DOSEQ(J1)=.TRUE.
         NSEQ=NSEQ+1
      ENDIF
!     WRITE(*,'(A,2I5,2G20.10,4X,L1)') 'J1,NSEQ,SMALLSTEP(J1),LARGESTEP(J1)=',J1,NSEQ,SMALLSTEP(J1),LARGESTEP(J1),SEQUENTIAL(J1)
      DELTA(J1)=SMALLSTEP(J1)
      SHORTER(J1)=.TRUE.
   ENDDO
ELSE
   DO J1=1,NCOORDS
      SEQUENTIAL(J1)=.FALSE.
      DOSEQ(J1)=.FALSE.
      IF (ABS(START(J1)-FINISH(J1)).GT.GTHRESHOLD) THEN
         SEQUENTIAL(J1)=.TRUE.
         DOSEQ(J1)=.TRUE.
         NSEQ=NSEQ+1
         WRITE(*,'(A,2I5,2G20.10)') 'J1,NSEQ,START,FINISH=',J1,NSEQ,START(J1),FINISH(J1)
      ENDIF
   ENDDO
ENDIF
IF (NSEQ.EQ.0) THEN
   WRITE(*,'(A,F12.3)') ' No sequential coordinates for a threshold of ',GTHRESHOLD
   NINTERP=0
   CALL FLUSH(6)
   RETURN
ENDIF
WRITE(*,'(2(A,I6))') ' Number of coordinates to be treated sequentially in GUESSPATH=',NSEQ, &
                     ' out of ', NCOORDS
NINTERP=NSEQ*GSTEPS-1
ALLOCATE(STEPORDER(NSEQ*GSTEPS),INTERPENERGY(NSEQ*GSTEPS),INTERPENERGYSAVE(NSEQ*GSTEPS))
NCOUNT=0
! initial order of steps for sequential coordinates
DO J1=1,NCOORDS
   IF (SEQUENTIAL(J1)) THEN
      NCOUNT=NCOUNT+1
      STEPORDER(GSTEPS*(NCOUNT-1)+1:GSTEPS*NCOUNT)=J1
   ENDIF
ENDDO
! to optimise the order of the steps we need to recalculate INTERPENERGY for all configurations
! between the elements of STEPORDER that have been changed, and save the current optimal values
! of the highest configuration and sum of energies. 
! 
! The energy at step NSEQ*GSTEPS is just the energy of the minimum in FINISH - should never change

CALL GENERGIES(1,NSEQ*GSTEPS,.FALSE.) ! fill the initial INTERPENERGY vector with interpolated energies
INTERPENERGYSAVE(1:NSEQ*GSTEPS)=INTERPENERGY(1:NSEQ*GSTEPS)

MAXE=MAXVAL(INTERPENERGY(1:NSEQ*GSTEPS-1))
JMAXE=MAXLOC(INTERPENERGY(1:NSEQ*GSTEPS-1)) ! note that MAXLOC returns a vector of dimension one, not a scalar
SUME=SUM(INTERPENERGY(1:NSEQ*GSTEPS))

WRITE(*,'(A,2G20.10)') 'maximum energy and  mean=',MAXE,SUME/(NSEQ*GSTEPS-1)
WRITE(*,'(A,I6)') 'maximum is at location ',JMAXE(1)
PRINT*,'order of coordinate steps, energies:'
WRITE(*,'(2I6,G20.10)') (J1,STEPORDER(J1),INTERPENERGY(J1),J1=1,NSEQ*GSTEPS)
! try changing the step in position JMAXE(1) with all the others in turn to find the lowest
! MAXE that results
DO 
   IMPROVED=.FALSE.
! For UNRES first try using LARGESTEP instead of SMALLSTEP for the degree of freedom
! corresponding to JMAXE(1). Or try SMALLSTEP instead of LARGESTEP if we've changed this
! an odd number of times already!
   IF (UNRST) THEN
      IF (SHORTER(STEPORDER(JMAXE(1)))) THEN
         WRITE(*,'(A,I5)') ' trying longer interpolation for coordinate ',STEPORDER(JMAXE(1))
         DELTA(STEPORDER(JMAXE(1)))=LARGESTEP(STEPORDER(JMAXE(1)))
         SHORTER(STEPORDER(JMAXE(1)))=.FALSE.
      ELSE
         WRITE(*,'(A,I5)') ' trying shorter interpolation for coordinate ',STEPORDER(JMAXE(1))
         DELTA(STEPORDER(JMAXE(1)))=SMALLSTEP(STEPORDER(JMAXE(1)))
         SHORTER(STEPORDER(JMAXE(1)))=.TRUE.
      ENDIF
      CALL GENERGIES(1,NSEQ*GSTEPS-1,.FALSE.) 
      NEWMAXE=MAXVAL(INTERPENERGY(1:NSEQ*GSTEPS-1))
      NEWJMAXE=MAXLOC(INTERPENERGY(1:NSEQ*GSTEPS-1))
      WRITE(*,'(A,G20.10,A,I6)') ' maximum energy is now ',NEWMAXE,' at position ',NEWJMAXE(1)
      IF (MAXE-NEWMAXE.GT.EDIFFTOLLOCAL) THEN ! accept change
         MAXE=NEWMAXE
         INTERPENERGYSAVE(1:NSEQ*GSTEPS-1)=INTERPENERGY(1:NSEQ*GSTEPS-1)
         IMPROVED=.TRUE.
         IF (JMAXE(1).NE.NEWJMAXE(1)) THEN
            JMAXE(1)=NEWJMAXE(1)
            CYCLE
         ENDIF
      ELSE
         INTERPENERGY(1:NSEQ*GSTEPS-1)=INTERPENERGYSAVE(1:NSEQ*GSTEPS-1)
         IF (SHORTER(STEPORDER(JMAXE(1)))) THEN
            DELTA(STEPORDER(JMAXE(1)))=LARGESTEP(STEPORDER(JMAXE(1)))
            SHORTER(STEPORDER(JMAXE(1)))=.FALSE.
         ELSE
            DELTA(STEPORDER(JMAXE(1)))=SMALLSTEP(STEPORDER(JMAXE(1)))
            SHORTER(STEPORDER(JMAXE(1)))=.TRUE.
         ENDIF
      ENDIF

      DO J1=1,NCOORDS
         IF (SEQUENTIAL(J1).AND.(J1.NE.STEPORDER(JMAXE(1)))) THEN
            IF (SHORTER(J1)) THEN
               WRITE(*,'(A,I5)') ' trying longer interpolation for coordinate ',J1
               DELTA(J1)=LARGESTEP(J1)
               SHORTER(J1)=.FALSE.
            ELSE
               WRITE(*,'(A,I5)') ' trying shorter interpolation for coordinate ',J1
               DELTA(J1)=SMALLSTEP(J1)
               SHORTER(J1)=.TRUE.
            ENDIF
            CALL GENERGIES(1,NSEQ*GSTEPS-1,.FALSE.) 
            NEWMAXE=MAXVAL(INTERPENERGY(1:NSEQ*GSTEPS-1))
            NEWJMAXE=MAXLOC(INTERPENERGY(1:NSEQ*GSTEPS-1))
            WRITE(*,'(A,G20.10,A,I6)') ' maximum energy is now ',NEWMAXE,' at position ',NEWJMAXE(1)
            IF (MAXE-NEWMAXE.GT.EDIFFTOLLOCAL) THEN ! accept change
               MAXE=NEWMAXE
               INTERPENERGYSAVE(1:NSEQ*GSTEPS-1)=INTERPENERGY(1:NSEQ*GSTEPS-1)
               IMPROVED=.TRUE.
               IF (JMAXE(1).NE.NEWJMAXE(1)) THEN
                  JMAXE(1)=NEWJMAXE(1)
                  CYCLE
               ENDIF
            ELSE
               INTERPENERGY(1:NSEQ*GSTEPS-1)=INTERPENERGYSAVE(1:NSEQ*GSTEPS-1)
               IF (SHORTER(J1)) THEN
                  DELTA(J1)=LARGESTEP(J1)
                  SHORTER(J1)=.FALSE.
               ELSE
                  DELTA(J1)=SMALLSTEP(J1)
                  SHORTER(J1)=.TRUE.
               ENDIF
            ENDIF
         ENDIF
      ENDDO
   ENDIF

   WRITE(*,'(A,I6,A,G20.10,A,L3,A,G20.10)') ' try toggling linear/sequential interpolation for eligible coordinates'
   DO J1=1,NCOORDS
      IF (SEQUENTIAL(J1)) THEN
         IF (DOSEQ(J1)) THEN
            WRITE(*,'(A,I5)') ' trying linear interpolation for coordinate ',J1
            DOSEQ(J1)=.FALSE.
         ELSE
            WRITE(*,'(A,I5)') ' trying sequential interpolation for coordinate ',J1
            DOSEQ(J1)=.TRUE.
         ENDIF
         CALL GENERGIES(1,NSEQ*GSTEPS-1,.FALSE.)
         NEWMAXE=MAXVAL(INTERPENERGY(1:NSEQ*GSTEPS-1))
         NEWJMAXE=MAXLOC(INTERPENERGY(1:NSEQ*GSTEPS-1))
         WRITE(*,'(A,G20.10,A,I6)') ' maximum energy is now ',NEWMAXE,' at position ',NEWJMAXE(1)
         IF (MAXE-NEWMAXE.GT.EDIFFTOLLOCAL) THEN ! accept change
            PRINT*,'accepting linear/sequential change'
            MAXE=NEWMAXE
            INTERPENERGYSAVE(1:NSEQ*GSTEPS-1)=INTERPENERGY(1:NSEQ*GSTEPS-1)
            IMPROVED=.TRUE.
            IF (JMAXE(1).NE.NEWJMAXE(1)) THEN
               JMAXE(1)=NEWJMAXE(1)
               CYCLE
            ENDIF
         ELSE
            PRINT*,'rejecting linear/sequential change'
            INTERPENERGY(1:NSEQ*GSTEPS-1)=INTERPENERGYSAVE(1:NSEQ*GSTEPS-1)
            IF (DOSEQ(J1)) THEN
               DOSEQ(J1)=.FALSE.
            ELSE
               DOSEQ(J1)=.TRUE.
            ENDIF
         ENDIF
      ENDIF
   ENDDO

   IF (UNRST) THEN
      WRITE(*,'(A,I6,A,G20.10,A,L3,A,G20.10)') ' starting series of exchanges for step in position ',JMAXE(1), &
                              ' highest energy=',MAXE,' SHORTER=',SHORTER(STEPORDER(JMAXE(1))),' step=',DELTA(STEPORDER(JMAXE(1)))
   ELSE
      WRITE(*,'(A,I6,A,G20.10,A,L3,A,G20.10)') ' starting series of exchanges for step in position ',JMAXE(1), &
                              ' highest energy=',MAXE
   ENDIF
   DO J1=1,NSEQ*GSTEPS
      IF (J1.EQ.JMAXE(1)) CYCLE
      NDUMMY=STEPORDER(J1)
      STEPORDER(J1)=STEPORDER(JMAXE(1))
      STEPORDER(JMAXE(1))=NDUMMY
      CALL GENERGIES(J1,JMAXE(1),.FALSE.) ! fill the initial INTERPENERGY vector with interpolated energies
      NEWMAXE=MAXVAL(INTERPENERGY(1:NSEQ*GSTEPS-1))
      NEWJMAXE=MAXLOC(INTERPENERGY(1:NSEQ*GSTEPS-1))
      WRITE(*,'(2(A,I6),A,G20.10,2(A,I6))') ' after swapping steps at locations ',J1,' and ',JMAXE(1), &
                            ' maximum energy=',NEWMAXE,' at position ',NEWJMAXE(1),' coordinate ',STEPORDER(NEWJMAXE(1))
      IF (MAXE-NEWMAXE.GT.EDIFFTOLLOCAL) THEN ! accept swap
         MAXE=NEWMAXE
         INTERPENERGYSAVE(1:NSEQ*GSTEPS-1)=INTERPENERGY(1:NSEQ*GSTEPS-1)
         IMPROVED=.TRUE.
         IF (JMAXE(1).NE.NEWJMAXE(1)) THEN
            JMAXE(1)=NEWJMAXE(1)
            EXIT
         ENDIF
      ELSE ! undo swap
         INTERPENERGY(1:NSEQ*GSTEPS-1)=INTERPENERGYSAVE(1:NSEQ*GSTEPS-1)
         NDUMMY=STEPORDER(J1)
         STEPORDER(J1)=STEPORDER(JMAXE(1))
         STEPORDER(JMAXE(1))=NDUMMY
      ENDIF
   ENDDO
   IF (.NOT.IMPROVED) EXIT
ENDDO 

SUME=SUM(INTERPENERGY(1:NSEQ*GSTEPS-1))
WRITE(*,'(A,2G20.10)') 'maximum energy and mean=',MAXE,SUME/(NSEQ*GSTEPS-1)
WRITE(*,'(A,I6)') 'maximum is at location ',JMAXE(1)
PRINT*,'order of coordinate steps, energies:'
WRITE(*,'(2I6,G20.10)') (J1,STEPORDER(J1),INTERPENERGY(J1),J1=1,NSEQ*GSTEPS)
! dump xyz path
CALL GENERGIES(1,NSEQ*GSTEPS,.TRUE.) 

DEALLOCATE(STEPORDER,INTERPENERGY,INTERPENERGYSAVE,SHORTER)
IF (ALLOCATED(PEPCOORDS)) DEALLOCATE(PEPCOORDS,DELTA,SMALLSTEP,LARGESTEP)
CALL FLUSH(6)

CONTAINS
   SUBROUTINE GENERGIES(BOTTOM,TOP,DUMPPATH)
   INTEGER, INTENT(IN) :: BOTTOM, TOP
   INTEGER :: K1, K2, NDUMMY, LBOTTOM, LTOP, K3
   LOGICAL, INTENT(IN) :: DUMPPATH
   DOUBLE PRECISION :: X(3*NATOMS), RMS, GRAD(3*NATOMS)
               ! there should be a way of not 
               ! declaring the Hessian unless it is actually needed? Use (*,*) in potential
               ! and routines called from it?
   DOUBLE PRECISION :: DPRAND, ENERGY


   IF (DUMPPATH) THEN
      IF (FILTH.EQ.0) THEN
         OPEN(UNIT=7,FILE='guess.xyz',STATUS='UNKNOWN')
      ELSE
         LFNAME='guess.xyz.'//TRIM(ADJUSTL(FILTHSTR))
         OPEN(UNIT=7,FILE=TRIM(ADJUSTL(LFNAME)),STATUS='UNKNOWN')
      ENDIF
      IF (UNRST) THEN
         IF (FILTH.EQ.0) THEN
            OPEN(UNIT=8,FILE='guess.unres.xyz',STATUS='UNKNOWN')
         ELSE
            LFNAME='guess.unres.xyz.'//TRIM(ADJUSTL(FILTHSTR))
            OPEN(UNIT=8,FILE=TRIM(ADJUSTL(LFNAME)),STATUS='UNKNOWN')
         ENDIF
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
      ELSE
         WRITE(7,'(I6)') NATOMS
         WRITE(7,'(A)') ' '
         WRITE(7,'(A2,3G20.10)') ('LA  ',START(3*K2-2),START(3*K2-1),START(3*K2),K2=1,NCOORDS/3)
      ENDIF
   ENDIF
   LBOTTOM=BOTTOM
   LTOP=TOP
   IF (LBOTTOM.GT.LTOP) THEN
      NDUMMY=LBOTTOM
      LBOTTOM=LTOP
      LTOP=NDUMMY
   ENDIF
   X(1:NCOORDS)=START(1:NCOORDS)
! Do the first LBOTTOM steps - the order only changes between LBOTTOM and LTOP
   IF (UNRST) THEN
      DO K1=1,LBOTTOM-1
         X(STEPORDER(K1))=X(STEPORDER(K1))+DELTA(STEPORDER(K1))/GSTEPS
      ENDDO
   ELSE
      DO K1=1,LBOTTOM-1
         X(STEPORDER(K1))=X(STEPORDER(K1))+(FINISH(STEPORDER(K1))-START(STEPORDER(K1)))/GSTEPS
      ENDDO
   ENDIF
! Remove exact atom superpositions using some random numerical noise with DPRAND
   DO K1=LBOTTOM,LTOP
      IF (DOSEQ(STEPORDER(K1))) THEN
         IF (UNRST) THEN
            X(STEPORDER(K1))=X(STEPORDER(K1))+DELTA(STEPORDER(K1))/GSTEPS
         ELSE
            X(STEPORDER(K1))=X(STEPORDER(K1))+ &
               (FINISH(STEPORDER(K1))-START(STEPORDER(K1)))/GSTEPS+(DPRAND()-1.0D0)*1.0D-20
         ENDIF
      ENDIF
      DO K2=1,NCOORDS ! linear interpolation
         IF (.NOT.SEQUENTIAL(K2).OR.(.NOT.DOSEQ(K2))) THEN
            IF (UNRST) THEN
               X(K2)=START(K2)+DELTA(K2)*K1/(NSEQ*GSTEPS)
            ELSE
               X(K2)=(START(K2)*(NSEQ*GSTEPS-K1)+FINISH(K2)*K1)/(NSEQ*GSTEPS) +(DPRAND()-1.0D0)*1.0D-20
            ENDIF
         ENDIF
      ENDDO
      IF (UNRST) THEN
         CALL var_to_geom(NCOORDS,X)
         CALL chainbuild
      ENDIF
      IF (DUMPPATH) THEN
         IF (UNRST) THEN
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
            WRITE(7,'(I6)') NATOMS
            WRITE(7,'(A)') ' '
            WRITE(7,'(A2,3G20.10)') ('LA  ',X(3*K2-2),X(3*K2-1),X(3*K2),K2=1,NCOORDS/3)
         ENDIF
      ELSE
         CALL POTENTIAL(X,ENERGY,GRAD,.FALSE.,.FALSE.,RMS,.FALSE.,.FALSE.)
         INTERPENERGY(K1)=ENERGY
      ENDIF
   ENDDO
   IF (DUMPPATH) CLOSE(7)
   IF (DUMPPATH.AND.UNRST) CLOSE(8)

   END SUBROUTINE GENERGIES

END SUBROUTINE GUESSPATH
