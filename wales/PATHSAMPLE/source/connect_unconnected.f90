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
!  CONNECTUNCONNECTED subroutine
!
SUBROUTINE GETCUPAIR(NAVAIL,NUSED,MINS,MINF,SPOINTS,FPOINTS)
USE COMMONS, ONLY: UMIN, NATOMS, DMIN1, DMIN2,  NMIN,  PAIR1, PAIR2, NPAIRDONE,MAXPAIRS
IMPLICIT NONE
INTEGER NUSED, MINS, MINF, NAVAIL, PAIRSTODO
DOUBLE PRECISION SPOINTS(3*NATOMS), FPOINTS(3*NATOMS)

10 CONTINUE
IF (NAVAIL.EQ.0) THEN
   CALL GETNCONN
   CALL CONNECTUNCONNECTED(NAVAIL)
   IF (NAVAIL.EQ.0) THEN
      PRINT '(A)','getcupair> No more candidate pairs of minima in getdpair - quit'
      STOP
   ENDIF
   NUSED=0
ENDIF
NUSED=NUSED+1
NAVAIL=NAVAIL-1
MINS=DMIN1(NUSED)
MINF=DMIN2(NUSED)
WRITE(*,'(4(A,I8))') 'getcupair> connecting minima ',MINS,' and ',MINF, ' pairs used=',NUSED,' remaining=',NAVAIL
NPAIRDONE=NPAIRDONE+1
IF (NPAIRDONE.GT.MAXPAIRS) CALL PAIRDOUBLE
PAIR1(NPAIRDONE)=DMIN1(NUSED)
PAIR2(NPAIRDONE)=DMIN2(NUSED)
READ(UMIN,REC=MINS) SPOINTS(1:3*NATOMS)
READ(UMIN,REC=MINF) FPOINTS(1:3*NATOMS)

END SUBROUTINE GETCUPAIR


SUBROUTINE CONNECTUNCONNECTED(NAVAIL)
USE PORFUNCS
USE COMMONS
!USE COMMONS, ONLY: UMIN, NMIN, NMINA, NMINB, NCPU, DMIN1, DMIN2, CONNECTLOWESTT, NATOMS, &
!&                  CONNECTETHRESHT, CONNECTDTHRESHT, REFMIN,  DEBUG, BOXLX, BOXLY, BOXLZ, &
!  &                BULKT, TWOD, RIGIDBODY, EMIN, NTS, LOCATIONA, LOCATIONB, NCONNMAX, NATT
IMPLICIT NONE

INTEGER :: J1, J2, J3, J4, J5, J6, J7, J8, NPUSED, CLOSESTMIN
LOGICAL :: ISA(NMIN), ISB(NMIN), CHANGED, DEADTS(NTS), UNCONN(NMIN), DISTCALC
INTEGER NOFF(NMIN)
INTEGER :: NCOL(NMIN), NVAL(2*NTS), NDISTA(NMIN), NDISTB(NMIN), NCYCLE, DMIN
INTEGER :: DMAX, NUNCONA, NUNCONB, NDEAD, NUNCONN, PAIRSTODO, NAVAIL
INTEGER, ALLOCATABLE :: UNCMIN(:),CONMIN(:)
DOUBLE PRECISION, ALLOCATABLE :: EUNCON(:), UNCONDIST(:), EUNDIFF(:), DISTCON(:)
DOUBLE PRECISION :: DMATMC(2*NTS), KSUM(NMIN), CUT_UNDERFLOW, LOWESTX(3*NATOMS),SETX(3*NATOMS)
DOUBLE PRECISION :: DISTMIN, EREF, DISTANCE, DIST2, RMAT(3,3), REFX(3*NATOMS)

!
! Set up connectivity analysis
!
UNCONN(1:NMIN)=.FALSE.
ISA(1:NMIN)=.FALSE.
ISB(1:NMIN)=.FALSE.
DO J1=1,NMINA
   ISA(LOCATIONA(J1))=.TRUE.
ENDDO
DO J1=1,NMINB
   ISB(LOCATIONB(J1))=.TRUE.
ENDDO


CUT_UNDERFLOW=-300.0D0
CALL RATECONST_SETUP(KSUM,DEADTS,NDEAD,.FALSE.,CUT_UNDERFLOW)


CALL MAKEDLIN(DMATMC,NCOL,NOFF,NVAL,DEADTS,.TRUE.,ISA,ISB,KSUM)


!
! Check connections to A set
!
NDISTA(1:NMIN)=1000000
DO J1=1,NMINA
   NDISTA(LOCATIONA(J1))=0
ENDDO 
NCYCLE=0
5    CHANGED=.FALSE.
NCYCLE=NCYCLE+1
DMIN=100000
DMAX=0
NUNCONA=0
DO J1=1,NMIN
   IF (NDISTA(J1).EQ.0) CYCLE ! A minimum
   DO J2=1,NCOL(J1)
      IF (NDISTA(NVAL(NOFF(J1)+J2))+1.LT.NDISTA(J1)) THEN
         CHANGED=.TRUE.
         NDISTA(J1)=NDISTA(NVAL(NOFF(J1)+J2))+1
      ENDIF
   ENDDO
   IF ((NDISTA(J1).GT.DMAX).AND.(NDISTA(J1).NE.1000000)) DMAX=NDISTA(J1)
   IF (NDISTA(J1).LT.DMIN) DMIN=NDISTA(J1)
   IF (NDISTA(J1).EQ.1000000) NUNCONA=NUNCONA+1
ENDDO
IF (CHANGED) GOTO 5
PRINT '(3(A,I8))','Dijkstra> steps to A region converged in ',NCYCLE-1, &
&                    ' cycles; maximum=',DMAX,' disconnected=',NUNCONA
!
!  Check connections to B set
!  
NDISTB(1:NMIN)=1000000
DO J1=1,NMINB
   NDISTB(LOCATIONB(J1))=0
ENDDO
NCYCLE=0
51    CHANGED=.FALSE.
NCYCLE=NCYCLE+1
DMIN=100000
DMAX=0
NUNCONB=0
DO J1=1,NMIN
   IF (NDISTB(J1).EQ.0) CYCLE ! B minimum
   DO J2=1,NCOL(J1)
      IF (NDISTB(NVAL(NOFF(J1)+J2))+1.LT.NDISTB(J1)) THEN
         CHANGED=.TRUE.
         NDISTB(J1)=NDISTB(NVAL(NOFF(J1)+J2))+1
      ENDIF
   ENDDO
   IF ((NDISTB(J1).GT.DMAX).AND.(NDISTB(J1).NE.1000000)) DMAX=NDISTB(J1)
   IF (NDISTB(J1).LT.DMIN) DMIN=NDISTB(J1)
   IF (NDISTB(J1).EQ.1000000) NUNCONB=NUNCONB+1
ENDDO 
IF (CHANGED) GOTO 51
PRINT '(3(A,I8))','Dijkstra> steps to B region converged in ',NCYCLE-1, &
&                    ' cycles; maximum=',DMAX,' disconnected=',NUNCONB
!
!  This could happen if disconnected minima lie in the A or B region
!
IF (NUNCONB.NE.NUNCONA) PRINT '(A)', &
&                   'Dijkstra> WARNING - number of disconnected minima from A and B is different'
!
!  Get logical array for connection to A or B sets (True means unconnected)
!

NUNCONN=0
DO J1=1,NMIN
   IF ((NDISTA(J1).EQ.1000000).OR.(NDISTB(J1).EQ.1000000)) THEN
      UNCONN(J1)=.TRUE.
      NUNCONN=NUNCONN+1
   ENDIF
ENDDO
PRINT '(A,I8)','connectunconnected> Number of unconnected minima:',NUNCONN

ALLOCATE(UNCMIN(NUNCONN)) !get an unsorted list of all unconnected minima
UNCMIN(1:NUNCONN)=0
ALLOCATE(CONMIN(NMIN-NUNCONN))   !get list of minima in connected set
CONMIN(1:(NMIN-NUNCONN))=0
ALLOCATE(EUNCON(NUNCONN))
J2=1
J3=1
DO J1=1,NMIN
   IF (UNCONN(J1)) THEN
      EUNCON(J2)=EMIN(J1)               !assign energy of minimum, minus as we want sorting lowest->highest
      UNCMIN(J2)=J1                      !assign minimum id (for sorting)
      J2=J2+1
   ELSE
      CONMIN(J3)=J1                      !similarly create list of minima in AB set
      J3=J3+1
   ENDIF
ENDDO

PAIRSTODO=NCPU                                     !recalculate for every cycle
IF (ALLOCATED(DMIN1)) DEALLOCATE(DMIN1,DMIN2)      
ALLOCATE(DMIN1(PAIRSTODO),DMIN2(PAIRSTODO))
IF (CONNECTLOWESTT) THEN   !connect lowest energy minima to closest from connected AB set
   !need to sort list by energies, lowest to highest
   !CALL SORT(NUNCONN,NUNCONN,EUNCON,UNCMIN)
   CALL DSORT(EUNCON,UNCMIN,NUNCONN,2)
   IF (DEBUG) THEN
      PRINT '(A)','connectunconnected> Sorted energies of unconnected minima'
      DO J7 = 1,NUNCONN,1
         PRINT '(2(A,I8),F20.10)' , 'Minimum',J7,':',UNCMIN(J7),EUNCON(J7)
      ENDDO
   ENDIF
!   !use closest distance and attempt one connection for each minimum
!   J7=0
!   DO J4=1,PAIRSTODO
!54    J7=J7+1
!      READ(UMIN,REC=UNCMIN(J7)) (LOWESTX(J6),J6=1,3*NATOMS)
!      DISTMIN = 1.0D100
!      CLOSESTMIN = 0
!      DO J5=1,(NMIN-NUNCONN),1
!         READ(UMIN,REC=CONMIN(J5)) (SETX(J6),J6=1,3*NATOMS)
!         CALL MINPERMDIST(SETX,LOWESTX,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE, &
!&                         DIST2,RIGIDBODY,RMAT,.FALSE.)
!         IF (DISTANCE.LT.DISTMIN) THEN
!            DISTMIN = DISTANCE
!            !IF (DEBUG) PRINT '(A,I8,A,F15.8)','connectlowest> New closest minimum:',CONMIN(J5),' distance:',DISTMIN
!            CLOSESTMIN = CONMIN(J5)
!         ENDIF
!      ENDDO
!      DO J3=1,NPAIRDONE,1
!         IF ((PAIR1(J3).EQ.UNCMIN(J7)).AND.(PAIR2(J3).EQ.CLOSESTMIN)) THEN
!            GOTO 54
!         ELSE IF ((PAIR1(J3).EQ.CLOSESTMIN).AND.(PAIR2(J3).EQ.UNCMIN(J7))) THEN
!            GOTO 54
!         ENDIF
!      ENDDO
!      DMIN1(J4) = UNCMIN(J7)
!      DMIN2(J4) = CLOSESTMIN
!      PRINT '(A,I8,A,F15.8)','connectlowest> Minimum 1:',DMIN1(J4),' Energy:',EMIN(DMIN1(J4))
!      PRINT '(A,I8,A,F15.8)','connectlowest> Minimum 2:',DMIN2(J4),' Energy:',EMIN(DMIN2(J4))
!      PRINT '(A,F15.8)','connectlowest> Distance:',DISTMIN
!   ENDDO

   !find n closest minima starting from the lowest energy minimum
   ALLOCATE(DISTCON(NMIN-NUNCONN))   !get list of distances to reference
   DISTCON(1:(NMIN-NUNCONN))=0.0     !initialise to 0.0
   DISTCALC = .FALSE.                !new minimum means we need to calculate the distances 
   J7 = 0                            !minimum under consideration
   DO J4=1,PAIRSTODO                 
54    IF (.NOT.DISTCALC) THEN        !only recalculate distance if we have a new minimum!
         J7 = J7 + 1                 !next minimum
         REFMIN = UNCMIN(J7)         !assign new reference minimum (starting from the lowest energy minimum)
         READ(UMIN,REC=UNCMIN(J7)) (LOWESTX(J6),J6=1,3*NATOMS)
         DO J5=1,(NMIN-NUNCONN),1   
            READ(UMIN,REC=CONMIN(J5)) (SETX(J6),J6=1,3*NATOMS)
            CALL MINPERMDIST(SETX,LOWESTX,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE, &
&                         DIST2,RIGIDBODY,RMAT,.FALSE.) 
            DISTCON(J5) = DISTANCE  !again use negative as we sort lowest->highest
            IF (DEBUG) PRINT '(A,I8,A,F15.8)', 'connectlowest> Minimum:',CONMIN(J5),' Distance:',DISTANCE
         ENDDO
         !CALL SORT((NMIN-NUNCONN),(NMIN-NUNCONN),DISTCON,CONMIN) !sort AB set by distance
         CALL DSORT(DISTCON,CONMIN,(NMIN-NUNCONN),2)
         IF (DEBUG) PRINT '(A)', 'connectlowest> Sorted distances'
         DISTCALC = .TRUE.  !have valid distances now
         NPUSED = 0 !numbers of connections used, reset for each new reference minmum
         J8 = 0 !actual minimum used from CONMIN (need to consider previously used pairs)
      ENDIF
      IF (NPUSED.LT.NATT) THEN    !if the number of pairs for the current reference is less than the user specified
57       J8 = J8 + 1              !next minimum in CONMIN
         DO J3=1,NPAIRDONE,1   
            IF ((PAIR1(J3).EQ.REFMIN).AND.(PAIR2(J3).EQ.CONMIN(J8))) THEN
               GOTO 57
            ELSE IF ((PAIR1(J3).EQ.CONMIN(J8)).AND.(PAIR2(J3).EQ.REFMIN)) THEN
               GOTO 57
            ENDIF
         ENDDO
         NPUSED = NPUSED + 1     !we have found a new pair
         DMIN1(J4) = REFMIN      !assign them accordingly 
         DMIN2(J4) = CONMIN(J8)
         PRINT '(A,I8,A,F15.8)','connectlowest> Minimum 1:',DMIN1(J4),' Energy:',EMIN(DMIN1(J4))
         PRINT '(A,I8,A,F15.8)','connectlowest> Minimum 2:',DMIN2(J4),' Energy:',EMIN(DMIN2(J4))
         PRINT '(A,F15.8)','connectlowest> Distance:',DISTCON(J8)
      ELSE                      !reached the number previously specified, go for the next lowest minimum in UNCONMIN
         DISTCALC = .FALSE.     !reset DISTCALC and return to 54 to reinitialise everything (as we assign a new distance
         GOTO 54                !to each minimum and resort, we dont need to reinitialise the distance list
      ENDIF
   ENDDO 

ELSEIF (CONNECTETHRESHT) THEN  !connect states closest in energy to a given min to set
   !need to sort list by relative energies
   EREF=EMIN(REFMIN)        
   ALLOCATE(EUNDIFF(NUNCONN))
   DO J1=1,NUNCONN,1
      EUNDIFF(J1)=ABS(EUNCON(J1)-EREF)  !use absolute difference
   ENDDO
   !CALL SORT(NUNCONN,NUNCONN,EUNDIFF,UNCMIN)  !sort UNCMIN by the energy differences
   CALL DSORT(EUNDIFF,UNCMIN,NUNCONN,2)
   !now find the pairs to connect
   J7=0
   DO J4=1,PAIRSTODO
55    IF (J7.LT.NUNCONN) THEN 
         J7=J7+1
      ELSE
         PRINT '(A)', 'connectclosestE> exhausted the list of unconnected minima for this reference'
         STOP
      ENDIF
      DO J3=1,NPAIRDONE,1
         IF ((PAIR1(J3).EQ.UNCMIN(J7)).AND.(PAIR2(J3).EQ.REFMIN)) THEN
            GOTO 55
         ELSE IF ((PAIR1(J3).EQ.REFMIN).AND.(PAIR2(J3).EQ.UNCMIN(J7))) THEN
            GOTO 55
         ENDIF
      ENDDO
      DMIN1(J4)=UNCMIN(J7)
      DMIN2(J4)=REFMIN
      PRINT '(A,I8,A,F15.8)','connectclosestE> Minimum 1:',DMIN1(J4),' Energy:',EMIN(DMIN1(J4))
      PRINT '(A,I8,A,F15.8)','connectclosestE> Minimum 2:',DMIN2(J4),' Energy:',EMIN(DMIN2(J4))
      PRINT '(A,F15.8)','connectlowestE> Energy difference:',EUNDIFF(J7)
   ENDDO

ELSEIF (CONNECTDTHRESHT) THEN  !connect states closest in energy to a given min to set
   !need to compute distances for all unconnected minima and sort it
   READ(UMIN,REC=REFMIN) (REFX(J6),J6=1,3*NATOMS)
   ALLOCATE(UNCONDIST(NUNCONN))
   DO J5=1,NUNCONN,1
      READ(UMIN,REC=UNCMIN(J5)) (SETX(J6),J6=1,3*NATOMS)
      CALL MINPERMDIST(SETX,REFX,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE, &
&                         DIST2,RIGIDBODY,RMAT,.FALSE.)
      UNCONDIST(J5) = DISTANCE                           
   ENDDO
  !CALL SORT(NUNCONN,NUNCONN,UNCONDIST,UNCMIN)
   CALL DSORT(UNCONDIST,UNCMIN,NUNCONN,2)
   !now find the pairs to connect
   J7=0
   DO J4=1,PAIRSTODO
56    IF (J7.LT.NUNCONN) THEN
         J7=J7+1
      ELSE
         PRINT '(A)', 'connectclosestE> exhausted the list of unconnected minima for this reference'
         STOP
      ENDIF
      DO J3=1,NPAIRDONE,1
         IF ((PAIR1(J3).EQ.UNCMIN(J7)).AND.(PAIR2(J3).EQ.REFMIN)) THEN
            GOTO 56
         ELSE IF ((PAIR1(J3).EQ.REFMIN).AND.(PAIR2(J3).EQ.UNCMIN(J7))) THEN
            GOTO 56
         ENDIF
      ENDDO
      DMIN1(J4)=UNCMIN(J7)
      DMIN2(J4)=REFMIN
      PRINT '(A,I8,A,F15.8)','connectclosestD> Minimum 1:',DMIN1(J4),' Energy:',EMIN(DMIN1(J4))
      PRINT '(A,I8,A,F15.8)','connectclosestD> Minimum 2:',DMIN2(J4),' Energy:',EMIN(DMIN2(J4))
      PRINT '(A,F15.8)','connectclosestD> Distance:', UNCONDIST(J7)
   ENDDO
ENDIF

NAVAIL = PAIRSTODO
DEALLOCATE(UNCMIN,CONMIN,EUNCON)
IF (ALLOCATED(UNCONDIST)) DEALLOCATE(UNCONDIST)
IF (ALLOCATED(EUNDIFF)) DEALLOCATE(EUNDIFF)
IF (ALLOCATED(DISTCON)) DEALLOCATE(DISTCON)
END SUBROUTINE CONNECTUNCONNECTED


