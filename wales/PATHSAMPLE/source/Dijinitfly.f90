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
!  Dijkstra connection algorithm for pathsample.
!  This version calculates the weights for the edges on the fly for
!  cases where there are too many minima to hold the distance matrix in memory.
!
SUBROUTINE DIJINITFLY(NWORST)
USE PORFUNCS
USE COMMONS
IMPLICIT NONE

INTEGER J1, J2, J4, PARENT(NMIN), JMINW, NPERM, J5, LJ1, LJ2, NWORST, NSTEPS, NMINSTART, NMINEND, J3, LNCONN, NPREV
INTEGER PREVCONN(NMIN)
INTEGER, ALLOCATABLE :: LOCATIONSTART(:), LOCATIONEND(:)
LOGICAL PERMANENT(NMIN), ISA(NMIN), ISB(NMIN), ISSTART(NMIN), CONNECTED, PREVTEST, NOTDONE, DEADTS
DOUBLE PRECISION MINWEIGHT, TMPWEIGHT, WEIGHT(NMIN), DUMMY, TNEW, ELAPSED, PFTOTALSTART, HUGESAVE
DOUBLE PRECISION LOCALPOINTS(3*NATOMS), LOCALPOINTS2(3*NATOMS), DISTANCE, DIST2, RMAT(3,3)
INTEGER LCONNECTION(NCONNMAX)
DOUBLE PRECISION MAXWEIGHT, SCALEFAC
DOUBLE PRECISION :: CUT_UNDERFLOW=-300.0D0

CALL CPU_TIME(ELAPSED)
DEALLOCATE(DMIN1,DMIN2)
ALLOCATE(DMIN1(10000),DMIN2(10000)) ! surely this will be enough room!
ISA(1:NMIN)=.FALSE.
ISB(1:NMIN)=.FALSE.
DO J1=1,NMINA
   ISA(LOCATIONA(J1))=.TRUE.
ENDDO
DO J1=1,NMINB
   ISB(LOCATIONB(J1))=.TRUE.
ENDDO

!!!!!!!!!!!!!!!!!!!   Dijkstra calculation    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Dijkstra connect process similar to that in OPTIM with a weight for missing
!  connections based on a distance metric.
!
IF (DIRECTION.EQ.'AB') THEN
   NMINSTART=NMINB
   NMINEND=NMINA
   ALLOCATE(LOCATIONSTART(NMINB),LOCATIONEND(NMINA))
   LOCATIONSTART(1:NMINB)=LOCATIONB(1:NMINB)
   LOCATIONEND(1:NMINA)=LOCATIONA(1:NMINA)
   ISSTART(1:NMIN)=ISB(1:NMIN)
   PFTOTALSTART=PFTOTALB
ELSEIF (DIRECTION.EQ.'BA') THEN
   NMINSTART=NMINA
   NMINEND=NMINB
   ALLOCATE(LOCATIONSTART(NMINA),LOCATIONEND(NMINB))
   LOCATIONSTART(1:NMINA)=LOCATIONA(1:NMINA)
   LOCATIONEND(1:NMINB)=LOCATIONB(1:NMINB)
   ISSTART(1:NMIN)=ISA(1:NMIN)
   PFTOTALSTART=PFTOTALA
ENDIF
!
!  Find largest weight for each B(A) minimum to all A(B) minima.
!
loopstart: DO J1=1,NMINSTART ! cycle over all minima in the starting state
   MAXWEIGHT=1.0D10
   SCALEFAC=1.0D0
222 LJ1=LOCATIONSTART(J1)
   WEIGHT(1:NMIN)=HUGE(1.0D0)
   HUGESAVE=WEIGHT(1)
   WEIGHT(LJ1)=0.0D0
   PERMANENT(1:NMIN)=.FALSE.
   PERMANENT(LJ1)=.TRUE.
   NPERM=1
   PARENT(1:NMIN)=0 ! parent is initially undefined
   J4=LJ1

   dijkstraloop: DO

      READ(UMIN,REC=J4) (LOCALPOINTS(J2),J2=1,3*NATOMS)
!
!  Save connections for minimum J4
!
      LNCONN=0
      DO J3=1,NTS 
!
! Disregard transition states that are too high or have two barriers > MAXBARRIER
! Need to match here the criteria applied in GETNCONN!
! JMC n.b. don't apply the nconnmin criteria at this point, hence the huge(1) 's in place of NCONN() for the plus and minus minima.
         CALL CHECKTS(ETS(J3),EMIN(PLUS(J3)),EMIN(MINUS(J3)),KPLUS(J3),KMINUS(J3),HUGE(1),HUGE(1), &
                      PLUS(J3),MINUS(J3),.FALSE.,CUT_UNDERFLOW,DEADTS)

         IF ((.NOT. DEADTS).AND.(PLUS(J3).NE.MINUS(J3))) THEN
            IF (PLUS(J3).EQ.J4) THEN
               LNCONN=LNCONN+1
               LCONNECTION(LNCONN)=MINUS(J3)
            ELSE IF (MINUS(J3).EQ.J4) THEN
               LNCONN=LNCONN+1
               LCONNECTION(LNCONN)=PLUS(J3)
            ENDIF
         ENDIF
      ENDDO
      IF (LNCONN.NE.NCONN(J4)) THEN
         PRINT '(A,3I8)',' Dijinitfly> ERROR - J4,LNCONN,NCONN=',J4,LNCONN,NCONN(J4)
         STOP
      ENDIF
!
!  Save previous connection attempts for minimum J4
!
      NPREV=0
      DO J3=1,NPAIRDONE
         IF (PAIR1(J3).EQ.J4) THEN
            NPREV=NPREV+1
            PREVCONN(NPREV)=PAIR2(J3)
         ELSE IF (PAIR2(J3).EQ.J4) THEN
            NPREV=NPREV+1
            PREVCONN(NPREV)=PAIR1(J3)
         ENDIF
      ENDDO 
      IF (NPREV.GE.NMIN-1) THEN
         PRINT '(A,I6)','Dijinitfly> WARNING *** All connections have been tried for parent minimum ',J4
         STOP
      ENDIF

      DO J2=1,NMIN
         IF (J2.EQ.J4) CYCLE
         IF (PERMANENT(J2)) CYCLE
         CONNECTED=.FALSE.
         DO J3=1,LNCONN
            IF (LCONNECTION(J3).EQ.J2) THEN
               TMPWEIGHT=0.0D0
               CONNECTED=.TRUE.
               EXIT
            ENDIF
         ENDDO
         PREVTEST=.FALSE.
         IF (.NOT.CONNECTED) THEN
            DO J3=1,NPREV
               IF (PREVCONN(J3).EQ.J2) THEN
                  TMPWEIGHT=1.0D100
                  PREVTEST=.TRUE.
                  EXIT
               ENDIF
            ENDDO
         ENDIF
         IF ((.NOT.CONNECTED).AND.(.NOT.PREVTEST)) THEN
            IF (INDEXCOSTFUNCTION) THEN
               TMPWEIGHT=1.0D0
            ELSE
               READ(UMIN,REC=J2) (LOCALPOINTS2(J3),J3=1,3*NATOMS)
!
! Calculate distances on the fly
!
               CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, &
  &                             RMAT,.FALSE.)
               IF (INTERPCOSTFUNCTION) CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD, &
  &                                                     DISTANCE,DIST2,RIGIDBODY,RMAT,INTERPCOSTFUNCTION)
               TMPWEIGHT=DISTANCE*SCALEFAC
            ENDIF
         ENDIF
         IF (TMPWEIGHT.LT.HUGE(1.0D0)/10.0D0) THEN ! don;t raise a huge number to any power!
            IF (INDEXCOSTFUNCTION) THEN
               IF (CONNECTED) THEN
                  TMPWEIGHT=0.0D0 
               ELSE
                  TMPWEIGHT=ABS(J4-J2)
                  IF (DIRECTION.EQ.'BA') THEN
                     IF (J4.LE.NMINA) TMPWEIGHT=NMIN+1-J2
                     IF (J2.LE.NMINA) TMPWEIGHT=NMIN+1-J4
                  ELSE
                     IF ((J4.LE.NMINA+NMINB).AND.(J4.GT.NMINA)) TMPWEIGHT=NMIN+1-J2
                     IF ((J2.LE.NMINA+NMINB).AND.(J2.GT.NMINA)) TMPWEIGHT=NMIN+1-J4
                  ENDIF
               ENDIF
            ELSEIF (EXPCOSTFUNCTION) THEN ! saves memory and CPU when endpoint separation is very large SAT
               IF (CONNECTED) THEN
                  TMPWEIGHT=0.0D0 ! not one !! DJW 22/7/08
               ELSEIF (TMPWEIGHT.GT.300.0D0) THEN
                  TMPWEIGHT=EXP(300.0D0)
               ELSE
                  TMPWEIGHT=EXP(TMPWEIGHT)
               ENDIF
            ELSE ! compare squares to favour more small jumps over big ones DJW
               IF (CONNECTED) THEN
                  TMPWEIGHT=0.0D0 
               ELSEIF (TMPWEIGHT.EQ.0.0D0) THEN
               ELSEIF (INTERPCOSTFUNCTION) THEN
                  TMPWEIGHT=TMPWEIGHT**COSTFUNCTIONPOWER 
               ELSEIF (COSTFUNCTIONPOWER.EQ.0) THEN
                  TMPWEIGHT=TMPWEIGHT+1.0D0
               ELSEIF (COSTFUNCTIONPOWER.EQ.-1) THEN
                  IF (TMPWEIGHT.GT.0.0D0) TMPWEIGHT=1.0D0/TMPWEIGHT
               ELSE
                  TMPWEIGHT=TMPWEIGHT**COSTFUNCTIONPOWER 
               ENDIF
            ENDIF
         ENDIF
         IF (TMPWEIGHT+WEIGHT(J4).LT.WEIGHT(J2)) THEN ! relax J2
            WEIGHT(J2)=TMPWEIGHT+WEIGHT(J4)
            PARENT(J2)=J4
         ENDIF
      ENDDO

      MINWEIGHT=HUGE(1.0D0)
      NOTDONE=.TRUE.
      DO J2=1,NMIN
         IF (.NOT.PERMANENT(J2)) THEN
            IF (WEIGHT(J2).LT.MINWEIGHT) THEN
               MINWEIGHT=WEIGHT(J2)
               JMINW=J2
               NOTDONE=.FALSE.
            ENDIF
         ENDIF
      ENDDO
      IF (NOTDONE) THEN
         PRINT '(A,I8,A,I8)','Dijinitfly> WARNING - JMINW not set - value=',JMINW,' J4=',J4
      ENDIF

      J4=JMINW
      PERMANENT(J4)=.TRUE.
      NPERM=NPERM+1
      IF (WEIGHT(J4).GT.MAXWEIGHT) THEN
         SCALEFAC=SCALEFAC/2.0D0
         PRINT '(A,G20.10)','Dijinitfly> Maximum weight is too large - scaling distances by ',SCALEFAC
         GOTO 222
      ENDIF

      IF (NPERM.EQ.NMIN) EXIT dijkstraloop

   ENDDO dijkstraloop

ENDDO loopstart
! 
!  Summarise the best path for any A(B) and any B(A)
!
LJ2=LOCATIONEND(1)
LJ1=LOCATIONSTART(1)
J5=LJ2
NWORST=0
NSTEPS=0
PRINT '(A)','Dijinit> Summary of best path based on missing connection metric - note distance scaling is removed'
PRINT '(A)','    min1          energy        min2          energy             metric          edge weight            weight'
DO 
   IF (PARENT(J5).EQ.0) THEN
      PRINT '(A,I6,A)','Dijinitfly> ERROR - parent for J5=',J5,' is zero'
      PRINT '(A)',     'Dijinitfly> Suggests all possible pairs have been tried!'
      STOP
   ENDIF
!
! Calculate distances on the fly
!
   CONNECTED=.FALSE.
   DO J3=1,NTS
      IF (((PLUS(J3).EQ.J5).AND.(MINUS(J3).EQ.PARENT(J5))).OR. &
      &   ((MINUS(J3).EQ.J5).AND.(PLUS(J3).EQ.PARENT(J5)))) THEN
         DUMMY=0.0D0
         CONNECTED=.TRUE.
         EXIT
      ENDIF
   ENDDO

   IF (.NOT.CONNECTED) THEN
      IF (INDEXCOSTFUNCTION) THEN
         DUMMY=1.0D0
      ELSE
         READ(UMIN,REC=J5) (LOCALPOINTS(J2),J2=1,3*NATOMS)
         READ(UMIN,REC=PARENT(J5)) (LOCALPOINTS2(J2),J2=1,3*NATOMS)
         CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY, &
  &                       RMAT,.FALSE.)
         IF (INTERPCOSTFUNCTION) CALL MINPERMDIST(LOCALPOINTS,LOCALPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD, &
  &                                               DISTANCE,DIST2,RIGIDBODY,RMAT,INTERPCOSTFUNCTION)
         DUMMY=DISTANCE*SCALEFAC
      ENDIF
   ENDIF
   IF (DUMMY.LT.HUGE(1.0D0)/10.0D0) THEN ! don;t raise a huge number to any power!
      IF (INDEXCOSTFUNCTION) THEN
         IF (DUMMY.EQ.0.0D0) THEN
            TMPWEIGHT=0.0D0
         ELSE
            TMPWEIGHT=ABS(J5-PARENT(J5))
            IF (DIRECTION.EQ.'AB') THEN
               IF (J5.LE.NMINA) TMPWEIGHT=NMIN+1-PARENT(J5)
               IF (PARENT(J5).LE.NMINA) TMPWEIGHT=NMIN+1-J5
            ELSE
               IF ((J5.LE.NMINA+NMINB).AND.(J5.GT.NMINA)) TMPWEIGHT=NMIN+1-PARENT(J5)
               IF ((PARENT(J5).LE.NMINA+NMINB).AND.(PARENT(J5).GT.NMINA)) TMPWEIGHT=NMIN+1-J5
            ENDIF
         ENDIF
      ELSEIF (EXPCOSTFUNCTION) THEN ! saves memory and CPU when endpoint separation is very large SAT
         IF (DUMMY.EQ.0.0D0) THEN
            TMPWEIGHT=0.0D0
         ELSEIF (DUMMY.GT.300.0D0) THEN
             TMPWEIGHT=EXP(300.0D0)
         ELSE  
            TMPWEIGHT=EXP(DUMMY)
         ENDIF 
      ELSE ! compare squares to favour more small jumps over big ones DJW
         IF (DUMMY.EQ.0.0D0) THEN
            TMPWEIGHT=0.0D0
         ELSEIF (INTERPCOSTFUNCTION) THEN
            TMPWEIGHT=DUMMY**COSTFUNCTIONPOWER
         ELSEIF (COSTFUNCTIONPOWER.EQ.0) THEN
            TMPWEIGHT=1.0D0
         ELSEIF (COSTFUNCTIONPOWER.EQ.-1) THEN
            TMPWEIGHT=1.0D0/DUMMY
         ELSE
            TMPWEIGHT=DUMMY**COSTFUNCTIONPOWER
         ENDIF
      ENDIF
   ELSE
      TMPWEIGHT=DUMMY
   ENDIF
   NSTEPS=NSTEPS+1
   
   PRINT '(2(I8,G20.10),3G20.10)',J5,EMIN(J5),parent(J5),EMIN(PARENT(J5)),DUMMY/SCALEFAC,TMPWEIGHT,WEIGHT(J5)
   IF ((DUMMY.GT.0.0D0).AND.(TMPWEIGHT.LT.HUGE(1.0D0)/10.0D0)) THEN
      NWORST=NWORST+1
      IF (NWORST.GT.10000) THEN
         PRINT '(A,I8)','ERROR in Dijinit, too many gaps, NWORST=',NWORST
         STOP
      ENDIF
      DMIN1(NWORST)=J5
      DMIN2(NWORST)=PARENT(J5)
!
!  On the first attempt there may only be two minima, and with more than one cpu available
!  we may need to choose the same best path several times before we can get
!  more than one distinct candidate.
!
!     IF (NMIN.GT.2) PAIRDIST(MAX(J5,PARENT(J5))*(MAX(J5,PARENT(J5))-1)/2+MIN(J5,PARENT(J5)))=HUGE(1.0D0)
   ENDIF
   J5=PARENT(J5)
   IF (J5.EQ.LJ1) EXIT
   IF (J5.EQ.0) EXIT
ENDDO
PRINT '(2(A,I8))','Dijinitfly> Number of steps=',NSTEPS,' number of missing connections=',NWORST
IF (NWORST.EQ.0) THEN
   PRINT '(A)','Dijinitfly> Connected path found'
   STOP
ENDIF
CALL CPU_TIME(TNEW)
TDIJKSTRA=TDIJKSTRA+TNEW-ELAPSED

RETURN
END
