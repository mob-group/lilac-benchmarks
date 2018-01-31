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
!
SUBROUTINE DIJINIT(NWORST)
USE PORFUNCS
USE COMMONS
USE UTILS,ONLY : GETUNIT
IMPLICIT NONE

INTEGER J1, J2, J4, PARENT(NMIN), JMINW, NPERM, J5, LJ1, LJ2, NWORST, NSTEPS, NMINSTART, NMINEND, J6, MUNIT, J7
INTEGER JN, JM, NPOSITION, LUNIT
INTEGER NMINGAP, NPRUNEDONE, NPRUNEPAIRS, NPRUNEMIN,NPRUNEPAIRSOLD, RELAXED(NMIN), NNEIGH
INTEGER, ALLOCATABLE :: LOCATIONSTART(:), LOCATIONEND(:), PRUNEPAIRS(:,:)
LOGICAL PERMANENT(NMIN), ISA(NMIN), ISB(NMIN), ISSTART(NMIN), NOTDONE, PRUNEMIN(NMIN), REDODIJKSTRA
DOUBLE PRECISION MINWEIGHT, DUMMY, TNEW, ELAPSED, PFTOTALSTART, HUGESAVE, THRESH, LPOINTS1(NOPT), LPOINTS2(NOPT)
DOUBLE PRECISION MAXWEIGHT, SCALEFAC, PDMAX, PD, MINGAPTHRESH, DIST2, DISTANCE, RMAT(3,3), MINNWEIGHT, MINNDIST
!
! KIND=16 is not supported by Portland. If you want extra precision, uncomment the following line
! and use NAG.
!
! REAL(KIND=16) :: TMPWEIGHT, WEIGHT(NMIN)
REAL(KIND=8) :: TMPWEIGHT, WEIGHT(NMIN)

IF (DIJPRUNET) THEN
   MUNIT=GETUNIT()
   IF (.NOT.(PRUNECYCLET)) NPRUNE=1
   NPRUNEDONE=0
   NPRUNEPAIRS=0
   ALLOCATE(PRUNEPAIRS(2,10000))
ENDIF
PRUNEMIN(1:NMIN)=.FALSE.

CALL CPU_TIME(ELAPSED)
121 DEALLOCATE(DMIN1,DMIN2)
ALLOCATE(DMIN1(10000),DMIN2(10000)) ! surely this will be enough room!
ISA(1:NMIN)=.FALSE.
ISB(1:NMIN)=.FALSE.
DO J1=1,NMINA
   ISA(LOCATIONA(J1))=.TRUE.
ENDDO
DO J1=1,NMINB
   ISB(LOCATIONB(J1))=.TRUE.
ENDDO
IF (DIJPRUNET) NPRUNEPAIRSOLD=NPRUNEPAIRS

!!!!!!!!!!!!!!!!!!!   Dijkstra calculation    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Dijkstra connect process similar to that in OPTIM with a weight for missing
!  connections based on a distance metric.
!
IF (DIRECTION.EQ.'AB') THEN
   NMINSTART=NMINB
   NMINEND=NMINA
   IF (ALLOCATED(LOCATIONSTART)) DEALLOCATE(LOCATIONSTART,LOCATIONEND)
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

642 CONTINUE ! return here for REDODIJKSTRA
PDMAX=-1.0D0
IF (INITIALDIST) THEN
   DO J2=1,(NMIN*(NMIN-1))/2
      IF (ABS(ALLPAIRS(J2)).GT.PDMAX) THEN
         PDMAX=ABS(ALLPAIRS(J2))
!        IF (DEBUG) PRINT '(A,G20.10)','Dijinit> maximum neighbour metric value increased to',PDMAX 
!        IF (DEBUG) PRINT '(A,I8,G20.10)','Dijinit> J2,ALLPAIRS=',J2,ALLPAIRS(J2)
      ENDIF
   ENDDO
ELSE
   DO J2=1,NMIN
      DO J5=1,PAIRDISTMAX
         IF (PAIRDIST(J2,J5).GT.PDMAX) THEN
            PDMAX=PAIRDIST(J2,J5)
            IF (DEBUG) PRINT '(A,G20.10)','Dijinit> maximum neighbour metric value increased to',PDMAX 
            IF (DEBUG) PRINT '(A,2I8,G20.10,I8)','Dijinit> J2,J5,PAIRDIST,PAIRLIST=',J2,J5,PAIRDIST(J2,J5),PAIRLIST(J2,J5)
         ENDIF
      ENDDO
   ENDDO
   PRINT '(A,G20.10)','Dijinit> maximum neighbour metric value=',PDMAX
ENDIF

IF (MINGAPT) THEN
   IF (MINGAPRATIOT) THEN
      MINGAPTHRESH=PDMAX*MINGAPINP
   ELSE
      MINGAPTHRESH=MINGAPINP
   ENDIF
ENDIF

!
!  Find largest weight for each B(A) minimum to all A(B) minima.
!
!  Added maximum weight condition via a scale factor.
!  Otherwise loss of precision can cause connections to be missed completely. DJW 29/7/08
!  PAIR1 and PAIR2 are connections from pairs.data that have already been tried
!
MAXWEIGHT=HUGE(1.0D0)/1.0D1
! MAXWEIGHT=1.0D6
loopstart: DO J1=1,NMINSTART ! cycle over all minima in the starting state
   SCALEFAC=1.0D0
222   LJ1=LOCATIONSTART(J1)
   WEIGHT(1:NMIN)=HUGE(1.0D0)
   HUGESAVE=WEIGHT(1)
   WEIGHT(LJ1)=0.0D0
   PERMANENT(1:NMIN)=.FALSE.
   PERMANENT(LJ1)=.TRUE.
   RELAXED(1:NMIN)=0
   NPERM=1
   PARENT(1:NMIN)=0 ! parent is initially undefined
   J4=LJ1
   dijkstraloop: DO
      NNEIGH=0
      MINNWEIGHT=HUGE(1.0D0)
      MINNDIST=HUGE(1.0D0)
      DO J2=1,NMIN
         IF (J2.EQ.J4) CYCLE
         IF (PERMANENT(J2)) CYCLE
         PD=1.0D4*PDMAX
         JM=MIN(J4,J2)
         JN=MAX(J4,J2)
         NPOSITION=((JN-2)*(JN-1))/2+JM
!        IF (INITIALDIST) PRINT '(A,I8,A,I8,A,I10,A,G20.10)','Dijinit> minima ',J2,' and ',J4,' position ',NPOSITION,' distance ',ALLPAIRS(NPOSITION)
         IF (.NOT.PAIRSIGNORET) THEN !for pruning the database all minima count not just the ones not searched yet
            DO J5=1,NPAIRDONE ! skip
               IF (INITIALDIST) THEN
                  IF ((PAIR1(J5).EQ.J4).AND.(PAIR2(J5).EQ.J2)) THEN 
                     PD=ABS(ALLPAIRS(NPOSITION))
                     IF (PD.EQ.0.0D0) THEN
                        GOTO 973
                     ELSE
                        PD=1.0D4*PDMAX
                        GOTO 973
                     ENDIF
                     GOTO 973
                  ENDIF
                  IF ((PAIR1(J5).EQ.J2).AND.(PAIR2(J5).EQ.J4)) THEN
                     PD=ABS(ALLPAIRS(NPOSITION))
                     IF (PD.EQ.0.0D0) THEN
                        GOTO 973
                     ELSE
                        PD=1.0D4*PDMAX
                        GOTO 973
                     ENDIF
                     GOTO 973
                  ENDIF
               ELSE
!kr366> check if pair has been searched before, if so enter DO loop to check if PAIRDIST
!is 0.0D0: If yes go to 973, else set PAIRDIST to 1.0D4*PDMAX
                  IF ((PAIR1(J5).EQ.J4).AND.(PAIR2(J5).EQ.J2)) THEN 
                     DO J6=1,PAIRDISTMAX
                        IF (PAIRLIST(J4,J6).EQ.J2) THEN
                           PD=PAIRDIST(J4,J6)
                           IF (PD.EQ.0.0D0) THEN
                              GOTO 973
                           ELSE
                              PD=1.0D4*PDMAX
                              GOTO 973
                           ENDIF
                        ENDIF
                     ENDDO
                     GOTO 973
                  ENDIF
                  IF ((PAIR1(J5).EQ.J2).AND.(PAIR2(J5).EQ.J4)) THEN
                     DO J6=1,PAIRDISTMAX
                        IF (PAIRLIST(J2,J6).EQ.J4) THEN
                           PD=PAIRDIST(J2,J6)
                           IF (PD.EQ.0.0D0) THEN
                              GOTO 973
                           ELSE
                              PD=1.0D4*PDMAX
                              GOTO 973
                           ENDIF
                        ENDIF
                     ENDDO
                     GOTO 973
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
         IF (PRUNECYCLET.AND.(NPRUNEPAIRSOLD.GT.0)) THEN
            DO J5=1,NPRUNEPAIRS
               IF ((PRUNEPAIRS(1,J5).EQ.J4).AND.(PRUNEPAIRS(2,J5).EQ.J2)) THEN
                  IF (DEBUG) PRINT '(A,2I8)','pruning> pair used: ',PRUNEPAIRS(1,J5),PRUNEPAIRS(2,J5)
                  PD=1.0D4*PDMAX
                  GOTO 973
               ENDIF
               IF ((PRUNEPAIRS(2,J5).EQ.J4).AND.(PRUNEPAIRS(1,J5).EQ.J2)) THEN
                  IF (DEBUG) PRINT '(A,2I8)','pruning> pair used: ',PRUNEPAIRS(1,J5),PRUNEPAIRS(2,J5)
                  PD=1.0D4*PDMAX
                  GOTO 973
               ENDIF
            ENDDO
         ENDIF 
         IF (INITIALDIST) THEN
            PD=ABS(ALLPAIRS(NPOSITION))
         ELSE
            DO J5=1,PAIRDISTMAX
               IF (PAIRLIST(J4,J5).EQ.J2) THEN
                  PD=PAIRDIST(J4,J5)
                  NNEIGH=NNEIGH+1
                  IF (WEIGHT(J2).LT.MINNWEIGHT) MINNWEIGHT=WEIGHT(J2)
                  IF (PD.LT.MINNDIST) MINNDIST=PD
                  GOTO 973
               ENDIF
            ENDDO
            DO J5=1,PAIRDISTMAX
               IF (PAIRLIST(J2,J5).EQ.J4) THEN
                  PD=PAIRDIST(J2,J5)
                  NNEIGH=NNEIGH+1
                  IF (WEIGHT(J2).LT.MINNWEIGHT) MINNWEIGHT=WEIGHT(J2)
                  IF (PD.LT.MINNDIST) MINNDIST=PD
                  GOTO 973
               ENDIF
            ENDDO
         ENDIF
973      CONTINUE
         TMPWEIGHT=PD*SCALEFAC 
         IF (TMPWEIGHT.LT.HUGE(1.0D0)/10.0D0) THEN ! don;t raise a huge number to any power!
            IF (INDEXCOSTFUNCTION) THEN 
               IF (TMPWEIGHT.EQ.0.0D0) THEN ! minima are connected!
               ELSE
                  TMPWEIGHT=ABS(J4-J2)
                  IF (DIRECTION.EQ.'BA') THEN
                     IF (J4.LE.NMINA) TMPWEIGHT=NMIN+1-J2 ! not sure that this really makes sense for A and B ! DJW
                     IF (J2.LE.NMINA) TMPWEIGHT=NMIN+1-J4
                  ELSE
                     IF ((J4.LE.NMINA+NMINB).AND.(J4.GT.NMINA)) TMPWEIGHT=NMIN+1-J2
                     IF ((J2.LE.NMINA+NMINB).AND.(J2.GT.NMINA)) TMPWEIGHT=NMIN+1-J4
                  ENDIF
                ENDIF
            ELSEIF (EXPCOSTFUNCTION) THEN 
               IF (TMPWEIGHT.EQ.0.0D0) THEN
                  ! do nothing - don;t set the weight to one !! DJW 22/7/08
               ELSEIF (TMPWEIGHT.GT.700.0D0) THEN
                  TMPWEIGHT=DEXP(700.0D0)
               ELSE
                  TMPWEIGHT=DEXP(TMPWEIGHT)
               ENDIF
            ELSE ! compare squares to favour more small jumps over big ones DJW
               IF (TMPWEIGHT.EQ.0.0D0) THEN
               ELSEIF (COSTFUNCTIONPOWER.EQ.0) THEN
                  TMPWEIGHT=TMPWEIGHT+1.0D0
               ELSEIF (COSTFUNCTIONPOWER.EQ.-1) THEN
                  TMPWEIGHT=1.0D0/TMPWEIGHT
               ELSE
                  TMPWEIGHT=TMPWEIGHT**COSTFUNCTIONPOWER 
               ENDIF
            ENDIF
         ENDIF
         
!        PRINT '(A,2I10,3G20.10)','J2,J4,TMPWEIGHT,WEIGHT(J4),WEIGHT(J2)=',J2,J4,TMPWEIGHT,WEIGHT(J4),WEIGHT(J2)
         IF (TMPWEIGHT+WEIGHT(J4).LT.WEIGHT(J2)) THEN ! relax J2
            RELAXED(J2)=RELAXED(J2)+1
            WEIGHT(J2)=WEIGHT(J4)+TMPWEIGHT
            PARENT(J2)=J4
!           PRINT '(A,2I10)','J2,PARENT=',J2,PARENT(J2)
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
!        PRINT '(A,I10,L5,2G20.10,I10)','J2,PERMANENT,WEIGHT,MINWEIGHT,JMINW=',J2,PERMANENT(J2),WEIGHT(J2),MINWEIGHT,JMINW
      ENDDO
      IF (NOTDONE) THEN
         PRINT '(A,I8,A,I8)','dijinit> WARNING - JMINW not set - value=',JMINW,' J4=',J4
         PRINT '(A,I8)','dijinit> NPERM=',NPERM
         DO J2=1,NMIN
            PRINT '(A,I8,L5,2G20.10)','J2,PERMANENT,WEIGHT,MINWEIGHT=',J2,PERMANENT(J2),WEIGHT(J2),MINWEIGHT
            IF (.NOT.PERMANENT(J2)) THEN
               IF (WEIGHT(J2).LT.MINWEIGHT) THEN
                  MINWEIGHT=WEIGHT(J2)
                  JMINW=J2
                  NOTDONE=.FALSE.
               ENDIF
            ENDIF
         ENDDO
         STOP !!! DJW
      ENDIF

      IF (.NOT.INITIALDIST) PRINT '(A,I10,A,I10,A,2G20.10)','Dijinit> Number of non-permanent nodes in neighbour list for ', &
  &                       J4,' is ',NNEIGH,' min weight and dist=',MINNWEIGHT,MINNDIST
      J4=JMINW
      PERMANENT(J4)=.TRUE.
      NPERM=NPERM+1
!     PRINT '(A,2I8,G20.10,I8)','permanent minimum J4,NPERM,WEIGHT,relaxations=',J4,NPERM,WEIGHT(J4),RELAXED(J4)
      IF (WEIGHT(J4).GT.MAXWEIGHT) THEN
         SCALEFAC=SCALEFAC/10.0D0
         PRINT '(A,G20.10)','dijinit> Maximum weight is too large - scaling by ',SCALEFAC
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
REDODIJKSTRA=.FALSE.
J5=LJ2
NWORST=0
NSTEPS=0
IF (MINGAPT) NMINGAP=0
PRINT '(A)','Dijinit> Summary of best path based on missing connection metric - note distance scaling is removed'
PRINT '(A)','    min1          energy        min2          energy             metric          edge weight            weight'
DO 
   IF (PARENT(J5).EQ.0) THEN
      PRINT '(A,I6,A)','Dijinit> ERROR - parent for J5=',J5,' is zero'
      PRINT '(A)',     'Dijinit> Suggests all possible pairs have been tried!'
      STOP
   ENDIF
   DUMMY=1.0D4*PDMAX*SCALEFAC
!  PRINT '(A,4G20.10)','PDMAX,SCALEFAC,DUMMY,DUMMY/SCALEFAC=',PDMAX,SCALEFAC,DUMMY,DUMMY/SCALEFAC
   IF (.NOT.PAIRSIGNORET) THEN
     DO J2=1,NPAIRDONE ! skip
       IF ((PAIR1(J2).EQ.J5).AND.(PAIR2(J2).EQ.PARENT(J5))) THEN
          IF (INITIALDIST) THEN
             JM=MIN(J5,PARENT(J5))
             JN=MAX(J5,PARENT(J5))
             NPOSITION=((JN-2)*(JN-1))/2+JM
             DUMMY=ABS(ALLPAIRS(NPOSITION))*SCALEFAC
             IF (DUMMY.EQ.0.0D0) THEN
                GOTO 864
             ELSE
                DUMMY=1.0D4*PDMAX*SCALEFAC
                GOTO 864
             ENDIF
          ELSE
             DO J6=1,PAIRDISTMAX
                IF (PAIRLIST(J5,J6).EQ.PARENT(J5)) THEN
                   DUMMY=PAIRDIST(J5,J6)*SCALEFAC
                   IF (DUMMY.EQ.0.0D0) THEN
                      GOTO 864
                   ELSE
                      DUMMY=1.0D4*PDMAX*SCALEFAC
                      GOTO 864
                   ENDIF
                ENDIF
             ENDDO
          ENDIF
          GOTO 864
       ENDIF
       IF ((PAIR1(J2).EQ.PARENT(J5)).AND.(PAIR2(J2).EQ.J5)) THEN
          IF (INITIALDIST) THEN
             JM=MIN(J5,PARENT(J5))
             JN=MAX(J5,PARENT(J5))
             NPOSITION=((JN-2)*(JN-1))/2+JM
             DUMMY=ABS(ALLPAIRS(NPOSITION))*SCALEFAC
             IF (DUMMY.EQ.0.0D0) THEN
                GOTO 864
             ELSE
                DUMMY=1.0D4*PDMAX*SCALEFAC
                GOTO 864
             ENDIF
          ELSE
             DO J6=1,PAIRDISTMAX
                IF (PAIRLIST(PARENT(J5),J6).EQ.J5) THEN
                   DUMMY=PAIRDIST(PARENT(J5),J6)*SCALEFAC
                   IF (DUMMY.EQ.0.0D0) THEN
                       GOTO 864
                   ELSE
                      DUMMY=1.0D4*PDMAX*SCALEFAC
                      GOTO 864
                   ENDIF
                ENDIF
             ENDDO
          ENDIF
          GOTO 864
       ENDIF
     ENDDO
   ENDIF
   IF (PRUNECYCLET.AND.(NPRUNEPAIRSOLD.GT.0)) THEN
      DO J2=1,NPRUNEPAIRSOLD
         IF ((PRUNEPAIRS(1,J2).EQ.J5).AND.(PRUNEPAIRS(2,J2).EQ.PARENT(J5))) THEN
            IF (DEBUG) PRINT '(A,2I8)','pruning> pair used: ',PRUNEPAIRS(1,J2),PRUNEPAIRS(2,J2)
            DUMMY=1.0D4*PDMAX*SCALEFAC
            GOTO 864
         ENDIF
         IF ((PRUNEPAIRS(2,J2).EQ.PARENT(J5)).AND.(PRUNEPAIRS(1,J2).EQ.J5)) THEN
            IF (DEBUG) PRINT '(A,2I8)','pruning> pair used: ',PRUNEPAIRS(1,J2),PRUNEPAIRS(2,J2)
            DUMMY=1.0D4*PDMAX*SCALEFAC
            GOTO 864
         ENDIF
      ENDDO
   ENDIF
   IF (INITIALDIST) THEN
      JM=MIN(J5,PARENT(J5))
      JN=MAX(J5,PARENT(J5))
      NPOSITION=((JN-2)*(JN-1))/2+JM
      DUMMY=ABS(ALLPAIRS(NPOSITION))*SCALEFAC
   ELSE
      DO J2=1,PAIRDISTMAX
         IF (PAIRLIST(J5,J2).EQ.PARENT(J5)) THEN
            DUMMY=PAIRDIST(J5,J2)*SCALEFAC
            GOTO 864
         ENDIF
      ENDDO
      DO J2=1,PAIRDISTMAX
         IF (PAIRLIST(PARENT(J5),J2).EQ.J5) THEN
            DUMMY=PAIRDIST(PARENT(J5),J2)*SCALEFAC
            GOTO 864
         ENDIF
      ENDDO
   ENDIF
864 CONTINUE
   IF (DUMMY.LT.HUGE(1.0D0)/10.0D0) THEN ! don;t raise a huge number to any power!
      IF (INDEXCOSTFUNCTION) THEN
         IF (DUMMY.EQ.0.0D0) THEN ! minima are connected!
            TMPWEIGHT=0.0D0
         ELSE
            TMPWEIGHT=ABS(J5-PARENT(J5))
            IF (DIRECTION.EQ.'AB') THEN
               IF (J5.LE.NMINA) TMPWEIGHT=NMIN+1-PARENT(J5)
               IF (PARENT(J5).LE.NMINA) TMPWEIGHT=NMIN+1-J5
            ELSE
               IF ((PARENT(J5).LE.NMINA+NMINB).AND.(PARENT(J5).GT.NMINA)) TMPWEIGHT=NMIN+1-J5
               IF ((J5.LE.NMINA+NMINB).AND.(J5.GT.NMINA)) TMPWEIGHT=NMIN+1-PARENT(J5)
            ENDIF 
          ENDIF
      ELSEIF (EXPCOSTFUNCTION) THEN ! saves memory and CPU when endpoint separation is very large SAT
         IF (DUMMY.EQ.0.0D0) THEN
            TMPWEIGHT=0.0D0
         ELSEIF (DUMMY.GT.700.0D0) THEN  !numerical limit doubles:10^308 KR
             TMPWEIGHT=DEXP(700.0D0)
         ELSE 
            TMPWEIGHT=DEXP(DUMMY)
         ENDIF
      ELSE ! compare higher powers to favour more small jumps over big ones DJW
         IF (DUMMY.EQ.0.0D0) THEN
            TMPWEIGHT=0.0D0
         ELSEIF (INTERPCOSTFUNCTION) THEN
            TMPWEIGHT=DUMMY**COSTFUNCTIONPOWER
         ELSEIF (COSTFUNCTIONPOWER.EQ.0) THEN
            TMPWEIGHT=1.0D0
         ELSEIF (COSTFUNCTIONPOWER.EQ.-1) THEN
            TMPWEIGHT=1.0D0/TMPWEIGHT
         ELSE
            TMPWEIGHT=DUMMY**COSTFUNCTIONPOWER
         ENDIF
      ENDIF
   ELSE
      TMPWEIGHT=DUMMY
   ENDIF
   NSTEPS=NSTEPS+1
   
   PRINT '(2(I8,G20.10),3G20.10)',J5,EMIN(J5),parent(J5),EMIN(PARENT(J5)),DUMMY/SCALEFAC,TMPWEIGHT,WEIGHT(J5)
   IF (INITIALDIST) THEN
      JM=MIN(J5,PARENT(J5))
      JN=MAX(J5,PARENT(J5))
      NPOSITION=((JN-2)*(JN-1))/2+JM
      IF (ALLPAIRS(NPOSITION).LT.0.0D0) THEN
         READ(UMIN,REC=J5) (LPOINTS1(J2),J2=1,NOPT)
         READ(UMIN,REC=PARENT(J5)) (LPOINTS2(J2),J2=1,NOPT)
         CALL MINPERMDIST(LPOINTS1,LPOINTS2,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE, &
  &                          DIST2,RIGIDBODY,RMAT,.FALSE.)
         PRINT '(A,I10,A,I10,A,G20.10)','Dijinit> true distance for minima ',J5,' and ',PARENT(J5),' is ',DISTANCE
         ALLPAIRS(NPOSITION)=DISTANCE
!
! Need to run Dikstra again and rewrite allpairs
!
         REDODIJKSTRA=.TRUE.
      ENDIF
   ENDIF
   IF (DIJPRUNET) PRUNEMIN(J5)=.TRUE.
   IF (DEBUG.AND.DIJPRUNET) PRINT '(A,I8)','pruning> minimum added to min.retain: ',J5
   THRESH=0.0D0
   IF (BHINTERPT) THRESH=BHDISTTHRESH ! for bhinterp runs raise the threshold to BHDISTTHRESH
   IF (BISECTT) THRESH=BISECTMINDIST ! for bisect runs raise the threshold to BISECTMINDIST
   IF ((DUMMY/SCALEFAC.GT.THRESH).AND.(TMPWEIGHT.LT.HUGE(1.0D0)/10.0D0)) THEN
      NWORST=NWORST+1
      IF (PRUNECYCLET) THEN
         NPRUNEPAIRS=NPRUNEPAIRS+1
         PRUNEPAIRS(1,NPRUNEPAIRS)=J5
         PRUNEPAIRS(2,NPRUNEPAIRS)=PARENT(J5)
      ENDIF
      IF (NWORST.GT.10000) THEN
         PRINT '(A,I8)','ERROR in Dijinit, too many gaps, NWORST=',NWORST
         STOP
      ENDIF
      IF (.NOT.MINGAPT) THEN
         DMIN1(NWORST)=J5
         DMIN2(NWORST)=PARENT(J5)
      ELSE
         IF (DUMMY/SCALEFAC.GT.MINGAPTHRESH) THEN
            NMINGAP=NMINGAP+1
            DMIN1(NMINGAP)=J5
            DMIN2(NMINGAP)=PARENT(J5)
            PRINT '(A,2I8)',' Dijinit> Pair with distance larger than minimum required:',J5,PARENT(J5)
         ELSE
            PRINT '(A,2I8)',' Dijinit> Pair with distance smaller than minimum required.'
         ENDIF
      ENDIF
!     
!  On the first attempt there may only be two minima, and with more than one cpu available
!  we may need to choose the same best path several times before we can get
!  more than one distinct candidate.
!
!      IF (NMIN.GT.2) THEN
!         DO J4=1,PAIRDISTMAX
!            IF (PAIRLIST(J5,J4).EQ.PARENT(J5)) THEN
!               PAIRDIST(J5,J4)=HUGE(1.0D0)
!               GOTO 753
!            ENDIF
!         ENDDO
!753      CONTINUE
!         DO J4=1,PAIRDISTMAX
!            IF (PAIRLIST(PARENT(J5),J4).EQ.J5) THEN
!               PAIRDIST(PARENT(J5),J4)=HUGE(1.0D0)
!               GOTO 751
!            ENDIF
!         ENDDO
!751      CONTINUE
!      ENDIF
   ENDIF
   J5=PARENT(J5)
   IF (J5.EQ.LJ1) EXIT
   IF (J5.EQ.0) EXIT
ENDDO
PRINT '(2(A,I8))','Dijinit> Number of steps=',NSTEPS,' number of missing connections=',NWORST
IF (REDODIJKSTRA) THEN
   LUNIT=GETUNIT()
   OPEN(UNIT=LUNIT,FILE='allpairs',STATUS='UNKNOWN')
   WRITE(LUNIT,'(G20.10)') ALLPAIRS(1:(NMIN*(NMIN-1))/2)
   CLOSE(LUNIT)
   GOTO 642
ENDIF
PRUNEMIN(J5)=.TRUE.
IF (PRUNECYCLET) THEN
   NPRUNEDONE=NPRUNEDONE+1
   PRINT '(A,I8)','Dijinit> Pruning cycle completed: ',NPRUNEDONE
   IF (NPRUNEDONE.LT.NPRUNE) THEN
      DEALLOCATE(LOCATIONSTART,LOCATIONEND)
      GOTO 121
   ENDIF
ENDIF
IF (MINGAPT) THEN
   IF (NMINGAP.GT.0) THEN
      PRINT '(A,I8)','Dijinit> Number of connections larger than minimum distance:',NMINGAP
      NWORST=NMINGAP
   ELSE
      PRINT '(A)','Dijinit> No missing connection is above the required distance cut off.'
      PRINT '(A)','Dijinit> Set MINGAP to false and redo analysis'
      MINGAPT=.FALSE.
      DEALLOCATE(LOCATIONSTART,LOCATIONEND)
      GOTO 121
   ENDIF
ENDIF
IF (DIJPRUNET) THEN !write the best path out to min.retain and then terminate
   OPEN(MUNIT,FILE='min.retain',POSITION='APPEND',ACTION='WRITE',STATUS='NEW')
   NPRUNEMIN=0
   DO J7=1,NMIN
      IF (PRUNEMIN(J7)) NPRUNEMIN=NPRUNEMIN+1
   ENDDO
   WRITE(MUNIT,'(I8)') NPRUNEMIN
   DO J7=1,NMIN
      IF (PRUNEMIN(J7)) THEN
         WRITE(MUNIT,'(I8)') J7
      ENDIF
   ENDDO
   CLOSE(MUNIT)
   PRINT '(A,I6,A)','Dijprune> Best ',NPRUNE,' paths written to min.retain'
   STOP
ENDIF
IF (NWORST.EQ.0) THEN
   PRINT '(A)','Dijinit> Connected path found'
!  IF (DIJCONT) THEN
!     DIJINITT=.FALSE.
!  ELSE
      STOP
!  ENDIF
ENDIF
CALL CPU_TIME(TNEW)
TDIJKSTRA=TDIJKSTRA+TNEW-ELAPSED

RETURN
END
