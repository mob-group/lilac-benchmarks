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

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      SUBROUTINE KMC
      USE PORFUNCS
      USE COMMONS
      IMPLICIT NONE
      INTEGER J1, J2, J3, NCYCLE, NKMC, ATMIN, MIN2, NDEAD, ISTAT, NSTATS, JSTART, JFINISH, LJ1, MODVAL
      LOGICAL CHANGED, DEADTS(NTS), DOLEAP(NMIN)
      DOUBLE PRECISION DMATMC(NCONNMAX,NMIN), RANDOM, WAIT, WAIT2, WAITMEAN, NSTEPMEAN, DUMMY, XNSTEP
      INTEGER MAXVALS, NAVAIL
      PARAMETER (MAXVALS=10000,NSTATS=100)
      INTEGER NCOL(NMIN), NVAL(NCONNMAX,NMIN), NUNCON, ATMINP, DMIN, DMAX
      DOUBLE PRECISION DENOM,NDIST(NMIN), DPRAND, P12P21, PFTOTALBCONN, EMKSUM(NMIN), P2(NMIN), T2(NMIN),
     &                 WAITAB(NMIN), LDUMMY, TDUMMY, PRETMAX, MEANTAU, MEANLENGTH, XSTEPS, RATE, KSUM(NMIN), ELAPSED, TNEW
      LOGICAL ISA(NMIN), ISB(NMIN)

      CALL CPU_TIME(ELAPSED)

      WRITE(*,'(A,G20.10)') 'KMC> Threshold for combining minima is p12*p21>',PAIRTHRESH

      ISA(1:NMIN)=.FALSE.
      ISB(1:NMIN)=.FALSE.
      DO J1=1,NMINA
         ISA(LOCATIONA(J1))=.TRUE.
      ENDDO
      DO J1=1,NMINB
         ISB(LOCATIONB(J1))=.TRUE.
      ENDDO

      CALL RATECONST_SETUP(KSUM,DEADTS,NDEAD,.TRUE.,-300.0D0)

      DO J1=1,NMIN
         EMKSUM(J1)=EXP(-KSUM(J1))
      ENDDO
!
!  Flag transition states that have the energy higher than the total for microcanonical ensemble.
!
      IF (ENSEMBLE.EQ.'E') THEN
         DO J1=1,NTS
            IF (ETS(J1).GT.TOTALE) DEADTS(J1)=.TRUE.
         ENDDO
      ENDIF

!!!!!!!!!!!!!!!!!!!   KMC calculation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  DIRECTION AB: calculate a mean transition time from each b minimum to any a minimum
!  DIRECTION BA: calculate a mean transition time from each a minimum to any b minimum
!
      CALL MAKEDMAT3(DMATMC,NCOL,NVAL,DEADTS,KSUM)
!
!  Check that the stationary point database is actually connected, and remove minima that lie in disjoint graphs.
!  Calculate minimum number of steps of each minimum from the A set.
!
      IF (DIRECTION.EQ.'AB') THEN
         NDIST(1:NMIN)=1000000
         DO J1=1,NMINA
            NDIST(LOCATIONA(J1))=0
         ENDDO
         PRINT *,'locationa(1)=',LOCATIONA(1)
         NCYCLE=0
5        CHANGED=.FALSE.
         NCYCLE=NCYCLE+1
         DMIN=100000
         DMAX=0
         NUNCON=0
         DO J1=1,NMIN
            IF (NDIST(J1).EQ.0) CYCLE ! A minimum
            DO J2=1,NCOL(J1)
               IF (NDIST(NVAL(J2,J1))+1.LT.NDIST(J1)) THEN
                  CHANGED=.TRUE.
                  NDIST(J1)=NDIST(NVAL(J2,J1))+1
               ENDIF
            ENDDO
            IF ((NDIST(J1).GT.DMAX).AND.(NDIST(J1).NE.1000000)) DMAX=NDIST(J1)
            IF (NDIST(J1).LT.DMIN) DMIN=NDIST(J1)
            IF (NDIST(J1).EQ.1000000) NUNCON=NUNCON+1
         ENDDO
         IF (CHANGED) GOTO 5
         PRINT '(3(A,I8))','KMC> steps to A region converged in ',NCYCLE-1,' cycles; maximum=',DMAX,' disconnected=',NUNCON
!
!  Might as well exclude disconnected minima.
!  
         DO J1=1,NMIN
            IF (NDIST(J1).EQ.1000000) NCOL(J1)=0
         ENDDO 
!
!  Calculate minimum number of steps of each minimum from the B set.
!  
      ELSE
         NDIST(1:NMIN)=1000000
         DO J1=1,NMINB
            NDIST(LOCATIONB(J1))=0
         ENDDO
         NCYCLE=0
51       CHANGED=.FALSE.
         NCYCLE=NCYCLE+1
         DMIN=100000
         DMAX=0
         NUNCON=0
         DO J1=1,NMIN
            IF (NDIST(J1).EQ.0) CYCLE ! B minimum
            DO J2=1,NCOL(J1)
               IF (NDIST(NVAL(J2,J1))+1.LT.NDIST(J1)) THEN
                  CHANGED=.TRUE.
                  NDIST(J1)=NDIST(NVAL(J2,J1))+1
               ENDIF
            ENDDO
            IF ((NDIST(J1).GT.DMAX).AND.(NDIST(J1).NE.1000000)) DMAX=NDIST(J1)
            IF (NDIST(J1).LT.DMIN) DMIN=NDIST(J1)
            IF (NDIST(J1).EQ.1000000) NUNCON=NUNCON+1
         ENDDO
         IF (CHANGED) GOTO 51
         PRINT '(3(A,I8))','KMC> steps to B region converged in ',NCYCLE-1,' cycles; maximum=',DMAX,' disconnected=',NUNCON
!
!  Might as well exclude disconnected minima.
!
         DO J1=1,NMIN
            IF (NDIST(J1).EQ.1000000) NCOL(J1)=0
         ENDDO
      ENDIF
C
C  Work out probabilities and times for leapfrog (second neighbour) moves,
C  which should help to stop the KMC getting stuck for pairs of low barrier
C  local minima.
C
      DOLEAP(1:NMIN)=.FALSE.
      P2(1:NMIN)=0.0D0
      T2(1:NMIN)=0.0D0
      IF (PAIRTHRESH.GT.1.0D0) GOTO 555
      DO J1=1,NMIN
         IF (ISA(J1).OR.ISB(J1)) CYCLE
         LDUMMY=1.0D0
         TDUMMY=0.0D0
         PRETMAX=-1.0D0
         DO J2=1,NCOL(J1)
            MIN2=NVAL(J2,J1)
            IF ((.NOT.((DIRECTION.EQ.'AB').AND.ISA(MIN2))).AND.(.NOT.((DIRECTION.EQ.'BA').AND.ISB(MIN2)))) THEN
               DO J3=1,NCOL(MIN2)
                  IF (NVAL(J3,MIN2).EQ.J1) THEN
                     LDUMMY=LDUMMY-DMATMC(J3,MIN2)*DMATMC(J2,J1)
                     TDUMMY=TDUMMY+DMATMC(J2,J1)*EMKSUM(MIN2)
                     IF (DMATMC(J3,MIN2).GT.PRETMAX) PRETMAX=DMATMC(J3,MIN2)
                     GOTO 43
                  ENDIF
               ENDDO
               PRINT*,'KMC> ERROR, J1,J2,NCOL,LDUMMY=',J1,J2,NCOL(ATMIN),LDUMMY
               STOP
43             CONTINUE
            ENDIF
         ENDDO
         IF (LDUMMY.EQ.0.0D0) THEN
            PRINT '(A,I7,A,G20.10)','KMC> ERROR - LDUMMY=0 for minimum ',J1,' PRETMAX=',PRETMAX
C           STOP
         ENDIF
         IF (LDUMMY.GT.0.0D0) THEN
            P2(J1)=1.0D0/LDUMMY
            T2(J1)=(EMKSUM(J1)+TDUMMY)*P2(J1)
            IF (PRETMAX.GT.PAIRTHRESH) DOLEAP(J1)=.TRUE.
            IF (DEBUG) PRINT '(A,I6,A,L5,2(A,G20.10))','KMC> minimum ',J1,' leapfrog moves ',DOLEAP(J1),
     &                 ' probability factor ',P2(J1),' waiting time ',T2(J1)
         ENDIF
C        IF (DOLEAP(J1).AND.DEBUG) THEN
!
!  There are two ways to calculate the leapfrog normalisation, gamma. One is to subtract products
!  of branching ratios from one (above). I suspect that the method below will be more accurate if
!  there are numerical problems. We do this alternative calculation as a check, and use it
!  if there is a discrepancy.
!
         IF (DOLEAP(J1)) THEN
            LDUMMY=0.0D0
            DO J2=1,NCOL(J1)
               MIN2=NVAL(J2,J1)
C              IF (ISA(MIN2).OR.ISB(MIN2)) THEN
               IF (((DIRECTION.EQ.'AB').AND.ISA(MIN2)).OR.((DIRECTION.EQ.'BA').AND.ISB(MIN2))) THEN
                  LDUMMY=LDUMMY+DMATMC(J2,J1)
               ELSE
                  DO J3=1,NCOL(MIN2) ! but not back to J1!
                     IF (NVAL(J3,MIN2).NE.J1) THEN
                        LDUMMY=LDUMMY+DMATMC(J2,J1)*DMATMC(J3,MIN2)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            IF (DEBUG) PRINT '(A,I6,A,G20.10)','KMC>  sum of leapfrog probabilities for minimum ',J1,' is ',LDUMMY
            IF (ABS(LDUMMY*P2(J1)-1.0D0).GT.1.0D-5) THEN
               PRINT '(A,I6,A,G20.10)','KMC> WARNING sum of leapfrog probabilities for minimum ',J1,' is ',LDUMMY*P2(J1)
C              IF (LDUMMY.NE.0.0D0) THEN
               IF (ABS(LDUMMY*P2(J1)-1.0D0).LT.1.0D-1) THEN
                  LDUMMY=1.0D0/LDUMMY
                  PRINT '(A,G20.10)','KMC>  changing waiting time to ',T2(J1)*LDUMMY/P2(J1)
                  PRINT '(A,G20.10)','KMC>  changing probability factor to ',LDUMMY
                  P2(J1)=LDUMMY
               ELSE
                  PRINT '(A,G20.10)','KMC>  turning off leapfrog for minimum ',J1
                  DOLEAP(J1)=.FALSE.
               ENDIF
            ENDIF
         ENDIF
      ENDDO
555   CONTINUE
      IF (DIRECTION.EQ.'AB') THEN
         JSTART=1
         JFINISH=NMINB
      ELSE
         JSTART=1
         JFINISH=NMINA
      ENDIF
      RATE=0.0D0
      ABloop: DO J1=JSTART,JFINISH
         IF (DIRECTION.EQ.'AB') THEN
            LJ1=LOCATIONB(J1)
         ELSE
            LJ1=LOCATIONA(J1)
         ENDIF
         IF (NCOL(LJ1).EQ.0) THEN
            IF (DEBUG) PRINT '(A,I5,A)','KMC> starting minimum ',LJ1,' unconnected - skipping KMC runs'
            WAITAB(LJ1)=0.0D0 ! actually, it should be infinity, but not to worry
            CYCLE ABloop
         ENDIF
         XSTEPS=0.0D0
         MEANTAU=0.0D0
         MEANLENGTH=0.0D0
10       XSTEPS=XSTEPS+1.0D0 ! KMC runs for minmum LJ1
         XNSTEP=0.0D0             ! KMC steps in current KMC run
         ATMIN=LJ1
         WAIT=0.0D0
20       XNSTEP=XNSTEP+1.0D0
C        PRINT '(A,F15.0,I10)','KMC> XNSTEP,ATMIN=',XNSTEP,ATMIN
C        CALL RANDOM_NUMBER(RANDOM)
         RANDOM=DPRAND()
         LDUMMY=0.0D0
C        PRINT '(A,I7,L5)','KMC> ATMIN,DOLEAP=',ATMIN,DOLEAP(ATMIN)
         IF (DOLEAP(ATMIN)) THEN ! move to second neighbour, or neighbouring A or B
C           PRINT '(A,I7,A,I7)','KMC> leapfrog move for minimum ',ATMIN
            DO J2=1,NCOL(ATMIN)
               MIN2=NVAL(J2,ATMIN)
               IF (((DIRECTION.EQ.'AB').AND.ISA(MIN2)).OR.((DIRECTION.EQ.'BA').AND.ISB(MIN2))) THEN
                  LDUMMY=LDUMMY+P2(ATMIN)*DMATMC(J2,ATMIN)
                  IF (LDUMMY.GT.RANDOM-1.0D-7) THEN ! the 0.0000001 is to deal with numerical errors
                     WAIT=WAIT+T2(ATMIN)
                     ATMIN=MIN2
                     GOTO 30 ! we have hit an A or B sink
                  ENDIF
               ELSE
                  DO J3=1,NCOL(MIN2) ! but not back to ATMIN!
                     IF (NVAL(J3,MIN2).NE.ATMIN) THEN
                        LDUMMY=LDUMMY+P2(ATMIN)*DMATMC(J2,ATMIN)*DMATMC(J3,MIN2)
                        IF (LDUMMY.GT.RANDOM-1.0D-7) THEN ! the 0.0000001 is to deal with numerical errors
                           WAIT=WAIT+T2(ATMIN)
                           ATMIN=NVAL(J3,MIN2)
C                          IF (ATMIN.LE.NMINA+NMINB) GOTO 30
                           IF (((DIRECTION.EQ.'AB').AND.ISA(ATMIN)).OR.((DIRECTION.EQ.'BA').AND.ISB(ATMIN))) GOTO 30
                           GOTO 20
                        ENDIF
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
         ELSE
            DO J2=1,NCOL(ATMIN)
               LDUMMY=LDUMMY+DMATMC(J2,ATMIN)
               IF (LDUMMY.GT.RANDOM) THEN
                  WAIT=WAIT+EMKSUM(ATMIN)
                  ATMIN=NVAL(J2,ATMIN)
                  IF (((DIRECTION.EQ.'AB').AND.ISA(ATMIN)).OR.((DIRECTION.EQ.'BA').AND.ISB(ATMIN))) GOTO 30
                  GOTO 20
               ENDIF
            ENDDO
         ENDIF
         PRINT *,'KMC> ERROR, ATMIN,DOLEAP,J2,NCOL,LDUMMY,RANDOM=',ATMIN,DOLEAP(ATMIN),J2,NCOL(ATMIN),LDUMMY,RANDOM
         STOP
30       IF (DEBUG) PRINT '(A,I6,A,F15.1,A,G20.10)','KMC> minimum ',ATMIN,' hit in ',XNSTEP,' KMC steps in time ',WAIT
C        PRINT '(A,I6,A,F15.1,A,G20.10)','KMC> minimum ',ATMIN,' hit in ',XNSTEP,' KMC steps in time ',WAIT
         MEANTAU=MEANTAU+WAIT
         MEANLENGTH=MEANLENGTH+XNSTEP
         IF (DEBUG) THEN
            MODVAL=INT(LOG(XSTEPS*1.0D0)/LOG(10*1.0D0))  ! this should give the integer part of log10(XSTEPS)
            IF (MOD(XSTEPS,10.0D0**MODVAL).EQ.0) PRINT '(A,F12.1,A,G20.10,A,G20.10,A,F12.5)',
     1         'KMC> After ',XSTEPS,' KMC runs <tau>=',MEANTAU/XSTEPS,' <length>=',MEANLENGTH/XSTEPS
         ENDIF
         IF (XSTEPS.LT.NKMCCYCLES) GOTO 10
         WAITAB(LJ1)=MEANTAU/XSTEPS
         IF (DIRECTION.EQ.'AB') THEN
            PRINT '(A,I6,A,G20.10,A,F12.1,A,G20.10,A,G10.3)',
     &            'KMC>  T_b min ',LJ1,' is ',WAITAB(LJ1),
     &            ' <KMC steps>=',MEANLENGTH/XSTEPS,' tau_b ',EMKSUM(LJ1), 
     &            ' ratio ',WAITAB(LJ1)/EMKSUM(LJ1)
            RATE=RATE+EXP(PFMIN(LJ1)-PFTOTALB)*XSTEPS/MEANTAU
         ELSE
            PRINT '(A,I6,A,G20.10,A,F12.1,A,G20.10,A,G10.3)',
     &            'KMC>  T_a min ',LJ1,' is ',WAITAB(LJ1),
     &            ' <KMC steps>=',MEANLENGTH/XSTEPS,' tau_a ',EMKSUM(LJ1), 
     &            ' ratio ',WAITAB(LJ1)/EMKSUM(LJ1)
            RATE=RATE+EXP(PFMIN(LJ1)-PFTOTALA)*XSTEPS/MEANTAU
         ENDIF

         CALL FLUSH(6)
      ENDDO ABloop
      IF (DIRECTION.EQ.'BA') THEN
         PRINT '(2(A,G20.10))','KMC>  from detailed balance: k(A<-B)=',RATE*EXP(PFTOTALA-PFTOTALB),' from KMC: k(B<-A)=',RATE
      ELSE
         PRINT '(2(A,G20.10))','KMC>  from KMC: k(A<-B)=',RATE,' from detailed balance: k(B<-A)=',RATE*EXP(PFTOTALB-PFTOTALA)
      ENDIF
!
!!!!!!!!!!!!!!!!!!!   end of KMC calculation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      CALL CPU_TIME(TNEW)
      WRITE(*,'(A,G15.5,A)') 'KMC> CPU time spent in KMC                          =',TNEW-ELAPSED,' s'

      STOP
      END
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C  Construct DMAT0.
C
C  Degenerate rearrangements are excluded.
C
      SUBROUTINE MAKEDMAT3(DMAT0,NCOL,NVAL,DEADTS,KSUM)
      USE COMMONS
      IMPLICIT NONE
      INTEGER NCOL(NMIN), J1, M1, M2, J2, NVAL(NCONNMAX,NMIN)
      DOUBLE PRECISION DMAT0(NCONNMAX,NMIN), TEMP0, DUMMY, KSUM(NMIN)
      LOGICAL NONZERO, DEADTS(NTS)

      NCOL(1:NMIN)=0
C
C  We are setting up DMAT(thing,M2)
C  Cycle over connected transition states using pointers for minimum M2.
C
      DO M2=1,NMIN   
         DO M1=1,NMIN  
            TEMP0=0.0D0
            NONZERO=.FALSE.
            J1=TOPPOINTER(M2)  !  sets J1 to the ts connected to minimum M2 with the highest value
            IF (J1.GT.0) THEN
11             IF ((.NOT.DEADTS(J1)).AND.(M2.NE.M1)) THEN
                  IF ((MINUS(J1).EQ.M1).AND.(PLUS(J1).EQ.M2)) THEN  !  M2 M1  
                     NONZERO=.TRUE.
                     TEMP0=TEMP0+EXP(KPLUS(J1)-KSUM(PLUS(J1)))
                  ENDIF
                  IF ((PLUS(J1).EQ.M1).AND.(MINUS(J1).EQ.M2)) THEN  !  M1 M2 
                     NONZERO=.TRUE.
                     TEMP0=TEMP0+EXP(KMINUS(J1)-KSUM(MINUS(J1)))
                  ENDIF
               ENDIF
               IF (PLUS(J1).EQ.M2) THEN
                  J1=POINTERP(J1)
               ELSE IF (MINUS(J1).EQ.M2) THEN
                  J1=POINTERM(J1)
               ENDIF
               IF (J1.GT.0) GOTO 11
            ENDIF
            IF (NONZERO) THEN
               NCOL(M2)=NCOL(M2)+1
               NVAL(NCOL(M2),M2)=M1
               DMAT0(NCOL(M2),M2)=TEMP0
            ENDIF
         ENDDO
      ENDDO
C
C  Check row normalisation. 
C
      DO J1=1,NMIN
         DUMMY=0.0D0
         DO J2=1,NCOL(J1)
            DUMMY=DUMMY+DMAT0(J2,J1)
            IF (DEBUG) WRITE(*,'(A,3I6,3G20.10)') 'makedmat3> J1,J2,NVAL,DMAT0,sum=',J1,J2,NVAL(J2,J1),DMAT0(J2,J1),DUMMY
         ENDDO
         IF (DEBUG) WRITE(*,'(A,2I6,3G20.10)') 'J1,ncol,sum=',J1,NCOL(J1),DUMMY
         IF ((NCOL(J1).GT.0).AND.(ABS(DUMMY-1.0D0).GT.1.0D-5)) THEN
            PRINT*,'ERROR - J1,NCOL(J1),DUMMY=',J1,NCOL(J1),DUMMY
            STOP
         ENDIF
      ENDDO

      RETURN
      END
