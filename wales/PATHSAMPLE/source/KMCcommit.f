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
!  Calculate rate constant by committor probability analysis combined with a waiting time.
!
      SUBROUTINE KMCCOMMIT
      USE COMMONS
      USE PORFUNCS
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4, MODVAL, JSTART, JFINISH, UNFROZEN, LJ1, LJ3, NCYCLE
      INTEGER MINMAP(NMIN), MIN2, NONZERO, DMIN, DMAX, NUNCON
      LOGICAL DEADTS(NTS), LFROZEN(NMIN), DOLEAP(NMIN), CHANGED
      INTEGER NCOL(NMIN), NVAL(NCONNMAX,NMIN), NSTEP, ATMIN, NDEAD, ISTAT, NCOUNT, NDIST(NMIN), NCOLSAVE(NMIN)
      DOUBLE PRECISION MEANTAU, MEANLENGTH, XSTEPS, RANDOM, PDUMMY, TOL, 
     &                 WAITAB(NMIN), EMKSUM(NMIN), LKAB, LKBA, LKTOTAL, PLKTOTAL, 
     &                 WAIT, TAUBAR, PFOLD(NMIN), NEWPFOLD(NMIN), TDUM(NMIN), 
     &                 P2(NMIN), T2(NMIN), PRETMAX, TDUMMY, TRUEK, SSK, DPRAND, KSUM(NMIN)
      DOUBLE PRECISION, ALLOCATABLE :: DVEC(:)
      INTEGER UNFROZENINDEX(NMIN), NCOLPREV, PUNFROZEN, ROW_PTR(NMIN+1)
      INTEGER, ALLOCATABLE :: COL_IND(:)
      DOUBLE PRECISION DMATMC(NCONNMAX,NMIN), LDUMMY
C     INTEGER, PARAMETER :: extra = selected_real_kind(20,200)
C     REAL(KIND=extra) :: LDUMMY, DMATMC(NCONNMAX,NMIN)
      INTEGER NMINASAVE, NMINBSAVE, LOCATIONASAVE(NMINA), LOCATIONBSAVE(NMINB)
      LOGICAL ISA(NMIN), ISB(NMIN)

      ISA(1:NMIN)=.FALSE.
      ISB(1:NMIN)=.FALSE.
      IF (REGROUPT) THEN
         NMINASAVE=NMINA
         NMINBSAVE=NMINB
         LOCATIONASAVE(1:NMINA)=LOCATIONA(1:NMINA)
         LOCATIONBSAVE(1:NMINB)=LOCATIONB(1:NMINB)
         CALL REGROUP(MINMAP) ! here NMINA and NMINB can change! Minima are also reordered!
         DEALLOCATE(LOCATIONA,LOCATIONB)
         ALLOCATE(LOCATIONA(NMINA),LOCATIONB(NMINB))
         DO J1=1,NMINA
            LOCATIONA(J1)=J1
         ENDDO
         DO J1=1,NMINB
            LOCATIONB(J1)=NMINA+J1
         ENDDO
      ELSE
         DO J1=1,NMIN
            MINMAP(J1)=J1
         ENDDO
      ENDIF
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
C
C  Flag transition states that have the energy higher than the total for microcanonical ensemble.
C  
      IF (ENSEMBLE.EQ.'E') THEN
         DO J1=1,NTS
            IF (ETS(J1).GT.TOTALE) DEADTS(J1)=.TRUE.
         ENDDO 
      ENDIF

!!!!!!!!!!!!!!!!!!!   KMC calculation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  DIRECTION AB: calculate a mean transition time from each b minimum to any other b or a minimum
!  DIRECTION BA: calculate a mean transition time from each a minimum to any other b or a minimum
!
!  If NKMCCYCLES = 0 use the waiting time for any escape from the a or b minimum as a
!  lower bound on the above mean transition time to A or B. If the steady state approximation
!  really holds for intervening minima then there should be no difference.
!
      IF (NKMCCYCLES.EQ.0) THEN
         WAITAB(1:NMIN)=EMKSUM(1:NMIN)
         GOTO 800
      ENDIF
!
!  We need to allow transitions to all A and B minima - DMATMC is different from P^fold part
!
      CALL MAKED(DMATMC,NCOL,NVAL,DEADTS,.TRUE.,ISA,ISB,KSUM) 
!
!  Check that the stationary point database is actually connected, and remove minima that lie in disjoint graphs.
!  Calculate minimum number of steps of each minimum from the A set.
!
      IF (DIRECTION.EQ.'AB') THEN
         NDIST(1:NMIN)=1000000
         DO J1=1,NMINA
            NDIST(LOCATIONA(J1))=0
         ENDDO
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
         PRINT '(3(A,I8))','KMCcommit> steps to A region converged in ',NCYCLE-1,' cycles; maximum=',DMAX,' disconnected=',NUNCON
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
            IF (J1.EQ.31) PRINT '(A,2I8,L5)','J1,NDIST,CHANGED=',J1,NDIST(J1),CHANGED
            IF ((NDIST(J1).GT.DMAX).AND.(NDIST(J1).NE.1000000)) DMAX=NDIST(J1)
            IF (NDIST(J1).LT.DMIN) DMIN=NDIST(J1)
            IF (NDIST(J1).EQ.1000000) NUNCON=NUNCON+1
         ENDDO 
         IF (CHANGED) GOTO 51
         PRINT '(3(A,I8))','KMCcommit> steps to B region converged in ',NCYCLE-1,' cycles; maximum=',DMAX,' disconnected=',NUNCON
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
C     DO J1=NMINA+NMINB+1,NMIN
      DO J1=1,NMIN
         IF (ISA(J1).OR.ISB(J1)) CYCLE
         LDUMMY=1.0D0
         TDUMMY=0.0D0
         PRETMAX=-1.0D0
         DO J2=1,NCOL(J1)
            MIN2=NVAL(J2,J1)
            IF ((.NOT.ISA(MIN2)).AND.(.NOT.ISB(MIN2))) THEN
               DO J3=1,NCOL(MIN2)
                  IF (NVAL(J3,MIN2).EQ.J1) THEN
                     LDUMMY=LDUMMY-DMATMC(J3,MIN2)*DMATMC(J2,J1)
                     TDUMMY=TDUMMY+DMATMC(J2,J1)*EMKSUM(MIN2)
                     IF (DMATMC(J3,MIN2).GT.PRETMAX) PRETMAX=DMATMC(J3,MIN2)
                     GOTO 43
                  ENDIF
               ENDDO 
               PRINT*,'KMCcommit> ERROR, J1,J2,NCOL,LDUMMY=',J1,J2,NCOL(ATMIN),LDUMMY
               STOP
43             CONTINUE
            ENDIF
         ENDDO
         IF (LDUMMY.EQ.0.0D0) THEN
            PRINT '(A,I7,A,G20.10)','KMCcommit> ERROR - LDUMMY=0 for minimum ',J1,' PRETMAX=',PRETMAX
C           STOP
         ENDIF
         IF (LDUMMY.GT.0.0D0) THEN
            P2(J1)=1.0D0/LDUMMY
            T2(J1)=(EMKSUM(J1)+TDUMMY)*P2(J1)
            IF (PRETMAX.GT.PAIRTHRESH) DOLEAP(J1)=.TRUE.
            IF (DEBUG) PRINT '(A,I6,A,L5,2(A,G20.10))','KMCcommit> minimum ',J1,' leapfrog moves ',DOLEAP(J1),
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
C              IF (MIN2.LE.NMINA+NMINB) THEN
               IF (ISA(MIN2).OR.ISB(MIN2)) THEN
                  LDUMMY=LDUMMY+DMATMC(J2,J1)
               ELSE
                  DO J3=1,NCOL(MIN2) ! but not back to J1!
                     IF (NVAL(J3,MIN2).NE.J1) THEN
                        LDUMMY=LDUMMY+DMATMC(J2,J1)*DMATMC(J3,MIN2)
                     ENDIF
                  ENDDO
               ENDIF
            ENDDO
            IF (DEBUG) PRINT '(A,I6,A,G20.10)','KMCcommit>  sum of leapfrog probabilities for minimum ',J1,' is ',LDUMMY
            IF (ABS(LDUMMY*P2(J1)-1.0D0).GT.1.0D-5) THEN
               PRINT '(A,I6,A,G20.10)','KMCcommit> WARNING sum of leapfrog probabilities for minimum ',J1,' is ',LDUMMY*P2(J1)
               IF (LDUMMY.NE.0.0D0) THEN
                  LDUMMY=1.0D0/LDUMMY
                  PRINT '(A,G20.10)','KMCcommit>  changing waiting time to ',T2(J1)*LDUMMY/P2(J1)
                  PRINT '(A,G20.10)','KMCcommit>  changing probability factor to ',LDUMMY
                  P2(J1)=LDUMMY
               ELSE
                  PRINT '(A,G20.10)','KMCcommit>  turning off leapfrog for minimum ',J1
                  DOLEAP(J1)=.FALSE.
               ENDIF
            ENDIF
         ENDIF
      ENDDO
555   CONTINUE

      IF (DIRECTION.EQ.'AB') THEN
C        JSTART=NMINA+1
C        JFINISH=NMINA+NMINB
         JSTART=1
         JFINISH=NMINB
      ELSE
C        JSTART=1
C        JFINISH=NMINA
         JSTART=1
         JFINISH=NMINA
      ENDIF
      ABloop: DO J1=JSTART,JFINISH
         IF (DIRECTION.EQ.'AB') THEN
            LJ1=LOCATIONB(J1)
         ELSE
            LJ1=LOCATIONA(J1)
         ENDIF
         IF (NCOL(LJ1).EQ.0) THEN
            IF (DEBUG) PRINT '(A,I5,A)','KMCcommit> minimum ',LJ1,' unconnected - skipping KMC runs'
            WAITAB(LJ1)=0.0D0 ! actually, it should be infinity, but not to worry
            CYCLE ABloop
         ENDIF
         XSTEPS=0.0D0
         MEANTAU=0.0D0
         MEANLENGTH=0.0D0
10       XSTEPS=XSTEPS+1.0D0 ! KMC runs for minmum LJ1
         NSTEP=0             ! KMC steps in current KMC run
         ATMIN=LJ1
         WAIT=0.0D0
20       NSTEP=NSTEP+1
C        PRINT '(A,3I10)','KMCcommit> NSTEP,ATMIN,formerly=',NSTEP,ATMIN,MINMAP(ATMIN) 
C        CALL RANDOM_NUMBER(RANDOM)
         RANDOM=DPRAND()
         LDUMMY=0.0D0
C        PRINT '(A,I7,L5)','KMCcommit> ATMIN,DOLEAP=',ATMIN,DOLEAP(ATMIN)
         IF (DOLEAP(ATMIN)) THEN ! move to second neighbour, or neighbouring A or B
C           PRINT '(A,I7,A,I7)','KMCcommit> leapfrog move for minimum ',ATMIN,' formerly ',MINMAP(ATMIN) 
            DO J2=1,NCOL(ATMIN)
               MIN2=NVAL(J2,ATMIN)
C              IF (MIN2.LE.NMINA+NMINB) THEN
               IF (ISA(MIN2).OR.ISB(MIN2)) THEN
                  LDUMMY=LDUMMY+P2(ATMIN)*DMATMC(J2,ATMIN)
                  IF (LDUMMY.GT.RANDOM-1.0D-7) THEN ! the 0.0000001 is to deal with numerical errors
                     WAIT=WAIT+T2(ATMIN)
                     ATMIN=MIN2
                     GOTO 30 ! we have hit a B or A minimum
                  ENDIF
               ELSE
                  DO J3=1,NCOL(MIN2) ! but not back to ATMIN!
                     IF (NVAL(J3,MIN2).NE.ATMIN) THEN
                        LDUMMY=LDUMMY+P2(ATMIN)*DMATMC(J2,ATMIN)*DMATMC(J3,MIN2)
                        IF (LDUMMY.GT.RANDOM-1.0D-7) THEN ! the 0.0000001 is to deal with numerical errors
                           WAIT=WAIT+T2(ATMIN)
                           ATMIN=NVAL(J3,MIN2)
C                          IF (ATMIN.LE.NMINA+NMINB) GOTO 30
                           IF (ISA(ATMIN).OR.ISB(ATMIN)) GOTO 30
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
C                 IF (ATMIN.LE.NMINA+NMINB) GOTO 30
                  IF (ISA(ATMIN).OR.ISB(ATMIN)) GOTO 30
                  GOTO 20
               ENDIF
            ENDDO
         ENDIF
         PRINT*,'KMCcommit> ERROR, ATMIN,J2,NCOL,LDUMMY,RANDOM=',ATMIN,J2,NCOL(ATMIN),LDUMMY,RANDOM
         STOP
30       IF (DEBUG) PRINT '(A,I6,A,I6,A,I10,A,G20.10)','KMCcommit> minimum ',ATMIN,' formerly ',MINMAP(ATMIN),
     1                    ' hit in ',NSTEP,' KMC steps in time ',WAIT
C        PRINT '(A,I6,A,I6,A,I10,A,G20.10)','KMCcommit> minimum ',ATMIN,' formerly ',MINMAP(ATMIN),
C    1                    ' hit in ',NSTEP,' KMC steps in time ',WAIT
         MEANTAU=MEANTAU+WAIT
         MEANLENGTH=MEANLENGTH+NSTEP
         IF (DEBUG) THEN
            MODVAL=INT(LOG(XSTEPS*1.0D0)/LOG(10*1.0D0))  ! this should give the integer part of log10(XSTEPS)
            IF (MOD(XSTEPS,10.0D0**MODVAL).EQ.0) PRINT '(A,F12.1,A,G20.10,A,G20.10,A,F12.5)',
     1         'KMCcommit> After ',XSTEPS,' KMC runs <tau>=',MEANTAU/XSTEPS,' <length>=',MEANLENGTH/XSTEPS
         ENDIF
         IF (XSTEPS.LT.NKMCCYCLES) GOTO 10
         WAITAB(LJ1)=MEANTAU/XSTEPS
         IF (DIRECTION.EQ.'AB') THEN
            PRINT '(A,I6,A,I6,A,G20.10,A,F12.1,A,G20.10,A,G10.3)',
     &            'KMCcommit>  t_b min ',LJ1,' formerly ',MINMAP(LJ1),' is ',WAITAB(LJ1),
     &            ' <KMC steps>=',MEANLENGTH/XSTEPS,' tau_b ',EMKSUM(LJ1),
     &            ' ratio ',WAITAB(LJ1)/EMKSUM(LJ1)
         ELSE
            PRINT '(A,I6,A,I6,A,G20.10,A,F12.1,A,G20.10,A,G10.3)',
     &            'KMCcommit>  t_a min ',LJ1,' formerly ',MINMAP(LJ1),' is ',WAITAB(LJ1),
     &            ' <KMC steps>=',MEANLENGTH/XSTEPS,' tau_a ',EMKSUM(LJ1),
     &            ' ratio ',WAITAB(LJ1)/EMKSUM(LJ1)
         ENDIF
         CALL FLUSH(6)
      ENDDO ABloop
!
!!!!!!!!!!!!!!!!!!!   end of KMC calculation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
800   CONTINUE

!!!!!!!!!!!!!!!!!!!   P^fold calculation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  IF DIRECTION is AB then we want P^fold->A and A minima are sinks
!  IF DIRECTION is BA then we want P^fold->B and B minima are sinks
!
      NCOLSAVE(1:NMIN)=NCOL(1:NMIN)
      CALL MAKED(DMATMC,NCOL,NVAL,DEADTS,.FALSE.,ISA,ISB,KSUM) 
C
C  Now iterate PFOLD's to convergence. 
C  For DIRECTION AB calculate PFA directly: initial values for A minima are 1, the rest are 0.
C  For DIRECTION BA calculate PFB directly: initial values for B minima are 1, the rest are 0.
C  Gauss-Seidel iteration if OMEGA=1: successive over-relaxation if 1<OMEGA<2.
C  We get the same answer whether we calculate PFA or PFB, but it may converge
C  faster one way.
C
      UNFROZEN=0
      DO J3=1,NMIN
         IF (NCOL(J3).EQ.0) THEN
            LFROZEN(J3)=.TRUE.
         ELSE
            LFROZEN(J3)=.FALSE.
            UNFROZEN=UNFROZEN+1
            UNFROZENINDEX(UNFROZEN)=J3
         ENDIF
      ENDDO
C
C  Pfold values were initialised in setup.f, and may have been read from commit.data.
C  We only want to reset the values that may not have been unity before regrouping.
C
      IF (DIRECTION.EQ.'AB') THEN ! PFA
C        PFOLD(1:NMIN)=GPFOLD(1:NMIN)
         PFOLD(1:NMIN)=0.0D0
         DO J1=1,NMINA
            PFOLD(LOCATIONA(J1))=1.0D0
         ENDDO
      ELSE ! PFB
C        PFOLD(1:NMIN)=GPFOLD(1:NMIN)
         PFOLD(1:NMIN)=0.0D0
         DO J1=1,NMINB
            PFOLD(LOCATIONB(J1))=1.0D0
         ENDDO
      ENDIF
      NEWPFOLD(1:NMIN)=PFOLD(1:NMIN)
      PDUMMY=1.0D100
      PLKTOTAL=1.0D100
C
C  Make compressed row storage for DMAT.
C  We are only going to save entries for unfrozen minima - data structure
C  must be updated when new minima join the frozen list.
C
      NONZERO=0
      DO J1=1,NMIN
         NONZERO=NONZERO+NCOL(J1)
      ENDDO
      ALLOCATE(DVEC(NONZERO),COL_IND(NONZERO))
      NCOUNT=0
      ROW_PTR(1)=1
      UNFROZEN=0
      DO J1=1,NMIN
         IF (.NOT.LFROZEN(J1)) THEN
            UNFROZEN=UNFROZEN+1
            IF (UNFROZEN.GT.1) ROW_PTR(UNFROZEN)=ROW_PTR(UNFROZEN-1)+NCOLPREV
            NCOLPREV=NCOL(J1)
            DO J2=1,NCOL(J1)
               NCOUNT=NCOUNT+1
               DVEC(NCOUNT)=DMATMC(J2,J1)
               COL_IND(NCOUNT)=NVAL(J2,J1)
            ENDDO
         ENDIF
      ENDDO
      ROW_PTR(UNFROZEN+1)=NONZERO+1

      J1=0
C
C  Main P^fold loop.
C
      itloop: DO
         J1=J1+1
C
C        OMEGA is the damping factor for successive overrelaxation method (SOR)
C        OMEGA=1 is pure Gauss-Seidel. OMEGA should be < 2
C
         DO J3=1,UNFROZEN
            LDUMMY=0.0D0
            DO J2=ROW_PTR(J3),ROW_PTR(J3+1)-1
C              LDUMMY=LDUMMY+DVEC(J2)*PFOLD(COL_IND(J2))    ! Jacobi
               LDUMMY=LDUMMY+DVEC(J2)*NEWPFOLD(COL_IND(J2)) ! Gauss-Seidel, a bit faster
            ENDDO 
            J4=UNFROZENINDEX(J3)
            NEWPFOLD(J4)=LDUMMY
            PFOLD(J4)=OMEGA*LDUMMY+(1.0D0-OMEGA)*PFOLD(J4)  ! SOR
C           PRINT '(A,4I6,2G20.10)','KMCcommit> J3,J4,limits,LDUMMY,PFOLD=',J3,J4,ROW_PTR(J3),ROW_PTR(J3+1)-1,LDUMMY,PFOLD(J4)
C           NEWPFOLD(J4)=MAX(MIN(LDUMMY,1.0D0),0.0D0)
C           PFOLD(J4)=MAX(MIN(OMEGA*LDUMMY+(1.0D0-OMEGA)*PFOLD(J4),1.0D0),0.0D0)  ! SOR
         ENDDO
C
C Print results for the following numbers of iterations:
C 2 3 4 5 6 7 8 9 10 20 30 40 50 60 70 80 90 100 200 300 400 500 600 700 800 900 1000 2000 3000 etc.
C
         MODVAL=INT(LOG(J1*1.0D0)/LOG(10*1.0D0))  ! this should give the integer part of log10(j1)
         IF ((MOD(J1,10**MODVAL).EQ.0).OR.(J1*1.0D0.GE.MAXBREAK)) THEN
            TOL=PABCONV
            TAUBAR=0.0D0
            IF (DEBUG) PRINT'(A)','   min       P_B                    P_A                    %'
            IF (DIRECTION.EQ.'BA') THEN  ! PFB calculated directly
               DO J3=1,NMINA
                  LJ3=LOCATIONA(J3)
C
C  We effectively assign zero rates to unconnected minima. Different from KMC routine.
C
                  IF (WAITAB(LJ3).GT.0) TAUBAR=TAUBAR+(PFOLD(LJ3)/WAITAB(LJ3))*EXP(PFMIN(LJ3)-PFTOTALA)
               ENDDO
               LKAB=EXP(PFTOTALA-PFTOTALB)*TAUBAR
               LKBA=TAUBAR
               IF (DEBUG) THEN
                  DO J3=1,NMIN
                  IF ((.NOT.LFROZEN(J3)).AND.NCOL(J3).GT.0) WRITE(*,'(I6,3G20.10)') 
     &                J3,PFOLD(J3),MAX(1.0D0-PFOLD(J3),0.0D0),ABS(PFOLD(J3)-NEWPFOLD(J3))*100.0/MAX(PFOLD(J3),1.0D-200)
                  ENDDO
               ENDIF
            ELSE  ! PFA calculated directly
C              DO J3=NMINA+1,NMINA+NMINB
               DO J3=1,NMINB
                  LJ3=LOCATIONB(J3)
                  IF (WAITAB(LJ3).GT.0) THEN
C
C  We effectively assign zero rates to unconnected minima. Different from KMC routine.
C
                     TAUBAR=TAUBAR+(PFOLD(LJ3)/WAITAB(LJ3))*EXP(PFMIN(LJ3)-PFTOTALB)
                  ENDIF
               ENDDO
               LKAB=TAUBAR
               LKBA=EXP(PFTOTALB-PFTOTALA)*TAUBAR
               IF (DEBUG) THEN
                  DO J3=1,NMIN
                     IF ((.NOT.LFROZEN(J3)).AND.NCOL(J3).GT.0) WRITE(*,'(I6,3G20.10)') 
     &                J3,MAX(1.0D0-PFOLD(J3),0.0D0),PFOLD(J3),ABS(PFOLD(J3)-NEWPFOLD(J3))*100.0/MAX(PFOLD(J3),1.0D-200)
                  ENDDO
               ENDIF
            ENDIF
            IF (J1.GT.10000) THEN
C           IF (.FALSE.) THEN
               PUNFROZEN=UNFROZEN
               UNFROZEN=0
               DO J3=1,NMIN
                  IF (.NOT.LFROZEN(J3)) THEN
                     IF ((ABS(PFOLD(J3)-NEWPFOLD(J3))/MAX(PFOLD(J3),1.0D-200).LT.TOL).AND.(PFOLD(J3).GT.0.0D0)) THEN
                        LFROZEN(J3)=.TRUE.
                        PFOLD(J3)=MAX(MIN(PFOLD(J3),1.0D0),0.0D0)
                        NEWPFOLD(J3)=MAX(MIN(NEWPFOLD(J3),1.0D0),0.0D0)
                        IF (DEBUG) PRINT '(A,I6,A,G20.10)','KMCcommit> freezing P^fold for minimum ',J3,' value=',PFOLD(J3)
C
C  Must update PUNFROZEN, UNFROZEN, NONZERO, DVEC, ROW_PTR, COL_IND to remove
C  the entries for this minimum that is now being frozen.
C
                        DO J4=ROW_PTR(UNFROZEN+1),NONZERO-NCOL(J3)
                           DVEC(J4)=DVEC(J4+NCOL(J3))
                           COL_IND(J4)=COL_IND(J4+NCOL(J3))
                        ENDDO
                        DO J4=UNFROZEN+2,PUNFROZEN-1
                           ROW_PTR(J4)=ROW_PTR(J4-1)+(ROW_PTR(J4+1)-ROW_PTR(J4))
                        ENDDO
                        DO J4=UNFROZEN+1,PUNFROZEN-1
                           UNFROZENINDEX(J4)=UNFROZENINDEX(J4+1)
                        ENDDO
                        NONZERO=NONZERO-NCOL(J3)
                        PUNFROZEN=PUNFROZEN-1
C                    ELSEIF (ABS(PFOLD(J3)-1.0D0).LT.1.0D-2) THEN
C                       LFROZEN(J3)=.TRUE.
C                       PRINT '(A,I6,A,G20.10)','KMCcommit> freezing P^fold for minimum ',J3,' value=',PFOLD(J3)
C                       IF (DEBUG) PRINT '(A,I6,A,G20.10)','KMCcommit> freezing P^fold for minimum ',J3,' value=',PFOLD(J3)
                     ELSE
                        UNFROZEN=UNFROZEN+1
                     ENDIF
                  ENDIF
               ENDDO
               IF (DEBUG.OR.(UNFROZEN.EQ.0)) PRINT '(A,I6)','KMCcommit> number of unfrozen minima=',UNFROZEN
               PRINT '(A,I6)','KMCcommit> number of unfrozen minima=',UNFROZEN
            ENDIF
            LKTOTAL=LKAB+LKBA
            LDUMMY=0.0D0
            IF ((PLKTOTAL.GT.0.0D0).AND.(LKTOTAL.GT.0.0D0)) LDUMMY=ABS(LKTOTAL-PLKTOTAL)/LKTOTAL 
            WRITE(*,'(A,I10,A,G13.3,A,G20.10,A,G20.10,A,G20.10)') 'P_fold iter ',J1,' conv=',LDUMMY,
     1               ' a<-b ',LKAB,' b<-a ',LKBA,' total ',LKTOTAL
            CALL FLUSH(6)
            IF ((PLKTOTAL.GT.0.0D0).AND.(LKTOTAL.GT.0.0D0)) THEN
               IF (LDUMMY.LT.TOL) EXIT itloop
            ENDIF
            PLKTOTAL=LKTOTAL
            IF (J1*1.0D0.GE.MAXBREAK) THEN
               WRITE(*,'(A,I10,A,G15.5)') 'WARNING - iterative P_(fold) calculation failed to converge in ',
     &                                    J1,' steps, RMS=',LDUMMY
               EXIT itloop
            ENDIF
         ENDIF
         PDUMMY=LDUMMY
         CYCLE itloop
      ENDDO itloop
!
!  check solution
!
      LDUMMY=0.0D0
      NCOUNT=0
      DO J1=1,NMIN
         TDUM(J1)=0.0D0
         DO J2=1,NCOL(J1)
            TDUM(J1)=TDUM(J1)+DMATMC(J2,J1)*NEWPFOLD(NVAL(J2,J1))
         ENDDO
         IF (NCOL(J1).GT.0) THEN
            NCOUNT=NCOUNT+1
            LDUMMY=LDUMMY+(TDUM(J1)-NEWPFOLD(J1))**2
         ENDIF
      ENDDO
      LDUMMY=SQRT(LDUMMY/NCOUNT)
      WRITE(*,'(A,G20.10)') 'RMS residual error=',LDUMMY
      CALL FLUSH(6)
      PRINT*
      IF (DEBUG) THEN
         PRINT'(A)','   min       P_B                    P_A'
         J2=0
         J3=0
         J4=0
         IF (DIRECTION.EQ.'BA') THEN
            DO J1=1,NMIN 
               IF (NCOL(J1).GT.0) THEN
                  WRITE(*,'(I6,3G20.10)') J1,PFOLD(J1),MAX(1.0D0-PFOLD(J1),0.0D0)
                  J2=J2+1
                  IF (DABS(PFOLD(J1)-1.0D0).LT.1.0D-2) J3=J3+1
                  IF (DABS(PFOLD(J1)).LT.1.0D-2) J4=J4+1
               ENDIF
            ENDDO
         ELSE
            DO J1=1,NMIN
               IF (NCOL(J1).GT.0) THEN
                  WRITE(*,'(I6,3G20.10)') J1,MAX(1.0D0-PFOLD(J1),0.0D0),PFOLD(J1)
                  J2=J2+1
                  IF (DABS(PFOLD(J1)-1.0D0).LT.1.0D-2) J3=J3+1
                  IF (DABS(PFOLD(J1)).LT.1.0D-2) J4=J4+1
               ENDIF
            ENDDO
         ENDIF
         PRINT '(A,I6,A,I6,A,I6,A)','KMCcommit> For ',J2,' active minima ',J3,' are close to 1 and ',J4,' are close to 0'
      ENDIF
!
! compare actual result versus steady-state approximation for I minima
!
      PRINT '(A)','KMCcommit> (minimum labels refer to before regrouping)'
      IF (NKMCCYCLES.EQ.0) PRINT '(A)',
     &     'KMCcommit> (no KMC was run, so all values are for steady-state treatment of intervening minima'
      TRUEK=0.0D0
      SSK=0.0D0
      IF (DIRECTION.EQ.'BA') THEN
         PRINT '(A)','KMCcommit>  A minimum   P_(Ba) committor prob        k(B<-a)          k(B<-a) SS           ratio'
         DO J1=1,NMINA
            LJ1=LOCATIONA(J1)
            IF (NCOLSAVE(LJ1).GT.0) THEN
               WRITE(*,'(A,I8,4X,4G20.10)') 'KMCcommit> ',MINMAP(LJ1),PFOLD(LJ1),PFOLD(LJ1)*EXP(PFMIN(LJ1)-PFTOTALA)/WAITAB(LJ1),
     &                                   PFOLD(LJ1)*EXP(PFMIN(LJ1)-PFTOTALA)/EMKSUM(LJ1),
     &                                   EMKSUM(LJ1)/WAITAB(LJ1)
               TRUEK=TRUEK+PFOLD(LJ1)*EXP(PFMIN(LJ1)-PFTOTALA)/WAITAB(LJ1)
               SSK=SSK+PFOLD(LJ1)*EXP(PFMIN(LJ1)-PFTOTALA)/EMKSUM(LJ1)
            ENDIF
         ENDDO
         PRINT *
         PRINT '(3(A,G20.10))','KMCcommit>  kNSS k(B<-A)=',TRUEK,' kSS k(B<-A)=',SSK,' ratio ',TRUEK/SSK
         PRINT '(2(A,G20.10))','KMCcommit>  from detailed balance kNSS: k(A<-B)=',TRUEK*EXP(PFTOTALA-PFTOTALB),
     &                         ' from KMC kNSS: k(B<-A)=',TRUEK
      ELSE
         PRINT'(A)','KMCcommit>  B minimum   P_(Ab) committor prob         k(A<-b)          k(A<-b) SS           ratio'
         DO J1=1,NMINB
            LJ1=LOCATIONB(J1)
            IF (NCOLSAVE(LJ1).GT.0) THEN
               WRITE(*,'(A,I8,4X,4G20.10)') 'KMCcommit> ',MINMAP(LJ1),PFOLD(LJ1),PFOLD(LJ1)*EXP(PFMIN(LJ1)-PFTOTALB)/WAITAB(LJ1),
     &                                      PFOLD(LJ1)*EXP(PFMIN(LJ1)-PFTOTALB)/EMKSUM(LJ1),
     &                                      EMKSUM(LJ1)/WAITAB(LJ1)
               TRUEK=TRUEK+PFOLD(LJ1)*EXP(PFMIN(LJ1)-PFTOTALB)/WAITAB(LJ1)
               SSK=SSK+PFOLD(LJ1)*EXP(PFMIN(LJ1)-PFTOTALB)/EMKSUM(LJ1)
            ENDIF
         ENDDO
         PRINT *
         PRINT '(3(A,G20.10))','KMCcommit>  kNSS k(A<-B)=',TRUEK,' kSS k(A<-B)=',SSK,' ratio ',TRUEK/SSK
         PRINT '(2(A,G20.10))','KMCcommit>  from KMC kNSS: k(A<-B)=',TRUEK,
     &                         ' from detailed balance kNSS: k(B<-A)=',TRUEK*EXP(PFTOTALB-PFTOTALA)
      ENDIF
!
!!!!!!!!!!!!!!!!!!!   end of P^fold calculation  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      DEALLOCATE(DVEC,COL_IND)
!
!  The minima may be in a different order if regrouping has occurred.
!
      OPEN(UNIT=1,FILE='commit.data.regrouped',STATUS='UNKNOWN')
      WRITE(1,'(G20.10)') PFOLD(1:NMIN)
      CLOSE(1)
      IF (REGROUPT) THEN
         NMINA=NMINASAVE
         NMINB=NMINBSAVE
         DEALLOCATE(LOCATIONA,LOCATIONB)
         ALLOCATE(LOCATIONA(NMINA),LOCATIONB(NMINB))
         LOCATIONA(1:NMINA)=LOCATIONASAVE(1:NMINA)
         LOCATIONB(1:NMINB)=LOCATIONBSAVE(1:NMINB)
      ENDIF
 
      RETURN
      END
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  Construct DMATMC
!
!  NCOL(M2)     = # connections for minimum M2
!  NVAL(J1,M2)  = index of minimum involved in connection J1 from minimum M2
!  DMAT0(J1,M2) = KMC-type probability of taking connection J1 from minimum 
!                 M2 to minimum NVAL(J1,M2)
!
!  Degenerate rearrangements are excluded.
!
      SUBROUTINE MAKED(DMAT0,NCOL,NVAL,DEADTS,SORTT,ISA,ISB,KSUM)
      USE COMMONS
      IMPLICIT NONE
      INTEGER NCOL(NMIN), J1, M1, M2, J2, NVAL(NCONNMAX,NMIN), NVALTMP(NCONNMAX), OTHER
      DOUBLE PRECISION DMATTMP(NCONNMAX), DUMMY, KSUM(NMIN)
      LOGICAL DEADTS(NTS), MATCHED, SORTT, ALLOWED, ISA(NMIN), ISB(NMIN)
C     INTEGER, PARAMETER :: extra = selected_real_kind(20,200)
C     REAL(KIND=extra) :: LDUMMY, DMAT0(NCONNMAX,NMIN)
      DOUBLE PRECISION :: DMAT0(NCONNMAX,NMIN)

!
!  We are setting up DMAT0(thing,M2)
!  Cycle over connected transition states using pointers for minimum M2.
!
      NCOL(1:NMIN)=0

      fromloop: DO M2=1,NMIN   
         IF ((DIRECTION.EQ.'AB').AND.ISA(M2).AND.(.NOT.SORTT)) THEN
            CYCLE fromloop ! no escape from A
         ELSEIF ((DIRECTION.EQ.'BA').AND.ISB(M2).AND.(.NOT.SORTT)) THEN
            CYCLE fromloop ! no escape from B
         ENDIF
         J1=TOPPOINTER(M2)  !  sets J1 to the ts connected to minimum M2 with the highest value
         IF (J1.LE.0) CYCLE fromloop
11       CONTINUE
         IF ((.NOT.DEADTS(J1)).AND.(PLUS(J1).NE.MINUS(J1))) THEN
            IF (PLUS(J1).EQ.M2) THEN
               OTHER=MINUS(J1)
            ELSE 
               OTHER=PLUS(J1)
            ENDIF
            MATCHED=.FALSE.
            ALLOWED=.TRUE.
            IF (.NOT.SORTT) THEN ! assume this is the real KMC run if SORTT is .TRUE.
               IF ((DIRECTION.EQ.'BA').AND.ISA(OTHER)) THEN
                  ALLOWED=.FALSE. ! no contribution from transitions into A minima
               ELSEIF ((DIRECTION.EQ.'AB').AND.ISB(OTHER)) THEN
                  ALLOWED=.FALSE. ! no contribution from transitions into B minima
               ENDIF
            ENDIF
            IF (ALLOWED) THEN
               IF (PLUS(J1).EQ.M2) THEN  !  M2 M1  
                  matchcolp: DO M1=1,NCOL(M2)
                     IF (NVAL(M1,M2).EQ.MINUS(J1)) THEN ! a previous ts also links this pair
                        DMAT0(M1,M2)=MIN(DMAT0(M1,M2)+EXP(KPLUS(J1)-KSUM(PLUS(J1))),1.0D0)
                        MATCHED=.TRUE.
                        EXIT matchcolp
                     ENDIF
                  ENDDO matchcolp
                  IF (.NOT.MATCHED) THEN ! this minimum has not been connected to from M1 before
                     NCOL(M2)=NCOL(M2)+1
                     NVAL(NCOL(M2),M2)=MINUS(J1)
                     DMAT0(NCOL(M2),M2)=MIN(EXP(KPLUS(J1)-KSUM(PLUS(J1))),1.0D0)
                  ENDIF
!                 IF (DEBUG) PRINT '(A,3I6,G20.10)','KMCcommit> M2,NCOL,NVAL,DMAT0=',
!    &                                               M2,NCOL(M2),NVAL(NCOL(M2),M2),DMAT0(NCOL(M2),M2)
               ELSE IF (MINUS(J1).EQ.M2) THEN  !  M1 M2 
                  matchcolm: DO M1=1,NCOL(M2)
                     IF (NVAL(M1,M2).EQ.PLUS(J1)) THEN ! a previous ts also links this pair
                        DMAT0(M1,M2)=MIN(DMAT0(M1,M2)+EXP(KMINUS(J1)-KSUM(MINUS(J1))),1.0D0)
                        MATCHED=.TRUE.
                        EXIT matchcolm
                     ENDIF
                  ENDDO matchcolm
                  IF (.NOT.MATCHED) THEN ! this minimum has not been connected to from M1 before
                     NCOL(M2)=NCOL(M2)+1
                     NVAL(NCOL(M2),M2)=PLUS(J1)
                     DMAT0(NCOL(M2),M2)=MIN(EXP(KMINUS(J1)-KSUM(MINUS(J1))),1.0D0)
                  ENDIF
!                 IF (DEBUG) PRINT '(A,3I6,G20.10)','KMCcommit> M2,NCOL,NVAL,DMAT0=',
!    &                                               M2,NCOL(M2),NVAL(NCOL(M2),M2),DMAT0(NCOL(M2),M2)
               ENDIF
            ENDIF
         ENDIF
         IF (PLUS(J1).EQ.M2) THEN
            J1=POINTERP(J1)
         ELSE IF (MINUS(J1).EQ.M2) THEN
            J1=POINTERM(J1)
         ENDIF
         IF (J1.GT.0) GOTO 11
      ENDDO fromloop
C
C  Reorder DMAT0 so that the largest entry appears first etc.
C  This should speed up KMC, since largest branching probabilities
C  are most likely to be accepted.
C
      IF (SORTT) THEN
         DO J1=1,NMIN
            IF (NCOL(J1).GT.1) THEN
               NVALTMP(1:NCOL(J1))=NVAL(1:NCOL(J1),J1)
               DMATTMP(1:NCOL(J1))=DMAT0(1:NCOL(J1),J1)
               CALL SORT(NCOL(J1),NCONNMAX,DMATTMP,NVALTMP)
               NVAL(1:NCOL(J1),J1)=NVALTMP(1:NCOL(J1))
               DMAT0(1:NCOL(J1),J1)=DMATTMP(1:NCOL(J1))
            ENDIF
         ENDDO
C
C  Check row normalisation.
C
         DO J1=1,NMIN
            DUMMY=0.0D0
            DO J2=1,NCOL(J1)
               DUMMY=DUMMY+DMAT0(J2,J1)
               IF (DEBUG) WRITE(*,'(A,3I6,3G20.10)') 'J1,J2,NVAL,DMAT0,sum=',J1,J2,NVAL(J2,J1),DMAT0(J2,J1),DUMMY
            ENDDO
            IF (DEBUG.AND.(NCOL(J1).GT.0)) WRITE(*,'(A,2I6,3G20.10)') 'J1,ncol,sum=',J1,NCOL(J1),DUMMY
            IF ((NCOL(J1).GT.0).AND.(ABS(DUMMY-1.0D0).GT.1.0D-10)) THEN
            PRINT*,'MAKED> ERROR - J1,NCOL(J1),DUMMY=',J1,NCOL(J1),DUMMY
                  STOP
            ENDIF
         ENDDO
      ELSE
C
C  Reorder DMAT0 so that the connected minima for minimum J1
C  appear in ascending order. Could speed up PFOLD calculation in 
C  new compressed row storage format. Doesn't seem to make much difference for LJ55.
C
         DO J1=1,NMIN
            IF (NCOL(J1).GT.1) THEN
               NVALTMP(1:NCOL(J1))=NVAL(1:NCOL(J1),J1)
               DMATTMP(1:NCOL(J1))=DMAT0(1:NCOL(J1),J1)
               CALL SORT4(NCOL(J1),NCONNMAX,DMATTMP,NVALTMP)
               NVAL(1:NCOL(J1),J1)=NVALTMP(1:NCOL(J1))
               DMAT0(1:NCOL(J1),J1)=DMATTMP(1:NCOL(J1))
            ENDIF
         ENDDO
      ENDIF

      RETURN
      END
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!  Construct DMATMC using an economical storage format
!
!  NCOL(M2)     = # connections for minimum M2
!  NVAL(J1,M2)  = index of minimum involved in connection J1 from minimum M2
!  DMAT0(J1,M2) = KMC-type probability of taking connection J1 from minimum 
!                 M2 to minimum NVAL(J1,M2)
!
!  Degenerate rearrangements are excluded.
!
      SUBROUTINE MAKEDLIN(DMAT0,NCOL,NOFF,NVAL,DEADTS,SORTT,ISA,ISB,KSUM)
      USE COMMONS
      IMPLICIT NONE
      INTEGER NCOL(NMIN), J1, M1, M2, J2, NVAL(2*NTS), NVALTMP(NCONNMAX), OTHER, NOFF(NMIN)
      DOUBLE PRECISION DMATTMP(NCONNMAX), DUMMY, KSUM(NMIN)
      LOGICAL DEADTS(NTS), MATCHED, SORTT, ALLOWED, ISA(NMIN), ISB(NMIN)
C     INTEGER, PARAMETER :: extra = selected_real_kind(20,200)
C     REAL(KIND=extra) :: LDUMMY, DMAT0(NCONNMAX,NMIN)
      DOUBLE PRECISION :: DMAT0(2*NTS)

!
!  We are setting up DMAT0(NOFF(M2),NOFF(M2)+NCOL(M2))
!  First entry for M2 at position NOFF(M2)+1
!  No entries if NCOL(M2)=0
!  Cycle over connected transition states using pointers for minimum M2.
!
      NCOL(1:NMIN)=0
      NOFF(1:NMIN)=0

      fromloop: DO M2=1,NMIN   
         IF (M2.GT.1) NOFF(M2)=NOFF(M2-1)+NCOL(M2-1)
         IF ((M2.GT.1).AND.DEBUG) PRINT '(A,4I10)','M2,NOFF(M2-1),NCOL(M2-1),NOFF(M2)=',M2,NOFF(M2-1),NCOL(M2-1),NOFF(M2)
         IF ((DIRECTION.EQ.'AB').AND.ISA(M2).AND.(.NOT.SORTT)) THEN
            CYCLE fromloop ! no escape from A
         ELSEIF ((DIRECTION.EQ.'BA').AND.ISB(M2).AND.(.NOT.SORTT)) THEN
            CYCLE fromloop ! no escape from B
         ENDIF
         J1=TOPPOINTER(M2)  !  sets J1 to the ts connected to minimum M2 with the highest value
         IF (J1.LE.0) CYCLE fromloop
11       CONTINUE
         IF ((.NOT.DEADTS(J1)).AND.(PLUS(J1).NE.MINUS(J1))) THEN
            IF (PLUS(J1).EQ.M2) THEN
               OTHER=MINUS(J1)
            ELSE 
               OTHER=PLUS(J1)
            ENDIF
            MATCHED=.FALSE.
            ALLOWED=.TRUE.
            IF (.NOT.SORTT) THEN ! assume this is the real KMC run if SORTT is .TRUE.
               IF ((DIRECTION.EQ.'BA').AND.ISA(OTHER)) THEN
                  ALLOWED=.FALSE. ! no contribution from transitions into A minima
               ELSEIF ((DIRECTION.EQ.'AB').AND.ISB(OTHER)) THEN
                  ALLOWED=.FALSE. ! no contribution from transitions into B minima
               ENDIF
            ENDIF
            IF (ALLOWED) THEN
               IF (PLUS(J1).EQ.M2) THEN  !  M2 M1  
                  matchcolp: DO M1=1,NCOL(M2)
                     IF (NVAL(NOFF(M2)+M1).EQ.MINUS(J1)) THEN ! a previous ts also links this pair
                        DMAT0(NOFF(M2)+M1)=MIN(DMAT0(NOFF(M2)+M1)+EXP(KPLUS(J1)-KSUM(PLUS(J1))),1.0D0)
                        MATCHED=.TRUE.
                        EXIT matchcolp
                     ENDIF
                  ENDDO matchcolp
                  IF (.NOT.MATCHED) THEN ! this minimum has not been connected to from M1 before
                     NCOL(M2)=NCOL(M2)+1
                     NVAL(NOFF(M2)+NCOL(M2))=MINUS(J1)
                     DMAT0(NOFF(M2)+NCOL(M2))=MIN(EXP(KPLUS(J1)-KSUM(PLUS(J1))),1.0D0)
                  ENDIF
!                 IF (DEBUG) PRINT '(A,3I6,G20.10)','KMCcommit> M2,NCOL,NVAL,DMAT0=',
!    &                                               M2,NCOL(M2),NVAL(NOFF(M2)+NCOL(M2)),DMAT0(NOFF(M2)+NCOL(M2)
               ELSE IF (MINUS(J1).EQ.M2) THEN  !  M1 M2 
                  matchcolm: DO M1=1,NCOL(M2)
                     IF (NVAL(NOFF(M2)+M1).EQ.PLUS(J1)) THEN ! a previous ts also links this pair
                        DMAT0(NOFF(M2)+M1)=MIN(DMAT0(NOFF(M2)+M1)+EXP(KMINUS(J1)-KSUM(MINUS(J1))),1.0D0)
                        MATCHED=.TRUE.
                        EXIT matchcolm
                     ENDIF
                  ENDDO matchcolm
                  IF (.NOT.MATCHED) THEN ! this minimum has not been connected to from M1 before
                     NCOL(M2)=NCOL(M2)+1
                     NVAL(NOFF(M2)+NCOL(M2))=PLUS(J1)
                     DMAT0(NOFF(M2)+NCOL(M2))=MIN(EXP(KMINUS(J1)-KSUM(MINUS(J1))),1.0D0)
                  ENDIF
!                 IF (DEBUG) PRINT '(A,3I6,G20.10)','KMCcommit> M2,NCOL,NVAL,DMAT0=',
!    &                                               M2,NCOL(M2),NVAL(NOFF(M2)+NCOL(M2)),DMAT0(NOFF(M2)+NCOL(M2))
               ENDIF
            ENDIF
         ENDIF
         IF (PLUS(J1).EQ.M2) THEN
            J1=POINTERP(J1)
         ELSE IF (MINUS(J1).EQ.M2) THEN
            J1=POINTERM(J1)
         ENDIF
         IF (J1.GT.0) GOTO 11
      ENDDO fromloop
C
C  Check row normalisation.
C
      DO J1=1,NMIN
         DUMMY=0.0D0
         DO J2=1,NCOL(J1)
            DUMMY=DUMMY+DMAT0(NOFF(J1)+J2)
            IF (DEBUG) WRITE(*,'(A,3I8,3G20.10)') 'J1,J2,NVAL,DMAT0,sum=',J1,J2,NVAL(NOFF(J1)+J2),DMAT0(NOFF(J1)+J2),DUMMY
         ENDDO
         IF (DEBUG.AND.(NCOL(J1).GT.0)) WRITE(*,'(A,2I8,3G20.10)') 'J1,ncol,sum=',J1,NCOL(J1),DUMMY
         IF ((NCOL(J1).GT.0).AND.(ABS(DUMMY-1.0D0).GT.1.0D-10)) THEN
         PRINT*,'MAKEDLIN> ERROR - J1,NCOL(J1),DUMMY=',J1,NCOL(J1),DUMMY
               STOP
         ENDIF
      ENDDO

      RETURN
      END
