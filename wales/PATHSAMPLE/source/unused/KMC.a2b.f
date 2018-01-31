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
C  A to B
C
      SUBROUTINE KMCA2B
      IMPLICIT NONE
      INCLUDE 'common.h'
      INTEGER J1, J2, J3, NCYCLE, NKMC, ATMIN, MIN2, NDEAD
      LOGICAL CHANGED, DEADTS(NTS)
      DOUBLE PRECISION DMAT0(NCONNMAX,NMIN), RANDOM, WAIT, WAIT2, WAITMEAN, NSTEP, NSTEPMEAN, DPRAND
      INTEGER NCOL(NMIN), NVAL(NCONNMAX,NMIN), NUNCON, ATMINP
      DOUBLE PRECISION TOTALAB, TOTALBA, DENOM,NDIST(NMIN),P12P21, PFTOTALACONN, EMKSUM(NMIN), DMIN, DMAX

      PRINT*,'WARNING - removing dead ends from the calculation'
      WRITE(*,'(A,1A,A,1A)') 'Running KMC in the direction ',DIRECTION(1:1),'<-',DIRECTION(2:2)
      WRITE(*,'(A,G20.10)') 'threshold for combining minima is p12*p21>',PAIRTHRESH
      WRITE(*,'(A,G20.10)') 'exponent factor for rewighting (0 means no reweighting): ',EXPFAC

      NUNCON=0
      PFTOTALACONN=0.0D0
      DO J1=1,NMINA
         IF (NCONN(J1).EQ.0) THEN
            NUNCON=NUNCON+1
         ELSE IF (NCONN(J1).GT.NCONNMIN) THEN
            PFTOTALACONN=PFTOTALACONN+EXP(PFMIN(J1)) ! PFMIN has already been adjusted by PFMEAN in setup.f
         ENDIF
      ENDDO
      PFTOTALACONN=LOG(PFTOTALACONN)
      PRINT*,NUNCON,' A minima are unconnected and will not be considered as starting points'
      NDEAD=0
      DO J1=1,NMINA
         IF (NCONN(J1).LE.NCONNMIN) THEN
            NDEAD=NDEAD+1
            IF (DEBUG) PRINT*,'discarding minimum ',J1,' with ',NCONN(J1),' connections'
         ENDIF
      ENDDO
      DO J1=NMINA+NMINB+1,NMIN
         IF (NCONN(J1).LE.NCONNMIN) THEN
            NDEAD=NDEAD+1
            IF (DEBUG) PRINT*,'discarding minimum ',J1,' with ',NCONN(J1),' connections'
         ENDIF
      ENDDO
      PRINT*,NDEAD,' minima with ',NCONNMIN,' connections or fewer will not be considered'
C
C  Flag transition states to dead-end minima as DEAD and correct KSUM values.
C  NCONN only counted non-degenerate rearrangements as connections.
C
      DO J1=1,NMIN
         KSUM(J1)=0.0D0
      ENDDO
      DO J1=1,NTS
         IF ((NCONN(PLUS(J1)).LE.NCONNMIN).OR.(NCONN(MINUS(J1)).LE.NCONNMIN)) THEN
            DEADTS(J1)=.TRUE.
!        ELSEIF (KPLUS(J1)-KMEAN.LT.-300.0D0) THEN  ! to prevent underflow for transition states with daft energies
!           DEADTS(J1)=.TRUE.
!        ELSEIF (KMINUS(J1)-KMEAN.LT.-300.0D0) THEN ! to prevent underflow for transition states with daft energies
!           DEADTS(J1)=.TRUE.
         ELSE
            DEADTS(J1)=.FALSE.
            IF (KPLUS(J1)-KMEAN.GT.-300.0D0) THEN
               IF (PLUS(J1).NE.MINUS(J1)) LKSUM(PLUS(J1))=LKSUM(PLUS(J1))+EXP(KPLUS(J1)-KMEAN)
            ENDIF 
            IF  (KMINUS(J1)-KMEAN.GT.-300.0D0) THEN
               IF (PLUS(J1).NE.MINUS(J1)) LKSUM(MINUS(J1))=LKSUM(MINUS(J1))+EXP(KMINUS(J1)-KMEAN)
            ENDIF 
         ENDIF
      ENDDO
      DO J1=1,NMIN
         IF (KSUM(J1).GT.0.0D0) THEN
            KSUM(J1)=LOG(KSUM(J1))+KMEAN
         ENDIF
      ENDDO

      DO J1=1,NMIN
         EMKSUM(J1)=EXP(-KSUM(J1))
      ENDDO

      CALL MAKEDMAT3(DMAT0,NCOL,NVAL,DEADTS)
C
C  Calculate minimum number of steps of each minimum from the B set.
C
      DO J1=1,NMIN
         NDIST(J1)=1000000.0D0
      ENDDO
      DO J1=NMINA+1,NMINA+NMINB
         NDIST(J1)=0.0D0
      ENDDO
      NCYCLE=0
5     CHANGED=.FALSE.
      NCYCLE=NCYCLE+1
      DMIN=100000.0D0
      DMAX=0.0D0
      NUNCON=0
      DO J1=1,NMINA
         DO J2=1,NCOL(J1)
            IF (NDIST(NVAL(J2,J1))+1.0D0.LT.NDIST(J1)) THEN
               CHANGED=.TRUE.
               NDIST(J1)=NDIST(NVAL(J2,J1))+1.0D0
            ENDIF
         ENDDO
         IF ((NDIST(J1).GT.DMAX).AND.(NCONN(J1).NE.0)) DMAX=NDIST(J1)
         IF (NDIST(J1).LT.DMIN) DMIN=NDIST(J1)
         IF (NDIST(J1).EQ.1000000.0) NUNCON=NUNCON+1
      ENDDO
      DO J1=NMINA+NMINB+1,NMIN
         DO J2=1,NCOL(J1)
            IF (NDIST(NVAL(J2,J1))+1.0D0.LT.NDIST(J1)) THEN
               CHANGED=.TRUE.
               NDIST(J1)=NDIST(NVAL(J2,J1))+1.0D0
            ENDIF
         ENDDO
         IF ((NDIST(J1).GT.DMAX).AND.(NCONN(J1).NE.0)) DMAX=NDIST(J1)
         IF (NDIST(J1).LT.DMIN) DMIN=NDIST(J1)
         IF (NDIST(J1).EQ.1000000.0) NUNCON=NUNCON+1
      ENDDO
      IF (CHANGED) GOTO 5
      PRINT*,'number of steps to B region calculation converged in ',NCYCLE-1,' cycles'
      PRINT*,'maximum number of steps=',DMAX
      PRINT*,'number of minima not connected to B region=',NUNCON
      IF (DEBUG) THEN
         DO J1=1,NMIN
            PRINT*,'J1,NDIST=',J1,NDIST(J1)
         ENDDO
      ENDIF

C     DO J1=NMINA+1,NMIN
C        WRITE(*,'(A,I5,3G20.10)') 'J1,NDIST,wait,prob=',J1,NDIST(J1),EMKSUM(J1),EXP(PFMIN(J1)-PFTOTALB)
C     ENDDO
C
C  KMC loop
C
      NKMC=0
      WAITMEAN=0.0D0
      NSTEPMEAN=0.0D0
      WAIT2=0.0D0
      DENOM=0.0D0
10    NKMC=NKMC+1
      WAIT=0.0D0
      NSTEP=0.0D0
C
C  Choose a connected A minimum with more than one connection. 
C
C     CALL RANDOM_NUMBER(RANDOM)
      RANDOM=DPRAND()
      DUMMY=0.0D0
      DO J1=1,NMINA
         IF (NCONN(J1).GT.NCONNMIN) THEN
            DUMMY=DUMMY+EXP(PFMIN(J1)-PFTOTALACONN)
C           PRINT*,'J1,NDIST,DUMMY,RANDOM,PFMIN,PFTOTALBCONN=',J1,NDIST(J1),DUMMY,RANDOM,PFMIN(J1),PFTOTALBCONN
            IF (DUMMY.GE.RANDOM) THEN
C              TRUEPROB=PFMIN(J1)-PFTOTALBCONN  ! new
               ATMIN=J1
               ATMINP=ATMIN
               PRINT *,'A minimum ',J1,'chosen'
               CALL FLUSH(6)
               GOTO 20
            ENDIF
         ENDIF
      ENDDO
      PRINT*,'ERROR, J1,DUMMY,RANDOM=',J1,DUMMY,RANDOM
      STOP
20    NSTEP=NSTEP+1.0D0
C     PRINT*,'NSTEP,ATMIN,WAIT=',NSTEP,ATMIN,WAIT
c      IF (WAIT.GT.1000*WAITMEAN/DENOM) THEN ! this test shouldn;t have exp(-weight) times wait
c         PRINT*,'contribution for KMC run ',NKMC,' step ',NSTEP,' is more than 1000 x greater than the mean'
c         PRINT*,'contribution is ',WAIT,' mean is ',WAITMEAN/DENOM
c         PRINT*,'NSTEP,ATMIN,WAIT,ATMINP=',NSTEP,ATMIN,WAIT,ATMINP
cC        STOP
c      ENDIF
      DUMMY=0.0D0
C     CALL RANDOM_NUMBER(RANDOM)
      RANDOM=DPRAND()
      DO J1=1,NCOL(ATMIN)
         DUMMY=DUMMY+DMAT0(J1,ATMIN)
         IF (DUMMY.GT.RANDOM) THEN    
            MIN2=NVAL(J1,ATMIN)
C
C  In this branch we catch pairs of minima separated by low barriers where the trajectory
C  tends to be caught.
C
            IF (((MIN2.LE.NMINA).OR.(MIN2.GT.NMINA+NMINB)).AND.(DMAT0(J1,ATMIN).GT.PAIRTHRESH)) THEN
C              PRINT *,'In this branch we catch pairs of minima separated by low barriers' ; CALL FLUSH(6) ! jmc testing
C
C  Find the probability of going back to ATMIN from the next minimum.
C
               DO J2=1,NCOL(MIN2)
C                 PRINT*,'ATMIN,MIN2,NVAL(J2,MIN2)=',ATMIN,MIN2,NVAL(J2,MIN2)
                  IF (NVAL(J2,MIN2).EQ.ATMIN) THEN
                     P12P21=DMAT0(J1,ATMIN)*DMAT0(J2,MIN2)
                     GOTO 25
                  ENDIF
               ENDDO
               PRINT*,'ERROR, ATMIN,J1,MIN2,NCOL(MIN2)=',ATMIN,J1,MIN2,NCOL(MIN2)
               STOP
C
C  Choose one of the minima connected to ATMIN or MIN2 with renormalised probabilities
C  and waiting times.
C
25             IF (P12P21.GT.PAIRTHRESH) THEN
                  PRINT *,'In this branch we catch pairs of minima separated by low barriers' ; CALL FLUSH(6) ! jmc testing
C                 CALL RANDOM_NUMBER(RANDOM)
                  RANDOM=DPRAND()
                  DUMMY=0.0D0
                  DO J3=1,NCOL(ATMIN)
                     IF (NVAL(J3,ATMIN).NE.MIN2) THEN
                        DUMMY=DUMMY+DMAT0(J3,ATMIN)/(1.0D0-P12P21)
C                       PRINT*,'J3,NCOL(ATMIN),DUMMY,RANDOM=',J3,NCOL(ATMIN),DUMMY,RANDOM
                        IF (DUMMY.GT.RANDOM) THEN
                           WAIT=WAIT+(EMKSUM(ATMIN)+P12P21*EMKSUM(MIN2))/(1.0D0-P12P21)
C                          WRITE(*,'(A,2I5,3G20.10)') 
C    1                      'ATMIN,MIN2,P12P21,exp1,exp2=',ATMIN,MIN2,P12P21,EMKSUM(ATMIN),EMKSUM(MIN2)
                           ATMINP=ATMIN
                           ATMIN=NVAL(J3,ATMIN)
                           IF ((ATMIN.LE.(NMINA+NMINB)) .AND. (ATMIN.GT.NMINA)) GOTO 30
                           GOTO 20
                        ENDIF
                     ENDIF
                  ENDDO
                  DO J3=1,NCOL(MIN2)
                     IF (NVAL(J3,MIN2).NE.ATMIN) THEN
                        DUMMY=DUMMY+DMAT0(J3,MIN2)*DMAT0(J1,ATMIN)/(1.0D0-P12P21)
C                       PRINT*,'J3,NCOL(MIN2),DUMMY,RANDOM=',J3,NCOL(MIN2),DUMMY,RANDOM
                        IF (DUMMY.GT.RANDOM) THEN
                           WAIT=WAIT+(EMKSUM(ATMIN)+EMKSUM(MIN2))/(1.0D0-P12P21)
C                          WRITE(*,'(A,2I5,3G20.10)') 
C    1                      'ATMIN,MIN2,P12P21,exp1,exp2=',ATMIN,MIN2,P12P21,EMKSUM(ATMIN),EMKSUM(MIN2)
                           ATMINP=ATMIN
                           ATMIN=NVAL(J3,MIN2)
                           IF ((ATMIN.LE.(NMINA+NMINB)) .AND. (ATMIN.GT.NMINA)) GOTO 30
                           GOTO 20
                        ENDIF
                     ENDIF
                  ENDDO
                  PRINT*,'ERROR, ATMIN,J1,NCOL,DUMMY,RANDOM=',ATMIN,J1,NCOL(ATMIN),DUMMY,RANDOM
               ELSE
                  WAIT=WAIT+EMKSUM(ATMIN)
                  ATMINP=ATMIN
                  ATMIN=NVAL(J1,ATMIN)
                  IF ((ATMIN.LE.(NMINA+NMINB)) .AND. (ATMIN.GT.NMINA)) GOTO 30
                  GOTO 20
               ENDIF
            ELSE
               WAIT=WAIT+EMKSUM(ATMIN) 
C              TRUEPROB=TRUEPROB+LOG(DMAT0(J1,ATMIN))  ! new
               ATMINP=ATMIN
               ATMIN=NVAL(J1,ATMIN)
               IF ((ATMIN.LE.(NMINA+NMINB)) .AND. (ATMIN.GT.NMINA)) GOTO 30
               GOTO 20
            ENDIF
         ENDIF
      ENDDO
      PRINT*,'ERROR, ATMIN,J1,NCOL,DUMMY,RANDOM=',ATMIN,J1,NCOL(ATMIN),DUMMY,RANDOM
      STOP
30    WAITMEAN=WAITMEAN+WAIT
      WAIT2=WAIT2+WAIT**2
      NSTEPMEAN=NSTEPMEAN+NSTEP
      DENOM=DENOM+1.0D0
      WRITE(*,'(A,G12.4,A,G20.10,A,2G15.5,A,G15.5)') 
     1         'B region hit in ',NSTEP,' steps, time =',WAIT,' mean and sd=',
     2         WAITMEAN/DENOM,SQRT(MAX(WAIT2/DENOM-(WAITMEAN/DENOM)**2,0.0D0))
      IF (MOD(NKMC,100).EQ.0) WRITE(*,'(A,G20.10,A,G20.10,A,G20.10))') 
     1    ' a<-b ',DENOM*EXP(PFTOTALA-PFTOTALB)/WAITMEAN,' b<-a ',DENOM/WAITMEAN,' total ',
     2    DENOM*(1.0D0+EXP(PFTOTALA-PFTOTALB))/WAITMEAN
      CALL FLUSH(6)
      IF (NKMC.LT.NKMCCYCLES) GOTO 10
      WRITE(*,'(A,I10,A)') 'statistics for ',NKMC,' KMC steps:'
      WRITE(*,'(A,2G20.10)') 'mean and sd of first encounter time=',WAITMEAN/DENOM,SQRT(MAX(WAIT2/DENOM-(WAITMEAN/DENOM)**2,0.0D0))
      WRITE(*,'(A,G20.10)') 'mean number of steps=',NSTEPMEAN/NKMC

      TOTALBA=DENOM/WAITMEAN
      TOTALAB=TOTALBA*EXP(PFTOTALA-PFTOTALB)
      PRINT*,'convert waiting times to rates:'
      WRITE(*,'(A,G20.10,A,G20.10,A,G20.10))') ' a<-b ',TOTALAB,' b<-a ',TOTALBA,' total ',TOTALAB+TOTALBA

      RETURN
      END
