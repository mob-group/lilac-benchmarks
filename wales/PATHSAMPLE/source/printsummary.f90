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
! Printing some summary info about the KTN, to standard out 

      SUBROUTINE PRINTSUMMARY
      USE COMMONS, ONLY : EMIN, ETS, NMIN, NTS, NCONN, PLUS, MINUS, NCONNMIN, TSTHRESH, MAXBARRIER, MINBARRIER, &
                          KPLUS, KMINUS
      IMPLICIT NONE

      INTEGER :: NDEAD, MINTARG, NLEFTTS, NCOUNT, NDONE, J1, NLEFTMIN, MAXTSINDEX, MINTSINDEX, CURRTARG
      INTEGER :: NDISTA(NMIN)
      DOUBLE PRECISION :: MAXTS, MINTS, GMIN, CURRENERGY, LOWBARRIER, HIGHBARRIER
      DOUBLE PRECISION :: BARRIER_MAX, BARRIER_MIN
      DOUBLE PRECISION :: CUT_UNDERFLOW=-300.0D0 ! set here but NOT used in the current coding of checkTS
      LOGICAL :: CHANGED
      LOGICAL :: DEADTS(NTS), INCC(NMIN)

      WRITE(*,'(A)')'*************************************************************************************************'
      WRITE(*,'(A)')'printsummary> Network summary information'
      WRITE(*,'(A)')'using thresholds:'
      WRITE(*,'(A,I8)')'NCONNMIN=',NCONNMIN
      WRITE(*,'(A,G20.10)')'MAXBARRIER=',MAXBARRIER
      WRITE(*,'(A,G20.10)')'MINBARRIER=',MINBARRIER
      WRITE(*,'(A,G20.10)')'TSTHRESH or MAXTSENERGY=',TSTHRESH
      WRITE(*,'(A)')'and NOT including degenerate rearrangements.'

      CALL GETNCONN ! GETNCONN only counts non-degenerate connections
      NDEAD=0
      IF (NCONNMIN.GE.0) THEN
         DO J1=1,NMIN
            IF (NCONN(J1).LE.NCONNMIN) THEN
               NDEAD=NDEAD+1
            ENDIF 
         ENDDO
      ENDIF
      WRITE(*,'(A,I8,A,I8)')'There are ',NMIN-NDEAD,' valid minima out of ',NMIN

      GMIN=HUGE(1.0D0)
      DO J1=1,NMIN
         IF (EMIN(J1).LT.GMIN .AND. (NCONN(J1).GT.NCONNMIN)) THEN
            GMIN=EMIN(J1)
            MINTARG=J1
         ENDIF
      ENDDO
      WRITE(*,'(A,I8,A,G20.10)')'Valid global minimum at index ',MINTARG,' with energy ',GMIN

      NLEFTTS=0
      DEADTS(1:NTS)=.FALSE.
      DO J1=1,NTS
         CALL CHECKTS(ETS(J1),EMIN(PLUS(J1)),EMIN(MINUS(J1)),KPLUS(J1),KMINUS(J1),NCONN(PLUS(J1)),NCONN(MINUS(J1)), &
                      PLUS(J1),MINUS(J1),.FALSE.,CUT_UNDERFLOW,DEADTS(J1))
         IF (.NOT.DEADTS(J1)) NLEFTTS=NLEFTTS+1
      ENDDO
      WRITE(*,'(A,I8,A,I8)')'There are ',NLEFTTS,' valid transition states out of ',NTS

      IF ((NMIN-NDEAD.EQ.0).OR.(NLEFTTS.EQ.0)) RETURN ! sanity check

      MINTS=HUGE(1.0D0)
      MAXTS=-HUGE(1.0D0)
      DO J1=1,NTS
         IF (DEADTS(J1)) CYCLE
         IF (ETS(J1).LE.MINTS) THEN
            MINTS=ETS(J1)
            MINTSINDEX=J1
         ENDIF
         IF (ETS(J1).GE.MAXTS) THEN
            MAXTS=ETS(J1)
            MAXTSINDEX=J1
         ENDIF
      ENDDO
      WRITE(*,'(A,I8,A,G20.10)')'Highest valid transition state at index ',MAXTSINDEX,' with energy ',MAXTS
      WRITE(*,'(A,I8,A,G20.10)')'Lowest valid transition state at index ',MINTSINDEX,' with energy ',MINTS

      LOWBARRIER=HUGE(1.0D0)
      HIGHBARRIER=-HUGE(1.0D0)
      DO J1=1,NTS
         IF (DEADTS(J1)) CYCLE
         BARRIER_MAX = MAX(ETS(J1)-EMIN(PLUS(J1)),ETS(J1)-EMIN(MINUS(J1)))
         BARRIER_MIN = MIN(ETS(J1)-EMIN(PLUS(J1)),ETS(J1)-EMIN(MINUS(J1)))
         IF (BARRIER_MIN.LE.LOWBARRIER) THEN
            LOWBARRIER=BARRIER_MIN
            MINTSINDEX=J1
         ENDIF
         IF (BARRIER_MAX.GE.HIGHBARRIER) THEN
            HIGHBARRIER=BARRIER_MAX
            MAXTSINDEX=J1
         ENDIF
      ENDDO
      WRITE(*,'(A,G20.10,A,I8)')'Highest valid barrier ',HIGHBARRIER,' for TS index ',MAXTSINDEX
      WRITE(*,'(A,G20.10,A,I8)')'Lowest valid barrier ',LOWBARRIER,' for TS index ',MINTSINDEX

      WRITE(*,'(A)')'Connected components for valid minima and TSs:'
      WRITE(*,'(A)')'CC #        # minima    index of lowest min    energy of lowest min'
      NCOUNT=0
      NDONE=0
      CURRTARG=MINTARG ! start with the valid global minimum 
      CURRENERGY=GMIN
      INCC(1:NMIN)=.FALSE.
      DO 
         NCOUNT=NCOUNT+1
         NDISTA(1:NMIN)=1000000
         NDISTA(CURRTARG)=0

5        CHANGED=.FALSE.
         DO J1=1,NTS
            IF (DEADTS(J1)) CYCLE
            IF (NDISTA(MINUS(J1))+1.LT.NDISTA(PLUS(J1))) THEN
               CHANGED=.TRUE.
               NDISTA(PLUS(J1))=NDISTA(MINUS(J1))+1
            ENDIF
            IF (NDISTA(PLUS(J1))+1.LT.NDISTA(MINUS(J1))) THEN
               CHANGED=.TRUE.
               NDISTA(MINUS(J1))=NDISTA(PLUS(J1))+1
            ENDIF
         ENDDO
         IF (CHANGED) GOTO 5

         NLEFTMIN=0
         DO J1=1,NMIN
            IF (INCC(J1)) CYCLE
            IF ((NDISTA(J1).LT.1000000).AND.(NCONN(J1).GT.NCONNMIN)) THEN
               NLEFTMIN=NLEFTMIN+1
               NDONE=NDONE+1
               INCC(J1)=.TRUE.
            ENDIF
         ENDDO

         WRITE(*,'(2(I8,4X),I8,15X,G20.10)')NCOUNT,NLEFTMIN,CURRTARG,CURRENERGY
         IF (NDONE.EQ.NMIN-NDEAD) EXIT

         ! update CURRTARG and CURRENERGY: find the lowest minimum not yet in a CC
         CURRENERGY=HUGE(1.0D0)
         DO J1=1,NMIN
            IF (INCC(J1)) CYCLE
            IF (EMIN(J1).LT.CURRENERGY .AND. (NCONN(J1).GT.NCONNMIN)) THEN
               CURRENERGY=EMIN(J1)
               CURRTARG=J1
            ENDIF
         ENDDO
      ENDDO

      WRITE(*,'(A)')'*************************************************************************************************'

      RETURN
      END SUBROUTINE PRINTSUMMARY
