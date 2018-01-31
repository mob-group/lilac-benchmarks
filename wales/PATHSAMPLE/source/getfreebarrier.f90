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

!
!  Find the approximate lowest barrier from each minimum to a product minimum
!  using a superbasin analysis. This routine is designed to work with free 
!  energy groups and transition states from regroupfree2.
!
SUBROUTINE GETFREEBARRIER(BARRIER,NGROUPS,GLOBALMIN,NEWEMIN,NTS,NEWETS,NCONNGROUP,GROUPTS,STARTGROUP,GROUPA,GROUPB,EINC,NMIN, &
  &                       DIRECTION,DEBUG,GROUPCONN)
IMPLICIT NONE
INTEGER NMIN, NGROUPS, NTS
DOUBLE PRECISION HIGHESTTS, ETHRESH, BARRIER(NMIN), NEWEMIN(NMIN), NEWETS(NTS), LOCALETS, EINC, GLOBALMIN
INTEGER NCONNGROUP(NMIN), LPLUS, LMINUS, GROUPCONN(2*NTS)
LOGICAL GROUPA(NGROUPS), GROUPB(NGROUPS)
INTEGER J1, BASIN(NGROUPS), NBASIN, J2, GROUPTS(2*NTS), STARTGROUP(NMIN), NGROUPA, NGROUPB, LOCATIONA(NGROUPS), LOCATIONB(NGROUPS)
LOGICAL CHANGED, BASINT(NGROUPS), DEBUG
CHARACTER(LEN=2) DIRECTION

NGROUPA=0
NGROUPB=0
DO J1=1,NGROUPS
   IF (GROUPA(J1)) THEN
      NGROUPA=NGROUPA+1
      LOCATIONA(NGROUPA)=J1
   ELSEIF (GROUPB(J1)) THEN
      NGROUPB=NGROUPB+1
      LOCATIONB(NGROUPB)=J1
   ENDIF
ENDDO

BARRIER(1:NGROUPS)=-1.0D0
!
!  Find the highest transition state.
!
PRINT '(A,G20.10)','getfreebarrier> lowest minimum in product set lies at ',GLOBALMIN
HIGHESTTS=-1.0D100
DO J1=1,NGROUPS
   IF (NCONNGROUP(J1).EQ.0) CYCLE
   DO J2=1,NCONNGROUP(J1)
      LOCALETS=NEWETS(GROUPTS(STARTGROUP(J1)+J2-1))
      IF (LOCALETS.GT.HIGHESTTS) THEN
!
! Avoid infinite superbasin analysis for rates that are zero in regroupfree2
! where the ts energy is set ot HUGE(1.0D0).
!
         IF (LOCALETS.LT.HUGE(1.0D0)/10.0D0) HIGHESTTS=LOCALETS
      ENDIF
   ENDDO
ENDDO
PRINT '(A,G20.10)','getfreebarrier> highest transition state lies at ',HIGHESTTS

ETHRESH=GLOBALMIN

DO
   BASIN(1:NGROUPS)=0
   NBASIN=0
   DO 
      CHANGED=.FALSE.
      DO J1=1,NGROUPS
         IF (NCONNGROUP(J1).EQ.0) CYCLE
         LPLUS=J1
         DO J2=1,NCONNGROUP(J1)
            IF (NEWETS(GROUPTS(STARTGROUP(J1)+J2-1)).LT.ETHRESH) THEN
               LMINUS=GROUPCONN(STARTGROUP(J1)+J2-1)
               IF ((BASIN(LPLUS).EQ.0).AND.(BASIN(LMINUS).EQ.0)) THEN
                  CHANGED=.TRUE.
                  NBASIN=NBASIN+1
                  BASIN(LPLUS)=NBASIN
                  BASIN(LMINUS)=NBASIN
               ELSEIF (BASIN(LPLUS).NE.BASIN(LMINUS)) THEN
                  CHANGED=.TRUE.
                  IF (BASIN(LPLUS).EQ.0) THEN
                     BASIN(LPLUS)=BASIN(LMINUS)
                  ELSEIF (BASIN(LMINUS).EQ.0) THEN
                     BASIN(LMINUS)=BASIN(LPLUS)
                  ELSE
                     BASIN(LPLUS)=MIN(BASIN(LPLUS),BASIN(LMINUS))
                     BASIN(LMINUS)=BASIN(LPLUS)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO
      IF (.NOT.CHANGED) EXIT
   ENDDO 
!
!  At this point all minima are assigned to superbasins.
!
   IF (DEBUG) PRINT '(A,I6)','superbasin analysis done NBASIN=',NBASIN
   BASINT(1:NBASIN)=.FALSE.
   IF (DIRECTION.EQ.'AB') THEN
      DO J1=1,NGROUPA
         IF (BASIN(LOCATIONA(J1)).GT.0) BASINT(BASIN(LOCATIONA(J1)))=.TRUE.
      ENDDO 
      DO J1=1,NGROUPS
         IF (BASIN(J1).EQ.0) CYCLE
         IF (BASINT(BASIN(J1)).AND.(.NOT.GROUPA(J1))) THEN
            IF (BARRIER(J1).LT.0.0D0) BARRIER(J1)=ETHRESH-NEWEMIN(J1)
         ENDIF
      ENDDO
   ELSE
      DO J1=1,NGROUPB
         IF (BASIN(LOCATIONB(J1)).GT.0) BASINT(BASIN(LOCATIONB(J1)))=.TRUE.
      ENDDO 
      DO J1=1,NGROUPS
         IF (BASIN(J1).EQ.0) CYCLE
         IF (BASINT(BASIN(J1)).AND.(.NOT.GROUPB(J1))) THEN
            IF (BARRIER(J1).LT.0.0D0) THEN
               BARRIER(J1)=ETHRESH-NEWEMIN(J1)
            ENDIF
         ENDIF
      ENDDO
   ENDIF
   IF (DEBUG) PRINT '(A,F15.5,A,I6)','superbasin analysis performed for threshold ',ETHRESH,' # basins=',NBASIN
   ETHRESH=ETHRESH+EINC
   IF (ETHRESH.GT.HIGHESTTS+EINC) EXIT 
ENDDO
IF (DEBUG) THEN
   DO J1=1,NGROUPS
!     IF (BARRIER(J1).GT.0.0D0) PRINT '(A,I6,G20.10,I6)','getfreebarrier> J1,BARRIER=',J1,BARRIER(J1)
      PRINT '(A,I6,G20.10,I6)','getfreebarrier> J1,BARRIER=',J1,BARRIER(J1)
   ENDDO
ENDIF

END SUBROUTINE GETFREEBARRIER

