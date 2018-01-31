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
!  Subroutine to provide candidate pairs of minima for untraping.
!  A metric is used to find the closest match in the product basin.
!
SUBROUTINE GETUPAIRMETRIC(UPMIN1,MINMET, METRICUPAIR, METRICM, METMATMAX)
IMPLICIT NONE
INTEGER METRICUPAIR, UPMIN1, UPMIN2, I, NUM, NUMA, METMATMAX, MATMID 
DOUBLE PRECISION ANGLE, ANGLE2, ANGLE2A, FCC, FCC2, HCP, HCP2
INTEGER MINMET(1:METMATMAX)
LOGICAL METRICM, TWOONLY
LOGICAL :: DEBUG=.TRUE.  

TWOONLY=.FALSE.
METRICM=.FALSE.
MATMID=INT(METMATMAX/2.0D0)
! IF (DEBUG) PRINT*, METRICUPAIR
OPEN(UNIT=21,FILE='metric',ACTION='READ',STATUS='OLD')
IF (METRICUPAIR.EQ.1) THEN
  DO I=1, UPMIN1
   ANGLE=-HUGE(1.0D0)
   READ(21,*, END=22) ANGLE
  ENDDO
  22 CLOSE(21) 
  IF (DEBUG) WRITE(6,'(A33, F11.5)') 'Metricupair=1, metric value is ', ANGLE 
ELSE IF (METRICUPAIR.EQ.2) THEN
  DO I=1, UPMIN1
   ANGLE=-HUGE(1.0D0)
   READ(21,*, END=32) ANGLE, NUM, FCC
  ENDDO
  32 CLOSE(21) 
  IF (DEBUG) WRITE(6,'(A33, 2F11.5)') 'Metricupair=2, metric values are ', ANGLE, FCC 
ELSE IF (METRICUPAIR.EQ.3) THEN
  DO I=1, UPMIN1
   ANGLE=-HUGE(1.0D0)
   READ(21,*, END=33) ANGLE, NUM, FCC, HCP
  ENDDO
  33 CLOSE(21) 
  IF (DEBUG) WRITE(6,'(A33, 3F11.5)') 'Metricupair=3, metric values are ', ANGLE, FCC, HCP 
ELSE
  PRINT*, 'METRICUPAIR must use 1, 2 or 3 metrics - using default, 1 metric'
  METRICUPAIR=1
  DO I=1, UPMIN1
   ANGLE=-HUGE(1.0D0)
   READ(21,*, END=23) ANGLE
  ENDDO
  23 CLOSE(21)
  IF (DEBUG) WRITE(6,'(A33, F11.5)') 'Metricupair=1, metric value is ', ANGLE
ENDIF  
IF (ANGLE.EQ.-HUGE(1.0D0)) THEN
  METRICM=.FALSE.
  RETURN
ENDIF  
OPEN(UNIT=23,FILE='metric.ordered',ACTION='READ',STATUS='OLD')
ANGLE2=-HUGE(1.0D0)
MINMET(:)=-1
IF (METRICUPAIR.EQ.1) THEN
  DO 
   ANGLE2=-HUGE(1.0D0)
   READ(23,*, END=24) ANGLE2, NUM
   MINMET(1:MATMID-1)=MINMET(2:MATMID)
   MINMET(MATMID)=NUM
   IF (ANGLE2.GT.ANGLE) THEN
     DO I=MATMID+1, METMATMAX
      READ(23,*, END=24) ANGLE, NUM
      MINMET(I)=NUM 
     ENDDO
   EXIT
   ENDIF
  ENDDO
ELSE IF (METRICUPAIR.EQ.2) THEN
  FCC2=-HUGE(1.0D0)
  DO 
   ANGLE2=-HUGE(1.0D0)
   FCC2=-HUGE(1.0D0)
   READ(23,*, END=24) ANGLE2, NUM, FCC2
   MINMET(1:MATMID-1)=MINMET(2:MATMID)
   MINMET(MATMID)=NUM
   IF (ANGLE2.GE.ANGLE.AND.FCC2.GE.FCC) THEN
     DO I=MATMID+1, METMATMAX
      READ(23,*, END=24) ANGLE, NUM, FCC
      MINMET(I)=NUM 
     ENDDO
     EXIT
   ENDIF
  ENDDO
ELSE IF (METRICUPAIR.EQ.3) THEN
  FCC2=-HUGE(1.0D0)
  HCP=-HUGE(1.0D0)
  DO 
   ANGLE2=-HUGE(1.0D0)
   FCC2=-HUGE(1.0D0)
   HCP=-HUGE(1.0D0)
   READ(23,*, END=24) ANGLE2, NUM, FCC2, HCP2
   MINMET(1:MATMID-1)=MINMET(2:MATMID)
   MINMET(MATMID)=NUM
   IF (ANGLE2.GE.ANGLE.AND.FCC2.GE.FCC.AND.HCP2.GE.HCP) THEN
     DO I=MATMID+1, METMATMAX
      READ(23,*, END=24) ANGLE, NUM, FCC, HCP
      MINMET(I)=NUM 
     ENDDO
     EXIT
   ENDIF
  ENDDO
ENDIF
24 IF (ANGLE2.GT.-HUGE(1.0D0)) METRICM=.TRUE.
METRICM=.TRUE.
! PRINT*, MINMET(MATMID)
CLOSE(23)
END SUBROUTINE GETUPAIRMETRIC
