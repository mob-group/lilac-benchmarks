!   CONNECT module is an implementation of a connection algorithm for finding rearrangement pathways.
!   Copyright (C) 2003-2006 Semen A. Trygubenko and David J. Wales
!   This file is part of CONNECT module. CONNECT module is part of OPTIM.
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

! DJW stationary point dump in PATHSAMPLE format, although with dummy entries
! for frequencies, point group orders and moments of inertia.

     SUBROUTINE DODUMPSP
     USE CONNECTDATA
     IMPLICIT NONE
     INTEGER NDUMMY, J1

! dump minimum energies in min.A.new, min.B.new and min.data.new

     OPEN(24,FILE='min.A.new',STATUS='UNKNOWN')
     WRITE(24,'(2F20.10,I5,4F20.10)') mi(1)%data%E,0.0D0,1,1.0D0,1.0D0,1.0, 1.0D0
     CLOSE(24)
     OPEN(24,FILE='min.B.new',STATUS='UNKNOWN')
     WRITE(24,'(2F20.10,I5,4F20.10)') mi(2)%data%E,0.0D0,1,1.0D0,1.0D0,1.0, 1.0D0
     CLOSE(24)
     OPEN(24,FILE='min.data.new',STATUS='UNKNOWN')
     DO J1=3,NMIN
        WRITE(24,'(2F20.10,I5,4F20.10)') mi(J1)%data%E,0.0D0,1,1.0D0,1.0D0,1.0, 1.0D0
     ENDDO
     CLOSE(24)

! dump ts energies in ts.data.new

     OPEN(24,FILE='ts.data.new',STATUS='UNKNOWN')
     DO J1=1,NTS
        WRITE(24,'(2F20.10,3I6,3F20.10)') TS(J1)%data%E,1.0D0,1,TS(J1)%DATA%P,TS(J1)%DATA%M,1.0D0,1.0D0,1.0D0
     ENDDO
     CLOSE(24)

! dump minimum coordinates in points.min.new

     INQUIRE(IOLENGTH=NDUMMY) MI(1)%DATA%X(1:NOPT)
     OPEN(13,FILE='points.min.new',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=NDUMMY)
     DO J1=1,NMIN
        WRITE(13,REC=J1) MI(J1)%DATA%X(1:NOPT)
     ENDDO
     CLOSE(13)

     OPEN(13,FILE='points.ts.new',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=NDUMMY)
     DO J1=1,NTS
        WRITE(13,REC=J1) TS(J1)%DATA%X(1:NOPT)
     ENDDO
     CLOSE(13)
 
     END SUBROUTINE DODUMPSP

