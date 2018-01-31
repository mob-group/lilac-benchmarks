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
MODULE KEYCONNECT
     USE KEY,ONLY: DEBUG
     USE KEYUTILS
     USE KEYDECIDE
     IMPLICIT NONE
     SAVE
     LOGICAL :: FCD             = .FALSE. 
     INTEGER :: NCONMAX         = 1
     INTEGER :: NTRIESMAX       = 1
     INTEGER :: IMAGEMAX        = 2
     DOUBLE PRECISION :: IMAGEINCR    = 0.5D0
     DOUBLE PRECISION :: IMAGEDENSITY = 10.0D0
     DOUBLE PRECISION :: ITERDENSITY  = 30.0D0
     !logical :: debug          = .False.
     CONTAINS

     SUBROUTINE KEYCONNECTPRINT
          USE CHARUTILS
          USE KEYNEB, ONLY: NIMAGE,NITERMAX,ALLKEYNEBPRINT
          USE GSDATA, ONLY: KEYGSPRINT, GSITERDENSITY
          USE KEY, ONLY : GROWSTRINGT
          IMPLICIT NONE
          
          INTSTR=WI(NCONMAX)
          INTSTR2=WI(IMAGEMAX)
          WRITE(*,'(a)') ' KeyConnect> Maximum cycles = '//trim(IntStr)//', maximum images = '//trim(IntStr2)
          INTSTR=WI(NTRIESMAX)
          REALSTR=WR(IMAGEINCR,2)
          WRITE(*,'(a)') ' KeyConnect> Maximum attempts per pair of minima = '//trim(IntStr)//&
                        &', with increment image density of '//trim(RealStr)
          REALSTR=WR(IMAGEDENSITY,2)
          REALSTR2=WR(ITERDENSITY,2)
          WRITE(*,'(a)') ' KeyConnect> Image density = '//trim(RealStr)//', iteration density = '//trim(RealStr2)
          IF (GROWSTRINGT) THEN
             CALL KEYGSPRINT(.TRUE.)
          ELSE
             CALL ALLKEYNEBPRINT(.TRUE.)
          ENDIF
          IF (FCD) THEN ! FIRST CYCLE DIFFERENT
               INTSTR=WI(NIMAGE)
               IF (GROWSTRINGT) THEN
                  INTSTR2=WI(INT(GSITERDENSITY*NIMAGE))
                  WRITE(*,'(a)')&
                       &' KeyConnect> Using '//trim(IntStr)//' images and '//trim(IntStr2)//' iterations in the first string run'
               ELSE
                  INTSTR2=WI(NITERMAX)
                  WRITE(*,'(a)')&
                       &' KeyConnect> Using '//trim(IntStr)//' images and '//trim(IntStr2)//' iterations in the first NEB run'
               ENDIF
          ENDIF
          IF (DEBUG) THEN
               WRITE(*,'(1x,a)') "KeyConnect> Verbose printing is on"
          ENDIF
     END SUBROUTINE KEYCONNECTPRINT

     SUBROUTINE ALLKEYCONNECTPRINT
          IMPLICIT NONE
          CALL KEYCONNECTPRINT
          CALL KEYDECIDEPRINT
     END SUBROUTINE ALLKEYCONNECTPRINT
END MODULE KEYCONNECT
