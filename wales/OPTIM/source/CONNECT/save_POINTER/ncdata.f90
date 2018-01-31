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
MODULE CONNECTDATA
     IMPLICIT NONE
     SAVE

     INTEGER :: NDIJPAIRS,listlength=0
     INTEGER, ALLOCATABLE :: DIJPAIR(:,:)
     DOUBLE PRECISION, ALLOCATABLE :: DIJPAIRDIST(:)
     LOGICAL :: FINISHED
!--------------------------------------------------------

     INTEGER :: NMIN=0,NTS=0 !NUMBER OF TS/MIN FOUND SO FAR
     INTEGER :: NTSOLD=0
     INTEGER :: NATOMS,NOPT,NCONDONE,DEPTH
     INTEGER :: TSRACKSIZE=100,MINRACKSIZE=150
     DOUBLE PRECISION :: RMS
     DOUBLE PRECISION,ALLOCATABLE :: G(:)
     DOUBLE PRECISION :: D

     LOGICAL :: MOREPRINTING,ENDREACHED
     LOGICAL :: RANDOM

     CHARACTER(LEN=7)  :: CHR,CHR2
     CHARACTER(LEN=20) :: STR1,STR2,STR3,STR4,STR5,STR6

!--------------------------------------------------------
     TYPE MINDATA
          DOUBLE PRECISION,POINTER :: E
          DOUBLE PRECISION,POINTER :: X(:)
          DOUBLE PRECISION,POINTER :: D(:)      ! HOLDS DISTANCES TO ALL MINIMA KNOWN AT THE MOMENT OF ADDITION OF THIS MINIMUM
          DOUBLE PRECISION,POINTER :: INTERP(:) ! HOLDS DISTANCES TO ALL MINIMA KNOWN AT THE MOMENT OF ADDITION OF THIS MINIMUM
          INTEGER,POINTER :: NTRIES(:) ! HOLDS NUMBER OF TRIES THIS MINIMA WAS TRIED TO BE CONNECTED WITH ANY OTHER MINIMA -//-;
                                       ! this affects number of images used in next try (if any)
          INTEGER,POINTER :: C(:)      ! TOTAL NUMBER OF CONNECTIONS TO THIS MINIMUM (SAME AS CURRENT SIZE OF CTS AND CMIN!)
          INTEGER,POINTER :: CTS(:)
          INTEGER,POINTER :: CMIN(:)
          LOGICAL         :: S,F
     END TYPE MINDATA

     TYPE TSDATA
          DOUBLE PRECISION,POINTER :: E
          DOUBLE PRECISION,POINTER :: X(:)
          DOUBLE PRECISION,POINTER :: EVALMIN
          DOUBLE PRECISION,POINTER :: VECS(:)
          INTEGER         :: P,M !NUMBERS OF MINIMA PLUS AND MINUS WHICH TS IS CONNECTING
          DOUBLE PRECISION         :: SLENGTH,SLP,DISP,GAMMA,NTILDE ! SLP - SLENGTHPLUS
          LOGICAL         :: BAD ! IS SET TO .TRUE. IF PATH RUN DID NOT CONVERGE
     END TYPE TSDATA

     TYPE MINTYPE
          TYPE(MINDATA) :: DATA
     END TYPE MINTYPE

     TYPE TSTYPE
          TYPE(TSDATA) :: DATA
     END TYPE TSTYPE

     TYPE(MINTYPE),POINTER :: MI(:),TEMPMINRACK(:)
     TYPE(TSTYPE),POINTER  ::  TS(:),TEMPTSRACK(:)
!--------------------------------------------------------
     TYPE CHAIN
          INTEGER :: I,J ! MINIMA I IS CONNECTED TO NEXT%I VIA TS J; J IS INDEX OF TS IN ARRAY CTS OF MINIMA I
          TYPE(CHAIN),POINTER :: NEXT
     END TYPE CHAIN

     TYPE(CHAIN),POINTER :: START,NEW,DUMMY
!--------------------------------------------------------
     TYPE CHAIN2
          INTEGER :: I
          TYPE(CHAIN2),POINTER :: NEXT
     END TYPE CHAIN2

     TYPE(CHAIN2),POINTER :: START2,NEW2,DUMMY2,ONEUP
!--------------------------------------------------------
     DOUBLE PRECISION :: EV
     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: FRQSTS, EVEC, FRQSPLUS, FRQSMINUS, GDUMMY
END MODULE CONNECTDATA
