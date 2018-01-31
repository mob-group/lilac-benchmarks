!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of OPTIM.
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
MODULE modcharmm
   IMPLICIT NONE
   SAVE

! DAE old CHARMM-related common blocks

      INTEGER, DIMENSION(:), ALLOCATABLE  :: NPHIPSI,NOMEGAC,NSIDECHAIN,NCHIRAL  ! CHDIHE
      INTEGER :: NPHIPSITOT,NOMEGACTOT,NSIDECHAINTOT,NCHIRALTOT
      INTEGER, DIMENSION(:), ALLOCATABLE :: PHIPSI,OMEGAC,SIDECHAIN,PSFPHI  ! CHDIHE
      INTEGER, ALLOCATABLE, DIMENSION(:) :: NSEGATOMS     ! dimension will be NSEG
      INTEGER :: CHARMMTYPE, CHSEEDI, CHSEEDJ, CHSEEDK    ! CHARMMWORD
      LOGICAL, DIMENSION(:), ALLOCATABLE :: ISSIDECHAIN
      LOGICAL :: TOMEGAC,TSIDECHAIN     ! CHARMMTWIST
      INTEGER, DIMENSION(:), ALLOCATABLE :: CHIRAL(:)
      LOGICAL :: CALCDIHE,OSASAT,ODIHET,CHECKCHIRALT,CHECKOMEGAT, ICINTERPT, CHBIT
      DOUBLE PRECISION :: MINOMEGA

      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: REFCOORD, REFPPSANGLE   ! CHREF
      INTEGER, DIMENSION(:), ALLOCATABLE :: REFIHB, REFJHB, REFKHB
      INTEGER  :: REFNHB

      INTEGER :: NMINB,NORDERB
      DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: EMINB, ITXMINB, ITYMINB, ITZMINB, OMEGASTART
      DOUBLE PRECISION, DIMENSION(:) :: ORDERBPARAM(2)
      LOGICAL :: ORDERBT
      CHARACTER(LEN=8) ORDERBLABEL(2)

      LOGICAL :: FAILT, REBUILDSCT   ! used to be passed to NEB
      INTEGER :: TWISTTYPE, NGUESS, NCHENCALLS, ACEUPSTEP
      LOGICAL :: TRYNEB, ACESOLV, CHDEBUG, EDEBUG
      DOUBLE PRECISION  :: GUESSTHRESH, REBUILDSC
      DOUBLE PRECISION :: RPRO,PTRANS,TRANSMAX,PROT,ROTMAX
      LOGICAL :: CHARMMNOTUPDATE

! formerly in key.h

      INTEGER :: NPERMDIHE
      LOGICAL :: CHRMMT,TWISTDIHET,PERTDIHET,GUESSTST,NOCISTRANS,NORANDOM,ENDHESS,ENDNUMHESS,PERMDIHET,INTMINT
      DOUBLE PRECISION :: MAXBFGSB,MAXBFGST,MAXBFGSP,MAXBFGSI,IMINCUT,RANDOMCUTOFF,PERMDIHE(20), EXTRASTEPS

END MODULE
