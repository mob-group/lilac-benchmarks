!   NEB module is an implementation of the nudged elastic band method for performing double-ended pathway searches.
!   Copyright (C) 2003-2006 Semen A. Trygubenko and David J. Wales
!   This file is part of NEB module. NEB module is part of OPTIM.
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
MODULE NEBDATA
     IMPLICIT NONE
     SAVE
     INTEGER :: NATOMS,NOPT,NINTS=0
     INTEGER :: NITERDONE,NITERDONESAVE=0,EXITSTATUS
     DOUBLE PRECISION :: STARTTIME,ENDTIME,AVDEV=0.0D0,AVDEVOLD,ETOTAL,SEPARATION,RMS,STEPTOT
     DOUBLE PRECISION, ALLOCATABLE :: TRUEGRAD(:)
     DOUBLE PRECISION,POINTER :: DVEC(:),DEVIATION(:),TANVEC(:,:),STEPIMAGE(:),X(:),EIMAGE(:),G(:),GSPR(:),RMSFIMAGE(:),&
 &                               XYZ(:),EEE(:),GGG(:),SSS(:),RRR(:), XCART(:), XYZCART(:), GCART(:), GGGCART(:)
     DOUBLE PRECISION, ALLOCATABLE :: NEWNEBK(:)
     LOGICAL, ALLOCATABLE :: BADIMAGE(:), BADPEPTIDE(:)
     LOGICAL :: BADTAU=.FALSE.
     ! efk: for freezenodes
     LOGICAL, ALLOCATABLE :: IMGFREEZE(:)
END MODULE NEBDATA
