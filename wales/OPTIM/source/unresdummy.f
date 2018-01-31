C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
      SUBROUTINE UNRSETZSYMATMASS
      IMPLICIT NONE

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END
 
      SUBROUTINE UENERGY(X,GRAD,etot,GRADT,SECT)
      IMPLICIT NONE
      DOUBLE PRECISION GRAD(*), X(*)
      DOUBLE PRECISION etot
      LOGICAL GRADT,SECT

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END 
 
      SUBROUTINE UPDATEDC
      IMPLICIT NONE

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END SUBROUTINE UPDATEDC

      SUBROUTINE chainbuild
      IMPLICIT NONE

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE int_from_cart(DUM1,DUM2)
      IMPLICIT NONE
      LOGICAL DUM1,DUM2

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE geom_to_var(DUM1,DUM2)
      IMPLICIT NONE
      INTEGER DUM1
      DOUBLE PRECISION DUM2(DUM1)

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE var_to_geom(DUM1,DUM2)
      IMPLICIT NONE
      INTEGER DUM1
      DOUBLE PRECISION DUM2(DUM1)

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE UNRSETDIHE
      IMPLICIT NONE

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE UNEWREAD(X,Y,Z,NATOMS,FILTH,FILTHSTR)
      IMPLICIT NONE
      INTEGER NATOMS, FILTH
      DOUBLE PRECISION X(NATOMS),Y(NATOMS),Z(NATOMS)
      CHARACTER(LEN=80) :: FILTHSTR

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE UNRESDUMP2(COORDS,IUNIT)
      IMPLICIT NONE
      INTEGER IUNIT
      DOUBLE PRECISION COORDS(*)

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE MYUNRESDUMP(COORDS,FNAMEF)
      IMPLICIT NONE
      CHARACTER(LEN=*)  FNAMEF 
      DOUBLE PRECISION COORDS(*)

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE UNRESDUMP3(COORDS,FNAMEF)
      IMPLICIT NONE
      CHARACTER(LEN=*)  FNAMEF
      DOUBLE PRECISION COORDS(*)

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE GETSTUFF(KD,NNZ,NINTB)
      IMPLICIT NONE
      INTEGER KD,NINTB,NNZ

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE INTSECDET(COORDS,NCART,KD,NNZ,NINTB,LAMBDA)
      IMPLICIT NONE
      INTEGER NNZ,NCART,KD,NINTB
      DOUBLE PRECISION COORDS(NCART),LAMBDA(NCART)

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE UNRESCALCDIHE(DIHE,ALLANG,QLOCAL,FIN)
      IMPLICIT NONE
      DOUBLE PRECISION DIHE,ALLANG,QLOCAL(*),FIN(*)

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE UNRESCALCDIHEREF(DIHE,ALLANG,QLOCAL)
      IMPLICIT NONE
      DOUBLE PRECISION DIHE,ALLANG,QLOCAL(*)

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE UNREADREF(NATOMS)
      IMPLICIT NONE
      INTEGER NATOMS

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE UNRESCALCRGYR(RGYR,QLOCAL)
      IMPLICIT NONE
      DOUBLE PRECISION RGYR,QLOCAL(*)

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END
 
      SUBROUTINE UNRSTWISTDIHE(X,Y,Z,DMODE,DPERT)
      IMPLICIT NONE
      DOUBLE PRECISION X(*),Y(*),Z(*),DPERT
      INTEGER DMODE

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE UNRSPERTDIHE(X,Y,Z,UNPMIN,UNPMAX,UNNMIN,UNNMAX,ISEED)
      IMPLICIT NONE
      INTEGER ISEED
      DOUBLE PRECISION UNPMIN,UNPMAX,UNNMIN,UNNMAX,X(*),Y(*),Z(*)

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END
 
      SUBROUTINE UNRSTWISTALL(X,Y,Z,DMODE,DPERT)
      IMPLICIT NONE
      INTEGER DMODE
      DOUBLE PRECISION X(*),Y(*),Z(*),DPERT

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE UNRESGUESSTS(Q,ITEST,PTEST,TWISTTYPE,TWISTFRAC,GUESSFAIL,DISTPF)
      IMPLICIT NONE
      DOUBLE PRECISION TWISTFRAC,Q(*)
      DOUBLE PRECISION DISTPF
      INTEGER TWISTTYPE
      LOGICAL PTEST,ITEST,GUESSFAIL

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE UNRESGUESSMIN(Q,PTEST,TWISTTYPE,NGUESS)
      IMPLICIT NONE
      DOUBLE PRECISION Q(*)
      LOGICAL PTEST
      INTEGER TWISTTYPE,NGUESS

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE UNRESGUESSTS3(Q,ITEST,PTEST,TWISTTYPE,TWISTFRAC,GUESSFAIL,DISTPF)
      IMPLICIT NONE
      DOUBLE PRECISION TWISTFRAC,Q(*)
      DOUBLE PRECISION DISTPF
      INTEGER TWISTTYPE
      LOGICAL PTEST,ITEST,GUESSFAIL

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE INTBFGSTS(ITMAX,COORDS,ENERGY,GRAD,MFLAG,RMS,EVALMIN,EVALMAX,VECS,ITER,POTCALL,PTEST)
      IMPLICIT NONE
      INTEGER ITER,ITMAX
      DOUBLE PRECISION GRAD(*),ENERGY,COORDS(*),RMS,EVALMIN,EVALMAX,VECS(*)
      LOGICAL MFLAG,PTEST,POTCALL

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE UNRESCALCDIHESEC(DIHE,ALLANG,QLOCAL,ORDERSTOP)
      IMPLICIT NONE
      DOUBLE PRECISION DIHE,ALLANG,QLOCAL(*)
      LOGICAL ORDERSTOP

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE UNRESGUESSTSSEC(Q,ITEST,PTEST,TWISTTYPE,TWISTFRAC,GUESSFAIL,DISTPF)
      IMPLICIT NONE
      DOUBLE PRECISION TWISTFRAC,Q(*),DISTPF
      INTEGER TWISTTYPE
      LOGICAL PTEST,ITEST,GUESSFAIL

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END 

      SUBROUTINE UNRESINIT
      IMPLICIT NONE

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE UNRESDIHENEB(Q,FINISH,POINTS)
      IMPLICIT NONE
      DOUBLE PRECISION Q(*), FINISH(*), POINTS(*)

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP


      RETURN
      END

      SUBROUTINE UNRESGETDIHEDIST(DIHEDIST,QA,QB,QC)
      IMPLICIT NONE

      DOUBLE PRECISION QA(*),QB(*),QC(*)
      DOUBLE PRECISION DIHEDIST

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END

      SUBROUTINE UNMECCANO(ITEST,PTEST,EMAX,PERM,Q,FIN,E1,E2,RMSINITIAL,RMSFINAL)
      IMPLICIT NONE
      DOUBLE PRECISION :: EMAX, Q(*), FIN(*), E1, E2, RMSINITIAL, RMSFINAL
      LOGICAL :: ITEST, PTEST, PERM

      PRINT*,'ERROR - this dummy routine should never be called - wrong executable?'
      STOP

      RETURN
      END
