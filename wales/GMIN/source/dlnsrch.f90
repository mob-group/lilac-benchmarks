!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
SUBROUTINE LNSRCH(N, XOLD, FOLD, G, P, X, F, STPMAX, CHECK, PROBLEM)
! This subroutine doesn't do anything at the moment, apart from stop the program
! because the user has requested the BFGS potential without implementing a 
! line search subroutine.
    USE PREC
    IMPLICIT NONE
! Arguments
    INTEGER(INT32)  :: N
    REAL(REAL64)    :: XOLD(N), FOLD, G(N), P(N), X(N), F, STPMAX
    LOGICAL         :: CHECK, PROBLEM

    WRITE(*, '(A)') 'Please supply the lnsrch routine!'
    STOP

END SUBROUTINE LNSRCH
