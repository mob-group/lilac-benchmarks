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
!
!*************************************************************************
!
!  Subroutine MORSE_bulk calculates the energy and gradient analytically for 
!  the Morse potential with periodic boundary conditions. The potential has
!  a cutoff and is shifted to make the potential continuous.  It is not smooth
!  though.
!
!  subroutine morse_bulk_wrapper makes calling it from potentials.f a bit
!  simpler.
!
!  use this potential by assigning the atoms the label "M" as with normal morse
!  potential, but additionally pass the keyword BULK
!
!  the options for the potential are passed using the keyword 
!
!  PARAMS rho boxlx boxly boxlz rcut
!  
!  where rho is the inverse width of the well.  boxlx, boxly, boxlz are the box
!  lengths.  And rcut is the cutoff.
!
!*************************************************************************
!
subroutine morse_bulk_wrapper(x, v, emorse, gtest, stest)
      use commons, only : natoms, param1, param2, param3, param4, param5, debug
      IMPLICIT NONE 
      LOGICAL, intent(IN) :: GTEST, stest
      DOUBLE PRECISION, intent(IN) :: X(3*NATOMS)
      DOUBLE PRECISION, intent(OUT) :: V(3*NATOMS), EMORSE
      logical periodic, use_cutoff
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION DIST, R, DUMMY, &
                       RR(NATOMS,NATOMS), &
                       XMUL2, iboxvec(3), dx(3), eshift
      DOUBLE PRECISION rho, R0, A, boxvec(3), rcut

      rho = param1
      boxvec(1) = param2
      boxvec(2) = param3
      boxvec(3) = param4
      rcut = param5

      R0 = 1.d0
      A = 1.d0
      !R0 = 2.8970
      !A = 0.7102
      !rho = 1.6047
      periodic = .true.
      use_cutoff = .true.

      if (rcut / R0 < 1.d-1) then
         write(*,*) "morse_bulk> warning the cutoff is very small", rcut
      endif

      if (stest) then
         write(*,*) "second derivatives are not yet implemented for morse bulk"
         write(*,*) "exiting"
!        call exit(1)
         STOP
      endif
      !if (debug) then
         !write(*,*) "rho, R0, A, rcut, boxl", rho, r0, a, rcut, boxvec(1)
      !endif

      call morse_bulk(X,V,EMORSE,GTEST, natoms, rho, R0, A, periodic, & 
         boxvec, use_cutoff, rcut)


end subroutine morse_bulk_wrapper

      SUBROUTINE MORSE_BULK(X,V,EMORSE,GTEST, natoms, rho, R0, A, periodic, &
         boxvec, use_cutoff, rcut)
      ! R0 is the position of the bottom of the well
      ! rho is the width of the well and has units of inverse length
      ! A is the energy scale
!      USE commons
      IMPLICIT NONE 
      LOGICAL, intent(IN) :: GTEST, periodic, use_cutoff
      integer, intent(IN) :: NATOMS
      DOUBLE PRECISION, intent(IN) :: X(3*NATOMS), rho, R0, A, boxvec(3), rcut
      DOUBLE PRECISION, intent(OUT) :: V(3*NATOMS), EMORSE
      INTEGER J1, J2, J3, J4
      DOUBLE PRECISION DIST, R, DUMMY, &
                       RR(NATOMS,NATOMS), &
                       XMUL2, iboxvec(3), dx(3), eshift
!     LOGICAL EVAP, evapreject
!     COMMON /EV/ EVAP, evapreject
      if (periodic) iboxvec(:) = 1.d0 / boxvec(:)

      if (use_cutoff) then
         Eshift = (1.d0 - exp(rho * (r0 - rcut)))**2 - 1.d0
         !write(*,*) "Eshift", eshift, rcut
      endif

!     EVAP=.FALSE.
      V(:) = 0.D0
      EMORSE=0.0D0
      DO J1=1,NATOMS
         J3=3*J1
         RR(J1,J1)=0.0D0
         DO J2=J1+1,NATOMS
            J4=3*J2
            dx(:) = X(J3-2:j3)-X(J4-2:j4)
            if (periodic) then
               dx = dx - boxvec * nint(dx * iboxvec)
            endif
            dist = max(sqrt(sum(dx**2)), 1.d-5)

            if (use_cutoff .and. dist.ge.rcut) cycle

            R=DEXP(RHO*R0-RHO*DIST)
            DUMMY=R*(R-2.0D0)
            EMORSE=EMORSE+DUMMY - Eshift

            if (gtest) then
               xmul2 = 2.0D0*R*(R-1.0D0)/DIST * A
               V(J3-2:j3) = V(j3-2:j3) - xmul2 * dx
               V(J4-2:j4) = V(j4-2:j4) + xmul2 * dx
            endif
         ENDDO
      ENDDO
      EMORSE = EMORSE * A

      RETURN
      END
