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
!  Add energy and gradient correction terms for compression
!
SUBROUTINE COMPRESS(X,V,ENERGY,GTEST)
    USE COMMONS
    USE GENRIGID, ONLY: NRIGIDBODY, RIGIDINIT
    USE OPP_MOD, ONLY: OPP_MODE
    IMPLICIT NONE

    LOGICAL GTEST
    INTEGER J1, J3, NMOL
    DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), ENERGY, XMASS, YMASS, ZMASS

    IF (PERIODIC) RETURN

! jwrm2> for rigid body angle axis system, only the first half of X is the body position
!        Note, as written this is a harmonic potential acting on the centre of the rigid
!        body, NOT on each site of the rigid body 
    NMOL = NATOMS
    IF (RIGID) NMOL = NATOMS/2
    IF (RIGIDINIT) NMOL = NRIGIDBODY

! jwrm2> Need to remove alchemical and unit cell degrees of freedom.
    IF (OPPT) THEN
        SELECT CASE (OPP_MODE)
        CASE(1)
            NMOL = NATOMS - 2
        CASE(2)
            NMOL = NATOMS - 3
        CASE(3)
            NMOL = NATOMS - 1
        END SELECT
    END IF

! Find centre of mass
    XMASS=0.0D0
    YMASS=0.0D0
    ZMASS=0.0D0
    DO J1=1,NMOL
        XMASS=XMASS+X(3*(J1-1)+1)
        YMASS=YMASS+X(3*(J1-1)+2)
        ZMASS=ZMASS+X(3*(J1-1)+3)
    ENDDO
    XMASS=XMASS/NMOL
    YMASS=YMASS/NMOL
    ZMASS=ZMASS/NMOL

    IF (DEBUG) WRITE (MYUNIT,*) 'compress> Compressing'

! Loop over all bodies
    DO J1=1, NMOL
        J3=3*J1
        DIST=(X(J3-2)-XMASS)**2+(X(J3-1)-YMASS)**2+(X(J3)-ZMASS)**2

! Add energy
        ENERGY=ENERGY+K_COMP*DIST/2.0D0

        IF (GTEST) THEN
! Add gradients
            V(J3-2)=V(J3-2)+K_COMP*(X(J3-2)-XMASS)
            V(J3-1)=V(J3-1)+K_COMP*(X(J3-1)-YMASS)
            V(J3)=  V(J3)  +K_COMP*(X(J3)  -ZMASS)
        ENDIF
    ENDDO

END
