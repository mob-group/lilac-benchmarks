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
   
SUBROUTINE AADISTANCE( V1, V2 )
!js850>
!minimize the "distance" between two angle axis vectors
!
!perform symmetry operations on angle axis vectors v1 and v2 such that 
!the cartesian distance norm(v1(:) - v2(:)) between them is minimized.
! 
!the possible symmetry operations are
! v1 -> v1 * ( 1 + n * 2*pi / norm(v1) ) where n is an integer
   IMPLICIT NONE
   
   DOUBLE PRECISION, INTENT(INOUT) :: V1(3), V2(3)
   DOUBLE PRECISION D, V1N, V2N
   DOUBLE PRECISION, PARAMETER :: PID = 3.14159265358973D0
   DOUBLE PRECISION, PARAMETER :: PID2 = 2.D0*PID
   DOUBLE PRECISION DOLD, DNEW, DNEWNEW
   
   
   
!dold = sqrt(sum( (v2 - v1)**2 )) !for debugging
   
!make the magnitude of v1 as close to zero as possible
   V1N = SQRT( SUM(V1(:)**2) )
   DO WHILE ( V1N .GE. PID )
      V1(:) = V1(:)/V1N*(V1N-PID2)
      V1N = SQRT( SUM(V1(:)**2) )
   ENDDO
   
!make the magnitude of v2 as close to zero as possible
   V2N = SQRT( SUM(V2(:)**2) )
   DO WHILE ( V2N .GE. PID )
      V2(:) = V2(:)/V2N*(V2N-PID2)
      V2N = SQRT( SUM(V2(:)**2) )
   ENDDO
   
! js850>
! We still need to check the two possibilities
! v1 -> v1 * (v1n - 2*pi)/v1n   
! v2 -> v2
!
! and
!
! v1 -> v1
! v2 -> v2 * (v2n - 2*pi)/v2n
!
! Other symmetry operations can be ruled out as options to give to closer vectors.
! We can simplify the calculation by noticing that the smaller of v1 and v2
! will always remain unchanged.  Finally, we can check if the larger vector
! needs to change by checking the condition (asuming v1 is the larger)
!
!   norm(v1) - dot_product(v1, v2) / norm(v1)  > Pi
   
!dnew = sqrt(sum( (v2 - v1)**2 )) ! for debugging
   
   V1N = SQRT(SUM( V1**2 ))
   V2N = SQRT(SUM( V2**2 ))
   IF ( V1N > V2N ) THEN
      D = V1N - DOT_PRODUCT(V1, V2) / V1N
      IF (D > PID) THEN
         V1 = V1 * (1.D0 - PID2/V1N)
      ENDIF
   ELSE
      D = V2N - DOT_PRODUCT(V1, V2) / V2N
      IF (D > PID) THEN
         V2 = V2 * (1.D0 - PID2/V2N)
      ENDIF
   ENDIF
!dnewnew = sqrt(sum( (v2 - v1)**2 )) !js850>
!if (dold > 2) write(*,*) "js850>", dold, dnew, dnewnew !js850>
!if (dnewnew > 3) write(*,*) "js850> something may have gone wrotg"
   RETURN
END SUBROUTINE AADISTANCE
