      SUBROUTINE FLATHELIX(R,C,T,S,PHI,FMAT,D)

! Computes the position of the points along a helix defined by r, c, t, s and 
! phi where r is the radius of the helix, c=L/2*pi is the pitch length per radian, 
! t is the conventional helix parameter in radians, which denotes the positions to 
! be computed on the helix, s is the length of arc (along the helix itself) by 
! which the helix is displaced in a screw like motion about its central axis 
! and phi is the phase of the helix about its central axis.
! rotation matrix is returned in which the radial vector is returned as the first row vector
! and the vector perpendicular to this in the xy plane as the second row vector.
! the z-axis is returned as the third row vector.

! this rotation matrix defines the rotation about the z-axis at each point on the helix.


      IMPLICIT NONE

      DOUBLE PRECISION :: D(1:3), FMAT(1:3,1:3)
      DOUBLE PRECISION :: R, C, T, S, PHI
      DOUBLE PRECISION :: RSQ, CSQ, TZERO, HARG, HSIZE

      RSQ=R*R
      CSQ=C*C
      HSIZE=SQRT(RSQ+CSQ)
      TZERO=S/HSIZE

!     HARG = t + t0 + phi = t + s/sqrt(r^2+c^2) + (phi)  
      HARG=T+TZERO+PHI

      FMAT(1,1)=cos(HARG)
      FMAT(1,2)=sin(HARG)
      FMAT(1,3)=0
      FMAT(2,1)=-sin(HARG)
      FMAT(2,2)=cos(HARG)
      FMAT(2,3)=0.0
      FMAT(3,1)=0
      FMAT(3,2)=0
      FMAT(3,3)=1

      D(1)=R*cos(HARG)
      D(2)=R*sin(HARG)
      D(3)=C*(T+TZERO)

      !write(6,*) 'S:',S,' HSIZE: ',HSIZE,' T:',T,' TZERO: ', TZERO,' Phi:', Phi,' HARG:',HARG,' R:', R,' C:', C,' D:',D

      return
      END SUBROUTINE FLATHELIX
