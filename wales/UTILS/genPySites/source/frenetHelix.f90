      SUBROUTINE FRENETHELIX(R,C,TH,S,PHI,T,N,B,D)

! Computes the Frenet frame of a helix with parameters: 
! r, c, t, s and phi. where r is the radius of the 
! helix, c=L/2*pi is the pitch length per radian, t is the conventional helix 
! parameter in radians, which denotes the position on the helix, s is the 
! length of arc (along the helix itself) by which the helix is displaced in a
! screw like motion about its central axis (ie small rotation and longitudinal 
! displacement) and phi is the phase of the helix about its central axis.
!The helix is aligned with the z-axis of the lab frame. The x-axis is defined by the point on the helix where t=0, s=0 and phi=0.
!z=0 is also defined by this point. 
     IMPLICIT NONE

      DOUBLE PRECISION :: D(1:3), T(1:3), N(1:3), B(1:3)
      DOUBLE PRECISION :: R,C,TH,S,PHI
      DOUBLE PRECISION :: RSQ, CSQ, TZERO, HARG, HSIZE

      RSQ=R*R
      CSQ=C*C
      HSIZE=SQRT(RSQ+CSQ)
      TZERO=S/HSIZE

!     HARG = t + t0 + phi = t + s/sqrt(r^2+c^2) + (phi)  
      HARG=TH+TZERO+PHI

      T(1)=-1.0*R*sin(HARG)/HSIZE
      T(2)=R*cos(HARG)/HSIZE
      T(3)=C/HSIZE
      N(1)=-1.0*cos(HARG)
      N(2)=-1.0*sin(HARG)
      N(3)=0.0
      B(1)=C*sin(HARG)/HSIZE
      B(2)=-1.0*C*cos(HARG)/HSIZE
      B(3)=R/HSIZE

      D(1)=R*cos(HARG)
      D(2)=R*sin(HARG)
      D(3)=C*(TH+TZERO)

      return
      END SUBROUTINE FRENETHELIX
