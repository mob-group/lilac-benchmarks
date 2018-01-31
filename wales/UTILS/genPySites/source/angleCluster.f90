      SUBROUTINE ANGLECLUSTER(a,b,g,d,lPOS,lRMAT)
      !Computes the position (d) and orientation (a,b,g) of an ellipsoid in the cluster frame

     IMPLICIT NONE

      DOUBLE PRECISION :: a,b,g,d,lPOS(1:3),lRMAT(1:3,1:3)
      DOUBLE PRECISION :: av(1:3),bv(1:3),cv(1:3)
      !set the position of the ellipsoid in the cluster frame at a point d along x-axis
      lPOS(1)=d
      lPOS(2)=DBLE(0.0)
      lPOS(3)=DBLE(0.0)

      CALL EULER2RMAT(a,b,g,lRMAT)

      !compute the unit vector for the a body axis (the cluster x-axis rotated anticlockwise about z an amount b)
      !av(1)=cos(b)
      !av(2)=sin(b)
      !av(3)=0

      !compute the unit vector for the b body axis (the y-axis rotated about the x'-axis an amount a and then the cluster frame z axis an amount b (spherical polars!))
      !bv(1)=-cos(a)*sin(b)
      !bv(2)=cos(a)*cos(b)
      !bv(3)=sin(a)

      !compute the unit vector for the c body axis (the z axis rotated about the x axis an amount a and then the cluster frame z axis an amount b)
      !cv(1)=sin(a)*sin(b)
      !cv(2)=-sin(a)*cos(b)
      !cv(3)=cos(a)

      !rotate the a and c body axes about the b-body axis by angle g. Write results as row vectors into the rotation matrix directly.
      !CALL ROTATEVEC(av(1:3),bv(1:3),g,lRMAT(1,1:3))
      !CALL ROTATEVEC(cv(1:3),bv(1:3),g,lRMAT(3,1:3))

      !copy the b-axis into the rotation matrix
      !lRMAT(2,1:3)=bv

      return
      END SUBROUTINE ANGLECLUSTER
