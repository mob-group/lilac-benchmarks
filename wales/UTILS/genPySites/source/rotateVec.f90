      SUBROUTINE ROTATEVEC(v,u,g,vNew)
!Rotates the vector v by g radians about the axis u. +ve rotation in anticlockwise direction
      IMPLICIT NONE

      DOUBLE PRECISION :: v(1:3),u(1:3),g,vNew(1:3)
      DOUBLE PRECISION :: uCrossv(1:3)

      uCrossv(1) = u(2)*v(3) - u(3)*v(2)
      uCrossv(2) = u(3)*v(1) - u(1)*v(3)
      uCrossv(3) = u(1)*v(2) - u(2)*v(1)

      vNew=v*cos(g)-uCrossV*sin(g)+u*DOT_PRODUCT(u,v)*(1-cos(g))

      return
      END SUBROUTINE ROTATEVEC
