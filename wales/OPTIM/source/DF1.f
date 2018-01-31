c     Routine to calculation forces and energy.
c     x and y are particle coordinates.
c     D1 is the diameter of small particles.
c     D is the diameter of particles in the unit of small particle diameter,
c     so D=1 for particles 1 to N / 2 and D=1.4 for particles N/2+1 to N.
c     N is the number of particles.  The particles are bidisperse with diameter
c     ratio of 1.4.
c     fx and fy are forces in x and y direction acting on particles.
c     V is the total potential energy.


      SUBROUTINE DF1grad(COORDS,N,VNEW,V,GTEST,SSTEST,BOXLX,BOXLY)
      USE KEY
      IMPLICIT NONE
      
      integer N, J1, I, J, J2
      LOGICAL GTEST, SSTEST
      double precision x(N), y(N), COORDS(3*N), VNEW(3*N), BOXLX, BOXLY
      double precision D(N), D1
      double precision fx(N), fy(N), V, fr
      double precision xij, yij, rij, dij

      D1=1.4D0
      DO J1=1,N/2
         D(J1)=1.0D0
      ENDDO
      DO J1=N/2+1,N
         D(J1)=D1
      ENDDO
C
C  Deal with any atoms that have left the box.
C
      IF ((.NOT.FIXIMAGE).AND.(.NOT.NORESET)) THEN
         DO J1=1,N
            J2 = 3*(J1-1)
            COORDS(J2+1)=COORDS(J2+1) - BOXLX*ANINT(COORDS(J2+1)/BOXLX)
            COORDS(J2+2)=COORDS(J2+2) - BOXLY*ANINT(COORDS(J2+2)/BOXLY)
         ENDDO
      ENDIF

      DO J1=1,N
         X(J1)=COORDS(3*(J1-1)+1)
         Y(J1)=COORDS(3*(J1-1)+2)
      ENDDO

      V = 0.D0

      do i = 1, N
         fx(i) = 0.D0
         fy(i) = 0.D0
      enddo
      
      do i = 1, N - 1
         do j = i + 1, N
            xij = x(i) - x(j)
            xij = xij - BOXLX*NINT(xij/BOXLX)
            yij = y(i) - y(j)
            yij = yij - BOXLY*NINT(yij/BOXLY)
            rij = dsqrt(xij * xij + yij * yij)
            dij = (D(i) + D(j)) / 2.D0
            if (rij .lt. dij) then
               fr = (1.D0 - rij / dij) / dij
               fx(i) = fx(i) + fr * xij / rij
               fx(j) = fx(j) - fr * xij / rij
               fy(i) = fy(i) + fr * yij / rij
               fy(j) = fy(j) - fr * yij / rij
              
               V = V + (1.D0 - rij / dij) ** 2 / 2.D0
            endif
         enddo
      enddo
      DO J1=1,N
         VNEW(3*(J1-1)+1)=-FX(J1)
         VNEW(3*(J1-1)+2)=-FY(J1)
      ENDDO

      IF (.NOT.SSTEST) RETURN
      CALL DF1hess(X,Y,D1,D,BOXLX,BOXLY,N)
 
      return
      end

c     Routine to calculation hessian matrix.
c     x and y are particle coordinates.
c     D1 is the diameter of small particles.
c     D is the diameter of particles in the unit of small particle diameter,
c     so D=1 for particles 1 to N / 2 and D=1.4 for particles N/2+1 to N.
c     N is the number of particles.  The particles are bidisperse with diameter
c     ratio of 1.4.
c     M is the returned hessian matrix.

      subroutine DF1hess(X, Y, D1, D, BOXLX,BOXLY,N)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N
      DOUBLE PRECISION X(N), Y(N), BOXLX, BOXLY
      double precision D(N), D1
      double precision xij, yij, rij, dij, tij, cij, coe
      integer i1, i2, j1, j2, I, J

      HESS(1:3*N,1:3*N)=0.0D0
      
      do i = 1, N - 1
         i1 = 3 * i - 2  ! x
         i2 = i1 + 1     ! y
         do j = i + 1, N
            xij = x(i) - x(j)
            xij = xij - BOXLX*NINT(xij/BOXLX)
            yij = y(i) - y(j)
            yij = yij - BOXLY*NINT(yij/BOXLY)
            rij = dsqrt(xij * xij + yij * yij)
            dij = (D(i) + D(j)) / 2.D0
            if(rij .lt. dij) then
               tij = - (1.D0 - rij / dij) / dij
               cij = 1.D0 / dij ** 2

               coe = - cij / rij ** 2 + tij / rij ** 3
               
               j1 = 3 * j - 2 ! x
               j2 = j1 + 1    ! y
               
               HESS(i1, j1) = coe * xij * xij - tij / rij
               HESS(i1, j2) = coe * xij * yij
               HESS(i2, j1) = HESS(i1, j2)
               HESS(i2, j2) = coe * yij * yij - tij / rij
               
               HESS(i1, i1) = HESS(i1, i1) - HESS(i1, j1)
               HESS(j1, j1) = HESS(j1, j1) - HESS(i1, j1)
               HESS(i1, i2) = HESS(i1, i2) - HESS(i1, j2)
               HESS(j1, j2) = HESS(j1, j2) - HESS(i1, j2)
               HESS(i2, i2) = HESS(i2, i2) - HESS(i2, j2)
               HESS(j2, j2) = HESS(j2, j2) - HESS(i2, j2)

            end if
         enddo
      enddo

      do i = 2, 3*N
         do j = 1, i - 1
            HESS(i, j) = HESS(j, i)
         enddo
      enddo
!     WRITE(*,'(2I6,G15.5)') ((I,J,HESS(I,J),I=3,3),J=3,3)
      
      return
      end
