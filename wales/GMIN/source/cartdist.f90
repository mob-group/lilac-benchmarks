module cartdist
use commons, only: natoms, box_params, box_paramsgrad

implicit none

! dj337: util modules for converting between absolute and fractional coordinate
! systems, building the H matrix, and computing the reciprocal lattice vectors

contains

! -----------------------------------------------------------------------------------
! VARIABLES
! x: atomic positions in absolute coordinates
! xfrac: atomic positions in fractional coordinates

! xr: absolute rigid body (COM+AA) coordinates
! xrfrac: fractional rigid body coordinates
!!! NOTE: AA coordinates are the same in absolute and fractional coordinate systems

! box_params: unit cell lengths and angles, in radians (a, b, c, alpha, beta, gamma)
! box_paramsgrad: gradient of energy wrt unit cell parameters

! H: matrix that transforms between fractional and absolute coordinates;
!    x = H*xfrac
!    first index corresponds to row, second to column
! H_grad: derivatives of the H matrix wrt cell parameters;
!    first index corresponds to row, second to column, third to cell parameter
! H_inverse: inverse of the H matrix
! -----------------------------------------------------------------------------------
! CONVERTS ATOMIC POSITIONS FROM ABSOLUTE TO FRACTIONAL COORDINATES.
! -----------------------------------------------------------------------------------
! Assumes atomistic coordinates and an orthorhombic unit cell.

       subroutine cart2frac_ortho(x, xfrac)

       implicit none

       integer                       :: j1, j3
       double precision, intent(in)  :: x(3*natoms)
       double precision, intent(out) :: xfrac(3*natoms)

       xfrac(:) = 0.0d0

       do j1 = 1, natoms
          j3 = 3*j1
          ! convert from absolute to fractional
          xfrac(j3-2:j3) = x(j3-2:j3)/box_params(1:3)
       enddo

       end subroutine cart2frac_ortho

! -----------------------------------------------------------------------------------
! Assumes atomistic coordinates and a tricilinic unit cell.

       subroutine cart2frac_tri(x, xfrac, H_inverse)

       implicit none

       integer                       :: j1, j3
       double precision, intent(in)  :: x(3*natoms), H_inverse(3,3)
       double precision, intent(out) :: xfrac(3*natoms)

       xfrac(:) = 0.0d0

       do j1 = 1, natoms
          j3 = 3*j1
          ! convert from absolute to fractional
          xfrac(j3-2:j3) = matmul(H_inverse, x(j3-2:j3))
       enddo

       end subroutine

! -----------------------------------------------------------------------------------
! Assumes rigid body coordinates and an orthorhombic unit cell.

       subroutine cart2frac_rb_ortho(nrigidbody, xr, xrfrac)

       implicit none

       integer                       :: j1, j3, degfreedoms
       integer, intent(in)           :: nrigidbody
       double precision, intent(in)  :: xr(3*natoms)
       double precision, intent(out) :: xrfrac(3*natoms)

       ! number of degrees of freedom
       degfreedoms = 6*nrigidbody
       xrfrac(:) = 0.0d0

       do j1 = 1, nrigidbody
          j3 = 3*j1
          ! convert COM coordinates from absolute to fractional
          xrfrac(j3-2:j3) = xr(j3-2:j3)/box_params(1:3)
       enddo

       ! AA coordinates are unchanged
       xrfrac(3*nrigidbody+1:degfreedoms) = xr(3*nrigidbody+1:degfreedoms)

       end subroutine cart2frac_rb_ortho

! -----------------------------------------------------------------------------------
! Assumes rigid body coordinates an a triclinic unit cell.

       subroutine cart2frac_rb_tri(nrigidbody, xr, xrfrac, H_inverse)

       implicit none

       integer                       :: j1, j3, degfreedoms
       integer, intent(in)           :: nrigidbody
       double precision, intent(in)  :: xr(3*natoms), H_inverse(3,3)
       double precision, intent(out) :: xrfrac(3*natoms)

       ! number of degrees of freedom
       degfreedoms = 6*nrigidbody
       xrfrac(:) = 0.0d0

       do j1 = 1, nrigidbody
          j3 = 3*j1
          ! convert COM coordinates from absolute to fractional
          xrfrac(j3-2:j3) = matmul(H_inverse, xr(j3-2:j3))
       enddo

       ! AA coordinates are unchanged
       xrfrac(3*nrigidbody+1:degfreedoms) = xr(3*nrigidbody+1:degfreedoms)

       end subroutine cart2frac_rb_tri

! -----------------------------------------------------------------------------------
! CONVERTS ATOMIC POSITIONS FROM FRACTIONAL TO ABSOLUTE COORDINATES.
! -----------------------------------------------------------------------------------
! Assumes atomistic coordinates and an orthorhombic unit cell.

       subroutine frac2cart_ortho(x, xfrac)

       implicit none

       integer                       :: j1, j3
       double precision, intent(in)  :: xfrac(3*natoms)
       double precision, intent(out) :: x(3*natoms)

       x(:) = 0.0d0

       do j1 = 1, natoms
          j3 = 3*j1
          ! convert from fractional to absolute
          x(j3-2:j3) = xfrac(j3-2:j3)*box_params(1:3)
       enddo

       end subroutine frac2cart_ortho

! -----------------------------------------------------------------------------------
! Assumes atomistic coordinates and a triclinic unit cell.

       subroutine frac2cart_tri(x, xfrac, H)

       implicit none

       integer :: j1, j3
       double precision, intent(in)  :: xfrac(3*natoms), H(3,3)
       double precision, intent(out) :: x(3*natoms)

       x(:) = 0.0d0

       do j1 = 1, natoms
          j3 = 3*j1
          ! convert from fractional to absolute
          x(j3-2:j3) = matmul(H, xfrac(j3-2:j3))
      enddo

      end subroutine frac2cart_tri

! -----------------------------------------------------------------------------------
! Assumes rigid body coordinates and an orthorhombic unit cell.

       subroutine frac2cart_rb_ortho(nrigidbody, xr, xrfrac)

       implicit none

       integer                       :: j1, j3, degfreedoms
       integer, intent(in)           :: nrigidbody
       double precision, intent(in)  :: xrfrac(3*natoms)
       double precision, intent(out) :: xr(3*natoms)

       ! number of degfrees of freedom
       degfreedoms = 6*nrigidbody
       xr(:) = 0.0d0

       do j1 = 1, nrigidbody
          j3 = 3*j1
          ! convert COM coordinates from fractional to absolute
          xr(j3-2:j3) = xrfrac(j3-2:j3)*box_params(1:3)
       enddo

       ! AA coordinates are unchanged
       xr(3*nrigidbody+1:degfreedoms) = xrfrac(3*nrigidbody+1:degfreedoms)

       end subroutine frac2cart_rb_ortho

! -----------------------------------------------------------------------------------
! Assumes rigid body coordinates and a triclinic unit cell.

       subroutine frac2cart_rb_tri(nrigidbody, xr, xrfrac, H)

       implicit none

       integer                       :: j1, j3, degfreedoms
       integer, intent(in)           :: nrigidbody
       double precision, intent(in)  :: xrfrac(3*natoms), H(3,3)
       double precision, intent(out) :: xr(3*natoms)

       ! number of degrees of freedom
       degfreedoms = 6*nrigidbody
       xr(:) = 0.0d0

       do j1 = 1, nrigidbody
          j3 = 3*j1
          ! convert COM coordinates from fractional to absolute
          xr(j3-2:j3) = matmul(H, xrfrac(j3-2:j3))
       enddo

       ! AA coordinates are unchanged
       xr(3*nrigidbody+1:degfreedoms) = xrfrac(3*nrigidbody+1:degfreedoms)

       end subroutine frac2cart_rb_tri

! -----------------------------------------------------------------------------------

! BUILDS THE H MATRIX that transforms between fractional and absolute coordinates.
! If GTEST is true, computes the six derivative matrices of the H matrix with respects
! to the six cell parameters. This works for any triclinic unit cell.

      subroutine build_H(H, H_grad, gtest)

      implicit none

      double precision, intent(out)          :: H(3,3), H_grad(3,3,6)
      logical, intent(in)                    :: gtest
      !double precision, intent(in), optional :: box_parameters
      double precision                       :: box_lengths(3), box_angles(3)
      double precision                       :: c(3), s(3), v

      H(:,:) = 0.0d0
      H_grad(:,:,:) = 0.0d0
      box_lengths(:) = box_params(1:3)
      box_angles(:) = box_params(4:6)

      ! cosine of the angles
      c(:) = dcos(box_angles(:))
      ! sine of the angles
      s(:) = dsin(box_angles(:))
      ! factor that is related to the volume (but not quite volume)
      v = dsqrt(1.0d0 - c(1)**2 - c(2)**2 - c(3)**2 + 2.0d0*c(1)*c(2)*c(3))

      ! define the H transformation matrix
      ! first row of matrix
      H(1,1) = box_lengths(1)
      H(1,2) = box_lengths(2)*c(3)
      H(1,3) = box_lengths(3)*c(2)
      ! second row
      H(2,2) = box_lengths(2)*s(3)
      H(2,3) = box_lengths(3)*(c(1) - c(2)*c(3))/s(3)
      ! third row
      H(3,3) = box_lengths(3)*v/s(3)

      ! compute derivatives of H matrix
      if (gtest) then
         ! wrt box length a
         H_grad(1,1,1) = 1.0d0
         ! wrt box length b
         H_grad(1,2,2) = c(3)
         H_grad(2,2,2) = s(3)
         ! wrt box length c
         H_grad(1,3,3) = c(2)
         H_grad(2,3,3) = (c(1) - c(2)*c(3))/s(3)
         H_grad(3,3,3) = v/s(3)
         ! wrt box angle alpha
         H_grad(2,3,4) = -box_lengths(3)*s(1)/s(3)
         H_grad(3,3,4) = box_lengths(3)*(c(1)*s(1) - s(1)*c(2)*c(3))/(s(3)*v)
         ! wrt box angle beta
         H_grad(1,3,5) = -box_lengths(3)*s(2)
         H_grad(2,3,5) = box_lengths(3)*s(2)*c(3)/s(3)
         H_grad(3,3,5) = box_lengths(3)*s(2)*(c(2) - c(1)*c(3))/(s(3)*v)
         ! wrt box angle gamma
         H_grad(1,2,6) = -box_lengths(2)*s(3)
         H_grad(2,2,6) = box_lengths(2)*c(3)
         H_grad(2,3,6) = box_lengths(3)*(c(2) - c(1)*c(3))/s(3)**2
         H_grad(3,3,6) = box_lengths(3)*((c(3) - c(1)*c(2))/v - v*c(3)/s(3)**2)
      endif

      return
      end subroutine build_H

! -----------------------------------------------------------------------------------

! COMPUTES CELL VOLUME. This works for any orthorhombic or triclinic unit cell.
! VARIABLES
! vol: cell volume

       subroutine get_volume(vol)

       use commons, only: ortho

       implicit none

       !double precision, intent(in), optional :: box_parameters(6)
       double precision, intent(out)          :: vol
       double precision                       :: box_lengths(3), box_angles(3), c(3)

       box_lengths(:) = box_params(1:3)
       box_angles(:) = box_params(4:6)

       vol = box_lengths(1)*box_lengths(2)*box_lengths(3)
       if (.not.ortho) then
          c(:) = dcos(box_angles(:))
          vol = vol * dsqrt(1.0d0 - c(1)**2 - c(2)**2 - c(3)**2 + 2.0d0*c(1)*c(2)*c(3))
       endif

       end subroutine get_volume

! -----------------------------------------------------------------------------------

! BUILDS THE K MATRIX whose columns are the reciprocal lattice vectors. If GTEST is
! true, computes the six derivative matrices of the K matrix with respects to the six
! cell parameters. This works for any triclinic unit cell.
! VARIABLES
! reciplatvec: matrix whose columns are the reciprocal lattice vectors
!    first index corresponds to row, second to column
! reciplatvec_grad: derivatives of the reciprocal lattice vector matrix wrt cell parameters
!    first index corresponds to row, second to column, third to cell parameter

      subroutine get_reciplatvec(reciplatvec, reciplatvec_grad, gtest)

      implicit none

      double precision, intent(out) :: reciplatvec(3,3), reciplatvec_grad(3,3,6)
      logical, intent(in)           :: gtest
      double precision              :: box_lengths(3), box_angles(3)
      double precision              :: c(3), s(3), v, dv(3), cfact
      double precision, parameter   :: pi = 3.141592654d0

      reciplatvec(:,:) = 0.0d0
      reciplatvec_grad(:,:,:) = 0.0d0
      box_lengths(:) = box_params(1:3)
      box_angles(:) = box_params(4:6)
       
      ! cosine of the angles
      c(:) = dcos(box_angles(:))
      ! sine of the angles
      s(:) = dsin(box_angles(:))
      ! factor that is related to the volume (but not quite volume)
      v = dsqrt(1.0d0 - c(1)**2 - c(2)**2 - c(3)**2 + 2.0d0*c(1)*c(2)*c(3))

      ! define the reciprocal lattice vector matrix
      ! first row of matrix
      reciplatvec(1,1) = 1.0d0/box_lengths(1)
      ! second row
      reciplatvec(2,1) = -c(3)/(box_lengths(1)*s(3))
      reciplatvec(2,2) = 1.0d0/(box_lengths(2)*s(3))
      ! third row
      reciplatvec(3,1) = (c(3)*(c(1) - c(2)*c(3)) - c(2)*s(3)**2)/(box_lengths(1)*v*s(3))
      reciplatvec(3,2) = -(c(1) - c(2)*c(3))/(box_lengths(2)*v*s(3))
      reciplatvec(3,3) = s(3)/(box_lengths(3)*v)
      ! multiply by 2*pi
      reciplatvec(:,:) = 2.0d0*pi*reciplatvec(:,:)

      ! compute derivatives of reciprocal lattice vector matrix
      if (gtest) then
         ! gradient of v wrt cell angles
         dv(1) = s(1)*(c(1)-c(2)*c(3))/v
         dv(2) = s(2)*(c(2)-c(1)*c(3))/v
         dv(3) = s(3)*(c(3)-c(1)*c(2))/v
         ! cosine factor: cos(alpha) - cos(beta)cos(gamma)
         cfact = c(1)-c(2)*c(3)

         ! wrt box length a
         reciplatvec_grad(1,1,1) = -1.0d0/box_lengths(1)**2 
         reciplatvec_grad(2,1,1) = c(3)/(box_lengths(1)**2*s(3))
         reciplatvec_grad(3,1,1) = (c(2)*s(3) - c(3)*cfact/s(3))/(box_lengths(1)**2*v)
         ! wrt box length b
         reciplatvec_grad(2,2,2) = -1.0d0/(box_lengths(2)**2*s(3))
         reciplatvec_grad(3,2,2) = cfact/(box_lengths(2)**2*s(3)*v)
         ! wrt box length c
         reciplatvec_grad(3,3,3) = -s(3)/(box_lengths(3)**2*v)
         ! wrt cell angle alpha
         reciplatvec_grad(3,1,4) = (-c(3)*s(1)/s(3) - c(3)*cfact*dv(1)/(s(3)*v) + c(2)*s(3)*dv(1)/v)/(box_lengths(1)*v)
         reciplatvec_grad(3,2,4) = (s(1) + cfact*dv(1)/v)/(box_lengths(2)*s(3)*v)
         reciplatvec_grad(3,3,4) = -s(3)*dv(1)/(box_lengths(3)*v**2)
         ! wrt cell angle beta
         reciplatvec_grad(3,1,5) = ((s(2)*c(3)**2)/s(3) - c(3)*cfact*dv(2)/(s(3)*v) + &
                                   s(2)*s(3) + c(2)*s(3)*dv(2)/v)/(box_lengths(1)*v)
         reciplatvec_grad(3,2,5) = (-s(2)*c(3) + (cfact*dv(2)/v))/(box_lengths(2)*s(3)*v)
         reciplatvec_grad(3,3,5) = -s(3)*dv(2)/(box_lengths(3)*v**2)
         ! wrt cell angle gamma
         reciplatvec_grad(2,1,6) = 1.0d0/(box_lengths(1)*s(3)**2)
         reciplatvec_grad(2,2,6) = -c(3)/(box_lengths(2)*s(3)**2)
         reciplatvec_grad(3,1,6) = (-cfact - c(3)**2*cfact/(s(3)**2) - c(3)*cfact*dv(3)/(s(3)*v) + &
                                   c(2)*s(3)*dv(3)/v)/(box_lengths(1)*v)
         reciplatvec_grad(3,2,6) = (-c(2)*s(3) + c(3)*cfact/s(3) + cfact*dv(3)/v)/(box_lengths(2)*s(3)*v)
         reciplatvec_grad(3,3,6) = (c(3) - s(3)*dv(3)/v)/(box_lengths(3)*v)

         ! multiply by 2*pi
         reciplatvec_grad(:,:,:) = 2.0d0*pi*reciplatvec_grad(:,:,:)
      endif

      end subroutine get_reciplatvec

! -----------------------------------------------------------------------------------

end module
