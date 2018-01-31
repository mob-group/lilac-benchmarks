module box_derivatives
use commons, only: natoms, box_params, box_paramsgrad
use genrigid, only: degfreedoms, nrigidbody
use cartdist

implicit none

public :: check_angles

! dj337: module to convert gradient from absolute to 
! fractional coordinates to perform basin-hopping on
! periodic systems with varying lattice parameters

contains

! -----------------------------------------------------------------------------------
! VARIABLES
! x: atomic positions in absolute coordinates
! xfrac: atomic positions in fractional coordinates
! grad: gradient of energy wrt absolute atomic positions
! gradfrac: gradient of energy wrt fractional atomic positions

! xr: absolute rigid body (COM+AA) coordinates
! xrfrac: fractional rigid body coordinates
! gradr: gradient of energy wrt absolute rigid body coordinates
! gradrfrac: gradient of energy wrt fractional rigid body coordinates
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

! CONVERTS ATOMIC POSITIONS FROM ABSOLUTE TO FRACTIONAL COORDINATES. Wrapper subroutine
! that calls appropriate subroutine depending on whether atomistic or rigidbody
! coordinates, orthorhombic or tricilinc unit cell.

      subroutine cart2frac(x, xfrac, H_inverse)

      use commons, only: ortho
      use genrigid, only: rigidinit, atomrigidcoordt, inversematrix

      implicit none

      double precision, intent(in)           :: x(3*natoms)
      double precision, intent(out)          :: xfrac(3*natoms)
      double precision, intent(in), optional :: H_inverse(3,3)
      double precision                       :: H(3,3), H_grad(3,3,6), H_inv(3,3)

      ! orthorhombic cells
      if (ortho) then
         if (rigidinit.and.(.not.atomrigidcoordt)) then
            call cart2frac_rb_ortho(nrigidbody, x, xfrac)
         else
            call cart2frac_ortho(x, xfrac)
         endif
      ! triclinic cells
      else
         ! compute inverse of H matrix if not given
         if (present(H_inverse)) then
            H_inv = H_inverse
         else
            call build_H(H, H_grad, .false.)
            call inversematrix(H, H_inv)
         endif
         ! convert coordinates
         if (rigidinit.and.(.not.atomrigidcoordt)) then
            call cart2frac_rb_tri(nrigidbody, x, xfrac, H_inv)
         else
            call cart2frac_tri(x, xfrac, H_inv)
         endif
      endif

      end subroutine cart2frac

! -----------------------------------------------------------------------------------

! CONVERTS ATOMIC POSITIONS FROM FRACTIONAL TO ABSOLUTE COORDINATES. Wrapper subroutine
! that calls appropriate subroutine depending on whether atomistic or rigidbody
! coordinates, orthorhombic or tricilinc unit cell.

      subroutine frac2cart(x, xfrac, H)

      use commons, only: ortho
      use genrigid, only: rigidinit, atomrigidcoordt

      implicit none

      double precision, intent(in)           :: xfrac(3*natoms)
      double precision, intent(out)          :: x(3*natoms)
      double precision, intent(in), optional :: H(3,3)
      double precision                       :: H_mat(3,3), H_grad(3,3,6)

      ! orthorhombic cell
      if (ortho) then
         if (rigidinit.and.(.not.atomrigidcoordt)) then
            call frac2cart_rb_ortho(nrigidbody, x, xfrac)
         else
            call frac2cart_ortho(x, xfrac)
         endif
      ! triclinic cell
      else
         ! compute H matrix if not given
         if (present(H)) then
            H_mat = H
         else
            call build_H(H_mat, H_grad, .false.)
         endif
         ! convert coordinates
         if (rigidinit.and.(.not.atomrigidcoordt)) then
            call frac2cart_rb_tri(nrigidbody, x, xfrac, H_mat)
         else
            call frac2cart_tri(x, xfrac, H_mat)
         endif
      endif

      end subroutine frac2cart
 
! -----------------------------------------------------------------------------------

! CONVERTS GRADIENT WRT ATOMIC POSITIONS FROM ABSOLUTE TO FRACTIONAL COORDINATES.
! Wrapper subroutine that calls appropriate subroutine depending on whether atomistic
! or rigidbody coordinates, orthorhombic or triclinic unit cell.

      subroutine boxderiv(grad, gradfrac, H)

      use commons, only: ortho
      use genrigid, only: rigidinit, atomrigidcoordt

      implicit none

      double precision, intent(in)           :: grad(3*natoms)
      double precision, intent(out)          :: gradfrac(3*natoms)
      double precision, intent(in), optional :: H(3,3)
      double precision                       :: H_mat(3,3), H_grad(3,3,6)

      ! orthorhombic cell
      if (ortho) then
         if (rigidinit.and.(.not.atomrigidcoordt)) then
            call boxderiv_rb_ortho(grad, gradfrac)
         else
            call boxderiv_ortho(grad, gradfrac)
         endif
      ! triclinic
      else
         ! compute H matrix if not given
         if (present(H)) then
            H_mat = H
         else
            call build_H(H_mat, H_grad, .false.)
         endif
         ! convert gradient
         if (rigidinit.and.(.not.atomrigidcoordt)) then
            call boxderiv_rb_tri(grad, gradfrac, H_mat)
         else
            call boxderiv_tri(grad, gradfrac, H_mat)
         endif
      endif

      end subroutine boxderiv

! -----------------------------------------------------------------------------------
! Assumes atomistic coordinates and an orthorhombic unit cell.

       subroutine boxderiv_ortho(grad, gradfrac)

       implicit none

       integer                       :: j1, j3
       double precision, intent(in)  :: grad(3*natoms)
       double precision, intent(out) :: gradfrac(3*natoms)

       gradfrac(:) = 0.0d0
       !box_paramsgrad(:) = 0.0d0

       ! iterate over atoms
       do j1 = 1, natoms
          j3 = 3*j1

          ! convert gradient wrt atom positions from absolute to fractional
          gradfrac(j3-2:j3) = grad(j3-2:j3)*box_params(1:3)

          ! TODO: trying to compute box derivatives generally (doesn't work!!)
          ! add contribution to gradient wrt cell lengths
          ! box_paramsgrad(1:3) = box_paramsgrad(1:3) + xfrac(j3-2:j3)*gradfrac(j3-2:j3)/box_params(1:3)
       enddo

       end subroutine boxderiv_ortho

! -----------------------------------------------------------------------------------
! Assumes atomistic coordinates and a triclinic unit cell.

       subroutine boxderiv_tri(grad, gradfrac, H)

       implicit none

       integer                       :: j1, j3
       double precision, intent(in)  :: grad(3*natoms), H(3,3)
       double precision, intent(out) :: gradfrac(3*natoms)

       gradfrac(:) = 0.0d0

       do j1 = 1, natoms
          j3 = 3*j1
          ! convert from absolute to fractional
          gradfrac(j3-2:j3) = matmul(grad(j3-2:j3), H)
       enddo

       end subroutine boxderiv_tri

! -----------------------------------------------------------------------------------
! Assumes rigid body coordinates and an orthorhombic unit cell.

       subroutine boxderiv_rb_ortho(gradr, gradrfrac)

       implicit none

       integer                       :: j1, j3
       double precision, intent(in)  :: gradr(3*natoms)
       double precision, intent(out) :: gradrfrac(3*natoms)

       gradrfrac(:) = 0.0d0

       do j1 = 1, nrigidbody
          j3 = 3*j1
          ! convert gradient wrt COM positions from absolute to fractional
          gradrfrac(j3-2:j3) = gradr(j3-2:j3)*box_params(1:3)
       enddo

       ! gradient wrt AA coordinates are unchanged
       gradrfrac(3*nrigidbody+1:degfreedoms) = gradr(3*nrigidbody+1:degfreedoms)

       end subroutine boxderiv_rb_ortho

! -----------------------------------------------------------------------------------
! Assumes rigid body coordinates and a triclinic unit cell.

       subroutine boxderiv_rb_tri(gradr, gradrfrac, H)

       implicit none

       integer :: j1, j3
       double precision, intent(in)  :: gradr(3*natoms), H(3,3)
       double precision, intent(out) :: gradrfrac(3*natoms)

       gradrfrac(:) = 0.0d0

       do j1 = 1, nrigidbody
          j3 = 3*j1
          ! convert gradient wrt COM positions from absolute to fractional
          gradrfrac(j3-2:j3) = matmul(gradr(j3-2:j3), H)
       enddo

       ! gradient wrt AA coordinates are unchanged
       gradrfrac(3*nrigidbody+1:degfreedoms) = gradr(3*nrigidbody+1:degfreedoms)

       end subroutine boxderiv_rb_tri

! -----------------------------------------------------------------------------------
! TAKES A BASIN-HOPPING STEP by making uniformly random changes to the cell lengths
! and angles

       subroutine bd_takestep(np)

       use commons, only: ortho, box_params
       use vec3, only: vec_random

       integer, intent(in)         :: np
       double precision            :: new_angles(3)
       double precision, parameter :: max_length_step = 0.3d0
       double precision, parameter :: max_angle_step = 0.1d0

       ! generate new box lengths
       box_params(1:3) = box_params(1:3) + vec_random() * max_length_step
       ! if triclinic
       if (.not.(ortho)) then
          new_angles(:) = box_params(4:6) + vec_random() * max_angle_step
          ! check to make sure combination of new angles is valid
          do while (.not.check_angles(new_angles(:)))
             new_angles(:) = box_params(4:6) + vec_random() * max_angle_step
          enddo
          box_params(4:6) = new_angles(:)
       endif ! triclinic

       end subroutine bd_takestep

! -----------------------------------------------------------------------------------
! CHECKS SET OF TRICLINIC CELL ANGLES to make sure they are valid. Non-valid set of angles
! corresponds to a cell with zero or imaginary volume. Function returns True if the set
! of angles meets the criteria for valid cell angles.

       pure logical function check_angles(angles)

       implicit none

       double precision, intent(in) :: angles(3)
       double precision             :: sums(4)
       double precision, parameter   :: pi = 3.141592654d0

       ! calculate necessary sums
       sums(1) =  angles(1) + angles(2) + angles(3)
       sums(2) = -angles(1) + angles(2) + angles(3)
       sums(3) =  angles(1) - angles(2) + angles(3)
       sums(4) =  angles(1) + angles(2) - angles(3)

       ! check all sums are between 0 and 2pi and all angles between 0 and pi
       check_angles = all(sums.gt.0.0d0).and.all(sums.lt.2*pi).and.all(angles.gt.0.0d0).and.all(angles.lt.pi)
       end function check_angles

! -----------------------------------------------------------------------------------
! REJECTS structures that are invalid. Tricks the LBFGS minimizer by returning a very
! large energy and a very small gradient, so will the quench will immediately be 
! considered converged but the structure will never be saved as a low-energy minimum.

       subroutine reject(energy, grad)

       implicit none

       double precision, intent(out) :: energy, grad(3*natoms)

       energy = 1.0d20
       grad(:) = 1.0d-20

       end subroutine reject

! -----------------------------------------------------------------------------------
! ADDS WCA-STYLE REPULSION term to the energy and gradients to repel the cell away
! from having zero volume.

       subroutine constrain_volume(v, v_deriv, energy, grad_angles, gtest)

       implicit none

       double precision, intent(in)    :: v, v_deriv(3)
       double precision, intent(inout) :: energy, grad_angles(3)
       logical, intent(in)             :: gtest
       double precision, parameter     :: v_sigma = 3.0d-1
       double precision, parameter     :: v_eps = 1.0d-3

       if (v.lt.v_sigma**(1.0d0/6.d0)) then
          ! add purely repulsive WCA energy term
          energy = energy + 4.0d0*v_eps*((v_sigma/v)**12 - (v_sigma/v)**6) + v_eps

          if (gtest) then
             ! add gradient contribution
             grad_angles(:) = grad_angles(:) + 24.0d0*v_eps/v_sigma*((v_sigma/v)**7 - 2.0d0*(v_sigma/v)**13)*v_deriv(:)
          endif
       endif

       end subroutine constrain_volume

! ---------------------------------------------------------------------------------
! 
! Rotates all rigid bodies after a step is taken in the cell parameters.
! VARIABLES
! box_paramsold: unit cell lengths and angles from before the step was taken
! TODO: not sure this is working properly (or that the equations are even valid...)


!       subroutine rotate_bodies(box_paramsold, xrfrac)
!
!       use genrigid, only: transformctorigid, sitesrigidbody, maxsite, nsiteperbody, inversematrix
!
!       integer                         :: j1, j3, j5, j2, j4, j6
!       double precision, intent(in)    :: box_paramsold(6)
!       double precision, intent(inout) :: xrfrac(3*natoms)
!       double precision                :: H_old(3,3), H_grad(3,3,6), H_oldinverse(3,3), H_new(3,3)
!       double precision                :: ri(3), p(3), xr(3*natoms), x(3*natoms), rot(3,3)
!       double precision                :: rmi(3,3), drmi1(3,3), drmi2(3,3), drmi3(3,3), vol_new, vol_old
!
!       print *, 'rotating!'
!
!       !print *, 'xrfrac old: ', xrfrac(:degfreedoms)
!       ! get H matrix from old box parameters
!       call build_H(H_old, H_grad, .false., box_paramsold)
!       !print *, 'H_old: ', H_old(:3,:3)
!       call inversematrix(H_old, H_oldinverse)
!       call get_volume(vol_old, box_paramsold)
!       ! get H matrix from current box parameters
!       call build_H(H_new, H_grad, .false.)
!       call get_volume(vol_new)
!
!       !print *, 'H_old: ', H_old(:3,:3)
!       !print *, 'H_new: ', H_new(:3,:3)
! 
!       call frac2cart_rb_tri(xr, xrfrac, H_old)
!
!       do j1 = 1, 1!nrigidbody
!          j3 = 3*j1
!          j5 = 3*nrigidbody + j3
!
!          ri(:) = xr(j3-2:j3)
!          p(:) = xr(j5-2:j5)
!
!          call rmdrvt(p, rmi, drmi1, drmi2, drmi3, .false.)
!          ! get new rotation matrix
!          rot(:,:) = (vol_old/vol_new)*matmul(H_new, matmul(H_oldinverse, rmi))
!          !print *, 'rmi: ', rmi(:3,:3)
!          !print *, 'rot: ', rot(:3,:3)
!
!          do j2 = 1, nsiteperbody(j1)
!             j4 = maxsite*(j1-1) + j2
!             j6 = 3*j4
!             x(j6-2:j6) = ri(:) + matmul(rot(:,:), sitesrigidbody(j2,:,j1))
!          enddo
!          
!          print *, 'reference geometry: '
!          do j2 = 1, nsiteperbody(j1)
!             print *, sitesrigidbody(j2,:,j1)
!          enddo
!          print *, 'beginning rotation: '
!          do j2 = 1, nsiteperbody(j1)
!             print *, (ri(:) + matmul(rmi, sitesrigidbody(j2,:,j1)))
!          enddo
!          print *, 'base rotation: '
!          do j2 = 1, nsiteperbody(j1)
!             print *, (ri(:) + vol_old*matmul(matmul(H_oldinverse, rmi), sitesrigidbody(j2,:,j1)))
!          enddo
!          print *, 'end rotation: '
!          do j2 = 1, nsiteperbody(j1)
!             print *, (ri(:) + matmul(rot, sitesrigidbody(j2,:,j1)))
!          enddo
!
!       enddo
!
!       !print *, 'x         : ', x(:3*natoms)
!       call transformctorigid(x, xrfrac)
!       !print *, 'xrfrac new: ', xrfrac(:degfreedoms)
!
!       end subroutine rotate_bodies

! -----------------------------------------------------------------------------------
end module
