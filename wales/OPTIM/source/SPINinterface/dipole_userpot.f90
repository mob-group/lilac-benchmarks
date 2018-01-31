!File containing all important routines for a system of dipole-coupled 
!nanoparticles.

! create the lattice and calculate the ewald-sums for the interaction
subroutine userpot_init
	use main_mod
	use list_mod
	use init_mod
	implicit none
	
	call init_lattice(.false.,'quad')
	return
end subroutine userpot_init	

! return the number of atoms
subroutine userpot_get_natoms(num_atoms)
	use main_mod
	implicit none
	integer, intent(out) :: num_atoms

	num_atoms = nats
	return
end subroutine

! copy all coordinates to gmin coords(:,1) array (in commons)
! TODO: what about MPI?
! for very simple standard cases this is just initial configuration,
! but can be also more complicated, e.g. initalizing generalized rigid
! body framework
subroutine userpot_initialize_gmin(dof,x)
	use init_mod
	implicit none
	integer, intent(in) :: dof
	double precision, intent(out) :: x(dof)
	
	call init_shiftspins(x)
	return
end subroutine

! called by gmin to calculate potential
subroutine userpot_potential(dof,X,GRAD,EREAL,GRADT)
	use main_mod
	use force_mod
	implicit none
	integer, intent(in) :: dof                 ! number of degrees of freedom
	double precision, intent(in) :: X(dof)     ! current coordinates
	double precision, intent(out) :: GRAD(dof) ! gradient
	double precision, intent(out) :: EREAL     ! energy
	logical, intent(in) :: gradt               ! is the gradient needed?
	
	if (dof .ne. n) then
		print *,'Failure in initializing userpot_potential!'
		stop
	endif
	call get_energy(x,ereal)
	print *,x
	print *,"energy", ereal
	if (gradt) then
		call get_gradient(x,grad)
	endif
	return
end subroutine

!called by optim to calculate potential
subroutine userpot_potentialhess(DOF,X,GRAD,EREAL,GRADT,SECT,HESS)
	use main_mod
	use force_mod
	implicit none
	integer, intent(in) :: DOF                 ! number of degrees of freedom
	double precision, intent(in) :: X(DOF)     ! current coordinates
	double precision, intent(out) :: GRAD(DOF) ! gradient
	double precision, intent(out) :: EREAL     ! energy
	logical, intent(in) :: GRADT               ! is the gradient needed?
	logical, intent(in) :: SECT		   ! is the hessian needed?
	double precision, intent(out) :: HESS(DOF,DOF)

	if (DOF .ne. n) then
		print *,'Failure in initializing userpot_potential!'
		stop
	endif
	call get_energy(x,ereal)
	if (gradt) then
		call get_gradient(x,grad)
	endif
	if (sect) then
		call get_hessian(x,hess)
	endif
	return
end subroutine

subroutine userpot_distance(vec1,vec2,dist)
	use main_mod
	use math_mod
	implicit none
	double precision, intent(in) :: vec1(n)
	double precision, intent(in) :: vec2(n)
	double precision, intent(out) :: dist
	double precision :: rvec1(n),rvec2(n)

	rvec1=vec1
	rvec2=vec2
	call math_rescale(rvec1)
	call math_rescale(rvec2)
	dist=math_dis(rvec1,rvec2)
	return
end subroutine userpot_distance

subroutine userpot_delta(x,y,delta)
	use main_mod
	use math_mod
	implicit none
	double precision, intent(in) :: x(n)
	double precision, intent(in) :: y(n)
	double precision :: delta(n),h1(n),h2(n)
	integer :: i
	
	h1=x
	h2=y
	
	call math_rescale(h1)
	call math_rescale(h2)

	do i=1,3*nats
		if (abs(h1(i)-h2(i)) .le. pi) then
			delta(i)=y(i)-x(i)
		elseif (h1(i)-h2(i) .lt. -pi) then
			delta(i)=2.d0*pi+(h1(i)-h2(i))
		elseif (h1(i)-h2(i) .gt. pi) then
			delta(i)=2.d0*pi-(h1(i)-h2(i))
		endif
	enddo
	return
end subroutine userpot_delta	

	


	

