!dummyfile for an external userpotential. If you want to implement a
!userpotential to OPTIM/GMIN, you have to provide the routines
!below.
!Some minor fixes have to be done to some routines:
!	- fetchz.f -> IF (USERPOTT) THEN LNATOMS=NATOMS
!This should be changed if possible!

subroutine userpot_init
    print *,'ERROR: you are using USERPOT with a standard GMIN binary'
    stop
end subroutine

! return the number of atoms
subroutine userpot_get_natoms(num_atoms)
    integer, intent(out) :: num_atoms
    print *,'ERROR: you are using USERPOT with a standard GMIN binary'
    stop
end subroutine

! copy all coordinates to gmin coords(:,1) array (in commons)
! TODO: what about MPI?
! for very simple standard cases this is just initial configuration,
! but can be also more complicated, e.g. initalizing generalized rigid
! body framework
subroutine userpot_initialize_gmin(dof, x)
    integer, intent(in) :: dof
    double precision, intent(out) :: x(dof)
    print *,'ERROR: you are using dmacrys with a non-dmacrys binary'
    stop
end subroutine

! called by gmin to calculate potential
subroutine userpot_potential(dof,X,GRAD,EREAL,GRADT)
    integer, intent(in) :: dof                 ! number of degrees of freedom
    double precision, intent(in) :: X(dof)     ! current coordinates
    double precision, intent(out) :: GRAD(dof) ! gradient
    double precision, intent(out) :: EREAL     ! energy
    logical, intent(in) :: gradt              ! is the gradient needed?
    print *,'ERROR: you are using dmacrys with a non-dmacrys binary'
    stop
end subroutine

subroutine userpot_potentialhess(dof,X,GRAD,EREAL,GRADT,SECT,HESS)
    integer, intent(in) :: dof			!degrees of freedom
    double precision, intent(in) :: X(dof)	!current coordinates
    double precision, intent(out) :: GRAD(dof)	!gradient
    double precision, intent(out) :: EREAL	!energy
    logical, intent(in) :: GRADT		!is the gradient needed?
    logical, intent(in) :: SECT 		!is the hessian needed?
    logical, intent(out) :: HESS(dof,dof)	!hessian

    print *,'ERROR: you are using dmacrys with a non-dmacrys binary'
    stop
end subroutine userpot_potentialhess

subroutine userpot_distance(dof,x,y,dist)
    implicit none
    integer, intent(in) :: dof
    double precision, intent(in) :: x(dof)
    double precision, intent(in) :: y(dof)
    double precision, intent(out) :: dist
    
    print *,'ERROR: you are using dmacrys with a non-dmacrys binary'
    stop
end subroutine userpot_distance

subroutine userpot_dump_configuration(filename, coords)
implicit none
double precision COORDS(*)
CHARACTER filename
    
end subroutine

subroutine userpot_dump_lowest
end subroutine
