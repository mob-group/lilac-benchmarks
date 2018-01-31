subroutine userpot_init
    call load_plugin
end subroutine

! return the number of atoms
subroutine userpot_get_natoms(num_atoms)
    call plugin_get_natoms
end subroutine

subroutine userpot_potentialhess(dof,X,GRAD,EREAL,GRADT,SECT,HESS)
    integer, intent(in) :: dof			!degrees of freedom
    double precision, intent(in) :: X(dof)	!current coordinates
    double precision, intent(out) :: GRAD(dof)	!gradient
    double precision, intent(out) :: EREAL	!energy
    logical, intent(in) :: GRADT		!is the gradient needed?
    logical, intent(in) :: SECT 		!is the hessian needed?
    logical, intent(out) :: HESS(dof,dof)	!hessian

    call plugin_potential(x, ereal, grad, hess, gradt, sect)

end subroutine userpot_potentialhess

subroutine userpot_distance(dof,x,y,dist)
    implicit none
    integer, intent(in) :: dof
    double precision, intent(in) :: x(dof)
    double precision, intent(in) :: y(dof)
    double precision, intent(out) :: dist
    
    print *,'USERPOT distance not supported yet'
    stop
end subroutine userpot_distance

