! Called from keywords.f to initialize DMACRYS when specified
! in input file
subroutine dmacrys_setup
    use dmacrys_interface
    use commons
    implicit none
    integer dmacrys_get_natoms
    open (5,FILE='data',STATUS='OLD')
    print *,'DMACRYS INTERFACE> Using dmacrys, setting up interface'
    CALL dmacrys_initialize
    natoms = dmacrys_get_natoms()+2
    ALLOCATE(COORDS(3*NATOMS, 1))
    call dmacrys_initialize_gmin()
end subroutine

subroutine dmacrys_get_coords(x)
    use genrigid
    use commons, only: coords, natoms
    double precision :: x(degfreedoms)

    call transformctorigid(coords(:,1), x)
end subroutine

function dmacrys_get_dof()
    use genrigid, only: degfreedoms
    integer dmacrys_get_dof
    dmacrys_get_dof = degfreedoms
end function

function dmacrys_get_natoms()
    use dmacrys_interface
    integer dmacrys_get_natoms
    dmacrys_get_natoms = num_atoms
end function

subroutine dmacrys_initialize_gmin
    use dmacrys_interface
    print *,'DMACRYS INTERFACE> Setting up GMIN from DMACRYS structures'
    CALL dmacrys_init_gmin
end subroutine

! calculate the potential using dmacrys
subroutine dmacrys_potential(X,GRAD,EREAL,GRADT)
    use dmacrys_interface
    use commons
    use genrigid, only: degfreedoms
    implicit none
    double precision x(1:degfreedoms)
    double precision grad(1:degfreedoms)
    double precision ereal
    logical GRADT

    !print *,"calc potential"
    !print *,X
    !print *,"end"
    grad(1:degfreedoms) = 0
    ereal = 0
!    call dmacrys_set_coords(X)
!    EREAL=dmagrys_calc_energy(X, GRAD, GRADT)
    call dmacrys_calc_potential(X,GRAD,EREAL,GRADT)
!    print *,"Energy",ereal
!    if(dmadidstep) then
!        call dmacrys_check_derivatives(3*NATOMS,x,grad)
!    endif
    dmadidstep=.false.
end subroutine

subroutine dmacrys_dump(filename, coords)
    USE DMACRYS_INTERFACE, only : dmacrys_get_lattice,DMACRYS_DUMP_CIF_NO_TRANSFORM
    USE COMMONS, ONLY : NATOMS
    use genrigid, only : get_lattice_matrix
    IMPLICIT NONE
    INTEGER I
    DOUBLE PRECISION coords(3*NATOMS)
    CHARACTER*(*) :: filename
    DOUBLE PRECISION MLATTICE(3,3)
    call dmacrys_reduce_cell(NATOMS, coords)
    call get_lattice_matrix(coords(3*NATOMS-5:3*NATOMS), MLATTICE)
    CALL DMACRYS_DUMP_CIF_NO_TRANSFORM(filename,&
            NATOMS-2,coords(1:3*NATOMS-6), MLATTICE, "noname")
end subroutine

subroutine dmacrys_reduce_cell_rigid(rcoords,cell_was_changed)
    use genrigid
    use cell_reduction
    use vec3
    use commons, only: natoms
    implicit none
    integer i
    double precision rcoords(DEGFREEDOMS)
    double precision coords(1:3*NATOMS)
    double precision lattice(3,3), tmp(3,3),transform(3,3),lattice_inv(3,3)
    logical cell_was_changed

    cell_was_changed = .false.
    call get_lattice_matrix(rcoords(DEGFREEDOMS-5:DEGFREEDOMS), lattice)

    call cell_set_lattice(lattice)
    call cell_minimum_reduction

    if(.not.cell_has_changed()) return
    cell_was_changed = .true.

    call invert3x3(lattice, lattice_inv)

    tmp = cell_get_transformation()
    call cell_get_lattice(lattice)

    transform = matmul(lattice,matmul(tmp, lattice_inv))
    call TRANSFORMRIGIDTOC (1, NRIGIDBODY, coords, rcoords(1:DEGFREEDOMS))
    do i=1,natoms
        coords(3*i-2:3*i) = matmul(transform, coords(3*i-2:3*i))
    end do
    coords(3*NATOMS-5) = lattice(1,1)
    coords(3*NATOMS-4) = lattice(2,1)
    coords(3*NATOMS-3) = lattice(3,1)
    coords(3*NATOMS-2) = lattice(2,2)
    coords(3*NATOMS-1) = lattice(3,2)
    coords(3*NATOMS-0) = lattice(3,3)

    call TRANSFORMCTORIGID (coords, rcoords(1:DEGFREEDOMS))

    do i=1,3*NRIGIDBODY
        rcoords(i) = rcoords(i) - entier(rcoords(i))
    end do
end subroutine

subroutine dmacrys_reduce_cell(natoms, atomcoords)
    use cell_reduction
    use genrigid
    use vec3
    implicit none
    integer natoms, i
    double precision atomcoords(1:3*NATOMS)
    double precision lattice(3,3)
    double precision transform(3,3)
    double precision lattice_inv(3,3)
    double precision tmp(3,3)
    double precision rcoords(DEGFREEDOMS)

    call get_lattice_matrix(atomcoords(3*natoms-5:3*natoms), lattice)

    call cell_set_lattice(lattice)
    call cell_minimum_reduction

    if(.not.cell_has_changed()) return

    call invert3x3(lattice, lattice_inv)

    tmp = cell_get_transformation()
    call cell_get_lattice(lattice)

    transform = matmul(lattice,matmul(tmp, lattice_inv))
    do i=1,natoms
        atomcoords(3*i-2:3*i) = matmul(transform, atomcoords(3*i-2:3*i))
    end do
    atomcoords(3*natoms-6+1) = lattice(1,1)
    atomcoords(3*natoms-6+2) = lattice(2,1)
    atomcoords(3*natoms-6+3) = lattice(3,1)
    atomcoords(3*natoms-6+4) = lattice(2,2)
    atomcoords(3*natoms-6+5) = lattice(3,2)
    atomcoords(3*natoms-6+6) = lattice(3,3)

    call TRANSFORMCTORIGID (atomcoords, rcoords(1:DEGFREEDOMS))
    do i=1,3*NRIGIDBODY
        rcoords(i) = rcoords(i) - entier(rcoords(i))
    end do
    call TRANSFORMRIGIDTOC (1, NRIGIDBODY, atomcoords, rcoords(1:DEGFREEDOMS))
end subroutine

subroutine transformrigidtoc_wrapper (CMIN, CMAX, XCOORDS, XRIGIDCOORDS)
      
  USE COMMONS, ONLY: NATOMS
  use genrigid, only:degfreedoms,transformrigidtoc
  
  IMPLICIT NONE
  
  INTEGER :: J1, J2, J5, J7, J9, NP        !NP = No of processors
  INTEGER :: CMIN, CMAX
  DOUBLE PRECISION :: P(3), RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
  DOUBLE PRECISION :: XRIGIDCOORDS(DEGFREEDOMS), XCOORDS(3*NATOMS)
  DOUBLE PRECISION :: COM(3) ! center of mass
  LOGICAL          :: GTEST !, ATOMTEST
  DOUBLE PRECISION :: MLATTICE(3,3)
  call transformrigidtoc(cmin, cmax, xcoords, xrigidcoords)
end subroutine

subroutine dmacrys_dump_cif_wrapper(filename, rigid_coords, cifname)
    use dmacrys_interface, only : dmacrys_dump_cif
    use commons, only : natoms
    implicit none
    character*(*) :: filename, cifname
    double precision, intent(in)  :: rigid_coords(1:3*NATOMS) ! that is the natoms hack
	call dmacrys_dump_cif(filename, rigid_coords, cifname)
end subroutine 

subroutine dmacrys_get_masses(masses)
    use commons, only: natoms
    use genrigid, only: gr_weights
    implicit none
    double precision, intent(out) :: masses(natoms - 2)
    masses(:) = gr_weights(1:natoms-2)
end subroutine

subroutine dmacrys_get_atomname(atom, name)
    use iso_c_binding
    use system, only: atoms
    implicit none
    character (kind=c_char, len=64) :: name
    integer, intent(in)  :: atom

    name = atoms(atom)%name
end subroutine

function dmacrys_rigidbody_natoms(rigidbody)
    use genrigid
    integer dmacrys_rigidbody_natoms, rigidbody
    dmacrys_rigidbody_natoms = NSITEPERBODY(rigidbody)
end function

subroutine dmacrys_rigidbody_atoms(rigidbody, atoms, coords)
    use genrigid
    integer, intent(in) :: rigidbody
    integer, intent(out) :: atoms(NSITEPERBODY(rigidbody))
    double precision, intent(out) :: coords(3*NSITEPERBODY(rigidbody))
    integer nsites, i
    nsites = NSITEPERBODY(rigidbody)
    atoms(1:nsites) = rigidgroups(1:nsites, rigidbody)
    do i=1,nsites
        coords(3*i-2:3*i) = SITESRIGIDBODY(i,:,rigidbody)
    end do
end subroutine

