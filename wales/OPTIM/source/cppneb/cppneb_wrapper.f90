module cpp_potential
contains
function get_energy_gradient(ncoords, coords, grad, userdata) result(energy) bind(c)
    use iso_c_binding
    implicit none
    integer(c_int), value, intent(in) :: ncoords
    real(c_double), intent(in) :: coords(ncoords)
    real(c_double), intent(out) :: grad(ncoords) 
    type(c_ptr) :: userdata
    real(c_double) :: energy
    real(c_double) :: rms
    real(c_double) :: dummy

    grad = 0.0
    energy = 0.0

    call potential(coords,energy,grad,.true.,.false.,rms,.false.,.false.)
    ! Warning. There seems to be a compiler bug here: the energy is calculated correctly but doesn't get passed
    ! down to the C++ end. I'll investigate stopping this by lowering the level of optimisation, but for the moment
    ! it seems that the following write statement is enough to correct the behaviour.
    write(*,*) "Called potential from C++, energy ", energy
    !                                 !gtest  stest       ptest  boxtest!!
    ! first, but not second derivatives; no printing; is boxtest, ever?

end function
end module cpp_potential

module cpp_distance
use genrigid, only: rigidinit, rb_distance, degfreedoms
contains
function get_distance(ncoords, left, right, gradient_left, gradient_right) result(distance) bind(c)
    use iso_c_binding
    implicit none
    real(c_double) :: distance
    integer(c_int), value, intent(in) :: ncoords
    real(c_double), intent(in) :: left(ncoords), right(ncoords)
    real(c_double), intent(out) :: gradient_left(ncoords), gradient_right(ncoords)

    distance=0

    if (rigidinit) then
        call rb_distance(distance, left, right, gradient_left, gradient_right, .true.)
        gradient_left(degfreedoms+1:) = 0
        gradient_right(degfreedoms+1:) = 0
    else
        call cartesian_distance(ncoords, left, right, gradient_left, gradient_right, distance)
    endif

end function

subroutine cartesian_distance(ncoords, left, right, gradient_left, gradient_right, distance)
    implicit none
    integer :: ncoords
    double precision, intent(in)  :: left(ncoords), right(ncoords)
    double precision, intent(out) :: gradient_left(ncoords), gradient_right(ncoords), distance
    distance=0
    gradient_left=-(right-left)
    gradient_right=(right-left)
    distance=sqrt(sum((right-left)*(right-left)))

end subroutine

end module cpp_distance

!subroutine neb_example(qm,qp,nimages,ncoords,xyz,eee,nitermax)  ! sn402: changed
subroutine neb_example(nimages,ncoords,xyz,eee,nitermax)
    use iso_c_binding
    use neb
    use cpp_distance
    use cpp_potential
    use keyneb, only: dumpnebeos, dumpnebeosfreq, dumpnebxyz, dumpnebxyzfreq, dumpnebpts, dumpnebptsfreq
    use nebdata, only: truegrad
    use keyminimizer, only: rmstol
    use key, only: maxnebbfgs, nebk, kadjusttol, kadjustfrq, kadjustfrac, nebmupdate, nebdguess
    use commons, only: debug
    use nebutils

    implicit none
    integer :: i, i2, j, nimages, ncoords, k
    integer, intent(in) :: nitermax
!    double precision, intent(in)  :: qm(ncoords),qp(ncoords)
    double precision :: energy, rms
!    double precision, intent(out) :: xyz(ncoords*nimages), eee(nimages)
    double precision, intent(out) :: xyz(ncoords*(nimages+2)), eee(nimages+2)
    double precision :: thisimage(ncoords), tempgrad(ncoords)
    integer :: verbosity = 0, iprint = -1

    if(debug) then ! Set the lbfgs parameters to maximum information printed out
        verbosity = 2
        iprint = 5
    endif

    ! set up the band
    call neb_setup(c_funloc(get_energy_gradient), c_funloc(get_distance), c_null_ptr)
    call neb_initialize_path(nimages+2, ncoords) !sn402 changed
!    call neb_initialize_path(nimages, ncoords)

    ! load images into band - including the fixed configurations at the start and end of the band
    do i=0,nimages+1  ! sn402 changed
!     do i=0,nimages-1
!        i2=ncoords*(i+1)  ! sn402 changed
        i2=ncoords*i
        call neb_set_image_coords(i,ncoords,xyz(i2+1:i2+ncoords))
    enddo

    ! begin
!    call neb_start()
    call neb_start_with_lbfgs(rmstol, nebmupdate, nebdguess)
    call neb_parameters(1, rmstol, nebk, kadjusttol, kadjustfrac, maxnebbfgs, 1, iprint, verbosity)
!    call lbfgs_parameters(nebmupdate, nebmaxerise, nebdguess)  ! Need to sort out initialising these.

    ! run neb until converged or nitermax steps
    do i=1,nitermax

        ! If required, print out data for this step: EofS, xyz coordinates or pts files
        if((dumpnebeos .and. mod(i-1,dumpnebeosfreq)==0) .or. &
            & (dumpnebxyz .and. mod(i-1,dumpnebxyzfreq)==0) .or. &
                & (dumpnebpts .and. mod(i-1,dumpnebptsfreq)==0)) then
            ! All three require knowing the current band coordinates
            do j=0,nimages+1
                i2=ncoords*j
                call neb_get_image_coords(j,ncoords,xyz(i2+1:i2+ncoords))
            enddo

            if(dumpnebeos .and. mod(i-1,dumpnebeosfreq)==0) then
                call neb_get_image_energies(nimages+2, eee)
                ! Compute the distances between images
                call imagedistribution
                ! Print the actual data
                call writeprofile(i-1)
            endif
            if(dumpnebxyz .and. mod(i-1,dumpnebxyzfreq)==0) call rwg("w",.false.,i-1)
            if(dumpnebpts .and. mod(i-1,dumpnebptsfreq)==0) call savebandcoord
        endif

        if (neb_step()) exit
        if (mod(i,kadjustfrq)==0) then
            call neb_adjust_k()
        endif

    enddo

    ! write band and energies to xyz and eee
     do i=0,nimages+1  !sn402 changed
!    do i=0,nimages-1
        i2=ncoords*i
!        i2=ncoords*(i+1) !sn402 changed
        thisimage(:) = xyz(i2+1:i2+ncoords)

        call neb_get_image_coords(i,ncoords,thisimage)
        call potential(thisimage,energy,tempgrad,.true.,.false.,rms,.false.,.false.)

        ! sn402: save the gradients for use in bfgsts.
        truegrad(i2+1:i2+ncoords) = tempgrad(:)
        eee(i+1)=energy
         enddo
    call imagedistribution

    call neb_cleanup()

end subroutine
