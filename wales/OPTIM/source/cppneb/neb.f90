!        typedef double potential_callback_t(int ncoords, double
!        *coords, double
!*grad, void
!                *userdata);
!
!
!        void neb_set_image_coords(int image, int ncoords, const double
!        *coords);
!        void neb_get_image_coords(int image, int ncoords, double
!        *coords);
!

module neb
   use iso_c_binding
   interface
       ! void neb_setup(potential_callback_t *potential, void
       ! *userdata);
       subroutine neb_setup(potential_callback, distance_callback, userdata) bind(c)
           use iso_c_binding
           implicit none
           type(c_funptr), value, intent(in) :: potential_callback 
           type(c_funptr), value, intent(in) :: distance_callback 
           type(c_ptr), intent(in) :: userdata
       end subroutine
 
       ! void neb_parameters(double tol, double maxstep, int maxiter, int iprint, int
       ! verbosity);
       subroutine neb_parameters(double_nudging, tol, k_initial, adjust_k_tol, adjust_k_factor, &
            & maxstep, maxiter, iprint, verbosity) bind(c)
           use iso_c_binding
           implicit none
           integer(c_int), value, intent(in) :: double_nudging
           real(c_double), value, intent(in) :: tol
           real(c_double), value, intent(in) :: k_initial
           real(c_double), value, intent(in) :: adjust_k_tol
           real(c_double), value, intent(in) :: adjust_k_factor
           real(c_double), value, intent(in) :: maxstep
           integer(c_int), value, intent(in) :: maxiter
           integer(c_int), value, intent(in) :: iprint
           integer(c_int), value, intent(in) :: verbosity
       end subroutine

       ! void neb_cleanup();
       subroutine neb_cleanup() bind(C)
       end subroutine 
    
       ! void neb_initialize_path(int nimages, int
       ! num_coords_per_image);
       subroutine neb_initialize_path(nimages, num_coords_per_image) bind(c)
           use iso_c_binding
           integer(c_int), value, intent(in) :: nimages
           integer(c_int), value, intent(in) :: num_coords_per_image
       end subroutine

       ! void neb_start();
       subroutine neb_start() bind(c)
       end subroutine

       ! sn402: this should replace the previous subroutine.
       ! void neb_start_with_lbfgs(double rmstol, int setM, double H0);
       subroutine neb_start_with_lbfgs(rmstol, setM, H0) bind(c)
           use iso_c_binding
           implicit none
           real(c_double), value, intent(in) :: rmstol
           integer(c_int), value, intent(in) :: setM
           real(c_double), value, intent(in) :: H0
       end subroutine
       
       ! bool neb_step();
       function neb_step() bind(c)
           use iso_c_binding
           logical(c_bool) :: neb_step
       end function

       ! void neb_adjust_k();
       subroutine neb_adjust_k() bind(c)
       end subroutine

       subroutine neb_set_image_coords(image, ncoords, coords) bind(c)
           use iso_c_binding
           integer(c_int), value, intent(in) :: image
           integer(c_int), value, intent(in) :: ncoords
           real(c_double), intent(in) :: coords(ncoords)
       end subroutine

       subroutine neb_get_image_coords(image, ncoords, coords) bind(c)
           use iso_c_binding
           integer(c_int), value, intent(in) :: image
           integer(c_int), value, intent(in) :: ncoords
           real(c_double), intent(in) :: coords(ncoords)
       end subroutine

       subroutine neb_get_image_energies(nimages, energies) bind(c)
           use iso_c_binding
           integer(c_int), value, intent(in) :: nimages
           real(c_double), intent(in) :: energies(nimages)
       end subroutine

       ! bool neb_step();
       function neb_rms() bind(c)
           use iso_c_binding
           real(c_double) :: neb_rms
       end function

    end interface
end module neb
