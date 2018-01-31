module bowmanwater

implicit none

contains

subroutine bowmaninit(nw, pes, dname)
use pes_shell, only : iwaterfcn, pes_init
integer, intent(in) :: nw, pes
character (len=*), intent(in) :: dname
call pes_init(nw, dname)
iwaterfcn = pes
end subroutine bowmaninit

double precision function bowmanpot(x)
use pes_shell, only : f, auang
double precision, dimension(:), intent(in) :: x
real, dimension(:,:), allocatable :: xx
integer :: natom, i
natom = size(x) / 3
allocate(xx(3,natom))
forall(i=1:natom) xx(:,i) = x(3*i-2:3*i) / auang
bowmanpot = f(xx)
deallocate(xx)
end function bowmanpot

end module  bowmanwater
