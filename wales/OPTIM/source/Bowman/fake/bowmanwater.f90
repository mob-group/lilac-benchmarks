module bowmanwater

implicit none

contains

subroutine bowmaninit(nw, pes, dname)
integer, intent(in) :: nw, pes
character (len=*), intent(in) :: dname
print *, "must compile with ifort in order to use Bowmans water potential"
stop
end subroutine bowmaninit

double precision function bowmanpot(x)
double precision, dimension(:), intent(in) :: x
print *, "must compile with ifort in order to use Bowmans water potential"
bowmanpot=0.0D0 ! to avoid compiler warning for unassigned function value.
stop
end function bowmanpot

end module bowmanwater
