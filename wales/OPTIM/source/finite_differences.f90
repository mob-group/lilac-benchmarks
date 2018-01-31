module finite_differences

implicit none

contains

!###############################

function FinDifGrad(x, pot, eps, double)
DOUBLE PRECISION, intent(in) :: x(:), eps
LOGICAL, intent(in) :: double ! better accuracy gradient
interface
	DOUBLE PRECISION function pot(x)
	DOUBLE PRECISION, intent(in), dimension(:) :: x
	end function pot
end interface
DOUBLE PRECISION, dimension(size(x)) :: FinDifGrad, ei
integer :: i

ei = 0.0
if (double) then
	do i = 1, size(x)
		ei(i) = eps
		FinDifGrad(i) = (-pot(x+2*ei) + 8*pot(x+ei) - 8*pot(x-ei) + pot(x-2*ei)) / (12*eps)
		ei(i) = 0.0
	end do
else
	do i = 1, size(x)
		ei(i) = eps
		FinDifGrad(i) = (pot(x+ei) - pot(x-ei)) / (2*eps)
		ei(i) = 0.0
	end do
end if
end function FinDifGrad

!###############################

function FinDifHess(x, g, grad, eps) result(hess)
DOUBLE PRECISION, intent(in) :: x(:), g(:), eps
interface
	function grad(x)
	DOUBLE PRECISION, intent(in), dimension(:) :: x
	DOUBLE PRECISION, dimension(size(x)) :: grad
	end function grad
end interface
DOUBLE PRECISION, dimension(size(x)) :: x0
DOUBLE PRECISION, dimension(size(x),size(x)) :: hess
integer :: i

x0 = x
do i = 1, size(x)
	x0(i) = x(i) + eps
	hess(i,:) = (grad(x0) - g) / eps
	x0(i) = x(i)
end do
call symmetrize(hess)
end function FinDifHess

!###############################

function FinDifHess_Pot(x, pot, eps) result(hess)
DOUBLE PRECISION, intent(in) :: x(:), eps
interface
	DOUBLE PRECISION function pot(x)
	DOUBLE PRECISION, intent(in), dimension(:) :: x
	end function pot
end interface
DOUBLE PRECISION :: V0
DOUBLE PRECISION, dimension(size(x)) :: ei, ej
DOUBLE PRECISION, dimension(size(x),size(x)) :: hess
integer :: i, j

V0 = pot(x)

ei = 0.0
ej = 0.0

do i = 1, size(x)
   ei(i) = eps
   hess(i,i) = (-pot(x+2*ei) + 16*pot(x+ei) - 30*V0 + 16*pot(x-ei) - pot(x-2*ei)) / (12*eps**2)
   do j = 1, i-1
      ej(j) = eps
!     hess(i,j) = (pot(x+ei+ej) - pot(x-ei+ej) - pot(x+ei-ej) + pot(x-ei-ej)) / (2*eps)**2
      hess(i,j) = (-pot(x+2*ei+2*ej) + 16*pot(x+ei+ej) - 30*V0 + 16*pot(x-ei-ej) - pot(x-2*ei-2*ej)) / &
  &   (24*eps**2) - (hess(i,i) + hess(j,j)) / 2
      hess(j,i) = hess(i,j)
      ej(j) = 0.0
   end do
   ei(i) = 0.0
end do
end function FinDifHess_Pot

!###############################

subroutine symmetrize(hess)
DOUBLE PRECISION, intent(inout) :: hess(:,:)
integer :: i, j, N
DOUBLE PRECISION :: tmp
N = size(hess,1)
do i = 1, N
	do j = i+1, N
		tmp = 0.5 * (hess(i,j) + hess(j,i))
		hess(i,j) = tmp
		hess(j,i) = tmp
	end do
end do
end subroutine symmetrize

!###############################

end module finite_differences
