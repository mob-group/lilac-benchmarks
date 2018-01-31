! Leforestier 2002 VRT(MCY-5f)
! x(i,j) = ith dimension of jth atom
! O is 1, H is 2 and 3
! distances in a.u., energies in a.u.

module MCY

implicit none

! parameters in a.u.
real*8, parameter :: Q2 = 0.583599d0,&
							AOO = 2376.94d0,&
							bOO = 2.71770d0,&
							AHH = 0.98642d0,&
							bHH = 1.56704d0,&
							AOH = 2.91318d0,&
							bOH = 1.50467d0,&
							AVH = 0.475543d0,&
							bVH = 1.18216d0,&
							dVC = 0.517424d0,&
							dVS = 0.03857d0

contains

!#####################################

real*8 function potential(coords)
real*8, intent(in) :: coords(:)
integer :: i
real*8 :: x(3,6), R1(3), R2(3), q1, q2, theta, pot
x = reshape(coords,(/3,6/))
potential = inter(x)
do i = 0, 1
	R1 = x(:,2+i*3) - x(:,1+i*3)
	R2 = x(:,3+i*3) - x(:,1+i*3)
	q1 = sqrt(sum(R1**2))
	q2 = sqrt(sum(R2**2))
	theta = acos(dot_product(R1, R2) / (q1*q2))
	call POTS(pot, q1, q2, theta)
	potential = potential + pot
end do
end function potential

!#####################################

real*8 function inter(x)
real*8, intent(in) :: x(3,6)
integer :: i, j
real*8 :: xn(3,2), xv(3,2), r, R1(3), R2(3)

inter = 0.0d0

!calculate position of negative/virtual sites

R1 = x(:,2) + x(:,3) - 2.0d0*x(:,1)
R1 = R1 / sqrt(sum(R1**2))
R2 = x(:,5) + x(:,6) - 2.0d0*x(:,4)
R2 = R2 / sqrt(sum(R2**2))

xn(:,1) = x(:,1) + dVC*R1
xn(:,2) = x(:,4) + dVC*R2
xv(:,1) = x(:,1) + dVS*R1
xv(:,2) = x(:,4) + dVS*R2

!sum over interatomic Coulombic interactions and dispersion terms

!between negative sites
r = sqrt(sum((xn(:,1) - xn(:,2))**2))
inter = inter + 4*Q2/r

!between oxygen atoms
r = sqrt(sum((x(:,1) - x(:,4))**2))
inter = inter + AOO * exp(-bOO*r)

do i = 2, 3 ! hydrogens on molecule 1
	do j = 5, 6 ! hydrogens on molecule 2
		r = sqrt(sum((x(:,i) - x(:,j))**2))
		inter = inter + Q2/r + AHH*exp(-bHH*r)
	end do
end do

do i = 1, 2 ! negative/virtual sites
	do j = 2, 3 ! hydrogens
		r = sqrt(sum((xn(:,i) - x(:,j+6-i*3))**2))
		inter = inter - 2*Q2/r
		r = sqrt(sum((xv(:,i) - x(:,j+6-i*3))**2))
		inter = inter - AVH*exp(-bVH*r)
	end do
end do

do i = 1, 4, 3 ! oxygen atoms
	do j = 2, 3 ! hydrogens
		r = sqrt(sum((x(:,i) - x(:,j+4-i))**2))
		inter = inter + AOH*exp(-bOH*r)
	end do
end do

end function inter

!#####################################

end module MCY
