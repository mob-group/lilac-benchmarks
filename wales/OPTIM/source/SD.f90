module SDwater

implicit none

real*8 :: alpha = 1.444 !Angstroms**3 !polarizability of O. (H+ can't be polarized at all)
real*8 :: re = 0.9584 !Angstroms
integer :: NO, Na
integer :: lwork
real*8, dimension(:), allocatable :: work, q

interface dist
   module procedure dist_x, dist_d
end interface

contains

!#######################################################################

!vector between ith particle and jth particle
pure function diff(x, j, i)
real*8, dimension(3) :: diff
real*8,  intent(IN), dimension(:) :: x
integer, intent(IN)               :: i, j
diff = x(i*3-2:i*3) - x(j*3-2:j*3)
end function diff

!distance between ith particle and jth particle in Angstrom
pure real*8 function dist_x(x, i, j)
real*8,  intent(IN), dimension(:) :: x
integer, intent(IN)               :: i, j
dist_x = sqrt(sum(diff(x, i, j)**2))
end function dist_x

pure real*8 function dist_d(diff)
real*8,  intent(IN), dimension(3) :: diff
dist_d = sqrt(sum(diff**2))
end function dist_d

!#######################################################################
! FUNCTIONS
!#######################################################################

pure real*8 function phiHH(r)
real*8, intent(in) :: r ! Angstrom
phiHH = 332.1669 / r
end function phiHH

!#######################################################################

pure real*8 function phiOH(r)
real*8, intent(in) :: r
phiOH = 332.1669 / r * (10*exp(-3.699392820*r)-2) &
        + (-184.6966743*(r-re) + 123.9762188*(r-re)**2) * exp(-8*(r-re)**2)
end function phiOH

!#######################################################################

pure real*8 function phiOO(r)
real*8, intent(in) :: r
phiOO = 1328.6676/r               &
        + 24/(1+exp(2.5*(r-2.9))) &
        + 90/(1+exp(8*(r-2.45)))  &
        + exp(-6*(r-2.7))
end function phiOO

!#######################################################################

pure real*8 function omKr3(r) ! (1 - K) / r**3
real*8, intent(in) :: r
omKr3 = 1 / ( r**3                                        &
               + 1.855785223*(r-re)**2 * exp(-8*(r-re)**2) &
               + 16.95145727 * exp(-2.702563425*r) )
end function omKr3

!#######################################################################

pure real*8 function omLr3(r) ! (1 - L) / r**3
real*8, intent(in) :: r
omLr3 = 1 - exp(-3.169888166*r) * (1 + (3.169888166 + (5.024095492 + (-17.99599078 + 23.92285*r)*r)*r)*r )
omLr3 = omLr3 / r**3
end function omLr3

!#######################################################################
! FUNCTION DERIVATIVES
!#######################################################################

pure real*8 function dphiHHdror(r) ! minus derivative of phiHH w.r.t r divided by r
real*8, intent(in) :: r ! Angstrom
dphiHHdror = 332.1669 / r**3
end function dphiHHdror

!#######################################################################

pure real*8 function dphiOHdror(r) ! minus derivative of phiOH w.r.t r divided by r
real*8, intent(in) :: r
real*8 :: tmp
tmp = 3.699392820*r
dphiOHdror = 332.1669 / r**2 * (10*exp(-tmp)*(1+tmp) - 2)
tmp = r - re
dphiOHdror = dphiOHdror + exp(-8*tmp**2) * (184.6966743 + tmp*(-2*123.9762188 + 16*tmp*(-184.6966743 + 123.9762188*tmp)))
dphiOHdror = dphiOHdror / r
end function dphiOHdror

!#######################################################################

pure real*8 function dphiOOdror(r) ! minus derivative of phiOO w.r.t r divided by r
real*8, intent(in) :: r
real*8 :: tmp
dphiOOdror = 1328.6676 / r**2
tmp = exp(2.5*(r-2.9))
dphiOOdror = dphiOOdror + 24 / (1+tmp)**2 * 2.5*tmp
tmp = exp(8*(r-2.45))
dphiOOdror = dphiOOdror + 90 / (1+tmp)**2 * 8 * tmp
dphiOOdror = dphiOOdror + 6 * exp(-6*(r-2.7))
dphiOOdror = dphiOOdror / r
end function dphiOOdror

!#######################################################################

pure real*8 function domKr3dr(r) ! d/dr (1 - K)/r**3
real*8, intent(in) :: r
real*8 :: rre2
rre2 = (r - re)**2
domKr3dr = r**3 + 1.855785223*rre2 * exp(-8*rre2) + 16.95145727 * exp(-2.702563425*r)
domKr3dr = - ( 3*r**2 + 1.855785223*(r-re)*(2 - 16.95145727*rre2)*exp(-8*rre2) - 2.702563425 * 16.95145727 * exp(-2.702563425*r) ) &
  &      / domKr3dr**2
end function domKr3dr

!#######################################################################

pure real*8 function domLr3dr(r) ! d/dr (1 - L)/r**3
real*8, intent(in) :: r
domLr3dr = r**2
domLr3dr = - (3 - exp(-3.169888166*r) * (3 + (2*3.169888166 + (5.024095492 - 23.92285*domLr3dr)*r)*r +  &
  &           3.169888166*(1 + (3.169888166 + (5.024095492 + (-17.99599078 + 23.92285*r)*r)*r)*r)*r ) ) / domLr3dr**2
end function domLr3dr

!#######################################################################
! DIPOLE
!#######################################################################

! calculates dipole mu(i,:) on the ith Oxygen atom

! mu = alpha * G
! where alpha is polarizability and G is modified field
! G = - sum( |r> * omK * q / r**3) - sum( <T|mu> * omK / r**3)
! T = 1 - 3 |r><r| / r**2
! phi = 0.5 * sum( <mu|r> * omL * q / r**3)

! linear coupled equations: mu=alpha*G, G=b+C.mu
! G = b + alpha*C.G
! I.G = b + D.G
! (I-D).G = b
! A.G = b !solve this with ?getrs

!######################################################################

function field1(x) ! at Oxygens = - sum( |r> * omK * q / r**3 )
real*8, intent(in) :: x(3*Na)
real*8 :: field1(3*NO)
real*8 :: rv(3), r
integer :: i, j
field1 = 0
do i = 1, NO !oxygen atoms
	do j = 1, Na
      if (j==i) cycle
      rv = diff(x,i,j)
      r = dist(rv)
		field1(3*i-2:3*i) = field1(3*i-2:3*i) - rv*q(j)*omKr3(r)
   end do
end do
end function field1

!######################################################################

function field2(x) ! at Oxygens = - sum( <T| * omK * q / r**3 )
real*8, intent(in) :: x(3*Na)
real*8 :: field2(3*NO,3*NO)
real*8 :: rv(3), r
integer :: i, j
field2 = 0
do i = 1, NO !oxygen atoms
	do j = 1, NO
      if (j==i) cycle
		rv = diff(x, i, j)
		r = dist(rv)
      field2(3*i-2:3*i,3*j-2:3*j) = field2(3*i-2:3*i,3*j-2:3*j) - T(r, rv)*omKr3(r)
   end do
end do
end function field2

!######################################################################

!field from an electric dipole
!pot(r) = mu.\hat{r} / r**2 / 4pie_0
!where \hat{r} is unit vector in direction of r
!field = (3*(mu.\hat{r})*\hat{r} - mu)/r**3 / 4pie_0
!T_ij = \delta_ij - 3*r_i*r_j
pure function T(r, rv)
real*8, dimension(3,3) :: T
real*8, intent(in) :: r, rv(3)
real*8 :: tmp
integer :: k, l
tmp = - 3 / r**2
do k = 1, 3
   do l = k, 3
      T(k,l) = tmp * rv(k) * rv(l)
   end do
end do
do k = 2, 3
   do l = 1, k-1
      T(k,l) = T(l,k)
   end do
end do
forall(k=1:3) T(k,k) = T(k,k) + 1
end function T

!#######################################################################

subroutine factorize(A, ipiv)
real*8, intent(inout) :: A(3*NO,3*NO)
integer, intent(out) :: ipiv(3*NO)
integer :: info

call DSYTRF('U', 3*NO, A, 3*NO, ipiv, work, lwork, info)

if (int(work(1)) /= lwork) print *, work(1), "is better than", lwork
if (info < 0) then
   print '("the ", I0, "th parameter had an illegal value")', -info
   stop
else if (info > 0) then
   print '("error: U(",I0,",",I0,") = 0")', info, info
   stop
end if
end subroutine factorize

!#######################################################################

subroutine solve(A, ipiv, mu)
real*8, intent(in) :: A(3*NO,3*NO)
integer, intent(in) :: ipiv(3*NO)
real*8, intent(inout) :: mu(NO,3)
integer :: info

call DSYTRS('U', 3*NO, 1, A, 3*NO, ipiv, mu, 3*NO, info)
if (info < 0) then
   print '("the ", I0, "th parameter had an illegal value")', -info
   stop
else if (info > 0) then
   print '("error: U(",I0,",",I0,") = 0")', info, info
   stop
end if
end subroutine solve

!#######################################################################
! DIPOLE DERIVATIVES
!######################################################################

function dTdx(r, rv, ndim, s)
real*8, dimension(3,3) :: dTdx
real*8, intent(in) :: r, rv(3)
integer, intent(in) :: ndim, s
real*8 :: tmp
integer :: k, l
forall(k=1:3,l=1:3) dTdx(k,l) = 2 * rv(k) * rv(l) * rv(ndim)
tmp = r**2
forall(l=1:3) dTdx(ndim,l) = dTdx(ndim,l) - tmp*rv(l)
forall(k=1:3) dTdx(k,ndim) = dTdx(k,ndim) - tmp*rv(k)
dTdx = dTdx * s * 3 / tmp**2
end function dTdx

!#######################################################################

function dbdx(x, n, mu)
real*8, intent(in) :: x(3*Na), mu(3*NO)
integer, intent(in) :: n
real*8 :: dbdx(3*NO)
integer :: i, j
real*8 :: r, rv(3), tmp
integer :: nat, ndim
nat = (n-1)/3 + 1
ndim = n - 3*(nat-1)

if (nat > NO) then ! .: nat != i
	do i = 1, NO
		rv = diff(x, i, nat)
		r = dist(rv)
		dbdx(3*i-2:3*i) = q(nat) * domKr3dr(r) * rv(ndim) / r * rv
		dbdx(3*(i-1)+ndim) = dbdx(3*(i-1)+ndim) + omKr3(r) * q(nat)
	end do
else ! nat <= NO
	! i != nat
	do i = 1, NO
		if (i==nat) cycle
		rv = diff(x, i, nat)
		r = dist(rv)
		dbdx(3*i-2:3*i) = domKr3dr(r)*rv(ndim)/r * (rv*q(nat) + matmul(T(r,rv),mu(3*nat-2:3*nat)))

		tmp = omKr3(r)
		dbdx(3*(i-1)+ndim) = dbdx(3*(i-1)+ndim) + tmp * q(nat)
		dbdx(3*i-2:3*i) = dbdx(3*i-2:3*i) + tmp * matmul(dTdx(r,rv,ndim,1), mu(3*nat-2:3*nat))
	end do

	! i = nat
	dbdx(3*nat-2:3*nat) = 0
	do j = 1, Na
		if (j==nat) cycle
		rv = diff(x, nat, j)
		r = dist(rv)
                dbdx(3*nat-2:3*nat) = dbdx(3*nat-2:3*nat) - domKr3dr(r) * rv(ndim)/r * q(j) * rv
		dbdx(n) = dbdx(n) - omKr3(r) * q(j)
	end do

	do j = 1, NO
		if (j==nat) cycle
		rv = diff(x, nat, j)
		r = dist(rv)
		dbdx(3*nat-2:3*nat) = dbdx(3*nat-2:3*nat) - domKr3dr(r) * rv(ndim)/r * matmul(T(r,rv),mu(3*j-2:3*j))
		dbdx(3*nat-2:3*nat) = dbdx(3*nat-2:3*nat) + omKr3(r) * matmul(dTdx(r,rv,ndim,-1),mu(3*j-2:3*j))
	end do
end if
dbdx = - alpha * dbdx
end function dbdx

!#######################################################################
! MAIN SUBROUTINES
!#######################################################################

subroutine sdinit(NO_, Np_)
integer, intent(in) :: NO_, Np_
integer :: info
real*8, dimension(3*NO_, 3*NO_) :: A
real*8, dimension(3*NO_) :: b
integer, dimension(3*NO_) :: ipiv

NO = NO_
Na = 3*NO_ + Np_
allocate(q(Na))
q(1:NO) = -2 * sqrt(332.1669)
q(NO+1:) = sqrt(332.1669)

allocate(work(1))
call DSYSV('U', size(A,dim=2), 1, A, size(A,dim=1), ipiv, b, size(b), work, -1, info)
if (info < 0) then
   print '("the ", I0, "th parameter had an illegal value")', -info
   stop
else if (info > 0) then
   print '("error: U(",I0,",",I0,") = 0")', info, info
   stop
end if
lwork = work(1)
!print *, 'best lwork =', lwork
deallocate(work)
allocate(work(lwork))
end subroutine sdinit

!#######################################################################

subroutine destruct()
deallocate(work, q)
end subroutine destruct

!#######################################################################

real*8 function sdpotential(x) result(V)
real*8, intent(in) :: x(:)
integer :: i, j
real*8 :: mu(3*NO), rv(3), r, phi2
real*8, dimension(3*NO,3*NO) :: A
integer :: ipiv(3*NO)

!sum over HH pairs
V = 0
do j = NO+2, Na
   do i = NO+1, j-1
      r = dist(x,i,j)
      V = V + phiHH(r)
   end do
end do

!sum over OH pairs
do j = 1, NO !O atoms
   do i = NO+1, Na !H atoms
      r = dist(x,i,j)
      V = V + phiOH(r)
   end do
end do 

!sum over OO pairs
do j = 2, NO
   do i = 1, j-1
      r = dist(x, i, j)
      V = V + phiOO(r)
   end do
end do

mu = field1(x) * alpha

A = field2(x)
A = A * alpha ! is D
A = - A
forall(i=1:3*NO) A(i,i) = 1 + A(i,i)

call factorize(A, ipiv)
call solve(A, ipiv, mu)

phi2 = 0
do i = 1, Na !all atoms
   do j = 1, NO !Oxygen atoms
      if (i==j) cycle
      rv = diff(x,j,i)
      r = dist(rv)
      phi2 = phi2 + dot_product(mu(3*j-2:3*j), rv) * q(i) * omLr3(r)
   end do
end do
phi2 = phi2 / 2
V = V + phi2

end function sdpotential

!######################################################################

function sdgrad(x) result(g)
real*8, intent(in) :: x(:)
real*8 :: g(size(x))
integer :: i, j, nat, ndim, n
real*8 :: mu(3*NO), rv(3), r, dmudx(3*NO)
real*8, dimension(3*NO,3*NO) :: A
integer :: ipiv(3*NO)

g = 0

! HH pairs
do j = NO+2, Na
   do i = NO+1, j-1
		rv = diff(x,i,j)
		r = dist(rv)
		rv = rv * dphiHHdror(r)
		g(3*j-2:3*j) = g(3*j-2:3*j) - rv
		g(3*i-2:3*i) = g(3*i-2:3*i) + rv
   end do
end do

! OH pairs
do j = 1, NO !O atoms
   do i = NO+1, Na !H atoms
  		rv = diff(x,i,j)
		r = dist(rv)
		rv = rv * dphiOHdror(r)
		g(3*j-2:3*j) = g(3*j-2:3*j) - rv
		g(3*i-2:3*i) = g(3*i-2:3*i) + rv
   end do
end do 

! OO pairs
do j = 2, NO
   do i = 1, j-1
  		rv = diff(x,i,j)
		r = dist(rv)
		rv = rv * dphiOOdror(r)
		g(3*j-2:3*j) = g(3*j-2:3*j) - rv
		g(3*i-2:3*i) = g(3*i-2:3*i) + rv
   end do
end do

mu = field1(x) * alpha

A = field2(x)
A = A * alpha ! is D
A = - A
forall(i=1:3*NO) A(i,i) = 1 + A(i,i)

call factorize(A, ipiv)
call solve(A, ipiv, mu)

do n = 1, 3*Na
	nat = (n-1)/3 + 1
	ndim = n - 3*(nat-1)
	dmudx = dbdx(x, n, mu)
	call solve(A, ipiv, dmudx)
do i = 1, Na
		do j = 1, NO
			if (i==j) cycle
			rv = diff(x,j,i)
			r = dist(rv)
			g(n) = g(n) + q(i) * omLr3(r) * dot_product(dmudx(3*j-2:3*j),rv) / 2
		end do
	end do

	do j = 1, NO
		if (j==nat) cycle
		rv = diff(x,j,nat)
		r = dist(rv)
		g(n) = g(n) + q(nat) / 2 * (mu(3*(j-1)+ndim) * omLr3(r) + dot_product(mu(3*j-2:3*j),rv) * domLr3dr(r) * rv(ndim)/r)
	end do

	if (nat <= NO) then
		do i = 1, Na
			if (i==nat) cycle
			rv = diff(x,nat,i)
			r = dist(rv)
			g(n) = g(n) - q(i) / 2 * (mu(n) * omLr3(r) + domLr3dr(r) * rv(ndim)/r * dot_product(mu(3*nat-2:3*nat),rv))
		end do
	end if
end do
end function sdgrad

function sdhess(x, g) result(hess)
real*8, intent(in) :: x(:), g(:)
real*8 :: hess(size(x),size(x)), x0(size(x))
integer :: i, j
real*8 :: eps = 1.0D-8, tmp

x0 = x
do i = 1, 3*Na
	x0(i) = x(i) + eps
	hess(i,:) = (sdgrad(x0) - g) / eps
	x0(i) = x(i)
end do
! symmetrize
do i = 1, 3*Na
	do j = 1, 3*Na
		tmp = 0.5 * (hess(i,j) + hess(j,i))
		hess(i,j) = tmp
		hess(j,i) = tmp
	end do
end do
end function sdhess

!######################################################################

end module sdwater
