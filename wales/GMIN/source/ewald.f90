module ewald
use commons, only: natoms, stchrg, ortho, boxderivt, box_params, box_paramsgrad, &
&                  ewaldalpha, ewaldrealc, ewaldrecipc

implicit none

contains

! -----------------------------------------------------------------------------------
! dj337

! COMPUTES ENERGY AND GRADIENT OF POTENTIALS USING EWALD SUMMATION.
! Usable for any potential that satifisfies the equation:
! U_n = (1/2)*sum_L(sum_i,j(B_ij/(rij+L)**n))
! where n is any integer and L are lattice vectors.
! A separate subroutine is used to calculate the special case for the
! Coulomb potential (when n=1).

! All equations for energy and gradient of Coulomb summation follow from:
! Karasawa, N. and Goddard III, W. A. J. Phys. Chem., 93, 7320-7327 (1989).
 
! All input / output are in absolute Cartesian coordinates.

! Assuming all units for length, charge, and energy are in atomic units.

! Works for either orthorhombic or triclinic unit cells. Computes energy gradient wrt
! cell parameters when BOXDERIVT keyword is true.
! -----------------------------------------------------------------------------------
      subroutine ewaldsum(n, x, g, etot, gtest)

      use cartdist, only: get_reciplatvec, build_H

      implicit none

      integer, intent(in)           :: n
      integer                       :: newaldreal(3), newaldrecip(3)
      logical, intent(in)           :: gtest
      double precision, intent(in)  :: x(3*natoms)
      double precision, intent(out) :: etot, g(3*natoms)
      double precision              :: H(3,3), H_grad(3,3,6)
      double precision              :: reciplatvec(3,3), reciplatvec_grad(3,3,6)
      double precision, parameter   :: pi = 3.141592654d0

      etot = 0.0d0
      g(:) = 0.0d0

      if (n > 1) then
         ! TODO: implement general Ewald summation
         print *, 'Ewald summation not yet implemented for n > 1!'
         return
      else
         ! orthorhombic unit cell
         if (ortho) then
            ! determine number of lattice vectors to sum over
            newaldreal(:) = floor(ewaldrealc/box_params(1:3) + 0.5d0)
            ! compute real-space contribution to energy
            call coulombreal_ortho(x, newaldreal, etot)

            ! determine number of reciprocal lattice vectors to sum over
            newaldrecip(:) = floor(ewaldrecipc*box_params(1:3)/(2.0d0*pi))
            ! compute reciprocal-space contribution to energy
            call coulombrecip_ortho(x, newaldrecip, etot)

            if (gtest) then
               ! compute real-space contribution to gradient
               call coulombrealgrad_ortho(x, newaldreal, g)

               ! compute reciprocal-space contribution to gradient
               call coulombrecipgrad_ortho(x, newaldrecip, g)
            endif
         ! triclinic unit cell
         else
            ! get reciprocal lattice vectors
            call get_reciplatvec(reciplatvec, reciplatvec_grad, .false.)
            ! determine number of lattice vectors to sum over
            newaldreal(1) = floor(ewaldrealc*dsqrt(sum(reciplatvec(1,:)**2))/(2.0d0*pi) + 0.5d0)
            newaldreal(2) = floor(ewaldrealc*dsqrt(sum(reciplatvec(2,:)**2))/(2.0d0*pi) + 0.5d0)
            newaldreal(3) = floor(ewaldrealc*dsqrt(sum(reciplatvec(3,:)**2))/(2.0d0*pi) + 0.5d0)
            ! compute real-space contribution to energy
            call coulombreal_tri(x, newaldreal, etot)

            ! get lattice vectors
            call build_H(H, H_grad, .false.)
            ! determine number of reciprocal lattice vectors to sum over
            newaldrecip(1) = floor(ewaldrecipc*dsqrt(sum(H(1,:)**2))/(2.0d0*pi))
            newaldrecip(2) = floor(ewaldrecipc*dsqrt(sum(H(2,:)**2))/(2.0d0*pi))
            newaldrecip(3) = floor(ewaldrecipc*dsqrt(sum(H(3,:)**2))/(2.0d0*pi))
            ! compute reciprocal-space contribution to energy
            call coulombrecip_tri(x, newaldrecip, etot)

            if (gtest) then
               ! compute real-space contribution to gradient
               call coulombrealgrad_tri(x, newaldreal, g)

               ! compute reciprocal-space contribution to gradient
               call coulombrecipgrad_tri(x, newaldrecip, g)
            endif
         endif ! ortho or triclinic
      endif ! n < 1

      return
      end subroutine ewaldsum

! -----------------------------------------------------------------------------------
! Calculates short-range contribution to Coulomb sum energy. Also includes the self-
! correction term and subtracts within-rigidbody interactions, if needed.

! Assumes orthorhombic unit cell.
! -----------------------------------------------------------------------------------
      subroutine coulombreal_ortho(x, newaldreal, ereal)

      use genrigid, only: rigidinit, nrigidbody, nsiteperbody

      implicit none

      integer                         :: j1, j3, j2, j4, l, m, n, i
      integer, intent(in)             :: newaldreal(3)
      double precision, intent(in)    :: x(3*natoms)
      double precision                :: rmin(3), r(3)
      double precision                :: q1, q2, sumq2, dist, dist2, ewaldrealc2
      double precision                :: vshift, esum, eself, ewrb
      double precision, intent(inout) :: ereal
      double precision, parameter     :: pi = 3.141592654D0

      ! real-space cutoff
      ewaldrealc2 = ewaldrealc**2
      esum = 0.0d0

      ! compute real-space sum
      ! U_real-space = sum_L,i>j(Qij*erfc(alpha*rij)/rij)
      ! iterate over atoms j
      do j1 = 1, natoms
         j3 = 3*j1
         q1 = stchrg(j1)

         ! iterate over atoms i > j
         do j2 = j1+1, natoms
            j4 = 3*j2
            q2 = stchrg(j2)

            ! get distance between atoms
            rmin(1) = x(j3-2)-x(j4-2)
            rmin(2) = x(j3-1)-x(j4-1)
            rmin(3) = x(j3)-x(j4)
            ! minimum image convention
            rmin(1) = rmin(1)-box_params(1)*anint(rmin(1)/box_params(1))
            rmin(2) = rmin(2)-box_params(2)*anint(rmin(2)/box_params(2))
            rmin(3) = rmin(3)-box_params(3)*anint(rmin(3)/box_params(3))

            ! calculate vertical shift
            vshift = q1*q2*erfc(ewaldalpha*ewaldrealc)/ewaldrealc

            ! iterate over boxes
            do l = -newaldreal(1), newaldreal(1)
               r(1) = rmin(1)+box_params(1)*l
               do m = -newaldreal(2), newaldreal(2)
                  r(2) = rmin(2)+box_params(2)*m
                  do n = -newaldreal(3), newaldreal(3)
                     r(3) = rmin(3)+box_params(3)*n
                     dist2 = r(1)**2 + r(2)**2 + r(3)**2
                     if (dist2 < ewaldrealc2) then
                        dist = dsqrt(dist2)
                        ! calculate short-range contribution
                        ! note: don't need factor of 1/2 bc summing over j,i>j
                        esum = esum + q1*q2*erfc(ewaldalpha*dist)/dist - vshift
                     endif ! within cutoff
                  enddo ! n
               enddo ! m
            enddo ! l
         enddo ! atoms j
      enddo ! atoms i

      ! include contribution due to interaction of j1 with periodic images of itself
      ! (separated due to efficiency)
      ! U_periodic-self = 0.5*sum_L(erfc(alpha*rL)/rL)*sum_i(Qi**2)
      sumq2 = 0.0d0
      do j1 = 1, natoms
        q1 = stchrg(j1)
        sumq2 = sumq2 + q1*q1
      enddo

      ! calculate vertical shift
      vshift = erfc(ewaldalpha*ewaldrealc)/(2*ewaldrealc)

      eself = 0.0d0
      ! iterate over boxes
      do l = -newaldreal(1), newaldreal(1)
         r(1) = box_params(1)*l
         do m = -newaldreal(2), newaldreal(2)
            r(2) = box_params(2)*m
            do n = -newaldreal(3), newaldreal(3)
               r(3) = box_params(3)*n
               ! check not in central box
               if (.not.(l.eq.0.and.m.eq.0.and.n.eq.0)) then
                  dist2 = r(1)**2 + r(2)**2 + r(3)**2
                  if (dist2 < ewaldrealc2) then
                     dist = dsqrt(dist2)
                     ! calculate short-range contribution
                     ! note: need factor of 1/2 to prevent double-counting
                     eself = eself + erfc(ewaldalpha*dist)/(2.0d0*dist) - vshift
                  endif ! within cutoff
               endif ! not in central box
            enddo ! n
         enddo ! m
      enddo ! l

      esum = esum + sumq2*eself

      ! compensate for within-rigidbody interactions
      ! calculate within-rigidbody energy using exact Coulomb sum
      ! U_wrb = sum_J(sum_i>j(Qij/rij))
      ! note: don't need factor of 1/2 because summing over i > j
      if (rigidinit) then
         ewrb = 0.0d0
         ! iterate over rigidbodies
         do i = 1, nrigidbody
   
            ! iterate over atoms i
            do j1 = 1, nsiteperbody(i)
               j3 = 3*j1
               q1 = stchrg(j1)
   
               ! iterate over atoms i > j
               do j2 = j1+1, nsiteperbody(i)
                  j4 = 3*j2
                  q2 = stchrg(j2)
   
                  ! calculate rij
                  r(1) = x(j3-2)-x(j4-2)
                  r(2) = x(j3-1)-x(j4-1)
                  r(3) = x(j3)-x(j4)
                  dist2 = r(1)**2 + r(2)**2 + r(3)**2
                  dist = dsqrt(dist2)
   
                  ! calculate within-rigidbody contribution
                  ewrb = ewrb + q1*q2/dist
               enddo ! sites j
            enddo ! sites i
         enddo ! rigid bodies
   
         ! subtract U_wrb
         esum = esum - ewrb
      endif ! rigidinit

      ! compensate for contribution due to self-interaction
      ! U_self-interaction = -alpha*sum_i(Qi**2)/sqrt(pi)
      esum = esum - sumq2*ewaldalpha/dsqrt(pi)

      ereal = ereal + esum

      return
      end subroutine coulombreal_ortho

! -----------------------------------------------------------------------------------
! Calculates short-range contribution to Coulomb sum energy. Also includes the self-
! correction term and subtracts within-rigidbody interactions, if needed.

! Assumes triclinic unit cell.
! -----------------------------------------------------------------------------------
      subroutine coulombreal_tri(x, newaldreal, ereal)

      use genrigid, only: rigidinit, nrigidbody, nsiteperbody, inversematrix
      use cartdist, only: build_H

      implicit none

      integer                         :: j1, j3, j2, j4, l, m, n, i
      integer, intent(in)             :: newaldreal(3)
      double precision, intent(in)    :: x(3*natoms)
      double precision                :: rr(3), rrfracmin(3), rfrac(3), r(3)
      double precision                :: q1, q2, sumq2, dist, dist2, ewaldrealc2
      double precision                :: vshift, esum, eself, ewrb
      double precision                :: H(3,3), H_grad(3,3,6), H_inverse(3,3)
      double precision, intent(inout) :: ereal
      double precision, parameter     :: pi = 3.141592654D0

      ! real-space cutoff
      ewaldrealc2 = ewaldrealc**2
      esum = 0.0d0

      ! get H matrix and inverse
      call build_H(H, H_grad, .false.)
      call inversematrix(H, H_inverse)

      ! compute real-space sum
      ! U_real-space = sum_L,i>j(Qij*erfc(alpha*rij)/rij)
      ! iterate over atoms j
      do j1 = 1, natoms
         j3 = 3*j1
         q1 = stchrg(j1)

         ! iterate over atoms i > j
         do j2 = j1+1, natoms
            j4 = 3*j2
            q2 = stchrg(j2)

            ! get distance between atoms
            rr(:) = x(j3-2:j3) - x(j4-2:j4)
            ! convert to fractional coordinates
            rrfracmin(:) = matmul(H_inverse, rr(:))
            ! minimum image convention
            rrfracmin(1) = rrfracmin(1) - anint(rrfracmin(1))
            rrfracmin(2) = rrfracmin(2) - anint(rrfracmin(2))
            rrfracmin(3) = rrfracmin(3) - anint(rrfracmin(3))

            ! calculate vertical shift
            vshift = q1*q2*erfc(ewaldalpha*ewaldrealc)/ewaldrealc

            ! iterate over boxes
            do l = -newaldreal(1), newaldreal(1)
               rfrac(1) = rrfracmin(1) + l
               do m = -newaldreal(2), newaldreal(2)
                  rfrac(2) = rrfracmin(2) + m
                  do n = -newaldreal(3), newaldreal(3)
                     rfrac(3) = rrfracmin(3) + n

                     ! convert to absolute coordinates
                     r(:) = matmul(H, rfrac(:))

                     dist2 = r(1)**2 + r(2)**2 + r(3)**2
                     if (dist2 < ewaldrealc2) then
                        dist = dsqrt(dist2)
                        ! calculate short-range contribution
                        ! note: don't need factor of 1/2 bc summing over j,i>j
                        esum = esum + q1*q2*erfc(ewaldalpha*dist)/dist - vshift
                     endif ! within cutoff
                  enddo ! n
               enddo ! m
            enddo ! l
         enddo ! atoms j
      enddo ! atoms i

      ! include contribution due to interaction of j1 with periodic images of itself
      ! (separated due to efficiency)
      ! U_periodic-self = 0.5*sum_L(erfc(alpha*rL)/rL)*sum_i(Qi**2)
      sumq2 = 0.0d0
      do j1 = 1, natoms
        q1 = stchrg(j1)
        sumq2 = sumq2 + q1*q1
      enddo

      ! calculate vertical shift
      vshift = erfc(ewaldalpha*ewaldrealc)/(2*ewaldrealc)

      eself = 0.0d0
      ! iterate over boxes
      do l = -newaldreal(1), newaldreal(1)
         rfrac(1) = l
         do m = -newaldreal(2), newaldreal(2)
            rfrac(2) = m
            do n = -newaldreal(3), newaldreal(3)
               rfrac(3) = n

               ! check not in central box
               if (.not.(l.eq.0.and.m.eq.0.and.n.eq.0)) then
                  ! convert from fractional to absolute
                  r(:) = matmul(H, rfrac(:))

                  dist2 = r(1)**2 + r(2)**2 + r(3)**2
                  if (dist2 < ewaldrealc2) then
                     dist = dsqrt(dist2)
                     ! calculate short-range contribution
                     ! note: need factor of 1/2 to prevent double-counting
                     eself = eself + erfc(ewaldalpha*dist)/(2.0d0*dist) - vshift
                  endif ! within cutoff
               endif ! not in central box
            enddo ! n
         enddo ! m
      enddo ! l

      esum = esum + sumq2*eself

      ! compensate for within-rigidbody interactions
      ! calculate within-rigidbody energy using exact Coulomb sum
      ! U_wrb = sum_J(sum_i>j(Qij/rij))
      ! note: don't need factor of 1/2 because summing over i > j
      if (rigidinit) then
         ewrb = 0.0d0
         ! iterate over rigidbodies
         do i = 1, nrigidbody
   
            ! iterate over atoms i
            do j1 = 1, nsiteperbody(i)
               j3 = 3*j1
               q1 = stchrg(j1)
   
               ! iterate over atoms i > j
               do j2 = j1+1, nsiteperbody(i)
                  j4 = 3*j2
                  q2 = stchrg(j2)
   
                  ! calculate rij
                  r(1) = x(j3-2)-x(j4-2)
                  r(2) = x(j3-1)-x(j4-1)
                  r(3) = x(j3)-x(j4)
                  dist2 = r(1)**2 + r(2)**2 + r(3)**2
                  dist = dsqrt(dist2)
   
                  ! calculate within-rigidbody contribution
                  ewrb = ewrb + q1*q2/dist
               enddo ! sites j
            enddo ! sites i
         enddo ! rigidbodies
   
         ! subtract U_wrb
         esum = esum - ewrb
      endif ! rigidinit

      ! compensate for contribution due to self-interaction
      ! U_self-interaction = -alpha*sum_i(Qi**2)/sqrt(pi)
      esum = esum - sumq2*ewaldalpha/dsqrt(pi)

      ereal = ereal + esum

      return
      end subroutine coulombreal_tri

! -----------------------------------------------------------------------------------
! Calculates and stores terms that are needed to calculate structure factors,
! S(k) and S(-k), to facilitate the computation of the reciprocal-space part of the 
! Ewald sum.

! Because the coefficient of the Coulomb term satisfies the geometric combination rule,
! Q_ij = sqrt(Q_ii*Q_jj), a summation over two indices can be converted to two
! summations over one index.

! Assumes orthorhombic unit cell.
! -----------------------------------------------------------------------------------
      subroutine ftdensity_ortho(x, newaldrecip)

      use commons, only: rerhoarray, imrhoarray 

      implicit none

      integer                      :: j1, j3, l, m, n, dims(3)
      integer, intent(in)          :: newaldrecip(3)
      double precision, intent(in) :: x(3*natoms)
      double precision             :: k(3), r(3)
      double precision             :: q1, k2, kdotr, rerho, imrho, ewaldrecipc2
      double precision, parameter  :: pi = 3.141592654D0

      ! reciprocal-space cutoff
      ewaldrecipc2 = ewaldrecipc**2

      ! make sure allocated arrays for structure factors are the correct size
      dims(:) = 2*newaldrecip(1:3)+1 
      if (.not.allocated(rerhoarray)) allocate(rerhoarray(dims(1), dims(2), dims(3)))
      if (.not.allocated(imrhoarray)) allocate(imrhoarray(dims(1), dims(2), dims(3)))

      if (.not.(size(rerhoarray,1).eq.dims(1).and.size(rerhoarray,2).eq.dims(2).and.size(rerhoarray,3).eq.dims(3))) then
         deallocate(rerhoarray) 
         deallocate(imrhoarray)
         allocate(rerhoarray(dims(1), dims(2), dims(3)))
         allocate(imrhoarray(dims(1), dims(2), dims(3)))
      endif

      ! iterate over boxes and calculate reciprocal lattice vectors
      ! note: because of anti/symmetry in sine and cosine functions,
      ! only need to calculate terms for half of the k-values
      do l = 0,newaldrecip(1)
         k(1) = 2*pi*l/box_params(1)
         do m = -newaldrecip(2), newaldrecip(2)
            k(2) = 2*pi*m/box_params(2)
            do n = -newaldrecip(3), newaldrecip(3)
               k(3) = 2*pi*n/box_params(3)
               ! check not in central box
               if (.not.(l.eq.0.and.m.eq.0.and.n.eq.0)) then
                  k2 = k(1)**2 + k(2)**2 + k(3)**2
                  rerho=0.0d0
                  imrho=0.0d0
                  if (k2 < ewaldrecipc2) then
                     ! iterate over atoms
                     do j1 = 1, natoms
                        j3 = 3*j1
                        q1 = stchrg(j1)
                        r(:) = x(j3-2:j3)
                        ! dot product of k and ri
                        kdotr = k(1)*r(1) + k(2)*r(2) + k(3)*r(3)
                        ! rerho = sum_i(Qi*cos(k*ri))
                        rerho = rerho + q1*dcos(kdotr)
                        ! imrho = sum_i(Qi*sin(k*ri))
                        imrho = imrho + q1*dsin(kdotr)
                     enddo ! atoms
                  endif ! within cutoff
                  ! store rerho and imrho values
                  rerhoarray(-l+newaldrecip(1)+1, -m+newaldrecip(2)+1, -n+newaldrecip(3)+1) = rerho
                  rerhoarray(l+newaldrecip(1)+1, m+newaldrecip(2)+1, n+newaldrecip(3)+1) = rerho
                  imrhoarray(-l+newaldrecip(1)+1, -m+newaldrecip(2)+1, -n+newaldrecip(3)+1) = -imrho
                  imrhoarray(l+newaldrecip(1)+1, m+newaldrecip(2)+1, n+newaldrecip(3)+1) = imrho
               endif ! not in central box
            enddo ! n
         enddo ! m
      enddo ! l

      return
      end subroutine ftdensity_ortho

! -----------------------------------------------------------------------------------
! Calculates and stores terms that are needed to calculate structure factors,
! S(k) and S(-k) to facilitate the computation of the reciprocal-space part of the 
! Ewald sum.

! Because the coefficient of the Coulomb term satisfies the geometric combination rule,
! Q_ij = sqrt(Q_ii*Q_jj), a summation over two indices can be converted to two
! summations over one index.

! Assumes triclinic unit cell.
! -----------------------------------------------------------------------------------
      subroutine ftdensity_tri(x, newaldrecip)

      use commons, only: rerhoarray, imrhoarray 
      use cartdist, only: get_reciplatvec

      implicit none

      integer                      :: j1, j3, l, m, n, dims(3)
      integer, intent(in)          :: newaldrecip(3)
      double precision, intent(in) :: x(3*natoms)
      double precision             :: k(3), r(3), reciplatvec(3,3), reciplatvec_grad(3,3,6)
      double precision             :: q1, k2, kdotr, rerho, imrho, ewaldrecipc2

      ! reciprocal-space cutoff
      ewaldrecipc2 = ewaldrecipc**2

      ! make sure allocated arrays for structure factors are the correct size
      dims(:) = 2*newaldrecip(1:3)+1 
      if (.not.allocated(rerhoarray)) allocate(rerhoarray(dims(1), dims(2), dims(3)))
      if (.not.allocated(imrhoarray)) allocate(imrhoarray(dims(1), dims(2), dims(3)))

      if (.not.(size(rerhoarray,1).eq.dims(1).and.size(rerhoarray,2).eq.dims(2).and.size(rerhoarray,3).eq.dims(3))) then
         deallocate(rerhoarray) 
         deallocate(imrhoarray)
         allocate(rerhoarray(dims(1), dims(2), dims(3)))
         allocate(imrhoarray(dims(1), dims(2), dims(3)))
      endif

      ! get reciprocal lattice vectors
      call get_reciplatvec(reciplatvec, reciplatvec_grad, .false.)

      ! iterate over boxes and calculate reciprocal lattice vectors
      ! note: because of anti/symmetry in sine and cosine functions,
      ! only need to calculate terms for half of the k-values
      do l = 0,newaldrecip(1)
         do m = -newaldrecip(2), newaldrecip(2)
            do n = -newaldrecip(3), newaldrecip(3)
               ! check not in central box
               if (.not.(l.eq.0.and.m.eq.0.and.n.eq.0)) then
                  k = l*reciplatvec(:,1) + m*reciplatvec(:,2) + n*reciplatvec(:,3)
                  k2 = k(1)**2 + k(2)**2 + k(3)**2
                  rerho=0.0d0
                  imrho=0.0d0
                  if (k2 < ewaldrecipc2) then
                     ! iterate over atoms
                     do j1 = 1, natoms
                        j3 = 3*j1
                        q1 = stchrg(j1)
                        r(:) = x(j3-2:j3)
                        ! dot product of k and ri
                        kdotr = k(1)*r(1) + k(2)*r(2) + k(3)*r(3)
                        ! rerho = sum_i(Qi*cos(k*ri))
                        rerho = rerho + q1*dcos(kdotr)
                        ! imrho = sum_i(Qi*sin(k*ri))
                        imrho = imrho + q1*dsin(kdotr)
                     enddo ! atoms
                  endif ! within cutoff
                  ! store rerho and imrho values
                  rerhoarray(-l+newaldrecip(1)+1, -m+newaldrecip(2)+1, -n+newaldrecip(3)+1) = rerho
                  rerhoarray(l+newaldrecip(1)+1, m+newaldrecip(2)+1, n+newaldrecip(3)+1) = rerho
                  imrhoarray(-l+newaldrecip(1)+1, -m+newaldrecip(2)+1, -n+newaldrecip(3)+1) = -imrho
                  imrhoarray(l+newaldrecip(1)+1, m+newaldrecip(2)+1, n+newaldrecip(3)+1) = imrho
               endif ! not in central box
            enddo ! n
         enddo ! m
      enddo ! l

      return
      end subroutine ftdensity_tri

! -----------------------------------------------------------------------------------
! Calculates long-range contribution to Coulomb energy. Uses terms calculated by
! ftdensity_ortho subroutine (structure factors) to simplify computation.

! Assumes orthorhombic unit cell.
! -----------------------------------------------------------------------------------
      subroutine coulombrecip_ortho(x, newaldrecip, erecip)

      use commons, only: rerhoarray, imrhoarray
      use cartdist, only: get_volume

      implicit none

      integer                         :: l, m, n
      integer, intent(in)             :: newaldrecip(3)
      double precision, intent(in)    :: x(3*natoms)
      double precision                :: vol, ewaldrecipc2, k(3)
      double precision                :: k2, rerho, imrho, esum
      double precision, intent(inout) :: erecip
      double precision, parameter     :: pi = 3.141592654D0

      ! cell volume
      call get_volume(vol)
      ! reciprocal-space cutoff
      ewaldrecipc2 = ewaldrecipc**2
      ! compute / store structure factors
      call ftdensity_ortho(x, newaldrecip)
      esum = 0.0d0

      ! compute reciprocal-space sum
      ! U_f = (2*pi/V)*(sum_k(exp(-k**2/4*alpha**2)*S(k)S(-k)/k**2)
      ! iterate over boxes and calculate reciprocal lattice vectors
      do l = -newaldrecip(1), newaldrecip(1)
         k(1) = 2*pi*l/box_params(1)
         do m = -newaldrecip(2), newaldrecip(2)
            k(2) = 2*pi*m/box_params(2)
            do n = -newaldrecip(3), newaldrecip(3)
               k(3) = 2*pi*n/box_params(3)
               ! check not in central box
               if (.not.(l.eq.0.and.m.eq.0.and.n.eq.0)) then
                  k2 = k(1)**2 + k(2)**2 + k(3)**2
                  if (k2 < ewaldrecipc2) then
                     ! get structure factors
                     rerho = rerhoarray(l+newaldrecip(1)+1, m+newaldrecip(2)+1, n+newaldrecip(3)+1)
                     imrho = imrhoarray(l+newaldrecip(1)+1, m+newaldrecip(2)+1, n+newaldrecip(3)+1)
                     ! calculate long-range contribution
                     esum = esum + dexp(-k2/(4.0d0*ewaldalpha**2))*(rerho**2+imrho**2)/k2
                  endif ! within cutoff
               endif ! not in central box
            enddo ! n
         enddo ! m
      enddo ! l

      ! multiply sum by factor of 2*pi/vol
      erecip = erecip + 2.0d0*pi*esum/vol

      return
      end subroutine coulombrecip_ortho

! -----------------------------------------------------------------------------------
! Calculates long-range contribution to Coulomb energy. Uses terms calculated by
! ftdensity_ortho subroutine (structure factors) to simplify computation.

! Assumes triclinic unit cell.
! -----------------------------------------------------------------------------------
      subroutine coulombrecip_tri(x, newaldrecip, erecip)

      use commons, only: rerhoarray, imrhoarray
      use cartdist, only: get_volume, get_reciplatvec

      implicit none

      integer                         :: l, m, n
      integer, intent(in)             :: newaldrecip(3)
      double precision, intent(in)    :: x(3*natoms)
      double precision                :: reciplatvec(3,3), reciplatvec_grad(3,3,6), k(3)
      double precision                :: vol, ewaldrecipc2, k2, rerho, imrho, esum
      double precision, intent(inout) :: erecip
      double precision, parameter     :: pi = 3.141592654D0

      ! cell volume
      call get_volume(vol)
      ! reciprocal lattice vectors
      call get_reciplatvec(reciplatvec, reciplatvec_grad, .false.)
      ! reciprocal-space cutoff
      ewaldrecipc2 = ewaldrecipc**2
      ! compute / store structure factors
      call ftdensity_tri(x, newaldrecip)
      esum = 0.0d0

      ! compute reciprocal-space sum
      ! U_f = (2*pi/V)*(sum_k(exp(-k**2/4*alpha**2)*S(k)S(-k)/k**2)
      ! iterate over boxes and calculate reciprocal lattice vectors
      do l = -newaldrecip(1), newaldrecip(1)
         do m = -newaldrecip(2), newaldrecip(2)
            do n = -newaldrecip(3), newaldrecip(3)
               ! check not in central box
               if (.not.(l.eq.0.and.m.eq.0.and.n.eq.0)) then
                  k = l*reciplatvec(:,1) + m*reciplatvec(:,2) + n*reciplatvec(:,3)
                  k2 = k(1)**2 + k(2)**2 + k(3)**2
                  if (k2 < ewaldrecipc2) then
                     ! get structure factors
                     rerho = rerhoarray(l+newaldrecip(1)+1, m+newaldrecip(2)+1, n+newaldrecip(3)+1)
                     imrho = imrhoarray(l+newaldrecip(1)+1, m+newaldrecip(2)+1, n+newaldrecip(3)+1)
                     ! calculate long-range contribution
                     esum = esum + dexp(-k2/(4.0d0*ewaldalpha**2))*(rerho**2+imrho**2)/k2
                  endif ! within cutoff
               endif ! not in central box
            enddo ! n
         enddo ! m
      enddo ! l

      ! multiply sum by factor of 2*pi/vol
      erecip = erecip + 2.0d0*pi*esum/vol

      return
      end subroutine coulombrecip_tri

! -----------------------------------------------------------------------------------
! Calculates the real-space contribution to the gradient with respects to atomic
! positions. Also calculates real-space contribution to the gradient wrt lattice
! vectors, if BOXDERIVT is true.

! Assumes orthorhombic unit cell.
! -----------------------------------------------------------------------------------
      subroutine coulombrealgrad_ortho(x, newaldreal, g)

      use genrigid, only: rigidinit, nrigidbody, nsiteperbody, rigidgroups, gr_weights

      implicit none

      integer                         :: j1, j3, j2, j4, l, m, n
      integer, intent(in)             :: newaldreal(3)
      double precision, intent(in)    :: x(3*natoms)
      double precision, intent(inout) :: g(3*natoms)
      double precision                :: com(3), mass, comcoords(3*natoms)
      double precision                :: rss(3), rmin(3), r(3), rcommin(3), rcom(3), f(3)
      double precision                :: ewaldrealc2, q1, q2, mul, dist, dist2
      double precision, parameter     :: pi = 3.141592654d0

      ! if rigid bodies, calculate COM coordinates
      ! to compute box derivatives
      if (rigidinit.and.boxderivt) then
         do j1 = 1, nrigidbody
            ! calculate COM
            com(:) = 0.0d0
            mass = 0.0d0
            do j2 = 1, nsiteperbody(j1)
               j3 = rigidgroups(j2, j1)
               com(1:3) = com(1:3) + x(3*j3-2:3*j3)*gr_weights(j3)
               mass = mass + gr_weights(j3)
            enddo
            com(1:3) = com(1:3) / mass
            ! store COM coords
            do j2 = 1, nsiteperbody(j1)
               j3 = rigidgroups(j2, j1)
               comcoords(3*j3-2:3*j3) = com(1:3)
            enddo
         enddo
      endif

      ! real-space cutoff
      ewaldrealc2 = ewaldrealc**2

      ! compute real-space contribution to gradient
      ! G_r = sum_L,i>j(-Qij*r*((erfc(alpha*rij)/(alpha*dist)**3) + 2*alpha*exp(-(alpha*rij)**2)/(sqrt(pi)*rij**2))
      ! iterate over atoms i
      do j1 = 1, natoms
         j3 = 3*j1
         q1 = stchrg(j1)

         ! iterate over atoms i > j
         do j2 = j1+1, natoms
            j4 = 3*j2
            q2 = stchrg(j2)

            ! get distance between atoms
            rss(1) = x(j3-2)-x(j4-2)
            rss(2) = x(j3-1)-x(j4-1)
            rss(3) = x(j3)-x(j4) 
            ! minimum image convention
            rmin(1) = rss(1) - box_params(1)*anint(rss(1)/box_params(1))
            rmin(2) = rss(2) - box_params(2)*anint(rss(2)/box_params(2))
            rmin(3) = rss(3) - box_params(3)*anint(rss(3)/box_params(3))

            ! get minimum distance between COM
            ! NOTE: use rss for minimum image convention to ensure COM corresponds to right atoms
            if (rigidinit.and.boxderivt) then
               rcommin(1) = comcoords(j3-2)-comcoords(j4-2) - box_params(1)*anint(rss(1)/box_params(1))
               rcommin(2) = comcoords(j3-1)-comcoords(j4-1) - box_params(2)*anint(rss(2)/box_params(2))
               rcommin(3) = comcoords(j3)-comcoords(j4) - box_params(3)*anint(rss(3)/box_params(3))
            endif

            ! get gradient contribution per box
            f(:) = 0.0d0

            ! iterate over boxes
            do l = -newaldreal(1), newaldreal(1)
               r(1) = rmin(1)+box_params(1)*l
               do m = -newaldreal(2), newaldreal(2)
                  r(2) = rmin(2)+box_params(2)*m
                  do n = -newaldreal(3), newaldreal(3)
                     r(3) = rmin(3)+box_params(3)*n

                     if (rigidinit.and.boxderivt) then
                        rcom(1) = rcommin(1)+box_params(1)*l
                        rcom(2) = rcommin(2)+box_params(2)*m
                        rcom(3) = rcommin(3)+box_params(3)*n
                     endif

                     dist2 = r(1)**2 + r(2)**2 + r(3)**2
                     if (dist2 < ewaldrealc2) then
                        dist = dsqrt(dist2)
                        ! calculate short-range gradient contribution per box
                        mul = q1*q2*(erfc(ewaldalpha*dist)/dist**3 + 2.0d0*ewaldalpha*dexp(-(ewaldalpha*dist)**2)/(dsqrt(pi)*dist2))
                        f(1) = f(1) + mul*r(1)
                        f(2) = f(2) + mul*r(2)
                        f(3) = f(3) + mul*r(3)

                        ! compute contribution to box derivatives
                        if (boxderivt) then
                           if (rigidinit) then
                              box_paramsgrad(1:3) = box_paramsgrad(1:3) - mul*r(1:3)*rcom(1:3)/box_params(1:3)
                           else ! not rigid bodies
                              box_paramsgrad(1:3) = box_paramsgrad(1:3) - mul*r(1:3)*r(1:3)/box_params(1:3)
                           endif 
                        endif 

                     endif ! within cutoff
                  enddo ! n
               enddo ! m
            enddo ! l

            ! add gradient contribution
            g(j3-2) = g(j3-2)-f(1)
            g(j3-1) = g(j3-1)-f(2)
            g(j3)   = g(j3)-f(3)
            g(j4-2) = g(j4-2)+f(1)
            g(j4-1) = g(j4-1)+f(2)
            g(j4)   = g(j4)+f(3)
         enddo ! atoms j
      enddo ! atoms i

      ! include contribution due to interaction of j1 with periodic images of itself
      ! (separated due to efficiency)
      ! G_periodic-self = sum_L(Qi**2*rL*(erfc(alpha*rL)/rL**3 + 2*alpha*exp(-(alpha*rL)**2)/(sqrt(pi)*rL**2)))
      ! iterate over boxes
      do l = -newaldreal(1), newaldreal(1)
         rmin(1) = box_params(1)*l
         do m = -newaldreal(2), newaldreal(2)
            rmin(2) = box_params(2)*m
            do n = -newaldreal(3), newaldreal(3)
               rmin(3) = box_params(3)*n
               ! check not in central box
               if (.not.(l.eq.0.and.m.eq.0.and.n.eq.0)) then
                  dist2 = rmin(1)**2 + rmin(2)**2 + rmin(3)**2
                  if (dist2 < ewaldrealc2) then
                     dist = dsqrt(dist2)

                     if (rigidinit.and.boxderivt) then
                        rcom(1) = box_params(1)*l
                        rcom(2) = box_params(2)*m
                        rcom(3) = box_params(3)*n
                     endif

                     mul = erfc(ewaldalpha*dist)/dist**3 + 2.0d0*ewaldalpha*dexp(-(ewaldalpha*dist)**2)/(dsqrt(pi)*dist2)
                     ! iterate over atoms and calculate gradient terms
                     do j1 = 1, natoms
                        j3 = 3*j1
                        q1 = stchrg(j1)
                        g(j3-2) = g(j3-2) - q1*q1*mul*rmin(1)
                        g(j3-1) = g(j3-1) - q1*q1*mul*rmin(2)
                        g(j3)   = g(j3)   - q1*q1*mul*rmin(3)

                        ! compute contribution to box derivatives
                        if (boxderivt) then
                           if (rigidinit) then
                              box_paramsgrad(1:3) = box_paramsgrad(1:3) - q1*q1*mul*rmin(1:3)*rcom(1:3)/box_params(1:3)
                           else ! not rigid bodies
                              box_paramsgrad(1:3) = box_paramsgrad(1:3) + q1*q1*mul*rmin(1:3)*rmin(1:3)/box_params(1:3)
                           endif
                        endif 

                     enddo ! atoms
                  endif ! within cutoff
               endif ! not in central box
            enddo ! n
         enddo ! m
      enddo ! l

      return
      end subroutine coulombrealgrad_ortho

! -----------------------------------------------------------------------------------
! Calculates the real-space contribution to the gradient with respects to atomic
! positions. Also calculates real-space contribution to the gradient wrt lattice
! vectors, if BOXDERIVT is true.

! Assumes triclinic unit cell.
! -----------------------------------------------------------------------------------
      subroutine coulombrealgrad_tri(x, newaldreal, g)

      use genrigid, only: rigidinit, nrigidbody, nsiteperbody, rigidgroups, &
      &                   gr_weights, inversematrix
      use cartdist, only: build_H

      implicit none

      integer                         :: j1, j3, j2, j4, l, m, n, idx
      integer, intent(in)             :: newaldreal(3)
      double precision, intent(in)    :: x(3*natoms)
      double precision, intent(inout) :: g(3*natoms)
      double precision                :: com(3), mass, comcoords(3*natoms)
      double precision                :: rr(3), rrfrac(3), rrfracmin(3), r(3), f(3)
      double precision                :: rcom(3), rcomfracmin(3), rcomfrac(3)
      double precision                :: H(3,3), H_grad(3,3,6), H_inverse(3,3)
      double precision                :: ewaldrealc2, q1, q2, mul, dist, dist2
      double precision, parameter     :: pi = 3.141592654d0

      ! if rigid bodies, calculate COM coordinates
      ! to compute box derivatives
      if (rigidinit.and.boxderivt) then
         do j1 = 1, nrigidbody
            ! calculate COM
            com(:) = 0.0d0
            mass = 0.0d0
            do j2 = 1, nsiteperbody(j1)
               j3 = rigidgroups(j2, j1)
               com(1:3) = com(1:3) + x(3*j3-2:3*j3)*gr_weights(j3)
               mass = mass + gr_weights(j3)
            enddo
            com(1:3) = com(1:3) / mass
            ! store COM coords
            do j2 = 1, nsiteperbody(j1)
               j3 = rigidgroups(j2, j1)
               comcoords(3*j3-2:3*j3) = com(1:3)
            enddo
         enddo
      endif

      ! real-space cutoff
      ewaldrealc2 = ewaldrealc**2

      ! get H matrix and inverse
      call build_H(H, H_grad, boxderivt)
      call inversematrix(H, H_inverse)

      ! compute real-space contribution to gradient
      ! G_r = sum_L,i>j(-Qij*r*((erfc(alpha*rij)/(alpha*dist)**3) + 2*alpha*exp(-(alpha*rij)**2)/(sqrt(pi)*rij**2))
      ! iterate over atoms i
      do j1 = 1, natoms
         j3 = 3*j1
         q1 = stchrg(j1)

         ! iterate over atoms i > j
         do j2 = j1+1, natoms
            j4 = 3*j2
            q2 = stchrg(j2)

            ! get distance between atoms
            rr(:) = x(j3-2:j3) - x(j4-2:j4)
            ! convert to fractional coordinates
            rrfrac(:) = matmul(H_inverse, rr(:))
            ! minimum image convention
            rrfracmin(1) = rrfrac(1) - anint(rrfrac(1))
            rrfracmin(2) = rrfrac(2) - anint(rrfrac(2))
            rrfracmin(3) = rrfrac(3) - anint(rrfrac(3))

            ! get minimum distance between COM
            if (rigidinit.and.boxderivt) then
               rcom(:) = comcoords(j3-2:j3) - comcoords(j4-2:j4)
               ! convert to fractional coords
               rcomfracmin(:) = matmul(H_inverse, rcom(:))
               ! minimum image convention
               ! NOTE: use rrfrac for minimum image convention to ensure COM corresponds to right atoms
               rcomfracmin(1) = rcomfracmin(1) - anint(rrfrac(1))
               rcomfracmin(2) = rcomfracmin(2) - anint(rrfrac(2))
               rcomfracmin(3) = rcomfracmin(3) - anint(rrfrac(3))
            endif

            ! get gradient contribution per box
            f(:) = 0.0d0

            ! iterate over boxes
            do l = -newaldreal(1), newaldreal(1)
               rrfrac(1) = rrfracmin(1) + l
               do m = -newaldreal(2), newaldreal(2)
                  rrfrac(2) = rrfracmin(2) + m
                  do n = -newaldreal(3), newaldreal(3)
                     rrfrac(3) = rrfracmin(3) + n

                     ! convert to absolute coordinates
                     r(:) = matmul(H, rrfrac(:))

                     if (rigidinit.and.boxderivt) then
                        rcomfrac(1) = rcomfracmin(1) + l
                        rcomfrac(2) = rcomfracmin(2) + m
                        rcomfrac(3) = rcomfracmin(3) + n
                     endif

                     dist2 = r(1)**2 + r(2)**2 + r(3)**2
                     if (dist2 < ewaldrealc2) then
                        dist = dsqrt(dist2)
                        ! calculate short-range gradient contribution per box
                        mul = q1*q2*(erfc(ewaldalpha*dist)/dist**3 + 2.0d0*ewaldalpha*dexp(-(ewaldalpha*dist)**2)/(dsqrt(pi)*dist2))
                        f(1) = f(1) + mul*r(1)
                        f(2) = f(2) + mul*r(2)
                        f(3) = f(3) + mul*r(3)

                        ! compute contribution to box derivatives
                        if (boxderivt) then
                           if (rigidinit) then
                              ! iterate over cell parameters
                              do idx = 1,6
                                 box_paramsgrad(idx) = box_paramsgrad(idx) - mul*dot_product(r(:), matmul(H_grad(:,:,idx),rcomfrac))
                              enddo
                           else ! not rigid bodies
                              ! iterate over cell parameters
                              do idx = 1, 6
                                 box_paramsgrad(idx) = box_paramsgrad(idx) - mul*dot_product(r(:), matmul(H_grad(:,:,idx), rrfrac))
                              enddo
                           endif 
                        endif 

                     endif ! within cutoff
                  enddo ! n
               enddo ! m
            enddo ! l

            ! add gradient contribution
            g(j3-2) = g(j3-2)-f(1)
            g(j3-1) = g(j3-1)-f(2)
            g(j3)   = g(j3)-f(3)
            g(j4-2) = g(j4-2)+f(1)
            g(j4-1) = g(j4-1)+f(2)
            g(j4)   = g(j4)+f(3)
         enddo ! atoms j
      enddo ! atoms i

      ! include contribution due to interaction of j1 with periodic images of itself
      ! (separated due to efficiency)
      ! G_periodic-self = sum_L(Qi**2*rL*(erfc(alpha*rL)/rL**3 + 2*alpha*exp(-(alpha*rL)**2)/(sqrt(pi)*rL**2)))
      ! iterate over boxes
      do l = -newaldreal(1), newaldreal(1)
         rrfrac(1) = l
         do m = -newaldreal(2), newaldreal(2)
            rrfrac(2) = m
            do n = -newaldreal(3), newaldreal(3)
               rrfrac(3) = n
               ! check not in central box
               if (.not.(l.eq.0.and.m.eq.0.and.n.eq.0)) then
                  ! convert from fractional to absolute
                  r(:) = matmul(H, rrfrac(:))

                  dist2 = r(1)**2 + r(2)**2 + r(3)**2
                  if (dist2 < ewaldrealc2) then
                     dist = dsqrt(dist2)

                     if (rigidinit.and.boxderivt) then
                        rcomfrac(1) = l
                        rcomfrac(2) = m
                        rcomfrac(3) = n
                     endif

                     mul = erfc(ewaldalpha*dist)/dist**3 + 2.0d0*ewaldalpha*dexp(-(ewaldalpha*dist)**2)/(dsqrt(pi)*dist2)
                     ! iterate over atoms and calculate gradient terms
                     do j1 = 1, natoms
                        j3 = 3*j1
                        q1 = stchrg(j1)
                        g(j3-2) = g(j3-2) - q1*q1*mul*r(1)
                        g(j3-1) = g(j3-1) - q1*q1*mul*r(2)
                        g(j3)   = g(j3)   - q1*q1*mul*r(3)

                        ! compute contribution to box derivatives
                        if (boxderivt) then
                           if (rigidinit) then
                              ! iterate over cell parameters
                              do idx = 1,6
                                 box_paramsgrad(idx) = box_paramsgrad(idx) - mul*dot_product(r(:), matmul(H_grad(:,:,idx),rcomfrac))
                              enddo
                           else ! not rigid bodies
                              ! iterate over cell parameters
                              do idx = 1, 6
                                 box_paramsgrad(idx) = box_paramsgrad(idx) - mul*dot_product(r(:), matmul(H_grad(:,:,idx), rrfrac))
                              enddo
                           endif
                        endif

                     enddo ! atoms
                  endif ! within cutoff
               endif ! not in central box
            enddo ! n
         enddo ! m
      enddo ! l

      return
      end subroutine coulombrealgrad_tri

! -----------------------------------------------------------------------------------
! Calculates the reciprocal-space contribution to the gradient with respects to atomic
! positions. Also calculates reciprocal-space contribution to the gradient wrt lattice
! vectors, if BOXDERIVT is true. Uses structure factors to simplify computation.

! Assumes orthorhombic unit cell.
! -----------------------------------------------------------------------------------
      subroutine coulombrecipgrad_ortho(x, newaldrecip, g)

      use commons, only: rerhoarray, imrhoarray
      use genrigid, only: rigidinit, nrigidbody, nsiteperbody, rigidgroups, gr_weights
      use cartdist, only: get_volume

      implicit none

      integer                         :: l, m, n, j1, j2, j3
      integer, intent(in)             :: newaldrecip(3)
      double precision, intent(in)    :: x(3*natoms)
      double precision, intent(inout) :: g(3*natoms)
      double precision                :: vol, ewaldrecipc2, k(3), r(3)
      double precision                :: k2, kdotr, rerho, imrho, q1, mul, mul2
      double precision                :: com(3), mass, comcoords(3*natoms)
      double precision, parameter     :: pi = 3.141592654D0

      ! cell volume
      call get_volume(vol)
      ! reciprocal-space cutoff
      ewaldrecipc2 = ewaldrecipc**2

      ! if rigid bodies, compute COM coords
      ! to compute box derivatives
      if (rigidinit.and.boxderivt) then
         do j1 = 1, nrigidbody
            com(:) = 0.0d0
            mass = 0.0d0
            ! compute COM
            do j2 = 1, nsiteperbody(j1)
               j3 = rigidgroups(j2, j1)
               com(1:3) = com(1:3) + x(3*j3-2:3*j3)*gr_weights(j3)
               mass = mass + gr_weights(j3)
            enddo
            com(1:3) = com(1:3) / mass
            ! store COM coords
            do j2 = 1, nsiteperbody(j1)
               j3 = rigidgroups(j2, j1)
               comcoords(3*j3-2:3*j3) = com(1:3)
            enddo
         enddo
      endif

      ! compute reciprocal-space gradient
      ! G_f = (-4*pi/vol)*q*sum_k((k/k2)*exp(-k2/4*alpha)*(sin(k*r)*sum_i(qi*cos(k*ri)) - cos(k*r)*sum_i(qi*sin(k*ri))))
      ! iterate over boxes and calculate repciprocal lattice vectors
      do l = -newaldrecip(1), newaldrecip(1)
         k(1) = 2*pi*l/box_params(1)
         do m = -newaldrecip(2), newaldrecip(2)
            k(2) = 2*pi*m/box_params(2)
            do n = -newaldrecip(3), newaldrecip(3)
               k(3) = 2*pi*n/box_params(3)
               ! check not in central box
               if (.not.(l.eq.0.and.m.eq.0.and.n.eq.0)) then
                  k2 = k(1)**2 + k(2)**2 + k(3)**2
                  if (k2 < ewaldrecipc2) then
                     ! calculate multiplicative factor
                     mul = -4.0d0*pi*dexp(-k2/(4.0d0*ewaldalpha**2))/(vol*k2)
                     ! get structure factors
                     rerho = rerhoarray(l+newaldrecip(1)+1, m+newaldrecip(2)+1, n+newaldrecip(3)+1)
                     imrho = imrhoarray(l+newaldrecip(1)+1, m+newaldrecip(2)+1, n+newaldrecip(3)+1)

                     ! add contribution to box derivatives
                     if (boxderivt) then
                        box_paramsgrad(1:3) = box_paramsgrad(1:3) + mul*(rerho**2+imrho**2)* &
                                              (1.0d0 - (k2 + 4.0d0*ewaldalpha**2)*k(1:3)*k(1:3)/ &
                                              (2.0d0*ewaldalpha**2*k2))/(2.0d0*box_params(1:3))
                     endif

                     ! iterate over atoms and calculate long-range gradient terms
                     do j1 = 1, natoms
                        j3 = 3*j1
                        r(:) = x(j3-2:j3)
                        kdotr = k(1)*r(1) + k(2)*r(2) + k(3)*r(3)
                        q1 = stchrg(j1)
                        mul2 = mul*q1*(dsin(kdotr)*rerho - dcos(kdotr)*imrho)
                        
                        ! add contribution to gradient 
                        g(j3-2) = g(j3-2) + mul2*k(1)
                        g(j3-1) = g(j3-1) + mul2*k(2)
                        g(j3)   = g(j3)   + mul2*k(3)

                        ! add contribution to box derivatives from rigid bodies
                        ! NOTE: no contribition if not using rigid bodies
                        if (rigidinit.and.boxderivt) then
                           box_paramsgrad(1:3) = box_paramsgrad(1:3) - mul2*k(1:3)*(x(j3-2:j3)-comcoords(j3-2:j3))/box_params(1:3)
                        endif

                     enddo ! atoms

                  endif ! within cutoff
               endif ! not in central box
            enddo ! n
         enddo ! m
      enddo ! l

      return
      end subroutine coulombrecipgrad_ortho

! -----------------------------------------------------------------------------------
! Calculates the reciprocal-space contribution to the gradient with respects to atomic
! positions. Also calculates reciprocal-space contribution to the gradient wrt lattice
! vectors, if BOXDERIVT is true. Uses structure factors to simplify computation.

! Assumes triclinic unit cell.
! -----------------------------------------------------------------------------------
      subroutine coulombrecipgrad_tri(x, newaldrecip, g)

      use commons, only: rerhoarray, imrhoarray
      use genrigid, only: rigidinit, nrigidbody, nsiteperbody, rigidgroups, gr_weights, inversematrix
      use cartdist, only: get_volume, get_reciplatvec, build_H, cart2frac_tri

      implicit none

      integer                         :: l, m, n, j1, j2, j3, idx
      integer, intent(in)             :: newaldrecip(3)
      double precision, intent(in)    :: x(3*natoms)
      double precision, intent(inout) :: g(3*natoms)
      double precision                :: vol, ewaldrecipc2, c(3), s(3), abc, vfact, dvol(6), r(3)
      double precision                :: reciplatvec(3,3), reciplatvec_grad(3,3,6), xfrac(3*natoms)
      double precision                :: H(3,3), H_grad(3,3,6), H_inverse(3,3), k(3), k_grad(3,6)
      double precision                :: k2, kdotr, rerho, imrho, q1, mul, mul2
      double precision                :: com(3), mass, comcoords(3*natoms), comcoordsfrac(3*natoms)
      double precision, parameter     :: pi = 3.141592654D0

      ! cell volume
      call get_volume(vol)
      ! gradient of volume wrt cell parameters
      if (boxderivt) then
         c(:) = dcos(box_params(4:6))
         s(:) = dsin(box_params(4:6))
         abc = box_params(1)*box_params(2)*box_params(3)
         vfact = vol/abc
         dvol(1) = vol/box_params(1)
         dvol(2) = vol/box_params(2)
         dvol(3) = vol/box_params(3)
         dvol(4) = s(1)*(c(1)-c(2)*c(3))
         dvol(5) = s(2)*(c(2)-c(1)*c(3))
         dvol(6) = s(3)*(c(3)-c(1)*c(2))
         dvol(4:6) = abc*dvol(4:6)/vfact
      endif

      ! reciprocal lattice vectors
      call get_reciplatvec(reciplatvec, reciplatvec_grad, boxderivt)
      ! get H matrix and inverse
      call build_H(H, H_grad, boxderivt)
      call inversematrix(H, H_inverse)
      ! get fractional coordinates
      if (boxderivt) call cart2frac_tri(x, xfrac, H_inverse)
      ! reciprocal-space cutoff
      ewaldrecipc2 = ewaldrecipc**2

      ! if rigid bodies, compute COM coords
      ! to compute box derivatives
      if (rigidinit.and.boxderivt) then
         do j1 = 1, nrigidbody
            com(:) = 0.0d0
            mass = 0.0d0
            ! compute COM
            do j2 = 1, nsiteperbody(j1)
               j3 = rigidgroups(j2, j1)
               com(1:3) = com(1:3) + x(3*j3-2:3*j3)*gr_weights(j3)
               mass = mass + gr_weights(j3)
            enddo
            com(1:3) = com(1:3) / mass
            ! store COM coords
            do j2 = 1, nsiteperbody(j1)
               j3 = rigidgroups(j2, j1)
               comcoords(3*j3-2:3*j3) = com(1:3)
            enddo
         enddo
         ! convert to fractional
         call cart2frac_tri(comcoords, comcoordsfrac, H_inverse)
      endif

      ! compute reciprocal-space gradient
      ! G_f = (-4*pi/vol)*q*sum_k((k/k2)*exp(-k2/4*alpha)*(sin(k*r)*sum_i(qi*cos(k*ri)) - cos(k*r)*sum_i(qi*sin(k*ri))))
      ! iterate over boxes and calculate repciprocal lattice vectors
      do l = -newaldrecip(1), newaldrecip(1)
         do m = -newaldrecip(2), newaldrecip(2)
            do n = -newaldrecip(3), newaldrecip(3)
               ! check not in central box
               if (.not.(l.eq.0.and.m.eq.0.and.n.eq.0)) then
                  k = l*reciplatvec(:,1) + m*reciplatvec(:,2) + n*reciplatvec(:,3)
                  k2 = k(1)**2 + k(2)**2 + k(3)**2
                  if (k2 < ewaldrecipc2) then

                     ! get gradient of reciprocal lattice vector wrt cell parameters
                     if (boxderivt) then
                        do idx = 1,6
                           k_grad(:,idx) = l*reciplatvec_grad(:,1,idx) + m*reciplatvec_grad(:,2,idx) + n*reciplatvec_grad(:,3,idx)
                        enddo
                     endif
                     
                     ! calculate multiplicative factor
                     mul = -4.0d0*pi*dexp(-k2/(4.0d0*ewaldalpha**2))/(vol*k2)
                     ! get structure factors
                     rerho = rerhoarray(l+newaldrecip(1)+1, m+newaldrecip(2)+1, n+newaldrecip(3)+1)
                     imrho = imrhoarray(l+newaldrecip(1)+1, m+newaldrecip(2)+1, n+newaldrecip(3)+1)

                     ! add contribution to box derivatives
                     if (boxderivt) then
                        ! iterate over cell parameters
                        do idx = 1, 6
                            box_paramsgrad(idx) = box_paramsgrad(idx) + &
                                                  mul*(rerho**2+imrho**2)*(dvol(idx)/(2.0d0*vol) + &
                                                  (k2 + 4.0d0*ewaldalpha**2)*dot_product(k, k_grad(:,idx))/ &
                                                  (4.0d0*ewaldalpha**2*k2))
                        enddo
                     endif

                     ! iterate over atoms and calculate long-range gradient terms
                     do j1 = 1, natoms
                        j3 = 3*j1
                        r(:) = x(j3-2:j3)
                        kdotr = k(1)*r(1) + k(2)*r(2) + k(3)*r(3)
                        q1 = stchrg(j1)
                        mul2 = mul*q1*(dsin(kdotr)*rerho - dcos(kdotr)*imrho)
                        
                        ! add contribution to gradient 
                        g(j3-2) = g(j3-2) + mul2*k(1)
                        g(j3-1) = g(j3-1) + mul2*k(2)
                        g(j3)   = g(j3)   + mul2*k(3)

                        ! add contribution to box derivatives
                        if (boxderivt) then
                           if (rigidinit) then
                              ! iterate over cell parameters
                              do idx = 1,6
                                 box_paramsgrad(idx) = box_paramsgrad(idx) + &
                                                       mul2*(dot_product(k_grad(:,idx), r) + &
                                                       dot_product(k, matmul(H_grad(:,:,idx), comcoordsfrac(j3-2:j3))))
                              enddo
                           else ! not rigid bodies
                              ! iterate over cell parameters
                              do idx = 1,6
                                 box_paramsgrad(idx) = box_paramsgrad(idx) + &
                                                       mul2*(dot_product(k_grad(:,idx), r) + &
                                                       dot_product(k, matmul(H_grad(:,:,idx), xfrac(j3-2:j3))))
                              enddo
                           endif
                        endif

                     enddo ! atoms

                  endif ! within cutoff
               endif ! not in central box
            enddo ! n
         enddo ! m
      enddo ! l

      return
      end subroutine coulombrecipgrad_tri

end module
