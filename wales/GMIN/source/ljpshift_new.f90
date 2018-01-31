! dj337: the following subroutine computes the LJ potential
! for periodic systems without using minimum image convention
! (i.e. uses neighbor lists).
! NOTE: only works for orthorhombic cells
! NOTE:  binary LJ not yet implemented

      subroutine ljpshift_new(x, potel, g, gtest)

      use commons, only: natoms, epsab, sigab, box_params, cutoff, &
                         boxderivt, ortho, box_params, box_paramsgrad

      implicit none
      integer                       :: j1, j3, j2, j4, l, m, n, i
      double precision, intent(in)  :: x(3*natoms)
      logical, intent(in)           :: gtest
      double precision, intent(out) :: potel, g(3*natoms)
      double precision              :: eps, sig, sig6, rcut, rcut2, sigrc6
      double precision              :: const, rconst, dist2, idist2, idist6
      double precision              :: sig12, idist8, idist14, dvdr, dist, temp1(3), temp(3)
      double precision              :: val, eself, sigrc12, xj(3), rmin(3), r(3), temp2(3)
      integer                       :: cell_range(3)

      ! confirm that cell is orthorhombic
      if (.not.(ortho)) then
         print *, 'LJPSHIFT_NEW not implemented for triclinic cells!'
         return
      endif

      ! define and calculate constants
      eps = epsab
      sig = sigab
      sig6 = sig**6
      sig12 = sig6**2
      rcut = cutoff*sig
      rcut2 = rcut**2
      sigrc6 = sig6/rcut**6
      sigrc12 = sigrc6**2
      const = 4.0d0*(sigrc6)-7.0d0*sigrc12
      rconst = (6.0d0*sigrc12-3.0d0*sigrc6)/rcut**2

      potel = 0.0d0
      if (gtest) g(:) = 0.0d0
      box_paramsgrad(:) = 0.0d0

      ! determine cell range needed for cutoff
      cell_range(:) = ceiling(rcut/box_params(1:3))

      ! iterate over atoms i
      do j1 = 1, natoms
         j3 = 3*j1

         ! iterate over atoms j > i
         do j2 = j1+1, natoms
            j4 = 3*j2

            rmin(1:3) = x(j3-2:j3) - x(j4-2:j4)
            rmin(1) = rmin(1) - box_params(1)*anint(rmin(1)/box_params(1))
            rmin(2) = rmin(2) - box_params(2)*anint(rmin(2)/box_params(2))
            rmin(3) = rmin(3) - box_params(3)*anint(rmin(3)/box_params(3))

            ! iterate over cells
            do l = -cell_range(1), cell_range(1)
               do m = -cell_range(2), cell_range(2)
                  do n = -cell_range(3), cell_range(3)

                     ! calculate atom-atom separation
                     r(1) = rmin(1)+box_params(1)*l
                     r(2) = rmin(2)+box_params(2)*m
                     r(3) = rmin(3)+box_params(3)*n
                     dist2 = r(1)**2 + r(2)**2 + r(3)**2

                     ! check that distance within cutoff
                     if (dist2 < rcut2) then
                        idist2 = 1.d0/dist2
                        idist6 = idist2**3

                        ! calculate contribution to energy
                        val = 4.d0*eps*(sig6*idist6*(sig6*idist6 - 1.0d0) + rconst*dist2 + const)
                        potel = potel + val

                        if (gtest) then
                           idist8 = idist2*idist6
                           idist14 = idist8*idist6

                           ! calculate partial wrt distance
                           dvdr = -8.0d0*eps*(3.0d0*(2.0d0*idist14*(sig12)-idist8*sig6)-rconst)

                           ! add contribution to gradient
                           g(j3-2:j3) = g(j3-2:j3) + dvdr*r(1:3)
                           g(j4-2:j4) = g(j4-2:j4) - dvdr*r(1:3)

                           if (boxderivt) then
                              if (ortho) then
                                 box_paramsgrad(1:3) = box_paramsgrad(1:3) + dvdr*r(1:3)*r(1:3)/box_params(1:3)
                              endif ! ortho
                           endif ! box derivatives

                        endif ! gradient

                     endif ! cutoff

                  enddo ! n
               enddo ! m
            enddo ! l

         enddo ! atoms j
      enddo ! atoms i

      ! add contribution of atoms with periodic images of itself
      eself = 0.0d0

      ! iterate over cells
      do l = -cell_range(1), cell_range(1)
         do m = -cell_range(2), cell_range(2)
            do n = -cell_range(3), cell_range(3)

               ! check that not in central cell
               if (.not.(l.eq.0.and.m.eq.0.and.n.eq.0)) then

                  ! calculate atom-atom separation
                  r(1) = box_params(1)*l
                  r(2) = box_params(2)*m
                  r(3) = box_params(3)*n
                  dist2 = r(1)**2 + r(2)**2 + r(3)**2

                  ! check that distance within cutoff
                  if (dist2 < rcut2) then
                     idist2 = 1.d0/dist2
                     idist6 = idist2**3

                     ! calculate contribution to energy
                     val = 2.d0*eps*(sig6*idist6*(sig6*idist6 - 1.0d0) + rconst*dist2 + const)
                     eself = eself + val

                     if (gtest) then
                        idist8 = idist2*idist6
                        idist14 = idist8*idist6

                        ! calculate partial wrt distance
                        dvdr = -4.0d0*eps*(3.0d0*(2.0d0*idist14*(sig12)-idist8*sig6)-rconst)

                        ! add contribution to gradient for each atom
                        do j1 = 1, natoms
                           j3 = 3*j1

                           g(j3-2:j3) = g(j3-2:j3) + dvdr*r(1:3)

                           if (boxderivt) then
                              if (ortho) then
                                 box_paramsgrad(1:3) = box_paramsgrad(1:3) + dvdr*r(1:3)*r(1:3)/box_params(1:3)
                              endif ! ortho
                           endif ! box derivatives

                        enddo ! atoms

                     endif ! gradient

                  endif ! cutoff
               endif ! central box
            enddo ! n
         enddo ! m
      enddo ! l

      ! multiply energy by number of atoms
      eself = eself*natoms
      ! total potential energy
      potel = potel + eself

      end subroutine ljpshift_new
