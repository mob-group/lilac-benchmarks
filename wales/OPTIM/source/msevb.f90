!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of OPTIM.
!   
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!   
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!   
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!###################################################################
!################# SUBROUTINE MSEVB ################################
!###################################################################

!###################################################################

! 02/04/04 Changed jacobi routine for LAPACK DSYEVR
! Individual results are consistent to 1E-12 but there seem to be differences in Markov chains none the less

     subroutine msevb(ncoords, systemCoords, assignVBstates, energy, onlyAssignStates)
     use commons
     use msevb_common

     implicit none

!###################################################################
!
! Parameter set used is from:
!    U. Schmitt and G. A. Voth, J. Chem. Phys. 111, 1999, 9361
! and not the earlier set, nor the later one, or even the buggy one.
!
! Calculate the MSEVB energy of a given protonated water cluster.
!
! Steps: (i) assign the possible EVB states in the cluster, (ii)
! transfer the EVB bonding arrangements into coordinate vectors
! xx(i,j) etc. (iii) calculate the Hamiltonian elements by means
! of only the nuclear positions, and finally diagonalizing the
! Hamiltonian to yield the Eigenvectors and Eigenvalues.
!
!###################################################################

! Subroutine arguments

     INTEGER, INTENT(IN) :: ncoords  
     DOUBLE PRECISION, INTENT(IN) :: SYSTEMCOORDS(NCOORDS)
     LOGICAL, INTENT(IN) :: assignVBstates
     DOUBLE PRECISION, INTENT(OUT) :: ENERGY
     LOGICAL, OPTIONAL, INTENT(IN) :: onlyAssignStates   ! Only assign the VB states, do nothing else

! Local variables

     integer i,j,nrot

       DOUBLE PRECISION MIN

! Variables for LAPACK routine

     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: D
     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: V, E
!     integer :: numEVfound
!     integer :: ISUPPZ(2)
!     integer :: localStatus
!     DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: WORK
!     integer, allocatable, dimension(:) :: IWORK

!     DOUBLE PRECISION :: DLAMCH ! EXTERNAL FUNCTION

! Store the atom coordinates in psix, psiy, psiz rather than the single array coords

     j = 0

     do i = 1,ncoords,3
        j = j + 1
        psix(j) = systemCoords(i)
        psiy(j) = systemCoords(i+1)
        psiz(j) = systemCoords(i+2)
     enddo

! Calculate the interatomic distances

     call calcInterAtomicDistances(natoms,psix,psiy,psiz,interAtomicR)

!###################################################################
! Check that there is no Roo vectors smaller in magnitude than
! the RCut cutoff value
! 
!      do i = 3,natoms-4,3
!         do j = i+3,natoms-1,3
!            dx = psix(i) - psix(j) 
!            dy = psiy(i) - psiy(j)   
!            dz = psiz(i) - psiz(j)   
!            R = DSQRT(dx**2 + dy**2 + dz**2)
!            if (R.lt.RCut) then
!            energy = 1.0d03
!            RETURN
!            endif
!         enddo
!      enddo
! 
!###################################################################
! Assign the VB states if required, else use the ones in memory
!
!###################################################################

      if (assignVBstates) then
        call assignvb3
        if (PRESENT(onlyAssignStates).and.onlyAssignStates) return
     endif

! Allocate the relevant arrays

     ALLOCATE(xx(reduced_num_eig,natoms), yy(reduced_num_eig,natoms), zz(reduced_num_eig,natoms))
     ALLOCATE(aij(reduced_num_eig,reduced_num_eig))
     ALLOCATE (d(reduced_num_eig), v(reduced_num_eig, reduced_num_eig), e(reduced_num_eig, reduced_num_eig))
!     ALLOCATE (d(reduced_num_eig), v(reduced_num_eig, 1), e(reduced_num_eig, reduced_num_eig))

! vijexch and grstwfu are the only arrays used externally, all the rest are deallocated at the end of the routine

     if (ALLOCATED(vijexch)) DEALLOCATE(vijexch)
     ALLOCATE(vijexch(reduced_num_eig,reduced_num_eig))

     if (ALLOCATED(grstwfu)) DEALLOCATE(grstwfu)
     ALLOCATE (grstwfu(reduced_num_eig))

! Sometimes want ham for debugging

!     if (ALLOCATED(ham)) DEALLOCATE(ham)
     ALLOCATE(ham(reduced_num_eig,reduced_num_eig))

!###################################################################
! Assign the configurations to xx, the coordinate vector used in
! the subsequent calculation of the Hamiltonian elements
!
!###################################################################
 
     do i = 1,reduced_num_eig
        do j = 1,natoms
           xx(i,j) = psix(atmpl(i,j))
           yy(i,j) = psiy(atmpl(i,j))
           zz(i,j) = psiz(atmpl(i,j))
        enddo
     enddo

!###################################################################
! Build the Hamiltonian: call the routines for constructing each 
! element Hij
!
!###################################################################

     call bildham

       do i = 1,reduced_num_eig
          do j = 1,reduced_num_eig
             e(i,j) = ham(i,j)
          enddo
       enddo

!###################################################################
! Diagonalize the Hamiltonian
!
! Find the Ground State energy and Eigenvector. Note the monitoring
! of which EVB state is occupied in this geometry
!
!###################################################################

      call jacobi(e,reduced_num_eig,d,v,nrot)

! d contains the eigenvalues
! v contains the eigenvectors

!----find the lowest eigenvalue and print out the geometry
      min=1.0d10
      ground_state = 3
      do i = 1,reduced_num_eig
         if (d(i).lt.min) then
            min = d(i)
            ground_state = i
         endif 
       enddo

!----assign the lowest energy from the lowest eigenvalue
      energy = min

! Run LAPACK routine which calculates only the lowest eigenvalue, rather than the full spectrum
! Only need one triangle
! Routine overwrites the Hamiltonian which is why we must pass it a copy rather than the original

!     do i = 1,reduced_num_eig
!        do j = 1,i
!           e(i,j) = ham(i,j)
!        enddo
!     enddo

!     allocate(WORK(33*reduced_num_eig), IWORK(33*reduced_num_eig))
     
!     call DSYEVR('V','I','L',reduced_num_eig,e,reduced_num_eig,0.0d0,0.0d0,1,1,DLAMCH('Safe  minimum'),numEVfound,
!    &              d,v,reduced_num_eig,ISUPPZ,WORK,SIZE(WORK),IWORK,SIZE(IWORK),localStatus)

!     deallocate(WORK,IWORK)

! Check for problems with LAPACK routine

!     if (localStatus.ne.0 .or. numEVfound.ne.1) then
!        print *, 'Error in DSYEVR routine: Return status:', localStatus, 'Num EVs found:', numEVfound
!        stop
!     endif

! Assign energy and ground state wavefunction

!     energy = d(1)
!     grstwfu(:) = v(:,1)

!###################################################################
! Assign the ground state wavefunction to variable
!###################################################################
       do i = 1,reduced_num_eig
          grstwfu(i) = v(i,ground_state)
       enddo      

! Deallocate memory

     DEALLOCATE(xx, yy, zz, aij)
     DEALLOCATE(d, v, e)
     DEALLOCATE(ham)

     RETURN
     END     

!###################################################################
!############### SUBROUTINE VH2OINTER #############################
!###################################################################

     subroutine vh2ointer(vh2o_inter)
     use commons
     use msevb_common

     implicit none

!################################################################### 
! Calculate the Water Intermolecular interaction
! Units are kcal/mole 
!
!################################################################### 

     DOUBLE PRECISION VH2O_INTER(REDUCED_NUM_EIG)
     DOUBLE PRECISION DX_OO,DY_OO,DZ_OO,ROO,ROO_6,ROO_12,LJ_TMP, &
                coulomb_tmp,coulomb_mono,lj_term,dx,dy,dz, &
                Rcoul,qk,qj,e_sq,coulomb,a,c
     integer ii,i,j,kas,kk,start,count,count_sub,jj,kaskas
     character*1 dummy

!###################################################################
! Begin loop over the Eigenstates
!
!################################################################### 
      do 699 ii = 1,reduced_num_eig
!----reinitialize the variables:
        lj_tmp = 0.0d0
        lj_term = 0.0d0
        coulomb = 0.0d0
        coulomb_tmp = 0.0d0
        coulomb_mono = 0.0d0

!################################################################### 
! Determine the LJ component of the Water intermolecular interaction
! for this Eigenstate.
!
!################################################################### 
!     if (QUICK) then
        do i = 6,(natoms-4),3
           jj = atmpl(ii,i) 
           do j = i+3,(natoms-1),3 
              kaskas = atmpl(ii,j)
              lj_tmp = ljr(jj,kaskas)
              lj_term = lj_tmp + lj_term
           enddo
        enddo

!     else
!
!      do i = 6,(natoms-4),3
!        do j = i+3,(natoms-1),3
!           dx_oo = xx(ii,i) - xx(ii,j)
!           dy_oo = yy(ii,i) - yy(ii,j)
!           dz_oo = zz(ii,i) - zz(ii,j)
!           Roo = DSQRT(dx_oo*dx_oo + dy_oo*dy_oo + dz_oo*dz_oo)
!            Roo = 3.1506d0/Roo
!           Roo_6 = Roo**6
!           Roo_12 = Roo_6*Roo_6
!              lj_tmp = (ALJ*Roo_12) - (CLJ*Roo_6)
!            lj_tmp = 4.0d0*(1.522d-01)*(Roo_12 - Roo_6)
!           lj_term = lj_tmp + lj_term
!        enddo
!      enddo
!
!     endif
!
!################################################################### 
! Calculate the coulombic component of the Water intermolecular
! interaction for this Eigenstate. 
!
!################################################################### 
     start = 5

       do j = start,(natoms-3) 
          if (j.eq.(start+3)) then
               start = start + 3
           endif 
          jj = atmpl(ii,j)
          do kas = start+3,natoms  
               kaskas = atmpl(ii,kas)
             coulomb_mono = water_inter_coulomb(jj,kaskas)
             coulomb_tmp = coulomb_tmp + coulomb_mono            
          enddo
       enddo
  
!     do j = start,(natoms-3)
!        if (mod(j,3).eq.0) then
!           qj = -0.834d0
!        else
!          qj = 0.417d0
!        endif
!        if (j.eq.(start+3)) then
!           start = start + 3
!        endif
!        do kas = start+3,natoms
!          dx = xx(ii,j) - xx(ii,kas)
!           dy = yy(ii,j) - yy(ii,kas)
!          dz = zz(ii,j) - zz(ii,kas)
!          Rcoul = DSQRT(dx*dx + dy*dy + dz*dz)
!          if (mod(kas,3).eq.0) then
!             qk = -0.834d0
!          else
!             qk = 0.417d0
!          endif
!              coulomb_mono = (qk*qj*fourpieo)/(Rcoul)
!             coulomb_tmp = coulomb_tmp + coulomb_mono
!        enddo
!     enddo
!
!###################################################################
! Finally work out the Water intermolecular interaction for this
! Eigenstate
!
!###################################################################

     vh2o_inter(ii) = coulomb_tmp + lj_term

!###################################################################
! End of Eigenstate

 699     enddo             

     RETURN 
     END 

!###################################################################
!################### VH2OINTRA SUBROUTINE ##########################
!###################################################################

     subroutine vh2ointra(vh2o_intra)
     use commons
     use msevb_common

     implicit none

!###################################################################
!
! DETERMINATION OF THE INTRAMOLECULAR POTENTIAL FOR H20
!
!################################################################### 

     DOUBLE PRECISION VH2O_INTRA(REDUCED_NUM_EIG)
     DOUBLE PRECISION DX_ROH1,DY_ROH1,DZ_ROH1,DX_ROH2,DY_ROH2,DZ_ROH2
     DOUBLE PRECISION ROH(3),ALPHA
     DOUBLE PRECISION VH2O_INTRA_TMP,VH2O_INTRA_R,VH2O_INTRA_ALPHA &
                ,vh2o_intra_rtmp,vh2o_intra_mono
     integer i,j,l,m,n
     character*4 form
     DOUBLE PRECISION CHECKED_ACOS

!----intialize variables
     do i = 1,reduced_num_eig
        vh2o_intra(i) = 0.0d0
     enddo
        
!################################################################### 
! Calculate the 2(O-H) bond lengths and H-O-H bond angle in each of
! the water molecules of the Eigenstate considered
!
!################################################################### 
     form = 'tip3'
!----Sum over all Eigenstates
     do 799 i = 1,reduced_num_eig
        vh2o_intra_mono = 0.0d00
!----Sum over all the water molecules
     do n = 5,natoms-2,3
!----intialize variables
           roh(1) = 0.0d00
           roh(2) = 0.0d00
           alpha = 0.0d00
!----determine the roh distances in water
        dx_roh1 = psix(atmpl(i,n)) - psix(atmpl(i,n+1))
        dy_roh1 = psiy(atmpl(i,n)) - psiy(atmpl(i,n+1))
           dz_roh1 = psiz(atmpl(i,n)) - psiz(atmpl(i,n+1))      
        roh(1) = dsqrt(dx_roh1*dx_roh1 + dy_roh1*dy_roh1 + dz_roh1*dz_roh1)

        dx_roh2 = psix(atmpl(i,n+2)) - psix(atmpl(i,n+1))
        dy_roh2 = psiy(atmpl(i,n+2)) - psiy(atmpl(i,n+1))   
        dz_roh2 = psiz(atmpl(i,n+2)) - psiz(atmpl(i,n+1))   
           roh(2) = dsqrt(dx_roh2*dx_roh2 + dy_roh2*dy_roh2 + dz_roh2*dz_roh2)

!----determine the bending angle, alpha
        alpha = (dx_roh1*dx_roh2 + dy_roh1*dy_roh2 + dz_roh1*dz_roh2)
        alpha = checked_acos(alpha/(roh(1)*roh(2)))

!###################################################################
! Now determine the Bending and Stretching terms together by looping 
! over the 2 stretched and one angle.
!
!###################################################################
        vh2o_intra_r = 0.0d00
        vh2o_intra_alpha = 0.0d00
        vh2o_intra_rtmp = 0.0d00

        do j=1,2
           vh2o_intra_r = roh(j) - h2oroheq
           vh2o_intra_r = 0.5d0*h2okb*(vh2o_intra_r*vh2o_intra_r) 
              vh2o_intra_rtmp = vh2o_intra_rtmp + vh2o_intra_r
        enddo

        vh2o_intra_alpha = alpha-thetaeq
        vh2o_intra_alpha = 0.5d0*h2oktheta*(vh2o_intra_alpha*vh2o_intra_alpha)
          
!################################################################### 
! Determine the contribution from this, single, water molecule
!
!################################################################### 

        vh2o_intra_mono = vh2o_intra_mono + vh2o_intra_rtmp + vh2o_intra_alpha
 
!----calc'ed the energy for 1 water in 1 eigenstate
!    Finished the loop over a water molecule
       enddo

!----finally assign the vh2o_intra for eigenstate (i) |i>
        vh2o_intra(i) = vh2o_intra_mono
 799     enddo

     return
     end 

!#####################################################################
!################## VH3OINTRA SUBROUTINE #############################
!#####################################################################

     subroutine vh3ointra(vh3o_intra)
     use commons
     use msevb_common

     implicit none

!###################################################################
! Determine the contribution to the Hamiltonian from the H3O+
! Intermolecular interaction
!
! Only requires the nuclear coordinate vector xx(i,j)
!
!###################################################################

     DOUBLE PRECISION VH3O_INTRA(REDUCED_NUM_EIG)
     DOUBLE PRECISION DX_ROH,DY_ROH,DZ_ROH,DX_HH,DY_HH,DZ_HH
     DOUBLE PRECISION ROH(3),ALPHA(3),HH(3)
     DOUBLE PRECISION VH3O_INTRA_TMP,VH3O_INTRA_R,VH3O_INTRA_ALPHA
     integer i,j,keg,l,m
     DOUBLE PRECISION CHECKED_ACOS

!###################################################################
! Determine the stretching (3) and bending (3) vectors
!
!###################################################################

     do 899 i = 1,reduced_num_eig
!----intialize variables
     vh3o_intra_tmp = 0.0d0
        do keg = 1,3
           roh(keg) = 0.0d0
           hh(keg) = 0.0d0
           alpha(keg) = 0.0d0
        enddo
!----determine the roh distances in hydronium
        dx_roh = xx(i,1) - xx(i,3)
        dy_roh = yy(i,1) - yy(i,3)
        dz_roh = zz(i,1) - zz(i,3)
        roh(1) = dsqrt(dx_roh*dx_roh + dy_roh*dy_roh + dz_roh*dz_roh)
        dx_roh = xx(i,2) - xx(i,3)
           dy_roh = yy(i,2) - yy(i,3)
        dz_roh = zz(i,2) - zz(i,3)
           roh(2) = dsqrt(dx_roh*dx_roh + dy_roh*dy_roh + dz_roh*dz_roh)
        dx_roh = xx(i,4) - xx(i,3)
           dy_roh = yy(i,4) - yy(i,3)
           dz_roh = zz(i,4) - zz(i,3)
           roh(3) = dsqrt(dx_roh*dx_roh + dy_roh*dy_roh + dz_roh*dz_roh)
!----determine the bending angles, alpha
!------(a)first determine the hh distances
        dx_hh = xx(i,1) - xx(i,2)
        dy_hh = yy(i,1) - yy(i,2)
           dz_hh = zz(i,1) - zz(i,2)
        hh(1) = dsqrt(dx_hh*dx_hh + dy_hh*dy_hh + dz_hh*dz_hh)
        dx_hh = xx(i,2) - xx(i,4)
           dy_hh = yy(i,2) - yy(i,4)
           dz_hh = zz(i,2) - zz(i,4)
           hh(2) = dsqrt(dx_hh*dx_hh + dy_hh*dy_hh + dz_hh*dz_hh)
        dx_hh = xx(i,4) - xx(i,1)
           dy_hh = yy(i,4) - yy(i,1)
           dz_hh = zz(i,4) - zz(i,1)
           hh(3) = dsqrt(dx_hh*dx_hh + dy_hh*dy_hh + dz_hh*dz_hh)
!------(b)now determine the angles alpha
        alpha(1) = (-(hh(1)*hh(1)) + (roh(1)*roh(1)) + (roh(2)*roh(2)))
        alpha(1) = checked_acos(alpha(1)/(2.0d0*roh(1)*roh(2)))

        alpha(2) = (-(hh(2)*hh(2)) + (roh(2)*roh(2)) + (roh(3)*roh(3)))
        alpha(2) =  checked_acos(alpha(2)/(2.0d0*roh(2)*roh(3)))

        alpha(3) = (-(hh(3)*hh(3)) + (roh(3)*roh(3)) + (roh(1)*roh(1)))
        alpha(3) = checked_acos(alpha(3)/(2.0d0*roh(3)*roh(1)))

!###################################################################
! Now determine the H3O+ intermolecular energy for this Eigenstate 
! from the bending and stretching terms.
!
!###################################################################
        vh3o_intra_r = 0.0d0
        vh3o_intra_alpha = 0.0d0
         vh3o_intra_tmp = 0.0d0
      
          do j=1,3
           vh3o_intra_r = (1.0d00 - DEXP(-1.0d0*aoheq*(roh(j)-roheq)))
           vh3o_intra_r = doh * (vh3o_intra_r*vh3o_intra_r)
           vh3o_intra_alpha = alpha(j)-alphaeq
           vh3o_intra_alpha = 0.5d0*kalpha*(vh3o_intra_alpha* vh3o_intra_alpha)
           vh3o_intra_tmp = vh3o_intra_tmp + vh3o_intra_r + vh3o_intra_alpha
        enddo
!###################################################################
! Finally assign the vh3o_intra for eigenstate (i) |i>
!
!###################################################################
        vh3o_intra(i) = vh3o_intra_tmp

!###################################################################
! End of loop over Eigenstates

 899     enddo

     RETURN
     END 

!###################################################################
!############### SUBROUTINE VINTER #################################
!###################################################################

     subroutine vinter(v_inter)
     use commons
     use msevb_common

     implicit none

!###################################################################
! Calculation of the Intermolecular interaction between H3O+ and H2O
!
! Coulomb, Repulsion and Lennard-Jones Terms
!
!###################################################################

     DOUBLE PRECISION V_INTER(REDUCED_NUM_EIG)
     DOUBLE PRECISION DX,DY,DZ,ROO,DXQQ,DYQQ,DZQQ,RQQ,QII,QJJ,ROO_6,ROO_12
     DOUBLE PRECISION LJ_TERM,REPULSE,COULOMB
     integer i,ii,jj,kk,ll,iii,jjj
     DOUBLE PRECISION LJ_TMP,REPULSE_TMP,COULOMB_TMP,QTIP3P_HY,QTIP3P_OX
     character*1 dummy

!###################################################################
! Start of loop over Eigenstates

     do 599 ll = 1,reduced_num_eig
        lj_term = 0.0d00
!        lj_tmp =0.0d0
        coulomb = 0.0d0
!        coulomb_tmp = 0.0d0
        repulse = 0.0d0
!        repulse_tmp = 0.0d0

!###################################################################
! Determination of the LJ and Replusive terms
!
!###################################################################
     
        ii = atmpl(ll,3)     ! The hydronium O atom in this eigenstate

!     if (QUICK) then
        do i = 6,(natoms-1),3
           jj = atmpl(ll,i)        ! The current H2O atom
           repulse = repulse + repulse_inter(ii,jj)
           lj_term = lj_term + lj_inter(ii,jj)
        enddo

!     else
!
!     do i = 6,(natoms-1),3 
!        dx = xx(ll,3) - xx(ll,i)
!         dy = yy(ll,3) - yy(ll,i)
!        dz = zz(ll,3) - zz(ll,i)
!        Roo = DSQRT(dx*dx + dy*dy + dz*dz)
!        repulse_tmp = big_b*(1.0d0-DTANH(small_b*(Roo - dooeq)))
!        Roo = sigma_mix/Roo
!         Roo = Roo * Roo
!         Roo_6 = Roo**6 
!         Roo_12 = Roo**12
!        lj_tmp = epsilon_mix*(Roo_12-Roo_6)
!        lj_term = lj_term + lj_tmp
!        repulse = repulse + repulse_tmp
!     enddo
!
!     endif

!###################################################################
! Determination of the Coulombic terms
!
!###################################################################

!      if (QUICK) then
            do ii = 1,4
             iii = atmpl(ll,ii)
             do jj = 5,natoms
                jjj = atmpl(ll,jj)
                coulomb = coulomb + inter_coulomb(iii,jjj)
             enddo
          enddo
 
!      else
!
!     do ii = 1,4
!        if (ii.eq.3) then
!            qii = qintero
!        else
!            qii = qinterh
!        endif
!        do jj = 5,natoms
!           dxqq = xx(ll,ii) - xx(ll,jj) 
!           dyqq = yy(ll,ii) - yy(ll,jj)
!           dzqq = zz(ll,ii) - zz(ll,jj)
!           Rqq = DSQRT(dxqq*dxqq+dyqq*dyqq+dzqq*dzqq)
!           if (mod(jj,3).eq.0) then
!             qjj = qtip3p_o
!           else
!             qjj = qtip3p_h
!           endif
!           coulomb_tmp = dampchge*(fourpieo*qii*qjj)/(Rqq)
!           coulomb = coulomb + coulomb_tmp
!        enddo
!     enddo
!
!      endif
!
!###################################################################
! Work out the v_inter() term for this Eigenstate from the sum of 
! coulomb, LJ and  repulsive terms

     v_inter(ll) = lj_term + coulomb + repulse

!     print *, 'Eigenstate:', ll
!     print *, 'Hydronium O:', atmpl(ll,3)
!     print *, 'LJ interaction:', lj_term
!     print *, 'Coulomb attraction:', coulomb
!     print *, 'Repulsive term:', repulse

!     LJ_values(ll)=lj_term
!     coulomb_values(ll)=coulomb
!     rep_values(ll)=repulse

!###################################################################
! End of loop over Eigenstates

 599     enddo

     RETURN
     END 

!#######################################################################
!################# SUBROUTINE BILDHAM ##################################
!####################################################################### 

     subroutine bildham
     use commons
     use msevb_common

     implicit none

!###################################################################
! Form the numerical Hamiltonian formed by the nuclear coordinates
!
!###################################################################

     integer i,j
     DOUBLE PRECISION VH2O_INTER(REDUCED_NUM_EIG),VH3O_INTRA(REDUCED_NUM_EIG)
        DOUBLE PRECISION VH2O_INTRA(REDUCED_NUM_EIG),V_INTER(REDUCED_NUM_EIG)
      
!###################################################################
! From only the nuclear coordinate vectors, xx(i,j) etc, calculate
! the diagonal and then the offdiagonal Hamiltonian elements
!
! Pass back from these routines the array of diagonal element 
! components.
!
!###################################################################

     do i = 1,reduced_num_eig
           vh3o_intra(i) = 0.0d0
        vh2o_intra(i) = 0.0d0
        v_inter(i) = 0.0d0
           vh2o_inter(i) = 0.0d0
     enddo

     call vh3ointra(vh3o_intra)
     call vh2ointra(vh2o_intra)

! Calculate the various energy components that are needed

     call calcAtomicInteractions
              
!----is there more than one water molecule?
        if (num_eig.gt.2) call vh2ointer(vh2o_inter)

     call voffdiag2
     call vinter(v_inter)

!###################################################################
! Build the Hamiltonian by using the element components just 
! determined
!
!###################################################################

     do i = 1,reduced_num_eig
         ham(i,i) = vh3o_intra(i) + vh2o_intra(i) + v_inter(i) + vh2o_inter(i)
        if (i.ne.reduced_num_eig) then
           do j = i+1,reduced_num_eig
              ham(i,j) = aij(i,j)
              ham(j,i) = aij(i,j)
           enddo
        endif

! Debugging

!        h3o_intra_energy(i)=vh3o_intra(i)
!        h2o_intra_energy(i)=vh2o_intra(i)
!        h3o_h2o_energy(i)=v_inter(i)
!        h2o_h2o_energy(i)=vh2o_inter(i)

     enddo

!     do i=1, reduced_num_eig-1
!        do j=i+1, reduced_num_eig
!           off_diag_energies(i,j)=ham(i,j)
!        enddo
!     enddo

     RETURN
     END

!#####################################################################

! Separate out the calculation of the energy terms because we are going reuse them
! in the force calculations

     subroutine calcAtomicInteractions ()

     use commons
     use msevb_common

     implicit none

! Local variables

     integer :: i, j
     DOUBLE PRECISION :: COULOMBPREFACTOR
     DOUBLE PRECISION :: QEXCH1, QWATER1, QH3O1, QEXCH2, QWATER2, QH3O2
     DOUBLE PRECISION :: R6,R12,R

     DOUBLE PRECISION, PARAMETER :: SIGMA6  = SIGMA_MIX**6
        DOUBLE PRECISION, PARAMETER :: TIP3P6  = H2OINTERSIGMA**6
     DOUBLE PRECISION, PARAMETER :: SIGMA12 = SIGMA6*SIGMA6
        DOUBLE PRECISION, PARAMETER :: TIP3P12 = TIP3P6*TIP3P6

! Calculate Coulombic interactions

     atom_coulomb = 0.0d0
     each_coulomb = 0.0d0
     water_inter_coulomb = 0.0d0
     inter_coulomb = 0.0d0
     
     do i = 1, natoms-1

        if (MOD(i,3).eq.0) then
           qexch1 = qexcho
           qwater1 = qtip3p_o
           qh3o1 = qintero
        else
           qexch1 = qexchh
           qwater1 = qtip3p_h
           qh3o1 = qinterh
        endif

        do j = i+1, natoms

           coulombPrefactor = fourpieo/interAtomicR(i,j)

           if (MOD(j,3).eq.0) then
           qexch2 = qexcho
           qwater2 = qtip3p_o
           qh3o2 = qintero
           else
           qexch2 = qexchh
           qwater2 = qtip3p_h
           qh3o2 = qinterh
           endif

           each_coulomb(i,j) = qexch1*qwater2*coulombPrefactor
           atom_coulomb(i) = atom_coulomb(i) + each_coulomb(i,j)
           each_coulomb(j,i) = qexch2*qwater1*coulombPrefactor
           atom_coulomb(j) = atom_coulomb(j) + each_coulomb(j,i)

           water_inter_coulomb(i,j) = qwater1*qwater2*coulombPrefactor
           water_inter_coulomb(j,i) = water_inter_coulomb(i,j)

           inter_coulomb(i,j) = dampchge*qh3o1*qwater2*coulombPrefactor
           inter_coulomb(j,i) = dampchge*qh3o2*qwater1*coulombPrefactor      
        enddo
     enddo

! O-O terms

     do i = 3, natoms-4, 3             ! Loop over O atoms
        do j = i+3, natoms-1, 3 
           R = interAtomicR(i,j) 
           r6 = 1.0d0/R**6
           r12 = r6*r6

! H2O-H2O O-O LJ term
           ljr(i,j) = 4.0d0*h2ointerEpsilon*((tip3p12*r12) - (tip3p6*r6))
           ljr(j,i) = ljr(i,j)
           
! H2O-H3O LJ term
           lj_inter(i,j) = epsilon_mix*((sigma12*r12) - (sigma6*r6))
           lj_inter(j,i) = lj_inter(i,j)

! H2O-H3O repulsion term
           repulse_inter(i,j) = big_b*(1.0d0-DTANH(small_b*(R - dooeq)))  
           repulse_inter(j,i) = repulse_inter(i,j)
        enddo
     enddo

     return

     end subroutine

!#####################################################################
!################# SUBROUTINE JACOBI #################################
!#####################################################################
     
      SUBROUTINE jacobi(a,n,d,v,nrot)
     implicit none

      INTEGER n,nrot,NMAX
      DOUBLE PRECISION A(N,N),D(N),V(N,N)
      PARAMETER (NMAX=500)
      INTEGER i,ip,iq,j,rac,mv
      DOUBLE PRECISION C,G,H,S,SM,T,TAU,THETA,TRESH,B(NMAX),Z(NMAX)

      do 12 ip=1,n
        do 11 iq=1,n
          v(ip,iq)=0.0d0
11      continue
        v(ip,ip)=1.0d0
12    continue
      do 13 ip=1,n
        b(ip)=a(ip,ip)
        d(ip)=b(ip)
        z(ip)=0.0d0
13    continue
      nrot=0
      do 24 i=1,50
        sm=0.0d0
        do 15 ip=1,n-1
          do 14 iq=ip+1,n
            sm=sm+abs(a(ip,iq))
14        continue
15      continue
        if(sm.eq.0.0d0)return
        if(i.lt.4)then
          tresh=0.2d0*sm/n**2
        else
          tresh=0.d0
        endif
        do 22 ip=1,n-1
          do 21 iq=ip+1,n
            g=100.d0*abs(a(ip,iq))
            if((i.gt.4).and.(abs(d(ip))+g.eq.abs(d(ip))).and.(abs(d(iq))+g.eq.abs(d(iq))))then
              a(ip,iq)=0.d0
            else if(abs(a(ip,iq)).gt.tresh)then
              h=d(iq)-d(ip)
              if(abs(h)+g.eq.abs(h))then
                t=a(ip,iq)/h
              else
                theta=0.5d0*h/a(ip,iq)
                t=1./(abs(theta)+sqrt(1.+theta**2))
                if(theta.lt.0.d0)t=-t
              endif
              c=1.d0/sqrt(1+t**2)
              s=t*c
              tau=s/(1.d0+c)
              h=t*a(ip,iq)
              z(ip)=z(ip)-h
              z(iq)=z(iq)+h
              d(ip)=d(ip)-h
              d(iq)=d(iq)+h
              a(ip,iq)=0.d0
              do 16 j=1,ip-1
                g=a(j,ip)
                h=a(j,iq)
                a(j,ip)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
16            continue
              do 17 j=ip+1,iq-1
                g=a(ip,j)
                h=a(j,iq)
                a(ip,j)=g-s*(h+g*tau)
                a(j,iq)=h+s*(g-h*tau)
17            continue
              do 18 j=iq+1,n
                g=a(ip,j)
                h=a(iq,j)
                a(ip,j)=g-s*(h+g*tau)
                a(iq,j)=h+s*(g-h*tau)
18            continue
              do 19 j=1,n
                g=v(j,ip)
                h=v(j,iq)
                v(j,ip)=g-s*(h+g*tau)
                v(j,iq)=h+s*(g-h*tau)
19            continue
              nrot=nrot+1
            endif
21        continue
22      continue
        do 23 ip=1,n
          b(ip)=b(ip)+z(ip)
          d(ip)=b(ip)
          z(ip)=0.
23      continue
24    continue
      print *, 'too many iterations in jacobi -- aborting '
     STOP

      return
      END

!#####################################################################
!################## SUBROUTINE VOFFDIAG ##############################
!#####################################################################

     subroutine voffdiag2
     use commons
     use msevb_common

     implicit none

!================================================================
!  r.christie   10/11/98
!
! DETERMINATION OF THE OFF-DIAGONAL ELEMENTS OF THE HAMILTONIAN
!
!================================================================

     integer i,j
     DOUBLE PRECISION ROH,ROO
     DOUBLE PRECISION ASYMM_COORD,F1,F2,F3

! Global storage of the Zundel terms for reuse in the force calculations

     zundel_f = 0.0d0
     zundel_g = 0.0d0

     do i = 1,reduced_num_eig
        do j = 1,reduced_num_eig
           vijexch(i,j) = 0.0d0
        enddo
     enddo

     if (num_eig.gt.2) call exchange2

!###################################################################
! Calculate simple Molecular Mechanics functions in off diagonal
! term. Only two variables: the asymmetric stretch coordinate, and
! (ii) the R(OO) vector in the H5O2+ dimer considered as 
! participating in the exchange charge interaction.

     do i = 1,(reduced_num_eig-1)
          do j = (i+1),reduced_num_eig

           if (statesInteract(i,j)) then

!           print *, 'States', i, 'and', j, 'interact'

! Already assigned a Zundel species for the interaction of these two 
! states 

           roh = interAtomicR(zundel_species(i,j,4), atmpl(i,3))
                                                                                       
!----Determine the R(OO) separation between the CENTRAL hydronium
!    in state |i> and that in state |2>, as is necessary in the
!    exchange term

           roo = interAtomicR(atmpl(i,3), atmpl(j,3))

!##################################################################
! Two versions of determining the asymmetric stretch Coordinate
! ----------------------------------------------------------------
! 1. Udo Schmitt Pre-Print version of calculating the asymmetric
!    stretch                                                    
!----determine the unit vector eoo
!             eoox = (1.0d0/roo)*oodx
!            eooy = (1.0d0/roo)*oody
!            eooz = (1.0d0/roo)*oodz
!----determine the asymmetric stretch terms qi
!             qx = ((oodx*0.5d0) - dx)*eoox 
!             qy = ((oody*0.5d0) - dy)*eooy
!             qz = ((oodz*0.5d0) - dz)*eooz
!----determine the magnitude of the asymmetric stretch
!             asymm_coord = (qx + qy + qz)
!-----------------------------------------------------------------
! 2. Plain old method for determining the asymmetric stretch
!    Coordinate                                              
           asymm_coord = 0.5d0*roo - roh
!#################################################################          

           f1 = (1.00d0 + P*DEXP(-1.0d00*kexp*(roo-doo)*(roo-doo)))
           f2 = 0.50d0*(1.0d0-DTANH(beta*(roo - big_rooeq )))
           f3 = 10.0d0*DEXP(-1.0d00*alpha_bleh*(roo - rooeq))
           zundel_f(i,j) = f1*(f2+f3)
           zundel_g(i,j) = DEXP(-1.0d00*gamma_msevb*asymm_coord*asymm_coord)

!           print *, 'State 1:', i, 'State 2:', j
!           print *, 'Components of the off diagonal contribution to the energy:'
!           print *, 'g:', g
!           print *, 'f=f1(f2+f3)'
!           print *, 'f1:', f1
!           print *, 'f2:', f2
!           print *, 'f3:', f3
!           print *, 'f:', f
!           print *, 'Roo:', roo
!           print *, 'Roh:', roh
!           print *, 'q:', asymm_coord
!           print *, 'vijexch:', vijexch(i,j)

!----assign the matrix element
           aij(i,j) = (vij+vijexch(i,j))*zundel_f(i,j)*zundel_g(i,j)
!----hermetian matrix mind
           aij(j,i) = aij(i,j)

           else
           aij(i,j) = 0.0d0
           aij(j,i) = 0.0d0
           endif
         enddo
      enddo

     return 
     end

! ##################################################################################################

        subroutine exchange2

     use msevb_common

     implicit none
!---------------------------------------------------------------------
!  calculation of the interaction of the exchange charge distribution 
!  for the H5O2+ dimers resulting from <i|j> with the rest of the 
!  TIP3P charges on the remaining water molecules
!
!---------------------------------------------------------------------

! At this point atom_coulomb(i) holds the sum of exchange interactions with all other atoms,
! where i holds the exchange charge and the other atoms are assigned TIP3 charges
! Need to therefore subtract the interactions with the other members of the hydronium species

     integer i,j,k,l

     do i = 1,(reduced_num_eig-1)
        do j = (i+1),reduced_num_eig

           vijexch(i,j) = 0.0D0

           if (statesInteract(i,j)) then

           do k=1,7

! Sum of exchange charge interactions with all other atoms
              vijexch(i,j) = vijexch(i,j) + atom_coulomb(zundel_species(i,j,k))

! Subtract the interactions with the other Zundel atoms
              do l=1,7
                 if (l.eq.k) cycle
                 vijexch(i,j) = vijexch(i,j) - each_coulomb(zundel_species(i,j,k),zundel_species(i,j,l))
              enddo
           enddo
           endif
        enddo
     enddo

        RETURN
     END

!#############################################################################################

! J Chem Phys 2002 version of this algorithm
! Including the H-bond angle cutoff that they say doesn`t make any difference
! 04/12/2003 Modified (corrected) so that each O atom can now form part of more than one
! hydronium species

! Experimental min roh rather than roh^2 code

        subroutine assignvb3

     use commons
     use msevb_common
     use msevb_interfaces

     implicit none

     DOUBLE PRECISION MINHBONDANGLERAD,HBONDANGLE

     DOUBLE PRECISION MINSUMROH,CURRENTSUMROH     
     logical :: assignedPivotState

     integer :: generatorState(maxNumVBstates)  ! VB state i was generated from generatorState(i)

! Holds the closest H atoms for the Os and the closest O atoms for the Hs

     integer proximalAtoms(natoms,num_hyd)   
     integer neighbourCounts(natoms)

     integer currentO, currentH

     integer firstVBstate, secondVBstate, currentShellHcount, nextShellHcount, currentShell
     integer currentExchangeH, newEigenO, oldVBstate
     integer, dimension(:,:), allocatable :: exchangeHs
     integer, allocatable :: fingerprints(:,:)
     logical newFingerprint
     integer O1pos,O2pos,centralHpos
     integer firstStateO, secondStateO

! Loop variables

     integer i,j,k,l,eigen_O

! Convert minimum H bond angle to radians

     minHbondAngleRad = (minHbondAngle/180.0d0)*pi

! Order the O atoms in terms of the distance from each H, only consider those which are less than
! the cutoff distance apart

     neighbourCounts = 0

     do i = 1,num_hyd

        currentH = i + INT((i-1)/2)

        do j = 1,num_eig
           currentO = 3*j
           if (interAtomicR(currentO,currentH).lt.maxHbondLength) then
           neighbourCounts(currentH) = neighbourCounts(currentH) + 1

           do k=1,neighbourCounts(currentH)
              if (k.eq.neighbourCounts(currentH)) then
                 proximalAtoms(currentH,k) = currentO
                 exit
              elseif (interAtomicR(currentO,currentH).lt.interAtomicR(proximalAtoms(currentH,k),currentH)) then
                 do l=neighbourCounts(currentH),k+1,-1
                 proximalAtoms(currentH,l)=proximalAtoms(currentH,l-1)
                 enddo
                 proximalAtoms(currentH,k) = currentO
                 exit
              endif
           enddo
           endif
        enddo

! Check whether we have assigned at least one O 

        if (neighbourCounts(currentH).eq.0) then
           print *, 'Within the cutoff of', maxHbondLength, 'program assigned zero O atoms to H', currentH
           print *, 'Each H must be assigned at least one O'
           print *, 'Increase cutoff or examine your geometry'

           do j=1, NATOMS
           if (MOD(j,3).eq.0) then
              print *, 'O', psix(j), psiy(j), psiz(j)
           else
              print *, 'H', psix(j), psiy(j), psiz(j)
           endif
           enddo
           stop
        endif

     enddo

!     print *, 'Assigned Os:'

!     do i=1, num_hyd
        
!        print *, 'H:', Hposition(i)
!        do j=1, assignedOs(i)
!           print *, j, 3*neighbouringOs(i,j)
!        enddo
!     enddo

! Try a simple VB assignment first, if that doesn`t work try the iterative procedure

     assignedPivotState = .FALSE.
     reduced_num_eig = 0

     call trySimpleVBAssignment(proximalAtoms, neighbourCounts, reduced_num_eig, atmpl(1:2,:))

     if (reduced_num_eig.eq.0) then

        call assignPivotVBstate(proximalAtoms, neighbourCounts, assignedPivotState, atmpl(1,:))

        if (assignedPivotState) then
           reduced_num_eig = 1
        else
           print *, 'Unable to assign pivot VB state'
           print *, 'System:'
           do i = 1, natoms
           if (MOD(i,3).eq.0) then
              print *, 'O', psix(i), psiy(i), psiz(i)
           else
              print *, 'H', psix(i), psiy(i), psiz(i)
           endif
           enddo
           stop   
        endif
     endif

     generatorState(1:reduced_num_eig) = -1    ! -1 indicates that the states are pivot states
     
!     print *, 'Number of pivot states assigned:', reduced_num_eig
     
!     do i=1,reduced_num_eig
!        print *, '*** State', i, '***'
!        do j=1,natoms
!           print *, j, atmpl(i,j)
!        enddo
!     enddo

!     stop

     allocate(fingerprints(maxNumVBstates,natoms))
     fingerprints = 0

     do i=1,reduced_num_eig

        fingerprints(i,atmpl(i,1)) = atmpl(i,3)
        fingerprints(i,atmpl(i,2)) = atmpl(i,3)
        fingerprints(i,atmpl(i,4)) = atmpl(i,3)
        
        do j=6, (3*num_eig), 3
           fingerprints(i,atmpl(i,j-1))=atmpl(i,j)
           fingerprints(i,atmpl(i,j+1))=atmpl(i,j)
        enddo
     enddo

! Initialise the VB state interaction matrix

     statesInteract = .FALSE.

! Allocate the exchange Hs array

     if (shellsToCount.lt.2) then
        ALLOCATE(exchangeHs((reduced_num_eig+2),2))
     else
        ALLOCATE(exchangeHs((reduced_num_eig+2)*3*(2**(shellsToCount-2)), 2))
     endif

! Run through the hydration shells

     if (reduced_num_eig.gt.1) then

        currentShellHcount = 0

        FirstPivotState: do i = 1, 4

           if (i.ne.3) then

! For second and subsequent pivot states the fourth entry is the symmetric H
           
           do secondVBstate = 2, reduced_num_eig

              if (atmpl(secondVBstate,4).eq.atmpl(1,i)) then

                 call assignZundelSpecies(1,secondVBstate,atmpl(1,i),.false.,minHbondAngleRad)

                 exchangeHs(currentShellHcount+1,1) = atmpl(secondVBstate,1)
                 exchangeHs(currentShellHcount+1,2) = secondVBstate
                 exchangeHs(currentShellHcount+2,1) = atmpl(secondVBstate,2)
                 exchangeHs(currentShellHcount+2,2) = secondVBstate

                 currentShellHcount = currentShellHcount+2
                 cycle FirstPivotState
              endif
           enddo

           currentShellHcount = currentShellHcount+1
           exchangeHs(currentShellHcount,1) = atmpl(1,i)
           exchangeHs(currentShellHcount,2) = 1           
           endif
        enddo FirstPivotState

     else
        currentShellHcount = 3
        exchangeHs(1,1) = atmpl(reduced_num_eig,1)
        exchangeHs(1,2) = reduced_num_eig
        exchangeHs(2,1) = atmpl(reduced_num_eig,2)
        exchangeHs(2,2) = reduced_num_eig
        exchangeHs(3,1) = atmpl(reduced_num_eig,4)
        exchangeHs(3,2) = reduced_num_eig
     endif
     
! Loop through hydration shells, 3 at present (set in msevb_common module)
! Only assign new VB states for the specified number of shells

     do currentShell = 1, shellsToCount

!          print *, 'Shell', currentShell, currentShellHcount

!          print *, 'Exchanging Hs in this shell'
!          do i = 1, currentShellHcount
!             print *, i, exchangeHs(i,1)
!          enddo

        nextShellHcount = 0

        currentShellLoop: do i = 1, currentShellHcount
           
           currentExchangeH = exchangeHs(i,1)

           if (neighbourCounts(currentExchangeH).gt.1) then

!             print *, 'Currently considering exchange H', currentExchangeH
!             print *, 'Assigned O count:', neighbourCounts(currentExchangeH)

! At most, each H can be assigned to two different hydronium species therefore check the closest two Os
! If the H was only assigned one O then this hydronium has already been taken care of
        
           firstVBstate = exchangeHs(i,2)

           if (atmpl(firstVBstate,3).eq.proximalAtoms(currentExchangeH,1)) then
              newEigenO = proximalAtoms(currentExchangeH,2)
           else
              newEigenO = proximalAtoms(currentExchangeH,1)
           endif

! Check the H bond angle

           O2pos = atmpl(firstVBstate,3)

           hBondAngle = (psix(newEigenO)-psix(currentExchangeH))*(psix(O2pos)-psix(currentExchangeH))+& 
                             (psiy(newEigenO)-psiy(currentExchangeH))*(psiy(O2pos)-psiy(currentExchangeH))+&
                          (psiz(newEigenO)-psiz(currentExchangeH))*(psiz(O2pos)-psiz(currentExchangeH))


           hbondAngle = DACOS(hbondAngle/(interAtomicR(newEigenO,currentExchangeH) &
                              * interAtomicR(O2pos,currentExchangeH)))
              
!              print *, 'Bond angle: ', O1pos, centralHpos, O2pos, ((hbondAngle/pi)*180d0)
              
           if (hbondAngle.ge.minHbondAngleRad) then

              reduced_num_eig = reduced_num_eig + 1

! Check whether we have seen this VB state before - can happen with cycles

              fingerprints(reduced_num_eig,:)=fingerprints(firstVBstate,:)
              fingerprints(reduced_num_eig,currentExchangeH)=newEigenO

              do j=1,reduced_num_eig-1
                 if (j.ne.firstVBstate) then

                 newFingerprint = .FALSE.

                 do k=1, natoms
                    if (fingerprints(j,k).ne.fingerprints(reduced_num_eig,k)) then
                    newFingerprint = .TRUE.
                    exit
                    endif
                 enddo

                 if (.not.newFingerprint) then
!                    print *, 'Seen this fingerprint before', reduced_num_eig, j
!                    do k=1, natoms
!                    print *, k, atmpl(j,k)
!                    enddo
                    reduced_num_eig = reduced_num_eig - 1
                    cycle currentShellLoop                   
                 endif
                 endif
              enddo
              
              atmpl(reduced_num_eig, 3) = newEigenO
              atmpl(reduced_num_eig, 4) = currentExchangeH                 

!                   print *, 'New eigen O:', newEigenO

              do j=6, natoms, 3
                 if (atmpl(firstVBstate,j).eq.newEigenO) then
                 atmpl(reduced_num_eig,1) = atmpl(firstVBstate,j-1)
                 atmpl(reduced_num_eig,2) = atmpl(firstVBstate,j+1)

! Add these protons to the potential exchange list for the next shell

                 if (currentShell.ne.shellsToCount) then

                    nextShellHcount = nextShellHcount + 1
                    if ((currentShellHcount+nextShellHcount).gt.SIZE(exchangeHs,1)) then
                    print *, 'Too many Hs in exchange H list'
                    print *, 'Increase size of array'
                    do k = 1, natoms
                       if (MOD(k,3).eq.0) then
                          print *, 'O', psix(k), psiy(k), psiz(k)
                       else
                          print *, 'H', psix(k), psiy(k), psiz(k)
                       endif
                    enddo
                    print *, 'Size of array:', SIZE(exchangeHs,1)
                    print *, 'Current shell:', currentShell
                    do k=1,currentShellHcount
                       print *, k, exchangeHs(k,1), exchangeHs(k,2)
                    enddo
                    print *, 'Next shell:', (currentShell+1)
                    do k=currentShellHcount+1, currentShellHcount+nextShellHcount-1
                       print *, k, exchangeHs(k,1), exchangeHs(k,2)
                    enddo
                    stop
                    endif
                    exchangeHs(currentShellHcount+nextShellHcount,1) = atmpl(reduced_num_eig,1)
                    exchangeHs(currentShellHcount+nextShellHcount,2) = reduced_num_eig
                 
                    nextShellHcount = nextShellHcount + 1
                    if ((currentShellHcount+nextShellHcount).gt.SIZE(exchangeHs,1)) then
                    print *, 'Too many Hs in exchange H list'
                    print *, 'Increase size of array'
                    do k = 1, natoms
                       if (MOD(k,3).eq.0) then
                          print *, 'O', psix(k), psiy(k), psiz(k)
                       else
                          print *, 'H', psix(k), psiy(k), psiz(k)
                       endif
                    enddo
                    print *, 'Size of array:', SIZE(exchangeHs,1)
                    print *, 'Current shell:', currentShell
                    do k=1,currentShellHcount
                       print *, k, exchangeHs(k,1), exchangeHs(k,2)
                    enddo
                    print *, 'Next shell:', (currentShell+1)
                    do k=currentShellHcount+1, currentShellHcount+nextShellHcount-1
                       print *, k, exchangeHs(k,1), exchangeHs(k,2)
                    enddo
                    stop
                    endif
                    exchangeHs(currentShellHcount+nextShellHcount,1) = atmpl(reduced_num_eig,2)
                    exchangeHs(currentShellHcount+nextShellHcount,2) = reduced_num_eig
                 endif
                 exit
                 endif
              enddo

! Add the old hydronium species, minus the exchange proton

              atmpl(reduced_num_eig,6) = atmpl(firstVBstate, 3)
           
              if (atmpl(firstVBstate,1).eq.currentExchangeH) then
                 atmpl(reduced_num_eig,5) = atmpl(firstVBstate,2)
                 atmpl(reduced_num_eig,7) = atmpl(firstVBstate,4)
              else
                 atmpl(reduced_num_eig,5) = atmpl(firstVBstate,1)            
                 
                 if (atmpl(firstVBstate,2).eq.currentExchangeH) then
                 atmpl(reduced_num_eig,7) = atmpl(firstVBstate,4)                    
                 else
                 atmpl(reduced_num_eig,7) = atmpl(firstVBstate,2)     
                 endif
              endif
              
! Copy across the remainder of the water species from the old VB state, these will not have changed

              k = 9

              do j = 6, natoms-1, 3
                 if (atmpl(firstVBstate, j).eq.atmpl(reduced_num_eig,3)) cycle
              
                 atmpl(reduced_num_eig, k-1) = atmpl(firstVBstate,j-1)
                 atmpl(reduced_num_eig, k) = atmpl(firstVBstate,j)
                 atmpl(reduced_num_eig, k+1) = atmpl(firstVBstate,j+1)
              
                 k = k + 3
              enddo

!              if (.not.newFingerprint) then
!                 print *, 'New VB state:'
!                 do k= 1, natoms
!                 print *, k, atmpl(reduced_num_eig,k)
!                 enddo
!              endif

! Assign the Zundel species

              call assignZundelSpecies (firstVBstate,reduced_num_eig,currentExchangeH,.TRUE.,minHbondAngleRad)

 !                 print *, 'New zundel species:'
 !                 do k = 1, 7
 !                 print *, k, zundel_species(firstVBstate, reduced_num_eig, k)
 !                 enddo

              generatorState(reduced_num_eig) = firstVBstate

! Check whether any other VB states use this first VB state O as the hydronium with currentExchangeH as one of the Hs
! Assign interactions if they do

              do j = 1, (reduced_num_eig-1)
                 if (j.ne.firstVBstate .and. atmpl(j,3).eq.O2pos .and. &
                            (currentExchangeH.eq.atmpl(j,1) & 
                            .or. currentExchangeH.eq.atmpl(j,2) .or. currentExchangeH.eq.atmpl(j,4))) then

                 call assignZundelSpecies(j,reduced_num_eig,currentExchangeH,.TRUE.,minHbondAngleRad)

!                     print *, 'State', j, 'also interacts with this new state'
!                     print *, 'New zundel species:'
!                     do k = 1, 7
!                     print *, k, zundel_species(j, reduced_num_eig, k)
!                     enddo               
                 endif
              enddo
              
           endif
 
           endif

        enddo currentShellLoop ! End of loop over current shell

! Update the exchange H list

        do i = 1, nextShellHcount
           exchangeHs(i,1) = exchangeHs(i+currentShellHcount,1)
           exchangeHs(i,2) = exchangeHs(i+currentShellHcount,2)      
        enddo
       
        currentShellHcount = nextShellHcount

     enddo

! Deallocate memory

     DEALLOCATE(exchangeHs,fingerprints)

! Check for unassigned interactions - this bit works in tandem with the in loop checking
! The in loop checking takes care of extra interactions with states which were assigned before the current one
! This bit takes care of extra interactions with states which were assigned after the current one
! See journal entry 06/04/04

     do i = 2, reduced_num_eig-1
        if (generatorState(i).ne.-1) then

           do j = i+1, reduced_num_eig
           if (j.ne.generatorState(i) .and. (.not.statesInteract(i,j)) .and. &
                     atmpl(j,3).eq.atmpl(generatorState(i),3)) then

!              print *, 'Missed interaction between states', i, 'and', j

! atmpl(i,4) holds the exchange H with which this state was initially generated (except for the pivot state(s))
              call assignZundelSpecies (i,j,atmpl(i,4),.TRUE.,minHbondAngleRad)

!              print *, 'New Zundel species'
!              do k=1,7
!                 print *, k, zundel_species(i,j,k)
!              enddo

           endif
           enddo
        endif
     enddo

     if (reduced_num_eig.gt.maxNumVBstates) then
        print *, 'Number of VB states generated exceeds hard limit'
        stop
     endif

!      print *, 'VB states assigned', reduced_num_eig

!       do i=1, reduced_num_eig

!          print *, '*** State', i, '***'

!          do j=1, natoms
!             print *, j, atmpl(i,j)
!          enddo
!       enddo

!     print *, 'Exchanging H list'

!     do i=1, (currentShellStartH + currentShellHcount - 1)
!        print *, i, Hposition(exchangeHs(i))
!     enddo

!       print *, 'Assigned Zundel species'

!       do i = 1, (reduced_num_eig-1)
!          do j = (i+1), reduced_num_eig
!             if (statesInteract(i,j)) then
!             print *, 'State 1:', i, 'State 2:', j
!             do k=1,7
!                print *, k, zundel_species(i,j,k)
!             enddo
!             endif
!          enddo
!       enddo

!     stop

     end subroutine

! ######################################################################################

! Assigns the exchange species between two hydronium species, checking the H bond angle
! if necessary

     subroutine assignZundelSpecies(firstVBstate, secondVBstate, exchangeH, checkedHbond, minHbondAngleRad)

     use commons
     use msevb_common

     implicit none

     integer, intent(IN) :: firstVBstate, secondVBstate, exchangeH
     logical, intent(IN) :: checkedHbond
     DOUBLE PRECISION, INTENT(IN) :: MINHBONDANGLERAD
     
     integer :: O1pos, O2pos, centralHpos
     DOUBLE PRECISION :: HBONDANGLE
     integer :: lowerIndex, higherIndex

! Make sure that the first VB state has the lower index - matrices are upper-triangular

     if (firstVBstate.lt.secondVBstate) then
        lowerIndex = firstVBstate
        higherIndex = secondVBstate
     else
        lowerIndex = secondVBstate
        higherIndex = firstVBstate
     endif

     if (.not.statesInteract(lowerIndex,higherIndex)) then

! Check the H bond angle if necessary

        if (.not.checkedHbond) then
           O1pos = atmpl(firstVBstate,3)
           O2pos = atmpl(secondVBstate,3)

           hBondAngle = (psix(O1pos)-psix(exchangeH))*(psix(O2pos)-psix(exchangeH)) + &
                          (psiy(O1pos)-psiy(exchangeH))*(psiy(O2pos)-psiy(exchangeH)) + &
                       (psiz(O1pos)-psiz(exchangeH))*(psiz(O2pos)-psiz(exchangeH))

           hbondAngle = DACOS(hbondAngle/(interAtomicR(O1pos,exchangeH)*interAtomicR(O2pos,exchangeH))) 
              
!              print *, 'Bond angle: ', O1pos, exchangeH, O2pos, ((hbondAngle/pi)*180d0)
              
           if (hbondAngle.lt.minHbondAngleRad) return
        endif

        zundel_species(lowerIndex,higherIndex,3) = atmpl(lowerIndex,3)  ! O
        zundel_species(lowerIndex,higherIndex,4) = exchangeH            ! Exchanging H
        zundel_species(lowerIndex,higherIndex,6) = atmpl(higherIndex,3) ! O
           
! Sort out the other protons from the first hydronium species

        if (atmpl(lowerIndex,1).eq.exchangeH) then
           zundel_species(lowerIndex,higherIndex,1) = atmpl(lowerIndex,2)         
           zundel_species(lowerIndex,higherIndex,2) = atmpl(lowerIndex,4)
        else
           zundel_species(lowerIndex,higherIndex,1) = atmpl(lowerIndex,1)
              
           if (atmpl(lowerIndex,2).eq.exchangeH) then
           zundel_species(lowerIndex,higherIndex,2) = atmpl(lowerIndex,4)
           else
           zundel_species(lowerIndex,higherIndex,2) = atmpl(lowerIndex,2)
           endif
        endif

! And those from the second hydronium species
              
        if (atmpl(higherIndex,1).eq.exchangeH) then
           zundel_species(lowerIndex,higherIndex,5) = atmpl(higherIndex,2)    
           zundel_species(lowerIndex,higherIndex,7) = atmpl(higherIndex,4)
        else
           zundel_species(lowerIndex,higherIndex,5) = atmpl(higherIndex,1)
           
           if (atmpl(higherIndex,2).eq.exchangeH) then
           zundel_species(lowerIndex,higherIndex,7) = atmpl(higherIndex,4)
           else
           zundel_species(lowerIndex,higherIndex,7) = atmpl(higherIndex,2)
           endif
        endif

! Set the interaction label

        statesInteract(lowerIndex,higherIndex) = .TRUE.

     endif

     return

     end subroutine

! ###############################################################################################

! Tries the simplest possible pivot VB assignment by assigning each H to the nearest O
! and checking whether or not that is viable
! Now considers the possibility that there could be two approximately equally good 'pivot' states
! in the case where the local geometry represents a Zundel species

     subroutine trySimpleVBAssignment(proximalAtoms, neighbourCounts, numStatesAssigned, pivotVBstates)

     use commons
     use msevb_common

     implicit none

! Subroutine arguments

     integer, intent(IN) :: proximalAtoms(natoms, num_hyd)
     integer, intent(IN) :: neighbourCounts(natoms)
     integer, intent(OUT) :: numStatesAssigned
     integer, intent(OUT) :: pivotVBstates(2,natoms)

! Local variables

     integer :: currentAtom, nearestOref, insertRef
     integer :: assignedHcounts(num_eig)
     integer :: assignedHs(num_eig,3)

     logical :: pivotOassigned

     logical :: symmetricH(natoms)
     integer :: checkH, sharedH, pivotO1ref, pivotO2ref
     DOUBLE PRECISION, PARAMETER :: EQUIDISTANCECUTOFF = 1.2D0
     DOUBLE PRECISION :: DISTANCERATIO

! See if we can identify a symmetric H atom

     symmetricH = .false.

     do currentAtom = 1, natoms
        if ((MOD(currentAtom,3).ne.0) .and. (neighbourCounts(currentAtom).gt.1)) then

           distanceRatio = interAtomicR(currentAtom,proximalAtoms(currentAtom,2))/ &
                                  interAtomicR(currentAtom,proximalAtoms(currentAtom,1))

           if (distanceRatio.le.equidistanceCutoff) symmetricH(currentAtom) = .true.
        endif
     enddo

! Find the nearest O atom to each H

     numStatesAssigned = 1
     pivotOassigned = .FALSE.
     assignedHcounts = 0
     pivotVBstates = 0

     do currentAtom = 1, natoms
        if (MOD(currentAtom,3).ne.0) then

           nearestOref = proximalAtoms(currentAtom,1)/3

           if (assignedHcounts(nearestOref).eq.3) then
           numStatesAssigned = 0
           exit
           else if (assignedHcounts(nearestOref).eq.2) then
           if (pivotOassigned) then
              numStatesAssigned = 0
              exit
           else
              assignedHcounts(nearestOref) = assignedHcounts(nearestOref) + 1
              assignedHs(nearestOref,assignedHcounts(nearestOref)) = currentAtom
              pivotOassigned = .TRUE.    
           endif                
           else
           assignedHcounts(nearestOref) = assignedHcounts(nearestOref) + 1
           assignedHs(nearestOref,assignedHcounts(nearestOref)) = currentAtom
           endif
        endif
     enddo

     if (numStatesAssigned.eq.1) then

! Fill the VB assignment

        insertRef = 6

        do currentAtom = 1, num_eig
           if (assignedHcounts(currentAtom).eq.3) then    ! H3O+ molecule
           pivotVBstates(numStatesAssigned,1) = assignedHs(currentAtom,1)
           pivotVBstates(numStatesAssigned,2) = assignedHs(currentAtom,2)
           pivotVBstates(numStatesAssigned,3) = currentAtom*3
           pivotVBstates(numStatesAssigned,4) = assignedHs(currentAtom,3)
           else    ! Normal H2O molecule
           pivotVBstates(numStatesAssigned,insertRef-1) = assignedHs(currentAtom,1)
           pivotVBstates(numStatesAssigned,insertRef) = currentAtom*3
           pivotVBstates(numStatesAssigned,insertRef+1) = assignedHs(currentAtom,2)

           insertRef = insertRef + 3
           endif
        enddo

! Check if there is another, almost equally good pivot state which needs assigning

        pivotO1ref = pivotVBstates(1,3)/3
        sharedH = -1

        do checkH = 1, 3
           if (symmetricH(assignedHs(pivotO1ref,checkH))) then
           if (sharedH.eq.-1) then
              sharedH = checkH
           else 
              sharedH = -1
              exit
           endif
           endif
        enddo
        
        if (sharedH.ne.-1) then

           pivotO2ref = proximalAtoms(assignedHs(pivotO1ref,sharedH),2)/3
           numStatesAssigned = numStatesAssigned + 1

! Assign another pivot state

           pivotVBstates(numStatesAssigned,1) = assignedHs(pivotO2ref,1)
           pivotVBstates(numStatesAssigned,2) = assignedHs(pivotO2ref,2)
           pivotVBstates(numStatesAssigned,3) = pivotO2ref*3
           pivotVBstates(numStatesAssigned,4) = assignedHs(pivotO1ref,sharedH)

           pivotVBstates(numStatesAssigned,6) = pivotVBstates(1,3)
           if (sharedH.eq.1) then
           pivotVBstates(numStatesAssigned,5) = assignedHs(pivotO1ref,2)
           pivotVBstates(numStatesAssigned,7) = assignedHs(pivotO1ref,3)
           else
           pivotVBstates(numStatesAssigned,5) = assignedHs(pivotO1ref,1)

           if (sharedH.eq.2) then
              pivotVBstates(numStatesAssigned,7) = assignedHs(pivotO1ref,3)
           else
              pivotVBstates(numStatesAssigned,7) = assignedHs(pivotO1ref,2)              
           endif
           endif

           insertRef = 9

           do currentAtom = 1, num_eig             
           if ((currentAtom.ne.pivotO1ref) .and. (currentAtom.ne.pivotO2ref)) then
              pivotVBstates(numStatesAssigned,insertRef-1) = assignedHs(currentAtom,1)
              pivotVBstates(numStatesAssigned,insertRef) = currentAtom*3
              pivotVBstates(numStatesAssigned,insertRef+1) = assignedHs(currentAtom,2)
              insertRef = insertRef + 3
           endif
           enddo
        endif

     endif

     return

     end subroutine

! ###############################################################################################

! Assign the pivot VB state by minimising the sum of bonded OH distances
! The pivot O is identified as the O with the 3 closest Hs

     subroutine assignPivotVBstate(proximalAtoms, neighbourCounts, success, pivotVBstate)

     use commons
     use msevb_common

     implicit none

! Subroutine arguments

     integer, intent(IN OUT) :: proximalAtoms(natoms, num_hyd)
     integer, intent(IN OUT) :: neighbourCounts(natoms)
     logical, intent(OUT) :: success
     integer, intent(OUT) :: pivotVBstate(natoms)

! Local variables

     integer :: i, j, k, l
     integer :: currentO, currentH, pivotOxygen
     integer :: assignedOcount
     integer :: Olist(num_eig)

     DOUBLE PRECISION :: MINSUMROH, CURRENTSUMROH

     logical :: usedH(natoms), goingUp, foundAnH
     integer :: Hcounters(num_eig,3)

     success = .TRUE.

! Order the H atoms in terms of the distance from each O, only consider those which are less than
! the cutoff distance apart

     do i=1,num_eig

        currentO = 3*i
        neighbourCounts(currentO) = 0   

        do j=1,num_hyd

           currentH = j + INT((j-1)/2)

           if (interAtomicR(currentO,currentH).lt.maxHbondLength) then
           neighbourCounts(currentO) = neighbourCounts(currentO) + 1

           do k=1,neighbourCounts(currentO)
              if (k.eq.neighbourCounts(currentO)) then
                 proximalAtoms(currentO,k) = currentH
                 exit
              elseif (interAtomicR(currentO,currentH).lt.interAtomicR(currentO,proximalAtoms(currentO,k))) then
                 do l=neighbourCounts(currentO),k+1,-1
                 proximalAtoms(currentO,l)=proximalAtoms(currentO,l-1)
                 enddo
                 proximalAtoms(currentO,k) = currentH
                 exit
              endif
           enddo
           endif
        enddo

! Check whether we have assigned sufficient Hs within the cutoff to form at least a water molecule

        if (neighbourCounts(currentO).lt.2) then
           print *, 'Within the cutoff of', maxHbondLength, 'program assigned', neighbourCounts(currentO), &
                       'as possibly bonded to O atom', currentO
           print *, 'Require a minimum of 2 assigned Hs to form a water species'
           print *, 'Increase cutoff or examine your geometry'
           success = .FALSE.
           return
        endif
     enddo

! Find the pivot oxygen, the O atom with the three closest Hs

     minSumRoh = HUGE(minSumRoh)
     pivotOxygen = -1

     do currentO = 3, (natoms-1), 3
        if (neighbourCounts(currentO).ge.3) then

           currentSumRoh = interAtomicR(currentO,proximalAtoms(currentO,1)) + &
                           interAtomicR(currentO,proximalAtoms(currentO,2)) + &
                              interAtomicR(currentO,proximalAtoms(currentO,3))

           if (currentSumRoh.lt.minSumRoh) then
           pivotOxygen = currentO
           minSumRoh = currentSumRoh
           endif
        endif
     enddo

! We must have assigned a pivot O because every H has been assigned at least one O
! Assign the rest of this VB state

     pivotVBstate = -1

! Order the other oxygens in decreasing proximity from the hydronium O
! This should be a more efficient way of finding the best overall bond assignment

     Olist = -1
     Olist(1) = pivotOxygen

     assignedOcount = 1

     do i=1,num_eig
        currentO = 3*i

        if (currentO.eq.pivotOxygen) cycle
           
        assignedOcount = assignedOcount + 1

        do j=2,assignedOcount
           
           if (Olist(j).eq.-1) then
           Olist(j) = currentO
           exit
           elseif (interAtomicR(pivotOxygen,currentO).lt.interAtomicR(pivotOxygen,Olist(j))) then
           
           do k=assignedOcount,j+1,-1
              Olist(k) = Olist(k-1)
           enddo

           Olist(j) = currentO
           exit                 
           endif
        enddo
     enddo

     usedH = .FALSE.

! Set the counters to their initial values

     Hcounters(1,1) = 1
     Hcounters(1,2) = 2
     Hcounters(1,3) = 3

     assignedOcount = 1
     goingUp = .TRUE.

! Assign the best VB species by minimising the sum of the distances between bonded atoms
 
! Three Hs will be associated with the pivot O and two with each of the others

     minSumRoh = HUGE(minSumRoh)
     currentSumRoh = 0.0D0
           
     InfiniteLoop: do
        
        currentO = Olist(assignedOcount)

        if (assignedOcount.eq.1) then ! The Eigen O in this state, which has 3 associated Hs

! Check if we arrived here as a result of increasing or decreasing assignOcount
              
           if (goingUp) then     ! Initialisation
           usedH(proximalAtoms(currentO,Hcounters(assignedOcount,1))) = .TRUE.
           usedH(proximalAtoms(currentO,Hcounters(assignedOcount,2))) = .TRUE.
           usedH(proximalAtoms(currentO,Hcounters(assignedOcount,3))) = .TRUE.

           currentSumRoh = interAtomicR(currentO,proximalAtoms(currentO,Hcounters(assignedOcount,1)))+&
                                    interAtomicR(currentO,proximalAtoms(currentO,Hcounters(assignedOcount,2)))+& 
                                 interAtomicR(currentO,proximalAtoms(currentO,Hcounters(assignedOcount,3)))

           assignedOcount = assignedOcount + 1

           else          ! Increment the counters if possible, exit if all iterations are done

           if (Hcounters(assignedOcount,3).eq.neighbourCounts(currentO)) then
              
              if (Hcounters(assignedOcount,2).eq.(neighbourCounts(currentO)-1)) then
                 
                 if (Hcounters(assignedOcount,1).eq.(neighbourCounts(currentO)-2)) then
                 exit InfiniteLoop
                 else
                 usedH(proximalAtoms(currentO,Hcounters(assignedOcount,1))) = .FALSE.
                 usedH(proximalAtoms(currentO,Hcounters(assignedOcount,2))) = .FALSE.
                 usedH(proximalAtoms(currentO,Hcounters(assignedOcount,3))) = .FALSE.
                 Hcounters(assignedOcount,1) = Hcounters(assignedOcount,1) + 1
                 Hcounters(assignedOcount,2) = Hcounters(assignedOcount,1) + 1
                 Hcounters(assignedOcount,3) = Hcounters(assignedOcount,2) + 1
                 usedH(proximalAtoms(currentO,Hcounters(assignedOcount,1))) = .TRUE.
                 usedH(proximalAtoms(currentO,Hcounters(assignedOcount,2))) = .TRUE.
                 usedH(proximalAtoms(currentO,Hcounters(assignedOcount,3))) = .TRUE.

                 currentSumRoh = &
                               interAtomicR(currentO,proximalAtoms(currentO,Hcounters(assignedOcount,1))) + &
                               interAtomicR(currentO,proximalAtoms(currentO,Hcounters(assignedOcount,2))) + &
                               interAtomicR(currentO,proximalAtoms(currentO,Hcounters(assignedOcount,3)))
       
                 if (currentSumRoh.ge.minSumRoh) then
                    goingUp = .FALSE.
                 else
                    goingUp = .TRUE.
                 endif
                 endif

              else 
                 usedH(proximalAtoms(currentO,Hcounters(assignedOcount,2))) = .FALSE.
                 usedH(proximalAtoms(currentO,Hcounters(assignedOcount,3))) = .FALSE.
                 Hcounters(assignedOcount,2) = Hcounters(assignedOcount,2) + 1
                 Hcounters(assignedOcount,3) = Hcounters(assignedOcount,2) + 1
                 usedH(proximalAtoms(currentO,Hcounters(assignedOcount,2))) = .TRUE.
                 usedH(proximalAtoms(currentO,Hcounters(assignedOcount,3))) = .TRUE.

                 currentSumRoh = &
                            interAtomicR(currentO,proximalAtoms(currentO,Hcounters(assignedOcount,1))) + &
                            interAtomicR(currentO,proximalAtoms(currentO,Hcounters(assignedOcount,2))) + &
                            interAtomicR(currentO,proximalAtoms(currentO,Hcounters(assignedOcount,3)))

                 if (currentSumRoh.ge.minSumRoh) then
                 goingUp = .FALSE.
                 else
                 goingUp = .TRUE.
                 endif
                 
              endif
           else 
              usedH(proximalAtoms(currentO,Hcounters(assignedOcount,3))) = .FALSE.
              Hcounters(assignedOcount,3) = Hcounters(assignedOcount,3) + 1
              usedH(proximalAtoms(currentO,Hcounters(assignedOcount,3))) = .TRUE.       

              currentSumRoh = &
                         interAtomicR(currentO,proximalAtoms(currentO,Hcounters(assignedOcount,1))) + &
                         interAtomicR(currentO,proximalAtoms(currentO,Hcounters(assignedOcount,2))) + &
                         interAtomicR(currentO,proximalAtoms(currentO,Hcounters(assignedOcount,3)))

              if (currentSumRoh.ge.minSumRoh) then
                 goingUp = .FALSE.
              else
                 goingUp = .TRUE.
              endif
              
           endif

           if (goingUp) assignedOcount = assignedOcount + 1

           endif

        else               ! A water O, has 2 associated Hs

! Check whether we have arrived here as a result of increasing or decreasing assignedOcount

           if (goingUp) then ! Initialisation

!           print *, 'Adding a new water O:', currentO
!                 print *, 'Assignment to date:'
!                 do i=1,assignedOcount-1
!                 print *, 'O:', (3*Olist(i))
!                 if (i.eq.1) then
!                    print *, 'H1:', Hposition(neighbouringHs(Olist(i),Hcounters(i,1))),
!    &                               'H2:', Hposition(neighbouringHs(Olist(i),Hcounters(i,2))),
!    &                             'H3:', Hposition(neighbouringHs(Olist(i),Hcounters(i,3)))
!                 else
!                    print *, 'H1:', Hposition(neighbouringHs(Olist(i),Hcounters(i,1))),
!    &                             'H2:', Hposition(neighbouringHs(Olist(i),Hcounters(i,2)))
!                 endif
!                 enddo

!                 print *, 'roh_sq to date:', currentSumRohSq

           foundAnH = .FALSE.

           FirstHLoop: do i=1, neighbourCounts(currentO)-1

              if (.not.usedH(proximalAtoms(currentO,i))) then

! If adding this OH distance brings the total to more than the current minimum then any subsequent H will also do
! so since they are arranged in ascending distance from the O

                 if ((currentSumRoh+interAtomicR(currentO,proximalAtoms(currentO,i))).ge.minSumRoh) then
                 exit FirstHLoop
                 else

                 do j = i+1, neighbourCounts(currentO)
                    if (.not.usedH(proximalAtoms(currentO,j))) then

                    if ((currentSumRoh + interAtomicR(currentO,proximalAtoms(currentO,i)) + &
                                     interAtomicR(currentO,proximalAtoms(currentO,j))).ge.minSumRoh) then
                       exit FirstHLoop
                    else
                       Hcounters(assignedOcount,1) = i
                       usedH(proximalAtoms(currentO,i)) = .TRUE.
                       Hcounters(assignedOcount,2) = j
                       usedH(proximalAtoms(currentO,j)) = .TRUE.
                       foundAnH = .TRUE.

                       currentSumRoh = currentSumRoh + &
                                        interAtomicR(currentO,proximalAtoms(currentO,i)) + &
                                        interAtomicR(currentO,proximalAtoms(currentO,j))
                       exit FirstHLoop
                    endif
                    endif
                 enddo
                 endif
              endif
           enddo FirstHLoop

           if (foundAnH) then

!              print *, 'Successfully allocated new Hs'
              if (assignedOcount.eq.num_eig) then

                 pivotVBstate(1) = proximalAtoms(pivotOxygen, Hcounters(1,1))
                 pivotVBstate(2) = proximalAtoms(pivotOxygen, Hcounters(1,2))
                 pivotVBstate(3) = pivotOxygen
                 pivotVBstate(4) = proximalAtoms(pivotOxygen, Hcounters(1,3))

                 i = 5

                 do j = 2, num_eig
                 pivotVBstate(i) = proximalAtoms(Olist(j),Hcounters(j,1))
                 pivotVBstate(i+1) = Olist(j)
                 pivotVBstate(i+2) = proximalAtoms(Olist(j),Hcounters(j,2))
                 i = i + 3
                 enddo

                 minSumRoh = currentSumRoh

!                 print *, 'New assignment:', minSumRoh
!                    do i = 1,num_eig
!                    print *, 'O:', (3*Olist(i))
!                    if (i.eq.1) then
!                       print *, 'H1:', Hposition(neighbouringHs(Olist(i),Hcounters(i,1))),
!    &                             'H2:', Hposition(neighbouringHs(Olist(i),Hcounters(i,2))),
!    &                             'H3:', Hposition(neighbouringHs(Olist(i),Hcounters(i,3)))
!                    else
!                       print *, 'H1:', Hposition(neighbouringHs(Olist(i),Hcounters(i,1))),
!    &                                'H2:', Hposition(neighbouringHs(Olist(i),Hcounters(i,2)))
!                    endif
!                    enddo

                 goingUp = .FALSE.
              else
                 assignedOcount = assignedOcount + 1
              endif
           else
!              print *, 'Could not allocate new Hs, dropping a level'
              assignedOcount = assignedOcount - 1
              goingUp = .FALSE.
           endif   
           else 

!           print *, 'Attempting to increment O:', currentO
                 
!                 print *, 'Assignment at this point:'
                 
!                 do i=1,assignedOcount
!                 print *, 'O:', (3*Olist(i))
!                 if (i.eq.1) then
!                    print *, 'H1:', Hposition(neighbouringHs(Olist(i),Hcounters(i,1))),
!    &                               'H2:', Hposition(neighbouringHs(Olist(i),Hcounters(i,2))),
!    &                             'H3:', Hposition(neighbouringHs(Olist(i),Hcounters(i,3)))
!                 else
!                    print *, 'H1:', Hposition(neighbouringHs(Olist(i),Hcounters(i,1))),
!    &                             'H2:', Hposition(neighbouringHs(Olist(i),Hcounters(i,2)))
!                 endif
!                 enddo

!                 print *, 'Current sum roh_sq:', currentSumRohSq

! Increment the current counter if possible, else drop down another level

           currentSumRoh = currentSumRoh - &
                      interAtomicR(currentO,proximalAtoms(currentO,Hcounters(assignedOcount,2)))

           usedH(proximalAtoms(currentO,Hcounters(assignedOcount,2))) = .FALSE.

           if (Hcounters(assignedOcount,2).eq.neighbourCounts(currentO)) then                 
              if (Hcounters(assignedOcount,1).eq.(neighbourCounts(currentO)-1)) then                 

                 currentSumRoh = currentSumRoh - & 
                            interAtomicR(currentO,proximalAtoms(currentO,Hcounters(assignedOcount,1)))

                 usedH(proximalAtoms(currentO,Hcounters(assignedOcount,1))) = .FALSE.
!                    Hcounters(assignedOcount,1) = 1
!                    Hcounters(assignedOcount,2) = 2
                 assignedOcount = assignedOcount - 1
              else
                 usedH(proximalAtoms(currentO,Hcounters(assignedOcount,1))) = .FALSE.
                 currentSumRoh = currentSumRoh - &
                            interAtomicR(currentO,proximalAtoms(currentO,Hcounters(assignedOcount,1)))

! Find a suitable H atom for the first H atom

                 foundAnH = .FALSE.

                 FirstHLoop2: do i = Hcounters(assignedOcount,1)+1, neighbourCounts(currentO)-1
                    if (.not.usedH(proximalAtoms(currentO,i))) then

                    if ((currentSumRoh + interAtomicR(currentO,proximalAtoms(currentO,i))).ge.minSumRoh) then
                    exit FirstHLoop2
                    else
                    do j = i+1, neighbourCounts(currentO) 
                       if (.not.usedH(proximalAtoms(currentO,j))) then
                         
                          if ((currentSumRoh + &
                                           interAtomicR(currentO,proximalAtoms(currentO,i)) + &
                                           interAtomicR(currentO,proximalAtoms(currentO,j))).ge.minSumRoh) then
                          exit FirstHLoop2
                          else
                          Hcounters(assignedOcount,1) = i
                          usedH(proximalAtoms(currentO,i)) = .TRUE.
                          Hcounters(assignedOcount,2) = j
                          usedH(proximalAtoms(currentO,j)) = .TRUE.
                          foundAnH = .TRUE.

                          currentSumRoh = currentSumRoh + &
                                              interAtomicR(currentO,proximalAtoms(currentO,i)) + &
                                              interAtomicR(currentO,proximalAtoms(currentO,j))

                          exit FirstHLoop2
                          endif
                       endif
                    enddo
                    endif
                 endif
                 enddo FirstHLoop2

                 if (foundAnH) then
                 if (assignedOcount.eq.num_eig) then

                    pivotVBstate(1) = proximalAtoms(pivotOxygen, Hcounters(1,1))
                    pivotVBstate(2) = proximalAtoms(pivotOxygen, Hcounters(1,2))
                    pivotVBstate(3) = pivotOxygen
                    pivotVBstate(4) = proximalAtoms(pivotOxygen, Hcounters(1,3))

                    i = 5
                       
                    do j = 2, num_eig
                    pivotVBstate(i) = proximalAtoms(Olist(j),Hcounters(j,1))
                    pivotVBstate(i+1) = Olist(j)
                    pivotVBstate(i+2) = proximalAtoms(Olist(j),Hcounters(j,2))
                    i = i + 3
                    enddo

                    minSumRoh = currentSumRoh

!                    print *, 'New assignment:', minSumRoh
!                       do i = 1,num_eig
!                          print *, 'O:', (3*Olist(i))
!                          if (i.eq.1) then
!                          print *, 'H1:', Hposition(neighbouringHs(Olist(i),Hcounters(i,1))),
!    &                                   'H2:', Hposition(neighbouringHs(Olist(i),Hcounters(i,2))),
!    &                                   'H3:', Hposition(neighbouringHs(Olist(i),Hcounters(i,3)))
!                          else
!                          print *, 'H1:', Hposition(neighbouringHs(Olist(i),Hcounters(i,1))),
!    &                                      'H2:', Hposition(neighbouringHs(Olist(i),Hcounters(i,2)))
!                          endif
!                       enddo
                 else
                    assignedOcount = assignedOcount + 1
                    goingUp = .TRUE.
                 endif
                 else 

! Drop down a level
                 assignedOcount = assignedOcount - 1
                 endif
              endif
           else 

! Check for an unused H atom
                 
              foundAnH = .FALSE.
              
              do i = Hcounters(assignedOcount,2)+1,neighbourCounts(currentO)
                 if (.not.usedH(proximalAtoms(currentO,i))) then

                 if ((currentSumRoh + interAtomicR(currentO,proximalAtoms(currentO,i))).ge.minSumRoh) then
                    exit
                 else
                    Hcounters(assignedOcount,2) = i
                    usedH(proximalAtoms(currentO,i)) = .TRUE.
                    foundAnH = .TRUE.
                    currentSumRoh = currentSumRoh + interAtomicR(currentO,proximalAtoms(currentO,i))
                    exit
                 endif
                 endif
              enddo

              if (foundAnH) then

                 if (assignedOcount.eq.num_eig) then

                 pivotVBstate(1) = proximalAtoms(pivotOxygen, Hcounters(1,1))
                 pivotVBstate(2) = proximalAtoms(pivotOxygen, Hcounters(1,2))
                 pivotVBstate(3) = pivotOxygen
                 pivotVBstate(4) = proximalAtoms(pivotOxygen, Hcounters(1,3))

                 i = 5

                 do j = 2, num_eig
                    pivotVBstate(i) = proximalAtoms(Olist(j),Hcounters(j,1))
                    pivotVBstate(i+1) = Olist(j)
                    pivotVBstate(i+2) = proximalAtoms(Olist(j),Hcounters(j,2))
                    i = i + 3
                 enddo

                 minSumRoh = currentSumRoh

!                 print *, 'New assignment:', minSumRoh
!                    do i = 1,num_eig
!                       print *, 'O:', (3*Olist(i))

!                       if (i.eq.1) then
!                          print *, 'H1:', Hposition(neighbouringHs(Olist(i),Hcounters(i,1))),
!    &                                         'H2:', Hposition(neighbouringHs(Olist(i),Hcounters(i,2))),
!    &                                         'H3:', Hposition(neighbouringHs(Olist(i),Hcounters(i,3)))
!                       else
!                          print *, 'H1:', Hposition(neighbouringHs(Olist(i),Hcounters(i,1))),
!    &                                         'H2:', Hposition(neighbouringHs(Olist(i),Hcounters(i,2)))
!                       endif
!                    enddo
                 else
                 assignedOcount = assignedOcount + 1
                 goingUp = .TRUE.
                 endif

! Could not find an acceptable second H, try running through the first H selection in tandem

              else     

                 usedH(proximalAtoms(currentO,Hcounters(assignedOcount,1))) = .FALSE.
                 currentSumRoh = currentSumRoh - &
                            interAtomicR(currentO,proximalAtoms(currentO,Hcounters(assignedOcount,1)))

                 foundAnH = .FALSE.

                 FirstHLoop3: do i = Hcounters(assignedOcount,1)+1, neighbourCounts(currentO)-1
                    if (.not.usedH(proximalAtoms(currentO,i))) then

                    if ((currentSumRoh + interAtomicR(currentO,proximalAtoms(currentO,i))).ge.minSumRoh) then
                    exit FirstHLoop3
                    else

                    do j = i+1, neighbourCounts(currentO) 
                       if (.not.usedH(proximalAtoms(currentO,j))) then

                          if ((currentSumRoh + &
                                           interAtomicR(currentO,proximalAtoms(currentO,i)) + &
                                           interAtomicR(currentO,proximalAtoms(currentO,j))).ge.minSumRoh) then
                          exit FirstHLoop3
                          else

                          Hcounters(assignedOcount,1) = i
                          usedH(proximalAtoms(currentO,i)) = .TRUE.
                          Hcounters(assignedOcount,2) = j
                          usedH(proximalAtoms(currentO,j)) = .TRUE.

                          currentSumRoh = currentSumRoh + &
                                              interAtomicR(currentO,proximalAtoms(currentO,i)) + &
                                              interAtomicR(currentO,proximalAtoms(currentO,j))
                          foundAnH = .TRUE.
                          exit FirstHLoop3
                          endif
                       endif
                    enddo
                    endif
                 endif
                 enddo FirstHLoop3

                 if (foundAnH) then
                 if (assignedOcount.eq.num_eig) then

                    pivotVBstate(1) = proximalAtoms(pivotOxygen, Hcounters(1,1))
                    pivotVBstate(2) = proximalAtoms(pivotOxygen, Hcounters(1,2))
                    pivotVBstate(3) = pivotOxygen
                    pivotVBstate(4) = proximalAtoms(pivotOxygen, Hcounters(1,3))

                    i = 5

                    do j = 2, num_eig
                    pivotVBstate(i) = proximalAtoms(Olist(j),Hcounters(j,1))
                    pivotVBstate(i+1) = Olist(j)
                    pivotVBstate(i+2) = proximalAtoms(Olist(j),Hcounters(j,2))
                    i = i + 3
                    enddo

                    minSumRoh = currentSumRoh

!                    print *, 'New assignment:', minSumRoh
!                       do i = 1,num_eig
!                          print *, 'O:', (3*Olist(i))
!                          if (i.eq.1) then
!                          print *, 'H1:', Hposition(neighbouringHs(Olist(i),Hcounters(i,1))),
!    &                                    'H2:', Hposition(neighbouringHs(Olist(i),Hcounters(i,2))),
!    &                                      'H3:', Hposition(neighbouringHs(Olist(i),Hcounters(i,3)))
!                          else
!                          print *, 'H1:', Hposition(neighbouringHs(Olist(i),Hcounters(i,1))),
!    &                                      'H2:', Hposition(neighbouringHs(Olist(i),Hcounters(i,2)))
!                          endif
!                       enddo
                 else
                    assignedOcount = assignedOcount + 1                 
                    goingUp = .TRUE.
                 endif
                 else 

! Drop down a level

                 assignedOcount = assignedOcount - 1
                 endif                 
              endif
           endif
           endif
        endif
     enddo InfiniteLoop

! Check also that it was possible to assign a complete VB state

     do i=1,natoms
        if (pivotVBstate(i).eq.-1) then
           print *, 'Unable to assign a complete VB state for pivot oxygen', pivotOxygen
           print *, 'Check geometry or cutoff'
           success = .FALSE.
           exit
        endif
     enddo

     return

     end subroutine

! ###############################################################################################

! Calculate interatomic distances - these will be used by the energy/force calculations

     subroutine calcInterAtomicDistances(natoms,xcoords,ycoords,zcoords,separations)
     implicit none

! Subroutine arguments

     integer, intent(IN) :: natoms
     DOUBLE PRECISION, INTENT (IN) :: XCOORDS(NATOMS), YCOORDS(NATOMS), ZCOORDS(NATOMS)
     DOUBLE PRECISION, DIMENSION(NATOMS,NATOMS), INTENT (OUT) :: SEPARATIONS

! Local variables

     integer atom1, atom2

! Initialise

     separations = 0.0d0

! Calculate the upper triangle but explicitly fill in the whole matrix to avoid in future having to work out
! which atom number is lower

     do atom1 = 1, (natoms-1)
        do atom2 = (atom1+1), natoms
           separations(atom1,atom2) = DSQRT((xcoords(atom1)-xcoords(atom2))**2 + &
                                               (ycoords(atom1)-ycoords(atom2))**2 + &
                                               (zcoords(atom1)-zcoords(atom2))**2)
           separations(atom2,atom1) = separations(atom1,atom2)
        enddo
     enddo

     return

     end subroutine

! ##################################################################################################     

