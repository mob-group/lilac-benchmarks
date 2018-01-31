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
subroutine drvmsevb(ncoord, Q, dQ, energy, computeGrad, computeHess)
	use commons
	use msevb_common
	use msevb_interfaces
        use MODHESS

	implicit none

!###################################################################

	integer ncoord
	DOUBLE PRECISION :: Q(NCOORD), DQ(NCOORD)
	DOUBLE PRECISION ENERGY
	logical computeGrad, computeHess

        DOUBLE PRECISION :: DELTA_COORD

! Debugging variables

	integer :: i,j
        DOUBLE PRECISION RMS_DHESS

!	DOUBLE PRECISION RMS_DQ, NUMERICAL_RMS_DQ
!	DOUBLE PRECISION DELTA_Q, UPPER_ENERGY, LOWER_ENERGY, UPPER_DIAG, LOWER_DIAG
!	DOUBLE PRECISION NUMERICAL_DQ(NCOORD)
!	DOUBLE PRECISION DQ_DIAGONAL(3*MXATMS), NUMERICAL_DQ_DIAGONAL(3*MXATMS)

!	DOUBLE PRECISION UPPER_Q_LJ(MXEIG), UPPER_Q_COULOMB(MXEIG)
!	DOUBLE PRECISION UPPER_Q_REP(MXEIG)

!	DOUBLE PRECISION UPPER_Q_H3O_INTRA(MXEIG)
!	DOUBLE PRECISION UPPER_Q_H2O_INTRA(MXEIG)
!	DOUBLE PRECISION UPPER_Q_H3O_H2O(MXEIG)
!	DOUBLE PRECISION UPPER_Q_H2O_H2O(MXEIG)

!	DOUBLE PRECISION LOWER_Q_H3O_INTRA(MXEIG)
!	DOUBLE PRECISION LOWER_Q_H2O_INTRA(MXEIG)
!	DOUBLE PRECISION LOWER_Q_H3O_H2O(MXEIG)
!	DOUBLE PRECISION LOWER_Q_H2O_H2O(MXEIG)

!	DOUBLE PRECISION UPPER_Q_OFF_DIAG(MXEIG,MXEIG)

!	DOUBLE PRECISION NUM_FORCE

!	integer requiredCoordcount
!	integer requiredCoord(3*MXATMS)
!	integer newCoord

! Calculate the energy

	call msevb(ncoord, Q, .TRUE., energy, .FALSE.)

        if (printCoefficients) then
           print *, '#########################'
           print *, 'MSEVB coefficients:'
           do i = 1, reduced_num_eig
              print *, i, grstwfu(i)
           enddo
           print *, '#########################'
        end if

!	print *, 'Energy:', energy

!	print *, 'Points:'
!	j = 1
!	do i = 1, natoms

!	   if (MOD(i,3).eq.0) then
!	      print *, 'O', Q(j), Q(j+1), Q(j+2)
!	   else
!	      print *, 'H', Q(j), Q(j+1), Q(j+2)
!	   endif

!	   j = j + 3
!	enddo

! 	print *, 'VB states:'

! 	do i=1,reduced_num_eig
! 	   print *, 'State:', i

! 	   do j=1, natoms
! 	      print *, j, atmpl(i,j)
!  	   enddo
!  	enddo

! 	print *, 'Ground state coefficients and diagonal energies:'
! 	do i=1, reduced_num_eig
! 	   print *, i, grstwfu(i), ham(i,i)
! 	enddo

! 	print *, 'Hamiltonian:'

! 	do i=1,reduced_num_eig
! 	   do j=i, reduced_num_eig
! 	      print *, i , j, ham(i,j)
! 	   enddo
! 	enddo

!	print *, 'Zundel species:'

! 	do i=1, reduced_num_eig-1
! 	   do j=i+1, reduced_num_eig

!	      if (statesInteract(i,j)) then
! 		 print *, 'State 1:', i, 'State 2:', j

!		 do k=1,7
!		    print *, k, zundel_species(i,j,k)
!		 enddo
!	      endif
! 	   enddo
!  	enddo

!	do i=1,num_eig
!	   print *, '*** VB state:', i, '***'
!	   print *, 'H3O intra energy:', h3o_intra_energy(i)
!	   print *, 'H2O intra energy:', h2o_intra_energy(i)
!	   print *, 'H3O-H2O inter energy:', h3o_h2O_energy(i)
!	   print *, 'H2O-H2O inter energy:', h2o_h2o_energy(i)
!	   print *, ' '
!	enddo

!	print *, 'Exchange interactions:'
!	do i=1, num_eig-1
!	   do j=i+1, num_eig
!	      print *, i, j, ham(i,j)
!	   enddo
!	enddo

!	stop

	if (computeGrad) call fmsevb(ncoord, Q, dQ, .FALSE.)

        if (computeHess) then
           delta_coord = 1.0d-6
           call MSEVB_hess (ncoord, Q, delta_coord)

!           rms_dHess=0.0d0

!           do i=1, ncoord-1
!              do j = i+1, ncoord
!                 rms_dHess = rms_dHess + (HESS(i,j)-HESS(j,i))**2
!              enddo
!           enddo

!           rms_dHess = dsqrt(2.0d0*rms_dHess/(ncoord**2 - ncoord))
           
!           print *, 'RMS dHess:', rms_dHess

!           stop
        endif

! Deallocate memory

	DEALLOCATE(vijexch, grstwfu)

!	endif

!	print *, 'Forces:'
!	j=1
!	do i=1,ncoord,3
!	   print *, j, dQ(i), dQ(i+1), dQ(i+2)
!	   j = j + 1
!	enddo

! Calculate the numerical gradients

!	delta_Q = 1D-7
!	do i=1, ncoord

!	   print *, '*** Coordinate:', i, '***'

!	   Q(i) = Q(i) + delta_Q
! 	   call msevb(ncoord, Q, .FALSE., upper_energy)

!	   do j=1,reduced_num_eig
!	      upper_Q_h3o_intra(j) = h3o_intra_energy(j)
!	      upper_Q_h2o_intra(j) = h2o_intra_energy(j)
!	      upper_Q_h3o_h2o(j) = h3o_h2o_energy(j)
!	      upper_Q_h2o_h2o(j) = h2o_h2o_energy(j)
!	   enddo

!  	   Q(i) = Q(i) - 2*delta_Q
! 	   call msevb(ncoord, Q, .FALSE., lower_energy)
!  	   numerical_dQ(i) = (upper_energy-lower_energy)/(2*delta_Q)

!	   do j=1,reduced_num_eig
!	      print *, 'VB state:', j

!	      num_force = (upper_Q_h3o_intra(j)-h3o_intra_energy(j))/(2*delta_Q)
!	      print *, 'H3O intra forces:', h3o_intra_forces(j,i), num_force, (h3o_intra_forces(j,i)-num_force)

!	      num_force = (upper_Q_h2o_intra(j)-h2o_intra_energy(j))/(2*delta_Q)
!	      print *, 'H2O intra forces:', h2o_intra_forces(j,i), num_force, (h2o_intra_forces(j,i)-num_force)

!	      num_force = (upper_Q_h3o_h2o(j)-h3o_h2o_energy(j))/(2*delta_Q)
!	      print *, 'H3O-H2O forces:', h3o_h2o_forces(j,i), num_force, (h3o_h2o_forces(j,i)-num_force)

!	      num_force = (upper_Q_h2o_h2o(j)-h2o_h2o_energy(j))/(2*delta_Q)
!	      print *, 'H2O-H2O forces:', h2o_h2o_forces(j,i), num_force, (h2o_h2o_forces(j,i)-num_force)

!	   enddo

!  	   Q(i) = Q(i) + delta_Q
!  	enddo

! Calculate the numerical gradients

!	delta_Q = 1D-7

!	requiredCoordcount = 6
!	requiredCoord(1)=52
!	requiredCoord(2)=53
!	requiredCoord(3)=54
!	requiredCoord(4)=61
!	requiredCoord(5)=62
!	requiredCoord(6)=63

!	do i=1, requiredCoordCount

!	   newCoord = 3*new_order(int((requiredCoord(i)-1)/3)+1) - 2 + MOD(requiredCoord(i)-1, 3)

!	   print *, ' '
!	   print *, '*** Coordinate:', requiredCoord(i), newCoord, '***'
!	   print *, ' '

!	   print *, '*** Calculating upper energy bound ***'
!	   Q(i) = Q(i) + delta_Q
!	   call msevb_ind(ncoord, Q, upper_energy, wf_buffer)

!	   do j=1, num_eig
!	      upper_Q_LJ(j)=LJ_values(j)
!	      upper_Q_coulomb(j)=coulomb_values(j)
!	      upper_Q_rep(j)=rep_values(j)
! 	   upper_Q_h3o_intra(j)=h3o_intra_energy(j)
!	      upper_Q_h2o_intra(j)=h2o_intra_energy(j)
!	      upper_Q_h3o_h2o(j)=h3o_h2o_energy(j)
!	      upper_Q_h2o_h2o(j)=h2o_h2o_energy(j)
! 	   enddo

! 	   do j=1,num_eig-1
! 	      do k=j+1, num_eig
! 		 upper_Q_off_diag(j,k)=off_diag_energies(j,k)
! 	      enddo
! 	   enddo

! 	   print *, '*** Calculating lower energy bound ***'
!	   Q(i) = Q(i) - 2*delta_Q
!	   call msevb_ind(ncoord, Q, lower_energy, wf_buffer)

!	   numerical_dQ(i) = (upper_energy-lower_energy)/(2*delta_Q)

!	   Q(i) = Q(i) + delta_Q

! 	   print *, 'Diagonal components:'

! 	   do j=1, num_eig
!	      print *, ' '
! 	      print *, 'VB state:', j

! 	      num_force=(upper_Q_h3o_intra(j)-h3o_intra_energy(j))/(2*delta_Q)
! 	      print *, 'H3O+ intramolecular contribution:'
! 	      print *, h3o_intra_forces(j,i), num_force, (h3o_intra_forces(j,i)-num_force)
!	      print *, h3o_intra_forces(j,newCoord)


! 	      num_force=(upper_Q_h2o_intra(j)-h2o_intra_energy(j))/(2*delta_Q)
! 	      print *, 'H2O intramolecular contribution:'
! 	      print *, h2o_intra_forces(j,i), num_force, (h2o_intra_forces(j,i)-num_force)
! 	      print *, h2o_intra_forces(j,newCoord)

! 	      num_force=(upper_Q_h3o_h2o(j)-h3o_h2o_energy(j))/(2*delta_Q)
! 	      print *, 'H3O+-H2O intermolecular contribution:'
! 	      print *, h3o_h2O_forces(j,i), num_force, (h3o_h2O_forces(j,i)-num_force)
! 	      print *, h3o_h2O_forces(j,newCoord)

! 	      num_force=(upper_Q_h2o_h2o(j)-h2o_h2o_energy(j))/(2*delta_Q)
! 	      print *, 'H2O-H2O intermolecular contribution:'
! 	      print *, h2o_h2o_forces(j,i), num_force, (h2o_h2o_forces(j,i)-num_force)
! 	      print *, h2o_h2o_forces(j,newCoord)

! 	      print *, 'H2O-H3O interaction:'
! 	      print *, 'Component, Analytic, Numerical, Difference'

! 	      num_force=(upper_Q_LJ(j)-LJ_values(j))/(2*delta_Q)
! 	      print *, 'LJ force:', LJ_forces(j,i), num_force, (LJ_forces(j,i)-num_force)
! 	      print *, 'LJ force:', LJ_forces(j,newCoord)

! 	      num_force=(upper_Q_coulomb(j)-coulomb_values(j))/(2*delta_Q)
! 	      print *, 'Coulomb force:', coulomb_forces(j,i), 
!     &                 num_force, (coulomb_forces(j,i)-num_force)
!	      print *, 'Coulomb force:', coulomb_forces(j,newCoord)

! 	      num_force=(upper_Q_rep(j)-rep_values(j))/(2*delta_Q)
! 	      print *, 'Repulsive force:', rep_forces(j,i), num_force, 
!     &                 (rep_forces(j,i)-num_force) 
! 	      print *, 'Repulsive force:', rep_forces(j,newCoord)

! 	   enddo

! 	   print *, ' '
! 	   print *, 'Off-diagonal components:'

! 	   do j=1, num_eig-1
! 	      do k=j+1,num_eig
! 		 num_force=(upper_Q_off_diag(j,k)-off_diag_energies(j,k))/(2*delta_Q)
! 		 print *, 'VB state 1, state 2, analytic force, numeric, difference'
! 		 print *, j, k, off_diag_forces(j,k,i),num_force,(off_diag_forces(j,k,i)-num_force)
! 		 print *, 'VB state 1, state 2, analytic force'
!		 print *, j, k, off_diag_forces(j,k,newCoord)

! 		 print *, 'Numerical energy bounds:'
! 		 print *, 'Upper:', upper_Q_off_diag(j,k), 'Lower:',off_diag_energies(j,k)
! 	      enddo
! 	   enddo

! 	   print *, ' '

!	enddo   

!  	print *, 'Derivatives: '
!  	print *, 'Coordinate Number'
!  	print *, 'Analytical, Numerical, Difference'

!   	rms_dQ = 0.0D0
!   	numerical_rms_dQ = 0.0D0

!   	do i=1,ncoord
!   	   rms_dQ = rms_dQ + (dQ(i)*dQ(i))
!   	   numerical_rms_dQ = numerical_rms_dQ + (numerical_dQ(i)*numerical_dQ(i))
!   	   print *, i, dQ(i), numerical_dQ(i), (dQ(i)-numerical_dQ(i))
!   	enddo

!   	rms_dQ = SQRT(rms_dQ/ncoord)
!   	numerical_rms_dQ = SQRT(numerical_rms_dQ/ncoord)

!   	print *, 'RMS derivatives: ', rms_dQ, ' ', numerical_rms_dQ, ' ', (rms_dQ-numerical_rms_dQ) 
!   	print *, 'Delta_Q: ', delta_Q

!   	stop

! End debugging

	return
	end

! ##################################################################################################
	
! Numerical Hessian routine for MSEVB potential
! Do not call this directly, only via drvmsevb

subroutine MSEVB_hess (ncoord, coords, delta_coord)
  use msevb_interfaces
  use MODHESS, hessian => HESS

  implicit none

! Subroutine arguments

  integer, intent(IN) :: ncoord
  DOUBLE PRECISION, INTENT(IN) :: COORDS(NCOORD)
  DOUBLE PRECISION, INTENT(IN) :: DELTA_COORD

! Local variables

  integer :: currentCoord, currentChange
  DOUBLE PRECISION :: LOCALCOORDS(NCOORD)
  DOUBLE PRECISION :: PLUSSIDEGRAD(NCOORD), MINUSSIDEGRAD(NCOORD)
  DOUBLE PRECISION :: ENERGY

! Try to avoid problems with arithmetic by referring always to the original copy of the coordinates

  localCoords = coords

  do currentCoord = 1, ncoord

     localCoords(currentCoord) = coords(currentCoord) + delta_coord

! We do not want to recompute the VB states because we assume that they don`t change

     call msevb(ncoord, localCoords, .FALSE., energy, .FALSE.)
     call fmsevb(ncoord, localCoords, plusSideGrad, .FALSE.)

     localCoords(currentCoord) = coords(currentCoord) - delta_coord     

     call msevb(ncoord, localCoords, .FALSE., energy, .FALSE.)
     call fmsevb(ncoord, localCoords, minusSideGrad, .FALSE.)     

     do currentChange = currentCoord, ncoord
        hessian(currentChange,currentCoord) = &
        (plusSideGrad(currentChange)-minusSideGrad(currentChange))/(2.0d0*delta_coord)

        if (currentChange.ne.currentCoord) then
           hessian(currentCoord,currentChange) = hessian(currentChange,currentCoord)
        endif
     enddo
     
     localCoords(currentCoord) = coords(currentCoord)      
  enddo

end subroutine MSEVB_hess	

! ##################################################################################################











