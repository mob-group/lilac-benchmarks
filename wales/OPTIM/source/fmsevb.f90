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
!###################################################################################
!############################ FMSEVB SUBROUTINE ####################################
!###################################################################################

	subroutine fmsevb(ncoord, Q, dQ, assignVBstates)

	use commons
	use msevb_common
	use msevb_interfaces

	implicit none

!#####################################################################
! 24-09-00  R. A. Christie
!
! Parameter set used is from:
!    U. Schmitt and G. A. Voth, J. Chem. Phys. 111, 1999, 9361
! and not the earlier set, nor the later one, or even the buggy one.
!
! Call the various routines necessary to build up the Hellman-Feynman
!
! *NB* This routine used to be called hellfeyn, but the name was 
!      altered on 3rd Aug 2002 when introduced into OPTIM 2.3
!
!#####################################################################  

	integer ncoord
	DOUBLE PRECISION DQ(NCOORD), Q(NCOORD)
	logical assignVBstates

	DOUBLE PRECISION XF(NATOMS),YF(NATOMS),ZF(NATOMS)
	integer atompos,atompos_temp,proton,atom_num
	integer i,jj,j,ii
	DOUBLE PRECISION FXH(REDUCED_NUM_EIG,REDUCED_NUM_EIG),FYH(REDUCED_NUM_EIG,REDUCED_NUM_EIG), &
               fzH(reduced_num_eig,reduced_num_eig),dx,dy, &
      	       dz,R1,R2,minimum_x(natoms),scale,minimum_y(natoms),minimum_z(natoms), &
               probable_state(reduced_num_eig),tempScale
	integer m,eigenstate,pivotOxygen,jk,posn_in_state1,position,position_temp
	DOUBLE PRECISION FXATOM1,FXATOM3,FXATOM4(NATOMS),OFFDIAGX,OFFDIAGY, &
      		offdiagz,Roo,Roh,Fyatom1,Fzatom1, &
      		Fyatom3,Fzatom4(natoms),Fyatom4(natoms),Fzatom3, &
      		force_x,force_y,force_z, &
      		dRoo(3),dRoh(3),startEnergy,maxforcemag,keech
	character*5 atomtype
	DOUBLE PRECISION DX1,DX2,DX3,DY1,DY2,DY3,DZ1,DZ2,DZ3,R12_1,R12_2,R12_3,ASYMM
	integer hmin,positioni,positionj,jl,atomorder(natoms)
	DOUBLE PRECISION ENEW
	DOUBLE PRECISION R,RCUT
	integer currentAtom,state1,positionInState1,state2,positionInState2
	DOUBLE PRECISION :: H2OINTERXFORCES(NATOMS,REDUCED_NUM_EIG)
	DOUBLE PRECISION :: H2OINTERYFORCES(NATOMS,REDUCED_NUM_EIG)
	DOUBLE PRECISION :: H2OINTERZFORCES(NATOMS,REDUCED_NUM_EIG)

!###################################################################
! Check that there is no Roo vectors smaller in magnitude than
! the RCut cutoff value
!  
! 	RCut = 2.0
!         do i = 3,num_atom-4,3
!            do j = i+3,num_atom-1,3
!               dx = psix(i) - psix(j)
!               dy = psiy(i) - psiy(j)
!               dz = psiz(i) - psiz(j)
!               R = DSQRT(dx**2 + dy**2 + dz**2)
!               if (R.lt.RCut) then
! 		  dQ(3*i+1) = -1.0d02
! 	          dQ(3*i+2) = -1.0d02
!                 dQ(3*i+3) = -1.0d02
!                 dQ(3*j+1) = -1.0d02
!                 dQ(3*j+2) = -1.0d02
!                 dQ(3*j+3) = -1.0d02
!                 RETURN
!               endif
!            enddo
!         enddo              
!
!#######################################################################################
! Evaluation of the force terms requires the coordinates to be in 
! vectors "psix,psiy" etc. However, Polak-Ribiere routines pass the
! vector with which the energy is to be minimized as the vector "p".
! 
!########################################################################################

! Not always necessary to recall the energy step
! The only things carried over are the VB assignments and the exchange energies

	if (assignVBstates) call msevb(ncoord,Q,assignVBstates,enew)

! Calculate H2O-H2O interactions all in one go

	call calculateH2OinterForces (h2oInterXforces, h2oInterYforces, h2oInterZforces)

!#######################################################################################
! Loop over all atoms
!
!####################################################################################### 
 	do currentAtom = 1,natoms

!#######################################################################################
! Loop Over the Eigenstates of System
!#######################################################################################

	   do state1 = 1, reduced_num_eig

! Which atom is this in the current VB state?

	      do j = 1, natoms
		 if (atmpl(state1,j).eq.currentAtom) then
		    positionInState1 = j
		    exit
		 endif
	      enddo

!###############################################################################################
! Determine the diagonal matrix element

	      call fdiag(positionInState1,state1,force_x,force_y,force_z)

	      fxH(state1,state1) = force_x + h2oInterXforces(currentAtom,state1)
	      fyH(state1,state1) = force_y + h2oInterYforces(currentAtom,state1)
	      fzH(state1,state1) = force_z + h2oInterZforces(currentAtom,state1)

!###############################################################################################
! Having calculated the diagonal element, now determine the off-diagonal elements

	      do state2 = (state1+1),reduced_num_eig

		 if (statesInteract(state1,state2)) then
		 
		    do j = 1, natoms
		       if (atmpl(state2,j).eq.currentAtom) then
			  positionInState2 = j
			  exit
		       endif
		    enddo  

! Already have determined the Zundel species for this combination

		    roh = interAtomicR(zundel_species(state1,state2,4), atmpl(state1,3))
    
		    if (currentAtom.eq.zundel_species(state1,state2,4)) then
		       dRoh(1) = psix(zundel_species(state1,state2,4)) - psix(atmpl(state1,3))
		       dRoh(2) = psiy(zundel_species(state1,state2,4)) - psiy(atmpl(state1,3))
		       dRoh(3) = psiz(zundel_species(state1,state2,4)) - psiz(atmpl(state1,3))
		    else
		       dRoh(1) = psix(atmpl(state1,3)) - psix(zundel_species(state1,state2,4))
		       dRoh(2) = psiy(atmpl(state1,3)) - psiy(zundel_species(state1,state2,4))
		       dRoh(3) = psiz(atmpl(state1,3)) - psiz(zundel_species(state1,state2,4))            
		    endif

		    if (currentAtom.eq.atmpl(state2,3)) then
		       dRoo(1) = psix(atmpl(state2,3)) - psix(atmpl(state1,3))
		       dRoo(2) = psiy(atmpl(state2,3)) - psiy(atmpl(state1,3))
		       dRoo(3) = psiz(atmpl(state2,3)) - psiz(atmpl(state1,3))  
		    else
		       dRoo(1) = psix(atmpl(state1,3)) - psix(atmpl(state2,3))
		       dRoo(2) = psiy(atmpl(state1,3)) - psiy(atmpl(state2,3))
		       dRoo(3) = psiz(atmpl(state1,3)) - psiz(atmpl(state2,3))  
		    endif

		    Roo = interAtomicR(atmpl(state1,3), atmpl(state2,3))

!###############################################################################################
! Determine the appropriate off-diagonal routine to call
!

		    if (MOD(currentAtom,3).eq.0) then
		       if ((positionInState1.eq.3).OR.(positionInState2.eq.3)) then

! Hydronium O atom in either of the VB states under consideration
			  call offd_atom3(currentAtom,Roo,Roh,dRoo,dRoh,Fxatom3,Fyatom3,Fzatom3,state1,state2)
			  offdiagx = Fxatom3
			  offdiagy = Fyatom3
			  offdiagz = Fzatom3 

		       else	
! Pure H2O O atom
			  call offd_atom1(currentAtom,Roo,Roh,Fxatom1,Fyatom1,Fzatom1,state1,state2)
			  offdiagx = Fxatom1
			  offdiagy = Fyatom1
			  offdiagz = Fzatom1
		       endif
!---		       elseif (atmpl(i,4).eq.4) then
!---		       elseif (position.eq.hmin) then
		    elseif (currentAtom.eq.zundel_species(state1,state2,4)) then

! Exchanging H atom
		       call offd_atom4(currentAtom,Roo,Roh,dRoo,dRoh,Fxatom1,Fyatom1,Fzatom1,state1,state2)
		       offdiagx = Fxatom1
		       offdiagy = Fyatom1
		       offdiagz = Fzatom1 
		    else

! Either a pure water H or a non-exchanging Zundel H

		       call offd_atom1(currentAtom,Roo,Roh,Fxatom1,Fyatom1,Fzatom1,state1,state2)
		       offdiagx = Fxatom1
		       offdiagy = Fyatom1
		       offdiagz = Fzatom1
		    endif
		    
		    fxH(state1,state2) = offdiagx
		    fyH(state1,state2) = offdiagy
		    fzH(state1,state2) = offdiagz
		    fxH(state2,state1) = offdiagx
		    fyH(state2,state1) = offdiagy
		    fzH(state2,state1) = offdiagz

! Debugging
!		       off_diag_forces(state1,state2,3*jj-2)=offdiagx
!		       off_diag_forces(state1,state2,3*jj-1)=offdiagy
!		       off_diag_forces(state1,state2,3*jj)=offdiagz
		 
		 else
		    fxH(state1,state2) = 0.0d0
		    fyH(state1,state2) = 0.0d0
		    fzH(state1,state2) = 0.0d0
		    fxH(state2,state1) = 0.0d0
		    fyH(state2,state1) = 0.0d0
		    fzH(state2,state1) = 0.0d0
		 endif
	      enddo

!----Sum up the components to make the element Hij

	   enddo

! Sum components

!#####################################################################
! Assign the ground state eigenvector as corresponding to the 
! lowest eigenvalue
!  *NB* No consideration of sign is observed here
!
!#####################################################################

	   minimum_x(currentAtom) = 0.0d0
	   minimum_y(currentAtom) = 0.0d0
	   minimum_z(currentAtom) = 0.0d0

	   do i = 1,reduced_num_eig
	      do j = 1,reduced_num_eig
		 minimum_x(currentAtom) = minimum_x(currentAtom) + grstwfu(i)*grstwfu(j)*fxH(i,j)
		 minimum_y(currentAtom) = minimum_y(currentAtom) + grstwfu(i)*grstwfu(j)*fyH(i,j)
		 minimum_z(currentAtom) = minimum_z(currentAtom) + grstwfu(i)*grstwfu(j)*fzH(i,j)
	      enddo
	   enddo

!	   print *, 'Force matrix: atom', currentAtom
!	   do i=1,reduced_num_eig
!	   do j=1,reduced_num_eig
!	      print *, '***', i, j, '***'
!	      print *, '***', i, '***'
!	      write(*, '(3(a3,F14.6))') 'x:', fxH(i,i), 'y:', fyH(i,i), 'z:', fzH(i,i) 
!	   enddo
!	   enddo

!#####################################################################
! End loop over Atoms
!
!#####################################################################
 	enddo

	j = 1
	do i = 1,natoms
	   dQ(j) = minimum_x(i)
	   dQ(j+1) = minimum_y(i)
	   dQ(j+2) = minimum_z(i)
	   j = j + 3
	enddo

!	stop

	RETURN 
        END

  !#####################################################################
!############### SUBROUTINE FDIAG ####################################
!#####################################################################

! H20-H2O intermolecular forces are calculated separately

	subroutine fdiag(atompos, vbState, force_x, force_y, force_z)

	use msevb_common

! atompos is the number of the atom in VBstate, the VB state of interest

	implicit none

!-----------------------------------------------------------------
! 23-7-00 R.A. Christie
! routine to call the diagonal Hellman-Feynman force calculations
!
! Includes the following routines:
! 	fh3ointra.f
!	finter.f
!	fh2ointra.f
!-----------------------------------------------------------------

!----Variables to be passed
        integer atompos,vbState
	DOUBLE PRECISION FH3OINTRAX,FH3OINTRAY,FH3OINTRAZ
	DOUBLE PRECISION FINTERX,FINTERY,FINTERZ
	DOUBLE PRECISION FORCE_X,FORCE_Y,FORCE_Z
	DOUBLE PRECISION FH2OINTRAX,FH2OINTRAY,FH2OINTRAZ

!-----------------------------------------------------------------
! Initialize variables
!
!----------------------------------------------------------------- 

	fh3ointrax = 0.0d0
	fh3ointray = 0.0d0
	fh3ointraz = 0.0d0
	finterx = 0.0d0
	fintery = 0.0d0
	finterz = 0.0d0
	fh2ointrax = 0.0d0
	fh2ointray = 0.0d0
	fh2ointraz = 0.0d0
	force_x = 0.0d0
	force_y = 0.0d0
	force_z = 0.0d0

!-----------------------------------------------------------------
! Call subroutines according to the atom in question
!-----------------------------------------------------------------

        if (atompos.le.4) then  ! Hydronium atom in this state

           call fh3ointra(atompos,vbState,fh3ointrax,fh3ointray,fh3ointraz)
           call finter(atompos,vbState,finterx,fintery,finterz)

	   force_x = fh3ointrax + finterx 
	   force_y = fh3ointray + fintery
           force_z = fh3ointraz + finterz 

        else  ! Water atom in this state

           call fh2ointra(atompos,vbState,fh2ointrax,fh2ointray,fh2ointraz)
           call finter(atompos,vbState,finterx,fintery,finterz)

	   force_x = fh2ointrax + finterx
	   force_y = fh2ointray + fintery
	   force_z = fh2ointraz + finterz
        endif

        RETURN
        END

!#####################################################################
!################### SUBROUTINE FH3OINTRA ############################
!#####################################################################

        subroutine fh3ointra(atompos,vbState,fh3ointrax,fh3ointray,fh3ointraz)

	use msevb_common

	implicit none

!-----------------------------------------------------------------
! 30-7-99 R.A. Christie
!
!   -- Update: 21-07-00
! Subroutine to determine the forces from the intramolecular bend
! and stretch of the H3O+ molecule
!
!-----------------------------------------------------------------

!----Incoming variables
        integer atompos,vbState
	DOUBLE PRECISION FH3OINTRAX,FH3OINTRAY,FH3OINTRAZ

!----Local variables
	integer atomNo
        DOUBLE PRECISION DX,DY,DZ,EXPTERM,VROH(10),RFH3OSTR,FH3OSTR,RFH3OSTRX, &
             rfh3ostry,rfh3ostrz,roh(3),theta(3),dxroh(10),dyroh(10), &
             dzroh(10),rfh3obendx1,rfh3obendy1,rfh3obendz1,rohstr, &
             rfh3obendx,rfh3obendy,rfh3obendz,rfh3obendx2,AdotA,AdotB, &
      	     AdotC,rfh3obendy2,rfh3obendz2,fh3ostr_temp,dxoh(10),dyoh(10), &
      	     dzoh(10),vectSq(3),angleDiff(3),part1x,part2x,part3x,part4x,part1y, &
      	     part2y,part3y,part4y,part1z,part2z,part3z,part4z,angleDiff2, &
      	     vectSq2,theta_temp
        integer i,j,otherb1,otherb2,bnd,hyd
!	DOUBLE PRECISION RADS
!----New variables
	DOUBLE PRECISION ATOMPOSVECTOR,OTHERB1VECTOR,OTHERB2VECTOR,RFH3OBENDX_TEMP1, &
      		rfh3oBendX_temp2,rfh3oBendX_temp3,rfh3oBendY_temp1, &
      		rfh3oBendY_temp2,rfh3oBendY_temp3,rfh3oBendZ_temp1, &
      		rfh3oBendZ_temp2,rfh3oBendZ_temp3,tempR
	integer pivotOxygen
	DOUBLE PRECISION CHECKED_ACOS

!-----Initialize variables
	expterm = 0.0d0
	fh3ostr = 0.0d0
	rfh3ostrx = 0.0d0
	fh3ostr_temp = 0.0d0
	rfh3ostry = 0.0d0
	rfh3ostrz = 0.0d0
	rfh3obendx1 = 0.0d0
	rfh3obendy1 = 0.0d0
	rfh3obendz1 = 0.0d0
	rfh3obendx = 0.0d0
        rfh3obendy = 0.0d0
        rfh3obendz = 0.0d0
	do i = 1,10
	   dxroh(i) = 0.0d0
	   dxoh(i) = 0.0d0
	   dyroh(i) = 0.0d0
	   dyoh(i) = 0.0d0
	   dzroh(i) = 0.0d0
	   dzoh(i) = 0.d0
	enddo

!----Initialize new variables
	atomposVector = 0.0d0
	otherb1Vector = 0.0d0
	otherb2Vector = 0.0d0
	rfh3oBendX_temp1 = 0.0d0
	rfh3oBendX_temp2 = 0.0d0
	rfh3oBendX_temp3 = 0.0d0
	rfh3oBendY_temp1 = 0.0d0
        rfh3oBendY_temp2 = 0.0d0
        rfh3oBendY_temp3 = 0.0d0
	rfh3oBendZ_temp1 = 0.0d0
        rfh3oBendZ_temp2 = 0.0d0
        rfh3oBendZ_temp3 = 0.0d0
	do i = 1,3
           vectSq(i) = 0.0d0
	   angleDiff(i) = 0.0d0
	enddo

! The position of the atom in the original input vectors

	atomNo = atmpl(vbState, atompos)

!-----------------------------------------------------------------------
! Calculate the stretch contribution to the force first.
!
!-----------------------------------------------------------------------

        if (atompos.ne.3) then    ! H atom
           dx = psix(atomNo) - psix(atmpl(vbState,3))
           dy = psiy(atomNo) - psiy(atmpl(vbState,3))
           dz = psiz(atomNo) - psiz(atmpl(vbState,3))
	   rohstr = interAtomicR(atomNo,atmpl(vbState,3))
           expterm = DEXP(-1.0d0*aoheq*(rohstr-roheq))
           fh3ostr = (2.0d0/rohstr)*expterm*aoheq*doh*(1.0d0-expterm)
	   rfh3ostrx = dx*fh3ostr
	   rfh3ostry = dy*fh3ostr
           rfh3ostrz = dz*fh3ostr
        else                    ! O atom
	   j = 1
           do i = 1,3
              dxoh(i) = psix(atomNo) - psix(atmpl(vbState,j))
              dyoh(i) = psiy(atomNo) - psiy(atmpl(vbState,j))
              dzoh(i) = psiz(atomNo) - psiz(atmpl(vbState,j))
	      vroh(i) = interAtomicR(atomNo,atmpl(vbState,j))
              expterm = DEXP(-1.0d0*aoheq*(vroh(i)-roheq))
              fh3ostr_temp = (2.0d0/vroh(i))*expterm*aoheq*doh*(1.0d0-expterm)
              rfh3ostrx = rfh3ostrx + fh3ostr_temp*dxoh(i)
 	      rfh3ostry = rfh3ostry + fh3ostr_temp*dyoh(i) 
 	      rfh3ostrz = rfh3ostrz + fh3ostr_temp*dzoh(i) 
              j = j+1
              if (j.eq.3) j = 4
           enddo
        endif

!---------------------------------------------------------------------
! Calculate the contribution to the force from the harmonic bending
! potential
!
! 1. Determine the "other atoms" which the atom in question is
!    interacting with.
! 2. Find the magnitude of the position vectors in a global coordinate
!    system.
! 3. Determine the magnitude of the O-H vectors in a local, molecule,
!    coordinate system.
! 4. Calculate the angle between the O and H vectors
! 5. Finally, calculate the bending component of force.
!
!---------------------------------------------------------------------

!----OXYGEN
        if (atompos.eq.3) then

	   dxoh(1) = psix(atmpl(vbState,1)) - psix(atomNo) 
	   dyoh(1) = psiy(atmpl(vbState,1)) - psiy(atomNo) 
	   dzoh(1) = psiz(atmpl(vbState,1)) - psiz(atomNo) 
	   vroh(1) = interAtomicR(atmpl(vbState,1),atomNo)

	   dxoh(2) = psix(atmpl(vbState,2)) - psix(atomNo) 
	   dyoh(2) = psiy(atmpl(vbState,2)) - psiy(atomNo) 
	   dzoh(2) = psiz(atmpl(vbState,2)) - psiz(atomNo) 
	   vroh(2) = interAtomicR(atmpl(vbState,2),atomNo)

	   dxoh(3) = psix(atmpl(vbState,4)) - psix(atomNo) 
	   dyoh(3) = psiy(atmpl(vbState,4)) - psiy(atomNo) 
	   dzoh(3) = psiz(atmpl(vbState,4)) - psiz(atomNo) 
	   vroh(3) = interAtomicR(atmpl(vbState,4),atomNo)

	   theta_temp = dxoh(1)*dxoh(2) + dyoh(1)*dyoh(2) &
              + dzoh(1)*dzoh(2)
	   theta(1) = checked_acos(theta_temp/(vroh(1)*vroh(2)))

	   theta_temp = dxoh(1)*dxoh(3) + dyoh(1)*dyoh(3) &
              + dzoh(1)*dzoh(3)
	   theta(2) = checked_acos(theta_temp/(vroh(1)*vroh(3)))

	   theta_temp = dxoh(3)*dxoh(2) + dyoh(3)*dyoh(2) &
              + dzoh(3)*dzoh(2)
	   theta(3) = checked_acos(theta_temp/(vroh(3)*vroh(2)))
 
	   angleDiff(1) = theta(1) - alphaeq
	   angleDiff(2) = theta(2) - alphaeq
	   angleDiff(3) = theta(3) - alphaeq

	   vectSq(1) = dxoh(1)*dxoh(2) + dyoh(1)*dyoh(2) + dzoh(1)*dzoh(2)
	   vectSq(2) = dxoh(1)*dxoh(3) + dyoh(1)*dyoh(3) + dzoh(1)*dzoh(3)
	   vectSq(3) = dxoh(2)*dxoh(3) + dyoh(2)*dyoh(3) + dzoh(2)*dzoh(3)
 
!---------------------------------------------------------------------
!  X
!---------------------------------------------------------------------
	   part1x = (dxoh(1)*vectSq(1))/(vroh(2)*vroh(1)**3) + (dxoh(2)*vectSq(1))/(vroh(1)*vroh(2)**3)
	   part2x = dxoh(1)/(vroh(1)*vroh(2)) + dxoh(2)/(vroh(1)*vroh(2))
	   part3x = -1.0d0*DSIN(theta(1))

	   if (DABS(part3x)<numericalZeroLimit) then
	      rfh3obendx_temp1 = large_force
	   else
	      rfh3obendx_temp1 = kalpha*(theta(1) - alphaeq)*(part1x - part2x)/part3x
	   endif

	   rfh3obendx = rfh3obendx + rfh3obendx_temp1
	     
	   part1x = (dxoh(1)*vectSq(2))/(vroh(3)*vroh(1)**3) + (dxoh(3)*vectSq(2))/(vroh(1)*vroh(3)**3)
	   part2x = dxoh(1)/(vroh(1)*vroh(3)) + dxoh(3)/(vroh(1)*vroh(3))
	   part3x = -1.0d0*DSIN(theta(2))

	   if (DABS(part3x)<numericalZeroLimit) then
	      rfh3obendx_temp1 = large_force
	   else 
	      rfh3obendx_temp1 = kalpha*(theta(2) - alphaeq)*(part1x - part2x)/part3x
	   endif

	   rfh3obendx = rfh3obendx + rfh3obendx_temp1 

	   part1x = (dxoh(2)*vectSq(3))/(vroh(3)*vroh(2)**3) + (dxoh(3)*vectSq(3))/(vroh(2)*vroh(3)**3)
	   part2x = dxoh(2)/(vroh(2)*vroh(3)) + dxoh(3)/(vroh(2)*vroh(3))
	   part3x = -1.0d0*DSIN(theta(3))

	   if (DABS(part3x)<numericalZeroLimit) then
	      rfh3obendx_temp1 = large_force
	   else     
	      rfh3obendx_temp1 = kalpha*(theta(3) - alphaeq)*(part1x - part2x)/part3x
	   endif
	   
	   rfh3obendx = rfh3obendx + rfh3obendx_temp1 

!---------------------------------------------------------------------
!  Y
!---------------------------------------------------------------------
	   part1y = -1.0d0*dyoh(1)*vectSq(1)/(vroh(2)*vroh(1)**3)
	   part2y = (dyoh(2)+dyoh(1))/(vroh(2)*vroh(1))
	   part3y = -1.0d0*dyoh(2)*vectSq(1)/(vroh(1)*vroh(2)**3)
	   part4y = 1.0d0 - vectSq(1)**2/((vroh(2)**2)*(vroh(1)**2))

	   if (DABS(part4y)<numericalZeroLimit) then
	      rfh3obendy_temp1 = large_force
	   else
	      rfh3obendy_temp1 = (kalpha*(part1y + part2y + part3y)*angleDiff(1))/DSQRT(part4y)
	   endif
		 
	   rfh3obendy = rfh3obendy + rfh3obendy_temp1
	   
	   part1y = -1.0d0*dyoh(2)*vectSq(3)/(vroh(3)*vroh(2)**3)
	   part2y = (dyoh(3)+dyoh(2))/(vroh(3)*vroh(2))
	   part3y = -1.0d0*dyoh(3)*vectSq(3)/(vroh(2)*vroh(3)**3)
	   part4y = 1.0d0 - vectSq(3)**2/((vroh(3)**2)*(vroh(2)**2))
	      
	   if (DABS(part4y)<numericalZeroLimit) then
	      rfh3obendy_temp1 = large_force
	   else
	      rfh3obendy_temp1 = (kalpha*(part1y + part2y + part3y)*angleDiff(3))/DSQRT(part4y)
	   endif
	   
	   rfh3obendy = rfh3obendy + rfh3obendy_temp1 

	   part1y = -1.0d0*dyoh(3)*vectSq(2)/(vroh(1)*vroh(3)**3)
	   part2y = (dyoh(1)+dyoh(3))/(vroh(1)*vroh(3))
	   part3y = -1.0d0*dyoh(1)*vectSq(2)/(vroh(3)*vroh(1)**3)
	   part4y = 1.0d0 - vectSq(2)**2/((vroh(1)**2)*(vroh(3)**2))

	   if (DABS(part4y)<numericalZeroLimit) then
	      rfh3obendy_temp1 = large_force
	   else
	      rfh3obendy_temp1 = (kalpha*(part1y + part2y + part3y)*angleDiff(2))/DSQRT(part4y)
	   endif

	   rfh3obendy = rfh3obendy + rfh3obendy_temp1 
!---------------------------------------------------------------------
!  Z
!---------------------------------------------------------------------
	   part1z = -1.0d0*dzoh(1)*vectSq(1)/(vroh(2)*vroh(1)**3)
	   part2z = (dzoh(2)+dzoh(1))/(vroh(2)*vroh(1))
	   part3z = -1.0d0*dzoh(2)*vectSq(1)/(vroh(1)*vroh(2)**3)
	   part4z = 1.0d0 - vectSq(1)**2/((vroh(2)**2)*(vroh(1)**2))

	   if (DABS(part4z)<numericalZeroLimit) then
	      rfh3obendz_temp1 = large_force
	   else
	      rfh3obendz_temp1 = (kalpha*(part1z + part2z + part3z)*angleDiff(1))/DSQRT(part4z)
	   endif
	   
	   rfh3obendz = rfh3obendz + rfh3obendz_temp1
	   
	   part1z = -1.0d0*dzoh(2)*vectSq(3)/(vroh(3)*vroh(2)**3)
	   part2z = (dzoh(3)+dzoh(2))/(vroh(3)*vroh(2))
	   part3z = -1.0d0*dzoh(3)*vectSq(3)/(vroh(2)*vroh(3)**3)
	   part4z = 1.0d0 - vectSq(3)**2/((vroh(3)**2)*(vroh(2)**2))

	   if (DABS(part4z)<numericalZeroLimit) then
	      rfh3obendz_temp1 = large_force
	   else
	      rfh3obendz_temp1 = (kalpha*(part1z + part2z + part3z)*angleDiff(3))/DSQRT(part4z)
	   endif
	      
	   rfh3obendz = rfh3obendz + rfh3obendz_temp1   

	   part1z = -1.0d0*dzoh(3)*vectSq(2)/(vroh(1)*vroh(3)**3)
	   part2z = (dzoh(1)+dzoh(3))/(vroh(1)*vroh(3))
	   part3z = -1.0d0*dzoh(1)*vectSq(2)/(vroh(3)*vroh(1)**3)
	   part4z = 1.0d0 - vectSq(2)**2/((vroh(1)**2)*(vroh(3)**2))

	   if (DABS(part4z)<numericalZeroLimit) then
	      rfh3obendz_temp1 = large_force
	   else
	      rfh3obendz_temp1 = (kalpha*(part1z + part2z + part3z)*angleDiff(2))/DSQRT(part4z)
	   endif

	   rfh3obendz = rfh3obendz + rfh3obendz_temp1   
 	     
	else

!---------------------------------------------------------------------
!  HYDROGEN
!---------------------------------------------------------------------   
           if (atompos.eq.1) then
              otherb1 = atmpl(vbState,2)
              otherb2 = atmpl(vbState,4)
           elseif (atompos.eq.2) then
              otherb1 = atmpl(vbState,1)
              otherb2 = atmpl(vbState,4)
	   else
	      otherb1 = atmpl(vbState,1)
	      otherb2 = atmpl(vbState,2)
           endif

           do i = 1,3
              roh(i) = 0.0d0
              theta(i) = 0.0d0
           enddo

!----Determine the O-H distances in hydronium
           dxroh(1) = psix(atomNo) - psix(atmpl(vbState,3))
           dyroh(1) = psiy(atomNo) - psiy(atmpl(vbState,3))
           dzroh(1) = psiz(atomNo) - psiz(atmpl(vbState,3))
	   roh(1) = interAtomicR(atomNo,atmpl(vbState,3))

           dxroh(2) = psix(otherb1) - psix(atmpl(vbState,3))
           dyroh(2) = psiy(otherb1) - psiy(atmpl(vbState,3))
           dzroh(2) = psiz(otherb1) - psiz(atmpl(vbState,3))
	   roh(2) = interAtomicR(otherb1,atmpl(vbState,3))

           dxroh(3) = psix(otherb2) - psix(atmpl(vbState,3))
           dyroh(3) = psiy(otherb2) - psiy(atmpl(vbState,3))
           dzroh(3) = psiz(otherb2) - psiz(atmpl(vbState,3))
	   roh(3) = interAtomicR(otherb2,atmpl(vbState,3))

!----Determine the bending angles, alpha
           theta(1) = dxroh(1)*dxroh(2)+dyroh(1)*dyroh(2) &
                        + dzroh(1)*dzroh(2)
	   theta(1) = checked_acos(theta(1)/(roh(1)*roh(2)))

           theta(2) = dxroh(3)*dxroh(1)+dyroh(3)*dyroh(1) &
                        + dzroh(3)*dzroh(1)
           theta(2) = checked_acos(theta(2)/(roh(3)*roh(1)))

	   vectSq(1) = dxroh(1)*dxroh(2) + dyroh(1)*dyroh(2) + dzroh(1)*dzroh(2)
	   vectSq(2) = dxroh(1)*dxroh(3) + dyroh(1)*dyroh(3) + dzroh(1)*dzroh(3)
	   angleDiff(1) = theta(1) - alphaeq
	   angleDiff(2) = theta(2) - alphaeq

!---------------------------------------------------------------------
!  X
!---------------------------------------------------------------------
!----Consider the first angle
	   part1x = dxroh(1)*vectSq(1)/(roh(2)*roh(1)**3)
	   part2x = dxroh(2)/(Roh(1)*Roh(2))
	   part3x = 1.0d0 - (vectSq(1)/(roh(1)*roh(2)))**2 ! = (sin(alpha))^2

! Check, there is a discontinuity in the potential/force when alpha = 180

	   if (DABS(part3x)<numericalZeroLimit) then
	      rfh3obendx_temp1 = large_force
	   else
	      rfh3obendx_temp1 = (kalpha*(part1x - part2x)*angleDiff(1))/DSQRT(part3x)
	   endif

!----Consider the Second angle
	   part1x = dxroh(1)*vectSq(2)/(roh(3)*roh(1)**3)
	   part2x = dxroh(3)/(Roh(1)*Roh(3))
	   part3x = 1.0d0 - vectSq(2)**2/((roh(1)**2)*(roh(3)**2))

	   if (DABS(part3x)<numericalZeroLimit) then
	      rfh3obendx_temp1 = large_force
	   else
	      rfh3obendx_temp2 = (kalpha*(part1x - part2x)*angleDiff(2))/DSQRT(part3x)
	   endif

!----Sum up the contributions
	   rfh3obendx = rfh3obendx_temp1 + rfh3obendx_temp2
!---------------------------------------------------------------------
!  Y
!---------------------------------------------------------------------
!----Consider the first angle
	   part1y = dyroh(1)*vectSq(1)/(roh(2)*roh(1)**3)
	   part2y = dyroh(2)/(Roh(1)*Roh(2))
	   part3y = 1.0d0 - vectSq(1)**2/((roh(1)**2)*(roh(2)**2))

	   if (DABS(part3y)<numericalZeroLimit) then
	      rfh3obendy_temp1 = large_force
	   else
	      rfh3obendy_temp1 = (kalpha*(part1y - part2y)*angleDiff(1))/DSQRT(part3y)
	   endif

!----Consider the Second angle
	   part1y = dyroh(1)*vectSq(2)/(roh(3)*roh(1)**3)
	   part2y = dyroh(3)/(Roh(1)*Roh(3))
	   part3y = 1.0d0 - vectSq(2)**2/((roh(1)**2)*(roh(3)**2))

	   if (DABS(part3y)<numericalZeroLimit) then
	      rfh3obendy_temp2 = large_force
	   else
	      rfh3obendy_temp2 = (kalpha*(part1y - part2y)*angleDiff(2))/DSQRT(part3y)
	   endif

!----Sum up the contributions
	   rfh3obendy = rfh3obendy_temp1 + rfh3obendy_temp2  

!---------------------------------------------------------------------
!  Z
!---------------------------------------------------------------------
!----Consider the first angle
	   part1z = dzroh(1)*vectSq(1)/(roh(2)*roh(1)**3)
	   part2z = dzroh(2)/(Roh(1)*Roh(2))
	   part3z = 1.0d0 - vectSq(1)**2/((roh(1)**2)*(roh(2)**2))

	   if (DABS(part3z)<numericalZeroLimit) then
	      rfh3obendz_temp1 = large_force
	   else
	      rfh3obendz_temp1 = (kalpha*(part1z - part2z)*angleDiff(1))/DSQRT(part3z)
	   endif

!----Consider the Second angle
	   part1z = dzroh(1)*vectSq(2)/(roh(3)*roh(1)**3)
	   part2z = dzroh(3)/(Roh(1)*Roh(3))
	   part3z = 1.0d0 - vectSq(2)**2/((roh(1)**2)*(roh(3)**2))

	   if (DABS(part3z)<numericalZeroLimit) then
	      rfh3obendz_temp2 = large_force
	   else
	      rfh3obendz_temp2 = (kalpha*(part1z - part2z)*angleDiff(2))/DSQRT(part3z)
	   endif

!----Sum up the contributions
	   rfh3obendz = rfh3obendz_temp1 + rfh3obendz_temp2  

!----end of the conditional statement (if(atompos.ne.3)---
        endif

!---------------------------------------------------------------------
! Finally sum the force components to get the total force
!
!---------------------------------------------------------------------

        fh3ointrax = rfh3obendx + rfh3ostrx
        fh3ointray = rfh3obendy + rfh3ostry
        fh3ointraz = rfh3obendz + rfh3ostrz

! Debugging

!	h3o_intra_forces(vbState,(3*atomNo)-2) = fh3ointrax
!	h3o_intra_forces(vbState,(3*atomNo)-1) = fh3ointray
!	h3o_intra_forces(vbState,(3*atomNo)) = fh3ointraz	

        RETURN
        END

!#####################################################################
!############## SUBROUTINE FINTER ####################################
!#####################################################################

	subroutine finter(atompos,vbState,finterx,fintery,finterz)

	use commons
	use msevb_common

	implicit none

!---------------------------------------------------------------------------
! 06/08/99 R.A. Christie
! routine to determine the forces resulting from the H3O+ - H2O interaction
!
! Need as input into this routine: 	atompos (atom # in question)
!
! Type scenarios:
!	1. Atom in Q is hydrogen atom; only has coulomb interaction.
!	2. Atom in Q is a oxygen atom; coulomb, repulse and LJ interaction.
!---------------------------------------------------------------------------

!----Incoming variables
	DOUBLE PRECISION FINTERX,FINTERY,FINTERZ
	integer atompos,vbState,atomNo
!----Local variables 
	DOUBLE PRECISION QI,QJ,RFCOULOMB,DX,DY,DZ,RIJ,TRIG,RFREPULSE,FREPULSE, &
      	  	rfLJ,fLJ,Roo,Roo2,Roo8,Roo14,fcoulombx,fcoulomby,fcoulombz	
	integer i,j
!----New variables:
	DOUBLE PRECISION FLJX,FLJY,FLJZ,FREPULSEX,FREPULSEY,FREPULSEZ,ROO_SIGMA
	integer wheech

! Debugging variables

	integer coord_pos

!----initialize variables
	rfcoulomb = 0.0d0
	fcoulombx = 0.0d0
	fcoulomby = 0.0d0
	fcoulombz = 0.0d0
	fLJx = 0.0d0
	fLJy = 0.0d0
	fLJz = 0.0d0
	rfLJ = 0.0d0
	frepulsex = 0.0d0
	frepulsey = 0.0d0
	frepulsez = 0.0d0
	dx = 0.0d0
	dy = 0.0d0
	dz = 0.0d0
	frepulse = 0.0d0
	wheech = 1
	roo_sigma = 0.0d0

! Number of the atom in the original input vectors

	atomNo = atmpl(vbState,atompos)

!---------------------------------------------------------------------
! Determine the Coulombic part of the interaction first
!
! Two possibilities:
!   1. Atom is in hydronium unit
!   2. Atom is in solvating water molecule
!
! This part of the code calculates the Coulombic interaction of both
! oxygen and hydrogen atoms.
!
! ** NB ** units for Coulombic  interaction should be kcal/angstroms
!
!---------------------------------------------------------------------

	if (atompos.le.4) then
!----case 1. When atom Q is in hydronium
	   if (atompos.eq.3) then
	      qi = qintero
	   else
	      qi = qinterh
	   endif
	   do j = 5,natoms
	      wheech = atmpl(vbState,j)
	      if (MOD(wheech,3).eq.0) then
                 qj = qtip3p_o
              else
                 qj = qtip3p_h
              endif
              dx = psix(atomNo) - psix(wheech)
              dy = psiy(atomNo) - psiy(wheech)
              dz = psiz(atomNo) - psiz(wheech)
	      rij = interAtomicR(atomNo,wheech)
 	      rfcoulomb = -1.0d0*dampchge*(qi*qj*fourpieo)/(rij**3)
	      fcoulombx = rfcoulomb*dx + fcoulombx
	      fcoulomby = rfcoulomb*dy + fcoulomby
	      fcoulombz = rfcoulomb*dz + fcoulombz
           enddo 
	else

!----case 2. When atom in Q is in water molecule
	   if (MOD(atompos,3).eq.0) then
              qi = qtip3p_o
           else
              qi = qtip3p_h
           endif
	   do j = 1,4
	      wheech = atmpl(vbState,j)
              if (MOD(j,3).eq.0) then
                 qj = qintero
              else
                 qj = qinterh
              endif
              dx = psix(atomNo) - psix(wheech)
              dy = psiy(atomNo) - psiy(wheech)
              dz = psiz(atomNo) - psiz(wheech)
	      rij = interAtomicR(atomNo,wheech)
   	      rfcoulomb = -1.0d0*dampchge*(qi*qj*fourpieo)/(rij**3)
              fcoulombx = rfcoulomb*dx + fcoulombx
	      fcoulomby = rfcoulomb*dy + fcoulomby
	      fcoulombz = rfcoulomb*dz + fcoulombz
	   enddo
	endif

!---------------------------------------------------------------------
! Calculate the LJ and Repulsive interactions
!
! This part of the code only applies to Oxygen atoms, and applies to
! both Hydronium and solvating Water molecule Oxygen atoms
!
!
!---------------------------------------------------------------------

	if (MOD(atompos,3).eq.0) then
!----------------Oxygen atom on water molecule
	   if (atompos.ne.3) then
	       dx = psix(atomNo) - psix(atmpl(vbState,3))
	       dy = psiy(atomNo) - psiy(atmpl(vbState,3))
	       dz = psiz(atomNo) - psiz(atmpl(vbState,3))
	       Roo = interAtomicR(atomNo,atmpl(vbState,3))
 	       Roo_sigma = sigma_mix/Roo
	       Roo14 = Roo_sigma**14
	       Roo8 = Roo_sigma**8
 	       rfLJ = (-6.0d0/sigma_mix**2)*epsilon_mix*(2.0d0*Roo14 - Roo8)
	       fLJx = rfLJ*dx
	       fLJy = rfLJ*dy
	       fLJz = rfLJ*dz	
	       trig = DTANH(small_b*(Roo-dooeq))
	       frepulse = (-1.0d00/Roo)*small_b*big_b*(1.0d0 - trig**2)
	       frepulsex = frepulse*dx
	       frepulsey = frepulse*dy
	       frepulsez = frepulse*dz
	   else
!----------------Oxygen atom on hydronium ion
	       do i = 6,(natoms-1),3
 		  dx = psix(atomNo) - psix(atmpl(vbState,i))
 	          dy = psiy(atomNo) - psiy(atmpl(vbState,i))
 	 	  dz = psiz(atomNo) - psiz(atmpl(vbState,i))
		  Roo = interAtomicR(atomNo,atmpl(vbState,i))
 	          Roo_sigma = sigma_mix/Roo
                  Roo8 = Roo_sigma**8
                  Roo14 = Roo_sigma**14 
                  rfLJ = (-6.0d0/sigma_mix**2)*epsilon_mix*(2.0d0*Roo14 - Roo8)
                  trig = DTANH(small_b*(Roo-dooeq))
		  rfrepulse = (-1.0d00/Roo)*small_b*big_b*(1.0d0 - trig**2)
	          fLJx = fLJx + rfLJ*dx
		  fLJy = fLJy + rfLJ*dy
                  fLJz = fLJz + rfLJ*dz
	          frepulsex = frepulsex + rfrepulse*dx
		  frepulsey = frepulsey + rfrepulse*dy
		  frepulsez = frepulsez + rfrepulse*dz
	      enddo
	   endif
!----End of Oxygen-Atom Case
	endif

 	if (MOD(atompos,3).eq.0) then
 	   finterx = 1.0d0*(fLJx+frepulsex+fcoulombx)
 	   fintery = 1.0d0*(fLJy+frepulsey+fcoulomby)
 	   finterz = 1.0d0*(fLJz+frepulsez+fcoulombz)
 	else
 	   finterx = 1.0d0*fcoulombx
 	   fintery = 1.0d0*fcoulomby
 	   finterz = 1.0d0*fcoulombz
 	endif

! Debugging

!	coord_pos=3*atomNo-2

!	LJ_forces(eigenstate,coord_pos) = fLJx
!	LJ_forces(eigenstate,coord_pos+1) = fLJy
!	LJ_forces(eigenstate,coord_pos+2) = fLJz
!	coulomb_forces(eigenstate,coord_pos) = fcoulombx
!	coulomb_forces(eigenstate,coord_pos+1) = fcoulomby
!	coulomb_forces(eigenstate,coord_pos+2) = fcoulombz
!	rep_forces(eigenstate,coord_pos) = frepulsex
!	rep_forces(eigenstate,coord_pos+1) = frepulsey
!	rep_forces(eigenstate,coord_pos+2) = frepulsez

!	h3o_h2o_forces(vbState,coord_pos) = finterx
!	h3o_h2o_forces(vbState,coord_pos+1) = fintery
!	h3o_h2o_forces(vbState,coord_pos+2) = finterz	

	RETURN
	END

!#####################################################################

! Calculate the H2O-H2O interaction forces all in one go

	subroutine calculateH2OinterForces (h2oInterXforces, h2oInterYforces, h2oInterZforces)

	use commons
	use msevb_common

	implicit none

! Subroutine arguments

	DOUBLE PRECISION, INTENT(OUT) :: H2OINTERXFORCES(NATOMS, REDUCED_NUM_EIG)
	DOUBLE PRECISION, INTENT(OUT) :: H2OINTERYFORCES(NATOMS, REDUCED_NUM_EIG)
	DOUBLE PRECISION, INTENT(OUT) :: H2OINTERZFORCES(NATOMS, REDUCED_NUM_EIG)

! Local variables

	DOUBLE PRECISION :: XFORCES(NATOMS, NATOMS)
	DOUBLE PRECISION :: YFORCES(NATOMS, NATOMS)
	DOUBLE PRECISION :: ZFORCES(NATOMS, NATOMS)
	
	integer :: atom1, atom2, currentOatom, Hatom1, Hatom2
	integer :: currentVBstate
	DOUBLE PRECISION :: CHARGEONATOM1, CHARGEONATOM2
	DOUBLE PRECISION :: DX, DY, DZ
	DOUBLE PRECISION :: COULOMBPREFACTOR
	DOUBLE PRECISION :: INVROO2, SIGMAINVROO6, LJPREFACTOR
	DOUBLE PRECISION, PARAMETER :: LJCONSTANT = -24.0D0*H2OINTEREPSILON
	DOUBLE PRECISION, PARAMETER :: H2OINTERSIGMA_SQ = H2OINTERSIGMA**2

! Initialise

	xForces = 0.0
	yForces = 0.0
	zForces = 0.0

	h2oInterXforces = 0.0
	h2oInterYforces = 0.0
	h2oInterZforces = 0.0

! Calculate the interactions

	do atom1 = 1, (natoms-1)

	   if (MOD(atom1,3).eq.0) then
	      chargeOnAtom1 = qtip3p_o
	   else
	      chargeOnAtom1 = qtip3p_h
	   endif

	   do atom2 = (atom1+1), natoms

	      if (MOD(atom2,3).eq.0) then
		 chargeOnAtom2 = qtip3p_o
	      else
		 chargeOnAtom2 = qtip3p_h
	      endif

! Determine the coulombic part of the force

	      dx = psix(atom1) - psix(atom2)   ! Defined this way around for consistency
	      dy = psiy(atom1) - psiy(atom2)
	      dz = psiz(atom1) - psiz(atom2)      

	      coulombPrefactor = -1.0d0*fourpieo*chargeOnAtom1*chargeOnAtom2/((interAtomicR(atom1,atom2))**3)

! xForces(i,j) is the force on atom i due to presence of atom j
! xForces(j,i) = -xForces(i,j) because Newton`s law holds

	      xForces(atom1,atom2) = coulombPrefactor * dx
	      yForces(atom1,atom2) = coulombPrefactor * dy
	      zForces(atom1,atom2) = coulombPrefactor * dz

! Determine the LJ part of the force if necessary

	      if (chargeOnAtom1.eq.qtip3p_o .and. chargeOnAtom2.eq.qtip3p_o) then
		 invRoo2 = 1.0d0/(interAtomicR(atom1,atom2)**2)
		 sigmaInvRoo6 = h2ointerSigma_sq*invRoo2
		 sigmaInvRoo6 = sigmaInvRoo6**3

		 LJprefactor = LJconstant*(2.0d0*(sigmaInvRoo6**2) - sigmaInvRoo6)*invRoo2

		 xForces(atom1,atom2) = xForces(atom1,atom2) + LJprefactor*dx
		 yForces(atom1,atom2) = yForces(atom1,atom2) + LJprefactor*dy
		 zForces(atom1,atom2) = zForces(atom1,atom2) + LJprefactor*dz
	      endif

! Fill in the lower triangle explicitly
! This makes it easier than trying to work out which number is lower all the time in the next section

	      xForces(atom2,atom1) = -1.0d0 *  xForces(atom1,atom2)
	      yForces(atom2,atom1) = -1.0d0 *  yForces(atom1,atom2)	      
	      zForces(atom2,atom1) = -1.0d0 *  zForces(atom1,atom2)	      

! Add the new forces to the accumulating totals

	      h2oInterXforces(atom1,1) = h2oInterXforces(atom1,1) + xForces(atom1,atom2)
	      h2oInterXforces(atom2,1) = h2oInterXforces(atom2,1) - xForces(atom1,atom2)	      

	      h2oInterYforces(atom1,1) = h2oInterYforces(atom1,1) + yForces(atom1,atom2)
	      h2oInterYforces(atom2,1) = h2oInterYforces(atom2,1) - yForces(atom1,atom2)

	      h2oInterZforces(atom1,1) = h2oInterZforces(atom1,1) + zForces(atom1,atom2)
	      h2oInterZforces(atom2,1) = h2oInterZforces(atom2,1) - zForces(atom1,atom2)
	   enddo
	enddo
	
! Copy across the accumulated totals

	do currentVBstate = 2, reduced_num_eig
	   h2oInterXforces(:,currentVBstate) =  h2oInterXforces(:,1)
	   h2oInterYforces(:,currentVBstate) =  h2oInterYforces(:,1)
	   h2oInterZforces(:,currentVBstate) =  h2oInterZforces(:,1)
	enddo

! Subtract for each VB state those interactions which are wrongly included in the totals

	do currentVBstate = 1, reduced_num_eig

! Firstly make zero the forces for those atoms in the hydronium species

	   do atom1 = 1, 4
	      h2oInterXforces(atmpl(currentVBstate,atom1),currentVBstate) = 0.0	      
	      h2oInterYforces(atmpl(currentVBstate,atom1),currentVBstate) = 0.0	      
	      h2oInterZforces(atmpl(currentVBstate,atom1),currentVBstate) = 0.0	      
	   enddo

! Other atoms

	   do atom1 = 6, (natoms-1), 3

	      currentOatom = atmpl(currentVBstate,atom1)
	      Hatom1 = atmpl(currentVBstate,atom1-1)
	      Hatom2 = atmpl(currentVBstate,atom1+1)

! Remove interactions between other atoms and the hydronium atom

	      do atom2 = 1,4
! O-atom
		 h2oInterXforces(currentOatom,currentVBstate) =   &
      		 h2oInterXforces(currentOatom,currentVBstate) - xForces(currentOatom, atmpl(currentVBstate,atom2))
  		 h2oInterYforces(currentOatom,currentVBstate) =   &
      		 h2oInterYforces(currentOatom,currentVBstate) - yForces(currentOatom, atmpl(currentVBstate,atom2))
		 h2oInterZforces(currentOatom,currentVBstate) =   &
      		 h2oInterZforces(currentOatom,currentVBstate) - zForces(currentOatom, atmpl(currentVBstate,atom2))

! H-atom 1
		 h2oInterXforces(Hatom1,currentVBstate) =   &
      		 h2oInterXforces(Hatom1,currentVBstate) - xForces(Hatom1, atmpl(currentVBstate,atom2))
		 h2oInterYforces(Hatom1,currentVBstate) =   &
      		 h2oInterYforces(Hatom1,currentVBstate) - yForces(Hatom1, atmpl(currentVBstate,atom2))
		 h2oInterZforces(Hatom1,currentVBstate) =   &
      		 h2oInterZforces(Hatom1,currentVBstate) - zForces(Hatom1, atmpl(currentVBstate,atom2))

! H-atom 2
		 h2oInterXforces(Hatom2,currentVBstate) =  &
      		 h2oInterXforces(Hatom2,currentVBstate) - xForces(Hatom2, atmpl(currentVBstate,atom2))
		 h2oInterYforces(Hatom2,currentVBstate) =  &
      		 h2oInterYforces(Hatom2,currentVBstate) - yForces(Hatom2, atmpl(currentVBstate,atom2))
		 h2oInterZforces(Hatom2,currentVBstate) =  &
      		 h2oInterZforces(Hatom2,currentVBstate) - zForces(Hatom2, atmpl(currentVBstate,atom2))
	      enddo

! Remove interactions between atoms in the same water molecule

! O-atom
		 h2oInterXforces(currentOatom,currentVBstate) = h2oInterXforces(currentOatom,currentVBstate) - &
      		 xForces(currentOatom, Hatom1) - xForces(currentOatom, Hatom2)
		 h2oInterYforces(currentOatom,currentVBstate) = h2oInterYforces(currentOatom,currentVBstate) - &
      		 yForces(currentOatom, Hatom1) - yForces(currentOatom, Hatom2)
		 h2oInterZforces(currentOatom,currentVBstate) = h2oInterZforces(currentOatom,currentVBstate) - &
      		 zForces(currentOatom, Hatom1) - zForces(currentOatom, Hatom2)

! H-atom 1
		 h2oInterXforces(Hatom1,currentVBstate) = h2oInterXforces(Hatom1,currentVBstate) - &
      		 xForces(Hatom1, currentOatom) - xForces(Hatom1, Hatom2)
		 h2oInterYforces(Hatom1,currentVBstate) = h2oInterYforces(Hatom1,currentVBstate) - &
      		 yForces(Hatom1, currentOatom) - yForces(Hatom1, Hatom2)
		 h2oInterZforces(Hatom1,currentVBstate) = h2oInterZforces(Hatom1,currentVBstate) - &
      		 zForces(Hatom1, currentOatom) - zForces(Hatom1, Hatom2)

! H-atom 2
		 h2oInterXforces(Hatom2,currentVBstate) = h2oInterXforces(Hatom2,currentVBstate) - &
      		 xForces(Hatom2, currentOatom) - xForces(Hatom2, Hatom1)
		 h2oInterYforces(Hatom2,currentVBstate) = h2oInterYforces(Hatom2,currentVBstate) - &
      		 yForces(Hatom2, currentOatom) - yForces(Hatom2, Hatom1)
		 h2oInterZforces(Hatom2,currentVBstate) = h2oInterZforces(Hatom2,currentVBstate) - &
      		 zForces(Hatom2, currentOatom) - zForces(Hatom2, Hatom1)
	   enddo
	enddo

	return

	end subroutine

!#####################################################################
!###################### SUBROUTINE FH2OINTER #########################
!#####################################################################

	subroutine fh2ointer(atompos,vbState,fh2ointerx,fh2ointery,fh2ointerz)

	use commons
	use msevb_common

	implicit none

!---------------------------------------------------------------------------
!  06/08/99 R.A. Christie
! routine for determining the intermolecular water-water forces
!
!---------------------------------------------------------------------------

	DOUBLE PRECISION FH2OINTERX,FH2OINTERY,FH2OINTERZ
	integer atompos,vbState,atomNo
	DOUBLE PRECISION RFLJ,RFCOULOMB,FCOULOMBX,FCOULOMBY,FCOULOMBZ,QI,QJ,DX_OO, &
                dy_oo,dz_oo, &
      		Roo,Roo2,Roo8,Roo14,rij,dx,dy,dz,fLJx,fLJy,fLJz
	integer i,j,avoid2,avoid3

	DOUBLE PRECISION, PARAMETER :: LJPREFACTOR = -2.4D01*1.522D-01
	DOUBLE PRECISION :: COULOMBPREFACTOR

!----initialize variables
	coulombPrefactor = -1.0d0*fourpieo

	rfLJ = 0.0d0
	fLJx = 0.0d0
	fLJy = 0.0d0
	fLJz = 0.0d0
	rfcoulomb = 0.0d0
	fcoulombx = 0.0d0
	fcoulomby = 0.0d0
	fcoulombz = 0.0d0
	
	avoid2 = 0
	avoid3 = 0
	dx_oo = 0.0d0
	dy_oo = 0.0d0
	dz_oo = 0.0d0

! Position of the atom in the original input vector

	atomNo = atmpl(vbState,atompos)

!----find out which type of atom is: 

	if (MOD(atompos,3).eq.0) then	
	   qi = qtip3p_o
	   avoid2 = atompos - 1 
	   avoid3 = atompos + 1

!----detemine LJ part of force if considering an oxygen atom

           do j = 6,(natoms-1),3

	      if (j.eq.atompos) cycle

              dx_oo = psix(atomNo) - psix(atmpl(vbState,j))
              dy_oo = psiy(atomNo) - psiy(atmpl(vbState,j))
              dz_oo = psiz(atomNo) - psiz(atmpl(vbState,j))
	      Roo = DSQRT(dx_oo**2 + dy_oo**2 + dz_oo**2)
	      Roo2 = Roo**2
	      Roo = 3.1506d0/Roo
	      Roo8 = Roo**6
	      Roo14 = Roo8**2
!	      rfLJ = -2.4d01*1.522d-01*(2.0d0*Roo14 - Roo8)*(1.0d0/Roo2)
	      rfLJ = LJprefactor*(2.0d0*Roo14 - Roo8)/Roo2
	      fLJx = fLJx + rfLJ*dx_oo
	      fLJy = fLJy + rfLJ*dy_oo
	      fLJz = fLJz + rfLJ*dz_oo
	   enddo

	elseif (MOD(atompos+1,3).eq.0) then
	   avoid2 = atompos + 1
	   avoid3 = atompos + 2
	   qi = qtip3p_h
	elseif (MOD(atompos-1,3).eq.0) then
	   qi = qtip3p_h
	   avoid2 = atompos-1
	   avoid3 = atompos-2
	else
	   print *, " Stopping here..........#$%^&***"
	   stop
	endif

!----determine coulomb part of TIP3P potential
	do i = 5,natoms
	   if ((i.eq.atompos).or.(i.eq.avoid2).or.(i.eq.avoid3)) cycle
	   if (MOD(i,3).eq.0) then
	      qj = qtip3p_o
	   else
	      qj = qtip3p_h
	   endif
	   dx = psix(atomNo) - psix(atmpl(vbState,i))
	   dy = psiy(atomNo) - psiy(atmpl(vbState,i))
	   dz = psiz(atomNo) - psiz(atmpl(vbState,i))
	   rij = DSQRT(dx**2+dy**2+dz**2)  
!	   rfcoulomb = -1.0d0*qi*qj*fourpieo/(rij**3)
	   rfcoulomb = coulombPrefactor*qi*qj/(rij**3)   
	   fcoulombx = rfcoulomb*dx + fcoulombx 
	   fcoulomby = rfcoulomb*dy + fcoulomby 
	   fcoulombz = rfcoulomb*dz + fcoulombz 
	enddo

	fh2ointerx = 1.0d0*(fLJx + fcoulombx)
	fh2ointery = 1.0d0*(fLJy + fcoulomby)
	fh2ointerz = 1.0d0*(fLJz + fcoulombz)

! Debugging

!	h2o_h2o_forces(vbState,(3*atomNo)-2) = fh2ointerx
!	h2o_h2o_forces(vbState,(3*atomNo)-1) = fh2ointery
!	h2o_h2o_forces(vbState,(3*atomNo)) = fh2ointerz	

	return
	end

!#####################################################################
!####################### SUBROUTINE FH2OINTRA ########################
!#####################################################################

	subroutine fh2ointra(atompos,vbState,fh2ointrax,fh2ointray,fh2ointraz)

	use msevb_common

	implicit none

!-------------------------------------------------------------------------
! 07/08/99 R.A. Christie
!
! routine for determining the intramolecular water forces
!
!-------------------------------------------------------------------------

!----Incoming variables
	integer atompos,vbState,atomNo
	DOUBLE PRECISION FH2OINTRAX,FH2OINTRAY,FH2OINTRAZ
!----Local
	DOUBLE PRECISION DXROH(2),DYROH(2),DZROH(2),ROH(2),RFH2OSTRX,RFH2OSTRY, &
      		rfh2ostrz,rfh2obendx,rfh2obendy,rfh2obendz,theta, &
      		rfh2obendx1,rfh2obendy1,rfh2obendz1,AdotA,AdotB,vectSq, &
      		angleDiff,part1x,part2x,part3x,part1y,part2y,part3y,part1z, &
      		part2z,part3z,part4x,part4y,part4z,theta_temp
	integer i,j,otherb1,otherb2
	character*2 type
	DOUBLE PRECISION CHECKED_ACOS

! Position of the atom in the original input vector

	atomNo = atmpl(vbState,atompos)

!#####################################################################
! Consider whether atom selected was originally in a hydronium
! position, or a water position
!##################################################################### 

	if (MOD(atompos,3).eq.0) then
           type = "ox"
	   otherb1 = atmpl(vbState,atompos+1)
           otherb2 = atmpl(vbState,atompos-1)
	elseif (MOD((atompos+1),3).eq.0) then
           type = "hg"
           otherb1 = atmpl(vbState,atompos+1)
           otherb2 = atmpl(vbState,atompos+2)
	elseif (MOD((atompos-1),3).eq.0) then
           type = "lw"
           otherb1 = atmpl(vbState,atompos-1)
           otherb2 = atmpl(vbState,atompos-2)
	endif
	
!----determine the Roh distances in water
	if (MOD(atompos,3).eq.0) then 
           dxroh(1) = psix(otherb1) - psix(atomNo)
           dyroh(1) = psiy(otherb1) - psiy(atomNo)
           dzroh(1) = psiz(otherb1) - psiz(atomNo)
	   roh(1) = interAtomicR(otherb1,atomNo)

           dxroh(2) = psix(otherb2) - psix(atomNo)
           dyroh(2) = psiy(otherb2) - psiy(atomNo)
           dzroh(2) = psiz(otherb2) - psiz(atomNo)
	   roh(2) = interAtomicR(otherb2,atomNo)
	else
	   dxroh(1) = psix(otherb1) - psix(atomNo)
           dyroh(1) = psiy(otherb1) - psiy(atomNo)
           dzroh(1) = psiz(otherb1) - psiz(atomNo)
	   roh(1) = interAtomicR(otherb1,atomNo)
           dxroh(2) = psix(otherb1) - psix(otherb2)
           dyroh(2) = psiy(otherb1) - psiy(otherb2)
           dzroh(2) = psiz(otherb1) - psiz(otherb2)
	   roh(2) = interAtomicR(otherb1,otherb2)
	endif

!---------------------------------------------------------------------
! Calculate the bending terms:
!
!---------------------------------------------------------------------

        theta_temp = dxroh(1)*dxroh(2)+dyroh(1)*dyroh(2) &
                        + dzroh(1)*dzroh(2)
	theta = checked_acos(theta_temp/(roh(1)*roh(2)))

        angleDiff = theta - thetaeq
        vectSq = dxroh(1)*dxroh(2) + dyroh(1)*dyroh(2) + dzroh(1)*dzroh(2)

!----Initially consider the Oxygen
        if (MOD(atompos,3).eq.0) then 
!---------------------------------------------------------------------
!  X
!---------------------------------------------------------------------
	   part1x = (dxroh(1)*vectSq)/(roh(2)*roh(1)**3) + (dxroh(2)*vectSq)/(roh(1)*roh(2)**3)
	   part2x = dxroh(1)/(roh(1)*roh(2)) + dxroh(2)/(roh(1)*roh(2))
	   part3x = -1.0d0*DSIN(theta)

	   if (DABS(part3x)<numericalZeroLimit) then
	      rfh2obendz = large_force
	   else
	      rfh2obendx = h2oktheta*(theta - thetaeq)*(part1x - part2x)/part3x 
	   endif
!---------------------------------------------------------------------
!  Y
!---------------------------------------------------------------------
            part1y = dyroh(2)*vectSq/(roh(1)*roh(2)**3)
            part2y = -1.0d0*(dyroh(1)+dyroh(2))/(roh(1)*roh(2))
            part3y = dyroh(1)*vectSq/(roh(2)*roh(1)**3)
            part4y = 1.0d0 - vectSq**2/((roh(1)**2)*(roh(2)**2))

	    if (DABS(part4y)<numericalZeroLimit) then
	       rfh2obendy = large_force
	    else
	       rfh2obendy = -1.0d0*(h2oktheta*(part1y + part2y + part3y)*angleDiff)/DSQRT(part4y)
	    endif
!---------------------------------------------------------------------
!  Z
!---------------------------------------------------------------------
            part1z = dzroh(2)*vectSq/(roh(1)*roh(2)**3)
            part2z = -1.0d0*(dzroh(1)+dzroh(2))/(roh(1)*roh(2))
            part3z = dzroh(1)*vectSq/(roh(2)*roh(1)**3)
            part4z = 1.0d0 - vectSq**2/((roh(1)**2)*(roh(2)**2))

	    if (DABS(part4z)<numericalZeroLimit) then
	       rfh2obendz = large_force
	    else
	       rfh2obendz = -1.0d0*(h2oktheta*(part1z + part2z + part3z)*angleDiff)/DSQRT(part4z)
	    endif

!----Now Consider Hydrogen
	 else
!---------------------------------------------------------------------
!  X
!--------------------------------------------------------------------- 
	    part1x = dxroh(1)*vectSq/(roh(2)*roh(1)**3) 
	    part2x = dxroh(2)/(Roh(1)*Roh(2))
	    part3x = 1.0d0 - vectSq**2/((roh(1)**2)*(roh(2)**2))

	    if (DABS(part3x)<numericalZeroLimit) then
	       rfh2obendx = large_force
	    else
	       rfh2obendx = -1.0d0*(h2oktheta*(part1x - part2x)*angleDiff)/DSQRT(part3x)
	    endif
!---------------------------------------------------------------------
!  Y
!---------------------------------------------------------------------
            part1y = dyroh(1)*vectSq/(roh(2)*roh(1)**3)
            part2y = dyroh(2)/(Roh(1)*Roh(2))
            part3y = 1.0d0 - vectSq**2/((roh(1)**2)*(roh(2)**2))
	    
	    if (DABS(part3y)<numericalZeroLimit) then
	       rfh2obendy = large_force
	    else
	       rfh2obendy = -1.0d0*(h2oktheta*(part1y - part2y)*angleDiff)/DSQRT(part3y) 
	    endif
!---------------------------------------------------------------------
!  Z
!---------------------------------------------------------------------
            part1z = dzroh(1)*vectSq/(roh(2)*roh(1)**3)
            part2z = dzroh(2)/(Roh(1)*Roh(2))
            part3z = 1.0d0 - vectSq**2/((roh(1)**2)*(roh(2)**2))

	    if (DABS(part3z)<numericalZeroLimit) then
	       rfh2obendz = large_force
	    else
	       rfh2obendz = -1.0d0*(h2oktheta*(part1z - part2z)*angleDiff)/DSQRT(part3z) 
	    endif
	 endif
!---------------------------------------------------------------------
! Calculate the stretching terms:
!
!---------------------------------------------------------------------

	 rfh2ostrx = -1.0d0*dxroh(1)*(h2okb/roh(1))*(roh(1)-h2oroheq) 
	 rfh2ostry = -1.0d0*dyroh(1)*(h2okb/roh(1))*(roh(1)-h2oroheq)
	 rfh2ostrz = -1.0d0*dzroh(1)*(h2okb/roh(1))*(roh(1)-h2oroheq)
	
	 if (type.eq."ox") then
	    rfh2ostrx = rfh2ostrx - 1.0d0*dxroh(2)*(h2okb/roh(2))*(roh(2)-h2oroheq) 
	    rfh2ostry = rfh2ostry - 1.0d0*dyroh(2)*(h2okb/roh(2))*(roh(2)-h2oroheq)
	    rfh2ostrz = rfh2ostrz - 1.0d0*dzroh(2)*(h2okb/roh(2))*(roh(2)-h2oroheq)
	 endif

	 fh2ointrax = 1.0d0*(rfh2ostrx + rfh2obendx)
	 fh2ointray = 1.0d0*(rfh2ostry + rfh2obendy)
	 fh2ointraz = 1.0d0*(rfh2ostrz + rfh2obendz)

! Debugging

!	h2o_intra_forces(vbState,(3*atomNo)-2) = fh2ointrax
!	h2o_intra_forces(vbState,(3*atomNo)-1) = fh2ointray
!	h2o_intra_forces(vbState,(3*atomNo)) = fh2ointraz

	RETURN 
	END

!#####################################################################
!##################### SUBROUTINE OFFD_ATOM1 #########################
!#####################################################################

! Calculate the off diagonal contributions to the force for an H2O atom (O or H)
! or a hydronium H which is not the exchanging H

	subroutine offd_atom1(atomNo,Roo,Roh,Fxatom1,Fyatom1,Fzatom1,vbState1,vbState2) 

	use commons
	use msevb_common
	
	implicit none

!========================================================================
! 29/7/99 R.A. Christie
!
! Routine for calculating the forces on atom type 1
! *NB* This routine also includes f2_force
! *NB* Roo takes on the distance of the H5O2+ dimer being considered
! INPUT:
!	1. atom number(position)
!
!  -- Update: 21-07-00
! Atom 1 undergoes interactions with the H3O+ atoms, the H2O atoms of a
! water molecule and the exchange-coupling terms.
! Thus the derivatives of the V(H3O+), V(H3O+,H2O) and V(OffDiag) are
! included in this routine.
!
!========================================================================

! Subroutine arguments

        integer, intent(IN) :: atomNo
        DOUBLE PRECISION, INTENT(IN) :: ROO,ROH
	DOUBLE PRECISION, INTENT(OUT) :: FXATOM1,FYATOM1,FZATOM1
        integer, intent(IN) :: vbState1, vbState2

! Local variables

	DOUBLE PRECISION DX,DY,DZ,Q	
	integer i,j,wheech
        DOUBLE PRECISION KROO,TANHROO,ALPHAROO
	DOUBLE PRECISION DHIJ_TEMP
	logical zundelH

!       DOUBLE PRECISION, INTENT(IN) :: DROO(3),DROH(3)
!       DOUBLE PRECISION F,G,RIJ,VEX,DVEXDRIJ
!	DOUBLE PRECISION SECHROO,DGDQ,DROODX,DROODY,DROODZ,DROHDX,DROHDY,DROHDZ,DQDX,DQDY,DQDZ,GAMMAROO

	Fxatom1 = 0.0d0
	Fyatom1 = 0.0d0
	Fzatom1 = 0.0d0
        
 	if (num_eig.gt.2) then
 
!---------------------------------------------------------------------
! Determine initial information to be utilised in the loop
! This is now stored from the energy calculations rather than recalculated
!---------------------------------------------------------------------  

!	   q = 5.0d-1*Roo - Roh
!	   alphaRoo = DEXP(-1.0d0*alpha_bleh*(Roo - rooeq))
!	   kRoo = DEXP(-1.0d0*kexp*(Roo - Doo)**2)
!	   tanhRoo = DTANH(beta*(Roo - big_rooeq))
!	   f = (1.0d0 + P*kRoo)*(5.0d-01*(1.0d0 - tanhRoo) + 1.0d01*(alphaRoo))
!	   g = DEXP(-1.0d0*gamma_msevb*(q)**2)

!	   gammaRoo = DEXP(-1.0d0*gamma_msevb*(q)**2)
!	   sechRoo = 1.0d0 - tanhRoo**2
!	   dgdq = -2.0d0*gammaRoo*gamma_msevb*q
!	   dRoodx = dRoo(1)/Roo
!	   dRoody = dRoo(2)/Roo
!	   dRoodz = dRoo(3)/Roo
!	   dRohdx = dRoh(1)/Roh
!	   dRohdy = dRoh(2)/Roh
!	   dRohdz = dRoh(3)/Roh
!	   dqdx = 5.0d-01*dRoodx - dRohdx
!	   dqdy = 5.0d-01*dRoody - dRohdy
!	   dqdz = 5.0d-01*dRoodz - dRohdz    

! Decide whether this is an H2O atom or a Zundel atom

	   zundelH = .FALSE.

 	   if (MOD(atomNo,3).ne.0) then
 	      do i=1,7
 		 if (zundel_species(vbState1,vbState2,i).eq.atomNo) then
 		    zundelH = .TRUE.
 		    exit
 		 endif
 	      enddo
 	   endif

 	   if (zundelH) then

 	      Interact: do i = 1, natoms
 		 do j = 1, 7
 		    if (atmpl(vbState1,i).eq.zundel_species(vbState1,vbState2,j)) cycle Interact
 		 enddo

 		 wheech = atmpl(vbState1,i)

 		 dx = psix(atomNo) - psix(wheech)
 		 dy = psiy(atomNo) - psiy(wheech)
 		 dz = psiz(atomNo) - psiz(wheech)

! 		 Rij = interAtomicR(atomNo,wheech)
! 		 Vex = each_coulomb(atomNo,wheech)     ! subscripts have to be this way around
! 		 dVexdRij = -1.0d0*Vex/Rij 
!		 dHij_temp = dVexdRij*f*g/Rij

		 dHij_temp = -1.0d0*each_coulomb(atomNo,wheech)*zundel_f(vbState1,vbState2)*zundel_g(vbState1,vbState2)/ &
                             (interAtomicR(atomNo,wheech)**2)

 		 Fxatom1 = Fxatom1 + dHij_temp*dx
 		 Fyatom1 = Fyatom1 + dHij_temp*dy
 		 Fzatom1 = Fzatom1 + dHij_temp*dz

 	      enddo Interact

 	   else  ! H2O atom

 	      do i = 1, 7
 		 wheech = zundel_species(vbState1,vbState2,i)

 		 dx = psix(atomNo) - psix(wheech)
 		 dy = psiy(atomNo) - psiy(wheech)
 		 dz = psiz(atomNo) - psiz(wheech)

!!		 Vex = (qexchi*qexchj*fourpieo)/rij
!!		 dVexdRij = -1.0d0*qexchi*qexchj*fourpieo/(Rij**2)

! 		 Rij = interAtomicR(atomNo,wheech)
! 		 Vex = each_coulomb(wheech,atomNo)
! 		 dVexdRij = -1.0d0*Vex/Rij
!		 dHij_temp = dVexdRij*f*g/Rij

		 dHij_temp = -1.0d0*each_coulomb(wheech,atomNo)*zundel_f(vbState1,vbState2)*zundel_g(vbState1,vbState2)/ &
                             (interAtomicR(atomNo,wheech)**2)

 		 Fxatom1 = Fxatom1 + dHij_temp*dx
 		 Fyatom1 = Fyatom1 + dHij_temp*dy
 		 Fzatom1 = Fzatom1 + dHij_temp*dz
 	      enddo
	   endif 
	endif

	RETURN
	END

!##############################################################################
!####################### SUBROUTINE OFFD_ATOM3 ################################
!##############################################################################

! Called when the atom in question is the hydronium O atom in one of the two
! interacting VB states

	subroutine offd_atom3(atomNo,Roo,Roh,dRoo,dRoh,Fxatom3,Fyatom3,Fzatom3, &
      		vbState1,vbState2)

	use commons
	use msevb_common

	implicit none

! atomNo is the position of the O atom in the initial coordinate vector

!===============================================================
! 23-7-00 R.A. Christie
! Routine to determine the atom-type-3 forces. 
!
! Need to determine the force on the particle in question 
! resulting from both the coulombic interaction, the asymmetric
! stretch and the Roo repulsion
!
!===============================================================

!----Variables to be passed
	DOUBLE PRECISION ROO,ROH,DROO(3),DROH(3),FXATOM3,FYATOM3,FZATOM3

	integer atomNo,vbState1,vbState2
!----Local Variables
	DOUBLE PRECISION DX,DY,DZ,RIJ,RCOULOMB,COULOMB,COULOMBSQ,QEXCHJ, &
      		part1f,part2f,part3fai,part3fa,coeff,derf2f, &
      		alphaRoo,f2_f,part3f,derf2_f,f2force,f2Roo,qexchi, &
      		rcoulombsq
	integer i,start,end,wheech,j,ninteractions,count
!----New method variables
	DOUBLE PRECISION KROO,GAMMAROO,TANHROO,SECHROO,DXOO,DYOO,DZOO,DXOH,DYOH, &
      		dzoh,q
!----Even more recent method variables
	DOUBLE PRECISION F,DFDROO,G,DGDQ,DROODX,DROODY,DROODZ,DROHDX,DROHDY,DROHDZ, &
      	   dqdx,dqdy,dqdz,dHijdx,dHijdy,dHijdz, &
           dVexdRij,Vex,f1,f2,df2dRoo,dVexdRijAdivRij
	logical inexchange

!	print *, 'Atom:', position, 'x-coord:', 

!----Initialize variables

	Fxatom3 = 0.0d0
	Fyatom3 = 0.0d0
	Fzatom3 = 0.0d0

	dHijdx = 0.0d0
	dHijdy = 0.0d0
	dHijdz = 0.0d0

!----determine the initial information 
	q = 5.0d-1*Roo - Roh
	alphaRoo = DEXP(-1.0d0*alpha_bleh*(Roo - rooeq))
	kRoo = DEXP(-1.0d0*kexp*(Roo - Doo)*(Roo - Doo))
	gammaRoo = DEXP(-1.0d0*gamma_msevb*q*q)
	tanhRoo = DTANH(beta*(Roo - big_rooeq))
	sechRoo = 1.0d0 - tanhRoo**2
	qexchi = qexcho
	dxoo = dRoo(1)
	dyoo = dRoo(2)
	dzoo = dRoo(3)
	dxoh = dRoh(1)
	dyoh = dRoh(2)
	dzoh = dRoh(3)

!#####################################################################
! New method variables:
!
 	f = (1.0d0 + P*kRoo)*(5.0d-01*(1.0d0 - tanhRoo) + 1.0d01*(alphaRoo))

	dfdroo = (1.0d0+ kRoo*P)*(-1.0d01*alpha_bleh*alphaRoo - 5.0d-01*beta*sechRoo) - &
      		2.0d0*kRoo*kexp*P*(Roo - Doo)*(1.0d01*alphaRoo + 5.0d-01*(1.0d0 - tanhroo))

!
!----AFTER REVIEW OF DERIVATIVES - 19/08/02 - atom3_offd_forces.nb:
!----THIS NOW AGREES WITH MATHEMATICA VERSION OF df/dRoo.
! trj25 22/10/03 - not according to Mathematica4 it doesn`t, previous expression was correct 
!	dfdroo = kRoo*P*(-1.0d01*alpha_bleh*alphaRoo - 5.0d-01*beta*sechRoo) -
!     &         2.0d0*kRoo*kexp*P*(Roo - Doo)*(1.0d01*alphaRoo + 5.0d-01*(1.0d0 - tanhroo)) 

!	f1 = (1.0d0 + P*kRoo)
!	f2 = 5.0d-01*(1.0d0 - tanhRoo) + 1.0d01*(alphaRoo)
!	df2droo = (-1.0d01*alpha_bleh*alphaRoo - 5.0d-01*beta*sechRoo)
	g = gammaRoo
	dgdq = -2.0d0*gammaRoo*gamma_msevb*q
	dRoodx = dxoo/Roo
	dRoody = dyoo/Roo 
	dRoodz = dzoo/Roo
	dRohdx = dxoh/Roh
	dRohdy = dyoh/Roh
	dRohdz = dzoh/Roh

! trj25 Change so that dq/dcoordinate now includes the influence of d(Roh)/dcoordinate
! in a manner consistent with the energy calculation

	if (atomNo.eq.zundel_species(vbState1,vbState2,3)) then
	   dqdx = 0.5d0*dRoodx - dRohdx 
	   dqdy = 0.5d0*dRoody - dRohdy
	   dqdz = 0.5d0*dRoodz - dRohdz
	else
	   dqdx = 0.5d0*dRoodx
	   dqdy = 0.5d0*dRoody
	   dqdz = 0.5d0*dRoodz
	endif

!---------------------------------------------------------------------
! Determine the Coulombic interaction arising from the Exchange term
! first. 
! Note that this is calculated as the sum of all interactions to be 
! used later in the full for this off-diagonal term
!
!---------------------------------------------------------------------

 	if (num_eig.gt.2) then

	   count = 0

	   Exchange: do i = 1,natoms
	      if (count.eq.(natoms-7)) exit

! Check that we are not calculating a hydronium-hydronium exchange interaction
	      do j = 1,7
		 if (atmpl(vbState2,i).eq.zundel_species(vbState1,vbState2,j)) cycle Exchange
	      enddo

	      wheech = atmpl(vbState2,i)
	      qexchi = qexcho
	      count = count + 1
 
	      if (MOD(wheech,3).eq.0) then
		 qexchj = qtip3p_o
	      else
		 qexchj = qtip3p_h
	      endif
	   
	      dx = psix(atomNo) - psix(wheech)
	      dy = psiy(atomNo) - psiy(wheech)
	      dz = psiz(atomNo) - psiz(wheech)
	      Rij = interAtomicR(atomNo,wheech)
	      Vex = (qexchi*qexchj*fourpieo)/Rij
	      dVexdRij = -1.0d0*Vex/Rij

	      dVexdRijAdivRij = dVexdRij*f*g/Rij

!	   if (inexchange) then
!	      dHijdx_temp = (gammaRoo*f1*(vij+Vex)*df2dRoo*(dxoo/Roo)) - 
!     &  (fourpieo*gammaRoo*f1*qexchi*qexchj*dx*f2)/(Rij)**3 - 
!     &  (2.0d0*gammaRoo*kRoo*kexp*P*dxoo*(vij+Vex)*(Roo-Doo)*f2)/Roo -
!     &  (2.0d0*gamma_msevb*f1*(vij+Vex)*(dxoo/(2.0d0*Roo)-dxoh/Roh)*(Roo/2.0d0 - Roh)*f2)
!	      dHijdy_temp = (gammaRoo*f1*(vij+Vex)*df2dRoo*(dyoo/Roo)) -
!     &  (fourpieo*gammaRoo*f1*qexchi*qexchj*dy*f2)/(Rij)**3 -
!     &  (2.0d0*gammaRoo*kRoo*kexp*P*dyoo*(vij+Vex)*(Roo-Doo)*f2)/Roo -
!     &  (2.0d0*gamma_msevb*f1*(vij+Vex)*(dyoo/(2.0d0*Roo)-dyoh/Roh)*(Roo/2.0d0 - Roh)*f2) 
!	      dHijdz_temp = (gammaRoo*f1*(vij+Vex)*df2dRoo*(dzoo/Roo)) -
!     &  (fourpieo*gammaRoo*f1*qexchi*qexchj*dz*f2)/(Rij)**3 -
!     &  (2.0d0*gammaRoo*kRoo*kexp*P*dzoo*(vij+Vex)*(Roo-Doo)*f2)/Roo -
!     &  (2.0d0*gamma_msevb*f1*(vij+Vex)*(dzoo/(2.0d0*Roo)-dzoh/Roh)*(Roo/2.0d0 - Roh)*f2) 
!	print *, "dHijdx_temp = ", dHijdx_temp

	      dHijdx = dHijdx + dVexdRijAdivRij*dx
	      dHijdy = dHijdy + dVexdRijAdivRij*dy
	      dHijdz = dHijdz + dVexdRijAdivRij*dz

!----End of Loop over the solvating atoms

	   enddo Exchange

!----END OF IF STATEMENT FOR > 2 EIGENSTATES

! Add the second part of the derivatives, use the total exchange energy between these two VB states	

	   dHijdx = dHijdx + (vij+vijexch(vbState1,vbState2))*(dfdroo*droodx*g + f*dgdq*dqdx)
	   dHijdy = dHijdy + (vij+vijexch(vbState1,vbState2))*(dfdroo*droody*g + f*dgdq*dqdy)
	   dHijdz = dHijdz + (vij+vijexch(vbState1,vbState2))*(dfdroo*droodz*g + f*dgdq*dqdz)

	else
!----DETERMINE THE OFF-DIAGONAL ELEMENT FOR 2 EIGENSTATE SYSTEM

	    dHijdx = vij*(dfdroo*g*droodx + f*dgdq*dqdx)
	    dHijdy = vij*(dfdroo*g*droody + f*dgdq*dqdy)  
	    dHijdz = vij*(dfdroo*g*droodz + f*dgdq*dqdz)  
	endif

!	    print *, 'At end:'
!	print *, 'vij:', vij
!	print *, 'Roo:', roo
!	print *, 'Roh:', roh
!	print *, 'f:', f
!	print *, 'df/droo:', dfdroo
!	print *, 'g:', g
!	print *, 'dg/dq:', dgdq
!	print *, 'q:', q
!	print *, 'dq/dx:', dqdx
!	print *, 'dq/dy:', dqdy
!	print *, 'dq/dz:', dqdz
!	print *, 'droo/dx:', droodx
!	print *, 'droo/dy:', droody	
!	print *, 'droo/dz:', droodz
!	print *, 'Total exchange energy: ', vijexch(eigenstate1,eigenstate2)

!	print *, 'dHij/dx:', dHijdx
!	print *, 'dHij/dy:', dHijdy
!	print *, 'dHij/dz:', dHijdz

	Fxatom3 = dHijdx
	Fyatom3 = dHijdy
	Fzatom3 = dHijdz

	RETURN
	END

!#####################################################################
!#################### SUBROUTINE OFFD_ATOM4 ##########################
!#####################################################################

	subroutine offd_atom4(atomNo,Roo,Roh,dRoo,dRoh,Fxatom4,Fyatom4,Fzatom4,vbState1,vbState2) 

	use commons
	use msevb_common
	  
	implicit none

! Calculate the off-diagonal forces for the exchanging H atom

!===============================================================
! 30/7/99 R.A. Christie
! Routine to calculate the atom-4-type forces. Applies only to
! atom4
!
! *NB* function f2force is in routine f_atom1
! *NB* Need two Roo distances as input
!
!===============================================================

!----Incoming variables
	DOUBLE PRECISION ROO,ROH,DROO(3),DROH(3),FXATOM4,FYATOM4,FZATOM4
        integer atomNo,vbState1,vbState2
!----Local variables
	DOUBLE PRECISION DX,DY,DZ,REXCH,COULOMB,TMPFX,TMPFY,TMPFZ,QEXCHJ, &
      		qexchi,part1f,part2f,alphaRoo,q,Rij,rcoulomb
	integer i,j,start,wheech,count
!----New method variables
        DOUBLE PRECISION KROO,GAMMAROO,TANHROO,SECHROO,DXOO,DYOO,DZOO,DXOH,DYOH, &
                dzoh
!----New method variables
        DOUBLE PRECISION F,G,DGDQ,DROHDX,DROHDY,DROHDZ,DQDX,DQDY,DQDZ,DHIJDX,DHIJDY, &
      		dHijdz,Vex,dVexdRij
	logical inexchange

!---------------------------------------------------------------------
! Initialize variables
!
!---------------------------------------------------------------------  

	Fxatom4 = 0.0d0
	Fyatom4 = 0.0d0
	Fzatom4 = 0.0d0

        dHijdx = 0.0d0
        dHijdy = 0.0d0
        dHijdz = 0.0d0  
	count = 0

!---------------------------------------------------------------------
! Determine the initial information
!
!---------------------------------------------------------------------
	q = Roo*5.0d-01 - Roh
	alphaRoo = DEXP(-1.0d0*alpha_bleh*(Roo - rooeq))
        kRoo = DEXP(-1.0d0*kexp*(Roo - Doo)**2)
        gammaRoo = DEXP(-1.0d0*gamma_msevb*(q)**2)
        tanhRoo = DTANH(beta*(Roo - big_rooeq))
        sechRoo = 1.0d0 - tanhRoo**2
        qexchi = qexchh
        dxoo = dRoo(1)
        dyoo = dRoo(2)
        dzoo = dRoo(3)
        dxoh = dRoh(1)
        dyoh = dRoh(2)
        dzoh = dRoh(3) 

!#####################################################################
! New method variables
!
	f = (1.0d0 + P*kRoo)*(5.0d-01*(1.0d0 - tanhRoo) + 1.0d01*(alphaRoo))
        g = gammaRoo
        dgdq = -2.0d0*gammaRoo*gamma_msevb*q
        dRohdx = dxoh/Roh
        dRohdy = dyoh/Roh
        dRohdz = dzoh/Roh
        dqdx = -1.0d0*dRohdx
        dqdy = -1.0d0*dRohdy
        dqdz = -1.0d0*dRohdz    

!################################################################################
! Calculate the Off-Diagonal Element.
! Two cases:
!   (1) Any system other than H5O2+
!   (2) H5O2+
!
! Loop over all atoms in this eigenstate, given by atmpl(eigenstate,j), which 
! are not part of the H5O2+ in this eigenstate, if the atom in question, atompos,
! is part of the H5O2+ in this eigenstate. If atompos isn't part of the H5O2+
! in this eigenstate, then the loop over atoms should only include the H5O2+
! atoms for this eigenstate.
!
!################################################################################ 

 	if (num_eig.gt.2) then

!----Determine which atoms atomNo should be interacting with

	   Exchange: do i = 1,natoms
	      if (count.eq.(natoms-7)) exit

!----Determine whether atom i is in H5O2+

	      do j = 1,7
		 if (atmpl(vbState1,i).eq.zundel_species(vbState1,vbState2,j)) cycle Exchange
	      enddo 
	      wheech = atmpl(vbState1,i) 
	      count = count + 1

	      if (MOD(wheech,3).eq.0) then 
		 qexchj = qtip3p_o
	      else 
		 qexchj = qtip3p_h
	      endif

	      dx = psix(atomNo) - psix(wheech)
	      dy = psiy(atomNo) - psiy(wheech)
	      dz = psiz(atomNo) - psiz(wheech)

	      Rij = interAtomicR(atomNo,wheech)
	      Vex = (qexchi*qexchj*fourpieo)/Rij 
	      dVexdRij = -1.0d0*qexchi*qexchj*fourpieo/(Rij**2)

	      dHijdx = dHijdx + dVexdRij*(dx/Rij)
	      dHijdy = dHijdy + dVexdRij*(dy/Rij)
	      dHijdz = dHijdz + dVexdRij*(dz/Rij)	      

	   enddo Exchange

! Complete the exchange terms

	   dHijdx = f*g*dHijdx + (vij+vijexch(vbState1,vbState2))*dgdq*f*dqdx 
	   dHijdy = f*g*dHijdy + (vij+vijexch(vbState1,vbState2))*dgdq*f*dqdy
	   dHijdz = f*g*dHijdz + (vij+vijexch(vbState1,vbState2))*dgdq*f*dqdz	

!----Consider now the case when have only H5O2+ i.e. no exchange term
	else

	   dHijdx = vij*f*dgdq*dqdx  
	   dHijdy = vij*f*dgdq*dqdy 
	   dHijdz = vij*f*dgdq*dqdz 

	endif

	Fxatom4 = dHijdx
	Fyatom4 = dHijdy
	Fzatom4 = dHijdz

	RETURN
	END

















