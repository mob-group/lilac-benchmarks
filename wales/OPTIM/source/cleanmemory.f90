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
! Cleans up all msevb allocated memory at the end of a run

SUBROUTINE cleanMemory ()

  USE msevb_common

  IMPLICIT NONE

! MSEVB storage

  DEALLOCATE(each_coulomb, water_inter_coulomb, ljr)
  DEALLOCATE(inter_coulomb, lj_inter, repulse_inter)
  DEALLOCATE(atom_coulomb)

  DEALLOCATE(atmpl, statesInteract)
  DEALLOCATE(zundel_species)
  DEALLOCATE(zundel_f,zundel_g)

  DEALLOCATE(psix, psiy, psiz)
  DEALLOCATE(interAtomicR)

  RETURN

END SUBROUTINE cleanMemory
