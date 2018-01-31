
SUBROUTINE NIMETP (XC, DX, ENERGY, GTEST, STEST)

! *****************************************************************************************
! This is the EAM 4 Ni(100)/H for the interaction of a single H atom with a Ni(100)
! surface by Yury Suleimanov and added by Judith Rommel, if used cite:

!S. M. Foiles, M. I. Baskes and M. S. Daw, Phys. Rev. B 33, 7983 (1986).
!M. S. Daw and M. I. Baskes, Phys. Rev. B 29, 6443 (1984).
!T. N. Truong, D. G. Truhlar and B. C. Garrett, J. Phys. Chem. 93, 8227 (1989).
!T. N. Truong and D. G. Truhlar, J. Phys. Chem. 94, 8262 (1990).
!REV B 51, 9985 (1995)
!and Yure's RPMD paper:
!Suleimanov, Y. V., J. Phys. Chem. C 116.20 (2012): 11141-11153.

! Input/output geometry in Angstrom, energies in eV, and gradients in eV/Angstrom!
! This version is for 1 hydrogen atom only!

      USE MODHESS
      USE COMMONS, ONLY: NATOMS, ZSYM
      USE KEY, ONLY: FROZEN, NFREEZE

      IMPLICIT NONE
      DOUBLE PRECISION :: XC(3*NATOMS), DX(3*NATOMS) !Cartesian coords, derivatives
      DOUBLE PRECISION :: ENERGY
!      CHARACTER(LEN=5)::ZSYM(393)
      integer          :: nclassic,nrigid,nb,i1
      integer :: iflag ! -1 - coordinates to Angs, +1 - coordinates, energy and forces to a.u
      double precision :: boxlx,boxly,boxlz
      double precision :: qh(3,1,1),dvdqh(3,1,1)
      double precision :: qmclassic(3,NATOMS-NFREEZE-1,1),dvdqmclassic(3,NATOMS-NFREEZE-1,1)
      double precision :: qmrigid(3,NFREEZE,1),dvdqmrigid(3,NFREEZE,1)
      character(len=6)      :: boundary
      logical :: h_included,GTEST, STEST
!Description of the variables
!nclassic - number of classical Ni atoms
!nrigid - number of rigid Ni atoms
!nb - number of the ring polymer beads
!boxlx - length of the simulation box in x direction
!boxly - length of the simulation box in y direction
!boxlz - length of the simulation box in z direction
!boundary - 'afixed', 'bfixed', or 'cfixed'. 'cfixed' or 'bfixed'- the lattices with 2D periodic boundary conditions
!h_included - '.true.' or  '.false.'  for hydrogen (whether it is on the surface or not)
!qmclassic  - array for classic Ni atoms coordinates
!qmrigid  - array for rigid Ni atoms coordinates
!qh - array for H atom coordinates
!v - interaction potential
!dvdqmclassic - array for classic Ni atoms forces
!dvdqmrigid  - array for rigid Ni atoms forces
!dvdqh  - array for H atom forces
! *****************************************************************************************

! PRINT*, 'Nimet potential'

nclassic=NATOMS-NFREEZE-1 !natoms-nfrozen !number of mobile atoms not including H
nrigid=NFREEZE!nfrozen ! number of frozen atoms
nb=1

!boxlength (caluclated based on approx 2 times the atom radius 4.7 of NI in a.u.)
boxlx = 24.64d0 ! 65.8d0
boxly = 24.64d0 !65.8d0
boxlz = 7.04d0 !14.1d0

!boundary - 'afixed', 'bfixed', or 'cfixed'. 'cfixed' or 'bfixed'- the lattices with 2D periodic boundary conditions
boundary='bfixed'
!h_included - '.true.' or  '.false.'  for hydrogen (whether it is on the surface or not)
h_included =.false.


!initialise coords and forces
qh(:,:,:) =0.d0
dvdqh(:,:,:) =0.d0
qmclassic(:,:,:)=0.D0
dvdqmclassic(:,:,:)=0.D0
qmrigid(:,:,:)=0.D0
dvdqmrigid(:,:,:)=0.D0
ENERGY = 0.0D0
DX(:)=0

! The geometry of the mobile classical atoms
DO  I1 = 1,nclassic+1
   IF ((ZSYM(I1).EQ.'NI').or.(ZSYM(I1).EQ.'Ni')) THEN
      qmclassic(1,I1-1,1)= XC(3*I1-2)!I1-1 because of H excluded
      qmclassic(2,I1-1,1)= XC(3*I1-1)
      qmclassic(3,I1-1,1)= XC(3*I1)
   ELSE IF (ZSYM(I1).EQ.'H') THEN
      qh(1,I1,1)= XC(3*I1-2)
      qh(2,I1,1)= XC(3*I1-1)
      qh(3,I1,1)= XC(3*I1)
      h_included =.true.
   END IF
ENDDO


! The geometry of the rigid atoms
IF (nclassic+2 .LE. NATOMS) THEN
   DO  I1 = nclassic+2,NATOMS
      IF (ZSYM(I1).EQ.'NI') THEN
         qmrigid(1,I1-nclassic-1,1)= XC(3*I1-2)!I1-nclassic-1 because of H excluded
         qmrigid(2,I1-nclassic-1,1)= XC(3*I1-1)
         qmrigid(3,I1-nclassic-1,1)= XC(3*I1)
      ELSE
         PRINT*, 'There are gas atoms in frozen zone'
         STOP
      END IF
   ENDDO
ENDIF


IF (GTEST) THEN
! PRINT*, 'Calculating energy and forces.  One moment please ...'

    !     Calculate the new forces and the interaction energy
    call eamnih_forces_umb(nclassic,nrigid,nb,boxlx,boxly,boxlz, &
            boundary,h_included,qmclassic,qmrigid,&
            qh,ENERGY,dvdqmclassic,dvdqmrigid,dvdqh)

    !print*,dvdqh,'dvdqh'
!    print*,dvdqmclassic,'dvdqmclassic'
!    print*,dvdqmclassic(1,1,1),dvdqmclassic(2,1,1),dvdqmclassic(3,1,1),'dvdqmclassic'

  !Gradient of Hydrogen
  DO I1=1,nb
     DX(3*I1-2)=dvdqh(1,I1,1)
     DX(3*I1-1)=dvdqh(2,I1,1)
     DX(3*I1)=dvdqh(3,I1,1)
  ENDDO
  !Gradient of mobile metal atoms
  DO I1 = 1,nclassic
     DX(3*(I1+nb) - 2) = dvdqmclassic(1,I1,1)
     DX(3*(I1+nb) - 1) = dvdqmclassic(2,I1,1)
     DX(3*(I1+nb))     = dvdqmclassic(3,I1,1)
  ENDDO

!    PRINT*, 'Gradient:'
!            DO I1 = 1,15,3!1,(nclassic+1)*3,3
!                PRINT*, DX(I1), DX(I1+1), DX(I1+2)
!            ENDDO

!           PRINT*, 'Energy = ', ENERGY, 'eV'
!           PRINT*, 'Done with ENE'

ELSE
!PRINT*, 'Calculating energy.  One moment please ...'

call eamnih_energy(nclassic,nrigid,nb,boxlx,boxly,boxlz, &
        boundary,h_included,qmclassic,qmrigid,qh,ENERGY)

!            PRINT*, 'Energy = ', ENERGY,'eV'
!            PRINT*, 'Done with ENE'

ENDIF

IF (STEST) THEN
    PRINT*, 'Warning: There is no analytical hessian in the nimet (nickel + h) potential.'
    STOP
    !CALL MAKENUMHESS(XC,NATOMS)
ENDIF

!call writespec_xyz(17,XC)

END SUBROUTINE NIMETP 

subroutine writespec_xyz(out_unit,COORDS)
      USE COMMONS, ONLY: NATOMS, ZSYM
  implicit none
  DOUBLE PRECISION, intent(in):: coords(3*NATOMS)
  INTEGER :: IAT
  INTEGER, INTENT(IN):: out_unit

! **********************************************************************
!  open (out_unit,file="Nimetpath.xyz",ACTION="WRITE")

  write(out_unit,*) NATOMS
  write(out_unit,*)
    DO  IAT = 1,NATOMS
          write(out_unit,'(a2,3f12.7)') ZSYM(IAT),coords(3*IAT-2),&
                coords(3*IAT-1),coords(3*IAT)
    ENDDO
!  call flush(out_unit)
!  close(out_unit)
end subroutine writespec_xyz



subroutine eamnih_forces_umb(nclassic,nrigid,nb,boxlx,boxly,boxlz, &
        boundary,h_included,qmclassic,qmrigid,qh,v, &
        dvdqmclassic,dvdqmrigid,dvdqh)
!
       Implicit None
!
!      EAM 4 Ni(100)/H. Input geometry in angstom Output energy, grad in eV
!
       integer :: i,j,k,ilist
       integer, intent(in) :: nclassic,nrigid, nb
       double precision :: dx,dy,dz
       double precision :: vmm,vmh,fm,fh !interaction potene mm, mh, embedding ene fm, fh
       double precision, intent(in) :: boxlx,boxly,boxlz
       double precision :: onboxlx,onboxly,onboxlz
       character(len=6),intent(in) :: boundary
       double precision :: rc1sq,rc1,drsq,dv,dvdr,dvdx,dvdy,dvdz,dr
       double precision :: pma,pha,ph,pm,df,dxr,dyr,dzr,dpadr_m,dpadr_h
       double precision :: dfidp,dfjdp,fjx,fjy,fjz,summ,fix,fiy,fiz
       double precision,intent(out) :: v
       logical, intent(in) :: h_included
!
!
       double precision :: f_mclh(1:3,1:1,nclassic,1:nb), &
        f_hmcl(1:3,1:1,nclassic,1:nb), &
        f_mrigh(1:3,1:1,nrigid,1:nb), &
        f_hmrig(1:3,1:1,nrigid,1:nb), &
        f_mclcl(1:3,nclassic,nclassic,1:1), &
        f_mclrig(1:3,nrigid,nclassic,1:1), &
        f_mrigcl(1:3,nrigid,nclassic,1:1)
!
       INTEGER :: nlist_mclh(nclassic,1:nb), &
        list_mclh(1,nclassic,1:nb),      &
        nlist_mrigh(nrigid,1:nb), &
        list_mrigh(1,nrigid,1:nb), &
        nlist_mclcl(nclassic,1:1), &
        list_mclcl(nclassic,nclassic,1:1), &
        nlist_mclrig(nclassic,1:1), &
        list_mclrig(nrigid,nclassic,1:1), &
        nlist_mrigrig(nrigid,1:1)
!
       double precision :: dfdp_mclassic(nclassic,1:nb), &
        dfdp_mrigid(nrigid,1:nb), &
        dfdp_h(1:1,1:nb)
!
       double precision,intent(in) :: qmclassic(1:3,nclassic,1:1), &
         qmrigid(1:3,nrigid,1:1)
!
       double precision,intent(in) :: qh(1:3,1:1,1:nb)
!
       double precision,intent(out) :: dvdqmclassic(1:3,nclassic,1:1), &
        dvdqmrigid(1:3,nrigid,1:1), &
        dvdqh(1:3,1:1,1:nb)
!
       double precision :: eden_mrigid(nrigid,1:nb), &
                          eden_mclassic(nclassic,1:nb), &
                          eden_h(1:1,1:nb)
!
        double precision :: smooth,dsmooth

!***************************************************************

!        nh = 1 ! this version is for 1 hydrogen atom only!
        vmm = 0.d0
        vmh = 0.d0
        v = 0.d0
        fm = 0.d0
        fh = 0.d0
        rc1 = 10.D0 ! 10.0d0  !!!  DEBUG ! fixed cut off for Ni-Ni and H-Ni interactions (Angstroms)
        rc1sq = rc1*rc1
!
!
!
!
        dvdqmclassic = 0.d0
        dvdqmrigid = 0.d0
        dvdqh = 0.d0
        dv = 0.d0
        dvdr = 0.d0

!      Cl-H
       nlist_mclh = 0
       list_mclh = 0
       f_mclh = 0.d0
       f_hmcl = 0.d0
!
!      RIG-H
       nlist_mrigh = 0
       list_mrigh = 0
       f_mrigh = 0.d0
       f_hmrig = 0.d0
!
!      Cl-Cl
       nlist_mclcl = 0
       list_mclcl = 0
       f_mclcl = 0.d0
 !
 !     Cl-RIG
       nlist_mclrig = 0
       list_mclrig = 0
       f_mclrig = 0.d0
       f_mrigcl = 0.d0
!
!      RIG-RIG
       nlist_mrigrig = 0
!       list_mrigrig = 0
!       f_mrigrig = 0.d0
!
!
       dfdp_mclassic = 0.d0
       dfdp_mrigid = 0.d0
       dfdp_h = 0.d0
!
       eden_mclassic = 0.d0
       eden_mrigid = 0.d0
       eden_h = 0.d0
!
!      ! Converting in eV/ANgs if input is not in Angstrom
!       call atomic_units(nclassic,nrigid,nb,boxlx,boxly,boxlz,qmclassic,&
!                    qmrigid,qh,v,dvdqmclassic,dvdqmrigid,dvdqh, -1)
        onboxlx = 1.d0/boxlx
        onboxly = 1.d0/boxly
        onboxlz = 1.d0/boxlz
!
    IF(nclassic.ne.0) then
        do i = 1,nb
        do j = 1, nclassic
!
              dx = qmclassic(1,j,1)-qh(1,1,i)
              dy = qmclassic(2,j,1)-qh(2,1,i)
              dz = qmclassic(3,j,1)-qh(3,1,i)
          if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
            dx = dx-boxlx*nint(onboxlx*dx)
            dy = dy-boxly*nint(onboxly*dy)
!            !dz = dz-boxlz*nint(onboxlz*dz)
           end if
              drsq = dx*dx + dy*dy + dz*dz
              if (drsq.lt.rc1sq) then
                 ! Add to neighbour list
!
                 nlist_mclh(j,i) = nlist_mclh(j,i)+1
                 list_mclh(nlist_mclh(j,i),j,i) = 1
!
                 dr = dsqrt(drsq)
                 CALL SMOOTHING(dR,Smooth,DSmooth)
                 call pair_NiH(dr,dv,dvdr,Smooth,DSmooth)
                 vmh = vmh + dv
                 dxr = dx/dr
                 dyr = dy/dr
                 dzr = dz/dr
                 dvdx = dvdr*dxr
                 dvdy = dvdr*dyr
                 dvdz = dvdr*dzr
!
                 dvdqmclassic(1,j,1) = dvdqmclassic(1,j,1) + dvdx
                 dvdqmclassic(2,j,1) = dvdqmclassic(2,j,1) + dvdy
                 dvdqmclassic(3,j,1) = dvdqmclassic(3,j,1) + dvdz
                 dvdqh(1,1,i) = dvdqh(1,1,i)  - dvdx
                 dvdqh(2,1,i) = dvdqh(2,1,i)  - dvdy
                 dvdqh(3,1,i) = dvdqh(3,1,i)  - dvdz
!
!
                 call rho_hnew(dr,pha,dpadr_h,Smooth,DSmooth)
                 call rho_Ni(dr,pma,dpadr_m,Smooth,DSmooth)
!
                 eden_mclassic(j,i) = eden_mclassic(j,i) + pha
                 eden_h(1,i) = eden_h(1,i) + pma
!
!       !       DERIVATIVES HERE
!
                 ilist = nlist_mclh(j,i)
                 f_mclh(1,ilist,j,i) = dpadr_m*dxr
                 f_mclh(2,ilist,j,i) = dpadr_m*dyr
                 f_mclh(3,ilist,j,i) = dpadr_m*dzr
                 f_hmcl(1,ilist,j,i) = dpadr_h*dxr
                 f_hmcl(2,ilist,j,i) = dpadr_h*dyr
                 f_hmcl(3,ilist,j,i) = dpadr_h*dzr
!
!
!
        end if
        end do
        end do
    END IF

    IF(nrigid.ne.0) then
        do i = 1,nb
        do j = 1, nrigid
!
              dx = qmrigid(1,j,1)-qh(1,1,i)
              dy = qmrigid(2,j,1)-qh(2,1,i)
              dz = qmrigid(3,j,1)-qh(3,1,i)
          if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
            dx = dx-boxlx*nint(onboxlx*dx)
            dy = dy-boxly*nint(onboxly*dy)
!            !dz = dz-boxlz*nint(onboxlz*dz)
          end if
              drsq = dx*dx + dy*dy + dz*dz
              if (drsq.lt.rc1sq) then
                 ! Add to neighbour list
!
                 nlist_mrigh(j,i) = nlist_mrigh(j,i)+1
                 list_mrigh(nlist_mrigh(j,i),j,i) = 1
!
                 dr = dsqrt(drsq)
                 CALL SMOOTHING(dR,Smooth,DSmooth)
                 call pair_NiH(dr,dv,dvdr,Smooth,DSmooth)
                 vmh = vmh + dv
                 dxr = dx/dr
                 dyr = dy/dr
                 dzr = dz/dr
                 dvdx = dvdr*dxr
                 dvdy = dvdr*dyr
                 dvdz = dvdr*dzr
!
                 dvdqh(1,1,i) = dvdqh(1,1,i)  - dvdx
                 dvdqh(2,1,i) = dvdqh(2,1,i)  - dvdy
                 dvdqh(3,1,i) = dvdqh(3,1,i)  - dvdz
!
!
!
                 call rho_hnew(dr,pha,dpadr_h,Smooth,DSmooth)
                 call rho_Ni(dr,pma,dpadr_m,Smooth,DSmooth)
!
                 eden_mrigid(j,i) = eden_mrigid(j,i) + pha
                 eden_h(1,i) = eden_h(1,i) + pma
!
!       !       DERIVATIVES HERE
!
                 ilist = nlist_mrigh(j,i)
                 f_mrigh(1,ilist,j,i) = dpadr_m*dxr
                 f_mrigh(2,ilist,j,i) = dpadr_m*dyr
                 f_mrigh(3,ilist,j,i) = dpadr_m*dzr
                 f_hmrig(1,ilist,j,i) = dpadr_h*dxr
                 f_hmrig(2,ilist,j,i) = dpadr_h*dyr
                 f_hmrig(3,ilist,j,i) = dpadr_h*dzr
!
!
!
        end if
        end do
        end do
    END IF

      ! Classic-Classic

        do i=1,nclassic
        do j=i+1,nclassic
              dx = qmclassic(1,i,1) - qmclassic(1,j,1)
              dy = qmclassic(2,i,1) - qmclassic(2,j,1)
              dz = qmclassic(3,i,1) - qmclassic(3,j,1)
          if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
            dx = dx-boxlx*nint(onboxlx*dx)
            dy = dy-boxly*nint(onboxly*dy)
!            !dz = dz-boxlz*nint(onboxlz*dz)
          end if
              drsq = dx*dx + dy*dy + dz*dz
              if (drsq.lt.rc1sq) then
                 ! Add to neighbour list
!
                 nlist_mclcl(i,1) = nlist_mclcl(i,1)+1
                 list_mclcl(nlist_mclcl(i,1),i,1) = j
!
                 dr = dsqrt(drsq)
!
                 CALL SMOOTHING(dR,Smooth,DSmooth)
                call pair_nini(dr,dv,dvdr,Smooth,DSmooth)
                 vmm = vmm  + dv
                 dxr = dx/dr
                 dyr = dy/dr
                 dzr = dz/dr
                 dvdx = dvdr*dxr
                 dvdy = dvdr*dyr
                 dvdz = dvdr*dzr
                 dvdqmclassic(1,i,1) = dvdqmclassic(1,i,1) + dvdx
                 dvdqmclassic(2,i,1) = dvdqmclassic(2,i,1) + dvdy
                 dvdqmclassic(3,i,1) = dvdqmclassic(3,i,1) + dvdz
!
                 dvdqmclassic(1,j,1) = dvdqmclassic(1,j,1) - dvdx
                 dvdqmclassic(2,j,1) = dvdqmclassic(2,j,1) - dvdy
                 dvdqmclassic(3,j,1) = dvdqmclassic(3,j,1) - dvdz
!

                 call rho_Ni(dr,pma,dpadr_m,Smooth,DSmooth)
!
                 do k = 1, nb
                 eden_mclassic(j,k) = eden_mclassic(j,k) + pma
                 eden_mclassic(i,k) = eden_mclassic(i,k) + pma
                 end do
!
!        !      DERIVATIVES HERE
!
                 ilist = nlist_mclcl(i,1)
                 f_mclcl(1,ilist,i,1) = dpadr_m*dxr
                 f_mclcl(2,ilist,i,1) = dpadr_m*dyr
                 f_mclcl(3,ilist,i,1) = dpadr_m*dzr
!
         end if
         end do
         end do

        do i = 1,nclassic
        do k = 1,nrigid
!
              dx = qmclassic(1,i,1) - qmrigid(1,k,1)
              dy = qmclassic(2,i,1) - qmrigid(2,k,1)
              dz = qmclassic(3,i,1) - qmrigid(3,k,1)
          if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
            dx = dx-boxlx*nint(onboxlx*dx)
            dy = dy-boxly*nint(onboxly*dy)
            !dz = dz-boxlz*nint(onboxlz*dz)

          end if
              drsq = dx*dx + dy*dy + dz*dz
              if (drsq.lt.rc1sq) then
                 ! Add to neighbour list
!
                 nlist_mclrig(i,1) = nlist_mclrig(i,1)+1
                 list_mclrig(nlist_mclrig(i,1),i,1) = k
                 dr = dsqrt(drsq)
!
                 CALL SMOOTHING(dR,Smooth,DSmooth)
!
                 call pair_nini(dr,dv,dvdr,Smooth,DSmooth)
                 vmm = vmm + dv
                 dxr = dx/dr
                 dyr = dy/dr
                 dzr = dz/dr
                 dvdx = dvdr*dxr
                 dvdy = dvdr*dyr
                 dvdz = dvdr*dzr
!
                 dvdqmclassic(1,i,1) = dvdqmclassic(1,i,1) + dvdx
                 dvdqmclassic(2,i,1) = dvdqmclassic(2,i,1) + dvdy
                 dvdqmclassic(3,i,1) = dvdqmclassic(3,i,1) + dvdz
!
!
                 call rho_ni(dr,pma,dpadr_m,Smooth,DSmooth)
!
                    do  j = 1,nb
                    eden_mclassic(i,j) = eden_mclassic(i,j) + pma
                    eden_mrigid(k,j) = eden_mrigid(k,j) + pma
                    end do
!
!
!
!     !     DERIVATIVES HERE
!
                 ilist = nlist_mclrig(i,1)
                 f_mclrig(1,ilist,i,1) = dpadr_m*dxr
                 f_mclrig(2,ilist,i,1) = dpadr_m*dyr
                 f_mclrig(3,ilist,i,1) = dpadr_m*dzr
                 f_mrigcl(1,ilist,i,1) = dpadr_m*dxr
                 f_mrigcl(2,ilist,i,1) = dpadr_m*dyr
                 f_mrigcl(3,ilist,i,1) = dpadr_m*dzr
!
!
!
        end if
    end do
    end do


        do i = 1,nrigid
        do j = i+1,nrigid
              dx = qmrigid(1,i,1) - qmrigid(1,j,1)
              dy = qmrigid(2,i,1) - qmrigid(2,j,1)
              dz = qmrigid(3,i,1) - qmrigid(3,j,1)
          if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
           dx = dx-boxlx*nint(onboxlx*dx)
           dy = dy-boxly*nint(onboxly*dy)
           !dz = dz-boxlz*nint(onboxlz*dz)
          end if
              drsq = dx*dx + dy*dy + dz*dz
              if (drsq.lt.rc1sq) then
                 ! Add to neighbour list
!
                 nlist_mrigrig(i,1) = nlist_mrigrig(i,1)+1
!
                 dr = dsqrt(drsq)

                 CALL SMOOTHING(dR,Smooth,DSmooth)
                 call rho_Ni(dr,pma,dpadr_m,Smooth,DSmooth)
 !!
                    do k = 1, nb
                    eden_mrigid(j,k) = eden_mrigid(j,k) + pma
                    eden_mrigid(i,k) = eden_mrigid(i,k) + pma
                    end do
 !!      !       DERIVATIVES HERE
 !
                 ilist = nlist_mrigrig(i,1)
!                 f_mrigrig(1,ilist,i,1) = dpadr_m*dxr
!                 f_mrigrig(2,ilist,i,1) = dpadr_m*dyr
!                 f_mrigrig(3,ilist,i,1) = dpadr_m*dzr
!
        end if
        end do
        end do
!
! Embedding Energies
!
      fm = 0.d0
      fh = 0.d0
      df = 0.d0
!
!!!

        if(nclassic.ne.0) then   !!!debug
      do k = 1,nb
      do j = 1,nclassic
        pm = eden_mclassic(j,k)
        call embeding_nickel(pm,df,dfdp_mclassic(j,k))
        fm = fm + df
      enddo
      enddo
      end if
!
!      do j = 1,15
!        print*,j,eden_mrigid(j,1),'rhonickel'
!      enddo
        if(nrigid.ne.0) then
      do k = 1,nb
      do j = 1,nrigid
        pm = eden_mrigid(j,k)
        call embeding_nickel(pm,df,dfdp_mrigid(j,k))
        fm = fm + df
      enddo
      enddo
    end if

!       Embedding Forces for Ni-Ni
!        Classic-Classic
!
        if(nclassic.ne.0) then  !!! DEBUG
      do k = 1,1
      do j = 1,nclassic
        dfjdp = dfdp_mclassic(j,k)
        do ilist = 1,nlist_mclcl(j,k)
     if(ilist.gt.nclassic) write(*,*) ' classic-classic listing error'
           i = list_mclcl(ilist,j,k)
           dfidp = dfdp_mclassic(i,k)
           summ = dfidp + dfjdp
           fjx = f_mclcl(1,ilist,j,k)
           fjy = f_mclcl(2,ilist,j,k)
           fjz = f_mclcl(3,ilist,j,k)
           dvdqmclassic(1,j,k) = dvdqmclassic(1,j,k) + summ*fjx
           dvdqmclassic(2,j,k) = dvdqmclassic(2,j,k) + summ*fjy
           dvdqmclassic(3,j,k) = dvdqmclassic(3,j,k) + summ*fjz
           dvdqmclassic(1,i,k) = dvdqmclassic(1,i,k) - summ*fjx
           dvdqmclassic(2,i,k) = dvdqmclassic(2,i,k) - summ*fjy
           dvdqmclassic(3,i,k) = dvdqmclassic(3,i,k) - summ*fjz
        enddo
        enddo
        enddo
       end if

!
!     Classic-RIGID

              if(nclassic.ne.0.and.nrigid.ne.0) then
      do k = 1, 1
        do j = 1, nclassic
           dfjdp = dfdp_mclassic(j,k)
           do ilist = 1, nlist_mclrig(j,k)
    if(ilist.gt.nrigid) write(*,*) ' classic-rigid listing error'
              i = list_mclrig(ilist,j,k)
              dfidp = dfdp_mrigid(i,k)
              fjx = f_mrigcl(1,ilist,j,k)
              fjy = f_mrigcl(2,ilist,j,k)
              fjz = f_mrigcl(3,ilist,j,k)
              fix = f_mclrig(1,ilist,j,k)
              fiy = f_mclrig(2,ilist,j,k)
              fiz = f_mclrig(3,ilist,j,k)
              dvdqmclassic(1,j,k)   =  dvdqmclassic(1,j,k) + (dfidp*fjx + dfjdp*fix)
              dvdqmclassic(2,j,k)   =  dvdqmclassic(2,j,k) + (dfidp*fjy + dfjdp*fiy)
              dvdqmclassic(3,j,k)   =  dvdqmclassic(3,j,k) + (dfidp*fjz + dfjdp*fiz)
           enddo
        enddo
      enddo
      end if

      if(h_included) then
       do k = 1,nb
        do j = 1,1
           ph = eden_h(j,k)
           call embeding_hydrogen(ph,df,dfdp_h(j,k))
           fh = fh + df
        enddo
       enddo
       end if


      if(h_included) then
!    !  Classic-H
      if(nclassic.ne.0) then
      do k = 1,nb
        do j = 1,nclassic
           dfjdp = dfdp_mclassic(j,k)
           do ilist = 1,nlist_mclh(j,k)
    if(ilist.gt.nclassic) write(*,*) ' classic-H listing error'
              i = list_mclh(ilist,j,k)
              dfidp = dfdp_h(i,k)
              fjx = f_mclh(1,ilist,j,k)
              fjy = f_mclh(2,ilist,j,k)
              fjz = f_mclh(3,ilist,j,k)
              fix = f_hmcl(1,ilist,j,k)
              fiy = f_hmcl(2,ilist,j,k)
              fiz = f_hmcl(3,ilist,j,k)
              dvdqmclassic(1,j,1) = dvdqmclassic(1,j,1) + (dfidp*fjx + dfjdp*fix)
              dvdqmclassic(2,j,1) = dvdqmclassic(2,j,1) + (dfidp*fjy + dfjdp*fiy)
              dvdqmclassic(3,j,1) = dvdqmclassic(3,j,1) + (dfidp*fjz + dfjdp*fiz)
              dvdqh(1,i,k) = dvdqh(1,i,k) - (dfidp*fjx + dfjdp*fix)
              dvdqh(2,i,k) = dvdqh(2,i,k) - (dfidp*fjy + dfjdp*fiy)
              dvdqh(3,i,k) = dvdqh(3,i,k) - (dfidp*fjz + dfjdp*fiz)
           enddo
        enddo
       enddo
       end if
!

!      !RIGID-H
        if(nrigid.ne.0) then
      do k = 1,nb
        do j = 1,nrigid
           dfjdp = dfdp_mrigid(j,k)
           do ilist = 1,nlist_mrigh(j,k)
    if(ilist.gt.nrigid) write(*,*) ' rigid-H listing error'
              i = list_mrigh(ilist,j,k)
              dfidp = dfdp_h(i,k)
              fjx = f_mrigh(1,ilist,j,k)
              fjy = f_mrigh(2,ilist,j,k)
              fjz = f_mrigh(3,ilist,j,k)
              fix = f_hmrig(1,ilist,j,k)
              fiy = f_hmrig(2,ilist,j,k)
              fiz = f_hmrig(3,ilist,j,k)
              dvdqh(1,i,k) = dvdqh(1,i,k) - (dfidp*fjx + dfjdp*fix)
              dvdqh(2,i,k) = dvdqh(2,i,k) - (dfidp*fjy + dfjdp*fiy)
              dvdqh(3,i,k) = dvdqh(3,i,k) - (dfidp*fjz + dfjdp*fiz)
           enddo
        enddo
       enddo
       end if
!
!
      end if  ! h_incl
!
!

        v = vmh + vmm + fh + fm
!        print*,v,vmh, vmm, fh, fm,'v,vmh, vmm, fh, fm'
!
      ! Converting forces, energy and distances in a.u.

!        call atomic_units(nclassic,nrigid,nb,boxlx,boxly,boxlz,qmclassic, &
!                       qmrigid,qh,v,dvdqmclassic,dvdqmrigid,dvdqh, 1)

!       v=v/dble(nb)
       dvdqmrigid   = 0.d0
!
END subroutine eamnih_forces_umb
!
!
!

subroutine eamnih_energy(nclassic,nrigid,nb,boxlx,boxly,boxlz, &
        boundary,h_included,qmclassic,qmrigid,qh,v)
!
       Implicit None
!
!      EAM 4 Ni(100)/H. Input in a.u. Output in a.u. thanks to atomic_units subroutine
!
       integer :: i,j,k,ilist
       integer, intent(in) :: nclassic,nrigid, nb
       double precision :: dx,dy,dz
       double precision :: vmm,vmh,fm,fh !interaction potene mm, mh, embedding ene fm, fh
       double precision, intent(in) :: boxlx,boxly,boxlz
       double precision :: onboxlx,onboxly,onboxlz
       character(len=6),intent(in) :: boundary
       double precision :: rc1sq,rc1,drsq,dv,dvdr,dvdx,dvdy,dvdz,dr
       double precision :: pma,pha,ph,pm,df,dxr,dyr,dzr,dpadr_m,dpadr_h
       double precision :: dfidp,dfjdp,fjx,fjy,fjz,summ,fix,fiy,fiz
       double precision,intent(out) :: v
       logical, intent(in) :: h_included
!
!
       double precision :: f_mclh(1:3,1:1,nclassic,1:nb), &
        f_hmcl(1:3,1:1,nclassic,1:nb), &
        f_mrigh(1:3,1:1,nrigid,1:nb), &
        f_hmrig(1:3,1:1,nrigid,1:nb), &
        f_mclcl(1:3,nclassic,nclassic,1:1), &
        f_mclrig(1:3,nrigid,nclassic,1:1), &
        f_mrigcl(1:3,nrigid,nclassic,1:1)
!
       INTEGER :: nlist_mclh(nclassic,1:nb), &
        list_mclh(1,nclassic,1:nb),      &
        nlist_mrigh(nrigid,1:nb), &
        list_mrigh(1,nrigid,1:nb), &
        nlist_mclcl(nclassic,1:1), &
        list_mclcl(nclassic,nclassic,1:1), &
        nlist_mclrig(nclassic,1:1), &
        list_mclrig(nrigid,nclassic,1:1), &
        nlist_mrigrig(nrigid,1:1)
!
       double precision :: dfdp_mclassic(nclassic,1:nb), &
        dfdp_mrigid(nrigid,1:nb), &
        dfdp_h(1:1,1:nb)
!
       double precision,intent(in) :: qmclassic(1:3,nclassic,1:1), &
         qmrigid(1:3,nrigid,1:1)
!
       double precision,intent(in) :: qh(1:3,1:1,1:nb)
!
       double precision :: eden_mrigid(nrigid,1:nb), &
                          eden_mclassic(nclassic,1:nb), &
                          eden_h(1:1,1:nb)
!
        double precision :: smooth,dsmooth

!***************************************************************

!        nh = 1 ! this version is for 1 hydrogen atom only!
        vmm = 0.d0
        vmh = 0.d0
        v = 0.d0
        fm = 0.d0
        fh = 0.d0
        rc1 = 10.d0  !!!  DEBUG ! fixed cut off for Ni-Ni and H-Ni interactions (Angstroms)
        rc1sq = rc1*rc1
!
!
!
!
        dv = 0.d0
        dvdr = 0.d0

!      Cl-H
       nlist_mclh = 0
       list_mclh = 0
       f_mclh = 0.d0
       f_hmcl = 0.d0
!
!      RIG-H
       nlist_mrigh = 0
       list_mrigh = 0
       f_mrigh = 0.d0
       f_hmrig = 0.d0
!
!      Cl-Cl
       nlist_mclcl = 0
       list_mclcl = 0
       f_mclcl = 0.d0
 !
 !     Cl-RIG
       nlist_mclrig = 0
       list_mclrig = 0
       f_mclrig = 0.d0
       f_mrigcl = 0.d0
!
!      RIG-RIG
       nlist_mrigrig = 0
!       list_mrigrig = 0
!       f_mrigrig = 0.d0
!
!
       dfdp_mclassic = 0.d0
       dfdp_mrigid = 0.d0
       dfdp_h = 0.d0
!
       eden_mclassic = 0.d0
       eden_mrigid = 0.d0
       eden_h = 0.d0
!
!      ! Converting in eV/ANgs if input is not in Angstrom
!       call atomic_units(nclassic,nrigid,nb,boxlx,boxly,boxlz,qmclassic,&
!                    qmrigid,qh,v,dvdqmclassic,dvdqmrigid,dvdqh, -1)
        onboxlx = 1.d0/boxlx
        onboxly = 1.d0/boxly
        onboxlz = 1.d0/boxlz
!
    IF(nclassic.ne.0) then
        do i = 1,nb
        do j = 1, nclassic
!
              dx = qmclassic(1,j,1)-qh(1,1,i)
              dy = qmclassic(2,j,1)-qh(2,1,i)
              dz = qmclassic(3,j,1)-qh(3,1,i)
          if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
            dx = dx-boxlx*nint(onboxlx*dx)
            dy = dy-boxly*nint(onboxly*dy)
!            !dz = dz-boxlz*nint(onboxlz*dz)
           end if
              drsq = dx*dx + dy*dy + dz*dz
              if (drsq.lt.rc1sq) then
                 ! Add to neighbour list
!
                 nlist_mclh(j,i) = nlist_mclh(j,i)+1
                 list_mclh(nlist_mclh(j,i),j,i) = 1
!
                 dr = dsqrt(drsq)
                 CALL SMOOTHING(dR,Smooth,DSmooth)
                 call pair_NiH(dr,dv,dvdr,Smooth,DSmooth)
                 vmh = vmh + dv
                 dxr = dx/dr
                 dyr = dy/dr
                 dzr = dz/dr

!
                 call rho_hnew(dr,pha,dpadr_h,Smooth,DSmooth)
                 call rho_Ni(dr,pma,dpadr_m,Smooth,DSmooth)
!
                 eden_mclassic(j,i) = eden_mclassic(j,i) + pha
                 eden_h(1,i) = eden_h(1,i) + pma
!

        end if
        end do
        end do
    END IF

    IF(nrigid.ne.0) then
        do i = 1,nb
        do j = 1, nrigid
!
              dx = qmrigid(1,j,1)-qh(1,1,i)
              dy = qmrigid(2,j,1)-qh(2,1,i)
              dz = qmrigid(3,j,1)-qh(3,1,i)
          if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
            dx = dx-boxlx*nint(onboxlx*dx)
            dy = dy-boxly*nint(onboxly*dy)
!            !dz = dz-boxlz*nint(onboxlz*dz)
          end if
              drsq = dx*dx + dy*dy + dz*dz
              if (drsq.lt.rc1sq) then
                 ! Add to neighbour list
!
                 nlist_mrigh(j,i) = nlist_mrigh(j,i)+1
                 list_mrigh(nlist_mrigh(j,i),j,i) = 1
!
                 dr = dsqrt(drsq)
                 CALL SMOOTHING(dR,Smooth,DSmooth)
                 call pair_NiH(dr,dv,dvdr,Smooth,DSmooth)
                 vmh = vmh + dv
                 dxr = dx/dr
                 dyr = dy/dr
                 dzr = dz/dr

!
                 call rho_hnew(dr,pha,dpadr_h,Smooth,DSmooth)
                 call rho_Ni(dr,pma,dpadr_m,Smooth,DSmooth)
!
                 eden_mrigid(j,i) = eden_mrigid(j,i) + pha
                 eden_h(1,i) = eden_h(1,i) + pma
!

        end if
        end do
        end do
    END IF

!


      ! Classic-Classic

        do i=1,nclassic
        do j=i+1,nclassic
              dx = qmclassic(1,i,1) - qmclassic(1,j,1)
              dy = qmclassic(2,i,1) - qmclassic(2,j,1)
              dz = qmclassic(3,i,1) - qmclassic(3,j,1)
          if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
            dx = dx-boxlx*nint(onboxlx*dx)
            dy = dy-boxly*nint(onboxly*dy)
!            !dz = dz-boxlz*nint(onboxlz*dz)
          end if
              drsq = dx*dx + dy*dy + dz*dz
              if (drsq.lt.rc1sq) then
                 ! Add to neighbour list
!
                 nlist_mclcl(i,1) = nlist_mclcl(i,1)+1
                 list_mclcl(nlist_mclcl(i,1),i,1) = j
!
                 dr = dsqrt(drsq)
!
                 CALL SMOOTHING(dR,Smooth,DSmooth)
                call pair_nini(dr,dv,dvdr,Smooth,DSmooth)
                 vmm = vmm  + dv
                 dxr = dx/dr
                 dyr = dy/dr
                 dzr = dz/dr

                 call rho_Ni(dr,pma,dpadr_m,Smooth,DSmooth)
!
                 do k = 1, nb
                 eden_mclassic(j,k) = eden_mclassic(j,k) + pma
                 eden_mclassic(i,k) = eden_mclassic(i,k) + pma
                 end do
!

         end if
         end do
         end do

        do i = 1,nclassic
        do k = 1,nrigid
!
              dx = qmclassic(1,i,1) - qmrigid(1,k,1)
              dy = qmclassic(2,i,1) - qmrigid(2,k,1)
              dz = qmclassic(3,i,1) - qmrigid(3,k,1)
          if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
            dx = dx-boxlx*nint(onboxlx*dx)
            dy = dy-boxly*nint(onboxly*dy)
            !dz = dz-boxlz*nint(onboxlz*dz)

          end if
              drsq = dx*dx + dy*dy + dz*dz
              if (drsq.lt.rc1sq) then
                 ! Add to neighbour list
!
                 nlist_mclrig(i,1) = nlist_mclrig(i,1)+1
                 list_mclrig(nlist_mclrig(i,1),i,1) = k
                 dr = dsqrt(drsq)
!
                 CALL SMOOTHING(dR,Smooth,DSmooth)
!
                 call pair_nini(dr,dv,dvdr,Smooth,DSmooth)
                 vmm = vmm + dv
                 dxr = dx/dr
                 dyr = dy/dr
                 dzr = dz/dr

!
                 call rho_ni(dr,pma,dpadr_m,Smooth,DSmooth)
!
                    do  j = 1,nb
                    eden_mclassic(i,j) = eden_mclassic(i,j) + pma
                    eden_mrigid(k,j) = eden_mrigid(k,j) + pma
                    end do
!
        end if
    end do
    end do


        do i = 1,nrigid
        do j = i+1,nrigid
              dx = qmrigid(1,i,1) - qmrigid(1,j,1)
              dy = qmrigid(2,i,1) - qmrigid(2,j,1)
              dz = qmrigid(3,i,1) - qmrigid(3,j,1)
          if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
           dx = dx-boxlx*nint(onboxlx*dx)
           dy = dy-boxly*nint(onboxly*dy)
           !dz = dz-boxlz*nint(onboxlz*dz)
          end if
              drsq = dx*dx + dy*dy + dz*dz
              if (drsq.lt.rc1sq) then
                 ! Add to neighbour list
!
                 nlist_mrigrig(i,1) = nlist_mrigrig(i,1)+1
!
                 dr = dsqrt(drsq)

                 CALL SMOOTHING(dR,Smooth,DSmooth)
                 call rho_Ni(dr,pma,dpadr_m,Smooth,DSmooth)
 !!
                    do k = 1, nb
                    eden_mrigid(j,k) = eden_mrigid(j,k) + pma
                    eden_mrigid(i,k) = eden_mrigid(i,k) + pma
                    end do
!
        end if
        end do
        end do
!
! Embedding Energies
!
      fm = 0.d0
      fh = 0.d0
      df = 0.d0
!
!!!

        if(nclassic.ne.0) then   !!!debug
      do k = 1,nb
      do j = 1,nclassic
        pm = eden_mclassic(j,k)
        call embeding_nickel(pm,df,dfdp_mclassic(j,k))
        fm = fm + df
      enddo
      enddo
      end if
!
!      do j = 1,15
!        print*,j,eden_mrigid(j,1),'rhonickel'
!      enddo
        if(nrigid.ne.0) then
      do k = 1,nb
      do j = 1,nrigid
        pm = eden_mrigid(j,k)
        call embeding_nickel(pm,df,dfdp_mrigid(j,k))
        fm = fm + df
      enddo
      enddo
    end if


      if(h_included) then
       do k = 1,nb
        do j = 1,1
           ph = eden_h(j,k)
!           print*,j,k,ph,'rhoh'
           call embeding_hydrogen(ph,df,dfdp_h(j,k))
           fh = fh + df
        enddo
       enddo
       end if


!
!
        v = vmh + vmm + fh + fm
        print*,v,vmh, vmm, fh, fm,'v,vmh, vmm, fh, fm'
!
      ! Converting forces, energy and distances in a.u.

!        call atomic_units(nclassic,nrigid,nb,boxlx,boxly,boxlz,qmclassic, &
!                       qmrigid,qh,v,dvdqmclassic,dvdqmrigid,dvdqh, 1)

!       v=v/dble(nb)

END subroutine eamnih_energy

!!   EAM 4 PARAMETRIZATION for NI(100)/H by Yury V. Suleymanov
!!  PARAMETERS FOR THE NI DENSITY BASED SPHERICALLY ADJUSTED S- AND D- LIKE ORBITALS
!! TAKEN FROM SINGLE DETERMINANT HARTREE-fOCK CALCULATIONS.
      SUBROUTINE rho_Ni(r,p,dpdr,S, DSDR)
      IMPLICIT NONE
!      DOUBLE PRECISION :: C(1:10)
      DOUBLE PRECISION :: EPS(1:10) !10, the total number of valence electrons in NI
      DOUBLE PRECISION :: NI(1:10)
      DOUBLE PRECISION :: R,P,DPDR,PS,PD,PPS,PPD,DPDRS,DPDRD
      INTEGER :: I
      DOUBLE PRECISION :: S, DSDR,PI
!      DOUBLE PRECISION :: FACTORIAL(10)
      DOUBLE PRECISION :: COEF(10), RNI(1:10),DEXPEPS(1:10)
!!
      PI = DACOS(-1.D0)!3.141592653589793115997963468544185162d0 !
!!
      NI = 0.D0
      NI(1) = 1.D0
      NI(2) = 1.D0
      NI(3) = 2.D0
      NI(4) = 2.D0
      NI(5) = 3.D0
      NI(6) = 3.D0
      NI(7) = 4.D0
      NI(8) = 4.D0
      NI(9) = 3.D0
      NI(10) = 3.D0
!!
      EPS = 0.D0
      EPS(1) = 54.88885D0
      EPS(2) = 38.48431D0
      EPS(3) = 27.42703D0
      EPS(4) = 20.88204D0
      EPS(5) = 10.95707D0
      EPS(6) = 7.31958D0
      EPS(7) = 3.92650D0
      EPS(8) = 2.15289D0
      EPS(9) = 12.67582D0
      EPS(10) = 5.43253D0
!!
!      C = 0.D0
!      C(1) = -0.00389D0
!      C(2) = -0.02991D0
!      C(3) = -0.03189D0
!      C(4) =  0.15289D0
!      C(5) = -0.20048D0
!      C(6) = -0.05423D0
!      C(7) =  0.49292D0
!      C(8) =  0.61875D0
!      C(9) =  0.42120D0
!      C(10) = 0.70658D0

!      FACTORIAL(1) = 2.D0
!      FACTORIAL(2) = 2.D0
!      FACTORIAL(3) = 24.D0
!      FACTORIAL(4) = 24.D0
!      FACTORIAL(5) = 720.D0
!      FACTORIAL(6) = 720.D0
!      FACTORIAL(7) = 40320.D0
!      FACTORIAL(8) = 40320.D0
!      FACTORIAL(9) = 720.D0
!      FACTORIAL(10) = 720.D0

!     DO I = 1, 10
!     COEF(I) = C(I)*(1.d0/DSQRT(FACTORIAL(I)))*(2.d0*EPS(I))**(NI(I)+0.5D0)
!     END DO
    COEF(1) =      -3.163776491313133654159628349589d0
    COEF(2) =     -14.281438867475564791220676852390d0
    COEF(3) =    -145.067738725436015556624624878168d0
    COEF(4) =     351.787786416577375803171889856458d0
    COEF(5) =    -368.078386826022267541702603921294d0
    COEF(6) =     -24.259391731366324762575459317304d0
    COEF(7) =      26.162326823407976661428619991057d0
    COEF(8) =       2.197800222225523736341301628272d0
    COEF(9) =    1287.784878683165061374893411993980d0
    COEF(10) =     111.328804019789757262515195179731d0
!!
        DO I = 1, 10
        RNI(I) = R**(NI(I)-1.D0)
        DEXPEPS(I) = DEXP(-EPS(I)*R)
        END DO
!!



!!
      PS = 0.D0
!!     ATOMIC NI ELECTRONIC DENSITY for S-state
      DO I = 1, 8
      PS = PS + COEF(I)*RNI(I)*DEXPEPS(I)
      END DO
!!
      PPS = PS
      PS = PS*PS/(4.d0*PI)
!!
      PD = 0.D0
!!     ATOMIC NI ELECTRONIC DENSITY for d-state
      DO I = 9,10
      PD = PD + COEF(I)*RNI(I)*DEXPEPS(I)
      END DO
!!
      PPD = PD
      PD = PD*PD/(4.d0*PI)
!!
      P = 2.D0*PS + 8.D0*PD
!!
!!     NOW DERIVATIVES
      DPDRS = 0.D0
!!
      DO I = 1,2
      DPDRS = DPDRS + COEF(I)*(-EPS(I)*DEXPEPS(I))
      END DO
!!     GOOD CHECK IS TO CYCLE 1,10 THIS ONE
      DO I = 3, 8
      DPDRS = DPDRS + COEF(I)*((NI(I)-1.D0)*R**(NI(I)-2.D0)-EPS(I)*R**(NI(I)-1.D0))*DEXPEPS(I)
      END DO
!!
      DPDRS = DPDRS*PPS/(2.D0*PI)
!!
!!
      DPDRD = 0.D0
!!
      DO I = 9, 10
      DPDRD = DPDRD + COEF(I)*((NI(I)-1.D0)*R**(NI(I)-2.D0)-EPS(I)*RNI(I))*DEXPEPS(I)
      END DO
!!
      DPDRD = DPDRD*PPD/(2.D0*PI)
!!
      DPDR = 2.D0*DPDRS + 8.D0*DPDRD
!!
!!   NOW ADDING SMOOTHING FUNCTION
 !     CALL SMOOTHING(R,S,DSDR)
!!
      DPDR = DPDR*S+P*DSDR
      P = P*S
!!
!!    CHECKED
      RETURN
      END
!!
      SUBROUTINE SMOOTHING(R,S,DSDR)
      IMPLICIT NONE
      DOUBLE PRECISION :: R, RC, DELTA,S,DSDR
      DOUBLE PRECISION :: PI
!
      PI = DACOS(-1.D0)
      RC = 5.D0
      DELTA = 5.D0
!
      S = 0.D0
      DSDR = 0.D0
      IF(R.LT.RC) THEN
      S = 1.D0
      DSDR = 0.D0
      END IF
      IF(R.GE.RC.AND.R.LE.(RC+DELTA)) THEN
      S = (RC+DELTA-R)/DELTA - 1.D0/(2.D0*PI)*DSIN(PI*(2.D0*R-2.D0*RC-DELTA)/DELTA)
      DSDR = (-1.D0 - DCOS(PI*(2.D0*R-2.D0*RC-DELTA)/DELTA))/DELTA
      END IF
      IF(R.GT.(RC+DELTA)) THEN
      S = 0.D0
      DSDR = 0.D0
      END IF
!
      RETURN
      END
!
      SUBROUTINE EMBEDING_HYDROGEN(P,FH,DFHDP)
      IMPLICIT NONE
      DOUBLE PRECISION :: P, FH, DFHDP
      DOUBLE PRECISION :: ALPHAH, BETAH, BETAHP
!
      ALPHAH = -70.5461D0 !eV*Angstrom^3
      BETAH  = 6.9507D0 !Angstrom^3
      BETAHP = BETAH*P
!
      FH = ALPHAH*P*DEXP(-BETAHP)
!
      DFHDP = ALPHAH*DEXP(-BETAHP)*(1.D0-BETAHP)
!     CHECKED
      RETURN
      END
!
      SUBROUTINE EMBEDING_NICKEL(P,FNI,DFNIDP)
      IMPLICIT NONE
      DOUBLE PRECISION :: P, FNI, DFNIDP
      DOUBLE PRECISION :: A,B,C,AS,BS,CS,DS,DELTA,PC
      DOUBLE PRECISION :: ALPHA, BETA, GAMMANI
      DOUBLE PRECISION :: ALPHAP,BETAP,GAMMAP
    DOUBLE PRECISION :: DEXPALPHAP,DEXPBETAP,DEXPGAMMAP, P3,P2,PmPC,PmPC3,PmPC4
!
      A = -121.9314D0 !eV*Angstrom^3
      ALPHA = 0.0877D0 !Angstrom^3
      B = 6033.2396902178215896128676831722259520d0 !18.2047D0**3.d0 ! in paper18.2047D0 eV^1/2*Angstrom^3
      BETA = 10.5889D0 !Angstrom^3
      C = -205.2505D0 !eV*Angstrom^3
      GAMMANI = 51.8831D0 !Angstrom^3
      PC = 0.21D0 !Angstrom^-3
      DELTA = 0.01D0 !Angstrom^-3
      AS = -28379000000.d0 ! -2.8379D0*10.D0**(10.D0)  !eV*Angstrom^15
      BS = -756260000.d0 !-7.5626D0*10.D0**(8.D0) !eV*Angstrom^12
      CS = -5660400.d0 !-5.6604D0*10.D0**(6.D0) !eV*Angstrom^9
      DS = -19.0930D0 !eV
!
!
      IF(P.GE.0.D0.AND.P.LE.(PC-DELTA)) THEN
!
!
      ALPHAP=alpha*p
      BETAP=beta*p
      GAMMAP=gammaNI*p
!
        DEXPALPHAP=A*DEXP(-ALPHAP)
        DEXPBETAP=B*DEXP(-BETAP)
        DEXPGAMMAP=C*DEXP(-GAMMAP)
!
        P3= P**3.D0
        P2 = P**2.d0
!
      FNI = P*DEXPALPHAP+P3*DEXPBETAP+P*DEXPGAMMAP
      DFNIDP = DEXPALPHAP*(1.D0-ALPHAP) + DEXPBETAP*(3.D0*P2 - BETA*P3) + DEXPGAMMAP*(1.D0-GAMMAP)

      END IF
!
      IF(P.GT.(PC-DELTA).AND.P.LE.PC) THEN
        PmPC = P-PC
        PmPC3 = PmPC**3.d0
        PmPC4 = PmPC**4.d0
      FNI = AS*PmPC**5.D0 + BS*PmPC4 + CS*PmPC3 + DS
      DFNIDP =5.D0*AS*PmPC4 + 4.D0*BS*PmPC3 + 3.D0*CS*PmPC**2.D0
      END IF
!
      IF(P.GT.PC) THEN
      FNI = DS
      DFNIDP = 0.D0
      END IF
!
!
!
      RETURN
      END
!
      SUBROUTINE PAIR_NIH(R,PHINIH,DPHINIHDR,S, DSDR)
!! Calculates the short range pauirwise repulsive function phi and its gradient for NIH
      IMPLICIT NONE
!!
      DOUBLE PRECISION :: C
      DOUBLE PRECISION :: Z0H, BETAH, ALPHAH
      DOUBLE PRECISION :: Z0NI, BETANI, ALPHANI
      DOUBLE PRECISION :: R,PHINIH,DPHINIHDR
      DOUBLE PRECISION :: DZHDR,DZNIDR
      DOUBLE PRECISION :: ZH, ZNI
      DOUBLE PRECISION :: S, DSDR
      DOUBLE PRECISION :: ZHCOEF,ZNICOEF, A1, A2
!!
      C = 14.3888D0 !eV*Angstrom
      Z0H = 0.1959D0 !Z effective nuclear charge of H atom at distance R
      BETAH = 3.2108D0 ! Angstrom^-1
      ALPHAH = 1.7957D0  ! Angstrom^-1
      Z0NI = 10.0D0 !Z effective nuclear charge of Ni atom at distance R
      BETANI = 0.8957D0  ! Angstrom^-1
      ALPHANI = 1.8633D0  ! Angstrom^-1
!!
    ZHCOEF = Z0H*DEXP(-ALPHAH*R)
    A1 = 1.D0+BETAH*R
        ZH = ZHCOEF*A1
!!
    ZNICOEF = Z0NI*DEXP(-ALPHANI*R)
    A2= 1.D0+BETANI*R
        ZNI = ZNICOEF*A2
!!
      DZHDR = ZHCOEF*(BETAH-ALPHAH*A1)
      DZNIDR = ZNICOEF*(BETANI-ALPHANI*A2)
!!    NOW ADDING SMOOTHING FUNCTION
!!
      DZHDR = DZHDR*S + ZH*DSDR
      DZNIDR = DZNIDR*S + ZNI*DSDR
!!
      ZH = ZH*S
      ZNI = ZNI*S
!!
      !PHINIH = 0.D0
!!
      PHINIH = C*ZH*ZNI/R
!!
      DPHINIHDR = C*(DZHDR*ZNI + DZNIDR*ZH-ZNI*ZH/R)/R
!!
      RETURN
      END
!
!
!
      SUBROUTINE PAIR_NINI(R,PHININI,DPHININIDR,S, DSDR)
! Calculates the short range pauirwise repulsive function phi and its gradient for NINI
      IMPLICIT NONE
!
      DOUBLE PRECISION :: C, S, DSDR
      DOUBLE PRECISION :: Z0NI, BETANI, ALPHANI
      DOUBLE PRECISION :: DZNIDR1!, DZNIDR2
      DOUBLE PRECISION :: ZNI1!, ZNI2
      DOUBLE PRECISION :: R,PHININI,DPHININIDR
      DOUBLE PRECISION :: ZNICOEF,A1,ZNI1S,ZNI1SR
!
      C = 14.3888D0 !eV*Angstrom
!
      Z0NI = 10.0D0 !Z effective nuclear charge of Ni atom at distance R
      BETANI = 0.8957D0 ! Angstrom^-1
      ALPHANI = 1.8633D0 ! Angstrom^-1
!
!
      ZNICOEF = Z0NI*DEXP(-ALPHANI*R)
      A1=1.D0+BETANI*R
      ZNI1 = ZNICOEF*A1
!
      DZNIDR1 = ZNICOEF*(BETANI-ALPHANI*A1)
!
!!    NOW ADDING SMOOTHING FUNCTION
 !     CALL SMOOTHING(R,S,DSDR)
!!
     ! DZNIDR1 = DZNIDR1*S+DSDR*ZNI1
!!
      !ZNI1 = ZNI1*S
 !
!!
!
      !PHININI = 0.D0
!
      ZNI1S = ZNI1*S
      ZNI1SR=ZNI1S/R
      PHININI = C*ZNI1S*ZNI1SR

      DPHININIDR = (C*ZNI1SR)*(2.d0*(DZNIDR1*S+DSDR*ZNI1)-ZNI1SR)
!
      RETURN
      END
!!
      subroutine rho_hnew(r,p,dpdr,S, DSDR)
      implicit none
    ! ------------------------------------------------------------------
    ! H Electron Density and derivatives for NI(100)/H EAM4 parametrization: REV B 51, 9985 (1995)
    ! ------------------------------------------------------------------
      double precision :: r,delh,ch,dpdr,p
      double precision :: s, dsdr
!
      delh = 3.779452267842503321304548080661334097d0 !2.d0/0.5291772083D0                    !1.30927d0   Tom Markland's values
      ch =  2.148061616605061452389691112330183387d0 !1.d0/(PI*0.5291772083D0**3.d0)          !11.0025d0
!
      p = ch*dexp(-delh*r)
      dpdr = -delh*p
!
 !    CALL SMOOTHING(R,S,DSDR)
!
      DPDR=DPDR*S+DSDR*P
!
      P=P*S
!
      RETURN
      END


