
SUBROUTINE NIH2LEPS(XC, DX, ENERGY, GTEST, STEST)

  USE MODHESS
  USE COMMONS, ONLY: NATOMS, ZSYM
  USE KEY, ONLY: FROZEN, NFREEZE

  IMPLICIT NONE

!  integer          :: Natoms
  integer          :: nclassic,nrigid,i1,nh
  DOUBLE PRECISION :: XC(3*NATOMS), DX(3*NATOMS) !Cartesian coords, derivatives
  DOUBLE PRECISION :: ENERGY
  double precision :: boxlx,boxly,boxlz
  double precision :: qh(3,2),dvdqh(3,2)
  double precision :: qmclassic(3,NATOMS-NFREEZE-2),dvdqmclassic(3,NATOMS-NFREEZE-2)
  double precision :: qmrigid(3,NFREEZE),dvdqmrigid(3,NFREEZE)
!  CHARACTER(len=5) :: ZSYM(394)
  character(len=6) :: boundary
  logical          :: GTEST, STEST

    DOUBLE PRECISION   :: alphaHH, re, DHH, deltaHH
    DOUBLE PRECISION   :: nv,var
    DOUBLE PRECISION   :: alpha2, D2, r2
    DOUBLE PRECISION   :: deltaHNi, alphaH, DH, rH
    DOUBLE PRECISION   :: alphaNi, DNi, rNi,alpha2Ni, D2Ni, r2Ni
    DOUBLE PRECISION   :: alphaHcor,DHcor,rHcor

  nv    = 0.335d0 !effective number of 4s metal electrons
  var   = 1.3407d0
    ! !values for 0.335 4s electrons


  alphaHH = 0.d0 !A = 1.04435 a.u.
  DHH     = 0.d0      !eV
  re      = 0.d0 !.725d0 !A = 1.40083 a.u.

!  D2     = 0.d0  !eV
!  alpha2 = 0.d0  !A-1
!  r2     = 0.d0  !A 3.7335d0 bohr
    alpha2 = 1.7957D0  ! Angstrom^-1
    D2 = 0.1959D0 !Z effective nuclear charge of H atom at distance R
    r2 = 3.2108D0 ! Angstrom^-1

!  alphaH = 0.d0  !A⁻1 = 0.31138 bohr⁻1
!  DH     = 0.d0  ! 2.051d0 !ev
  rH     = 0.d0   !A   = 4.51141 bohr
!      alphaH  = 4.03690
!      DH      = -11.7724
!      rH      =  0.232963
alphaH  = 7.8549199D0
DH = -60.59949628D0

alphaHcor = 0.d0
DHcor = 0.d0
rHcor = 0.d0

!   r2ni    = 0.d0 !A   - equilibrium distance of two ni atoms at zero temperature (rT for tmp. chance)
!   alpha2Ni = 0.d0  !1/A - width of potential curve
!   D2Ni    = 0.d0  !eV  - depth of potential curve (DT for tmp. chance)
    alpha2Ni = 1.8633D0 ! Angstrom^-1
    D2Ni = 10.0D0 !Z effective nuclear charge of Ni atom at distance R
    r2Ni = 0.8957D0 ! Angstrom^-1


!   alphaNi = 0.d0
!   DNi    = 0.d0
!   rni    = 0.d0
    alphaNi =    3.1095475751174053D0
    DNi    =   -373.74637149155762D0
    rni    = -1.0250842641622047D0





!quantity of different kinds of atoms
nclassic = NATOMS-NFREEZE-2 !natoms-nfrozen = number of mobile atoms not including H
nrigid   =  NFREEZE !number of frozen/rigid atoms

!boxlength of the bulk (lattice constant = 3.52 A)
boxlx = 24.64d0 !length of the simulation box in x direction
boxly = 24.64d0 !length of the simulation box in y direction
boxlz =  7.04d0 !length of the simulation box in z direction

boundary   = 'bfixed' !'afixed', 'bfixed', or 'cfixed'. 'cfixed' or 'bfixed'- the lattices with 2D periodic boundary conditions

!initialise coords and forces
qh           = 0.d0  !position hydrogen
dvdqh        = 0.d0  !array for H atom forces
qmclassic    = 0.D0  !position metal classical
dvdqmclassic = 0.D0  !array for classic Ni atoms forces
qmrigid      = 0.D0  !position metal rigid
dvdqmrigid   = 0.D0  !array for rigid Ni atoms forces
dx           = 0.D0
nh = 0.D0

! The geometry of the mobile classical atoms
DO  I1 = 1,nclassic+2
   IF (ZSYM(I1).EQ.'Ni'.or.(ZSYM(I1).EQ.'NI')) THEN
      qmclassic(1,I1-2) = XC(3*I1-2) !I1-1 because of H excluded
      qmclassic(2,I1-2) = XC(3*I1-1)
      qmclassic(3,I1-2) = XC(3*I1)
   ELSE IF (ZSYM(I1).EQ.'H') THEN
      qh(1,I1) = XC(3*I1-2)
      qh(2,I1) = XC(3*I1-1)
      qh(3,I1) = XC(3*I1)
      nh=nh+1
   END IF
ENDDO

!Secondary zone
! The geometry of the rigid atoms
IF (nclassic+2 .LE. NATOMS) THEN
   DO  I1 = nclassic+3,NATOMS
      IF ((ZSYM(I1).EQ.'Ni').or.(ZSYM(I1).EQ.'NI')) THEN
         qmrigid(1,I1-nclassic-2) = XC(3*I1-2) !I1-nclassic-1 because of H excluded
         qmrigid(2,I1-nclassic-2) = XC(3*I1-1)
         qmrigid(3,I1-nclassic-2) = XC(3*I1)
      ELSE
         PRINT*, 'There are gas atoms in frozen zone'
         STOP
      END IF
   ENDDO
ENDIF

ENERGY = 0.0D0

call nih2_forces_umb(boxlx,boxly,boxlz, &
            boundary,qmclassic,qmrigid,qh,nclassic,nrigid,nh,&
            alphaHH, re, DHH, nv,var,alpha2, D2, r2, alpha2Ni, D2Ni, r2Ni,&
            alphaNi, DNi, rNi, alphaH, DH, rH,alphaHcor,DHcor,rHcor,&
            energy,dvdqmclassic,dvdqmrigid,dvdqh)


  !Gradient of Hydrogen
  DO I1=1,nh
     DX(3*I1-2)=dvdqh(1,I1)
     DX(3*I1-1)=dvdqh(2,I1)
     DX(3*I1)=dvdqh(3,I1)
  ENDDO
  !Gradient of mobile metal atoms
  DO I1 = 1,nclassic
     DX(3*(I1+nh) - 2) = dvdqmclassic(1,I1)
     DX(3*(I1+nh) - 1) = dvdqmclassic(2,I1)
     DX(3*(I1+nh))     = dvdqmclassic(3,I1)
  ENDDO

!    PRINT*, 'Gradient:'
!
!    DO I1 = 1,15,3 !(nclassic+1)*3,3
!        PRINT*, DX(I1), DX(I1+1), DX(I1+2)
!    ENDDO
!write(11,*) DX(6)
!  write(11,*) 'Gradients 2H,3Ni'
!  DO I1 = 1,5*3,3
!    write(11,*) DX(I1), DX(I1+1), DX(I1+2)
!  ENDDO
!  write(11,*) 'Energy = ', ENERGY,'eV'
!close(11)
!
!PRINT*, 'Energy = ', ENERGY,'eV'
!PRINT*, 'Done with ENE'

call writespec_xyz(17,XC)

END SUBROUTINE NIH2LEPS


subroutine nih2_forces_umb(boxlx,boxly,boxlz, &
     boundary,qmclassic,qmrigid,qh,nclassic,nrigid,nh,&
     alphaHH, re, DHH, nv,var,alpha2, D2, r2, alpha2Ni, D2Ni, r2Ni,&
     alphaNi, DNi, rNi, alphaH, DH, rH,alphaHcor,DHcor,rHcor,&
     v,dvdqmclassic,dvdqmrigid,dvdqh)


  !
  !      Mod leps eam Ni(100)/H2. Input in Output in Angstrom and eV
  !
  Implicit None

  DOUBLE PRECISION, intent(in)   :: alphaHH, re, DHH
  DOUBLE PRECISION, intent(in)   :: nv,var
  DOUBLE PRECISION, intent(in)   :: alpha2, D2, r2, alphaNi, DNi, rNi
  DOUBLE PRECISION, intent(in)   :: alphaH, DH, rH,alphaHcor,DHcor,rHcor
  DOUBLE PRECISION, intent(in)   :: alpha2Ni, D2Ni, r2Ni

  integer, intent(in) :: nclassic,nrigid, nh
  character(len=6),intent(in) :: boundary

  DOUBLE PRECISION, dimension(3*(nclassic+nrigid+nh)) :: Dx0,DX1,DX2,DX3,DX4,DX5
  DOUBLE PRECISION, dimension(nclassic+nrigid+nh)   :: edens

  integer :: i,j,k,ilist, I1
  DOUBLE PRECISION :: dx,dy,dz
  DOUBLE PRECISION :: vmm,vmh,fac,faci,fm,fh
  DOUBLE PRECISION :: boxlx,boxly,boxlz,onboxlx,onboxly,onboxlz
  DOUBLE PRECISION :: rc1sq,rc1,drsq,dv,dvdr,dvdx,dvdy,dvdz,dr
  DOUBLE PRECISION :: pma,pha,ph,pm,df,dxr,dyr,dzr,dpadr_m,dpadr_h
  DOUBLE PRECISION :: dfidp,dfjdp,fjx,fjy,fjz,sum,fix,fiy,fiz
  DOUBLE PRECISION :: v,vhh
  DOUBLE PRECISION :: rc2sq, rc2

  !
  !
  real(8) :: f_mclh(1:3,1:nh,nclassic), &
       f_hmcl(1:3,1:nh,nclassic), &
       f_mrigh(1:3,1:nh,nrigid), &
       f_hmrig(1:3,1:nh,nrigid), &
       f_mclcl(1:3,nclassic,nclassic), &
       f_mclrig(1:3,nrigid,nclassic), &
       f_mrigcl(1:3,nrigid,nclassic), &
       f_mrigrig(1:3,nrigid,nrigid)
  !
  !
  INTEGER :: nlist_mclh(nclassic), &
       list_mclh(1:nh,nclassic),      &
       nlist_mrigh(nrigid), &
       list_mrigh(1:nh,nrigid), &
       nlist_mclcl(nclassic), &
       list_mclcl(nclassic,nclassic), &
       nlist_mclrig(nclassic), &
       list_mclrig(nrigid,nclassic), &
       nlist_mrigrig(nrigid), &
       list_mrigrig(nrigid,nrigid)
  !
  real(8) :: dfdp_mclassic(nclassic), &
       dfdp_mrigid(nrigid), &
       dfdp_h(nh)
  !
  real(8) :: qmclassic(1:3,nclassic), &
       qmrigid(1:3,nrigid)
  !
  real(8) ::  qh(1:3,1:nh)
  !
  real(8) :: dvdqmclassic(1:3,nclassic), &
       dvdqmrigid(1:3,nrigid), &
       dvdqh(1:3,1:nh)
  !
  real(8) :: eden_mrigid(nrigid), &
       eden_mclassic(nclassic), &
       eden_h(nh)
  !
  !
  REAL(8) :: fs, sq, s, halffs,ds,vb
  REAL(8) :: sxc, syc
  REAL(8) :: smooth,dsmooth
  !
  !********************
  !From new LEPS
  double precision :: JHH, QHH,rHNi,QHNi_HNi,JHNi(1:nh),NIJH(1:nh)
  double precision :: dJHHdr, dJHH(1:3)
  double precision :: dQHHdr,dJHNidHNi, dQHNidHNi

  QHH =0.D0
  JHH =0.D0
  dQHHdr =0.D0
  dJHHdr =0.D0

  !********************
  Dx0(:)=0.D0
  DX1(:)=0.D0
  DX2(:)=0.D0
  DX3(:)=0.D0
  DX4(:)=0.D0
  DX5(:)=0.D0

  vmm = 0.d0
  vmh = 0.d0
  vhh = 0.d0
  v = 0.d0
  fm = 0.d0
  fh = 0.d0
  dv = 0.d0
  dvdr = 0.d0

  rc1 = 2.6d0  !!!  DEBUG ! fixed cut off for Ni-Ni and H-Ni interactions (Angstroms)
  rc1sq = rc1*rc1
  rc2 = 2.6d0  !!!  DEBUG  ! fixed cut off for Ni-Ni interactions (Angstroms)
  rc2sq = rc2*rc2



  dvdqmclassic(:,:) = 0.d0
  dvdqmrigid(:,:) = 0.d0
  dvdqh(:,:) = 0.d0
  !
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
  list_mrigrig = 0
  f_mrigrig = 0.d0
  !
  !
  dfdp_mclassic = 0.d0
  dfdp_mrigid = 0.d0
  !       dfdp_h = 0.d0
  !
  eden_mclassic = 0.d0
  eden_mrigid = 0.d0
  eden_h = 0.d0

  onboxlx = 1.d0/boxlx
  onboxly = 1.d0/boxly
  onboxlz = 1.d0/boxlz
  !


  !        !calculation of the distance of the two H2 atoms
  !            IF(nh.ne.0) then
  !                do i = 1,nh
  !                do j = i+1,nh
  !                    dx = qh(1,j)-qh(1,i)
  !                    dy = qh(2,j)-qh(2,i)
  !                    dz = qh(3,j)-qh(3,i)
  !                    if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
  !                        dx = dx-boxlx*nint(onboxlx*dx)
  !                        dy = dy-boxly*nint(onboxly*dy)
  !                        !dz = dz-boxlz*nint(onboxlz*dz)  !DEBUG
  !                    end if
  !                    drsq = dx*dx + dy*dy + dz*dz
  !
  !                 if (drsq.lt.rc1sq) then
  !                         ! Add to neighbour list
  !        !
  !        !            nlist_mclh(j,i) = nlist_mclh(j,i)+1
  !        !            list_mclh(nlist_mclh(j,i),j,i) = j
  !        !
  !
  !                    dr = dsqrt(drsq)
  !
  !                    call QHH_JHH(dr, QHH, dQHHdr,alphaHH, re, DHH, rc1)!SAME CUTOFF AS HNI FOR HH
  !                    vhh = v + QHH
  !                    dxr = dx/dr
  !                    dyr = dy/dr
  !                    dzr = dz/dr
  !
  !                    dvdqh(1,j)=dvdqh(1,j)+dQHHdr*dxr
  !                    dvdqh(2,j)=dvdqh(2,j)+dQHHdr*dyr
  !                    dvdqh(3,j)=dvdqh(3,j)+dQHHdr*dzr
  !
  !                    dvdqh(1,i)=dvdqh(1,i)-dQHHdr*dxr
  !                    dvdqh(2,i)=dvdqh(2,i)-dQHHdr*dyr
  !                    dvdqh(3,i)=dvdqh(3,i)-dQHHdr*dzr
  !
  !                 end if
  !                end do
  !                end do
  !            endif
  !
  !
  !
  !calculation of all distances between H and Ni(move)
  IF(nclassic.ne.0) then
     do i = 1,nh
        do j = 1, nclassic
           !
           dx = qmclassic(1,j)-qh(1,i)
           dy = qmclassic(2,j)-qh(2,i)
           dz = qmclassic(3,j)-qh(3,i)
           if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
              dx = dx-boxlx*nint(onboxlx*dx)
              dy = dy-boxly*nint(onboxly*dy)
              !dz = dz-boxlz*nint(onboxlz*dz)  !DEBUG
           end if
           drsq = dx*dx + dy*dy + dz*dz
           if (drsq.lt.rc1sq) then
              ! Add to neighbour list
              !
              nlist_mclh(j) = nlist_mclh(j)+1
              list_mclh(nlist_mclh(j),j) = i
              !
              dr = dsqrt(drsq)

              call QHNi_JHNi(dr,dv,dvdr,&
                   alpha2, D2, r2,alpha2Ni, D2Ni, r2Ni)

              vmh = vmh + dv
              dxr = dx/dr
              dyr = dy/dr
              dzr = dz/dr
              dvdx = dvdr*dxr
              dvdy = dvdr*dyr
              dvdz = dvdr*dzr
              !
              dvdqmclassic(1,j) = dvdqmclassic(1,j) + dvdx
              dvdqmclassic(2,j) = dvdqmclassic(2,j) + dvdy
              dvdqmclassic(3,j) = dvdqmclassic(3,j) + dvdz
              dvdqh(1,i) = dvdqh(1,i)  - dvdx
              dvdqh(2,i) = dvdqh(2,i)  - dvdy
              dvdqh(3,i) = dvdqh(3,i)  - dvdz
              !
              call rho_hydrogen(dr,pha,dpadr_h)
              call rho_Nickel(dr,pma,dpadr_m)
              !
              eden_mclassic(j) = eden_mclassic(j) + pha
              eden_h(i) = eden_h(i) + pma

              !              DERIVATIVES HERE
              !
              ilist = nlist_mclh(j)
              f_mclh(1,ilist,j) = dpadr_m*dxr
              f_mclh(2,ilist,j) = dpadr_m*dyr
              f_mclh(3,ilist,j) = dpadr_m*dzr
              f_hmcl(1,ilist,j) = dpadr_h*dxr
              f_hmcl(2,ilist,j) = dpadr_h*dyr
              f_hmcl(3,ilist,j) = dpadr_h*dzr
              !
           end if
        end do
     end do
  END IF

  !Gradient of Hydrogen
  DO I1=1,nh
     DX0(3*I1-2)=dvdqh(1,I1)
     DX0(3*I1-1)=dvdqh(2,I1)
     DX0(3*I1)=dvdqh(3,I1)
  ENDDO
  !Gradient of mobile metal atoms
  DO I1 = 1,nclassic
     DX0(3*(I1+nh) - 2) = dvdqmclassic(1,I1)
     DX0(3*(I1+nh) - 1) = dvdqmclassic(2,I1)
     DX0(3*(I1+nh))     = dvdqmclassic(3,I1)
  ENDDO

  !calculation of all distances between H and Ni(no_move)
  IF(nrigid.ne.0) then
     do i = 1,nh
        do j = 1, nrigid
           !
           dx = qmrigid(1,j)-qh(1,i)
           dy = qmrigid(2,j)-qh(2,i)
           dz = qmrigid(3,j)-qh(3,i)
           if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
              dx = dx-boxlx*nint(onboxlx*dx)
              dy = dy-boxly*nint(onboxly*dy)
              !dz = dz-boxlz*nint(onboxlz*dz)  !DEBUG
           end if
           drsq = dx*dx + dy*dy + dz*dz
           ! Add to neighbour list
           !
           if (drsq.lt.rc1sq) then
              nlist_mrigh(j) = nlist_mrigh(j)+1
              list_mrigh(nlist_mrigh(j),j) = i
              !
              dr = dsqrt(drsq)

              call QHNi_JHNi(dr,dv,dvdr,&
                   alpha2, D2, r2,alpha2Ni, D2Ni, r2Ni)

              vmh = vmh + dv
              dxr = dx/dr
              dyr = dy/dr
              dzr = dz/dr
              dvdx = dvdr*dxr
              dvdy = dvdr*dyr
              dvdz = dvdr*dzr
              !
              dvdqh(1,i) = dvdqh(1,i)  - dvdx
              dvdqh(2,i) = dvdqh(2,i)  - dvdy
              dvdqh(3,i) = dvdqh(3,i)  - dvdz

              call rho_hydrogen(dr,pha,dpadr_h)
              call rho_Nickel(dr,pma,dpadr_m)
              !
              eden_mrigid(j) = eden_mrigid(j) + pha
              eden_h(i) = eden_h(i) + pma

              !              DERIVATIVES HERE
              !
              ilist = nlist_mrigh(j)
              f_mrigh(1,ilist,j) = dpadr_m*dxr
              f_mrigh(2,ilist,j) = dpadr_m*dyr
              f_mrigh(3,ilist,j) = dpadr_m*dzr
              f_hmrig(1,ilist,j) = dpadr_h*dxr
              f_hmrig(2,ilist,j) = dpadr_h*dyr
              f_hmrig(3,ilist,j) = dpadr_h*dzr
              !
              !
              !
           end if
        end do
     end do
  END IF

  !Gradient of Hydrogen
  DO I1=1,nh
     DX1(3*I1-2)=dvdqh(1,I1)
     DX1(3*I1-1)=dvdqh(2,I1)
     DX1(3*I1)=dvdqh(3,I1)
  ENDDO
  !Gradient of mobile metal atoms
  DO I1 = 1,nclassic
     DX1(3*(I1+nh) - 2) = dvdqmclassic(1,I1)
     DX1(3*(I1+nh) - 1) = dvdqmclassic(2,I1)
     DX1(3*(I1+nh))     = dvdqmclassic(3,I1)
  ENDDO


  !
  !       Classic-Classic nickel
  !
  if(nclassic.ne.0) then
     do i=1,nclassic
        do j=i+1,nclassic
           dx = qmclassic(1,i) - qmclassic(1,j)
           dy = qmclassic(2,i) - qmclassic(2,j)
           dz = qmclassic(3,i) - qmclassic(3,j)
           if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
              dx = dx-boxlx*nint(onboxlx*dx)
              dy = dy-boxly*nint(onboxly*dy)
              !dz = dz-boxlz*nint(onboxlz*dz)  !DEBUG
           end if
           drsq = dx*dx + dy*dy + dz*dz
           if (drsq.lt.rc2sq) then
              ! Add to neighbour list
              !
              nlist_mclcl(i) = nlist_mclcl(i)+1
              list_mclcl(nlist_mclcl(i),i) = j
              !
              dr = dsqrt(drsq)
              !
              !CALL SMOOTHING(dR,Smooth,DSmooth)
              !call pair_nini(dr,dv,dvdr,Smooth,DSmooth)

              call QNiNi_JNiNi(dr,dv,dvdr,&
                   alpha2Ni, D2Ni, R2NI)

              vmm = vmm  + dv
              dxr = dx/dr
              dyr = dy/dr
              dzr = dz/dr
              dvdx = dvdr*dxr
              dvdy = dvdr*dyr
              dvdz = dvdr*dzr
              dvdqmclassic(1,i) = dvdqmclassic(1,i) + dvdx
              dvdqmclassic(2,i) = dvdqmclassic(2,i) + dvdy
              dvdqmclassic(3,i) = dvdqmclassic(3,i) + dvdz
              !
              dvdqmclassic(1,j) = dvdqmclassic(1,j) - dvdx
              dvdqmclassic(2,j) = dvdqmclassic(2,j) - dvdy
              dvdqmclassic(3,j) = dvdqmclassic(3,j) - dvdz

              call rho_Nickel(dr,pma,dpadr_m)
              !
              eden_mclassic(j) = eden_mclassic(j) + pma
              eden_mclassic(i) = eden_mclassic(i) + pma
              !
              !              DERIVATIVES HERE
              !
              ilist = nlist_mclcl(i)
              f_mclcl(1,ilist,i) = dpadr_m*dxr
              f_mclcl(2,ilist,i) = dpadr_m*dyr
              f_mclcl(3,ilist,i) = dpadr_m*dzr
              !
           end if
        end do
     end do
  endif

  !Gradient of Hydrogen
  DO I1=1,nh
     DX2(3*I1-2)=dvdqh(1,I1)
     DX2(3*I1-1)=dvdqh(2,I1)
     DX2(3*I1)=dvdqh(3,I1)
  ENDDO
  !Gradient of mobile metal atoms
  DO I1 = 1,nclassic
     DX2(3*(I1+nh) - 2) = dvdqmclassic(1,I1)
     DX2(3*(I1+nh) - 1) = dvdqmclassic(2,I1)
     DX2(3*(I1+nh))     = dvdqmclassic(3,I1)
  ENDDO

  !
  !       Classic-rigid nickel
  !
  do i = 1,nclassic
     do k = 1,nrigid
        !
        dx = qmclassic(1,i) - qmrigid(1,k)
        dy = qmclassic(2,i) - qmrigid(2,k)
        dz = qmclassic(3,i) - qmrigid(3,k)
        if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
           dx = dx-boxlx*nint(onboxlx*dx)
           dy = dy-boxly*nint(onboxly*dy)
           !dz = dz-boxlz*nint(onboxlz*dz)  !DEBUG
        end if
        drsq = dx*dx + dy*dy + dz*dz
        if (drsq.lt.rc2sq) then
           ! Add to neighbour list
           !
           nlist_mclrig(i) = nlist_mclrig(i)+1
           list_mclrig(nlist_mclrig(i),i) = k
           dr = dsqrt(drsq)
           !
           !CALL SMOOTHING(dR,Smooth,DSmooth)
           !
           !call pair_nini(dr,dv,dvdr,Smooth,DSmooth)
           call QNiNi_JNiNi(dr,dv,dvdr,&
                alpha2NI, D2Ni, r2Ni)

           vmm = vmm + dv
           dxr = dx/dr
           dyr = dy/dr
           dzr = dz/dr
           dvdx = dvdr*dxr
           dvdy = dvdr*dyr
           dvdz = dvdr*dzr
           !
           dvdqmclassic(1,i) = dvdqmclassic(1,i) + dvdx
           dvdqmclassic(2,i) = dvdqmclassic(2,i) + dvdy
           dvdqmclassic(3,i) = dvdqmclassic(3,i) + dvdz

           call rho_Nickel(dr,pma,dpadr_m)
           !
           eden_mclassic(i) = eden_mclassic(i) + pma
           eden_mrigid(k) = eden_mrigid(k) + pma
           !
           !
           !
           !          DERIVATIVES HERE
           !
           ilist = nlist_mclrig(i)
           f_mclrig(1,ilist,i) = dpadr_m*dxr
           f_mclrig(2,ilist,i) = dpadr_m*dyr
           f_mclrig(3,ilist,i) = dpadr_m*dzr
           f_mrigcl(1,ilist,i) = dpadr_m*dxr
           f_mrigcl(2,ilist,i) = dpadr_m*dyr
           f_mrigcl(3,ilist,i) = dpadr_m*dzr
           !
           !
           !
        end if
     end do
  end do

  !Gradient of Hydrogen
  DO I1=1,nh
     DX3(3*I1-2)=dvdqh(1,I1)
     DX3(3*I1-1)=dvdqh(2,I1)
     DX3(3*I1)=dvdqh(3,I1)
  ENDDO
  !Gradient of mobile metal atoms
  DO I1 = 1,nclassic
     DX3(3*(I1+nh) - 2) = dvdqmclassic(1,I1)
     DX3(3*(I1+nh) - 1) = dvdqmclassic(2,I1)
     DX3(3*(I1+nh))     = dvdqmclassic(3,I1)
  ENDDO

  !
  !       rigid-rigid nickel
  !
  IF(nrigid.ne.0) then
     do i = 1,nrigid
        do j = i+1,nrigid
           dx = qmrigid(1,i) - qmrigid(1,j)
           dy = qmrigid(2,i) - qmrigid(2,j)
           dz = qmrigid(3,i) - qmrigid(3,j)
           if(boundary.eq.'bfixed'.or.boundary.eq.'cfixed') then
              dx = dx-boxlx*nint(onboxlx*dx)
              dy = dy-boxly*nint(onboxly*dy)
              ! dz = dz-boxlz*nint(onboxlz*dz)  !DEBUG
           end if
           drsq = dx*dx + dy*dy + dz*dz
           if (drsq.lt.rc2sq) then
              ! Add to neighbour list
              !
              nlist_mrigrig(i) = nlist_mrigrig(i)+1
              list_mrigrig(nlist_mrigrig(i),i) = j
              !
              dr = dsqrt(drsq)

              dxr = dx/dr
              dyr = dy/dr
              dzr = dz/dr

              call RHO_NICKEL(dr,pma,dpadr_m)
              !
              eden_mrigid(j) = eden_mrigid(j) + pma
              eden_mrigid(i) = eden_mrigid(i) + pma
              !
              !
              !
              !              DERIVATIVES HERE
              !
              ilist = nlist_mrigrig(i)
              f_mrigrig(1,ilist,i) = dpadr_m*dxr
              f_mrigrig(2,ilist,i) = dpadr_m*dyr
              f_mrigrig(3,ilist,i) = dpadr_m*dzr
              !
           end if
        end do
     end do
  endif


  !Gradient of Hydrogen
  DO I1=1,nh
      DX4(3*I1-2)=dvdqh(1,I1)
      DX4(3*I1-1)=dvdqh(2,I1)
      DX4(3*I1)=dvdqh(3,I1)
  ENDDO
  !Gradient of mobile metal atoms
  DO I1 = 1,nclassic
      DX4(3*(I1+nh) - 2) = dvdqmclassic(1,I1)
      DX4(3*(I1+nh) - 1) = dvdqmclassic(2,I1)
      DX4(3*(I1+nh))     = dvdqmclassic(3,I1)
  ENDDO

  !  fjx = 0.D0
  !  fjy = 0.D0
  !  fjz = 0.D0
  !  fix = 0.D0
  !  fiy = 0.D0
  !  fiz = 0.D0

  !dvdqmclassic(:,:) = 0.d0
  !dvdqh(:,:) = 0.d0

  !
  ! Embedding Energies
  !
  fm = 0.d0
  fh = 0.d0
  df = 0.d0
  !
  if(nclassic.ne.0) then
      do j = 1,nclassic
          pm = eden_mclassic(j)
          call QNiNi_JNiNi_RHO(pm,df,dfdp_mclassic(j),alphaNi, DNi, rNi)
          fm = fm + df
          edens(nh+j)=eden_mclassic(j)
      enddo
  end if
  !
  !
  if(nrigid.ne.0) then
      do j = 1,nrigid
          pm = eden_mrigid(j)
          call QNiNi_JNiNi_RHO(pm,df,dfdp_mrigid(j),alphaNi, DNi, rNi)
          fm = fm + df
          edens(nh+nclassic+j)=eden_mrigid(j)
      enddo
  end if
  !
  !        Classic-Classic
  !
  if(nclassic.ne.0) then
      do j = 1,nclassic
          dfjdp = dfdp_mclassic(j)
          do ilist = 1,nlist_mclcl(j)
              if(ilist.gt.nclassic) write(*,*) ' classic-classic listing error'
              i = list_mclcl(ilist,j)
              dfidp = dfdp_mclassic(i)
              sum = dfidp + dfjdp
              fjx = f_mclcl(1,ilist,j)
              fjy = f_mclcl(2,ilist,j)
              fjz = f_mclcl(3,ilist,j)
              dvdqmclassic(1,j) = dvdqmclassic(1,j) + sum*fjx
              dvdqmclassic(2,j) = dvdqmclassic(2,j) + sum*fjy
              dvdqmclassic(3,j) = dvdqmclassic(3,j) + sum*fjz
              dvdqmclassic(1,i) = dvdqmclassic(1,i) - sum*fjx
              dvdqmclassic(2,i) = dvdqmclassic(2,i) - sum*fjy
              dvdqmclassic(3,i) = dvdqmclassic(3,i) - sum*fjz
          enddo
      enddo
  end if
  !
  !     Classic-RIGID
  !
  if(nclassic.ne.0.and.nrigid.ne.0) then
      do j = 1, nclassic
          dfjdp = dfdp_mclassic(j)
          do ilist = 1, nlist_mclrig(j)
              if(ilist.gt.nrigid) write(*,*) ' classic-rigid listing error'
              i = list_mclrig(ilist,j)
              dfidp = dfdp_mrigid(i)
              fjx = f_mrigcl(1,ilist,j)
              fjy = f_mrigcl(2,ilist,j)
              fjz = f_mrigcl(3,ilist,j)
              fix = f_mclrig(1,ilist,j)
              fiy = f_mclrig(2,ilist,j)
              fiz = f_mclrig(3,ilist,j)
              dvdqmclassic(1,j)   =  dvdqmclassic(1,j) + (dfidp*fjx + dfjdp*fix)
              dvdqmclassic(2,j)   =  dvdqmclassic(2,j) + (dfidp*fjy + dfjdp*fiy)
              dvdqmclassic(3,j)   =  dvdqmclassic(3,j) + (dfidp*fjz + dfjdp*fiz)
          enddo
      enddo
  end if


  !Gradient of Hydrogen
  DO I1=1,nh
      DX5(3*I1-2)=dvdqh(1,I1)
      DX5(3*I1-1)=dvdqh(2,I1)
      DX5(3*I1)=dvdqh(3,I1)
  ENDDO
  !Gradient of mobile metal atoms
  DO I1 = 1,nclassic
      DX5(3*(I1+nh) - 2) = dvdqmclassic(1,I1)
      DX5(3*(I1+nh) - 1) = dvdqmclassic(2,I1)
      DX5(3*(I1+nh))     = dvdqmclassic(3,I1)
  ENDDO

  !  fjx = 1.D0
  !  fjy = 1.D0
  !  fjz = 1.D0
  !  fix = 0.D0
  !  fiy = 0.D0
  !  fiz = 0.D0

  !  dvdqmclassic(:,:)=0.D0
  !  dvdqh(:,:)=0.D0

  !The density dependent parts
  !hydrogen density
  do k = 1,nh
      ph = eden_h(k)
      call QHNi_JHNi_RHO(ph,df,dfdp_h(k),&
      alphaH, DH,rH,alphaHcor,DHcor,rHcor)
      fh = fh + df
      edens(k)=eden_h(k)
     !dfdp_h(k)=1.D0
     !dfdp_h(k)=sign(1.D0,dfdp_h(k))
  enddo
  !
  !Classic-H
  !
  if(nclassic.ne.0) then
      do j = 1, nclassic
          dfjdp = dfdp_mclassic(j)
          do ilist = 1,nlist_mclh(j)
              if(ilist.gt.nclassic) write(*,*) ' classic-H listing error'
              i = list_mclh(ilist,j)
              dfidp =  dfdp_h(i)
              fjx = f_mclh(1,ilist,j)
              fjy = f_mclh(2,ilist,j)
              fjz = f_mclh(3,ilist,j)
              fix = f_hmcl(1,ilist,j)
              fiy = f_hmcl(2,ilist,j)
              fiz = f_hmcl(3,ilist,j)
              dvdqmclassic(1,j)   =  dvdqmclassic(1,j) + (dfidp*fjx + dfjdp*fix)
              dvdqmclassic(2,j)   =  dvdqmclassic(2,j) + (dfidp*fjy + dfjdp*fiy)
              dvdqmclassic(3,j)   =  dvdqmclassic(3,j) + (dfidp*fjz + dfjdp*fiz)
              dvdqh(1,i) = dvdqh(1,i) - (dfidp*fjx + dfjdp*fix)
              dvdqh(2,i) = dvdqh(2,i) - (dfidp*fjy + dfjdp*fiy)
              dvdqh(3,i) = dvdqh(3,i) - (dfidp*fjz + dfjdp*fiz)
             !              write(13,*) j, dvdqh(3,i)
          enddo
      enddo
  end if

  !      RIGID-H
  if(nrigid.ne.0) then
      do j = 1, nrigid
          dfjdp = dfdp_mrigid(j)
          do ilist = 1,nlist_mrigh(j)
              if(ilist.gt.nrigid) write(*,*) ' rigid-H listing error'
              i = list_mrigh(ilist,j)
              dfidp =  dfdp_h(i)
              fjx = f_mclh(1,ilist,j)
              fjy = f_mclh(2,ilist,j)
              fjz = f_mclh(3,ilist,j)
              fix = f_hmcl(1,ilist,j)
              fiy = f_hmcl(2,ilist,j)
              fiz = f_hmcl(3,ilist,j)
              dvdqh(1,i) = dvdqh(1,i) - (dfidp*fjx + dfjdp*fix)
              dvdqh(2,i) = dvdqh(2,i) - (dfidp*fjy + dfjdp*fiy)
              dvdqh(3,i) = dvdqh(3,i) - (dfidp*fjz + dfjdp*fiz)
             !                    write(13,*) j, dvdqh(3,i)

          enddo
      enddo
  end if

!
  !
  v = vmh + vmm + vhh + fh + fm
  !
  !dvdqmquantum = 0.d0
  !dvdqmclassic = 0.d0
  dvdqmrigid   = 0.d0
  !dvdqh = 0.d0
  !v = 0.d0




END subroutine nih2_forces_umb





! **********************************************************************
! ***                       Utilities                                ***
! **********************************************************************

! **********************************************************************
!subroutine for calculation of potential and derivatives of for HH
subroutine QHH_JHH(rHH, QHH, dQHHdr, &
     alphaHH, re, DHH, rcut)

  implicit none

  DOUBLE PRECISION,intent(in)  :: rHH,alphaHH, re, DHH, rcut
  DOUBLE PRECISION             :: expo2,rcutdiff
  DOUBLE PRECISION,intent(out) :: QHH, dQHHdr

  expo2   = dexp(alphaHH*(re-rHH))
  rcutdiff = dexp(alphaHH*(re-rcut))

  QHH  = DHH*expo2*(expo2-2.d0)-&
       DHH*rcutdiff*(rcutdiff-2.d0)+(rHH-rcut)*2.d0*alphaHH*DHH*rcutdiff*(rcutdiff-1.d0)

  dQHHdr = -2.d0*alphaHH*DHH*expo2*(expo2-1.d0)+2.d0*alphaHH*DHH*rcutdiff*(rcutdiff-1.d0)


end subroutine QHH_JHH



! **********************************************************************
subroutine QHNi_JHNi(rHNi,QHNi_HNi,dQHNidHNi,&
     alpha2, D2, r2,alpha2Ni, D2Ni, r2Ni)
  implicit none

  DOUBLE PRECISION,intent(in)  :: rHNi,alpha2, D2, r2,alpha2Ni, D2Ni, r2Ni
  DOUBLE PRECISION             :: expo2, sato,rcutdiff
  DOUBLE PRECISION :: DZHDR,DZNIDR
  DOUBLE PRECISION :: ZH, ZNI
  DOUBLE PRECISION :: ZHCOEF,ZNICOEF, A1, A2,C


  DOUBLE PRECISION,intent(out) :: QHNi_HNi,dQHNidHNi

  C = 14.3888D0 !eV*Angstrom

  ZHCOEF = D2*DEXP(-alpha2*rHNi)
  A1 = 1.D0+r2*rHNi
  ZH = ZHCOEF*A1
  !!
  ZNICOEF = D2Ni*DEXP(-alpha2Ni*rHNi)
  A2= 1.D0+r2Ni*rHNi
  ZNI = ZNICOEF*A2
  !!
  DZHDR = ZHCOEF*(r2-alpha2*A1)
  DZNIDR = ZNICOEF*(r2Ni-alpha2Ni*A2)
  !!
  QHNi_HNi = C*ZH*ZNI/rHNi
  !!
  dQHNidHNi = C*(DZHDR*ZNI + DZNIDR*ZH-ZNI*ZH/rHNi)/rHNi


  !  rcutdiff = dexp(alpha2*(r2-rcut))
  !  expo2  = dexp(alpha2*(r2-rHNI))
  !
  !
  !  QHNi_HNi  = D2*expo2*(expo2-2.d0)-&
  !             D2*rcutdiff*(rcutdiff-2.d0)+(rHNI-rcut)*2.d0*alpha2*D2*rcutdiff*(rcutdiff-1.d0)
  !
  !  dQHNidHNi = -2.d0*alpha2*D2*expo2*(expo2-1.d0)+2.d0*alpha2*D2*rcutdiff*(rcutdiff-1.d0)


end subroutine QHNi_JHNi

! **********************************************************************
subroutine QNiNi_JNiNi(rNiNi,PhiT,FT,alpha2Ni, D2Ni, r2Ni)

  IMPLICIT NONE

  DOUBLE PRECISION,intent(in)  :: alpha2Ni, D2Ni, r2Ni!,rcut

  DOUBLE PRECISION :: rNiNi, C,ZNI1
  DOUBLE PRECISION ::   ZNICOEF,A1,DZNIDR1

  !   DOUBLE PRECISION :: rNiNi, r0, D0, alpha
  !   DOUBLE PRECISION :: DT, rT, betaT,rcutdiff
  !   DOUBLE PRECISION :: rdiff, alpha2, betaTfac
  DOUBLE PRECISION,intent(out)  :: PhiT, FT
  !   DOUBLE PRECISION, dimension(14) :: epsilonT  !(epsilonT - linear termal expansion)
  !   INTEGER :: tmp

  !*************************************************
  ! Non-zero temperature Morse-Potential
  ! in molecular statics for Ni-crystals
  ! Computational Material Science 62 (2012) 126-130
  !*************************************************

  !   tmp   = 6 !1=5;2=25;3=50;4=100;5=200;6=293;7=400;8=500;9=600;10=700;11=800;12=900;13=1000;14=1100
  !   epsilonT = (/ -0.00231d0, -0.0023d0, -0.00229d0, 0.00208d0, -0.00114d0, 0.d0, 0.0015d0, 0.0029d0, 0.00455d0, 0.00617d0,&
  !               0.00783d0, 0.00953d0, 0.01126d0, 0.01302d0 /)
  !   alpha = alpha2Ni
  !   D0 = D2Ni
  !   r0 = r2Ni
  !
  !   alpha2   = 2.d0*alpha
  !   betaTfac = alpha2*epsilonT(tmp)*r0
  !
  !   betaT = dexp(-betaTfac)
  !   DT    = betaT * D0
  !   rT    = r0*(1.d0 + epsilonT(tmp))
  !   rcutdiff = rT - rcut
  !
  !   rdiff = rT - rNiNi
  !
  !   PhiT  = DT * (dexp(alpha2*(rdiff))-2.d0*dexp(alpha*(rdiff)))&
  !             -DT * (dexp(alpha2*(rcutdiff))-2.d0*dexp(alpha*(rcutdiff)))-&
  !             (rNINI-10.D0)*alpha2*DT * (dexp(alpha*(rcutdiff))-dexp(alpha2*(rcutdiff)))
  !
  !   FT    = alpha2*DT * (dexp(alpha*(rdiff))-dexp(alpha2*(rdiff)))&
  !           -alpha2*DT * (dexp(alpha*(rcutdiff))-dexp(alpha2*(rcutdiff))) !derivative of PhiT with repect to r

  C = 14.3888D0 !eV*Angstrom

  ZNICOEF = D2Ni*DEXP(-alpha2Ni*rNINI)
  A1=1.D0+r2Ni*rNINI
  ZNI1 = ZNICOEF*A1
  !
  DZNIDR1 = ZNICOEF*(r2Ni-alpha2Ni*A1)
  !!
  PhiT = C*ZNI1*ZNI1/rNINI
  FT = C*ZNI1/rNINI*(2.d0*DZNIDR1-ZNI1/rNINI)



END subroutine QNiNi_JNiNi



! **********************************************************************
!subroutine for calculation of potential and derivatives for HNi
subroutine QHNi_JHNi_RHO(RHO,QHNi_RHO,dQHNidrho,&
     alphaH, DH,rH,alphaHcor,DHcor,rHcor)

  implicit none

  DOUBLE PRECISION,intent(in)  :: RHO,alphaH, DH,rH,alphaHcor,DHcor,rHcor
  DOUBLE PRECISION             :: dEhdrho, drsdrho, EH
  DOUBLE PRECISION             :: alphaHP,alphaHPcut,alphaHPcor
  DOUBLE PRECISION             :: expo, expo0
  DOUBLE PRECISION,intent(out) :: QHNi_RHO,dQHNidrho


      call EMBEDING_HYDROGEN(RHO,QHNi_RHO,dQHNidrho)




!if(RHO .lt. 0.1D0) then
!QHNi_RHO=0.D0*rho
!dQHNidrho=0.D0
!
!  elseif((RHO .lt. 0.11D0).and.(RHO .ge. 0.1D0)) then

!      alphaHP = alphaH*(RHO-rh)
!      alphaHPcor = alphaHcor*(RHO-rHcor)
!
!
!      dQHNidrho = DH*(RHO-rH)*DEXP(-alphaHP)!+ DHcor*(RHO-rHcor)*DEXP(-alphaHPcor)
!      dQHNidrho = DH*(DEXP(-alphaHP)*(1.D0-alphaHP))!*DHcor*DEXP(-alphaHPcor)*(1.D0-alphaHPcor)


!else
!QHNi_RHO=0.D0*rho
!dQHNidrho=0.D0
!  endif

!  if(RHO <= 0.01D0) then

!  else !if(RHO <= 0.3D0) then

!
!  else
!      QHNi_RHO=1.D0
!      dQHNidrho=1.D0
!  endif



!
!  !      QHNi_RHO = DH*RHO*DEXP(-alphaHP)
!  !    dQHNidrho = DH*DEXP(-alphaHP)*(1.D0-alphaH*RHO)
!

  !
  !  expo   = dexp(-alphaH*(RHO-rH))
  !
  !Parts of exchange and coulomb interaction only depending on NiNi distance
  !  QHNi_RHO  = DH*expo*(expo-2.d0)
  !  dQHNidrho = alphaH*DH*expo*((-2.D0)*expo+2.d0)
  !      dQHNidrho=-6.95D0

  !print*,RHO,'RHO',QHNi_RHO,dQHNidrho

end subroutine QHNi_JHNi_RHO



! **********************************************************************
!subroutine for calculation of potential and derivatives of rho for NINI
subroutine QNiNi_JNiNi_RHO(rsd,QNiNi_RHO,dQNiNidrho,&
     alphaNi, DNi, rNi)

  implicit none

  DOUBLE PRECISION,intent(in)  :: rsd, alphaNi, DNi, rNi
  DOUBLE PRECISION             :: const,expo,sato,numberPI,expo0
  DOUBLE PRECISION,intent(out) :: QNiNi_RHO,dQNiNidrho


    call EMBEDING_NICKEL(rsd,QNiNi_RHO,dQNiNidrho)

!  expo   = dexp(alphaNi*(-rsd+rNi))
!  expo0   = dexp(alphaNi*(rNi))
!
!  QNiNi_RHO  = DNi*expo*(expo-2.d0)-DNi*expo0*(expo0-2.d0)
!  dQNiNidrho = alphaNi*DNi*expo*(-2.D0*expo+2.d0)

!expo=alphaNi*(rsd-rNi)

!QNiNi_RHO = DNi*(rsd-rNi)*DEXP(-expo)
!dQNiNidrho = DNi*DEXP(-expo)*(1.D0-expo)

end subroutine QNiNi_JNiNi_RHO



! **********************************************************************
!subroutine to calculate valence electron densty at one h atom
subroutine rho_hydrogen(r,p,dpdr)
  implicit none
  ! ------------------------------------------------------------------
  ! H Electron Density and derivatives for NI(100)/H EAM4 parametrization: REV B 51, 9985 (1995)
  ! ------------------------------------------------------------------
  double precision :: r,delh,ch,dpdr,p,DDPD2R
  double precision :: s, dsdr
  !
  delh = 3.779452267842503321304548080661334097d0 !2.d0/0.5291772083D0                    !1.30927d0   Tom Markland's values !2.53355563336653493078
  ch =  2.148061616605061452389691112330183387d0 !1.d0/(PI*0.5291772083D0**3.d0)          !11.0025d0
  !
  p = ch*dexp(-delh*r)
  dpdr = -delh*p
  !
  RETURN
END subroutine rho_hydrogen


!!  PARAMETERS FOR THE NI DENSITY BASED SPHERICALLY ADJUSTED S- AND D- LIKE ORBITALS
!! TAKEN FROM SINGLE DETERMINANT HARTREE-fOCK CALCULATIONS.
SUBROUTINE rho_Nickel(r,p,dpdr)
  IMPLICIT NONE
  !      DOUBLE PRECISION :: C(1:10)
  DOUBLE PRECISION :: EPS(1:10) !10, the total number of valence electrons in NI
  DOUBLE PRECISION :: NI(1:10)
  DOUBLE PRECISION :: R,P,DPDR,PS,PD,PPS,PPD,DPDRS,DPDRD
  INTEGER :: I
  DOUBLE PRECISION ::PI
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
  RETURN
END SUBROUTINE rho_Nickel
