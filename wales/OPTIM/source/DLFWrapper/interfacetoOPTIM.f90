! **********************************************************************
! **       Interface file to a modified version of dl_find            **
! **               by Dr. Judith B. Rommel (jbr36)                    **
! **             the potential and the keywords to do                 **
! **     quantum rate calculations are provided by OPTIM              **
! **                                                                  **
! **********************************************************************

module driver_module
  use dlf_parameter_module, only: rk

! variables for non-continuous differentiable MEP potential
  real(rk) :: smallbarvar !parameter to control the hight of the small barrier
end module driver_module


subroutine dlf_interfaceOPTIM(NATOMS)
! jbr36 DL-find modules
  use driver_module
  implicit none
  integer :: ivar,NATOMS

  call dlf_output(6,0)

  !nspec should be: nat + nz + 5*ncons + 2*nconn
  !nvarin should be 3* number of atoms
  !nvarin2 should be: nframe*nat*3 + nweight + nmass + n_po_scaling
  !master 1 if this task is the master of a parallel run, 0 otherwise (for parallel runs)
  !master 1 for serial runs
  !dl_find(nvarin,nvarin2,nspec,master)
  !nweight should be either 0 or nat
  !nmass should be either 0 or nat


  print*,'Calling the routines for instanton rate calculations in DL-find'
  call dl_find(3*natoms,4*natoms,2*natoms,1)

 ! Note on units:   DLFIND
 ! masses           atomic mass units
 ! coordinates      a.u. (internal) Angstrom (to read in and print out)
 ! temperature      Kelvin
 ! energie          hartree (atomic units) (to read in)

end subroutine dlf_interfaceOPTIM

! **********************************************************************
! subroutines that have to be provided to dl_find from outside
! **********************************************************************

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_get_params(nvar,nvar2,nspec,coords,coords2,spec,ierr, &
    tolerance,printl,maxcycle,maxene,tatoms,icoord, &
    iopt,iline,maxstep,scalestep,newlbfgs_mem,nimage,nebk,dump,restart,&
    nz,ncons,nconn,update,maxupd,delta,soft,inithessian,carthessian,tsrel, &
    maxrot,tolrot,nframe,nmass,nweight,timestep,fric0,fricfac,fricp, &
    imultistate, state_i,state_j,pf_c1,pf_c2,gp_c3,gp_c4,ln_t1,ln_t2, &
    printf,tolerance_e,distort,massweight,minstep,maxdump,task,temperature, &
    po_pop_size,po_radius,po_contraction,po_tolerance_r,po_tolerance_g, &
    po_distribution,po_maxcycle,po_init_pop_size,po_reset,po_mutation_rate, &
    po_death_rate,po_scalefac,po_nsave,ntasks,tdlf_farm,n_po_scaling, &
    neb_climb_test, neb_freeze_test, &
    nzero, coupled_states, qtsflag )
! jbr36 OPTIM modules
  use KEY, only: NFREEZE,INSTANTONSTARTDUMPT,INSTANTONOPTT, &
                  INSTANTONRATET,TSPLITTINGT,CLASSICALRATEST,VARSTEPOPTT,&
                  TEMPERATURE1,DISTORTINST,NIMAGEINST,DELTAINST
  use COMMONS, only: ATMASS
! jbr36 DL-find modules
  use dlf_parameter_module, only: rk
  use dlf_constants, only: dlf_constants_get
  use driver_module
  implicit none
  integer   ,intent(in)      :: nvar 
  integer   ,intent(in)      :: nvar2
  integer   ,intent(in)      :: nspec
  real(rk)  ,intent(inout)   :: coords(nvar) ! start coordinates
  real(rk)  ,intent(inout)   :: coords2(nvar2) ! a real array that can be used
                                               ! depending on the calculation
                                               ! e.g. a second set of coordinates
  integer   ,intent(inout)   :: spec(nspec)  ! specifications like fragment or frozen
  integer   ,intent(out)     :: ierr
  real(rk)  ,intent(inout)   :: tolerance
  real(rk)  ,intent(inout)   :: tolerance_e
  integer   ,intent(inout)   :: printl
  integer   ,intent(inout)   :: maxcycle
  integer   ,intent(inout)   :: maxene
  integer   ,intent(inout)   :: tatoms
  integer   ,intent(inout)   :: icoord
  integer   ,intent(inout)   :: iopt
  integer   ,intent(inout)   :: iline
  real(rk)  ,intent(inout)   :: maxstep
  real(rk)  ,intent(inout)   :: scalestep
  integer   ,intent(inout)   :: newlbfgs_mem
  integer   ,intent(inout)   :: nimage
  real(rk)  ,intent(inout)   :: nebk
  integer   ,intent(inout)   :: dump
  integer   ,intent(inout)   :: restart
  integer   ,intent(inout)   :: nz
  integer   ,intent(inout)   :: ncons
  integer   ,intent(inout)   :: nconn
  integer   ,intent(inout)   :: update
  integer   ,intent(inout)   :: maxupd
  real(rk)  ,intent(inout)   :: delta
  real(rk)  ,intent(inout)   :: soft
  integer   ,intent(inout)   :: inithessian
  integer   ,intent(inout)   :: carthessian
  integer   ,intent(inout)   :: tsrel
  integer   ,intent(inout)   :: maxrot
  real(rk)  ,intent(inout)   :: tolrot
  integer   ,intent(inout)   :: nframe
  integer   ,intent(inout)   :: nmass
  integer   ,intent(inout)   :: nweight
  real(rk)  ,intent(inout)   :: timestep
  real(rk)  ,intent(inout)   :: fric0
  real(rk)  ,intent(inout)   :: fricfac
  real(rk)  ,intent(inout)   :: fricp
  integer   ,intent(inout)   :: imultistate
  integer   ,intent(inout)   :: state_i
  integer   ,intent(inout)   :: state_j
  real(rk)  ,intent(inout)   :: pf_c1  
  real(rk)  ,intent(inout)   :: pf_c2  
  real(rk)  ,intent(inout)   :: gp_c3  
  real(rk)  ,intent(inout)   :: gp_c4
  real(rk)  ,intent(inout)   :: ln_t1  
  real(rk)  ,intent(inout)   :: ln_t2  
  integer   ,intent(inout)   :: printf
  real(rk)  ,intent(inout)   :: distort
  integer   ,intent(inout)   :: massweight
  real(rk)  ,intent(inout)   :: minstep
  integer   ,intent(inout)   :: maxdump
  integer   ,intent(inout)   :: task
  real(rk)  ,intent(inout)   :: temperature
  integer   ,intent(inout)   :: po_pop_size
  real(rk)  ,intent(inout)   :: po_radius
  real(rk)  ,intent(inout)   :: po_contraction
  real(rk)  ,intent(inout)   :: po_tolerance_r
  real(rk)  ,intent(inout)   :: po_tolerance_g
  integer   ,intent(inout)   :: po_distribution
  integer   ,intent(inout)   :: po_maxcycle
  integer   ,intent(inout)   :: po_init_pop_size
  integer   ,intent(inout)   :: po_reset
  real(rk)  ,intent(inout)   :: po_mutation_rate
  real(rk)  ,intent(inout)   :: po_death_rate
  real(rk)  ,intent(inout)   :: po_scalefac
  integer   ,intent(inout)   :: po_nsave
  integer   ,intent(inout)   :: ntasks
  integer   ,intent(inout)   :: tdlf_farm
  integer   ,intent(inout)   :: n_po_scaling
  real(rk)  ,intent(inout)   :: neb_climb_test
  real(rk)  ,intent(inout)   :: neb_freeze_test
  integer   ,intent(inout)   :: nzero
  integer   ,intent(inout)   :: coupled_states
  integer   ,intent(inout)   :: qtsflag
  ! local variables
  real(rk)                   :: svar
  integer                    :: iat,jat
  interface
    subroutine read_rand(arr)
      use dlf_parameter_module, only: rk
      use dlf_global, only : glob
      real(rk)  :: arr(:)
    end subroutine read_rand
  end interface
  integer  :: natoms !natoms read in from input.xyz
  character(LEN=5)  :: atomsymbol(nvar/3) !atomsymbol read in from input.xyz
  real(rk)  :: inputcoords(3,nvar/3) !coordinates read in from input.xyz
  integer  :: i1
  integer, external :: get_atom_charge
  real(rk) :: ang_au
! **********************************************************************
  ierr=0
  tsrel=1
  nweight=0

  print*,"External sizes:",nvar,nvar2,nspec

  open(501,file="input.xyz")
  call read_xyz(501,nvar/3,natoms,atomsymbol,inputcoords)
  close(501)
  if(nvar/=natoms*3) then
    write(*,'(A,A)') ' dlf_interface> Number of atoms in coords inconsistent with natoms in odata'
    STOP
  endif
  call dlf_constants_get("ANG_AU",ang_au)
  DO  I1 = 1,natoms
      coords(3*I1-2)=inputcoords(1,I1)/ang_au
      coords(3*I1-1)=inputcoords(2,I1)/ang_au
      coords(3*I1)=inputcoords(3,I1)/ang_au
  ENDDO


    nframe=1
    nzero=0!3*NFREEZE+6
    nmass=natoms
    coords2(1:nvar)=coords(:)
    coords2(nvar+1:nvar+natoms)=ATMASS(1:natoms) ! in a.m.u.
!
    spec(:)=0 ! all atoms active in internal coordinates
    spec(natoms-nfreeze+1:natoms)=-1 ! frozen atoms (the last nfreeze one)
    nz=natoms
    !spec(natoms+1)=1 !H
    !spec(natoms+2:natoms+nz)=28 !Ni
    !Adapt for other systems than nickel + H
    do I1=1,nz
        spec(natoms+i1)=get_atom_charge(atomsymbol(i1))
    enddo


!****************************************************
!*************END reading in coords ******************
! todo change maxcycle maxene tolerance

  tolerance=4.5D-11 ! negative: default settings
  printl=4
  printf=4
  maxcycle=100 !200
  maxene=100000

  tolrot=1.D2
!  tolrot=0.1D0
!  maxrot=100 !was 100

  task=0 !1011

  distort=DISTORTINST !0.D0 !0.4 
  tatoms=0
  icoord=190 !0 cartesian coord !210 Dimer !190 qts search !120 NEB frozen endpoint
  massweight=0
  iline=0
  maxstep=2.2D0
  scalestep=1.0D0
  newlbfgs_mem=100
  nimage=NIMAGEINST !*k-k+1
  temperature=TEMPERATURE1  !K

  ! Hessian
  delta=DELTAINST!1.D-2
  soft=-6.D-4
  update=2
  maxupd=0

  minstep=1.D0 **2 ! 1.D-5
  minstep=1.D-5

  nebk=0.D0
!****************** BEGIN of INSTANTON and classical rate OPTIONS********************
  qtsflag=0 ! 1: Tunneling splittings, 11/10: read image hessians
  if(INSTANTONSTARTDUMPT) then
   iopt=11
   inithessian=2
  elseif(CLASSICALRATEST) then
   iopt=13
   inithessian=2
  elseif(INSTANTONOPTT) then
    iopt=20
    inithessian=2
    if(TSPLITTINGT) qtsflag=1
    if(VARSTEPOPTT) nebk=1.0D0 ! for QTS calculations with variable tau,
                               ! nebk transports the parameter alpha (0=equidist. in tau)
  elseif(INSTANTONRATET) then
    iopt=12
    inithessian=5
    if(TSPLITTINGT) qtsflag=1
  endif


!***********************END of INSTANTON  adn classical rate OPTIONS********************
  ! damped dynamics
  fric0=0.1D0
  fricfac=1.0D0
  fricp=0.1D0

  dump=0
  restart=0

  ! Parallel optimization

  po_pop_size=25
  po_radius=0.5D0
  po_init_pop_size=50
  po_contraction=0.95D0
  po_tolerance_r=1.0D-8
  po_tolerance_g=1.0D-6
  po_distribution=3
  po_maxcycle=100000
  po_reset=500
  po_mutation_rate=0.15D0
  po_death_rate=0.5D0
  po_scalefac=10.0D0
  po_nsave=10
  n_po_scaling=0 ! meaning that the base radii values for the sampling and tolerance 
                 ! are used for all components of the working coordinate vector.
                 ! Remember to change the second arg in the call to dl_find if a 
                 ! non-zero value for n_po_scaling is desired, and also to add the 
                 ! necessary values to the coords2 array...

  ! Taskfarming 
  ! (the two lines below, with any values assigned, 
  ! may be safely left in place for a serial build) 
  ntasks = 1
  tdlf_farm = 1

  tatoms=1
  
!  call test_ene
   imultistate=0
  
end subroutine dlf_get_params

subroutine test_update
  use dlf_parameter_module, only: rk
  use dlf_allocate, only: allocate,deallocate
  use dlf_global, only: glob,printl
  implicit none
  integer(4) :: varperimage,nimage,iimage,ivar4
  real(rk), allocatable:: coords(:,:),grad(:,:) ! varperimage,nimage
  real(rk), allocatable:: hess(:,:,:) ! varperimage,varperimage,nimage
  real(rk), allocatable:: fhess(:,:) ! 2*varperimage*nimage,2*varperimage*nimage
  real(rk), allocatable:: eigval(:),eigvec(:,:)
  real(rk), allocatable:: vec0(:)
  real(rk), allocatable:: determinant(:)
  real(rk), allocatable:: tmphess(:,:)
  real(rk) :: svar
  integer :: ivar,jvar,vp8,target_image,step,lastimage,jimage,turnimg
  logical :: havehessian,fracrecalc

  open(unit=102,file="grad_coor.bin",form="unformatted")
  read(102) varperimage,nimage
  print*,"varperimage,nimage",varperimage,nimage
  vp8=varperimage

!!$  call allocate(coords,int(varperimage,kind=8),int(nimage,kind=8))
!!$  call allocate(grad,int(varperimage,kind=8),int(nimage,kind=8))
!!$  call allocate(determinant,int(nimage,kind=8))

  do iimage=1,nimage
    read(102) ivar4
    print*,"Reading coords/grad of image",ivar4
    read(102) grad(:,iimage)
    read(102) coords(:,iimage)
  end do
  close(102)
  print*,"Coords sucessfully read"

  ! print coords and grad
  iimage=2
  print*,"Coords and grad for image ",iimage
  do ivar=1,vp8
    write(6,"(i6,1x,2es18.9)") &
        ivar,coords(ivar,iimage),grad(ivar,iimage)
  end do

  open(unit=101,file="hessian.bin",form="unformatted")
  ! Hessian in mass-weighted coordinates on the diagnal blocks - everything else should be zero
  read(101) iimage,ivar4!neb%varperimage,neb%nimage
  if(iimage/=varperimage.or.ivar4/=nimage) then
    print*,"Dimensions read",iimage,ivar4
    call dlf_fail("ERROR: wrong dimensions in hessian.bin!")
  end if
  ivar=2*varperimage*nimage
  print*,"File Hessian size",ivar
  call allocate(fhess,ivar,ivar)
  read(101) fhess
  close(101)
  print*,"Hessian sucessfully read"

  ! map hessian to different array:
!!$  call allocate(hess,int(varperimage,kind=8),int(varperimage,kind=8),int(nimage,kind=8))
  do iimage=1,nimage
    print*,"Image",iimage,"Hessian positions",(iimage-1)*varperimage+1,iimage*varperimage
    hess(:,:,iimage)=fhess((iimage-1)*varperimage+1:iimage*varperimage,(iimage-1)*varperimage+1:iimage*varperimage)
  end do
  call deallocate(fhess)

  call allocate(eigval,vp8)
  call allocate(eigvec,vp8,vp8)
  ! now we have all we need

  print*,"# Distance from previous image"
  do iimage=2,nimage
    print*,iimage,sqrt(sum( (coords(:,iimage)-coords(:,iimage-1))**2))
  end do

  print*,"# Determinant of Hessian"
  do iimage=1,nimage
    do ivar=1,vp8
      do jvar=ivar+1,vp8
        if(abs(hess(ivar,jvar,iimage)-hess(jvar,ivar,iimage))>1.D-20) &
            print*,"Unsymmetric:",ivar,jvar,iimage,hess(ivar,jvar,iimage),hess(jvar,ivar,iimage)
      end do
    end do
    call dlf_matrix_diagonalise(vp8,hess(:,:,iimage),eigval,eigvec)
    !write(6,"(i6,1x,9es12.3)") iimage,product(eigval(8:vp8)),eigval(1:8)
    do ivar=1,6
      eigval(minloc(abs(eigval)))=1.D0
    end do
    determinant(iimage)=product(eigval)
    write(6,"(i6,1x,9es12.3)") iimage,product(eigval),eigval(1:8)
  end do

  print*,"maxval(hess(:,:,nimage-1)-hess(:,:,nimage))",maxval(hess(:,:,nimage-1)-hess(:,:,nimage))

!!$  ! Richtungsableitung
!!$  call allocate(vec0,vp8)!int(varperimage,kind=8))
!!$  do iimage=2,nimage-1
!!$    ! eigval is Vector along which the derivative is taken
!!$    eigval=coords(:,iimage+1)-coords(:,iimage-1)
!!$    svar=sqrt(sum(eigval**2))
!!$    eigval=eigval/svar
!!$    vec0=matmul(hess(:,:,iimage),eigval)
!!$    do ivar=1,vp8
!!$      write(6,"(2i6,1x,2es18.9,1x,f10.5)") &
!!$          iimage,ivar,(grad(ivar,iimage+1)-grad(ivar,iimage-1))/svar,vec0(ivar),&
!!$          vec0(ivar)/((grad(ivar,iimage+1)-grad(ivar,iimage-1))/svar)
!!$    end do
!!$  end do

  !
  ! now test updates
  !
  call allocate(tmphess,vp8,vp8)
  havehessian=.true.
  fracrecalc=.false.
  printl=2
  glob%maxupd=30000

  target_image=1 !nimage
  ! update hessians to the one of the first image
  print*,"Updating Hessians to that of image",target_image
  print*,"Sum-of-squares difference"
  do iimage=1,nimage
    tmphess(:,:)=hess(:,:,iimage)
    call dlf_hessian_update(vp8, &
        coords(:,target_image),coords(:,iimage),&
        grad(:,target_image),grad(:,iimage), &
        tmphess, havehessian, fracrecalc)
    if(.not.havehessian) then
      print*,"Problem with hessian update, image",iimage
      havehessian=.true.
    end if
    print*,iimage,sum( (tmphess-hess(:,:,target_image))**2),&
        sum( (hess(:,:,iimage)-hess(:,:,target_image))**2)
  end do

  print*,"Minstep",glob%minstep
  print*,"Updating Hessians to that of image",target_image
  print*,"Determinant"
  open(file="determinant",unit=10)
  do iimage=1,nimage
    tmphess(:,:)=hess(:,:,iimage)
    call dlf_hessian_update(vp8, &
        coords(:,target_image),coords(:,iimage),&
        grad(:,target_image),grad(:,iimage), &
        tmphess, havehessian, fracrecalc)
    if(.not.havehessian) then
      print*,"Problem with hessian update, image",iimage
      havehessian=.true.
    end if
    call dlf_matrix_diagonalise(vp8,tmphess,eigval,eigvec)
    !write(6,"(i6,1x,9es12.3)") iimage,product(eigval(8:vp8)),eigval(1:8)
    do ivar=1,6
      eigval(minloc(abs(eigval)))=1.D0
    end do

    print*,iimage,product(eigval),determinant(iimage),determinant(target_image)
    write(10,*) iimage,product(eigval),determinant(iimage),determinant(target_image)
  end do
  close(10)

  print*,"Updating Hessians to that of image",target_image
  print*,"Determinant - incremental"
  open(file="determinant_incr",unit=10)
  do iimage=1,nimage
    tmphess(:,:)=hess(:,:,iimage)
    step=1
    if(iimage>target_image) step=-1
    lastimage=iimage
    do jimage=iimage+step,target_image,step
      !print*,"updating",lastimage," to ",jimage
      call dlf_hessian_update(vp8, &
          coords(:,jimage),coords(:,lastimage),&
          grad(:,jimage),grad(:,lastimage), &
          tmphess, havehessian, fracrecalc)

      if(.not.havehessian) then
        print*,"Problem with hessian update, image",iimage
        havehessian=.true.
      end if
      lastimage=jimage
    end do
    call dlf_matrix_diagonalise(vp8,tmphess,eigval,eigvec)
    !write(6,"(i6,1x,9es12.3)") iimage,product(eigval(8:vp8)),eigval(1:8)
    do ivar=1,6
      eigval(minloc(abs(eigval)))=1.D0
    end do

    print*,iimage,product(eigval),determinant(iimage),determinant(target_image)
    write(10,*) iimage,product(eigval),determinant(iimage),determinant(target_image)
  end do
  close(10)

  print*,"Updating Hessians to that of image",target_image
  print*,"Determinant - incremental turning around"
  turnimg=20
  open(file="determinant_turn",unit=10)
  do iimage=1,nimage
    tmphess(:,:)=hess(:,:,iimage)
    lastimage=iimage
    ! first upwards to turnimg
    do jimage=iimage+1,turnimg
      print*,"updating",lastimage," to ",jimage
      call dlf_hessian_update(vp8, &
          coords(:,jimage),coords(:,lastimage),&
          grad(:,jimage),grad(:,lastimage), &
          tmphess, havehessian, fracrecalc)

      if(.not.havehessian) then
        print*,"Problem with hessian update, image",iimage
        havehessian=.true.
      end if
      lastimage=jimage
    end do
    step=1
    if(lastimage>target_image) step=-1
    do jimage=lastimage+step,target_image,step
      print*,"updating",lastimage," to ",jimage
      call dlf_hessian_update(vp8, &
          coords(:,jimage),coords(:,lastimage),&
          grad(:,jimage),grad(:,lastimage), &
          tmphess, havehessian, fracrecalc)

      if(.not.havehessian) then
        print*,"Problem with hessian update, image",iimage
        havehessian=.true.
      end if
      lastimage=jimage
    end do
    call dlf_matrix_diagonalise(vp8,tmphess,eigval,eigvec)
    !write(6,"(i6,1x,9es12.3)") iimage,product(eigval(8:vp8)),eigval(1:8)
    do ivar=1,6
      eigval(minloc(abs(eigval)))=1.D0
    end do

    print*,iimage,product(eigval),determinant(iimage),determinant(target_image)
    write(10,*) iimage,product(eigval),determinant(iimage),determinant(target_image)
  end do
  close(10)

!!$  do ivar=1,vp8
!!$    WRITE(*,'(33f10.5)') hess(ivar,:,1)*1.D6
!!$  end do
!!$
!!$  do ivar=1,vp8
!!$    write(6,"(i6,1x,es18.9)") &
!!$            ivar,eigval(ivar)
!!$  end do

  call deallocate(coords)
  call deallocate(grad)

  call dlf_fail("stop in test_update")
end subroutine test_update

! just a test routine
subroutine test_ene
  use dlf_parameter_module, only: rk
  implicit none
  integer :: ivar,ivar2,status
  integer :: samples
  real(rk) :: halfSamples
  real(rk) :: coords(3),grad(3),hess(3,3),ene
  coords(:)=0.D0
!  open(file="energy",unit=13)
  open(file="energy-2d.dat",unit=13)
  samples = 100
  halfSamples = dble(samples) * 0.25D0
  do ivar2=1,samples
    do ivar=1,samples
      coords(1)=dble(ivar-samples/2)/halfSamples
      coords(2)=dble(ivar2-samples/2)/halfSamples
      call dlf_get_gradient(3,coords,ene,grad,1,status)
      call dlf_get_hessian(3,coords,hess,status)
     write(13,*) coords(1),coords(2),ene,grad(1),grad(2),hess(1,1),hess(2,2),hess(2,1)
!     write(13,*) coords(1),coords(2),ene,grad(1)
    end do
    write(13,*) ""
  end do
  close(13)
  call dlf_fail("stop in test_ene")
end subroutine test_ene

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_get_gradient(nvar,coords,energy,gradient,iimage,status)
  use key, only: TTM3T
  use dlf_parameter_module, only: rk
  use dlf_constants, only: dlf_constants_get
  use driver_module
  implicit none
  integer   ,intent(in)    :: nvar
!  real(rk)  ,intent(in)    :: coords(nvar)
  real(rk)                 :: coords(nvar)
  real(rk)  ,intent(out)   :: energy
  real(rk)  ,intent(out)   :: gradient(nvar)
  integer   ,intent(in)    :: iimage
  integer   ,intent(out)   :: status

  real(rk) :: RMS

! **********************************************************************
!  call test_update
  status=1
  if(ttm3t) then
    call atomic_unitskcal(coords(:),gradient(:),energy,nvar,-1) !conversion: Hartree (a.u.) to kcal/Mol Bohr to Angstrom

    call  potential(coords,energy,gradient,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)

    call atomic_unitskcal(coords(:),gradient(:),energy,nvar,+1) !conversion from kcal/mol Hartree (a.u.)

  else
    call atomic_units(coords(:),gradient(:),energy,nvar,-1) !conversion: Hartree (a.u.) to eV, Bohr to Angstrom

    call  potential(coords,energy,gradient,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)

    call atomic_units(coords(:),gradient(:),energy,nvar,+1) !conversion from eV to Hartree (a.u.)
  endif

  status=0
end subroutine dlf_get_gradient




! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_get_hessian(nvar,coords,hessian,status)
  !  get the hessian at a given geometry
  use dlf_parameter_module
  use driver_module
  implicit none
  integer   ,intent(in)    :: nvar
  real(rk)  ,intent(in)    :: coords(nvar)
  real(rk)  ,intent(out)   :: hessian(nvar,nvar)
  integer   ,intent(out)   :: status

! **********************************************************************
  hessian(:,:)=0.D0
  status=1

  print*, 'WARNING: No external Hessian'

   status=0
end subroutine dlf_get_hessian




! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_put_coords(nvar,mode,energy,coords,iam)
  use dlf_parameter_module
  implicit none
  integer   ,intent(in)    :: nvar
  integer   ,intent(in)    :: mode
  integer   ,intent(in)    :: iam
  real(rk)  ,intent(in)    :: energy
  real(rk)  ,intent(in)    :: coords(nvar)
  integer                  :: iat
! **********************************************************************

! Only do this writing of files if I am the rank-zero processor
  if (iam /= 0) return

  if(mod(nvar,3)==0) then
    !assume coords are atoms
    if(mode==2) then
      open(unit=20,file="tsmode.xyz")
    else
      open(unit=20,file="coords.xyz")
    end if
    write(20,*) nvar/3
    write(20,*) 
    do iat=1,nvar/3
      write(20,'("H ",3f12.7)') coords((iat-1)*3+1:(iat-1)*3+3)
    end do
    close(20)
  else
    !print*,"Coords in put_coords: ",coords
  end if
end subroutine dlf_put_coords

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_error()
  implicit none
! **********************************************************************
  call dlf_mpi_abort() ! only necessary for a parallel build;
                       ! can be present for a serial build
  STOP !call exit(1)
end subroutine dlf_error

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_update()
  implicit none
! **********************************************************************
  ! only a dummy routine here.
end subroutine dlf_update


subroutine dlf_get_multistate_gradients(nvar,coords,energy,gradient,&
    mscoupling,needcoupling,iimage,status)
  ! only a dummy routine up to now
  ! for conical intersection search
  use dlf_parameter_module
  implicit none
  integer   ,intent(in)    :: nvar
  integer   ,intent(in)    :: coords(nvar)
  real(rk)  ,intent(in)    :: energy(2)
  real(rk)  ,intent(in)    :: gradient(nvar,2)
  integer   ,intent(in)    :: iimage
  integer   ,intent(in)    :: status
  real(rk) :: mscoupling(:,:)
  integer  :: needcoupling
end subroutine dlf_get_multistate_gradients


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_put_procinfo(dlf_nprocs, dlf_iam, dlf_global_comm)

  implicit none

  integer, intent(in) :: dlf_nprocs ! total number of processors
  integer, intent(in) :: dlf_iam ! my rank, from 0, in mpi_comm_world
  integer, intent(in) :: dlf_global_comm ! world-wide communicator
! **********************************************************************

!!! variable in the calling program = corresponding dummy argument

end subroutine dlf_put_procinfo


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_get_procinfo(dlf_nprocs, dlf_iam, dlf_global_comm)

  implicit none

  integer :: dlf_nprocs ! total number of processors
  integer :: dlf_iam ! my rank, from 0, in mpi_comm_world
  integer :: dlf_global_comm ! world-wide communicator
! **********************************************************************

!!! dummy argument = corresponding variable in the calling program

end subroutine dlf_get_procinfo


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_put_taskfarm(dlf_ntasks, dlf_nprocs_per_task, dlf_iam_in_task, &
                        dlf_mytask, dlf_task_comm, dlf_ax_tasks_comm)

  implicit none

  integer, intent(in) :: dlf_ntasks          ! number of taskfarms
  integer, intent(in) :: dlf_nprocs_per_task ! no of procs per farm
  integer, intent(in) :: dlf_iam_in_task     ! my rank, from 0, in my farm
  integer, intent(in) :: dlf_mytask          ! rank of my farm, from 0
  integer, intent(in) :: dlf_task_comm       ! communicator within each farm
  integer, intent(in) :: dlf_ax_tasks_comm   ! communicator involving the 
                                             ! i-th proc from each farm
! **********************************************************************

!!! variable in the calling program = corresponding dummy argument 

end subroutine dlf_put_taskfarm


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_get_taskfarm(dlf_ntasks, dlf_nprocs_per_task, dlf_iam_in_task, &
                        dlf_mytask, dlf_task_comm, dlf_ax_tasks_comm)

  implicit none

  integer :: dlf_ntasks          ! number of taskfarms
  integer :: dlf_nprocs_per_task ! no of procs per farm
  integer :: dlf_iam_in_task     ! my rank, from 0, in my farm
  integer :: dlf_mytask          ! rank of my farm, from 0
  integer :: dlf_task_comm       ! communicator within each farm
  integer :: dlf_ax_tasks_comm   ! communicator involving the
                                 ! i-th proc from each farm
! **********************************************************************

!!! dummy argument = corresponding variable in the calling program

end subroutine dlf_get_taskfarm


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_output(dum_stdout, dum_stderr)
  use dlf_parameter_module, only: rk
  use dlf_global, only: glob,stderr,stdout,keep_alloutput
  implicit none
  integer :: dum_stdout
  integer :: dum_stderr
  integer :: ierr
  logical :: topened
  character(len=10) :: suffix

! sort out output units; particularly important on multiple processors
 
! set unit numbers for main output and error messages
  if (dum_stdout >= 0) stdout = dum_stdout 
  if (dum_stderr >= 0) stderr = dum_stderr

  if (glob%iam /= 0) then
     inquire(unit=stdout, opened=topened, iostat=ierr)
     if (topened .and. ierr == 0) close(stdout)
     if (keep_alloutput) then ! hardwired in dlf_global_module.f90
        write(suffix,'(i10)') glob%iam
        open(unit=stdout,file='output.proc'//trim(adjustl(suffix)))
     else
        open(unit=stdout,file='/dev/null')
     end if
  endif

  if (glob%nprocs > 1) then
     ! write some info on the parallelization
     write(stdout,'(1x,a,i10,a)')"I have rank ",glob%iam," in mpi_comm_world"
     write(stdout,'(1x,a,i10)')"Total number of processors = ",glob%nprocs
     if (keep_alloutput) then
        write(stdout,'(1x,a)')"Keeping output from all processors"
     else
        write(stdout,'(1x,a)')"Not keeping output from processors /= 0"
     end if
  end if

end subroutine dlf_output


! **********************************************************************
! **********************************************************************
! The following routine either writes random numbers to a file, or reads
! them. This is to have equal starting conditions for different compilers
subroutine read_rand(arr)
  use dlf_parameter_module, only: rk
  use dlf_global, only : glob
  real(rk)  :: arr(:)
  integer, parameter :: si1=12
  integer, parameter :: si2=3000
  logical,parameter :: readf=.true.
  real(rk) :: ar(si2)
  integer :: l(1),length
  l=ubound(arr)
  length=l(1)
  if(readf) then
    if(length<=si1) then
      open(unit=201,file="random1.bin",form="unformatted")
    else if(length<=si2) then
      open(unit=201,file="random2.bin",form="unformatted")
    else
      call dlf_mpi_finalize() ! only necessary for a parallel build;
                              ! can be present for a serial build
      stop "Too many coordinates to be read from random.bin file"
    end if
    read(201) ar(1:length)
    close(201)
    arr=ar(1:length)
  else
    if (glob%iam == 0) then
       call random_number(ar)
       open(unit=201,file="random1.bin",form="unformatted")
       write(201) ar(1:si1)
       close(201)
       open(unit=201,file="random2.bin",form="unformatted")
       write(201) ar(1:si2)
       close(201)
    end if
  end if
end subroutine read_rand

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     read_xyz accepts
!     coords.xyz input files (cartesian angstrom)
!         header           : # of atoms
!         with atom records: tag x y z

subroutine read_xyz(unit,nat,natoms,atomsymbol,coords)
  use dlf_parameter_module, only: rk
  use porfuncs

  implicit none
  integer,intent(in) :: unit
  integer,intent(in) :: nat
  integer,intent(out) :: natoms
  character(LEN=5),intent(out) :: atomsymbol(nat)
  real(rk),intent(out):: coords(3,nat)
  integer            :: iat

! **********************************************************************
  read(unit,*) natoms
  read(unit,*)
  do iat=1,natoms
    read(unit,*) atomsymbol(iat),coords(1,iat),coords(2,iat),coords(3,iat)
  end do
  call flush(unit)
end subroutine read_xyz


! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
integer function get_atom_charge(atomsymbol)
  implicit none
  character(2), intent(in) :: atomsymbol
  integer iat
  character(2), parameter :: elements(111) = &
       (/ 'H ','He', &
          'Li','Be','B ','C ','N ','O ','F ','Ne', &
          'Na','Mg','Al','Si','P ','S ','Cl','Ar', &
          'K ','Ca','Sc','Ti','V ','Cr','Mn','Fe','Co','Ni','Cu', &
          'Zn','Ga','Ge','As','Se','Br','Kr', &
          'Rb','Sr','Y ','Zr','Nb','Mo','Tc','Ru','Rh','Pd','Ag', &
          'Cd','In','Sn','Sb','Te','I ','Xe', &
          'Cs','Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy', &
          'Ho','Er','Tm','Yb','Lu','Hf','Ta','W ','Re','Os','Ir','Pt', &
          'Au','Hg','Tl','Pb','Bi','Po','At','Rn', &
          'Fr','Ra','Ac','Th','Pa','U ','Np','Pu','Am','Cm','Bk','Cf', &
          'Es','Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','Ds', &
          'Rg' /)
! **********************************************************************
  get_atom_charge=1
  do iat=1,111
    if(elements(iat).eq.atomsymbol) then
        get_atom_charge=iat
        return
    endif
  enddo
end function get_atom_charge


  subroutine atomic_units(xcoord,dxcoord,energyconv,nvar,iflag)
  ! Input/output geometry in Bohr (a.u.), energies in Hartree (a.u.), and gradients in Hartree/Bohr.
  ! Internal geometry in Angstrom, energies in eV, and gradients in eV/Angstrom!
   implicit none
   integer, intent(in) :: nvar
   double precision ::xcoord(nvar),dxcoord(nvar)
   double precision :: energyconv
   double precision :: ang_au, ev_hartree, factor2
   integer :: iflag  ! -1 - coordinates to Angs, +1 - coordinates, energy and forces to a.u
!Constants taken from NIST
   ang_au=5.2917720810086E-01
   ev_hartree=3.67493237981D-2

   if(iflag.eq.-1) then !conversion: Hartree (a.u.) to eV, Bohr to Angstrom
!
   xcoord(:)=xcoord(:)*ang_au
!
   else if(iflag.eq.1) then !conversion: eV to Hartree, Angstrom to Bohr
    xcoord(:)=xcoord(:)/ang_au
    energyconv=energyconv*ev_hartree!conversion from eV to Hartree (a.u.)
    dxcoord(:)=dxcoord(:)*ev_hartree*ang_au ! conversion from eV/Angstrom to Hartree/Bohr
   else
   write(6,*) '* Convert units - iflag invalid'
   stop
   endif
!
   return
   end

     subroutine atomic_unitskcal(xcoord,dxcoord,energyconv,nvar,iflag)
  ! Input/output geometry in Bohr (a.u.), energies in Hartree (a.u.), and gradients in Hartree/Bohr.
  ! Internal geometry in Angstrom, energies in kcal/mol, and gradients in kcal/mol/Angstrom!
   implicit none
   integer, intent(in) :: nvar
   double precision ::xcoord(nvar),dxcoord(nvar)
   double precision :: energyconv
   double precision :: ang_au, kcalmol_hartree, factor2
   integer :: iflag  ! -1 - coordinates to Angs, +1 - coordinates, energy and forces to a.u
!Constants taken from NIST
   ang_au=5.2917720810086E-01
   kcalmol_hartree=1.5936D-3
   if(iflag.eq.-1) then !conversion: Hartree (a.u.) to kcal, Bohr to Angstrom
!
   xcoord(:)=xcoord(:)*ang_au
!
   else if(iflag.eq.1) then !conversion: kcal to Hartree, Angstrom to Bohr
    xcoord(:)=xcoord(:)/ang_au
    energyconv=energyconv*kcalmol_hartree!conversion from kcal to Hartree (a.u.)
    dxcoord(:)=dxcoord(:)*kcalmol_hartree*ang_au ! conversion from kcal/Angstrom to Hartree/Bohr
   else
   write(6,*) '* Convert units - iflag invalid'
   stop
   endif
!
   return
   end
