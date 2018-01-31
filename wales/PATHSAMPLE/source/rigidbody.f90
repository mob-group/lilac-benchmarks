!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999-2005 David J. Wales
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

module rigidbody

  use mathsconstants

  implicit none

! #####################################################################################################

! Define a new structure type
! This will make it easier if we want to have more than one type of rigid body
! in the same system

  type rigidbodyPotential

     integer :: nSites
     integer :: nPhysicalSites        ! For the purposes of calculating path lengths etc

! Reference geometry
     real (kind=kind(0.0d0)), allocatable :: SITE(:,:)

! Site symbols for output

     character(5), allocatable :: siteLabel(:)

! Masses for calculating rotational constants etc

     real (kind=kind(0.0d0)), allocatable :: mass(:)

! Atomic numbers - not sure why we need these and only valid for atomic models anyway

     integer, allocatable :: atomicNumber(:)

! Partial charges
     real (kind=kind(0.0d0)), allocatable :: CHARGE(:)

! LJ parameters
     real (kind=kind(0.0d0)), allocatable :: C6(:)
     real (kind=kind(0.0d0)), allocatable :: C12(:)

! Morse parameters - epsilon, rho, Re

     real (kind=kind(0.0d0)), allocatable :: morseParameters(:,:)    

! Conversion factor for normal mode frequencies

     real (kind=kind(0.0d0)) :: freqConvFactor = 0.0d0

! Conversion factor for Coulomb energies 

     real (kind=kind(0.0d0)) :: coulombConversionFactor = 1.0d0

! Conversion factor for electric field interactions

     real (kind=kind(0.0d0)) :: efieldConversionFactor = 1.0d0     

  end type rigidbodyPotential

! #####################################################################################################

! Store the definition of a virus capsomer

  type capsomer
     integer :: nBasalSites
     real (kind=kind(0.0d0)) :: rho, radius, height
     real (kind=kind(0.0d0)) :: epsilonRep   ! Repulsive as a fraction of epsilon_Morse
  end type capsomer

! #####################################################################################################

  save

  integer :: numRBtypes
  type (capsomer), allocatable :: capsomerDefs(:)

  type(rigidbodyPotential) :: rbPotential
  logical, allocatable :: sitesInteract(:,:,:)

! Small angle below which approximations to trigonometric functions are employed

  real (kind=kind(0.0d0)), parameter :: smallAngle = 1D-5 
  
contains

! #####################################################################################################

! Initialise known rigid body models

  subroutine initialiseRigidBody (atomType, nbodies, systemMasses, systemAtNumbers)

    implicit none
    character(5), intent(IN) :: atomType
    integer, intent(IN), optional :: nbodies
    real (kind=kind(0.0d0)), pointer, dimension(:), optional :: systemMasses
    integer, pointer, dimension(:), optional :: systemAtNumbers

    integer :: currentBody, currentSite
    real (kind=kind(0.0d0)) :: currentAngle

    select case (atomType(1:1))

! TIPnP family of potentials
    case ('W')
       rbPotential%nPhysicalSites = 3

       allocate(rbPotential%siteLabel(rbPotential%nPhysicalSites))
       allocate(rbPotential%mass(rbPotential%nPhysicalSites))
       allocate(rbPotential%atomicNumber(rbPotential%nPhysicalSites))

       rbPotential%siteLabel = (/'O', 'H', 'H'/)
       rbPotential%mass = (/16.0d0, 1.0d0, 1.0d0/)       ! Pure H2O, not isotopic abundance weighted
       rbPotential%atomicNumber = (/8, 1, 1/)

       rbPotential%freqConvFactor = 1d3/(2.99792458d0*2.0d0*pi)    ! Converts frequencies to wavenumbers
                                                                   ! See journal reference 01/07/05

       rbPotential%coulombConversionFactor = 1389.354848D0
       rbPotential%efieldConversionFactor = 96.4847D0 ! eV to kJ/mol conversion, input field in V/A

! Initialise known system types

       select case (atomType(2:2))
       case ('1')
          rbPotential%nSites = 3
          allocate(rbPotential%site(rbPotential%nSites,3),rbPotential%CHARGE(rbPotential%nSites))
          allocate(rbPotential%C6(rbPotential%nSites),rbPotential%C12(rbPotential%nSites))
          
          rbPotential%site(1,:) = 0.0d0
          rbPotential%site(2,:) = (/0.2785156560837814d0, -0.6990486725944498d0, 0.5916010671560346d0/)
          rbPotential%site(3,:) = (/-0.0640734628674161d0, -0.4219928632410997D0, -0.8567662777734411D0/)

          rbPotential%CHARGE = (/-0.8D0, 0.4D0, 0.4D0/)

! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
          rbPotential%C6 = (/2510.4D0, 0.0d0, 0.0d0/) 
          rbPotential%C12 = (/2426720.0D0, 0.0d0, 0.0d0/)

       case ('2')
          rbPotential%nSites = 4
          allocate(rbPotential%site(rbPotential%nSites,3),rbPotential%CHARGE(rbPotential%nSites))
          allocate(rbPotential%C6(rbPotential%nSites),rbPotential%C12(rbPotential%nSites))
          
          rbPotential%site(1,:) = 0.0d0
          rbPotential%site(2,:) = (/0.2785156560837814d0, -0.6990486725944498d0, 0.5916010671560346d0/)
          rbPotential%site(3,:) = (/-0.0640734628674161d0, -0.4219928632410997D0, -0.8567662777734411D0/)
          rbPotential%site(4,:) = (/0.0274511879486426D0, -0.1435068418061116D0, -0.0339443461425310D0/)

          rbPotential%CHARGE = (/0.0d0, 0.535D0, 0.535D0, -1.07D0/)

! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
          rbPotential%C6 = (/2510.4D0, 0.0d0, 0.0d0, 0.0d0/) 
          rbPotential%C12 = (/2907880.0D0, 0.0d0, 0.0d0, 0.0d0/)

       case ('3')
          rbPotential%nSites = 3
          allocate(rbPotential%site(rbPotential%nSites,3),rbPotential%CHARGE(rbPotential%nSites))
          allocate(rbPotential%C6(rbPotential%nSites),rbPotential%C12(rbPotential%nSites))
          
          rbPotential%site(1,:) = 0.0d0
          rbPotential%site(2,:) = (/0.2785156560837814d0, -0.6990486725944498d0, 0.5916010671560346d0/)
          rbPotential%site(3,:) = (/-0.0640734628674161d0, -0.4219928632410997D0, -0.8567662777734411D0/)

          rbPotential%CHARGE = (/-0.834D0, 0.417D0, 0.417D0/)

! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
          rbPotential%C6 = (/2489.48D0, 0.0d0, 0.0d0/)             ! 2489.71000D0 
          rbPotential%C12 = (/2435088.0D0, 0.0d0, 0.0d0/)          ! 2435099.13639D0

       case ('4')
          rbPotential%nSites = 4
          allocate(rbPotential%site(rbPotential%nSites,3),rbPotential%CHARGE(rbPotential%nSites))
          allocate(rbPotential%C6(rbPotential%nSites),rbPotential%C12(rbPotential%nSites))
          
          rbPotential%site(1,:) = 0.0d0
          rbPotential%site(2,:) = (/0.2785156560837814d0, -0.6990486725944498d0, 0.5916010671560346d0/)
          rbPotential%site(3,:) = (/-0.0640734628674161d0, -0.4219928632410997D0, -0.8567662777734411D0/)
          rbPotential%site(4,:) = (/0.0274511879486426D0, -0.1435068418061116D0, -0.0339443461425310D0/)

          rbPotential%CHARGE = (/0.0d0, 0.52D0, 0.52D0, -1.04D0/)

! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
          rbPotential%C6 = (/2552.24D0, 0.0d0, 0.0d0, 0.0d0/)       ! 2551.90393D0
          rbPotential%C12 = (/2510.4D3, 0.0d0, 0.0d0, 0.0d0/)       ! 2510413.58406D0      

       case ('5')
          rbPotential%nSites = 5
          allocate(rbPotential%site(rbPotential%nSites,3),rbPotential%CHARGE(rbPotential%nSites))
          allocate(rbPotential%C6(rbPotential%nSites),rbPotential%C12(rbPotential%nSites))
          
          rbPotential%site(1,:) = 0.0d0
!          rbPotential%site(2,:) = (/0.75695032726366118236d0,0.0d0,-0.58588227661829495041d0/)
!          rbPotential%site(3,:) = (/-0.75695032726366118236d0,0.0d0,-0.58588227661829495041d0/)
!          rbPotential%site(4,:) = (/0.0d0,0.57154330164408195802d0,0.40415127656087129759d0/)
!          rbPotential%site(5,:) = (/0.0d0,-0.57154330164408195802d0,0.40415127656087129759d0/)
         
          rbPotential%site(2,:) = (/0.2785156560837814d0, -0.6990486725944498d0, 0.5916010671560346d0/)
          rbPotential%site(3,:) = (/-0.0640734628674161d0, -0.4219928632410997D0, -0.8567662777734411D0/)
          rbPotential%site(4,:) = (/0.4728396101454916D0, 0.5159942465174050D0, -0.0131392784579428d0/)
          rbPotential%site(5,:) = (/-0.6207653788462422D0, 0.2573187309647155D0, 0.1960546227983159d0/)

          rbPotential%CHARGE = (/0.0d0, 0.241D0, 0.241D0, -0.241D0, -0.241D0/)

! LJ coefficients in kJ/mol Angstrom**6 or Angstrom**12
          rbPotential%C6 = (/2470.012857D0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/) 
          rbPotential%C12 = (/2278383.244D0, 0.0d0, 0.0d0, 0.0d0, 0.0d0/)

       case default
          print *, 'Unrecognised atom type passed to rigid body initialisation:', atomType
          STOP ! call exit(1)

       end select

! PAH molecule
    case ('P')
       rbPotential%nSites = 36
       ALLOCATE(rbPotential%site(rbPotential%nSites,3), rbPotential%charge(rbPotential%nSites))

       rbPotential%site(:,3) = 0.0d0

       rbPotential%site( 1,1)=1.42d0
       rbPotential%site( 1,2)=0.0d0
       rbPotential%CHARGE( 1)=-0.0066d0

       rbPotential%site( 2,1)=0.71d0
       rbPotential%site( 2,2)=1.229755d0
       rbPotential%CHARGE( 2)=-0.0066d0

       rbPotential%site( 3,1)=-0.71d0
       rbPotential%site( 3,2)=1.229755d0
       rbPotential%CHARGE( 3)=-0.0066d0

       rbPotential%site( 4,1)=-1.42d0
       rbPotential%site( 4,2)=0.0d0
       rbPotential%CHARGE( 4)=-0.0066D0

       rbPotential%site( 5,1)=-0.71d0
       rbPotential%site( 5,2)=-1.229755d0
       rbPotential%CHARGE( 5)=-0.0066D0

       rbPotential%site( 6,1)=0.71d0
       rbPotential%site( 6,2)=-1.229755d0
       rbPotential%CHARGE( 6)=-0.0066D0

       rbPotential%site( 7,1)=2.838d0
       rbPotential%site( 7,2)=0.0d0
       rbPotential%CHARGE( 7)=0.0158d0

       rbPotential%site( 8,1)=1.419d0
       rbPotential%site( 8,2)=2.457779d0
       rbPotential%CHARGE( 8)=0.0158d0

       rbPotential%site( 9,1)=-1.419d0
       rbPotential%site( 9,2)=2.457779d0
       rbPotential%CHARGE( 9)=0.0158d0

       rbPotential%site(10,1)=-2.838d0
       rbPotential%site(10,2)=0.0d0
       rbPotential%CHARGE(10)=0.0158d0

       rbPotential%site(11,1)=-1.419d0
       rbPotential%site(11,2)=-2.45779d0
       rbPotential%CHARGE(11)=0.0158d0

       rbPotential%site(12,1)=1.419d0
       rbPotential%site(12,2)=-2.45779d0
       rbPotential%CHARGE(12)=0.0158d0

       rbPotential%site(13,1)=3.523366d0
       rbPotential%site(13,2)=1.248219d0
       rbPotential%CHARGE(13)=-0.0986d0

       rbPotential%site(14,1)=2.842673d0
       rbPotential%site(14,2)=2.427211d0
       rbPotential%CHARGE(14)=-0.0986d0

       rbPotential%site(15,1)=0.680706d0
       rbPotential%site(15,2)=3.675430d0
       rbPotential%CHARGE(15)=-0.0986d0

       rbPotential%site(16,1)=-0.680677d0
       rbPotential%site(16,2)=3.675437d0
       rbPotential%CHARGE(16)=-0.0986d0

       rbPotential%site(17,1)=-2.842673d0
       rbPotential%site(17,2)=2.427211d0
       rbPotential%CHARGE(17)=-0.0986D0

       rbPotential%site(18,1)=-3.523366D0
       rbPotential%site(18,2)=1.248219D0
       rbPotential%CHARGE(18)=-0.0986D0

       rbPotential%site(19,1)=-3.523366D0
       rbPotential%site(19,2)=-1.248219D0
       rbPotential%CHARGE(19)=-0.0986D0

       rbPotential%site(20,1)=-2.842673D0
       rbPotential%site(20,2)=-2.427211D0
       rbPotential%CHARGE(20)=-0.0986D0

       rbPotential%site(21,1)=-0.680677D0
       rbPotential%site(21,2)=-3.675437D0
       rbPotential%CHARGE(21)=-0.0986D0

       rbPotential%site(22,1)=0.680706D0
       rbPotential%site(22,2)=-3.675437D0
       rbPotential%CHARGE(22)=-0.0986D0

       rbPotential%site(23,1)=2.842673D0
       rbPotential%site(23,2)=-2.427211D0
       rbPotential%CHARGE(23)=-0.0986D0

       rbPotential%site(24,1)=3.523366D0
       rbPotential%site(24,2)=-1.248219D0
       rbPotential%CHARGE(24)=-0.0986D0

       rbPotential%site(25,1)=4.602048D0
       rbPotential%site(25,2)=1.274393D0
       rbPotential%CHARGE(25)=0.094D0

       rbPotential%site(26,1)=3.404682D0
       rbPotential%site(26,2)=3.348290D0
       rbPotential%CHARGE(26)=0.094D0

       rbPotential%site(27,1)=1.197383D0
       rbPotential%site(27,2)=4.622681D0
       rbPotential%CHARGE(27)=0.094D0

       rbPotential%site(28,1)=-1.197347D0
       rbPotential%site(28,2)=4.622693D0
       rbPotential%CHARGE(28)=0.094D0

       rbPotential%site(29,1)=-3.404682D0
       rbPotential%site(29,2)=3.348290D0
       rbPotential%CHARGE(29)=0.094D0

       rbPotential%site(30,1)=-4.602048D0
       rbPotential%site(30,2)=1.274393D0
       rbPotential%CHARGE(30)=0.094D0

       rbPotential%site(31,1)=-4.602048D0
       rbPotential%site(31,2)=-1.274393D0
       rbPotential%CHARGE(31)=0.094D0

       rbPotential%site(32,1)=-3.404682D0
       rbPotential%site(32,2)=-3.348290D0
       rbPotential%CHARGE(32)=0.094D0

       rbPotential%site(33,1)=-1.197347D0
       rbPotential%site(33,2)=-4.622693D0
       rbPotential%CHARGE(33)=0.094D0

       rbPotential%site(34,1)=1.197383D0
       rbPotential%site(34,2)=-4.622681D0
       rbPotential%CHARGE(34)=0.094D0

       rbPotential%site(35,1)=3.404682D0
       rbPotential%site(35,2)=-3.348290D0
       rbPotential%CHARGE(35)=0.094D0

       rbPotential%site(36,1)=4.602048D0
       rbPotential%site(36,2)=-1.274393D0
       rbPotential%CHARGE(36)=0.094D0

    case ('C')
       if (numRBtypes.ne.1) then
          print *, 'Multiple capsid types not yet initialised'
          stop
       endif

       rbPotential%nSites = capsomerDefs(1)%nBasalSites + 1
       rbPotential%nPhysicalSites = rbPotential%nSites
       
       allocate(rbPotential%site(rbPotential%nSites,3), rbPotential%siteLabel(rbPotential%nSites))
       allocate(rbPotential%C6(rbPotential%nSites), rbPotential%C12(rbPotential%nSites))
       allocate(rbPotential%morseParameters(rbPotential%nSites,3))
       allocate(rbPotential%mass(rbPotential%nSites))

       rbPotential%C6 = 0.0d0
       rbPotential%C12 = 0.0d0

! Sites around the base
! Morse parameters are in the order epsilon, rho, Re
! At the moment epsilon and Re are implicitly the units of energy and distance, respectively

       do currentSite = 1, capsomerDefs(1)%nBasalSites
          currentAngle = (currentSite-1)*twopi/(1.0d0*capsomerDefs(1)%nBasalSites)

          rbPotential%site(currentSite,:) = (/capsomerDefs(1)%radius*cos(currentAngle), &
                                              capsomerDefs(1)%radius*sin(currentAngle), &
                                              0.0d0/)
          
          rbPotential%siteLabel(currentSite) = 'C';
          rbPotential%morseParameters(currentSite,:) = (/1.0d0,capsomerDefs(1)%rho,1.0d0/)
          rbPotential%mass(currentSite) = 0.6d0
       enddo

! Then the repulsive site at the apex
! epsilon_repulsive is specified as a fraction of epsilon_Morse

       rbPotential%site(rbPotential%nSites,:) = (/0.0d0, 0.0d0, capsomerDefs(1)%height/)
       rbPotential%siteLabel(rbPotential%nSites) = 'O';
       rbPotential%C12(rbPotential%nSites) = &
         capsomerDefs(1)%epsilonRep*1.0d0*(1.0d0+capsomerDefs(1)%radius*sqrt((sqrt(5.0d0)+5.0d0)/2.0d0))**12
       rbPotential%morseParameters(rbPotential%nSites,:) = 0.0d0
       rbPotential%mass(rbPotential%nSites) = 3d0

       rbPotential%freqConvFactor = 1d0    ! No frequency conversion just yet

!       do currentSite = 1, rbPotential%nSites
!          print *, rbPotential%site(currentSite,1), rbPotential%site(currentSite,2), rbPotential%site(currentSite,3)
!       enddo
!       stop
    case default
       print *, 'Unrecognised atom type passed to rigid body initialisation:', atomType
       STOP ! call exit(1)

    end select

! Optional specification of system masses and atomic numbers
! OPTIM uses this but GMIN doesn`t

    if (.not.allocated(rbPotential%mass)) then
       print *, 'initialiseRigidBody> No masses specified so initialising to default value (1.0)'
       allocate(rbPotential%mass(rbPotential%nPhysicalSites))
       rbPotential%mass = 1.0d0
    endif
 
    if (.not.allocated(rbPotential%atomicNumber)) then
       print *, 'initialiseRigidBody> No atomic numbers specified so initialising to default value (0)'
       allocate(rbPotential%atomicNumber(rbPotential%nPhysicalSites))
       rbPotential%atomicNumber = 0
    endif
      
    if (present(systemMasses)) then
       if (.not.present(nbodies)) then
          print *, 'initialiseRigidBody> Must specify number of bodies if you wish to initialise system masses'
          STOP ! call exit(1)
       endif

       allocate(systemMasses(nbodies*rbPotential%nPhysicalSites))

       do currentBody = 1, nbodies
          systemMasses((currentBody-1)*rbPotential%nPhysicalSites+1:currentBody*rbPotential%nPhysicalSites) = &
               rbPotential%mass
       enddo
    endif

    if (present(systemAtNumbers)) then
       if (.not.present(nbodies)) then
          print *, 'initialiseRigidBody> Must specify number of bodies if you wish to initialise system atomic numbers'
          STOP ! call exit(1)
       endif

       allocate(systemAtNumbers(nbodies*rbPotential%nPhysicalSites))
       do currentBody = 1, nbodies
          systemAtNumbers((currentBody-1)*rbPotential%nPhysicalSites+1:currentBody*rbPotential%nPhysicalSites) = &
               rbPotential%atomicNumber
       enddo
    endif

    call initialiseInteractionMap

  end subroutine initialiseRigidBody

! #####################################################################################################

! Saves us having to check explicitly in the potential evaluation routine

  subroutine initialiseInteractionMap
    implicit none

    integer :: site1, site2

    if (allocated(sitesInteract)) then
       if (size(sitesInteract,1).ne.rbPotential%nSites .or. size(sitesInteract,2).ne.rbPotential%nSites) then
          deallocate(sitesInteract)
          allocate(sitesInteract(rbPotential%nSites, rbPotential%nSites, 4))
       endif
    else
       allocate(sitesInteract(rbPotential%nSites,rbPotential%nSites, 4))
    endif

    sitesInteract = .false.

    do site1 = 1, rbPotential%nSites
       do site2 = 1, rbPotential%nSites

! Coulomb interaction

          if (allocated(rbPotential%charge)) then
             if (rbPotential%charge(site1).ne.0.0d0 .and. rbPotential%charge(site2).ne.0.0d0) then
                sitesInteract(site1, site2, 1) = .true.
                sitesInteract(site1, site2, 2) = .true.   
             endif
          endif

! LJ interaction

          if (allocated(rbPotential%C6) .and. allocated(rbPotential%C12)) then
              if ((rbPotential%C6(site1).ne.0.0d0 .or. rbPotential%C12(site1).ne.0.0d0) .and. &
                   (rbPotential%C6(site2).ne.0.0d0 .or. rbPotential%C12(site2).ne.0.0d0)) then
                 sitesInteract(site1, site2, 1) = .true.
                 sitesInteract(site1, site2, 3) = .true.
              endif
          endif

! Morse interaction

          if (allocated(rbPotential%morseParameters)) then
              if (rbPotential%morseParameters(site1,1).ne.0.0d0 .and. &
                   rbPotential%morseParameters(site2,1).ne.0.0d0) then
                 sitesInteract(site1, site2, 1) = .true.
                 sitesInteract(site1, site2, 4) = .true.
              endif
          endif
       enddo
    enddo

  end subroutine initialiseInteractionMap

! #####################################################################################################

! Return allocated memory on request

  subroutine cleanRigidBodies
    implicit none
    
    if (allocated(rbPotential%site)) deallocate(rbPotential%site)
    if (allocated(rbPotential%siteLabel)) deallocate(rbPotential%siteLabel)
    if (allocated(rbPotential%mass)) deallocate(rbPotential%mass)
    if (allocated(rbPotential%atomicNumber)) deallocate(rbPotential%atomicNumber)
    if (allocated(rbPotential%CHARGE)) deallocate(rbPotential%CHARGE)
    if (allocated(rbPotential%C6)) deallocate(rbPotential%C6)
    if (allocated(rbPotential%C12)) deallocate(rbPotential%C12)
    if (allocated(rbPotential%morseParameters)) deallocate(rbPotential%morseParameters)

    if (allocated(sitesInteract)) deallocate(sitesInteract)

  end subroutine cleanRigidBodies

! #####################################################################################################

! Moves the center of mass of the reference geometry to the origin
! This is necessary for some implementations of the normal mode analysis

  subroutine CoMtoOrigin ()

    implicit none

    real (kind=kind(0.0d0)) :: initialCoM(3)
    integer :: currentSite

! Check we have relevant data

    if (.not.allocated(rbPotential%site)) then
       print *, 'Cannot perform CoM calculation without specified sites'
       STOP ! call exit(1)
    else if (.not.allocated(rbPotential%mass)) then
       print *, 'Cannot perform CoM calculation without specified masses'
       STOP ! call exit(1)
    endif

! Calculate initial CoM - use only physical sites

    initialCoM = 0.0d0

    do currentSite = 1, rbPotential%nPhysicalSites
       initialCoM = initialCoM + rbPotential%mass(currentSite)*rbPotential%site(currentSite,:)
    enddo

    initialCoM = initialCoM/monomerMass()

! Correct all sites

    do currentSite = 1, rbPotential%nSites
       rbPotential%site(currentSite,:) = rbPotential%site(currentSite,:) - initialCoM
    enddo

    return

  end subroutine CoMtoOrigin

! #####################################################################################################

! Cartesian generation routines

  function cartesianX (refSite,xCoC,px,py,pz)
    implicit none
    real (kind=kind(0.0d0)) :: cartesianX
    real (kind=kind(0.0d0)), intent(IN) :: refSite(3),xCoC,px,py,pz 

   real (kind=kind(0.0d0)) :: angle, cos_angle, cos_factor, sin_factor 

! Normalise the rotational variables

    angle = DSQRT(px**2 + py**2 + pz**2)
    cos_angle = cos(angle)

    if (angle.lt.smallAngle) then
!       print *, 'Small angle approximation in cartesianX'
       cos_factor = 0.5d0 - (angle**2)/24.0d0
       sin_factor = 1.0d0 - (angle**2)/6.0d0
    else
       cos_factor = (1.0d0-cos_angle)/(angle**2)
       sin_factor = sin(angle)/angle
    endif

    cartesianX =  xCoC + refSite(1)*cos_angle + &
                  px*(px*refSite(1)+py*refSite(2)+pz*refSite(3))*cos_factor + &
                  (refSite(2)*pz - refSite(3)*py)*sin_factor
    return
  end function cartesianX

! #####################################################################################################

  function cartesianY (refSite,yCoC,px,py,pz)
    implicit none
    real (kind=kind(0.0d0)) :: cartesianY
    real (kind=kind(0.0d0)), intent(IN) :: refSite(3),yCoC,px,py,pz 

    real (kind=kind(0.0d0)) :: angle, cos_angle, cos_factor, sin_factor

! Normalise the rotational variables

    angle = DSQRT(px**2 + py**2 + pz**2)
    cos_angle = cos(angle)   

    if (angle.lt.smallAngle) then
!       print *, 'Small angle approximation in cartesianY'
       cos_factor = 0.5d0 - (angle**2)/24.0d0
       sin_factor = 1.0d0 - (angle**2)/6.0d0
    else
       cos_factor = (1.0d0-cos_angle)/(angle**2)
       sin_factor = sin(angle)/angle
    endif

    cartesianY = yCoC + refSite(2)*cos_angle + &
                 py*(px*refSite(1)+py*refSite(2)+pz*refSite(3))*cos_factor + &
                 (refSite(3)*px - refSite(1)*pz)*sin_factor

    return
  end function cartesianY

! #####################################################################################################

  function cartesianZ (refSite,zCoC,px,py,pz)
    implicit none
    real (kind=kind(0.0d0)) :: cartesianZ
    real (kind=kind(0.0d0)), intent(IN) :: refSite(3),zCoC,px,py,pz 

    real (kind=kind(0.0d0)) :: angle, cos_angle, cos_factor, sin_factor

! Normalise the rotational variables

    angle = DSQRT(px**2 + py**2 + pz**2)
    cos_angle = cos(angle)   

    if (angle.lt.smallAngle) then
!       print *, 'Small angle approximation in cartesianZ'
       cos_factor = 0.5d0 - (angle**2)/24.0d0
       sin_factor = 1.0d0 - (angle**2)/6.0d0
    else
       cos_factor = (1.0d0-cos_angle)/(angle**2)
       sin_factor = sin(angle)/angle
    endif

    cartesianZ = zCoC + refSite(3)*cos_angle + &
                 pz*(px*refSite(1)+py*refSite(2)+pz*refSite(3))*cos_factor + &
                 (refSite(1)*py - refSite(2)*px)*sin_factor
    return
  end function cartesianZ

! #####################################################################################################

  function cartesianSeparation(refSite1,xCoC1,yCoC1,zCoC1,px1,py1,pz1, &
                               refSite2,xCoC2,yCoC2,zCoC2,px2,py2,pz2)

    implicit none

    real (kind=kind(0.0d0)) :: cartesianSeparation
    real (kind=kind(0.0d0)), intent(IN) :: refSite1(3),xCoC1,yCoC1,zCoC1,px1,py1,pz1
    real (kind=kind(0.0d0)), intent(IN) :: refSite2(3),xCoC2,yCoC2,zCoC2,px2,py2,pz2
    
    cartesianSeparation = DSQRT((cartesianX(refSite1,xCoC1,px1,py1,pz1) -     &
                                 cartesianX(refSite2,xCoC2,px2,py2,pz2))**2 + &
                                (cartesianY(refSite1,yCoC1,px1,py1,pz1) -     &
                                 cartesianY(refSite2,yCoC2,px2,py2,pz2))**2 + &
                                (cartesianZ(refSite1,zCoC1,px1,py1,pz1) -     &
                                 cartesianZ(refSite2,zCoC2,px2,py2,pz2))**2)
    
    return
  end function cartesianSeparation

! #####################################################################################################

! Convert to Cartesians for the purposes of calculating path lengths etc
! Assumes that the physical sites are specified first in the list

  subroutine monomerToCartesians (xCoC,yCoC,zCoC,px,py,pz,cartesians)
    implicit none

! Subroutine arguments

    real (kind=kind(0.0d0)), intent(IN) :: xCoC,yCoC,zCoC,px,py,pz
    real (kind=kind(0.0d0)), intent(OUT) :: cartesians (3*rbPotential%nPhysicalSites)
    
! Local variables

    integer :: currentSite
    real (kind=kind(0.0d0)) :: angle, pdotX0, cos_angle, cos_factor, sin_factor  

! Work out the angle, this only needs to be done once

    angle = DSQRT(px**2 + py**2 + pz**2)
    cos_angle = cos(angle)

    if (angle.lt.smallAngle) then
!       print *, 'Small angle approximation in monomerToCartesians'
       cos_factor = 0.5d0 - (angle**2)/24.0d0
       sin_factor = 1.0d0 - (angle**2)/6.0d0
    else
       cos_factor = (1.0d0-cos_angle)/(angle**2)
       sin_factor = sin(angle)/angle
    endif

! Run through the sites

    do currentSite=1, rbPotential%nPhysicalSites
       pdotX0 = px*rbPotential%site(currentSite,1)+py*rbPotential%site(currentSite,2)+pz*rbPotential%site(currentSite,3)

       cartesians(3*currentSite-2) =  xCoC + rbPotential%site(currentSite,1)*cos_angle + &
            px*pdotX0*cos_factor + &
            (rbPotential%site(currentSite,2)*pz - rbPotential%site(currentSite,3)*py)*sin_factor

       cartesians(3*currentSite-1) = yCoC + rbPotential%site(currentSite,2)*cos_angle + &
                 py*pdotX0*cos_factor + &
                 (rbPotential%site(currentSite,3)*px - rbPotential%site(currentSite,1)*pz)*sin_factor

       cartesians(3*currentSite) = zCoC + rbPotential%site(currentSite,3)*cos_angle + &
                 pz*pdotX0*cos_factor + &
                 (rbPotential%site(currentSite,1)*py - rbPotential%site(currentSite,2)*px)*sin_factor

!       cartesians(3*currentSite-2) = cartesianX(rbPotential%site(currentSite,:),xCoC,px,py,pz) 
!       cartesians(3*currentSite-1) = cartesianY(rbPotential%site(currentSite,:),yCoC,px,py,pz)
!       cartesians(3*currentSite) =   cartesianZ(rbPotential%site(currentSite,:),zCoC,px,py,pz)
    enddo

  end subroutine monomerToCartesians

! #####################################################################################################

  subroutine systemToCartesians (nbodies, angleAxisCoords, cartesianCoords)
    implicit none

! Subroutine arguments

    integer, intent(IN) :: nbodies
    real (kind=kind(0.0d0)), intent(IN) :: angleAxisCoords(6*nbodies)
    real (kind=kind(0.0d0)), intent(OUT) :: cartesianCoords(nbodies*rbPotential%nPhysicalSites*3)  
    
! Local variables

    integer :: currentBody

    cartesianCoords = 0.0d0

    do currentBody = 1, nbodies
       call monomerToCartesians(angleAxisCoords(3*currentBody-2), &
                                angleAxisCoords(3*currentBody-1), &
                                angleAxisCoords(3*currentBody), &
                                angleAxisCoords(3*(currentBody+nbodies)-2), &
                                angleAxisCoords(3*(currentBody+nbodies)-1), &
                                angleAxisCoords(3*(currentBody+nbodies)), &
         cartesianCoords((currentBody-1)*rbPotential%nPhysicalSites*3+1:currentBody*rbPotential%nPhysicalSites*3))
    enddo

  end subroutine systemToCartesians

! #####################################################################################################

  subroutine writeToCartesianStream (writeUnit, nbodies, angleAxisCoords, comment)

    implicit none

! Subroutine arguments

    integer, intent(IN) :: writeUnit
    integer, intent(IN) :: nbodies
    real (kind=kind(0.0d0)), intent(IN) :: angleAxisCoords(6*nbodies)
    character(*), intent(IN) :: comment

! Local variables

    real (kind=kind(0.0d0)) :: cartesians(3*nbodies*rbPotential%nPhysicalSites)
    integer :: nCartPoints, i

! Convert to Cartesians and write to specified unit

    call systemToCartesians(nbodies, angleAxisCoords, cartesians)

    nCartPoints = nbodies*rbPotential%nPhysicalSites
    write(writeUnit, *) nCartPoints
    write(writeUnit, *) trim(comment)

    do i = 1, nCartPoints
       WRITE(writeUnit,'(A5,4X,3g20.10)') rbPotential%siteLabel(MOD(i-1,rbPotential%nPhysicalSites)+1), &
            cartesians(3*(i-1)+1), cartesians(3*(i-1)+2), cartesians(3*(i-1)+3)
    enddo

  end subroutine writeToCartesianStream

! #####################################################################################################

! Converts a single monomer from Cartesians to angle-axis coordinates

  subroutine monomerToAA (cartesians, AAtrans, AArot)

    implicit none

    real (kind=kind(0.0d0)), intent(IN) :: cartesians(rbPotential%nPhysicalSites,3)
    real (kind=kind(0.0d0)), intent(OUT) :: AAtrans(3), AArot(3)    

    real (kind=kind(0.0d0)) :: fittingRMSD    
    real (kind=kind(0.0d0)) :: rmsdTol = 1D-4
    integer :: i, j  

! Simply uses the quaternion fitting method to determine the necessary parameters, not very efficient

    call quaternionMatch(rbPotential%nPhysicalSites, cartesians, rbPotential%site(1:rbPotential%nPhysicalSites,:), &
                         fittingRMSD, AArot, AAtrans)

    if (fittingRMSD.gt.rmsdTol) then
       print *, 'Error in attempting to convert Cartesians to angle-axis coordinates: RMSD is too high', fittingRMSD
!       STOP ! call exit(1)
    endif

  end subroutine monomerToAA

! #####################################################################################################

  subroutine systemToAA (nbodies, cartesians, AAcoords)

    implicit none
    
! Subroutine arguments

    integer, intent(IN) :: nbodies
    real (kind=kind(0.0d0)), intent(IN) :: cartesians(rbPotential%nPhysicalSites*nbodies,3)
    real (kind=kind(0.0d0)), intent(OUT) :: AAcoords(6*nbodies)

! Local variables

    integer :: currentBody

    do currentBody=1,nbodies
       call monomerToAA(cartesians(rbPotential%nPhysicalSites*(currentBody-1)+1:rbPotential%nPhysicalSites*currentBody,:), &
                        AAcoords(3*currentBody-2:3*currentBody), &
                        AAcoords(3*(currentBody+nbodies)-2:3*(currentBody+nbodies)))
    enddo
    
  end subroutine systemToAA

! #####################################################################################################

  function monomerMass ()

    implicit none

    real (kind=kind(0.0d0)) :: monomerMass
 
    monomerMass = SUM(rbPotential%mass)
    return

  end function monomerMass

! #####################################################################################################

! Interaction functions and derivatives, note no unit conversion is attempted

! #####################################################################################################

! Lennard-Jones interactions

  function fLJ(k12,k6,distance)
    implicit none
    real (kind=kind(0.0d0)) :: fLJ
    real (kind=kind(0.0d0)), intent(IN) :: k12, k6    ! C12 and C6 coefficients 
    real (kind=kind(0.0d0)), intent(IN) :: distance

    real (kind=kind(0.0d0)) :: invR6

    invR6 = 1.0d0/(distance**6)
    fLJ = invR6*(k12*invR6 - k6)

    return

  end function fLJ

! #####################################################################################################

  function dfLJdR(k12,k6,distance)
    implicit none
    real (kind=kind(0.0d0)) :: dfLJdR
    real (kind=kind(0.0d0)), intent(IN) :: k12, k6    ! C12 and C6 coefficients 
    real (kind=kind(0.0d0)), intent(IN) :: distance

    real (kind=kind(0.0d0)) :: invR, invR6

    invR = 1.0d0/distance
    invR6 = invR**6

    dfLJdR = -6.0D0*(-k6+2.0D0*k12*invR6)*invR*invR6
!    dfLJdR = -6.0D0*(-k6+2.0D0*k12/(distance**6))/(distance**7)

    return

  end function dfLJdR

! #####################################################################################################

  function d2fLJdR2(k12,k6,distance)
    implicit none
    real (kind=kind(0.0d0)) :: d2fLJdR2
    real (kind=kind(0.0d0)), intent(IN) :: k12, k6    ! C12 and C6 coefficients 
    real (kind=kind(0.0d0)), intent(IN) :: distance

    real (kind=kind(0.0d0)) :: invR2, invR6

    invR2 = 1.0d0/(distance**2)
    invR6 = invR2**3

    d2fLJdR2 = 6.0D0*(26.0d0*k12*invR6 - 7.0d0*k6)*invR2*invR6
!    d2fLJdR2 = 6.0D0*(26.0d0*k12/distance**6 - 7.0d0*k6)/distance**8
    return

  end function d2fLJdR2

! #####################################################################################################

! Coulombic interactions

  function fCoulomb (q1, q2, distance)
    implicit none
    real (kind=kind(0.0d0)) :: fCoulomb
    real (kind=kind(0.0d0)), intent(IN) :: q1, q2
    real (kind=kind(0.0d0)), intent(IN) :: distance

    fCoulomb = q1*q2/distance

    return
  end function fCoulomb

! #####################################################################################################

  function dfCoulombdR (q1, q2, distance)
    implicit none
    real (kind=kind(0.0d0)) :: dfCoulombdR
    real (kind=kind(0.0d0)), intent(IN) :: q1, q2
    real (kind=kind(0.0d0)), intent(IN) :: distance

    dfCoulombdR = -q1*q2/(distance**2)

    return
  end function dfCoulombdR

! #####################################################################################################

  function d2fCoulombdR2 (q1, q2, distance)
    implicit none
    real (kind=kind(0.0d0)) :: d2fCoulombdR2
    real (kind=kind(0.0d0)), intent(IN) :: q1, q2
    real (kind=kind(0.0d0)), intent(IN) :: distance

    d2fCoulombdR2 = 2.0d0*q1*q2/(distance**3)

    return
  end function d2fCoulombdR2

! #####################################################################################################

! Note this Morse function is shifted down by epsilon

  function fMorse(epsilon, rho, Re, distance)
    implicit none
    real (kind=kind(0.0d0)) :: fMorse
    real (kind=kind(0.0d0)), intent(IN) :: epsilon, rho, Re   ! Morse parameters
    real (kind=kind(0.0d0)), intent(IN) :: distance    

    fMorse = epsilon*((1.0d0-exp(rho*(1.0d0 - distance/Re)))**2 - 1.0d0)

  end function fMorse

! #####################################################################################################

  function dfMorsedR(epsilon, rho, Re, distance)
    implicit none
    real (kind=kind(0.0d0)) :: dfMorsedR
    real (kind=kind(0.0d0)), intent(IN) :: epsilon, rho, Re   ! Morse parameters
    real (kind=kind(0.0d0)), intent(IN) :: distance 

    real (kind=kind(0.0d0)) :: expFactor

    expFactor = exp(rho*(1.0d0-distance/Re))

    dfMorsedR = 2.0d0*epsilon*rho*(1.0d0 - expFactor)*expFactor/Re

  end function dfMorsedR

! #####################################################################################################

  function d2fMorsedR2(epsilon, rho, Re, distance)
    implicit none
    real (kind=kind(0.0d0)) :: d2fMorsedR2
    real (kind=kind(0.0d0)), intent(IN) :: epsilon, rho, Re   ! Morse parameters
    real (kind=kind(0.0d0)), intent(IN) :: distance 

    real (kind=kind(0.0d0)) :: expFactor

    expFactor = exp(rho*(1.0d0-distance/Re))

    d2fMorsedR2 = 2.0d0*epsilon*(rho**2)*expFactor*(2.0d0*expFactor - 1.0d0)/(Re**2)

  end function d2fMorsedR2

! #####################################################################################################

  subroutine updateC12(newValue) 

    implicit none

    real (kind=kind(0.0d0)), intent(IN) :: newValue
    integer :: currentSite
    
    if (allocated(rbPotential%C12)) then
       do currentSite = 1, rbPotential%nSites
          if (rbPotential%C12(currentSite).ne.0.0d0) rbPotential%C12(currentSite) = newValue
       enddo
    endif

    call initialiseInteractionMap

  end subroutine updateC12

! #####################################################################################################

  subroutine updateC6(newValue) 

    implicit none

    real (kind=kind(0.0d0)), intent(IN) :: newValue
    integer :: currentSite
    
    if (.not.allocated(rbPotential%C6)) then
       allocate(rbPotential%C6(rbPotential%nSites))
       rbPotential%C6 = 0.0d0
    endif

    rbPotential%C6(rbPotential%nSites) = newValue

    call initialiseInteractionMap

  end subroutine updateC6

! #####################################################################################################

end module rigidbody
