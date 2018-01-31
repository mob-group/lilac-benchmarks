!!****h* formstep/newlbfgs
!!
!! NAME
!! L-BFGS
!!
!! FUNCTION
!! Optimisation algorithms: determine a search direction. 
!!
!!
!!       LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
!!                         JORGE NOCEDAL
!!                       *** July 1990 ***
!!
!!
!!    This subroutine solves the unconstrained minimization problem
!!
!!                     min F(x),    x= (x1,x2,...,xN),
!!
!!     using the limited memory BFGS method. The routine is especially
!!     effective on problems involving a large number of variables. In
!!     a typical iteration of this method an approximation Hk to the
!!     inverse of the Hessian is obtained by applying M BFGS updates to
!!     a diagonal matrix Hk0, using information from the previous M steps.
!!     The user specifies the number M, which determines the amount of
!!     storage required by the routine. 
!!     The algorithm is described in "On the limited memory BFGS method
!!     for large scale optimization", by D. Liu and J. Nocedal,
!!     Mathematical Programming B 45 (1989) 503-528.
!!
!!    M (Nmem) is an INTEGER variable that must be set by the user to
!!            the number of corrections used in the BFGS update. It
!!            is not altered by the routine. Values of M less than 3 are
!!            not recommended; large values of M will result in excessive
!!            computing time. 3<= M <=7 is recommended. Restriction: M>0.
!!
!! An f77 version of this file was originally obtained from
!! http://www.ece.northwestern.edu/~nocedal/lbfgs.html
!! Condition for Use: This software is freely available for educational
!! or commercial purposes. We expect that all publications describing 
!! work using this software quote at least one of the references:
!!
!! J. Nocedal. Updating Quasi-Newton Matrices with Limited Storage (1980),
!! Mathematics of Computation 35, pp. 773-782.
!!
!! D.C. Liu and J. Nocedal. On the Limited Memory Method for Large Scale
!! Optimization (1989), Mathematical Programming B, 45, 3, pp. 503-528.
!!
!! For normal use in one instance, use dlf_newlbfgs_init, _step, and _destroy.
!! To restart the optimiser, use dlf_newlbfgs_restart.
!!
!! If two overlapping optimisations should use L-BFGS, the module can be started
!! in more instances: use dlf_newlbfgs_select to select an instance. Always the current
!! instance will be affected by dlf_newlbfgs_step, dlf_newlbfgs_restart and dlf_newlbfgs_destroy. Use
!! dlf_newlbfgs_deselect to select the main instance (first instance). To invoke a new instance,
!! use dlf_newlbfgs_deselect("newname",.true.) and then dlf_newlbfgs_init.
!!
!!**********************************************************************

!! THE VARIABLES IN THIS FILE WERE ADAPTED DURING THE INTERFACING PROCESS TO OPTIM BY Judith Rommel jbr36
!! DATA
!! $Date: 2011-01-20 13:43:17 +0000 (Thu, 20 Jan 2011) $
!! $Rev: 465 $
!! $Author: twk $
!! $URL: http://ccpforge.cse.rl.ac.uk/svn/dl-find/branches/release_chemsh3.5/dlf_lbfgs.f90 $
!! $Id: dlf_lbfgs.f90 465 2011-01-20 13:43:17Z twk $
!!
!! COPYRIGHT
!!
!!  Copyright 2007 Johannes Kaestner (kaestner@theochem.uni-stuttgart.de),
!!  Tom Keal (thomas.keal@stfc.ac.uk)
!!
!!  This file is part of DL-FIND.
!!
!!  DL-FIND is free software: you can redistribute it and/or modify
!!  it under the terms of the GNU Lesser General Public License as 
!!  published by the Free Software Foundation, either version 3 of the 
!!  License, or (at your option) any later version.
!!
!!  DL-FIND is distributed in the hope that it will be useful,
!!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!  GNU Lesser General Public License for more details.
!!
!!  You should have received a copy of the GNU Lesser General Public 
!!  License along with DL-FIND.  If not, see 
!!  <http://www.gnu.org/licenses/>.
!!
!!****

MODULE newlbfgs_MODULE
  USE dlf_parameter_module, only: rk 
  type newlbfgs_type
    integer                 :: N ! number of variables
    integer                 :: M ! number steps to remember
    REAL(RK), ALLOCATABLE   :: store(:)  ! N WORK 1 - N
    REAL(RK), ALLOCATABLE   :: store2(:) ! N old coords
    REAL(RK), ALLOCATABLE   :: rho(:)    ! M WORK N+1 - N+M
    REAL(RK), ALLOCATABLE   :: alpha(:)  ! M WORK N+M+1 - N+2M
    REAL(RK), ALLOCATABLE   :: step(:,:) ! M,N WORK N+2M+1 - N+2M+NM
    REAL(RK), ALLOCATABLE   :: dgrad(:,:)! M,N WORK N+2M+NM+1 - N+2M+2NM
    logical                 :: tprecon ! is precon set (and allocated)?
    real(rk), allocatable   :: precon(:,:) ! (N,N) guess for the initial inverse hessian
                                ! full matrix, not just the diagonal. If this is used,
                                ! the algorithm gets order N^2 !
    
    integer                 :: point ! CURRENT POSITION IN THE WORK ARRAY
    INTEGER                 :: iter ! number of iteration
    logical                 :: tinit
    character(40)           :: tag
    type(newlbfgs_type),pointer  :: next
  end type newlbfgs_type
  type(newlbfgs_type),pointer,save :: newlbfgs
  type(newlbfgs_type),pointer,save :: newlbfgs_first
  logical, save               :: tinit=.false.
  LOGICAL,PARAMETER :: dbg=.false. ! write debug info 
  character(40),save     :: newtag="none"
END MODULE newlbfgs_MODULE

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* newlbfgs/dlf_newlbfgs_step
!!
!! FUNCTION
!!
!! Form an L-BFGS step from the input geometry, gradient, and the history
!!
!! SYNOPSIS
SUBROUTINE dlf_newlbfgs_step(x,g,step_)
!! SOURCE
  USE dlf_parameter_module, only: rk
  use dlf_global, only: stdout,printl
  USE newlbfgs_module
  IMPLICIT NONE
  !
  ! Dummy arguments
  !
  REAL(RK), intent(in)    :: x(newlbfgs%n) ! position (coords)
  REAL(RK), intent(in)    :: g(newlbfgs%n) ! gradient
  REAL(RK), intent(out)   :: step_(newlbfgs%n) ! Step
  !
  ! Local variables
  !
  REAL(RK) :: diag(newlbfgs%n)
  REAL(RK) :: beta , gnorm , sq , stp , yr , ys , yy 
  INTEGER  :: bound , cp , i 
  INTEGER  :: oldpoint,ivar
  real(RK) ,external :: ddot
! **********************************************************************
  if(.not.tinit) call dlf_fail("LBFGS not initialised!")
  if(.not.newlbfgs%tinit) then
    print*,"Instance of L-BFGS:",trim(newlbfgs%tag)
    call dlf_fail("This instance of LBFGS not initialised!")
  end if

  if(newlbfgs%iter==0) then ! first iteration, steepest descent!
    newlbfgs%point = 1
    oldpoint = 1

    newlbfgs%iter=1
    ! if the fist step should be smaller, include a factor here 
    step_(:) = -g(:)*0.02D0/dsqrt(sum(g**2))
    
    ! Store old gradient and coordinates
    newlbfgs%store(:) = g(:)
    newlbfgs%store2(:)= x(:)

    return

  end if

  ! ====================================================================
  ! All steps but first: calculate L-BFGS Step
  ! ====================================================================

  ! COMPUTE THE NEW STEP AND GRADIENT CHANGE
  newlbfgs%step(newlbfgs%point,:) = x(:) - newlbfgs%store2(:)
  newlbfgs%dgrad(newlbfgs%point,:) = g(:) - newlbfgs%store(:)

  oldpoint=newlbfgs%point

  ! take next point
  newlbfgs%point = newlbfgs%point + 1
  IF ( newlbfgs%point>newlbfgs%m ) newlbfgs%point = 1

  newlbfgs%iter = newlbfgs%iter + 1

  if(dbg) print*,"@100 newlbfgs%point=",newlbfgs%point
  if(dbg) print*,"newlbfgs%dgrad",newlbfgs%dgrad
  if(dbg) print*,"newlbfgs%step",newlbfgs%step

  bound = newlbfgs%iter - 1
  IF ( newlbfgs%iter>newlbfgs%m ) bound = newlbfgs%m
  ys = DDOT(newlbfgs%n,newlbfgs%dgrad(oldpoint,:),1,newlbfgs%step(oldpoint,:),1)
  yy = DDOT(newlbfgs%n,newlbfgs%dgrad(oldpoint,:),1,newlbfgs%dgrad(oldpoint,:),1)
  diag(:) = ys/yy ! default guess for diag
  !print*,"JK precon scale factor:",ys/yy

  if(dbg) print*,"Before 200: ys,yy",ys,yy

  ! ====================================================================
  ! COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
  ! "Updating quasi-Newton matrices with limited storage",
  ! Mathematics of Computation, Vol.24, No.151, pp. 773-782.
  ! ====================================================================
  !
  cp = newlbfgs%point
  newlbfgs%rho(oldpoint) = 1.D0/ys
  if(dbg) print*,"newlbfgs%rho",newlbfgs%rho
  if(dbg) print*,"bound",bound

  newlbfgs%store(:) = -g(:)
  cp = newlbfgs%point
  DO i = 1 , bound
    cp = cp - 1
    IF ( cp==0 ) cp = newlbfgs%m
    sq = DDOT(newlbfgs%n,newlbfgs%step(cp,:),1,newlbfgs%store(:),1)
    newlbfgs%alpha(cp) = newlbfgs%rho(cp)*sq
    !CALL DAXPY(newlbfgs%n,-newlbfgs%alpha(cp),newlbfgs%dgrad(cp,:),1,newlbfgs%store(:),1)
    newlbfgs%store(:)=newlbfgs%store(:)-newlbfgs%alpha(cp)*newlbfgs%dgrad(cp,:)
  END DO
  if(dbg) print*,"removed DAXPY"
  if(dbg) print*,"AFTER CALL TO DAXPY"
  !
  if(newlbfgs%tprecon) then
    ! diag= newlbfgs%precon * newlbfgs%store manually :-(
    do ivar=1,newlbfgs%n
      diag(ivar)=sum(newlbfgs%precon(:,ivar)*newlbfgs%store(:))
    end do
    !diag = matmul(newlbfgs%precon,newlbfgs%store)
    newlbfgs%store(:) = diag !* ys/yy
    !newlbfgs%store(:) = newlbfgs%store(:)*ys/yy
  else
    newlbfgs%store(:) = diag(:)*newlbfgs%store(:)
  end if
  !
  DO i = 1 , bound
    if(dbg) print*,"cp",cp
    yr = DDOT(newlbfgs%n,newlbfgs%dgrad(cp,:),1,newlbfgs%store(:),1)
    beta = newlbfgs%rho(cp)*yr
    beta = newlbfgs%alpha(cp) - beta
    !CALL DAXPY(newlbfgs%n,beta,newlbfgs%step(cp,:),1,newlbfgs%store(:),1)
    newlbfgs%store(:)=newlbfgs%store(:)+beta*newlbfgs%step(cp,:)
    cp = cp + 1
    IF ( cp>newlbfgs%m ) cp = 1
  END DO
  !
  ! STORE THE NEW SEARCH DIRECTION
  newlbfgs%step(newlbfgs%point,:) = newlbfgs%store(:)

  DO i = 1 , newlbfgs%n
    newlbfgs%store(i) = g(i)
  END DO

  stp = 1.D0

  ! ====================================================================
  ! check that function is descending
  ! ====================================================================
  yr= ddot (newlbfgs%n,newlbfgs%step(newlbfgs%point,:),1,g(:),1)
  if(yr>0.D0) then
    if(printl>=4) write(stdout,*) "Inverting BFGS direction"
    stp=-stp
  end if

  if(dbg) then
    print*,"stp",stp
    print*,"newlbfgs%step(:,:)",newlbfgs%step(:,:)
    print*,"newlbfgs%point",newlbfgs%point
    print*,"newlbfgs%step(newlbfgs%point,:)",newlbfgs%step(newlbfgs%point,:)
  end if

  ! set step (return value)
  step_(:)=stp*newlbfgs%step(newlbfgs%point,:)
  
  ! Store old gradient and coordinates
  newlbfgs%store(:) = g(:)
  newlbfgs%store2(:)= x(:)

END SUBROUTINE DLF_newlbfgs_STEP
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* newlbfgs/dlf_newlbfgs_restart
!!
!! FUNCTION
!!
!! Restart the L-BFGS algorithm, reset the memory
!!
!! SYNOPSIS
subroutine dlf_newlbfgs_restart
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: stdout, printl
  USE newlbfgs_module
  implicit none
  ! **********************************************************************
  if(.not.tinit) call dlf_fail("LBFGS not initialised!")
  if(.not.newlbfgs%tinit) then
    print*,"Instance of L-BFGS:",trim(newlbfgs%tag)
    call dlf_fail("This instance of LBFGS not initialised!")
  end if

  if(printl>=4) then
    if(trim(newlbfgs%tag)=="main") then
      write(stdout,'("Restarting L-BFGS optimiser")')
    else
      write(stdout,'("Restarting L-BFGS optimiser, instance: ",a)') &
          trim(newlbfgs%tag)
    end if
  end if
  newlbfgs%iter=0
end subroutine dlf_newlbfgs_restart
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* newlbfgs/dlf_newlbfgs_select
!!
!! FUNCTION
!!
!! Select an instance of L-BFGS 
!!
!! To be called before initialisation of new instances, but not before
!! the first instance, which is automaticall called "main".
!!
!! SYNOPSIS
subroutine dlf_newlbfgs_select(tag,newinstance)
!! SOURCE
  USE newlbfgs_module
  implicit none
  character(*), intent(in)  :: tag
  logical     , intent(in)  :: newinstance
  ! **********************************************************************
  if(.not.tinit) then
    !call dlf_fail("newlbfgs not initialised in newlbfgs_select!")
    call dlf_newlbfgs_init(1,1) ! initialise a dummy first instance
    newlbfgs%tinit=.false.
  end if
    

  ! try to selcet instance with %tag=tag
  newlbfgs=>newlbfgs_first
  do while (associated(newlbfgs%next))
    if(trim(newlbfgs%tag)==trim(tag)) exit
    newlbfgs=>newlbfgs%next
  end do

  if(newinstance) then
    if(trim(newlbfgs%tag)==trim(tag)) then
      print*,"Error, instance ",trim(tag)," already exists and selected &
          &with flag 'new'"
      call dlf_fail("Error selecting new hdlcopt instance")
    end if
    newtag=tag
    ! last instance selected, newtag contains name of new instance
  else
    if(newlbfgs%tag/=tag) then
      print*,"Error, instance ",trim(tag)," does not exist"
      print*,"Existing inctances:"
      newlbfgs=>newlbfgs_first
      do while (associated(newlbfgs))
        print*,"--",trim(newlbfgs%tag),"--"
        newlbfgs=>newlbfgs%next
      end do
      call dlf_fail("Error selecting new hdlcopt instance")

    end if
    ! instance with %tag = tag selected
  end if
  if(dbg) PRINT*,"SELECTED -",trim(tag),"-"
end subroutine dlf_newlbfgs_select
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* newlbfgs/dlf_newlbfgs_deselect
!!
!! FUNCTION
!!
!! Select first instance of L-BFGS 
!!
!! SYNOPSIS
subroutine dlf_newlbfgs_deselect
!! SOURCE
  USE newlbfgs_module
  implicit none
  ! **********************************************************************
  ! do nothing in case newlbfgs does not exist
  if(.not.tinit) return !call dlf_fail("newlbfgs not initialised!")
  call dlf_newlbfgs_select("main",.false.)
  newtag="main"
end subroutine dlf_newlbfgs_deselect
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* newlbfgs/dlf_newlbfgs_exists
!!
!! FUNCTION
!!
!! Find out if an instance with a specified tag is initialised
!!
!! SYNOPSIS
subroutine dlf_newlbfgs_exists(tag,exists)
!! SOURCE
  USE newlbfgs_module
  implicit none
  character(*),   intent(in)  :: tag
  logical     ,   intent(out) :: exists
  type(newlbfgs_type),pointer    :: newlbfgs_search
  ! **********************************************************************
  if(.not.tinit) then
    exists=.false.
    return
  end if

  newlbfgs_search=>newlbfgs_first
  exists=.false.
  do while (associated(newlbfgs_search))
    if(trim(newlbfgs_search%tag)==trim(tag)) then
      if(newlbfgs_search%tinit) then
        exists=.true.
        return
      end if
    end if
    newlbfgs_search => newlbfgs_search%next
  end do
end subroutine dlf_newlbfgs_exists
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* newlbfgs/dlf_newlbfgs_init
!!
!! FUNCTION
!!
!! Initialise the L-BFGS routines, allocate memory
!!
!! SYNOPSIS
subroutine dlf_newlbfgs_init(nvar,nmem)
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: stderr
  USE newlbfgs_module
  use dlf_allocate, only: allocate, deallocate
  implicit none
  integer,  intent(in)    :: Nvar ! number of variables
  integer,  intent(in)    :: Nmem ! number steps to remember
  ! **********************************************************************

  if(tinit) then
     ! this is initialisation of a new (not first) instance

    ! trap initialisation of "main" if another instance is allready there
    if(trim(newtag)=="main") then
      if(trim(newlbfgs%tag)/="main".or.newlbfgs%n/=1) then
        call dlf_fail("L-BFGS main instance is allready initialised")
      end if
      ! deallocate 
      call deallocate(newlbfgs%store)
      call deallocate(newlbfgs%store2)
      call deallocate(newlbfgs%rho)
      call deallocate(newlbfgs%alpha)
      call deallocate(newlbfgs%step)
      call deallocate(newlbfgs%dgrad)
      
    else

      ! newlbfgs should now point to the last existing instance
      if(dbg) print*,"Current lbfgs instance: ",trim(newlbfgs%tag)
      
      ! check that no instance with tag=newtag exists
      newlbfgs=>newlbfgs_first
      do while (associated(newlbfgs%next))
        newlbfgs=>newlbfgs%next
        if(trim(newlbfgs%tag)==trim(newtag)) then
          print*,"Instance with name ",trim(newtag)," already initialised"
          call dlf_fail("Instance with name already initialised")
        end if
      end do
      
      allocate(newlbfgs%next)
      newlbfgs=>newlbfgs%next
      nullify(newlbfgs%next)
      
    end if ! (trim(newtag)=="main")

  else
    ! this is the initialisation of the first instance
    tinit=.true.
    newtag="main"
    !if(associated(newlbfgs)) call dlf_fail("This instance of newlbfgs has already been initialised")
    ! allocate the newlbfgs pointer
    allocate(newlbfgs)

    nullify(newlbfgs%next)
    newlbfgs_first => newlbfgs

  end if


  newlbfgs%tag=newtag
  newlbfgs%tinit=.true.

  if(dbg) print*,"Allocating ",trim(newlbfgs%tag)
  newlbfgs%n=nvar
  newlbfgs%m=nmem
  newlbfgs%tprecon=.false.
  if(newlbfgs%n<=0) call dlf_fail("nvar in L-BFGS has to be > 0")
  if(newlbfgs%m<=0) call dlf_fail("Nmem in L-BFGS has to be > 0")

  ! allocate memory
  call allocate(newlbfgs%store,newlbfgs%N)
  call allocate(newlbfgs%store2,newlbfgs%N)
  call allocate(newlbfgs%rho,newlbfgs%M)
  call allocate(newlbfgs%alpha,newlbfgs%M)
  call allocate(newlbfgs%step,newlbfgs%M,newlbfgs%N)
  call allocate(newlbfgs%dgrad,newlbfgs%M,newlbfgs%N)

  ! variables to set at the beginning
  newlbfgs%iter = 0

  ! initialise (mainly to avoid NaNs in checkpointing)
  newlbfgs%store(:)=0.D0
  newlbfgs%store2(:)=0.D0
  newlbfgs%rho(:)=0.D0
  newlbfgs%alpha(:)=0.D0
  newlbfgs%step(:,:)=0.D0
  newlbfgs%dgrad(:,:)=0.D0

end subroutine dlf_newlbfgs_init
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* newlbfgs/dlf_newlbfgs_destroy
!!
!! FUNCTION
!!
!! Deallocate memory
!!
!! SYNOPSIS
subroutine dlf_newlbfgs_destroy
!! SOURCE
  use dlf_parameter_module, only: rk
  use dlf_global, only: stderr
  USE newlbfgs_module
  use dlf_allocate, only: deallocate
  implicit none
  logical         :: allgone
  ! **********************************************************************
  if(.not.tinit) call dlf_fail("LBFGS not initialised in dlf_newlbfgs_destroy!")

  !deallocate memory of the present instance
  call deallocate(newlbfgs%store)
  call deallocate(newlbfgs%store2)
  call deallocate(newlbfgs%rho)
  call deallocate(newlbfgs%alpha)
  call deallocate(newlbfgs%step)
  call deallocate(newlbfgs%dgrad)
  if(newlbfgs%tprecon) then
    deallocate(newlbfgs%precon)
  end if
  newlbfgs%tinit=.false.
!print*,"Destroying ",trim(newlbfgs%tag)

  ! check if all instances are deleted
  allgone=.true.
  newlbfgs=>newlbfgs_first
  do while (associated(newlbfgs))
    if(newlbfgs%tinit) then
      allgone=.false.
      exit
    end if
    newlbfgs => newlbfgs%next
  end do
  ! this may leave newlbfgs pointing nowhere. This is fine as one cannot expect
  ! it to point somewhere usefull after dlf_newlbfgs_destroy


  if(allgone) then

    ! if only a dummy has been initialised for the instance MAIN, it may still be
    ! allocated, even though %tinit would be false. Deallocate in this case
    newlbfgs=>newlbfgs_first
    if(allocated(newlbfgs%store)) call deallocate(newlbfgs%store)
    if(allocated(newlbfgs%store2)) call deallocate(newlbfgs%store2)
    if(allocated(newlbfgs%rho)) call deallocate(newlbfgs%rho)
    if(allocated(newlbfgs%alpha)) call deallocate(newlbfgs%alpha)
    if(allocated(newlbfgs%step)) call deallocate(newlbfgs%step)
    if(allocated(newlbfgs%dgrad)) call deallocate(newlbfgs%dgrad)

    ! deallocate everything
    tinit=.false.
    newlbfgs=>newlbfgs_first
    do while (associated(newlbfgs%next))
      newlbfgs_first => newlbfgs
      newlbfgs => newlbfgs%next
!print*,"Finally deallocating ",trim(newlbfgs_first%tag)
      deallocate(newlbfgs_first)
    end do
!print*,"Finally deallocating ",trim(newlbfgs%tag)
    deallocate(newlbfgs)
    nullify(newlbfgs)
    nullify(newlbfgs_first)
  end if

end subroutine dlf_newlbfgs_destroy
!!****

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!!****f* newlbfgs/dlf_newlbfgs_precon
!!
!! FUNCTION
!!
!! Supply a precondition matrix, an estimate of the inverse Hessian
!! to the module
!! The space for the matrix is only allocated here.
!!
!! SYNOPSIS
subroutine dlf_newlbfgs_precon(precon)
!! SOURCE
  use dlf_parameter_module, only: rk
  !use dlf_global, only: glob,stderr
  USE newlbfgs_module
  use dlf_allocate, only: allocate
  implicit none
  real(rk)  ,intent(in):: precon(newlbfgs%N,newlbfgs%N)
  ! **********************************************************************
  if(.not.tinit) call dlf_fail("LBFGS not initialised in newlbfgs_precon!")
  if(.not.newlbfgs%tinit) then
    print*,"Instance of L-BFGS:",trim(newlbfgs%tag)
    call dlf_fail("This instance of newlbfgs not initialised!")
  end if
  if(.not.newlbfgs%tprecon) then
    newlbfgs%tprecon=.true.
    call allocate(newlbfgs%precon, newlbfgs%N, newlbfgs%N)
  end if
  newlbfgs%precon=precon
end subroutine dlf_newlbfgs_precon
!!****

!   ----------------------------------------------------------
!   local routine, only to be used if no external ddot is
!   available (which is not recommended!)
FUNCTION DDOT_internal(n,dx,incx,dy,incy)
  USE dlf_parameter_module, only: rk                        
  IMPLICIT NONE
  !
  ! Dummy arguments
  !
  INTEGER :: incx , incy , n
  REAL(RK) :: DDOT_internal
  REAL(RK) , DIMENSION(n) :: dx , dy
  INTENT (IN) dx , dy , incx , incy , n
  !
  ! Local variables
  !
  REAL(RK) :: dtemp
  INTEGER :: i , ix , iy , m , mp1
  !
  !     forms the dot product of two vectors.
  !     uses unrolled loops for increments equal to one.
  !     jack dongarra, linpack, 3/11/78.
  !
  !
  DDOT_internal = 0.0D0
  dtemp = 0.0D0
  IF ( n<=0 ) RETURN
  IF ( incx==1 .AND. incy==1 ) THEN
    !
    !        code for both increments equal to 1
    !
    !
    !        clean-up loop
    !
    m = MOD(n,5)
    IF ( m/=0 ) THEN
      DO i = 1 , m
        dtemp = dtemp + dx(i)*dy(i)
      END DO
      IF ( n<5 ) THEN
        DDOT_internal = dtemp
        return
      END IF
    END IF
    mp1 = m + 1
    DO i = mp1 , n , 5
      dtemp = dtemp + dx(i)*dy(i) + dx(i+1)*dy(i+1) + dx(i+2)     &
          & *dy(i+2) + dx(i+3)*dy(i+3) + dx(i+4)*dy(i+4)
    END DO
    DDOT_internal = dtemp
  ELSE
    !
    !        code for unequal increments or equal increments
    !          not equal to 1
    !
    ix = 1
    iy = 1
    IF ( incx<0 ) ix = (-n+1)*incx + 1
    IF ( incy<0 ) iy = (-n+1)*incy + 1
    DO i = 1 , n
      dtemp = dtemp + dx(ix)*dy(iy)
      ix = ix + incx
      iy = iy + incy
    END DO
    DDOT_internal = dtemp
    RETURN
  END IF
END FUNCTION DDOT_INTERNAL

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_checkpoint_newlbfgs_write
  use dlf_parameter_module, only: rk
  use dlf_global, only: stderr
  USE newlbfgs_module
  use dlf_checkpoint, only: tchkform,write_separator
  implicit none
  type(newlbfgs_type),pointer :: newlbfgs_current
  ! **********************************************************************
  if(.not.tinit) call dlf_fail("LBFGS not initialised! (in checkpoint write)")

  newlbfgs_current => newlbfgs
  newlbfgs => newlbfgs_first

  if(tchkform) then
    open(unit=100,file="dlf_newlbfgs.chk",form="formatted")
    call write_separator(100,"current")
    write(100,*) newlbfgs_current%tag
    do while (associated(newlbfgs))
      call write_separator(100,"NM")
      write(100,*) newlbfgs%n,newlbfgs%m
      call write_separator(100,"Arrays")
      write(100,*) newlbfgs%store,newlbfgs%store2,newlbfgs%rho,newlbfgs%alpha,newlbfgs%step,newlbfgs%dgrad
      call write_separator(100,"Position")
      write(100,*) newlbfgs%point,newlbfgs%iter
      newlbfgs=>newlbfgs%next
    end do
    call write_separator(100,"END")
  else
    open(unit=100,file="dlf_newlbfgs.chk",form="unformatted")
    call write_separator(100,"current")
    write(100) newlbfgs_current%tag
    do while (associated(newlbfgs))
      call write_separator(100,"NM")
      write(100) newlbfgs%n,newlbfgs%m
      call write_separator(100,"Arrays")
      write(100) newlbfgs%store,newlbfgs%store2,newlbfgs%rho,newlbfgs%alpha,newlbfgs%step,newlbfgs%dgrad
      call write_separator(100,"Position")
      write(100) newlbfgs%point,newlbfgs%iter
      newlbfgs=>newlbfgs%next
    end do
    call write_separator(100,"END")
  end if
  close(100)
  newlbfgs => newlbfgs_current

end subroutine dlf_checkpoint_newlbfgs_write

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subroutine dlf_checkpoint_newlbfgs_read(tok)
  use dlf_global, only: stdout,printl
  USE newlbfgs_module
  use dlf_checkpoint, only: tchkform, read_separator
  implicit none
  logical,intent(out) :: tok
  logical             :: tchk
  integer             :: n_f,m_f
  character(40)       :: tag_read
  ! **********************************************************************
  tok=.false.
  if(.not.tinit) call dlf_fail("LBFGS not initialised! (in checkpoint read)")

  ! check if checkpoint file exists
  INQUIRE(FILE="dlf_newlbfgs.chk",EXIST=tchk)
  if(.not.tchk) then
    write(stdout,10) "File dlf_newlbfgs.chk not found"
    return
  end if

  if(tchkform) then
    open(unit=100,file="dlf_newlbfgs.chk",form="formatted")
  else
    open(unit=100,file="dlf_newlbfgs.chk",form="unformatted")
  end if

  newlbfgs => newlbfgs_first

  call read_separator(100,"current",tchk)
  if(.not.tchk) return    

  if(tchkform) then
    read(100,*,end=201,err=200) tag_read
  else
    read(100,end=201,err=200) tag_read
  end if

  do while (associated(newlbfgs))


    call read_separator(100,"NM",tchk)
    if(.not.tchk) return    

    if(tchkform) then
      read(100,*,end=201,err=200) n_f,m_f
    else
      read(100,end=201,err=200) n_f,m_f
    end if
    
    if(n_f/=newlbfgs%n) then
      write(stdout,10) "Different L-BFGS system size"
      close(100)
      return
    end if
    if(m_f/=newlbfgs%m) then
      write(stdout,10) "Different L-BFGS memory size"
      close(100)
      return
    end if
    
    call read_separator(100,"Arrays",tchk)
    if(.not.tchk) return 
    
    if(tchkform) then
      read(100,*,end=201,err=200) newlbfgs%store,newlbfgs%store2,newlbfgs%rho,newlbfgs%alpha,newlbfgs%step,newlbfgs%dgrad
    else
      read(100,end=201,err=200) newlbfgs%store,newlbfgs%store2,newlbfgs%rho,newlbfgs%alpha,newlbfgs%step,newlbfgs%dgrad
    end if
    
    call read_separator(100,"Position",tchk)
    if(.not.tchk) return 

    if(tchkform) then
      read(100,*,end=201,err=200) newlbfgs%point,newlbfgs%iter
    else
      read(100,end=201,err=200) newlbfgs%point,newlbfgs%iter
    end if

    newlbfgs=>newlbfgs%next

  end do ! while (associated(newlbfgs))

  call read_separator(100,"END",tchk)
  if(.not.tchk) return 
    
  ! now make sure that the instance of the checkpoint is selected
  call dlf_newlbfgs_select(trim(tag_read),.false.)

  if(printl >= 6) write(stdout,"('LBFGS checkpoint file sucessfully read')")
  close(100)
  tok=.true.

  return

  ! return on error
  close(100)
200 continue
  write(stdout,10) "Error reading LBFGS checkpoint file"
  return
201 continue
  write(stdout,10) "Error (EOF) reading file"
  return

10 format("Checkpoint reading WARNING: ",a)

end subroutine dlf_checkpoint_newlbfgs_read
