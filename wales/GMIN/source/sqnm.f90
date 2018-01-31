!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2017 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!                   SUBROUTINE SLBFGS
!Written by jk669 Feb 2017, based on code from mymylbfgs.f
!
!   Stabilized L-BFGS algorithm (Schaefer et al) added 2017:
!   Stabilized quasi-Newton optimization of noisy potential energy surfaces
!	J. Chem. Phys. 142. 034112 (2015) doi:10.1063/1.4905665
!      
!This subroutine is called from mylbfgs.f90 only.
!Requires LAPACK functions: dsyev, dsysv
!Dependencies: Requires calls to potential and topology_reader, and uses some global variables from Commons
!
module sqnm_func 
    implicit none
contains
    !dot product
    pure double precision function dot(numcoords, v1, v2)
        implicit none
        integer, intent(in) :: numcoords
        double precision, intent(in) :: v1(numcoords), v2(numcoords)
        double precision :: temp
        integer :: i
        temp=0.0d0
        do i=1, numcoords
            temp=temp+v1(i)*v2(i)
        end do 
        dot=temp
    end function dot

    !vector norm
    pure double precision function norm(numcoords, v1)
        implicit none
        integer, intent(in) :: numcoords
        double precision, intent(in) :: v1(numcoords)
        norm = sqrt(dot(numcoords,v1,v1))
    end function norm
end module sqnm_func

subroutine sqnm(numcoords,xcoords,maxrmsgrad,nhistmax,converget,energy,itermax,iter)
    use COMMONS, only: MAXERISE,SQNM_DEBUGT,AMBERT,AMBER12T,BONDS,SQNM_DEBUGRUN,SQNM_DEBUGLEVEL,SQNM_WRITEMAX,SQNM_BIOT,RMS
    use sqnm_func
    implicit none
    !arguments
    integer, intent(in)             :: numcoords 
    double precision, intent(inout) :: xcoords(numcoords) 
    double precision, intent(in)    :: maxrmsgrad
    integer, intent(in)             :: nhistmax 
    logical, intent(inout)          :: converget 
    double precision, intent(inout) :: energy 
    integer, intent(in)             :: itermax 
    integer, intent(inout)          :: iter

!The subroutine is invoked as 
!   slbfgs(numcoords,xcoords,maxrmsgrad,converget,energy,itermax,iter)
!where
!   numcoords is a positive integer equal to the number of coordinates. This should be equal to NATOMS*3. It is the responisibility
!           of the calling function to verify this. 
!   xcoords is a vector of length numcoords containing the coordinates of the initial position on input, and containing the 
!           optimized position upon exit if converget==.true., or the best positon obtained thus far if converget==.false.
!   maxrmsgrad is a cut-off for the gradient. The rms of the gradient is compared to maxrmsgrad to determine if the algorithm
!           has found a local minimum. The algorithm is converged if rms(grad)<maxrmsgrad.
!   nhistmax is a positive integer containing the maximum number of historic position and gradient vectors that will be stored for  
!           determination of teh significant subspace.
!   converget is a logical which should be false on input (since the system is not converged) and is false on output UNLESS the 
!           gradient at the new position is less than maxrmsgrad.
!   energy is a single double precision value that contains the energy of the function at the initial position on input, and the 
!           energy of the optimal position on output. 
!   itermax is a positive integer containing the maximum number of minimization steps that may be performed during the algorithm.
!   iter is the actual number of steps performed during the minimization algorithm. Note that any preliminary steepest descent
!           steps used to move the initial position closer to a minimum are NOT included in this count. 
!
    !**********parameters**********
    double precision, parameter :: alphascale=1.1d0 !preconditioned gradient scaling
    double precision, parameter :: minalpha=1.0d-15 !used to prevent alpha from getting too small, which may cause underflow
    double precision, parameter :: initsteepstepsize=1.0d-3 !initial steepest descent stepsize
    double precision, parameter :: cutoffratio=1.0d-4 !cutoff ratio for small eigenvalues for finding significant subspace
    double precision, parameter :: anglecutoff=0.2d0 !cos(angle) cutoff for changing preconditioning with factor alpha 
    integer, parameter :: maxdescents=15 !maximum number of steepest descent steps allowed
    integer, parameter :: maxeeval=100 !maximum number potential calls allowed in a line search
    !**********local variables**********
    integer :: i, j !counters
    integer :: nhist=0; !current history length
    double precision :: grad(numcoords) !current gradient
    double precision :: ediff !difference between new and old energy
    double precision :: pregrad(numcoords) !preconditioned gradient in subspace G
    double precision :: gradhist(numcoords,nhistmax) !history of gradients
    double precision :: xhist(numcoords,nhistmax) !history of gradients
    integer :: lwork !used for LAPACK function dsyev
    integer :: ndim !number of dimensions of the significant subspace. 
    integer, save :: nbonds !number of bonds in the molecule (Used for biomolecules)
    integer :: eeval ! number of potential calls in the linesearch 
    double precision :: maxlambda !the value of 1/(maximum eigenvalue of the overlap matrix)
    double precision :: tempval !exactly what you think.
    double precision :: gradperp(numcoords) !portion of gradient perpendicular to preconditioned gradient
    double precision :: newxcoords(numcoords) !newly calculated x-coordinates
    double precision :: oldxcoords(numcoords) !previous x-coordinates
    double precision :: oldgrad(numcoords) !previous gradient.
    double precision :: newgrad(numcoords) ! gradient of energy at newly calculated x-coordinates
    double precision :: residuetemp(numcoords) !temporary vector for calculating residues
    double precision :: invproj(numcoords,numcoords) !inverse projection from subspace to find orthogonal complement
    double precision :: enew !energy at newly calculated x-coordinates
    double precision :: ethresh !Energy threshhold to determine if a step is accepted. 
    double precision :: alpha !preconditioned gradient scaling
    double precision :: angle !angle between gradient and preconditioned gradient
    double precision :: gradrms !Root mean squared of the gradient
    double precision :: steepstepsize !steepest descent step size
    logical :: acceptstep !Set based on whether the last step was accepted or not.
    logical :: bio = .false. !Whether the molecule being minimized is a biomolecule (true) or not (false).
    double precision :: gradstr(numcoords) !the bond-stretching component of the gradient for biomolecules
    integer :: atom1, atom2 !the atom numbers used for finding atom displacements for biomolecules 
    double precision :: atom_disp(3) !Displacement between atom1 and atom2 for biomolecules
    integer :: signflip !number of projections which have flipped since the last iteration for biomolecules
    double precision :: alpha_bio !steepest descent stepsize for biomolecule gradients 
    integer :: infoc !information integer from subroutines and LAPACK functions. 
    integer :: totalfunceval ! number of calls to the potential energy function.
    integer, save :: runs=0 !keeps track of how many times this subroutine is called during the program.
    logical :: debug_run !If set to true, this prints debug information for the run. 
    integer :: writemax !maximum number of lines of something to pring when debugging

    !allocatable local variables. All variables here have at least one dimension of length nhist, ndim, or nbond, which 
    !change or are unknown.
    double precision, allocatable, dimension(:,:) :: bondvect !sparce vectors containing atom distances of bonded atoms.
    double precision, allocatable, dimension(:)   :: oldcoeff !sparce vectors containing atom distances of bonded atoms.  
    double precision, allocatable, dimension(:)   :: bondcoeff !vector with coefficients of bond vectors to find gradstr
    double precision, allocatable, dimension(:,:) :: bondmat !Matrix containging bondvect^t bondvect
    double precision, allocatable, dimension(:)   :: ipiv !vector of length nbonds for dsysv
    double precision, allocatable, dimension(:,:) :: graddisp !gradient displacement vectors
    double precision, allocatable, dimension(:,:) :: xdisp !xcoord displacement vectors
    double precision, allocatable, dimension(:)   :: norms !the norms of the xdisp vectors
    double precision, allocatable, dimension(:,:) :: overlap !Overlap matrix of normalized displacement vectors
    double precision, allocatable, dimension(:)   :: work !used for LAPACK function dsyev
    double precision, allocatable, dimension(:)   :: eigenvalues !contains eigenvalues of overlap matrix or projected hessian
    double precision, allocatable, dimension(:,:) :: subspan !matrix with subspace-spanning vectors as columns
    double precision, allocatable, dimension(:,:) :: gradspan !matrix with curvature-describing vectors as columns
    double precision, allocatable, dimension(:,:) :: sigsubspan !matrix with significant subspace-spanning vectors as columns
    double precision, allocatable, dimension(:,:) :: siggradspan !matrix with significant curvature-describing vectors as columns
    double precision, allocatable, dimension(:,:) :: SO !approximation of shape operator
    double precision, allocatable, dimension(:,:) :: principals !holds principal directions of the projected hessian
    double precision, allocatable, dimension(:)   :: residues !contains residues calculated by Weinstein's criterion

    !increment number of runs to keep track of calls to the function during GMIN.
    runs=runs+1

    !Determine if run is being debugged
    debug_run=.false.
    if(SQNM_DEBUGT .and. (SQNM_DEBUGRUN==0 .or. SQNM_DEBUGRUN+1==runs)) then 
        debug_run=.true.
        writemax=SQNM_WRITEMAX
    end if 
    if((runs > SQNM_DEBUGRUN+1) .and. SQNM_DEBUGT ) STOP !stop after debug info is printed

    !check input for any potential issues (responisbility for correctness here lies with calling subroutine). 
    if (numcoords<=0 .or.itermax<=0 .or. maxrmsgrad<=0.0d0 ) then 
        if(debug_run .and. SQNM_DEBUGLEVEL>=0) write(*, '(a,i5,a)') &
            'WARNING: subroutine slbfgs in slbfgs.f90 received bad input on run ', runs, &
            '. No minimization occurred.'
        return
    end if

    if(debug_run .and. SQNM_DEBUGLEVEL>=0) write(*, '(a,i5,a)') 'STARTING RUN :', runs, '.'

    !initialize variables.
    converget=.false.
    iter=1
    nhist=0
    ethresh=MAXERISE !from COMMONS
    lwork=4*numcoords+64 !size of work vector for LAPACK function dsyev
    totalfunceval=0
    call potential(xcoords,grad,energy,numcoords,.FALSE.)
    if(energy==1.0d6) return
    allocate(work(lwork))
    steepstepsize=min(initsteepstepsize,1.0d0/maxval(abs(grad)))
    totalfunceval=totalfunceval+1 !update number of energy function calls
    if (debug_run .and. SQNM_DEBUGLEVEL>=0) then 
        write(*, '(A)') 'xcoords and grad before steepest descent steps'
        do i=1,min(writemax,numcoords)
            write(*, '(2(f20.15))') xcoords(i), grad(i)
        end do 
    end if
    !In GMIN, atoms are moved randomly, and initial gradients may be very large. Do steepest descent steps to reduce energy.
    !perform at least one steepest descent step, and store more recent previous position and gradient in the histories.
    !Steepest descent is performed with a constant stepsize to prevent numerous function calls for linesearch. 
    steepest_descent: do while (iter<maxdescents)
        if (debug_run .and. SQNM_DEBUGLEVEL>=1) then
            print *, 'Iter: ', iter, ' Init E: ', energy, '. Init gradrms: ', norm(numcoords, grad), 'stepsize', steepstepsize 
        end if
        eeval=0
        ! store position and gradient in history before linsearch, since this subroutine changes values. 
        xhist(:,1)=xcoords(:)
        gradhist(:,1)=grad(:)
        !update the stepsize. Try a larger step if possible, or way larger if step is very small.
        if (steepstepsize<1.0d-8 .and. iter>=3) then 
            steepstepsize=5.0D-2
        else 
            steepstepsize=alphascale*steepstepsize 
        end if 
        tempval=energy !store energy value to ensure it's lower
        call linesearch(numcoords,xcoords,energy,grad,-grad,steepstepsize,ethresh,maxeeval,infoc,eeval)
        totalfunceval=totalfunceval+eeval !update number of energy function calls
        gradrms=norm(numcoords,grad)
        if (RMS<maxrmsgrad) then !we already converged! Great work, everyone.
            converget=.true.
            if(debug_run .and. SQNM_DEBUGLEVEL>=0 ) then
                write(*, '(a, i4,a,i4)') 'Run ', runs, ' has converged during steepest descent, iteration ', iter
                write(*,'(a,f12.8,a,f12.6)') 'Final rms of gradient = ', gradrms, ', Final energy: ', energy
                print *, 'Total potential energy function evaluations: ', totalfunceval
            end if
            return
        else if (tempval+ethresh<energy.or.(gradrms>1.0d10.and.iter>=itermax/2).or.(gradrms>1.0d15.and.iter>=2)) then 
            !bad things have happened; the energy has increased after the stepsearch.
            !the gradient was probably way too massive coming in, and a small enough step could not be taken.
            !Escape to prevent overflow. 
            if(debug_run .and. SQNM_DEBUGLEVEL>=0 ) then
                print *, 'ERROR: Run did not converge due to stepsearch being unable to locate a minimum. '
                write(*, '(a, i4,a,i4)') ' Stepsearch terminated on run ', runs, ', stepsearch iteration ', iter
                write(*,'(a,f12.8,a,f12.6)') 'RMS of gradient = ', gradrms, ', Energy: ', energy
                print *, 'Total potential energy function evaluations: ', totalfunceval
            end if
            return
        end if
        iter=iter+1
    end do steepest_descent

    acceptstep=.true. !since we did a steepest descent step, we had a successful step, and it is stored in history
    nhist=1 !since previous position stored in steepest descent portion
    alpha=min(initsteepstepsize,1.0d0/maxval(abs(grad))) !works well emperically for LJ clusters and oligopeptides
    if (debug_run .and. SQNM_DEBUGLEVEL>=1) then
        write(*, '(a, i4, a, f15.8, a, f15.8)') 'Number of steepest descent steps taken: ', iter, & 
            '. Energy after Steepest Descents: ', energy, '. Gradient RMS after  Steepest Descents: ', gradrms
    end if
    if (converget) then
        if(debug_run .and. SQNM_DEBUGLEVEL>=0) write(*, '(a,i4,a,i3,a)') 'Run ', runs, ' converged after ', iter,&
         ' steps. No slbfgs needed.'
        return
    end if 
    if (debug_run .and. SQNM_DEBUGLEVEL>=1) then 
        write(*, '(A)') 'xcoords and grad after steepest descent steps'
        do i=1,min(writemax,numcoords)
            write(*, '(2(f20.15))') xcoords(i), grad(i)
        end do 
    end if
    !if we are using an AMBER force field, we will assume that we are dealing with a biomolecule, and use the 
    !procedure for biomolecules. For this, we need the bond information. Also make sure that the bio portion was not 
    ! turned off. 
    biocheck: if((AMBERT .or. AMBER12T) .and. (SQNM_BIOT)) then 
        !this function call sets nbonds to the number of bonds, and then the global variable BONDS 
        !is a nbonds*2 integer array containing pairs of bonded atoms. 
        bio = .true.
        if (.not. allocated(BONDS)) then 
            call topology_reader(nbonds)
        end if 
        if (.not. allocated(bondvect)) then
            allocate(bondvect(numcoords,nbonds))
            allocate(bondcoeff(nbonds))
            allocate(oldcoeff(nbonds))
            allocate(bondmat(nbonds,nbonds))
        end if 
        bondcoeff(:)=0.0d0
        !set bond steepest descent stepsize. 
        alpha_bio = steepstepsize 
        if (debug_run .and. SQNM_DEBUGLEVEL>=1) then 
            write(*,'(a,i5,a)') 'The molecule has been identified as a biomolecule. ', nbonds, ' bonds were found.'
            write(*,'(a,i5,a)') 'The molecule has has ', numcoords/3, ' atoms.'
        end if 
    else if ((AMBERT .or. AMBER12T) .and. .not. (SQNM_BIOT)) then biocheck
        if (debug_run .and. SQNM_DEBUGLEVEL>=1) then 
            write(*,'(2a)') 'The minimization of this molecule is using the AMBER potential but is not being treated as a ', &
            'biomolcule due to the presence of the BIOOFF keyword in the data file.'
        end if 
    end if biocheck

!********************************************************************************************************
!Main slbfgs do loop. This is broken down into numerous parts:                                          *
!    If we have a biomolecule,                                                                          *
!       Compute bond stretching components.                                                             *
!       Find bond vector coefficients                                                                   *
!       Adjust alpha_bio                                                                                *
!       Adjust gradient and position based on the stretching components.                                *
!    Find preconditioned gradient.                                                                      *
!       find x-displacements and grad displacements from histories.                                     *
!           find x-displacement norms and normalize x-displacement vectors.                             *
!       find overlap matrix.                                                                            *
!           get eigenvectors and eigenvalues of overlap matrix                                          *
!           find number of dimensions of significant subspace based on eigenvalues                      *
!       find significant subspace of tangent space on energy manifold at position                       *
!           calculate spanning vectors of subspace                                                      *
!           determine significance of vectors from overlap eigenvalues                                  *
!       obtain curvature information from significant subspace                                          *
!           find 'normalized' gradients in tangent subspace                                             *
!           produce symmetrized shape operator                                                          *
!           project principal directions of shape operator onto energy manifold                         *
!       correct curvatures due to deviations of calculated shape operator from second fundamental form  *
!           calculate Weinstein residues                                                                *
!           adjust curvatures                                                                           *
!       precondition the gradient                                                                       *
!           find rejection of gradient from projection onto preconditioned gradient vector              *
!           adjust rejection by factor of alpha                                                         *
!    x(iter+1)=x-preconditioned gradient.                                                               *
!    Find energy at x(iter+1), and determine if minimization occurred.                                  *
!       If energy decreased or alpha < alphastart/10, accept step                                       *
!           update x-coords and gradient histories                                                      *
!           if nhist<nhistmax, nhist++                                                                  *
!           adjust alpha                                                                                *
!       Else, reject step (x(iter+1)=x(iter)                                                            *
!           alpha=alpha/2                                                                               *
!           purge history                                                                               *
!   iter++                                                                                              *
!Loop until converged or maxiter is reached.                                                            *
!********************************************************************************************************
    mainloop: do while (.not. converget .and. iter .le. itermax)
        if(debug_run .and. SQNM_DEBUGLEVEL>=0) then 
            write(*, '(a, i4, a, f15.8, a, f15.8)') 'Number of iterations done: ', iter, & 
            '. Energy: ', energy, '. Gradient RMS: ', gradrms
            print *, 'Bio On:',SQNM_BIOT,'alpha_bio: ', alpha_bio, ' bio: ',bio,' acceptstep: ',acceptstep,' nhist : ',nhist
            print *, 'run: ',runs,' iter: ',iter,' nhist: ',nhist,' Energy: ',energy,' RMSgrad: ',gradrms,' alpha: ',alpha 
            write(*, '(A)') 'xcoords and grad at loop start'
            do i=1,min(writemax,numcoords)
                write(*, '(2(f24.20))') xcoords(i), grad(i)
            end do 
        end if 

        oldgrad(:)=grad(:)
        oldxcoords(:)=xcoords(:)
        !FOR BIOMOLECULES, WE MUST FIRST FIND BOND SRETCHING COMPONENTS
        bio_optim: if (bio .and. acceptstep) then
            !save old positions in case new ones are higher in energy. 
            !FIND ATOM DISTANCES FOR BOND VECTORS 
            do i=1,nbonds
                !get bonded atom distances.
                !bondvect is sparce with six elements in each column, corresponding to the bonded atom displacements. 
                !Because displacement is antisymmetric, each column sum should be 0.
                atom1=BONDS(i,1)
                atom2=BONDS(i,2)
                !displacement is xyz(atom2)-xyz(atom1)
                atom_disp(:)=xcoords(3*atom2-2:3*atom2)-xcoords(3*atom1-2:3*atom1)
                bondvect(3*atom1-2:3*atom1,i)=atom_disp(:)
                bondvect(3*atom2-2:3*atom2,i)=-atom_disp(:) !by antisymmetry
            end do
            if(debug_run .and. SQNM_DEBUGLEVEL>=1) then
                print *, 'bond vectors'
                do i=1,min(writemax,nbonds) 
                    write(*, '(f13.9)', advance='no') (bondvect(i,j),j=1,min(writemax,nbonds))
                    write(*,*) ""
                end do
            end if 
            !FIND BOND VECTOR COEFFICIENTS 
            !This is done by solving (bondvect^t * bondvect)*coeffs = bondvect^t * grad
            !find bondgrad=bondvect^t * grad and bondmat=(bondvect^t * bondvect)
            do i=1,nbonds
                bondcoeff(i)=dot(numcoords,bondvect(:,i),grad(:))
                bondmat(i,i)=dot(numcoords,bondvect(:,i),bondvect(:,i))
                do j=i+1,nbonds
                    !compute only lower half since matrix is symmetric.
                    bondmat(j,i)=dot(numcoords,bondvect(:,i),bondvect(:,j))
                end do
            end do
            !count number of projections bondcoeff which have not changed sign from last time by comparing to the old coefficient
            signflip=0
            do i=1,nbonds
                if ((oldcoeff(i)>=0.0d0 .and. bondcoeff(i)>=0.0d0) .or. (oldcoeff(i)<0.0d0 .and. bondcoeff(i)<0.0d0)) then
                    signflip=signflip+1
                end if 
            end do
            oldcoeff(:)=bondcoeff(:)

            if(debug_run .and. SQNM_DEBUGLEVEL>=1) then
                print *, 'bond matrix'
                do i=1,min(writemax,nbonds) 
                    write(*, '(f13.9)', advance='no') (bondmat(i,j),j=1,min(writemax,nbonds))
                    write(*,*) ""
                end do
                print *, 'bond coeffs'
                do i=1,min(writemax,nbonds)
                    write(*,'(f16.10)') bondcoeff(i)
                end do
                print *, 'Number of sign flips: ', signflip, '; alpha_bio = ', alpha_bio
            end if 


            !Solve linear system Bondmat*coeff=bondvect. 
            !Note that a LAPACK function for symmetric matrices is called since we calculated the lower
            !triangular portion of bondmat. Solution should now be in bondcoeff. 
            allocate(ipiv(nbonds))
            call dsysv('L',nbonds,1,bondmat,nbonds,ipiv,bondcoeff,nbonds,work,lwork,infoc)
            deallocate(ipiv)
            if(debug_run .and. SQNM_DEBUGLEVEL>=1) then
                print *, 'bond coeffs AFTER solving'
                do i=1,min(writemax,nbonds) 
                    write(*,'(f16.10)') bondcoeff(i)
                end do
            end if 


            !FIND BOND STRETCHING COMPONENTS OF THE GRDAIENT, gradstr
            gradstr(:)=0
            do i=1,nbonds
                gradstr(:)=gradstr(:)+bondcoeff(i)*bondvect(:,i)
            end do

            !MINIMIZE ALONG BOND-STRETCHING COMPONENTS AND SEND REST OF GREADIENT TO PRECONDITIONING
            !adjust gradient to non-bonded portions 
            grad(:)=grad(:)-gradstr(:)
            !scale bonded gradient component steepest descent stepsize based on number of sign flips 
            if(signflip > 2*nbonds/3) then 
                alpha_bio=min(alpha_bio*alphascale,1.0d0)
            else
                alpha_bio=alpha_bio/alphascale
            end if
            !perform simple steepest descent step along bonded portions of gradient 
            xcoords(:)=xcoords(:)-alpha_bio*gradstr(:)
            if (debug_run .and. SQNM_DEBUGLEVEL>=1) then
                print *, 'New alpha_bio: ', alpha_bio
                write(*, '(A)') '   old xcoords         xcoords             old grad            new grad            grad_str'
                do i=1,min(writemax,numcoords)
                    write(*, '(5(f20.15))') oldxcoords(i), xcoords(i), oldgrad(i), grad(i), gradstr(i)
                end do 
            end if
        end if bio_optim


        !**********FINDING PRECONDITIONED GRADIENT**********
        !FIND COORDINATE AND GRADIANT DISPLACEMENTS  
        preconditioning: if (nhist==0) then !last step failed; perform a steepest descent-like step
            pregrad(:)=alpha*grad
        else preconditioning
            !FIND POSITION AND GRADIENT FINITE DIFFERENCES
            allocate(xdisp(numcoords,nhist))
            allocate(graddisp(numcoords,nhist)) 
            !for first column of displacements, use current position and then first historical position
            xdisp(:,1)=xcoords(:)-xhist(:,1)
            graddisp(:,1)=grad(:)-gradhist(:,1)
            !for remaining columns, just use historical positions.
            do i=2,nhist
                xdisp(:,i)=xhist(:,i-1)-xhist(:,i)
                graddisp(:,i)=gradhist(:,i-1)-gradhist(:,i)
            end do  

            if (debug_run .and. SQNM_DEBUGLEVEL>=2) then 
                write(*, '(A)') 'xhist'
                do i=1,min(writemax,numcoords)
                    write(*, '(f20.15)', advance='no') (xhist(i,j),j=1,nhist)
                    write(*,*) ""
                end do 
                write(*, '(A)') 'gradhist'
                do i=1,min(writemax,numcoords)
                    write(*, '(f20.15)', advance='no') (gradhist(i,j),j=1,nhist)
                    write(*,*) ""
                end do 
                write(*, '(A)') 'xdisp'
                do i=1,min(writemax,numcoords)
                    write(*, '(f16.10)', advance='no') (xdisp(i,j),j=1,nhist)
                    write(*,*) ""
                end do
                write(*, '(A)') 'graddisp'
                do i=1,min(writemax,numcoords)
                    write(*, '(f16.10)', advance='no') (graddisp(i,j),j=1,nhist)
                    write(*,*) ""
                end do   
            end if 

            !FIND NORMS AND NORMALIZE X-DISPLACEMENT VECTORS
            allocate(norms(nhist))
            do i=1, nhist 
                norms(i)=norm(numcoords,xdisp(:,i))
                if (norms(i)<1.0d-50) then !floating point errors caused underflow. Return to avoid division by 0.
                    if (debug_run  .and. SQNM_DEBUGLEVEL>=0) print *,'Norm too small for division on run ',runs,', iter,',iter,&
                        'hist ', nhist, ' Bombing out...'
                    converget=.false.
                    return
                end if
                !non-normalized x-displacements will not be needed again, so we can normalize in original matrix.
                xdisp(:,i)=xdisp(:,i)/norms(i)
            end do

            if (debug_run .and. SQNM_DEBUGLEVEL>=2) then 
                write(*, '(A)') 'displacement norms'
                do i=1,min(writemax,nhist)
                    write(*, '(f16.10)', advance='no') (norms(i))
                    write(*,*) ""
                end do 
                write(*, '(A)') 'xdisp - normalized'
                do i=1,min(writemax,numcoords)
                    write(*, '(f16.10)', advance='no') (xdisp(i,j),j=1,min(writemax,nhist))
                    write(*,*) ""
                end do  
            end if

            !FIND OVERLAP MATRIX
            allocate(overlap(nhist, nhist))
            do i=1, nhist
                !diagonal elements
                overlap(i,i)=dot(numcoords,xdisp(:,i),xdisp(:,i))
                do j=i+1, nhist
                    !compute only lower half since matrix is symmetric.
                    overlap(j,i)=dot(numcoords,xdisp(:,i),xdisp(:,j))
                end do
            end do

            if (debug_run .and. SQNM_DEBUGLEVEL>=3) then 
                write(*, '(A)') 'overlap'
                do i=1,nhist
                    write(*, '(f16.10)', advance='no') (overlap(i,j),j=1,nhist)
                    write(*,*) ""
                end do
            end if

            !GET EIGENVALUES AND EIGENVECTORS OF OVERLAP MATRIX
            !allocate vectors for eigenvectors and for work for dsyev
            allocate(eigenvalues(nhist))
            !dsyev will place eigenvalues of overlap into the eigenvalues vector in order of ascending value, 
            !and the columns of overlap will be replaces with the corresponding eigenvectors.
            !Note that the LAPACK function for symmetric matrices is called since we calculated the lower
            !triangular portion of overlap.
            call dsyev('V',"L",nhist,overlap,nhist,eigenvalues,work,lwork,infoc) 
            if (debug_run .and. SQNM_DEBUGLEVEL>=4) then 
                write(*, '(A)') 'eigenvalues'
                do i=1,nhist
                    write(*, '(f16.10)', advance='no') (eigenvalues(i))
                    write(*,*) ""
                end do
                write(*, '(A)') 'overlap matrix eigenvectors'
                do i=1,nhist
                    write(*, '(f16.10)', advance='no') (overlap(i,j),j=1,nhist)
                    write(*,*) ""
                end do  
            end if 

            !FIND THE DIMENSIONS OF THE SIGNIFICANT SUBSPACE. This is done by comparing the ratio of all eigenvalues to 
            !the largest eigenvalue. This is used to find the dimension of the significant subspace:
            ndim=1 !minimum dimension if nhist>0
            maxlambda=1.0D0/eigenvalues(nhist) !eigenvalues are in ascending order, so eigenvalues(nhist) is the largest.
            do i=1, nhist-1
                if (maxlambda*eigenvalues(i)>=cutoffratio) then !equivalent to checking eigenvalues(i)/eigenvalues(nhist)
                    ndim=ndim+1
                end if 
            end do

            !FIND SPANNING VECTORS OF SUBSPACE. ALSO FIND NORMALIZED GRADIENTS
            allocate(subspan(numcoords,nhist))
            allocate(gradspan(numcoords,nhist))
            do i = 1, nhist
                subspan(:,i)=0.0D0
                gradspan(:,i)=0.0D0
                do j=1, nhist
                    subspan(:,i)=subspan(:,i)+overlap(j,i)*xdisp(:,j)
                    gradspan(:,i)=gradspan(:,i)+overlap(j,i)*graddisp(:,j)/norms(j)
                end do
                !note that any eigenvalues should be at least zero since the overlap matrix is positive semi-definite; 
                !however, due to floating point error, we may have very slightly negative or zero eigenvalues. Take max(). 
                tempval=1.0D0/sqrt(max(eigenvalues(i),1.0d-10))
                subspan(:,i)=subspan(:,i)*tempval
                gradspan(:,i)=gradspan(:,i)*tempval
            end do

            if (debug_run .and. SQNM_DEBUGLEVEL>=4) then 
                write(*, '(A)') 'subspance span'
                do i=1,min(writemax,numcoords)
                   write(*, '(f16.10)', advance='no') (subspan(i,j),j=1,min(writemax,nhist))
                    write(*,*) ""
                end do 
                write(*, '(A)') 'subspance gradients'
                do i=1,min(writemax,numcoords)
                    write(*, '(f16.10)', advance='no') (gradspan(i,j),j=1,min(writemax,nhist))
                    write(*,*) ""
                end do 
            end if 

            !SELECT VECTORS AND GRADIENTS IN THE SIGNIFICANT SUBSPACE
            !there will be ndim of these vectors, out of nhist, with ndim<=nhist
            allocate(sigsubspan(numcoords,ndim))
            allocate(siggradspan(numcoords,ndim))
            !Note that eigenvalues were listed in ascending order before, so we reverse the direction of vector selection.
            !The first vectors correspond with the largest eigenvalue and so on.
            do i=1, ndim
                sigsubspan(:,i)=subspan(:,nhist-i+1)
                siggradspan(:,i)=gradspan(:,nhist-i+1)
            end do

            if (debug_run .and. SQNM_DEBUGLEVEL>=4) then 
                write(*, '(A)') 'significant subspance span'
                do i=1,min(writemax,numcoords)
                    write(*, '(f16.10)', advance='no') (sigsubspan(i,j),j=1,min(writemax,ndim))
                    write(*,*) ""
                end do 
                write(*, '(A)') 'significant subspance gradients'
                do i=1,min(writemax,numcoords)
                    write(*, '(f16.10)', advance='no') (siggradspan(i,j),j=1,min(writemax,ndim))
                    write(*,*) ""
                end do 
            end if 

            !PRODUCE SYMMETRIC SHAPE OPERATOR
            allocate(SO(ndim,ndim))
            SO(:,:)=0.0D0
            do i=1, ndim
                SO(i,i)=dot(numcoords,siggradspan(:,i),sigsubspan(:,i)) !diagonals of shape operator
                do j=i+1, ndim
                    !matrix is symmetric; find only lower triangular 
                    SO(j,i)=0.5D0*(dot(numcoords,siggradspan(:,i),sigsubspan(:,j))+ &
                        dot(numcoords,siggradspan(:,j),sigsubspan(:,i))) 
                end do
            end do

            if (debug_run .and. SQNM_DEBUGLEVEL>=5) then 
                write(*, '(A)') 'shape operator matrix'
                do i=1,min(writemax,ndim)
                    write(*, '(f16.10)', advance='no') (SO(i,j),j=1,min(writemax,ndim))
                    write(*,*) ""
                end do   
            end if 

            !PERFORM EIGENVALUE DECOMFC=ifort cmake -DCMAKE_BUILD_TYPE=Debug -DWITH_AMBER=1 ../../sourcePOSITION OF SHAPE OPERATOR
            !we no longer need the eigenvectors from the overlap matrix, so we reuse the vector.
            deallocate(eigenvalues)
            allocate(eigenvalues(ndim))
            !Note that the LAPACK function for symmetric matrices is called since we calculated the lower
            !triangular portion of the shape operator.
            call dsyev('V',"L",ndim,SO,ndim,eigenvalues,work,lwork,infoc)  

            if (debug_run .and. SQNM_DEBUGLEVEL>=5) then 
                write(*, '(A)') 'shape operator eigenvalues'
                do i=1,min(writemax,ndim)
                    write(*, '(f16.10)') (eigenvalues(i))
                end do
                write(*, '(A)') 'shape operator matrix eigenvectors'
                do i=1,min(writemax,ndim)
                    write(*, '(f16.10)', advance='no') (SO(i,j),j=1,min(writemax,ndim))
                    write(*,*) ""
                end do   
            end if 

            !FIND PROJECTED PRINCIPAL DIRECTIONS
            allocate(principals(numcoords,ndim))
            do i=1, ndim
                principals(:,i)=0.0D0
                do j=1, ndim
                    principals(:,i)=principals(:,i)+SO(j,i)*sigsubspan(:,j)
                end do
            end do 

            if (debug_run .and. SQNM_DEBUGLEVEL>=4) then 
                write(*, '(A)') 'principals'
                do i=1,min(writemax,numcoords)
                    write(*, '(f16.10)', advance='no') (principals(i,j),j=1,min(writemax,ndim))
                    write(*,*) ""
                end do  
            end if

            !ADJUST CURVATURES BY FINDING RESIDUES
            allocate(residues(ndim))
            !calculate residues
            do i=1, ndim
                residues(i)=0.0D0
                residuetemp(:)=0.0D0
                do j =1, ndim
                    residuetemp(:)=residuetemp(:)+SO(j,i)*siggradspan(:,j)
                end do 
                residuetemp(:)=residuetemp(:)-eigenvalues(i)*principals(:,i)
                residues(i)=norm(numcoords,residuetemp)
            end do 
            !adjust egienvalues
            do i=1,ndim
                eigenvalues(i)=sqrt(eigenvalues(i)**2+residues(i)**2)
            end do

            if (debug_run .and. SQNM_DEBUGLEVEL>=4) then 
                write(*, '(2A)') 'residues     ', ' Adjusted eigenvalues'
                do i=1,min(writemax,ndim)
                    write(*, '(i3,a,f12.5,A,g15.7)') i,': ', residues(i), ' ', eigenvalues(i)
                end do 
            end if


            !FORM SUBSPACE PORTION OF PRECONDITIONED GRADIENT
            pregrad(:)=0.0D0
            do i=1, ndim
                if (eigenvalues(i)<1.0d-50) return
                tempval=dot(numcoords,grad,principals(:,i))/eigenvalues(i)
                pregrad(:)=pregrad(:)+tempval*principals(:,i)
            end do
            if (debug_run .and. SQNM_DEBUGLEVEL>=3) then 
                write(*, '(A)') 'Preconditioned Gradient in Subspace'
                do i=1,min(writemax,numcoords)
                    write(*, '(f16.10)', advance='no') pregrad(i)
                    write(*,*) ""
                end do
            end if
            !FIND COS(ANGLE) BETWEEN GRADIENT AND PRECONDITIONED GRADIENT TO ADJUST ALPHA LATER.
            !ALSO FIND ORTHOGONAL COMPLEMENT OF GRADIENT: gradperp=(I-P')grad; invproj=(I-P')
            angle=dot(numcoords,pregrad,grad)/(norm(numcoords,pregrad)*norm(numcoords,grad))
            !Note P'=principals*principals^T; can use dot products for each entry of invproj; add 1 to diagonals to get (I-P')
            do i=1, numcoords
                invproj(i,i)=1.0d0-dot(ndim,principals(i,:),principals(i,:))
                do j=i+1, numcoords
                    invproj(j,i)=-dot(ndim,principals(i,:),principals(j,:))
                    invproj(i,j)=invproj(j,i)
                end do
            end do

            !Find rejection of the gradient projected onto the significant subspace 
            do i=1, numcoords
                gradperp(i)=0.0d0
                do j=1,numcoords
                gradperp(i)=gradperp(i)+invproj(j,i)*grad(j)
                end do
            end do
            pregrad(:)=pregrad(:)+alpha*gradperp(:) !full preconditioned gradient

            if (debug_run .and. SQNM_DEBUGLEVEL>=2) then 
                print *, 'alpha: ', alpha
                print *, 'Preconditioned Gradient, ', 'grad, ', 'gradperp'
                do i=1,min(writemax,numcoords)
                    print *, pregrad(i), grad(i), gradperp(i)
                end do
            end if

            !Iteration Cleanup
            deallocate(xdisp)
            deallocate(graddisp)
            deallocate(overlap)
            deallocate(subspan)
            deallocate(gradspan)
            deallocate(norms)
            deallocate(eigenvalues)
            deallocate(residues)
            deallocate(siggradspan)
            deallocate(sigsubspan)
            deallocate(SO)
            deallocate(principals)
        end if preconditioning
        !**********DONE FINDING PRECONDITIONED GRADIENT**********

        if (nhist==0) then !step was a pure steepest step. Do a line search.
            newxcoords(:)=xcoords(:)
            newgrad(:)=grad(:)
            enew=energy
            !try making alpha bigger
            alpha=max(1.0d-3,alpha)
            eeval=0
            call linesearch(numcoords,newxcoords,enew,newgrad,-grad,alpha,ethresh,maxeeval,infoc,eeval)
            totalfunceval=totalfunceval+eeval !update number of energy function calls
        else
            newxcoords(:)=xcoords(:)-pregrad(:)
            call potential(newxcoords,newgrad,enew,numcoords,.FALSE.)
            totalfunceval=totalfunceval+1 !update number of energy function calls
        end if

        !**********CHECK FOR STEP ACCEPTANCE**********
        accept_check: if(enew>energy+ethresh) then
            acceptstep=.false. !step rejected
            if (bio) then 
                grad(:)=oldgrad(:)
                xcoords(:)=oldxcoords(:)
            end if 
            nhist=0
            if (alpha>steepstepsize/1.0d1) alpha=max(alpha/1.5d0,minalpha)

            if (debug_run .and. SQNM_DEBUGLEVEL>=1) then 
                write(*,'(a,g12.5,a,f15.7,a,f12.6,a,f15.12,a,f15.7)') 'Step REJECTED. New angle:', angle, &
                    ', rms of gradient: ', RMS, ', Energy: ', energy, ', alpha: ', alpha, ', newE: ', enew
            end if
        else accept_check
            ediff=enew-energy
            if (debug_run .and. SQNM_DEBUGLEVEL>=1) then 
                write(*,'(a,g12.5,a,f15.7,a,f12.6,a,f15.12,a,f14.10)') 'Step ACCEPTED. New angle:', angle, &
                    ', rms of gradient: ', RMS, ',New Energy: ', enew, ', alpha: ', alpha, ', E decrease: ', ediff
            end if
            acceptstep=.true. !step accepted!
            !adjust histories, and current coordinates, gradient, and energy.
            energy=enew
            !first cycle old coordinates by shifting columns. Note that if there is an xcoord that is 
            !nhistmax old, it is deleted. Same for gradients
            do i=nhistmax-1,1,-1
                xhist(:,i+1)=xhist(:,i)
                gradhist(:,i+1)=gradhist(:,i)
            end do
            !now insert last previous x-coord and grads:
            xhist(:,1)=xcoords(:)
            gradhist(:,1)=grad(:)
            
            grad(:)=newgrad(:)
            xcoords(:)=newxcoords(:)
            !increment history length if allowed 
            if (nhist<nhistmax) nhist=nhist+1

            !algorithm seems to occasionally get stuck with small steps when nhist is high; this resets it
            if ((nhist==nhistmax .and. 1.0d3*ediff/RMS < maxrmsgrad) .or. maxrmsgrad<0.9d-4) nhist=nhist/2

            !Adjust alpha
            !alpha = alpha * alphascale
            if (angle>anglecutoff) then
                alpha=alpha*alphascale
            else
                alpha=max(alpha/alphascale, minalpha)
            end if
        end if accept_check
        !**********DONE CHECKING FOR STEP ACCEPTANCE**********

        !check for convergence. 
        gradrms=norm(numcoords,grad)
        if (RMS < maxrmsgrad) then 
            !converged! Yay! 
        if(debug_run .and. SQNM_DEBUGLEVEL>=0) then
            write(*, '(a, i4,a,i4)') 'Run ', runs, ' has converged! num iterations: ', iter
            write(*,'(a,f12.8,a,f12.6)') 'Final rms of gradient = ', gradrms, ', Final energy: ', energy
            print *, 'Total potential energy function evaluations: ', totalfunceval
        end if
            converget=.true.
            return
        end if
        !increment iterations. If not converged, loop on. 
        iter=iter+1
    end do mainloop


    !algorithm has not converged.
   if(debug_run .and. SQNM_DEBUGLEVEL>=0) then
        write(*, '(a, i4,a,i4)') 'Run ', runs, ' failed to converge.'
        write(*,'(a,f12.8,a,f12.6)') 'Final rms of gradient = ', gradrms, ', Final energy: ', energy
        print *, 'Total potential energy function evaluations: ', totalfunceval
   end if
end subroutine sqnm

!================================================================================================================================
!==================================================================================================================================
!     ******************************
!     LINE SEARCH ROUTINE linesearch
!     ******************************

SUBROUTINE linesearch(N,X,F,G,S,STP,FTOL,MAXFEV,INFO,NFEV)
    use sqnm_func, only: dot
    implicit none
    integer, intent(in)             :: N,MAXFEV
    integer, intent(out)            :: INFO
    integer, intent(inout)          :: NFEV
    double precision,intent(inout)  :: F,STP
    double precision, intent(in)    :: FTOL
    double precision, intent(inout) :: X(N),G(N)
    double precision, intent(in)    :: S(N)

!          SUBROUTINE linesearch 
! Taken from Nocedal's L-BFGS algorith (MCSRCH) at http://users.iems.northwestern.edu/~nocedal/lbfgs.html.
! J. Nocedal. Updating Quasi-Newton Matrices with Limited Storage (1980), Mathematics of Computation 35, pp. 773-782.
! Note in Nocedal's implemntation, it is called MCSRCH, not linesearch. This does not include any of the actual LBFGS algorithm
! 
! Modifications were made to bring up to Fortran 90 standards.  GOTO/CONTINUE commands were removed, and no returns to the 
! parent function are present since the energies and gradients are obtained from potential.f90, and several do loops were 
! removed in favour of operations on arrays. No interaction with the calling function occurs. Typically, things in lowercase
! were added by James Klimavicz March 2017
!
! THE PURPOSE OF linesearch IS TO FIND A STEP WHICH SATISFIES A SUFFICIENT DECREASE CONDITION AND A CURVATURE CONDITION.
! AT EACH STAGE THE SUBROUTINE UPDATES AN INTERVAL OF UNCERTAINTY WITH ENDPOINTS STX AND STY. THE INTERVAL OF UNCERTAINTY 
! IS INITIALLY CHOSEN SO THAT IT CONTAINS A MINIMIZER OF THE MODIFIED FUNCTION
!     F(X+STP*S) - F(X) - FTOL*STP*(GRADF(X)'S).
! IF A STEP IS OBTAINED FOR WHICH THE MODIFIED FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE, THEN THE 
! INTERVAL OF UNCERTAINTY IS CHOSEN SO THAT IT CONTAINS A MINIMIZER OF F(X+STP*S).
! THE ALGORITHM IS DESIGNED TO FIND A STEP WHICH SATISFIES THE SUFFICIENT DECREASE CONDITION
!     F(X+STP*S) .LE. F(X) + FTOL*STP*(GRADF(X)'S),
! AND THE CURVATURE CONDITION
!     ABS(GRADF(X+STP*S)'S)) .LE. GTOL*ABS(GRADF(X)'S).
! IF FTOL IS LESS THAN GTOL AND IF, FOR EXAMPLE, THE FUNCTION IS BOUNDED BELOW, THEN THERE IS ALWAYS A STEP WHICH SATISFIES BOTH
! CONDITIONS. IF NO STEP CAN BE FOUND WHICH SATISFIES BOTH CONDITIONS, THEN THE ALGORITHM USUALLY STOPS WHEN ROUNDING ERRORS PREVENT
! FURTHER PROGRESS. IN THIS CASE STP ONLY SATISFIES THE SUFFICIENT DECREASE CONDITION.
! 
! THE SUBROUTINE STATEMENT IS SUBROUTINE 
!       MCSRCH(N,X,F,G,S,STP,FTOL,XTOL, MAXFEV,INFO,NFEV,WA) WHERE
! N IS A POSITIVE INTEGER INPUT VARIABLE SET TO THE NUMBER OF VARIABLES.
! X IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE BASE POINT FOR THE LINE SEARCH. ON OUTPUT IT CONTAINS X + STP*S.
! F IS A VARIABLE. ON INPUT IT MUST CONTAIN THE VALUE OF F AT X. ON OUTPUT IT CONTAINS THE VALUE OF F AT X + STP*S.
! G IS AN ARRAY OF LENGTH N. ON INPUT IT MUST CONTAIN THE GRADIENT OF F AT X. ON OUTPUT IT CONTAINS THE GRADIENT OF F AT X + STP*S.
! S IS AN INPUT ARRAY OF LENGTH N WHICH SPECIFIES THE SEARCH DIRECTION.
! STP IS A NONNEGATIVE VARIABLE. ON INPUT STP CONTAINS AN INITIAL ESTIMATE OF A SATISFACTORY STEP. ON OUTPUT STP CONTAINS 
!     THE FINAL ESTIMATE.
! FTOL AND GTOL ARE NONNEGATIVE INPUT VARIABLES. TERMINATION OCCURS WHEN THE SUFFICIENT DECREASE CONDITION AND THE DIRECTIONAL 
!     DERIVATIVE CONDITION ARE SATISFIED.
! STPMIN AND STPMAX ARE NONNEGATIVE INPUT VARIABLES WHICH SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP.
! MAXFEV IS A POSITIVE INTEGER INPUT VARIABLE. TERMINATION OCCURS WHEN THE NUMBER OF CALLS TO FCN IS AT LEAST MAXFEV BY THE END OF 
! AN ITERATION.
! INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS:
!   INFO = 0  IMPROPER INPUT PARAMETERS.
!   INFO = 1  THE SUFFICIENT DECREASE CONDITION AND THE DIRECTIONAL DERIVATIVE CONDITION HOLD.
!   INFO = 2  RELATIVE WIDTH OF THE INTERVAL OF UNCERTAINTY IS AT MOST XTOL.
!   INFO = 3  NUMBER OF CALLS TO FCN HAS REACHED MAXFEV.
!   INFO = 4  THE STEP IS AT THE LOWER BOUND STPMIN.
!   INFO = 5  THE STEP IS AT THE UPPER BOUND STPMAX.
!   INFO = 6  ROUNDING ERRORS PREVENT FURTHER PROGRESS. THERE MAY NOT BE A STEP WHICH SATISFIES THE
!          SUFFICIENT DECREASE AND CURVATURE CONDITIONS. TOLERANCES MAY BE TOO SMALL.
!   info = 7  search direction was not a descent direction
! NFEV IS AN INTEGER OUTPUT VARIABLE SET TO THE NUMBER OF CALLS TO FCN.
! WA IS A WORK ARRAY OF LENGTH N.
! SUBPROGRAMS CALLED: MYMCSTEP
! ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983 JORGE J. MORE', DAVID J. THUENTE

    double precision, parameter :: GTOL = 1.0d-1 !controls the accuracy of the line search routine MCSRCH.
    double precision, parameter :: STPMIN = 1.0d-20
    double precision, parameter :: STPMAX = 1.0d5
    double precision, parameter :: XTOL = epsilon(1.0d0) 

    INTEGER :: INFOC=6
    LOGICAL :: BRACKT,STAGE1
    DOUBLE PRECISION :: DG,DGM,DGINIT,DGTEST,DGX,DGXM,DGY,DGYM,FINIT,FTEST1,FM,FX,FXM,FY,FYM,STX,STY
    DOUBLE PRECISION :: STMIN,STMAX,WIDTH,WIDTH1, WA(N)
    double precision, parameter :: P5 = 0.5d0
    double precision, parameter :: P66 = 2.0d0/3.0d0
    double precision, parameter :: XTRAPF=4.0d0
    double precision, parameter :: zero=0.0d0

    INFO=0

    ! CHECK THE INPUT PARAMETERS FOR ERRORS.

    IF (N<= 0.OR.STP<=ZERO.OR.FTOL<ZERO.OR.GTOL<ZERO.OR.XTOL<ZERO.OR.STPMIN<ZERO.OR.STPMAX<STPMIN.OR.MAXFEV<0) RETURN

    ! COMPUTE THE INITIAL GRADIENT IN THE SEARCH DIRECTION AND CHECK THAT S IS A DESCENT DIRECTION.
    DGINIT=dot(N,G,S)
        IF (DGINIT .GE. ZERO) then
            info = 8        
            RETURN
        ENDIF

    !INITIALIZE LOCAL VARIABLES.
    BRACKT = .FALSE.
    STAGE1 = .TRUE.
    NFEV = 0
    FINIT = F
    DGTEST = FTOL*DGINIT
    WIDTH = STPMAX - STPMIN
    WIDTH1 = WIDTH/P5
    WA(:) = X(:)

    ! THE VARIABLES STX, FX, DGX CONTAIN THE VALUES OF THE STEP, FUNCTION, AND DIRECTIONAL DERIVATIVE AT THE BEST STEP.
    ! THE VARIABLES STY, FY, DGY CONTAIN THE VALUE OF THE STEP, FUNCTION, AND DERIVATIVE AT THE OTHER ENDPOINT OF THE INTERVAL 
    ! OF UNCERTAINTY. THE VARIABLES STP, F, DG CONTAIN THE VALUES OF THE STEP, FUNCTION, AND DERIVATIVE AT THE CURRENT STEP.
    STX = ZERO
    FX = FINIT
    DGX = DGINIT
    STY = ZERO
    FY = FINIT
    DGY = DGINIT

    ! START OF ITERATION.
    do while (info==0)
        ! SET THE MINIMUM AND MAXIMUM STEPS TO CORRESPOND TO THE PRESENT INTERVAL OF UNCERTAINTY.
        IF (BRACKT) THEN
            STMIN = MIN(STX,STY)
            STMAX = MAX(STX,STY)
        ELSE
            STMIN = STX
            STMAX = STP + XTRAPF*(STP - STX)
        END IF

        ! FORCE THE STEP TO BE WITHIN THE BOUNDS STPMAX AND STPMIN.
         STP = MAX(STP,STPMIN)
         STP = MIN(STP,STPMAX)

        ! IF AN UNUSUAL TERMINATION IS TO OCCUR THEN LET STP BE THE LOWEST POINT OBTAINED SO FAR.
        IF ((BRACKT .AND. (STP <= STMIN .OR. STP >= STMAX)) .OR. NFEV >= MAXFEV-1 .OR. INFOC == 0 .OR. &
            (BRACKT .AND. STMAX-STMIN <= XTOL*STMAX)) STP = STX

        ! EVALUATE THE FUNCTION AND GRADIENT AT STP AND COMPUTE THE DIRECTIONAL DERIVATIVE.
        X(:) = WA(:)+STP*S(:)
        call potential(X,G,F,N,.FALSE.)
        NFEV = NFEV + 1
        DG = dot(N,G,S)

        FTEST1 = FINIT + STP*DGTEST

        ! TEST FOR CONVERGENCE.
        IF ((BRACKT .AND. (STP .LE. STMIN .OR. STP .GE. STMAX)) .OR. INFOC .EQ. 0) INFO = 6
        IF (STP .EQ. STPMAX .AND. F .LE. FTEST1 .AND. DG .LE. DGTEST) INFO = 5
        IF (STP .EQ. STPMIN .AND. (F .GT. FTEST1 .OR. DG .GE. DGTEST)) INFO = 4
        IF (NFEV .GE. MAXFEV) INFO = 3
        IF (BRACKT .AND. STMAX-STMIN .LE. XTOL*STMAX) INFO = 2
        IF (F .LE. FTEST1 .AND. ABS(DG) .LE. GTOL*(-DGINIT)) INFO = 1

        ! CHECK FOR TERMINATION.
        IF (INFO .NE. 0) RETURN

        ! IN THE FIRST STAGE WE SEEK A STEP FOR WHICH THE MODIFIED FUNCTION HAS A NONPOSITIVE VALUE AND NONNEGATIVE DERIVATIVE.
        IF (STAGE1 .AND. F .LE. FTEST1 .AND. DG .GE. MIN(FTOL,GTOL)*DGINIT) STAGE1 = .FALSE.

        ! A MODIFIED FUNCTION IS USED TO PREDICT THE STEP ONLY IF WE HAVE NOT OBTAINED A STEP FOR WHICH THE MODIFIED
        ! FUNCTION HAS A NONPOSITIVE FUNCTION VALUE AND NONNEGATIVE DERIVATIVE, AND IF A LOWER FUNCTION VALUE HAS BEEN
        ! OBTAINED BUT THE DECREASE IS NOT SUFFICIENT.
        IF (STAGE1 .AND. F .LE. FX .AND. F .GT. FTEST1) THEN
            ! DEFINE THE MODIFIED FUNCTION AND DERIVATIVE VALUES.
            FM = F - STP*DGTEST
            FXM = FX - STX*DGTEST
            FYM = FY - STY*DGTEST
            DGM = DG - DGTEST
            DGXM = DGX - DGTEST
            DGYM = DGY - DGTEST
            ! CALL MYMCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY AND TO COMPUTE THE NEW STEP.
            CALL MYMCSTEP(STX,FXM,DGXM,STY,FYM,DGYM,STP,FM,DGM,BRACKT,STMIN,STMAX,INFOC)

            ! RESET THE FUNCTION AND GRADIENT VALUES FOR F.
            FX = FXM + STX*DGTEST
            FY = FYM + STY*DGTEST
            DGX = DGXM + DGTEST
            DGY = DGYM + DGTEST
        ELSE
            ! CALL MYMCSTEP TO UPDATE THE INTERVAL OF UNCERTAINTY AND TO COMPUTE THE NEW STEP.
            CALL MYMCSTEP(STX,FX,DGX,STY,FY,DGY,STP,F,DG,BRACKT,STMIN,STMAX,INFOC)
        END IF

        ! FORCE A SUFFICIENT DECREASE IN THE SIZE OF THE INTERVAL OF UNCERTAINTY.
        IF (BRACKT) THEN
            IF (ABS(STY-STX) .GE. P66*WIDTH1) STP = STX + P5*(STY - STX)
            WIDTH1 = WIDTH
            WIDTH = ABS(STY-STX)
        END IF
        
    end do
end subroutine linesearch

!===========================================================================================================
SUBROUTINE MYMCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,STPMIN,STPMAX,INFO)
implicit none
    INTEGER INFO
    DOUBLE PRECISION STX,FX,DX,STY,FY,DY,STP,FP,DP,STPMIN,STPMAX
    LOGICAL BRACKT,BOUND

! THE PURPOSE OF MYMCSTEP IS TO COMPUTE A SAFEGUARDED STEP FOR A LINESEARCH AND TO UPDATE AN INTERVAL OF UNCERTAINTY FOR
! A MINIMIZER OF THE FUNCTION.
! THE PARAMETER STX CONTAINS THE STEP WITH THE LEAST FUNCTION VALUE. THE PARAMETER STP CONTAINS THE CURRENT STEP. IT IS
! ASSUMED THAT THE DERIVATIVE AT STX IS NEGATIVE IN THE DIRECTION OF THE STEP. IF BRACKT IS SET TRUE THEN A
! MINIMIZER HAS BEEN BRACKETED IN AN INTERVAL OF UNCERTAINTY WITH ENDPOINTS STX AND STY.
! THE SUBROUTINE STATEMENT IS SUBROUTINE MYMCSTEP(STX,FX,DX,STY,FY,DY,STP,FP,DP,BRACKT,STPMIN,STPMAX,INFO) WHERE
! STX, FX, AND DX ARE VARIABLES WHICH SPECIFY THE STEP, THE FUNCTION, AND THE DERIVATIVE AT THE BEST STEP OBTAINED
!     SO FAR. THE DERIVATIVE MUST BE NEGATIVE IN THE DIRECTION OF THE STEP, THAT IS, DX AND STP-STX MUST HAVE OPPOSITE
!     SIGNS. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY.
! STY, FY, AND DY ARE VARIABLES WHICH SPECIFY THE STEP, THE FUNCTION, AND THE DERIVATIVE AT THE OTHER ENDPOINT OF
!     THE INTERVAL OF UNCERTAINTY. ON OUTPUT THESE PARAMETERS ARE UPDATED APPROPRIATELY.
! STP, FP, AND DP ARE VARIABLES WHICH SPECIFY THE STEP, THE FUNCTION, AND THE DERIVATIVE AT THE CURRENT STEP.
!     IF BRACKT IS SET TRUE THEN ON INPUT STP MUST BE BETWEEN STX AND STY. ON OUTPUT STP IS SET TO THE NEW STEP.
!     BRACKT IS A LOGICAL VARIABLE WHICH SPECIFIES IF A MINIMIZER HAS BEEN BRACKETED. IF THE MINIMIZER HAS NOT BEEN BRACKETED
!     THEN ON INPUT BRACKT MUST BE SET FALSE. IF THE MINIMIZER IS BRACKETED THEN ON OUTPUT BRACKT IS SET TRUE.
! STPMIN AND STPMAX ARE INPUT VARIABLES WHICH SPECIFY LOWER AND UPPER BOUNDS FOR THE STEP.
! INFO IS AN INTEGER OUTPUT VARIABLE SET AS FOLLOWS: IF INFO = 1,2,3,4,5, THEN THE STEP HAS BEEN COMPUTED
! ACCORDING TO ONE OF THE FIVE CASES BELOW. OTHERWISE INFO = 0, AND THIS INDICATES IMPROPER INPUT PARAMETERS.
! ARGONNE NATIONAL LABORATORY. MINPACK PROJECT. JUNE 1983 JORGE J. MORE', DAVID J. THUENTE

    DOUBLE PRECISION GAMMA,P,Q,R,S,SGND,STPC,STPF,STPQ,THETA
    INFO = 0

    ! CHECK THE INPUT PARAMETERS FOR ERRORS.
    IF ((BRACKT .AND. (STP <= MIN(STX,STY) .OR. STP >= MAX(STX,STY))) .OR. DX*(STP-STX) >= 0.0 .OR. STPMAX < STPMIN) RETURN
    
    ! DETERMINE IF THE DERIVATIVES HAVE OPPOSITE SIGN.
    SGND = DP*(DX/ABS(DX))

    ! FIRST CASE. A HIGHER FUNCTION VALUE. THE MINIMUM IS BRACKETED. IF THE CUBIC STEP IS CLOSER TO STX THAN THE QUADRATIC STEP, 
    ! THE CUBIC STEP IS TAKEN, ELSE THE AVERAGE OF THE CUBIC AND QUADRATIC STEPS IS TAKEN.
    cases: IF (FP .GT. FX) THEN
        if ((STP - STX)==0.0d0) then 
            info=0
            return
        endif
        INFO = 1
        BOUND = .TRUE.
        THETA = 3*(FX - FP)/(STP - STX) + DX + DP
        S = MAX(ABS(THETA),ABS(DX),ABS(DP))
        if (s==0.0d0) then 
            info=0
            return
        endif
        GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
        IF (STP .LT. STX) GAMMA = -GAMMA
        P = (GAMMA - DX) + THETA
        Q = ((GAMMA - DX) + GAMMA) + DP
        if (q==0.0d0) then 
            info=0
            return
        endif
        R = P/Q
        STPC = STX + R*(STP - STX)
        STPQ = STX + ((DX/((FX-FP)/(STP-STX)+DX))/2)*(STP - STX)
        IF (ABS(STPC-STX) .LT. ABS(STPQ-STX)) THEN
            STPF = STPC
        ELSE
            STPF = STPC + (STPQ - STPC)/2
        END IF
        BRACKT = .TRUE.

    ! SECOND CASE. A LOWER FUNCTION VALUE AND DERIVATIVES OF OPPOSITE SIGN. THE MINIMUM IS BRACKETED. IF THE CUBIC
    ! STEP IS CLOSER TO STX THAN THE QUADRATIC (SECANT) STEP, THE CUBIC STEP IS TAKEN, ELSE THE QUADRATIC STEP IS TAKEN.
    ELSE IF (SGND .LT. 0.0) THEN cases
        INFO = 2
        if ((STP - STX)==0.0d0) then 
            info=0
            return
        endif
        BOUND = .FALSE.
        THETA = 3*(FX - FP)/(STP - STX) + DX + DP
        S = MAX(ABS(THETA),ABS(DX),ABS(DP))
        if (s==0.0d0) then 
            info=0
            return
        endif
        GAMMA = S*SQRT((THETA/S)**2 - (DX/S)*(DP/S))
        IF (STP .GT. STX) GAMMA = -GAMMA
        P = (GAMMA - DP) + THETA
        Q = ((GAMMA - DP) + GAMMA) + DX
        if (q==0.0d0) then 
            info=0
            return
        endif
        R = P/Q
        if ((DP-DX)==0.0d0) then 
            info=0
            return
        endif
        STPC = STP + R*(STX - STP)
        STPQ = STP + (DP/(DP-DX))*(STX - STP)
        IF (ABS(STPC-STP) .GT. ABS(STPQ-STP)) THEN
            STPF = STPC
        ELSE
            STPF = STPQ
        END IF
        BRACKT = .TRUE.

    ! THIRD CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DECREASES. THE CUBIC
    ! STEP IS ONLY USED IF THE CUBIC TENDS TO INFINITY IN THE DIRECTION OF THE STEP OR IF THE MINIMUM OF THE CUBIC IS BEYOND STP. 
    ! OTHERWISE THE CUBIC STEP IS DEFINED TO BE EITHER STPMIN OR STPMAX. THE QUADRATIC (SECANT) STEP IS ALSO COMPUTED AND IF THE 
    ! MINIMUM IS BRACKETED THEN THE THE STEP CLOSEST TO STX IS TAKEN, ELSE THE STEP FARTHEST AWAY IS TAKEN.
    ELSE IF (ABS(DP) .LT. ABS(DX)) THEN cases
        INFO = 3
        BOUND = .TRUE.
        if ((STP - STX)==0.0d0) then 
            info=0
            return
        endif
        THETA = 3*(FX - FP)/(STP - STX) + DX + DP
        S = MAX(ABS(THETA),ABS(DX),ABS(DP))
        if (s==0.0d0) then 
            info=0
            return
        endif
        !THE CASE GAMMA = 0 ONLY ARISES IF THE CUBIC DOES NOT TEND TO INFINITY IN THE DIRECTION OF THE STEP.
        GAMMA = S*SQRT(MAX(0.0D0,(THETA/S)**2 - (DX/S)*(DP/S)))
        IF (STP .GT. STX) GAMMA = -GAMMA
        P = (GAMMA - DP) + THETA
        Q = (GAMMA + (DX - DP)) + GAMMA
        if (q==0.0d0) then 
            info=0
            return
        endif
        R = P/Q
        IF (R .LT. 0.0d0 .AND. GAMMA .NE. 0.0d0) THEN
            STPC = STP + R*(STX - STP)
        ELSE IF (STP .GT. STX) THEN
            STPC = STPMAX
        ELSE
            STPC = STPMIN
        END IF
        if ((DP-DX)==0.0d0) then 
            info=0
            return
        endif
        STPQ = STP + (DP/(DP-DX))*(STX - STP)
        IF (BRACKT) THEN
            IF (ABS(STP-STPC) .LT. ABS(STP-STPQ)) THEN
                STPF = STPC
            ELSE
                STPF = STPQ
            END IF
        ELSE
            IF (ABS(STP-STPC) .GT. ABS(STP-STPQ)) THEN
                STPF = STPC
            ELSE
                STPF = STPQ
            END IF
        END IF

    ! FOURTH CASE. A LOWER FUNCTION VALUE, DERIVATIVES OF THE SAME SIGN, AND THE MAGNITUDE OF THE DERIVATIVE DOES
    ! NOT DECREASE. IF THE MINIMUM IS NOT BRACKETED, THE STEP IS EITHER STPMIN OR STPMAX, ELSE THE CUBIC STEP IS TAKEN.
    ELSE cases
        INFO = 4
        BOUND = .FALSE.
        IF (BRACKT) THEN
            if ((STY - STP)==0.0d0) then 
                info=0
                return
            endif
            THETA = 3*(FP - FY)/(STY - STP) + DY + DP
            S = MAX(ABS(THETA),ABS(DY),ABS(DP))
            if (s==0.0d0) then 
                info=0
                return
            endif
            GAMMA = S*SQRT((THETA/S)**2 - (DY/S)*(DP/S))
            IF (STP .GT. STY) GAMMA = -GAMMA
            P = (GAMMA - DP) + THETA
            Q = ((GAMMA - DP) + GAMMA) + DY
            if (q==0.0d0) then 
                info=0
                return
            endif
            R = P/Q
            STPC = STP + R*(STY - STP)
            STPF = STPC
        ELSE IF (STP .GT. STX) THEN
            STPF = STPMAX
        ELSE
            STPF = STPMIN
        END IF
    END IF cases

    ! UPDATE THE INTERVAL OF UNCERTAINTY. THIS UPDATE DOES NOT DEPEND ON THE NEW STEP OR THE CASE ANALYSIS ABOVE.
    IF (FP .GT. FX) THEN
        STY = STP
        FY = FP
        DY = DP
    ELSE
        IF (SGND .LT. 0.0d0) THEN
            STY = STX
            FY = FX
            DY = DX
            END IF
        STX = STP
        FX = FP
        DX = DP
    END IF

    ! COMPUTE THE NEW STEP AND SAFEGUARD IT.
    STPF = MIN(STPMAX,STPF)
    STPF = MAX(STPMIN,STPF)
    STP = STPF
    IF (BRACKT .AND. BOUND) THEN
        IF (STY .GT. STX) THEN
            STP = MIN(STX+0.66*(STY-STX),STP)
        ELSE
            STP = MAX(STX+0.66*(STY-STX),STP)
        END IF
    END IF
    RETURN
END subroutine mymcstep 
