MODULE MOVES
! This module is for the magical new way of doing custom movesets and also fixing our
! horrible old way of doing moves with 15 different implementations of Cartesian steps.
! Please code nicely in here...GOTOs will be removed, as will unclear variable names.
! Consistent indentation is mandatory!

CONTAINS

SUBROUTINE CARTESIAN_SPHERE(XYZ, MAX_STEP, ATOM_LIST)
! Add a random spherically symmetric displacement of up to MAXSTEP to each atom
! in the ATOM_LIST array if present, or all atoms if not.
!
! Arguments
! ---------
!
! Required: 
! XYZ(in/out): coordinates array from GMIN, in Cartesian coordinates
! MAX_STEP(in): the maximum step size
!
! Optional:
! ATOM_LIST(in): list of atoms to be moved - if omitted, all are moved
 
   USE COMMONS, ONLY: AMBERMUTATIONT !need to reallocate atom_mask size!

! The VEC3 module (vec3.f90) contains helper functions for handling vectors and matricies
   USE VEC3
! The SANITY module contains sanity check functions
   USE SANITY
   IMPLICIT NONE
   INTEGER                                       :: I 
   INTEGER                                       :: NUM_ATOMS
   INTEGER, OPTIONAL, DIMENSION(:), INTENT(IN)   :: ATOM_LIST
   DOUBLE PRECISION                              :: DPRAND
   DOUBLE PRECISION, INTENT(IN)                  :: MAX_STEP
   DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: XYZ
   LOGICAL, ALLOCATABLE, DIMENSION(:)            :: ATOM_MASK
   LOGICAL                                       :: TEST

! Sanity check - are the coordinates in XYZ Cartesian? 
! Check if the SIZE is a multiple of 3
   TEST=.FALSE.
   TEST=CHECK_DIMENSION(SIZE(XYZ),3)
   IF (.NOT.TEST) THEN
      STOP 'Coordinates in a non-Cartesian basis passed to CARTESIAN_SPHERE'
   ENDIF

! Set NUM_ATOMS
   NUM_ATOMS = SIZE(XYZ) / 3

! Set up ATOM_MASK
   IF (.NOT. ALLOCATED(ATOM_MASK)) ALLOCATE(ATOM_MASK(NUM_ATOMS))
   IF ((.NOT.(SIZE(ATOM_MASK).EQ.NUM_ATOMS)).AND.AMBERMUTATIONT) THEN
      DEALLOCATE(ATOM_MASK)
      ALLOCATE(ATOM_MASK(NUM_ATOMS))
   ENDIF
   ATOM_MASK = .FALSE.

! Check to see if an ATOM_LIST was provided
   IF (PRESENT(ATOM_LIST)) THEN
! If so, determine which atoms the move applies to and set up ATOM_MASK
      DO I = 1, SIZE(ATOM_LIST)
         ATOM_MASK(ATOM_LIST(I)) = .TRUE.
      END DO
   ELSE
! Otherwise, apply the move to all atoms
      ATOM_MASK = .TRUE.
   ENDIF

! Apply the move to the atoms specified 
   DO I = 1, NUM_ATOMS
! Skip atoms we do not want to move
      IF (.NOT. ATOM_MASK(I)) CYCLE
! Otherwise apply the move
      XYZ(3*I-2:3*I)=XYZ(3*I-2:3*I)+VEC_RANDOM()*(DPRAND()**(1.0D0/3.0D0))*MAX_STEP
   ENDDO

END SUBROUTINE CARTESIAN_SPHERE

SUBROUTINE CARTESIAN_SIMPLE(XYZ, MAX_STEP, ATOM_LIST)
! Add a random displacement of up to MAXSTEP to each atom
! in the ATOM_LIST array if present, or all atoms if not.
!
! Arguments
! ---------
!
! Required: 
! XYZ(in/out): coordinates array from GMIN, in Cartesian coordinates
! MAX_STEP(in): the maximum step size
!
! Optional:
! ATOM_LIST(in): list of atoms to be moved - if omitted, all are moved

! The SANITY module contains sanity check functions
   USE SANITY 
   IMPLICIT NONE
   INTEGER                                       :: I 
   INTEGER                                       :: NUM_ATOMS
   INTEGER, OPTIONAL, DIMENSION(:), INTENT(IN)   :: ATOM_LIST
   DOUBLE PRECISION                              :: DPRAND
   DOUBLE PRECISION                              :: RANDOMX
   DOUBLE PRECISION                              :: RANDOMY
   DOUBLE PRECISION                              :: RANDOMZ
   DOUBLE PRECISION, INTENT(IN)                  :: MAX_STEP
   DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT) :: XYZ
   LOGICAL, ALLOCATABLE, DIMENSION(:)            :: ATOM_MASK
   LOGICAL                                       :: TEST

! Sanity check - are the coordinates in XYZ Cartesian? 
! Check if the SIZE is a multiple of 3
   TEST=.FALSE.
   TEST=CHECK_DIMENSION(SIZE(XYZ),3)
   IF (.NOT.TEST) THEN
      STOP 'Coordinates in a non-Cartesian basis passed to CARTESIAN_SIMPLE'
   ENDIF

! Set NUM_ATOMS
   NUM_ATOMS = SIZE(XYZ) / 3

! Set up ATOM_MASK
   IF (.NOT. ALLOCATED(ATOM_MASK)) ALLOCATE(ATOM_MASK(NUM_ATOMS))
   ATOM_MASK = .FALSE.

! Check to see if an ATOM_LIST was provided
   IF (PRESENT(ATOM_LIST)) THEN
! If so, determine which atoms the move applies to and set up ATOM_MASK
      DO I = 1, SIZE(ATOM_LIST)
         ATOM_MASK(ATOM_LIST(I)) = .TRUE.
      END DO
   ELSE
! Otherwise, apply the move to all atoms
      ATOM_MASK = .TRUE.
   ENDIF

! Apply the move to the atoms specified 
   DO I = 1, NUM_ATOMS
! Skip atoms we do not want to move
      IF (.NOT. ATOM_MASK(I)) CYCLE
! Otherwise apply the move
! Draw a random number between -1 and 1 for each coordinate
      RANDOMX=(DPRAND()-0.5D0)*2.0D0
      RANDOMY=(DPRAND()-0.5D0)*2.0D0
      RANDOMZ=(DPRAND()-0.5D0)*2.0D0
! Displace each coordinate
      XYZ(3*I-2)=XYZ(3*I-2)+MAX_STEP*RANDOMX
      XYZ(3*I-1)=XYZ(3*I-1)+MAX_STEP*RANDOMY
      XYZ(3*I  )=XYZ(3*I  )+MAX_STEP*RANDOMZ
   ENDDO

END SUBROUTINE CARTESIAN_SIMPLE

SUBROUTINE ROTATION_ABOUT_AXIS(XYZ, VECTOR_START_XYZ, &
                               VECTOR_END_XYZ, ANGLE_DEGS, ATOM_LIST)
!
! Rotate the coordinates of the atoms in ATOM_LIST about the line from VECTOR_START_XYZ
! to VECTOR_END_XYZ through ANGLE degrees.
!
! Arguments
! ---------
!
! Required:
! XYZ(in/out): coordinates array from GMIN, in Cartesian coordinates
! VECTOR_START_XYZ(in): start of the line about which to rotate, in Cartesian coords
! VECTOR_END_XYZ(in): end of the line about which to rotate, in Cartesian coords
! ANGLE_DEGS(in): angle through which to rotate, in degrees
!
! Optional:
! ATOM_LIST(in): list of atoms to be rotated
!
! The SANITY module contains sanity check functions
   USE SANITY 
   IMPLICIT NONE
! Arguments
   DOUBLE PRECISION, DIMENSION(:), INTENT(INOUT)   :: XYZ
   INTEGER, OPTIONAL, DIMENSION(:), INTENT(IN)     :: ATOM_LIST
   DOUBLE PRECISION, DIMENSION(3), INTENT(IN)      :: VECTOR_START_XYZ
   DOUBLE PRECISION, DIMENSION(3), INTENT(IN)      :: VECTOR_END_XYZ
   DOUBLE PRECISION, INTENT(IN)                    :: ANGLE_DEGS
! Constants
! Variables
   DOUBLE PRECISION                                :: PI
   DOUBLE PRECISION                                :: DEGS_OVER_RADS
   DOUBLE PRECISION                                :: COS_THETA, SIN_THETA
   DOUBLE PRECISION                                :: A, B, C
   DOUBLE PRECISION                                :: D, E, F
   DOUBLE PRECISION                                :: U, V, W
   DOUBLE PRECISION                                :: X, Y, Z
   DOUBLE PRECISION                                :: VECTOR_MAG
   INTEGER                                         :: NUM_ATOMS
   INTEGER                                         :: I
   LOGICAL, ALLOCATABLE, DIMENSION(:)              :: ATOM_MASK
   LOGICAL                                         :: TEST
! Function declarations
   DOUBLE PRECISION                                :: DNRM2

! Define PI using ATAN and the conversion factor between degrees and radians.
! x degrees = x * DEGS_OVER_RADS radians
   PI = 4.0D0 * ATAN(1.0D0)
   DEGS_OVER_RADS = PI / 180.0D0

! Calculate cos and sin of the angle.
   SIN_THETA = SIN(ANGLE_DEGS * DEGS_OVER_RADS)
   COS_THETA = COS(ANGLE_DEGS * DEGS_OVER_RADS)

! Assign variables according to the functional form described on:
! http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html

   A = VECTOR_START_XYZ(1)
   B = VECTOR_START_XYZ(2)
   C = VECTOR_START_XYZ(3)
   
   D = VECTOR_END_XYZ(1)
   E = VECTOR_END_XYZ(2)
   F = VECTOR_END_XYZ(3)

! DNRM2(INTEGER N, REAL X(N), INTEGER INCX)
   VECTOR_MAG = DNRM2(3, VECTOR_END_XYZ - VECTOR_START_XYZ, 1)

   U = (D - A) / VECTOR_MAG
   V = (E - B) / VECTOR_MAG
   W = (F - C) / VECTOR_MAG

! Work out how many atoms the coordinates passed describe
   NUM_ATOMS = SIZE(XYZ) / 3 
! Sanity check - are the coordinates in XYZ Cartesian? 
! Check if the SIZE is a multiple of 3
   TEST=.FALSE.
   TEST=CHECK_DIMENSION(SIZE(XYZ),3)
   IF (.NOT.TEST) THEN
      STOP 'Coordinates in a non-Cartesian basis passed to ROTATION_ABOUT_AXIS'
   ENDIF

! Convert ATOM_LIST (a list of atoms to be rotated) into ATOM_MASK, which can
! applied in a WHERE loop.
   IF (.NOT. ALLOCATED(ATOM_MASK)) ALLOCATE(ATOM_MASK(NUM_ATOMS))
   ATOM_MASK = .FALSE.

! Check to see if an ATOM_LIST was provided
   IF (PRESENT(ATOM_LIST)) THEN
! If so, determine which atoms the move applies to and set up ATOM_MASK
      DO I = 1, SIZE(ATOM_LIST)
         ATOM_MASK(ATOM_LIST(I)) = .TRUE.
      END DO
   ELSE
! Otherwise, apply the move to all atoms
      ATOM_MASK = .TRUE.
   ENDIF

! Loop through and apply the formula if the atom described is in ATOM_LIST
! N.B. I've rearranged the formula so that the variables cycle, since I think it
! makes checking a bit easier.
   DO I = 1, NUM_ATOMS
      IF (.NOT. ATOM_MASK(I)) CYCLE
      X = XYZ(3 * I - 2)
      Y = XYZ(3 * I - 1)
      Z = XYZ(3 * I    )
      XYZ(3 * I - 2) = (A * (V**2 + W**2) - U * (B*V + C*W - U*X - V*Y - W*Z)) * (1 - COS_THETA) + &
                       X * COS_THETA + &
                       (-C*V + B*W + V*Z - W*Y) * SIN_THETA
      XYZ(3 * I - 1) = (B * (W**2 + U**2) - V * (C*W + A*U - U*X - V*Y - W*Z)) * (1 - COS_THETA) + &
                       Y * COS_THETA + &
                       (-A*W + C*U + W*X - U*Z) * SIN_THETA
      XYZ(3 * I    ) = (C * (U**2 + V**2) - W * (A*U + B*V - U*X - V*Y - W*Z)) * (1 - COS_THETA) + &
                       Z * COS_THETA + &
                       (-B*U + A*V + U*Y - V*X) * SIN_THETA
   END DO

END SUBROUTINE ROTATION_ABOUT_AXIS

SUBROUTINE MACROION_MOVES(j1,y,movableatomlist1,nmovableatoms,ligmovet,blockmovet,nblocks,atomsinblock1,LOCALSTEP)
! sf344> this is a cooperative move set for the MACROION model. It could be adapted for other systems, 
! where dynamically determined groups of atoms have to be moved together.

   use modamber9 ! macroion related keywords are stored in this module
   use commons, only : natoms, change_temp, newres_temp, perct, perccut, frozen
!   use nblist, only : nbflag, skinnb
   use porfuncs

   implicit none


   integer itime1, now(3), nmovableatoms,i,j,movableatomlist1(nmovableatoms),nblocks,atomsinblock1(nblocks),offset1,k,j1
   double precision        :: y(3*natoms), grad(3*natoms), ereal,ligandcentre(3),ligandcoords(3),ligandcoordsrotated(3)
   double precision        :: twopi,pi,DPRAND,randomphi,randompsi,randomtheta,DIST,dummyz,mindistance
   double precision        :: st,ct,sph,cph,sps,cps, VECBAX, VECBAY, VECBAZ, randomx, randomy, randomz, dummyx, dummyy
   double precision        :: distancematrix(natoms,natoms), ysave(3*natoms),cartstepsave,transstepsave,localstep
   logical                 :: ligmovet,blockmovet,overlapt
   logical                 :: includedatom(natoms),enabledmove(natoms)
   character(len=10)       :: datechar,timechar,zonechar
   integer                 :: values(8),iostatus,restemp,i1,currentresidue,resnumber,overlaparray(natoms,natoms)
   integer                 :: overlaparraysave(natoms,natoms),nretries,loopcounter
   character(len=10)       :: rotmaxchangestr,rotpselectstr,rotcutoffstr,rotcentrestr,rotoccuwstr
   integer, allocatable    :: movableatomlist(:), atomsinblock(:) 
   double precision, allocatable :: atomblockcentre(:,:), rotationmatrix(:,:,:)

! sf344> for compatibility with dynamically determining macroion nearest neighbours and moving them with LIGMOVE
   if(allocated(movableatomlist)) deallocate(movableatomlist)
   allocate(movableatomlist(nmovableatoms))
   if(allocated(atomsinblock)) deallocate(atomsinblock)
   allocate(atomsinblock(nblocks))
   if(allocated(rotationmatrix)) deallocate(rotationmatrix)
   allocate(rotationmatrix(nblocks,3,3))
   if(allocated(atomblockcentre)) deallocate(atomblockcentre)
   allocate(atomblockcentre(nblocks,3))
   movableatomlist(:)=movableatomlist1(:)
   atomsinblock(:)=atomsinblock1(:)

   cartstepsave=ligcartstep
   transstepsave=ligtransstep

! parameters
   pi=ATAN(1.0D0)*4
   twopi=2.0D0*pi

! 1) RIGID BODY ROTATION

! sf344> work out for binary systems (macroion model) the bound counterions to each macroion, and move them cooperatively
!  with BLOCKMOVE.


! Macroions are the first in the atom list.
          nblocks=nmacroions
          if(allocated(movableatomlist)) deallocate(movableatomlist)
          if(allocated(atomsinblock)) deallocate(atomsinblock)
          if(allocated(atomblockcentre)) deallocate(atomblockcentre)
          if(allocated(rotationmatrix)) deallocate(rotationmatrix)
          allocate(movableatomlist(natoms),atomsinblock(nblocks),atomblockcentre(nblocks,3),rotationmatrix(nblocks,3,3))
          movableatomlist(:)=0
          atomsinblock(:)=0
          atomblockcentre(:,:)=0.0D0
          rotationmatrix(:,:,:)=0.0D0
          offset1=1
          includedatom(1:natoms)=.false.
          do i=1,nmacroions
               movableatomlist(offset1)=i
               offset1=offset1+1
               atomsinblock(i)=1
             do j=nmacroions+1,natoms
               distancematrix(i,j)=sqrt((y(3*i-2)-y(3*j-2))**2+(y(3*i-1)-y(3*j-1))**2+(y(3*i)-y(3*j))**2)
               if(distancematrix(i,j)<macroiondist.and..not.includedatom(j))then
                includedatom(j)=.true. ! include one counterion in only one block!
                movableatomlist(offset1)=j
                offset1=offset1+1
                atomsinblock(i)=atomsinblock(i)+1
               end if
             end do
          end do
!          write(MYUNITNEW,*) 'dynamically determining bound counterions to each macroion and moving them together'
!          write(MYUNITNEW,*) 'natoms, nblocks', natoms, nblocks
!          write(MYUNITNEW,*) 'atomsinblock(:)', atomsinblock(:)
!          write(MYUNITNEW,*) 'movableatomlist(:)', movableatomlist(:)
        ysave(:)=y(:)

   IF (LIGMOVET.AND.MOD(J1,ligmovefreq).EQ.0) THEN
        overlapt=.true.
        overlaparraysave(:,:)=0
        overlaparray(:,:)=0
        distancematrix(:,:)=0.0D0
        ysave(:)=y(:)
        IF(BLOCKMOVET) THEN
         ! Try to avoid moves that cause extremely high initial forces (or cold fusion),
         ! by redoing the rotation moves whenever two atoms get closer than 0.5 Angstroms.
          do i=1,natoms
           do j=1,natoms
              if(i==j) cycle
              distancematrix(i,j)=sqrt((y(3*i-2)-y(3*j-2))**2+(y(3*i-1)-y(3*j-1))**2+(y(3*i)-y(3*j))**2)
              if((MACROIONT.and.distancematrix(i,j)<1.0D0).or.distancematrix(i,j)<0.7D0) then 
                overlaparraysave(i,j)=1  !there's an initial overlap
!                write(*,*) 'atoms close together before rotating ', i,j,distancematrix(i,j)
              end if
           end do
          end do
        END IF

     nretries=0
     loopcounter=0
     overlapt=.true.
     DO WHILE(OVERLAPT)
       IF(nretries>10) then ! scale ligcartstep and ligtransstep
          nretries=0
!          ligcartstep=ligcartstep*0.5
!          ligtransstep=ligtransstep*0.5
!         WRITE(MYUNITNEW,*) 'scaling down ligcartstep and ligtransstep', ligcartstep, ligtransstep
       END IF
       IF(loopcounter>1000) exit !give up taking non-overlapping moves
       overlapt=.false.
       DO i=1,nblocks ! nblocks=1 if BLOCKMOVET=.FALSE. 
          randomphi=(DPRAND()-0.5)*twopi*ligrotscale*LOCALSTEP
          randomtheta=(DPRAND()-0.5)*pi*ligrotscale*LOCALSTEP
          randompsi=(DPRAND()-0.5)*twopi*ligrotscale*LOCALSTEP
          st=sin(randomtheta)
          ct=cos(randomtheta)
          sph=sin(randomphi)
          cph=cos(randomphi)
          sps=sin(randompsi)
          cps=cos(randompsi)

          rotationmatrix(i,1,1)=cps*cph-ct*sph*sps
          rotationmatrix(i,2,1)=cps*sph+ct*cph*sps
          rotationmatrix(i,3,1)=sps*st
          rotationmatrix(i,1,2)=-sps*cph-ct*sph*cps
          rotationmatrix(i,2,2)=-sps*sph+ct*cph*cps
          rotationmatrix(i,3,2)=cps*st
          rotationmatrix(i,1,3)=st*sph
          rotationmatrix(i,2,3)=-st*cph
          rotationmatrix(i,3,3)=ct
!        WRITE(*,*) 'rotation matrix: ', rotationmatrix(i,:,:)
       END DO

! work out the centre of coordinates for the ligand
       atomblockcentre(:,:)=0.0D0
       offset1=0
!      overlapt=.false.
       do i=1,nblocks
!          if(.not.enabledmove(i)) cycle
          if(i>1) offset1=offset1+atomsinblock(i-1)
          do k=1,atomsinblock(i)
             j=movableatomlist(k+offset1)
!            WRITE(*,*) k, j
             atomblockcentre(i,1)=atomblockcentre(i,1)+y(3*j-2)
             atomblockcentre(i,2)=atomblockcentre(i,2)+y(3*j-1)
             atomblockcentre(i,3)=atomblockcentre(i,3)+y(3*j  )
          end do
          atomblockcentre(i,:)=atomblockcentre(i,:)/atomsinblock(i)
          do k=1,atomsinblock(i)
             j=movableatomlist(k+offset1)
! move each block to the origin
             ligandcoords(1)=y(3*j-2)-atomblockcentre(i,1)
             ligandcoords(2)=y(3*j-1)-atomblockcentre(i,2)
             ligandcoords(3)=y(3*j  )-atomblockcentre(i,3)
! rotate the block of atoms with random rotation matrix
             ligandcoordsrotated=MATMUL(rotationmatrix(i,:,:),ligandcoords)
! translate back the block of atoms to the original centre of their coordinates
             y(3*j-2)=ligandcoordsrotated(1)+atomblockcentre(i,1)
             y(3*j-1)=ligandcoordsrotated(2)+atomblockcentre(i,2)
             y(3*j  )=ligandcoordsrotated(3)+atomblockcentre(i,3)
          end do            
       end do

! 2) RANDOM CARTESIAN PERTURBATION
       if((LIGMOVET).and.(ligcartstep.gt.0)) then
          do i=1,nmovableatoms
            if(frozen(i)) cycle
!            if(.not.enabledmove(i)) cycle
            j=movableatomlist(i)
            randomx=2*(DPRAND()-0.5D0)
            randomy=2*(DPRAND()-0.5D0)
            randomz=2*(DPRAND()-0.5D0)
            y(3*j-2)=y(3*j-2)+randomx*ligcartstep*LOCALSTEP
            y(3*j-1)=y(3*j-1)+randomy*ligcartstep*LOCALSTEP
            y(3*j  )=y(3*j  )+randomz*ligcartstep*LOCALSTEP
          enddo
       endif
! 3) RIGID TRANSLATION OF THE BLOCK OF ATOMS
       if((LIGMOVET).and.(ligtransstep.gt.0).and..not.BLOCKMOVET) then
          randomx=2*(DPRAND()-0.5D0)
          randomy=2*(DPRAND()-0.5D0)
          randomz=2*(DPRAND()-0.5D0)
          do i=1,nmovableatoms
            if(frozen(i)) cycle
!             if(.not.enabledmove(i)) cycle
             j=movableatomlist(i)
             y(3*j-2)=y(3*j-2)+randomx*ligtransstep*LOCALSTEP
             y(3*j-1)=y(3*j-1)+randomy*ligtransstep*LOCALSTEP
             y(3*j  )=y(3*j  )+randomz*ligtransstep*LOCALSTEP
          enddo
       endif
       if((LIGMOVET).and.(ligtransstep.gt.0).and.BLOCKMOVET) then
          offset1=0
!          WRITE(*,*) 'nblocks ', nblocks
          do i=1,nblocks
            if(frozen(i)) cycle
!             if(.not.enabledmove(i)) cycle
             if(i>1) offset1=offset1+atomsinblock(i-1)
             randomx=2*(DPRAND()-0.5D0)
             randomy=2*(DPRAND()-0.5D0)
             randomz=2*(DPRAND()-0.5D0)
             do k=1,atomsinblock(i)
                j=movableatomlist(k+offset1)
                y(3*j-2)=y(3*j-2)+randomx*ligtransstep*LOCALSTEP
                y(3*j-1)=y(3*j-1)+randomy*ligtransstep*LOCALSTEP
                y(3*j  )=y(3*j  )+randomz*ligtransstep*LOCALSTEP
             end do
          end do
       end if
       ! compute the new distances between all atoms and
       ! determine whether the overlap array has changed.
       do i=1,natoms
!          if(overlapt) exit
          do j=i+1,natoms
             distancematrix(i,j)=sqrt((y(3*i-2)-y(3*j-2))**2+(y(3*i-1)-y(3*j-1))**2+(y(3*i)-y(3*j))**2)
             if(i<=nblocks.and.j<=nblocks) then ! distance between two macroions
                mindistance=macroiondist/2
             else if((i<=nblocks.and.j>nblocks).or.(i>nblocks.and.j<=nblocks)) then ! distance between a macroion and a counterion
                mindistance=macroiondist/4
             else if(i>nblocks.and.j>nblocks) then  ! distance between two counterions
                mindistance=1.0D0
             end if
             if(distancematrix(i,j)<mindistance)then
                overlaparray(i,j)=1
!                overlapt=.true.
!                enabledmove(j)=.true.
             else
!                enabledmove(j)=.false. 
                overlaparray(i,j)=0 
             end if
             if(overlaparraysave(i,j)-overlaparray(i,j)<0) then 
                overlapt=.true.
!                write(*,*) 'atoms close together after rotation, random perturbation and group translation: ', i,j,distancematrix(i,j)
!                y(:)=ysave(:)
                nretries=nretries+1
                loopcounter=loopcounter+1
!                exit 
             else
                overlapt=.false. 
             end if
          end do
       end do

       if(overlapt) cycle
 ! now we finally have no overlap, so check for percolation, and take a new step until we start off from a percolated structure
       PERCT=.TRUE.
       CALL PERC(y,NATOMS,PERCCUT,PERCT,.FALSE.,MYUNITNEW,.FALSE.)
!      WRITE(*,*) 'perccut=', perccut
       IF(.NOT.PERCT) then
!          write(MYUNITNEW,*) ' macroion_moves> disconnected structure detected, undoing step'
          y(:)=ysave(:)
          overlapt=.true.
          nretries=nretries+1  ! make sure we are not getting stuck in an infinite loop
          loopcounter=loopcounter+1
          cycle
       ELSE
          overlapt=.false.
       END IF
     END DO  ! while(overlapt)
   END IF
! reset ligcartstep and ligtransstep
   ligcartstep=cartstepsave
   ligtransstep=transstepsave

END SUBROUTINE MACROION_MOVES

END MODULE MOVES
