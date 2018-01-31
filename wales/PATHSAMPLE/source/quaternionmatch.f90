!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999-2009 David J. Wales
!   This file is part of PATHSAMPLE.
!
!   PATHSAMPLE is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   PATHSAMPLE is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!

! Superposes structure B onto structure A using the quaternion method
! Neither structure is changed
! Also returns the RMSD, the rotation and the translation to perform (in that order) to achieve 
! the superposition

subroutine quaternionMatch (natoms, structureA, structureB, minRMSD, rotation, translation)

  implicit none

! Subroutine arguments

  integer, intent(IN) :: natoms
  real (kind=kind(0.0d0)), intent(IN) :: structureA(NATOMS,3)
  real (kind=kind(0.0d0)), intent(IN) :: structureB(NATOMS,3)
  real (kind=kind(0.0d0)), intent(OUT) :: minRMSD
  real (kind=kind(0.0d0)), intent(OUT) :: translation(3), rotation(3)
  
! Local variables

  real (kind=kind(0.0d0)) :: eigenValues(4)
  real (kind=kind(0.0d0)) :: centroid1(3), centroid2(3), quaternionMatrix(4,4)
  real (kind=kind(0.0d0)) :: differenceMatrix(natoms,3), sumMatrix(natoms,3)
  integer :: currentAtom
  real (kind=kind(0.0d0)), parameter :: zeroThreshold = 1d-20   ! For identifying identical structures  

! DSYEV stuff
  
  integer :: INFO
  real (kind=kind(0.0d0)) :: WORK(12) 

! Calculate centroids

  centroid1 = 0.0d0
  centroid2 = 0.0d0
  
  do currentAtom = 1, natoms
     centroid1 = centroid1 + structureA(currentAtom,:)
     centroid2 = centroid2 + structureB(currentAtom,:)
  enddo

  centroid1 = centroid1/natoms
  centroid2 = centroid2/natoms 

!  print *, 'Centroids:'
!  print *, 'System 1:', (centroid1(currentAtom),currentAtom=1,3)
!  print *, 'System 2:', (centroid2(currentAtom),currentAtom=1,3)

! Calculate matrices which will provide the elements of the quaternions

  differenceMatrix = 0.0d0
  sumMatrix = 0.0d0

  do currentAtom = 1, natoms
     sumMatrix(currentAtom,:) = structureA(currentAtom,:) - centroid1 + &
                                structureB(currentAtom,:) - centroid2

     differenceMatrix(currentAtom,:) =  structureA(currentAtom,:) - centroid1 - &
                                        (structureB(currentAtom,:) - centroid2)    
  enddo

! Calculate quaternion matrix  

  quaternionMatrix = 0.0d0

  do currentAtom = 1, natoms
    quaternionMatrix(1,1) = quaternionMatrix(1,1) + differenceMatrix(currentAtom,1)**2 + &
                             differenceMatrix(currentAtom,2)**2 + differenceMatrix(currentAtom,3)**2

     quaternionMatrix(2,1) = quaternionMatrix(2,1) + sumMatrix(currentAtom,2)*differenceMatrix(currentAtom,3) - &
                             differenceMatrix(currentAtom,2)*sumMatrix(currentAtom,3)

     quaternionMatrix(2,2) = quaternionMatrix(2,2) + differenceMatrix(currentAtom,1)**2 + &
                             sumMatrix(currentAtom,2)**2 + sumMatrix(currentAtom,3)**2  

     quaternionMatrix(3,1) = quaternionMatrix(3,1) + sumMatrix(currentAtom,3)*differenceMatrix(currentAtom,1) - &
                             differenceMatrix(currentAtom,3)*sumMatrix(currentAtom,1)     

     quaternionMatrix(3,2) = quaternionMatrix(3,2) + differenceMatrix(currentAtom,1)*differenceMatrix(currentAtom,2) - &
                             sumMatrix(currentAtom,1)*sumMatrix(currentAtom,2)     

     quaternionMatrix(3,3) = quaternionMatrix(3,3) + sumMatrix(currentAtom,1)**2 + &
                             differenceMatrix(currentAtom,2)**2 + sumMatrix(currentAtom,3)**2 

     quaternionMatrix(4,1) = quaternionMatrix(4,1) + sumMatrix(currentAtom,1)*differenceMatrix(currentAtom,2) - &
                             differenceMatrix(currentAtom,1)*sumMatrix(currentAtom,2)

     quaternionMatrix(4,2) = quaternionMatrix(4,2) + differenceMatrix(currentAtom,1)*differenceMatrix(currentAtom,3) - &
                             sumMatrix(currentAtom,1)*sumMatrix(currentAtom,3)

     quaternionMatrix(4,3) = quaternionMatrix(4,3) + differenceMatrix(currentAtom,2)*differenceMatrix(currentAtom,3) - &
                             sumMatrix(currentAtom,2)*sumMatrix(currentAtom,3)     

     quaternionMatrix(4,4) = quaternionMatrix(4,4) + sumMatrix(currentAtom,1)**2 + &
                             sumMatrix(currentAtom,2)**2 + differenceMatrix(currentAtom,3)**2
  enddo

! Diagonalise then sort eigenvectors/values

  call DSYEV('V','L',4,quaternionMatrix,4,eigenValues,WORK,size(WORK),INFO)

  if (INFO.ne.0) then
     print *, 'Problem occurred in DSYEV'
     stop
  endif

! Despite the name actually sorts into descending order
  call eigensort_val_asc(eigenValues,quaternionMatrix,4,4)

  if (eigenValues(4).lt.-1.0d-12) then
     print *, 'Problem in quaternion match routine: significantly negative RMSD:', eigenValues(4)
     stop
  end if

! RMSD and rotation vector are straightforward

  minRMSD = sqrt(max(0.0d0,eigenValues(4))/natoms)

  call quaternionToAngleAxis(quaternionMatrix(:,4),rotation)

! Work out effect of rotating about reference origin on reference centroid

  call rotateAngleAxis(centroid2, rotation)
  translation = centroid1 - centroid2


  return

end subroutine quaternionMatch


! this returns the square of the euclidean distance, not the rmsd!
! Edist**2 =  RMSD**2 * natoms
subroutine quaternionMatch2 (natoms, structureA, structureB, distance)

  implicit none

! Subroutine arguments

  integer, intent(IN) :: natoms
  real (kind=kind(0.0d0)), intent(IN) :: structureA(3,NATOMS)
  real (kind=kind(0.0d0)), intent(IN) :: structureB(3,NATOMS)
  real (kind=kind(0.0d0)), intent(OUT) :: distance
  
! Local variables

  real (kind=kind(0.0d0)) :: eigenValues(4), Rmatrix(3,3)
  real (kind=kind(0.0d0)) :: centroid1(3), centroid2(3), quaternionMatrix(4,4)
  real (kind=kind(0.0d0)) :: sA(3,natoms), sB(3,natoms)
  integer :: currentAtom, i, j
  real (kind=kind(0.0d0)), parameter :: zeroThreshold = 1d-20   ! For identifying identical structures  

! DSYEV stuff
  
  integer :: INFO
  real (kind=kind(0.0d0)) :: WORK(12) 

! Calculate centroids

  centroid1 = sum(structureA, 2)/natoms
  centroid2 = sum(structureB, 2)/natoms

!  print *, 'Centroids:'
!  print *, 'System 1:', (centroid1(currentAtom),currentAtom=1,3)
!  print *, 'System 2:', (centroid2(currentAtom),currentAtom=1,3)


  do currentAtom = 1, natoms
     sA(:, currentAtom) = structureA(:, currentAtom) - centroid1
     sB(:, currentAtom) = structureB(:, currentAtom) - centroid2
  enddo

  Rmatrix = matmul(sA, transpose(sB))

! Calculate quaternion matrix  

  quaternionMatrix = 0.0d0

  quaternionMatrix(1, 1) = Rmatrix(1, 1) + Rmatrix(2, 2) + Rmatrix(3, 3)
  quaternionMatrix(2, 1) = Rmatrix(2, 3) - Rmatrix(3, 2)
  quaternionMatrix(3, 1) = Rmatrix(3, 1) - Rmatrix(1, 3)
  quaternionMatrix(4, 1) = Rmatrix(1, 2) - Rmatrix(2, 1)

  quaternionMatrix(1, 2) = quaternionMatrix(2, 1)
  quaternionMatrix(2, 2) = Rmatrix(1, 1) - Rmatrix(2, 2) - Rmatrix(3, 3)
  quaternionMatrix(3, 2) = Rmatrix(1, 2) + Rmatrix(2, 1)
  quaternionMatrix(4, 2) = Rmatrix(1, 3) + Rmatrix(3, 1)

  quaternionMatrix(1, 3) = quaternionMatrix(3, 1)
  quaternionMatrix(2, 3) = quaternionMatrix(3, 2)
  quaternionMatrix(3, 3) =-Rmatrix(1, 1) + Rmatrix(2, 2) - Rmatrix(3, 3)
  quaternionMatrix(4, 3) = Rmatrix(2, 3) + Rmatrix(3, 2)

  quaternionMatrix(1, 4) = quaternionMatrix(4, 1)
  quaternionMatrix(2, 4) = quaternionMatrix(4, 2)
  quaternionMatrix(3, 4) = quaternionMatrix(4, 3)
  quaternionMatrix(4, 4) =-Rmatrix(1, 1) - Rmatrix(2, 2) + Rmatrix(3, 3) 

! Diagonalise then sort eigenvectors/values

  call DSYEV('V','L',4,quaternionMatrix,4,eigenValues,WORK,size(WORK),INFO)

  if (INFO.ne.0) then
     print *, 'Problem occurred in DSYEV'
     stop
  endif

! Despite the name actually sorts into descending order
  call eigensort_val_asc(eigenValues,quaternionMatrix,4,4)

  if (eigenValues(1).lt.-1.0d-12) then
     print *, 'Problem in quaternion match routine: significantly negative RMSD:', eigenValues(4)
     stop
  end if

  distance= max(0.0d0, sum(sA**2)+sum(sB**2)-2*eigenValues(1))

  return

end subroutine quaternionMatch2



! #########################################################################################

! Converts a quaternion into a rotation vector

subroutine quaternionToAngleAxis (quaternion,rotationVector)

  implicit none
  
  real (kind=kind(0.0d0)), intent(IN) :: quaternion(4)
  real (kind=kind(0.0d0)), intent(OUT) :: rotationVector(3)
    
  real (kind=kind(0.0d0)) :: cosAngle, sinAngle, angle, initialModulus

! Angle stuff first

  cosAngle = (3.0d0*quaternion(1)**2 - quaternion(2)**2 - &
              quaternion(3)**2 - quaternion(4)**2 - 1.0d0)/2.0d0

  if (cosAngle > 1.0d0) then
     cosAngle = 1.0d0
  elseif (cosAngle < -1.0d0) then
     cosAngle = -1.0d0
  endif

  angle = acos(cosAngle)
  sinAngle = dsqrt(1.0d0 - cosAngle**2)
! Unit vector
! sf344> sinAngle very rarely can become zero, and we don't want to divide by that! (quick dirty fix)
if(sinAngle/=0.0) then
  rotationVector(1) = 2.0d0*quaternion(1)*quaternion(2)/sinAngle
  rotationVector(2) = 2.0d0*quaternion(1)*quaternion(3)/sinAngle
  rotationVector(3) = 2.0d0*quaternion(1)*quaternion(4)/sinAngle
else 
  rotationVector(1) = 2.0d0*quaternion(1)*quaternion(2)
  rotationVector(2) = 2.0d0*quaternion(1)*quaternion(3)
  rotationVector(3) = 2.0d0*quaternion(1)*quaternion(4)
end if  
  initialModulus = dsqrt(rotationVector(1)**2 + rotationVector(2)**2 + rotationVector(3)**2)
  
! Combine the two

  rotationVector = rotationVector*(angle/initialModulus)
    
  return

end subroutine quaternionToAngleAxis

! #########################################################################################

subroutine eigensort_val_asc(eval, evec, dof1, dof2)
      implicit none

      integer,intent(in)         :: dof1, dof2
      DOUBLE PRECISION,INTENT(INOUT)  :: EVAL(DOF1), EVEC(DOF1,DOF2)

      integer :: i,j,k
      DOUBLE PRECISION :: TMP 

      do i=1,dof1-1
         k=i
         tmp=eval(i)
         do j=i+1, dof1
            if (eval(j)>=tmp) then
               k=j
               tmp=eval(j)
            endif
         enddo
         if (.not.k==i) then
            eval(k)=eval(i)
            eval(i)=tmp
            do j=1,dof1
               tmp = evec(j,i)
               evec(j,i) = evec(j,k)
               evec(j,k) = tmp
            enddo
         endif
      enddo
 end subroutine eigensort_val_asc

! #########################################################################################


 subroutine rotateAngleAxis(point, rotationAxis)

   implicit none

! Subroutine arguments

   real (kind=kind(0.0d0)), intent(IN OUT) :: point(3)
   real (kind=kind(0.0d0)), intent(IN) :: rotationAxis(3)

! Local variables
  
   real (kind=kind(0.0d0)) :: angle, cosAngle, cosFactor, sinFactor, pdotX, rotatedPoint(3)

   real (kind=kind(0.0d0)) :: vectorModulus  

   angle = vectorModulus(3,rotationAxis)
   cosAngle = cos(angle)
   cosFactor = (1.0d0-cosAngle)/(angle**2)
   sinFactor = sin(angle)/angle

   pdotX = dot_product(rotationAxis, point)            
   
   rotatedPoint(1) = point(1)*cosAngle + &
                     rotationAxis(1)*pdotX*cosFactor + &
                     (point(2)*rotationAxis(3) - &
                      point(3)*rotationAxis(2))*sinFactor

   rotatedPoint(2) = point(2)*cosAngle + &
                     rotationAxis(2)*pdotX*cosFactor + &
                     (point(3)*rotationAxis(1) - &
                      point(1)*rotationAxis(3))*sinFactor

   rotatedPoint(3) = point(3)*cosAngle + &
                     rotationAxis(3)*pdotX*cosFactor + &
                     (point(1)*rotationAxis(2) - &
                      point(2)*rotationAxis(1))*sinFactor

   point = rotatedPoint
     
   return

 end subroutine rotateAngleAxis

! #########################################################################################

 function vectorModulus(nElements, vector)
   implicit none
      
   real (kind=kind(0.0d0)) :: vectorModulus   
   integer, intent(IN) :: nElements
   real (kind=kind(0.0d0)), intent(IN) :: vector(nElements)
   integer :: currentElement

   vectorModulus = 0.0d0
   do currentElement = 1, nElements
      vectorModulus = vectorModulus + vector(currentElement)**2.0d0
   enddo

   vectorModulus = sqrt(vectorModulus)
   return
 end function vectorModulus

! #########################################################################################
