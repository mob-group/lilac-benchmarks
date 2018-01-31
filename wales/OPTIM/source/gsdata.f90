MODULE GSDATA
  
  USE COMMONS
  USE INTCOMMONS, ONLY : NINTC, DESMINT

  IMPLICIT NONE
  SAVE

  TYPE IMGNODE
     DOUBLE PRECISION, POINTER :: XYZ(:), GRAD(:), TGT(:), XINT(:)
     DOUBLE PRECISION, POINTER :: PREVGRAD(:), DIFF(:), PREVXYZ(:)
     DOUBLE PRECISION, POINTER :: XYZCART(:), GCART(:), PREVXCART(:)
     DOUBLE PRECISION :: TNORM, FNORM2, E !energy
     DOUBLE PRECISION :: CHORD, ARC, ARCDIFF
     LOGICAL :: LF, RF ! is this left or right front?
     LOGICAL :: FREEZE
     INTEGER :: IND ! index

     TYPE(IMGNODE), POINTER :: NEXT
     TYPE(IMGNODE), POINTER :: PREV
  ENDTYPE IMGNODE  

  DOUBLE PRECISION, PARAMETER :: TINY=1.0D-8

  ! keywords
  INTEGER :: GSUPDATE, MAXGROWSTEPS
  LOGICAL :: EVOLVESTRINGT
  LOGICAL :: CUBSPLT,FIXATMS, PREROTATE, DUMPGSALL = .FALSE.
  INTEGER :: TANTYPE
  DOUBLE PRECISION :: GSMXSTP, GSMAXTOTITD
  DOUBLE PRECISION :: GSGROWTOL, GSCONV, REPARAMTOL
  DOUBLE PRECISION :: GSIMGDENSITY, GSITERDENSITY, GSDGUESS
  DOUBLE PRECISION :: MAXLENPERIM
  LOGICAL :: NOLBFGS

  INTEGER :: NC ! number of coordinates
  LOGICAL :: PTEST, JOINED, HESSGRAD, INTPTEST
  DOUBLE PRECISION :: STRINGLEN ! total string length
  DOUBLE PRECISION :: MAXLEN
  INTEGER :: NIM, NLIM, NRIM ! number of images on left and right sides
  INTEGER :: TIM ! total number of images
  INTEGER :: MAXTOTSTEPS, TOTSTEPS ! total steps taken so far

  ! pointers to the first node, last node, left and right frontiers
  ! first and last point to actual images not terminal configurations
  ! (ie: there's actually one node before first and one after last)
  TYPE(IMGNODE), POINTER :: FIRST, LAST, LEFTFRONT, RIGHTFRONT

  DOUBLE PRECISION :: IMGD, ITERD ! current image and iteration densities
  INTEGER :: MAXDROPS, M ! M is the lbfgs memory
  DOUBLE PRECISION :: MAXERISE  

  INTEGER :: DUMPIND

CONTAINS
  SUBROUTINE KEYGSPRINT(VARIABLE)
    USE KEY, ONLY : DESMAXAVGE, DESMAXEJUMP
    IMPLICIT NONE
    LOGICAL :: VARIABLE ! is the number of images variable?

    IF(VARIABLE) THEN
       WRITE(*,'(1x,a)') 'KeyGS> Number of images will vary depending on the separation of the endpoints'
       WRITE(*,'(1x,a,F10.5)') 'KeyGS> Starting iteration density per image: ', GSITERDENSITY
    ELSE
       WRITE(*,'(1x,a,F10.5)') 'KeyGS> iteration density per image: ', ITERD
       WRITE(*,'(1x,a,I4,a)') 'KeyGS> Using ', TIM, ' images'
    ENDIF

    IF (NOLBFGS) WRITE(*,'(1x,a)') 'KeyGS> No LBFGS optimization. Move images according to FPERP'
    IF (EVOLVESTRINGT) WRITE(*,'(1x,a)') 'KeyGS> Using evolving string method.&
         & String will start out populated with images'
    IF (CUBSPLT) WRITE(*,'(1x,a)') 'KeyGS> Using cubic spline interpolation between images'
    IF (FIXATMS) WRITE(*,'(1x,a)') 'KeyGS> Projecting out overall translation and rotations &
         & by fixing coordinates 1-4 and 7-8'
    IF(HESSGRAD) WRITE(*,'(1x,a)') 'KeyGS> Using method described in Appendix of &
         & Peters et al to generate the Newton-Raphston step; Hessian determined by &
         & changing gradient rather than changing perpendicular force; tangential &
         & component projected after multiplying by Hessian'
    WRITE(*,'(1x,a,I2)') 'KeyGS> Tangent type: ', TANTYPE
    WRITE(*,'(1x,a,2F10.5)') 'KeyGS> Reparametrization tolerance, growth tolerance: ', &
         & REPARAMTOL, GSGROWTOL
    WRITE(*,'(1x,a,F10.5)') 'KeyGS> Convergence tolerance: ', GSCONV
    WRITE(*,'(1x,a,F10.5,I4)') 'KeyGS> Max step size, LBFGS memory: ', GSMXSTP, GSUPDATE
    WRITE(*,'(1x,a,I6)') 'KeyGS> Max growth steps: ', MAXGROWSTEPS 
    IF (DESMAXAVGE.LT.0.99*HUGE(1.0D0)) WRITE(*,'(1x,a,G20.10)') 'KeyGS> Max average energy: ', DESMAXAVGE
    IF (DESMAXEJUMP.LT.0.99*HUGE(1.0D0)) WRITE(*,'(1x,a,G20.10)') 'KeyGS> Max energy jump per image: ', DESMAXEJUMP
    IF(GSMAXTOTITD.GE.0) WRITE(*,'(1x,a,F10.1)') 'KeyGS> Max total iteration density: ', GSMAXTOTITD    
  END SUBROUTINE KEYGSPRINT

  SUBROUTINE NEWNODE(P)
    
    IMPLICIT NONE    
    TYPE(IMGNODE), POINTER :: P
    
    ! make a new imgnode and return a pointer to that node
    ! assumes previous association of P is checked elsewhere and it's ok to allocate it
    
    ALLOCATE(P)
    NULLIFY(P%NEXT, P%PREV, P%XYZ, P%GRAD, P%TGT, P%PREVGRAD, P%DIFF)
    
    ALLOCATE(P%XYZ(NC), P%GRAD(NC), P%TGT(NC), &
            & P%PREVGRAD(NC), P%DIFF(NC), P%PREVXYZ(NC))
    IF (DESMINT) ALLOCATE(P%XYZCART(3*NATOMS), P%GCART(3*NATOMS), P%PREVXCART(3*NATOMS))
       
    P%ARC = -1; P%CHORD = -1
    P%RF = .FALSE.; P%LF = .FALSE.
    P%FREEZE = .FALSE.

  END SUBROUTINE NEWNODE
  
  SUBROUTINE DELETENODE(P)
    IMPLICIT NONE
    
    TYPE(IMGNODE), POINTER :: P
    
    ! deallocate the imagenode to which P is pointing and nullify p
    
    DEALLOCATE(P%XYZ, P%GRAD, P%TGT, P%PREVGRAD, P%DIFF, P%PREVXYZ)
    IF (DESMINT) DEALLOCATE(P%XYZCART, P%GCART, P%PREVXCART)
    DEALLOCATE(P)
    NULLIFY(P)
    
  END SUBROUTINE DELETENODE
  
  SUBROUTINE EXTENDLEFT
    
    ! extend left side of string    

    LEFTFRONT%LF = .FALSE.

    CALL NEWNODE(LEFTFRONT%NEXT)    
    LEFTFRONT%NEXT%PREV => LEFTFRONT
    LEFTFRONT => LEFTFRONT%NEXT
    LEFTFRONT%NEXT => RIGHTFRONT
    RIGHTFRONT%PREV => LEFTFRONT

    LEFTFRONT%LF = .TRUE.
  END SUBROUTINE EXTENDLEFT
  
  SUBROUTINE EXTENDRIGHT   
    ! extend right side of string    

    RIGHTFRONT%RF = .FALSE.

    CALL NEWNODE(RIGHTFRONT%PREV)
    RIGHTFRONT%PREV%NEXT => RIGHTFRONT
    RIGHTFRONT => RIGHTFRONT%PREV
    RIGHTFRONT%PREV => LEFTFRONT
    LEFTFRONT%NEXT => RIGHTFRONT

    RIGHTFRONT%RF = .TRUE.
  END SUBROUTINE EXTENDRIGHT

  SUBROUTINE DELETESTRING
    ! delete entire string, from FIRST%PREV to LAST%NEXT

    IMPLICIT NONE
    TYPE(IMGNODE), POINTER :: DUMMYP

    DUMMYP => FIRST%PREV
    DO WHILE (ASSOCIATED(DUMMYP%NEXT))
       DUMMYP => DUMMYP%NEXT
       CALL DELETENODE(DUMMYP%PREV)
    ENDDO
    CALL DELETENODE(DUMMYP)

    NULLIFY(FIRST, LAST, RIGHTFRONT, LEFTFRONT)
  END SUBROUTINE DELETESTRING
END MODULE GSDATA

