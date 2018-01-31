MODULE GENRIGID

      INTEGER :: NRIGIDBODY, DEGFREEDOMS, MAXSITE, NRELAXRIGIDR, NRELAXRIGIDA
      INTEGER :: XNATOMS
      INTEGER, ALLOCATABLE :: NSITEPERBODY(:), REFVECTOR(:), RIGIDSINGLES(:)
      INTEGER, ALLOCATABLE, DIMENSION (:,:) :: RIGIDGROUPS
      DOUBLE PRECISION, ALLOCATABLE :: RIGIDCOORDS(:)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: SITESRIGIDBODY
      DOUBLE PRECISION, ALLOCATABLE :: GR_WEIGHTS(:) ! weights for com calculation, e.g. masses
      LOGICAL :: RIGIDINIT, ATOMRIGIDCOORDT, RELAXRIGIDT
      LOGICAL :: GENRIGIDT

      LOGICAL :: RIGIDOPTIMROTAT, FREEZERIGIDBODYT
      DOUBLE PRECISION :: OPTIMROTAVALUES(3)

!   vr274:  added lattice coordinates
!           if HAS_LATTICE_COORDS is true, the last two atoms are treated
!           as lattice coordintes and rigidcoords is in reduced lattice units
      LOGICAL HAS_LATTICE_COORDS

!-----------------------------------------------------------------------------------!
! NRIGIDBODY  = number of rigid bodies
! DEGFREEDOMS = number of degrees of freedom = 6 * NRIGIDBODY + 3 * ADDITIONAL ATOMS
! MAXSITE     = maximum number of sites in a rigid body
! NRELAXRIGIDR = rigid body minimisation for this number of steps
! NRELAXRIGIDA = atom minimisation for this number of steps
! NSITEPERBODY= number of rigid body sites, no need to be the same for all bodies
! REFVECTOR   = reference vector for the atomistic to rigic coordinate transformation
! RIGIDSINGLES= list of atoms not in rigid bodies
! RIGIDGROUPS = list of atoms in rigid bodies, need a file called rbodyconfig
! RIGIDCOORDS = 6 * NRIGIDBODY + 3 * ADDITIONAL ATOMS coordinates
! SITESRIGIDBODY = coordinates of the rigid body sites
! RIGIDINIT   = logical variable for generalised rigid body
! ATOMRIGIDCOORDT, .TRUE. = atom coords active, .FALSE. = rigid coords active, used in mylbfgs & potential
! GENRIGIDT = generalised rigid body takestep taken if .TRUE.
!-----------------------------------------------------------------------------------!


CONTAINS
! vr274> TODO: better initialization of SITESRIDIGBODY, why use maxsite?

!-------------------------------------------
! vr274> Initializes basic structures
! hast to be the first call to GENRIGID function in order to setup basic structures.
! After that, the array which defines the sites can be filled. Then GENRIGID_INITIALIZE
! completes the initialization of rigid bodies.
!-------------------------------------------
SUBROUTINE GENRIGID_ALLOCATE(NEW_NRIGIDBODY,NEW_MAXSITE)
  USE COMMONS, only: NATOMS
  IMPLICIT NONE
  INTEGER, intent(in) :: NEW_NRIGIDBODY, NEW_MAXSITE
  
  NRIGIDBODY = NEW_NRIGIDBODY
  MAXSITE = NEW_MAXSITE

  ! hk286 > Allocate NSITEPERBODY
  ALLOCATE (NSITEPERBODY(NRIGIDBODY))
  ALLOCATE (SITESRIGIDBODY(MAXSITE,3,NRIGIDBODY))
  ALLOCATE (RIGIDGROUPS(MAXSITE,NRIGIDBODY))
  ALLOCATE (REFVECTOR(NRIGIDBODY))
  ALLOCATE (GR_WEIGHTS(NATOMS))

  ! by default use center of geometry
  GR_WEIGHTS=1d0
END SUBROUTINE

!-------------------------------------------
! vr274> Setup rigid body stuff after site definitions are done
!-------------------------------------------
SUBROUTINE GENRIGID_INITIALISE(INICOORDS)
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  INTEGER :: J1, J2, J3, DUMMY
  DOUBLE PRECISION :: XMASS, YMASS, ZMASS, PNORM, MASS
  LOGICAL :: SATOMT, RTEST
  DOUBLE PRECISION :: P(3), RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
  DOUBLE PRECISION INICOORDS(3*NATOMS)
  
! vr275> initialize coordinates for rigid bodies
  DO J1 = 1, NRIGIDBODY
     DO J2 = 1, NSITEPERBODY(J1)
        DUMMY=RIGIDGROUPS(J2,J1)
        SITESRIGIDBODY(J2,:,J1) = INICOORDS(3*DUMMY-2:3*DUMMY)
     ENDDO
  ENDDO


  ! hk286 > determine number of degrees of freedom
  DEGFREEDOMS = 0
  DO J1 = 1, NRIGIDBODY
     DEGFREEDOMS = DEGFREEDOMS + NSITEPERBODY(J1)
  ENDDO
  DEGFREEDOMS = 6 * NRIGIDBODY + 3 * (NATOMS - DEGFREEDOMS)

! hk286 > Allocate further data 
  ALLOCATE (RIGIDSINGLES((DEGFREEDOMS/3 - 2 * NRIGIDBODY)))
  ALLOCATE (RIGIDCOORDS(DEGFREEDOMS))

  DUMMY = 0
  DO J1 = 1, NATOMS
     SATOMT = .TRUE.
     DO J2 = 1, NRIGIDBODY
        DO J3 = 1, NSITEPERBODY(J2)
           IF (J1 == RIGIDGROUPS(J3,J2)) SATOMT = .FALSE.
        ENDDO
     ENDDO
     IF (SATOMT) THEN
        DUMMY = DUMMY + 1
        RIGIDSINGLES(DUMMY) = J1
     ENDIF
  ENDDO
  
  DO J1 = 1, NRIGIDBODY
     XMASS = 0.0D0
     YMASS = 0.0D0
     ZMASS = 0.0D0
     MASS = 0.0d0
     DO J2 = 1, NSITEPERBODY(J1)
        XMASS = XMASS + SITESRIGIDBODY(J2,1,J1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        YMASS = YMASS + SITESRIGIDBODY(J2,2,J1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        ZMASS = ZMASS + SITESRIGIDBODY(J2,3,J1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        MASS  = MASS + GR_WEIGHTS(RIGIDGROUPS(J2,J1))
     ENDDO
     XMASS = XMASS / MASS
     YMASS = YMASS / MASS
     ZMASS = ZMASS / MASS
     DO J2 = 1, NSITEPERBODY(J1)
        SITESRIGIDBODY(J2,1,J1) = SITESRIGIDBODY(J2,1,J1) - XMASS
        SITESRIGIDBODY(J2,2,J1) = SITESRIGIDBODY(J2,2,J1) - YMASS
        SITESRIGIDBODY(J2,3,J1) = SITESRIGIDBODY(J2,3,J1) - ZMASS
     ENDDO

  ENDDO

!  PRINT *, SITESRIGIDBODY
  
! hk286 > make sure the two atoms used as reference for rigid bodies are suitable
! Checks: (1) Atoms 1 and 2 do not sit on COM, and (2) Vector 1 and 2 are not parallel
  
  DO J1 = 1, NRIGIDBODY
     REFVECTOR(J1) = 1
     RTEST = .TRUE.
     DO WHILE (RTEST)
        RTEST = .FALSE.
        DO J2 = REFVECTOR(J1), REFVECTOR(J1) + 1 
           PNORM = SQRT(DOT_PRODUCT(SITESRIGIDBODY(J2,:,J1),SITESRIGIDBODY(J2,:,J1)))
           IF ( (PNORM  < 0.001) .AND. (PNORM > -0.001)) THEN
              RTEST = .TRUE.
           ENDIF
        ENDDO
        PNORM = DOT_PRODUCT(SITESRIGIDBODY(REFVECTOR(J1),:,J1),SITESRIGIDBODY(REFVECTOR(J1)+1,:,J1)) 
        PNORM = PNORM / SQRT(DOT_PRODUCT(SITESRIGIDBODY(REFVECTOR(J1),:,J1),SITESRIGIDBODY(REFVECTOR(J1),:,J1))) 
        PNORM = PNORM / SQRT(DOT_PRODUCT(SITESRIGIDBODY(REFVECTOR(J1)+1,:,J1),SITESRIGIDBODY(REFVECTOR(J1)+1,:,J1)))
        IF (PNORM < 0.0) PNORM = -1.0D0 * PNORM
        IF ( (PNORM < 1.0 + 0.001) .AND. (PNORM > 1.0 - 0.001) ) THEN
           RTEST = .TRUE.
        ENDIF
        IF (RTEST) THEN
           REFVECTOR(J1) = REFVECTOR(J1) + 1               
        ENDIF
     ENDDO
  ENDDO


! hk286 - new 7/1/12
  IF (RIGIDOPTIMROTAT .EQV. .TRUE.) THEN
     CALL ROTATEINITIALREF ()
  ENDIF

! list of useful checks
!  P(1) = 0.0001D0
!  P(2) = 0.0003D0
!  P(3) = 0.0002D0
!  CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, .TRUE.)

  !PRINT *, "TEST"
  !ATOMRIGIDCOORDT = .FALSE.
  !AMBERT = .TRUE.
  !CALL TRANSFORMCTORIGID(COORDS, RIGIDCOORDS)
  !COORDS(1:DEGFREEDOMS,1) = RIGIDCOORDS(1:DEGFREEDOMS)
  !COORDS(DEGFREEDOMS+1:3*NATOMS,1) = 0.0D0
  !CALL POTENTIAL(COORDS,GRAD,ENERGY,.TRUE.,.FALSE.)
  !PRINT *, GRAD
  !CALL TRANSFORMGRAD(GRAD,RIGIDCOORDS,GRADR)
  !PRINT*, ENERGY
 
  
!      PRINT *, REFVECTOR
!      CHECKS  
!      PRINT*, COORDS(:,1)
!      PRINT*, " "
!      CALL TRANSFORMCTORIGID(COORDS(:,1), RIGIDCOORDS)
!      PRINT*, RIGIDCOORDS
!      PRINT*, " "
!      CALL TRANSFORMRIGIDTOC(1,NRIGIDBODY, COORDS(:,1), RIGIDCOORDS)
!      PRINT*, COORDS(:,1)
!      PRINT*, " "
!      READ *, DUMMY
!      PRINT*, NRIGIDBODY
!      PRINT*, NSITEPERBODY
!      PRINT*, NATOMS

!      DO J1 = 1, 100
!         CALL TAKESTEPGENRIGID ()
!         PRINT*, " "
!         PRINT*, COORDS(:,1)
!         CALL TRANSFORMCTORIGID()
!         CALL TRANSFORMRIGIDTOC(1,2)
!         PRINT*, COORDS(:,1)
!      ENDDO
!      NATOMS = 4
!      CALL NEWCAPSID(RIGIDCOORDS,GRADR,ENERGY, .TRUE.)
!      NATOMS = 12
!      PRINT*, ENERGY
!      PRINT*, GRADR(1), GRADR(2), GRADR(3), GRADR(4), GRADR(5), GRADR(6)
!      PRINT*, GRADR(7), GRADR(8), GRADR(9), GRADR(10), GRADR(11), GRADR(12)

!      CALL GENCAPSID(COORDS(:,1),GRAD,ENERGY,.TRUE.)
!      CALL TRANSFORMGRAD(GRAD,RIGIDCOORDS,GRADR)
!      PRINT*, ENERGY
!      PRINT*, GRADR
!      PRINT*, " "
!      PRINT*, " "
!      PRINT*, NRIGIDBODY
!      PRINT*, NSITEPERBODY
!      PRINT*, NATOMS

      !READ*, DUMMY
    
!      CALL POTENTIALLBFGS(RIGIDCOORDS,GRADR,ENERGY,.TRUE.,.FALSE.,NATOMS)
!      PRINT*, ENERGY
!      PRINT*, RIGIDCOORDS
!      PRINT*, GRADR
!      PRINT*, " "

      !CALL VIEWGENCAPSID()

      !READ *, DUMMY
END SUBROUTINE

!-----------------------------------------------------------

SUBROUTINE TRANSFORMRIGIDTOC (CMIN, CMAX, XCOORDS, XRIGIDCOORDS)
      
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  
  INTEGER :: J1, J2, J5, J7, J9, NP        !NP = No of processors
  INTEGER :: CMIN, CMAX
  DOUBLE PRECISION :: P(3), RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
  DOUBLE PRECISION :: XRIGIDCOORDS(DEGFREEDOMS), XCOORDS(3*NATOMS)
  DOUBLE PRECISION :: COM(3) ! center of mass
  LOGICAL          :: GTEST !, ATOMTEST
  DOUBLE PRECISION :: MLATTICE(3,3)
  
  GTEST = .FALSE.
  NP = 1

! vr274 > are there additional lattice coordinates? If yes, setup transformation matrix
  IF(HAS_LATTICE_COORDS) THEN
    CALL GET_LATTICE_MATRIX(XRIGIDCOORDS(DEGFREEDOMS-5:DEGFREEDOMS), MLATTICE)
  ELSE ! vr274 > otherwise identity matrix
    MLATTICE = 0D0
    MLATTICE(1,1)=1d0
    MLATTICE(2,2)=1D0
    MLATTICE(3,3)=1D0
  ENDIF


  ! hk286 > coord transformations for rigid bodies CMIN to CMAX
  DO J1 = CMIN, CMAX
     J5   = 3*J1
     J7   = 3*NRIGIDBODY + J5
     P(:) = XRIGIDCOORDS(J7-2:J7)
     CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, GTEST)

! vr274 > MLATTICE can have lattice transformation or be identity matrix
     COM = matmul(MLATTICE, XRIGIDCOORDS(J5-2:J5))
     DO J2 = 1,  NSITEPERBODY(J1)
        J9 = RIGIDGROUPS(J2, J1)
        XCOORDS(3*J9-2:3*J9) = COM + MATMUL(RMI(:,:),SITESRIGIDBODY(J2,:,J1))
     ENDDO
     
  ENDDO
  
! hk286 > now the single atoms
! vr274 > this copies lattice coordinates as well which is stored in last 2 atoms
  IF (DEGFREEDOMS > 6 * NRIGIDBODY) THEN
     DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
        J9 = RIGIDSINGLES(J1)
        XCOORDS(3*J9-2:3*J9) = XRIGIDCOORDS(6*NRIGIDBODY + 3*J1-2:6*NRIGIDBODY + 3*J1)
     ENDDO
  ENDIF
      
END SUBROUTINE TRANSFORMRIGIDTOC

!----------------------------------------------------------

SUBROUTINE ROTATEINITIALREF ()
IMPLICIT NONE
DOUBLE PRECISION :: P(3)
INTEGER J1

! hk286 - rotate the system - new
  P(:) = OPTIMROTAVALUES(:)
!  P(1) = -1.0D0 * 8.0D0 * ATAN(1.0D0) 
!  P(2) = 0.0D0 !4.0D0 * ATAN(1.0D0) !-(8*ATAN(1.0D0) - 5.0D0)/DSQRT(2.0D0)
!  P(3) = 0.0D0 !4.0D0 * ATAN(1.0D0)
  DO J1 = 1, NRIGIDBODY
     CALL REDEFINERIGIDREF (J1,P)
  ENDDO

END SUBROUTINE ROTATEINITIALREF

!----------------------------------------------------------

SUBROUTINE REDEFINERIGIDREF (J1,P)

  IMPLICIT NONE
  
  INTEGER :: J1, J2     !No of processor
  DOUBLE PRECISION :: P(3), RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)

!  PRINT *, "REDEFINE ", J1
  CALL RMDRVT(P, RMI, DRMI1, DRMI2, DRMI3, .FALSE.)  
  DO J2 = 1, NSITEPERBODY(J1)
     SITESRIGIDBODY(J2,:,J1) = MATMUL(RMI(:,:),SITESRIGIDBODY(J2,:,J1))
  ENDDO

END SUBROUTINE REDEFINERIGIDREF

!----------------------------------------------------------

SUBROUTINE TRANSFORMCTORIGID (XCOORDS, XRIGIDCOORDS)
  USE COMMONS, ONLY: NATOMS
  USE VEC3
  USE ROTATIONS
  IMPLICIT NONE
  
  INTEGER :: J1, J2, J9     !No of processor
  DOUBLE PRECISION :: P(3)
  DOUBLE PRECISION :: COM(3), PNORM, PT(3,3), PI(3,3), MASS
  DOUBLE PRECISION :: XRIGIDCOORDS (DEGFREEDOMS), XCOORDS(3*NATOMS)

! vr274 > lattice matrix and inverse
  DOUBLE PRECISION MLATTICE(3,3), MLATTICEINV(3,3)
  INTEGER NLATTICECOORDS

! vr274 > if has lattice coordinates, setup matrices
  IF(HAS_LATTICE_COORDS) THEN
    NLATTICECOORDS=6
    CALL GET_LATTICE_MATRIX(XCOORDS(3*NATOMS-5:3*NATOMS),MLATTICE)
  ELSE
    NLATTICECOORDS=0
    MLATTICE=0
    MLATTICE(1,1)=1
    MLATTICE(2,2)=1
    MLATTICE(3,3)=1
  ENDIF
  CALL INVERT3X3(MLATTICE, MLATTICEINV)

! loop over all rigid bodies
  DO J1 = 1, NRIGIDBODY
     COM = 0.0D0
     MASS = 0.0D0
     ! calculate center of mass
     DO J2 = 1, NSITEPERBODY(J1)
        J9 = RIGIDGROUPS(J2, J1)
        COM = COM + XCOORDS(3*J9-2:3*J9)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        MASS = MASS + GR_WEIGHTS(RIGIDGROUPS(J2,J1))
     ENDDO
     COM = COM / MASS
     XRIGIDCOORDS(3*J1-2:3*J1) = COM

     ! determine 3 points + 3 reference points
     DO J2 = 1, 3
        J9 = RIGIDGROUPS(J2+REFVECTOR(J1)-1, J1)
        PT(:,J2) = XCOORDS(3*J9-2:3*J9) - COM
        PI(:,J2) = SITESRIGIDBODY(J2+REFVECTOR(J1)-1,:,J1)
! hk286 - seems to work better
        PNORM = DSQRT(DOT_PRODUCT(PT(:,J2),PT(:,J2)))
        PT(:,J2) = PT(:,J2) / PNORM
        PNORM = DSQRT(DOT_PRODUCT(PI(:,J2),PI(:,J2)))
        PI(:,J2) = PI(:,J2) / PNORM
! hk286 - seems to work better
     ENDDO

     ! determine orientation based on these points
     XRIGIDCOORDS(3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1) = ROT_GET_ORIENTATION_AA(PT, PI)
  ENDDO

! vr274> now translate everything to reduced units
  DO J1 = 1, NRIGIDBODY
    XRIGIDCOORDS(3*J1-2:3*J1) = MATMUL(MLATTICEINV, XRIGIDCOORDS(3*J1-2:3*J1))
  END DO

! hk286 > now the single atoms
  IF (DEGFREEDOMS > 6 * NRIGIDBODY) THEN
     DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY - NLATTICECOORDS)/3
        J9 = RIGIDSINGLES(J1)
        ! vr274 > added lattice stuff
        XRIGIDCOORDS(6*NRIGIDBODY + 3*J1-2:6*NRIGIDBODY + 3*J1) = MATMUL(MLATTICEINV, XCOORDS(3*J9-2:3*J9))
     ENDDO
  ENDIF

! vr274 > copy lattice coords
  IF(HAS_LATTICE_COORDS) THEN
    XRIGIDCOORDS(DEGFREEDOMS - 5:DEGFREEDOMS) =  XCOORDS(3*NATOMS-5:3*NATOMS)
  ENDIF
END SUBROUTINE TRANSFORMCTORIGID

!-----------------------------------------------------------

SUBROUTINE TRANSFORMCTORIGID_OLD (XCOORDS, XRIGIDCOORDS)

  USE COMMONS, ONLY: NATOMS
  USE VEC3
  IMPLICIT NONE

  INTEGER :: J1, J2, J9, NP     !No of processor
  DOUBLE PRECISION :: P(3), RMI(3,3), RMITEMP(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
  DOUBLE PRECISION :: TEMPP(3)
  DOUBLE PRECISION :: XMASS, YMASS, ZMASS, PNORM, PT(3,3), PI(3,3), MASS
  DOUBLE PRECISION :: QR1, QR2, QR3, PCROSS(3), PCROSS21(3), PCROSS22(3), P2R(3)
  DOUBLE PRECISION :: PTEMP(3), PTEMP1(3), PTEMP2(3), DUMMY, DUMMY2
  DOUBLE PRECISION :: XRIGIDCOORDS (DEGFREEDOMS), XCOORDS(3*NATOMS)
  DOUBLE PRECISION :: THETA, THETA2, PID, FRPISQ
  LOGICAL          :: GTEST, REDEFINET

! vr274 > lattice matrix and inverse
  DOUBLE PRECISION MLATTICE(3,3), MLATTICEINV(3,3)
  INTEGER NLATTICECOORDS

  PID        = 4.D0*ATAN(1.0D0)
  FRPISQ     = 4.D0*PID*PID
  GTEST = .FALSE.
  NP = 1

! vr274 > if has lattice coordinates, setup matrices
  IF(HAS_LATTICE_COORDS) THEN
    NLATTICECOORDS=6
    CALL GET_LATTICE_MATRIX(XCOORDS(3*NATOMS-5:3*NATOMS),MLATTICE)
  ELSE
    NLATTICECOORDS=0
    MLATTICE=0
    MLATTICE(1,1)=1
    MLATTICE(2,2)=1
    MLATTICE(3,3)=1
  ENDIF
  CALL INVERT3X3(MLATTICE, MLATTICEINV)

  DO J1 = 1, NRIGIDBODY

     REDEFINET = .FALSE.
10   IF (REDEFINET) THEN
        PRINT *, "REDEFINE CALLED"
!        TEMPP = (/0.2D0,0.2D0,0.2D0/)
        CALL REDEFINERIGIDREF(J1,OPTIMROTAVALUES)
        REDEFINET = .FALSE.
     ENDIF

! hk286 > compute the rotated reference vectors, save in PT(:,:)
     XMASS = 0.0D0
     YMASS = 0.0D0
     ZMASS = 0.0D0
     MASS = 0.0D0
     DO J2 = 1, NSITEPERBODY(J1)
        J9 = RIGIDGROUPS(J2, J1)
        IF( (J2-REFVECTOR(J1) .GE. 0) .AND. (J2-REFVECTOR(J1) .LE. 2) ) THEN
           PT(J2-REFVECTOR(J1)+1,:) = XCOORDS(3*J9-2:3*J9);
        ENDIF
        XMASS = XMASS + XCOORDS(3*J9-2)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        YMASS = YMASS + XCOORDS(3*J9-1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        ZMASS = ZMASS + XCOORDS(3*J9)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        MASS = MASS + GR_WEIGHTS(RIGIDGROUPS(J2,J1))
     ENDDO
     XMASS = XMASS / MASS
     YMASS = YMASS / MASS
     ZMASS = ZMASS / MASS
     XRIGIDCOORDS(3*J1-2) = XMASS
     XRIGIDCOORDS(3*J1-1) = YMASS
     XRIGIDCOORDS(3*J1)   = ZMASS         
     DO J2 = 1, 3
        PT(J2,1) = PT(J2,1) - XMASS
        PT(J2,2) = PT(J2,2) - YMASS
        PT(J2,3) = PT(J2,3) - ZMASS
     ENDDO
     DO J2 = 1, 3
        PNORM   = DSQRT(DOT_PRODUCT(PT(J2,:),PT(J2,:)))
        PT(J2,1) = PT(J2,1) / PNORM
        PT(J2,2) = PT(J2,2) / PNORM
        PT(J2,3) = PT(J2,3) / PNORM
        
! hk286 > the original, unrotated reference vectors PI
        PNORM   = DSQRT(DOT_PRODUCT(SITESRIGIDBODY(J2+REFVECTOR(J1)-1,:,J1),SITESRIGIDBODY(J2+REFVECTOR(J1)-1,:,J1)))
        PI(J2,1) = SITESRIGIDBODY(J2+REFVECTOR(J1)-1,1,J1) / PNORM
        PI(J2,2) = SITESRIGIDBODY(J2+REFVECTOR(J1)-1,2,J1) / PNORM
        PI(J2,3) = SITESRIGIDBODY(J2+REFVECTOR(J1)-1,3,J1) / PNORM           
     ENDDO

! hk286 > compute rotation around axis perpendicular to PI(1,:) and PT(1,:)
! hk286 > QR1 rotation around PCROSS as the axis 
     QR1 = DOT_PRODUCT(PI(1,:),PT(1,:))
     IF ( ((QR1 < 0.01D0 .AND. QR1 > -0.01D0) .OR. &
         (QR1 < 2.0D0*PID + 0.01D0 .AND. QR1 > 2.0D0 * PID -0.01D0 )) &
         .AND. (RIGIDOPTIMROTAT .EQV. .TRUE.) ) THEN
        REDEFINET = .TRUE.
        GOTO 10
     ENDIF
     IF (QR1 > 1.0D0) QR1 = 1.0D0
     IF (QR1 < -1.0D0) QR1 = -1.0D0
     QR1 = ACOS(QR1)
     PCROSS(1) = PI(1,2)*PT(1,3) - PT(1,2)*PI(1,3)
     PCROSS(2) = PI(1,3)*PT(1,1) - PT(1,3)*PI(1,1)
     PCROSS(3) = PI(1,1)*PT(1,2) - PT(1,1)*PI(1,2)
     PNORM = DSQRT(DOT_PRODUCT(PCROSS,PCROSS))
     IF ((QR1 < 1D-6 .AND. QR1 > -1D-6).OR.(PNORM < 1D-6)) THEN
        QR1 = QR1 + 4.0D0*ASIN(1.0D0)
        PCROSS(:) = PT(1,:) 
     ELSE
        PCROSS(1) = PCROSS(1) / PNORM
        PCROSS(2) = PCROSS(2) / PNORM
        PCROSS(3) = PCROSS(3) / PNORM
     ENDIF
     CALL RMDRVT(QR1*PCROSS, RMI, DRMI1, DRMI2, DRMI3, GTEST)

! hk286 > any additional rotation not accounted by the above operation must be
! hk286 > rotation around PT(1,:), now compute this rotation
! hk286 > rotate the second reference vector PI(2,:) around PCROSS by QR1
     P2R = MATMUL(RMI(:,:),PI(2,:))
! hk286 > take suitable projections
! hk286 > QR2 is the amount of rotation around PT(1,:)
     PCROSS21(1) = P2R(2)*PT(1,3) - PT(1,2)*P2R(3)
     PCROSS21(2) = P2R(3)*PT(1,1) - PT(1,3)*P2R(1)
     PCROSS21(3) = P2R(1)*PT(1,2) - PT(1,1)*P2R(2)
     PCROSS22(1) = PT(2,2)*PT(1,3) - PT(1,2)*PT(2,3)
     PCROSS22(2) = PT(2,3)*PT(1,1) - PT(1,3)*PT(2,1)
     PCROSS22(3) = PT(2,1)*PT(1,2) - PT(1,1)*PT(2,2)

     PNORM = DSQRT(DOT_PRODUCT(PCROSS21,PCROSS21))
     PCROSS21 = PCROSS21 / PNORM
     PNORM = DSQRT(DOT_PRODUCT(PCROSS22,PCROSS22))
     PCROSS22 = PCROSS22 / PNORM
     QR2 = DOT_PRODUCT(PCROSS21,PCROSS22)
     IF ( ((QR2 < 0.01D0 .AND. QR2 > -0.01D0) .OR. &
        (QR2 < 2.0D0*PID + 0.01D0 .AND. QR2 > 2.0D0 * PID -0.01D0)) .AND. &
        (RIGIDOPTIMROTAT .EQV. .TRUE.) ) THEN
        REDEFINET = .TRUE.
        GOTO 10
     ENDIF     
     IF (QR2 > 1.0D0) QR2 = 1.0D0
     IF (QR2 < -1.0D0) QR2 = -1.0D0
     QR2 = ACOS(QR2)
     CALL RMDRVT(QR2*PT(1,:), RMI, DRMI1, DRMI2, DRMI3, GTEST)
     PTEMP = MATMUL(RMI(:,:),P2R(:))
     DUMMY = (PTEMP(1)-PT(2,1))**2 + (PTEMP(2)-PT(2,2))**2 + (PTEMP(3)-PT(2,3))**2
     CALL RMDRVT((2.0D0*PID-QR2)*PT(1,:), RMI, DRMI1, DRMI2, DRMI3, GTEST)
     PTEMP = MATMUL(RMI(:,:),P2R(:))
     DUMMY2 = (PTEMP(1)-PT(2,1))**2 + (PTEMP(2)-PT(2,2))**2 + (PTEMP(3)-PT(2,3))**2
     IF (DUMMY2 < DUMMY) THEN
        QR2 = 2.0D0 * PID - QR2
     ENDIF

! hk286 > Construct quarternions
!         Quarternion1 = [COS(QR1/2) PCROSS(:)*SIN(QR1/2)]
!         Quarternion2 = [COS(QR2/2) PT(1,:)*SIN(QR2/2)]
     PTEMP1 = PCROSS(:)*SIN(QR1/2)
     PTEMP2 = PT(1,:)*SIN(QR2/2)
     QR3 = COS(QR1/2)*COS(QR2/2) - DOT_PRODUCT(PTEMP1,PTEMP2)
     IF (QR3 > 1.0D0) QR3 = 1.0D0
     IF (QR3 < -1.0D0) QR3 = -1.0D0    
     QR3 = 2*ACOS(QR3)
     IF ( ((QR3 < 0.01D0 .AND. QR3 > -0.01D0) .OR. &
        (QR3 < 2.0D0*PID + 0.01D0 .AND. QR3 > 2.0D0 * PID -0.01D0)) .AND. &
        (RIGIDOPTIMROTAT .EQV. .TRUE.) ) THEN
        REDEFINET = .TRUE.
        GOTO 10
     ENDIF          
     P(1) = (PT(1,2)*PCROSS(3) - PCROSS(2)*PT(1,3))*SIN(QR1/2)*SIN(QR2/2)
     P(2) = (PT(1,3)*PCROSS(1) - PCROSS(3)*PT(1,1))*SIN(QR1/2)*SIN(QR2/2)
     P(3) = (PT(1,1)*PCROSS(2) - PCROSS(1)*PT(1,2))*SIN(QR1/2)*SIN(QR2/2)
     P(:) = P(:) + COS(QR1/2)*PT(1,:)*SIN(QR2/2) + COS(QR2/2)*PCROSS(:)*SIN(QR1/2)
     IF (abs(sin(QR3)) < 1D-6 ) THEN
        PNORM = SQRT(DOT_PRODUCT(P,P))
        P(:) = P(:) * 2.0D0
     ELSE
        P(:) = P(:) / SIN(QR3/2) * QR3
     ENDIF
     XRIGIDCOORDS(3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1) = P(:)

  ENDDO
  
! vr274> now translate everything to reduced units
  DO J1 = 1, NRIGIDBODY
    XRIGIDCOORDS(3*J1-2:3*J1) = MATMUL(MLATTICEINV, XRIGIDCOORDS(3*J1-2:3*J1))
  END DO

! hk286 > now the single atoms
  IF (DEGFREEDOMS > 6 * NRIGIDBODY) THEN
     DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY - NLATTICECOORDS)/3
        J9 = RIGIDSINGLES(J1)
        ! vr274 > added lattice stuff
        XRIGIDCOORDS(6*NRIGIDBODY + 3*J1-2:6*NRIGIDBODY + 3*J1) = MATMUL(MLATTICEINV, XCOORDS(3*J9-2:3*J9))
     ENDDO
  ENDIF

! vr274 > copy lattice coords
  IF(HAS_LATTICE_COORDS) THEN
    XRIGIDCOORDS(DEGFREEDOMS - 5:DEGFREEDOMS) =  XCOORDS(3*NATOMS-5:3*NATOMS)
  ENDIF


! hk286 > reduce the amount of rotation between 2*PI and 4*PI
  DO J1 = 1, NRIGIDBODY

     J2       = 3*NRIGIDBODY + 3*J1
     THETA2  = DOT_PRODUCT(XRIGIDCOORDS(J2-2:J2),XRIGIDCOORDS(J2-2:J2))
     THETA   = DSQRT(THETA2)
     IF (THETA2 > FRPISQ) THEN
        THETA2   = DSQRT(THETA2)
        THETA    = THETA2 - INT(THETA2/(2.D0*PID))*2.D0*PID
        XRIGIDCOORDS(J2-2:J2) = XRIGIDCOORDS(J2-2:J2)/THETA2 * THETA
     ENDIF
!     IF (THETA > PID) THEN
!        XRIGIDCOORDS(J2-2:J2) = XRIGIDCOORDS(J2-2:J2)/THETA * (-2.0D0*PID + THETA)
!        THETA = 2.0D0 * PID - THETA
!     ENDIF
     IF (RIGIDOPTIMROTAT .EQV. .TRUE.) THEN
        XRIGIDCOORDS(J2-2:J2) = XRIGIDCOORDS(J2-2:J2)/THETA * (2.0D0*PID + THETA)
        THETA = 2.0D0*PID + THETA
     ENDIF

  ENDDO
  
END SUBROUTINE TRANSFORMCTORIGID_OLD

!-----------------------------------------------------------


SUBROUTINE TRANSFORMGRAD (G, XR, GR)
  
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  
  INTEGER          :: J1, J2, J9
  DOUBLE PRECISION :: G(3*NATOMS), XR(DEGFREEDOMS), GR(DEGFREEDOMS)
  DOUBLE PRECISION :: PI(3)
  DOUBLE PRECISION :: DR1(3),DR2(3),DR3(3) 
  DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
  LOGICAL :: GTEST
  INTEGER :: NLATTICECOORDS
  DOUBLE PRECISION :: MLATTICE(3,3)
  
  NLATTICECOORDS=0
  IF(HAS_LATTICE_COORDS) THEN
      NLATTICECOORDS=6
  ENDIF

  GTEST = .TRUE.
  GR(:) = 0.0D0
  
  DO J1 = 1, NRIGIDBODY
     
     PI = XR(3*NRIGIDBODY+3*J1-2 : 3*NRIGIDBODY+3*J1)
     CALL RMDRVT(PI, RMI, DRMI1, DRMI2, DRMI3, GTEST)
     
     DO J2 = 1, NSITEPERBODY(J1)
        J9 = RIGIDGROUPS(J2, J1)

! hk286 > translation
        GR(3*J1-2:3*J1) = GR(3*J1-2:3*J1) + G(3*J9-2:3*J9)
        
! hk286 > rotation
        DR1(:) = MATMUL(DRMI1,SITESRIGIDBODY(J2,:,J1))
        DR2(:) = MATMUL(DRMI2,SITESRIGIDBODY(J2,:,J1))
        DR3(:) = MATMUL(DRMI3,SITESRIGIDBODY(J2,:,J1))
        GR(3*NRIGIDBODY+3*J1-2) = GR(3*NRIGIDBODY+3*J1-2) + DOT_PRODUCT(G(3*J9-2:3*J9),DR1(:))
        GR(3*NRIGIDBODY+3*J1-1) = GR(3*NRIGIDBODY+3*J1-1) + DOT_PRODUCT(G(3*J9-2:3*J9),DR2(:))
        GR(3*NRIGIDBODY+3*J1)   = GR(3*NRIGIDBODY+3*J1)   + DOT_PRODUCT(G(3*J9-2:3*J9),DR3(:))
     ENDDO
  ENDDO

! hk286 - testing 6/6/12
  IF (FREEZERIGIDBODYT .EQV. .TRUE.) THEN
     GR(3*NRIGIDBODY-2:3*NRIGIDBODY) = 0.0D0
     GR(6*NRIGIDBODY-2:6*NRIGIDBODY) = 0.0D0
  ENDIF

! hk286 > single atoms
! vr274 > and lattice
  IF (DEGFREEDOMS > 6 * NRIGIDBODY - NLATTICECOORDS) THEN
     DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY - NLATTICECOORDS)/3
        J9 = RIGIDSINGLES(J1)
        GR(6*NRIGIDBODY + 3*J1-2:6*NRIGIDBODY + 3*J1) = G(3*J9-2:3*J9)
     ENDDO
  ENDIF

  IF(HAS_LATTICE_COORDS) THEN
      CALL GET_LATTICE_MATRIX(XR(DEGFREEDOMS-5:DEGFREEDOMS),MLATTICE)

      ! vr274> for lattice, go to reduced coordinates
      DO J1 = 1, NRIGIDBODY
          GR(3*J1-2:3*J1) =  matmul(transpose(mlattice), GR(3*J1-2:3*J1))
      ENDDO
      ! vr274> and single atoms
      IF (DEGFREEDOMS > 6 * NRIGIDBODY + NLATTICECOORDS) THEN
          DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY - NLATTICECOORDS)/3
              J2 = 6*NRIGIDBODY + 3*J1
              GR(J2-2:J2) = matmul(transpose(mlattice), GR(J2-2:J2))
          ENDDO
      ENDIF
      ! copy lattice gradient
      GR(DEGFREEDOMS-5:DEGFREEDOMS) = G(3*NATOMS-5:3*NATOMS)
  ENDIF

END SUBROUTINE TRANSFORMGRAD

!--------------------------------------------------------------

! hk286 > Often we want to check if the atoms grouped in a rigid body has moved or not
! hk286 > They should not if everything is done correctly
! hk286 > REDEFINESITEST = .FALSE. then it prints to standard output
! hk286 > REDEFINESITEST = .TRUE. then regroup atoms, SITESRIGIDBODY rewritten

SUBROUTINE CHECKSITES (REDEFINESITEST, COORDS)
      
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE

  INTEGER :: J1, J2, J3, DUMMY
  DOUBLE PRECISION :: XMASS, YMASS, ZMASS, PNORM, MASS
  DOUBLE PRECISION :: XSITESRIGIDBODY(MAXSITE,3,NRIGIDBODY)
  DOUBLE PRECISION :: COORDS(3*NATOMS)
  LOGICAL :: RTEST, REDEFINESITEST
  

  DO J1 = 1, NRIGIDBODY
     DO J2 = 1, NSITEPERBODY(J1)
        DUMMY = RIGIDGROUPS(J2,J1)
        XSITESRIGIDBODY(J2,:,J1) = COORDS(3*DUMMY-2:3*DUMMY)
     ENDDO
  ENDDO

  DO J1 = 1, NRIGIDBODY
     XMASS = 0.0D0
     YMASS = 0.0D0
     ZMASS = 0.0D0
     MASS = 0.0D0
     DO J2 = 1, NSITEPERBODY(J1)
        XMASS = XMASS + XSITESRIGIDBODY(J2,1,J1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        YMASS = YMASS + XSITESRIGIDBODY(J2,2,J1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        ZMASS = ZMASS + XSITESRIGIDBODY(J2,3,J1)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        MASS = MASS + GR_WEIGHTS(RIGIDGROUPS(J2,J1))
     ENDDO
     XMASS = XMASS / MASS
     YMASS = YMASS / MASS
     ZMASS = ZMASS / MASS
     DO J2 = 1, NSITEPERBODY(J1)
        XSITESRIGIDBODY(J2,1,J1) = XSITESRIGIDBODY(J2,1,J1) - XMASS
        XSITESRIGIDBODY(J2,2,J1) = XSITESRIGIDBODY(J2,2,J1) - YMASS
        XSITESRIGIDBODY(J2,3,J1) = XSITESRIGIDBODY(J2,3,J1) - ZMASS
     ENDDO
  ENDDO
  

  IF (REDEFINESITEST) THEN
!     PRINT *, " SITES REDEFINED "
     SITESRIGIDBODY(:,:,:) = XSITESRIGIDBODY(:,:,:)

!Checks: (1) Atoms 1 and 2 do not sit on COM, and (2) Vector 1 and 2 are not parallel
  
     DO J1 = 1, NRIGIDBODY
        REFVECTOR(J1) = 1
        RTEST = .TRUE.
        DO WHILE (RTEST)
           RTEST = .FALSE.
           DO J2 = REFVECTOR(J1), REFVECTOR(J1) + 1 
              PNORM = SQRT(DOT_PRODUCT(SITESRIGIDBODY(J2,:,J1),SITESRIGIDBODY(J2,:,J1)))
              IF ( (PNORM  < 0.001) .AND. (PNORM > -0.001)) THEN
                 RTEST = .TRUE.
              ENDIF
           ENDDO
           PNORM = DOT_PRODUCT(SITESRIGIDBODY(REFVECTOR(J1),:,J1),SITESRIGIDBODY(REFVECTOR(J1)+1,:,J1)) 
           PNORM = PNORM / SQRT(DOT_PRODUCT(SITESRIGIDBODY(REFVECTOR(J1),:,J1),SITESRIGIDBODY(REFVECTOR(J1),:,J1))) 
           PNORM = PNORM / SQRT(DOT_PRODUCT(SITESRIGIDBODY(REFVECTOR(J1)+1,:,J1),SITESRIGIDBODY(REFVECTOR(J1)+1,:,J1)))
           IF (PNORM < 0.0) PNORM = -1.0D0 * PNORM
           IF ( (PNORM < 1.0 + 0.001) .AND. (PNORM > 1.0 - 0.001) ) THEN
              RTEST = .TRUE.
           ENDIF
           IF (RTEST) THEN
              REFVECTOR(J1) = REFVECTOR(J1) + 1               
           ENDIF
        ENDDO
     ENDDO
  ELSE
!     PRINT *, XSITESRIGIDBODY
  ENDIF

END SUBROUTINE CHECKSITES

!--------------------------------------------------------------

! vr274 > build the lattice matrix.
!         The matrix is an triangular matrix,  c vector is always perpendicular to z
SUBROUTINE GET_LATTICE_MATRIX(LATTICE_COORDS, M)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN) :: LATTICE_COORDS(6)
    DOUBLE PRECISION, INTENT(OUT) :: M(3,3)
    M=0
    M(1,1) = LATTICE_COORDS(1)
    M(2,1) = LATTICE_COORDS(2)
    M(3,1) = LATTICE_COORDS(3)
    M(2,2) = LATTICE_COORDS(4)
    M(3,2) = LATTICE_COORDS(5)
    M(3,3) = LATTICE_COORDS(6)
END SUBROUTINE

! vr274 > set lattice coordinates from lattice matrix.
!         The matrix is an triangular matrix,  c vector is always perpendicular to z
SUBROUTINE SET_LATTICE_MATRIX(LATTICE_COORDS, M)
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(OUT) :: LATTICE_COORDS(6)
    DOUBLE PRECISION, INTENT(IN) :: M(3,3)

    LATTICE_COORDS(1) = M(1,1)
    LATTICE_COORDS(2) = M(2,1)
    LATTICE_COORDS(3) = M(3,1)
    LATTICE_COORDS(4) = M(2,2)
    LATTICE_COORDS(5) = M(3,2)
    LATTICE_COORDS(6) = M(3,3)
END SUBROUTINE

END MODULE GENRIGID
