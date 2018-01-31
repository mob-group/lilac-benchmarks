MODULE GENRIGID

      INTEGER :: NRIGIDBODY, DEGFREEDOMS, MAXSITE, NRELAXRIGIDR, NRELAXRIGIDA
      INTEGER :: XNATOMS
      INTEGER, ALLOCATABLE :: NSITEPERBODY(:), REFVECTOR(:), RIGIDSINGLES(:)
      INTEGER, ALLOCATABLE, DIMENSION (:,:) :: RIGIDGROUPS
      DOUBLE PRECISION, ALLOCATABLE :: RIGIDCOORDS(:)
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:,:) :: SITESRIGIDBODY
      DOUBLE PRECISION, ALLOCATABLE :: GR_WEIGHTS(:) ! weights for com calculation, e.g. masses

! Sm, S, RBMASS and COG are lists of: the mass-weighted and non-mass-weighted tensors
! of gyration, the total mass and the centre of geometry, for all rigid bodies in the system.
! RCOM is a list of the centre of mass for each rigid body in the reference structure.
      DOUBLE PRECISION, ALLOCATABLE :: S(:,:), Sm(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: RBMASS(:)
      DOUBLE PRECISION, ALLOCATABLE :: COG(:,:), RCOM(:,:)

      LOGICAL :: MASSES ! Used in setup to indicate a user-specified mass file.

      LOGICAL :: RIGIDINIT         ! The keyword used to activate genrigid-specific blocks of code in OPTIM and GMIN
      LOGICAL :: ATOMRIGIDCOORDT   ! TRUE when we are working in atomistic coordinates, FALSE otherwise.
      LOGICAL :: RELAXRIGIDT, GENRIGIDT
      LOGICAL :: RIGIDMOLECULEST   ! Specifies a system of fully rigid molecules

      DOUBLE PRECISION, ALLOCATABLE :: IINVERSE(:,:,:)

      LOGICAL :: RIGIDOPTIMROTAT, FREEZERIGIDBODYT
      DOUBLE PRECISION :: OPTIMROTAVALUES(3)
      LOGICAL :: AACONVERGENCET
      LOGICAL :: DUMMYLOGICAL = .FALSE.
      INTEGER, ALLOCATABLE :: LRBNPERMGROUP(:), LRBNPERMSIZE(:,:), LRBPERMGROUP(:,:), LRBNSETS(:,:), LRBSETS(:,:,:)

!   vr274:  added lattice coordinates
!           if HAS_LATTICE_COORDS is true, the last two atoms are treated
!           as lattice coordintes and rigidcoords is in reduced lattice units
      LOGICAL HAS_LATTICE_COORDS
!   csw34>  RIGIDISRIGID = logical array, size NATOMS, TRUE if atom is part of
!           RB
      LOGICAL, ALLOCATABLE :: RIGIDISRIGID(:)
      INTEGER, ALLOCATABLE :: RB_BY_ATOM(:) ! sn402: records which RB an atom belongs to

!-----------------------------------------------------------------------------------!
! NRIGIDBODY  = number of rigid bodies
! DEGFREEDOMS = number of degrees of freedom = 6 * NRIGIDBODY + 3 * ADDITIONAL ATOMS
! MAXSITE     = maximum number of sites in a rigid body
! NRELAXRIGIDR = rigid body minimisation for this number of steps
! NRELAXRIGIDA = atom minimisation for this number of steps
! NSITEPERBODY= number of rigid body sites, no need to be the same for all bodies
! REFVECTOR   = reference vector for the atomistic to rigic coordinate transformation (used for checking only)
! RIGIDSINGLES= list of atoms not in rigid bodies
! RIGIDGROUPS = list of atoms in rigid bodies, need a file called rbodyconfig
! RIGIDCOORDS = 6 * NRIGIDBODY + 3 * ADDITIONAL ATOMS coordinates
! SITESRIGIDBODY = coordinates of the rigid body sites
! RIGIDINIT   = logical variable for generalised rigid body
! ATOMRIGIDCOORDT, .TRUE. = atom coords active, .FALSE. = rigid coords active, used in mylbfgs & potential
! GENRIGIDT = generalised rigid body takestep taken if .TRUE. (GMIN only)
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
  USE COMMONS, only: NATOMS, ATMASS
! hk286
  USE MODAMBER9, ONLY : ATMASS1
  USE AMBER12_INTERFACE_MOD, ONLY: AMBER12_MASSES
  USE KEY, ONLY: AMBER12T

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
  ALLOCATE (IINVERSE(NRIGIDBODY,3,3))
! csw34> 
  ALLOCATE (RIGIDISRIGID(NATOMS))
  ALLOCATE (RB_BY_ATOM(NATOMS))
  RIGIDISRIGID=.FALSE.
  RB_BY_ATOM(:) = -1

  ! by default use center of geometry - check if ok for AMBER 9/8/12
  IF ( ALLOCATED(ATMASS1) ) THEN
     GR_WEIGHTS=ATMASS1
     MASSES=.TRUE.  ! We need to have this set to deactivate the safety check in RB_DISTANCE
                    ! (which doesn't apply if the atoms have non-uniform mass)
  ELSE IF ( ALLOCATED(AMBER12_MASSES) .AND. AMBER12T ) THEN
     ! AMBER 12 atom masses
     GR_WEIGHTS = AMBER12_MASSES
     MASSES=.TRUE. ! We need to have this set to deactivate the safety check in RB_DISTANCE
                   ! (which doesn't apply if the atoms have non-uniform mass)
  ELSE IF (MASSES) THEN
     GR_WEIGHTS=ATMASS
  ELSE
     GR_WEIGHTS=1.0D0
  ENDIF
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

  INTEGER          :: INFO
  INTEGER, PARAMETER :: LWORK = 1000000 ! the dimension is set arbitrarily
  INTEGER :: I, J
  DOUBLE PRECISION :: DR(3), KBLOCK(3,3), KBEGNV(3)
  DOUBLE PRECISION :: WORK(LWORK)

  DOUBLE PRECISION :: DET

  
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
  ! Here, DEGFREEDOMS is just being used as a dummy to count the number of atoms which do belong to rigid bodies
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
        IF (DUMMY.EQ.DEGFREEDOMS/3-2*NRIGIDBODY) THEN
           WRITE(*,*) "genrigid> Error. More free atoms than expected."
           WRITE(*,*) "Likely problem with rbodyconfig."
           STOP
        ENDIF
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

  DO J1 = 1, NRIGIDBODY
     IINVERSE(J1,:,:) = 0.0D0
  END DO

  IF (AACONVERGENCET .EQV. .TRUE.) THEN
     DO J1 = 1, NRIGIDBODY
        KBLOCK(:,:) = 0.0D0
        DO J2 = 1, NSITEPERBODY(J1)
           DR(:)  = SITESRIGIDBODY(J2,:,J1)
           DO I = 1, 3
              ! KBLOCK is the unweighted tensor of gyration
              KBLOCK(I,I) = KBLOCK(I,I) + (DR(1)*DR(1) + DR(2)*DR(2) + DR(3)*DR(3))
              DO J = 1, 3    ! could have been J = 1, I; KBLOCK is a symmetric matrix
                 KBLOCK(I,J) = KBLOCK(I,J) - DR(I)*DR(J)
              ENDDO
           ENDDO
        ENDDO
        CALL DSYEV('V','U',3,KBLOCK,3,KBEGNV,WORK,LWORK,INFO)
        CALL RBDET(KBLOCK, DET)
        IF (DET < 0.0D0) THEN
           KBLOCK(:,3) = -KBLOCK(:,3)
           CALL RBDET(KBLOCK, DET)
           IF (DET < 0.0D0) THEN
              PRINT *, "GENRIGID> BAD ALIGNMENT", J1
              STOP
           ENDIF
        ENDIF
        KBLOCK = TRANSPOSE(KBLOCK)
!        PRINT *, KBEGNV
        IINVERSE(J1,1,1) = 1.0D0/KBEGNV(1)
        IINVERSE(J1,2,2) = 1.0D0/KBEGNV(2)
        IINVERSE(J1,3,3) = 1.0D0/KBEGNV(3)
        DO J2 = 1, NSITEPERBODY(J1)
           SITESRIGIDBODY(J2,:,J1) = MATMUL(KBLOCK,SITESRIGIDBODY(J2,:,J1))
        ENDDO
     ENDDO
  ENDIF
!  PRINT *, SITESRIGIDBODY(1,:,1)
!  PRINT *, SITESRIGIDBODY(2,:,1)

! hk286 > make sure the two atoms used as reference for rigid bodies are suitable
! Checks: (1) Atoms 1 and 2 do not sit on COM, and (2) Vector 1 and 2 are not parallel
! sn402: this seems like a very restrictive test - for instance, it means we can't work
! with fixed-length diatomics using this framework. Is there any way of avoiding this?
  
  DO J1 = 1, NRIGIDBODY
     REFVECTOR(J1) = 1
     RTEST = .TRUE.
     DO WHILE (RTEST)
        RTEST = .FALSE.
        DO J2 = REFVECTOR(J1), REFVECTOR(J1) + 1 
            IF (J2.GT.NSITEPERBODY(J1)) THEN
                WRITE(*,*) "genrigid.f90> There is no pair of atoms which makes an appropriate rigid body reference structure"
                STOP
            ENDIF
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

! Compute the tensors of gyration and centres of mass for all bodies.
  CALL SETUP_TENSORS()

END SUBROUTINE

!-----------------------------------------------------------
SUBROUTINE GENRIGID_READ_FROM_FILE ()
      
  USE COMMONS, ONLY: NATOMS, ATMASS, DEBUG
  USE KEY, ONLY: AMBER12T
  USE MODAMBER9, ONLY: ATMASS1
  IMPLICIT NONE

  CHARACTER(LEN=10) CHECK1, CHECK2
  INTEGER :: J1, J2, DUMMY, iostatus
  DOUBLE PRECISION :: INICOORDS(3*NATOMS)

! hk286 > read atomistic coordinates
! hk286 > in future, no need for separate coordsinirigid
! hk286 > currently the input coords files vary for CHARMM, AMBER, and RIGID BODIES
  IF (NATOMS == 0) THEN
     PRINT *, "ERROR STOP NOW > During generalised rigid body initialisation NATOMS = 0"
     STOP
  ENDIF

! vr274 > by standard don't do lattice coordinates
  HAS_LATTICE_COORDS = .FALSE.
  OPEN(UNIT = 28, FILE = 'coordsinirigid', STATUS = 'OLD')
  DO J1 = 1, NATOMS
     READ(28, *) INICOORDS(3*J1-2), INICOORDS(3*J1-1), INICOORDS(3*J1)
  ENDDO
  CLOSE(UNIT = 28)
! hk286 > determine no of rigid bodies
  NRIGIDBODY=0
  OPEN(UNIT=222,FILE='rbodyconfig',status='old')
  DO
     READ(222,*,IOSTAT=iostatus) CHECK1
     CALL UPPERCASE(CHECK1)
     IF (iostatus<0) THEN
        CLOSE(222)
        EXIT
     ELSE IF (TRIM(ADJUSTL(CHECK1)).EQ.'GROUP') then
        NRIGIDBODY=NRIGIDBODY+1
     ENDIF
  END DO
  CLOSE(222)
! hk286 > determine maximum no of rigid body sites
  MAXSITE = 0
  OPEN(UNIT=222,FILE='rbodyconfig',status='old')
  DO J1 = 1, NRIGIDBODY
     READ(222,*) CHECK1,DUMMY
     IF (MAXSITE < DUMMY) MAXSITE = DUMMY
     DO J2 = 1, DUMMY
        READ(222,*) CHECK1
     ENDDO
  ENDDO
  CLOSE(222)
  ! Insert mass read from file here.
  INQUIRE(FILE='mass', EXIST=MASSES)
  IF(MASSES) THEN
    IF(DEBUG) WRITE(*,*) "genrigid> Reading masses from file 'mass'"
    OPEN(UNIT=222,FILE='mass',status='old')
    IF(.NOT. ALLOCATED(ATMASS)) THEN
        ALLOCATE(ATMASS(NATOMS))
    ELSE IF (SIZE(ATMASS) .NE. NATOMS) THEN
        WRITE(*,*) "genrigid> problem: ATMASS has been allocated with the wrong number of atoms: ", &
            & SIZE(ATMASS)
        STOP
    ENDIF
    DO J1 = 1,NATOMS
        READ(222,*) CHECK1, ATMASS(J1)
    ENDDO
    CLOSE(222)
  ELSE
    IF(DEBUG .AND. .NOT.(ALLOCATED(ATMASS1).OR.AMBER12T)) WRITE(*,*) "genrigid> No file 'mass': assuming all atoms have mass = 1."
  ENDIF

!  vr274> Calling function for allocation to make more general (setup rigid bodies from code)
  CALL GENRIGID_ALLOCATE(NRIGIDBODY,MAXSITE)
! hk286 > initialise SITESRIGIDBODY, RIGIDGROUPS, RIGIDSINGLES 
  OPEN(UNIT=222,FILE='rbodyconfig',status='unknown')
  DO J1 = 1, NRIGIDBODY
     READ(222,*) CHECK1, NSITEPERBODY(J1)
     DO J2 = 1, NSITEPERBODY(J1)
        READ(222,*) RIGIDGROUPS(J2,J1)
! csw34> check to make sure the current atom is not already in a rigid body
        IF (RIGIDISRIGID(RIGIDGROUPS(J2,J1))) THEN
            PRINT *," genrigid> ERROR: atom ",RIGIDGROUPS(J2,J1)," is in multiple rigid bodies! Stopping."
            STOP
        ELSE       
! csw34> if not, flag the current atom
            RIGIDISRIGID(RIGIDGROUPS(J2,J1))=.TRUE.
            RB_BY_ATOM(RIGIDGROUPS(J2,J1)) = J1
        ENDIF
! vr274> Moved initialization of coordinates to GENRIGID_INITIALISE, here only read the setup
!        SITESRIGIDBODY(J2,:,J1) = COORDS(3*DUMMY-2:3*DUMMY,1)
     ENDDO
  ENDDO
  CLOSE(222)
  CALL GENRIGID_INITIALISE(INICOORDS)
END SUBROUTINE GENRIGID_READ_FROM_FILE

!-----------------------------------------------------------
SUBROUTINE SETUP_TENSORS()
! Some information regarding the rigid bodies does not change through a simulation -
! e.g. the tensor of gyration. So these can all be computed once for the reference
! orientations of each rigid body, then re-used later.
USE COMMONS, ONLY: DEBUG
IMPLICIT NONE
INTEGER :: J1, J2, J5, J6, J7, ATOMINDEX, DUMMY

! Sm and S are lists of the mass-weighted and non-mass-weighted tensors of gyration for
! the rigid bodies. These consist of NRIGIDBODY blocks of 3x3 matrices, one for each body.
! RBMASS is simply a list of the total masses for each rigid body
! COG is a list of 3-tuples holding the reduced coordinates of the centre of geometry for
! each rigid body. These are obtained from coordsinirigid, and the molecule-frame coordinate
! system is usually chosen such that the origin is on the centre of geometry for every
! rigid body. If this is the case, COG(:,:) = 0 (approximately).
! Note that if GR_WEIGHTS are not all 1 then COG no longer represents the unweighted
! centre of geometry. If GR_WEIGHTS are set equal to the masses of the atoms, for instance,
! then COG becomes a list of the centres of mass.
ALLOCATE(S(3*NRIGIDBODY,3), Sm(3*NRIGIDBODY,3),RBMASS(NRIGIDBODY),COG(NRIGIDBODY,3),RCOM(NRIGIDBODY,3))
S(:,:) = 0.0D0
Sm(:,:) = 0.0D0
RBMASS(:) = 0.0D0
COG(:,:) = 0.0D0
RCOM(:,:) = 0.0D0

! All quantities are calculated according to the definitions in Ruehle, Kusumaatmaja,
! Chakrabarti and Wales, JCTC (2013)

DUMMY = 0
DO J1 = 1, NRIGIDBODY
    DO J6 = 1,3
        DO J7 = 1,3
            DO J5 = 1, NSITEPERBODY(J1)
                ATOMINDEX = RIGIDGROUPS(J5,J1)
                S(J6+DUMMY,J7) = S(J6+DUMMY,J7) + SITESRIGIDBODY(J5,J6,J1)*SITESRIGIDBODY(J5,J7,J1)
                Sm(J6+DUMMY,J7) = Sm(J6+DUMMY,J7) + GR_WEIGHTS(ATOMINDEX)*SITESRIGIDBODY(J5,J6,J1)*SITESRIGIDBODY(J5,J7,J1)
            ENDDO
        ENDDO
    ENDDO
    DUMMY = DUMMY + 3

    ! Compute total mass and centre of geometry of the rigid body. This mostly exists to check that the
    ! reference-frame coordinates have been correctly centred.
    DO J2 = 1, NSITEPERBODY(J1)
        RBMASS(J1) = RBMASS(J1) + GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        COG(J1,:) = COG(J1,:) + SITESRIGIDBODY(J2,:,J1) ! sn402: problem here!
        RCOM(J1,:) = RCOM(J1,:) + GR_WEIGHTS(RIGIDGROUPS(J2,J1))*SITESRIGIDBODY(J2,:,J1)
    ENDDO
    RCOM(J1,:) = RCOM(J1,:)/RBMASS(J1)      ! centre of mass
    COG(J1,:) = COG(J1,:)/NSITEPERBODY(J1) ! centre of geometry

ENDDO

RETURN
END SUBROUTINE SETUP_TENSORS
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
  USE COMMONS, ONLY: NATOMS, PARAM1,PARAM2,PARAM3, DEBUG ! hk286
  USE KEY, ONLY : PERMDIST, NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, SETS, BULKT, BULK_BOXVEC, AMBER12T, AMBERT, NABT
  USE VEC3
  USE ROTATIONS
  IMPLICIT NONE
  
  INTEGER :: J1, J2, J9     !No of processor
  DOUBLE PRECISION :: P(3), XSAVE(3)
  DOUBLE PRECISION :: COM(3), MASS
  DOUBLE PRECISION :: XRIGIDCOORDS (DEGFREEDOMS), XCOORDS(3*NATOMS), TEMPXCOORDS(3*NATOMS), SAVEXCOORDS(3*NATOMS)

! vr274 > lattice matrix and inverse
  DOUBLE PRECISION MLATTICE(3,3), MLATTICEINV(3,3)
  INTEGER NLATTICECOORDS

! hk286 - extra variables for minpermdist
  DOUBLE PRECISION :: D, DIST2, RMAT(3,3), RMATINV(3,3)
  DOUBLE PRECISION :: PP1(3*NATOMS), PP2(3*NATOMS)
  DOUBLE PRECISION :: NEWP(3), PNORM
  LOGICAL :: TEMPPERMDIST  ! sn402
  INTEGER :: TEMPNPERMGROUP, TEMPNPERMSIZE(3*NATOMS), TEMPPERMGROUP(3*NATOMS), TEMPNSETS(3*NATOMS), TEMPSETS(NATOMS,3)
  INTEGER :: FAILCOUNT = 0
  SAVEXCOORDS(:) = XCOORDS(:)

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

  ! sn402: added. There seems to be a bug (which appears to be compiler-dependent) with the coordinate transform when
  ! PERMDIST is set. So we temporarily turn it off. This will hopefully be removed when we figure out what's going on.
  IF(RIGIDMOLECULEST) THEN
    TEMPPERMDIST = PERMDIST
    PERMDIST = .FALSE.
  ENDIF

  IF (PERMDIST) THEN
     TEMPNPERMGROUP = NPERMGROUP
     TEMPNPERMSIZE(:) = NPERMSIZE(:)
     TEMPPERMGROUP(:) = PERMGROUP(:)
     TEMPNSETS(:) = NSETS(:)
     TEMPSETS(:,:) = SETS(:,:)
  ENDIF
! loop over all rigid bodies
  DO J1 = 1, NRIGIDBODY
     COM = 0.0D0
     MASS = 0.0D0
     XSAVE = XCOORDS(3*RIGIDGROUPS(1,J1)-2:3*RIGIDGROUPS(1,J1))
     ! calculate center of mass
     DO J2 = 1, NSITEPERBODY(J1)
        J9 = RIGIDGROUPS(J2, J1)

        IF(BULKT) THEN ! sn402. If a rigid body has got split up over the boundary
        ! of the box, the CoM calculated normally will be incorrect. We must
        ! ensure that the correct image of each atom is used.
        ! For each atom, measure its displacement from the previously considered atom and
        ! replace it with the coordinates of the image closest to that atom.
        ! This method should work fine for most systems, but if the rigid body has a length
        ! comparable with the box size (e.g. for large biological molecules) then it may have
        ! problems.
            IF(DEBUG .AND. ANY(NINT((XSAVE(:) - XCOORDS(3*J9-2:3*J9))/BULK_BOXVEC(:)) .NE. 0.0d0)) THEN
                WRITE(*,*) "In TRANSFORMCTORIGID, rigid body ", J1, "moving atom", J2,"due to PBCs."
            ENDIF
            XCOORDS(3*J9-2:3*J9) = XCOORDS(3*J9-2:3*J9) + &
                        & NINT((-XCOORDS(3*J9-2:3*J9) + XSAVE(:))/BULK_BOXVEC(:))*BULK_BOXVEC(:)
            XSAVE(:) = XCOORDS(3*J9-2:3*J9)
        ENDIF
        COM = COM + XCOORDS(3*J9-2:3*J9)*GR_WEIGHTS(RIGIDGROUPS(J2,J1))
        MASS = MASS + GR_WEIGHTS(RIGIDGROUPS(J2,J1))
     ENDDO
     COM = COM / MASS
     XRIGIDCOORDS(3*J1-2:3*J1) = COM

     DO J2 = 1, NSITEPERBODY(J1)
        J9 = RIGIDGROUPS(J2, J1)
        PP1(3*J2-2:3*J2) = XCOORDS(3*J9-2:3*J9) - COM
        PP2(3*J2-2:3*J2) = SITESRIGIDBODY(J2,:,J1)
     ENDDO

     IF (PERMDIST) THEN
        NPERMGROUP = LRBNPERMGROUP(J1)
        NPERMSIZE(:) = LRBNPERMSIZE(J1,:)
        PERMGROUP(:) = LRBPERMGROUP(J1,:)
        NSETS(:) = LRBNSETS(J1,:)
        SETS(:,:) = LRBSETS(J1,:,:)
     ENDIF

    ! PP2 are the com-frame coordinates for the reference structure of this rigid fragment
    ! PP1 are the com-frame coordinates for the corresponding fragment in XCOORDS.

     RIGIDINIT = .FALSE. ! sn402
     CALL MINPERMDIST(PP2(1:3*NSITEPERBODY(J1)),PP1(1:3*NSITEPERBODY(J1)),NSITEPERBODY(J1),.FALSE., &
          PARAM1,PARAM2,PARAM3,.FALSE.,.FALSE.,D,DIST2,.FALSE.,RMAT)
     RIGIDINIT = .TRUE. ! sn402

     RMATINV(1,1:3) = RMAT(1:3,1)
     RMATINV(2,1:3) = RMAT(1:3,2)
     RMATINV(3,1:3) = RMAT(1:3,3)
!     XRIGIDCOORDS(3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1) = rot_mx2aa(RMATINV)
     NEWP = rot_mx2aa(RMATINV)

!     PNORM = NORM2(NEWP)  ! NAG doesn't like this.
     PNORM = 0.0D0
     DO J2 = 1,3
        PNORM = PNORM + NEWP(J2)*NEWP(J2)
     ENDDO
     PNORM = SQRT(PNORM)
     IF(PNORM .GT. 8*ATAN(1.0D0)) THEN
        IF(DEBUG) WRITE(*,*) "Subtracting off rotation of 2PI from AA vector"
        NEWP(:) = NEWP(:)/PNORM
        NEWP(:) = NEWP(:)*(PNORM-8*ATAN(1.0D0))
     ELSE IF(PNORM .GT. 4*ATAN(1.0D0)) THEN
        IF(DEBUG) WRITE(*,*) "Moving to a symmetry-related AA vector with minimum magnitude"
        NEWP(:) = -(NEWP(:)/PNORM)
        NEWP(:) = NEWP(:)*(8*ATAN(1.0D0)-PNORM)
     ENDIF
     XRIGIDCOORDS(3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1) = NEWP(:)

     IF ( D/NSITEPERBODY(J1) > 0.1D0 ) THEN
        PRINT *,  'Warning: Genrigid > mapping looks bad for RB no ', J1 
        PRINT *,  'Warning: Genrigid > average deviation per particle is', D/NSITEPERBODY(J1)
        PRINT *,  'Warning: Genrigid >  Often it is the permutation of the RB members, e.g. Hs in NH3'
        PRINT *,  'Warning: Genrigid >  Use perm.allow'
     ENDIF
  ENDDO

  IF (PERMDIST) THEN
     NPERMGROUP = TEMPNPERMGROUP
     NPERMSIZE(:) = TEMPNPERMSIZE(:)
     PERMGROUP(:) = TEMPPERMGROUP(:)
     NSETS(:) = TEMPNSETS(:)
     SETS(:,:) = TEMPSETS(:,:)
  ENDIF

  ! sn402: see previous comment
  IF(RIGIDMOLECULEST) PERMDIST = TEMPPERMDIST

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

! sn402: the remainder of the subroutine is a safety check.
  IF(DEBUG .AND. .NOT. (AMBERT .OR. AMBER12T .OR. NABT)) THEN
    CALL TRANSFORMRIGIDTOC(1,NRIGIDBODY,TEMPXCOORDS,XRIGIDCOORDS)
    DO J1 = 1, 3*NATOMS
        IF(ABS(TEMPXCOORDS(J1)-SAVEXCOORDS(J1)) .GT. 1.0E-7) THEN
            IF(.NOT.(BULKT) .OR. (ABS(ABS(TEMPXCOORDS(J1)-SAVEXCOORDS(J1))-BULK_BOXVEC(MOD(J1-1,3)+1)) .GT. 1.0E-7)) THEN
                WRITE(*,*) "Warning, coordinate transform may have failed for coordinate ", J1
                WRITE(*,*) "Original coordinate:", SAVEXCOORDS(J1)
                WRITE(*,*) "New coordinate:", TEMPXCOORDS(J1)
                WRITE(*,*) "Difference:", TEMPXCOORDS(J1)-SAVEXCOORDS(J1)
                FAILCOUNT = FAILCOUNT+1
            ENDIF
        ENDIF
!        IF(BULKT) THEN
!            IF(ABS(TEMPXCOORDS(J1)-SAVEXCOORDS(J1)-BULK_BOXVEC(MOD(J1-1,3)+1)) .GT. 1.0E-7) THEN
!                WRITE(*,*) "Warning, coordinate transform may have failed."
!                WRITE(*,*) "Original coordinate:", SAVEXCOORDS(J1)
!                WRITE(*,*) "New coordinate:", TEMPXCOORDS(J1)
!                WRITE(*,*) "Difference:", TEMPXCOORDS(J1)-SAVEXCOORDS(J1)
!                FAILCOUNT = FAILCOUNT+1
!            ENDIF
!        ELSE
!            IF(ABS(TEMPXCOORDS(J1)-SAVEXCOORDS(J1)) .GT. 1.0E-7) THEN
!                WRITE(*,*) "Warning, coordinate transform may have failed."
!                WRITE(*,*) "Original coordinate:", SAVEXCOORDS(J1)
!                WRITE(*,*) "New coordinate:", TEMPXCOORDS(J1)
!                WRITE(*,*) "Difference:", TEMPXCOORDS(J1)-SAVEXCOORDS(J1)
!                FAILCOUNT = FAILCOUNT+1
!            ENDIF
!        ENDIF
    ENDDO
    IF (FAILCOUNT > 0) THEN
        WRITE(*,*) "Transformation failed ", FAILCOUNT, "times. Aborting."
        STOP 99
    ENDIF
  ENDIF

END SUBROUTINE TRANSFORMCTORIGID

!-----------------------------------------------------------

SUBROUTINE TRANSFORMGRAD (G, XR, GR)
  
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE
  
  INTEGER          :: J1, J2, J9
  DOUBLE PRECISION, INTENT(IN) :: G(3*NATOMS), XR(DEGFREEDOMS)
  DOUBLE PRECISION, INTENT(OUT) :: GR(DEGFREEDOMS)
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

SUBROUTINE AACONVERGENCE (G, XR, GR, RMS)
! Atomistic gradient, rigid coords, rigid gradient

  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE

  INTEGER          :: J1, J2, J9
  DOUBLE PRECISION :: G(3*NATOMS), XR(DEGFREEDOMS), GR(DEGFREEDOMS)
  DOUBLE PRECISION :: PI(3)
  DOUBLE PRECISION :: RMI(3,3), DRMI1(3,3), DRMI2(3,3), DRMI3(3,3)
  DOUBLE PRECISION :: TORQUE(3)
  DOUBLE PRECISION :: RMI0(3,3), DRMI10(3,3), DRMI20(3,3), DRMI30(3,3)
  DOUBLE PRECISION :: DR1(3),DR2(3),DR3(3), RMI3(3,3), RMS

  RMS = 0.0D0
  PI = (/0.0D0, 0.0D0, 0.0D0/)
  CALL RMDRVT(PI, RMI0, DRMI10, DRMI20, DRMI30, .TRUE.)

  ! For each rigid body
  DO J1 = 1, NRIGIDBODY
     ! Get the AA vector and its derivatives
     PI = XR(3*NRIGIDBODY+3*J1-2 : 3*NRIGIDBODY+3*J1)
     CALL RMDRVT(PI, RMI, DRMI1, DRMI2, DRMI3, .FALSE.)

     TORQUE(:) = 0.0D0
     DO J2 = 1, NSITEPERBODY(J1)
        J9 = RIGIDGROUPS(J2, J1)
        DR1(:) = MATMUL(DRMI10,MATMUL(RMI,SITESRIGIDBODY(J2,:,J1)))
        DR2(:) = MATMUL(DRMI20,MATMUL(RMI,SITESRIGIDBODY(J2,:,J1)))
        DR3(:) = MATMUL(DRMI30,MATMUL(RMI,SITESRIGIDBODY(J2,:,J1)))
        TORQUE(1) = TORQUE(1) + DOT_PRODUCT(G(3*J9-2:3*J9),DR1(:))
        TORQUE(2) = TORQUE(2) + DOT_PRODUCT(G(3*J9-2:3*J9),DR2(:))
        TORQUE(3) = TORQUE(3) + DOT_PRODUCT(G(3*J9-2:3*J9),DR3(:))
     ENDDO
     TORQUE = MATMUL(TRANSPOSE(RMI), TORQUE)
     ! IINVERSE contains the inverse eigenvalues of the unweighted tensor of gyration for this rigid body.
     RMS = RMS + DOT_PRODUCT(TORQUE, MATMUL(TRANSPOSE(IINVERSE(J1,:,:)),TORQUE))
!     RMI3 = MATMUL( RMI, IINVERSE(J1,:,:))
!     RMI3 = MATMUL( RMI3, TRANSPOSE(RMI) )
!     RMS = RMS + DOT_PRODUCT(TORQUE, MATMUL(TRANSPOSE(RMI3),TORQUE))

     ! Add the component of the rms force due to translational degrees of freedom
     RMS = RMS + 1.0D0/NSITEPERBODY(J1) * DOT_PRODUCT(GR(3*J1-2:3*J1),GR(3*J1-2:3*J1))
  ENDDO

  ! Deal with free atoms (translational component only)
  IF (DEGFREEDOMS > 6 * NRIGIDBODY) THEN
     DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
        RMS = RMS + DOT_PRODUCT(GR(6*NRIGIDBODY + 3*J1-2:6*NRIGIDBODY + 3*J1),GR(6*NRIGIDBODY + 3*J1-2:6*NRIGIDBODY + 3*J1))
     ENDDO
  ENDIF

  RMS=MAX(DSQRT(RMS/(3*NATOMS)),1.0D-100)

END SUBROUTINE AACONVERGENCE

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

!--------------------------------------------------------------

SUBROUTINE RBDET(A, DET)

  IMPLICIT NONE
  DOUBLE PRECISION :: A (3,3), DET

  DET = A(1,1)*(A(2,2)*A(3,3)-A(2,3)*A(3,2)) - A(1,2)*(A(2,1)*A(3,3)-A(2,3)*A(3,1)) + A(1,3)*(A(2,1)*A(3,2)-A(2,2)*A(3,1))

END SUBROUTINE RBDET

SUBROUTINE GENRIGID_IMAGE_CTORIGID(NIMAGE, XYZ)

  USE COMMONS, only: NATOMS
  IMPLICIT NONE
  
  INTEGER :: I, NIMAGE, NOPT 
  DOUBLE PRECISION :: XCOORDS(3*NATOMS), XRIGIDCOORDS (DEGFREEDOMS)
  DOUBLE PRECISION :: XYZ(3*NATOMS*(NIMAGE+2))

! XYZ holds the atomistic coordinates of the NIMAGE images and the two end points.
  NOPT = 3*NATOMS ! The number of (atomistic) coordinates to be optimised.
  DO I=1,NIMAGE+2
     XCOORDS(1:NOPT)=XYZ(NOPT*(I-1)+1:NOPT*I)
     CALL TRANSFORMCTORIGID (XCOORDS, XRIGIDCOORDS)
     ! Fill up the first DEGFREEDOMS coordinates in this subsection of XYZ
     ! with the correct rigid-body coordinates.
     XYZ(NOPT*(I-1)+1:NOPT*(I-1)+DEGFREEDOMS) = XRIGIDCOORDS(1:DEGFREEDOMS)
     ! Fill up the excess with 0s.
     XYZ(NOPT*(I-1)+DEGFREEDOMS+1:NOPT*(I-1)+NOPT) = 0.0D0
  ENDDO

END SUBROUTINE

SUBROUTINE GENRIGID_IMAGE_RIGIDTOC(NIMAGE, XYZ)

  USE COMMONS, only: NATOMS
  IMPLICIT NONE
  
  INTEGER :: I, NIMAGE, NOPT
  DOUBLE PRECISION :: XCOORDS(3*NATOMS), XRIGIDCOORDS (DEGFREEDOMS)
  DOUBLE PRECISION :: XYZ(3*NATOMS*(NIMAGE+2))

  NOPT = 3*NATOMS
  DO I=1,NIMAGE+2
     XRIGIDCOORDS(1:DEGFREEDOMS) = XYZ( NOPT*(I-1)+1:NOPT*(I-1)+DEGFREEDOMS )
     CALL TRANSFORMRIGIDTOC (1, NRIGIDBODY, XCOORDS, XRIGIDCOORDS)
     XYZ(NOPT*(I-1)+1:NOPT*I) = XCOORDS(1:NOPT)
  ENDDO

END SUBROUTINE

!-----------------------------------------------------------
! Calculate the Hessian in rigid body coordinates from the Gradient and Hessian in Cartesians.
! Follows the procedure outlined in Mochizuki et al, PCCP 16 (2014)
SUBROUTINE TRANSFORMHESSIAN (H, G, XR, HR, RBAANORMALMODET)
  
  USE COMMONS, ONLY: NATOMS, DEBUG
  IMPLICIT NONE
  
  INTEGER          :: J1, J2, J3, J4, J8, J9, K, L
  DOUBLE PRECISION, INTENT(IN) :: G(3*NATOMS), H(3*NATOMS,3*NATOMS), XR(DEGFREEDOMS)
  DOUBLE PRECISION, INTENT(OUT) :: HR(DEGFREEDOMS,DEGFREEDOMS)
  DOUBLE PRECISION :: PI(3)
  DOUBLE PRECISION :: AD2R11(3),AD2R22(3),AD2R33(3),AD2R12(3),AD2R23(3),AD2R31(3) 
  DOUBLE PRECISION :: ADR1(3),ADR2(3),ADR3(3) 
  DOUBLE PRECISION :: ARMI(3,3), ADRMI1(3,3), ADRMI2(3,3), ADRMI3(3,3)
  DOUBLE PRECISION :: AD2RMI11(3,3), AD2RMI22(3,3), AD2RMI33(3,3)
  DOUBLE PRECISION :: AD2RMI12(3,3), AD2RMI23(3,3), AD2RMI31(3,3)
  DOUBLE PRECISION :: BDR1(3),BDR2(3),BDR3(3) 
  DOUBLE PRECISION :: BRMI(3,3), BDRMI1(3,3), BDRMI2(3,3), BDRMI3(3,3)
  DOUBLE PRECISION :: BD2RMI11(3,3), BD2RMI22(3,3), BD2RMI33(3,3)
  DOUBLE PRECISION :: BD2RMI12(3,3), BD2RMI23(3,3), BD2RMI31(3,3)
  LOGICAL :: GTEST, STEST, RBAANORMALMODET
  DOUBLE PRECISION :: RMI0(3,3), DRMI10(3,3), DRMI20(3,3), DRMI30(3,3)
  DOUBLE PRECISION :: D2RMI10(3,3), D2RMI20(3,3), D2RMI30(3,3), D2RMI120(3,3), D2RMI230(3,3), D2RMI310(3,3)
  
  GTEST = .TRUE.
  STEST = .TRUE.
  HR(:,:) = 0.0D0

  IF(DEBUG) THEN
     IF (.NOT.(ANY(ABS(G(DEGFREEDOMS+1:3*NATOMS)) .GT. 1.0E-8))) THEN
           WRITE(*,*) "genrigid> Error: Gradient passed into TRANSFORM_HESSIAN seems to be in rigid body coordinates."
           WRITE(*,*) G(DEGFREEDOMS+1:)
           STOP
     ENDIF
  ENDIF

  IF ( RBAANORMALMODET ) THEN
     PI = (/0.0D0, 0.0D0, 0.0D0/)
     CALL RMDFAS(PI, RMI0, DRMI10, DRMI20, DRMI30, D2RMI10, D2RMI20, D2RMI30, D2RMI120, D2RMI230, D2RMI310, GTEST, STEST)
  ENDIF

  DO J1 = 1, NRIGIDBODY

     PI = XR(3*NRIGIDBODY+3*J1-2 : 3*NRIGIDBODY+3*J1)
     CALL RMDFAS(PI, ARMI, ADRMI1, ADRMI2, ADRMI3, AD2RMI11, AD2RMI22, AD2RMI33, AD2RMI12, AD2RMI23, AD2RMI31, GTEST, STEST)
    
     DO J2 = J1, NRIGIDBODY

        PI = XR(3*NRIGIDBODY+3*J2-2 : 3*NRIGIDBODY+3*J2)
        CALL RMDFAS(PI, BRMI, BDRMI1, BDRMI2, BDRMI3, BD2RMI11, BD2RMI22, BD2RMI33, BD2RMI12, BD2RMI23, BD2RMI31, GTEST, STEST)
             
        DO J3 = 1, NSITEPERBODY(J1)
           J8 = RIGIDGROUPS(J3, J1)

           DO J4 = 1, NSITEPERBODY(J2)
              J9 = RIGIDGROUPS(J4, J2)

! hk286 > translation
              HR(3*J1-2:3*J1, 3*J2-2:3*J2) = HR(3*J1-2:3*J1, 3*J2-2:3*J2) + H(3*J8-2:3*J8, 3*J9-2:3*J9)

! hk286 > rotations
              IF ( RBAANORMALMODET ) THEN
                 ADR1(:) = MATMUL(DRMI10,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 ADR2(:) = MATMUL(DRMI20,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 ADR3(:) = MATMUL(DRMI30,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 BDR1(:) = MATMUL(DRMI10,MATMUL(BRMI,SITESRIGIDBODY(J4,:,J2)))
                 BDR2(:) = MATMUL(DRMI20,MATMUL(BRMI,SITESRIGIDBODY(J4,:,J2)))
                 BDR3(:) = MATMUL(DRMI30,MATMUL(BRMI,SITESRIGIDBODY(J4,:,J2)))
              ELSE
                 ADR1(:) = MATMUL(ADRMI1,SITESRIGIDBODY(J3,:,J1))
                 ADR2(:) = MATMUL(ADRMI2,SITESRIGIDBODY(J3,:,J1))
                 ADR3(:) = MATMUL(ADRMI3,SITESRIGIDBODY(J3,:,J1))
                 BDR1(:) = MATMUL(BDRMI1,SITESRIGIDBODY(J4,:,J2))
                 BDR2(:) = MATMUL(BDRMI2,SITESRIGIDBODY(J4,:,J2))
                 BDR3(:) = MATMUL(BDRMI3,SITESRIGIDBODY(J4,:,J2))
              ENDIF

! hk286 - mixed translation rotation
              HR(3*J1-2, 3*NRIGIDBODY+3*J2-2)=HR(3*J1-2, 3*NRIGIDBODY+3*J2-2)+DOT_PRODUCT(H(3*J8-2, 3*J9-2:3*J9),BDR1(:))
              HR(3*J1-1, 3*NRIGIDBODY+3*J2-2)=HR(3*J1-1, 3*NRIGIDBODY+3*J2-2)+DOT_PRODUCT(H(3*J8-1, 3*J9-2:3*J9),BDR1(:))
              HR(3*J1  , 3*NRIGIDBODY+3*J2-2)=HR(3*J1  , 3*NRIGIDBODY+3*J2-2)+DOT_PRODUCT(H(3*J8  , 3*J9-2:3*J9),BDR1(:))
              HR(3*J1-2, 3*NRIGIDBODY+3*J2-1)=HR(3*J1-2, 3*NRIGIDBODY+3*J2-1)+DOT_PRODUCT(H(3*J8-2, 3*J9-2:3*J9),BDR2(:))
              HR(3*J1-1, 3*NRIGIDBODY+3*J2-1)=HR(3*J1-1, 3*NRIGIDBODY+3*J2-1)+DOT_PRODUCT(H(3*J8-1, 3*J9-2:3*J9),BDR2(:))
              HR(3*J1  , 3*NRIGIDBODY+3*J2-1)=HR(3*J1  , 3*NRIGIDBODY+3*J2-1)+DOT_PRODUCT(H(3*J8  , 3*J9-2:3*J9),BDR2(:))
              HR(3*J1-2, 3*NRIGIDBODY+3*J2  )=HR(3*J1-2, 3*NRIGIDBODY+3*J2  )+DOT_PRODUCT(H(3*J8-2, 3*J9-2:3*J9),BDR3(:))
              HR(3*J1-1, 3*NRIGIDBODY+3*J2  )=HR(3*J1-1, 3*NRIGIDBODY+3*J2  )+DOT_PRODUCT(H(3*J8-1, 3*J9-2:3*J9),BDR3(:))
              HR(3*J1  , 3*NRIGIDBODY+3*J2  )=HR(3*J1  , 3*NRIGIDBODY+3*J2  )+DOT_PRODUCT(H(3*J8  , 3*J9-2:3*J9),BDR3(:))        
              
              IF (J2 > J1) THEN
                 HR(3*J2-2, 3*NRIGIDBODY+3*J1-2) = HR(3*J2-2, 3*NRIGIDBODY+3*J1-2) &
                      + DOT_PRODUCT(H(3*J9-2, 3*J8-2:3*J8),ADR1(:))
                 HR(3*J2-1, 3*NRIGIDBODY+3*J1-2) = HR(3*J2-1, 3*NRIGIDBODY+3*J1-2) &
                      + DOT_PRODUCT(H(3*J9-1, 3*J8-2:3*J8),ADR1(:))
                 HR(3*J2  , 3*NRIGIDBODY+3*J1-2) = HR(3*J2  , 3*NRIGIDBODY+3*J1-2) &
                      + DOT_PRODUCT(H(3*J9  , 3*J8-2:3*J8),ADR1(:))
                 HR(3*J2-2, 3*NRIGIDBODY+3*J1-1) = HR(3*J2-2, 3*NRIGIDBODY+3*J1-1) &
                      + DOT_PRODUCT(H(3*J9-2, 3*J8-2:3*J8),ADR2(:))
                 HR(3*J2-1, 3*NRIGIDBODY+3*J1-1) = HR(3*J2-1, 3*NRIGIDBODY+3*J1-1) &
                      + DOT_PRODUCT(H(3*J9-1, 3*J8-2:3*J8),ADR2(:))
                 HR(3*J2  , 3*NRIGIDBODY+3*J1-1) = HR(3*J2  , 3*NRIGIDBODY+3*J1-1) &
                      + DOT_PRODUCT(H(3*J9  , 3*J8-2:3*J8),ADR2(:))
                 HR(3*J2-2, 3*NRIGIDBODY+3*J1  ) = HR(3*J2-2, 3*NRIGIDBODY+3*J1  ) &
                      + DOT_PRODUCT(H(3*J9-2, 3*J8-2:3*J8),ADR3(:))
                 HR(3*J2-1, 3*NRIGIDBODY+3*J1  ) = HR(3*J2-1, 3*NRIGIDBODY+3*J1  ) &
                      + DOT_PRODUCT(H(3*J9-1, 3*J8-2:3*J8),ADR3(:))
                 HR(3*J2  , 3*NRIGIDBODY+3*J1  ) = HR(3*J2  , 3*NRIGIDBODY+3*J1  ) &
                      + DOT_PRODUCT(H(3*J9  , 3*J8-2:3*J8),ADR3(:))        
              ENDIF
! hk286 - double rotation
              DO K = 1, 3
                 DO L = 1, 3
                    HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-2)=HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-2)+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR1(K) * BDR1(L)
                    HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2-2)=HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2-2)+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR2(K) * BDR1(L)
                    HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2-2)=HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2-2)+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR3(K) * BDR1(L)
                    HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-1)=HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-1)+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR1(K) * BDR2(L)
                    HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2-1)=HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2-1)+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR2(K) * BDR2(L)
                    HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2-1)=HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2-1)+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR3(K) * BDR2(L)
                    HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2  )=HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2  )+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR1(K) * BDR3(L)
                    HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2  )=HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2  )+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR2(K) * BDR3(L)
                    HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2  )=HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2  )+&
                    & H(3*J8-3+K, 3*J9-3+L) * ADR3(K) * BDR3(L)   
                 ENDDO
              ENDDO
           ENDDO ! Loop over J4

           IF (J1 .EQ. J2) THEN
              IF ( RBAANORMALMODET ) THEN
                 AD2R11(:) = MATMUL(D2RMI10, MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 AD2R22(:) = MATMUL(D2RMI20, MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 AD2R33(:) = MATMUL(D2RMI30, MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 AD2R12(:) = MATMUL(D2RMI120,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 AD2R23(:) = MATMUL(D2RMI230,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
                 AD2R31(:) = MATMUL(D2RMI310,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
              ELSE
                 AD2R11(:) = MATMUL(AD2RMI11,SITESRIGIDBODY(J3,:,J1))
                 AD2R22(:) = MATMUL(AD2RMI22,SITESRIGIDBODY(J3,:,J1))
                 AD2R33(:) = MATMUL(AD2RMI33,SITESRIGIDBODY(J3,:,J1))
                 AD2R12(:) = MATMUL(AD2RMI12,SITESRIGIDBODY(J3,:,J1))
                 AD2R23(:) = MATMUL(AD2RMI23,SITESRIGIDBODY(J3,:,J1))
                 AD2R31(:) = MATMUL(AD2RMI31,SITESRIGIDBODY(J3,:,J1))
              ENDIF

              ! p_x, p_x
              HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-2) = HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-2) &
                   + DOT_PRODUCT(G(3*J8-2:3*J8),AD2R11(:))
              ! p_y, p_y
              HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2-1) = HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2-1) &
                   + DOT_PRODUCT(G(3*J8-2:3*J8),AD2R22(:))
              ! p_z, p_z
              HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2  ) = HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2  ) &
                   + DOT_PRODUCT(G(3*J8-2:3*J8),AD2R33(:))
              ! p_x, p_y
              HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-1) = HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J2-1) &
                   + DOT_PRODUCT(G(3*J8-2:3*J8),AD2R12(:))
              ! p_y, p_z
              HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2  ) = HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J2  ) &
                   + DOT_PRODUCT(G(3*J8-2:3*J8),AD2R23(:))
              ! p_z, p_x
              HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2-2) = HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J2-2) &
                   + DOT_PRODUCT(G(3*J8-2:3*J8),AD2R31(:))
           ENDIF
        ENDDO ! Loop over J3 (i.e. J8)

        IF (J1 .EQ. J2) THEN
           HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J1-2) = HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J1-1) 
           HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J1-1) = HR(3*NRIGIDBODY+3*J1-1, 3*NRIGIDBODY+3*J1  )
           HR(3*NRIGIDBODY+3*J1-2, 3*NRIGIDBODY+3*J1  ) = HR(3*NRIGIDBODY+3*J1  , 3*NRIGIDBODY+3*J1-2)
           HR(3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1, 3*J1-2:3*J1) = &
                TRANSPOSE(HR(3*J1-2:3*J1, 3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1))
        ELSE
           HR(3*J2-2:3*J2, 3*J1-2:3*J1) = TRANSPOSE(HR(3*J1-2:3*J1, 3*J2-2:3*J2))
           HR(3*NRIGIDBODY+3*J2-2:3*NRIGIDBODY+3*J2, 3*J1-2:3*J1) = &
                TRANSPOSE(HR(3*J1-2:3*J1, 3*NRIGIDBODY+3*J2-2:3*NRIGIDBODY+3*J2))
           HR(3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1, 3*J2-2:3*J2) = &
                TRANSPOSE(HR(3*J2-2:3*J2, 3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1))
           HR(3*NRIGIDBODY+3*J2-2:3*NRIGIDBODY+3*J2, 3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1) = &
                TRANSPOSE(HR(3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1, 3*NRIGIDBODY+3*J2-2:3*NRIGIDBODY+3*J2))
        ENDIF

     ENDDO
  ENDDO

  DO J1 = 1, NRIGIDBODY

     DO J2 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
            
        J9 = RIGIDSINGLES(J2)
        PI = XR(3*NRIGIDBODY+3*J1-2 : 3*NRIGIDBODY+3*J1)
        CALL RMDFAS(PI, ARMI, ADRMI1, ADRMI2, ADRMI3, AD2RMI11, AD2RMI22, AD2RMI33, AD2RMI12, AD2RMI23, AD2RMI31, GTEST, STEST)
        
        DO J3 = 1, NSITEPERBODY(J1)
           J8 = RIGIDGROUPS(J3, J1)
                     
! hk286 > translation
           HR(3*J1-2:3*J1, 6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2) = HR(3*J1-2:3*J1, 6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2) &
                + H(3*J8-2:3*J8, 3*J9-2:3*J9)

! hk286 > rotations
           IF ( RBAANORMALMODET ) THEN
              ADR1(:) = MATMUL(DRMI10,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
              ADR2(:) = MATMUL(DRMI20,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
              ADR3(:) = MATMUL(DRMI30,MATMUL(ARMI,SITESRIGIDBODY(J3,:,J1)))
           ELSE
              ADR1(:) = MATMUL(ADRMI1,SITESRIGIDBODY(J3,:,J1))
              ADR2(:) = MATMUL(ADRMI2,SITESRIGIDBODY(J3,:,J1))
              ADR3(:) = MATMUL(ADRMI3,SITESRIGIDBODY(J3,:,J1))
           ENDIF
           HR(3*NRIGIDBODY+3*J1-2, 6*NRIGIDBODY+3*J2-2) = HR(3*NRIGIDBODY+3*J1-2, 6*NRIGIDBODY+3*J2-2) &
                + DOT_PRODUCT(H(3*J9-2, 3*J8-2:3*J8),ADR1(:))
           HR(3*NRIGIDBODY+3*J1-2, 6*NRIGIDBODY+3*J2-1) = HR(3*NRIGIDBODY+3*J1-2, 6*NRIGIDBODY+3*J2-1) &
                + DOT_PRODUCT(H(3*J9-1, 3*J8-2:3*J8),ADR1(:))
           HR(3*NRIGIDBODY+3*J1-2, 6*NRIGIDBODY+3*J2  ) = HR(3*NRIGIDBODY+3*J1-2, 6*NRIGIDBODY+3*J2  ) &
                + DOT_PRODUCT(H(3*J9  , 3*J8-2:3*J8),ADR1(:))
           HR(3*NRIGIDBODY+3*J1-1, 6*NRIGIDBODY+3*J2-2) = HR(3*NRIGIDBODY+3*J1-1, 6*NRIGIDBODY+3*J2-2) &
                + DOT_PRODUCT(H(3*J9-2, 3*J8-2:3*J8),ADR2(:))
           HR(3*NRIGIDBODY+3*J1-1, 6*NRIGIDBODY+3*J2-1) = HR(3*NRIGIDBODY+3*J1-1, 6*NRIGIDBODY+3*J2-1) &
                + DOT_PRODUCT(H(3*J9-1, 3*J8-2:3*J8),ADR2(:))
           HR(3*NRIGIDBODY+3*J1-1, 6*NRIGIDBODY+3*J2  ) = HR(3*NRIGIDBODY+3*J1-1, 6*NRIGIDBODY+3*J2  ) &
                + DOT_PRODUCT(H(3*J9  , 3*J8-2:3*J8),ADR2(:))
           HR(3*NRIGIDBODY+3*J1  , 6*NRIGIDBODY+3*J2-2) = HR(3*NRIGIDBODY+3*J1  , 6*NRIGIDBODY+3*J2-2) &
                + DOT_PRODUCT(H(3*J9-2, 3*J8-2:3*J8),ADR3(:))
           HR(3*NRIGIDBODY+3*J1  , 6*NRIGIDBODY+3*J2-1) = HR(3*NRIGIDBODY+3*J1  , 6*NRIGIDBODY+3*J2-1) &
                + DOT_PRODUCT(H(3*J9-1, 3*J8-2:3*J8),ADR3(:))
           HR(3*NRIGIDBODY+3*J1  , 6*NRIGIDBODY+3*J2  ) = HR(3*NRIGIDBODY+3*J1  , 6*NRIGIDBODY+3*J2  ) &
                + DOT_PRODUCT(H(3*J9  , 3*J8-2:3*J8),ADR3(:))        
        ENDDO
        
        HR(6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2, 3*J1-2:3*J1) = &
             TRANSPOSE(HR(3*J1-2:3*J1, 6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2))
        HR(6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2, 3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1) = &
             TRANSPOSE(HR(3*NRIGIDBODY+3*J1-2:3*NRIGIDBODY+3*J1, 6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2))
     ENDDO
  ENDDO

  DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
     J8 = RIGIDSINGLES(J1)
     DO J2 = J1, (DEGFREEDOMS - 6*NRIGIDBODY)/3            
        J9 = RIGIDSINGLES(J2)                     
        HR(6*NRIGIDBODY+3*J1-2:6*NRIGIDBODY+3*J1, 6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2) = H(3*J8-2:3*J8, 3*J9-2:3*J9)
        HR(6*NRIGIDBODY+3*J2-2:6*NRIGIDBODY+3*J2, 6*NRIGIDBODY+3*J1-2:6*NRIGIDBODY+3*J1) = H(3*J9-2:3*J9, 3*J8-2:3*J8)
     ENDDO
  ENDDO

  ! sn402: safety check
  IF(DEBUG) THEN
      DO J1 = 1,DEGFREEDOMS
         DO J2 = 1,DEGFREEDOMS
             IF(ABS(HR(J1,J2)-HR(J2,J1)) .GT. 1.0E-7) THEN
                 write(*,*) "transformHessian> Asymmetric Hessian, coords ", J1, J2
             ENDIF
         ENDDO
      ENDDO
  ENDIF

END SUBROUTINE TRANSFORMHESSIAN


!SUBROUTINE NABSECDER(OLDX,STEST)
!use modhess
!USE COMMONS, ONLY: NATOMS

!implicit none

!DOUBLE PRECISION :: V(3*NATOMS),ENERGY,XDUMMY(3*NATOMS),OLDX(3*NATOMS),VTEMP(2,3*NATOMS),ksi
!DOUBLE PRECISION :: XRIGIDCOORDS(DEGFREEDOMS), XCOORDS(3*NATOMS), XRIGIDGRAD(DEGFREEDOMS)
!INTEGER    :: i,j,k
!LOGICAL    :: GTEST,STEST

!ksi=0.00001D0
!XDUMMY(:)=OLDX(:)

!IF (.NOT.ALLOCATED(HESS)) ALLOCATE(HESS(3*NATOMS,3*NATOMS))
!DO i=1,DEGFREEDOMS

!   XDUMMY(i)=XDUMMY(i)-ksi 
!   XRIGIDCOORDS(1:DEGFREEDOMS) = XDUMMY(1:DEGFREEDOMS)
!   CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, XCOORDS, XRIGIDCOORDS)               
!   CALL MME(ENERGY,XCOORDS,V,1)
!   CALL TRANSFORMGRAD(V, XRIGIDCOORDS, XRIGIDGRAD)
!   VTEMP(1,1:DEGFREEDOMS)=XRIGIDGRAD(1:DEGFREEDOMS)
!!   CALL POTENTIAL (XDUMMY,ENERGY,V,.TRUE.,.FALSE.,1.0D-6,.FALSE.,.FALSE.)
!!   VTEMP(1,:)=V(:)
   
!   XDUMMY(i)=XDUMMY(i)+2.0D0*ksi   
!   XRIGIDCOORDS(1:DEGFREEDOMS) = XDUMMY(1:DEGFREEDOMS)
!   CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, XCOORDS, XRIGIDCOORDS)               
!   CALL MME(ENERGY,XCOORDS,V,1)
!   CALL TRANSFORMGRAD(V, XRIGIDCOORDS, XRIGIDGRAD)
!   VTEMP(2,1:DEGFREEDOMS)=XRIGIDGRAD(1:DEGFREEDOMS)
!!   CALL POTENTIAL (XDUMMY,ENERGY,V,.TRUE.,.FALSE.,1.0D-6,.FALSE.,.FALSE.)
!!   VTEMP(2,:)=V(:)

!   XDUMMY(i)=XDUMMY(i)-ksi 
   
!   DO j=i,DEGFREEDOMS
!      HESS(i,j)=(VTEMP(2,j)-VTEMP(1,j))/(2.0D0*ksi)
!      HESS(j,i)=HESS(i,j)
!!      PRINT *, i, j, HESS(i,j)
!!      PRINT *, VTEMP(2,j), VTEMP(1,j)
!   END DO

!END DO

!END SUBROUTINE NABSECDER


!SUBROUTINE ATOMNABSECDER(OLDX,V,STEST)
!use modhess
!USE COMMONS, ONLY: NATOMS

!implicit none

!DOUBLE PRECISION :: V(3*NATOMS),ENERGY,XDUMMY(3*NATOMS),OLDX(3*NATOMS),VTEMP(2,3*NATOMS),ksi
!INTEGER    :: i,j,k
!LOGICAL    :: GTEST,STEST

!ksi=0.0000001D0
!XDUMMY(:)=OLDX(:)

!IF (.NOT.ALLOCATED(HESS)) ALLOCATE(HESS(3*NATOMS,3*NATOMS))
!DO i=1,3*NATOMS

!   XDUMMY(i)=XDUMMY(i)-ksi 
!   CALL MME(ENERGY,XDUMMY,V,1)
!   VTEMP(1,:)=V(:)
   
!   XDUMMY(i)=XDUMMY(i)+2.0D0*ksi   
!   CALL MME(ENERGY,XDUMMY,V,1)
!   VTEMP(2,:)=V(:)

!   XDUMMY(i)=XDUMMY(i)-ksi 
   
!   DO j=i,3*NATOMS
!      HESS(i,j)=(VTEMP(2,j)-VTEMP(1,j))/(2.0D0*ksi)
!      HESS(j,i)=HESS(i,j)
!!      PRINT *, i, j, HESS(i,j)
!!      PRINT *, VTEMP(2,j), VTEMP(1,j)
!   END DO

!END DO

!XDUMMY(:)=OLDX(:)
!CALL MME(ENERGY,XDUMMY,V,1)

!END SUBROUTINE ATOMNABSECDER

! sn402: New version of eigenvalue calculation using the metric tensor
! Following dx.doi.org/10.1021/ct400403y | J. Chem. Theory Comput. 2013, 9, 4026−4034
!----------------------------------------------------------------------------------------------
SUBROUTINE GENRIGID_NORMALMODES(X, ATOMMASS, DIAG, INFO)
  USE COMMONS
  USE MODHESS     ! For the Hessian
  USE MODCHARMM
  USE KEY
  IMPLICIT NONE

  ! Communicating with the rest of OPTIM
  DOUBLE PRECISION, INTENT(IN)  :: X(3*NATOMS), ATOMMASS(3*NATOMS) ! Input coordinates and atom masses
                                                ! Note, ATOMMASS is currently not being used.
  DOUBLE PRECISION, INTENT(OUT) :: DIAG(3*NATOMS) ! Output: will contain the normal mode frequencies
  LOGICAL                       :: ART  ! Temporary variable to save the value of ATOMRIGIDCOORDT
  LOGICAL                       :: RBT  ! Temporary variable to save the value of RBAANORMALMODET

  ! Constructing the metric tensor
  DOUBLE PRECISION :: P(3) ! Angle-axis vector of a particular rigid body
  DOUBLE PRECISION :: W ! Total weight of a particular rigid body
  DOUBLE PRECISION :: COM(3) ! Centre of mass for a particular rigid body
  DOUBLE PRECISION :: R(3,3), DR1(3,3), DR2(3,3), DR3(3,3) ! Rotation matrix and its derivatives for an RB
  DOUBLE PRECISION :: RB_IMT(6,6)  ! Inverse Metric Tensor (IMT) for a particular rigid body
  DOUBLE PRECISION :: LR(3,3) ! A temporary matrix used to invert the rotational component of the metric tensor
  DOUBLE PRECISION :: DET ! Determinant of the rotational component of the metric tensor
  DOUBLE PRECISION :: INVERSE_METRIC_TENSOR(DEGFREEDOMS,DEGFREEDOMS) ! For the whole system

  ! Constructing the Hessian
  DOUBLE PRECISION :: XRIGIDCOORDS(DEGFREEDOMS), XCOORDS(3*NATOMS) ! Temporary coords arrays
  DOUBLE PRECISION :: ENERGY, G(3*NATOMS), RMS  ! Effectively dummies to pass in to the POTENTIAL call.

  ! Determining the normal modes
  DOUBLE PRECISION :: METRICHESS(DEGFREEDOMS,DEGFREEDOMS) ! The inverse MT multiplied by the Hessian
  DOUBLE PRECISION :: EIG_REAL(DEGFREEDOMS), EIG_IM(DEGFREEDOMS) ! Real and Imaginary parts of the eigenvalues
  DOUBLE PRECISION :: DUMMY_EVECS(DEGFREEDOMS,DEGFREEDOMS) ! A dummy array for DGEEV
  DOUBLE PRECISION :: EVECS(DEGFREEDOMS,DEGFREEDOMS) ! To hold the normal mode eigenvectors

  ! Temporary variables and dummys
  INTEGER          :: DUMMY  ! Counts the number of atoms considered in previous RBs when assembling the IMT.
  INTEGER          :: J1, J2, J3
  INTEGER          :: NFAILS ! Counts the number of times a failure is recorded when assembling the IMT.

! The following required to call the LAPACK routines DGETRF, DGETRI and DGEEV
  INTEGER, INTENT(OUT) :: INFO  ! An exit code for the LAPACK functions
  INTEGER, PARAMETER   :: LWORK = 1000000 ! Arbitrary dimension for WORK
  DOUBLE PRECISION     :: WORK(LWORK) ! Internal arrays for the LAPACK functions
  INTEGER              :: PIVOTS(3)

! Initialize variables
  DIAG(:) = 0.0D0
  ART = ATOMRIGIDCOORDT
  RBT = RBAANORMALMODET

  INVERSE_METRIC_TENSOR(:,:) = 0.0D0
  EIG_REAL(:) = 0.0D0
  EIG_IM(:) = 0.0D0
  DUMMY_EVECS(:,:) = 0.0D0

  DUMMY = 0
  NFAILS = 0

! Ensure we are in rigid body coordinates
  IF ( ATOMRIGIDCOORDT .EQV. .TRUE. ) THEN
     CALL TRANSFORMCTORIGID (X, XRIGIDCOORDS)
     ATOMRIGIDCOORDT = .FALSE.
  ELSE
     XRIGIDCOORDS(1:DEGFREEDOMS) = X(1:DEGFREEDOMS)
  ENDIF

! Construct the inverse metric tensor for each body, then we will combine them to get the overall IMT.
  DO J1 = 1, NRIGIDBODY

    W = 0.0D0
    COM(:) = RCOM(J1,:)
    RB_IMT(:,:) = 0.0D0

    DO J2 = 1, NSITEPERBODY(J1)
        ! DUMMY holds the total number of atoms belonging to rigid bodies that have already been considered
!        W = W + ATOMMASS(DUMMY+J2) ! It would be better to be able to pass in masses from outside but for the moment
!                                     the ATMASS vector seems to be getting reset to 0 (probably in fetchz.f). Until
!                                     I fix that, we'll just use the weight values saved in genrigid.
        W = W + GR_WEIGHTS(DUMMY+J2)
    ENDDO

    ! Translational components first (eq. 5 in the paper). These are the components of the inverse matrix.
    DO J2 = 1, 3
        RB_IMT(J2, J2) = 1.0D0/W
    ENDDO

    ! Rotational components. We start by constructing the metric tensor for each body exactly as given in
    ! eq. 6 of the paper, then we invert it.
    !
    ! Note, this procedure of inverting the two submatrices separately only works if the metric tensor for each body
    ! is block-diagonal i.e. the mixing terms are zero (which they should be)

    P  = XRIGIDCOORDS(3*(NRIGIDBODY+J1)-2:3*(NRIGIDBODY+J1))
    CALL RMDRVT(P, R, DR1, DR2, DR3, .TRUE.)

    RB_IMT(4,4) = TRACE3(MATMUL(DR1,MATMUL(Sm,TRANSPOSE(DR1)))) ! Using the mass-weighted gyration tensor Sm
    RB_IMT(4,5) = TRACE3(MATMUL(DR1,MATMUL(Sm,TRANSPOSE(DR2)))) ! which was computed in SETUP_TENSORS.
    RB_IMT(4,6) = TRACE3(MATMUL(DR1,MATMUL(Sm,TRANSPOSE(DR3))))
    RB_IMT(5,5) = TRACE3(MATMUL(DR2,MATMUL(Sm,TRANSPOSE(DR2))))
    RB_IMT(5,6) = TRACE3(MATMUL(DR2,MATMUL(Sm,TRANSPOSE(DR3))))
    RB_IMT(6,6) = TRACE3(MATMUL(DR3,MATMUL(Sm,TRANSPOSE(DR3))))

    ! The rotational part of the metric tensor should be symmetric, so we usually just assume that it is.
    IF(DEBUG) THEN
        RB_IMT(5,4) = TRACE3(MATMUL(DR2,MATMUL(Sm,TRANSPOSE(DR1))))
        RB_IMT(6,4) = TRACE3(MATMUL(DR3,MATMUL(Sm,TRANSPOSE(DR1))))
        RB_IMT(6,5) = TRACE3(MATMUL(DR3,MATMUL(Sm,TRANSPOSE(DR2))))
        IF(ABS(RB_IMT(5,4)-RB_IMT(4,5)).GT.1.0e-8) WRITE(*,*) "Warning: Metric tensor appears to be unsymmetric for body",J1
        IF(ABS(RB_IMT(6,4)-RB_IMT(4,6)).GT.1.0e-8) WRITE(*,*) "Warning: Metric tensor appears to be unsymmetric for body",J1
        IF(ABS(RB_IMT(6,5)-RB_IMT(5,6)).GT.1.0e-8) WRITE(*,*) "Warning: Metric tensor appears to be unsymmetric for body",J1
    ELSE
        RB_IMT(5,4) = RB_IMT(4,5)
        RB_IMT(6,4) = RB_IMT(4,6)
        RB_IMT(6,5) = RB_IMT(5,6)
    ENDIF

    ! Something has gone badly wrong if the mixing terms are non-zero: this implies that COM is non-zero but it ought
    ! to be (0,0,0) by construction.
    IF(DEBUG) THEN
        RB_IMT(1:3,4) = 2.0D0*W*MATMUL(DR1,COM)
        RB_IMT(1:3,5) = 2.0D0*W*MATMUL(DR2,COM)
        RB_IMT(1:3,6) = 2.0D0*W*MATMUL(DR3,COM)
        RB_IMT(4,1:3) = 2.0D0*W*MATMUL(DR1,COM)
        RB_IMT(5,1:3) = 2.0D0*W*MATMUL(DR2,COM)
        RB_IMT(6,1:3) = 2.0D0*W*MATMUL(DR3,COM)
        IF(ANY(ABS(RB_IMT(1:3,4:6)).GT.1.0e-8) .OR. ANY(ABS(RB_IMT(4:6,1:3)).GT.1.0e-8)) THEN
            WRITE(*,*) "Error: mixing terms in Metric Tensor are non-zero for body", J1
            WRITE(*,*) "Upper-right block:"
            DO J3=1,3
                WRITE(*,*) RB_IMT(J3,4:6)
            ENDDO
            WRITE(*,*) "Lower-left block:"
            DO J3=4,6
                WRITE(*,*) RB_IMT(J3,1:3)
            ENDDO
            NFAILS = NFAILS+1  ! We will now terminate at the end of the tensor-construction step.
            ! We don't terminate immediately in order to count the number of failures that occur.
        ENDIF
    ENDIF

! We have now constructed the rotational part of the metric tensor.
!
! Next we need to invert the matrix. I've considered using a LAPACK routine, but we don't seem to have DGETRI
! in our library, and in any case it's probably faster to code it explicitly for a 3x3 matrix.

    ! Compute the determinant (use LR to save the original values in the matrix we're inverting)
    LR(:,:) = RB_IMT(4:6,4:6)
    DET = LR(1,1)*(LR(2,2)*LR(3,3)-LR(2,3)*LR(3,2)) - LR(1,2)*(LR(2,1)*LR(3,3)-LR(2,3)*LR(3,1)) + &
        & LR(1,3)*(LR(2,1)*LR(3,2)-LR(2,2)*LR(3,1))
    IF(ABS(DET).LT.1e-8) THEN
        WRITE(*,*) "Error: Metric Tensor is singular. At this point we ought to calculate the pseudoinverse &
                 & and use that instead, but it's not yet implemented."
        STOP
    ENDIF

    ! Perform the inversion (we're assuming a symmetric metric tensor but that ought to be what we have)
    RB_IMT(4,4) = (LR(2,2)*LR(3,3)-LR(2,3)*LR(3,2))/DET
    RB_IMT(4,5) = -(LR(2,1)*LR(3,3)-LR(2,3)*LR(3,1))/DET
    RB_IMT(4,6) = (LR(2,1)*LR(3,2)-LR(2,2)*LR(3,1))/DET

    ! Copy the previous result where possible, to save on calculations
    RB_IMT(5,4) = RB_IMT(4,5) ! = -(LR(1,2)*LR(3,3)-LR(1,3)*LR(3,2))/DET
    RB_IMT(5,5) = (LR(1,1)*LR(3,3)-LR(1,3)*LR(3,1))/DET
    RB_IMT(5,6) = -(LR(1,1)*LR(3,2)-LR(1,2)*LR(3,1))/DET

    RB_IMT(6,4) = RB_IMT(4,6) ! = (LR(1,2)*LR(2,3)-LR(1,3)*LR(2,2))/DET
    RB_IMT(6,5) = RB_IMT(5,6) ! = -(LR(1,1)*LR(2,3)-LR(1,3)*LR(2,1))/DET
    RB_IMT(6,6) = (LR(1,1)*LR(2,2)-LR(1,2)*LR(2,1))/DET

!    write(*,*) "After inversion, RB_IMT for rigid body ", J1
!    write(*,*) "UL:"
!    DO J3=1,3
!        write(*,*) RB_IMT(J3,1:3)
!    ENDDO
!    write(*,*) "LR:"
!    DO J3=4,6
!        write(*,*) RB_IMT(J3,4:6)
!    ENDDO

    ! We finally have the IMT for this particular body!
    ! Now put it into the whole system IMT (which is block-diagonal)
    ! Translational block first
    INVERSE_METRIC_TENSOR(3*(J1-1)+1:3*J1,3*(J1-1)+1:3*J1) = RB_IMT(1:3,1:3)
    ! Then rotational block
    INVERSE_METRIC_TENSOR(3*NRIGIDBODY+3*(J1-1)+1:3*NRIGIDBODY+3*J1,3*NRIGIDBODY+3*(J1-1)+1:3*NRIGIDBODY+3*J1) = RB_IMT(4:6,4:6)

    ! Update the number of atoms which have already been considered
    DUMMY = DUMMY + NSITEPERBODY(J1)
  ENDDO ! End loop over rigid bodies

  ! If we had any serious errors detected during the construction of the inverse metric tensor, stop now.
  IF(DEBUG .AND. NFAILS .GT. 0) THEN
    WRITE(*,*) "GENRIGID_NORMALMODES> Construction of the IMT failed ", NFAILS, "times"
    STOP
  ENDIF

  ! Now we need to compute the correct rigid body hessian
  RBAANORMALMODET = .FALSE.
  XCOORDS(1:DEGFREEDOMS) = XRIGIDCOORDS(1:DEGFREEDOMS)
  XCOORDS(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
  IF (ENDNUMHESS .OR. AMBERT .OR. AMBER12T) THEN
     CALL GENRIGID_MAKENUMHESS(XCOORDS,NATOMS,DEGFREEDOMS)
  ELSE
     CALL POTENTIAL(XCOORDS,ENERGY,G,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
  ENDIF

  ! We need to left-multiply the Hessian by the inverse MT we calculated earlier, then find its eigenvalues and
  ! eigenvectors. See eq.31 in the paper.
  METRICHESS(:,:) = MATMUL(INVERSE_METRIC_TENSOR,HESS(1:DEGFREEDOMS,1:DEGFREEDOMS))

  IF (DUMPV) THEN ! Save right-eigenvectors as well.
     CALL DGEEV('N','V', DEGFREEDOMS, METRICHESS(1:DEGFREEDOMS,1:DEGFREEDOMS), DEGFREEDOMS, EIG_REAL, EIG_IM, &
                & DUMMY_EVECS, DEGFREEDOMS, EVECS(1:DEGFREEDOMS,1:DEGFREEDOMS), DEGFREEDOMS, WORK, LWORK, INFO)
  ELSE
     CALL DGEEV('N','N', DEGFREEDOMS, METRICHESS(1:DEGFREEDOMS,1:DEGFREEDOMS), DEGFREEDOMS, EIG_REAL, EIG_IM, &
                & DUMMY_EVECS, DEGFREEDOMS, EVECS(1:DEGFREEDOMS,1:DEGFREEDOMS), DEGFREEDOMS, WORK, LWORK, INFO)
  ENDIF

  IF(DEBUG .AND. ANY(ABS(EIG_IM(:)).GT.1.0e-8)) THEN
    WRITE(*,*) "GENRIGID_NORMALMODES> Error: complex eigenvalues returned from metric-tensor Hessian."
    STOP
  ENDIF

 CALL EIGENSORT_VAL_ASC(EIG_REAL,METRICHESS,DEGFREEDOMS,DEGFREEDOMS)

  IF (MWVECTORS .AND. (AMBERT .OR. AMBER12T)) THEN
     DO J1 = 1, DEGFREEDOMS
        IF (EIG_REAL(J1) > 0.0D0) THEN
           EIG_REAL(J1) = SQRT((EIG_REAL(J1)))*108.52
        ELSE
           EIG_REAL(J1) = -SQRT((-EIG_REAL(J1)))*108.52
        ENDIF
     ENDDO
  ENDIF

  ! Copy the eigenvalues to DIAG
  DIAG(1:DEGFREEDOMS) = EIG_REAL(1:DEGFREEDOMS)
  ! Fill in the remaining elements with a very large number so that after sorting, these get ignored.
  DIAG(DEGFREEDOMS+1:3*NATOMS) = 1.0D10

! Set the Hessian matrix to its eigenvector matrix
  HESS(1:DEGFREEDOMS,1:DEGFREEDOMS) = EVECS(1:DEGFREEDOMS,1:DEGFREEDOMS)

  ! Reset this variable to its original value
  ATOMRIGIDCOORDT = ART
  RBAANORMALMODET = RBT

  OPEN(UNIT = 28, FILE = 'LRBNORMALMODES')
  WRITE(28, *) ENERGY, 3*NATOMS, DEGFREEDOMS
  DO J1 = 3, DEGFREEDOMS/3
     WRITE(28, *) DIAG(3*J1-2), DIAG(3*J1-1), DIAG(3*J1)
  ENDDO
  CLOSE(UNIT = 28)

  IF (DUMPV) THEN
     CALL GENRIGID_EIGENVECTORS(METRICHESS, XRIGIDCOORDS, ATOMMASS)
  ENDIF

END SUBROUTINE GENRIGID_NORMALMODES
!----------------------------------------------------------------------------------------------
! See section 2.2.1 of Mochizuki et al, PCCP 16 (2014), 2842-53  DOI 10.1039/C3CP53537A
! However, the notation used is closer to Wales and Ohmine, JCP 98 (1993), 7257-7268 (KBLOCK, U)
! See also rigidb.f90>NRMLMD, an obsolete version of the same procedure which has comments from the
! original authors. The comments on this routine were mostly added by sn402.
SUBROUTINE GENRIGID_EIGENVALUES(X, ATOMMASS, DIAG, INFO)
  
  USE COMMONS
  USE MODHESS
  USE MODCHARMM
  USE KEY
  IMPLICIT NONE

  INTEGER          :: IR, IC, OFFSET, NDIM
  INTEGER          :: I, J, J1, J2, J3, J5, J8, K1, K2, ISTART, IFINISH, JSTART, JFINISH 
  DOUBLE PRECISION :: U(DEGFREEDOMS,DEGFREEDOMS), KBLOCK(3,3), KBEGNV(3), TMASS, ENERGY
  DOUBLE PRECISION :: XRIGIDCOORDS(DEGFREEDOMS), XCOORDS(3*NATOMS), G(3*NATOMS)
  DOUBLE PRECISION :: X(3*NATOMS), FRQN(DEGFREEDOMS), ATOMMASS(3*NATOMS), DIAG(3*NATOMS)
  DOUBLE PRECISION :: P(3), RMI(3,3), DRMI(3,3), DR(3)
  DOUBLE PRECISION :: KD(DEGFREEDOMS), AP(DEGFREEDOMS,DEGFREEDOMS), RMS
! the following required to call the LAPACK routine DSYEV
  INTEGER          :: INFO
  INTEGER, PARAMETER :: LWORK = 1000000 ! the dimension is set arbitrarily
  DOUBLE PRECISION :: WORK(LWORK)
  LOGICAL          :: ART

  ! The METRICTENSOR keyword instructs us to override all calls to this subroutine with a call
  ! to the newer subroutine GENRIGID_NORMALMODES. They should have exactly the same behaviour.
  IF(METRICTENSOR) THEN
     CALL GENRIGID_NORMALMODES(X, ATOMMASS, DIAG, INFO)
     RETURN
  ENDIF

! Initialize
! hk286
  OFFSET = 3*NRIGIDBODY
  U(:,:) = 0.D0
  IR     = 0  ! IR+1 is the index of the first coordinate for the current rigid body
  IC     = 0
  ART = ATOMRIGIDCOORDT

! Transform coordinates
  IF ( ATOMRIGIDCOORDT .EQV. .TRUE. ) THEN
     CALL TRANSFORMCTORIGID (X, XRIGIDCOORDS)
  ELSE
     XRIGIDCOORDS(1:DEGFREEDOMS) = X(1:DEGFREEDOMS)
  ENDIF
  
  DO J1 = 1, NRIGIDBODY

     J3 = 3*J1
     J5 = OFFSET + J3
     P  = XRIGIDCOORDS(J5-2:J5)
     ! KBLOCK is the moment of inertia matrix for this rigid body.
     KBLOCK(:,:) = 0.D0     
     ! Get the rotation matrix RMI that corresponds to P (doesn't calculate any derivatives)
     CALL RMDFAS(P, RMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, .FALSE., .FALSE.)

     TMASS = 0.0D0 ! The total mass of the rigid body
     DO J2 = 1, NSITEPERBODY(J1)
        J8 = RIGIDGROUPS(J2, J1)
        ! DR is the coordinates of this atom relative to the rigid body CoM.
        DR(:)  = MATMUL(RMI(:,:),SITESRIGIDBODY(J2,:,J1))
        TMASS = TMASS + ATOMMASS(J8)
        DO I = 1, 3
           ! Diagonal terms in the moment of inertia about the rigid-body CoM
           ! sn402: Is this right? Shouldn't Ixx = sum_i m_i * (y_i^2+z_i^2) ? Here we have Ixx = sum_i m_i * (x_i^2+y_i^2+z_i^2)
           KBLOCK(I,I) = KBLOCK(I,I) + ATOMMASS(J8)*(DR(1)*DR(1) + DR(2)*DR(2) + DR(3)*DR(3))
           DO J = 1, 3    ! could have been J = 1, I; KBLOCK is a symmetric matrix
              KBLOCK(I,J) = KBLOCK(I,J) - ATOMMASS(J8)*DR(I)*DR(J)
           ENDDO
        ENDDO
     ENDDO
     ! Diagonalise KBLOCK
     ! KBEGNV are the KBLOCK eigenvalues, KBLOCK now contains the eigenvectors
     CALL DSYEV('V','L',3,KBLOCK,3,KBEGNV,WORK,LWORK,INFO)
     IF (INFO /= 0) THEN
        WRITE(*,*) 'NRMLMD > Error in DSYEV with KBLOCK, INFO =', INFO
        STOP
     ENDIF

!    Construction of the matrix U, which diagonalises the coordinate vector x = (r,p) for this rigid body.
!    The reason we want this is that making the transform w = matmul(U,p) gives us coordinates
!    which are diagonal in the kinetic energy.

!    First: translation coordinates (recall IR+1 is the first coordinate for this rigid body)
!    No diagonalisation required. U is identity.
     U(IR+1,IC+1) = 1.D0; U(IR+2,IC+2) = 1.D0; U(IR+3,IC+3) = 1.D0
     ! KD is the diagonalised KBLOCK. Only diagonal elements are recorded.
     ! Scaling by sqrt(TMASS) is required to move to the CoM frame, as usual in a normal mode calculation.
     KD(IC+1:IC+3) = 1.D0/SQRT(TMASS)

!    Now rotational coordinates.
!    Recall that KBLOCK now contains the eigenvectors of the original matrix. These are required to diagonalise
!    the rotational component of the coordinate vector.
     U(OFFSET+IR+1:OFFSET+IR+3,OFFSET+IC+1:OFFSET+IC+3) = KBLOCK(:,:)
     KD(OFFSET+IC+1:OFFSET+IC+3) = 1.D0/SQRT(KBEGNV(:))               ! See Mochizuki et al

     ! IR = IC = 3*(J1-1) - these two variables indicate the first coordinate position of the current rigid body.
     IR = IR + 3
     IC = IC + 3 

  ENDDO ! Loop over rigid bodies (J1)

  ! Deal with the free atoms. No off-diagonal terms in the kinetic energy here, so they are automatically diagonalised and only
  ! require mass-rescaling.
  DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
     J8 = RIGIDSINGLES(J1)
     U(6*NRIGIDBODY+3*J1-2,6*NRIGIDBODY+3*J1-2) = 1.D0
     U(6*NRIGIDBODY+3*J1-1,6*NRIGIDBODY+3*J1-1) = 1.D0
     U(6*NRIGIDBODY+3*J1  ,6*NRIGIDBODY+3*J1  ) = 1.D0
     KD(6*NRIGIDBODY+3*J1-2:6*NRIGIDBODY+3*J1 ) = 1.D0/SQRT(ATOMMASS(J8))     
  ENDDO

  ! We have now obtained the required diagonal coordinates for the total kinetic energy.
  ! We next need to obtain the Hessian in terms of these coordinates.

  RBAANORMALMODET = .TRUE.
  ATOMRIGIDCOORDT = .FALSE.
  XCOORDS(1:DEGFREEDOMS) = XRIGIDCOORDS(1:DEGFREEDOMS)
  XCOORDS(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
  IF (ENDNUMHESS .OR. AMBERT .OR. AMBER12T) THEN
     CALL GENRIGID_MAKENUMHESS(XCOORDS,NATOMS,DEGFREEDOMS)
  ELSE
     ! STEST is set to TRUE so we get the Hessian calculated.
     ! When TRANSFORM_HESSIAN is called from POTENTIAL, it will behave differently because of RBAANORMALMODET=TRUE.
     ! See TRANSFORM_HESSIAN for comments.
     CALL POTENTIAL(XCOORDS,ENERGY,G,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)     
  ENDIF
  RBAANORMALMODET = .FALSE.
!  ATOMRIGIDCOORDT = .TRUE.

  NDIM = DEGFREEDOMS  
  AP(:,:) = 0.D0  ! This will be our Hessian in diagonalised general coordinates.
  ! Fill in the lower-left triangle of the transformed Hessian
  DO I = 1, NDIM
     DO J = 1, I
        
        ! sn402: I'm not quite sure what's going on here.
        ! This block seems to suggest that the block AP(1:3,1:3) gets written to more than everything else.
        ! Wouldn't it be easier to use the upper-right triangle and use
        !DO I=1,NDIM
        !   DO J=I,NDIM
        !      DO K1=1,NDIM
        !         DO K2=1,NDIM
        !            ?
        ! Possibly this way is more efficient because U is 0 away from the main diagonal?
        IF (I .LE. 2) THEN
           ISTART = 1
        ELSE
           ISTART = I-2
        ENDIF
        IF (I .GE. NDIM -2) THEN
           IFINISH = NDIM
        ELSE
           IFINISH = I+2
        ENDIF
        IF (J .LE. 2) THEN
           JSTART = 1
        ELSE
           JSTART = J-2
        ENDIF
        IF (J .GE. NDIM-2) THEN
           JFINISH = NDIM
        ELSE
           JFINISH = J+2
        ENDIF

        DO K1 = ISTART, IFINISH
           DO K2 = JSTART, JFINISH
              ! The actual coordinate transform. Recall that U is diagonal for translational coordinates and single
              ! atoms (the elements are square roots of the corresponding masses)
              ! But U and HESS are both different for the rotational coordinates.
              AP(I,J) = AP(I,J) + U(K1,I)*HESS(K1,K2)*U(K2,J)
           ENDDO
        ENDDO
        AP(I,J) = KD(I)*AP(I,J)*KD(J) ! Weight the coordinates by the eigenvalues of the inertia tensor (see above)
     ENDDO
  ENDDO

  ! Now diagonalise the modified Hessian. The eigenvalues (returned as FRQN) are the normal mode frequencies.
  ! If DUMPV is set, then on return AP contains the eigenvectors of the transformed Hessian - i.e. the normal modes.
  IF (DUMPV) THEN
     CALL DSYEV('V','L',NDIM,AP(1:NDIM,1:NDIM),NDIM,FRQN,WORK,LWORK,INFO)
  ELSE
     CALL DSYEV('N','L',NDIM,AP(1:NDIM,1:NDIM),NDIM,FRQN,WORK,LWORK,INFO)
  ENDIF

!  call eigensort_val_asc(FRQN,AP,NDIM,NDIM)
! Mass-weight the normal mode frequencies here?
! Looks like someone has hard-coded the units in. Need to change that.
  IF (MWVECTORS .AND. (AMBERT .OR. AMBER12T)) THEN ! .OR. CHARMMT)) THEN  ! I need to get this included at some point
     ! Hessian diagonalisation returns square frequencies in internal units. For AMBER/CHARMM, these units are
     ! (kCal mol^-1)/(amu Angstrom^2) = 4.184E26 s^-2.
     ! So to convert to rad s^-1, we want sqrt(FRQN(I)*4.184D26) = 2.045E13*sqrt(FRQN(I))
     ! For Hz: sqrt(FRQN(I)*4.184D26)/2*pi = 3.255E12*sqrt(FRQN(I))
     ! For cm^-1: sqrt(FRQN(I)*4.184D26)/(2*pi*c*100) = 108.52*sqrt(FRQN(I))
     DO I = 1, NDIM
        IF (FRQN(I) > 0.0D0) THEN
           FRQN(I) = SQRT((FRQN(I)))*108.52
        ELSE
           FRQN(I) = -SQRT((-FRQN(I)))*108.52
        ENDIF
     ENDDO
  ENDIF

  DIAG(1:DEGFREEDOMS) = FRQN
  DIAG(DEGFREEDOMS+1:3*NATOMS) = 1.0D10

! hk286 - set the Hessian matrix to its eigenvector matrix
  HESS(1:DEGFREEDOMS,1:DEGFREEDOMS) = AP
! Restore this variable to its saved value.
  ATOMRIGIDCOORDT = ART

!  OPEN(UNIT = 28, FILE = 'LRBNORMALMODES')
!  WRITE(28, *) ENERGY, 3*NATOMS, DEGFREEDOMS
!  DO J1 = 3, DEGFREEDOMS/3
!     WRITE(28, *) FRQN(3*J1-2), FRQN(3*J1-1), FRQN(3*J1)
!  ENDDO
!  CLOSE(UNIT = 28)

  ! hk286
  IF (DUMPV) THEN      
     CALL GENRIGID_EIGENVECTORS(AP, XRIGIDCOORDS, ATOMMASS)
  ENDIF
  
END SUBROUTINE GENRIGID_EIGENVALUES

!----------------------------------------------------------------------------------------------

SUBROUTINE GENRIGID_EIGENVECTORS(AP, XRIGIDCOORDS, ATOMMASS)

  USE COMMONS
  USE MODHESS
  IMPLICIT NONE
        
  INTEGER :: OFFSET, SMODE, I, J, J1, J2, J3, J5, J8, J9
  DOUBLE PRECISION :: XTEMP(DEGFREEDOMS), XXTEMP(DEGFREEDOMS), AP(DEGFREEDOMS,DEGFREEDOMS)
  DOUBLE PRECISION :: X(3*NATOMS), XCOORDS(3*NATOMS), XRIGIDCOORDS(3*NATOMS)
  DOUBLE PRECISION :: P(3), KBLOCK(3,3), KBEGNV(3), RMI(3,3), RRMI(3,3), DRMI(3,3), DR(3)
  DOUBLE PRECISION :: TMASS, ATOMMASS(3*NATOMS)
  INTEGER, PARAMETER :: LWORK = 10000 ! the dimension is set arbitrarily
  DOUBLE PRECISION :: WORK(LWORK)
  INTEGER          :: INFO, TFRAME

  OFFSET = 3*NRIGIDBODY
  CALL TRANSFORMRIGIDTOC (1, NRIGIDBODY, X, XRIGIDCOORDS)

! Here are the normal modes you want to be computed - the first six are zero
  DO SMODE = 1, DEGFREEDOMS
     DO I = 1, DEGFREEDOMS
! This gives you DX where DX is taken along the direction of an eigenvector
        XTEMP(I) = AP(I,SMODE) * 0.05
     ENDDO
     DO J1 = 1, NRIGIDBODY
        J3 = 3*J1
! The rotational part - undiagonalized it
        J5 = OFFSET + J3
        P  = XRIGIDCOORDS(J5-2:J5)
        KBLOCK(:,:) = 0.D0                 
! Computes the rotation matrix which brings the reference geometry in stationary frame 
! to the atom positions in current minimum geometry
        CALL RMDFAS(P, RMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, .FALSE., .FALSE.)
                 
! Computing inertia matrix in the moving frame, i.e. current minimum geometry
        TMASS = 0.0D0
        DO J2 = 1, NSITEPERBODY(J1)
           J8 = RIGIDGROUPS(J2, J1)        
           DR(:)  = MATMUL(RMI(:,:),SITESRIGIDBODY(J2,:,J1))
           TMASS = TMASS + ATOMMASS(J8)
           DO I = 1, 3
              KBLOCK(I,I) = KBLOCK(I,I) + ATOMMASS(J8)*(DR(1)*DR(1) + DR(2)*DR(2) + DR(3)*DR(3))
              DO J = 1, 3    ! could have been J = 1, I; KBLOCK is a symmetric matrix
                 KBLOCK(I,J) = KBLOCK(I,J) - ATOMMASS(J8)*DR(I)*DR(J)
              ENDDO
           ENDDO
        ENDDO
! Diagonalise inertia matrix
        CALL DSYEV('V','L',3,KBLOCK,3,KBEGNV,WORK,LWORK,INFO)

        XTEMP(J5-2) = XTEMP(J5-2)/SQRT(KBEGNV(1)) 
        XTEMP(J5-1) = XTEMP(J5-1)/SQRT(KBEGNV(2)) 
        XTEMP(J5  ) = XTEMP(J5  )/SQRT(KBEGNV(3)) 
                 
! Going from the diagonalised rotation coordinates to per rigid body angle-axis coordinates in the moving frame
        XXTEMP(J5-2) = KBLOCK(1,1)*XTEMP(J5-2) + KBLOCK(1,2)*XTEMP(J5-1) + KBLOCK(1,3)*XTEMP(J5  )
        XXTEMP(J5-1) = KBLOCK(2,1)*XTEMP(J5-2) + KBLOCK(2,2)*XTEMP(J5-1) + KBLOCK(2,3)*XTEMP(J5  )
        XXTEMP(J5  ) = KBLOCK(3,1)*XTEMP(J5-2) + KBLOCK(3,2)*XTEMP(J5-1) + KBLOCK(3,3)*XTEMP(J5  )
                 
! Computing the rotation matrix which rotates the current minimum geometry due to its normal mode displacements
        P = XXTEMP(J5-2:J5)
        CALL RMDFAS(P, RRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, .FALSE., .FALSE.)

! Compute the displaced positions of the atoms
        DO J2 = 1, NSITEPERBODY(J1)
           J8 = RIGIDGROUPS(J2, J1)        
           XCOORDS(3*J8-2:3*J8) = XRIGIDCOORDS(J3-2:J3) + XTEMP(J3-2:J3)/SQRT(TMASS) &
                + MATMUL(RRMI(:,:),MATMUL(RMI(:,:),SITESRIGIDBODY(J2,:,J1)))
        ENDDO
     ENDDO

     DO J1 = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
        J8 = RIGIDSINGLES(J1)
        J5 = 2*OFFSET + 3*J1
        XCOORDS(3*J8-2:3*J8) = XRIGIDCOORDS(J5-2:J5) + XTEMP(J5-2:J5)/SQRT(ATOMMASS(J8))
     ENDDO

! Computes the atoms eigenvectors
     DO I = 1, 3*NATOMS
        HESS(I,SMODE) = (XCOORDS(I)-X(I)) / 0.05
     ENDDO

  ENDDO

END SUBROUTINE GENRIGID_EIGENVECTORS

!----------------------------------------------------------------------------------------------

SUBROUTINE GENRIGID_VDUMP(DIAG,ZT,N,M)
  USE KEY
  USE MODHESS
  USE MODCHARMM, ONLY: CHRMMT
  USE PORFUNCS
  IMPLICIT NONE
  INTEGER M, N, J1, J2, ISTAT, MCOUNT
  DOUBLE PRECISION DIAG(M)
  LOGICAL ZT(M)

!
!  dump the eigenvectors which correspond to non-zero eigenvalues
!  in file vectors.dump
!
  IF (.NOT.ALLSTEPS) REWIND(44)
  IF (ALLVECTORS) THEN
     MCOUNT=0
     DO J1=1,M
        IF (ZT(J1)) MCOUNT=MCOUNT+1
     ENDDO
     OPEN(UNIT=499,FILE='nmodes.dat',STATUS='UNKNOWN')
     WRITE(499,'(I6)') MCOUNT
     CLOSE(499)
     DO J1=1,M
        IF (ZT(J1)) THEN
! If printing the mass weighted vectors (normal modes), convert omega^2
! into the vibrational frequency in the specified unit system using FRQCONV.
! Normally, this will be either internal units or rad/s, for compatibility with PATHSAMPLE.
! Other unit systems can be specified using the FRQCONV keyword.
! NOTE: This behaviour has changed as of 23/9/16. Until now, the frequencies were always
! multiplied by 108.52, which is the conversion factor from kCal mol^-1 and Angstrom units
! to cm^-1 frequencies. However, this is obviously not appropriate for all systems. If
! you wish to retrieve the old behaviour, simply add FRQCONV 108.52 to your odata file.
! But note that the square frequencies used for the log product in path.info will then be given
! in units of cm^-2, rather than the default which is s^-1 for AMBER and CHARMM, internal units
! for most other potentials.
           IF (MWVECTORS) THEN
              WRITE(44,'(F20.10)') DSQRT(DIAG(J1))*FRQCONV
           ELSE
              WRITE(44,'(F20.10)') DIAG(J1)
           ENDIF
           WRITE(44,'(3F20.10)') (HESS(J2,J1),J2=1,N)
        ENDIF
     ENDDO
  ELSE
     DO J1=M,1,-1
        IF (ZT(J1)) THEN
! As above
           IF (MWVECTORS) THEN
              WRITE(44,'(F20.10)') DSQRT(DIAG(J1))*FRQCONV
           ELSE
              WRITE(44,'(F20.10)') DIAG(J1)
           ENDIF
           WRITE(44,'(3F20.10)') (HESS(J2,J1),J2=1,N)
           CALL FLUSH(44)
           RETURN
        ENDIF
     ENDDO
  ENDIF
  CALL FLUSH(44)
  RETURN

END SUBROUTINE GENRIGID_VDUMP

!----------------------------------------------------------------------------------------------

SUBROUTINE GENRIGID_MAKENUMHESS(X,NATOMS,DEGFREEDOMS)
  
  USE MODHESS
  USE MODCHARMM
  USE COMMONS,ONLY : DEBUG
  use porfuncs
  IMPLICIT NONE
  LOGICAL KNOWE, KNOWG, KNOWH
  COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
  
  INTEGER I1,J1,NATOMS,DEGFREEDOMS,ISTAT
  DOUBLE PRECISION X(3*NATOMS)
  DOUBLE PRECISION DUM(3*NATOMS),GRAD1(3*NATOMS),GRAD2(3*NATOMS),DELTA,RMS,ENERGY
  
  IF (DEBUG) WRITE(*,'(A)') ' makenumhess> Making numerical hessian for local rigid body'
  DO I1=1,DEGFREEDOMS
     DUM(I1)=X(I1)
  ENDDO
  
  DELTA=1.0D-6
  IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
  
  IF (.NOT.ALLOCATED(HESS)) ALLOCATE(HESS(3*NATOMS,3*NATOMS))
  DO I1=1,DEGFREEDOMS
     DUM(I1)=X(I1)-DELTA
     CALL POTENTIAL(DUM,ENERGY,GRAD1,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
     DUM(I1)=X(I1)+DELTA
     CALL POTENTIAL(DUM,ENERGY,GRAD2,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
     DUM(I1)=X(I1)
     DO J1=I1,DEGFREEDOMS
        HESS(I1,J1)=(GRAD2(J1)-GRAD1(J1))/(2.0D0*DELTA)
        HESS(J1,I1)=HESS(I1,J1)
     ENDDO
     
  ENDDO
  
  HESS(DEGFREEDOMS+1:3*NATOMS,:) = 0.0D0
  HESS(:,DEGFREEDOMS+1:3*NATOMS) = 0.0D0

  IF (DEBUG) WRITE(*,'(A)') ' makenumhess> Hessian made for local rigid body'
  KNOWH=.TRUE.

  RETURN
END SUBROUTINE GENRIGID_MAKENUMHESS

!----------------------------------------------------------------------------------------------
! hk286 - this is still not working
SUBROUTINE SHIFTLOCALRIGID (Q, NATOMS)

!     THIS SUBROUTINE SHIFTS THE 'ZERO' EIGENVALUES CORRESPONDING TO OVERALL TRANSLATION AND
!     ROTATION OF A SYSTEM OF (IDENTICAL) RIGID BODIES WITH C-O-M & ANGLE-AXIS COORDINATES.

  USE KEY
  USE MODHESS

  IMPLICIT NONE
  
  INTEGER            :: NATOMS, I, J, J1, J2
  DOUBLE PRECISION   :: Q(3*NATOMS), EV(3*NATOMS,6), NRMFCT(6)
  DOUBLE PRECISION   :: CMX, CMY, CMZ, THETA, THETA2, THETAH, DUMMY
  
!     INITIALIZE

  EV(:,:)   = 0.D0
  NRMFCT(:) = 0.D0
  CMX       = 0.D0
  CMY       = 0.D0
  CMZ       = 0.D0
  
  SHIFTED = .TRUE.
  NZERO = 6

  DO I = 1, NRIGIDBODY
     J = 3*I     
     CMX = CMX + Q(J-2) !* NSITEPERBODY(I)
     CMY = CMY + Q(J-1) !* NSITEPERBODY(I)
     CMZ = CMZ + Q(J)   !* NSITEPERBODY(I)     
  ENDDO
  DO I = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
     J = 6*NRIGIDBODY + 3*I
     CMX = CMX + Q(J-2)
     CMY = CMY + Q(J-1)
     CMZ = CMZ + Q(J)
  ENDDO
  CMX = CMX / FLOAT(NATOMS)
  CMY = CMY / FLOAT(NATOMS)
  CMZ = CMZ / FLOAT(NATOMS)

  CMX = 0.0D0
  CMY = 0.0D0
  CMZ = 0.0D0

  DO I = 1, NRIGIDBODY
     J  = 3*I
     J1 = 3*NRIGIDBODY + J

     THETA2 = DOT_PRODUCT(Q(J1-2:J1), Q(J1-2:J1))
     THETA  = DSQRT(THETA2)
     THETAH = 0.5D0*THETA

!     TRANSLATION ALONG X
     EV(J-2,1) = 1.D0
     NRMFCT(1) = NRMFCT(1) + 1.D0

!     TRANSLATION ALONG Y
     EV(J-1,2) = 1.D0
     NRMFCT(2) = NRMFCT(2) + 1.D0

!     TRANSLATION ALONG Z
     EV(J,3) = 1.D0
     NRMFCT(3) = NRMFCT(3) + 1.D0

     IF (THETA == 0.D0) THEN

!     ROTATION ABOUT Z
        EV(J-2,4)  = - Q(J-1) + CMY
        EV(J-1,4)  = Q(J-2) - CMX
        EV(J1,4)   = 1.D0
        NRMFCT(4)  = NRMFCT(4) + EV(J-2,4)**2 + EV(J-1,4)**2 + EV(J1-2,4)**2 + EV(J1-1,4)**2 + EV(J1,4)**2
        
!     ROTATION ABOUT X
        EV(J-1,5)  = - Q(J) + CMZ
        EV(J,5)    = Q(J-1) - CMY
        EV(J1-2,5) = 1.D0
        NRMFCT(5)  = NRMFCT(5) + EV(J-1,5)**2 + EV(J,5)**2 + EV(J1-2,5)**2 + EV(J1-1,5)**2 + EV(J1,5)**2
            
!     ROTATION ABOUT Y
        EV(J-2,6)  = Q(J) - CMZ
        EV(J,6)    = - Q(J-2) + CMX
        EV(J1-1,6) = 1.D0
        NRMFCT(6)  = NRMFCT(6) + EV(J-2,6)**2 + EV(J,6)**2 + EV(J1-2,6)**2 + EV(J1-1,6)**2 + EV(J1,6)**2
        
     ELSE

!     ROTATION ABOUT Z
        EV(J-2,4)  = - Q(J-1) + CMY
        EV(J-1,4)  = Q(J-2) - CMX
        EV(J1-2,4) = - 0.5D0*Q(J1-1) + Q(J1)*Q(J1-2)/THETA2 - 0.5D0*Q(J1)*Q(J1-2)/(THETA*TAN(THETAH))
        EV(J1-1,4) = 0.5D0*Q(J1-2) + Q(J1)*Q(J1-1)/THETA2 - 0.5D0*Q(J1)*Q(J1-1)/(THETA*TAN(THETAH))
        EV(J1,4)   = THETAH/TAN(THETAH) + Q(J1)**2/THETA2 - 0.5D0*Q(J1)**2/(THETA*TAN(THETAH))
        NRMFCT(4)  = NRMFCT(4) + EV(J-2,4)**2 + EV(J-1,4)**2 + EV(J1-2,4)**2 + EV(J1-1,4)**2 + EV(J1,4)**2
             
!     ROTATION ABOUT X
        EV(J-1,5)  = - Q(J) + CMZ
        EV(J,5)    = Q(J-1) - CMY
        EV(J1-2,5) = THETAH/TAN(THETAH) + Q(J1-2)**2/THETA2 - 0.5D0*Q(J1-2)**2/(THETA*TAN(THETAH))
        EV(J1-1,5) = - 0.5D0*Q(J1) + Q(J1-2)*Q(J1-1)/THETA2 - 0.5D0*Q(J1-2)*Q(J1-1)/(THETA*TAN(THETAH))
        EV(J1,5)   = 0.5D0*Q(J1-1) + Q(J1-2)*Q(J1)/THETA2 - 0.5D0*Q(J1-2)*Q(J1)/(THETA*TAN(THETAH))
        NRMFCT(5)  = NRMFCT(5) + EV(J-1,5)**2 + EV(J,5)**2 + EV(J1-2,5)**2 + EV(J1-1,5)**2 + EV(J1,5)**2
     
!     ROTATION ABOUT Y
        EV(J-2,6)  = Q(J) - CMZ
        EV(J,6)    = - Q(J-2) + CMX
        EV(J1-2,6) = 0.5D0*Q(J1) + Q(J1-1)*Q(J1-2)/THETA2 - 0.5D0*Q(J1-1)*Q(J1-2)/(THETA*TAN(THETAH))
        EV(J1-1,6) = THETAH/TAN(THETAH) + Q(J1-1)**2/THETA2 - 0.5D0*Q(J1-1)**2/(THETA*TAN(THETAH))
        EV(J1,6)   = - 0.5D0*Q(J1-2) + Q(J1-1)*Q(J1)/THETA2 - 0.5D0*Q(J1-1)*Q(J1)/(THETA*TAN(THETAH))
        NRMFCT(6)  = NRMFCT(6) + EV(J-2,6)**2 + EV(J,6)**2 + EV(J1-2,6)**2 + EV(J1-1,6)**2 + EV(J1,6)**2
     
     ENDIF     

  ENDDO

  DO I = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
     J = 6*NRIGIDBODY + 3*I
!     TRANSLATION ALONG X
     EV(J-2,1) = 1.D0
     NRMFCT(1) = NRMFCT(1) + 1.D0
!     TRANSLATION ALONG Y
     EV(J-1,2) = 1.D0
     NRMFCT(2) = NRMFCT(2) + 1.D0
!     TRANSLATION ALONG Z
     EV(J,3) = 1.D0
     NRMFCT(3) = NRMFCT(3) + 1.D0
!     ROTATION ABOUT Z
     EV(J-2,4)  = - Q(J-1) + CMY
     EV(J-1,4)  = Q(J-2) - CMX
     NRMFCT(4)  = NRMFCT(4) + EV(J-2,4)**2 + EV(J-1,4)**2        
!     ROTATION ABOUT X
     EV(J-1,5)  = - Q(J) + CMZ
     EV(J,5)    = Q(J-1) - CMY
     NRMFCT(5)  = NRMFCT(5) + EV(J-1,5)**2 + EV(J,5)**2            
!     ROTATION ABOUT Y
     EV(J-2,6)  = Q(J) - CMZ
     EV(J,6)    = - Q(J-2) + CMX
     NRMFCT(6)  = NRMFCT(6) + EV(J-2,6)**2 + EV(J,6)**2        
  ENDDO

  DO J = 1, NZERO

     NRMFCT(J) = DSQRT(NRMFCT(J))
     EV(:,J)   = EV(:,J)/NRMFCT(J)
   
  ENDDO

!     GRAM-SCHMIDT ORTHOGONALISATION TO OBTAIN ORTHONORMAL ROTATIONAL EIGENVECTORS

  DO J = 4, NZERO
     
     DO J1 = 4, J-1
        
        EV(:,J) = EV(:,J) - DOT_PRODUCT(EV(:,J),EV(:,J1))*EV(:,J1)
        
     ENDDO
     
     EV(:,J) = EV(:,J) / DSQRT(DOT_PRODUCT(EV(:,J),EV(:,J)))
     
  ENDDO

  DO J1 = 1, DEGFREEDOMS
     
     DO J2 = 1, DEGFREEDOMS
        
        DO J = 1, NZERO 
           
           HESS(J2,J1) = HESS(J2,J1) + SHIFTV*EV(J2,J)*EV(J1,J)
           
        ENDDO
        
     ENDDO
     
  ENDDO

END SUBROUTINE SHIFTLOCALRIGID

!----------------------------------------------------------------------------------------------

SUBROUTINE ORTHOGLOCALRIGID (VEC1, Q, OTEST)

!     THIS SUBROUTINE ORTHOGONALISES VEC1 TO THE EIGENVECTORS CORRESPONDING TO OVERALL TRANSLATION
!     AND ROTATION OF A SYSTEM OF (IDENTICAL) RIGID BODIES WITH C-O-M & ANGLE-AXIS COORDINATES.

  USE COMMONS
  USE KEY
  
  IMPLICIT NONE
  
! hk286 - use 3*NATOMS or DEGFREEDOMS
  INTEGER            :: I, J, J1
  DOUBLE PRECISION   :: Q(3*NATOMS), EV(3*NATOMS,6), NRMFCT(6), VEC1(3*NATOMS)
  DOUBLE PRECISION   :: CMX, CMY, CMZ, THETA, THETA2, THETAH, DUMMY
  LOGICAL            :: OTEST

  !     INITIALIZE

  EV(:,:)   = 0.D0
  NRMFCT(:) = 0.D0
  CMX       = 0.D0
  CMY       = 0.D0
  CMZ       = 0.D0
  
  NZERO = 6

  CMX = 0.0D0
  CMY = 0.0D0
  CMZ = 0.0D0

  DO I = 1, NRIGIDBODY
     J  = 3*I
     J1 = 3*NRIGIDBODY + J

     THETA2 = DOT_PRODUCT(Q(J1-2:J1), Q(J1-2:J1))
     THETA  = DSQRT(THETA2)
     THETAH = 0.5D0*THETA

!     TRANSLATION ALONG X
     EV(J-2,1) = 1.D0
     NRMFCT(1) = NRMFCT(1) + 1.D0

!     TRANSLATION ALONG Y
     EV(J-1,2) = 1.D0
     NRMFCT(2) = NRMFCT(2) + 1.D0

!     TRANSLATION ALONG Z
     EV(J,3) = 1.D0
     NRMFCT(3) = NRMFCT(3) + 1.D0

     IF (THETA == 0.D0) THEN

!     ROTATION ABOUT Z
        EV(J-2,4)  = - Q(J-1) + CMY
        EV(J-1,4)  = Q(J-2) - CMX
        EV(J1,4)   = 1.D0
        NRMFCT(4)  = NRMFCT(4) + EV(J-2,4)**2 + EV(J-1,4)**2 + EV(J1-2,4)**2 + EV(J1-1,4)**2 + EV(J1,4)**2
        
!     ROTATION ABOUT X
        EV(J-1,5)  = - Q(J) + CMZ
        EV(J,5)    = Q(J-1) - CMY
        EV(J1-2,5) = 1.D0
        NRMFCT(5)  = NRMFCT(5) + EV(J-1,5)**2 + EV(J,5)**2 + EV(J1-2,5)**2 + EV(J1-1,5)**2 + EV(J1,5)**2
            
!     ROTATION ABOUT Y
        EV(J-2,6)  = Q(J) - CMZ
        EV(J,6)    = - Q(J-2) + CMX
        EV(J1-1,6) = 1.D0
        NRMFCT(6)  = NRMFCT(6) + EV(J-2,6)**2 + EV(J,6)**2 + EV(J1-2,6)**2 + EV(J1-1,6)**2 + EV(J1,6)**2
        
     ELSE

!     ROTATION ABOUT Z
        EV(J-2,4)  = - Q(J-1) + CMY
        EV(J-1,4)  = Q(J-2) - CMX
        EV(J1-2,4) = - 0.5D0*Q(J1-1) + Q(J1)*Q(J1-2)/THETA2 - 0.5D0*Q(J1)*Q(J1-2)/(THETA*TAN(THETAH))
        EV(J1-1,4) = 0.5D0*Q(J1-2) + Q(J1)*Q(J1-1)/THETA2 - 0.5D0*Q(J1)*Q(J1-1)/(THETA*TAN(THETAH))
        EV(J1,4)   = THETAH/TAN(THETAH) + Q(J1)**2/THETA2 - 0.5D0*Q(J1)**2/(THETA*TAN(THETAH))
        NRMFCT(4)  = NRMFCT(4) + EV(J-2,4)**2 + EV(J-1,4)**2 + EV(J1-2,4)**2 + EV(J1-1,4)**2 + EV(J1,4)**2
             
!     ROTATION ABOUT X
        EV(J-1,5)  = - Q(J) + CMZ
        EV(J,5)    = Q(J-1) - CMY
        EV(J1-2,5) = THETAH/TAN(THETAH) + Q(J1-2)**2/THETA2 - 0.5D0*Q(J1-2)**2/(THETA*TAN(THETAH))
        EV(J1-1,5) = - 0.5D0*Q(J1) + Q(J1-2)*Q(J1-1)/THETA2 - 0.5D0*Q(J1-2)*Q(J1-1)/(THETA*TAN(THETAH))
        EV(J1,5)   = 0.5D0*Q(J1-1) + Q(J1-2)*Q(J1)/THETA2 - 0.5D0*Q(J1-2)*Q(J1)/(THETA*TAN(THETAH))
        NRMFCT(5)  = NRMFCT(5) + EV(J-1,5)**2 + EV(J,5)**2 + EV(J1-2,5)**2 + EV(J1-1,5)**2 + EV(J1,5)**2
     
!     ROTATION ABOUT Y
        EV(J-2,6)  = Q(J) - CMZ
        EV(J,6)    = - Q(J-2) + CMX
        EV(J1-2,6) = 0.5D0*Q(J1) + Q(J1-1)*Q(J1-2)/THETA2 - 0.5D0*Q(J1-1)*Q(J1-2)/(THETA*TAN(THETAH))
        EV(J1-1,6) = THETAH/TAN(THETAH) + Q(J1-1)**2/THETA2 - 0.5D0*Q(J1-1)**2/(THETA*TAN(THETAH))
        EV(J1,6)   = - 0.5D0*Q(J1-2) + Q(J1-1)*Q(J1)/THETA2 - 0.5D0*Q(J1-1)*Q(J1)/(THETA*TAN(THETAH))
        NRMFCT(6)  = NRMFCT(6) + EV(J-2,6)**2 + EV(J,6)**2 + EV(J1-2,6)**2 + EV(J1-1,6)**2 + EV(J1,6)**2
     
     ENDIF     

  ENDDO
  
  DO I = 1, (DEGFREEDOMS - 6*NRIGIDBODY)/3
     J = 6*NRIGIDBODY + 3*I
!     TRANSLATION ALONG X
     EV(J-2,1) = 1.D0
     NRMFCT(1) = NRMFCT(1) + 1.D0
!     TRANSLATION ALONG Y
     EV(J-1,2) = 1.D0
     NRMFCT(2) = NRMFCT(2) + 1.D0
!     TRANSLATION ALONG Z
     EV(J,3) = 1.D0
     NRMFCT(3) = NRMFCT(3) + 1.D0
!     ROTATION ABOUT Z
     EV(J-2,4)  = - Q(J-1) + CMY
     EV(J-1,4)  = Q(J-2) - CMX
     NRMFCT(4)  = NRMFCT(4) + EV(J-2,4)**2 + EV(J-1,4)**2        
!     ROTATION ABOUT X
     EV(J-1,5)  = - Q(J) + CMZ
     EV(J,5)    = Q(J-1) - CMY
     NRMFCT(5)  = NRMFCT(5) + EV(J-1,5)**2 + EV(J,5)**2            
!     ROTATION ABOUT Y
     EV(J-2,6)  = Q(J) - CMZ
     EV(J,6)    = - Q(J-2) + CMX
     NRMFCT(6)  = NRMFCT(6) + EV(J-2,6)**2 + EV(J,6)**2        
  ENDDO

  DO J = 1, NZERO

     NRMFCT(J) = DSQRT(NRMFCT(J))
     EV(:,J)   = EV(:,J)/NRMFCT(J)
   
  ENDDO


!     GRAM-SCHMIDT ORTHOGONALISATION TO OBTAIN ORTHONORMAL ROTATIONAL EIGENVECTORS

  DO J = 4, NZERO
     
     DO J1 = 4, J-1
        
        EV(:,J) = EV(:,J) - DOT_PRODUCT(EV(:,J),EV(:,J1))*EV(:,J1)
        
     ENDDO
     
     EV(:,J) = EV(:,J) / DSQRT(DOT_PRODUCT(EV(:,J),EV(:,J)))
     
  ENDDO

!     PROJECT TRANS/ROT SET OUT OF VEC1

  DO J = 1, NZERO 

     DUMMY   = DOT_PRODUCT(VEC1(:),EV(:,J))
     VEC1(:) = VEC1(:) - DUMMY*EV(:,J)

  ENDDO

  IF (OTEST) CALL VECNORM(VEC1,DEGFREEDOMS) ! NORMALIZE VEC1

END SUBROUTINE ORTHOGLOCALRIGID


SUBROUTINE BEIGLOCALRIGID(ITER,COORDS,ENERGY,VECS,EVALMIN,NS,SOVER,PTEST,CONVERGED)

  USE COMMONS
  USE KEY
  USE MODNEB
  USE MODTWOEND
  use porfuncs
  IMPLICIT NONE

  INTEGER, INTENT(OUT) ::            NS ! the number of iterations spent in LBFGS
  INTEGER, INTENT(IN) ::             ITER ! the number of iterations of Hybrid eigenvector following we have already done
  DOUBLE PRECISION, INTENT(OUT) ::   EVALMIN ! the eigenvalue
  DOUBLE PRECISION, INTENT(OUT) ::   SOVER ! the overlap with the previous eigenvector
  LOGICAL, INTENT(IN) ::             PTEST ! Flag for printing status info
  DOUBLE PRECISION, INTENT(IN) ::    COORDS(3*NATOMS) ! the coordinates at which the eigenvector is being calculated
  DOUBLE PRECISION, INTENT(INOUT) :: ENERGY ! the energy at COORDS
  INTEGER J1, ISEED, NITS
  DOUBLE PRECISION VEC(3*NATOMS),VECS(3*NATOMS),DPRAND,EPS,FRET
  DOUBLE PRECISION DIAG(3*NATOMS),WORK(3*NATOMS*(2*XMUPDATE+1)+2*XMUPDATE)
  PARAMETER (EPS=1.D-6)
  LOGICAL CONVERGED
  COMMON /IS/ ISEED
  SAVE
  
  IF ((ITER.EQ.1).AND.(.NOT.READV).AND.(.NOT.TWOENDS).AND.(.NOT.TTDONE).AND.(.NOT.(NEWCONNECTT.OR.NEWNEBT))) THEN
     DO J1=1,DEGFREEDOMS
        VEC(J1)=2*(DPRAND()-0.5D0)
     ENDDO
     DO J1=DEGFREEDOMS+1,NOPT
        VEC(J1)=0.0D0
     ENDDO
  ELSE
     DO J1=1,NOPT
        VEC(J1)=VECS(J1)
     ENDDO
  ENDIF

  CALL XMYLBFGS(NOPT,XMUPDATE,VEC,.FALSE.,DIAG,CEIG,WORK,ENERGY,COORDS,NITS,NEVS,FRET,PTEST,CONVERGED)
  CALL VECNORM(VEC,NOPT) 

  EVALMIN=FRET
  NS=NITS
  SOVER=0.0D0
  DO J1=1,NOPT
     SOVER=SOVER+VECS(J1)*VEC(J1)
     VECS(J1)=VEC(J1)
  ENDDO
  IF (PTEST) WRITE(*,'(A,F15.7)') 'beig> Overlap with previous vector=',SOVER

  FIXIMAGE=.FALSE.
  
  RETURN
END SUBROUTINE BEIGLOCALRIGID




SUBROUTINE GENRIGID_POTENTIAL(COORDS,ENERGY,VNEW,GTEST,STEST,RMS,PTEST,BOXTEST)

  USE COMMONS, only: NATOMS
  IMPLICIT NONE

  DOUBLE PRECISION, DIMENSION(3*NATOMS) :: COORDS, XRIGIDCOORDS
  DOUBLE PRECISION ENERGY
  DOUBLE PRECISION, DIMENSION(3*NATOMS) :: VNEW 
  LOGICAL GTEST, STEST
  DOUBLE PRECISION RMS
  LOGICAL PTEST, BOXTEST, TEMPATOMRIGIDCOORDT

  XRIGIDCOORDS(:) = 0.0D0
  CALL TRANSFORMCTORIGID(COORDS, XRIGIDCOORDS(1:DEGFREEDOMS))
  TEMPATOMRIGIDCOORDT = ATOMRIGIDCOORDT
  ATOMRIGIDCOORDT = .FALSE.
  CALL POTENTIAL(XRIGIDCOORDS,ENERGY,VNEW,GTEST,STEST,RMS,PTEST,BOXTEST)
  ATOMRIGIDCOORDT = TEMPATOMRIGIDCOORDT

END SUBROUTINE GENRIGID_POTENTIAL


SUBROUTINE GENRIGID_PERMDIST()
  USE KEY,ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, SETS
  USE COMMONS, ONLY: NATOMS
  IMPLICIT NONE

  INTEGER J1, J2, J3, J4, J5, J6,  NDUMMY, NDUMMY2
  INTEGER, ALLOCATABLE :: LRBNDUMMY(:)
  LOGICAL :: INTRAGROUPPERM

  ALLOCATE (LRBNPERMGROUP(NRIGIDBODY))
  ALLOCATE (LRBNPERMSIZE(NRIGIDBODY,3*NATOMS))
  ALLOCATE (LRBPERMGROUP(NRIGIDBODY,3*NATOMS))
  ALLOCATE (LRBNSETS(NRIGIDBODY,3*NATOMS))
  ALLOCATE (LRBSETS(NRIGIDBODY,NATOMS,3))
  ALLOCATE (LRBNDUMMY(NRIGIDBODY))
  
  LRBNPERMGROUP(:) = 0
  LRBPERMGROUP(:,:)=0
  LRBNSETS(:,:)=0
  LRBSETS(:,:,:)=0
  NDUMMY = 1
  LRBNDUMMY(:) = 1

  DO J1 = 1, NPERMGROUP
!    write(*,*) "Looking at permgroup ", J1, "which is"
!    write(*,*) permgroup(NDUMMY:NDUMMY+NPERMSIZE(J1)-1)
     DO J2 = 1, NRIGIDBODY
        DO J3 = 1, NSITEPERBODY(J2)
!           write(*,*) "looking at site", J3, "in body", J2, "which is atom", rigidgroups(J3,J2)
           IF ( PERMGROUP(NDUMMY) .EQ. RIGIDGROUPS(J3,J2) ) THEN
              NDUMMY2 = NDUMMY
!jdf43> check that permutations are restricted to a single RB
              INTRAGROUPPERM=.TRUE.
!                WRITE(*,*) "Testing the remaining atoms in this group against the rest of this permutation"
              DO J4 = 1, NPERMSIZE(J1)
!                WRITE(*,*) "Testing the group", RIGIDGROUPS(:,J2), "against atom", PERMGROUP(NDUMMY2+J4-1)
!                 IF (.NOT.ANY(RIGIDGROUPS(:,J2)==PERMGROUP(NDUMMY2+J4-1))) INTRAGROUPPERM=.FALSE.
                  ! sn402: changed the following line because it seemed to be wrong.
!                 INTRAGROUPPERM = INTRAGROUPPERM .AND. (.NOT.ANY(RIGIDGROUPS(:,J2)==PERMGROUP(NDUMMY2+J4-1)))
                 INTRAGROUPPERM = INTRAGROUPPERM .AND. (ANY(RIGIDGROUPS(:,J2)==PERMGROUP(NDUMMY2+J4-1)))
              ENDDO
              IF (INTRAGROUPPERM) THEN
!jdf43>
!		        WRITE(*,*) "Setting LRBPERMGROUP number", J2 !sn402
                 LRBNPERMGROUP(J2) = LRBNPERMGROUP(J2) + 1
!                 WRITE(*,*) "Setting group ", J2, "to have size", NPERMSIZE(J1) !sn402
                 LRBNPERMSIZE(J2,LRBNPERMGROUP(J2)) = NPERMSIZE(J1)
                 LRBNSETS(J2,LRBNPERMGROUP(J2)) = NSETS(J1)
                 DO J4 = 1, NPERMSIZE(J1)
                    DO J5 = 1, NSITEPERBODY(J2)
                       IF (PERMGROUP(NDUMMY2+J4-1) .EQ. RIGIDGROUPS(J5,J2)) THEN
!                       WRITE(*,*) "Adding atom", J5, "to position (",J2,",",LRBNDUMMY(J2)+J4-1,")"
                       ! sn402: I think this line was wrong as well?
                          LRBPERMGROUP(J2,LRBNDUMMY(J2)+J4-1) = J5
!                           LRBPERMGROUP(J2,LRBNDUMMY(J2)+J4-1) = RIGIDGROUPS(J5,J2)
                       ENDIF
                    ENDDO
                       
                    DO J6 = 1, NSETS(J1)
                       DO J5 = 1, NSITEPERBODY(J2)
                          IF (SETS(PERMGROUP(NDUMMY2+J4-1),J6) .EQ. RIGIDGROUPS(J5,J2)) THEN
                             LRBSETS(J2,LRBPERMGROUP(J2,LRBNDUMMY(J2)+J4-1),J6) = J5
!                             PRINT *, J2,LRBPERMGROUP(J2,LRBNDUMMY(J2)+J4-1), LRBNPERMGROUP(J2)
                          ENDIF
                       ENDDO
                    ENDDO
                 ENDDO
                 LRBNDUMMY(J2) = LRBNDUMMY(J2) + NPERMSIZE(J1)
!jdf43>
              ENDIF
!jdf43>
           ENDIF
           
        ENDDO
        
     ENDDO

     NDUMMY = NDUMMY + NPERMSIZE(J1)
  ENDDO

END SUBROUTINE GENRIGID_PERMDIST

! ---------------------------------------------------------------------------------
! sn402: modified version of minpermdist which aligns only a subset of rigid bodies in the system
SUBROUTINE ALIGN_RBS(COORDSA, COORDSB, DEBUG, BULKT, TWOD, DISTANCE, DIST2, RMATBEST)
! NOTE: I haven't implemented BULKT or TWOD yet, but they should be easy to do.

USE KEY, ONLY: ALIGNRBST, N_TO_ALIGN, TO_ALIGN, BULK_BOXVEC
USE COMMONS, ONLY: NATOMS

IMPLICIT NONE

DOUBLE PRECISION, INTENT(INOUT) :: COORDSA(:)
DOUBLE PRECISION, INTENT(IN) :: COORDSB(:)
LOGICAL, INTENT(IN) :: DEBUG, BULKT, TWOD
DOUBLE PRECISION, INTENT(OUT) :: DISTANCE, DIST2
DOUBLE PRECISION, INTENT(OUT) :: RMATBEST(3,3)

INTEGER :: DUMMY_NATOMS
DOUBLE PRECISION, ALLOCATABLE :: DUMMY_COORDSA(:), DUMMY_COORDSB(:)
DOUBLE PRECISION :: TRANSLATION(3)
INTEGER :: J1, J2, OFFSET

IF(DEBUG) WRITE(*,*) "genrigid> aligning subset of rigid bodies"

!DISTANCE=0.0D0
!DO J1=1,3*NATOMS
!     DISTANCE=DISTANCE+(COORDSA(J1)-COORDSB(J1))**2
!ENDDO
!IF(DEBUG) WRITE(*,*) "genrigid> Total distance before alignment=", SQRT(DISTANCE)
!write(*,*) "COORDSA before alignment", COORDSA
!write(*,*) "COORDSB before alignment", COORDSB

DUMMY_NATOMS = 0
DO J1 = 1, N_TO_ALIGN
    DUMMY_NATOMS = DUMMY_NATOMS + NSITEPERBODY(TO_ALIGN(J1))
ENDDO

ALLOCATE(DUMMY_COORDSA(3*DUMMY_NATOMS),DUMMY_COORDSB(3*DUMMY_NATOMS))

OFFSET = 0
DO J1 = 1, N_TO_ALIGN
    DO J2 = 1, NSITEPERBODY(TO_ALIGN(J1))
        DUMMY_COORDSA(OFFSET+3*J2-2:3*J2) = COORDSA(3*RIGIDGROUPS(J2,TO_ALIGN(J1))-2:3*RIGIDGROUPS(J2,TO_ALIGN(J1)))
        DUMMY_COORDSB(OFFSET+3*J2-2:3*J2) = COORDSB(3*RIGIDGROUPS(J2,TO_ALIGN(J1))-2:3*RIGIDGROUPS(J2,TO_ALIGN(J1)))
    ENDDO
    OFFSET = OFFSET+3*NSITEPERBODY(TO_ALIGN(J1))
ENDDO

ALIGNRBST = .FALSE. ! Switch this off temporarily so we don't end up back at this routine (infinite loop)

! Use MINPERMDIST to obtain RMAT, the rotation matrix that puts DUMMYCOORDSA into the best alignment
! with DUMMYCOORDSB
CALL MINPERMDIST(DUMMY_COORDSB(:),DUMMY_COORDSA(:),DUMMY_NATOMS,.FALSE., &
          BULK_BOXVEC(1),BULK_BOXVEC(2),BULK_BOXVEC(3),.FALSE.,.FALSE.,DISTANCE,DIST2,.FALSE.,RMATBEST)

IF(DEBUG) WRITE(*,*) "genrigid> Distance for aligned subset=", SQRT(DISTANCE)

ALIGNRBST = .TRUE.

! MINPERMDIST should also have aligned the centres of mass, so we can extract the required translation from DUMMY_COORDSA
TRANSLATION(:) = DUMMY_COORDSA(1:3)-COORDSA(3*RIGIDGROUPS(1,TO_ALIGN(1))-2:3*RIGIDGROUPS(1,TO_ALIGN(1)))

DIST2 = 0.0D0
! Apply the translation and rotation returned from MINPERMDIST, and calculate the optimised distance.
DO J1 = 1,NATOMS
    COORDSA(3*J1-2:3*J1) = COORDSA(3*J1-2:3*J1)+TRANSLATION(:)
    COORDSA(3*J1-2:3*J1) = MATMUL(RMATBEST,COORDSA(3*J1-2:3*J1))
    DIST2 = DIST2 + (COORDSA(3*J1-2)-COORDSB(3*J1-2))**2 + &
                    (COORDSA(3*J1-1)-COORDSB(3*J1-1))**2 + &
                    (COORDSA(3*J1)-COORDSB(3*J1))**2
ENDDO

DISTANCE=SQRT(DIST2)

IF(DEBUG) WRITE(*,*) "genrigid> Total distance after alignment=", DISTANCE

RETURN

END SUBROUTINE ALIGN_RBS
! ---------------------------------------------------------------------------------
! Pass in a list of coords arrays for configurations for the entire system, and return an equivalent list
! containing only the coordinates for a subset of the rigid bodies. Currently, this is coded to assume that the subset of RBs
! will always be the first MAX_RB bodies stored in RIGIDGROUPS.
! NOTE: although this subroutine does seem to work, I wasn't able to use it for its intended purpose so it is currently redundant.
! It turns out that once you've extracted the coordinates of a subsystem, there isn't much you can really do with them.
! I leave the code here in case it proves useful (maybe if you want to print out coordinates for part of your system only
! in order to visualise them?)
SUBROUTINE EXTRACT_SUBSET(BAND, SUBSET_BAND, MAX_RB, NCONFIGS, NOPT)
IMPLICIT NONE
INTEGER, INTENT(IN) :: MAX_RB, NCONFIGS, NOPT
DOUBLE PRECISION, INTENT(IN) :: BAND(:)
DOUBLE PRECISION, INTENT(OUT) :: SUBSET_BAND(:)
INTEGER :: J1, J2, J3, NEW_NATOMS, OFFSET

! May eventually change this so we can specify a list of RBs rather than just using the ones which come first in the list
DO J1 = 1, MAX_RB
    ! Add up how many atoms are in the subset
    NEW_NATOMS = NEW_NATOMS + NSITEPERBODY(J1)
ENDDO

IF (ATOMRIGIDCOORDT) THEN
    DO J1 = 1, NCONFIGS
        OFFSET = 0
        DO J2 = 1, MAX_RB
            DO J3 = 1, NSITEPERBODY(J2)
                SUBSET_BAND(3*NEW_NATOMS*(J1-1)+OFFSET+3*J3-2:3*NEW_NATOMS*J1+OFFSET+3*J3) = &
                    BAND(NOPT*(J1-1)+3*RIGIDGROUPS(J3,J2)-2:NOPT*(J1-1)+3*RIGIDGROUPS(J3,J2))
            ENDDO
            OFFSET = OFFSET + NSITEPERBODY(J2)  ! Keep a record of how many atoms we've already added in for this configuration
        ENDDO
    ENDDO

ELSE
    SUBSET_BAND(:) = 0.0D0
    DO J1 = 1, NCONFIGS
        DO J2 = 1, MAX_RB
            SUBSET_BAND(3*NEW_NATOMS*(J1-1)+3*J2-2:3*NEW_NATOMS*(J1-1)+3*J2) = & ! centre-of-mass coordinates
                    BAND(NOPT*(J1-1)+3*J2-2:NOPT*(J1-1)+3*J2)
            SUBSET_BAND(3*NEW_NATOMS*(J1-1)+3*MAX_RB+3*J2-2:3*NEW_NATOMS*(J1-1)+3*MAX_RB+3*J2) = & 
                    BAND(NOPT*(J1-1)+3*NRIGIDBODY+3*J2-2:NOPT*(J1-1)+3*NRIGIDBODY+3*J2)  ! AA coordinates
        ENDDO
    ENDDO
ENDIF

RETURN

END SUBROUTINE EXTRACT_SUBSET
! ---------------------------------------------------------------------------------
! sn402: angle-axis distance measure for rigid body systems
SUBROUTINE RB_DISTANCE(RB_DIST, RBCOORDS1, RBCOORDS2, GRAD1, GRAD2, GRADT)
USE COMMONS, ONLY: NATOMS, DEBUG
USE KEY, ONLY: PERMDIST
IMPLICIT NONE

! RBCOORDS1, RBCOORDS2 are angle-axis coordinates for two poses of the entire rigid-body
! system.
! GRAD1 and GRAD2 are the gradient components of the square distance with respect to
! the first and second sets of coordinates, respectively.
! SQ_DIST is the square distance measure between the two poses. RB_DIST = SQRT(SQ_DIST)
! GRADT = FALSE if we don't want to compute the distance gradient this time round.
DOUBLE PRECISION, INTENT(IN) :: RBCOORDS1(DEGFREEDOMS)
DOUBLE PRECISION, INTENT(IN) :: RBCOORDS2(DEGFREEDOMS)
DOUBLE PRECISION, INTENT(OUT) :: RB_DIST
DOUBLE PRECISION, INTENT(OUT) :: GRAD1(DEGFREEDOMS)
DOUBLE PRECISION, INTENT(OUT) :: GRAD2(DEGFREEDOMS)
LOGICAL, INTENT(IN) :: GRADT
! Temporary variables. The TEMPGRADs hold gradients of the distance measure with respect
! to the coordinates of a particular fragment. P1 and P2 are used to pass angle-axis vectors
! into the site_dist functions. drij is the shortest vector connecting two rigid-body centres
! of mass.
DOUBLE PRECISION :: TEMPGRAD1(3)
DOUBLE PRECISION :: TEMPGRAD2(3)
DOUBLE PRECISION :: P1(3)
DOUBLE PRECISION :: P2(3)
DOUBLE PRECISION :: drij(3)
DOUBLE PRECISION :: SQ_DIST
! temporary arrays to stop the ifort compiler complaining when using DEBUG mode.
DOUBLE PRECISION :: thisS(3,3), thisCOG(3)

INTEGER :: J1, J2, J5, J6, J7, J9

!IF (DEBUG) THEN
!    WRITE(*, *) "Rigid Body coords input to RB_DISTANCE"
!    DO J1 = 1, DEGFREEDOMS
!      WRITE(*, *) RBCOORDS1(J1), RBCOORDS2(J1)
!    ENDDO
!ENDIF

SQ_DIST = 0.D0
! Compute the square distance measure between each pair of rigid body fragments and add it to
! the total


! Start of big loop over all rigid fragments: compute the square distance between the two
! poses of each rigid body.
DO J1 = 1, NRIGIDBODY

    ! Need to define the smallest rij.
    drij = smallest_rij(RBCOORDS1(3*J1-2:3*J1), RBCOORDS2(3*J1-2:3*J1))

    ! Extract the angle-axis vectors of each pose
    P1 = RBCOORDS1(3*NRIGIDBODY+3*J1-2 : 3*NRIGIDBODY+3*J1)
    P2 = RBCOORDS2(3*NRIGIDBODY+3*J1-2 : 3*NRIGIDBODY+3*J1)

    ! S and COG are, respectively, lists of the non-mass-weighted tensor of gyration
    ! and the centre of geometry for all rigid bodies in the system.
    ! These are computed in SETUP_TENSORS.
    ! Note, we are using unit weights for every atom regardless of the actual atom type. This
    ! is desired for the NEB because all atoms should affect the pathway roughly equally.
    ! First, extract the S matrix and COG vector for this particular rigid body
    thisS = S(3*J1-2:3*J1,:)
    thisCOG = COG(J1,:)
    ! Compute the square distance:
!    SQ_DIST = SQ_DIST + sitedist(drij, P1, P2, thisS, RBMASS(J1), thisCOG)
    SQ_DIST = SQ_DIST + sitedist(drij, P1, P2, thisS, DBLE(NSITEPERBODY(J1)), thisCOG)

    IF(GRADT) THEN
        ! Compute the gradient of the square distance with respect to COORDS1
        ! for each pair of sites
!        CALL sitedist_grad(drij, P1, P2, thisS, RBMASS(J1), thisCOG, TEMPGRAD1, TEMPGRAD2)
        CALL sitedist_grad(drij, P1, P2, thisS, DBLE(NSITEPERBODY(J1)), thisCOG, TEMPGRAD1, TEMPGRAD2)
        GRAD1(3*J1-2:3*J1) = TEMPGRAD1
        GRAD1(3*NRIGIDBODY+3*J1-2 : 3*NRIGIDBODY+3*J1) = TEMPGRAD2
        ! Compute the gradient with respect to COORDS2 for each pair of sites
        drij = -1.0*drij

!        CALL sitedist_grad(drij, P2, P1, thisS, RBMASS(J1), thisCOG, TEMPGRAD1, TEMPGRAD2)
        CALL sitedist_grad(drij, P2, P1, thisS, DBLE(NSITEPERBODY(J1)), thisCOG, TEMPGRAD1, TEMPGRAD2)
        GRAD2(3*J1-2:3*J1) = TEMPGRAD1
        GRAD2(3*NRIGIDBODY+3*J1-2 : 3*NRIGIDBODY+3*J1) = TEMPGRAD2
    ENDIF

ENDDO ! End of big loop over rigid bodies

! Compute the contributions to the distance and gradient of any loose atoms.
DO J9 = (6*NRIGIDBODY+1), DEGFREEDOMS
    ! Straightforward cartesian distance and derivatives
    SQ_DIST = SQ_DIST + (RBCOORDS2(J9)-RBCOORDS1(J9))**2
    GRAD1(J9) = RBCOORDS2(J9) - RBCOORDS1(J9)
    GRAD2(J9) = -(RBCOORDS2(J9) - RBCOORDS1(J9))
ENDDO

!IF (DEBUG) THEN
!    WRITE(*, *) "GRAD1"
!    DO J9 = 1, DEGFREEDOMS
!        WRITE(*, *) GRAD1(J9)
!    ENDDO
!    WRITE(*, *) "GRAD2"
!    DO J9 = 1, DEGFREEDOMS
!        WRITE(*, *) GRAD2(J9)
!    ENDDO
!    WRITE(*, *) " RB Distance", SQ_DIST, SQRT(SQ_DIST)
!ENDIF
RB_DIST = SQRT(SQ_DIST)
RETURN

END SUBROUTINE RB_DISTANCE

! -----------------------------------------------------------------------------

function smallest_rij(r1, r2) result(rij)
    ! Calculate the shortest vector between two rigid-body centres of mass, with or without
    ! periodic boundary conditions
    USE KEY, only: BULKT, BULK_BOXVEC

    implicit none
    double precision, intent(in) :: r1(3), r2(3)
    double precision rij(3), ibulk_boxvec(3)
    integer j1

    do j1 = 1,3

        rij(j1) = r2(j1) - r1(j1)
        if(BULKT) then
            ibulk_boxvec(j1) = 1.d0/BULK_BOXVEC(j1)
            ! Subtract off any whole-number multiples of the box length in this dimension
            ! so that the vector connects r1 with the closest image of r2.
            rij(j1) = rij(j1) - nint(rij(j1)*ibulk_boxvec(j1))*BULK_BOXVEC(j1)
        endif
    enddo
end function smallest_rij

! ---------------------------------------------------------------------------------------------

function sitedist(drij, p1, p2, S, W, cog) result(dist)
    use commons, only: debug
    implicit none
    ! drij is a vector between the centres of mass of the two rigid body poses
    ! p1 and p2 are their angle-axis vectors
    ! DRk is the kth derivative of the rotation matrix, as computed by RMDRVT
    ! S is the non-mass-weighted tensor of gyration
    ! W is the total mass of the rigid body
    ! cog is the coordinates of the centre of mass
    ! dist is the square distance measure between the two poses
    double precision, intent(in) :: drij(3), p1(3), p2(3)
    double precision :: DR1(3,3), DR2(3,3), DR3(3,3)
    double precision, intent(in) :: S(3,3), W, cog(3)
    double precision dist
    ! R1, R2 and dR are, respectively, the rotation matrices of the two poses
    ! and the difference between these two matrices
    ! d_M, d_P, d_mix are the translational, rotational and mixed components of
    ! the distance measure (See http://dx.doi.org/10.1021/ct400403y)
    double precision R1(3,3), R2(3,3), dR(3,3)
    double precision d_M, d_P, d_mix

    ! Call RMDRVT to obtain the rotation matrices via Rodrigues' formula
    ! Note, we don't need to compute the RM derivatives this time.
    call RMDRVT(p1, R1, DR1, DR2, DR3, .FALSE.)
    call RMDRVT(p2, R2, DR1, DR2, DR3, .FALSE.)

    dR = R2 - R1

    ! Compute the various terms of the distance measure as described in the paper.
    ! Translational component of the distance measure
    d_M = W*sum((drij)**2)
    DR1 = matmul(dR, matmul(S, transpose(dR)))
    ! Rotational component
    d_P = DR1(1,1) + DR1(2,2) + DR1(3,3)
    ! Mixing component - this should be zero for a consistent choice of GR_WEIGHTS
    ! However, if the atoms are not equally weighted, then the reference coordinates for
    ! the rigid body are centred relative to the centre of mass, not the centre of geometry.
    ! The distance measure (which is used for the NEB etc.) still uses the non-mass-weighted
    ! description and hence the centre of geometry is not necessarily at (0,0,0) in reduced
    ! coordinates. In this case, we do not expect d_mix to be 0.
    d_mix = 2. * W * dot_product(drij, matmul(dR, cog))

    if (debug .and. .not. masses) then
        if (d_mix .gt. 1.0E-12) then
            print *, "ERROR STOP NOW genrigid.f90> Mixing term in square distance is > 0."
            print *, "Have you specified your masses correctly?"
            print *, "p1", p1
            print *, "p2", p2
            print *, "drij", drij
            print *, "cog", cog
            print *, "d_mix", d_mix
            stop
        endif
    endif

    dist = d_M + d_P + d_mix

end function sitedist

! ---------------------------------------------------------------------------------------------

subroutine sitedist_grad(drij, p1, p2, S, W, cog, g_M, g_P)
    use commons, only: debug
    implicit none
    ! See function sitedist, above, for meanings of variables
    double precision, intent(in) :: drij(3), p1(3), p2(3)
    double precision :: R11(3,3), R12(3,3), R13(3,3)
    double precision :: g_X(3)
    double precision, intent(in) :: S(3,3), W, cog(3)
    ! g_M is the translational part of the gradient, g_P is the rotational part
    ! and g_X is the mixing part (which should be 0 for consistent choice of weights)
    double precision, intent(out) :: g_M(3), g_P(3)
    double precision :: temp(3,3)

    double precision R1(3,3), R2(3,3), dR(3,3)

    ! Get the rotation matrices for each pose of the fragment, and the derivatives of the p1 pose.
    call RMDRVT(p2, R2, R11, R12, R13, .FALSE.)
    call RMDRVT(p1, R1, R11, R12, R13, .TRUE.)

    dR = R2 - R1

    ! Translational component of the gradient
    g_M = -2.*W*(drij)

    ! Rotational component of the gradient. The ugly use of "temp" here stops the ifort compiler in
    ! debug mode from printing out thousands of warning messages about creating an array temporary
    ! every time trace3 is called.
    temp = matmul(R11, matmul(S, transpose(dR)))
    g_P(1) = -2.*trace3(temp)
    temp = matmul(R12, matmul(S, transpose(dR)))
    g_P(2) = -2.*trace3(temp)
    temp = matmul(R13, matmul(S, transpose(dR)))
    g_P(3) = -2.*trace3(temp)

    ! Mixing component of the gradient (should be effectively zero if atoms have equal mass -
    ! see comment in the previous subroutine)
    g_X = 2.*W * matmul(dR, cog)

    if (debug .and. .not. masses) then
        if ((g_X(1)+g_X(2)+g_X(3)) .gt. 1.e-12) then
            print *, "Warning, genrigid.f90> Mixing term in distance gradient is > 0."
            print *, "Have you specified your masses correctly?"
            print *, "Mixing term = ",g_X
            stop
        endif
    endif

    g_M = g_M - g_X
    g_P(1) = g_P(1) - 2.*W * dot_product(drij, matmul(R11, cog))
    g_P(2) = g_P(2) - 2.*W * dot_product(drij, matmul(R12, cog))
    g_P(3) = g_P(3) - 2.*W * dot_product(drij, matmul(R13, cog))

end subroutine

!------------------------------------------------------------------------------------------

function trace3(M) result(tr)
    implicit none
    double precision, intent(in) :: M(3,3)
    double precision tr
    tr = M(1,1) + M(2,2) + M(3,3)
end function

!-------------------------------------------------------------------------------------------

END MODULE GENRIGID

