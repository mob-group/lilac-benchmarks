MODULE CONFIGS
USE PRECISION, ONLY: REAL64
#ifdef AMBER_12
USE AMBER12_INTERFACE_MOD, ONLY : POT_ENE_REC_C
#endif
IMPLICIT NONE

TYPE CONFIG
! Type for storing configurations in GMIN or OPTIM.

! Coords and information about our coordinate system. COORD_SYSTEM refers to a
! series of named parameters, described below.
   REAL(REAL64), ALLOCATABLE  :: COORDS(:)
   INTEGER                    :: COORD_SYSTEM
   INTEGER                    :: NUM_DIMS

! Energy terms
   REAL(REAL64)               :: POT_ENERGY
   REAL(REAL64)               :: FREE_ENERGY
#ifdef AMBER12
   TYPE(POT_ENE_REC_C)        :: POT_ENERGY_DECOMP
#endif
   REAL(REAL64), ALLOCATABLE  :: CUSTOM_ENERGY_DECOMP(:)

! First derivative terms (i.e. gradient/force)
   REAL(REAL64)               :: RMS_FORCE
   REAL(REAL64), ALLOCATABLE  :: GRADIENT(:)

! Second derivative terms, eigenvectors and eigenvalues
   REAL(REAL64), ALLOCATABLE  :: HESSIAN(:, :)
   REAL(REAL64), ALLOCATABLE  :: EIGENVALUES(:)
   REAL(REAL64), ALLOCATABLE  :: EIGENVECTORS(:, :)
   INTEGER                    :: NUM_ZEROS
   REAL(REAL64)               :: LOG_PROD_EVALUES

! Mass related stuff
   REAL(REAL64), ALLOCATABLE  :: MASSES(:)
   LOGICAL                    :: MASS_WEIGHTED

! Type of configuration (see parameters below).
   INTEGER                    :: CONFIG_TYPE

! Information about when the configuration was generated or found.
   INTEGER                    :: STEP
   INTEGER                    :: POT_CALLS
   REAL(REAL64)               :: TIME

! Order of the point group
   INTEGER                    :: ORDER

! Related configurations (e.g. sloppy minimum corresponding to tight minimum)
   TYPE(CONFIG), POINTER      :: ASSOC_CONFIGS(:) 
END TYPE CONFIG

! Named parameters for CONFIG_TYPE
INTEGER, PARAMETER :: MINIMUM       = 1
INTEGER, PARAMETER :: TRANS_STATE   = 2
INTEGER, PARAMETER :: AFTER_STEP    = 3

! Named parameters for COORD_SYSTEM
INTEGER, PARAMETER :: CARTESIAN     = 1
INTEGER, PARAMETER :: RIGID_BODY    = 2
INTEGER, PARAMETER :: ANGLE_AXIS    = 3
INTEGER, PARAMETER :: INTERNALS     = 4

END MODULE CONFIGS
