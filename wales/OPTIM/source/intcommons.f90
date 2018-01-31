MODULE INTCOMMONS
  USE modamber9, only : INTMINPERMT
  ! some global variables for internals
  IMPLICIT NONE  

  ! keywords
  LOGICAL :: NATINT = .FALSE., INTNEWT = .TRUE.
  LOGICAL :: OLDINTMINPERMT=.FALSE. !this is BAD for amber
  LOGICAL :: BBCART = .FALSE., INTINTERPT = .FALSE., INTERPSIMPLE = .FALSE.
  LOGICAL :: RIGIDBONDS = .FALSE. ! freeze bonds when using internal coordinates
  LOGICAL :: BONDSFROMFILE = .FALSE. ! read bondlength from file
  LOGICAL :: INTERPCHOICE = .FALSE.
  INTEGER :: NINTIM = 100 ! number of images to define path in internals interpolation
  ! starting from this residue use cartesians
  ! if cartresstart is less than zero, it will be set to NRES+1
  INTEGER :: CARTRESSTART = -1, CARTATMSTART = -1
  DOUBLE PRECISION :: MINBACKTCUT = 1.0D-4, INTERPBACKTCUT = 1.0D-8
  LOGICAL :: PRINTCOORDS
  LOGICAL :: BACKTPRINTERR = .FALSE.
  LOGICAL :: DESMINT = .FALSE., CHICDNEB=.FALSE.
  INTEGER :: NURINGS=0, URINGS(100,0:6)
  INTEGER :: NUBONDS=0, UBONDS(100,2)
  LOGICAL :: USEPARFILE = .FALSE.
  CHARACTER*100 :: INTPARFILE = 'naturals-noH.par'
  DOUBLE PRECISION :: GMAXINT = -1  ! convergence check in internals. a negative value means its disabled,

  DOUBLE PRECISION, PARAMETER :: PI=3.141592653589793D0, TWOPI = 2*PI

  ! -----------------------
  ! Natural internals stuff
  ! -----------------------
  ! number of bonds, centers, rings, fused rings, linear dihedrals, 
  ! improper dihedrals, cartesian atoms, and individual dihedral angles
  INTEGER :: NBDS, NCNT, NRNG, NFRG, NLDH, NIMP, NCRT, NDIH
  ! angle and torsion coefficients defining the natural coordinates
  DOUBLE PRECISION :: COEFF(30,6,6) 
  !     CENTERS lists type of center, then center atom followed by nonterminal atoms
  !     followed by terminal atoms; RINGS lists atoms in ring in order; 
  !     FRINGS lists atoms involved in a 
  !     fused ring, the two bridge atoms first; LINDIH lists m, n, A, B, 
  !     X_i, Y_j for X_m -- A -- B -- Y_n; IMPDIH lists atoms in an improper dihedral 
  !     ('wag') coordinate; starting with central atom and ending with the wagging
  !     atom.
  !     ATOMxPOSINC gives, for a given center type, the position of the x coordinate 
  !     of an atom when for each natural angle-sum coordinate all the atoms involved 
  !     are listed in order with triplets for the x,y,z coordinates
  INTEGER, ALLOCATABLE :: CENTERS(:,:), RINGS(:,:), FRINGS(:,:), LINDIH(:,:), &
       & IMPDIH(:,:), CARTATMS(:)
  INTEGER, ALLOCATABLE :: CENTER2(:) !msb50 - lists center nos if centertyp is 2
  INTEGER :: NCNT2 
  INTEGER :: ATOMxPOSINC(10,5,5)

  ! -------------------
  DOUBLE PRECISION, ALLOCATABLE, TARGET :: DIHINFOSINGLE(:), DIHINFO(:,:)
  DOUBLE PRECISION, POINTER :: PREVDIH(:)
  LOGICAL :: ALIGNDIR = .FALSE.
  LOGICAL, ALLOCATABLE :: BACKBONE(:), USECART(:)
  INTEGER :: NINTC = 0
  INTEGER :: KD, NNZ
  DOUBLE PRECISION :: BACKTCUTOFF

  ! copied from charmm common blocks
  INTEGER :: TOTBOND, TOTANG, TOTDIH, TOTRES! total bonds, angles, dihedrals, residues in molecule
  INTEGER, ALLOCATABLE :: RESSTARTS(:) ! starting atoms (n-1) for each residue
  CHARACTER*8, ALLOCATABLE :: ATMTYPES(:) ! list of atom types
  CHARACTER*8, ALLOCATABLE :: RESLIST(:) ! list of residues
  INTEGER:: STARTLINDH
  INTEGER, ALLOCATABLE :: MBBONDNUM(:), MBADJACENT(:,:)
  INTEGER, ALLOCATABLE :: BONDNUM(:), ADJACENT(:,:)
  INTEGER, ALLOCATABLE :: PERMNEIGHBOURS(:,:), PERMCHAIN(:,:)
  DOUBLE PRECISION, ALLOCATABLE :: PERMNANGLE(:,:)
  INTEGER, ALLOCATABLE :: GLYDIH(:,:)
  INTEGER ::NGLYDIH
  LOGICAL :: GLYCART = .FALSE.
  LOGICAL :: INTDISTANCET = .FALSE. !distances from minpermdist in dihedral angle metric

END MODULE INTCOMMONS
