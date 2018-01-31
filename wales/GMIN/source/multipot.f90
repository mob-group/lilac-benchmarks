! Module that allows different potential types to be specified for different atoms in a single system

! Checklist for adding a new potential:
!   a) Write a subroutine to calculate the energy (preferably saving it in isotropic_potentials.f90 or pairwise_potentials.f90)
!   b) Make sure the signature of your subroutine matches the interface in COMPUTE_PAIRWISE_POTENTIAL or COMPUTE_ISOTROPIC_POTENTIAL
!      (see below). Is neither of these is appropriate for your potential, you will need to write a new COMPUTE subroutine.
!   c) Add your potential to the SELECT CASE statement in MULTIPOT_CALL
!   d) If your potential requires atom indices to be input in any format other than a simple list, add your requirements to the
!      SELECT CASE statement in MULTIPOT_INITIALISE
!   e) Check that you've made equivalent changes in OPTIM so that the two systems stay in sync.
!   f) Enjoy!
MODULE MULTIPOT

    USE PAIRWISE_POTENTIALS
    USE ISOTROPIC_POTENTIALS
    USE COMMONS, ONLY: DEBUG, MYUNIT

    IMPLICIT NONE

    ! We pre-calculate as much as possible, hence the need for lots of different storage arrays, documented here.
    INTEGER :: NPOTTYPES = 0                        ! The number of different potential types being used
    CHARACTER(LEN=10), ALLOCATABLE :: POTTYPES(:)   ! An array of identifiers for the different types of potential being
                                                    ! used for this particular system
    DOUBLE PRECISION, ALLOCATABLE :: POTSCALES(:)   ! A list of the energy scale factors required to match the energy units of the
                                                    ! different potentials
    DOUBLE PRECISION, ALLOCATABLE :: POTPARAMS(:,:) ! Holds the parameters required for each potential type
    INTEGER, ALLOCATABLE :: NATOM_BY_POT(:)         ! Holds the number of atoms using each potential


                                                    ! The next two variables are for use by pairwise potentials only (the entries in
                                                    ! these arrays for non-pairwise potentials must be filled by dummies)
    INTEGER, ALLOCATABLE :: N_ATOM_PARTNERS(:,:)    ! Holds the number of interaction partners for each atom, sorted by potential type
    INTEGER, ALLOCATABLE :: POTLISTS(:,:,:)         ! Holds the indices of all atoms belonging to each potential type, and
                                                    ! the indices of the partners which interact with each atom
                                                    ! POTLISTS(m,n,1) contains the index of the nth atom which has an attraction of
                                                    ! type POTTYPES(m), and POTLISTS(m,n,2:N_ATOM_PARTNERS(m,n)) are the indices of
                                                    ! its partners.

    INTEGER, ALLOCATABLE :: DEGS_BY_POT(:,:)        ! For use by isotropic potentials, where all atoms interact with all other atoms.
                                                    ! Contains a list of the degrees of freedom belonging to all the atoms using this
                                                    ! potential.


CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Parse the input file which is required to specify the different types of interaction we have in the system
SUBROUTINE MULTIPOT_INITIALISE

    USE COMMONS, ONLY: NATOMS
    USE GENRIGID, ONLY: RIGIDINIT, RB_BY_ATOM

    IMPLICIT NONE

    INTEGER NPARAMS, ATOM1, ATOM2
    INTEGER :: ATOMLIST(NATOMS)
    CHARACTER(LEN=1000) :: DUMMYCHAR
    LOGICAL :: END
    INTEGER :: J1, J2, J3, J4, iostatus, COUNTER

    ! dj337: Variables needed for binary potentials
    INTEGER :: NATOMS1, NATOMS2
    INTEGER, ALLOCATABLE :: ATOM1LIST(:), ATOM2LIST(:)

    ! Variables needed to read in exclusion lists
    INTEGER :: N_EXCLUDE_LINES, MAX_LINE_LENGTH
    INTEGER, ALLOCATABLE :: LINE_LEN(:)
    INTEGER, ALLOCATABLE :: EXCLUSIONS(:,:)  ! Holds the indices of excluded interactions. EXCLUSIONS(l,m) contains
                                             ! the m'th atom in the l'th set of excluded interactions for the current potential.

    ! Input file: multipotconfig
    ! Format as follows.
    !
    ! Introduce each type of interaction with a line: 'POT'
    ! To comment out a potential type, simply insert a '#' at the start of the 'POT' line.
    ! The next line has the form: POTTYPE n scale nparams
    !   POTTYPE is a string specifying the type of potential
    !   n is the number of atoms which will use this type of potential
    !   scale is the factor by which the energy of this interaction must be scaled to match the others
    !   nparams is the number of arguments to the pairwise potential, which will be potential specific. nparams<=10, currently.
    ! The next line is the list of input parameters (up to 10 of them), separated by spaces. 
    !   Note that nparams needs to be at least 1. If your potential has no extra parameters, set nparams to 1
    !   and have this line contain a single "0", which you then don't need to use.
    ! The next line(s) contain lists of the atoms which will interact according to the current potential.
    !   See below for example formats which may be used.
    ! Some formats (notably EWCA, ELJ) may require additional input lines in multipotconfig.

!   Read multipotconfig once to count the number of types of potential
    NPOTTYPES=0
    OPEN(UNIT=22,FILE='multipotconfig',status='old')
    DO
       READ(22,*,IOSTAT=iostatus) DUMMYCHAR
       IF (iostatus<0) THEN
          CLOSE(22)
          EXIT
       ELSE IF (DUMMYCHAR(1:3).EQ.'POT') THEN
          NPOTTYPES = NPOTTYPES + 1
       ENDIF
    END DO
    CLOSE(22)

    ALLOCATE(POTTYPES(NPOTTYPES))
    ALLOCATE(NATOM_BY_POT(NPOTTYPES))
    ALLOCATE(DEGS_BY_POT(NPOTTYPES,3*NATOMS))
    ALLOCATE(POTSCALES(NPOTTYPES))
    ALLOCATE(POTPARAMS(NPOTTYPES,10))  ! Currently each potential is allowed to have up to 10 parameters

!   Determine the names of the potentials and the number of atoms using each potential
    COUNTER = 1 ! Counts the number of types of potential that have been read thus far
    OPEN(UNIT=22,FILE='multipotconfig',status='old')
    DO
        READ(22,*,IOSTAT=iostatus) DUMMYCHAR
        IF (iostatus<0) THEN
            CLOSE(22)
            EXIT
        ELSE IF (DUMMYCHAR(1:3).EQ.'POT') THEN ! Read in information and parameters associated with this potential
            ! NPARAMS is the number of parameters to be read from the next line
            READ(22,*) POTTYPES(COUNTER), NATOM_BY_POT(COUNTER), POTSCALES(COUNTER), NPARAMS
            READ(22,*) POTPARAMS(COUNTER,:NPARAMS)
            IF(NATOM_BY_POT(COUNTER).GT.NATOMS) WRITE(MYUNIT,*) "WARNING: NATOM_BY_POT for potential ", POTTYPES(COUNTER), &
                "is larger than NATOMS"

            COUNTER = COUNTER + 1
        ENDIF
    END DO
    CLOSE(22)

    ALLOCATE(N_ATOM_PARTNERS(NPOTTYPES,MAXVAL(NATOM_BY_POT)))
    ALLOCATE(POTLISTS(NPOTTYPES,MAXVAL(NATOM_BY_POT),NATOMS))

    ! Now the important bit: read in the actual atom lists
    OPEN(UNIT=22,FILE='multipotconfig',status='unknown')
    DO J1 = 1, NPOTTYPES
        ! Read the three header lines
        READ(22,*) DUMMYCHAR

        IF (DUMMYCHAR(1:1).EQ.'#') THEN  ! This potential is commented out: skip it.
            IF(DEBUG) WRITE(MYUNIT,*) "Skipping commented block in multipotconfig"
            DO
                READ(22,*) DUMMYCHAR
                IF (DUMMYCHAR(1:3).EQ.'POT') EXIT  ! We have found the start of the next uncommented header
            ENDDO
        ENDIF

        READ(22,*) DUMMYCHAR
        READ(22,*) DUMMYCHAR


        ! Now read in the atom numbers from the following line(s).
        ! The required input format will vary between potentials, although most will use the format described under CASE DEFAULT
        SELECT CASE(POTTYPES(J1))

        ! For potentials which only operate between specific pairs of atoms, such as harmonic springs, each line should contain an interacting pair.
        ! This method of reading input should extend simply to 3- and 4-body potentials (e.g. dihedral terms).
        CASE('HSPR', 'PWCA', 'PLJ')
            DO J2 = 1, NATOM_BY_POT(J1)/2  ! J2 ranges over the number of pairs (Note, NATOM_BY_POT contains the actual no of atoms)
                ! We only want to compute the potential once for each pair.

                READ(22,*) ATOM1, ATOM2

                POTLISTS(J1,2*J2-1,1) = ATOM1  ! Make an entry for each of these atoms in POTLISTS
                POTLISTS(J1,2*J2,1) = ATOM2

                N_ATOM_PARTNERS(J1,2*J2-1)=1 ! ATOM1 partners with ATOM2...
                POTLISTS(J1,2*J2-1,2) = ATOM2
                N_ATOM_PARTNERS(J1,2*J2)=0   ! ...but atom B doesn't partner with any (we've already included its interaction with A)
            ENDDO

            IF(DEBUG) THEN
                WRITE(MYUNIT,*) "Potential:", POTTYPES(J1)
                WRITE(MYUNIT,*) "N_ATOM_PARTNERS:"
                WRITE(MYUNIT,*) N_ATOM_PARTNERS(J1,:NATOM_BY_POT(J1))
                WRITE(MYUNIT,*) "POTLISTS:"
                DO J2 = 1,NATOM_BY_POT(J1)
                    WRITE(MYUNIT,*) "Atom number", J2, "in this potential"
                    WRITE(MYUNIT,*) POTLISTS(J1,J2,:N_ATOM_PARTNERS(J1,J2))
                ENDDO
            ENDIF

        ! dj337: For pairwise potentials which operate between two types of atoms. All the atoms of type 1 interact with the atoms
        ! of type 2 but do not interact with themselves. This is essentially a binary potential by creating several 
        ! pairwise interactions. The input format is:
        !
        ! type1_natoms type2_natoms
        ! t1_atom1 t1_atom2 t1_atom3 ...
        ! t2_atom1 t2_atom2 t2_atom3 ...
        !
        ! where type1_natoms is the number of atoms of the first type and type2_natoms is the number of atoms of the second type.
        ! In the second line, list all the type1 atoms using this potential on a single line, separated by spaces.
        ! In the third line, list all the type2 atoms.
        CASE('BLJ')
            READ(22,*) DUMMYCHAR
            READ(22,*) NATOMS1, NATOMS2

            ! allocate atom lists for type1 and type2 atoms
            IF (ALLOCATED(ATOM1LIST)) DEALLOCATE(ATOM1LIST)
            IF (ALLOCATED(ATOM2LIST)) DEALLOCATE(ATOM2LIST)
            ALLOCATE(ATOM1LIST(NATOMS1))
            ALLOCATE(ATOM2LIST(NATOMS2))

            ! read atom lists
            READ(22,*) ATOM1LIST(:NATOMS1)
            READ(22,*) ATOM2LIST(:NATOMS2)

            DO J2 = 1, NATOMS1
                ! create entries for all type1 atoms
                POTLISTS(J1,J2,1) = ATOM1LIST(J2)
                N_ATOM_PARTNERS(J1,J2) = NATOMS2
                ! populate partners with type2 atoms
                DO J3 = 1, NATOMS2
                    POTLISTS(J1,J2,1+J3) = ATOM2LIST(J3)
                ENDDO
            ENDDO

            ! create entries for all type2 atoms and set n_partners to 0
            DO J3 = 1, NATOMS2
                POTLISTS(J1,NATOMS1+J3,1) = ATOM2LIST(J3)
                N_ATOM_PARTNERS(J1,NATOMS1+J3) = 0
            ENDDO

        ! dj337: Coulomb potential where all atoms interact with one another
        ! The input format is the same as for ILJ and IWCA (defined below)
        CASE('ICOU')
            READ(22,*) ATOMLIST(:NATOM_BY_POT(J1))

            ! create entries for all atoms
            DO J2=1, NATOM_BY_POT(J1)
               N_ATOM_PARTNERS(J1,J2) = 0
               POTLISTS(J1,J2,1) = ATOMLIST(J2)
               ! populate entries, preventing double counting
               DO J3=J2+1, NATOM_BY_POT(J1)
                  N_ATOM_PARTNERS(J1,J2) = N_ATOM_PARTNERS(J1,J2) + 1
                  POTLISTS(J1,J2,1+N_ATOM_PARTNERS(J1,J2)) = ATOMLIST(J3)
               ENDDO
           ENDDO

        ! dj337: Coulomb potential where all type1 atoms interact with all type2 atoms
        ! The input format is the same as for BLJ (defined above)
        CASE('BCOU')
            READ(22,*) DUMMYCHAR
            READ(22,*) NATOMS1, NATOMS2

            ! allocate atom lists for type1 and type2 atoms
            IF (ALLOCATED(ATOM1LIST)) DEALLOCATE(ATOM1LIST)
            IF (ALLOCATED(ATOM2LIST)) DEALLOCATE(ATOM2LIST)
            ALLOCATE(ATOM1LIST(NATOMS1))
            ALLOCATE(ATOM2LIST(NATOMS2))

            ! read atom lists
            READ(22,*) ATOM1LIST(:NATOMS1)
            READ(22,*) ATOM2LIST(:NATOMS2)

            DO J2 = 1, NATOMS1
                ! create entries for all type1 atoms
                POTLISTS(J1,J2,1) = ATOM1LIST(J2)
                N_ATOM_PARTNERS(J1,J2) = NATOMS2
                ! populate partners with type2 atoms
                DO J3 = 1, NATOMS2
                    POTLISTS(J1,J2,1+J3) = ATOM2LIST(J3)
                ENDDO
            ENDDO

            ! create entries for all type2 atoms and set n_partners to 0
            DO J3 = 1, NATOMS2
                POTLISTS(J1,NATOMS1+J3,1) = ATOM2LIST(J3)
                N_ATOM_PARTNERS(J1,NATOMS1+J3) = 0
            ENDDO

        ! For "iso_" potentials. Every specified atom with this potential type interacts with all the others.
        ! The input format is simple: list all the atoms using this potential on a single line, separated by spaces.
        ! They don't have to be in index order, but everything will make more sense if they are!
        ! There is currently no facility to exclude the interactions within a rigid body. Use EWCA/ELJ if you need that.
        CASE('ILJ','IWCA')
            READ(22,*) ATOMLIST(:NATOM_BY_POT(J1))

            ! All we need to save for this type of potential is the list of whole-system degrees of freedom on which
            ! the potential will depend.
            DO J2=1,NATOM_BY_POT(J1)
                J3 = ATOMLIST(J2)
                DEGS_BY_POT(J1,3*J2-2) = 3*J3-2
                DEGS_BY_POT(J1,3*J2-1) = 3*J3-1
                DEGS_BY_POT(J1,3*J2)   = 3*J3
            ENDDO

            IF(DEBUG) THEN
                write(MYUNIT,*) "Potential:", POTTYPES(J1)
                write(MYUNIT,*) "DEGS_BY_POT:"
                write(MYUNIT,*) DEGS_BY_POT(J1,:NATOM_BY_POT(J1)*3)
            ENDIF

        CASE('EWCA', 'ELJ')

        ! An isotropic potential with the option to switch off particular interactions.
        !
        ! Specify the allowed interactions as normal for isotropic potentials: with a list of atom indices on a single line.
        ! If RIGIDINIT is set, then interactions within a rigid body are not calculated.
        ! Additional excluded interactions are specified by appending the following to this multipotconfig entry:
        !
        ! EXCLUDE nlines max_length
        ! line_len_1
        ! index1 index2 ...
        ! line_len_2
        ! index1 index2 ...
        ! ...
        !
        ! nlines is the number of sets of interactions you wish to exclude
        ! max_length is the largest number of atoms involved in an excluded set (interactions will not be calculated for any
        ! pair of atoms within the set)
        ! The next nlines pairs of lines specify each excluded set. line_len_1 specifies the number of atoms in excluded set 1,
        ! and the following line contains the indices of these atoms.
        !
        ! If this potential type is the last one in your MULTIPOTCONFIG file, then you must include a final line 

            ! Read the list of interacting atoms as normal
            READ(22,*) ATOMLIST(:NATOM_BY_POT(J1))

            ! Parse the exclusions lists. Read the header line once to check that it is an EXCLUDE header line as expected.
            READ(22,*) DUMMYCHAR
            IF (DUMMYCHAR(1:7) .EQ. 'EXCLUDE') THEN
                ! Re-read this line to extract the parameters that we want
                BACKSPACE(22)
                READ(22,*) DUMMYCHAR, N_EXCLUDE_LINES, MAX_LINE_LENGTH

                ALLOCATE(LINE_LEN(N_EXCLUDE_LINES), EXCLUSIONS(N_EXCLUDE_LINES,MAX_LINE_LENGTH))
           		
                ! Read in the lists of atoms for which interactions will be excluded
            	DO J2 = 1, N_EXCLUDE_LINES
                    READ(22,*) LINE_LEN(J2)  ! This line tells us how many indices to read from the next line
                    READ(22,*) EXCLUSIONS(J2,:LINE_LEN(J2))
                ENDDO

            ELSE
                WRITE(MYUNIT,*) "multipot> Warning. An exclusion potential is being used with no EXCLUDE lists specified."
                WRITE(MYUNIT,*) "multipot> It's probably more efficient to use the ISO- version of this potential instead."
                BACKSPACE(22)  ! Step back one record so that we can move on to the next potential
            ENDIF

            ! We do not calculate interactions which have been explicitly excluded, or interactions within rigid bodies.
            ! If there are no rigid bodies and no exclusions have been specified, we can use a simpler method to construct
            ! POTLISTS (see the ELSE block)
            IF(RIGIDINIT .OR. ALLOCATED(EXCLUSIONS)) THEN  

                ! We need the list of which rigid body each atom belongs to. GENRIGID_READ_FROM_FILE must already have been
                ! called. For GMIN, this IF block should never be satisfied but for OPTIM the RIGIDINIT keyword needs to be
                ! placed higher in the odata file than MULTIPOT.
                IF (.NOT. ALLOCATED(RB_BY_ATOM)) THEN
                    WRITE(*,*) "multipot> Error: RB_BY_ATOM is not allocated by the time MULTIPOT_INITIALISE is called."
                    WRITE(*,*) "multipot> The RIGIDINIT keyword must precede the MULTIPOT keyword in your odata file."
                    STOP 98
                ENDIF
				
                ! Work out which interactions are still allowed after accounting for exclusions, and put them into POTLISTS.
                DO J2=1,NATOM_BY_POT(J1)
                    POTLISTS(J1,J2,1) = ATOMLIST(J2)   ! The first element in POTLISTS(x,y,:) is the index of atom y in the list
                                                       ! for potential type x. The remaining elements are the partners of atom y.
                    N_ATOM_PARTNERS(J1,J2) = 0

                    ! Work out which other atoms in this potential will interact with POTLISTS(J1,J2,1).
                    ! Only include each pair once to avoid double counting, so the following loop runs over indices higher than J2.
                    inner: DO J3=J2+1,NATOM_BY_POT(J1)
                    
                        IF ((.NOT. RIGIDINIT) .OR. (RB_BY_ATOM(ATOMLIST(J3)) .NE. RB_BY_ATOM(ATOMLIST(J2)))) THEN
                            ! Only add partners that belong to different rigid bodies 
                            ! (always fulfilled if there are no rigid bodies!)

                            IF(ALLOCATED(EXCLUSIONS)) THEN
                                ! Check to see if this interaction was excluded.
                                DO J4=1,N_EXCLUDE_LINES
                                    IF ((ANY(EXCLUSIONS(J4,:)==ATOMLIST(J2))) .AND. (ANY(EXCLUSIONS(J4,:)==ATOMLIST(J3)))) THEN
                                        CYCLE inner  ! Don't incude this interaction.
                                    ENDIF
                                ENDDO
                            ENDIF

                            ! If we have got this far, this pair must be a genuine interacting pair.
                            ! Avoid double counting: each interacting pair is only included once. 
                            ! (The fact that J3 starts from J2+1 ensures this)
                            N_ATOM_PARTNERS(J1,J2) = N_ATOM_PARTNERS(J1,J2) + 1
                            POTLISTS(J1,J2,1+N_ATOM_PARTNERS(J1,J2)) = ATOMLIST(J3)
                        ENDIF
                    ENDDO inner
                ENDDO
            ELSE   
               ! No rigid bodies or exclusions so all atoms in ATOMLIST should partner with all others.
                DO J2 = 1,NATOM_BY_POT(J1)  ! We still only add each pair once.
                    N_ATOM_PARTNERS(J1,J2) = NATOM_BY_POT(J1)-J2  ! So each atom partners with every atom that comes later in the list
                    POTLISTS(J1,J2,:N_ATOM_PARTNERS(J1,J2)+1)=ATOMLIST(J2:NATOM_BY_POT(J1))
                ENDDO
            ENDIF        

            IF(DEBUG) THEN
                WRITE(MYUNIT,*) "Potential:", POTTYPES(J1)
                WRITE(MYUNIT,*) "Number of atoms:", NATOM_BY_POT(J1)
                WRITE(MYUNIT,*) "N_ATOM_PARTNERS:"
                WRITE(MYUNIT,*) N_ATOM_PARTNERS(J1,:NATOM_BY_POT(J1))
                WRITE(MYUNIT,*) "POTLISTS:"
                DO J2 = 1,NATOM_BY_POT(J1)
                    WRITE(MYUNIT,*) "Atom number", J2, "in this potential"
                    WRITE(MYUNIT,*) POTLISTS(J1,J2,:1+N_ATOM_PARTNERS(J1,J2))
                ENDDO
            ENDIF

        CASE DEFAULT
           WRITE(MYUNIT,*) "multipot> Error: no rule to read input for this type of potential."
           STOP
        END SELECT
    ENDDO
    CLOSE(22)

END SUBROUTINE MULTIPOT_INITIALISE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculate the potential energy and the gradient (if GTEST is set) and hessian (if STEST)
SUBROUTINE MULTIPOT_CALL (X, G, ENERGY, GTEST, STEST)
    USE COMMONS, ONLY: NATOMS
    USE MODHESS

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN)  :: X(3*NATOMS)
    DOUBLE PRECISION, INTENT(OUT) :: ENERGY, G(3*NATOMS)
    LOGICAL, INTENT(IN)           :: GTEST, STEST
    INTEGER                       :: J1

    INTEGER                       :: THISN
    DOUBLE PRECISION              :: THESEPARAMS(10)

    ENERGY  = 0.D0
    IF (GTEST) G(:) = 0.D0
    IF (STEST) THEN
        HESS(:,:) = 0.D0
    ENDIF

    ! Cycle through the different types of interaction and perform the potential call for all the corresponding atoms
    DO J1 = 1, NPOTTYPES

        ! Extract a few sub-arrays manually to avoid ifort's irritating "array temporary" warning message.
        THISN = NATOM_BY_POT(J1)
        THESEPARAMS = POTPARAMS(J1,:)

        SELECT CASE(POTTYPES(J1))
            CASE('ILJ')
                ! Only one parameter for this potential: the particle radius, sigma
                CALL COMPUTE_ISOTROPIC_POTENTIAL(X, G, ENERGY, GTEST, STEST, ISO_LJ, J1, THISN)
            CASE('IWCA')
                ! Only one parameter for this potential: the particle radius, sigma
                CALL COMPUTE_ISOTROPIC_POTENTIAL(X, G, ENERGY, GTEST, STEST, ISO_WCA, J1, THISN)
            CASE('PLJ')
                ! Only one parameter for this potential: the particle radius, sigma
                CALL COMPUTE_PAIRWISE_POTENTIAL(X, G, ENERGY, GTEST, STEST, PAIRWISE_LJ, J1)
            CASE('PWCA')
                ! Only one parameter for this potential: the particle radius, sigma
                CALL COMPUTE_PAIRWISE_POTENTIAL(X, G, ENERGY, GTEST, STEST, PAIRWISE_WCA, J1)
            CASE('BLJ')
                CALL COMPUTE_PAIRWISE_POTENTIAL(X, G, ENERGY, GTEST, STEST, PAIRWISE_LJ, J1)
            CASE('HSPR')
                ! Parameter is the equilibrium bond length, R0. Energy is returned in units of k_spr.
                CALL COMPUTE_PAIRWISE_POTENTIAL(X, G, ENERGY, GTEST, STEST, HARMONIC_SPRINGS, J1)
            ! For exclusion potentials, we must also pass in a list specifying the number of interacting partners possessed
            ! by each atom, and a list specifying which atoms make up these partners (from POTLISTS)
            CASE('ELJ')
                ! Only one parameter for this potential: the particle radius, sigma
                CALL EXCLUDE_ISO_LJ(X, POTLISTS(J1,:THISN,:), N_ATOM_PARTNERS(J1,:THISN), POTSCALES(J1), THESEPARAMS, ENERGY, G, GTEST, STEST)
            CASE('EWCA')
                ! Only one parameter for this potential: the particle radius, sigma
                CALL EXCLUDE_ISO_WCA(X, POTLISTS(J1,:THISN,:), N_ATOM_PARTNERS(J1,:THISN), POTSCALES(J1), THESEPARAMS, ENERGY, G, GTEST, STEST)
            CASE('ICOU')
                ! Two parameters for this potential: the charges for the atoms (these numbers should be the same)
                CALL ISO_COULOMB(X, POTLISTS(J1,:THISN,:), N_ATOM_PARTNERS(J1,:THISN), POTSCALES(J1), THESEPARAMS, ENERGY, G, GTEST, STEST)
            CASE('BCOU')
                ! Two parameters for this potential: the charge for type1 atom and for type2 atom (they should be different)
                CALL ISO_COULOMB(X, POTLISTS(J1,:THISN,:), N_ATOM_PARTNERS(J1,:THISN), POTSCALES(J1), THESEPARAMS, ENERGY, G, GTEST, STEST)
            CASE DEFAULT
                ! We shouldn't every get here, unless you have added a new type of potential to MULTIPOT_INITIALISE and forgotten
                ! to add it to this SELECT CASE.
                WRITE(MYUNIT,*) "multipot> Error: unspecified potential"
                STOP
        END SELECT
    ENDDO


END SUBROUTINE MULTIPOT_CALL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Evaluate the energy(+gradient(+hessian)) for a set of atoms interacting according to this particular potential
SUBROUTINE COMPUTE_ISOTROPIC_POTENTIAL(X, G, ENERGY, GTEST, STEST, POT, POTID, TMP_NATOMS)
    USE COMMONS, ONLY: NATOMS
    USE MODHESS
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN)    :: X(3*NATOMS)
    DOUBLE PRECISION, INTENT(INOUT) :: G(3*NATOMS)
    DOUBLE PRECISION, INTENT(INOUT) :: ENERGY
    LOGICAL, INTENT(IN)             :: GTEST, STEST
    INTEGER, INTENT(IN)             :: POTID, TMP_NATOMS
    DOUBLE PRECISION                :: TMP_X(3*TMP_NATOMS), TMP_G(3*TMP_NATOMS), TMP_HESS(3*TMP_NATOMS,3*TMP_NATOMS), TMP_E
    INTEGER                         :: J1, J2, NDEGS

    ! Interface to the potential type which is passed in
    INTERFACE
        SUBROUTINE POT(TMP_NATOMS, X, PARAMS, TMP_ENERGY, TMP_G, TMP_HESS, GTEST, STEST)
            INTEGER, INTENT(IN)           :: TMP_NATOMS
            DOUBLE PRECISION, INTENT(IN)  :: X(3*TMP_NATOMS)
            DOUBLE PRECISION, INTENT(IN)  :: PARAMS(10)  ! Maximum number of parameters is hardcoded here
            DOUBLE PRECISION, INTENT(OUT) :: TMP_ENERGY
            DOUBLE PRECISION, INTENT(OUT) :: TMP_G(3*TMP_NATOMS), TMP_HESS(3*TMP_NATOMS,3*TMP_NATOMS)
            LOGICAL, INTENT(IN)           :: GTEST, STEST
        END SUBROUTINE POT
    END INTERFACE


    NDEGS = 3*TMP_NATOMS

    DO J1 = 1, NDEGS  ! Loop over all the atoms with this kind of potential
        TMP_X(J1) = X(DEGS_BY_POT(POTID,J1))
    ENDDO

!    IF(DEBUG) THEN
!        WRITE(MYUNIT,*) "Calling potential", POTTYPES(J1)
!        WRITE(MYUNIT,*) "Degrees of freedom being used:"
!        WRITE(MYUNIT,*) DEGS_BY_POT(POTID,:NDEGS)
!    ENDIF

    ! The advantage of this slightly tortuous method is that we only need to call POT once. Compare with the pairwise
    ! routine, where we need to make a function call for every pair of atoms in the system.
    CALL POT(TMP_NATOMS, TMP_X, POTPARAMS(POTID,:), TMP_E, TMP_G, TMP_HESS, GTEST, STEST)

    ENERGY = ENERGY + TMP_E*POTSCALES(POTID)

    ! Unpack the gradient, if it is being calculated
    IF(GTEST) THEN
        DO J1 = 1,NDEGS
            G(DEGS_BY_POT(POTID,J1)) = G(DEGS_BY_POT(POTID,J1)) + TMP_G(J1)*POTSCALES(POTID)
        ENDDO
    ENDIF

    ! Unpack the Hessian, if it is being calculated
    IF(STEST) THEN
        DO J1 = 1,NDEGS
            DO J2 = 1,NDEGS
                HESS(DEGS_BY_POT(POTID,J1), DEGS_BY_POT(POTID,J2)) = HESS(DEGS_BY_POT(POTID,J1), DEGS_BY_POT(POTID,J2)) &
                                                                               + TMP_HESS(J1,J2)*POTSCALES(POTID)
            ENDDO
        ENDDO
    ENDIF

END SUBROUTINE COMPUTE_ISOTROPIC_POTENTIAL

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!   Evaluate the energy(+gradient(+hessian)) for each pair of atoms listed as interacting according to this particular potential
SUBROUTINE COMPUTE_PAIRWISE_POTENTIAL(X, G, ENERGY, GTEST, STEST, POT, POTID)
    USE COMMONS, ONLY: NATOMS
    USE MODHESS
    IMPLICIT NONE
    DOUBLE PRECISION, INTENT(IN)    :: X(3*NATOMS)
    DOUBLE PRECISION, INTENT(INOUT) :: G(3*NATOMS)
    DOUBLE PRECISION, INTENT(INOUT) :: ENERGY
    LOGICAL, INTENT(IN)             :: GTEST, STEST
    INTEGER, INTENT(IN)             :: POTID
    DOUBLE PRECISION                :: TMP_E, TMP_G(6), TMP_HESS(6,6)
    INTEGER                         :: J1, J2, J3, ATOM1, ATOM2

    ! Interface to the potential type which is passed in
    INTERFACE
        SUBROUTINE POT(X1, X2, PARAMS, PG, PAIR_ENERGY, P_HESS, GTEST, STEST)
            DOUBLE PRECISION, INTENT(IN)  :: X1(3), X2(3)
            DOUBLE PRECISION, INTENT(IN)  :: PARAMS(10)  ! Maximum number of parameters is hardcoded here
            DOUBLE PRECISION, INTENT(OUT) :: PAIR_ENERGY
            DOUBLE PRECISION, INTENT(OUT) :: PG(6), P_HESS(6,6)
            LOGICAL, INTENT(IN)           :: GTEST, STEST
        END SUBROUTINE POT
    END INTERFACE

    DO J1 = 1, NATOM_BY_POT(POTID)  ! Loop over all the atoms with this kind of potential
        ATOM1 = POTLISTS(POTID,J1,1)
        DO J2 = 1, N_ATOM_PARTNERS(POTID,J1)  ! Loop over every partner with which this atom must interact
            ATOM2 = POTLISTS(POTID,J1,1+J2)

            !!!!!! The actual potential call! !!!!!
            CALL POT(X(3*ATOM1-2:3*ATOM1),X(3*ATOM2-2:3*ATOM2),POTPARAMS(POTID,:),TMP_G,TMP_E,TMP_HESS,GTEST,STEST)
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ENERGY = ENERGY + TMP_E*POTSCALES(POTID)

            ! Unpack the pairwise gradient, if required
            IF(GTEST) THEN
                G(3*ATOM1-2:3*ATOM1) = G(3*ATOM1-2:3*ATOM1) + TMP_G(:3)*POTSCALES(POTID)
                G(3*ATOM2-2:3*ATOM2) = G(3*ATOM2-2:3*ATOM2) + TMP_G(4:)*POTSCALES(POTID)

                ! Unpack the pairwise Hessian, if required
                IF(STEST) THEN
                    DO J3 = 1, 3
                        ! On-diagonal block for atom 1
                        HESS(3*(ATOM1-1)+J3,3*ATOM1-2:3*ATOM1) = HESS(3*(ATOM1-1)+J3,3*ATOM1-2:3*ATOM1) + &
                                                                        TMP_HESS(J3,1:3)*POTSCALES(POTID)
                        ! On-diagonal block for atom 2
                        HESS(3*(ATOM2-1)+J3,3*ATOM2-2:3*ATOM2) = HESS(3*(ATOM2-1)+J3,3*ATOM2-2:3*ATOM2) + &
                                                                        TMP_HESS(J3+3,4:6)*POTSCALES(POTID)
                        ! Top-right off-diagonal block
                        HESS(3*(ATOM1-1)+J3,3*ATOM2-2:3*ATOM2) = HESS(3*(ATOM1-1)+J3,3*ATOM2-2:3*ATOM2) + &
                                                                        TMP_HESS(J3,4:6)*POTSCALES(POTID)
                        ! Bottom-left off-diagonal block
                        HESS(3*(ATOM2-1)+J3,3*ATOM1-2:3*ATOM1) = HESS(3*(ATOM2-1)+J3,3*ATOM1-2:3*ATOM1) + &
                                                                        TMP_HESS(J3+3,1:3)*POTSCALES(POTID)
                    ENDDO
                ENDIF

            ENDIF

        ENDDO
    ENDDO

END SUBROUTINE COMPUTE_PAIRWISE_POTENTIAL

END MODULE MULTIPOT
