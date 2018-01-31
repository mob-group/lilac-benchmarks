!Subroutines to set up and use mutations as moves in a general BH run
MODULE AMBER12_MUTATIONS
  USE COMMONS
  USE PORFUNCS
  USE CHIRALITY, ONLY: DEALLOC_STATES_MUTATION
  USE AMBER12_INTERFACE_MOD
  USE QMODULE
  IMPLICIT NONE


!******************************************************************************
! Types to represent mutation information.
  TYPE RESIDUE_MUTATION
     INTEGER                           :: RESNUM          !residue number
     INTEGER                           :: NMUTATIONS      !number of mutations so far
     INTEGER                           :: NENTRIES        !number of possible residues
     CHARACTER(LEN=4)                  :: CURRENT_RES     !current residue
     CHARACTER(LEN=4) , DIMENSION(:) , ALLOCATABLE :: RESCHOICE       !residues to choose for mutations
     DOUBLE PRECISION , DIMENSION(:) , ALLOCATABLE :: PROBABILITIES   !selection probability for selection
  END TYPE RESIDUE_MUTATION

  CHARACTER(LEN=4), ALLOCATABLE , SAVE :: AMBER12_RESNAME(:)
  INTEGER , SAVE :: NRESIDUES , NRESMUT , MUNIT, NRESRB, NCYXBONDS
  INTEGER , ALLOCATABLE , SAVE :: TERMINI_RES(:), AMBER12_RESSTART(:), AMBER12_RESEND(:), AMBER12_RESNATOM(:)
  INTEGER , ALLOCATABLE , SAVE :: RESFORRB(:,:), CYX_BONDS(:,:)
  TYPE(RESIDUE_MUTATION) , DIMENSION(:) , ALLOCATABLE ,SAVE :: MUTATION_INFO , PREVIOUS_MUTATION
  DOUBLE PRECISION, ALLOCATABLE, SAVE :: CA_REFERENCE(:)
  LOGICAL :: CYX_BONDT

  CONTAINS
  !setup mutational system, initialise coordinates correctly, and according to the right size
  SUBROUTINE AMBERMUTATION_SETUP()  
     IMPLICIT NONE
     INTEGER :: MIUNIT,GETUNIT,J1,J2,NENTRIES,NTERMINI,TESTINT,CYXUNIT
     LOGICAL :: YESNO , NTERT
     CHARACTER(200) :: ENTRY_
     CHARACTER(25) , DIMENSION(:) , ALLOCATABLE :: ENTRIES
     INTEGER, PARAMETER :: NMUTOCC = 1000 

     !check there is a file contianing the mutational information
     YESNO = .FALSE.
     INQUIRE(FILE='amber_mutations',EXIST=YESNO)
     IF (.NOT.YESNO) THEN
        WRITE(MYUNIT,'(A)') ' ambermut> No mutation information given'
        STOP
     ENDIF
     !get the number of residues, and their atom positions
     CALL TOPOLOGY_READER()
     !open the mutation information
     MIUNIT = GETUNIT()
     OPEN(UNIT=MIUNIT,FILE='amber_mutations',status='unknown')
     WRITE(MYUNIT,*) 'ambermut> Reading in mutations allowed'
     READ(MIUNIT,*) NRESMUT
     WRITE(MYUNIT,*) 'ambermut> ',NRESIDUES,' residues, of which ',NRESMUT,'can be mutated'
     !all allocations and reallocations are taken care of here, as we will not chnage the number of residues
     IF (ALLOCATED(MUTATION_INFO)) DEALLOCATE(MUTATION_INFO)
     ALLOCATE(MUTATION_INFO(NRESMUT))
     IF (ALLOCATED(PREVIOUS_MUTATION)) DEALLOCATE(PREVIOUS_MUTATION)
     ALLOCATE(PREVIOUS_MUTATION(NRESMUT))
     IF (ALLOCATED(TERMINI_RES)) DEALLOCATE(TERMINI_RES)
     ALLOCATE(TERMINI_RES(NRESIDUES))
     TERMINI_RES(:) = 0
     !best sequence
     IF (ALLOCATED(BESTMUTSEQ)) DEALLOCATE(BESTMUTSEQ)
     ALLOCATE(BESTMUTSEQ(NRESIDUES))
     BESTMUTSEQ(:) = "    "
     !all sequences including rejected ones
     IF (ALLOCATED(MUTSEQ)) DEALLOCATE(MUTSEQ)
     ALLOCATE(MUTSEQ(NMUTOCC,NRESIDUES))
     MUTSEQ(:,:) = "    "
     !storing information about rejecting sequences
     IF (ALLOCATED(SEQREJECTED)) DEALLOCATE(SEQREJECTED)
     ALLOCATE(SEQREJECTED(NMUTOCC))
     SEQREJECTED(:) = .FALSE.
     NSEQSTORED = 0
     !first sequence to be stored is initial configuration
     MUTSEQ(1,:) = AMBER12_RESNAME(:)
     NSEQSTORED = NSEQSTORED + 1
     !read next line, this contains the terminal residues
     READ(MIUNIT,*) NTERMINI
     READ(MIUNIT,'(A)',END=101) ENTRY_               !read line
     ALLOCATE(ENTRIES(NTERMINI))
     ENTRIES(:)=''
     CALL READ_LINE(ENTRY_,NTERMINI,ENTRIES)
     !entries now contains the information about the termini, we need to differentiate between C and N termini
     !the first entry ought to be the N terminus and then we switch between then, this might cause problems for ACE, NME, NHE!!!
     NTERT = .TRUE.
     DO J1=1,NRESIDUES
        DO J2=1,NTERMINI
           READ(ENTRIES(J2),'(I8)') TESTINT
           IF (TESTINT.EQ.J1) THEN
              IF (NTERT) THEN
                 TERMINI_RES(J1)=1
                 NTERT = .FALSE.
              ELSE
                 TERMINI_RES(J1)=2
                 NTERT = .TRUE.
              ENDIF
           ENDIF
        ENDDO
     ENDDO
     !now we get to actual mutation information
     !line 1: RESNUM NENTRIES CURRENT_RES
     !line 2: RESNAME1 RESNAME2 RESNAME3 ...
     !line 3: PROB1 PROB2 PROB3 ...
     ! we can give the probabilities as any series of numbers, they are normalised later, as we need to discount the residue that we currently have
     !i.e. if we try to mutate we will mutate, and then check after some more group rotation steps in mc.F whether the energy is lower or not 
     DO J1=1,NRESMUT
        READ(MIUNIT,'(A)',END=101) ENTRY_               !read line
        !reallocate the length of the entries list
        DEALLOCATE(ENTRIES)
        ALLOCATE(ENTRIES(4))
        ENTRIES(:)=''
        CALL READ_LINE(ENTRY_,4,ENTRIES)
        READ(ENTRIES(1),'(I8)') MUTATION_INFO(J1)%RESNUM
        READ(ENTRIES(2),'(I8)') MUTATION_INFO(J1)%NENTRIES
        READ(ENTRIES(3),'(A)') MUTATION_INFO(J1)%CURRENT_RES
        READ(ENTRIES(4),'(I8)') MUTATION_INFO(J1)%NMUTATIONS
        NENTRIES=MUTATION_INFO(J1)%NENTRIES
        ALLOCATE(MUTATION_INFO(J1)%RESCHOICE(NENTRIES))
        ALLOCATE(MUTATION_INFO(J1)%PROBABILITIES(NENTRIES))
        ALLOCATE(PREVIOUS_MUTATION(J1)%RESCHOICE(NENTRIES))
        ALLOCATE(PREVIOUS_MUTATION(J1)%PROBABILITIES(NENTRIES))
        !reallocate the length of the entries list
        READ(MIUNIT,'(A)',END=101) ENTRY_
        DEALLOCATE(ENTRIES)
        ALLOCATE(ENTRIES(MUTATION_INFO(J1)%NENTRIES))
        ENTRIES(:)=''
        CALL READ_LINE(ENTRY_,MUTATION_INFO(J1)%NENTRIES,ENTRIES)
        DO J2=1,MUTATION_INFO(J1)%NENTRIES
           MUTATION_INFO(J1)%RESCHOICE(J2)=ENTRIES(J2)
        ENDDO
        READ(MIUNIT,'(A)',END=101) ENTRY_
        ENTRIES(:)=''
        CALL READ_LINE(ENTRY_,MUTATION_INFO(J1)%NENTRIES,ENTRIES)
        DO J2=1,MUTATION_INFO(J1)%NENTRIES
           READ(ENTRIES(J2),*) MUTATION_INFO(J1)%PROBABILITIES(J2)
        ENDDO
     ENDDO
101  CONTINUE
     CLOSE(MIUNIT)
     !check if we have Cysteine disulfide bonds!
     CYX_BONDT=.FALSE.
     INQUIRE(FILE='amber_mut_cyx',EXIST=YESNO)
     IF (YESNO) THEN
        CYX_BONDT=.TRUE.
        CYXUNIT=GETUNIT()
        OPEN(UNIT=CYXUNIT,FILE='amber_mut_cyx',STATUS='unknown')
        READ(CYXUNIT,*) NCYXBONDS
        IF (ALLOCATED(CYX_BONDS)) DEALLOCATE(CYX_BONDS)
        ALLOCATE(CYX_BONDS(NCYXBONDS,2))
        DO J1=1,NCYXBONDS
           READ(MIUNIT,'(A)',END=101) ENTRY_
           ENTRIES(:)=''
           CALL READ_LINE(ENTRY_,2,ENTRIES)
           READ(ENTRIES(1),'(I8)') CYX_BONDS(J1,1)
           READ(ENTRIES(2),'(I8)') CYX_BONDS(J1,2)
        ENDDO
        CLOSE(CYXUNIT)
     ENDIF
     !call the grouprotation set up here (not in keywords)
     CALL MUT_SETUP_GROUPROTATION(1,.FALSE.,.FALSE.,0)
!old implementation --> use file to specify residues to be rigidified 
!     IF (AMBERMUTRIGIDT) THEN
!        CALL CREATE_RIGID_FILES()
!     END IF    
     RETURN
  END SUBROUTINE AMBERMUTATION_SETUP

  !setup rigid bodies for AMBERMUTATIONS
  SUBROUTINE AMBERMUTRB_SETUP()
     IMPLICIT NONE
     INTEGER :: GETUNIT, RBUNIT, J1
     LOGICAL :: YESNO
     CHARACTER(200) :: ENTRY_
     CHARACTER(25) :: ENTRIES(2)

     !check there is a file contianing the mutational information
     YESNO = .FALSE.
     INQUIRE(FILE='residues_rb',EXIST=YESNO)
     IF (.NOT.YESNO) THEN
        WRITE(MYUNIT,'(A)') ' ambermut> No information to recreate rigid body input'
        STOP
     ENDIF
     RBUNIT = GETUNIT()
     OPEN(UNIT=RBUNIT,FILE='residues_rb',status='unknown')
     WRITE(MYUNIT,*) 'ambermut> Reading in information to create new rbodyconfig files'
     READ(RBUNIT,*) NRESRB !number of rigid bodies
     ALLOCATE(RESFORRB(NRESRB,2))
     DO J1=1,NRESRB
        READ(RBUNIT,'(A)',END=159) ENTRY_
        ENTRIES(:)=''
        CALL READ_LINE(ENTRY_,2,ENTRIES)
        READ(ENTRIES(1),'(I8)') RESFORRB(J1,1)
        READ(ENTRIES(2),'(I8)') RESFORRB(J1,2)
     ENDDO
159  CONTINUE
     CLOSE(RBUNIT)
     RETURN
  END SUBROUTINE AMBERMUTRB_SETUP
  
  !mutate protein
  SUBROUTINE AMBERMUT_STEP(COORDINATES , RESNUMBER)
     INTEGER , INTENT(OUT) :: RESNUMBER
     CHARACTER(LEN=4) :: OLDRES , OLDRES1 , NEWRES , NEWRES1
     CHARACTER(LEN=6) :: NMUT_STRING , STARTINDEX_STRING, START_STR, SHIFT_STR
     CHARACTER(LEN=50) :: OPTION_STRING, PATH_STRING
     DOUBLE PRECISION :: COORDINATES(3*NATOMS)
     INTEGER ::  SHIFT, START, J1
 
     !let's store all information first in case we have to go back!
     PREVIOUS_MUTATION = MUTATION_INFO
     !we have a new mutation
     NMUTATION = NMUTATION + 1
     WRITE(NMUT_STRING,'(I6)') NMUTATION - 1 
     !before we do anything, we save the old lowest minima
     CALL AMBERMUT_CURR_LOWEST()
     !select a residue to mutate
     CALL SELECT_MUTATION(RESNUMBER , OLDRES1 , NEWRES1)
     !if it is a terminal residue, we need to go for a different set of atoms and coordinates in the coordinate creation script
     IF (TERMINI_RES(RESNUMBER).EQ.1) THEN
        OLDRES = "C" // OLDRES1
        NEWRES = "C" // NEWRES1
     ELSE IF (TERMINI_RES(RESNUMBER).EQ.2) THEN
        OLDRES = "N" // OLDRES1
        NEWRES = "N" // NEWRES1
     ELSE
        OLDRES = OLDRES1
        NEWRES = NEWRES1
     ENDIF
     WRITE(MUTUNIT,'(A)') 'Currently the best match to the objective function is:'
     WRITE(MUTUNIT,'(A,F20.10)') 'Mutation score: ' , BESTMUTSCORE
     WRITE(MUTUNIT,'(A)',ADVANCE='NO') 'Sequence: '
     DO J1=1,NRESIDUES-1
        WRITE(MUTUNIT,'(A)',ADVANCE='NO') BESTMUTSEQ(J1) // "  "
     ENDDO
     WRITE(MUTUNIT,'(A)') BESTMUTSEQ(NRESIDUES)
     WRITE(MUTUNIT,'(A)') '=============================='
     WRITE(MUTUNIT,'(A,I6,4A)') 'Mutate residue ' , RESNUMBER , ' from ' , OLDRES , ' to ' , NEWRES
     WRITE(STARTINDEX_STRING,'(I6)') AMBER12_RESSTART(RESNUMBER)

     !dump the coordinates for the old residue, and move things to safety
     CALL DUMP_RESIDUE_COORDS(RESNUMBER , COORDINATES)
     CALL SYSTEM('mv coords.prmtop coords.prmtop.'//TRIM(ADJUSTL(NMUT_STRING)))
     CALL SYSTEM('mv coords.inpcrd coords.inpcrd.'//TRIM(ADJUSTL(NMUT_STRING)))
     CALL SYSTEM('mv start start.'//TRIM(ADJUSTL(NMUT_STRING)))
     CALL SYSTEM('mv atomgroups atomgroups.'//TRIM(ADJUSTL(NMUT_STRING)))
     !create mutated coordinates and a new perm.allow file
     OPTION_STRING=OLDRES//' '//NEWRES//' '//STARTINDEX_STRING
#ifdef _SVN_ROOT_
     CALL SYSTEM('python '//_SVN_ROOT_//'/SCRIPTS/AMBER/BHmutation_steps/mutate_aa.py '//OLDRES//' '//NEWRES)
     CALL SYSTEM('python '//_SVN_ROOT_//'/SCRIPTS/AMBER/BHmutation_steps/perm_allow.py '//OPTION_STRING)
#else
     CALL SYSTEM('python ' // mutation_script // OLDRES // ' ' // NEWRES )
     CALL SYSTEM('python ' // perm_allow_script //OPTION_STRING)
#endif
     CALL SYSTEM('mv perm.allow perm.allow.'//TRIM(ADJUSTL(NMUT_STRING)))
     CALL SYSTEM('mv perm.allow.new perm.allow')
     !create a new topology, update the residue information and adjust coordinates for unchanged residues
     CALL CREATE_NEW_TOPOLOGY(RESNUMBER ,  NEWRES , COORDS, START, SHIFT)
     WRITE(START_STR,'(I6)') START
     WRITE(SHIFT_STR,'(I6)') SHIFT
     !create new atom groups
     IF (.NOT.AMBERMUTRIGIDT) THEN
#ifdef _SVN_ROOT_
        CALL SYSTEM('python ' // _SVN_ROOT_ // '/SCRIPTS/AMBER/BHmutation_steps/grouprotations.py tmp.pdb')
#else
        CALL SYSTEM('python ' // grouprotation_script // ' tmp.pdb')
#endif
     ELSE
        OPTION_STRING=' rbodyconfig '//SHIFT_STR//' '//START_STR//' tmp.pdb'
        PATH_STRING='/SCRIPTS/AMBER/BHmutation_steps/rbody_grot.py'
#ifdef _SVN_ROOT_
        
        CALL SYSTEM('python '// _SVN_ROOT_ //PATH_STRING//OPTION_STRING)
        CALL SYSTEM('mv rbodyconfig rbodyconfig.'//TRIM(ADJUSTL(NMUT_STRING)))
#else
        CALL SYSTEM('python ' // grouprotation_script //OPTION_STRING)      
#endif
        CALL SYSTEM('mv rbodyconfig.new rbodyconfig')
        CALL SYSTEM('mv coordsinirigid coordsinirigid.'//TRIM(ADJUSTL(NMUT_STRING)))
        CALL SYSTEM('cp start coordsinirigid')
     ENDIF
     CALL SYSTEM('rm tmp.pdb')
     !finally reinitialise AMBER with new groups, coordinates and topology
     CALL REINITIALISE_AMBER()
     !now remove old chiral states used for checking (the rest is done when we initialise the chirality in mc.F)
     CALL DEALLOC_STATES_MUTATION()
     !finally store the new sequence
     NSEQSTORED = NSEQSTORED + 1
     MUTSEQ(NSEQSTORED, :) = AMBER12_RESNAME(:)
     RETURN     
  END SUBROUTINE AMBERMUT_STEP

  SUBROUTINE SELECT_MUTATION(RESNUMBER , OLDRES , NEWRES)
     INTEGER , INTENT(OUT) :: RESNUMBER
     CHARACTER(LEN=4) , INTENT(OUT) :: OLDRES , NEWRES 
     CHARACTER(LEN=4) :: SELECTED_MUT
     INTEGER :: ENTRIES , NCURR , J1 , SELECTED_ID , SELECTED_RES
     DOUBLE PRECISION :: PROB_RES_SELECT(NRESMUT,2) , NMUTATED , PROB , PROBTOT , RANDOM, DPRAND
     DOUBLE PRECISION , ALLOCATABLE :: PROB_MUT_SELECT(:,:)
     !create probability array to select residue id
     NMUTATED = 0.0
     DO J1 = 1,NRESMUT
        !We take the number of previous mutations plus 1 (otherwise we are at zeros to start with ...)
        NMUTATED = NMUTATED + 1.0/((MUTATION_INFO(J1)%NMUTATIONS) + 1)
     ENDDO
     DO J1 = 1,NRESMUT
        PROB = 1.0/(NMUTATED * ((MUTATION_INFO(J1)%NMUTATIONS) + 1))
        IF (J1.EQ.1) THEN
           !for the first choice we go from zero to prob
           PROB_RES_SELECT(J1,1) = 0.0
           PROB_RES_SELECT(J1,2) = PROB
        ELSE IF (J1.LT.NRESMUT) THEN
           !then we go in intervalls
           PROB_RES_SELECT(J1,1) = PROB_RES_SELECT((J1-1),2)
           PROB_RES_SELECT(J1,2) = PROB_RES_SELECT(J1,1) + PROB
        ELSE
           !finally making sure the array stretches to 1.0
           PROB_RES_SELECT(J1,1) = PROB_RES_SELECT((J1-1),2)
           PROB_RES_SELECT(J1,2) = 1.0
        ENDIF
     ENDDO
     !select residue
     RANDOM=DPRAND()
     DO J1 = 1,NRESMUT
        IF ((PROB_RES_SELECT(J1,1).LT.RANDOM).AND.(RANDOM.LE.PROB_RES_SELECT(J1,2))) THEN
           SELECTED_RES = MUTATION_INFO(J1) % RESNUM
           SELECTED_ID = J1
           WRITE(MYUNIT,'(A,I6)') ' ambermut> Selected residue for mutation: ' , SELECTED_RES
           GOTO 20
        ENDIF
     ENDDO
     !independent of whether we accept or reject the mutation attempt later, we store that it has occured
20   CONTINUE
     MUTATION_INFO(J1)%NMUTATIONS = (MUTATION_INFO(J1)%NMUTATIONS)+1
     PREVIOUS_MUTATION(J1)%NMUTATIONS = (PREVIOUS_MUTATION(J1)%NMUTATIONS)+1 

     !create normalisation for probabilities, same procedure as for the residue
     ENTRIES = MUTATION_INFO(SELECTED_ID)%NENTRIES
     IF (ALLOCATED(PROB_MUT_SELECT)) DEALLOCATE(PROB_MUT_SELECT)
     ALLOCATE(PROB_MUT_SELECT(ENTRIES,2))
     PROB_MUT_SELECT(:,:) = 0.0D0
     PROBTOT = 0.0
     DO J1 = 1,ENTRIES
        IF (.NOT.((MUTATION_INFO(SELECTED_ID)%CURRENT_RES).EQ.(MUTATION_INFO(SELECTED_ID)%RESCHOICE(J1)))) THEN
           PROBTOT = PROBTOT + MUTATION_INFO(SELECTED_ID)%PROBABILITIES(J1)
        ELSE
           NCURR = J1
        ENDIF
     ENDDO
     !create probabilities (making sure we actually mutate)
     DO J1 = 1,ENTRIES
        PROB = (MUTATION_INFO(SELECTED_ID)%PROBABILITIES(J1))/PROBTOT
        IF (J1.EQ.1) THEN
           PROB_MUT_SELECT(J1,1) = 0.0
           IF (J1.EQ.NCURR) THEN
              PROB_MUT_SELECT(J1,2) = 0.0
           ELSE
              PROB_MUT_SELECT(J1,2) = PROB
           ENDIF
        ELSE IF (J1.LT.ENTRIES) THEN
           PROB_MUT_SELECT(J1,1) = PROB_MUT_SELECT((J1-1),2)
           IF (J1.EQ.NCURR) THEN
              PROB_MUT_SELECT(J1,2) = PROB_MUT_SELECT(J1,1)
           ELSE
              PROB_MUT_SELECT(J1,2) = PROB_MUT_SELECT(J1,1) + PROB
           ENDIF
        ELSE
           IF (J1.EQ.NCURR) THEN
              PROB_MUT_SELECT(J1,1) = 1.0
              PROB_MUT_SELECT(J1-1,2) = 1.0
           ELSE
              PROB_MUT_SELECT(J1,1) = PROB_MUT_SELECT(J1-1,2)
           ENDIF
           PROB_MUT_SELECT(J1,2) = 1.0
        ENDIF
     ENDDO
     PROB_MUT_SELECT(NCURR,1) = -1.0
     PROB_MUT_SELECT(NCURR,2) = -1.0
     !select mutation
     RANDOM=DPRAND()
     DO J1 = 1,ENTRIES
     IF ((PROB_MUT_SELECT(J1,1).LT.RANDOM).AND.(RANDOM.LE.PROB_MUT_SELECT(J1,2))) THEN
           SELECTED_MUT = MUTATION_INFO(SELECTED_ID)%RESCHOICE(J1)
           WRITE(MYUNIT,'(A,A)') ' ambermut> Mutate to: ' , SELECTED_MUT
           GOTO 30
        ENDIF
     ENDDO       
30   CONTINUE
     !assign everything to our intent out variables
     RESNUMBER = SELECTED_RES
     OLDRES = MUTATION_INFO(SELECTED_ID)%CURRENT_RES
     NEWRES = SELECTED_MUT  
     MUTATION_INFO(SELECTED_ID)%CURRENT_RES =  NEWRES
  END SUBROUTINE SELECT_MUTATION

  SUBROUTINE DUMP_RESIDUE_COORDS(RESNUMBER , COORD)
     INTEGER , INTENT(IN) :: RESNUMBER
     DOUBLE PRECISION , INTENT(IN) :: COORD(3*NATOMS)
     INTEGER :: STARTATOM , NRESATOM , CUNIT , GETUNIT , J1
     
     !simply dump the coordinates of the residue we want to mutate
     STARTATOM = AMBER12_RESSTART(RESNUMBER)
     NRESATOM = AMBER12_RESNATOM(RESNUMBER)
     CUNIT = GETUNIT()
     OPEN(UNIT=CUNIT , FILE='coords.oldres' , STATUS='NEW')
     DO J1 = 1,NRESATOM
        WRITE(CUNIT,'(3F20.10)') COORD(3*(STARTATOM+J1-1)-2),COORD(3*(STARTATOM+J1-1)-1),COORD(3*(STARTATOM+J1-1))
     ENDDO
     CLOSE(CUNIT)
  END SUBROUTINE DUMP_RESIDUE_COORDS

  SUBROUTINE CREATE_NEW_TOPOLOGY(RESNUMBER , NEWRES , COORDS_OLD, STARTATOM, SHIFT)
     INTEGER , INTENT(IN) :: RESNUMBER
     CHARACTER(LEN=4) , INTENT(IN) :: NEWRES
     DOUBLE PRECISION , INTENT(IN) :: COORDS_OLD(3*NATOMS)
     INTEGER, INTENT(OUT) :: STARTATOM , SHIFT
     DOUBLE PRECISION , ALLOCATABLE :: COORDS_NEW(:) , COORDS_NEWRES(:,:)
     DOUBLE PRECISION , ALLOCATABLE :: COORDS_RES(:)
     INTEGER :: J1 , TUNIT , CUNIT , CUNIT2 , GETUNIT ,  FINALATOM_OLD , FINALATOM_NEW 
     CHARACTER(LEN=4) :: RESNAMES(NRESIDUES)
     CHARACTER(LEN=6) :: CYX_STRING

     TUNIT = GETUNIT()
     DO J1=1,NRESIDUES
        RESNAMES(J1) = AMBER12_RESNAME(J1)
     ENDDO
     !create a leap.in file
     OPEN(TUNIT , FILE='leap.in' , STATUS='NEW')
     !currently we either go for ff14SB or ff99SB
     IF (AMBERMUTFF.EQ.14) THEN
        WRITE(TUNIT,'(A)') 'source leaprc.ff14SB'
     ELSE IF (AMBERMUTFF.EQ.99) THEN
        WRITE(TUNIT,'(A)') 'source oldff/leaprc.ff99SB'
     ENDIF
     !make sure we use the correct adjustment of radii for the solvent model used
     IF (AMBERMUTIGB.EQ.2) THEN
        WRITE(TUNIT,'(A)') 'set default PBradii mbondi2'
     ELSE IF (AMBERMUTIGB.EQ.8) THEN
        WRITE(TUNIT,'(A)') 'set default PBradii mbondi3'
     ENDIF
     !write the sequence including the correct termini (all stored residues have len=3, but the newres is already adjusted to len=4!)
     WRITE(TUNIT,'(A)',ADVANCE='NO') 'mol = sequence {'
     DO J1=1,NRESIDUES
        IF (J1.EQ.RESNUMBER) THEN
           WRITE(TUNIT,'(A)',ADVANCE='NO') NEWRES // " "
        ELSE IF (TERMINI_RES(J1).EQ.2) THEN
           WRITE(TUNIT,'(A)',ADVANCE='NO') "C" // RESNAMES(J1)
        ELSE IF (TERMINI_RES(J1).EQ.1) THEN
           WRITE(TUNIT,'(A)',ADVANCE='NO') "N" // RESNAMES(J1)
        ELSE
           WRITE(TUNIT,'(A)',ADVANCE='NO') RESNAMES(J1)
        ENDIF 
     ENDDO
     WRITE(TUNIT,'(A)') '}'
     !enter bonding for cysteine bonds
     IF (CYX_BONDT) THEN
        DO J1=1,NCYXBONDS
           WRITE(CYX_STRING,'(I6)') CYX_BONDS(J1,1)
           WRITE(TUNIT,'(A)',ADVANCE='NO') "bond mol."//TRIM(ADJUSTL(CYX_STRING))//".SG "
           WRITE(CYX_STRING,'(I6)') CYX_BONDS(J1,2)
           WRITE(TUNIT,'(A)') "mol."//TRIM(ADJUSTL(CYX_STRING))//".SG "
        ENDDO
     ENDIF
     WRITE(TUNIT,'(A)') 'saveamberparm mol coords.prmtop tmp.inpcrd'
     WRITE(TUNIT,'(A)') 'savepdb mol tmp.pdb'
     WRITE(TUNIT,'(A)') 'quit'
     CLOSE(TUNIT)
     !finished creating leap input, now run leap and get the right coordinates     
     CALL SYSTEM('tleap -f leap.in >> output')
     !save the old information
     STARTATOM = AMBER12_RESSTART(RESNUMBER)
     FINALATOM_OLD = AMBER12_RESEND(RESNUMBER)
     CALL TOPOLOGY_READER()
     !get the changed number of atoms
     FINALATOM_NEW = AMBER12_RESEND(RESNUMBER)
     SHIFT = FINALATOM_NEW - FINALATOM_OLD
     !correct wrong information (we havent reinitialised yet, so NATOMS is still wrong)
     AMBER12_RESEND(NRESIDUES) = AMBER12_RESEND(NRESIDUES) + SHIFT
     AMBER12_RESNATOM(NRESIDUES) = AMBER12_RESNATOM(NRESIDUES) + SHIFT
     !create final input files needed
     IF (ALLOCATED(COORDS_NEW)) DEALLOCATE(COORDS_NEW)
     ALLOCATE(COORDS_NEW(3*(NATOMS+SHIFT)))
     IF (ALLOCATED(COORDS_NEWRES)) DEALLOCATE(COORDS_NEWRES)
     ALLOCATE(COORDS_NEWRES(AMBER12_RESNATOM(RESNUMBER),3))
     IF (ALLOCATED(COORDS_RES)) DEALLOCATE(COORDS_RES)
     ALLOCATE(COORDS_RES(3*AMBER12_RESNATOM(RESNUMBER)))
     !fill the new coordinates array 
     COORDS_NEW(:) = 0.0D0
     DO J1 = 1,3*(STARTATOM-1)
        COORDS_NEW(J1) = COORDS_OLD(J1)
     ENDDO

     CUNIT = GETUNIT()
     OPEN(CUNIT , FILE='coords.newres' , STATUS='OLD')
     READ(CUNIT,*) COORDS_NEWRES
     COORDS_RES(:) = RESHAPE(COORDS_NEWRES,(/3*AMBER12_RESNATOM(RESNUMBER)/))
     DO J1 = 1,AMBER12_RESNATOM(RESNUMBER)
        COORDS_NEW(3*(STARTATOM+(J1-1))-2) = COORDS_RES(3*J1-2)
        COORDS_NEW(3*(STARTATOM+(J1-1))-1) = COORDS_RES(3*J1-1)
        COORDS_NEW(3*(STARTATOM+(J1-1))) = COORDS_RES(3*J1)
     ENDDO
40   CLOSE(CUNIT)
   
     DO J1 = 1,(3*(NATOMS-FINALATOM_OLD))
        COORDS_NEW(3*FINALATOM_NEW + J1) = COORDS_OLD(3*FINALATOM_OLD + J1)
     ENDDO
     !we can't write to coords.inpcrd (as this is protected by the interface)
     !hence we trick the program by writing it to a different name and using a system call to move it
     CALL AMBER12_WRITE_RESTART_MUT(COORDS_NEW, AMBER12_RESEND(NRESIDUES),&
                  &'coords.inpcrd.xxx',LEN('coords.inpcrd.xxx'))
     CALL SYSTEM('mv coords.inpcrd.xxx coords.inpcrd')
     CUNIT2 = GETUNIT()
     !create a start file (format specifications are less strict here)
     OPEN(CUNIT2 , FILE='start' , STATUS='NEW')
     DO J1 = 1,NATOMS+SHIFT
        WRITE(CUNIT2 , '(3f12.7)') COORDS_NEW(3*J1-2) , COORDS_NEW(3*J1-1) , COORDS_NEW(3*J1)
     ENDDO
     CLOSE(CUNIT2)
     !finally remove the files we dont need, except if we are in DEBUG mode
     IF (.NOT.DEBUG) CALL SYSTEM('rm coords.newres coords.oldres leap.in leap.log tmp.inpcrd')
  END SUBROUTINE CREATE_NEW_TOPOLOGY

  SUBROUTINE CREATE_RIGID_FILES()
     INTEGER :: J1, J2, GETUNIT, RBCONFUNIT, ID_START, ID_END, ATOMSINGROUP

     RBCONFUNIT = GETUNIT()
     OPEN(RBCONFUNIT , FILE='rbodyconfig' , STATUS='NEW')
     DO J1=1,NRESRB
        ID_START = AMBER12_RESSTART(RESFORRB(J1,1))
        ID_END = AMBER12_RESEND(RESFORRB(J1,2))
        ATOMSINGROUP = ID_END - ID_START + 1
        WRITE(RBCONFUNIT, '(A,I8)') 'GROUP ',ATOMSINGROUP
        DO J2=ID_START,ID_END
           WRITE(RBCONFUNIT,'(I8)') J2
        END DO
     END DO
     CLOSE(RBCONFUNIT)
     CALL SYSTEM('cp start coordsinirigid')
     RETURN
  END SUBROUTINE CREATE_RIGID_FILES

  SUBROUTINE REINITIALISE_AMBER()
        USE GENRIGID, ONLY: GENRIGID_READ_FROM_FILE, DEALLOCATE_GENRIGID
        INTEGER :: NUMBER_OF_ATOMS , J1
        CHARACTER(LEN=20) OSTRING
        DOUBLE PRECISION , ALLOCATABLE :: COORDS1(:)

        !first of all we close all open AMBER files, deallocate all internal arrays, and remove traces from the previous initialisation
        CALL AMBER12_MUT_FINISH()
        !new number of atoms and amber setup
        NUMBER_OF_ATOMS=AMBER12_RESEND(NRESIDUES)
        WRITE(OSTRING,'(A)') 'coords.inpcrd'
        !reinitialise AMBER with the new information
        CALL AMBER12_SETUP(NUMBER_OF_ATOMS, OSTRING, LEN(OSTRING))
        NATOMS = NUMBER_OF_ATOMS
        NATOMSALLOC = NUMBER_OF_ATOMS
        WRITE(MYUNIT,'(A,I8)') ' ambermut> new number of atoms: ',NATOMS
        !new coordinates
        IF(ALLOCATED(COORDS1)) DEALLOCATE(COORDS1)
        ALLOCATE(COORDS1(3*NATOMS))
        IF(ALLOCATED(COORDS)) DEALLOCATE(COORDS)
        ! Read the coords from AMBER12 into COORDS1(:)
        CALL AMBER12_GET_COORDS(NATOMS, COORDS1(:))
        ALLOCATE(COORDS(3*NATOMS,NPAR))
        DO J1=1,NPAR
           COORDS(:,J1) = COORDS1(:)
        END DO
        !setup the new group rotation information
        CALL MUT_SETUP_GROUPROTATION(1,.FALSE.,.FALSE.,0)
        !call reinitialisation for rigid bodies
        IF (AMBERMUTRIGIDT) THEN
           CALL DEALLOCATE_GENRIGID()
           CALL GENRIGID_READ_FROM_FILE()
        ENDIF
        !deallocate, reallocate and initialise a bunch of globals that we need to reset
        DEALLOCATE(QMINP)
        ALLOCATE(QMINP(NSAVE,3*NATOMS))
        DEALLOCATE(QMINT)
        ALLOCATE(QMINT(NSAVE,NATOMS))
        DEALLOCATE(COORDSO)
        ALLOCATE(COORDSO(3*NATOMS,NPAR))
        DEALLOCATE(VT)
        ALLOCATE(VT(NATOMS))
        DEALLOCATE(VAT)
        ALLOCATE(VAT(NATOMS,NPAR))
        DEALLOCATE(VATO)
        ALLOCATE(VATO(NATOMS,NPAR))
        DEALLOCATE(LABELS)
        ALLOCATE(LABELS(NATOMS,NPAR))
        DEALLOCATE(LABELSO)
        ALLOCATE(LABELSO(NATOMS,NPAR))
        QMINP(1:NSAVE,1:3*NATOMS)=0.0D0 ! to prevent reading from uninitialised memory
        QMINT(1:NSAVE,1:NATOMS)=1 ! to prevent reading from uninitialised memory
        QMINNATOMS(1:NSAVE)=NATOMS ! to prevent reading from uninitialised memory
        COORDSO(1:3*NATOMS,1:NPAR)=0.0D0
        VT(1:NATOMS)=0.0D0
        VAT(1:NATOMS,1:NPAR)=0.0D0
        DO J1=1,NSAVE
           QMIN(J1)=1.0D10
           NPCALL_QMIN(J1)=0
        ENDDO
        !this is assuming we do not use frozen atoms, but the genrigid
        !framework
        IF (ALLOCATED(FROZEN)) DEALLOCATE(FROZEN)
        ALLOCATE(FROZEN(NATOMS))
        FROZEN(:) = .FALSE.
  END SUBROUTINE REINITIALISE_AMBER

  SUBROUTINE REVERSE_MUTATION(RESNUMBER)
     CHARACTER(LEN=6) :: NMUT_STRING
     INTEGER :: STARTATOM, FINALATOM_OLD,FINALATOM_NEW,SHIFT
     INTEGER , INTENT(IN) :: RESNUMBER    

     
     !save structures and reload the correct information into MUTATION_INFO
     CALL AMBERMUT_REJ_LOWEST()
     MUTATION_INFO = PREVIOUS_MUTATION 
     WRITE(NMUT_STRING,'(I6)') NMUTATION - 1
     STARTATOM = AMBER12_RESSTART(RESNUMBER)
     FINALATOM_OLD = AMBER12_RESEND(RESNUMBER)
     !move all the files we need back into place (we use the lowest previous minimum to restart)
     CALL SYSTEM('cp coords.prmtop.'//TRIM(ADJUSTL(NMUT_STRING))//' coords.prmtop')
     CALL TOPOLOGY_READER()
     !get the changed number of atoms
     FINALATOM_NEW = AMBER12_RESEND(RESNUMBER)
     SHIFT = FINALATOM_NEW - FINALATOM_OLD
     !correct wrong information (we havent reinitialised yet, so NATOMS is still wrong)
     AMBER12_RESEND(NRESIDUES) = AMBER12_RESEND(NRESIDUES) + SHIFT
     AMBER12_RESNATOM(NRESIDUES) = AMBER12_RESNATOM(NRESIDUES) + SHIFT
     !create final input files needed
     CALL SYSTEM('cp coords.'//TRIM(ADJUSTL(NMUT_STRING))//'.1.rst coords.inpcrd')
     CALL SYSTEM('cp start.'//TRIM(ADJUSTL(NMUT_STRING))//' start')  !this one is the wrong file?
     CALL SYSTEM('cp atomgroups.'//TRIM(ADJUSTL(NMUT_STRING))//' atomgroups')
     CALL SYSTEM('cp perm.allow.'//TRIM(ADJUSTL(NMUT_STRING))//' perm.allow')
     IF (AMBERMUTRIGIDT) THEN
        CALL SYSTEM('cp rbodyconfig.'//TRIM(ADJUSTL(NMUT_STRING))//' rbodyconfig') 
        CALL SYSTEM('cp start coordsinirigid')     
     ENDIF
     !now reinitialise once more
     CALL REINITIALISE_AMBER()
     !reset chirality
     CALL DEALLOC_STATES_MUTATION()
     !store new (old) sequence and note the previous one was rejected
     SEQREJECTED(NSEQSTORED) = .TRUE.
     NSEQSTORED = NSEQSTORED + 1
     MUTSEQ(NSEQSTORED, :) = AMBER12_RESNAME(:)
     RETURN
  END SUBROUTINE REVERSE_MUTATION

  !scoring function to judge how good mutation is
  SUBROUTINE MUTATION_E(SCORE,COORDS,MODE,TERMID)
     DOUBLE PRECISION, INTENT(OUT) :: SCORE
     DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS)
     INTEGER, INTENT(IN) :: MODE, TERMID
     DOUBLE PRECISION :: DPRAND, EREAL, NORMDIFF, ESEP
     DOUBLE PRECISION :: GRADATOMS(3*NATOMS), XSEP(3*NATOMS), XCA_1(3), XCA_2(3), DIFF12(3)
     DOUBLE PRECISION, PARAMETER :: DSEP=1.5D2
     INTEGER :: ATOMID, PARMEDUNIT, GETUNIT, J1, NDUMMY
     TYPE(POT_ENE_REC_C) :: DECOMPOSED_E
     TYPE(AMBER12_ATOM) :: ATOMDATA(NATOMS)
     CHARACTER(200) ENTRY_
     INTEGER , PARAMETER :: NWORDS=20
     CHARACTER(25) :: ENTRIES(NWORDS)=''
     CHARACTER(LEN=6) :: J1_STRING
   
     !random number as penalty function
     IF (MODE.EQ.1) THEN
        SCORE=DPRAND()
     !use contribution form AMBER energies
     ELSE IF (MODE.EQ.2) THEN
        CALL AMBER12_ENERGY_AND_GRADIENT(NATOMS, COORDS, EREAL, GRADATOMS, DECOMPOSED_E)
        WRITE(MUTUNIT,'(A)') 'Energy decomposition'
        WRITE(MUTUNIT,'(A,F20.10)') 'Total energy:        ', DECOMPOSED_E % TOTAL
        WRITE(MUTUNIT,'(A,F20.10)') 'Total van der Waals: ', DECOMPOSED_E % VDW_TOT
        WRITE(MUTUNIT,'(A,F20.10)') 'Total electronic:    ', DECOMPOSED_E % ELEC_TOT
        WRITE(MUTUNIT,'(A,F20.10)') 'Generalised Born:    ', DECOMPOSED_E % GB
        WRITE(MUTUNIT,'(A,F20.10)') 'Surface energy:      ', DECOMPOSED_E % SURF
        WRITE(MUTUNIT,'(A,F20.10)') 'Bond energy:         ', DECOMPOSED_E % BOND
        WRITE(MUTUNIT,'(A,F20.10)') 'Angular term:        ', DECOMPOSED_E % ANGLE
        WRITE(MUTUNIT,'(A,F20.10)') 'Dihedral term:       ', DECOMPOSED_E % DIHEDRAL
        WRITE(MUTUNIT,'(A,F20.10)') 'vdW 1-4 term:        ', DECOMPOSED_E % VDW_14
        WRITE(MUTUNIT,'(A,F20.10)') 'Electronic 1-4:      ', DECOMPOSED_E % ELEC_14
        WRITE(MUTUNIT,'(A,F20.10)') 'Restraints:          ', DECOMPOSED_E % RESTRAINT
        WRITE(MUTUNIT,'(A,F20.10)') 'Urey Bradley angle:  ', DECOMPOSED_E % ANGLE_UB
        WRITE(MUTUNIT,'(A,F20.10)') 'Improper energy:     ', DECOMPOSED_E % IMP
        WRITE(MUTUNIT,'(A,F20.10)') 'CMAP:                ', DECOMPOSED_E % CMAP
        IF (TERMID.EQ.0) THEN
           SCORE = DECOMPOSED_E % TOTAL
        ELSE IF (TERMID.EQ.1) THEN
           SCORE = DECOMPOSED_E % VDW_TOT
        ELSE IF (TERMID.EQ.2) THEN
           SCORE = DECOMPOSED_E % ELEC_TOT
        ELSE IF (TERMID.EQ.3) THEN
           SCORE = DECOMPOSED_E % GB
        ELSE IF (TERMID.EQ.4) THEN
           SCORE = DECOMPOSED_E % SURF
        ELSE IF (TERMID.EQ.5) THEN
           SCORE = DECOMPOSED_E % BOND
        ELSE IF (TERMID.EQ.6) THEN
           SCORE = DECOMPOSED_E % ANGLE
        ELSE IF (TERMID.EQ.7) THEN
           SCORE = DECOMPOSED_E % DIHEDRAL
        ELSE IF (TERMID.EQ.8) THEN
           SCORE = DECOMPOSED_E % VDW_14
        ELSE IF (TERMID.EQ.9) THEN
           SCORE = DECOMPOSED_E % ELEC_14
        ELSE IF (TERMID.EQ.10) THEN
           SCORE = DECOMPOSED_E % RESTRAINT
        ELSE IF (TERMID.EQ.11) THEN
           SCORE = DECOMPOSED_E % ANGLE_UB
        ELSE IF (TERMID.EQ.12) THEN
           SCORE = DECOMPOSED_E % IMP
        ELSE IF (TERMID.EQ.13) THEN
           SCORE = DECOMPOSED_E % CMAP
        ENDIF
     !approximate brinding interactions
     ELSE IF (MODE.EQ.3) THEN
        ATOMID = AMBER12_RESEND(TERMID) !get last atom in first group
        !old way - very slow as it uses external script
        !open new file to write parmed input
        !PARMEDUNIT = GETUNIT()
        !OPEN(PARMEDUNIT,FILE='parmed_in',STATUS='NEW')
        !WRITE(PARMEDUNIT,'(A)',ADVANCE='NO') 'addExclusions @1'        
        !DO J1=2,ATOMID
        !   WRITE(J1_STRING,'(I6)') J1
        !   WRITE(PARMEDUNIT,'(A)',ADVANCE='NO') ','//TRIM(ADJUSTL(J1_STRING))
        !ENDDO
        !WRITE(PARMEDUNIT,'(A)',ADVANCE='NO') ' @1'       
        !DO J1=2,ATOMID-1
        !   WRITE(J1_STRING,'(I6)') J1
        !   WRITE(PARMEDUNIT,'(A)',ADVANCE='NO') ','//TRIM(ADJUSTL(J1_STRING))
        !ENDDO
        !WRITE(J1_STRING,'(I6)') ATOMID
        !WRITE(PARMEDUNIT,'(A)') ','//TRIM(ADJUSTL(J1_STRING))
        !WRITE(J1_STRING,'(I6)') ATOMID+1
        !WRITE(PARMEDUNIT,'(A)',ADVANCE='NO') 'addExclusions @'//TRIM(ADJUSTL(J1_STRING))        
        !DO J1=ATOMID+2,NATOMS
        !   WRITE(J1_STRING,'(I6)') J1
        !   WRITE(PARMEDUNIT,'(A)',ADVANCE='NO') ','//TRIM(ADJUSTL(J1_STRING))
        !ENDDO
        !WRITE(J1_STRING,'(I6)') ATOMID+1
        !WRITE(PARMEDUNIT,'(A)',ADVANCE='NO') ' @'//TRIM(ADJUSTL(J1_STRING))       
        !DO J1=ATOMID+2,NATOMS-1
        !   WRITE(J1_STRING,'(I6)') J1
        !   WRITE(PARMEDUNIT,'(A)',ADVANCE='NO') ','//TRIM(ADJUSTL(J1_STRING))
        !ENDDO
        !WRITE(J1_STRING,'(I6)') NATOMS-1
        !WRITE(PARMEDUNIT,'(A)') ','//TRIM(ADJUSTL(J1_STRING))
        !WRITE(PARMEDUNIT,'(A)') 'loadRestrt current.inpcrd'
        !WRITE(J1_STRING,'(I6)') AMBERMUTIGB
        !WRITE(PARMEDUNIT,'(A)') 'energy cutoff 15.0 igb '//TRIM(ADJUSTL(J1_STRING))//' saltcon 0.1'
        !WRITE(PARMEDUNIT,'(A)') 'quit'
        !CLOSE(PARMEDUNIT)
        !CALL AMBER12_WRITE_RESTART(COORDS, 'current.inpcrd',LEN('current.inpcrd'))
        !!create new topology without interactions and calculate energy
        !CALL SYSTEM('parmed.py -n coords.prmtop parmed_in > parmed_out')
        !OPEN(PARMEDUNIT,FILE='parmed_out',STATUS='OLD')
        !DO
        !  ENTRIES(:)=''
        !  READ(PARMEDUNIT,'(A)',END=588) ENTRY_
        !  CALL READ_LINE(ENTRY_,NWORDS,ENTRIES)
        !  IF (ENTRIES(1).EQ.'TOTAL') THEN
        !     READ(ENTRIES(3),'(F20.10)') SCORE
        !     GOTO 588
        !  ENDIF
        !ENDDO
!588     CONTINUE
        !CLOSE(PARMEDUNIT)
        !CALL SYSTEM('rm current.inpcrd parmed_in parmed_out')


        !new way - simply move first and second half apart and compute amber energy
        CALL AMBER12_GET_ATOMDATA(ATOMDATA,NATOMS)
        !get the positions of the C alpha atoms and then assign the centre of these coordinates
        NDUMMY = 0
        XCA_1(:) = 0.0D0
        XCA_2(:) = 0.0D0 
        DO J1=1,NATOMS
            IF (ATOMDATA(J1)%NAME.EQ.'C') THEN
                NDUMMY = NDUMMY + 1
                IF (J1.LE.ATOMID) THEN
                    XCA_1(1) = XCA_1(1) + COORDS(3*J1-2)
                    XCA_1(2) = XCA_1(2) + COORDS(3*J1-1)
                    XCA_1(3) = XCA_1(3) + COORDS(3*J1)
                    IF (NDUMMY.EQ.TERMID) THEN
                        XCA_1(1) = XCA_1(1)/TERMID
                        XCA_1(2) = XCA_1(2)/TERMID
                        XCA_1(3) = XCA_1(3)/TERMID
                    ENDIF
                ELSE
                    XCA_2(1) = XCA_2(1) + COORDS(3*J1-2)
                    XCA_2(2) = XCA_2(2) + COORDS(3*J1-1)
                    XCA_2(3) = XCA_2(3) + COORDS(3*J1)
                ENDIF
            ENDIF
            IF (NDUMMY.EQ.NRESIDUES) GOTO 760
        END DO
760     CONTINUE
        XCA_1(1) = XCA_1(1)/(NRESIDUES-TERMID)
        XCA_1(2) = XCA_1(2)/(NRESIDUES-TERMID)
        XCA_1(3) = XCA_1(3)/(NRESIDUES-TERMID)
        !now define the vector between the two centres
        DIFF12(1) = XCA_2(1) - XCA_1(1)
        DIFF12(2) = XCA_2(2) - XCA_1(2)
        DIFF12(3) = XCA_2(3) - XCA_1(3)
        NORMDIFF = 1.0D0/SQRT(DIFF12(1)**2 + DIFF12(2)**2 + DIFF12(3)**2)
        DIFF12(1) = NORMDIFF * DIFF12(1)
        DIFF12(2) = NORMDIFF * DIFF12(2)
        DIFF12(3) = NORMDIFF * DIFF12(3)
        !now we have the normalised vector pointing along the Calpha centres
        XSEP(:) = COORDS(:)
        !first recentre coordinates
        DO J1=1,ATOMID
           XSEP(3*J1-2) = COORDS(3*J1-2) - XCA_1(1)
           XSEP(3*J1-1) = COORDS(3*J1-1) - XCA_1(2)
           XSEP(3*J1)   = COORDS(3*J1)   - XCA_1(3)
        END DO
        DO J1=(ATOMID+1),NATOMS
           XSEP(3*J1-2) = COORDS(3*J1-2) - XCA_1(1) + DSEP * DIFF12(1)
           XSEP(3*J1-1) = COORDS(3*J1-1) - XCA_1(2) + DSEP * DIFF12(2)
           XSEP(3*J1)   = COORDS(3*J1)   - XCA_1(3) + DSEP * DIFF12(3)
        END DO
        CALL AMBER12_ENERGY_AND_GRADIENT(NATOMS, XSEP, ESEP, GRADATOMS, DECOMPOSED_E)
        CALL AMBER12_ENERGY_AND_GRADIENT(NATOMS, COORDS, EREAL, GRADATOMS, DECOMPOSED_E)
        SCORE = EREAL - ESEP
     ELSE IF (MODE.EQ.4) THEN
         IF (.NOT.(ALLOCATED(CA_REFERENCE))) THEN
            ALLOCATE(CA_REFERENCE(3*NRESIDUES))
            CALL CREATE_CA_REF(TERMID)
         END IF
         CALL CALPHA_RMSD(SCORE,COORDS)
     !helical optimisation for alpha helix
     ELSE IF (MODE.EQ.5) THEN
         CALL HELICAL_DSSP(SCORE,COORDS,1)
     !helical optimisation for 3-10 helix
     ELSE IF (MODE.EQ.6) THEN
         CALL HELICAL_DSSP(SCORE,COORDS,2)
     !helical optimisation for pi helix
     ELSE IF (MODE.EQ.7) THEN
         CALL HELICAL_DSSP(SCORE,COORDS,3)
     ENDIF
  END SUBROUTINE MUTATION_E

  SUBROUTINE AMBERMUTDUMP(MUTATEDT)
     LOGICAL, INTENT(IN) :: MUTATEDT
     INTEGER :: MIUNIT, GETUNIT, NTERMINI, J1, NDUMMY

     !save a copy of all the restart files needed
     CALL SYSTEM('cp coords.prmtop DUMP.coords.prmtop')
     CALL SYSTEM('cp coords.inpcrd DUMP.coords.inpcrd')
     CALL SYSTEM('cp start DUMP.start')
     CALL SYSTEM('cp atomgroups DUMP.atomgroups')
     CALL SYSTEM('cp perm.allow DUMP.perm.allow')
     IF (AMBERMUTRIGIDT) THEN
        CALL SYSTEM('cp rbodyconfig DUMP.rbodyconfig') 
        CALL SYSTEM('cp coordsinirigid DUMP.coordsinirigid')     
     ENDIF
     !create amber_mutation file
     MIUNIT = GETUNIT()
     OPEN(MIUNIT,FILE='DUMP.amber_mutations',STATUS='UNKNOWN')
     !first line is the number of residues to mutate
     WRITE(MIUNIT,'(I8)') NRESMUT
     !number of termini
     NTERMINI = 0
     DO J1=1,NRESIDUES
        IF ((TERMINI_RES(J1).EQ.1).OR.(TERMINI_RES(J1).EQ.2)) NTERMINI = NTERMINI + 1
     ENDDO
     WRITE(MIUNIT,'(I8)') NTERMINI
     !list of termini id
     NDUMMY = 0
     DO J1=1,NRESIDUES
        IF ((TERMINI_RES(J1).EQ.1).OR.(TERMINI_RES(J1).EQ.2)) THEN
           NDUMMY = NDUMMY + 1
           IF (NDUMMY.EQ.NTERMINI) THEN
              WRITE(MIUNIT,'(I8)') J1
           ELSE
              WRITE(MIUNIT,'(I8)',ADVANCE='NO') J1
           ENDIF
        ENDIF
     ENDDO
     !now we get to actual mutation information
     !line 1: RESNUM NENTRIES CURRENT_RES NMUTATIONS
     !line 2: RESNAME1 RESNAME2 RESNAME3 ...
     !line 3: PROB1 PROB2 PROB3 ...
     DO J1=1,NRESMUT
        WRITE(MIUNIT,'(2(I8),A6,I8)') MUTATION_INFO(J1)%RESNUM, MUTATION_INFO(J1)%NENTRIES, &
       &                              MUTATION_INFO(J1)%CURRENT_RES, MUTATION_INFO(J1)%NMUTATIONS
        WRITE(MIUNIT, *) MUTATION_INFO(J1)%RESCHOICE(:)
        WRITE(MIUNIT, *) MUTATION_INFO(J1)%PROBABILITIES(:)       
     ENDDO
     CLOSE(MIUNIT)
  END SUBROUTINE AMBERMUTDUMP

  SUBROUTINE AMBERMUT_CURR_LOWEST()
     CHARACTER(LEN=6) :: J1_STRING, NMUT_STRING
     INTEGER :: J1    
 
     !to save the previous best strutures 
     DO J1=1,NSAVE
        WRITE(J1_STRING,'(I6)') J1
        WRITE(NMUT_STRING,'(I6)') NMUTATION - 1
        CALL AMBER12_WRITE_RESTART(QMINP(J1,:), 'coords.'//TRIM(ADJUSTL(NMUT_STRING))//'.'//&
                   &TRIM(ADJUSTL(J1_STRING))//'.rst', &
                   & LEN('coords.'//TRIM(ADJUSTL(NMUT_STRING))//'.'//TRIM(ADJUSTL(J1_STRING))//'.rst'))
        CALL AMBER12_WRITE_PDB(QMINP(J1,:), 'coords.'//TRIM(ADJUSTL(NMUT_STRING))//'.'//&
                   &TRIM(ADJUSTL(J1_STRING))//'.pdb', &
                   & LEN('coords.'//TRIM(ADJUSTL(NMUT_STRING))//'.'//TRIM(ADJUSTL(J1_STRING))//'.pdb'))
     ENDDO
  END SUBROUTINE AMBERMUT_CURR_LOWEST

  SUBROUTINE AMBERMUT_REJ_LOWEST()
     CHARACTER(LEN=6) :: J1_STRING, NMUT_STRING
     INTEGER :: J1

     !to save the previous best strutures 
     DO J1=1,NSAVE
        WRITE(J1_STRING,'(I6)') J1
        WRITE(NMUT_STRING,'(I6)') NMUTATION - 1
        CALL AMBER12_WRITE_RESTART(QMINP(J1,:), 'rej.coords.'//TRIM(ADJUSTL(NMUT_STRING))//'.'//&
                   &TRIM(ADJUSTL(J1_STRING))//'.rst', &
                   & LEN('rej.coords.'//TRIM(ADJUSTL(NMUT_STRING))//'.'//TRIM(ADJUSTL(J1_STRING))//'.rst'))
        CALL AMBER12_WRITE_PDB(QMINP(J1,:), 'rej.coords.'//TRIM(ADJUSTL(NMUT_STRING))//'.'//&
                   &TRIM(ADJUSTL(J1_STRING))//'.pdb', &
                   & LEN('rej.coords.'//TRIM(ADJUSTL(NMUT_STRING))//'.'//TRIM(ADJUSTL(J1_STRING))//'.pdb'))
     ENDDO
  END SUBROUTINE AMBERMUT_REJ_LOWEST

  
  !tidy up after run is complete
  SUBROUTINE FINISH_AMBERMUT()
     CHARACTER(LEN=6) :: J1_STRING, NMUT_STRING
     INTEGER :: J1    
 
     !to save the previous best strutures 
     DO J1=1,NSAVE
        WRITE(J1_STRING,'(I6)') J1
        WRITE(NMUT_STRING,'(I6)') NMUTATION
        CALL AMBER12_WRITE_RESTART(QMINP(J1,:), 'coords.'//TRIM(ADJUSTL(NMUT_STRING))//'.'//&
                   &TRIM(ADJUSTL(J1_STRING))//'.rst', &
                   & LEN('coords.'//TRIM(ADJUSTL(NMUT_STRING))//'.'//TRIM(ADJUSTL(J1_STRING))//'.rst'))
        CALL AMBER12_WRITE_PDB(QMINP(J1,:), 'coords.'//TRIM(ADJUSTL(NMUT_STRING))//'.'//&
                   &TRIM(ADJUSTL(J1_STRING))//'.pdb', &
                   & LEN('coords.'//TRIM(ADJUSTL(NMUT_STRING))//'.'//TRIM(ADJUSTL(J1_STRING))//'.pdb'))
     ENDDO

     IF (ALLOCATED(AMBER12_RESNAME)) DEALLOCATE(AMBER12_RESNAME)
     IF (ALLOCATED(AMBER12_RESSTART)) DEALLOCATE(AMBER12_RESSTART)
     IF (ALLOCATED(AMBER12_RESEND)) DEALLOCATE(AMBER12_RESEND)
     IF (ALLOCATED(AMBER12_RESNATOM)) DEALLOCATE(AMBER12_RESNATOM)
     IF (ALLOCATED(TERMINI_RES)) DEALLOCATE(TERMINI_RES)
     IF (ALLOCATED(MUTATION_INFO)) DEALLOCATE(MUTATION_INFO)
     IF (ALLOCATED(PREVIOUS_MUTATION)) DEALLOCATE(PREVIOUS_MUTATION)
     CLOSE(MUTUNIT)
  END SUBROUTINE FINISH_AMBERMUT

  SUBROUTINE PRINT_CURRENT_SEQ()
     INTEGER J
     
     DO J=1,NRESIDUES-1
        WRITE(MYUNIT, '(A)',ADVANCE='NO') AMBER12_RESNAME(J)
     END DO
     WRITE(MYUNIT, '(A)')  AMBER12_RESNAME(NRESIDUES)
  END SUBROUTINE PRINT_CURRENT_SEQ

  SUBROUTINE CREATE_CA_REF(HELTYPE)
    INTEGER, INTENT(IN) :: HELTYPE
    DOUBLE PRECISION, PARAMETER :: HBOND=5.40D0
    DOUBLE PRECISION, PARAMETER :: PI=3.14159265358979D0
    DOUBLE PRECISION :: CA_DIST, CA_AXIS, CA_ANGLE, XYZ_OLD(3), XYZ_NEW(3), RMAT(3,3)
    INTEGER :: J1, REFUNIT, GETUNIT

    IF (HELTYPE.EQ.1) THEN !alpha helix
      CA_ANGLE = 100.0D0
    ELSE IF (HELTYPE.EQ.2) THEN !3-10 helix
      CA_ANGLE = 120.0D0
    ELSE IF (HELTYPE.EQ.3) THEN !pi helix
      CA_ANGLE = 87.0D0
    END IF
    CA_DIST = 2.5D0
    CA_AXIS = HBOND/(360.0D0/CA_ANGLE)
    CA_ANGLE = (CA_ANGLE/180.0D0)*PI
    !set firsat atom on x-axis, and have helical axis along z-axis, and set RMAT 
    CA_REFERENCE(1) = CA_DIST
    CA_REFERENCE(2) = 0.0D0
    CA_REFERENCE(3) = 0.0D0
    RMAT = RESHAPE((/ COS(CA_ANGLE),SIN(CA_ANGLE),0.0D0,-SIN(CA_ANGLE),COS(CA_ANGLE),0.0D0,0.0D0,0.0D0,1.0D0 /),(/ 3,3 /))
    !now create the new coordinates
    DO J1=2,NRESIDUES
      XYZ_OLD=(/ CA_REFERENCE(3*(J1-1)-2),CA_REFERENCE(3*(J1-1)-1),CA_REFERENCE(3*(J1-1)) /)
      XYZ_NEW=MATMUL(RMAT,XYZ_OLD)
      CA_REFERENCE(3*J1-2) = XYZ_NEW(1) 
      CA_REFERENCE(3*J1-1) = XYZ_NEW(2)
      CA_REFERENCE(3*J1)   = XYZ_NEW(3) + CA_AXIS
    ENDDO
    REFUNIT = GETUNIT()
    OPEN(UNIT=REFUNIT , FILE='reference.CA' , STATUS='NEW')
     DO J1 = 1,NRESIDUES
        WRITE(REFUNIT,'(3F20.10)') CA_REFERENCE(3*J1-2),CA_REFERENCE(3*J1-1),CA_REFERENCE(3*J1)
     ENDDO
     CLOSE(REFUNIT)

  END SUBROUTINE CREATE_CA_REF

  SUBROUTINE CALPHA_RMSD(SCORE,COORDS)
    DOUBLE PRECISION, INTENT(OUT) :: SCORE
    DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS)
    DOUBLE PRECISION :: CA_COORDS(3*NRESIDUES), DIST, RMAT(3,3)
    TYPE(AMBER12_ATOM) :: ATOMDATA(NATOMS)
    INTEGER :: J1, NDUMMY, CA_POS(NRESIDUES)
 
    CALL AMBER12_GET_ATOMDATA(ATOMDATA,NATOMS)
    !get the positions of the C alpha atoms (as they change with the mutations this is the easiest way)
    NDUMMY = 0
    DO J1=1,NATOMS
       IF (ATOMDATA(J1)%NAME.EQ.'C') THEN
          NDUMMY = NDUMMY + 1
          CA_POS(NDUMMY) = J1
       END IF
       IF (NDUMMY.EQ.NRESIDUES) GOTO 713
    END DO
    !now create a new array of coordinates just containing the C alpha coordinates
713 CONTINUE
    DO J1=1,NRESIDUES
       CA_COORDS(3*J1-2) = COORDS(3*CA_POS(J1)-2)
       CA_COORDS(3*J1-1) = COORDS(3*CA_POS(J1)-1)
       CA_COORDS(3*J1)   = COORDS(3*CA_POS(J1))
    END DO
    !now compare CA_COORDS and CA_REFERENCE
    CALL NEWMINDIST(CA_REFERENCE,CA_COORDS,NRESIDUES,DIST,.FALSE.,.FALSE.,'     ',.FALSE.,.FALSE.,DEBUG,RMAT)
    SCORE = DIST
  END SUBROUTINE CALPHA_RMSD

  SUBROUTINE HELICAL_DSSP(SCORE,COORDS,MODE)
     INTEGER, INTENT(IN) :: MODE                     !1 - alpha helix, 2 - 3-10 helix, 3 -pi
     DOUBLE PRECISION, INTENT(IN) :: COORDS(3*NATOMS)
     DOUBLE PRECISION, INTENT(OUT) :: SCORE
     LOGICAL :: HBOND(NRESIDUES,NRESIDUES)
     TYPE(AMBER12_ATOM) :: ATOMDATA(NATOMS)
     TYPE(AMBER12_RESIDUE) :: RESDATA(NRESIDUES)
     INTEGER :: J1, J2, J3, N_ID, H_ID, C_ID, O_ID
     DOUBLE PRECISION :: R_ON, R_CH, R_OH, R_CN, NX(3), HX(3), CX(3), OX(3), EHBOND
     INTEGER :: NIDEAL, NACTUAL, NOTHER

     HBOND(:,:) = .FALSE.
     CALL AMBER12_GET_ATOMDATA(ATOMDATA,NATOMS)
     CALL AMBER12_GET_RESDATA(RESDATA,NRESIDUES)
     !test for the existance of H bonds
     !ignore terminal residues NH3+ and COO-
     DO J1=2,NRESIDUES           !N-H
       DO J2=1,NRESIDUES-1       !C-O
         IF (.NOT.(J1.EQ.J2)) THEN
           DO J3=RESDATA(J1)%START_INDEX,RESDATA(J1)%END_INDEX
              IF (ATOMDATA(J3)%NAME.EQ.'N') N_ID = J3
              IF (ATOMDATA(J3)%NAME.EQ.'H') H_ID = J3
           END DO
           DO J3=RESDATA(J1)%START_INDEX,RESDATA(J1)%END_INDEX
              IF (ATOMDATA(J3)%NAME.EQ.'C') C_ID = J3
              IF (ATOMDATA(J3)%NAME.EQ.'O') O_ID = J3
           END DO
           NX = (/ COORDS(3*N_ID-2),COORDS(3*N_ID-1),COORDS(3*N_ID) /)
           HX = (/ COORDS(3*H_ID-2),COORDS(3*H_ID-1),COORDS(3*H_ID) /)
           CX = (/ COORDS(3*C_ID-2),COORDS(3*C_ID-1),COORDS(3*C_ID) /)
           OX = (/ COORDS(3*O_ID-2),COORDS(3*O_ID-1),COORDS(3*O_ID) /)
           CALL DISTANCE(OX,NX,R_ON)
           CALL DISTANCE(CX,HX,R_CH)
           CALL DISTANCE(OX,HX,R_OH)
           CALL DISTANCE(CX,NX,R_CN)
           EHBOND = 0.084*332*(1/R_ON + 1/R_CH - 1/R_OH - 1/R_CH)
           IF (EHBOND.LT.-0.5) HBOND(J1,J2) = .TRUE.
         ENDIF
       ENDDO
     ENDDO
     !check all entries to identify secondary structure patterns
     NIDEAL = 0
     NACTUAL = 0
     NOTHER = 0
     IF (MODE.EQ.1) THEN
        NIDEAL = NRESIDUES - 4
     ELSE IF (MODE.EQ.2) THEN
        NIDEAL = NRESIDUES - 3
     ELSE IF (MODE.EQ.3) THEN
        NIDEAL = NRESIDUES - 5
     ENDIF
     DO J1=1,NRESIDUES
       DO J2=1,NRESIDUES
         IF (HBOND(J1,J2)) THEN
           IF (MODE.EQ.1) THEN
             IF ((J1-4).EQ.J2) THEN
                NACTUAL = NACTUAL + 1
             ELSE
                NOTHER = NOTHER + 1
             END IF
           ELSE IF (MODE.EQ.2) THEN
             IF ((J1-3).EQ.J2) THEN
                NACTUAL = NACTUAL + 1
             ELSE
                NOTHER = NOTHER + 1
             END IF
           ELSE IF (MODE.EQ.3) THEN
             IF ((J1-5).EQ.J2) THEN
                NACTUAL = NACTUAL + 1
             ELSE
                NOTHER = NOTHER
             END IF
           END IF
         END IF
       ENDDO
     ENDDO
     SCORE = (NACTUAL - NOTHER)/NIDEAL 
  END SUBROUTINE HELICAL_DSSP

  SUBROUTINE MUT_SETUP_GROUPROTATION(GROUPROTFREQ,GR_SCALEROT,GR_SCALEPROB,GROUPOFFSET)
     INTEGER, INTENT(IN) :: GROUPROTFREQ , GROUPOFFSET
     LOGICAL, INTENT(IN) :: GR_SCALEROT , GR_SCALEPROB
     INTEGER ::  GROUPSIZE , GROUPATOM , AXIS1 , AXIS2 , IOSTATUS, J1,J2
     CHARACTER(LEN=10) :: CHECK1
     LOGICAL :: YESNO

     !check we actually have a grouprotation file!
     YESNO=.FALSE.
     INQUIRE(FILE='atomgroups',EXIST=YESNO)
     IF (YESNO) THEN
        GROUPROTT=.TRUE.
        WRITE(MYUNIT,'(A)') ' ambermut> AMBER group rotation moves enabled for new sequence'
     ELSE
        WRITE(MYUNIT,'(A)') ' keyword> ERROR: atom groups must be defined in atomgroups file'
        STOP
     ENDIF
     !check the grouprotation frequency
     IF(GROUPROTFREQ.EQ.0) THEN
        GROUPROTT=.FALSE.
        WRITE(MYUNIT,'(A)') ' keyword> WARNING: frequency of GROUPROTATION moves set to 0 - moves DISABLED!'
     ENDIF
     !kr366> copy ffrom keywords.f
     !csw34> Figure out how many atom groups have been defined
     NGROUPS=0
     OPEN(UNIT=222,FILE='atomgroups',status='old')
     DO
        READ(222,*,IOSTAT=iostatus) CHECK1
        IF (iostatus<0) THEN
           CLOSE(222)
           EXIT
        ELSE IF (TRIM(ADJUSTL(check1)).EQ.'GROUP') then
           NGROUPS=NGROUPS+1
        ENDIF
     END DO
     CLOSE(222)
     !DEALLOCATE old arrays first
     DEALLOCATE(ATOMGROUPNAMES)
     DEALLOCATE(ATOMGROUPAXIS)
     DEALLOCATE(ATOMGROUPPSELECT)
     DEALLOCATE(ATOMGROUPSCALING)
     DEALLOCATE(ATOMGROUPS)
     !Allocate atom group info arrays appropriately
     ALLOCATE(ATOMGROUPNAMES(NGROUPS))
     ALLOCATE(ATOMGROUPAXIS(NGROUPS,2))
     ALLOCATE(ATOMGROUPPSELECT(NGROUPS))
     ALLOCATE(ATOMGROUPSCALING(NGROUPS))
     ALLOCATE(ATOMGROUPS(NGROUPS,NATOMSALLOC))
     !Set safe defaults
     ATOMGROUPS(:,:)=.FALSE.
     ATOMGROUPNAMES(:)='EMPTY'
     ATOMGROUPAXIS(:,:)=0
     ATOMGROUPSCALING(:)=1.0D0
     ATOMGROUPPSELECT(:)=1.0D0
     ! Read in group info
     ! Here is an example entry:
     ! GROUP OME 6 5 4 1.0
     ! 1
     ! 2
     ! 3
     ! 4
     ! This says that group OME is to be rotated about the bond from atom 6->5.
     ! There are 4 atoms in the OME group. Rotations of -pi->+pi are to be scaled by 1.0.
     ! Finally, the group members are specified one per line
     OPEN(UNIT=222,FILE='atomgroups',status='unknown')
     WRITE(MYUNIT,*) 'keyword> Reading in atom groups for GROUPROTATION'
     IF(GROUPOFFSET.NE.0) WRITE(MYUNIT,*) 'keyword> Group atom numbering offset by ',GROUPOFFSET
     DO J1=1,NGROUPS
        READ(222,*) CHECK1,ATOMGROUPNAMES(J1),AXIS1,AXIS2,GROUPSIZE,ATOMGROUPSCALING(J1),&
     &                ATOMGROUPPSELECT(J1)
        ATOMGROUPAXIS(J1,1)=AXIS1+GROUPOFFSET
        ATOMGROUPAXIS(J1,2)=AXIS2+GROUPOFFSET
        CALL FLUSH(MYUNIT)
        IF (TRIM(ADJUSTL(CHECK1)).EQ.'GROUP') THEN
           DO J2=1,GROUPSIZE
              READ(222,*) GROUPATOM
              IF(GROUPOFFSET.GT.0) GROUPATOM=GROUPATOM+GROUPOFFSET
              !add bound checks
              IF (GROUPATOM > NATOMSALLOC) THEN
                WRITE(MYUNIT,'(A)') 'ambermut> ERROR! GROUPATOM > NATOMSALLOC'
              ENDIF
              ATOMGROUPS(J1,GROUPATOM)=.TRUE.
           END DO
        ELSE
           WRITE(MYUNIT,'(A)') ' keyword: ERROR! Group file not formatted correctly!'
           STOP
        ENDIF
        WRITE(MYUNIT,'(3A)') '<GROUP ',TRIM(ADJUSTL(ATOMGROUPNAMES(J1))),'>'
        WRITE(MYUNIT,'(A,I3)') 'Index: ',J1
        WRITE(MYUNIT,'(A,I4)') 'Size: ',GROUPSIZE
        WRITE(MYUNIT,'(A,2I6)') 'Atoms defining axis: ',ATOMGROUPAXIS(J1,1),ATOMGROUPAXIS(J1,2)
        WRITE(MYUNIT,'(A,F5.2)') 'Rotation scaling: ',ATOMGROUPSCALING(J1)
        WRITE(MYUNIT,'(A,F5.4)') 'Selection probablity: ',ATOMGROUPPSELECT(J1)
        WRITE(MYUNIT,'(A)') 'Members:'
        DO J2=1,NATOMSALLOC
           IF(ATOMGROUPS(J1,J2)) WRITE(MYUNIT,*) J2
        ENDDO
     ENDDO
     CLOSE(222)
  END SUBROUTINE MUT_SETUP_GROUPROTATION 

  SUBROUTINE TOPOLOGY_READER()

     IMPLICIT NONE
     CHARACTER(200) ENTRY_
     INTEGER :: MYUNIT2,GETUNIT
     INTEGER :: HENTRIES,J3,J4,NDUMMY
     INTEGER , PARAMETER :: NWORDS=20
     CHARACTER(25) :: ENTRIES(NWORDS)=''
     CHARACTER(LEN=4) :: WORD

     MYUNIT2=GETUNIT()
     OPEN(MYUNIT2,FILE='coords.prmtop',STATUS='OLD')
     reading:DO
98      ENTRIES(:)=''
        READ(MYUNIT2,'(A)',END=99) ENTRY_
        CALL READ_LINE(ENTRY_,NWORDS,ENTRIES)      !get all words in line
        IF (ENTRIES(2).EQ.'POINTERS') THEN        !get number of residues here
           READ(MYUNIT2,*)                          !ignore format identifier after flag
           READ(MYUNIT2,*)                          !ignore first line, no information we need in here
           READ(MYUNIT2,'(A)',END=99) ENTRY_
           ENTRIES(:)=''
           CALL READ_LINE(ENTRY_,NWORDS,ENTRIES)
           READ(ENTRIES(2),'(I8)') NRESIDUES
           WRITE(MYUNIT,'(A,I8)') 'ambermut> reading topology - Number of residues:' , NRESIDUES
           IF (ALLOCATED(AMBER12_RESNAME)) DEALLOCATE(AMBER12_RESNAME)
           ALLOCATE(AMBER12_RESNAME(NRESIDUES))
           AMBER12_RESNAME(:) = "    "
           IF (ALLOCATED(AMBER12_RESSTART)) DEALLOCATE(AMBER12_RESSTART)
           ALLOCATE(AMBER12_RESSTART(NRESIDUES))
           AMBER12_RESSTART(:) = 0
           IF (ALLOCATED(AMBER12_RESEND)) DEALLOCATE(AMBER12_RESEND)
           ALLOCATE(AMBER12_RESEND(NRESIDUES))
           AMBER12_RESEND(:) = 0
           IF (ALLOCATED(AMBER12_RESNATOM)) DEALLOCATE(AMBER12_RESNATOM)
           ALLOCATE(AMBER12_RESNATOM(NRESIDUES))
           AMBER12_RESNATOM(:) = 0
        ENDIF
        IF (ENTRIES(2).EQ. 'RESIDUE_LABEL') THEN
           READ(MYUNIT2,*)                        !ignore format identifier after flag
           IF (MOD(NRESIDUES,20).EQ.0) THEN       !get the number of lines (20 entries per line!)
              HENTRIES=NRESIDUES/20
           ELSE
              HENTRIES=NRESIDUES/20 + 1
           ENDIF
           !We leave th complication of terminal residues out here and take care of it in the atom mapping when taking a step 
           NDUMMY=1
           DO J3=1,HENTRIES                             !go through all lines
              READ(MYUNIT2,'(A)',END=99) ENTRY_               !read line
              ENTRIES(:)=''
              CALL READ_LINE(ENTRY_,NWORDS,ENTRIES)
              J4=1
              DO WHILE(J4.LE.20)
                 IF (NDUMMY.LE.NRESIDUES) THEN
                    WORD = ENTRIES(J4)(1:3)
                    AMBER12_RESNAME(NDUMMY) = WORD
                    NDUMMY = NDUMMY + 1
                 ELSE
                    GOTO 98
                 ENDIF
                 J4=J4+1
              ENDDO
           ENDDO
        ENDIF
        IF (ENTRIES(2).EQ. 'RESIDUE_POINTER') THEN
           READ(MYUNIT2,*)                             !ignore format identifier after flag
           IF (MOD(NRESIDUES,10).EQ.0) THEN       !get the number of lines (10 entries per line!)
              HENTRIES=NRESIDUES/10
           ELSE
              HENTRIES=NRESIDUES/10 + 1
           ENDIF
           NDUMMY=1
           DO J3=1,HENTRIES                             !go through all lines
              READ(MYUNIT2,'(A)',END=99) ENTRY_               !read line
              CALL READ_LINE(ENTRY_,NWORDS,ENTRIES)
              J4=1
              DO WHILE(J4.LE.10)
                 IF (NDUMMY.LE.NRESIDUES) THEN
                    READ(ENTRIES(J4),'(I8)') AMBER12_RESSTART(NDUMMY)
                    NDUMMY = NDUMMY + 1
                 ELSE
                    GOTO 98
                 ENDIF
                 J4=J4+1
              ENDDO
           ENDDO
        ENDIF
     ENDDO reading
99   CLOSE(MYUNIT2)
     DO J4=1,NRESIDUES-1
        AMBER12_RESEND(J4) = AMBER12_RESSTART(J4+1) - 1
        AMBER12_RESNATOM(J4) = AMBER12_RESEND(J4) - AMBER12_RESSTART(J4) + 1
     ENDDO
     AMBER12_RESEND(NRESIDUES) = NATOMS
     AMBER12_RESNATOM(NRESIDUES) = AMBER12_RESEND(NRESIDUES) - AMBER12_RESSTART(NRESIDUES) + 1
     IF (DEBUG) THEN
        WRITE(MUTUNIT,'(A)') 'Residue names, start index, end index and number of atoms'
        WRITE(MUTUNIT,*) AMBER12_RESNAME
        WRITE(MUTUNIT,*) AMBER12_RESSTART
        WRITE(MUTUNIT,*) AMBER12_RESEND
        WRITE(MUTUNIT,*) AMBER12_RESNATOM
    ENDIF
  END SUBROUTINE TOPOLOGY_READER
  
  SUBROUTINE READ_LINE(LINE,NWORDS,WORDSOUT)
      CHARACTER(*), INTENT(IN) :: LINE
      INTEGER, INTENT(IN) :: NWORDS
      CHARACTER(*), DIMENSION(NWORDS), INTENT(OUT) :: WORDSOUT
      INTEGER:: J1,START_IND,END_IND,J2
      CHARACTER(25) :: WORD
      START_IND=0
      END_IND=0
      J1=1
      J2=0
      DO WHILE(J1.LE.LEN(LINE))
          IF ((START_IND.EQ.0).AND.(LINE(J1:J1).NE.' ')) THEN
             START_IND=J1
          ENDIF
          IF (START_IND.GT.0) THEN
             IF (LINE(J1:J1).EQ.' ') END_IND=J1-1
             IF (J1.EQ.LEN(LINE)) END_IND=J1
             IF (END_IND.GT.0) THEN
                J2=J2+1
                WORD=LINE(START_IND:END_IND)
                WORDSOUT(J2)=TRIM(WORD)
                START_IND=0
                END_IND=0
             ENDIF
          ENDIF
          J1=J1+1
      ENDDO
  END SUBROUTINE READ_LINE
 
  SUBROUTINE DISTANCE(A1,A2,DIST)
     DOUBLE PRECISION, INTENT(IN) :: A1(3), A2(3)
     DOUBLE PRECISION, INTENT(OUT) :: DIST

     DIST = SQRT((A1(1)-A2(1))**2+(A1(2)-A2(2))**2+(A1(3)-A2(3))**2)
  END SUBROUTINE DISTANCE

  SUBROUTINE AMBER12DECOMP(NA,X,DECOMPOSED)
     INTEGER, INTENT(IN) :: NA
     DOUBLE PRECISION, INTENT(IN) :: X(3*NA)
     DOUBLE PRECISION :: DUMMYE
     DOUBLE PRECISION :: DUMMYGRAD(3*NA)
     TYPE(POT_ENE_REC_C) ,INTENT(OUT) :: DECOMPOSED

     CALL AMBER12_ENERGY_AND_GRADIENT(NATOMS, COORDS, DUMMYE, DUMMYGRAD, DECOMPOSED)
  END SUBROUTINE AMBER12DECOMP

END MODULE AMBER12_MUTATIONS
