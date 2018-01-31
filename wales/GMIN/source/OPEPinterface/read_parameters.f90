!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Reads the parameter file 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE definitions()
    USE COMMONS, ONLY: MYUNIT
    USE defs
    USE RANDOM
    USE md_defs
    USE ion_pair
    USE PORFUNCS
  
    IMPLICIT NONE
    INTEGER :: i, PARAMS_UNIT, GETUNIT
    CHARACTER(8) :: date
    CHARACTER(10) :: time
    CHARACTER(5) :: zone
    INTEGER, DIMENSION(8) :: value
    LOGICAL :: exists_already, end
    CHARACTER WORD*16
    INTEGER :: fchain
    fchain = GETUNIT() 
    use_qbug = .FALSE.             !debug force field
    force_scaling_factor = 1.0D0   !scaling of potential
    ion_pair_control = .FALSE.     !ion pair potential
    ion_pair_scaling = 1.0D0       !ion potential scaling
    PBC = .FALSE.                  !periodic boundary conditions
    BL = 1.0                       !box size
    C_M = .FALSE.                  !centre of mass for pdb output
    idum = 0                       !random number seed
    usextc = .FALSE.              !use xtc and not pdb

    PARAMS_UNIT = GETUNIT()
    CALL FILE_OPEN('OPEP_params', PARAMS_UNIT, .FALSE.)
10  CALL INPUT(END,PARAMS_UNIT)
    IF (.NOT. END) THEN
       CALL READU(WORD)
    ELSE
       GOTO 20
    ENDIF
    IF (WORD.EQ.'    '.OR.WORD.EQ.'NOTE'.OR.WORD.EQ.'COMMENT' &
   &                  .OR.WORD.EQ.'\\'.OR.WORD.EQ."!"         &
   &                  .OR.WORD.EQ."#") THEN 
         GOTO 10
    !force field debug option
    ELSEIF (WORD .EQ. 'use_qbug') THEN
       use_qbug = .TRUE.     
    
    !scaling factor for potential
    ELSEIF (WORD .EQ. 'Potential_Scaling_Factor') THEN
       CALL READF(force_scaling_factor)

    !ion pair potential
    ELSEIF (WORD .EQ. 'Ion_Pair_Potential') THEN
       ion_pair_control = .TRUE.
    ELSEIF (WORD .EQ. 'Ion_Pair_Scaling') THEN
       CALL READF(ion_pair_scaling)

    !periodic boundary condition
    ELSEIF (WORD .EQ. 'Periodic_Boundary_Condition') THEN
       PBC = .TRUE.
       CALL READF(BL)

    !centre of mass for pdb files for periodic boundary conditions
    ELSEIF (WORD .EQ. 'PDB_center_of_mass') THEN
       C_M = .TRUE.

    !random number seed
    ELSEIF (WORD .EQ. 'RANDOM_SEED') THEN
       CALL READI(idum)
       !use the clock to generate a random number
       IF (idum .EQ. 0) THEN                         
          CALL DATE_AND_TIME(date,time,zone,value)
          idum = -1 * MOD( (1000 * value(7) + value(8)),1024)
       ENDIF

    ELSEIF (WORD .EQ. 'usextc') THEN
       usextc = .TRUE.
    ENDIF
    GOTO 10

  ! We now read the number of fragments after all parameters are set
20  INQUIRE(FILE='ichain.dat',EXIST=exists_already)
    IF (exists_already) THEN
       OPEN(UNIT=fchain,FILE='ichain.dat',STATUS='UNKNOWN',ACTION='READ')
       READ(fchain,*) nfrag
       WRITE(MYUNIT,'(A,I6)') 'read ichain for restriction',nfrag     
       ALLOCATE(list_fragments(nfrag,2))   
       DO i=1, nfrag
          READ(fchain,*) list_fragments(i,1),list_fragments(i,2)
       ENDDO
       CLOSE(fchain)
    ELSE
       WRITE(MYUNIT,'(A)') 'ichain.dat does not exist'
       STOP
    ENDIF

    VECSIZE = 3 * NATOMS
    VECSIZE1 =  VECSIZE / NFRAG


!=========================================================================================================

    WRITE(MYUNIT,'(A39,I12)')   ' Number of atoms                     : ', NATOMS
    WRITE(MYUNIT,'(A39,I12)')   ' Number of fragments                 : ', NFRAG
    WRITE(MYUNIT,'(A39,I12)')   ' Random seed                         : ', idum
    WRITE(MYUNIT,'(A39,F12.6)') ' Potential scaling factor            : ', force_scaling_factor

    WRITE(MYUNIT,'(A39,L12)')   ' Periodic Boundary Condition         : ', PBC
    WRITE(MYUNIT,'(A39,F12.6)') ' Box Length                          : ', BL
    WRITE(MYUNIT,'(A39,L12)')   ' center of mass for pdb              : ', C_M


    !CALL read_parameters_md()
    ALLOCATE(pos(VECSIZE))       
    ALLOCATE(posref(VECSIZE))       
    ALLOCATE(force(VECSIZE))       
    ALLOCATE(atomic_type(NATOMS))
    ALLOCATE(mass(vecsize))

  ! We first set-up pointers for the x, y, z components in the position and
  ! forces

    x    => pos(1:3*natoms:3)
    y    => pos(2:3*natoms:3)
    z    => pos(3:3*natoms:3)

    xref => posref(1:3*natoms:3)
    yref => posref(2:3*natoms:3)
    zref => posref(3:3*natoms:3)

    fx   => force(1:3*natoms:3)
    fy   => force(2:3*natoms:3)
    fz   => force(3:3*natoms:3)

   RETURN
END SUBROUTINE definitions

!clean-up routine for end of run
SUBROUTINE end_definitions()
    USE defs
    USE RANDOM
    USE md_defs
    USE ion_pair
    IF (ALLOCATED(pos)) DEALLOCATE(pos)
    IF (ALLOCATED(posref)) DEALLOCATE(posref)
    IF (ALLOCATED(force)) DEALLOCATE(force)
    IF (ALLOCATED(atomic_type)) DEALLOCATE(atomic_type)
    IF (ALLOCATED(mass)) DEALLOCATE(mass)
    RETURN
END SUBROUTINE end_definitions


