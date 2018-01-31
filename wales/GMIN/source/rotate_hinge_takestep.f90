! A takestep type for the plate-folding potential.
MODULE HINGE_MOVES

    USE COMMONS, ONLY: DEBUG, MYUNIT

    IMPLICIT NONE

    INTEGER :: NHINGES                       ! The number of hinges in the system
    INTEGER, ALLOCATABLE :: REFATOMS(:,:)    ! Will hold a set of 4 reference atoms for each hinge
    INTEGER, ALLOCATABLE :: SET_SIZES(:)     ! A list of the number of plates which need to move for each hinge.
    INTEGER, ALLOCATABLE :: PLATE_SETS(:,:)  ! For each hinge, this will hold the list of plates which 
                                             ! need to be moved when the hinge angle is changed. 

CONTAINS

SUBROUTINE HINGE_INITIALISE

    IMPLICIT NONE

    INTEGER :: J1, HINGEUNIT, JUNK
    INTEGER :: GETUNIT

    ! Input file: hingeconfig
    ! Format as follows.
    !
    !
    ! The first line contains the number of (flexible) hinges in the system.
    ! Each hinge is specified by a set of 3 lines. Sets may be separated by a 
    ! blank line, for clarity, but this is not compulsory.
    !
    ! We need to know the following:
    ! - where the hinge is located (we get this by specifying 4 atom indices, 2
    !   of which lie on either side of the hinge)
    ! - which plates need to move (one side of the hinge is the "moving" side;
    !   the plate on that side needs to move, along with all other plates hinged
    !   to it)
    !
    ! The first line of each set contains a number, corresponding to the number
    ! of plates on the "moving" side of the hinge. It will be more efficient to
    ! use the side of the hinge with the smallest number of connected plates as
    ! the "moving side"
    !
    ! The second line of each set contains a list of indices specifying which
    ! plates belong to the moving side. These indices must match the order in
    ! which the plates are specified in rbodyconfig.
    !
    ! The third line contains a list of 4 atom indices. These are the reference
    ! atoms for the hinge. They should fall into 2 pairs (listed one after the
    ! other), with each pair straddling the hinge and joined by a stiff spring.
    ! The two atoms on each side of the hinge should be connected by a vector
    ! parallel to the hinge.

    HINGEUNIT = GETUNIT()
    OPEN(UNIT=HINGEUNIT,FILE='hingeconfig',status='old')
    READ(HINGEUNIT,*) NHINGES
    IF(DEBUG) WRITE(MYUNIT,*) "Reading hingeconfig. NHINGES:", NHINGES
    ! We may eventually want to move to 6 reference atoms so as to compute the
    ! angle between the two hinges as well.
    ALLOCATE(REFATOMS(NHINGES,4))
    ALLOCATE(SET_SIZES(NHINGES))

! jwrm2> First go through the whole file to find the maximum number of plates in
!        a set
    DO J1 = 1, NHINGES
        ! Number of plates this hinge moves
        READ(HINGEUNIT, *) SET_SIZES(J1)
        ! Identity of plates to move, gobble this time sice the array to store
        ! them isn't allocated
        READ(HINGEUNIT, *) JUNK
        ! Identity of reference atoms
        READ(HINGEUNIT, *) REFATOMS(J1, :)
        IF(DEBUG) THEN
            WRITE(MYUNIT,*) "Reference atoms:"
            WRITE(MYUNIT,*) REFATOMS(J1,:)
        ENDIF
    END DO

    IF (DEBUG) THEN
        WRITE (MYUNIT, *) "SET_SIZES = ", SET_SIZES(:)
    END IF

! jwrm2> Now we can allocate enough space for the plate sets
    ALLOCATE(PLATE_SETS(NHINGES, MAXVAL(SET_SIZES)))

! jwrm2> Read the file again to get the plate sets
    REWIND (UNIT=HINGEUNIT)

! jwrm2> Gobble number of hinges
    READ (HINGEUNIT, *) JUNK

    DO J1 = 1, NHINGES
        ! Gobble set size
        READ(HINGEUNIT,*) JUNK

        ! Get plate set
        IF(DEBUG) WRITE(MYUNIT,'(A,I2,A)') "Hinge ", J1, ". Moving side:"
        READ(HINGEUNIT,*) PLATE_SETS(J1,1:SET_SIZES(J1))
        IF(DEBUG) WRITE(MYUNIT,*) PLATE_SETS(J1,1:SET_SIZES(J1))

        ! Gobble reference atoms
        READ(HINGEUNIT,*) JUNK
    END DO

    CLOSE(HINGEUNIT)

END SUBROUTINE HINGE_INITIALISE

SUBROUTINE HINGE_ROTATE(XCOORDS, ROTATEFACTOR)

    USE COMMONS, ONLY: NATOMS
    USE ROTATIONS, ONLY: rot_aa2mx
    USE GENRIGID, ONLY: NSITEPERBODY, RIGIDGROUPS

    IMPLICIT NONE

    INTEGER :: J1, J2, J3, ATOM, PLATE_INDEX
    DOUBLE PRECISION :: REFPOINT(3), PI, TWOPI, DPRAND
    DOUBLE PRECISION :: AXIS(3), MATRIX(3,3), TOROTATE(3), ATOMROTATED(3)
    DOUBLE PRECISION :: ANGLE, NORM
    DOUBLE PRECISION, INTENT(INOUT) :: XCOORDS(3*NATOMS)
    DOUBLE PRECISION, INTENT(IN) :: ROTATEFACTOR

    ! Current implementation - somewhat crude - is to change every hinge every
    ! time the subroutine is called. We pick a random rotation angle between
    ! -pi*rotatefactor and pi*rotatefactor

!    write(*,*) "Before hinge takestep"
!    write(*,*) XCOORDS

    MATRIX(:,:) = 0.0D0
    TOROTATE(:) = 0.0D0
    ! Define some constants
    PI=ATAN(1.0D0)*4
    TWOPI=2.0D0*PI

    ! Loop over all hinges
    DO J1 = 1, NHINGES

        ! Calculate the centrepoint of the hinge
        REFPOINT(:) = 0.0D0
        DO J2 = 1,4
            ATOM = REFATOMS(J1,J2)
            REFPOINT(:) = REFPOINT(:) + XCOORDS(3*ATOM-2:3*ATOM)
        ENDDO
        REFPOINT(:) = REFPOINT(:)/4.0D0

        ! Calculate the hinge axis (from the centrepoint of reference atoms 1 and 2 to the centrepoint of 3 and 4)
        AXIS(:) = XCOORDS(3*REFATOMS(J1,3)-2:3*REFATOMS(J1,3)) + XCOORDS(3*REFATOMS(J1,4)-2:3*REFATOMS(J1,4)) - &
                  XCOORDS(3*REFATOMS(J1,1)-2:3*REFATOMS(J1,1)) - XCOORDS(3*REFATOMS(J1,2)-2:3*REFATOMS(J1,2))
!       AXIS(:) = AXIS(:) * 0.5D0   ! This line is needed for the definition of the centrepoints, but since
                                    ! we are about to normalise anyway, we don't really need it.
        ! Normalise the axis
        NORM = SQRT(DOT_PRODUCT(AXIS,AXIS))
        AXIS(:) = AXIS(:)/NORM

        ! Centre the system about the hinge centrepoint
        DO J2 = 1,NATOMS
            XCOORDS(3*J2-2:3*J2) = XCOORDS(3*J2-2:3*J2) - REFPOINT(:)
        ENDDO

!        write(*,*) "After first translation"
!        write(*,*) XCOORDS

        ! Choose a random angle by which to rotate
        ANGLE = (DPRAND()-0.5)*TWOPI*ROTATEFACTOR

        ! Create an angle-axis vector corresponding to the desired rotation, and convert it to a rotation matrix
        AXIS(:) = AXIS(:)*ANGLE
        !write(*,*) "Angle Axis:", AXIS(:)
        MATRIX = rot_aa2mx(AXIS)
        !write(*,*) "Matrix"
        !DO J2=1,3
        !    write(*,*) MATRIX(J2,:)
        !ENDDO

        ! Apply the rotation matrix to the atoms in the rigid bodies indicated by PLATE_SETS(J1)
        DO J2 = 1, SET_SIZES(J1)
            PLATE_INDEX = PLATE_SETS(J1,J2)
            DO J3 = 1, NSITEPERBODY(PLATE_INDEX)
                ATOM = RIGIDGROUPS(J3, PLATE_INDEX)
                !write(*,'(A,I3,A,I1,A,I1)') "Rotating atom ", ATOM, " in RB ", PLATE_INDEX, " for hinge ", J1
                TOROTATE = XCOORDS(3*ATOM-2:3*ATOM)
                !write(*,'(A,3F20.10)') "Before rotation: ", TOROTATE(1), TOROTATE(2), TOROTATE(3)
                ATOMROTATED=MATMUL(MATRIX,TOROTATE)
                !write(*,'(A,3F20.10)') "After rotation: ", ATOMROTATED(1), ATOMROTATED(2), ATOMROTATED(3)
                XCOORDS(3*ATOM-2:3*ATOM) = ATOMROTATED
            ENDDO
        ENDDO

!        write(*,*) "After rotation, before back-translation"
!        write(*,*) XCOORDS

        ! Translate the system back to the old centre of mass
        DO J2 = 1, NATOMS
            XCOORDS(3*J2-2:3*J2) = XCOORDS(3*J2-2:3*J2) + REFPOINT(:)
        ENDDO

!        write(*,*) "After hinge", J1
!        write(*,*) XCOORDS


    ENDDO  ! End of loop over hinges

!    write(*,*) "After hinge takestep"
!    write(*,*) XCOORDS
END SUBROUTINE HINGE_ROTATE

END MODULE

