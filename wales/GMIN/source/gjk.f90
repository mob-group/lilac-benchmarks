!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2017 David J. Wales
!   This file is part of GMIN.
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
! An implementation of the Gilbert–Johnson–Keerthi algorithm to test whether
! convex shapes overlap, and the expanding polytope algorithm to find the
! minimum distance move to get rid of overlap. See
! http://www.dyn4j.org/2010/04/gjk-gilbert-johnson-keerthi/#gjk-support
! (particularly the video tutorial) for an explanation of GJK and
! http://www.dyn4j.org/2010/05/epa-expanding-polytope-algorithm/ for EPA.
! Questions to jwrm2

MODULE GJK_MODULE

    ! DIMS: Working in 3 Dimensions
    INTEGER, PARAMETER :: DIMS = 3

    ! NVERTS: number of vertices in each polyhedron
    INTEGER :: NVERTS

    ! Stores a vertex in Minkowski difference space
    TYPE :: SUPPORT_POINT
        ! V: the point in difference space
        ! S1, S2: the points of the individual support functions
        DOUBLE PRECISION :: V(DIMS), S1(DIMS), S2(DIMS)
    END TYPE SUPPORT_POINT

    ! Stores a single polyhedron
    TYPE :: POLYHEDRON
        ! R: Coordinates of the centre of the particle
        ! P: Angle-axis style orientation of the particle
        ! RM: Rotation matrix, derived from P
        ! INVRM: Inverse of RM
        ! RMD: Derivatives of the rotation matrix
        DOUBLE PRECISION :: R(DIMS), P(DIMS), RM(DIMS, DIMS), INVRM(DIMS, DIMS)
        DOUBLE PRECISION :: RMD(DIMS, DIMS, DIMS)
        ! VERTS: Vertex positions
        DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: VERTS
    END TYPE POLYHEDRON

    ! Stores a triangular face for the EPA algorithm
    TYPE :: FACE
        ! A, B, C: vertex positions of the faces
        ! NORMAL: normal to the face (AB x AC)
        ! ORIG_DIST: distance to the origin
        TYPE(SUPPORT_POINT) :: A, B, C
        DOUBLE PRECISION :: NORMAL(DIMS), ORIG_DIST
    END TYPE FACE

    ! Stores an edge for the EPA algorithm
    TYPE :: EDGE
        ! A, B: end points of the edge
        TYPE(SUPPORT_POINT) :: A, B
    END TYPE EDGE

    PRIVATE :: NEAREST_SIMPLEX, CHECK_FACE, SAME_DIRECTION, CROSS_121
    PRIVATE :: FACE, NEW_FACE, CHECK_EDGE, BARYCENTRIC

CONTAINS
!-------------------------------------------------------------------------------
! Performs the GJK intersection test
! SHAPE1: The first shape of the test
! SHAPE2: The second shape of the test
! INIT_AXIS: The initial guess, just needs to be a point on the surface of the
!            Minkowski difference of the shapes
! RESULT: whether or not the shapes intersect
! RESULT_SIMPLEX: if the shapes intersect, the simplex containing the origin
!-------------------------------------------------------------------------------
    SUBROUTINE GJK_INTERSECTION(SHAPE1, SHAPE2, INIT_AXIS, RESULT, &
                                RESULT_SIMPLEX)

        ! MYUNIT: File unit of GMIN_out
        ! SUPPORT_FUNC: Support function for convex polyhedra
        ! POLYHEDRON: Type to store a polyhedron
        USE COMMONS, ONLY: MYUNIT

        IMPLICIT NONE

        TYPE(POLYHEDRON), INTENT(IN) :: SHAPE1, SHAPE2
        DOUBLE PRECISION, INTENT(IN) :: INIT_AXIS(DIMS)
        TYPE(SUPPORT_POINT), INTENT(OUT) :: RESULT_SIMPLEX(DIMS+1)
        LOGICAL, INTENT(OUT) :: RESULT

        ! WORK_SIMPLEX: the working simplex
        ! NEW_POINT: the next point to be considered for addition to the simplex
        ! LEN_SIMPLEX: how many vertices the simplex currently contains
        TYPE(SUPPORT_POINT) :: WORK_SIMPLEX(DIMS+1), NEW_POINT

        ! SEARCH_DIRECTION: the current search direction
        DOUBLE PRECISION :: SEARCH_DIRECTION(DIMS)

        ! LEN_SIMPLEX: how many vertices the simplex currently contains
        INTEGER :: LEN_SIMPLEX

        ! Initialise the output variables
        RESULT = .FALSE.

        ! Initialise the working simplex
        WORK_SIMPLEX(1)%S1(:) = SUPPORT_FUNC(SHAPE1,  INIT_AXIS(:))
        WORK_SIMPLEX(1)%S2(:) = SUPPORT_FUNC(SHAPE2, -INIT_AXIS(:))
        WORK_SIMPLEX(1)%V(:) = WORK_SIMPLEX(1)%S1(:) - WORK_SIMPLEX(1)%S2(:)
        SEARCH_DIRECTION(:) = - WORK_SIMPLEX(1)%V(:)
        LEN_SIMPLEX = 1

        ! Start the GJK refinement loop
        DO
            ! Go as far as we can within the Minkowski difference in the search
            ! direction.
            NEW_POINT%S1(:) = SUPPORT_FUNC(SHAPE1,  SEARCH_DIRECTION(:))
            NEW_POINT%S2(:) = SUPPORT_FUNC(SHAPE2, -SEARCH_DIRECTION(:))
            NEW_POINT%V(:) = NEW_POINT%S1(:) - NEW_POINT%S2(:)

            ! If we haven't gone as far as the origin, the shapes definitely
            ! don't intersect
            IF (DOT_PRODUCT(SEARCH_DIRECTION(:), NEW_POINT%V(:)) < 0.D0) THEN
                RESULT = .FALSE.
                EXIT
            END IF

            ! Add the new point to the simplex
            LEN_SIMPLEX = LEN_SIMPLEX + 1
            WORK_SIMPLEX(LEN_SIMPLEX) = NEW_POINT

            ! Check if the new simplex contains the origin
            ! Update the simplex and the search direction
            CALL NEAREST_SIMPLEX(WORK_SIMPLEX, LEN_SIMPLEX, &
                                 SEARCH_DIRECTION, RESULT)

            ! If the simplex contains the origin, we're done
            IF (RESULT) THEN
                ! The simplex length should be DIMS+1 (ie a tetrahedron in 3D),
                ! or something horrible has happened.
                IF (LEN_SIMPLEX .NE. DIMS + 1) THEN
                    WRITE(MYUNIT, *) 'GJK_INTERSECTION> ERROR, length of ', &
                    ' simpled is ', LEN_SIMPLEX, ' but dimension is ', DIMS
                    CLOSE(MYUNIT)
                    STOP
                END IF

                RESULT_SIMPLEX(:) = WORK_SIMPLEX(:)
                EXIT
            END IF

        END DO ! GJK refinement loop

    END SUBROUTINE GJK_INTERSECTION

!-------------------------------------------------------------------------------
! Support function for the GJK algorithm. Returns the vertex furthest along
! the direction DIR
! POLY: the polyhedron of interest
! DIR: the search direction
!-------------------------------------------------------------------------------
        PURE FUNCTION SUPPORT_FUNC(POLY, DIR)

            IMPLICIT NONE

            DOUBLE PRECISION :: SUPPORT_FUNC(DIMS)
            DOUBLE PRECISION, INTENT(IN) :: DIR(:)
            TYPE(POLYHEDRON), INTENT(IN) :: POLY

            ! BEST_INDEX: index of the best vertex
            INTEGER :: I, BEST_INDEX
            ! DOT: temporary dot product
            ! MAX_DOT: largest dot product between DIR and a vertex
            DOUBLE PRECISION :: DOT, MAX_DOT

            BEST_INDEX = 1
            MAX_DOT = -1.0D10

            ! Loop over the vertices
            DO I = 1, NVERTS
                DOT = DOT_PRODUCT(POLY%VERTS(:,I), DIR(:))
                IF (DOT .GT. MAX_DOT) THEN
                    MAX_DOT = DOT
                    BEST_INDEX = I
                END IF
            END DO ! Loop over vertices

            SUPPORT_FUNC(:) = POLY%VERTS(:,BEST_INDEX)
!            WRITE (*,*) 'SEARCH_DIRECTION = ', DIR(:)
!            WRITE (*,*) 'SUPPORT_FUNC = ', SUPPORT_FUNC(:)

        END FUNCTION SUPPORT_FUNC

!-------------------------------------------------------------------------------
! Finds the nearest simplex to the origin of rank lower than the supplied
! simplex. If the simplex is at the maximum size, check if the origin lies
! within the simplex.
! Provides the new search direction based on the closest simplex.
! TEST_SIMPLEX: the simplex to check
! LEN_SIMPLEX: current length of the simplex
! SEARCH_DIRECTION: the new search direction
! RESULT: whether the origin lies within the simplex
!-------------------------------------------------------------------------------
    SUBROUTINE NEAREST_SIMPLEX(TEST_SIMPLEX, LEN_SIMPLEX, &
                               SEARCH_DIRECTION, RESULT)

        ! MYUNIT: File unit of GMIN_out
        ! VEC_CROSS: vector cross product
        ! VEC_RANDOM: a vector with random direction
        USE COMMONS, ONLY: MYUNIT
        USE VEC3, ONLY: VEC_CROSS, VEC_RANDOM

        IMPLICIT NONE

        TYPE(SUPPORT_POINT), INTENT(INOUT) :: TEST_SIMPLEX(DIMS+1)
        DOUBLE PRECISION, INTENT(OUT) :: SEARCH_DIRECTION(DIMS)
        INTEGER, INTENT(INOUT) :: LEN_SIMPLEX
        LOGICAL, INTENT(OUT) :: RESULT

        ! O: the origin
        ! AB, AC, AD, AO, ABPERP, ACPERP, ADPERP, ABCPERP, ACDPERP, ADBPERP:
        !     storage for plane tests and search direction calculation
        ! TOL: For checking whether things are zero
        DOUBLE PRECISION :: AB(DIMS), AC(DIMS), AD(DIMS), AO(DIMS)
        DOUBLE PRECISION :: ABPERP(DIMS), ACPERP(DIMS), ADPERP(DIMS)
        DOUBLE PRECISION :: ABCPERP(DIMS), ACDPERP(DIMS), ADBPERP(DIMS)
        DOUBLE PRECISION, PARAMETER :: TOL = 1.D-06
        ! OVERABC, OVERACD, OVERADB: Results of plane tests
        LOGICAL :: OVERABC, OVERACD, OVERADB

        ! Initialise output
        SEARCH_DIRECTION(:) = 1.D100
        RESULT = .FALSE.

        ! The point referred to as A is always the mose recent added to the
        ! simplex, B the second most recent, etc.
        ! No point can be the closest to the origin. In the GJK loop we've
        ! checked to make sure the origin has been passed. Since it has, there's
        ! no way the points can be closer than the lines connecting them to
        ! other points.
        ! Likewise, only lines and faces involving A are candidates.

        ! Now assuming that DIMS is 3. May need to rewrite.
        IF (LEN_SIMPLEX .EQ. 2) THEN
            ! Simplest case of a line. The nearest simplex can only be the line.
            RESULT = .FALSE.
            AO(:) = - TEST_SIMPLEX(2)%V(:)
            AB(:) = TEST_SIMPLEX(1)%V(:) - TEST_SIMPLEX(2)%V(:)
            SEARCH_DIRECTION(:) = CROSS_121(AB(:), AO(:))

            ! Check if the origin lies on the line, indicated by a zero vector
            IF (ABS(SEARCH_DIRECTION(1)) .LT. TOL .AND. &
                ABS(SEARCH_DIRECTION(2)) .LT. TOL .AND. &
                ABS(SEARCH_DIRECTION(3)) .LT. TOL) THEN
                ! It's on the line. Search direction doesn't really matter.
                SEARCH_DIRECTION(:) = VEC_RANDOM()
            END IF

        ELSE IF (LEN_SIMPLEX .EQ. 3) THEN
            ! Simplex is a triangle. Nearest simplex still cannot be a point.
            ! It could be either of the lines with the latest point, but not the
            ! third.
            ! Origin could lie above or below the triangle.
            RESULT = .FALSE.

            ! Generate the necessary vectors
            AB(:) = TEST_SIMPLEX(2)%V(:) - TEST_SIMPLEX(3)%V(:)
            AC(:) = TEST_SIMPLEX(1)%V(:) - TEST_SIMPLEX(3)%V(:)
            ABCPERP(:) = VEC_CROSS(AB(:), AC(:))

            ! Check whether the origin is outside the triangle on the AB side
            CALL CHECK_FACE(TEST_SIMPLEX, LEN_SIMPLEX, SEARCH_DIRECTION, AB, &
                            AC, ABCPERP)

            ! Check if the origin lies on the line or face, indicated by a zero
            ! vector
            IF (ABS(SEARCH_DIRECTION(1)) .LT. TOL .AND. &
                ABS(SEARCH_DIRECTION(2)) .LT. TOL .AND. &
                ABS(SEARCH_DIRECTION(3)) .LT. TOL) THEN
                ! It's on the line or face. Search direction doesn't really
                ! matter.
                SEARCH_DIRECTION(:) = VEC_RANDOM()
                ! Except it should be 'above' the triangle
                IF (DOT_PRODUCT(SEARCH_DIRECTION(:), ABCPERP(:)) .LT. 0.D0) &
                THEN
                    SEARCH_DIRECTION(:) = -SEARCH_DIRECTION(:)
                END IF
            END IF

        ELSE IF (LEN_SIMPLEX .EQ. 4) THEN
            ! Simplex is a tetrahedron. It might enclose the origin.
            ! Only lines and faces involving A need to be considered.
            ! See http://vec3.ca/gjk/implementation/ for an explanation

            ! Generate the necessary vectors
            AO(:) = - TEST_SIMPLEX(4)%V(:)
            ! Edge vectors
            AB(:) = TEST_SIMPLEX(3)%V(:) - TEST_SIMPLEX(4)%V(:)
            AC(:) = TEST_SIMPLEX(2)%V(:) - TEST_SIMPLEX(4)%V(:)
            AD(:) = TEST_SIMPLEX(1)%V(:) - TEST_SIMPLEX(4)%V(:)
            ! Perpendiculars to faces (pointing outwards)
            ABCPERP(:) = VEC_CROSS(AB(:), AC(:))
            ACDPERP(:) = VEC_CROSS(AC(:), AD(:))
            ADBPERP(:) = VEC_CROSS(AD(:), AB(:))

            ! Store the results of plane tests, so we only have to do them once
            OVERABC = SAME_DIRECTION(ABCPERP(:), AO(:), MYUNIT)
            OVERACD = SAME_DIRECTION(ACDPERP(:), AO(:), MYUNIT)
            OVERADB = SAME_DIRECTION(ADBPERP(:), AO(:), MYUNIT)

            IF (.NOT. (OVERABC .OR. OVERACD .OR. OVERADB)) THEN
                ! Origin is inside the simplex
                RESULT = .TRUE.

            ! Deal with the two faces case
            ELSE IF (OVERABC .AND. OVERACD .AND. .NOT. OVERADB) THEN
                ! Check the common edge, relative to the first face
                ACPERP(:) = VEC_CROSS(ABCPERP(:), AC(:))
                IF (SAME_DIRECTION(ACPERP(:), AO(:), MYUNIT)) THEN
                    ! Now need to do the full check for the second face.
                    ! Fiddle the logicals and let it fall through.
                    OVERABC = .FALSE.
                ELSE
                    ! Just need to decide between the first face and the other
                    ! edge. However, it's simpler (if slightly less efficient)
                    ! to do the full check for the first face.
                    OVERACD = .FALSE.
                ENDIF
            ELSE IF (OVERACD .AND. OVERADB .AND. .NOT. OVERABC) THEN
                ! Check the common edge, relative to the first face
                ADPERP(:) = VEC_CROSS(ACDPERP(:), AD(:))
                IF (SAME_DIRECTION(ADPERP(:), AO(:), MYUNIT)) THEN
                    ! Now need to do the full check for the second face.
                    ! Fiddle the logicals and let it fall through.
                    OVERACD = .FALSE.
                ELSE
                    ! Just need to decide between the first face and the other
                    ! edge. However, it's simpler (if slightly less efficient)
                    ! to do the full check for the first face.
                    OVERADB = .FALSE.
                ENDIF
            ELSE IF (OVERADB .AND. OVERABC .AND. .NOT. OVERACD) THEN
                ! Check the common edge, relative to the first face
                ABPERP(:) = VEC_CROSS(ADBPERP(:), AB(:))
                IF (SAME_DIRECTION(ABPERP(:), AO(:), MYUNIT)) THEN
                    ! Now need to do the full check for the second face.
                    ! Fiddle the logicals and let it fall through.
                    OVERADB = .FALSE.
                ELSE
                    ! Just need to decide between the first face and the other
                    ! edge. However, it's simpler (if slightly less efficient)
                    ! to do the full check for the first face.
                    OVERABC = .FALSE.
                ENDIF
            ENDIF ! Two faces case

            ! Deal with the only one face cases.
            ! Need to check the bordering lines.
            IF (OVERABC .AND. .NOT. OVERACD .AND. .NOT. OVERADB) THEN
                ! Adjust TEXT_SIMPLEX for the ABC face
                TEST_SIMPLEX(1) = TEST_SIMPLEX(2)
                TEST_SIMPLEX(2) = TEST_SIMPLEX(3)
                TEST_SIMPLEX(3) = TEST_SIMPLEX(4)
                ! Check the face and it's edges
                CALL CHECK_FACE(TEST_SIMPLEX, LEN_SIMPLEX, SEARCH_DIRECTION, &
                                AB, AC, ABCPERP)
            ELSE IF (OVERACD .AND. .NOT. OVERADB .AND. .NOT. OVERABC) THEN
                ! Adjust TEXT_SIMPLEX for the ACD face
                TEST_SIMPLEX(3) = TEST_SIMPLEX(4)
                ! Check the face and it's edges
                CALL CHECK_FACE(TEST_SIMPLEX, LEN_SIMPLEX, SEARCH_DIRECTION, &
                                AC, AD, ACDPERP)
            ELSE IF (OVERADB .AND. .NOT. OVERABC .AND. .NOT. OVERACD) THEN
                ! Adjust TEXT_SIMPLEX for the ADB face
                TEST_SIMPLEX(2) = TEST_SIMPLEX(1)
                TEST_SIMPLEX(1) = TEST_SIMPLEX(3)
                TEST_SIMPLEX(3) = TEST_SIMPLEX(4)
                ! Check the face and it's edges
                CALL CHECK_FACE(TEST_SIMPLEX, LEN_SIMPLEX, SEARCH_DIRECTION, &
                                AD, AB, ADBPERP)
            END IF ! One face case
        ELSE
            ! Simplex is not the right size, something horrible has happened
            WRITE (MYUNIT, *) 'NEAREST_SIMPLEX> ERROR, simplex size is ', &
                LEN_SIMPLEX, ' This should not happen.'
            CLOSE(MYUNIT)
            STOP
        END IF ! Simplex size

    END SUBROUTINE NEAREST_SIMPLEX

!-------------------------------------------------------------------------------
! Assuming we know a given face or it's edges are the nearest simplex, find
! which of the two edges or face it is, and generate the result simplex and
! search direction.
! TEST_SIMPLEX: the simplex we're checking. The third point is 'A'
! LEN_SIMPLEX: length of the simplex, 3 at the start
! SEARCH_DIRECTION: new direction to search next
! AB, AC: edge vectors
! ABCPERP: perpendicular to the face
!-------------------------------------------------------------------------------
    SUBROUTINE CHECK_FACE(TEST_SIMPLEX, LEN_SIMPLEX, SEARCH_DIRECTION, AB, AC, &
                          ABCPERP)
        USE VEC3, ONLY: VEC_CROSS
        USE COMMONS, ONLY: MYUNIT

        IMPLICIT NONE

        TYPE(SUPPORT_POINT), INTENT(INOUT) :: TEST_SIMPLEX(DIMS+1)
        INTEGER, INTENT(INOUT) :: LEN_SIMPLEX
        DOUBLE PRECISION, INTENT(IN) :: AB(DIMS), AC(DIMS), ABCPERP(DIMS)
        DOUBLE PRECISION, INTENT(OUT) :: SEARCH_DIRECTION(DIMS)

        ! ABPERP, ACPERP: vectors perpendicular to the lines in the plane of the
        !                 triangle
        ! AO: edge from A to the origin
        DOUBLE PRECISION :: ABPERP(DIMS), ACPERP(DIMS), AO(DIMS)

        AO(:) = - TEST_SIMPLEX(3)%V(:)
        ! Perpendiculars to the edges, in plane of this face
        ABPERP(:) = VEC_CROSS(AB(:), ABCPERP(:))
        ACPERP(:) = VEC_CROSS(ABCPERP(:), AC(:))

        IF (SAME_DIRECTION(ABPERP(:), AO(:), MYUNIT)) THEN
            ! In the AB line region
            TEST_SIMPLEX(1) = TEST_SIMPLEX(2)
            TEST_SIMPLEX(2) = TEST_SIMPLEX(3)
            LEN_SIMPLEX = 2
            SEARCH_DIRECTION(:) = CROSS_121(AB(:), AO(:))
        ELSE IF (SAME_DIRECTION(ACPERP(:), AO(:), MYUNIT)) THEN
            ! In the AC line region
            TEST_SIMPLEX(2) = TEST_SIMPLEX(3)
            LEN_SIMPLEX = 2
            SEARCH_DIRECTION(:) = CROSS_121(AC(:), AO(:))
        ELSE
            ! In the ABC region
            ! Make sure we get the right side
            LEN_SIMPLEX = 3

            IF (SAME_DIRECTION(ABCPERP(:), AO(:), MYUNIT)) THEN
                SEARCH_DIRECTION(:) = ABCPERP(:)
            ELSE
                SEARCH_DIRECTION(:) = -ABCPERP(:)
                ! Swap the simplex ordering, so the normal is 'outside'. Makes
                ! it easier for tetrahedra if this is guaranteed.
                TEST_SIMPLEX(4) = TEST_SIMPLEX(1)
                TEST_SIMPLEX(1) = TEST_SIMPLEX(2)
                TEST_SIMPLEX(2) = TEST_SIMPLEX(4)
            END IF
        END IF
    END SUBROUTINE CHECK_FACE

!-------------------------------------------------------------------------------
! Performs the Expanding Polytope Algorithm to generate overlap information
! SHAPE1, SHAPE2: the two polyhedra to test
! START_SIMPLEX: a simplex containing the origin in the Minkowski difference
!                space. Get from the output of GJK_INTERSECTION
!                We'll assume it's 3D, so simplex is of length 4
! OVERLAP_DIST: the overlap, or penetration distance. Defined as the smallest
!               distance a polyhedron must move to remove overlap
! OVERLAP_VECT: overlap vector, required for derivatives
! WITNESS1, WITNESS2: the points on the shapes at the ends of OVERLAP_VECT
!-------------------------------------------------------------------------------
    SUBROUTINE EXPANDING_POLYTOPE(SHAPE1, SHAPE2, START_SIMPLEX, OVERLAP_DIST, &
                                  OVERLAP_VECT, WITNESS1, WITNESS2)

        ! MYUNIT: File unit for GMIN_out
        USE COMMONS, ONLY: MYUNIT

        IMPLICIT NONE

        TYPE(POLYHEDRON), INTENT(IN) :: SHAPE1, SHAPE2
        TYPE(SUPPORT_POINT), INTENT(IN) :: START_SIMPLEX(DIMS+1)
        DOUBLE PRECISION, INTENT(OUT) :: OVERLAP_DIST, OVERLAP_VECT(DIMS)
        DOUBLE PRECISION, INTENT(OUT) :: WITNESS1(DIMS), WITNESS2(DIMS)

        ! FACES: Stores all the current faces
        ! TEMP_FACES: Used for reallocating the array if necessary
        ! TEMP_FACE: Used for checking the normal
        TYPE(FACE), ALLOCATABLE :: FACES(:), TEMP_FACES(:)
        TYPE(FACE) :: TEMP_FACE

        ! CULL_EDGES: Stores edges while checking for culling
        TYPE(EDGE), ALLOCATABLE :: CULL_EDGES(:)

        ! NEW_POINT: the new point on the surface of the Minkowski difference
        TYPE(SUPPORT_POINT) :: NEW_POINT

        ! CLOSEST_BARY: Barycentric coordinates of the origin projected onto the
        !               closest triangle
        ! MIN_DISTANCE: Keep track of the closest distance to the origin
        ! TOL: Tolerance used for deciding if we've moved closer to the origin
        DOUBLE PRECISION :: CLOSEST_BARY(DIMS), MIN_DISTANCE
        DOUBLE PRECISION, PARAMETER :: TOL = 1.D-06

        ! CLOSEST_INDEX: Index of the closest face to the origin
        ! LEN_CULL: The number of faces to cull
        ! MAX_FACES: Current maximum number of faces in array
        ! NUM_FACES: Current number of faces in the expanding simplex
        ! NUM_NEW: The number of new faces to add
        ! CULL_FACES: Which faces need to be removed at each addition
        INTEGER :: CLOSEST_INDEX, LEN_CULL, MAX_FACES, NUM_FACES, NUM_NEW
        INTEGER, ALLOCATABLE :: CULL_FACES(:)
        INTEGER :: I, J, K

        ! INCLUDE_EDGES: Which edges to make new faces from
        ! DELETE_FACE: Used for a sanity check
        LOGICAL, ALLOCATABLE :: INCLUDE_EDGES(:)
        LOGICAL :: DELETE_FACE

        ! GETUNIT: for finding an unused file unit, from utils.f
        ! TEST_UNIT, TEST_UNIT_2: file unit for testing purposes
!        INTEGER :: GETUNIT, TEST_UNIT, TEST_UNIT_2

        ! Get a file unit for test_vertices.xyz
!        TEST_UNIT = GETUNIT()
!        OPEN(UNIT=TEST_UNIT, FILE="test_vertices.xyz", STATUS="REPLACE")
!        TEST_UNIT_2 = GETUNIT()
!        OPEN(UNIT=TEST_UNIT_2, FILE="test_output.txt", STATUS="REPLACE")

!        WRITE (TEST_UNIT_2,*) 'Starting EXPANDING_POLYTOPE'

        ! Start by allocating 20 faces. Will be doubled if required
        MAX_FACES = 20
        ALLOCATE(FACES(MAX_FACES))

        ! Add the initial simplex to the array
        NUM_FACES = SIZE(START_SIMPLEX)
        IF (NUM_FACES .NE. 4) THEN
            WRITE (MYUNIT, *) 'EXPANDING_POLYTOPE> ERROR, START_SIMPLEX size ',&
                              'is ', NUM_FACES, ' but should be 4'
            CLOSE(MYUNIT)
            STOP
        END IF

        ! ABC
        FACES(1) = NEW_FACE(START_SIMPLEX(4),START_SIMPLEX(3),START_SIMPLEX(2))
        ! ACD
        FACES(2) = NEW_FACE(START_SIMPLEX(4),START_SIMPLEX(2),START_SIMPLEX(1))
        ! ADB
        FACES(3) = NEW_FACE(START_SIMPLEX(4),START_SIMPLEX(1),START_SIMPLEX(3))
        ! CBD
        FACES(4) = NEW_FACE(START_SIMPLEX(2),START_SIMPLEX(3),START_SIMPLEX(1))

!        WRITE (TEST_UNIT_2,*) 'Entering EPA loop'
        ! Now start the refinement loop
        DO
            ! Write all the vertices to a file, for testing
!            WRITE (TEST_UNIT, *) (NUM_FACES*3)
!            WRITE (TEST_UNIT, *) 'Whatever'
!            DO I = 1, NUM_FACES
!                WRITE (TEST_UNIT, *) 'O ', FACES(I)%A%V(1), FACES(I)%A%V(2), &
!                                           FACES(I)%A%V(3)
!                WRITE (TEST_UNIT, *) 'O ', FACES(I)%B%V(1), FACES(I)%B%V(2), &
!                                           FACES(I)%B%V(3)
!                WRITE (TEST_UNIT, *) 'O ', FACES(I)%C%V(1), FACES(I)%C%V(2), &
!                                           FACES(I)%C%V(3)
!            END DO

            ! First we have to identify the closest face to the origin
            ! Fortunately, we already have all the distances
            MIN_DISTANCE = 1.D10
            CLOSEST_INDEX = -1
            DO I = 1, NUM_FACES
!                WRITE (TEST_UNIT_2,*) '*** I = ', I, ' *** ORIG_DIST = ', FACES(I)%ORIG_DIST
!                WRITE (TEST_UNIT_2,*) 'Normal = ', FACES(I)%NORMAL(:)
!                WRITE (TEST_UNIT_2,*) 'A = ', FACES(I)%A%V(:)
!                WRITE (TEST_UNIT_2,*) 'B = ', FACES(I)%B%V(:)
!                WRITE (TEST_UNIT_2,*) 'C = ', FACES(I)%C%V(:)
                IF (FACES(I)%ORIG_DIST .LT. MIN_DISTANCE) THEN
                    MIN_DISTANCE = FACES(I)%ORIG_DIST
                    CLOSEST_INDEX = I
                END IF
            END DO ! Loop over faces

            ! Check that we've actually found something
            IF (CLOSEST_INDEX .EQ. -1) THEN
                WRITE (MYUNIT, *) 'EXPANDING_POLYTOPE> ERROR, failed to find ',&
                                  'closest distance'
                CLOSE(MYUNIT)
                STOP
            END IF

!            WRITE (TEST_UNIT_2,*) 'Found minimum distance ', MIN_DISTANCE, ' Closest face index is ', CLOSEST_INDEX
            ! Generate the new point on the surface of the Minkowski difference
            NEW_POINT%S1 = SUPPORT_FUNC(SHAPE1,  FACES(CLOSEST_INDEX)%NORMAL(:))
            NEW_POINT%S2 = SUPPORT_FUNC(SHAPE2, -FACES(CLOSEST_INDEX)%NORMAL(:))
            NEW_POINT%V(:) = NEW_POINT%S1(:) - NEW_POINT%S2(:)
!            WRITE (TEST_UNIT_2,*) 'NEW_POINT = ', NEW_POINT%V(:)

!            WRITE (TEST_UNIT_2,*) 'Found NEW_POINT, distance check: ', &
!                         DOT_PRODUCT(NEW_POINT%V(:), &
!                         FACES(CLOSEST_INDEX)%NORMAL(:)) &
!                         - FACES(CLOSEST_INDEX)%ORIG_DIST

            ! Check the distance from the origin to the support point against
            ! the distance from the origin to the face
            IF (DOT_PRODUCT(NEW_POINT%V(:), FACES(CLOSEST_INDEX)%NORMAL(:)) &
                - FACES(CLOSEST_INDEX)%ORIG_DIST .LT. TOL) THEN
                ! We've found the closest point
                OVERLAP_DIST = FACES(CLOSEST_INDEX)%ORIG_DIST
                OVERLAP_VECT = FACES(CLOSEST_INDEX)%NORMAL(:)

                ! Find the barycentric coordinates of the origin projected onto
                ! the triangle
                CLOSEST_BARY(:) = BARYCENTRIC(FACES(CLOSEST_INDEX))
!                WRITE(*,*) 'CLOSEST_BARY = ', CLOSEST_BARY(:)
!                WRITE(*,*) 'A%S1 = ', FACES(CLOSEST_INDEX)%A%S1(:)
!                WRITE(*,*) 'A%S1 = ', FACES(CLOSEST_INDEX)%B%S1(:)
!                WRITE(*,*) 'A%S1 = ', FACES(CLOSEST_INDEX)%C%S1(:)

                ! Calculate the witness points
                WITNESS1(:) = CLOSEST_BARY(1)*FACES(CLOSEST_INDEX)%A%S1(:) + &
                              CLOSEST_BARY(2)*FACES(CLOSEST_INDEX)%B%S1(:) + &
                              CLOSEST_BARY(3)*FACES(CLOSEST_INDEX)%C%S1(:)

                WITNESS2(:) = CLOSEST_BARY(1)*FACES(CLOSEST_INDEX)%A%S2(:) + &
                              CLOSEST_BARY(2)*FACES(CLOSEST_INDEX)%B%S2(:) + &
                              CLOSEST_BARY(3)*FACES(CLOSEST_INDEX)%C%S2(:)
                EXIT ! Break out of the main EPA loop
            END IF

            ! We haven't found the answer

            ! Something horrible can happen. It's possible to generate a support
            ! point already in the list of vertices. If everything was working
            ! well, then that means termination. However, due to numerical
            ! inaccuracies, and the possibility of faces being coplanar,
            ! spurious faces can be added to the list. If we have a repeated
            ! vertex as the support point, we want to remove the face that
            ! generated it and try again.
            DELETE_FACE = .FALSE.
            DO I = 1, NUM_FACES
                IF (ALL(ABS(NEW_POINT%V(:) - FACES(I)%A%V(:)) .LT. TOL) .OR. &
                    ALL(ABS(NEW_POINT%V(:) - FACES(I)%B%V(:)) .LT. TOL) .OR. &
                    ALL(ABS(NEW_POINT%V(:) - FACES(I)%C%V(:)) .LT. TOL)) THEN
                    DELETE_FACE = .TRUE.
                    EXIT
                END IF
            END DO

            IF (DELETE_FACE) THEN
                ! We have to delete the face
                ALLOCATE(TEMP_FACES(NUM_FACES - 1))
                ! Use J to keep track of the new index
                J = 1
                DO I = 1, NUM_FACES
                    IF (I .NE. CLOSEST_INDEX) THEN
                        TEMP_FACES(J) = FACES(I)
                        J = J + 1
                    END IF
                END DO
                NUM_FACES = NUM_FACES - 1
                FACES(1:NUM_FACES) = TEMP_FACES(:)
                DEALLOCATE(TEMP_FACES)
                CYCLE ! Cycle main EPA loop
            END IF ! End deleting a face

            ! Do a winding check to see which faces need to be removed
            ! We remove all faces that can be 'seen' from the new point
            ! It's possible that due to numerical rounding errors, this check
            ! will be incorrect.
            ALLOCATE(CULL_FACES(NUM_FACES))
            LEN_CULL = 0
            DO I = 1, NUM_FACES
!                WRITE (TEST_UNIT_2,*) "I = ", I, " new point distance = ", &
!                                      DOT_PRODUCT(FACES(I)%NORMAL(:), &
!                                      NEW_POINT%V(:) - FACES(I)%A%V(:))
                IF (DOT_PRODUCT(FACES(I)%NORMAL(:), &
                    NEW_POINT%V(:) - FACES(I)%A%V(:)) .GT. 0.D0) THEN
                    LEN_CULL = LEN_CULL + 1
                    CULL_FACES(LEN_CULL) = I
                END IF
            END DO

            ! Testing only
!            WRITE(TEST_UNIT_2,*) 'LEN_CULL = ', LEN_CULL
!            DO I = 1, LEN_CULL
!                WRITE (TEST_UNIT_2,*) 'Cull face ', CULL_FACES(I)
!            END DO

            ! Now we loop through the culled edges. We only want to include
            ! unique edges, which gives us the boundary of the 'hole' left by
            ! the culled faces.
            ALLOCATE(CULL_EDGES(LEN_CULL*3))
            ALLOCATE(INCLUDE_EDGES(LEN_CULL*3))
            INCLUDE_EDGES(:) = .TRUE.
            DO I = 1, LEN_CULL ! Loop over the culled faces
                ! Add the new edges to the list
                CULL_EDGES(3*I-2)%A = FACES(CULL_FACES(I))%A
                CULL_EDGES(3*I-2)%B = FACES(CULL_FACES(I))%B
                CULL_EDGES(3*I-1)%A = FACES(CULL_FACES(I))%B
                CULL_EDGES(3*I-1)%B = FACES(CULL_FACES(I))%C
                CULL_EDGES(3*I  )%A = FACES(CULL_FACES(I))%C
                CULL_EDGES(3*I  )%B = FACES(CULL_FACES(I))%A

                DO J = 1, 3 ! Loop over the new edges
                    DO K = 1, 3*(I-1) ! Loop over the old edges
                        IF (INCLUDE_EDGES(K)) THEN
                            IF (CHECK_EDGE(CULL_EDGES(K), &
                                           CULL_EDGES(3*(I-1) + J))) THEN
                                INCLUDE_EDGES(K) = .FALSE.
                                INCLUDE_EDGES(3*(I-1) + J) = .FALSE.
                                EXIT ! Don't need to check any more old edges
                            END IF
                        END IF
                    END DO ! Loop over the old edges
                END DO ! Loop over the new edges
            END DO ! Loop over the culled faces

!            DO I = 1, LEN_CULL*3
!                WRITE (*,*) 'Added edge: ', INCLUDE_EDGES(I), &
!                CULL_EDGES(I)%A%V(1), CULL_EDGES(I)%A%V(2), &
!                CULL_EDGES(I)%A%V(3), ' to ', &
!                CULL_EDGES(I)%B%V(1), CULL_EDGES(I)%B%V(2), &
!                CULL_EDGES(I)%B%V(3)
!            END DO

            ! Get the number of new faces
            NUM_NEW = 0
            DO I = 1, LEN_CULL*3 ! Loop over new edges
                IF (INCLUDE_EDGES(I)) NUM_NEW = NUM_NEW + 1
            END DO

            ! Check whether we're about to exceed the FACES array bounds
            IF (NUM_FACES - LEN_CULL + NUM_NEW .GT. MAX_FACES) THEN
                ! Reallocate the array
                ALLOCATE(TEMP_FACES(NUM_FACES))
                TEMP_FACES(:) = FACES(1:NUM_FACES)
                DEALLOCATE(FACES)
                MAX_FACES = MAX_FACES*2
                ALLOCATE(FACES(MAX_FACES))
                FACES(1:NUM_FACES) = TEMP_FACES(:)
                DEALLOCATE(TEMP_FACES)
            END IF

            ! Testing
            IF (NUM_FACES .GT. 100) THEN
!                TEST_UNIT = GETUNIT()
!                OPEN(UNIT=TEST_UNIT, FILE="test_vertices.xyz", STATUS="REPLACE")
!                WRITE (TEST_UNIT, *) (NUM_FACES*3)
!                WRITE (TEST_UNIT, *) 'Whatever'
!                DO I = 1, NUM_FACES
!                    WRITE(TEST_UNIT, *) 'O ', FACES(I)%A%V(:)
!                    WRITE(TEST_UNIT, *) 'O ', FACES(I)%B%V(:)
!                    WRITE(TEST_UNIT, *) 'O ', FACES(I)%C%V(:)
!                END DO
!                CLOSE(TEST_UNIT)
                WRITE (*,*) 'EPA: High number of faces ', NUM_FACES
                STOP
            END IF

!            WRITE(*,*) 'INCLUDE_EDGES = ', INCLUDE_EDGES(:)
!            DO I = 1, SIZE(INCLUDE_EDGES)
!                WRITE(*, *) 'Edge I A = ', CULL_EDGES(I)%A%V(:), ' B = ', CULL_EDGES(I)%B%V(:)
!            END DO

            ! Add all the new faces
            ! First delete the old ones and compress the array
            ! Use J to keep track of the new index and K to keep track of the
            ! cull index
            J = 1
            K = 1
            ALLOCATE(TEMP_FACES(NUM_FACES - LEN_CULL))
            DO I = 1, NUM_FACES
                IF (CULL_FACES(K) .EQ. I .AND. K .LE. LEN_CULL) THEN
                    K = K + 1
                ELSE
                    TEMP_FACES(J) = FACES(I)
                    J = J + 1
                END IF
            END DO
            NUM_FACES = NUM_FACES - LEN_CULL
            FACES(1:NUM_FACES) = TEMP_FACES(:)
            DEALLOCATE(TEMP_FACES)

            ! Now add the new faces to the end
            DO I = 1, LEN_CULL*3
                IF (INCLUDE_EDGES(I)) THEN
                    TEMP_FACE = NEW_FACE(CULL_EDGES(I)%A, CULL_EDGES(I)%B, &
                                         NEW_POINT)
                    ! It's possible numerical stability errors will have got an
                    ! incorrect face in the mix. If the normal is pointing the
                    ! wrong way, we know this has happened.
                    IF (DOT_PRODUCT(TEMP_FACE%A%V(:), TEMP_FACE%NORMAL(:)) .LT.&
                        0.D0) THEN
                        CYCLE
                    ELSE
                        NUM_FACES = NUM_FACES + 1
                        FACES(NUM_FACES) = TEMP_FACE
                    END IF
                END IF
            END DO

            ! Deallocate the arrays for next time round
            DEALLOCATE(CULL_FACES)
            DEALLOCATE(CULL_EDGES)
            DEALLOCATE(INCLUDE_EDGES)

!            WRITE (TEST_UNIT_2, *) 'End of EPA loop: NUM_FACES = ', NUM_FACES
        END DO ! EPA refinement loop

        ! Deallocate the arrays for next time round
        DEALLOCATE(FACES)

!        WRITE (TEST_UNIT_2, *) 'Ending EPA> NUM_FACES = ', NUM_FACES

!        CLOSE(TEST_UNIT)
!        CLOSE(TEST_UNIT_2)

    END SUBROUTINE EXPANDING_POLYTOPE

!-------------------------------------------------------------------------------
! Utility function for dot products.
! All we care about is whether or not the dot is positive or negative
! Print a warning if it's zero
! VEC1, VEC2: the vectors to dot
! OUT_UNIT: File unit for warnings
!-------------------------------------------------------------------------------
    LOGICAL FUNCTION SAME_DIRECTION(VEC1, VEC2, OUT_UNIT)

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: OUT_UNIT
        DOUBLE PRECISION, INTENT(IN) :: VEC1(DIMS), VEC2(DIMS)

        ! DOT: store the dot product
        DOUBLE PRECISION :: DOT

        DOT = DOT_PRODUCT(VEC1(:), VEC2(:))

        IF (DOT .EQ. 0.D0) THEN
            WRITE(OUT_UNIT, *) 'NEAREST_SIMPLEX> WARNING: exact overlap'
            SAME_DIRECTION = .FALSE.
        ELSE
            SAME_DIRECTION = (DOT .GT. 0.D0)
        END IF
    END FUNCTION SAME_DIRECTION

!-------------------------------------------------------------------------------
! Utility function for cross products.
! Carries out VEC1 x VEC2 x VEC1
! VEC1, VEC2: the vectors to cross
!-------------------------------------------------------------------------------
    PURE FUNCTION CROSS_121(VEC1, VEC2)

        ! VEC_CROSS: cross producr
        USE VEC3, ONLY: VEC_CROSS

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(IN) :: VEC1(DIMS), VEC2(DIMS)
        DOUBLE PRECISION :: CROSS_121(DIMS)

        CROSS_121(:) = VEC_CROSS(VEC_CROSS(VEC1(:), VEC2(:)), VEC1(:))
    END FUNCTION CROSS_121

!-------------------------------------------------------------------------------
! Makes a new FACE. Always use this instead of generating by hand, to make sure
! the normal and origin distance are always calculated.
! A, B, C: Vertices of the face
!-------------------------------------------------------------------------------
    TYPE(FACE) FUNCTION NEW_FACE(A, B, C)

        ! VEC_CROSS: Vector cross product
        ! VEC_LEN: Length of a vector
        ! MYUNIT: File handle for GMIN output
        USE VEC3, ONLY: VEC_CROSS, VEC_LEN
        USE COMMONS, ONLY: MYUNIT

        IMPLICIT NONE

        TYPE(SUPPORT_POINT), INTENT(IN) :: A, B, C

        ! AB, AC: edge vectors
        DOUBLE PRECISION :: AB(DIMS), AC(DIMS)

        ! Set the vertices
        NEW_FACE%A = A
        NEW_FACE%B = B
        NEW_FACE%C = C

        AB(:) = B%V(:) - A%V(:)
        AC(:) = C%V(:) - A%V(:)

        ! Calculate the normal (should be pointing away from the origin)
        NEW_FACE%NORMAL(:) = VEC_CROSS(AB(:), AC(:))

        ! Normalise
        NEW_FACE%NORMAL(:) = NEW_FACE%NORMAL(:) / VEC_LEN(NEW_FACE%NORMAL(:))

        ! Calculate the distance to the origin
        NEW_FACE%ORIG_DIST = ABS(DOT_PRODUCT(NEW_FACE%NORMAL(:), &
                                             NEW_FACE%A%V(:)))

        IF (ALL(AB .EQ. 0.D0) .OR. ALL(AC .EQ. 0.D0)) THEN
            WRITE (MYUNIT,*) 'NEW_FACE> NaN'
            WRITE (MYUNIT,*) 'A = ', A%V(:)
            WRITE (MYUNIT,*) 'B = ', B%V(:)
            WRITE (MYUNIT,*) 'C = ', C%V(:)
            WRITE (MYUNIT,*) 'AB = ', AB(:)
            WRITE (MYUNIT,*) 'AC = ', AB(:)
            WRITE (MYUNIT,*) 'NORMAL = ', NEW_FACE%NORMAL(:)
            STOP
        END IF

    END FUNCTION NEW_FACE

!-------------------------------------------------------------------------------
! Checks whether two edges are the same as each other (regardless of direction)
! EDGE1, EDGE2: The two edges to check
!-------------------------------------------------------------------------------
    PURE LOGICAL FUNCTION CHECK_EDGE(EDGE1, EDGE2)

        IMPLICIT NONE

        TYPE(EDGE), INTENT(IN) :: EDGE1, EDGE2

        ! TOL: numerical tolerance for testing equality
        DOUBLE PRECISION, PARAMETER :: TOL = 1.D-06

        ! DIFFAA, DIFFBB, DIFFAB, DIFFBA: differences between the vertices
        DOUBLE PRECISION :: DIFFAA(DIMS), DIFFBB(DIMS)
        DOUBLE PRECISION :: DIFFAB(DIMS), DIFFBA(DIMS)

        DIFFAA(:) = ABS(EDGE1%A%V(:) - EDGE2%A%V(:))
        DIFFBB(:) = ABS(EDGE1%B%V(:) - EDGE2%B%V(:))
        DIFFAB(:) = ABS(EDGE1%A%V(:) - EDGE2%B%V(:))
        DIFFBA(:) = ABS(EDGE1%B%V(:) - EDGE2%A%V(:))

!        WRITE(*, *) 'Checking edge: ',EDGE1%A%V(1),EDGE1%A%V(2),EDGE1%A%V(3), &
!                               ' to ',EDGE1%B%V(1),EDGE1%B%V(2),EDGE1%B%V(3)
!        WRITE(*, *) 'Against: ', EDGE2%A%V(1), EDGE2%A%V(2), EDGE2%A%V(3), &
!                         ' to ', EDGE2%B%V(1), EDGE2%B%V(2), EDGE2%B%V(3)


        CHECK_EDGE = ((ALL(DIFFAA .LT. TOL) .AND. ALL(DIFFBB .LT. TOL)) .OR. &
                      (ALL(DIFFAB .LT. TOL) .AND. ALL(DIFFBA .LT. TOL)))

!        WRITE(*, *) 'Result is: ', CHECK_EDGE

    END FUNCTION CHECK_EDGE

!-------------------------------------------------------------------------------
! Calculates the projection of the origin onto a triangle in barycentric
! coordinates
! See the bottom of http://hacktank.net/blog/?p=119
! TRIANG: the face we are calculating for
!-------------------------------------------------------------------------------
    PURE FUNCTION BARYCENTRIC(TRIANG)

        IMPLICIT NONE

        TYPE(FACE), INTENT(IN) :: TRIANG
        DOUBLE PRECISION :: BARYCENTRIC(DIMS)

        ! V0, V1, V2, D00, D01, D11, D20, D21, DENOM: temporary quantities
        DOUBLE PRECISION :: V0(DIMS), V1(DIMS), V2(DIMS)
        DOUBLE PRECISION :: D00, D01, D11, D20, D21, DENOM

        V0(:) = TRIANG%B%V(:) - TRIANG%A%V(:)
        V1(:) = TRIANG%C%V(:) - TRIANG%A%V(:)
        V2(:) = TRIANG%ORIG_DIST*TRIANG%NORMAL(:) - TRIANG%A%V(:)

        D00 = DOT_PRODUCT(V0(:), V0(:))
        D01 = DOT_PRODUCT(V0(:), V1(:))
        D11 = DOT_PRODUCT(V1(:), V1(:))
        D20 = DOT_PRODUCT(V2(:), V0(:))
        D21 = DOT_PRODUCT(V2(:), V1(:))

        DENOM = D00*D11 - D01*D01
        BARYCENTRIC(2) = (D11*D20 - D01*D21) / DENOM
        BARYCENTRIC(3) = (D00*D21 - D01*D20) / DENOM
        BARYCENTRIC(1) = 1.D0 - BARYCENTRIC(2) - BARYCENTRIC(3)

    END FUNCTION BARYCENTRIC

END MODULE GJK_MODULE
