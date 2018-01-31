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
! Potential for convex polyhedra. The polyhedra only interact when they are. In
! this case, there is a harmonic repulsion based on the smallest distance
! required to make them not overlap. Additionally, there is a harmonic
! compression towards the origin on all particles, to promote dense packing.
! Questions to jwrm2.

MODULE CONVEX_POLYHEDRA_MODULE

! POLYHEDRON: type to store a polyhedron
    USE GJK_MODULE, ONLY: POLYHEDRON, DIMS

    IMPLICIT none

    ! NPOLY: Number of polyhedra
    ! X_SIZE: Size of the coordinates array
    INTEGER :: NPOLY, X_SIZE

    ! The array of polyhedra
    TYPE(POLYHEDRON), ALLOCATABLE :: POLYHEDRA(:)

    ! K_COMPRESS: Potential parameter, compression towards the centre
    ! K_OVERLAP: Potential parameter, cost of overlap
    ! MAX_VERT: distance squared of the furthest vertex from the origin
    DOUBLE PRECISION :: K_COMPRESS, K_OVERLAP, MAX_VERT_SQ

    ! Array of the reference vertex positions
    DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:,:) :: VERTS_0

    PRIVATE :: VERTS_0, NPOLY, X_SIZE, POLYHEDRA, UPDATE_POLYHEDRON

CONTAINS
!-------------------------------------------------------------------------------
! Initialise the system, including reading the vertex list from the file
! "polyhedron_vertices.dat"
!-------------------------------------------------------------------------------
    SUBROUTINE INITIALISE_POLYHEDRA()

        ! MYUNIT: File unit of GMIN_out
        ! NATOMS: the number of coordinate lines in the coords file, so twice
        !         the actual number of bodies
        ! DIMS: Number of dimensions of the system
        ! NVERTS: the number of vertices per polyhedron
        ! VEC_CROSS: finds the cross product of two vectors
        USE COMMONS, ONLY: MYUNIT, NATOMS
        USE GJK_MODULE, ONLY: DIMS, NVERTS
        USE VEC3, ONLY: VEC_CROSS

        IMPLICIT NONE

        ! GETUNIT: Function for finding an unused file unit, from utils.f
        ! VERT_UNIT: File unit for "polyhedron_vertices.dat"
        INTEGER :: GETUNIT, VERT_UNIT, I

        ! NORMAL: Normal to the first face read
        ! TEMP_DIST: Distance squared of a vertex from the origin
        ! TOLERANCE: Numerical tolerance for coplanar test
        DOUBLE PRECISION :: NORMAL(DIMS), TEMP_DIST
        DOUBLE PRECISION, PARAMETER :: TOL = 1.0D-6

        ! COPLANAR: whether all the vertices are coplanar
        LOGICAL :: COPLANAR

        ! Set global module parameters
        NPOLY = NATOMS / 2
        X_SIZE = NATOMS * 3

        WRITE (MYUNIT, *) 'INITIALISE_POLYHEDRA> ', NPOLY, ' polyhedra'

        ! Allocate array of polyhedra
        ! We don't need to worry about the positions and orientations here,
        ! they'll be dealt with the first time the potential is called.
        ALLOCATE(POLYHEDRA(NPOLY))

        ! Read the reference vertex positions from file.
        ! Format is
        !   NVERT
        !   VERT(X, 1) VERT(Y, 1) VERT(Z, 1)
        !   VERT(X, 2) VERT(Y, 2) VERT(Z, 2)
        !   ...
        VERT_UNIT = GETUNIT()
        OPEN(UNIT=VERT_UNIT, FILE="polyhedron_vertices.dat", STATUS="old")
        READ(VERT_UNIT, *) NVERTS

        ! Check that we've got at least four vertices
        IF (NVERTS .LT. 4) THEN
            WRITE (MYUNIT, *) &
                'INITIALISE_POLYHEDA> ERROR: require at least 4 vertices'
            CLOSE(MYUNIT)
            STOP
        END IF

        ! Read the reference vertices
        ALLOCATE(VERTS_0(DIMS, NVERTS))
        DO I = 1, NVERTS
            READ(VERT_UNIT, *) VERTS_0(1,I), VERTS_0(2,I), VERTS_0(3,I)
        END DO

        ! Check that the vertices are not all coplanar
        NORMAL(:) = VEC_CROSS((VERTS_0(:, 1) - VERTS_0(:, 2)), &
                               (VERTS_0(:, 1) - VERTS_0(:, 3)))
        COPLANAR = .TRUE.
        DO I = 4, NVERTS ! Loop over vertices
            IF (ABS(DOT_PRODUCT((VERTS_0(:, I) - VERTS_0(:, 1)), NORMAL)) &
                .GT. TOL) THEN
                COPLANAR = .FALSE.
                EXIT
            END IF
        END DO ! Loop over vertices

        ! Error out if the vertices are all coplanar
        IF (COPLANAR) THEN
            WRITE (MYUNIT, *) &
                'INITIALISE_POLYHEDRA> ERROR: all vertices are coplanar'
            CLOSE(MYUNIT)
            STOP
        END IF

        ! Work out the maximum vertex distance
        MAX_VERT_SQ = 0.D0
        DO I = 1, NVERTS
            TEMP_DIST = DOT_PRODUCT(VERTS_0(:, I), VERTS_0(:, I))
            MAX_VERT_SQ = MAX(MAX_VERT_SQ, TEMP_DIST)
        END DO

        ! Copy the reference vertices to the polyhedra
        DO I = 1, NPOLY
            ALLOCATE(POLYHEDRA(I)%VERTS(DIMS, NVERTS))
            POLYHEDRA(I)%VERTS(:, :) = VERTS_0(:, :)
        END DO

    END SUBROUTINE INITIALISE_POLYHEDRA

!-------------------------------------------------------------------------------
! Updates a polyhedron with the rotation matrix and derivatives,
! and the position of the vertices.
! POLY: the polyhedron to operate on
! R: Coordinates of the centre of the particle
! P: Angle-axis style orientation of the particle
!    INTENT(INOUT) because orientations are reranged
! GTEST: Whether to computer derivatives
!-------------------------------------------------------------------------------
    SUBROUTINE UPDATE_POLYHEDRON(POLY, R, P, GTEST)

        ! DIMS: Dimension of the space
        ! NVERTS: Number of vertices per polyhedron
        ! POLYHEDRON: Type for storing a polyhedron
        ! INVERT3X3: Matrix inversion
        USE GJK_MODULE, ONLY: DIMS, NVERTS, POLYHEDRON
        USE VEC3, ONLY: INVERT3X3

        IMPLICIT NONE

        TYPE(POLYHEDRON), INTENT(OUT) :: POLY
        DOUBLE PRECISION, INTENT(IN) :: R(DIMS)
        DOUBLE PRECISION, INTENT(INOUT) :: P(DIMS)
        LOGICAL, INTENT(IN) :: GTEST

        INTEGER :: I
        DOUBLE PRECISION :: PI = 4.0D0 * ATAN(1.0D0)
        ! PMOD: angle of the angle axis orientation
        DOUBLE PRECISION :: PMOD

        POLY%R(:) = R(:)

        ! Make sure that 0 < |p| < 2*pi
        PMOD = sqrt(dot_product(p, p))
        IF(PMOD > 2 * PI) P(:) = P(:) / (PMOD * MOD(PMOD, 2 * PI))
        POLY%P(:) = P(:)

        ! Compute the rotation matrices and derivatives
        CALL RMDRVT(POLY%P, POLY%RM, POLY%RMD(:,:,1), POLY%RMD(:,:,2), &
            POLY%RMD(:,:,3), GTEST)

        ! We need the inverse of the rotation matrix for derivatives
        IF (GTEST) CALL INVERT3X3(POLY%RM(:,:), POLY%INVRM(:,:))

        ! Apply rotation matrix and translational offset to all the vertices
        ALLOCATE(POLY%VERTS(DIMS,NVERTS))
        DO I = 1, NVERTS
            POLY%VERTS(:, I) = POLY%R(:) + MATMUL(POLY%RM, VERTS_0(:, I))
        END DO

    END SUBROUTINE UPDATE_POLYHEDRON

!-------------------------------------------------------------------------------
! Calculates the potential, and optionally the gradients
! X: array of coordinates, all translations, then all orientations
!    INTENT(INOUT) because orientation may be reranged
! GRAD: array of gradients, all translations, then all orientations
! ENERGY: the energy of the system
! GTEST: whether gradients are required
!-------------------------------------------------------------------------------
    SUBROUTINE CONVEX_POLYHEDRA(X, GRAD, ENERGY, GTEST)

        ! GJK_INTERSECTION: Tests whether two shapes overlap
        ! EXPANDING_POLYTOP: Finds the intersection information of two shapes
        ! DIMS: Dimension of the space
        ! SUPPORT_POINT: Type for storing points in Minkowski difference space
        ! NVERTS: testing only, the number of vertices
        USE GJK_MODULE, ONLY: GJK_INTERSECTION, EXPANDING_POLYTOPE, DIMS, &
                              SUPPORT_POINT !, NVERTS

        IMPLICIT NONE

        DOUBLE PRECISION, INTENT(INOUT) :: X(X_SIZE)
        DOUBLE PRECISION, INTENT(OUT) :: GRAD(X_SIZE), ENERGY
        LOGICAL, INTENT(IN) :: GTEST

        ! MOLI_IND, MOLJ_IND: indices for looping over molecules
        ! OFFSET: The end of the position coordinates in X
        INTEGER :: MOLI_IND, MOLJ_IND, OFFSET

        ! RESULT_SIMPLEX: A simplex enclosing the origin if GJK finds an
        !                 intersection
        TYPE(SUPPORT_POINT) :: RESULT_SIMPLEX(DIMS+1)

        ! INIT_AXIS: Initial direction guess for GJK, just use the difference of
        !            the centres
        ! CENTRE_DIST: distance squared from the centre
        ! OVERLAP_VECT: Normalised vector in the overlap direction, from EPA
        ! OVERLAP_DIST: Penetration distance, from EPA
        ! WITNESS1, WITNESS2: The witness points, from EPA
        ! MOLI_REF, MOLJ_REF: the witness points in the reference geometries
        DOUBLE PRECISION :: INIT_AXIS(DIMS)
        DOUBLE PRECISION :: CENTRE_DIST, OVERLAP_VECT(DIMS), OVERLAP_DIST
        DOUBLE PRECISION :: WITNESS1(DIMS), WITNESS2(DIMS)
        DOUBLE PRECISION :: MOLI_REF(DIMS), MOLJ_REF(DIMS)

        ! RESULT: whether or not GJK finds an intersection
        LOGICAL :: RESULT

        ! Testing
!        INTEGER :: GETUNIT, TEST_UNIT, I
!        TEST_UNIT = GETUNIT()
!        OPEN(UNIT=TEST_UNIT, FILE="test_vertices.xyz", STATUS="REPLACE")

        ! Initialise ouput variables
        ENERGY = 0.D0
        IF (GTEST) GRAD(:) = 0.D0

        ! Print out all coordinates, for testing
!        DO MOLI_IND = 1, NPOLY
!            WRITE (*,*) 'Particle ', MOLI_IND, ', trans ', X(MOLI_IND*3-2:MOLI_IND*3)
!            WRITE (*,*) 'Particle ', MOLI_IND, ', rot ', X(MOLI_IND*3-2+OFFSET:MOLI_IND*3+OFFSET)
!        END DO

        ! Calculate offset from the number of particles
        OFFSET = 3*NPOLY

        ! Update all the particle positions and rotation matrices
        DO MOLI_IND = 1, NPOLY
            CALL UPDATE_POLYHEDRON(POLYHEDRA(MOLI_IND), &
                                   X(MOLI_IND*3-2:MOLI_IND*3), &
                                   X(OFFSET+MOLI_IND*3-2:OFFSET+MOLI_IND*3), &
                                   GTEST)
        END DO

        ! Just looking to test at the moment
        ! Testing with two particles, just print whether they overlap or not
!        INIT_AXIS(:) = POLYHEDRA(1)%R(:) - POLYHEDRA(2)%R(:)
!        CALL GJK_INTERSECTION(POLYHEDRA(1), POLYHEDRA(2), INIT_AXIS(:), &
!                              RESULT, RESULT_SIMPLEX(:))

!        WRITE(*,*) 'GJK returned ', RESULT
!        IF (RESULT) THEN
!            WRITE(*,*) 'Result simplex is:'
!            DO MOLJ_IND = 1, 4
!                WRITE(*,*) RESULT_SIMPLEX(MOLJ_IND)%V(1), &
!                           RESULT_SIMPLEX(MOLJ_IND)%V(2), &
!                           RESULT_SIMPLEX(MOLJ_IND)%V(3)
!            END DO
!        END IF

!        IF (RESULT) THEN
!            WRITE(*,*) 'Calling EXPANDING_POLYTOPE'
!            CALL EXPANDING_POLYTOPE(POLYHEDRA(1), POLYHEDRA(2), &
!                                    RESULT_SIMPLEX(:), OVERLAP_DIST, &
!                                    OVERLAP_VECT(:), WITNESS1(:), WITNESS2(:))
!            WRITE (*,*) 'Overlap distance is ', OVERLAP_DIST
!            WRITE (*,*) 'Overlap vector is ', OVERLAP_VECT(:)
!            WRITE (*,*) 'Witness point on shape 1 is ', WITNESS1(:)
!            WRITE (*,*) 'Witness point on shape 2 is ', WITNESS2(:)
!        END IF

!        WRITE (*,*) 'Exiting now'
!        STOP

    ! Harmonic compression towards the origin
    DO MOLI_IND = 1, NPOLY
        CENTRE_DIST = DOT_PRODUCT(POLYHEDRA(MOLI_IND)%R, POLYHEDRA(MOLI_IND)%R)
        ENERGY = ENERGY + 0.5 * K_COMPRESS * CENTRE_DIST

        IF (GTEST) THEN
            GRAD(MOLI_IND*3-2:MOLI_IND*3) = GRAD(MOLI_IND*3-2:MOLI_IND*3) + &
                                           K_COMPRESS*X(MOLI_IND*3-2:MOLI_IND*3)
        END IF
    END DO ! Harmonic compression loop

    ! Loop over all pairs of particles
    DO MOLI_IND = 1, NPOLY-1
        DO MOLJ_IND = MOLI_IND+1, NPOLY

            ! Get the vector between the body centres
            INIT_AXIS(:) = POLYHEDRA(MOLI_IND)%R(:) - POLYHEDRA(MOLJ_IND)%R(:)

            ! If the distance betweent the body centres is greater than twice
            ! the max vertex distance, we can completely skip this pair
            IF (DOT_PRODUCT(INIT_AXIS(:), INIT_AXIS(:)) .GT. 4.D0*MAX_VERT_SQ) &
                CYCLE

            ! Otherwise we have to call GJK to see if the shapes intersect
            CALL GJK_INTERSECTION(POLYHEDRA(MOLI_IND), POLYHEDRA(MOLJ_IND), &
                                  INIT_AXIS(:), RESULT, RESULT_SIMPLEX(:))

            ! If they don't intersect, we can skip this pair
            IF (.NOT. RESULT) CYCLE

            ! If they do intersect, we run EPA to get the penetration distance
            ! and witness points.
            CALL EXPANDING_POLYTOPE(POLYHEDRA(1), POLYHEDRA(2), &
                                    RESULT_SIMPLEX(:), OVERLAP_DIST, &
                                    OVERLAP_VECT(:), WITNESS1(:), WITNESS2(:))

!            WRITE (*,*) 'Overlap distance is ', OVERLAP_DIST
!            WRITE (*,*) 'Overlap vector is ', OVERLAP_VECT(:)
!            WRITE (*,*) 'Witness point on shape 1 is ', WITNESS1(:)
!            WRITE (*,*) 'Witness point on shape 2 is ', WITNESS2(:)

!            WRITE (TEST_UNIT, *) (2*NVERTS + 2)
!            WRITE (TEST_UNIT, *) 'Whatever'
!            DO I = 1, NVERTS
!                WRITE (TEST_UNIT, *) 'C ', POLYHEDRA(MOLI_IND)%VERTS(:, I)
!            END DO
!            DO I = 1, NVERTS
!                WRITE (TEST_UNIT, *) 'N ', POLYHEDRA(MOLJ_IND)%VERTS(:, I)
!            END DO
!            WRITE (TEST_UNIT, *) 'O ', WITNESS1(:)
!            WRITE (TEST_UNIT, *) 'F ', WITNESS2(:)

            ! Unnormalised overlap vector is more useful
            OVERLAP_VECT(:) = OVERLAP_VECT(:) * OVERLAP_DIST

            ! Add the harmonic overlap repulsion
            ENERGY = ENERGY + 0.5D0 * K_OVERLAP * OVERLAP_DIST**2

            IF (GTEST) THEN
                ! Translational derivatives
                GRAD(3*MOLI_IND-2:3*MOLI_IND) = GRAD(3*MOLI_IND-2:3*MOLI_IND) &
                                                + K_OVERLAP * OVERLAP_VECT(:)
                GRAD(3*MOLJ_IND-2:3*MOLJ_IND) = GRAD(3*MOLJ_IND-2:3*MOLJ_IND) &
                                                - K_OVERLAP * OVERLAP_VECT(:)

                ! Rotational derivatives. We need to get the witness points in
                ! the reference geometries
                MOLI_REF(:) = MATMUL(POLYHEDRA(MOLI_IND)%INVRM, WITNESS1(:) - &
                                     POLYHEDRA(MOLI_IND)%R(:))
                MOLJ_REF(:) = MATMUL(POLYHEDRA(MOLJ_IND)%INVRM, WITNESS1(:) - &
                                     POLYHEDRA(MOLJ_IND)%R(:))

!                WRITE (*,*) 'MOLI_REF = ', MOLI_REF(:)
!                WRITE (*,*) 'MOLJ_REF = ', MOLJ_REF(:)

                ! Get the rotational derivatives by taking the dot product of
                ! the overlap vector with the derivative of the rotation matrix
                ! acting on the witness point in the reference geometry
                GRAD(3*MOLI_IND-2+OFFSET) = GRAD(3*MOLI_IND-2+OFFSET) + &
                    K_OVERLAP*DOT_PRODUCT(MATMUL(POLYHEDRA(MOLI_IND)%RMD(:,:,1)&
                    , MOLI_REF(:)), OVERLAP_VECT(:))
                GRAD(3*MOLI_IND-1+OFFSET) = GRAD(3*MOLI_IND-1+OFFSET) + &
                    K_OVERLAP*DOT_PRODUCT(MATMUL(POLYHEDRA(MOLI_IND)%RMD(:,:,2)&
                    , MOLI_REF(:)), OVERLAP_VECT(:))
                GRAD(3*MOLI_IND + OFFSET) = GRAD(3*MOLI_IND + OFFSET) + &
                    K_OVERLAP*DOT_PRODUCT(MATMUL(POLYHEDRA(MOLI_IND)%RMD(:,:,3)&
                    , MOLI_REF(:)), OVERLAP_VECT(:))
                GRAD(3*MOLJ_IND-2+OFFSET) = GRAD(3*MOLJ_IND-2+OFFSET) - &
                    K_OVERLAP*DOT_PRODUCT(MATMUL(POLYHEDRA(MOLJ_IND)%RMD(:,:,1)&
                    , MOLJ_REF(:)), OVERLAP_VECT(:))
                GRAD(3*MOLJ_IND-1+OFFSET) = GRAD(3*MOLJ_IND-1+OFFSET) - &
                    K_OVERLAP*DOT_PRODUCT(MATMUL(POLYHEDRA(MOLJ_IND)%RMD(:,:,2)&
                    , MOLJ_REF(:)), OVERLAP_VECT(:))
                GRAD(3*MOLJ_IND + OFFSET) = GRAD(3*MOLJ_IND + OFFSET) - &
                    K_OVERLAP*DOT_PRODUCT(MATMUL(POLYHEDRA(MOLJ_IND)%RMD(:,:,3)&
                    , MOLJ_REF(:)), OVERLAP_VECT(:))
            END IF

        END DO ! Inner particle loop
    END DO ! Outer particle loop

    ! Testing
!    CLOSE(TEST_UNIT)

    END SUBROUTINE CONVEX_POLYHEDRA

!------------------------------------------------------------------------------
! Prints out the vertex information to the file poly.xyz
!------------------------------------------------------------------------------
    SUBROUTINE VIEW_POLYHEDRA()

        IMPLICIT NONE

        ! Doesn't actually do anything yet

    END SUBROUTINE VIEW_POLYHEDRA

END MODULE CONVEX_POLYHEDRA_MODULE
