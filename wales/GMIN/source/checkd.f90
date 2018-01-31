      SUBROUTINE CHECKD(X)

      USE COMMONS, ONLY: NATOMS, COMPRESST, PERCOLATET, CHECKDID, GTHOMSONT, &
                         BOXDERIVT, ORTHO, BOX_PARAMS, BOX_PARAMSGRAD
      USE GENRIGID, ONLY: RIGIDINIT, ATOMRIGIDCOORDT, DEGFREEDOMS, TRANSFORMCTORIGID

      USE MODHESS
      IMPLICIT NONE

      INTEGER          :: IVRNO, IVRNO1, IVRNO2, J1, J3, dof, doff
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), ENERGY, FM, FP, DFA, DFN, TMPCOORDS(3*NATOMS)
      double precision :: box_paramsold(6)
      LOGICAL          :: GTEST, STEST, COMPON
      DOUBLE PRECISION, PARAMETER :: ERRLIM = 1.D-04, DELX = 1.0D-6
      COMMON /CO/ COMPON

      ! dj337: allow for rigid bodies
      if (rigidinit) then
         dof = degfreedoms
      else
         dof = 3*natoms
      endif

      print *, 'DELX: ', DELX

! jwrm2> Turning compression on, if required
      IF (COMPRESST .OR. PERCOLATET) COMPON = .TRUE.

! jwrm2> Converting GTHOMSON coordinates to polars
      IF (GTHOMSONT) THEN
        CALL GTHOMSONCTOANG(X(1:3*NATOMS), TMPCOORDS(1:3*NATOMS), NATOMS)
        X(1:3*NATOMS) = TMPCOORDS(1:3*NATOMS)
      END IF

      STEST = .FALSE.

      IF (CHECKDID == 0) THEN
         GTEST = .FALSE.
         CALL POTENTIAL (X, G, ENERGY, GTEST, STEST)
         WRITE(*, *) 'Energy  = ', ENERGY

      ELSEIF (CHECKDID == 1) THEN

!     Checks gradients

      ! check derivatives wrt atomic positions
      DO IVRNO = 1, DOF
         WRITE(*, *) IVRNO

         ! dj337: make sure coordinates rigid body
         if (rigidinit.and.atomrigidcoordt) then
            call transformctorigid(x, tmpcoords)
            x(1:degfreedoms) = tmpcoords(1:degfreedoms)
            x(degfreedoms+1:3*natoms) = 0.0d0
            atomrigidcoordt = .false.
         endif

         !print *, 'x: ', x(ivrno)

         GTEST = .FALSE.
         X(IVRNO) = X(IVRNO) - DELX
         !print *, 'xplus1: ', x(ivrno)
         CALL POTENTIAL(X, G, FM, GTEST, STEST)
         !print *, 'xplus2: ', x(ivrno)

         if (rigidinit.and.atomrigidcoordt) then
            call transformctorigid(x, tmpcoords)
            x(1:degfreedoms) = tmpcoords(1:degfreedoms)
            x(degfreedoms+1:3*natoms) = 0.0d0
            atomrigidcoordt = .false.
         endif

         X(IVRNO) = X(IVRNO) + 2.D0*DELX
         !print *, 'xminus1: ', x(ivrno)
         CALL POTENTIAL(X, G, FP, GTEST, STEST)
         !print *, 'xminus2: ', x(ivrno)
     
         if (rigidinit.and.atomrigidcoordt) then
            call transformctorigid(x, tmpcoords)
            x(1:degfreedoms) = tmpcoords(1:degfreedoms)
            x(degfreedoms+1:3*natoms) = 0.0d0
            atomrigidcoordt = .false. 
         endif
 
         GTEST = .TRUE.
         X(IVRNO) = X(IVRNO) - DELX
         CALL POTENTIAL(X, G, ENERGY, GTEST, STEST)
         DFN = (FP - FM) / (2.D0*DELX)
         IF (ABS(DFN).LT.1.0D-10) DFN = 0.D0
         DFA = G(IVRNO)

         WRITE(*, *) 'Gradient numerical  = ', DFN
         WRITE(*, *) 'Gradient analytical = ', DFA

         IF (ABS(DFN - DFA) > ERRLIM) WRITE(*, *) IVRNO, DFN, DFA, ABS(DFN-DFA)
      ENDDO

      ! dj337: check lattice derivatives
      if (boxderivt) then
         doff = 3
         if (.not.ortho) doff = 6
         DO IVRNO = 1, doff

            if (rigidinit.and.atomrigidcoordt) then
               call transformctorigid(x, tmpcoords)
               x(1:degfreedoms) = tmpcoords(1:degfreedoms)
               atomrigidcoordt = .false.
            endif

            WRITE(*, *) 'Box parameter ', IVRNO

            GTEST = .FALSE.
            BOX_PARAMS(IVRNO) = BOX_PARAMS(IVRNO) - DELX
            CALL POTENTIAL(X, G, FM, GTEST, STEST)
            WRITE(*, *) 'Energy minus = ', FM

            BOX_PARAMS(IVRNO) = BOX_PARAMS(IVRNO) + 2.D0*DELX
            CALL POTENTIAL(X, G, FP, GTEST, STEST)
            WRITE(*, *) 'Energy plus  = ', FP

            GTEST = .TRUE.
            BOX_PARAMS(IVRNO) = BOX_PARAMS(IVRNO) - DELX
            CALL POTENTIAL(X, G, ENERGY, GTEST, STEST)
            DFN = (FP - FM) / (2.D0*DELX)
            IF (ABS(DFN).LT.1.0D-10) DFN = 0.D0
            DFA = BOX_PARAMSGRAD(IVRNO)

            WRITE(*, *) 'Box gradient numerical  = ', DFN
            WRITE(*, *) 'Box gradient analytical = ', DFA

            IF (ABS(DFN - DFA) > ERRLIM) WRITE(*, *) IVRNO, DFN, DFA, ABS(DFN-DFA)

         ENDDO
      endif

      ELSE IF (CHECKDID == 2) THEN

         IF (.NOT. ALLOCATED(HESS)) ALLOCATE(HESS(3*NATOMS,3*NATOMS))

         DO IVRNO1 = 1, 3*NATOMS
            DO IVRNO2 = 1, 3*NATOMS
               WRITE(*,*) IVRNO1, IVRNO2
               X(IVRNO1) = X(IVRNO1) - DELX
               CALL POTENTIAL (X,G,ENERGY,.TRUE.,.FALSE.)
               FM   = G(IVRNO2)
!              WRITE(*, *) 'Gradient minus = ', FM

               X(IVRNO1) = X(IVRNO1) + 2.D0*DELX
               CALL POTENTIAL (X,G,ENERGY,.TRUE.,.FALSE.)
               FP   = G(IVRNO2)
!              WRITE(*, *) 'Gradient plus = ', FP

               X(IVRNO1) = X(IVRNO1) - DELX
               CALL POTENTIAL (X,G,ENERGY,.TRUE.,.TRUE.)
               DFN  = (FP - FM) / (2.D0*DELX)
               DFA  = HESS(IVRNO1,IVRNO2)

               WRITE(*, *) 'Hessian numerical  = ', DFN
               WRITE(*, *) 'Hessian analytical = ', DFA

               IF (ABS(DFN - DFA) > ERRLIM) WRITE(*,*) 'Error:', IVRNO1, IVRNO2, DFN, DFA, ABS(DFN-DFA)
            ENDDO
         ENDDO
      ENDIF

      STOP

      END SUBROUTINE CHECKD
