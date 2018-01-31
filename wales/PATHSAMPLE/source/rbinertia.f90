      SUBROUTINE RBINERTIA(Q, ITX, ITY, ITZ)

      USE COMMONS, ONLY: NATOMS, NRBSITES, NTSITES, RBSITE, SITEMASS

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J4, J5
      DOUBLE PRECISION :: Q(3*NATOMS), ITX, ITY, ITZ, MSRBST(NTSITES)
      DOUBLE PRECISION :: X(3*NTSITES), R(3), P(3), RM(3,3), IT(3,3), CMX, CMY, CMZ, VEC(3,3), MASST

      J5 = 0

      DO J1 = 1, NATOMS/2

         J2   = 3*J1
         R(:) = Q(J2-2:J2)
         J2   = 3*NATOMS/2 + J2
         P(:) = Q(J2-2:J2)

         CALL ROTMAT(P, RM)

         DO J3 = 1, NRBSITES
            J4          = 3*((J1-1)*NRBSITES + J3)
            X(J4-2:J4)  = R(:) + MATMUL(RM(:,:), RBSITE(J3,:))
            J5 = J5 + 1
            MSRBST(J5) = SITEMASS(J3)
         ENDDO

      ENDDO

      CMX   = 0.0D0
      CMY   = 0.0D0
      CMZ   = 0.0D0
      MASST = 0.0D0

      DO J1 = 1, NTSITES
         CMX = CMX + X(3*(J1-1)+1)*MSRBST(J1)
         CMY = CMY + X(3*(J1-1)+2)*MSRBST(J1)
         CMZ = CMZ + X(3*(J1-1)+3)*MSRBST(J1)
         MASST = MASST + MSRBST(J1)
      ENDDO
      CMX = CMX/MASST
      CMY = CMY/MASST
      CMZ = CMZ/MASST
      DO J1 = 1, NTSITES
         X(3*(J1-1)+1) = X(3*(J1-1)+1) - CMX
         X(3*(J1-1)+2) = X(3*(J1-1)+2) - CMY
         X(3*(J1-1)+3) = X(3*(J1-1)+3) - CMZ
      ENDDO
      DO J1 = 1, 3
         DO J2 =1, 3
            IT(J1,J2) = 0.0D0
            DO J3 = 1, NTSITES
               IT(J1,J2) = IT(J1,J2) - X(3*(J3-1)+J1)*X(3*(J3-1)+J2)*MSRBST(J3)
            ENDDO
            IF (J1 == J2) THEN
               DO J3 = 1, NTSITES
                  IT(J1,J2) = IT(J1,J2) + (X(3*(J3-1)+1)**2 + X(3*(J3-1)+2)**2 + X(3*(J3-1)+3)**2)*MSRBST(J3)
               ENDDO
            ENDIF
         ENDDO
      ENDDO

      CALL EIG(IT,VEC,3,3,0)

      ITX = IT(1,1)
      ITY = IT(2,2)
      ITZ = IT(3,3)

      DO J1 = 1, NTSITES
         X(3*(J1-1)+1) = X(3*(J1-1)+1) + CMX
         X(3*(J1-1)+2) = X(3*(J1-1)+2) + CMY
         X(3*(J1-1)+3) = X(3*(J1-1)+3) + CMZ
      ENDDO

      END SUBROUTINE RBINERTIA
