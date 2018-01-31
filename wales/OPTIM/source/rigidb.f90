!
!     OPTIM is free software; you can redistribute it and/or modify
!     it under the terms of the GNU General Public License as published by
!     the Free Software Foundation; either version 2 of the License, or
!     (at your option) any later version.
!
!     OPTIM is distributed in the hope that it will be useful,
!     but WITHOUT ANY WARRANTY; without even the implied warranty of
!     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!     GNU General Public License for more details.
!
!     You should have received a copy of the GNU General Public License
!     along with this program; if not, write to the Free Software
!     Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!     ----------------------------------------------------------------------------------------------

      SUBROUTINE ROTMAT(P, RM)

      INTEGER          :: K
      DOUBLE PRECISION :: P(3), PN(3), THETA, THETA2, CT, ST, I3(3,3), E(3,3), RM(3,3)

      I3(:,:) = 0.D0
      DO K = 1, 3
         I3(K,K) = 1.D0
      ENDDO

      THETA2 = DOT_PRODUCT(P,P)

      IF (THETA2 < 1.D-12) THEN

         RM(:,:)     = I3(:,:)
! vr274> first order corrections to rotation matrix
! 11/1/12
         RM(1,2) = -P(3)
         RM(2,1) = P(3)
         RM(1,3) = P(2)
         RM(3,1) = -P(2)
         RM(2,3) = -P(1)
         RM(3,2) = P(1)

      ELSE
         THETA   = SQRT(THETA2)
         CT      = COS(THETA)
         ST      = SIN(THETA)
         THETA   = 1.D0/THETA
         PN(:)   = THETA*P(:)
         E(:,:)  = 0.D0
         E(1,2)  = -PN(3)
         E(1,3)  =  PN(2)
         E(2,3)  = -PN(1)
         E(2,1)  = -E(1,2)
         E(3,1)  = -E(1,3)
         E(3,2)  = -E(2,3)
         RM      = I3(:,:) + (1.D0-CT)*MATMUL(E(:,:),E(:,:)) + ST*E(:,:)
      ENDIF

      END SUBROUTINE ROTMAT

!     --------------------------------------------------------------------------

      SUBROUTINE RMDRVT(P, RM, DRM1, DRM2, DRM3, GTEST)

      IMPLICIT NONE

      DOUBLE PRECISION :: P(3), PN(3), THETA, THETA2, THETA3, CT, ST, I3(3,3), E(3,3), ESQ(3,3)
      DOUBLE PRECISION :: DE1(3,3), DE2(3,3), DE3(3,3), RM(3,3), DRM1(3,3), DRM2(3,3), DRM3(3,3)
      LOGICAL          :: GTEST

      I3(:,:) = 0.D0
      I3(1,1) = 1.D0; I3(2,2) = 1.D0; I3(3,3) = 1.D0

      THETA2  = DOT_PRODUCT(P,P)

      IF (THETA2 == 0.D0) THEN

         RM(:,:) = I3(:,:)

! vr274> first order corrections to rotation matrix
! 11/1/12
         RM(1,2) = -P(3)
         RM(2,1) = P(3)
         RM(1,3) = P(2)
         RM(3,1) = -P(2)
         RM(2,3) = -P(1)
         RM(3,2) = P(1)

!         IF (.NOT. GTEST) RETURN

! hk286 - now up to the linear order in theta
! 11/1/12
         E(:,:)    = 0.D0
         E(1,1)    = 0.0D0
         E(1,2)    = P(2)
         E(1,3)    = P(3)
         E(2,1)    = P(2)
         E(2,2)    = -2.0D0*P(1)
         E(2,3)    = -2.0D0
         E(3,1)    = P(3)
         E(3,2)    = 2.0D0
         E(3,3)    = -2.0D0*P(1)
         DRM1(:,:) = 0.5D0*E(:,:)

         E(:,:)    = 0.D0
         E(1,1)    = -2.0D0*P(2)
         E(1,2)    = P(1)
         E(1,3)    = 2.0D0
         E(2,1)    = P(1)
         E(2,2)    = 0.0D0
         E(2,3)    = P(3)
         E(3,1)    = -2.0D0
         E(3,2)    = P(3)
         E(3,3)    = -2.0D0*P(2)
         DRM2(:,:) = 0.5D0*E(:,:)

         E(:,:)    = 0.D0
         E(1,1)    = -2.0D0*P(3)
         E(1,2)    = -2.0D0
         E(1,3)    = P(1)
         E(2,1)    = 2.0D0
         E(2,2)    = -2.0D0*P(3)
         E(2,3)    = P(2)
         E(3,1)    = P(1)
         E(3,2)    = P(2)
         E(3,3)    = 0.0D0
         DRM3(:,:) = 0.5D0*E(:,:)

      ELSE
        
         THETA   = SQRT(THETA2)
         CT      = COS(THETA)
         ST      = SIN(THETA)
         THETA3  = 1.D0/(THETA2*THETA)
         THETA   = 1.D0/THETA
         PN(:)   = THETA*P(:)
         E(:,:)  = 0.D0
         E(1,2)  = -PN(3)
         E(1,3)  =  PN(2)
         E(2,3)  = -PN(1) 
         E(2,1)  = -E(1,2)
         E(3,1)  = -E(1,3)
         E(3,2)  = -E(2,3)

         ESQ     = MATMUL(E,E)
         RM      = I3(:,:) + (1.D0-CT)*ESQ(:,:) + ST*E(:,:)

!         IF (.NOT. GTEST) RETURN

         DE1(:,:) = 0.D0  
         DE1(1,2) = P(3)*P(1)*THETA3
         DE1(1,3) = -P(2)*P(1)*THETA3
         DE1(2,3) = -(THETA - P(1)*P(1)*THETA3)
         DE1(2,1) = -DE1(1,2)
         DE1(3,1) = -DE1(1,3)
         DE1(3,2) = -DE1(2,3)

         DE2(:,:) = 0.D0
         DE2(1,2) = P(3)*P(2)*THETA3
         DE2(1,3) = THETA - P(2)*P(2)*THETA3
         DE2(2,3) = P(1)*P(2)*THETA3
         DE2(2,1) = -DE2(1,2)
         DE2(3,1) = -DE2(1,3)
         DE2(3,2) = -DE2(2,3)

         DE3(:,:) = 0.D0
         DE3(1,2) = -(THETA - P(3)*P(3)*THETA3)
         DE3(1,3) = -P(2)*P(3)*THETA3
         DE3(2,3) = P(1)*P(3)*THETA3
         DE3(2,1) = -DE3(1,2)
         DE3(3,1) = -DE3(1,3)
         DE3(3,2) = -DE3(2,3)

         DRM1(:,:) = ST*PN(1)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE1,E) + MATMUL(E,DE1)) &
                   + CT*PN(1)*E(:,:) + ST*DE1(:,:)

         DRM2(:,:) = ST*PN(2)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE2,E) + MATMUL(E,DE2)) &
                   + CT*PN(2)*E(:,:) + ST*DE2(:,:)

         DRM3(:,:) = ST*PN(3)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE3,E) + MATMUL(E,DE3)) &
                   + CT*PN(3)*E(:,:) + ST*DE3(:,:)

      ENDIF

      END SUBROUTINE RMDRVT

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RMDFAS(P, RM, DRM1, DRM2, DRM3, D2RM1, D2RM2, D2RM3, D2RI12, D2RI23, D2RI31, GTEST, STEST)

      IMPLICIT NONE

      INTEGER          :: K
      DOUBLE PRECISION :: P(3), PN(3), THETA, THETA2, THETA3, CT, ST, E(3,3), ESQ(3,3), I3(3,3) 
      DOUBLE PRECISION :: DE1(3,3), DE2(3,3), DE3(3,3), RM(3,3), DRM1(3,3), DRM2(3,3), DRM3(3,3)
      DOUBLE PRECISION :: D2E1(3,3), D2E2(3,3), D2E3(3,3), D2E12(3,3), D2E23(3,3), D2E31(3,3)
      DOUBLE PRECISION :: D2RM1(3,3), D2RM2(3,3), D2RM3(3,3)
      DOUBLE PRECISION :: D2RI12(3,3), D2RI23(3,3), D2RI31(3,3)
      DOUBLE PRECISION :: FCTR, FCTRSQ1, FCTRSQ2, FCTRSQ3
      LOGICAL          :: GTEST, STEST

      I3(:,:) = 0.D0
      I3(1,1) = 1.D0; I3(2,2) = 1.D0; I3(3,3) = 1.D0

!      RM(3,3) = 0.D0; DRM1(3,3) = 0.D0; DRM2(3,3) = 0.D0; DRM3(3,3) = 0.D0 
!      D2RM1(3,3) = 0.D0; D2RM2(3,3) = 0.D0; D2RM3(3,3) = 0.D0
!      D2RI12(3,3) = 0.D0; D2RI23(3,3) = 0.D0; D2RI31(3,3) = 0.D0

      THETA2  = DOT_PRODUCT(P,P)

      IF (THETA2 < 1.D-12) THEN

         RM(:,:)     = I3(:,:)
! vr274> first order corrections to rotation matrix
! 11/1/12
         RM(1,2) = -P(3)
         RM(2,1) = P(3)
         RM(1,3) = P(2)
         RM(3,1) = -P(2)
         RM(2,3) = -P(1)
         RM(3,2) = P(1)

         IF (.NOT. GTEST .AND. .NOT. STEST) RETURN

! hk286 - now up to the linear order in theta
! 11/1/12
         E(:,:)    = 0.D0
         E(1,1)    = 0.0D0
         E(1,2)    = P(2)
         E(1,3)    = P(3)
         E(2,1)    = P(2)
         E(2,2)    = -2.0D0*P(1)
         E(2,3)    = -2.0D0
         E(3,1)    = P(3)
         E(3,2)    = 2.0D0
         E(3,3)    = -2.0D0*P(1)
         DRM1(:,:) = 0.5D0*E(:,:)

         E(:,:)    = 0.D0
         E(1,1)    = -2.0D0*P(2)
         E(1,2)    = P(1)
         E(1,3)    = 2.0D0
         E(2,1)    = P(1)
         E(2,2)    = 0.0D0
         E(2,3)    = P(3)
         E(3,1)    = -2.0D0
         E(3,2)    = P(3)
         E(3,3)    = -2.0D0*P(2)
         DRM2(:,:) = 0.5D0*E(:,:)

         E(:,:)    = 0.D0
         E(1,1)    = -2.0D0*P(3)
         E(1,2)    = -2.0D0
         E(1,3)    = P(1)
         E(2,1)    = 2.0D0
         E(2,2)    = -2.0D0*P(3)
         E(2,3)    = P(2)
         E(3,1)    = P(1)
         E(3,2)    = P(2)
         E(3,3)    = 0.0D0
         DRM3(:,:) = 0.5D0*E(:,:)

         IF (.NOT. STEST) RETURN

! Fix bug and now to linear order in theta
! 11/1/12
         D2RM1(:,:) = 0.0D0
         D2RM1(1,1) = 0.0D0
         D2RM1(1,2) =  P(3)
         D2RM1(1,3) = -P(2)
         D2RM1(2,1) = -P(3)
         D2RM1(2,2) = -3.0D0
         D2RM1(2,3) =  3.0D0*P(1)
         D2RM1(3,1) =  P(2)
         D2RM1(3,2) = -3.0D0*P(1)
         D2RM1(3,3) = -3.0D0
         D2RM1(:,:) = D2RM1(:,:)/3.0D0

         D2RM2(:,:) = 0.0D0
         D2RM2(1,1) = -3.0D0
         D2RM2(1,2) =  P(3)
         D2RM2(1,3) = -3.0D0*P(2)
         D2RM2(2,1) = -P(3)
         D2RM2(2,2) =  0.0D0
         D2RM2(2,3) =  P(1)
         D2RM2(3,1) =  3.0D0*P(2)
         D2RM2(3,2) = -P(1)
         D2RM2(3,3) = -3.0D0
         D2RM2(:,:) = D2RM2(:,:)/3.0D0

         D2RM3(:,:) = 0.0D0
         D2RM3(1,1) = -3.0D0
         D2RM3(1,2) =  3.0D0*P(3)
         D2RM3(1,3) = -P(2)
         D2RM3(2,1) = -3.0D0*P(3)
         D2RM3(2,2) = -3.0D0
         D2RM3(2,3) =  P(1)
         D2RM3(3,1) =  P(2)
         D2RM3(3,2) = -P(1)
         D2RM3(3,3) =  0.0D0
         D2RM3(:,:) = D2RM3(:,:)/3.0D0

         D2RI12(:,:) = 0.0D0
         D2RI12(1,1) = 0.0D0
         D2RI12(1,2) = 0.5D0
         D2RI12(1,3) = -P(1)/3.0D0
         D2RI12(2,1) = 0.5D0
         D2RI12(2,2) = 0.0D0
         D2RI12(2,3) = P(2)/3.0D0
         D2RI12(3,1) = P(1)/3.0D0
         D2RI12(3,2) = -P(2)/3.0D0
         D2RI12(3,3) = 0.0D0

         D2RI23(:,:) = 0.0D0
         D2RI23(1,1) = 0.0D0
         D2RI23(1,2) = P(2)/3.0D0
         D2RI23(1,3) = -P(3)/3.0D0
         D2RI23(2,1) = -P(2)/3.0D0
         D2RI23(2,2) = 0.0D0
         D2RI23(2,3) = 0.5D0
         D2RI23(3,1) = P(3)/3.0D0
         D2RI23(3,2) = 0.5D0
         D2RI23(3,3) = 0.0D0

         D2RI31(:,:) = 0.0D0
         D2RI31(1,1) = 0.0D0
         D2RI31(1,2) = P(1)/3.0D0
         D2RI31(1,3) = 0.5D0
         D2RI31(2,1) = -P(1)/3.0D0
         D2RI31(2,2) = 0.0D0
         D2RI31(2,3) = P(3)/3.0D0
         D2RI31(3,1) = 0.5D0
         D2RI31(3,2) = -P(3)/3.0D0
         D2RI31(3,3) = 0.0D0

      ELSE

         THETA   = SQRT(THETA2)
         CT      = COS(THETA)
         ST      = SIN(THETA)
         THETA3  = 1.D0/(THETA2*THETA)
         THETA   = 1.D0/THETA
         PN(:)   = THETA*P(:)
         E(:,:)  = 0.D0
         E(1,2)  = -PN(3)
         E(1,3)  =  PN(2)
         E(2,3)  = -PN(1)
         E(2,1)  = -E(1,2)
         E(3,1)  = -E(1,3)
         E(3,2)  = -E(2,3)

         ESQ     = MATMUL(E(:,:),E(:,:))
         RM      = I3(:,:) + (1.D0-CT)*ESQ(:,:) + ST*E(:,:)

         IF (.NOT. GTEST .AND. .NOT. STEST) RETURN

         DE1(:,:) = 0.D0
         DE1(1,2) = P(3)*P(1)*THETA3
         DE1(1,3) = -P(2)*P(1)*THETA3
         DE1(2,3) = -(THETA - P(1)*P(1)*THETA3)
         DE1(2,1) = -DE1(1,2)
         DE1(3,1) = -DE1(1,3)
         DE1(3,2) = -DE1(2,3)

         DE2(:,:) = 0.D0
         DE2(1,2) = P(3)*P(2)*THETA3
         DE2(1,3) = THETA - P(2)*P(2)*THETA3
         DE2(2,3) = P(1)*P(2)*THETA3
         DE2(2,1) = -DE2(1,2)
         DE2(3,1) = -DE2(1,3)
         DE2(3,2) = -DE2(2,3)

         DE3(:,:) = 0.D0
         DE3(1,2) = -(THETA - P(3)*P(3)*THETA3)
         DE3(1,3) = -P(2)*P(3)*THETA3
         DE3(2,3) = P(1)*P(3)*THETA3
         DE3(2,1) = -DE3(1,2)
         DE3(3,1) = -DE3(1,3)
         DE3(3,2) = -DE3(2,3)

         DRM1(:,:) = ST*PN(1)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE1,E) + MATMUL(E,DE1)) &
                   + CT*PN(1)*E(:,:) + ST*DE1(:,:)

         DRM2(:,:) = ST*PN(2)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE2,E) + MATMUL(E,DE2)) &
                   + CT*PN(2)*E(:,:) + ST*DE2(:,:)

         DRM3(:,:) = ST*PN(3)*ESQ(:,:) + (1.D0-CT)*(MATMUL(DE3,E) + MATMUL(E,DE3)) &
                   + CT*PN(3)*E(:,:) + ST*DE3(:,:)

         IF (.NOT. STEST) RETURN

         FCTRSQ1 = PN(1)*PN(1) 
         FCTRSQ2 = PN(2)*PN(2)
         FCTRSQ3 = PN(3)*PN(3)

         D2E1(:,:) =  0.D0
         FCTR      =  (1.D0 - 3.D0*PN(1)*PN(1))*THETA3
         D2E1(1,2) =  P(3)*FCTR
         D2E1(1,3) = -P(2)*FCTR
         D2E1(2,3) =  P(1)*FCTR + 2.D0*P(1)*THETA3
         D2E1(2,1) = -D2E1(1,2)
         D2E1(3,1) = -D2E1(1,3)
         D2E1(3,2) = -D2E1(2,3)

         D2E2(:,:) =  0.D0
         FCTR      =  (1.D0 - 3.D0*PN(2)*PN(2))*THETA3
         D2E2(1,2) =  P(3)*FCTR
         D2E2(1,3) = -P(2)*FCTR - 2.D0*P(2)*THETA3
         D2E2(2,3) =  P(1)*FCTR 
         D2E2(2,1) = -D2E2(1,2)
         D2E2(3,1) = -D2E2(1,3)
         D2E2(3,2) = -D2E2(2,3)

         D2E3(:,:) =  0.D0
         FCTR      =  (1.D0 - 3.D0*PN(3)*PN(3))*THETA3
         D2E3(1,2) =  P(3)*FCTR + 2.D0*P(3)*THETA3
         D2E3(1,3) = -P(2)*FCTR
         D2E3(2,3) =  P(1)*FCTR 
         D2E3(2,1) = -D2E3(1,2)
         D2E3(3,1) = -D2E3(1,3)
         D2E3(3,2) = -D2E3(2,3)

         D2E12(:,:) =  0.D0
         D2E12(1,2) = -3.D0*PN(1)*PN(2)*PN(3)/THETA2 
         D2E12(1,3) = -PN(1)*(1.D0 - 3.D0*FCTRSQ2)/THETA2
         D2E12(2,3) =  PN(2)*(1.D0 - 3.D0*FCTRSQ1)/THETA2
         D2E12(2,1) = -D2E12(1,2)
         D2E12(3,1) = -D2E12(1,3)
         D2E12(3,2) = -D2E12(2,3)

         D2E23(:,:) =  0.D0
         D2E23(1,2) =  PN(2)*(1.D0 - 3.D0*FCTRSQ3)/THETA2
         D2E23(1,3) = -PN(3)*(1.D0 - 3.D0*FCTRSQ2)/THETA2
         D2E23(2,3) = -3.D0*PN(1)*PN(2)*PN(3)/THETA2
         D2E23(2,1) = -D2E23(1,2)
         D2E23(3,1) = -D2E23(1,3)
         D2E23(3,2) = -D2E23(2,3)

         D2E31(:,:) =  0.D0
         D2E31(1,2) =  PN(1)*(1.D0 - 3.D0*FCTRSQ3)/THETA2
         D2E31(1,3) =  3.D0*PN(1)*PN(2)*PN(3)/THETA2
         D2E31(2,3) =  PN(3)*(1.D0 - 3.D0*FCTRSQ1)/THETA2
         D2E31(2,1) = -D2E31(1,2)
         D2E31(3,1) = -D2E31(1,3)
         D2E31(3,2) = -D2E31(2,3)

         D2RM1(:,:)  = ST*PN(1)*(MATMUL(DE1(:,:),E(:,:)) + MATMUL(E(:,:),DE1(:,:))) &
                     + (CT*FCTRSQ1 - ST*FCTRSQ1*THETA + ST*THETA)*ESQ(:,:) &
                     + ST*PN(1)*(MATMUL(DE1(:,:),E(:,:)) + MATMUL(E(:,:),DE1(:,:))) &
                     + (1.D0-CT)*(2.D0*MATMUL(DE1(:,:),DE1(:,:)) + MATMUL(D2E1(:,:),E(:,:)) &
                     + MATMUL(E(:,:),D2E1(:,:))) + (- ST*FCTRSQ1 - CT*FCTRSQ1*THETA + CT*THETA)*E(:,:) &
                     + 2.D0*CT*PN(1)*DE1(:,:) + ST*D2E1(:,:)
         
         D2RM2(:,:)  = ST*PN(2)*(MATMUL(DE2(:,:),E(:,:)) + MATMUL(E(:,:),DE2(:,:))) &
                     + (CT*FCTRSQ2 - ST*FCTRSQ2*THETA + ST*THETA)*ESQ(:,:) &
                     + ST*PN(2)*(MATMUL(DE2(:,:),E(:,:)) + MATMUL(E(:,:),DE2(:,:))) &
                     + (1.D0-CT)*(2.D0*MATMUL(DE2(:,:),DE2(:,:)) + MATMUL(D2E2(:,:),E(:,:)) & 
                     + MATMUL(E(:,:),D2E2(:,:))) + (- ST*FCTRSQ2 - CT*FCTRSQ2*THETA + CT*THETA)*E(:,:) &
                     + 2.D0*CT*PN(2)*DE2(:,:) + ST*D2E2(:,:)

         D2RM3(:,:)  = ST*PN(3)*(MATMUL(DE3(:,:),E(:,:)) + MATMUL(E(:,:),DE3(:,:)))   &
                     + (CT*FCTRSQ3 - ST*FCTRSQ3*THETA + ST*THETA)*ESQ(:,:) &
                     + ST*PN(3)*(MATMUL(DE3(:,:),E(:,:)) + MATMUL(E(:,:),DE3(:,:)))   &
                     + (1.D0-CT)*(2.D0*MATMUL(DE3(:,:),DE3(:,:)) + MATMUL(D2E3(:,:),E(:,:)) &
                     + MATMUL(E(:,:),D2E3(:,:))) + (- ST*FCTRSQ3 - CT*FCTRSQ3*THETA + CT*THETA)*E(:,:) &
                     + 2.D0*CT*PN(3)*DE3(:,:) + ST*D2E3(:,:)

         D2RI12(:,:) = ST*PN(1)*(MATMUL(E(:,:),DE2(:,:)) + MATMUL(DE2(:,:),E(:,:))) &
                     + (CT*PN(1)*PN(2) - ST*PN(1)*PN(2)*THETA)*ESQ(:,:) &
                     + ST*PN(2)*(MATMUL(DE1(:,:),E(:,:)) + MATMUL(E(:,:),DE1(:,:))) &
                     + (1.D0 - CT)*(MATMUL(D2E12(:,:),E(:,:)) + MATMUL(DE1(:,:),DE2(:,:)) & 
                     + MATMUL(DE2(:,:),DE1(:,:)) + MATMUL(E(:,:),D2E12(:,:))) & 
                     - (ST*PN(1)*PN(2) + CT*PN(1)*PN(2)*THETA)*E(:,:) + CT*PN(1)*DE2(:,:) + CT*PN(2)*DE1(:,:) &
                     + ST*D2E12(:,:)

         D2RI23(:,:) = ST*PN(2)*(MATMUL(E(:,:),DE3(:,:)) + MATMUL(DE3(:,:),E(:,:))) &
                     + (CT*PN(2)*PN(3) - ST*PN(2)*PN(3)*THETA)*ESQ(:,:) &
                     + ST*PN(3)*(MATMUL(DE2(:,:),E(:,:)) + MATMUL(E(:,:),DE2(:,:))) &
                     + (1.D0 - CT)*(MATMUL(D2E23(:,:),E(:,:)) + MATMUL(DE2(:,:),DE3(:,:)) &
                     + MATMUL(DE3(:,:),DE2(:,:)) + MATMUL(E(:,:),D2E23(:,:))) &
                     - (ST*PN(2)*PN(3) + CT*PN(2)*PN(3)*THETA)*E(:,:) + CT*PN(2)*DE3(:,:) + CT*PN(3)*DE2(:,:) &
                     + ST*D2E23(:,:)

         D2RI31(:,:) = ST*PN(3)*(MATMUL(E(:,:),DE1(:,:)) + MATMUL(DE1(:,:),E(:,:))) &
                     + (CT*PN(3)*PN(1) - ST*PN(3)*PN(1)*THETA)*ESQ(:,:) &
                     + ST*PN(1)*(MATMUL(DE3(:,:),E(:,:)) + MATMUL(E(:,:),DE3(:,:))) &
                     + (1.D0 - CT)*(MATMUL(D2E31(:,:),E(:,:)) + MATMUL(DE3(:,:),DE1(:,:)) &
                     + MATMUL(DE1(:,:),DE3(:,:)) + MATMUL(E(:,:),D2E31(:,:))) &
                     - (ST*PN(3)*PN(1) + CT*PN(3)*PN(1)*THETA)*E(:,:) + CT*PN(3)*DE1(:,:) + CT*PN(1)*DE3(:,:) &
                     + ST*D2E31(:,:)

      ENDIF
 
      END SUBROUTINE RMDFAS

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE SHIFTRIGID (Q, NATOMS)

!     THIS SUBROUTINE SHIFTS THE 'ZERO' EIGENVALUES CORRESPONDING TO OVERALL TRANSLATION AND
!     ROTATION OF A SYSTEM OF (IDENTICAL) RIGID BODIES WITH C-O-M & ANGLE-AXIS COORDINATES.

      USE KEY
      USE MODHESS

      IMPLICIT NONE

      INTEGER            :: NATOMS, I, J, J1, J2
      DOUBLE PRECISION   :: Q(3*NATOMS), EV(3*NATOMS,6+NATOMS/2), NRMFCT(6+NATOMS/2)
      DOUBLE PRECISION   :: CMX, CMY, CMZ, THETA, THETA2, THETAH, DUMMY

!     INITIALIZE
      EV(:,:)   = 0.D0
      NRMFCT(:) = 0.D0
      CMX       = 0.D0
      CMY       = 0.D0
      CMZ       = 0.D0

      SHIFTED = .TRUE.
      IF (EFIELDT) THEN
          NZERO = 4
      ELSE
          NZERO = 6
      ENDIF
! bug fix for pgi compiler
      IF(.NOT.ALLOCATED(UNIAXARRAY)) ALLOCATE(UNIAXARRAY(NATOMS/2))

      DO I = 1, NATOMS/2

         J = 3*I

         CMX = CMX + Q(J-2)
         CMY = CMY + Q(J-1)
         CMZ = CMZ + Q(J)

      ENDDO

      CMX = CMX / FLOAT(NATOMS/2)
      CMY = CMY / FLOAT(NATOMS/2)
      CMZ = CMZ / FLOAT(NATOMS/2)

      DO I = 1, NATOMS/2

         J  = 3*I
         J1 = 3*NATOMS/2 + J

         THETA2 = DOT_PRODUCT(Q(J1-2:J1), Q(J1-2:J1))
         THETA  = DSQRT(THETA2)
         THETAH = 0.5D0*THETA

!     TRANSLATION ALONG X
         EV(J-2,1) = 1.D0
         NRMFCT(1) = NRMFCT(1) + 1.D0

!     TRANSLATION ALONG Y
         EV(J-1,2) = 1.D0
         NRMFCT(2) = NRMFCT(2) + 1.D0

!     TRANSLATION ALONG Z
         EV(J,3) = 1.D0
         NRMFCT(3) = NRMFCT(3) + 1.D0

!     IF EFIELD IS PRESENT, IT HAS TO BE ALONG THE Z-AXIS IN THE LAB FRAME BY CONVENTION.

         IF (THETA == 0.D0) THEN

!     ROTATION ABOUT Z
            EV(J-2,4)  = - Q(J-1) + CMY
            EV(J-1,4)  = Q(J-2) - CMX
            EV(J1,4)   = 1.D0
            NRMFCT(4)  = NRMFCT(4) + EV(J-2,4)**2 + EV(J-1,4)**2 + EV(J1-2,4)**2 + EV(J1-1,4)**2 + EV(J1,4)**2

            IF (.NOT. EFIELDT) THEN

!     ROTATION ABOUT X
               EV(J-1,5)  = - Q(J) + CMZ
               EV(J,5)    = Q(J-1) - CMY
               EV(J1-2,5) = 1.D0
               NRMFCT(5)  = NRMFCT(5) + EV(J-1,5)**2 + EV(J,5)**2 + EV(J1-2,5)**2 + EV(J1-1,5)**2 + EV(J1,5)**2

!     ROTATION ABOUT Y
               EV(J-2,6)  = Q(J) - CMZ
               EV(J,6)    = - Q(J-2) + CMX
               EV(J1-1,6) = 1.D0
               NRMFCT(6)  = NRMFCT(6) + EV(J-2,6)**2 + EV(J,6)**2 + EV(J1-2,6)**2 + EV(J1-1,6)**2 + EV(J1,6)**2

            ENDIF

         ELSE

!     ROTATION ABOUT Z
           EV(J-2,4)  = - Q(J-1) + CMY
           EV(J-1,4)  = Q(J-2) - CMX
           EV(J1-2,4) = - 0.5D0*Q(J1-1) + Q(J1)*Q(J1-2)/THETA2 - 0.5D0*Q(J1)*Q(J1-2)/(THETA*TAN(THETAH))
           EV(J1-1,4) = 0.5D0*Q(J1-2) + Q(J1)*Q(J1-1)/THETA2 - 0.5D0*Q(J1)*Q(J1-1)/(THETA*TAN(THETAH))
           EV(J1,4)   = THETAH/TAN(THETAH) + Q(J1)**2/THETA2 - 0.5D0*Q(J1)**2/(THETA*TAN(THETAH))
           NRMFCT(4)  = NRMFCT(4) + EV(J-2,4)**2 + EV(J-1,4)**2 + EV(J1-2,4)**2 + EV(J1-1,4)**2 + EV(J1,4)**2
             
           IF (.NOT. EFIELDT) THEN

!     ROTATION ABOUT X
              EV(J-1,5)  = - Q(J) + CMZ
              EV(J,5)    = Q(J-1) - CMY
              EV(J1-2,5) = THETAH/TAN(THETAH) + Q(J1-2)**2/THETA2 - 0.5D0*Q(J1-2)**2/(THETA*TAN(THETAH))
              EV(J1-1,5) = - 0.5D0*Q(J1) + Q(J1-2)*Q(J1-1)/THETA2 - 0.5D0*Q(J1-2)*Q(J1-1)/(THETA*TAN(THETAH))
              EV(J1,5)   = 0.5D0*Q(J1-1) + Q(J1-2)*Q(J1)/THETA2 - 0.5D0*Q(J1-2)*Q(J1)/(THETA*TAN(THETAH))
              NRMFCT(5)  = NRMFCT(5) + EV(J-1,5)**2 + EV(J,5)**2 + EV(J1-2,5)**2 + EV(J1-1,5)**2 + EV(J1,5)**2

!     ROTATION ABOUT Y
              EV(J-2,6)  = Q(J) - CMZ
              EV(J,6)    = - Q(J-2) + CMX
              EV(J1-2,6) = 0.5D0*Q(J1) + Q(J1-1)*Q(J1-2)/THETA2 - 0.5D0*Q(J1-1)*Q(J1-2)/(THETA*TAN(THETAH))
              EV(J1-1,6) = THETAH/TAN(THETAH) + Q(J1-1)**2/THETA2 - 0.5D0*Q(J1-1)**2/(THETA*TAN(THETAH))
              EV(J1,6)   = - 0.5D0*Q(J1-2) + Q(J1-1)*Q(J1)/THETA2 - 0.5D0*Q(J1-1)*Q(J1)/(THETA*TAN(THETAH))
              NRMFCT(6)  = NRMFCT(6) + EV(J-2,6)**2 + EV(J,6)**2 + EV(J1-2,6)**2 + EV(J1-1,6)**2 + EV(J1,6)**2

            ENDIF

         ENDIF

         IF (GBT.OR.GBDT.OR.(PYGT.AND.UNIAXT).OR.STOCKAAT.OR.((PYT.OR.PYGPERIODICT.OR.PYBINARYT.OR.MULTISITEPYT).AND.UNIAXT)) THEN

!     ROTATION ABOUT THE SYMMETRY AXIS
            IF (THETA == 0.D0) THEN

               EV(J1,NZERO+I)  = 1.D0 
               NRMFCT(NZERO+I) = NRMFCT(NZERO+I) + EV(J1,NZERO+I)**2
            
            ELSE 
               EV(J1-2,NZERO+I) = 0.5D0*Q(J1-1) + Q(J1-2)*Q(J1)/THETA2 - 0.5D0*Q(J1-2)*Q(J1)/(THETA*TAN(THETAH))
               EV(J1-1,NZERO+I) = - 0.5D0*Q(J1-2) + Q(J1-1)*Q(J1)/THETA2 - 0.5D0*Q(J1-1)*Q(J1)/(THETA*TAN(THETAH))
               EV(J1,NZERO+I)   = THETAH*SIN(THETA) + Q(J1)*Q(J1)/THETA2 &
                                + 0.5D0*(THETA*COS(THETA)-Q(J1)*Q(J1)/THETA)/TAN(THETAH)
               NRMFCT(NZERO+I)  = NRMFCT(NZERO+I) + EV(J1-2,NZERO+I)**2 + EV(J1-1,NZERO+I)**2 + EV(J1,NZERO+I)**2  
            ENDIF
         ELSE IF (SANDBOXT.AND.UNIAXARRAY(I)) THEN   !sandbox potential with arbitrary indices of rotationally symmetric rigid bodies
            IF (THETA == 0.D0) THEN

               EV(J1,NZERO+I)  = 1.D0
               NRMFCT(NZERO+I) = NRMFCT(NZERO+I) + EV(J1,NZERO+I)**2

            ELSE
               EV(J1-2,NZERO+I) = 0.5D0*Q(J1-1) + Q(J1-2)*Q(J1)/THETA2 - 0.5D0*Q(J1-2)*Q(J1)/(THETA*TAN(THETAH))
               EV(J1-1,NZERO+I) = - 0.5D0*Q(J1-2) + Q(J1-1)*Q(J1)/THETA2 - 0.5D0*Q(J1-1)*Q(J1)/(THETA*TAN(THETAH))
               EV(J1,NZERO+I)   = THETAH*SIN(THETA) + Q(J1)*Q(J1)/THETA2 &
                                + 0.5D0*(THETA*COS(THETA)-Q(J1)*Q(J1)/THETA)/TAN(THETAH)
               NRMFCT(NZERO+I)  = NRMFCT(NZERO+I) + EV(J1-2,NZERO+I)**2 + EV(J1-1,NZERO+I)**2 + EV(J1,NZERO+I)**2
            ENDIF
         ENDIF

      ENDDO
! sf344> set the number of zeros properly for SANDBOX, if only some of the building blocks have rotational symmetry
      IF(UNIAXT) THEN 
         DO I=1,NATOMS/2
           IF (SANDBOXT.AND.UNIAXARRAY(I)) THEN
             NZERO = NZERO + 1
           END IF
         END DO
      END IF
      IF (GBT.OR.GBDT.OR.(PYGT.AND.UNIAXT).OR.STOCKAAT.OR.((PYT.OR.PYGPERIODICT.OR.PYBINARYT.OR.MULTISITEPYT).AND.UNIAXT)) THEN
         NZERO = NZERO + NATOMS/2
      ELSE
         NZERO = NZERO
      ENDIF

      DO J = 1, NZERO

         NRMFCT(J) = DSQRT(NRMFCT(J))
         EV(:,J)   = EV(:,J)/NRMFCT(J)

      ENDDO

!     GRAM-SCHMIDT ORTHOGONALISATION TO OBTAIN ORTHONORMAL ROTATIONAL EIGENVECTORS

      DO J = 4, NZERO

         DO J1 = 4, J-1

            EV(:,J) = EV(:,J) - DOT_PRODUCT(EV(:,J),EV(:,J1))*EV(:,J1)

         ENDDO

         EV(:,J) = EV(:,J) / DSQRT(DOT_PRODUCT(EV(:,J),EV(:,J)))

      ENDDO

      DO J1 = 1, 3*NATOMS

         DO J2 = 1, 3*NATOMS

            DO J = 1, NZERO 

               HESS(J2,J1) = HESS(J2,J1) + SHIFTV*EV(J2,J)*EV(J1,J)

            ENDDO

         ENDDO

      ENDDO

      END SUBROUTINE SHIFTRIGID
 
!     ----------------------------------------------------------------------------------------------

      SUBROUTINE ORTHOGRIGID (VEC1, Q, OTEST)

!     THIS SUBROUTINE ORTHOGONALISES VEC1 TO THE EIGENVECTORS CORRESPONDING TO OVERALL TRANSLATION
!     AND ROTATION OF A SYSTEM OF (IDENTICAL) RIGID BODIES WITH C-O-M & ANGLE-AXIS COORDINATES.

      USE COMMONS
      USE KEY

      IMPLICIT NONE

      INTEGER            :: I, J, J1
      DOUBLE PRECISION   :: Q(3*NATOMS), EV(3*NATOMS,6+NATOMS/2), NRMFCT(6+NATOMS/2), VEC1(3*NATOMS)
      DOUBLE PRECISION   :: CMX, CMY, CMZ, THETA, THETA2, THETAH, DUMMY
      LOGICAL            :: OTEST

!     INITIALIZE

      EV(:,:)   = 0.D0
      NRMFCT(:) = 0.D0
      CMX       = 0.D0
      CMY       = 0.D0
      CMZ       = 0.D0
      
      IF (EFIELDT) THEN
          NZERO = 4
      ELSE
          NZERO = 6
      ENDIF

      DO I = 1, NATOMS/2

         J = 3*I

         CMX = CMX + Q(J-2)
         CMY = CMY + Q(J-1)
         CMZ = CMZ + Q(J)

      ENDDO

      CMX = CMX / FLOAT(NATOMS/2)
      CMY = CMY / FLOAT(NATOMS/2)
      CMZ = CMZ / FLOAT(NATOMS/2)

      DO I = 1, NATOMS/2

         J  = 3*I   
         J1 = 3*NATOMS/2 + J

         THETA2 = DOT_PRODUCT(Q(J1-2:J1), Q(J1-2:J1))
         THETA  = DSQRT(THETA2)
         THETAH = 0.5D0*THETA

!     TRANSLATION ALONG X
         EV(J-2,1) = 1.D0
         NRMFCT(1) = NRMFCT(1) + 1.D0

!     TRANSLATION ALONG X
         EV(J-1,2) = 1.D0
         NRMFCT(2) = NRMFCT(2) + 1.D0

!     TRANSLATION ALONG X
         EV(J,3) = 1.D0
         NRMFCT(3) = NRMFCT(3) + 1.D0

!     IF EFIELD IS PRESENT, IT HAS TO BE ALONG THE Z-AXIS IN THE LAB FRAME BY CONVENTION.

         IF (THETA == 0.D0) THEN

!     ROTATION ABOUT Z
            EV(J-2,4)  = - Q(J-1) + CMY
            EV(J-1,4)  = Q(J-2) - CMX
            EV(J1,4)   = 1.D0
            NRMFCT(4)  = NRMFCT(4) + EV(J-2,4)**2 + EV(J-1,4)**2 + EV(J1-2,4)**2 + EV(J1-1,4)**2 + EV(J1,4)**2

            IF (.NOT. EFIELDT) THEN

!     ROTATION ABOUT X
              EV(J-1,5)  = - Q(J) + CMZ
              EV(J,5)    = Q(J-1) - CMY
              EV(J1-2,5) = 1.D0
              NRMFCT(5)  = NRMFCT(5) + EV(J-1,5)**2 + EV(J,5)**2 + EV(J1-2,5)**2 + EV(J1-1,5)**2 + EV(J1,5)**2

!     ROTATION ABOUT Y
              EV(J-2,6)  = Q(J) - CMZ
              EV(J,6)    = - Q(J-2) + CMX
              EV(J1-1,6) = 1.D0
              NRMFCT(6)  = NRMFCT(6) + EV(J-2,6)**2 + EV(J,6)**2 + EV(J1-2,6)**2 + EV(J1-1,6)**2 + EV(J1,6)**2

            ENDIF

         ELSE

!     ROTATION ABOUT Z
            EV(J-2,4)  = - Q(J-1) + CMY
            EV(J-1,4)  = Q(J-2) - CMX
            EV(J1-2,4) = - 0.5D0*Q(J1-1) + Q(J1)*Q(J1-2)/THETA2 - 0.5D0*Q(J1)*Q(J1-2)/(THETA*TAN(THETAH))
            EV(J1-1,4) = 0.5D0*Q(J1-2) + Q(J1)*Q(J1-1)/THETA2 - 0.5D0*Q(J1)*Q(J1-1)/(THETA*TAN(THETAH))
            EV(J1,4)   = THETAH/TAN(THETAH) + Q(J1)**2/THETA2 - 0.5D0*Q(J1)**2/(THETA*TAN(THETAH))
            NRMFCT(4)  = NRMFCT(4) + EV(J-2,4)**2 + EV(J-1,4)**2 + EV(J1-2,4)**2 + EV(J1-1,4)**2 + EV(J1,4)**2

            IF (.NOT. EFIELDT) THEN

!     ROTATION ABOUT X
               EV(J-1,5)  = - Q(J) + CMZ
               EV(J,5)    = Q(J-1) - CMY
               EV(J1-2,5) = THETAH/TAN(THETAH) + Q(J1-2)**2/THETA2 - 0.5D0*Q(J1-2)**2/(THETA*TAN(THETAH))
               EV(J1-1,5) = - 0.5D0*Q(J1) + Q(J1-2)*Q(J1-1)/THETA2 - 0.5D0*Q(J1-2)*Q(J1-1)/(THETA*TAN(THETAH))
               EV(J1,5)   = 0.5D0*Q(J1-1) + Q(J1-2)*Q(J1)/THETA2 - 0.5D0*Q(J1-2)*Q(J1)/(THETA*TAN(THETAH))
               NRMFCT(5)  = NRMFCT(5) + EV(J-1,5)**2 + EV(J,5)**2 + EV(J1-2,5)**2 + EV(J1-1,5)**2 + EV(J1,5)**2

!     ROTATION ABOUT Y
               EV(J-2,6)  = Q(J) - CMZ
               EV(J,6)    = - Q(J-2) + CMX
               EV(J1-2,6) = 0.5D0*Q(J1) + Q(J1-1)*Q(J1-2)/THETA2 - 0.5D0*Q(J1-1)*Q(J1-2)/(THETA*TAN(THETAH))
               EV(J1-1,6) = THETAH/TAN(THETAH) + Q(J1-1)**2/THETA2 - 0.5D0*Q(J1-1)**2/(THETA*TAN(THETAH))
               EV(J1,6)   = - 0.5D0*Q(J1-2) + Q(J1-1)*Q(J1)/THETA2 - 0.5D0*Q(J1-1)*Q(J1)/(THETA*TAN(THETAH))
               NRMFCT(6)  = NRMFCT(6) + EV(J-2,6)**2 + EV(J,6)**2 + EV(J1-2,6)**2 + EV(J1-1,6)**2 + EV(J1,6)**2

             ENDIF

         ENDIF

         IF (GBT.OR.GBDT.OR.(PYGT.AND.UNIAXT).OR.STOCKAAT.OR.((PYT.OR.PYGPERIODICT.OR.PYBINARYT.OR.MULTISITEPYT).AND.UNIAXT)) THEN

!     ROTATION ABOUT THE SYMMETRY AXIS
            IF (THETA == 0.D0) THEN

               EV(J1,NZERO+I)  = 1.D0
               NRMFCT(NZERO+I) = NRMFCT(NZERO+I) + EV(J1,NZERO+I)**2

            ELSE
               EV(J1-2,NZERO+I) = 0.5D0*Q(J1-1) + Q(J1-2)*Q(J1)/THETA2 - 0.5D0*Q(J1-2)*Q(J1)/(THETA*TAN(THETAH))
               EV(J1-1,NZERO+I) = - 0.5D0*Q(J1-2) + Q(J1-1)*Q(J1)/THETA2 - 0.5D0*Q(J1-1)*Q(J1)/(THETA*TAN(THETAH))
               EV(J1,NZERO+I)   = THETAH*SIN(THETA) + Q(J1)*Q(J1)/THETA2 &
                                + 0.5D0*(THETA*COS(THETA)-Q(J1)*Q(J1)/THETA)/TAN(THETAH)
               NRMFCT(NZERO+I)  = NRMFCT(NZERO+I) + EV(J1-2,NZERO+I)**2 + EV(J1-1,NZERO+I)**2 + EV(J1,NZERO+I)**2

            ENDIF

         ENDIF

      ENDDO

      IF (GBT.OR.GBDT.OR.(PYGT.AND.UNIAXT).OR.STOCKAAT.OR.((PYT.OR.PYGPERIODICT.OR.PYBINARYT.OR.MULTISITEPYT).AND.UNIAXT)) THEN
         NZERO = NZERO + NATOMS/2
      ELSE
         NZERO = NZERO
      ENDIF

      DO J = 1, NZERO  

         NRMFCT(J) = DSQRT(NRMFCT(J))
         EV(:,J)   = EV(:,J)/NRMFCT(J)

      ENDDO

!     GRAM-SCHMIDT ORTHOGONALISATION TO OBTAIN ORTHONORMAL ROTATIONAL EIGENVECTORS

      DO J = 4, NZERO  

         DO J1 = 4, J-1

            EV(:,J) = EV(:,J) - DOT_PRODUCT(EV(:,J),EV(:,J1))*EV(:,J1)

         ENDDO

         EV(:,J) = EV(:,J) / DSQRT(DOT_PRODUCT(EV(:,J),EV(:,J)))

      ENDDO

!     PROJECT TRANS/ROT SET OUT OF VEC1

      DO J = 1, NZERO 

         DUMMY   = DOT_PRODUCT(VEC1(:),EV(:,J))
         VEC1(:) = VEC1(:) - DUMMY*EV(:,J)

      ENDDO

      IF (OTEST) CALL VECNORM(VEC1,NOPT) ! NORMALIZE VEC1

      END SUBROUTINE ORTHOGRIGID

!     ----------------------------------------------------------------------------------------------
      SUBROUTINE GENRIGID_MINDIST_BULK(XA,XB,BOXLX,BOXLY,BOXLZ,DIST,DEBUG)
      ! Puts XA into its best alignment with XB with respect to translation and permutation
      ! assuming a system of rigid bodies described using the GENRIGID module, and using
      ! periodic boundary conditions.
      ! This is currently a rather crude approach to the problem.
      ! Alignment with respect to rotation is not yet implemented.

      ! There are two options for handling the permutational alignment, which are selected
      ! between using the following keywords.
      ! The default should be to use PERMDIST only, using a perm.allow file to specify
      ! atomic permutations caused by inner-molecular symmetry operations and/or
      ! permutations of rigid bodies. Note that the latter option may be very slow.
      ! If the system contains only one type of rigid body and no free atoms, then RBSYMT
      ! may optionally be used as well (see OPTIM documentation). This may be faster than
      ! the previous option. In this case, PERMDIST and a perm.allow file must be provided
      ! but will not be used: only inner-molecular permutations are considered.
      USE KEY, ONLY: RBSYMT, PERMDIST
      USE COMMONS, ONLY: NATOMS
      USE GENRIGID, ONLY: RB_DISTANCE, DEGFREEDOMS, NRIGIDBODY, TRANSFORMCTORIGID, TRANSFORMRIGIDTOC

      IMPLICIT NONE

      ! XA and XB are atomistic (cartesian) coordinates vectors, RA and RB are corresponding
      ! angle-axis coordinates vectors.
      ! DIST is the distance measure between the two structures.
      ! DUMMYGRAD1 and 2 are dummies for the RB_DIST subroutine - their values are never used.
      ! AVE and TEMP will be used to compute the average displacement of rigid bodies between
      ! the two structures, so XA/RA can be translated by AVE to bring them into alignment
      ! with the centre of mass for XB/RB.
      DOUBLE PRECISION, INTENT(INOUT) :: XA(3*NATOMS), XB(3*NATOMS)
      DOUBLE PRECISION, INTENT(OUT) :: DIST
      LOGICAL, INTENT(IN) :: DEBUG
      DOUBLE PRECISION, INTENT(IN) :: BOXLX, BOXLY, BOXLZ
      DOUBLE PRECISION :: RA(DEGFREEDOMS), RB(DEGFREEDOMS), DUMMYGRAD1(DEGFREEDOMS), DUMMYGRAD2(DEGFREEDOMS)
      DOUBLE PRECISION :: BOXL(3), AVE(3), TEMP(3)
      INTEGER J1, J2, NRB

      ! Periodic boundary box lengths
      BOXL(1) = BOXLX
      BOXL(2) = BOXLY
      BOXL(3) = BOXLZ
      NRB = NRIGIDBODY

      ! The general procedure:
      ! - Optimise the distance measure with respect to permutations
      ! - Compute the average vector between a molecule in RA and the corresponding molecule
      ! in RB. Translate RA along this vector to align them
      ! - If any molecules have crossed a periodic boundary, replace them with their image
      ! on the original side of the boundary. This will mean that at the end of the minimisation,
      ! some molecules sit outside the periodic box. But this is required for NEB to work correctly.
      ! - Optimise again with respect to permutations. This will usually not change
      ! the structure very much.

      ! Permutational alignment.
      IF (RBSYMT) THEN
!       Convert to rigid-body coordinates (needed by ALIGN_INNER)
        CALL TRANSFORMCTORIGID(XA, RA)
        CALL TRANSFORMCTORIGID(XB, RB)
!       Compute the initial distance between the two structures (needed by ALIGN_INNER)
        CALL RB_DISTANCE(DIST, RA, RB, DUMMYGRAD1, DUMMYGRAD2, .FALSE.)
        CALL ALIGN_INNER(NRB, RA, RB, DIST)
      ELSE IF (PERMDIST) THEN
        CALL GENRIGID_MINPERM(XA, XB)
!       Convert to rigid-body coordinates (needed for the next bit)
!       TRANSFORMCTORIGID involves a call to MINPERMDIST, which can change the angle-axis
!       vectors if PERMDIST=TRUE. So we temporarily set it to false.
        PERMDIST = .FALSE.
        CALL TRANSFORMCTORIGID(XA, RA)
        CALL TRANSFORMCTORIGID(XB, RB)
        PERMDIST = .TRUE.
      ENDIF

      ! Compute the average displacement between the structures
      AVE(:) = 0.D0
      DO J1=1,NRB
        ! Vector between equivalent molecules in the two structures
        TEMP(:) = RB(3*J1-2:3*J1) - RA(3*J1-2:3*J1)
        ! Correct for PBCs
        DO J2 = 1,3
            TEMP(J2) = TEMP(J2) - BOXL(J2)*ANINT(TEMP(J2)/BOXL(J2))
        ENDDO
        AVE(:) = AVE(:) + TEMP(:)
      ENDDO
      AVE(:) = AVE(:)/NRB

      DO J1=1,NRB
        ! Perform translation
        RA(3*J1-2:3*J1) = RA(3*J1-2:3*J1) + AVE(:)
        DO J2 = 1,3
            ! Put all molecules in the box
            RA(3*(J1-1)+J2) = RA(3*(J1-1)+J2) - BOXL(J2)*ANINT(RA(3*(J1-1)+J2)/BOXL(J2))
            ! If any molecules are now on different sides of a periodic boundary in the
            ! two different structures, move them back. This is needed for the NEB.
            RA(3*(J1-1)+J2) = RA(3*(J1-1)+J2) + BOXL(J2)*ANINT((RB(3*(J1-1)+J2)-RA(3*(J1-1)+J2))/BOXL(J2))
        ENDDO
      ENDDO


      IF (RBSYMT) THEN
        CALL ALIGN_INNER(NRB, RA, RB, DIST)
      ELSE IF (PERMDIST) THEN
        ! GENRIGID_MINPERM needs atomistic coordinates
        CALL TRANSFORMRIGIDTOC(1,NRIGIDBODY,XA,RA)
        CALL TRANSFORMRIGIDTOC(1,NRIGIDBODY,XB,RB)
        CALL GENRIGID_MINPERM(XA, XB)

        ! RB_DISTANCE needs angle-axis coordinates
        PERMDIST = .FALSE.
        CALL TRANSFORMCTORIGID(XA,RA)
        CALL TRANSFORMCTORIGID(XB,RB)
        PERMDIST = .TRUE.
      ENDIF

      ! Compute optimised distance
      CALL RB_DISTANCE(DIST, RA, RB, DUMMYGRAD1, DUMMYGRAD2, .FALSE.)

      RETURN

      END SUBROUTINE GENRIGID_MINDIST_BULK

!     ----------------------------------------------------------------------------------------------
    SUBROUTINE GENRIGID_MINPERM(XA, XB)
    ! Optimise the distance measure with respect to all the permutations in perm.allow
    USE COMMONS, ONLY: NATOMS
    ! smallest_rij returns the shortest (periodic) vector between two points.
    USE GENRIGID, ONLY: smallest_rij
    USE KEY, ONLY: NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, SETS

    IMPLICIT NONE

    ! XA, XB are atomistic coordinates vectors. TO_SWAP is a list of the atoms which need
    ! to be swapped by a particular permutation.
    DOUBLE PRECISION, INTENT(INOUT) :: XA(3*NATOMS), XB(3*NATOMS)
    DOUBLE PRECISION :: OLD_DIST, NEW_DIST, TEMP(3)
    INTEGER, DIMENSION(:,:), ALLOCATABLE :: TO_SWAP
    INTEGER :: DUMMY, J1, J2, J3, J4

    DUMMY = 0
    ! Step through all possible permutations
    DO J1 = 1, NPERMGROUP
        ! J2 and J3 both run over all the atoms involved in this permutation (avoiding
        ! double counting)
        DO J2 = DUMMY+1, DUMMY+NPERMSIZE(J1)
            DO J3 = J2+1, DUMMY+NPERMSIZE(J1)
                ! TO_SWAP has two columns, with corresponding elements in each column
                ! containing the indices of atoms which are swapped by this permutation.
                ALLOCATE(TO_SWAP(1:1+NSETS(J1),1:2))
                OLD_DIST = 0
                NEW_DIST = 0

                ! Put the main atom of the permutation into TO_SWAP
                TO_SWAP(1,1) = PERMGROUP(J2)
                TO_SWAP(1,2) = PERMGROUP(J3)
                ! Now put in any secondary atoms that also swap over.
                DO J4 = 1,NSETS(J1)
                    TO_SWAP(1+J4,1) = SETS(PERMGROUP(J2),J4)
                    TO_SWAP(1+J4,2) = SETS(PERMGROUP(J3),J4)
                ENDDO
                ! Step through all the swapping atoms and compute the old and new
                ! contributions of each pair to the distance measure (which in cartesian
                ! coordinates is simply the sum of the square displacements)
                DO J4 = 1,NSETS(J1)+1
                    TEMP(:) = smallest_rij(XA(3*TO_SWAP(J4,1)-2:3*TO_SWAP(J4,1)), &
                                         & XB(3*TO_SWAP(J4,1)-2:3*TO_SWAP(J4,1)))
                    OLD_DIST = OLD_DIST + TEMP(1)**2 + TEMP(2)**2 + TEMP(3)**2

                    TEMP(:) = smallest_rij(XA(3*TO_SWAP(J4,2)-2:3*TO_SWAP(J4,2)), &
                                         & XB(3*TO_SWAP(J4,2)-2:3*TO_SWAP(J4,2)))
                    OLD_DIST = OLD_DIST + TEMP(1)**2 + TEMP(2)**2 + TEMP(3)**2


                    TEMP(:) = smallest_rij(XA(3*TO_SWAP(J4,2)-2:3*TO_SWAP(J4,2)), &
                                         & XB(3*TO_SWAP(J4,1)-2:3*TO_SWAP(J4,1)))
                    NEW_DIST = NEW_DIST + TEMP(1)**2 + TEMP(2)**2 + TEMP(3)**2

                    TEMP(:) = smallest_rij(XA(3*TO_SWAP(J4,1)-2:3*TO_SWAP(J4,1)), &
                                         & XB(3*TO_SWAP(J4,2)-2:3*TO_SWAP(J4,2)))
                    NEW_DIST = NEW_DIST + TEMP(1)**2 + TEMP(2)**2 + TEMP(3)**2
                ENDDO
                ! If performing the permutation reduces the distance, keep it.
                IF (NEW_DIST .LT. OLD_DIST) THEN
                    DO J4 = 1,NSETS(J1)+1
                        ! Perform the swaps.
                        TEMP(:) = XA(3*TO_SWAP(J4,1)-2:3*TO_SWAP(J4,1))
                        XA(3*TO_SWAP(J4,1)-2:3*TO_SWAP(J4,1)) = XA(3*TO_SWAP(J4,2)-2:3*TO_SWAP(J4,2))
                        XA(3*TO_SWAP(J4,2)-2:3*TO_SWAP(J4,2)) = TEMP(:)
                    ENDDO
                ENDIF
                ! Deallocate ready for the next permutation group.
                DEALLOCATE(TO_SWAP)
            ENDDO
        ENDDO
        DUMMY = DUMMY + NPERMSIZE(J1)
    ENDDO

    RETURN
    END SUBROUTINE GENRIGID_MINPERM
!--------------------------------------------------------------------------------------------------
      ! An alternative permutational alignment using the keyword RBSYMT which is designed
      ! to go with the fully-rigid framework (RBAAT==.TRUE.)
      ! See the description of GENRIGID_MINDIST_BULK for a description.
      ! This routine has not been extensively tested.
      SUBROUTINE ALIGN_INNER(NRB, RA, RB, DIST)

      USE COMMONS, ONLY: NATOMS
      USE KEY, ONLY: NRBGROUP, RBOPS
      ! Respectively: Distance measure for angle-axis coordinates and convert AA
      ! coordinates to atomistic.
      USE GENRIGID, ONLY: RB_DISTANCE

      IMPLICIT NONE

      ! NRB: the number of rigid bodies in the system
      ! RA and RB: the angle-axis coordinate vectors of the two structures being aligned
      ! XB: an atomistic cartesian coordinate array for structure B
      ! DIST: the distance measure between the two poses
      ! THETAH: Half the rotation angle for a symmetry operation belonging to the point group of
      ! the rigid body (which has order NRBGROUP). The quaternions for the symmetry operations are
      ! stored in RBOPS.
      ! P and PVEC are temporary variables used to hold angle-axis vectors for a particular rigid body
      ! QTMP is a temporary quaternion, RTEMP1 is a temporary rotation matrix
      ! TRIAL and BESTA are temporary AA coordinate vectors. TRIAL_AT is an atomistic version of TRIAL
      ! TRIALDIST is a temporary distance variable used to compare trial permutations against the best current distance.
      ! DUMMYGRAD1 and 2 are just to pass in to RB_DISTANCE. They are never used.
      INTEGER, INTENT (IN) :: NRB
      DOUBLE PRECISION, INTENT (INOUT) :: RA(6*NRB), RB(6*NRB)
      DOUBLE PRECISION, INTENT (INOUT) :: DIST
      DOUBLE PRECISION :: THETAH, TRIALDIST
      DOUBLE PRECISION :: P(3), PVEC(3), QTMP(4)
      DOUBLE PRECISION :: TRIAL(6*NRB), TRIAL_AT(3*NATOMS), BESTA(6*NRB),  XB(6*NRB)
      DOUBLE PRECISION :: DUMMYGRAD1(6*NRB), DUMMYGRAD2(6*NRB)
      DOUBLE PRECISION :: RTEMP1(3,3)
      INTEGER :: J1, J2

      BESTA(1:6*NRB) = RA(1:6*NRB)

      DO J1=1,NRB
          ! Step through all the inner-molecular symmetry operations provided in rbsymops
          DO J2=1,NRBGROUP
              ! P is the (normalised) angle-axis vector for body J1.
              P(:) = RA(3*(NRB+J1-1)+1:3*(NRB+J1-1)+3)

              ! J2 = 1 corresponds to the identity. If J2 /= 1, perform the operation
              ! being considered.
              IF (J2 /= 1) THEN
                 ! Convert P to a rotation matrix
                 CALL ROTMAT(P,RTEMP1)
                 ! PVEC is a new angle-axis vector generated from P by rotation by one of
                 ! the symmetry operations of the rigid body.
                 PVEC(1:3) = MATMUL(RTEMP1,RBOPS(1:3,J2))
                 PVEC(1:3) = PVEC(1:3)/SQRT(DOT_PRODUCT(PVEC(1:3),PVEC(1:3)))

                 ! Construct the quaternion for symmetry operation being considered
                 THETAH    = 0.5D0*RBOPS(4,J2)
                 QTMP(1)   = COS(THETAH)
                 QTMP(2:4) = SIN(THETAH)*PVEC(1:3)
                 ! Convert this quaternion into the new angle-axis vector P.
                 CALL QROTAA(QTMP,P)
              ENDIF
              ! Copy the current coordinates and update the corresponding AA vector
              TRIAL(1:6*NRB) = RA(1:6*NRB)
              TRIAL(3*(NRB+J1-1)+1:3*(NRB+J1-1)+3) = P(1:3)

              ! Compute the distance between RB and the trial version of RA.
              CALL RB_DISTANCE(TRIALDIST, TRIAL, RB, DUMMYGRAD1, DUMMYGRAD2, .FALSE.)

              ! If the trial symmetry operation is an improvement, keep it.
              IF (TRIALDIST.LT.DIST) THEN
                  DIST=TRIALDIST
                  BESTA(1:6*NRB)=TRIAL(1:6*NRB)
              ENDIF
          ENDDO
          ! After looping through all symmetry operations, keep the one which gave the best distance
          ! (Which could be the identity operation)
          RA(1:6*NRB)=BESTA(1:6*NRB)
      ENDDO

      ! Modified values to be returned are Dist, RA (and technically RB, although this hasn't actually changed)

      END SUBROUTINE ALIGN_INNER

! --------------------------------------------------------------------------------------------------------
      SUBROUTINE RBMINDIST_BULK(RA,RB,NATOMS,DIST,Q2,DEBUG)
      !js850> This will find the distance between two structures RA and RB.
      !This routine will be with a bulk system of rigid bodies.
      !The distance will be the r.m.s. distance between the atoms.


      USE COMMONS, ONLY: NRBSITES, PARAM1, PARAM2, PARAM3
      USE KEY, ONLY: NTSITES, DBPT, DBPTDT, MSSTOCKT, STOCKAAT, EFIELDT, BULKT

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NATOMS
      DOUBLE PRECISION, INTENT(INOUT) :: RA(3*NATOMS)
      DOUBLE PRECISION, INTENT(INOUT) :: RB(3*NATOMS)
      DOUBLE PRECISION, INTENT(OUT) :: DIST, Q2(4)
      LOGICAL, INTENT(IN)          :: DEBUG
      DOUBLE PRECISION, ALLOCATABLE :: XA(:), XB(:)
      DOUBLE PRECISION dr, boxl(3)
      INTEGER J1, J2, NSIZE

      !write(*,*) "RBMINDIST_BULK>"

      BOXL(1) = PARAM1
      BOXL(2) = PARAM2
      BOXL(3) = PARAM3
! sf344 modification
!      BOXL(:) = 100.0D0
      NSIZE = NTSITES
      ALLOCATE(XA(3*NSIZE),XB(3*NSIZE))

!     CONVERT TO THE CARTESIAN CORDINATES FOR THE SITES

      CALL SITEPOS(RA,XA)

      CALL SITEPOS(RB,XB)

!     apply periodic boundary conditions

      DIST = 0.D0
      DO J1 = 1, NSIZE
         DO J2 = 1, 3
            DR = ( XA(3*(J1-1)+J2) - XB(3*(J1-1)+J2) )
            DR = DR - BOXL(J2)*ANINT( DR / BOXL(J2) )
            DIST = DIST + DR**2
         ENDDO
      ENDDO
      DIST = SQRT(DIST)

!     apply symmetry opperations on the angle axis coordinates so that the
!     cartesian distance between them is minimized.  This distance will be used
!     for the nudged elastic band
      do J1 = NATOMS/2+1,NATOMS
         J2 = 3*(J1-1)
         CALL AADISTANCE ( RA(J2+1:J2+3) , RB(J2+1:J2+3) )
      enddo


      DEALLOCATE(XA,XB)

      END SUBROUTINE RBMINDIST_BULK

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBMINDIST(RA,RB,NATOMS,DIST,Q2,DEBUG)

!     Follows the prescription of Kearsley, Acta Cryst. A, 45, 208-210, 1989, making necessary changes 
!     to conform to right-handed rotation in the right-handed coordinate system.
!     Brings RB to the best alignment with RA at the centre of mass of RA
!     Returns DIST as the actual distance, rather than the squared distance 
!     http://dx.doi.org/10.1107/S0108767388010128

      USE COMMONS, ONLY: NRBSITES
      USE KEY, ONLY: NTSITES, DBPT, DBPTDT, DMBLPYT, MSSTOCKT, STOCKAAT, EFIELDT, BULKT, MULTISITEPYT

      IMPLICIT NONE

      INTEGER :: NATOMS
      DOUBLE PRECISION :: RA(3*NATOMS)
      DOUBLE PRECISION :: RB(3*NATOMS)
      DOUBLE PRECISION :: DIST, Q2(4)
      LOGICAL          :: DEBUG
      INTEGER          :: J1, J2, NSIZE, JMIN, INFO
      DOUBLE PRECISION :: QMAT(4,4), TEMPA(9*NATOMS), XM, YM, ZM, XP, YP, ZP
      DOUBLE PRECISION :: DIAG(4), MINV, CMXA, CMYA, CMZA, CMXB, CMYB, CMZB
      DOUBLE PRECISION :: R(3), P(3), RM(3,3)
      DOUBLE PRECISION, ALLOCATABLE :: XA(:), XB(:)
      DOUBLE PRECISION :: ENERGY, VNEW(3*NATOMS), RMS

      IF ((DBPT.AND.EFIELDT).OR.(DBPTDT.AND.EFIELDT).OR.(DMBLPYT.AND.EFIELDT)  &
     &    .OR.(MSSTOCKT.AND.EFIELDT).OR.(STOCKAAT.AND.EFIELDT).OR.(MULTISITEPYT.AND.EFIELDT)) THEN

         CALL FLDMINDIST(RA,RB,NATOMS,DIST,DEBUG,Q2)
         RETURN
    
      ENDIF 

      IF (BULKT) THEN
         CALL RBMINDIST_BULK(RA,RB,NATOMS,DIST,Q2,DEBUG)
         RETURN
      ENDIF

      NSIZE = NTSITES
      ALLOCATE(XA(3*NSIZE),XB(3*NSIZE))

!     MOVE CENTRE OF COORDINATES OF XA AND XB TO THE ORIGIN

      CMXA = 0.0D0; CMYA = 0.0D0; CMZA = 0.0D0
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         CMXA = CMXA + RA(J2+1)
         CMYA = CMYA + RA(J2+2)
         CMZA = CMZA + RA(J2+3)
      ENDDO
      CMXA = 2.D0*CMXA/NATOMS; CMYA = 2.D0*CMYA/NATOMS; CMZA = 2.D0*CMZA/NATOMS
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         RA(J2+1) = RA(J2+1) - CMXA
         RA(J2+2) = RA(J2+2) - CMYA
         RA(J2+3) = RA(J2+3) - CMZA
      ENDDO

      CMXB = 0.0D0; CMYB = 0.0D0; CMZB = 0.0D0
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         CMXB = CMXB + RB(J2+1)
         CMYB = CMYB + RB(J2+2)
         CMZB = CMZB + RB(J2+3)
      ENDDO
      CMXB = 2.D0*CMXB/NATOMS; CMYB = 2.D0*CMYB/NATOMS; CMZB = 2.D0*CMZB/NATOMS
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         RB(J2+1) = RB(J2+1) - CMXB
         RB(J2+2) = RB(J2+2) - CMYB
!sf344 modification         
         RB(J2+3) = RB(J2+3) - CMZB
      ENDDO

!     CONVERT TO THE CARTESIAN CORDINATES FOR THE SITES

      CALL SITEPOS(RA,XA)
      CALL SITEPOS(RB,XB)
 
!     Create matrix QMAT

      QMAT(1:4,1:4) = 0.0D0

      DO J1 = 1, NSIZE
         J2 = 3*(J1-1) 
         XM = XA(J2+1) - XB(J2+1)
         YM = XA(J2+2) - XB(J2+2)
         ZM = XA(J2+3) - XB(J2+3)
         XP = XA(J2+1) + XB(J2+1)
         YP = XA(J2+2) + XB(J2+2)
         ZP = XA(J2+3) + XB(J2+3)
         QMAT(1,1) = QMAT(1,1) + XM**2 + YM**2 + ZM**2
         QMAT(1,2) = QMAT(1,2) - YP*ZM + YM*ZP
         QMAT(1,3) = QMAT(1,3) - XM*ZP + XP*ZM
         QMAT(1,4) = QMAT(1,4) - XP*YM + XM*YP
         QMAT(2,2) = QMAT(2,2) + YP**2 + ZP**2 + XM**2
         QMAT(2,3) = QMAT(2,3) + XM*YM - XP*YP
         QMAT(2,4) = QMAT(2,4) + XM*ZM - XP*ZP
         QMAT(3,3) = QMAT(3,3) + XP**2 + ZP**2 + YM**2
         QMAT(3,4) = QMAT(3,4) + YM*ZM - YP*ZP
         QMAT(4,4) = QMAT(4,4) + XP**2 + YP**2 + ZM**2
      ENDDO

      QMAT(2,1) = QMAT(1,2); QMAT(3,1) = QMAT(1,3); QMAT(3,2) = QMAT(2,3); QMAT(4,1) = QMAT(1,4) 
      QMAT(4,2) = QMAT(2,4); QMAT(4,3) = QMAT(3,4)

!     Diagonalize QMAT
      CALL DSYEV('V','U',4,QMAT,4,DIAG,TEMPA,9*NATOMS,INFO)

      IF (INFO /= 0) PRINT '(A,I6,A)','rbmindist> WARNING - INFO=',INFO,' in DSYEV'

!     Calculate the distance DIST from DIAG (from QMAT)
!     MINV = the lowest eigenvalue ?
!     JMIN = the index of the lowest eigenvalue ?
      MINV = 1.0D100
      DO J1 = 1,4
         IF (DIAG(J1).LT.MINV) THEN
            JMIN = J1
            MINV = DIAG(J1)
         ENDIF
      ENDDO
      IF (MINV < 0.0D0) THEN
         IF (ABS(MINV)< 1.0D-6) THEN
            MINV = 0.0D0
         ELSE
            PRINT '(A,G20.10,A)','newmindist> WARNING MINV is ',MINV,' change to absolute value'
            MINV = -MINV
         ENDIF
      ENDIF
      DIST = SQRT(MINV)

!     Calculate Q2 from QMAT and JMIN
      Q2(1) = QMAT(1,JMIN); Q2(2) = QMAT(2,JMIN); Q2(3) = QMAT(3,JMIN); Q2(4) = QMAT(4,JMIN)

!     Update RB based on center of mass
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         RB(J2+1) = RB(J2+1) - CMXB
         RB(J2+2) = RB(J2+2) - CMYB
! sf344 modification         
         RB(J2+3) = RB(J2+3) - CMZB
      ENDDO

!     rotate RB based on Q2 and center of mass
      CALL RBNEWROTGEOM(NATOMS,RB,Q2,RM,CMXA,CMYA,CMZA)

      !write(*,*) "RBMINDIST> after  ", RA(1), RA(2), RA(3), RB(1), RB(2), RB(3)
      DEALLOCATE(XA,XB)

      END SUBROUTINE RBMINDIST

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBNEWROTGEOM(NATOMS,COORDS,Q2,RM,CX,CY,CZ)

      IMPLICIT NONE

      INTEGER, INTENT(IN) :: NATOMS
      DOUBLE PRECISION, INTENT(INOUT) :: COORDS(3*NATOMS)
      DOUBLE PRECISION, INTENT(IN) :: CX, CY, CZ, Q2(4)
      DOUBLE PRECISION, INTENT(OUT) :: RM(3,3)
      INTEGER          :: I, J
      DOUBLE PRECISION :: R(3), P(3)

!     RMAT CONTAINS THE MATRIX THAT MAPS RB ONTO THE BEST CORRESPONDENCE WITH RA

      CALL QROTMAT(Q2,RM)

      DO I = 1, NATOMS/2
  
         J    = 3*(I-1)
         R(:) = MATMUL(RM(:,:), COORDS(J+1:J+3))

         COORDS(J+1) = R(1) + CX
         COORDS(J+2) = R(2) + CY
         COORDS(J+3) = R(3) !+ CZ
      
!     CONVERT THE ANGLE-AXIS COORDINATES

         J      = 3*NATOMS/2 + J
         P(:)   = COORDS(J+1:J+3)

         CALL QROTAA(Q2,P)

         COORDS(J+1:J+3) = P(1:3)

      ENDDO

      END SUBROUTINE RBNEWROTGEOM

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE FLDMINDIST(RA,RB,NATOMS,DIST,DEBUG,Q2)

!     returns DIST as the actual distance, rather than the squared distance

      USE COMMONS, ONLY: NRBSITES
      USE KEY, ONLY: NTSITES, STOCKAAT, MULTISITEPYT

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J4, NATOMS, NSIZE, JMIN, INFO
      DOUBLE PRECISION :: RA(3*NTSITES), RB(3*NTSITES), DIST, QMAT(2,2), XM, YM, ZM, XP, YP, ZP
      DOUBLE PRECISION :: MINV, Q2(4), CMXA, CMYA, CMZA, CMXB, CMYB, CMZB
      DOUBLE PRECISION :: R(3), P(3), RM(3,3) 
      DOUBLE PRECISION, ALLOCATABLE :: XA(:), XB(:)
      DOUBLE PRECISION :: ENERGY, VNEW(3*NATOMS), RMS, DUMMY
      LOGICAL          :: DEBUG

      NSIZE = NTSITES
      ALLOCATE(XA(3*NSIZE),XB(3*NSIZE))

!     MOVE CENTRE OF COORDINATES TO THE ORIGIN

      CMXA = 0.0D0; CMYA = 0.0D0; CMZA = 0.0D0
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         CMXA = CMXA + RA(J2+1)
         CMYA = CMYA + RA(J2+2)
         CMZA = CMZA + RA(J2+3)
      ENDDO
      CMXA = 2.D0*CMXA/NATOMS; CMYA = 2.D0*CMYA/NATOMS; CMZA = 2.D0*CMZA/NATOMS
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         RA(J2+1) = RA(J2+1) - CMXA
         RA(J2+2) = RA(J2+2) - CMYA
! do not move along the z axis in case of a gravity field (currently implemented for MULTISITEPY)
         IF(MULTISITEPYT) CMZA=0.0D0 
         RA(J2+3) = RA(J2+3) - CMZA 
      ENDDO

      CMXB = 0.0D0; CMYB = 0.0D0; CMZB = 0.0D0
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         CMXB = CMXB + RB(J2+1)
         CMYB = CMYB + RB(J2+2)
         CMZB = CMZB + RB(J2+3)
      ENDDO
      CMXB = 2.D0*CMXB/NATOMS; CMYB = 2.D0*CMYB/NATOMS; CMZB = 2.D0*CMZB/NATOMS
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         RB(J2+1) = RB(J2+1) - CMXB
         RB(J2+2) = RB(J2+2) - CMYB
! do not move along the z axis in case of a gravity field (currently implemented for MULTISITEPY)
         IF(MULTISITEPYT) CMZB=0.0D0
         RB(J2+3) = RB(J2+3) - CMZB
      ENDDO

!     CONVERT TO THE CARTESIAN CORDINATES FOR THE SITES

      IF (STOCKAAT) THEN
         XA(:) = RA(:)
         XB(:) = RB(:)
      ELSE
         CALL SITEPOS(RA,XA)
         CALL SITEPOS(RB,XB)
      ENDIF

      QMAT(1:2,1:2) = 0.0D0

      DO J1 = 1, NSIZE
         J2 = 3*(J1-1) 
         XM = XA(J2+1) - XB(J2+1)
         YM = XA(J2+2) - XB(J2+2)
         ZM = XA(J2+3) - XB(J2+3)
         XP = XA(J2+1) + XB(J2+1)
         YP = XA(J2+2) + XB(J2+2)
         ZP = XA(J2+3) + XB(J2+3)
         QMAT(1,1) = QMAT(1,1) + XM**2 + YM**2 + ZM**2
         QMAT(1,2) = QMAT(1,2) - XP*YM + XM*YP
         QMAT(2,2) = QMAT(2,2) + XP**2 + YP**2 + ZM**2
      ENDDO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      MINV=MIN(QMAT(1,1),QMAT(2,2))
      Q2(1)=1.D0
      Q2(2:4)=0.D0
      IF (QMAT(1,2)**2>1.D-8) THEN
         MINV = 0.5D0*(QMAT(1,1) + QMAT(2,2) - SQRT(4.D0*QMAT(1,2)*QMAT(1,2) + (QMAT(1,1) - QMAT(2,2))**2.D0))
         Q2(1) = SQRT((MINV-QMAT(2,2))**2.D0/(QMAT(1,2)*QMAT(1,2) + (MINV-QMAT(2,2))**2.D0))
         Q2(4) = QMAT(1,2)*Q2(1)/(MINV - QMAT(2,2))
      ENDIF
      IF (DEBUG) WRITE(*,'(3F20.10)')Q2(1),Q2(4),ABS(Q2(1)**2+Q2(4)**2-1.D0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!     QMAT IS SYMMETRIC; QMAT(2,1) = QMAT(1,2)
!
!      MINV = 0.5D0*(QMAT(1,1) + QMAT(2,2) - SQRT(4.D0*QMAT(1,2)*QMAT(1,2) + (QMAT(1,1) - QMAT(2,2))**2.D0))
!
!      Q2(2) = 0.D0; Q2(3) = 0.D0
!
!      IF ((QMAT(1,1)<1.D-12).AND.(QMAT(1,2)<1.D-12).AND.(QMAT(2,2)<1.D-12)) THEN
!         Q2(1) = 1.D0
!         Q2(4) = 0.D0
!      ELSEIF (((MINV-QMAT(1,1))**2.D0 + QMAT(1,2)*QMAT(1,2)) < 1.D-12 .OR. QMAT(1,2) < 1.D-12) THEN
!         Q2(1) = SQRT((MINV-QMAT(2,2))**2.D0/(QMAT(1,2)*QMAT(1,2) + (MINV-QMAT(2,2))**2.D0))
!         Q2(4) = QMAT(1,2)*Q2(1)/(MINV - QMAT(2,2))
!      ELSE
!         Q2(1) = SQRT(QMAT(1,2)*QMAT(1,2)/((MINV-QMAT(1,1))**2.D0 + QMAT(1,2)*QMAT(1,2)))
!         Q2(4) = (MINV - QMAT(1,1))*Q2(1)/QMAT(1,2)
!      ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      IF (MINV < 0.0D0) THEN
         IF (ABS(MINV)< 1.0D-6) THEN
            MINV = 0.0D0
         ELSE
            PRINT '(A,G20.10,A)','newmindist> WARNING MINV is ',MINV,' change to absolute value'
            MINV = -MINV
         ENDIF
      ENDIF
      DIST = SQRT(MINV)
  
      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         RB(J2+1) = RB(J2+1) - CMXB
         RB(J2+2) = RB(J2+2) - CMYB
! do not move along the z axis in case of a gravity field (currently implemented for MULTISITEPY)
         IF(MULTISITEPYT) CMZB=0.0D0 
         RB(J2+3) = RB(J2+3) - CMZB
      ENDDO

      CALL RBNEWROTGEOM(NATOMS,RB,Q2,RM,CMXA,CMYA,CMZA)

      DEALLOCATE(XA,XB)

      END SUBROUTINE FLDMINDIST

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBNEWROTGEOMMYORIENT(NATOMS,COORDS,RM,CX,CY,CZ)

      IMPLICIT NONE

      INTEGER          :: I, J, NATOMS
      DOUBLE PRECISION :: COORDS(3*NATOMS), RM(3,3), CX, CY, CZ, R(3), P(3), Q1(4), Q2(4), Q(4)
      DOUBLE PRECISION :: THETA, THETAH, ST, FCT

      DOUBLE PRECISION :: r1, diag1, diag2, diag3
      INTEGER          :: u,v,w
!     RMAT CONTAINS THE MATRIX THAT MAPS RB ONTO THE BEST CORRESPONDENCE WITH RA
!     sf344> when called from minpermdist, RM is already given from putting the centres of rigid bodies 
!            into standard orientation, so this subroutine simply rotates all coordinates with that rotation
!            matrix

!      RM(1,1) = Q2(1)**2 + Q2(2)**2 - Q2(3)**2 - Q2(4)**2
!      RM(1,2) = 2.D0*(Q2(2)*Q2(3) - Q2(1)*Q2(4))
!      RM(1,3) = 2.D0*(Q2(2)*Q2(4) + Q2(1)*Q2(3))
!      RM(2,1) = 2.D0*(Q2(2)*Q2(3) + Q2(1)*Q2(4))
!      RM(2,2) = Q2(1)**2 + Q2(3)**2 -Q2(2)**2-Q2(4)**2
!      RM(2,3) = 2.D0*(Q2(3)*Q2(4) - Q2(1)*Q2(2))
!      RM(3,1) = 2.D0*(Q2(2)*Q2(4) - Q2(1)*Q2(3))
!      RM(3,2) = 2.D0*(Q2(3)*Q2(4) + Q2(1)*Q2(2))
!      RM(3,3) = Q2(1)**2 +Q2(4)**2 -Q2(2)**2 - Q2(3)**2

      diag1=RM(1,1)
      diag2=RM(2,2)
      diag3=RM(3,3)

       
!     if the rotation matrix is the identity matrix, then return
      IF(abs(diag1-1.0D0)<1.0D-8.AND.abs(diag2-1.0D0)<1.0D-8.AND.abs(diag3-1.0D0)<1.0D-8) THEN
        RETURN
      END IF
!     otherwise figure out the quaternion from the rotation matrix
!      WRITE(*,'(A)') 'coords before rotation'
!      WRITE(*,'(3F13.8)') COORDS(:)

      IF(ABS(RM(1,1))>=ABS(RM(2,2)).AND.ABS(RM(1,1))>=ABS(RM(3,3))) THEN
          diag1=RM(1,1)
          u=1
          v=2
          w=3
      ELSE IF (ABS(RM(2,2))>=ABS(RM(1,1)).AND.ABS(RM(2,2))>=ABS(RM(3,3))) THEN
          diag1=RM(2,2)
          diag2=RM(3,3)
          diag3=RM(1,1)
          u=2
          v=3
          w=1
      ELSE IF (ABS(RM(3,3))>=ABS(RM(1,1)).AND.ABS(RM(3,3))>=ABS(RM(2,2))) THEN
          diag1=RM(3,3)
          diag2=RM(1,1)
          diag3=RM(2,2)
          u=3
          v=1
          w=2
      END IF

      r1=SQRT(1+RM(u,u)-RM(v,v)-RM(w,w))
      
      Q2(1)=(RM(w,v)-RM(v,w))/(2.0D0*r1)
      Q2(u+1)=r1/2.0D0
      Q2(v+1)=(RM(u,v)+RM(v,u))/(2.0D0*r1)
      Q2(w+1)=(RM(w,u)+RM(u,w))/(2.0D0*r1)

      DO I = 1, NATOMS/2
  
         J    = 3*(I-1)
         R(:) = MATMUL(RM(:,:), COORDS(J+1:J+3))

         COORDS(J+1) = R(1) + CX
         COORDS(J+2) = R(2) + CY
         COORDS(J+3) = R(3) + CZ
      
!     CONVERT THE ANGLE-AXIS COORDINATES

         J      = 3*NATOMS/2 + J
         P(:)   = COORDS(J+1:J+3)
         THETA  = DSQRT(DOT_PRODUCT(P,P))
         THETAH = 0.5D0*THETA
         ST     = SIN(THETAH)
!        WRITE(*,*) 'st=', THETAH, ST, P(:)
         Q1(1)  = COS(THETAH)
         Q1(2)  = P(1)*ST/THETA
         Q1(3)  = P(2)*ST/THETA
         Q1(4)  = P(3)*ST/THETA

         Q(1)   = Q2(1)*Q1(1) - Q2(2)*Q1(2) - Q2(3)*Q1(3) - Q2(4)*Q1(4)
         Q(2)   = Q2(1)*Q1(2) + Q2(2)*Q1(1) + Q2(3)*Q1(4) - Q2(4)*Q1(3)
         Q(3)   = Q2(1)*Q1(3) + Q2(3)*Q1(1) + Q2(4)*Q1(2) - Q2(2)*Q1(4)
         Q(4)   = Q2(1)*Q1(4) + Q2(4)*Q1(1) + Q2(2)*Q1(3) - Q2(3)*Q1(2)

         THETA  = 2.D0*ACOS(Q(1))

         IF (THETA == 0.D0) THEN
            COORDS (J+1:J+3) = 0.D0
         ELSE
            FCT    = DSQRT(DOT_PRODUCT(Q(2:4),Q(2:4)))
            COORDS(J+1) = THETA*Q(2)/FCT 
            COORDS(J+2) = THETA*Q(3)/FCT
            COORDS(J+3) = THETA*Q(4)/FCT
         ENDIF

      ENDDO
!      WRITE(*,'(A)') 'coords after rotation'
!      WRITE(*,'(3F13.8)') COORDS(:)

      END SUBROUTINE RBNEWROTGEOMMYORIENT

!     ----------------------------------------------------------------------------------------------


SUBROUTINE UNIAXGETPATHLENGTH(RA,RB,TEMP)
use key, only : PYA1
use commons, only : NATOMS
implicit none

integer          :: NSIZE, J1, J2, J3, J4, NRBSITES
double precision :: XA(9*NATOMS/2), XB(9*NATOMS/2), RA(3*NATOMS), RB(3*NATOMS), P(3), RM(3,3), RBSITE(3,3), TEMP
double precision :: R(3),rbdistance

      NRBSITES = 3
      NSIZE = NRBSITES*NATOMS/2
! sf344> try to generalize this for ellipsoids, using sites displaced by 1 along the symmetry axis      
      IF(PYA1(3)==0.0D0) THEN
        rbdistance= 1.0D0
      ELSE
        rbdistance=PYA1(3)
      END IF

      RBSITE(1,1) = 0.0D0
      RBSITE(1,2) = 0.0D0
      RBSITE(1,3) = 0.0D0
      RBSITE(2,1) = 0.0D0
      RBSITE(2,2) = 0.0D0
      RBSITE(2,3) = rbdistance
      RBSITE(3,1) = 0.0D0
      RBSITE(3,2) = 0.0D0
      RBSITE(3,3) = -rbdistance
!      ALLOCATE(XA(3*NSIZE),XB(3*NSIZE))

!      write(*,*) 'ra'
!      write(*,*) ra(:)
!     ----------------------------------------------------------------------------------------------
!      CALL LWOTPGH(RB,VNEW,ENERGY,.TRUE.,.FALSE.)
!      WRITE(*,*) ENERGY, SQRT(DOT_PRODUCT(VNEW(:),VNEW(:)))
!      CALL DUMBBELLP(RB,VNEW,ENERGY,.TRUE.,.FALSE.)
!      WRITE(*,*) ENERGY, SQRT(DOT_PRODUCT(VNEW(:),VNEW(:)))
!     ----------------------------------------------------------------------------------------------
!     CONVERT TO THE CARTESIAN CORDINATES FOR THE SITES

      DO J1 = 1, NATOMS/2

         J2   = 3*J1
         R(:) = RA(J2-2:J2)
         J2   = 3*NATOMS/2 + J2
         P(:) = RA(J2-2:J2)

         CALL ROTMAT(P, RM)
 
         DO J3 = 1, NRBSITES

            J4          = 3*((J1-1)*NRBSITES + J3)
            XA(J4-2:J4) = R(:) + MATMUL(RM(:,:), RBSITE(J3,:))

         ENDDO

         J2   = 3*J1
         R(:) = RB(J2-2:J2)
         J2   = 3*NATOMS/2 + J2
         P(:) = RB(J2-2:J2)
 
         CALL ROTMAT(P, RM)

         DO J3 = 1, NRBSITES

            J4          = 3*((J1-1)*NRBSITES + J3)
            XB(J4-2:J4) = R(:) + MATMUL(RM(:,:), RBSITE(J3,:))

         ENDDO

      ENDDO

      TEMP=0.0D0
      DO J2=1,NSIZE
            TEMP=TEMP+(XB(J2)-XA(J2))**2
      ENDDO

      END SUBROUTINE UNIAXGETPATHLENGTH

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBCOMMINDIST(RA,RB,NATOMS,DIST,RM,DEBUG)

!     returns squared distance DIST

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3, J4, NATOMS, NSIZE, JMIN, INFO
      DOUBLE PRECISION :: RA(3*NATOMS), RB(3*NATOMS), DIST, QMAT(4,4), TEMPA(9*NATOMS), XM, YM, ZM, XP, YP, ZP
      DOUBLE PRECISION :: DIAG(4), MINV, Q2(4), CMXA, CMYA, CMZA, CMXB, CMYB, CMZB
      DOUBLE PRECISION :: R(3), P(3), RM(3,3) 
      DOUBLE PRECISION, ALLOCATABLE :: XA(:), XB(:)
      DOUBLE PRECISION :: ENERGY, VNEW(3*NATOMS), RMS, DUMMY
      LOGICAL          :: BULKT, PRESERVET, DEBUG 

      NSIZE = NATOMS/2
      ALLOCATE(XA(3*NSIZE),XB(3*NSIZE))
      XA(1:3*NSIZE) = RA(1:3*NSIZE); XB(1:3*NSIZE) = RB(1:3*NSIZE)

!     MOVE CENTRE OF COORDINATES OF XA AND XB TO THE ORIGIN

      CMXA = 0.0D0; CMYA = 0.0D0; CMZA = 0.0D0
      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         CMXA = CMXA + XA(J2+1)
         CMYA = CMYA + XA(J2+2)
         CMZA = CMZA + XA(J2+3)
      ENDDO
      CMXA = CMXA/NSIZE; CMYA = CMYA/NSIZE; CMZA = CMZA/NSIZE
      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         XA(J2+1) = XA(J2+1) - CMXA
         XA(J2+2) = XA(J2+2) - CMYA
         XA(J2+3) = XA(J2+3) - CMZA
      ENDDO

      CMXB = 0.0D0; CMYB = 0.0D0; CMZB = 0.0D0
      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         CMXB = CMXB + XB(J2+1)
         CMYB = CMYB + XB(J2+2)
         CMZB = CMZB + XB(J2+3)
      ENDDO
      CMXB = CMXB/NSIZE; CMYB = CMYB/NSIZE; CMZB = CMZB/NSIZE
      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         XB(J2+1) = XB(J2+1) - CMXB
         XB(J2+2) = XB(J2+2) - CMYB
         XB(J2+3) = XB(J2+3) - CMZB
      ENDDO

      QMAT(1:4,1:4) = 0.0D0

      DO J1 = 1, NSIZE
         J2 = 3*(J1-1)
         XM = XA(J2+1) - XB(J2+1)
         YM = XA(J2+2) - XB(J2+2)
         ZM = XA(J2+3) - XB(J2+3)
         XP = XA(J2+1) + XB(J2+1)
         YP = XA(J2+2) + XB(J2+2)
         ZP = XA(J2+3) + XB(J2+3)
         QMAT(1,1) = QMAT(1,1) + XM**2 + YM**2 + ZM**2
         QMAT(1,2) = QMAT(1,2) - YP*ZM + YM*ZP 
         QMAT(1,3) = QMAT(1,3) - XM*ZP + XP*ZM
         QMAT(1,4) = QMAT(1,4) - XP*YM + XM*YP
         QMAT(2,2) = QMAT(2,2) + YP**2 + ZP**2 + XM**2
         QMAT(2,3) = QMAT(2,3) + XM*YM - XP*YP 
         QMAT(2,4) = QMAT(2,4) + XM*ZM - XP*ZP
         QMAT(3,3) = QMAT(3,3) + XP**2 + ZP**2 + YM**2
         QMAT(3,4) = QMAT(3,4) + YM*ZM - YP*ZP 
         QMAT(4,4) = QMAT(4,4) + XP**2 + YP**2 + ZM**2
      ENDDO

      QMAT(2,1) = QMAT(1,2); QMAT(3,1) = QMAT(1,3); QMAT(3,2) = QMAT(2,3); QMAT(4,1) = QMAT(1,4)
      QMAT(4,2) = QMAT(2,4); QMAT(4,3) = QMAT(3,4)
      CALL DSYEV('V','U',4,QMAT,4,DIAG,TEMPA,9*NATOMS,INFO)

      IF (INFO /= 0) PRINT '(A,I6,A)','newmindist> WARNING - INFO=',INFO,' in DSYEV'

      MINV = 1.0D100
      DO J1 = 1,4
         IF (DIAG(J1).LT.MINV) THEN
            JMIN = J1
            MINV = DIAG(J1)
         ENDIF
      ENDDO
      IF (MINV < 0.0D0) THEN
         IF (ABS(MINV)< 1.0D-6) THEN
            MINV = 0.0D0
         ELSE
            PRINT '(A,G20.10,A)','newmindist> WARNING MINV is ',MINV,' change to absolute value'
            MINV = -MINV
         ENDIF
      ENDIF
      DIST = MINV

      IF (DEBUG) PRINT '(A,G20.10,A,I6)',' rbmindist2> minimum residual is ',DIAG(JMIN),' for eigenvector ',JMIN

      Q2(1) = QMAT(1,JMIN); Q2(2) = QMAT(2,JMIN); Q2(3) = QMAT(3,JMIN); Q2(4) = QMAT(4,JMIN)

      DO J1 = 1, NATOMS/2
         J2 = 3*(J1-1)
         RB(J2+1) = RB(J2+1) - CMXB
         RB(J2+2) = RB(J2+2) - CMYB
         RB(J2+3) = RB(J2+3) - CMZB
      ENDDO

      CALL RBNEWROTGEOM(NATOMS,RB,Q2,RM,CMXA,CMYA,CMZA)

      DEALLOCATE(XA,XB)

      END SUBROUTINE RBCOMMINDIST

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE CHECKDRVTS(X)

      USE COMMONS
      USE KEY
      USE MODHESS

      IMPLICIT NONE

      INTEGER          :: IVRNO1, IVRNO2
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), ENERGY, RMS, FM, FP, DF, DFN
      LOGICAL          :: GTEST, STEST
! jwrm2> Added this to make changing the error limit and coordinate offset easier
      DOUBLE PRECISION, PARAMETER :: ERRLIM = 1.D-06, DELX = 1.D-06

      IF (CHECKDID == 0) THEN 
!     ENERGY
         CALL POTENTIAL(X, ENERGY, G, .FALSE., .FALSE., RMS, .FALSE., .FALSE.)
         WRITE(*,*) 'ENERGY =', ENERGY
      ELSEIF (CHECKDID == 1) THEN
!     GRADIENTS
         DO IVRNO1 = 1, 3*NATOMS

! jwrm2> Skip z coordinates in 2D systems
            IF (TWOD .AND. (MODULO(IVRNO1, 3) .EQ. 0)) CYCLE 

            WRITE(*, *) IVRNO1
            X(IVRNO1) = X(IVRNO1) - DELX
            CALL POTENTIAL(X, ENERGY, G, .FALSE., .FALSE., RMS, .FALSE., .FALSE.)
            FM   = ENERGY
!            WRITE(*, *) FM
            X(IVRNO1) = X(IVRNO1) + 2.D0*DELX
            CALL POTENTIAL(X, ENERGY, G, .FALSE., .FALSE., RMS, .FALSE., .FALSE.)
            FP   = ENERGY
!            WRITE(*, *) FP
            X(IVRNO1) = X(IVRNO1) - DELX
            CALL POTENTIAL(X, ENERGY, G, .TRUE., .FALSE., RMS, .FALSE., .FALSE.)
            DFN = (FP - FM) / (2.D0*DELX)
            DF   = G(IVRNO1)

            WRITE(*, *) 'DFN=', DFN
            WRITE(*, *) 'DFA=', DF

            IF (ABS(DFN - DF) > ERRLIM) WRITE(*, *) IVRNO1, DFN, DF
         ENDDO
      ELSEIF (CHECKDID == 2) THEN
!     HESSIAN
         ALLOCATE(HESS(3*NATOMS,3*NATOMS))
         DO IVRNO1 = 1, 3*NATOMS

! jwrm2> Skip z coordinates in 2D systems
            IF (TWOD .AND. (MODULO(IVRNO1, 3) .EQ. 0)) CYCLE

            DO IVRNO2 = 1, 3*NATOMS

! jwrm2> Skip z coordinates in 2D systems
               IF (TWOD .AND. (MODULO(IVRNO2, 3) .EQ. 0)) CYCLE

               WRITE(*,*) IVRNO1, IVRNO2
               X(IVRNO1) = X(IVRNO1) - DELX
               CALL POTENTIAL (X, ENERGY, G, .TRUE., .FALSE., RMS, .FALSE., .FALSE.)
               FM   = G(IVRNO2)
!               WRITE(*, *) FM
               X(IVRNO1) = X(IVRNO1) + 2.D0*DELX
               CALL POTENTIAL (X, ENERGY, G, .TRUE., .FALSE., RMS, .FALSE., .FALSE.)
               FP   = G(IVRNO2)
!               WRITE(*, *) FP
               X(IVRNO1) = X(IVRNO1) - DELX
               CALL POTENTIAL (X, ENERGY, G, .TRUE., .TRUE., RMS, .FALSE., .FALSE.)
               DFN  = (FP - FM) / (2.D0*DELX)
               DF   = HESS(IVRNO1,IVRNO2)
               WRITE(*, *) 'DFN=', DFN
               WRITE(*, *) 'DFA=', DF
               IF (ABS(DFN - DF) > ERRLIM) WRITE(*,*) IVRNO1, IVRNO2, DFN, DF, ABS(DFN - DF)

            ENDDO

         ENDDO

      ENDIF

      STOP

      END SUBROUTINE CHECKDRVTS

!     ---------------------------------------------------------------------------------------------

      SUBROUTINE RMDFASN(P, RM, DRM1, DRM2, DRM3, D2RM1, D2RM2, D2RM3, D2RI12, D2RI23, D2RI31, GTEST, STEST)

      IMPLICIT NONE

      DOUBLE PRECISION :: P(3), PN(3), THETA, THETA2, THETA3, CT, ST, E(3,3), ESQ(3,3), I3(3,3) 
      DOUBLE PRECISION :: DE1(3,3), DE2(3,3), DE3(3,3), RM(3,3), DRM1(3,3), DRM2(3,3), DRM3(3,3)
      DOUBLE PRECISION :: D2E1(3,3), D2E2(3,3), D2E3(3,3), D2E12(3,3), D2E23(3,3), D2E31(3,3)
      DOUBLE PRECISION :: D2RM1(3,3), D2RM2(3,3), D2RM3(3,3)
      DOUBLE PRECISION :: D2RI12(3,3), D2RI23(3,3), D2RI31(3,3)
      DOUBLE PRECISION :: FCTR, FCTRSQ1, FCTRSQ2, FCTRSQ3
      DOUBLE PRECISION :: COSOMG, SINOMG, COSTHT, SINTHT, COSPSI, SINPSI
      LOGICAL          :: GTEST, STEST

      I3(:,:) = 0.D0
      I3(1,1) = 1.D0; I3(2,2) = 1.D0; I3(3,3) = 1.D0
!     P(1) = \THETA, P(2) = \PSI, P(3) = \OMEGA

      COSTHT  = COS(P(1))
      SINTHT  = SIN(P(1))
      COSPSI  = COS(P(2))
      SINPSI  = SIN(P(2))
      COSOMG  = COS(P(3))
      SINOMG  = SIN(P(3))
      PN(1)   = SINTHT*COSPSI
      PN(2)   = SINTHT*SINPSI
      PN(3)   = COSTHT
      E(:,:)  = 0.D0
      E(1,2)  = -PN(3)
      E(1,3)  =  PN(2)
      E(2,3)  = -PN(1)
      E(2,1)  = -E(1,2)
      E(3,1)  = -E(1,3)
      E(3,2)  = -E(2,3)

      ESQ(:,:) = MATMUL(E(:,:),E(:,:))
      RM      = I3(:,:) + (1.D0-COSOMG)*ESQ(:,:) + SINOMG*E(:,:)

      IF (.NOT. GTEST .AND. .NOT. STEST) RETURN

      DE1(:,:) = 0.D0
      DE1(1,2) = SINTHT
      DE1(1,3) = COSTHT*SINPSI
      DE1(2,3) =-COSTHT*COSPSI
      DE1(2,1) =-DE1(1,2)
      DE1(3,1) =-DE1(1,3)
      DE1(3,2) =-DE1(2,3)

      DE2(:,:) = 0.D0
      DE2(1,3) = PN(1)
      DE2(2,3) = PN(2)
      DE2(3,1) =-DE2(1,3)
      DE2(3,2) =-DE2(2,3)

      DRM1(:,:) = (1.D0-COSOMG)*(MATMUL(DE1,E) + MATMUL(E,DE1)) + SINOMG*DE1(:,:)
      DRM2(:,:) = (1.D0-COSOMG)*(MATMUL(DE2,E) + MATMUL(E,DE2)) + SINOMG*DE2(:,:)
      DRM3(:,:) = SINOMG*ESQ(:,:) + COSOMG*E(:,:)

      IF (.NOT. STEST) RETURN


      D2E1(:,:) = 0.D0
      D2E1(1,2) = COSTHT
      D2E1(1,3) =-PN(2)
      D2E1(2,3) = PN(1)
      D2E1(2,1) =-D2E1(1,2)
      D2E1(3,1) =-D2E1(1,3)
      D2E1(3,2) =-D2E1(2,3)

      D2E2(:,:) = 0.D0
      D2E2(1,3) =-SINTHT*SINPSI
      D2E2(2,3) = SINTHT*COSPSI
      D2E2(3,1) =-D2E2(1,3)
      D2E2(3,2) =-D2E2(2,3)

      D2E12(:,:) = 0.D0
      D2E12(1,3) = COSTHT*COSPSI
      D2E12(2,3) = COSTHT*SINPSI
      D2E12(3,1) =-D2E12(1,3)
      D2E12(3,2) =-D2E12(2,3)

      D2RM1(:,:)  = (1.D0-COSOMG)*(MATMUL(D2E1,E)+2.D0*MATMUL(DE1,DE1)+MATMUL(E,D2E1))+SINOMG*D2E1(:,:)

      D2RM2(:,:)  = (1.D0-COSOMG)*(MATMUL(D2E2,E)+2.D0*MATMUL(DE2,DE2)+MATMUL(E,D2E2))+SINOMG*D2E2(:,:)

      D2RM3(:,:)  = COSOMG*ESQ(:,:) - SINOMG*E(:,:)

      D2RI12(:,:) = (1.D0-COSOMG)*(MATMUL(D2E12,E)+MATMUL(DE1,DE2)+MATMUL(DE2,DE1)+MATMUL(E,D2E12)) &
                  & +SINOMG*D2E12(:,:)

      D2RI23(:,:) = SINOMG*(MATMUL(DE2(:,:),E(:,:))+MATMUL(E(:,:),DE2(:,:))) + COSOMG*DE2(:,:)

      D2RI31(:,:) = SINOMG*(MATMUL(DE1(:,:),E(:,:))+MATMUL(E(:,:),DE1(:,:))) + COSOMG*DE1(:,:)

      END SUBROUTINE RMDFASN

!     ----------------------------------------------------------------------------------------------

! ---------------------------------------------------------
! hk286 - computing normal modes for rigid body angle-axis
! ---------------------------------------------------------

      SUBROUTINE NRMLMD (X, FRQN, EIGENVECTORT)

      USE COMMONS
      USE MODHESS

      IMPLICIT NONE

      INTEGER          :: I, J, K, J1, J2, J3, J5, OFFSET, NDIM, IR, IC, K1, K2 
      DOUBLE PRECISION :: X(3*NATOMS), FRQN(3*NATOMS), G(3*NATOMS)
      DOUBLE PRECISION :: TMASS, FRQCNV, ENERGY, RMS
      DOUBLE PRECISION :: KBLOCK(3,3), KBEGNV(3)
      DOUBLE PRECISION :: P(3), RMI(3,3), DRMI(3,3), DR(3)
      DOUBLE PRECISION :: KD(3*NATOMS), U(3*NATOMS,3*NATOMS)
      DOUBLE PRECISION :: HUK(3*NATOMS,3*NATOMS), AP(3*NATOMS,3*NATOMS)
      LOGICAL          :: GTEST, STEST
!     the following required to call the LAPACK routine DSYEV
      INTEGER          :: INFO
      INTEGER, PARAMETER :: LWORK = 10000 ! the dimension is set arbitrarily
      DOUBLE PRECISION :: WORK(LWORK)
! hk286
      LOGICAL          :: EIGENVECTORT

!     Computes the normal modes and frequencies
!     following Pohorille et al. JCP 87, 6070 (1987)
!     INPUT:
!     ndim   : number of degrees of freedom
!     second : second derivatives (ORIENT convention) from derivs
!     common "sites" : coordinates of sites in global coordinates
!     OUTPUT:
!     freq   : eigenfrequencies of normal modes

!     We adopt Pohorille's notation for clarity:
!     K matrix : block diagonal kinetic energy matrix (6N x 6N)
!     kblock : for each rigid body we have one 6x6 matrix which consists
!     of two blocks : left upper or "mass" diagonal submatrix
!     and right lower or inertia tensor matrix. Here kblock
!     is the 3 x 3 inertia tensor.
!     KD     : using the diagonalized kinetic energy tensor
!     we keep track of the diagonal of KD only
!     U      : eigenvector matrix (Pohorille's S)

!     This has now been deprecated. Use the FRQCONV keyword to set your unit conversions. For the TIP4P potential,
!     the value given here is set as the default value of FRQCONV, so frequencies will be given in cm^-1.
!     Frequency conversion factor: Energy in KJ/mol and length in angstrom
!     to get frequencies in cm^{-1}      
!      FRQCNV = 1.D03/(2.D0*4.D0*DATAN(1.D0)*2.998D0)

!     Initialize
      U(:,:) = 0.D0
      IR     = 0
      IC     = 0

!     Get the site positions
      OFFSET = 3*NATOMS/2
      GTEST = .FALSE.; STEST = .FALSE.

      DO J1 = 1, NATOMS/2

         J3 = 3*J1
         J5 = OFFSET + J3
         P  = X(J5-2:J5)
         KBLOCK(:,:) = 0.D0
         CALL RMDFAS(P, RMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, DRMI, GTEST, STEST)
         CALL COMPUTEINERTIA(RMI, KBLOCK, TMASS) 
!     Diagonalise KBLOCK using LAPACK routine DSYEV
!     DSYEV computes all eigenvalues and, optionally, eigenvectors of a real symmetric matrix KBLOCK 
         CALL DSYEV('V','L',3,KBLOCK,3,KBEGNV,WORK,LWORK,INFO)
!     The character 'V' instructs to return the eigenvector as well
!     On exit, if INFO = 0, KBLOCK contains the orthonormal eigenvectors of the matrix KBLOCK in columns
!     The character 'L' tells that the Lower triangle of KBLOCK is stored 
!     Next is the order of the matrix KBLOCK
!     The integer after KBLOCK is the leading dimension of the array KBLOCK
!     KEGNV holds the eigenvalues in ascending order if INFO = 0

         IF (INFO /= 0) THEN
            WRITE(*,*) 'NRMLMD > Error in DSYEV with KBLOCK, INFO =', INFO
            STOP
         ENDIF

!     Construction of the matrix U
!     First: translation coordinates
         U(IR+1,IC+1) = 1.D0; U(IR+2,IC+2) = 1.D0; U(IR+3,IC+3) = 1.D0
         KD(IC+1:IC+3) = 1.D0/SQRT(TMASS)            
!     Now rotational coordinates
         U(OFFSET+IR+1:OFFSET+IR+3,OFFSET+IC+1:OFFSET+IC+3) = KBLOCK(:,:) 
         KD(OFFSET+IC+1:OFFSET+IC+3) = 1.D0/SQRT(KBEGNV(:))
         IR = IR + 3
         IC = IC + 3 

      ENDDO

      NDIM = 3*NATOMS

! hk286 - you want to compute HESSIAN differently for the normal modes

      CALL POTENTIAL(X,ENERGY,G,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)

      AP(:,:) = 0.D0
      DO I = 1, NDIM
         DO J = 1, I
            DO K1 = 1, NDIM
               DO K2 = 1, NDIM
                  AP(I,J) = AP(I,J) + U(K1,I)*HESS(K1,K2)*U(K2,J)
               ENDDO
            ENDDO
            AP(I,J) = KD(I)*AP(I,J)*KD(J)
         ENDDO
      ENDDO

      IF (EIGENVECTORT) THEN
         CALL DSYEV('V','L',NDIM,AP,NDIM,FRQN,WORK,LWORK,INFO)
      ELSE
         CALL DSYEV('N','L',NDIM,AP,NDIM,FRQN,WORK,LWORK,INFO)
      ENDIF

      ! sn402: commented this out: it doesn't play well with geopt, because this subroutine actually sorts
      ! eigenvalues into descending order but geopt expects ascending order.
!      call eigensort_val_asc(FRQN,AP,NDIM,3*NATOMS)



!  sn402: As of 23/9/16, I'm rationalising the frequency unit conversions. These are now controlled by a single keyword,
! FRQCONV. Check the comments on that keyword in keywords.f for more information.
!      DO I = 1, NDIM
!         IF (FRQN(I) > 0.0D0) THEN
!            FRQN(I) = FRQCNV * SQRT((FRQN(I)))
!         ELSE
!            FRQN(I) = -FRQCNV * SQRT((-FRQN(I)))
!         ENDIF
!      ENDDO
!      FRQN(:) = FRQCNV*SQRT(ABS(FRQN(:)))

!      IF(.FALSE.) THEN
!      FRQN(:) = FRQN(:) * 1.0D26
!      ENDIF
!      print *, 'FRQN'
!      print *, FRQN

!      PRINT *, "CALLED FREQ"
!      PRINT *, FRQN

      IF (EIGENVECTORT) THEN      
         HESS = AP
      ENDIF

      END SUBROUTINE NRMLMD 


      SUBROUTINE COMPUTEINERTIA(RMI, KBLOCK, TMASS) 

        USE KEY, ONLY: NTIPT, PAPT, MULTISITEPYT
        IMPLICIT NONE
        DOUBLE PRECISION :: TMASS, KBLOCK(3,3), RMI(3,3)       

        IF (NTIPT) THEN
           
           CALL INERTIANTIP(RMI, KBLOCK, TMASS)
           
        ELSE IF (PAPT) THEN
!jwrm2> added PAP potential
           CALL INERTIAPAP(RMI, KBLOCK, TMASS)
        ELSE IF (MULTISITEPYT) THEN
!sf344> added PY potential (not generic yet)
           CALL INERTIAPY(RMI, KBLOCK, TMASS)
       
        ENDIF

      END SUBROUTINE COMPUTEINERTIA
