SUBROUTINE STEALTHY(XS,TIDU,ENERGY,GTEST,STEST)
    USE COMMONS
    USE MODHESS
    IMPLICIT NONE
    LOGICAL GTEST,STEST
    INTEGER M, M1, M2, M3, J1, J2, I1, I2, L1, L2, KI!, KNUM, GNUM
    DOUBLE PRECISION XS(3*NATOMS), KX(3), KY(3), KZ(3), KN(3), TIDU(3*NATOMS), &
                     CKR(NATOMS), SKR(NATOMS)
    DOUBLE PRECISION PI, KR ,C, S, ENERGY, VF, IM, KM

!    IF (ALLOCATED(HESS)) DEALLOCATE(HESS)
!    ALLOCATE(HESS(3*NATOMS, 3*NATOMS))

!   initialize
    PI=3.141592653589793
    KM=KLIM*KLIM
    KI=DINT(KLIM)
    ENERGY=0

    IF (GTEST) THEN
       TIDU=0
    END IF

    IF (STEST) THEN
       HESS=0
    END IF

!    KNUM=0
!    GNUM=0

!   get the initial box

!   calculate the wave vector
    KX(1)=2*PI/BOXLX
    KX(2)=0
    KX(3)=0
    KY(1)=0
    KY(2)=2*PI/BOXLY
    KY(3)=0
    KZ(1)=0
    KZ(2)=0
    KZ(3)=2*PI/BOXLZ
    VF=BOXLX*BOXLY*BOXLZ
!   WRITE(*,*) KI(1), KI(2), KI(3)

!   periodic boundary condition
    IF (PERIODIC) THEN
       DO J1=1, NATOMS
          J2=3*(J1-1)
          XS(J2+1)=XS(J2+1) - BOXLX*DNINT(XS(J2+1)/BOXLX)
          XS(J2+2)=XS(J2+2) - BOXLY*DNINT(XS(J2+2)/BOXLY)
          XS(J2+3)=XS(J2+3) - BOXLZ*DNINT(XS(J2+3)/BOXLZ)
       END DO
    END IF

!   sum of k
    DO M1=0, KI
       DO M2=-KI, KI
          DO M3=-KI, KI

             IF (M1==0.AND.M2>=0.AND.M3<0) CYCLE
             IF (M1==0.AND.M2<=0.AND.M3<=0) CYCLE

             M=M1*M1+M2*M2+M3*M3
             IF (M>KM.OR.M==0) CYCLE

!             KNUM=KNUM+1
             C=0
             S=0
             KN(1)=KX(1)*M1 + KY(1)*M2 + KZ(1)*M3
             KN(2)=KX(2)*M1 + KY(2)*M2 + KZ(2)*M3
             KN(3)=KX(3)*M1 + KY(3)*M2 + KZ(3)*M3

!            calculation of Re[n(k)] and Im[n(k)]
             DO I1=1, NATOMS
                I2=3*(I1-1)
                KR=KN(1)*XS(I2+1) + KN(2)*XS(I2+2) + KN(3)*XS(I2+3)
                CKR(I1)=COS(KR)
                SKR(I1)=SIN(KR)
                C=C+CKR(I1)
                S=S+SKR(I1)
             END DO

!            E=n(k)^2/Vf
             ENERGY=ENERGY+ SCA*(C*C+S*S)/VF

!            gradient
             IF (GTEST) THEN
!                GNUM=GNUM+1
                DO L1=1, NATOMS
                   L2=3*(L1-1)
                   IM=( C*SKR(L1) - S*CKR(L1) )*SCA
                   TIDU(L2+1)=TIDU(L2+1) - 2*KN(1)*IM/VF
                   TIDU(L2+2)=TIDU(L2+2) - 2*KN(2)*IM/VF
                   TIDU(L2+3)=TIDU(L2+3) - 2*KN(3)*IM/VF
                END DO
             END IF

!            hessian matrix
             IF (STEST) THEN
                CALL STHESS(XS, KN, VF, CKR, C, SKR, S)
             END IF

          END DO
       END DO
    END DO

!    WRITE(*,*) KNUM
!    WRITE(*,*) GNUM
!   change the box back to cubic

END SUBROUTINE

SUBROUTINE STHESS(X, K, V, CK, CSUM, SK, SSUM)
   USE COMMONS
   USE MODHESS
   IMPLICIT NONE
   INTEGER I1, I2, J1, J2
   DOUBLE PRECISION X(3*NATOMS), K(3), RIJ(3), CK(NATOMS), SK(NATOMS)
   DOUBLE PRECISION KRIJ, V, CSUM, SSUM, CI, SI, CF, CSKR

   DO I1=1, NATOMS
      I2=3*(I1-1)

      CI=CSUM-CK(I1)
      SI=SSUM-SK(I1)
      CF=( CK(I1)*CI + SK(I1)*SI )*SCA

      HESS(I2+1, I2+1)=HESS(I2+1, I2+1) + 2*K(1)*K(1)/V * CF
      HESS(I2+1, I2+2)=HESS(I2+1, I2+2) + 2*K(1)*K(2)/V * CF
      HESS(I2+1, I2+3)=HESS(I2+1, I2+3) + 2*K(1)*K(3)/V * CF
      HESS(I2+2, I2+2)=HESS(I2+2, I2+2) + 2*K(2)*K(2)/V * CF
      HESS(I2+2, I2+3)=HESS(I2+2, I2+3) + 2*K(2)*K(3)/V * CF
      HESS(I2+3, I2+3)=HESS(I2+3, I2+3) + 2*K(3)*K(3)/V * CF

      DO J1=I1+1, NATOMS
         J2=3*(J1-1)

         RIJ(1)=X(I2+1) - X(J2+1)
         RIJ(2)=X(I2+2) - X(J2+2)
         RIJ(3)=X(I2+3) - X(J2+3)

         KRIJ=K(1)*RIJ(1) + K(2)*RIJ(2) + K(3)*RIJ(3)
         CSKR= COS(KRIJ)*SCA

         HESS(I2+1, J2+1)=HESS(I2+1, J2+1) + 2*K(1)*K(1)/V * CSKR
         HESS(I2+1, J2+2)=HESS(I2+1, J2+2) + 2*K(1)*K(2)/V * CSKR
         HESS(I2+1, J2+3)=HESS(I2+1, J2+3) + 2*K(1)*K(3)/V * CSKR

         HESS(I2+2, J2+1)=HESS(I2+2, J2+1) + 2*K(2)*K(1)/V * CSKR
         HESS(I2+2, J2+2)=HESS(I2+2, J2+2) + 2*K(2)*K(2)/V * CSKR
         HESS(I2+2, J2+3)=HESS(I2+2, J2+3) + 2*K(2)*K(3)/V * CSKR

         HESS(I2+3, J2+1)=HESS(I2+3, J2+1) + 2*K(3)*K(1)/V * CSKR
         HESS(I2+3, J2+2)=HESS(I2+3, J2+2) + 2*K(3)*K(2)/V * CSKR
         HESS(I2+3, J2+3)=HESS(I2+3, J2+3) + 2*K(3)*K(3)/V * CSKR

      END DO
   END DO

   DO I1=1, NATOMS
      I2=3*(I1-1)

      HESS(I2+2, I2+1)=HESS(I2+1, I2+2)
      HESS(I2+3, I2+1)=HESS(I2+1, I2+3)
      HESS(I2+3, I2+2)=HESS(I2+2, I2+3)

      DO J1=I1+1, NATOMS
         J2=3*(J1-1)

         HESS(J2+1, I2+1)=HESS(I2+1, J2+1)
         HESS(J2+2, I2+1)=HESS(I2+1, J2+2)
         HESS(J2+3, I2+1)=HESS(I2+1, J2+3)

         HESS(J2+1, I2+2)=HESS(I2+2, J2+1)
         HESS(J2+2, I2+2)=HESS(I2+2, J2+2)
         HESS(J2+3, I2+2)=HESS(I2+2, J2+3)

         HESS(J2+1, I2+3)=HESS(I2+3, J2+1)
         HESS(J2+2, I2+3)=HESS(I2+3, J2+2)
         HESS(J2+3, I2+3)=HESS(I2+3, J2+3)

      END DO
   END DO

END SUBROUTINE
