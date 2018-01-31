SUBROUTINE STEALTHY(XS,TIDU,ENERGY,GTEST,STEST)
    USE COMMONS
    USE KEY
    USE MODHESS
    IMPLICIT NONE
    LOGICAL GTEST,STEST
    INTEGER M, M1, M2, M3, J1, J2, I1, I2, L1, L2, KI, I
    DOUBLE PRECISION XS(3*NATOMS), KX(3), KY(3), KZ(3), KN(3), TIDU(3*NATOMS), &
                     CKR(NATOMS), SKR(NATOMS)
    DOUBLE PRECISION PI, KR ,C, S, ENERGY, VF, IM, KM, FORPI, STRMS

    IF (.NOT.ALLOCATED(HESS)) ALLOCATE(HESS(3*NATOMS, 3*NATOMS))

!   initialize
    FORPI=1.0
    PI=4*DATAN(FORPI)
    KM=KLIM*KLIM
    KI=INT(KLIM)
    ENERGY=0
    STM=0

    IF (GTEST) THEN
       TIDU=0
    END IF

    IF (STEST) THEN
       HESS=0
    END IF

!   calculate the wave vector
    KX(1)=2*PI/BULK_BOXVEC(1)
    KX(2)=0
    KX(3)=0
    KY(1)=0
    KY(2)=2*PI/BULK_BOXVEC(2)
    KY(3)=0
    KZ(1)=0
    KZ(2)=0
    KZ(3)=2*PI/BULK_BOXVEC(3)
    VF=BULK_BOXVEC(1)*BULK_BOXVEC(2)*BULK_BOXVEC(3)
!   WRITE(*,*) KI(1), KI(2), KI(3)

!   periodic boundary condition
    IF (BULKT.AND. .NOT. BULKBOXT) THEN
       DO J1=1, NATOMS
          J2=3*(J1-1)
          XS(J2+1)=XS(J2+1) - BULK_BOXVEC(1)*DNINT(XS(J2+1)/BULK_BOXVEC(1))
          XS(J2+2)=XS(J2+2) - BULK_BOXVEC(2)*DNINT(XS(J2+2)/BULK_BOXVEC(2))
          XS(J2+3)=XS(J2+3) - BULK_BOXVEC(3)*DNINT(XS(J2+3)/BULK_BOXVEC(3))
       END DO
    END IF
!    write(*,*) XS
 
!   sum of k
    DO M1=0, KI
       DO M2=-KI, KI
          DO M3=-KI, KI

             IF (M1.EQ.0.AND.M2.GE.0.AND.M3.LT.0) CYCLE
             IF (M1.EQ.0.AND.M2.LE.0.AND.M3.LE.0) CYCLE

             M=M1*M1+M2*M2+M3*M3
             IF (M.GT.KM.OR.M.EQ.0) CYCLE

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
                !if (M2==5) write(*,'(T1,4(I,1X),F,1X,F)') M1, M2, M3, I1, CKR(I1), SKR(I1)

                C=C+CKR(I1)
                S=S+SKR(I1)
!                if (M1==0.AND.M2==-3) write(*,'(T1,F,1X,F)') C, S
             END DO
             !write(*,*) M1, M2, M3
             !write(*,*) CKR
             !write(*,*) SKR
             !if (M1==0.AND.M2==-8) then
             !   write(*,*) CKR(1:6)
             !   write(*,*) SKR(1:6)
             !end if
             !write(*,'(T1,3(I,1X),F,1X,F)') M1, M2, M3, C, S

!            E=n(k)^2/Vf
             ENERGY=ENERGY+ SCA*(C*C+S*S)/VF

!            gradient
             IF (GTEST) THEN
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

             STM=STM+1

          END DO
       END DO
    END DO
    !write(*,*) "STM", STM

    !DO I=1, NATOMS
    !   STRMS=STRMS+TIDU(I)**2
    !ENDDO

    !STRMS=DSQRT(STRMS/NATOMS)

    !write(*,*) "STRMS", STRMS

    IF (STEST.AND.STEALTV) THEN
       CALL STEALTHY_EIGEN_TEST(ENERGY, TIDU)
    ENDIF

END SUBROUTINE

SUBROUTINE STHESS(X, K, V, CK, CSUM, SK, SSUM)
   USE COMMONS
   USE KEY
   USE MODHESS
   IMPLICIT NONE
   INTEGER I1, I2, J1, J2
   DOUBLE PRECISION X(3*NATOMS), K(3), RIJ(3), CK(NATOMS), SK(NATOMS)
   DOUBLE PRECISION KRIJ, V, CSUM, SSUM, CI, SI, CF, CSKR

   DO I1=1, NATOMS
      I2=3*(I1-1)

      CI=CSUM-CK(I1)
      SI=SSUM-SK(I1)
      CF=(CK(I1)*CI + SK(I1)*SI)*SCA

      HESS(I2+1, I2+1)=HESS(I2+1, I2+1) - 2*K(1)*K(1)/V * CF
      HESS(I2+1, I2+2)=HESS(I2+1, I2+2) - 2*K(1)*K(2)/V * CF
      HESS(I2+1, I2+3)=HESS(I2+1, I2+3) - 2*K(1)*K(3)/V * CF
      HESS(I2+2, I2+2)=HESS(I2+2, I2+2) - 2*K(2)*K(2)/V * CF
      HESS(I2+2, I2+3)=HESS(I2+2, I2+3) - 2*K(2)*K(3)/V * CF
      HESS(I2+3, I2+3)=HESS(I2+3, I2+3) - 2*K(3)*K(3)/V * CF

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

SUBROUTINE STEALTHY_EIGEN_TEST(EN, GRAD)
   USE COMMONS
   USE KEY
   USE MODHESS
   IMPLICIT NONE
   INTEGER I, NZEV, INFO, LSTECUT, HSTECUT
   DOUBLE PRECISION GRAD(3*NATOMS), EIGEN(3*NATOMS), STTEMP(9*NATOMS)
   DOUBLE PRECISION STRMS, EN

   CALL DSYEV('N', 'U', 3*NATOMS, HESS, 3*NATOMS, EIGEN, STTEMP, 9*NATOMS, INFO)
   STRMS=0

   !write(*,*) EIGEN
   DO I=1, NATOMS
      STRMS=STRMS+GRAD(I)**2
   ENDDO

   STRMS=DSQRT(STRMS/NATOMS)
   !write(*,'(T2,A,G20.10)') 'stealthy> Calculated RMS',STRMS

   IF (STRMS.LE.GMAX) THEN
      IF (EN.LT.1.0D-8) THEN
         NZEV=NATOMS*3-STM*2
         PRINT '(A)',' stealthy> Ground state'
         !write(*,'(T2,A,I,G20.10)') 'stealthy> Test highest zero eigenvalue', NZEV, EIGEN(NZEV)
         !write(*,'(T2,A,I,G20.10)') 'stealthy> Test lowest non-zero eigenvalue', NZEV+1, EIGEN(NZEV+1)
         IF (EIGEN(1).LT.-1.0D-8) THEN
            PRINT '(A,G20.10)',' stealthy> Error with lowest eigenvalue', EIGEN(1)
         ENDIF
         IF (EIGEN(NZEV+1).LT.1.0D-8) THEN
            PRINT '(A,G20.10)',' stealthy> Error with lowest non-zero eigenvalue',  EIGEN(NZEV+1)
         ENDIF
         IF (EIGEN(NZEV).GT.1.0D-8) THEN
            PRINT '(A,G20.10)',' stealthy> Error with highest zero eigenvalue', EIGEN(NZEV)
         ENDIF
      ELSE
         NZEV=3
         PRINT '(A)',' stealthy> Local minimum'
         CALL ST_EIGEN_SORT(EIGEN, NZEV, LSTECUT, HSTECUT)
         !write(*,'(T2,A,I,G20.10)') 'stealthy> Test highest zero eigenvalue', LSTECUT, EIGEN(LSTECUT)
         !write(*,'(T2,A,I,G20.10)') 'stealthy> Test lowest non-zero eigenvalue', HSTECUT, EIGEN(HSTECUT)
         IF (EIGEN(1).LT.-1.0D-8) THEN
            PRINT '(A,G20.10)',' stealthy> Error with lowest eigenvalue', EIGEN(1)
         ENDIF
         IF (EIGEN(HSTECUT).LT.1.0D-8) THEN
            PRINT '(A,G20.10)',' stealthy> Error with lowest non-zero eigenvalue', EIGEN(HSTECUT)
         ENDIF
         IF (EIGEN(LSTECUT).GT.1.0D-8) THEN
            PRINT '(A,G20.10)',' stealthy> Error with highest zero eigenvalue', EIGEN(LSTECUT)
         ENDIF
      ENDIF
   ELSE
      NZEV=3
      PRINT '(A)',' stealthy> Non-minimum configuration'
      CALL ST_EIGEN_SORT(EIGEN, NZEV, LSTECUT, HSTECUT)
      !write(*,'(T2,A,I,G20.10)') 'stealthy> Test highest zero eigenvalue', LSTECUT, EIGEN(LSTECUT)
      !write(*,'(T2,A,I,G20.10)') 'stealthy> Test lowest non-zero eigenvalue', HSTECUT, EIGEN(HSTECUT+1)
      IF (EIGEN(HSTECUT).LT.1.0D-8) THEN
         PRINT '(A,G20.10)',' stealthy> Error with lowest non-zero eigenvalue', EIGEN(HSTECUT)
      ENDIF
      IF (EIGEN(LSTECUT).GT.1.0D-8) THEN
         PRINT '(A,G20.10)',' stealthy> Error with highest zero eigenvalue', EIGEN(LSTECUT)
      ENDIF
   ENDIF

END SUBROUTINE

SUBROUTINE ST_EIGEN_SORT(EVARY, NLING, ELCUT, EHCUT)
   USE COMMONS
   IMPLICIT NONE
   LOGICAL EQYN
   INTEGER I, J, K, L, ELCUT, EHCUT, NLING
   INTEGER,ALLOCATABLE :: BNNUM(:)
   DOUBLE PRECISION EVARY(3*NATOMS), EVABS(3*NATOMS)
   DOUBLE PRECISION MYMIN

   ALLOCATE(BNNUM(NLING+1))
   EVABS=DABS(EVARY)
   BNNUM=3*NATOMS
   EQYN=.FALSE.
   
   DO I=1, NLING+1
      MYMIN=EVABS(3*NATOMS)
      DO J=1, 3*NATOMS-1
         DO K=1, I
            IF (J.EQ.BNNUM(K)) EQYN=.TRUE.
         ENDDO
         IF (EQYN) THEN
            EQYN=.FALSE.
            CYCLE
         ENDIF
         IF (EVABS(J).LT.MYMIN) THEN
            MYMIN=EVABS(J)
            BNNUM(I)=J
         ENDIF
      ENDDO
      IF (I.EQ.NLING) ELCUT=BNNUM(I)
      IF (I.EQ.NLING+1) EHCUT=BNNUM(I)
   ENDDO

END SUBROUTINE
