SUBROUTINE MLP3(X,V,ENERGY,GTEST,SECT)
USE COMMONS
USE MODHESS
IMPLICIT NONE
LOGICAL GTEST,SECT
DOUBLE PRECISION X(NMLP), V(NMLP), ENERGY, DUMMY1, DUMMY2, DUMMY3, DUMMY4
DOUBLE PRECISION Y(MLPOUT), PROB(MLPOUT), PMLPOUTJ1
DOUBLE PRECISION DYW1G(MLPHIDDEN), DPCW1BG(MLPOUT,MLPHIDDEN)
DOUBLE PRECISION DYW2G(MLPOUT,MLPHIDDEN,MLPIN), DPCW2BG(MLPHIDDEN,MLPIN), TANHSUM(MLPHIDDEN), SECH2(MLPHIDDEN)
INTEGER MLPOUTJ1, MLPOFFSET
INTEGER GETUNIT, LUNIT, J1, J2, J3, J4, K4, K2, K3, J5

!
! Variables are ordered 
! w^2_{jk} at (j-1)*MLPIN+k
!   up to MLPHIDDEN*MLPIN, then
! w^1_{ij} at MLPHIDDEN*MLPIN + (i-1)*MLPHIDDEN+j
!   up to MLPHIDDEN*MLPIN + MLPOUT*MLPHIDDEN
!

MLPOFFSET=MLPHIDDEN*MLPIN
ENERGY=0.0D0
V(1:NMLP)=0.0D0
! HESS(1:NMLP,1:NMLP)=0.0D0
DO J1=1,MLPDATA
   MLPOUTJ1=MLPOUTCOME(J1)
!  WRITE(MYUNIT,'(A,3I10)') 'J1,MLPDATA,MLPOUTCOME(J1)=',J1,MLPDATA,MLPOUTCOME(J1)
   DO J2=1,MLPHIDDEN
      DUMMY1=0.0D0
      DO J3=1,MLPIN
         DUMMY1=DUMMY1+X((J2-1)*MLPIN+J3)*MLPDAT(J1,J3)
      ENDDO
      TANHSUM(J2)=TANH(DUMMY1) 
      DYW1G(J2)=TANHSUM(J2)
      SECH2(J2)=1.0D0/COSH(DUMMY1)**2 
   ENDDO
   DUMMY3=0.0D0
   DO J4=1,MLPOUT
      DUMMY2=0.0D0
      DO J2=1,MLPHIDDEN
         DO J3=1,MLPIN
            DYW2G(J4,J2,J3)=X( MLPOFFSET + (J4-1)*MLPHIDDEN + J2 ) * MLPDAT(J1,J3)*SECH2(J2)
         ENDDO
         DUMMY2=DUMMY2+X(MLPOFFSET+(J4-1)*MLPHIDDEN+J2)*TANHSUM(J2)
      ENDDO
      Y(J4)=DUMMY2
      DUMMY3=DUMMY3+EXP(DUMMY2)
   ENDDO  
   DO J4=1,MLPOUT
      PROB(J4)=EXP(Y(J4))/DUMMY3
!     WRITE(MYUNIT,'(A,8G15.5)') 'J4,Y(J4),DUMMY3,PROB(J4)=',J4,Y(J4),DUMMY3,PROB(J4)
   ENDDO
!  WRITE(MYUNIT,'(A,I10)') 'MLPOUTJ1=',MLPOUTJ1
   PMLPOUTJ1=PROB(MLPOUTJ1)
!  WRITE(MYUNIT,'(A,G20.10)') 'PMLPOUTJ1',PMLPOUTJ1
   IF (DEBUG) THEN
      WRITE(MYUNIT,'(A,I8,A)') 'MLP3> data point ',J1,' outputs and probabilities:'
      WRITE(MYUNIT,'(8G15.5)') Y(1:MLPOUT),PROB(1:MLPOUT)
   ENDIF
   ENERGY=ENERGY-LOG(PMLPOUTJ1)
   IF (GTEST) THEN
!
! We only need the probability derivative for the probability corresponding to the correct outcome for this data point
!
      DPCW1BG(1:MLPOUT,1:MLPHIDDEN)=0.0D0
      DO J2=1,MLPHIDDEN
         DO J4=1,MLPOUT
            DPCW1BG(J4,J2)=DPCW1BG(J4,J2)-PMLPOUTJ1*PROB(J4)*DYW1G(J2)
         ENDDO
         DPCW1BG(MLPOUTJ1,J2)=DPCW1BG(MLPOUTJ1,J2)+PMLPOUTJ1*DYW1G(J2)
      ENDDO

      DO J3=1,MLPIN
         DO J2=1,MLPHIDDEN
            DUMMY3=0.0D0
            DO J4=1,MLPOUT
               DUMMY3=DUMMY3+PROB(J4)*DYW2G(J4,J2,J3)
            ENDDO
            DPCW2BG(J2,J3)=PMLPOUTJ1*(DYW2G(MLPOUTJ1,J2,J3)-DUMMY3)
         ENDDO
      ENDDO

      DO J4=1,MLPOUT
         DO J2=1,MLPHIDDEN
            V(MLPOFFSET+(J4-1)*MLPHIDDEN+J2)=V(MLPOFFSET+(J4-1)*MLPHIDDEN+J2)-DPCW1BG(J4,J2)/PMLPOUTJ1
         ENDDO
      ENDDO

      DO J3=1,MLPIN
         DO J2=1,MLPHIDDEN
            V((J2-1)*MLPIN+J3)=V((J2-1)*MLPIN+J3)-DPCW2BG(J2,J3)/PMLPOUTJ1
         ENDDO
      ENDDO
   ENDIF
   IF (SECT) THEN
!
! This block w^1 with w^1 is locally symmetric
!
      DO J4=1,MLPOUT ! J4 is beta 
         DO J2=1,MLPHIDDEN ! J2 is gamma
            DO K4=1,J4 ! K4 is alpha
               DO K2=1,MLPHIDDEN ! K2 is epsilon
                  DUMMY1=0.0D0
                  IF ((J4.EQ.MLPOUTJ1).AND.(K4.EQ.MLPOUTJ1)) DUMMY1=1.0D0
                  IF (J4.EQ.MLPOUTJ1) DUMMY1=DUMMY1-PROB(K4)
                  IF (K4.EQ.MLPOUTJ1) DUMMY1=DUMMY1-PROB(J4)
                  IF (K4.EQ.J4) DUMMY1=DUMMY1-PROB(J4)
                  DUMMY1=DUMMY1+2.0D0*PROB(J4)*PROB(K4)
                  HESS(MLPOFFSET+(J4-1)*MLPHIDDEN+J2,MLPOFFSET+(K4-1)*MLPHIDDEN+K2)= &
  &               HESS(MLPOFFSET+(J4-1)*MLPHIDDEN+J2,MLPOFFSET+(K4-1)*MLPHIDDEN+K2) & ! sum over data points
  &               +DPCW1BG(J4,J2)*DPCW1BG(K4,K2)/PMLPOUTJ1**2 &
  &               -DUMMY1*TANHSUM(J2)*TANHSUM(K2)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
      DO J4=1,MLPOUT ! J4 is beta 
         DO J2=1,MLPHIDDEN ! J2 is gamma
            DO K4=1,MLPOUT ! K4 is alpha
               DO K2=1,MLPHIDDEN ! K2 is epsilon
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
! Off-diagonal w^1 with w^2 blocks
!
      DO J3=1,MLPOUT ! J3 is beta for w^1 outputs
         DO J2=1,MLPHIDDEN ! J2 is gamma for w^1 hidden
            DO K4=1,MLPHIDDEN ! K4 is alpha for w^2 hidden
               DUMMY3=0.0D0
               DO J5=1,MLPOUT
                  DUMMY3=DUMMY3+PROB(J5)*X(MLPOFFSET + (J5-1)*MLPHIDDEN + K4) 
               ENDDO
               DO K2=1,MLPIN ! K2 is epsilon for w^2 inputs
                  DUMMY1=0.0D0
                  IF (K4.EQ.J2) DUMMY1=DUMMY1-PMLPOUTJ1*PROB(J3)*MLPDAT(J1,K2)*SECH2(J2)
                  IF ((K4.EQ.J2).AND.(J3.EQ.MLPOUTJ1)) DUMMY1=DUMMY1+PMLPOUTJ1*MLPDAT(J1,K2)*SECH2(J2)

                  DUMMY2=TANHSUM(J2)*PMLPOUTJ1*MLPDAT(J1,K2)*SECH2(K4) &
  &                      *(X(MLPOFFSET+(MLPOUTJ1-1)*MLPHIDDEN+K4)-DUMMY3)
                  DUMMY1=DUMMY1-PROB(J3)*DUMMY2
                  IF (MLPOUTJ1.EQ.J3) DUMMY1=DUMMY1+DUMMY2
                  DUMMY1=DUMMY1-PMLPOUTJ1*PROB(J3)*MLPDAT(J1,K2)*SECH2(K4)*TANHSUM(J2) &
  &                             *(X(MLPOFFSET + (J3-1)*MLPHIDDEN + K4)-DUMMY3)
                  HESS(MLPOFFSET+(J3-1)*MLPHIDDEN+J2,(K4-1)*MLPIN+K2)= &
  &               HESS(MLPOFFSET+(J3-1)*MLPHIDDEN+J2,(K4-1)*MLPIN+K2) & ! sum over data points
  &               +DPCW1BG(J3,J2)*DPCW2BG(K4,K2)/PMLPOUTJ1**2 &
  &               -DUMMY1/PMLPOUTJ1
               ENDDO
            ENDDO
         ENDDO
      ENDDO
!
! diagonal w^2 with w^2 
!
      DO J3=1,MLPIN ! J3 is gamma for w^2 inputs
         DO J2=1,MLPHIDDEN ! J2 is beta for w^2 hidden
            DUMMY2=0.0D0
            DO J5=1,MLPOUT
               DUMMY2=DUMMY2+PROB(J5)*X(MLPOFFSET + (J5-1)*MLPHIDDEN + J2) 
            ENDDO
            DO K4=1,MLPIN ! K4 is epsilon for w^2 inputs
               DO K2=1,J2 ! K2 is alpha for w^2 hidden
                  DUMMY3=0.0D0
                  DO J5=1,MLPOUT
                     DUMMY3=DUMMY3+PROB(J5)*X(MLPOFFSET + (J5-1)*MLPHIDDEN + K2) 
                  ENDDO
                  DUMMY4=0.0D0
                  DO J5=1,MLPOUT
                     DUMMY4=DUMMY4+PROB(J5)*X(MLPOFFSET + (J5-1)*MLPHIDDEN + K2)  & ! take out of loops
  &                                        *X(MLPOFFSET + (J5-1)*MLPHIDDEN + J2)
                  ENDDO
                  DUMMY1=DPCW2BG(K2,K4)*MLPDAT(J1,J3)*SECH2(J2) &
  &                      *(X(MLPOFFSET+(MLPOUTJ1-1)*MLPHIDDEN+J2)-DUMMY2) &
  &                      -PMLPOUTJ1*MLPDAT(J1,J3)*SECH2(J2)*MLPDAT(J1,K4)*SECH2(K2) &
  &                      *(DUMMY4-DUMMY2*DUMMY3)
                  IF (K2.EQ.J2) DUMMY1=DUMMY1-2.0D0*PMLPOUTJ1*MLPDAT(J1,K4)*MLPDAT(J1,J3) &
  &                                     *SECH2(J2)*TANHSUM(J2)*(X(MLPOFFSET + (MLPOUTJ1-1)*MLPHIDDEN + J2)-DUMMY2) 

                  HESS((J2-1)*MLPIN+J3,(K2-1)*MLPIN+K4)= &
  &               HESS((J2-1)*MLPIN+J3,(K2-1)*MLPIN+K4) & ! sum over data points
  &               +DPCW2BG(J2,J3)*DPCW2BG(K2,K4)/PMLPOUTJ1**2 - DUMMY1/PMLPOUTJ1
               ENDDO
            ENDDO
         ENDDO
      ENDDO
   ENDIF
ENDDO

DUMMY1=0.0D0
DO J1=1,NMLP
   DUMMY1=DUMMY1+X(J1)**2
ENDDO

ENERGY=ENERGY/MLPDATA + MLPLAMBDA*DUMMY1
! IF (DEBUG) WRITE(MYUNIT,'(A,G20.10)') 'MLP3> objective function=',ENERGY

IF (GTEST) V(1:NMLP)=V(1:NMLP)/MLPDATA + 2.0D0*MLPLAMBDA*X(1:NMLP)
!
! Symmetrise Hessian here
!
IF (SECT) HESS(1:NMLP,1:NMLP)=HESS(1:NMLP,1:NMLP)/MLPDATA
IF (SECT) THEN
   DO J1=1,NMLP
      HESS(J1,J1)=HESS(J1,J1)+2*MLPLAMBDA
      DO J2=1,J1-1
         HESS(J2,J1)=HESS(J1,J2)
      ENDDO
   ENDDO
ENDIF

END SUBROUTINE MLP3
