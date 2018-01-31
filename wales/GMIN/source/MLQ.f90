SUBROUTINE MLQ(X,V,ENERGY,GTEST,SECT)
USE MODHESS
USE COMMONS
IMPLICIT NONE
LOGICAL GTEST,SECT
DOUBLE PRECISION X(NMLQ), V(NMLQ), ENERGY, DUMMY1, DUMMY2, DUMMY3
DOUBLE PRECISION Y(MLQOUT), PROB(MLQOUT), DMAX, PMLQOUTJ1
DOUBLE PRECISION DYDW(NMLQ,MLQOUT), DPDW(NMLQ,MLQOUT), PSUM(NMLQ)
INTEGER MLQOUTJ1, MLQOFFSET, NDUMMY
INTEGER J1, J2, J3, J4
!
! Variables are ordered 
! w^0 at 1
! w^1_j at j+1
! w^2_jk with k >= j at 1+NMLQIN+(j*(2*NMLQIN-j+1))/2+1-j+k = 1+k+(j*(2*NMLQIN-j+1))/2
! for input 1, then 2 offset by NMLQOUT*(1+(NMLQIN*(NMLQIN+3))/2), then 3, ...
!

MLQOFFSET=1+(MLQIN*(MLQIN+3))/2
ENERGY=0.0D0
V(1:NMLQ)=0.0D0
IF (SECT) HESS(1:NMLQ,1:NMLQ)=0.0D0
!
! J1 is data item from 1 to MLQDATA
! MLQOUTJ1=MLQOUTCOME(J1) is the outcome for data entry J1
!
DO J1=1,MLQDATA
   DYDW(1:NMLQ,1:MLQOUT)=0.0D0
   MLQOUTJ1=MLQOUTCOME(J1)
!  WRITE(MYUNIT,'(A,2I10)') 'J1,outcome=',J1,MLQOUTCOME(J1)
   DMAX=-1.0D100
   NDUMMY=0
   DO J4=1,MLQOUT
      DUMMY2=X(MLQOFFSET*(J4-1)+1) ! initialises the output value J4 for data point J1
      NDUMMY=NDUMMY+1
      DYDW(MLQOFFSET*(J4-1)+1,J4)=1.0D0 ! not summed over data points J1
      DO J2=1,MLQIN
         DUMMY2=DUMMY2+X(MLQOFFSET*(J4-1)+1+J2)*MLQDAT(J1,J2)
         NDUMMY=NDUMMY+1
         DYDW(MLQOFFSET*(J4-1)+1+J2,J4)=MLQDAT(J1,J2)
!        PRINT '(A,2I8)','w1 loop NDUMMY and index are ',NDUMMY,MLQOFFSET*(J4-1)+1+J2
      ENDDO
      DO J2=1,MLQIN      ! j
         DO J3=J2,MLQIN  ! k >= j
            DUMMY2=DUMMY2+MLQDAT(J1,J2)*MLQDAT(J1,J3)*X(MLQOFFSET*(J4-1)+1+J3+(J2*(2*MLQIN-J2+1))/2)
            NDUMMY=NDUMMY+1
            DYDW(MLQOFFSET*(J4-1)+1+J3+(J2*(2*MLQIN-J2+1))/2,J4)=MLQDAT(J1,J2)*MLQDAT(J1,J3)
!           PRINT '(A,2I8)','w2 loop NDUMMY and index are ',NDUMMY,MLQOFFSET*(J4-1)+1+J3+(J2*(2*MLQIN-J2+1))/2
         ENDDO
      ENDDO
      IF (DUMMY2.GT.DMAX) DMAX=DUMMY2
      Y(J4)=DUMMY2 ! output J4 for data point J1 
   ENDDO  
!  PRINT '(A,2I10)','NDUMMY should equal MLQOUT*MLQOFFSET=',NDUMMY,MLQOUT*MLQOFFSET
!
! Factor DMAX out of the exponentials to prevent overflow.
!
   DUMMY3=0.0D0
   DO J4=1,MLQOUT
      Y(J4)=Y(J4)-DMAX
      DUMMY3=DUMMY3+EXP(Y(J4))
   ENDDO
   DO J4=1,MLQOUT
      PROB(J4)=EXP(Y(J4))/DUMMY3
   ENDDO
   PMLQOUTJ1=PROB(MLQOUTJ1) ! this is p_c(d), the probability predicted for the known outcome for data point d
   IF (PMLQOUTJ1.LE.0.0D0) THEN
      IF (DEBUG) WRITE(MYUNIT,'(A,G20.10)') 'WARNING PROB=',PROB(MLQOUTJ1)
      PMLQOUTJ1=1.0D-100
!     WRITE(MYUNIT,'(I10,G20.10)') (J4,PROB(J4),J4=1,MLQOUT)
!     STOP
   ENDIF
   ENERGY=ENERGY-LOG(PMLQOUTJ1) ! summing over contributions for data points, indexed by J1
!  WRITE(MYUNIT,'(A,I8,2G20.10)') 'J1,PMLQOUTJ1,ENERGY=',J1,PMLQOUTJ1,ENERGY

   IF (GTEST) THEN
!
! We only need the probability derivatives for the probability corresponding to the correct outcome for this data point
!
      DO J3=1,NMLQ ! for each variable w
         DUMMY3=0.0D0
         DO J4=1,MLQOUT ! the only non-zero terms in this sum will correspond to variables associated with the given output
            DUMMY3=DUMMY3+PROB(J4)*DYDW(J3,J4)
         ENDDO
         PSUM(J3)=DUMMY3 
         V(J3)=V(J3)-(DYDW(J3,MLQOUTJ1)-DUMMY3)
      ENDDO

   ENDIF

   IF (SECT) THEN

      IF (.NOT.GTEST) THEN ! need PSUM
         DO J3=1,NMLQ ! for each variable w
            DUMMY3=0.0D0
            DO J4=1,MLQOUT ! the only non-zero terms in this sum will correspond to variables associated with the given output
               DUMMY3=DUMMY3+PROB(J4)*DYDW(J3,J4)
            ENDDO
            PSUM(J3)=DUMMY3 
         ENDDO
      ENDIF

      DO J3=1,NMLQ
         DO J4=1,MLQOUT 
            DPDW(J3,J4)=PROB(J4)*(DYDW(J3,J4)-PSUM(J3))
         ENDDO
      ENDDO

      DO J2=1,NMLQ  ! alpha
         DO J3=J2,NMLQ  ! beta
            DUMMY3=0.0D0
            DO J4=1,MLQOUT
               DUMMY3=DUMMY3+DPDW(J3,J4)*DYDW(J2,J4)
            ENDDO
            HESS(J3,J2)=HESS(J3,J2) + DUMMY3
! &                                +DPDW(J3,MLQOUTJ1)*DPDW(J2,MLQOUTJ1)/PMLQOUTJ1**2 &
! &                                -(DYDW(J3,MLQOUTJ1)-PSUM(J3))*(DYDW(J2,MLQOUTJ1)-PSUM(J2)) &
! &                                -DYDW(J3,MLQOUTJ1)*DYDW(J2,MLQOUTJ1) &
! &                                +DYDW(J3,MLQOUTJ1)*PSUM(J2)+DYDW(J2,MLQOUTJ1)*PSUM(J3)-PSUM(J3)*PSUM(J2) &
         ENDDO
      ENDDO

   ENDIF

ENDDO ! end of sum over data points, counter J1, limit MLQDATA

DUMMY1=0.0D0
DO J1=1,NMLQ
   DUMMY1=DUMMY1+X(J1)**2
ENDDO

ENERGY=ENERGY/MLQDATA + MLQLAMBDA*DUMMY1
! IF (DEBUG) WRITE(*,'(A,G20.10)') 'MLQ> objective function=',ENERGY

IF (GTEST) V(1:NMLQ)=V(1:NMLQ)/MLQDATA 
IF (GTEST) V(1:NMLQ)=V(1:NMLQ) + 2.0D0*MLQLAMBDA*X(1:NMLQ)
!
! Symmetrise Hessian here
!
IF (SECT) HESS(1:NMLQ,1:NMLQ)=HESS(1:NMLQ,1:NMLQ)/MLQDATA
IF (SECT) THEN
   DO J1=1,NMLQ
      HESS(J1,J1)=HESS(J1,J1)+2*MLQLAMBDA
   ENDDO
   DO J1=1,NMLQ
      DO J2=1,J1-1
         HESS(J2,J1)=HESS(J1,J2)
      ENDDO
   ENDDO
ENDIF

END SUBROUTINE MLQ

