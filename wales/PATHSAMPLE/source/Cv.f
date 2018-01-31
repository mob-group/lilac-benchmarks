      SUBROUTINE CV
      USE COMMONS
      IMPLICIT NONE
      INTEGER J1
      DOUBLE PRECISION DUMMY, CVTEMP, Z0, Z1, Z2, NUMIN, Z0R, Z1R, Z2R


      PRINT '(A,2F15.5,A,F20.10)','Cv> Calculating Cv for temperature range ',CVTMIN,CVTMAX,' increment ',CVTINC
      NUMIN=1.0D100
      DO J1=1,NMIN
         IF (FVIBMIN(J1).LT.NUMIN) NUMIN=FVIBMIN(J1)
      ENDDO

      CVTEMP=CVTMIN
      IF (.NOT.CVMINIMAT) THEN
         OPEN(UNIT=1,FILE='Cv.out',STATUS='UNKNOWN')
         OPEN(UNIT=2,FILE='Cv.config.out',STATUS='UNKNOWN')
      ENDIF
20    CONTINUE     
      Z0=0.0D0
      Z1=0.0D0
      Z2=0.0D0
      Z0R=0.0D0
      Z1R=0.0D0
      Z2R=0.0D0
      DO J1=1,NMIN
         DUMMY=EXP(-(EMIN(J1)-EMIN(1))/CVTEMP+(NUMIN-FVIBMIN(J1))/2)/HORDERMIN(J1)
         Z0=Z0+DUMMY
         Z1=Z1+EMIN(J1)*DUMMY
         Z2=Z2+EMIN(J1)**2*DUMMY
         DUMMY=DUMMY/SQRT(IXMIN(J1)*IYMIN(J1)*IZMIN(J1))
         Z0R=Z0R+DUMMY
         Z1R=Z1R+EMIN(J1)*DUMMY
         Z2R=Z2R+EMIN(J1)**2*DUMMY
      ENDDO
      WRITE(1,'(3G20.10)') CVTEMP,KAPPA-Z1**2/(CVTEMP*Z0)**2+Z2/(Z0*CVTEMP**2),
     &                            KAPPA-Z1R**2/(CVTEMP*Z0R)**2+Z2R/(Z0R*CVTEMP**2)   
      WRITE(2,'(3G20.10)') CVTEMP,KAPPA/2.0D0-Z1**2/(CVTEMP*Z0)**2+Z2/(Z0*CVTEMP**2),
     &                            KAPPA/2.0D0-Z1R**2/(CVTEMP*Z0R)**2+Z2R/(Z0R*CVTEMP**2)   
 
      CVTEMP=CVTEMP+CVTINC
      IF (CVTEMP.GT.CVTMAX) THEN
         IF (.NOT.CVMINIMAT) THEN
            CLOSE(1)
            CLOSE(2)
         ENDIF
         RETURN
      ENDIF
      GOTO 20
    
      RETURN
      END
