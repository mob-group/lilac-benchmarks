      SUBROUTINE DOS
      USE COMMONS
      IMPLICIT NONE
      INTEGER J1
      DOUBLE PRECISION DUMMY, DOSE

      PRINT '(A,2F15.5,A,F20.10)','DOS> Calculating total harmonic energy density of states for energy range ',DOSEMIN,DOSEMAX,
     &                            ' increment ',DOSEINC
      DOSE=DOSEMIN
      OPEN(UNIT=1,FILE='DOS.out',STATUS='UNKNOWN')
20    CONTINUE     
      DUMMY=0.0D0
      DO J1=1,NMIN
         IF (DOSE.GT.EMIN(J1)) DUMMY=DUMMY+EXP((KAPPA-1)*LOG(DOSE-EMIN(J1)) - FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1)))
      ENDDO
      IF (DUMMY.GT.0.0D0) WRITE(1,'(3G20.10)') DOSE,DUMMY,LOG(DUMMY)
 
      DOSE=DOSE+DOSEINC
      IF (DOSE.GT.DOSEMAX) THEN
         CLOSE(1)
         STOP
      ENDIF
      GOTO 20
    
      RETURN
      END
