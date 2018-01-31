      PROGRAM GMINCONV3
      IMPLICIT NONE
      INTEGER J1, MMAX, N, NPOW
      PARAMETER(MMAX=10000,NPOW=10)
      DOUBLE PRECISION X(MMAX), XMAX, XMIN, XSIG(NPOW), XMEAN(NPOW)

      XMAX=-1.0D100
      XMIN=1.0D100
      XMEAN(1)=0.0D0
      XSIG(1)=0.0D0
      
      DO J1=1,MMAX
         READ(*,*,END=20) X(J1)
         XMEAN(1)=XMEAN(1)+ X(J1)
         XSIG(1)=XSIG(1)+ X(J1)**2
         IF (X(J1).GT.XMAX) XMAX=X(J1)
         IF (X(J1).LT.XMIN) XMIN=X(J1)
      ENDDO
20    N=J1-1
      XMEAN(1)=XMEAN(1)/N
      XSIG(1)=DSQRT((XSIG(1)-N*XMEAN(1)**2)/(N-1))
      WRITE(*,'(A,4F13.2)') 'lowest, highest, mean and sd of energies: ',XMIN,XMAX,XMEAN(1),XSIG(1)

      STOP
      END
