!
!  Persistence analysis using getbarrier2 routine.
!  Persistence diagram is a plot of local minimum energy
!  versus transition state energy where that minimum joins a
!  superbasin with a lower minimum in it. The corresponding downhill 
!  barrier height is called the culminance of the minimum in question.
!
! Approximate method:
!  Usual superbasin analysis at fine grained total energies.
!  Check if minimum has a threshold assigned.
!  If not, check if it is the lowest minimum in the superbasin.
!  If not, assign that threshold as its culminance, and mark it done.
!
! Exact method uses a sorted list of ts energies and a Union-Find algorithm 
! (see Introduction to Algorithms; Cormen, Leiserson, Rivest and Stein; Chapter 21 Data structures for disjoint sets)
!
SUBROUTINE PERSISTENCE
USE COMMONS,ONLY : NMIN, NTS, EMIN, DEBUG, TSTHRESH, FVIBMIN, HORDERMIN, TEMPERATURE, TSTHRESH, &
  &                PEQTHRESH, PERTHRESH, PERSISTAPPROXT
USE UTILS,ONLY : GETUNIT
IMPLICIT NONE
DOUBLE PRECISION HIGHESTTS, LOWESTTARG, ETHRESH, BARRIER(NMIN), GMIN, SLOPE, INTERCEPT, CORRELATION, EDISC(NMIN)
DOUBLE PRECISION SUMX, SUMY, SUMX2, SUMY2, SUMXY, XTEMP(NMIN), PEQ(NMIN), PFNORM, DUMMY, LOCPFMIN(NMIN)
INTEGER J1, BASIN(NMIN), NBASIN, MINTARG, NGMIN, NCONNECTED, MINID(NMIN), LUNIT
LOGICAL CHANGED, BASINT(NMIN)

IF (PERSISTAPPROXT) THEN
  CALL GETBARRIER2(BARRIER,NGMIN,GMIN)
ELSE
  CALL GETEXACTBARRIER(BARRIER,NGMIN,GMIN)
ENDIF
PRINT '(A,I6,G20.10)','persistence> global minimum position and energy: ',NGMIN,GMIN
SUMX=0.0D0 
SUMX2=0.0D0 
SUMXY=0.0D0 
SUMY=0.0D0 
SUMY2=0.0D0 

NCONNECTED=0
EDISC(1:NMIN)=-1.0D0
XTEMP(1:NMIN)=-1.0D0
DO J1=1,NMIN
   IF (BARRIER(J1).GT.0.0D0) THEN
      NCONNECTED=NCONNECTED+1
      MINID(NCONNECTED)=J1
      XTEMP(NCONNECTED)=BARRIER(J1)/SQRT(2.0D0) ! minimum distance of point to diagonal line
      EDISC(J1)=EMIN(J1)+BARRIER(J1)
      SUMX=SUMX+EMIN(J1)
      SUMX2=SUMX2+EMIN(J1)**2
      SUMXY=SUMXY+EMIN(J1)*EDISC(J1)
      SUMY=SUMY+EDISC(J1)
      SUMY2=SUMY2+EDISC(J1)**2
!     PRINT '(A,I6,3G20.10)','J1,barrier,sumx,sumy=',J1,barrier(J1),sumx,sumy
   ENDIF
ENDDO
!
! Regression analysis - best straight line fit.
!
SLOPE=(NCONNECTED * sumxy  -  sumx * sumy) / (NCONNECTED * sumx2 - sumx**2)
INTERCEPT=(sumy * sumx2  -  sumx * sumxy) / (NCONNECTED * sumx2  -  sumx**2)
CORRELATION=(sumxy - sumx * sumy / NCONNECTED) /sqrt(MAX(1.0D-100,(sumx2 - sumx**2/NCONNECTED) * (sumy2 - sumy**2/NCONNECTED)))
PRINT '(A,G20.10,A,G20.10)','persistence> best fit straight line is EDISC=',slope,'* EMIN + ',intercept
PRINT '(A,G20.10,A,I8,A)','persistence> correlation coefficient=',CORRELATION,' for ',NCONNECTED,' connected minima'

PRINT '(A)','rank (barrier)      energy       energy+barrier           barrier       min           deviation  '
CALL SORT(NCONNECTED,NCONNECTED,XTEMP,MINID)
LUNIT=GETUNIT()
OPEN(LUNIT,FILE='min.persist',STATUS='UNKNOWN')
DO J1=1,NCONNECTED
   PRINT '(I8,3F20.10,I8,F20.10)',J1,EMIN(MINID(J1)),EMIN(MINID(J1))+BARRIER(MINID(J1)),BARRIER(MINID(J1)),MINID(J1), &
  &                        XTEMP(J1)
   WRITE(LUNIT,'(I8,3F20.10,I8,F20.10)') J1,EMIN(MINID(J1)),EMIN(MINID(J1))+BARRIER(MINID(J1)),BARRIER(MINID(J1)),MINID(J1), &
  &                        XTEMP(J1)
ENDDO
CLOSE(LUNIT)

PFNORM=0.0D0
DO J1=1,NMIN
   LOCPFMIN(J1) = -EMIN(J1)/TEMPERATURE - FVIBMIN(J1)/2.0D0 - LOG(1.0D0*HORDERMIN(J1))
   PFNORM=PFNORM+EXP(LOCPFMIN(J1)-LOCPFMIN(1))
ENDDO
PFNORM=LOG(PFNORM)+LOCPFMIN(1)
DUMMY=0.0D0
DO J1=1,NMIN
   PEQ(J1)=EXP(LOCPFMIN(J1)-PFNORM)
   DUMMY=DUMMY+PEQ(J1)
ENDDO
PRINT '(A,G20.10)','sum of equilibrium occupation probabilities=',DUMMY

PRINT '(3(A,G20.10))','data for connected minima with Peq > ',PEQTHRESH,' energy+barrier < ',TSTHRESH,&
                     & ' and deviation > ',PERTHRESH
LUNIT=GETUNIT()
OPEN(LUNIT,FILE='min.persist.select',STATUS='UNKNOWN')
DO J1=1,NCONNECTED
   IF ((PEQ(MINID(J1)).GT.PEQTHRESH).AND.(EMIN(MINID(J1))+BARRIER(MINID(J1)).LT.TSTHRESH).AND.(XTEMP(J1).GT.PERTHRESH)) THEN
      PRINT '(I8,3F20.10,I8,F20.10)',J1,EMIN(MINID(J1)),EMIN(MINID(J1))+BARRIER(MINID(J1)),BARRIER(MINID(J1)),MINID(J1), &
  &                        XTEMP(J1)
      WRITE(LUNIT,'(I8,3F20.10,I8,F20.10)') J1,EMIN(MINID(J1)),EMIN(MINID(J1))+BARRIER(MINID(J1)),BARRIER(MINID(J1)),MINID(J1), &
  &                        XTEMP(J1)
   ENDIF
ENDDO
CLOSE(LUNIT)
LUNIT=GETUNIT()
OPEN(LUNIT,FILE='min.include',STATUS='UNKNOWN')
WRITE(LUNIT,'(I6)') NGMIN
DO J1=1,NCONNECTED
   IF ((PEQ(MINID(J1)).GT.PEQTHRESH).AND.(EMIN(MINID(J1))+BARRIER(MINID(J1)).LT.TSTHRESH).AND.(XTEMP(J1).GT.PERTHRESH)) THEN
      WRITE(LUNIT,'(I6)') MINID(J1)
   ENDIF
ENDDO
CLOSE(LUNIT)
LUNIT=GETUNIT()
OPEN(LUNIT,FILE='persist.trvalfile',STATUS='UNKNOWN')
DO J1=1,NMIN
   IF ((PEQ(J1).GT.PEQTHRESH).AND.(BARRIER(J1).GT.(PERTHRESH*SQRT(2.0D0))).AND.(EMIN(J1)+BARRIER(J1).LT.TSTHRESH)) THEN 
     write(LUNIT,*) BARRIER(J1)
   else
     write(LUNIT,*) "-1.0"
   endif
ENDDO
CLOSE(LUNIT)

END SUBROUTINE PERSISTENCE
