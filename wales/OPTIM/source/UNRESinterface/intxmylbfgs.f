C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
C
C        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
C                          JORGE NOCEDAL
C                        *** July 1990 ***
C
C        Line search removed plus small modifications, DJW 2001
C        This routine is for minimising a Rayleigh-Ritz ratio
C        in BFGS/BFGS transition state searches.
C        Called if INTMINT is true, COORDS passed from intbfgsts 
C        via intbeig contains internals, N is NINTS. jmc
C
      SUBROUTINE INTXMYLBFGS(N,M,X,DIAGCO,DIAG,EPS,W,EENERGY,COORDS,ITDONE,ITMAX,ENERGY,PTEST,CONVERGED)
      USE COMMONS
      USE KEY
      USE MODTWOEND
      IMPLICIT NONE

      INTEGER N,M,J1,J2,ITMAX,ITDONE,NFAIL
C JMC      DOUBLE PRECISION X(N),G(3*NATOMS),DIAG(N),W(N*(2*M+1)+2*M),SLENGTH,DDOT,DUMMY2
      DOUBLE PRECISION X(3*NATOMS),G(3*NATOMS),DIAG(N),W(N*(2*M+1)+2*M),SLENGTH,DDOT,DUMMY2
      DOUBLE PRECISION EPS,DUMMY1,ENERGY,ENEW,GNEW(3*NATOMS),RMS,EENERGY
      LOGICAL DIAGCO, PTEST
      DOUBLE PRECISION GNORM,STP,YS,YY,SQ,YR,BETA,COORDS(3*NATOMS),OVERLAP,DOT1,DOT2
      INTEGER ITER,POINT,ISPT,IYPT,BOUND,NPT,CP,I,INMC,IYCN,ISCN,NDECREASE, ISTAT
      LOGICAL MFLAG, FAILED, CONVERGED

      ITER=0
      ITDONE=0
      NFAIL=0
      FAILED=.FALSE.

C jmc X is displacement vector (converging to eigenvector), COORDS is coords, eenergy is real energy, 
C G is gradient of lambda wrt x, energy is (almost) lambda, rms is rms G

      CALL INTSECDIAG(X,COORDS,EENERGY,G,ENERGY,.TRUE.,RMS)

      IF (PTEST) WRITE(*,'(A,2G20.10,A,I6,A)') 'intxmylbfgs> Eigenvalue and RMS force=',ENERGY,RMS,' after ',ITDONE,' steps'

10    CALL FLUSH(6)
      MFLAG=.FALSE.
      IF (RMS.LE.EPS) THEN
         MFLAG=.TRUE.
         IF (MFLAG) THEN
            FIXIMAGE=.FALSE.
            IF (PTEST) WRITE(*,'(A,I4,A,F15.7,A,F15.7)') 'intxmylbfgs> Smallest eigenvalue converged in ',ITDONE,
     1              ' steps. Eigenvalue=',ENERGY,' RMS force=',RMS
            CALL FLUSH(6)
            IF (PTEST.AND.(ITER.GT.0))
     1          WRITE(*,'(A,F20.10)') 'intxmylbfgs> Diagonal inverse Hessian elements are now ',DIAG(1)
            CONVERGED=.TRUE.
            RETURN
         ENDIF
      ENDIF

      IF ((ITDONE.EQ.ITMAX).OR.FAILED) THEN
         FIXIMAGE=.FALSE.
         IF (PTEST) WRITE(*,'(A,G15.7,A,G15.7)') 'intxmylbfgs> **WARNING - Smallest eigenvalue did not converge, value=',ENERGY,
     1           ' RMS force=',RMS
         IF (PTEST) WRITE(*,'(A,F20.10)') 'intxmylbfgs> Diagonal inverse Hessian elements are now ',DIAG(1)
         CONVERGED=.FALSE.
         RETURN
      ENDIF

      IF (ITER.EQ.0) THEN
         IF (N.LE.0.OR.M.LE.0) THEN
            WRITE(*,240)
 240        FORMAT('intxmylbfgs> IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
            STOP
         ENDIF
         POINT=0
         MFLAG=.FALSE.
         IF (DIAGCO) THEN
            PRINT*,'using estimate of the inverse diagonal elements'
            DO I=1,N
               IF (DIAG(I).LE.0.0D0) THEN
                  WRITE(*,235) I
 235              FORMAT(' THE',I5,'-TH DIAGONAL ELEMENT OF THE',/,
     1                   ' INVERSE HESSIAN APPROXIMATION IS NOT POSITIVE')
                  STOP
               ENDIF
            ENDDO
         ELSE
C           INQUIRE(FILE='diag',EXIST=YESNO)
C           IF (YESNO) THEN
C              OPEN(UNIT=34,FILE='diag',STATUS='OLD')
C              READ(34,*) (DIAG(I),I=1,N)
C              PRINT*,'intxmylbfgs> diag read in LBFGS'
C              WRITE(*,'(6F15.5)') (DIAG(I),I=1,N)
C           ELSE
            DO I=1,N
               DIAG(I)=XDGUESS
            ENDDO
         ENDIF
C
C     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
C     ---------------------------------------
C     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
C         OTHER TEMPORARY INFORMATION.
C     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
C     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
C         IN THE FORMULA THAT COMPUTES H*G.
C     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
C         STEPS.
C     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
C         GRADIENT DIFFERENCES.
C
C     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
C     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
C
         ISPT= N+2*M
         IYPT= ISPT+N*M
C
C  NR step for diagonal inverse Hessian
C
         DO I=1,N
            W(ISPT+I)= -G(I)*DIAG(I)
            W(I)= -G(I)*DIAG(I)
         ENDDO
         GNORM= DSQRT(DDOT(N,G,1,G,1))
C
C  Make the first guess for the step length cautious.
C
         IF (GNORM.EQ.0.0D0) THEN
            GNORM=1.0D0 ! exact zero is presumably wrong!
            PRINT '(A)','WARNING - GNORM was zero in xmylbfgs, resetting to one'
         ENDIF
         STP=MIN(GNORM,1.0D0/GNORM)
      ELSE 
         BOUND=ITER
         IF (ITER.GT.M) BOUND=M
         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
         IF (YS.EQ.0.0D0) YS=1.0D0
C
C  Update estimate of diagonal inverse Hessian elements
C  We divide by both YS and YY at different points, so
C  they had better not be zero!
C
         IF (.NOT.DIAGCO) THEN
            YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
            IF (YY.EQ.0.0D0) YY=1.0D0
            DO I=1,N
               DIAG(I)= YS/YY
C              DIAG(I)= ABS(YS/YY)
            ENDDO
         ELSE
            PRINT*,'using estimate of the inverse diagonal elements'
            DO I=1,N
               IF (DIAG(I).LE.0.0D0) THEN
                  WRITE(*,235) I
                  STOP
               ENDIF
            ENDDO
         ENDIF
C
C     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
C     "Updating quasi-Newton matrices with limited storage",
C     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
C     ---------------------------------------------------------
C
         CP= POINT
         IF (POINT.EQ.0) CP=M
         W(N+CP)= 1.0D0/YS
         DO I=1,N
            W(I)= -G(I)
         ENDDO
         CP= POINT
         DO I= 1,BOUND
            CP=CP-1
            IF (CP.EQ. -1)CP=M-1
            SQ= DDOT(N,W(ISPT+CP*N+1),1,W,1)
            INMC=N+M+CP+1
            IYCN=IYPT+CP*N
            W(INMC)= W(N+CP+1)*SQ
            CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
         ENDDO
        
         DO I=1,N
            W(I)=DIAG(I)*W(I)
         ENDDO

         DO I=1,BOUND
            YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
            BETA= W(N+CP+1)*YR
            INMC=N+M+CP+1
            BETA= W(INMC)-BETA
            ISCN=ISPT+CP*N
            CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
            CP=CP+1
            IF (CP.EQ.M) CP=0
         ENDDO
         STP=1.0D0
      ENDIF
C
C  Orthogonalisation of the search direction to translations and rotations
C  appears to make no difference.
C
C      CALL ORTHOGOPT(W,COORDS,.TRUE.)

C
C  Store the new search direction
C
C     IF (ITER.NE.0) THEN
         DO I=1,N
            W(ISPT+POINT*N+I)= W(I)
         ENDDO
C     ENDIF

C     OVERLAP=DDOT(N,G,1,W,1)/SQRT(DDOT(N,G,1,G,1)*DDOT(N,W,1,W,1))
      DOT1=SQRT(DDOT(N,G,1,G,1))
      DOT2=SQRT(DDOT(N,W,1,W,1))
      OVERLAP=0.0D0
      IF (DOT1*DOT2.NE.0.0D0) OVERLAP=DDOT(N,G,1,W,1)/(DOT1*DOT2)

      IF (OVERLAP.GT.0.0D0) THEN
         IF (PTEST) PRINT*,'Search direction has positive projection onto gradient - reversing step'
         DO I=1,N
            W(ISPT+POINT*N+I)= -W(I)
         ENDDO
C        ITER=0
C        GOTO 10
      ENDIF
      
      DO I=1,N
         W(I)=G(I)
      ENDDO
      SLENGTH=0.0D0
      DO J1=1,N
         SLENGTH=SLENGTH+W(ISPT+POINT*N+J1)**2
      ENDDO
      SLENGTH=SQRT(SLENGTH)
      IF (STP*SLENGTH.GT.MAXXBFGS) STP=MAXXBFGS/SLENGTH
C
C  We now have the proposed step.
C
      DO J1=1,N
         X(J1)=X(J1)+STP*W(ISPT+POINT*N+J1)
      ENDDO 
      NDECREASE=0
20    CONTINUE
      CALL INTSECDIAG(X,COORDS,EENERGY,GNEW,ENEW,.TRUE.,RMS)

      IF ((ENEW-ENERGY)/MAX(ABS(ENEW),1.0D-100).LE.1.0D-3) THEN
C     IF ((ENEW-ENERGY)/MAX(ABS(ENEW),1.0D-100).LE.1.0D-8) THEN
         ITER=ITER+1
         ITDONE=ITDONE+1
         ENERGY=ENEW
         DO J1=1,N
            G(J1)=GNEW(J1)
         ENDDO
         IF (PTEST) WRITE(*,'(A,2G20.10,A,I6,A,F13.10)') 
     1          ' Eigenvalue and RMS force=',ENERGY,RMS,' after ',ITDONE,' steps, step:',STP*SLENGTH
      ELSE
C
C  Eigenvalue increased - try again with a smaller step size?
C
         IF (NDECREASE.GT.5) THEN
C jmc         IF (NDECREASE.GT.10) THEN
C           PRINT*,' LBFGS step cannot find a lower non-zero eigenvalue'
            NFAIL=NFAIL+1
            ITER=0  !  try resetting
            DO J1=1,N
               X(J1)=X(J1)-STP*W(ISPT+POINT*N+J1)
            ENDDO
            IF (NFAIL.GT.2) FAILED=.TRUE. ! allow one restart DJW
            GOTO 10
         ENDIF
         DO J1=1,N
            X(J1)=X(J1)-0.9*STP*W(ISPT+POINT*N+J1)
         ENDDO 
         NDECREASE=NDECREASE+1
         STP=STP/10.0D0
         IF (PTEST) WRITE(*,'(A,F19.10,A,G16.10,A,F15.8)') 
     1                  ' eigenvalue increased from ',ENERGY,' to ',ENEW,' decreasing step to ',STP*SLENGTH
         FIXIMAGE=.TRUE.
         GOTO 20
      ENDIF
C
C     Compute the new step and gradient change
C
      NPT=POINT*N
      DO I=1,N
         W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
         W(IYPT+NPT+I)= G(I)-W(I)
      ENDDO
C
C  Does it help to keep X normalised?
C
      CALL VECNORM(X,N)

      POINT=POINT+1
      IF (POINT.EQ.M) POINT=0
      FIXIMAGE=.FALSE.
      IF ((FIXAFTER.GT.0).AND.(ITDONE.GE.FIXAFTER)) FIXIMAGE=.TRUE.
      GOTO 10

      RETURN
      END
