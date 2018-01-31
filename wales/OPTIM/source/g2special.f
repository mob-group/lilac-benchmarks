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
      SUBROUTINE G2SPECIAL(NATOMS,XOLD,GRAD,G2GRAD,G2ENERGY,G2RMS,EREAL,RMSREAL,MFLAG)
      USE MODHESS
      USE KEY
      use porfuncs
      IMPLICIT NONE 
      INTEGER NATOMS, INFO, J1, JMIN, JZERO, NTRANS, J2, NCOUNT, NINFO
      LOGICAL MFLAG
      DOUBLE PRECISION GRAD(3*NATOMS), TEMPA(9*NATOMS), EVMIN, ZERO, STEP,
     1                 G2GRAD(3*NATOMS), G2ENERGY, DDOT, G2RMSP, SCALE, EREAL, RMSREAL, XNEW(3*NATOMS),
     2                 SVEC(3*NATOMS), GOVER, GMOD, GSTEP, SIXC, Q1MQ0, ZEROOLD, ZEROOLD2, G2RMS, GSTHRESH,
     3                 XOLD(3*NATOMS), XOLD2(3*NATOMS)
      DOUBLE PRECISION DIAG(3*NATOMS)
      INTEGER NSPECIAL, NALLOW, ISTAT
      COMMON /G2/ GSTHRESH, SCALE, NSPECIAL, NALLOW, NINFO

      MFLAG=.FALSE.
      G2RMSP=G2RMS

      CALL DSYEV('V','U',3*NATOMS,HESS,3*NATOMS,DIAG,TEMPA,9*NATOMS,INFO)
      IF (INFO.NE.0) PRINT*,'WARNING - INFO=',INFO,' in DSYEV'
c     call eigensort_val_asc(diag,hess,3*natoms,3*natoms)

      EVMIN=1.0D6
      ZERO=1.0D6
      NTRANS=0
      DO J1=1,3*NATOMS
         IF ((ABS(DIAG(J1)).LT.1.0D-10).AND.(NTRANS.LT.3)) THEN
            NTRANS=NTRANS+1
         ELSE
            IF (DIAG(J1).LT.EVMIN) THEN
               EVMIN=DIAG(J1)
               JMIN=J1
            ENDIF
            IF (ABS(DIAG(J1)).LT.ZERO) THEN
               ZERO=ABS(DIAG(J1))
               JZERO=J1
            ENDIF
         ENDIF
      ENDDO
      ZERO=DIAG(JZERO)
      ZEROOLD=ZERO
      ZEROOLD2=ZERO
C     WRITE(*,'(6F20.10)') (DIAG(J1),J1=1,3*NATOMS)
C     WRITE(*,'(A,i4)') 'Number of translations=',NTRANS
      WRITE(*,'(A,2F20.10)') 'Lowest eigenvalue and closest to zero= ',EVMIN,ZERO
      WRITE(*,'(A,4F20.10)') 'G^2, RMS force and real energy and RMS=',G2ENERGY,G2RMS,EREAL,RMSREAL
      GMOD=SQRT(G2ENERGY)
      GOVER=DDOT(3*NATOMS,GRAD,1,HESS(1,JZERO),1)/GMOD
      WRITE(*,'(A,F20.10)') 'Overlap between normalised gradient vector and hessian eigenvector=',GOVER

      IF (ABS(ABS(GOVER)-1.0D0).LT.1.0D-2) THEN

         WRITE(*,'(A)') 'Trying special steps for solution with an extra zero eigenvalue'
         FIXIMAGE=.TRUE.

         WRITE(*,'(A,I3,A,4F20.10)') 'At step    ',0,' zero, G2energy, G2RMS=',ZERO,G2ENERGY,G2RMS
         NCOUNT=0

10       NCOUNT=NCOUNT+1
C
C  Combined tangent and gradient step. If this is the first step, then no move
C  along the gradient.
C
         DO J2=1,3*NATOMS
            SVEC(J2)=0.0D0
         ENDDO
         DO J1=1,3*NATOMS
            IF ((J1.NE.JZERO).AND.(ABS(DIAG(J1)).GT.1.0D-10)) THEN
               STEP=0.0D0
               DO J2=1,3*NATOMS
                  STEP=STEP+G2GRAD(J2)*HESS(J2,J1)
               ENDDO
               STEP=STEP/(2*DIAG(J1)**2)
C              WRITE(*,'(A,I4,2F20.10)') 'J1,DIAG,STEP=',J1,DIAG(J1),STEP
               DO J2=1,3*NATOMS
                  SVEC(J2)=SVEC(J2)+STEP*HESS(J2,J1)
               ENDDO
            ENDIF
         ENDDO

30       CONTINUE 
C
C  Step along the gradient. 
C
         IF (NCOUNT.GT.1) THEN
            GSTEP=0.0D0
            DO J1=1,3*NATOMS
               GSTEP=GSTEP+(XOLD(J1)-XOLD2(J1))*GRAD(J1)/GMOD
            ENDDO
            SIXC=(ZEROOLD-ZEROOLD2)/GSTEP
            Q1MQ0=ZEROOLD2/SIXC

            WRITE(*,'(A,4G20.10)') '6*C, Q1-Q0 a, Q1-Q0 b, STEP=',SIXC,Q1MQ0,(ZEROOLD-SIXC*GSTEP)/SIXC,-Q1MQ0
         ELSE
            Q1MQ0=0.0D0
         ENDIF

         DO J2=1,3*NATOMS
            XNEW(J2)=XOLD(J2)-(SVEC(J2)+GRAD(J2)*Q1MQ0/GMOD)*SCALE
         ENDDO
C
C  Calculate properties for the new geometry.
C
         CALL POTENTIAL(XNEW,EREAL,GRAD,.TRUE.,.TRUE.,RMSREAL,.FALSE.,.FALSE.)
         CALL DSYMV('U',3*NATOMS,2.0D0,HESS,3*NATOMS,GRAD,1,0.0D0,G2GRAD,1)
         G2RMS= DSQRT(DDOT(3*NATOMS,G2GRAD,1,G2GRAD,1)/(3*NATOMS))
         CALL DSYEV('V','U',3*NATOMS,HESS,3*NATOMS,DIAG,TEMPA,9*NATOMS,INFO)
         IF (INFO.NE.0) PRINT*,'WARNING - INFO=',INFO,' in DSYEV'

         call eigensort_val_asc(diag,hess,3*natoms,3*natoms)
         WRITE(*,'(6F20.10)') (DIAG(J1),J1=1,3*NATOMS)

         G2ENERGY=DDOT(3*NATOMS,GRAD,1,GRAD,1)
         GMOD=SQRT(G2ENERGY)

         ZERO=1.0D6
         NTRANS=0
         DO J1=1,3*NATOMS
            IF ((ABS(DIAG(J1)).LT.1.0D-10).AND.(NTRANS.LT.3)) THEN
               NTRANS=NTRANS+1
            ELSE IF (ABS(DIAG(J1)).LT.ZERO) THEN
               ZERO=ABS(DIAG(J1))
               JZERO=J1
            ENDIF
         ENDDO
         ZERO=DIAG(JZERO)
         WRITE(*,'(A,I3,A,4F20.10)') 'After step ',NCOUNT,' zero, G2energy, G2RMS=',ZERO,G2ENERGY,G2RMS,SCALE
         CALL FLUSH(6)
         IF (G2RMS.GT.G2RMSP) THEN
            IF (SCALE.LT.0.000001D0) THEN
               WRITE(*,'(A,F15.5,A)') 'G2 RMS increased for step of ',SCALE,' return'
               RETURN
            ELSE
C              WRITE(*,'(A,F15.5,A)') 'G2 RMS increased for step of ',SCALE,' try half the step'
               SCALE=SCALE/2.0D0
               GOTO 30
            ENDIF
         ELSE
            SCALE=MIN(1.0D0,SCALE*1.1D0)
         ENDIF
         IF (G2RMS.LT.GMAX) THEN
            MFLAG=.TRUE.
            FIXIMAGE=.FALSE.
            RETURN
         ENDIF
         DO J1=1,3*NATOMS
            XOLD2(J1)=XOLD(J1)
            XOLD(J1)=XNEW(J1)
         ENDDO
         ZEROOLD2=ZEROOLD
         ZEROOLD=ZERO
         G2RMSP=G2RMS
         
         PRINT*

         IF (NCOUNT.LT.NALLOW) GOTO 10

      ELSE IF ((G2ENERGY.LT.1.0D-1).AND.(RMSREAL.LT.1.0D-1)) THEN

         WRITE(*,'(A)') 'Trying special steps for solution at a true stationary point'
         FIXIMAGE=.TRUE.

20       NCOUNT=NCOUNT+1

         DO J2=1,3*NATOMS
            SVEC(J2)=0.0D0
         ENDDO
         DO J1=1,3*NATOMS
            IF (ABS(DIAG(J1)).GT.1.0D-10) THEN
               STEP=0.0D0
               DO J2=1,3*NATOMS
                  STEP=STEP-GRAD(J2)*HESS(J2,J1)
               ENDDO
               STEP=STEP/DIAG(J1)
C              WRITE(*,'(A,I4,2F20.10)') 'J1,DIAG,STEP=',J1,DIAG2(J1),STEP
               DO J2=1,3*NATOMS
                  SVEC(J2)=SVEC(J2)+STEP*HESS(J2,J1)
               ENDDO
            ENDIF
         ENDDO

40       DO J2=1,3*NATOMS
            XOLD(J2)=XOLD(J2)+SVEC(J2)*SCALE
         ENDDO

         CALL POTENTIAL(XOLD,EREAL,GRAD,.TRUE.,.TRUE.,RMSREAL,.FALSE.,.FALSE.)
         CALL DSYMV('U',3*NATOMS,2.0D0,HESS,3*NATOMS,GRAD,1,0.0D0,G2GRAD,1)
         G2RMS= DSQRT(DDOT(3*NATOMS,G2GRAD,1,G2GRAD,1)/(3*NATOMS))
         CALL DSYEV('V','U',3*NATOMS,HESS,3*NATOMS,DIAG,TEMPA,9*NATOMS,INFO)
         IF (INFO.NE.0) PRINT*,'WARNING - INFO=',INFO,' in DSYEV'

         call eigensort_val_asc(diag,hess,3*natoms,3*natoms)
         WRITE(*,'(6F20.10)') (DIAG(J1),J1=1,3*NATOMS)

         G2ENERGY=DDOT(3*NATOMS,GRAD,1,GRAD,1)
         GMOD=SQRT(G2ENERGY)

         ZERO=1.0D6
         NTRANS=0
         DO J1=1,3*NATOMS
            IF ((ABS(DIAG(J1)).LT.1.0D-10).AND.(NTRANS.LT.3)) THEN
               NTRANS=NTRANS+1
            ELSE IF (ABS(DIAG(J1)).LT.ZERO) THEN
               ZERO=ABS(DIAG(J1))
               JZERO=J1
            ENDIF
         ENDDO
         ZERO=DIAG(JZERO)
         WRITE(*,'(A,I3,A,4F20.10)') 'After NR/EF   step ',NCOUNT,' zero, G2energy, G2RMS=',ZERO,G2ENERGY,G2RMS
         IF (G2RMS.GT.G2RMSP) THEN
            DO J2=1,3*NATOMS
               XOLD(J2)=XOLD(J2)-SVEC(J2)*SCALE
            ENDDO
            IF (SCALE.LT.0.000001D0) THEN
               WRITE(*,'(A,F15.5,A)') 'G2 RMS increased for step of ',SCALE,' return'
               RETURN
            ELSE
C              WRITE(*,'(A,F15.5,A)') 'G2 RMS increased for step of ',SCALE,' try half the step'
               SCALE=SCALE/2.0D0
               GOTO 40
            ENDIF
         ELSE
            SCALE=SCALE*1.1D0
         ENDIF
         IF (G2RMS.LT.GMAX) THEN
            MFLAG=.TRUE.
            FIXIMAGE=.FALSE.
            RETURN
         ENDIF
         
         PRINT*

         G2RMSP=G2RMS
         IF (NCOUNT.LT.NALLOW) GOTO 20

      ENDIF

      FIXIMAGE=.FALSE.

      RETURN
      END
