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
      SUBROUTINE DIIS(NOPT,SVEC,ITER,QSAVE,VNEW,Q,PRMS,ENERGY,DONE)
      USE KEY
      USE GDIIS
      IMPLICIT NONE
      INTEGER J1, J2, NOPT, INDX(NDIIS+1), J, ITER
      DOUBLE PRECISION SVEC(NOPT+1,NDIIS), QD(NOPT), 
     1                 BD(NDIIS+1,NDIIS+1), CVEC(NDIIS+1), TEMP,
     2                 QSAVE(NOPT,NDIIS), D, SINV(NDIIS+1,NDIIS+1),
     3                 Q(NOPT), VNEW(NOPT), PRMS,
     4                 ESAVE(NDIIS), ENERGY
      LOGICAL SFLAG, DONE
C
C  NDIIS defines the maximum dimension
C  NDIIA is the number of steps to use
C  SVEC  contains saved steps
C  QSAVE contains saved points
C  QD    contains the extrapolated position
C  BD    is the B matrix of direct products
C
C  Dimension is so small that I might as well recompute the whole
C  B matrix each time.
C
      DO J1=1,NDIIA-1
         ESAVE(J1)=ESAVE(J1+1)
         DO J2=1,NOPT+1
            SVEC(J2,J1)=SVEC(J2,J1+1)
            QSAVE(J2,J1)=QSAVE(J2,J1+1)
         ENDDO
      ENDDO
      ESAVE(NDIIA)=ENERGY
      DO J2=1,NOPT
         SVEC(J2,NDIIA)=VNEW(J2)
         QSAVE(J2,NDIIA)=Q(J2)
      ENDDO

      IF ((MOD(ITER,NINTV).EQ.0).AND.(PRMS.LT.PCUT).AND.(ITER.GT.1).AND.(ITER.GE.NDIIA)) THEN
         WRITE(*,'(A)') ' Taking a DIIS step'
         DONE=.TRUE.
      ELSE
         RETURN
      ENDIF

      DO J1=1,NDIIA
         DO J2=J1,NDIIA
            TEMP=0.0D0
C           DO J3=1,NOPT
C              TEMP=TEMP+SVEC(J3,J2)*SVEC(J3,J1)
C           ENDDO
            TEMP=(ESAVE(J1)-ESAVE(NDIIA)+(ESAVE(NDIIA-1)-ESAVE(NDIIA)))*
     1           (ESAVE(J2)-ESAVE(NDIIA)+(ESAVE(NDIIA-1)-ESAVE(NDIIA)))/(ESAVE(NDIIA-1)-ESAVE(NDIIA))**2
            BD(J2,J1)=TEMP
            BD(J1,J2)=TEMP
         ENDDO
      ENDDO

      DO J1=1,NDIIA+1
         BD(NDIIA+1,J1)=1.0D0
         BD(J1,NDIIA+1)=1.0D0
      ENDDO
      BD(NDIIA+1,NDIIA+1)=0.0D0

       PRINT*,'B matrix:'
       DO J1=1,NDIIA+1
          WRITE(*,98) (BD(J1,J2),J2=1,NDIIA+1)
 98       FORMAT(10G13.6)
       ENDDO

C      PRINT*,'QSAVE matrix:'
C      DO J1=1,NOPT
C         WRITE(*,96) (QSAVE(J1,J2),J2=1,NDIIA)
C96       FORMAT(10F13.7)
C      ENDDO

      DO J1=1,NDIIA+1
         DO 52 J2=1,NDIIA+1
            SINV(J1,J2)=0.0D0
52       CONTINUE
         SINV(J1,J1)=1.0D0
      ENDDO
      SFLAG=.FALSE.
      CALL LUDCMP(BD,NDIIA+1,NDIIS+1,INDX,D,SFLAG)
      IF (SFLAG) PRINT*,'Singular matrix'
      DO J=1,NDIIA+1
         D=D*BD(J,J)
      ENDDO
      WRITE(*,'(A,G15.8)') ' DIIS determinant=',D
      DO J1=1,NDIIA+1
         CALL LUBKSB(BD,NDIIA+1,NDIIS+1,INDX,SINV(1,J1))
      ENDDO
C
C  SINV is symmetric, so no problem here.
C
      DO J1=1,NDIIA
         CVEC(J1)=SINV(J1,NDIIA+1)
      ENDDO
      WRITE(*,'(A)') ' Vector of DIIS expansion coefficients:'
      WRITE(*,93) (CVEC(J1),J1=1,NDIIA)
93    FORMAT(6(2X,E12.5))
C     PRINT*,'DIIS residual=',CVEC(NDIIA+1)

      DO J1=1,NOPT
         QD(J1)=0.0D0
      ENDDO
      DO J1=1,NDIIA
         DO J2=1,NOPT
            QD(J2)=QD(J2)+CVEC(J1)*QSAVE(J2,J1)
         ENDDO
      ENDDO
      DO J1=1,NOPT
         Q(J1)=QD(J1)
      ENDDO

      RETURN
      END
