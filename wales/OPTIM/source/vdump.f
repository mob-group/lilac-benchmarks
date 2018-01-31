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
      SUBROUTINE VDUMP(DIAG,ZT,N,M)
      USE KEY
      USE MODHESS
      USE MODCHARMM, ONLY: CHRMMT
      USE PORFUNCS
      IMPLICIT NONE
      INTEGER M, N, J1, J2, ISTAT, MCOUNT
      DOUBLE PRECISION DIAG(M)
      LOGICAL ZT(M)
  
C
C  dump the eigenvectors which correspond to non-zero eigenvalues
C  in file vectors.dump
C
      IF (.NOT.ALLSTEPS) REWIND(44)
      IF (ALLVECTORS) THEN
!       csw34> Find how many modes are going to be printed and echo this
!       to a file for the analysis script (this could be removed when all the
!       dumping is done in OPTIM). This is a problem with FREEZE where
!       some non-real modes have zero eiganvalue.
            MCOUNT=0
            DO J1=1,N
                IF (ZT(J1)) MCOUNT=MCOUNT+1
            ENDDO
            OPEN(UNIT=499,FILE='nmodes.dat',STATUS='UNKNOWN')
            WRITE(499,'(I6)') MCOUNT
            CLOSE(499)
         DO J1=1,N
            IF (ZT(J1)) THEN

! If printing the mass weighted vectors (normal modes), convert omega^2
! into the vibrational frequency in the specified unit system using FRQCONV.
! Normally, this will be either internal units or rad/s, for compatibility with PATHSAMPLE.
! Other unit systems can be specified using the FRQCONV keyword.
! NOTE: This behaviour has changed as of 23/9/16. Until now, the frequencies were always
! multiplied by 108.52, which is the conversion factor from kCal mol^-1 and Angstrom units
! to cm^-1 frequencies. However, this is obviously not appropriate for all systems. If
! you wish to retrieve the old behaviour, simply add FRQCONV 108.52 to your odata file.
! But note that the square frequencies used for the log product in path.info will then be given
! in units of cm^-2
              IF (MWVECTORS) THEN
                        WRITE(44,'(F20.10)') DSQRT(DIAG(J1))*FRQCONV
              ELSE
                        WRITE(44,'(F20.10)') DIAG(J1)
              ENDIF
!              IF (ANGLEAXIS2) N=N/2
               WRITE(44,'(3F20.10)') (HESS(J2,J1),J2=1,N)
!              IF (ANGLEAXIS2) N=N*2
            ENDIF
         ENDDO
      ELSE
         DO J1=N,1,-1
            IF (ZT(J1)) THEN
! As above
               IF (MWVECTORS) THEN
                        WRITE(44,'(F20.10)') DSQRT(DIAG(J1))*FRQCONV
               ELSE
                        WRITE(44,'(F20.10)') DIAG(J1)
               ENDIF
               WRITE(44,'(3F20.10)') (HESS(J2,J1),J2=1,N)
               CALL FLUSH(44)
               RETURN
            ENDIF
         ENDDO
      ENDIF
      CALL FLUSH(44)
      RETURN
      END

