!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999-2009 David J. Wales
!   This file is part of PATHSAMPLE.
!
!   PATHSAMPLE is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   PATHSAMPLE is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!

c     SUBROUTINE UNRESDUMP(COORDS,FNAMEF)
c     IMPLICIT NONE
c     INCLUDE 'common.h'
c     INCLUDE 'unres.h'
c     CHARACTER*15  FNAMEF 
c     INTEGER IUNIT,I,IRES
c     DOUBLE PRECISION COORDS(3*NATOMS),X(NATOMS),Y(NATOMS),Z(NATOMS)

C jmc
c     INTEGER nnt,nct,nint_gr,istart,iend,itype,maxint_gr,expon,expon2
c     PARAMETER (maxint_gr=2)
c     DOUBLE PRECISION aa,bb,augm,aad,bad,appu,bpp,ael6,ael3
c     COMMON /interact/aa(ntyp,ntyp),bb(ntyp,ntyp),augm(ntyp,ntyp),
c    & aad(ntyp,2),bad(ntyp,2),appu(2,2),bpp(2,2),ael6(2,2),ael3(2,2),
c    & expon,expon2,nnt,nct,nint_gr(maxres),istart(maxres,maxint_gr),
c    & iend(maxres,maxint_gr),itype(maxres)

c     DO I=1,NATOMS
c        X(I)=COORDS(3*(I-1)+1)
c        Y(I)=COORDS(3*(I-1)+2)
c        Z(I)=COORDS(3*(I-1)+3)
c     ENDDO

c     IUNIT = 19
c     OPEN(UNIT=IUNIT,FILE=FNAMEF,STATUS='UNKNOWN')
c       WRITE(IUNIT,'(I5)') NATOMS
c       DO I=1,NRES
c          WRITE(IUNIT,'(2I5,1X,A3,2X,A3,3F15.10)')
c    &     2*I-1,I,RESTYP(ITYPE(I)),'CA ',X(2*I-1),Y(2*I-1),Z(2*I-1)
c          WRITE(IUNIT,'(2I5,1X,A3,2X,A3,3F15.10)')
c    &     2*I,I,RESTYP(ITYPE(I)),'CB ',X(2*I),Y(2*I),Z(2*I)
c       ENDDO

c     CLOSE (IUNIT)


c     RETURN

c     END

C dae 
C get coordinates from external file 'coords'
C which is CHARMM format but use standard fortran
C reading commands, unlike CHARMM which has a limit
C on total line length and does involved procedure with strings
C
C jmc      SUBROUTINE NEWREAD(X,Y,Z,NATOMS)
      SUBROUTINE NEWREAD(X,Y,Z)
      USE COMMONS
      IMPLICIT NONE
  
      CHARACTER*4 A1
      INTEGER I
      DOUBLE PRECISION X(NATOMS),Y(NATOMS),Z(NATOMS)
 
      OPEN(UNIT=19,FILE='coords',STATUS='UNKNOWN')
      READ (19,*)
      DO I=1,NATOMS
         READ (19,'(A,3F20.10)') A1,X(I),Y(I),Z(I)
      ENDDO
  
      CLOSE(19)

      RETURN
      END

      SUBROUTINE UNRESREAD(COORDS,FNAME)
      USE COMMONS
      IMPLICIT NONE

      CHARACTER*4 A1
C     CHARACTER*72 FNAME
      CHARACTER(LEN=*) FNAME
      INTEGER I
      DOUBLE PRECISION COORDS(3*NATOMS)

      OPEN(UNIT=19,FILE=FNAME,STATUS='UNKNOWN')
      READ (19,*,END=26)
      DO I=1,NATOMS
         READ (19,*,END=26) A1,COORDS(3*(I-1)+1),COORDS(3*(I-1)+2),COORDS(3*(I-1)+3)
C        READ (19,'(A,3F20.10)') A1,X(I),Y(I),Z(I)
      ENDDO

26    CONTINUE
      CLOSE(19)

      RETURN
      END

      SUBROUTINE UNRESDUMP2(COORDS,IUNIT)
      USE COMMONS
      IMPLICIT NONE
c     CHARACTER*15  FNAMEF 
      INTEGER IUNIT,J,K
      DOUBLE PRECISION COORDS(3*NATOMS)
      DOUBLE PRECISION PEPCOORDS(3*NATOMS)

c     DO I=1,NATOMS
c        X(I)=COORDS(3*(I-1)+1)
c        Y(I)=COORDS(3*(I-1)+2)
c        Z(I)=COORDS(3*(I-1)+3)
c     ENDDO

C also writes coordinates for dummy peptide atoms (O [representing C=O],N), for
C visualisation purposes.
      DO J=1,(NATOMS/2)-1
         DO K=1,3
            PEPCOORDS(6*(J-1)+K)=(2.0D0*COORDS(6*(J-1)+K)+COORDS(6*J+K))/3.0D0
            PEPCOORDS(6*(J-1)+K+3)=(COORDS(6*(J-1)+K)+2.0D0*COORDS(6*J+K))/3.0D0
         END DO
      END DO

      WRITE(IUNIT,'(I5)') 2*NATOMS-2
      WRITE(IUNIT,'(A)') ' ' 

      DO J=1,NATOMS/2
C Calpha
         WRITE(IUNIT,'(A2,3F20.10)') 'C ',(COORDS(6*(J-1)+K),K=1,3)
C side chain
         WRITE(IUNIT,'(A2,3F20.10)') 'C ',(COORDS(6*(J-1)+K),K=4,6)
C peptide atoms
         IF (J.LT.(NATOMS/2)) THEN
             WRITE(IUNIT,'(A2,3F20.10)') 'O ',(PEPCOORDS(6*(J-1)+K),K=1,3)
             WRITE(IUNIT,'(A2,3F20.10)') 'N ',(PEPCOORDS(6*(J-1)+K),K=4,6)
         END IF
      ENDDO

c     DO I=1,NATOMS
c        WRITE(IUNIT,'(A2,3F20.10)') 'C ',X(I),Y(I),Z(I)
c     ENDDO

      RETURN

      END

      SUBROUTINE MYUNRESDUMP(COORDS,FNAMEF)
      USE COMMONS
      IMPLICIT NONE

c     CHARACTER*15  FNAMEF 
      CHARACTER(LEN=*)  FNAMEF 
      INTEGER IUNIT,J,K
      DOUBLE PRECISION COORDS(3*NATOMS)

      IUNIT=19
      OPEN(UNIT=IUNIT,FILE=FNAMEF,STATUS='UNKNOWN')
      WRITE(IUNIT,'(I5)') NATOMS

      DO J=1,NATOMS/2
C Calpha
         WRITE(IUNIT,'(A2,3F20.10)') 'C ',(COORDS(6*(J-1)+K),K=1,3)
C side chain
         WRITE(IUNIT,'(A2,3F20.10)') 'C ',(COORDS(6*(J-1)+K),K=4,6)
      END DO

      CLOSE(IUNIT)

      RETURN

      END
