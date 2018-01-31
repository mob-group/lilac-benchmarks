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

      SUBROUTINE CHARMMDUMP(COORDS,FNAME)
      USE COMMONS
      IMPLICIT NONE
      CHARACTER(LEN=4) SID
      CHARACTER(LEN=*)  FNAME
      INTEGER IUNIT,I,RID
      DOUBLE PRECISION COORDS(3*NATOMS)

C     CHARACTER(LEN=4) RESLABEL(NATOMS),ATOMLABEL(NATOMS)
C     INTEGER RESNUMBER(NATOMS)
C     COMMON /CHARMMLABEL/ RESLABEL, ATOMLABEL, RESNUMBER

      iunit = 19
      if (machine) then
           OPEN(IUNIT,FILE=FNAME,access='direct',form='unformatted',status='unknown',recl=3*8*Natoms)
           write(IUNIT,rec=1) (COORDS(i),i=1,3*Natoms)
      else
           SID=' 1'
           RID=0 ! was uninitialised!
           OPEN(UNIT=IUNIT,FILE=FNAME,STATUS='UNKNOWN')
           WRITE(IUNIT,'(I5)') NATOMS
           DO I=1,NATOMS
              WRITE(IUNIT,'(2I5,1X,A4,1X,A4,3F20.12,2X,A2,2X,I2)')
     &        I,RESNUMBER(I),RESLABEL(I),ATOMLABEL(I),
     &        COORDS(3*(I-1)+1),COORDS(3*(I-1)+2),COORDS(3*(I-1)+3)
           ENDDO
      endif
      close (iunit)
      end

C dae 
C get coordinates from external file 'coords'
C which is CHARMM format but use standard fortran
C reading commands, unlike CHARMM which has a limit
C on total line length and does involved procedure with strings
C
      SUBROUTINE READREF(FNAME)
      use utils
      USE COMMONS
      IMPLICIT NONE
  
      CHARACTER(LEN=*) FNAME
      INTEGER I,I1
 
C     CHARACTER(LEN=4) RESLABEL(NATOMS),ATOMLABEL(NATOMS)
C     INTEGER RESNUMBER(NATOMS)
C     COMMON /CHARMMLABEL/ RESLABEL, ATOMLABEL, RESNUMBER

      OPEN(UNIT=19,FILE=FNAME,STATUS='OLD',iostat=I)
      call openiostat(I,fname,'readref')
      READ (19,*)
      DO I=1,NATOMS
         READ (19,*) I1,RESNUMBER(I),RESLABEL(I),ATOMLABEL(I)
      ENDDO
  
      CLOSE(19)

      RETURN
      END


C get coordinates from external file 
C which is CHARMM format but use standard fortran
C reading commands, unlike CHARMM which has a limit
C on total line length and does involved procedure with strings
C
      SUBROUTINE CHARMMREAD(COORDS,FNAME)
      USE COMMONS
      IMPLICIT NONE

      CHARACTER(LEN=4) A1
      CHARACTER(LEN=*)  FNAME
      INTEGER I,I1
      DOUBLE PRECISION COORDS(3*NATOMS)

      if (machine) then
           OPEN(19,FILE=FNAME,access='direct',form='unformatted',status='old',recl=3*8*Natoms,action='read')
           read(19,rec=1) (COORDS(i),i=1,3*Natoms)
      else
           OPEN(UNIT=19,FILE=FNAME,STATUS='UNKNOWN')
           READ (19,*)
           DO I=1,NATOMS
              READ (19,*) I1,I1,A1,A1,COORDS(3*(I-1)+1),COORDS(3*(I-1)+2),COORDS(3*(I-1)+3),I1,I1
           ENDDO
      endif
      CLOSE(19)

      RETURN
      END

