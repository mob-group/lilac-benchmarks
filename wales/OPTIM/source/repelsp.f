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
      SUBROUTINE REPELSP(COORDS,RMS,INEG,MFLAG)
      USE COMMONS
      USE KEY
      IMPLICIT NONE
      DOUBLE PRECISION COORDS(3*NATOMS), DUM(3*NATOMS), DIST, RMAT(3,3)
      DOUBLE PRECISION RMS, DUMMY
      INTEGER INEG
      LOGICAL NOSHIFTSAVE, NOITSAVE, MFLAG

C     DOUBLE PRECISION REPELPUSH
C     DOUBLE PRECISION REPELTS(3*NATOMS,100)
C     INTEGER NREPELTS, REPELFROM
C     LOGICAL REPELTST, REPEL
C     COMMON /OTS/ NREPELTS, REPELTST, REPELPUSH, REPEL, REPELFROM

      INTEGER J1, J2
      SAVE NOSHIFTSAVE, NOITSAVE

C     PRINT*,'NOIT=',NOIT
C     PRINT*,'REPEL,REPELFROM,RMS=',REPEL,REPELFROM,RMS
C     PRINT*,'INEG,NOSHIFT,NOIT=',INEG,NOSHIFT,NOIT

      IF (RMS.LT.1.0D-4) THEN
         IF (REPEL.AND.(INEG.EQ.2)) THEN
            DO J2=1,NOPT
               DUM(J2)=REPELTS(J2,REPELFROM)
            ENDDO
            CALL NEWMINDIST(COORDS,DUM,NATOMS,DIST,BULKT,TWOD,ZSYM(1),.FALSE.,RIGIDBODY,DEBUG,RMAT)
            DUMMY=0.0D0
            DO J1=1,NOPT
               DUM(J1)=COORDS(J1)-DUM(J1)
               DUMMY=DUMMY+DUM(J1)**2
            ENDDO
            DUMMY=1.0D0/SQRT(DUMMY)
            DO J1=1,NOPT
               COORDS(J1)=COORDS(J1)+REPELPUSH*DUM(J1)*DUMMY
            ENDDO
            REPEL=.FALSE.
            INDEXT=.FALSE.
            HINDEX=1
            NOIT=NOITSAVE
C           NOSHIFT=NOSHIFTSAVE
         ELSE IF ((.NOT.REPEL).AND.(INEG.EQ.1)) THEN
            DO J1=1,NREPELTS
               DO J2=1,NOPT
                  DUM(J2)=REPELTS(J2,J1)
               ENDDO
               CALL NEWMINDIST(COORDS,DUM,NATOMS,DIST,BULKT,TWOD,ZSYM(1),.FALSE.,RIGIDBODY,DEBUG,RMAT)
               IF (DIST.LT.0.05D0) THEN
                  WRITE(*,'(A,F15.5,A,I5)') ' Search has moved to a distance of ',DIST,' from transition state ',J1
                  REPEL=.TRUE.
                  MFLAG=.FALSE.
                  REPELFROM=J1
                  INDEXT=.TRUE.
                  HINDEX=2
                  NOITSAVE=NOIT
                  NOIT=.TRUE.
C                 NOSHIFTSAVE=NOSHIFT
C
C     Eigenvalue shifting
C
C                 IF (NOSHIFT) THEN
C                    NOSHIFT=.FALSE.
                     DO J2=1,6
                        SHIFTL(J2)=0.0D0
                     ENDDO
                     IF (RTEST) THEN
                        IF (JZ.NE.0.0D0) THEN
                           WRITE(*,'(A,G20.10)') ' Using trans/z rotation shift=',SHIFTV
                           SHIFTL(3)=SHIFTV
                           SHIFTL(6)=SHIFTV
                        ELSE
                           WRITE(*,'(A,G20.10)') ' Using z rotation shift=',SHIFTV
                           SHIFTL(6)=SHIFTV
                        ENDIF
C                    ELSE IF (NOSHIFT) THEN
C                       WRITE(*,'(A)') 'No shifting'
                     ELSE IF (TWOD) THEN
                        IF (.NOT.BULKT) THEN
                           SHIFTL(1)=SHIFTV
                           SHIFTL(2)=SHIFTV
                           SHIFTL(6)=SHIFTV
                           WRITE(*,'(A)') ' x,y translational and z rotational shifting'
                        ELSE
                           SHIFTL(1)=SHIFTV
                        SHIFTL(2)=SHIFTV
                           WRITE(*,'(A)') ' x,y translational shifting'
                        ENDIF
                     ELSE IF (MKTRAPT) THEN
                        WRITE(*,'(A,G20.10)') ' no shifting'
                     ELSE IF (EYTRAPT.OR.(ZSYM(NATOMS).EQ.'BE')) THEN
                        SHIFTL(4)=SHIFTV
                        SHIFTL(5)=SHIFTV
                        SHIFTL(6)=SHIFTV
                        WRITE(*,'(A,G20.10)') ' Using rotational shift=',SHIFTV
                     ELSE IF (BULKT) THEN
                        WRITE(*,'(A,G20.10)') ' Using translational shift=',SHIFTV
                        DO J2=1,3
                           SHIFTL(J2)=SHIFTV
                        ENDDO
                     ELSE
                        WRITE(*,'(A,G20.10)') ' Using translational/rotational shift=',SHIFTV
                        DO J2=1,6
                           SHIFTL(J2)=SHIFTV
                        ENDDO
                     ENDIF
C                 ENDIF

                  RETURN
               ENDIF
C           DO J2=1,NOPT
C              REPELTS(J2,J1)=DUM(J2)
C           ENDDO
C           IF (DIST.LT.DMIN) THEN
C              DMIN=DIST
C              JMIN=J1
C           ENDIF
C           COMP=0.0D0
C           DO J2=1,NOPT
C              COMP=COMP+VEC(J2)*(COORDS(J2)-DUM(J2))/DIST
C           ENDDO
C           DO J2=1,NOPT
C              VEC(J2)=VEC(J2)-(1.0D0-AMOUNT)*COMP*(COORDS(J2)-DUM(J2))/DIST
C           ENDDO
C           WRITE(*,'(A,I6,F20.10)') 'J1,COMP=',J1,COMP
            ENDDO
         ENDIF
      ENDIF

      RETURN
      END

