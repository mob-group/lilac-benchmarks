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
      SUBROUTINE SIAZ(VEC,ROT,IAXIS)
      IMPLICIT NONE
      DOUBLE PRECISION VEC(3),ROT(3,3),RM(3,3),TEMP(3,3),DIS,ARG1,AZIM,DISP,ARG2,ANG1,DOTOPT
      INTEGER J1,J2,IAXIS
C     Write (6,*) 'SIAZ: vec, iaxis'
C     Write (6,'(a,3f12.6,i3)') 'SIAZ: ',vec, iaxis
C 
C DETERMINE LENGTH OF VECTOR.
C
      ROT(1,1)=1.0D0
      ROT(2,2)=1.0D0
      ROT(3,3)=1.0D0
      ROT(1,2)=0.0D0
      ROT(1,3)=0.0D0
      ROT(2,1)=0.0D0
      ROT(2,3)=0.0D0
      ROT(3,1)=0.0D0
      ROT(3,2)=0.0D0
      DIS=DSQRT(DOTOPT(VEC,VEC,3))
      IF (DABS(DIS).LT.1.0D-10) RETURN
C     Write (6,'(a,f12.6)') 'SIAZ: dis', dis
C
C   RETURNS ROTATION MATRIX NEEDED TO PLACE VECTOR
C   VEC ALONG THE X (IAXIS=1), Y (IAXIS=2) OR Z (IAXIS=3)
C   AXIS.  TAIL OF VECTOR IMPLICITLY ASSUMED TO BE AT 
C   THE ORIGIN.
C
      IF(IAXIS.EQ.3)THEN
        ARG1=VEC(3)/DIS
        IF(DABS(ARG1).GT.1.D0)ARG1=ARG1/DABS(ARG1)
        AZIM=-1.D0*DACOS(ARG1)
        DISP=DSQRT(DOTOPT(VEC,VEC,2))
        IF (DABS(DISP).LT.1.0D-10) RETURN
        ARG2=VEC(1)/DISP
        IF(DABS(ARG2).GT.1.D0)ARG2=ARG2/DABS(ARG2)
        ANG1=DACOS(ARG2)
        ANG1=SIGN(ANG1,-1.D0*VEC(2))
        CALL ROTM(3,ANG1,0,RM)
        CALL ROTM(2,AZIM,0,ROT)
C       CALL MATMUL(TEMP,RM,ROT,3,3,3,3,3,3)
        CALL MUL3(TEMP,RM,ROT)
        DO 17 J1=1,3
           DO 18 J2=1,3
              ROT(J2,J1)=TEMP(J2,J1)
18         CONTINUE
17      CONTINUE
      ELSEIF(IAXIS.EQ.1)THEN
        ARG1=VEC(1)/DIS
        IF(DABS(ARG1).GT.1.D0)ARG1=ARG1/DABS(ARG1)
        AZIM=-1.D0*DACOS(ARG1)
        DISP=DSQRT(DOTOPT(VEC(2),VEC(2),2))
        IF (DABS(DISP).LT.1.0D-10) RETURN
        ARG2=VEC(2)/DISP
        IF(DABS(ARG2).GT.1.D0)ARG2=ARG2/DABS(ARG2)
        ANG1=DACOS(ARG2)
        ANG1=SIGN(ANG1,-1.D0*VEC(3))
        CALL ROTM(1,ANG1,0,RM)
        CALL ROTM(3,AZIM,0,ROT)
C       CALL MATMUL(TEMP,RM,ROT,3,3,3,3,3,3)
        CALL MUL3(TEMP,RM,ROT)
        DO 27 J1=1,3
           DO 28 J2=1,3
              ROT(J2,J1)=TEMP(J2,J1)
28         CONTINUE
27      CONTINUE
      ELSEIF(IAXIS.EQ.2)THEN
        ARG1=VEC(2)/DIS
        IF(DABS(ARG1).GT.1.D0)ARG1=ARG1/DABS(ARG1)
        AZIM=-1.D0*DACOS(ARG1)
        DISP=DSQRT(DIS**2-VEC(2)**2)
        IF (DABS(DISP).LT.1.0D-10) RETURN
        ARG2=VEC(3)/DISP
        IF(DABS(ARG2).GT.1.D0)ARG2=ARG2/DABS(ARG2)
        ANG1=DACOS(ARG2)
        ANG1=SIGN(ANG1,-1.D0*VEC(1))
        CALL ROTM(2,ANG1,0,RM)
        CALL ROTM(1,AZIM,0,ROT)
C       CALL MATMUL(TEMP,RM,ROT,3,3,3,3,3,3)
        CALL MUL3(TEMP,RM,ROT)
        DO 37 J1=1,3
           DO 38 J2=1,3
              ROT(J2,J1)=TEMP(J2,J1)
38         CONTINUE
37      CONTINUE
       ENDIF
       RETURN
       END

      SUBROUTINE MUL3(TEMP,RM,ROT)
      IMPLICIT NONE
 
      DOUBLE PRECISION TEMP(3,3), RM(3,3), ROT(3,3)

      TEMP(1,1)=RM(1,1)*ROT(1,1)+RM(1,2)*ROT(2,1)+RM(1,3)*ROT(3,1)
      TEMP(2,1)=RM(2,1)*ROT(1,1)+RM(2,2)*ROT(2,1)+RM(2,3)*ROT(3,1)
      TEMP(3,1)=RM(3,1)*ROT(1,1)+RM(3,2)*ROT(2,1)+RM(3,3)*ROT(3,1)
      TEMP(1,2)=RM(1,1)*ROT(1,2)+RM(1,2)*ROT(2,2)+RM(1,3)*ROT(3,2)
      TEMP(2,2)=RM(2,1)*ROT(1,2)+RM(2,2)*ROT(2,2)+RM(2,3)*ROT(3,2)
      TEMP(3,2)=RM(3,1)*ROT(1,2)+RM(3,2)*ROT(2,2)+RM(3,3)*ROT(3,2)
      TEMP(1,3)=RM(1,1)*ROT(1,3)+RM(1,2)*ROT(2,3)+RM(1,3)*ROT(3,3)
      TEMP(2,3)=RM(2,1)*ROT(1,3)+RM(2,2)*ROT(2,3)+RM(2,3)*ROT(3,3)
      TEMP(3,3)=RM(3,1)*ROT(1,3)+RM(3,2)*ROT(2,3)+RM(3,3)*ROT(3,3)

      RETURN
      END
      
