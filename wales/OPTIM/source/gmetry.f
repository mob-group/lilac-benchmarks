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
C     ROUTINE TO CALCULATE CARTESIAN COORDINATES FROM INTERNAL
C     COORDINATE REPRESENTATION.  SOME OF THIS HAS BEEN LIFTED FROM
C     PRDDO, ALTHOUGH SOME IMPROVEMENTS HAVE BEEN MADE.
C     CONNECTIVITY OF FIRST THREE MUST BE 1-2-3 IN
C     INTERNAL COORDINATE REP.
C
      SUBROUTINE GMETRY(ITER,VEC,Q)
      USE KEY
      USE COMMONS
      IMPLICIT NONE
      INTEGER ITER, NUMAT, NDIS, I, J1, J2, J, NCOUNT
      DOUBLE PRECISION VEC(3*NATOMS), Q(3*NATOMS), SCALE, TX, TY, TZ
C
      NUMAT=0
      NDIS=1
      SCALE=1.0D0
      IF (IPRNT .GE. 20) then
         WRITE (*, *) 'Cartesian coordinates before scaling.'
         WRITE (*, '(I3,3f12.6)') (i,Q(3*i-2),Q(3*i-1),Q(3*i),
     1      i=1,NATOMS)
      ENDIF
C     IF (NDIS.NE.0) SCALE=0.5291772D0
      IF (NDIS.NE.0) SCALE=1.0D0
      DO I=1,NATOMS
         Q(3*I-2)=Q(3*I-2)/SCALE
         Q(3*I-1)=Q(3*I-1)/SCALE
         Q(3*I)=Q(3*I)/SCALE
      ENDDO
      IF (IPRNT .GE. 20) then
         WRITE (*, *) 'Cartesian coordinates after scaling.'
         WRITE (*, '(I3,3F12.6)') (I,Q(3*i-2),Q(3*i-1),Q(3*i),
     1      i=1,NATOMS)
      ENDIF
C
C Permute the coordinates if necessary - only for FH
C
      IF (ZSYM(NATOMS)(1:2).EQ.'FH') THEN
      NCOUNT=0
130   CONTINUE
      DO 150 J1=1,NATOMS-1
         TX=Q(3*(J1-1)+1)-Q(3*J1+1)
         TY=Q(3*(J1-1)+2)-Q(3*J1+2)
         TZ=Q(3*(J1-1)+3)-Q(3*J1+3)
         IF (DABS(TZ*TZ+TY*TY).LT.1.0D-10) THEN
            NCOUNT=NCOUNT+1
            IF (NCOUNT.GT.2) THEN
               PRINT*,'*** WARNING ***'
               PRINT*,'*** Cannot permute the coordinates in gmetry'
               GOTO 160
c              STOP
            ENDIF
            DO 140 J2=1,NATOMS
               TX=Q(3*(J2-1)+1)
               TY=Q(3*(J2-1)+2)
               TZ=Q(3*(J2-1)+3)
               Q(3*(J2-1)+1)=TY
               Q(3*(J2-1)+2)=TZ
               Q(3*(J2-1)+3)=TX
140         CONTINUE
            GOTO 130
         ENDIF
150   CONTINUE
      ENDIF

160   IF (NATOMS.EQ.NUMAT) RETURN
      J=0
      DO 170 I=1,NATOMS
         IF (NR(I) .NE. 99) THEN
            J=J+1
C
C  Turn off the axis permutations -
C  this is the only way I can find of preserving the
C  original Cartesian coordinates.
C
            Q(3*J-2)=Q(3*I-2)
            Q(3*J-1)=Q(3*I-1)
            Q(3*J)=Q(3*I)
         ENDIF
170   CONTINUE
      If (IPRNT .GE. 20) then
         WRITE (*,*) 'Cartesian coordinates from GMETRY.'
C
         WRITE (*,'(I3,3F12.6)') (i,Q(3*i-2),Q(3*i-1),Q(3*i),i=1,NATOMS)
      ENDIF
C
C Check the theta angles if we are doing water.
C
      IF ((ZSYM(NATOMS)(1:1).EQ.'W').AND.(.NOT.ANGLEAXIS)) CALL TCHECK(VEC,Q)
      RETURN
      END
