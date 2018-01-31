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
C     GIVES ATOMIC NUMBER AND ATOMIC MASS INFORMATION
C
      SUBROUTINE PERTABLE
      USE COMMONS
      USE KEY, ONLY: NTIPT, PAHAT, STOCKAAT, NCARBON, QUIPATOMMASS

      IMPLICIT NONE
C
C     Setup the periodic table - we have entered MXTABL elements here
C
      INTEGER MXTABL,I,J
      PARAMETER (MXTABL = 110)
      DOUBLE PRECISION ATMSS(MXTABL)
      CHARACTER(LEN=2) ATSYM(MXTABL)
C
      DATA ATSYM /'H ','HE','LI','BE','B ','C ','N ','O ','F ','NE',
     &     'NA','MG','AL','SI','P ','S ','CL','AR','K ','CA',
     &     'SC','TI','V ','CR','MN','FE','CO','NI','CU','ZN',
     &     'GA','GE','AS','SE','BR','KR','AX','M1','M2','M3','M4',
     &     'AU','AG','NI','FH','W2','W4','W3','W1','ME','JM','SC',
     &     'M ','JC','CC','TB','C6','OR','P6',
     &     'W5','TE','GV','MS','MP','MV','AZ','PL','IN','LP','TT',
     &     'PL','MI','LS','GL','SW','SM','CP','DZ','LC','LK','PR',
     &     'BC','Z1','Z2','CD','C1','CK','LM','BL','TH','GO','SY',
     &     'ZF','DS','GP','AM','DB', 'ST', 'PT', 'PD', 'AK', 'DF', 
     &     'S4', 'DT', 'P4', 'AO', 'UP', 'G1', 'G2', 'RM' /
      DATA ATMSS / 1.007825D+00 , 4.00260D+00 , 7.01600D+00 ,
     &     9.01218D+00  , 11.00931D+00 , 12.00000D+00 ,
     &     14.00307D+00 , 15.99491D+00 , 18.99840D+00 ,
     &     19.99244D+00 , 22.989800D+00, 23.98504D+00 ,
     &     26.98153D+00 , 27.97693D+00 , 30.97376D+00 ,
     &     31.97207D+00 , 34.95885D+00 , 39.94800D+00  ,
     &     39.09830D+00 , 40.08000D+00 , 44.95592D+00  ,
     &     47.90000D+00 , 50.49200D+00 , 51.99600D+00  ,
     &     54.93810D+00 , 55.84700D+00 , 58.93320D+00  ,
     &     58.71000D+00 , 63.54000D+00 , 63.57000D+00  ,
     &     69.72000D+00 , 72.59000D+00 , 74.92160D+00  ,
     &     78.96000D+00 , 79.90900D+00 , 83.80000D+00  ,
     &     1.000000D+00 , 39.94800D+00 ,  1.00000D+00  ,
     &     39.94800D+00 , 39.94800D+00 ,196.96665D0    ,
     &    107.868D0     , 58.69D0      , 39.94800D0    ,
     &     18.0D0       , 18.0D0       , 18.0D0        ,
     &     18.0D0       , 39.94800D+00 , 27.97693D+00  ,
     &     39.948D0     , 39.948D0     , 39.948D0      ,
     &     12.0D0       ,  1.0D0       , 1.0D0         ,
     &      1.0D0       ,  1.0D0       , 1.0D0         ,
     &     18.0D0       , 127.60D0     , 1.0D0         ,
     &      1.0D0       ,  1.0D0       , 1.0D0         ,
     &      1.0D0       ,  1.0D0       , 1.0D0         ,
     &      1.0D0       ,  1.0D0       , 1.0D0         ,
     &      1.0D0       ,  1.0D0       , 1.0D0         ,
     &      1.0D0       ,  1.0D0       , 1.0D0         ,
     &      1.0D0       ,  1.0D0       , 1.0D0         ,
     &      1.0D0       ,  1.0D0       , 1.0D0         ,
     &      1.0D0       ,  1.0D0       , 1.0D0         ,
     &      1.0D0       ,  1.0D0       , 1.0D0         ,
     &      1.0D0       ,  1.0D0       , 1.0D0         ,
     &      1.0D0       ,  1.0D0       , 1.0D0         ,
     &      1.D0        ,  1.0D0       , 195.08D0      ,
     &    106.42D0      ,  1.0D0       , 1.0D0         ,
     &      1.0D0       ,  1.D0        , 1.D0          ,
     &      1.0D0       ,  1.0D0       , 1.0D0         ,
     &      1.0D0       , 1.0D0 /

C 
C     dg413: added a new symbol for userpotential called UP
C            at the moment its mass is set to 1.D0     
C
C     For each entry in the Z-matrix, get a mass
C
      IF (IPRNT .GE. 10) WRITE (*,9000)
 9000 FORMAT ('PERTABLE: Mass and at. nr. lookup'/
     1     'Line Symbol AtNr   At. Mass')
      DO 10 I = 1, NATOMS
         IF (ZSYM(I) .EQ. 'X ') THEN
            IATNUM(I) = 0
            ATMASS(I) = 0.0
         ELSE
            DO 20 J = 1, MXTABL
               IF (ZSYM(I) .EQ. ATSYM(J)) THEN
                  IF (ZSYM(I)(1:1).EQ.'W') THEN ! to generalise this to other rigid 
                     IF (I.LE.NATOMS/2) THEN
                        IATNUM(3*(I-1)+1)=8        ! bodies will require site/mass specification !
                        IATNUM(3*(I-1)+2)=1
                        IATNUM(3*(I-1)+3)=1
                        ATMASS(3*(I-1)+1)=16.0D0
                        ATMASS(3*(I-1)+2)=1.0D0
                        ATMASS(3*(I-1)+3)=1.0D0
                     ENDIF
                  ELSE
                     IATNUM(I) = J
                     ATMASS(I) = ATMSS(J)
                  ENDIF
               ENDIF
 20         CONTINUE
         ENDIF
         IF (IPRNT .GE.10) WRITE (*,9010) I,ZSYM(I),IATNUM(I),ATMASS(I)
 9010    FORMAT(I3,2X,A5,2X,I5,2X,F10.5)
 10   CONTINUE

! RAC: very, very lazy..

      IF (ZSYM(NATOMS).EQ.'SV') THEN
         print *, " Assigning the atomic masses for MSEVB system "
         do i = 1,NATOMS
            if (MOD(i,3).eq.0) then
               ATMASS(i) = 15.9994d0
            else
               ATMASS(i) = 1.00794
            endif
         enddo
      ENDIF

!Assign atomic masses to QUIP atoms
      IF (ZSYM(NATOMS).EQ.'QI') ATMASS = QUIPATOMMASS

      IF (NTIPT) THEN
 
         IF(ALLOCATED(IATNUM)) DEALLOCATE(IATNUM) 
         IF(ALLOCATED(ATMASS)) DEALLOCATE(ATMASS) 
           
         ALLOCATE(IATNUM(3*NATOMS/2))
         ALLOCATE(ATMASS(3*NATOMS/2))

         DO I = 1, NATOMS/2
 
            IATNUM(3*(I-1)+1)=8
            IATNUM(3*(I-1)+2)=1
            IATNUM(3*(I-1)+3)=1
            ATMASS(3*(I-1)+1)=16.0D0
            ATMASS(3*(I-1)+2)=1.0D0
            ATMASS(3*(I-1)+3)=1.0D0

         ENDDO

      ELSE IF (PAHAT) THEN

         IF(ALLOCATED(IATNUM)) DEALLOCATE(IATNUM) 
         IF(ALLOCATED(ATMASS)) DEALLOCATE(ATMASS)
           
         ALLOCATE(IATNUM(NRBSITES*NATOMS/2))
         ALLOCATE(ATMASS(NRBSITES*NATOMS/2))
 
         DO I = 1, NATOMS/2

            DO J = 1, NRBSITES

               IF ( J <= NCARBON) THEN
                  IATNUM(NRBSITES*(I-1)+J)=6 
                  ATMASS(NRBSITES*(I-1)+J)=12.0D0
               ELSE
                  IATNUM(NRBSITES*(I-1)+J)=1 
                  ATMASS(NRBSITES*(I-1)+J)=1.0D0
               ENDIF

            ENDDO

         ENDDO

      ENDIF

      RETURN
      END
