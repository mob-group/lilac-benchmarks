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
      SUBROUTINE TCHECK(VEC,Q)
      USE COMMONS
      IMPLICIT NONE

      DOUBLE PRECISION SIZE(64), ROT(3,3)
      LOGICAL WTEST, TEST
      DOUBLE PRECISION X(3*NATOMS),Y(3*NATOMS),Z(3*NATOMS),
     1                 AT(3*NATOMS),BT(3*NATOMS),CT(3*NATOMS),
     2                 XX1(3*NATOMS),YY1(3*NATOMS),ZZ1(3*NATOMS),
     3                 XX2(3*NATOMS),YY2(3*NATOMS),ZZ2(3*NATOMS),
     4                 XX3(3*NATOMS),YY3(3*NATOMS),ZZ3(3*NATOMS),
     5                 VEC(3*NATOMS),Q(3*NATOMS),
     6                 DX(3*NATOMS),DY(3*NATOMS),DZ(3*NATOMS),
     7                 DAT(3*NATOMS),DBT(3*NATOMS),DCT(3*NATOMS),
     8                 DXX1(3*NATOMS),DYY1(3*NATOMS),DZZ1(3*NATOMS),
     9                 DXX2(3*NATOMS),DYY2(3*NATOMS),DZZ2(3*NATOMS),
     A                 DXX3(3*NATOMS),DYY3(3*NATOMS),DZZ3(3*NATOMS), ARGUMENT
      DOUBLE PRECISION ANGLE, HOLEN, PI, PI2, RANGLE, OHZ, HYL, HZL, OL,
     1  TH, PH, PS, SINA, SINB, SINC, COSA, COSB, COSC, SP12, SP13,
     1  SP22, SP23, SP32, SP33, SY12, SY22, SY32, SZ13, SZ23, SZ33, TEMPY,
     1  TEMPX, SP13T, SP12T, SP22T, TEMP, SBEST, TEMPZ

      INTEGER J6, JBEST, J1, J2, J5, I, NMOLS, LBOUND, UBOUND

      DATA ANGLE,HOLEN/104.52D0,0.9572D0/
C
      WTEST=.FALSE.
C
C  Turn off rotation for reaction paths.
C
      IF ((INR.EQ.3).OR.REDOPATH) GOTO 440
      NMOLS=NATOMS/2  !  NMOLS is the number of molecules
      DO 190 J1=1,NMOLS
         J2=3*NMOLS+3*(J1-1)+1
         IF (ABS(DSIN(Q(J2))).LT.0.1D0) THEN
            WTEST=.TRUE.
            TEST=.FALSE.
            PRINT*,'Euler angle going bad - evasive action'
            GOTO 200
          ENDIF
190   CONTINUE
C
C  If it;s the first step find the best orientation anyway.
C
200   IF (ISTATUS.LT.0) WTEST=.TRUE.
C     PRINT*,'WTEST=',WTEST
      IF (WTEST) THEN
         PI=4.0D0*ATAN(1.0D0)
         PI2=2.0D0*PI
         RANGLE=PI*ANGLE/360.0D0
         OHZ=HOLEN*COS(RANGLE)
         HYL=HOLEN*SIN(RANGLE)
         HZL=16.0D0*OHZ/18.0D0
         OL=-OHZ+HZL
C        PRINT*,'Original angles'
C        WRITE(*,210) (Q(3*NMOLS+J1),J1=1,3*(NMOLS))
C210      FORMAT(3F20.16)
         IF ((NMOLS.EQ.64).OR.(NMOLS.EQ.67)) THEN
         LBOUND=1
         UBOUND=3
220      DO 280 J5=LBOUND,UBOUND
            DO 240 I=1,NMOLS
               AT(I)=Q(3*NMOLS+3*(I-1)+1)
               BT(I)=Q(3*NMOLS+3*(I-1)+2)
               CT(I)=Q(3*NMOLS+3*(I-1)+3)
               TH=AT(I)
               PH=BT(I)
               PS=CT(I)
               SINA=SIN(TH)
               SINB=SIN(PH)
               SINC=SIN(PS)
               COSA=COS(TH)
               COSB=COS(PH)
               COSC=COS(PS)
               SP12=-(SINC*COSB+COSA*SINB*COSC)
               SP13=SINA*SINB
               SP22=-SINC*SINB+COSA*COSB*COSC
               SP23=-SINA*COSB
               SP32=SINA*COSC
               SP33=COSA
               SY12=SP12*HYL
               SY22=SP22*HYL
               SY32=SP32*HYL
               SZ13=SP13*HZL
               SZ23=SP23*HZL
               SZ33=SP33*HZL
               X(I)=Q(3*(I-1)+1)
               Y(I)=Q(3*(I-1)+2)
               Z(I)=Q(3*(I-1)+3)
               XX1(I)=SY12+SZ13+X(I)
               YY1(I)=SY22+SZ23+Y(I)
               ZZ1(I)=SY32+SZ33+Z(I)
               XX2(I)=-SY12+SZ13+X(I)
               YY2(I)=-SY22+SZ23+Y(I)
               ZZ2(I)=-SY32+SZ33+Z(I)
               XX3(I)=SP13*OL   +X(I)
               YY3(I)=SP23*OL   +Y(I)
               ZZ3(I)=SP33*OL   +Z(I)
C              IF (TEST) THEN
C                 WRITE(*,15) X(I),Y(I),Z(I),TH,PH,PS
C              ENDIF
240         CONTINUE
            DO 250 J1=1,NMOLS
               IF (MOD(J5-1,4).EQ.1) THEN
                  TEMPY=YY1(J1)
                  YY1(J1)=ZZ1(J1)
                  ZZ1(J1)=-TEMPY
                  TEMPY=YY2(J1)
                  YY2(J1)=ZZ2(J1)
                  ZZ2(J1)=-TEMPY
                  TEMPY=YY3(J1)
                  YY3(J1)=ZZ3(J1)
                  ZZ3(J1)=-TEMPY
               ELSE IF (MOD(J5-1,4).EQ.2) THEN
                  YY1(J1)=-YY1(J1)
                  ZZ1(J1)=-ZZ1(J1)
                  YY2(J1)=-YY2(J1)
                  ZZ2(J1)=-ZZ2(J1)
                  YY3(J1)=-YY3(J1)
                  ZZ3(J1)=-ZZ3(J1)
               ELSE IF (MOD(J5-1,4).EQ.3) THEN
                  TEMPY=YY1(J1)
                  YY1(J1)=-ZZ1(J1)
                  ZZ1(J1)=TEMPY
                  TEMPY=YY2(J1)
                  YY2(J1)=-ZZ2(J1)
                  ZZ2(J1)=TEMPY
                  TEMPY=YY3(J1)
                  YY3(J1)=-ZZ3(J1)
                  ZZ3(J1)=TEMPY
               ENDIF
               IF (MOD((J5-1)/4,4).EQ.1) THEN
                  TEMPX=XX1(J1)
                  XX1(J1)=ZZ1(J1)
                  ZZ1(J1)=-TEMPX
                  TEMPX=XX2(J1)
                  XX2(J1)=ZZ2(J1)
                  ZZ2(J1)=-TEMPX
                  TEMPX=XX3(J1)
                  XX3(J1)=ZZ3(J1)
                  ZZ3(J1)=-TEMPX
               ELSE IF (MOD((J5-1)/4,4).EQ.2) THEN
                  XX1(J1)=-XX1(J1)
                  ZZ1(J1)=-ZZ1(J1)
                  XX2(J1)=-XX2(J1)
                  ZZ2(J1)=-ZZ2(J1)
                  XX3(J1)=-XX3(J1)
                  ZZ3(J1)=-ZZ3(J1)
               ELSE IF (MOD((J5-1)/4,4).EQ.3) THEN
                  TEMPX=XX1(J1)
                  XX1(J1)=-ZZ1(J1)
                  ZZ1(J1)=TEMPX
                  TEMPX=XX2(J1)
                  XX2(J1)=-ZZ2(J1)
                  ZZ2(J1)=TEMPX
                  TEMPX=XX3(J1)
                  XX3(J1)=-ZZ3(J1)
                  ZZ3(J1)=TEMPX
               ENDIF
               IF (MOD((J5-1)/16,4).EQ.1) THEN
                  TEMPX=XX1(J1)
                  XX1(J1)=YY1(J1)
                  YY1(J1)=-TEMPX
                  TEMPX=XX2(J1)
                  XX2(J1)=YY2(J1)
                  YY2(J1)=-TEMPX
                  TEMPX=XX3(J1)
                  XX3(J1)=YY3(J1)
                  YY3(J1)=-TEMPX
               ELSE IF (MOD((J5-1)/16,4).EQ.2) THEN
                  XX1(J1)=-XX1(J1)
                  YY1(J1)=-YY1(J1)
                  XX2(J1)=-XX2(J1)
                  YY2(J1)=-YY2(J1)
                  XX3(J1)=-XX3(J1)
                  YY3(J1)=-YY3(J1)
               ELSE IF (MOD((J5-1)/16,4).EQ.3) THEN
                  TEMPX=XX1(J1)
                  XX1(J1)=-YY1(J1)
                  YY1(J1)=TEMPX
                  TEMPX=XX2(J1)
                  XX2(J1)=-YY2(J1)
                  YY2(J1)=TEMPX
                  TEMPX=XX3(J1)
                  XX3(J1)=-YY3(J1)
                  YY3(J1)=TEMPX
               ENDIF
250         CONTINUE

            DO 260 I=1,NMOLS
               X(I)=(XX1(I)+XX2(I)+16.0D0*XX3(I))/18.0D0
               Y(I)=(YY1(I)+YY2(I)+16.0D0*YY3(I))/18.0D0
               Z(I)=(ZZ1(I)+ZZ2(I)+16.0D0*ZZ3(I))/18.0D0

               AT(I)=DACOS((ZZ3(I)-Z(I))/OL)
               COSA=((ZZ3(I)-Z(I))/OL)
               SINA=DSIN(AT(I))

               BT(I)=DACOS(-(YY3(I)-Y(I))/(OL*SINA))
               COSB=-(YY3(I)-Y(I))/(OL*SINA)
               SINB=DSIN(BT(I))

               CT(I)=DACOS((ZZ1(I)-ZZ2(I))/(2.0D0*HYL*SINA))
               COSC=DCOS(CT(I))
               SINC=DSIN(CT(I))

               SP13T=(XX3(I)-X(I))/OL
               SP12T=(XX1(I)-XX2(I))/(2.0D0*HYL)
               SP22T=(YY1(I)-YY2(I))/(2.0D0*HYL)
               SP13=SINA*SINB
               SP12=-(SINC*COSB+COSA*SINB*COSC)
               SP22=-SINC*SINB+COSA*COSB*COSC
               IF (DABS(SP12T/SP12+1.0D0).LT.1.0D-4) THEN
               TEMP=AT(I)
                  AT(I)=PI2-AT(I)
                  BT(I)=PI-BT(I)
                  CT(I)=PI-CT(I)
                  SINA=-SINA
                  COSB=-COSB
                  COSC=-COSC
                  GOTO 260
               ENDIF
               IF (DABS(SP12T/SP12-1.0D0).GT.1.0D-4) THEN
                  IF (DABS(SP22T/SP22-1.0D0).LT.1.0D-4) THEN
C                    PRINT*,'Inconsistent errors - QUIT'
c                    STOP
                  ENDIF
                  IF (DABS(SP13T/SP13+1.0D0).LT.1.0D-4) THEN
                     TEMP=BT(I)
                     BT(I)=PI2-BT(I)
                     SINB=-SINB
                  ELSE
                     TEMP=CT(I)
                     CT(I)=PI2-CT(I)
                     SINC=-SINC
                  ENDIF
               ENDIF
260         CONTINUE
            IF (TEST) GOTO 400
            SIZE(J5)=1.0D6
            DO 270 J6=1,NMOLS
               IF (DABS(DSIN(AT(J6))).LT.SIZE(J5)) THEN
                  SIZE(J5)=DABS(DSIN(AT(J6)))
               ENDIF
270         CONTINUE
280      CONTINUE
         JBEST=1
         SBEST=-1.0D6
         PRINT*,' Smallest sin theta values are'
         DO 300 J6=1,3
            WRITE(*,290) J6,SIZE(J6)
290         FORMAT(I3,F20.10)
            IF (SIZE(J6).GT.SBEST) THEN
               JBEST=J6
               SBEST=SIZE(J6)
            ENDIF
300      CONTINUE
         TEST=.TRUE.
         J5=JBEST
         PRINT*
         PRINT*,' Best is',J5
         LBOUND=JBEST
         UBOUND=JBEST
         GOTO 220
C
C  This branch for clusters
C
         ELSE
         DO 320 I=1,NMOLS
            AT(I)=Q(3*NMOLS+3*(I-1)+1)
            BT(I)=Q(3*NMOLS+3*(I-1)+2)
            CT(I)=Q(3*NMOLS+3*(I-1)+3)
            TH=AT(I)
            PH=BT(I)
            PS=CT(I)
            SINA=SIN(TH)
            SINB=SIN(PH)
            SINC=SIN(PS)
            COSA=COS(TH)
            COSB=COS(PH)
            COSC=COS(PS)
            SP12=-(SINC*COSB+COSA*SINB*COSC)
            SP13=SINA*SINB
            SP22=-SINC*SINB+COSA*COSB*COSC
            SP23=-SINA*COSB
            SP32=SINA*COSC
            SP33=COSA
            SY12=SP12*HYL
            SY22=SP22*HYL
            SY32=SP32*HYL
            SZ13=SP13*HZL
            SZ23=SP23*HZL
            SZ33=SP33*HZL
            X(I)=Q(3*(I-1)+1)
            Y(I)=Q(3*(I-1)+2)
            Z(I)=Q(3*(I-1)+3)
            XX1(I)=SY12+SZ13+X(I)
            YY1(I)=SY22+SZ23+Y(I)
            ZZ1(I)=SY32+SZ33+Z(I)
            XX2(I)=-SY12+SZ13+X(I)
            YY2(I)=-SY22+SZ23+Y(I)
            ZZ2(I)=-SY32+SZ33+Z(I)
            XX3(I)=SP13*OL   +X(I)
            YY3(I)=SP23*OL   +Y(I)
            ZZ3(I)=SP33*OL   +Z(I)
C
C  This part is to rotate the vector so that eigenvector-following
C  can continue
C
            DAT(I)=Q(3*NMOLS+3*(I-1)+1)+VEC(3*NMOLS+3*(I-1)+1)/1000 
            DBT(I)=Q(3*NMOLS+3*(I-1)+2)+VEC(3*NMOLS+3*(I-1)+2)/1000 
            DCT(I)=Q(3*NMOLS+3*(I-1)+3)+VEC(3*NMOLS+3*(I-1)+3)/1000 
            TH=DAT(I)
            PH=DBT(I)
            PS=DCT(I)
            SINA=SIN(TH)
            SINB=SIN(PH)
            SINC=SIN(PS)
            COSA=COS(TH)
            COSB=COS(PH)
            COSC=COS(PS)
            SP12=-(SINC*COSB+COSA*SINB*COSC)
            SP13=SINA*SINB
            SP22=-SINC*SINB+COSA*COSB*COSC
            SP23=-SINA*COSB
            SP32=SINA*COSC
            SP33=COSA
            SY12=SP12*HYL
            SY22=SP22*HYL
            SY32=SP32*HYL
            SZ13=SP13*HZL
            SZ23=SP23*HZL
            SZ33=SP33*HZL
            DX(I)=Q(3*(I-1)+1)+VEC(3*(I-1)+1)/1000
            DY(I)=Q(3*(I-1)+2)+VEC(3*(I-1)+2)/1000
            DZ(I)=Q(3*(I-1)+3)+VEC(3*(I-1)+3)/1000
            DXX1(I)=SY12+SZ13+DX(I)
            DYY1(I)=SY22+SZ23+DY(I)
            DZZ1(I)=SY32+SZ33+DZ(I)
            DXX2(I)=-SY12+SZ13+DX(I)
            DYY2(I)=-SY22+SZ23+DY(I)
            DZZ2(I)=-SY32+SZ33+DZ(I)
            DXX3(I)=SP13*OL   +DX(I)
            DYY3(I)=SP23*OL   +DY(I)
            DZZ3(I)=SP33*OL   +DZ(I)
320      CONTINUE
C        PRINT*,'Rotating coordinates'
         ROT(2,2)= DCOS(10.0D0*PI/180.0D0)
         ROT(2,3)=-DSIN(10.0D0*PI/180.0D0)
         ROT(3,2)= DSIN(10.0D0*PI/180.0D0)
         ROT(3,3)= DCOS(10.0D0*PI/180.0D0)
         DO 360 J5=1,36
            DO 330 J1=1,NMOLS
               TEMPY=ROT(2,2)*YY1(J1)+ROT(2,3)*ZZ1(J1)
               TEMPZ=ROT(3,2)*YY1(J1)+ROT(3,3)*ZZ1(J1)
               YY1(J1)=TEMPY
               ZZ1(J1)=TEMPZ
               TEMPY=ROT(2,2)*YY2(J1)+ROT(2,3)*ZZ2(J1)
               TEMPZ=ROT(3,2)*YY2(J1)+ROT(3,3)*ZZ2(J1)
               YY2(J1)=TEMPY
               ZZ2(J1)=TEMPZ
               TEMPY=ROT(2,2)*YY3(J1)+ROT(2,3)*ZZ3(J1)
               TEMPZ=ROT(3,2)*YY3(J1)+ROT(3,3)*ZZ3(J1)
               YY3(J1)=TEMPY
               ZZ3(J1)=TEMPZ
               TEMPY=ROT(2,2)*DYY1(J1)+ROT(2,3)*DZZ1(J1)
               TEMPZ=ROT(3,2)*DYY1(J1)+ROT(3,3)*DZZ1(J1)
               DYY1(J1)=TEMPY
               DZZ1(J1)=TEMPZ
               TEMPY=ROT(2,2)*DYY2(J1)+ROT(2,3)*DZZ2(J1)
               TEMPZ=ROT(3,2)*DYY2(J1)+ROT(3,3)*DZZ2(J1)
               DYY2(J1)=TEMPY
               DZZ2(J1)=TEMPZ
               TEMPY=ROT(2,2)*DYY3(J1)+ROT(2,3)*DZZ3(J1)
               TEMPZ=ROT(3,2)*DYY3(J1)+ROT(3,3)*DZZ3(J1)
               DYY3(J1)=TEMPY
               DZZ3(J1)=TEMPZ
330         CONTINUE

            DO 340 I=1,NMOLS
               X(I)=(XX1(I)+XX2(I)+16.0D0*XX3(I))/18.0D0
               Y(I)=(YY1(I)+YY2(I)+16.0D0*YY3(I))/18.0D0
               Z(I)=(ZZ1(I)+ZZ2(I)+16.0D0*ZZ3(I))/18.0D0

               ARGUMENT=MIN((ZZ3(I)-Z(I))/OL,1.0D0)
               ARGUMENT=MAX(ARGUMENT,-1.0D0)
               AT(I)=DACOS(ARGUMENT)
               COSA=((ZZ3(I)-Z(I))/OL)
               SINA=DSIN(AT(I))
               IF (SINA.EQ.0.0D0) SINA=1.0D-100

               ARGUMENT=MIN(-(YY3(I)-Y(I))/(OL*SINA),1.0D0)
               ARGUMENT=MAX(ARGUMENT,-1.0D0)
               BT(I)=DACOS(ARGUMENT)
               COSB=-(YY3(I)-Y(I))/(OL*SINA)
               SINB=DSIN(BT(I))

               ARGUMENT=MIN((ZZ1(I)-ZZ2(I))/(2.0D0*HYL*SINA),1.0D0)
               ARGUMENT=MAX(ARGUMENT,-1.0D0)
               CT(I)=DACOS(ARGUMENT)
               COSC=DCOS(CT(I))
               SINC=DSIN(CT(I))

               SP13T=(XX3(I)-X(I))/OL
               SP12T=(XX1(I)-XX2(I))/(2.0D0*HYL)
               SP22T=(YY1(I)-YY2(I))/(2.0D0*HYL)
               SP13=SINA*SINB
               SP12=-(SINC*COSB+COSA*SINB*COSC)
               SP22=-SINC*SINB+COSA*COSB*COSC
               IF (SP12.NE.0.0D0) THEN
                  IF (DABS(SP12T/SP12+1.0D0).LT.1.0D-4) THEN
                     TEMP=AT(I)
                     AT(I)=PI2-AT(I)
                     BT(I)=PI-BT(I)
                     CT(I)=PI-CT(I)
                     SINA=-SINA
                     COSB=-COSB
                     COSC=-COSC
                     GOTO 345
                  ENDIF
                  IF (DABS(SP12T/SP12-1.0D0).GT.1.0D-4) THEN
                     IF (DABS(SP22T/SP22-1.0D0).LT.1.0D-4) THEN
C                       PRINT*,'Inconsistent errors - QUIT'
c                       STOP
                     ENDIF
                     IF (SP13.NE.0.0D0) THEN
                        IF (DABS(SP13T/SP13+1.0D0).LT.1.0D-4) THEN
                           TEMP=BT(I)
                           BT(I)=PI2-BT(I)
                           SINB=-SINB
                        ELSE
                           TEMP=CT(I)
                           CT(I)=PI2-CT(I)
                           SINC=-SINC
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
345            DX(I)=(DXX1(I)+DXX2(I)+16.0D0*DXX3(I))/18.0D0
               DY(I)=(DYY1(I)+DYY2(I)+16.0D0*DYY3(I))/18.0D0
               DZ(I)=(DZZ1(I)+DZZ2(I)+16.0D0*DZZ3(I))/18.0D0

               ARGUMENT=MIN((DZZ3(I)-DZ(I))/OL,1.0D0)
               ARGUMENT=MAX(ARGUMENT,-1.0D0)
               DAT(I)=DACOS(ARGUMENT)
               COSA=((DZZ3(I)-DZ(I))/OL)
               SINA=DSIN(DAT(I))
               IF (SINA.EQ.0.0D0) SINA=1.0D-100

               ARGUMENT=MIN(-(DYY3(I)-DY(I))/(OL*SINA),1.0D0)
               ARGUMENT=MAX(ARGUMENT,-1.0D0)
               DBT(I)=DACOS(ARGUMENT)
               COSB=-(DYY3(I)-DY(I))/(OL*SINA)
               SINB=DSIN(DBT(I))

               ARGUMENT=MIN((DZZ1(I)-DZZ2(I))/(2.0D0*HYL*SINA),1.0D0)
               ARGUMENT=MAX(ARGUMENT,-1.0D0)
               DCT(I)=DACOS(ARGUMENT)
               COSC=DCOS(DCT(I))
               SINC=DSIN(DCT(I))

               SP13T=(DXX3(I)-DX(I))/OL
               SP12T=(DXX1(I)-DXX2(I))/(2.0D0*HYL)
               SP22T=(DYY1(I)-DYY2(I))/(2.0D0*HYL)
               SP13=SINA*SINB
               SP12=-(SINC*COSB+COSA*SINB*COSC)
               SP22=-SINC*SINB+COSA*COSB*COSC
               IF (SP12.NE.0.0D0) THEN
                  IF (DABS(SP12T/SP12+1.0D0).LT.1.0D-4) THEN
                     TEMP=DAT(I)
                     DAT(I)=PI2-DAT(I)
                     DBT(I)=PI-DBT(I)
                     DCT(I)=PI-DCT(I)
                     SINA=-SINA
                     COSB=-COSB
                     COSC=-COSC
                     GOTO 340
                  ENDIF
                  IF (DABS(SP12T/SP12-1.0D0).GT.1.0D-4) THEN
                     IF (DABS(SP22T/SP22-1.0D0).LT.1.0D-4) THEN
C                    PRINT*,'Inconsistent errors - QUIT'
c                    STOP
                     ENDIF
                     IF (SP13.NE.0.0D0) THEN
                        IF (DABS(SP13T/SP13+1.0D0).LT.1.0D-4) THEN
                           TEMP=DBT(I)
                           DBT(I)=PI2-DBT(I)
                           SINB=-SINB
                        ELSE
                           TEMP=DCT(I)
                           DCT(I)=PI2-DCT(I)
                           SINC=-SINC
                        ENDIF
                     ENDIF
                  ENDIF
               ENDIF
340         CONTINUE
            SIZE(J5)=1.0D6
            DO 350 J6=1,NMOLS
               IF (DABS(DSIN(AT(J6))).LT.SIZE(J5)) THEN
                  SIZE(J5)=DABS(DSIN(AT(J6)))
               ENDIF
350         CONTINUE
360      CONTINUE
         JBEST=1
         SBEST=-1.0D6
         DO 370 J6=1,36
            IF (SIZE(J6).GT.SBEST) THEN
               JBEST=J6
               SBEST=SIZE(J6)
            ENDIF
370      CONTINUE
         ROT(2,2)= DCOS(JBEST*10.0D0*PI/180.0D0)
         ROT(2,3)=-DSIN(JBEST*10.0D0*PI/180.0D0)
         ROT(3,2)= DSIN(JBEST*10.0D0*PI/180.0D0)
         ROT(3,3)= DCOS(JBEST*10.0D0*PI/180.0D0)
         DO 380 J1=1,NMOLS
            TEMPY=ROT(2,2)*YY1(J1)+ROT(2,3)*ZZ1(J1)
            TEMPZ=ROT(3,2)*YY1(J1)+ROT(3,3)*ZZ1(J1)
            YY1(J1)=TEMPY
            ZZ1(J1)=TEMPZ
            TEMPY=ROT(2,2)*YY2(J1)+ROT(2,3)*ZZ2(J1)
            TEMPZ=ROT(3,2)*YY2(J1)+ROT(3,3)*ZZ2(J1)
            YY2(J1)=TEMPY
            ZZ2(J1)=TEMPZ
            TEMPY=ROT(2,2)*YY3(J1)+ROT(2,3)*ZZ3(J1)
            TEMPZ=ROT(3,2)*YY3(J1)+ROT(3,3)*ZZ3(J1)
            YY3(J1)=TEMPY
            ZZ3(J1)=TEMPZ
            TEMPY=ROT(2,2)*DYY1(J1)+ROT(2,3)*DZZ1(J1)
            TEMPZ=ROT(3,2)*DYY1(J1)+ROT(3,3)*DZZ1(J1)
            DYY1(J1)=TEMPY
            DZZ1(J1)=TEMPZ
            TEMPY=ROT(2,2)*DYY2(J1)+ROT(2,3)*DZZ2(J1)
            TEMPZ=ROT(3,2)*DYY2(J1)+ROT(3,3)*DZZ2(J1)
            DYY2(J1)=TEMPY
            DZZ2(J1)=TEMPZ
            TEMPY=ROT(2,2)*DYY3(J1)+ROT(2,3)*DZZ3(J1)
            TEMPZ=ROT(3,2)*DYY3(J1)+ROT(3,3)*DZZ3(J1)
            DYY3(J1)=TEMPY
            DZZ3(J1)=TEMPZ
380      CONTINUE

         DO 390 I=1,NMOLS
            X(I)=(XX1(I)+XX2(I)+16.0D0*XX3(I))/18.0D0
            Y(I)=(YY1(I)+YY2(I)+16.0D0*YY3(I))/18.0D0
            Z(I)=(ZZ1(I)+ZZ2(I)+16.0D0*ZZ3(I))/18.0D0

            AT(I)=DACOS((ZZ3(I)-Z(I))/OL)
            COSA=((ZZ3(I)-Z(I))/OL)
            SINA=DSIN(AT(I))

            ARGUMENT=MIN(-(YY3(I)-Y(I))/(OL*SINA),1.0D0)
            ARGUMENT=MAX(ARGUMENT,-1.0D0)
            BT(I)=DACOS(ARGUMENT)
            COSB=-(YY3(I)-Y(I))/(OL*SINA)
            SINB=DSIN(BT(I))

            ARGUMENT=MIN((ZZ1(I)-ZZ2(I))/(2.0D0*HYL*SINA),1.0D0)
            ARGUMENT=MAX(ARGUMENT,-1.0D0)
            CT(I)=DACOS(ARGUMENT)
            COSC=DCOS(CT(I))
            SINC=DSIN(CT(I))

            SP13T=(XX3(I)-X(I))/OL
            SP12T=(XX1(I)-XX2(I))/(2.0D0*HYL)
            SP22T=(YY1(I)-YY2(I))/(2.0D0*HYL)
            SP13=SINA*SINB
            SP12=-(SINC*COSB+COSA*SINB*COSC)
            SP22=-SINC*SINB+COSA*COSB*COSC
            IF (SP12.NE.0.0D0) THEN
               IF (DABS(SP12T/SP12+1.0D0).LT.1.0D-4) THEN
                  TEMP=AT(I)
                  AT(I)=PI2-AT(I)
                  BT(I)=PI-BT(I)
                  CT(I)=PI-CT(I)
                  SINA=-SINA
                  COSB=-COSB
                  COSC=-COSC
                  GOTO 395
               ENDIF
               IF (DABS(SP12T/SP12-1.0D0).GT.1.0D-4) THEN
                  IF (DABS(SP22T/SP22-1.0D0).LT.1.0D-4) THEN
C                 PRINT*,'Inconsitent errors - QUIT'
c                 STOP
               ENDIF
            ENDIF
            IF (SP13.NE.0.0D0) THEN
               IF (DABS(SP13T/SP13+1.0D0).LT.1.0D-4) THEN
                  TEMP=BT(I)
                     BT(I)=PI2-BT(I)
                     SINB=-SINB
                  ELSE
                     TEMP=CT(I)
                     CT(I)=PI2-CT(I)
                     SINC=-SINC
                  ENDIF
               ENDIF
            ENDIF
395         DX(I)=(DXX1(I)+DXX2(I)+16.0D0*DXX3(I))/18.0D0
            DY(I)=(DYY1(I)+DYY2(I)+16.0D0*DYY3(I))/18.0D0
            DZ(I)=(DZZ1(I)+DZZ2(I)+16.0D0*DZZ3(I))/18.0D0

            DAT(I)=DACOS((DZZ3(I)-DZ(I))/OL)
            COSA=((DZZ3(I)-DZ(I))/OL)
            SINA=DSIN(DAT(I))

            ARGUMENT=MIN(-(DYY3(I)-DY(I))/(OL*SINA),1.0D0)
            ARGUMENT=MAX(ARGUMENT,-1.0D0)
            DBT(I)=DACOS(ARGUMENT)
            COSB=-(DYY3(I)-DY(I))/(OL*SINA)
            SINB=DSIN(DBT(I))

            ARGUMENT=MIN((DZZ1(I)-DZZ2(I))/(2.0D0*HYL*SINA),1.0D0)
            ARGUMENT=MAX(ARGUMENT,-1.0D0)
            DCT(I)=DACOS(ARGUMENT)
            COSC=DCOS(DCT(I))
            SINC=DSIN(DCT(I))

            SP13T=(DXX3(I)-DX(I))/OL
            SP12T=(DXX1(I)-DXX2(I))/(2.0D0*HYL)
            SP22T=(DYY1(I)-DYY2(I))/(2.0D0*HYL)
            SP13=SINA*SINB
            SP12=-(SINC*COSB+COSA*SINB*COSC)
            SP22=-SINC*SINB+COSA*COSB*COSC
            IF (SP12.NE.0.0D0) THEN
               IF (DABS(SP12T/SP12+1.0D0).LT.1.0D-4) THEN
                  TEMP=DAT(I)
                  DAT(I)=PI2-DAT(I)
                  DBT(I)=PI-DBT(I)
                  DCT(I)=PI-DCT(I)
                  SINA=-SINA
                  COSB=-COSB
                  COSC=-COSC
                  GOTO 390
               ENDIF
               IF (DABS(SP12T/SP12-1.0D0).GT.1.0D-4) THEN
                  IF (DABS(SP22T/SP22-1.0D0).LT.1.0D-4) THEN
C                 PRINT*,'Inconsistent errors - QUIT'
c                 STOP
                  ENDIF
                  IF (SP13.NE.0.0D0) THEN
                     IF (DABS(SP13T/SP13+1.0D0).LT.1.0D-4) THEN
                        TEMP=DBT(I)
                        DBT(I)=PI2-DBT(I)
                        SINB=-SINB
                     ELSE
                        TEMP=DCT(I)
                        DCT(I)=PI2-DCT(I)
                        SINC=-SINC
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
390      CONTINUE
         ENDIF
400      CONTINUE
C        REWIND(18)
         TEMP=0.0D0
         DO 430 I=1,NMOLS
C           WRITE(*,410) X(I),Y(I),Z(I),AT(I),BT(I),CT(I)
C410         FORMAT(6F13.6)
            Q(3*(I-1)+1)=X(I)
            Q(3*(I-1)+2)=Y(I)
            Q(3*(I-1)+3)=Z(I)
            Q(3*NMOLS+3*(I-1)+1)=AT(I)
            Q(3*NMOLS+3*(I-1)+2)=BT(I)
            Q(3*NMOLS+3*(I-1)+3)=CT(I)
C           WRITE(18,420) AT(I),BT(I),CT(I)
420         FORMAT(3F20.16)
            VEC(3*(I-1)+1)=X(I)-DX(I)
            VEC(3*(I-1)+2)=Y(I)-DY(I)
            VEC(3*(I-1)+3)=Z(I)-DZ(I)
            VEC(3*NMOLS+3*(I-1)+1)=AT(I)-DAT(I)
            VEC(3*NMOLS+3*(I-1)+2)=BT(I)-DBT(I)
            VEC(3*NMOLS+3*(I-1)+3)=CT(I)-DCT(I)
            TEMP=TEMP+VEC(3*(I-1)+1)**2
            TEMP=TEMP+VEC(3*(I-1)+2)**2
            TEMP=TEMP+VEC(3*(I-1)+3)**2
            TEMP=TEMP+VEC(3*NMOLS+3*(I-1)+1)**2
            TEMP=TEMP+VEC(3*NMOLS+3*(I-1)+2)**2
            TEMP=TEMP+VEC(3*NMOLS+3*(I-1)+3)**2
430      CONTINUE
         TEMP=DSQRT(TEMP)
         IF (TEMP.GT.0.0D0) THEN
         DO 435 I=1,NMOLS
            VEC(3*(I-1)+1)=VEC(3*(I-1)+1)/TEMP 
            VEC(3*(I-1)+2)=VEC(3*(I-1)+2)/TEMP 
            VEC(3*(I-1)+3)=VEC(3*(I-1)+3)/TEMP 
            VEC(3*NMOLS+3*(I-1)+1)=VEC(3*NMOLS+3*(I-1)+1)/TEMP
            VEC(3*NMOLS+3*(I-1)+2)=VEC(3*NMOLS+3*(I-1)+2)/TEMP
            VEC(3*NMOLS+3*(I-1)+3)=VEC(3*NMOLS+3*(I-1)+3)/TEMP
435      CONTINUE
         ENDIF
      ENDIF
440   RETURN
      END
