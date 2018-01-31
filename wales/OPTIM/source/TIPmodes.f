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
      SUBROUTINE H2OMODES(NMOLS,IPOT,COORDS,DIAG)
      USE MODHESS
      IMPLICIT NONE
C ----------------------------------------------------------------------
C 
C       NORMAL MODE ANALYSIS FOR WATER MOLECULES 
C       TO THREE TRANSLATIONAL MODES AND THREE ROTATINAL MODES 
C                          CODED BY H. TANAKA      1988, 01, 08 
C                               SEE J.CHEM.PHYS. 87, 6070 (1987) 
C       FIRST AND SECOND DERIVATIVES ARE EXPRESSED IN ANALYTICAL FORM 
C       EWALD SUM IS NOT TAKEN INTO ACCOUNT 
C ----------------------------------------------------------------------
      INTEGER NMOLS, IPOT, N, NN, N3, N6, NTHREE, M, L, K, L3, L2, L1, K3, K2, K1, I3, I2, I1,
     1        J1, J2, I, J3, J
      DOUBLE PRECISION RC, BXL, RC2, RL, RQ, RQ3, RQ6, EP, EMUPID, P1, P2, P3, P4, P5, ANUM,
     1                QQ, AD1, AD2, ANGLE, HOLEN, COLEN, QE, UJ, SG, WM, CSINTH, PI, PID,
     2                QS, SW, RANGLE, OHZ, FZ4(NMOLS), SS, CC33O, CC23O, CC13O, BC23O, BC13O, BB23O,
     3                BB13O, AC33O, AC23O, AC13O, AB23O, AB13O, AA33O, AA23O, AA13O, CC33C,
     4                CC23C, CC13C, BC23C, BC13C, BB23C,
     1 P6, HYL, HZL, OL, CL, ENU, RIX, RIY, RIZ, PN, PNI, EMU, RDENS, DENS, VOL,
     1 TH, PH, PS, SINA, SINB, SINC, COSA, COSB, COSC, COSAP, AP,
     1 SINAP, COSCP, CP, SINCP, COSBP, BP, SINBP, SP12, DDA12, DDB12, DDC12, DAA12,
     1 DAB12, DBB12, DAC12, DBC12, DCC12, SP13, DDA13, DDB13, DDC13, DAA13, DAB13, DAC13,
     1 DBB13, DBC13, DCC13, SP22, DDA22, DDB22, DDC22, DAA22, DAB22, DAC22, DBB22, DCC22,
     1 DBC22, SP23, DDA23, DDB23, DDC23, DAA23, DAB23, DAC23, DBB23, DBC23, DCC23, SP32, DDA32,
     1 DDC32, DAA32, DAC32, DCC32, SP33, DDA33, DDC33, DAA33, DAC33, DCC33, SY12, SY22, SY32,
     1 SZ13, SZ23, SZ33, A12H, A22H, A32H, A13H, A23H, A33H, B12H, B22H, B13H, B23H, C12H, C22H,
     1 C32H, C13H, C23H, C33H, A13C, A23C, A33C, B13C, B23C, C13C, C23C, C33C, A13O, A23O, A33O,
     1 B13O, B23O, C13O, C23O, C33O, AA12H, AA22H, AA32H, AA13H, AA23H, AA33H, AB12H, AB22H, AB13H,
     1 AB23H, AC12H, AC22H, AC32H, AC13H, AC23H, AC33H, BB12H, BB22H, BB13H, BB23H, BC12H, BC22H, BC13H,
     1 BC23H, CC12H, CC22H, CC32H, CC13H, CC23H, CC33H, AA13C, AA23C, AA33C, AB13C, AB23C, AC13C, AC23C,
     1 AC33C, BB13C, FX(NMOLS), FY(NMOLS), FZ(NMOLS), FA(NMOLS), FB(NMOLS), FC(NMOLS), FX0(NMOLS), FY0(NMOLS), 
     & FZ0(NMOLS), FX1(NMOLS), FY1(NMOLS), FZ1(NMOLS), FX2(NMOLS), FY2(NMOLS),
     1 FZ2(NMOLS), FX3(NMOLS), FY3(NMOLS), FZ3(NMOLS), FX4(NMOLS), FY4(NMOLS)
      LOGICAL GTEST,SSTEST
      DOUBLE PRECISION E(3),EV(3,3),T(3,3), 
     *            EKV(6*NMOLS),EKM(3,3,NMOLS),S(3,3) 
      DOUBLE PRECISION X(NMOLS+1),Y(NMOLS+1),Z(NMOLS+1), 
     *            A(NMOLS),B(NMOLS),C(NMOLS),D(NMOLS), 
     *            EA(NMOLS+1),EB(NMOLS+1),EC(NMOLS+1), 
     *            XX1(NMOLS),YY1(NMOLS),ZZ1(NMOLS), 
     *            XX2(NMOLS),YY2(NMOLS),ZZ2(NMOLS), 
     *            XX3(NMOLS),YY3(NMOLS),ZZ3(NMOLS), 
     *            XX4(NMOLS),YY4(NMOLS),ZZ4(NMOLS) 
      DOUBLE PRECISION AX1(NMOLS),AX2(NMOLS),AX3(NMOLS),AX4(NMOLS), 
     *            AY1(NMOLS),AY2(NMOLS),AY3(NMOLS),AY4(NMOLS), 
     *            AZ1(NMOLS),AZ2(NMOLS),AZ3(NMOLS),AZ4(NMOLS), 
     *            BX1(NMOLS),BX2(NMOLS),BX3(NMOLS),BX4(NMOLS), 
     *            BY1(NMOLS),BY2(NMOLS),BY3(NMOLS),BY4(NMOLS), 
     *            CX1(NMOLS),CX2(NMOLS),CX3(NMOLS),CX4(NMOLS), 
     *            CY1(NMOLS),CY2(NMOLS),CY3(NMOLS),CY4(NMOLS), 
     *            CZ1(NMOLS),CZ2(NMOLS),CZ3(NMOLS),CZ4(NMOLS) 
      DOUBLE PRECISION AAX1(NMOLS),AAX2(NMOLS),AAX3(NMOLS),AAX4(NMOLS), 
     *            ABX1(NMOLS),ABX2(NMOLS),ABX3(NMOLS),ABX4(NMOLS), 
     *            ACX1(NMOLS),ACX2(NMOLS),ACX3(NMOLS),ACX4(NMOLS), 
     *            BBX1(NMOLS),BBX2(NMOLS),BBX3(NMOLS),BBX4(NMOLS), 
     *            BCX1(NMOLS),BCX2(NMOLS),BCX3(NMOLS),BCX4(NMOLS), 
     *            CCX1(NMOLS),CCX2(NMOLS),CCX3(NMOLS),CCX4(NMOLS) 
      DOUBLE PRECISION AAY1(NMOLS),AAY2(NMOLS),AAY3(NMOLS),AAY4(NMOLS), 
     *            ABY1(NMOLS),ABY2(NMOLS),ABY3(NMOLS),ABY4(NMOLS), 
     *            ACY1(NMOLS),ACY2(NMOLS),ACY3(NMOLS),ACY4(NMOLS), 
     *            BBY1(NMOLS),BBY2(NMOLS),BBY3(NMOLS),BBY4(NMOLS), 
     *            BCY1(NMOLS),BCY2(NMOLS),BCY3(NMOLS),BCY4(NMOLS), 
     *            CCY1(NMOLS),CCY2(NMOLS),CCY3(NMOLS),CCY4(NMOLS) 
      DOUBLE PRECISION AAZ1(NMOLS),AAZ2(NMOLS),AAZ3(NMOLS),AAZ4(NMOLS), 
     *            ACZ1(NMOLS),ACZ2(NMOLS),ACZ3(NMOLS),ACZ4(NMOLS), 
     *            CCZ1(NMOLS),CCZ2(NMOLS),CCZ3(NMOLS),CCZ4(NMOLS) 
      DOUBLE PRECISION COORDS(6*NMOLS), DIAG(6*NMOLS)
      COMMON /DR/ RC,BXL,RC2,RL,RQ,RQ3,RQ6,EP,EMUPID,P1,P2,P3,P4,P5,P6 
      INTEGER LC(NMOLS)
      COMMON /IN/ N,NN,N3,N6 
      DATA ANUM/6.0225D+23/ 
      DATA QQ,AD1,AD2/0.535D0,6.95D+5,6.00D+2/ 
      DATA ANGLE,HOLEN,COLEN/104.52D0,0.9572D0,0.15D0/ 
C     DATA QE/332.17752D0/ 
      DATA QE/332.0637782D0/  ! hartree to kJ/mol divided by 4.184 times A to bohr
      DATA UJ/4.184D3/ 
      DATA SG,WM/1.0D-8,18.0D0/ 
      INTEGER INFO
      DOUBLE PRECISION TEMPA(18*NMOLS)

      N=NMOLS


!      OPEN(UNIT=1,FILE='INPUTCOORDS',STATUS='UNKNOWN')
!      DO J1 = 1, 2*NMOLS
!         READ(1,*) COORDS(3*J1-2:3*J1)
!      ENDDO

C
C  CSINTH is a fudge to muck around with the Euler angles if
C  one of them is small. Swaps the y and z axes for the molecule in question.
C  IPOT defines the TIPS type as usual. 
C
      CSINTH =0.001D0
C
C  This is the only time IPOT is used!
C
      IF (IPOT.EQ.4) THEN 
         QQ=0.52D0 
         AD1=6.0D+5 
         AD2=6.10D+2 
      ENDIF 
      IF (IPOT.EQ.3) THEN
         QQ=0.417D0 
         AD1=5.82D+5 
         AD2=5.950D+2 
      ENDIF
      IF (IPOT.EQ.1) THEN
         QQ=0.4D0 
         AD1=5.8D+5 
         AD2=6.0D+2 
      ENDIF
      PI=4.0D0*ATAN(1.0D0) 
      PID=2.0D0*PI
      NTHREE=3 
      QS=QQ**2*QE 
      SW=SQRT(QS*UJ*1.0D3/18.0D0)/(6.0D0*PI) 
      P1=AD1/QS 
      P2=-AD2/QS 
      P3=-P1*12.0D0 
      P4=-P2*6.0D0 
      P5=-P3*13.0D0 
      P6=-P4*7.0D0 
      RANGLE=PI*ANGLE/360.0D0 
      OHZ=HOLEN*COS(RANGLE) 
      HYL=HOLEN*SIN(RANGLE) 
      HZL=16.0D0*OHZ/18.0D0 
      OL=-OHZ+HZL 
      IF ((IPOT.EQ.2).OR.(IPOT.EQ.4)) CL=COLEN+OL 
      IF ((IPOT.EQ.1).OR.(IPOT.EQ.3)) CL=OL 
      ENU=QS/ANUM*UJ*1.0D+7 
      NN=N-1 
      N3=N*3 
      N6=N3+N3 
      RIX=(16.0D0*OL*OL+2.0D0*(HYL*HYL+HZL*HZL))/WM 
      RIY=(16.0D0*OL*OL+2.0D0*HZL*HZL)/WM 
      RIZ=2.0D0*HYL*HYL/WM 
      PN=FLOAT(N) 
      PNI=1.0D0/PN 
      EMU=1.0D-10*ENU*ANUM*PNI 
      EMUPID=EMU*PN 
      IF (N.EQ.64) THEN
         IPOT=2
         RDENS=0.9970710D+00          
         DENS=ANUM*SG*SG*SG*RDENS/WM
         VOL=PN/DENS
         BXL=VOL**(1.0D0/3.0D0)
         RC=0.5D0*BXL
      ELSE
         BXL=100000.0 
         RC=10000.0 
      ENDIF
      RC2=RC*RC 
      RL=RC-2.0D0 
      RQ=1.0D0/(RL-RC)**5 
      RQ3=30.0D0/(RL-RC)**5 
      RQ6=60.0D0/(RL-RC)**5 
      DO J1=1,N
         J2=3*(J1-1)
         X(J1)=COORDS(J2+1)
         Y(J1)=COORDS(J2+2)
         Z(J1)=COORDS(J2+3)
         EA(J1)=COORDS(3*N+J2+1)
         EB(J1)=COORDS(3*N+J2+2)
         EC(J1)=COORDS(3*N+J2+3)
      ENDDO

      DO 48 I=1,N 
         TH=EA(I) 
         PH=EB(I) 
         PS=EC(I) 
         SINA=SIN(TH) 
         SINB=SIN(PH) 
         SINC=SIN(PS) 
         COSA=COS(TH) 
         COSB=COS(PH) 
         COSC=COS(PS) 
         COSAP=SINA*SINC
C        PRINT*,'I,ABS(SINA),test=',I,ABS(SINA),ABS(SINA).GT.CSINTH
         IF (ABS(SINA).GT.CSINTH) THEN
            LC(I)=1 
            CALL LC1(NMOLS,SINA,SINB,SINC,COSA,COSB,COSC,HYL,HZL,OL,CL,I,
     1            X,Y,Z,A,B,C,D,EA,EB,EC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,XX3,YY3,ZZ3,XX4,YY4,ZZ4,
     2            AX1,AX2,AX3,AX4,AY1,AY2,AY3,AY4,AZ1,AZ2,AZ3,AZ4,BX1,BX2,BX3,BX4,BY1,BY2,BY3,BY4,CX1,CX2,CX3,CX4,
     3            CY1,CY2,CY3,CY4,CZ1,CZ2,CZ3,CZ4,AAX1,AAX2,AAX3,AAX4,ABX1,ABX2,ABX3,ABX4,ACX1,ACX2,ACX3,ACX4,BBX1,BBX2,
     4            BBX3,BBX4,BCX1,BCX2,BCX3,BCX4,CCX1,CCX2,CCX3,CCX4,AAY1,AAY2,AAY3,AAY4,ABY1,ABY2,ABY3,ABY4,
     5            ACY1,ACY2,ACY3,ACY4,BBY1,BBY2,BBY3,BBY4,BCY1,BCY2,BCY3,BCY4,CCY1,CCY2,CCY3,CCY4,
     6            AAZ1,AAZ2,AAZ3,AAZ4,ACZ1,ACZ2,ACZ3,ACZ4,CCZ1,CCZ2,CCZ3,CCZ4,LC)
         ELSE
            LC(I)=2
            AP=ACOS(COSAP)
            SINAP=SIN(AP)
            COSCP=COSA/SINAP
            CP=ACOS(COSCP)
            SINCP=SINA*COSC/SINAP
            IF(SINCP.LT.0.0D0) THEN
               CP=PID-CP
            ENDIF
            COSBP=-(COSC*SINB+COSA*COSB*SINC)/SINAP
            BP=ACOS(COSBP)
            SINBP=(COSC*COSB-COSA*SINB*SINC)/SINAP
            IF(SINBP.LT.0.0D0) THEN
               BP=PID-BP
            ENDIF
            EA(I)=AP
            EB(I)=BP
            EC(I)=CP
C           CALL LC2(SINAP,SINBP,SINCP,COSAP,COSBP,COSCP,HYL,HZL,OL,CL,I)
            SP12=COSCP*COSBP-COSAP*SINBP*SINCP 
            DDA12=SINAP*SINBP*SINCP 
            DDB12=-COSCP*SINBP-COSAP*COSBP*SINCP 
            DDC12=-SINCP*COSBP-COSAP*SINBP*COSCP 
            DAA12=COSAP*SINBP*SINCP 
            DAB12=SINAP*COSBP*SINCP 
            DBB12=-SP12 
            DAC12=SINAP*SINBP*COSCP 
            DBC12=SINCP*SINBP-COSAP*COSBP*COSCP 
            DCC12=-SP12 
            SP13=-(SINCP*COSBP+COSAP*SINBP*COSCP) 
            DDA13=SINAP*SINBP*COSCP 
            DDB13=SINCP*SINBP-COSAP*COSBP*COSCP 
            DDC13=-COSCP*COSBP+COSAP*SINBP*SINCP 
            DAA13=COSAP*SINBP*COSCP 
            DAB13=SINAP*COSBP*COSCP 
            DAC13=-SINAP*SINBP*SINCP 
            DBB13=-SP13 
            DBC13=COSCP*SINBP+COSAP*COSBP*SINCP 
            DCC13=-SP13 
            SP22=COSCP*SINBP+COSAP*COSBP*SINCP 
            DDA22=-SINAP*COSBP*SINCP 
            DDB22=COSCP*COSBP-COSAP*SINBP*SINCP 
            DDC22=-SINCP*SINBP+COSAP*COSBP*COSCP 
            DAA22=-COSAP*COSBP*SINCP 
            DAB22=SINAP*SINBP*SINCP 
            DAC22=-SINAP*COSBP*COSCP 
            DBB22=-SP22 
            DCC22=-SP22 
            DBC22=-SINCP*COSBP-COSAP*SINBP*COSCP 
            SP23=-SINCP*SINBP+COSAP*COSBP*COSCP 
            DDA23=-SINAP*COSBP*COSCP 
            DDB23=-SINCP*COSBP-COSAP*SINBP*COSCP 
            DDC23=-COSCP*SINBP-COSAP*COSBP*SINCP 
            DAA23=-COSAP*COSBP*COSCP 
            DAB23=SINAP*SINBP*COSCP 
            DAC23=SINAP*COSBP*SINCP 
            DBB23=-SP23 
            DBC23=-COSCP*COSBP+COSAP*SINBP*SINCP 
            DCC23=-SP23 
            SP32=SINAP*SINCP 
            DDA32=COSAP*SINCP 
            DDC32=SINAP*COSCP 
            DAA32=-SP32 
            DAC32=COSAP*COSCP 
            DCC32=-SP32 
            SP33=SINAP*COSCP 
            DDA33=COSAP*COSCP 
            DDC33=-SINAP*SINCP 
            DAA33=-SP33 
            DAC33=-COSAP*SINCP 
            DCC33=-SP33 
            SY12=SP12*HYL 
            SY22=SP22*HYL 
            SY32=SP32*HYL 
            SZ13=SP13*HZL 
            SZ23=SP23*HZL 
            SZ33=SP33*HZL 
            XX1(I)=SY12+SZ13 
            YY1(I)=SY22+SZ23 
            ZZ1(I)=SY32+SZ33 
            XX2(I)=-SY12+SZ13 
            YY2(I)=-SY22+SZ23 
            ZZ2(I)=-SY32+SZ33 
            XX3(I)=SP13*OL 
            YY3(I)=SP23*OL 
            ZZ3(I)=SP33*OL 
            XX4(I)=SP13*CL 
            YY4(I)=SP23*CL 
            ZZ4(I)=SP33*CL 
            A12H=DDA12*HYL 
            A22H=DDA22*HYL 
            A32H=DDA32*HYL 
            A13H=DDA13*HZL 
            A23H=DDA23*HZL 
            A33H=DDA33*HZL 
            B12H=DDB12*HYL 
            B22H=DDB22*HYL 
            B13H=DDB13*HZL 
            B23H=DDB23*HZL 
            C12H=DDC12*HYL 
            C22H=DDC22*HYL 
            C32H=DDC32*HYL 
            C13H=DDC13*HZL 
            C23H=DDC23*HZL 
            C33H=DDC33*HZL 
            A13C=DDA13*CL 
            A23C=DDA23*CL 
            A33C=DDA33*CL 
            B13C=DDB13*CL 
            B23C=DDB23*CL 
            C13C=DDC13*CL 
            C23C=DDC23*CL 
            C33C=DDC33*CL 
            A13O=DDA13*OL 
            A23O=DDA23*OL 
            A33O=DDA33*OL 
            B13O=DDB13*OL 
            B23O=DDB23*OL 
            C13O=DDC13*OL 
            C23O=DDC23*OL 
            C33O=DDC33*OL 
            AA12H=DAA12*HYL 
            AA22H=DAA22*HYL 
            AA32H=DAA32*HYL 
            AA13H=DAA13*HZL 
            AA23H=DAA23*HZL 
            AA33H=DAA33*HZL 
            AB12H=DAB12*HYL 
            AB22H=DAB22*HYL 
            AB13H=DAB13*HZL 
            AB23H=DAB23*HZL 
            AC12H=DAC12*HYL 
            AC22H=DAC22*HYL 
            AC32H=DAC32*HYL 
            AC13H=DAC13*HZL 
            AC23H=DAC23*HZL 
            AC33H=DAC33*HZL 
            BB12H=DBB12*HYL 
            BB22H=DBB22*HYL 
            BB13H=DBB13*HZL 
            BB23H=DBB23*HZL 
            BC12H=DBC12*HYL 
            BC22H=DBC22*HYL 
            BC13H=DBC13*HZL 
            BC23H=DBC23*HZL 
            CC12H=DCC12*HYL 
            CC22H=DCC22*HYL 
            CC32H=DCC32*HYL 
            CC13H=DCC13*HZL 
            CC23H=DCC23*HZL 
            CC33H=DCC33*HZL 
            AA13C=DAA13*CL 
            AA23C=DAA23*CL 
            AA33C=DAA33*CL 
            AB13C=DAB13*CL 
            AB23C=DAB23*CL 
            AC13C=DAC13*CL 
            AC23C=DAC23*CL 
            AC33C=DAC33*CL 
            BB13C=DBB13*CL 
            BB23C=DBB23*CL 
            BC13C=DBC13*CL 
            BC23C=DBC23*CL 
            CC13C=DCC13*CL 
            CC23C=DCC23*CL 
            CC33C=DCC33*CL 
            AA13O=DAA13*OL 
            AA23O=DAA23*OL 
            AA33O=DAA33*OL 
            AB13O=DAB13*OL 
            AB23O=DAB23*OL 
            AC13O=DAC13*OL 
            AC23O=DAC23*OL 
            AC33O=DAC33*OL 
            BB13O=DBB13*OL 
            BB23O=DBB23*OL 
            BC13O=DBC13*OL 
            BC23O=DBC23*OL 
            CC13O=DCC13*OL 
            CC23O=DCC23*OL 
            CC33O=DCC33*OL 
            AX1(I)=A12H+A13H 
            AY1(I)=A22H+A23H 
            AZ1(I)=A32H+A33H 
            AX2(I)=-A12H+A13H 
            AY2(I)=-A22H+A23H 
            AZ2(I)=-A32H+A33H 
            AX3(I)=A13O 
            AY3(I)=A23O 
            AZ3(I)=A33O 
            AX4(I)=A13C 
            AY4(I)=A23C 
            AZ4(I)=A33C 
            BX1(I)=B12H+B13H 
            BY1(I)=B22H+B23H 
            BX2(I)=-B12H+B13H 
            BY2(I)=-B22H+B23H 
            BX3(I)=B13O 
            BY3(I)=B23O 
            BX4(I)=B13C 
            BY4(I)=B23C 
            CX1(I)=C12H+C13H 
            CY1(I)=C22H+C23H 
            CZ1(I)=C32H+C33H 
            CX2(I)=-C12H+C13H 
            CY2(I)=-C22H+C23H 
            CZ2(I)=-C32H+C33H 
            CX3(I)=C13O 
            CY3(I)=C23O 
            CZ3(I)=C33O 
            CX4(I)=C13C 
            CY4(I)=C23C 
            CZ4(I)=C33C 
            AAX1(I)=AA12H+AA13H 
            AAY1(I)=AA22H+AA23H 
            AAZ1(I)=AA32H+AA33H 
            AAX2(I)=-AA12H+AA13H 
            AAY2(I)=-AA22H+AA23H 
            AAZ2(I)=-AA32H+AA33H 
            AAX3(I)=AA13O 
            AAY3(I)=AA23O 
            AAZ3(I)=AA33O 
            AAX4(I)=AA13C 
            AAY4(I)=AA23C 
            AAZ4(I)=AA33C 
            BBX1(I)=BB12H+BB13H 
            BBY1(I)=BB22H+BB23H 
            BBX2(I)=-BB12H+BB13H 
            BBY2(I)=-BB22H+BB23H 
            BBX3(I)=BB13O 
            BBY3(I)=BB23O 
            BBX4(I)=BB13C 
            BBY4(I)=BB23C 
            CCX1(I)=CC12H+CC13H 
            CCY1(I)=CC22H+CC23H 
            CCZ1(I)=CC32H+CC33H 
            CCX2(I)=-CC12H+CC13H 
            CCY2(I)=-CC22H+CC23H 
            CCZ2(I)=-CC32H+CC33H 
            CCX3(I)=CC13O 
            CCY3(I)=CC23O 
            CCZ3(I)=CC33O 
            CCX4(I)=CC13C 
            CCY4(I)=CC23C 
            CCZ4(I)=CC33C 
            ABX1(I)=AB12H+AB13H 
            ABY1(I)=AB22H+AB23H 
            ABX2(I)=-AB12H+AB13H 
            ABY2(I)=-AB22H+AB23H 
            ABX3(I)=AB13O 
            ABY3(I)=AB23O 
            ABX4(I)=AB13C 
            ABY4(I)=AB23C 
            ACX1(I)=AC12H+AC13H 
            ACY1(I)=AC22H+AC23H 
            ACZ1(I)=AC32H+AC33H 
            ACX2(I)=-AC12H+AC13H 
            ACY2(I)=-AC22H+AC23H 
            ACZ2(I)=-AC32H+AC33H 
            ACX3(I)=AC13O 
            ACY3(I)=AC23O 
            ACZ3(I)=AC33O 
            ACX4(I)=AC13C 
            ACY4(I)=AC23C 
            ACZ4(I)=AC33C 
            BCX1(I)=BC12H+BC13H 
            BCY1(I)=BC22H+BC23H 
            BCX2(I)=-BC12H+BC13H 
            BCY2(I)=-BC22H+BC23H 
            BCX3(I)=BC13O 
            BCY3(I)=BC23O 
            BCX4(I)=BC13C 
            BCY4(I)=BC23C 
         ENDIF
48    CONTINUE 
      DO 148 I=1,N 
         TH=EA(I) 
         PH=EB(I) 
         PS=EC(I) 
         SINA=SIN(TH) 
         SINB=SIN(PH) 
         SINC=SIN(PS) 
         COSA=COS(TH) 
         COSB=COS(PH) 
         COSC=COS(PS) 
         EV(1,1)=RIX*COSC**2+RIY*SINC**2 
         EV(1,2)=(RIX-RIY)*SINA*SINC*COSC 
         EV(1,3)=0.0D0 
         EV(2,2)=(RIX*SINC**2+RIY*COSC**2)*SINA**2+RIZ*COSA**2 
         EV(2,3)=RIZ*COSA 
         EV(3,3)=RIZ 
         EV(2,1)=EV(1,2) 
         EV(3,1)=EV(1,3) 
         EV(3,2)=EV(2,3) 
         CALL DSYEV('V','U',NTHREE,EV,3,E,TEMPA,18*NMOLS,INFO)
         J1=3*I-2 
         J2=J1+1 
         J3=J2+1 
         EKV(J1)=1.0D0/SQRT(E(1)) 
         EKV(J2)=1.0D0/SQRT(E(2)) 
         EKV(J3)=1.0D0/SQRT(E(3)) 
         EKM(1,1,I)=EV(1,1) 
         EKM(2,1,I)=EV(2,1) 
         EKM(3,1,I)=EV(3,1) 
         EKM(1,2,I)=EV(1,2) 
         EKM(2,2,I)=EV(2,2) 
         EKM(3,2,I)=EV(3,2) 
         EKM(1,3,I)=EV(1,3) 
         EKM(2,3,I)=EV(2,3) 
         EKM(3,3,I)=EV(3,3) 
  148 CONTINUE 
      DO I=1,6*N
         DO J=1,6*N
            HESS(J,I)=0.0D0 
         ENDDO
      ENDDO
      GTEST=.TRUE.
      SSTEST=.TRUE.
      CALL DRVTV(NMOLS,GTEST,SSTEST,
     1           FX,FY,FZ,FA,FB,FC,FX0,FY0,FZ0,FX1,FY1,FZ1,FX2,FY2,FZ2,FX3,FY3,FZ3,FX4,FY4,FZ4,
     2           X,Y,Z,EA,EB,EC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,XX3,YY3,ZZ3,XX4,YY4,ZZ4,
     3           AX1,AX2,AX3,AX4,
     3           AY1,AY2,AY3,AY4,AZ1,AZ2,AZ3,AZ4,BX1,BX2,BX3,BX4,BY1,BY2,BY3,BY4,
     4           CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,CZ1,CZ2,CZ3,CZ4,
     5           AAX1,AAX2,AAX3,AAX4,ABX1,ABX2,ABX3,ABX4,ACX1,ACX2,ACX3,ACX4,BBX1,BBX2,BBX3,BBX4,BCX1,BCX2,BCX3,BCX4,CCX1,
     6           CCX2,CCX3,CCX4,AAY1,AAY2,AAY3,AAY4,ABY1,ABY2,ABY3,ABY4,ACY1,ACY2,ACY3,ACY4,BBY1,BBY2,BBY3,BBY4,
     7           BCY1,BCY2,BCY3,BCY4,CCY1,CCY2,CCY3,CCY4,AAZ1,AAZ2,AAZ3,AAZ4,ACZ1,ACZ2,ACZ3,ACZ4,CCZ1,CCZ2,CCZ3,CCZ4,
     8           LC)

      DO I=1,N 
         I1=3*I-2 
         I2=I1+1 
         I3=I2+1 
         K1=N3+I1 
         K2=K1+1 
         K3=K2+1 
         DO J=1,N 
            J1=3*J-2 
            J2=J1+1 
            J3=J2+1 
            L1=N3+J1 
            L2=L1+1 
            L3=L2+1 
            T(1,1)=HESS(L1,K1) 
            T(2,1)=HESS(L2,K1) 
            T(3,1)=HESS(L3,K1) 
            T(1,2)=HESS(L1,K2) 
            T(2,2)=HESS(L2,K2) 
            T(3,2)=HESS(L3,K2) 
            T(1,3)=HESS(L1,K3) 
            T(2,3)=HESS(L2,K3) 
            T(3,3)=HESS(L3,K3) 
            DO K=1,3 
               DO L=1,3 
                  SS=0.0D0 
                  DO M=1,3 
                     SS=SS+T(L,M)*EKM(M,K,I) 
                  ENDDO
                  S(L,K)=SS 
               ENDDO
            ENDDO
            DO K=1,3 
               DO L=1,3 
                  SS=0.0D0 
                  DO M=1,3 
                     SS=SS+S(M,K)*EKM(M,L,J) 
                  ENDDO
                  T(L,K)=SS 
               ENDDO
            ENDDO
            HESS(L1,K1)=T(1,1)*EKV(J1)*EKV(I1) 
            HESS(L2,K1)=T(2,1)*EKV(J2)*EKV(I1) 
            HESS(L3,K1)=T(3,1)*EKV(J3)*EKV(I1) 
            HESS(L1,K2)=T(1,2)*EKV(J1)*EKV(I2) 
            HESS(L2,K2)=T(2,2)*EKV(J2)*EKV(I2) 
            HESS(L3,K2)=T(3,2)*EKV(J3)*EKV(I2) 
            HESS(L1,K3)=T(1,3)*EKV(J1)*EKV(I3) 
            HESS(L2,K3)=T(2,3)*EKV(J2)*EKV(I3) 
            HESS(L3,K3)=T(3,3)*EKV(J3)*EKV(I3) 
         ENDDO
      ENDDO
      DO I=1,N 
         I1=3*I-2 
         I2=I1+1 
         I3=I2+1 
         DO J=1,N 
            J1=3*J-2 
            J2=J1+1 
            J3=J2+1 
            K1=N3+J1 
            K2=K1+1 
            K3=K2+1 
            T(1,1)=HESS(K1,I1) 
            T(2,1)=HESS(K2,I1) 
            T(3,1)=HESS(K3,I1) 
            T(1,2)=HESS(K1,I2) 
            T(2,2)=HESS(K2,I2) 
            T(3,2)=HESS(K3,I2) 
            T(1,3)=HESS(K1,I3) 
            T(2,3)=HESS(K2,I3) 
            T(3,3)=HESS(K3,I3) 
            DO K=1,3 
               DO L=1,3 
                  SS=0.0D0 
                  DO M=1,3 
                     SS=SS+EKM(M,L,J)*T(M,K) 
                  ENDDO
                  S(L,K)=SS 
               ENDDO
            ENDDO
            HESS(K1,I1)=S(1,1)*EKV(J1) 
            HESS(K2,I1)=S(2,1)*EKV(J2) 
            HESS(K3,I1)=S(3,1)*EKV(J3) 
            HESS(K1,I2)=S(1,2)*EKV(J1) 
            HESS(K2,I2)=S(2,2)*EKV(J2) 
            HESS(K3,I2)=S(3,2)*EKV(J3) 
            HESS(K1,I3)=S(1,3)*EKV(J1) 
            HESS(K2,I3)=S(2,3)*EKV(J2) 
            HESS(K3,I3)=S(3,3)*EKV(J3) 
         ENDDO
      ENDDO
      DO I=1,N6 
         DO J=1,I 
            HESS(J,I)=HESS(I,J) 
         ENDDO
      ENDDO

      CALL DSYEV('V','U',N6,HESS,SIZE(HESS,1),DIAG,TEMPA,18*NMOLS,INFO)
      call eigensort_val_asc(diag,hess,n6,6*nmols)

      DO 772 I=1,N6 
         IF(DIAG(I).GT.0.0D0) THEN 
            DIAG(I)=SQRT(DIAG(I))*SW 
         ELSE 
            DIAG(I)=-SQRT(-DIAG(I))*SW 
         ENDIF 
772   CONTINUE 

      RETURN
      END 
