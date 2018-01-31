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
      SUBROUTINE H2O(NMOLS,IPOT,COORDS,GRAD,ENERGY,GTEST,SSTEST)
C ----------------------------------------------------------------------
C 
C       NORMAL MODE ANALYSIS FOR WATER MOLECULES 
C       TO THREE TRANSLATIONAL MODES AND THREE ROTATINAL MODES 
C                          CODED BY H. TANAKA      1988, 01, 08 
C                               SEE J.CHEM.PHYS. 87, 6070 (1987) 
C                          Modified by me 1992
C       FIRST AND SECOND DERIVATIVES ARE EXPRESSED IN ANALYTICAL FORM 
C       EWALD SUM IS NOT TAKEN INTO ACCOUNT 
C       Writing to unit 13 suppressed
C
C       This is the subroutine designed to return the energy and
C       first and second derivatives with respect to the centre
C       of mass and the Euler angles.
C
C ----------------------------------------------------------------------
      USE MODHESS
      IMPLICIT NONE
      LOGICAL GTEST, SSTEST
      INTEGER NMOLS, IPOT, N, NN, N3, N6, J1, J2, I, J
      DOUBLE PRECISION ENERGY, RC, BXL, RC2, RL, RQ, RQ3, RQ6, EP, EMUPID, P1, P2, P3,
     1                 P4, P5, P6, ANUM, QQ, AD1, AD2, ANGLE, HOLEN, COLEN, QE, UJ, SG,
     1                 WM, PI, QS, RANGLE, OHZ, HYL, HZL, OL, CL, ENU, PN, PNI, EMU, RDENS,
     1                 DENS, VOL, TH, PH, PS, SINA, SINB, SINC, COSA, COSB, COSC, EMUP
      DOUBLE PRECISION COORDS(6*NMOLS), GRAD(6*NMOLS)
      DOUBLE PRECISION FX(NMOLS),FY(NMOLS),FZ(NMOLS), 
     *               FA(NMOLS),FB(NMOLS),FC(NMOLS), 
     *               FX0(NMOLS),FY0(NMOLS),FZ0(NMOLS), 
     *               FX1(NMOLS),FY1(NMOLS),FZ1(NMOLS), 
     *               FX2(NMOLS),FY2(NMOLS),FZ2(NMOLS), 
     *               FX3(NMOLS),FY3(NMOLS),FZ3(NMOLS), 
     *               FX4(NMOLS),FY4(NMOLS),FZ4(NMOLS)
      DOUBLE PRECISION X(NMOLS+1),Y(NMOLS+1),Z(NMOLS+1), ! THE +1 IS SUPPOSED TO SPEED UP EXECUTION?!
     *            A(NMOLS),B(NMOLS),C(NMOLS),D(NMOLS), 
     *            EA(NMOLS+1),EB(NMOLS+1),EC(NMOLS+1), 
     *            XX1(NMOLS),YY1(NMOLS),ZZ1(NMOLS), 
     *            XX2(NMOLS),YY2(NMOLS),ZZ2(NMOLS), 
     *            XX3(NMOLS),YY3(NMOLS),ZZ3(NMOLS), 
     *            XX4(NMOLS),YY4(NMOLS),ZZ4(NMOLS) 
      DOUBLE PRECISION AX1(NMOLS),AX2(NMOLS),AX3(NMOLS),AX4(NMOLS), 
     *          AY1(NMOLS),AY2(NMOLS),AY3(NMOLS),AY4(NMOLS), 
     *          AZ1(NMOLS),AZ2(NMOLS),AZ3(NMOLS),AZ4(NMOLS), 
     *          BX1(NMOLS),BX2(NMOLS),BX3(NMOLS),BX4(NMOLS), 
     *          BY1(NMOLS),BY2(NMOLS),BY3(NMOLS),BY4(NMOLS), 
     *          CX1(NMOLS),CX2(NMOLS),CX3(NMOLS),CX4(NMOLS), 
     *          CY1(NMOLS),CY2(NMOLS),CY3(NMOLS),CY4(NMOLS), 
     *          CZ1(NMOLS),CZ2(NMOLS),CZ3(NMOLS),CZ4(NMOLS) 
      DOUBLE PRECISION AAX1(NMOLS),AAX2(NMOLS),AAX3(NMOLS),AAX4(NMOLS), 
     *          ABX1(NMOLS),ABX2(NMOLS),ABX3(NMOLS),ABX4(NMOLS), 
     *          ACX1(NMOLS),ACX2(NMOLS),ACX3(NMOLS),ACX4(NMOLS), 
     *          BBX1(NMOLS),BBX2(NMOLS),BBX3(NMOLS),BBX4(NMOLS), 
     *          BCX1(NMOLS),BCX2(NMOLS),BCX3(NMOLS),BCX4(NMOLS), 
     *          CCX1(NMOLS),CCX2(NMOLS),CCX3(NMOLS),CCX4(NMOLS) 
      DOUBLE PRECISION AAY1(NMOLS),AAY2(NMOLS),AAY3(NMOLS),AAY4(NMOLS), 
     *          ABY1(NMOLS),ABY2(NMOLS),ABY3(NMOLS),ABY4(NMOLS), 
     *          ACY1(NMOLS),ACY2(NMOLS),ACY3(NMOLS),ACY4(NMOLS), 
     *          BBY1(NMOLS),BBY2(NMOLS),BBY3(NMOLS),BBY4(NMOLS), 
     *          BCY1(NMOLS),BCY2(NMOLS),BCY3(NMOLS),BCY4(NMOLS), 
     *          CCY1(NMOLS),CCY2(NMOLS),CCY3(NMOLS),CCY4(NMOLS) 
      DOUBLE PRECISION AAZ1(NMOLS),AAZ2(NMOLS),AAZ3(NMOLS),AAZ4(NMOLS), 
C    *          ABZ1(NMOLS),ABZ2(NMOLS),ABZ3(NMOLS),ABZ4(NMOLS), 
     *          ACZ1(NMOLS),ACZ2(NMOLS),ACZ3(NMOLS),ACZ4(NMOLS), 
C    *          BBZ1(NMOLS),BBZ2(NMOLS),BBZ3(NMOLS),BBZ4(NMOLS), 
C    *          BCZ1(NMOLS),BCZ2(NMOLS),BCZ3(NMOLS),BCZ4(NMOLS), 
     *          CCZ1(NMOLS),CCZ2(NMOLS),CCZ3(NMOLS),CCZ4(NMOLS) 
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
      N=NMOLS
C
C  IPOT defines the TIPS type as usual.
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
      QS=QQ**2*QE 
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
      PN=FLOAT(N) 
      PNI=1.0D0/PN 
      EMU=1.0D-10*ENU*ANUM*PNI 
      EMUPID=EMU*PN 
      IF (N.EQ.64) THEN
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
C
C  COORDS contains the coordinates as X1, Y1, Z1, X2, Y2,...
C  ...,Z3N,E11,E12,E13, etc.
C
      DO 666 J1=1,N
         J2=3*(J1-1)
         X(J1)=COORDS(J2+1)
         Y(J1)=COORDS(J2+2)
         Z(J1)=COORDS(J2+3)
C        READ(18,*) EA(J1),EB(J1),EC(J1)
C        WRITE(*,*) J1,EA(J1),EB(J1),EC(J1)
         EA(J1)=COORDS(3*N+J2+1)
         EB(J1)=COORDS(3*N+J2+2)
         EC(J1)=COORDS(3*N+J2+3)
C        WRITE(*,'(I6,6F15.5)') J1,X(J1),Y(J1),Z(J1),EA(J1),EB(J1),EC(J1)
C        COORDS(3*N+J2+1)=EA(J1)
C        COORDS(3*N+J2+2)=EB(J1)
C        COORDS(3*N+J2+3)=EC(J1)
666   CONTINUE
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
         LC(I)=1 
         CALL LC1(NMOLS,SINA,SINB,SINC,COSA,COSB,COSC,HYL,HZL,OL,CL,I,
     1            X,Y,Z,A,B,C,D,EA,EB,EC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,XX3,YY3,ZZ3,XX4,YY4,ZZ4,
     2            AX1,AX2,AX3,AX4,AY1,AY2,AY3,AY4,AZ1,AZ2,AZ3,AZ4,BX1,BX2,BX3,BX4,BY1,BY2,BY3,BY4,CX1,CX2,CX3,CX4,
     3            CY1,CY2,CY3,CY4,CZ1,CZ2,CZ3,CZ4,AAX1,AAX2,AAX3,AAX4,ABX1,ABX2,ABX3,ABX4,ACX1,ACX2,ACX3,ACX4,BBX1,BBX2,
     4            BBX3,BBX4,BCX1,BCX2,BCX3,BCX4,CCX1,CCX2,CCX3,CCX4,AAY1,AAY2,AAY3,AAY4,ABY1,ABY2,ABY3,ABY4,
     5            ACY1,ACY2,ACY3,ACY4,BBY1,BBY2,BBY3,BBY4,BCY1,BCY2,BCY3,BCY4,CCY1,CCY2,CCY3,CCY4,
     6            AAZ1,AAZ2,AAZ3,AAZ4,ACZ1,ACZ2,ACZ3,ACZ4,CCZ1,CCZ2,CCZ3,CCZ4,LC) 
48    CONTINUE 
      DO I=1,N
         FA(I)=0.0D0
         FB(I)=0.0D0
         FC(I)=0.0D0
         FX(I)=0.0D0
         FY(I)=0.0D0
         FZ(I)=0.0D0
      ENDDO
      IF (SSTEST) THEN
         DO I=1,6*N ! WCOMMENT
            DO J=1,6*N
               HESS(J,I)=0.0D0 
            ENDDO
         ENDDO
      ENDIF
C
C  DRVTV does most of the work:
C    EP is the energy in kJ per mole / 2
C    SD are the second derivatives
C
C    N3=3*N, so I1(J1) = 1, 4, 7,...
C               I2(J2) = 2, 5, 8,...
C               I3(J3) = 3, 6, 9,...
C               K1(L1) = 3*N + I1(J1)
C               K2(L2) = 3*N + I2(J2)
C               K3(L3) = 3*N + I3(J3)
C
C  The second derivative matrix therefore appears to be in exactly
C  the form that I need it, with all the centre of mass derivatives
C  and then all the Euler angle derivatives.
C
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
C
C  EMUP is an energy conversion factor that puts everything
C  in kJ per mol per molecule
C
      EMUP=49.72553808D0
      DO 754 J1=1,N
         J2=3*(J1-1)
         GRAD(J2+1)=FX(J1)*EMUP
         GRAD(J2+2)=FY(J1)*EMUP
         GRAD(J2+3)=FZ(J1)*EMUP
         GRAD(3*N+J2+1)=FA(J1)*EMUP
         GRAD(3*N+J2+2)=FB(J1)*EMUP
         GRAD(3*N+J2+3)=FC(J1)*EMUP
754   CONTINUE
C     EP=EP/N 
      EP=EP
      ENERGY=EP 
      IF (.NOT.SSTEST) RETURN
      DO 113 J1=1,6*N
         DO 112 J2=1,6*N
            HESS(J2,J1)=HESS(J2,J1)*EMUP
C           WRITE(*,'(A,2I6,F20.10)') 'J1,J2,HESS=',J1,J2,HESS(J2,J1)
112      CONTINUE 
113   CONTINUE

      RETURN
      END 
C
C********************************************************************
C
      SUBROUTINE LC1(NMOLS,SINA,SINB,SINC,COSA,COSB,COSC,HYL,HZL,OL,CL,I,
     1            X,Y,Z,A,B,C,D,EA,EB,EC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,XX3,YY3,ZZ3,XX4,YY4,ZZ4,
     2            AX1,AX2,AX3,AX4,AY1,AY2,AY3,AY4,AZ1,AZ2,AZ3,AZ4,BX1,BX2,BX3,BX4,BY1,BY2,BY3,BY4,CX1,CX2,CX3,CX4,
     3            CY1,CY2,CY3,CY4,CZ1,CZ2,CZ3,CZ4,AAX1,AAX2,AAX3,AAX4,ABX1,ABX2,ABX3,ABX4,ACX1,ACX2,ACX3,ACX4,BBX1,BBX2,
     4            BBX3,BBX4,BCX1,BCX2,BCX3,BCX4,CCX1,CCX2,CCX3,CCX4,AAY1,AAY2,AAY3,AAY4,ABY1,ABY2,ABY3,ABY4,
     5            ACY1,ACY2,ACY3,ACY4,BBY1,BBY2,BBY3,BBY4,BCY1,BCY2,BCY3,BCY4,CCY1,CCY2,CCY3,CCY4,
     6            AAZ1,AAZ2,AAZ3,AAZ4,ACZ1,ACZ2,ACZ3,ACZ4,CCZ1,CCZ2,CCZ3,CCZ4,LC) 
      IMPLICIT NONE
      INTEGER N, I, NN, N3, N6, NMOLS
      DOUBLE PRECISION SINA, SINB, SINC, COSA, COSB, COSC, HYL, HZL, OL, CL, RC, BXL,
     1  RC2, RL, RQ, RQ3, RQ6, EP, EMUPID, P1, P2, P3, P4, P5, P6, SP12,
     1  DDA12, DDB12, DDC12, DAA12, DAB12, DAC12, DBB12, DBC12, DCC12, SP13,
     1  DDA13, DDB13, DAA13, DAB13, DBB13, SP22, DDA22, DDB22, DDC22, DAA22,
     1  DAB22, DAC22, DBB22, DBC22, DCC22, SP23, DDA23, DDB23, DAA23, DAB23,
     1  DBB23, SP32, DDA32, DDC32, DAA32, DAC32, DCC32, SP33, DDA33, DAA33, SY12,
     1  SY22, SY32, SZ13, SZ23, SZ33, A12H, A22H, A32H, A13H, A23H, A33H, B12H,
     1  B22H, B13H, B23H, C12H, C22H, C32H, A13C, A23C, A33C, B13C, B23C, A13O,
     1  A23O, A33O, B13O, B23O, AA12H, AA22H, AA32H, AA13H, AA23H, AA33H, AB12H,
     1  AB22H, AB13H, AB23H, AC12H, AC22H, AC32H, BB12H, BB22H, BB13H, BB23H, BC12H,
     1  BC22H, CC12H, CC22H, CC32H, AA13C, AA23C, AA33C, AB13C, AB23C, BB13C, BB23C,
     1  AA13O, AA23O, AA33O, AB13O, AB23O, BB13O, BB23O
      DOUBLE PRECISION X(NMOLS+1),Y(NMOLS+1),Z(NMOLS+1), 
     *            A(NMOLS),B(NMOLS),C(NMOLS),D(NMOLS), 
     *            EA(NMOLS+1),EB(NMOLS+1),EC(NMOLS+1), 
     *            XX1(NMOLS),YY1(NMOLS),ZZ1(NMOLS), 
     *            XX2(NMOLS),YY2(NMOLS),ZZ2(NMOLS), 
     *            XX3(NMOLS),YY3(NMOLS),ZZ3(NMOLS), 
     *            XX4(NMOLS),YY4(NMOLS),ZZ4(NMOLS) 
      DOUBLE PRECISION AX1(NMOLS),AX2(NMOLS),AX3(NMOLS),AX4(NMOLS), 
     *          AY1(NMOLS),AY2(NMOLS),AY3(NMOLS),AY4(NMOLS), 
     *          AZ1(NMOLS),AZ2(NMOLS),AZ3(NMOLS),AZ4(NMOLS), 
     *          BX1(NMOLS),BX2(NMOLS),BX3(NMOLS),BX4(NMOLS), 
     *          BY1(NMOLS),BY2(NMOLS),BY3(NMOLS),BY4(NMOLS), 
     *          CX1(NMOLS),CX2(NMOLS),CX3(NMOLS),CX4(NMOLS), 
     *          CY1(NMOLS),CY2(NMOLS),CY3(NMOLS),CY4(NMOLS), 
     *          CZ1(NMOLS),CZ2(NMOLS),CZ3(NMOLS),CZ4(NMOLS) 
      DOUBLE PRECISION AAX1(NMOLS),AAX2(NMOLS),AAX3(NMOLS),AAX4(NMOLS), 
     *          ABX1(NMOLS),ABX2(NMOLS),ABX3(NMOLS),ABX4(NMOLS), 
     *          ACX1(NMOLS),ACX2(NMOLS),ACX3(NMOLS),ACX4(NMOLS), 
     *          BBX1(NMOLS),BBX2(NMOLS),BBX3(NMOLS),BBX4(NMOLS), 
     *          BCX1(NMOLS),BCX2(NMOLS),BCX3(NMOLS),BCX4(NMOLS), 
     *          CCX1(NMOLS),CCX2(NMOLS),CCX3(NMOLS),CCX4(NMOLS) 
      DOUBLE PRECISION AAY1(NMOLS),AAY2(NMOLS),AAY3(NMOLS),AAY4(NMOLS), 
     *          ABY1(NMOLS),ABY2(NMOLS),ABY3(NMOLS),ABY4(NMOLS), 
     *          ACY1(NMOLS),ACY2(NMOLS),ACY3(NMOLS),ACY4(NMOLS), 
     *          BBY1(NMOLS),BBY2(NMOLS),BBY3(NMOLS),BBY4(NMOLS), 
     *          BCY1(NMOLS),BCY2(NMOLS),BCY3(NMOLS),BCY4(NMOLS), 
     *          CCY1(NMOLS),CCY2(NMOLS),CCY3(NMOLS),CCY4(NMOLS) 
      DOUBLE PRECISION AAZ1(NMOLS),AAZ2(NMOLS),AAZ3(NMOLS),AAZ4(NMOLS), 
     *          ACZ1(NMOLS),ACZ2(NMOLS),ACZ3(NMOLS),ACZ4(NMOLS), 
     *          CCZ1(NMOLS),CCZ2(NMOLS),CCZ3(NMOLS),CCZ4(NMOLS) 
      INTEGER LC(NMOLS)
      COMMON /IN/ N,NN,N3,N6 
      COMMON /DR/ RC,BXL,RC2,RL,RQ,RQ3,RQ6,EP,EMUPID,P1,P2,P3,P4,P5,P6 
      SP12=-(SINC*COSB+COSA*SINB*COSC) 
      DDA12=SINA*SINB*COSC 
      DDB12=SINC*SINB-COSA*COSB*COSC 
      DDC12=-COSC*COSB+COSA*SINB*SINC 
      DAA12=COSA*SINB*COSC 
      DAB12=SINA*COSB*COSC 
      DAC12=-SINA*SINB*SINC 
      DBB12=-SP12 
      DBC12=COSC*SINB+COSA*COSB*SINC 
      DCC12=-SP12 
      SP13=SINA*SINB 
      DDA13=COSA*SINB 
      DDB13=SINA*COSB 
      DAA13=-SP13 
      DAB13=COSA*COSB 
      DBB13=-SP13 
      SP22=-SINC*SINB+COSA*COSB*COSC 
      DDA22=-SINA*COSB*COSC 
      DDB22=-SINC*COSB-COSA*SINB*COSC 
      DDC22=-COSC*SINB-COSA*COSB*SINC 
      DAA22=-COSA*COSB*COSC 
      DAB22=SINA*SINB*COSC 
      DAC22=SINA*COSB*SINC 
      DBB22=-SP22 
      DBC22=-COSC*COSB+COSA*SINB*SINC 
      DCC22=-SP22 
      SP23=-SINA*COSB 
      DDA23=-COSA*COSB 
      DDB23=SINA*SINB 
      DAA23=-SP23 
      DAB23=COSA*SINB 
      DBB23=-SP23 
      SP32=SINA*COSC 
      DDA32=COSA*COSC 
      DDC32=-SINA*SINC 
      DAA32=-SP32 
      DAC32=-COSA*SINC 
      DCC32=-SP32 
      SP33=COSA 
      DDA33=-SINA 
      DAA33=-COSA 
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
      A13C=DDA13*CL 
      A23C=DDA23*CL 
      A33C=DDA33*CL 
      B13C=DDB13*CL 
      B23C=DDB23*CL 
      A13O=DDA13*OL 
      A23O=DDA23*OL 
      A33O=DDA33*OL 
      B13O=DDB13*OL 
      B23O=DDB23*OL 
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
      BB12H=DBB12*HYL 
      BB22H=DBB22*HYL 
      BB13H=DBB13*HZL 
      BB23H=DBB23*HZL 
      BC12H=DBC12*HYL 
      BC22H=DBC22*HYL 
      CC12H=DCC12*HYL 
      CC22H=DCC22*HYL 
      CC32H=DCC32*HYL 
      AA13C=DAA13*CL 
      AA23C=DAA23*CL 
      AA33C=DAA33*CL 
      AB13C=DAB13*CL 
      AB23C=DAB23*CL 
      BB13C=DBB13*CL 
      BB23C=DBB23*CL 
      AA13O=DAA13*OL 
      AA23O=DAA23*OL 
      AA33O=DAA33*OL 
      AB13O=DAB13*OL 
      AB23O=DAB23*OL 
      BB13O=DBB13*OL 
      BB23O=DBB23*OL 
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
      CX1(I)=C12H 
      CY1(I)=C22H 
      CZ1(I)=C32H 
      CX2(I)=-C12H 
      CY2(I)=-C22H 
      CZ2(I)=-C32H 
      CX3(I)=0.0D0 
      CY3(I)=0.0D0 
      CZ3(I)=0.0D0 
      CX4(I)=0.0D0 
      CY4(I)=0.0D0 
      CZ4(I)=0.0D0 
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
      CCX1(I)=CC12H 
      CCY1(I)=CC22H 
      CCZ1(I)=CC32H 
      CCX2(I)=-CC12H 
      CCY2(I)=-CC22H 
      CCZ2(I)=-CC32H 
      CCX3(I)=0.0D0 
      CCY3(I)=0.0D0 
      CCZ3(I)=0.0D0 
      CCX4(I)=0.0D0 
      CCY4(I)=0.0D0 
      CCZ4(I)=0.0D0 
      ABX1(I)=AB12H+AB13H 
      ABY1(I)=AB22H+AB23H 
      ABX2(I)=-AB12H+AB13H 
      ABY2(I)=-AB22H+AB23H 
      ABX3(I)=AB13O 
      ABY3(I)=AB23O 
      ABX4(I)=AB13C 
      ABY4(I)=AB23C 
      ACX1(I)=AC12H 
      ACY1(I)=AC22H 
      ACZ1(I)=AC32H 
      ACX2(I)=-AC12H 
      ACY2(I)=-AC22H 
      ACZ2(I)=-AC32H 
      ACX3(I)=0.0D0 
      ACY3(I)=0.0D0 
      ACZ3(I)=0.0D0 
      ACX4(I)=0.0D0 
      ACY4(I)=0.0D0 
      ACZ4(I)=0.0D0 
      BCX1(I)=BC12H 
      BCY1(I)=BC22H 
      BCX2(I)=-BC12H 
      BCY2(I)=-BC22H 
      BCX3(I)=0.0D0 
      BCY3(I)=0.0D0 
      BCX4(I)=0.0D0 
      BCY4(I)=0.0D0 
      RETURN 
      END 
C
C******************************************************************
C
C  This is where the derivatives are calculated.
C  The total EP = sum V(i,j) / 2 kJ per mole
C            FDX is the gradient with respect to the X coordinate
C                 of the centre of mass on a given site
C            FDA is the derivative of the energy with respect to the
C                 first Euler angle on a given site
C            SDXX, SDXA, SDAA are the analogous second derivatives
C
      SUBROUTINE DRVTV(NMOLS,GTEST,SSTEST,
     1                 FX,FY,FZ,FA,FB,FC,FX0,FY0,FZ0,FX1,FY1,FZ1,FX2,FY2,FZ2,FX3,FY3,FZ3,FX4,FY4,FZ4,
     2                 X,Y,Z,EA,EB,EC,XX1,YY1,ZZ1,XX2,YY2,ZZ2,XX3,YY3,ZZ3,XX4,YY4,ZZ4,
     3                 AX1,AX2,AX3,AX4,
     3                 AY1,AY2,AY3,AY4,AZ1,AZ2,AZ3,AZ4,BX1,BX2,BX3,BX4,BY1,BY2,BY3,BY4,
     4                 CX1,CX2,CX3,CX4,CY1,CY2,CY3,CY4,CZ1,CZ2,CZ3,CZ4,
     5                 AAX1,AAX2,AAX3,AAX4,ABX1,ABX2,ABX3,ABX4,ACX1,ACX2,ACX3,ACX4,BBX1,BBX2,BBX3,BBX4,BCX1,BCX2,BCX3,BCX4,CCX1,
     6                 CCX2,CCX3,CCX4,AAY1,AAY2,AAY3,AAY4,ABY1,ABY2,ABY3,ABY4,ACY1,ACY2,ACY3,ACY4,BBY1,BBY2,BBY3,BBY4,
     7                 BCY1,BCY2,BCY3,BCY4,CCY1,CCY2,CCY3,CCY4,AAZ1,AAZ2,AAZ3,AAZ4,ACZ1,ACZ2,ACZ3,ACZ4,CCZ1,CCZ2,CCZ3,CCZ4,
     8                 LC)
      USE MODHESS
      IMPLICIT NONE
      LOGICAL GTEST,SSTEST
      INTEGER NMOLS, N, NN, N3, N6, I, NTI1, NRI1, NRI2, NRI3, NTI2, NTI3, J, NTJ1, NRJ1, NRJ2, 
     1        NRJ3, NTJ2, NTJ3
      DOUBLE PRECISION RC, BXL, RC2, RL, RQ, RQ3, RQ6, EP, EMUPID, P1, P2, P3, P4, P5, P6,
     1 XX, YY, ZZ, XL1, YL1, ZL1, XL2, YL2, ZL2, XL3, YL3, ZL3, XL4, YL4,
     1 ZL4, DX, DY, DZ, SX, SY, SZ, R2, XK1, YK1, ZK1, XK2, YK2, ZK2, XK3, YK3,
     1 ZK3, XK4, YK4, ZK4, DX1, DY1, DZ1, DX2, DY2, DZ2, DX3, DY3, DZ3, DX4, DY4,
     1 DZ4, RR, RSI, RI, RRL, RRS, RRL2, RRS2, SF, DF, DS, XY, YZ, ZX, DFRI, SDX,
     1 SDY, SDZ, DSF, DSFRSI, SXX, SYY, SZZ, SXY, SYZ, SZX, RSIM, RTI, RFI, RTIM, DDX, DDY, DDZ,
     1 RAI, DDAI, RAJ, DDAJ, RBI, DDBI, RBJ, DDBJ, RCI, DDCI, RCJ, DDCJ, RFI3, DXX,
     1 DYY, DZZ, DXY, DYZ, DZX, RI3, RI3X, RI3Y, RI3Z, DXAI, DXAJ, DXBI, DXBJ, DXCI,
     1 DXCJ, DYAI, DYAJ, DYBI, DYBJ, DYCI, DYCJ, DZAI, DZAJ, DZBI, DZBJ, DZCI, DZCJ, RTI3,
     1 DAAI, DAAJ, DABI, DABJ, DACI, DACJ, DBBI, DBBJ, DBCI, DBCJ, DCCI, DCCJ, DAAK, DABK, DBAK,
     1 DACK, DCAK, DBBK, DBCK, DCBK, DCCK, FDX, FDY, FDZ, FDAI, FDBI, FDCI, FDAJ, FDBJ, FDCJ,
     1 SDXX, SDYY, SDZZ, SDXY, SDYZ, SDZX, SDXAI, SDYAI, SDZAI, SDXBI, SDYBI, SDZBI, SDXCI, SDYCI,
     1 SDZCI, SDXAJ, SDYAJ, SDZAJ, SDXBJ, SDYBJ, SDZBJ, SDXCJ, SDYCJ, SDZCJ, SDAAI, SDABI, SDACI,
     1 SDBBI, SDBCI, SDCCI, SDAAJ, SDABJ, SDACJ, SDBBJ, SDBCJ, SDCCJ, SDAAK, SDABK, SDBAK, SDACK,
     1 SDCAK, SDBBK, SDBCK, SDCBK, SDCCK, EPI, RHI, RDI, DE, DSFRI, DXDSF, DYDSF, DZDSF,
     1 U, VDX, VDY, VDZ, EPID
      DOUBLE PRECISION FX(NMOLS),FY(NMOLS),FZ(NMOLS), 
     *               FA(NMOLS),FB(NMOLS),FC(NMOLS), 
     *               FX0(NMOLS),FY0(NMOLS),FZ0(NMOLS), 
     *               FX1(NMOLS),FY1(NMOLS),FZ1(NMOLS), 
     *               FX2(NMOLS),FY2(NMOLS),FZ2(NMOLS), 
     *               FX3(NMOLS),FY3(NMOLS),FZ3(NMOLS), 
     *               FX4(NMOLS),FY4(NMOLS),FZ4(NMOLS)
      DOUBLE PRECISION X(NMOLS+1),Y(NMOLS+1),Z(NMOLS+1), 
     *            A(NMOLS),B(NMOLS),C(NMOLS),D(NMOLS), 
     *            EA(NMOLS+1),EB(NMOLS+1),EC(NMOLS+1), 
     *            XX1(NMOLS),YY1(NMOLS),ZZ1(NMOLS), 
     *            XX2(NMOLS),YY2(NMOLS),ZZ2(NMOLS), 
     *            XX3(NMOLS),YY3(NMOLS),ZZ3(NMOLS), 
     *            XX4(NMOLS),YY4(NMOLS),ZZ4(NMOLS) 
      DOUBLE PRECISION AX1(NMOLS),AX2(NMOLS),AX3(NMOLS),AX4(NMOLS), 
     *          AY1(NMOLS),AY2(NMOLS),AY3(NMOLS),AY4(NMOLS), 
     *          AZ1(NMOLS),AZ2(NMOLS),AZ3(NMOLS),AZ4(NMOLS), 
     *          BX1(NMOLS),BX2(NMOLS),BX3(NMOLS),BX4(NMOLS), 
     *          BY1(NMOLS),BY2(NMOLS),BY3(NMOLS),BY4(NMOLS), 
     *          CX1(NMOLS),CX2(NMOLS),CX3(NMOLS),CX4(NMOLS), 
     *          CY1(NMOLS),CY2(NMOLS),CY3(NMOLS),CY4(NMOLS), 
     *          CZ1(NMOLS),CZ2(NMOLS),CZ3(NMOLS),CZ4(NMOLS) 
      DOUBLE PRECISION AAX1(NMOLS),AAX2(NMOLS),AAX3(NMOLS),AAX4(NMOLS), 
     *          ABX1(NMOLS),ABX2(NMOLS),ABX3(NMOLS),ABX4(NMOLS), 
     *          ACX1(NMOLS),ACX2(NMOLS),ACX3(NMOLS),ACX4(NMOLS), 
     *          BBX1(NMOLS),BBX2(NMOLS),BBX3(NMOLS),BBX4(NMOLS), 
     *          BCX1(NMOLS),BCX2(NMOLS),BCX3(NMOLS),BCX4(NMOLS), 
     *          CCX1(NMOLS),CCX2(NMOLS),CCX3(NMOLS),CCX4(NMOLS) 
      DOUBLE PRECISION AAY1(NMOLS),AAY2(NMOLS),AAY3(NMOLS),AAY4(NMOLS), 
     *          ABY1(NMOLS),ABY2(NMOLS),ABY3(NMOLS),ABY4(NMOLS), 
     *          ACY1(NMOLS),ACY2(NMOLS),ACY3(NMOLS),ACY4(NMOLS), 
     *          BBY1(NMOLS),BBY2(NMOLS),BBY3(NMOLS),BBY4(NMOLS), 
     *          BCY1(NMOLS),BCY2(NMOLS),BCY3(NMOLS),BCY4(NMOLS), 
     *          CCY1(NMOLS),CCY2(NMOLS),CCY3(NMOLS),CCY4(NMOLS) 
      DOUBLE PRECISION AAZ1(NMOLS),AAZ2(NMOLS),AAZ3(NMOLS),AAZ4(NMOLS), 
     *          ACZ1(NMOLS),ACZ2(NMOLS),ACZ3(NMOLS),ACZ4(NMOLS), 
     *          CCZ1(NMOLS),CCZ2(NMOLS),CCZ3(NMOLS),CCZ4(NMOLS) 
      INTEGER LC(NMOLS)
      COMMON /IN/ N,NN,N3,N6 
      COMMON /DR/ RC,BXL,RC2,RL,RQ,RQ3,RQ6,EP,EMUPID,P1,P2,P3,P4,P5,P6 
      EP=0.0D0 
C
C  Loop over all water molecules, except the last
C
      DO 60 I=1,NN 
      XX=X(I) 
      YY=Y(I) 
      ZZ=Z(I) 
      XL1=XX1(I) 
      YL1=YY1(I) 
      ZL1=ZZ1(I) 
      XL2=XX2(I) 
      YL2=YY2(I) 
      ZL2=ZZ2(I) 
      XL3=XX3(I) 
      YL3=YY3(I) 
      ZL3=ZZ3(I) 
      XL4=XX4(I) 
      YL4=YY4(I) 
      ZL4=ZZ4(I) 
C
C        T for translation R for rotation
C        I index
C  NTI1 = 1, 4, 7,...
C  NTI2 = 2, 5, 8,...
C  NTI3 = 3, 6, 9,...
C  NRI1 = 3*N + 1, 4, 7,...
C  NRI2 = 3*N + 2, 5, 8,...
C  NRI3 = 3*N + 3, 6, 9,...
C        J index
C  NTJ1 = 1, 4, 7,...
C  NTJ2 = 2, 5, 8,...
C  NTJ3 = 3, 6, 9,...
C  NRJ1 = 3*N + 1, 4, 7,...
C  NRJ2 = 3*N + 2, 5, 8,...
C  NRJ3 = 3*N + 3, 6, 9,...
C
      NTI1=3*I-2 
      NRI1=N3+NTI1 
      NRI2=NRI1+1 
      NRI3=NRI2+1 
      NTI2=NTI1+1 
      NTI3=NTI2+1 
C
C  Loop over all the other water molecules
C
*VOPTION NOVEC 
      DO 70 J=I+1,N 
      DX=XX-X(J) 
      IF (N.EQ.64) THEN
         IF(ABS(DX).LE.RC) GOTO 35 
         DX=DX-SIGN(BXL,DX) 
      ENDIF
   35 DY=YY-Y(J) 
      IF (N.EQ.64) THEN
         IF(ABS(DY).LE.RC) GOTO 45 
         DY=DY-SIGN(BXL,DY) 
      ENDIF
   45 DZ=ZZ-Z(J) 
      IF (N.EQ.64) THEN
         IF(ABS(DZ).LE.RC) GOTO 55 
         DZ=DZ-SIGN(BXL,DZ) 
      ENDIF
   55 SX=DX*DX 
      SY=DY*DY 
      SZ=DZ*DZ 
      R2=SX+SY+SZ 
      IF(R2.GT.RC2) GOTO 70 
      XK1=XX1(J) 
      YK1=YY1(J) 
      ZK1=ZZ1(J) 
      XK2=XX2(J) 
      YK2=YY2(J) 
      ZK2=ZZ2(J) 
      XK3=XX3(J) 
      YK3=YY3(J) 
      ZK3=ZZ3(J) 
      XK4=XX4(J) 
      YK4=YY4(J) 
      ZK4=ZZ4(J) 
      DX1=DX+XL1 
      DY1=DY+YL1 
      DZ1=DZ+ZL1 
      DX2=DX+XL2 
      DY2=DY+YL2 
      DZ2=DZ+ZL2 
      DX3=DX+XL3 
      DY3=DY+YL3 
      DZ3=DZ+ZL3 
      DX4=DX+XL4 
      DY4=DY+YL4 
      DZ4=DZ+ZL4 
      RR=SQRT(R2) 
      IF(RR.GT.RL) THEN 
      RSI=1.0D0/R2 
      RI=1.0D0/RR 
      RRL=RR-RC 
      RRS=RR-RL 
      RRL2=RRL*RRL 
      RRS2=RRS*RRS 
      SF=RRL2*RRL*RQ*(10.0D0*RRS2-5.0D0*RRL*RRS+RRL2) 
      DF=RQ3*RRL2*RRS2 
      DS=RQ6*RRL*RRS*(RRL+RRS) 
      XY=DX*DY 
      YZ=DY*DZ 
      ZX=DZ*DX 
      DFRI=DF*RI 
      SDX=DFRI*DX 
      SDY=DFRI*DY 
      SDZ=DFRI*DZ 
      DSF=DS-DFRI 
      DSFRSI=DSF*RSI 
      SXX=DFRI+DSFRSI*SX 
      SYY=DFRI+DSFRSI*SY 
      SZZ=DFRI+DSFRSI*SZ 
      SXY=DSFRSI*XY 
      SYZ=DSFRSI*YZ 
      SZX=DSFRSI*ZX 
      END IF 
C ::::::: FIRST AND SECOND DERIVATIVES FOR EACH SITE 
      DX=DX1-XK1 
      DY=DY1-YK1 
      DZ=DZ1-ZK1 
      SX=DX*DX 
      SY=DY*DY 
      SZ=DZ*DZ 
      RSI=1.0D0/(SX+SY+SZ) 
      RI=SQRT(RSI) 
      IF (GTEST.OR.SSTEST) THEN
      RSIM=-RSI 
      RTI=RI*RSI 
      RFI=RSI*RTI 
      XY=DX*DY 
      YZ=DY*DZ 
      ZX=DZ*DX 
      RTIM=-RTI 
      DDX=RTIM*DX 
      DDY=RTIM*DY 
      DDZ=RTIM*DZ 
      RAI=RI*(AX1(I)*DX+AY1(I)*DY+AZ1(I)*DZ) 
      DDAI=RSIM*RAI 
      RAJ=RI*(AX1(J)*DX+AY1(J)*DY+AZ1(J)*DZ) 
      DDAJ=RSIM*RAJ 
      RBI=RI*(BX1(I)*DX+BY1(I)*DY) 
      DDBI=RSIM*RBI 
      RBJ=RI*(BX1(J)*DX+BY1(J)*DY) 
      DDBJ=RSIM*RBJ 
      RCI=RI*(CX1(I)*DX+CY1(I)*DY+CZ1(I)*DZ) 
      DDCI=RSIM*RCI 
      RCJ=RI*(CX1(J)*DX+CY1(J)*DY+CZ1(J)*DZ) 
      DDCJ=RSIM*RCJ 
      RFI3=RFI*3.0D0 
      DXX=RFI3*SX+RTIM 
      DYY=RFI3*SY+RTIM 
      DZZ=RFI3*SZ+RTIM 
      DXY=RFI3*XY 
      DYZ=RFI3*YZ 
      DZX=RFI3*ZX 
      RI3=3.0D0*RI 
      RI3X=DX*RI3 
      RI3Y=DY*RI3 
      RI3Z=DZ*RI3 
      ENDIF
      IF (SSTEST) THEN
      DXAI=(RI3X*RAI-AX1(I))*RTI 
      DXAJ=(RI3X*RAJ-AX1(J))*RTI 
      DXBI=(RI3X*RBI-BX1(I))*RTI 
      DXBJ=(RI3X*RBJ-BX1(J))*RTI 
      DXCI=(RI3X*RCI-CX1(I))*RTI 
      DXCJ=(RI3X*RCJ-CX1(J))*RTI 
      DYAI=(RI3Y*RAI-AY1(I))*RTI 
      DYAJ=(RI3Y*RAJ-AY1(J))*RTI 
      DYBI=(RI3Y*RBI-BY1(I))*RTI 
      DYBJ=(RI3Y*RBJ-BY1(J))*RTI 
      DYCI=(RI3Y*RCI-CY1(I))*RTI 
      DYCJ=(RI3Y*RCJ-CY1(J))*RTI 
      DZAI=(RI3Z*RAI-AZ1(I))*RTI 
      DZAJ=(RI3Z*RAJ-AZ1(J))*RTI 
      DZBI=RI3Z*RBI*RTI 
      DZBJ=RI3Z*RBJ*RTI 
      DZCI=(RI3Z*RCI-CZ1(I))*RTI 
      DZCJ=(RI3Z*RCJ-CZ1(J))*RTI 
      RTI3=3.0D0*RTI 
      DAAI=RTIM*(AX1(I)*AX1(I)+AY1(I)*AY1(I)+AZ1(I)*AZ1(I) 
     *         +DX*AAX1(I)+DY*AAY1(I)+DZ*AAZ1(I)) 
     *    +RTI3*RAI*RAI 
      DAAJ=RTIM*(AX1(J)*AX1(J)+AY1(J)*AY1(J)+AZ1(J)*AZ1(J) 
     *         -DX*AAX1(J)-DY*AAY1(J)-DZ*AAZ1(J)) 
     *    +RTI3*RAJ*RAJ 
      DABI=RTIM*(AX1(I)*BX1(I)+AY1(I)*BY1(I) 
     *         +DX*ABX1(I)+DY*ABY1(I)) 
     *    +RTI3*RAI*RBI 
      DABJ=RTIM*(AX1(J)*BX1(J)+AY1(J)*BY1(J) 
     *         -DX*ABX1(J)-DY*ABY1(J)) 
     *    +RTI3*RAJ*RBJ 
      DACI=RTIM*(AX1(I)*CX1(I)+AY1(I)*CY1(I)+AZ1(I)*CZ1(I) 
     *         +DX*ACX1(I)+DY*ACY1(I)+DZ*ACZ1(I)) 
     *    +RTI3*RAI*RCI 
      DACJ=RTIM*(AX1(J)*CX1(J)+AY1(J)*CY1(J)+AZ1(J)*CZ1(J) 
     *         -DX*ACX1(J)-DY*ACY1(J)-DZ*ACZ1(J)) 
     *    +RTI3*RAJ*RCJ 
      DBBI=RTIM*(BX1(I)*BX1(I)+BY1(I)*BY1(I) 
     *         +DX*BBX1(I)+DY*BBY1(I)) 
     *    +RTI3*RBI*RBI 
      DBBJ=RTIM*(BX1(J)*BX1(J)+BY1(J)*BY1(J) 
     *         -DX*BBX1(J)-DY*BBY1(J)) 
     *    +RTI3*RBJ*RBJ 
      DBCI=RTIM*(BX1(I)*CX1(I)+BY1(I)*CY1(I) 
     *         +DX*BCX1(I)+DY*BCY1(I)) 
     *    +RTI3*RBI*RCI 
      DBCJ=RTIM*(BX1(J)*CX1(J)+BY1(J)*CY1(J) 
     *         -DX*BCX1(J)-DY*BCY1(J)) 
     *    +RTI3*RBJ*RCJ 
      DCCI=RTIM*(CX1(I)*CX1(I)+CY1(I)*CY1(I)+CZ1(I)*CZ1(I) 
     *         +DX*CCX1(I)+DY*CCY1(I)+DZ*CCZ1(I)) 
     *    +RTI3*RCI*RCI 
      DCCJ=RTIM*(CX1(J)*CX1(J)+CY1(J)*CY1(J)+CZ1(J)*CZ1(J) 
     *         -DX*CCX1(J)-DY*CCY1(J)-DZ*CCZ1(J)) 
     *    +RTI3*RCJ*RCJ 
      DAAK=RTIM*(AX1(I)*AX1(J)+AY1(I)*AY1(J)+AZ1(I)*AZ1(J)) 
     *    +RTI3*RAI*RAJ 
      DABK=RTIM*(AX1(I)*BX1(J)+AY1(I)*BY1(J)) 
     *    +RTI3*RAI*RBJ 
      DBAK=RTIM*(BX1(I)*AX1(J)+BY1(I)*AY1(J)) 
     *    +RTI3*RAJ*RBI 
      DACK=RTIM*(AX1(I)*CX1(J)+AY1(I)*CY1(J)+AZ1(I)*CZ1(J)) 
     *    +RTI3*RAI*RCJ 
      DCAK=RTIM*(CX1(I)*AX1(J)+CY1(I)*AY1(J)+CZ1(I)*AZ1(J)) 
     *    +RTI3*RAJ*RCI 
      DBBK=RTIM*(BX1(I)*BX1(J)+BY1(I)*BY1(J)) 
     *    +RTI3*RBI*RBJ 
      DBCK=RTIM*(BX1(I)*CX1(J)+BY1(I)*CY1(J)) 
     *    +RTI3*RBI*RCJ 
      DCBK=RTIM*(CX1(I)*BX1(J)+CY1(I)*BY1(J)) 
     *    +RTI3*RBJ*RCI 
      DCCK=RTIM*(CX1(I)*CX1(J)+CY1(I)*CY1(J)+CZ1(I)*CZ1(J)) 
     *    +RTI3*RCI*RCJ 
      ENDIF
      IF (GTEST) THEN
      FDX=DDX 
      FDY=DDY 
      FDZ=DDZ 
      FDAI=DDAI 
      FDBI=DDBI 
      FDCI=DDCI 
      FDAJ=DDAJ 
      FDBJ=DDBJ 
      FDCJ=DDCJ 
      ENDIF
      IF (SSTEST) THEN
      SDXX=DXX 
      SDYY=DYY 
      SDZZ=DZZ 
      SDXY=DXY 
      SDYZ=DYZ 
      SDZX=DZX 
      SDXAI=DXAI 
      SDYAI=DYAI 
      SDZAI=DZAI 
      SDXBI=DXBI 
      SDYBI=DYBI 
      SDZBI=DZBI 
      SDXCI=DXCI 
      SDYCI=DYCI 
      SDZCI=DZCI 
      SDXAJ=DXAJ 
      SDYAJ=DYAJ 
      SDZAJ=DZAJ 
      SDXBJ=DXBJ 
      SDYBJ=DYBJ 
      SDZBJ=DZBJ 
      SDXCJ=DXCJ 
      SDYCJ=DYCJ 
      SDZCJ=DZCJ 
      SDAAI=DAAI 
      SDABI=DABI 
      SDACI=DACI 
      SDBBI=DBBI 
      SDBCI=DBCI 
      SDCCI=DCCI 
      SDAAJ=DAAJ 
      SDABJ=DABJ 
      SDACJ=DACJ 
      SDBBJ=DBBJ 
      SDBCJ=DBCJ 
      SDCCJ=DCCJ 
      SDAAK=DAAK 
      SDABK=DABK 
      SDBAK=DBAK 
      SDACK=DACK 
      SDCAK=DCAK 
      SDBBK=DBBK 
      SDBCK=DBCK 
      SDCBK=DCBK 
      SDCCK=DCCK 
      ENDIF
      EPI=RI 
C ::::: @@@ 11 :::::::::::::::::::::::::::::::::::::::::::::::::::: 
      DX=DX1-XK2 
      DY=DY1-YK2 
      DZ=DZ1-ZK2 
      SX=DX*DX 
      SY=DY*DY 
      SZ=DZ*DZ 
      RSI=1.0D0/(SX+SY+SZ) 
      RI=SQRT(RSI) 
      IF (GTEST.OR.SSTEST) THEN
      RSIM=-RSI 
      RTI=RI*RSI 
      RFI=RSI*RTI 
      XY=DX*DY 
      YZ=DY*DZ 
      ZX=DZ*DX 
      RTIM=-RTI 
      DDX=RTIM*DX 
      DDY=RTIM*DY 
      DDZ=RTIM*DZ 
      RAI=RI*(AX1(I)*DX+AY1(I)*DY+AZ1(I)*DZ) 
      DDAI=RSIM*RAI 
      RAJ=RI*(AX2(J)*DX+AY2(J)*DY+AZ2(J)*DZ) 
      DDAJ=RSIM*RAJ 
      RBI=RI*(BX1(I)*DX+BY1(I)*DY) 
      DDBI=RSIM*RBI 
      RBJ=RI*(BX2(J)*DX+BY2(J)*DY) 
      DDBJ=RSIM*RBJ 
      RCI=RI*(CX1(I)*DX+CY1(I)*DY+CZ1(I)*DZ) 
      DDCI=RSIM*RCI 
      RCJ=RI*(CX2(J)*DX+CY2(J)*DY+CZ2(J)*DZ) 
      DDCJ=RSIM*RCJ 
      ENDIF
      IF (SSTEST) THEN
      RFI3=RFI*3.0D0 
      DXX=RFI3*SX+RTIM 
      DYY=RFI3*SY+RTIM 
      DZZ=RFI3*SZ+RTIM 
      DXY=RFI3*XY 
      DYZ=RFI3*YZ 
      DZX=RFI3*ZX 
      RI3=3.0D0*RI 
      RI3X=DX*RI3 
      RI3Y=DY*RI3 
      RI3Z=DZ*RI3 
      DXAI=(RI3X*RAI-AX1(I))*RTI 
      DXAJ=(RI3X*RAJ-AX2(J))*RTI 
      DXBI=(RI3X*RBI-BX1(I))*RTI 
      DXBJ=(RI3X*RBJ-BX2(J))*RTI 
      DXCI=(RI3X*RCI-CX1(I))*RTI 
      DXCJ=(RI3X*RCJ-CX2(J))*RTI 
      DYAI=(RI3Y*RAI-AY1(I))*RTI 
      DYAJ=(RI3Y*RAJ-AY2(J))*RTI 
      DYBI=(RI3Y*RBI-BY1(I))*RTI 
      DYBJ=(RI3Y*RBJ-BY2(J))*RTI 
      DYCI=(RI3Y*RCI-CY1(I))*RTI 
      DYCJ=(RI3Y*RCJ-CY2(J))*RTI 
      DZAI=(RI3Z*RAI-AZ1(I))*RTI 
      DZAJ=(RI3Z*RAJ-AZ2(J))*RTI 
      DZBI=RI3Z*RBI*RTI 
      DZBJ=RI3Z*RBJ*RTI 
      DZCI=(RI3Z*RCI-CZ1(I))*RTI 
      DZCJ=(RI3Z*RCJ-CZ2(J))*RTI 
      RTI3=3.0D0*RTI 
      DAAI=RTIM*(AX1(I)*AX1(I)+AY1(I)*AY1(I)+AZ1(I)*AZ1(I) 
     *         +DX*AAX1(I)+DY*AAY1(I)+DZ*AAZ1(I)) 
     *    +RTI3*RAI*RAI 
      DAAJ=RTIM*(AX2(J)*AX2(J)+AY2(J)*AY2(J)+AZ2(J)*AZ2(J) 
     *         -DX*AAX2(J)-DY*AAY2(J)-DZ*AAZ2(J)) 
     *    +RTI3*RAJ*RAJ 
      DABI=RTIM*(AX1(I)*BX1(I)+AY1(I)*BY1(I) 
     *         +DX*ABX1(I)+DY*ABY1(I)) 
     *    +RTI3*RAI*RBI 
      DABJ=RTIM*(AX2(J)*BX2(J)+AY2(J)*BY2(J) 
     *         -DX*ABX2(J)-DY*ABY2(J)) 
     *    +RTI3*RAJ*RBJ 
      DACI=RTIM*(AX1(I)*CX1(I)+AY1(I)*CY1(I)+AZ1(I)*CZ1(I) 
     *         +DX*ACX1(I)+DY*ACY1(I)+DZ*ACZ1(I)) 
     *    +RTI3*RAI*RCI 
      DACJ=RTIM*(AX2(J)*CX2(J)+AY2(J)*CY2(J)+AZ2(J)*CZ2(J) 
     *         -DX*ACX2(J)-DY*ACY2(J)-DZ*ACZ2(J)) 
     *    +RTI3*RAJ*RCJ 
      DBBI=RTIM*(BX1(I)*BX1(I)+BY1(I)*BY1(I) 
     *         +DX*BBX1(I)+DY*BBY1(I)) 
     *    +RTI3*RBI*RBI 
      DBBJ=RTIM*(BX2(J)*BX2(J)+BY2(J)*BY2(J) 
     *         -DX*BBX2(J)-DY*BBY2(J)) 
     *    +RTI3*RBJ*RBJ 
      DBCI=RTIM*(BX1(I)*CX1(I)+BY1(I)*CY1(I) 
     *         +DX*BCX1(I)+DY*BCY1(I)) 
     *    +RTI3*RBI*RCI 
      DBCJ=RTIM*(BX2(J)*CX2(J)+BY2(J)*CY2(J) 
     *         -DX*BCX2(J)-DY*BCY2(J)) 
     *    +RTI3*RBJ*RCJ 
      DCCI=RTIM*(CX1(I)*CX1(I)+CY1(I)*CY1(I)+CZ1(I)*CZ1(I) 
     *         +DX*CCX1(I)+DY*CCY1(I)+DZ*CCZ1(I)) 
     *    +RTI3*RCI*RCI 
      DCCJ=RTIM*(CX2(J)*CX2(J)+CY2(J)*CY2(J)+CZ2(J)*CZ2(J) 
     *         -DX*CCX2(J)-DY*CCY2(J)-DZ*CCZ2(J)) 
     *    +RTI3*RCJ*RCJ 
      DAAK=RTIM*(AX1(I)*AX2(J)+AY1(I)*AY2(J)+AZ1(I)*AZ2(J)) 
     *    +RTI3*RAI*RAJ 
      DABK=RTIM*(AX1(I)*BX2(J)+AY1(I)*BY2(J)) 
     *    +RTI3*RAI*RBJ 
      DBAK=RTIM*(BX1(I)*AX2(J)+BY1(I)*AY2(J)) 
     *    +RTI3*RAJ*RBI 
      DACK=RTIM*(AX1(I)*CX2(J)+AY1(I)*CY2(J)+AZ1(I)*CZ2(J)) 
     *    +RTI3*RAI*RCJ 
      DCAK=RTIM*(CX1(I)*AX2(J)+CY1(I)*AY2(J)+CZ1(I)*AZ2(J)) 
     *    +RTI3*RAJ*RCI 
      DBBK=RTIM*(BX1(I)*BX2(J)+BY1(I)*BY2(J)) 
     *    +RTI3*RBI*RBJ 
      DBCK=RTIM*(BX1(I)*CX2(J)+BY1(I)*CY2(J)) 
     *    +RTI3*RBI*RCJ 
      DCBK=RTIM*(CX1(I)*BX2(J)+CY1(I)*BY2(J)) 
     *    +RTI3*RBJ*RCI 
      DCCK=RTIM*(CX1(I)*CX2(J)+CY1(I)*CY2(J)+CZ1(I)*CZ2(J)) 
     *    +RTI3*RCI*RCJ 
      ENDIF
      IF (GTEST) THEN
      FDX=DDX+FDX 
      FDY=DDY+FDY 
      FDZ=DDZ+FDZ 
      FDAI=FDAI+DDAI 
      FDBI=FDBI+DDBI 
      FDCI=FDCI+DDCI 
      FDAJ=FDAJ+DDAJ 
      FDBJ=FDBJ+DDBJ 
      FDCJ=FDCJ+DDCJ 
      ENDIF
      IF (SSTEST) THEN
      SDXX=DXX+SDXX 
      SDYY=DYY+SDYY 
      SDZZ=DZZ+SDZZ 
      SDXY=DXY+SDXY 
      SDYZ=DYZ+SDYZ 
      SDZX=DZX+SDZX 
      SDXAI=DXAI+SDXAI 
      SDYAI=DYAI+SDYAI 
      SDZAI=DZAI+SDZAI 
      SDXBI=DXBI+SDXBI 
      SDYBI=DYBI+SDYBI 
      SDZBI=DZBI+SDZBI 
      SDXCI=DXCI+SDXCI 
      SDYCI=DYCI+SDYCI 
      SDZCI=DZCI+SDZCI 
      SDXAJ=DXAJ+SDXAJ 
      SDYAJ=DYAJ+SDYAJ 
      SDZAJ=DZAJ+SDZAJ 
      SDXBJ=DXBJ+SDXBJ 
      SDYBJ=DYBJ+SDYBJ 
      SDZBJ=DZBJ+SDZBJ 
      SDXCJ=DXCJ+SDXCJ 
      SDYCJ=DYCJ+SDYCJ 
      SDZCJ=DZCJ+SDZCJ 
      SDAAI=DAAI+SDAAI 
      SDABI=DABI+SDABI 
      SDACI=DACI+SDACI 
      SDBBI=DBBI+SDBBI 
      SDBCI=DBCI+SDBCI 
      SDCCI=DCCI+SDCCI 
      SDAAJ=DAAJ+SDAAJ 
      SDABJ=DABJ+SDABJ 
      SDACJ=DACJ+SDACJ 
      SDBBJ=DBBJ+SDBBJ 
      SDBCJ=DBCJ+SDBCJ 
      SDCCJ=DCCJ+SDCCJ 
      SDAAK=DAAK+SDAAK 
      SDABK=DABK+SDABK 
      SDBAK=DBAK+SDBAK 
      SDACK=DACK+SDACK 
      SDCAK=DCAK+SDCAK 
      SDBBK=DBBK+SDBBK 
      SDBCK=DBCK+SDBCK 
      SDCBK=DCBK+SDCBK 
      SDCCK=DCCK+SDCCK 
      ENDIF
      EPI=EPI+RI 
C ::::: @@@ 12 :::::::::::::::::::::::::::::::::::::::::::::::::::: 
      DX=DX2-XK1 
      DY=DY2-YK1 
      DZ=DZ2-ZK1 
      SX=DX*DX 
      SY=DY*DY 
      SZ=DZ*DZ 
      RSI=1.0D0/(SX+SY+SZ) 
      RI=SQRT(RSI) 
      IF (GTEST.OR.SSTEST) THEN
      RSIM=-RSI 
      RTI=RI*RSI 
      RFI=RSI*RTI 
      XY=DX*DY 
      YZ=DY*DZ 
      ZX=DZ*DX 
      RTIM=-RTI 
      DDX=RTIM*DX 
      DDY=RTIM*DY 
      DDZ=RTIM*DZ 
      RAI=RI*(AX2(I)*DX+AY2(I)*DY+AZ2(I)*DZ) 
      DDAI=RSIM*RAI 
      RAJ=RI*(AX1(J)*DX+AY1(J)*DY+AZ1(J)*DZ) 
      DDAJ=RSIM*RAJ 
      RBI=RI*(BX2(I)*DX+BY2(I)*DY) 
      DDBI=RSIM*RBI 
      RBJ=RI*(BX1(J)*DX+BY1(J)*DY) 
      DDBJ=RSIM*RBJ 
      RCI=RI*(CX2(I)*DX+CY2(I)*DY+CZ2(I)*DZ) 
      DDCI=RSIM*RCI 
      RCJ=RI*(CX1(J)*DX+CY1(J)*DY+CZ1(J)*DZ) 
      DDCJ=RSIM*RCJ 
      ENDIF
      IF (SSTEST) THEN
      RFI3=RFI*3.0D0 
      DXX=RFI3*SX+RTIM 
      DYY=RFI3*SY+RTIM 
      DZZ=RFI3*SZ+RTIM 
      DXY=RFI3*XY 
      DYZ=RFI3*YZ 
      DZX=RFI3*ZX 
      RI3=3.0D0*RI 
      RI3X=DX*RI3 
      RI3Y=DY*RI3 
      RI3Z=DZ*RI3 
      DXAI=(RI3X*RAI-AX2(I))*RTI 
      DXAJ=(RI3X*RAJ-AX1(J))*RTI 
      DXBI=(RI3X*RBI-BX2(I))*RTI 
      DXBJ=(RI3X*RBJ-BX1(J))*RTI 
      DXCI=(RI3X*RCI-CX2(I))*RTI 
      DXCJ=(RI3X*RCJ-CX1(J))*RTI 
      DYAI=(RI3Y*RAI-AY2(I))*RTI 
      DYAJ=(RI3Y*RAJ-AY1(J))*RTI 
      DYBI=(RI3Y*RBI-BY2(I))*RTI 
      DYBJ=(RI3Y*RBJ-BY1(J))*RTI 
      DYCI=(RI3Y*RCI-CY2(I))*RTI 
      DYCJ=(RI3Y*RCJ-CY1(J))*RTI 
      DZAI=(RI3Z*RAI-AZ2(I))*RTI 
      DZAJ=(RI3Z*RAJ-AZ1(J))*RTI 
      DZBI=RI3Z*RBI*RTI 
      DZBJ=RI3Z*RBJ*RTI 
      DZCI=(RI3Z*RCI-CZ2(I))*RTI 
      DZCJ=(RI3Z*RCJ-CZ1(J))*RTI 
      RTI3=3.0D0*RTI 
      DAAI=RTIM*(AX2(I)*AX2(I)+AY2(I)*AY2(I)+AZ2(I)*AZ2(I) 
     *         +DX*AAX2(I)+DY*AAY2(I)+DZ*AAZ2(I)) 
     *    +RTI3*RAI*RAI 
      DAAJ=RTIM*(AX1(J)*AX1(J)+AY1(J)*AY1(J)+AZ1(J)*AZ1(J) 
     *         -DX*AAX1(J)-DY*AAY1(J)-DZ*AAZ1(J)) 
     *    +RTI3*RAJ*RAJ 
      DABI=RTIM*(AX2(I)*BX2(I)+AY2(I)*BY2(I) 
     *         +DX*ABX2(I)+DY*ABY2(I)) 
     *    +RTI3*RAI*RBI 
      DABJ=RTIM*(AX1(J)*BX1(J)+AY1(J)*BY1(J) 
     *         -DX*ABX1(J)-DY*ABY1(J)) 
     *    +RTI3*RAJ*RBJ 
      DACI=RTIM*(AX2(I)*CX2(I)+AY2(I)*CY2(I)+AZ2(I)*CZ2(I) 
     *         +DX*ACX2(I)+DY*ACY2(I)+DZ*ACZ2(I)) 
     *    +RTI3*RAI*RCI 
      DACJ=RTIM*(AX1(J)*CX1(J)+AY1(J)*CY1(J)+AZ1(J)*CZ1(J) 
     *         -DX*ACX1(J)-DY*ACY1(J)-DZ*ACZ1(J)) 
     *    +RTI3*RAJ*RCJ 
      DBBI=RTIM*(BX2(I)*BX2(I)+BY2(I)*BY2(I) 
     *         +DX*BBX2(I)+DY*BBY2(I)) 
     *    +RTI3*RBI*RBI 
      DBBJ=RTIM*(BX1(J)*BX1(J)+BY1(J)*BY1(J) 
     *         -DX*BBX1(J)-DY*BBY1(J)) 
     *    +RTI3*RBJ*RBJ 
      DBCI=RTIM*(BX2(I)*CX2(I)+BY2(I)*CY2(I) 
     *         +DX*BCX2(I)+DY*BCY2(I)) 
     *    +RTI3*RBI*RCI 
      DBCJ=RTIM*(BX1(J)*CX1(J)+BY1(J)*CY1(J) 
     *         -DX*BCX1(J)-DY*BCY1(J)) 
     *    +RTI3*RBJ*RCJ 
      DCCI=RTIM*(CX2(I)*CX2(I)+CY2(I)*CY2(I)+CZ2(I)*CZ2(I) 
     *         +DX*CCX2(I)+DY*CCY2(I)+DZ*CCZ2(I)) 
     *    +RTI3*RCI*RCI 
      DCCJ=RTIM*(CX1(J)*CX1(J)+CY1(J)*CY1(J)+CZ1(J)*CZ1(J) 
     *         -DX*CCX1(J)-DY*CCY1(J)-DZ*CCZ1(J)) 
     *    +RTI3*RCJ*RCJ 
      DAAK=RTIM*(AX2(I)*AX1(J)+AY2(I)*AY1(J)+AZ2(I)*AZ1(J)) 
     *    +RTI3*RAI*RAJ 
      DABK=RTIM*(AX2(I)*BX1(J)+AY2(I)*BY1(J)) 
     *    +RTI3*RAI*RBJ 
      DBAK=RTIM*(BX2(I)*AX1(J)+BY2(I)*AY1(J)) 
     *    +RTI3*RAJ*RBI 
      DACK=RTIM*(AX2(I)*CX1(J)+AY2(I)*CY1(J)+AZ2(I)*CZ1(J)) 
     *    +RTI3*RAI*RCJ 
      DCAK=RTIM*(CX2(I)*AX1(J)+CY2(I)*AY1(J)+CZ2(I)*AZ1(J)) 
     *    +RTI3*RAJ*RCI 
      DBBK=RTIM*(BX2(I)*BX1(J)+BY2(I)*BY1(J)) 
     *    +RTI3*RBI*RBJ 
      DBCK=RTIM*(BX2(I)*CX1(J)+BY2(I)*CY1(J)) 
     *    +RTI3*RBI*RCJ 
      DCBK=RTIM*(CX2(I)*BX1(J)+CY2(I)*BY1(J)) 
     *    +RTI3*RBJ*RCI 
      DCCK=RTIM*(CX2(I)*CX1(J)+CY2(I)*CY1(J)+CZ2(I)*CZ1(J)) 
     *    +RTI3*RCI*RCJ 
      ENDIF
      IF (GTEST) THEN
      FDX=DDX+FDX 
      FDY=DDY+FDY 
      FDZ=DDZ+FDZ 
      FDAI=FDAI+DDAI 
      FDBI=FDBI+DDBI 
      FDCI=FDCI+DDCI 
      FDAJ=FDAJ+DDAJ 
      FDBJ=FDBJ+DDBJ 
      FDCJ=FDCJ+DDCJ 
      ENDIF
      IF (SSTEST) THEN
      SDXX=DXX+SDXX 
      SDYY=DYY+SDYY 
      SDZZ=DZZ+SDZZ 
      SDXY=DXY+SDXY 
      SDYZ=DYZ+SDYZ 
      SDZX=DZX+SDZX 
      SDXAI=DXAI+SDXAI 
      SDYAI=DYAI+SDYAI 
      SDZAI=DZAI+SDZAI 
      SDXBI=DXBI+SDXBI 
      SDYBI=DYBI+SDYBI 
      SDZBI=DZBI+SDZBI 
      SDXCI=DXCI+SDXCI 
      SDYCI=DYCI+SDYCI 
      SDZCI=DZCI+SDZCI 
      SDXAJ=DXAJ+SDXAJ 
      SDYAJ=DYAJ+SDYAJ 
      SDZAJ=DZAJ+SDZAJ 
      SDXBJ=DXBJ+SDXBJ 
      SDYBJ=DYBJ+SDYBJ 
      SDZBJ=DZBJ+SDZBJ 
      SDXCJ=DXCJ+SDXCJ 
      SDYCJ=DYCJ+SDYCJ 
      SDZCJ=DZCJ+SDZCJ 
      SDAAI=DAAI+SDAAI 
      SDABI=DABI+SDABI 
      SDACI=DACI+SDACI 
      SDBBI=DBBI+SDBBI 
      SDBCI=DBCI+SDBCI 
      SDCCI=DCCI+SDCCI 
      SDAAJ=DAAJ+SDAAJ 
      SDABJ=DABJ+SDABJ 
      SDACJ=DACJ+SDACJ 
      SDBBJ=DBBJ+SDBBJ 
      SDBCJ=DBCJ+SDBCJ 
      SDCCJ=DCCJ+SDCCJ 
      SDAAK=DAAK+SDAAK 
      SDABK=DABK+SDABK 
      SDBAK=DBAK+SDBAK 
      SDACK=DACK+SDACK 
      SDCAK=DCAK+SDCAK 
      SDBBK=DBBK+SDBBK 
      SDBCK=DBCK+SDBCK 
      SDCBK=DCBK+SDCBK 
      SDCCK=DCCK+SDCCK 
      ENDIF
      EPI=EPI+RI 
C ::::: @@@ 21 :::::::::::::::::::::::::::::::::::::::::::::::::::: 
      DX=DX2-XK2 
      DY=DY2-YK2 
      DZ=DZ2-ZK2 
      SX=DX*DX 
      SY=DY*DY 
      SZ=DZ*DZ 
      RSI=1.0D0/(SX+SY+SZ) 
      RI=SQRT(RSI) 
      IF (GTEST.OR.SSTEST) THEN
      RSIM=-RSI 
      RTI=RI*RSI 
      RFI=RSI*RTI 
      XY=DX*DY 
      YZ=DY*DZ 
      ZX=DZ*DX 
      RTIM=-RTI 
      DDX=RTIM*DX 
      DDY=RTIM*DY 
      DDZ=RTIM*DZ 
      RAI=RI*(AX2(I)*DX+AY2(I)*DY+AZ2(I)*DZ) 
      DDAI=RSIM*RAI 
      RAJ=RI*(AX2(J)*DX+AY2(J)*DY+AZ2(J)*DZ) 
      DDAJ=RSIM*RAJ 
      RBI=RI*(BX2(I)*DX+BY2(I)*DY) 
      DDBI=RSIM*RBI 
      RBJ=RI*(BX2(J)*DX+BY2(J)*DY) 
      DDBJ=RSIM*RBJ 
      RCI=RI*(CX2(I)*DX+CY2(I)*DY+CZ2(I)*DZ) 
      DDCI=RSIM*RCI 
      RCJ=RI*(CX2(J)*DX+CY2(J)*DY+CZ2(J)*DZ) 
      DDCJ=RSIM*RCJ 
      ENDIF
      IF (SSTEST) THEN
      RFI3=RFI*3.0D0 
      DXX=RFI3*SX+RTIM 
      DYY=RFI3*SY+RTIM 
      DZZ=RFI3*SZ+RTIM 
      DXY=RFI3*XY 
      DYZ=RFI3*YZ 
      DZX=RFI3*ZX 
      RI3=3.0D0*RI 
      RI3X=DX*RI3 
      RI3Y=DY*RI3 
      RI3Z=DZ*RI3 
      DXAI=(RI3X*RAI-AX2(I))*RTI 
      DXAJ=(RI3X*RAJ-AX2(J))*RTI 
      DXBI=(RI3X*RBI-BX2(I))*RTI 
      DXBJ=(RI3X*RBJ-BX2(J))*RTI 
      DXCI=(RI3X*RCI-CX2(I))*RTI 
      DXCJ=(RI3X*RCJ-CX2(J))*RTI 
      DYAI=(RI3Y*RAI-AY2(I))*RTI 
      DYAJ=(RI3Y*RAJ-AY2(J))*RTI 
      DYBI=(RI3Y*RBI-BY2(I))*RTI 
      DYBJ=(RI3Y*RBJ-BY2(J))*RTI 
      DYCI=(RI3Y*RCI-CY2(I))*RTI 
      DYCJ=(RI3Y*RCJ-CY2(J))*RTI 
      DZAI=(RI3Z*RAI-AZ2(I))*RTI 
      DZAJ=(RI3Z*RAJ-AZ2(J))*RTI 
      DZBI=RI3Z*RBI*RTI 
      DZBJ=RI3Z*RBJ*RTI 
      DZCI=(RI3Z*RCI-CZ2(I))*RTI 
      DZCJ=(RI3Z*RCJ-CZ2(J))*RTI 
      RTI3=3.0D0*RTI 
      DAAI=RTIM*(AX2(I)*AX2(I)+AY2(I)*AY2(I)+AZ2(I)*AZ2(I) 
     *         +DX*AAX2(I)+DY*AAY2(I)+DZ*AAZ2(I)) 
     *    +RTI3*RAI*RAI 
      DAAJ=RTIM*(AX2(J)*AX2(J)+AY2(J)*AY2(J)+AZ2(J)*AZ2(J) 
     *         -DX*AAX2(J)-DY*AAY2(J)-DZ*AAZ2(J)) 
     *    +RTI3*RAJ*RAJ 
      DABI=RTIM*(AX2(I)*BX2(I)+AY2(I)*BY2(I) 
     *         +DX*ABX2(I)+DY*ABY2(I)) 
     *    +RTI3*RAI*RBI 
      DABJ=RTIM*(AX2(J)*BX2(J)+AY2(J)*BY2(J) 
     *         -DX*ABX2(J)-DY*ABY2(J)) 
     *    +RTI3*RAJ*RBJ 
      DACI=RTIM*(AX2(I)*CX2(I)+AY2(I)*CY2(I)+AZ2(I)*CZ2(I) 
     *         +DX*ACX2(I)+DY*ACY2(I)+DZ*ACZ2(I)) 
     *    +RTI3*RAI*RCI 
      DACJ=RTIM*(AX2(J)*CX2(J)+AY2(J)*CY2(J)+AZ2(J)*CZ2(J) 
     *         -DX*ACX2(J)-DY*ACY2(J)-DZ*ACZ2(J)) 
     *    +RTI3*RAJ*RCJ 
      DBBI=RTIM*(BX2(I)*BX2(I)+BY2(I)*BY2(I) 
     *         +DX*BBX2(I)+DY*BBY2(I)) 
     *    +RTI3*RBI*RBI 
      DBBJ=RTIM*(BX2(J)*BX2(J)+BY2(J)*BY2(J) 
     *         -DX*BBX2(J)-DY*BBY2(J)) 
     *    +RTI3*RBJ*RBJ 
      DBCI=RTIM*(BX2(I)*CX2(I)+BY2(I)*CY2(I) 
     *         +DX*BCX2(I)+DY*BCY2(I)) 
     *    +RTI3*RBI*RCI 
      DBCJ=RTIM*(BX2(J)*CX2(J)+BY2(J)*CY2(J) 
     *         -DX*BCX2(J)-DY*BCY2(J)) 
     *    +RTI3*RBJ*RCJ 
      DCCI=RTIM*(CX2(I)*CX2(I)+CY2(I)*CY2(I)+CZ2(I)*CZ2(I) 
     *         +DX*CCX2(I)+DY*CCY2(I)+DZ*CCZ2(I)) 
     *    +RTI3*RCI*RCI 
      DCCJ=RTIM*(CX2(J)*CX2(J)+CY2(J)*CY2(J)+CZ2(J)*CZ2(J) 
     *         -DX*CCX2(J)-DY*CCY2(J)-DZ*CCZ2(J)) 
     *    +RTI3*RCJ*RCJ 
      DAAK=RTIM*(AX2(I)*AX2(J)+AY2(I)*AY2(J)+AZ2(I)*AZ2(J)) 
     *    +RTI3*RAI*RAJ 
      DABK=RTIM*(AX2(I)*BX2(J)+AY2(I)*BY2(J)) 
     *    +RTI3*RAI*RBJ 
      DBAK=RTIM*(BX2(I)*AX2(J)+BY2(I)*AY2(J)) 
     *    +RTI3*RAJ*RBI 
      DACK=RTIM*(AX2(I)*CX2(J)+AY2(I)*CY2(J)+AZ2(I)*CZ2(J)) 
     *    +RTI3*RAI*RCJ 
      DCAK=RTIM*(CX2(I)*AX2(J)+CY2(I)*AY2(J)+CZ2(I)*AZ2(J)) 
     *    +RTI3*RAJ*RCI 
      DBBK=RTIM*(BX2(I)*BX2(J)+BY2(I)*BY2(J)) 
     *    +RTI3*RBI*RBJ 
      DBCK=RTIM*(BX2(I)*CX2(J)+BY2(I)*CY2(J)) 
     *    +RTI3*RBI*RCJ 
      DCBK=RTIM*(CX2(I)*BX2(J)+CY2(I)*BY2(J)) 
     *    +RTI3*RBJ*RCI 
      DCCK=RTIM*(CX2(I)*CX2(J)+CY2(I)*CY2(J)+CZ2(I)*CZ2(J)) 
     *    +RTI3*RCI*RCJ 
      ENDIF
      IF (GTEST) THEN
      FDX=DDX+FDX 
      FDY=DDY+FDY 
      FDZ=DDZ+FDZ 
      FDAI=FDAI+DDAI 
      FDBI=FDBI+DDBI 
      FDCI=FDCI+DDCI 
      FDAJ=FDAJ+DDAJ 
      FDBJ=FDBJ+DDBJ 
      FDCJ=FDCJ+DDCJ 
      ENDIF
      IF (SSTEST) THEN
      SDXX=DXX+SDXX 
      SDYY=DYY+SDYY 
      SDZZ=DZZ+SDZZ 
      SDXY=DXY+SDXY 
      SDYZ=DYZ+SDYZ 
      SDZX=DZX+SDZX 
      SDXAI=DXAI+SDXAI 
      SDYAI=DYAI+SDYAI 
      SDZAI=DZAI+SDZAI 
      SDXBI=DXBI+SDXBI 
      SDYBI=DYBI+SDYBI 
      SDZBI=DZBI+SDZBI 
      SDXCI=DXCI+SDXCI 
      SDYCI=DYCI+SDYCI 
      SDZCI=DZCI+SDZCI 
      SDXAJ=DXAJ+SDXAJ 
      SDYAJ=DYAJ+SDYAJ 
      SDZAJ=DZAJ+SDZAJ 
      SDXBJ=DXBJ+SDXBJ 
      SDYBJ=DYBJ+SDYBJ 
      SDZBJ=DZBJ+SDZBJ 
      SDXCJ=DXCJ+SDXCJ 
      SDYCJ=DYCJ+SDYCJ 
      SDZCJ=DZCJ+SDZCJ 
      SDAAI=DAAI+SDAAI 
      SDABI=DABI+SDABI 
      SDACI=DACI+SDACI 
      SDBBI=DBBI+SDBBI 
      SDBCI=DBCI+SDBCI 
      SDCCI=DCCI+SDCCI 
      SDAAJ=DAAJ+SDAAJ 
      SDABJ=DABJ+SDABJ 
      SDACJ=DACJ+SDACJ 
      SDBBJ=DBBJ+SDBBJ 
      SDBCJ=DBCJ+SDBCJ 
      SDCCJ=DCCJ+SDCCJ 
      SDAAK=DAAK+SDAAK 
      SDABK=DABK+SDABK 
      SDBAK=DBAK+SDBAK 
      SDACK=DACK+SDACK 
      SDCAK=DCAK+SDCAK 
      SDBBK=DBBK+SDBBK 
      SDBCK=DBCK+SDBCK 
      SDCBK=DCBK+SDCBK 
      SDCCK=DCCK+SDCCK 
      ENDIF
      EPI=EPI+RI 
C ::::: @@@ 22 :::::::::::::::::::::::::::::::::::::::::::::::::::: 
      DX=DX1-XK4 
      DY=DY1-YK4 
      DZ=DZ1-ZK4 
      SX=DX*DX 
      SY=DY*DY 
      SZ=DZ*DZ 
      RSI=1.0D0/(SX+SY+SZ) 
      RI=SQRT(RSI) 
      IF (GTEST.OR.SSTEST) THEN
      RSIM=2.0D0*RSI 
      RTI=-RI*RSIM 
      RFI=RSI*RTI 
      XY=DX*DY 
      YZ=DY*DZ 
      ZX=DZ*DX 
      RTIM=-RTI 
      DDX=RTIM*DX 
      DDY=RTIM*DY 
      DDZ=RTIM*DZ 
      RAI=RI*(AX1(I)*DX+AY1(I)*DY+AZ1(I)*DZ) 
      DDAI=RSIM*RAI 
      RAJ=RI*(AX4(J)*DX+AY4(J)*DY+AZ4(J)*DZ) 
      DDAJ=RSIM*RAJ 
      RBI=RI*(BX1(I)*DX+BY1(I)*DY) 
      DDBI=RSIM*RBI 
      RBJ=RI*(BX4(J)*DX+BY4(J)*DY) 
      DDBJ=RSIM*RBJ 
      RCI=RI*(CX1(I)*DX+CY1(I)*DY+CZ1(I)*DZ) 
      DDCI=RSIM*RCI 
      RCJ=RI*(CX4(J)*DX+CY4(J)*DY+CZ4(J)*DZ) 
      DDCJ=RSIM*RCJ 
      ENDIF
      IF (SSTEST) THEN
      RFI3=RFI*3.0D0 
      DXX=RFI3*SX+RTIM 
      DYY=RFI3*SY+RTIM 
      DZZ=RFI3*SZ+RTIM 
      DXY=RFI3*XY 
      DYZ=RFI3*YZ 
      DZX=RFI3*ZX 
      RI3=3.0D0*RI 
      RI3X=DX*RI3 
      RI3Y=DY*RI3 
      RI3Z=DZ*RI3 
      DXAI=(RI3X*RAI-AX1(I))*RTI 
      DXAJ=(RI3X*RAJ-AX4(J))*RTI 
      DXBI=(RI3X*RBI-BX1(I))*RTI 
      DXBJ=(RI3X*RBJ-BX4(J))*RTI 
      DXCI=(RI3X*RCI-CX1(I))*RTI 
      DXCJ=(RI3X*RCJ-CX4(J))*RTI 
      DYAI=(RI3Y*RAI-AY1(I))*RTI 
      DYAJ=(RI3Y*RAJ-AY4(J))*RTI 
      DYBI=(RI3Y*RBI-BY1(I))*RTI 
      DYBJ=(RI3Y*RBJ-BY4(J))*RTI 
      DYCI=(RI3Y*RCI-CY1(I))*RTI 
      DYCJ=(RI3Y*RCJ-CY4(J))*RTI 
      DZAI=(RI3Z*RAI-AZ1(I))*RTI 
      DZAJ=(RI3Z*RAJ-AZ4(J))*RTI 
      DZBI=RI3Z*RBI*RTI 
      DZBJ=RI3Z*RBJ*RTI 
      DZCI=(RI3Z*RCI-CZ1(I))*RTI 
      DZCJ=(RI3Z*RCJ-CZ4(J))*RTI 
      RTI3=3.0D0*RTI 
      DAAI=RTIM*(AX1(I)*AX1(I)+AY1(I)*AY1(I)+AZ1(I)*AZ1(I) 
     *         +DX*AAX1(I)+DY*AAY1(I)+DZ*AAZ1(I)) 
     *    +RTI3*RAI*RAI 
      DAAJ=RTIM*(AX4(J)*AX4(J)+AY4(J)*AY4(J)+AZ4(J)*AZ4(J) 
     *         -DX*AAX4(J)-DY*AAY4(J)-DZ*AAZ4(J)) 
     *    +RTI3*RAJ*RAJ 
      DABI=RTIM*(AX1(I)*BX1(I)+AY1(I)*BY1(I) 
     *         +DX*ABX1(I)+DY*ABY1(I)) 
     *    +RTI3*RAI*RBI 
      DABJ=RTIM*(AX4(J)*BX4(J)+AY4(J)*BY4(J) 
     *         -DX*ABX4(J)-DY*ABY4(J)) 
     *    +RTI3*RAJ*RBJ 
      DACI=RTIM*(AX1(I)*CX1(I)+AY1(I)*CY1(I)+AZ1(I)*CZ1(I) 
     *         +DX*ACX1(I)+DY*ACY1(I)+DZ*ACZ1(I)) 
     *    +RTI3*RAI*RCI 
      DACJ=RTIM*(AX4(J)*CX4(J)+AY4(J)*CY4(J)+AZ4(J)*CZ4(J) 
     *         -DX*ACX4(J)-DY*ACY4(J)-DZ*ACZ4(J)) 
     *    +RTI3*RAJ*RCJ 
      DBBI=RTIM*(BX1(I)*BX1(I)+BY1(I)*BY1(I) 
     *         +DX*BBX1(I)+DY*BBY1(I)) 
     *    +RTI3*RBI*RBI 
      DBBJ=RTIM*(BX4(J)*BX4(J)+BY4(J)*BY4(J) 
     *         -DX*BBX4(J)-DY*BBY4(J)) 
     *    +RTI3*RBJ*RBJ 
      DBCI=RTIM*(BX1(I)*CX1(I)+BY1(I)*CY1(I) 
     *         +DX*BCX1(I)+DY*BCY1(I)) 
     *    +RTI3*RBI*RCI 
      DBCJ=RTIM*(BX4(J)*CX4(J)+BY4(J)*CY4(J) 
     *         -DX*BCX4(J)-DY*BCY4(J)) 
     *    +RTI3*RBJ*RCJ 
      DCCI=RTIM*(CX1(I)*CX1(I)+CY1(I)*CY1(I)+CZ1(I)*CZ1(I) 
     *         +DX*CCX1(I)+DY*CCY1(I)+DZ*CCZ1(I)) 
     *    +RTI3*RCI*RCI 
      DCCJ=RTIM*(CX4(J)*CX4(J)+CY4(J)*CY4(J)+CZ4(J)*CZ4(J) 
     *         -DX*CCX4(J)-DY*CCY4(J)-DZ*CCZ4(J)) 
     *    +RTI3*RCJ*RCJ 
      DAAK=RTIM*(AX1(I)*AX4(J)+AY1(I)*AY4(J)+AZ1(I)*AZ4(J)) 
     *    +RTI3*RAI*RAJ 
      DABK=RTIM*(AX1(I)*BX4(J)+AY1(I)*BY4(J)) 
     *    +RTI3*RAI*RBJ 
      DBAK=RTIM*(BX1(I)*AX4(J)+BY1(I)*AY4(J)) 
     *    +RTI3*RAJ*RBI 
      DACK=RTIM*(AX1(I)*CX4(J)+AY1(I)*CY4(J)+AZ1(I)*CZ4(J)) 
     *    +RTI3*RAI*RCJ 
      DCAK=RTIM*(CX1(I)*AX4(J)+CY1(I)*AY4(J)+CZ1(I)*AZ4(J)) 
     *    +RTI3*RAJ*RCI 
      DBBK=RTIM*(BX1(I)*BX4(J)+BY1(I)*BY4(J)) 
     *    +RTI3*RBI*RBJ 
      DBCK=RTIM*(BX1(I)*CX4(J)+BY1(I)*CY4(J)) 
     *    +RTI3*RBI*RCJ 
      DCBK=RTIM*(CX1(I)*BX4(J)+CY1(I)*BY4(J)) 
     *    +RTI3*RBJ*RCI 
      DCCK=RTIM*(CX1(I)*CX4(J)+CY1(I)*CY4(J)+CZ1(I)*CZ4(J)) 
     *    +RTI3*RCI*RCJ 
      ENDIF
      IF (GTEST) THEN
      FDX=DDX+FDX 
      FDY=DDY+FDY 
      FDZ=DDZ+FDZ 
      FDAI=FDAI+DDAI 
      FDBI=FDBI+DDBI 
      FDCI=FDCI+DDCI 
      FDAJ=FDAJ+DDAJ 
      FDBJ=FDBJ+DDBJ 
      FDCJ=FDCJ+DDCJ 
      ENDIF
      IF (SSTEST) THEN
      SDXX=DXX+SDXX 
      SDYY=DYY+SDYY 
      SDZZ=DZZ+SDZZ 
      SDXY=DXY+SDXY 
      SDYZ=DYZ+SDYZ 
      SDZX=DZX+SDZX 
      SDXAI=DXAI+SDXAI 
      SDYAI=DYAI+SDYAI 
      SDZAI=DZAI+SDZAI 
      SDXBI=DXBI+SDXBI 
      SDYBI=DYBI+SDYBI 
      SDZBI=DZBI+SDZBI 
      SDXCI=DXCI+SDXCI 
      SDYCI=DYCI+SDYCI 
      SDZCI=DZCI+SDZCI 
      SDXAJ=DXAJ+SDXAJ 
      SDYAJ=DYAJ+SDYAJ 
      SDZAJ=DZAJ+SDZAJ 
      SDXBJ=DXBJ+SDXBJ 
      SDYBJ=DYBJ+SDYBJ 
      SDZBJ=DZBJ+SDZBJ 
      SDXCJ=DXCJ+SDXCJ 
      SDYCJ=DYCJ+SDYCJ 
      SDZCJ=DZCJ+SDZCJ 
      SDAAI=DAAI+SDAAI 
      SDABI=DABI+SDABI 
      SDACI=DACI+SDACI 
      SDBBI=DBBI+SDBBI 
      SDBCI=DBCI+SDBCI 
      SDCCI=DCCI+SDCCI 
      SDAAJ=DAAJ+SDAAJ 
      SDABJ=DABJ+SDABJ 
      SDACJ=DACJ+SDACJ 
      SDBBJ=DBBJ+SDBBJ 
      SDBCJ=DBCJ+SDBCJ 
      SDCCJ=DCCJ+SDCCJ 
      SDAAK=DAAK+SDAAK 
      SDABK=DABK+SDABK 
      SDBAK=DBAK+SDBAK 
      SDACK=DACK+SDACK 
      SDCAK=DCAK+SDCAK 
      SDBBK=DBBK+SDBBK 
      SDBCK=DBCK+SDBCK 
      SDCBK=DCBK+SDCBK 
      SDCCK=DCCK+SDCCK 
      ENDIF
      EPI=EPI-2.0D0*RI 
C ::::: @@@ 14 :::::::::::::::::::::::::::::::::::::::::::::::::::: 
      DX=DX2-XK4 
      DY=DY2-YK4 
      DZ=DZ2-ZK4 
      SX=DX*DX 
      SY=DY*DY 
      SZ=DZ*DZ 
      RSI=1.0D0/(SX+SY+SZ) 
      RI=SQRT(RSI) 
      IF (GTEST.OR.SSTEST) THEN
      RSIM=2.0D0*RSI 
      RTI=-RI*RSIM 
      RFI=RSI*RTI 
      XY=DX*DY 
      YZ=DY*DZ 
      ZX=DZ*DX 
      RTIM=-RTI 
      DDX=RTIM*DX 
      DDY=RTIM*DY 
      DDZ=RTIM*DZ 
      RAI=RI*(AX2(I)*DX+AY2(I)*DY+AZ2(I)*DZ) 
      DDAI=RSIM*RAI 
      RAJ=RI*(AX4(J)*DX+AY4(J)*DY+AZ4(J)*DZ) 
      DDAJ=RSIM*RAJ 
      RBI=RI*(BX2(I)*DX+BY2(I)*DY) 
      DDBI=RSIM*RBI 
      RBJ=RI*(BX4(J)*DX+BY4(J)*DY) 
      DDBJ=RSIM*RBJ 
      RCI=RI*(CX2(I)*DX+CY2(I)*DY+CZ2(I)*DZ) 
      DDCI=RSIM*RCI 
      RCJ=RI*(CX4(J)*DX+CY4(J)*DY+CZ4(J)*DZ) 
      DDCJ=RSIM*RCJ 
      ENDIF
      IF (SSTEST) THEN
      RFI3=RFI*3.0D0 
      DXX=RFI3*SX+RTIM 
      DYY=RFI3*SY+RTIM 
      DZZ=RFI3*SZ+RTIM 
      DXY=RFI3*XY 
      DYZ=RFI3*YZ 
      DZX=RFI3*ZX 
      RI3=3.0D0*RI 
      RI3X=DX*RI3 
      RI3Y=DY*RI3 
      RI3Z=DZ*RI3 
      DXAI=(RI3X*RAI-AX2(I))*RTI 
      DXAJ=(RI3X*RAJ-AX4(J))*RTI 
      DXBI=(RI3X*RBI-BX2(I))*RTI 
      DXBJ=(RI3X*RBJ-BX4(J))*RTI 
      DXCI=(RI3X*RCI-CX2(I))*RTI 
      DXCJ=(RI3X*RCJ-CX4(J))*RTI 
      DYAI=(RI3Y*RAI-AY2(I))*RTI 
      DYAJ=(RI3Y*RAJ-AY4(J))*RTI 
      DYBI=(RI3Y*RBI-BY2(I))*RTI 
      DYBJ=(RI3Y*RBJ-BY4(J))*RTI 
      DYCI=(RI3Y*RCI-CY2(I))*RTI 
      DYCJ=(RI3Y*RCJ-CY4(J))*RTI 
      DZAI=(RI3Z*RAI-AZ2(I))*RTI 
      DZAJ=(RI3Z*RAJ-AZ4(J))*RTI 
      DZBI=RI3Z*RBI*RTI 
      DZBJ=RI3Z*RBJ*RTI 
      DZCI=(RI3Z*RCI-CZ2(I))*RTI 
      DZCJ=(RI3Z*RCJ-CZ4(J))*RTI 
      RTI3=3.0D0*RTI 
      DAAI=RTIM*(AX2(I)*AX2(I)+AY2(I)*AY2(I)+AZ2(I)*AZ2(I) 
     *         +DX*AAX2(I)+DY*AAY2(I)+DZ*AAZ2(I)) 
     *    +RTI3*RAI*RAI 
      DAAJ=RTIM*(AX4(J)*AX4(J)+AY4(J)*AY4(J)+AZ4(J)*AZ4(J) 
     *         -DX*AAX4(J)-DY*AAY4(J)-DZ*AAZ4(J)) 
     *    +RTI3*RAJ*RAJ 
      DABI=RTIM*(AX2(I)*BX2(I)+AY2(I)*BY2(I) 
     *         +DX*ABX2(I)+DY*ABY2(I)) 
     *    +RTI3*RAI*RBI 
      DABJ=RTIM*(AX4(J)*BX4(J)+AY4(J)*BY4(J) 
     *         -DX*ABX4(J)-DY*ABY4(J)) 
     *    +RTI3*RAJ*RBJ 
      DACI=RTIM*(AX2(I)*CX2(I)+AY2(I)*CY2(I)+AZ2(I)*CZ2(I) 
     *         +DX*ACX2(I)+DY*ACY2(I)+DZ*ACZ2(I)) 
     *    +RTI3*RAI*RCI 
      DACJ=RTIM*(AX4(J)*CX4(J)+AY4(J)*CY4(J)+AZ4(J)*CZ4(J) 
     *         -DX*ACX4(J)-DY*ACY4(J)-DZ*ACZ4(J)) 
     *    +RTI3*RAJ*RCJ 
      DBBI=RTIM*(BX2(I)*BX2(I)+BY2(I)*BY2(I) 
     *         +DX*BBX2(I)+DY*BBY2(I)) 
     *    +RTI3*RBI*RBI 
      DBBJ=RTIM*(BX4(J)*BX4(J)+BY4(J)*BY4(J) 
     *         -DX*BBX4(J)-DY*BBY4(J)) 
     *    +RTI3*RBJ*RBJ 
      DBCI=RTIM*(BX2(I)*CX2(I)+BY2(I)*CY2(I) 
     *         +DX*BCX2(I)+DY*BCY2(I)) 
     *    +RTI3*RBI*RCI 
      DBCJ=RTIM*(BX4(J)*CX4(J)+BY4(J)*CY4(J) 
     *         -DX*BCX4(J)-DY*BCY4(J)) 
     *    +RTI3*RBJ*RCJ 
      DCCI=RTIM*(CX2(I)*CX2(I)+CY2(I)*CY2(I)+CZ2(I)*CZ2(I) 
     *         +DX*CCX2(I)+DY*CCY2(I)+DZ*CCZ2(I)) 
     *    +RTI3*RCI*RCI 
      DCCJ=RTIM*(CX4(J)*CX4(J)+CY4(J)*CY4(J)+CZ4(J)*CZ4(J) 
     *         -DX*CCX4(J)-DY*CCY4(J)-DZ*CCZ4(J)) 
     *    +RTI3*RCJ*RCJ 
      DAAK=RTIM*(AX2(I)*AX4(J)+AY2(I)*AY4(J)+AZ2(I)*AZ4(J)) 
     *    +RTI3*RAI*RAJ 
      DABK=RTIM*(AX2(I)*BX4(J)+AY2(I)*BY4(J)) 
     *    +RTI3*RAI*RBJ 
      DBAK=RTIM*(BX2(I)*AX4(J)+BY2(I)*AY4(J)) 
     *    +RTI3*RAJ*RBI 
      DACK=RTIM*(AX2(I)*CX4(J)+AY2(I)*CY4(J)+AZ2(I)*CZ4(J)) 
     *    +RTI3*RAI*RCJ 
      DCAK=RTIM*(CX2(I)*AX4(J)+CY2(I)*AY4(J)+CZ2(I)*AZ4(J)) 
     *    +RTI3*RAJ*RCI 
      DBBK=RTIM*(BX2(I)*BX4(J)+BY2(I)*BY4(J)) 
     *    +RTI3*RBI*RBJ 
      DBCK=RTIM*(BX2(I)*CX4(J)+BY2(I)*CY4(J)) 
     *    +RTI3*RBI*RCJ 
      DCBK=RTIM*(CX2(I)*BX4(J)+CY2(I)*BY4(J)) 
     *    +RTI3*RBJ*RCI 
      DCCK=RTIM*(CX2(I)*CX4(J)+CY2(I)*CY4(J)+CZ2(I)*CZ4(J)) 
     *    +RTI3*RCI*RCJ 
      ENDIF
      IF (GTEST) THEN
      FDX=DDX+FDX 
      FDY=DDY+FDY 
      FDZ=DDZ+FDZ 
      FDAI=FDAI+DDAI 
      FDBI=FDBI+DDBI 
      FDCI=FDCI+DDCI 
      FDAJ=FDAJ+DDAJ 
      FDBJ=FDBJ+DDBJ 
      FDCJ=FDCJ+DDCJ 
      ENDIF
      IF (SSTEST) THEN
      SDXX=DXX+SDXX 
      SDYY=DYY+SDYY 
      SDZZ=DZZ+SDZZ 
      SDXY=DXY+SDXY 
      SDYZ=DYZ+SDYZ 
      SDZX=DZX+SDZX 
      SDXAI=DXAI+SDXAI 
      SDYAI=DYAI+SDYAI 
      SDZAI=DZAI+SDZAI 
      SDXBI=DXBI+SDXBI 
      SDYBI=DYBI+SDYBI 
      SDZBI=DZBI+SDZBI 
      SDXCI=DXCI+SDXCI 
      SDYCI=DYCI+SDYCI 
      SDZCI=DZCI+SDZCI 
      SDXAJ=DXAJ+SDXAJ 
      SDYAJ=DYAJ+SDYAJ 
      SDZAJ=DZAJ+SDZAJ 
      SDXBJ=DXBJ+SDXBJ 
      SDYBJ=DYBJ+SDYBJ 
      SDZBJ=DZBJ+SDZBJ 
      SDXCJ=DXCJ+SDXCJ 
      SDYCJ=DYCJ+SDYCJ 
      SDZCJ=DZCJ+SDZCJ 
      SDAAI=DAAI+SDAAI 
      SDABI=DABI+SDABI 
      SDACI=DACI+SDACI 
      SDBBI=DBBI+SDBBI 
      SDBCI=DBCI+SDBCI 
      SDCCI=DCCI+SDCCI 
      SDAAJ=DAAJ+SDAAJ 
      SDABJ=DABJ+SDABJ 
      SDACJ=DACJ+SDACJ 
      SDBBJ=DBBJ+SDBBJ 
      SDBCJ=DBCJ+SDBCJ 
      SDCCJ=DCCJ+SDCCJ 
      SDAAK=DAAK+SDAAK 
      SDABK=DABK+SDABK 
      SDBAK=DBAK+SDBAK 
      SDACK=DACK+SDACK 
      SDCAK=DCAK+SDCAK 
      SDBBK=DBBK+SDBBK 
      SDBCK=DBCK+SDBCK 
      SDCBK=DCBK+SDCBK 
      SDCCK=DCCK+SDCCK 
      ENDIF
      EPI=EPI-2.0D0*RI 
C ::::: @@@ 24 :::::::::::::::::::::::::::::::::::::::::::::::::::: 
      DX=DX4-XK1 
      DY=DY4-YK1 
      DZ=DZ4-ZK1 
      SX=DX*DX 
      SY=DY*DY 
      SZ=DZ*DZ 
      IF (SX+SY+SZ.EQ.0.0D0) THEN
         RSI=1.0D100
      ELSE
         RSI=1.0D0/(SX+SY+SZ) 
      ENDIF
      RI=SQRT(RSI) 
      IF (GTEST.OR.SSTEST) THEN
      RSIM=2.0D0*RSI 
      RTI=-RI*RSIM 
      RFI=RSI*RTI 
      XY=DX*DY 
      YZ=DY*DZ 
      ZX=DZ*DX 
      RTIM=-RTI 
      DDX=RTIM*DX 
      DDY=RTIM*DY 
      DDZ=RTIM*DZ 
      RAI=RI*(AX4(I)*DX+AY4(I)*DY+AZ4(I)*DZ) 
      DDAI=RSIM*RAI 
      RAJ=RI*(AX1(J)*DX+AY1(J)*DY+AZ1(J)*DZ) 
      DDAJ=RSIM*RAJ 
      RBI=RI*(BX4(I)*DX+BY4(I)*DY) 
      DDBI=RSIM*RBI 
      RBJ=RI*(BX1(J)*DX+BY1(J)*DY) 
      DDBJ=RSIM*RBJ 
      RCI=RI*(CX4(I)*DX+CY4(I)*DY+CZ4(I)*DZ) 
      DDCI=RSIM*RCI 
      RCJ=RI*(CX1(J)*DX+CY1(J)*DY+CZ1(J)*DZ) 
      DDCJ=RSIM*RCJ 
      ENDIF
      IF (SSTEST) THEN
      RFI3=RFI*3.0D0 
      DXX=RFI3*SX+RTIM 
      DYY=RFI3*SY+RTIM 
      DZZ=RFI3*SZ+RTIM 
      DXY=RFI3*XY 
      DYZ=RFI3*YZ 
      DZX=RFI3*ZX 
      RI3=3.0D0*RI 
      RI3X=DX*RI3 
      RI3Y=DY*RI3 
      RI3Z=DZ*RI3 
      DXAI=(RI3X*RAI-AX4(I))*RTI 
      DXAJ=(RI3X*RAJ-AX1(J))*RTI 
      DXBI=(RI3X*RBI-BX4(I))*RTI 
      DXBJ=(RI3X*RBJ-BX1(J))*RTI 
      DXCI=(RI3X*RCI-CX4(I))*RTI 
      DXCJ=(RI3X*RCJ-CX1(J))*RTI 
      DYAI=(RI3Y*RAI-AY4(I))*RTI 
      DYAJ=(RI3Y*RAJ-AY1(J))*RTI 
      DYBI=(RI3Y*RBI-BY4(I))*RTI 
      DYBJ=(RI3Y*RBJ-BY1(J))*RTI 
      DYCI=(RI3Y*RCI-CY4(I))*RTI 
      DYCJ=(RI3Y*RCJ-CY1(J))*RTI 
      DZAI=(RI3Z*RAI-AZ4(I))*RTI 
      DZAJ=(RI3Z*RAJ-AZ1(J))*RTI 
      DZBI=RI3Z*RBI*RTI 
      DZBJ=RI3Z*RBJ*RTI 
      DZCI=(RI3Z*RCI-CZ4(I))*RTI 
      DZCJ=(RI3Z*RCJ-CZ1(J))*RTI 
      RTI3=3.0D0*RTI 
      DAAI=RTIM*(AX4(I)*AX4(I)+AY4(I)*AY4(I)+AZ4(I)*AZ4(I) 
     *         +DX*AAX4(I)+DY*AAY4(I)+DZ*AAZ4(I)) 
     *    +RTI3*RAI*RAI 
      DAAJ=RTIM*(AX1(J)*AX1(J)+AY1(J)*AY1(J)+AZ1(J)*AZ1(J) 
     *         -DX*AAX1(J)-DY*AAY1(J)-DZ*AAZ1(J)) 
     *    +RTI3*RAJ*RAJ 
      DABI=RTIM*(AX4(I)*BX4(I)+AY4(I)*BY4(I) 
     *         +DX*ABX4(I)+DY*ABY4(I)) 
     *    +RTI3*RAI*RBI 
      DABJ=RTIM*(AX1(J)*BX1(J)+AY1(J)*BY1(J) 
     *         -DX*ABX1(J)-DY*ABY1(J)) 
     *    +RTI3*RAJ*RBJ 
      DACI=RTIM*(AX4(I)*CX4(I)+AY4(I)*CY4(I)+AZ4(I)*CZ4(I) 
     *         +DX*ACX4(I)+DY*ACY4(I)+DZ*ACZ4(I)) 
     *    +RTI3*RAI*RCI 
      DACJ=RTIM*(AX1(J)*CX1(J)+AY1(J)*CY1(J)+AZ1(J)*CZ1(J) 
     *         -DX*ACX1(J)-DY*ACY1(J)-DZ*ACZ1(J)) 
     *    +RTI3*RAJ*RCJ 
      DBBI=RTIM*(BX4(I)*BX4(I)+BY4(I)*BY4(I) 
     *         +DX*BBX4(I)+DY*BBY4(I)) 
     *    +RTI3*RBI*RBI 
      DBBJ=RTIM*(BX1(J)*BX1(J)+BY1(J)*BY1(J) 
     *         -DX*BBX1(J)-DY*BBY1(J)) 
     *    +RTI3*RBJ*RBJ 
      DBCI=RTIM*(BX4(I)*CX4(I)+BY4(I)*CY4(I) 
     *         +DX*BCX4(I)+DY*BCY4(I)) 
     *    +RTI3*RBI*RCI 
      DBCJ=RTIM*(BX1(J)*CX1(J)+BY1(J)*CY1(J) 
     *         -DX*BCX1(J)-DY*BCY1(J)) 
     *    +RTI3*RBJ*RCJ 
      DCCI=RTIM*(CX4(I)*CX4(I)+CY4(I)*CY4(I)+CZ4(I)*CZ4(I) 
     *         +DX*CCX4(I)+DY*CCY4(I)+DZ*CCZ4(I)) 
     *    +RTI3*RCI*RCI 
      DCCJ=RTIM*(CX1(J)*CX1(J)+CY1(J)*CY1(J)+CZ1(J)*CZ1(J) 
     *         -DX*CCX1(J)-DY*CCY1(J)-DZ*CCZ1(J)) 
     *    +RTI3*RCJ*RCJ 
      DAAK=RTIM*(AX4(I)*AX1(J)+AY4(I)*AY1(J)+AZ4(I)*AZ1(J)) 
     *    +RTI3*RAI*RAJ 
      DABK=RTIM*(AX4(I)*BX1(J)+AY4(I)*BY1(J)) 
     *    +RTI3*RAI*RBJ 
      DBAK=RTIM*(BX4(I)*AX1(J)+BY4(I)*AY1(J)) 
     *    +RTI3*RAJ*RBI 
      DACK=RTIM*(AX4(I)*CX1(J)+AY4(I)*CY1(J)+AZ4(I)*CZ1(J)) 
     *    +RTI3*RAI*RCJ 
      DCAK=RTIM*(CX4(I)*AX1(J)+CY4(I)*AY1(J)+CZ4(I)*AZ1(J)) 
     *    +RTI3*RAJ*RCI 
      DBBK=RTIM*(BX4(I)*BX1(J)+BY4(I)*BY1(J)) 
     *    +RTI3*RBI*RBJ 
      DBCK=RTIM*(BX4(I)*CX1(J)+BY4(I)*CY1(J)) 
     *    +RTI3*RBI*RCJ 
      DCBK=RTIM*(CX4(I)*BX1(J)+CY4(I)*BY1(J)) 
     *    +RTI3*RBJ*RCI 
      DCCK=RTIM*(CX4(I)*CX1(J)+CY4(I)*CY1(J)+CZ4(I)*CZ1(J)) 
     *    +RTI3*RCI*RCJ 
      ENDIF
      IF (GTEST) THEN
      FDX=DDX+FDX 
      FDY=DDY+FDY 
      FDZ=DDZ+FDZ 
      FDAI=FDAI+DDAI 
      FDBI=FDBI+DDBI 
      FDCI=FDCI+DDCI 
      FDAJ=FDAJ+DDAJ 
      FDBJ=FDBJ+DDBJ 
      FDCJ=FDCJ+DDCJ 
      ENDIF
      IF (SSTEST) THEN
      SDXX=DXX+SDXX 
      SDYY=DYY+SDYY 
      SDZZ=DZZ+SDZZ 
      SDXY=DXY+SDXY 
      SDYZ=DYZ+SDYZ 
      SDZX=DZX+SDZX 
      SDXAI=DXAI+SDXAI 
      SDYAI=DYAI+SDYAI 
      SDZAI=DZAI+SDZAI 
      SDXBI=DXBI+SDXBI 
      SDYBI=DYBI+SDYBI 
      SDZBI=DZBI+SDZBI 
      SDXCI=DXCI+SDXCI 
      SDYCI=DYCI+SDYCI 
      SDZCI=DZCI+SDZCI 
      SDXAJ=DXAJ+SDXAJ 
      SDYAJ=DYAJ+SDYAJ 
      SDZAJ=DZAJ+SDZAJ 
      SDXBJ=DXBJ+SDXBJ 
      SDYBJ=DYBJ+SDYBJ 
      SDZBJ=DZBJ+SDZBJ 
      SDXCJ=DXCJ+SDXCJ 
      SDYCJ=DYCJ+SDYCJ 
      SDZCJ=DZCJ+SDZCJ 
      SDAAI=DAAI+SDAAI 
      SDABI=DABI+SDABI 
      SDACI=DACI+SDACI 
      SDBBI=DBBI+SDBBI 
      SDBCI=DBCI+SDBCI 
      SDCCI=DCCI+SDCCI 
      SDAAJ=DAAJ+SDAAJ 
      SDABJ=DABJ+SDABJ 
      SDACJ=DACJ+SDACJ 
      SDBBJ=DBBJ+SDBBJ 
      SDBCJ=DBCJ+SDBCJ 
      SDCCJ=DCCJ+SDCCJ 
      SDAAK=DAAK+SDAAK 
      SDABK=DABK+SDABK 
      SDBAK=DBAK+SDBAK 
      SDACK=DACK+SDACK 
      SDCAK=DCAK+SDCAK 
      SDBBK=DBBK+SDBBK 
      SDBCK=DBCK+SDBCK 
      SDCBK=DCBK+SDCBK 
      SDCCK=DCCK+SDCCK 
      ENDIF
      EPI=EPI-2.0D0*RI 
C ::::: @@@ 41 :::::::::::::::::::::::::::::::::::::::::::::::::::: 
      DX=DX4-XK2 
      DY=DY4-YK2 
      DZ=DZ4-ZK2 
      SX=DX*DX 
      SY=DY*DY 
      SZ=DZ*DZ 
      IF (SX+SY+SZ.EQ.0.0D0) THEN
         RSI=1.0D100
      ELSE
         RSI=1.0D0/(SX+SY+SZ) 
      ENDIF
      RI=SQRT(RSI) 
      IF (GTEST.OR.SSTEST) THEN
      RSIM=2.0D0*RSI 
      RTI=-RI*RSIM 
      RFI=RSI*RTI 
      XY=DX*DY 
      YZ=DY*DZ 
      ZX=DZ*DX 
      RTIM=-RTI 
      DDX=RTIM*DX 
      DDY=RTIM*DY 
      DDZ=RTIM*DZ 
      RAI=RI*(AX4(I)*DX+AY4(I)*DY+AZ4(I)*DZ) 
      DDAI=RSIM*RAI 
      RAJ=RI*(AX2(J)*DX+AY2(J)*DY+AZ2(J)*DZ) 
      DDAJ=RSIM*RAJ 
      RBI=RI*(BX4(I)*DX+BY4(I)*DY) 
      DDBI=RSIM*RBI 
      RBJ=RI*(BX2(J)*DX+BY2(J)*DY) 
      DDBJ=RSIM*RBJ 
      RCI=RI*(CX4(I)*DX+CY4(I)*DY+CZ4(I)*DZ) 
      DDCI=RSIM*RCI 
      RCJ=RI*(CX2(J)*DX+CY2(J)*DY+CZ2(J)*DZ) 
      DDCJ=RSIM*RCJ 
      ENDIF
      IF (SSTEST) THEN
      RFI3=RFI*3.0D0 
      DXX=RFI3*SX+RTIM 
      DYY=RFI3*SY+RTIM 
      DZZ=RFI3*SZ+RTIM 
      DXY=RFI3*XY 
      DYZ=RFI3*YZ 
      DZX=RFI3*ZX 
      RI3=3.0D0*RI 
      RI3X=DX*RI3 
      RI3Y=DY*RI3 
      RI3Z=DZ*RI3 
      DXAI=(RI3X*RAI-AX4(I))*RTI 
      DXAJ=(RI3X*RAJ-AX2(J))*RTI 
      DXBI=(RI3X*RBI-BX4(I))*RTI 
      DXBJ=(RI3X*RBJ-BX2(J))*RTI 
      DXCI=(RI3X*RCI-CX4(I))*RTI 
      DXCJ=(RI3X*RCJ-CX2(J))*RTI 
      DYAI=(RI3Y*RAI-AY4(I))*RTI 
      DYAJ=(RI3Y*RAJ-AY2(J))*RTI 
      DYBI=(RI3Y*RBI-BY4(I))*RTI 
      DYBJ=(RI3Y*RBJ-BY2(J))*RTI 
      DYCI=(RI3Y*RCI-CY4(I))*RTI 
      DYCJ=(RI3Y*RCJ-CY2(J))*RTI 
      DZAI=(RI3Z*RAI-AZ4(I))*RTI 
      DZAJ=(RI3Z*RAJ-AZ2(J))*RTI 
      DZBI=RI3Z*RBI*RTI 
      DZBJ=RI3Z*RBJ*RTI 
      DZCI=(RI3Z*RCI-CZ4(I))*RTI 
      DZCJ=(RI3Z*RCJ-CZ2(J))*RTI 
      RTI3=3.0D0*RTI 
      DAAI=RTIM*(AX4(I)*AX4(I)+AY4(I)*AY4(I)+AZ4(I)*AZ4(I) 
     *         +DX*AAX4(I)+DY*AAY4(I)+DZ*AAZ4(I)) 
     *    +RTI3*RAI*RAI 
      DAAJ=RTIM*(AX2(J)*AX2(J)+AY2(J)*AY2(J)+AZ2(J)*AZ2(J) 
     *         -DX*AAX2(J)-DY*AAY2(J)-DZ*AAZ2(J)) 
     *    +RTI3*RAJ*RAJ 
      DABI=RTIM*(AX4(I)*BX4(I)+AY4(I)*BY4(I) 
     *         +DX*ABX4(I)+DY*ABY4(I)) 
     *    +RTI3*RAI*RBI 
      DABJ=RTIM*(AX2(J)*BX2(J)+AY2(J)*BY2(J) 
     *         -DX*ABX2(J)-DY*ABY2(J)) 
     *    +RTI3*RAJ*RBJ 
      DACI=RTIM*(AX4(I)*CX4(I)+AY4(I)*CY4(I)+AZ4(I)*CZ4(I) 
     *         +DX*ACX4(I)+DY*ACY4(I)+DZ*ACZ4(I)) 
     *    +RTI3*RAI*RCI 
      DACJ=RTIM*(AX2(J)*CX2(J)+AY2(J)*CY2(J)+AZ2(J)*CZ2(J) 
     *         -DX*ACX2(J)-DY*ACY2(J)-DZ*ACZ2(J)) 
     *    +RTI3*RAJ*RCJ 
      DBBI=RTIM*(BX4(I)*BX4(I)+BY4(I)*BY4(I) 
     *         +DX*BBX4(I)+DY*BBY4(I)) 
     *    +RTI3*RBI*RBI 
      DBBJ=RTIM*(BX2(J)*BX2(J)+BY2(J)*BY2(J) 
     *         -DX*BBX2(J)-DY*BBY2(J)) 
     *    +RTI3*RBJ*RBJ 
      DBCI=RTIM*(BX4(I)*CX4(I)+BY4(I)*CY4(I) 
     *         +DX*BCX4(I)+DY*BCY4(I)) 
     *    +RTI3*RBI*RCI 
      DBCJ=RTIM*(BX2(J)*CX2(J)+BY2(J)*CY2(J) 
     *         -DX*BCX2(J)-DY*BCY2(J)) 
     *    +RTI3*RBJ*RCJ 
      DCCI=RTIM*(CX4(I)*CX4(I)+CY4(I)*CY4(I)+CZ4(I)*CZ4(I) 
     *         +DX*CCX4(I)+DY*CCY4(I)+DZ*CCZ4(I)) 
     *    +RTI3*RCI*RCI 
      DCCJ=RTIM*(CX2(J)*CX2(J)+CY2(J)*CY2(J)+CZ2(J)*CZ2(J) 
     *         -DX*CCX2(J)-DY*CCY2(J)-DZ*CCZ2(J)) 
     *    +RTI3*RCJ*RCJ 
      DAAK=RTIM*(AX4(I)*AX2(J)+AY4(I)*AY2(J)+AZ4(I)*AZ2(J)) 
     *    +RTI3*RAI*RAJ 
      DABK=RTIM*(AX4(I)*BX2(J)+AY4(I)*BY2(J)) 
     *    +RTI3*RAI*RBJ 
      DBAK=RTIM*(BX4(I)*AX2(J)+BY4(I)*AY2(J)) 
     *    +RTI3*RAJ*RBI 
      DACK=RTIM*(AX4(I)*CX2(J)+AY4(I)*CY2(J)+AZ4(I)*CZ2(J)) 
     *    +RTI3*RAI*RCJ 
      DCAK=RTIM*(CX4(I)*AX2(J)+CY4(I)*AY2(J)+CZ4(I)*AZ2(J)) 
     *    +RTI3*RAJ*RCI 
      DBBK=RTIM*(BX4(I)*BX2(J)+BY4(I)*BY2(J)) 
     *    +RTI3*RBI*RBJ 
      DBCK=RTIM*(BX4(I)*CX2(J)+BY4(I)*CY2(J)) 
     *    +RTI3*RBI*RCJ 
      DCBK=RTIM*(CX4(I)*BX2(J)+CY4(I)*BY2(J)) 
     *    +RTI3*RBJ*RCI 
      DCCK=RTIM*(CX4(I)*CX2(J)+CY4(I)*CY2(J)+CZ4(I)*CZ2(J)) 
     *    +RTI3*RCI*RCJ 
      ENDIF
      IF (GTEST) THEN
      FDX=DDX+FDX 
      FDY=DDY+FDY 
      FDZ=DDZ+FDZ 
      FDAI=FDAI+DDAI 
      FDBI=FDBI+DDBI 
      FDCI=FDCI+DDCI 
      FDAJ=FDAJ+DDAJ 
      FDBJ=FDBJ+DDBJ 
      FDCJ=FDCJ+DDCJ 
      ENDIF
      IF (SSTEST) THEN
      SDXX=DXX+SDXX 
      SDYY=DYY+SDYY 
      SDZZ=DZZ+SDZZ 
      SDXY=DXY+SDXY 
      SDYZ=DYZ+SDYZ 
      SDZX=DZX+SDZX 
      SDXAI=DXAI+SDXAI 
      SDYAI=DYAI+SDYAI 
      SDZAI=DZAI+SDZAI 
      SDXBI=DXBI+SDXBI 
      SDYBI=DYBI+SDYBI 
      SDZBI=DZBI+SDZBI 
      SDXCI=DXCI+SDXCI 
      SDYCI=DYCI+SDYCI 
      SDZCI=DZCI+SDZCI 
      SDXAJ=DXAJ+SDXAJ 
      SDYAJ=DYAJ+SDYAJ 
      SDZAJ=DZAJ+SDZAJ 
      SDXBJ=DXBJ+SDXBJ 
      SDYBJ=DYBJ+SDYBJ 
      SDZBJ=DZBJ+SDZBJ 
      SDXCJ=DXCJ+SDXCJ 
      SDYCJ=DYCJ+SDYCJ 
      SDZCJ=DZCJ+SDZCJ 
      SDAAI=DAAI+SDAAI 
      SDABI=DABI+SDABI 
      SDACI=DACI+SDACI 
      SDBBI=DBBI+SDBBI 
      SDBCI=DBCI+SDBCI 
      SDCCI=DCCI+SDCCI 
      SDAAJ=DAAJ+SDAAJ 
      SDABJ=DABJ+SDABJ 
      SDACJ=DACJ+SDACJ 
      SDBBJ=DBBJ+SDBBJ 
      SDBCJ=DBCJ+SDBCJ 
      SDCCJ=DCCJ+SDCCJ 
      SDAAK=DAAK+SDAAK 
      SDABK=DABK+SDABK 
      SDBAK=DBAK+SDBAK 
      SDACK=DACK+SDACK 
      SDCAK=DCAK+SDCAK 
      SDBBK=DBBK+SDBBK 
      SDBCK=DBCK+SDBCK 
      SDCBK=DCBK+SDCBK 
      SDCCK=DCCK+SDCCK 
      ENDIF
      EPI=EPI-2.0D0*RI 
C ::::: @@@ 42 :::::::::::::::::::::::::::::::::::::::::::::::::::: 
      DX=DX4-XK4 
      DY=DY4-YK4 
      DZ=DZ4-ZK4 
      SX=DX*DX 
      SY=DY*DY 
      SZ=DZ*DZ 
      RSI=1.0D0/(SX+SY+SZ) 
      RI=SQRT(RSI) 
      IF (GTEST.OR.SSTEST) THEN
      RSIM=-4.0D0*RSI 
      RTI=-RI*RSIM 
      RFI=RSI*RTI 
      XY=DX*DY 
      YZ=DY*DZ 
      ZX=DZ*DX 
      RTIM=-RTI 
      DDX=RTIM*DX 
      DDY=RTIM*DY 
      DDZ=RTIM*DZ 
      RAI=RI*(AX4(I)*DX+AY4(I)*DY+AZ4(I)*DZ) 
      DDAI=RSIM*RAI 
      RAJ=RI*(AX4(J)*DX+AY4(J)*DY+AZ4(J)*DZ) 
      DDAJ=RSIM*RAJ 
      RBI=RI*(BX4(I)*DX+BY4(I)*DY) 
      DDBI=RSIM*RBI 
      RBJ=RI*(BX4(J)*DX+BY4(J)*DY) 
      DDBJ=RSIM*RBJ 
      RCI=RI*(CX4(I)*DX+CY4(I)*DY+CZ4(I)*DZ) 
      DDCI=RSIM*RCI 
      RCJ=RI*(CX4(J)*DX+CY4(J)*DY+CZ4(J)*DZ) 
      DDCJ=RSIM*RCJ 
      ENDIF
      IF (SSTEST) THEN
      RFI3=RFI*3.0D0 
      DXX=RFI3*SX+RTIM 
      DYY=RFI3*SY+RTIM 
      DZZ=RFI3*SZ+RTIM 
      DXY=RFI3*XY 
      DYZ=RFI3*YZ 
      DZX=RFI3*ZX 
      RI3=3.0D0*RI 
      RI3X=DX*RI3 
      RI3Y=DY*RI3 
      RI3Z=DZ*RI3 
      DXAI=(RI3X*RAI-AX4(I))*RTI 
      DXAJ=(RI3X*RAJ-AX4(J))*RTI 
      DXBI=(RI3X*RBI-BX4(I))*RTI 
      DXBJ=(RI3X*RBJ-BX4(J))*RTI 
      DXCI=(RI3X*RCI-CX4(I))*RTI 
      DXCJ=(RI3X*RCJ-CX4(J))*RTI 
      DYAI=(RI3Y*RAI-AY4(I))*RTI 
      DYAJ=(RI3Y*RAJ-AY4(J))*RTI 
      DYBI=(RI3Y*RBI-BY4(I))*RTI 
      DYBJ=(RI3Y*RBJ-BY4(J))*RTI 
      DYCI=(RI3Y*RCI-CY4(I))*RTI 
      DYCJ=(RI3Y*RCJ-CY4(J))*RTI 
      DZAI=(RI3Z*RAI-AZ4(I))*RTI 
      DZAJ=(RI3Z*RAJ-AZ4(J))*RTI 
      DZBI=RI3Z*RBI*RTI 
      DZBJ=RI3Z*RBJ*RTI 
      DZCI=(RI3Z*RCI-CZ4(I))*RTI 
      DZCJ=(RI3Z*RCJ-CZ4(J))*RTI 
      RTI3=3.0D0*RTI 
      DAAI=RTIM*(AX4(I)*AX4(I)+AY4(I)*AY4(I)+AZ4(I)*AZ4(I) 
     *         +DX*AAX4(I)+DY*AAY4(I)+DZ*AAZ4(I)) 
     *    +RTI3*RAI*RAI 
      DAAJ=RTIM*(AX4(J)*AX4(J)+AY4(J)*AY4(J)+AZ4(J)*AZ4(J) 
     *         -DX*AAX4(J)-DY*AAY4(J)-DZ*AAZ4(J)) 
     *    +RTI3*RAJ*RAJ 
      DABI=RTIM*(AX4(I)*BX4(I)+AY4(I)*BY4(I) 
     *         +DX*ABX4(I)+DY*ABY4(I)) 
     *    +RTI3*RAI*RBI 
      DABJ=RTIM*(AX4(J)*BX4(J)+AY4(J)*BY4(J) 
     *         -DX*ABX4(J)-DY*ABY4(J)) 
     *    +RTI3*RAJ*RBJ 
      DACI=RTIM*(AX4(I)*CX4(I)+AY4(I)*CY4(I)+AZ4(I)*CZ4(I) 
     *         +DX*ACX4(I)+DY*ACY4(I)+DZ*ACZ4(I)) 
     *    +RTI3*RAI*RCI 
      DACJ=RTIM*(AX4(J)*CX4(J)+AY4(J)*CY4(J)+AZ4(J)*CZ4(J) 
     *         -DX*ACX4(J)-DY*ACY4(J)-DZ*ACZ4(J)) 
     *    +RTI3*RAJ*RCJ 
      DBBI=RTIM*(BX4(I)*BX4(I)+BY4(I)*BY4(I) 
     *         +DX*BBX4(I)+DY*BBY4(I)) 
     *    +RTI3*RBI*RBI 
      DBBJ=RTIM*(BX4(J)*BX4(J)+BY4(J)*BY4(J) 
     *         -DX*BBX4(J)-DY*BBY4(J)) 
     *    +RTI3*RBJ*RBJ 
      DBCI=RTIM*(BX4(I)*CX4(I)+BY4(I)*CY4(I) 
     *         +DX*BCX4(I)+DY*BCY4(I)) 
     *    +RTI3*RBI*RCI 
      DBCJ=RTIM*(BX4(J)*CX4(J)+BY4(J)*CY4(J) 
     *         -DX*BCX4(J)-DY*BCY4(J)) 
     *    +RTI3*RBJ*RCJ 
      DCCI=RTIM*(CX4(I)*CX4(I)+CY4(I)*CY4(I)+CZ4(I)*CZ4(I) 
     *         +DX*CCX4(I)+DY*CCY4(I)+DZ*CCZ4(I)) 
     *    +RTI3*RCI*RCI 
      DCCJ=RTIM*(CX4(J)*CX4(J)+CY4(J)*CY4(J)+CZ4(J)*CZ4(J) 
     *         -DX*CCX4(J)-DY*CCY4(J)-DZ*CCZ4(J)) 
     *    +RTI3*RCJ*RCJ 
      DAAK=RTIM*(AX4(I)*AX4(J)+AY4(I)*AY4(J)+AZ4(I)*AZ4(J)) 
     *    +RTI3*RAI*RAJ 
      DABK=RTIM*(AX4(I)*BX4(J)+AY4(I)*BY4(J)) 
     *    +RTI3*RAI*RBJ 
      DBAK=RTIM*(BX4(I)*AX4(J)+BY4(I)*AY4(J)) 
     *    +RTI3*RAJ*RBI 
      DACK=RTIM*(AX4(I)*CX4(J)+AY4(I)*CY4(J)+AZ4(I)*CZ4(J)) 
     *    +RTI3*RAI*RCJ 
      DCAK=RTIM*(CX4(I)*AX4(J)+CY4(I)*AY4(J)+CZ4(I)*AZ4(J)) 
     *    +RTI3*RAJ*RCI 
      DBBK=RTIM*(BX4(I)*BX4(J)+BY4(I)*BY4(J)) 
     *    +RTI3*RBI*RBJ 
      DBCK=RTIM*(BX4(I)*CX4(J)+BY4(I)*CY4(J)) 
     *    +RTI3*RBI*RCJ 
      DCBK=RTIM*(CX4(I)*BX4(J)+CY4(I)*BY4(J)) 
     *    +RTI3*RBJ*RCI 
      DCCK=RTIM*(CX4(I)*CX4(J)+CY4(I)*CY4(J)+CZ4(I)*CZ4(J)) 
     *    +RTI3*RCI*RCJ 
      ENDIF
      IF (GTEST) THEN
      FDX=DDX+FDX 
      FDY=DDY+FDY 
      FDZ=DDZ+FDZ 
      FDAI=FDAI+DDAI 
      FDBI=FDBI+DDBI 
      FDCI=FDCI+DDCI 
      FDAJ=FDAJ+DDAJ 
      FDBJ=FDBJ+DDBJ 
      FDCJ=FDCJ+DDCJ 
      ENDIF
      IF (SSTEST) THEN
      SDXX=DXX+SDXX 
      SDYY=DYY+SDYY 
      SDZZ=DZZ+SDZZ 
      SDXY=DXY+SDXY 
      SDYZ=DYZ+SDYZ 
      SDZX=DZX+SDZX 
      SDXAI=DXAI+SDXAI 
      SDYAI=DYAI+SDYAI 
      SDZAI=DZAI+SDZAI 
      SDXBI=DXBI+SDXBI 
      SDYBI=DYBI+SDYBI 
      SDZBI=DZBI+SDZBI 
      SDXCI=DXCI+SDXCI 
      SDYCI=DYCI+SDYCI 
      SDZCI=DZCI+SDZCI 
      SDXAJ=DXAJ+SDXAJ 
      SDYAJ=DYAJ+SDYAJ 
      SDZAJ=DZAJ+SDZAJ 
      SDXBJ=DXBJ+SDXBJ 
      SDYBJ=DYBJ+SDYBJ 
      SDZBJ=DZBJ+SDZBJ 
      SDXCJ=DXCJ+SDXCJ 
      SDYCJ=DYCJ+SDYCJ 
      SDZCJ=DZCJ+SDZCJ 
      SDAAI=DAAI+SDAAI 
      SDABI=DABI+SDABI 
      SDACI=DACI+SDACI 
      SDBBI=DBBI+SDBBI 
      SDBCI=DBCI+SDBCI 
      SDCCI=DCCI+SDCCI 
      SDAAJ=DAAJ+SDAAJ 
      SDABJ=DABJ+SDABJ 
      SDACJ=DACJ+SDACJ 
      SDBBJ=DBBJ+SDBBJ 
      SDBCJ=DBCJ+SDBCJ 
      SDCCJ=DCCJ+SDCCJ 
      SDAAK=DAAK+SDAAK 
      SDABK=DABK+SDABK 
      SDBAK=DBAK+SDBAK 
      SDACK=DACK+SDACK 
      SDCAK=DCAK+SDCAK 
      SDBBK=DBBK+SDBBK 
      SDBCK=DBCK+SDBCK 
      SDCBK=DCBK+SDCBK 
      SDCCK=DCCK+SDCCK 
      ENDIF
      EPI=EPI+4.0D0*RI 
C ::::: @@@ 44 :::::::::::::::::::::::::::::::::::::::::::::::::::: 
      DX=DX3-XK3 
      DY=DY3-YK3 
      DZ=DZ3-ZK3 
      SX=DX*DX 
      SY=DY*DY 
      SZ=DZ*DZ 
      RSI=1.0D0/(SX+SY+SZ) 
      RI=SQRT(RSI) 
      RHI=RSI*RSI*RSI 
      RDI=RHI*RHI 
      DE=P1*RDI+P2*RHI 
      IF (GTEST.OR.SSTEST) THEN
      DF=(P3*RDI+P4*RHI)*RI 
      DS=(P5*RDI+P6*RHI)*RSI 
      XY=DX*DY 
      YZ=DY*DZ 
      ZX=DZ*DX 
      DFRI=DF*RI 
      DDX=DFRI*DX 
      DDY=DFRI*DY 
      DDZ=DFRI*DZ 
      RAI=RI*(AX3(I)*DX+AY3(I)*DY+AZ3(I)*DZ) 
      DDAI=DF*RAI 
      RAJ=RI*(AX3(J)*DX+AY3(J)*DY+AZ3(J)*DZ) 
      DDAJ=DF*RAJ 
      RBI=RI*(BX3(I)*DX+BY3(I)*DY) 
      DDBI=DF*RBI 
      RBJ=RI*(BX3(J)*DX+BY3(J)*DY) 
      DDBJ=DF*RBJ 
      RCI=RI*(CX3(I)*DX+CY3(I)*DY+CZ3(I)*DZ) 
      DDCI=DF*RCI 
      RCJ=RI*(CX3(J)*DX+CY3(J)*DY+CZ3(J)*DZ) 
      DDCJ=DF*RCJ 
      ENDIF
      IF (SSTEST) THEN
      DSF=DS-DFRI 
      DSFRSI=DSF*RSI 
      DSFRI=DSF*RI 
      DXX=DFRI+DSFRSI*SX 
      DYY=DFRI+DSFRSI*SY 
      DZZ=DFRI+DSFRSI*SZ 
      DXY=DSFRSI*XY 
      DYZ=DSFRSI*YZ 
      DZX=DSFRSI*ZX 
      DXDSF=DX*DSFRI 
      DYDSF=DY*DSFRI 
      DZDSF=DZ*DSFRI 
      DXAI=DFRI*AX3(I)+RAI*DXDSF 
      DXAJ=DFRI*AX3(J)+RAJ*DXDSF 
      DXBI=DFRI*BX3(I)+RBI*DXDSF 
      DXBJ=DFRI*BX3(J)+RBJ*DXDSF 
      DXCI=DFRI*CX3(I)+RCI*DXDSF 
      DXCJ=DFRI*CX3(J)+RCJ*DXDSF 
      DYAI=DFRI*AY3(I)+RAI*DYDSF 
      DYAJ=DFRI*AY3(J)+RAJ*DYDSF 
      DYBI=DFRI*BY3(I)+RBI*DYDSF 
      DYBJ=DFRI*BY3(J)+RBJ*DYDSF 
      DYCI=DFRI*CY3(I)+RCI*DYDSF 
      DYCJ=DFRI*CY3(J)+RCJ*DYDSF 
      DZAI=DFRI*AZ3(I)+RAI*DZDSF 
      DZAJ=DFRI*AZ3(J)+RAJ*DZDSF 
      DZBI=RBI*DZDSF 
      DZBJ=RBJ*DZDSF 
      DZCI=DFRI*CZ3(I)+RCI*DZDSF 
      DZCJ=DFRI*CZ3(J)+RCJ*DZDSF 
      DAAI=DFRI*(AX3(I)*AX3(I)+AY3(I)*AY3(I)+AZ3(I)*AZ3(I) 
     *         +DX*AAX3(I)+DY*AAY3(I)+DZ*AAZ3(I)) 
     *    +DSF*RAI*RAI 
      DAAJ=DFRI*(AX3(J)*AX3(J)+AY3(J)*AY3(J)+AZ3(J)*AZ3(J) 
     *         -DX*AAX3(J)-DY*AAY3(J)-DZ*AAZ3(J)) 
     *    +DSF*RAJ*RAJ 
      DABI=DFRI*(AX3(I)*BX3(I)+AY3(I)*BY3(I) 
     *         +DX*ABX3(I)+DY*ABY3(I)) 
     *    +DSF*RAI*RBI 
      DABJ=DFRI*(AX3(J)*BX3(J)+AY3(J)*BY3(J) 
     *         -DX*ABX3(J)-DY*ABY3(J)) 
     *    +DSF*RAJ*RBJ 
      DBBI=DFRI*(BX3(I)*BX3(I)+BY3(I)*BY3(I) 
     *         +DX*BBX3(I)+DY*BBY3(I)) 
     *    +DSF*RBI*RBI 
      DBBJ=DFRI*(BX3(J)*BX3(J)+BY3(J)*BY3(J) 
     *         -DX*BBX3(J)-DY*BBY3(J)) 
     *    +DSF*RBJ*RBJ 
      DCCI=DFRI*(CX3(I)*CX3(I)+CY3(I)*CY3(I)+CZ3(I)*CZ3(I) 
     *         +DX*CCX3(I)+DY*CCY3(I)+DZ*CCZ3(I)) 
     *    +DSF*RCI*RCI 
      DCCJ=DFRI*(CX3(J)*CX3(J)+CY3(J)*CY3(J)+CZ3(J)*CZ3(J) 
     *         -DX*CCX3(J)-DY*CCY3(J)-DZ*CCZ3(J)) 
     *    +DSF*RCJ*RCJ 
      DACI=DFRI*(AX3(I)*CX3(I)+AY3(I)*CY3(I)+AZ3(I)*CZ3(I) 
     *         +DX*ACX3(I)+DY*ACY3(I)+DZ*ACZ3(I)) 
     *    +DSF*RAI*RCI 
      DACJ=DFRI*(AX3(J)*CX3(J)+AY3(J)*CY3(J)+AZ3(J)*CZ3(J) 
     *         -DX*ACX3(J)-DY*ACY3(J)-DZ*ACZ3(J)) 
     *    +DSF*RAJ*RCJ 
      DBCI=DFRI*(BX3(I)*CX3(I)+BY3(I)*CY3(I) 
     *         +DX*BCX3(I)+DY*BCY3(I)) 
     *    +DSF*RBI*RCI 
      DBCJ=DFRI*(BX3(J)*CX3(J)+BY3(J)*CY3(J) 
     *         -DX*BCX3(J)-DY*BCY3(J)) 
     *    +DSF*RBJ*RCJ 
      DAAK=DFRI*(AX3(I)*AX3(J)+AY3(I)*AY3(J)+AZ3(I)*AZ3(J)) 
     *    +DSF*RAI*RAJ 
      DABK=DFRI*(AX3(I)*BX3(J)+AY3(I)*BY3(J)) 
     *    +DSF*RAI*RBJ 
      DBAK=DFRI*(BX3(I)*AX3(J)+BY3(I)*AY3(J)) 
     *    +DSF*RAJ*RBI 
      DBBK=DFRI*(BX3(I)*BX3(J)+BY3(I)*BY3(J)) 
     *    +DSF*RBI*RBJ 
      DCCK=DFRI*(CX3(I)*CX3(J)+CY3(I)*CY3(J)+CZ3(I)*CZ3(J)) 
     *    +DSF*RCI*RCJ 
      DACK=DFRI*(AX3(I)*CX3(J)+AY3(I)*CY3(J)+AZ3(I)*CZ3(J)) 
     *    +DSF*RAI*RCJ 
      DCAK=DFRI*(CX3(I)*AX3(J)+CY3(I)*AY3(J)+CZ3(I)*AZ3(J)) 
     *    +DSF*RAJ*RCI 
      DBCK=DFRI*(BX3(I)*CX3(J)+BY3(I)*CY3(J)) 
     *    +DSF*RBI*RCJ 
      DCBK=DFRI*(CX3(I)*BX3(J)+CY3(I)*BY3(J)) 
     *    +DSF*RBJ*RCI 
      ENDIF
      IF (GTEST) THEN
      FDX=DDX+FDX 
      FDY=DDY+FDY 
      FDZ=DDZ+FDZ 
      FDAI=FDAI+DDAI 
      FDBI=FDBI+DDBI 
      FDCI=FDCI+DDCI 
      FDAJ=FDAJ+DDAJ 
      FDBJ=FDBJ+DDBJ 
      FDCJ=FDCJ+DDCJ 
      ENDIF
      IF (SSTEST) THEN
      SDXX=DXX+SDXX 
      SDYY=DYY+SDYY 
      SDZZ=DZZ+SDZZ 
      SDXY=DXY+SDXY 
      SDYZ=DYZ+SDYZ 
      SDZX=DZX+SDZX 
      SDXAI=DXAI+SDXAI 
      SDYAI=DYAI+SDYAI 
      SDZAI=DZAI+SDZAI 
      SDXBI=DXBI+SDXBI 
      SDYBI=DYBI+SDYBI 
      SDZBI=DZBI+SDZBI 
      SDXCI=DXCI+SDXCI 
      SDYCI=DYCI+SDYCI 
      SDZCI=DZCI+SDZCI 
      SDXAJ=DXAJ+SDXAJ 
      SDYAJ=DYAJ+SDYAJ 
      SDZAJ=DZAJ+SDZAJ 
      SDXBJ=DXBJ+SDXBJ 
      SDYBJ=DYBJ+SDYBJ 
      SDZBJ=DZBJ+SDZBJ 
      SDXCJ=DXCJ+SDXCJ 
      SDYCJ=DYCJ+SDYCJ 
      SDZCJ=DZCJ+SDZCJ 
      SDAAI=DAAI+SDAAI 
      SDABI=DABI+SDABI 
      SDACI=DACI+SDACI 
      SDBBI=DBBI+SDBBI 
      SDBCI=DBCI+SDBCI 
      SDCCI=DCCI+SDCCI 
      SDAAJ=DAAJ+SDAAJ 
      SDABJ=DABJ+SDABJ 
      SDACJ=DACJ+SDACJ 
      SDBBJ=DBBJ+SDBBJ 
      SDBCJ=DBCJ+SDBCJ 
      SDCCJ=DCCJ+SDCCJ 
      SDAAK=DAAK+SDAAK 
      SDABK=DABK+SDABK 
      SDBAK=DBAK+SDBAK 
      SDACK=DACK+SDACK 
      SDCAK=DCAK+SDCAK 
      SDBBK=DBBK+SDBBK 
      SDBCK=DBCK+SDBCK 
      SDCBK=DCBK+SDCBK 
      SDCCK=DCCK+SDCCK 
      ENDIF
      EPI=EPI+DE 
C ::::: @@@ 33 :::::::::::::::::::::::::::::::::::::::::::::::::::: 
      NTJ1=3*J-2 
      NRJ1=N3+NTJ1 
      NRJ2=NRJ1+1 
      NRJ3=NRJ2+1 
      NTJ2=NTJ1+1 
      NTJ3=NTJ2+1 
      IF(RR.GT.RL) THEN 
         U=EPI 
         EPI=EPI*SF 
         IF (GTEST) THEN
         VDX=FDX 
         VDY=FDY 
         VDZ=FDZ 
         FDX=VDX*SF+SDX*U 
         FDY=VDY*SF+SDY*U 
         FDZ=VDZ*SF+SDZ*U 
         FDAI=FDAI*SF 
         FDBI=FDBI*SF 
         FDCI=FDCI*SF 
         FDAJ=FDAJ*SF 
         FDBJ=FDBJ*SF 
         FDCJ=FDCJ*SF 
         ENDIF
         IF (SSTEST) THEN
         SDXX=SDXX*SF+2.0D0*VDX*SDX+SXX*U 
         SDYY=SDYY*SF+2.0D0*VDY*SDY+SYY*U 
         SDZZ=SDZZ*SF+2.0D0*VDZ*SDZ+SZZ*U 
         SDXY=SDXY*SF+VDX*SDY+VDY*SDX+SXY*U 
         SDYZ=SDYZ*SF+VDY*SDZ+VDZ*SDY+SYZ*U 
         SDZX=SDZX*SF+VDZ*SDX+VDX*SDZ+SZX*U 
         SDXAI=SDX*FDAI+SDXAI*SF 
         SDXBI=SDX*FDBI+SDXBI*SF 
         SDXCI=SDX*FDCI+SDXCI*SF 
         SDYAI=SDY*FDAI+SDYAI*SF 
         SDYBI=SDY*FDBI+SDYBI*SF 
         SDYCI=SDY*FDCI+SDYCI*SF 
         SDZAI=SDZ*FDAI+SDZAI*SF 
         SDZBI=SDZ*FDBI+SDZBI*SF 
         SDZCI=SDZ*FDCI+SDZCI*SF 
         SDXAJ=SDX*FDAJ+SDXAJ*SF 
         SDXBJ=SDX*FDBJ+SDXBJ*SF 
         SDXCJ=SDX*FDCJ+SDXCJ*SF 
         SDYAJ=SDY*FDAJ+SDYAJ*SF 
         SDYBJ=SDY*FDBJ+SDYBJ*SF 
         SDYCJ=SDY*FDCJ+SDYCJ*SF 
         SDZAJ=SDZ*FDAJ+SDZAJ*SF 
         SDZBJ=SDZ*FDBJ+SDZBJ*SF 
         SDZCJ=SDZ*FDCJ+SDZCJ*SF 
         SDAAI=SDAAI*SF 
         SDABI=SDABI*SF 
         SDACI=SDACI*SF 
         SDBBI=SDBBI*SF 
         SDBCI=SDBCI*SF 
         SDCCI=SDCCI*SF 
         SDAAJ=SDAAJ*SF 
         SDABJ=SDABJ*SF 
         SDACJ=SDACJ*SF 
         SDBBJ=SDBBJ*SF 
         SDBCJ=SDBCJ*SF 
         SDCCJ=SDCCJ*SF 
         SDAAK=SDAAK*SF 
         SDABK=SDABK*SF 
         SDACK=SDACK*SF 
         SDBBK=SDBBK*SF 
         SDBCK=SDBCK*SF 
         SDCCK=SDCCK*SF 
         ENDIF
      ENDIF
      IF (GTEST) THEN
      FX(I)=FX(I)+FDX
      FY(I)=FY(I)+FDY
      FZ(I)=FZ(I)+FDZ
      FX(J)=FX(J)-FDX
      FY(J)=FY(J)-FDY
      FZ(J)=FZ(J)-FDZ
      FA(I)=FA(I)+FDAI
      FB(I)=FB(I)+FDBI
      FC(I)=FC(I)+FDCI
      FA(J)=FA(J)-FDAJ
      FB(J)=FB(J)-FDBJ
      FC(J)=FC(J)-FDCJ
      ENDIF
      EPID=EPI*EMUPID 
C
C  Total Energy
C
      EP=EP+EPID 
C
C translation - translation for molecule I, diagonal block
C
      IF (SSTEST) THEN
      HESS(NTI1,NTI1)=SDXX+HESS(NTI1,NTI1) 
      HESS(NTI1,NTI2)=SDXY+HESS(NTI1,NTI2) 
      HESS(NTI2,NTI1)=SDXY+HESS(NTI2,NTI1) 
      HESS(NTI1,NTI3)=SDZX+HESS(NTI1,NTI3) 
      HESS(NTI3,NTI1)=SDZX+HESS(NTI3,NTI1) 
      HESS(NTI2,NTI2)=SDYY+HESS(NTI2,NTI2) 
      HESS(NTI2,NTI3)=SDYZ+HESS(NTI2,NTI3) 
      HESS(NTI3,NTI2)=SDYZ+HESS(NTI3,NTI2) 
      HESS(NTI3,NTI3)=SDZZ+HESS(NTI3,NTI3) 
C
C translation - translation for molecule J, diagonal block
C
      HESS(NTJ1,NTJ1)=SDXX+HESS(NTJ1,NTJ1) 
      HESS(NTJ1,NTJ2)=SDXY+HESS(NTJ1,NTJ2) 
      HESS(NTJ2,NTJ1)=SDXY+HESS(NTJ2,NTJ1) 
      HESS(NTJ1,NTJ3)=SDZX+HESS(NTJ1,NTJ3) 
      HESS(NTJ3,NTJ1)=SDZX+HESS(NTJ3,NTJ1) 
      HESS(NTJ2,NTJ2)=SDYY+HESS(NTJ2,NTJ2) 
      HESS(NTJ2,NTJ3)=SDYZ+HESS(NTJ2,NTJ3) 
      HESS(NTJ3,NTJ2)=SDYZ+HESS(NTJ3,NTJ2) 
      HESS(NTJ3,NTJ3)=SDZZ+HESS(NTJ3,NTJ3) 
C
C  Translation on I with translation on J, off-diagonal
C
      HESS(NTJ1,NTI1)=-SDXX+HESS(NTJ1,NTI1) 
      HESS(NTJ1,NTI2)=-SDXY+HESS(NTJ1,NTI2) 
      HESS(NTJ2,NTI1)=-SDXY+HESS(NTJ2,NTI1) 
      HESS(NTJ1,NTI3)=-SDZX+HESS(NTJ1,NTI3) 
      HESS(NTJ3,NTI1)=-SDZX+HESS(NTJ3,NTI1) 
      HESS(NTJ2,NTI2)=-SDYY+HESS(NTJ2,NTI2) 
      HESS(NTJ2,NTI3)=-SDYZ+HESS(NTJ2,NTI3) 
      HESS(NTJ3,NTI2)=-SDYZ+HESS(NTJ3,NTI2) 
      HESS(NTJ3,NTI3)=-SDZZ+HESS(NTJ3,NTI3) 
C
C  Translation on J with translation on I, off-diagonal
C
      HESS(NTI1,NTJ1)=-SDXX+HESS(NTI1,NTJ1) 
      HESS(NTI1,NTJ2)=-SDXY+HESS(NTI1,NTJ2) 
      HESS(NTI2,NTJ1)=-SDXY+HESS(NTI2,NTJ1) 
      HESS(NTI1,NTJ3)=-SDZX+HESS(NTI1,NTJ3) 
      HESS(NTI3,NTJ1)=-SDZX+HESS(NTI3,NTJ1) 
      HESS(NTI2,NTJ2)=-SDYY+HESS(NTI2,NTJ2) 
      HESS(NTI2,NTJ3)=-SDYZ+HESS(NTI2,NTJ3) 
      HESS(NTI3,NTJ2)=-SDYZ+HESS(NTI3,NTJ2) 
      HESS(NTI3,NTJ3)=-SDZZ+HESS(NTI3,NTJ3) 
C
C  Rotation on I with rotation on I, diagonal
C
      HESS(NRI1,NRI1)=SDAAI+HESS(NRI1,NRI1) 
      HESS(NRI1,NRI2)=SDABI+HESS(NRI1,NRI2) 
      HESS(NRI2,NRI1)=SDABI+HESS(NRI2,NRI1) 
      HESS(NRI1,NRI3)=SDACI+HESS(NRI1,NRI3) 
      HESS(NRI3,NRI1)=SDACI+HESS(NRI3,NRI1) 
      HESS(NRI2,NRI2)=SDBBI+HESS(NRI2,NRI2) 
      HESS(NRI2,NRI3)=SDBCI+HESS(NRI2,NRI3) 
      HESS(NRI3,NRI2)=SDBCI+HESS(NRI3,NRI2) 
      HESS(NRI3,NRI3)=SDCCI+HESS(NRI3,NRI3) 
C
C  Rotation on J with rotation on J, diagonal
C
      HESS(NRJ1,NRJ1)=SDAAJ+HESS(NRJ1,NRJ1) 
      HESS(NRJ1,NRJ2)=SDABJ+HESS(NRJ1,NRJ2) 
      HESS(NRJ2,NRJ1)=SDABJ+HESS(NRJ2,NRJ1) 
      HESS(NRJ1,NRJ3)=SDACJ+HESS(NRJ1,NRJ3) 
      HESS(NRJ3,NRJ1)=SDACJ+HESS(NRJ3,NRJ1) 
      HESS(NRJ2,NRJ2)=SDBBJ+HESS(NRJ2,NRJ2) 
      HESS(NRJ2,NRJ3)=SDBCJ+HESS(NRJ2,NRJ3) 
      HESS(NRJ3,NRJ2)=SDBCJ+HESS(NRJ3,NRJ2) 
      HESS(NRJ3,NRJ3)=SDCCJ+HESS(NRJ3,NRJ3) 
C
C  Rotation on I with rotation on J, off-diagonal
C
      HESS(NRI1,NRJ1)=-SDAAK+HESS(NRI1,NRJ1) 
      HESS(NRJ1,NRI1)=-SDAAK+HESS(NRJ1,NRI1) 
      HESS(NRI1,NRJ2)=-SDABK+HESS(NRI1,NRJ2) 
      HESS(NRJ2,NRI1)=-SDABK+HESS(NRJ2,NRI1) 
      HESS(NRI2,NRJ1)=-SDBAK+HESS(NRI2,NRJ1) 
      HESS(NRJ1,NRI2)=-SDBAK+HESS(NRJ1,NRI2) 
      HESS(NRI1,NRJ3)=-SDACK+HESS(NRI1,NRJ3) 
      HESS(NRJ3,NRI1)=-SDACK+HESS(NRJ3,NRI1) 
      HESS(NRI3,NRJ1)=-SDCAK+HESS(NRI3,NRJ1) 
C
C  Rotation on J with rotation on I, off-diagonal
C
      HESS(NRJ1,NRI3)=-SDCAK+HESS(NRJ1,NRI3) 
      HESS(NRI2,NRJ2)=-SDBBK+HESS(NRI2,NRJ2) 
      HESS(NRJ2,NRI2)=-SDBBK+HESS(NRJ2,NRI2) 
      HESS(NRI2,NRJ3)=-SDBCK+HESS(NRI2,NRJ3) 
      HESS(NRJ3,NRI2)=-SDBCK+HESS(NRJ3,NRI2) 
      HESS(NRI3,NRJ2)=-SDCBK+HESS(NRI3,NRJ2) 
      HESS(NRJ2,NRI3)=-SDCBK+HESS(NRJ2,NRI3) 
      HESS(NRI3,NRJ3)=-SDCCK+HESS(NRI3,NRJ3) 
      HESS(NRJ3,NRI3)=-SDCCK+HESS(NRJ3,NRI3) 
C
C  Translation on I with rotation on I, off-diagonal
C
      HESS(NTI1,NRI1)=SDXAI+HESS(NTI1,NRI1) 
      HESS(NRI1,NTI1)=SDXAI+HESS(NRI1,NTI1) 
      HESS(NTI1,NRI2)=SDXBI+HESS(NTI1,NRI2) 
      HESS(NRI2,NTI1)=SDXBI+HESS(NRI2,NTI1) 
      HESS(NTI1,NRI3)=SDXCI+HESS(NTI1,NRI3) 
      HESS(NRI3,NTI1)=SDXCI+HESS(NRI3,NTI1) 
      HESS(NTI2,NRI1)=SDYAI+HESS(NTI2,NRI1) 
      HESS(NRI1,NTI2)=SDYAI+HESS(NRI1,NTI2) 
      HESS(NTI2,NRI2)=SDYBI+HESS(NTI2,NRI2) 
      HESS(NRI2,NTI2)=SDYBI+HESS(NRI2,NTI2) 
      HESS(NTI2,NRI3)=SDYCI+HESS(NTI2,NRI3) 
      HESS(NRI3,NTI2)=SDYCI+HESS(NRI3,NTI2) 
      HESS(NTI3,NRI1)=SDZAI+HESS(NTI3,NRI1) 
      HESS(NRI1,NTI3)=SDZAI+HESS(NRI1,NTI3) 
      HESS(NTI3,NRI2)=SDZBI+HESS(NTI3,NRI2) 
      HESS(NRI2,NTI3)=SDZBI+HESS(NRI2,NTI3) 
      HESS(NTI3,NRI3)=SDZCI+HESS(NTI3,NRI3) 
      HESS(NRI3,NTI3)=SDZCI+HESS(NRI3,NTI3) 
C
C  Translation on J with rotation on J, off-diagonal
C
      HESS(NTJ1,NRJ1)=SDXAJ+HESS(NTJ1,NRJ1) 
      HESS(NRJ1,NTJ1)=SDXAJ+HESS(NRJ1,NTJ1) 
      HESS(NTJ1,NRJ2)=SDXBJ+HESS(NTJ1,NRJ2) 
      HESS(NRJ2,NTJ1)=SDXBJ+HESS(NRJ2,NTJ1) 
      HESS(NTJ1,NRJ3)=SDXCJ+HESS(NTJ1,NRJ3) 
      HESS(NRJ3,NTJ1)=SDXCJ+HESS(NRJ3,NTJ1) 
      HESS(NTJ2,NRJ1)=SDYAJ+HESS(NTJ2,NRJ1) 
      HESS(NRJ1,NTJ2)=SDYAJ+HESS(NRJ1,NTJ2) 
      HESS(NTJ2,NRJ2)=SDYBJ+HESS(NTJ2,NRJ2) 
      HESS(NRJ2,NTJ2)=SDYBJ+HESS(NRJ2,NTJ2) 
      HESS(NTJ2,NRJ3)=SDYCJ+HESS(NTJ2,NRJ3) 
      HESS(NRJ3,NTJ2)=SDYCJ+HESS(NRJ3,NTJ2) 
      HESS(NTJ3,NRJ1)=SDZAJ+HESS(NTJ3,NRJ1) 
      HESS(NRJ1,NTJ3)=SDZAJ+HESS(NRJ1,NTJ3) 
      HESS(NTJ3,NRJ2)=SDZBJ+HESS(NTJ3,NRJ2) 
      HESS(NRJ2,NTJ3)=SDZBJ+HESS(NRJ2,NTJ3) 
      HESS(NTJ3,NRJ3)=SDZCJ+HESS(NTJ3,NRJ3) 
      HESS(NRJ3,NTJ3)=SDZCJ+HESS(NRJ3,NTJ3) 
C
C  Translation on J with rotation on I, off-diagonal
C
      HESS(NTJ1,NRI1)=-SDXAI+HESS(NTJ1,NRI1) 
      HESS(NRI1,NTJ1)=-SDXAI+HESS(NRI1,NTJ1) 
      HESS(NTJ1,NRI2)=-SDXBI+HESS(NTJ1,NRI2) 
      HESS(NRI2,NTJ1)=-SDXBI+HESS(NRI2,NTJ1) 
      HESS(NTJ1,NRI3)=-SDXCI+HESS(NTJ1,NRI3) 
      HESS(NRI3,NTJ1)=-SDXCI+HESS(NRI3,NTJ1) 
      HESS(NTJ2,NRI1)=-SDYAI+HESS(NTJ2,NRI1) 
      HESS(NRI1,NTJ2)=-SDYAI+HESS(NRI1,NTJ2) 
      HESS(NTJ2,NRI2)=-SDYBI+HESS(NTJ2,NRI2) 
      HESS(NRI2,NTJ2)=-SDYBI+HESS(NRI2,NTJ2) 
      HESS(NTJ2,NRI3)=-SDYCI+HESS(NTJ2,NRI3) 
      HESS(NRI3,NTJ2)=-SDYCI+HESS(NRI3,NTJ2) 
      HESS(NTJ3,NRI1)=-SDZAI+HESS(NTJ3,NRI1) 
      HESS(NRI1,NTJ3)=-SDZAI+HESS(NRI1,NTJ3) 
      HESS(NTJ3,NRI2)=-SDZBI+HESS(NTJ3,NRI2) 
      HESS(NRI2,NTJ3)=-SDZBI+HESS(NRI2,NTJ3) 
      HESS(NTJ3,NRI3)=-SDZCI+HESS(NTJ3,NRI3) 
      HESS(NRI3,NTJ3)=-SDZCI+HESS(NRI3,NTJ3) 
C
C  Translation on I with rotation on J, off-diagonal
C
      HESS(NTI1,NRJ1)=-SDXAJ+HESS(NTI1,NRJ1) 
      HESS(NRJ1,NTI1)=-SDXAJ+HESS(NRJ1,NTI1) 
      HESS(NTI1,NRJ2)=-SDXBJ+HESS(NTI1,NRJ2) 
      HESS(NRJ2,NTI1)=-SDXBJ+HESS(NRJ2,NTI1) 
      HESS(NTI1,NRJ3)=-SDXCJ+HESS(NTI1,NRJ3) 
      HESS(NRJ3,NTI1)=-SDXCJ+HESS(NRJ3,NTI1) 
      HESS(NTI2,NRJ1)=-SDYAJ+HESS(NTI2,NRJ1) 
      HESS(NRJ1,NTI2)=-SDYAJ+HESS(NRJ1,NTI2) 
      HESS(NTI2,NRJ2)=-SDYBJ+HESS(NTI2,NRJ2) 
      HESS(NRJ2,NTI2)=-SDYBJ+HESS(NRJ2,NTI2) 
      HESS(NTI2,NRJ3)=-SDYCJ+HESS(NTI2,NRJ3) 
      HESS(NRJ3,NTI2)=-SDYCJ+HESS(NRJ3,NTI2) 
      HESS(NTI3,NRJ1)=-SDZAJ+HESS(NTI3,NRJ1) 
      HESS(NRJ1,NTI3)=-SDZAJ+HESS(NRJ1,NTI3) 
      HESS(NTI3,NRJ2)=-SDZBJ+HESS(NTI3,NRJ2) 
      HESS(NRJ2,NTI3)=-SDZBJ+HESS(NRJ2,NTI3) 
      HESS(NTI3,NRJ3)=-SDZCJ+HESS(NTI3,NRJ3) 
      HESS(NRJ3,NTI3)=-SDZCJ+HESS(NRJ3,NTI3) 
      ENDIF
   70 CONTINUE 
60    CONTINUE 
      RETURN 
      END 
