C Subroutine to combine GETKD,GETNNZ and GETNINT
C KD is width of sparse band in G matrix
C needed to dimension storage arrays GCS and GCS2
C NNZ is number of non-zero elements in B-matrix
C NINTB is number of internal coordinates including bond lengths (therefore /= NINTS!!!)
C (changed from NINT as this is an intrinsic fortran function!!)

      SUBROUTINE GETSTUFF(KD,NNZ,NINTB)
      USE MODUNRES
      IMPLICIT NONE

      INTEGER KD,NINTB,NNZ
C
C jmc want difference between two atom numbers of dihedral angle.
C in unres, c contains all the calpha coords than all the side chain centroid coords
C main chain dihedrals will always have diff of 3 (4-1,etc) but side chain dihedrals 
C will not be so straightforward.
C SC dihedral is defined as rotation of the side chain centroid about the bisector of the 
C angle defined by Ca(i-1), Ca(i) and Ca(i+1)
C BUT COORDS contains ca, sc, ca, sc etc...
C so max diff will be for main chain dihedral angles: diff = 6 (from 7-1)
 
      KD=3*6+2
C
C Gradient and variable internal coordinate arrays do not contain entries for non-capping glycine
C side chains, so use nside not nres-2.
C Note this is only true if use GLY caps (Ac and NHMe) not NH3+ and COO- (as the latter don't have 
C side chains or even Calphas, whereas the former do.
      NINTB=nphi+ntheta+3*nside+nres-1
C
C 3*number of atoms in each type of int. coor.
C 12 each for alpha and 'omega' (aka beta) - the side chain polar angle and dihedral angle.
      NNZ=6*(nres-1)+6*nside+9*ntheta+12*nphi+24*nside
C again not true if 'D' caps used.
C
      RETURN
      END

c subroutine to calculate B-matrix
C calculation of internal coordinates copied from eintern.src <- charmm only
C NEW now also calculates upper triangle of symmetric G-matrix BEEtBEE on the fly
C stored in sparse format GCS as required for dbptrf routine
C
C Problem is that some of the coords (for glycine side chains) don't actually have 
C internals associated with them in mylbfgs.  Bond lengths have been added to DQ in transback.
C 1. Could have new coords array, excluding the gly side chain coords.
C But, having coords for calpha then side chain for all residues in the chain makes the 
C numbering easier...
C So, coords passed to here still contains zeroes for glycine side chains.
C      SUBROUTINE GETBEE(COORDS,VAL,INDX,JNDX,NINT,NCART,NNZ,GCS,NOCOOR,KD)
      SUBROUTINE UNRSGETBEE(BEEMATRIX,COORDS,VAL,INDX,JNDX,NINTB,NCART,NNZ,GCS,NOCOOR,KD,INTCOORDS)
      USE MODUNRES
      IMPLICIT NONE
C
C NNZ is number of non-zero elements in B-matrix
C NINTB is number of internal coordinates
C NCART is number of cartesians (excluding glycine side chains [including caps])
C
      INTEGER NNZ,NINTB,NCART
C
C jmc      REAL*8 COORDS(6*NATOMS)
      REAL*8 COORDS(NCART)
      REAL*8 X(NCART),Y(NCART),Z(NCART) ! jmc note that NCART=3*NATOMS...
      REAL*8 VAL(NNZ),GCS(2*KD+1,NCART)
      INTEGER INDX(NNZ),JNDX(NNZ)
      REAL*8 RX,RY,RZ,S2,S,SR,RXSR,RYSR,RZSR
      REAL*8 DXI,DYI,DZI,DXJ,DYJ,DZJ,RI2,RJ2,RIR,RJR,RI,RJ
      REAL*8 DXIR,DYIR,DZIR,DXJR,DYJR,DZJR,CST,ISNT
      REAL*8 FX,FY,FZ,GX,GY,GZ,HX,HY,HZ,AX,AY,AZ,BX,BY,BZ
      REAL*8 RF2,RG2,RH2,RF,RG,RH,RFR,RGR,RHR,RA2,RB2,RA2R,RB2R,RABR,CP
      REAL*8 CSTTWO,SNTTWO2,SNTTWO2R,CSTTHREE,SNTTHREE2,SNTTHREE2R,DUMMY,DUMMY2,DUMMY3
      REAL*8 VXI,VYI,VZI,VXJ,VYJ,VZJ,VXK,VYK,VZK,VXL,VYL,VZL
      INTEGER MM,ITH,IPHI,I,J,K,L,IC,I1,SUMPHI,DIFF,ATOMWIDTH
      INTEGER MMM,ITHH,IPHII,II,JJ,KK,LL,IIKD,JJKD,KKKD,LLKD,KD1,KD
      INTEGER II1,JJ1,KK1,LL1,IIKD1,JJKD1,KKKD1,LLKD1
      INTEGER II2,JJ2,KK2,LL2,IIKD2,JJKD2,KKKD2,LLKD2
      LOGICAL NOCOOR
      INTEGER J1,IPASS

      REAL*8 SCALAR,TX,TY,TZ

      REAL*8 BEEMATRIX(NINTB,NCART),INTCOORDS(NINTB)
C
C put cartesians in COORDS in angstroms into X,Y,Z in metres
C
      DO I1=1,nres
        X(2*I1-1)=COORDS(6*(I1-1)+1)*1.0D-10
        Y(2*I1-1)=COORDS(6*(I1-1)+2)*1.0D-10
        Z(2*I1-1)=COORDS(6*(I1-1)+3)*1.0D-10
        X(2*I1)=COORDS(6*(I1-1)+4)*1.0D-10
        Y(2*I1)=COORDS(6*(I1-1)+5)*1.0D-10
        Z(2*I1)=COORDS(6*(I1-1)+6)*1.0D-10
      ENDDO
 
C IPASS is a counter for the number of times the alpha and beta loops are passed (to account for when 
C the loop is skipped due to no side chain angle being present.
      IPASS=0
C
C set ATOMWIDTH to 10 for sparse storage; KD1 is dimension of GCS sparse GC matrix
C that we need to call dbptrf in transform
C in fact using matrix double this size so can calculate upper and lower triangle separately
C will combine the two at the end
C     ATOMWIDTH=6
C     KD=3*ATOMWIDTH+2
      KD1=KD+1
C 
C zero GCS here
C
      DO J1=1,NINTB
        DO I1=1,NCART
           BEEMATRIX(J1,I1)=0.0D0
        END DO
      END DO

       DO J1=1,KD
         DO I1=1,J1
           GCS(KD+1+I1-J1,J1)=0.0D0
           GCS(KD+1+J1-I1,I1)=0.0D0
         ENDDO
       ENDDO
       DO J1=KD+1,NCART
         DO I1=J1-KD,J1
           GCS(KD+1+I1-J1,J1)=0.0D0
           GCS(KD+1+J1-I1,I1)=0.0D0
         ENDDO
       ENDDO
C
C the following loops will put internals in COORDS and pass them back
C if the flag NOCOOR is false
C
C main chain dihedrals
      DO 10 IPHI=1,NPHI
        I=2*IPHI-1
        J=I+2
        K=I+4
        L=I+6
        II=3*(I-1)+1
        II1=3*(I-1)+2
        II2=3*(I-1)+3
        JJ=3*(J-1)+1
        JJ1=3*(J-1)+2
        JJ2=3*(J-1)+3
        KK=3*(K-1)+1
        KK1=3*(K-1)+2
        KK2=3*(K-1)+3
        LL=3*(L-1)+1
        LL1=3*(L-1)+2
        LL2=3*(L-1)+3
C      print *,'pdihe I J K L',I,J,K,L
C       IC=ICP(IPHI)
C F=Ri-Rj, G=Rj-Rk, H-Rl-Rk.
        FX=X(I)-X(J)
        FY=Y(I)-Y(J)
        FZ=Z(I)-Z(J)
        GX=X(J)-X(K)
        GY=Y(J)-Y(K)
        GZ=Z(J)-Z(K)
        HX=X(L)-X(K)
        HY=Y(L)-Y(K)
        HZ=Z(L)-Z(K)
C A=F^G, B=H^G
        AX=FY*GZ-FZ*GY
        AY=FZ*GX-FX*GZ
        AZ=FX*GY-FY*GX
        BX=HY*GZ-HZ*GY
        BY=HZ*GX-HX*GZ
        BZ=HX*GY-HY*GX
C RG=|G|, RGR=1/|G|
        RG2=GX*GX+GY*GY+GZ*GZ
        RG=SQRT(RG2)
        RGR=1/RG
C dae for use in evaluating B-matrix
        RF2=FX*FX+FY*FY+FZ*FZ
        RF=SQRT(RF2)
        RFR=1/RF
        RH2=HX*HX+HY*HY+HZ*HZ
        RH=SQRT(RH2)
        RHR=1/RH
C jmc        CSTTWO=(FX*GX+FY*GY+FZ*GZ)*RFR*RGR
        CSTTWO=-(FX*GX+FY*GY+FZ*GZ)*RFR*RGR
        SNTTWO2=1-CSTTWO*CSTTWO
        SNTTWO2R=1/SNTTWO2
        CSTTHREE=(HX*GX+HY*GY+HZ*GZ)*RHR*RGR
        SNTTHREE2=1-CSTTHREE*CSTTHREE
        SNTTHREE2R=1/SNTTHREE2
        IPHII=12*(IPHI-1)+1
C
        IF (.NOT.NOCOOR) THEN
          RA2=AX*AX+AY*AY+AZ*AZ
          RB2=BX*BX+BY*BY+BZ*BZ
          RA2R=1/RA2
          RB2R=1/RB2
          RABR=SQRT(RA2R*RB2R)
C CP=cos(phi)
          CP=(AX*BX+AY*BY+AZ*BZ)*RABR
          INTCOORDS(IPHI)=ACOS(CP)
          IF (CP.GT.1.0D0) INTCOORDS(IPHI)=0.0D0
          IF (CP.LT.-1.0D0) INTCOORDS(IPHI)=3.141592653589793D0
C jmc 
C Test to determine whether the dihedral angle should be > 0 or < 0, from unres (intcor.f) 
C > 0 if cw rotation from 2-> 1 to 3-> 4.
C ?? Does this affect the B-matrix elements ??
          TX=AY*BZ-BY*AZ
          TY=-AX*BZ+BX*AZ
          TZ=AX*BY-BX*AY
          SCALAR=TX*GX+TY*GY+TZ*GZ
C         PRINT *,'SCALAR ',SCALAR
          IF (SCALAR.GT.0.0D0) INTCOORDS(IPHI)=-INTCOORDS(IPHI)
C         PRINT *,'COORDS = ',INTCOORDS(IPHI)
        ENDIF
C dae calculate B-matrix elements
        DUMMY=RFR*RFR*RGR*SNTTWO2R
        VAL(IPHII)=-AX*DUMMY
        INDX(IPHII)=IPHI
        JNDX(IPHII)=II
        BEEMATRIX(IPHI,II)=-AX*DUMMY
C
        VAL(IPHII+1)=-AY*DUMMY
        INDX(IPHII+1)=IPHI
        JNDX(IPHII+1)=II+1
        BEEMATRIX(IPHI,II1)=-AY*DUMMY
C
        VAL(IPHII+2)=-AZ*DUMMY
        INDX(IPHII+2)=IPHI
        JNDX(IPHII+2)=II+2
        BEEMATRIX(IPHI,II2)=-AZ*DUMMY

        DUMMY=RFR*RFR*RGR*RGR*SNTTWO2R*(RG-RF*CSTTWO)
        DUMMY2=RHR*RGR*RGR*SNTTHREE2R*CSTTHREE
C
        VAL(IPHII+3)=AX*DUMMY-BX*DUMMY2
        INDX(IPHII+3)=IPHI
        JNDX(IPHII+3)=JJ
        BEEMATRIX(IPHI,JJ)=AX*DUMMY-BX*DUMMY2
C
        VAL(IPHII+4)=AY*DUMMY-BY*DUMMY2
        INDX(IPHII+4)=IPHI
        JNDX(IPHII+4)=JJ+1
        BEEMATRIX(IPHI,JJ1)=AY*DUMMY-BY*DUMMY2
C
        VAL(IPHII+5)=AZ*DUMMY-BZ*DUMMY2
        INDX(IPHII+5)=IPHI
        JNDX(IPHII+5)=JJ+2
        BEEMATRIX(IPHI,JJ2)=AZ*DUMMY-BZ*DUMMY2
C
        DUMMY=RHR*RHR*RGR*RGR*SNTTHREE2R*(RG-RH*CSTTHREE)
        DUMMY2=RFR*RGR*RGR*SNTTWO2R*CSTTWO
        VAL(IPHII+6)=-BX*DUMMY+AX*DUMMY2
        INDX(IPHII+6)=IPHI
        JNDX(IPHII+6)=KK
        BEEMATRIX(IPHI,KK)=-BX*DUMMY+AX*DUMMY2
C
        VAL(IPHII+7)=-BY*DUMMY+AY*DUMMY2
        INDX(IPHII+7)=IPHI
        JNDX(IPHII+7)=KK+1
        BEEMATRIX(IPHI,KK1)=-BY*DUMMY+AY*DUMMY2
C
        VAL(IPHII+8)=-BZ*DUMMY+AZ*DUMMY2
        INDX(IPHII+8)=IPHI
        JNDX(IPHII+8)=KK+2
        BEEMATRIX(IPHI,KK2)=-BZ*DUMMY+AZ*DUMMY2
C
        DUMMY=RHR*RHR*RGR*SNTTHREE2R
        VAL(IPHII+9)=BX*DUMMY
        INDX(IPHII+9)=IPHI
        JNDX(IPHII+9)=LL
        BEEMATRIX(IPHI,LL)=BX*DUMMY
C
        VAL(IPHII+10)=BY*DUMMY
        INDX(IPHII+10)=IPHI
        JNDX(IPHII+10)=LL+1
        BEEMATRIX(IPHI,LL1)=BY*DUMMY
C
        VAL(IPHII+11)=BZ*DUMMY
        INDX(IPHII+11)=IPHI
        JNDX(IPHII+11)=LL+2
        BEEMATRIX(IPHI,LL2)=BZ*DUMMY
C
C calculate GEE for this ic
        IIKD=II+KD1
        IIKD1=II1+KD1
        IIKD2=II2+KD1
        JJKD=JJ+KD1
        JJKD1=JJ1+KD1
        JJKD2=JJ2+KD1
        KKKD=KK+KD1
        KKKD1=KK1+KD1
        KKKD2=KK2+KD1
        LLKD=LL+KD1
        LLKD1=LL1+KD1
        LLKD2=LL2+KD1
C
        VXI=VAL(IPHII)
        VYI=VAL(IPHII+1)
        VZI=VAL(IPHII+2)
        GCS(KD1,II)=GCS(KD1,II)+VXI**2
        GCS(KD1+II-II1,II1)=GCS(KD1+II-II1,II1)+VXI*VYI
        GCS(KD1+II-II2,II2)=GCS(KD1+II-II2,II2)+VXI*VZI
        GCS(KD1,II1)=GCS(KD1,II1)+VYI**2
        GCS(KD1+II1-II2,II2)=GCS(KD1+II1-II2,II2)+VYI*VZI
        GCS(KD1,II2)=GCS(KD1,II2)+VZI**2
C
        VXJ=VAL(IPHII+3)
        VYJ=VAL(IPHII+4)
        VZJ=VAL(IPHII+5)
        GCS(KD1,JJ)=GCS(KD1,JJ)+VXJ**2
        GCS(KD1+JJ-JJ1,JJ1)=GCS(KD1+II-II1,JJ1)+VXJ*VYJ
        GCS(KD1+JJ-JJ2,JJ2)=GCS(KD1+II-II2,JJ2)+VXJ*VZJ
        GCS(KD1,JJ1)=GCS(KD1,JJ1)+VYJ**2
        GCS(KD1+JJ1-JJ2,JJ2)=GCS(KD1+JJ1-JJ2,JJ2)+VYJ*VZJ
        GCS(KD1,JJ2)=GCS(KD1,JJ2)+VZJ**2
C
        VXK=VAL(IPHII+6)
        VYK=VAL(IPHII+7)
        VZK=VAL(IPHII+8)
        GCS(KD1,KK)=GCS(KD1,KK)+VXK**2
        GCS(KD1+KK-KK1,KK1)=GCS(KD1+KK-KK1,KK1)+VXK*VYK
        GCS(KD1+KK-KK2,KK2)=GCS(KD1+KK-KK2,KK2)+VXK*VZK
        GCS(KD1,KK1)=GCS(KD1,KK1)+VYK**2
        GCS(KD1+KK1-KK2,KK2)=GCS(KD1+KK1-KK2,KK2)+VYK*VZK
        GCS(KD1,KK2)=GCS(KD1,KK2)+VZK**2
C
        VXL=VAL(IPHII+9)
        VYL=VAL(IPHII+10)
        VZL=VAL(IPHII+11)
        GCS(KD1,LL)=GCS(KD1,LL)+VXL**2
        GCS(KD1+LL-LL1,LL1)=GCS(KD1+LL-LL1,LL1)+VXL*VYL
        GCS(KD1+LL-LL2,LL2)=GCS(KD1+LL-LL2,LL2)+VXL*VZL
        GCS(KD1,LL1)=GCS(KD1,LL1)+VYL**2
        GCS(KD1+LL1-LL2,LL2)=GCS(KD1+LL1-LL2,LL2)+VYL*VZL
        GCS(KD1,LL2)=GCS(KD1,LL2)+VZL**2
C
        IIKD=II+KD1
        IIKD1=II1+KD1
        IIKD2=II2+KD1
        JJKD=JJ+KD1
        JJKD1=JJ1+KD1
        JJKD2=JJ2+KD1
        KKKD=KK+KD1
        KKKD1=KK1+KD1
        KKKD2=KK2+KD1
        LLKD=LL+KD1
        LLKD1=LL1+KD1
        LLKD2=LL2+KD1
C
C       IF (II.LT.JJ) THEN
        GCS(IIKD-JJ,JJ)=GCS(IIKD-JJ,JJ)+VXI*VXJ
        GCS(IIKD-JJ1,JJ1)=GCS(IIKD-JJ1,JJ1)+VXI*VYJ
        GCS(IIKD-JJ2,JJ2)=GCS(IIKD-JJ2,JJ2)+VXI*VZJ
        GCS(IIKD1-JJ,JJ)=GCS(IIKD1-JJ,JJ)+VYI*VXJ
        GCS(IIKD1-JJ1,JJ1)=GCS(IIKD1-JJ1,JJ1)+VYI*VYJ
        GCS(IIKD1-JJ2,JJ2)=GCS(IIKD1-JJ2,JJ2)+VYI*VZJ
        GCS(IIKD2-JJ,JJ)=GCS(IIKD2-JJ,JJ)+VZI*VXJ
        GCS(IIKD2-JJ1,JJ1)=GCS(IIKD2-JJ1,JJ1)+VZI*VYJ
        GCS(IIKD2-JJ2,JJ2)=GCS(IIKD2-JJ2,JJ2)+VZI*VZJ
C       ELSE
C       GCS(JJKD-II,II)=GCS(IIKD-JJ,JJ)
C       GCS(JJKD1-II,II)=GCS(IIKD-JJ1,JJ1)
C       GCS(JJKD2-II,II)=GCS(IIKD-JJ2,JJ2)
C       GCS(JJKD-II1,II1)=GCS(IIKD1-JJ,JJ)
C       GCS(JJKD1-II1,II1)=GCS(IIKD1-JJ1,JJ1)
C       GCS(JJKD2-II1,II1)=GCS(IIKD1-JJ2,JJ2)
C       GCS(JJKD-II2,II2)=GCS(IIKD2-JJ,JJ)
C       GCS(JJKD1-II2,II2)=GCS(IIKD2-JJ1,JJ1)
C       GCS(JJKD2-II2,II2)=GCS(IIKD2-JJ2,JJ2)
C       ENDIF
C
C       IF (II.LT.KK) THEN
        GCS(IIKD-KK,KK)=GCS(IIKD-KK,KK)+VXI*VXK
        GCS(IIKD-KK1,KK1)=GCS(IIKD-KK1,KK1)+VXI*VYK
        GCS(IIKD-KK2,KK2)=GCS(IIKD-KK2,KK2)+VXI*VZK
        GCS(IIKD1-KK,KK)=GCS(IIKD1-KK,KK)+VYI*VXK
        GCS(IIKD1-KK1,KK1)=GCS(IIKD1-KK1,KK1)+VYI*VYK
        GCS(IIKD1-KK2,KK2)=GCS(IIKD1-KK2,KK2)+VYI*VZK
        GCS(IIKD2-KK,KK)=GCS(IIKD2-KK,KK)+VZI*VXK
        GCS(IIKD2-KK1,KK1)=GCS(IIKD2-KK1,KK1)+VZI*VYK
        GCS(IIKD2-KK2,KK2)=GCS(IIKD2-KK2,KK2)+VZI*VZK
C       ELSE
C       GCS(KKKD-II,II)=GCS(IIKD-KK,KK)
C       GCS(KKKD1-II,II)=GCS(IIKD-KK1,KK1)
C       GCS(KKKD2-II,II)=GCS(IIKD-KK2,KK2)
C       GCS(KKKD-II1,II1)=GCS(IIKD1-KK,KK)
C       GCS(KKKD1-II1,II1)=GCS(IIKD1-KK1,KK1)
C       GCS(KKKD2-II1,II1)=GCS(IIKD1-KK2,KK2)
C       GCS(KKKD-II2,II2)=GCS(IIKD2-KK,KK)
C       GCS(KKKD1-II2,II2)=GCS(IIKD2-KK1,KK1)
C       GCS(KKKD2-II2,II2)=GCS(IIKD2-KK2,KK2)
C       ENDIF
C
C       IF (II.LT.LL) THEN
        GCS(IIKD-LL,LL)=GCS(IIKD-LL,LL)+VXI*VXL
        GCS(IIKD-LL1,LL1)=GCS(IIKD-LL1,LL1)+VXI*VYL
        GCS(IIKD-LL2,LL2)=GCS(IIKD-LL2,LL2)+VXI*VZL
        GCS(IIKD1-LL,LL)=GCS(IIKD1-LL,LL)+VYI*VXL
        GCS(IIKD1-LL1,LL1)=GCS(IIKD1-LL1,LL1)+VYI*VYL
        GCS(IIKD1-LL2,LL2)=GCS(IIKD1-LL2,LL2)+VYI*VZL
        GCS(IIKD2-LL,LL)=GCS(IIKD2-LL,LL)+VZI*VXL
        GCS(IIKD2-LL1,LL1)=GCS(IIKD2-LL1,LL1)+VZI*VYL
        GCS(IIKD2-LL2,LL2)=GCS(IIKD2-LL2,LL2)+VZI*VZL
C       ELSE
C       GCS(LLKD-II,II)=GCS(IIKD-LL,LL)
C       GCS(LLKD1-II,II)=GCS(IIKD-LL1,LL1)
C       GCS(LLKD2-II,II)=GCS(IIKD-LL2,LL2)
C       GCS(LLKD-II1,II1)=GCS(IIKD1-LL,LL)
C       GCS(LLKD1-II1,II1)=GCS(IIKD1-LL1,LL1)
C       GCS(LLKD2-II1,II1)=GCS(IIKD1-LL2,LL2)
C       GCS(LLKD-II2,II2)=GCS(IIKD2-LL,LL)
C       GCS(LLKD1-II2,II2)=GCS(IIKD2-LL1,LL1)
C       GCS(LLKD2-II2,II2)=GCS(IIKD2-LL2,LL2)
C       ENDIF
C
C       IF (JJ.LT.KK) THEN
        GCS(JJKD-KK,KK)=GCS(JJKD-KK,KK)+VXJ*VXK
        GCS(JJKD-KK1,KK1)=GCS(JJKD-KK1,KK1)+VXJ*VYK
        GCS(JJKD-KK2,KK2)=GCS(JJKD-KK2,KK2)+VXJ*VZK
        GCS(JJKD1-KK,KK)=GCS(JJKD1-KK,KK)+VYJ*VXK
        GCS(JJKD1-KK1,KK1)=GCS(JJKD1-KK1,KK1)+VYJ*VYK
        GCS(JJKD1-KK2,KK2)=GCS(JJKD1-KK2,KK2)+VYJ*VZK
        GCS(JJKD2-KK,KK)=GCS(JJKD2-KK,KK)+VZJ*VXK
        GCS(JJKD2-KK1,KK1)=GCS(JJKD2-KK1,KK1)+VZJ*VYK
        GCS(JJKD2-KK2,KK2)=GCS(JJKD2-KK2,KK2)+VZJ*VZK
C       ELSE
C       GCS(KKKD-JJ,JJ)=GCS(JJKD-KK,KK)+1000
C       GCS(KKKD1-JJ,JJ)=GCS(JJKD-KK1,KK1)+2000
C       GCS(KKKD2-JJ,JJ)=GCS(JJKD-KK2,KK2)+3000
C       GCS(KKKD-JJ1,JJ1)=GCS(JJKD1-KK,KK)+4000
C       GCS(KKKD1-JJ1,JJ1)=GCS(JJKD1-KK1,KK1)+5000
C       GCS(KKKD2-JJ1,JJ1)=GCS(JJKD1-KK2,KK2)+6000
C       GCS(KKKD-JJ2,JJ2)=GCS(JJKD2-KK,KK)
C       GCS(KKKD1-JJ2,JJ2)=GCS(JJKD2-KK1,KK1)+8000
C       GCS(KKKD2-JJ2,JJ2)=GCS(JJKD2-KK2,KK2)+9000
C       ENDIF
C
C       IF (JJ.LT.LL) THEN
        GCS(JJKD-LL,LL)=GCS(JJKD-LL,LL)+VXJ*VXL
        GCS(JJKD-LL1,LL1)=GCS(JJKD-LL1,LL1)+VXJ*VYL
        GCS(JJKD-LL2,LL2)=GCS(JJKD-LL2,LL2)+VXJ*VZL
        GCS(JJKD1-LL,LL)=GCS(JJKD1-LL,LL)+VYJ*VXL
        GCS(JJKD1-LL1,LL1)=GCS(JJKD1-LL1,LL1)+VYJ*VYL
        GCS(JJKD1-LL2,LL2)=GCS(JJKD1-LL2,LL2)+VYJ*VZL
        GCS(JJKD2-LL,LL)=GCS(JJKD2-LL,LL)+VZJ*VXL
        GCS(JJKD2-LL1,LL1)=GCS(JJKD2-LL1,LL1)+VZJ*VYL
        GCS(JJKD2-LL2,LL2)=GCS(JJKD2-LL2,LL2)+VZJ*VZL
C       ELSE
C       GCS(LLKD-JJ,JJ)=GCS(JJKD-LL,LL)
C       GCS(LLKD1-JJ,JJ)=GCS(JJKD-LL1,LL1)
C       GCS(LLKD2-JJ,JJ)=GCS(JJKD-LL2,LL2)
C       GCS(LLKD-JJ1,JJ1)=GCS(JJKD1-LL,LL)
C       GCS(LLKD1-JJ1,JJ1)=GCS(JJKD1-LL1,LL1)
C       GCS(LLKD2-JJ1,JJ1)=GCS(JJKD1-LL2,LL2)
C       GCS(LLKD-JJ2,JJ2)=GCS(JJKD2-LL,LL)
C       GCS(LLKD1-JJ2,JJ2)=GCS(JJKD2-LL1,LL1)
C       GCS(LLKD2-JJ2,JJ2)=GCS(JJKD2-LL2,LL2)
C       ENDIF
C
C       IF (KK.LT.LL) THEN
        GCS(KKKD-LL,LL)=GCS(KKKD-LL,LL)+VXK*VXL
        GCS(KKKD-LL1,LL1)=GCS(KKKD-LL1,LL1)+VXK*VYL
        GCS(KKKD-LL2,LL2)=GCS(KKKD-LL2,LL2)+VXK*VZL
        GCS(KKKD1-LL,LL)=GCS(KKKD1-LL,LL)+VYK*VXL
        GCS(KKKD1-LL1,LL1)=GCS(KKKD1-LL1,LL1)+VYK*VYL
        GCS(KKKD1-LL2,LL2)=GCS(KKKD1-LL2,LL2)+VYK*VZL
        GCS(KKKD2-LL,LL)=GCS(KKKD2-LL,LL)+VZK*VXL
        GCS(KKKD2-LL1,LL1)=GCS(KKKD2-LL1,LL1)+VZK*VYL
        GCS(KKKD2-LL2,LL2)=GCS(KKKD2-LL2,LL2)+VZK*VZL
C       ELSE
C       GCS(LLKD-KK,KK)=GCS(KKKD-LL,LL)
C       GCS(LLKD1-KK,KK)=GCS(KKKD-LL1,LL1)
C       GCS(LLKD2-KK,KK)=GCS(KKKD-LL2,LL2)
C       GCS(LLKD-KK1,KK1)=GCS(KKKD1-LL,LL)
C       GCS(LLKD1-KK1,KK1)=GCS(KKKD1-LL1,LL1)
C       GCS(LLKD2-KK1,KK1)=GCS(KKKD1-LL2,LL2)
C       GCS(LLKD-KK2,KK2)=GCS(KKKD2-LL,LL)
C       GCS(LLKD1-KK2,KK2)=GCS(KKKD2-LL1,LL1)
C       GCS(LLKD2-KK2,KK2)=GCS(KKKD2-LL2,LL2)
C       ENDIF
C
C end dae
C
   10 CONTINUE
C      DO I1=1,NPHI
C        print *,'phi I J K L',IP(I1),JP(I1),KP(I1),LP(I1)
C      ENDDO

      DO 20 ITH=1,NTHETA
C
C theta i is angle between Cai, Cai+1 and Cai+2
        I=2*ITH-1
        J=I+2
        K=I+4
        II=3*(I-1)+1
        II1=3*(I-1)+2
        II2=3*(I-1)+3
        JJ=3*(J-1)+1
        JJ1=3*(J-1)+2
        JJ2=3*(J-1)+3
        KK=3*(K-1)+1
        KK1=3*(K-1)+2
        KK2=3*(K-1)+3
C      print *,'angle I J K',I,J,K
C       IC=ICT(ITH)
C       IF(IC.EQ.0) GOTO 20
        DXI=X(I)-X(J)
        DYI=Y(I)-Y(J)
        DZI=Z(I)-Z(J)
        DXJ=X(K)-X(J)
        DYJ=Y(K)-Y(J)
        DZJ=Z(K)-Z(J)
        RI2=DXI*DXI+DYI*DYI+DZI*DZI
        RJ2=DXJ*DXJ+DYJ*DYJ+DZJ*DZJ
        RI=SQRT(RI2)
        RJ=SQRT(RJ2)
        RIR=1/RI
        RJR=1/RJ
        DXIR=DXI*RIR
        DYIR=DYI*RIR
        DZIR=DZI*RIR
        DXJR=DXJ*RJR
        DYJR=DYJ*RJR
        DZJR=DZJ*RJR
        CST=DXIR*DXJR+DYIR*DYJR+DZIR*DZJR
C
        IF (.NOT.NOCOOR) THEN
C dae store theta (in radians) in COORDS
          INTCOORDS(nphi+ITH)=ACOS(CST)
c         PRINT *,'COORDS = ',COORDS(nphi+ITH)
        ENDIF
C dae calculate B-matrix elements for these atoms (Wilson Decius + Cross ch.4)
        ISNT=1/SQRT(1-CST*CST)
        ITHH=12*nphi+9*(ITH-1)+1
        VAL(ITHH)=(DXIR*CST-DXJR)*RIR*ISNT
        INDX(ITHH)=nphi+ITH
        JNDX(ITHH)=II
        BEEMATRIX(nphi+ITH,II)=(DXIR*CST-DXJR)*RIR*ISNT
C
        VAL(ITHH+1)=(DYIR*CST-DYJR)*RIR*ISNT
        INDX(ITHH+1)=nphi+ITH
        JNDX(ITHH+1)=II1
        BEEMATRIX(nphi+ITH,II1)=(DYIR*CST-DYJR)*RIR*ISNT
C
        VAL(ITHH+2)=(DZIR*CST-DZJR)*RIR*ISNT
        INDX(ITHH+2)=nphi+ITH
        JNDX(ITHH+2)=II2
        BEEMATRIX(nphi+ITH,II2)=(DZIR*CST-DZJR)*RIR*ISNT
C
        VAL(ITHH+3)=(DXIR*(RI-RJ*CST)+DXJR*(RJ-RI*CST))*RIR*RJR*ISNT
        INDX(ITHH+3)=nphi+ITH
        JNDX(ITHH+3)=JJ
        BEEMATRIX(nphi+ITH,JJ)=(DXIR*(RI-RJ*CST)+DXJR*(RJ-RI*CST))*RIR*RJR*ISNT
C
        VAL(ITHH+4)=(DYIR*(RI-RJ*CST)+DYJR*(RJ-RI*CST))*RIR*RJR*ISNT
        INDX(ITHH+4)=nphi+ITH
        JNDX(ITHH+4)=JJ1
        BEEMATRIX(nphi+ITH,JJ1)=(DYIR*(RI-RJ*CST)+DYJR*(RJ-RI*CST))*RIR*RJR*ISNT
C
        VAL(ITHH+5)=(DZIR*(RI-RJ*CST)+DZJR*(RJ-RI*CST))*RIR*RJR*ISNT
        INDX(ITHH+5)=nphi+ITH
        JNDX(ITHH+5)=JJ2
        BEEMATRIX(nphi+ITH,JJ2)=(DZIR*(RI-RJ*CST)+DZJR*(RJ-RI*CST))*RIR*RJR*ISNT
C
        VAL(ITHH+6)=(DXJR*CST-DXIR)*RJR*ISNT
        INDX(ITHH+6)=nphi+ITH
        JNDX(ITHH+6)=KK
        BEEMATRIX(nphi+ITH,KK)=(DXJR*CST-DXIR)*RJR*ISNT
C
        VAL(ITHH+7)=(DYJR*CST-DYIR)*RJR*ISNT
        INDX(ITHH+7)=nphi+ITH
        JNDX(ITHH+7)=KK1
        BEEMATRIX(nphi+ITH,KK1)=(DYJR*CST-DYIR)*RJR*ISNT
C
        VAL(ITHH+8)=(DZJR*CST-DZIR)*RJR*ISNT
        INDX(ITHH+8)=nphi+ITH
        JNDX(ITHH+8)=KK2
        BEEMATRIX(nphi+ITH,KK2)=(DZJR*CST-DZIR)*RJR*ISNT
C
C calculate GEE for this ic
        VXI=VAL(ITHH)
        VYI=VAL(ITHH+1)
        VZI=VAL(ITHH+2)
        GCS(KD1,II)=GCS(KD1,II)+VXI**2
        GCS(KD1+II-II1,II1)=GCS(KD1+II-II1,II1)+VXI*VYI
        GCS(KD1+II-II2,II2)=GCS(KD1+II-II2,II2)+VXI*VZI
        GCS(KD1,II1)=GCS(KD1,II1)+VYI**2
        GCS(KD1+II1-II2,II2)=GCS(KD1+II1-II2,II2)+VYI*VZI
        GCS(KD1,II2)=GCS(KD1,II2)+VZI**2
C
        VXJ=VAL(ITHH+3)
        VYJ=VAL(ITHH+4)
        VZJ=VAL(ITHH+5)
        GCS(KD1,JJ)=GCS(KD1,JJ)+VXJ**2
        GCS(KD1+JJ-JJ1,JJ1)=GCS(KD1+JJ-JJ1,JJ1)+VXJ*VYJ
        GCS(KD1+JJ-JJ2,JJ2)=GCS(KD1+JJ-JJ2,JJ2)+VXJ*VZJ
        GCS(KD1,JJ1)=GCS(KD1,JJ1)+VYJ**2
        GCS(KD1+JJ1-JJ2,JJ2)=GCS(KD1+JJ1-JJ2,JJ2)+VYJ*VZJ
        GCS(KD1,JJ2)=GCS(KD1,JJ2)+VZJ**2
C
        VXK=VAL(ITHH+6)
        VYK=VAL(ITHH+7)
        VZK=VAL(ITHH+8)
        GCS(KD1,KK)=GCS(KD1,KK)+VXK**2
        GCS(KD1+KK-KK1,KK1)=GCS(KD1+KK-KK1,KK1)+VXK*VYK
        GCS(KD1+KK-KK2,KK2)=GCS(KD1+KK-KK2,KK2)+VXK*VZK
        GCS(KD1,KK1)=GCS(KD1,KK1)+VYK**2
        GCS(KD1+KK1-KK2,KK2)=GCS(KD1+KK1-KK2,KK2)+VYK*VZK
        GCS(KD1,KK2)=GCS(KD1,KK2)+VZK**2
C
        IIKD=II+KD1
        IIKD1=II1+KD1
        IIKD2=II2+KD1
        JJKD=JJ+KD1
        JJKD1=JJ1+KD1
        JJKD2=JJ2+KD1
        KKKD=KK+KD1
        KKKD1=KK1+KD1
        KKKD2=KK2+KD1
C
C       IF (II.LT.JJ) THEN
        GCS(IIKD-JJ,JJ)=GCS(IIKD-JJ,JJ)+VXI*VXJ
        GCS(IIKD-JJ1,JJ1)=GCS(IIKD-JJ1,JJ1)+VXI*VYJ
        GCS(IIKD-JJ2,JJ2)=GCS(IIKD-JJ2,JJ2)+VXI*VZJ
        GCS(IIKD1-JJ,JJ)=GCS(IIKD1-JJ,JJ)+VYI*VXJ
        GCS(IIKD1-JJ1,JJ1)=GCS(IIKD1-JJ1,JJ1)+VYI*VYJ
        GCS(IIKD1-JJ2,JJ2)=GCS(IIKD1-JJ2,JJ2)+VYI*VZJ
        GCS(IIKD2-JJ,JJ)=GCS(IIKD2-JJ,JJ)+VZI*VXJ
        GCS(IIKD2-JJ1,JJ1)=GCS(IIKD2-JJ1,JJ1)+VZI*VYJ
        GCS(IIKD2-JJ2,JJ2)=GCS(IIKD2-JJ2,JJ2)+VZI*VZJ
C       ELSE
C       GCS(JJKD-II,II)=GCS(IIKD-JJ,JJ)
C       GCS(JJKD1-II,II)=GCS(IIKD-JJ1,JJ1)
C       GCS(JJKD2-II,II)=GCS(IIKD-JJ2,JJ2)
C       GCS(JJKD-II1,II1)=GCS(IIKD1-JJ,JJ)
C       GCS(JJKD1-II1,II1)=GCS(IIKD1-JJ1,JJ1)
C       GCS(JJKD2-II1,II1)=GCS(IIKD1-JJ2,JJ2)
C       GCS(JJKD-II2,II2)=GCS(IIKD2-JJ,JJ)
C       GCS(JJKD1-II2,II2)=GCS(IIKD2-JJ1,JJ1)
C       GCS(JJKD2-II2,II2)=GCS(IIKD2-JJ2,JJ2)
C       ENDIF
C
C       IF (II.LT.KK) THEN
        GCS(IIKD-KK,KK)=GCS(IIKD-KK,KK)+VXI*VXK
        GCS(IIKD-KK1,KK1)=GCS(IIKD-KK1,KK1)+VXI*VYK
        GCS(IIKD-KK2,KK2)=GCS(IIKD-KK2,KK2)+VXI*VZK
        GCS(IIKD1-KK,KK)=GCS(IIKD1-KK,KK)+VYI*VXK
        GCS(IIKD1-KK1,KK1)=GCS(IIKD1-KK1,KK1)+VYI*VYK
        GCS(IIKD1-KK2,KK2)=GCS(IIKD1-KK2,KK2)+VYI*VZK
        GCS(IIKD2-KK,KK)=GCS(IIKD2-KK,KK)+VZI*VXK
        GCS(IIKD2-KK1,KK1)=GCS(IIKD2-KK1,KK1)+VZI*VYK
        GCS(IIKD2-KK2,KK2)=GCS(IIKD2-KK2,KK2)+VZI*VZK
C       ELSE
C       GCS(KKKD-II,II)=GCS(IIKD-KK,KK)
C       GCS(KKKD1-II,II)=GCS(IIKD-KK1,KK1)
C       GCS(KKKD2-II,II)=GCS(IIKD-KK2,KK2)
C       GCS(KKKD-II1,II1)=GCS(IIKD1-KK,KK)
C       GCS(KKKD1-II1,II1)=GCS(IIKD1-KK1,KK1)
C       GCS(KKKD2-II1,II1)=GCS(IIKD1-KK2,KK2)
C       GCS(KKKD-II2,II2)=GCS(IIKD2-KK,KK)
C       GCS(KKKD1-II2,II2)=GCS(IIKD2-KK1,KK1)
C       GCS(KKKD2-II2,II2)=GCS(IIKD2-KK2,KK2)
C       ENDIF
C
C       IF (JJ.LT.KK) THEN
        GCS(JJKD-KK,KK)=GCS(JJKD-KK,KK)+VXJ*VXK
        GCS(JJKD-KK1,KK1)=GCS(JJKD-KK1,KK1)+VXJ*VYK
        GCS(JJKD-KK2,KK2)=GCS(JJKD-KK2,KK2)+VXJ*VZK
        GCS(JJKD1-KK,KK)=GCS(JJKD1-KK,KK)+VYJ*VXK
        GCS(JJKD1-KK1,KK1)=GCS(JJKD1-KK1,KK1)+VYJ*VYK
        GCS(JJKD1-KK2,KK2)=GCS(JJKD1-KK2,KK2)+VYJ*VZK
        GCS(JJKD2-KK,KK)=GCS(JJKD2-KK,KK)+VZJ*VXK
        GCS(JJKD2-KK1,KK1)=GCS(JJKD2-KK1,KK1)+VZJ*VYK
        GCS(JJKD2-KK2,KK2)=GCS(JJKD2-KK2,KK2)+VZJ*VZK
C       ELSE
C       GCS(KKKD-JJ,JJ)=GCS(JJKD-KK,KK)
C       GCS(KKKD1-JJ,JJ)=GCS(JJKD-KK1,KK1)
C       GCS(KKKD2-JJ,JJ)=GCS(JJKD-KK2,KK2)
C       GCS(KKKD-JJ1,JJ1)=GCS(JJKD1-KK,KK)
C       GCS(KKKD1-JJ1,JJ1)=GCS(JJKD1-KK1,KK1)
C       GCS(KKKD2-JJ1,JJ1)=GCS(JJKD1-KK2,KK2)
C       GCS(KKKD-JJ2,JJ2)=GCS(JJKD2-KK,KK)
C       GCS(KKKD1-JJ2,JJ2)=GCS(JJKD2-KK1,KK1)
C       GCS(KKKD2-JJ2,JJ2)=GCS(JJKD2-KK2,KK2)
C       ENDIF
C
C
C end dae
C
   20  CONTINUE

       IPASS=0
C side chain polar angle alpha
C alpha is angle between the vector SCi-Cai, and the bisector of theta(i-1) in the plane of the 
C three Ca's. Depends therefore on the coordinates of 4 atoms
C ***********************************************************************************************
C NOTE that glycine residues (whether capping or normal) DO NOT have alphas or dihedral angles...
C therefore can't simply use nside to work out I,J,K,L...
C ***********************************************************************************************
      DO 30 ITH=1,nres-2
C
C skip this pass if the 'itype' of the Ca(i) is glycine, 
        IF (ITYPE(ITH+1).EQ.10) GOTO 30
        IPASS=IPASS+1
        I=2*ITH-1
        J=I+2
C K is the side chain centroid.
        K=I+3
        L=I+4
        II=3*(I-1)+1
        II1=3*(I-1)+2
        II2=3*(I-1)+3
        JJ=3*(J-1)+1
        JJ1=3*(J-1)+2
        JJ2=3*(J-1)+3
        KK=3*(K-1)+1
        KK1=3*(K-1)+2
        KK2=3*(K-1)+3
        LL=3*(L-1)+1
        LL1=3*(L-1)+2
        LL2=3*(L-1)+3
C      print *,'angle I J K L',I,J,K,L
C       IC=ICT(ITH)
C       IF(IC.EQ.0) GOTO 30
C DXI is vector Cai -> SCi
        DXI=X(K)-X(J)
        DYI=Y(K)-Y(J)
        DZI=Z(K)-Z(J)
C DXJ is bisector of I->L
        DXJ=X(J)-(X(I)+X(L))/2.0D0
        DYJ=Y(J)-(Y(I)+Y(L))/2.0D0
        DZJ=Z(J)-(Z(I)+Z(L))/2.0D0
        RI2=DXI*DXI+DYI*DYI+DZI*DZI
        RJ2=DXJ*DXJ+DYJ*DYJ+DZJ*DZJ
        RI=SQRT(RI2)
        RJ=SQRT(RJ2)
        RIR=1/RI
        RJR=1/RJ
        DXIR=DXI*RIR
        DYIR=DYI*RIR
        DZIR=DZI*RIR
        DXJR=DXJ*RJR
        DYJR=DYJ*RJR
        DZJR=DZJ*RJR
        CST=DXIR*DXJR+DYIR*DYJR+DZIR*DZJR
C
        IF (.NOT.NOCOOR) THEN
C dae store theta (in radians) in COORDS
C         COORDS(nphi+ntheta+ITH)=3.141592653589793D0-ACOS(CST)
          INTCOORDS(nphi+ntheta+IPASS)=3.141592653589793D0-ACOS(CST)
C         PRINT *,'COORDS = ',COORDS(nphi+ntheta+IPASS)
        ENDIF
C dae calculate B-matrix elements for these atoms (Wilson Decius + Cross ch.4)
C jmc testing!!!        ISNT=1/SQRT(1-CST*CST)
        ISNT=-1/SQRT(1-CST*CST)
        ITHH=12*nphi+9*ntheta+12*(IPASS-1)+1
        VAL(ITHH)=-(DXJR*CST-DXIR)*RJR*ISNT/2.0D0
        INDX(ITHH)=nphi+ntheta+IPASS
        JNDX(ITHH)=II
        BEEMATRIX(nphi+ntheta+IPASS,II)=-(DXJR*CST-DXIR)*RJR*ISNT/2.0D0
C
        VAL(ITHH+1)=-(DYJR*CST-DYIR)*RJR*ISNT/2.0D0
        INDX(ITHH+1)=nphi+ntheta+IPASS
        JNDX(ITHH+1)=II1
        BEEMATRIX(nphi+ntheta+IPASS,II1)=-(DYJR*CST-DYIR)*RJR*ISNT/2.0D0
C
        VAL(ITHH+2)=-(DZJR*CST-DZIR)*RJR*ISNT/2.0D0
        INDX(ITHH+2)=nphi+ntheta+IPASS
        JNDX(ITHH+2)=II2
        BEEMATRIX(nphi+ntheta+IPASS,II2)=-(DZJR*CST-DZIR)*RJR*ISNT/2.0D0
C
        VAL(ITHH+3)=(DXIR*(RI-RJ*CST)+DXJR*(RJ-RI*CST))*RIR*RJR*ISNT
     &                                   +2.0D0*(DXJR*CST-DXIR)*RJR*ISNT
        INDX(ITHH+3)=nphi+ntheta+IPASS
        JNDX(ITHH+3)=JJ
        BEEMATRIX(nphi+ntheta+IPASS,JJ)=(DXIR*(RI-RJ*CST)+DXJR*(RJ-RI*CST))*RIR*RJR*ISNT
     &                                   +2.0D0*(DXJR*CST-DXIR)*RJR*ISNT
C
        VAL(ITHH+4)=(DYIR*(RI-RJ*CST)+DYJR*(RJ-RI*CST))*RIR*RJR*ISNT
     &                                   +2.0D0*(DYJR*CST-DYIR)*RJR*ISNT
        INDX(ITHH+4)=nphi+ntheta+IPASS
        JNDX(ITHH+4)=JJ1
        BEEMATRIX(nphi+ntheta+IPASS,JJ1)=(DYIR*(RI-RJ*CST)+DYJR*(RJ-RI*CST))*RIR*RJR*ISNT
     &                                   +2.0D0*(DYJR*CST-DYIR)*RJR*ISNT
C
        VAL(ITHH+5)=(DZIR*(RI-RJ*CST)+DZJR*(RJ-RI*CST))*RIR*RJR*ISNT
     &                                   +2.0D0*(DZJR*CST-DZIR)*RJR*ISNT
        INDX(ITHH+5)=nphi+ntheta+IPASS
        JNDX(ITHH+5)=JJ2
        BEEMATRIX(nphi+ntheta+IPASS,JJ2)=(DZIR*(RI-RJ*CST)+DZJR*(RJ-RI*CST))*RIR*RJR*ISNT
     &                                   +2.0D0*(DZJR*CST-DZIR)*RJR*ISNT
C
        VAL(ITHH+6)=(DXIR*CST-DXJR)*RIR*ISNT
        INDX(ITHH+6)=nphi+ntheta+IPASS
        JNDX(ITHH+6)=KK
        BEEMATRIX(nphi+ntheta+IPASS,KK)=(DXIR*CST-DXJR)*RIR*ISNT
C
        VAL(ITHH+7)=(DYIR*CST-DYJR)*RIR*ISNT
        INDX(ITHH+7)=nphi+ntheta+IPASS
        JNDX(ITHH+7)=KK1
        BEEMATRIX(nphi+ntheta+IPASS,KK1)=(DYIR*CST-DYJR)*RIR*ISNT
C
        VAL(ITHH+8)=(DZIR*CST-DZJR)*RIR*ISNT
        INDX(ITHH+8)=nphi+ntheta+IPASS
        JNDX(ITHH+8)=KK2
        BEEMATRIX(nphi+ntheta+IPASS,KK2)=(DZIR*CST-DZJR)*RIR*ISNT
C
        VAL(ITHH+9)=-(DXJR*CST-DXIR)*RJR*ISNT/2.0D0
        INDX(ITHH+9)=nphi+ntheta+IPASS
        JNDX(ITHH+9)=LL
        BEEMATRIX(nphi+ntheta+IPASS,LL)=-(DXJR*CST-DXIR)*RJR*ISNT/2.0D0
C
        VAL(ITHH+10)=-(DYJR*CST-DYIR)*RJR*ISNT/2.0D0
        INDX(ITHH+10)=nphi+ntheta+IPASS
        JNDX(ITHH+10)=LL1
        BEEMATRIX(nphi+ntheta+IPASS,LL1)=-(DYJR*CST-DYIR)*RJR*ISNT/2.0D0
C
        VAL(ITHH+11)=-(DZJR*CST-DZIR)*RJR*ISNT/2.0D0
        INDX(ITHH+11)=nphi+ntheta+IPASS
        JNDX(ITHH+11)=LL2
        BEEMATRIX(nphi+ntheta+IPASS,LL2)=-(DZJR*CST-DZIR)*RJR*ISNT/2.0D0
         
C
C calculate GEE for this ic
        VXI=VAL(ITHH)
        VYI=VAL(ITHH+1)
        VZI=VAL(ITHH+2)
        GCS(KD1,II)=GCS(KD1,II)+VXI**2
        GCS(KD1+II-II1,II1)=GCS(KD1+II-II1,II1)+VXI*VYI
        GCS(KD1+II-II2,II2)=GCS(KD1+II-II2,II2)+VXI*VZI
        GCS(KD1,II1)=GCS(KD1,II1)+VYI**2
        GCS(KD1+II1-II2,II2)=GCS(KD1+II1-II2,II2)+VYI*VZI
        GCS(KD1,II2)=GCS(KD1,II2)+VZI**2
C
        VXJ=VAL(ITHH+3)
        VYJ=VAL(ITHH+4)
        VZJ=VAL(ITHH+5)
        GCS(KD1,JJ)=GCS(KD1,JJ)+VXJ**2
        GCS(KD1+JJ-JJ1,JJ1)=GCS(KD1+JJ-JJ1,JJ1)+VXJ*VYJ
        GCS(KD1+JJ-JJ2,JJ2)=GCS(KD1+JJ-JJ2,JJ2)+VXJ*VZJ
        GCS(KD1,JJ1)=GCS(KD1,JJ1)+VYJ**2
        GCS(KD1+JJ1-JJ2,JJ2)=GCS(KD1+JJ1-JJ2,JJ2)+VYJ*VZJ
        GCS(KD1,JJ2)=GCS(KD1,JJ2)+VZJ**2
C
        VXK=VAL(ITHH+6)
        VYK=VAL(ITHH+7)
        VZK=VAL(ITHH+8)
        GCS(KD1,KK)=GCS(KD1,KK)+VXK**2
        GCS(KD1+KK-KK1,KK1)=GCS(KD1+KK-KK1,KK1)+VXK*VYK
        GCS(KD1+KK-KK2,KK2)=GCS(KD1+KK-KK2,KK2)+VXK*VZK
        GCS(KD1,KK1)=GCS(KD1,KK1)+VYK**2
        GCS(KD1+KK1-KK2,KK2)=GCS(KD1+KK1-KK2,KK2)+VYK*VZK
        GCS(KD1,KK2)=GCS(KD1,KK2)+VZK**2
C
        VXL=VAL(ITHH+9)
        VYL=VAL(ITHH+10)
        VZL=VAL(ITHH+11)
        GCS(KD1,LL)=GCS(KD1,LL)+VXL**2
        GCS(KD1+LL-LL1,LL1)=GCS(KD1+LL-LL1,LL1)+VXL*VYL
        GCS(KD1+LL-LL2,LL2)=GCS(KD1+LL-LL2,LL2)+VXL*VZL
        GCS(KD1,LL1)=GCS(KD1,LL1)+VYL**2
        GCS(KD1+LL1-LL2,LL2)=GCS(KD1+LL1-LL2,LL2)+VYL*VZL
        GCS(KD1,LL2)=GCS(KD1,LL2)+VZL**2
C
        IIKD=II+KD1
        IIKD1=II1+KD1
        IIKD2=II2+KD1
        JJKD=JJ+KD1
        JJKD1=JJ1+KD1
        JJKD2=JJ2+KD1
        KKKD=KK+KD1
        KKKD1=KK1+KD1
        KKKD2=KK2+KD1
        LLKD=LL+KD1
        LLKD1=LL1+KD1
        LLKD2=LL2+KD1
C
C       IF (II.LT.JJ) THEN
        GCS(IIKD-JJ,JJ)=GCS(IIKD-JJ,JJ)+VXI*VXJ
        GCS(IIKD-JJ1,JJ1)=GCS(IIKD-JJ1,JJ1)+VXI*VYJ
        GCS(IIKD-JJ2,JJ2)=GCS(IIKD-JJ2,JJ2)+VXI*VZJ
        GCS(IIKD1-JJ,JJ)=GCS(IIKD1-JJ,JJ)+VYI*VXJ
        GCS(IIKD1-JJ1,JJ1)=GCS(IIKD1-JJ1,JJ1)+VYI*VYJ
        GCS(IIKD1-JJ2,JJ2)=GCS(IIKD1-JJ2,JJ2)+VYI*VZJ
        GCS(IIKD2-JJ,JJ)=GCS(IIKD2-JJ,JJ)+VZI*VXJ
        GCS(IIKD2-JJ1,JJ1)=GCS(IIKD2-JJ1,JJ1)+VZI*VYJ
        GCS(IIKD2-JJ2,JJ2)=GCS(IIKD2-JJ2,JJ2)+VZI*VZJ

C II < KK
        GCS(IIKD-KK,KK)=GCS(IIKD-KK,KK)+VXI*VXK
        GCS(IIKD-KK1,KK1)=GCS(IIKD-KK1,KK1)+VXI*VYK
        GCS(IIKD-KK2,KK2)=GCS(IIKD-KK2,KK2)+VXI*VZK
        GCS(IIKD1-KK,KK)=GCS(IIKD1-KK,KK)+VYI*VXK
        GCS(IIKD1-KK1,KK1)=GCS(IIKD1-KK1,KK1)+VYI*VYK
        GCS(IIKD1-KK2,KK2)=GCS(IIKD1-KK2,KK2)+VYI*VZK
        GCS(IIKD2-KK,KK)=GCS(IIKD2-KK,KK)+VZI*VXK
        GCS(IIKD2-KK1,KK1)=GCS(IIKD2-KK1,KK1)+VZI*VYK
        GCS(IIKD2-KK2,KK2)=GCS(IIKD2-KK2,KK2)+VZI*VZK

C II < LL
        GCS(IIKD-LL,LL)=GCS(IIKD-LL,LL)+VXI*VXL
        GCS(IIKD-LL1,LL1)=GCS(IIKD-LL1,LL1)+VXI*VYL
        GCS(IIKD-LL2,LL2)=GCS(IIKD-LL2,LL2)+VXI*VZL
        GCS(IIKD1-LL,LL)=GCS(IIKD1-LL,LL)+VYI*VXL
        GCS(IIKD1-LL1,LL1)=GCS(IIKD1-LL1,LL1)+VYI*VYL
        GCS(IIKD1-LL2,LL2)=GCS(IIKD1-LL2,LL2)+VYI*VZL
        GCS(IIKD2-LL,LL)=GCS(IIKD2-LL,LL)+VZI*VXL
        GCS(IIKD2-LL1,LL1)=GCS(IIKD2-LL1,LL1)+VZI*VYL
        GCS(IIKD2-LL2,LL2)=GCS(IIKD2-LL2,LL2)+VZI*VZL

C JJ < KK
        GCS(JJKD-KK,KK)=GCS(JJKD-KK,KK)+VXJ*VXK
        GCS(JJKD-KK1,KK1)=GCS(JJKD-KK1,KK1)+VXJ*VYK
        GCS(JJKD-KK2,KK2)=GCS(JJKD-KK2,KK2)+VXJ*VZK
        GCS(JJKD1-KK,KK)=GCS(JJKD1-KK,KK)+VYJ*VXK
        GCS(JJKD1-KK1,KK1)=GCS(JJKD1-KK1,KK1)+VYJ*VYK
        GCS(JJKD1-KK2,KK2)=GCS(JJKD1-KK2,KK2)+VYJ*VZK
        GCS(JJKD2-KK,KK)=GCS(JJKD2-KK,KK)+VZJ*VXK
        GCS(JJKD2-KK1,KK1)=GCS(JJKD2-KK1,KK1)+VZJ*VYK
        GCS(JJKD2-KK2,KK2)=GCS(JJKD2-KK2,KK2)+VZJ*VZK

C       IF (JJ.LT.LL) THEN
        GCS(JJKD-LL,LL)=GCS(JJKD-LL,LL)+VXJ*VXL
        GCS(JJKD-LL1,LL1)=GCS(JJKD-LL1,LL1)+VXJ*VYL
        GCS(JJKD-LL2,LL2)=GCS(JJKD-LL2,LL2)+VXJ*VZL
        GCS(JJKD1-LL,LL)=GCS(JJKD1-LL,LL)+VYJ*VXL
        GCS(JJKD1-LL1,LL1)=GCS(JJKD1-LL1,LL1)+VYJ*VYL
        GCS(JJKD1-LL2,LL2)=GCS(JJKD1-LL2,LL2)+VYJ*VZL
        GCS(JJKD2-LL,LL)=GCS(JJKD2-LL,LL)+VZJ*VXL
        GCS(JJKD2-LL1,LL1)=GCS(JJKD2-LL1,LL1)+VZJ*VYL
        GCS(JJKD2-LL2,LL2)=GCS(JJKD2-LL2,LL2)+VZJ*VZL

C       IF (KK.LT.LL) THEN
        GCS(KKKD-LL,LL)=GCS(KKKD-LL,LL)+VXK*VXL
        GCS(KKKD-LL1,LL1)=GCS(KKKD-LL1,LL1)+VXK*VYL
        GCS(KKKD-LL2,LL2)=GCS(KKKD-LL2,LL2)+VXK*VZL
        GCS(KKKD1-LL,LL)=GCS(KKKD1-LL,LL)+VYK*VXL
        GCS(KKKD1-LL1,LL1)=GCS(KKKD1-LL1,LL1)+VYK*VYL
        GCS(KKKD1-LL2,LL2)=GCS(KKKD1-LL2,LL2)+VYK*VZL
        GCS(KKKD2-LL,LL)=GCS(KKKD2-LL,LL)+VZK*VXL
        GCS(KKKD2-LL1,LL1)=GCS(KKKD2-LL1,LL1)+VZK*VYL
        GCS(KKKD2-LL2,LL2)=GCS(KKKD2-LL2,LL2)+VZK*VZL

C
C end dae
C
   30  CONTINUE

C side chain dihedrals
      IPASS=0

      DO 40 IPHI=1,nres-2
C I= Ca(i-1), J = Cai, K = SC, L = Ca(i+1)
C for glycine
        IF (ITYPE(IPHI+1).EQ.10) GOTO 40
        IPASS=IPASS+1
        I=2*IPHI-1
        J=I+2
        K=I+3
        L=I+4
        II=3*(I-1)+1
        II1=3*(I-1)+2
        II2=3*(I-1)+3
        JJ=3*(J-1)+1
        JJ1=3*(J-1)+2
        JJ2=3*(J-1)+3
        KK=3*(K-1)+1
        KK1=3*(K-1)+2
        KK2=3*(K-1)+3
        LL=3*(L-1)+1
        LL1=3*(L-1)+2
        LL2=3*(L-1)+3
C       IC=ICI(IPHI)
C point BIS is the point where the bisector of theta meets the vector Ca (i-1) -> Ca(i+1)
C (since all virtual bond lengths are equal, this point lies at the midpoint of Ca (i-1) -> Ca(i+1).
C F=Rk-Rj, G=Rj-Rbis, H-Rl-Rbis.
        FX=X(K)-X(J)
        FY=Y(K)-Y(J)
        FZ=Z(K)-Z(J)
        GX=X(J)-(X(I)+X(L))/2.0D0
        GY=Y(J)-(Y(I)+Y(L))/2.0D0
        GZ=Z(J)-(Z(I)+Z(L))/2.0D0
        HX=(X(L)-X(I))/2.0D0
        HY=(Y(L)-Y(I))/2.0D0
        HZ=(Z(L)-Z(I))/2.0D0
C A=F^G, B=H^G
        AX=FY*GZ-FZ*GY
        AY=FZ*GX-FX*GZ
        AZ=FX*GY-FY*GX
        BX=HY*GZ-HZ*GY
        BY=HZ*GX-HX*GZ
        BZ=HX*GY-HY*GX
C RG=|G|, RGR=1/|G|
        RG2=GX*GX+GY*GY+GZ*GZ
        RG=SQRT(RG2)
        RGR=1/RG
C dae for use in evaluating B-matrix
        RF2=FX*FX+FY*FY+FZ*FZ
        RF=SQRT(RF2)
        RFR=1/RF
        RH2=HX*HX+HY*HY+HZ*HZ
        RH=SQRT(RH2)
        RHR=1/RH
C jmc        CSTTWO=(FX*GX+FY*GY+FZ*GZ)*RFR*RGR
        CSTTWO=-(FX*GX+FY*GY+FZ*GZ)*RFR*RGR
        SNTTWO2=1-CSTTWO*CSTTWO
        SNTTWO2R=1/SNTTWO2
        CSTTHREE=(HX*GX+HY*GY+HZ*GZ)*RHR*RGR
C       PRINT *,'COSTHREE = ',CSTTHREE
        SNTTHREE2=1-CSTTHREE*CSTTHREE
C       PRINT *,'SINTHREE2 = ',SNTTHREE2
        SNTTHREE2R=1/SNTTHREE2
        IPHII=12*nphi+9*NTHETA+12*nside+12*(IPASS-1)+1
        SUMPHI=NPHI+NTHETA+nside+IPASS
C
        IF (.NOT.NOCOOR) THEN
          RA2=AX*AX+AY*AY+AZ*AZ
          RB2=BX*BX+BY*BY+BZ*BZ
          RA2R=1/RA2
          RB2R=1/RB2
          RABR=SQRT(RA2R*RB2R)
C CP=cos(phi)
          CP=(AX*BX+AY*BY+AZ*BZ)*RABR
          INTCOORDS(SUMPHI)=ACOS(CP)
          IF (CP.GT.1.0D0) INTCOORDS(SUMPHI)=0.0D0
C         PRINT *,'cp ',CP
C         PRINT *,'acos ',ACOS(CP)
C jmc 
C Test to determine whether the dihedral angle should be > 0 or < 0, from unres (intcor.f) 
          TX=AY*BZ-BY*AZ
          TY=-AX*BZ+BX*AZ
          TZ=AX*BY-BX*AY
          SCALAR=TX*GX+TY*GY+TZ*GZ
C         PRINT *,'SCALAR',SCALAR
          IF (SCALAR.GT.0.0D0) INTCOORDS(SUMPHI)=-INTCOORDS(SUMPHI)
C         PRINT *,'COORDS = ',COORDS(SUMPHI)
        ENDIF
C dae calculate B-matrix elements

C I
        DUMMY=RHR*RHR*RGR*RGR*SNTTHREE2R*(RG-RH*CSTTHREE)
        DUMMY2=RFR*RGR*RGR*SNTTWO2R*CSTTWO
        VAL(IPHII)=(-BX*DUMMY+AX*DUMMY2)/2.0D0
        INDX(IPHII)=SUMPHI
        JNDX(IPHII)=II
        BEEMATRIX(SUMPHI,II)=(-BX*DUMMY+AX*DUMMY2)/2.0D0
C
        VAL(IPHII+1)=(-BY*DUMMY+AY*DUMMY2)/2.0D0
        INDX(IPHII+1)=SUMPHI
        JNDX(IPHII+1)=II+1
        BEEMATRIX(SUMPHI,II+1)=(-BY*DUMMY+AY*DUMMY2)/2.0D0
C
        VAL(IPHII+2)=(-BZ*DUMMY+AZ*DUMMY2)/2.0D0
        INDX(IPHII+2)=SUMPHI
        JNDX(IPHII+2)=II+2
        BEEMATRIX(SUMPHI,II+2)=(-BZ*DUMMY+AZ*DUMMY2)/2.0D0

C J
        DUMMY=RFR*RFR*RGR*RGR*SNTTWO2R*(RG-RF*CSTTWO)
        DUMMY2=RHR*RGR*RGR*SNTTHREE2R*CSTTHREE
C
        VAL(IPHII+3)=AX*DUMMY-BX*DUMMY2
        INDX(IPHII+3)=SUMPHI
        JNDX(IPHII+3)=JJ
        BEEMATRIX(SUMPHI,JJ)=AX*DUMMY-BX*DUMMY2
C
        VAL(IPHII+4)=AY*DUMMY-BY*DUMMY2
        INDX(IPHII+4)=SUMPHI
        JNDX(IPHII+4)=JJ+1
        BEEMATRIX(SUMPHI,JJ+1)=AY*DUMMY-BY*DUMMY2
C
        VAL(IPHII+5)=AZ*DUMMY-BZ*DUMMY2
        INDX(IPHII+5)=SUMPHI
        JNDX(IPHII+5)=JJ+2
        BEEMATRIX(SUMPHI,JJ+2)=AZ*DUMMY-BZ*DUMMY2

C K: side chain centroid.
        DUMMY=RFR*RFR*RGR*SNTTWO2R
        VAL(IPHII+6)=-AX*DUMMY
        INDX(IPHII+6)=SUMPHI
        JNDX(IPHII+6)=KK
        BEEMATRIX(SUMPHI,KK)=-AX*DUMMY
C
        VAL(IPHII+7)=-AY*DUMMY
        INDX(IPHII+7)=SUMPHI
        JNDX(IPHII+7)=KK+1
        BEEMATRIX(SUMPHI,KK+1)=-AY*DUMMY
C
        VAL(IPHII+8)=-AZ*DUMMY
        INDX(IPHII+8)=SUMPHI
        JNDX(IPHII+8)=KK+2
        BEEMATRIX(SUMPHI,KK+2)=-AZ*DUMMY

C L
        DUMMY=RHR*RHR*RGR*SNTTHREE2R
        DUMMY2=RHR*RHR*RGR*RGR*SNTTHREE2R*(RG-RH*CSTTHREE)
        DUMMY3=RFR*RGR*RGR*SNTTWO2R*CSTTWO
        VAL(IPHII+9)=BX*DUMMY+(-BX*DUMMY2+AX*DUMMY3)/2.0D0
        INDX(IPHII+9)=SUMPHI
        JNDX(IPHII+9)=LL
        BEEMATRIX(SUMPHI,LL)=BX*DUMMY+(-BX*DUMMY2+AX*DUMMY3)/2.0D0
C
        VAL(IPHII+10)=BY*DUMMY+(-BY*DUMMY2+AY*DUMMY3)/2.0D0
        INDX(IPHII+10)=SUMPHI
        JNDX(IPHII+10)=LL+1
        BEEMATRIX(SUMPHI,LL+1)=BY*DUMMY+(-BY*DUMMY2+AY*DUMMY3)/2.0D0
C
        VAL(IPHII+11)=BZ*DUMMY+(-BZ*DUMMY2+AZ*DUMMY3)/2.0D0
        INDX(IPHII+11)=SUMPHI
        JNDX(IPHII+11)=LL+2
        BEEMATRIX(SUMPHI,LL+2)=BZ*DUMMY+(-BZ*DUMMY2+AZ*DUMMY3)/2.0D0

C calculate GEE for this ic
        IIKD=II+KD1
        IIKD1=II1+KD1
        IIKD2=II2+KD1
        JJKD=JJ+KD1
        JJKD1=JJ1+KD1
        JJKD2=JJ2+KD1
        KKKD=KK+KD1
        KKKD1=KK1+KD1
        KKKD2=KK2+KD1
        LLKD=LL+KD1
        LLKD1=LL1+KD1
        LLKD2=LL2+KD1
C
        VXI=VAL(IPHII)
        VYI=VAL(IPHII+1)
        VZI=VAL(IPHII+2)
        GCS(KD1,II)=GCS(KD1,II)+VXI**2
        GCS(KD1+II-II1,II1)=GCS(KD1+II-II1,II1)+VXI*VYI
        GCS(KD1+II-II2,II2)=GCS(KD1+II-II2,II2)+VXI*VZI
        GCS(KD1,II1)=GCS(KD1,II1)+VYI**2
        GCS(KD1+II1-II2,II2)=GCS(KD1+II1-II2,II2)+VYI*VZI
        GCS(KD1,II2)=GCS(KD1,II2)+VZI**2
C
        VXJ=VAL(IPHII+3)
        VYJ=VAL(IPHII+4)
        VZJ=VAL(IPHII+5)
        GCS(KD1,JJ)=GCS(KD1,JJ)+VXJ**2
        GCS(KD1+JJ-JJ1,JJ1)=GCS(KD1+JJ-JJ1,JJ1)+VXJ*VYJ
        GCS(KD1+JJ-JJ2,JJ2)=GCS(KD1+JJ-JJ2,JJ2)+VXJ*VZJ
        GCS(KD1,JJ1)=GCS(KD1,JJ1)+VYJ**2
        GCS(KD1+JJ1-JJ2,JJ2)=GCS(KD1+JJ1-JJ2,JJ2)+VYJ*VZJ
        GCS(KD1,JJ2)=GCS(KD1,JJ2)+VZJ**2
C
        VXK=VAL(IPHII+6)
        VYK=VAL(IPHII+7)
        VZK=VAL(IPHII+8)
        GCS(KD1,KK)=GCS(KD1,KK)+VXK**2
        GCS(KD1+KK-KK1,KK1)=GCS(KD1+KK-KK1,KK1)+VXK*VYK
        GCS(KD1+KK-KK2,KK2)=GCS(KD1+KK-KK2,KK2)+VXK*VZK
        GCS(KD1,KK1)=GCS(KD1,KK1)+VYK**2
        GCS(KD1+KK1-KK2,KK2)=GCS(KD1+KK1-KK2,KK2)+VYK*VZK
        GCS(KD1,KK2)=GCS(KD1,KK2)+VZK**2
C
        VXL=VAL(IPHII+9)
        VYL=VAL(IPHII+10)
        VZL=VAL(IPHII+11)
        GCS(KD1,LL)=GCS(KD1,LL)+VXL**2
        GCS(KD1+LL-LL1,LL1)=GCS(KD1+LL-LL1,LL1)+VXL*VYL
        GCS(KD1+LL-LL2,LL2)=GCS(KD1+LL-LL2,LL2)+VXL*VZL
        GCS(KD1,LL1)=GCS(KD1,LL1)+VYL**2
        GCS(KD1+LL1-LL2,LL2)=GCS(KD1+LL1-LL2,LL2)+VYL*VZL
        GCS(KD1,LL2)=GCS(KD1,LL2)+VZL**2
C
C       IF (II.LT.JJ) THEN
        GCS(IIKD-JJ,JJ)=GCS(IIKD-JJ,JJ)+VXI*VXJ
        GCS(IIKD-JJ1,JJ1)=GCS(IIKD-JJ1,JJ1)+VXI*VYJ
        GCS(IIKD-JJ2,JJ2)=GCS(IIKD-JJ2,JJ2)+VXI*VZJ
        GCS(IIKD1-JJ,JJ)=GCS(IIKD1-JJ,JJ)+VYI*VXJ
        GCS(IIKD1-JJ1,JJ1)=GCS(IIKD1-JJ1,JJ1)+VYI*VYJ
        GCS(IIKD1-JJ2,JJ2)=GCS(IIKD1-JJ2,JJ2)+VYI*VZJ
        GCS(IIKD2-JJ,JJ)=GCS(IIKD2-JJ,JJ)+VZI*VXJ
        GCS(IIKD2-JJ1,JJ1)=GCS(IIKD2-JJ1,JJ1)+VZI*VYJ
        GCS(IIKD2-JJ2,JJ2)=GCS(IIKD2-JJ2,JJ2)+VZI*VZJ
C       ELSE
C       GCS(JJKD-II,II)=GCS(IIKD-JJ,JJ)
C       GCS(JJKD1-II,II)=GCS(IIKD-JJ1,JJ1)
C       GCS(JJKD2-II,II)=GCS(IIKD-JJ2,JJ2)
C       GCS(JJKD-II1,II1)=GCS(IIKD1-JJ,JJ)
C       GCS(JJKD1-II1,II1)=GCS(IIKD1-JJ1,JJ1)
C       GCS(JJKD2-II1,II1)=GCS(IIKD1-JJ2,JJ2)
C       GCS(JJKD-II2,II2)=GCS(IIKD2-JJ,JJ)
C       GCS(JJKD1-II2,II2)=GCS(IIKD2-JJ1,JJ1)
C       GCS(JJKD2-II2,II2)=GCS(IIKD2-JJ2,JJ2)
C       ENDIF
C
C       IF (II.LT.KK) THEN
        GCS(IIKD-KK,KK)=GCS(IIKD-KK,KK)+VXI*VXK
        GCS(IIKD-KK1,KK1)=GCS(IIKD-KK1,KK1)+VXI*VYK
        GCS(IIKD-KK2,KK2)=GCS(IIKD-KK2,KK2)+VXI*VZK
        GCS(IIKD1-KK,KK)=GCS(IIKD1-KK,KK)+VYI*VXK
        GCS(IIKD1-KK1,KK1)=GCS(IIKD1-KK1,KK1)+VYI*VYK
        GCS(IIKD1-KK2,KK2)=GCS(IIKD1-KK2,KK2)+VYI*VZK
        GCS(IIKD2-KK,KK)=GCS(IIKD2-KK,KK)+VZI*VXK
        GCS(IIKD2-KK1,KK1)=GCS(IIKD2-KK1,KK1)+VZI*VYK
        GCS(IIKD2-KK2,KK2)=GCS(IIKD2-KK2,KK2)+VZI*VZK
C       ELSE
C       GCS(KKKD-II,II)=GCS(IIKD-KK,KK)
C       GCS(KKKD1-II,II)=GCS(IIKD-KK1,KK1)
C       GCS(KKKD2-II,II)=GCS(IIKD-KK2,KK2)
C       GCS(KKKD-II1,II1)=GCS(IIKD1-KK,KK)
C       GCS(KKKD1-II1,II1)=GCS(IIKD1-KK1,KK1)
C       GCS(KKKD2-II1,II1)=GCS(IIKD1-KK2,KK2)
C       GCS(KKKD-II2,II2)=GCS(IIKD2-KK,KK)
C       GCS(KKKD1-II2,II2)=GCS(IIKD2-KK1,KK1)
C       GCS(KKKD2-II2,II2)=GCS(IIKD2-KK2,KK2)
C       ENDIF
C
C       IF (II.LT.LL) THEN
        GCS(IIKD-LL,LL)=GCS(IIKD-LL,LL)+VXI*VXL
        GCS(IIKD-LL1,LL1)=GCS(IIKD-LL1,LL1)+VXI*VYL
        GCS(IIKD-LL2,LL2)=GCS(IIKD-LL2,LL2)+VXI*VZL
        GCS(IIKD1-LL,LL)=GCS(IIKD1-LL,LL)+VYI*VXL
        GCS(IIKD1-LL1,LL1)=GCS(IIKD1-LL1,LL1)+VYI*VYL
        GCS(IIKD1-LL2,LL2)=GCS(IIKD1-LL2,LL2)+VYI*VZL
        GCS(IIKD2-LL,LL)=GCS(IIKD2-LL,LL)+VZI*VXL
        GCS(IIKD2-LL1,LL1)=GCS(IIKD2-LL1,LL1)+VZI*VYL
        GCS(IIKD2-LL2,LL2)=GCS(IIKD2-LL2,LL2)+VZI*VZL
C       ELSE
C       GCS(LLKD-II,II)=GCS(IIKD-LL,LL)
C       GCS(LLKD1-II,II)=GCS(IIKD-LL1,LL1)
C       GCS(LLKD2-II,II)=GCS(IIKD-LL2,LL2)
C       GCS(LLKD-II1,II1)=GCS(IIKD1-LL,LL)
C       GCS(LLKD1-II1,II1)=GCS(IIKD1-LL1,LL1)
C       GCS(LLKD2-II1,II1)=GCS(IIKD1-LL2,LL2)
C       GCS(LLKD-II2,II2)=GCS(IIKD2-LL,LL)
C       GCS(LLKD1-II2,II2)=GCS(IIKD2-LL1,LL1)
C       GCS(LLKD2-II2,II2)=GCS(IIKD2-LL2,LL2)
C       ENDIF
C
C       IF (JJ.LT.KK) THEN
        GCS(JJKD-KK,KK)=GCS(JJKD-KK,KK)+VXJ*VXK
        GCS(JJKD-KK1,KK1)=GCS(JJKD-KK1,KK1)+VXJ*VYK
        GCS(JJKD-KK2,KK2)=GCS(JJKD-KK2,KK2)+VXJ*VZK
        GCS(JJKD1-KK,KK)=GCS(JJKD1-KK,KK)+VYJ*VXK
        GCS(JJKD1-KK1,KK1)=GCS(JJKD1-KK1,KK1)+VYJ*VYK
        GCS(JJKD1-KK2,KK2)=GCS(JJKD1-KK2,KK2)+VYJ*VZK
        GCS(JJKD2-KK,KK)=GCS(JJKD2-KK,KK)+VZJ*VXK
        GCS(JJKD2-KK1,KK1)=GCS(JJKD2-KK1,KK1)+VZJ*VYK
        GCS(JJKD2-KK2,KK2)=GCS(JJKD2-KK2,KK2)+VZJ*VZK
C       ELSE
C       GCS(KKKD-JJ,JJ)=GCS(JJKD-KK,KK)
C       GCS(KKKD1-JJ,JJ)=GCS(JJKD-KK1,KK1)
C       GCS(KKKD2-JJ,JJ)=GCS(JJKD-KK2,KK2)
C       GCS(KKKD-JJ1,JJ1)=GCS(JJKD1-KK,KK)
C       GCS(KKKD1-JJ1,JJ1)=GCS(JJKD1-KK1,KK1)
C       GCS(KKKD2-JJ1,JJ1)=GCS(JJKD1-KK2,KK2)
C       GCS(KKKD-JJ2,JJ2)=GCS(JJKD2-KK,KK)
C       GCS(KKKD1-JJ2,JJ2)=GCS(JJKD2-KK1,KK1)
C       GCS(KKKD2-JJ2,JJ2)=GCS(JJKD2-KK2,KK2)
C       ENDIF
C
C       IF (JJ.LT.LL) THEN
        GCS(JJKD-LL,LL)=GCS(JJKD-LL,LL)+VXJ*VXL
        GCS(JJKD-LL1,LL1)=GCS(JJKD-LL1,LL1)+VXJ*VYL
        GCS(JJKD-LL2,LL2)=GCS(JJKD-LL2,LL2)+VXJ*VZL
        GCS(JJKD1-LL,LL)=GCS(JJKD1-LL,LL)+VYJ*VXL
        GCS(JJKD1-LL1,LL1)=GCS(JJKD1-LL1,LL1)+VYJ*VYL
        GCS(JJKD1-LL2,LL2)=GCS(JJKD1-LL2,LL2)+VYJ*VZL
        GCS(JJKD2-LL,LL)=GCS(JJKD2-LL,LL)+VZJ*VXL
        GCS(JJKD2-LL1,LL1)=GCS(JJKD2-LL1,LL1)+VZJ*VYL
        GCS(JJKD2-LL2,LL2)=GCS(JJKD2-LL2,LL2)+VZJ*VZL
C       ELSE
C       GCS(LLKD-JJ,JJ)=GCS(JJKD-LL,LL)
C       GCS(LLKD1-JJ,JJ)=GCS(JJKD-LL1,LL1)
C       GCS(LLKD2-JJ,JJ)=GCS(JJKD-LL2,LL2)
C       GCS(LLKD-JJ1,JJ1)=GCS(JJKD1-LL,LL)
C       GCS(LLKD1-JJ1,JJ1)=GCS(JJKD1-LL1,LL1)
C       GCS(LLKD2-JJ1,JJ1)=GCS(JJKD1-LL2,LL2)
C       GCS(LLKD-JJ2,JJ2)=GCS(JJKD2-LL,LL)
C       GCS(LLKD1-JJ2,JJ2)=GCS(JJKD2-LL1,LL1)
C       GCS(LLKD2-JJ2,JJ2)=GCS(JJKD2-LL2,LL2)
C       ENDIF
C
C       IF (KK.LT.LL) THEN
        GCS(KKKD-LL,LL)=GCS(KKKD-LL,LL)+VXK*VXL
        GCS(KKKD-LL1,LL1)=GCS(KKKD-LL1,LL1)+VXK*VYL
        GCS(KKKD-LL2,LL2)=GCS(KKKD-LL2,LL2)+VXK*VZL
        GCS(KKKD1-LL,LL)=GCS(KKKD1-LL,LL)+VYK*VXL
        GCS(KKKD1-LL1,LL1)=GCS(KKKD1-LL1,LL1)+VYK*VYL
        GCS(KKKD1-LL2,LL2)=GCS(KKKD1-LL2,LL2)+VYK*VZL
        GCS(KKKD2-LL,LL)=GCS(KKKD2-LL,LL)+VZK*VXL
        GCS(KKKD2-LL1,LL1)=GCS(KKKD2-LL1,LL1)+VZK*VYL
        GCS(KKKD2-LL2,LL2)=GCS(KKKD2-LL2,LL2)+VZK*VZL
C       ELSE
C       GCS(LLKD-KK,KK)=GCS(KKKD-LL,LL)
C       GCS(LLKD1-KK,KK)=GCS(KKKD-LL1,LL1)
C       GCS(LLKD2-KK,KK)=GCS(KKKD-LL2,LL2)
C       GCS(LLKD-KK1,KK1)=GCS(KKKD1-LL,LL)
C       GCS(LLKD1-KK1,KK1)=GCS(KKKD1-LL1,LL1)
C       GCS(LLKD2-KK1,KK1)=GCS(KKKD1-LL2,LL2)
C       GCS(LLKD-KK2,KK2)=GCS(KKKD2-LL,LL)
C       GCS(LLKD1-KK2,KK2)=GCS(KKKD2-LL1,LL1)
C       GCS(LLKD2-KK2,KK2)=GCS(KKKD2-LL2,LL2)
C       ENDIF
C
C end dae
C
   40 CONTINUE

C for main chain virtual bonds:
      DO 50 MM=1,nres-1
        I=2*MM-1
        J=2*MM+1
C MM+2 from Ca, sc, Ca, sc ... in X,Y,Z
        II=3*(I-1)+1
        II1=3*(I-1)+2
        II2=3*(I-1)+3
        JJ=3*(J-1)+1
        JJ1=3*(J-1)+2
        JJ2=3*(J-1)+3
C      print *,'bond I J',I,J
C ? charmm      IC=ICB(MM)
C ICB array contains pointers for bond parameters
C       IF(IC.EQ.0) GOTO 50
C ? charmm      IF(CBC(IC).EQ.0) GOTO 50
C CBC is force constant of bond
        RX=X(I)-X(J)
        RY=Y(I)-Y(J)
        RZ=Z(I)-Z(J)
        S2=RX*RX + RY*RY + RZ*RZ
        S=SQRT(S2)
        IF (.NOT.NOCOOR) THEN
C dae put S into COORDS: it is internal coordinate
          INTCOORDS(MM+nphi+ntheta+2*nside)=S
        ENDIF
C dae calculate B-matrix elements for this bond
C dae VAL is 1-D array of non-zero B elements
C INDX is array of the ic(row) indices of VAL
C JNDX is array of the cart(column) indices of VAL
        SR=1/S
        MMM=12*nphi+9*ntheta+24*nside+6*(MM-1)+1
        SUMPHI=nphi+ntheta+2*nside+MM
        RXSR=RX*SR
        RYSR=RY*SR
        RZSR=RZ*SR
        VAL(MMM)=RXSR
        INDX(MMM)=SUMPHI
        JNDX(MMM)=II
        BEEMATRIX(SUMPHI,II)=RXSR
C
        VAL(MMM+1)=RYSR
        INDX(MMM+1)=SUMPHI
        JNDX(MMM+1)=II1
        BEEMATRIX(SUMPHI,II1)=RYSR
C       
        VAL(MMM+2)=RZSR
        INDX(MMM+2)=SUMPHI
        JNDX(MMM+2)=II2
        BEEMATRIX(SUMPHI,II2)=RZSR
C      
        VAL(MMM+3)=-RXSR
        INDX(MMM+3)=SUMPHI
        JNDX(MMM+3)=JJ
        BEEMATRIX(SUMPHI,JJ)=-RXSR
C     
        VAL(MMM+4)=-RYSR
        INDX(MMM+4)=SUMPHI
        JNDX(MMM+4)=JJ1
        BEEMATRIX(SUMPHI,JJ1)=-RYSR
C    
        VAL(MMM+5)=-RZSR
        INDX(MMM+5)=SUMPHI
        JNDX(MMM+5)=JJ2
        BEEMATRIX(SUMPHI,JJ2)=-RZSR
C calculate GEE for this ic
        VXI=VAL(MMM)
        VYI=VAL(MMM+1)
        VZI=VAL(MMM+2)
        GCS(KD1,II)=GCS(KD1,II)+VXI**2
        GCS(KD1+II-II1,II1)=GCS(KD1+II-II1,II1)+VXI*VYI
        GCS(KD1+II-II2,II2)=GCS(KD1+II-II2,II2)+VXI*VZI
        GCS(KD1,II1)=GCS(KD1,II1)+VYI**2
        GCS(KD1+II1-II2,II2)=GCS(KD1+II1-II2,II2)+VYI*VZI
        GCS(KD1,II2)=GCS(KD1,II2)+VZI**2
C
        VXJ=VAL(MMM+3)
        VYJ=VAL(MMM+4)
        VZJ=VAL(MMM+5)
        GCS(KD1,JJ)=GCS(KD1,JJ)+VXJ**2
        GCS(KD1+JJ-JJ1,JJ1)=GCS(KD1+JJ-JJ1,JJ1)+VXJ*VYJ
        GCS(KD1+JJ-JJ2,JJ2)=GCS(KD1+JJ-JJ2,JJ2)+VXJ*VZJ
        GCS(KD1,JJ1)=GCS(KD1,JJ1)+VYJ**2
        GCS(KD1+JJ1-JJ2,JJ2)=GCS(KD1+JJ1-JJ2,JJ2)+VYJ*VZJ
        GCS(KD1,JJ2)=GCS(KD1,JJ2)+VZJ**2
C
        IIKD=II+KD1
        IIKD1=II1+KD1
        IIKD2=II2+KD1
        JJKD=JJ+KD1
        JJKD1=JJ1+KD1
        JJKD2=JJ2+KD1
C
C       IF (II.LT.JJ) THEN
        GCS(IIKD-JJ,JJ)=GCS(IIKD-JJ,JJ)+VXI*VXJ
        GCS(IIKD-JJ1,JJ1)=GCS(IIKD-JJ1,JJ1)+VXI*VYJ
        GCS(IIKD-JJ2,JJ2)=GCS(IIKD-JJ2,JJ2)+VXI*VZJ
        GCS(IIKD1-JJ,JJ)=GCS(IIKD1-JJ,JJ)+VYI*VXJ
        GCS(IIKD1-JJ1,JJ1)=GCS(IIKD1-JJ1,JJ1)+VYI*VYJ
        GCS(IIKD1-JJ2,JJ2)=GCS(IIKD1-JJ2,JJ2)+VYI*VZJ
        GCS(IIKD2-JJ,JJ)=GCS(IIKD2-JJ,JJ)+VZI*VXJ
        GCS(IIKD2-JJ1,JJ1)=GCS(IIKD2-JJ1,JJ1)+VZI*VYJ
        GCS(IIKD2-JJ2,JJ2)=GCS(IIKD2-JJ2,JJ2)+VZI*VZJ
C end dae
   50  CONTINUE
       
      IPASS=0
C for side chain bonds:
      DO 60 MM=1,nres
        IF (ITYPE(MM).EQ.10) GOTO 60
        IPASS=IPASS+1
        I=2*MM-1
        J=2*MM
        II=3*(I-1)+1
        II1=3*(I-1)+2
        II2=3*(I-1)+3
        JJ=3*(J-1)+1
        JJ1=3*(J-1)+2
        JJ2=3*(J-1)+3
C      print *,'bond I J',I,J
C ? charmm      IC=ICB(MM)
C ICB array contains pointers for bond parameters
C       IF(IC.EQ.0) GOTO 10
C ? charmm      IF(CBC(IC).EQ.0) GOTO 10
C CBC is force constant of bond
        RX=X(I)-X(J)
        RY=Y(I)-Y(J)
        RZ=Z(I)-Z(J)
        S2=RX*RX + RY*RY + RZ*RZ
        S=SQRT(S2)
        IF (.NOT.NOCOOR) THEN
C dae put S into COORDS: it is internal coordinate
          INTCOORDS(nphi+ntheta+2*nside+nres-1+IPASS)=S
        ENDIF
C dae calculate B-matrix elements for this bond
C dae VAL is 1-D array of non-zero B elements
C INDX is array of the ic(row) indices of VAL
C JNDX is array of the cart(column) indices of VAL
        SR=1/S
        MMM=12*nphi+9*ntheta+24*nside+6*(nres-1)+6*(IPASS-1)+1
        SUMPHI=nphi+ntheta+2*nside+nres-1+IPASS
        RXSR=RX*SR
        RYSR=RY*SR
        RZSR=RZ*SR
        VAL(MMM)=RXSR
        INDX(MMM)=SUMPHI
        JNDX(MMM)=II
        BEEMATRIX(SUMPHI,II)=RXSR
C
        VAL(MMM+1)=RYSR
        INDX(MMM+1)=SUMPHI
        JNDX(MMM+1)=II1
        BEEMATRIX(SUMPHI,II1)=RYSR
C       
        VAL(MMM+2)=RZSR
        INDX(MMM+2)=SUMPHI
        JNDX(MMM+2)=II2
        BEEMATRIX(SUMPHI,II2)=RZSR
C      
        VAL(MMM+3)=-RXSR
        INDX(MMM+3)=SUMPHI
        JNDX(MMM+3)=JJ
        BEEMATRIX(SUMPHI,JJ)=-RXSR
C     
        VAL(MMM+4)=-RYSR
        INDX(MMM+4)=SUMPHI
        JNDX(MMM+4)=JJ1
        BEEMATRIX(SUMPHI,JJ1)=-RYSR
C    
        VAL(MMM+5)=-RZSR
        INDX(MMM+5)=SUMPHI
        JNDX(MMM+5)=JJ2
        BEEMATRIX(SUMPHI,JJ2)=-RZSR
C calculate GEE for this ic
        VXI=VAL(MMM)
        VYI=VAL(MMM+1)
        VZI=VAL(MMM+2)
        GCS(KD1,II)=GCS(KD1,II)+VXI**2
        GCS(KD1+II-II1,II1)=GCS(KD1+II-II1,II1)+VXI*VYI
        GCS(KD1+II-II2,II2)=GCS(KD1+II-II2,II2)+VXI*VZI
        GCS(KD1,II1)=GCS(KD1,II1)+VYI**2
        GCS(KD1+II1-II2,II2)=GCS(KD1+II1-II2,II2)+VYI*VZI
        GCS(KD1,II2)=GCS(KD1,II2)+VZI**2
C
        VXJ=VAL(MMM+3)
        VYJ=VAL(MMM+4)
        VZJ=VAL(MMM+5)
        GCS(KD1,JJ)=GCS(KD1,JJ)+VXJ**2
        GCS(KD1+JJ-JJ1,JJ1)=GCS(KD1+JJ-JJ1,JJ1)+VXJ*VYJ
        GCS(KD1+JJ-JJ2,JJ2)=GCS(KD1+JJ-JJ2,JJ2)+VXJ*VZJ
        GCS(KD1,JJ1)=GCS(KD1,JJ1)+VYJ**2
        GCS(KD1+JJ1-JJ2,JJ2)=GCS(KD1+JJ1-JJ2,JJ2)+VYJ*VZJ
        GCS(KD1,JJ2)=GCS(KD1,JJ2)+VZJ**2
C
        IIKD=II+KD1
        IIKD1=II1+KD1
        IIKD2=II2+KD1
        JJKD=JJ+KD1
        JJKD1=JJ1+KD1
        JJKD2=JJ2+KD1
C
C       IF (II.LT.JJ) THEN
        GCS(IIKD-JJ,JJ)=GCS(IIKD-JJ,JJ)+VXI*VXJ
        GCS(IIKD-JJ1,JJ1)=GCS(IIKD-JJ1,JJ1)+VXI*VYJ
        GCS(IIKD-JJ2,JJ2)=GCS(IIKD-JJ2,JJ2)+VXI*VZJ
        GCS(IIKD1-JJ,JJ)=GCS(IIKD1-JJ,JJ)+VYI*VXJ
        GCS(IIKD1-JJ1,JJ1)=GCS(IIKD1-JJ1,JJ1)+VYI*VYJ
        GCS(IIKD1-JJ2,JJ2)=GCS(IIKD1-JJ2,JJ2)+VYI*VZJ
        GCS(IIKD2-JJ,JJ)=GCS(IIKD2-JJ,JJ)+VZI*VXJ
        GCS(IIKD2-JJ1,JJ1)=GCS(IIKD2-JJ1,JJ1)+VZI*VYJ
        GCS(IIKD2-JJ2,JJ2)=GCS(IIKD2-JJ2,JJ2)+VZI*VZJ
C end dae
   60  CONTINUE
C
C add lower diagonal elements to upper diagonal. Avoid doubling diagonal
C by looping only to J1-1
       DO J1=1,KD
         DO I1=1,J1-1
           GCS(KD+1+I1-J1,J1)=GCS(KD+1+I1-J1,J1)+GCS(KD+1+J1-I1,I1)
         ENDDO
       ENDDO
       DO J1=KD+1,NCART
         DO I1=J1-KD,J1-1
           GCS(KD+1+I1-J1,J1)=GCS(KD+1+I1-J1,J1)+GCS(KD+1+J1-I1,I1)
         ENDDO
       ENDDO
       RETURN
       END
C
C subroutine to do loads of coordinate transformations to obtain 
C normal mode frequencies via internal coordinate representation.
C calls GETBEE (renamed UNRSGETBEE) to calculate B-matrix
C
      SUBROUTINE INTSECDET(COORDS,NCART,KD,NNZ,NINTB,LAMBDA)
      USE COMMONS
      USE MODHESS
      USE MODUNRES
      IMPLICIT NONE
C
C NNZ is number of non-zero elements in B-matrix
C NINTB is number of internal coordinates (including bond lengths)
C NINTS in commons.f90 is number of main chain dihedrals and angles, and side chain dihedrals & angles 
C not including glycines
C NCART is number of cartesians
C
      INTEGER NNZ,NCART,KD,K,NINTB
C
C     REAL*8 COORDS(6*NATOMS),INTCOORDS(NINTB),VAL(NNZ),GCS2(2*KD+1,NCART)
      REAL*8 COORDS(NCART),INTCOORDS(NINTB),VAL(NNZ),GCS2(2*KD+1,NCART)
      INTEGER INDX(NNZ),JNDX(NNZ),I1,J1,K1,INFO
      REAL*8 BEEMATRIX(NINTB,NCART)
C     REAL*8 BEENUMMATRIX(NINTB,NCART)
      REAL*8 GEE(NINTS,NINTS)
C     REAL*8 WORK(4*NINTS),LAMBDA(NCART),EVECIM(NINTS),CRAP(NINTS),CRAP2(NINTS)
      REAL*8 LAMBDA(NCART),EVECIM(NINTS),CRAP(NINTS),CRAP2(NINTS)
      REAL*8 KVEC(NINTS),WORK1(3*NINTS),A(NINTS,NINTS),TMPA(NINTS,NINTS)
      REAL*8 TMPG(NINTS,NINTS)
      REAL*8 NEWATMASS(NATOMS),P
C
C  Assign enough memory to WORK for a blocksize of 32 to be possible.
C  This is for DSYEVR.
C
      INTEGER ILWORK, LWORK, NFOUND, ISUPPZ(2*NINTS), IWORK(10*NINTS), INFO
      DOUBLE PRECISION WORK(33*NINTS), ABSTOL, DLAMCH, ZZWORK(NINTS,NINTS)
      ABSTOL=DLAMCH( 'Safe  minimum' )
      LWORK=33*NINTS
c     ILWORK=33*NINTS
      ILWORK=10*NINTS
C
C     print *,'transform NINTB NCART NNZ',NINTB,NCART,NNZ
C
C Initialisation
C
      GEE=0.0D0
      TMPA=0.0D0
      A=0.0D0
      LAMBDA=0.0D0
C
C get B-matrix
C
      CALL UNRSGETBEE(BEEMATRIX,COORDS,VAL,INDX,JNDX,NINTB,NCART,NNZ,GCS2,.TRUE.,KD,INTCOORDS)

C jmc for testing
C      CALL BEENUM(NINTB,NCART,COORDS,BEENUMMATRIX,NNZ,KD)
C      DO I1=1,NINTB-nres+1-nside
C         DO J1=1,NCART
C            IF (ABS(BEEMATRIX(I1,J1)-BEENUMMATRIX(I1,J1)).GE.0.1D-5) PRINT *,'Matrix elements',
c    &     I1,J1,'differ (anal, num):',BEEMATRIX(I1,J1),BEENUMMATRIX(I1,J1)
C         END DO
C      END DO

      DO I1=1,NINTS
C jmc         DO J1=1,NINTS ! since only need upper triangle of GEE to pass to DSYEV
         DO J1=I1,NINTS
            DO K1=1,NATOMS
               GEE(I1,J1)=GEE(I1,J1)+BEEMATRIX(I1,3*(K1-1)+1)*BEEMATRIX(J1,3*(K1-1)+1)
     &         *1.0D3/ATMASS(K1) ! NOTE g -> kg !!!
               GEE(I1,J1)=GEE(I1,J1)+BEEMATRIX(I1,3*(K1-1)+2)*BEEMATRIX(J1,3*(K1-1)+2)
     &         *1.0D3/ATMASS(K1)
               GEE(I1,J1)=GEE(I1,J1)+BEEMATRIX(I1,3*(K1-1)+3)*BEEMATRIX(J1,3*(K1-1)+3)
     &         *1.0D3/ATMASS(K1)
            END DO
            IF (I1.EQ.1.AND.J1.EQ.1) PRINT *,'printing to make matrix multiplication work...'
         END DO
      END DO

C jmc test1 - is GEE symmetric? Commented out since only calculate upper triangle of GEE now.  Tested, OK.
C     DO I1=1,NINTS
C         DO J1=1,NINTS
C            IF (ABS(GEE(I1,J1)-GEE(J1,I1)).GE.1.0D-5) PRINT *,'Gee non symmetric ',I1,J1,GEE(I1,J1),GEE(J1,I1)
C         END DO
C     END DO
 
C
C Wilson FG method: entirely equivalent, negligible time difference compared to doing the transformations explicitly, 
C step-by-step.  ***Need to make the full GEE matrix above, not just the upper triangle***
C
C     DO I1=1,NINTS
C        DO J1=1,NINTS
C           HESS(I1,J1)=4.184D0*1.0D3*HESS(I1,J1)
C        ENDDO
C     ENDDO
C     A=matmul(GEE,HESS)
C     CALL DGEEV('N','N',NINTS,A,NINTS,LAMBDA(1:NINTS),EVECIM,CRAP,NINTS,CRAP2,NINTS,WORK,4*NINTS,INFO)
C     IF (INFO.NE.0) PRINT*,'WARNING - INFO=',INFO,' in DGEEV'
C     GOTO 112
C Wilson GF method

C
C Find eigenvalues and eigenvectors of GEE, which are eigenvalues^-1 and the eigenvectors themselves of GEE^-1.
C Previously used DSYEV (below) but switched to DSYEVR as the former gave substantially different results with
C the sunperf libraries and lapack.f/blas.f compiled in.  DSYEVR is not perfect, but better...
C
       CALL DSYEVR('V','A','U',NINTS,GEE,NINTS,0.0D0,1.0D0,1,2,ABSTOL,NFOUND,KVEC,ZZWORK,NINTS,ISUPPZ,WORK,
     1            LWORK, IWORK, ILWORK, INFO )

       IF (INFO.NE.0) PRINT*,'WARNING - INFO=',INFO,' in DSYEVR'
C
C NOTE now ZZWORK is the matrix of eigenvectors, was GEE itself from DSYEV!!!
C
c     CALL DSYEV('V','U',NINTS,GEE,NINTS,KVEC,WORK1,3*NINTS,INFO)
c     IF (INFO.NE.0) PRINT*,'WARNING - INFO=',INFO,' in DSYEV'

c     DO J1=1,NINTS
c        PRINT *,'K ',KVEC(J1)
c     END DO

C now do A = LUHU^TL
C where L_ij = delta_ij/sqrt(k_i), where k_i are the eigenvalues of GEE^-1, diagonal)
 
c     do i1=1,nints
c        ll(i1,i1)=dsqrt(kvec(i1))
c     enddo
c     tmp1=matmul(transpose(gee),ll)
c     tmp2=matmul(hess(:nints,:nints),tmp1)
c     tmp1=matmul(gee,tmp2)
c     A=matmul(ll,tmp1)
c     A=4.184D0*A

      DO I1=1,NINTS
         DO J1=1,NINTS
            DO K1=1,NINTS
C              TMPA(I1,J1)=TMPA(I1,J1)+GEE(I1,K1)*HESS(K1,J1)*4.184D0*DSQRT(KVEC(I1))
               TMPA(I1,J1)=TMPA(I1,J1)+ZZWORK(I1,K1)*HESS(K1,J1)*4.184D0*1.0D3*DSQRT(KVEC(I1))
C factor of 4.184D0*1.0D3 from kcal -> J in HESS
            END DO
         END DO
      END DO

      DO I1=1,NINTS
         DO J1=1,NINTS
            DO K1=1,NINTS
C              A(I1,J1)=A(I1,J1)+TMPA(I1,K1)*GEE(J1,K1)*DSQRT(KVEC(J1))
               A(I1,J1)=A(I1,J1)+TMPA(I1,K1)*ZZWORK(J1,K1)*DSQRT(KVEC(J1))
            END DO
         END DO
      END DO

      CALL DGEEV('N','N',NINTS,A,NINTS,LAMBDA(1:NINTS),EVECIM,CRAP,NINTS,CRAP2,NINTS,WORK,4*NINTS,INFO)
      IF (INFO.NE.0) PRINT*,'WARNING - INFO=',INFO,' in DGEEV'

112   CONTINUE
C     DO I1=1,NINTS
C        PRINT *,'EVALMIN ',LAMBDA(I1),I1
C        PRINT *,'EVALMIN ',EVECIM(I1),I1
C     END DO

C jmc sort lambdas
      DO 30 I1=1,NINTS-1
         K1=I1
         P=LAMBDA(I1)
         DO 10 J1=I1+1,NINTS
            IF (LAMBDA(J1).GE.P) THEN
               K1=J1
               P=LAMBDA(J1)
            ENDIF
10       CONTINUE
         IF (K1.NE.I1) THEN
            LAMBDA(K1)=LAMBDA(I1)
            LAMBDA(I1)=P
         ENDIF
30    CONTINUE

      RETURN
      END

      SUBROUTINE BEENUM(NINTB,NCART,COORDS,BEENUMMATRIX,NNZ,KD)
      USE MODUNRES
      IMPLICIT NONE

      INTEGER NNZ,NINTB,NCART,KD,K
C
C jmc      REAL*8 COORDS(6*NATOMS)
      REAL*8 COORDS(3*NCART)
      REAL*8 VAL(NNZ)
      INTEGER INDX(NNZ),JNDX(NNZ)
      INTEGER I1,J1
      REAL*8 GCS(KD+1,NCART),GCS2(2*KD+1,NCART)
C
      REAL*8 BEEMATRIX(NINTB,NCART),BEENUMMATRIX(NINTB,NCART)
      REAL*8 INTCOORDS(NINTB),NEWINTCOORDS(NINTB)

      INTCOORDS=0.0D0
      NEWINTCOORDS=0.0D0

      DO I1=1,NCART
         COORDS(I1)=COORDS(I1)+0.1D-5
         CALL UNRSGETBEE(BEEMATRIX,COORDS,VAL,INDX,JNDX,NINTB,NCART,NNZ,GCS2,.FALSE.,KD,INTCOORDS)

         COORDS(I1)=COORDS(I1)-0.2D-5
         CALL UNRSGETBEE(BEEMATRIX,COORDS,VAL,INDX,JNDX,NINTB,NCART,NNZ,GCS2,.FALSE.,KD,NEWINTCOORDS)

         COORDS(I1)=COORDS(I1)+0.1D-5

         DO J1=1,NINTB-nres+1-nside
c           PRINT *,'INTCOORDS1 ',INTCOORDS(J1),J1
c           PRINT *,'INTCOORDS2 ',NEWINTCOORDS(J1),J1
            BEENUMMATRIX(J1,I1)=(INTCOORDS(J1)-NEWINTCOORDS(J1))/0.2D-5
            IF ((INTCOORDS(J1)-NEWINTCOORDS(J1)).GE.3.1D0) THEN
              BEENUMMATRIX(J1,I1)=(INTCOORDS(J1)-NEWINTCOORDS(J1)-(2*3.141592653589793D0))/0.2D-5
            ELSE IF ((INTCOORDS(J1)-NEWINTCOORDS(J1)).LE.-3.1D0) THEN
              BEENUMMATRIX(J1,I1)=(INTCOORDS(J1)-NEWINTCOORDS(J1)+(2*3.141592653589793D0))/0.2D-5
            END IF
         END DO
      END DO

      RETURN
      END SUBROUTINE

