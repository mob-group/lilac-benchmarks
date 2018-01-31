!
! Calculate the Forces and energies for a dihedral angle potential
! Based on Paul Whitford's routine from SBM.f.
!
      SUBROUTINE TWIST(COORDS,NATOMS,GRAD,ENERGY,GTEST)
      USE KEY, ONLY : TWISTF, ITWIST, JTWIST, KTWIST, LTWIST, TWISTREF
      IMPLICIT NONE
      INTEGER I, J, NATOMS
      DOUBLE PRECISION X(NATOMS),Y(NATOMS),Z(NATOMS),GRAD(3*NATOMS),ENERGY,COORDS(3*NATOMS)
      LOGICAL GTEST
      DOUBLE PRECISION PK
      INTEGER IP, JP, KP, LP 
      INTEGER I3, J3, K3, L3
      DOUBLE PRECISION  XIJ,YIJ,ZIJ,XKJ,YKJ,ZKJ,XKL,YKL,ZKL,RKJ,DX,DY,
     & DZ, GX,GY,GZ,CT,CPHI,SPHI,FXI,FYI,
     & DF,Z1,Z2,Z12,CT0,CT1,AP0,AP1,S,HGoverG,FGoverG,A1

      DOUBLE PRECISION ZERO,ONE
      DOUBLE PRECISION TT1, TT2, TT3, TT4, TT1X,TT1Y,TT1Z,TT2X,TT2Y,
     & TT2Z, TT3X, TT3Y, TT3Z, TT4X, TT4Y, TT4Z
      DATA ZERO,ONE /0.d0,1.d0/

      DOUBLE PRECISION PI

      PI = 3.14159265358979323846264338327950288419716939937510

      DO I = 1, NATOMS
         J = (I-1)*3
         X(i) = COORDS(J+1)
         Y(i) = COORDS(J+2)
         Z(i) = COORDS(J+3)
      ENDDO

      PK=TWISTF
      IP=ITWIST
      JP=JTWIST
      KP=KTWIST
      LP=LTWIST

      I3 = IP
      J3 = JP
      K3 = KP
      L3 = LP

      XIJ = X(I3)-X(J3)
      YIJ = Y(I3)-Y(J3)
      ZIJ = Z(I3)-Z(J3)
      XKJ = X(K3)-X(J3)
      YKJ = Y(K3)-Y(J3)
      ZKJ = Z(K3)-Z(J3)
! RKJ = size of kj vector
      RKJ = SQRT(XKJ**2+YKJ**2+ZKJ**2)
      XKL = X(K3)-X(L3)
      YKL = Y(K3)-Y(L3)
      ZKL = Z(K3)-Z(L3)                                  

! r_ij.r_kj / |r_kj|
      FGoverG=-(XIJ*XKJ+YIJ*YKJ+ZIJ*ZKJ)/RKJ
! r_kl.r_kj / |r_kj|
      HGoverG=(XKL*XKJ+YKL*YKJ+ZKL*ZKJ)/RKJ

! DX is the M vector and G is the N vector

! D = r_ij x r_kj
      DX = YIJ*ZKJ-ZIJ*YKJ
      DY = ZIJ*XKJ-XIJ*ZKJ
      DZ = XIJ*YKJ-YIJ*XKJ
! G = r_kl x r_kj
      GX = ZKJ*YKL-YKJ*ZKL
      GY = XKJ*ZKL-ZKJ*XKL
      GZ = YKJ*XKL-XKJ*YKL

      FXI = SQRT(DX*DX+DY*DY+DZ*DZ)
      FYI = SQRT(GX*GX+GY*GY+GZ*GZ)
! CT = D x G
      CT = DX*GX+DY*GY+DZ*GZ

! Z1 = 1/|D|
      Z1 = 1.0D0/FXI
! Z2 = 1/|G|
      Z2 = 1.0D0/FYI
! Z12 = 1/|D||G|
      Z12 = Z1*Z2
! CT0 = min(1, D x G / |D||G|)
      CT0 = MIN(one,CT*Z12)
! CT1 = max(-1, CT0)
      CT1 = MAX(-one,CT0)
! S = r_kj . (G x D)
      S = XKJ*(DZ*GY-DY*GZ)+YKJ*(DX*GZ-DZ*GX)+ZKJ*(DY*GX-DX*GY)
! AP0 = acos(CT1)
      AP0 = ACOS(CT1)
! AP1 = pi - AP0 * sign of S
      AP1 = PI-SIGN(AP0,S)

! CT = AP1
      CT = AP1
! CPHI = COS(AP1)
      CPHI = COS(AP1)
! SPHI = SIN(AP1)
      SPHI = SIN(AP1)

! Here is the energy part

      IF (CT.GT.PI) THEN
         CT=CT-2*PI
      ELSEIF (CT.LT.-PI) THEN
         CT=CT+2*PI
      ENDIF

      A1=CT-TWISTREF
      IF (A1.GT.PI) THEN
         A1=A1-2*PI
      ELSEIF (A1.LT.-PI) THEN
         A1=A1+2*PI
      ENDIF

!     Energy=Energy+PK*A1
      Energy=Energy+PK*SIN(A1)
!     PRINT '(A,5G20.10)',' twist> f, energy, de, dphi, angle=',PK,energy,PK*A1,A1,CT
! dE/dPHI
!     DF=PK
      DF=PK*COS(A1)

      IF (.NOT.GTEST) RETURN

! now, do dPhi/dX

! TT1 = r_kj/|D|^2 * dE/dphi
      TT1 = Z1*Z1*RKJ*DF
! TT2 = (r_ij.r_kj / |r_kj||D|^2) * dE/dphi
      TT2 = FGoverG*Z1*Z1*DF
! TT3 = (r_kl.r_kj / |r_kj||G|^2) * dE/dphi
      TT3 = HGoverG*Z2*Z2*DF
! TT4 = r_kj/|G|^2 * dE/dphi
      TT4 = Z2*Z2*RKJ*DF

! note: negatives are flipped from paper because A=-DX
      TT1X=TT1*DX
      TT1Y=TT1*DY
      TT1Z=TT1*DZ

      TT2X=TT2*DX
      TT2Y=TT2*DY
      TT2Z=TT2*DZ


      TT3X=TT3*GX
      TT3Y=TT3*GY
      TT3Z=TT3*GZ


      TT4X=TT4*GX
      TT4Y=TT4*GY
      TT4Z=TT4*GZ

      I3 = IP
      J3 = JP
      K3 = KP
      L3 = LP

      grad(I3*3-2) =  grad(I3*3-2)  + TT1X  
      grad(I3*3-1) =  grad(I3*3-1)  + TT1Y
      grad(I3*3)   =  grad(I3*3)    + TT1Z
      grad(J3*3-2) =  grad(J3*3-2)  - TT1X - TT2X - TT3X
      grad(J3*3-1) =  grad(J3*3-1)  - TT1Y - TT2Y - TT3Y
      grad(J3*3)   =  grad(J3*3)    - TT1Z - TT2Z - TT3Z
      grad(K3*3-2) =  grad(K3*3-2)  + TT2X + TT3X - TT4X
      grad(K3*3-1) =  grad(K3*3-1)  + TT2Y + TT3Y - TT4Y
      grad(K3*3)   =  grad(K3*3)    + TT2Z + TT3Z - TT4Z
      grad(L3*3-2) =  grad(L3*3-2)  + TT4X
      grad(L3*3-1) =  grad(L3*3-1)  + TT4Y
      grad(L3*3)   =  grad(L3*3)    + TT4Z

      RETURN
      END

