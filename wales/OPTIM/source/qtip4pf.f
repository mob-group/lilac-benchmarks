C
C*************************************************************************
C
C  Subroutine qTIP4PF calculates the energy and cartesian gradient analytically. 
C  This corresponds to a flexible water model introduced 
C  by Habershon et al. (JCP 131, 024501 (2009))
C
C*************************************************************************
C
      SUBROUTINE qtip4pf(nmol,X,V,ENERGY,GTEST)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6, J7, J8, J9, J10, nmol
      LOGICAL GTEST
C      INCLUDE 'param_atoms.h'
      DOUBLE PRECISION X(9*nmol), V(9*nmol), ENERGY, epsilon, sigma, qm,
     1                 gamma, Dr, alphar, req, ktheta, thetaeq, r6, r12, DIST, 
     2                 ENERGY_LJ, ENERGY_Q, Q(4), HTOKCAL, rmx, rmy, rmz,
     3                 FLJ, DFLJ, FC, DFC, C12, C6, Q1, Q2, ENERGY_OH, ENERGY_THETA,
     4                 XDUMM, FOH, DFOH, FTHETA, DFTHETA, X0(12*nmol), THETA,
     5                 OH1VECX,OH1VECY,OH1VECZ,OH2VECX,OH2VECY,OH2VECZ,COSTHETA,
     6                 RDIST, DUMMY, DRO(4,4),DRH1(4,4),DRH2(4,4),RDIST1,RDIST2,
     7                 DtxH1,DtyH1,DtzH1,DtxH2,DtyH2,DtzH2,dotnorma,
     8                 Gtx1,Gty1,Gtz1,Gtx2,Gty2,Gtz2

      PARAMETER (HTOKCAL=331.952261d0)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C     Statement functions
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      FLJ(XDUMM)=(C12/XDUMM**6-C6)/XDUMM**6
      DFLJ(XDUMM)=-6.0D0*(-C6+2.0D0*C12/XDUMM**6)/XDUMM**7
      FC(Q1,Q2,XDUMM)=Q1*Q2/XDUMM
      DFC(Q1,Q2,XDUMM)=-Q1*Q2/XDUMM**2
      FOH(XDUMM)=Dr*(alphar**2*(XDUMM-req)**2-alphar**3*(XDUMM-req)**3+
     &               7.d0*alphar**4*(XDUMM-req)**4/12.d0)
      DFOH(XDUMM)=Dr*(2.d0*alphar**2*(XDUMM-req)-3.d0*alphar**3*(XDUMM-req)**2+
     &               7.d0*alphar**4*(XDUMM-req)**3/3.d0)
      FTHETA(XDUMM)=0.5d0*ktheta*(XDUMM-thetaeq)**2
      DFTHETA(XDUMM)=1.d0*ktheta*(XDUMM-thetaeq)

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
C    Define the parameters of the potential
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      epsilon=0.1852d0
      sigma=3.1589d0
      qm=1.1128d0
      gamma=0.73612d0
      Dr=116.09d0
      alphar=2.287d0
      req=0.9419d0
      ktheta=87.85d0
      thetaeq=107.4d0*dacos(-1.d0)/180.d0

      C6=4.d0*epsilon*(sigma)**6
      C12=C6*(sigma)**6
      Q(1)=0.d0
      Q(2)=qm/2.d0
      Q(3)=qm/2.d0
      Q(4)=-qm

      N=nmol

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C      Define the M point for all water molecules
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       
      DO J1=1,N
       J2=9*(J1-1)
       J3=12*(J1-1)
       X0(J3+1)=X(J2+1)
       X0(J3+2)=X(J2+2)
       X0(J3+3)=X(J2+3)
       X0(J3+4)=X(J2+4)
       X0(J3+5)=X(J2+5)
       X0(J3+6)=X(J2+6)
       X0(J3+7)=X(J2+7)
       X0(J3+8)=X(J2+8)
       X0(J3+9)=X(J2+9)
       rmx=gamma*X(J2+1)+(1.d0-gamma)*(X(J2+4)+X(J2+7))/2.d0
       rmy=gamma*X(J2+2)+(1.d0-gamma)*(X(J2+5)+X(J2+8))/2.d0  ! cartesian coordinates of the M point
       rmz=gamma*X(J2+3)+(1.d0-gamma)*(X(J2+6)+X(J2+9))/2.d0
       X0(J3+10)=rmx
       X0(J3+11)=rmy
       X0(J3+12)=rmz
      END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      DO J1=1,9*N
       V(J1)=0.d0         ! Three points per molecule
      END DO

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       Intermolecular potential
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      ENERGY_LJ=0.0D0
      ENERGY_Q=0.0D0

      DO J1=1,4
       DO J2=1,4
        DRO(J1,J2)=0.d0
        DRH1(J1,J2)=0.d0
        DRH2(J1,J2)=0.d0
       END DO
      END DO

      DO J2=1,4
       DRO(1,J2)=0.d0
       DRH1(1,J2)=0.d0
       DRH2(1,J2)=0.d0
       DRO(2,J2)=0.d0
       DRH1(2,J2)=1.d0
       DRH2(2,J2)=0.d0
       DRO(3,J2)=0.d0
       DRH1(3,J2)=0.d0
       DRH2(3,J2)=1.d0
       DRO(4,J2)=gamma
       DRH1(4,J2)=0.5d0*(1.d0-gamma)
       DRH2(4,J2)=0.5d0*(1.d0-gamma)
      END DO

      DO J1=1,N
        J3=12*(J1-1)
        J9=9*(J1-1)
        DO J2=1,N
          IF (J1.NE.J2) THEN
          J4=12*(J2-1)
          J10=9*(J1-1)
           DIST=(X0(J3+1)-X0(J4+1))**2
     1         +(X0(J3+2)-X0(J4+2))**2
     2         +(X0(J3+3)-X0(J4+3))**2
           DIST=DSQRT(DIST)
C           IF (DIST.LT.1.8d0) THEN
C             ENERGY=1.d10
C             RETURN
C           END IF
           ENERGY_LJ=ENERGY_LJ+FLJ(DIST)   ! LJ contribution	

           IF (GTEST) THEN
            RDIST=1.d0/DIST
            V(J9+1)=V(J9+1)+DFLJ(DIST)*RDIST*(X0(J3+1)-X0(J4+1))
            V(J9+2)=V(J9+2)+DFLJ(DIST)*RDIST*(X0(J3+2)-X0(J4+2))
            V(J9+3)=V(J9+3)+DFLJ(DIST)*RDIST*(X0(J3+3)-X0(J4+3))
           END IF

           DO J5=2,4
            J7=J3+3*J5-3
            DO J6=2,4
             J8=J4+3*J6-3
             DUMMY=(X0(J7+1)-X0(J8+1))**2
     1            +(X0(J7+2)-X0(J8+2))**2
     2            +(X0(J7+3)-X0(J8+3))**2
             DUMMY=DSQRT(DUMMY)
             ENERGY_Q=ENERGY_Q+HTOKCAL*FC(Q(J5),Q(J6),DUMMY)  ! Coulomb contribution

            IF (GTEST) THEN
              RDIST=1.d0/DUMMY
              V(J9+1)=V(J9+1)+HTOKCAL*DFC(Q(J5),Q(J6),DUMMY)*RDIST*(X0(J7+1)-X0(J8+1))*DRO(J5,J6)
              V(J9+2)=V(J9+2)+HTOKCAL*DFC(Q(J5),Q(J6),DUMMY)*RDIST*(X0(J7+2)-X0(J8+2))*DRO(J5,J6)
              V(J9+3)=V(J9+3)+HTOKCAL*DFC(Q(J5),Q(J6),DUMMY)*RDIST*(X0(J7+3)-X0(J8+3))*DRO(J5,J6)
              V(J9+4)=V(J9+4)+HTOKCAL*DFC(Q(J5),Q(J6),DUMMY)*RDIST*(X0(J7+1)-X0(J8+1))*DRH1(J5,J6)
              V(J9+5)=V(J9+5)+HTOKCAL*DFC(Q(J5),Q(J6),DUMMY)*RDIST*(X0(J7+2)-X0(J8+2))*DRH1(J5,J6)
              V(J9+6)=V(J9+6)+HTOKCAL*DFC(Q(J5),Q(J6),DUMMY)*RDIST*(X0(J7+3)-X0(J8+3))*DRH1(J5,J6)
              V(J9+7)=V(J9+7)+HTOKCAL*DFC(Q(J5),Q(J6),DUMMY)*RDIST*(X0(J7+1)-X0(J8+1))*DRH2(J5,J6)
              V(J9+8)=V(J9+8)+HTOKCAL*DFC(Q(J5),Q(J6),DUMMY)*RDIST*(X0(J7+2)-X0(J8+2))*DRH2(J5,J6)
              V(J9+9)=V(J9+9)+HTOKCAL*DFC(Q(J5),Q(J6),DUMMY)*RDIST*(X0(J7+3)-X0(J8+3))*DRH2(J5,J6)
             END IF

            END DO
           END DO
         END IF
        ENDDO 
      ENDDO 


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C       Intramolecular potential
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      ENERGY_OH=0.0D0
      ENERGY_THETA=0.0D0

      DO J1=1,N
       J2=9*(J1-1)
       OH1VECX=X(J2+4)-X(J2+1)
       OH1VECY=X(J2+5)-X(J2+2)
       OH1VECZ=X(J2+6)-X(J2+3)
       DIST=OH1VECX**2+OH1VECY**2+OH1VECZ**2
       DIST=DSQRT(DIST)                        !OH1 DISTANCE
       RDIST1=1.d0/DIST

       OH2VECX=X(J2+7)-X(J2+1)
       OH2VECY=X(J2+8)-X(J2+2)
       OH2VECZ=X(J2+9)-X(J2+3)
       DUMMY=OH2VECX**2+OH2VECY**2+OH2VECZ**2
       DUMMY=DSQRT(DUMMY)                     !OH2 DISTANCE
       RDIST2=1.d0/DUMMY

       ENERGY_OH=ENERGY_OH+FOH(DIST)+FOH(DUMMY)

       COSTHETA=(OH1VECX*OH2VECX+OH1VECY*OH2VECY+OH1VECZ*OH2VECZ)/(DIST*DUMMY)
       dotnorma=COSTHETA
       THETA=DACOS(COSTHETA)
       ENERGY_THETA=ENERGY_THETA+FTHETA(THETA)   ! HOH 

       IF (GTEST) THEN
        DtxH1=-(1.d0/DSQRT(1.d0-dotnorma*dotnorma))*
     &          (OH2VECX/(DIST*DUMMY)-dotnorma*OH1VECX/DIST**2)
        DtyH1=-(1.d0/DSQRT(1.d0-dotnorma*dotnorma))*
     &          (OH2VECY/(DIST*DUMMY)-dotnorma*OH1VECY/DIST**2)
        DtzH1=-(1.d0/DSQRT(1.d0-dotnorma*dotnorma))*
     &          (OH2VECZ/(DIST*DUMMY)-dotnorma*OH1VECZ/DIST**2)
        DtxH2=-(1.d0/DSQRT(1.d0-dotnorma*dotnorma))*
     &          (OH1VECX/(DIST*DUMMY)-dotnorma*OH2VECX/DUMMY**2)
        DtyH2=-(1.d0/DSQRT(1.d0-dotnorma*dotnorma))*
     &          (OH1VECY/(DIST*DUMMY)-dotnorma*OH2VECY/DUMMY**2)
        DtzH2=-(1.d0/DSQRT(1.d0-dotnorma*dotnorma))*
     &          (OH1VECZ/(DIST*DUMMY)-dotnorma*OH2VECZ/DUMMY**2)

        Gtx1=DFOH(DIST)*RDIST1*OH1VECX+DFTHETA(THETA)*DtxH1
        Gty1=DFOH(DIST)*RDIST1*OH1VECY+DFTHETA(THETA)*DtyH1
        Gtz1=DFOH(DIST)*RDIST1*OH1VECZ+DFTHETA(THETA)*DtzH1
        Gtx2=DFOH(DUMMY)*RDIST2*OH2VECX+DFTHETA(THETA)*DtxH2
        Gty2=DFOH(DUMMY)*RDIST2*OH2VECY+DFTHETA(THETA)*DtyH2
        Gtz2=DFOH(DUMMY)*RDIST2*OH2VECZ+DFTHETA(THETA)*DtzH2
        V(J2+4)=V(J2+4)+Gtx1
        V(J2+5)=V(J2+5)+Gty1
        V(J2+6)=V(J2+6)+Gtz1
        V(J2+7)=V(J2+7)+Gtx2
        V(J2+8)=V(J2+8)+Gty2
        V(J2+9)=V(J2+9)+Gtz2
        V(J2+1)=V(J2+1)-Gtx1-Gtx2
        V(J2+2)=V(J2+2)-Gty1-Gty2
        V(J2+3)=V(J2+3)-Gtz1-Gtz2
       END IF

      END DO


      ENERGY=0.5d0*ENERGY_LJ+0.5d0*ENERGY_Q+ENERGY_OH+ENERGY_THETA

C      WRITE (*,*) 'LJ ',ENERGY_LJ
c      WRITE (*,*) 'Q ',ENERGY_Q
c      WRITE (*,*) 'OH ',ENERGY_OH
C      WRITE (*,*) 'THETA ',ENERGY_THETA
       

      RETURN
      END
