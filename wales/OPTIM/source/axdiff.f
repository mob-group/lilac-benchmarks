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
C*************************************************************************
C
C  Here we calculate the analytic gradient and Hessian for the 
C  Axilrod-Teller term. These are added to the existing derivatives.
C                                        
C*************************************************************************
C
      SUBROUTINE AXDIFF(N, X, V, ZSTAR, GTEST, STEST)
      USE MODHESS
      IMPLICIT NONE 
      INTEGER N, J1, J2, J3, J4, J5
      DOUBLE PRECISION X(3*N), ZSTAR, ABBC, ABAC, ACBC, ABI, ACI, BCI,
     1                 R2(N,N), V(3*N), RX, RY, RZ, 
     2                 VEC(N,N,3), TEMP, RR2(N,N), TDOT,
     3                 RAB, RRAB, RAC, RRAC, RBC, RRBC,
     4                 ABJ, ACJ, BCJ
      LOGICAL GTEST, STEST
      DO 20 J1=1,N
         R2(J1,J1)=0.0D0
         RR2(J1,J1)=0.0D0
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
         DO 10 J2=J1+1,N
            R2(J2,J1)=
     1                ( X(3*(J2-1)+1)-X(3*(J1-1)+1) )**2 +
     2                ( X(3*(J2-1)+2)-X(3*(J1-1)+2) )**2 +
     3                ( X(3*(J2-1)+3)-X(3*(J1-1)+3) )**2 
            RR2(J2,J1)=1.0D0/R2(J2,J1)
            VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
            VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
            VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            R2(J1,J2)=R2(J2,J1) 
            RR2(J1,J2)=RR2(J2,J1) 
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
10       CONTINUE 
20    CONTINUE 
      IF (GTEST) THEN
C
C  First the gradient.
C
      DO 120 J1=1,N
         DO 110 J2=1,3
            TEMP=0.0D0
            DO 100 J3=1,N
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               ABI=VEC(J3,J1,J2)
               RX=VEC(J3,J1,1)
               RY=VEC(J3,J1,2)
               RZ=VEC(J3,J1,3)
               DO 95 J4=J3+1,N
                  ABAC=RX*VEC(J4,J1,1)+RY*VEC(J4,J1,2)+RZ*VEC(J4,J1,3)
                  ABBC=RX*VEC(J4,J3,1)+RY*VEC(J4,J3,2)+RZ*VEC(J4,J3,3)
                  ACBC=VEC(J4,J1,1)*VEC(J4,J3,1)
     1                +VEC(J4,J1,2)*VEC(J4,J3,2)
     2                +VEC(J4,J1,3)*VEC(J4,J3,3)
                  TDOT=ABAC*ACBC*ABBC
                  BCI=VEC(J4,J3,J2)
                  ACI=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1)
C
                  TEMP=TEMP+ DSQRT(RRAB*RRAC*RRBC)**5 *          (
     1    3*(ABAC*ACBC*BCI + ABBC*(ACBC*(ABI + ACI) + ABAC*BCI) + 
     2                            (ACI*RAB + ABI*RAC)*RBC) - 
     3                         15*(ABI*RRAB + ACI*RRAC)*TDOT   )

C
95             CONTINUE
100         CONTINUE
            V(3*(J1-1)+J2)=ZSTAR*TEMP
     1                     + V(3*(J1-1)+J2)
C           PRINT*,'K2,V=',3*(J1-1)+J2,V(3*(J1-1)+J2)
110      CONTINUE
120   CONTINUE
      ENDIF
      
      IF (STEST) THEN
C
C  Diagonal bits of the Hessian.
C
      DO 160 J1=1,N
         DO 150 J2=1,3
            TEMP=0.0D0
            DO 140 J3=1,N
               RAB=R2(J3,J1)
               RRAB=RR2(J3,J1)
               RX=VEC(J3,J1,1)
               RY=VEC(J3,J1,2)
               RZ=VEC(J3,J1,3)
               ABI=VEC(J3,J1,J2)
               DO 130 J4=J3+1,N
                  BCI=VEC(J4,J3,J2)
                  ACI=VEC(J4,J1,J2)
                  RBC=R2(J4,J3)
                  RRBC=RR2(J4,J3)
                  RAC=R2(J4,J1)
                  RRAC=RR2(J4,J1)
                  ABAC=RX*VEC(J4,J1,1)+RY*VEC(J4,J1,2)+RZ*VEC(J4,J1,3)
                  ABBC=RX*VEC(J4,J3,1)+RY*VEC(J4,J3,2)+RZ*VEC(J4,J3,3)
                  ACBC=VEC(J4,J1,1)*VEC(J4,J3,1)
     1                +VEC(J4,J1,2)*VEC(J4,J3,2)
     2                +VEC(J4,J1,3)*VEC(J4,J3,3)
                  TDOT=ABAC*ACBC*ABBC
C
                  TEMP=TEMP + DSQRT(RRAB*RRAC*RRBC)**5 * (
     1    -6*((ABBC*ACI + ACBC*(ABI + ACI))*BCI + ABAC*BCI**2 + 
     2     ABBC*(ACBC + ABI*BCI)) + 
     3  30*(ABBC*(ABAC*BCI*(ABI*RRAB + ACI*RRAC) + 
     4        ACBC*((ABI**2 + ABI*ACI)*RRAB + ABI*ACI*RRAC)) + 
     5     ACBC*(ABAC*BCI*(ABI*RRAB + ACI*RRAC) + ABBC*ACI**2*RRAC)) +  
     6  (18*ABI*ACI + RAB*(-3.0 + 15*ACI**2*RRAC) + 
     7     (-3. + 15*ABI**2*RRAB)*RAC)*RBC + 
     8  (-105*(ABI**2*RRAB**2 + ACI**2*RRAC**2) + 15*(RRAB + RRAC) - 
     9     150*ABI*ACI*RRAB*RRAC)*TDOT               )
C
130            CONTINUE
140         CONTINUE
            HESS(3*(J1-1)+J2,3*(J1-1)+J2)=ZSTAR*TEMP 
     1                               + HESS(3*(J1-1)+J2,3*(J1-1)+J2) 
C           PRINT*,'K2,A=',3*(J1-1)+J2,HESS(3*(J1-1)+J2,3*(J1-1)+J2)
150      CONTINUE
160   CONTINUE
C
C  Same atom, different component.
C
      DO 210 J1=1,N
        DO 200 J2=1,3
          DO 190 J5=J2+1,3 
            TEMP=0.0D0
            DO 180 J3=1,N
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              RX=VEC(J3,J1,1)
              RY=VEC(J3,J1,2)
              RZ=VEC(J3,J1,3)
              ABI=VEC(J3,J1,J2)
              ABJ=VEC(J3,J1,J5)
              DO 170 J4=J3+1,N
                BCI=VEC(J4,J3,J2)
                ACI=VEC(J4,J1,J2)
                BCJ=VEC(J4,J3,J5)
                ACJ=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                ABAC=RX*VEC(J4,J1,1)+RY*VEC(J4,J1,2)+RZ*VEC(J4,J1,3)
                ABBC=RX*VEC(J4,J3,1)+RY*VEC(J4,J3,2)+RZ*VEC(J4,J3,3) 
                ACBC=VEC(J4,J1,1)*VEC(J4,J3,1)
     1              +VEC(J4,J1,2)*VEC(J4,J3,2)
     2              +VEC(J4,J1,3)*VEC(J4,J3,3)
                TDOT=ABAC*ACBC*ABBC

                TEMP=TEMP + DSQRT(RRAB*RRAC*RRBC)**5 *              (
     1   -3*(ABBC*ACI*BCJ + (ABBC + ACBC)*((ABJ + ACJ)*BCI + ABI*BCJ)) + 
     2  BCJ*(ACBC*(15*ABAC*ABI*RRAB + ACI*(-3.0 + 15*ABAC*RRAC)) + 
     3     ABAC*(-6*BCI + 15*ABBC*(ABI*RRAB + ACI*RRAC))) + 
     4  ABBC*(ABJ*(15*ABAC*BCI*RRAB + 
     5        ACBC*((30*ABI + 15*ACI)*RRAB + 15*ACI*RRAC)) + 
     6     ACBC*ACJ*(15*ABI*RRAB + 30*ACI*RRAC)) + 
     7  15*(ACBC*(ABAC*BCI*(ABJ*RRAB + ACJ*RRAC) + ABBC*ABI*ACJ*RRAC) + 
     8     ABI*ABJ*RAC*RBC*RRAB + 
     9     ACJ*(ABAC*ABBC*BCI + ACI*RAB*RBC)*RRAC) - 
     A  105*(ABI*ABJ*RRAB**2 + ACI*ACJ*RRAC**2)*TDOT + 
     B  (ABJ*ACI + ABI*ACJ)*(9*RBC - 75*TDOT*RRAB*RRAC)          )

170            CONTINUE
180         CONTINUE
            HESS(3*(J1-1)+J5,3*(J1-1)+J2)=ZSTAR*TEMP
     1                                 + HESS(3*(J1-1)+J5,3*(J1-1)+J2)
190      CONTINUE
200   CONTINUE
210   CONTINUE
C
C  Different atoms, same component.
C
      DO 260 J1=1,N
        DO 250 J2=1,3
          DO 230 J3=J1+1,N
            RAB=R2(J3,J1)
            RRAB=RR2(J3,J1)
            RX=VEC(J3,J1,1)
            RY=VEC(J3,J1,2)
            RZ=VEC(J3,J1,3)
            ABI=VEC(J3,J1,J2)
            TEMP=0.0D0
            DO 220 J4=1,N
              BCI=VEC(J4,J3,J2)
              ACI=VEC(J4,J1,J2)
              RBC=R2(J4,J3)
              RRBC=RR2(J4,J3)
              RAC=R2(J4,J1)
              RRAC=RR2(J4,J1)
              ABAC=RX*VEC(J4,J1,1)+RY*VEC(J4,J1,2)+RZ*VEC(J4,J1,3)
              ABBC=RX*VEC(J4,J3,1)+RY*VEC(J4,J3,2)+RZ*VEC(J4,J3,3) 
              ACBC=VEC(J4,J1,1)*VEC(J4,J3,1)
     1            +VEC(J4,J1,2)*VEC(J4,J3,2)
     2            +VEC(J4,J1,3)*VEC(J4,J3,3)
              TDOT=ABAC*ACBC*ABBC

              TEMP=TEMP + DSQRT(RRAB*RRAC*RRBC)**5 *            ( 
     1  -3*(ABBC*(ABI*ACI + ACI**2) + ABAC*(ABBC + ACBC + ABI*BCI)) + 
     2  ABAC*((15*ABI**2*ACBC - 15*ABBC*ABI*BCI)*RRAB + 
     3     15*ABBC*ACI**2*RRAC + BCI**2*(3.0 + 15*ABBC*RRBC)) + 
     4  ACBC*(-3*(ABI**2 + ABI*ACI) + 3*(ABBC + ABI*BCI) + 
     5     (-15*ABBC*ABI**2 - 30*ABAC*ABI*BCI)*RRAB + 
     6     (15*ABAC*ABI*ACI - 15*ABBC*ACI**2)*RRAC + 
     7     15*(ABBC*ABI*BCI + ABAC*BCI**2)*RRBC) + 
     8  (3 - 15*ABI**2*RRAB)*RAC*RBC + 
     9  (105*ABI**2*RRAB**2 - 15*RRAB - 75*ACI*BCI*RRAC*RRBC)*
     A   TDOT + ACI*(BCI*(3*ABBC + 6*ACBC + 9*RAB + 
     B        ABAC*(-3.0 - 15*ACBC*RRAC)) + 
     C     ABBC*(ABI*(15*ABAC - 30*ACBC)*RRAB + 15*ACBC*BCI*RRBC) + 
     D     75*ABI*TDOT*RRAB*RRAC) + 
     E  ABI*(-9*ACI*RBC + BCI*(9*RAC - 75*TDOT*RRAB*RRBC))  )

220         CONTINUE
            HESS(3*(J3-1)+J2,3*(J1-1)+J2)=ZSTAR*TEMP
     1                                 + HESS(3*(J3-1)+J2,3*(J1-1)+J2)
230       CONTINUE
250     CONTINUE
260   CONTINUE
C
C  Different atoms and different components
C
      DO 310 J1=1,N
        DO 300 J2=1,3
          DO 280 J3=J1+1,N
            DO 290 J5=1,J2-1 
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              RX=VEC(J3,J1,1)
              RY=VEC(J3,J1,2)
              RZ=VEC(J3,J1,3)
              ABI=VEC(J3,J1,J2)
              ABJ=VEC(J3,J1,J5)
              TEMP=0.0D0
              DO 270 J4=1,N
                BCI=VEC(J4,J3,J2)
                ACI=VEC(J4,J1,J2)
                BCJ=VEC(J4,J3,J5)
                ACJ=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                ABAC=RX*VEC(J4,J1,1)+RY*VEC(J4,J1,2)+RZ*VEC(J4,J1,3)
                ABBC=RX*VEC(J4,J3,1)+RY*VEC(J4,J3,2)+RZ*VEC(J4,J3,3) 
                ACBC=VEC(J4,J1,1)*VEC(J4,J3,1)
     1              +VEC(J4,J1,2)*VEC(J4,J3,2)
     2              +VEC(J4,J1,3)*VEC(J4,J3,3)
                TDOT=ABAC*ACBC*ABBC
                TEMP=TEMP + DSQRT(RRAB*RRAC*RRBC)**5 * (
     1  -3*(ABBC*(ABI + ACI)*ACJ + ABJ*(ACBC*(ABI + ACI) + ABAC*BCI)) + 
     2  ACJ*((-3*ABAC + 3*(ABBC + ACBC))*BCI + 
     3     ABBC*(ABI*(15*ABAC - 15*ACBC)*RRAB + 15*ABAC*ACI*RRAC)) - 
     4  15*(ABBC*(ABJ*(ACBC*ACI + ABAC*BCI)*RRAB + ACBC*ACI*ACJ*RRAC) + 
     5     ABJ*(ABAC*ACBC*BCI + ABI*RAC*RBC)*RRAB) + 
     6  ABI*(BCJ*(9*RAC + ACBC*(-15*ABAC*RRAB + 15*ABBC*RRBC)) + 
     7     ABJ*((15*ABAC - 15*ABBC)*ACBC*RRAB + 105*TDOT*RRAB**2)) + 
     8  BCJ*(3*(ACBC*(ABI + ACI) + ABAC*BCI) + 
     9     (15*ABAC*(ABBC + ACBC)*BCI - 
     A        75*(ABI*RRAB + ACI*RRAC)*TDOT)*RRBC) + 
     B  ACI*(BCJ*(9*RAB + ACBC*(-15*ABAC*RRAC + 15*ABBC*RRBC)) + 
     C     ABJ*(-9*RBC + (15*ABAC*ACBC + 75*TDOT*RRAB)*RRAC)) )

270           CONTINUE
              HESS(3*(J3-1)+J5,3*(J1-1)+J2)=ZSTAR*TEMP
     1                                   + HESS(3*(J3-1)+J5,3*(J1-1)+J2)
290         CONTINUE
            DO 295 J5=J2+1,3
              RAB=R2(J3,J1)
              RRAB=RR2(J3,J1)
              RX=VEC(J3,J1,1)
              RY=VEC(J3,J1,2)
              RZ=VEC(J3,J1,3)
              ABI=VEC(J3,J1,J2)
              ABJ=VEC(J3,J1,J5)
              TEMP=0.0D0
              DO 275 J4=1,N
                BCI=VEC(J4,J3,J2)
                ACI=VEC(J4,J1,J2)
                BCJ=VEC(J4,J3,J5)
                ACJ=VEC(J4,J1,J5)
                RBC=R2(J4,J3)
                RRBC=RR2(J4,J3)
                RAC=R2(J4,J1)
                RRAC=RR2(J4,J1)
                ABAC=RX*VEC(J4,J1,1)+RY*VEC(J4,J1,2)+RZ*VEC(J4,J1,3)
                ABBC=RX*VEC(J4,J3,1)+RY*VEC(J4,J3,2)+RZ*VEC(J4,J3,3) 
                ACBC=VEC(J4,J1,1)*VEC(J4,J3,1)
     1              +VEC(J4,J1,2)*VEC(J4,J3,2)
     2              +VEC(J4,J1,3)*VEC(J4,J3,3)
                TDOT=ABAC*ACBC*ABBC
                TEMP=TEMP + DSQRT(RRAB*RRAC*RRBC)**5 * (
     1  -3*(ABBC*(ABI + ACI)*ACJ + ABJ*(ACBC*(ABI + ACI) + ABAC*BCI)) + 
     2  ACJ*((-3*ABAC + 3*(ABBC + ACBC))*BCI + 
     3     ABBC*(ABI*(15*ABAC - 15*ACBC)*RRAB + 15*ABAC*ACI*RRAC)) - 
     4  15*(ABBC*(ABJ*(ACBC*ACI + ABAC*BCI)*RRAB + ACBC*ACI*ACJ*RRAC) + 
     5     ABJ*(ABAC*ACBC*BCI + ABI*RAC*RBC)*RRAB) + 
     6  ABI*(BCJ*(9*RAC + ACBC*(-15*ABAC*RRAB + 15*ABBC*RRBC)) + 
     7     ABJ*((15*ABAC - 15*ABBC)*ACBC*RRAB + 105*TDOT*RRAB**2)) + 
     8  BCJ*(3*(ACBC*(ABI + ACI) + ABAC*BCI) + 
     9     (15*ABAC*(ABBC + ACBC)*BCI - 
     A        75*(ABI*RRAB + ACI*RRAC)*TDOT)*RRBC) + 
     B  ACI*(BCJ*(9*RAB + ACBC*(-15*ABAC*RRAC + 15*ABBC*RRBC)) + 
     C     ABJ*(-9*RBC + (15*ABAC*ACBC + 75*TDOT*RRAB)*RRAC)) )

275           CONTINUE
              HESS(3*(J3-1)+J5,3*(J1-1)+J2)=ZSTAR*TEMP
     1                                   + HESS(3*(J3-1)+J5,3*(J1-1)+J2)
295         CONTINUE
280       CONTINUE
300     CONTINUE
310   CONTINUE
C
C  Symmetrise
C
      DO 1000 J1=1,3*N
         DO 1010 J2=J1+1,3*N
C           PRINT*,'J1,J2,A=',J1,J2,HESS(J2,J1)
            HESS(J1,J2)=HESS(J2,J1)
1010     CONTINUE
1000  CONTINUE
      ENDIF

      RETURN
      END
