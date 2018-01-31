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
C  Subroutine MDIFF calculates the cartesian gradient and second
C  derivative matrix analytically for the Morse potential.
C
C*************************************************************************
C
      SUBROUTINE MPDIFF(N,X,V,P2,RHO,BOXLX,BOXLY,BOXLZ,CUTOFF,GTEST,STEST)
      USE MODHESS
      USE KEY,ONLY : NORESET,FIXIMAGE
      IMPLICIT NONE 
      INTEGER N, J1, J2, J3, J4, J5, J6, I, J, MX, MY, MZ
      DOUBLE PRECISION X(3*N), RHO, 
     1                 V(3*N), DIST, R(N,N), RR(N,N), P2, MEXP(N,N),
     2                 BOXLX, BOXLY, BOXLZ, VEC(N,N,3), 
     3                 CUTOFF, RCUT, RRCUT
      LOGICAL GTEST, STEST

      IF ((BOXLX.LT.BOXLY).AND.(BOXLX.LT.BOXLZ)) THEN
        RCUT=BOXLX*CUTOFF
      ELSE IF (BOXLY.LT.BOXLZ) THEN
        RCUT=BOXLY*CUTOFF
      ELSE
        RCUT=BOXLZ*CUTOFF
      ENDIF
!     PRINT*,'Cutoff used = ',RCUT

C
C  Deal with atoms leaving the box:
C
      IF ((.NOT.FIXIMAGE).AND.(.NOT.NORESET)) THEN
         DO 41 J1=1,N
            X(3*(J1-1)+1)=X(3*(J1-1)+1) - BOXLX*DNINT(X(3*(J1-1)+1)/BOXLX) 
            X(3*(J1-1)+2)=X(3*(J1-1)+2) - BOXLY*DNINT(X(3*(J1-1)+2)/BOXLY)
            X(3*(J1-1)+3)=X(3*(J1-1)+3) - BOXLZ*DNINT(X(3*(J1-1)+3)/BOXLZ)
41       CONTINUE
      ENDIF
C
C  Calculation of connecting vectors; to implement the periodic
C  boundary conditions, the shortest vector between two atoms is
C  used:
C
      DO 25 J1=1,N
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
         DO 15 J2=J1+1,N
            VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
            VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
            VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            MX=NINT(VEC(J2,J1,1)/BOXLX)
            MY=NINT(VEC(J2,J1,2)/BOXLY)
            MZ=NINT(VEC(J2,J1,3)/BOXLZ)
            VEC(J2,J1,1)=VEC(J2,J1,1) - BOXLX * MX
            VEC(J2,J1,2)=VEC(J2,J1,2) - BOXLY * MY
            VEC(J2,J1,3)=VEC(J2,J1,3) - BOXLZ * MZ
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
15       CONTINUE
25    CONTINUE
C
C  Store distance matrices.
C
      P2=0.0D0
      DO 20 J1=1,N
         R(J1,J1)=0.0D0
         RR(J1,J1)=0.0D0
         DO 10 J2=1,J1-1
            DIST=VEC(J1,J2,1)**2 + VEC(J1,J2,2)**2 +
     1           VEC(J1,J2,3)**2
            R(J2,J1)=RHO*DSQRT(DIST)
            R(J1,J2)=R(J2,J1)
            RR(J2,J1)=1.0D0/R(J2,J1)
            RR(J1,J2)=RR(J2,J1)
            MEXP(J1,J2)=DEXP(RHO-R(J1,J2))
            MEXP(J2,J1)=MEXP(J1,J2)
10       CONTINUE
20    CONTINUE
C
C  Here we calculate the energy
C                                        
      RRCUT=RHO*RCUT
      DO 22 I=1,N
         DO 23 J=I+1,N
            IF ((I.NE.J).AND.(R(I,J).LT.RRCUT)) THEN
               P2=P2+MEXP(I,J)*(MEXP(I,J)-2.0D0)
            ENDIF
23       CONTINUE
22    CONTINUE
      IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN
C
C  First calculate the gradient analytically.
C
      DO 50 J1=1,N
         DO 40 J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO 30 J4=1,N
              IF (R(J4,J1).LT.RRCUT) THEN 
                V(J3)=V(J3)-2.0*RHO*RHO*VEC(J1,J4,J2)*
     1               MEXP(J4,J1)* 
     2              (MEXP(J4,J1)-1.0D0)*RR(J4,J1)
              ENDIF 
30          CONTINUE
C           PRINT*,'J3,V(J3)=',J3,V(J3)
40       CONTINUE
50    CONTINUE
      IF (.NOT.STEST) RETURN
C
C  Now do the hessian. First are the entirely diagonal terms.
C
      DO 80 J1=1,N
         DO 70 J2=1,3
            J3=3*(J1-1)+J2
            HESS(J3,J3)=0.0D0
            DO 60 J4=1,N
              IF (R(J4,J1).LT.RRCUT) THEN 
               HESS(J3,J3)=HESS(J3,J3)+2.0*RHO*RHO*MEXP(J4,J1)
     1         *(  (MEXP(J4,J1)-1.0D0) 
     2         * ( (RHO * VEC(J1,J4,J2) *RR(J4,J1) )**2-1.0D0)  
     4         + ( RHO*VEC(J1,J4,J2) )**2*RR(J4,J1)*
     5             (2.0*MEXP(J4,J1)-1.0D0) ) * RR(J4,J1)
              ENDIF 
60          CONTINUE
C           PRINT*,'HESS(',J3,',',J3,')=',HESS(J3,J3)
70       CONTINUE
80    CONTINUE
C
C  Next are the terms where x_i and x_j are on the same atom
C  but are different, e.g. y and z.
C
      DO 120 J1=1,N
         DO 110 J2=1,3
            J3=3*(J1-1)+J2
            DO 100 J4=J2+1,3
               HESS(J3,3*(J1-1)+J4)=0.0D0
               DO 90 J5=1,N
                IF (R(J5,J1).LT.RRCUT) THEN 
                 HESS(J3,3*(J1-1)+J4)=HESS(J3,3*(J1-1)+J4)
     1           +2.0*RHO**4*MEXP(J5,J1)
     2           *VEC(J1,J5,J2)*VEC(J1,J5,J4)
     3           *( (MEXP(J5,J1)-1.0D0)*RR(J5,J1)+
     4              2.0*MEXP(J5,J1)-1.0D0) * RR(J5,J1)**2
                ENDIF 
90             CONTINUE
               HESS(3*(J1-1)+J4,J3)=HESS(J3,3*(J1-1)+J4)
C              PRINT*,'HESS(',J3,',',3*(J1-1)+J4,')=',HESS(J3,3*(J1-1)+J4)
100         CONTINUE
110      CONTINUE
120   CONTINUE
C
C  Case III, different atoms, same cartesian coordinate.
C
      DO 150 J1=1,N
         DO 140 J2=1,3
            J3=3*(J1-1)+J2
            DO 130 J4=J1+1,N
              IF (R(J4,J1).LT.RRCUT) THEN 
               HESS(J3,3*(J4-1)+J2)=2.0*RHO**2*MEXP(J4,J1)
     1      *(  (-RHO**2*VEC(J1,J4,J2)**2*RR(J4,J1)**2+1.0) 
     2           *   (MEXP(J4,J1)-1.0D0)
     3           -   VEC(J1,J4,J2)**2 * RR(J4,J1) * RHO**2
     4           *   (2.0*MEXP(J4,J1)-1.0D0)   )  * RR(J4,J1)
              ENDIF 
C              PRINT*,'HESS(',J3,',',3*(J4-1)+J2,')=',HESS(J3,3*(J4-1)+J2)
               HESS(3*(J4-1)+J2,J3)=HESS(J3,3*(J4-1)+J2)
130         CONTINUE
140      CONTINUE
150   CONTINUE
C
C  Case IV: different atoms and different cartesian coordinates.
C
      DO 180 J1=1,N
         DO 170 J2=1,3
            J3=3*(J1-1)+J2
            DO 160 J4=J1+1,N
               DO 155 J5=1,J2-1
                  J6=3*(J4-1)+J5
                IF (R(J4,J1).LT.RRCUT) THEN
                  HESS(J3,J6)=
     1           -2.0*RHO**4*MEXP(J4,J1)
     2           *VEC(J1,J4,J2)*VEC(J1,J4,J5)
     3           *( (MEXP(J4,J1)-1.0D0)*RR(J4,J1)+
     4              2.0*MEXP(J4,J1)-1.0D0) * RR(J4,J1)**2
                ENDIF
                  HESS(J6,J3)=HESS(J3,J6)
C                 PRINT*,'HESS(',J3,',',J6,')=',HESS(J3,J6)
155            CONTINUE
               DO 156 J5=J2+1,3
                 J6=3*(J4-1)+J5
                IF (R(J4,J1).LT.RRCUT) THEN
                  HESS(J3,J6)=
     1           -2.0*RHO**4*MEXP(J4,J1)
     2           *VEC(J1,J4,J2)*VEC(J1,J4,J5)
     3           *( (MEXP(J4,J1)-1.0D0)*RR(J4,J1)+
     4              2.0*MEXP(J4,J1)-1.0D0) * RR(J4,J1)**2
                ENDIF
                  HESS(J6,J3)=HESS(J3,J6)
C                 PRINT*,'HESS(',J3,',',J6,')=',HESS(J3,J6)
156            CONTINUE
160         CONTINUE
170      CONTINUE
180   CONTINUE
      RETURN
      END
C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2005 David J. Wales
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
C  Subroutine MDIFF calculates the cartesian gradient and second
C  derivative matrix analytically for the Morse potential.
C  This is Dan Sheppard's version with different units for the potential
C  and a shift at the cutoff distance.
C
C*************************************************************************
C
      SUBROUTINE MPDIFFDS(N,X,V,P2,RHO,BOXLX,BOXLY,BOXLZ,CUTOFF,GTEST,STEST)
      USE KEY,ONLY : NORESET,FIXIMAGE
      USE MODHESS
      IMPLICIT NONE 
      INTEGER N, J1, J2, J3, J4, J5, J6, I, J, MX, MY, MZ
      DOUBLE PRECISION X(3*N), RHO, 
     1                 V(3*N), DIST, R(N,N), RR(N,N), P2, MEXP(N,N),
     2                 BOXLX, BOXLY, BOXLZ, VEC(N,N,3), 
     3                 CUTOFF, RCUT, RRCUT, U0, R0, UC
      LOGICAL GTEST, STEST

CGH: added to match Pt parameters
      R0 = 2.8970D0
      U0 = 0.7102D0

      IF ((BOXLX.LT.BOXLY).AND.(BOXLX.LT.BOXLZ)) THEN
        RCUT=BOXLX*CUTOFF
      ELSE IF (BOXLY.LT.BOXLZ) THEN
        RCUT=BOXLY*CUTOFF
      ELSE
        RCUT=BOXLZ*CUTOFF
      ENDIF
!     PRINT*,'Cutoff used = ',RCUT

CGH: cutoff energy shift
      UC = U0*DEXP(RHO-RHO*RCUT/R0)*(DEXP(RHO-RHO*RCUT/R0)-2.0D0)
C
C  Deal with atoms leaving the box:
C
      IF ((.NOT.FIXIMAGE).AND.(.NOT.NORESET)) THEN
         DO 41 J1=1,N
            X(3*(J1-1)+1)=X(3*(J1-1)+1) - BOXLX*DNINT(X(3*(J1-1)+1)/BOXLX) 
            X(3*(J1-1)+2)=X(3*(J1-1)+2) - BOXLY*DNINT(X(3*(J1-1)+2)/BOXLY)
            X(3*(J1-1)+3)=X(3*(J1-1)+3) - BOXLZ*DNINT(X(3*(J1-1)+3)/BOXLZ)
41       CONTINUE
      ENDIF
C
C  Calculation of connecting vectors; to implement the periodic
C  boundary conditions, the shortest vector between two atoms is
C  used:
C
      DO 25 J1=1,N
         VEC(J1,J1,1)=0.0D0
         VEC(J1,J1,2)=0.0D0
         VEC(J1,J1,3)=0.0D0
         DO 15 J2=J1+1,N
            VEC(J2,J1,1)=X(3*(J2-1)+1)-X(3*(J1-1)+1)
            VEC(J2,J1,2)=X(3*(J2-1)+2)-X(3*(J1-1)+2)
            VEC(J2,J1,3)=X(3*(J2-1)+3)-X(3*(J1-1)+3)
            MX=NINT(VEC(J2,J1,1)/BOXLX)
            MY=NINT(VEC(J2,J1,2)/BOXLY)
            MZ=NINT(VEC(J2,J1,3)/BOXLZ)
            VEC(J2,J1,1)=VEC(J2,J1,1) - BOXLX * MX
            VEC(J2,J1,2)=VEC(J2,J1,2) - BOXLY * MY
            VEC(J2,J1,3)=VEC(J2,J1,3) - BOXLZ * MZ
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
15       CONTINUE
25    CONTINUE
C
C  Store distance matrices.
C
      P2=0.0D0
      DO 20 J1=1,N
         R(J1,J1)=0.0D0
         RR(J1,J1)=0.0D0
         DO 10 J2=1,J1-1
            DIST=VEC(J1,J2,1)**2 + VEC(J1,J2,2)**2 +
     1           VEC(J1,J2,3)**2
CGH            R(J2,J1)=RHO*DSQRT(DIST)/R0
            R(J2,J1)=RHO*DSQRT(DIST)/R0
            R(J1,J2)=R(J2,J1)
            RR(J2,J1)=1.0D0/R(J2,J1)
            RR(J1,J2)=RR(J2,J1)
            MEXP(J1,J2)=DEXP(RHO-R(J1,J2))
            MEXP(J2,J1)=MEXP(J1,J2)
10       CONTINUE
20    CONTINUE
C
C  Here we calculate the energy
C
CGH      RRCUT=RHO*RCUT
      RRCUT=RHO*RCUT/R0
      DO 22 I=1,N
         DO 23 J=I+1,N
            IF ((I.NE.J).AND.(R(I,J).LT.RRCUT)) THEN
CGH               P2=P2+MEXP(I,J)*(MEXP(I,J)-2.0D0)
               P2=P2+U0*MEXP(I,J)*(MEXP(I,J)-2.0D0)-UC
            ENDIF
23       CONTINUE
22    CONTINUE
      IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN
C
C  First calculate the gradient analytically.
C
      DO 50 J1=1,N
         DO 40 J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO 30 J4=1,N
              IF (R(J4,J1).LT.RRCUT) THEN 
CGH                V(J3)=V(J3)-2.0*RHO*RHO*VEC(J1,J4,J2)*
                V(J3)=V(J3)-2.0*U0*RHO*RHO*VEC(J1,J4,J2)*
     1               MEXP(J4,J1)* 
     2              (MEXP(J4,J1)-1.0D0)*RR(J4,J1)/(R0*R0)
              ENDIF 
30          CONTINUE
C           PRINT*,'J3,V(J3)=',J3,V(J3)
40       CONTINUE
50    CONTINUE

      IF (.NOT.STEST) RETURN
C
C  Now do the hessian. First are the entirely diagonal terms.
C
      DO 80 J1=1,N
         DO 70 J2=1,3
            J3=3*(J1-1)+J2
            HESS(J3,J3)=0.0D0
            DO 60 J4=1,N
              IF (R(J4,J1).LT.RRCUT) THEN 
               HESS(J3,J3)=HESS(J3,J3)+2.0*RHO*RHO*MEXP(J4,J1)
     1         *(  (MEXP(J4,J1)-1.0D0) 
     2         * ( (RHO * VEC(J1,J4,J2) *RR(J4,J1) )**2-1.0D0)  
     4         + ( RHO*VEC(J1,J4,J2) )**2*RR(J4,J1)*
     5             (2.0*MEXP(J4,J1)-1.0D0) ) * RR(J4,J1)
              ENDIF 
60          CONTINUE
C           PRINT*,'HESS(',J3,',',J3,')=',HESS(J3,J3)
70       CONTINUE
80    CONTINUE
C
C  Next are the terms where x_i and x_j are on the same atom
C  but are different, e.g. y and z.
C
      DO 120 J1=1,N
         DO 110 J2=1,3
            J3=3*(J1-1)+J2
            DO 100 J4=J2+1,3
               HESS(J3,3*(J1-1)+J4)=0.0D0
               DO 90 J5=1,N
                IF (R(J5,J1).LT.RRCUT) THEN 
                 HESS(J3,3*(J1-1)+J4)=HESS(J3,3*(J1-1)+J4)
     1           +2.0*RHO**4*MEXP(J5,J1)
     2           *VEC(J1,J5,J2)*VEC(J1,J5,J4)
     3           *( (MEXP(J5,J1)-1.0D0)*RR(J5,J1)+
     4              2.0*MEXP(J5,J1)-1.0D0) * RR(J5,J1)**2
                ENDIF 
90             CONTINUE
               HESS(3*(J1-1)+J4,J3)=HESS(J3,3*(J1-1)+J4)
C              PRINT*,'HESS(',J3,',',3*(J1-1)+J4,')=',HESS(J3,3*(J1-1)+J4)
100         CONTINUE
110      CONTINUE
120   CONTINUE
C
C  Case III, different atoms, same cartesian coordinate.
C
      DO 150 J1=1,N
         DO 140 J2=1,3
            J3=3*(J1-1)+J2
            DO 130 J4=J1+1,N
              IF (R(J4,J1).LT.RRCUT) THEN 
               HESS(J3,3*(J4-1)+J2)=2.0*RHO**2*MEXP(J4,J1)
     1      *(  (-RHO**2*VEC(J1,J4,J2)**2*RR(J4,J1)**2+1.0) 
     2           *   (MEXP(J4,J1)-1.0D0)
     3           -   VEC(J1,J4,J2)**2 * RR(J4,J1) * RHO**2
     4           *   (2.0*MEXP(J4,J1)-1.0D0)   )  * RR(J4,J1)
              ENDIF 
C              PRINT*,'HESS(',J3,',',3*(J4-1)+J2,')=',HESS(J3,3*(J4-1)+J2)
               HESS(3*(J4-1)+J2,J3)=HESS(J3,3*(J4-1)+J2)
130         CONTINUE
140      CONTINUE
150   CONTINUE
C
C  Case IV: different atoms and different cartesian coordinates.
C
      DO 180 J1=1,N
         DO 170 J2=1,3
            J3=3*(J1-1)+J2
            DO 160 J4=J1+1,N
               DO 155 J5=1,J2-1
                  J6=3*(J4-1)+J5
                IF (R(J4,J1).LT.RRCUT) THEN
                  HESS(J3,J6)=
     1           -2.0*RHO**4*MEXP(J4,J1)
     2           *VEC(J1,J4,J2)*VEC(J1,J4,J5)
     3           *( (MEXP(J4,J1)-1.0D0)*RR(J4,J1)+
     4              2.0*MEXP(J4,J1)-1.0D0) * RR(J4,J1)**2
                ENDIF
                  HESS(J6,J3)=HESS(J3,J6)
C                 PRINT*,'HESS(',J3,',',J6,')=',HESS(J3,J6)
155            CONTINUE
               DO 156 J5=J2+1,3
                 J6=3*(J4-1)+J5
                IF (R(J4,J1).LT.RRCUT) THEN
                  HESS(J3,J6)=
     1           -2.0*RHO**4*MEXP(J4,J1)
     2           *VEC(J1,J4,J2)*VEC(J1,J4,J5)
     3           *( (MEXP(J4,J1)-1.0D0)*RR(J4,J1)+
     4              2.0*MEXP(J4,J1)-1.0D0) * RR(J4,J1)**2
                ENDIF
                  HESS(J6,J3)=HESS(J3,J6)
C                 PRINT*,'HESS(',J3,',',J6,')=',HESS(J3,J6)
156            CONTINUE
160         CONTINUE
170      CONTINUE
180   CONTINUE
      RETURN
      END
