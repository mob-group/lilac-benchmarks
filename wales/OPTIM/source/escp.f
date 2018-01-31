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
C********************************************************************
C
C Subroutine ENERGY calculates the SC energy:
C
C********************************************************************
C
      SUBROUTINE OESCP(N,X,PSC,VNEW,BOXLX,BOXLY,BOXLZ,CUTOFF,PRESSURE,GTEST,STEST)
      IMPLICIT NONE
      INTEGER N, J1, J2, J, I, NN, MM, MZ, MY, MX 
      LOGICAL PRESSURE, YESNO
      DOUBLE PRECISION X(3*N), POTA, POTB, DIST, SIG, C, EPS, CUTOFF, BOXLZ,
     1                 BOXLY, BOXLX, PSC, VNEW(3*N)
      COMMON /POWERS/ NN, MM
      DOUBLE PRECISION RHO(3*N)
      COMMON /INT1/ EPS, C, SIG
      DOUBLE PRECISION RM(N,N), RN(N,N),VEC(N,N,3)
      LOGICAL GTEST,STEST

      INQUIRE(FILE='SCparams',EXIST=YESNO)
      IF (.NOT.YESNO) THEN
C        PRINT*,'No SCparams file - 12,6 assumed'
         NN=12
         MM=6
         EPS=2.5415D-03
         C=144.41
         SIG=1.414
C
C traditional 10,8 params
C
C        EPS=1.2793D-02
C        C=34.408
C        SIG=1.414
      ELSE
         OPEN(UNIT=33,FILE='SCparams',STATUS='OLD')
         READ(33,*) NN,MM,EPS,C,SIG
C        WRITE(*,'(A,I4,A,I4,A,F15.5,A,F15.5,A,F15.5)')  ' n=',NN,' m=',MM,' epsilon=',EPS,' c=',C,' sigma=',SIG
         CLOSE(33)
      ENDIF
C
C  Deal with atoms leaving the box:
C
      DO 21 J1=1,N
         X(3*(J1-1)+1)=X(3*(J1-1)+1) - BOXLX *
     1                  DNINT(X(3*(J1-1)+1)/BOXLX)
         X(3*(J1-1)+2)=X(3*(J1-1)+2) -  BOXLY *
     1                  DNINT(X(3*(J1-1)+2)/BOXLY)
         X(3*(J1-1)+3)=X(3*(J1-1)+3) -  BOXLZ *
     1                  DNINT(X(3*(J1-1)+3)/BOXLZ)
21    CONTINUE
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
            VEC(J2,J1,1)=VEC(J2,J1,1) - BOXLX * FLOAT(MX)
            VEC(J2,J1,2)=VEC(J2,J1,2) - BOXLY * FLOAT(MY)
            VEC(J2,J1,3)=VEC(J2,J1,3) - BOXLZ * FLOAT(MZ)
            VEC(J1,J2,1)=-VEC(J2,J1,1)
            VEC(J1,J2,2)=-VEC(J2,J1,2)
            VEC(J1,J2,3)=-VEC(J2,J1,3)
15       CONTINUE
25    CONTINUE
C
C  Store distance matrices.
C
      DO 20 J1=1,N
         RM(J1,J1)=0.0D0
         RN(J1,J1)=0.0D0
         DO 10 J2=1,J1-1
            DIST=VEC(J1,J2,1)**2 + VEC(J1,J2,2)**2 + VEC(J1,J2,3)**2
            IF (DSQRT(DIST).LT.CUTOFF) THEN
               RM(J2,J1)=DIST**(-FLOAT(MM)/2.0D0)
               RN(J2,J1)=DIST**(-FLOAT(NN)/2.0D0)
               RM(J1,J2)=RM(J2,J1)
               RN(J1,J2)=RN(J2,J1)
            ELSE
               RM(J2,J1)=0.0D0
               RN(J2,J1)=0.0D0
               RM(J1,J2)=RM(J2,J1)
               RN(J1,J2)=RN(J2,J1)
            ENDIF
10       CONTINUE
20    CONTINUE
C
C Call scl.f for lattice constant optimisation if required:
C
      IF (PRESSURE) THEN
         CALL SCL(N,X,EPS,C,SIG,BOXLX,BOXLY,BOXLZ,CUTOFF)
         PRINT*,'Energy minimised with respect to lattice constants' 
C        CUTOFF=BOXLX/2.0D0
         PRINT*,'New box length and cutoff=',BOXLX,CUTOFF
      ENDIF
C
C Store density matrix: In the case of the perfect fcc lattice,
C the infinitely extended crystal implies that every RHO(J) is
C equal to RHO(1).
C
      DO 11 I=1,N
         RHO(I)=0.0D00
         DO 122 J=1,N
            RHO(I)=RHO(I) + SIG**MM*RM(I,J)
122      CONTINUE
11    CONTINUE
C
C calculate the potential energy:
C
      POTA=0.0D0
      POTB=0.0D0
      DO 13 I=1,N
         DO 14 J=1,N
            POTA=POTA + 0.50D00*EPS*SIG**NN*RN(I,J)
14       CONTINUE
         POTB=POTB + EPS*DSQRT(RHO(I))*C
13    CONTINUE
      PSC=POTA - POTB

      IF (GTEST.OR.STEST) CALL DSCP(N,X,VNEW,BOXLX,BOXLY,BOXLZ,CUTOFF,GTEST,STEST,RHO,VEC)

      RETURN
      END
C
C********************************************************************
C
C Subroutine DSCP calculates first and second derivative 
C analytically for periodic boundary conditions. 
C
C********************************************************************
C
      SUBROUTINE DSCP(N,X,V,BOXLX,BOXLY,BOXLZ,CUTOFF,GTEST,SSTEST,RHO,VEC)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6, J7, NN, MM
      COMMON /POWERS/ NN, MM
      DOUBLE PRECISION RHO(3*N)
      COMMON /INT1/ EPS, C, SIG
      DOUBLE PRECISION VEC(N,N,3),
     1                 RN2(N,N), RN4(N,N), RM4(N,N),
     2                 RM2(N,N), DIST

      DOUBLE PRECISION X(3*N), V(3*N), RHOP(N),
     1                 RHOH(N), TEMP, TEMP1, TEMP2, TEMP3, SIGNN, SIG2MM, SIGMM,
     2                 CUTOFF, EPS, C, SIG, BOXLX, BOXLY, BOXLZ
      LOGICAL GTEST, SSTEST

      SIGMM=SIG**MM
      SIG2MM=SIG**(2*MM)
      SIGNN=SIG**NN
C
C  Store distance matrices.
C
      DO 20 J1=1,N
         RM2(J1,J1)=0.0D0
         RN2(J1,J1)=0.0D0
         RM4(J1,J1)=0.0D0
         RN4(J1,J1)=0.0D0
         RHOP(J1)=RHO(J1)**(1.5)
         RHOH(J1)=1.0D0/DSQRT(RHO(J1))
         DO 10 J2=1,J1-1
            DIST=VEC(J1,J2,1)**2 + VEC(J1,J2,2)**2 +
     1           VEC(J1,J2,3)**2
            IF (DSQRT(DIST).LT.CUTOFF) THEN
               RM2(J2,J1)=DIST**(-(FLOAT(MM)+2.0D0)/2.0D0)
               RN2(J2,J1)=DIST**(-(FLOAT(NN)+2.0D0)/2.0D0)
               RM2(J1,J2)=RM2(J2,J1)
               RN2(J1,J2)=RN2(J2,J1)
               RM4(J2,J1)=DIST**(-(FLOAT(MM)+4.0D0)/2.0D0)
               RN4(J2,J1)=DIST**(-(FLOAT(NN)+4.0D0)/2.0D0)
               RM4(J1,J2)=RM4(J2,J1)
               RN4(J1,J2)=RN4(J2,J1)
            ELSE
               RM2(J2,J1)=0.0D0
               RN2(J2,J1)=0.0D0
               RM4(J2,J1)=0.0D0
               RN4(J2,J1)=0.0D0
               RM2(J1,J2)=RM2(J2,J1)
               RN2(J1,J2)=RN2(J2,J1)
               RM4(J1,J2)=RM4(J2,J1)
               RN4(J1,J2)=RN4(J2,J1)
            ENDIF
10       CONTINUE
20    CONTINUE
C
C Now calculate the gradient analytically.
C
      DO 50 J1=1,N
         DO 40 J2=1,3
            J3=3*(J1-1)+J2
            TEMP=0.0D0
            DO 30 J4=1,N
               TEMP=TEMP+( -NN*EPS*SIGNN*RN2(J1,J4)
     1                    +FLOAT(MM)/2*EPS*C*SIGMM*RM2(J1,J4)*
     2                     (RHOH(J1) + RHOH(J4)))
     3              *VEC(J1,J4,J2)
30          CONTINUE
            V(J3)=TEMP
C           PRINT*,'V(',J3,')=',V(J3)
40       CONTINUE
50    CONTINUE

      IF (.NOT.SSTEST) RETURN
C
C  Now do the hessian. First are the entirely diagonal terms.
C
      DO 80 J1=1,N
         DO 70 J2=1,3
            J3=3*(J1-1)+J2
            HESS(J3,J3)=0.0D0
            TEMP=0.0D0
            DO 60 J4=1,N
               HESS(J3,J3)=HESS(J3,J3)
     1             - NN*EPS*(
     2                     SIGNN*RN2(J4,J1)
     3                   - SIGNN*RN4(J4,J1)*(NN+2)
     4                   * VEC(J1,J4,J2)**2 )
     5             + FLOAT(MM)*EPS*C/2*(
     6                   ( SIGMM*RM2(J4,J1) 
     7                   - SIGMM*RM4(J4,J1)*(MM+2)
     8                   * VEC(J1,J4,J2)**2 )
     9                   *(RHOH(J1) + RHOH(J4))
     A                   + 0.5*SIG2MM*MM*RM2(J4,J1)**2
     B                   * VEC(J1,J4,J2)**2
     C                   /RHOP(J4) )
               TEMP=TEMP+VEC(J1,J4,J2)*RM2(J4,J1)
60          CONTINUE
            TEMP=MM*MM*EPS*C*SIG2MM*TEMP*TEMP/(4.0*RHOP(J1)) 
            HESS(J3,J3)=HESS(J3,J3)+TEMP
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
            DO 100 J4=1,J2-1
               HESS(J3,3*(J1-1)+J4)=0.0D0
               TEMP1=0.0D0
               TEMP2=0.0D0
               DO 90 J5=1,N
                  HESS(J3,3*(J1-1)+J4)=HESS(J3,3*(J1-1)+J4)
     1           + NN*EPS*(SIGNN*RN4(J1,J5)*(NN+2)
     2          * VEC(J1,J5,J2) * VEC(J1,J5,J4) )
     3          +  FLOAT(MM)*EPS*C/2*(
     4            -(RHOH(J1) + RHOH(J5))
     5            * SIGMM*RM4(J5,J1)*(MM+2)
     6          * VEC(J1,J5,J2) * VEC(J1,J5,J4) 
     7           + 0.5*SIG2MM*MM*RM2(J5,J1)**2
     8          * VEC(J1,J5,J2) * VEC(J1,J5,J4)
     9           /RHOP(J5) ) 
                  TEMP1=TEMP1+VEC(J1,J5,J2)*RM2(J5,J1)
                  TEMP2=TEMP2+VEC(J1,J5,J4)*RM2(J5,J1)
90             CONTINUE
               TEMP=MM*MM*EPS*C*SIG2MM*TEMP1*TEMP2/
     1              (4.0*RHOP(J1))  
               HESS(J3,3*(J1-1)+J4)=HESS(J3,3*(J1-1)+J4)+TEMP
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
               HESS(J3,3*(J4-1)+J2)=
     1               NN*EPS*(
     2                     SIGNN*RN2(J4,J1)
     3                   - SIGNN*RN4(J4,J1)*(NN+2)
     4                   * VEC(J1,J4,J2)**2 )
     5             - FLOAT(MM)*EPS*C/2*(
     6                   ( SIGMM*RM2(J4,J1) 
     7                   - SIGMM*RM4(J4,J1)*(MM+2)
     8                   * VEC(J1,J4,J2)**2 )
     9                   *(RHOH(J1) + RHOH(J4)))
               TEMP1=0.0D0
               TEMP2=0.0D0
               TEMP3=0.0D0
               DO 125 J5=1,N
                  TEMP1=TEMP1 + VEC(J4,J5,J2)*RM2(J5,J4)
                  TEMP2=TEMP2 + VEC(J1,J5,J2)*RM2(J5,J1)
                  TEMP3=TEMP3 + VEC(J4,J5,J2)*RM2(J5,J4)
     1                        * VEC(J1,J5,J2)*RM2(J5,J1)/
     2                        RHOP(J5)
125            CONTINUE
              TEMP1=MM*MM*EPS*C*SIG2MM*TEMP1*VEC(J1,J4,J2)  
     1              *RM2(J4,J1)/(4.0*RHOP(J4))  
              TEMP2=MM*MM*EPS*C*SIG2MM*TEMP2*VEC(J1,J4,J2) 
     1              *RM2(J4,J1)/(4.0*RHOP(J1))  
               TEMP3=MM*MM*EPS*C*SIG2MM*TEMP3/4.0
               HESS(J3,3*(J4-1)+J2)=HESS(J3,3*(J4-1)+J2)+TEMP1-TEMP2+TEMP3
               HESS(3*(J4-1)+J2,J3)=HESS(J3,3*(J4-1)+J2)
C              PRINT*,'HESS(',J3,',',3*(J4-1)+J2,')=',HESS(J3,3*(J4-1)+J2)
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
                  HESS(J3,J6)=
     1          -  NN*EPS*SIGNN*RN4(J4,J1)*(NN+2)
     2           * VEC(J1,J4,J2)*VEC(J1,J4,J5)
     3          +  FLOAT(MM)*EPS*C/2*
     4             (RHOH(J1) + RHOH(J4))
     5            * SIGMM*RM4(J4,J1)*(MM+2)
     6           * VEC(J1,J4,J2)*VEC(J1,J4,J5)
                  TEMP1=0.0D0
                  TEMP2=0.0D0
                  TEMP3=0.0D0
                  DO 151 J7=1,N
                     TEMP1=TEMP1+VEC(J4,J7,J5)*RM2(J7,J4)
                     TEMP2=TEMP2+VEC(J1,J7,J2)*RM2(J7,J1)
                     TEMP3=TEMP3+VEC(J1,J7,J2)*RM2(J7,J1)*
     1                           VEC(J4,J7,J5)*RM2(J7,J4)/
     2                            RHOP(J7)
151               CONTINUE
                  TEMP1=MM*MM*EPS*C*SIG2MM*TEMP1*
     1            VEC(J1,J4,J2)*RM2(J4,J1)/(4.0*RHOP(J4)) 
                  TEMP2=MM*MM*EPS*C*SIG2MM*TEMP2*
     1            VEC(J1,J4,J5)*RM2(J4,J1)/(4.0*RHOP(J1)) 
                  TEMP3=MM*MM*EPS*C*SIG2MM*TEMP3/4.0
                  HESS(J3,J6)=HESS(J3,J6)+TEMP1-TEMP2+TEMP3
                  HESS(J6,J3)=HESS(J3,J6)
C              PRINT*,'HESS(',J3,',',J6,')=',HESS(J3,J6)
155            CONTINUE
               DO 156 J5=J2+1,3
                  J6=3*(J4-1)+J5
                  HESS(J3,J6)=
     1          -  NN*EPS*(SIGNN*RN4(J1,J4)*(NN+2)
     2            * VEC(J1,J4,J2)*VEC(J1,J4,J5) )
     3          -  FLOAT(MM)*EPS*C/2*(
     4            -(RHOH(J1) + RHOH(J4))
     5            * SIGMM*RM4(J4,J1)*(MM+2)
     6             * VEC(J1,J4,J2)*VEC(J1,J4,J5) )
                  TEMP1=0.0D0
                  TEMP2=0.0D0
                  TEMP3=0.0D0
                  DO 152 J7=1,N
                     TEMP1=TEMP1+VEC(J4,J7,J5)*RM2(J7,J4)
                     TEMP2=TEMP2+VEC(J1,J7,J2)*RM2(J7,J1)
                     TEMP3=TEMP3+VEC(J1,J7,J2)*RM2(J7,J1)*
     1                           VEC(J4,J7,J5)*RM2(J7,J4)/
     2                            RHOP(J7)
152               CONTINUE
                  TEMP1=MM*MM*EPS*C*SIG2MM*TEMP1*
     1            VEC(J1,J4,J2)*RM2(J4,J1)/(4.0*RHOP(J4))
                  TEMP2=MM*MM*EPS*C*SIG2MM*TEMP2*
     1            VEC(J1,J4,J5)*RM2(J4,J1)/(4.0*RHOP(J1))
                  TEMP3=MM*MM*EPS*C*SIG2MM*TEMP3/4.0
                  HESS(J3,J6)=HESS(J3,J6)+TEMP1-TEMP2+TEMP3
                  HESS(J6,J3)=HESS(J3,J6)
C              PRINT*,'HESS(',J3,',',J6,')=',HESS(J3,J6)
156            CONTINUE
160         CONTINUE
170      CONTINUE
180   CONTINUE
      RETURN
      END
