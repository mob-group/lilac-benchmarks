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
C Subroutine SCDERIV calculates first and second derivative 
C analytically:
C
C********************************************************************
C
      SUBROUTINE SCDIFF(N,X,V,EPS,C,SIG,NN,MM,PSC)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6, J7, NN, MM, I, J
      
      DOUBLE PRECISION X(3*N), 
     1                 V(3*N), RHO(N), RM(N,N),
     2                 RN(N,N),RN2(N,N),
     3                 RM2(N,N), RN4(N,N),EPS,C,
     4                 RM4(N,N), SIG,PSC,DIST,POTA,POTB,TEMP,TEMP1,TEMP2,TEMP3
C
C
C  Store distance matrices.
C
      DO 20 J1=1,N
         RM(J1,J1)=0.0D0
         RN(J1,J1)=0.0D0
         RM2(J1,J1)=0.0D0
         RM4(J1,J1)=0.0D0
         RN2(J1,J1)=0.0D0
         RN4(J1,J1)=0.0D0
         DO 10 J2=1,J1-1
            DIST=(X(3*(J1-1)+1)-X(3*(J2-1)+1))**2
     1          +(X(3*(J1-1)+2)-X(3*(J2-1)+2))**2
     2          +(X(3*(J1-1)+3)-X(3*(J2-1)+3))**2
            RM(J2,J1)=DIST**(-FLOAT(MM)/2.0D0)
            RN(J2,J1)=DIST**(-FLOAT(NN)/2.0D0)
            RM2(J2,J1)=DIST**(-(FLOAT(MM)+2.0D0)/2.0D0)
            RM4(J2,J1)=DIST**(-(FLOAT(MM)+4.0D0)/2.0D0)
            RN2(J2,J1)=DIST**(-(FLOAT(NN)+2.0D0)/2.0D0)
            RN4(J2,J1)=DIST**(-(FLOAT(NN)+4.0D0)/2.0D0)
            RM(J1,J2)=RM(J2,J1)
            RN(J1,J2)=RN(J2,J1)
            RM2(J1,J2)=RM2(J2,J1)
            RM4(J1,J2)=RM4(J2,J1)
            RN2(J1,J2)=RN2(J2,J1)
            RN4(J1,J2)=RN4(J2,J1)
10       CONTINUE
20    CONTINUE
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
C First calculate the potential energy:
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
      PRINT*, ' Sutton-Chen n=', NN,' m=', MM
      PRINT*,'Two-body contribution=',POTA,' eV'
      PRINT*,'Many-body contribution=',-POTB,' eV'
C     PRINT*, 'Total Energy for last step=', PSC,' eV'
C
C Now calculate the gradient analytically.
C
      DO 50 J1=1,N
         DO 40 J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO 30 J4=1,N
               V(J3)=V(J3)+( -NN*EPS*SIG**NN*RN2(J1,J4)
     1                    +FLOAT(MM)/2*EPS*C*SIG**MM*RM2(J1,J4)*
     2                     (1/DSQRT(RHO(J1)) + 1/DSQRT(RHO(J4))))
     3              *(X(J3)-X(3*(J4-1)+J2))
30          CONTINUE
C           PRINT*,'V(',J3,')=',V(J3)
40       CONTINUE
50    CONTINUE
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
     2                     SIG**NN*RN2(J4,J1)
     3                   - SIG**NN*RN4(J4,J1)*(NN+2)
     4                   * (X(J3)-X(3*(J4-1)+J2))**2 )
     5             + FLOAT(MM)*EPS*C/2*(
     6                   ( SIG**MM*RM2(J4,J1) 
     7                   - SIG**MM*RM4(J4,J1)*(MM+2)
     8                   * (X(J3)-X(3*(J4-1)+J2))**2 )
     9                   *(1/DSQRT(RHO(J4)) + 1/DSQRT(RHO(J1)))
     A                   + 0.5*SIG**(2*MM)*MM*RM2(J4,J1)**2
     B                   *(X(J3)-X(3*(J4-1)+J2))**2
     C                   *RHO(J4)**(-1.5) )
               TEMP=TEMP+(X(J3)-X(3*(J4-1)+J2))*RM2(J4,J1)
60          CONTINUE
            TEMP=MM*MM*EPS*C*SIG**(2*MM)*TEMP*TEMP/(4.0*RHO(J1)**(1.5)) 
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
     1           + NN*EPS*(SIG**NN*RN4(J1,J5)*(NN+2)
     2          *(X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4)))
     3          +  FLOAT(MM)*EPS*C/2*(
     4            -(1/DSQRT(RHO(J5)) + 1/DSQRT(RHO(J1)))
     5            * SIG**MM*RM4(J5,J1)*(MM+2)
     6           *(X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4))
     7           + 0.5*SIG**(2*MM)*MM*RM2(J5,J1)**2
     8           *(X(J3)-X(3*(J5-1)+J2))*(X(3*(J1-1)+J4)-X(3*(J5-1)+J4))
     9           *RHO(J5)**(-1.5)  )
                  TEMP1=TEMP1+(X(J3)-X(3*(J5-1)+J2))*RM2(J5,J1)
                  TEMP2=TEMP2+(X(3*(J1-1)+J4)-X(3*(J5-1)+J4))*RM2(J5,J1)
90             CONTINUE
               TEMP=MM*MM*EPS*C*SIG**(2*MM)*TEMP1*TEMP2/
     1              (4.0*RHO(J1)**(1.5))  
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
     2                     SIG**NN*RN2(J4,J1)
     3                   - SIG**NN*RN4(J4,J1)*(NN+2)
     4                   * (X(J3)-X(3*(J4-1)+J2))**2 )
     5             - FLOAT(MM)*EPS*C/2*(
     6                   ( SIG**MM*RM2(J4,J1) 
     7                   - SIG**MM*RM4(J4,J1)*(MM+2)
     8                   * (X(J3)-X(3*(J4-1)+J2))**2 )
     9                   *(1/DSQRT(RHO(J4)) + 1/DSQRT(RHO(J1))))
               TEMP1=0.0D0
               TEMP2=0.0D0
               TEMP3=0.0D0
               DO 125 J5=1,N
                  TEMP1=TEMP1+(X(3*(J4-1)+J2)-X(3*(J5-1)+J2))*RM2(J5,J4) 
                  TEMP2=TEMP2+(X(J3)-X(3*(J5-1)+J2))*RM2(J5,J1) 
                  TEMP3=TEMP3+(X(3*(J4-1)+J2)-X(3*(J5-1)+J2))*RM2(J5,J4) 
     1                       *(X(J3)-X(3*(J5-1)+J2))*RM2(J5,J1)*   
     2                        RHO(J5)**(-1.5)
125            CONTINUE
              TEMP1=MM*MM*EPS*C*SIG**(2*MM)*TEMP1*(X(J3)-X(3*(J4-1)+J2))  
     1              *RM2(J4,J1)/(4.0*RHO(J4)**(1.5))  
              TEMP2=MM*MM*EPS*C*SIG**(2*MM)*TEMP2*(X(J3)-X(3*(J4-1)+J2)) 
     1              *RM2(J4,J1)/(4.0*RHO(J1)**(1.5))  
               TEMP3=MM*MM*EPS*C*SIG**(2*MM)*TEMP3/4.0
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
     1          -  NN*EPS*SIG**NN*RN4(J4,J1)*(NN+2)
     2           *(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6))
     3          +  FLOAT(MM)*EPS*C/2*
     4             (1/DSQRT(RHO(J4)) + 1/DSQRT(RHO(J1)))
     5            * SIG**MM*RM4(J4,J1)*(MM+2)
     6           *(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6))
                  TEMP1=0.0D0
                  TEMP2=0.0D0
                  TEMP3=0.0D0
                  DO 151 J7=1,N
                     TEMP1=TEMP1+(X(J6)-X(3*(J7-1)+J5))*RM2(J7,J4)
                     TEMP2=TEMP2+(X(J3)-X(3*(J7-1)+J2))*RM2(J7,J1)
                     TEMP3=TEMP3+(X(J3)-X(3*(J7-1)+J2))*RM2(J7,J1)*
     1                           (X(J6)-X(3*(J7-1)+J5))*RM2(J7,J4)*   
     2                            RHO(J7)**(-1.5)
151               CONTINUE
                  TEMP1=MM*MM*EPS*C*SIG**(2*MM)*TEMP1*
     1            (X(J3)-X(3*(J4-1)+J2))*RM2(J4,J1)/(4.0*RHO(J4)**(1.5)) 
                  TEMP2=MM*MM*EPS*C*SIG**(2*MM)*TEMP2*
     1            (X(3*(J1-1)+J5)-X(J6))*RM2(J4,J1)/(4.0*RHO(J1)**(1.5)) 
                  TEMP3=MM*MM*EPS*C*SIG**(2*MM)*TEMP3/4.0
                  HESS(J3,J6)=HESS(J3,J6)+TEMP1-TEMP2+TEMP3
                  HESS(J6,J3)=HESS(J3,J6)
C              PRINT*,'HESS(',J3,',',J6,')=',HESS(J3,J6)
155            CONTINUE
               DO 156 J5=J2+1,3
                  J6=3*(J4-1)+J5
                  HESS(J3,J6)=
     1          -  NN*EPS*(SIG**NN*RN4(J1,J4)*(NN+2)
     2           *(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6)))
     3          -  FLOAT(MM)*EPS*C/2*(
     4            -(1/DSQRT(RHO(J4)) + 1/DSQRT(RHO(J1)))
     5            * SIG**MM*RM4(J4,J1)*(MM+2)
     6           *(X(J3)-X(3*(J4-1)+J2))*(X(3*(J1-1)+J5)-X(J6)))
                  TEMP1=0.0D0
                  TEMP2=0.0D0
                  TEMP3=0.0D0
                  DO 152 J7=1,N
                     TEMP1=TEMP1+(X(J6)-X(3*(J7-1)+J5))*RM2(J7,J4)
                     TEMP2=TEMP2+(X(J3)-X(3*(J7-1)+J2))*RM2(J7,J1)
                     TEMP3=TEMP3+(X(J3)-X(3*(J7-1)+J2))*RM2(J7,J1)*
     1                           (X(J6)-X(3*(J7-1)+J5))*RM2(J7,J4)*  
     1                            RHO(J7)**(-1.5)
152               CONTINUE
                  TEMP1=MM*MM*EPS*C*SIG**(2*MM)*TEMP1*
     1            (X(J3)-X(3*(J4-1)+J2))*RM2(J4,J1)/(4.0*RHO(J4)**(1.5)) 
                  TEMP2=MM*MM*EPS*C*SIG**(2*MM)*TEMP2*
     1            (X(3*(J1-1)+J5)-X(J6))*RM2(J4,J1)/(4.0*RHO(J1)**(1.5)) 
                  TEMP3=MM*MM*EPS*C*SIG**(2*MM)*TEMP3/4.0
                  HESS(J3,J6)=HESS(J3,J6)+TEMP1-TEMP2+TEMP3
                  HESS(J6,J3)=HESS(J3,J6)
156            CONTINUE
160         CONTINUE
170      CONTINUE
180   CONTINUE
      RETURN
      END
