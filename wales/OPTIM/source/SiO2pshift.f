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
C  Subroutine Si02PSHIFT calculates the energy, cartesian gradient and second
C  derivative matrix analytically for Si02 using a shifted, truncated potential.
C
C  Adapted for the Si02 glass described by van Beest, Kramer, van Santen, Physical Review Letters, 64, 1990. 
C  
C  Atom types are A(Si) and B(O). The first
C  NTYPEA are A, the next NTYPEB=NATOMS-NTYPEA are B.  
C  (input must be in this order!)
C
C
C  The shifting and truncation for the short-range part of the potential is 
C  as described by Stoddard and Ford, Phys. Rev. A, 8, 1504, 1973. 
C  The Coulomb-part is shifted and truncated via the Wolf-method, 
C  (see Carre et al., Jour. Chem. Phys., 127, 114512, 2007 and 
C  Fennell and Gezelter, Jour. Chem. Phys., 124, 234104, 2006) 
C  
C
C*************************************************************************
C
      SUBROUTINE SIO2PSHIFT(N, X, V, POTEL, BOXLX, BOXLY, BOXLZ, CUTOFF, GTEST, STEST, PTEST, BOXTEST)
      USE KEY
      USE MODHESS 
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, NTYPEA, ANV(N,N,3)
      DOUBLE PRECISION X(3*N), VEC1, VEC2, VEC3, 
     1                 V(3*N), R6, R2DUM,
     2                 Q_A, Q_B, QQ_AA, QQ_AB, QQ_BB, A_AA, A_AB, A_BB,
     3                 B_AA, B_AB, B_BB, C_AA, C_AB, C_BB, EPS_AA, SIG_AA, 
     4                 EPS_AB, EPS_BB, SIG_AB, SIG_BB, POTEL, 
     5                 BOXLX, BOXLY, BOXLZ, IRCUT2AA, IRCUT2AB, IRCUT2BB,   
     6                 CUTOFF, CONSTAA, CONSTBB, CONSTAB,
     7                 RCONSTAA, RCONSTAB, RCONSTBB, CUTAA, CUTAB, CUTBB,
     8                 SIG6_AB, SIG30_AB, SIG6_BB, SIG30_BB
      LOGICAL GTEST, STEST, PTEST, BOXTEST
      COMMON /BIN/ NTYPEA
      COMMON /RCONSTBIN/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB
      parameter(QQ_AA=82.9428D0,QQ_AB=-41.4714D0, QQ_BB=20.7357D0)
C      parameter(A_AA=0.0D0,A_AB=661.62554D0,A_BB=51.03644D0) ! [A] = Eh  
C      parameter(B_AA=0.0D0,B_AB=2.57877477D0,B_BB=1.4605285D0) ! [B] = 1/a
C      parameter(C_AA=0.0D0,C_AB=223.485086D0,C_BB=292.8744D0) ! [C] = Eh a^6 
      parameter(A_AA=0.0D0,A_AB=18003.7572D0,A_BB=1388.7730D0) ! [A] = eV  
      parameter(B_AA=0.0D0,B_AB=4.87318D0,B_BB=2.76000D0) ! [B] = 1/Angstrom
      parameter(C_AA=0.0D0,C_AB=133.5381D0,C_BB=175.0000D0) ! [C] = eV Angstrom^6
      parameter(EPS_AA=0.0D0, EPS_AB=0.003097948D0, EPS_BB=0.00105105048D0)
      parameter(SIG_AA=0.0D0, SIG_AB=1.313635D0, SIG_BB=1.779239D0)

      IF(MOD(N,3).EQ.0) THEN
        NTYPEA = N/3
      ELSE
        WRITE(*,*) "SiO2pshift.f> warning: number of atoms must correspond to stoichiometry of SiO2!"
        STOP
      ENDIF
   
C      IF (PTEST) THEN
C         WRITE(102,*) "SiO2pshift.f> warning: PTEST not implemented"
C         STOP
C      ENDIF
      IF (BOXTEST) THEN
         WRITE(102,*) "SiO2pshift.f> warning: box derivatives not implemented"
         STOP
      ENDIF

      SIG30_AB = SIG_AB**30
      SIG6_AB = SIG_AB**6
      SIG30_BB = SIG_BB**30
      SIG6_BB = SIG_BB**6      

      CUTAA=CUTOFF 
      CUTAB=CUTOFF ! change? and set CUTOFF!
      CUTBB=CUTOFF !
      IRCUT2AA = 1.D0/CUTAA**2
      IRCUT2AB = 1.D0/CUTAB**2
      IRCUT2BB = 1.D0/CUTBB**2
      
      CONSTAA = 0.0D0
      CONSTAB = -A_AB*EXP(-B_AB*CUTAB)*(1.0D0+B_AB*CUTAB/2.0D0)
     &          +4.0D0*C_AB*IRCUT2AB**3
     &          +4.0D0*EPS_AB*(-16.D0*SIG30_AB/CUTAB**30+4.0D0*SIG6_AB/CUTAB**6)
      CONSTBB = -A_BB*EXP(-B_BB*CUTBB)*(1.0D0+B_BB*CUTBB/2.0D0)
     &          +4.0D0*C_BB*IRCUT2BB**3
     &          +4.0D0*EPS_BB*(-16.D0*SIG30_BB/CUTBB**30+4.0D0*SIG6_BB/CUTBB**6)
      RCONSTAA = 0.0D0
      RCONSTAB = A_AB*B_AB/(2.0D0*CUTAB)*EXP(-B_AB*CUTAB)
     &          -3.0D0*C_AB/CUTAB**8
     &          +4.0D0*EPS_AB*(15.0D0*SIG30_AB/CUTAB**32-3.0D0*SIG6_AB/CUTAB**8)
      RCONSTBB = A_BB*B_BB/(2.0D0*CUTBB)*EXP(-B_BB*CUTBB)
     &          -3.0D0*C_BB/CUTBB**8
     &          +4.0D0*EPS_BB*(15.0D0*SIG30_BB/CUTBB**32-3.0D0*SIG6_BB/CUTBB**8)


C  Work out cutoff for potential. Two particles interact if r<c, but
C  we will use the equivalent condition 1/r^2 > 1/c^2.
C
C  Deal with any atoms that have left the box.
C
      IF ((.NOT.FIXIMAGE).AND.(.NOT.NORESET)) THEN
         DO J1=1,N
            J2 = 3*(J1-1)
C           IF (ANINT(X(J2+1)/BOXLX).NE.0) PRINT*,'resetting X for J1,X,Xnew=',J1,X(J2+1),X(J2+1) - BOXLX*ANINT(X(J2+1)/BOXLX)
C           IF (ANINT(X(J2+2)/BOXLX).NE.0) PRINT*,'resetting Y for J1,X,Xnew=',J1,X(J2+2),X(J2+2) - BOXLY*ANINT(X(J2+2)/BOXLY)
C           IF (ANINT(X(J2+3)/BOXLX).NE.0) PRINT*,'resetting Z for J1,X,Xnew=',J1,X(J2+3),X(J2+3) - BOXLZ*ANINT(X(J2+3)/BOXLZ)
            X(J2+1)=X(J2+1) - BOXLX*ANINT(X(J2+1)/BOXLX)
            X(J2+2)=X(J2+2) - BOXLY*ANINT(X(J2+2)/BOXLY)
            X(J2+3)=X(J2+3) - BOXLZ*ANINT(X(J2+3)/BOXLZ)
         ENDDO
      ENDIF
C
C  Calculate interatomic vectors using the minimum image convention.
C  VEC(i,j,alpha) is the alpha (x, y or z) component of the vector between
C  atoms i and j.
C

      POTEL=0.0D0

      ! sn402: We will profile this at some point and decide whether it needs to be rewritten
      ! in a more efficient manner.
      IF (STEST) THEN
         ! calculate first and second derivatives of potential (hessian)
         V(1:3*N)=0.0D0 ! not necessary if .not.GTEST .and. .not.STEST
         HESS(:,:)=0.0D0 

         DO J1=1,NTYPEA
            DO J2=J1+1,NTYPEA
               CALL SIO2PSHIFT_UPDATE_PAIRS(N, X, QQ_AA, A_AA, B_AA, C_AA,
     &           EPS_AA, SIG_AA, RCONSTAA, CONSTAA, CUTAA, IRCUT2AA,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, V)
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            DO J2=NTYPEA+1,N
               CALL SIO2PSHIFT_UPDATE_PAIRS(N, X, QQ_AB, A_AB, B_AB, C_AB,
     &           EPS_AB, SIG_AB, RCONSTAB, CONSTAB, CUTAB, IRCUT2AB,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, V)
            ENDDO
         ENDDO
         DO J1=NTYPEA+1,N
            DO J2=J1+1,N
              CALL SIO2PSHIFT_UPDATE_PAIRS(N, X, QQ_BB, A_BB, B_BB, C_BB,
     &           EPS_BB, SIG_BB, RCONSTBB, CONSTBB, CUTBB, IRCUT2BB,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, V)
            ENDDO
         ENDDO

       ELSEIF (GTEST) THEN 
         ! calculate first derivatives
         V(1:3*N)=0.0D0 ! not necessary if .not.GTEST .and. .not.STEST

         DO J1=1,NTYPEA
            DO J2=J1+1,NTYPEA
               CALL SIO2PSHIFT_UPDATE_PAIRG(N, X, QQ_AA, A_AA, B_AA, C_AA,
     &           EPS_AA, SIG_AA, RCONSTAA, CONSTAA, CUTAA, IRCUT2AA,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, V)
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            DO J2=NTYPEA+1,N
               CALL SIO2PSHIFT_UPDATE_PAIRG(N, X, QQ_AB, A_AB, B_AB, C_AB,
     &           EPS_AB, SIG_AB, RCONSTAB, CONSTAB, CUTAB, IRCUT2AB,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, V)
            ENDDO
         ENDDO
         DO J1=NTYPEA+1,N
            DO J2=J1+1,N
               CALL SIO2PSHIFT_UPDATE_PAIRG(N, X, QQ_BB, A_BB, B_BB, C_BB,
     &           EPS_BB, SIG_BB, RCONSTBB, CONSTBB, CUTBB, IRCUT2BB,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, V)
            ENDDO
         ENDDO

       ELSE 
         ! calculate potential only

         DO J1=1,NTYPEA
            DO J2=J1+1,NTYPEA
               ! A-A interaction
               CALL SIO2PSHIFT_UPDATE_PAIR(N, X, QQ_AA, A_AA, B_AA, C_AA,
     &           EPS_AA, SIG_AA, RCONSTAA, CONSTAA, CUTAA, IRCUT2AA,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE)
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            DO J2=NTYPEA+1,N
               ! A-B interaction
               CALL SIO2PSHIFT_UPDATE_PAIR(N, X, QQ_AB, A_AB, B_AB, C_AB,
     &           EPS_AB, SIG_AB, RCONSTAB, CONSTAB, CUTAB, IRCUT2AB,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE)
            ENDDO
         ENDDO
         DO J1=NTYPEA+1,N
            DO J2=J1+1,N
               ! B-B interaction
               CALL SIO2PSHIFT_UPDATE_PAIR(N, X, QQ_BB, A_BB, B_BB, C_BB,
     &           EPS_BB, SIG_BB, RCONSTBB, CONSTBB, CUTBB, IRCUT2BB,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE)
            ENDDO
         ENDDO

      ENDIF

C      PRINT *,'mb2098: POTEL=',POTEL
      RETURN
      END



C*****************************************************************************

      SUBROUTINE SIO2PSHIFT_UPDATE_PAIRS(N, X, QQ, A, B, C,
     &           EPS, SIG, RCONST, CONST, RCUT, IRCUT2,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, V)
        !Calculate the potential energy and gradient and hessian between atoms
        !j1 and j2.
        !Update also ANV
        USE MODHESS

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: J1, J2, N
        DOUBLE PRECISION, INTENT(IN) :: QQ, A, B, C, RCONST, CONST, RCUT, IRCUT2
        DOUBLE PRECISION, INTENT(IN) :: EPS, SIG
        DOUBLE PRECISION, INTENT(IN) :: X(3*N)
        DOUBLE PRECISION, INTENT(IN) :: BOXLX, BOXLY, BOXLZ
        LOGICAL, INTENT(IN) :: FIXIMAGE
        DOUBLE PRECISION, INTENT(OUT) :: POTEL
        INTEGER, INTENT(INOUT) :: ANV(N,N,3)
        DOUBLE PRECISION, INTENT(INOUT) :: V(3*N)

        DOUBLE PRECISION :: R, R2, R6, VEC(3), G, R8, R30, R32, FR2, TEMP, COULPOT
        DOUBLE PRECISION :: SIG6, SIG30
        INTEGER :: J3, J4, J5, J6

        SIG6 = SIG**6
        SIG30 = SIG**30

        J3=3*(J1-1)
        J4=3*(J2-1)
        !find particle separation with suitable boundary conditions
        VEC(1)=X(J3+1)-X(J4+1)
        VEC(2)=X(J3+2)-X(J4+2)
        VEC(3)=X(J3+3)-X(J4+3)
        IF (.NOT.FIXIMAGE) THEN
          ANV(J2,J1,1)=NINT(VEC(1)/BOXLX)
          ANV(J2,J1,2)=NINT(VEC(2)/BOXLY)
          ANV(J2,J1,3)=NINT(VEC(3)/BOXLZ)
          ANV(J1,J2,1)=-ANV(J2,J1,1)
          ANV(J1,J2,2)=-ANV(J2,J1,2)
          ANV(J1,J2,3)=-ANV(J2,J1,3)
        ENDIF
        VEC(1)=VEC(1)-BOXLX*ANV(J2,J1,1)
        VEC(2)=VEC(2)-BOXLY*ANV(J2,J1,2)
        VEC(3)=VEC(3)-BOXLZ*ANV(J2,J1,3)


        R2=VEC(1)**2+VEC(2)**2+VEC(3)**2
        R = SQRT(R2)
        R2=1.0D0/R2
        IF (R2.GT.IRCUT2) THEN
          !update potential energy
          R6=R2**3
          R30=R2**15
          COULPOT = QQ*(1.0D0/R - 1.0D0/RCUT + IRCUT2*(R-RCUT))
          POTEL = POTEL + COULPOT+A*EXP(-B*R)-C*R6+RCONST/R2+CONST
     &            +4.0D0*EPS*(SIG30*R30-SIG6*R6) 
          !update derivative of potential energy
          R8=R6*R2
          R32=R8**4
          G = QQ*(-1.0D0/R**3 + 1.0D0/R*IRCUT2) - A*B/R*EXP(-B*R) + 6.0D0*C*R8 + 2.0D0*RCONST
     &        +4.0D0*EPS*(6.0D0*SIG6*R8-30.0D0*SIG30*R32)
          DO J5=1,3
            !careful with signs: both VEC and V are antisymmetric in J1<->J2
            V(J3+J5)=V(J3+J5)+G*VEC(J5)
            V(J4+J5)=V(J4+J5)-G*VEC(J5)
          END DO

          !calculate the second derivative, and update the components of the
          !hessian
          FR2 = (QQ*(3.0D0/R**3 - 1.0D0/R*IRCUT2) +A*B*EXP(-B*R)*(1.0D0/R+B)-48.0D0*C*R8
     &          +4.0D0*EPS*(960.0D0*SIG30*R32-48.0D0*SIG6*R8))*R2
C  js850> same cartesian coordinate
          DO J5=1,3
            TEMP = FR2 * VEC(J5)**2 + G
C  Now do the hessian. First are the entirely diagonal terms.
            HESS(J3+J5,J3+J5) = HESS(J3+J5,J3+J5) + TEMP
            HESS(J4+J5,J4+J5) = HESS(J4+J5,J4+J5) + TEMP
C  Case III, different atoms, same cartesian coordinate.
            HESS(J4+J5,J3+J5) = - TEMP
            HESS(J3+J5,J4+J5) = - TEMP
          END DO
C  js850> different cartesian coordinates
          DO J5=1,3
            DO J6=J5+1,3
              TEMP = FR2 * VEC(J5) * VEC(J6) 
C  Next are the terms where x_i and x_j are on the same atom
C  but are different, e.g. y and z.
              HESS(J3+J5,J3+J6) = HESS(J3+J5,J3+J6) + TEMP
              HESS(J4+J5,J4+J6) = HESS(J4+J5,J4+J6) + TEMP
              ! J5 <-> J6
              HESS(J3+J6,J3+J5) = HESS(J3+J6,J3+J5) + TEMP
              HESS(J4+J6,J4+J5) = HESS(J4+J6,J4+J5) + TEMP
C  Case IV: different atoms and different cartesian coordinates.
              HESS(J4+J5,J3+J6) = -TEMP
              HESS(J3+J5,J4+J6) = -TEMP
              ! J5 <-> J6
              HESS(J4+J6,J3+J5) = -TEMP
              HESS(J3+J6,J4+J5) = -TEMP
            END DO
          END DO

        ELSE
          !set hessian to zero

          DO J5=1,3
            DO J6=J5,3
              HESS(J4+J5,J3+J6) = 0.D0
              HESS(J3+J5,J4+J6) = 0.D0
              ! J5 <-> J6
              HESS(J4+J6,J3+J5) = 0.D0
              HESS(J3+J6,J4+J5) = 0.D0
            END DO
          END DO

        ENDIF

      END SUBROUTINE SIO2PSHIFT_UPDATE_PAIRS



C*****************************************************************************

      SUBROUTINE SIO2PSHIFT_UPDATE_PAIRG(N, X, QQ, A, B, C, 
     &           EPS, SIG, RCONST, CONST, RCUT, IRCUT2,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, V)
        !Calculate the potential energy and gradient between atoms j1 and j2.
        !Update also ANV

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: J1, J2, N
        DOUBLE PRECISION, INTENT(IN) :: QQ, A, B, C, RCONST, CONST, RCUT, IRCUT2
        DOUBLE PRECISION, INTENT(IN) :: EPS, SIG
        DOUBLE PRECISION, INTENT(IN) :: X(3*N)
        DOUBLE PRECISION, INTENT(IN) :: BOXLX, BOXLY, BOXLZ
        LOGICAL, INTENT(IN) :: FIXIMAGE
        DOUBLE PRECISION, INTENT(OUT) :: POTEL
        INTEGER, INTENT(INOUT) :: ANV(N,N,3)
        DOUBLE PRECISION, INTENT(INOUT) :: V(3*N)

        DOUBLE PRECISION :: R, R2, R6, VEC(3), G, R8, R30, R32, COULPOT
        DOUBLE PRECISION :: SIG6, SIG30
        INTEGER :: J3, J4, J5

        SIG6 = SIG**6
        SIG30 = SIG**30

        J3=3*(J1-1)
        J4=3*(J2-1)
        !find particle separation with suitable boundary conditions
        VEC(1)=X(J3+1)-X(J4+1)
        VEC(2)=X(J3+2)-X(J4+2)
        VEC(3)=X(J3+3)-X(J4+3)
        IF (.NOT.FIXIMAGE) THEN
          ANV(J2,J1,1)=NINT(VEC(1)/BOXLX)
          ANV(J2,J1,2)=NINT(VEC(2)/BOXLY)
          ANV(J2,J1,3)=NINT(VEC(3)/BOXLZ)
          ANV(J1,J2,1)=-ANV(J2,J1,1) 
          ANV(J1,J2,2)=-ANV(J2,J1,2)
          ANV(J1,J2,3)=-ANV(J2,J1,3)
        ENDIF
        VEC(1)=VEC(1)-BOXLX*ANV(J2,J1,1)
        VEC(2)=VEC(2)-BOXLY*ANV(J2,J1,2)
        VEC(3)=VEC(3)-BOXLZ*ANV(J2,J1,3)

        !calculate the potential
        R2=VEC(1)**2+VEC(2)**2+VEC(3)**2
        R=SQRT(R2)
        R2=1.0D0/R2
        IF (R2.GT.IRCUT2) THEN
          R30=R2**15
          R6=R2**3
          COULPOT = QQ*(1.0D0/R - 1.0D0/RCUT + IRCUT2*(R-RCUT))
          POTEL = POTEL + COULPOT+A*EXP(-B*R)-C*R6+RCONST/R2+CONST
     &            +4.0D0*EPS*(SIG30*R30-SIG6*R6)
          !update derivative of potential energy
          R8=R6*R2
          R32=R8**4
          G = QQ*(-1.0D0/R**3 + 1.0D0/R*IRCUT2) - A*B/R*EXP(-B*R) + 6.0D0*C*R8 + 2.0D0*RCONST
     &        +4.0D0*EPS*(6.0D0*SIG6*R8-30.0D0*SIG30*R32)
          DO J5=1,3
            V(J3+J5)=V(J3+J5)+G*VEC(J5)
            V(J4+J5)=V(J4+J5)-G*VEC(J5)
          END DO  
        ENDIF

      END SUBROUTINE SIO2PSHIFT_UPDATE_PAIRG


C*****************************************************************************

      SUBROUTINE SIO2PSHIFT_UPDATE_PAIR(N, X, QQ, A, B, C, EPS, SIG, RCONST, 
     &           CONST, RCUT, IRCUT2, POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE)
        !calculate the potential energy between atoms j1 and j2 
        !update ANV

        IMPLICIT NONE

        INTEGER, INTENT(IN) :: J1, J2, N
        DOUBLE PRECISION, INTENT(IN) :: QQ, A, B, C, RCONST, CONST, RCUT, IRCUT2
        DOUBLE PRECISION, INTENT(IN) :: EPS, SIG
        DOUBLE PRECISION, INTENT(IN) :: X(3*N)
        DOUBLE PRECISION, INTENT(IN) :: BOXLX, BOXLY, BOXLZ
        LOGICAL, INTENT(IN) :: FIXIMAGE
        DOUBLE PRECISION, INTENT(OUT) :: POTEL
        INTEGER, INTENT(INOUT) :: ANV(N,N,3)

        DOUBLE PRECISION :: R, R2, R6, R30, VEC1, VEC2, VEC3, COULPOT
        DOUBLE PRECISION :: SIG6, SIG30
        INTEGER :: J3, J4

        SIG6 = SIG**6
        SIG30 = SIG**30

        J3=3*(J1-1)
        J4=3*(J2-1)
        !find particle separation with suitable boundary conditions
        VEC1=X(J3+1)-X(J4+1)
        VEC2=X(J3+2)-X(J4+2)
        VEC3=X(J3+3)-X(J4+3)
        IF (.NOT.FIXIMAGE) THEN
          ANV(J2,J1,1)=NINT(VEC1/BOXLX)
          ANV(J2,J1,2)=NINT(VEC2/BOXLY)
          ANV(J2,J1,3)=NINT(VEC3/BOXLZ)
          ANV(J1,J2,1)=-ANV(J2,J1,1) 
          ANV(J1,J2,2)=-ANV(J2,J1,2)
          ANV(J1,J2,3)=-ANV(J2,J1,3)
        ENDIF
        VEC1=VEC1-BOXLX*ANV(J2,J1,1)
        VEC2=VEC2-BOXLY*ANV(J2,J1,2)
        VEC3=VEC3-BOXLZ*ANV(J2,J1,3)
        !calculate the potential
        R2=VEC1**2+VEC2**2+VEC3**2

        R=SQRT(R2)
        R2=1.0D0/R2
        R30=R2**15

        IF (R2.GT.IRCUT2) THEN
          R6=R2**3
          COULPOT = QQ*(1.0D0/R - 1.0D0/RCUT + IRCUT2*(R-RCUT))
          POTEL = POTEL + COULPOT+A*EXP(-B*R)-C*R6+RCONST/R2+CONST
     &            +4.0D0*EPS*(SIG30*R30-SIG6*R6)
        ENDIF

      END SUBROUTINE SIO2PSHIFT_UPDATE_PAIR
