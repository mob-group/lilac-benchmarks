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
C  Subroutine LJPSHIFTBIN calculates the energy, cartesian gradient and second
C  derivative matrix analytically for Lennard-Jones in reduced units
C  (epsilon=sigma=1) using a shifted, truncated potential.
C
C  Adapted for the binary LJ glass described by Sastry, Debenetti and
C  Stillinger, Nature, 393, 554, 1998. Atom types are A and B. The first
C  NTYPEA are A, the next NBTYPE=NATOMS-NTYPEA are B. epsilon and sigma for A are the
C  units of energy and distance, so we also need EPSAB, EPSAB, SIGAB and
C  SIGAA in these units. Sastry et al. density is 1.2 i.e. a box length
C  of 5.975206 for 256 atoms. 
C
C  The shifting and truncation is as described by Stoddard and Ford, Phys.
C  Rev. A, 8, 1504, 1973. Atom type 'LS'.
C
C
C*************************************************************************
C
      SUBROUTINE LJPSHIFTBIN(N, X, V, POTEL, BOXLX, BOXLY, BOXLZ, CUTOFF, GTEST, STEST, PTEST, BOXTEST)
      USE KEY
      USE MODHESS 
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, NTYPEA, ANV(N,N,3)
      DOUBLE PRECISION X(3*N), VEC1, VEC2, VEC3, VBOX(3), SBOX(3),
     1                 V(3*N), R6, R2DUM,
     2                 EPSAB, EPSBB, SIGAB, SIGBB,
     4                 POTEL, BOXLX, BOXLY, BOXLZ, SIGAB6, SIGAB12, SIGBB6, SIGBB12,
     5                 IRCUT2AA, IRCUT2AB, IRCUT2BB, CUTOFF, SIGRCAA6, SIGRCAA12,
     6                 SIGRCAB6, SIGRCAB12, SIGRCBB6, SIGRCBB12, CONSTAA, CONSTBB, CONSTAB,
     7                 RCONSTAA, RCONSTAB, RCONSTBB, CUTAA, CUTAB, CUTBB
      LOGICAL GTEST, STEST, PTEST, BOXTEST
      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      COMMON /RCONSTBIN/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB
      COMMON /BOXDERIVS/ VBOX, SBOX
C     COMMON /ANVCOMMON/ ANV

      !WRITE(*,*) "GTEST, STEST, PTEST, BOXTEST ", GTEST, STEST, PTEST, BOXTEST
      CUTAA=CUTOFF
      CUTAB=CUTOFF*SIGAB
      CUTBB=CUTOFF*SIGBB
      IRCUT2AA = 1.D0/CUTAA**2
      IRCUT2AB = 1.D0/CUTAB**2
      IRCUT2BB = 1.D0/CUTBB**2
      SIGAB6=SIGAB**6
      SIGAB12=SIGAB6**2
      SIGBB6=SIGBB**6
      SIGBB12=SIGBB6**2
      SIGRCAA6= 1.0D0/CUTAA**6
      SIGRCAA12=SIGRCAA6**2
      SIGRCAB6=SIGAB6/CUTAB**6
      SIGRCAB12=SIGRCAB6**2
      SIGRCBB6=SIGBB6/CUTBB**6
      SIGRCBB12=SIGRCBB6**2
      CONSTAA=4.0D0*SIGRCAA6-7.0D0*SIGRCAA12
      CONSTAB=4.0D0*SIGRCAB6-7.0D0*SIGRCAB12
      CONSTBB=4.0D0*SIGRCBB6-7.0D0*SIGRCBB12
      RCONSTAA=(6.0D0*SIGRCAA12-3.0D0*SIGRCAA6)/CUTAA**2
      RCONSTAB=(6.0D0*SIGRCAB12-3.0D0*SIGRCAB6)/CUTAB**2
      RCONSTBB=(6.0D0*SIGRCBB12-3.0D0*SIGRCBB6)/CUTBB**2
C
C  Work out cutoff for potential. Two particles interact if r<c, but
C  we will use the equivalent condition 1/r^2 > 1/c^2.
C
C     IF (PTEST) PRINT*,'Cutoff used = ',CUTOFF
C
C  Deal with any atoms that have left the box.
C
C     PRINT*,'In ljpshift FIXIMAGE,NORESET,GTEST=',FIXIMAGE,NORESET,GTEST
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

      IF (PTEST) WRITE(102,*) "ljpshift.f> warning: PTEST not implimented"
      !js850> the box derivatives and the potential derivatives (also hessian)
      !are not calculated together.  They could easily be (and it would be much
      !faster to do them together), but I'm maintaining the structure that
      !already existed.
      IF (BOXTEST .AND. (GTEST .OR. STEST) ) 
     &   WRITE(102,*) "ljpshift.f> warning: box derivatives and potential derivatives are implimented seperately"

      POTEL=0.0D0
      IF (BOXTEST) THEN
         ! calculte box derivatives only
         VBOX(1:3) = 0.D0
         SBOX(1:3) = 0.D0

         DO J1=1,NTYPEA
            DO J2=J1+1,NTYPEA
               CALL LJPSHIFT_UPDATE_PAIRBOX(N, X,  1.D0,   1.D0, RCONSTAA,
     &           CONSTAA, IRCUT2AA,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, VBOX,
     &           SBOX, FRACTIONAL )
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            DO J2=NTYPEA+1,N
               CALL LJPSHIFT_UPDATE_PAIRBOX(N, X, EPSAB, SIGAB6, RCONSTAB,
     &           CONSTAB, IRCUT2AB,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, VBOX,
     &           SBOX, FRACTIONAL )
            ENDDO
         ENDDO
         DO J1=NTYPEA+1,N
            DO J2=J1+1,N
               CALL LJPSHIFT_UPDATE_PAIRBOX(N, X, EPSBB, SIGBB6, RCONSTBB,
     &           CONSTBB, IRCUT2BB,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, VBOX,
     &           SBOX, FRACTIONAL )
            ENDDO
         ENDDO

      ELSEIF (STEST) THEN
         ! calculate first and second derivatives of potential (hessian)
         V(1:3*N)=0.0D0 ! not necessary if .not.GTEST .and. .not.STEST
         HESS(:,:)=0.0D0 

         DO J1=1,NTYPEA
            DO J2=J1+1,NTYPEA
               CALL LJPSHIFT_UPDATE_PAIRS(N, X,  1.D0,   1.D0, RCONSTAA,
     &           CONSTAA, IRCUT2AA,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, V)
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            DO J2=NTYPEA+1,N
               CALL LJPSHIFT_UPDATE_PAIRS(N, X, EPSAB, SIGAB6, RCONSTAB,
     &           CONSTAB, IRCUT2AB,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, V)
            ENDDO
         ENDDO
         DO J1=NTYPEA+1,N
            DO J2=J1+1,N
               CALL LJPSHIFT_UPDATE_PAIRS(N, X, EPSBB, SIGBB6, RCONSTBB,
     &           CONSTBB, IRCUT2BB,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, V)
            ENDDO
         ENDDO

       ELSEIF (GTEST) THEN 
         ! calculate first derivatives
         V(1:3*N)=0.0D0 ! not necessary if .not.GTEST .and. .not.STEST

         DO J1=1,NTYPEA
            DO J2=J1+1,NTYPEA
               CALL LJPSHIFT_UPDATE_PAIRG(N, X,  1.D0,   1.D0, RCONSTAA,
     &           CONSTAA, IRCUT2AA,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, V)
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            DO J2=NTYPEA+1,N
               CALL LJPSHIFT_UPDATE_PAIRG(N, X, EPSAB, SIGAB6, RCONSTAB,
     &           CONSTAB, IRCUT2AB,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, V)
            ENDDO
         ENDDO
         DO J1=NTYPEA+1,N
            DO J2=J1+1,N
               CALL LJPSHIFT_UPDATE_PAIRG(N, X, EPSBB, SIGBB6, RCONSTBB,
     &           CONSTBB, IRCUT2BB,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, V)
            ENDDO
         ENDDO

       ELSE 
         ! calculate potential only

         DO J1=1,NTYPEA
            DO J2=J1+1,NTYPEA
               CALL LJPSHIFT_UPDATE_PAIR(N, X,  1.D0,   1.D0, RCONSTAA,
     &           CONSTAA, IRCUT2AA,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE)
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            DO J2=NTYPEA+1,N
               CALL LJPSHIFT_UPDATE_PAIR(N, X, EPSAB, SIGAB6, RCONSTAB,
     &           CONSTAB, IRCUT2AB,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE)
            ENDDO
         ENDDO
         DO J1=NTYPEA+1,N
            DO J2=J1+1,N
               CALL LJPSHIFT_UPDATE_PAIR(N, X, EPSBB, SIGBB6, RCONSTBB,
     &           CONSTBB, IRCUT2BB,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE)
            ENDDO
         ENDDO

      ENDIF
      POTEL = POTEL * 4.D0

      IF ((.NOT.GTEST).AND.(.NOT.STEST).AND.(.NOT.BOXTEST)) RETURN
      IF (BOXTEST) THEN
         IF (FRACTIONAL) THEN
            VBOX(1)=VBOX(1)/BOXLX+PRESS*BOXLY*BOXLZ
            VBOX(2)=VBOX(2)/BOXLY+PRESS*BOXLX*BOXLZ
            VBOX(3)=VBOX(3)/BOXLZ+PRESS*BOXLY*BOXLX
         ELSE
            VBOX(1)=VBOX(1)+PRESS*BOXLY*BOXLZ
            VBOX(2)=VBOX(2)+PRESS*BOXLX*BOXLZ
            VBOX(3)=VBOX(3)+PRESS*BOXLY*BOXLX
         ENDIF
      ENDIF

      IF ((.NOT.STEST).AND.(.NOT.BOXTEST)) RETURN
      IF (BOXTEST) THEN
         IF (FRACTIONAL) THEN
            SBOX(1)=SBOX(1)/BOXLX**2
            SBOX(2)=SBOX(2)/BOXLY**2
            SBOX(3)=SBOX(3)/BOXLZ**2
         ENDIF
      ENDIF

      RETURN
      END

C*****************************************************************************
  

C*****************************************************************************

C
C  Same as above without periodic boundaries.
C
      SUBROUTINE LJPSHIFTBINC(N, X, V, POTEL, CUTOFF, GTEST, STEST, PTEST)
      USE KEY
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, NTYPEA
      DOUBLE PRECISION X(3*N), VEC1, VEC2, VEC3, 
     1                 V(3*N), R2(N,N),  R6, R2DUM,
     2                 R8(N,N), G(N,N), EPSAB, EPSBB, SIGAB, SIGBB,
     3                 R14(N,N), F(N,N), 
     4                 POTEL, SIGAB6, SIGAB12, SIGBB6, SIGBB12,
     5                 VEC(N,N,3), IRCUT2AA, IRCUT2AB, IRCUT2BB, CUTOFF, SIGRCAA6, SIGRCAA12,
     6                 SIGRCAB6, SIGRCAB12, SIGRCBB6, SIGRCBB12, CONSTAA, CONSTBB, CONSTAB,
     7                 RCONSTAA, RCONSTAB, RCONSTBB, CUTAA, CUTAB, CUTBB
      LOGICAL GTEST, STEST, PTEST
      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      COMMON /RCONSTBIN/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB

      CUTAA=CUTOFF
      CUTAB=CUTOFF*SIGAB
      CUTBB=CUTOFF*SIGBB
      IRCUT2AA = 1.D0/CUTAA**2
      IRCUT2AB = 1.D0/CUTAB**2
      IRCUT2BB = 1.D0/CUTBB**2
      SIGAB6=SIGAB**6
      SIGAB12=SIGAB6**2
      SIGBB6=SIGBB**6
      SIGBB12=SIGBB6**2
      SIGRCAA6= 1.0D0/CUTAA**6
      SIGRCAA12=SIGRCAA6**2
      SIGRCAB6=SIGAB6/CUTAB**6
      SIGRCAB12=SIGRCAB6**2
      SIGRCBB6=SIGBB6/CUTBB**6
      SIGRCBB12=SIGRCBB6**2
      CONSTAA=4.0D0*SIGRCAA6-7.0D0*SIGRCAA12
      CONSTAB=4.0D0*SIGRCAB6-7.0D0*SIGRCAB12
      CONSTBB=4.0D0*SIGRCBB6-7.0D0*SIGRCBB12
      RCONSTAA=(6.0D0*SIGRCAA12-3.0D0*SIGRCAA6)/CUTAA**2
      RCONSTAB=(6.0D0*SIGRCAB12-3.0D0*SIGRCAB6)/CUTAB**2
      RCONSTBB=(6.0D0*SIGRCBB12-3.0D0*SIGRCBB6)/CUTBB**2
C
C  Work out cutoff for potential. Two particles interact if r<c, but
C  we will use the equivalent condition 1/r^2 > 1/c^2.
C
C     IF (PTEST) PRINT*,'Cutoff used = ',CUTOFF
C
C  Calculate interatomic vectors 
C  VEC(i,j,alpha) is the alpha (x, y or z) component of the vector between
C  atoms i and j.
C
      POTEL=0.0D0
      IF (GTEST) THEN
         DO J1=1, N
            VEC(J1,J1,1)=0.0D0
            VEC(J1,J1,2)=0.0D0
            VEC(J1,J1,3)=0.0D0
            J3=3*(J1-1)
            DO J2=J1+1, N
               J4=3*(J2-1)
               VEC(J2,J1,1)=X(J3+1)-X(J4+1)
               VEC(J2,J1,2)=X(J3+2)-X(J4+2)
               VEC(J2,J1,3)=X(J3+3)-X(J4+3)
               VEC(J1,J2,1)=-VEC(J2,J1,1)
               VEC(J1,J2,2)=-VEC(J2,J1,2)
               VEC(J1,J2,3)=-VEC(J2,J1,3)
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            R2(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,NTYPEA
               R2DUM=1.0D0/(VEC(J1,J2,1)**2+VEC(J1,J2,2)**2+VEC(J1,J2,3)**2)
               R2(J2,J1)=R2DUM
               R6=R2(J2,J1)**3
               IF (R2DUM.GT.IRCUT2AA) THEN
                  POTEL=POTEL + R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA  ! AA
C                 WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA,POTEL
C                 IF (R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA.GT.10.0D0) 
C    1            WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA,POTEL
               ENDIF
               R8(J2,J1)=R6*R2(J2,J1)
               R14(J2,J1)=R8(J2,J1)*R6
               R2(J1,J2)=R2(J2,J1)
               R8(J1,J2)=R8(J2,J1)
               R14(J1,J2)=R14(J2,J1)
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            DO J2=NTYPEA+1,N
               R2DUM=1.0D0/(VEC(J1,J2,1)**2+VEC(J1,J2,2)**2+VEC(J1,J2,3)**2)
               R2(J2,J1)=R2DUM
               R6=R2(J2,J1)**3
               IF (R2DUM.GT.IRCUT2AB) THEN
                  POTEL=POTEL + EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB) ! AB
C                 WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                         EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB),POTEL
C                 IF (EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB).GT.10.0D0) 
C    1            WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                         EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB),POTEL
               ENDIF
               R8(J2,J1)=R6*R2(J2,J1)
               R14(J2,J1)=R8(J2,J1)*R6
               R2(J1,J2)=R2(J2,J1)
               R8(J1,J2)=R8(J2,J1)
               R14(J1,J2)=R14(J2,J1)
            ENDDO
         ENDDO
         DO J1=NTYPEA+1,N
            R2(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,N
               R2DUM=1.0D0/(VEC(J1,J2,1)**2+VEC(J1,J2,2)**2+VEC(J1,J2,3)**2)
               R2(J2,J1)=R2DUM
               R6=R2(J2,J1)**3
               IF (R2DUM.GT.IRCUT2BB) THEN
                  POTEL=POTEL + EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB) ! BB
C                 WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB),POTEL
C                 IF (EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB).GT.10.0D0) 
C    1            WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB),POTEL
               ENDIF
               R8(J2,J1)=R6*R2(J2,J1)
               R14(J2,J1)=R8(J2,J1)*R6
               R2(J1,J2)=R2(J2,J1)
               R8(J1,J2)=R8(J2,J1)
               R14(J1,J2)=R14(J2,J1)
            ENDDO
         ENDDO
C        WRITE(*,'(A,6F20.10)') 'X53,X93=,',X(157),X(158),X(159),X(277),X(278),X(279)
      ELSE
         DO J1=1,NTYPEA
            J3=3*(J1-1)
            DO J2=J1+1,NTYPEA
               J4=3*(J2-1)
               VEC1=X(J3+1)-X(J4+1)
               VEC2=X(J3+2)-X(J4+2)
               VEC3=X(J3+3)-X(J4+3)
               R2DUM=VEC1**2+VEC2**2+VEC3**2
               R2DUM=1.0D0/R2DUM
               R6=R2DUM**3
               IF (R2DUM.GT.IRCUT2AA) THEN
                  POTEL=POTEL + R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA  ! AA
C                 WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA,POTEL
C                 IF (R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA.GT.10.0D0) 
C    1            WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA,POTEL
               ENDIF
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            J3=3*(J1-1)
            DO J2=NTYPEA+1, N
               J4=3*(J2-1)
               VEC1=X(J3+1)-X(J4+1)
               VEC2=X(J3+2)-X(J4+2)
               VEC3=X(J3+3)-X(J4+3)
               R2DUM=VEC1**2+VEC2**2+VEC3**2
C              IF (J2.EQ.211) WRITE(*,'(A,I4,4F20.10)') 'J1,VEC12,VEC22,VEC32,R2DUM=',J1,VEC1**2,VEC2**2,VEC3**2,R2DUM
               R2DUM=1.0D0/R2DUM
               R6=R2DUM**3
               IF (R2DUM.GT.IRCUT2AB) THEN
                  POTEL=POTEL + EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB)  ! AB
C                 WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                       EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB),POTEL
C                 IF (EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB).GT.10.0D0)
C    1            WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                       EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB),POTEL
               ENDIF
            ENDDO
         ENDDO
         DO J1=NTYPEA+1, N
            J3=3*(J1-1)
            DO J2=J1+1, N
               J4=3*(J2-1)
               VEC1=X(J3+1)-X(J4+1)
               VEC2=X(J3+2)-X(J4+2)
               VEC3=X(J3+3)-X(J4+3)
               R2DUM=VEC1**2+VEC2**2+VEC3**2
C              IF (J2.EQ.211) WRITE(*,'(A,I4,4F20.10)') 'J1,VEC12,VEC22,VEC32,R2DUM=',J1,VEC1**2,VEC2**2,VEC3**2,R2DUM
               R2DUM=1.0D0/R2DUM
               R6=R2DUM**3
               IF (R2DUM.GT.IRCUT2BB) THEN
                  POTEL=POTEL + EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB) ! BB
C                 IF (EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB).GT.10.0D0)
C    1            WRITE(*,'(A,2I4,3F20.10)') 'J1,J2,R2DUM,PAIR,POTEL=',J1,J2,R2DUM,
C    1                 EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB),POTEL
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      POTEL = POTEL * 4.D0
     
      IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN
      CALL LJPSHIFTGBINC(G,R2,R14,R8,V,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12)
      
      IF (.NOT.STEST) RETURN
      CALL LJPSHIFTSBINC(G,F,R2,R14,R8,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12)

      RETURN
      END

C*****************************************************************************
  
      SUBROUTINE LJPSHIFTGBINC(G,R2,R14,R8,V,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, NTYPEA
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N), 
     1                 VEC(N,N,3), V(3*N), R2(N,N),
     2                 EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12,
     3                 RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB
      COMMON /RCONSTBIN/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB
C
C  Calculate the g tensor.
C
      DO J1=1,NTYPEA
         G(J1,J1)=0.0D0
         DO J2=J1+1,NTYPEA 
            G(J2,J1)=0.0D0
            IF (R2(J2,J1).GT.IRCUT2AA) G(J2,J1)=-8.0D0      *(3.0D0*(2.0D0*R14(J2,J1)        -R8(J2,J1)       )-RCONSTAA) 
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
      DO J1=1,NTYPEA
         DO J2=NTYPEA+1,N 
            G(J2,J1)=0.0D0
            IF (R2(J2,J1).GT.IRCUT2AB) G(J2,J1)=-8.0D0*EPSAB*(3.0D0*(2.0D0*R14(J2,J1)*SIGAB12-R8(J2,J1)*SIGAB6)-RCONSTAB)
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
      DO J1=NTYPEA+1,N
         G(J1,J1)=0.0D0
         DO J2=J1+1,N 
            G(J2,J1)=0.0D0
            IF (R2(J2,J1).GT.IRCUT2BB) G(J2,J1)=-8.0D0*EPSBB*(3.0D0*(2.0D0*R14(J2,J1)*SIGBB12-R8(J2,J1)*SIGBB6)-RCONSTBB)
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO J4=1,N
               IF (J1.LE.NTYPEA) THEN
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AA) V(J3)=V(J3)+G(J4,J1)*VEC(J4,J1,J2)
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2AB) V(J3)=V(J3)+G(J4,J1)*VEC(J4,J1,J2)
                  ENDIF
               ELSE
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AB) V(J3)=V(J3)+G(J4,J1)*VEC(J4,J1,J2)
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2BB) V(J3)=V(J3)+G(J4,J1)*VEC(J4,J1,J2)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
      RETURN
      END

C*****************************************************************************

      SUBROUTINE LJPSHIFTSBINC(G,F,R2,R14,R8,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6, NTYPEA
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N), IRCUT2AA, IRCUT2AB, IRCUT2BB, 
     1                 VEC(N,N,3), R2(N,N), RCONSTAA, RCONSTAB, RCONSTBB,
     2                 F(N,N), EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12
      COMMON /RCONSTBIN/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB

      DO J1=1,N
         F(J1,J1)=0.0D0
         DO J2=J1+1,N 
            F(J2,J1)=0.0D0
            IF (J1.LE.NTYPEA) THEN
               IF (J2.LE.NTYPEA) THEN
                  IF (R2(J2,J1).GT.IRCUT2AA) F(J2,J1)=8.0D0*      (84.0D0*R14(J2,J1)-24.0D0*R8(J2,J1))  
               ELSE
                  IF (R2(J2,J1).GT.IRCUT2AB) F(J2,J1)=8.0D0*EPSAB*(84.0D0*R14(J2,J1)*SIGAB12-24.0D0*R8(J2,J1)*SIGAB6)
               ENDIF
            ELSE
               IF (R2(J2,J1).GT.IRCUT2BB) F(J2,J1)=   8.0D0*EPSBB*(84.0D0*R14(J2,J1)*SIGBB12-24.0D0*R8(J2,J1)*SIGBB6)
            ENDIF
C           F(J2,J1)=672.0D0*R14(J2,J1)-192.0D0*R8(J2,J1)
            F(J1,J2)=F(J2,J1)
         ENDDO
      ENDDO
C
C  Now do the hessian. First are the entirely diagonal terms.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            HESS(J3,J3)=0.0D0
            DO J4=1,N
               IF (J1.LE.NTYPEA) THEN
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AA) HESS(J3,J3)=HESS(J3,J3)+F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2 + G(J4,J1)
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2AB) HESS(J3,J3)=HESS(J3,J3)+F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2 + G(J4,J1)
                  ENDIF
               ELSE
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AB) HESS(J3,J3)=HESS(J3,J3)+F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2 + G(J4,J1)
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2BB) HESS(J3,J3)=HESS(J3,J3)+F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2 + G(J4,J1)
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C
C  Next are the terms where x_i and x_j are on the same atom
C  but are different, e.g. y and z.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J2+1,3
               HESS(3*(J1-1)+J4,J3)=0.0D0
               DO J5=1,N
                  IF (J1.LE.NTYPEA) THEN
                     IF (J5.LE.NTYPEA) THEN
                        IF (R2(J5,J1).GT.IRCUT2AA) 
     1                      HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*VEC(J5,J1,J2)*VEC(J5,J1,J4)
                     ELSE
                        IF (R2(J5,J1).GT.IRCUT2AB) 
     1                      HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*VEC(J5,J1,J2)*VEC(J5,J1,J4)
                     ENDIF
                  ELSE
                     IF (J5.LE.NTYPEA) THEN
                        IF (R2(J5,J1).GT.IRCUT2AB)
     1                      HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*VEC(J5,J1,J2)*VEC(J5,J1,J4)
                     ELSE
                        IF (R2(J5,J1).GT.IRCUT2BB) 
     1                      HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*VEC(J5,J1,J2)*VEC(J5,J1,J4)
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C  Case III, different atoms, same cartesian coordinate.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               IF (J1.LE.NTYPEA) THEN
                  IF (J4.LE.NTYPEA) THEN
                     IF (R2(J4,J1).GT.IRCUT2AA) THEN
                        HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2-G(J4,J1)
                     ELSE
                        HESS(3*(J4-1)+J2,J3)=0.0D0
                     ENDIF
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2AB) THEN
                        HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2-G(J4,J1)
                     ELSE
                        HESS(3*(J4-1)+J2,J3)=0.0D0
                     ENDIF
                  ENDIF
               ELSE
                  IF (R2(J4,J1).GT.IRCUT2BB) THEN
                     HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2-G(J4,J1)
                  ELSE
                     HESS(3*(J4-1)+J2,J3)=0.0D0
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C
C  Case IV: different atoms and different cartesian coordinates.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               DO J5=1,J2-1
                  J6=3*(J4-1)+J5
                  IF (J1.LE.NTYPEA) THEN
                     IF (J4.LE.NTYPEA) THEN
                        IF (R2(J4,J1).GT.IRCUT2AA) THEN
                           HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)*VEC(J4,J1,J5)
                        ELSE
                           HESS(J6,J3)=0.0D0
                        ENDIF
                     ELSE
                        IF (R2(J4,J1).GT.IRCUT2AB) THEN
                           HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)*VEC(J4,J1,J5)
                        ELSE
                           HESS(J6,J3)=0.0D0
                        ENDIF
                     ENDIF
                  ELSE
                     IF (R2(J4,J1).GT.IRCUT2BB) THEN
                        HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)*VEC(J4,J1,J5)
                     ELSE
                        HESS(J6,J3)=0.0D0
                     ENDIF
                  ENDIF
                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C  Symmetrise Hessian
C
      DO J1=1,3*N
         DO J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
         ENDDO
      ENDDO

      RETURN
      END
C
C*************************************************************************
C
C  Subroutine LJPSHIFT calculates the energy, cartesian gradient and second
C  derivative matrix analytically for Lennard-Jones in reduced units
C  (epsilon=sigma=1) using a shifted, truncated potential.
C
C  Adapted for the binary LJ glass described by Sastry, Debenetti and
C  Stillinger, Nature, 393, 554, 1998. Atom types are A and B. The first
C  NTYPEA are A, the next NBTYPE=NATOMS-NTYPEA are B. epsilon and sigma for A are the
C  units of energy and distance, so we also need EPSAB, EPSAB, SIGAB and
C  SIGAA in these units. Sastry et al. density is 1.2 i.e. a box length
C  of 5.975206 for 256 atoms. 
C
C
C*************************************************************************
C
      SUBROUTINE LJPSHIFT(N, X, V, POTEL, BOXLX, BOXLY, BOXLZ, CUTOFF, GTEST, STEST, PTEST)
      USE KEY
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, ANV(N,N,3)
      DOUBLE PRECISION X(3*N), VEC1, VEC2, VEC3,
     1                 V(3*N), R2(N,N),  R6, R2DUM,
     2                 R8(N,N), G(N,N), 
     3                 R14(N,N), F(N,N), 
     4                 POTEL, BOXLX, BOXLY, BOXLZ, 
     5                 VEC(N,N,3), IRCUT2AA, CUTOFF, SIGRCAA6, SIGRCAA12,
     6                 CONSTAA, RCONSTAA, CUTAA
      LOGICAL GTEST, STEST, PTEST
      COMMON /RCONST/ RCONSTAA, IRCUT2AA
C     COMMON /ANVCOMMON/ ANV

      CUTAA=CUTOFF
      IRCUT2AA = 1.D0/CUTAA**2
      SIGRCAA6= 1.0D0/CUTAA**6
      SIGRCAA12=SIGRCAA6**2
      CONSTAA=4.0D0*SIGRCAA6-7.0D0*SIGRCAA12
      RCONSTAA=(6.0D0*SIGRCAA12-3.0D0*SIGRCAA6)/CUTAA**2
C
C  Work out cutoff for potential. Two particles interact if r<c, but
C  we will use the equivalent condition 1/r^2 > 1/c^2.
C
C     IF (PTEST) PRINT*,'Cutoff used = ',CUTOFF
C
C  Deal with any atoms that have left the box.
C
      IF ((.NOT.FIXIMAGE).AND.(.NOT.NORESET)) THEN
         DO J1=1,N
            J2 = 3*(J1-1)
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
      IF (GTEST) THEN
         DO J1=1, N
            VEC(J1,J1,1)=0.0D0
            VEC(J1,J1,2)=0.0D0
            VEC(J1,J1,3)=0.0D0
            J3=3*(J1-1)
            DO J2=J1+1, N
               J4=3*(J2-1)
               VEC(J2,J1,1)=X(J3+1)-X(J4+1)
               VEC(J2,J1,2)=X(J3+2)-X(J4+2)
               VEC(J2,J1,3)=X(J3+3)-X(J4+3)
               IF (.NOT.FIXIMAGE) THEN
                  ANV(J2,J1,1)=NINT(VEC(J2,J1,1)/BOXLX)
                  ANV(J2,J1,2)=NINT(VEC(J2,J1,2)/BOXLY)
                  ANV(J2,J1,3)=NINT(VEC(J2,J1,3)/BOXLZ)
               ENDIF
               VEC(J2,J1,1)=VEC(J2,J1,1)-BOXLX*ANV(J2,J1,1)
               VEC(J2,J1,2)=VEC(J2,J1,2)-BOXLY*ANV(J2,J1,2)
               VEC(J2,J1,3)=VEC(J2,J1,3)-BOXLZ*ANV(J2,J1,3)
               VEC(J1,J2,1)=-VEC(J2,J1,1)
               VEC(J1,J2,2)=-VEC(J2,J1,2)
               VEC(J1,J2,3)=-VEC(J2,J1,3)
            ENDDO
         ENDDO
         DO J1=1,N
            R2(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,N
               IF (DOUBLET.AND.(J1.EQ.1).AND.(J2.EQ.2)) THEN
                  R2(J2,J1)=0.0D0 
                  R6=0.0D0
               ELSE
                  R2DUM=1.0D0/(VEC(J1,J2,1)**2+VEC(J1,J2,2)**2+VEC(J1,J2,3)**2)
                  R2(J2,J1)=R2DUM
                  R6=R2(J2,J1)**3
                  IF (R2DUM.GT.IRCUT2AA) THEN
                     POTEL=POTEL + R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA
                  ENDIF
               ENDIF
               R8(J2,J1)=R6*R2(J2,J1)
               R14(J2,J1)=R8(J2,J1)*R6
               R2(J1,J2)=R2(J2,J1)
               R8(J1,J2)=R8(J2,J1)
               R14(J1,J2)=R14(J2,J1)
            ENDDO
         ENDDO
      ELSE
         DO J1=1,N
            J3=3*(J1-1)
            DO J2=J1+1,N
               J4=3*(J2-1)
               VEC1=X(J3+1)-X(J4+1)
               VEC2=X(J3+2)-X(J4+2)
               VEC3=X(J3+3)-X(J4+3)
               IF (.NOT.FIXIMAGE) THEN
                  ANV(J2,J1,1)=NINT(VEC1/BOXLX)
                  ANV(J2,J1,2)=NINT(VEC2/BOXLY)
                  ANV(J2,J1,3)=NINT(VEC3/BOXLZ)
               ENDIF
               VEC1=VEC1-BOXLX*ANV(J2,J1,1)
               VEC2=VEC2-BOXLY*ANV(J2,J1,2)
               VEC3=VEC3-BOXLZ*ANV(J2,J1,3)
               R2DUM=VEC1**2+VEC2**2+VEC3**2
               R2DUM=1.0D0/R2DUM
               R6=R2DUM**3
               IF (.NOT.(DOUBLET.AND.(J1.EQ.1).AND.(J2.EQ.2))) THEN
                  IF (R2DUM.GT.IRCUT2AA) THEN
                     POTEL=POTEL + R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA 
                  ENDIF
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      POTEL = POTEL * 4.D0
     
      IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN
      CALL LJPSHIFTG(G,R2,R14,R8,V,VEC,N,DOUBLET)
      
      IF (.NOT.STEST) RETURN
      CALL LJPSHIFTS(G,F,R2,R14,R8,VEC,N)

      RETURN
      END

C*****************************************************************************
  
      SUBROUTINE LJPSHIFTG(G,R2,R14,R8,V,VEC,N,DOUBLET)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N),
     1                 VEC(N,N,3), V(3*N), R2(N,N),
     3                 RCONSTAA, IRCUT2AA
      LOGICAL DOUBLET
      COMMON /RCONST/ RCONSTAA, IRCUT2AA
C
C  Calculate the g tensor.
C
      DO J1=1,N
         G(J1,J1)=0.0D0
         DO J2=J1+1,N
            G(J2,J1)=0.0D0
            IF (R2(J2,J1).GT.IRCUT2AA) G(J2,J1)=-8.0D0      *(3.0D0*(2.0D0*R14(J2,J1)        -R8(J2,J1)       )-RCONSTAA) 
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
      IF (DOUBLET) THEN
         G(1,2)=0.0D0
         G(2,1)=0.0D0
      ENDIF
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO J4=1,N
               IF (R2(J4,J1).GT.IRCUT2AA) V(J3)=V(J3)+G(J4,J1)*VEC(J4,J1,J2)
            ENDDO
         ENDDO
      ENDDO

      RETURN
      END

C*****************************************************************************

      SUBROUTINE LJPSHIFTS(G,F,R2,R14,R8,VEC,N)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N), IRCUT2AA, 
     1                 VEC(N,N,3), R2(N,N), RCONSTAA, 
     2                 F(N,N)
      COMMON /RCONST/ RCONSTAA, IRCUT2AA

      DO J1=1,N
         F(J1,J1)=0.0D0
         DO J2=J1+1,N 
            F(J2,J1)=0.0D0
            IF (R2(J2,J1).GT.IRCUT2AA) F(J2,J1)=8.0D0*      (84.0D0*R14(J2,J1)-24.0D0*R8(J2,J1))  
            F(J1,J2)=F(J2,J1)
         ENDDO
      ENDDO
C
C  Now do the hessian. First are the entirely diagonal terms.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            HESS(J3,J3)=0.0D0
            DO J4=1,N
               IF (R2(J4,J1).GT.IRCUT2AA) HESS(J3,J3)=HESS(J3,J3)+F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2 + G(J4,J1)
            ENDDO
         ENDDO
      ENDDO
C
C  Next are the terms where x_i and x_j are on the same atom
C  but are different, e.g. y and z.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J2+1,3
               HESS(3*(J1-1)+J4,J3)=0.0D0
               DO J5=1,N
                  IF (R2(J5,J1).GT.IRCUT2AA) 
     1          HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*VEC(J5,J1,J2)*VEC(J5,J1,J4)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C  Case III, different atoms, same cartesian coordinate.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               IF (R2(J4,J1).GT.IRCUT2AA) THEN
                  HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2-G(J4,J1)
               ELSE
                  HESS(3*(J4-1)+J2,J3)=0.0D0
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C
C  Case IV: different atoms and different cartesian coordinates.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               DO J5=1,J2-1
                  J6=3*(J4-1)+J5
                  IF (R2(J4,J1).GT.IRCUT2AA) THEN
                     HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)*VEC(J4,J1,J5)
                  ELSE
                     HESS(J6,J3)=0.0D0
                  ENDIF
                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C  Symmetrise Hessian
C
      DO J1=1,3*N
         DO J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
         ENDDO
      ENDDO

      RETURN
      END

      SUBROUTINE LJPSLATMIN(N,XC,BOXLX,BOXLY,BOXLZ,CUTOFF,V)
      IMPLICIT NONE
      INTEGER N, J1, NCOUNT
      DOUBLE PRECISION XC(3*N),F1,F2,F3,GRAD,SECOND,
     1                 XSAVE(3*N),BSAVE,V(3*N),
     2                 BOXLX,BOXLY,BOXLZ,CUTOFF,DIF,TEMP1

C
C  Value of DIF is the order of magnitude to which the lattice
C  constant can be optimised. Setting it smaller than 10^(-7)
C  causes numerical problems on the DEC.
C
      DIF=1.0D-7
      BSAVE=BOXLX
      NCOUNT=1
      DO 20 J1=1,3*N
         XSAVE(J1)=XC(J1)
20    CONTINUE

10    TEMP1=BOXLX
      BOXLX=TEMP1+DIF
      BOXLY=TEMP1+DIF
      BOXLZ=TEMP1+DIF

      CALL LJPSHIFT(N, XC,  V, F1,   BOXLX, BOXLY, BOXLZ, CUTOFF, .FALSE., .FALSE., .FALSE.)

      BOXLX=TEMP1-DIF
      BOXLY=TEMP1-DIF
      BOXLZ=TEMP1-DIF
      DO 30 J1=1,3*N
         XC(J1)=XSAVE(J1)
30    CONTINUE
      CALL LJPSHIFT(N, XC,  V, F2,   BOXLX, BOXLY, BOXLZ, CUTOFF, .FALSE., .FALSE., .FALSE.)

      GRAD=(F1-F2)/(2.0D0*DIF)

      BOXLX=TEMP1
      BOXLY=TEMP1
      BOXLZ=TEMP1
      DO 40 J1=1,3*N
         XC(J1)=XSAVE(J1)
40    CONTINUE
      CALL LJPSHIFT(N, XC,  V, F1,   BOXLX, BOXLY, BOXLZ, CUTOFF, .FALSE., .FALSE., .FALSE.)

      BOXLX=TEMP1+2.0D0*DIF
      BOXLY=TEMP1+2.0D0*DIF
      BOXLZ=TEMP1+2.0D0*DIF
      DO 50 J1=1,3*N
         XC(J1)=XSAVE(J1)
50    CONTINUE
      CALL LJPSHIFT(N, XC,  V, F2,   BOXLX, BOXLY, BOXLZ, CUTOFF, .FALSE., .FALSE., .FALSE.)

      BOXLX=TEMP1-2.0D0*DIF
      BOXLY=TEMP1-2.0D0*DIF
      BOXLZ=TEMP1-2.0D0*DIF
      DO 60 J1=1,3*N
         XC(J1)=XSAVE(J1)
60    CONTINUE
      CALL LJPSHIFT(N, XC,  V, F3,   BOXLX, BOXLY, BOXLZ, CUTOFF, .FALSE., .FALSE., .FALSE.)

      SECOND=(F3+F2-2.0D0*F1)/(4.0D0*DIF*DIF)
      PRINT*,'Energy for lattice cycle ',NCOUNT,' is ',F1
      PRINT*,'Gradient wrt box length=',GRAD
      PRINT*,'Second derivative wrt box length=',SECOND
      NCOUNT=NCOUNT+1
      IF (DABS(GRAD/SECOND).GT.1.0D-4) THEN
         BOXLX=BOXLX-0.0001D0*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
         BOXLY=BOXLY-0.0001D0*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
         BOXLZ=BOXLZ-0.0001D0*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
         PRINT*,'Step=',-0.0001D0*GRAD*DABS(SECOND)/(SECOND*DABS(GRAD))
         GOTO 10
      ELSE
         BOXLX=BOXLX-GRAD/SECOND
         BOXLY=BOXLY-GRAD/SECOND
         BOXLZ=BOXLZ-GRAD/SECOND
         PRINT*,'Step=',-GRAD/SECOND
         IF (DABS(GRAD/SECOND).GT.1.0D-6) GOTO 10
      ENDIF
      RETURN
      END
C
C*************************************************************************
C
C  Subroutine LJPSHIFTBIN2 calculates the energy, cartesian gradient and second
C  derivative matrix analytically for Lennard-Jones in reduced units
C  (epsilon=sigma=1) using a shifted, truncated potential.
C
C  Adapted for the binary LJ glass described by Sastry, Debenetti and
C  Stillinger, Nature, 393, 554, 1998. Atom types are A and B. The first
C  NTYPEA are A, the next NBTYPE=NATOMS-NTYPEA are B. epsilon and sigma for A are the
C  units of energy and distance, so we also need EPSAB, EPSAB, SIGAB and
C  SIGAA in these units. Sastry et al. density is 1.2 i.e. a box length
C  of 5.975206 for 256 atoms. 
C
C  This routine uses the same neighbours when FIXIMAGE is set to try and avoid
C  discontinuities in the second derivatives.
C
C
C*************************************************************************
C
      SUBROUTINE LJPSHIFTBIN2(N, X, V, POTEL, BOXLX, BOXLY, BOXLZ, CUTOFF, GTEST, STEST, PTEST)
      USE KEY
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, NTYPEA, ANV(N,N,3)
      DOUBLE PRECISION X(3*N), VEC1, VEC2, VEC3,
     1                 V(3*N), R2(N,N),  R6, R2DUM,
     2                 R8(N,N), G(N,N), EPSAB, EPSBB, SIGAB, SIGBB,
     3                 R14(N,N), F(N,N), 
     4                 POTEL, BOXLX, BOXLY, BOXLZ, SIGAB6, SIGAB12, SIGBB6, SIGBB12,
     5                 VEC(N,N,3), IRCUT2AA, IRCUT2AB, IRCUT2BB, CUTOFF, SIGRCAA6, SIGRCAA12,
     6                 SIGRCAB6, SIGRCAB12, SIGRCBB6, SIGRCBB12, CONSTAA, CONSTBB, CONSTAB,
     7                 RCONSTAA, RCONSTAB, RCONSTBB, CUTAA, CUTAB, CUTBB
      LOGICAL GTEST, STEST, PTEST, CHECKCUT(N,N)
      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      COMMON /RCONSTBIN/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB
C     COMMON /ANVCOMMON/ ANV
C     SAVE CHECKCUT

      CUTAA=CUTOFF
      CUTAB=CUTOFF*SIGAB
      CUTBB=CUTOFF*SIGBB
      IRCUT2AA = 1.D0/CUTAA**2
      IRCUT2AB = 1.D0/CUTAB**2
      IRCUT2BB = 1.D0/CUTBB**2
      SIGAB6=SIGAB**6
      SIGAB12=SIGAB6**2
      SIGBB6=SIGBB**6
      SIGBB12=SIGBB6**2
      SIGRCAA6= 1.0D0/CUTAA**6
      SIGRCAA12=SIGRCAA6**2
      SIGRCAB6=SIGAB6/CUTAB**6
      SIGRCAB12=SIGRCAB6**2
      SIGRCBB6=SIGBB6/CUTBB**6
      SIGRCBB12=SIGRCBB6**2
      CONSTAA=4.0D0*SIGRCAA6-7.0D0*SIGRCAA12
      CONSTAB=4.0D0*SIGRCAB6-7.0D0*SIGRCAB12
      CONSTBB=4.0D0*SIGRCBB6-7.0D0*SIGRCBB12
      RCONSTAA=(6.0D0*SIGRCAA12-3.0D0*SIGRCAA6)/CUTAA**2
      RCONSTAB=(6.0D0*SIGRCAB12-3.0D0*SIGRCAB6)/CUTAB**2
      RCONSTBB=(6.0D0*SIGRCBB12-3.0D0*SIGRCBB6)/CUTBB**2
      DO J1=1,N
         CHECKCUT(J1,J1)=.FALSE.
      ENDDO
C
C  Work out cutoff for potential. Two particles interact if r<c, but
C  we will use the equivalent condition 1/r^2 > 1/c^2.
C
C     IF (PTEST) PRINT*,'Cutoff used = ',CUTOFF
C
C  Deal with any atoms that have left the box.
C
C     PRINT*,'in LJPSHIFTBIN2 FIXIMAGE,NORESET=',FIXIMAGE,NORESET
      IF ((.NOT.FIXIMAGE).AND.(.NOT.NORESET)) THEN
         DO J1=1,N
            J2 = 3*(J1-1)
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
      IF (GTEST) THEN
         DO J1=1, N
            VEC(J1,J1,1)=0.0D0
            VEC(J1,J1,2)=0.0D0
            VEC(J1,J1,3)=0.0D0
            J3=3*(J1-1)
            DO J2=J1+1, N
               J4=3*(J2-1)
               VEC(J2,J1,1)=X(J3+1)-X(J4+1)
               VEC(J2,J1,2)=X(J3+2)-X(J4+2)
               VEC(J2,J1,3)=X(J3+3)-X(J4+3)
               IF (.NOT.FIXIMAGE) THEN
                  ANV(J2,J1,1)=NINT(VEC(J2,J1,1)/BOXLX)
                  ANV(J2,J1,2)=NINT(VEC(J2,J1,2)/BOXLY)
                  ANV(J2,J1,3)=NINT(VEC(J2,J1,3)/BOXLZ)
               ENDIF
               VEC(J2,J1,1)=VEC(J2,J1,1)-BOXLX*ANV(J2,J1,1)
               VEC(J2,J1,2)=VEC(J2,J1,2)-BOXLY*ANV(J2,J1,2)
               VEC(J2,J1,3)=VEC(J2,J1,3)-BOXLZ*ANV(J2,J1,3)
               VEC(J1,J2,1)=-VEC(J2,J1,1)
               VEC(J1,J2,2)=-VEC(J2,J1,2)
               VEC(J1,J2,3)=-VEC(J2,J1,3)
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            R2(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,NTYPEA
               R2DUM=1.0D0/(VEC(J1,J2,1)**2+VEC(J1,J2,2)**2+VEC(J1,J2,3)**2)
               IF (.NOT.FIXIMAGE) THEN
                  CHECKCUT(J2,J1)=.FALSE.
                  IF (R2DUM.GT.IRCUT2AA) CHECKCUT(J2,J1)=.TRUE.
                  CHECKCUT(J1,J2)=CHECKCUT(J2,J1)
               ENDIF
               R2(J2,J1)=R2DUM
               R6=R2(J2,J1)**3
               IF (CHECKCUT(J2,J1)) POTEL=POTEL + R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA  ! AA
               R8(J2,J1)=R6*R2(J2,J1)
               R14(J2,J1)=R8(J2,J1)*R6
               R2(J1,J2)=R2(J2,J1)
               R8(J1,J2)=R8(J2,J1)
               R14(J1,J2)=R14(J2,J1)
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            DO J2=NTYPEA+1,N
               R2DUM=1.0D0/(VEC(J1,J2,1)**2+VEC(J1,J2,2)**2+VEC(J1,J2,3)**2)
               IF (.NOT.FIXIMAGE) THEN
                  CHECKCUT(J2,J1)=.FALSE.
                  IF (R2DUM.GT.IRCUT2AB) CHECKCUT(J2,J1)=.TRUE.
                  CHECKCUT(J1,J2)=CHECKCUT(J2,J1)
               ENDIF
               R2(J2,J1)=R2DUM
               R6=R2(J2,J1)**3
               IF (CHECKCUT(J2,J1)) POTEL=POTEL + EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB) ! AB
               R8(J2,J1)=R6*R2(J2,J1)
               R14(J2,J1)=R8(J2,J1)*R6
               R2(J1,J2)=R2(J2,J1)
               R8(J1,J2)=R8(J2,J1)
               R14(J1,J2)=R14(J2,J1)
            ENDDO
         ENDDO
         DO J1=NTYPEA+1,N
            R2(J1,J1)=0.0D0
            R8(J1,J1)=0.0D0
            R14(J1,J1)=0.0D0
            DO J2=J1+1,N
               R2DUM=1.0D0/(VEC(J1,J2,1)**2+VEC(J1,J2,2)**2+VEC(J1,J2,3)**2)
               IF (.NOT.FIXIMAGE) THEN
                  CHECKCUT(J2,J1)=.FALSE.
                  IF (R2DUM.GT.IRCUT2BB) CHECKCUT(J2,J1)=.TRUE.
                  CHECKCUT(J1,J2)=CHECKCUT(J2,J1)
               ENDIF
               R2(J2,J1)=R2DUM
               R6=R2(J2,J1)**3
               IF (CHECKCUT(J2,J1)) POTEL=POTEL + EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB) ! BB
               R8(J2,J1)=R6*R2(J2,J1)
               R14(J2,J1)=R8(J2,J1)*R6
               R2(J1,J2)=R2(J2,J1)
               R8(J1,J2)=R8(J2,J1)
               R14(J1,J2)=R14(J2,J1)
            ENDDO
         ENDDO
      ELSE
         DO J1=1,NTYPEA
            J3=3*(J1-1)
            DO J2=J1+1,NTYPEA
               J4=3*(J2-1)
               VEC1=X(J3+1)-X(J4+1)
               VEC2=X(J3+2)-X(J4+2)
               VEC3=X(J3+3)-X(J4+3)
               IF (.NOT.FIXIMAGE) THEN
                  ANV(J2,J1,1)=NINT(VEC1/BOXLX)
                  ANV(J2,J1,2)=NINT(VEC2/BOXLY)
                  ANV(J2,J1,3)=NINT(VEC3/BOXLZ)
               ENDIF
               VEC1=VEC1-BOXLX*ANV(J2,J1,1)
               VEC2=VEC2-BOXLY*ANV(J2,J1,2)
               VEC3=VEC3-BOXLZ*ANV(J2,J1,3)
               R2DUM=VEC1**2+VEC2**2+VEC3**2
               R2DUM=1.0D0/R2DUM
               IF (.NOT.FIXIMAGE) THEN
                  CHECKCUT(J2,J1)=.FALSE.
                  IF (R2DUM.GT.IRCUT2AA) CHECKCUT(J2,J1)=.TRUE.
                  CHECKCUT(J1,J2)=CHECKCUT(J2,J1)
               ENDIF
               R6=R2DUM**3
               IF (CHECKCUT(J2,J1)) POTEL=POTEL + R6*(R6-1.0D0) + RCONSTAA/R2DUM + CONSTAA  ! AA
            ENDDO
         ENDDO
         DO J1=1,NTYPEA
            J3=3*(J1-1)
            DO J2=NTYPEA+1, N
               J4=3*(J2-1)
               VEC1=X(J3+1)-X(J4+1)
               VEC2=X(J3+2)-X(J4+2)
               VEC3=X(J3+3)-X(J4+3)
               IF (.NOT.FIXIMAGE) THEN
                  ANV(J2,J1,1)=NINT(VEC1/BOXLX)
                  ANV(J2,J1,2)=NINT(VEC2/BOXLY)
                  ANV(J2,J1,3)=NINT(VEC3/BOXLZ)
               ENDIF
               VEC1=VEC1-BOXLX*ANV(J2,J1,1)
               VEC2=VEC2-BOXLY*ANV(J2,J1,2)
               VEC3=VEC3-BOXLZ*ANV(J2,J1,3)
               R2DUM=VEC1**2+VEC2**2+VEC3**2
               R2DUM=1.0D0/R2DUM
               IF (.NOT.FIXIMAGE) THEN
                  CHECKCUT(J2,J1)=.FALSE.
                  IF (R2DUM.GT.IRCUT2AB) CHECKCUT(J2,J1)=.TRUE.
                  CHECKCUT(J1,J2)=CHECKCUT(J2,J1)
               ENDIF
               R6=R2DUM**3
               IF (CHECKCUT(J2,J1)) POTEL=POTEL + EPSAB*(SIGAB6*R6*(SIGAB6*R6-1.0D0) + RCONSTAB/R2DUM + CONSTAB)  ! AB
            ENDDO
         ENDDO
         DO J1=NTYPEA+1, N
            J3=3*(J1-1)
            DO J2=J1+1, N
               J4=3*(J2-1)
               VEC1=X(J3+1)-X(J4+1)
               VEC2=X(J3+2)-X(J4+2)
               VEC3=X(J3+3)-X(J4+3)
               IF (.NOT.FIXIMAGE) THEN
                  ANV(J2,J1,1)=NINT(VEC1/BOXLX)
                  ANV(J2,J1,2)=NINT(VEC2/BOXLY)
                  ANV(J2,J1,3)=NINT(VEC3/BOXLZ)
               ENDIF
               VEC1=VEC1-BOXLX*ANV(J2,J1,1)
               VEC2=VEC2-BOXLY*ANV(J2,J1,2)
               VEC3=VEC3-BOXLZ*ANV(J2,J1,3)
               R2DUM=VEC1**2+VEC2**2+VEC3**2
               R2DUM=1.0D0/R2DUM
               IF (.NOT.FIXIMAGE) THEN
                  CHECKCUT(J2,J1)=.FALSE.
                  IF (R2DUM.GT.IRCUT2BB) CHECKCUT(J2,J1)=.TRUE.
                  CHECKCUT(J1,J2)=CHECKCUT(J2,J1)
               ENDIF
               R6=R2DUM**3
               IF (CHECKCUT(J2,J1)) POTEL=POTEL + EPSBB*(SIGBB6*R6*(SIGBB6*R6-1.0D0) + RCONSTBB/R2DUM + CONSTBB) ! BB
            ENDDO
         ENDDO
      ENDIF
      POTEL = POTEL * 4.D0
     
      IF ((.NOT.GTEST).AND.(.NOT.STEST)) RETURN
      CALL LJPSHIFTGBIN2(G,R2,R14,R8,V,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12,CHECKCUT)
      
      IF (.NOT.STEST) RETURN
      CALL LJPSHIFTSBIN2(G,F,R2,R14,R8,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12,CHECKCUT)

      RETURN
      END

C*****************************************************************************
  
      SUBROUTINE LJPSHIFTGBIN2(G,R2,R14,R8,V,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12,CHECKCUT)
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, NTYPEA
      LOGICAL CHECKCUT(N,N)
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N),
     1                 VEC(N,N,3), V(3*N), R2(N,N),
     2                 EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12,
     3                 RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB
      COMMON /RCONSTBIN/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB
C
C  Calculate the g tensor.
C
      DO J1=1,NTYPEA
         G(J1,J1)=0.0D0
         DO J2=J1+1,NTYPEA 
            G(J2,J1)=0.0D0
            IF (CHECKCUT(J2,J1)) G(J2,J1)=-8.0D0      *(3.0D0*(2.0D0*R14(J2,J1)        -R8(J2,J1)       )-RCONSTAA) 
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
      DO J1=1,NTYPEA
         DO J2=NTYPEA+1,N 
            G(J2,J1)=0.0D0
            IF (CHECKCUT(J2,J1)) G(J2,J1)=-8.0D0*EPSAB*(3.0D0*(2.0D0*R14(J2,J1)*SIGAB12-R8(J2,J1)*SIGAB6)-RCONSTAB)
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
      DO J1=NTYPEA+1,N
         G(J1,J1)=0.0D0
         DO J2=J1+1,N 
            G(J2,J1)=0.0D0
            IF (CHECKCUT(J2,J1)) G(J2,J1)=-8.0D0*EPSBB*(3.0D0*(2.0D0*R14(J2,J1)*SIGBB12-R8(J2,J1)*SIGBB6)-RCONSTBB)
            G(J1,J2)=G(J2,J1)
         ENDDO
      ENDDO
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            V(J3)=0.0D0
            DO J4=1,N
               IF (CHECKCUT(J4,J1)) V(J3)=V(J3)+G(J4,J1)*VEC(J4,J1,J2)
            ENDDO
         ENDDO
      ENDDO
C     DO J1=1,10
C        WRITE(*,'(A,I6,F20.10)') 'J1,V=',J1,V(J1)
C     ENDDO

      RETURN
      END

C*****************************************************************************

      SUBROUTINE LJPSHIFTSBIN2(G,F,R2,R14,R8,VEC,N,NTYPEA,EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12,CHECKCUT)
      USE MODHESS
      IMPLICIT NONE
      INTEGER N, J1, J2, J3, J4, J5, J6, NTYPEA
      LOGICAL CHECKCUT(N,N)
      DOUBLE PRECISION G(N,N), R14(N,N), R8(N,N), IRCUT2AA, IRCUT2AB, IRCUT2BB,
     1                 VEC(N,N,3), R2(N,N), RCONSTAA, RCONSTAB, RCONSTBB,
     2                 F(N,N), EPSAB,EPSBB,SIGAB6,SIGBB6,SIGAB12,SIGBB12
      COMMON /RCONSTBIN/ RCONSTAA, RCONSTAB, RCONSTBB, IRCUT2AA, IRCUT2AB, IRCUT2BB

C     DO J1=1,N
C        DO J2=1,N
C           WRITE(*,'(A,2I6,L)') 'J1,J2,CHECK=',J1,J2,CHECKCUT(J1,J2)
C        ENDDO
C     ENDDO

      DO J1=1,N
         F(J1,J1)=0.0D0
         DO J2=J1+1,N 
            F(J2,J1)=0.0D0
            IF (J1.LE.NTYPEA) THEN
               IF (J2.LE.NTYPEA) THEN
                  IF (CHECKCUT(J2,J1)) F(J2,J1)=8.0D0*      (84.0D0*R14(J2,J1)-24.0D0*R8(J2,J1))  
               ELSE
                  IF (CHECKCUT(J2,J1)) F(J2,J1)=8.0D0*EPSAB*(84.0D0*R14(J2,J1)*SIGAB12-24.0D0*R8(J2,J1)*SIGAB6)
               ENDIF
            ELSE
               IF (CHECKCUT(J2,J1)) F(J2,J1)=   8.0D0*EPSBB*(84.0D0*R14(J2,J1)*SIGBB12-24.0D0*R8(J2,J1)*SIGBB6)
            ENDIF
C           F(J2,J1)=672.0D0*R14(J2,J1)-192.0D0*R8(J2,J1)
            F(J1,J2)=F(J2,J1)
         ENDDO
      ENDDO
C
C  Now do the hessian. First are the entirely diagonal terms.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            HESS(J3,J3)=0.0D0
            DO J4=1,N
               IF (CHECKCUT(J4,J1)) HESS(J3,J3)=HESS(J3,J3)+F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2 + G(J4,J1)
            ENDDO
         ENDDO
      ENDDO
C
C  Next are the terms where x_i and x_j are on the same atom
C  but are different, e.g. y and z.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J2+1,3
               HESS(3*(J1-1)+J4,J3)=0.0D0
               DO J5=1,N
                  IF (CHECKCUT(J5,J1)) HESS(3*(J1-1)+J4,J3)=HESS(3*(J1-1)+J4,J3)+F(J5,J1)*R2(J5,J1)*VEC(J5,J1,J2)*VEC(J5,J1,J4)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C  Case III, different atoms, same cartesian coordinate.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               IF (CHECKCUT(J4,J1)) THEN
                  HESS(3*(J4-1)+J2,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)**2-G(J4,J1)
               ELSE
                  HESS(3*(J4-1)+J2,J3)=0.0D0
               ENDIF
            ENDDO
         ENDDO
      ENDDO
C
C  Case IV: different atoms and different cartesian coordinates.
C
      DO J1=1,N
         DO J2=1,3
            J3=3*(J1-1)+J2
            DO J4=J1+1,N
               DO J5=1,J2-1
                  J6=3*(J4-1)+J5
                  IF (CHECKCUT(J4,J1)) THEN
                     HESS(J6,J3)=-F(J4,J1)*R2(J4,J1)*VEC(J4,J1,J2)*VEC(J4,J1,J5)
                  ELSE
                     HESS(J6,J3)=0.0D0
                  ENDIF
                  HESS(3*(J4-1)+J2,3*(J1-1)+J5)=HESS(J6,J3)
               ENDDO
            ENDDO
         ENDDO
      ENDDO
C
C  Symmetrise Hessian
C
      DO J1=1,3*N
         DO J2=J1+1,3*N
            HESS(J1,J2)=HESS(J2,J1)
         ENDDO
      ENDDO

C     DO J1=1,10
C        DO J2=J1,10
C           WRITE(*,'(A,2I6,F20.10)') 'J1,J2,A=',J1,J2,HESS(J1,J2)
C        ENDDO
C     ENDDO


      RETURN
      END

C*****************************************************************************

      SUBROUTINE LJPSHIFT_UPDATE_PAIRBOX(N, X, EPS, SIG6, RCONST, CONST, IRCUT2,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, VBOX,
     &           SBOX, FRACTIONAL )
        !Calculate the potential energy and box derivatives between atoms j1 and j2.
        !USE MODHESS


        IMPLICIT NONE

        INTEGER, INTENT(IN) :: J1, J2, N
        DOUBLE PRECISION, INTENT(IN) :: EPS, SIG6, RCONST, CONST, IRCUT2
        DOUBLE PRECISION, INTENT(IN) :: X(3*N)
        DOUBLE PRECISION, INTENT(IN) :: BOXLX, BOXLY, BOXLZ
        LOGICAL, INTENT(IN) :: FIXIMAGE, FRACTIONAL
        DOUBLE PRECISION, INTENT(INOUT) :: POTEL
        INTEGER, INTENT(INOUT) :: ANV(N,N,3)
        DOUBLE PRECISION, INTENT(INOUT) :: VBOX(3), SBOX(3)

        DOUBLE PRECISION :: R2, R6, VEC(3), G, R8, R14, FR2, TEMP
        INTEGER :: J3, J4, J5, J6

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
        R2=1.0D0/R2
        IF (R2.GT.IRCUT2) THEN
          !update potential energy
          R6=R2**3
          POTEL = POTEL + EPS*(SIG6*R6*(SIG6*R6-1.0D0) + RCONST/R2 + CONST)  ! AB

          !the derivative of potential energy
          R8=R6*R2
          R14=R8*R6
          G=-8.0D0*EPS*(3.0D0*(2.0D0*R14*(SIG6*SIG6)-R8*SIG6)-RCONST)

          !calculate the second derivative
          FR2 = 8.0D0*EPS * (84.0D0 * R14 * (SIG6*SIG6) - 24.0D0*R8 * SIG6) * R2

C  Box length derivatives
          IF (FRACTIONAL) THEN
            DO J5=1,3
              VBOX(J5) = VBOX(J5) + G*VEC(J5)**2
              SBOX(J5) = SBOX(J5) + (FR2 * VEC(J5)**2 + G) * VEC(J5)**2
            END DO
          ELSE
            DO J5=1,3
              !careful with signs: both VEC and ANV are antisymmetric in J1<->J2
              VBOX(J5) = VBOX(J5) - G * VEC(J5) * ANV(J2,J1,J5)
              SBOX(J5) = SBOX(J5) + (FR2 * VEC(J5)**2 + G) * ANV(J2,J1,J5)**2
            END DO
          ENDIF

        ENDIF

      END SUBROUTINE LJPSHIFT_UPDATE_PAIRBOX

C*****************************************************************************

      SUBROUTINE LJPSHIFT_UPDATE_PAIRS(N, X, EPS, SIG6, RCONST, CONST, IRCUT2,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, V)
        !Calculate the potential energy and gradient and hessian between atoms
        !j1 and j2.
        !Update also ANV
        USE MODHESS


        IMPLICIT NONE

        INTEGER, INTENT(IN) :: J1, J2, N
        DOUBLE PRECISION, INTENT(IN) :: EPS, SIG6, RCONST, CONST, IRCUT2
        DOUBLE PRECISION, INTENT(IN) :: X(3*N)
        DOUBLE PRECISION, INTENT(IN) :: BOXLX, BOXLY, BOXLZ
        LOGICAL, INTENT(IN) :: FIXIMAGE
        DOUBLE PRECISION, INTENT(INOUT) :: POTEL
        INTEGER, INTENT(INOUT) :: ANV(N,N,3)
        DOUBLE PRECISION, INTENT(INOUT) :: V(3*N)

        DOUBLE PRECISION :: R2, R6, VEC(3), G, R8, R14, FR2, TEMP
        INTEGER :: J3, J4, J5, J6

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
        R2=1.0D0/R2
        IF (R2.GT.IRCUT2) THEN
          !update potential energy
          R6=R2**3
          POTEL = POTEL + EPS*(SIG6*R6*(SIG6*R6-1.0D0) + RCONST/R2 + CONST)  ! AB

          !update derivative of potential energy
          R8=R6*R2
          R14=R8*R6
          G=-8.0D0*EPS*(3.0D0*(2.0D0*R14*(SIG6*SIG6)-R8*SIG6)-RCONST)
          DO J5=1,3
            !careful with signs: both VEC and V are antisymmetric in J1<->J2
            V(J3+J5)=V(J3+J5)+G*VEC(J5)
            V(J4+J5)=V(J4+J5)-G*VEC(J5)
          END DO


          !calculate the second derivative, and update the compoents of the
          !hessian
          FR2 = 8.0D0*EPS * (84.0D0 * R14 * (SIG6*SIG6) - 24.0D0*R8 * SIG6) * R2

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

      END SUBROUTINE LJPSHIFT_UPDATE_PAIRS

C*****************************************************************************

      SUBROUTINE LJPSHIFT_UPDATE_PAIRG(N, X, EPS, SIG6, RCONST, CONST, IRCUT2,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE, V)
        !Calculate the potential energy and gradient between atoms j1 and j2.
        !Update also ANV


        IMPLICIT NONE

        INTEGER, INTENT(IN) :: J1, J2, N
        DOUBLE PRECISION, INTENT(IN) :: EPS, SIG6, RCONST, CONST, IRCUT2
        DOUBLE PRECISION, INTENT(IN) :: X(3*N)
        DOUBLE PRECISION, INTENT(IN) :: BOXLX, BOXLY, BOXLZ
        LOGICAL, INTENT(IN) :: FIXIMAGE
        DOUBLE PRECISION, INTENT(INOUT) :: POTEL
        INTEGER, INTENT(INOUT) :: ANV(N,N,3)
        DOUBLE PRECISION, INTENT(INOUT) :: V(3*N)

        DOUBLE PRECISION :: R2, R6, VEC(3), G, R8, R14
        INTEGER :: J3, J4, J5

C       POTELDUM = 0.D0;

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
        R2=1.0D0/R2
        IF (R2.GT.IRCUT2) THEN
          R6=R2**3
          POTEL = POTEL + EPS*(SIG6*R6*(SIG6*R6-1.0D0) + RCONST/R2 + CONST)  ! AB

          !update derivative of potential energy
          R8=R6*R2
          R14=R8*R6
          G=-8.0D0*EPS*(3.0D0*(2.0D0*R14*(SIG6*SIG6)-R8*SIG6)-RCONST)
          DO J5=1,3
            V(J3+J5)=V(J3+J5)+G*VEC(J5)
            V(J4+J5)=V(J4+J5)-G*VEC(J5)
          END DO
        ENDIF

      END SUBROUTINE LJPSHIFT_UPDATE_PAIRG

C*****************************************************************************

      SUBROUTINE LJPSHIFT_UPDATE_PAIR(N, X, EPS, SIG6, RCONST, CONST, IRCUT2,
     &           POTEL, J1, J2, BOXLX, BOXLY, BOXLZ, ANV, FIXIMAGE )
        !calculate the potential energy between atoms j1 and j2 
        !update ANV


        IMPLICIT NONE

        INTEGER, INTENT(IN) :: J1, J2, N
        DOUBLE PRECISION, INTENT(IN) :: EPS, SIG6, RCONST, CONST, IRCUT2
        DOUBLE PRECISION, INTENT(IN) :: X(3*N)
        DOUBLE PRECISION, INTENT(IN) :: BOXLX, BOXLY, BOXLZ
        LOGICAL, INTENT(IN) :: FIXIMAGE
        DOUBLE PRECISION, INTENT(INOUT) :: POTEL
        INTEGER, INTENT(INOUT) :: ANV(N,N,3)

        DOUBLE PRECISION :: R2DUM, R6, VEC1, VEC2, VEC3
        INTEGER :: J3, J4

C       POTELDUM = 0.D0;

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
        R2DUM=VEC1**2+VEC2**2+VEC3**2
        R2DUM=1.0D0/R2DUM
        R6=R2DUM**3
        IF (R2DUM.GT.IRCUT2) THEN
          POTEL = POTEL + EPS*(SIG6*R6*(SIG6*R6-1.0D0) + RCONST/R2DUM + CONST)  ! AB
        ENDIF

      END SUBROUTINE LJPSHIFT_UPDATE_PAIR
