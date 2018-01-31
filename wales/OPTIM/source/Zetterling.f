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
C From fzet@nada.kth.se Thu Aug 24 08:56:18 2000
C Hello!
C I am sorry for the delay. 
C I have attached two potentials in the desired format in 
C this email, the first one is the gamma-brassy one, the 
C second one is the "do-not-crystallize" one. It is of 
C course the same potential, but with different values of 
C the parameters. The analytical expression is:
C 
C V=A*exp(alpha*r)*cos(2*K_F*r)/r**3 + B*(r/sigma)**pow + Const.
C 
C The constant is needed to have V=0 at the cutoff.
C (The cutoff is of course another parameter...)
C 
C Here is the first one:

!-----------------------------------------------------------------------*
C
C  Energy and gradient for the modified Dzugutov(Zetterling) potential-1.
C
C  23/08-2000
C  Fredrik Zetterling <fzet@pdc.kth.se>
C
! correspondence with the parameters in the file "potential.f"
! V     = GRAD
! EDZ   = EREAL
! GTEST = GRADT
!
! meaning of the variables
! DIST - distance
! EDZ  - total interaction energy
! VT   - (3*N) vector of interaction energies of each particle
! G    - (3*N,3*N) matrix of "forces" (first derivatives of the potential)/r
! V    - (3*N) vector - gradient (force acting on a particle)
! 
!-----------------------------------------------------------------------*
      SUBROUTINE Z1(NATOMS,X,V,EDZ,GTEST,STEST,BOXLX,BOXLY,BOXLZ)
      use porfuncs
      USE MODHESS
      IMPLICIT NONE
      INTEGER J1, J2, J3, J4, NATOMS
      LOGICAL GTEST, STEST
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), G(NATOMS,NATOMS), 
     1                 EDZ, DUMMYX, DUMMYY, DUMMYZ, XMUL2, DUMMY,
     2                 DDUMMY, DDDUMMY, RCUT2, BOXLX,BOXLY,BOXLZ
      PARAMETER (RCUT2=7.0176544700968D0)

      IF (STEST) THEN
         PRINT*,'second derivatives not implemented for Z1'
         STOP
      ENDIF
      EDZ=0.0D0
      IF (GTEST) THEN
         DO J1=1,NATOMS
            J3=3*J1
            G(J1,J1)=0.0D0
            DO J2=J1+1,NATOMS
               J4=3*J2
                     DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               G(J2,J1)=0.0D0
               IF (DIST.LT.RCUT2) THEN 
                  call derphiz1(DIST,DUMMY,DDUMMY,DDDUMMY,.true.,.false.)
                  EDZ=EDZ+DUMMY
                  G(J2,J1)=DDUMMY
               ENDIF
               G(J1,J2)=G(J2,J1)
            ENDDO
         ENDDO
      ELSE
         DO J1=1,NATOMS
            J3=3*J1
            DO J2=J1+1,NATOMS
               J4=3*J2
               DIST=(X(J3-2)-X(J4-2))**2+(X(J3-1)-X(J4-1))**2+(X(J3)-X(J4))**2
               IF (DIST.LT.RCUT2) THEN 
                  call derphiz1(DIST,DUMMY,DDUMMY,DDDUMMY,.false.,.false.)
                  EDZ=EDZ+DUMMY
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      IF (.NOT.GTEST) RETURN

      DO J1=1,NATOMS
         J3=3*J1
         DUMMYX=0.0D0
         DUMMYY=0.0D0
         DUMMYZ=0.0D0
         DO J4=1,NATOMS
            J2=3*J4
            XMUL2=G(J4,J1)
            DUMMYX=DUMMYX+XMUL2*(X(J3-2)-X(J2-2))
            DUMMYY=DUMMYY+XMUL2*(X(J3-1)-X(J2-1))
            DUMMYZ=DUMMYZ+XMUL2*(X(J3)  -X(J2))
         ENDDO
         V(J3-2)=DUMMYX
         V(J3-1)=DUMMYY
         V(J3)=DUMMYZ
      ENDDO

      RETURN
      END

!%%%%%%%
! values of the potential, first, and second derivatives
! if secder = .false. ddphi = 0 on output
!%%%%%%%
C  ***  The derivative is divided by R !
!-----------------------------------------------------------------------*
      subroutine derphiz1(r2,phi,dphi,ddphi,gtest,secder)
      implicit none
      DOUBLE PRECISION R, R2, PHI, DPHI, DDPHI
      logical secder, gtest

      DOUBLE PRECISION A,B,KF,ALPH,SIG,POW,ST,RCUT ! PAIR POTENTIAL PARAMETERS
      DOUBLE PRECISION DUMMY1,DUMMY2,DUMMY3

      parameter ( A=1.58 )
      parameter ( B=420000000.0 )
      parameter ( KF=4.12 )
      parameter ( ALPH=-0.22 )
      parameter ( SIG=0.331 )
      parameter ( POW=-18 )
C     parameter ( ST=0.046822 )
C     parameter ( RCUT=2.65d0 )
      parameter ( ST=0.04682632412414856D0 )
      parameter ( RCUT=2.649085591311991D0 )

      R = SQRT(R2)

      DUMMY1=B*(r/SIG)**POW
      DUMMY2=A*DEXP(ALPH*r)/r**3
      DUMMY3=DCOS(2*KF*r)
      phi=ST+DUMMY2*DUMMY3+DUMMY1

      IF (GTEST) dphi=(DUMMY2*(DUMMY3*(ALPH-3/r)-2*KF*DSIN(2*KF*r))+POW*DUMMY1/r)/r
            
      IF (SECDER) THEN 
         ddphi=DBLE(POW)*DBLE(POW-1)*B*((r/SIG)**POW)/(r**2)
     1      +A*(DEXP(ALPH*r)/(r**3))*((ALPH-3/r)*((ALPH-3/r)
     2      *DCOS(2*KF*r)-4*KF*DSIN(2*KF*r))+(3/(r**2)-4*KF*KF)
     3      *DCOS(2*KF*r))
      ENDIF

      end                       ! subroutine derphi
!-----------------------------------------------------------------------*

C And here is the second one:

!-----------------------------------------------------------------------*
C
C  Energy and gradient for the modified Dzugutov(Zetterling) potential-2.
C
C  23/08-2000
C  Fredrik Zetterling <fzet@pdc.kth.se>
C
! correspondence with the parameters in the file "potential.f"
! V     = GRAD
! EDZ   = EREAL
! GTEST = GRADT
!
! meaning of the variables
! DIST - distance
! EDZ  - total interaction energy
! VT   - (3*N) vector of interaction energies of each particle
! G    - (3*N,3*N) matrix of "forces" (first derivatives of the potential)/r
! V    - (3*N) vector - gradient (force acting on a particle)
! 
!-----------------------------------------------------------------------*
      SUBROUTINE Z2(NATOMS,X,V,EDZ,GTEST,STEST,BOXLX,BOXLY,BOXLZ)
      use porfuncs
      USE KEY
      IMPLICIT NONE
      LOGICAL GTEST, STEST
      INTEGER J1, J2, J3, J4, NATOMS, ANV(3)
      DOUBLE PRECISION X(3*NATOMS), DIST, V(3*NATOMS), ! G(NATOMS,NATOMS), 
     1                 EDZ, DUMMYX, DUMMYY, DUMMYZ, XMUL2, DUMMY, VEC(3),
     2                 DDUMMY, DDDUMMY, RCUT2, BOXLX,BOXLY,BOXLZ
      PARAMETER (RCUT2=6.995378792429925D0)

      IF (STEST) THEN
         PRINT*,'second derivatives not implemented for Z2'
         STOP
      ENDIF
C
C  Deal with any atoms that have left the box.
C
C     PRINT*,'In Zetterling FIXIMAGE=',FIXIMAGE
      IF (BULKT.AND.(.NOT.FIXIMAGE).AND.(.NOT.NORESET)) THEN
         DO J1=1,NATOMS
            J2 = 3*(J1-1)
            X(J2+1)=X(J2+1) - BOXLX*ANINT(X(J2+1)/BOXLX)
            X(J2+2)=X(J2+2) - BOXLY*ANINT(X(J2+2)/BOXLY)
            X(J2+3)=X(J2+3) - BOXLZ*ANINT(X(J2+3)/BOXLZ)
         ENDDO
      ENDIF

      EDZ=0.0D0
      IF (GTEST) THEN
         V(1:3*NATOMS)=0.0D0
         DO J1=1,NATOMS
            J3=3*(J1-1)
            DO J2=J1+1,NATOMS
               J4=3*(J2-1)
               VEC(1)=X(J3+1)-X(J4+1)
               VEC(2)=X(J3+2)-X(J4+2)
               VEC(3)=X(J3+3)-X(J4+3)
               ANV(1:3)=0
               IF (BULKT.AND.(.NOT.FIXIMAGE)) THEN
                  ANV(1)=NINT(VEC(1)/BOXLX)
                  ANV(2)=NINT(VEC(2)/BOXLY)
                  ANV(3)=NINT(VEC(3)/BOXLZ)
               ENDIF

               VEC(1)=VEC(1)-BOXLX*ANV(1)
               VEC(2)=VEC(2)-BOXLY*ANV(2)
               VEC(3)=VEC(3)-BOXLZ*ANV(3)

               DIST=VEC(1)**2+VEC(2)**2+VEC(3)**2

               IF (DIST.LT.RCUT2) THEN
                  CALL DERPHIZ2(DIST,DUMMY,DDUMMY,DDDUMMY,.TRUE.,.FALSE.)
                  EDZ=EDZ+DUMMY

                  V(J3+1)=V(J3+1)+DDUMMY*VEC(1) 
                  V(J3+2)=V(J3+2)+DDUMMY*VEC(2)
                  V(J3+3)=V(J3+3)+DDUMMY*VEC(3)

                  V(J4+1)=V(J4+1)-DDUMMY*VEC(1) 
                  V(J4+2)=V(J4+2)-DDUMMY*VEC(2)
                  V(J4+3)=V(J4+3)-DDUMMY*VEC(3)

               ENDIF
            ENDDO
         ENDDO
      ELSE
         DO J1=1,NATOMS
            J3=3*(J1-1)
            DO J2=J1+1,NATOMS
               J4=3*(J2-1)
               VEC(1)=X(J3+1)-X(J4+1)
               VEC(2)=X(J3+2)-X(J4+2)
               VEC(3)=X(J3+3)-X(J4+3)
               ANV(1:3)=0
               IF (BULKT.AND.(.NOT.FIXIMAGE)) THEN
                  ANV(1)=NINT(VEC(1)/BOXLX)
                  ANV(2)=NINT(VEC(2)/BOXLY)
                  ANV(3)=NINT(VEC(3)/BOXLZ)
               ENDIF
               VEC(1)=VEC(1)-BOXLX*ANV(1)
               VEC(2)=VEC(2)-BOXLY*ANV(2)
               VEC(3)=VEC(3)-BOXLZ*ANV(3)

               DIST=VEC(1)**2+VEC(2)**2+VEC(3)**2

               IF (DIST.LT.RCUT2) THEN
                  CALL DERPHIZ2(DIST,DUMMY,DDUMMY,DDDUMMY,.FALSE.,.FALSE.)
                  EDZ=EDZ+DUMMY
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      RETURN
      END

!%%%%%%%
! values of the potential, first, and second derivatives
! if secder = .false. ddphi = 0 on output
!%%%%%%%
!-----------------------------------------------------------------------*
      subroutine derphiz2(r2,phi,dphi,ddphi,gtest,secder)
      implicit none
      DOUBLE PRECISION R, R2, PHI, DPHI, DDPHI
      logical secder, gtest

      DOUBLE PRECISION A,B,KF,ALPH,SIG,POW,ST,RCUT ! PAIR POTENTIAL PARAMETERS
      DOUBLE PRECISION DUMMY1,DUMMY2,DUMMY3

      parameter ( A=1.04D0 )
      parameter ( B=4200000.0D0 )
      parameter ( KF=4.139D0 )
      parameter ( ALPH=0.33D0 )
      parameter ( SIG=0.348D0 )
      parameter ( POW=-14.5D0 )
C     parameter ( ST=0.133915D0 )
C     parameter ( RCUT=2.645D0 )
      parameter ( ST=0.1339154253770228D0 )
      parameter ( RCUT=2.644877840738571D0 )

      r = sqrt(r2)

      DUMMY1=B*(r/SIG)**POW
      DUMMY2=A*DEXP(ALPH*r)/r**3
      DUMMY3=DCOS(2*KF*r)
      phi=ST+DUMMY2*DUMMY3+DUMMY1

      IF (GTEST) dphi=(DUMMY2*(DUMMY3*(ALPH-3/r)-2*KF*DSIN(2*KF*r))+POW*DUMMY1/r)/R

      if (secder) then 
            
         ddphi=DBLE(POW)*DBLE(POW-1)*B*((r/SIG)**POW)/(r**2)
     1      +A*(DEXP(ALPH*r)/(r**3))*((ALPH-3/r)*((ALPH-3/r)
     2      *DCOS(2*KF*r)-4*KF*DSIN(2*KF*r))+(3/(r**2)-4*KF*KF)
     3      *DCOS(2*KF*r))

      endif
  
      end                       ! subroutine derphi
!-----------------------------------------------------------------------*
