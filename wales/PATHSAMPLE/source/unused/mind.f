!  License info {{{
!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999-2009 David J. Wales
!   This file is part of PATHSAMPLE.
!
!   PATHSAMPLE is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   PATHSAMPLE is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
! }}}

C>     finds the minimum distance between two geometries
C>     minimization by LBFGS
C
      SUBROUTINE MIND(RA,RB,NATOMS,DIST,BULKT,TWOD,ZSYM)
! Declarations {{{
      USE PORFUNCS
      IMPLICIT NONE
      INTEGER J1, IG, I, ITER, NCOUNT, NATOMS
      DOUBLE PRECISION r(3,NATOMS),r0(3,NATOMS),R1(3,NATOMS),RA(3*NATOMS),RB(3*NATOMS)
      DOUBLE PRECISION P(3),DIST, DIST0, RG, RG0, OMEGA(3,3), RANDOM, R1SAVE(3,NATOMS), DSAVE,
     &                 OVEC(3), H1VEC(3), H2VEC(3), DPRAND
      INTEGER NSIZE
      LOGICAL BULKT, MFLAG, AGAIN, TWOD
C     common /geom /r,r0,r1,nsize
      COMMON /MINDOM/ OMEGA
      CHARACTER*2 ZSYM
! }}}
! Subroutine body {{{

      AGAIN=.FALSE.
      NCOUNT=0
      IF (ZSYM(1:1).EQ.'W') THEN
        ! RA => (via CONVERT) OVEC,H1VEC,H2VEC => R 
        ! RB => (via CONVERT) OVEC,H1VEC,H2VEC => R0 
        !{{{
         NSIZE=3*NATOMS
         DO J1=1,NATOMS
            CALL CONVERT(RA(3*(J1-1)+1),RA(3*(J1-1)+2),RA(3*(J1-1)+3),
     1                 RA(3*(NATOMS+J1-1)+1),RA(3*(NATOMS+J1-1)+2),RA(3*(NATOMS+J1-1)+3),OVEC,H1VEC,H2VEC)
            R(1,(J1-1)*3+1)=OVEC(1)
            R(2,(J1-1)*3+1)=OVEC(2)
            R(3,(J1-1)*3+1)=OVEC(3)
            R(1,(J1-1)*3+2)=H1VEC(1)
            R(2,(J1-1)*3+2)=H1VEC(2)
            R(3,(J1-1)*3+2)=H1VEC(3)
            R(1,(J1-1)*3+3)=H2VEC(1)
            R(2,(J1-1)*3+3)=H2VEC(2)
            R(3,(J1-1)*3+3)=H2VEC(3)
            CALL CONVERT(RB(3*(J1-1)+1),RB(3*(J1-1)+2),RB(3*(J1-1)+3),
     1                 RB(3*(NATOMS+J1-1)+1),RB(3*(NATOMS+J1-1)+2),RB(3*(NATOMS+J1-1)+3),OVEC,H1VEC,H2VEC)
            R0(1,(J1-1)*3+1)=OVEC(1)
            R0(2,(J1-1)*3+1)=OVEC(2)
            R0(3,(J1-1)*3+1)=OVEC(3)
            R0(1,(J1-1)*3+2)=H1VEC(1)
            R0(2,(J1-1)*3+2)=H1VEC(2)
            R0(3,(J1-1)*3+2)=H1VEC(3)
            R0(1,(J1-1)*3+3)=H2VEC(1)
            R0(2,(J1-1)*3+3)=H2VEC(2)
            R0(3,(J1-1)*3+3)=H2VEC(3)
         ENDDO
         ! }}}
      ELSE
        ! RA => R, RB => R0
        ! {{{
         nsize=NATOMS
         DO J1=1,NSIZE
            R(1,J1)=RA(3*(J1-1)+1)
            R(2,J1)=RA(3*(J1-1)+2)
            R(3,J1)=RA(3*(J1-1)+3)
            R0(1,J1)=RB(3*(J1-1)+1)
            R0(2,J1)=RB(3*(J1-1)+2)
            R0(3,J1)=RB(3*(J1-1)+3)
         ENDDO
         ! }}}
      ENDIF
C
C     put c.o.m. to origin {{{
C
      do ig=1,3
         rg=0.d0
         rg0=0.d0
         do i=1,nsize
            rg=rg+r(ig,i)
            rg0=rg0+r0(ig,i)
         enddo
         rg=rg/nsize
         rg0=rg0/nsize
         do i=1,nsize
            r(ig,i)=r(ig,i)-rg
            r0(ig,i)=r0(ig,i)-rg0
         enddo
      enddo

! }}}
C
C     initial angles
C
11    CONTINUE
      P(1)=0.0D0
      P(2)=0.0D0
      P(3)=0.0D0
C
C     calculate initial distance
C
      NCOUNT=0
10    CALL DISTFUNC(P,R,R0,R1,NSIZE,DIST0,NATOMS)
      DIST0=SQRT(DIST0)
      IF (BULKT) THEN
         DIST=DIST0
         RETURN
      ENDIF

      CALL MMYLBFGS(P,1.0D-6,MFLAG,DIST,200,ITER,TWOD,R,R0,R1,NSIZE,NATOMS)
      CALL DISTFUNC(P,R,R0,R1,NSIZE,DIST,NATOMS)
      DIST=SQRT(DIST)
      IF (MFLAG) THEN
C        WRITE(*,'(A,2F15.5,A,I6)') 'Initial and final distances:',DIST0,DIST,' iterations=',ITER
C        PRINT*
      ELSE
         NCOUNT=NCOUNT+1
         IF (NCOUNT.GT.0) THEN
            PRINT*,'WARNING - convergence failure in mind'
C           STOP
         ELSE
C           WRITE(*,'(A,2F15.5,A,I6,A,I6)') 'Initial and final distances:',DIST0,DIST,' iterations=',ITER,' NCOUNT=',NCOUNT
            DO J1=1,NSIZE
               R0(1,J1)=R1(1,J1)
               R0(2,J1)=R1(2,J1)
               R0(3,J1)=R1(3,J1)
            ENDDO
            IF (.NOT.TWOD) THEN
C              CALL RANDOM_NUMBER(RANDOM)
               RANDOM=DPRAND()
               P(1)=2*(RANDOM-0.5D0)/10.0D0
C              CALL RANDOM_NUMBER(RANDOM)
               RANDOM=DPRAND()
               P(2)=2*(RANDOM-0.5D0)/10.0D0
            ENDIF
C           CALL RANDOM_NUMBER(RANDOM)
            RANDOM=DPRAND()
            P(3)=2*(RANDOM-0.5D0)/10.0D0
            GOTO 10
         ENDIF
      ENDIF

C
C  This block allows the second geometry to rotate out of the xy plane. This is allowed
C  in mind for pathsample, where we need to align pathways using ORIENT. It is not 
C  allowed in OPTIM when we are checking for the minimum distance between permutational
C  isomers.
C
      IF (TWOD) THEN
         IF (AGAIN) THEN
C           WRITE(*,'(A,2F20.10)') 'DIST,DSAVE=',DIST,DSAVE
            IF (DIST.GT.DSAVE) THEN
               DIST=DSAVE
               DO J1=1,NSIZE
                  R1(1,J1)=R1SAVE(1,J1)
                  R1(2,J1)=R1SAVE(2,J1)
                  R1(3,J1)=R1SAVE(3,J1)
                  R0(1,J1)=-R0(1,J1)
               ENDDO
            ENDIF
         ELSE
            AGAIN=.TRUE.
            DSAVE=DIST
            DO J1=1,NSIZE
               R0(1,J1)=-R0(1,J1)
               R1SAVE(1,J1)=R1(1,J1)
               R1SAVE(2,J1)=R1(2,J1)
               R1SAVE(3,J1)=R1(3,J1)
            ENDDO
            GOTO 10
         ENDIF
      ENDIF


      IF (ZSYM(1:1).EQ.'W') THEN
         DO J1=1,NATOMS
            OVEC(1)=R(1,(J1-1)*3+1)
            OVEC(2)=R(2,(J1-1)*3+1)
            OVEC(3)=R(3,(J1-1)*3+1)
            H1VEC(1)=R(1,(J1-1)*3+2)
            H1VEC(2)=R(2,(J1-1)*3+2)
            H1VEC(3)=R(3,(J1-1)*3+2)
            H2VEC(1)=R(1,(J1-1)*3+3)
            H2VEC(2)=R(2,(J1-1)*3+3)
            H2VEC(3)=R(3,(J1-1)*3+3)
            CALL CONVERT2(OVEC,H1VEC,H2VEC,RA(3*(J1-1)+1),RA(3*(J1-1)+2),RA(3*(J1-1)+3),
     1                 RA(3*(NATOMS+J1-1)+1),RA(3*(NATOMS+J1-1)+2),RA(3*(NATOMS+J1-1)+3))
            OVEC(1)=R1(1,(J1-1)*3+1)
            OVEC(2)=R1(2,(J1-1)*3+1)
            OVEC(3)=R1(3,(J1-1)*3+1)
            H1VEC(1)=R1(1,(J1-1)*3+2)
            H1VEC(2)=R1(2,(J1-1)*3+2)
            H1VEC(3)=R1(3,(J1-1)*3+2)
            H2VEC(1)=R1(1,(J1-1)*3+3)
            H2VEC(2)=R1(2,(J1-1)*3+3)
            H2VEC(3)=R1(3,(J1-1)*3+3)
            CALL CONVERT2(OVEC,H1VEC,H2VEC,RB(3*(J1-1)+1),RB(3*(J1-1)+2),RB(3*(J1-1)+3),
     1                 RB(3*(NATOMS+J1-1)+1),RB(3*(NATOMS+J1-1)+2),RB(3*(NATOMS+J1-1)+3))
         ENDDO
      ELSE
         DO J1=1,NSIZE
            RB(3*(J1-1)+1)=R1(1,J1)
            RB(3*(J1-1)+2)=R1(2,J1)
            RB(3*(J1-1)+3)=R1(3,J1)
            RA(3*(J1-1)+1)=R(1,J1)
            RA(3*(J1-1)+2)=R(2,J1)
            RA(3*(J1-1)+3)=R(3,J1)
         ENDDO
      ENDIF
! }}}
      RETURN
      END

c_____________________________________________________________________

      subroutine prod(A,B,C)
! Declarations & sub. body {{{
      implicit DOUBLE PRECISION(a-h,o-z)
      dimension A(3,3),B(3,3),C(3,3)

      do i=1,3
         do j=1,3
            DUMMY=0.d0
            do k=1,3
               DUMMY=DUMMY+A(i,k)*B(k,j)
            enddo
            C(i,j)=DUMMY
         enddo
      enddo
! }}}
      return
      end
c_____________________________________________________________________

      SUBROUTINE DISTFUNC(P,R,R0,R1,NSIZE,DIST,NATOMS)
! Declarations {{{
      IMPLICIT NONE

      INTEGER NSIZE, NATOMS, I, J
      DOUBLE PRECISION R(3,NATOMS),r0(3,NATOMS),A(3,3),R1(3,NATOMS),DIST

      DOUBLE PRECISION omega(3,3),rotx(3,3),roty(3,3),rotz(3,3)
      DOUBLE PRECISION P(3)
      COMMON /MINDOM/ OMEGA
! }}}
! Subroutine body {{{

C     common/geom/r,r0,r1,nsize

c     calculate rotation matrices

      rotx(1,1)=1.0D0
      rotx(1,2)=0.0D0
      rotx(2,1)=0.0D0
      rotx(1,3)=0.0D0
      rotx(3,1)=0.0D0
      rotx(2,2)=dcos(p(1))
      rotx(3,3)=rotx(2,2)
      rotx(2,3)=dsin(p(1))
      rotx(3,2)=-rotx(2,3)

      roty(1,2)=0.0D0
      roty(2,1)=0.0D0
      roty(2,3)=0.0D0
      roty(3,2)=0.0D0
      roty(2,2)=1.0D0
      roty(1,1)=dcos(p(2))
      roty(3,3)=roty(1,1)
      roty(1,3)=-dsin(p(2))
      roty(3,1)=-roty(1,3)

      rotz(1,3)=0.0D0
      rotz(3,1)=0.0D0
      rotz(2,3)=0.0D0
      rotz(3,2)=0.0D0
      rotz(3,3)=1.0D0
      rotz(1,1)=dcos(p(3))
      rotz(2,2)=rotz(1,1)
      rotz(1,2)=dsin(p(3))
      rotz(2,1)=-rotz(1,2)

      A(1,1)=ROTY(1,1)*ROTZ(1,1)+ROTY(1,2)*ROTZ(2,1)+ROTY(1,3)*ROTZ(3,1)
      A(1,2)=ROTY(1,1)*ROTZ(1,2)+ROTY(1,2)*ROTZ(2,2)+ROTY(1,3)*ROTZ(3,2)
      A(1,3)=ROTY(1,1)*ROTZ(1,3)+ROTY(1,2)*ROTZ(2,3)+ROTY(1,3)*ROTZ(3,3)
      A(2,1)=ROTY(2,1)*ROTZ(1,1)+ROTY(2,2)*ROTZ(2,1)+ROTY(2,3)*ROTZ(3,1)
      A(2,2)=ROTY(2,1)*ROTZ(1,2)+ROTY(2,2)*ROTZ(2,2)+ROTY(2,3)*ROTZ(3,2)
      A(2,3)=ROTY(2,1)*ROTZ(1,3)+ROTY(2,2)*ROTZ(2,3)+ROTY(2,3)*ROTZ(3,3)
      A(3,1)=ROTY(3,1)*ROTZ(1,1)+ROTY(3,2)*ROTZ(2,1)+ROTY(3,3)*ROTZ(3,1)
      A(3,2)=ROTY(3,1)*ROTZ(1,2)+ROTY(3,2)*ROTZ(2,2)+ROTY(3,3)*ROTZ(3,2)
      A(3,3)=ROTY(3,1)*ROTZ(1,3)+ROTY(3,2)*ROTZ(2,3)+ROTY(3,3)*ROTZ(3,3)

      OMEGA(1,1)=ROTX(1,1)*A(1,1)+ROTX(1,2)*A(2,1)+ROTX(1,3)*A(3,1)
      OMEGA(1,2)=ROTX(1,1)*A(1,2)+ROTX(1,2)*A(2,2)+ROTX(1,3)*A(3,2)
      OMEGA(1,3)=ROTX(1,1)*A(1,3)+ROTX(1,2)*A(2,3)+ROTX(1,3)*A(3,3)
      OMEGA(2,1)=ROTX(2,1)*A(1,1)+ROTX(2,2)*A(2,1)+ROTX(2,3)*A(3,1)
      OMEGA(2,2)=ROTX(2,1)*A(1,2)+ROTX(2,2)*A(2,2)+ROTX(2,3)*A(3,2)
      OMEGA(2,3)=ROTX(2,1)*A(1,3)+ROTX(2,2)*A(2,3)+ROTX(2,3)*A(3,3)
      OMEGA(3,1)=ROTX(3,1)*A(1,1)+ROTX(3,2)*A(2,1)+ROTX(3,3)*A(3,1)
      OMEGA(3,2)=ROTX(3,1)*A(1,2)+ROTX(3,2)*A(2,2)+ROTX(3,3)*A(3,2)
      OMEGA(3,3)=ROTX(3,1)*A(1,3)+ROTX(3,2)*A(2,3)+ROTX(3,3)*A(3,3)

C     call prod(roty,rotz,A)
C     call prod(rotx,A,omega)

      do i=1,nsize
         do j=1,3
            R1(j,i)=omega(j,1)*R0(1,i)+omega(j,2)*R0(2,i)+omega(j,3)*R0(3,i)
         enddo
      enddo

      dist=0.d0
      do i=1,nsize
         dist=dist+(r(1,i)-r1(1,i))**2+(r(2,i)-r1(2,i))**2+(r(3,i)-r1(3,i))**2
      enddo
! }}}
      return
      end
c_____________________________________________________________________

!> \param[in] P,R,R0,R1,NSIZE,NATOMS
!> \param[out] F

      SUBROUTINE DFUNC(P,F,R,R0,R1,NSIZE,NATOMS)
! {{{
      IMPLICIT NONE

      INTEGER NSIZE, NATOMS, I
      DOUBLE PRECISION r(3,NATOMS),r0(3,NATOMS),r1(3,NATOMS)
      DOUBLE PRECISION P(3),F(3)

      F(1)=0.d0
      F(2)=0.d0
      F(3)=0.d0
      DO I=1,NSIZE
         F(1)=F(1)+R(2,i)*R1(3,i)-R(3,i)*R1(2,i)
         F(2)=F(2)+R(3,i)*R1(1,i)-R(1,i)*R1(3,i)
         F(3)=F(3)+R(1,i)*R1(2,i)-R(2,i)*R1(1,i)
      ENDDO
      F(1)=-2.d0*F(1)
      F(2)=-2.d0*F(2)
      F(3)=-2.d0*F(3)
! }}}
      RETURN
      END
c_____________________________________________________________________

      SUBROUTINE DDFUNC(P,H,R,R0,R1,NSIZE,NATOMS)
! {{{
      IMPLICIT NONE

      INTEGER NSIZE, I, J, NATOMS
      DOUBLE PRECISION r(3,NATOMS),r0(3,NATOMS),r1(3,NATOMS)
      DOUBLE PRECISION P(3),H(3,3)

      do i=1,3
         do j=1,3
            H(i,j)=0.d0
            H(i,j)=0.d0
            H(i,j)=0.d0
         enddo
      enddo
      do i=1,nsize
         H(1,1)=H(1,1)+R(2,i)*R1(2,i)+R(3,i)*R1(3,i)
         H(2,2)=H(2,2)+R(1,i)*R1(1,i)+R(3,i)*R1(3,i)
         H(3,3)=H(3,3)+R(1,i)*R1(1,i)+R(2,i)*R1(2,i)
         H(2,1)=H(2,1)-R(2,i)*R1(1,i)
         H(3,1)=H(3,1)-R(3,i)*R1(1,i)
         H(3,2)=H(3,2)-R(3,i)*R1(2,i)
      enddo
      H(1,2)=H(2,1)
      H(1,3)=H(3,1)
      H(2,3)=H(3,2)

      do i=1,3
         do j=1,3
            H(i,j)=2.d0*H(i,j)
         enddo
      enddo
! }}}
      return

      end

! Doxygen {{{
!>
!> \brief Given cartesian coordinates (X,Y,Z) and Euler angles (A,B,C), convert to HHO positions  \n
!>
!> \param[in] X,Y,Z Cartesian coordinates
!> \param[in] A,B,C Euler angles 
!> \param[out] H1VEC, H2VEC Hydrogen positions 
!> \param[out] OVEC Oxygen positions 
!>  
! }}}
      SUBROUTINE CONVERT(X,Y,Z,A,B,C,OVEC,H1VEC,H2VEC)
      ! Declarations {{{
      IMPLICIT NONE

      DOUBLE PRECISION X,Y,Z,A,B,C

      DOUBLE PRECISION OVEC(3),H1VEC(3),H2VEC(3)

      DOUBLE PRECISION PI,OL,TH,PH,PS,SINA,SINB,SINC,COSA,COSB,COSC,SP12,SP13,SP22,
     1                 SP23,SP32,SP33,SY12,SY22,SY32,SZ13,SZ23,SZ33,ANGLE,HOLEN,RANGLE,OHZ,HYL,HZL
      DATA ANGLE,HOLEN/104.52D0,0.9572D0/

      PI=4.0D0*ATAN(1.0D0)
      RANGLE=PI*ANGLE/360.0D0
      OHZ=HOLEN*COS(RANGLE)
      HYL=HOLEN*SIN(RANGLE)
      HZL=16.0D0*OHZ/18.0D0
      OL=-OHZ+HZL
C}}}
! Subroutine body {{{
CCCCCC    A,B,C ARE EULER ANGLES {{{
C
      TH=A
      PH=B
      PS=C
      SINA=SIN(TH)
      SINB=SIN(PH)
      SINC=SIN(PS)
      COSA=COS(TH)
      COSB=COS(PH)
      COSC=COS(PS)
      SP12=-(SINC*COSB+COSA*SINB*COSC)
      SP13=SINA*SINB
      SP22=-SINC*SINB+COSA*COSB*COSC
      SP23=-SINA*COSB
      SP32=SINA*COSC
      SP33=COSA
      SY12=SP12*HYL
      SY22=SP22*HYL
      SY32=SP32*HYL
      SZ13=SP13*HZL
      SZ23=SP23*HZL
      SZ33=SP33*HZL
C}}}
CCCCC HYDROGEN POSITIONS {{{
C
      H1VEC(1)=SY12+SZ13+X
      H1VEC(2)=SY22+SZ23+Y
      H1VEC(3)=SY32+SZ33+Z

      H2VEC(1)=-SY12+SZ13+X
      H2VEC(2)=-SY22+SZ23+Y
      H2VEC(3)=-SY32+SZ33+Z
C}}}
CCCC OXYGEN POSITION    {{{
C
      OVEC(1)=SP13*OL   +X
      OVEC(2)=SP23*OL   +Y
      OVEC(3)=SP33*OL   +Z
C }}}
! }}}
      RETURN
      END
! Doxygen {{{
!
!> \brief convert HHO positions to Euler angles (inverse to CONVERT)
!>
!>  Convert H, H, O positions to Eulers - needs H, H, O
!>  positions and centres of mass. Here we allow for the possibility that the
!>  ideal rigid water geometry may be broken and take the best fit.
!>
!> \param[out] A,B,C
!>
!> \param[in] X,Y,Z  center of mass coordinates 
!> \param[in] OVEC   oxygen positions 
!> \param[in] H1VEC  hydrogen positions
!> \param[in] H2VEC      --//--
!>
! }}}
       SUBROUTINE CONVERT2(OVEC,H1VEC,H2VEC,X,Y,Z,A,B,C)
! Declarations {{{
     IMPLICIT NONE
! output
      DOUBLE PRECISION A,B,C
! input
      DOUBLE PRECISION X,Y,Z,OVEC(3),H1VEC(3),H2VEC(3)
! other parameters  
      DOUBLE PRECISION PI,PI2, OL

      DOUBLE PRECISION SINA,SINB,SINC,COSA,COSB,COSC

      DOUBLE PRECISION SP12,SP13,SP22,SP23,SP32,SP33

      DOUBLE PRECISION SY12,SY22,SY32,SZ13,SZ23,SZ33

      DOUBLE PRECISION ANGLE,HOLEN,RANGLE,OHZ,HYL,HZL

      DOUBLE PRECISION SP13T,SP23T,SP33T

      DOUBLE PRECISION ASAVE,BSAVE,CSAVE

      DOUBLE PRECISION SY12T,SY22T,SY32T,SZ13T,SZ23T,SZ33T

      DOUBLE PRECISION DMIN

      DOUBLE PRECISION ABEST,BBEST,CBEST
      
      INTEGER JA,JB,JC
      DATA ANGLE,HOLEN/104.52D0,0.9572D0/
! }}}
! subroutine body  {{{
      PI=4.0D0*ATAN(1.0D0)
      PI2=2.0D0*PI
      RANGLE=PI*ANGLE/360.0D0
      OHZ=HOLEN*COS(RANGLE)
      HYL=HOLEN*SIN(RANGLE)
      HZL=16.0D0*OHZ/18.0D0
      OL=-OHZ+HZL
C
      X=(H1VEC(1)+H2VEC(1)+16.0D0*OVEC(1))/18.0D0
      Y=(H1VEC(2)+H2VEC(2)+16.0D0*OVEC(2))/18.0D0
      Z=(H1VEC(3)+H2VEC(3)+16.0D0*OVEC(3))/18.0D0


      COSA=((OVEC(3)-Z)/OL)
      IF (ABS(COSA).GT.1.0D0) THEN
         COSA=ABS(COSA)/COSA
         A=DACOS(-COSA)
      ELSE
         A=DACOS((OVEC(3)-Z)/OL)
      ENDIF
      SINA=DSIN(A)
      IF (SINA.EQ.0.0D0) SINA=1.0D-10

      COSB=-(OVEC(2)-Y)/(OL*SINA)
      IF (ABS(COSB).GT.1.0D0) THEN
         COSB=ABS(COSB)/COSB
         B=DACOS(-COSB)
      ELSE
         B=DACOS(-(OVEC(2)-Y)/(OL*SINA))
      ENDIF
      SINB=DSIN(B)

      COSC=(H1VEC(3)-H2VEC(3))/(2.0D0*HYL*SINA)
      IF (ABS(COSC).GT.1.0D0) THEN
         COSC=ABS(COSC)/COSC
         C=DACOS(-COSC)
      ELSE
         C=DACOS((H1VEC(3)-H2VEC(3))/(2.0D0*HYL*SINA))
      ENDIF
      SINC=DSIN(C)

      SP13T=(OVEC(1)-X)/OL
      SP23T=(OVEC(2)-Y)/OL
      SP33T=(OVEC(3)-Z)/OL

      SY12T=(H1VEC(1)-H2VEC(1))/(2.0D0)
      SY22T=(H1VEC(2)-H2VEC(2))/(2.0D0)
      SY32T=(H1VEC(3)-H2VEC(3))/(2.0D0)

      SZ13T=(H1VEC(1)+H2VEC(1)-2.0D0*X)/(2.0D0)
      SZ23T=(H1VEC(2)+H2VEC(2)-2.0D0*Y)/(2.0D0)
      SZ33T=(H1VEC(3)+H2VEC(3)-2.0D0*Z)/(2.0D0)

      ASAVE=A
      BSAVE=B
      CSAVE=C

      DMIN=1.0D100
      DO JA=1,4
         IF (JA.EQ.1) A=ASAVE
         IF (JA.EQ.2) A=PI2-ASAVE
         IF (JA.EQ.3) A=PI-ASAVE
         IF (JA.EQ.4) A=PI+ASAVE
         SINA=SIN(A)
         COSA=COS(A)
         DO JB=1,4
            IF (JB.EQ.1) B=BSAVE
            IF (JB.EQ.2) B=PI2-BSAVE
            IF (JB.EQ.3) B=PI-BSAVE
            IF (JB.EQ.4) B=PI+BSAVE
            SINB=SIN(B)
            COSB=COS(B)
            DO JC=1,4
               IF (JC.EQ.1) C=CSAVE
               IF (JC.EQ.2) C=PI2-CSAVE
               IF (JC.EQ.3) C=PI-CSAVE
               IF (JC.EQ.4) C=PI+CSAVE
               SINC=SIN(C)
               COSC=COS(C)

               SP12=-(SINC*COSB+COSA*SINB*COSC)
               SP13=SINA*SINB
               SP22=-SINC*SINB+COSA*COSB*COSC
               SP23=-SINA*COSB
               SP32=SINA*COSC
               SP33=COSA
               SY12=SP12*HYL
               SY22=SP22*HYL
               SY32=SP32*HYL
               SZ13=SP13*HZL
               SZ23=SP23*HZL
               SZ33=SP33*HZL

               IF ( DABS(SP13T-SP13)+
     1              DABS(SP23T-SP23)+
     1              DABS(SP33T-SP33)+
     1              DABS(SY12T-SY12)+
     1              DABS(SY22T-SY22)+
     1              DABS(SY32T-SY32)+
     1              DABS(SZ13T-SZ13)+
     1              DABS(SZ23T-SZ23)+
     1              DABS(SZ33T-SZ33) .LT. DMIN) THEN
                   DMIN= DABS(SP13T-SP13)+
     1              DABS(SP23T-SP23)+
     1              DABS(SP33T-SP33)+
     1              DABS(SY12T-SY12)+
     1              DABS(SY22T-SY22)+
     1              DABS(SY32T-SY32)+
     1              DABS(SZ13T-SZ13)+
     1              DABS(SZ23T-SZ23)+
     1              DABS(SZ33T-SZ33)
                  ABEST=A
                  BBEST=B
                  CBEST=C
                  IF (DMIN.LT.1.0D-10) GOTO 20
               ENDIF
            ENDDO
         ENDDO
      ENDDO

20    A=ABEST
      B=BBEST
      C=CBEST
C     IF (DMIN.GT.0.1D0) WRITE(*,'(A,F15.5)') 'WARNING, deviation from rigid body geometry detected, best fit is ',DMIN

      RETURN
! }}}
      END
! Doxygen {{{
C>        LIMITED MEMORY BFGS METHOD FOR LARGE SCALE OPTIMIZATION
C>                          JORGE NOCEDAL
C>                        *** July 1990 ***
C>
C>        Line search removed plus small modifications, DJW 2001
C
!> \param X
! }}}
      SUBROUTINE MMYLBFGS(X,EPS,MFLAG,ENERGY,ITMAX,ITDONE,TWOD,R,R0,R1,NSIZE,NATOMS)
! Declarations {{{
      USE PORFUNCS
    
      IMPLICIT NONE
     
      INTEGER N,M
      PARAMETER (N=3,M=5)

! in/out parameters {{{

! input

      DOUBLE PRECISION X(N)
      DOUBLE PRECISION EPS,ENERGY

      INTEGER ITMAX,ITDONE,NSIZE,NATOMS
      LOGICAL MFLAG,TWOD
      
! output 

      DOUBLE PRECISION R(3,NATOMS),R0(3,NATOMS),R1(3,NATOMS)

! }}}
     
! other parameters {{{ 

      INTEGER J1,NFAIL,J,ISTAT

      DOUBLE PRECISION G(3),DIAG(N),W(N*(2*M+1)+2*M),SLENGTH,DDOT,OVERLAP
      DOUBLE PRECISION ENEW,RMS,GSAVE(3)
      DOUBLE PRECISION GNORM,STP,YS,YY,SQ,YR,BETA,DUMMY1
      INTEGER ITER,POINT,ISPT,IYPT,BOUND,NPT,CP,I,INMC,IYCN,ISCN,NDECREASE
      
      LOGICAL PTEST
! }}}
! }}}
! Subroutine body {{{
C
C  SGI appears to need this SAVE statement!
C
      SAVE

      PTEST=.TRUE.
      PTEST=.FALSE.
      NFAIL=0
      ITER=0
      ITDONE=0
      CALL DISTFUNC(X,R,R0,R1,NSIZE,ENERGY,NATOMS)
      CALL DFUNC(X,GSAVE,R,R0,R1,NSIZE,NATOMS)
      IF (TWOD) GSAVE(1)=0.0D0
      IF (TWOD) GSAVE(2)=0.0D0
      RMS=SQRT((GSAVE(1)**2+GSAVE(2)**2+GSAVE(3)**2)/3)
      DO J1=1,N
         G(J1)=GSAVE(J1)
      ENDDO

      IF (PTEST) WRITE(*,'(A,2F20.10,A,I6,A)') ' Distance and RMS force=',ENERGY,RMS,' after ',ITDONE,' LBFGS steps'
16    FORMAT(A,27X,F20.10,A)

10    CALL FLUSH(6,ISTAT)
      MFLAG=.FALSE.
      IF (RMS.LE.EPS) THEN
         MFLAG=.TRUE.
         IF (MFLAG) THEN
C           WRITE(*,'(A,F20.10)') ' Diagonal inverse Hessian elements are now ',DIAG(1)
            RETURN
         ENDIF
      ENDIF

      IF (ITDONE.EQ.ITMAX) THEN
C        WRITE(*,'(A,F20.10)') ' Diagonal inverse Hessian elements are now ',DIAG(1)
         RETURN
      ENDIF

      IF (ITER.EQ.0) THEN
         IF (N.LE.0.OR.M.LE.0) THEN
            WRITE(*,240)
 240        FORMAT(' IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
            STOP
         ENDIF
         POINT=0
         MFLAG=.FALSE.
         DO I=1,N
            DIAG(I)=1.0D0
         ENDDO
C
C     THE WORK VECTOR W IS DIVIDED AS FOLLOWS:
C     ---------------------------------------
C     THE FIRST N LOCATIONS ARE USED TO STORE THE GRADIENT AND
C         OTHER TEMPORARY INFORMATION.
C     LOCATIONS (N+1)...(N+M) STORE THE SCALARS RHO.
C     LOCATIONS (N+M+1)...(N+2M) STORE THE NUMBERS ALPHA USED
C         IN THE FORMULA THAT COMPUTES H*G.
C     LOCATIONS (N+2M+1)...(N+2M+NM) STORE THE LAST M SEARCH
C         STEPS.
C     LOCATIONS (N+2M+NM+1)...(N+2M+2NM) STORE THE LAST M
C         GRADIENT DIFFERENCES.
C
C     THE SEARCH STEPS AND GRADIENT DIFFERENCES ARE STORED IN A
C     CIRCULAR ORDER CONTROLLED BY THE PARAMETER POINT.
C
         ISPT= N+2*M
         IYPT= ISPT+N*M
C
C  NR step for diagonal inverse Hessian
C
         DO I=1,N
            W(ISPT+I)= -G(I)*DIAG(I)
            W(I)= -G(I)*DIAG(I)
         ENDDO
         GNORM= DSQRT(DDOT(N,G,1,G,1))
C
C  Make the first guess for the step length cautious.
C
         STP=MIN(1.0D0/GNORM,GNORM)
C        STP=1.0D0
      ELSE 
         BOUND=ITER
         IF (ITER.GT.M) BOUND=M
         YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
         IF (YS.EQ.0.0D0) YS=1.0D0
C
C  Update estimate of diagonal inverse Hessian elements
C  We have to divide by YS and YY, so make sure they
C  aren't zero!
C
         YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
         IF (YY.EQ.0.0D0) YY=1.0D0
         DUMMY1= YS/YY
C        DUMMY1=ABS(YS/YY)
         DO I=1,N
            DIAG(I)=DUMMY1
         ENDDO
C
C     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
C     "Updating quasi-Newton matrices with limited storage",
C     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
C     ---------------------------------------------------------
C
         CP= POINT
         IF (POINT.EQ.0) CP=M
         W(N+CP)= 1.0D0/YS
         DO I=1,N
            W(I)= -G(I)
         ENDDO
         CP= POINT
         DO I= 1,BOUND
            CP=CP-1
            IF (CP.EQ. -1)CP=M-1
            SQ= DDOT(N,W(ISPT+CP*N+1),1,W,1)
            INMC=N+M+CP+1
            IYCN=IYPT+CP*N
            W(INMC)= W(N+CP+1)*SQ
            CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
         ENDDO
        
         DO I=1,N
            W(I)=DIAG(I)*W(I)
         ENDDO

         DO I=1,BOUND
            YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
            BETA= W(N+CP+1)*YR
            INMC=N+M+CP+1
            BETA= W(INMC)-BETA
            ISCN=ISPT+CP*N
            CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
            CP=CP+1
            IF (CP.EQ.M) CP=0
         ENDDO
         STP=1.0D0
      ENDIF
C
C  Store the new search direction
C
      IF (ITER.GT.0) THEN
         DO I=1,N
            W(ISPT+POINT*N+I)= W(I)
         ENDDO
      ENDIF
      IF (TWOD) W(ISPT+POINT*N+1)=0.0D0
      IF (TWOD) W(ISPT+POINT*N+2)=0.0D0

      OVERLAP=DDOT(N,G,1,W,1)/SQRT(DDOT(N,G,1,G,1)*DDOT(N,W,1,W,1))
C     PRINT*,'OVERLAP,DIAG(1)=',OVERLAP,DIAG(1)
C     PRINT*,'G . G=',DDOT(N,G,1,G,1)
C     PRINT*,'W . W=',DDOT(N,W,1,W,1)
      IF (OVERLAP.GT.0.0D0) THEN
         IF (PTEST) PRINT*,'Search direction has positive projection onto gradient - reversing step'
C        IF (PTEST) PRINT*,'Search direction has positive projection onto gradient - reset'
         DO I=1,N
            W(ISPT+POINT*N+I)= -W(I)
         ENDDO
C        ITER=0
C        GOTO 10
      ENDIF

      DO I=1,N
         W(I)=G(I)
      ENDDO
      SLENGTH=0.0D0
      DO J1=1,N
         SLENGTH=SLENGTH+W(ISPT+POINT*N+J1)**2
      ENDDO
      SLENGTH=SQRT(SLENGTH)
C     IF (STP*SLENGTH.GT.MAXMBFGS) STP=MAXMBFGS/SLENGTH
      IF (STP*SLENGTH.GT.0.2D0) STP=0.2D0/SLENGTH
C
C  We now have the proposed step.
C
      DO J1=1,N
         X(J1)=X(J1)+STP*W(ISPT+POINT*N+J1)
      ENDDO 
      NDECREASE=0
20    CALL DISTFUNC(X,R,R0,R1,NSIZE,ENEW,NATOMS)

      IF (ENEW-ENERGY.LE.1.0D-3) THEN
         ITER=ITER+1
         ITDONE=ITDONE+1
         ENERGY=ENEW
         IF (PTEST) WRITE(*,'(A,2F20.10,A,I6,A,F13.10)') ' Distance and RMS force=',ENERGY,RMS,' after ',ITDONE,
     1           ' LBFGS steps, step:',STP*SLENGTH
      ELSE 
C
C  Energy increased - try again with a smaller step size
C
C        IF (STP*SLENGTH.LT.1.0D-10) THEN
         IF (NDECREASE.GT.5) THEN
            PRINT*,' in mmylbfgs LBFGS step cannot find a lower distance try again'
            NFAIL=NFAIL+1
            IF (NFAIL.GT.20) THEN
               PRINT*,' Too many failures - give up'
               RETURN
            ENDIF
            DO J1=1,N
               X(J1)=X(J1)-STP*W(ISPT+POINT*N+J1)
            ENDDO 
            GOTO 30
         ENDIF
         DO J1=1,N
            X(J1)=X(J1)-0.9*STP*W(ISPT+POINT*N+J1)
         ENDDO 
         NDECREASE=NDECREASE+1
         STP=STP/10.0D0
         IF (PTEST) 
     1    WRITE(*,'(A,F19.10,A,F16.10,A,F15.8)') ' Distance increased from ',ENERGY,' to ',ENEW,' decreasing step to ',STP*SLENGTH
         GOTO 20
      ENDIF
C
C  Reset comparison geometry to R1 and rotation angles to zero.
C
      X(1)=0.0D0
      X(2)=0.0D0
      X(3)=0.0D0
      DO I=1,NSIZE
         DO J=1,3
            R0(J,I)=R1(J,I)
         ENDDO
      ENDDO

      CALL DFUNC(X,GSAVE,R,R0,R1,NSIZE,NATOMS)
      IF (TWOD) GSAVE(1)=0.0D0
      IF (TWOD) GSAVE(2)=0.0D0
      RMS=SQRT((GSAVE(1)**2+GSAVE(2)**2+GSAVE(3)**2)/3)
      DO J1=1,N
         G(J1)=GSAVE(J1)
      ENDDO
C
C     Compute the new step and gradient change
C
30    NPT=POINT*N
      DO I=1,N
         W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
         W(IYPT+NPT+I)= G(I)-W(I)
      ENDDO
      POINT=POINT+1
      IF (POINT.EQ.M) POINT=0
      GOTO 10
! }}}
      RETURN
      END
