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
      SUBROUTINE UNRSETZSYMATMASS
      USE COMMONS
      USE MODUNRES
      IMPLICIT NONE

      INTEGER I1

      ALLOCATE(ATMASS(NATOMS))
      DO I1=1,nres
         ZSYM(2*I1-1)='C'
         ZSYM(2*I1)='C'
         ATMASS(2*I1-1)=MASSES(1)
         ATMASS(2*I1)=MASSES(itype(I1)+1)
      ENDDO
C jmc make a masses file for pathsample...
c     OPEN (UNIT=78,FILE='mass',STATUS='UNKNOWN')
c     DO I1=1,nres
c        WRITE(78,'(A3,F5.1)') ZSYM(2*I1-1),ATMASS(2*I1-1)
c        WRITE(78,'(A3,F5.1)') ZSYM(2*I1),ATMASS(2*I1)
c     END DO
c     CLOSE (UNIT=78)
 
      RETURN
      END


      SUBROUTINE UENERGY(X,GRAD,etot,GRADT,SECT)
      USE COMMONS
      USE MODUNRES
      IMPLICIT NONE

      DOUBLE PRECISION GRAD(3*NATOMS), X(3*NATOMS)
      REAL*8 evdw,evdw1,evdw2,ees,ebe,esc,etors,ehpb,edihcnstr,ecorr,etot
      LOGICAL GRADT,SECT
      INTEGER icall
      common /srutu/ icall

      nfl=0
      icg=1
      icall=1

      CALL zerograd
      CALL etotal(evdw,evdw1,evdw2,ees,ebe,esc,etors,ehpb,edihcnstr,ecorr,etot)

      CALL INTGRAD(GRAD)
C jmc INTGRAD takes grad in 'carts' (i.e. in DC,DX: Ca-Ca and sc-sc vectors) and converts 
C it to a gradient in internals, which is passed back to POTENTIAL etc.

C jmc testing stuff!!!
C     CALL zerograd
C     CALL check_cartgrad
C     CALL check_ecart
C     CALL check_eint

C     PRINT *,'UNRES energy :'
C     PRINT *,etot

C     CALL cartprint
C     CALL int_from_cart(.true.,.true.)

      RETURN
      END SUBROUTINE UENERGY

      SUBROUTINE INTGRAD(gint)
      USE COMMONS
      USE MODUNRES
      IMPLICIT NONE ! note not present in UNOPTIM.2.3

      DOUBLE PRECISION X(NINTS),G(NINTS)
      INTEGER uiparm(1),icall,nf
      DOUBLE PRECISION URPARM(1)
      EXTERNAL fdum
      DOUBLE PRECISION GINT(3*NATOMS)
      INTEGER J1
      COMMON /srutu/ icall

      call geom_to_var(nvaru,x)
C     call var_to_geom(nvaru,x)
C     call chainbuild

      nf=1
      nfl=3
      call gradient(nvaru,x,nf,g,uiparm,urparm,fdum)

      gint=0.0D0

      gint(1:nvaru)=g(1:nvaru)

      RETURN
      END SUBROUTINE INTGRAD

      SUBROUTINE UPDATEDC
      USE COMMONS
      USE MODUNRES
      IMPLICIT NONE

      INTEGER J1

      DC=0.0D0
      DC_NORM=0.0D0

      DO J1=1,nres-1
         DC(1,J1)=c(1,J1+1)-c(1,J1)
         DC(2,J1)=c(2,J1+1)-c(2,J1)
         DC(3,J1)=c(3,J1+1)-c(3,J1)
         DC_NORM(1,J1)=DC(1,J1)/vbl
         DC_NORM(2,J1)=DC(2,J1)/vbl
         DC_NORM(3,J1)=DC(3,J1)/vbl
C        PRINT *,MYDC(1,J1)/DC(1,J1),MYDC(2,J1)/DC(2,J1),MYDC(3,J1)/DC(3,J1)
C        PRINT *,MYDCNORM(1,J1)/dc_norm(1,J1),MYDCNORM(2,J1)/dc_norm(2,J1),MYDCNORM(3,J1)/dc_norm(3,J1)
      END DO

      DO J1=nres+1,2*nres
         DC(1,J1)=c(1,J1)-c(1,J1-nres)
         DC(2,J1)=c(2,J1)-c(2,J1-nres)
         DC(3,J1)=c(3,J1)-c(3,J1-nres)
         IF (itype(J1-nres).EQ.10) THEN   ! i.e. glycine
            DC_NORM(1,J1)=DC(1,J1)
            DC_NORM(2,J1)=DC(2,J1)
            DC_NORM(3,J1)=DC(3,J1)
         ELSE
            DC_NORM(1,J1)=DC(1,J1)/dsc(itype(J1-nres))
            DC_NORM(2,J1)=DC(2,J1)/dsc(itype(J1-nres))
            DC_NORM(3,J1)=DC(3,J1)/dsc(itype(J1-nres))
         ENDIF
C        PRINT *,MYDC(1,J1)/DC(1,J1),MYDC(2,J1)/DC(2,J1),MYDC(3,J1)/DC(3,J1)
C        PRINT *,MYDCNORM(1,J1)/dc_norm(1,J1),MYDCNORM(2,J1)/dc_norm(2,J1),MYDCNORM(3,J1)/dc_norm(3,J1)
      END DO
 
      RETURN
      END SUBROUTINE UPDATEDC
