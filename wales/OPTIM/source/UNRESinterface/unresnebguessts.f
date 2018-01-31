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
C Routine to generate neb images for unres by interpolating
C between different internal coordinates 
C Changed to be called from sat's newneb -- make finish passed array
C
      SUBROUTINE UNRESDIHENEB(Q,FINISH,POINTS)
      USE COMMONS
      USE KEY
C     USE MODNEB
      USE MODUNRES
      use KeyNEB,only: Nimage
      IMPLICIT NONE
C
      DOUBLE PRECISION ANGLE,TWISTFRAC,POINTS(3*NATOMS*NIMAGE),PI
      INTEGER I1,J1,J2,NM
      REAL*8 DUMMY(3*NATOMS),Q(3*NATOMS),FINISH(3*NATOMS),DIFFPP,DIST,FINPPSANGLE(NINTS),QPPSANGLE(NINTS)
      PARAMETER (PI=3.141592653589793D0)

      DO I1=1,nres
         c(1,I1)=Q(6*(I1-1)+1)
         c(2,I1)=Q(6*(I1-1)+2)
         c(3,I1)=Q(6*(I1-1)+3)
         c(1,I1+nres)=Q(6*(I1-1)+4)
         c(2,I1+nres)=Q(6*(I1-1)+5)
         c(3,I1+nres)=Q(6*(I1-1)+6)
      END DO
      CALL UPDATEDC
      CALL int_from_cart(.true.,.false.)
      CALL geom_to_var(nvaru,QPPSANGLE)

      DO I1=1,nres
         c(1,I1)=FINISH(6*(I1-1)+1)
         c(2,I1)=FINISH(6*(I1-1)+2)
         c(3,I1)=FINISH(6*(I1-1)+3)
         c(1,I1+nres)=FINISH(6*(I1-1)+4)
         c(2,I1+nres)=FINISH(6*(I1-1)+5)
         c(3,I1+nres)=FINISH(6*(I1-1)+6)
      END DO
      CALL UPDATEDC
      CALL int_from_cart(.true.,.false.)
      CALL geom_to_var(nvaru,FINPPSANGLE)

C DAE's comment...not relevant for unres all-internal implementation.
C We need to keep our images along the path between q and fin. Their alignment will
C be messed up by rebuilding. Therefore rebuild starting point q (before any dihes have
C been changed), and put these rebuilt coordinates in maximal alignment with the proper q.
C using mind. This will store the rotation matrix required, which we then apply to all the 
C rebuilt images

      TWISTFRAC=1.D0/(Nimage+1)   ! Incremental, so TWISTFRAC is always the same

      DO J1=1,Nimage

         DO I1=1,nvaru
            DIFFPP = FINPPSANGLE(I1) - QPPSANGLE(I1)
C
C next two lines are meant to ensure that you always interpolate
C along the shortest distance between the dihedral angles.
C
            IF (DIFFPP.LT.-PI) DIFFPP = DIFFPP+2.0D0*PI
            IF (DIFFPP.GT.PI) DIFFPP = DIFFPP-2.0D0*PI
            ANGLE=TWISTFRAC*DIFFPP
            QPPSANGLE(I1)=QPPSANGLE(I1)+ANGLE 
         ENDDO
 
         CALL var_to_geom(nvaru,QPPSANGLE)
         CALL chainbuild
         DO I1=1,nres
            DUMMY(6*(I1-1)+1)=c(1,I1)
            DUMMY(6*(I1-1)+2)=c(2,I1)
            DUMMY(6*(I1-1)+3)=c(3,I1)
            DUMMY(6*(I1-1)+4)=c(1,I1+nres)
            DUMMY(6*(I1-1)+5)=c(2,I1+nres)
            DUMMY(6*(I1-1)+6)=c(3,I1+nres)
         ENDDO

         DO I1=1,3*NATOMS
            POINTS(I1+NOPT*(J1-1))=DUMMY(I1)
         END DO
         
      ENDDO

      RETURN
      END

C
C Routine to calculate distance in spring force for neb, when using the distance in dihedral angle space
C DISTDIHE = |QC - QB| - |QB - QA|
C
      SUBROUTINE UNRESGETDIHEDIST(DIHEDIST,QA,QB,QC)
      USE COMMONS
      USE MODUNRES
      IMPLICIT NONE

      INTEGER I1
      REAL*8 QA(3*NATOMS),QB(3*NATOMS),QC(3*NATOMS)
      REAL*8 DIFFPP,DIHEDIST,DISTCB,DISTBA
      REAL*8 QAPPSANGLE(NINTS),QBPPSANGLE(NINTS),QCPPSANGLE(NINTS)
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)

      DO I1=1,nres
         c(1,I1)=QA(6*(I1-1)+1)
         c(2,I1)=QA(6*(I1-1)+2)
         c(3,I1)=QA(6*(I1-1)+3)
         c(1,I1+nres)=QA(6*(I1-1)+4)
         c(2,I1+nres)=QA(6*(I1-1)+5)
         c(3,I1+nres)=QA(6*(I1-1)+6)
      END DO
      CALL UPDATEDC
      CALL int_from_cart(.true.,.false.)
c     CALL chainbuild
      CALL geom_to_var(nvaru,QAPPSANGLE)

      DO I1=1,nres
         c(1,I1)=QB(6*(I1-1)+1)
         c(2,I1)=QB(6*(I1-1)+2)
         c(3,I1)=QB(6*(I1-1)+3)
         c(1,I1+nres)=QB(6*(I1-1)+4)
         c(2,I1+nres)=QB(6*(I1-1)+5)
         c(3,I1+nres)=QB(6*(I1-1)+6)
      END DO
      CALL UPDATEDC
      CALL int_from_cart(.true.,.false.)
c     CALL chainbuild
      CALL geom_to_var(nvaru,QBPPSANGLE)

      DO I1=1,nres
         c(1,I1)=QC(6*(I1-1)+1)
         c(2,I1)=QC(6*(I1-1)+2)
         c(3,I1)=QC(6*(I1-1)+3)
         c(1,I1+nres)=QC(6*(I1-1)+4)
         c(2,I1+nres)=QC(6*(I1-1)+5)
         c(3,I1+nres)=QC(6*(I1-1)+6)
      END DO
      CALL UPDATEDC
      CALL int_from_cart(.true.,.false.)
c     CALL chainbuild
      CALL geom_to_var(nvaru,QCPPSANGLE)

      DISTCB=0.D0
      DISTBA=0.D0
      DO I1=1,nvaru
C jmc why not adjust if > pi etc?
         DIFFPP=QCPPSANGLE(I1)-QBPPSANGLE(I1)
         IF (DIFFPP.LT.-PI) DIFFPP = DIFFPP+2.0D0*PI
         IF (DIFFPP.GT.PI) DIFFPP = DIFFPP-2.0D0*PI

         DISTCB=DISTCB+ABS(DIFFPP)

         DIFFPP=QBPPSANGLE(I1)-QAPPSANGLE(I1)
         IF (DIFFPP.LT.-PI) DIFFPP = DIFFPP+2.0D0*PI
         IF (DIFFPP.GT.PI) DIFFPP = DIFFPP-2.0D0*PI
         DISTBA=DISTBA+ABS(DIFFPP)
      ENDDO

      DISTCB=DISTCB/nvaru
      DISTBA=DISTBA/nvaru

      DIHEDIST=DISTCB-DISTBA

      RETURN
      END

