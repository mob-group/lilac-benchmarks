!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of OPTIM.
!
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!   Finds the minimum distance between two geometries.
!   Geometry in RA should not change. RB is returned as the
!   closest geometry to RA if PRESERVET is .FALSE.
!
!   New analytic method based on quaterions from
!   Kearsley, Acta Cryst. A, 45, 208-210, 1989.
!
! jmc As long as zsym isn't 'W' (in which case mind does something special) mind
! doesn't care what atomic symbol we give it.
!
SUBROUTINE NEWMINDIST(RA,RB,NATOMS,DIST,BULKT,TWOD,ZUSE,PRESERVET,ANGLEAXIS,DEBUG,RMAT)
use rigidbodymod
USE COMMONS,ONLY : FREEZE, BOXLX, BOXLY, BOXLZ, RBAAT, PULLT, TRAPT, USERPOTT, GTHOMSONT, GTHOMMET, MKTRAPT, NOTRANSROTT ! hk286
IMPLICIT NONE
INTEGER J1, NATOMS, NSIZE, INFO, JINFO, JMIN
DOUBLE PRECISION RA(3*NATOMS), RB(3*NATOMS), DIST, QMAT(4,4), XM, YM, ZM, XP, YP, ZP, OVEC(3), H1VEC(3), H2VEC(3), &
  &              DIAG(4), TEMPA(9*NATOMS), RMAT(3,3), MINV, Q1, Q2, Q3, Q4, CMXA, CMYA, CMZA, CMXB, CMYB, CMZB, &
  &              NCMXB, NCMYB, NCMZB, ROTMAT(3,3), OMEGATOT(3,3)
DOUBLE PRECISION :: COCA(3), COCB(3)
DOUBLE PRECISION, ALLOCATABLE :: XA(:), XB(:), TEMPCOORDS(:,:)
COMMON /MINDOM/ ROTMAT, OMEGATOT
LOGICAL, INTENT(IN) :: BULKT, TWOD, PRESERVET, DEBUG, ANGLEAXIS
LOGICAL UPRETURN
CHARACTER(LEN=5) ZUSE
INTEGER NCIT
DOUBLE PRECISION XSHIFT, YSHIFT, ZSHIFT, XSHIFTNEW, YSHIFTNEW, ZSHIFTNEW

IF (RBAAT) THEN

   CALL RBMINDIST(RA,RB,NATOMS,DIST,RMAT,DEBUG)
   RETURN

ELSEIF (MKTRAPT) THEN
   DIST=0.0D0
   DO J1=1,3*NATOMS
      DIST=DIST+(RA(J1)-RB(J1))**2
   ENDDO
   DIST=SQRT(DIST)
   RETURN

ELSEIF (USERPOTT) THEN
!
! Substitute a real routine for calculating the distance between two
! configurations.
!
   call quaternionMatch2 (NATOMS, RA, RB, dist)
   dist = sqrt(dist)
   RMAT(1:3,1:3)=0.0D0 ! rotation matrix is the identity
   RMAT(1,1)=1.0D0; RMAT(2,2)=1.0D0; RMAT(3,3)=1.0D0
   UPRETURN=.TRUE.
   If (UPRETURN) RETURN
ENDIF



! Convert rigid body coordinates to Cartesians for rigid bodies.

IF (ANGLEAXIS) THEN
   NSIZE = RBPOTENTIAL%NPHYSICALSITES*(NATOMS/2)
   ALLOCATE(XA(3*NSIZE),XB(3*NSIZE))

   CALL SYSTEMTOCARTESIANS(NATOMS/2,RA,XA)
   CALL SYSTEMTOCARTESIANS(NATOMS/2,RB,XB)

ELSE IF (ZUSE(1:1).EQ.'W') THEN
   NSIZE=3*(NATOMS/2)
   ALLOCATE(XA(3*NSIZE),XB(3*NSIZE))

   DO J1=1,NATOMS/2
      CALL CONVERT(RA(3*J1-2),RA(3*J1-1),RA(3*J1), &
      RA(3*(NATOMS/2+J1)-2),RA(3*(NATOMS/2+J1)-1),RA(3*(NATOMS/2+J1)), &
      XA(9*J1-8:9*J1-6),XA(9*J1-5:9*J1-3),XA(9*J1-2:9*J1))

      CALL CONVERT(RB(3*J1-2),RB(3*J1-1),RB(3*J1), &
      RB(3*(NATOMS/2+J1)-2),RB(3*(NATOMS/2+J1)-1),RB(3*(NATOMS/2+J1)), &
      XB(9*J1-8:9*J1-6),XB(9*J1-5:9*J1-3),XB(9*J1-2:9*J1))
   ENDDO
ELSEIF ((PULLT.OR.(GTHOMSONT.AND.(GTHOMMET<5))).AND.(.NOT.BULKT)) THEN
!  ALLOCATE(XA(3*(NATOMS/2)*number of sites,XB(3*(NATOMS/2)*number of sites))
!  NSIZE=(NATOMS/2)*number of sites
!  PRINT '(A)',' newmindist> New quaternion procedure not yet coded for flatland'
! There is one unknown angle, so this should be trivial!'
   CALL MINDIST(RA,RB,NATOMS,DIST,BULKT,TWOD,ZUSE,PRESERVET)
   RMAT(1:3,1:3)=OMEGATOT(1:3,1:3)
   RETURN
ELSE
   ALLOCATE(XA(3*NATOMS),XB(3*NATOMS))
   NSIZE=NATOMS
   XA(1:3*NATOMS)=RA(1:3*NATOMS)
   XB(1:3*NATOMS)=RB(1:3*NATOMS)
ENDIF
IF (FREEZE) THEN
   DIST=0.0D0
   DO J1=1,3*NSIZE
      DIST=DIST+(RA(J1)-RB(J1))**2
   ENDDO
   DIST=SQRT(DIST)
   IF (ALLOCATED(XA)) DEALLOCATE(XA,XB)
   RETURN
ENDIF
IF (BULKT) THEN
   DO J1=1,NSIZE
      XA(3*(J1-1)+1)=XA(3*(J1-1)+1)-BOXLX*NINT(XA(3*(J1-1)+1)/BOXLX)
      XA(3*(J1-1)+2)=XA(3*(J1-1)+2)-BOXLY*NINT(XA(3*(J1-1)+2)/BOXLY)
      IF (.NOT.TWOD) XA(3*(J1-1)+3)=XA(3*(J1-1)+3)-BOXLZ*NINT(XA(3*(J1-1)+3)/BOXLZ)
   ENDDO
   DO J1=1,NSIZE
      XB(3*(J1-1)+1)=XB(3*(J1-1)+1)-BOXLX*NINT(XB(3*(J1-1)+1)/BOXLX)
      XB(3*(J1-1)+2)=XB(3*(J1-1)+2)-BOXLY*NINT(XB(3*(J1-1)+2)/BOXLY)
      IF (.NOT.TWOD) XB(3*(J1-1)+3)=XB(3*(J1-1)+3)-BOXLZ*NINT(XB(3*(J1-1)+3)/BOXLZ)
   ENDDO
ENDIF
!
! Move centre of coordinates of XA and XB to the origin.
!
CMXA=0.0D0; CMYA=0.0D0; CMZA=0.0D0
! hk286
IF (.NOT. GTHOMSONT) THEN
   DO J1=1,NSIZE
      CMXA=CMXA+XA(3*(J1-1)+1)
      CMYA=CMYA+XA(3*(J1-1)+2)
      CMZA=CMZA+XA(3*(J1-1)+3)
   ENDDO
   CMXA=CMXA/NSIZE; CMYA=CMYA/NSIZE; CMZA=CMZA/NSIZE
   DO J1=1,NSIZE
      XA(3*(J1-1)+1)=XA(3*(J1-1)+1)-CMXA
      XA(3*(J1-1)+2)=XA(3*(J1-1)+2)-CMYA
      XA(3*(J1-1)+3)=XA(3*(J1-1)+3)-CMZA
   ENDDO
   CMXB=0.0D0; CMYB=0.0D0; CMZB=0.0D0
   DO J1=1,NSIZE
      CMXB=CMXB+XB(3*(J1-1)+1)
      CMYB=CMYB+XB(3*(J1-1)+2)
      CMZB=CMZB+XB(3*(J1-1)+3)
   ENDDO
   CMXB=CMXB/NSIZE; CMYB=CMYB/NSIZE; CMZB=CMZB/NSIZE
   DO J1=1,NSIZE
      XB(3*(J1-1)+1)=XB(3*(J1-1)+1)-CMXB
      XB(3*(J1-1)+2)=XB(3*(J1-1)+2)-CMYB
      XB(3*(J1-1)+3)=XB(3*(J1-1)+3)-CMZB
   ENDDO
ENDIF
!
! Deal with TWO2 non-bulk here.
!
IF (TWOD.AND.(.NOT.BULKT)) THEN
   CALL MINDIST(XA,XB,NATOMS,DIST,BULKT,TWOD,ZUSE,PRESERVET)
   RMAT(1:3,1:3)=OMEGATOT(1:3,1:3)
!
!  To align the structures for minpermdist we use the orientation in XA
!  and translate the B coordinates to the centre of coordinates of RA.
!
   IF (.NOT.PRESERVET) THEN
!      IF (DEBUG) PRINT '(3(A,G20.10), A)', 'XSHIFT = ', XSHIFT, '; YSHIFT = ', YSHIFT, '; ZSHIFT = ', ZSHIFT, ';'
! jwrm2> It seems that the intention should be to return the rotated RB at the centre of mass of RA. 
!        This would not seem to do that.
!      DO J1=1,NATOMS
!         RB(3*(J1-1)+1)=XB(3*(J1-1)+1)-CMXB+CMXA+XSHIFT
!         RB(3*(J1-1)+2)=XB(3*(J1-1)+2)-CMYB+CMYA+YSHIFT
!         RB(3*(J1-1)+3)=XB(3*(J1-1)+3)-CMZB+CMZA+ZSHIFT
!      ENDDO
! jwrm2> Try this instead
      DO J1=1,NATOMS
         RB(3*(J1-1)+1)=XB(3*(J1-1)+1)+CMXA
         RB(3*(J1-1)+2)=XB(3*(J1-1)+2)+CMYA
         RB(3*(J1-1)+3)=XB(3*(J1-1)+3)+CMZA
      ENDDO
   ENDIF
   RETURN
ENDIF

IF (NOTRANSROTT) THEN
   DIST=0.0D0
   DO J1=1,NSIZE
      DIST=DIST + (XA(3*(J1-1)+1)-XB(3*(J1-1)+1))**2 &
   &            + (XA(3*(J1-1)+2)-XB(3*(J1-1)+2))**2 &
   &            + (XA(3*(J1-1)+3)-XB(3*(J1-1)+3))**2
   ENDDO
   DIST=SQRT(DIST)

   RMAT(1:3,1:3)=0.0D0 ! rotation matrix is the identity
   RMAT(1,1)=1.0D0; RMAT(2,2)=1.0D0; RMAT(3,3)=1.0D0
   RETURN
ENDIF

XSHIFT=0.0D0; YSHIFT=0.0D0; ZSHIFT=0.0D0
NCIT=0
IF (BULKT) THEN 
! 1  NCIT=NCIT+1
!    IF (NCIT.GT.1000) THEN
!       PRINT '(A)','inertia> WARNING - iterative calculation of centre of mass shift did not converge'
!    ENDIF
!    XSHIFTNEW=0.0D0
!    YSHIFTNEW=0.0D0
!    ZSHIFTNEW=0.0D0
!    DO J1=1,NSIZE
!       XSHIFTNEW=XSHIFTNEW + XA(3*(J1-1)+1)-XB(3*(J1-1)+1) - BOXLX*NINT((XA(3*(J1-1)+1)-XB(3*(J1-1)+1)-XSHIFT)/BOXLX)
!       YSHIFTNEW=YSHIFTNEW + XA(3*(J1-1)+2)-XB(3*(J1-1)+2) - BOXLY*NINT((XA(3*(J1-1)+2)-XB(3*(J1-1)+2)-YSHIFT)/BOXLY)
!       IF (.NOT.TWOD) ZSHIFTNEW=ZSHIFTNEW + XA(3*(J1-1)+3)-XB(3*(J1-1)+3) - BOXLZ*NINT((XA(3*(J1-1)+3)-XB(3*(J1-1)+3)-ZSHIFT)/BOXLZ)
!    ENDDO
!    XSHIFTNEW=XSHIFTNEW/NSIZE; YSHIFTNEW=YSHIFTNEW/NSIZE; ZSHIFTNEW=ZSHIFTNEW/NSIZE
!    IF ((ABS(XSHIFTNEW-XSHIFT).GT.1.0D-6).OR.(ABS(YSHIFTNEW-YSHIFT).GT.1.0D-6).OR.(ABS(ZSHIFTNEW-ZSHIFT).GT.1.0D-6)) THEN
! !     IF (DEBUG) PRINT '(A,I6,6F15.7)',' newmindist> ',NCIT,XSHIFTNEW,YSHIFTNEW,ZSHIFTNEW,XSHIFT,YSHIFT,ZSHIFT
!       XSHIFT=0.05D0*XSHIFT+0.95D0*XSHIFTNEW
!       YSHIFT=0.05D0*YSHIFT+0.95D0*YSHIFTNEW
!       ZSHIFT=0.05D0*ZSHIFT+0.95D0*ZSHIFTNEW
!       GOTO 1
!    ENDIF
! !  IF (DEBUG) PRINT '(A,I6,3F15.7)',' newmindist> coordinate shift converged. Cycles and values: ',NCIT,XSHIFT,YSHIFT,ZSHIFT
!    XSHIFT=XSHIFTNEW
!    YSHIFT=YSHIFTNEW
!    ZSHIFT=ZSHIFTNEW
! Actually, the iterative solution seems to be worse than simply putting the centre of mass
! at the origin.
!
   DIST=0.0D0   
   DO J1=1,NSIZE
      DIST=DIST + (XA(3*(J1-1)+1)-XB(3*(J1-1)+1)-XSHIFT - BOXLX*NINT((XA(3*(J1-1)+1)-XB(3*(J1-1)+1)-XSHIFT)/BOXLX))**2 &
   &            + (XA(3*(J1-1)+2)-XB(3*(J1-1)+2)-YSHIFT - BOXLY*NINT((XA(3*(J1-1)+2)-XB(3*(J1-1)+2)-YSHIFT)/BOXLY))**2 
      IF (.NOT.TWOD) DIST=DIST &
   &            + (XA(3*(J1-1)+3)-XB(3*(J1-1)+3)-ZSHIFT - BOXLZ*NINT((XA(3*(J1-1)+3)-XB(3*(J1-1)+3)-ZSHIFT)/BOXLZ))**2
   ENDDO
   DIST=SQRT(DIST)

   RMAT(1:3,1:3)=0.0D0 ! rotation matrix is the identity
   RMAT(1,1)=1.0D0; RMAT(2,2)=1.0D0; RMAT(3,3)=1.0D0
ELSE
!
!  The formula below is not invariant to overall translation because XP, YP, ZP
!  involve a sum of coordinates! We need to have XA and XB coordinate centres both
!  at the origin!!
!
   QMAT(1:4,1:4)=0.0D0
   DO J1=1,NSIZE
      XM=XA(3*(J1-1)+1)-XB(3*(J1-1)+1)
      YM=XA(3*(J1-1)+2)-XB(3*(J1-1)+2)
      ZM=XA(3*(J1-1)+3)-XB(3*(J1-1)+3)
      XP=XA(3*(J1-1)+1)+XB(3*(J1-1)+1)
      YP=XA(3*(J1-1)+2)+XB(3*(J1-1)+2)
      ZP=XA(3*(J1-1)+3)+XB(3*(J1-1)+3)
      QMAT(1,1)=QMAT(1,1)+XM**2+YM**2+ZM**2
      QMAT(1,2)=QMAT(1,2)+YP*ZM-YM*ZP
      QMAT(1,3)=QMAT(1,3)+XM*ZP-XP*ZM
      QMAT(1,4)=QMAT(1,4)+XP*YM-XM*YP
      QMAT(2,2)=QMAT(2,2)+YP**2+ZP**2+XM**2
      QMAT(2,3)=QMAT(2,3)+XM*YM-XP*YP
      QMAT(2,4)=QMAT(2,4)+XM*ZM-XP*ZP
      QMAT(3,3)=QMAT(3,3)+XP**2+ZP**2+YM**2
      QMAT(3,4)=QMAT(3,4)+YM*ZM-YP*ZP
      QMAT(4,4)=QMAT(4,4)+XP**2+YP**2+ZM**2
   ENDDO
   QMAT(2,1)=QMAT(1,2); QMAT(3,1)=QMAT(1,3); QMAT(3,2)=QMAT(2,3); QMAT(4,1)=QMAT(1,4); QMAT(4,2)=QMAT(2,4); QMAT(4,3)=QMAT(3,4)

   CALL DSYEV('V','U',4,QMAT,4,DIAG,TEMPA,9*NATOMS,INFO)
   IF (INFO.NE.0) PRINT '(A,I6,A)','newmindist> WARNING - INFO=',INFO,' in DSYEV'

   MINV=1.0D100
   DO J1=1,4
      IF (DIAG(J1).LT.MINV) THEN
         JMIN=J1
         MINV=DIAG(J1)
      ENDIF
   ENDDO
   IF (MINV.LT.0.0D0) THEN
      IF (ABS(MINV).LT.1.0D-6) THEN
         MINV=0.0D0
      ELSE
         PRINT '(A,G20.10,A)','newmindist> WARNING MINV is ',MINV,' change to absolute value'
         MINV=-MINV
      ENDIF
   ENDIF
   DIST=SQRT(MINV)

! IF (DEBUG) PRINT '(A,G20.10,A,I6)','newmindist> minimum residual is ',DIAG(JMIN),' for eigenvector ',JMIN
   Q1=QMAT(1,JMIN); Q2=QMAT(2,JMIN); Q3=QMAT(3,JMIN); Q4=QMAT(4,JMIN)
!
! RMAT will contain the matrix that maps RB onto the best correspondence with RA
!
   RMAT(1,1)=Q1**2+Q2**2-Q3**2-Q4**2
   RMAT(1,2)=2*(Q2*Q3+Q1*Q4)
   RMAT(1,3)=2*(Q2*Q4-Q1*Q3)
   RMAT(2,1)=2*(Q2*Q3-Q1*Q4)
   RMAT(2,2)=Q1**2+Q3**2-Q2**2-Q4**2
   RMAT(2,3)=2*(Q3*Q4+Q1*Q2)
   RMAT(3,1)=2*(Q2*Q4+Q1*Q3)
   RMAT(3,2)=2*(Q3*Q4-Q1*Q2)
   RMAT(3,3)=Q1**2+Q4**2-Q2**2-Q3**2
ENDIF
!
! Superpose RB if required
! For the angle-axis formulation we are not translating RB but simply
! setting RB equal to XB, so the origin will be the centre of coordinates
! including all the rigid body sites. This might need to be changed.
!
IF (.NOT.PRESERVET.OR.NOTRANSROTT) THEN
   CALL NEWROTGEOM(NSIZE,XB,RMAT,CMXA,CMYA,CMZA)

   IF (ANGLEAXIS) THEN
      ALLOCATE(TEMPCOORDS(NSIZE,3))
      TEMPCOORDS = RESHAPE(XB,SHAPE=(/NSIZE,3/),ORDER=(/2,1/))
      CALL SYSTEMTOAA(NATOMS/2,TEMPCOORDS,RB)
      DEALLOCATE(TEMPCOORDS)
   ELSE IF (ZUSE(1:1).EQ.'W') THEN
      DO J1=1,NATOMS/2
         CALL CONVERT2(XB(9*J1-8:9*J1-6),XB(9*J1-5:9*J1-3),XB(9*J1-2:9*J1), &
                       RB(3*J1-2),RB(3*J1-1),RB(3*J1), &
                       RB(3*(J1+NATOMS/2)-2),RB(3*(J1+NATOMS/2)-1),RB(3*(J1+NATOMS/2)))
      ENDDO
!
!  Preserve centre of mass for TRAPT, but not orientatation.
!
   ELSEIF (TRAPT) THEN
      CALL NEWROTGEOM(NSIZE,RB,RMAT,0.0D0,0.0D0,0.0D0)
   ELSE
!
!  Translate the RB coordinates to the centre of coordinates of RA.
!
   DO J1=1,NATOMS
      RB(3*(J1-1)+1)=RB(3*(J1-1)+1)-CMXB+CMXA+XSHIFT
      RB(3*(J1-1)+2)=RB(3*(J1-1)+2)-CMYB+CMYA+YSHIFT
      RB(3*(J1-1)+3)=RB(3*(J1-1)+3)-CMZB+CMZA+ZSHIFT
   ENDDO
   IF (.NOT.BULKT) CALL NEWROTGEOM(NSIZE,RB,RMAT,CMXA,CMYA,CMZA)
   ENDIF
ENDIF

DEALLOCATE(XA,XB)

END SUBROUTINE NEWMINDIST

SUBROUTINE NEWROTGEOM(NATOMS,COORDS,ROTMAT,CX,CY,CZ)
IMPLICIT NONE
INTEGER I, J, K, NATOMS
DOUBLE PRECISION COORDS(*), R1, R0(3), ROTMAT(3,3), CX, CY, CZ

DO I=1,NATOMS
   R0(1)=COORDS(3*(I-1)+1)-CX
   R0(2)=COORDS(3*(I-1)+2)-CY
   R0(3)=COORDS(3*(I-1)+3)-CZ
   DO J=1,3
      R1=0.0D0
      DO K=1,3
         R1=R1+ROTMAT(J,K)*R0(K)
      ENDDO
      IF (J.EQ.1) COORDS(3*(I-1)+J)=R1+CX
      IF (J.EQ.2) COORDS(3*(I-1)+J)=R1+CY
      IF (J.EQ.3) COORDS(3*(I-1)+J)=R1+CZ
   ENDDO
ENDDO

RETURN
END SUBROUTINE NEWROTGEOM

