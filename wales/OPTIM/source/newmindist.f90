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
!  If PRESERVET is false we put RB into best correspondence with RA. This involves
!  a translation to the same centre of coordinates, followed by a rotation about that
!  centre.
!
SUBROUTINE NEWMINDIST(RA,RB,NATOMS,DIST,BULKT,TWOD,ZUSE,PRESERVET,RIGIDBODY,DEBUG,RMAT)
USE COMMONS,ONLY : PARAM1, PARAM2, PARAM3, NOPT
USE KEY,ONLY : STOCKT, NFREEZE, MIEFT, NOTRANSROTT, RBAAT, PULLT, EFIELDT, EYTRAPT, MKTRAPT, USERPOTT, &
  &            GTHOMSONT, VARIABLES, BULKBOXT! hk286

IMPLICIT NONE

INTEGER J1, NATOMS, NSIZE, INFO, JINFO, JMIN
DOUBLE PRECISION RA(NOPT), RB(NOPT), DIST, QMAT(4,4), XM, YM, ZM, XP, YP, ZP, OVEC(3), H1VEC(3), H2VEC(3), &
  &              DIAG(4), TEMPA(3*NOPT), RMAT(3,3), MINV, Q1, Q2, Q3, Q4, CMXA, CMYA, CMZA, CMXB, CMYB, CMZB, &
  &              MYROTMAT(3,3), OMEGATOT(3,3)
DOUBLE PRECISION, ALLOCATABLE :: XA(:), XB(:)
LOGICAL BULKT, TWOD, RIGIDBODY, PRESERVET, DEBUG, UPRETURN
CHARACTER(LEN=5) ZUSE
COMMON /MINDOM/ MYROTMAT, OMEGATOT
INTEGER NCIT
DOUBLE PRECISION XSHIFT, YSHIFT, ZSHIFT, XSHIFTNEW, YSHIFTNEW, ZSHIFTNEW, BOXLX, BOXLY, BOXLZ
DOUBLE PRECISION ENERGY, VNEW(3*NATOMS), RMS

! hk286
IF (GTHOMSONT) THEN
   CALL GTHOMSONNEWMINDIST(RA,RB,NATOMS,DIST,BULKT,TWOD,ZUSE,PRESERVET,RIGIDBODY,DEBUG,RMAT)
   RETURN
ENDIF

! DJW
IF (VARIABLES.OR.MKTRAPT) THEN
   DIST=0.0D0
   DO J1=1,NOPT
      DIST=DIST+(RA(J1)-RB(J1))**2
   ENDDO
   DIST=SQRT(DIST)
   RETURN
ENDIF

IF (RBAAT) THEN
   CALL RBMINDIST(RA,RB,NATOMS,DIST,RMAT,DEBUG)
   RETURN

ELSEIF (USERPOTT) THEN
   CALL USERPOT_DISTANCE(3*NATOMS, RA,RB,DIST)
   RMAT(1:3,1:3)=0.0D0 ! rotation matrix is the identity
   RMAT(1,1)=1.0D0; RMAT(2,2)=1.0D0; RMAT(3,3)=1.0D0
   UPRETURN=.TRUE.
   If (UPRETURN) RETURN
ENDIF

! CALL POTENTIAL(RA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
! PRINT '(2(A,F25.15))','newmindist> Initial RA energy=',ENERGY,' RMS=',RMS
! CALL POTENTIAL(RB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
! PRINT '(2(A,F25.15))','newmindist> Initial RB energy=',ENERGY,' RMS=',RMS

! WRITE(*,*) NATOMS
! WRITE(*,*) 'RA starting geometry'
! WRITE(*,'(A3,3G20.10)') ('LA ',RA(3*(J1-1)+1),RA(3*(J1-1)+2),RA(3*(J1-1)+3),J1=1,NATOMS)
! WRITE(*,*) NATOMS
! WRITE(*,*) 'RB starting geometry'
! WRITE(*,'(A3,3G20.10)') ('LA ',RB(3*(J1-1)+1),RB(3*(J1-1)+2),RB(3*(J1-1)+3),J1=1,NATOMS)
! 
! Convert rigid body coordinates to Cartesians for rigid bodies. 
!
IF (ZUSE(1:1).EQ.'W') THEN
   ALLOCATE(XA(3*3*(NATOMS/2)),XB(3*3*(NATOMS/2)))
   NSIZE=3*(NATOMS/2)
   DO J1=1,NATOMS/2
      CALL CONVERT(RA(3*(J1-1)+1),RA(3*(J1-1)+2),RA(3*(J1-1)+3), &
     &        RA(3*(NATOMS/2+J1-1)+1),RA(3*(NATOMS/2+J1-1)+2),RA(3*(NATOMS/2+J1-1)+3),OVEC,H1VEC,H2VEC)
      XA(9*(J1-1)+0+1)=OVEC(1)
      XA(9*(J1-1)+0+2)=OVEC(2)
      XA(9*(J1-1)+0+3)=OVEC(3)
      XA(9*(J1-1)+3+1)=H1VEC(1)
      XA(9*(J1-1)+3+2)=H1VEC(2)
      XA(9*(J1-1)+3+3)=H1VEC(3)
      XA(9*(J1-1)+6+1)=H2VEC(1)
      XA(9*(J1-1)+6+2)=H2VEC(2)
      XA(9*(J1-1)+6+3)=H2VEC(3)
      CALL CONVERT(RB(3*(J1-1)+1),RB(3*(J1-1)+2),RB(3*(J1-1)+3), &
     &      RB(3*(NATOMS/2+J1-1)+1),RB(3*(NATOMS/2+J1-1)+2),RB(3*(NATOMS/2+J1-1)+3),OVEC,H1VEC,H2VEC)
      XB(9*(J1-1)+0+1)=OVEC(1)
      XB(9*(J1-1)+0+2)=OVEC(2)
      XB(9*(J1-1)+0+3)=OVEC(3)
      XB(9*(J1-1)+3+1)=H1VEC(1)
      XB(9*(J1-1)+3+2)=H1VEC(2)
      XB(9*(J1-1)+3+3)=H1VEC(3)
      XB(9*(J1-1)+6+1)=H2VEC(1)
      XB(9*(J1-1)+6+2)=H2VEC(2)
      XB(9*(J1-1)+6+3)=H2VEC(3)
   ENDDO
ELSEIF (RIGIDBODY) THEN
   PRINT '(A)',' newmindist> New quaternion procedure not yet coded for general angle-axis variables'
   STOP
ELSEIF ((PULLT.OR.EFIELDT).AND.(.NOT.BULKT)) THEN
!  ALLOCATE(XA(3*(NATOMS/2)*number of sites,XB(3*(NATOMS/2)*number of sites))
!  NSIZE=(NATOMS/2)*number of sites
!  PRINT '(A)',' newmindist> New quaternion procedure not yet coded for flatland'
! There is one unknown angle, so this should be trivial!'
   CALL MINDIST(RA,RB,NATOMS,DIST,BULKT,TWOD,ZUSE,PRESERVET)
   RMAT(1:3,1:3)=OMEGATOT(1:3,1:3)
   RETURN
!  STOP
ELSEIF (STOCKT) THEN
   ALLOCATE(XA(3*(NATOMS/2)),XB(3*(NATOMS/2)))
   NSIZE=(NATOMS/2)
!  PRINT *,'newmindist> WARNING *** for STOCKT only CofM coordinates loaded into XA and XB???'
   XA(1:3*NSIZE)=RA(1:3*NSIZE)
   XB(1:3*NSIZE)=RB(1:3*NSIZE)
ELSE
   ALLOCATE(XA(3*NATOMS),XB(3*NATOMS))
   NSIZE=NATOMS
   XA(1:3*NATOMS)=RA(1:3*NATOMS)
   XB(1:3*NATOMS)=RB(1:3*NATOMS)
ENDIF
!
! If there are frozen atoms then just calculate the distance and return.
!
IF (NFREEZE.GT.0 .OR. MIEFT .OR. NOTRANSROTT) THEN
   DIST=0.0D0
   IF (BULKT) THEN
      BOXLX=PARAM1; BOXLY=PARAM2; BOXLZ=PARAM3
      DO J1=1,NSIZE
         DIST=DIST + (XA(3*(J1-1)+1)-XB(3*(J1-1)+1) - BOXLX*NINT((XA(3*(J1-1)+1)-XB(3*(J1-1)+1))/BOXLX))**2 &
   &               + (XA(3*(J1-1)+2)-XB(3*(J1-1)+2) - BOXLY*NINT((XA(3*(J1-1)+2)-XB(3*(J1-1)+2))/BOXLY))**2 &
   &               + (XA(3*(J1-1)+3)-XB(3*(J1-1)+3) - BOXLZ*NINT((XA(3*(J1-1)+3)-XB(3*(J1-1)+3))/BOXLZ))**2
      ENDDO
   ELSE
      ! Note: this is the site-by-site cartesian distance, which is equivalent to the angle-axis
      ! distance measure and is therefore valid for rigid bodies.
      DO J1=1,NSIZE
         DIST=DIST + (XA(3*(J1-1)+1)-XB(3*(J1-1)+1))**2 &
   &               + (XA(3*(J1-1)+2)-XB(3*(J1-1)+2))**2 &
   &               + (XA(3*(J1-1)+3)-XB(3*(J1-1)+3))**2
      ENDDO
   ENDIF
   DIST=SQRT(DIST)
   
   RMAT(1:3,1:3)=0.0D0 ! rotation matrix is the identity
   RMAT(1,1)=1.0D0; RMAT(2,2)=1.0D0; RMAT(3,3)=1.0D0
   RETURN
ENDIF
!
! If bulkboxt is on, the particles are allow to move out of the box, but
! structure still remains a box.
!
IF (BULKT.AND.(.NOT.BULKBOXT)) THEN
   BOXLX=PARAM1; BOXLY=PARAM2; BOXLZ=PARAM3
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
! ds656> Is this supposed to be done with substrate?!?
!
CMXA=0.0D0; CMYA=0.0D0; CMZA=0.0D0
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
!ENDIF
! PRINT '(A,6F15.7)','CMA,CMB=',CMXA,CMYA,CMZA,CMXB,CMYB,CMZB

!  
! Deal with TWOD non-bulk here.
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

XSHIFT=0.0D0; YSHIFT=0.0D0; ZSHIFT=0.0D0
NCIT=0
IF (BULKT) THEN 
   BOXLX=PARAM1; BOXLY=PARAM2; BOXLZ=PARAM3
! Iterative solution
! 1  NCIT=NCIT+1
!    IF (NCIT.GT.1000) THEN
!       PRINT '(A)','inertia> WARNING - iterative calculation of centre of mass shift did not converge'
!    ENDIF
!    XSHIFTNEW=0.0D0
!    YSHIFTNEW=0.0D0
!    ZSHIFTNEW=0.0D0
!    DO J1=1,NSIZE
!       XSHIFTNEW= SHIFTNEW + XA(3*(J1-1)+1)-XB(3*(J1-1)+1) - BOXLX*NINT((XA(3*(J1-1)+1)-XB(3*(J1-1)+1)-XSHIFT)/BOXLX)
!       YSHIFTNEW=YSHIFTNEW + XA(3*(J1-1)+2)-XB(3*(J1-1)+2) - BOXLY*NINT((XA(3*(J1-1)+2)-XB(3*(J1-1)+2)-YSHIFT)/BOXLY)
!       IF (.NOT.TWOD) ZSHIFTNEW=ZSHIFTNEW + XA(3*(J1-1)+3)-XB(3*(J1-1)+3) - BOXLZ*NINT((XA(3*(J1-1)+3)-XB(3*(J1-1)+3)-ZSHIFT)/BOXLZ)
!    ENDDO
!    XSHIFTNEW=XSHIFTNEW/NSIZE; YSHIFTNEW=YSHIFTNEW/NSIZE; ZSHIFTNEW=ZSHIFTNEW/NSIZE
!    IF ((ABS(XSHIFTNEW-XSHIFT).GT.1.0D-6).OR.(ABS(YSHIFTNEW-YSHIFT).GT.1.0D-6).OR.(ABS(ZSHIFTNEW-ZSHIFT).GT.1.0D-6)) THEN
! !     IF (DEBUG) PRINT '(A,I6,6F15.7)',' newmindist> ',NCIT,XSHIFTNEW,YSHIFTNEW,ZSHIFTNEW,XSHIFT,YSHIFT,ZSHIFT
!       XSHIFT=0.05D0*XSHIFT+0.95D0*XSHIFTNEW
!       YSHIFT=0.05D0*YSHIFT+0.95D0*YSHIFTNEW
!       IF (.NOT.TWOD) ZSHIFT=0.05D0*ZSHIFT+0.95D0*ZSHIFTNEW
!       GOTO 1
!    ENDIF
! !  IF (DEBUG) PRINT '(A,I6,3F15.7)',' newmindist> coordinate shift converged. Cycles and values: ',NCIT,XSHIFT,YSHIFT,ZSHIFT
!    XSHIFT=XSHIFTNEW
!    YSHIFT=YSHIFTNEW 
!    ZSHIFT=ZSHIFTNEW
!
! Actually, the iterative solution seems to be worse than simply putting the centre of mass
! at the origin. 
!
   DIST=0.0D0
   DO J1=1,NSIZE
    DIST=DIST + (XA(3*(J1-1)+1)-XB(3*(J1-1)+1)-XSHIFT - BOXLX*NINT((XA(3*(J1-1)+1)-XB(3*(J1-1)+1)-XSHIFT)/BOXLX))**2 &
   &          + (XA(3*(J1-1)+2)-XB(3*(J1-1)+2)-YSHIFT - BOXLY*NINT((XA(3*(J1-1)+2)-XB(3*(J1-1)+2)-YSHIFT)/BOXLY))**2 
      IF (.NOT.TWOD) DIST=DIST &
   &          + (XA(3*(J1-1)+3)-XB(3*(J1-1)+3)-ZSHIFT - BOXLZ*NINT((XA(3*(J1-1)+3)-XB(3*(J1-1)+3)-ZSHIFT)/BOXLZ))**2
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
!     PRINT '(A,I8,6G18.8)','J1,XM,YM,ZM,XP,YP,ZP=',J1,XM,YM,ZM,XP,YP,ZP
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
!  PRINT '(A,G20.10)','QMAT(1,1)=',QMAT(1,1)
!  PRINT '(A,G20.10)','QMAT(1,2)=',QMAT(1,2)
!  PRINT '(A,G20.10)','QMAT(1,3)=',QMAT(1,3)
!  PRINT '(A,G20.10)','QMAT(1,4)=',QMAT(1,4)
!  PRINT '(A,G20.10)','QMAT(2,2)=',QMAT(2,2)
!  PRINT '(A,G20.10)','QMAT(2,3)=',QMAT(2,3)
!  PRINT '(A,G20.10)','QMAT(2,4)=',QMAT(2,4)
!  PRINT '(A,G20.10)','QMAT(3,3)=',QMAT(3,3)
!  PRINT '(A,G20.10)','QMAT(3,4)=',QMAT(3,4)
!  PRINT '(A,G20.10)','QMAT(4,4)=',QMAT(4,4)

   CALL DSYEV('V','U',4,QMAT,4,DIAG,TEMPA,9*NATOMS,INFO)
   IF (INFO.NE.0) PRINT '(A,I6,A)',' newmindist> WARNING - INFO=',INFO,' in DSYEV'

   MINV=1.0D200
   JMIN=1
   DO J1=1,4
!     PRINT '(A,I8,G20.10)','newmindist> J1,DIAG=',J1,DIAG(J1)
      IF (DIAG(J1).LT.MINV) THEN
         JMIN=J1
         MINV=DIAG(J1)
      ENDIF
   ENDDO
   IF (MINV.LT.0.0D0) THEN
      IF (ABS(MINV).LT.1.0D-6) THEN
         MINV=0.0D0
      ELSE
         PRINT '(A,G20.10,A)',' newmindist> WARNING MINV is ',MINV,' change to absolute value'
         MINV=-MINV
      ENDIF
   ENDIF
   DIST=SQRT(MINV)

!  IF (DEBUG) PRINT '(A,G20.10,A,I6)',' newmindist> minimum residual is ',DIAG(JMIN),' for eigenvector ',JMIN
   Q1=QMAT(1,JMIN); Q2=QMAT(2,JMIN); Q3=QMAT(3,JMIN); Q4=QMAT(4,JMIN)
!
! RMAT will contain the matrix that maps XB onto the best correspondence with XA
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
IF (.NOT.(PRESERVET.OR.MIEFT.OR.NOTRANSROTT)) THEN ! never translate with substrate!
   IF (ZUSE(1:1).EQ.'W') THEN
!
!  Translate the XB coordinates to the centre of coordinates of XA.
!
      DO J1=1,NSIZE
         XB(3*(J1-1)+1)=XB(3*(J1-1)+1)+CMXA+XSHIFT
         XB(3*(J1-1)+2)=XB(3*(J1-1)+2)+CMYA+YSHIFT
         XB(3*(J1-1)+3)=XB(3*(J1-1)+3)+CMZA+ZSHIFT
      ENDDO
!
!  Rotate XB coordinates about new centre of mass
!
      CALL NEWROTGEOM(NSIZE,XB,RMAT,CMXA,CMYA,CMZA)
      DO J1=1,NATOMS/2
         OVEC(1)=XB(1+(J1-1)*9+0)
         OVEC(2)=XB(2+(J1-1)*9+0)
         OVEC(3)=XB(3+(J1-1)*9+0)
         H1VEC(1)=XB(1+(J1-1)*9+3)
         H1VEC(2)=XB(2+(J1-1)*9+3)
         H1VEC(3)=XB(3+(J1-1)*9+3)
         H2VEC(1)=XB(1+(J1-1)*9+6)
         H2VEC(2)=XB(2+(J1-1)*9+6)
         H2VEC(3)=XB(3+(J1-1)*9+6)
         CALL CONVERT2(OVEC,H1VEC,H2VEC,RB(3*(J1-1)+1),RB(3*(J1-1)+2),RB(3*(J1-1)+3), &
  &                    RB(3*(NATOMS/2+J1-1)+1),RB(3*(NATOMS/2+J1-1)+2),RB(3*(NATOMS/2+J1-1)+3))
      ENDDO
   ELSEIF (RIGIDBODY) THEN
!
!  Needs some thought for the angle/axis rigid body formulation.
!
      PRINT '(A)',' newmindist> WARNING *** back transformation not programmed yet for rigid bodies'
!
!  Preserve centre of mass for EYTRAP, but not orientatation.
!
!  ELSEIF (EYTRAPT.OR.MKTRAPT) THEN
   ELSEIF (EYTRAPT) THEN
      CALL NEWROTGEOM(NSIZE,RB,RMAT,0.0D0,0.0D0,0.0D0)
   !ELSEIF (BULKT) THEN
   !   DO J1=1,NSIZE
   !      RB(3*(J1-1)+1)=RB(3*(J1-1)+1)+VPBX+ZSHIFT
   !      RB(3*(J1-1)+2)=RB(3*(J1-1)+2)+VPBY+YSHIFT
   !      RB(3*(J1-1)+3)=RB(3*(J1-1)+3)+VPBZ+ZSHIFT
   !      RB(3*(J1-1)+1)=RB(3*(J1-1)+1)+BOXLX*NINT(RB(3*(J1-1)+1)/BOXLX)
   !      RB(3*(J1-1)+2)=RB(3*(J1-1)+2)+BOXLY*NINT(RB(3*(J1-1)+2)/BOXLY)
   !      RB(3*(J1-1)+3)=RB(3*(J1-1)+3)+BOXLZ*NINT(RB(3*(J1-1)+3)/BOXLZ)
   !   ENDDO
   ELSE
!
!  Translate the RB coordinates to the centre of coordinates of RA.
!
      DO J1=1,NSIZE
         RB(3*(J1-1)+1)=RB(3*(J1-1)+1)-CMXB+CMXA+XSHIFT
         RB(3*(J1-1)+2)=RB(3*(J1-1)+2)-CMYB+CMYA+YSHIFT
         RB(3*(J1-1)+3)=RB(3*(J1-1)+3)-CMZB+CMZA+ZSHIFT
      ENDDO

      !IF (BULKT) THEN
      !   BOXLX=PARAM1; BOXLY=PARAM2; BOXLZ=PARAM3
      !   DO J1=1,NSIZE
      !      RB(3*(J1-1)+1)=RB(3*(J1-1)+1)-BOXLX*NINT(RB(3*(J1-1)+1)/BOXLX)
      !      RB(3*(J1-1)+2)=RB(3*(J1-1)+2)-BOXLY*NINT(RB(3*(J1-1)+2)/BOXLY)
      !      IF (.NOT.TWOD) RB(3*(J1-1)+3)=RB(3*(J1-1)+3)-BOXLZ*NINT(RB(3*(J1-1)+3)/BOXLZ)
      !   ENDDO
      !ENDIF


      IF (.NOT.BULKT) THEN
         IF (STOCKT) THEN
            CALL NEWROTGEOMSTOCK(NATOMS,RB,RMAT,CMXA,CMYA,CMZA)
         ELSE
            CALL NEWROTGEOM(NSIZE,RB,RMAT,CMXA,CMYA,CMZA)
         ENDIF
      ENDIF
   ENDIF
ENDIF

DEALLOCATE(XA,XB)

! WRITE(*,*) NATOMS
! WRITE(*,*) 'RB finishing geometry'
! WRITE(*,'(A3,3G20.10)') ('LA ',RB(3*(J1-1)+1),RB(3*(J1-1)+2),RB(3*(J1-1)+3),J1=1,NATOMS)

RETURN

END SUBROUTINE NEWMINDIST

SUBROUTINE NEWROTGEOM(NATOMS,COORDS,MYROTMAT,CX,CY,CZ)
IMPLICIT NONE
INTEGER I, J, K, NATOMS
DOUBLE PRECISION COORDS(*), R1, R0(3), MYROTMAT(3,3), CX, CY, CZ

DO I=1,NATOMS
   R0(1)=COORDS(3*(I-1)+1)-CX
   R0(2)=COORDS(3*(I-1)+2)-CY
   R0(3)=COORDS(3*(I-1)+3)-CZ
   DO J=1,3
      R1=0.0D0
      DO K=1,3
         R1=R1+MYROTMAT(J,K)*R0(K)
      ENDDO
      IF (J.EQ.1) COORDS(3*(I-1)+J)=R1+CX
      IF (J.EQ.2) COORDS(3*(I-1)+J)=R1+CY
      IF (J.EQ.3) COORDS(3*(I-1)+J)=R1+CZ
   ENDDO
ENDDO

RETURN
END SUBROUTINE NEWROTGEOM

SUBROUTINE OLDROTGEOMSTOCK(NATOMS,COORDS,MYROTMAT,CX,CY,CZ)
IMPLICIT NONE
INTEGER I, J, K, NATOMS, NREALATOMS, J3, J1, OFFSET
DOUBLE PRECISION COORDS(*), R1, R0(3), MYROTMAT(3,3), CX, CY, CZ, X1, Y1, Z1, X2, Y2, Z2, CT1, ST1, P1, CP1, SP1, T1B, P1B, T1
DOUBLE PRECISION START(3), FINISH(3), DIFF, DIFFBEST, DUMMY
DOUBLE PRECISION, PARAMETER ::  PI=3.141592654D0

NREALATOMS=(NATOMS/2)
OFFSET = 3*NREALATOMS
!
! First rotate the dipoles.
!
DO J1=1,NREALATOMS
   J3=3*J1
   X1=COORDS(J3-2)-CX
   Y1=COORDS(J3-1)-CY
   Z1=COORDS(J3)-CZ
   T1=COORDS(OFFSET+J3-2)
   CT1=COS(T1)
   ST1=SIN(T1)
   P1=COORDS(OFFSET+J3-1)
   CP1=COS(P1)
   SP1=SIN(P1)
   X2=X1+ST1*CP1
   Y2=Y1+ST1*SP1
   Z2=Z1+CT1
   START(1)=MYROTMAT(1,1)*X1+MYROTMAT(1,2)*Y1+MYROTMAT(1,3)*Z1
   START(2)=MYROTMAT(2,1)*X1+MYROTMAT(2,2)*Y1+MYROTMAT(2,3)*Z1
   START(3)=MYROTMAT(3,1)*X1+MYROTMAT(3,2)*Y1+MYROTMAT(3,3)*Z1
   FINISH(1)=MYROTMAT(1,1)*X2+MYROTMAT(1,2)*Y2+MYROTMAT(1,3)*Z2
   FINISH(2)=MYROTMAT(2,1)*X2+MYROTMAT(2,2)*Y2+MYROTMAT(2,3)*Z2
   FINISH(3)=MYROTMAT(3,1)*X2+MYROTMAT(3,2)*Y2+MYROTMAT(3,3)*Z2
   DUMMY=FINISH(3)-START(3)
   IF (DUMMY.GT.1.0D0) DUMMY=1.0D0; IF (DUMMY.LT.-1.0D0) DUMMY=-1.0D0
   T1=ACOS(DUMMY)
   DIFFBEST=1.0D100
   IF (SIN(T1).NE.0.0D0) THEN
      DUMMY=(FINISH(1)-START(1))/SIN(T1)
      IF (DUMMY.GT.1.0D0) DUMMY=1.0D0; IF (DUMMY.LT.-1.0D0) DUMMY=-1.0D0
      P1=ACOS(DUMMY)
   ELSE
      P1=1.0D0
   ENDIF
   DIFFBEST=(FINISH(1)-START(1)-SIN(T1)*COS(P1))**2+(FINISH(2)-START(2)-SIN(T1)*SIN(P1))**2+(FINISH(3)-START(3)-COS(T1))**2
   T1B=T1; P1B=P1
   DIFF=(FINISH(1)-START(1)-SIN(2*PI-T1)*COS(P1))**2+(FINISH(2)-START(2)-SIN(2*PI-T1)*SIN(P1))**2+ &
  &          (FINISH(3)-START(3)-COS(2*PI-T1))**2
   IF (DIFF.LT.DIFFBEST) THEN
      T1B=2*PI-T1; P1B=P1
      DIFFBEST=DIFF
   ENDIF
   DIFF=(FINISH(1)-START(1)-SIN(2*PI-T1)*COS(2*PI-P1))**2+(FINISH(2)-START(2)-SIN(2*PI-T1)*SIN(2*PI-P1))**2+ &
  &           (FINISH(3)-START(3)-COS(2*PI-T1))**2
   IF (DIFF.LT.DIFFBEST) THEN
      T1B=2*PI-T1; P1B=2*PI-P1
      DIFFBEST=DIFF
   ENDIF
   DIFF=(FINISH(1)-START(1)-SIN(T1)*COS(2*PI-P1))**2+(FINISH(2)-START(2)-SIN(T1)*SIN(2*PI-P1))**2+(FINISH(3)-START(3)-COS(T1))**2
   IF (DIFF.LT.DIFFBEST) THEN
      T1B=T1; P1B=2*PI-P1
      DIFFBEST=DIFF
   ENDIF
   IF (DIFFBEST.GT.1.0D-5) THEN
      PRINT '(A,G20.10)','newrotgeomstock> WARNING - angle rotation failed - DIFFBEST=',DIFFBEST
   ENDIF
!  PRINT '(A,G20.10)','newrotgeomstock> DIFFBEST=',DIFFBEST
!
!  Inverse cos gives us an angle between 0 and pi. However, 2*pi - angle gives
!  the same cos. There are therefore two possibilities for theta and two for phi,
!  and only one should regenerate the correct displacements. Find it!
!
   COORDS(OFFSET+J3-2)=T1B; COORDS(OFFSET+J3-1)=P1B
ENDDO

DO I=1,(NATOMS/2)
   R0(1)=COORDS(3*(I-1)+1)-CX
   R0(2)=COORDS(3*(I-1)+2)-CY
   R0(3)=COORDS(3*(I-1)+3)-CZ
   DO J=1,3
      R1=0.0D0
      DO K=1,3
         R1=R1+MYROTMAT(J,K)*R0(K)
      ENDDO
      IF (J.EQ.1) COORDS(3*(I-1)+J)=R1+CX
      IF (J.EQ.2) COORDS(3*(I-1)+J)=R1+CY
      IF (J.EQ.3) COORDS(3*(I-1)+J)=R1+CZ
   ENDDO
ENDDO

RETURN
END SUBROUTINE OLDROTGEOMSTOCK

SUBROUTINE NEWROTGEOMSTOCK(NATOMS,COORDS,MYROTMAT,CX,CY,CZ)
IMPLICIT NONE
INTEGER I, J, K, NATOMS, NREALATOMS, J3, J1, OFFSET
DOUBLE PRECISION COORDS(*), BEFORE(3), AFTER(3), MYROTMAT(3,3)
DOUBLE PRECISION CX, CY, CZ, THETA, PHI
DOUBLE PRECISION START(3), FINISH(3), DIFF, DIFFBEST, DUMMY
DOUBLE PRECISION, PARAMETER ::  PI=3.141592654D0

NREALATOMS=(NATOMS/2)
OFFSET = 3*NREALATOMS
DO J1=1,NREALATOMS
   J3 = J1*3
   THETA = COORDS(OFFSET+J3-2)
   PHI = COORDS(OFFSET+J3-1)
!  Make a unit vector pointing along the dipole.
   BEFORE(1) = SIN(THETA) * COS(PHI)
   BEFORE(2) = SIN(THETA) * SIN(PHI)
   BEFORE(3) = COS(THETA)
!  Rotate the unit vector using the rotation matrix.
   AFTER(1) = MYROTMAT(1,1)*BEFORE(1) + MYROTMAT(1,2)*BEFORE(2) + MYROTMAT(1,3)*BEFORE(3)
   AFTER(2) = MYROTMAT(2,1)*BEFORE(1) + MYROTMAT(2,2)*BEFORE(2) + MYROTMAT(2,3)*BEFORE(3)
   AFTER(3) = MYROTMAT(3,1)*BEFORE(1) + MYROTMAT(3,2)*BEFORE(2) + MYROTMAT(3,3)*BEFORE(3)
!  Convert the unit vector back to spherical polars.
   IF (AFTER(3) > 1.0D0) AFTER(3)=1.0D0
   IF (AFTER(3) < -1.0D0) AFTER(3)=-1.0D0
   COORDS(OFFSET+J3-2) = ACOS(AFTER(3))
   COORDS(OFFSET+J3-1) = ATAN2(AFTER(2), AFTER(1))
ENDDO
!
! Now rotate the particle positions.
! 
DO I=1,(NATOMS/2)
   BEFORE(1)=COORDS(3*(I-1)+1)-CX
   BEFORE(2)=COORDS(3*(I-1)+2)-CY
   BEFORE(3)=COORDS(3*(I-1)+3)-CZ
   AFTER(1) = MYROTMAT(1,1)*BEFORE(1) + MYROTMAT(1,2)*BEFORE(2) + MYROTMAT(1,3)*BEFORE(3)
   AFTER(2) = MYROTMAT(2,1)*BEFORE(1) + MYROTMAT(2,2)*BEFORE(2) + MYROTMAT(2,3)*BEFORE(3)
   AFTER(3) = MYROTMAT(3,1)*BEFORE(1) + MYROTMAT(3,2)*BEFORE(2) + MYROTMAT(3,3)*BEFORE(3)
   COORDS(3*(I-1)+1) = AFTER(1) + CX
   COORDS(3*(I-1)+2) = AFTER(2) + CY
   COORDS(3*(I-1)+3) = AFTER(3) + CZ
ENDDO

RETURN
END SUBROUTINE NEWROTGEOMSTOCK
