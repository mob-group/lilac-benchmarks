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
C jmc May need to use difference between polar angles as well in the order parameter for unres
C like with twisttype10...
C
C Routine to calculate difference between different internal coordinate dihedral values for QLOCAL and FIN (Cartesian) coords
C
      SUBROUTINE UNRESCALCDIHE(DIHE,ALLANG,QLOCAL,FIN)
      USE COMMONS
c     USE KEY
      USE MODUNRES
      IMPLICIT NONE
C
      INTEGER I1,J1
      REAL*8 DIFFPP,DIHE, SUMD2,ALLDIFFPP, ALLANG, ALLSUMD2,FINPPSANGLE(4*nres-9) ! dim of *ppsangle is nphi+ntheta+2(nres-2)
      REAL*8 QPPSANGLE(4*nres-9),DIFFARRAY(4*nres-9),QLOCAL(3*NATOMS),FIN(3*NATOMS)
      INTEGER IICD,TWISTMODE,TWISTTYPE,NM,NWRONG
      LOGICAL LINTCOOR,PTEST
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)

      DO I1=1,nres
         c(1,I1)=FIN(6*(I1-1)+1)
         c(2,I1)=FIN(6*(I1-1)+2)
         c(3,I1)=FIN(6*(I1-1)+3)
         c(1,I1+nres)=FIN(6*(I1-1)+4)
         c(2,I1+nres)=FIN(6*(I1-1)+5)
         c(3,I1+nres)=FIN(6*(I1-1)+6)
      END DO
      CALL UPDATEDC
      CALL int_from_cart(.true.,.false.)

C use unres geometry arrays phi (bb dihedrals) and omeg (sc dihedrals)
C Take care with numbering - see /unres/src/readpdb.f (subroutine int_from_cart)
C For side chain dihedrals, have zero elements for proper (i.e. not capping) glycines.
C Need to remember not to try to twist them though!
C No entries in QPPSANGLE for capping 'residues'.
      DO I1=1,nres-3
        FINPPSANGLE(I1)=phi(I1+3)
      END DO
      DO I1=1,nres-2
        FINPPSANGLE(I1+nres-3)=omeg(I1+1)
C jmc 30/4/03 try adding backbone and side chain polar angles
C This should be more important for unres than for charmm...
C Order is bb dihedrals, sc dihedrals, bb polars, sc polars.
        FINPPSANGLE(I1+2*nres-5)=theta(I1+2)
        FINPPSANGLE(I1+3*nres-7)=alph(I1+1)
      END DO

      DO I1=1,nres
         c(1,I1)=QLOCAL(6*(I1-1)+1)
         c(2,I1)=QLOCAL(6*(I1-1)+2)
         c(3,I1)=QLOCAL(6*(I1-1)+3)
         c(1,I1+nres)=QLOCAL(6*(I1-1)+4)
         c(2,I1+nres)=QLOCAL(6*(I1-1)+5)
         c(3,I1+nres)=QLOCAL(6*(I1-1)+6)
      END DO
      CALL UPDATEDC
      CALL int_from_cart(.true.,.false.)

      DO I1=1,nres-3
        QPPSANGLE(I1)=phi(I1+3)
      END DO
      DO I1=1,nres-2
        QPPSANGLE(I1+nres-3)=omeg(I1+1)
C jmc 30/4/03 try adding backbone and side chain polar angles
C This should be more important for unres than for charmm...
C Order is bb dihedrals, sc dihedrals, bb polars, sc polars.
        QPPSANGLE(I1+2*nres-5)=theta(I1+2)
        QPPSANGLE(I1+3*nres-7)=alph(I1+1)
      END DO
C
      SUMD2=0.D0
      DO I1=1,NPHI+NRES-2
         DIFFPP = QPPSANGLE(I1) - FINPPSANGLE(I1)
C
C next two lines are meant to ensure that you always interpolate
C along the shortest distance between the dihedral angles.
C
C jmc         IF (DIFFPP.GT.180.0) DIFFPP = DIFFPP-360.D0
C jmc         IF (DIFFPP.GT.180.0) DIFFPP = DIFFPP-360.D0
         IF (DIFFPP.LT.-PI) DIFFPP = DIFFPP+2.0D0*PI
         IF (DIFFPP.GT.PI) DIFFPP = DIFFPP-2.0D0*PI
c        WRITE(*,'(A,I6,3F15.10)') 'FINS QS',I1,FINPPSANGLE(I1),QPPSANGLE(I1),DIFFPP
         SUMD2=SUMD2+DIFFPP*DIFFPP
      ENDDO

      ALLSUMD2=SUMD2
      DO I1=NPHI+NRES-1,NPHI+ntheta+2*nres-4
         DIFFPP = QPPSANGLE(I1) - FINPPSANGLE(I1)
c        WRITE(*,'(A,I6,3F15.10)') 'FINS QS',I1,FINPPSANGLE(I1),QPPSANGLE(I1),DIFFPP
         ALLSUMD2=ALLSUMD2+DIFFPP*DIFFPP
      ENDDO
 
C jmc testing
c     PRINT *,'DIHE,ALLANG ',SUMD2,ALLSUMD2
C jmc REMEMBER if sumd is in radians, then dihe will have different
C range of values for charmm vs unres...
      DIHE=DSQRT(SUMD2/(NPHI+nside))
      ALLANG=DSQRT(ALLSUMD2/(NPHI+ntheta+2.0D0*nside))

      RETURN

      END

C jmc May need to use difference between polar angles as well in the order parameter for unres
C like with twisttype10...
C
C Routine to calculate difference between
C different internal coordinate dihedral values for passed coordinate array QLOCAL and a stored reference array
C (see keywords.f)
C
      SUBROUTINE UNRESCALCDIHEREF(DIHE,ALLANG,QLOCAL)
      USE COMMONS
c     USE KEY
      USE MODUNRES
      IMPLICIT NONE

      INTEGER I1,J1
C
      REAL*8 X(NATOMS),Y(NATOMS),Z(NATOMS)
      REAL*8 DIFFPP,DIHE, SUMD2,ALLDIFFPP, ALLANG, ALLSUMD2
      REAL*8 QPPSANGLE(4*nres-9),DIFFARRAY(4*nres-9),QLOCAL(3*NATOMS)
      INTEGER IICD,TWISTMODE,TWISTTYPE,NM,NWRONG
      LOGICAL LINTCOOR,PTEST
      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)

      DO I1=1,nres
         c(1,I1)=QLOCAL(6*(I1-1)+1)
         c(2,I1)=QLOCAL(6*(I1-1)+2)
         c(3,I1)=QLOCAL(6*(I1-1)+3)
         c(1,I1+nres)=QLOCAL(6*(I1-1)+4)
         c(2,I1+nres)=QLOCAL(6*(I1-1)+5)
         c(3,I1+nres)=QLOCAL(6*(I1-1)+6)
c     PRINT *,'QLOCAL in unrescalcdihe: ',QLOCAL(6*(I1-1)+1),QLOCAL(6*(I1-1)+2),QLOCAL(6*(I1-1)+3)
c     PRINT *,'QLOCAL in unrescalcdihe: ',QLOCAL(6*(I1-1)+4),QLOCAL(6*(I1-1)+5),QLOCAL(6*(I1-1)+6)
      END DO
      CALL UPDATEDC
      CALL int_from_cart(.true.,.false.)

C use unres geometry arrays phi (bb dihedrals) and omeg (sc dihedrals)
C Take care with numbering - see /unres/src/readpdb.f (subroutine int_from_cart)
C For side chain dihedrals, have zero elements for proper (i.e. not capping) glycines.
C Need to remember not to try to twist them though!
C No entries in QPPSANGLE for capping 'residues'.
      DO I1=1,nres-3
        QPPSANGLE(I1)=phi(I1+3)
      END DO
      DO I1=1,nres-2
        QPPSANGLE(I1+nres-3)=omeg(I1+1)
C jmc 30/4/03 try adding backbone and side chain polar angles
C This should be more important for unres than for charmm...
C Order is bb dihedrals, sc dihedrals, bb polars, sc polars.
        QPPSANGLE(I1+2*nres-5)=theta(I1+2)
        QPPSANGLE(I1+3*nres-7)=alph(I1+1)
      END DO
C
      SUMD2=0.D0
      DO I1=1,NPHI+NRES-2
         DIFFPP = QPPSANGLE(I1) - UREFPPSANGLE(I1)

C
C next two lines are meant to ensure that you always interpolate
C along the shortest distance between the dihedral angles.
C
C jmc         IF (DIFFPP.GT.180.0) DIFFPP = DIFFPP-360.D0
C jmc         IF (DIFFPP.GT.180.0) DIFFPP = DIFFPP-360.D0
         IF (DIFFPP.LT.-PI) DIFFPP = DIFFPP+2.0D0*PI
         IF (DIFFPP.GT.PI) DIFFPP = DIFFPP-2.0D0*PI
         SUMD2=SUMD2+DIFFPP*DIFFPP
      ENDDO

      ALLSUMD2=SUMD2
      DO I1=NPHI+NRES-1,NPHI+ntheta+2*nres-4
         DIFFPP = QPPSANGLE(I1) - UREFPPSANGLE(I1)
C
C next two lines are meant to ensure that you always interpolate
C along the shortest distance between the angles.
C
C jmc         IF (DIFFPP.GT.180.0) DIFFPP = DIFFPP-360.D0
C jmc         IF (DIFFPP.GT.180.0) DIFFPP = DIFFPP-360.D0
c        IF (DIFFPP.LT.-PI) DIFFPP = DIFFPP+2.0D0*PI
c        IF (DIFFPP.GT.PI) DIFFPP = DIFFPP-2.0D0*PI
         ALLSUMD2=ALLSUMD2+DIFFPP*DIFFPP
      ENDDO
 
C jmc testing
c     PRINT *,'DIHE,ALLANG ',SUMD2,ALLSUMD2
C jmc REMEMBER if sumd is in radians, then dihe will have different
C range of values for charmm vs unres...
      DIHE=DSQRT(SUMD2/(NPHI+nside))
      ALLANG=DSQRT(ALLSUMD2/(NPHI+ntheta+2.0D0*nside))

      RETURN

      END

C
C routine to read in reference coordinates in plain xyz format (unres) for mind/rmsd comparison
C stored in CHREF common block
C
      SUBROUTINE UNREADREF(NATOMS)
      USE MODUNRES
      IMPLICIT NONE
       
      INTEGER I1,NATOMS,J1
      REAL*8 X(NATOMS), Y(NATOMS), Z(NATOMS)
      CHARACTER*2 DUM

      OPEN (UNIT=10,FILE='ref.crd',STATUS='OLD')

      READ(10,*)
      DO I1=1,NATOMS
         READ(10,*) DUM,UREFCOORD(3*(I1-1)+1),UREFCOORD(3*(I1-1)+2),UREFCOORD(3*(I1-1)+3)
      ENDDO
      CLOSE(10)
C
      DO I1=1,nres
         c(1,I1)=UREFCOORD(6*(I1-1)+1)
         c(2,I1)=UREFCOORD(6*(I1-1)+2)
         c(3,I1)=UREFCOORD(6*(I1-1)+3)
         c(1,I1+nres)=UREFCOORD(6*(I1-1)+4)
         c(2,I1+nres)=UREFCOORD(6*(I1-1)+5)
         c(3,I1+nres)=UREFCOORD(6*(I1-1)+6)
      END DO
      CALL UPDATEDC
      CALL int_from_cart(.true.,.false.)

C use unres geometry arrays phi (bb dihedrals) and omeg (sc dihedrals)
C Take care with numbering - see /unres/src/readpdb.f (subroutine int_from_cart)
C For side chain dihedrals, have zero elements for proper (i.e. not capping) glycines.
C Need to remember not to try to twist them though!
C No entries in QPPSANGLE for capping 'residues'.
      DO I1=1,nres-3
        UREFPPSANGLE(I1)=phi(I1+3)
      END DO
      DO I1=1,nres-2
        UREFPPSANGLE(I1+nres-3)=omeg(I1+1)
C jmc 30/4/03 try adding backbone and side chain polar angles
C This should be more important for unres than for charmm...
C Order is bb dihedrals, sc dihedrals, bb polars, sc polars.
        UREFPPSANGLE(I1+2*nres-5)=theta(I1+2)
        UREFPPSANGLE(I1+3*nres-7)=alph(I1+1)
      END DO

      RETURN
      END

C
C Routine to calculate radius of gyration
C
      SUBROUTINE UNRESCALCRGYR(RGYR,QLOCAL)
      USE COMMONS
      IMPLICIT NONE

      INTEGER I1
      REAL*8 X(NATOMS),Y(NATOMS),Z(NATOMS),QLOCAL(3*NATOMS)
      REAL*8 RGYR
      REAL*8 AM(NATOMS),W(NATOMS)
      INTEGER ISLCT(NATOMS)
      REAL*8 FACT
      LOGICAL LMASS,LWEIG
 
      DO I1=1,NATOMS
        X(I1)=QLOCAL(3*(I1-1)+1)
        Y(I1)=QLOCAL(3*(I1-1)+2)
        Z(I1)=QLOCAL(3*(I1-1)+3)
        ISLCT(I1)=1
      ENDDO

      LMASS=.FALSE.
      LWEIG=.FALSE.
      FACT=0.0D0

C jmc W is weighting array.
C jmc we're always calculating the geometric rgy, since lmass and lweig are false.
      CALL UEDITRGYR(RGYR,NATOMS,X,Y,Z,W,AM,ISLCT,FACT,LMASS,LWEIG)

      RETURN
      END
C
C DAE modified solely to pass RG back to OPTIM
C

CHARMM Element source/manip/rgyr.src 1.1

      SUBROUTINE UEDITRGYR(RG,NATOMS,X,Y,Z,W,AM,ISLCT,FACT,LMASS,LWEIG)

      IMPLICIT NONE
C-----------------------------------------------------------------------
C     Compute the radius of gyration, center of mass,
C     and total mass of the selected subset of either the main or
C     the comparison structure. The results are output to unit OUTU.
C     If keyword WEIG is given the weighting array, which the user
C     must fill with CORMAN or SCALAR commands, is used for the
C     weighting. This is indicated by LWEIG=.TRUE. and is taken care
C     of in CORMAN's call to RGYR.
C     If keyword MASS is given, then the mass-weighted radius of
C     gyration, etc are computed, otherwise unit weighting per
C     atom is used (giving rms distance from the geometric center).
C     The default option is to do e geometric RGYR calculation.
C     A constant offset to be subtracted from the weights can be
C     specified with keyword FACT.
C
C     SYNTAX:
C
C     COOR RGYR  [FACT <real>] {[MASS]} [COMP]  [<atom-selection>]
C     {[WEIG]}
C
C     1983-09-01/LN
C     1985-01-05/ Default revised /LN
C
C##INCLUDE '/export/home/dae22/charmmcode/fcm/impnon.fcm'
C##INCLUDE '/export/home/dae22/charmmcode/fcm/number.fcm'
C##INCLUDE '/export/home/dae22/charmmcode/fcm/stream.fcm'
C     IMPLICIT REAL*8(A-H,O-Z)

      INTEGER NATOMS
      REAL*8 X(*),Y(*),Z(*),W(*)
      REAL*8 AM(*)
      INTEGER ISLCT(*)
      REAL*8 FACT,ANUM
      LOGICAL LMASS,LWEIG
      REAL*8 XCM,YCM,ZCM
      INTEGER NMISS,I
      REAL*8 TMASS,WW,RG,AR

C jmc
      INTEGER*1 OUTU
      OUTU=6

      ANUM=9999.0D0
C
C     Center-of-mass:
C
      XCM=0.0D0
      YCM=0.0D0
      ZCM=0.0D0
      TMASS=0.0D0
      NMISS=0
C
      DO 10 I=1,NATOMS
        IF (ISLCT(I).NE.1) THEN
        ELSE IF (X(I).EQ.ANUM) THEN
          NMISS=NMISS+1
        ELSE
          IF (LMASS.AND.LWEIG) THEN
            WW=W(I)*AM(I)
          ELSE IF (LMASS) THEN
            WW=AM(I)
          ELSE IF (LWEIG) THEN
            WW=W(I)
          ELSE
            WW=1.0D0
          ENDIF
          WW=WW-FACT
          XCM=WW*X(I)+XCM
          YCM=WW*Y(I)+YCM
          ZCM=WW*Z(I)+ZCM
          TMASS=WW+TMASS
        ENDIF
 10   CONTINUE
C
      IF(NMISS.EQ.NATOMS) THEN
C       IF(WRNLEV.GE.2) WRITE(OUTU,25)
        WRITE(OUTU,25)
C       RETURN
      ENDIF
 25   FORMAT(/' RGYR: *** ERROR ***  All coordinates were missing'/)
C
      IF(TMASS.LE.0.0) THEN
C       IF(PRNLEV.GE.2) WRITE(OUTU,35) TMASS
        WRITE(OUTU,35) TMASS
        IF(TMASS.EQ.0.0) RETURN
      ENDIF
 35   FORMAT(/' RGYR: *** WARNING *** Net "mass"=',F12.5/)
C
      XCM=XCM/TMASS
      YCM=YCM/TMASS
      ZCM=ZCM/TMASS
C     IF(NMISS.NE.0 .AND. WRNLEV.GE.2) WRITE(OUTU,45) NMISS
      IF(NMISS.NE.0) WRITE(OUTU,45) NMISS
 45   FORMAT(/' RGYR:   There were',I5,' missing coordinates.'/)
C
C     Radius of gyration:
C
      RG=0.0D0
      DO 60 I=1,NATOMS
        IF (ISLCT(I).NE.1) THEN
        ELSE IF (X(I).EQ.ANUM) THEN
        ELSE
          IF (LMASS.AND.LWEIG) THEN
            WW=W(I)*AM(I)
          ELSE IF (LMASS) THEN
            WW=AM(I)
          ELSE IF (LWEIG) THEN
            WW=W(I)
          ELSE
            WW=1.0D0
          ENDIF
          WW=WW-FACT
          RG=RG+WW*((X(I)-XCM)**2+(Y(I)-YCM)**2+(Z(I)-ZCM)**2)
        ENDIF
 60   CONTINUE
      AR=ABS(RG/TMASS)
C
C     Compute an RG with the same sign as RG/TMASS:
C
      RG=SQRT(AR)*RG/(TMASS*AR)

C     CALL SETMSR('RGYR',RG)
c     CALL SETMSR('MASS',TMASS)
c     CALL SETMSR('XCM ',XCM)
c     CALL SETMSR('YCM ',YCM)
c     CALL SETMSR('ZCM ',ZCM)
C
C     IF(PRNLEV.GE.2) THEN
        IF(LWEIG) WRITE(OUTU,71)
        IF(LMASS) WRITE(OUTU,72)
        IF(.NOT. (LMASS.OR.LWEIG)) WRITE(OUTU,73)
        WRITE(OUTU,74) RG,TMASS,XCM,YCM,ZCM
C     ENDIF
 71   FORMAT(/' RGYR:'/)
 72   FORMAT(/' RGYR:  Mass weighted results:'/)
 73   FORMAT(/' RGYR:  Geometric results:'/)
 74   FORMAT(/'       Radius of gyration=',F12.5,5X,
     $     'Net "mass"=',F12.3/
     $     '       Center-of-"mass" = ' ,3F12.5/)
C
      RETURN
      END
