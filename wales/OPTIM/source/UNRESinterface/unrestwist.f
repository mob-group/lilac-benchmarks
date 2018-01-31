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
C Perturbs one angle as defined by keyword TWISTDIHE nmode dpert
C DMODE refers to dihedrals in order backbone, sidechain.
C
      SUBROUTINE UNRSTWISTDIHE(X,Y,Z,DMODE,DPERT)
      USE COMMONS 
      USE MODUNRES
      IMPLICIT NONE

      DOUBLE PRECISION P,ANGLE,PINT(NVARU)
      REAL*8 X(NATOMS),Y(NATOMS),Z(NATOMS)
      REAL*8 DPERT
      INTEGER DMODE,DMODEMOD,NUMMODES,IICD,I1,J1
      LOGICAL L1,L2,SMARTMOVE
C
      NUMMODES=NPHI+NSIDE

      PRINT *,'NUMMODES,DMODE,DPERT',NUMMODES,DMODE,DPERT
      IF ((DMODE.GT.NUMMODES).OR.(DMODE.LT.-NUMMODES)) THEN
         PRINT *,'Index of mode to be perturbed is larger than number of dihedral angles in protein'
         STOP
      ENDIF
C     
C DPERT should be entered in degrees in odata
C
      IF (DMODE.LT.0) THEN
         ANGLE=-1.0D0*DPERT/57.29577951D0
         DMODEMOD=-DMODE
      ELSE
         ANGLE=DPERT/57.29577951D0
         DMODEMOD=DMODE
      ENDIF

      CALL geom_to_var(nvaru,PINT)
C so PINT now contains old unres internal coordinates (in radians!!). Order is backbone dihedrals, 
C backbone bond angles, sc polar angles then sc dihedral angles (called omega in unres!!!). 
C
C     print *,'tc ANGLE ',ANGLE
      IF(DMODEMOD.LE.NPHI) IICD=DMODEMOD
      IF(DMODEMOD.GT.NPHI) IICD=DMODEMOD+NRES-2+NSIDE

      PINT(IICD)=PINT(IICD)+ANGLE
C     IF (PINT(IICD).LE.-3.141592654D0) PINT(IICD)=PINT(IICD)+2.0D0*3.141592654D0
C     IF (PINT(IICD).GT.3.141592654D0) PINT(IICD)=PINT(IICD)-2.0D0*3.141592654D0

      WRITE (*,'(A20,I3,A2,1X,F10.5,1X,A8)') 'Twisting dihedral ',DMODEMOD,'by',ANGLE,'radians'

C now update the stored internal coordinate arrays and the Cartesian coordinates.
      CALL var_to_geom(nvaru,PINT)
      CALL chainbuild
      DO J1=1,nres
C jmc coords contains x,y,z for all the Calphas
         X(2*J1-1)=c(1,J1)
         Y(2*J1-1)=c(2,J1)
         Z(2*J1-1)=c(3,J1)
C jmc then x,y,z for the side chain centroids
         X(2*J1)=c(1,J1+nres)
         Y(2*J1)=c(2,J1+nres)
         Z(2*J1)=c(3,J1+nres)
      END DO

      RETURN

      END


      SUBROUTINE UNRSPERTDIHE(X,Y,Z,UNPMIN,UNPMAX,UNNMIN,UNNMAX,ISEED)
C adapted from TAKESTEPCH in gmin
      USE COMMONS
      USE MODUNRES
      IMPLICIT NONE

      REAL*8 UNPMIN,UNPMAX,UNNMIN,UNNMAX
      REAL*8 P,ANGLE,DPRAND,RANDOM,STEP,PINT(nvaru)
      INTEGER          a,b,d,NP,NTEST1,NTEST2,I1,J1,RESNUM,ISEED
C
      REAL*8 X(NATOMS),Y(NATOMS),Z(NATOMS)
      REAL*8 XVEC(3*NATOMS),G(3*NATOMS),TWISTARRAY(NATOMS)
      REAL*8 EP1,EP2,E,EMIN1,EMIN2
      INTEGER IICD
      LOGICAL L1,L2,CHIRALFAIL
      LOGICAL TPP(NATOMS),TS(NATOMS)

C set STEP (maximum twist in any one dihedral) to 180.d0
C this could be an input parameter

c     print *,'UNPMIN,UNPMAX,UNNMIN,UNNMAX in unrspertdihe',UNPMIN,UNPMAX,UNNMIN,UNNMAX,ISEED
      STEP=180.d0

C initialise random number generator with input seed
      CALL SDPRND(ISEED)

      CALL geom_to_var(nvaru,PINT)
C
C will be sent back to 192 if too many or too few dihedrals are altered
C as determined by UNNMIN and UNNMAX

192   CONTINUE
C
      b=0
      DO a=1,NPHI
        TPP(a)=.FALSE.
C
C  Calculate P, the probability of twisting
C
        P=UNPMAX
C jmc testing        IF (REAL(a).LE.(0.5*NPHI)) THEN
C jmc testing          P=UNPMAX-a*((UNPMAX-UNPMIN)/(NPHI*0.5))
C jmc testing        ELSE
C jmc testing          P=UNPMIN+(a-0.5*NPHI)*((UNPMAX-UNPMIN)/(NPHI*0.5))
C jmc testing        END IF

        RANDOM=DPRAND()
c       print *,'P, RANDOM',P,RANDOM
        IF (RANDOM.LT.P) THEN
           WRITE (*,'(A,I3)') 'Twisting dihedral ',a
           TPP(a)=.TRUE.
           b=b+1
        END IF
      END DO

C jmc note loop over nres not nside, to take into account glycines and capping groups.
      DO a=1,NRES
        TS(a)=.FALSE.
C Skip if glycine ! jmc can use CYCLE (f90)
        IF (itype(a).EQ.10) GOTO 100

C  Calculate P, the probability of twisting
         P=UNPMAX
C        IF (REAL(a).LE.(0.5E0*(NRES+1))) THEN
C          P=UNPMAX-(REAL(a)-1.E0)*(2.0E0*(UNPMAX-UNPMIN)/(NRES-1.E0))
C        ELSE
C          P=UNPMIN+(REAL(a)-(0.5E0*(NRES+1.E0)))*(2.0E0*(UNPMAX-UNPMIN)/(NRES-1.E0))
C        END IF

C          IF (REAL(RESNUM).LE.(0.5*NRES)) THEN
C            P=CHPMAX-a*((CHPMAX-CHPMIN)/(NRES*0.5))
c          ELSE
c            P=CHPMIN+(a-0.5*NRES)*((CHPMAX-CHPMIN)/(NRES*0.5))
c          END IF
C
        RANDOM=DPRAND()
        IF (RANDOM.LT.P) THEN
           WRITE (*,'(A,I3)') 'Twisting dihedral ',a
           TS(a)=.TRUE.
           b=b+1
        END IF
100   CONTINUE
      END DO
C
c        print *,'UNNMIN UNNMAX',UNNMIN,UNNMAX
C        print *,'NPHIPSI NSIDECHAIN',NPHIPSI,NSIDECHAIN

      NTEST1=INT(UNNMIN*(NPHI+NSIDE))
      IF (NTEST1.LT.1) NTEST1=1
      NTEST2=INT(UNNMAX*(NPHI+NSIDE))

c     WRITE (*,'(A,I3,A,I3,A)') 'Must shift between ',NTEST1,' and ',NTEST2,' dihedrals'
c     WRITE (*,'(A,I3)') 'Attempting to shift ',b

      IF (b.LT.NTEST1 .OR. b.GT.NTEST2) THEN
         WRITE (*,'(A)') 'Too many dihedrals shifted - retrying'
        GOTO 192
      END IF

      DO a=1,NPHI
         IF (TPP(a)) THEN
            ANGLE=(DPRAND()-0.5D0)*2.0D0*STEP/57.29577951D0
            PINT(a)=PINT(a)+ANGLE
            print *,'twisting dihe ',a,' by ',ANGLE
c              IF (PINT(a).LE.-3.141592654D0) PINT(a)=PINT(a)+2.0D0*3.141592654D0
c              IF (PINT(a).GT.3.141592654D0) PINT(a)=PINT(a)-2.0D0*3.141592654D0
         ENDIF
      ENDDO

      DO a=1,NRES
         IF (TS(a)) THEN
            IICD=a+2*NRES+NSIDE-5
            ANGLE=(DPRAND()-0.5D0)*2.0D0*STEP/57.29577951D0
            PINT(IICD)=PINT(IICD)+ANGLE
            print *,'twisting dihe ',a,' by ',ANGLE
c              IF (PINT(IICD).LE.-3.141592654D0) PINT(IICD)=PINT(IICD)+2.0D0*3.141592654D0
c              IF (PINT(IICD).GT.3.141592654D0) PINT(IICD)=PINT(IICD)-2.0D0*3.141592654D0
         ENDIF
      ENDDO

C now update the stored internal coordinate arrays and the Cartesian coordinates.
      CALL var_to_geom(nvaru,PINT)
      CALL chainbuild
      DO J1=1,nres
C jmc coords contains x,y,z for all the Calphas
         X(2*J1-1)=c(1,J1)
         Y(2*J1-1)=c(2,J1)
         Z(2*J1-1)=c(3,J1)
C jmc then x,y,z for the side chain centroids
         X(2*J1)=c(1,J1+nres)
         Y(2*J1)=c(2,J1+nres)
         Z(2*J1)=c(3,J1+nres)
      END DO

C jmc commented these lines out...
C     DO J1=1,NATOMS
C        print *,'pd X',X(J1)
C     ENDDO

      RETURN

      END

C Perturbs one angle as defined by keyword TWISTALL nmode dpert
C DMODE refers to angles in order backbone dihedrals, backbone bond angles, sc polar angles then sc dihedral angles
C
      SUBROUTINE UNRSTWISTALL(X,Y,Z,DMODE,DPERT)
      USE COMMONS
      USE MODUNRES
      IMPLICIT NONE

      DOUBLE PRECISION P,ANGLE,PINT(NVARU)
      REAL*8 X(NATOMS),Y(NATOMS),Z(NATOMS)
      REAL*8 DPERT
      INTEGER DMODE,DMODEMOD,NUMMODES,IICD,I1,J1
      LOGICAL L1,L2,SMARTMOVE
C
      NUMMODES=NPHI+NTHETA+2*NSIDE

      PRINT *,'NUMMODES,DMODE,DPERT',NUMMODES,DMODE,DPERT
      IF ((DMODE.GT.NUMMODES).OR.(DMODE.LT.-NUMMODES)) THEN
         PRINT *,'Index of mode to be perturbed is larger than number of angles in protein'
         STOP
      ENDIF
C     
C DPERT should be entered in degrees in odata
C
      IF (DMODE.LT.0) THEN
         ANGLE=-1.0D0*DPERT/57.29577951D0
         DMODEMOD=-DMODE
      ELSE
         ANGLE=DPERT/57.29577951D0
         DMODEMOD=DMODE
      ENDIF

      CALL geom_to_var(nvaru,PINT)
C so PINT now contains old unres internal coordinates (in radians!!). Order is backbone dihedrals, 
C backbone bond angles, sc polar angles then sc dihedral angles (called omega in unres!!!). 
C
C     print *,'tc ANGLE ',ANGLE

      PINT(DMODEMOD)=PINT(DMODEMOD)+ANGLE
C jmc call to var_to_geom sorts out periodicities for us, so no need to impose restrictions here.

      WRITE (*,'(A20,I3,A2,1X,F10.5,1X,A8)') 'Twisting angle ',DMODEMOD,' by ',ANGLE,' radians'

C now update the stored internal coordinate arrays and the Cartesian coordinates.
      CALL var_to_geom(nvaru,PINT)
      CALL chainbuild
      DO J1=1,nres
C jmc coords contains x,y,z for all the Calphas
         X(2*J1-1)=c(1,J1)
         Y(2*J1-1)=c(2,J1)
         Z(2*J1-1)=c(3,J1)
C jmc then x,y,z for the side chain centroids
         X(2*J1)=c(1,J1+nres)
         Y(2*J1)=c(2,J1+nres)
         Z(2*J1)=c(3,J1+nres)
      END DO

      RETURN

      END
