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
C Routine to guess transition states for unres (despite it's name...) 
C by interpolating between different internal coordinate dihedral values
C Designed to replace neb routine (which is called from connect).
C Equivalent to dae's chguessts for charmm!
C
      SUBROUTINE UNRESGUESSTS(Q,ITEST,PTEST,TWISTTYPE,TWISTFRAC,GUESSFAIL,DISTPF)
      USE COMMONS
      USE KEY
      USE MODTWOEND
      USE MODUNRES
      IMPLICIT NONE
C
      DOUBLE PRECISION ANGLE,TWISTFRAC,Q(3*NATOMS)
C
      REAL*8 DIFFPP,SAVEDIFFPP,MAXDIFF2, DUMMYA, RAND, SUMDIFF, DPRAND
C jmc changed dimension of the following three arrays... Was mxatms.
      REAL*8 FINPPSANGLE(4*nres-9),QPPSANGLE(4*nres-9),DIFFARRAY(4*nres-9),DISTPF
      INTEGER IMIN1,IMIN2,IICD,TWISTMODE,TWISTTYPE,NM,NWRONG,I1,J1
      LOGICAL PTEST,ITEST,RANDOM,NORANDOM,GUESSFAIL
      CHARACTER(LEN=18) GUESSFNAME

      DOUBLE PRECISION PI
      PARAMETER (PI=3.141592653589793D0)
      CHARACTER*15 DIHENAME
      DOUBLE PRECISION ENERGY,RMS,GRAD(3*NATOMS)
      INTEGER NWRONGPOL,TWISTMODEPOL,SAVEDIFFPPPOL,NANGLE

      IF (FILTH.EQ.0) THEN
         GUESSFNAME='unguessts.xyz'
      ELSE
         WRITE(GUESSFNAME,'(A)') 'unguessts.xyz.'//TRIM(ADJUSTL(FILTHSTR))
      ENDIF

C     OPEN(78,FILE='chguessts.xyz',STATUS='UNKNOWN')
      OPEN(78,FILE=GUESSFNAME,STATUS='UNKNOWN')
      CALL UNRESDUMP2(Q,78)
C
      DIFFARRAY=0.0D0 ! jmc initialising

      DO I1=1,nres
         c(1,I1)=FIN(6*(I1-1)+1)
         c(2,I1)=FIN(6*(I1-1)+2)
         c(3,I1)=FIN(6*(I1-1)+3)
         c(1,I1+nres)=FIN(6*(I1-1)+4)
         c(2,I1+nres)=FIN(6*(I1-1)+5)
         c(3,I1+nres)=FIN(6*(I1-1)+6)
c     PRINT *,'FIN in unresguessts: ',FIN(6*(I1-1)+1),FIN(6*(I1-1)+2),FIN(6*(I1-1)+3)
c     PRINT *,'FIN in unresguessts: ',FIN(6*(I1-1)+4),FIN(6*(I1-1)+5),FIN(6*(I1-1)+6)
      END DO
      CALL UPDATEDC
      CALL int_from_cart(.true.,.false.)

      DO I1=1,nres-3
        FINPPSANGLE(I1)=phi(I1+3)
      END DO
      DO I1=1,nres-2
        FINPPSANGLE(I1+nres-3)=omeg(I1+1)
C jmc 30/4/03 try adding backbone and side chain polar angles to the interpolation procedure...
C This should be more important for unres than for charmm...
C Order is bb dihedrals, sc dihedrals, bb bond angles, sc polars.
        FINPPSANGLE(I1+2*nres-5)=theta(I1+2)
        FINPPSANGLE(I1+3*nres-7)=alph(I1+1)
      END DO

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

C use unres geometry arrays phi (bb dihedrals) and omeg (sc dihedrals)
C NOTE THAT ANGLES ARE IN RADIANS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C Take care with numbering - see /unres/src/readpdb.f (subroutine int_from_cart)
C For side chain dihedrals, the actual stored arrays (alpha and omeg) contain zero elements for
C proper (i.e. not capping) glycines but the variable array from a call to geom_to_var does not
C contain these elements.
C Need to remember not to try to twist them though!
C No entries in QPPSANGLE for capping 'residues'.
      DO I1=1,nres-3
        QPPSANGLE(I1)=phi(I1+3)
      END DO
      DO I1=1,nres-2
        QPPSANGLE(I1+nres-3)=omeg(I1+1)
        QPPSANGLE(I1+2*nres-5)=theta(I1+2)
        QPPSANGLE(I1+3*nres-7)=alph(I1+1)
      END DO
C jmc note that the Q internal coord set is now saved in the unres int coor common block...
C Put the TS guess coords into common block before exiting this subroutine.

C Now decide which phi/psi or sidechain angle to twist
C
C Based on TWISTTYPE
C TWISTTYPE = 1  means take one with biggest difference and interpolate between
C      the two values using TWISTFRAC as the fraction
C TWISTTYPE = 2; interpolates like 1 but sets chosen angle to the nearest of -120, 0, 120 degrees
C i.e. maxima of the dihedral potential (which is k(1+cos(3*phi)) for phi and psi angles. ! charmm
C In fact k = 0 for psi, so this method may be a bit futile for psi angles, but it may give sensible ! charmm
C geometries anyway) ! charmm
C
C TWISTTYPE =3; like 1, but also interpolates the dihedral either side of the maximum
C
C TWISTTYPE =4; does on one dihedral, chosen with probability based on size of displacement
C
C TWISTTYPE =5; If only one dihedral differs by >60deg then interpolates on one dihedral,
C               If more than ones differs then proceeds like random mode (TWISTTYPE =4)
C
C TWISTTYPE =6; If only one dihedral differs by >60deg then interpolates that dihedral,
C               and the dihedrals either side (like TT=3)
C               If more than ones differs then proceeds like random mode (TWISTTYPE =4)
C
C TWISTTYPE =7; Just interpolate all dihedrals
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C jmc NOTE for now, only these options below are working!!!
C !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
C TWISTTYPE =8; Interpolate all backbone angles
C TWISTTYPE =9; Interpolate largest dihedral and largest polar angle
C TWISTTYPE =10; Just interpolate all angles
C
      RANDOM=.FALSE.
      MAXDIFF2 = 0.0D0
      NWRONG=0
      NWRONGPOL=0
C
C turn random displacements off once two minima are close enough together for neb
C to be successful
C
c     IF (DISTPF.LT.RANDOMCUTOFF) THEN
c        NORANDOM=.TRUE.
c     ELSE
         NORANDOM=.FALSE.
c     ENDIF

      DO I1=1,nphi+nres-2
         IF (I1.GT.nphi) THEN
            IF (itype(I1-nphi+1).EQ.10) GOTO 100 ! glycine
         END IF
c        WRITE(*,'(A,I6,2F15.10)') 'FINS QS',I1,FINPPSANGLE(I1),QPPSANGLE(I1)
         DIFFPP = FINPPSANGLE(I1) - QPPSANGLE(I1)
C
C next two lines are meant to ensure that you always interpolate
C along the shortest distance between the dihedral angles.
C
         IF (DIFFPP.GT.PI) DIFFPP = DIFFPP-2.0D0*PI
         IF (DIFFPP.LT.-PI) DIFFPP = DIFFPP+2.0D0*PI

         DIFFARRAY(I1)=DIFFPP

         IF ((DIFFPP*DIFFPP).GT.MAXDIFF2) THEN
            MAXDIFF2=DIFFPP*DIFFPP
            TWISTMODE=I1
            SAVEDIFFPP=DIFFPP
         ENDIF

C jmc         IF (ABS(DIFFPP).GT.60.0D0) NWRONG=NWRONG+1
         IF (ABS(DIFFPP).GT.PI/3.0D0) NWRONG=NWRONG+1

100   CONTINUE 
      ENDDO

C jmc don't duplicate work from above do loop...
C Remember the polar angles run from 0 to pi, whereas dihedrals go from -pi to pi.
      DO I1=nphi+nres-1,nphi+ntheta+2*nres-4
         IF (I1.GT.nphi+nres-2+ntheta) THEN
            IF (itype(I1-nphi-nres+2-ntheta+1).EQ.10) GOTO 200 ! glycine
         END IF
c        WRITE(*,'(A,I6,2F15.10)') 'FINS QS',I1,FINPPSANGLE(I1),QPPSANGLE(I1)
         DIFFPP = FINPPSANGLE(I1) - QPPSANGLE(I1)
         DIFFARRAY(I1)=DIFFPP
         IF ((DIFFPP*DIFFPP).GT.MAXDIFF2) THEN
            MAXDIFF2=DIFFPP*DIFFPP
            TWISTMODEPOL=I1
            SAVEDIFFPPPOL=DIFFPP
         ENDIF

C jmc         IF (ABS(DIFFPP).GT.60.0D0) NWRONG=NWRONG+1
         IF (ABS(DIFFPP).GT.PI/3.0D0) NWRONGPOL=NWRONGPOL+1

200   CONTINUE
      ENDDO

C
C Now do twisting 
C

      IF (TWISTTYPE.EQ.7) THEN
         DO I1=1,NPHI
            ANGLE=TWISTFRAC*DIFFARRAY(I1)
C jmc            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*(2.0D0*PI-DIFFARRAY(I1))
            IF (TWISTFRAC.LT.0.D0) THEN
               ANGLE = 0.5D0*(2.0D0*PI-ABS(DIFFARRAY(I1))) ! so angle will always be between pi/2 and pi (always > 0)
               IF (DIFFARRAY(I1).GT.0.0D0) ANGLE = -ANGLE ! jmc need to test this!! or do we need -1.0D0*ANGLE??
            ENDIF
            phi(I1+3)=phi(I1+3)+ANGLE
         ENDDO
         DO I1=NPHI+1,NPHI+NRES-2
            ANGLE=TWISTFRAC*DIFFARRAY(I1)
            IF (TWISTFRAC.LT.0.D0) THEN
               ANGLE = 0.5D0*(2.0D0*PI-ABS(DIFFARRAY(I1))) ! so angle will always be between pi/2 and pi (always > 0)
               IF (DIFFARRAY(I1).GT.0.0D0) ANGLE = -ANGLE ! jmc need to test this!! or do we need -1.0D0*ANGLE??
            ENDIF
C jmc            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*(2.0D0*PI-DIFFARRAY(I1))
            omeg(I1+1-nphi)=omeg(I1+1-nphi)+ANGLE
         ENDDO
         GOTO 20
      ENDIF

      IF (TWISTTYPE.EQ.10) THEN
         DO I1=1,NPHI
            ANGLE=TWISTFRAC*DIFFARRAY(I1)
C jmc            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*(2.0D0*PI-DIFFARRAY(I1))
            IF (TWISTFRAC.LT.0.D0) THEN
               ANGLE = 0.5D0*(2.0D0*PI-ABS(DIFFARRAY(I1))) ! so angle will always be between pi/2 and pi (always > 0)
               IF (DIFFARRAY(I1).GT.0.0D0) ANGLE = -ANGLE ! jmc need to test this!! or do we need -1.0D0*ANGLE??
            ENDIF
            phi(I1+3)=phi(I1+3)+ANGLE
         ENDDO
         DO I1=NPHI+1,NPHI+NRES-2
            ANGLE=TWISTFRAC*DIFFARRAY(I1)
C jmc            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*(2.0D0*PI-DIFFARRAY(I1))
            IF (TWISTFRAC.LT.0.D0) THEN
               ANGLE = 0.5D0*(2.0D0*PI-ABS(DIFFARRAY(I1))) ! so angle will always be between pi/2 and pi (always > 0)
               IF (DIFFARRAY(I1).GT.0.0D0) ANGLE = -ANGLE ! jmc need to test this!! or do we need -1.0D0*ANGLE??
            ENDIF
            omeg(I1+1-nphi)=omeg(I1+1-nphi)+ANGLE
         ENDDO
         DO I1=NPHI+NRES-1,NPHI+NRES-2+NTHETA
            ANGLE=TWISTFRAC*DIFFARRAY(I1)
            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*DIFFARRAY(I1) ! 'going the long way round' doesn't apply for bond angles
            theta(I1-NPHI-NRES+4)=theta(I1-NPHI-NRES+4)+ANGLE
         ENDDO
         DO I1=NPHI+NRES-1+NTHETA,NPHI+2*NRES-4+NTHETA
            ANGLE=TWISTFRAC*DIFFARRAY(I1)
            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*DIFFARRAY(I1) ! 'going the long way round' doesn't apply for bond angles
            alph(I1-NPHI-NTHETA-NRES+3)=alph(I1-NPHI-NTHETA-NRES+3)+ANGLE
         ENDDO
         GOTO 20
      ENDIF

      IF (TWISTTYPE.EQ.8) THEN
C jmc backbone angles only
         DO I1=1,NPHI
            ANGLE=TWISTFRAC*DIFFARRAY(I1)
C jmc            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*(2.0D0*PI-DIFFARRAY(I1))
            IF (TWISTFRAC.LT.0.D0) THEN
               ANGLE = 0.5D0*(2.0D0*PI-ABS(DIFFARRAY(I1))) ! so angle will always be between pi/2 and pi (always > 0)
               IF (DIFFARRAY(I1).GT.0.0D0) ANGLE = -ANGLE ! jmc need to test this!! or do we need -1.0D0*ANGLE??
            ENDIF
            phi(I1+3)=phi(I1+3)+ANGLE
         ENDDO
         DO I1=NPHI+NRES-1,NPHI+NRES-2+NTHETA
            ANGLE=TWISTFRAC*DIFFARRAY(I1)
            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*DIFFARRAY(I1) ! 'going the long way round' doesn't apply for bond angles
            theta(I1-NPHI-NRES+4)=theta(I1-NPHI-NRES+4)+ANGLE
         ENDDO
         GOTO 20
      ENDIF

      IF (TWISTTYPE.EQ.9) THEN
            ANGLE=TWISTFRAC*DIFFARRAY(TWISTMODE)
C jmc            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*(2.0D0*PI-DIFFARRAY(TWISTMODE))
            IF (TWISTFRAC.LT.0.D0) THEN
               ANGLE = 0.5D0*(2.0D0*PI-ABS(DIFFARRAY(TWISTMODE))) ! so angle will always be between pi/2 and pi (always > 0)
               IF (DIFFARRAY(TWISTMODE).GT.0.0D0) ANGLE = -ANGLE ! jmc need to test this!! or do we need -1.0D0*ANGLE??
            ENDIF
            phi(TWISTMODE+3)=phi(TWISTMODE+3)+ANGLE

            ANGLE=TWISTFRAC*DIFFARRAY(TWISTMODE+nphi+nres-2)
            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*DIFFARRAY(TWISTMODE+nphi+nres-2) ! 'going the long way round' doesn't apply for bond angles
            theta(TWISTMODE-NPHI-NRES+4)=theta(TWISTMODE-NPHI-NRES+4)+ANGLE
            ANGLE=TWISTFRAC*DIFFARRAY(TWISTMODE+1+nphi+nres-2)
            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*DIFFARRAY(TWISTMODE+1+nphi+nres-2) ! 'going the long way round' doesn't apply for bond angles
            theta(TWISTMODE+1-NPHI-NRES+4)=theta(TWISTMODE+1-NPHI-NRES+4)+ANGLE
         GOTO 20
      ENDIF

      IF ((TWISTTYPE.EQ.5).OR.(TWISTTYPE.EQ.6)) THEN
         IF (NWRONG.GT.2) THEN
C            WRITE (*,'(A)') 'More than one dihedral displaced - unlikely to be a direct connection'
            WRITE (*,'(A)') 'More than two dihedrals displaced - unlikely to be a direct connection'
            IF (NORANDOM) THEN
               WRITE (*,'(A)') 'Switching to neb'
               GUESSFAIL=.TRUE.
               RETURN
            ELSE
               WRITE (*,'(A)') 'Choosing a mode to twist at random'
               RANDOM=.TRUE.
c              STOP
            ENDIF
         ENDIF
      ENDIF

      ANGLE=TWISTFRAC*SAVEDIFFPP
C jmc what if savediffpp is lt 0?
      IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*(2.0D0*PI-SAVEDIFFPP) ! this will not work properly - don't use!!!
      IF (TWISTTYPE.EQ.2) THEN
         DUMMYA=QPPSANGLE(TWISTMODE)+ANGLE
C jmc         IF ((DUMMYA.GT.-180.0D0).AND.(DUMMYA.LT.-60.0D0)) DUMMYA=-120.0D0
C jmc         IF ((DUMMYA.GT.-60.0D0).AND.(DUMMYA.LT.60.0D0)) DUMMYA=0.0D0
C jmc         IF ((DUMMYA.GT.60.0D0).AND.(DUMMYA.LT.180.0D0)) DUMMYA=120.0D0
         IF ((DUMMYA.GT.-PI).AND.(DUMMYA.LT.-PI/3.0D0)) DUMMYA=-2.0D0*PI/3.0D0
         IF ((DUMMYA.GE.-PI/3.0D0).AND.(DUMMYA.LT.PI/3.0D0)) DUMMYA=0.0D0
         IF ((DUMMYA.GE.PI/3.0D0).AND.(DUMMYA.LE.PI)) DUMMYA=2.0D0*PI/3.0D0
         ANGLE=DUMMYA-QPPSANGLE(TWISTMODE)
      ENDIF
      
      IF ((TWISTTYPE.EQ.4).OR.RANDOM) THEN
         SUMDIFF=0.D0
         DO I1=1,nphi+nres-2
            SUMDIFF=SUMDIFF+ABS(DIFFARRAY(I1))
         ENDDO
         RAND=DPRAND()*SUMDIFF
         print *,'RAND',RAND
         SUMDIFF=0.D0
         DO I1=1,NPHI+nres-2
            SUMDIFF=SUMDIFF+ABS(DIFFARRAY(I1))
c              PRINT *,'DIFFARRAY ',DIFFARRAY(I1)
            IF (SUMDIFF.GT.RAND) THEN 
               TWISTMODE=I1
               ANGLE=TWISTFRAC*DIFFARRAY(I1)
               PRINT *,'ANGLE ',ANGLE
               IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5*(2.0D0*PI-SAVEDIFFPP)
               GOTO 10
            ENDIF
         ENDDO
10      CONTINUE
C jmc huh?
C       IF (RANDOM) THEN
C jmc           ANGLE=DPRAND()*60.D0
c          ANGLE=DPRAND()*PI/3.0D0
c          IF (DIFFARRAY(TWISTMODE).LT.0.D0) ANGLE=-1.D0*ANGLE
c       ENDIF
      ENDIF

C jmc      WRITE (*,'(A20,I3,A2,1X,F10.5,1X,A8)') 'Twisting phi/psi dihedral ',TWISTMODE,' by ',ANGLE,' degrees'
      WRITE (*,'(A20,I3,A2,1X,F10.5,1X,A8)') 'Twisting phi/psi dihedral ',TWISTMODE,' by ',ANGLE,' radians'

      IF (TWISTMODE.LE.NPHI) THEN
         phi(TWISTMODE+3)=phi(TWISTMODE+3)+ANGLE
      ELSE
         omeg(TWISTMODE+1-nphi)=omeg(TWISTMODE+1-nphi)+ANGLE
      END IF

      IF ((TWISTTYPE.EQ.3).OR.((TWISTTYPE.EQ.6).AND.(.NOT.RANDOM))) THEN
         NM=TWISTMODE-1
         IF (NM.GE.1) THEN
            DIFFPP = FINPPSANGLE(NM) - QPPSANGLE(NM)
C jmc            IF (DIFFPP.GT.180.0) DIFFPP = DIFFPP-360.D0
C jmc            IF (DIFFPP.GT.180.0) DIFFPP = DIFFPP-360.D0
            IF (DIFFPP.LT.-PI) DIFFPP = DIFFPP+2.0D0*PI
            IF (DIFFPP.GT.PI) DIFFPP = DIFFPP-2.0D0*PI
            ANGLE=TWISTFRAC*DIFFPP
            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*(2.0D0*PI-SAVEDIFFPP)
C jmc            IICD=PHIPSI(NM)
C jmc            CALL TWISTCH(IICD,ANGLE)
            phi(NM+3)=phi(NM+3)+ANGLE
         ENDIF
         NM=TWISTMODE+1
C jmc         IF (NM.LE.NPHIPSI) THEN
         IF (NM.LE.NPHI) THEN
            DIFFPP = FINPPSANGLE(NM) - QPPSANGLE(NM)
            IF (DIFFPP.GT.PI) DIFFPP = DIFFPP-2.0D0*PI
            IF (DIFFPP.LT.-PI) DIFFPP = DIFFPP+2.0D0*PI
            ANGLE=TWISTFRAC*DIFFPP
            IF (TWISTFRAC.LT.0.D0) ANGLE = 0.5D0*(2.0D0*PI-SAVEDIFFPP)
C jmc            IICD=PHIPSI(NM)
C jmc            CALL TWISTCH(IICD,ANGLE)
            phi(NM+3)=phi(NM+3)+ANGLE
         ENDIF
       ENDIF
C
20    CONTINUE

      CALL chainbuild

      DO J1=1,nres
         Q(6*(J1-1)+1)=c(1,J1)
         Q(6*(J1-1)+2)=c(2,J1)
         Q(6*(J1-1)+3)=c(3,J1)
         Q(6*(J1-1)+4)=c(1,J1+nres)
         Q(6*(J1-1)+5)=c(2,J1+nres)
         Q(6*(J1-1)+6)=c(3,J1+nres)
      END DO

      CALL UNRESDUMP2(Q,78)
      CALL UNRESDUMP2(FIN,78)

      CLOSE(78)

C jmc testing
C     CALL POTENTIAL(Q,ENERGY,GRAD,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
C     PRINT *,'ENERGY,RMS IN CHGUESSTS ',ENERGY,RMS
C     DIHENAME='QTS'
C     CALL PRINTDIHE(Q,QPPSANGLE,NANGLE,DIHENAME)

      RETURN

      END

      SUBROUTINE UNRESGUESSMIN(Q,PTEST,TWISTTYPE,NGUESS)
      USE COMMONS
      USE MODTWOEND
      USE MODUNRES
      IMPLICIT NONE

      DOUBLE PRECISION Q(3*NATOMS),TWISTFRAC,DISTPF,TMPQ(3*NATOMS)
      LOGICAL PTEST,GUESSFAIL
      INTEGER TWISTTYPE,I1,NGUESS,J1,K1
      DOUBLE PRECISION Q1(3*NATOMS),Q2(3*NATOMS),Q3(3*NATOMS)
C     COMMON /MINARRAY/ MYQMINSAVE(3*NATOMS,100),COUNTER   ! note arbitrarily choosing 100, could have as a parameter but can't be bothered now...
      CHARACTER*5 ZSYMSAVE
      COMMON /SYS/ ZSYMSAVE

      CALL NEWMINDIST(Q,FIN,NATOMS,DISTPF,.FALSE.,.FALSE.,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
      PRINT *,'Q, FIN ',DISTPF

      TMPQ=Q

      COUNTER=0
      DO I1=1,NGUESS
         TWISTFRAC=1.0D0*I1/(NGUESS+1)
         CALL UNRESGUESSTS(Q,.FALSE.,PTEST,TWISTTYPE,TWISTFRAC,GUESSFAIL,DISTPF)
         Q=TMPQ
         PRINT *,'GUESSFAIL ',GUESSFAIL
      END DO

C now do some funky mind stuff on MYQMINSAVE...
      DO I1=1,COUNTER
         DO K1=1,3*NATOMS
            Q1(K1)=MYQMINSAVE(K1,I1)
            Q2(K1)=TMPQ(K1)
            Q3(K1)=FIN(K1)
         END DO
         CALL NEWMINDIST(Q1,Q2,NATOMS,DISTPF,.FALSE.,.FALSE.,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
         PRINT *,'DISTPF ',DISTPF,I1,' START'
         CALL NEWMINDIST(Q1,Q3,NATOMS,DISTPF,.FALSE.,.FALSE.,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
         PRINT *,'DISTPF ',DISTPF,I1,' FIN'
         DO J1=I1+1,COUNTER
            DO K1=1,3*NATOMS
               Q1(K1)=MYQMINSAVE(K1,I1)
               Q2(K1)=MYQMINSAVE(K1,J1)
            END DO
            CALL NEWMINDIST(Q1,Q2,NATOMS,DISTPF,.FALSE.,.FALSE.,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
            PRINT *,'DISTPF ',DISTPF,I1,J1
         END DO
      END DO

      STOP

      END SUBROUTINE UNRESGUESSMIN
