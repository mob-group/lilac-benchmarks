! 
!  Copyright (C) 1999-2008 David J. Wales                       
!  This file is part of OPTIM.
!
!  OPTIM is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  OPTIM is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program; if not, write to the Free Software
!  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  COORDSA becomes the optimal alignment of the optimal permutation(-inversion)
!  isomer, but without the permutations. DISTANCE is the residual square distance
!  for the best alignment with respect to permutation(-inversion)s as well as
!  orientation and centre of mass.
!
!  MYORIENT is called first for both COORDSA and COORDSB to put them into
!  a standard orientation in DUMMYA and DUMMYB (which both have the centre of
!  coordinates at the origin). 
!  The objective is to identify permutation-inversion isomers without fail. 
!  However, we have to cycle over all equivalent atoms in two particular orbits for DUMMYA
!  to achieve this.
!  We iterate permutations and newmindist minimisations up to a maximum number or
!  until no more permutations are required for each instance of DUMMYA aligned 
!  according to NCHOOSE1 and NCHOOSE2 by MYORIENT. The cumulative rotation
!  matrix that takes the initial DUMMYA to the one that aligns best with DUMMYB
!  is saved in RMATCUMUL.
!  Then, if we've not going BULK, AMBER, or CHARMM, we try again for the inverted
!  version of COORDSA. The transformation corresponding to the minimum distance
!  is saved whenever it is improved - the best alignment including permutations
!  is saved in XBEST, and the last step is to rotate this back to coincide best
!  with COORDSB (rather than DUMMYB) using ROTINVBBEST. This gives suitable
!  fixed end points for DNEB.
!  Finally, we transform COORDSA to be in optimal alignment, with the
!  permutations in XBEST. The overall transformation is
!  COORDSA -> +/- ROTINVB RMATCUMUL ROTA (permutation(COORDSA) - CMA) 
!
!  The correspondence between COORDSA and DUMMYA after DUMMYA has been aligned by
!  newmindist is
!  +/- RMATCUMUL ROTA (COORDSA - CMA) = permutation(DUMMYA)
!  where +/- is given by the value of INVERT.
!  The centres of coordinates for COORDSA and COORDSB can be anywhere. On return, the
!  centre of coordinates of COORDSA will be the same as for COORDSB, unless we
!  are doing an ion trap potential.
!
SUBROUTINE MINPERMDIST(COORDSB,COORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGID,RMATBEST,USEINT)
USE COMMONS,ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, SETS, GEOMDIFFTOL, AMBERT, NFREEZE, CHARMMT, RBAAT, PULLT, &
  &               ANGLEAXIS, PERMISOMER, PERMDIST, ZSYM, INTCONSTRAINTT, INTLJT, OHCELLT, ATOMMATCHDIST, LPERMDIST, &
  &               TRAPT, NRANROT, MACROCYCLET, MCYCLEPERIOD, MCYCLEREPEATS, NOINVERSION, GTHOMSONT, &
  &               RATETARGETT, MLP3T, NOPT, MKTRAPT, MLPB3T, NOTRANSROTT, MACROIONT
USE PORFUNCS 
USE UTILS,ONLY : GETUNIT

IMPLICIT NONE

INTEGER :: MAXIMUMTRIES=10
INTEGER NATOMS, NPERM, PATOMS, NTRIES, BESTINVERT, I, LUNIT
INTEGER INVERT, NORBIT1, NORBIT2, NCHOOSE2, NDUMMY, LPERM(NATOMS), J1, J2, NCHOOSE1, OPNUM, J3, NROTDONE
INTEGER NORBITB1, NORBITB2, NCHOOSEB1, NCHOOSEB2

DOUBLE PRECISION DIST2, COORDSA(3*NATOMS), COORDSB(3*NATOMS), DISTANCE, DUMMYA(3*NATOMS), DUMMYB(3*NATOMS), DUMMY(3*NATOMS)
DOUBLE PRECISION BOXLX,BOXLY,BOXLZ,WORSTRAD,RMAT(3,3),ENERGY, VNEW(3*NATOMS), DX, DY, DZ, RMS, DBEST, XBEST(3*NATOMS)
DOUBLE PRECISION CMXA, CMXB, CMXC, CMX, CMY, CMZ
DOUBLE PRECISION ROTA(3,3), ROTINVA(3,3), ROTB(3,3), ROTINVB(3,3), ROTINVBBEST(3,3), ROTABEST(3,3), RMATBEST(3,3), TMAT(3,3)
DOUBLE PRECISION CMAX, CMAY, CMAZ, CMBX, CMBY, CMBZ, RMATCUMUL(3,3)
DOUBLE PRECISION REFXZ(3,3), QBEST(4)
DOUBLE PRECISION BMDIST
DOUBLE PRECISION, ALLOCATABLE :: BMCOORDS(:), BMCOORDSSV(:)
LOGICAL DEBUG, TWOD, RIGID, BULKT, PITEST, BMTEST, TNMATCH
LOGICAL, INTENT(IN) :: USEINT
DOUBLE PRECISION PDUMMYA(3*NATOMS), PDUMMYB(3*NATOMS), LDISTANCE, DUMMYC(3*NATOMS), XDUMMY
INTEGER NEWPERM(NATOMS), ALLPERM(NATOMS), SAVEPERM(NATOMS), NMAXINT, NMININT, BESTPERM(NATOMS)
DOUBLE PRECISION CONSTRAINTE, XYZLOCAL(6*NATOMS), LMINCOORDS(2,3*NATOMS)
INTEGER :: MCYCLESTEP =1,MCYCLESHIFT

! hk286
IF (GTHOMSONT) THEN
   CALL GTHOMSONMINPERMDIST(COORDSB,COORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGID,RMATBEST,USEINT)
   RETURN
ELSEIF (MLP3T.OR.MLPB3T) THEN
   DISTANCE=0.0D0
   DO J1=1,NOPT
      DISTANCE=DISTANCE+(COORDSA(J1)-COORDSB(J1))**2
   ENDDO
   DISTANCE=SQRT(DISTANCE)
   RETURN
ENDIF

NROTDONE=-1
MAXIMUMTRIES=MAX(MAXIMUMTRIES,NRANROT+1)
IF (RATETARGETT) CALL SDPRND(0) ! to make demos deterministic

BMTEST=.FALSE.
IF (ATOMMATCHDIST) BMTEST=.TRUE. 
IF (BMTEST) THEN 
 ALLOCATE(BMCOORDS(3*NATOMS))
 ALLOCATE(BMCOORDSSV(3*NATOMS))
ENDIF


IF (USEINT) THEN
   LMINCOORDS(1,1:3*NATOMS)=COORDSB(1:3*NATOMS)
   LMINCOORDS(2,1:3*NATOMS)=COORDSA(1:3*NATOMS)
   XYZLOCAL(1:3*NATOMS)=LMINCOORDS(1,1:3*NATOMS)
   XYZLOCAL(3*NATOMS+1:6*NATOMS)=LMINCOORDS(2,1:3*NATOMS)
   IF (INTCONSTRAINTT) THEN
!     PRINT '(A)','minpermdist> calling make_conpot'
      CALL MAKE_CONPOT(2,LMINCOORDS)
      CALL CONGRAD2(NMAXINT,NMININT,CONSTRAINTE,XYZLOCAL)
   ELSE IF (INTLJT) THEN
      CALL INTGRADLJ(CONSTRAINTE,XYZLOCAL)
   ENDIF
   DISTANCE=CONSTRAINTE/2.0D0
   IF (DEBUG) PRINT '(A,2G20.10)','minpermdist> interpolation metric, constraint energy ',DISTANCE,CONSTRAINTE
!  IF (DEBUG) PRINT '(A,G20.10)','minpermdist> interpolation metric ',DISTANCE
   RETURN
ENDIF

IF (RBAAT) THEN
   IF (PERMDIST) THEN
      CALL RBMINPERMDIST(COORDSB,COORDSA,DISTANCE,DIST2,QBEST,RMATBEST,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT)
   ELSE
      CALL RBMINDIST(COORDSB,COORDSA,NATOMS,DISTANCE,QBEST,DEBUG)
   ENDIF
   RETURN
ELSEIF (LPERMDIST) THEN
   CALL LOPERMDIST(COORDSB,COORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGID,RMATBEST)
   RETURN
ELSEIF (.NOT.(PERMDIST.OR.PERMISOMER)) THEN
   CALL NEWMINDIST(COORDSB,COORDSA,NATOMS,DISTANCE,BULKT,TWOD,ZSYM,.FALSE.,ANGLEAXIS,DEBUG,RMAT)
   RETURN
ENDIF

REFXZ(1:3,1:3)=0.0D0
REFXZ(1,1)=1.0D0; REFXZ(2,2)=-1.0D0; REFXZ(3,3)=1.0D0

! 
! DUMMYB(1:3*NATOMS)=COORDSB(1:3*NATOMS)
! DUMMYA(1:3*NATOMS)=COORDSA(1:3*NATOMS)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CALL OCHARMM(DUMMYA,VNEW,ENERGY,.FALSE.,.FALSE.)
! CALL POTENTIAL(DUMMYA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
! PRINT '(2(A,F25.15))',' Initial energy=',ENERGY,' RMS=',RMS
! PRINT '(2(A,F25.15))',' for coordinates:'
! PRINT '(3F25.15)',DUMMYA(1:3*NATOMS)
! PRINT '(A,F25.15,A)',' Initial energy=',ENERGY,' kcal/mol'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!  LUNIT=GETUNIT()
!  OPEN(UNIT=LUNIT,FILE='minpermdist.xyz',STATUS='UNKNOWN')
!  WRITE(LUNIT,'(I6)') NATOMS/2
!  WRITE(LUNIT,'(A)') 'A initial'
!  DO J3=1,NATOMS/2
!      WRITE(LUNIT,'(A2,2X,3F20.10)') 'LA',COORDSA(3*(J3-1)+1),COORDSA(3*(J3-1)+2),COORDSA(3*(J3-1)+3)
!  ENDDO
!  WRITE(LUNIT,'(I6)') NATOMS/2
!  WRITE(LUNIT,'(A)') 'B initial'
!  DO J3=1,NATOMS/2
!      WRITE(LUNIT,'(A2,2X,3F20.10)') 'LA',COORDSB(3*(J3-1)+1),COORDSB(3*(J3-1)+2),COORDSB(3*(J3-1)+3)
!  ENDDO
!  CLOSE(LUNIT)
!
!  Calculate original centres of mass.
!
CMAX=0.0D0; CMAY=0.0D0; CMAZ=0.0D0
IF ((NFREEZE.LE.0).AND.(.NOT.MKTRAPT).AND.(.NOT. NOTRANSROTT)) THEN
   IF (RBAAT) THEN
      DO J1=1,NATOMS/2
         CMAX=CMAX+COORDSA(3*(J1-1)+1)
         CMAY=CMAY+COORDSA(3*(J1-1)+2)
         CMAZ=CMAZ+COORDSA(3*(J1-1)+3)
      ENDDO
      CMAX=2*CMAX/NATOMS; CMAY=2*CMAY/NATOMS; CMAZ=2*CMAZ/NATOMS
      CMBX=0.0D0; CMBY=0.0D0; CMBZ=0.0D0
      DO J1=1,NATOMS/2
         CMBX=CMBX+COORDSB(3*(J1-1)+1)
         CMBY=CMBY+COORDSB(3*(J1-1)+2)
         CMBZ=CMBZ+COORDSB(3*(J1-1)+3)
      ENDDO
      CMBX=2*CMBX/NATOMS; CMBY=2*CMBY/NATOMS; CMBZ=2*CMBZ/NATOMS
   ELSE
      DO J1=1,NATOMS
         CMAX=CMAX+COORDSA(3*(J1-1)+1)
         CMAY=CMAY+COORDSA(3*(J1-1)+2)
         CMAZ=CMAZ+COORDSA(3*(J1-1)+3)
      ENDDO
      CMAX=CMAX/NATOMS; CMAY=CMAY/NATOMS; CMAZ=CMAZ/NATOMS
      CMBX=0.0D0; CMBY=0.0D0; CMBZ=0.0D0
      DO J1=1,NATOMS
         CMBX=CMBX+COORDSB(3*(J1-1)+1)
         CMBY=CMBY+COORDSB(3*(J1-1)+2)
         CMBZ=CMBZ+COORDSB(3*(J1-1)+3)
      ENDDO
      CMBX=CMBX/NATOMS; CMBY=CMBY/NATOMS; CMBZ=CMBZ/NATOMS
   ENDIF
ENDIF
!
! It is possible for the standard orientation to result in a distance that is worse than
! the starting distance. Hence we need to set XBEST here.
!
DUMMYA(1:3*NATOMS)=COORDSA(1:3*NATOMS)
DUMMYB(1:3*NATOMS)=COORDSB(1:3*NATOMS)
DBEST=1.0D100
CALL NEWMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,BULKT,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT)
DBEST=DISTANCE**2
! IF (DEBUG) PRINT '(A,G20.10)',' minpermdist> Initial distance before standard orientation=',DISTANCE

! If global translation is possible (i.e. no frozen atoms and not periodic boundary conditions) then
! translate the origin of structure A to the CoM of structure B.
IF ((NFREEZE.GT.0).OR.BULKT.OR.MKTRAPT.OR.NOTRANSROTT) THEN
   XBEST(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
ELSE
   DO J1=1,NATOMS
      XBEST(3*(J1-1)+1)=DUMMYA(3*(J1-1)+1)-CMBX
      XBEST(3*(J1-1)+2)=DUMMYA(3*(J1-1)+2)-CMBY
      XBEST(3*(J1-1)+3)=DUMMYA(3*(J1-1)+3)-CMBZ
   ENDDO
ENDIF
BESTINVERT=1
DO J1=1,NATOMS
   BESTPERM(J1)=J1
ENDDO
RMATBEST(1:3,1:3)=RMAT(1:3,1:3)
ROTINVBBEST(1:3,1:3)=0.0D0
ROTINVBBEST(1,1)=1.0D0;ROTINVBBEST(2,2)=1.0D0;ROTINVBBEST(3,3)=1.0D0;
ROTABEST(1:3,1:3)=0.0D0
ROTABEST(1,1)=1.0D0;ROTABEST(2,2)=1.0D0;ROTABEST(3,3)=1.0D0;
!
! End of XBEST associated initialisation.
! 
NROTDONE=-1
11 CONTINUE
NROTDONE=NROTDONE+1

MCYCLESTEP=1
70 MCYCLESHIFT=(MCYCLESTEP-1)*MCYCLEPERIOD !Macrocycle loop point
INVERT=1
60 CONTINUE ! jump back here if INVERT changes sign.
   NCHOOSEB1=0
66 NCHOOSEB1=NCHOOSEB1+1
   NCHOOSEB2=0
31 NCHOOSEB2=NCHOOSEB2+1
   NCHOOSE1=0
65 NCHOOSE1=NCHOOSE1+1
40 NCHOOSE2=0
30 NCHOOSE2=NCHOOSE2+1
OPNUM=0
IF (BMTEST) THEN
  BMDIST=HUGE(1.0D0)
  BMCOORDS(1:3*NATOMS)=COORDSA(1:3*NATOMS)
ENDIF
TNMATCH=.FALSE.
25 OPNUM=OPNUM+1 ! Point group operation counter for Oh supercell if OHCELLT is true.
DUMMYB(1:3*NATOMS)=COORDSB(1:3*NATOMS)
IF (OHCELLT) THEN
!  IF (DEBUG) PRINT '(A,I8)','minpermdist> Trying Oh symmetry operation number ',OPNUM
   CALL OHOPS(COORDSA,DUMMYA,OPNUM,NATOMS)
ELSE
   DUMMYA(1:3*NATOMS)=COORDSA(1:3*NATOMS)
ENDIF

DO J1=1,NATOMS
   ALLPERM(J1)=J1
ENDDO

! The optimal alignment returned by minpermdist is a local minimum, but may not  
! be the global minimum. Calling MYORIENT first should put permutational isomers
! into a standard alignment and spot the global minimum zedro distance in one
! go. However, we also need to cycle over equivalent atoms in orbits using NCHOOSE2.
!
! Problems can occur if we don't use all the atoms specified by NORBIT1 and NORBIT2
! because of the numerical cutoffs employed in MYORIENT. We could miss the
! right orientation! 
!
! If we use MYORIENT to produce particular orientations then we end up aligning 
! COORDSA not with COORDSB but with the standard orientation of COORDSB in DUMMYB.
! We now deal with this by tracking the complete transformation, including the
! contribution of MYORIENT using ROTB and ROTINVB.
!
DISTANCE=0.0D0 
IF (NFREEZE.LE.0..AND. .NOT. NOTRANSROTT) THEN
   IF (BULKT) THEN
      IF (BMTEST) BMCOORDSSV(1:3*NATOMS)=BMCOORDS(1:3*NATOMS)
      NORBIT1=1; NORBIT2=1; NORBITB1=1; NORBITB2=1;
      CALL BULKMINDIST(DUMMYB,DUMMYA,BMCOORDS,NATOMS,DISTANCE,TWOD,DEBUG,BOXLX,BOXLY,BOXLZ,PITEST,.TRUE.,TNMATCH, BMTEST)
      IF (PITEST) THEN
         COORDSA(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
         RMATBEST(1:3,1:3)=0.0D0
         RMATBEST(1,1)=1.0D0; RMATBEST(2,2)=1.0D0; RMATBEST(3,3)=1.0D0
         DISTANCE=SQRT(DISTANCE)
         RETURN
       ELSE IF (BMTEST) THEN
          IF (SQRT(DISTANCE).LT.BMDIST) THEN
             BMDIST=SQRT(DISTANCE)
             IF (DEBUG) PRINT '(A,G20.10)',' minpermdist> after initial call to BULKMINDIST distance=', BMDIST
          ELSE
             BMCOORDS(1:3*NATOMS)=BMCOORDSSV(1:3*NATOMS)
             IF (DEBUG) PRINT '(A,G20.10)',' minpermdist> after initial call to BULKMINDIST distance has not improved'
          ENDIF
!         IF ((OHCELLT.AND.OPNUM.EQ.0).OR.(.NOT.OHCELLT)) THEN
!           CALL NEWMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,BULKT,TWOD,'AX',.FALSE.,RIGID,DEBUG,RMAT)
!         IF (BMDIST.GT.DISTANCE) THEN
!            IF (DEBUG) PRINT '(A,G20.10)',' minpermdist> shortest distance so far from NEWMINDIST distance=', DISTANCE 
!          BMTEST=.FALSE. 
!          OPNUM=48 
!         ENDIF
!         ENDIF
          IF (BMDIST.LT.DISTANCE) THEN         
             DISTANCE=BMDIST
             IF (DEBUG) PRINT '(A,G20.10)',' minpermdist> shortest distance so far from BULKMINDIST distance=', DISTANCE 
             IF (OHCELLT.AND.(OPNUM.LT.48)) THEN
              GOTO 25
             ELSE
              DUMMYA(1:3*NATOMS)=BMCOORDS(1:3*NATOMS)
              COORDSA(1:3*NATOMS)=BMCOORDS(1:3*NATOMS)
              DISTANCE=0.0D0
              DO J3=1,NATOMS
                  DISTANCE=DISTANCE + (COORDSA(3*(J3-1)+1)-DUMMYB(3*(J3-1)+1) &
  &                       -BOXLX*NINT((COORDSA(3*(J3-1)+1)-DUMMYB(3*(J3-1)+1))/BOXLX))**2 &
  &                                 + (COORDSA(3*(J3-1)+2)-DUMMYB(3*(J3-1)+2) &
  &                       -BOXLY*NINT((COORDSA(3*(J3-1)+2)-DUMMYB(3*(J3-1)+2))/BOXLY))**2
                  IF (.NOT.TWOD) DISTANCE=DISTANCE+(COORDSA(3*(J3-1)+3)-DUMMYB(3*(J3-1)+3)&
  &                       -BOXLZ*NINT((COORDSA(3*(J3-1)+3)-DUMMYB(3*(J3-1)+3))/BOXLZ))**2
              ENDDO
              DISTANCE=SQRT(DISTANCE)  
              IF (DEBUG) PRINT*, ' minpermdist> Recalculated distance=',DISTANCE
              RMATBEST(1:3,1:3)=0.0D0
              RMATBEST(1,1)=1.0D0; RMATBEST(2,2)=1.0D0; RMATBEST(3,3)=1.0D0
              RETURN  
             ENDIF 
          ENDIF 
       ELSE  
         CALL NEWMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,BULKT,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT)
       !  IF (DEBUG) PRINT '(A,G20.10)','minpermdist> after initial call to BULK/NEWMINDIST distance=',DISTANCE
         DISTANCE=DISTANCE**2 ! minperdist returns the distance squared for historical reasons
      ENDIF
   ELSEIF (MKTRAPT) THEN
      TMAT(1:3,1:3)=0.0D0
      TMAT(1,1)=INVERT*1.0D0; TMAT(2,2)=INVERT*1.0D0; TMAT(3,3)=INVERT*1.0D0
      NORBIT1=1; NORBIT2=1; NORBITB1=1; NORBITB2=1;
      ROTB(1:3,1:3)=0.0D0
      ROTB(1,1)=1.0D0; ROTB(2,2)=1.0D0; ROTB(3,3)=1.0D0
      ROTINVB(1:3,1:3)=0.0D0
      ROTINVB(1,1)=1.0D0; ROTINVB(2,2)=1.0D0; ROTINVB(3,3)=1.0D0
      ROTA(1:3,1:3)=0.0D0
      ROTA(1,1)=1.0D0; ROTA(2,2)=1.0D0; ROTA(3,3)=1.0D0
      ROTINVA(1:3,1:3)=0.0D0
      ROTINVA(1,1)=1.0D0; ROTINVA(2,2)=1.0D0; ROTINVA(3,3)=1.0D0
      RMAT(1:3,1:3)=0.0D0
      RMAT(1,1)=1.0D0; RMAT(2,2)=1.0D0; RMAT(3,3)=1.0D0
      CMX=0.0D0; CMY=0.0D0; CMZ=0.0D0
      DUMMYA(1:3*NATOMS)=INVERT*COORDSA(1:3*NATOMS)
      DISTANCE=0.0D0
      DO J1=1,3*NATOMS
         DISTANCE=DISTANCE+(DUMMYA(J1)-DUMMYB(J1))**2      
      ENDDO
   ELSE
      IF ((TWOD.OR.PULLT).AND.(INVERT.EQ.-1)) THEN ! reflect in xz plane
         DO J1=1,NATOMS
            DUMMYC(3*(J1-1)+1)=DUMMYA(3*(J1-1)+1)
            DUMMYC(3*(J1-1)+2)=-DUMMYA(3*(J1-1)+2)
            DUMMYC(3*(J1-1)+3)=DUMMYA(3*(J1-1)+3)
         ENDDO
      ELSE IF (MACROCYCLET) THEN !Permutations for macrocycles
      ! We can't use MACROCYCLE and PULL at the same time, but it's unlikely anyone will want to.
         DO J1=1,NATOMS
            J2=J1+MCYCLESHIFT
            IF (J2.GT.NATOMS) THEN
               J2=J2-NATOMS
            ENDIF
            DUMMYC(3*J1-2)=INVERT*DUMMYA(3*J2-2)
            DUMMYC(3*J1-1)=INVERT*DUMMYA(3*J2-1)
            DUMMYC(3*J1  )=INVERT*DUMMYA(3*J2  )
         ENDDO
      ELSE
         DUMMYC(1:3*NATOMS)=INVERT*DUMMYA(1:3*NATOMS)
      ENDIF
      IF ((NRANROT.GT.0).AND.(NROTDONE.LE.NRANROT).AND.(NROTDONE.GT.0).AND.(.NOT.MKTRAPT)) THEN
!        IF (DEBUG) PRINT '(A,I6,A,G20.10)',' minpermdist> Trying random starting orientation number ',NROTDONE, &
! &                                         ' minimum distance=',SQRT(DBEST)
         NORBIT1=1; NORBIT2=1; NORBITB1=1; NORBITB2=1;
         ROTB(1:3,1:3)=0.0D0
         ROTB(1,1)=1.0D0; ROTB(2,2)=1.0D0; ROTB(3,3)=1.0D0
         ROTINVB(1:3,1:3)=0.0D0
         ROTINVB(1,1)=1.0D0; ROTINVB(2,2)=1.0D0; ROTINVB(3,3)=1.0D0
         ROTA(1:3,1:3)=0.0D0
         ROTA(1,1)=1.0D0; ROTA(2,2)=1.0D0; ROTA(3,3)=1.0D0
         ROTINVA(1:3,1:3)=0.0D0
         ROTINVA(1,1)=1.0D0; ROTINVA(2,2)=1.0D0; ROTINVA(3,3)=1.0D0
         RMAT(1:3,1:3)=0.0D0
         RMAT(1,1)=1.0D0; RMAT(2,2)=1.0D0; RMAT(3,3)=1.0D0
         CMX=0.0D0; CMY=0.0D0; CMZ=0.0D0
         DO I=1,NATOMS
            CMX=CMX+DUMMYC(3*(I-1)+1)
            CMY=CMY+DUMMYC(3*(I-1)+2)
            CMZ=CMZ+DUMMYC(3*(I-1)+3)
         ENDDO
         CMX=CMX/NATOMS; CMY=CMY/NATOMS; CMZ=CMZ/NATOMS
         DO I=1,NATOMS
            DUMMYC(3*(I-1)+1)=DUMMYC(3*(I-1)+1)-CMX
            DUMMYC(3*(I-1)+2)=DUMMYC(3*(I-1)+2)-CMY
            DUMMYC(3*(I-1)+3)=DUMMYC(3*(I-1)+3)-CMZ
         ENDDO
         CMX=0.0D0; CMY=0.0D0; CMZ=0.0D0
         DO I=1,NATOMS
            CMX=CMX+DUMMYB(3*(I-1)+1)
            CMY=CMY+DUMMYB(3*(I-1)+2)
            CMZ=CMZ+DUMMYB(3*(I-1)+3)
         ENDDO
         CMX=CMX/NATOMS; CMY=CMY/NATOMS; CMZ=CMZ/NATOMS
         DO I=1,NATOMS
            DUMMYB(3*(I-1)+1)=DUMMYB(3*(I-1)+1)-CMX
            DUMMYB(3*(I-1)+2)=DUMMYB(3*(I-1)+2)-CMY
            DUMMYB(3*(I-1)+3)=DUMMYB(3*(I-1)+3)-CMZ
         ENDDO
         CALL RANROT(DUMMYC,ROTA,ROTINVA,NATOMS)
         DUMMYA(1:3*NATOMS)=DUMMYC(1:3*NATOMS)
      ELSE
         CALL MYORIENT(DUMMYC,DUMMY,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,NATOMS,DEBUG,ROTA,ROTINVA,PULLT)
         DUMMYA(1:3*NATOMS)=DUMMY(1:3*NATOMS)
         CALL MYORIENT(DUMMYB,DUMMY,NORBITB1,NCHOOSEB1,NORBITB2,NCHOOSEB2,NATOMS,DEBUG,ROTB,ROTINVB,PULLT)
         DUMMYB(1:3*NATOMS)=DUMMY(1:3*NATOMS)
      ENDIF
      DISTANCE=0.0D0
      DO J1=1,3*NATOMS
         DISTANCE=DISTANCE+(DUMMYA(J1)-DUMMYB(J1))**2
      ENDDO
   ENDIF
!  IF (DEBUG) PRINT '(A,G20.10,A,I6,A)', &
! &       ' minpermdist> after initial call to MYORIENT distance=',SQRT(DISTANCE), ' for ',NATOMS,' atoms'
!  IF (DEBUG) PRINT '(A,6I8)',' minpermdist> size of orbits, selected atoms, random rotations, invert: ', &
! &       NORBIT1,NORBIT2,NCHOOSE1,NCHOOSE2,NROTDONE,INVERT
ELSE
  NORBIT1=1; NORBIT2=1; NORBITB1=1; NORBITB2=1
  CALL NEWMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,BULKT,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT)
  IF (DEBUG) PRINT '(A,G20.10)','minpermdist> after initial call to NEWMINDIST distance=',DISTANCE
  DISTANCE=DISTANCE**2
ENDIF
!
!  Bipartite matching routine for permutations. Coordinates in DUMMYB do not change
!  but the coordinates in DUMMYA do. DISTANCE is the distance^2 in this case.
!  We return to label 10 after every round of permutational/orientational alignment
!  unless we have converged to the identity permutation.
!
!  Atoms are not allowed to appear in more than one group.
!  The maximum number of pair exchanges associated with a group is two.
!
NTRIES=0
!
!  RMATCUMUL contains the accumulated rotation matrix that relates the original 
!  DUMMYA obtained from COORDSA to the final one.
!
RMATCUMUL(1:3,1:3)=0.0D0
RMATCUMUL(1,1)=1.0D0; RMATCUMUL(2,2)=1.0D0; RMATCUMUL(3,3)=1.0D0
10 CONTINUE
NTRIES=NTRIES+1

NDUMMY=1
DO J1=1,NATOMS
   NEWPERM(J1)=J1
ENDDO
!
! ALLPERM saves the permutation from the previous cycle.
! NEWPERM contains the permutation for this cycle, relative to the identity.
! SAVEPERM is temporary storage for NEWPERM.
! NEWPERM must be applied to ALLPERM after the loop over NPERMGROUP and
! corresponding swaps.
!
! New version allows for overlapping atoms in NPERMGROUP, so that atoms
! can appear in more thsan one group. This was needed to flexible water potentials.
!

DO J1=1,NPERMGROUP
   PATOMS=NPERMSIZE(J1)
   DO J2=1,PATOMS
      PDUMMYA(3*(J2-1)+1)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+1)
      PDUMMYA(3*(J2-1)+2)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+2)
      PDUMMYA(3*(J2-1)+3)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+3)
      PDUMMYB(3*(J2-1)+1)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+1)
      PDUMMYB(3*(J2-1)+2)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+2)
      PDUMMYB(3*(J2-1)+3)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+3)
   ENDDO
   CALL MINPERM(PATOMS, PDUMMYB, PDUMMYA, BOXLX, BOXLY, BOXLZ, BULKT, LPERM, LDISTANCE, DIST2, WORSTRAD)
   SAVEPERM(1:NATOMS)=NEWPERM(1:NATOMS)
   DO J2=1,PATOMS
      SAVEPERM(PERMGROUP(NDUMMY+J2-1))=NEWPERM(PERMGROUP(NDUMMY+LPERM(J2)-1))
   ENDDO
!
! Update permutation of associated atoms, if any.
! We must do this as we go along, because these atoms could move in more than
! one permutational group now.
!
   IF (NSETS(J1).GT.0) THEN
      DO J2=1,PATOMS
         DO J3=1,NSETS(J1)
            SAVEPERM(SETS(PERMGROUP(NDUMMY+J2-1),J3))=SETS(NEWPERM(PERMGROUP(NDUMMY+LPERM(J2)-1)),J3)
         ENDDO
      ENDDO
   ENDIF
   NDUMMY=NDUMMY+NPERMSIZE(J1)
   NEWPERM(1:NATOMS)=SAVEPERM(1:NATOMS)
ENDDO
!
! Update the overall permutation here.
! The latest NEWPERM(J1) tells us which position moves to J1 in the latest
! permutation, relative to the identity.
! ALLPERM(J2) tells us which atom has moved to position J2.
! So, the new overall permutation, i.e. the atoms that moves to position J1
! After ALLPERM followed by NEWPERM is ALLPERM(NEWPERM(J1))
! 
DO J1=1,NATOMS
   SAVEPERM(J1)=ALLPERM(NEWPERM(J1))
ENDDO
ALLPERM(1:NATOMS)=SAVEPERM(1:NATOMS)

DUMMY(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
NPERM=0
DISTANCE=0.0D0
IF (RBAAT) THEN ! additional permutation of angle-axis coordinates
   DO J1=(NATOMS/2)+1,NATOMS
      ALLPERM(J1)=ALLPERM(J1-(NATOMS/2))+(NATOMS/2)
      NEWPERM(J1)=NEWPERM(J1-(NATOMS/2))+(NATOMS/2)
   ENDDO
ENDIF
!
! Update coordinates in DUMMYA to overall permutation using NEWPERM.
!
DO J3=1,NATOMS
   DUMMYA(3*(J3-1)+1)=DUMMY(3*(NEWPERM(J3)-1)+1)
   DUMMYA(3*(J3-1)+2)=DUMMY(3*(NEWPERM(J3)-1)+2)
   DUMMYA(3*(J3-1)+3)=DUMMY(3*(NEWPERM(J3)-1)+3)
   IF (J3.NE.NEWPERM(J3)) THEN
!     IF (DEBUG) WRITE(*,'(A,I5,A,I5)') ' minpermdist> move position ',NEWPERM(J3),' to ',J3
      NPERM=NPERM+1
   ENDIF
   IF (.NOT.BULKT) THEN
      DISTANCE=DISTANCE+(DUMMYA(3*(J3-1)+1)-DUMMYB(3*(J3-1)+1))**2 &
  &                    +(DUMMYA(3*(J3-1)+2)-DUMMYB(3*(J3-1)+2))**2 &
  &                    +(DUMMYA(3*(J3-1)+3)-DUMMYB(3*(J3-1)+3))**2
   ELSE
      DISTANCE=DISTANCE + (DUMMYA(3*(J3-1)+1)-DUMMYB(3*(J3-1)+1)- BOXLX*NINT((DUMMYA(3*(J3-1)+1)-DUMMYB(3*(J3-1)+1))/BOXLX))**2 &
  &                     + (DUMMYA(3*(J3-1)+2)-DUMMYB(3*(J3-1)+2)- BOXLY*NINT((DUMMYA(3*(J3-1)+2)-DUMMYB(3*(J3-1)+2))/BOXLY))**2 
      IF (.NOT.TWOD) DISTANCE=DISTANCE &
  &                     + (DUMMYA(3*(J3-1)+3)-DUMMYB(3*(J3-1)+3)- BOXLZ*NINT((DUMMYA(3*(J3-1)+3)-DUMMYB(3*(J3-1)+3))/BOXLZ))**2
   ENDIF
ENDDO

! IF (DEBUG) WRITE(*,'(A,I6,A,G20.10)') 'minpermdist> distance after permuting ',NPERM,' pairs of atoms=',SQRT(DISTANCE)

! CALL OCHARMM(DUMMYA,VNEW,ENERGY,.FALSE.,.FALSE.)
! PRINT '(A,F25.15,A)',' Energy for last cycle=',ENERGY,' kcal/mol'
! CALL UPDATENBONDS(DUMMYA)
! PRINT '(A,F25.15,A)',' Energy for last cycle=',ENERGY,' kcal/mol after update'
! CALL POTENTIAL(DUMMYA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
! PRINT '(2(A,F25.15))',' Energy for last cycle=',ENERGY,' RMS=',RMS
!
!  Optimal alignment. Coordinates in DUMMYA are reset by NEWMINDIST (second argument).
!  Must allow at least one call to NEWMINDIST in case the MYORIENT result is terrible
!  but gives zero permutations!
!  
IF ((NPERM.NE.0).OR.(NTRIES.EQ.1)) THEN 
   CALL NEWMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,BULKT,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT)
   RMATCUMUL=MATMUL(RMAT,RMATCUMUL)
   DISTANCE=DISTANCE**2 ! we are using DISTANCE^2 further down
!  IF (DEBUG) WRITE(*,'(A,G20.10)') ' minpermdist> distance after NEWMINDIST=                     ', &
! &                                    SQRT(DISTANCE) 
   IF (NTRIES.LT.MAXIMUMTRIES) THEN
      GOTO 10
   ELSE ! prevent infinite loop
      IF (DEBUG) PRINT '(A)','minpermdist> WARNING - number of tries exceeded, giving up'
   ENDIF
ENDIF

IF (DISTANCE.LT.DBEST) THEN
   DBEST=DISTANCE
   XBEST(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
   BESTPERM(1:NATOMS)=ALLPERM(1:NATOMS)
   BESTINVERT=INVERT
   RMATBEST(1:3,1:3)=RMATCUMUL(1:3,1:3)
   ROTINVBBEST(1:3,1:3)=ROTINVB(1:3,1:3) 
   ROTABEST(1:3,1:3)=ROTA(1:3,1:3)      
   RMATBEST=MATMUL(RMATBEST,ROTABEST)
   IF (INVERT.EQ.-1) THEN
      IF (PULLT.OR.TWOD) THEN ! reflect in xz plane rather than invert!
         RMATBEST(1:3,1:3)=MATMUL(RMATBEST,REFXZ)
      ELSE
         RMATBEST(1:3,1:3)=-RMATBEST(1:3,1:3)
      ENDIF
   ENDIF
ENDIF
!
! If GEOMDIFFTOL is set too large we could miss the best solution by exiting prematurely.
! Turn off the next line?!
! This greatly increases the time taken to parse path.info files - not recommended!
! 
IF (SQRT(DBEST).LT.GEOMDIFFTOL) GOTO 50
!
!IF(BULKT) GOTO 50 ! sn402

IF (OHCELLT.AND.(OPNUM.LT.48)) GOTO 25
IF (NCHOOSE2.LT.NORBIT2) GOTO 30
IF (NCHOOSE1.LT.NORBIT1) GOTO 65
IF (NCHOOSEB2.LT.NORBITB2) GOTO 31
IF (NCHOOSEB1.LT.NORBITB1) GOTO 66
!
!  Now try the enantiomer (or xz reflected structure for PULLT.OR.EFIELDT.OR.TWOD).
!  The tests for NCHOOSE1 and NCHOOSE2 appear to be redundant!
!
IF ((NCHOOSE2.EQ.NORBIT2).AND.(NCHOOSE1.EQ.NORBIT1).AND.(INVERT.EQ.1)) THEN
!
! don't try inversion for bulk or charmm or amber or frozen atoms
! Inversion is allowed for all macrocycles in AMBER/CHARMM
!
!  IF ((BULKT.OR.CHARMMT.OR.AMBERT.OR.(NFREEZE.GT.0)).AND.MACROCYCLET.EQV..FALSE.) GOTO 50 
! for the macroion model, allow inversion (we use AMBER which does not allow it by default)
   IF (MACROIONT.AND.NFREEZE.EQ.0) THEN
    INVERT=-1
    GOTO 60
   END IF
!
!  Bug fix DJW 12/7/12. When MACROCYCLET was false all inversion tests were turned off!!
!
   IF (NOINVERSION.OR.BULKT.OR.(CHARMMT.AND.(.NOT.MACROCYCLET)).OR.(AMBERT.AND.(.NOT.MACROCYCLET)) &
  &     .OR.(NFREEZE.GT.0).OR.NOTRANSROTT) GOTO 50 
!  IF (DEBUG) PRINT '(A)','minpermdist> inverting geometry for comparison with target'
   INVERT=-1
   GOTO 60
ENDIF
MCYCLESTEP=MCYCLESTEP+1 
IF (MACROCYCLET.AND.(MCYCLESTEP.LE.MCYCLEREPEATS)) GOTO 70 !End of macrocycle loop
IF (NROTDONE.LT.NRANROT) GOTO 11

50 DISTANCE=DBEST
!
!  XBEST contains the best alignment of A coordinates for the orientation of B coordinates in DUMMYB.
!  Rotate XBEST by ROTINVBBEST to put in best correspondence with COORDSB, 
!  undoing the reorientation to DUMMYB from MYORIENT. 
!  We should get the same result for ROTINVBBEST * RMATBEST * (COORDSA-CMA)
!  where RMATBEST = +/- RMATCUMUL * ROTA for the best alignment 
!  (aside from a possible permutation of the atom ordering)
!  

   IF (BULKT) THEN
      XDUMMY=0.0D0
      DO J1=1,NATOMS
         XDUMMY=XDUMMY+(COORDSB(3*(J1-1)+1)-XBEST(3*(J1-1)+1) - BOXLX*NINT((COORDSB(3*(J1-1)+1)-XBEST(3*(J1-1)+1))/BOXLX))**2+ &
  &                    (COORDSB(3*(J1-1)+2)-XBEST(3*(J1-1)+2) - BOXLY*NINT((COORDSB(3*(J1-1)+2)-XBEST(3*(J1-1)+2))/BOXLY))**2
         IF (.NOT.TWOD) XDUMMY=XDUMMY &
  &                   +(COORDSB(3*(J1-1)+3)-XBEST(3*(J1-1)+3) - BOXLZ*NINT((COORDSB(3*(J1-1)+3)-XBEST(3*(J1-1)+3))/BOXLZ))**2
      ENDDO
   ELSEIF (NFREEZE.GT.0.OR.RBAAT.OR.MKTRAPT .OR. NOTRANSROTT) THEN
      XDUMMY=0.0D0
      DO J1=1,NATOMS
         XDUMMY=XDUMMY+(COORDSB(3*(J1-1)+1)-XBEST(3*(J1-1)+1))**2+ &
  &                    (COORDSB(3*(J1-1)+2)-XBEST(3*(J1-1)+2))**2+ &
  &                    (COORDSB(3*(J1-1)+3)-XBEST(3*(J1-1)+3))**2
      ENDDO
   ELSE
      XDUMMY=0.0D0
      DO J1=1,NATOMS
         XBEST(3*(J1-1)+1:3*(J1-1)+3)=MATMUL(ROTINVBBEST,XBEST(3*(J1-1)+1:3*(J1-1)+3))
         XBEST(3*(J1-1)+1)=XBEST(3*(J1-1)+1)+CMBX
         XBEST(3*(J1-1)+2)=XBEST(3*(J1-1)+2)+CMBY
         XBEST(3*(J1-1)+3)=XBEST(3*(J1-1)+3)+CMBZ
         XDUMMY=XDUMMY+(COORDSB(3*(J1-1)+1)-XBEST(3*(J1-1)+1))**2+ &
  &                    (COORDSB(3*(J1-1)+2)-XBEST(3*(J1-1)+2))**2+ &
  &                    (COORDSB(3*(J1-1)+3)-XBEST(3*(J1-1)+3))**2
      ENDDO
   ENDIF
   IF (ABS(SQRT(XDUMMY)-SQRT(DISTANCE)).GT.GEOMDIFFTOL .AND. (.NOT. RBAAT)) THEN
      PRINT '(2(A,G20.10))','minpermdist> ERROR *** distance between transformed XBEST and COORDSB=',SQRT(XDUMMY), &
  &                         ' should be ',SQRT(DISTANCE)
      PRINT '(A)','transformed XBEST:'
      PRINT '(3F20.10)',XBEST(1:3*NATOMS)
      PRINT '(A)','COORDSB:'
      PRINT '(3F20.10)',COORDSB(1:3*NATOMS)
      STOP
   ENDIF

   IF ((NFREEZE.GT.0).OR.MKTRAPT.OR.NOTRANSROTT) THEN
!no rotation for NFREEZE .gt. 0
      RMATBEST(1:3,1:3)=0.0D0
      RMATBEST(1,1)=1.0D0; RMATBEST(2,2)=1.0D0; RMATBEST(3,3)=1.0D0
   ELSE
      RMATBEST=MATMUL(ROTINVBBEST,RMATBEST)
   ENDIF
!
! For TRAP potentials we need to preserve the centre of mass because
! the origin is a fixed point. Try using RMAT best and rotating COORDSA
! around the origin.
!
   IF (TRAPT) THEN
      XDUMMY=0.0D0
      DO J3=1,NATOMS
         DUMMYC(3*(J3-1)+1)=COORDSA(3*(BESTPERM(J3)-1)+1)
         DUMMYC(3*(J3-1)+2)=COORDSA(3*(BESTPERM(J3)-1)+2)
         DUMMYC(3*(J3-1)+3)=COORDSA(3*(BESTPERM(J3)-1)+3)
         XDUMMY=XDUMMY+(COORDSB(3*(J3-1)+1)-DUMMYC(3*(J3-1)+1))**2+ &
  &                    (COORDSB(3*(J3-1)+2)-DUMMYC(3*(J3-1)+2))**2+ &
  &                    (COORDSB(3*(J3-1)+3)-DUMMYC(3*(J3-1)+3))**2
      ENDDO
!     PRINT '(A,G20.10)',' minpermdist> with permutations distance is ',SQRT(XDUMMY)
      CALL NEWROTGEOM(NATOMS,DUMMYC,RMATBEST,0.0D0,0.0D0,0.0D0)
!     XDUMMY=0.0D0
!     DO J3=1,NATOMS
!        XDUMMY=XDUMMY+(COORDSB(3*(J3-1)+1)-DUMMYC(3*(J3-1)+1))**2+ &
! &                    (COORDSB(3*(J3-1)+2)-DUMMYC(3*(J3-1)+2))**2+ &
! &                    (COORDSB(3*(J3-1)+3)-DUMMYC(3*(J3-1)+3))**2
!     ENDDO
!     PRINT '(A,G20.10)',' minpermdist> with permutations after rotation distance is ',SQRT(XDUMMY)
      XBEST(1:3*NATOMS)=DUMMYC(1:3*NATOMS)
   ENDIF

   COORDSA(1:3*NATOMS)=XBEST(1:3*NATOMS) ! finally, best COORDSA should include permutations for DNEB input!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  DO J1=1,(NATOMS/2)
!     XDUMMY=XDUMMY+(COORDSB(3*(J1-1)+1)-COORDSA(3*(J1-1)+1))**2+ &
! &                 (COORDSB(3*(J1-1)+2)-COORDSA(3*(J1-1)+2))**2+ &
! &                 (COORDSB(3*(J1-1)+3)-COORDSA(3*(J1-1)+3))**2
!  ENDDO
!  PRINT '(A,F20.10)','XDUMMY=',XDUMMY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CLOSE(10)

DISTANCE=SQRT(DISTANCE)

RETURN
END SUBROUTINE MINPERMDIST

SUBROUTINE RANROT(COORDS,ROT,ROTINV,NATOMS)
! rewritten by js850 to use the new unbiased rotations
!
! subtract the center of mass from coords
! rotate coords by a uniform random rotation.
! return the rotation ROT and the inverse ROTINV
USE COMMONS, ONLY : GTHOMSONT, GTHOMMET
USE ROTATIONS, ONLY : ROT_RANDOM_Q, ROT_Q2MX
IMPLICIT NONE
INTEGER, INTENT(IN) :: NATOMS
DOUBLE PRECISION, INTENT(INOUT) :: COORDS(3*NATOMS)
DOUBLE PRECISION, INTENT(OUT) :: ROTINV(3,3), ROT(3,3)
INTEGER J1
DOUBLE PRECISION ANGLE, COST, SINT, RANDOM, PI, TX, TY, DPRAND
DOUBLE PRECISION CMX, CMY, CMZ
DOUBLE PRECISION T1(3*NATOMS)
DOUBLE PRECISION qrot(4)

! hk286
IF (GTHOMSONT .AND. (GTHOMMET < 5) ) THEN
   T1(1:3*NATOMS)=COORDS(1:3*NATOMS)

   PI=3.14159D0
   RANDOM=(DPRAND()-0.5D0)*2.0D0
   ANGLE=RANDOM*PI ! in radians
   COST=COS(ANGLE)
   SINT=SIN(ANGLE)
   DO J1=1,NATOMS
      TX= COST*T1(3*(J1-1)+1)+SINT*T1(3*(J1-1)+2)
      TY=-SINT*T1(3*(J1-1)+1)+COST*T1(3*(J1-1)+2)
      T1(3*(J1-1)+1)=TX
      T1(3*(J1-1)+2)=TY
   ENDDO
   
   ROT(1,1)=COST;  ROT(1,2)=-SINT; ROT(1,3)=0.0D0
   ROT(2,1)=SINT;  ROT(2,2)=COST;  ROT(2,3)=0.0D0
   ROT(3,1)=0.0D0; ROT(3,2)=0.0D0; ROT(3,3)=1.0D0

   rotinv(1,1:3) = ROT(1:3,1); rotinv(2,1:3)=ROT(1:3,2); rotinv(3,1:3)=ROT(1:3,3) ! transpose

ELSEIF (GTHOMSONT .AND. (GTHOMMET .EQ. 5) ) THEN
   
   !get the quaternion giving the random rotation in 3 dimensions
   qrot = rot_random_q()
   !convert it to matrix format
   rot = rot_q2mx(qrot)
   !get the transpose
   rotinv(1,1:3) = ROT(1:3,1); rotinv(2,1:3)=ROT(1:3,2); rotinv(3,1:3)=ROT(1:3,3) ! transpose
   !rotate coords by rot
   do j1=1,natoms
      T1(3*(J1-1)+1:3*(J1-1)+3)=MATMUL(ROT,coords(3*(J1-1)+1:3*(J1-1)+3))
   enddo

ELSE

!subtract the center of mass from the coordinates
   CMX = sum(COORDS(1:3*NATOMS:3)) / NATOMS
   CMY = sum(COORDS(2:3*NATOMS:3)) / NATOMS
   CMZ = sum(COORDS(3:3*NATOMS:3)) / NATOMS
   COORDS(1:3*NATOMS:3) = COORDS(1:3*NATOMS:3) - CMX
   COORDS(2:3*NATOMS:3) = COORDS(2:3*NATOMS:3) - CMY
   COORDS(3:3*NATOMS:3) = COORDS(3:3*NATOMS:3) - CMZ
   
!get the quaternion giving the random rotation in 3 dimensions
   qrot = rot_random_q()
!convert it to matrix format
   rot = rot_q2mx(qrot)
!get the transpose
   rotinv(1,1:3) = ROT(1:3,1); rotinv(2,1:3)=ROT(1:3,2); rotinv(3,1:3)=ROT(1:3,3) ! transpose
!rotate coords by rot
   do j1=1,natoms
      T1(3*(J1-1)+1:3*(J1-1)+3)=MATMUL(ROT,coords(3*(J1-1)+1:3*(J1-1)+3))
   enddo

endif

!copy T1 into coords
COORDS(1:3*NATOMS)=T1(1:3*NATOMS)

END SUBROUTINE RANROT


!SUBROUTINE RANROT(COORDS,ROT,ROTINV,NATOMS)
!IMPLICIT NONE
!INTEGER J1, NATOMS, I, J2, J3
!DOUBLE PRECISION ANGLE, COST, SINT, RANDOM, PI, TX, TY, TZ, COORDS(3*NATOMS), DPRAND, TEMP(3,3)
!DOUBLE PRECISION CMX, CMY, CMZ, ROT(3,3), ROTINV(3,3), TVEC(3*NATOMS)
!DOUBLE PRECISION T1(3*NATOMS), DUMMY, R1(3,3), R2(3,3)
!
!CMX=0.0D0
!CMY=0.0D0
!CMZ=0.0D0
!DO I=1,NATOMS
!   CMX=CMX+COORDS(3*(I-1)+1)
!   CMY=CMY+COORDS(3*(I-1)+2)
!   CMZ=CMZ+COORDS(3*(I-1)+3)
!ENDDO
!CMX=CMX/NATOMS
!CMY=CMY/NATOMS
!CMZ=CMZ/NATOMS
!DO I=1,NATOMS
!   COORDS(3*(I-1)+1)=COORDS(3*(I-1)+1)-CMX
!   COORDS(3*(I-1)+2)=COORDS(3*(I-1)+2)-CMY
!   COORDS(3*(I-1)+3)=COORDS(3*(I-1)+3)-CMZ
!ENDDO
!T1(1:3*NATOMS)=COORDS(1:3*NATOMS)
!
!PI=3.14159D0
!
!RANDOM=(DPRAND()-0.5D0)*PI
!ANGLE=RANDOM ! in radians
!COST=COS(ANGLE)
!SINT=SIN(ANGLE)
!!
!! Random rotation about the x axis.
!!
!DO J1=1,NATOMS
!   TY= COST*T1(3*(J1-1)+2)+SINT*T1(3*(J1-1)+3)
!   TZ=-SINT*T1(3*(J1-1)+2)+COST*T1(3*(J1-1)+3)
!   T1(3*(J1-1)+2)=TY
!   T1(3*(J1-1)+3)=TZ
!ENDDO
!
!R1(1,1)=1.0D0; R1(1,2)=0.0D0; R1(1,3)=0.0D0
!R1(2,1)=0.0D0; R1(2,2)=COST;  R1(2,3)=SINT
!R1(3,1)=0.0D0; R1(3,2)=-SINT; R1(3,3)=COST
!
!RANDOM=(DPRAND()-0.5D0)*PI
!ANGLE=RANDOM*PI ! in radians
!COST=COS(ANGLE)
!SINT=SIN(ANGLE)
!!
!! Random rotation about the y axis.
!!
!DO J1=1,NATOMS
!   TX= COST*T1(3*(J1-1)+1)+SINT*T1(3*(J1-1)+3)
!   TZ=-SINT*T1(3*(J1-1)+1)+COST*T1(3*(J1-1)+3)
!   T1(3*(J1-1)+1)=TX
!   T1(3*(J1-1)+3)=TZ
!ENDDO
!
!R2(1,1)=COST; R2(1,2)=0.0D0; R2(1,3)=-SINT
!R2(2,1)=0.0D0; R2(2,2)=1.0D0; R2(2,3)=0.0D0
!R2(3,1)=SINT; R2(3,2)=0.0D0; R2(3,3)=COST
!!
!! MATMULV returns A^T= B^T C. Hence R2 above needs to be the
!! transpose.
!!
!CALL MATMULV(ROT,R2,R1,3,3,3)
!R1(1,1:3)=ROT(1:3,1); R1(2,1:3)=ROT(1:3,2); R1(3,1:3)=ROT(1:3,3)
!
!!
!! Random rotation about the z axis.
!!
!RANDOM=(DPRAND()-0.5D0)*2.0D0
!ANGLE=RANDOM*PI ! in radians
!COST=COS(ANGLE)
!SINT=SIN(ANGLE)
!DO J1=1,NATOMS
!   TX= COST*T1(3*(J1-1)+1)+SINT*T1(3*(J1-1)+2)
!   TY=-SINT*T1(3*(J1-1)+1)+COST*T1(3*(J1-1)+2)
!   T1(3*(J1-1)+1)=TX
!   T1(3*(J1-1)+2)=TY
!ENDDO
!
!R2(1,1)=COST; R2(1,2)=-SINT; R2(1,3)=0.0D0
!R2(2,1)=SINT; R2(2,2)=COST; R2(2,3)=0.0D0
!R2(3,1)=0.0D0; R2(3,2)=0.0D0; R2(3,3)=1.0D0
!
!CALL MATMULV(ROTINV,R2,R1,3,3,3)
!!
!! MATMULV returns A^T= B^T C. Hence R2 above needs to be the
!! transpose.
!!
!
!ROT(1,1:3)=ROTINV(1:3,1); ROT(2,1:3)=ROTINV(1:3,2); ROT(3,1:3)=ROTINV(1:3,3) ! transpose
!
!
!COORDS(1:3*NATOMS)=T1(1:3*NATOMS)
!
!END SUBROUTINE RANROT

SUBROUTINE MATMULV(A,B,C,NA,NB,NC)
IMPLICIT NONE
INTEGER, INTENT(IN) :: NA, NB, NC
DOUBLE PRECISION, INTENT(IN) :: B(NB,NA), C(NB,NC)
DOUBLE PRECISION, INTENT(OUT) :: A(NC,NA)
INTEGER I,K,J
DOUBLE PRECISION Z

DO I=1,NA
   DO K=1,NC
      Z=0.D0
      DO J=1,NB
         Z=Z+B(J,I)*C(J,K)
      ENDDO
      A(K,I)=Z
   ENDDO
ENDDO

RETURN
END SUBROUTINE MATMULV


