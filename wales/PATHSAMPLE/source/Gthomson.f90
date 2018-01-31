SUBROUTINE HKMINPERMDIST(COORDSB,COORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGID,RMATBEST,USEINT)
USE COMMONS,ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, SETS, GEOMDIFFTOL, AMBERT, NFREEZE, CHARMMT, RBAAT, PULLT, &
  &               ANGLEAXIS, PERMISOMER, PERMDIST, ZSYM, INTCONSTRAINTT, INTLJT, OHCELLT, ATOMMATCHDIST, LPERMDIST, &
  &               TRAPT, NRANROT, MACROCYCLET, MCYCLEPERIOD, MCYCLEREPEATS, NOINVERSION, GTHOMSONT, GTHOMMET
USE PORFUNCS 
USE UTILS,ONLY : GETUNIT

IMPLICIT NONE

INTEGER :: MAXIMUMTRIES=100
INTEGER NATOMS, NPERM, PATOMS, NTRIES, I, LUNIT
INTEGER INVERT, NORBIT1, NORBIT2, NCHOOSE2, NDUMMY, LPERM(NATOMS), J1, J2, NCHOOSE1, OPNUM, J3, NROTDONE
INTEGER NORBITB1, NORBITB2, NCHOOSEB1, NCHOOSEB2

DOUBLE PRECISION DIST2, COORDSA(3*NATOMS), COORDSB(3*NATOMS), DISTANCE, DUMMYA(3*NATOMS), DUMMYB(3*NATOMS), DUMMY(3*NATOMS)
DOUBLE PRECISION BOXLX,BOXLY,BOXLZ,WORSTRAD,RMAT(3,3),ENERGY, VNEW(3*NATOMS), DX, DY, DZ, RMS, DBEST, XBEST(3*NATOMS)
DOUBLE PRECISION CMXA, CMXB, CMXC, CMX, CMY, CMZ
DOUBLE PRECISION ROTA(3,3), ROTINVA(3,3), ROTB(3,3), ROTINVB(3,3), ROTINVBBEST(3,3), ROTABEST(3,3), RMATBEST(3,3), TMAT(3,3)
DOUBLE PRECISION RMATCUMUL(3,3)
DOUBLE PRECISION REFXZ(3,3), QBEST(4)
DOUBLE PRECISION BMDIST
DOUBLE PRECISION, ALLOCATABLE :: BMCOORDS(:), BMCOORDSSV(:)
LOGICAL DEBUG, TWOD, RIGID, BULKT, PITEST, BMTEST, TNMATCH
LOGICAL, INTENT(IN) :: USEINT
DOUBLE PRECISION PDUMMYA(3*NATOMS), PDUMMYB(3*NATOMS), LDISTANCE, DUMMYC(3*NATOMS), XDUMMY
INTEGER NEWPERM(NATOMS), ALLPERM(NATOMS), SAVEPERM(NATOMS), NMAXINT, NMININT, BESTPERM(NATOMS)
DOUBLE PRECISION CONSTRAINTE, XYZLOCAL(6*NATOMS), LMINCOORDS(2,3*NATOMS)
INTEGER :: MCYCLESTEP =1,MCYCLESHIFT

NROTDONE=-1
MAXIMUMTRIES=MAX(MAXIMUMTRIES,NRANROT+1)

REFXZ(1:3,1:3)=0.0D0
REFXZ(1,1)=1.0D0; REFXZ(2,2)=-1.0D0; REFXZ(3,3)=1.0D0
!
! It is possible for the standard orientation to result in a distance that is worse than
! the starting distance. Hence we need to set XBEST here.
!
DUMMYA(1:3*NATOMS) = COORDSA(1:3*NATOMS)
DUMMYB(1:3*NATOMS) = COORDSB(1:3*NATOMS)
DBEST = 1.0D100

XDUMMY = 0.0D0
XBEST(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
DO J1 = 1,NATOMS
   XDUMMY = XDUMMY+(COORDSB(3*(J1-1)+1)-XBEST(3*(J1-1)+1))**2+ &
&                  (COORDSB(3*(J1-1)+2)-XBEST(3*(J1-1)+2))**2+ &
&                  (COORDSB(3*(J1-1)+3)-XBEST(3*(J1-1)+3))**2
ENDDO
DBEST = XDUMMY

DO J1=1,NATOMS
   BESTPERM(1:NATOMS)=J1
ENDDO
RMATBEST(1:3,1:3)=RMAT(1:3,1:3)
ROTINVBBEST(1:3,1:3)=0.0D0
ROTINVBBEST(1,1)=1.0D0;ROTINVBBEST(2,2)=1.0D0;ROTINVBBEST(3,3)=1.0D0;
ROTABEST(1:3,1:3)=0.0D0
ROTABEST(1,1)=1.0D0;ROTABEST(2,2)=1.0D0;ROTABEST(3,3)=1.0D0;

NROTDONE=-1
11 CONTINUE
NROTDONE=NROTDONE+1

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
DUMMYB(1:3*NATOMS)=COORDSB(1:3*NATOMS)
DUMMYA(1:3*NATOMS)=COORDSA(1:3*NATOMS)

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
IF ((NFREEZE.LE.0).AND.(.NOT.RBAAT)) THEN
   IF ((TWOD.OR.PULLT.OR.(GTHOMSONT .AND.(GTHOMMET < 5))).AND.(INVERT.EQ.-1)) THEN ! reflect in xz plane
      DO J1=1,NATOMS
         DUMMYC(3*(J1-1)+1)=DUMMYA(3*(J1-1)+1)
         DUMMYC(3*(J1-1)+2)=-DUMMYA(3*(J1-1)+2)
         DUMMYC(3*(J1-1)+3)=DUMMYA(3*(J1-1)+3)
      ENDDO
   ELSE
      DUMMYC(1:3*NATOMS)=INVERT*DUMMYA(1:3*NATOMS)
   ENDIF
   IF ((NRANROT.GT.0).AND.(NROTDONE.LE.NRANROT).AND.(NROTDONE.GT.0)) THEN
!     IF (DEBUG) PRINT '(A,I6,A,G20.10)',' minpermdist> Trying random starting orientation number ',NROTDONE, &
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
DO J1=1,NATOMS
!  SAVEPERM(ALLPERM(J1))=ALLPERM(NEWPERM(J1))
   SAVEPERM(J1)=ALLPERM(NEWPERM(J1))
ENDDO
ALLPERM(1:NATOMS)=SAVEPERM(1:NATOMS)

DUMMY(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
NPERM=0
DISTANCE=0.0D0

DO J3=1,NATOMS
   DUMMYA(3*(J3-1)+1)=DUMMY(3*(NEWPERM(J3)-1)+1)
   DUMMYA(3*(J3-1)+2)=DUMMY(3*(NEWPERM(J3)-1)+2)
   DUMMYA(3*(J3-1)+3)=DUMMY(3*(NEWPERM(J3)-1)+3)
   IF (J3.NE.NEWPERM(J3)) THEN
!     IF (DEBUG) WRITE(*,'(A,I5,A,I5)') ' minpermdist> move position ',NEWPERM(J3),' to ',J3
      NPERM=NPERM+1
   ENDIF
   DISTANCE=DISTANCE+(DUMMYA(3*(J3-1)+1)-DUMMYB(3*(J3-1)+1))**2 &
  &                    +(DUMMYA(3*(J3-1)+2)-DUMMYB(3*(J3-1)+2))**2 &
  &                    +(DUMMYA(3*(J3-1)+3)-DUMMYB(3*(J3-1)+3))**2
ENDDO

IF ((NPERM.NE.0).OR.(NTRIES.EQ.1)) THEN 

   CALL NEWMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,BULKT,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT)
   RMATCUMUL=MATMUL(RMAT,RMATCUMUL)
   DISTANCE=DISTANCE**2 ! we are using DISTANCE^2 further down

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
   RMATBEST(1:3,1:3)=RMATCUMUL(1:3,1:3)
   ROTINVBBEST(1:3,1:3)=ROTINVB(1:3,1:3) 
   ROTABEST(1:3,1:3)=ROTA(1:3,1:3)      
   RMATBEST=MATMUL(RMATBEST,ROTABEST)
   IF (INVERT.EQ.-1) THEN
      IF (PULLT.OR.TWOD.OR.(GTHOMSONT.AND.(GTHOMMET < 5))) THEN ! reflect in xz plane rather than invert!
         RMATBEST(1:3,1:3)=MATMUL(RMATBEST,REFXZ)
      ELSE
         RMATBEST(1:3,1:3)=-RMATBEST(1:3,1:3)
      ENDIF
   ENDIF
ENDIF

IF (DISTANCE.LT.GEOMDIFFTOL) GOTO 50
IF (NCHOOSE2.LT.NORBIT2) GOTO 30
IF (NCHOOSE1.LT.NORBIT1) GOTO 65
IF (NCHOOSEB2.LT.NORBITB2) GOTO 31
IF (NCHOOSEB1.LT.NORBITB1) GOTO 66

IF ((NCHOOSE2.EQ.NORBIT2).AND.(NCHOOSE1.EQ.NORBIT1).AND.(INVERT.EQ.1)) THEN

   IF (NOINVERSION.OR.BULKT.OR.(CHARMMT.AND.(.NOT.MACROCYCLET)).OR.(AMBERT.AND.(.NOT.MACROCYCLET)) &
  &     .OR.(NFREEZE.GT.0)) GOTO 50 
!  IF (DEBUG) PRINT '(A)','minpermdist> inverting geometry for comparison with target'
   INVERT=-1
   GOTO 60
ENDIF

IF (NROTDONE.LT.NRANROT) GOTO 11

50 DISTANCE=DBEST
!
!  XBEST contains the best alignment of A coordinates for the orientation of B coordinates in DUMMYB.
!  Rotate XBEST by ROTINVB to put in best correspondence with COORDSB, undoing the reorientation to DUMMYB from MYORIENT. 
!  We should get the same result for ROTINVB * RMATBEST * (COORDSA-CMA) 
!  where RMATBEST = +/- RMATCUMUL * ROTA for the best alignment 
!  (aside from a possible permutation of the atom ordering)
!
XDUMMY=0.0D0
DO J1=1,NATOMS
   XBEST(3*(J1-1)+1:3*(J1-1)+3)=MATMUL(ROTINVBBEST,XBEST(3*(J1-1)+1:3*(J1-1)+3))
   XDUMMY=XDUMMY+(COORDSB(3*(J1-1)+1)-XBEST(3*(J1-1)+1))**2+ &
  &                    (COORDSB(3*(J1-1)+2)-XBEST(3*(J1-1)+2))**2+ &
  &                    (COORDSB(3*(J1-1)+3)-XBEST(3*(J1-1)+3))**2
ENDDO
IF (ABS(SQRT(XDUMMY)-SQRT(DISTANCE)).GT.GEOMDIFFTOL .AND. (.NOT. RBAAT)) THEN
   PRINT '(2(A,G20.10))','minpermdist> ERROR *** distance between transformed XBEST and COORDSB=',SQRT(XDUMMY), &
  &                         ' should be ',SQRT(DISTANCE)
   PRINT '(A)','transformed XBEST:'
   PRINT '(3F20.10)',XBEST(1:3*NATOMS)
   PRINT '(A)','COORDSB:'
   PRINT '(3F20.10)',COORDSB(1:3*NATOMS)
ENDIF

RMATBEST=MATMUL(ROTINVBBEST,RMATBEST)
COORDSA(1:3*NATOMS)=XBEST(1:3*NATOMS) ! finally, best COORDSA should include permutations for DNEB input!
DISTANCE=SQRT(DISTANCE)

RETURN
END SUBROUTINE HKMINPERMDIST


SUBROUTINE GTHOMSONMINPERMDIST(COORDSB,COORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGID,RMATBEST,USEINT)

  USE COMMONS,ONLY : GTHOMMET
  IMPLICIT NONE
  INTEGER NATOMS, JJ
  DOUBLE PRECISION XCOORDSA (3*NATOMS), XCOORDSB(3*NATOMS), DD, DD2
  DOUBLE PRECISION DIST2, COORDSA(3*NATOMS), COORDSB(3*NATOMS), DISTANCE
  DOUBLE PRECISION BOXLX,BOXLY,BOXLZ,RMATBEST(3,3),RMATBEST2(3,3),REFXY(3,3)
  LOGICAL DEBUG, TWOD, RIGID, BULKT, USEINT

  REFXY(:,:) = 0.0D0
  REFXY(1,1) = 1.0D0; REFXY(2,2) = 1.0D0; REFXY(3,3) = -1.0D0

  IF (GTHOMMET .EQ. 5) THEN
    CALL GTHOMSONANGTOC(XCOORDSA,COORDSA,NATOMS)
    CALL GTHOMSONANGTOC(XCOORDSB,COORDSB,NATOMS)
  ELSE
    XCOORDSA = COORDSA
    XCOORDSB = COORDSB
  END IF

  DISTANCE = 1.0D10

  CALL HKMINPERMDIST(XCOORDSB,XCOORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DD,DD2,RIGID,RMATBEST,USEINT)
     
  IF (DD < DISTANCE) THEN
     DISTANCE = DD
     DIST2 = DD2
     COORDSA = XCOORDSA
     COORDSB = XCOORDSB
  ENDIF
!  PRINT *, 1, DD, DD2

  IF (GTHOMMET < 5) THEN
     DO JJ = 1, NATOMS
        XCOORDSA(3*JJ) = -XCOORDSA(3*JJ)
     ENDDO
     CALL HKMINPERMDIST(XCOORDSB,XCOORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DD,DD2,RIGID,RMATBEST2,USEINT)
     IF (DD < DISTANCE) THEN
        DISTANCE = DD
        DIST2 = DD2
        COORDSA = XCOORDSA
        COORDSB = XCOORDSB
        RMATBEST(:,:)=MATMUL(REFXY,RMATBEST)
        RMATBEST(:,:)=MATMUL(RMATBEST2,RMATBEST)
     ENDIF
  ENDIF
!  PRINT *, 2, DD, DD2


  RETURN
  
END SUBROUTINE GTHOMSONMINPERMDIST

!----------------------------------------------------------------------------------------------------------------------------------
!  Subroutine to convert Cartesians to theta, phi.
!  Order of output coordinates is phi1, theta1, r1, phi2, theta2, r2,...
!  where phi is the azimuthal angle and theta the vertical angle.
!  r is set to zero for all atoms and recalculated from the surface parameterisation when converting back to Cartesians
!
!  COORDS: the input coordinates, which should be Cartesians
!  P: the output coordinates in polars
!  NATOMS: the number of bodies
!  MYUNIT: file handle for debug printing. 6 is standard out.

      SUBROUTINE GTHOMSONCTOANG(COORDS, P, NATOMS, MYUNIT)
      USE COMMONS, ONLY: DEBUG, GTHOMMET, GTHOMSONZ
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NATOMS, MYUNIT
      INTEGER J1
      DOUBLE PRECISION, INTENT(INOUT)  :: COORDS(*)
      DOUBLE PRECISION, INTENT(OUT) :: P(*)
      DOUBLE PRECISION :: PI, HALFPI, phi, u, DIST
      DOUBLE PRECISION RADIUS
      DOUBLE PRECISION, PARAMETER :: TOLERANCE = 1.D-8
      LOGICAL POLAR

      PI = 4.0D0*ATAN(1.0D0)
      HALFPI = 2.0D0*ATAN(1.0D0)

      IF (GTHOMMET .NE. 5) THEN
        WRITE (*,*) 'GTHOMSONCTOANG> ERROR, conversion only implemented for spherical system'
        STOP
      END IF

! jwrm2> Attempt to detect if the coordintes are already polar. We'll assume they are if all the x coordinates
!        are between 0 and 2pi, all the y coordinates are between 0 and pi and all the z coordinates are zero, 
!        This is not guaranteed to mean they're polar, but should be ok since it shouldn't be a 2D system.
      POLAR = .TRUE.
      DO J1 = 3, 3*NATOMS, 3
        IF (COORDS(J1-2) .GT. 2*PI   .OR. COORDS(J1-2) .LT. 0.D0) POLAR = .FALSE. ! x coord
        IF (COORDS(J1-1) .GT. PI     .OR. COORDS(J1-1) .LT. 0.D0) POLAR = .FALSE. ! y coord
        IF (COORDS(J1) .GT. TOLERANCE .OR. COORDS(J1) .LT. -TOLERANCE) POLAR = .FALSE. ! z coord
      END DO
      IF (POLAR) THEN
        IF (DEBUG) WRITE (MYUNIT,*) 'GTHOMSONCTOANG> warning - incoming coordinates are already polar, skipping conversion'
        P(1:3*NATOMS) = COORDS(1:3*NATOMS)
        RETURN
      END IF

      IF (GTHOMMET .EQ. 5) THEN
         RADIUS = GTHOMSONZ
      ENDIF

      DO J1 = 1, NATOMS

         IF ( (COORDS(3*J1-2) .GE. 0.0D0) .AND. (COORDS(3*J1-1) .GE. 0.0D0) ) THEN
            IF ( ABS(COORDS(3*J1-2)) < 1.0D-5) THEN
               P(3*J1-2) = HALFPI
            ELSE IF ( ABS(COORDS(3*J1-1)) < 1.0D-5) THEN
               P(3*J1-2) = 0.0D0
            ELSE
               P(3*J1-2) = ATAN(COORDS(3*J1-1)/COORDS(3*J1-2))
            ENDIF
         ELSEIF ( (COORDS(3*J1-2) < 0.0D0) .AND. (COORDS(3*J1-1) .GE. 0.0D0) ) THEN
            IF ( ABS(COORDS(3*J1-1)) < 1.0D-5) THEN
               P(3*J1-2) = 2*HALFPI
            ELSE
               P(3*J1-2) = 2*HALFPI - ATAN(COORDS(3*J1-1)/(-COORDS(3*J1-2)))
            ENDIF
         ELSEIF ( (COORDS(3*J1-2) < 0.0D0) .AND. (COORDS(3*J1-1) < 0.0D0) ) THEN
               P(3*J1-2) = 2*HALFPI + ATAN(COORDS(3*J1-1)/COORDS(3*J1-2))
         ELSEIF ( (COORDS(3*J1-2) .GE. 0.0D0) .AND. (COORDS(3*J1-1) < 0.0D0) ) THEN
            IF ( ABS(COORDS(3*J1-2)) < 1.0D-5) THEN
               P(3*J1-2) = 3*HALFPI
            ELSE
               P(3*J1-2) = 4*HALFPI - ATAN(-COORDS(3*J1-1)/COORDS(3*J1-2))
            ENDIF
         ENDIF

         IF ( GTHOMMET .EQ. 5 )  THEN
            
            IF ( COORDS(3*J1) < 0.0D0 ) THEN               
               P(3*J1-1) = 2*HALFPI - ACOS(-COORDS(3*J1)/RADIUS)
            ELSE
               IF ( COORDS(3*J1)/RADIUS .GE. 1.0D0 ) THEN
                  P(3*J1-1) = 0.0D0
               ELSE IF ( COORDS(3*J1)/RADIUS .LE. -1.0D0 ) THEN
                  P(3*J1-1) = 2*HALFPI
               ELSE
                  P(3*J1-1) = ACOS(COORDS(3*J1)/RADIUS)
               ENDIF
            ENDIF

         ENDIF
         
         IF ( P(3*J1-2) > 4*HALFPI ) P(3*J1-2) = P(3*J1-2) - 4*HALFPI 
         IF ( P(3*J1-2) < 0.0D0    ) P(3*J1-2) = P(3*J1-2) + 4*HALFPI 
         IF ( P(3*J1-1) > 4*HALFPI ) P(3*J1-1) = P(3*J1-1) - 4*HALFPI 
         IF ( P(3*J1-1) < 0.0D0    ) P(3*J1-1) = P(3*J1-1) + 4*HALFPI 

      ENDDO

      DO J1 = 3, 3*NATOMS, 3
         P(J1) = 0.0D0
      END DO

      RETURN
      
    END SUBROUTINE GTHOMSONCTOANG

!----------------------------------------------------------------------------------------------------------------------------------
!  Subroutine to convert theta, phi to Cartesians.
!  Coordinates come in as phi1, theta1, r1, phi2, theta2, r2,...
!  where phi is the azimuthal angle and theta the vertical angle.
!  All rs are zero and the radius is calculated here from the surface parameterisation
!
!  COORDS: output coordinates, in Cartesians
!  P: input coordinates in polars
!  NATOMS: the number of bodies

      SUBROUTINE GTHOMSONANGTOC(COORDS,P,NATOMS)
      USE COMMONS, ONLY: DEBUG, GTHOMMET, GTHOMSONZ
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: NATOMS
      INTEGER J1
      DOUBLE PRECISION, INTENT(OUT)  :: COORDS(*)
      DOUBLE PRECISION, INTENT(IN)   :: P(*)
      DOUBLE PRECISION :: PI, FELINT, SELINT, u
      DOUBLE PRECISION :: ARG1, COSHARG
      DOUBLE PRECISION :: RADIUS
      DOUBLE PRECISION, PARAMETER :: TOLERANCE = 1.D-8
      LOGICAL POLAR

      PI = 4.0D0*ATAN(1.0D0)

      IF (GTHOMMET .NE. 5) THEN
        WRITE (*,*) 'GTHOMSONCTOANG> ERROR, conversion only implemented for spherical system'
        STOP
      END IF

! jwrm2> Detect if the coordintes are already Cartesian. If the r coordinates are not zero, then
!        the coordinates are already Cartesian. They could still be 
!        Cartesian if they are all zero, but this it unlikely since it shouldn't be a flat 2D system.
      POLAR = .TRUE.
      DO J1 = 3, 3*NATOMS, 3
        IF (P(J1) .GT. TOLERANCE .OR. P(J1) .LT. -TOLERANCE) POLAR = .FALSE. ! r coord
      END DO
      IF (.NOT. POLAR) THEN
        IF (DEBUG) PRINT *, 'GTHOMSONANGTOC> warning - incoming coordinates are already Cartesian, skipping conversion'
        COORDS(1:3*NATOMS) = P(1:3*NATOMS)
        RETURN
      END IF

      IF (GTHOMMET .EQ. 5) THEN
         RADIUS = GTHOMSONZ
      ENDIF

      DO J1 = 1, NATOMS
         IF ( GTHOMMET .EQ. 5) THEN

            ARG1 = RADIUS * SIN(P(3*J1-1))
            COORDS(3*J1-2)= ARG1 * COS(P(3*J1-2))
            COORDS(3*J1-1)= ARG1 * SIN(P(3*J1-2))
            COORDS(3*J1  )= RADIUS * COS(P(3*J1-1))

         ENDIF
      ENDDO
      RETURN
    END SUBROUTINE GTHOMSONANGTOC
