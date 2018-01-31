SUBROUTINE ALIGN_DECIDE(COORDSB,COORDSA,NATOMS,DEBUG,NBOXLX,NBOXLY,NBOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGID,RMATBEST)

USE COMMONS, ONLY: FASTOVERLAPT, BNB_ALIGNT, &    ! Logicals to determine which alignment routine to use
                   KERNELWIDTH,NDISPLACEMENTS, &  ! Parameters for the Bulk FASTOVERLAP routine
                   MAX_ANGMOM, NROTATIONS, &      ! Parameters for the Cluster FASTOVERLAP routine
                   BNB_NSTEPS, &                  ! Parameter for the BNB align routine    
                   MYUNIT, BOXLX, BOXLY, BOXLZ, & ! Misc variables from the main program
                   NSETS, PERMOPT, PERMDIST, LOCALPERMDIST, PERMINVOPT, NOINVERSION

USE GENRIGID, ONLY: RIGIDINIT, ATOMRIGIDCOORDT    ! Keywords that need checking for compatibility
USE BULKFASTOVERLAP, ONLY: FOM_ALIGN_BULK
USE CLUSTERFASTOVERLAP, ONLY: FOM_ALIGN_CLUSTERS, ALIGNHARM
USE GOPERMDIST, ONLY: BNB_ALIGN

IMPLICIT NONE

INTEGER NATOMS
DOUBLE PRECISION DIST2, COORDSA(3*NATOMS), COORDSB(3*NATOMS), DISTANCE, RMATBEST(3,3)
LOGICAL DEBUG, TWOD, RIGID, BULKT, SAVEPERMOPT, SAVEPERMINVOPT
DOUBLE PRECISION NBOXLX,NBOXLY,NBOXLZ

! Start by performing some sanity checks, to make sure the keywords being used are compatible with the requested alignment method.

IF (DEBUG .AND. BULKT .AND. ((ABS(NBOXLX-BOXLX).GT.1.0D-8) .OR. (ABS(NBOXLY-BOXLY).GT.1.0D-8) .OR. (ABS(NBOXLZ-BOXLZ).GT.1.0D-8))) THEN
   WRITE(MYUNIT,*) "align_decide> ERROR: Box parameters passed in as arguments differ to those USEd from COMMONS."
   WRITE(MYUNIT,*) "Passed in: ", NBOXLX,NBOXLY,NBOXLZ
   WRITE(MYUNIT,*) "USEd: ", BOXLX, BOXLY, BOXLZ
   STOP 1
ENDIF  

IF (FASTOVERLAPT .OR. BNB_ALIGNT) THEN
   IF ((RIGIDINIT .AND. (.NOT.ATOMRIGIDCOORDT)) .OR. RIGID) THEN
      WRITE(MYUNIT,'(A)') "align_decide> fastoverlap and BNB methods do not work in rigid body coordinates. Use cartesians instead."
      STOP
   ELSEIF (ANY(NSETS(:).GT.0)) THEN
      WRITE(MYUNIT,'(A)') "align_decide> fastoverlap and BNB methods is not tested for secondary permutable sets, and probably doesn't work. Stopping now."
      STOP
   ENDIF

   ! In order to ensure that the correct logic is followed in MINPERMDIST, we need to make sure that PERMOPT and PERMINVOPT
   ! have the correct values for the system type. Unfortunately, we need PERMOPT set for clusters even if it's not being used
   ! for the rest of the program.
   SAVEPERMOPT = PERMOPT; SAVEPERMINVOPT = PERMINVOPT
   IF(BULKT) THEN
      PERMOPT = .FALSE.
      PERMINVOPT = .FALSE.
   ELSE
      PERMOPT = .TRUE.
      IF(.NOT.NOINVERSION) PERMINVOPT=.TRUE.  ! Inversion isomers will be identified unless you have NOINVERSION set!
   ENDIF

ENDIF

! Now perform the actual alignment call.

! FASTOVERLAP and BNB are not designed to work with LOCALPERMDIST.
! Using these new routines for non-permutable systems (without PERMDIST) is unnecessary and likely to be inefficient.
! So in both of these cases, we bypass the new alignment routines and go straight to MINPERMDIST.
IF (PERMDIST .AND. (.NOT. LOCALPERMDIST)) THEN

   IF (FASTOVERLAPT) THEN
   
      IF(BULKT) THEN
         IF (DEBUG) WRITE(MYUNIT,*) "align_decide> using fastoverlap periodic alignment"
         CALL FOM_ALIGN_BULK(COORDSB,COORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,KERNELWIDTH,NDISPLACEMENTS,DISTANCE,DIST2)
      ELSE
         IF (DEBUG) WRITE(MYUNIT,*) "align_decide> using fastoverlap cluster alignment"
         CALL FOM_ALIGN_CLUSTERS(COORDSB,COORDSA,NATOMS,DEBUG,MAX_ANGMOM,KERNELWIDTH,DISTANCE,DIST2,RMATBEST,NROTATIONS)
      ENDIF

   ELSE IF (BNB_ALIGNT) THEN

      IF(DEBUG) WRITE(MYUNIT,*) "align_decide> using BNB align"
      CALL BNB_ALIGN(COORDSB,COORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,DISTANCE,DIST2,RMATBEST,BNB_NSTEPS)

   ELSE
      IF(DEBUG) WRITE(MYUNIT,*) "align_decide> using original MINPERMDIST routine"
      CALL MINPERMDIST(COORDSB,COORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGID,RMATBEST)

   ENDIF

ELSE  ! i.e. LOCALPERMDIST or NOT PERMDIST.

   IF (DEBUG .AND. (FASTOVERLAPT .OR. BNB_ALIGNT)) THEN
      WRITE(MYUNIT,*) "Warning: Specified new ALIGN routines without PERMDIST or with LOCALPERMDIST. Using MINPERMDIST instead."
   ENDIF
   IF(DEBUG) WRITE(MYUNIT,'(A)') "align_decide> using original MINPERMDIST routine"
   CALL MINPERMDIST(COORDSB,COORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGID,RMATBEST)
ENDIF

IF (FASTOVERLAPT .OR. BNB_ALIGNT) THEN
   PERMOPT = SAVEPERMOPT; PERMINVOPT = SAVEPERMINVOPT
ENDIF

END SUBROUTINE
