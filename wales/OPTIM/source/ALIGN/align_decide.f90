SUBROUTINE ALIGN_DECIDE(COORDSB,COORDSA,NATOMS,DEBUG,NBOXLX,NBOXLY,NBOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGID,RMATBEST)

USE KEY, ONLY: FASTOVERLAPT, BNB_ALIGNT, &    ! Logicals to determine which alignment routine to use
               KERNELWIDTH,NDISPLACEMENTS, &  ! Parameters for the Bulk FASTOVERLAP routine
               MAX_ANGMOM, NROTATIONS, &      ! Parameters for the Cluster FASTOVERLAP routine
               BNB_NSTEPS, &                  ! Parameter for the BNB align routine    
               BULK_BOXVEC, &                 ! Misc variables from the main program
               NSETS, PERMDIST, LOCALPERMDIST, NOINVERSION

USE GENRIGID, ONLY: RIGIDINIT, ATOMRIGIDCOORDT    ! Keywords that need checking for compatibility
USE BULKFASTOVERLAP, ONLY: FOM_ALIGN_BULK
USE CLUSTERFASTOVERLAP, ONLY: FOM_ALIGN_CLUSTERS, ALIGNHARM
USE GOPERMDIST, ONLY: BNB_ALIGN

IMPLICIT NONE

INTEGER NATOMS
DOUBLE PRECISION DIST2, COORDSA(3*NATOMS), COORDSB(3*NATOMS), DISTANCE, RMATBEST(3,3)
LOGICAL DEBUG, TWOD, RIGID, BULKT
DOUBLE PRECISION NBOXLX,NBOXLY,NBOXLZ

! Start by performing some sanity checks, to make sure the keywords being used are compatible with the requested alignment method.

IF (DEBUG .AND. BULKT .AND. ((ABS(NBOXLX-BULK_BOXVEC(1)).GT.1.0D-8) .OR. (ABS(NBOXLY-BULK_BOXVEC(2)).GT.1.0D-8) .OR. (ABS(NBOXLZ-BULK_BOXVEC(3)).GT.1.0D-8))) THEN
   WRITE(*,*) "align_decide> ERROR: Box parameters passed in as arguments differ to those USEd from COMMONS."
   WRITE(*,*) "Passed in: ", NBOXLX,NBOXLY,NBOXLZ
   WRITE(*,*) "USEd: ", BULK_BOXVEC(:)
   STOP 1
ENDIF  

IF (FASTOVERLAPT .OR. BNB_ALIGNT) THEN
   IF ((RIGIDINIT .AND. (.NOT.ATOMRIGIDCOORDT)) .OR. RIGID) THEN
      WRITE(*,*) "align_decide> fastoverlap and BNB methods do not work in rigid body coordinates. Use cartesians instead."
      STOP
   ELSEIF (ALLOCATED(NSETS)) THEN
      IF (ANY(NSETS(:).GT.0)) THEN
         WRITE(*,*) "align_decide> fastoverlap and BNB methods are not tested for secondary permutable sets, and probably don't work. Stopping now."
         STOP
      ENDIF
   ENDIF
ENDIF

! Now perform the actual alignment call.

! FASTOVERLAP and BNB are not designed to work with LOCALPERMDIST.
! Using these new routines for non-permutable systems (without PERMDIST) is unnecessary and likely to be inefficient.
! So in both of these cases, we bypass the new alignment routines and go straight to MINPERMDIST.
IF (PERMDIST .AND. (.NOT. LOCALPERMDIST)) THEN

   IF (FASTOVERLAPT) THEN  ! Without PERMDIST, we definitely don't need to call the ALIGN routines.
   
      IF(BULKT) THEN

         IF (.NOT. NOINVERSION) THEN
            IF (DEBUG) THEN
               WRITE(*,*) "align_decide> Warning: Bulk FASTOVERLAP does not support checking for inversion symmetry only. Setting NOINVERSION=.TRUE."
               WRITE(*,*) "align_decide> Use the OHCELL keyword to account for symmetries of a cubic box."
            ENDIF
            NOINVERSION = .TRUE.
         ENDIF

         IF (DEBUG) WRITE(*,*) "align_decide> using fastoverlap periodic alignment"
         CALL FOM_ALIGN_BULK(COORDSB,COORDSA,NATOMS,DEBUG,NBOXLX,NBOXLY,NBOXLZ,KERNELWIDTH,NDISPLACEMENTS,DISTANCE,DIST2)
      ELSE
         IF (DEBUG) WRITE(*,*) "align_decide> using fastoverlap cluster alignment"
         CALL FOM_ALIGN_CLUSTERS(COORDSB,COORDSA,NATOMS,DEBUG,MAX_ANGMOM,KERNELWIDTH,DISTANCE,DIST2,RMATBEST,NROTATIONS)
      ENDIF

   ELSE IF (BNB_ALIGNT) THEN

      IF(DEBUG) WRITE(*,*) "align_decide> using BNB align"
      CALL BNB_ALIGN(COORDSB,COORDSA,NATOMS,DEBUG,NBOXLX,NBOXLY,NBOXLZ,BULKT,DISTANCE,DIST2,RMATBEST,BNB_NSTEPS)

   ELSE
 
      IF(DEBUG) WRITE(*,*) "align_decide> using original MINPERMDIST routine"
      CALL MINPERMDIST(COORDSB,COORDSA,NATOMS,DEBUG,NBOXLX,NBOXLY,NBOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGID,RMATBEST)

   ENDIF

ELSE 

  IF (DEBUG .AND. (FASTOVERLAPT .OR. BNB_ALIGNT)) THEN
      WRITE(*,*) "Warning: Specified new ALIGN routines without PERMDIST or with LOCALPERMDIST. Using MINPERMDIST instead."
   ENDIF
   IF(DEBUG) WRITE(*,*) "align_decide> using original MINPERMDIST routine"
   CALL MINPERMDIST(COORDSB,COORDSA,NATOMS,DEBUG,NBOXLX,NBOXLY,NBOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGID,RMATBEST)

 ENDIF


END SUBROUTINE
