MODULE INTCUTILS
  USE MODCHARMM
  USE COMMONS
  USE INTCOMMONS
  USE KEY, ONLY : INTEPSILON
  USE SPFUNCTS, ONLY : DUMPCOORDS

  IMPLICIT NONE

CONTAINS
    SUBROUTINE INTINTERPOLATE(RS,RF,NIMI,NIMC, XINTERP, PTEST,OUTFAILED)
    ! do a linear interpolation between start and end 
    ! using internal coordinates, with NIMI images
    ! then place NIMC images at equidistant cartesian
    ! points along the resulting interpolated path    
    ! results go in XCART

    USE KEY, ONLY : RIGIDBODY, TWOD, BULKT, AMBERT, NABT

    IMPLICIT NONE

    DOUBLE PRECISION, INTENT(IN) :: RS(3*NATOMS), RF(3*NATOMS)
    INTEGER, INTENT(IN) :: NIMC
    INTEGER, INTENT(INOUT) :: NIMI
    DOUBLE PRECISION, INTENT(OUT) :: XINTERP(3*NATOMS*NIMC)
    LOGICAL, INTENT(IN) :: PTEST
    LOGICAL, INTENT(OUT) :: OUTFAILED

    DOUBLE PRECISION :: RSINT(NINTC), RFINT(NINTC)
    DOUBLE PRECISION :: XINT(NINTC), INTDIFF(NINTC)
    DOUBLE PRECISION,ALLOCATABLE :: XHELP(:,:), DIFFS2(:,:)
    INTEGER :: IM, I2, NC, I3
    DOUBLE PRECISION, ALLOCATABLE :: XCART(:,:), DIFFS(:,:), DISTS(:), ARCS(:)
    DOUBLE PRECISION :: WANTLEN, EXTRALEN, TOTLEN

    DOUBLE PRECISION :: EMAX, EMAXC, ENERGY, GDUMMY(3*NATOMS), TMPRMS
    DOUBLE PRECISION :: XCART2(3*NATOMS*NIMC), CARTDIFF(3*NATOMS), GINTDUMMY(NINTC)
    INTEGER :: ISTAT
    LOGICAL :: FAILED
    DOUBLE PRECISION :: RFORIG(3*NATOMS)

    ! alignment stuff
    DOUBLE PRECISION :: DISTF, DIST, DIST2, RMAT(3,3)
    CHARACTER(LEN=5) :: ZSYMSAVE
    COMMON /SYS/ ZSYMSAVE
    DOUBLE PRECISION TEST(3*NATOMS)

    IF (.NOT.NATINT) THEN
       print*, 'intinterpolate >> INTINTERPOLATE is not set up to work without NATINT'
       STOP
    ENDIF    

    OUTFAILED = .FALSE.

    BACKTCUTOFF = INTERPBACKTCUT

    IF (INTERPSIMPLE) NIMI = NIMC
!    NIMI = NIMC

    ALLOCATE(XCART(0:NIMI+1,3*NATOMS), DIFFS(0:NIMI+1,3*NATOMS), DISTS(0:NIMI+1), ARCS(0:NIMI+1))
    ALLOCATE(DIHINFO(NIMI+2,NDIH))
    ALLOCATE(DIFFS2(0:600,NINTC)) !debug
    ALLOCATE(XHELP(0:600,NINTC))
    DIHINFO(:,:) = 0.0D0

    IF (PTEST) print*, 'intinterpolate>> interpolating with internals: ', &
         & NIMI, NIMC
              
    NC = 3*NATOMS

    IF (INTERPCHOICE) THEN
       EMAXC = -9999
       CARTDIFF(:) = RF(:)-RS(:)
       DO IM = 1,NIMC
          XCART2(NC*(IM-1)+1:NC*IM) = RS(:) + CARTDIFF*DBLE(IM)/(NIMC+1)
          CALL POTENTIAL(XCART2(NC*(IM-1)+1:NC*IM), ENERGY, GDUMMY, .FALSE., .FALSE., &
               & TMPRMS, .FALSE., .FALSE.)
          
          IF (ENERGY.GT.EMAXC) EMAXC = ENERGY
       ENDDO
    ENDIF

    IF (NATINT) PREVDIH => DIHINFO(1,:)
    ALIGNDIR = .FALSE.
    CALL CART2INT(RS,RSINT)
   
    ! align individual dihedrals of all other images and endpoint to start
    DO IM = 2,NIMI+2
       DIHINFO(IM,:) = DIHINFO(1,:)
    ENDDO
    PREVDIH => DIHINFO(NIMI+2,:)
    ALIGNDIR = .TRUE.
    CALL CART2INT(RF,RFINT)

    ALIGNDIR = .FALSE.
    
    INTDIFF = RFINT(:) - RSINT(:)

    ! place the images in internals and convert to cartesians
    IF (PTEST) THEN
       print*, 'intinterpolate>> step size: ', SQRT(DOT_PRODUCT(INTDIFF,INTDIFF)) / (NIMI+1)
!      CALL FLUSH(6)
    ENDIF

    XHELP(0,1:NINTC) = RSINT(:)
    XCART(0,1:NC) = RS(:)
    ARCS(0) = 0.0D0
    DO IM = 1,NIMI
       PREVDIH => DIHINFO(IM+1,:)
       IF (AMBERT.OR.NABT) THEN !msb50
          CALL AMB_TRANSBACKDELTA(INTDIFF(:) * 1.0D0 / (NIMI+1),DIFFS(IM,:), &
   &         XCART(IM-1,:),NINTC,3*NATOMS,NNZ,KD,FAILED,.FALSE.,INTEPSILON)
       ELSE
       CALL TRANSBACKDELTA(INTDIFF(:) * 1.0D0 / (NIMI+1),DIFFS(IM,:),XCART(IM-1,:),&
            & NINTC,3*NATOMS,NNZ,KD,FAILED,.FALSE.,INTEPSILON)       
       ENDIF
       XCART(IM,:) = XCART(IM-1,:) + DIFFS(IM,:)                       

       ! line up structure with start point
       CALL NEWMINDIST(RS(:),XCART(IM,:),NATOMS,DIST,.FALSE.,.FALSE.,ZSYM(1),.FALSE.,.FALSE.,.FALSE.,RMAT)

       DIFFS(IM,:) = XCART(IM,:) - XCART(IM-1,:)
       DISTS(IM) = SQRT(DOT_PRODUCT(DIFFS(IM,:), DIFFS(IM,:)))
       ARCS(IM) = ARCS(IM-1) + DISTS(IM)

       IF (INTERPSIMPLE) XINTERP(NC*(IM-1)+1:NC*IM) = XCART(IM,:)

    ENDDO   

    DIFFS(NIMI+1,:) = RF(:) - XCART(NIMI,:)
    DISTS(NIMI+1) = SQRT(DOT_PRODUCT(DIFFS(NIMI+1,:), DIFFS(NIMI+1,:)))
    ARCS(NIMI+1) = ARCS(NIMI) + DISTS(NIMI+1)        
    
    IF (PTEST) THEN
       CALL DUMPCOORDS(XCART(0,:),'interp.xyz', .FALSE.)
       DO IM = 1,NIMI
          CALL DUMPCOORDS(XCART(IM,:),'interp.xyz', .TRUE.)
       ENDDO
       CALL DUMPCOORDS(RF(:),'interp.xyz', .TRUE.)
    ENDIF

    IF (DISTS(NIMI+1)/DISTS(NIMI).GT.10) THEN
       print*, 'Something wrong! Interpolation discontinuous at end.'
       print*, "nimi", nimi+1, nimi
       print*, DISTS(NIMI-6:NIMI+1)
       OUTFAILED = .TRUE.
!       IF (PTEST) THEN
          CALL DUMPCOORDS(RS,'badstart.xyz', .FALSE.)
          CALL DUMPCOORDS(RF, 'badfin.xyz', .FALSE.)
!       ENDIF
    ENDIF

    BACKTCUTOFF = MINBACKTCUT

    IF (INTERPSIMPLE) THEN
       DEALLOCATE(XCART, DIFFS, DISTS, ARCS)
       print*, 'Finished simple initial interpolation'
       RETURN
    ENDIF

    TOTLEN = ARCS(NIMI+1)

    EMAX = -9999

    ! get cartesians at equidistant intervals along interpolated path
    I2 = 0
    DO IM = 1,NIMC
       WANTLEN = TOTLEN*DBLE(IM)/(NIMC+1) !desired total length from start
       DO WHILE (ARCS(I2).LE.WANTLEN)
          I2 = I2 + 1
       ENDDO

       EXTRALEN = WANTLEN - ARCS(I2-1)

       XINTERP(NC*(IM-1)+1:NC*IM) = XCART(I2-1,:) + DIFFS(I2,:)*EXTRALEN/DISTS(I2)
       IF (INTERPCHOICE) THEN
          CALL POTENTIAL(XINTERP(NC*(IM-1)+1:NC*IM), ENERGY, GDUMMY, .FALSE., .FALSE., &
               & TMPRMS, .FALSE., .FALSE.)

          IF (ENERGY.GT.EMAX) EMAX = ENERGY
       ENDIF
       I2 = I2 - 1
    ENDDO

    IF (PTEST) THEN
       CALL DUMPCOORDS(RS(:),'interpfin.xyz', .FALSE.)
       DO IM = 1,NIMC
          CALL DUMPCOORDS(XINTERP(NC*(IM-1)+1:NC*IM),'interpfin.xyz', .TRUE.)
       ENDDO
       CALL DUMPCOORDS(RF(:),'interpfin.xyz', .TRUE.)
    ENDIF

    IF (INTERPCHOICE) THEN
       IF (EMAX.GT.EMAXC) THEN
         IF (PTEST) print '(A,2G20.10)', 'intinterpolate>> using cartesian interpolation', EMAXC, EMAX
          XINTERP(:) = XCART2(:)
       ELSE
          IF (PTEST) print '(A,2G20.10)', 'intinterpolate>> using internal interpolation', EMAXC, EMAX
       ENDIF
    ENDIF
!   CALL FLUSH(6)

!    PRINT*, "intinterp> deall"
    DEALLOCATE(XCART, DIFFS, DISTS, ARCS, DIHINFO)
    DEALLOCATE(DIFFS2,XHELP)


    IF (PTEST) print*, 'Finished initial interpolation'
    RETURN

  END SUBROUTINE INTINTERPOLATE

  SUBROUTINE INTMINPERM(RS,RF,DISTANCE,RMAT,PTEST)
    ! minimise Cartesian distance without permuting then
    ! permute the permutable groups to minimize distance in single torsions
    ! changes RF to resulting best alignment with RS
    ! output distance is cartesian distance **squared**
!msb50 - notes on intminperm:
! still expect problems when dragging permutable atoms along in a swap, as 
! technically to see whether this gives a smaller atoms the dragged atoms have
! to be permuted at the same time (i.e. two CH_3 groups as in val, leu)
! I currently don't do this - primarily because swaps for val, leu are re-checked
! in intminperm, so should be fine - but then alignment cartesian


      USE KEY, ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, SETS, AMBERT, NABT
      USE MODAMBER9, ONLY : ih, ix, NRES, i02,m04,m02
      USE intcommons, only: NDIH, PERMNEIGHBOURS, PERMCHAIN, NGLYDIH, GLYDIH, GLYCART !msb50
      USE commons, only: NATOMS !msb50
      USE PORFUNCS
    
      IMPLICIT NONE
    
      DOUBLE PRECISION, INTENT(IN)    :: RS(3*NATOMS)
      DOUBLE PRECISION, INTENT(INOUT) :: RF(3*NATOMS), RMAT(3,3)
      DOUBLE PRECISION, INTENT(OUT)   :: DISTANCE
      LOGICAL, INTENT(IN)             :: PTEST
    
      DOUBLE PRECISION :: RF2(3*NATOMS), RFNEW(3*NATOMS)
      INTEGER :: PG, SW, I, J, K, A1, A2, A3, B1, B2, B3, START, MINSWAP(3)
      DOUBLE PRECISION :: RSINT(NINTC), RFINT(NINTC), INTDIFF(NDIH), INTDIFFMSB50(NDIH), INTDISTMSB50
      DOUBLE PRECISION :: INTDIST, CURINTDIST, INTDIST2
      DOUBLE PRECISION :: GD_PREVDIH(NDIH)
      DOUBLE PRECISION :: TMP(3), TMP2(3), TMP3(3), STARTDIH(NDIH), RFTMP(3*NATOMS)
      DOUBLE PRECISION :: RFTMP2(3*NATOMS), MININTDIST, RFTMP3(3*NATOMS)
      DOUBLE PRECISION :: TMP_AR(NATOMS,3)
      INTEGER          :: MINSWAP1(3), MINSWAP2(3)
! al  ignment stuff
      CHARACTER(LEN=5) :: ZSYMSAVE
      COMMON /SYS/ ZSYMSAVE

      DOUBLE PRECISION :: DISTANCE2
      CHARACTER(LEN=4) :: resname
      INTEGER          :: III, ISTAT,JJJ
      INTEGER          :: nglyc, nprol, glyc(NRES), prol(NRES)
      LOGICAL          :: glyswap, prolswap, SKIPPG(NPERMGROUP)
      INTEGER          :: DUMMYPG1(3), DUMMYST1(3)
      INTEGER          :: SWAPARRAY(12,6), SWAPARRAY3(24,8), SWAPARRAY5(96,12)
      !contains all possibilities of having neighbouring groups of 
      !permable atoms permuted (i.e. a, b, c, d; b, a, c, d; a, b, d, c etc)
      DOUBLE PRECISION :: SWAPDIST(12), SWAPDIST3(24), SWAPDIST5(96)
      DOUBLE PRECISION :: GLYSTART(4), GLYFIN(4), GLYDIFF(4), GLYDIST1, GLYDIST
      DOUBLE PRECISION :: dihed1
      DOUBLE PRECISION :: PI
      INTEGER          :: i1, i2, i3, i4
      DOUBLE PRECISION :: ANGLES(NDIH)
      data DUMMYPG1 /0,0,0/
      data DUMMYST1 /0,0,0/ 
      data PI /3.1415926535897931/    

      SKIPPG(:)=.FALSE.
      IF (.NOT.ALLOCATED(PERMNEIGHBOURS)) CALL FINDPERMDIH
     
      print*, "msb50 just checking whether this is called in ambtransbackde" 

      prol(:) = 0;nprol=0
      glyc(:) = 0;nglyc=0
      glyswap = .FALSE.;prolswap=.FALSE.
      IF (AMBERT.OR.NABT) THEN !msb50 - for amber9.ff03 and glycine this doesn't work
          DO III = 1,NRES!as swapping the hydrogens makes no difference in natural internals in this case
            IF (ih(m02+III-1).EQ.'GLY'.OR.ih(m02+III-1).EQ.'NGLY'.OR.ih(m02+III-1).EQ.'CGLY') THEN
               nglyc = nglyc+1
               glyc(nglyc) = III
            ENDIF       
            IF (ih(m02+III-1).EQ.'PRO'.OR.ih(m02+III-1).EQ.'NPRO'.OR.ih(m02+III-1).EQ.'CPRO') THEN
               nprol = nprol+1 
               prol(nprol) =III
            ENDIF
          ENDDO
      ENDIF

      DISTANCE2 = DOT_PRODUCT(RF-RS,RF-RS)
      IF (PTEST) PRINT*, "intminperm> initial cart dist", DISTANCE2

! align without permuting
! msb50 - done at the end in minpermdist now!
!      CALL NEWMINDIST(RS(:),RF(:),NATOMS,DISTANCE,.FALSE.,.FALSE.,ZSYM(1),.FALSE.,.FALSE.,.FALSE.,RMAT)

      RFTMP3(:) = RF(:)

      DISTANCE2 = DOT_PRODUCT(RF-RS,RF-RS)
      IF (PTEST) PRINT*, "intminperm> initial cart dist", DISTANCE2

      PREVDIH => DIHINFOSINGLE(:)

      PREVDIH(:) = 0.0D0
      CALL CART2INT(RS, RSINT)

    STARTDIH(:) = PREVDIH(:)
    ALIGNDIR = .FALSE.

     CALL GETDIHONLY(RF) 

    INTDIFF(:) = PREVDIH(:)-STARTDIH(:)
! just for the moment

    INTDIFFMSB50(:) = INTDIFF(:)
    CURINTDIST = SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))
    IF (PTEST) print*, 'Starting intdistance: ', CURINTDIST

    RFTMP(:)=RF(:)
    IF (nglyc.NE.0.AND..NOT.GLYCART) THEN
        CALL GLYINTPERM(RS, RF, PTEST)
    ENDIF

    START = 1
    DO PG = 1,NPERMGROUP

    IF (NPERMSIZE(PG).EQ.2) THEN
      RFTMP(:) = RF(:)

      ! A1 and A2 are the atom numbers to permute
      A1 = PERMGROUP(START)
      A2 = PERMGROUP(START+1)

      IF (PTEST) PRINT '(A,2I6)', "intminperm> which PERMGROUP", PG, START
      IF (PTEST) PRINT '(A,2I6)','intminperm> Swapping atoms: ', A1, A2

!glycine can't be done with Lena's natural internals as they internals don't 
! change when the glycine H's are swapped
      IF (nglyc.NE.0) THEN
         DO III=1, nglyc
         IF (ix(i02+glyc(III)-1).GT.A1) EXIT
         glyswap = .FALSE.
            IF (ix(i02+glyc(III)-1).LE.A1 .AND.(ix(i02+glyc(III)).GT.A1).AND.&
       &          (ih(m04+A1-1).EQ.'HA2 '.OR.ih(m04+A1-1).EQ.'HA3 ')) THEN
               IF (PTEST) PRINT*, "glyc", glyc(III)
               IF (.NOT.GLYCART) THEN
                  glyswap=.TRUE.
                  CYCLE
               ENDIF 
              !the second atom has to be A2 in gly
               IF (ih(m04+A2-1).EQ.'HA3 '.OR. ih(m04+A2-1).EQ.'HA2 ') THEN
                  glyswap = .TRUE.
                  IF (PTEST) PRINT*, "glycine"
                  RF(3*(A1-1)+1:3*A1) = RFTMP(3*(A2-1)+1:3*A2)
                  RF(3*(A2-1)+1:3*A2) = RFTMP(3*(A1-1)+1:3*A1)
              !move any other groups that must be dragged along
                  IF (NSETS(PG).NE.0) THEN
                    PRINT*, "intminperm>this is not glycine. you DON'T know& 
        &                        what you are doing"
                    STOP
                  ENDIF
                  DISTANCE = DOT_PRODUCT(RFTMP-RS,RFTMP-RS)
                  DISTANCE2 = DOT_PRODUCT(RF-RS,RF-RS)
                  !IF (PTEST) PRINT*, "dist old", DISTANCE, "dist new", DISTANCE2
                  IF (DISTANCE2.LT.DISTANCE) THEN
                    IF (PTEST) print '(a40, f17.5,a10,f17.5)',"keep permutation&
        & aaccrding to cartesians (gly): dist1",DISTANCE, "dist_new", DISTANCE2
                  ELSE
                     RF(:) = RFTMP(:) !unddo
                  ENDIF
               ELSE
                  PRINT*, "intminperm> sth wrong with GLY"
                  STOP
               ENDIF
           ENDIF
           IF (glyswap) EXIT
         ENDDO
           IF (glyswap) THEN
              IF (PTEST) PRINT*, "glyswap true", PG
              START = START + NPERMSIZE(PG)
             CYCLE
           ENDIF
      ENDIF

      IF (nprol.NE.0) THEN
         prolswap=.FALSE.
         DO III=1, nprol

            IF (ix(i02+prol(III)-1).LT.A1 .AND. ix(i02+prol(III)).GT.A1) THEN
               IF (PTEST) PRINT*, "pro", prol(III)
               CALL PROL_PERMUTE(A1,A2,PTEST, RS, RF, DISTANCE)
               prolswap =.TRUE.
               EXIT
            ENDIF
         ENDDO
         IF (prolswap) THEN
            START= START+NPERMSIZE(PG)
            CYCLE
         ENDIF
      ENDIF

! now start the real thing:     
      IF (PERMNEIGHBOURS(PG,1).EQ.0) THEN
         RF(3*(A1-1)+1:3*A1) = RFTMP(3*(A2-1)+1:3*A2)
         RF(3*(A2-1)+1:3*A2) = RFTMP(3*(A1-1)+1:3*A1)
        
         !move any other groups that must be dragged along
         DO SW = 1,NSETS(PG)
!           A1 = SWAP1(PG,SW)
!           A2 = SWAP2(PG,SW)
            A1 = SETS(PERMGROUP(START),SW)
            A2 = SETS(PERMGROUP(START+1),SW)
        
            IF (PTEST) print*, 'Dragging swap: ', A1, A2
        
            RF(3*(A1-1)+1:3*A1) = RFTMP(3*(A2-1)+1:3*A2)
            RF(3*(A2-1)+1:3*A2) = RFTMP(3*(A1-1)+1:3*A1)
         ENDDO
       
         PREVDIH(:) = STARTDIH(:)
         CALL GETDIHONLY(RF) !msb50
   
         INTDIFF(:) = PREVDIH(:) - STARTDIH(:)
!!!!
         INTDIST = SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))
         IF (PTEST) PRINT '(2i4,a6,2f15.10)', A1, A2,"dists", INTDIST, CURINTDIST
         
         IF (PTEST) print*, 'torsion distance (not yet accepted!): ', INTDIST
    
    
         PREVDIH(:)=STARTDIH(:)

         
         RFNEW(:)=RFTMP(:)
         A1 = PERMGROUP(START)
         A2 = PERMGROUP(START+1)
         CALL DISTANCEPAIRSWAP(STARTDIH, RFTMP,A1,A2, PG, INTDIST2,PTEST,RFNEW)

         IF (ABS(INTDIST-INTDIST2).GT. 1e-6) THEN
             PRINT*, "not equal"
             STOP
         ENDIF

         IF (INTDIST.GE.CURINTDIST) THEN
            ! undo permutation
          RF(:) = RFTMP(:)
         ELSE
             IF (PTEST) print*, 'keep permutation'
   
             CURINTDIST = INTDIST
   
         ENDIF

      ELSE
         IF (PERMCHAIN(PG,1).EQ.2) THEN
            IF (SKIPPG(PG)) THEN
                 START = START + NPERMSIZE(PG)
                CYCLE
            ENDIF
            IF (PTEST) PRINT*, "coupled pg's!"
            CALL SWAP2ATONCE(STARTDIH,RF,PG,START,0,PTEST,SKIPPG,RFNEW,MINSWAP1,MINSWAP2,MININTDIST, SWAPARRAY, SWAPDIST)
            DO I =1, NATOMS
               CALL FLUSH(6)
               IF (RF(I).NE.RFNEW(I)) THEN
                 PRINT*, I, "not eq after swap2atonce", RF(I), RFNEW(I);
               !STOP
               ENDIF
            ENDDO
            RF(:)=RFNEW(:)
            IF (MININTDIST.LT.CURINTDIST) CURINTDIST=MININTDIST
            IF (PTEST) PRINT*, "MININTDIST after swap2atonce", CURINTDIST
          ELSEIF (PERMCHAIN(PG,1).EQ.3) THEN
            CALL GETRESID(A1,resname, I)
            IF (resname.EQ.'ARG '.OR.resname.EQ.'NARG'.OR.resname.EQ.'CARG') THEN
               IF (NSETS(PG)==2) THEN
                  PRINT*, "call argswap"
                  CALL ARG_SWAP(STARTDIH, RF,PG,START,PTEST,SKIPPG,RFNEW,CURINTDIST)
                  DO I=1,3*NATOMS 
                     IF (RFNEW(I).NE.RF(I)) THEN; PRINT*, "arg ne",I,RFNEW(I),RF(I); STOP;ENDIF
                  ENDDO
                  START=START+NPERMSIZE(PG)
               CYCLE
               ENDIF
            ENDIF
            IF (PTEST) PRINT*, "triple pg's"
            CALL SWAP3ATONCE(STARTDIH,RF,PG,START,DUMMYPG1,DUMMYST1,PTEST,SKIPPG,RFNEW,MININTDIST, SWAPARRAY3, SWAPDIST3)
            !DO I =1, NATOMS
            !   CALL FLUSH(6)
            !   IF (RF(I).NE.RFNEW(I)) THEN
            !     PRINT*, I, "not eq after swap3atonce", RF(I), RFNEW(I);
            !   !STOP
            !   ENDIF
            !ENDDO
            RF(:)=RFNEW(:)
            IF (MININTDIST.LT.CURINTDIST) CURINTDIST=MININTDIST
            IF (PTEST) PRINT*, "minintdist after swap3atonce", MININTDIST
          ELSEIF (PERMCHAIN(PG,1).EQ.4) THEN
            PRINT*, "intminperm>currently no case known. Code this yourself"
            STOP
          ELSEIF (PERMCHAIN(PG,1).EQ.5) THEN
             IF (PTEST) PRINT*, "this is lysine - 5pg's"
             CALL SWAP5ATONCE(STARTDIH,RF,PG,START,PTEST,SKIPPG,RFNEW,MININTDIST,SWAPARRAY5,SWAPDIST5)
             DO I =1, NATOMS
                CALL FLUSH(6)
                IF (RF(I).NE.RFNEW(I)) THEN
                  PRINT*, I, "not eq after swap5atonce"
                !STOP
                ENDIF
             ENDDO
             RF(:)=RFNEW(:)
             IF (MININTDIST.LT.CURINTDIST) THEN
                CURINTDIST=MININTDIST
                IF (PTEST) PRINT*, "MININTDIST", MININTDIST
             ENDIF
          ELSE
            PRINT*, "this should never ever happen. Better go back to bed. "
            STOP
          ENDIF
         !RF(:)=RFNEW(:)
       ENDIF
      
      
      ELSE IF (NPERMSIZE(PG).EQ.3) THEN
       ! this is an inefficient way of listing permutations
       ! and should be fixed and generalized at some point
        IF (PERMNEIGHBOURS(PG,1).EQ.1) THEN
           !next 5 lines to the bin once this is working 
           IF (PERMCHAIN(PG,1).GT.3) THEN
      !       PRINT*, "currently skipping", PG, PERMCHAIN(PG,1) 
             START = START + NPERMSIZE(PG)
             CYCLE
           ENDIF
             IF (SKIPPG(PG)) THEN
               START = START + NPERMSIZE(PG)
      !         PRINT*, "cycling", PG
                CYCLE
             ENDIF
             CALL SWAP2ATONCE(STARTDIH,RF,PG,START,0,PTEST,SKIPPG,RFNEW,MINSWAP1,MINSWAP2,MININTDIST,SWAPARRAY, SWAPDIST)
             RF(:)=RFNEW(:)
             IF (MININTDIST.LT.CURINTDIST) CURINTDIST=MININTDIST
       ELSE
         A1 = PERMGROUP(START)
         A2 = PERMGROUP(START+1)
         A3 = PERMGROUP(START+2)
        
         TMP(:) = RF(3*(A1-1)+1:3*A1)
         TMP2(:) = RF(3*(A2-1)+1:3*A2)
         TMP3(:) = RF(3*(A3-1)+1:3*A3)
        
         RFTMP2(:) = RF(:)
         RFTMP(:) = RF(:)
         !DO I = 0,2
         !   DO J = 0,2
         !      IF (J.EQ.I.OR.(I.EQ.0.AND.J.EQ.1)) CYCLE
         !      DO K = 0,2
         !         IF (K.EQ.I.OR.K.EQ.J) CYCLE
         ! 
         !         B1 = PERMGROUP(START+I)
         !         B2 = PERMGROUP(START+J)
         !         B3 = PERMGROUP(START+K)
         ! 
         !         IF (PTEST) print*, 'permuting: ', B1, B2, B3
        ! 
        !          RF(3*(B1-1)+1:3*B1) = TMP(:)
        !          RF(3*(B2-1)+1:3*B2) = TMP2(:)
        !          RF(3*(B3-1)+1:3*B3) = TMP3(:)     
        !       
        !          PREVDIH(:) = STARTDIH(:)
        !          CALL CART2INT(RF,RFINT)
         ! 
         !         INTDIFF(:) = PREVDIH(:) - STARTDIH(:)
         !          INTDIST = SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))
         ! 
         !          IF (PTEST) print*, 'dihedral distance: ', INTDIST
         ! 
         !          !PRINT*, "msb50 intdist", B1,B2,B3 ,intdist
         ! 
         !          IF (INTDIST.LT.CURINTDIST) THEN
         !             RFTMP(:) = RF(:)
         !             IF (PTEST) print '(2i5,a20,2f12.7)',A1,A2, 'keep permutation',INTDIST, CURINTDIST
         !             CURINTDIST = INTDIST
         !          ENDIF
         !       ENDDO
         !    ENDDO
         ! ENDDO
       !   RF(:) = RFTMP(:)
        
          RFTMP(:)= RFTMP2(:)
          MININTDIST=CURINTDIST
          MINSWAP(1)=A1; MINSWAP(2)=A2; MINSWAP(3)=A3 
          TMP_AR(A1,:)= RFTMP2(3*(A1-1)+1:3*A1)
          TMP_AR(A2,:)= RFTMP2(3*(A2-1)+1:3*A2)
          TMP_AR(A3,:)= RFTMP2(3*(A3-1)+1:3*A3)
          DO I=0,2
             RFTMP(3*(A1-1)+1:3*A1)=TMP_AR(A1,:)
             RFTMP(3*(A2-1)+1:3*A2)=TMP_AR(A2,:)
             RFTMP(3*(A3-1)+1:3*A3)=TMP_AR(A3,:)!initialise back to beginning
             !swap first two atoms
             B1=PERMGROUP(START+I)
             !find other two atoms
             IF (I.EQ.0) THEN !i.e. B1 = A1
                B2=PERMGROUP(START+1); B3= PERMGROUP(START+2)
             ELSEIF (I.EQ.1) THEN 
                B2=A1; B3= PERMGROUP(START+2) 
             ELSE
                B2=PERMGROUP(START+1); B3=A1
             ENDIF
         
! msb50 - difference in the way Lena prints her permutation and I do 
! when it prints here: B1, B2, B3 I mean the coordinates of original atom
! B1 are now on the first place etc.
! Lena means, the B1 has now the coordinates of whatever used to
! be B1 originally. 
! Therefore, I always swap A2 and A3 in DISTANCEPAIRSWAP, not B2 and B3
! as whatever is printed as B2, even though on the second place here, could
! be the first or third atoms
! otherwise, the way the coordinates are written is wrong

             RFTMP(3*(A1-1)+1:3*A1)=TMP_AR(B1,:)
             RFTMP(3*(B1-1)+1:3*B1)=TMP_AR(A1,:)
             PREVDIH(:) = STARTDIH(:)
             CALL GETDIHONLY(RFTMP)

             INTDIFF= PREVDIH(:) - STARTDIH(:)
             INTDIST= SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))
             IF (PTEST) PRINT*,"swap", B1,B2,B3 ,"intdist", intdist

             PREVDIH(:) = STARTDIH(:) 
             CALL DISTANCEPAIRSWAP(STARTDIH, RFTMP, A2,A3,PG, INTDIST2,PTEST,RFNEW)
             !PRINT*, "minintdist", MININTDIST
             IF (PTEST) PRINT*, "swap",B1, B3,B2,"intdist2", intdist2
             
             IF (INTDIST.LT.MININTDIST) THEN 
                 MININTDIST=INTDIST; MINSWAP(1)=B1;MINSWAP(2)=B2;MINSWAP(3)=B3
                 IF (PTEST) PRINT*, MINSWAP(:), MININTDIST
             ENDIF
             IF (INTDIST2.LT.MININTDIST) THEN
                 MININTDIST=INTDIST2; MINSWAP(1)=B1;MINSWAP(2)=B3;MINSWAP(3)=B2
                 IF (PTEST) PRINT*, MINSWAP(:), MININTDIST
             ENDIF
          ENDDO
          !PRINT*,A1,A2,A3, "minswap",MINSWAP(1), MINSWAP(2), MINSWAP(3)
          RFTMP(3*(A1-1)+1:3*A1)=RFTMP2(3*(MINSWAP(1)-1)+1:3*MINSWAP(1))
          RFTMP(3*(A2-1)+1:3*A2)=RFTMP2(3*(MINSWAP(2)-1)+1:3*MINSWAP(2))
          RFTMP(3*(A3-1)+1:3*A3)=RFTMP2(3*(MINSWAP(3)-1)+1:3*MINSWAP(3))
          
          CURINTDIST = MININTDIST   
          !  DO I=1, 3*NATOMS
          !    IF (RF(I).NE.RFTMP(I)) THEN 
          !     PRINT*, I, "coords not equal", RF(I), RFTMP(I)
          !     !STOP
          !    ENDIF
          !  ENDDO
          RF(:)=RFTMP(:)
      ENDIF !line 460 (473 atm) IF (PERMCHAIN...)
      
      ELSE
        print*, 'Error! INTMINPERM is only set up for permutation groups &
            &          of size 2 & 3 for now', PG, NPERMSIZE(PG)
        STOP
      ENDIF
      
      START = START + NPERMSIZE(PG)
      END DO
      
      IF (PTEST) print*, 'Final int distance: ', CURINTDIST    
      
      !PRINT*, "msb50 test"
      !RF(3*28+1:3*29)=RFTMP(3*29+1:3*30)
      !RF(3*29+1:3*30)=RFTMP(3*28+1:3*29)
      DISTANCE = DOT_PRODUCT(RF-RS,RF-RS)
      IF (PTEST) PRINT*, 'Final cart distance', DISTANCE    
 
      !STOP
      RETURN
      END SUBROUTINE INTMINPERM
      
   SUBROUTINE CART2INT(XCART, XINT) 
     ! converts cartesians to internals
     ! using parameters from intcommons
     ! no derivatives involved
     USE KEY, ONLY  : AMBERT, NABT
     
     IMPLICIT NONE
     DOUBLE PRECISION, INTENT(IN) :: XCART(3*NATOMS)
     DOUBLE PRECISION, INTENT(OUT) :: XINT(NINTC)
     DOUBLE PRECISION :: GDUMMY(3*NATOMS), GINT(NINTC)
     
          IF (AMBERT.OR.NABT) THEN
          CALL AMBTRANSFORM(XCART,GDUMMY,XINT,GINT,NINTC,3*NATOMS,NNZ,.FALSE.,.TRUE.,KD,INTEPSILON)
          ELSE
          CALL TRANSFORM(XCART,GDUMMY,XINT,GINT,NINTC,3*NATOMS,NNZ,.FALSE.,.TRUE.,KD,INTEPSILON)
          ENDIF
         
          RETURN
   END SUBROUTINE CART2INT
    

! ************************************************************* 
     SUBROUTINE INTSETUP
         USE MODAMBER9
         USE KEY, ONLY : AMBERT, NABT
         USE COMMONS, ONLY : DEBUG
         IMPLICIT NONE
         ! set up some of the initial global arrays and constants for working with internals
         ! interface to CHARMM    
         
         ALLOCATE(USECART(NATOMS))
         USECART(:) = .FALSE.
         
         IF (CHRMMT) THEN
            CALL GETCHRMINTPARAM
            ALLOCATE(BACKBONE(NATOMS))
            CALL GETBACKBONE
         ENDIF
         ! NINTS is the number of internals in commons file; NINTC should be same number in intcommons file
         IF (NATINT) THEN
           IF (AMBERT.OR.NABT) THEN
             CALL AMB_NATINTSETUP 
           ELSE 
             CALL NATINTSETUP
           ENDIF
           ALLOCATE(BACKBONE(NATOMS)) 
           CALL SETCARTATMS
           IF (USEPARFILE) THEN
             CALL GETNATINTERNFILE
         ELSE
           IF (AMBERT.OR.NABT) THEN
              IF (.NOT.ALLOCATED(PERMNEIGHBOURS)) CALL AMBGETNATINTERN
              IF (DEBUG) PRINT*, "intcoords> after ambgetnatintern"
           ELSE
               CALL GETNATINTERN
           ENDIF
         ENDIF
         NINTS = NINTC
         
         ALLOCATE(DIHINFOSINGLE(NDIH))
         ALIGNDIR = .FALSE.
         DIHINFOSINGLE(:) = 0.0D0
         PREVDIH => DIHINFOSINGLE(:)
         ENDIF
         
         IF( AMBERT.OR.NABT) THEN
            IF (NATINT) THEN 
              CALL AMBGETNNZNAT(NNZ)
              CALL AMB_GETKDNAT(KD)
            ELSE
              PRINT*, "error in intcoords - natint not set for AMBER"
              STOP
            ENDIF
         ELSE
            CALL GETNNZ(NNZ)
            CALL GETKD(KD)
         ENDIF
         BACKTCUTOFF = MINBACKTCUT
         
         CALL KEYINTPRINT
      END SUBROUTINE INTSETUP
         
      SUBROUTINE SETCARTATMS
      ! pick out the atoms for which to use cartesian coordinates
      ! make sure NATINTSETUP is done before this!
      
      IMPLICIT NONE
      INTEGER :: A
     
         NCRT = 0
         IF (.NOT. ALLOCATED(CARTATMS)) ALLOCATE(CARTATMS(NATOMS))
         
         DO A = 1,NATOMS
         IF (A.GE.CARTATMSTART.OR.(CHRMMT.AND.BBCART.AND.BACKBONE(A))) THEN
           NCRT = NCRT + 1
           CARTATMS(NCRT) = A
           USECART(A) = .TRUE.
         ENDIF       
         ENDDO
         
     END SUBROUTINE SETCARTATMS
     
     SUBROUTINE NATINTSETUP
     ! setup up various arrays necessary for natural internals
     ! this is the routine that defines the natural internal coordinates in terms of angle and torsion coefficients
     
     IMPLICIT NONE
     INTEGER, PARAMETER :: LARGENEG = -99999999
     DOUBLE PRECISION :: S6R, S2R, S26R, S18R, S12R, A, B, SABR, SABR2
     
     
     PRINT*, "natintsetup, NATOMS", NATOMS
     NBDS = 0; NCNT = 0; NRNG = 0; NFRG = 0; NLDH = 0; NIMP = 0
     NCRT = 0; NDIH = 0
     
     ALLOCATE(CENTERS(NATOMS,0:5), RINGS(TOTRES,6), FRINGS(TOTRES,6), LINDIH(NATOMS,-1:8), IMPDIH(NATOMS,4), CARTATMS(NATOMS))
     
     IF (CARTRESSTART.GT.TOTRES+1.OR.CARTRESSTART.LT.0) THEN
        CARTRESSTART = TOTRES + 1
        CARTATMSTART = NATOMS + 1
     ELSE
        CARTATMSTART = RESSTARTS(CARTRESSTART)+1
     ENDIF
 
     !     Set up COEFF table
     S6R = 1/SQRT(6.0)
     S2R = 1/SQRT(2.0)
     S26R = 1/SQRT(26.0)
     S18R = 1/SQRT(18.0)
     S12R = 1/SQRT(12.0)
     
     COEFF(1,1,1:6) = (/-1,-1,-1,1,1,1/)*S6R
     COEFF(1,2,1:6) = (/0,0,0,2,-1,-1/)*S6R
     COEFF(1,3,1:6) = (/0,0,0,0,1,-1/)*S2R
     COEFF(1,4,1:6) = (/2,-1,-1,0,0,0/)*S6R
     COEFF(1,5,1:6) = (/0,1,-1,0,0,0/)*S2R     
     
     COEFF(2,1,1:6) = (/0,1,-1,1,-1,0/)*0.5
     COEFF(2,2,1:6) = (/0,-1,-1,1,1,0/)*0.5
     COEFF(2,3,1:6) = (/0,-1,1,1,-1,0/)*0.5
     COEFF(2,4,1:6) = (/1,0,0,0,0,5/)*S26R
     COEFF(2,5,1:6) = (/5,0,0,0,0,1/)*S26R
     
     COEFF(3,1,1:6) = (/0,0,2,0,-1,-1/)*S6R
     COEFF(3,2,1:6) = (/0,0,0,0,1,-1/)*S2R
     COEFF(3,3,1:6) = (/1,4,0,1,0,0/)*S18R
     COEFF(3,4,1:6) = (/1,1,0,4,0,0/)*S18R
     COEFF(3,5,1:6) = (/4,1,0,1,0,0/)*S18R
     
     COEFF(4,1,1:3)=(/-1,-1,2/)*S6R
     COEFF(4,2,1:3)=(/1,-1,0/)*S2R

     COEFF(5,1,1:3)=(/0,1,-1/)*S2R
     COEFF(5,2,1:3)=(/2,-1,-1/)*S6R
     
     COEFF(6,1,1:1) = (/1/)
     
     COEFF(7,1,1:2) = (/1,-1/)*S2R
     
     COEFF(8,1,1:5) = (/1,-1,1,-1,0/)*0.5
     COEFF(8,2,1:5) = (/-1,-1,1,1,0/)*0.5
     COEFF(8,3,1:5) = (/-1,1,1,-1,0/)*0.5
     COEFF(8,4,1:5) = (/0,0,0,0,1/)
     
     A = COS(144.0*PI/180.0)
     B = COS(72.0*PI/180.0)
     SABR = 1/SQRT(1+2*A*A+2*B*B)
     SABR2 = 1/SQRT(2*(A-B)**2+2*(1-A)**2)
     COEFF(25,1,1:5) = (/1.0D0,A,B,B,A/)*SABR
     COEFF(25,2,1:5) = (/0.0D0,A-B,1.0D0-A,A-1.0D0,B-A/)*SABR2
     COEFF(25,3,1:5) = (/B,A,1.0D0,A,B/)*SABR
     COEFF(25,4,1:5) = (/A-1.0D0,B-A,0.0D0,A-B,1.0D0-A/)*SABR2
     
     COEFF(26,1,1:6) = (/1,-1,1,-1,1,-1/)*S6R
     COEFF(26,2,1:6) = (/2,-1,-1,2,-1,-1/)*S6R
     COEFF(26,3,1:6) = (/0,1,-1,0,1,-1/)*0.5D0
     COEFF(26,4,1:6) = (/1,-1,1,-1,1,-1/)*S6R
     COEFF(26,5,1:6) = (/1,0,-1,1,0,-1/)*0.5D0
     COEFF(26,6,1:6) = (/-1,2,-1,-1,2,-1/)*S12R

!     Set up ATOMxPOSINC table
!     large negative number listed for those atoms not involved
     ATOMxPOSinC(1,1,1:5) = (/0,3,6,9,12/)
     ATOMxPOSinC(1,2,1:5) = (/15,LARGENEG,18,21,24/)
     ATOMxPOSinC(1,3,1:5) = (/27,LARGENEG,30,33,36/)
     ATOMxPOSinC(1,4,1:5) = (/39,42,45,48,51/)
     ATOMxPOSinC(1,5,1:5) = (/54,57,LARGENEG,60,63/)
     
     ATOMxPOSinC(2,1,1:5) = (/0,3,6,9,12/)
     ATOMxPOSinC(2,2,1:5) = (/15,18,21,24,27/)
     ATOMxPOSinC(2,3,1:5) = (/30,33,36,39,42/)
     ATOMxPOSinC(2,4,1:5) = (/45,48,51,54,57/)
     ATOMxPOSinC(2,5,1:5) = (/60,63,66,69,72/)
     
     ATOMxPOSinC(3,1,1:5) = (/0,3,6,9,12/)
     ATOMxPOSinC(3,2,1:5) = (/15,LARGENEG,18,21,24/)
     ATOMxPOSinC(3,3,1:5) = (/27,30,33,36,LARGENEG/)
     ATOMxPOSinC(3,4,1:5) = (/39,42,45,48,LARGENEG/)
     ATOmxPOSinC(3,5,1:5) = (/51,54,57,60,LARGENEG/)
     
     ATOMxPOSinC(4,1,1:4)=(/0,3,6,9/)
     ATOMxPOSinC(4,2,1:4)=(/12,15,18,21/)
     
     ATOMxPOSinC(5,1,1:4)=(/0,3,6,9/)
     ATOMxPOSinC(5,2,1:4)=(/12,15,18,21/)
     
     ATOMxPOSinC(6,1,1:3) = (/0,3,6/)
     
     ATOMxPOSinC(7,1,1:4) = (/0,3,6,9/)      

     ATOMxPOSinC(8,1,1:4)=(/0,3,6,9/)
     ATOMxPOSinC(8,2,1:4)=(/12,15,18,21/)
     ATOMxPOSinC(8,3,1:4)=(/24,27,30,33/)
    END SUBROUTINE NATINTSETUP

    SUBROUTINE GETBACKBONE
    ! which of the atoms are part of the backbone?
    
    IMPLICIT NONE
    
    INTEGER :: A, R
    
      BACKBONE(1:NATOMS) = .FALSE.
      DO A = 1,NATOMS
        IF (ATMTYPES(A).EQ.'C'.OR.ATMTYPES(A).EQ.'CA'.OR.ATMTYPES(A).EQ.'N') BACKBONE(A) = .TRUE.
      ENDDO
      
      ! label all proline atoms (except O) as backbone
      DO R = 1,TOTRES
        IF (RESLIST(R).EQ.'PRO') THEN
          DO A = RESSTARTS(R)+1,RESSTARTS(R+1)
             IF (ATMTYPES(A)(1:1).NE.'O') BACKBONE(A) = .TRUE.
          ENDDO
        ENDIF
      ENDDO
      RETURN
      
   END SUBROUTINE GETBACKBONE
    
   SUBROUTINE INTCLEANUP
    ! cleanup global arrays for internals
   IMPLICIT NONE

    DEALLOCATE(USECART)
    IF (CHRMMT) DEALLOCATE(BACKBONE, RESLIST, RESSTARTS, ATMTYPES)
    IF (NATINT) THEN
    DEALLOCATE(DIHINFOSINGLE,CENTERS, RINGS, FRINGS, LINDIH, IMPDIH, CARTATMS)
    NULLIFY(PREVDIH)
    ENDIF
    
   END SUBROUTINE INTCLEANUP
    
   SUBROUTINE KEYINTPRINT
    ! print out key information related to internals
    
   IMPLICIT NONE
    
    IF (NATINT) THEN
      WRITE(*,'(1x,a,I10)') 'KeyInt>> Using natural internal coordinates. # of coords = ', NINTC
      WRITE(*,'(1x,a,a)') 'KeyInt>> Parameter file for internal coordinates is ', INTPARFILE
      WRITE(*,'(1x,a,2I10)') 'KeyInt>> KD, NNZ are ', KD, NNZ
    ELSE
      WRITE(*,'(1x,a,I10)') 'KeyInt>> Using primitive internal coordinates. # of coords = ', NINTC
    ENDIF
    
    IF (INTINTERPT.OR.DESMINT) THEN
      WRITE(*,'(1x,a,G20.10)') 'KeyInt>> Back transform convergence (for interpolation) ', INTERPBACKTCUT
    ENDIF
    WRITE(*,'(1x,a,G20.10)') &
     & 'KeyInt>> Back transform convergence criterion (except for interpolation) ', MINBACKTCUT
    
    IF (BBCART) WRITE(*,'(1x,a)') 'KeyInt>> Using cartesian coordinates for the backbone and prolines' 
    
    IF (CARTRESSTART.LE.TOTRES) &
     &   WRITE(*,'(1x,a,I6)') 'KeyInt>> using Cartesian coords for all residues starting with ', CARTRESSTART
    
    IF (.NOT.INTNEWT) WRITE(*,'(1x,a)') "KeyInt>> Linear back-transforms only. No Newton's iterations"
    
    IF (DESMINT) WRITE(*,'(1x,a)') 'KeyInt>> Using internal coords for double-ended search methods'
    IF (INTINTERPT) THEN 
       WRITE(*,'(1x,a)') 'KeyInt>> Using internal coordinates for interpolation.'
       IF (INTERPCHOICE) WRITE(*,'(1x,a)') &
     &   'KeyInt>> Choosing internal or Cartesian coordinates for interpolation, based on max energy.'
       IF (INTERPSIMPLE) WRITE(*,'(1x,a)') &
      &   'KeyInt>> using simple internal interpolation. &
      &   Image points will not be evenly distributed in Cartesians'
    ENDIF
    
    IF (INTMINPERMT) WRITE(*,'(1x,a)') 'KeyInt>> Permuting endpoint atoms to minimise torsion distance'
    
    IF (PRINTCOORDS) WRITE(*,'(1x,a)') 'KeyInt>> Printing coordinate information only.'

    RETURN
    END SUBROUTINE KEYINTPRINT


! ********************************************************************************
! msb50 to help with intminperm **************************************************

    SUBROUTINE DISTANCEPAIRSWAP(STARTDIH, RF, A1, A2, PG, INTDIST, PTEST,RFNEW)
    !gives you internal distance after swapping two atoms
    !NOTE - when swapping NPERMSIZE(PG).EQ.4 need to call this 6 times
        use key, only: NSETS, SETS
        use commons, only: NATOMS
        use intcommons, only: NDIH
        IMPLICIT NONE
         
        DOUBLE PRECISION, INTENT(IN) :: STARTDIH(NDIH)
        DOUBLE PRECISION, INTENT(IN) :: RF(NATOMS*3)
        INTEGER, INTENT(IN) :: PG !which permgroup (number)
        INTEGER, INTENT(IN) :: A1, A2 ! which atoms
        LOGICAL, INTENT(IN) :: PTEST
        DOUBLE PRECISION, INTENT(OUT) :: INTDIST
        DOUBLE PRECISION, INTENT(OUT) ::RFNEW(NATOMS*3)
        INTEGER :: A3, A4
        INTEGER :: SW
        DOUBLE PRECISION :: INTDIFF(NDIH), GD_PREVDIH(NDIH)
        DOUBLE PRECISION :: RFTMP(NATOMS*3)
        !note rftmp and rf wrong way round as rf intent(in), cannot change
        
          RFTMP(:) = RF(:)
          RFTMP(3*(A1-1)+1:3*A1) = RF(3*(A2-1)+1:3*A2)
          RFTMP(3*(A2-1)+1:3*A2) = RF(3*(A1-1)+1:3*A1)
          DO SW = 1,NSETS(PG)
!            A3 = SWAP1(PG,SW)
!            A4 = SWAP2(PG,SW)
             A3 = SETS(A1,SW)
             A4 = SETS(A2,SW)
        
             IF (PTEST) print*, 'Dragging swap: ', A3, A4
        
             RFTMP(3*(A3-1)+1:3*A3) = RF(3*(A4-1)+1:3*A4)
             RFTMP(3*(A4-1)+1:3*A4) = RF(3*(A3-1)+1:3*A3)
          ENDDO
          PREVDIH(:) = STARTDIH(:)
          CALL GETDIHONLY(RFTMP) !msb50
          GD_PREVDIH(:) = PREVDIH(:)
          PREVDIH(:) = STARTDIH(:)

          INTDIFF(:) = GD_PREVDIH(:) - STARTDIH(:)
          INTDIST = SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))

          RFNEW(:) = RFTMP(:)

       END SUBROUTINE DISTANCEPAIRSWAP

!*************************************************************************
!msb50 this subroutine is for intminperm, where we had the problem that if your metric
! is based on dihedrals the dihedrals must be independent. In case of connected (neighboured)
! permgroups, this is not the case. So here is the function if you have 2 neighbouring
! permgroups
      SUBROUTINE SWAP2ATONCE(STARTDIH,RFTMP,PERMG,START,PG2DUMMY,PTEST,&
     & SKIPPG,RFNEW,MINSWAP1,MINSWAP2,MININTDIST, SWAPARRAY, SWAPDIST)

         USE KEY, ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP,NSETS, SETS
!         USE MODAMBER9, ONLY : ih, ix, NRES, i02,m04,m02
         USE intcommons, only: NDIH, PERMNEIGHBOURS ,PREVDIH
         USE commons, only: NATOMS !msb50
         USE porfuncs
         IMPLICIT NONE

         DOUBLE PRECISION, INTENT(IN) :: STARTDIH(NDIH)
         DOUBLE PRECISION, INTENT(IN) :: RFTMP(3*NATOMS)
         INTEGER, INTENT(IN) ::PERMG, START !where atoms are in pg array
         INTEGER, INTENT(IN) :: PG2DUMMY
         LOGICAL, INTENT(IN) :: PTEST
         LOGICAL, INTENT(INOUT) :: SKIPPG(NPERMGROUP)
         INTEGER :: ISTAT
         INTEGER, INTENT(OUT) :: MINSWAP1(3), MINSWAP2(3)
         DOUBLE PRECISION, INTENT(OUT) :: MININTDIST
         
         INTEGER, INTENT(OUT) :: SWAPARRAY(12,6) !this thing stores possible atom combinations
         DOUBLE PRECISION, INTENT(OUT) :: SWAPDIST(12) !and their distances 

         INTEGER :: MINSWAP(2,3)!-order of permable atoms in which you want them
         DOUBLE PRECISION :: RFNEW(3*NATOMS),  INTDIFF(NDIH)
         DOUBLE PRECISION :: RFNEW2(3*NATOMS),RF2(3*NATOMS)
         DOUBLE PRECISION :: INTDIST2, INTDIST
         DOUBLE PRECISION :: TMP_AR(NATOMS,3)
         INTEGER:: A1, A2, A3, B1, B2, B3, START2,pg1,pg2, C1,C2,C3,ST,AN1,AN2
         INTEGER :: PERMG2, JJ,I, indx, A4
         INTEGER:: II, START_AR(2), swI, SW
         
         
         A1=0;A2=0;A3=0;B1=0;B2=0;B3=0
         IF (PG2DUMMY.EQ.0) THEN !we call it directly -only has 1neighbour
             PERMG2=PERMNEIGHBOURS(PERMG,2)
         ELSE
             PERMG2=PG2DUMMY
         ENDIF
         IF (SKIPPG(PERMG).AND.SKIPPG(PERMG2)) THEN
           IF (PTEST) PRINT*, "swap2atonce: groups already done:", PERMG, PERMG2
           RETURN
         ENDIF
         START2= START
         RF2(:)=RFTMP(:)
         RFNEW(:)=RFTMP(:)
         RFNEW2(:)=RFTMP(:)
         SWAPARRAY(:,:)=0
         SWAPDIST(:)=10000.0
   !      PRINT*, "this was at the begininng"
         IF (PERMG2.LT.PERMG) THEN
            PRINT*, "PERMG2 LT PERMG", PERMG2, "PERMG", PERMG
            STOP
         ENDIF
   !      PRINT*, "swap2, pg, pg2", PERMG,PERMG2
         DO II=PERMG,NPERMGROUP
          IF (II.EQ.PERMG2) EXIT
          START2 = START2 + NPERMSIZE(II)
         ENDDO
         START_AR(1)=START; START_AR(2)=START2
         
         PREVDIH(:)=STARTDIH(:)
         CALL GETDIHONLY(RFTMP)
         INTDIFF(:)= PREVDIH(:)-STARTDIH(:)      
         MININTDIST=SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))
         
      !the idea of this is that you start with a pg that has size 2
      !but you still want the swaparray in the initial order   
      DO II=1,2
       IF (II.EQ.1) THEN
          pg1=PERMG; pg2=PERMG2
       ELSE 
          pg1=PERMG2;pg2=PERMG
          IF (NPERMSIZE(pg2).EQ.2) CYCLE
       ENDIF
       IF (NPERMSIZE(pg1).EQ.3) THEN
           CYCLE !IF both 3 go for sth else
       ELSE
           A1 = PERMGROUP(START_AR(II))
           A2 = PERMGROUP(START_AR(II)+1)
   !        print*, II, "pg",pg1,A1,a2
           MINSWAP(II,1)=A1;MINSWAP(II,2)=A2;MINSWAP(II,3)=0
           SWAPDIST(1)=MININTDIST
           ! start with 1 if II==1, with 4 if II==2
           SWAPARRAY(1,3*(II-1)+1:3*(II-1)+3)=MINSWAP(II,:)
           SWAPARRAY(2,3*(II-1)+1:3*(II-1)+3)=MINSWAP(II,:)
           SWAPARRAY(3,3*(II-1)+1:3*(II-1)+3)=MINSWAP(II,:)
           SWAPARRAY(4,3*(II-1)+1:3*(II-1)+3)=MINSWAP(II,:)
           SWAPARRAY(5,3*(II-1)+1:3*(II-1)+3)=MINSWAP(II,:)
           SWAPARRAY(6,3*(II-1)+1:3*(II-1)+3)=MINSWAP(II,:)
           DO JJ=1,2  
              indx=MOD(II,2)+1
              !indx: 2 if II=1 (i.e. pg2=PERMG2, want A first)
              !      1 if II=2 (i.e. pg2=PERMG, want PERMG first still -so B first)
              IF (NPERMSIZE(pg2).EQ.2) THEN
                  C1=PERMGROUP(START_AR(indx))
                  C2=PERMGROUP(START_AR(indx)+1)
                  B1=C1
                  B2=C2
   !               print*, II, "pg2", pg2, B1,B2
                  IF (JJ.EQ.1) THEN !initialise minswap
                      MINSWAP(indx,1)=B1
                      MINSWAP(indx,2)=B2
                      MINSWAP(indx,3)=0
                  SWAPARRAY(1,3*(indx-1)+1:3*indx)=MINSWAP(indx,:)
                  ENDIF
                  PREVDIH(:) = STARTDIH(:) 
                  CALL DISTANCEPAIRSWAP(STARTDIH, RF2, B1,B2,pg2, &
&                          INTDIST,PTEST,RFNEW2)
                 !RFNEW2 going to be overwritten, don't care
                  SWAPARRAY(6*(JJ-1)+2,3*(indx-1)+1)=B2
                  SWAPARRAY(6*(JJ-1)+2,3*(indx-1)+2)=B1
                  SWAPDIST(6*(JJ-1)+2)=INTDIST
                  IF (INTDIST.LT.MININTDIST) THEN
                     MININTDIST=INTDIST; MINSWAP(indx,1)=B2
                     MINSWAP(indx,2)=B1
  !                   PRINT*, "swapping", indx, MINSWAP(indx,1)
                  ENDIF
  !                PRINT*, "minswap indx", indx, MINSWAP(indx,:)
                  IF (JJ.EQ.1) THEN
                  ! as coords for next run are still for this arrangement
                     SWAPARRAY(7,3*(indx-1)+1)=B1
                     SWAPARRAY(7,3*(indx-1)+2)=B2
                     SWAPARRAY(7,3*(indx-1)+3)=0
                  ELSE
                    CYCLE
                  ENDIF
                  
                  CALL FLUSH(6)         
               ELSEIF (NPERMSIZE(pg2).EQ.3) THEN
                  C1=PERMGROUP(START_AR(indx))
                  C2= PERMGROUP(START_AR(indx)+1)
                  C3=PERMGROUP(START_AR(indx)+2)
                  IF (JJ.EQ.1) THEN !initialise minswap
                         MINSWAP(indx,1)=C1
                         MINSWAP(indx,2)=C2
                         MINSWAP(indx,3)=C3
                  SWAPARRAY(1,3*(indx-1)+1:3*indx)=MINSWAP(indx,:) 
                  ENDIF
                  TMP_AR(C1,:)= RF2(3*(C1-1)+1:3*C1)
                  TMP_AR(C2,:)= RF2(3*(C2-1)+1:3*C2)
                  TMP_AR(C3,:)= RF2(3*(C3-1)+1:3*C3)
                  ST=START_AR(indx)
                  DO I=0,2
                     RF2(3*(C1-1)+1:3*C1)=TMP_AR(C1,:)
                     RF2(3*(C2-1)+1:3*C2)=TMP_AR(C2,:)
                     RF2(3*(C3-1)+1:3*C3)=TMP_AR(C3,:)!initialise back to beginning
                     !swap first two atoms
                     B1=PERMGROUP(ST+I)
                     !find other two atoms
                     IF (I.EQ.0) THEN !i.e. B1 = A1
                        B2=PERMGROUP(ST+1); B3= PERMGROUP(ST+2)
                     ELSEIF (I.EQ.1) THEN
                        B2=C1; B3= PERMGROUP(ST+2)
                     ELSE
                        B2=PERMGROUP(ST+1); B3=C1
                     ENDIF
             
                     RF2(3*(C1-1)+1:3*C1)=TMP_AR(B1,:)
                     RF2(3*(B1-1)+1:3*B1)=TMP_AR(C1,:)
                     PREVDIH(:) = STARTDIH(:)
                     CALL GETDIHONLY(RF2)
                     INTDIFF= PREVDIH(:) - STARTDIH(:)
                     INTDIST= SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))
                     swI=6*(JJ-1)+2*I
                     SWAPARRAY(swI+1,3*(indx-1)+1)=B1
                     SWAPARRAY(swI+1,3*(indx-1)+2)=B2
                     SWAPARRAY(swI+1,3*(indx-1)+3)=B3
                     SWAPDIST(swI+1)=INTDIST
                     PREVDIH(:) = STARTDIH(:)
!note always swap last two, not B2 and B3!!! check intminperm for comment
                     CALL DISTANCEPAIRSWAP(STARTDIH,RF2, C2,C3,pg2, &
     &                      INTDIST2,PTEST,RFNEW2)
                     !RFNEW2- going to be overwritten, don't care
                     RF2(:)=RFNEW(:) 
                     !coordinates have to be same for next loop!
!                     PRINT*, "swap",B1, B3,B2,"intdist2", intdist2
                     SWAPARRAY(swI+2,3*(indx-1)+1)=B1
                     SWAPARRAY(swI+2,3*(indx-1)+2)=B3
                     SWAPARRAY(swI+2,3*(indx-1)+3)=B2
                     SWAPDIST(swI+2)=INTDIST2
    !                 PRINT*, "swI=2",swI+2
    !                 PRINT '(a10,6i5,f11.7)',"swaparry",SWAPARRAY(swI+2,1),&
    ! &SWAPARRAY(swI+2,2),SWAPARRAY(swI+2,3),SWAPARRAY(swI+2,4),&
    ! &SWAPARRAY(swI+2,5),SWAPARRAY(swI+2,6), SWAPDIST(swI+2)
                     IF (INTDIST.LT.MININTDIST) THEN
                         MININTDIST=INTDIST
                         MINSWAP(indx,1)=B1
                         MINSWAP(indx,2)=B2;MINSWAP(indx,3)=B3
                     ENDIF
                     IF (INTDIST2.LT.MININTDIST) THEN
                         MININTDIST=INTDIST2
                         MINSWAP(indx,1)=B1
                         MINSWAP(indx,2)=B3;MINSWAP(indx,3)=B2
                     ENDIF
                  ENDDO
                  !because the coordinates were not changed, so we're back to this
                  IF (JJ.EQ.1) THEN 
                    SWAPARRAY(7,3*(indx-1)+1)=C1
                    SWAPARRAY(7,3*(indx-1)+2)=C2
                    SWAPARRAY(7,3*(indx-1)+3)=C3
                  ELSE 
                    CYCLE
                  ENDIF
              ENDIF
              PREVDIH(:) = STARTDIH(:)
              CALL DISTANCEPAIRSWAP(STARTDIH, RF2, A1,A2,pg1,INTDIST,&
     &              PTEST,RFNEW)
!              PRINT*, "swap", A1, A2, "dist", INTDIST
    !          PRINT*, "end, swaparray"
              SWAPARRAY(7,3*(II-1)+1)=A2;SWAPARRAY(7,3*(II-1)+2)=A1
              SWAPARRAY(8,3*(II-1)+1)=A2;SWAPARRAY(8,3*(II-1)+2)=A1
              SWAPARRAY(9,3*(II-1)+1)=A2;SWAPARRAY(9,3*(II-1)+2)=A1
              SWAPARRAY(10,3*(II-1)+1)=A2;SWAPARRAY(10,3*(II-1)+2)=A1
              SWAPARRAY(11,3*(II-1)+1)=A2;SWAPARRAY(11,3*(II-1)+2)=A1
              SWAPARRAY(12,3*(II-1)+1)=A2;SWAPARRAY(12,3*(II-1)+2)=A1
              SWAPDIST(7)=INTDIST
              RF2(:)=RFNEW(:)
     
  !            PRINT*, "minswap", minswap(1,:)
  !            PRINT*, "minswap", minswap(2,:)         
              IF (INTDIST.LT.MININTDIST) THEN 
                  IF (NPERMSIZE(pg2) == 2) THEN 
                     MINSWAP(indx, 1) = B1; MINSWAP(indx,2) = B2
                  ELSE
                     MINSWAP(indx,1) = C1; MINSWAP(indx,2)=C2; MINSWAP(indx,3)=C3
                  ENDIF
                  MININTDIST=INTDIST;MINSWAP(II,1)=A2; MINSWAP(II,2)=A1
              ENDIF
           ENDDO !jj - now do with changed RF2   
         ENDIF 
      ENDDO     

      IF (NPERMSIZE(PERMG).EQ.3.AND.NPERMSIZE(PERMG2).EQ.3) THEN
         PRINT*, "currently no case known where two permgroups&
     &      are neighbours with NPERMSIZE3 - have fun coding this"
         STOP
      ENDIF   
     
      !PRINT*, "in swap2atonce, final array"
      !DO I=1,12
      !  PRINT*, "swaps"
      !  PRINT*, SWAPARRAY(I,:), SWAPDIST(I)
      !ENDDO

      IF (PTEST) PRINT*, "minintdist in swap2atonce", minintdist
      !bookkeeping
      SKIPPG(PERMG)=.TRUE.; SKIPPG(PERMG2)=.TRUE.
      IF (NPERMSIZE(PERMG).EQ.2) THEN
         MINSWAP1(:)=MINSWAP(1,:); MINSWAP1(3)=0
         MINSWAP2(:)=MINSWAP(2,:)
      ELSE
         MINSWAP1(:)=MINSWAP(2,:)
         MINSWAP2(:)=MINSWAP(1,:); MINSWAP2(3)=0
      ENDIF
      IF (PTEST) PRINT*, "minswap1", MINSWAP(1,:)
      IF (PTEST) PRINT*, "minswap2", MINSWAP(2,:)
      AN1=MINSWAP(1,1); AN2=MINSWAP(1,2)
!      PRINT*, "ans", an1, an2
      RFNEW(:)=RFTMP(:)
      RFNEW(3*(A1-1)+1:3*A1)= RFTMP(3*(AN1-1)+1:3*AN1)
      RFNEW(3*(A2-1)+1:3*A2)= RFTMP(3*(AN2-1)+1:3*AN2) 
      RFNEW(3*(C1-1)+1:3*C1)= RFTMP(3*(MINSWAP(2,1)-1)+1:3*MINSWAP(2,1))
      RFNEW(3*(C2-1)+1:3*C2)= RFTMP(3*(MINSWAP(2,2)-1)+1:3*MINSWAP(2,2))
      IF (MINSWAP(2,3).NE.0) THEN !ie not both pgsizes where 2
            RFNEW(3*(C3-1)+1:3*C3)= RFTMP(3*(MINSWAP(2,3)-1)+1:3*MINSWAP(2,3))
      ENDIF
      !note that I never have a constellation where npermsize(1stPG)=3 
      ! and nswap(2ndPG)!=0 as 3 and sth draggable are only coupled in 
      ! val and leu for which there is an exception
      ! Val and leu aren't even recognised as having permneighbouring atoms
      ! (apart from the other case in leucine)
      IF (NSETS(PERMG).NE.0) THEN!automatically npermsize(permg)=2 i.e.permg->A1
         IF (A1.NE.AN1) THEN
           DO SW = 1,NSETS(PERMG)
!              A3 = SWAP1(PERMG,SW)
!              A4 = SWAP2(PERMG,SW)
               A3 = SETS(AN1,SW) ! is AN1 the right argument under the new scheme? DJW
               A4 = SETS(AN2,SW) ! is AN2 the right argument under the new scheme? DJW
               IF (PTEST) PRINT*, "dragging swap in swap2atonce", A3, A4
               RFNEW(3*(A3-1)+1:3*A3) = RFTMP(3*(A4-1)+1:3*A4)
               RFNEW(3*(A4-1)+1:3*A4) = RFTMP(3*(A3-1)+1:3*A3)
            ENDDO
          ENDIF
      ELSEIF (NSETS(PERMG2).NE.0) THEN !both is never the case
          IF (C1.NE.MINSWAP(2,1)) THEN !C1 must correspond to 2nd pg then
            DO SW = 1,NSETS(PERMG)
!              A3 = SWAP1(PERMG,SW)
!              A4 = SWAP2(PERMG,SW)
               A3 = SETS(MINSWAP(2,1),SW) ! is MINSWAP(2,1) the right argument under the new scheme? DJW
               A4 = SETS(MINSWAP(2,2),SW) ! is MINSWAP(2,2) the right argument under the new scheme? DJW
               IF (PTEST) PRINT*, "dragging swap in swap2atonce", A3, A4
               RFNEW(3*(A3-1)+1:3*A3) = RFTMP(3*(A4-1)+1:3*A4)
               RFNEW(3*(A4-1)+1:3*A4) = RFTMP(3*(A3-1)+1:3*A3)
            ENDDO
          ENDIF
      ENDIF      
 
      END SUBROUTINE SWAP2ATONCE

!************************************************************************
      SUBROUTINE ARG_SWAP(STARTDIH, RF, PERMG, START,PTEST,SKIPPG,RFNEW,CURINTDIST)
      USE KEY, ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP
!      USE MODAMBER9, ONLY : ih, ix, NRES, i02,m04
      USE intcommons, only: NDIH, PERMNEIGHBOURS,PREVDIH
      USE commons, only: NATOMS !msb50
      USE porfuncs 
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: STARTDIH(NDIH)
      DOUBLE PRECISION, INTENT(IN) :: RF(3*NATOMS)
      INTEGER, INTENT(IN) ::PERMG, START !where atoms are in pg array
      DOUBLE PRECISION, INTENT(INOUT) :: CURINTDIST
      LOGICAL, INTENT(IN) :: PTEST
      LOGICAL, INTENT(INOUT) :: SKIPPG(NPERMGROUP)
      DOUBLE PRECISION, INTENT(OUT) ::RFNEW(3*NATOMS)
      DOUBLE PRECISION :: INTDIST2, INTDIST,CURINTDIST2
      DOUBLE PRECISION :: RF2(3*NATOMS), RF3(3*NATOMS)
      INTEGER:: PERMG2, PG3
      INTEGER:: II, START_AR(2), START2, START3
      INTEGER:: A1, A2, B1,B2, C1, C2
      

      A1=0;A2=0;B1=0;B2=0;C1=0;C2=0
      START2=START; START3=START
      CURINTDIST2=CURINTDIST
      PERMG2=PERMNEIGHBOURS(PERMG,2); PG3=PERMNEIGHBOURS(PERMG,3)
      IF (NPERMSIZE(PERMG).NE.2.AND.NPERMSIZE(PERMG2).NE.2) THEN
         PRINT*, "ARG_SWAP error - wrong permsize"
         STOP
      ENDIF
      RF2(:)=RF(:)
      RF3(:)=RF(:)
      DO II=PERMG,NPERMGROUP
        IF (II.EQ.PERMG2) EXIT
        START2 = START2 + NPERMSIZE(II)
      ENDDO
      DO II=PERMG,NPERMGROUP
        IF (II.EQ.PG3) EXIT
        START3 = START3 + NPERMSIZE(II)
      ENDDO

      A1=PERMGROUP(START); A2=PERMGROUP(START+1)
      B1=PERMGROUP(START2); B2=PERMGROUP(START2+1)
      C1=PERMGROUP(START3); C2=PERMGROUP(START3+1)
      !swap nitrogens and H's along with them
      CALL DISTANCEPAIRSWAP(STARTDIH,RF2,A1,A2,PERMG,INTDIST,PTEST,RFNEW)
      RF2(:)=RFNEW(:)
      CURINTDIST=INTDIST
      IF (PTEST) PRINT*, "swap", a1,a2,intdist
      !swap H's 
      CALL DISTANCEPAIRSWAP(STARTDIH,RF2,B1,B2,PERMG2,INTDIST,PTEST,RFNEW)
      IF (INTDIST.LT.CURINTDIST) THEN
         RF2(:)=RFNEW(:)
         CURINTDIST=INTDIST
      ENDIF
      IF (PTEST) PRINT*, "swap",b1,b2,intdist
      !swap H's
      CALL DISTANCEPAIRSWAP(STARTDIH,RF2,C1,C2,PG3,INTDIST,PTEST,RFNEW)
      IF (INTDIST.LT.CURINTDIST) THEN
         RF2(:)=RFNEW(:)
         CURINTDIST=INTDIST
      ENDIF
      IF (PTEST) PRINT*, "swap",c1,c2,intdist

      !leave N's and swap H's
      CALL DISTANCEPAIRSWAP(STARTDIH,RF3,B1,B2,PERMG2,INTDIST2,PTEST,RFNEW)
      IF (INTDIST2.LT.CURINTDIST2) THEN
         RF3(:)=RFNEW(:)
         CURINTDIST2=INTDIST2
      ENDIF
      IF (PTEST) PRINT*, "swap",b1,b2,intdist2  
      CALL DISTANCEPAIRSWAP(STARTDIH,RF3,C1,C2,PG3,INTDIST2,PTEST,RFNEW)
      IF (INTDIST2.LT.CURINTDIST2) THEN
         RF3(:)=RFNEW(:)
         CURINTDIST2=INTDIST2
      ENDIF
      IF (PTEST) PRINT*, "swap",c1,c2,intdist2
      IF (CURINTDIST2.LT.CURINTDIST) THEN
         RFNEW(:)=RF3(:) !which is still the same as the orig one if no changes were better
         CURINTDIST=CURINTDIST2
      ELSE
         RFNEW(:)=RF2(:)
      ENDIF
      SKIPPG(PERMG)=.TRUE.;SKIPPG(PERMG2)=.TRUE.;SKIPPG(PG3)=.TRUE.

      RETURN
      END SUBROUTINE ARG_SWAP

! *******************************************************************************
      SUBROUTINE LEU_SWAP(STARTDIH, RF, PERMG, START,PTEST,SKIPPG,RFNEW,CURINTDIST)
      USE KEY, ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP
!      USE MODAMBER9, ONLY : ih, ix, NRES, i02,m04
      USE intcommons, only: NDIH,PREVDIH
      USE commons, only: NATOMS !msb50
      USE porfuncs
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: STARTDIH(NDIH)
      DOUBLE PRECISION, INTENT(IN) :: RF(3*NATOMS)
      INTEGER, INTENT(IN) ::PERMG, START !where atoms are in pg array
      DOUBLE PRECISION, INTENT(INOUT) :: CURINTDIST
      LOGICAL, INTENT(IN) :: PTEST
      LOGICAL, INTENT(INOUT) :: SKIPPG(NPERMGROUP)
      DOUBLE PRECISION, INTENT(OUT) ::RFNEW(3*NATOMS)
      DOUBLE PRECISION :: INTDIST2, INTDIST,CURINTDIST2
      DOUBLE PRECISION :: RF2(3*NATOMS), RF3(3*NATOMS)
      INTEGER:: PERMG2, PG3
      INTEGER:: II,START_AR(2),  START2, START3,ST
      INTEGER:: A1, A2, C1, C2, C3,B1, B2, B3, JJ,I, pg2
      DOUBLE PRECISION :: INTDIFF(NDIH)   
      
      INTEGER :: MINSWAP(2,3)
      DOUBLE PRECISION :: MININTDIST, TMP_AR(NATOMS,4)

      INTEGER :: PERM_AR(2)
 
      
      A1=0;A2=0;B1=0;B2=0;C1=0;C2=0; C3=0;B3=0
      START2=START; START3=START
      CURINTDIST2=CURINTDIST
!      PERMG2=PERMNEIGHBOURS(PERMG,2); PG3=PERMNEIGHBOURS(PERMG,3)
      IF (NPERMSIZE(PERMG).NE.2.AND.NPERMSIZE(PERMG2).NE.2) THEN
         PRINT*, "ARG_SWAP error - wrong permsize"
         STOP
      ENDIF
      RFNEW(:)=RF(:)
      RF2(:)=RF(:)
      RF3(:)=RF(:)
!      DO II=PERMG,NPERMGROUP
!        IF (II.EQ.PERMG2) EXIT
!        START2 = START2 + NPERMSIZE(II)
!      ENDDO
!      DO II=PERMG,NPERMGROUP
!        IF (II.EQ.PG3) EXIT
!        START3 = START3 + NPERMSIZE(II)
!      ENDDO
      PERMG2=PERMG+1 
      PERM_AR(1)=PERMG2; PERM_AR(2)=PERMG+2
      PRINT*, "permg", PERMG2, PERM_AR(2)
      START2=START+2; START3=START2+3
      START_AR(1)=START2;START_AR(2)=START3
      A1=PERMGROUP(START); A2=PERMGROUP(START+1)
      PRINT*, "as", A1, A2
      PRINT*, "curintidst", curintdist
      MININTDIST=CURINTDIST

      !find minimal hydrogen alignment first

      DO JJ =1,2
         PG2= PERM_AR(JJ)
         C1=PERMGROUP(START_AR(JJ))
         C2= PERMGROUP(START_AR(JJ)+1)
         C3=PERMGROUP(START_AR(JJ)+2)
         IF (PTEST) PRINT*,pg2, "pg2", c1,c2, c3
         MINSWAP(JJ,1)=C1; MINSWAP(JJ,2)=C2; MINSWAP(JJ,3)=C3
         TMP_AR(C1,:)= RF2(3*(C1-1)+1:3*C1)
         TMP_AR(C2,:)= RF2(3*(C2-1)+1:3*C2)
         TMP_AR(C3,:)= RF2(3*(C3-1)+1:3*C3)
         ST=START_AR(JJ)
!        PRINT*, "sec group3, start", st
         DO I=0,2
            RF2(3*(C1-1)+1:3*C1)=TMP_AR(C1,:)
            RF2(3*(C2-1)+1:3*C2)=TMP_AR(C2,:)
            RF2(3*(C3-1)+1:3*C3)=TMP_AR(C3,:)!initialise back to beginning
            !swap first two atoms
            B1=PERMGROUP(ST+I)
            !find other two atoms
            IF (I.EQ.0) THEN !i.e. B1 = A1
               B2=PERMGROUP(ST+1); B3= PERMGROUP(ST+2)
            ELSEIF (I.EQ.1) THEN
               B2=C1; B3= PERMGROUP(ST+2)
            ELSE
               B2=PERMGROUP(ST+1); B3=C1
            ENDIF

            RF2(3*(C1-1)+1:3*C1)=TMP_AR(B1,:)
            RF2(3*(B1-1)+1:3*B1)=TMP_AR(C1,:)
!            PRINT*, "swap", B1, C1
            PREVDIH(:) = STARTDIH(:)
            CALL GETDIHONLY(RF2)
            INTDIFF= PREVDIH(:) - STARTDIH(:)
            INTDIST= SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))
            IF (PTEST) PRINT*, B1, C1,"intdist",intdist

            CALL DISTANCEPAIRSWAP(STARTDIH,RF2,B2,B3,PG2,INTDIST2,PTEST,RFNEW)
            RFNEW(:)=RF2(:)
            !PRINT*, "swap",B1, B3,B2,"intdist2", intdist2

            IF (INTDIST.LT.MININTDIST) THEN
               MININTDIST=INTDIST; MINSWAP(JJ,1)=B1;MINSWAP(jj,2)=B2;MINSWAP(jj,3)=B3
               !PRINT*, MINSWAP(JJ,:), MININTDIST
            ENDIF
            IF (INTDIST2.LT.MININTDIST) THEN
               MININTDIST=INTDIST2; MINSWAP(jj,1)=B1;MINSWAP(jj,2)=B3;MINSWAP(jj,3)=B2
               !PRINT*, MINSWAP(JJ,:), MININTDIST
            ENDIF
         ENDDO
         RF2(:)=RF(:)
         !PRINT*, "minswap",MINSWAP(JJ,:)
         RF2(3*(C1-1)+1:3*C1)=RF(3*(MINSWAP(JJ,1)-1)+1:3*MINSWAP(JJ,1))
         RF2(3*(C2-1)+1:3*C2)=RF(3*(MINSWAP(JJ,2)-1)+1:3*MINSWAP(JJ,2))
         RF2(3*(C3-1)+1:3*C3)=RF(3*(MINSWAP(JJ,3)-1)+1:3*MINSWAP(JJ,3))
         !PRINT*, "new curintdist", minintdist
         CURINTDIST = MININTDIST
      ENDDO
      PRINT*, "Curintdist after 1st bit", curintdist     

      !swap nitrogens and H's along with them
      CALL DISTANCEPAIRSWAP(STARTDIH,RF3,A1,A2,PERMG,INTDIST,PTEST,RFNEW)
      RF3(:)=RFNEW(:)
      CURINTDIST2=INTDIST
      IF (PTEST) PRINT*, "swap", a1,a2,intdist
      !swap H's
      
      MININTDIST=CURINTDIST2
      DO JJ =1,2
         PG2= PERM_AR(JJ)
         C1=PERMGROUP(START_AR(JJ))
         C2=PERMGROUP(START_AR(JJ)+1)
         C3=PERMGROUP(START_AR(JJ)+2)
         IF (PTEST) PRINT*,pg2, "pg2", c1,c2, c3
         MINSWAP(JJ,1)=C1; MINSWAP(JJ,2)=C2; MINSWAP(JJ,3)=C3
         TMP_AR(C1,:)= RF3(3*(C1-1)+1:3*C1)
         TMP_AR(C2,:)= RF3(3*(C2-1)+1:3*C2)
         TMP_AR(C3,:)= RF3(3*(C3-1)+1:3*C3)
         ST=START_AR(JJ)
!        PRINT*, "sec group3, start", st
         DO I=0,2
            RF3(3*(C1-1)+1:3*C1)=TMP_AR(C1,:)
            RF3(3*(C2-1)+1:3*C2)=TMP_AR(C2,:)
            RF3(3*(C3-1)+1:3*C3)=TMP_AR(C3,:)!initialise back to beginning
            !swap first two atoms
            B1=PERMGROUP(ST+I)
            !find other two atoms
            IF (I.EQ.0) THEN !i.e. B1 = A1
               B2=PERMGROUP(ST+1); B3= PERMGROUP(ST+2)
            ELSEIF (I.EQ.1) THEN
               B2=C1; B3= PERMGROUP(ST+2)
            ELSE
               B2=PERMGROUP(ST+1); B3=C1
            ENDIF

            RF3(3*(C1-1)+1:3*C1)=TMP_AR(B1,:)
            RF3(3*(B1-1)+1:3*B1)=TMP_AR(C1,:)
!            PRINT*, "swap", B1, C1
            PREVDIH(:) = STARTDIH(:)
            CALL GETDIHONLY(RF3)
            INTDIFF= PREVDIH(:) - STARTDIH(:)
            INTDIST= SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))
 
            CALL DISTANCEPAIRSWAP(STARTDIH,RF3,B2,B3,PERMG2,INTDIST2,PTEST,RFNEW)
            IF (PTEST) PRINT*, "minintdist", MININTDIST
            IF (PTEST) PRINT*, "swap",B1, B3,B2,"intdist2", intdist2

            IF (INTDIST.LT.MININTDIST) THEN
                MININTDIST=INTDIST; MINSWAP(JJ,1)=B1;MINSWAP(JJ,2)=B2;MINSWAP(JJ,3)=B3
                IF (PTEST) PRINT*, "minswap",  MINSWAP(JJ,:), MININTDIST
            ENDIF
            IF (INTDIST2.LT.MININTDIST) THEN
                MININTDIST=INTDIST2; MINSWAP(JJ,1)=B1;MINSWAP(JJ,2)=B3;MINSWAP(JJ,3)=B2
                IF (PTEST) PRINT*, "minswap", MINSWAP(JJ,:), MININTDIST
            ENDIF
         ENDDO
         IF (PTEST) PRINT*, "minswap",MINSWAP(JJ,:) 
         RF3(3*(C1-1)+1:3*C1)=RF(3*(MINSWAP(JJ,1)-1)+1:3*MINSWAP(JJ,1))
         RF3(3*(C2-1)+1:3*C2)=RF(3*(MINSWAP(JJ,2)-1)+1:3*MINSWAP(JJ,2))
         RF3(3*(C3-1)+1:3*C3)=RF(3*(MINSWAP(JJ,3)-1)+1:3*MINSWAP(JJ,3))
         CURINTDIST2= MININTDIST
      ENDDO


      IF (CURINTDIST2.LT.CURINTDIST) THEN
         RFNEW(:)=RF3(:) !which is still the same as the orig one if no changes were better
         CURINTDIST=CURINTDIST2
      ELSE
         RFNEW(:)=RF2(:)
      ENDIF
      IF (PTEST) PRINT*, "PERMG", PERMG, PERM_AR(1), PERM_AR(2)
      SKIPPG(PERMG)=.TRUE.;SKIPPG(PERM_AR(1))=.TRUE.;SKIPPG(PERM_AR(2))=.TRUE.

      DO JJ=1,NATOMS
         PRINT '(3f12.7)', RFNEW(3*(jj-1)+1:3*jj)
      ENDDO

      RETURN

      END SUBROUTINE LEU_SWAP


! *****************************************************************************

SUBROUTINE SWAP3ATONCE(STARTDIH,RF,PERMG,START, PERM_ARIN, START_ARIN,PTEST,SKIPPG,RFNEW,MININTDIST, SWAPARRAY, SWAPDIST)
      USE KEY, ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP,NSETS, SETS
!         USE MODAMBER9, ONLY : ih, ix, NRES, i02,m04,m02
      USE intcommons, only: NDIH, PERMNEIGHBOURS, PERMCHAIN,PREVDIH
      USE commons, only: NATOMS !msb50
      USE porfuncs
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: STARTDIH(NDIH)
      DOUBLE PRECISION, INTENT(IN) :: RF(3*NATOMS)
      INTEGER, INTENT(IN) ::PERMG, START !where atoms are in pg array
      INTEGER , INTENT(IN):: PERM_ARIN(3),START_ARIN(3)
      LOGICAL, INTENT(IN) :: PTEST
      LOGICAL, INTENT(INOUT) :: SKIPPG(NPERMGROUP)
       DOUBLE PRECISION, INTENT(OUT) :: RFNEW(3*NATOMS)
      DOUBLE PRECISION, INTENT(OUT) :: MININTDIST
      INTEGER, INTENT(OUT) :: SWAPARRAY(24,8) !this thing stores possible atom combinations
      DOUBLE PRECISION, INTENT(OUT) :: SWAPDIST(24) !and their distances

      INTEGER :: ISTAT
      INTEGER :: PG,PG1, PG2, PG3, START2, A1, A2,II, MINDISTLOC,JJ,SW,A3,A4
      INTEGER :: START_AR(3), PERM_AR(3)
      DOUBLE PRECISION :: CURINTDIST, DIST, RF2(3*NATOMS), INTDIST
      INTEGER :: SWAPARRAY1(12,6), MINSWAP1(6), MINSWAP2(6)
      DOUBLE PRECISION :: SWAPDIST1(12)

      
      IF (SKIPPG(PERMG)) THEN
         !PRINT*, "swap3atonce, cycle"
         RETURN
      ENDIF
      SWAPARRAY(:,:)=0; SWAPDIST1=100000.00

      IF (PERM_ARIN(1).NE.0) THEN
         IF (PTEST) PRINT*, "PERM_AR coming in"
         PERM_AR(:)=PERM_ARIN(:)
         START_AR(:)= START_ARIN(:)
         PG1=PERM_AR(1); PG2=PERM_AR(2); PG3=PERM_AR(3)
         IF (PG2.GT.PG3) THEN
            PRINT*, "swap3atonce ERROR: wrong incoming PERM_AR"
            STOP
         ENDIF
      ELSE
         IF (NPERMSIZE(PERMCHAIN(PERMG,2))==3) THEN
            PG1=PERMCHAIN(PERMG,4); PERM_AR(1)=PG1 !start backwards
            PG2=PERMCHAIN(PERMG,3); PERM_AR(2)=PG2
            PG3=PERMCHAIN(PERMG,2); PERM_AR(3)=PG3
         ELSE
            PG1=PERMCHAIN(PERMG,2); PERM_AR(1)=PG1
            PG2=PERMCHAIN(PERMG,3); PERM_AR(2)=PG2
            PG3=PERMCHAIN(PERMG,4); PERM_AR(3)=PG3
         ENDIF
     
     !   PRINT*, "pgs in swap3", PG1, PG2, PG3
         START_AR(:)=START
         DO JJ=1,3
            IF (PERM_AR(JJ)==PERMG) THEN
               START_AR(JJ)=START
            ELSE !it has to be larger as otherwise SKIPPG(PERMG) would be true
               DO II=PERMG,NPERMGROUP
                  IF (II.EQ.PERM_AR(JJ)) EXIT
                  START_AR(JJ) = START_AR(JJ)+ NPERMSIZE(II)
               ENDDO
            ENDIF
         ENDDO
         IF (PG2.GT.PG3) THEN
            START2=START_AR(3)
            PERM_AR(3)=PG2;START_AR(3)=START_AR(2)
            PERM_AR(2)=PG3; START_AR(2)=START2
            PG2=PERM_AR(2); PG3=PERM_AR(3)
         ENDIF
      ENDIF
      START2=START_AR(2)
      !PRINT*, "STARTAR", START_AR(:)
      IF (PTEST) PRINT*, "swap3atonce PERMAR", PERM_AR(:)
      !PRINT*, "STARTAR", START_AR(:)
      !either PG2 or PG3 have to be in the middle, i.e. involve 2h only
      IF (NSETS(PG2).GT.0.AND.NSETS(PG3).GT.0) THEN
         PRINT*, "swap3atonce died because of wrong nswap"
         STOP
      ENDIF
 
      !PRINT*, "npermsize(Pg1)", NPERMSIZE(PG1)
      IF (NPERMSIZE(PG1).NE.2) THEN
         PRINT*, "swap3atonce died because of wrong npermsize"
         STOP
      ENDIF
       
      A1=PERMGROUP(START_AR(1))
      A2=PERMGROUP(START_AR(1)+1)
      !PRINT*, "A1, a2", A1, A2

      SWAPARRAY(1:12,1)=A1; SWAPARRAY(1:12,2)=A2
      CALL SWAP2ATONCE(STARTDIH,RF,PG2,START2,PERM_AR(3),PTEST,SKIPPG,RFNEW,MINSWAP1,MINSWAP2,MININTDIST, SWAPARRAY1, SWAPDIST1)
      SWAPARRAY(1:12,3:8)=SWAPARRAY1(:,:)
      SWAPDIST(1:12)=SWAPDIST1(:)
   !   DO II=1,12
   !      PRINT '(8i5,f12.7)', SWAPARRAY(II,:), SWAPDIST(II)
   !   ENDDO

      CALL DISTANCEPAIRSWAP(STARTDIH, RF, A1, A2, PG1, INTDIST, PTEST,RF2)
      !PRINT*, "intdist", INTDIST
      SWAPARRAY(13:24,1)=A2; SWAPARRAY(13:,2)=A1
   !   PRINT*, "swaparray test before second 2 at once"

      SKIPPG(PG2)=.FALSE.; SKIPPG(PG3)=.FALSE.!otherwise it doesn't do it
      CALL  SWAP2ATONCE(STARTDIH,RF2,PG2,START2, PERM_AR(3),PTEST,SKIPPG,RFNEW,MINSWAP1,MINSWAP2,MININTDIST, SWAPARRAY1, SWAPDIST1)
      SWAPARRAY(13:,3:)=SWAPARRAY1(:,:)
      SWAPDIST(13:)=SWAPDIST1(:)
      IF (PTEST.AND.(PERM_ARIN(1).EQ.0)) THEN
      PRINT*, "the array at the end"
      DO II=1,24
         PRINT '(8i5,f12.7)', SWAPARRAY(II,:), SWAPDIST(II)
      ENDDO
      ENDIF

      CURINTDIST=1000000.00
      MINDISTLOC=1; MININTDIST=SWAPDIST(1)
      DO II=2,24
        IF (SWAPDIST(II).LT.MININTDIST) THEN
          MININTDIST=SWAPDIST(II)
          MINDISTLOC=II
        ENDIF
      ENDDO
      IF (PTEST) PRINT*, "minintdist", minintdist, "mindistloc", mindistloc
      !PRINT*, "minswap", SWAPARRAY(MINDISTLOC,:)

      RFNEW(:)=RF(:)

      !coordinate bookkeeping
      DO II=1,8
        A1=SWAPARRAY(1,II)
        IF (A1.EQ.0) THEN 
         !  PRINT*, II,"a1 is zero"
           CYCLE
        ENDIF
        A2=SWAPARRAY(MINDISTLOC,II)
        IF (PTEST) PRINT*, "rewriting coords", A1, A2
        RFNEW(3*(A1-1)+1:3*A1)=RF(3*(A2-1)+1:3*A2)

        !rest for dragging groups - always when you're done with this group
        IF (II.EQ.NPERMSIZE(PG1)) THEN !i.e. 2
           JJ=1; PG1=PERM_AR(JJ)
           IF (NSETS(PG1).EQ.0) CYCLE
        !know npermsize(pg2)=2, nswap(pg2)=0
        ELSEIF (II.EQ.2+NPERMSIZE(PG2)) THEN
           JJ=2; PG2=PERM_AR(JJ)
           IF (NSETS(PG2).EQ.0) CYCLE
        ELSEIF (II.EQ.2+NPERMSIZE(PG2)+NPERMSIZE(PG3)) THEN
           JJ=3; PG3=PERM_AR(JJ)
           IF (NSETS(PG3).EQ.0) CYCLE
        ELSE
           JJ=0
        ENDIF
        IF (JJ.NE.0) THEN
           IF (NSETS(PERM_AR(JJ)).NE.0) THEN
              IF (A1.NE.A2) THEN
                 DO SW = 1,NSETS(PERMG)
!                   A3 = SWAP1(PERM_AR(JJ), SW)
!                   A4 = SWAP2(PERM_AR(JJ),SW)
                    A3 = SETS(A1,SW) ! check arguments DJW
                    A4 = SETS(A2,SW) ! check arguments DJW
                    IF (PTEST) PRINT*, "dragging swap in swap3tonce", A3, A4
                    RFNEW(3*(A3-1)+1:3*A3) = RF(3*(A4-1)+1:3*A4)
                    RFNEW(3*(A4-1)+1:3*A4) = RF(3*(A3-1)+1:3*A3)
                 ENDDO
              ENDIF
           ENDIF
        ENDIF
      ENDDO
     
!      DO II=1, 3*NATOMS
!        IF (RFNEW(II).NE.RF(II)) PRINT '(a13,i5,a7,2f12.7)', "have swapped", II, "coords", RFNEW(II), RF(II)
!      ENDDO
 
      SKIPPG(PG1)=.TRUE.;SKIPPG(PG2)=.TRUE.; SKIPPG(PG3)=.TRUE.

      RETURN
      END SUBROUTINE SWAP3ATONCE

! *********************************************************************************
    SUBROUTINE OLD_INTMINPERM(RS,RF,DISTANCE,RMAT,PTEST)
    ! minimise Cartesian distance without permuting then
    ! permute the permutable groups to minimize distance in single torsions
    ! changes RF to resulting best alignment with RS
    ! output distance is cartesian distance **squared**

    USE KEY, ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, SETS, NABT, AMBERT
    USE modamber9, only:nres, ih, m02,  ix, i02, m04

    IMPLICIT NONE
    
    DOUBLE PRECISION, INTENT(IN) :: RS(3*NATOMS)
    DOUBLE PRECISION, INTENT(INOUT) :: RF(3*NATOMS), RMAT(3,3)
    DOUBLE PRECISION, INTENT(OUT) :: DISTANCE
    LOGICAL, INTENT(IN) :: PTEST

    INTEGER :: PG, SW, I, J, K, A1, A2, A3, B1, B2, B3, START,III
    DOUBLE PRECISION :: RSINT(NINTC), RFINT(NINTC), INTDIFF(NDIH), INTDIFF2(NDIH)
    DOUBLE PRECISION :: INTDIST, CURINTDIST, DISTANCE2
    DOUBLE PRECISION :: TMP(3), TMP2(3), TMP3(3), STARTDIH(NDIH), RFTMP(3*NATOMS)
    INTEGER :: prol(NRES), glyc(NRES), nprol, nglyc,II
    LOGICAL :: prolswap, glyswap

    ! alignment stuff
    DOUBLE PRECISION :: DISTF
    CHARACTER(LEN=5) :: ZSYMSAVE
    COMMON /SYS/ ZSYMSAVE
    

      prol(:) = 0;nprol=0
      glyc(:) = 0;nglyc=0
      glyswap = .FALSE.;prolswap=.FALSE.
      IF (AMBERT.OR.NABT) THEN !msb50 - for amber9.ff03 and glycine this doesn't work
          DO III = 1,NRES!as swapping the hydrogens makes no difference in natural internals in this case
            IF (ih(m02+III-1).EQ.'GLY'.OR.ih(m02+III-1).EQ.'NGLY'.OR.ih(m02+III-1).EQ.'CGLY') THEN
               nglyc = nglyc+1
               glyc(nglyc) = III
            ENDIF
            IF (ih(m02+III-1).EQ.'PRO'.OR.ih(m02+III-1).EQ.'NPRO'.OR.ih(m02+III-1).EQ.'CPRO') THEN
               nprol = nprol+1
               prol(nprol) =III
            ENDIF
          ENDDO
      ENDIF

    ! align without permuting
    !msb50 - do at the end in minpermdist now
   !CALL NEWMINDIST(RS(:),RF(:),NATOMS,DISTANCE,.FALSE.,.FALSE.,ZSYM(1),.FALSE.,.FALSE.,.FALSE.,RMAT)


    PREVDIH => DIHINFOSINGLE(:)

    PREVDIH(:) = 0.0D0
    CALL CART2INT(RS, RSINT)

    STARTDIH(:) = PREVDIH(:)
    ALIGNDIR = .FALSE.

    CALL CART2INT(RF, RFINT)

    INTDIFF(:) = PREVDIH(:)-STARTDIH(:)
    CURINTDIST = SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))
    IF (PTEST) print*, 'Starting intdistance: ', CURINTDIST

    START = 1
    DO PG = 1,NPERMGROUP

       IF (NPERMSIZE(PG).EQ.2) THEN
          RFTMP(:) = RF(:)

          ! A1 and A2 are the atom numbers to permute
          A1 = PERMGROUP(START)
          A2 = PERMGROUP(START+1)

          IF (PTEST) PRINT '(A,2I6)','old_intminperm> Swapping atoms: ', A1, A2

           !could put in ngly bit too, but unnecessary for CHARMM 19
           !and not set up for CHARMM 22
            IF (nprol.NE.0) THEN
              prolswap=.FALSE.
              DO III=1, nprol
                 IF (ix(i02+prol(III)-1).LT.A1 .AND. ix(i02+prol(III)).GT.A1) THEN
                    IF (PTEST) PRINT*, "pro", prol(III)
                    CALL PROL_PERMUTE(A1,A2,PTEST, RS, RF, DISTANCE)
                    prolswap =.TRUE.
                    EXIT
                 ENDIF
              ENDDO
              IF (prolswap) THEN
                 START= START+NPERMSIZE(PG)
                 CYCLE
              ENDIF
           ENDIF



          RF(3*(A1-1)+1:3*A1) = RFTMP(3*(A2-1)+1:3*A2)
          RF(3*(A2-1)+1:3*A2) = RFTMP(3*(A1-1)+1:3*A1)

          !move any other groups that must be dragged along
          DO SW = 1,NSETS(PG)
!            A1 = SWAP1(PG,SW)
!            A2 = SWAP2(PG,SW)
             A1 = SETS(PERMGROUP(START),SW)
             A2 = SETS(PERMGROUP(START+1),SW)

             IF (PTEST) print*, 'Dragging swap: ', A1, A2

             RF(3*(A1-1)+1:3*A1) = RFTMP(3*(A2-1)+1:3*A2)
             RF(3*(A2-1)+1:3*A2) = RFTMP(3*(A1-1)+1:3*A1)
          ENDDO

          PREVDIH(:) = STARTDIH(:)
          CALL CART2INT(RF,RFINT)

          INTDIFF(:) = PREVDIH(:) - STARTDIH(:)
          INTDIST = SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))

          IF (PTEST) PRINT '(A,G20.10)','old_intminperm> torsion distance: ', INTDIST

          IF (INTDIST.GT.CURINTDIST) THEN
             ! undo permutation
             RF(:) = RFTMP(:)
          ELSE
             IF (PTEST) print*, 'keep permutation'
             CURINTDIST = INTDIST
          ENDIF
       ELSE IF (NPERMSIZE(PG).EQ.3) THEN
          ! this is an inefficient way of listing permutations
          ! and should be fixed and generalized at some point
          A1 = PERMGROUP(START)
          A2 = PERMGROUP(START+1)
          A3 = PERMGROUP(START+2)

          TMP(:) = RF(3*(A1-1)+1:3*A1)
          TMP2(:) = RF(3*(A2-1)+1:3*A2)
          TMP3(:) = RF(3*(A3-1)+1:3*A3)

          RFTMP(:) = RF(:)
          DO I = 0,2
             DO J = 0,2
                IF (J.EQ.I.OR.(I.EQ.0.AND.J.EQ.1)) CYCLE
                DO K = 0,2
                   IF (K.EQ.I.OR.K.EQ.J) CYCLE

                   B1 = PERMGROUP(START+I)
                   B2 = PERMGROUP(START+J)
                   B3 = PERMGROUP(START+K)

                   IF (PTEST) PRINT '(A,3I5)','old_intminperm> permuting: ', B1, B2, B3

                   RF(3*(B1-1)+1:3*B1) = TMP(:)
                   RF(3*(B2-1)+1:3*B2) = TMP2(:)
                   RF(3*(B3-1)+1:3*B3) = TMP3(:)
 

                   PREVDIH(:) = STARTDIH(:)
                   CALL CART2INT(RF,RFINT)

                   INTDIFF(:) = PREVDIH(:) - STARTDIH(:)
                   INTDIST = SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))

                   IF (PTEST) PRINT '(A,G20.10)','old_intminperm> dihedral distance: ', INTDIST

                   IF (INTDIST.LT.CURINTDIST) THEN
                      RFTMP(:) = RF(:)
                      CURINTDIST = INTDIST
                      IF (PTEST) print*, 'keep permutation'
                   ENDIF
                ENDDO
             ENDDO
          ENDDO
          RF(:) = RFTMP(:)
       ELSE
          print*, 'Error! INTMINPERM is only set up for permutation groups &
               &          of size 2 & 3 for now', PG, NPERMSIZE(PG)
          STOP
       ENDIF

       START = START + NPERMSIZE(PG)
    END DO

    IF (PTEST) print '(A,G20.10)','old_intminperm> Final int distance: ', CURINTDIST

!msb50 consider calling my orient and check then

    DISTANCE = DOT_PRODUCT(RF-RS,RF-RS)

!    CALL CART2INT(RF,RFINT)
!    PRINT*, "finish rfint"
!    PRINT*, RFINT(:)
    
    IF (PTEST) PRINT '(A,G20.10)',"old_intminperm> dist in old_intminperm", distance
    RETURN
  END SUBROUTINE OLD_INTMINPERM
!**********************************************************************************
SUBROUTINE SWAP4ATONCE(STARTDIH,RF,PERMG,START,PERM_AR, START_AR,PTEST,SKIPPG, SWAPARRAY, SWAPDIST)
      USE KEY, ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP,NSETS, SETS
!         USE MODAMBER9, ONLY : ih, ix, NRES, i02,m04,m02
      USE intcommons, only: NDIH, PERMNEIGHBOURS, PERMCHAIN,PREVDIH
      USE commons, only: NATOMS !msb50
      USE porfuncs
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: STARTDIH(NDIH)
      DOUBLE PRECISION, INTENT(IN) :: RF(3*NATOMS)
      INTEGER, INTENT(IN) ::PERMG, START !where atoms are in pg array
      LOGICAL, INTENT(IN) :: PTEST
      LOGICAL, INTENT(INOUT) :: SKIPPG(NPERMGROUP)
      INTEGER, INTENT(OUT) :: SWAPARRAY(48,10) !this thing stores possible atom combinations
      DOUBLE PRECISION, INTENT(OUT) :: SWAPDIST(48) !and their distances
      INTEGER, INTENT(IN) :: PERM_AR(4),START_AR(4)

      INTEGER :: ISTAT
      INTEGER :: PG4,PG1, PG2, PG3, START2, A1, A2,II, MINDISTLOC,JJ,SW,A3,A4
      DOUBLE PRECISION :: CURINTDIST, DIST, RF2(3*NATOMS), INTDIST
      DOUBLE PRECISION :: RFNEW(3*NATOMS)
      INTEGER :: SWAPARRAY1(24,8)
      DOUBLE PRECISION :: SWAPDIST1(24), MININTDIST

      IF (SKIPPG(PERMG)) THEN
         IF (PTEST) PRINT*, "swap4atonce, cycle"
         RETURN
      ENDIF
      SWAPARRAY(:,:)=0; SWAPDIST1=100000.00

      PG1=PERM_AR(1);PG2=PERM_AR(2);PG3=PERM_AR(3);PG4=PERM_AR(4)
      IF (NPERMSIZE(PG1)==3) THEN
         PRINT*, "wrong permsize of pg in swap4atonce"
         STOP
      ENDIF

      START2=START_AR(2)
      !PRINT*, "4,STARTAR", START_AR(:)
      !PRINT*, "4,PERMAR", PERM_AR(:)

      IF (NPERMSIZE(PG1).NE.2) THEN
         PRINT*, "swap4atonce died because of wrong npermsize"
         STOP
      ENDIF

      A1=PERMGROUP(START_AR(1))
      A2=PERMGROUP(START_AR(1)+1)
      IF (PTEST) PRINT*, "A1, a2", A1, A2

      SWAPARRAY(1:24,1)=A1; SWAPARRAY(1:24,2)=A2
      CALL SWAP3ATONCE(STARTDIH,RF,PG2,START2,PERM_AR(2:), START_AR(2:),PTEST,SKIPPG,RFNEW,MININTDIST, SWAPARRAY1, SWAPDIST1)
      SWAPARRAY(1:24,3:10)=SWAPARRAY1(:,:)
      SWAPDIST(1:24)=SWAPDIST1(:)

      CALL DISTANCEPAIRSWAP(STARTDIH, RF, A1, A2, PG1, INTDIST, PTEST,RF2)
      SWAPARRAY(25:48,1)=A2; SWAPARRAY(25:48,2)=A1

      SKIPPG(PG2)=.FALSE.; SKIPPG(PG3)=.FALSE.; SKIPPG(PG4)=.FALSE.
      CALL  SWAP3ATONCE(STARTDIH,RF2,PG2,START2,PERM_AR(2:), START_AR(2:),PTEST,SKIPPG,RFNEW,MININTDIST, SWAPARRAY1, SWAPDIST1)
      SWAPARRAY(25:,3:)=SWAPARRAY1(:,:)
      SWAPDIST(25:)=SWAPDIST1(:)
      !PRINT*, "swap4atonce,the array at the end"
      !DO II=1,48
      !   PRINT '(10i5,f12.7)', SWAPARRAY(II,:), SWAPDIST(II)
      !ENDDO

      RETURN
      END SUBROUTINE SWAP4ATONCE

! ********************************************************************************
      SUBROUTINE SWAP5ATONCE(STARTDIH,RF,PERMG,START,PTEST,SKIPPG,RFNEW,MININTDIST, SWAPARRAY, SWAPDIST)
! this is basically just for lysine - swap4atonce is never actually used as it doesn't exist
      USE KEY, ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP,NSETS, SETS
!         USE MODAMBER9, ONLY : ih, ix, NRES, i02,m04,m02
      USE intcommons, only: NDIH, PERMNEIGHBOURS, PERMCHAIN,PREVDIH
      USE commons, only: NATOMS !msb50
      USE porfuncs
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(IN) :: STARTDIH(NDIH)
      DOUBLE PRECISION, INTENT(IN) :: RF(3*NATOMS)
      INTEGER, INTENT(IN) ::PERMG, START !where atoms are in pg array
      LOGICAL, INTENT(IN) :: PTEST
      LOGICAL, INTENT(INOUT) :: SKIPPG(NPERMGROUP)
      DOUBLE PRECISION, INTENT(OUT) :: RFNEW(3*NATOMS)
      DOUBLE PRECISION, INTENT(OUT) :: MININTDIST
      INTEGER, INTENT(OUT) :: SWAPARRAY(96,12) !this thing stores possible atom combinations
      DOUBLE PRECISION, INTENT(OUT) :: SWAPDIST(96) !and their distances

      INTEGER :: ISTAT
      INTEGER :: PG4,PG1, PG2,PG3,PG5,START2, A1, A2,II, MINDISTLOC,JJ,SW,A3,A4,KK
      INTEGER :: START_AR(5), PERM_AR(5), SORT_PG(5), SORT_IND(5)
      INTEGER :: PERMCHSTARTS(2:6), START_SORT_AR(5)
      INTEGER :: THREEPG
      LOGICAL :: already_in
      DOUBLE PRECISION :: CURINTDIST, DIST, RF2(3*NATOMS), INTDIST
      INTEGER :: SWAPARRAY1(48,10)
      DOUBLE PRECISION :: SWAPDIST1(48)

      IF (SKIPPG(PERMG)) THEN
         IF (PTEST) PRINT*, "swap4atonce, cycle"
         RETURN
      ENDIF
      SWAPARRAY(:,:)=0; SWAPDIST1=100000.00
 
      !PRINT*, "PERMG", PERMG
      SORT_IND(:)=0; SORT_PG(:)=NPERMGROUP+2
      DO JJ=2,6
        IF (PERMCHAIN(PERMG,JJ).LT. SORT_PG(1)) THEN
           SORT_PG(1)=PERMCHAIN(PERMG,JJ)
           SORT_IND(JJ-1)=1
!           SORT_PG(2)=PERMCHAIN(PERMG, MOD(JJ-1,5)+2)!next one in the row
!           PRINT*, JJ, MOD(JJ-1,5)+2, "sug pg", SORT_PG(2)
        ENDIF
      ENDDO
      DO KK=2,5
        DO JJ=2,6
           already_in=.FALSE.
           DO II=1,KK-1
              IF (PERMCHAIN(PERMG,JJ).EQ.SORT_PG(II)) THEN
                 already_in=.TRUE. !i.e. the ones that are smaller already
              ENDIF
           ENDDO
           IF (already_in) CYCLE
           IF (PERMCHAIN(PERMG,JJ).LT.SORT_PG(KK)) THEN
              SORT_PG(KK)=PERMCHAIN(PERMG,JJ)
              SORT_IND(JJ-1)=KK!where in sort_pg do i find PERMCHAIN(PERMG,JJ)
           ENDIF
         ENDDO
      ENDDO
      !PRINT*, "sort pg", SORT_PG(:), SORT_IND(:)
      !PRINT*, "sort ind", SORT_IND(:)

      START_SORT_AR(:)=START
      DO JJ=1,5
        IF (SORT_PG(JJ).EQ.PERMG) THEN
           START_SORT_AR(JJ)=START
        ELSEIF (SORT_PG(JJ).GT.PERMG) THEN
           DO II=PERMG,NPERMGROUP
              IF (II.EQ.SORT_PG(JJ)) EXIT
              START_SORT_AR(JJ)=START_SORT_AR(JJ)+NPERMSIZE(II)
           ENDDO
        ELSEIF (SORT_PG(JJ).LT.PERMG) THEN
           DO II=PERMG, 2, -1
              IF (II.EQ.SORT_PG(JJ)) EXIT
              START_SORT_AR(JJ)=START_SORT_AR(JJ)-NPERMSIZE(II-1)
           ENDDO
        ENDIF
      ENDDO
      IF (PTEST) PRINT*, "swap5@once: SORT_START", START_SORT_AR(:)
      !in order of permgroups in permchain - starts
      DO JJ=2,6
        PERMCHSTARTS(JJ)=START_SORT_AR(SORT_IND(JJ-1))
      ENDDO
      IF (PTEST) PRINT*, "swap5@once: PERMCHSTARTS", PERMCHSTARTS(2:6)

      THREEPG=0
      DO II=2,6
        IF (NPERMSIZE(PERMCHAIN(PERMG,II))==3) THEN
          THREEPG=PERMCHAIN(PERMG,II)
        ENDIF
      ENDDO
      IF (THREEPG==0) THEN
         PRINT*, "error in swap5atonce, no three permgroups, no lysine"
         !setup should work for this too, just need to figure how to order the permgroups
      ENDIF

      !swap2atonce has to be called with the smallest of two neighbouring pg's
      !pgs have to be in order according to their lineup, only
      !threepg has to wait for swapatonce (i.e. be on place 4 or 5)
      IF (THREEPG==SORT_PG(5)) THEN
         PG5=THREEPG; PERM_AR(5)=THREEPG
         IF (THREEPG.EQ.PERMCHAIN(PERMG,6)) THEN
            PG1=PERMCHAIN(PERMG,2); PERM_AR(1)=PG1
            PG2=PERMCHAIN(PERMG,3); PERM_AR(2)=PG2
            PG3=PERMCHAIN(PERMG,4); PERM_AR(3)=PG3
            PG4=PERMCHAIN(PERMG,5); PERM_AR(4)=PG4
            DO JJ=2,6
              START_AR(JJ-1)=PERMCHSTARTS(JJ)
            ENDDO
         ELSEIF (THREEPG.EQ.PERMCHAIN(PERMG,2)) THEN
            PG1=PERMCHAIN(PERMG,6); PERM_AR(1)=PG1
            PG2=PERMCHAIN(PERMG,5); PERM_AR(2)=PG2
            PG3=PERMCHAIN(PERMG,4); PERM_AR(3)=PG3
            PG4=PERMCHAIN(PERMG,3); PERM_AR(4)=PG4
            DO JJ=2,6
              START_AR(JJ-1)=PERMCHSTARTS(8-JJ)
            ENDDO
         ELSE
            PRINT*, "ERROR in swap5atonce, npermsize 3 and not at end of permchain!"
           STOP
         ENDIF
      ELSE
         IF (THREEPG.EQ.PERMCHAIN(PERMG,6)) THEN
         !swap2atonce has to be called with the smallest of two neighbouring pg's
            IF (PERMCHAIN(PERMG,5).LT.THREEPG) THEN
               PG5=THREEPG; PERM_AR(5)=THREEPG
               PG4=PERMCHAIN(PERMG,5); PERM_AR(4)=PG4
               START_AR(5)=PERMCHSTARTS(6)
               START_AR(4)=PERMCHSTARTS(5)
            ELSE
               PG4=THREEPG; PERM_AR(4)=THREEPG
               PG5=PERMCHAIN(PERMG,5); PERM_AR(5)=PG5
               START_AR(4)=PERMCHSTARTS(6)
               START_AR(5)=PERMCHSTARTS(5)
            ENDIF
            PG1=PERMCHAIN(PERMG,2); PERM_AR(1)=PG1
            !sort according to neighbours, as START_AR will be given to the subsequent swapXatonce
            PG2=PERMCHAIN(PERMG,3);PERM_AR(2)=PG2
            PG3=PERMCHAIN(PERMG,4); PERM_AR(3)=PG3
            DO JJ=2,4
               START_AR(JJ-1)=PERMCHSTARTS(JJ)
            ENDDO
         ELSEIF (THREEPG.EQ.PERMCHAIN(PERMG,2)) THEN
            IF (PERMCHAIN(PERMG,3).LT.THREEPG) THEN
               PG5=THREEPG;PERM_AR(5)=THREEPG
               PG4=PERMCHAIN(PERMG,3); PERM_AR(4)=PG4
               START_AR(5)=PERMCHSTARTS(2)
               START_AR(4)=PERMCHSTARTS(3)
            ELSE
               PG4=THREEPG; PERM_AR(4)=THREEPG
               PG5=PERMCHAIN(PERMG,3); PERM_AR(5)=PG5
               START_AR(4)=PERMCHSTARTS(2)
               START_AR(5)=PERMCHSTARTS(3)
            ENDIF
            PG1=PERMCHAIN(PERMG,6); PERM_AR(1)=PG1
            PG2=PERMCHAIN(PERMG,5); PERM_AR(2)=PG2
            PG3=PERMCHAIN(PERMG,4); PERM_AR(3)=PG3
            DO JJ=2,4
               START_AR(JJ-1)=PERMCHSTARTS(8-JJ)
            ENDDO
         ELSE
           PRINT*, "ERROR in swap5atonce, npermsize 3 and not at end of permchain!"
           STOP
         ENDIF
      ENDIF

      START2=START_AR(2)
      IF (PTEST) PRINT*,"PERMAR", PERM_AR(:)
      IF (PTEST) PRINT*, "STARTAR", START_AR(:)

      A1=PERMGROUP(START_AR(1))
      A2=PERMGROUP(START_AR(1)+1)
      IF (PTEST) PRINT*, "swap A1, a2", A1, A2

      SWAPARRAY(1:48,1)=A1; SWAPARRAY(1:48,2)=A2
      !always call with smallest possible pg apart from pg1  as this is expected and
      !the subroutine will reorder itself
      CALL SWAP4ATONCE(STARTDIH,RF,PG2,START2,PERM_AR(2:), START_AR(2:),PTEST,SKIPPG, SWAPARRAY1, SWAPDIST1)
      SWAPARRAY(1:48,3:12)=SWAPARRAY1(:,:)
      SWAPDIST(1:48)=SWAPDIST1(:)

      !now swap the two PG1 H's, and go on
      CALL DISTANCEPAIRSWAP(STARTDIH, RF, A1, A2, PG1, INTDIST, PTEST,RF2)
      !PRINT*, "intdist", INTDIST
      SWAPARRAY(49:96,1)=A2; SWAPARRAY(49:96,2)=A1

      SKIPPG(PG2)=.FALSE.; SKIPPG(PG3)=.FALSE.; SKIPPG(PG3)=.FALSE.;SKIPPG(PG4)=.FALSE.
      CALL SWAP4ATONCE(STARTDIH,RF2,PG2,START2,PERM_AR(2:), START_AR(2:),PTEST,SKIPPG, SWAPARRAY1, SWAPDIST1)
      SWAPARRAY(49:96,3:12)=SWAPARRAY1(:,:)
      SWAPDIST(49:96)=SWAPDIST1(:)
      IF (PTEST) THEN
        PRINT*, "swap5atonce, final array"
        DO II=1,96
           PRINT '(12i5,f12.7)', SWAPARRAY(II,:), SWAPDIST(II)
        ENDDO
      ENDIF

      CURINTDIST=1000000.00
      MINDISTLOC=1; MININTDIST=SWAPDIST(1)
      DO II=2,96
        IF (SWAPDIST(II).LE.MININTDIST) THEN
          MINDISTLOC=II; MININTDIST=SWAPDIST(II)
        ENDIF
      ENDDO
      IF (PTEST) PRINT*, "minintdist", minintdist
      IF (PTEST) PRINT*, "minswap", SWAPARRAY(MINDISTLOC,:)

      SKIPPG(PG2)=.TRUE.;SKIPPG(PG3)=.TRUE.;SKIPPG(PG3)=.TRUE.;SKIPPG(PG4)=.TRUE.
      !now the coordinates
      RFNEW(:)=RF(:)
      DO II=1,12
        A1=SWAPARRAY(1,II)
        IF (A1.EQ.0) THEN
         !  PRINT*, II,"a1 is zero"
           CYCLE
        ENDIF
        A2=SWAPARRAY(MINDISTLOC,II)
        !PRINT*, "rewriting coords", A1, A2
        RFNEW(3*(A1-1)+1:3*A1)=RF(3*(A2-1)+1:3*A2)

        !rest for dragging groups - always when you're done with this group
        IF (II.EQ.NPERMSIZE(PG1)) THEN !i.e. 2
           JJ=1; PG1=PERM_AR(JJ)
           IF (NSETS(PG1).EQ.0) CYCLE
        !know npermsize(pg2)=2, nswap(pg2)=0
        ELSEIF (II.EQ.2+NPERMSIZE(PG2)) THEN
           JJ=2; PG2=PERM_AR(JJ)
           IF (NSETS(PG2).EQ.0) CYCLE
        ELSEIF (II.EQ.2+NPERMSIZE(PG2)+NPERMSIZE(PG3)) THEN
           JJ=3; PG3=PERM_AR(JJ)
           IF (NSETS(PG3).EQ.0) CYCLE
        ELSEIF (II.EQ.2+NPERMSIZE(PG2)+NPERMSIZE(PG3)+NPERMSIZE(PG4)) THEN
           JJ=4
           IF (NSETS(PG4).EQ.0) CYCLE
        ELSEIF (II.EQ.2+NPERMSIZE(PG2)+NPERMSIZE(PG3)+NPERMSIZE(PG4)+NPERMSIZE(PG5)) THEN
           JJ=5
           IF (NSETS(PG5).EQ.0) CYCLE
        ELSE
           JJ=0
        ENDIF
        IF (JJ.NE.0) THEN
           IF (NSETS(PERM_AR(JJ)).NE.0) THEN
              IF (A1.NE.A2) THEN
                 DO SW = 1,NSETS(PERMG)
!                   A3 = SWAP1(PERM_AR(JJ), SW)
!                   A4 = SWAP2(PERM_AR(JJ),SW)
                    A3 = SETS(A1,SW) ! ??? DJW
                    A4 = SETS(A2,SW) ! ??? DJW
                    IF (PTEST) PRINT*, "dragging swap in swap3tonce", A3, A4
                    RFNEW(3*(A3-1)+1:3*A3) = RF(3*(A4-1)+1:3*A4)
                    RFNEW(3*(A4-1)+1:3*A4) = RF(3*(A3-1)+1:3*A3)
                 ENDDO
              ENDIF
           ENDIF
        ENDIF
      ENDDO

      END SUBROUTINE SWAP5ATONCE


      
      SUBROUTINE INTMINPERM_CHIRAL(RS,RF,DISTANCE,RMAT,PTEST)
      !see intcoords.f90
      use key, only:  NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, SETS, AMBERT, NABT
      use modamber9, only: PROCHIRALH,ih, ix, NRES, i02,m04,m02
      use commons, only: natoms
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)    :: RS(3*NATOMS)
      DOUBLE PRECISION, INTENT(INOUT) :: RF(3*NATOMS), RMAT(3,3)
      DOUBLE PRECISION, INTENT(OUT)   :: DISTANCE
      LOGICAL, INTENT(IN)             :: PTEST 

      DOUBLE PRECISION :: RFTMP(3*NATOMS), RFTMP2(3*NATOMS), RFNEW(3*NATOMS)
      INTEGER :: PG, SW, I, J, K, A1, A2, A3, B1, B2, B3, START, III
      DOUBLE PRECISION :: RSINT(NINTC), RFINT(NINTC), INTDIFF(NDIH)
      DOUBLE PRECISION :: INTDIST, CURINTDIST, INTDIST2
      DOUBLE PRECISION :: STARTDIH(NDIH), GD_PREVDIH(NDIH)
      DOUBLE PRECISION :: TMP_AR(NATOMS,3)
      DOUBLE PRECISION :: MININTDIST, DISTANCE2
      CHARACTER(LEN=4) :: resname
      INTEGER          :: nglyc, nprol, glyc(NRES), prol(NRES)
      LOGICAL          :: glyswap, prolswap
      INTEGER          :: MINSWAP(3)
      DOUBLE PRECISION :: ANGLES(NDIH), INTDIFFMSB50(NDIH), INTDISTMSB50

      DISTANCE = DOT_PRODUCT(RF-RS,RF-RS)
      IF (PTEST) PRINT*, "intminperm> initial cart dist", DISTANCE
 
!glycine and proline exceptions
      prol(:) = 0;nprol=0
      glyc(:) = 0;nglyc=0
      glyswap = .FALSE.;prolswap=.FALSE.
      IF (AMBERT.OR.NABT) THEN !msb50 - for amber9.ff03 and glycine this doesn't work
          DO III = 1,NRES!as swapping the hydrogens makes no difference in natural internals in this case
            IF (ih(m02+III-1).EQ.'GLY'.OR.ih(m02+III-1).EQ.'NGLY'.OR.ih(m02+III-1).EQ.'CGLY') THEN
               nglyc = nglyc+1
               glyc(nglyc) = III
            ENDIF
            IF (ih(m02+III-1).EQ.'PRO'.OR.ih(m02+III-1).EQ.'NPRO'.OR.ih(m02+III-1).EQ.'CPRO') THEN
               nprol = nprol+1
               prol(nprol) =III
            ENDIF
          ENDDO
      ENDIF
!

      PREVDIH => DIHINFOSINGLE(:)
      PREVDIH(:) = 0.0D0
      CALL CART2INT(RS, RSINT)
      
      STARTDIH(:) = PREVDIH(:)
      ALIGNDIR = .FALSE.

      CALL GETDIHONLY(RF)  
      INTDIFF(:) = PREVDIH(:)-STARTDIH(:)
      CURINTDIST = SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))
      IF (PTEST) print*, 'Starting intdistance: ', CURINTDIST

      RFTMP(:) = RF(:)
      IF (nglyc.NE.0.AND..NOT.GLYCART) THEN
         CALL GLYINTPERM(RS, RF, PTEST)
      ENDIF

      PREVDIH(:) = STARTDIH(:)
      CALL GETDIHONLY(RF)
      INTDIFF(:) = PREVDIH(:)-STARTDIH(:)
      CURINTDIST = SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))
      IF (PTEST) print*, 'after intgly intdistance: ', CURINTDIST   

      START = 1
      DO PG = 1,NPERMGROUP

      IF (NPERMSIZE(PG).EQ.2) THEN
        RFTMP(:) = RF(:)
       
        ! A1 and A2 are the atom numbers to permute
        A1 = PERMGROUP(START)
        A2 = PERMGROUP(START+1)
       
        IF (PTEST) PRINT '(A,2I5)',"intminperm_chiral> which PERMGROUP", PG, START
        IF (PTEST) PRINT '(A,2I5)','intminperm_chiral> Swapping atoms: ', A1, A2
       
        IF (nglyc.NE.0) THEN
          DO III=1, nglyc
          IF (ix(i02+glyc(III)-1).GT.A1) EXIT
          glyswap = .FALSE.
             IF (ix(i02+glyc(III)-1).LE.A1 .AND.(ix(i02+glyc(III)).GT.A1).AND.&
        &          (ih(m04+A1-1).EQ.'HA2 '.OR.ih(m04+A1-1).EQ.'HA3 ')) THEN
                IF (PTEST) PRINT*, "glyc", glyc(III)
                IF (.NOT.GLYCART) THEN
                   glyswap=.TRUE.
                   CYCLE
                ENDIF
               !the second atom has to be A2 in gly
                IF (ih(m04+A2-1).EQ.'HA3 '.OR. ih(m04+A2-1).EQ.'HA2 ') THEN
                   glyswap = .TRUE.
                   IF (PTEST) PRINT*, "glycine"
                   RF(3*(A1-1)+1:3*A1) = RFTMP(3*(A2-1)+1:3*A2)
                   RF(3*(A2-1)+1:3*A2) = RFTMP(3*(A1-1)+1:3*A1)
               !move any other groups that must be dragged along
                   IF (NSETS(PG).NE.0) THEN
                     PRINT*, "intminperm>this is not glycine. you DON'T know& 
         &                        what you are doing"
                     STOP
                   ENDIF
                   DISTANCE = DOT_PRODUCT(RFTMP-RS,RFTMP-RS)
                   DISTANCE2 = DOT_PRODUCT(RF-RS,RF-RS)
                   !IF (PTEST) PRINT*, "dist old", DISTANCE, "dist new", DISTANCE2
                   IF (DISTANCE2.LT.DISTANCE) THEN
                     IF (PTEST) print '(a40, f17.5,a10,f17.5)',"keep permutation&
         & aaccrding to cartesians (gly): dist1",DISTANCE, "dist_new", DISTANCE2
                   ELSE
                      RF(:) = RFTMP(:) !unddo
                   ENDIF
                ELSE
                   PRINT*, "intminperm> sth wrong with GLY"
                   STOP
                ENDIF
            ENDIF
            IF (glyswap) EXIT
           ENDDO
             IF (glyswap) THEN
                IF (PTEST) PRINT*, "glyswap true", PG
                START = START + NPERMSIZE(PG)
               CYCLE
             ENDIF
        ENDIF
       
        IF (nprol.NE.0) THEN
           prolswap=.FALSE.
           DO III=1, nprol
       
              IF (ix(i02+prol(III)-1).LT.A1 .AND. ix(i02+prol(III)).GT.A1) THEN
                 IF (PTEST) PRINT*, "pro", prol(III)
                 CALL PROL_PERMUTE(A1,A2,PTEST, RS, RF, DISTANCE)
                 prolswap =.TRUE.
                 EXIT
              ENDIF
           ENDDO
           IF (prolswap) THEN
              START= START+NPERMSIZE(PG)
              CYCLE
           ENDIF
        ENDIF
       
! now start the real thing:
        IF (.NOT.PROCHIRALH(A1)) THEN
           RF(3*(A1-1)+1:3*A1) = RFTMP(3*(A2-1)+1:3*A2)
           RF(3*(A2-1)+1:3*A2) = RFTMP(3*(A1-1)+1:3*A1)
         
           !move any other groups that must be dragged along
           DO SW = 1,NSETS(PG)
!             A1 = SWAP1(PG,SW)
!             A2 = SWAP2(PG,SW)
              A1 = SETS(PERMGROUP(START),SW)
              A2 = SETS(PERMGROUP(START+1),SW)
         
              IF (PTEST) print*, 'Dragging swap: ', A1, A2
         
              RF(3*(A1-1)+1:3*A1) = RFTMP(3*(A2-1)+1:3*A2)
              RF(3*(A2-1)+1:3*A2) = RFTMP(3*(A1-1)+1:3*A1)
           ENDDO
        
        ELSE
           CALL ALLCHIRALH_ALIGN(RS,RF,A1,A2,PTEST)
        ENDIF
       
        PREVDIH(:) = STARTDIH(:)
        CALL GETDIHONLY(RF)
      
        INTDIFF(:) = PREVDIH(:) - STARTDIH(:)
        INTDIST = SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))

        IF (PTEST) PRINT '(2i4,a6,2f15.10)', A1, A2,"dists", INTDIST, CURINTDIST
        IF (INTDIST.GE.CURINTDIST.AND..NOT.PROCHIRALH(A1)) THEN
           ! undo permutation
           RF(:) = RFTMP(:)
        ELSEIF (INTDIST.LT.CURINTDIST) THEN
           IF (PTEST) print*, 'keep permutation'
           CURINTDIST = INTDIST
        ELSEIF (PROCHIRALH(A1)) THEN
           CURINTDIST = INTDIST
        ENDIF

      ELSEIF (NPERMSIZE(PG).EQ.3) THEN
         A1 = PERMGROUP(START)
         A2 = PERMGROUP(START+1)
         A3 = PERMGROUP(START+2)

         RFTMP2(:) = RF(:)
         RFTMP(:) = RF(:)
         MININTDIST=CURINTDIST
         MINSWAP(1)=A1; MINSWAP(2)=A2; MINSWAP(3)=A3
         TMP_AR(A1,:)= RFTMP2(3*(A1-1)+1:3*A1)
         TMP_AR(A2,:)= RFTMP2(3*(A2-1)+1:3*A2)
         TMP_AR(A3,:)= RFTMP2(3*(A3-1)+1:3*A3)
         DO I=0,2
            RFTMP(3*(A1-1)+1:3*A1)=TMP_AR(A1,:)
            RFTMP(3*(A2-1)+1:3*A2)=TMP_AR(A2,:)
            RFTMP(3*(A3-1)+1:3*A3)=TMP_AR(A3,:)!initialise back to beginning
            !swap first two atoms
            B1=PERMGROUP(START+I)
            !find other two atoms
            IF (I.EQ.0) THEN !i.e. B1 = A1
               B2=PERMGROUP(START+1); B3= PERMGROUP(START+2)
            ELSEIF (I.EQ.1) THEN
               B2=A1; B3= PERMGROUP(START+2)
            ELSE
               B2=PERMGROUP(START+1); B3=A1
            ENDIF
            RFTMP(3*(A1-1)+1:3*A1)=TMP_AR(B1,:)
            RFTMP(3*(B1-1)+1:3*B1)=TMP_AR(A1,:)
            PREVDIH(:) = STARTDIH(:)
            CALL GETDIHONLY(RFTMP)
            
            INTDIFF= PREVDIH(:) - STARTDIH(:)
            INTDIST= SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))
            IF (PTEST) PRINT*,"swap", B1,B2,B3 ,"intdist", intdist

            PREVDIH(:) = STARTDIH(:)
            CALL DISTANCEPAIRSWAP(STARTDIH, RFTMP, A2,A3,PG, INTDIST2,PTEST,RFNEW)
            !PRINT*, "minintdist", MININTDIST
            IF (PTEST) PRINT*, "swap",B1, B3,B2,"intdist2", intdist2

            IF (INTDIST.LT.MININTDIST) THEN
                MININTDIST=INTDIST; MINSWAP(1)=B1;MINSWAP(2)=B2;MINSWAP(3)=B3
                IF (PTEST) PRINT*, MINSWAP(:), MININTDIST
            ENDIF
            IF (INTDIST2.LT.MININTDIST) THEN
                MININTDIST=INTDIST2; MINSWAP(1)=B1;MINSWAP(2)=B3;MINSWAP(3)=B2
                IF (PTEST) PRINT*, MINSWAP(:), MININTDIST
            ENDIF
          ENDDO
          !PRINT*,A1,A2,A3, "minswap",MINSWAP(1), MINSWAP(2), MINSWAP(3)
          RFTMP(3*(A1-1)+1:3*A1)=RFTMP2(3*(MINSWAP(1)-1)+1:3*MINSWAP(1))
          RFTMP(3*(A2-1)+1:3*A2)=RFTMP2(3*(MINSWAP(2)-1)+1:3*MINSWAP(2))
          RFTMP(3*(A3-1)+1:3*A3)=RFTMP2(3*(MINSWAP(3)-1)+1:3*MINSWAP(3))
          PREVDIH(:)=STARTDIH(:)

          CURINTDIST = MININTDIST
          RF(:) = RFTMP(:)

      ELSE
          print*, 'Error! INTMINPERM is only set up for permutation groups &
            &          of size 2 & 3 for now', PG, NPERMSIZE(PG)
          STOP
      ENDIF

      START = START + NPERMSIZE(PG)
      END DO

      IF (PTEST) print*, 'Final int distance: ', CURINTDIST
      DISTANCE = DOT_PRODUCT(RF-RS,RF-RS)
      IF (PTEST) PRINT*, 'Final cart distance', DISTANCE

      RETURN
      END SUBROUTINE INTMINPERM_CHIRAL

!****************************************************************************************
      SUBROUTINE INTDISTANCE(RS, RF, DISTANCE, PTEST)      
! if INTDISTANCET - calculates distance for minpermdist in internals for dijkstra
! this should be better when using internal coordinate interpolation
      use commons, only: natoms
      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN)    :: RS(3*NATOMS)
      DOUBLE PRECISION, INTENT(INOUT) :: RF(3*NATOMS)
      DOUBLE PRECISION, INTENT(OUT)   :: DISTANCE
      LOGICAL, INTENT(IN)             :: PTEST

      DOUBLE PRECISION :: RSINT(NINTC), RFINT(NINTC), INTDIFF(NDIH)
      DOUBLE PRECISION :: STARTDIH(NDIH)

      PREVDIH => DIHINFOSINGLE(:)
      PREVDIH(:) = 0.0D0
      CALL CART2INT(RS, RSINT)

      STARTDIH(:) = PREVDIH(:)
      ALIGNDIR = .FALSE.

      CALL GETDIHONLY(RF)
      INTDIFF(:) = PREVDIH(:)-STARTDIH(:)
      DISTANCE = SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))
      IF (PTEST) print*, 'intdistance: ', DISTANCE

      RETURN
      END SUBROUTINE INTDISTANCE

END MODULE INTCUTILS

     

! **********************************************************************************
     SUBROUTINE INTH_ALIGN(COORDSB, COORDSA, centre, CURINTDIST)
     use intcutils
     USE intcommons, only:ndih
     USE commons, only: NATOMS,DEBUG
     IMPLICIT NONE
!aligns hydrogens in valine and leucine groups according to internals
     DOUBLE PRECISION, INTENT(IN) :: COORDSB(3*NATOMS)
     DOUBLE PRECISION, INTENT(INOUT) :: COORDSA(3*NATOMS)
     INTEGER, INTENT(IN) ::centre
     DOUBLE PRECISION, INTENT(OUT):: CURINTDIST
     DOUBLE PRECISION :: TMP(3), TMP2(3), TMP3(3), STARTDIH(NDIH), INTDIFF(NDIH)
     DOUBLE PRECISION :: COORDS_copy(3*NATOMS)
     DOUBLE PRECISION :: INTDIST

     !PREVDIH => DIHINFOSINGLE(:)

!    PRINT*, "start"
!    PREVDIH(:) = 0.0D0
!    CALL GETDIHONLY(RS) !msb50
!    GD_PREVDIH(:) = PREVDIH

     !PRINT*, "centre", centre, centre+1, centre+2, centre+3
     COORDS_copy(:) = COORDSA(:)
     PREVDIH(:) = 0.0D0
     CALL GETDIHONLY(COORDSB)
     STARTDIH(:)=PREVDIH(:)
     CALL GETDIHONLY(COORDSA)
     INTDIFF(:)=PREVDIH(:)-STARTDIH(:)
     CURINTDIST = SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))
     IF (DEBUG) PRINT*, "curintdist", CURINTDIST

     TMP(:)= COORDSA(3*(centre)+1:3*(centre+1)) !for centre+1
     TMP2(:)= COORDSA(3*(centre+1)+1:3*(centre+2))
     TMP3(:)= COORDSA(3*(centre+2)+1:3*(centre+3))
     COORDSA(3*(centre)+1:3*(centre+1))=TMP2(:)
     COORDSA(3*(centre+1)+1:3*(centre+2)) =TMP3(:)
     COORDSA(3*(centre+2)+1:3*(centre+3)) =TMP(:)
     PREVDIH(:)=STARTDIH(:)
     CALL GETDIHONLY(COORDSA)
     INTDIFF(:)=PREVDIH(:)-STARTDIH(:)
     INTDIST = SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))
     IF (INTDIST.LT.CURINTDIST) THEN
        CURINTDIST=INTDIST
        COORDS_copy(:)=COORDSA(:)
        IF (DEBUG) PRINT*, "swap accepted", centre+2, centre+3, centre+1, curintdist
     ELSE
        COORDSA(:)=COORDS_copy(:)
     ENDIF
     COORDSA(3*(centre)+1:3*(centre+1))=TMP3(:)
     COORDSA(3*(centre+1)+1:3*(centre+2)) =TMP(:)
     COORDSA(3*(centre+2)+1:3*(centre+3)) =TMP2(:)
     PREVDIH(:)=STARTDIH(:)
     CALL GETDIHONLY(COORDSA)
     INTDIFF(:)=PREVDIH(:)-STARTDIH(:)
     INTDIST = SQRT(DOT_PRODUCT(INTDIFF, INTDIFF))
     IF (INTDIST.LT.CURINTDIST) THEN
        CURINTDIST=INTDIST
        COORDS_copy(:)=COORDSA(:)
        IF (DEBUG) PRINT*, "swap accepted", centre+3, centre+1, centre+2, curintdist
     ELSE
        COORDSA(:)=COORDS_copy(:)
     ENDIF
     RETURN
END SUBROUTINE INTH_ALIGN


