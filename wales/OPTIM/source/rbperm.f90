!
!     COORDSA becomes the optimal alignment of the optimal permutation(-inversion) isomer, but 
!     without the permutations. DISTANCE is the residual square distance for the best alignment with
!     respect to permutation(-inversion)s as well as orientation and centre of mass.
!
!     RBSITESORIENT is called first for both COORDSA and COORDSB to put them into a standard
!     orientation in DUMMYA and DUMMYB (which both have the centre of coordinates at the origin)
!     with respect to the sites via rotations corresponding to the quaternions QA and QB, 
!     respectively. The objective is to identify permutation-inversion isomers without fail.
!     However, we have to cycle over all equivalent rigid bodies in two particular orbits for DUMMYA
!     to achieve this.
!     We iterate permutations and rbmindist minimisations up to a maximum number or until no more
!     permutations are required for each instance of DUMMYA aligned according to NCHOOSE1 and 
!     NCHOOSE2 by RBSITESORIENT. The cumulative rotation that takes the initial DUMMYA to the one
!     that aligns best with DUMMYB is saved in the quaternion QCUMUL.

!     Then, if we've not going with BULK, we try again for the inverted
!     version of COORDSA. The transformation corresponding to the minimum distance
!     is saved whenever it is improved - the best alignment including permutations
!     is saved in XBEST, and the last step is to rotate this back to coincide best
!     with COORDSB (rather than DUMMYB) using QBINV which is obtained from QB. This gives suitable
!     fixed end points for DNEB.

!     Finally, we transform COORDSA to be in optimal alignment, but without the
!     permutations in XBEST. The overall transformation is
!     COORDSA -> +/- QBINV QCUMUL QA (COORDSA - CMA) 
!
!     The correspondence between COORDSA and DUMMYA after DUMMYA has been aligned by rbmindist is
!     +/- RMATCUMUL ROTA (COORDSA - CMA) = permutation(DUMMYA)
!     where +/- is given by the value of INVERT which can assume +1 or -1.
!     The centres of coordinates for COORDSA and COORDSB can be anywhere. On return, the
!     centre of coordinates of COORDSA will be the same as for COORDSB.
!
!     ----------------------------------------------------------------------------------------------
 
      SUBROUTINE RBMINPERMDIST(COORDSB,COORDSA,DISTANCE,DIST2,QBEST,RMATBEST,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,SITESB,SITESA)

!     returns DISTANCE as the actual distance, rather than the squared distance

      USE COMMONS, ONLY : NATOMS, NRBSITES 
      USE KEY,ONLY : NTSITES, NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, GEOMDIFFTOL, BESTPERM, EFIELDT, PAPT, PTSTSTT, &
     &               RBOPS, NRBGROUP, RBSYMT, DBPTDT, MULTISITEPYT, PYGPERIODICT, SANDBOXT

      IMPLICIT NONE

      INTEGER, PARAMETER :: MAXIMUMTRIES = 100
      INTEGER            :: NPERM, PATOMS, NTRIES, NSIZE, J1, J2, J3, NRB !jwrm2> unused, INFO, JMAX, LOCMAX(1)
      INTEGER            :: INVERT, NORBIT1, NORBIT2, PERM(NATOMS), NCHOOSE2, NDUMMY, LPERM(NATOMS), NCHOOSE1
      INTEGER            :: NEWPERM(NATOMS), ALLPERM(NATOMS)
      DOUBLE PRECISION   :: COORDSA(3*NATOMS), COORDSB(3*NATOMS), DISTANCE, DIST2 !jwrm2> unused, DISTWP, TEMPA(9*NATOMS)  
      DOUBLE PRECISION   :: DUMMYA(3*NATOMS), DUMMYB(3*NATOMS), DUMMY(3*NATOMS) !jwrm2> unused, DUMMYWP(3*NATOMS)
      DOUBLE PRECISION   :: XA(3*NTSITES),  XB(3*NTSITES), XTMP(3*NATOMS) !jwrm2> unused, XBS(3*NTSITES)
      DOUBLE PRECISION   :: SITESA(3*NTSITES),  SITESB(3*NTSITES)
      DOUBLE PRECISION   :: DBEST, XBEST(3*NATOMS) !jwrm2> unused, DX, DY, DZ, ENERGY, RMAT(3,3), RMATI(3,3), RMS, VNEW(3*NATOMS) 
      DOUBLE PRECISION   :: ROTA(3,3), ROTB(3,3) !jwrm2> unused, ROTINVA(3,3), ROTINVB(3,3), RMATCUMUL(3,3), 
      DOUBLE PRECISION   :: RMATBEST(3,3) !jwrm2> unused, RMATWP(3,3), ROTABEST(3,3), ROTINVBBEST(3,3)
      DOUBLE PRECISION   :: CMAX, CMAY, CMAZ, CMBX, CMBY, CMBZ
      DOUBLE PRECISION   :: PDUMMYA(3*NATOMS), PDUMMYB(3*NATOMS), LDISTANCE, XDUMMY, BOXLX, BOXLY, BOXLZ, WORSTRAD
      DOUBLE PRECISION   :: Q2(4), P(3) !jwrm2> unused, AMAT(4,4), BMAT(4,4), DIAG(4), Q(4), Q1(4)
      DOUBLE PRECISION   :: THETAH, DUMMYC(3*NATOMS) !jwrm2> unused, FCT, DUMMYD(3*NATOMS), ST, THETA
      LOGICAL            :: DEBUG, BULKT
      DOUBLE PRECISION   :: BESTA(3*NATOMS), RBDISTANCE, PVEC(3), RTEMP1(3,3), SITEDIST !jwrm2> unused, RTEMP2(3,3)
      DOUBLE PRECISION   :: QCUMUL(4), QBEST(4), QA(4), QB(4), QBINV(4), QTMP(4) ! jwrm2> unused, QI(4)
      DOUBLE PRECISION   :: COORDSBS(3*NATOMS), COORDSAS(3*NATOMS), COORDSCOMA(3*NATOMS/2), COORDSCOMB(3*NATOMS/2)
      DOUBLE PRECISION   :: ENERGY, VNEW(3*NATOMS), RMS

      IF (PAPT .OR. PTSTSTT) GOTO 100
      COORDSCOMA(:) = COORDSA(1:3*NATOMS/2)
      COORDSCOMB(:) = COORDSB(1:3*NATOMS/2)

      IF (DEBUG) THEN
         WRITE(*,'(3F20.10)')COORDSA(1:3*NATOMS)
         WRITE(*,*)
         WRITE(*,'(3F20.10)')COORDSB(1:3*NATOMS)
      ENDIF

!     do not call minpermdistrbCOM for multisite PY 
      IF(MULTISITEPYT.OR.PYGPERIODICT) GOTO 100
        NATOMS = NATOMS/2

        CALL MINPERMDISTRBCOM(COORDSCOMB,COORDSCOMA,DISTANCE,DIST2,QBEST,RMATBEST,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT)

        NATOMS = 2*NATOMS

        IF (SQRT(DISTANCE) <= GEOMDIFFTOL) THEN
         DISTANCE = SQRT(DISTANCE)
         IF (DEBUG) PRINT '(A)',' rbpermdist> minpermdistrbcom suggests identical'
         CMAX = 0.0D0; CMAY = 0.0D0; CMAZ = 0.0D0
         CMBX = 0.0D0; CMBY = 0.0D0; CMBZ = 0.0D0
         DO J1 = 1, NATOMS/2
            CMAX = CMAX + COORDSA(3*(J1-1)+1)
            CMAY = CMAY + COORDSA(3*(J1-1)+2)
            CMAZ = CMAZ + COORDSA(3*(J1-1)+3)
            CMBX = CMBX + COORDSB(3*(J1-1)+1)
            CMBY = CMBY + COORDSB(3*(J1-1)+2)
            CMBZ = CMBZ + COORDSB(3*(J1-1)+3)
         ENDDO
         CMAX = 2*CMAX/NATOMS; CMAY = 2*CMAY/NATOMS; CMAZ = 2*CMAZ/NATOMS
         CMBX = 2*CMBX/NATOMS; CMBY = 2*CMBY/NATOMS; CMBZ = 2*CMBZ/NATOMS
         DO J1 = 1, NATOMS/2
            COORDSA(3*(J1-1)+1) = COORDSA(3*(J1-1)+1) - CMAX
            COORDSA(3*(J1-1)+2) = COORDSA(3*(J1-1)+2) - CMAY
            COORDSA(3*(J1-1)+3) = COORDSA(3*(J1-1)+3) - CMAZ
            COORDSB(3*(J1-1)+1) = COORDSB(3*(J1-1)+1) - CMBX
            COORDSB(3*(J1-1)+2) = COORDSB(3*(J1-1)+2) - CMBY
            COORDSB(3*(J1-1)+3) = COORDSB(3*(J1-1)+3) - CMBZ
         ENDDO
         CALL SITEPOS(COORDSB,SITESB)
         CALL SITEPOS(COORDSA,SITESA)
         DO J1 = 1, NTSITES
            J2 = 3*J1
            SITESA(J2-2:J2) = MATMUL(RMATBEST,SITESA(J2-2:J2))
         ENDDO
        RETURN
        ENDIF

100   CONTINUE

      COORDSBS(1:3*NATOMS) = COORDSB(1:3*NATOMS) ! to trace error, see at the end
      COORDSAS(1:3*NATOMS) = COORDSA(1:3*NATOMS) ! to trace error, see at the end
!      print *, 'COORDSB'
!      print *, COORDSB
!      print *, 'COORDSA'
!      print *, COORDSA

      NRB   = (NATOMS/2)
      IF (DBPTDT.OR.SANDBOXT) THEN
         NSIZE = NTSITES
      ELSE 
         NSIZE = NRB*NRBSITES
      ENDIF
       CMAX = 0.0D0; CMAY = 0.0D0; CMAZ = 0.0D0
!      DO J1 = 1, NATOMS/2
!         CMAX = CMAX + COORDSA(3*(J1-1)+1)
!         CMAY = CMAY + COORDSA(3*(J1-1)+2)
!         CMAZ = CMAZ + COORDSA(3*(J1-1)+3)
!       ENDDO
      CMAX = 2*CMAX/NATOMS; CMAY = 2*CMAY/NATOMS; CMAZ = 2*CMAZ/NATOMS
      CMBX = 0.0D0; CMBY = 0.0D0; CMBZ = 0.0D0
      DO J1 = 1, NATOMS/2
         CMBX = CMBX + COORDSB(3*(J1-1)+1)
         CMBY = CMBY + COORDSB(3*(J1-1)+2)
         CMBZ = CMBZ + COORDSB(3*(J1-1)+3)
       ENDDO
      CMBX = 2*CMBX/NATOMS; CMBY = 2*CMBY/NATOMS; CMBZ = 2*CMBZ/NATOMS
! do not move along the z direction for gravity field (currently implemented for MULTISITEPY)
      IF(EFIELDT.AND.MULTISITEPYT) CMBZ=0.0D0
!
!     Bring COORDSB into standard orientation with respect to the site position
!     The standard orientation needs to be done for the sites if we are going to identify
!     permutation-inversion isomers with respect to the sites metric!
!
!     DUMMY(1:3*NATOMS)=DUMMYB(1:3*NATOMS)
!     CALL RBSITESORIENT(DUMMY,DUMMYB,NORBIT1,1,NORBIT2,1,NATOMS,DEBUG,ROTB,ROTINVB)
!
!     ----------------------------------------------------------------------------------------------
!     !
!     CALL POTENTIAL(COORDSB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!     PRINT '(2(A,F25.15))','before ENERGYB =',ENERGY,' RMS=',RMS
!     CALL POTENTIAL(COORDSA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!     PRINT '(2(A,F25.15))','before ENERGYA =',ENERGY,' RMS=',RMS
!!     CALL POTENTIAL(DUMMYB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!     PRINT '(2(A,F25.15))',' after ENERGYB =',ENERGY,' RMS=',RMS
!     !
!     ---------------------------------------------------------------------------------------------- 
!
!     CALL SITEPOS(DUMMYB,XB)
!     WRITE(975,'(I6)') NRB*NRBSITES
!     WRITE(975,'(A,F20.10)') 'B sites 1'
!     WRITE(975,'(A,3F20.10)') ('LA ',XB(3*(J3-1)+1:3*(J3-1)+3),J3=1,NRB*NRBSITES)

!     Coordinates in DUMMYB do not change but the coordinates in DUMMYA do.

      DUMMYB(1:3*NATOMS) = COORDSB(1:3*NATOMS)

      DUMMYC(1:3*NATOMS) = DUMMYB(1:3*NATOMS)
      CALL RBSITESORIENT(DUMMYC,DUMMY,NORBIT1,1,NORBIT2,1,NATOMS,QB,DEBUG)
      CALL QROTMAT(QB,ROTB)
      DUMMYB(1:3*NATOMS)= DUMMY(1:3*NATOMS)
      CALL SITEPOS(DUMMYB,XB)
!
      INVERT = 1

      DBEST    = 1.0D100
60    NCHOOSE1 = 0
65    NCHOOSE1 = NCHOOSE1+1
40    NCHOOSE2 = 0
30    NCHOOSE2 = NCHOOSE2+1

      DO J1 = 1, NATOMS
         ALLPERM(J1) = J1
         NEWPERM(J1) = J1
      ENDDO
 
!     The optimal alignment returned by minperdist is a local minimum, but may not
!     be the global minimum. Calling RBORIENT first should put permutational isomers
!     into a standard alignment and spot the global minimum zedro distance in one
!     go. However, we also need to cycle over equivalent atoms in orbits using NCHOOSE2.
!
!     Problems can occur if we don't use all the atoms specified by NORBIT1 and NORBIT2
!     because of the numerical cutoffs employed in MYORIENT. We could miss the
!     right orientation! 
!
!     If we use RBORIENT to produce particular orientations then we end up aligning 
!     COORDSA not with COORDSB but with the standard orientation of COORDSB in DUMMYB.
!     We now deal with this by tracking the complete transformation, including the
!     contribution of MYORIENT using ROTB and ROTINVB.
!

!      IF (INVERT.NE.1) THEN
!         PRINT '(A)',' rbminpermdist> ERROR *** inversion not yet programmed'
!         STOP
!      ENDIF


!     CALL POTENTIAL(DUMMYB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!     PRINT '(2(A,F25.15))','before DUMMYB =',ENERGY,' RMS=',RMS

!     ----------------------------------------------------------------------------------------------
!     !
!      CALL POTENTIAL(DUMMYB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!      PRINT '(2(A,F25.15))', ' ENERGYBD =',ENERGY,' RMS=',RMS
!      CALL POTENTIAL(DUMMYA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!      PRINT '(2(A,F25.15))', ' ENERGYAD =',ENERGY,' RMS=',RMS
!
!      DISTANCE = 0.D0
!      DO J1 = 1, 3*NRB
!         DISTANCE = DISTANCE + (DUMMYA(J1) - DUMMYB(J1))**2
!      ENDDO
!     !
!     ----------------------------------------------------------------------------------------------
      
!      IF (INVERT == -1) ROTA(:,:) = MATMUL(ROTA(:,:),RMATI(:,:))  
!     RMATI was returned on inversion followed by rotation to bring about best superposition on DUMMYB in one go

      DUMMYA(1:3*NATOMS) = COORDSA(1:3*NATOMS)
      DUMMYC(1:3*NATOMS) = DUMMYA(1:3*NATOMS)

      CALL RBSITESORIENT(DUMMYC,DUMMY,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,NATOMS,QA,DEBUG)
      CALL QROTMAT(QA,ROTA)
      DUMMYA(1:3*NATOMS) = DUMMY(1:3*NATOMS)
      CALL SITEPOS(DUMMYA,XA)

!      IF (INVERT == -1) THEN
!         QTMP(1:4) = QI(1:4)
!         CALL QROTQ(QA,QTMP)
!         QA(1:4) = QTMP(1:4)
!      ENDIF  

      DISTANCE = 0.D0
      DO J1 = 1, 3*NSIZE
         DISTANCE = DISTANCE + (XA(J1) - XB(J1))**2
      ENDDO

!     ----------------------------------------------------------------------------------------------
!     !
!      WRITE(975,'(I6)') NRB*NRBSITES
!      WRITE(975,'(A,2I6,A,F20.10)') 'A sites before rotation'
!      WRITE(975,'(A,3F20.10)') ('LA ',XA(3*(J3-1)+1:3*(J3-1)+3),J3=1,NRB*NRBSITES)
!      WRITE(975,'(I6)') NRB*NRBSITES
!      WRITE(975,'(A,2I6,A,F20.10)') 'B sites before rotation'
!      WRITE(975,'(A,3F20.10)') ('LA ',XB(3*(J3-1)+1:3*(J3-1)+3),J3=1,NRB*NRBSITES)
!     !
!     ----------------------------------------------------------------------------------------------

      IF (DEBUG) PRINT '(A,G20.10)',' rbminpermdist> after initial call to RBSITESORIENT distance=',SQRT(DISTANCE)
!
!     Bipartite matching routine for permutations. Coordinates in DUMMYB do not change
!     but the coordinates in DUMMYA do. DISTANCE is the distance in this case.
!     We return to label 10 after every round of permutational/orientational alignment
!     unless we have converged to the identity permutation.
!
!     Atoms are not allowed to appear in more than one group.
!     The maximum number of pair exchanges associated with a group is two.
!
      NTRIES = 0
!
!     QCUMUL is a quaternion containing the information of accumulated rotation that relates the  
!     original DUMMYA obtained from COORDSA to the final one. Initialize QCUMUL.
!
      QCUMUL(1:4) = (/1.D0, 0.D0, 0.D0, 0.D0/)

10    CONTINUE

      NTRIES = NTRIES + 1
      NDUMMY = 1

      DO J1 = 1, NATOMS
         PERM(J1) = J1
      ENDDO

!     ALLPERM saves the permutation from the previous cycle.
!     NEWPERM contains the permutation for this cycle, relative to the identity.
!     SAVEPERM is temporary storage for NEWPERM.
!     NEWPERM must be applied to ALLPERM after the loop over NPERMGROUP and corresponding swaps.

      DO J1 = 1, NPERMGROUP

         PATOMS = NPERMSIZE(J1)
         IF (PATOMS.GT.NRB) THEN
             PRINT '(A,I6,A,I6,A,I6)',' rbminpermdist> ERROR *** number of permutable sites in group ',J1,' is ',PATOMS, &
  &                                   ' which exceeds the number of rigid bodies ',NRB
             STOP
         ENDIF

         DO J2 = 1, PATOMS
            PDUMMYA(3*(J2-1)+1) = DUMMYA(3*(PERMGROUP(NDUMMY+J2-1)-1)+1)
            PDUMMYA(3*(J2-1)+2) = DUMMYA(3*(PERMGROUP(NDUMMY+J2-1)-1)+2)
            PDUMMYA(3*(J2-1)+3) = DUMMYA(3*(PERMGROUP(NDUMMY+J2-1)-1)+3)
            PDUMMYB(3*(J2-1)+1) = DUMMYB(3*(PERMGROUP(NDUMMY+J2-1)-1)+1)
            PDUMMYB(3*(J2-1)+2) = DUMMYB(3*(PERMGROUP(NDUMMY+J2-1)-1)+2)
            PDUMMYB(3*(J2-1)+3) = DUMMYB(3*(PERMGROUP(NDUMMY+J2-1)-1)+3)
         ENDDO
!
!     All permutations within this group of size NPERMSIZE(J1) are now tried.
!     Note that we are just using a metric based on the rigid body centre of mass coordinates here!
!
         CALL MINPERM(PATOMS, PDUMMYB, PDUMMYA, BOXLX, BOXLY, BOXLZ, BULKT, LPERM, LDISTANCE, DIST2, WORSTRAD)

         DO J2 = 1, PATOMS
            PERM(PERMGROUP(NDUMMY+J2-1)) = PERMGROUP(NDUMMY+LPERM(J2)-1)
         ENDDO
!
!  Will this branch ever be executed for rigid bodies?!
!
         IF (NSETS(J1).GT.0) THEN
            PRINT '(A)',' rbperm> ERROR *** associated atom swaps not allowed for rigid bodies'
            STOP
         ENDIF

         NDUMMY = NDUMMY + NPERMSIZE(J1)

      ENDDO
!
!     Due to possible swaps above, we need to update the overall permutation here.
!
      ALLPERM(1:NATOMS) = NEWPERM(1:NATOMS)
      DUMMY(1:3*NATOMS) = DUMMYA(1:3*NATOMS)
      NPERM    = 0
      DISTANCE = 0.0D0
!
!     Now permute the rotational coordinates to correspond to the centre of mass coordinate permutation.
!
      DO J1 = (NATOMS/2)+1, NATOMS
         IF (PERM(J1).NE.J1) THEN
            PRINT '(A,2I8)',' rbminpermdist> ERROR - J1,PERM(J1) should be equal', J1,PERM(J1)
         ENDIF
         PERM(J1) = PERM(J1-(NATOMS/2)) + (NATOMS/2)
      ENDDO

      DO J3 = 1, NATOMS

         DUMMYA(3*(J3-1)+1) = DUMMY(3*(PERM(J3)-1)+1)
         DUMMYA(3*(J3-1)+2) = DUMMY(3*(PERM(J3)-1)+2)
         DUMMYA(3*(J3-1)+3) = DUMMY(3*(PERM(J3)-1)+3)

         IF (J3.NE.PERM(J3)) THEN

            IF (DEBUG) WRITE(*,'(A,I5,A,I5)') ' rbpermdist> permute rigid bodies ',J3,' and ',PERM(J3)
            NPERM = NPERM + 1
            NEWPERM(J3) = ALLPERM(PERM(J3))

         ENDIF

      ENDDO

!     Site-site distance

      CALL SITEPOS(DUMMYA,XA)      ! DUMMYB remained unaltered, so as XB.

!     WRITE(975,'(I6)') NRB*NRBSITES
!     WRITE(975,'(A,F20.10)') 'A sites after MINPERM'
!     WRITE(975,'(A,3F20.10)') ('LA ',XA(3*(J3-1)+1:3*(J3-1)+3),J3=1,NRB*NRBSITES)

      DO J1 = 1, 3*NSIZE
         DISTANCE = DISTANCE + (XA(J1) - XB(J1))**2
      ENDDO
!
!     Update ALLPERM again. NEWPERM is the same as ALLPERM at the start of the next pass.
!
      ALLPERM(1:NATOMS) = NEWPERM(1:NATOMS)

      IF (DEBUG) WRITE(*,'(A,I6,A,G20.10)') ' rbminpermdist> sites distance after permuting ',NPERM,'    &
     &                                        pairs of rigid bodies=', SQRT(DISTANCE)
!
!     Further alignment. Coordinates in DUMMYA are reset by RBMINDIST (second argument).
!     Must allow at least one call to RBMINDIST in case the standard orientation result is terrible
!     but gives zero permutations!
!     We try internal symmetry operations first for each rigid body, then minimise the
!     distance further (if possible) .
!  
      IF ((NPERM.NE.0) .OR. (NTRIES.EQ.1)) THEN 
!
!     Now if RBSYMT is .TRUE. we should minimise the distance metric for all the rigid bodies
!     by considering all the allowed symmetry operations for each rigid body in turn.
!     Do this by altering the orientational coordinates in DUMMYA for each rigid body in turn,
!     after saving them, calculate the new distance using the all sites metric, and accept the
!     new orientation if the distance is lower.
!     Could generalise to more than one sort of rigid body as well.
! 
         IF (RBSYMT) THEN
!
!     This call block does not do any further alignment so as not to break the standard reference geometry.
!
            CALL SITEPOS(DUMMYB,SITESB) 
            RBDISTANCE=1.0D10
            DO J1=1,NRB
               DO J2=1,NRBGROUP
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!debug
!                  CALL SITEPOS(DUMMYA,SITESA)
!                  WRITE(975,'(I6)') NRB*NRBSITES
!                  WRITE(975,'(A,2I6,A,F20.10)') 'A sites before rotation'
!                  WRITE(975,'(A,3F20.10)') ('LA ',SITESA(3*(J3-1)+1:3*(J3-1)+3),J3=1,NRB*NRBSITES)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!debug

                  P(:)      = DUMMYA(3*(NRB+J1-1)+1:3*(NRB+J1-1)+3)
                  IF (J2 /= 1) THEN
                     CALL ROTMAT(P,RTEMP1)
                     PVEC(1:3) = MATMUL(RTEMP1,RBOPS(1:3,J2))
                     PVEC(1:3) = PVEC(1:3)/SQRT(DOT_PRODUCT(PVEC(1:3),PVEC(1:3)))
                     THETAH    = 0.5D0*RBOPS(4,J2) 
                     QTMP(1)   = COS(THETAH)
                     QTMP(2:4) = SIN(THETAH)*PVEC(1:3)
                     CALL QROTAA(QTMP,P)
                  ENDIF
                  DUMMYC(1:6*NRB) = DUMMYA(1:6*NRB)
                  DUMMYC(3*(NRB+J1-1)+1:3*(NRB+J1-1)+3) = P(1:3)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!debug
!                  CALL POTENTIAL(DUMMYC,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!                  PRINT '(A,I6,2(A,F25.15))',' rbminpermdist> after RBROT following internal symmetry  & 
!     &                                          operation ',J2,' energy for A=',  ENERGY,' RMS=',RMS
!
!     Calculate new distance based on all sites metric rather than just centre-of-mass rigid body coordinates.
! 
                  CALL SITEPOS(DUMMYC,SITESA)
                  SITEDIST=0.D0
                  DO J3=1,3*NRB*NRBSITES
                     SITEDIST=SITEDIST+(SITESA(J3)-SITESB(J3))**2
                  ENDDO
!                  WRITE(975,'(I6)') NRB*NRBSITES
!                  WRITE(975,'(A,2I6,A,F20.10)') 'C sites for J1,J2=',J1,J2,' distance=',SQRT(SITEDIST)
!                  WRITE(975,'(A,3F20.10)') ('LA ',SITESA(3*(J3-1)+1:3*(J3-1)+3),J3=1,NRB*NRBSITES)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!debug

                  IF (SITEDIST.LT.RBDISTANCE) THEN
                     RBDISTANCE=SITEDIST
                     BESTA(1:6*NRB)=DUMMYC(1:6*NRB)
                  ENDIF
!                  PRINT '(A,2F20.10)',' rbpermdist> sites distance, RBDISTANCE, =',SQRT(SITEDIST), SQRT(RBDISTANCE)
               ENDDO
               DUMMYA(1:6*NRB)=BESTA(1:6*NRB)
            ENDDO
!
!     This call aligns the overall orientation with respect to the sites metric. We should already
!     have identified permutation-inversion isomers using the alignment based on centres of mass
!     above, if applicable. We now repeat but with an orientational alignment with respect to sites.
!
            CALL RBMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,Q2,DEBUG)
 
            CALL QROTQ(Q2,QCUMUL)
            DO J1=1,NRB
               DO J2=1,NRBGROUP
!
!     MUST NOT USE RMAT HERE - because RMAT is accumulated below.
! 
                  P(:)      = DUMMYA(3*(NRB+J1-1)+1:3*(NRB+J1-1)+3)
                  IF (J2 /= 1) THEN
                     CALL ROTMAT(P,RTEMP1)
                     PVEC(1:3) = MATMUL(RTEMP1,RBOPS(1:3,J2))
                     PVEC(1:3) = PVEC(1:3)/SQRT(DOT_PRODUCT(PVEC(1:3),PVEC(1:3)))
                     THETAH    = 0.5D0*RBOPS(4,J2)
                     QTMP(1)   = COS(THETAH)
                     QTMP(2:4) = SIN(THETAH)*PVEC(1:3)
                     CALL QROTAA(QTMP,P)
                  ENDIF
                  DUMMYC(1:6*NRB) = DUMMYA(1:6*NRB)
                  DUMMYC(3*(NRB+J1-1)+1:3*(NRB+J1-1)+3) = P(1:3)
!
!                 CALL POTENTIAL(DUMMYC,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!                 PRINT '(A,I6,2(A,F25.15))',' rbminpermdist> after RBROT following internal symmetry 
!                                              operation ',J2,' energy for A=', ENERGY,' RMS=',RMS
!
!     Calculate new distance based on all sites metric rather than just centre-of-mass rigid body coordinates.
! 
                  CALL SITEPOS(DUMMYC,SITESA)
                  SITEDIST=0.D0
                  DO J3=1,3*NRB*NRBSITES
                       SITEDIST=SITEDIST+(SITESA(J3)-SITESB(J3))**2
                    ENDDO
  !                 PRINT '(A,2F20.10)',' rbminpermdist> sites distance^2,RBDISTANCE=',SITEDIST,RBDISTANCE
                    IF (SITEDIST.LT.RBDISTANCE) THEN
                       RBDISTANCE=SITEDIST
                       BESTA(1:6*NRB)=DUMMYC(1:6*NRB)
                    ENDIF
                 ENDDO
                 DUMMYA(1:6*NRB)=BESTA(1:6*NRB)
              ENDDO
           ELSE
  !
  !     This call aligns the overall orientation with respect to the sites metric.
  !     Internal symmetry operations of the rigid bodies are not considered.
  !
  !sf344 debug
  !     CALL POTENTIAL(DUMMYB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
  !     PRINT '(2(A,F25.15))','DUMMYB before RBMINDIST =',ENERGY,' RMS=',RMS
  !     CALL POTENTIAL(DUMMYA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
  !     PRINT '(2(A,F25.15))','DUMMYA before RBMINDIST 1 =',ENERGY,' RMS=',RMS

              CALL RBMINDIST(DUMMYB,DUMMYA,NATOMS,DISTANCE,Q2,DEBUG)

  !     accumulate rotation
              CALL QROTQ(Q2,QCUMUL)

           ENDIF

           DISTANCE  = DISTANCE*DISTANCE

           IF (NTRIES .LT. MAXIMUMTRIES) THEN
              GOTO 10
  !         ELSE ! prevent infinite loop
  !            IF (INVERT == -1) THEN
  !               IF(DEBUG) PRINT '(A)',' rbminpermdist> WARNING - number of tries exceeded, giving up'
  !            ELSE
  !               PRINT '(A)',' rbminpermdist> WARNING - number of tries exceeded, giving up'
  !            ENDIF
           ENDIF

        ENDIF

        IF (DISTANCE .LT. DBEST) THEN

           DBEST                = DISTANCE
           XBEST(1:3*NATOMS)    = DUMMYA(1:3*NATOMS)
           BESTPERM(1:NATOMS)   = ALLPERM(1:NATOMS)
           QTMP(1:4)  = QA(1:4)
           CALL QROTQ(QCUMUL,QTMP)
           QBEST(1:4) = QTMP(1:4)
           QBINV(1:4) = (/QB(1), -QB(2:4)/)

        ENDIF
  !
  ! If GEOMDIFFTOLL is set too sloppy we could miss the best solution by exiting via the line
  ! below without trying other orbits. Turn off this escape?!
  ! The turn off seems to be a bug for non-RBPERM runs. LJ38 tests fail with all compilers, producing
  ! a "Distance is zero: this should not happen" error!!! DJW
  ! Put the escape back in for now.

        IF (SQRT(DBEST).LT.GEOMDIFFTOL) GOTO 50
        IF (NCHOOSE2.LT.NORBIT2) GOTO 30
        IF (NCHOOSE1.LT.NORBIT1) GOTO 65
        IF (EFIELDT) GOTO 50
  !
  50    DISTANCE = DBEST       ! squared distance
  !
  !     XBEST contains the best alignment of A coordinates for the orientation of B coordinates in DUMMYB.
  !     Rotate XBEST by ROTINVBBEST to put in best correspondence with COORDSB, undoing the reorientation to  
  !     DUMMYB from RBORIENT. 
  !     We should get the same result for ROTINVBBEST * RMATBEST * (COORDSA-CMA), 
  !     where RMATBEST = +/- RMATCUMUL * ROTA for the best alignment 
  !     (aside from a possible permutation of the atom ordering)

        CALL RBROT(XBEST, XTMP, QBINV, NATOMS)
        XBEST(1:3*NATOMS) = XTMP(1:3*NATOMS)

        DO J1 = 1, (NATOMS/2)
           XBEST(3*(J1-1)+1) = XBEST(3*(J1-1)+1) + CMBX
           XBEST(3*(J1-1)+2) = XBEST(3*(J1-1)+2) + CMBY
           XBEST(3*(J1-1)+3) = XBEST(3*(J1-1)+3) + CMBZ
        ENDDO

!       CALL POTENTIAL(XBEST,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
   
  !     Site-site distance

        CALL SITEPOS(XBEST,XA)

        CALL SITEPOS(COORDSB,XB)

        XDUMMY = 0.D0
        DO J1 = 1, 3*NSIZE
           XDUMMY = XDUMMY + (XA(J1) - XB(J1))**2
        ENDDO

        IF (ABS(SQRT(XDUMMY)-SQRT(DISTANCE)).GT.GEOMDIFFTOL) THEN
            PRINT '(2(A,F20.10))','rbminpermdist> ERROR *** distance between transformed XBEST and COORDSB=',  &
       &    SQRT(XDUMMY),  ' should be ', SQRT(DISTANCE)
            PRINT *, 'COORDSB'
            DO J1 = 1, NATOMS
               J2 = 3*J1
               PRINT *, COORDSBS(J2-2), COORDSBS(J2-1), COORDSBS(J2)
            ENDDO 
            PRINT *, 'COORDSA'
            DO J1 = 1, NATOMS
               J2 = 3*J1
               PRINT *, COORDSAS(J2-2), COORDSAS(J2-1), COORDSAS(J2)
            ENDDO 
  !          PRINT *, 'XBEST'
  !          DO J1 = 1, NATOMS
  !             J2 = 3*J1
  !             PRINT *, XBEST(J2-2), XBEST(J2-1), XBEST(J2)
  !          ENDDO 

  !          PRINT *, 'RMSD=', SQRT(DISTANCE/DFLOAT(NATOMS*NRBSITES/2))
           STOP ! commented to allow newtip to continue as a temporary workaround DJW 2/12/12

        ENDIF

        CALL QROTQ(QBINV,QBEST)

        CALL QROTMAT(QBEST,RMATBEST)

        COORDSA(1:3*NATOMS)=XBEST(1:3*NATOMS) ! finally, best COORDSA should include permutations for DNEB input!

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!debug
!        CALL POTENTIAL(COORDSA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!        PRINT '(2(A,F25.15))',' finally A energy=',ENERGY,' RMS=',RMS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!debug

      DISTANCE=SQRT(DISTANCE)

      IF (DISTANCE <= GEOMDIFFTOL .AND. (PAPT .OR. PTSTSTT)) THEN
         CMAX = 0.0D0; CMAY = 0.0D0; CMAZ = 0.0D0
         CMBX = 0.0D0; CMBY = 0.0D0; CMBZ = 0.0D0
         DO J1 = 1, NATOMS/2
            CMAX = CMAX + COORDSAS(3*(J1-1)+1)
            CMAY = CMAY + COORDSAS(3*(J1-1)+2)
            CMAZ = CMAZ + COORDSAS(3*(J1-1)+3)
            CMBX = CMBX + COORDSBS(3*(J1-1)+1)
            CMBY = CMBY + COORDSBS(3*(J1-1)+2)
            CMBZ = CMBZ + COORDSBS(3*(J1-1)+3)
         ENDDO
         CMAX = 2*CMAX/NATOMS; CMAY = 2*CMAY/NATOMS; CMAZ = 2*CMAZ/NATOMS
         CMBX = 2*CMBX/NATOMS; CMBY = 2*CMBY/NATOMS; CMBZ = 2*CMBZ/NATOMS
         DO J1 = 1, NATOMS/2
            COORDSAS(3*(J1-1)+1) = COORDSAS(3*(J1-1)+1) - CMAX
            COORDSAS(3*(J1-1)+2) = COORDSAS(3*(J1-1)+2) - CMAY
            COORDSAS(3*(J1-1)+3) = COORDSAS(3*(J1-1)+3) - CMAZ
            COORDSBS(3*(J1-1)+1) = COORDSBS(3*(J1-1)+1) - CMBX
            COORDSBS(3*(J1-1)+2) = COORDSBS(3*(J1-1)+2) - CMBY
            COORDSBS(3*(J1-1)+3) = COORDSBS(3*(J1-1)+3) - CMBZ
         ENDDO
         CALL SITEPOS(COORDSBS,SITESB)
         CALL SITEPOS(COORDSAS,SITESA)
         DO J1 = 1, NTSITES
            J2 = 3*J1
            SITESA(J2-2:J2) = MATMUL(RMATBEST,SITESA(J2-2:J2))
         ENDDO
      ENDIF

      END SUBROUTINE RBMINPERMDIST

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE SITEPOS(Q, X)

!     GENERATE SITE POSITIONS FROM THE RIGID-BODY COORDINATES

      USE COMMONS, ONLY: NATOMS, NRBSITES, RBSITE
      USE KEY, ONLY: NTSITES, DBPTDT, SANDBOXT
      use sandbox_module

      IMPLICIT NONE

      INTEGER :: J1, J2, J3, J4
      DOUBLE PRECISION :: Q(3*NATOMS), X(3*NTSITES), R(3), P(3), RM(3,3)
      DOUBLE PRECISION :: RBSITETD(4,3), SIGTD(4), SIGDB(2)
      type(sandbox_molecule), pointer :: mol
      type(sandbox_site), pointer :: site
      integer :: atom_num, atom_index, site_num, mol_num
      double precision :: dummy(3, 3)

      IF (DBPTDT) CALL DEFTD(RBSITETD,SIGTD,SIGDB)

      IF (SANDBOXT) THEN
        J1=0    
        do mol_num = 1, size(molecules)
            mol => molecules(mol_num)

            call rmdrvt(q(mol%p_index: mol%p_index + 2), RM, dummy, dummy, dummy, .false.)
!            WRITE(*,*) 'size(molecules)=', size(molecules(mol_num)%sites)

            call update_sandbox_molecule(.false., mol, q(mol%r_index: mol%r_index + 2),q(mol%p_index: mol%p_index + 2))

            do atom_num=1,size(molecules(mol_num)%sites)
               J1=J1+1
               X(3*J1-2:3*J1)=molecules(mol_num)%sites(atom_num)%r
            enddo
        end do
!        WRITE(*,*) X(:)
        RETURN 

      END IF

      DO J1 = 1, NATOMS/2

         J2   = 3*J1
         R(:) = Q(J2-2:J2)
         J2   = 3*NATOMS/2 + J2
         P(:) = Q(J2-2:J2)

         CALL ROTMAT(P, RM)

         IF (DBPTDT) THEN
            IF (J1 < NATOMS/2) THEN
               DO J3 = 1, NRBSITES
                  J4          = 3*((J1-1)*NRBSITES + J3)
                  X(J4-2:J4)  = R(:) + MATMUL(RM(:,:), RBSITE(J3,:))
               ENDDO
            ELSE
               DO J3 = 1, 4
                  J4          = 3*((J1-1)*NRBSITES + J3)
                  X(J4-2:J4)  = R(:) + MATMUL(RM(:,:), RBSITETD(J3,:))
               ENDDO
            ENDIF
         ELSE
            DO J3 = 1, NRBSITES
               J4          = 3*((J1-1)*NRBSITES + J3)
               X(J4-2:J4)  = R(:) + MATMUL(RM(:,:), RBSITE(J3,:))
            ENDDO
         ENDIF
      ENDDO
     
      END SUBROUTINE SITEPOS

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBSITESORIENT(X, T1, NORBIT1, NCHOOSE1, NORBIT2, NCHOOSE2, NATOMS, Q2, DEBUG)

      USE KEY, ONLY : NTSITES, EFIELDT, MULTISITEPYT 
!                        
!     This subroutine puts permutational isomers of a rigid-body system into a standard orientation
!     with respect to the sites
!
      IMPLICIT NONE
      INTEGER          :: NATOMS, I, J, J1, J2, JMAX1, JMAX2, NORBIT1, NCHOOSE1, NORBIT2, NCHOOSE2, NSIZE !jwrm2> unused, OFFSET
      DOUBLE PRECISION :: X(3*NATOMS), T1(3*NATOMS)
!      DOUBLE PRECISION :: XS(3*NATOMS*NRBSITES/2), T1S(3*NATOMS*NRBSITES/2), T2S(3*NATOMS*NRBSITES/2), DIST(NATOMS*NRBSITES/2)
      DOUBLE PRECISION :: XS(3*NTSITES), T1S(3*NTSITES), T2S(3*NTSITES), DIST(NTSITES)
      DOUBLE PRECISION :: AX(3), Q2(4), ROTM(3,3) !jwrm2> unused, P(3), ROTMINV(3,3)
      DOUBLE PRECISION :: THETA, THETAH, FCT !jwrm2> unused, COST, COSTH, SINT, SINTH, ST
      DOUBLE PRECISION :: CMX, CMY, CMZ, DMAX, PROJ, DMAX2, CUTOFF1, DTEMP !jwrm2> unused, DUMMY
      LOGICAL          :: DEBUG
      DOUBLE PRECISION :: ENERGY, VNEW(3*NATOMS), RMS

      NSIZE = NTSITES
      !jwrm2> unused: OFFSET  = 3*NATOMS/2
      CUTOFF1 = 1.D-03

      CALL SITEPOS(X, XS)
!
!     Move centre of mass to the origin.
!
      CMX = 0.D0; CMY = 0.D0; CMZ = 0.D0
      DO I = 1, NSIZE
         J = 3*I
         CMX = CMX + XS(J-2)
         CMY = CMY + XS(J-1)
         CMZ = CMZ + XS(J)
      ENDDO

      CMX = CMX/NSIZE; CMY = CMY/NSIZE; CMZ = CMZ/NSIZE
      DO I = 1, NSIZE
         J = 3*I
         XS(J-2) = XS(J-2) - CMX
         XS(J-1) = XS(J-1) - CMY
! sf344> modification
!        IF(.NOT.EFIELDT.AND..NOT.MULTISITEPYT) 
         XS(J)   = XS(J)   - CMZ
      ENDDO

      DMAX    = -1.D0
      NORBIT1 = 1

      IF (EFIELDT) THEN
         T1S(:) = XS(:)
        
         T1(:) = X(:)
         GOTO 100
      ENDIF
!
!     Find the rbsite which is at the largest distance from the centre of mass
!
      DO J1 = 1, NSIZE

         J = 3*J1
         DIST(J1) = SQRT(XS(J-2)**2 + XS(J-1)**2 + XS(J)**2)

         IF (ABS(DIST(J1) - DMAX) < CUTOFF1) THEN

            NORBIT1 = NORBIT1 + 1
            IF (NORBIT1 == NCHOOSE1) JMAX1 = J1

         ELSE IF (DIST(J1) > DMAX) THEN

            DMAX = DIST(J1)
            NORBIT1 = 1
            JMAX1 = J1

         ENDIF

      ENDDO

!     For tagged rbsites, the choice of the first rbsite matters if it belongs to an orbit of size > 1.
!
      IF ((ABS(XS(3*JMAX1-2)) < 1.D-08) .AND. (ABS(XS(3*JMAX1-1)) < 1.D-08)) THEN 
!
!     rb JMAX1 is already on the z axis!
!
         IF (XS(3*(JMAX1-1)+3) > 0.D0) THEN

            T1S(1:3*NSIZE) = XS(1:3*NSIZE)
            Q2(1:4)   = (/1.D0, 0.D0, 0.D0, 0.D0/)   ! Identity operation       
                  
         ELSE  ! rotate about the x-axis by \pi, DO NOT INVERT!

!            P(:) = (/4.D0*DATAN(1.D0), 0.D0, 0.D0/)
!            CALL ROTMAT(P(:), ROTM(:,:))
            Q2(1:4) = (/0.D0, 1.D0, 0.D0, 0.D0/)      ! the corresponding quaternion

            DO J1 = 1, NATOMS/2
               J       = 3*J1
               T1S(J-2) = XS(J-2)
               T1S(J-1) =-XS(J-1)
               T1S(J)   =-XS(J)
            ENDDO

         ENDIF

      ELSE
!
!     Rotate all rb's through an angle THETA about the axis P(:), so as to rotate the rbsite JMAX1 onto the z axis 
!
!       THETA = DACOS(XS(3*JMAX1)/DMAX)
!     For sloppy cutoffs we cannot assume that DMAX is exactly the same for members of the same orbit!

        THETA = DACOS(XS(3*JMAX1)/DIST(JMAX1))

!     The axis is on the xy-plane and perpendicular to the projection of the translation vector of the rb JMAX1
!     on the xy-plane. The axis is so chosen that a right-handed rotation in the right-handed coordinate system
!     will result in the required transformation.

         AX(3) = 0.D0
         FCT   = DSQRT(XS(3*JMAX1-2)**2 + XS(3*JMAX1-1)**2)
         AX(1) = THETA*XS(3*JMAX1-1)/FCT
         AX(2) =-THETA*XS(3*JMAX1-2)/FCT

         CALL ROTMAT(AX(:), ROTM(:,:))
         THETAH  = 0.5D0*THETA
         AX(1:3) = AX(1:3)/MAX(SQRT(DOT_PRODUCT(AX(1:3),AX(1:3))),1.0D-10)
         Q2(1:4) = (/COS(THETAH), SIN(THETAH)*AX(1:3)/)          ! the corresponding quaternion 
          
         DO J1 = 1, NSIZE
            J2 = 3*J1
            T1S(J2-2:J2) = MATMUL(ROTM(:,:),XS(J2-2:J2))
         ENDDO   

      ENDIF

!
!     Now find the rbsite with the largest distance from the z-axis 
!
100   DMAX = -1.0D0

      DO J1 = 1, NSIZE
         J2 = 3*J1 
         DIST(J1) = SQRT(T1S(J2-2)**2 + T1S(J2-1)**2)
         IF (DIST(J1) > DMAX) DMAX = DIST(J1)
      ENDDO

      DMAX2 = -1.0D100
!
!     PROJ is the sum of the x components. Use T2S as a dummy in order not to 
!     change T1 until we have decided which atom to put in the xz plane.
!
      DO J1 = 1, NSIZE

         IF (ABS(DIST(J1) - DMAX) < CUTOFF1) THEN

            T2S(1:3*NSIZE) = T1S(1:3*NSIZE)

            CALL RBSITESROTXZ(NSIZE, J1, T2S, PROJ, DIST, Q2, .FALSE.)

            IF (ABS(PROJ - DMAX2) < CUTOFF1) THEN
               NORBIT2 = NORBIT2+1
               IF (NORBIT2 == NCHOOSE2) THEN
                  JMAX2 = J1
                  DTEMP = PROJ
               ENDIF
            ELSE IF (PROJ > DMAX2) THEN
               NORBIT2 = 1
               DMAX2   = PROJ
               DTEMP   = PROJ
               JMAX2   = J1
            ENDIF

         ENDIF

      ENDDO

!      PRINT '(A,I6,F20.10)', 'JMAX2,DMAX2=', JMAX2,DMAX2
!
!     and now rotate it into the xz plane.
!
      CALL RBSITESROTXZ(NSIZE, JMAX2, T1S, DTEMP, DIST, Q2, .TRUE.)

!     Now change the coordinates by the rotation matrix corresponding to the net transformation

      CALL RBROT(X, T1, Q2, NATOMS)

      END SUBROUTINE RBSITESORIENT

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBSITESROTXZ(NSIZE, JDO, T1S, PROJ, DIST, Q1, ROTT)

      USE KEY, ONLY : EFIELDT

      IMPLICIT NONE
      INTEGER          :: NSIZE, JDO, I, J, J2
      DOUBLE PRECISION :: T1S(3*NSIZE), PROJ, DIST(NSIZE), THETA, THETAH, COST, SINT !jwrm2> unused, ST
      DOUBLE PRECISION :: Q2(4), Q1(4), RM(3,3) !jwrm2> unused, COSTH, SINTH, FCT, P(3), ROTM(3,3)
      LOGICAL          :: ROTT

      J2     = 3*JDO
 
      IF (ABS(T1S(J2-1)) < 1.0D-8) THEN            ! already on the xz plane

         Q2(1:4) = (/1.D0, 0.D0, 0.D0, 0.D0/) 

         IF (T1S(J2-2) < 0.D0) THEN                ! rotate about the z axis by \pi, DO NOT INVERT!!

            Q2(1:4) = (/0.D0, 0.D0, 0.D0, 1.D0/) 

            DO I = 1, NSIZE
               J = 3*I
               T1S(J-2) =-T1S(J-2)
               T1S(J-1) =-T1S(J-1)
            ENDDO

         ENDIF

      ELSE

         COST   = T1S(J2-2)/DIST(JDO)       
         SINT   = T1S(J2-1)/DIST(JDO)
         THETA  =-ATAN2(SINT,COST)              ! the negative sign appears as the rotation is reversed  
         THETAH = 0.5D0*THETA
         Q2(1:4) = (/COS(THETAH), 0.D0, 0.D0, SIN(THETAH)/)     ! the quaternion corresponding to R_{z}(-\theta)
         RM(1:3,1:3) = 0.D0                           ! the rotation matrix for R_{z}(-\theta)
         RM(1,1) = COST
         RM(2,2) = COST
         RM(3,3) = 1.D0
         RM(1,2) = SINT
         RM(2,1) =-SINT

         CALL QROTMAT(Q2,RM)

         DO I = 1, NSIZE
            IF (DIST(I) /= 0.0D0) THEN
               J = 3*I
               T1S(J-2:J) = MATMUL(RM(:,:),T1S(J-2:J))
            ENDIF  
         ENDDO

      ENDIF

      IF (ROTT) THEN

         IF (EFIELDT) THEN
            Q1(1:4) = Q2(1:4)
         ELSE
!            ROTM(:,:) = MATMUL(RM(:,:),ROTM(:,:))
            CALL QROTQ(Q2,Q1)
         ENDIF

      ELSE

         PROJ = 0.D0

         DO I = 1, NSIZE

            J = 3*I
            IF (T1S(J) > 1.0D-2) PROJ = PROJ + T1S(J-2)

         ENDDO

      ENDIF

      END SUBROUTINE RBSITESROTXZ

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE INVMTRX(S, SINV)
!     inverts a 3x3 matrix S

      IMPLICIT NONE

      DOUBLE PRECISION   :: S(3,3), SINV(3,3), A(3,3), DET

!     ADJOINT

      A(1,1) = S(2,2)*S(3,3) - S(2,3)*S(3,2)
      A(1,2) = S(1,3)*S(3,2) - S(1,2)*S(3,3)
      A(1,3) = S(1,2)*S(2,3) - S(1,3)*S(2,2)
      A(2,1) = S(3,1)*S(2,3) - S(2,1)*S(3,3)
      A(2,2) = S(1,1)*S(3,3) - S(1,3)*S(3,1)
      A(2,3) = S(1,3)*S(2,1) - S(1,1)*S(2,3)
      A(3,1) = S(2,1)*S(3,2) - S(3,1)*S(2,2)
      A(3,2) = S(1,2)*S(3,1) - S(1,1)*S(3,2)
      A(3,3) = S(1,1)*S(2,2) - S(1,2)*S(2,1)

      DET    =  S(1,1)* A(1,1) + S(1,2) * A(2,1) + S(1,3) * A(3,1)

      SINV(:,:) = A(:,:)/DET

      END SUBROUTINE INVMTRX

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBROT(T, X, Q2, NATOMS)
!     takes the set of rigid-body coordinates T and returns X after rotation via the quaternion Q2 
!     about the origin
      USE KEY, ONLY : EFIELDT, MULTISITEPYT
      IMPLICIT NONE

      INTEGER          :: I, J, NATOMS
      DOUBLE PRECISION :: T(3*NATOMS), X(3*NATOMS), T1(1:3), Q2(4), P(3), RM(3,3) 
      DOUBLE PRECISION :: CMX, CMY, CMZ
!
!     Move centre of mass to the origin.
!
      CMX = 0.D0; CMY = 0.D0; CMZ = 0.D0
      DO I = 1, NATOMS/2
         J = 3*I
         CMX = CMX + T(J-2)
         CMY = CMY + T(J-1)
         CMZ = CMZ + T(J)
      ENDDO
      CMX = 2*CMX/NATOMS; CMY = 2*CMY/NATOMS; CMZ = 2*CMZ/NATOMS
      DO I = 1, NATOMS/2
         J      = 3*I
         X(J-2) = T(J-2) - CMX
         X(J-1) = T(J-1) - CMY
         IF(EFIELDT.AND.MULTISITEPYT) CMZ = 0.0D0 
         X(J)   = T(J)   - CMZ
      ENDDO

!     extract the rotation matrix corresponding to rotation via Q

      CALL QROTMAT(Q2,RM)

      DO I = 1, NATOMS/2

         J        = 3*I
         T1(1:3)  = MATMUL(RM(:,:), X(J-2:J))
         X(J-2:J) = T1(1:3)  
!     CONVERT THE ANGLE-AXIS COORDINATES

         J        = 3*NATOMS/2 + J
         P(:)     = T(J-2:J)
         CALL QROTAA(Q2,P)
         X(J-2:J) = P(1:3)

      ENDDO

      END SUBROUTINE RBROT

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE QROTAA(Q2,P)

!     transforms the angle-axis variables P under rotation via quaternion Q2 and returns 
!     the transformed variables in P

      IMPLICIT NONE

      DOUBLE PRECISION :: Q(4), Q1(4), Q2(4), P(3), THETA, THETAH, ST, FCT

      THETA   = DSQRT(DOT_PRODUCT(P,P))
!     jwrm2> We mustn't divide by THETA when the angle of rotation is zero, in which case the axis is arbitrary
      IF (THETA <= 1.D-12) THEN
         Q1 = (/1.D0, 0.D0, 0.D0, 0.D0/)
      ELSE
         THETAH  = 0.5D0*THETA
         ST      = SIN(THETAH)
         Q1(1)   = COS(THETAH)
         Q1(2:4) = P(:)*ST/THETA
      ENDIF

      Q(1)   = Q2(1)*Q1(1) - Q2(2)*Q1(2) - Q2(3)*Q1(3) - Q2(4)*Q1(4)
      Q(2)   = Q2(1)*Q1(2) + Q2(2)*Q1(1) + Q2(3)*Q1(4) - Q2(4)*Q1(3)
      Q(3)   = Q2(1)*Q1(3) + Q2(3)*Q1(1) + Q2(4)*Q1(2) - Q2(2)*Q1(4)
      Q(4)   = Q2(1)*Q1(4) + Q2(4)*Q1(1) + Q2(2)*Q1(3) - Q2(3)*Q1(2)

      THETA  = 2.D0*ACOS(Q(1))

      IF (THETA <= 1.D-12) THEN
         P (1:3) = 0.D0
      ELSE
         FCT = DSQRT(DOT_PRODUCT(Q(2:4),Q(2:4)))
         P(1:3) = THETA*Q(2:4)/FCT
      ENDIF

      END SUBROUTINE QROTAA

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE QROTQ(Q2,Q1)

!     transforms the quaternion Q1 corresponding to rotation via quaternion Q2 and returns the
!     transformed variables in Q1

      IMPLICIT NONE

      DOUBLE PRECISION :: Q(4), Q1(4), Q2(4)

      Q(1)   = Q2(1)*Q1(1) - Q2(2)*Q1(2) - Q2(3)*Q1(3) - Q2(4)*Q1(4)
      Q(2)   = Q2(1)*Q1(2) + Q2(2)*Q1(1) + Q2(3)*Q1(4) - Q2(4)*Q1(3)
      Q(3)   = Q2(1)*Q1(3) + Q2(3)*Q1(1) + Q2(4)*Q1(2) - Q2(2)*Q1(4)
      Q(4)   = Q2(1)*Q1(4) + Q2(4)*Q1(1) + Q2(2)*Q1(3) - Q2(3)*Q1(2)

      Q1(1:4) = Q(1:4)

      END SUBROUTINE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE QROTMAT(Q,RM)

!     provides the rotation matrix RM corresponding to quaternion Q
!     right-handed rotation in the right-handed coordinate system

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: Q(4)
      DOUBLE PRECISION, INTENT(OUT) :: RM(3,3)

      RM(1,1) = Q(1)**2 + Q(2)**2 - Q(3)**2 - Q(4)**2
      RM(1,2) = 2.D0*(Q(2)*Q(3) - Q(1)*Q(4))
      RM(1,3) = 2.D0*(Q(2)*Q(4) + Q(1)*Q(3))
      RM(2,1) = 2.D0*(Q(2)*Q(3) + Q(1)*Q(4))
      RM(2,2) = Q(1)**2 + Q(3)**2 - Q(2)**2 - Q(4)**2
      RM(2,3) = 2.D0*(Q(3)*Q(4) - Q(1)*Q(2))
      RM(3,1) = 2.D0*(Q(2)*Q(4) - Q(1)*Q(3))
      RM(3,2) = 2.D0*(Q(3)*Q(4) + Q(1)*Q(2))
      RM(3,3) = Q(1)**2 + Q(4)**2 - Q(2)**2 - Q(3)**2

      END SUBROUTINE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBORIENT(X, T1, NORBIT1, NCHOOSE1, NORBIT2, NCHOOSE2, NATOMS, Q2, DEBUG)
      
      USE KEY, ONLY : EFIELDT,MULTISITEPYT 
!                        
!     This subroutine puts permutational isomers of a rigid-body system into a standard orientation.
!
      IMPLICIT NONE
      INTEGER          :: NATOMS, I, J, J1, J2, JMAX1, JMAX2, NORBIT1, NCHOOSE1, NORBIT2, NCHOOSE2 !jwrm2> unused, OFFSET
      DOUBLE PRECISION :: X(3*NATOMS), DIST(NATOMS), T1(3*NATOMS), T2(3*NATOMS)
      DOUBLE PRECISION :: AX(3), Q2(4), ROTM(3,3) !jwrm2> unused, P(3), Q(4), Q1(4), RM(3,3), ROTMINV(3,3)
      DOUBLE PRECISION :: THETA, THETAH, FCT !jwrm2> unused, COST, COSTH, SINT, SINTH, ST
      DOUBLE PRECISION :: CMX, CMY, CMZ, DMAX, PROJ, DMAX2, CUTOFF1, DTEMP !jwrm2> unused, DUMMY
      LOGICAL          :: DEBUG

      !jwrm2> unused: OFFSET  = 3*NATOMS/2
      CUTOFF1 = 1.D-03
!
!     Move centre of mass to the origin.
!
      CMX = 0.D0; CMY = 0.D0; CMZ = 0.D0
      DO I = 1, NATOMS/2
         J = 3*I
         CMX = CMX + X(J-2)
         CMY = CMY + X(J-1)
         CMZ = CMZ + X(J)
      ENDDO

      CMX = 2*CMX/NATOMS; CMY = 2*CMY/NATOMS; CMZ = 2*CMZ/NATOMS
      DO I = 1, NATOMS/2
         J = 3*I
         X(J-2) = X(J-2) - CMX
         X(J-1) = X(J-1) - CMY
! do not move along the z direction for gravity field (currently implemented in MULTISITEPY)
        IF(EFIELDT.AND.MULTISITEPYT) CMZ=0.0D0
        X(J)   = X(J)   - CMZ
      ENDDO

      DMAX    = -1.D0
      NORBIT1 = 1

      IF (EFIELDT) THEN

         T1(:) = X(:)
         GOTO 100

      ENDIF
!
!     Find the rb which is at the largest distance from the centre of mass
!
      DO J1 = 1, NATOMS/2

         J = 3*J1
         DIST(J1) = SQRT(X(J-2)**2 + X(J-1)**2 + X(J)**2)

         IF (ABS(DIST(J1) - DMAX) < CUTOFF1) THEN

            NORBIT1 = NORBIT1 + 1
            IF (NORBIT1 == NCHOOSE1) JMAX1 = J1

         ELSE IF (DIST(J1) > DMAX) THEN

            DMAX = DIST(J1)
            NORBIT1 = 1
            JMAX1 = J1

         ENDIF

      ENDDO

!     For tagged rbs, the choice of the first rb matters if it belongs to an orbit of size > 1.
!
      IF ((ABS(X(3*JMAX1-2)) < 1.D-08) .AND. (ABS(X(3*JMAX1-1)) < 1.D-08)) THEN 
!
!     rb JMAX1 is already on the z axis!
!
         IF (X(3*(JMAX1-1)+3) > 0.D0) THEN

            T1(1:3*NATOMS) = X(1:3*NATOMS)
            Q2(1:4) = (/1.D0, 0.D0, 0.D0, 0.D0/)
 
         ELSE  ! rotate about the x-axis by \pi, DO NOT INVERT!

            Q2(1:4) = (/0.D0, 1.D0, 0.D0, 0.D0/)

            DO J1 = 1, NATOMS/2
               J       = 3*J1
               T1(J-2) = X(J-2)
               T1(J-1) =-X(J-1)
               T1(J)   =-X(J)
            ENDDO

         ENDIF

      ELSE
!
!     Rotate all rb's through an angle THETA about the axis P(:), so as to rotate the rb JMAX1 onto the z axis 
!
!     For sloppy cutoffs we cannot assume that DMAX is exactly the same for members of the same orbit!

        THETA = DACOS(X(3*JMAX1)/DIST(JMAX1))

!     The axis is on the xy-plane and perpendicular to the projection of the translation vector of the rb JMAX1
!     on the xy-plane. The axis is so chosen that a right-handed rotation in the right-handed coordinate system
!     will result in the required transformation.

         AX(3) = 0.D0
         FCT   = DSQRT(X(3*JMAX1-2)**2 + X(3*JMAX1-1)**2)
         AX(1) = THETA*X(3*JMAX1-1)/FCT
         AX(2) =-THETA*X(3*JMAX1-2)/FCT

         CALL ROTMAT(AX(:), ROTM(:,:))

         THETAH  = 0.5D0*THETA
         AX(1:3) = AX(1:3)/SQRT(DOT_PRODUCT(AX(1:3),AX(1:3)))
         Q2(1:4) = (/COS(THETAH), SIN(THETAH)*AX(1:3)/)          ! the corresponding quaternion

         CALL RBROT(X, T1, Q2, NATOMS)

      ENDIF

!
!     Now find the rb with the largest distance from the z-axis 
!
100   DMAX = -1.0D0

      DO J1 = 1, NATOMS/2

         J2 = 3*J1 
         DIST(J1) = SQRT(T1(J2-2)**2 + T1(J2-1)**2)
         IF (DIST(J1) > DMAX) DMAX = DIST(J1)
         
      ENDDO

      DMAX2 = -1.0D100
!
!     PROJ is the sum of the x components. Use T2 as a dummy in order not to 
!     change T1 until we have decided which atom to put in the xz plane.
!
      DO J1 = 1, NATOMS/2

         IF (ABS(DIST(J1) - DMAX) < CUTOFF1) THEN

            T2(1:3*NATOMS) = T1(1:3*NATOMS)
            CALL RBROTXZ(NATOMS, J1, T2, PROJ, DIST, Q2, .FALSE.)

            IF (ABS(PROJ - DMAX2) < CUTOFF1) THEN
               NORBIT2 = NORBIT2+1
               IF (NORBIT2 == NCHOOSE2) THEN
                  JMAX2 = J1
                  DTEMP = PROJ
               ENDIF
            ELSE IF (PROJ > DMAX2) THEN
               NORBIT2 = 1
               DMAX2   = PROJ
               DTEMP   = PROJ
               JMAX2   = J1
            ENDIF

         ENDIF

      ENDDO
!
!     and now rotate it into the xz plane.
!
!      CALL RBROTXZ(NATOMS, JMAX2, T1, DMAX2, DIST, ROTM, .TRUE.)

      CALL RBROTXZ(NATOMS, JMAX2, T1, DTEMP, DIST, Q2, .TRUE.)

!      CALL INVMTRX(ROTM(:,:), ROTMINV(:,:))

      END SUBROUTINE RBORIENT

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBROTXZ(NATOMS, JDO, T1, PROJ, DIST, Q1, ROTT)

      !jwmr2> unused: USE KEY, ONLY : EFIELDT

      IMPLICIT NONE
      INTEGER          :: NATOMS, JDO, I, J, J2 !jwrm2> unused, OFFSET
      DOUBLE PRECISION :: T1(3*NATOMS), X(3*NATOMS), PROJ, DIST(NATOMS/2), THETA, THETAH, COST, SINT !jwrm2> unused, ST
      DOUBLE PRECISION :: P(3), Q1(4), Q2(4), RM(3,3) !jwrm2> unused, CMX, CMY, CMZ, COSTH, SINTH, FCT, ROTM(3,3)
      LOGICAL          :: ROTT

      J2     = 3*JDO
      !jwrm2> unused: OFFSET = 3*NATOMS/2

      X(1:3*NATOMS) = T1(1:3*NATOMS)

      IF (ABS(T1(J2-1)) < 1.0D-8) THEN            ! already on the xz plane

         RM(:,:) = 0.D0
         RM(1,1) = 1.D0; RM(2,2) = 1.D0; RM(3,3) = 1.D0
         Q2(1:4) = (/1.D0, 0.D0, 0.D0, 0.D0/)

         IF (T1(J2-2) < 0.D0) THEN                ! rotate about the z axis by \pi, DO NOT INVERT!!

            P(:) = (/0.D0, 0.D0, 4.D0*DATAN(1.D0)/)
            CALL ROTMAT(P(:), RM(:,:))
            Q2(1:4) = (/0.D0, 0.D0, 0.D0, 1.D0/)

            DO I = 1, NATOMS/2
               J = 3*I
               X(J-2) = -T1(J-2)
               X(J-1) = -T1(J-1)
            ENDDO

         ENDIF

      ELSE

         COST   = T1(J2-2)/DIST(JDO)
         SINT   = T1(J2-1)/DIST(JDO)
         THETA  =-ATAN2(SINT,COST)              ! the negative sign appears as the rotation is reversed
         THETAH = 0.5D0*THETA
         Q2(1:4) = (/COS(THETAH), 0.D0, 0.D0, SIN(THETAH)/)
         RM(:,:) = 0.D0                           ! the rotation matrix for R_{z}(-\theta)
         RM(1,1) = COST
         RM(2,2) = COST
         RM(3,3) = 1.D0
         RM(1,2) = SINT
         RM(2,1) =-SINT

         DO I = 1, NATOMS/2
            IF (DIST(I) /= 0.0D0) THEN
               J = 3*I
               X(J-2:J) = MATMUL(RM(:,:),T1(J-2:J))
            ENDIF
         ENDDO

      ENDIF

      IF (ROTT) THEN

!jwrm2> unused: 
!         IF (EFIELDT) THEN
!            ROTM(:,:) = RM(:,:)
!         ENDIF

!     Now change the angle-axis coordinates via quaternion multiplication equivalent to rotation by ROTM

         CALL RBROT(T1, X, Q2, NATOMS)
         T1(1:3*NATOMS) = X(1:3*NATOMS)

      ELSE

         PROJ = 0.D0
         DO I = 1, NATOMS/2
            J = 3*I
            IF (X(J) > 1.0D-2) PROJ = PROJ + X(J-2)
         ENDDO

      ENDIF

      T1(1:3*NATOMS) = X(1:3*NATOMS)

      END SUBROUTINE RBROTXZ

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBFRAME(C,XS,Q,RMATBEST)

      USE COMMONS, ONLY: NATOMS, NRBSITES
      USE KEY, ONLY    : DBPT, DBPTDT, DMBLPYT, NTIPT, MSSTOCKT, PAHAT, PATCHYDT, PAPT, PTSTSTT, NTSITES, NCARBON, &
&                        NPYSITE, MULTISITEPYT, PYT, SILANET, PYGPERIODICT, SANDBOXT !jwrm2> unused, STOCKAAT
      USE pymodule, ONLY : SITECOORDS,ELLST1,ELLMAT

      IMPLICIT NONE
      
      INTEGER :: J1, J2, J3
      DOUBLE PRECISION :: XS(3*NTSITES)
      CHARACTER(LEN=132),INTENT(IN) :: C
      DOUBLE PRECISION :: Q(3*NATOMS),RMATBEST(3,3)

      IF (DBPTDT ) THEN
         WRITE(55,'(i6)') (NATOMS/2-1)*(NRBSITES-1) + 4
         WRITE(55,'(1x,a)') trim(adjustl(c))
         DO J1 = 1, NATOMS/2
            IF (J1 < NATOMS/2) THEN
               DO J2 = 1, NRBSITES - 1
                  J3 = ((J1-1)*NRBSITES + (J2-1))*3
                  IF (J2 == 1) THEN
                     WRITE(55,'(A5,1X,3F20.10)') 'O ', XS(J3+1), XS(J3+2), XS(J3+3)
                  ELSE
                     WRITE(55,'(A5,1X,3F20.10)') 'C ', XS(J3+1), XS(J3+2), XS(J3+3)
                  ENDIF
               ENDDO
            ELSE
               DO J2 = 1, 4
                  J3 = ((J1-1)*NRBSITES + (J2-1))*3
                  IF (J2 == 1) THEN
                     WRITE(55,'(A5,1X,3F20.10)') 'H ', XS(J3+1), XS(J3+2), XS(J3+3)
                  ELSEIF (J2 == 2) THEN
                     WRITE(55,'(A5,1X,3F20.10)') 'N ', XS(J3+1), XS(J3+2), XS(J3+3)
                  ELSEIF (J2 == 3) THEN
                     WRITE(55,'(A5,1X,3F20.10)') 'B ', XS(J3+1), XS(J3+2), XS(J3+3)
                  ELSE
                     WRITE(55,'(A5,1X,3F20.10)') 'S ', XS(J3+1), XS(J3+2), XS(J3+3)
                  ENDIF
               ENDDO
            ENDIF
         ENDDO

      ELSEIF (DMBLPYT .OR. DBPT) THEN

         WRITE(55,'(i6)') (NATOMS/2)*(NRBSITES-1)
         WRITE(55,'(1x,a)') trim(adjustl(c))
         DO J1 = 1, NATOMS/2
            DO J2 = 1, NRBSITES - 1
               J3 = ((J1-1)*NRBSITES + (J2-1))*3
               IF (J2 == 1) THEN
                  WRITE(55,'(A5,1X,3F20.10)') 'O ', XS(J3+1), XS(J3+2), XS(J3+3)
               ELSE
                  WRITE(55,'(A5,1X,3F20.10)') 'C ', XS(J3+1), XS(J3+2), XS(J3+3)
               ENDIF
            ENDDO
         ENDDO

      ELSEIF (NTIPT) THEN

         WRITE(55,'(i6)') (NATOMS/2)*(NRBSITES-1)
         WRITE(55,'(1x,a)') trim(adjustl(c))
         DO J1 = 1, NATOMS/2
            J2 = 4*(J1-1)*3 
            WRITE(55,'(A5,1X,3F20.10)') 'O ', XS(J2+1), XS(J2+2), XS(J2+3)
            WRITE(55,'(A5,1X,3F20.10)') 'H ', XS(J2+4), XS(J2+5), XS(J2+6)
            WRITE(55,'(A5,1X,3F20.10)') 'H ', XS(J2+7), XS(J2+8), XS(J2+9)
         ENDDO

      ELSEIF (MSSTOCKT) THEN

         WRITE(55,'(i6)') (NATOMS/2)*NRBSITES
         WRITE(55,'(1x,a)') trim(adjustl(c))
         DO J1 = 1, NATOMS/2
            DO J2 = 1, NRBSITES
               J3 = ((J1-1)*NRBSITES + (J2-1))*3
               WRITE(55,'(A5,1X,3F20.10)') 'O ', XS(J3+1), XS(J3+2), XS(J3+3)
            ENDDO
         ENDDO

      ELSEIF (PAHAT) THEN 

         WRITE(55,'(i6)') (NATOMS/2)*NRBSITES
         WRITE(55,'(1x,a)') trim(adjustl(c))
         DO J1 = 1, NATOMS/2
            DO J2 = 1, NRBSITES
               J3 = ((J1-1)*NRBSITES + (J2-1))*3
               IF (J2 <= NCARBON) THEN
                  WRITE(55,'(A5,1X,3F20.10)') 'C ', XS(J3+1), XS(J3+2), XS(J3+3)
               ELSE
                  WRITE(55,'(A5,1X,3F20.10)') 'H ', XS(J3+1), XS(J3+2), XS(J3+3)
               ENDIF
            ENDDO
         ENDDO

! jwrm2> added PAP
      ELSEIF (PAPT) THEN

         WRITE(55,'(i6)') (NATOMS/2)*(NRBSITES)
         WRITE(55,'(1x,a)') trim(adjustl(c))
         DO J1 = 1, NATOMS/2
           DO J2 = 1, NRBSITES
            J3 = ((J1-1)*NRBSITES + (J2-1))*3
            IF (J2 .EQ. 1) THEN
              WRITE(55,'(A5,1X,3F20.10)') 'O ', XS(J3+1), XS(J3+2), XS(J3+3)
            ELSEIF (J2 .LE. (NRBSITES+1)/2) THEN
              WRITE(55,'(A5,1X,3F20.10)') 'C ', XS(J3+1), XS(J3+2), XS(J3+3)
            ELSE 
              WRITE(55,'(A5,1X,3F20.10)') 'F ', XS(J3+1), XS(J3+2), XS(J3+3)
            ENDIF
          ENDDO
         ENDDO
! jwrm2> end

     ELSEIF (PTSTSTT) THEN

         WRITE(55,'(i6)') (NATOMS/2)*(NRBSITES)
         WRITE(55,'(1x,a)') trim(adjustl(c))
         DO J1 = 1, NATOMS/2
           DO J2 = 1, NRBSITES
            J3 = ((J1-1)*NRBSITES + (J2-1))*3
            IF (J2 .EQ. NRBSITES) THEN
              WRITE(55,'(A5,1X,3F20.10)') 'O ', XS(J3+1), XS(J3+2), XS(J3+3)
            ELSEIF (J2 .LE. (NRBSITES-1)/2) THEN
              WRITE(55,'(A5,1X,3F20.10)') 'C ', XS(J3+1), XS(J3+2), XS(J3+3)
            ELSE
              WRITE(55,'(A5,1X,3F20.10)') 'F ', XS(J3+1), XS(J3+2), XS(J3+3)
            ENDIF
          ENDDO
         ENDDO


      ELSEIF (PATCHYDT) THEN

         WRITE(55,'(i6)') (NATOMS/2)*(NRBSITES+1)
         WRITE(55,'(1x,a)') trim(adjustl(c))
         DO J1 = 1, NATOMS/2
            J3 = (J1-1)*3
            WRITE(55,'(A5,1X,3F20.10)') 'O ', Q(J3+1), Q(J3+2), Q(J3+3)
            DO J2 = 1, NRBSITES 
               J3 = ((J1-1)*NRBSITES + (J2-1))*3
               WRITE(55,'(A5,1X,3F20.10)') 'H ', XS(J3+1), XS(J3+2), XS(J3+3)
            ENDDO
         ENDDO

      ELSEIF (MULTISITEPYT.OR.PYGPERIODICT) THEN
         WRITE(55,'(i6)') NATOMS*NPYSITE/2
         WRITE(55,'(1x,a)') trim(adjustl(c))
             DO J2=1,NATOMS/2
                   CALL AAtoSites(Q(3*(J2-1)+1:3*(J2-1)+3),Q(3*NATOMS/2+3*(J2-1)+1:3*NATOMS/2+3*(J2-1)+3),SITECOORDS,RMATBEST)
                    DO J3=1,NPYSITE
                        WRITE(55,'(a5,2x,3f20.10,2x,a8,12f15.8,2x,a11,3f15.8)') 'O',&
                        SITECOORDS(J3,1),SITECOORDS(J3,2),SITECOORDS(J3,3),&
                        'ellipse ',ELLST1(J3,1)*2.0D0,ELLST1(J3,2)*2.0D0,ELLST1(J3,3)*2.0D0,&
                        ELLMAT(J3,1,1),ELLMAT(J3,1,2),ELLMAT(J3,1,3),ELLMAT(J3,2,1),ELLMAT(J3,2,2),ELLMAT(J3,2,3),&
                        ELLMAT(J3,3,1),ELLMAT(J3,3,2),ELLMAT(J3,3,3),&
                        'atom_vector',Q(3*NATOMS/2+3*(J2-1)+1),&
                        Q(3*NATOMS/2+3*(J2-1)+2),Q(3*NATOMS/2+3*(J2-1)+3)
                    END DO
             END DO
      ELSEIF (PYT) THEN ! sf344> implementation still buggy, use MULTISITEPY instead!
         WRITE(55,'(i6)') NATOMS*NRBSITES/2
         WRITE(55,'(1x,a)') trim(adjustl(c))
!        WRITE(*,*) 'in here'
         CALL PY_OUTPUT(55,Q)
!      ELSEIF (STOCKAAT) THEN

!         WRITE(55,'(i6)') (NATOMS/2)*NRBSITES
!         WRITE(55,'(1x,a)') trim(adjustl(c))
!         DO J1 = 1, NATOMS/2
!            J2 = (J1-1)*3
!            WRITE(55,'(A5,1X,3F20.10)') 'O ', XS(J2+1), XS(J2+2), XS(J2+3)
!         ENDDO
      ELSEIF (SANDBOXT) THEN
         WRITE(55,'(i6)') NTSITES
         WRITE(55,'(1x,a)') trim(adjustl(c))
         CALL SANDBOX_OUTPUT(55,XS)

      ELSEIF (SILANET) THEN

         WRITE(55,'(i6)') (NATOMS/2)*NRBSITES
         WRITE(55,'(1x,a)') trim(adjustl(c))
         DO J1 = 1, NATOMS/2
            DO J2 = 1, NRBSITES
               J3 = ((J1-1)*NRBSITES + (J2-1))*3
               IF (J2 == 1) THEN
                  WRITE(55,'(A5,1X,3F20.10)') 'Si ', XS(J3+1), XS(J3+2), XS(J3+3)
               ELSE
                  WRITE(55,'(A5,1X,3F20.10)') 'H ', XS(J3+1), XS(J3+2), XS(J3+3)
               ENDIF
            ENDDO
         ENDDO


      ENDIF

      END SUBROUTINE RBFRAME

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE STFRAME(C,XS,XV)

      USE KEY, ONLY    : NTSITES

      IMPLICIT NONE
      
      INTEGER :: J1, J2
      CHARACTER(LEN=132),INTENT(IN) :: C
      DOUBLE PRECISION, INTENT(IN)  :: XS(3*NTSITES), XV(3*NTSITES)

      WRITE(55,'(i6)') NTSITES
      WRITE(55,'(1x,a)') trim(adjustl(c))

      DO J1 = 1, NTSITES
         J2 = J1*3
         WRITE(55,'(A5,1X,3F16.10,2X, A11,2X,3F16.10)') 'O ', XS(J2-2), XS(J2-1), XS(J2), 'atom_vector', XV(J2-2), XV(J2-1), XV(J2)

      ENDDO

      END SUBROUTINE STFRAME

!     ----------------------------------------------------------------------------------------------
      SUBROUTINE PYFRAME(C,XS,XV)

      USE KEY, ONLY    : NTSITES,PYA1BIN
      USE COMMONS, ONLY : NATOMS

      IMPLICIT NONE
      
      INTEGER :: J1, J2
      CHARACTER(LEN=132),INTENT(IN) :: C
      DOUBLE PRECISION, INTENT(IN)  :: XS(3*NATOMS), XV(3*NATOMS/2)
      DOUBLE PRECISION :: EulerPhiDeg, EulerPsiDeg, EulerThetaDeg

      WRITE(50,'(I6)') (NATOMS/2)
      WRITE(50,'(1X,A)') TRIM(ADJUSTL(C))

      DO J1 = 1, NATOMS/2
         J2 = J1*3
!             CALL AAtoEuler(XS(J2-2), XS(J2-1), XS(J2),EulerPhiDeg,EulerPsiDeg,EulerThetaDeg)
                 WRITE(50,'(a5,2x,3f20.10,2x,a8,6f15.8,2x,a11,3f15.8)') 'O',XS(J2-2), XS(J2-1), XS(J2),  'ellipse ',&
                        PYA1BIN(J1,1)*2.0D0,PYA1BIN(J1,2)*2.0D0,PYA1BIN(J1,3)*2.0D0,0.0,0.0,0.0,&
                      'atom_vector',XS(3*NATOMS/2 + J2-2), XS(3*NATOMS/2 + J2-1), XS(3*NATOMS/2 + J2)
      ENDDO

      END SUBROUTINE PYFRAME

      SUBROUTINE RBROTMAT(T, X, ROTMAT, NTSITES)
!     takes the set of rigid-body coordinates T and returns X after rotation via the quaternion Q2
!     about the origin
      USE KEY, ONLY : EFIELDT, MULTISITEPYT
      IMPLICIT NONE

      INTEGER          :: I, J, NTSITES
      DOUBLE PRECISION :: T(3*NTSITES), X(3*NTSITES), T1(1:3), ROTMAT(3,3) !jwmr2> unused, Q2(4)
      DOUBLE PRECISION :: CMX, CMY, CMZ
!
!     Move centre of mass to the origin.
!
      CMX = 0.D0; CMY = 0.D0; CMZ = 0.D0
      DO I = 1, NTSITES
         J = 3*I
         CMX = CMX + T(J-2)
         CMY = CMY + T(J-1)
         CMZ = CMZ + T(J)
      ENDDO
      CMX = CMX/NTSITES; CMY = CMY/NTSITES; CMZ = CMZ/NTSITES
      DO I = 1, NTSITES
         J      = 3*I
         X(J-2) = T(J-2) - CMX
         X(J-1) = T(J-1) - CMY
! do not move in the z direction for gravity field (currently implemented in MULTISITEPY)
        IF(EFIELDT.AND.MULTISITEPYT) CMZ=0.0D0 
        X(J)   = T(J)   - CMZ
      ENDDO

!     extract the rotation matrix corresponding to rotation via Q

!      CALL QROTMAT(Q2,RM)

      DO I = 1, NTSITES

         J        = 3*I
         T1(1:3)  = MATMUL(ROTMAT(:,:), X(J-2:J))
         X(J-2:J) = T1(1:3)

      ENDDO

      END SUBROUTINE RBROTMAT

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBREFLECTYZ(X)
      USE COMMONS, ONLY : NATOMS
      USE VEC3
      USE ROTATIONS
      IMPLICIT NONE
      INTEGER          :: J1,J2
      DOUBLE PRECISION :: X(3*NATOMS),M(3,3),I(3,3),A(3),P(3)

      CALL IDENTITY3X3(I(:,:))

      DO J1=1,NATOMS/2
! reflect COM coordinates
         J2=3*J1
         X(J2-2)=(-1.D0)*X(J2-2)
! reflect AA coordinates
         J2=3*(J1+NATOMS/2)
         M(:,:)=ROT_AA2MX(X(J2-2:J2))
         A(:)=VEC_CROSS( MATMUL(M(:,:),I(2,:)), MATMUL(M(:,:),I(3,:)) )
         P(:)=2.D0*ACOS(A(1)/VEC_LEN(A(:)))*(/0.D0,A(3),-1.D0*A(2)/)/SQRT(A(2)**2+A(3)**2)
         X(J2-2:J2)=ROT_ROTATE_AA(X(J2-2:J2),P(:))
      ENDDO

      END SUBROUTINE RBREFLECTYZ

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBPATHFRAME(Q,QE)

      USE COMMONS, ONLY: NATOMS, NRBSITES
      USE KEY, ONLY    : PAHAT, NTSITES, NCARBON

      IMPLICIT NONE

      INTEGER          :: J1, J2, J3
      DOUBLE PRECISION :: Q(3*NATOMS), XS(3*NTSITES), QE    
!      CHARACTER(LEN=*) ITSTRING


!      OPEN(UNIT=3,FILE=ITSTRING,STATUS='UNKNOWN')
      WRITE(3,'(I6)') NTSITES
      WRITE(3,*) 'Energy =', QE
      CALL SITEPOS(Q,XS)

      IF (PAHAT) THEN
         DO J1 = 1, NATOMS/2
            DO J2 = 1, NRBSITES
               J3 = ((J1-1)*NRBSITES + (J2-1))*3
               IF (J2 <= NCARBON) THEN
                  WRITE(3,'(A5,1X,3F20.10)') 'C ', XS(J3+1), XS(J3+2), XS(J3+3)
               ELSE
                  WRITE(3,'(A5,1X,3F20.10)') 'H ', XS(J3+1), XS(J3+2), XS(J3+3)
               ENDIF
            ENDDO
         ENDDO
      ENDIF

      END SUBROUTINE RBPATHFRAME
