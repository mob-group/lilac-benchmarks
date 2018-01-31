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

! Atom-matching method for aligning periodic structures in cartesian coordinates, and identifying permutational isomers.
! Input: DUMMYB and DUMMYA are the two structures to be aligned
! Output: 
!         DUMMYB will be left unchanged. If a permutational isomer is identified, it is returned in DUMMYA and PITEST=TRUE
!         Otherwise, DUMMYA is left unchanged and PITEST=FALSE
!         If BMTEST is true, we are using the atom-matching method to find the best alignment of the two structures. In this
!         case XBEST will contain structure A with the best-identified translation and permutation applied.
SUBROUTINE BULKMINDIST(DUMMYB,DUMMYA,XBEST, NATOMS,DISTANCE,TWOD,DEBUG,BOXLX,BOXLY,BOXLZ,PITEST,RESETA, TNMATCH, BMTEST)
USE KEY,ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, GEOMDIFFTOL, ATOMMATCHFULL, BULKBOXT, GDSQ, GDSQT

IMPLICIT NONE
INTEGER J1, NATOMS, NPMIN, NGMIN, J2, PERM(NATOMS), PBEST(NATOMS), NDUMMY, NMATCHED, PATOMS, J3, J4, NMBEST, ND1, ND2
INTEGER NMATCHSV, J5, J6, LPERM(NATOMS), NREP, NREP2, J7, J8, DEGS(3*NATOMS), TEMP_PERM(NATOMS)
DOUBLE PRECISION DUMMYB(3*NATOMS),DUMMYA(3*NATOMS),DISTANCE,BOXLX,BOXLY,BOXLZ,XSHIFT,YSHIFT,ZSHIFT,XTEMP(3*NATOMS)
DOUBLE PRECISION XBEST(3*NATOMS), DMIN, DTOTAL, DIST, TEMP_DIST
DOUBLE PRECISION DIST2, LDISTANCE, WORSTRAD , DUMMYE, VNEW(3*NATOMS), RMS
LOGICAL TWOD,DEBUG,PITEST,SAMEMIN,RESETA, TNMATCH, BMTEST, BMTESTLOCAL 
COMMON /BULKSHIFT/ XSHIFT,YSHIFT,ZSHIFT
SAVE NMATCHSV

IF(BMTEST .AND. ANY(NSETS(:).GT.0)) THEN
   WRITE(*,*) "bulkmindist> Error. Atom matching method does not currently support secondary sets of permutable atoms."
   STOP
ENDIF

IF (.NOT.TNMATCH) NMATCHSV=0 
!
! Find smallest group of permutable atoms.
! Translate first atom of group to all positions and then find nearest atom within
! the same group for every other atom.
! Keep the best translation/permutation, which corresponds to the smallest
! minimum image distance.
!
!PRINT*, NMATCHSV
TNMATCH=.TRUE.
DISTANCE=1.0D100
PITEST=.FALSE.
SAMEMIN=.TRUE.
BMTESTLOCAL=BMTEST
IF (.NOT.GDSQT) THEN
   GDSQ=GEOMDIFFTOL**2/NATOMS ! because GEOMDIFFTOL is used for the total distance elsewhere in the program
ENDIF
NPMIN=HUGE(1)

DO J1=1,NPERMGROUP
   IF (NPERMSIZE(J1).LT.NPMIN) THEN
      NPMIN=NPERMSIZE(J1)
      NGMIN=J1
   ENDIF
ENDDO
ND1=0
ND2=1
DO J1=1,NGMIN-1
   ND1=ND1+NPERMSIZE(J1)  ! ND1 is an offset index for the PERMGROUP array, indexing the start of the smallest permutational group.
ENDDO
 IF (DEBUG) PRINT '(3(A,I6))',' bulkmindist> Smallest group of permutable atoms is number ',NGMIN,' with ',NPMIN,' members'

NREP=0
NREP2=0
outer: DO J1=ND1+1,ND1+NPMIN  ! Loop over the smallest permutable group
   ! Consider renaming J2 in the next line - it's confusing to use the same name for a loop index and a normal variable.
   J2=PERMGROUP(J1)  ! Get the atom index for this cycle (PERMGROUP contains the indices of the atoms in all permutable groups, in group order)
   IF(DEBUG) WRITE(*,*) "Mapping atom ", J2, "in DUMMYA onto atom ", PERMGROUP(ND1+1), "in DUMMYB "
   ! Map the current atom (J2) onto the first atom in this group (PERMGROUP(ND1+1). First, calculate the (negative of the) required shift...
   XSHIFT=DUMMYA(3*(J2-1)+1)-DUMMYB(3*(PERMGROUP(ND1+1)-1)+1)-BOXLX*NINT((DUMMYA(3*(J2-1)+1)-DUMMYB(3*(PERMGROUP(ND1+1)-1)+1))/BOXLX)
   YSHIFT=DUMMYA(3*(J2-1)+2)-DUMMYB(3*(PERMGROUP(ND1+1)-1)+2)-BOXLY*NINT((DUMMYA(3*(J2-1)+2)-DUMMYB(3*(PERMGROUP(ND1+1)-1)+2))/BOXLY)
   IF (.NOT.TWOD) ZSHIFT=DUMMYA(3*(J2-1)+3)-DUMMYB(3*(PERMGROUP(ND1+1)-1)+3)-BOXLZ*NINT((DUMMYA(3*(J2-1)+3)-DUMMYB(3*(PERMGROUP(ND1+1)-1)+3))/BOXLZ)
   ! ...and then apply this shift to all atoms in DUMMYA.
   DO J2=1,NATOMS
      XTEMP(3*(J2-1)+1)=DUMMYA(3*(J2-1)+1)-XSHIFT
      XTEMP(3*(J2-1)+2)=DUMMYA(3*(J2-1)+2)-YSHIFT
      IF (.NOT.TWOD) XTEMP(3*(J2-1)+3)=DUMMYA(3*(J2-1)+3)-ZSHIFT
   ENDDO

   ! Now, we will find the best matches between the remaining atoms in the shifted DUMMYA and the unshifted DUMMYB. If all atoms match closely, we
   ! have a translation-permutation isomer. Otherwise, we continue superimposing atoms to find the translation which gives the greatest number
   ! of exact matches of atoms. This is returned as XBEST.
   NDUMMY=1
   PERM(1:NATOMS)=-1
   NMATCHED=0
   DTOTAL=0.0D0
   permgroups: DO J2=1,NPERMGROUP  ! For each permgroup (not just the smallest)  WARNING: J2 is redefined here!
      IF(J2.GT.1) NDUMMY=NDUMMY+NPERMSIZE(J2-1)
      PATOMS=NPERMSIZE(J2)
      loop1: DO J3=1,PATOMS    ! for each atom in fixed structure B in group J2
         DMIN=1.0D100
         loop2: DO J4=1,PATOMS ! which is the closest atom in the same group for the structure in XTEMP (shifted A)?
            DIST=(XTEMP(3*(PERMGROUP(NDUMMY+J4-1)-1)+1)-DUMMYB(3*(PERMGROUP(NDUMMY+J3-1)-1)+1) &
 &  - BOXLX*NINT((XTEMP(3*(PERMGROUP(NDUMMY+J4-1)-1)+1)-DUMMYB(3*(PERMGROUP(NDUMMY+J3-1)-1)+1))/BOXLX))**2 &
            &  + (XTEMP(3*(PERMGROUP(NDUMMY+J4-1)-1)+2)-DUMMYB(3*(PERMGROUP(NDUMMY+J3-1)-1)+2) &
 &  - BOXLY*NINT((XTEMP(3*(PERMGROUP(NDUMMY+J4-1)-1)+2)-DUMMYB(3*(PERMGROUP(NDUMMY+J3-1)-1)+2))/BOXLY))**2 
            IF (.NOT.TWOD) DIST=DIST+(XTEMP(3*(PERMGROUP(NDUMMY+J4-1)-1)+3)-DUMMYB(3*(PERMGROUP(NDUMMY+J3-1)-1)+3) &
 &  - BOXLZ*NINT((XTEMP(3*(PERMGROUP(NDUMMY+J4-1)-1)+3)-DUMMYB(3*(PERMGROUP(NDUMMY+J3-1)-1)+3))/BOXLZ))**2
            IF (DIST.LT.DMIN) THEN  ! Found a new atom which is the best match in the current cycle
               DMIN=DIST
               PERM(PERMGROUP(NDUMMY+J3-1))=PERMGROUP(NDUMMY+J4-1)
               IF (DIST.LT.GDSQ) THEN  ! Is it a genuine match as defined by GEOMDIFFTOL? (See above)
                 PRINT '(A,I6,A,I6,A,G20.10)',' match found between atom ',PERMGROUP(NDUMMY+J3-1), &
 &                                            ' and ',PERMGROUP(NDUMMY+J4-1),' DIST=',DIST
                  NMATCHED=NMATCHED+1
                  DTOTAL=DTOTAL+DMIN
                  IF (PERM(PERMGROUP(NDUMMY+J3-1)).NE.PERMGROUP(NDUMMY+J3-1)) SAMEMIN=.FALSE.
                  CYCLE loop1  ! We've matched this atom in the unshifted structure, so move on to the next.
               ENDIF
            ENDIF
         ENDDO loop2  ! We've tried matching atom J3 with every permutationally related atom in the shifted structure.
!       DTOTAL=DTOTAL+DMIN
       ! If we've got this far, it means we haven't hit the CYCLE loop1 statement, and hence we did not find a match.
!       PRINT '(A,I6,A,G20.10,A,I6)',' match failed for atom ',PERMGROUP(NDUMMY+J3-1),' DMIN=',DMIN,' J1=',J1

      ! If we reached here then the atom specified in the J3 loop does not have a partner, so this choice of atoms to superimpose does not
      ! give us exactly-aligned structures.
      ! If BMTESTLOCAL is FALSE, we are not trying to find the best alignment of the two structures by atom matching (either because 
      ! ATOMMATCHDIST is not set, or has already failed). In this case, we are only interested in whether the two structures are 
      ! permutational isomers so we can skip the next section. Instead, we try the next pair of atoms to superimpose.
      IF (.NOT.BMTESTLOCAL) CYCLE outer 

      IF ((J3-NMATCHED).GT.(NATOMS-NMATCHSV)) CYCLE outer ! Match cannot be better than the previous best
      IF ((.NOT.ATOMMATCHFULL).AND.NMATCHSV.GT.0.AND.(J3-NMATCHED).GT.INT(NATOMS/2.0D0)) THEN  ! Conditions to give up on ATOMMATCHDIST
          NREP2=NREP2+1
          IF (NREP2.GE.5) THEN  
            BMTESTLOCAL=.FALSE.
            CYCLE outer ! Match is not good enough - give up
          ENDIF 
      ELSE
          NREP2=0
      ENDIF
      ENDDO loop1 ! We have finished determining the number of matches associated with this choice of atoms to superimpose
      IF (DEBUG) PRINT '(A, I6, A, I6)',' bulkmindist> number of matching atoms=', NMATCHED, ' in permgroups ', J2   
      IF (BMTESTLOCAL.AND.((NATOMS-NMATCHED).GT.0)) THEN  ! We are still looking for better matches (i.e. translation vectors which match more
                                                          ! atoms). There are still atoms left to match.
         IF (J2.EQ.NPERMGROUP.AND.NMATCHED.GE.NMATCHSV) THEN ! We have checked all permutational groups for matches, and found an improvement
                                                             ! compared to the previous best.

            ! See whether a permutational alignment improves the match and record the new best match.

            ! This obsolete call to MINPERM will find the best permutational alignment of ALL atoms, rather than each of the
            ! permutational groups.
!            CALL MINPERM(NATOMS, DUMMYB, XTEMP, BOXLX, BOXLY, BOXLZ, .TRUE., LPERM, LDISTANCE, DIST2, WORSTRAD)

            ! Initialise some variables for the permutational alignment
            LDISTANCE = 0.0D0
            TEMP_DIST = 0.0D0
            DO J7=1, NATOMS
               LPERM(J7) = J7
            ENDDO
            DO J7=1,NPERMGROUP 
               ! ND2 is the index of the start of this permutation group.
               ! PERMGROUP(ND2:ND2+NPERMSIZE(J7)-1) is an array containing the indices of all atoms belonging to this group.
               ! We build up a new array DEGS(:NPERMSIZE(J7)) which contains the indices of the corresponding degrees of freedom
               ! (i.e. indices for the coordinate arrays that correspond to this permutational group)
               DO J8=1,NPERMSIZE(J7)
                  DEGS(3*J8-2) = 3*(PERMGROUP(ND2+J8-1))-2
                  DEGS(3*J8-1) = 3*(PERMGROUP(ND2+J8-1))-1
                  DEGS(3*J8) = 3*(PERMGROUP(ND2+J8-1))
               ENDDO
 
               ! Now do the actual permutational alignment
               ! DUMMYB(DEGS(1:NPERMSIZE(J7))) extracts from DUMMYB the coordinates of this group only.
               CALL MINPERM(NPERMSIZE(J7), DUMMYB(DEGS(1:3*NPERMSIZE(J7))), XTEMP(DEGS(1:3*NPERMSIZE(J7))), &
&                           BOXLX, BOXLY, BOXLZ, .TRUE., TEMP_PERM(:NPERMSIZE(J7)), TEMP_DIST, DIST2, WORSTRAD)
               DO J8=1,NPERMSIZE(J7)
                  LPERM(PERMGROUP(ND2+J8-1)) = PERMGROUP(ND2+TEMP_PERM(J8)-1)
               ENDDO
               LDISTANCE = LDISTANCE + TEMP_DIST
               ND2 = ND2+NPERMSIZE(J7)
            ENDDO
            ND2 = 1 ! Reset for next time we go through this block

            IF(DEBUG) THEN
               TEMP_DIST = 0.0D0
               DO J7=1,NATOMS     
                  TEMP_DIST=TEMP_DIST+ (XTEMP(3*LPERM(J7)-2)-DUMMYB(3*J7-2) - BOXLX*NINT((XTEMP(3*LPERM(J7)-2)-DUMMYB(3*J7-2))/BOXLX))**2 &
&                                    + (XTEMP(3*LPERM(J7)-1)-DUMMYB(3*J7-1) - BOXLY*NINT((XTEMP(3*LPERM(J7)-1)-DUMMYB(3*J7-1))/BOXLY))**2 
                  IF (.NOT.TWOD) TEMP_DIST=TEMP_DIST+(XTEMP(3*LPERM(J7))-DUMMYB(3*J7) - BOXLZ*NINT((XTEMP(3*LPERM(J7))-DUMMYB(3*J7))/BOXLZ))**2
               ENDDO
               IF(ABS(TEMP_DIST-LDISTANCE)/LDISTANCE.GT.1.0D-4) THEN
                  WRITE(*,*) "bulkmindist> Warning: distance from MINPERM and recalculated distance do not agree:"
                  WRITE(*,*) LDISTANCE, TEMP_DIST
               ENDIF
            ENDIF

            IF (LDISTANCE.LT.DISTANCE) THEN  ! Performing the alignment improved the distance. Update the best match.
               NREP=0
               ! We used to print DIST2 as well, but it is now meaningless (it only refers to one of the permutational groups)
               IF(DEBUG) PRINT*, 'bulkmindist> Distance from minperm', LDISTANCE!, DIST2
               DISTANCE=LDISTANCE
               NMATCHSV=NMATCHED
               IF (RESETA) THEN  ! (I think this is always TRUE for this subroutine)
                  ! XBEST contains the best-so-far aligned structure (translation+permutation)
                  DO J6=1,NATOMS
                     ! XTEMP is the most recent shifted structure, LPERM is the best permutation as obtained by MINPERM
                     XBEST(3*(J6-1)+1)=XTEMP(3*(LPERM(J6)-1)+1)-BOXLX*NINT(XTEMP(3*(LPERM(J6)-1)+1)/BOXLX)
                     XBEST(3*(J6-1)+2)=XTEMP(3*(LPERM(J6)-1)+2)-BOXLY*NINT(XTEMP(3*(LPERM(J6)-1)+2)/BOXLY)
                     IF (.NOT.TWOD) XBEST(3*(J6-1)+3)=XTEMP(3*(LPERM(J6)-1)+3)-BOXLZ*NINT(XTEMP(3*(LPERM(J6)-1)+3)/BOXLZ)
                  ENDDO
               ENDIF

            ELSE IF (NMATCHED.EQ.NMATCHSV) THEN  ! We have been through a whole cycle of outer (i.e. tried a different superposition)
                                                 ! without any change in the number of matches
               NREP=NREP+1  ! Increment the fail counter. Give up when it passes 10, unless ATOMMATCHFULL is set.
               IF (NREP.GT.10) THEN
                  IF (.NOT.ATOMMATCHFULL) BMTESTLOCAL=.FALSE. 
                  write(*,*) "Giving up now; ATOMMATCHDIST just set BMTESTLOCAL to FALSE"
                  CYCLE outer ! Match is always the same - give up, still do the PI test
               ENDIF
            ENDIF
         ELSE IF (J2.EQ.NPERMGROUP) THEN  ! we have reached the end of this permutational group, no improvement on previous best
            NREP=0 
         ENDIF  
         IF (J2.EQ.NPERMGROUP.AND.NMATCHED.NE.NATOMS) CYCLE outer ! We have finished this group and not yet found a full match.
                                                                  ! Move on to the next 
      ENDIF  
   ENDDO  permgroups
   IF (SQRT(DTOTAL).GE.GEOMDIFFTOL) CYCLE outer  ! We don't have an overall permutational isomer. Try next superposition.
   IF (SAMEMIN) THEN  ! no permutation was necessary to identify the isomers: the two structures are translationally related.
      IF (DEBUG) PRINT '(A,G20.10)',' bulkmindist> identical isomers identified for distance ',SQRT(DTOTAL)
   ELSE
      IF (DEBUG) PRINT '(A,G20.10)',' bulkmindist> permutational isomers identified for distance ',SQRT(DTOTAL)
   ENDIF
   PITEST=.TRUE.   ! We have a permutational isomer. This variable gets used by MINPERMDIST on return.
   DISTANCE=DTOTAL ! Update the distance (should be 0)
   ! Now if RESETA is TRUE (which it always is), update DUMMYA to contain the coordinates of the isomer
   IF (RESETA .AND. BULKBOXT) THEN
      ! Don't put the coordinates back in the box first
      DO J2=1,NATOMS
         DUMMYA(3*(J2-1)+1)=XTEMP(3*(PERM(J2)-1)+1)
         DUMMYA(3*(J2-1)+2)=XTEMP(3*(PERM(J2)-1)+2)
         IF (.NOT.TWOD) DUMMYA(3*(J2-1)+3)=XTEMP(3*(PERM(J2)-1)+3)
      ENDDO
   ELSE IF (RESETA) THEN
      ! Do put the coordinates back in the box.
      DO J2=1,NATOMS
         DUMMYA(3*(J2-1)+1)=XTEMP(3*(PERM(J2)-1)+1)-BOXLX*NINT(XTEMP(3*(PERM(J2)-1)+1)/BOXLX)
         DUMMYA(3*(J2-1)+2)=XTEMP(3*(PERM(J2)-1)+2)-BOXLY*NINT(XTEMP(3*(PERM(J2)-1)+2)/BOXLY)
         IF (.NOT.TWOD) DUMMYA(3*(J2-1)+3)=XTEMP(3*(PERM(J2)-1)+3)-BOXLZ*NINT(XTEMP(3*(PERM(J2)-1)+3)/BOXLZ)
      ENDDO
   ENDIF

   RETURN ! Don't need to do any more if we found an isomer!
ENDDO outer  

! If we got this far, we tried all the possible superpositions and didn't find a PI isomer. XBEST should contain the best translation-permutational
! alignment of DUMMYA that we found. DUMMYB is unchanged, as usual. DUMMYA is also unchanged, in this case.
IF (DEBUG) PRINT '(A, G20.10)',' bulkmindist> distance=', DISTANCE
IF (DEBUG) PRINT '(A)',' bulkmindist> structures are not permutational isomers'

RETURN

END SUBROUTINE BULKMINDIST
!
! Apply Oh point group operation number OPNUM to coordinates in
! vector X of dimension 3*NLOCAL, returning the result in 
! vector Y.
!
SUBROUTINE OHOPS(X,Y,OPNUM,NLOCAL)
IMPLICIT NONE
INTEGER OPNUM, J2, J3, NLOCAL
DOUBLE PRECISION RMAT(3,3,48), X(3*NLOCAL), Y(3*NLOCAL)
DATA RMAT / &
 & 1.00000000000,  0,  0,   & 
 & 0,  1.00000000000,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  -1.00000000000,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & 0,  0,  1.00000000000,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  1.00000000000,  0,   & 
 & 0,  -1.00000000000,  0,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  -1.00000000000,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & 0,  0,  -1.00000000000,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  1.00000000000,  0,   & 
 & 0,  1.00000000000,  0,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  1.00000000000,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & 0,  0,  1.00000000000,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  -1.00000000000,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  -1.00000000000,  0,   & 
 & 0,  1.00000000000,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  -1.00000000000,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & 0,  -1.00000000000,  0,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & 0,  1.00000000000,  0,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  -1.00000000000,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  -1.00000000000,  0,   & 
 & 0,  1.00000000000,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  -1.00000000000,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & 0,  1.00000000000,  0,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & 0,  -1.00000000000,  0,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  1.00000000000,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  1.00000000000,  0,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  -1.00000000000,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  1.00000000000,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & 0,  -1.00000000000,  0,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & 0,  -1.00000000000,  0,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & 0,  1.00000000000,  0,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & 0,  1.00000000000,  0,   & 
 & 0,  -1.00000000000,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  1.00000000000,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  -1.00000000000,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  1.00000000000,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  -1.00000000000,  0,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & 0,  1.00000000000,  0,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  1.00000000000,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  -1.00000000000,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & 0,  1.00000000000,  0,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & 0,  -1.00000000000,  0,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & 0,  0,  -1.00000000000,   & 
 & 0,  1.00000000000,  0,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & 0,  -1.00000000000,  0,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & 0,  -1.00000000000,  0,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & 0,  1.00000000000,  0,   & 
 & 1.00000000000,  0,  0,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & 0,  1.00000000000,  0,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & 0,  1.00000000000,  0,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & 0,  -1.00000000000,  0,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  0,  -1.00000000000,   & 
 & 0,  -1.00000000000,  0,   & 
 & 0,  1.00000000000,  0,   & 
 & 1.00000000000,  0,  0,   & 
 & 0,  0,  1.00000000000,   & 
 & 0,  -1.00000000000,  0,   & 
 & -1.00000000000,  0,  0,   & 
 & 0,  0,  1.00000000000 /

IF (OPNUM.EQ.0) THEN 
   Y(1:3*NLOCAL)=X(1:3*NLOCAL)
   RETURN
ENDIF

DO J2=1,NLOCAL
   J3=3*(J2-1)
   Y(J3+1)=RMAT(1,1,OPNUM)*X(J3+1)+RMAT(1,2,OPNUM)*X(J3+2)+RMAT(1,3,OPNUM)*X(J3+3)
   Y(J3+2)=RMAT(2,1,OPNUM)*X(J3+1)+RMAT(2,2,OPNUM)*X(J3+2)+RMAT(2,3,OPNUM)*X(J3+3)
   Y(J3+3)=RMAT(3,1,OPNUM)*X(J3+1)+RMAT(3,2,OPNUM)*X(J3+2)+RMAT(3,3,OPNUM)*X(J3+3)
ENDDO

END SUBROUTINE OHOPS
