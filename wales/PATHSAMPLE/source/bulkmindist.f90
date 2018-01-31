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
!
! This subroutine should match permutational isomers that correspond to the same orientation
! of the unit cell. We would have to consider all the cell symmetries as well to
! identify all possible equivalent permutation-inversion isomers.
!
SUBROUTINE BULKMINDIST(DUMMYB,DUMMYA,XBEST, NATOMS,DISTANCE,TWOD,DEBUG,BOXLX,BOXLY,BOXLZ,PITEST,RESETA, TNMATCH, BMTEST)
USE COMMONS,ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP, SETS, NSETS, GEOMDIFFTOL, ATOMMATCHFULL

IMPLICIT NONE
INTEGER J1, NATOMS, NPMIN, NGMIN, J2, PERM(NATOMS), PBEST(NATOMS), NDUMMY, NMATCHED, PATOMS, J3, J4, NMBEST, ND1
INTEGER NMATCHSV, J5, J6, LPERM(NATOMS), SAMEPERM(NATOMS), NEWPERM(NATOMS), NREP, NREP2
DOUBLE PRECISION DUMMYB(3*NATOMS),DUMMYA(3*NATOMS),DISTANCE,BOXLX,BOXLY,BOXLZ,XSHIFT,YSHIFT,ZSHIFT,XTEMP(3*NATOMS)
DOUBLE PRECISION XBEST(3*NATOMS), DMIN, DTOTAL, DIST, GDSQ, SAVECOORDS(3)
DOUBLE PRECISION DIST2, LDISTANCE, WORSTRAD 
LOGICAL TWOD,DEBUG,PITEST,SAMEMIN,RESETA, TNMATCH, BMTEST, BMTESTLOCAL 
COMMON /BULKSHIFT/ XSHIFT,YSHIFT,ZSHIFT
SAVE NMATCHSV

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
!GDSQ=GEOMDIFFTOL**2
GDSQ=GEOMDIFFTOL**2/NATOMS ! because GEOMDIFFTOL is usually used for the 3N dimensional space
BMTESTLOCAL=BMTEST

NPMIN=HUGE(1)
! PRINT *,'DUMMYA in bulkmindist:'
! PRINT '(3G20.10)',DUMMYA(1:3*NATOMS)
! PRINT *,'DUMMYB in bulkmindist:'
! PRINT '(3G20.10)',DUMMYB(1:3*NATOMS)
! Identify the smallest group of permutable atoms
DO J1=1,NPERMGROUP
   IF (NPERMSIZE(J1).LT.NPMIN) THEN
      NPMIN=NPERMSIZE(J1)
      NGMIN=J1
   ENDIF
ENDDO
ND1=0
! ND1 is an offset index for the start of this permutational group
DO J1=1,NGMIN-1
   ND1=ND1+NPERMSIZE(J1)
ENDDO
IF (DEBUG) PRINT '(3(A,I6))',' bulkmindist> Smallest group of permutable atoms is number ',NGMIN,' with ',NPMIN,' members'
NREP=0
NREP2=0

! Cycle through all atoms in the smallest group and superimpose their position in structure A over the position of the
! first atom of the same group in structure B.
outer: DO J1=ND1+1,ND1+NPMIN
!  USED(1:NATOMS)=.FALSE.
   J2=PERMGROUP(J1)

   ! Perform the superposition
   XSHIFT=DUMMYA(3*(J2-1)+1)-DUMMYB(3*(PERMGROUP(ND1+1))+1)-BOXLX*NINT((DUMMYA(3*(J2-1)+1)-DUMMYB(3*(PERMGROUP(ND1+1))+1))/BOXLX)
   YSHIFT=DUMMYA(3*(J2-1)+2)-DUMMYB(3*(PERMGROUP(ND1+1))+2)-BOXLY*NINT((DUMMYA(3*(J2-1)+2)-DUMMYB(3*(PERMGROUP(ND1+1))+2))/BOXLY)
   IF (.NOT.TWOD) ZSHIFT=DUMMYA(3*(J2-1)+3)-DUMMYB(3*(PERMGROUP(ND1+1))+3)-BOXLZ*NINT((DUMMYA(3*(J2-1)+3)-DUMMYB(3*(PERMGROUP(ND1+1))+3))/BOXLZ)
   DO J2=1,NATOMS
      XTEMP(3*(J2-1)+1)=DUMMYA(3*(J2-1)+1)-XSHIFT
      XTEMP(3*(J2-1)+2)=DUMMYA(3*(J2-1)+2)-YSHIFT
      IF (.NOT.TWOD) XTEMP(3*(J2-1)+3)=DUMMYA(3*(J2-1)+3)-ZSHIFT
   ENDDO
   NDUMMY=1

   ! Three arrays are used to keep track of the permutations.
   ! SAMEPERM is used to identify when two structures are identical up to a simple translation. Initially, all elements
   ! are set to -1, then each element is filled in with the index of the atom in XTEMP which best matches the corresponding
   ! atom in DUMMYB. If any non-trivial permutation is identified (i.e. SAMEPERM(J).NE.J for any J) then SAMEMIN is set
   ! to FALSE, i.e. we do not have identical structures.
   SAMEPERM(1:NATOMS)=-1
   ! NEWPERM is changed whenever a permutational match is found between an atom in DUMMYB and an atom in XTEMP, and unlike
   ! SAMEPERM it also takes into account secondary sets of permutable atoms which must be permuted every time the primary
   ! atom is permuted.
   ! PERM is used to calculate distances between atoms (see below). It is copied from NEWPERM at the end of each loop
   ! over a permutational group, because otherwise some atom indices get completely overwritten which then makes matching
   ! those atoms impossible.
   DO J2 = 1,NATOMS
       NEWPERM(J2) = J2
       PERM(J2) = J2
   ENDDO
   NMATCHED=0
   DTOTAL=0.0D0

   DO J2=1,NPERMGROUP  ! Cycle through permutational groups
      IF(J2.GT.1) NDUMMY=NDUMMY+NPERMSIZE(J2-1)
      PATOMS=NPERMSIZE(J2)
      loop1: DO J3=1,PATOMS    ! for each atom in fixed structure B that is in group J2
         DMIN=1.0D100
         loop2: DO J4=1,PATOMS ! which is the closest atom in the same group for the structure in XTEMP (shifted A)?
!           IF (USED(PERMGROUP(NDUMMY+J4-1)) CYCLE loop2 ! to prevent false matches with large distance cutoffs

! When we have secondary sets, we need to ensure that these atoms are permuted whenever the corresponding primary atom
! is permuted, and that when we try to match these secondary atoms between the two structures we compare their
! permuted positions rather than their original positions otherwise we won't find the correct match.
! So we apply PERM to the XTEMP indices.
            DIST=(XTEMP(3*(PERM(PERMGROUP(NDUMMY+J4-1))-1)+1)-DUMMYB(3*(PERMGROUP(NDUMMY+J3-1)-1)+1) &
 &  - BOXLX*NINT((XTEMP(3*(PERM(PERMGROUP(NDUMMY+J4-1))-1)+1)-DUMMYB(3*(PERMGROUP(NDUMMY+J3-1)-1)+1))/BOXLX))**2 &
            &  + (XTEMP(3*(PERM(PERMGROUP(NDUMMY+J4-1))-1)+2)-DUMMYB(3*(PERMGROUP(NDUMMY+J3-1)-1)+2) &
 &  - BOXLY*NINT((XTEMP(3*(PERM(PERMGROUP(NDUMMY+J4-1))-1)+2)-DUMMYB(3*(PERMGROUP(NDUMMY+J3-1)-1)+2))/BOXLY))**2
            IF (.NOT.TWOD) DIST=DIST+(XTEMP(3*(PERM(PERMGROUP(NDUMMY+J4-1))-1)+3)-DUMMYB(3*(PERMGROUP(NDUMMY+J3-1)-1)+3) &
 &  - BOXLZ*NINT((XTEMP(3*(PERM(PERMGROUP(NDUMMY+J4-1))-1)+3)-DUMMYB(3*(PERMGROUP(NDUMMY+J3-1)-1)+3))/BOXLZ))**2

            IF (DIST.LT.DMIN) THEN
            ! Remember, only if PERMGROUP(NDUMMY+J3-1)=PERMGROUP(NDUMMY+J4-1) for all J3, J4 do we have identical isomers
                 SAMEPERM(PERMGROUP(NDUMMY+J3-1))=PERMGROUP(NDUMMY+J4-1)

                 DMIN=DIST

                 IF (DIST.LT.GDSQ) THEN
!                    PRINT '(A,I6,A,I6,A,G20.10)',' match found between atom ',PERMGROUP(NDUMMY+J3-1), &
! &                                            ' and ',PERMGROUP(NDUMMY+J4-1),' DIST=',DIST

                    ! NEWPERM only holds "exact" matches
                    NEWPERM(PERMGROUP(NDUMMY+J3-1))=PERMGROUP(NDUMMY+J4-1)
                    ! sn402: If there are any secondary groups attached to the atom we have just permuted, we should
                    ! permute these as well.
                    IF (NSETS(J2).GT.0) THEN
                        DO J5=1,NSETS(J2)
                            NEWPERM(SETS(PERMGROUP(NDUMMY+J3-1),J5))=SETS(PERMGROUP(NDUMMY+J4-1),J5)
                        ENDDO
                    ENDIF

                   NMATCHED=NMATCHED+1

                   DTOTAL=DTOTAL+DMIN
                   IF (SAMEPERM(PERMGROUP(NDUMMY+J3-1)).NE.PERMGROUP(NDUMMY+J3-1)) SAMEMIN=.FALSE.
                   CYCLE loop1
                 ENDIF
            ENDIF
         ENDDO loop2
!       DTOTAL=DTOTAL+DMIN
!       PRINT '(A,I6,A,G20.10,A,I6)',' match failed for atom ',PERMGROUP(NDUMMY+J3-1),' DMIN=',DMIN,' J1=',J1

      IF (.NOT.BMTESTLOCAL) CYCLE outer ! If we reached here then we don't have a permutational isomer because
                      ! the atom specified in the J3 loop does not have a partner.

      IF ((J3-NMATCHED).GT.(NATOMS-NMATCHSV)) CYCLE outer ! Match cannot be better than the previous best
      IF ((.NOT.ATOMMATCHFULL).AND.NMATCHSV.GT.0.AND.(J3-NMATCHED).GT.INT(NATOMS/2.0D0)) THEN
          NREP2=NREP2+1
          IF (NREP2.GE.5) THEN  
            BMTESTLOCAL=.FALSE.
            CYCLE outer ! Match is not good enough - give up
          ENDIF 
      ELSE
          NREP2=0
      ENDIF
      ENDDO loop1 

      IF (DEBUG) PRINT '(A, I6)',' bulkmindist> number of matching atoms=', NMATCHED   
      IF (BMTESTLOCAL) THEN
        IF (J2.EQ.NPERMGROUP.AND.NMATCHED.GE.NMATCHSV) THEN ! We have cycled over all atoms. Record best match.
            CALL MINPERM(NATOMS, DUMMYB, XTEMP, BOXLX, BOXLY, BOXLZ, .TRUE., LPERM, LDISTANCE, DIST2, WORSTRAD)
            IF (LDISTANCE.LT.DISTANCE) THEN
                NREP=0
!               PRINT*, 'From minperm', LDISTANCE, DIST2
                DISTANCE=LDISTANCE
                NMATCHSV=NMATCHED
                IF (RESETA) THEN
                    DO J6=1,NATOMS
                        XBEST(3*(J6-1)+1)=XTEMP(3*(LPERM(J6)-1)+1)-BOXLX*NINT(XTEMP(3*(LPERM(J6)-1)+1)/BOXLX)
                        XBEST(3*(J6-1)+2)=XTEMP(3*(LPERM(J6)-1)+2)-BOXLY*NINT(XTEMP(3*(LPERM(J6)-1)+2)/BOXLY)
                        IF (.NOT.TWOD) XBEST(3*(J6-1)+3)=XTEMP(3*(LPERM(J6)-1)+3)-BOXLZ*NINT(XTEMP(3*(LPERM(J6)-1)+3)/BOXLZ)
                    ENDDO
                ENDIF
            ELSE IF (NMATCHED.EQ.NMATCHSV) THEN
                NREP=NREP+1
                IF (NREP.GT.10) THEN
                    IF (.NOT.ATOMMATCHFULL) BMTESTLOCAL=.FALSE.
                    CYCLE outer ! Match is always the same - give up, still do the PI test
                ENDIF
            ENDIF
        ELSE IF (J2.EQ.NPERMGROUP) THEN
            NREP=0
        ENDIF
        IF (J2.EQ.NPERMGROUP.AND.NMATCHED.NE.NATOMS) CYCLE outer
      ENDIF  

      ! sn402: Update the PERM arrays to hold the permutations identified for this group. NEWPERM should now contain
      ! every atom index once and once only (we can't guarantee that that is true before this point) and therefore is
      ! ready to start being used for distance calculations.
      SAMEPERM(:) = NEWPERM(:)
      PERM(:) = NEWPERM(:)

   ENDDO ! End of loop over J2
   IF (SAMEMIN) THEN
      IF (DEBUG) PRINT '(A,G20.10)',' bulkmindist> identical isomers identified, distance=',SQRT(DTOTAL)
   ELSE
      IF (DEBUG) PRINT '(A,G20.10)',' bulkmindist> permutational isomers identified, distance=',SQRT(DTOTAL)
   ENDIF
   PITEST=.TRUE.
   DISTANCE=DTOTAL
   IF (RESETA) THEN
      DO J2=1,NATOMS
         DUMMYA(3*(J2-1)+1)=XTEMP(3*(PERM(J2)-1)+1)-BOXLX*NINT(XTEMP(3*(PERM(J2)-1)+1)/BOXLX)
         DUMMYA(3*(J2-1)+2)=XTEMP(3*(PERM(J2)-1)+2)-BOXLY*NINT(XTEMP(3*(PERM(J2)-1)+2)/BOXLY)
         IF (.NOT.TWOD) DUMMYA(3*(J2-1)+3)=XTEMP(3*(PERM(J2)-1)+3)-BOXLZ*NINT(XTEMP(3*(PERM(J2)-1)+3)/BOXLZ)
      ENDDO
   ENDIF

   RETURN
ENDDO outer

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
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  -1.00000000000D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & 0.0D0,  1.00000000000D0,  0.0D0,   & 
 & 1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0,   & 
 & 0.0D0,  -1.00000000000D0,  0.0D0,   & 
 & -1.00000000000D0,  0.0D0,  0.0D0,   & 
 & 0.0D0,  0.0D0,  1.00000000000D0 /

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
