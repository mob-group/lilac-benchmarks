!     Copyright (C) 1999-2008 David J. Wales
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
!  This routine uses local optimal alignment for each group of permutable atoms.
!  It is intended for use with CHARMM and AMBER.
!  Overall alignment is based on the transformation for the best preserved local group.
!
!  COORDSA becomes the optimal alignment of the optimal permutation
!  isomer, but without the permutations. DISTANCE is the residual square distance
!  for the best alignment with respect to permutation as well as
!  orientation and centre of mass.
!
!  The centres of coordinates for COORDSA and COORDSB can be anywhere. On return, the
!  centre of coordinates of COORDSA will be the same as for COORDSB.
!
SUBROUTINE LOPERMDIST(COORDSB,COORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGID,RMATBEST)
USE COMMONS,ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, SETS, PULLT, &
  &            LOCALPERMCUT, LOCALPERMNEIGH, LOCALPERMCUT2
IMPLICIT NONE

INTEGER, PARAMETER :: MAXIMUMTRIES=10
INTEGER NATOMS, NPERM, PATOMS, NORBIT1, NORBIT2, NCHOOSE2, NCHOOSE1, NORBITB1, NORBITB2, BESTPERM(NATOMS), NDMEAN
INTEGER J3, J4, NDUMMY, LPERM(NATOMS), J1, J2, NOTHER, LPERMBEST(NATOMS), NCHOOSEB1, NCHOOSEB2, &
        LPERMBESTATOM(NATOMS)
DOUBLE PRECISION DIST2, COORDSA(3*NATOMS), COORDSB(3*NATOMS), DISTANCE, DUMMYA(3*NATOMS), &
  &              DUMMYB(3*NATOMS), DUMMY(3*NATOMS), DSUM
DOUBLE PRECISION BOXLX,BOXLY,BOXLZ,WORSTRAD,RMAT(3,3), DBEST, XBEST(3*NATOMS)
DOUBLE PRECISION ROTA(3,3), ROTINVA(3,3), ROTB(3,3), ROTINVB(3,3), RMATBEST(3,3)
LOGICAL DEBUG, TWOD, RIGID, BULKT, ADDED, PERMUTABLE(NATOMS)
DOUBLE PRECISION PDUMMYA(3*NATOMS), PDUMMYB(3*NATOMS), LDISTANCE, XDUMMY, &
   &             LDBEST(NPERMGROUP), LDBESTATOM
DOUBLE PRECISION SPDUMMYA(3*NATOMS), SPDUMMYB(3*NATOMS)
INTEGER NEWPERM(NATOMS), ALLPERM(NATOMS), SAVEPERM(NATOMS)
CHARACTER(LEN=5) ZSYMSAVE
COMMON /SYS/ ZSYMSAVE
DOUBLE PRECISION XA, XB, YA, YB, ZA, ZB, DMEAN(NATOMS), DA, DB
INTEGER TRIED(NATOMS), DLIST(NATOMS), SORTLIST(NATOMS), NDUMMY2, INGROUP(NATOMS), NADDED

DBEST=1.0D100
PERMUTABLE(1:NATOMS)=.FALSE.
NDUMMY=1
DO J1=1,NPERMGROUP
   DO J2=1,NPERMSIZE(J1)
      PERMUTABLE(PERMGROUP(NDUMMY+J2-1))=.TRUE.
      INGROUP(PERMGROUP(NDUMMY+J2-1))=J1
   ENDDO
   NDUMMY=NDUMMY+NPERMSIZE(J1)
ENDDO

DUMMYB(1:3*NATOMS)=COORDSB(1:3*NATOMS)
DUMMYA(1:3*NATOMS)=COORDSA(1:3*NATOMS)
!
!  Bipartite matching routine for permutations. Coordinates in DUMMYB do not change
!  but the coordinates in DUMMYA do. DISTANCE is the distance^2 in this case,
!  and is evaluated as a sum of local distances squared for permutable groups.
!  We return to label 10 after every round of permutational/orientational alignment
!  unless we have converged to the identity permutation.
!
!  The maximum number of pair exchanges associated with a group is two.
! 
DO J1=1,NATOMS
   NEWPERM(J1)=J1
ENDDO
DSUM=0.0D0
LOCALPERMNEIGH=MIN(LOCALPERMNEIGH,NATOMS)

NDUMMY=1
DO J1=1,NPERMGROUP
   PATOMS=NPERMSIZE(J1)
!  PRINT '(A,I6,A,I6)','group ',J1,' size=',NPERMSIZE(J1)
!  PRINT '(A)','members:'
!  DO J2=1,PATOMS
!     PRINT '(20I6)',NEWPERM(PERMGROUP(NDUMMY+J2-1))
!  ENDDO
   LDBEST(J1)=1.0D100
   TRIED(1:NATOMS)=0
   DO J2=1,PATOMS
      LPERMBEST(J2)=J2
   ENDDO
   XA=0.0D0; YA=0.0D0; ZA=0.0D0
   XB=0.0D0; YB=0.0D0; ZB=0.0D0
   DMEAN(1:LOCALPERMNEIGH)=1.0D100
   DO J2=1,PATOMS
!     TRIED(NEWPERM(PERMGROUP(NDUMMY+J2-1)))=-1
      TRIED(PERMGROUP(NDUMMY+J2-1))=-1
      PDUMMYA(3*(J2-1)+1)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+1)
      PDUMMYA(3*(J2-1)+2)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+2)
      PDUMMYA(3*(J2-1)+3)=DUMMYA(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+3)
!     PDUMMYB(3*(J2-1)+1)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+1)
!     PDUMMYB(3*(J2-1)+2)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+2)
!     PDUMMYB(3*(J2-1)+3)=DUMMYB(3*(NEWPERM(PERMGROUP(NDUMMY+J2-1))-1)+3)
      PDUMMYB(3*(J2-1)+1)=DUMMYB(3*(PERMGROUP(NDUMMY+J2-1)-1)+1)
      PDUMMYB(3*(J2-1)+2)=DUMMYB(3*(PERMGROUP(NDUMMY+J2-1)-1)+2)
      PDUMMYB(3*(J2-1)+3)=DUMMYB(3*(PERMGROUP(NDUMMY+J2-1)-1)+3)
      XA=XA+PDUMMYA(3*(J2-1)+1)
      YA=YA+PDUMMYA(3*(J2-1)+2)
      ZA=ZA+PDUMMYA(3*(J2-1)+3)
      XB=XB+PDUMMYB(3*(J2-1)+1)
      YB=YB+PDUMMYB(3*(J2-1)+2)
      ZB=ZB+PDUMMYB(3*(J2-1)+3)
   ENDDO
   XA=XA/PATOMS; YA=YA/PATOMS; ZA=ZA/PATOMS
   XB=XB/PATOMS; YB=YB/PATOMS; ZB=ZB/PATOMS
   SPDUMMYA(1:3*PATOMS)=PDUMMYA(1:3*PATOMS)
   SPDUMMYB(1:3*PATOMS)=PDUMMYB(1:3*PATOMS)
!
! TRIED(J2) is 0 if atom J2 is eligible to be a neighbour, but has not
! yet been tried. It is -1 if it is ineligible, or has been tried and
! broke the alignment. It is +1 if it has been tried and did not break
! the alignment. It is -1 for atoms already in the set of permutable
! atoms in question. We add neighbours one at a time in order of 
! increasing distance from primary permutable set
! and test whether they break the alignment.
!
   DMEAN(1:NATOMS)=1.0D10
   NDMEAN=0
!
! Make a sorted list of distance from the permuting atoms.
! DMEAN, SORTLIST, TRIED, PERMUTABLE, and DLIST entries refer to original
! atom labels. Use NEWPERM to find where they are in coordinate lists.
!
   outer1: DO J2=1,NATOMS
!
! Don't allow members of the same permutational group 
! to appear as reference neighbours.
!
      IF (TRIED(J2).EQ.-1) THEN
         XDUMMY=1.0D9
      ELSE
         DA=(XA-DUMMYA(3*(NEWPERM(J2)-1)+1))**2 &
  &        +(YA-DUMMYA(3*(NEWPERM(J2)-1)+2))**2 &
  &        +(ZA-DUMMYA(3*(NEWPERM(J2)-1)+3))**2
!        DB=(XB-DUMMYB(3*(NEWPERM(J2)-1)+1))**2 &
! &        +(YB-DUMMYB(3*(NEWPERM(J2)-1)+2))**2 &
! &        +(ZB-DUMMYB(3*(NEWPERM(J2)-1)+3))**2
         DB=(XB-DUMMYB(3*(J2-1)+1))**2 &
  &        +(YB-DUMMYB(3*(J2-1)+2))**2 &
  &        +(ZB-DUMMYB(3*(J2-1)+3))**2
         XDUMMY=(SQRT(DA)+SQRT(DB))/2.0D0
         IF (XDUMMY.GT.LOCALPERMCUT2) CYCLE outer1
      ENDIF

      NDMEAN=NDMEAN+1
      loop1: DO J3=1,NDMEAN ! J2
         IF (XDUMMY.LT.DMEAN(J3)) THEN
!
! Move the rest down.
!
            DO J4=NDMEAN,J3+1,-1   !    J2,J3+1,-1
               DMEAN(J4)=DMEAN(J4-1)
               SORTLIST(J4)=SORTLIST(J4-1)
            ENDDO
            DMEAN(J3)=XDUMMY
            SORTLIST(J3)=J2
            EXIT loop1
         ENDIF
      ENDDO loop1
   ENDDO outer1

71 CONTINUE
   PDUMMYA(1:3*PATOMS)=SPDUMMYA(1:3*PATOMS)
   PDUMMYB(1:3*PATOMS)=SPDUMMYB(1:3*PATOMS)

   LDBESTATOM=1.0D100
   NOTHER=0
   DO J2=1,NATOMS
      IF (TRIED(J2).EQ.1) THEN
         NOTHER=NOTHER+1
         DLIST(NOTHER)=J2
      ENDIF
   ENDDO
   ADDED=.FALSE.
   outer2: DO J2=1,NATOMS
      IF (DMEAN(J2).GT.LOCALPERMCUT2) THEN
!        PRINT '(A)',' lopermdist> No more atoms within cutoff'
         GOTO 91
      ENDIF
      IF (TRIED(SORTLIST(J2)).EQ.0) THEN
         ADDED=.TRUE.
         NOTHER=NOTHER+1
         IF (NOTHER+PATOMS.GT.NATOMS) THEN
            PRINT '(A,I6)', &
  & ' lopermdist> ERROR *** number of neighbours plus number of permutable atoms exceeds total for group ',J1
            STOP
         ENDIF
         DLIST(NOTHER)=SORTLIST(J2)
         EXIT outer2
      ENDIF
   ENDDO outer2

   NADDED=1
   IF (PERMUTABLE(DLIST(NOTHER))) THEN
!     IF (DEBUG) PRINT '(2(A,I6))',' lopermdist> Atom ',DLIST(NOTHER),' belongs to permutable set ', &
!  &                                INGROUP(DLIST(NOTHER))
      NDUMMY2=1
      DO J2=1,INGROUP(DLIST(NOTHER))-1
         NDUMMY2=NDUMMY2+NPERMSIZE(J2)
      ENDDO
      DO J2=1,NPERMSIZE(INGROUP(DLIST(NOTHER)))
         IF (PERMGROUP(NDUMMY2+J2-1).EQ.DLIST(NOTHER-NADDED+1)) CYCLE
         IF (TRIED(PERMGROUP(NDUMMY2+J2-1)).EQ.0) THEN
            NOTHER=NOTHER+1
            NADDED=NADDED+1
            IF (NOTHER+PATOMS.GT.NATOMS) THEN
               PRINT '(A,I6)', &
     ' lopermdist> ERROR *** number of neighbours plus number of permutable atoms exceeds total for group ',J1
               STOP
            ENDIF
            DLIST(NOTHER)=PERMGROUP(NDUMMY2+J2-1)
!           IF (DEBUG) PRINT '(A,I6)',' lopermdist> Adding partner atom ',DLIST(NOTHER)
         ELSE
            PRINT '(A,I6,A)',' lopermdist> ERROR *** Partner atom ',DLIST(NOTHER),' has already been tried'
            STOP
         ENDIF
      ENDDO
   ENDIF
   
   DO J2=1,NOTHER
      PDUMMYA(3*(PATOMS+J2-1)+1)=DUMMYA(3*(NEWPERM(DLIST(J2))-1)+1)
      PDUMMYA(3*(PATOMS+J2-1)+2)=DUMMYA(3*(NEWPERM(DLIST(J2))-1)+2)
      PDUMMYA(3*(PATOMS+J2-1)+3)=DUMMYA(3*(NEWPERM(DLIST(J2))-1)+3)
!     PDUMMYB(3*(PATOMS+J2-1)+1)=DUMMYB(3*(NEWPERM(DLIST(J2))-1)+1)
!     PDUMMYB(3*(PATOMS+J2-1)+2)=DUMMYB(3*(NEWPERM(DLIST(J2))-1)+2)
!     PDUMMYB(3*(PATOMS+J2-1)+3)=DUMMYB(3*(NEWPERM(DLIST(J2))-1)+3)
      PDUMMYB(3*(PATOMS+J2-1)+1)=DUMMYB(3*(DLIST(J2)-1)+1)
      PDUMMYB(3*(PATOMS+J2-1)+2)=DUMMYB(3*(DLIST(J2)-1)+2)
      PDUMMYB(3*(PATOMS+J2-1)+3)=DUMMYB(3*(DLIST(J2)-1)+3)
   ENDDO
!
! Save PDUMMYA and PDUMMYB for cycling over possible orbits in MYORIENT alignment.
!
   SPDUMMYA(3*PATOMS+1:3*(PATOMS+NOTHER))=PDUMMYA(3*PATOMS+1:3*(PATOMS+NOTHER))
   SPDUMMYB(3*PATOMS+1:3*(PATOMS+NOTHER))=PDUMMYB(3*PATOMS+1:3*(PATOMS+NOTHER))
   NCHOOSEB1=0
66 NCHOOSEB1=NCHOOSEB1+1
   NCHOOSEB2=0
31 NCHOOSEB2=NCHOOSEB2+1
   NCHOOSE1=0
65 NCHOOSE1=NCHOOSE1+1
   NCHOOSE2=0
30 NCHOOSE2=NCHOOSE2+1
!
! Reset the coordinates of the PATOMS+NOTHER atoms in PDUMMYA and PDUMMYB
! to the subset of atoms from COORDSA and COORDSB.
!
   PDUMMYA(1:3*(PATOMS+NOTHER))=SPDUMMYA(1:3*(PATOMS+NOTHER))
   PDUMMYB(1:3*(PATOMS+NOTHER))=SPDUMMYB(1:3*(PATOMS+NOTHER))

   CALL MYORIENT(PDUMMYA,DUMMY,NORBIT1,NCHOOSE1,NORBIT2,NCHOOSE2,PATOMS+NOTHER,DEBUG,ROTA,ROTINVA,PULLT)
   PDUMMYA(1:3*(PATOMS+NOTHER))=DUMMY(1:3*(PATOMS+NOTHER))
   CALL MYORIENT(PDUMMYB,DUMMY,NORBITB1,NCHOOSEB1,NORBITB2,NCHOOSEB2,PATOMS+NOTHER,DEBUG,ROTB,ROTINVB,PULLT)
   PDUMMYB(1:3*(PATOMS+NOTHER))=DUMMY(1:3*(PATOMS+NOTHER))
!
! Optimimise permutational isomer for the standard orientation for the
! current choice of atoms from the possible orbits.
!
! MINPERM does not change PDUMMYB and PDUMMYA.
!
! Note that LDISTANCE is actually the distance squared. LDBEST also has dimensions of
! length squared.
!
   LDISTANCE=0.0D0
   CALL MINPERM(PATOMS+NOTHER, PDUMMYB, PDUMMYA, BOXLX, BOXLY, BOXLZ, BULKT, LPERM, LDISTANCE, DIST2, WORSTRAD) 
!  PRINT '(A,I6,A,I6)','for group ',J1,' size=',NPERMSIZE(J1)
!  PRINT '(A)','original and new atom labels:'
!  DO J2=1,PATOMS
!     PRINT '(2I6)',PERMGROUP(NDUMMY+J2-1),NEWPERM(PERMGROUP(NDUMMY+J2-1))
!  ENDDO

!  PRINT '(I6,A)',NOTHER,' other atoms:'
!  PRINT '(A)','original and new other atom labels:'
!  DO J2=1,NOTHER
!     PRINT '(2I6)',DLIST(J2),NEWPERM(DLIST(J2))
!  ENDDO

!  PRINT '(A,3I6,3G20.10)','J1,PATOMS,NOTHER,LDBEST(J1),LDISTANCE=',J1,PATOMS,NOTHER,LDBEST(J1),LDISTANCE
!  PRINT '(A,20I6)','LPERM after MINPERM: ',LPERM(1:PATOMS+NOTHER)

   LDISTANCE=LDISTANCE
   DO J2=1,PATOMS
      IF (LPERM(J2).GT.PATOMS) THEN
         LDISTANCE=1.0D300
!        IF (DEBUG) PRINT '(A,I6,A,I6,A)',' lopermdist> For group ',J1,' with ',NOTHER,' neighbours - neighbours mix in' 
!        IF (DEBUG) PRINT '(A,I6,A,I6)',' lopermdist> atom ',J2,' lperm value is ',LPERM(J2)
         EXIT
      ENDIF
   ENDDO

   DO J2=1,NOTHER
      IF (LPERM(PATOMS+J2).NE.PATOMS+J2) THEN
!        IF (DEBUG) PRINT '(A,I6,A,I6)',' lopermdist> Atom ',DLIST(J2),' also needs to permute to ',LPERM(PATOMS+J2)
         IF (PERMUTABLE(DLIST(J2))) THEN
!           IF (DEBUG) PRINT '(2(A,I6))',' lopermdist> Atom ',DLIST(J2),' belongs to permutable set ', &
!  &                                INGROUP(DLIST(J2))
         ELSE
!           IF (DEBUG) PRINT '(2(A,I6))',' lopermdist> Atom ',DLIST(J2),' is NOT permutable!'
            LDISTANCE=1.0D300
         ENDIF
      ENDIF
   ENDDO
!
! Save the best permutation and local distance for this subset of atoms.
! NEWPERM and coordinates are only reset after all the cycles over orbits and NEWMINDIST.
! Hence we need to track a cumulative permutation and save the best current values.
!
   IF (LDISTANCE.LT.LDBESTATOM) THEN
      LDBESTATOM=LDISTANCE
      LPERMBESTATOM(1:PATOMS)=LPERM(1:PATOMS)
   ENDIF
!  PRINT '(A,2G20.10)','LDISTANCE,LDBESTATOM=',LDISTANCE,LDBESTATOM

!  PRINT '(A,4I6,2G20.10)','NORBIT1,NORBIT2,NCHOOSE1,NCHOOSE2,LDISTANCE,LDBEST=', &
! &                         NORBIT1,NORBIT2,NCHOOSE1,NCHOOSE2,LDISTANCE,LDBEST(J1)

   IF (NCHOOSE2.LT.NORBIT2) GOTO 30
   IF (NCHOOSE1.LT.NORBIT1) GOTO 65
   IF (NCHOOSEB2.LT.NORBITB2) GOTO 31
   IF (NCHOOSEB1.LT.NORBITB1) GOTO 66

!  PRINT '(A,2G20.10)','LDBESTATOM,LOCALPERMCUT=',LDBESTATOM,LOCALPERMCUT
   IF (SQRT(LDBESTATOM).GT.LOCALPERMCUT) THEN
!     IF (DEBUG) THEN
!        PRINT '(A,G15.5,A,I6)',' lopermdist> Best distance ',SQRT(LDBESTATOM), &
! &                                     ' is too large for atom ',DLIST(NOTHER)
!     ENDIF
      TRIED(DLIST(NOTHER))=-1
      IF (NADDED.GT.1) THEN
!        IF (DEBUG) THEN
!           PRINT '(A)',' lopermdist> and partner atoms:'
!           PRINT '(20I5)',DLIST(NOTHER-NADDED+1:NOTHER-1)
!        ENDIF
         TRIED(DLIST(NOTHER-NADDED+1:NOTHER-1))=-1
      ENDIF
      GOTO 71
   ELSE
!     IF (DEBUG) PRINT '(A,G20.10,3(A,I6))',' lopermdist> Best distance ',SQRT(LDBESTATOM), &
! &                    ' is OK for myorient with atom ',DLIST(NOTHER),' and ',NOTHER,' neighbours' 
      TRIED(DLIST(NOTHER))=1
      IF (NADDED.GT.1) THEN
!        IF (DEBUG) THEN
!           PRINT '(A)',' lopermdist> and partner atoms:'
!           PRINT '(20I5)',DLIST(NOTHER-NADDED+1:NOTHER-1)
!        ENDIF
         TRIED(DLIST(NOTHER-NADDED+1:NOTHER-1))=1
      ENDIF
      LDBEST(J1)=LDBESTATOM
      LPERMBEST(1:PATOMS)=LPERMBESTATOM(1:PATOMS)
!     PRINT '(A,2G20.10)','Updating permutation: sqrt(LDBEST)=',SQRT(LDBEST(J1))
!     PRINT '(A,10I6)','LPERMBEST: ',LPERMBEST(1:PATOMS)
   ENDIF
!
! Add the next eligible atom and try alignment again.
! Stop if we already have LOCALPERMNEIGH neighbours.
!
   IF (NOTHER.LT.LOCALPERMNEIGH) GOTO 71

91 CONTINUE ! jump here when there are no atoms left to try.

!  IF (DEBUG) PRINT '(2(A,I6),A,G15.5)',' lopermdist> For group ',J1,' maximum neighbours=', &
! &                                      NOTHER,' distance=',SQRT(LDBEST(J1))
!
! We now have the best permutation for group J1 and standard orientations
! based upon all atoms belonging to the two possible orbits that appear
! for the standard alignment.
!
   LPERM(1:PATOMS)=LPERMBEST(1:PATOMS)
!
! Fill SAVEPERM with NEWPERM, which contains the current best permutation
! after the previous pass through J1
!
   SAVEPERM(1:NATOMS)=NEWPERM(1:NATOMS)
!
! Update best permutation for atoms in subset J1, specified by PERMGROUP
! with offset NDUMMY (updated below after each pass through J1)
!
   DO J2=1,PATOMS
      SAVEPERM(PERMGROUP(NDUMMY+J2-1))=NEWPERM(PERMGROUP(NDUMMY+LPERMBEST(J2)-1))
!     PRINT '(2(A,I6))',' lopermdist> Atom ',NEWPERM(PERMGROUP(NDUMMY+LPERMBEST(J2)-1)), &
! &                     ' moves to position ',PERMGROUP(NDUMMY+LPERMBEST(J2)-1)
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
!
! Save current optimal permutation in NEWPERM
!
   NEWPERM(1:NATOMS)=SAVEPERM(1:NATOMS)
   DSUM=DSUM+SQRT(LDBEST(J1))
!  PRINT '(A,I6,2(A,F20.10))',' lopermdist> For group ',J1,' best distance=',SQRT(LDBEST(J1)),' total=',DSUM
!  PRINT '(A)','best permutation is now'
!  PRINT '(20I6)',NEWPERM(1:NATOMS)

!
! Update NDUMMY, the cumulative offset for PERMGROUP
!
   NDUMMY=NDUMMY+NPERMSIZE(J1)
ENDDO  !  end of loop over groups of permutable atoms
!
! NEWPERM(J1) is the atom that moves to position J1 to map COORDSA
! to the current best alignment. 
! This loop just appears to set SAVEPERM and ALLPERM equal to the current
! NEWPERM.
!
!
! Putting the ALLPERM(J1)=J1 into the second loop causes pgf90 to miscompile!!
!
DO J1=1,NATOMS
   ALLPERM(J1)=J1
ENDDO
DO J1=1,NATOMS
   SAVEPERM(J1)=ALLPERM(NEWPERM(J1))
ENDDO
ALLPERM(1:NATOMS)=SAVEPERM(1:NATOMS)
!
! At this point DUMMYA should not have changed from COORDSA, so we are
! putting COORDSA in DUMMY
!
DUMMY(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
NPERM=0
!
! Update coordinates in DUMMYA to current best overall permutation using NEWPERM.
! We are doing this to operate with NEWPERMDIST in the next block.
!
DO J3=1,NATOMS
   DUMMYA(3*(J3-1)+1)=DUMMY(3*(NEWPERM(J3)-1)+1)
   DUMMYA(3*(J3-1)+2)=DUMMY(3*(NEWPERM(J3)-1)+2)
   DUMMYA(3*(J3-1)+3)=DUMMY(3*(NEWPERM(J3)-1)+3)

!  IF (DEBUG) WRITE(*,'(A,I5,A,I5)') ' lopermdist> Overall permutations after MYORIENT alignment:'
   IF (J3.NE.NEWPERM(J3)) THEN
!     IF (DEBUG) WRITE(*,'(A,I5,A,I5)') ' lopermdist> Moving position ',NEWPERM(J3),' to ',J3
      NPERM=NPERM+1
   ENDIF
ENDDO

DISTANCE=DSUM
IF (DEBUG) WRITE(*,'(A,G20.10)') ' lopermdist> After myorient block sum of distances=',DISTANCE
!
! Save current best overall distance, permuted version of COORDSA, and permutation.
!
DBEST=DISTANCE
XBEST(1:3*NATOMS)=DUMMYA(1:3*NATOMS)
BESTPERM(1:NATOMS)=ALLPERM(1:NATOMS)
!
! At this point NEWPERM, ALLPERM, SAVEPERM, BESTPERM
! are all the same!
!
! PRINT '(A)',' lopermdist> NEWPERM, ALLPERM, SAVEPERM, BESTPERM:'
! PRINT '(4I6)',(NEWPERM(J1),ALLPERM(J1),SAVEPERM(J1),BESTPERM(J1),J1=1,NATOMS)

!!!!!!!!!!!!!!!!!!!!!!! DEBUG
!
! Test distance for COORDSA with permutation applied in BESTPERM
!
!  DO J1=1,NATOMS
!     DUMMYA(3*(J1-1)+1)=COORDSA(3*(BESTPERM(J1)-1)+1)
!     DUMMYA(3*(J1-1)+2)=COORDSA(3*(BESTPERM(J1)-1)+2)
!     DUMMYA(3*(J1-1)+3)=COORDSA(3*(BESTPERM(J1)-1)+3)
!  ENDDO

!  CALL NEWMINDIST(COORDSB,DUMMYA,NATOMS,DISTANCE,BULKT,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT)
!  CALL NEWMINDIST(COORDSB,XBEST,NATOMS,XDUMMY,BULKT,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT)
!  IF (DEBUG) WRITE(*,'(A,2G20.10)') & 
! &   ' lopermdist> distance check for permuted COORDSA and original COORDSB=',XDUMMY,DISTANCE
!!!!!!!!!!!!!!!!!!!!!!! DEBUG
!
! Now align and reorient the permuted coordinates in COORDSA 
! Try using the best locally aligned group of atoms
!
CALL NEWMINDIST(DUMMYB,XBEST,NATOMS,DISTANCE,BULKT,TWOD,'AX    ',.FALSE.,RIGID,DEBUG,RMAT)
IF (DEBUG) PRINT '(A,G20.10)',' lopermdist> after overall alignment distance=',DISTANCE
RMATBEST(1:3,1:3)=RMAT(1:3,1:3)

COORDSA(1:3*NATOMS)=XBEST(1:3*NATOMS) ! finally, best COORDSA should include permutations for DNEB input!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! IF (DEBUG) PRINT '(A)',' lopermdist> Overall permutation for COORDSA (second argument):'
! IF (DEBUG) PRINT '(20I6)',BESTPERM(1:NATOMS)

RETURN
END SUBROUTINE LOPERMDIST
