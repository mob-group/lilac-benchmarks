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
SUBROUTINE MINPERMRB(COORDSB,COORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGID,RMATBEST)
USE KEY,ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, SETS, STOCKT, GEOMDIFFTOL, &
  &            NFREEZE, RBAAT, ANGLEAXIS2, LPDGEOMDIFFTOL, BESTPERM, RBCUTOFF, NRBGROUP, RBGROUP, RBNINGROUP, &
  &            RBATOMSMAX
USE MODCHARMM,ONLY : CHRMMT
USE INTCOMMONS, ONLY : INTMINPERMT, INTINTERPT, DESMINT
USE COMMONS, ONLY : ZSYM
USE INTCUTILS, ONLY : INTMINPERM
IMPLICIT NONE

INTEGER NATOMS, NPERM, PATOMS, NCOUNT, NRUNNING1, NRUNNING2
INTEGER J3, PERM(NATOMS),NDUMMY, LPERM(NATOMS), J1, J2, NRUNNING
DOUBLE PRECISION DIST2, COORDSA(3*NATOMS), COORDSB(3*NATOMS), DISTANCE, DUMMYA(3*NATOMS), DUMMYB(3*NATOMS), DUMMY(3*NATOMS)
DOUBLE PRECISION BOXLX,BOXLY,BOXLZ,WORSTRAD,RMAT(3,3),ENERGY, VNEW(3*NATOMS), DX, DY, DZ, RMS, DBEST, XBEST(3*NATOMS)
DOUBLE PRECISION CMXA, CMXB, CMXC
DOUBLE PRECISION ROTA(3,3), ROTINVA(3,3), ROTB(3,3), ROTINVB(3,3), ROTINVBBEST(3,3), ROTABEST(3,3), RMATBEST(3,3), TMAT(3,3)
DOUBLE PRECISION CMAX, CMAY, CMAZ, CMBX, CMBY, CMBZ, RMATCUMUL(3,3)
LOGICAL DEBUG, TWOD, RIGID, BULKT, DONEGROUP(NATOMS), RBASSIGNED(NATOMS), TRYMERGE
DOUBLE PRECISION PDUMMYA(3*NATOMS), PDUMMYB(3*NATOMS), LDISTANCE, DUMMYC(3*NATOMS), XDUMMY, DISTA(NATOMS), DISTB(NATOMS), LDUMMY
DOUBLE PRECISION SAVEDIST(NATOMS), PDUMMYASAVE(3*NATOMS)

INTEGER SAVENPERMGROUP, SAVENPERMSIZE(NATOMS), SAVEPERMGROUP(NATOMS), SAVENSETS(NATOMS)
INTEGER SAVENATOMS, ATOMMAP(NATOMS), INVMAP(NATOMS), NNDUMMY, NNNDUMMY, J4, J5, NEXTRA, NSAVE, NPERMTOTAL, TEMPRBNINGROUP(NATOMS)
INTEGER, ALLOCATABLE :: RBTEMP(:), TEMPRBGROUP(:)

RBASSIGNED(1:NATOMS)=.FALSE.
NRBGROUP=0
SAVENPERMGROUP=NPERMGROUP
SAVENPERMSIZE(1:NATOMS)=NPERMSIZE(1:NATOMS)
SAVEPERMGROUP(1:NATOMS)=PERMGROUP(1:NATOMS)
SAVENSETS(1:NATOMS)=NSETS(1:NATOMS)
SAVENATOMS=NATOMS
DONEGROUP(1:NATOMS)=.FALSE.

! OPEN(UNIT=876,FILE='coordsa.xyz',STATUS='UNKNOWN')
! WRITE(876,'(i4)') SAVENATOMS
! WRITE(876,'(A)') 'coordsa initially'
! WRITE(876,'(A1,1X,3F20.10)') (ZSYM((j1+2)/3),COORDSA(J1),COORDSA(J1+1),COORDSA(J1+2),J1=1,3*NATOMS,3)
! OPEN(UNIT=774,FILE='coordsb.xyz',STATUS='UNKNOWN')
! WRITE(774,'(i4)') SAVENATOMS
! WRITE(774,'(A)') 'coordsb initially'
! WRITE(774,'(A1,1X,3F20.10)') (ZSYM((j1+2)/3),COORDSB(J1),COORDSB(J1+1),COORDSB(J1+2),J1=1,3*NATOMS,3)

! CALL POTENTIAL(COORDSA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
! PRINT '(A,2F20.10)',' minpermrb> at the start for min A energy and RMS are: ',ENERGY, RMS
! CALL POTENTIAL(COORDSB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
! PRINT '(A,2F20.10)',' minpermrb> at the start for min B energy and RMS are: ',ENERGY, RMS

!
! Trying to merge groups with common atoms may fix some problem interpolations.
! Set TRYMERGE true to use this.
!
TRYMERGE=.TRUE.
NDUMMY=1
NRUNNING=0
mainloop: DO J1=1,SAVENPERMGROUP
   IF (DONEGROUP(J1)) THEN
      NDUMMY=NDUMMY+SAVENPERMSIZE(J1)
      CYCLE mainloop
   ENDIF
!
! Permutable atoms are only allowed to be part of one rigid body. Other atoms can belong to
! more than one. Hence we fix on the first optimised permutation, but allow other atoms to
! overlap.
!
   DO J2=1,SAVENPERMSIZE(J1)
      IF (RBASSIGNED(SAVEPERMGROUP(NDUMMY+J2-1))) THEN
         DONEGROUP(J1)=.TRUE.
         IF (DEBUG) PRINT '(A,I6,A)',' minpermrb> Atom ',SAVEPERMGROUP(NDUMMY+J2-1),' is already in an RB group'
         NDUMMY=NDUMMY+SAVENPERMSIZE(J1)
         CYCLE mainloop
      ENDIF
   ENDDO
   DONEGROUP(J1)=.TRUE.
   NPERMGROUP=1
   PATOMS=SAVENPERMSIZE(J1)
   NPERMSIZE(1)=SAVENPERMSIZE(J1)
   NSETS(1)=SAVENSETS(J1)
   INVMAP(1:SAVENATOMS)=-1

   DO J2=1,SAVENPERMSIZE(J1)
      PERMGROUP(J2)=J2
      ATOMMAP(J2)=SAVEPERMGROUP(NDUMMY+J2-1)
      INVMAP(SAVEPERMGROUP(NDUMMY+J2-1))=J2
      PDUMMYA(3*(J2-1)+1)=COORDSA(3*(SAVEPERMGROUP(NDUMMY+J2-1)-1)+1)
      PDUMMYA(3*(J2-1)+2)=COORDSA(3*(SAVEPERMGROUP(NDUMMY+J2-1)-1)+2)
      PDUMMYA(3*(J2-1)+3)=COORDSA(3*(SAVEPERMGROUP(NDUMMY+J2-1)-1)+3)
      PDUMMYB(3*(J2-1)+1)=COORDSB(3*(SAVEPERMGROUP(NDUMMY+J2-1)-1)+1)
      PDUMMYB(3*(J2-1)+2)=COORDSB(3*(SAVEPERMGROUP(NDUMMY+J2-1)-1)+2)
      PDUMMYB(3*(J2-1)+3)=COORDSB(3*(SAVEPERMGROUP(NDUMMY+J2-1)-1)+3)
   ENDDO
   
   NCOUNT=SAVENPERMSIZE(J1)
   IF (SAVENSETS(J1).GT.0) THEN
      PRINT '(A)',' minpermrb> ERROR *** additional atom swaps not allowed for rigid bodies'
      STOP
   ENDIF

   NATOMS=PATOMS+2*SAVENSETS(J1)
   NSAVE=NATOMS
!
! Now check if any of the atoms in question lie in additional groups that can be permuted,
! such as the two sets of methyl hydrogens in valine. If so, we need to add this permutation
! information. This situation can only occur if NSETS(1) > 0.
!
   IF (NSETS(1).GT.0) THEN
      NNDUMMY=0
      NNNDUMMY=NPERMSIZE(1)
      outer: DO J2=1,SAVENPERMGROUP
         IF (J2.EQ.J1) THEN
            NNDUMMY=NNDUMMY+SAVENPERMSIZE(J2)
            CYCLE outer
         ENDIF
         DO J3=1,SAVENPERMSIZE(J2)
            DO J4=1,NATOMS
               IF (ATOMMAP(J4).EQ.SAVEPERMGROUP(NNDUMMY+J3)) THEN
                  NPERMGROUP=NPERMGROUP+1
                  NPERMSIZE(NPERMGROUP)=SAVENPERMSIZE(J2)
!
!  There can;t be any swaps associated with this group.
!
                  IF (SAVENSETS(J2).GT.0) THEN
                     IF (DEBUG) PRINT '(A,I6,A,I6)','minpermrb> ERROR *** number of sets associated with group ', & 
  &                                                  J2,' is ',SAVENSETS(J2)
                     STOP
                  ENDIF
                  NSETS(NPERMGROUP)=0
                  DO J5=1,SAVENPERMSIZE(J2)
                     PERMGROUP(NNNDUMMY+J5)=INVMAP(SAVEPERMGROUP(NNDUMMY+J5))
                  ENDDO
                  NNNDUMMY=NNNDUMMY+SAVENPERMSIZE(J2)
                  NNDUMMY=NNDUMMY+SAVENPERMSIZE(J2)
                  DONEGROUP(J2)=.TRUE.
                  CYCLE outer
               ENDIF
            ENDDO
         ENDDO
         NNDUMMY=NNDUMMY+SAVENPERMSIZE(J2)
      ENDDO outer
   ENDIF

   IF (DEBUG) THEN
      PRINT '(A,I6,A,I6)',' minpermrb> Number of sets of permutable atoms for group ',J1, ' is ',NPERMGROUP
      NNDUMMY=1
      DO J2=1,NPERMGROUP
         IF (DEBUG) THEN
            DO J3=1,NPERMSIZE(J2)
               WRITE(*,'(I6)',ADVANCE='NO') ATOMMAP(PERMGROUP(NNDUMMY-1+J3))
            ENDDO
            PRINT '(A)',' '
         ENDIF
         NNDUMMY=NNDUMMY+NPERMSIZE(J2)
      ENDDO
   ENDIF
!
! Add any atoms that are outside the current permutable groups but are equidistant from the 
! atoms that can be permuted.
! Forbid members of other permutational groups to avoid a subset of such atoms being included,
! which makes the bookkeeping tricky.
!
   NEXTRA=0
   SAVEDIST(1:SAVENATOMS)=-1.0D0
   outer2: DO J2=1,SAVENATOMS
      NPERMTOTAL=0
      DO J4=1,SAVENPERMGROUP
         DO J5=1,SAVENPERMSIZE(J4)
            IF (SAVEPERMGROUP(NPERMTOTAL+J5).EQ.J2) CYCLE outer2
         ENDDO
         NPERMTOTAL=NPERMTOTAL+SAVENPERMSIZE(J4)
      ENDDO
      DO J4=1,NPERMSIZE(1)
         DISTA(J4)=SQRT( &
  &                (COORDSA(3*(J2-1)+1)-PDUMMYA(3*(J4-1)+1))**2 &
  &               +(COORDSA(3*(J2-1)+2)-PDUMMYA(3*(J4-1)+2))**2 &
  &               +(COORDSA(3*(J2-1)+3)-PDUMMYA(3*(J4-1)+3))**2 )
         IF (DISTA(J4).GT.RBCUTOFF) CYCLE outer2
         DISTB(J4)=SQRT( &
  &                (COORDSB(3*(J2-1)+1)-PDUMMYB(3*(J4-1)+1))**2 &
  &               +(COORDSB(3*(J2-1)+2)-PDUMMYB(3*(J4-1)+2))**2 &
  &               +(COORDSB(3*(J2-1)+3)-PDUMMYB(3*(J4-1)+3))**2 )
         IF (DISTB(J4).GT.RBCUTOFF) CYCLE outer2
      ENDDO
!
! Sort the distances for comparison. They do not have to be equidistant from the
! permutable atoms in the two end points. This allows rotation of a whole arginine.
! To reduce the allowed RB size just reduce the cutoff radius RBCUTOFF.
!
      CALL RBSORT(NPERMSIZE(1),DISTA)
      CALL RBSORT(NPERMSIZE(1),DISTB)
      DO J4=1,NPERMSIZE(1)
         IF (ABS(DISTA(J4)-DISTB(J4)).GT.LPDGEOMDIFFTOL) CYCLE outer2
      ENDDO
      
      XDUMMY=0.5D0*(SUM(DISTA(1:NPERMSIZE(1)))+SUM(DISTB(1:NPERMSIZE(1))))/NPERMSIZE(1)
!     IF (DEBUG) PRINT '(A,I6,A,F15.5)',' minpermrb> atom ',J2,' satisfies equidistance condition with mean distance ',XDUMMY
!     IF (DEBUG) PRINT '(6F15.5)',DISTA(1:NPERMSIZE(1)),DISTB(1:NPERMSIZE(1))
!
!  Make an ordered list with furthest last so that we can remove the furthest atoms
!  one at a time of the alignment isn;t good enough.
!
write(*,*) "savedist 0" !sn402
      ! sn402: seems to be a (compiler-dependent?) bug here: if NEXTRA = 0 and the compiler
      ! doesn't short-circuit the following .or. statement, we get an error.
      ! So this line doesn't work for ifort in debug mode, at the very least.
      IF ((NEXTRA.EQ.0).OR.(XDUMMY.GT.SAVEDIST(NEXTRA))) THEN
         NATOMS=NATOMS+1
         ATOMMAP(NATOMS)=J2
         PDUMMYA(3*(NATOMS-1)+1:3*(NATOMS-1)+3)=COORDSA(3*(J2-1)+1:3*(J2-1)+3)
         PDUMMYB(3*(NATOMS-1)+1:3*(NATOMS-1)+3)=COORDSB(3*(J2-1)+1:3*(J2-1)+3)
         NEXTRA=NEXTRA+1
         write(*,*) "savedist 1"!sn402
         SAVEDIST(NEXTRA)=XDUMMY
      ELSE
         reorder: DO J4=NSAVE+1,NSAVE+NEXTRA
            IF (XDUMMY.LT.SAVEDIST(J4-NSAVE)) THEN
               DO J5=NATOMS+1,J4+1,-1 ! move the other entries up one
                  ATOMMAP(J5)=ATOMMAP(J5-1)
                  PDUMMYA(3*(J5-1)+1:3*(J5-1)+3)=PDUMMYA(3*(J5-2)+1:3*(J5-2)+3)
                  PDUMMYB(3*(J5-1)+1:3*(J5-1)+3)=PDUMMYB(3*(J5-2)+1:3*(J5-2)+3)
                    write(*,*) "savedist 2"!sn402
                  SAVEDIST(J5-NSAVE)=SAVEDIST(J5-NSAVE-1)
               ENDDO
               ATOMMAP(J4)=J2
               PDUMMYA(3*(J4-1)+1:3*(J4-1)+3)=COORDSA(3*(J2-1)+1:3*(J2-1)+3)
               PDUMMYB(3*(J4-1)+1:3*(J4-1)+3)=COORDSB(3*(J2-1)+1:3*(J2-1)+3)
               write(*,*) "savedist 3"!sn402
               SAVEDIST(J4-NSAVE)=XDUMMY
               EXIT reorder
            ENDIF
         ENDDO reorder
         NATOMS=NATOMS+1
         NEXTRA=NEXTRA+1
      ENDIF
   ENDDO outer2

20 CONTINUE
   XDUMMY=GEOMDIFFTOL
   GEOMDIFFTOL=NPERMSIZE(1)*LPDGEOMDIFFTOL 
!  GEOMDIFFTOL=LPDGEOMDIFFTOL
   PDUMMYASAVE(1:3*NATOMS)=PDUMMYA(1:3*NATOMS)
   IF (DEBUG) THEN
      PRINT '(A,I6,A)',' minpermrb> Calling minpermdist for ',NATOMS,' atoms:'
      DO J5=1,NATOMS
         WRITE(*,'(I6)',ADVANCE='NO') ATOMMAP(J5)
      ENDDO
      WRITE(*,'(A)') ' '
      PRINT '(A,I6,A)',' minpermrb> Distances of ',NEXTRA,' extra atoms:'
      DO J5=1,NEXTRA
         WRITE(*,'(F15.5)',ADVANCE='NO') SAVEDIST(J5)
      ENDDO
      WRITE(*,'(A)') ' '
   ENDIF
   CALL MINPERMDIST(PDUMMYB,PDUMMYA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGID,RMATBEST)
   GEOMDIFFTOL=XDUMMY
!
! LPDGEOMDIFFTOL is compared to the largest difference for a single atom in the aligned structures.
!
   XDUMMY=-1.0D0
   DO J2=1,NATOMS
      LDUMMY= (PDUMMYB(3*(J2-1)+1)-PDUMMYA(3*(J2-1)+1))**2 &
  &          +(PDUMMYB(3*(J2-1)+2)-PDUMMYA(3*(J2-1)+2))**2 &
  &          +(PDUMMYB(3*(J2-1)+3)-PDUMMYA(3*(J2-1)+3))**2 
      IF (LDUMMY.GT.XDUMMY) XDUMMY=LDUMMY
   ENDDO
   LDUMMY=SQRT(LDUMMY)

   IF (DEBUG) PRINT '(A,I6,2(A,G14.4))',' minpermrb> Distance for permutable group ',J1,' is ',DISTANCE, &
  &                                  ' largest displacement ',LDUMMY
!  IF ((DISTANCE.GT.NPERMSIZE(1)*LPDGEOMDIFFTOL).AND.(NATOMS.GT.PATOMS+2*SAVENSETS(J1))) THEN
   IF ((LDUMMY.GT.LPDGEOMDIFFTOL).AND.(NATOMS.GT.PATOMS+2*SAVENSETS(J1))) THEN
      NATOMS=NATOMS-1
      IF (DEBUG) PRINT '(A)',' minpermrb> Distance is too large - reduce number of equidistant atoms by one and try again'
      PDUMMYA(1:3*NATOMS)=PDUMMYASAVE(1:3*NATOMS)
      GOTO 20
   ENDIF
   PDUMMYA(1:3*SAVENATOMS)=COORDSA(1:3*SAVENATOMS)
   NRBGROUP=NRBGROUP+1
   IF (DEBUG) PRINT '(A,I6,A,I6)',' minpermrb> ',NATOMS,' atoms assigned to rigid body group ',NRBGROUP
   RBNINGROUP(NRBGROUP)=NATOMS
   DO J2=1,NATOMS
      RBASSIGNED(ATOMMAP(J2))=.TRUE.
      IF (NRUNNING+1.GT.RBATOMSMAX) THEN
         ALLOCATE(RBTEMP(NRUNNING))
         RBTEMP(1:NRUNNING)=RBGROUP(1:NRUNNING)
         DEALLOCATE(RBGROUP)
         ALLOCATE(RBGROUP(2*RBATOMSMAX))
         RBGROUP(1:NRUNNING)=RBTEMP(1:NRUNNING)
         DEALLOCATE(RBTEMP)
         RBATOMSMAX=2*RBATOMSMAX
      ENDIF
      NRUNNING=NRUNNING+1
      RBGROUP(NRUNNING)=ATOMMAP(J2)
      PDUMMYA(3*(ATOMMAP(J2)-1)+1)=COORDSA(3*(ATOMMAP(BESTPERM(J2))-1)+1)
      PDUMMYA(3*(ATOMMAP(J2)-1)+2)=COORDSA(3*(ATOMMAP(BESTPERM(J2))-1)+2)
      PDUMMYA(3*(ATOMMAP(J2)-1)+3)=COORDSA(3*(ATOMMAP(BESTPERM(J2))-1)+3)
   ENDDO
   COORDSA(1:3*SAVENATOMS)=PDUMMYA(1:3*SAVENATOMS)

   NDUMMY=NDUMMY+SAVENPERMSIZE(J1)
ENDDO mainloop

IF (.NOT.TRYMERGE) GOTO 111
!
! Now check if we can merge any of these groups. This could happen because only one permutable 
! set is initially allowed in each rigid body. Note that there could be overalpping atoms between
! rigid body groups as well.
! Only allow merging for overlapping groups.
!
!       IF (DEBUG) THEN
!           PRINT '(3(A,I6))',' minpermrb> Original RBGROUP entries:'
!           NDUMMY=0
!           DO J3=1,NRBGROUP
!              PRINT '(A,I6)',' group ',J3 
!              PRINT '(22I6)',RBGROUP(NDUMMY+1:NDUMMY+RBNINGROUP(J3))
!              NDUMMY=NDUMMY+RBNINGROUP(J3)
!           ENDDO
!       ENDIF
30 CONTINUE
NRUNNING1=0
DO J1=1,NRBGROUP
   NRUNNING2=NRUNNING1+RBNINGROUP(J1)
   DO J2=1,RBNINGROUP(J1)
      J3=RBGROUP(J2+NRUNNING1)
      ATOMMAP(J2)=J3
      PDUMMYA(3*(J2-1)+1:3*(J2-1)+3)=COORDSA(3*(J3-1)+1:3*(J3-1)+3)
      PDUMMYB(3*(J2-1)+1:3*(J2-1)+3)=COORDSB(3*(J3-1)+1:3*(J3-1)+3)
   ENDDO

   DO J2=J1+1,NRBGROUP
      NATOMS=RBNINGROUP(J1)
!
! Check for unique atoms in second group.
!
      g2atoms: DO J5=1,RBNINGROUP(J2)
         J3=RBGROUP(J5+NRUNNING2)
         DO J4=1,RBNINGROUP(J1)
            IF (J3.EQ.RBGROUP(J4+NRUNNING1)) CYCLE g2atoms
         ENDDO
         NATOMS=NATOMS+1
         ATOMMAP(NATOMS)=J3
         PDUMMYA(3*(NATOMS-1)+1:3*(NATOMS-1)+3)=COORDSA(3*(J3-1)+1:3*(J3-1)+3)
         PDUMMYB(3*(NATOMS-1)+1:3*(NATOMS-1)+3)=COORDSB(3*(J3-1)+1:3*(J3-1)+3)
      ENDDO g2atoms
!
! Only try to merge groups that have a common atom.
!
      IF (NATOMS.EQ.RBNINGROUP(J1)+RBNINGROUP(J2)) THEN
         NRUNNING2=NRUNNING2+RBNINGROUP(J2)
         CYCLE
      ENDIF

      XDUMMY=GEOMDIFFTOL
      GEOMDIFFTOL=NPERMSIZE(1)*LPDGEOMDIFFTOL
!     IF (DEBUG) THEN
!        PRINT '(3(A,I6))',' minpermrb> Calling minpermdist for ',NATOMS,' atoms merged from groups ',J1,' and ',J2
!        DO J5=1,NATOMS
!           WRITE(*,'(I6)',ADVANCE='NO') ATOMMAP(J5)
!        ENDDO
!        WRITE(*,'(A)') ' '
!     ENDIF
      CALL MINPERMDIST(PDUMMYB,PDUMMYA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGID,RMATBEST)
      GEOMDIFFTOL=XDUMMY
!
! LPDGEOMDIFFTOL is compared to the largest difference for a single atom in the aligned structures.
!
      XDUMMY=-1.0D0
      DO J3=1,NATOMS
         LDUMMY= (PDUMMYB(3*(J3-1)+1)-PDUMMYA(3*(J3-1)+1))**2 &
  &             +(PDUMMYB(3*(J3-1)+2)-PDUMMYA(3*(J3-1)+2))**2 &
  &             +(PDUMMYB(3*(J3-1)+3)-PDUMMYA(3*(J3-1)+3))**2
         IF (LDUMMY.GT.XDUMMY) XDUMMY=LDUMMY
      ENDDO
      LDUMMY=SQRT(XDUMMY) ! bug fix 30/7/09 DJW
!     IF (DEBUG) PRINT '(2(A,I6),2(A,G14.4))',' minpermrb> Distance for merged groups ',J1,' and ',J2,' is ',DISTANCE, &
! &                                           ' largest displacement ',LDUMMY
!
! If allowed, add non-duplicate entries for J2 group onto the end of J1 group, moving the rest of entries up
! appropriately. Reduced NRBGROUP by one, and start all over again.
! Otherwise, increment NRUNNING2 and try merging with the next group.
!
      IF (LDUMMY.LT.LPDGEOMDIFFTOL) THEN
         ALLOCATE(TEMPRBGROUP(RBATOMSMAX))
         NCOUNT=0
         DO J3=1,J1-1
            TEMPRBNINGROUP(J3)=RBNINGROUP(J3)
            DO J4=1,RBNINGROUP(J3)
               TEMPRBGROUP(NCOUNT+J4)=RBGROUP(NCOUNT+J4)
            ENDDO
            NCOUNT=NCOUNT+RBNINGROUP(J3)
         ENDDO
         TEMPRBGROUP(NCOUNT+1:NCOUNT+NATOMS)=ATOMMAP(1:NATOMS)
         NDUMMY=NCOUNT+RBNINGROUP(J1)
         NCOUNT=NCOUNT+NATOMS
         RBNINGROUP(J1)=NATOMS
         TEMPRBNINGROUP(J1)=NATOMS
         j3loop: DO J3=J1+1,NRBGROUP
            IF (J3.EQ.J2) THEN
               NDUMMY=NDUMMY+RBNINGROUP(J3)
               CYCLE j3loop
            ENDIF
            IF (J3.LT.J2) TEMPRBNINGROUP(J3)=RBNINGROUP(J3)
            IF (J3.GT.J2) TEMPRBNINGROUP(J3-1)=RBNINGROUP(J3)
            DO J4=1,RBNINGROUP(J3)
               TEMPRBGROUP(NCOUNT+J4)=RBGROUP(NDUMMY+J4)
            ENDDO
            NCOUNT=NCOUNT+RBNINGROUP(J3)
            NDUMMY=NDUMMY+RBNINGROUP(J3)
         ENDDO j3loop
         NRBGROUP=NRBGROUP-1
         RBGROUP(1:NCOUNT)=TEMPRBGROUP(1:NCOUNT)
         RBNINGROUP(1:NRBGROUP)=TEMPRBNINGROUP(1:NRBGROUP)
         DEALLOCATE(TEMPRBGROUP)
         IF (DEBUG) THEN
            PRINT '(3(A,I6))',' minpermrb> After merging groups ',J1,' and ',J2,' remaining groups=',NRBGROUP
            PRINT '(3(A,I6))',' minpermrb> New RBGROUP entries:'
            NDUMMY=0
            DO J3=1,NRBGROUP
               PRINT '(A,I6)',' group ',J3 
               PRINT '(22I6)',RBGROUP(NDUMMY+1:NDUMMY+RBNINGROUP(J3))
               NDUMMY=NDUMMY+RBNINGROUP(J3)
            ENDDO
         ENDIF
         GOTO 30
      ELSE
         NRUNNING2=NRUNNING2+RBNINGROUP(J2)
      ENDIF
   ENDDO
   NRUNNING1=NRUNNING1+RBNINGROUP(J1)
ENDDO

111 CONTINUE

NPERMGROUP=SAVENPERMGROUP
NATOMS=SAVENATOMS
NPERMSIZE(1:NATOMS)=SAVENPERMSIZE(1:NATOMS)
PERMGROUP(1:NATOMS)=SAVEPERMGROUP(1:NATOMS)
NSETS(1:NATOMS)=SAVENSETS(1:NATOMS)

!
! Routines like checkpair are expecting the minimum distance to be returned.
! Here we fix the permutations and optimise with respect to the overall orientation using
! newmindist.
!
CALL NEWMINDIST(COORDSB,COORDSA,NATOMS,DISTANCE,BULKT,TWOD,'AX    ',.TRUE.,RIGID,DEBUG,RMAT)

! WRITE(876,'(i4)') SAVENATOMS
! WRITE(876,'(A)') 'coordsa finally'
! WRITE(876,'(A1,1X,3F20.10)') (ZSYM((J1+2)/3),COORDSA(J1),COORDSA(J1+1),COORDSA(J1+2),J1=1,3*NATOMS,3)
! WRITE(774,'(i4)') SAVENATOMS
! WRITE(774,'(A)') 'coordsb finally '
! WRITE(774,'(A1,1X,3F20.10)') (ZSYM((J1+2)/3),COORDSB(J1),COORDSB(J1+1),COORDSB(J1+2),J1=1,3*NATOMS,3)
! CLOSE(876)
! CLOSE(774)

! STOP

! CALL POTENTIAL(COORDSA,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
! PRINT '(A,2F20.10)',' minpermrb> at the end for min A energy and RMS are: ',ENERGY, RMS
! CALL POTENTIAL(COORDSB,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
! PRINT '(A,2F20.10)',' minpermrb> at the end for min B energy and RMS are: ',ENERGY, RMS

DISTANCE=SQRT(DISTANCE)

RETURN
END SUBROUTINE MINPERMRB

!
!     This subprogram performs a sort on the input data and
!     arranges it from smallest to biggest. The exchange-sort
!     algorithm is used.
!
SUBROUTINE RBSORT(N,A)
IMPLICIT NONE
INTEGER J1, L, N, J2
DOUBLE PRECISION A(*), TEMP

DO J1=1,N-1
   L=J1
   DO J2=J1+1,N
      IF (A(L).GT.A(J2)) L=J2
   ENDDO
   TEMP=A(L)
   A(L)=A(J1)
   A(J1)=TEMP
ENDDO

RETURN
END SUBROUTINE RBSORT

