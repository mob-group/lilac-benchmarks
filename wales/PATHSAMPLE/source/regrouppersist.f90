!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using discrete path sampling and perform kinetic analysis
!   Copyright (C) 1999-2009 David J. Wales
!   This file is part of PATHSAMPLE.
!
!   PATHSAMPLE is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   PATHSAMPLE is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Regroup on the basis of persistence analysis
!
!  This routine performs a self-consistent regrouping based on a persistence
!  analysis. It is like regroupfree2, 
!  with groups containing the minima identified in persitence analysis
!  not allowed to merge, unless NOMERGEAB is false, which is set by the
!  ALLOWAB keyword, as for normal regroupfree2. Instead of just two groups,
!  A and B, we have a number of groups specified by the number of minima in 
!  persistent minima identified in min.include.
!
SUBROUTINE REGROUPPERSIST(GETPAIRST,PAIRSTODO,FREEMINLIST,FREEMINPOINT,NAVAIL)
USE COMMONS
USE UTILS,ONLY : GETUNIT
IMPLICIT NONE
INTEGER J1, J2, NGROUPS, PAIRSTODO, J3, J4, NEWA, NEWB
INTEGER NEWNMINA, NEWNMINB, NCOUNT, NDUMMY, GROUPMAP(NMIN), LOCALP, LOCALM, NMINCONNECTED
INTEGER NEWNMIN, NEWNTS, NEWPLUS(NTS), NEWMINUS(NTS), NMINGROUP(NMIN), NTSGROUP(NTS), CURRENTINDEX(NMIN)
INTEGER GROUPCONN(2*NTS), GROUPTS(2*NTS), STARTGROUP(NMIN)
INTEGER FREEMINLIST(NMIN), FREEMINPOINT(0:NMIN+1), FREETSLIST(NTS), FREETSPOINT(0:NTS), NAVAIL
INTEGER NCOL(NMIN), NVAL(NCONNMAX,NMIN), NDISTA(NMIN), NDISTB(NMIN), NCYCLE, DMIN, TSMAP(NTS), NEWJ1
INTEGER DMAX, NUNCONA, NUNCONB, NDEAD, LUNIT, NPERSIST
LOGICAL DEADTS(NTS), LREJECTTS, NOMERGEAB, GETPAIRST, FIRSTPASS, CHECKCONN
LOGICAL ISA(NMIN), ISB(NMIN), CHANGED, GROUPA(NMIN), GROUPB(NMIN), TESTIT
DOUBLE PRECISION NEWEMIN(NMIN), NEWETS(NTS), NEWPFMIN(NMIN), NEWKPLUS(NTS), NEWKMINUS(NTS)
DOUBLE PRECISION KSUM(NMIN), HBARRIER(NTS)
DOUBLE PRECISION DMATMC(NCONNMAX,NMIN)
DOUBLE PRECISION LNPROD, PFTOTAL
DOUBLE PRECISION :: CUT_UNDERFLOW=-300.0D0
CHARACTER(LEN=20) TSTRING
CHARACTER(LEN=80) FNAME
INTEGER, ALLOCATABLE :: NPGROUP(:), NEWNMINP(:), NMINP(:), LOCATIONP(:,:)
LOGICAL, ALLOCATABLE :: ISP(:,:), GROUPP(:,:)

!
! NMIN may have changed since last call. Cannot deallocate nconngroup in regrouppersist
! because it may be used elsewhere.
!
WRITE(TSTRING,'(F20.10)') TEMPERATURE
IF (ALLOCATED(NCONNGROUP)) DEALLOCATE(NCONNGROUP)
ALLOCATE(NCONNGROUP(NMIN))
FIRSTPASS=.TRUE.
NOMERGEAB=.TRUE.
IF (ALLOWABT) NOMERGEAB=.FALSE.
IF (ENSEMBLE.EQ.'E') THEN
   PRINT '(A)','regrouppersist> Regrouped entropies not yet coded for microcanonical ensemble'
   STOP
ENDIF
!
!  Remove potential energy minima that are not connected to product or reactant at this
!  stage. Why bother including them as free energy groups?
!  we have to set up the NCOL and NVAL arrays with MAKED first to do this.
!  We therefore have to set KSUM and DEADTS as well.
!
CALL RATECONST_SETUP(KSUM,DEADTS,NDEAD,.TRUE.,-300.0D0)

CALL MAKED(DMATMC,NCOL,NVAL,DEADTS,.TRUE.,ISA,ISB,KSUM)
!
!  Calculate minimum number of steps of each minimum from the A set.
!  
NDISTA(1:NMIN)=1000000
DO J1=1,NMINA
   NDISTA(LOCATIONA(J1))=0
ENDDO 
NCYCLE=0
5    CHANGED=.FALSE.
NCYCLE=NCYCLE+1
DMIN=100000
DMAX=0
NUNCONA=0
DO J1=1,NMIN
   IF (NDISTA(J1).EQ.0) CYCLE ! A minimum
   DO J2=1,NCOL(J1)
      IF (NDISTA(NVAL(J2,J1))+1.LT.NDISTA(J1)) THEN
         CHANGED=.TRUE.
         NDISTA(J1)=NDISTA(NVAL(J2,J1))+1
      ENDIF
   ENDDO
   IF ((NDISTA(J1).GT.DMAX).AND.(NDISTA(J1).NE.1000000)) DMAX=NDISTA(J1)
   IF (NDISTA(J1).LT.DMIN) DMIN=NDISTA(J1)
   IF (NDISTA(J1).EQ.1000000) NUNCONA=NUNCONA+1
ENDDO 
IF (CHANGED) GOTO 5
PRINT '(3(A,I8))','regrouppersist> steps to A region converged in ',NCYCLE-1, &
&                    ' cycles; maximum=',DMAX,' disconnected=',NUNCONA
!
!  Calculate minimum number of steps of each minimum from the B set.
!
NDISTB(1:NMIN)=1000000
DO J1=1,NMINB
   NDISTB(LOCATIONB(J1))=0
ENDDO
NCYCLE=0
51    CHANGED=.FALSE.
NCYCLE=NCYCLE+1
DMIN=100000
DMAX=0
NUNCONB=0
DO J1=1,NMIN
   IF (NDISTB(J1).EQ.0) CYCLE ! B minimum
   DO J2=1,NCOL(J1)
      IF (NDISTB(NVAL(J2,J1))+1.LT.NDISTB(J1)) THEN
         CHANGED=.TRUE.
         NDISTB(J1)=NDISTB(NVAL(J2,J1))+1
      ENDIF
   ENDDO
   IF ((NDISTB(J1).GT.DMAX).AND.(NDISTB(J1).NE.1000000)) DMAX=NDISTB(J1)
   IF (NDISTB(J1).LT.DMIN) DMIN=NDISTB(J1)
   IF (NDISTB(J1).EQ.1000000) NUNCONB=NUNCONB+1
ENDDO
IF (CHANGED) GOTO 51
PRINT '(3(A,I8))','regrouppersist> steps to B region converged in ',NCYCLE-1, &
&                    ' cycles; maximum=',DMAX,' disconnected=',NUNCONB
!
!  This could happen if disconnected minima lie in the A or B region.
!
IF (NUNCONB.NE.NUNCONA) PRINT '(A)', &
&                   'regrouppersist> WARNING - number of disconnected minima from A and B is different'
!
!  Check that we actually have a connection between the A and B regions.
!  If not, STOP.
!
CHECKCONN=.FALSE.
IF (DIRECTION.EQ.'AB') THEN
   DO J1=1,NMINB
      IF (NDISTA(LOCATIONB(J1)).LT.1000000) THEN
         CHECKCONN=.TRUE.
         EXIT
      ENDIF
   ENDDO
ELSE
   DO J1=1,NMINA
      IF (NDISTB(LOCATIONA(J1)).LT.1000000) THEN
         CHECKCONN=.TRUE.
         EXIT
      ENDIF
   ENDDO
ENDIF
IF (.NOT.CHECKCONN) THEN
   PRINT '(A)','regrouppersist> There is no connection between the A and B regions'
!  OPEN(UNIT=1,FILE='ts.attempts',STATUS='UNKNOWN')
!  WRITE(1,'(I8)') TSATTEMPT(1:NTS)
!  CLOSE(1)
   STOP
ENDIF

NMINCONNECTED=0
DO J1=1,NMIN
   IF ((NDISTA(J1).EQ.1000000).OR.(NDISTB(J1).EQ.1000000)) NCONN(J1)=0 ! exclude this pe minimum in everything that follows
   IF (NCONN(J1).GT.NCONNMIN) NMINCONNECTED=NMINCONNECTED+1
ENDDO
PRINT '(A,I8)','regrouppersist> Number of minima remaining after allowing for minimum connectivity and ts threshold=',NMINCONNECTED

ISA(1:NMIN)=.FALSE.
ISB(1:NMIN)=.FALSE.
DO J1=1,NMINA
   ISA(LOCATIONA(J1))=.TRUE.
ENDDO
DO J1=1,NMINB
   ISB(LOCATIONB(J1))=.TRUE.
ENDDO

LUNIT=GETUNIT()
OPEN(UNIT=LUNIT,FILE='min.include',STATUS='OLD')
NPERSIST=0
DO 
   READ(LUNIT,*,END=864) NDUMMY
   NPERSIST=NPERSIST+1
ENDDO
864 CONTINUE
ALLOCATE(NPGROUP(NPERSIST),NEWNMINP(NPERSIST),GROUPP(NPERSIST,NMIN),ISP(NPERSIST,NMIN),NMINP(NPERSIST),LOCATIONP(NPERSIST,NMIN))
!
! Set up initial ISP array with one minimum for each persistent minimum.
!
ISP(1:NPERSIST,1:NMIN)=.FALSE.
REWIND(LUNIT)
DO J1=1,NPERSIST
   READ(LUNIT,*,END=864) LOCATIONP(J1,1)
   PRINT '(A,I6,A,I6)','regrouppersist> Persistent minimum ',J1,' number ',LOCATIONP(J1,1)
   ISP(J1,LOCATIONP(J1,1))=.TRUE.
ENDDO
CLOSE(LUNIT)

NGROUPS=0
!
!  Assign minima to new groups. Initially, each group consists of a single connected minimum,
!  which constitutes its own free energy minimum. Minima with .LE. NCONNMIN connections are ignored.
!  NEWEMIN(J1) contains the free energy of group J1
!  NEWPFMIN(J1) contains the superposition partition function for group J1
!  NEWNMIN is the number of free energy minima
!  NEWKPLUS(J1) is the plus rate for inter-group ts J1
!  NEWKMINUS(J1) is the minus rate for inter-group ts J1
!  NEWNTS is the number of inter-group transition states
!  NTSGROUP(J1) is the number of transition states in inter-group J1
!  NEWETS(J1) is the effective free energy for the inter-group transition state J1
!
!  MINGROUP(J1) is the index of the group containing minimum J1
!
IF (ALLOCATED(MINGROUP)) DEALLOCATE(MINGROUP)
ALLOCATE(MINGROUP(NMIN)) 

MINGROUP(1:NMIN)=0 

NEWEMIN(1:NMIN)=HUGE(1.0D0)
NEWPFMIN(1:NMIN)=0.0D0
NEWNMIN=0

DO J1=1,NMIN
   IF (NCONN(J1).LE.NCONNMIN) CYCLE
   NEWNMIN=NEWNMIN+1
   MINGROUP(J1)=NEWNMIN
   NEWPFMIN(NEWNMIN)=EXP(PFMIN(J1))
   NMINGROUP(NEWNMIN)=1
ENDDO
PRINT '(A,G20.10)','regrouppersist> PFMEAN=',PFMEAN
DO J1=1,NEWNMIN
   IF (NEWPFMIN(J1).EQ.0.0D0) NEWPFMIN(J1)=1.0D-10 ! to avoid underflow
!
! NEWEMIN values are shifted by the missing term TEMPERATURE*LOG(PFMEAN)
!
   NEWEMIN(J1)=-TEMPERATURE*LOG(NEWPFMIN(J1))
   IF (DEBUG) PRINT '(A,I7,A,G20.10,A,I7,A,G20.10)','regrouppersist> For initial group ',J1,' Z(T)=',NEWPFMIN(J1), &
  &                  ' size ',NMINGROUP(J1),' free energy=',NEWEMIN(J1)
ENDDO
NEWKPLUS(1:NTS)=0.0D0
NEWKMINUS(1:NTS)=0.0D0
NTSGROUP(1:NTS)=0
NEWNTS=0
DO J1=1,NTS
   CALL CHECKTS(ETS(J1),EMIN(PLUS(J1)),EMIN(MINUS(J1)),KPLUS(J1),KMINUS(J1),NCONN(PLUS(J1)),NCONN(MINUS(J1)), &
                PLUS(J1),MINUS(J1),.FALSE.,CUT_UNDERFLOW,LREJECTTS)
   IF (LREJECTTS) CYCLE
   NEWNTS=NEWNTS+1
   NEWKPLUS(NEWNTS)=EXP(KPLUS(J1))
   NEWKMINUS(NEWNTS)=EXP(KMINUS(J1))
   NTSGROUP(NEWNTS)=1
   NEWPLUS(NEWNTS)=MINGROUP(PLUS(J1))
   NEWMINUS(NEWNTS)=MINGROUP(MINUS(J1))
ENDDO 

PRINT '(A,I7)','regrouppersist> Number of intergroup transition states=',NEWNTS
DO J1=1,NEWNTS
   IF (DEBUG) PRINT '(3(A,I7),2(A,G20.10),A,I7)','regrouppersist> Grouped ts ',J1,' between minima groups ',NEWPLUS(J1), &
  &    ' and ',NEWMINUS(J1), &
  &    ' k+=',NEWKPLUS(J1),' k-=',NEWKMINUS(J1),' members=',NTSGROUP(J1)
   IF (DEBUG) PRINT '(A,2G20.10)','regrouppersist> detailed balance - these numbers should be equal: ', &
  &      NEWKPLUS(J1)*NEWPFMIN(NEWPLUS(J1)), NEWKMINUS(J1)*NEWPFMIN(NEWMINUS(J1))

   IF ((NEWKPLUS(J1).EQ.0.0D0).OR.(NEWKMINUS(J1).EQ.0.0D0)) THEN
      IF (DEBUG) PRINT '(A,I7,2G20.10)','regrouppersist> WARNING - J1,NEWKPLUS,NEWKMINUS=',J1,NEWKPLUS(J1),NEWKMINUS(J1)
! 
! Setting this value to huge can cause an effectively infinite cycle in getfreebarrier, where the
! threshold increases to the huge value.
!
      NEWETS(J1)=HUGE(1.0D0)
   ELSE
!
! NEWEMIN has been shifted by T*LOG(PFMEAN), but NEWETS is shifted by exactly the same amount.
!
      NEWETS(J1)=NEWEMIN(NEWPLUS(J1))-TEMPERATURE*(LOG(NEWKPLUS(J1))+LOG(PLANCK/TEMPERATURE))
      IF (DEBUG) PRINT '(3(A,G20.10))','regrouppersist> Grouped ts free energy=', NEWETS(J1), &
  &              ' or ',NEWEMIN(NEWMINUS(J1))-TEMPERATURE*(LOG(NEWKMINUS(J1))+LOG(PLANCK/TEMPERATURE)), &
  &              ' or ',NEWEMIN(NEWPLUS(J1))-TEMPERATURE*(LOG(NEWKPLUS(J1))+LOG(PLANCK/TEMPERATURE))
      IF (NEWETS(J1).NE.0.0D0) THEN ! Check for consistency
         IF (ABS((NEWETS(J1)-NEWEMIN(NEWMINUS(J1))+TEMPERATURE*(LOG(NEWKMINUS(J1))+ &
  &                LOG(PLANCK/TEMPERATURE)))/NEWETS(J1)).GT.0.01D0) THEN
            PRINT '(A,I7,A,3G20.10)','regrouppersist> WARNING - free energies for ts group ',J1,' are ',  &
  &                  NEWETS(J1),NEWEMIN(NEWMINUS(J1))-TEMPERATURE*(LOG(NEWKMINUS(J1))+LOG(PLANCK/TEMPERATURE)), &
  &                             NEWEMIN(NEWPLUS(J1))-TEMPERATURE*(LOG(NEWKPLUS(J1))+LOG(PLANCK/TEMPERATURE))
            STOP
         ENDIF
      ENDIF
   ENDIF
ENDDO
!
!  A set.
!
NEWNMINA=0
GROUPA(1:NMIN)=.FALSE.
DO J1=1,NMIN
   IF (NCONN(J1).LE.NCONNMIN) CYCLE
   IF (ISA(J1)) GROUPA(MINGROUP(J1))=.TRUE.
ENDDO
DO J1=1,NMIN
   IF (NCONN(J1).LE.NCONNMIN) CYCLE ! MINGROUP(J1) is zero in this case!
   IF (GROUPA(MINGROUP(J1))) THEN
      ISA(J1)=.TRUE.
      NEWNMINA=NEWNMINA+1
      IF (DEBUG) PRINT '(A,I7,A)','regroup> potential energy minimum ',J1,' is in a free energy group that contains an A minimum'
   ENDIF
ENDDO
!
!  B set.
!  
NEWNMINB=0
GROUPB(1:NMIN)=.FALSE.
DO J1=1,NMIN
   IF (NCONN(J1).LE.NCONNMIN) CYCLE
   IF (ISB(J1)) GROUPB(MINGROUP(J1))=.TRUE.
ENDDO 
DO J1=1,NMIN
   IF (NCONN(J1).LE.NCONNMIN) CYCLE ! MINGROUP(J1) is zero in this case!
   IF (GROUPB(MINGROUP(J1))) THEN
      ISB(J1)=.TRUE.
      NEWNMINB=NEWNMINB+1
      IF (DEBUG) PRINT '(A,I7,A)','regroup> potential energy minimum ',J1,' is in a free energy group that contains a  B minimum'
   ENDIF 
ENDDO
!
!  Persistent sets.
!
NEWNMINP(1:NPERSIST)=0
GROUPP(1:NPERSIST,1:NMIN)=.FALSE.
DO J2=1,NPERSIST
   DO J1=1,NMIN
      IF (NCONN(J1).LE.NCONNMIN) CYCLE   
      IF (ISP(J2,J1)) GROUPP(J2,MINGROUP(J1))=.TRUE.
   ENDDO
   DO J1=1,NMIN
      IF (NCONN(J1).LE.NCONNMIN) CYCLE ! MINGROUP(J1) is zero in this case!
      IF (GROUPP(J2,MINGROUP(J1))) THEN
         ISP(J2,J1)=.TRUE.
         NEWNMINP(J2)=NEWNMINP(J2)+1
         IF (DEBUG) PRINT '(A,I7,A,I6)','regroup> potential energy minimum ',J1, &
  &                                     ' is in a free energy group that contains a persistent minimum ',J2
      ENDIF
   ENDDO
ENDDO
!
!  Initial setup is complete. No grouping has occurred yet, but some minima and transition states 
!  have been removed through TSTHRESH, MAXBARRIER and minimum connection conditions.
!
888 CONTINUE ! Top of iterative regrouping loop !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
CHANGED=.FALSE.
DO J1=1,NEWNMIN
   CURRENTINDEX(J1)=J1 ! currentindex tracks where we move the groups to once they have merged
ENDDO

! Find the higher of the two barrier heights, in free energy, for each TS.
DO J1=1,NEWNTS
   HBARRIER(J1)=MAX( (NEWETS(J1)-NEWEMIN(NEWPLUS(J1))), (NEWETS(J1)-NEWEMIN(NEWMINUS(J1))) )
   TSMAP(J1)=J1 ! also set this index-tracking array
END DO

! Order the TSs on the value of the higher barrier height, from small to large
! (The higher barrier height is the controlling one for each TS, because the test applied below 
! is that both barriers must lie below the threshold for grouping).
CALL SORT5(NEWNTS,NEWNTS,HBARRIER(1:NEWNTS),TSMAP(1:NEWNTS))

! We complete a sweep over all transition states using the free energies from the
! previous cycle, before updating the energies for the next cycles.
! Hence it is possible for the energies printed in the new min.data and ts.data
! files to give barriers below a given threshold, when that threshold will
! actually cause further regrouping.
!
! Alternatively we could jump out of the following loop after every regrouping even and
! start again with the new groups and energies.
! Hopefully these alternatives should not produce significantly different results.
!
DO J1=1,NEWNTS
!  PRINT '(A,2I7)','J1,NEWNTS=',J1,NEWNTS
!  PRINT '(A,2I7)','minima and current energies:'
!  PRINT '(I7,G20.10)',(J2,NEWEMIN(J2),J2=1,NEWNMIN)
!  PRINT '(A,2I7)','transition states and current energies:'
!  PRINT '(I7,G20.10)',(J2,NEWETS(J2),J2=1,NEWNTS)
!
! The order in which we encounter TSs and test them as below can affect the outcome of grouping.
! Now assessing TSs in order of low to high barrier height.
!
   NEWJ1=TSMAP(J1) ! NEWJ1 is the original index of the TS before the barrier-height reordering
   LOCALP=CURRENTINDEX(NEWPLUS(NEWJ1))
   LOCALM=CURRENTINDEX(NEWMINUS(NEWJ1))
   IF (LOCALP.EQ.LOCALM) CYCLE ! ignore intra-group rates
   TESTIT=(NEWETS(NEWJ1)-NEWEMIN(NEWPLUS(NEWJ1)).LE.REGROUPFREETHRESH).AND. &
         &(NEWETS(NEWJ1)-NEWEMIN(NEWMINUS(NEWJ1)).LE.REGROUPFREETHRESH)
!
!  PRINT '(A,I10,L5)','NEWJ1,TESTIT=',NEWJ1,TESTIT
!  PRINT '(A,4G20.10)','NEWETS(NEWJ1),NEWEMIN(NEWPLUS(NEWJ1)),NEWEMIN(NEWMINUS(NEWJ1)),REGROUPFREETHRESH=', &
! &            NEWETS(NEWJ1),NEWEMIN(NEWPLUS(NEWJ1)),NEWEMIN(NEWMINUS(NEWJ1)),REGROUPFREETHRESH
   IF (NOMERGEAB) THEN
      DO J3=1,NPERSIST
         DO J2=J3+1,NPERSIST
            IF (GROUPP(J3,LOCALP).AND.GROUPP(J2,LOCALM)) TESTIT=.FALSE.
         ENDDO
      ENDDO
      IF (GROUPA(LOCALP).AND.GROUPB(LOCALM)) TESTIT=.FALSE.
      IF (GROUPA(LOCALM).AND.GROUPB(LOCALP)) TESTIT=.FALSE.
   ENDIF
   IF (TESTIT) THEN
      CHANGED=.TRUE.
      IF (DEBUG) PRINT '(4(A,I7),A,3G15.5)','regrouppersist> merging groups ',NEWPLUS(NEWJ1),' and ',NEWMINUS(NEWJ1),' now at ', &
   &         LOCALP,' and ',LOCALM,' E ts,+,- ',NEWETS(NEWJ1), &
   &         NEWEMIN(NEWPLUS(NEWJ1)),NEWEMIN(NEWMINUS(NEWJ1))
!
!  Move minima from group with higher index to group with lower index.
!
      IF (LOCALP.LT.LOCALM) THEN
         NMINGROUP(LOCALP)=NMINGROUP(LOCALP)+NMINGROUP(LOCALM)
         NMINGROUP(LOCALM)=0
         DO J2=1,NPERSIST
            IF (GROUPP(J2,LOCALM)) GROUPP(J2,LOCALP)=.TRUE.
         ENDDO
         DO J2=1,NPERSIST
            GROUPP(J2,LOCALM)=.FALSE.
         ENDDO
         DO J2=1,NMIN
            IF (MINGROUP(J2).EQ.LOCALM) MINGROUP(J2)=LOCALP
         ENDDO
         NDUMMY=CURRENTINDEX(NEWMINUS(NEWJ1))
         DO J2=1,NEWNMIN
            IF (CURRENTINDEX(J2).EQ.NDUMMY) CURRENTINDEX(J2)=LOCALP
         ENDDO
!        CURRENTINDEX(NEWMINUS(NEWJ1))=LOCALP ! Any CURRENTINDEX(J2) = CURRENTINDEX(NEWMINUS(NEWJ1))
!                                          ! should also change to LOCALP ! DJW 11/1/08
      ELSE
         NMINGROUP(LOCALM)=NMINGROUP(LOCALM)+NMINGROUP(LOCALP)
         NMINGROUP(LOCALP)=0
         DO J2=1,NPERSIST
            IF (GROUPP(J2,LOCALP)) GROUPP(J2,LOCALM)=.TRUE.
         ENDDO
         DO J2=1,NPERSIST
            GROUPP(J2,LOCALP)=.FALSE.
         ENDDO
         DO J2=1,NMIN
            IF (MINGROUP(J2).EQ.LOCALP) MINGROUP(J2)=LOCALM
         ENDDO
         NDUMMY=CURRENTINDEX(NEWPLUS(NEWJ1))
         DO J2=1,NEWNMIN
            IF (CURRENTINDEX(J2).EQ.NDUMMY) CURRENTINDEX(J2)=LOCALM
         ENDDO
!        CURRENTINDEX(NEWPLUS(NEWJ1))=LOCALM ! Any CURRENTINDEX(J2) = CURRENTINDEX(NEWPLUS(NEWJ1))
!                                         ! should also change to LOCALM ! DJW 11/1/08
      ENDIF
!
! MINGROUP changes on merger to the lower group index.
! So MINGROUP(x) changes whenever the group containing pe min x changes.
! CURRENTINDEX(original group index) tells us where the original group maps to.
! So CURRENTINDEX(original MINGROUP(pe number)) should be current MINGROUP(pe number).
   ENDIF
ENDDO
! If there was no regrouping we have finished if CHANGED is false.
! However, we must do at least one pass to set up some of the arrays for
! free energy groups, even if they are exactly the same at the potential 
! energy groups!
IF ((.NOT.CHANGED).AND.(.NOT.FIRSTPASS)) GOTO 777 
FIRSTPASS=.FALSE.
!
!  Renumber groups of free energy minima. The free energy transition states
!  for inter-group rates are done from scratch each time. The free energy of
!  the groups and superposition partition functions are also recalculated from 
!  scratch.
!
NGROUPS=0
NDUMMY=0
DO J1=1,NEWNMIN
   IF (NMINGROUP(J1).GT.0) THEN 
      NGROUPS=NGROUPS+1
      DO J2=1,NMIN
         IF (MINGROUP(J2).EQ.J1) MINGROUP(J2)=NGROUPS
      ENDDO
      NMINGROUP(NGROUPS)=NMINGROUP(J1)
      DO J2=1,NPERSIST
         GROUPP(J2,NGROUPS)=GROUPP(J2,J1)
      ENDDO
      GROUPA(NGROUPS)=GROUPA(J1)
      GROUPB(NGROUPS)=GROUPB(J1)
      NDUMMY=NDUMMY+NMINGROUP(NGROUPS)
   ENDIF
ENDDO

PRINT '(4(A,I7))','regrouppersist> Number of free energy groups is now=',NGROUPS,' total PE minima=',NDUMMY
NEWNMIN=NGROUPS
!
!  A set.
!
NEWNMINA=0
DO J1=1,NMIN
   IF (NCONN(J1).LE.NCONNMIN) CYCLE
   IF (ISA(J1).AND.(.NOT.GROUPA(MINGROUP(J1)))) THEN
      PRINT '(2(A,I8),A,L5)','regrouppersist> ERROR - A minimum ',J1,' in group ',MINGROUP(J1),' where GROUPA=',GROUPA(MINGROUP(J1))
      STOP
   ENDIF
ENDDO
DO J1=1,NMIN
   IF (NCONN(J1).LE.NCONNMIN) CYCLE
   IF (GROUPA(MINGROUP(J1))) THEN
      ISA(J1)=.TRUE.
      NEWNMINA=NEWNMINA+1
      IF (DEBUG) PRINT '(A,I7,A)','regroup> potential energy minimum ',J1,' is in a free energy group that contains an A minimum'
   ENDIF
ENDDO
!
!  B set.
!
NEWNMINB=0
DO J1=1,NMIN
   IF (NCONN(J1).LE.NCONNMIN) CYCLE
   IF (ISB(J1).AND.(.NOT.GROUPB(MINGROUP(J1)))) THEN
      PRINT '(2(A,I8),A,L5)','regrouppersist> ERROR - B minimum ',J1,' in group ',MINGROUP(J1),' where GROUPB=',GROUPB(MINGROUP(J1))
      STOP
   ENDIF
ENDDO
DO J1=1,NMIN
   IF (NCONN(J1).LE.NCONNMIN) CYCLE
   IF (GROUPB(MINGROUP(J1))) THEN
      ISB(J1)=.TRUE.
      NEWNMINB=NEWNMINB+1
      IF (DEBUG) PRINT '(A,I7,A)','regroup> potential energy minimum ',J1,' is in a free energy group that contains a  B minimum'
   ENDIF
ENDDO
PRINT '(2(A,I7))','regrouppersist> After regrouping number of A PE minima=',NEWNMINA,' number of B PE minima=',NEWNMINB
IF ((NEWNMINA.EQ.0).OR.(NEWNMINB.EQ.0)) THEN
   PRINT '(A)','regroupfree2> ERROR - one or more of the A and B sets is empty!'
   STOP
ENDIF

IF (NOMERGEAB) THEN
   DO J1=1,NMIN
      IF (ISA(J1).AND.ISB(J1)) THEN
         PRINT '(3(A,I7))','regroup> WARNING - minimum ',J1,' belongs to A and B sets, MINGROUP=',MINGROUP(J1)
      ENDIF
   ENDDO
ENDIF
!
! Persistent sets.
!
DO J2=1,NPERSIST
   NEWNMINP(J2)=0

   DO J1=1,NMIN
      IF (NCONN(J1).LE.NCONNMIN) CYCLE
      IF (ISP(J2,J1).AND.(.NOT.GROUPP(J2,MINGROUP(J1)))) THEN
         PRINT '(2(A,I8),A,L5)','regrouppersist> ERROR - persistent minimum ',J1,' in group ',MINGROUP(J1), &
  &                             ' where GROUPP=',GROUPP(J2,MINGROUP(J1))
         STOP
      ENDIF
   ENDDO

   DO J1=1,NMIN
      IF (NCONN(J1).LE.NCONNMIN) CYCLE
      IF (GROUPP(J2,MINGROUP(J1))) THEN
         ISP(J2,J1)=.TRUE.
         NEWNMINP(J2)=NEWNMINP(J2)+1
         IF (DEBUG) PRINT '(A,I7,A)','regroup> potential energy minimum ',J1, &
  &                                  ' is in a free energy group that contains persistent minimum ',J2 
      ENDIF 
   ENDDO
   PRINT '(2(A,I7))','regrouppersist> After regrouping number of PE minima in persistent minimum ',J2,' group=',NEWNMINP(J2)
   IF (NEWNMINP(J2).EQ.0) THEN
      PRINT '(A,I6)','regrouppersist> ERROR - empty group for persistent minimum ',J2
      STOP
   ENDIF

ENDDO

IF (NOMERGEAB) THEN
   DO J1=1,NMIN
      DO J2=1,NPERSIST
         DO J3=J2+1,NPERSIST
            IF (ISP(J2,J1).AND.ISP(J3,J1)) THEN
               PRINT '(A,I7,A,2I7)','regrouppersist> ERROR - minimum ',J1,' belongs to persistent sets ',J2,J3
               STOP
            ENDIF
         ENDDO
      ENDDO
   ENDDO
ENDIF
!
!  Need to reset NEWNMIN, NEWNTS, NEWEMIN, NEWETS, NEWPLUS, NEWMINUS, NEWKPLUS and NEWKMINUS
!  to the corresponding grouped quantities and free energies.
!  Note that some of the groups of minima are actually empty (NEMPTY).
!  LOCATIONA and LOCATIONB are used in GT, so we need to reset them.
!  Probably also need to redo the TOPPOINTER and POINTER stuff for the
!  regrouped database, but only once the iterative regrouping has finished.
!  We can't renumber everything until we have calculated the free energy
!  of the grouped transition states.
!
!  Only the odd factor of Planck's constant that shifts transition state
!  free energies from minima is included.
!
NEWEMIN(1:NMIN)=HUGE(1.0D0)
NEWPFMIN(1:NMIN)=0.0D0
DO J1=1,NMIN
   IF (NCONN(J1).LE.NCONNMIN) CYCLE
   NEWPFMIN(MINGROUP(J1))=NEWPFMIN(MINGROUP(J1))+EXP(PFMIN(J1))
ENDDO
PFTOTAL=0.0D0
DO J1=1,NEWNMIN
   PFTOTAL=PFTOTAL+ NEWPFMIN(J1)
ENDDO
DO J1=1,NEWNMIN
   IF (NEWPFMIN(J1).EQ.0.0D0) NEWPFMIN(J1)=1.0D-10
   NEWEMIN(J1)=-TEMPERATURE*LOG(NEWPFMIN(J1))
!  IF (DEBUG) PRINT '(A,I7,A,G20.10,A,I7,2(A,G20.10))','regrouppersist> For group ',J1,' Z(T)=',NEWPFMIN(J1),' size ', &
! & NMINGROUP(J1), &
! &                           ' free energy=',NEWEMIN(J1),' Peq=',NEWPFMIN(J1)/PFTOTAL
   IF (DEBUG) PRINT '(A,I7,A,G20.10,A,I7,2(A,G20.10))','regrouppersist> For group ',J1,' Z(T)=',NEWPFMIN(J1),' size ', &
  &                           NMINGROUP(J1), &
  &                           ' free energy=',NEWEMIN(J1),' Peq=',NEWPFMIN(J1)/PFTOTAL
ENDDO
!
!  Store the connections for each group and the corresponding inter-group transition state
!  in a linear array.
!
!  NMINGROUP is the number of PE minima in each group.
!  NCONNGROUP(J1) is the number of other groups J1 is connected to
!  STARTGROUP(J1) is the position in array GROUPCONN(:) that the connections
!                 of group J1 start from
!  GROUPCONN(STARTGROUP(J1)+N-1) is the index of the Nth group connected to group J1
!  GROUPTS(STARTGROUP(J1)+N-1)   is the index of the inter-group transition state that connects
!                                group J1 to group GROUPCONN(STARTGROUP(J1)+N-1)
!
!  We first find an upper bound for the number of connections per group (including duplicates)
!  Then we calculate NCONNGROUP again removing the duplicates
!
NCONNGROUP(1:NEWNMIN)=0
DO J1=1,NTS
   CALL CHECKTS(ETS(J1),EMIN(PLUS(J1)),EMIN(MINUS(J1)),KPLUS(J1),KMINUS(J1),NCONN(PLUS(J1)),NCONN(MINUS(J1)), &
                PLUS(J1),MINUS(J1),.FALSE.,CUT_UNDERFLOW,LREJECTTS)
   IF (LREJECTTS) CYCLE
   IF (MINGROUP(PLUS(J1)).EQ.MINGROUP(MINUS(J1))) CYCLE ! Ignore intragroup rates
   NCONNGROUP(MINGROUP(PLUS(J1)))=NCONNGROUP(MINGROUP(PLUS(J1)))+1
   NCONNGROUP(MINGROUP(MINUS(J1)))=NCONNGROUP(MINGROUP(MINUS(J1)))+1
ENDDO
STARTGROUP(1)=1
DO J1=2,NEWNMIN
   STARTGROUP(J1)=STARTGROUP(J1-1)+NCONNGROUP(J1-1)
ENDDO
NCONNGROUP(1:NEWNMIN)=0
NEWNTS=0
tsloopagain: DO J1=1,NTS
   CALL CHECKTS(ETS(J1),EMIN(PLUS(J1)),EMIN(MINUS(J1)),KPLUS(J1),KMINUS(J1),NCONN(PLUS(J1)),NCONN(MINUS(J1)), &
                PLUS(J1),MINUS(J1),.FALSE.,CUT_UNDERFLOW,LREJECTTS)
   IF (LREJECTTS) CYCLE
   LOCALP=MINGROUP(PLUS(J1))
   LOCALM=MINGROUP(MINUS(J1))
   IF (LOCALP.EQ.LOCALM) CYCLE ! Ignore intragroup rates
   IF (NCONNGROUP(LOCALP).LT.NCONNGROUP(LOCALM)) THEN ! search the shorter list of connections
      DO J2=1,NCONNGROUP(LOCALP)
         IF (GROUPCONN(STARTGROUP(LOCALP)+J2-1).EQ.LOCALM) CYCLE tsloopagain ! we already have this connection
      ENDDO
   ELSE
      DO J2=1,NCONNGROUP(LOCALM)
         IF (GROUPCONN(STARTGROUP(LOCALM)+J2-1).EQ.LOCALP) CYCLE tsloopagain ! we already have this connection
      ENDDO
   ENDIF
   NEWNTS=NEWNTS+1
   NCONNGROUP(LOCALP)=NCONNGROUP(LOCALP)+1
   NCONNGROUP(LOCALM)=NCONNGROUP(LOCALM)+1
   GROUPCONN(STARTGROUP(LOCALP)+NCONNGROUP(LOCALP)-1)=LOCALM
   GROUPCONN(STARTGROUP(LOCALM)+NCONNGROUP(LOCALM)-1)=LOCALP
   GROUPTS(STARTGROUP(LOCALP)+NCONNGROUP(LOCALP)-1)=NEWNTS
   GROUPTS(STARTGROUP(LOCALM)+NCONNGROUP(LOCALM)-1)=NEWNTS
   NEWPLUS(NEWNTS)=LOCALP
   NEWMINUS(NEWNTS)=LOCALM
ENDDO tsloopagain
!
!  Change the double loop over NTS and NEWNTS to a single loop
!
NTSGROUP(1:NTS)=0
NEWKPLUS(1:NTS)=0.0D0
NEWKMINUS(1:NTS)=0.0D0
tsloop: DO J1=1,NTS
   CALL CHECKTS(ETS(J1),EMIN(PLUS(J1)),EMIN(MINUS(J1)),KPLUS(J1),KMINUS(J1),NCONN(PLUS(J1)),NCONN(MINUS(J1)), &
                PLUS(J1),MINUS(J1),.FALSE.,CUT_UNDERFLOW,LREJECTTS)
   IF (LREJECTTS) CYCLE
   LOCALP=MINGROUP(PLUS(J1))
   LOCALM=MINGROUP(MINUS(J1))
   IF (LOCALP.EQ.LOCALM) CYCLE ! Ignore intragroup rates
!
!  We know which groups this PE ts links. We just need to know which intra-group ts it
!  contributes to so we can add on the contribution.
!
   IF (NCONNGROUP(LOCALP).LT.NCONNGROUP(LOCALM)) THEN ! search the shorter list of connections
      DO J2=1,NCONNGROUP(LOCALP)
         IF (GROUPCONN(STARTGROUP(LOCALP)+J2-1).EQ.LOCALM) THEN
            NDUMMY=GROUPTS(STARTGROUP(LOCALP)+J2-1)
            GOTO 999
         ENDIF
      ENDDO
   ELSE
      DO J2=1,NCONNGROUP(LOCALM)
         IF (GROUPCONN(STARTGROUP(LOCALM)+J2-1).EQ.LOCALP) THEN
            NDUMMY=GROUPTS(STARTGROUP(LOCALM)+J2-1)
            GOTO 999
         ENDIF
      ENDDO
   ENDIF
   PRINT '(A,I8,A,2I8)','regrouppersist> ERROR - unmatched partner for pe ts ',J1,' groups linked: ', &
  &                LOCALP,LOCALM
   PRINT '(A,I8)','regrouppersist> NCONNGROUP for first of these groups:',NCONNGROUP(LOCALP)
   PRINT '(A)','regrouppersist> connected min:'
   PRINT '(16I8)',(GROUPCONN(STARTGROUP(LOCALP)+J2),J2=0,NCONNGROUP(LOCALP)-1)
   PRINT '(A)','regrouppersist> corresponding inter-group ts:'
   PRINT '(16I8)',(GROUPTS(STARTGROUP(LOCALP)+J2),J2=0,NCONNGROUP(LOCALP)-1)
   PRINT '(A,I8)','regrouppersist> NCONNGROUP for second of these groups:',NCONNGROUP(LOCALM)
   PRINT '(A)','regrouppersist> connected min:'
   PRINT '(16I8)',(GROUPCONN(STARTGROUP(LOCALM)+J2),J2=0,NCONNGROUP(LOCALM)-1)
   PRINT '(A)','regrouppersist> corresponding inter-group ts:'
   PRINT '(16I8)',(GROUPTS(STARTGROUP(LOCALM)+J2),J2=0,NCONNGROUP(LOCALM)-1)
   STOP
999 CONTINUE

!
! The PFMEAN terms cancel in the calculation of NEWKPLUS and NEWKMINUS
!
   IF ((LOCALP.EQ.NEWPLUS(NDUMMY)).AND.(LOCALM.EQ.NEWMINUS(NDUMMY))) THEN
      NEWKPLUS(NDUMMY)=NEWKPLUS(NDUMMY)+EXP(PFMIN(PLUS(J1))+KPLUS(J1))/NEWPFMIN(LOCALP)
      NEWKMINUS(NDUMMY)=NEWKMINUS(NDUMMY)+EXP(PFMIN(MINUS(J1))+KMINUS(J1))/NEWPFMIN(LOCALM)
      NTSGROUP(NDUMMY)=NTSGROUP(NDUMMY)+1 
   ELSEIF ((LOCALP.EQ.NEWMINUS(NDUMMY)).AND.(LOCALM.EQ.NEWPLUS(NDUMMY))) THEN
      NEWKPLUS(NDUMMY)=NEWKPLUS(NDUMMY)+EXP(PFMIN(MINUS(J1))+KMINUS(J1))/NEWPFMIN(LOCALM)
      NEWKMINUS(NDUMMY)=NEWKMINUS(NDUMMY)+EXP(PFMIN(PLUS(J1))+KPLUS(J1))/NEWPFMIN(LOCALP)
      NTSGROUP(NDUMMY)=NTSGROUP(NDUMMY)+1 
   ELSE
      PRINT '(A)','regrouppersist> ERROR - one of the two branches above should match'
      STOP
   ENDIF
ENDDO tsloop

PRINT '(A,I7)','regrouppersist> Number of intergroup transition states=',NEWNTS
DO J1=1,NEWNTS
   DO J2=1,NPERSIST
      DO J3=J2+1,NPERSIST
         IF (GROUPP(J2,NEWPLUS(J1)).AND.GROUPP(J3,NEWMINUS(J1))) THEN
            PRINT '(A,I7,A,2I7,2(A,G20.10))','regrouppersist> intergroup ts ',J1,' links persistent groups ',J2,J3, &
  &                      ' with k(B<-A)=', &
  &                      NEWKPLUS(J1),' k(A<-B)=',NEWKMINUS(J1)
            PRINT '(A,I7,A,I7,A)','regrouppersist> persistent group ',J2,' is number ',NEWPLUS(J1),' containing PE minima:'
            DO J4=1,NMIN
               IF (MINGROUP(J4).EQ.NEWPLUS(J1)) PRINT '(I7)',J4
            ENDDO
            PRINT '(A,I7,A,I7,A)','regrouppersist> persistent group ',J3,' is number ',NEWMINUS(J1),' containing PE minima:'
            DO J4=1,NMIN
               IF (MINGROUP(J4).EQ.NEWMINUS(J1)) PRINT '(I7)',J4
            ENDDO
         ENDIF
      ENDDO
   ENDDO
   IF (DEBUG) PRINT '(3(A,I7),2(A,G20.10),A,I7)','regrouppersist> Grouped ts ',J1,' between minima groups ',NEWPLUS(J1), &
  &    ' and ',NEWMINUS(J1), &
  &    ' k+=',NEWKPLUS(J1),' k-=',NEWKMINUS(J1),' members=',NTSGROUP(J1)
!
! PFMEAN is a constant factor here
!
!  PRINT '(A,2G20.10)','regrouppersist> detailed balance - these numbers should be equal: ',NEWKPLUS(J1)*NEWPFMIN(NEWPLUS (J1)), &
! &                                                                                      NEWKMINUS(J1)*NEWPFMIN(NEWMINUS(J1))

   IF ((NEWKPLUS(J1).EQ.0.0D0).OR.(NEWKMINUS(J1).EQ.0.0D0)) THEN
      IF (DEBUG) PRINT '(A,I7,2G20.10)','regrouppersist> WARNING - J1,NEWKPLUS,NEWKMINUS=',J1,NEWKPLUS(J1),NEWKMINUS(J1)
      NEWETS(J1)=HUGE(1.0D0)
   ELSE
      NEWETS(J1)=NEWEMIN(NEWPLUS(J1))-TEMPERATURE*(LOG(NEWKPLUS(J1))+LOG(PLANCK/TEMPERATURE))
      IF (DEBUG) PRINT '(3(A,G20.10))','regrouppersist> Grouped ts free energy=', &
  &                     NEWETS(J1), &
  &              ' or ',NEWEMIN(NEWMINUS(J1))-TEMPERATURE*(LOG(NEWKMINUS(J1))+LOG(PLANCK/TEMPERATURE)), &
  &              ' or ',NEWEMIN(NEWPLUS(J1))-TEMPERATURE*(LOG(NEWKPLUS(J1))+LOG(PLANCK/TEMPERATURE))
      IF (NEWETS(J1).NE.0.0D0) THEN ! Check for consistency
         IF (ABS((NEWETS(J1)-NEWEMIN(NEWMINUS(J1))+TEMPERATURE*(LOG(NEWKMINUS(J1))+ &
  &               LOG(PLANCK/TEMPERATURE)))/NEWETS(J1)).GT.0.01D0) THEN
            PRINT '(A,I7,A,2G20.10)','regrouppersist> WARNING - free energies for ts group ',J1,' are ',  &
  &                NEWETS(J1),NEWEMIN(NEWMINUS(J1))-TEMPERATURE*(LOG(NEWKMINUS(J1))+LOG(PLANCK/TEMPERATURE))
            STOP
         ENDIF
      ENDIF
   ENDIF
ENDDO

IF (.NOT.ONEREGROUPT) GOTO 888 ! free energies etc. have been updated - go back to see if we can regroup further
777 CONTINUE ! End of iterative regrouping loop.
PRINT '(A,I7,A,F15.5)','regrouppersist> Final number of free energy groups is ',NGROUPS,' at T=',TEMPERATURE
PRINT '(A)',' '
!
! We now have all the free energies for minima and transition state groups. 
! Some of the original groups for the minima will generally be empty, so 
! now we renumber.
!
! POINTERS are renumbered by calling REGROUP if required in GT, not here
!
NCOUNT=0
NDUMMY=0
NEWA=0; NEWB=0
DO J1=1,NEWNMIN
   IF (NMINGROUP(J1).GT.0) THEN
      NCOUNT=NCOUNT+1
      GROUPMAP(J1)=NCOUNT
      NCONNGROUP(NCOUNT)=NCONNGROUP(J1)
      STARTGROUP(NCOUNT)=STARTGROUP(J1)
      NMINGROUP(NCOUNT)=NMINGROUP(J1)
      NEWEMIN(NCOUNT)=NEWEMIN(J1)
      NEWPFMIN(NCOUNT)=NEWPFMIN(J1)
      GROUPA(NCOUNT)=GROUPA(J1)
      GROUPB(NCOUNT)=GROUPB(J1)
      IF (GROUPA(J1)) NEWA=NEWA+1  
      IF (GROUPB(J1)) NEWB=NEWB+1  
      NDUMMY=NDUMMY+NMINGROUP(NCOUNT)
   ENDIF
ENDDO
NGROUPS=NCOUNT
PRINT '(4(A,I7))','regrouppersist> Number of groups after removing empty sets=',NCOUNT, &
  &         ' total PE minima=',NDUMMY,' # A: ',NEWA,' # B: ',NEWB
IF (NDUMMY.NE.NMINCONNECTED) THEN
   PRINT '(A,I7)','regrouppersist> ERROR - number of minima in groups should be ',NMINCONNECTED
   STOP
ENDIF
PRINT '(A)','regrouppersist> Renumbering free energy minima and ts to remove empty sets'
!
! Excluded minima belong to group 0.
!
DO J1=1,NMIN
   IF (MINGROUP(J1).EQ.0) CYCLE
   MINGROUP(J1)=GROUPMAP(MINGROUP(J1))
ENDDO
!
!  Dump the members of the free energy groups in terms of pe stationary points
!
IF (GETPAIRST.OR.DUMPGROUPST) THEN
   IF (DUMPGROUPST) OPEN(UNIT=1,FILE='minima_groups.'//ADJUSTL(TRIM(TSTRING)),STATUS='UNKNOWN')
   FREEMINPOINT(0)=1
   NDUMMY=1
   DO J1=1,NGROUPS
      NCOUNT=0
      DO J2=1,NMIN
         IF (MINGROUP(J2).EQ.J1) THEN
            NCOUNT=NCOUNT+1
            FREEMINLIST(NDUMMY)=J2
            NDUMMY=NDUMMY+1
!           IF (DUMPGROUPST) WRITE(1,'(I8)',ADVANCE='NO') J2
            IF (DUMPGROUPST) WRITE(1,'(I8)',ADVANCE='YES') J2
         ENDIF
      ENDDO
      IF (DUMPGROUPST) WRITE(1,'(A,I8,A,G20.10,A,I7)') 'group ',J1,' free energy=',NEWEMIN(J1),' pe minima=',NCOUNT
      FREEMINPOINT(J1)=NDUMMY-NCOUNT
!     PRINT '(A,I8,A,I8)','regrouppersist> free energy group ',J1,' starts at FREEMINLIST entry ',FREEMINPOINT(J1)
      IF (DUMPGROUPST) WRITE(1,'(A)') ' '
   ENDDO
   FREEMINPOINT(NGROUPS+1)=NDUMMY ! needed to define the last entry for group NGROUPS below
   IF (DUMPGROUPST) CLOSE(1)
   
   FREETSPOINT(0)=1
   NDUMMY=1
   IF (DUMPGROUPST) OPEN(UNIT=1,FILE='ts_groups.'//ADJUSTL(TRIM(TSTRING)),STATUS='UNKNOWN')
   DO J1=1,NEWNTS
      IF (DUMPGROUPST) WRITE(1,'(A,I8,A,G20.10,A,2I8)') 'ts group ',J1,' free energy=',NEWETS(J1), &
  &             ' links groups: ',NEWPLUS(J1),NEWMINUS(J1)
      NCOUNT=0
      DO J2=1,NTS
         IF (((MINGROUP(PLUS(J2)).EQ.NEWPLUS(J1)).AND.(MINGROUP(MINUS(J2)).EQ.NEWMINUS(J1))).OR. &
     &       ((MINGROUP(PLUS(J2)).EQ.NEWMINUS(J1)).AND.(MINGROUP(MINUS(J2)).EQ.NEWPLUS(J1))))  THEN
!           IF (DUMPGROUPST) WRITE(1,'(I8)',ADVANCE='NO') J2
            IF (DUMPGROUPST) WRITE(1,'(I8)',ADVANCE='YES') J2
            NCOUNT=NCOUNT+1
            FREETSLIST(NDUMMY)=J2
            NDUMMY=NDUMMY+1
         ENDIF
      ENDDO
      FREETSPOINT(J1)=NDUMMY-NCOUNT
!     PRINT '(A,I8,A,I8)','regrouppersist> free energy ts ',J1,' starts at FREETSLIST entry ',FREETSPOINT(J1)
      IF (DUMPGROUPST) WRITE(1,'(A)') ' '
   ENDDO 
   IF (DUMPGROUPST) CLOSE(1)
ENDIF
!
!  Dump the members of the free energy groups separately for A, B, and the
!  persistent minima.
!
IF (DUMPGROUPST) THEN
   LUNIT=GETUNIT()
   OPEN(LUNIT,FILE='mrp.A',STATUS='UNKNOWN')
   DO J1=1,NMIN
      IF (MINGROUP(J1).EQ.0) CYCLE
      IF (GROUPA(MINGROUP(J1))) WRITE(LUNIT,'(I8)') J1
   ENDDO
   CLOSE(LUNIT)

   LUNIT=GETUNIT()
   OPEN(LUNIT,FILE='mrp.B',STATUS='UNKNOWN')
   DO J1=1,NMIN
      IF (MINGROUP(J1).EQ.0) CYCLE
      IF (GROUPB(MINGROUP(J1))) WRITE(LUNIT,'(I8)') J1
   ENDDO
   CLOSE(LUNIT)

   DO J3=1,NPERSIST
      LUNIT=GETUNIT()
      WRITE(FNAME,'(I6)') J3
      WRITE(FNAME,'(A)') 'mrp.' // TRIM(ADJUSTL(FNAME))
      OPEN(LUNIT,FILE=FNAME,STATUS='UNKNOWN')
      DO J1=1,NMIN
         IF (MINGROUP(J1).EQ.0) CYCLE
         IF (GROUPP(J3,MINGROUP(J1))) WRITE(LUNIT,'(I8)') J1
      ENDDO
      CLOSE(LUNIT)
   ENDDO
ENDIF
!
! From here on down we overwrite the PE groups with free energy groups.
! Everything goes over to the free energy group scenario, including rate
! constants, etc. We cannot continue growing a database after this, but we
! can calculate global rate constants or run DIJKSTRA or KSHORTEST paths based
! on free energy rather than pe groups.
!

NMINA=0; NMINB=0
NMINP(1:NPERSIST)=0

DO J2=1,NPERSIST
   DO J1=1,NEWNMIN
      PFMIN(J1)=LOG(NEWPFMIN(J1))  
      EMIN(J1)=NEWEMIN(J1)  
      IF (GROUPP(J2,J1)) THEN
         NMINP(J2)=NMINP(J2)+1  
         LOCATIONP(NMINP(J2),J2)=J1  
      ENDIF
   ENDDO
ENDDO
DO J1=1,NEWNMIN
   PFMIN(J1)=LOG(NEWPFMIN(J1))  
   EMIN(J1)=NEWEMIN(J1)  
   IF (GROUPA(J1)) THEN
      NMINA=NMINA+1
      LOCATIONA(NMINA)=J1
      ENDIF
   IF (GROUPB(J1)) THEN
      NMINB=NMINB+1
      LOCATIONB(NMINB)=J1
   ENDIF
ENDDO

NTS=NEWNTS

DO J1=1,NEWNTS
   PLUS(J1)=GROUPMAP(NEWPLUS(J1))
   MINUS(J1)=GROUPMAP(NEWMINUS(J1))
   IF (NEWKPLUS(J1).GT.0.0D0) THEN
      KPLUS(J1)=LOG(NEWKPLUS(J1))
   ELSE
      KPLUS(J1)=-HUGE(1.0D0)
   ENDIF
   IF (NEWKMINUS(J1).GT.0.0D0) THEN
      KMINUS(J1)=LOG(NEWKMINUS(J1))
   ELSE
      KMINUS(J1)=-HUGE(1.0D0)
   ENDIF
   ETS(J1)=NEWETS(J1)
ENDDO
PFTOTALA=0.0D0
PFTOTALB=0.0D0
DO J1=1,NMINB
   PFTOTALB=PFTOTALB+EXP(PFMIN(LOCATIONB(J1)))
ENDDO
PFTOTALB=LOG(PFTOTALB)
 DO J1=1,NMINA
   PFTOTALA=PFTOTALA+EXP(PFMIN(LOCATIONA(J1)))
ENDDO
PFTOTALA=LOG(PFTOTALA)
NMIN=NGROUPS
!
!  If we are going to analyse the min.data.regrouped.resorted and ts.data.regrouped.resorted
!  files for rates subsequently, then we have to arrange for the ln products of frequencies
!  to give us a factor of (kT/h). This can be done by setting the ln product equal to zero
!  for the transition state and 2 * ln(2*Pi*k*T/h) for all the minima. We already have h in the
!  units of kT, so this is easy. The 2*Pi factor occurs because the frequencies are assumed to be
!  angular normal mmode frequencies, and the factor of two occurs because they are assumed
!  to be squared.
!
LNPROD=2.0D0*LOG(2.0D0*3.141592654D0*TEMPERATURE/PLANCK)

OPEN(UNIT=1,FILE='min.data.regrouped.'//ADJUSTL(TRIM(TSTRING)),STATUS='UNKNOWN')
DO J1=1,NMIN
   WRITE(1,'(2G20.10,I7,4F20.10)') EMIN(J1),LNPROD,1,1.0,1.0,1.0,0.0
ENDDO
CLOSE(1)
OPEN(UNIT=1,FILE='ts.data.regrouped.'//ADJUSTL(TRIM(TSTRING)),STATUS='UNKNOWN')
DO J1=1,NTS
   WRITE(1,'(2G20.10,3I10,3F20.10)') ETS(J1),0.0,1,PLUS(J1),MINUS(J1),1.0,1.0,1.0
ENDDO
CLOSE(1)
OPEN(UNIT=1,FILE='min.A.regrouped.'//ADJUSTL(TRIM(TSTRING)),STATUS='UNKNOWN')
WRITE(1,'(I7)') NMINA
DO J1=1,NMINA
   WRITE(1,'(I7)') LOCATIONA(J1)
ENDDO
CLOSE(1)
OPEN(UNIT=1,FILE='min.B.regrouped.'//ADJUSTL(TRIM(TSTRING)),STATUS='UNKNOWN')
WRITE(1,'(I7)') NMINB
DO J1=1,NMINB
   WRITE(1,'(I7)') LOCATIONB(J1)
ENDDO
CLOSE(1)

PRINT '(A)','regrouppersist>  NOTE: from here on down min and ts refer to the new groups!'

DEALLOCATE(NPGROUP,NEWNMINP,GROUPP,ISP)

RETURN

END SUBROUTINE REGROUPPERSIST
