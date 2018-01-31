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

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
! "kth shortest paths" analysis subroutine via the Recursive Enumeration Algorithm of Jimenez and Marzal
!
!-------------------------------------------------------------------------------
! Module containing the network information
!
MODULE GRAPH
IMPLICIT NONE
SAVE

TYPE edge
   INTEGER :: index                ! index of ts, corresponding to line number in ts.data
   DOUBLE PRECISION :: l           ! weight of edge
   TYPE(node), POINTER :: from     ! edge going from "from" to "to"
   TYPE(node), POINTER :: to
   TYPE(edge), POINTER :: nextpred ! next in list of incoming edges associated with "to"
   TYPE(edge), POINTER :: nextsucc ! next in list of outgoing edges associated with "from"
END TYPE edge

TYPE path
   TYPE(node), POINTER :: node     ! node
   DOUBLE PRECISION :: l           ! accumulated weight of path
   INTEGER :: k                    ! if defined, this path is the k-th shortest path
   INTEGER :: i                    ! index in the now-linear candidate_paths array 
   TYPE(path), POINTER :: pred     ! predecessor in path
   TYPE(path), POINTER :: next     ! next in list of k-th shortest path
                                   ! or next in list of candidates
END TYPE path

TYPE node
   INTEGER :: index                ! index of minimum, corresponding to line number in min.data
   LOGICAL :: permanent            ! used in Dijkstra
   TYPE(edge), POINTER :: pred     ! top of list of incoming edges
   TYPE(edge), POINTER :: succ     ! top of list of outgoing edges
   TYPE(path), POINTER :: shortest ! top of list of k-th shortest paths
   TYPE(path), POINTER :: cand     ! top of list of candidates
END TYPE node

TYPE(edge), TARGET, ALLOCATABLE :: tsedge(:)
TYPE(node), TARGET, ALLOCATABLE :: minnode(:)
TYPE(path), TARGET, ALLOCATABLE :: shortest_paths(:,:), candidate_paths(:)

INTEGER, ALLOCATABLE :: index_next(:), index_start(:)

DOUBLE PRECISION :: start_shift_factor

END MODULE GRAPH

!-------------------------------------------------------------------------------
SUBROUTINE KSHORTESTPATHS(NWORST,GETCANDIDATES,PAIRSTODO,NMINSAVE,MINMAP)
USE PORFUNCS
USE COMMONS
USE GRAPH
IMPLICIT NONE

INTERFACE
   RECURSIVE SUBROUTINE print_path(p)
   USE graph, ONLY : path, edge, node
   IMPLICIT NONE
   TYPE(path), POINTER :: p
   END SUBROUTINE print_path
END INTERFACE

TYPE(path), POINTER :: pp

LOGICAL, INTENT(IN) :: GETCANDIDATES
LOGICAL DEADTS(NTS), ISA(NMIN), ISB(NMIN), CHANGED, ISSTART(NMIN), CHECKCONN, ADDMIN(NMIN), ADDTS(NTS)
LOGICAL SAVE_DEADTS(NTS)

INTEGER, PARAMETER :: MAXPATH=10000 ! longest path length allowed - surely 10000 is enough?!
INTEGER J1, J2, J3, PARENT(NMIN), J5, LJ1, LJ2, VERYBESTSTART, ENDMIN
INTEGER K1, K2, K3, KSHORT, NCOUNT
INTEGER VERYBESTEND, J6, NWORST, PAIRSTODO, NDUMMY, NDUMMY2
INTEGER VERYBESTPARENT(NMIN), NMINSTART, NMINEND, NDISTEND(NMIN), NSTEPS, TSPATHID(MAXPATH), MINPATHID(MAXPATH), NTSMIN
INTEGER NCOL(NMIN), NVAL(NCONNMAX,NMIN), NDEAD, NDISTA(NMIN), NDISTB(NMIN), NCYCLE, DMIN, ISTAT
INTEGER INDEX_TS(NCONNMAX,NMIN), INDEX_NEW_MIN(NMIN)
INTEGER DMAX, NUNCONA, NUNCONB, COUNTER
INTEGER, ALLOCATABLE :: LOCATIONSTART(:), LOCATIONEND(:)

DOUBLE PRECISION EMKSUM(NMIN), WEIGHT(NMIN), VERYBEST, ETSMIN, TNEW, ELAPSED, WAIT, WAITAB, WAITBA
DOUBLE PRECISION DMATMC(NTS,2), DUMMY, PFTOTALSTART, HUGESAVE, LOCALPOINTS(3*NATOMS), KSUM(NMIN)

INTEGER NMINSAVE, MINMAP(NMIN), DHORDERTS, DPLUS, DMINUS, NNMINA, NNMINB
INTEGER, ALLOCATABLE :: NEWINDEX(:)
LOGICAL, ALLOCATABLE :: INCLUDEMIN(:)
CHARACTER(LEN=130) DUMMYSTRING
DOUBLE PRECISION DETS,DFVIBTS,DIXTS,DIYTS,DIZTS,DNEGEIG, CUT_UNDERFLOW

CALL CPU_TIME(ELAPSED)
ISA(1:NMIN)=.FALSE.
ISB(1:NMIN)=.FALSE.
ADDMIN(1:NMIN)=.FALSE.
ADDTS(1:NTS)=.FALSE.
DO J1=1,NMINA
   ISA(LOCATIONA(J1))=.TRUE.
ENDDO
DO J1=1,NMINB
   ISB(LOCATIONB(J1))=.TRUE.
ENDDO

!  Original KSUM values are restored before return so that DIJKSTRA can be
!  called more than once.

IF (REGROUPFREET .OR. REGROUPFREEABT .OR. REGROUPPERSISTT) THEN
   TSTHRESH = HUGE(0.0D0)
   MAXBARRIER=HUGE(1.0D0)
END IF
CUT_UNDERFLOW=-300.0D0
CALL RATECONST_SETUP(KSUM,DEADTS,NDEAD,.FALSE.,CUT_UNDERFLOW)

DO J1=1,NMIN
   EMKSUM(J1)=EXP(-KSUM(J1))
ENDDO

IF (REGROUPFREET .OR. REGROUPPERSISTT) THEN
   IF (.NOT. ALLOCATED(NCONNGROUP)) THEN
      PRINT *,'kshortestpaths> NCONNGROUP is not allocated, yet REGROUPFREET is .TRUE.; stopping.'
      STOP
   ENDIF
   NCOUNT = 0
   DO k1 = 1, NMIN
      NCOUNT = NCOUNT + NCONNGROUP(k1) ! doesn't depend upon the order
   ENDDO
   IF (NCOUNT /= 2*NTS) THEN
      PRINT *,'kshortestpaths> stopping due to unexpected error: NCOUNT = ',NCOUNT,' but 2*NTS = ',2*NTS
      STOP
   ENDIF
ENDIF

IF (ALLOCATED(minnode)) DEALLOCATE(MINNODE)
IF (ALLOCATED(index_next)) DEALLOCATE(index_next)
IF (ALLOCATED(index_start)) DEALLOCATE(index_start)
ALLOCATE(minnode(NMIN), index_next(NMIN), index_start(NMIN))
IF (ALLOCATED(tsedge)) DEALLOCATE(tsedge)
ALLOCATE(tsedge(2*NTS))
IF (ALLOCATED(shortest_paths)) DEALLOCATE(shortest_paths)
IF (ALLOCATED(candidate_paths)) DEALLOCATE(candidate_paths)
ALLOCATE(shortest_paths(NMIN, NPATHS), candidate_paths(2*NTS))

VERYBEST=-1.0D0

!!!!!!!!!!!!!!!!!!!   Dijkstra calculation    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Find the single paths from a to b or b to a with the largest product of
!  branching probabilities. 
!
!  We need to allow transitions to all A and B minima - DMATMC is different from P^fold part
!
SAVE_DEADTS(:) = DEADTS(:)
! Added the above line as deadts will be used in MAKED4 to arrange that only one TS appears to
! connect a pair of minima, and that TS carries the sum of the rate consts for all the TS's
! that directly connect the relevant minima.  Need to revert to the normal deadts when we're
! reconstructing the paths with the lowest TS between each pair of minima.

CALL MAKED4(DMATMC,NCOL,NVAL,DEADTS,KSUM,INDEX_TS)
!
!  Check that the stationary point database is actually connected, 
!  and remove minima that lie in disjoint graphs.
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
PRINT '(3(A,I8))','kshortestpaths> steps to A region converged in ',NCYCLE-1, &
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
PRINT '(3(A,I8))','kshortestpaths> steps to B region converged in ',NCYCLE-1, &
&                    ' cycles; maximum=',DMAX,' disconnected=',NUNCONB
!
!  This could happen if disconnected minima lie in the A or B region
!
IF (NUNCONB.NE.NUNCONA) PRINT '(A)', &
&                   'kshortestpaths> WARNING - number of disconnected minima from A and B is different'
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
   PRINT '(A)','kshortestpaths> There is no connection between the A and B regions'
!  OPEN(UNIT=1,FILE='ts.attempts',STATUS='UNKNOWN')
!  WRITE(1,'(I8)') TSATTEMPT(1:NTS)
!  CLOSE(1)
   STOP
ENDIF

IF (DIRECTION.EQ.'AB') THEN
   NMINSTART=NMINB
   NMINEND=NMINA
   IF (ALLOCATED(LOCATIONSTART)) DEALLOCATE(LOCATIONSTART)
   IF (ALLOCATED(LOCATIONEND)) DEALLOCATE(LOCATIONEND)
   ALLOCATE(LOCATIONSTART(NMINB),LOCATIONEND(NMINA))
   LOCATIONSTART(1:NMINB)=LOCATIONB(1:NMINB)
   LOCATIONEND(1:NMINA)=LOCATIONA(1:NMINA)
   NDISTEND(1:NMIN)=NDISTA(1:NMIN)
   ISSTART(1:NMIN)=ISB(1:NMIN)
   PFTOTALSTART=PFTOTALB
ELSEIF (DIRECTION.EQ.'BA') THEN
   NMINSTART=NMINA
   NMINEND=NMINB
   IF (ALLOCATED(LOCATIONSTART)) DEALLOCATE(LOCATIONSTART)
   IF (ALLOCATED(LOCATIONEND)) DEALLOCATE(LOCATIONEND)
   ALLOCATE(LOCATIONSTART(NMINA),LOCATIONEND(NMINB))
   LOCATIONSTART(1:NMINA)=LOCATIONA(1:NMINA)
   LOCATIONEND(1:NMINB)=LOCATIONB(1:NMINB)
   NDISTEND(1:NMIN)=NDISTB(1:NMIN)
   ISSTART(1:NMIN)=ISA(1:NMIN)
   PFTOTALSTART=PFTOTALA
ENDIF
!
!  Find largest weight for each B(A) minimum to all A(B) minima.
!
ENDMIN=NMINSTART
NWORST=0
loopstart: DO J1=1,ENDMIN
   LJ1=LOCATIONSTART(J1)
   IF (NCOL(LJ1).LE.0) CYCLE loopstart
   IF (NDISTEND(LJ1).EQ.1000000) CYCLE loopstart
   IF (NDISTEND(LJ1).EQ.0) CYCLE loopstart

   loopend: DO K2 = 1, NMINEND ! cycle over end minima.

      ! jmc start of REA stuff here. 
      ! In order to ensure that the edge weights for TS's involving "start" are always negative, we include the constant
      ! factor of minim(start)%adjsum as for all other edges.  Need to remove this though in the final shortest-path rate
      ! constants.  N.B. start = LJ1 and end = LJ2

      LJ2=LOCATIONEND(K2)
      IF (ISSTART(LJ2)) CYCLE ! end point is also a start minimum...
      IF (NCOL(LJ2).EQ.0) CYCLE ! end point has no connections.
!
! jmc test for and remove any minima that are not accessible from the reactants 
! APART from through the product set (which is bad for the k-th shortest paths algorithm...)
!
! Calculate the minimum number of steps from every minimum to the reactant, in the absence of the product.
      NDISTA(1:NMIN)=1000000
      NDISTA(LJ1)=0
      NCYCLE=0
511   CHANGED=.FALSE.
      DO K1=1,NMIN
         IF (K1 == LJ1) CYCLE ! reactant minimum
         IF (K1 == LJ2) CYCLE ! product minimum
         innerloop: DO K3=1,NCOL(K1)
            IF (NVAL(K3,K1) == LJ2) CYCLE innerloop ! product min
            IF (NDISTA(NVAL(K3,K1))+1.LT.NDISTA(K1)) THEN
               CHANGED=.TRUE.
               NDISTA(K1)=NDISTA(NVAL(K3,K1))+1
            ENDIF
         ENDDO innerloop
      ENDDO
      IF (CHANGED) GOTO 511

      DO K1=1,NMIN
         IF (K1 == LJ2) CYCLE
         IF ((NDISTA(K1).EQ.1000000)) THEN
            IF (DEBUG) PRINT*,'kshortestpaths> minimum ',K1,' can only be reached from reactant via product; deleting.'
            DO K3=1,NCOL(K1)
               IF (NVAL(K3,K1) == LJ2) THEN
                  DEADTS(INDEX_TS(K3,K1)) = .TRUE. ! remove all connections to the product set.
                  ! We don't care if connections within this strange region remain, as we will never be able 
                  ! to reach them from either reactant or product after we disconnect them from product.
                  ! No need to change the elements of DMATMC, as with DEADTS set appropriately, they should never be 
                  ! used during the Dijkstra or REA algorithms.
                  ! Also note that no edges leave the product, and for the edges entering the product the DMATMC entry is 
                  ! k_out of non-prod-minimum / sum of k's out of non-product-minimum.
               ENDIF
            ENDDO
         ENDIF
      ENDDO
!
! end jmc extra stuff
!
      start_shift_factor = KSUM(LJ1)

! Now convert the minim and ts arrays into the pointer-based derived data types.

      DO k1 = 1, NMIN
         NULLIFY(minnode(k1)%pred,minnode(k1)%succ,minnode(k1)%shortest,minnode(k1)%cand)
      ENDDO
      DO k1 = 1, 2*NTS
         NULLIFY(tsedge(k1)%from,tsedge(k1)%to,tsedge(k1)%nextpred,tsedge(k1)%nextsucc)
      ENDDO

      DO k1 = 1, NMIN
         minnode(k1)%index = 0
         IF (NCOL(K1) == 0) CYCLE
         minnode(k1)%index = k1 ! Set the index of this minimum to its position in the minim array
      ENDDO
 
      DO k1 = 1, NTS
         tsedge(k1)%index = 0
         IF (DEADTS(K1)) CYCLE
         !tsedge(k1)%index = k1

         IF ( (PLUS(K1) == LJ1) .OR. (MINUS(K1) == LJ2) ) THEN
            ! the forwards edge corresponding to TS k1
            tsedge(k1)%l = -LOG(DMATMC(K1,1))
            !tsedge(k1)%l = DMATMC(K1,1)
            tsedge(k1)%from => minnode(PLUS(K1))
            tsedge(k1)%to => minnode(MINUS(K1))
         ELSEIF ( (MINUS(K1) == LJ1) .OR. (PLUS(K1) == LJ2) ) THEN
            ! the backwards edge corresponding to TS k1, now at k1 + nts in tsedge array.
            tsedge(nts + k1)%l = -LOG(DMATMC(K1,2))
            !tsedge(nts + k1)%l = DMATMC(K1,2)
            tsedge(nts + k1)%from => minnode(MINUS(K1))
            tsedge(nts + k1)%to => minnode(PLUS(K1))
         ELSE
            ! the forwards edge corresponding to TS k1
            tsedge(k1)%l = -LOG(DMATMC(K1,1))
            !tsedge(k1)%l = DMATMC(K1,1)
            tsedge(k1)%from => minnode(PLUS(K1))
            tsedge(k1)%to => minnode(MINUS(K1))

            ! the backwards edge corresponding to TS k1, now at k1 + nts in tsedge array.
            tsedge(nts + k1)%l = -LOG(DMATMC(K1,2))
            !tsedge(nts + k1)%l = DMATMC(K1,2)
            tsedge(nts + k1)%from => minnode(MINUS(K1))
            tsedge(nts + k1)%to => minnode(PLUS(K1))
         ENDIF
         IF (tsedge(k1)%l < 0.0D0) tsedge(k1)%l = 0.0D0
         IF (tsedge(nts + k1)%l < 0.0D0) tsedge(nts + k1)%l = 0.0D0
      ENDDO

      DO k1 = 1, NMIN
         !minnode(k1)%index = 0
         IF (NCOL(K1) == 0) CYCLE
         !minnode(k1)%index = k1 ! Set the index of this minimum to its position in the minim array
         DO k3 = 1, NCOL(K1)
            IF (ASSOCIATED(tsedge(INDEX_TS(K3,K1))%to)) THEN
               IF (tsedge(INDEX_TS(K3,K1))%to%index == k1) THEN
                  CALL add_pred(k1, INDEX_TS(K3,K1))
               ELSEIF (tsedge(INDEX_TS(K3,K1))%from%index == k1) THEN
                  CALL add_succ(k1, INDEX_TS(K3,K1))
               ELSE
                  PRINT *,'Something sadly awry here1'
               ENDIF
            ENDIF
            IF (ASSOCIATED(tsedge(INDEX_TS(K3,K1) + nts)%to)) THEN
               IF (tsedge(INDEX_TS(K3,K1) + nts)%to%index == k1) THEN
                  CALL add_pred(k1, INDEX_TS(K3,K1) + nts)
               ELSEIF (tsedge(INDEX_TS(K3,K1) + nts)%from%index == k1) THEN
                  CALL add_succ(k1, INDEX_TS(K3,K1) + nts)
               ELSE
                  PRINT *,'Something sadly awry here2'
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      HUGESAVE = HUGE(0.0D0)

      ! Do we need to reinitialize everything in NEXT_PATH at each pass through this do loop or not? 
      ! Yes, as the pattern of edges has changed....
      index_next(1) = 1
      DO k3 = 1, NPATHS
         NULLIFY(shortest_paths(1,k3)%pred,shortest_paths(1,k3)%next,shortest_paths(1,k3)%node)
         shortest_paths(1,k3)%l = HUGESAVE
      ENDDO
      IF (REGROUPFREET .OR. REGROUPPERSISTT) THEN
         DO k1 = 2, NMIN
            index_next(k1) = index_next(k1-1) + NCONNGROUP(MINMAP(k1-1))
            DO k3 = 1, NPATHS
               NULLIFY(shortest_paths(k1,k3)%pred,shortest_paths(k1,k3)%next,shortest_paths(k1,k3)%node)
               shortest_paths(k1,k3)%l = HUGESAVE
            ENDDO
         ENDDO
      ELSE
         ! To be safe, I'm including all minima here, and not excluding any on the grounds of nconn <= nconnmin, etc.
         ! Also note that because some of the TS's will probably have been excluded (because they mediate 
         ! degenerate rearrangements, or have too high energies), the sum of the elements of NCONN will often be less 
         ! than 2*NTS, the dimension of candidate_paths(). 
         DO k1 = 2, NMIN
            index_next(k1) = index_next(k1-1) + NCONN(k1-1)
            DO k3 = 1, NPATHS
               NULLIFY(shortest_paths(k1,k3)%pred,shortest_paths(k1,k3)%next,shortest_paths(k1,k3)%node)
               shortest_paths(k1,k3)%l = HUGESAVE
            ENDDO
         ENDDO
      ENDIF

      index_start = index_next

      DO k1 = 1, 2*NTS
         NULLIFY(candidate_paths(k1)%pred,candidate_paths(k1)%next,candidate_paths(k1)%node)
      ENDDO

      ! Use Dijkstra's algorithm to find the shortest path from 'start' to every accessible node in the graph.
      CALL MY_DIJKSTRA(LJ1, NMIN)

      DO K1 = 2, NPATHS
         CALL NEXT_PATH(LJ2, K1, LJ1) ! index of end min, which shortest path, index of start min.
      ENDDO

      ! Print information about the k-th shortest paths
      !CALL PRINT_PATH(minnode(LJ2)%shortest)

      NULLIFY(pp)
      pp => minnode(LJ2)%shortest
      DO
         IF (.NOT. ASSOCIATED(pp%next)) EXIT
         pp => pp%next
      ENDDO

      IF (pp%k /= 1) THEN 
         print *,'something wrong with the shortest shortest path'
         STOP
      ENDIF

      IF (pp%l == HUGESAVE) CYCLE ! end min LJ2 is disconnected from start min LJ1, hence from other end minima.
!
!  Now we loop over the NPATHS fastest paths.
!
      DO KSHORT=1,NPATHS
         DUMMY=1.0D0
         J5=LJ2
         NSTEPS=0
         PARENT = 0

!  Make the array PARENT for this particular path, then we don`t need to change it below!
         ! Descend the shortest paths list for end point LJ2 to find the KSHORT'th.  Note that the 
         ! shortest shortest path is at the bottom of this list.
         NULLIFY(pp)
         pp => minnode(LJ2)%shortest
         DO
            IF (pp%k == KSHORT) EXIT
            IF (.NOT. ASSOCIATED(pp%next)) THEN
               PRINT *,'kshortestpaths> something wrong with minnode(LJ2)%shortest ',LJ2
               STOP
            ENDIF
            pp => pp%next
         ENDDO

         PRINT *,'kshortestpaths> shortest path: ',pp%k,' for start min ',LJ1,' and end min ',LJ2
         WEIGHT(LJ2) = pp%l
         IF (pp%l == HUGESAVE) CYCLE loopend

         ! Now follow the required kth shortest path to find the predecessors (= parents).
         DO
            IF (.NOT. ASSOCIATED(pp%pred)) EXIT
            PARENT(pp%node%index) = pp%pred%node%index

            ! Identify the branching probability to calculate the weight of this path.
            DO J3=1,NCOL(PARENT(pp%node%index))
               IF (NVAL(J3,PARENT(pp%node%index)).EQ.pp%node%index) THEN
                   IF (PLUS(INDEX_TS(J3,PARENT(pp%node%index))) == PARENT(pp%node%index)) THEN
                      DUMMY=DUMMY*DMATMC(INDEX_TS(J3,PARENT(pp%node%index)),1)
                   ELSEIF (MINUS(INDEX_TS(J3,PARENT(pp%node%index))) == PARENT(pp%node%index)) THEN
                      DUMMY=DUMMY*DMATMC(INDEX_TS(J3,PARENT(pp%node%index)),2)
                   ELSE
                      PRINT *,'kshortestpaths> error with INDEX_TS or PARENT'
                      STOP
                   ENDIF
                   EXIT
               ENDIF
            ENDDO

            pp => pp%pred
         ENDDO

         DO 
!
!  Now identify the lowest transition state between the two minima.
!
            NSTEPS=NSTEPS+1
            ETSMIN=1.0D100
            J6=TOPPOINTER(J5) ! sets J6 to the transition state connected to J5 with the
                              ! highest number that isn't DEADTS
            pointa: DO
               IF (PLUS(J6).EQ.J5) THEN
                  IF ((MINUS(J6).EQ.PARENT(J5)).AND.(.NOT.SAVE_DEADTS(J6))) THEN 
                     IF (ETS(J6).LT.ETSMIN) THEN
                        ETSMIN=ETS(J6)
                        NTSMIN=J6
                     ENDIF
                  ENDIF
                  J6=POINTERP(J6)
               ELSE IF (MINUS(J6).EQ.J5) THEN
                  IF ((PLUS(J6).EQ.PARENT(J5)).AND.(.NOT.SAVE_DEADTS(J6))) THEN 
                     IF (ETS(J6).LT.ETSMIN) THEN
                        ETSMIN=ETS(J6)
                        NTSMIN=J6
                     ENDIF
                  ENDIF
                  J6=POINTERM(J6)
               ENDIF
               IF (J6.LT.0) EXIT pointa
            ENDDO pointa
   
            IF (NSTEPS.GT.MAXPATH) THEN
               PRINT '(A,I8)','ERROR in Dijkstra, path is longer than MAXPATH=',MAXPATH
               STOP
            ENDIF
   
            TSPATHID(NSTEPS)=NTSMIN
            IF (SAVE_DEADTS(NTSMIN)) THEN
               PRINT '(A,2I8,L5)','NTSMIN,TSATTEMPT,SAVE_DEADTS=',NTSMIN,TSATTEMPT(NTSMIN),SAVE_DEADTS(NTSMIN)
               PRINT '(A)','kshortestpaths> ERROR - this should never happen'
               STOP
            ENDIF
            MINPATHID(NSTEPS)=J5
!
!  End of TS identification.
!
            J5=PARENT(J5) 
            IF (J5.EQ.LJ1) EXIT
            IF (J5.EQ.0) EXIT
         ENDDO

         MINPATHID(NSTEPS+1)=LJ1

         IF (ALLOCATED(BESTPATH)) DEALLOCATE(BESTPATH)
         ALLOCATE(BESTPATH(2*NSTEPS+1))
         J5=LJ2
         DO J3=1,NSTEPS
            BESTPATH(2*J3-1)=J5         ! min
            ADDMIN(J5) = .TRUE.
            BESTPATH(2*J3)=TSPATHID(J3) ! ts
            ADDTS(TSPATHID(J3)) = .TRUE.
            IF (J3.EQ.NSTEPS) BESTPATH(2*J3+1)=PARENT(J5) ! start min
            J5=PARENT(J5)
         ENDDO
         ADDMIN(LJ1) = .TRUE.

!         IF (DEBUG) PRINT '(A,I6,A,A1,A,I6,A,G20.10,A,I8)', &
         PRINT '(A,I6,A,A1,A,I6,A,G20.10,A,I8)', &
 &         'kshortestpaths> kth best path for minimum ',LJ2,' and ', DIRECTION(2:2),' minimum ',LJ1,' k^SS=', &
 &            EXP(-WEIGHT(LJ2))*EXP(PFMIN(LJ1)-PFTOTALSTART)/EMKSUM(LJ1),' steps=',NSTEPS

!
! Calculate and return the forwards and backwards MFPTs for this unbranched kth fastest path, using GT.
!
         CALL MFPT(NSTEPS, TSPATHID, WAITAB, WAITBA)

         IF (DIRECTION.EQ.'AB') THEN ! B to A
            WRITE(*,'(A,G20.10)') 'kshortestpaths> detailed balance, ratio should be one if SS applies: ', &
   &                         WAITAB*EXP(PFMIN(LJ2)-PFMIN(LJ1))/WAITBA
            PRINT '(2(A,G20.10))','kshortestpaths> MFPT for this path: A<-B value=',WAITAB,' B<-A value=',WAITBA
            PRINT '(2(A,G20.10))','kshortestpaths> rates without conditional probability: A<-B value=', &
   &                        1.0D0/WAITAB,' B<-A value=',1.0D0/WAITBA
            PRINT '(2(A,G20.10))','kshortestpaths> rates with    conditional probability: A<-B value=', &
   &                        EXP(PFMIN(LJ1)-PFTOTALB)/WAITAB,' B<-A value=',EXP(PFMIN(LJ2)-PFTOTALA)/WAITBA
         ELSE ! A to B
            WRITE(*,'(A,G20.10)') 'kshortestpaths> detailed balance, ratio should be one if SS applies: ', &
   &                         WAITBA*EXP(PFMIN(LJ2)-PFMIN(LJ1))/WAITAB
            PRINT '(2(A,G20.10))','kshortestpaths> MFPT for this path: A<-B value=',WAITAB,' B<-A value=',WAITBA 
            PRINT '(2(A,G20.10))','kshortestpaths> rates without conditional probability: A<-B value=', &
   &                        1.0D0/WAITAB,' B<-A value=',1.0D0/WAITBA
            PRINT '(2(A,G20.10))','kshortestpaths> rates with    conditional probability: A<-B value=', &
   &                        EXP(PFMIN(LJ2)-PFTOTALB)/WAITAB,' B<-A value=',EXP(PFMIN(LJ1)-PFTOTALA)/WAITBA
         ENDIF
!
! Check that the product of DMAT entries gives the same path weight as we get from kshortestpaths
!
         IF (DUMMY.NE.0.0D0) THEN
            IF (DABS((EXP(-WEIGHT(LJ2))-DUMMY)/DUMMY).GT.1.0D-2) THEN
               PRINT '(A,2G20.10)','WARNING - alternative product of branching probabilities=',DUMMY, EXP(-WEIGHT(LJ2))
            ENDIF
         ENDIF
!
!  Note that VERYBEST contains the conditional probability and is divided by the waiting time.
!  We could get the contribution to NSS rate constants if we did some short KMC runs to get
!  the appropriate waiting time.
!
! We don't want DUMMY here but the reciprocal of the MFPT.
         IF (DIRECTION.EQ.'AB') THEN ! B to A
            WAIT = WAITAB
         ELSE
            WAIT = WAITBA
         ENDIF
         IF (EXP(PFMIN(LJ1)-PFTOTALSTART)/WAIT .GT. VERYBEST) THEN
            VERYBEST=EXP(PFMIN(LJ1)-PFTOTALSTART)/WAIT
            VERYBESTEND=LJ2
            VERYBESTSTART=LJ1
            VERYBESTPARENT(1:NMIN)=PARENT(1:NMIN)
         ENDIF
!
! Output printing for each fastest path
!
         IF (KSHORT_FULL_PRINTT) CALL KSHORT_PRINTING(NSTEPS, KSHORT, LJ1, LJ2, PARENT, TSPATHID)

      ENDDO ! end of k fastest paths loop for given starting minimum and endpoint
   ENDDO loopend! end of loop over endpoints
ENDDO loopstart ! end of loop over starting points
! 
!  Summarise the best path for any A(B) and any B(A)
!
PRINT *
PRINT '(2(A,A1),A)','kshortestpaths> Best path (chosen on GT rate constant weighted by conditional occupation &
                    &probability) between any ',DIRECTION(2:2),' minimum and any ',DIRECTION(1:1),' minimum:'
J5=VERYBESTEND
WRITE(*,'(I8)',ADVANCE='NO') J5
CALL FLUSH(6)
NSTEPS=0
DO
  NSTEPS=NSTEPS+1
!
!  Find the transition state with the lowest barrier between minima
!  J5 and VERYBESTPARENT(J5).
!
  J1=TOPPOINTER(J5) ! sets J1 to the transition state connected to J5 with the
                    ! highest number
  ETSMIN=1.0D100
  point: DO
     IF (PLUS(J1).EQ.J5) THEN
        IF (MINUS(J1).EQ.VERYBESTPARENT(J5)) THEN
           IF (ETS(J1).LT.ETSMIN) THEN
              ETSMIN=ETS(J1)
              NTSMIN=J1
           ENDIF
        ENDIF
        J1=POINTERP(J1)
     ELSE IF (MINUS(J1).EQ.J5) THEN
        IF (PLUS(J1).EQ.VERYBESTPARENT(J5)) THEN
           IF (ETS(J1).LT.ETSMIN) THEN
              ETSMIN=ETS(J1)
              NTSMIN=J1
           ENDIF
        ENDIF
        J1=POINTERM(J1)
     ENDIF
    IF (J1.LT.0) EXIT point
  ENDDO point
  TSPATHID(NSTEPS)=NTSMIN

  J5=VERYBESTPARENT(J5)
  WRITE(*,'(I8)',ADVANCE='NO') J5
  IF (J5.EQ.VERYBESTSTART) EXIT
  IF (J5.EQ.0) EXIT
ENDDO
PRINT*
PRINT '(2(A,A1),A,G20.10,A,I6,A)','kshortestpaths> Largest contribution to rate constant ', &
&  DIRECTION(1:1),'<-',DIRECTION(2:2), &
&       ' for any A and B is ',VERYBEST,' for ',NSTEPS,' transition states'

!
!  Dump potential energy stationary points corresponding to fastest path. Must allow for
!  regrouping and sorting of free energy groups using the MINGROUP and MINMAP arrays.
!  MINGROUP(J1) is the original free energy group to which potential energy minimum J1
!  belongs. MINMAP(J2) tells us which original free energy group the sorted free energy
!  group J2 corresponds to. This is enough information!
!
   IF (REGROUPFREET.OR.REGROUPFREEABT.OR.REGROUPRATET.OR.REGROUPPET.OR.REGROUPPERSISTT) THEN ! we have to deal with free energy groups
      OPEN(UNIT=2,FILE='min.data.fastest',STATUS='UNKNOWN')
      OPEN(UNIT=4,FILE='points.min.fastest',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
      IF (ALLOCATED(INCLUDEMIN)) DEALLOCATE(INCLUDEMIN)
      IF (ALLOCATED(NEWINDEX)) DEALLOCATE(NEWINDEX)
      ALLOCATE(INCLUDEMIN(NMINSAVE),NEWINDEX(NMINSAVE))
      INCLUDEMIN(1:NMINSAVE)=.FALSE.
      NEWINDEX(1:NMINSAVE)=0
      DO J1=1,NSTEPS+1 ! these are the free energy minima
         J5=MINMAP(BESTPATH(2*J1-1))
         DO J2=1,NMINSAVE
            IF (MINGROUP(J2).EQ.J5) THEN
               IF (.NOT.INCLUDEMIN(J2)) THEN
                  INCLUDEMIN(J2)=.TRUE.
               ENDIF
               IF (DEBUG) PRINT '(4(A,I8))','Dijkstra> Including PE minimum ',J2,' from group ',BESTPATH(2*J1-1), &
  &                             ' originally ',J5
            ENDIF
         ENDDO
      ENDDO
      OPEN(UNIT=7,FILE='min.data',ACTION='READ',STATUS='OLD')
      NDUMMY=0
      DO J2=1,NMINSAVE
         READ(7,'(A)') DUMMYSTRING
         IF (.NOT.INCLUDEMIN(J2)) CYCLE
         NDUMMY=NDUMMY+1
         NEWINDEX(J2)=NDUMMY
         WRITE(2,'(A,2X,I8)') TRIM(ADJUSTL(DUMMYSTRING)),J2
         READ(UMIN,REC=J2) (LOCALPOINTS(J3),J3=1,3*NATOMS)
         WRITE(4,REC=NDUMMY) (LOCALPOINTS(J3),J3=1,3*NATOMS)
      ENDDO
      CLOSE(2); CLOSE(4); CLOSE(7)
      OPEN(UNIT=7,FILE='ts.data',ACTION='READ',STATUS='OLD')
      OPEN(UNIT=3,FILE='ts.data.fastest',STATUS='UNKNOWN')
      OPEN(UNIT=5,FILE='points.ts.fastest',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
      NDUMMY=0
      NDUMMY2=0
      DO 
         NDUMMY=NDUMMY+1
         IF (IMFRQT) THEN
            READ(7,*,END=40) DETS,DFVIBTS,DHORDERTS,DPLUS,DMINUS,DIXTS,DIYTS,DIZTS,DNEGEIG
         ELSE
            READ(7,*,END=40) DETS,DFVIBTS,DHORDERTS,DPLUS,DMINUS,DIXTS,DIYTS,DIZTS
         ENDIF
         IF (INCLUDEMIN(DPLUS).AND.INCLUDEMIN(DMINUS)) THEN
            IF (DEBUG) PRINT '(3(A,I8))','Dijkstra> Including PE ts ',NDUMMY
            NDUMMY2=NDUMMY2+1
            IF (IMFRQT) THEN
               WRITE(3,'(2F20.10,3I8,4F20.10,2X,I8)')  &
  &                DETS,DFVIBTS,DHORDERTS,NEWINDEX(DPLUS),NEWINDEX(DMINUS),DIXTS,DIYTS,DIZTS,DNEGEIG,NDUMMY
            ELSE
               WRITE(3,'(2F20.10,3I8,3F20.10,2X,I8)')  &
  &                DETS,DFVIBTS,DHORDERTS,NEWINDEX(DPLUS),NEWINDEX(DMINUS),DIXTS,DIYTS,DIZTS,NDUMMY
            ENDIF
            IF (.NOT.DUMMYTST) THEN
               READ(UTS,REC=NDUMMY) (LOCALPOINTS(J2),J2=1,3*NATOMS)
               WRITE(5,REC=NDUMMY2) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            ENDIF
         ENDIF
      ENDDO
40    CONTINUE
      CLOSE(3); CLOSE(5); CLOSE(7)
      OPEN(UNIT=2,FILE='min.A.fastest',STATUS='UNKNOWN')
      OPEN(UNIT=3,FILE='min.A',ACTION='READ',STATUS='OLD')
      READ(3,*) NNMINA
      NDUMMY2=0
      DO J1=1,NNMINA
         READ(3,*) NDUMMY
         IF (INCLUDEMIN(NDUMMY)) NDUMMY2=NDUMMY2+1
      ENDDO
      WRITE(2,*) NDUMMY2
      REWIND(3)
      READ(3,*) NNMINA
      DO J1=1,NNMINA
         READ(3,*) NDUMMY
         IF (INCLUDEMIN(NDUMMY)) WRITE(2,*) NEWINDEX(NDUMMY)
      ENDDO
      CLOSE(2); CLOSE(3)

      OPEN(UNIT=4,FILE='min.B.fastest',STATUS='UNKNOWN')
      OPEN(UNIT=5,FILE='min.B',ACTION='READ',STATUS='OLD')
      READ(5,*) NNMINB
      NDUMMY2=0
      DO J1=1,NNMINB
         READ(5,*) NDUMMY
         IF (INCLUDEMIN(NDUMMY)) NDUMMY2=NDUMMY2+1
      ENDDO
      WRITE(4,*) NDUMMY2
      REWIND(5)
      READ(5,*) NNMINB
      DO J1=1,NNMINB
         READ(5,*) NDUMMY
         IF (INCLUDEMIN(NDUMMY)) WRITE(4,*) NEWINDEX(NDUMMY)
      ENDDO
      CLOSE(4); CLOSE(5)
      DEALLOCATE(INCLUDEMIN,NEWINDEX)
   ELSE ! the PE minima were not regrouped
      OPEN(UNIT=2,FILE='min.data.fastest',STATUS='UNKNOWN')
      OPEN(UNIT=4,FILE='points.min.fastest',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)

      J5 = VERYBESTEND
      WRITE(2,'(2F20.10,I6,3F20.10)') EMIN(J5), FVIBMIN(J5), HORDERMIN(J5), IXMIN(J5), IYMIN(J5), IZMIN(J5)
      DO J1=1,NSTEPS
         J5=VERYBESTPARENT(J5)
         WRITE(2,'(2F20.10,I6,3F20.10)') EMIN(J5), FVIBMIN(J5), HORDERMIN(J5), IXMIN(J5), IYMIN(J5), IZMIN(J5)
         READ(UMIN,REC=J5) (LOCALPOINTS(J2),J2=1,3*NATOMS)
         WRITE(4,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
      ENDDO
      CLOSE(2); CLOSE(4)
      OPEN(UNIT=3,FILE='ts.data.fastest',STATUS='UNKNOWN')
      OPEN(UNIT=5,FILE='points.ts.fastest',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
      DO J1=1,NSTEPS
         J5=TSPATHID(J1)
         IF (IMFRQT) THEN
            WRITE(3,'(2F20.10,3I10,4F20.10)') ETS(J5),FVIBTS(J5),HORDERTS(J5),J1,J1+1,IXTS(J5),IYTS(J5),IZTS(J5),NEGEIG(J5)
         ELSE
            WRITE(3,'(2F20.10,3I10,3F20.10)') ETS(J5),FVIBTS(J5),HORDERTS(J5),J1,J1+1,IXTS(J5),IYTS(J5),IZTS(J5)
         ENDIF
         IF (.NOT.DUMMYTST) THEN
            READ(UTS,REC=J5) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            WRITE(5,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
         ENDIF
      ENDDO
      CLOSE(3); CLOSE(5)
      OPEN(UNIT=2,FILE='min.A.fastest',STATUS='UNKNOWN')
      OPEN(UNIT=4,FILE='min.B.fastest',STATUS='UNKNOWN')
      WRITE(2,'(I6)') 1
      WRITE(4,'(I6)') 1
      IF (DIRECTION.EQ.'AB') THEN
!
! jmc swapped units - check!!!
!
         WRITE(2,'(I6)') 1
         WRITE(4,'(I6)') NSTEPS+1
      ELSE
         WRITE(4,'(I6)') 1
         WRITE(2,'(I6)') NSTEPS+1
      ENDIF
      CLOSE(2); CLOSE(4)
   ENDIF

PRINT '(A)','kshortestpaths> Writing [min.data,points.min,ts.data,points.ts].fastest.all for all the kth fastest paths found:'

OPEN(UNIT=2,FILE='min.data.fastest.all',STATUS='UNKNOWN')
OPEN(UNIT=4,FILE='points.min.fastest.all',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
COUNTER = 0
DO J5=1,NMIN
   IF (ADDMIN(J5)) THEN
      COUNTER = COUNTER + 1
      WRITE(2,'(2F20.10,I6,3F20.10)') EMIN(J5), FVIBMIN(J5), HORDERMIN(J5), IXMIN(J5), IYMIN(J5), IZMIN(J5)
      READ(UMIN,REC=J5) (LOCALPOINTS(J2),J2=1,3*NATOMS)
      WRITE(4,REC=COUNTER) (LOCALPOINTS(J2),J2=1,3*NATOMS)
      INDEX_NEW_MIN(J5) = COUNTER
   ENDIF
ENDDO
CLOSE(2); CLOSE(4)
OPEN(UNIT=3,FILE='ts.data.fastest.all',STATUS='UNKNOWN')
OPEN(UNIT=5,FILE='points.ts.fastest.all',ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='UNKNOWN',RECL=8*3*NATOMS)
COUNTER = 0
DO J5=1,NTS
   IF (ADDTS(J5)) THEN
      COUNTER = COUNTER + 1
      IF (IMFRQT) THEN
         WRITE(3,'(2F20.10,3I10,3F20.10)') ETS(J5),FVIBTS(J5),HORDERTS(J5),INDEX_NEW_MIN(PLUS(J5)),INDEX_NEW_MIN(MINUS(J5)), &
  &                                        IXTS(J5),IYTS(J5),IZTS(J5),NEGEIG(J5)
      ELSE
         WRITE(3,'(2F20.10,3I10,3F20.10)') ETS(J5),FVIBTS(J5),HORDERTS(J5),INDEX_NEW_MIN(PLUS(J5)),INDEX_NEW_MIN(MINUS(J5)), &
  &                                        IXTS(J5),IYTS(J5),IZTS(J5)
      ENDIF
      IF (.NOT.DUMMYTST) THEN
         READ(UTS,REC=J5) (LOCALPOINTS(J2),J2=1,3*NATOMS)
         WRITE(5,REC=COUNTER) (LOCALPOINTS(J2),J2=1,3*NATOMS)
      ENDIF
   ENDIF
ENDDO
CLOSE(3); CLOSE(5)

OPEN(UNIT=2,FILE='min.A.fastest.all',STATUS='UNKNOWN')
OPEN(UNIT=4,FILE='min.B.fastest.all',STATUS='UNKNOWN')
COUNTER = 0
DO J5 = 1, NMINA
   IF (ADDMIN(LOCATIONA(J5))) COUNTER = COUNTER + 1
ENDDO
WRITE(2,'(I6)') COUNTER
DO J5 = 1, NMINA
   IF (ADDMIN(LOCATIONA(J5))) WRITE(2,'(I6)') INDEX_NEW_MIN(LOCATIONA(J5))
ENDDO
COUNTER = 0
DO J5 = 1, NMINB
   IF (ADDMIN(LOCATIONB(J5))) COUNTER = COUNTER + 1
ENDDO
WRITE(4,'(I6)') COUNTER
DO J5 = 1, NMINB
   IF (ADDMIN(LOCATIONB(J5))) WRITE(4,'(I6)') INDEX_NEW_MIN(LOCATIONB(J5))
ENDDO
CLOSE(2); CLOSE(4)

CALL CPU_TIME(TNEW)
TKSHORTESTPATHS = TKSHORTESTPATHS + TNEW - ELAPSED

RETURN
END SUBROUTINE KSHORTESTPATHS


SUBROUTINE KSHORT_PRINTING(NSTEPS, KSHORT, START, END, PARENT, TSPATHID)
USE PORFUNCS
USE COMMONS
IMPLICIT NONE

INTEGER, PARAMETER :: MAXPATH=10000 ! longest path length allowed - surely 10000 is enough?!
INTEGER :: NSTEPS, J5, KSHORT, END, J1, PARENT(NMIN), TSPATHID(MAXPATH), J2, NDUMMY, START
INTEGER IRES(NATOMS)

DOUBLE PRECISION ALIGNPOINTS(3*NATOMS), DIST2, DISTANCES, RMAT(3,3), LOCALPOINTS(3*NATOMS), DUMMY

CHARACTER(LEN=4) SEGID(NATOMS), TYPE(NATOMS), RES(NATOMS)
CHARACTER(LEN=10) CONNSTR, SSTR, ESTR
CHARACTER(LEN=80) FPOO, FPOO2
!
! File name indices are x.start.end.k
!
IF (DIRECTION.EQ.'AB') PRINT '(A)','kshortestpaths> Note that path is printed backwards starting with A, ending with B'
IF (DIRECTION.EQ.'BA') PRINT '(A)','kshortestpaths> Note that path is printed backwards starting with B, ending with A'

WRITE(CONNSTR,'(I10)') KSHORT
WRITE(SSTR,'(I10)') START
WRITE(ESTR,'(I10)') END
FPOO='Epath.fastest.'//TRIM(ADJUSTL(SSTR))//'.'//TRIM(ADJUSTL(ESTR))//'.'//TRIM(ADJUSTL(CONNSTR)) ! workaround for Sun compiler bug
OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FPOO)),STATUS='UNKNOWN')
J5=END
DO J1=1,NSTEPS
   WRITE(1,'(I8,G20.10,I8)') 2*J1-1,EMIN(J5),J5
   WRITE(1,'(I8,G20.10,I8)') 2*J1,ETS(TSPATHID(J1)),TSPATHID(J1)
   IF (J1.EQ.NSTEPS) WRITE(1,'(I8,G20.10,I8)') 2*NSTEPS+1,EMIN(PARENT(J5)),PARENT(J5)
   J5=PARENT(J5)
ENDDO
CLOSE(1)

IF (.NOT.(REGROUPRATET.OR.REGROUPPET.OR.REGROUPFREET.OR.DUMMYTST .OR. REGROUPPERSISTT)) PRINT '(A,I6)', &
&       'kshortestpaths> dumping minima and transition state coordinates to file redopoints.fastest.xxx, steps=',NSTEPS

FPOO='stationary.points.pdb.fastest.'//TRIM(ADJUSTL(SSTR))//'.'//TRIM(ADJUSTL(ESTR))//'.'//TRIM(ADJUSTL(CONNSTR))
IF (CHARMMT) THEN
   OPEN(UNIT=4,FILE='input.crd',STATUS='OLD')
   READ(4,*) NDUMMY
   DO J1=1,NATOMS
      READ(4,*) NDUMMY,IRES(J1),RES(J1),TYPE(J1),DUMMY,DUMMY,DUMMY,SEGID(J1)
   ENDDO
   CLOSE(4)
   IF (.NOT.(REGROUPRATET.OR.REGROUPPET.OR.REGROUPFREET.OR.DUMMYTST .OR. REGROUPPERSISTT)) OPEN(UNIT=3,FILE=FPOO,STATUS='UNKNOWN')
!
!  Use our free format .crd file to get the necessary data to make a fixed format pdb file
!
   J5=END
   DO J1=1,NSTEPS
      IF (.NOT.(REGROUPRATET.OR.REGROUPPET.OR.REGROUPFREET.OR.DUMMYTST .OR. REGROUPPERSISTT)) THEN
         READ(UMIN,REC=J5) (LOCALPOINTS(J2),J2=1,3*NATOMS)
         IF (J1.GT.1) THEN
            CALL MINPERMDIST(ALIGNPOINTS,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCES,DIST2,RIGIDBODY, &
  &                          RMAT,.FALSE.)
         ENDIF
         DO J2=1,NATOMS
!           WRITE(3,'(a4,2x,i5,2x,a4,a4,1x,i4,4x,3f8.3,2f6.2,6x,a4)') &
            WRITE(3,'(a4,2x,i5,1x,a4,1x,a4,1x,i4,4x,3f8.3,2f6.2,6x,a4)') &
&           'ATOM',J2,TYPE(J2),RES(J2),IRES(J2),LOCALPOINTS(3*(J2-1)+1),LOCALPOINTS(3*(J2-1)+2),LOCALPOINTS(3*(J2-1)+3), &
&           1.00,0.00,SEGID(J2)
         ENDDO
         WRITE(3,'(A)') 'END'
         ALIGNPOINTS(1:3*NATOMS)=LOCALPOINTS(1:3*NATOMS)
         READ(UTS,REC=TSPATHID(J1)) (LOCALPOINTS(J2),J2=1,3*NATOMS)
         CALL MINPERMDIST(ALIGNPOINTS,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCES,DIST2,RIGIDBODY, &
  &                       RMAT,.FALSE.)
         DO J2=1,NATOMS
!           WRITE(3,'(a4,2x,i5,2x,a4,a4,1x,i4,4x,3f8.3,2f6.2,6x,a4)') &
            WRITE(3,'(a4,2x,i5,1x,a4,1x,a4,1x,i4,4x,3f8.3,2f6.2,6x,a4)') &
&           'ATOM',J2,TYPE(J2),RES(J2),IRES(J2),LOCALPOINTS(3*(J2-1)+1),LOCALPOINTS(3*(J2-1)+2),LOCALPOINTS(3*(J2-1)+3), &
&            1.00,0.00,SEGID(J2)
         ENDDO
         WRITE(3,'(A)') 'END'
         ALIGNPOINTS(1:3*NATOMS)=LOCALPOINTS(1:3*NATOMS)
         IF (J1.EQ.NSTEPS) THEN
            READ(UMIN,REC=PARENT(J5)) (LOCALPOINTS(J2),J2=1,3*NATOMS)
            CALL MINPERMDIST(ALIGNPOINTS,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCES,DIST2,RIGIDBODY, &
  &                          RMAT,.FALSE.)
            DO J2=1,NATOMS
               WRITE(3,'(a4,2x,i5,2x,a4,a4,1x,i4,4x,3f8.3,2f6.2,6x,a4)') &
&              'ATOM',J2,TYPE(J2),RES(J2),IRES(J2),LOCALPOINTS(3*(J2-1)+1),LOCALPOINTS(3*(J2-1)+2),LOCALPOINTS(3*(J2-1)+3), &
&              1.00,0.00,SEGID(J2)
            ENDDO
            WRITE(3,'(A)') 'END'
         ENDIF
      ENDIF
      J5=PARENT(J5)
   ENDDO
   IF (.NOT.(REGROUPRATET.OR.REGROUPPET.OR.REGROUPFREET.OR.DUMMYTST .OR. REGROUPPERSISTT)) CLOSE(3)
ENDIF

IF (.NOT.(REGROUPRATET.OR.REGROUPPET.OR.REGROUPFREET.OR.DUMMYTST .OR. REGROUPPERSISTT)) THEN
   FPOO2='redopoints.fastest.'//TRIM(ADJUSTL(SSTR))//'.'//TRIM(ADJUSTL(ESTR))//'.'//TRIM(ADJUSTL(CONNSTR))
   OPEN(UNIT=1,FILE=FPOO2,STATUS='UNKNOWN')
   J5=END
   DO J1=1,NSTEPS
      READ(UMIN,REC=J5) (LOCALPOINTS(J2),J2=1,3*NATOMS)
      IF (J1.GT.1) THEN
         CALL MINPERMDIST(ALIGNPOINTS,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCES,DIST2,RIGIDBODY, &
  &                       RMAT,.FALSE.)
      ENDIF
!
! We MUST use fixed format here in case the coordinates are read back into CHARMM.
! CHARMM cannot read exponential format, even if the free keyword is specified!!
!
      WRITE(1,'(3F25.15)') (LOCALPOINTS(J2),J2=1,3*NATOMS)
      ALIGNPOINTS(1:3*NATOMS)=LOCALPOINTS(1:3*NATOMS)
      READ(UTS,REC=TSPATHID(J1)) (LOCALPOINTS(J2),J2=1,3*NATOMS)
      CALL MINPERMDIST(ALIGNPOINTS,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCES,DIST2,RIGIDBODY, &
  &                    RMAT,.FALSE.)
      WRITE(1,'(3F25.15)') (LOCALPOINTS(J2),J2=1,3*NATOMS)
      ALIGNPOINTS(1:3*NATOMS)=LOCALPOINTS(1:3*NATOMS)
      IF (J1.EQ.NSTEPS) THEN ! last minimum
!        READ(UMIN,REC=VERYBESTPARENT(J5)) (LOCALPOINTS(J2),J2=1,3*NATOMS)
         READ(UMIN,REC=PARENT(J5)) (LOCALPOINTS(J2),J2=1,3*NATOMS)
         CALL MINPERMDIST(ALIGNPOINTS,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCES,DIST2,RIGIDBODY, &
  &                       RMAT,.FALSE.)
         WRITE(1,'(3F25.15)') (LOCALPOINTS(J2),J2=1,3*NATOMS)
      ENDIF
      J5=PARENT(J5)
   ENDDO
   CLOSE(1)
ENDIF

RETURN
END SUBROUTINE KSHORT_PRINTING


!-------------------------------------------------------------------------------

SUBROUTINE my_dijkstra(start, nmin)
USE GRAPH
IMPLICIT NONE

INTEGER :: start, i, j, minadj, n, m, nmin
DOUBLE PRECISION :: minweight
TYPE(edge), POINTER :: p

! Initialization
    DO i = 1, nmin
       shortest_paths(i,1)%l = HUGE(0.0D0)
       shortest_paths(i,1)%k = 1
       IF (minnode(i)%index == 0) CYCLE ! i.e. index is used to flag an unphysical minimum.
       minnode(i)%permanent = .FALSE.
       shortest_paths(i,1)%node => minnode(i)
    ENDDO

! Start algorithm at 'start' minimum
    i = start 
    minnode(i)%permanent = .TRUE.
    shortest_paths(i,1)%l = 0.0D0
  
    DO n = 1, nmin                         ! Iterate nmin times
 
       ! Loop over all minima connected to i
       p => minnode(i)%succ ! p points to an edge
       DO
          IF (.NOT. ASSOCIATED(p)) EXIT
          m = p%to%index                   ! m is connected minimum
          IF (minnode(m)%permanent) THEN   ! Only consider temporary minima
             p => p%nextsucc
             CYCLE
          ENDIF
 
          ! Update current estimate of shortest-path weight for each connected minimum
          IF ((shortest_paths(i,1)%l + p%l) < shortest_paths(m,1)%l) THEN
             shortest_paths(m,1)%l = shortest_paths(i,1)%l + p%l
             shortest_paths(m,1)%pred => shortest_paths(i,1)
          ENDIF
          p => p%nextsucc
       ENDDO

       ! Find temporary minimum with lowest current shortest-path weight
       minweight = HUGE(minweight)
       minadj = 0
       DO j = 1, nmin
          IF (minnode(j)%permanent .OR. (minnode(j)%index == 0)) CYCLE ! Only consider temporary, connected minima
          IF (shortest_paths(j,1)%l < minweight) THEN
             minweight = shortest_paths(j,1)%l
             minadj = j
          ENDIF
       ENDDO
       IF (minadj == 0) EXIT
       i = minadj
       minnode(i)%permanent = .TRUE.

    ENDDO
    
    DO i = 1, nmin
       IF (minnode(i)%index == 0) CYCLE
       minnode(i)%shortest => shortest_paths(i,1)
    ENDDO
 
END SUBROUTINE my_dijkstra


SUBROUTINE add_pred(i, j)
USE GRAPH
IMPLICIT NONE

INTEGER :: i, j

IF (ASSOCIATED(minnode(i)%pred)) THEN
   tsedge(j)%nextpred => minnode(i)%pred
   minnode(i)%pred => tsedge(j)
ELSE
   minnode(i)%pred => tsedge(j)
ENDIF

END SUBROUTINE add_pred


SUBROUTINE add_succ(i, j)
USE GRAPH
IMPLICIT NONE

INTEGER :: i, j

IF (ASSOCIATED(minnode(i)%succ)) THEN
   tsedge(j)%nextsucc => minnode(i)%succ
   minnode(i)%succ => tsedge(j)
ELSE
   minnode(i)%succ => tsedge(j)
ENDIF

END SUBROUTINE add_succ


SUBROUTINE add_candidate(i, j, k, weight) ! (node, predecessor, which path the predecessor is from, weight of j -> i)
USE GRAPH
IMPLICIT NONE

TYPE(path), POINTER :: p

INTEGER :: i, j, k, m, n1, n2
DOUBLE PRECISION :: weight

m = index_next(i)

! First add the j -> i connection to the kth shortest path to node j.
! The resulting path is at position m (=index_next(i)) in the candidate_paths array.
candidate_paths(m)%l = shortest_paths(j,k)%l + weight
candidate_paths(m)%pred => shortest_paths(j,k)
candidate_paths(m)%node => minnode(i)
candidate_paths(m)%i = m ! use %i to store the index for this path in candidate_paths
candidate_paths(m)%k = m - index_start(i) + 1

! Arrange the linked list of candidates in order of decreasing weight, i.e. tail of the list has the smallest weight.
IF (ASSOCIATED(minnode(i)%cand)) THEN ! there is at least one entry in the list already.
   NULLIFY(p)
   p => minnode(i)%cand
   IF (.NOT. ASSOCIATED(p%next)) THEN ! there is only one item in the list.
      !n1 = p%k
      n1 = p%i
      IF (p%l >= candidate_paths(m)%l) THEN ! add new path to the end of the list.
         candidate_paths(n1)%next => candidate_paths(m)
      ELSE ! add new path to the start of the list.
         candidate_paths(m)%next => candidate_paths(n1)
         minnode(i)%cand => candidate_paths(m)
      ENDIF
   ELSE
      IF (candidate_paths(m)%l >= p%l) THEN ! add the new path to the top of the list
         n1 = p%i
         candidate_paths(m)%next => candidate_paths(n1)
         minnode(i)%cand => candidate_paths(m)
      ELSE
         DO
            IF (.NOT. ASSOCIATED(p%next)) THEN ! this item is the last in the list, so add the new path after it.
               !n1 = p%k
               n1 = p%i
               candidate_paths(n1)%next => candidate_paths(m)
               EXIT
            ELSE ! Insert this candidate between the current item and the next.
               IF ((p%next%l < candidate_paths(m)%l) .AND. (p%l >= candidate_paths(m)%l)) THEN
                  !n1 = p%next%k
                  !n2 = p%k
                  n1 = p%next%i
                  n2 = p%i
                  candidate_paths(m)%next => candidate_paths(n1)
                  candidate_paths(n2)%next => candidate_paths(m)
                  EXIT
               ENDIF
            ENDIF
            p => p%next
         ENDDO
      ENDIF
   ENDIF
ELSE
   minnode(i)%cand => candidate_paths(m)
ENDIF

index_next(i) = index_next(i) + 1
!print *,'testing: at end of add_candidate, index_next(i), i = ',index_next(i),i

END SUBROUTINE add_candidate


SUBROUTINE delete_candidate(j, i, k) ! j is the index of the node; i is the index of the candidate to be removed; 
                                     ! k is the index of the predecessor of path i in the list.
USE GRAPH
IMPLICIT NONE

INTEGER :: i, j, k

! We're always removing the tail entry of the linked list, but 
! we need to be careful with its index in the candidate_paths array.
IF (.NOT. ASSOCIATED(minnode(j)%cand%next)) THEN
   NULLIFY(minnode(j)%cand)
ELSE
   NULLIFY(candidate_paths(k)%next)
ENDIF
candidate_paths(i)%l = 0.0D0
candidate_paths(i)%k = 0
candidate_paths(i)%i = 0
NULLIFY(candidate_paths(i)%pred,candidate_paths(i)%next,candidate_paths(i)%node)

index_next(j) = i
!print *,'testing: at end of delete_candidate, index_next(j), j = ',index_next(j),j

END SUBROUTINE delete_candidate


RECURSIVE SUBROUTINE next_path(v, k, start)
USE GRAPH
IMPLICIT NONE

INTEGER :: v, k, start, u, k1, w
DOUBLE PRECISION :: weight

TYPE(edge), POINTER :: p
TYPE(path), POINTER :: pp, newpp, n1, n2
NULLIFY(p, pp, newpp)

! For situations where the source has no incoming edges (e.g. analysis of DPS databases)
IF (v == start .AND. (.NOT. ASSOCIATED(minnode(v)%pred))) RETURN

! Initialization: add entries to the candidate list for node v using its predecessors 
! and the shortest shortest paths from start -> pred, excluding the pred of v 
! that occurs in its shortest shortest path.
IF (k == 2) THEN
   p => minnode(v)%pred ! There will be at least one incoming connection to every node that is considered at this point.
   DO
      IF (.NOT. ASSOCIATED(p)) EXIT
      u = p%from%index
      IF (shortest_paths(v,1)%pred%node%index /= u) THEN
         weight = p%l
         CALL add_candidate(v, u, 1, weight)
      ENDIF
      p => p%nextpred
   ENDDO

   IF (v == start) GOTO 300
ENDIF

! Find the predecessor of v in the k-1 th shortest path to v.
u = shortest_paths(v,k-1)%pred%node%index

! Get the weight of the u -> v connection
NULLIFY(p)
p => minnode(v)%pred
DO
   IF (p%from%index == u) THEN
      weight = p%l ! needed in the call to add_candidate
      EXIT
   ENDIF
   p => p%nextpred
ENDDO

newpp => minnode(u)%shortest
outer: DO
   IF (.NOT. ASSOCIATED(newpp)) THEN
      PRINT *,'Cannot find the kprime th shortest path for node ',u
      STOP
   ENDIF
   NULLIFY(n1, n2)
   n1 => newpp ! giving the nodes in the set of shortest paths to u
   n2 => shortest_paths(v, k-1)%pred ! giving the nodes in the k-1 th shortest path to v, stopping at u.
   inner: DO
      IF (.NOT. (ASSOCIATED(n1) .AND. ASSOCIATED(n2))) EXIT inner
      IF (n1%node%index /= n2%node%index) THEN
         newpp => newpp%next
         CYCLE outer
      ENDIF
      n1 => n1%pred
      n2 => n2%pred
   ENDDO inner
   k1 = newpp%k
   EXIT outer
ENDDO outer

! If the k1 + 1 th shortest path from node u does not exist, call next_path(u,k1+1) to try and find it.
IF (.NOT. ASSOCIATED(shortest_paths(u,k1 + 1)%node)) THEN
   !PRINT *,'Calling next_path to find the k1+1 th shortest path for node ',u,k1 + 1
   CALL next_path(u, k1 + 1, start)
ENDIF

! If it now exists, add the k1 + 1 th shortest path to u, plus the u-v edge, to the list of candidates for node v...
IF (ASSOCIATED(shortest_paths(u,k1 + 1)%node)) THEN
   CALL add_candidate(v, u, k1 + 1, weight) ! (node, predecessor, which path the predecessor is from, weight)
ENDIF

300 CONTINUE

! If we have some candidates for the kth shortest path to node v, then 
! add the one with the smallest weight to the top of the list of shortest paths.
! That path will be at the bottom of the candidate list, and must then be removed 
! from that list.  The entry in the candidate_paths array is also deleted, ready to 
! be refilled with the next item added to the candidate list.
NULLIFY(pp)
IF (ASSOCIATED(minnode(v)%cand)) THEN
   ! Find the bottom of the candidate linked list.
   pp => minnode(v)%cand
   DO
      !PRINT *,'candidate list before ',v,pp%l,pp%k,pp%i,pp%pred%node%index
      IF (.NOT. ASSOCIATED(pp%next)) EXIT
      !w = pp%k
      w = pp%i
      pp => pp%next
   ENDDO
   !u = pp%k
   u = pp%i
   shortest_paths(v, k) = candidate_paths(u)
   shortest_paths(v, k)%k = k
   shortest_paths(v, k)%next => minnode(v)%shortest
   minnode(v)%shortest => shortest_paths(v, k)
   CALL delete_candidate(v, u, w)
! testing...
!   NULLIFY(pp)
!   pp => minnode(v)%cand
!   DO
!      IF (.NOT. ASSOCIATED(pp)) EXIT
!      PRINT *,'candidate list after ',v,pp%l,pp%i,pp%pred%node%index
!      pp => pp%next
!   ENDDO
! end testing...
ELSE
   PRINT *,'Shortest path ',k,' for node ',v,' does not exist.'
   RETURN
ENDIF

END SUBROUTINE next_path


RECURSIVE SUBROUTINE print_path(p)
USE GRAPH, ONLY : path, edge, node, start_shift_factor
IMPLICIT NONE

INTERFACE
   RECURSIVE SUBROUTINE print_nodes_in_path(pj)
   USE GRAPH, ONLY : path, edge, node
   IMPLICIT NONE
   TYPE(path), POINTER :: pj
   END SUBROUTINE print_nodes_in_path
END INTERFACE

TYPE(path), POINTER :: p

IF (ASSOCIATED(p)) then
   CALL print_path(p%next)
   !PRINT*,'Shortest path ',p%k,' length ',p%l,' rate constant ',EXP(-p%l)
   PRINT*,'Shortest path ',p%k,' length ',p%l - start_shift_factor,' rate constant ',EXP(-(p%l - start_shift_factor))
   CALL print_nodes_in_path(p)
ENDIF

END SUBROUTINE print_path


RECURSIVE SUBROUTINE print_nodes_in_path(pj)
USE GRAPH, ONLY : path, edge, node
IMPLICIT NONE

TYPE(path), POINTER :: pj

IF (ASSOCIATED(pj)) then
   CALL print_nodes_in_path(pj%pred)
   PRINT *,pj%node%index
ENDIF

END SUBROUTINE print_nodes_in_path

!
!  Construct DMAT0
!
!  NCOL(M2)     = # connections for minimum M2
!  NVAL(J1,M2)  = index of minimum involved in connection J1 from minimum M2
!  INDEX_TS(J1,M2) = index of TS involved in connection J1 from minimum M2
!  DMAT0(J2,x)  = KMC-type probability of taking TS J2 in the forwards (J2,1) and backwards (J2,2) directions.
!
!  Degenerate rearrangements are excluded.
!
      SUBROUTINE MAKED4(DMAT0,NCOL,NVAL,DEADTS,KSUM,INDEX_TS)
      USE COMMONS
      IMPLICIT NONE

      INTEGER :: NCOL(NMIN), J1, M1, M2, NVAL(NCONNMAX,NMIN), OTHER, INDEX_TS(NCONNMAX,NMIN)
      DOUBLE PRECISION :: DMAT0(NTS,2), KSUM(NMIN)
      LOGICAL :: DEADTS(NTS), MATCHED, ALLOWED, DTS(NTS)
!
!  Cycle over connected transition states using pointers for minimum M2.
!
      DMAT0 = 0.0D0
      NCOL = 0
      NVAL = 0
      INDEX_TS = 0
      DTS = DEADTS

      fromloop: DO M2=1,NMIN   
         J1=TOPPOINTER(M2)  !  sets J1 to the ts connected to minimum M2 with the highest value
         IF (J1.LE.0) CYCLE fromloop
11       CONTINUE
         IF ((.NOT.DEADTS(J1)).AND.(PLUS(J1).NE.MINUS(J1))) THEN
            IF (PLUS(J1).EQ.M2) THEN
               OTHER=MINUS(J1)
            ELSE 
               OTHER=PLUS(J1)
            ENDIF
            MATCHED=.FALSE.
            ALLOWED=.TRUE.
            IF (ALLOWED) THEN
               IF (PLUS(J1).EQ.M2) THEN  !  M2 ->
                  matchcolp: DO M1=1,NCOL(M2)
                     IF (NVAL(M1,M2).EQ.MINUS(J1)) THEN ! a previous ts also links this pair
                        IF (PLUS(INDEX_TS(M1,M2)) == M2) THEN ! which way round was the original connecting ts?
                           DMAT0(INDEX_TS(M1,M2),1) = MIN(DMAT0(INDEX_TS(M1,M2),1)+EXP(KPLUS(J1)-KSUM(PLUS(J1))),1.0D0)
                        ELSEIF (MINUS(INDEX_TS(M1,M2)) == M2) THEN
                           DMAT0(INDEX_TS(M1,M2),2) = MIN(DMAT0(INDEX_TS(M1,M2),2)+EXP(KPLUS(J1)-KSUM(PLUS(J1))),1.0D0)
                        ELSE
                           PRINT *,'kshortestpaths> error in INDEX_TS connectivities'
                           STOP
                        ENDIF
                        DTS(J1) = .TRUE.
                        MATCHED=.TRUE.
                        EXIT matchcolp
                     ENDIF
                  ENDDO matchcolp
                  IF (.NOT.MATCHED) THEN ! this minimum has not been connected to from M1 before
                     NCOL(M2)=NCOL(M2)+1
                     NVAL(NCOL(M2),M2)=MINUS(J1)
                     INDEX_TS(NCOL(M2),M2) = J1
                     DMAT0(J1,1) = MIN(EXP(KPLUS(J1)-KSUM(PLUS(J1))),1.0D0)
                  ENDIF
               ELSE IF (MINUS(J1).EQ.M2) THEN  !   -> M2 
                  matchcolm: DO M1=1,NCOL(M2)
                     IF (NVAL(M1,M2).EQ.PLUS(J1)) THEN ! a previous ts also links this pair
                        IF (PLUS(INDEX_TS(M1,M2)) == M2) THEN ! which way round was the original connecting ts?
                           DMAT0(INDEX_TS(M1,M2),1) = MIN(DMAT0(INDEX_TS(M1,M2),1)+EXP(KMINUS(J1)-KSUM(MINUS(J1))),1.0D0)
                        ELSEIF (MINUS(INDEX_TS(M1,M2)) == M2) THEN
                           DMAT0(INDEX_TS(M1,M2),2) = MIN(DMAT0(INDEX_TS(M1,M2),2)+EXP(KMINUS(J1)-KSUM(MINUS(J1))),1.0D0)
                        ELSE
                           PRINT *,'kshortestpaths> error in INDEX_TS connectivities'
                           STOP
                        ENDIF
                        DTS(J1) = .TRUE.
                        MATCHED=.TRUE.
                        EXIT matchcolm
                     ENDIF
                  ENDDO matchcolm
                  IF (.NOT.MATCHED) THEN ! this minimum has not been connected to from M1 before
                     NCOL(M2)=NCOL(M2)+1
                     NVAL(NCOL(M2),M2)=PLUS(J1)
                     INDEX_TS(NCOL(M2),M2) = J1
                     DMAT0(J1,2) = MIN(EXP(KMINUS(J1)-KSUM(MINUS(J1))),1.0D0)
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
         IF (PLUS(J1).EQ.M2) THEN
            J1=POINTERP(J1)
         ELSE IF (MINUS(J1).EQ.M2) THEN
            J1=POINTERM(J1)
         ENDIF
         IF (J1.GT.0) GOTO 11
      ENDDO fromloop

      DEADTS = DTS 

      RETURN
      END
