!
!  Copyright (C) 2005-2006 Semen Trygubenko
!  This file is part of GT.
!
!  GT is free software; you can redistribute it and/or modify it
!  under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2, or (at your option)
!  any later version.
!
!  GT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with GT; see the file COPYING. If not, write to the
!  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
!  MA 02111-1307 USA
!
MODULE GRAPHTRANSFORMATIONMODULE
     IMPLICIT NONE
     CONTAINS
     SUBROUTINE SPARSEGRAPHTRANSFORMATION
          USE COMMONS, ONLY:DEBUG,SWITCH=>GT2SWITCH,RSWITCH=>GT2RSWITCH
          USE DATAMODULE, ONLY: G,NNODES,NODE
          USE FIBONACCIHEAPMODULE, ONLY: FHINIT,FHINSERT,FHEXTRACTMIN,N
          IMPLICIT NONE
          INTEGER :: I,NSAVE
          DOUBLE PRECISION :: TS,TF
          TYPE(NODE),POINTER :: Z
          LOGICAL :: SGTDONE
          NULLIFY(Z)
          SGTDONE=.FALSE.
          CALL CPU_TIME(TS)
          PRINT *, 'SGT> Start'
          CALL FHINIT
          DO I=1,NNODES
               IF (G(I)%P%T=='i') then
                    CALL FHINSERT(G(I)%P)
                    IF (DEBUG) PRINT *, 'SGT> Put i node',i,'on a heap'
               ENDIF
          ENDDO
          NSAVE=N
          I=N
          DO
               CALL FHEXTRACTMIN(Z)
               IF (ASSOCIATED(Z)) THEN
                    I=I-1
                    CALL SGDETACHNODE(Z)
                    IF (DEBUG) PRINT '(A,I6,A,I6,A,I6)', 'SGT> Iteration',i,': detached i node',z%index,'with d =',z%d
                    IF (N==0) THEN
                         SGTDONE=.TRUE.
                    ELSE IF (SWITCH.AND.REAL(Z%D)/I>RSWITCH) THEN
                         SGTDONE=.TRUE.
                    ENDIF
                    IF (SGTDONE) THEN
                         PRINT *, 'SGT> Detached',nsave-n,'out of',nsave,'intermediate nodes'
                         CALL CPU_TIME(TF)
                         PRINT *, 'SGT> Done. CPU time =',Tf-Ts
                         CALL COMPLETEGRAPHTRANSFORMATION
                         EXIT
                    ENDIF
               ELSE
                    IF (.NOT.N==0) THEN
                         PRINT *, 'SGT> ERROR: not all the nodes extracted from the heap'; stop
                    ELSE !AT THE END OF A CALCULATION SGT MUST CALL CGT TO DISCONNECT SOURCES, EVEN IF SWITCH=F.
                         PRINT *, 'SGT> ERROR: exiting SGT without a call to CGT'; stop
                    ENDIF
               ENDIF
          ENDDO
     END SUBROUTINE SPARSEGRAPHTRANSFORMATION
     SUBROUTINE SGDETACHNODE(C)
     USE COMMONS, ONLY:ALTPBB=>GT2ALTPBB,RESCALE=>GT2RESCALE
     USE DATAMODULE, ONLY: NODE
     USE DLLMODULE, ONLY: NODELIST,REALLIST,DLLDELR,DLLDELN,DLLADDR,DLLADDN,DLLDIVR,DLLSUMR,DLLMULR
     USE FIBONACCIHEAPMODULE, ONLY: FHDECREASEKEY,FHINCREASEKEY
     IMPLICIT NONE
     INTEGER :: DOLD,DNEW,FIRSTINDEX1,FIRSTINDEX2,FIRSTINDEX3,FIRSTINDEX4,FIRSTINDEX5,FIRSTINDEX6
     DOUBLE PRECISION :: PCB_C,PBC_B,PCA_C,PBB,ONEMPCB_C,ONEMPBC_B,PNORM
     TYPE(NODE),POINTER :: A,B,C
     TYPE(NODELIST),POINTER :: BA1, BA2, CA1, CA2
     TYPE(REALLIST),POINTER :: BP1, BP2, CP1, CP2
     LOGICAL :: FIRSTPASS1,FIRSTPASS2,FIRSTPASS3,FIRSTPASS4,FIRSTPASS5,FIRSTPASS6
     IF (ASSOCIATED(C%A)) THEN
          CA1=>C%A
          FIRSTINDEX1=CA1%VALUE%INDEX
          FIRSTPASS1=.TRUE.
          CP1=>C%P
          DO ! THIS DO LOOP IS OVER ALL B_C \IN ADJ[C]
               IF (.NOT.FIRSTINDEX1==CA1%VALUE%INDEX.OR.FIRSTPASS1) THEN
                    FIRSTPASS1=.FALSE.
                    B=>CA1%VALUE
                    PCB_C=CP1%VALUE
                    IF (ASSOCIATED(B%A)) THEN ! IT IS NOT A SINK
                         BA1=>B%A
                         FIRSTINDEX4=BA1%VALUE%INDEX
                         FIRSTPASS4=.TRUE.
                         BP1=>B%P
                         DO ! THIS DO LOOP IS OVER ALL C_B \IN ADJ[B]
                              IF (.NOT.FIRSTINDEX4==BA1%VALUE%INDEX.OR.FIRSTPASS4) THEN
                                   FIRSTPASS4=.FALSE.
                                   IF (BA1%VALUE%INDEX==C%INDEX) THEN ! WE HAVE FOUND C IN ADJACENCY LIST OF B
                                        PBC_B=BP1%VALUE
                                        IF (ALTPBB) THEN
                                             CALL CALCULATEPBB
                                        ELSE
                                             PBB=1.0D0-PBC_B*PCB_C
                                        ENDIF
                                        IF (PBB==0.0D0) THEN ! PRINT INFO HERE RATHER THAN BELOW AFTER ADJACENCY
                                             CALL PRINTNODEINFO(B) ! LISTS OF NODE C ARE MODIFIED, WHICH PRODUCES
                                             CALL PRINTNODEINFO(C) ! SOMEWHAT LESS STARTLING OUTPUT :-)
                                        ENDIF
                                        B%A=>BA1
                                        B%P=>BP1
                                        DOLD=B%D
                                        CALL DLLDELN(B%A)
                                        CALL DLLDELR(B%P)
                                        DNEW=DOLD-1
                                        CA2=>CA1
                                        FIRSTINDEX2=CA2%VALUE%INDEX
                                        FIRSTPASS2=.TRUE.
                                        CP2=>CP1
                                        DO
                                             IF (.NOT.CA2%VALUE%INDEX==FIRSTINDEX2.OR.FIRSTPASS2) THEN
                                                  FIRSTPASS2=.FALSE.
                                                  A=>CA2%VALUE
                                                  IF (.NOT.A%INDEX==B%INDEX) THEN
                                                       PCA_C=CP2%VALUE
                                                       IF (ASSOCIATED(B%A)) THEN
                                                            BA1=>B%A
                                                            FIRSTPASS5=.TRUE.
                                                            FIRSTINDEX5=BA1%VALUE%INDEX
                                                            BP1=>B%P
                                                            DO
                                                                 IF (.NOT.FIRSTINDEX5==BA1%VALUE%INDEX.OR.FIRSTPASS5) THEN
                                                                      FIRSTPASS5=.FALSE.
                                                                      IF (BA1%VALUE%INDEX==A%INDEX) THEN
                                                                           BP1%VALUE=BP1%VALUE+PBC_B*PCA_C
                                                                           NULLIFY(BA1,BP1)
                                                                           EXIT
                                                                      ELSE
                                                                           BA1=>BA1%F
                                                                           BP1=>BP1%F
                                                                      ENDIF
                                                                 ELSE
                                                                      NULLIFY(BA1,BP1)
                                                                      CALL DLLADDN(B%A,A)
                                                                      CALL DLLADDR(B%P,PBC_B*PCA_C)
                                                                      DNEW=DNEW+1
                                                                      EXIT
                                                                 ENDIF
                                                            ENDDO
                                                       ELSE
                                                            CALL DLLADDN(B%A,A)
                                                            CALL DLLADDR(B%P,PBC_B*PCA_C)
                                                            DNEW=DNEW+1
                                                       ENDIF
                                                  ENDIF
                                                  CA2=>CA2%F
                                                  CP2=>CP2%F
                                             ELSE
                                                  NULLIFY(CA2,CP2)
                                                  EXIT
                                             ENDIF
                                        ENDDO
                                        IF (PBB/=0.0D0) THEN
                                             IF (PBB<0.0D0) THEN
                                                  PRINT *,'SGDetachNode> WARNING: 0>Pbb^-1=',Pbb
                                             ENDIF
                                             CALL DLLDIVR(B%P,PBB)
                                             IF (RESCALE) THEN
                                                  CALL DLLSUMR(B%P,PNORM)
                                                  CALL DLLDIVR(B%P,PNORM)
                                             ENDIF
                                             B%TAU=(B%TAU+C%TAU*PBC_B)/PBB
                                             IF (B%T/='s') then
                                                  IF (DNEW>DOLD) THEN
                                                       CALL FHINCREASEKEY(B,DNEW)
                                                  ELSEIF (DNEW<DOLD) THEN
                                                       CALL FHDECREASEKEY(B,DNEW)
                                                  ENDIF
                                             ENDIF
                                        ELSE
                                             IF (DEADEND(B).AND.DEADEND(C)) THEN
                                                  PRINT *,'SGDetachNode> Disconnected subgraph detected (effectively)'
                                                  PRINT *,'SGDetachNode> Scheduling node',b%index,' for deletion'
                                                  IF (B%T/='s') CALL FHDECREASEKEY(B,0) ! This ensures b will be detached next !
                                                  CALL DLLMULR(B%P,0.0D0) ! BP's do not add up to 1 anymore !
                                                  B%TAU=0.0D0 ! THIS SHOULD BE INFINITY INSTEAD !
                                             ELSE
                                                 PRINT *,'SGDetachNode> ERROR: Pbb is infinite'; stop
                                             ENDIF
                                        ENDIF
                                        EXIT
                                   ENDIF
                              ENDIF
                              BA1=>BA1%F
                              BP1=>BP1%F
                         ENDDO
                    ENDIF
                    CA1=>CA1%F
                    CP1=>CP1%F
               ELSE
                    NULLIFY(CA1,CP1)
                    EXIT
               ENDIF
          ENDDO
     ELSE
          PRINT *,'SGDetachNode> ERROR: node to be detached has zero degree';stop
     ENDIF
     C%T='d' ! the node is now detached
     CONTAINS
     SUBROUTINE CALCULATEPBB
          IMPLICIT NONE
          ! we evaluate Pbb differently to maintain precision
          ! brings down the abs(1-Ptot) by up to an order of magnitude at extra computational cost
          !print *, 'Pcb_c=',Pcb_c
          ONEMPCB_C=0.0D0
          CA2=>CA1
          FIRSTPASS3=.TRUE.
          FIRSTINDEX3=CA2%VALUE%INDEX
          CP2=>CP1
          DO
               IF (.NOT.FIRSTINDEX3==CA2%VALUE%INDEX.OR.FIRSTPASS3) THEN
                    FIRSTPASS3=.FALSE.
                    IF (CA2%VALUE%INDEX/=CA1%VALUE%INDEX) THEN
                         ONEMPCB_C=ONEMPCB_C+CP2%VALUE
                    ENDIF 
                    CA2=>CA2%F
                    CP2=>CP2%F
               ELSE
                    NULLIFY(CA2,CP2)
                    EXIT
               ENDIF
          ENDDO
          !print *, 'Pcb_c=',1.0d0-onemPcb_c
          !print *, 'Pbc_b=',Pbc_b
          ONEMPBC_B=0.0D0
          BA2=>BA1
          FIRSTPASS6=.TRUE.
          FIRSTINDEX6=BA2%VALUE%INDEX
          BP2=>BP1
          DO
               IF (.NOT.FIRSTINDEX6==BA2%VALUE%INDEX.OR.FIRSTPASS6) THEN
                    FIRSTPASS6=.FALSE.
                    IF (BA2%VALUE%INDEX/=BA1%VALUE%INDEX) THEN
                         ONEMPBC_B=ONEMPBC_B+BP2%VALUE
                    ENDIF 
                    BA2=>BA2%F
                    BP2=>BP2%F
               ELSE
                    NULLIFY(BA2,BP2)
                    EXIT
               ENDIF
          ENDDO
          !print *, 'Pbc_b=',1.0d0-onemPbc_b
          PBB=ONEMPBC_B+ONEMPCB_C-ONEMPBC_B*ONEMPCB_C
     END SUBROUTINE CALCULATEPBB
     END SUBROUTINE SGDETACHNODE
     SUBROUTINE COMPLETEGRAPHTRANSFORMATION
          USE DATAMODULE,ONLY: NUNPROCESSEDNODES,G,P,TAU
          USE COMMONS, ONLY:DEBUG
          IMPLICIT NONE
          INTEGER :: I,J,IFIRST,SFIRST
          DOUBLE PRECISION :: TS,TF
          CALL CHANGEDS(IFIRST,SFIRST)
          IF (DEBUG.AND.NUNPROCESSEDNODES<9) CALL PRINTMATRIX(SIZE(P,1),SIZE(P,2),P)
          CALL CPU_TIME(TS)
          PRINT *, 'CGT> Start'
          IF (DEBUG) PRINT *, 'ifirst,sfirst=',ifirst,sfirst
          DO I=IFIRST,SFIRST-1
               CALL CGDETACHNODE(I)
               IF (DEBUG) PRINT *, 'CGT> Detached intermediate node',G(i)%p%index,'( old index is',G(i)%p%indexorig,')'
          ENDDO
          PRINT *, 'CGT> Detached',sfirst-ifirst,'intermediate nodes'
          IF (DEBUG.AND.NUNPROCESSEDNODES<9) CALL PRINTMATRIX(SIZE(P,1),SIZE(P,2),P)
          DO I=SFIRST,NUNPROCESSEDNODES
               DO J=SFIRST,NUNPROCESSEDNODES
                    IF (I/=J) THEN
                         IF (P(I,J)>0.0D0) THEN
                              CALL DISCONNECT(I,J,IFIRST,SFIRST)
                              IF (DEBUG) PRINT *, 'CGT> Disconnected nodes',i,'and',j
                         ENDIF
                    ENDIF
               ENDDO
          ENDDO
          PRINT *, 'CGT> Disconnected',NunprocessedNodes-sfirst+1,'source nodes'
          CALL CPU_TIME(TF)
          PRINT *, 'CGT> Done. CPU time =',Tf-Ts
          IF (DEBUG.AND.NUNPROCESSEDNODES<9) CALL PRINTMATRIX(SIZE(P,1),SIZE(P,2),P)
          CALL OUTPUT(IFIRST,SFIRST)
          DEALLOCATE(P,TAU)
     END SUBROUTINE COMPLETEGRAPHTRANSFORMATION
     SUBROUTINE CHANGEDS(IFIRST,SFIRST,SOURCES)
          USE DATAMODULE,ONLY: NNODES,NSOURCES,P,TAU,G,NODELIST,REALLIST,NODEARRAY,NUNPROCESSEDNODES
          USE COMMONS, ONLY:DEBUG
          IMPLICIT NONE
          INTEGER,INTENT(OUT) :: IFIRST
          INTEGER,INTENT(OUT),OPTIONAL :: SFIRST
          INTEGER,POINTER,DIMENSION(:),OPTIONAL :: SOURCES
          INTEGER :: I,J,K,LASTINDEX,NCON,SOURCEINDEX
          DOUBLE PRECISION :: TS,TF
          TYPE(REALLIST),POINTER :: BP
          TYPE(NODELIST),POINTER :: BA
          TYPE(NODEARRAY),ALLOCATABLE :: GTMP(:)
          IF (PRESENT(SFIRST).AND.PRESENT(SOURCES)) THEN
               PRINT*,'ChangeDS> ERROR: Arguments sfirst and sources are mutually exclusive';stop
          ENDIF
          ALLOCATE(GTMP(NNODES))
          DO I=1,NNODES
               NULLIFY(GTMP(I)%P)
          ENDDO
          CALL CPU_TIME(TS)
          PRINT *, 'ChangeDS> Start'
          ! sink nodes go first
          K=0
          DO I=1,NNODES
               IF (G(I)%P%T=='t') then
                    K=K+1
                    GTMP(K)%P=>G(I)%P
               ENDIF
          ENDDO
          IFIRST=K+1 ! INDEX OF THE FIRST INTERMEDIATE NODE
          ! next block of code orders intermediate and source nodes; order is different if this is a chain graph
          IF (PRESENT(SOURCES)) THEN
          ! block of Pij for i,j > ifirst: nodes are ordered so that only the elements of upper and lower diagonal are non-zero
               ALLOCATE(SOURCES(NSOURCES)); SOURCES=0; SOURCEINDEX=0
               DO I=1,NNODES ! FIND A NODE THAT IS CONNECTED TO ONLY ONE INTERMEDIATE NODE; IGNORE CONNECTIONS TO SINKS
                    IF (G(I)%P%T/='t') then
                         LASTINDEX=G(I)%P%A%VALUE%INDEX
                         BA=>G(I)%P%A%F
                         NCON=0 
                         DO
                              IF (BA%VALUE%T/='t') ncon=ncon+1
                              IF (BA%VALUE%INDEX==LASTINDEX) EXIT
                              BA=>BA%F
                         ENDDO
                         NULLIFY(BA)
                         IF (NCON==1) THEN ! VERIFY THAT THIS IS INDEED A CHAIN GRAPH AND INITIALISE TERMINAL NODES
                              IF (.NOT.ASSOCIATED(GTMP(IFIRST)%P)) THEN
                                   K=K+1
                                   GTMP(IFIRST)%P=>G(I)%P
                                   IF (GTMP(IFIRST)%P%T=='s') then
                                        SOURCEINDEX=SOURCEINDEX+1
                                        SOURCES(SOURCEINDEX)=IFIRST
                                   ENDIF
                                   GTMP(IFIRST)%P%T='p'
                              ELSEIF (.NOT.ASSOCIATED(GTMP(NNODES)%P)) THEN
                                   GTMP(NNODES)%P=>G(I)%P
                                   IF (GTMP(NNODES)%P%T=='s') then
                                        SOURCEINDEX=SOURCEINDEX+1
                                        SOURCES(SOURCEINDEX)=NNODES
                                   ENDIF
                                   GTMP(NNODES)%P%T='p'
                              ELSE
                                   PRINT*,'ChangeDS> ERROR: too many terminal nodes'; stop
                              ENDIF
                         ELSEIF (NCON>2) THEN
                              PRINT*,'ChangeDS> ERROR: a node with more than two connections to intermediate nodes encountered';stop
                         ELSEIF (NCON<1) THEN
                              PRINT*,'ChangeDS> ERROR: unconnected graph'; stop
                         ENDIF
                    ENDIF
               ENDDO
               MAIN: DO ! ORDER NON-TERMINAL NODES APPROPRIATELY
                    LASTINDEX=GTMP(K)%P%A%VALUE%INDEX
                    BA=>GTMP(K)%P%A%F
                    DO
                         IF (BA%VALUE%T/='t'.and.bA%value%t/='p') then
                              K=K+1
                              GTMP(K)%P=>G(BA%VALUE%INDEX)%P
                              IF (GTMP(K)%P%T=='s') then
                                   SOURCEINDEX=SOURCEINDEX+1
                                   SOURCES(SOURCEINDEX)=K
                              ENDIF
                              GTMP(K)%P%T='p'
                              NULLIFY(BA)
                              EXIT
                         ENDIF
                         IF (BA%VALUE%INDEX==LASTINDEX) THEN
                              NULLIFY(BA)
                              IF (K/=NNODES-1) THEN
                                   PRINT*,'ChangeDS> ERROR: unspecified';stop
                              ELSE
                                   EXIT MAIN
                              ENDIF
                         ENDIF
                         BA=>BA%F
                    ENDDO
               ENDDO MAIN
               K=K+1
          ELSE ! SOURCE NODES GO LAST
               DO I=1,NNODES
                    IF (G(I)%P%T=='i') then
                         K=K+1
                         GTMP(K)%P=>G(I)%P
                    ENDIF
               ENDDO
               SFIRST=K+1 ! INDEX OF THE FIRST SOURCE NODE
               DO I=1,NNODES
                    IF (G(I)%P%T=='s') then
                         K=K+1
                         GTMP(K)%P=>G(I)%P
                    ENDIF
               ENDDO
          ENDIF
          NUNPROCESSEDNODES=K!=NSOURCES+NSINKS+INTERMEDIATE NODES THAT WERE NOT DETACHED YET
          IF (DEBUG) PRINT*,'ChangeDS> NunprocessedNodes=',NunprocessedNodes
          DO I=1,NUNPROCESSEDNODES
               GTMP(I)%P%INDEXORIG=GTMP(I)%P%INDEX
               GTMP(I)%P%INDEX=I
          ENDDO
          K=NUNPROCESSEDNODES
          DO I=1,NNODES ! THIS DO LOOP IS NECESSARY ONLY FOR PROPER MEMORY MANAGEMENT
               IF (G(I)%P%T=='d') then
                    K=K+1
                    GTMP(K)%P=>G(I)%P
                    NULLIFY(G(I)%P)
               ENDIF
          ENDDO
          IF (K/=NNODES) THEN
               PRINT*,'ChangeDS> ERROR: k,Nnodes=',k,Nnodes; stop
          ENDIF
          DO I=1,NNODES
               G(I)%P=>GTMP(I)%P
               NULLIFY(GTMP(I)%P)
          ENDDO
          DEALLOCATE(GTMP)
          ALLOCATE(P(NUNPROCESSEDNODES,IFIRST:NUNPROCESSEDNODES),TAU(NUNPROCESSEDNODES))
          P=0.0D0 ! INITIALISE TRANSITION PROBABILITY MATRIX P
          DO J=IFIRST,NUNPROCESSEDNODES
               IF (ASSOCIATED(G(J)%P%P)) THEN
                    BP=>G(J)%P%P%F
                    BA=>G(J)%P%A%F
                    NULLIFY(BP%B)
                    I=-1
                    DO
                         IF (ASSOCIATED(BP%B).OR.I==-1) THEN
                              I=I+1
                              P(BA%VALUE%INDEX,J)=BP%VALUE
                              BP=>BP%F
                              BA=>BA%F
                         ELSE
                              G(J)%P%P%F%B=>G(J)%P%P
                              NULLIFY(BP,BA)
                              EXIT
                         ENDIF
                    ENDDO
               ENDIF
               TAU(J)=G(J)%P%TAU
          ENDDO
          CALL CPU_TIME(TF)
          PRINT *, 'ChangeDS> Done. CPU time =',Tf-Ts
     END SUBROUTINE CHANGEDS
     SUBROUTINE OUTPUT(IFIRST,SFIRST)
          USE DATAMODULE,ONLY:G,P,TAU,S1,S2,NUNPROCESSEDNODES
          USE COMMONS, ONLY:DEBUG,DIRECTION
          IMPLICIT NONE
          INTEGER,INTENT(IN) :: IFIRST,SFIRST
          INTEGER :: I,J
          DOUBLE PRECISION :: TAUTOT, KTOT
          TAUTOT=0.D0
          KTOT=0.D0
          DO I=SFIRST,NUNPROCESSEDNODES
               WRITE(S1,'(i20)') i
               WRITE(S2,'(i20)') G(i)%p%indexorig
               PRINT *, 'Output> Source '//trim(adjustl(s1))//' (old index '//trim(adjustl(s2))//'):'
               PRINT *, 'Output>    Tau =',tau(i)
               TAUTOT=TAUTOT+TAU(I)*G(I)%P%PROB
               KTOT=KTOT+G(I)%P%PROB/TAU(I)
               IF (DEBUG) THEN
                    DO J=1,IFIRST-1
                         PRINT *, 'Output>    P(',j,'<---',i,')=',P(j,i)
                    ENDDO
               ENDIF
               PRINT *, 'Output>    P sum =',sum(P(1:sfirst-1,i))
          ENDDO
          PRINT *,'Output> Tau total =',TauTot
          IF (DIRECTION.EQ.'AB') PRINT '(A,G20.10)','Output> rate constant k(A<-B) =',KTot
          IF (DIRECTION.EQ.'BA') PRINT '(A,G20.10)','Output> rate constant k(B<-A) =',KTot
     END SUBROUTINE OUTPUT
     SUBROUTINE PRINTMATRIX(NROW,NCOL,P,NROWSTART)
          IMPLICIT NONE
          INTEGER,INTENT(IN) :: NROW,NCOL
          INTEGER,INTENT(IN),OPTIONAL :: NROWSTART
          DOUBLE PRECISION,DIMENSION(NROW,NCOL),INTENT(IN) :: P
          INTEGER :: I,ISTART
          CHARACTER(LEN=100) :: F
          WRITE(F,'(i20)') ncol
          F='(1x,1a,'//trim(adjustl(f))//'g25.15,1x,1a)'
          IF (PRESENT(NROWSTART)) THEN
               ISTART=NROWSTART
          ELSE
               ISTART=1
          ENDIF
          DO I=ISTART,NROW
               WRITE(*,TRIM(ADJUSTL(F))) '|',P(i,:),'|'
          ENDDO
     END SUBROUTINE PRINTMATRIX
     SUBROUTINE CGDETACHNODE(C)
          USE COMMONS, ONLY:ALTPBB=>GT2ALTPBB,RESCALE=>GT2RESCALE
          USE DATAMODULE,ONLY:NUNPROCESSEDNODES,P,TAU
          IMPLICIT NONE
          INTEGER,INTENT(IN) :: C
          INTEGER :: A,B
          DOUBLE PRECISION :: D,ONEMPBC,ONEMPCB,PNORM
          DO B=C+1,NUNPROCESSEDNODES
               IF (P(C,B)>0.0D0) THEN
                    IF (ALTPBB) THEN
                         ONEMPBC=0.0D0
                         DO A=1,B-1
                              ONEMPBC=ONEMPBC+P(A,C)
                         ENDDO
                         DO A=B+1,NUNPROCESSEDNODES
                              ONEMPBC=ONEMPBC+P(A,C)
                         ENDDO
                         !print *, P(b,c),1.0d0-onemPbc
                         ONEMPCB=0.0D0
                         DO A=1,C-1
                              ONEMPCB=ONEMPCB+P(A,B)
                         ENDDO
                         DO A=C+1,NUNPROCESSEDNODES
                              ONEMPCB=ONEMPCB+P(A,B)
                         ENDDO
                         !print *, P(c,b),1.0d0-onemPcb
                         D=ONEMPBC+ONEMPCB-ONEMPBC*ONEMPCB
                    ELSE
                         D=1.0D0-P(B,C)*P(C,B)
                    ENDIF
                    TAU(B)=(TAU(B)+TAU(C)*P(C,B))/D
                    DO A=1,NUNPROCESSEDNODES
                         IF (A/=B.AND.A/=C) THEN
                              P(A,B)=(P(A,B)+P(A,C)*P(C,B))/D
                         ENDIF
                    ENDDO
                    P(C,B)=0.0D0
               ENDIF
               IF (RESCALE) THEN
                    PNORM=0.0D0
                    DO A=1,NUNPROCESSEDNODES
                         PNORM=PNORM+P(A,B)
                    ENDDO
                    !print *, 'Pnorm for node',b,'=',Pnorm
                    DO A=1,NUNPROCESSEDNODES
                         P(A,B)=P(A,B)/PNORM
                    ENDDO
               ENDIF
          ENDDO
          P(:,C)=0.0D0 ! NOT ABSOLUTELY NECESSARY, BUT THE FINAL MATRIX LOOKS NEATER
     END SUBROUTINE CGDETACHNODE
     SUBROUTINE DISCONNECT(A,B,IFIRST,SFIRST)
          USE DATAMODULE,ONLY:NUNPROCESSEDNODES,P,TAU
          IMPLICIT NONE
          INTEGER,INTENT(IN) :: A,B,IFIRST,SFIRST
          INTEGER :: I
          DOUBLE PRECISION :: PAB,PBA,D,TMP
          PAB=P(A,B)
          PBA=P(B,A)
          D=1.0D0-PAB*PBA
          TMP=TAU(A)
          TAU(A)=(TAU(A)+TAU(B)*PBA)/D
          TAU(B)=(TAU(B)+TMP*PAB)/D
          DO I=1,IFIRST-1
               TMP=P(I,A)
               P(I,A)=(P(I,A)+P(I,B)*PBA)/D
               P(I,B)=(P(I,B)+TMP*PAB)/D
          ENDDO                                
          DO I=SFIRST,NUNPROCESSEDNODES
               TMP=P(I,A)                     
               P(I,A)=(P(I,A)+P(I,B)*PBA)/D
               P(I,B)=(P(I,B)+TMP*PAB)/D
          ENDDO
          P(A,B)=0.0D0
          P(B,A)=0.0D0
     END SUBROUTINE DISCONNECT
     SUBROUTINE PRINTNODEINFO(A)
          USE DATAMODULE,ONLY:NODE
          USE DLLMODULE,ONLY:DLLPRINTR,DLLPRINTN
          IMPLICIT NONE
          TYPE(NODE),POINTER :: A
          PRINT *
          PRINT *,'PrintNodeInfo> Information for node with index =',a%index
          PRINT *,'PrintNodeInfo> b%d=',a%d
          PRINT *,'PrintNodeInfo> b%t=',a%t
          CALL DLLPRINTR(A%P)
          CALL DLLPRINTN(A%A)
          PRINT *
     END SUBROUTINE PRINTNODEINFO
     FUNCTION DEADEND(A) ! NODE A IS A DEAD END IF THERE's only one non-zero probability exit from it
          USE DATAMODULE,ONLY:NODE,REALLIST,NODELIST
          USE COMMONS, ONLY:PTOL=>GT2PTOL
          IMPLICIT NONE
          LOGICAL :: DEADEND
          TYPE(NODE),POINTER :: A
          TYPE(REALLIST),POINTER :: SP
          TYPE(NODELIST),POINTER :: SA
          INTEGER :: I,FIRSTINDEX
          LOGICAL :: FIRSTPASS
          NULLIFY(SP,SA)
          DEADEND=.TRUE.
          IF (ASSOCIATED(A%A)) THEN
               PRINT *,'DeadEnd> Checking node',a%index
               SP=>A%P
               SA=>A%A
               I=SA%VALUE%INDEX
               FIRSTPASS=.TRUE.
               FIRSTINDEX=0
               DO
                    IF (FIRSTPASS) THEN
                         FIRSTPASS=.FALSE.
                    ELSE IF (SA%VALUE%INDEX==I) THEN
                         EXIT
                    ENDIF
                    FIRSTINDEX=FIRSTINDEX+1
                    WRITE(*,'(1x,a,i7,a,i7,g20.10)') 'DeadEnd> Exit',firstindex,&
                    &' index,prob=',sa%value%index,sp%value
                    IF ( SP%VALUE<PTOL .OR. ABS(SP%VALUE-1.0D0)<PTOL ) THEN
                         CONTINUE
                    ELSE
                         DEADEND=.FALSE.
                         EXIT
                    ENDIF
                    SP=>SP%F
                    SA=>SA%F
               ENDDO
               NULLIFY(SP,SA)
          ELSE
               PRINT *,'DeadEnd> a%index=',a%index
               PRINT *,'DeadEnd> ERROR: Adjacency list pointer is disassociated'; stop
          ENDIF
     END FUNCTION DEADEND
END MODULE GRAPHTRANSFORMATIONMODULE
