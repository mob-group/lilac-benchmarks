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
MODULE INPUTMODULE
     IMPLICIT NONE
     CONTAINS
     ! SAT 19:46:13 BST 2006
     ! Modified GT's ReadInput routine,
     ! this version reads the data in from PATHSAMPLE's data structures
     SUBROUTINE READINPUT(NDISTA,NDISTB,NCOL,NVAL,PBRANCH,EMKSUM)
          USE COMMONS,ONLY:DEBUG,NORMALISE=>GT2NORMALISE,&
          &NMIN,NMINB,NMINA,NCONNMIN,NCONNMAX,NCONN,PFMIN,LOCATIONB,LOCATIONA,DIRECTION
          USE DATAMODULE,ONLY:NNODES,NSOURCES,NSINKS,NINODES,NUNPROCESSEDNODES,G
          USE DLLMODULE,ONLY:DLLINITN,DLLINITR,DLLADDR,DLLADDN,DLLDIVR
          IMPLICIT NONE
          INTEGER,INTENT(IN),DIMENSION(NMIN)::NCOL,NDISTA,NDISTB
          INTEGER,INTENT(IN),DIMENSION(NCONNMAX,NMIN)::NVAL
          DOUBLE PRECISION,INTENT(IN),DIMENSION(NMIN)::EMKSUM
          DOUBLE PRECISION,INTENT(IN),DIMENSION(NCONNMAX,NMIN)::PBRANCH
          
          INTEGER :: I,J,OLDINDEX,INDEX
          INTEGER,ALLOCATABLE,DIMENSION(:)::I2J,J2I
          DOUBLE PRECISION :: PROB,PSUM,PFTOTALB,PFTOTALA
          PRINT *, "ReadInput> Start"
          ! SAT 12:24:33 BST 2006; determine the number of nodes first
          NNODES=0
!         PRINT *,'readinput> size(nconn)=',size(nconn)
          DO I=1,NMIN
               IF (NCONN(I).LE.NCONNMIN) THEN
                    CYCLE
               ELSE
                    NNODES=NNODES+1
               ENDIF
          ENDDO
          ALLOCATE(G(NNODES),I2J(NMIN),J2I(NNODES))
          DO I=1,NNODES
               ALLOCATE(G(I)%P)
          ENDDO
          PRINT *, "ReadInput> Number of nodes =",Nnodes
          ! i \in [1,Nmin], j \in [1,Nnodes]; fill in index conversion arrays here
          J=0
          DO I=1,NMIN
               IF (NCONN(I).LE.NCONNMIN) THEN
                    CYCLE
               ELSE
                    J=J+1
                    I2J(I)=J
                    J2I(J)=I
               ENDIF
          ENDDO
          ! determine PFTOTALA and PFTOTALB
          PFTOTALB=0.0D0
          DO I=1,NMINB
               IF (NDISTB(LOCATIONB(I)).EQ.0) THEN
                    PFTOTALB=PFTOTALB+EXP(PFMIN(LOCATIONB(I)))
               ENDIF
          ENDDO
          PFTOTALB=LOG(PFTOTALB)
          PFTOTALA=0.0D0
          DO I=1,NMINA
               IF (NDISTA(LOCATIONA(I)).EQ.0) THEN
                    PFTOTALA=PFTOTALA+EXP(PFMIN(LOCATIONA(I)))
               ENDIF
          ENDDO
          PFTOTALA=LOG(PFTOTALA)
          NSOURCES=0; NSINKS=0
          DO I=1,NNODES
               G(I)%P%INDEX = I
               OLDINDEX=J2I(I)
               ! SAT Сбт Авг 19 12:53:20 BST 2006
               ! source/sink allocation determined by DIRECTION from pathdata file DJW
               IF (DIRECTION.EQ.'AB') THEN
                  IF (NDISTA(OLDINDEX).EQ.0) THEN
                       G(I)%P%D = 0
                  ELSE
                       G(I)%P%D = NCOL(OLDINDEX)
                  ENDIF
                  G(I)%P%TAU = EMKSUM(OLDINDEX)
                  IF (NDISTB(OLDINDEX).EQ.0) THEN
                       G(I)%P%PROB = EXP(PFMIN(OLDINDEX)-PFTOTALB)
                  ELSE
                       G(I)%P%PROB = 0.0D0
                  ENDIF
               ELSE
                  IF (NDISTB(OLDINDEX).EQ.0) THEN
                       G(I)%P%D = 0
                  ELSE
                       G(I)%P%D = NCOL(OLDINDEX)
                  ENDIF
                  G(I)%P%TAU = EMKSUM(OLDINDEX)
                  IF (NDISTA(OLDINDEX).EQ.0) THEN
                       G(I)%P%PROB = EXP(PFMIN(OLDINDEX)-PFTOTALA)
                  ELSE
                       G(I)%P%PROB = 0.0D0
                  ENDIF
               ENDIF
               IF (G(I)%P%D==0.AND.G(I)%P%PROB>0.0D0) THEN
                    PRINT *, 'ReadInput> ERROR: node with degree 0 is marked as source'; stop
               ENDIF
               NULLIFY(G(I)%P%P,G(I)%P%A,G(I)%P%LEFT,G(I)%P%RIGHT,G(I)%P%PARENT,G(I)%P%CHILD)
               IF (G(I)%P%D>0) THEN
                    IF (G(I)%P%PROB>0.0D0) THEN
                         IF (DEBUG) PRINT *, "ReadInput> Node",i,"is a source, weight =",G(i)%p%prob
                         G(I)%P%T='s'
                         NSOURCES=NSOURCES+1
                    ELSE
                         IF (DEBUG) PRINT *, "ReadInput> Node",i,"is an intermediate node"
                         G(I)%P%T='i'
                    ENDIF
                    PROB=PBRANCH(1,OLDINDEX)
                    INDEX=I2J(NVAL(1,OLDINDEX))
                    IF (.NOT.(INDEX<=NNODES.AND.INDEX>0)) THEN
                         PRINT *, "ReadInput> ERROR: index out of range 0...Nnodes, index=",index; stop
                    ENDIF
                    CALL DLLINITR(G(I)%P%P,PROB)
                    CALL DLLINITN(G(I)%P%A,G(INDEX)%P)
                    PSUM=PROB
                    IF (G(I)%P%D>1) THEN
                         DO J=2,G(I)%P%D
                              PROB=PBRANCH(J,OLDINDEX)
                              INDEX=I2J(NVAL(J,OLDINDEX))
                              IF (.NOT.(INDEX<=NNODES.AND.INDEX>0)) THEN
                                   PRINT *, "ReadInput> ERROR: index out of range 0...Nnodes, index=",index; stop
                              ENDIF
                              CALL DLLADDR(G(I)%P%P,PROB)
                              PSUM=PSUM+PROB
                              CALL DLLADDN(G(I)%P%A,G(INDEX)%P)
                         ENDDO
                    ENDIF
                    IF (NORMALISE) CALL DLLDIVR(G(I)%P%P,PSUM)
               ELSE
                    IF (DEBUG) PRINT *, "ReadInput> Node",i,"is a sink"
                    G(I)%P%T='t'
                    NSINKS=NSINKS+1
               ENDIF
          ENDDO
          PRINT *, "ReadInput> Number of source nodes =",Nsources
          PRINT *, "ReadInput> Number of sink nodes =",Nsinks
          NINODES=NNODES-NSOURCES-NSINKS
          NUNPROCESSEDNODES=NNODES
          PRINT *, "ReadInput> Number of intermediate nodes =",Ninodes
          CLOSE(UNIT=10)

          DEALLOCATE(I2J,J2I)
          PRINT *, "ReadInput> Done"
     END SUBROUTINE READINPUT
END MODULE INPUTMODULE
