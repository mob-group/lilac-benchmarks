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
! For use with order-parameter-generating protocols, ref.crd and native.xyz both correspond to min 2, i.e., the minimized PDB structure.
! Need an odata file "odata.order" set up as for a minimization, with the CALCDIHE keyword.
! Uses ~/OPTIM/source/CHARMM/calcdihe.src order param routines, hand-crafted to deal with GB1 specificities.
! 
! dbase files have the format:
! minimum       RMSD     dist from min 2      R_g    R_g of heavy atoms of H core  # native cross-chain H-bonds
! 1         9.74705      123.29152       11.49696       11.12788         0
! 2         0.00000        0.00000        7.68081        4.90476         5
! 3         2.83561       35.86794        7.86131        4.81227         2
! 
!
      SUBROUTINE CALCORDER(NATOMS,NMIN,NTS,UMIN,UTS,DEBUG)
      USE COMMONS,ONLY : RIGIDBODY,TWOD,BULKT,BOXLZ,BOXLY,BOXLX, EXEC
      IMPLICIT NONE
      INTEGER J1, J2, NMIN, NATOMS, NTS, UMIN, UTS, STATUS, NNHB
      DOUBLE PRECISION LOCALPOINTS(3*NATOMS), NATIVEPOINTS(3*NATOMS), RMSD, RMAT(3,3), DISTANCE, DIST2, &
     &                 RGYR, RGYRH
      LOGICAL DEBUG, YESNO

      WRITE(*,*) 'calcorder> opening the GB1 database file' 
      OPEN(UNIT=10,FILE='points.min',ACCESS='DIRECT', &
     & FORM='UNFORMATTED',STATUS='OLD',RECL=8*3*NATOMS)
      OPEN(UNIT=91,FILE='unfolded.dbase',STATUS='UNKNOWN')
      OPEN(UNIT=92,FILE='folded.dbase',STATUS='UNKNOWN')
      OPEN(UNIT=93,FILE='all.dbase',STATUS='UNKNOWN')
      OPEN(UNIT=100,FILE='native.xyz',STATUS='OLD')
      READ(100,*) (NATIVEPOINTS(J2),J2=1,3*NATOMS)

      DO J1=1,NMIN
        READ(10,REC=J1) (LOCALPOINTS(J2),J2=1,3*NATOMS)
        CALL MINPERMDIST(NATIVEPOINTS,LOCALPOINTS,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,RIGIDBODY,RMAT,.FALSE.)
        RMSD=0.0D0
        DO J2=1,NATOMS
          RMSD = RMSD + (LOCALPOINTS(3*(J2-1)+1) - NATIVEPOINTS(3*(J2-1)+1))**2 &
     &                + (LOCALPOINTS(3*(J2-1)+2) - NATIVEPOINTS(3*(J2-1)+2))**2 &
     &                + (LOCALPOINTS(3*(J2-1)+3) - NATIVEPOINTS(3*(J2-1)+3))**2
        ENDDO

        RMSD=DSQRT(RMSD/(1.0D0*NATOMS))

        CALL ORDERODATA(1,LOCALPOINTS)
        CALL MYSYSTEM(STATUS,DEBUG,EXEC//' 1 > output.order.1')
        INQUIRE(FILE='order.out',EXIST=YESNO)
        IF (YESNO) THEN
           OPEN(UNIT=1978,FILE='order.out',STATUS='OLD')
           READ(1978,*) RGYR, RGYRH, NNHB
           CLOSE(1978)
        ENDIF

! RMSD from minimized PDB, something out of minpermdist, all-atom radius of gyration, Rgyr for heavy atoms of 
! the hydrophobic core (residues trp, tyr, phe, val), number of native cross-chain H-bonds (in both PDB and mini-PDB).
        WRITE(93,'(I10,1X,4F15.5,I10)') J1, RMSD, DISTANCE, RGYR, RGYRH, NNHB

! These cutoffs may be total bolox.
        IF (NNHB < 4 .AND. RGYRH <= 9.0 .AND. RGYRH >= 5.5) THEN ! unfolded
           WRITE(91,'(I10,1X,4F15.5,I10)') J1, RMSD, DISTANCE, RGYR, RGYRH, NNHB
! These cutoffs have some basis from simulated FES (e.g. Zhou and Berne 2002).
        ELSEIF (NNHB == 5 .AND. RGYRH <= 6.0 .AND. RGYRH >= 5.0) THEN ! folded
           WRITE(92,'(I10,1X,4F15.5,I10)') J1, RMSD, DISTANCE, RGYR, RGYRH, NNHB
        ENDIF
      ENDDO

      CLOSE(10)
      CLOSE(100)
      CLOSE(91)
      CLOSE(92)
      CLOSE(93)
      RETURN 
      END
