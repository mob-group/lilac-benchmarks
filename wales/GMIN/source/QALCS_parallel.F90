SUBROUTINE QALCS_PARALLEL(ITER,TIME,BRUN,QDONE,SCREENC)
  !
  USE COMMONS, ONLY : MYUNIT, NATOMS
  !
  IMPLICIT NONE
  !
  ! Passed variables for QUENCH
  INTEGER, INTENT(INOUT) :: ITER, BRUN, QDONE
  DOUBLE PRECISION, INTENT(INOUT) :: TIME, SCREENC(3*NATOMS)
  !
  LOGICAL :: SWAP
  INTEGER :: I
  DOUBLE PRECISION :: POTEL
  COMMON /MYPOT/ POTEL  
  !
  I=0
  SWAP=.TRUE.
  DO WHILE(SWAP)
     CALL SCAN_ALL_SWAPS_PARALLEL(ITER,TIME,BRUN,QDONE,SCREENC,SWAP)
     IF(SWAP) WRITE(MYUNIT,'(A, G20.10)') &
          'QALCS_parallel> New E= ', POTEL
     I=I+1
  ENDDO
  WRITE(MYUNIT,'(A,I6,A)') 'QALCS_parallel> Converged in ', I, &
       ' swaps.'
  !
  RETURN
  !
END SUBROUTINE QALCS_PARALLEL

!SUBROUTINE EXE_BEST_SWAP(ITER,TIME,BRUN,QDONE,SCREENC,SWAP)
  !
  !
!END SUBROUTINE EXE_BEST_SWAP

SUBROUTINE SCAN_ALL_SWAPS_PARALLEL(ITER,TIME,BRUN,QDONE,SCREENC,SWAP)
  !
  USE COMMONS, ONLY : NATOMS,COORDS,COORDSO,LABELS,LABELSO, &
       MYUNIT, MYNODE, QALCS_NBRHD, NPAR_GBH, ECONV, NQ, RMS, &
       QALCSV, TARGET, HIT
  !
  IMPLICIT NONE
  !
#ifdef MPI
  INCLUDE 'mpif.h'
  INTEGER MPIERR
#endif
  !
  ! Passed variables for QUENCH
  INTEGER, INTENT(INOUT) :: ITER, BRUN, QDONE
  DOUBLE PRECISION, INTENT(INOUT) :: TIME, SCREENC(3*NATOMS)
  ! Indicator of whether an improvement has been made
  LOGICAL, INTENT(OUT) :: SWAP
  !
  INTEGER :: I, J, K, NBRHOOD(2*QALCS_NBRHD), NPPN, MYI1, MYI2, &
       LMIN(NATOMS), IPROC
  DOUBLE PRECISION :: POTEL, POTELO, EMIN, XMIN(3*NATOMS), &
       POTEL_LIST(NPAR_GBH), ELOWEST
  !
  ! Energy of COORDS from last quench. Common block in QUENCH.
  COMMON /MYPOT/ POTEL  
  !
  IF(QALCSV) WRITE(MYUNIT,'(A)') &
       'scan_all_swaps_parallel> Building neighbourhood list...'
  ! First build the list of atom pairs to be swapped
  K=1
  DO I=1,NATOMS-1
     DO J=I+1, NATOMS
        IF(LABELS(I,1) /= LABELS(J,1)) THEN
           IF(K+1<=2*QALCS_NBRHD) THEN
              NBRHOOD(K) = I; NBRHOOD(K+1) = J
              K=K+2
           ELSE
              WRITE(MYUNIT, '(A)') &
                   'qalcs_parallel> Neighbourhood size overflow!'
              
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  !
  IF(K-1 == 2*QALCS_NBRHD) THEN
     IF(QALCSV) WRITE(MYUNIT,'(A, I6)') &
          'scan_all_swaps_parallel> Built neighbourhood of size ',&
          QALCS_NBRHD
  ELSE
     WRITE(MYUNIT,'(A, I6, A, I6)') &
          'scan_all_swaps_parallel> ERROR: neighbourhood size ',&
          K,' not equal to the expected ', 2*QALCS_NBRHD
  ENDIF
  !
  COORDSO(1:3*NATOMS,1) = COORDS(1:3*NATOMS,1)
  LABELSO(1:NATOMS,1) = LABELS(1:NATOMS,1)
  POTELO = POTEL
  XMIN(1:3*NATOMS) = COORDS(1:3*NATOMS,1)
  LMIN(1:NATOMS) = LABELS(1:NATOMS,1)
  EMIN = POTEL
  !
  IF(QALCSV) WRITE(MYUNIT,'(A)') &
       'scan_all_swaps_parallel> Parallel scan of neighbourhood...'
  !
  IF(MOD(QALCS_NBRHD,NPAR_GBH)==0) THEN
     NPPN = QALCS_NBRHD/NPAR_GBH
  ELSE
     NPPN = QALCS_NBRHD/NPAR_GBH + 1 ! integer division!!!
  ENDIF
  MYI1 = 2*MYNODE*NPPN+1
  MYI2 = MIN(MYI1-1+2*NPPN, 2*QALCS_NBRHD)
  DO I = MYI1, MYI2, 2
     !
     !WRITE(MYUNIT,'(A,I6,I6)') &
     !     'scan_all_swaps_parallel> Swapping atoms', &
     !     NBRHOOD(I),NBRHOOD(I+1)
     ! Perform swap and quench
     J=LABELS(NBRHOOD(I),1)
     LABELS(NBRHOOD(I),1) = LABELS(NBRHOOD(I+1),1)
     LABELS(NBRHOOD(I+1),1) = J
     !
     NQ(1) = NQ(1) + 1
     CALL QUENCH(.FALSE.,1,ITER,TIME,BRUN,QDONE,SCREENC)
     IF(QALCSV) WRITE(MYUNIT,'(A,I10,A,F20.10,A,I5,A,G12.5)') &
          'Qu ',NQ(1),' E= ',POTEL,' steps= ',ITER,' RMS= ',RMS
     !
     ! Check for lower energy
     IF(POTEL < EMIN - ECONV) THEN
        EMIN = POTEL
        XMIN(:) = COORDS(:,1)
        LMIN(:) = LABELS(:,1)
     ENDIF
     !
     ! Undo swap and revert state
     !CALL SWAP_LABELS(NBRHOOD(I),NBRHOOD(I+1), 1)
     POTEL = POTELO
     COORDS(:,1) = COORDSO(:,1)
     LABELS(:,1) = LABELSO(:,1)
     !
  ENDDO
  !
#ifdef MPI
  ! Gather all parallel quench energies on master node.             
  CALL MPI_GATHER(EMIN, 1, MPI_DOUBLE_PRECISION, &
       POTEL_LIST, 1, MPI_DOUBLE_PRECISION, 0, &
       MPI_COMM_WORLD, MPIERR)
  !
  IF(MYNODE==0) THEN
     ELOWEST=MINVAL(POTEL_LIST)
     IPROC=MINLOC(POTEL_LIST, 1)-1
  ENDIF
  !
  ! Broadcast IPROC from master
  CALL MPI_BCAST(IPROC,1,MPI_INTEGER,&
       0,MPI_COMM_WORLD,MPIERR)
  ! Broadcast state from IPROC
  CALL MPI_BCAST(XMIN,3*NATOMS,&
       MPI_DOUBLE_PRECISION,IPROC,MPI_COMM_WORLD,MPIERR)
  CALL MPI_BCAST(LMIN,NATOMS,MPI_INTEGER,&
       IPROC,MPI_COMM_WORLD,MPIERR)
  CALL MPI_BCAST(EMIN,1,MPI_DOUBLE_PRECISION,&
       IPROC,MPI_COMM_WORLD,MPIERR)
  IF(TARGET) THEN 
     CALL MPI_BCAST(HIT,1,MPI_LOGICAL,&
          IPROC,MPI_COMM_WORLD,MPIERR)
     IF(HIT) WRITE(MYUNIT,'(A,I5)') &
          'QALCS_parallel> Target hit deterministically by node', &
          IPROC+1
  ENDIF
#endif
  !
  IF(EMIN < POTELO - ECONV) THEN
     COORDS(1:3*NATOMS,1) = XMIN(1:3*NATOMS)
     LABELS(1:NATOMS,1) = LMIN(1:NATOMS)
     POTEL=EMIN
     SWAP=.TRUE.
     IF(QALCSV) WRITE(MYUNIT,'(A, I5)') &
          'scan_all_swaps_parallel> Best improvement on node ',IPROC+1
  ELSE
     SWAP=.FALSE.
     IF(QALCSV) WRITE(MYUNIT,'(A)') &
          'scan_all_swaps_parallel> No improvement found!'
  ENDIF
  !
  RETURN
  !
END SUBROUTINE SCAN_ALL_SWAPS_PARALLEL
