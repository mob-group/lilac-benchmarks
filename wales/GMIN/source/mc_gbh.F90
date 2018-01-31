!   GMIN: A program for finding global minima
!   Copyright (C) 1999-2006 David J. Wales
!   This file is part of GMIN and was writtedn by ds656.   
!
!   GMIN is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or   
!   (at your option) any later version.
!
!   GMIN is distributed in the hope that it will be useful,   
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!                                               
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software    
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
SUBROUTINE MC_GBH(NSTEPS)
  !
  USE PORFUNCS
  USE COMMONS, ONLY : MYNODE, MYUNIT, DEBUG, NATOMSALLOC, NQ, &
       TSTART, COORDS, LABELS, ECONV, BOXCENTROIDT,RANDMULTIPERMT, &
       NPAR_GBH, TEMP, RMS, TARGET, HIT, QALCST, QALCS_NBRHD, &
       QALCSV, GBH_RESTART, GBH_NREJMAX, GBH_NAVOID, &
       QALCS_SURFT, QALCS_SURF_MODE
  !
  IMPLICIT NONE
  !
#ifdef MPI 
  INCLUDE 'mpif.h'
  INTEGER MPIERR
#endif
  !
  INTEGER, INTENT(IN) :: NSTEPS
  !
  LOGICAL :: RESET, DUPLICATE
  INTEGER :: I, ITRIAL, ITERNS, BRUN, QDONE, NQTOT, IPROC, IPROCLO, &
       NREJSTREAK, L0(NATOMSALLOC), LDUM(NATOMSALLOC)
  DOUBLE PRECISION :: TIME, SCREENC(3*NATOMSALLOC), POTEL, R, &
       POTEL_LIST(NPAR_GBH+1), DPRAND, X0(3*NATOMSALLOC), E0, &
       XDUM(3*NATOMSALLOC), EDUM
  DOUBLE PRECISION, ALLOCATABLE :: AVOIDLIST(:)
  !
  COMMON /MYPOT/ POTEL
  !
  WRITE(MYUNIT, '(A)')  'mc_gbh> Calculating initial energy'
  !WRITE(MYUNIT, *)  'mc_gbh> NATOMSALLOC=', NATOMSALLOC
  CALL FLUSH(MYUNIT)
  CALL QUENCH(.FALSE.,1,ITERNS,TIME,BRUN,QDONE,SCREENC)
  NQ(1) = 0
  WRITE(MYUNIT,&
       '(A,I10,A,G20.10,A,I5,A,G12.5,A,F10.1)') &
       'Qu ',NQ(1),' E= ',POTEL,' steps= ',ITERNS,' RMS= ',RMS,&
       ' t= ',TIME-TSTART
  CALL FLUSH(MYUNIT)
  !
  IF(GBH_NAVOID > 0) THEN
     ALLOCATE(AVOIDLIST(GBH_NAVOID))
     AVOIDLIST(:) = 0.0D0
  ENDIF
  !
  WRITE(MYUNIT, '(A,I6)') &
       'mc_gbh> Starting GBH loop of length ',NSTEPS
  !
  X0(1:3*NATOMSALLOC) = COORDS(1:3*NATOMSALLOC,1)
  L0(1:NATOMSALLOC) = LABELS(1:NATOMSALLOC,1)
  E0 = POTEL
  !
  NREJSTREAK=0
  gbh_loop: DO ITRIAL=1,NSTEPS
     !
     WRITE(MYUNIT, '(A,I10)') &
          'mc_gbh> Start of iteration ', ITRIAL     
     !
     RESET=.FALSE.
     IF(GBH_RESTART) THEN
        IF(NREJSTREAK >= GBH_NREJMAX) RESET=.TRUE.
        IF(GBH_NAVOID>0.AND.(RESET.OR.NREJSTREAK==1)) THEN
           CALL PROCESS_AVOID_LIST(E0,GBH_NAVOID,AVOIDLIST,RESET)
        ENDIF
        IF(RESET) THEN
           E0=1.0D+99 ! Guarantee move by inflating Markov energy
           NREJSTREAK=0           
        ENDIF
     ENDIF
     !
     ! --- Stochastic cartesian move -------
     CALL RANDOM_MOVE(RESET,3,NATOMSALLOC,COORDS(:,1))     
     !
     !IF (RANDMULTIPERMT.AND.MOD(J1,RANDMULTIPERM_STEP)==0) &
     !     CALL RANDMULTIPERM(1) ! re-write this routine!!!!!
     !IF(BOXCENTROIDT) CALL BOXCENTROID(COORDS(:,1))
     !
     ! Quench perturbed state
     NQ(1) = NQ(1) + 1
     CALL QUENCH(.FALSE.,1,ITERNS,TIME,BRUN,QDONE,SCREENC)
     WRITE(MYUNIT,&
          '(A,I10,A,G20.10,A,I5,A,G12.5,A,G20.10,A,F10.1)') &
          'Qu ',NQ(1),' E= ',POTEL,' steps= ',ITERNS,' RMS= ',&
          RMS,' E0=',E0,' t= ',TIME-TSTART
     CALL FLUSH(MYUNIT)
     !
#ifdef MPI
     IF (QALCS_SURF_MODE == -1) THEN
        !
        ! Back-up randomly perturbd state on each PROC
        XDUM(1:3*NATOMSALLOC) = COORDS(1:3*NATOMSALLOC,1)
        LDUM(1:NATOMSALLOC) = LABELS(1:NATOMSALLOC,1)
        EDUM = POTEL
        !
        proc_loop: DO IPROC=0,NPAR_GBH-1
           !
           WRITE(MYUNIT,'(A,I3)') &
                'mc_gbh> Refining state from IPROC= ',IPROC+1
           !
           ! Reinstate the random state on IPROC
           IF(MYNODE == IPROC) THEN
              COORDS(1:3*NATOMSALLOC,1) = XDUM(1:3*NATOMSALLOC)
              LABELS(1:NATOMSALLOC,1) = LDUM(1:NATOMSALLOC)
              POTEL = EDUM
           ENDIF
           !
           ! Broadcast state from IPROC
           CALL MPI_BCAST(COORDS(:,1),3*NATOMSALLOC,&
                MPI_DOUBLE_PRECISION,IPROC,MPI_COMM_WORLD,MPIERR)
           CALL MPI_BCAST(LABELS(:,1),NATOMSALLOC,MPI_INTEGER,&
                IPROC,MPI_COMM_WORLD,MPIERR)
           CALL MPI_BCAST(POTEL,1,MPI_DOUBLE_PRECISION,&
                IPROC,MPI_COMM_WORLD,MPIERR)
           !
           ! Accumulate POTEL_LIST with randomly perturbed 
           ! energies and check for duplicates
           POTEL_LIST(IPROC+1) = POTEL
           DUPLICATE=.FALSE.
           IF(IPROC>0) THEN
              DUPLICATE=.FALSE.
              check_dupes: DO I=IPROC,1,-1
                 IF(ABS(POTEL - POTEL_LIST(I)) < ECONV) THEN
                    DUPLICATE=.TRUE.
                    EXIT check_dupes
                 ENDIF
              ENDDO check_dupes
           ENDIF
           !
           IF(.NOT.DUPLICATE .AND. &
                ABS(POTEL - E0) > ECONV) THEN
              ! Do surface refinement (parallelised!)
              CALL QALCS_SURF(1, ITERNS, TIME, BRUN, QDONE, SCREENC)
              ! Update back-up on IPROC
              IF(MYNODE == IPROC) THEN
                 XDUM(1:3*NATOMSALLOC) = COORDS(1:3*NATOMSALLOC,1)
                 LDUM(1:NATOMSALLOC) = LABELS(1:NATOMSALLOC,1)
                 EDUM = POTEL
              ENDIF
           ELSE ! Just inflate energy to prohibit selection
              WRITE(MYUNIT,'(A)') &
                   'mc_gbh> Duplicate skipped.'
              IF(MYNODE == IPROC) EDUM = 1.0D+99 
           ENDIF
           !
        ENDDO proc_loop
        !
        ! Reset to back-up on each PROC
        COORDS(1:3*NATOMSALLOC,1) = XDUM(1:3*NATOMSALLOC)
        LABELS(1:NATOMSALLOC,1) = LDUM(1:NATOMSALLOC)
        POTEL = EDUM
        !
     ENDIF
     !
     ! Gather all parallel quench energies on master node.
     CALL MPI_GATHER(POTEL, 1, MPI_DOUBLE_PRECISION, &
          POTEL_LIST(1:NPAR_GBH), 1, MPI_DOUBLE_PRECISION, 0, &
          MPI_COMM_WORLD, MPIERR)
#else
     POTEL_LIST(1) = POTEL
#endif
     !
     IF(MYNODE==0) THEN
        !
        POTEL_LIST(NPAR_GBH+1) = E0 ! Allows rejection
        IF (TARGET) IPROCLO=MINLOC(POTEL_LIST(1:NPAR_GBH), 1)-1
        !
        ! Choose an IPROC.
        CALL CHOOSE_FROM_LIST(NPAR_GBH+1,POTEL_LIST,IPROC)
        !
        IF(IPROC>NPAR_GBH) THEN ! Reject all
           COORDS(1:3*NATOMSALLOC,1) = X0(1:3*NATOMSALLOC)
           LABELS(1:NATOMSALLOC,1) = L0(1:NATOMSALLOC) 
           POTEL = E0
           IPROC = 0
           WRITE(MYUNIT, '(A)') 'mc_gbh> Rejected trial(s)'
        ELSE
           WRITE(MYUNIT, '(A,I3)') &
                'mc_gbh> Accepted IPROC= ', IPROC    
           IPROC=IPROC-1
        ENDIF
        !
     ENDIF
     !
#ifdef MPI 
     ! Broadcast IPROC from master
     CALL MPI_BCAST(IPROC,1,MPI_INTEGER,&
          0,MPI_COMM_WORLD,MPIERR)
     ! Broadcast state from IPROC
     CALL MPI_BCAST(COORDS(:,1),3*NATOMSALLOC,&
          MPI_DOUBLE_PRECISION,IPROC,MPI_COMM_WORLD,MPIERR)
     CALL MPI_BCAST(LABELS(:,1),NATOMSALLOC,MPI_INTEGER,&
          IPROC,MPI_COMM_WORLD,MPIERR)
     CALL MPI_BCAST(POTEL,1,MPI_DOUBLE_PRECISION,&
          IPROC,MPI_COMM_WORLD,MPIERR)
     IF(QALCSV) WRITE(MYUNIT, '(A,F20.10,A,I6)') &
          'mc_gbh> Use E= ', POTEL,' from node ',IPROC+1
     ! Broadcast HIT from IPROCLO if targetting...
     IF(TARGET) THEN
        ! Broadcast IPROCLO from master
        CALL MPI_BCAST(IPROCLO,1,MPI_INTEGER,&
             0,MPI_COMM_WORLD,MPIERR)
        CALL MPI_BCAST(HIT,1,MPI_LOGICAL,&
             IPROCLO,MPI_COMM_WORLD,MPIERR)        
     ENDIF
#endif                     
     !
     !IF(POTEL < E0 - ECONV) THEN
     IF(ABS(POTEL-E0) > ECONV) THEN
        NREJSTREAK=0
     ELSE
        NREJSTREAK = NREJSTREAK+1
     ENDIF
     !
     IF(HIT) THEN
        WRITE(MYUNIT,'(A,I3,A,I6)') &
             'mc_gbh> Target hit stochastically by node ',IPROCLO+1,&
             ' on trial ',ITRIAL
        EXIT gbh_loop
     ENDIF
     !
     ! --- Deterministic refinement --------------------------
     IF(.NOT.HIT .AND. ABS(POTEL-E0)>ECONV) THEN
        IF(QALCS_SURFT .AND. QALCS_SURF_MODE /= -1) THEN
           CALL QALCS_SURF(1, ITERNS, TIME, BRUN, QDONE, SCREENC)
           WRITE(MYUNIT,'(A,G20.10,A,F11.1)') &
           'mc_gbh> Refined E= ',POTEL,' t= ', TIME-TSTART
        ENDIF
        IF( QALCST .AND. QALCS_NBRHD>0) THEN
           ! perform a variable neighbourhood search (in parallel)
           CALL QALCS_PARALLEL(ITERNS, TIME, BRUN, QDONE, SCREENC)
           WRITE(MYUNIT,'(A,G20.10,A,F11.1)') &
           'mc_gbh> Biminimum E= ',POTEL,' t= ', TIME-TSTART
        ENDIF
     ENDIF
     !
     IF(HIT) EXIT gbh_loop
     !
     X0(1:3*NATOMSALLOC) = COORDS(1:3*NATOMSALLOC,1)
     L0(1:NATOMSALLOC) = LABELS(1:NATOMSALLOC,1)
     E0 = POTEL
     !
  ENDDO gbh_loop
  !
  IF(ALLOCATED(AVOIDLIST)) DEALLOCATE(AVOIDLIST)
  !
  RETURN
  !
END SUBROUTINE MC_GBH

SUBROUTINE CHOOSE_FROM_LIST(N,VALUES,I)
  !
  ! Choose an element of VALUES with Boltzmann probability
  !
  USE COMMONS, ONLY : TEMP, QALCSV, MYUNIT, ECONV
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: N
  DOUBLE PRECISION, INTENT(IN) :: VALUES(N)
  INTEGER, INTENT(OUT) :: I
  !
  LOGICAL :: NOREJECT
  DOUBLE PRECISION :: PSUM(N), X, Y, DPRAND, ELOWEST
  !
  !X=0.0D0 ! initialise total sum
  !DO I=1,N-1
  !   Y = MIN(DEXP(-(VALUES(I)-VALUES(N))/TEMP(1)), 1.0D0)
  !   X = X + Y
  !   PSUM(I) = X ! store partial sum
  !ENDDO
  !PSUM(N) = DBLE(N-1)
  !
  ELOWEST=MINVAL(VALUES)
  IF(ABS(ELOWEST - VALUES(N)) < ECONV) THEN
     NOREJECT = .FALSE.
  ELSE
     NOREJECT = .TRUE.
  ENDIF
  !
  X=0.0D0 ! initialise total sum
  DO I=1,N-1
     IF( (NOREJECT .AND. VALUES(I) > VALUES(N)+ECONV) .OR. &
          ABS(VALUES(I) - VALUES(N)) < ECONV) THEN
        Y=0.0D0
     ELSE
        Y=DEXP(-(VALUES(I)-ELOWEST)/TEMP(1))
     ENDIF
     X = X + Y
     PSUM(I) = X ! store partial sum
  ENDDO
  IF(VALUES(N) > ELOWEST + ECONV) THEN
     PSUM(N) = PSUM(N-1) ! Disallow rejection
  ELSE
     PSUM(N) = DBLE(N-1) ! Allow rejection
  ENDIF
  !
  X=PSUM(N)*DPRAND()
  !
  I=1
  DO WHILE (X > PSUM(I))
     I=I+1
  ENDDO
  !
  !IF(QALCSV) THEN
  !   WRITE(MYUNIT, *) 'choose_from_list> vals=', VALUES
  !   WRITE(MYUNIT, *) 'choose_from_list> psum=', PSUM
  !   WRITE(MYUNIT, *) 'choose_from_list> choice=', I
  !ENDIF
  !
  RETURN
  !
END SUBROUTINE CHOOSE_FROM_LIST

SUBROUTINE RANDOM_MOVE(RESET,DIM,N,X)
  !
  USE COMMONS, ONLY : STEP, ASTEP, MYUNIT, MIEFT
  !
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: RESET
  INTEGER, INTENT(IN) :: DIM, N
  DOUBLE PRECISION, INTENT(INOUT) :: X(DIM*N)
  !
  INTEGER :: I, J, K, N1, N2
  DOUBLE PRECISION :: CNTR(DIM), RMIN, RMAX, RALL(N), R, &
       THETA, PHI, DPRAND, TWOPI
  !
  DOUBLE PRECISION, PARAMETER :: PI=3.141592654D0
  TWOPI=2.0D0*PI
  !
  CALL FIND_CENTROID(DIM,N,X,CNTR)
  !
  RMAX=0.0D0
  DO I=1,N
     J=DIM*(I-1)
     R=0.0D0
     DO K=1,DIM
        R = R + (X(I+K)-CNTR(K))**2
     ENDDO
     RALL(I) = DSQRT(R)
     IF(RMAX < RALL(I)) RMAX = RALL(I)
  ENDDO
  !
  R=(DPRAND()-0.5D0)*2.0D0 ! 50% chance of no angular move
  RMIN = RMAX - R*ASTEP(1)
  N1=0; N2=0
  !
  DO I=1,N
     J=DIM*(I-1)
     THETA = ACOS(2*DPRAND()-1) !THETA = DPRAND()*PI
     PHI = DPRAND()*TWOPI
     IF(RESET) THEN ! Reset to random coordinates
        !        
        R=RMAX*DPRAND()**(1.0D0/3.0D0) !R = DPRAND()*RMAX
        X(J+1) = CNTR(1) + R*DSIN(THETA)*DCOS(PHI)
        X(J+2) = CNTR(2) + R*DSIN(THETA)*DSIN(PHI)
        IF(DIM>2) THEN
           IF(MIEFT) THEN !<ds656 specific to z=0 substrate!!!
              X(J+3) = RMAX + R*DCOS(THETA)
           ELSE
              X(J+3) = CNTR(3) + R*DCOS(THETA)
           ENDIF
        ENDIF
        !
     ELSEIF(RALL(I) > RMIN) THEN ! Surface angular move
        !
        X(J+1) = CNTR(1) + RMAX*DSIN(THETA)*DCOS(PHI)
        X(J+2) = CNTR(2) + RMAX*DSIN(THETA)*DSIN(PHI)
        IF(DIM>2) X(J+3) = CNTR(3) + RMAX*DCOS(THETA)
        N2=N2+1
        !
     ELSE ! Interior shake move
        !
        R=STEP(1)*DPRAND()**(1.0D0/3.0D0) !R = DPRAND()*RMAX
        X(J+1) = X(J+1) + R*DSIN(THETA)*DCOS(PHI)
        X(J+2) = X(J+2) + R*DSIN(THETA)*DSIN(PHI)
        IF(DIM>2) X(J+3) = X(J+3) + R*DCOS(THETA)
        !
        ! old way:
        !DO K=1,DIM
           !R=(DPRAND()-0.5D0)*2.0D0
           !X(J+K)=X(J+K) + STEP(1)*R
        !ENDDO
        N1=N1+1
        !
     ENDIF
  ENDDO
  !
  IF(RESET) THEN
     WRITE(MYUNIT,'(A)') 'mc_gbh> Configuration reset!'
  !ELSEIF(N2>0) THEN
  !   WRITE(MYUNIT, '(A,I6)') 'mc_gbh> Atoms rotated: ', N2
  ENDIF
  !
  RETURN
  !
END SUBROUTINE RANDOM_MOVE
!
SUBROUTINE FIND_CENTROID(DIM,N,X,CNTR)
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: DIM, N
  DOUBLE PRECISION, INTENT(IN) :: X(DIM*N)
  DOUBLE PRECISION, INTENT(OUT) :: CNTR(DIM)
  !
  INTEGER :: I, J, K
  !
  CNTR(:) = 0.0D0
  DO I=1,N
     J=3*(I-1)
     DO K=1,DIM
        CNTR(K) = CNTR(K) + X(J+K)
     ENDDO
  ENDDO
  CNTR(:) = CNTR(:)/DBLE(N)
  !
  RETURN
  !
END SUBROUTINE FIND_CENTROID

SUBROUTINE PROCESS_AVOID_LIST(E,N,ELIST,SWITCH)
  !
  USE COMMONS, ONLY : ECONV, MYUNIT
  !
  IMPLICIT NONE
  !
  DOUBLE PRECISION, INTENT(IN) :: E
  INTEGER, INTENT(IN) :: N
  DOUBLE PRECISION, INTENT(INOUT) :: ELIST(N)
  LOGICAL, INTENT(INOUT) :: SWITCH
  !
  INTEGER :: I, J
  DOUBLE PRECISION :: EDUM
  !
  IF(SWITCH) THEN ! Append E to end of ELIST and then shuffle down.
     !
     ELIST(N) = E
     DO I=N-1,1,-1
        EDUM = ELIST(I)
        ELIST(I) = ELIST(I+1)
        ELIST(I+1) = EDUM
     ENDDO
     WRITE(MYUNIT,'(A,G20.10)') &
          'mc_gbh> Added to AVOIDLIST E= ',E 
     !
  ELSE ! Check if E is present in ELIST
     !
     elistscan: DO I=1,N
        IF(ABS(E-ELIST(I)) < ECONV) THEN
           SWITCH=.TRUE.
           WRITE(MYUNIT,'(A,G20.10)') &
                'mc_gbh> Found in AVOIDLIST E= ',E
           ! Shuffle match down to 1st place
           DO J=I,2,-1
              EDUM = ELIST(J)
              ELIST(J) = ELIST(J-1)
              ELIST(J-1) = EDUM
           ENDDO
           EXIT elistscan
        ENDIF
     ENDDO elistscan
     !
  ENDIF
  !
  RETURN
  !
END SUBROUTINE PROCESS_AVOID_LIST
