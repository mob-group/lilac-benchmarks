SUBROUTINE MIEF_INI()
  !
  USE KEY, ONLY : MIEF_FILENAME, MIEF_EPS, MIEF_SIG, &
       MIEF_NSITES, MIEF_N, MIEF_M, MIEF_SITES, MYUNIT, &
       MIEF_CUTT, MIEF_RCUT, MIEF_U_RCUT, MIEF_DUDR_RCUT
  USE COMMONS, ONLY: NATOMS
  !
  IMPLICIT NONE
  !
  LOGICAL :: YESNO
  INTEGER :: I, J, EOF, LUNIT, GETUNIT, NSPECIES, NTYPEA
  DOUBLE PRECISION :: ATT_TERM, REP_TERM, DUMMY, PREF
  !
  DOUBLE PRECISION :: EPSAB,EPSBB,SIGAB,SIGBB ! Fucking common block
  COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB,NTYPEA
  !
  IF(NTYPEA < NATOMS) THEN
     NSPECIES = 2
  ELSE
     NSPECIES = 1
  ENDIF
  ALLOCATE(MIEF_SIG(NSPECIES))
  ALLOCATE(MIEF_EPS(NSPECIES))
  !
  INQUIRE(FILE=MIEF_FILENAME, EXIST=YESNO)
  !
  IF(YESNO) THEN
     LUNIT=GETUNIT()
     OPEN (LUNIT, FILE=MIEF_FILENAME, STATUS='OLD')
     READ(LUNIT,*) MIEF_NSITES, MIEF_N, MIEF_M, &
          (MIEF_EPS(I), I=1,NSPECIES), &
          (MIEF_SIG(I), I=1,NSPECIES)
     !
     ALLOCATE(MIEF_SITES(MIEF_NSITES,3))
     !
     WRITE(*,'(A,I2,A,I1,A,I3,A)') 'mief_ini> ',&
          MIEF_N,'-',MIEF_M,' Mie surface with ', &
          MIEF_NSITES,' sites.'
     WRITE(*,'(A)', ADVANCE='NO') 'mief_ini> Epsilon value(s): '
     DO I=1,NSPECIES
        WRITE(*,'(F8.4)', ADVANCE='NO') MIEF_EPS(I)
     ENDDO
     WRITE(*,*)
     WRITE(*,'(A)', ADVANCE='NO') 'mief_ini> Sigma value(s): '
     DO I=1,NSPECIES
        WRITE(*,'(F8.4)', ADVANCE='NO') MIEF_SIG(I)
     ENDDO
     WRITE(*,*)
     !
     DO I=1,MIEF_NSITES
        READ(LUNIT,*,IOSTAT=EOF) ( MIEF_SITES(I,J), J=1,3 )
        IF(EOF.LT.0) THEN
           WRITE(*,*) "mief_ini> Premature end of file ", &
                TRIM(MIEF_FILENAME)
           STOP
        ENDIF
     ENDDO
     !
     CLOSE(LUNIT)
     !
  ELSE
     WRITE(*,*) "mief_ini> Missing file ", TRIM(MIEF_FILENAME)
     STOP
  ENDIF
  !
  IF (MIEF_CUTT) THEN
     ALLOCATE(MIEF_U_RCUT(NSPECIES))
     ALLOCATE(MIEF_DUDR_RCUT(NSPECIES))
     PREF = DBLE(MIEF_N)/DBLE(MIEF_N-MIEF_M)*&
          (DBLE(MIEF_N)/DBLE(MIEF_M))**&
          (DBLE(MIEF_M)/DBLE(MIEF_N-MIEF_M))
     DO I=1,NSPECIES
        DUMMY = MIEF_SIG(I)/MIEF_RCUT
        REP_TERM=DUMMY**MIEF_N
        ATT_TERM=DUMMY**MIEF_M
        MIEF_U_RCUT(I)=PREF*MIEF_EPS(I)*(REP_TERM-ATT_TERM)
        MIEF_DUDR_RCUT(I)=PREF*MIEF_EPS(I)*( -DBLE(MIEF_N)*REP_TERM + &
             DBLE(MIEF_M)*ATT_TERM ) / MIEF_RCUT
     ENDDO
  ENDIF
  !
  RETURN
  !
END SUBROUTINE MIEF_INI
!
SUBROUTINE MIEF(NATOMS,X,GRAD,EREAL,GRADT,HESST)
  !
  USE KEY, ONLY : MIEF_N, MIEF_M, MIEF_NSITES, MIEF_SITES, &
       MIEF_EPS, MIEF_SIG, MIEF_PBCT, MIEF_BOX, MIEF_CUTT, &
       MIEF_RCUT, MIEF_U_RCUT, MIEF_DUDR_RCUT
  USE MODHESS, ONLY : HESS
  !
  IMPLICIT NONE
  !
  INTEGER, INTENT(IN) :: NATOMS
  DOUBLE PRECISION, INTENT(IN) :: X(3*NATOMS)
  DOUBLE PRECISION, INTENT(INOUT) :: GRAD(3*NATOMS)
  DOUBLE PRECISION, INTENT(INOUT) :: EREAL
  LOGICAL, INTENT(IN) :: GRADT,HESST
  !
  INTEGER :: I, J, K, L, I3, NTYPEA, ITYPE
  DOUBLE PRECISION :: PREF, DX(3),DIST2,IDIST,SIG,EPS,&
       DUMMY, ATT_TERM, REP_TERM, DVDR, D2VDR2, DIST
  !
  DOUBLE PRECISION :: EPSAB,EPSBB,SIGAB,SIGBB ! Fucking common block
  COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB,NTYPEA
  !
  PREF = DBLE(MIEF_N)/DBLE(MIEF_N-MIEF_M)*&
       (DBLE(MIEF_N)/DBLE(MIEF_M))**&
       (DBLE(MIEF_M)/DBLE(MIEF_N-MIEF_M))
  !
  DO I=1,NATOMS ! Loop over atoms
     !
     I3=3*(I-1)
     IF(I > NTYPEA) THEN
        ITYPE=2
        EPS=MIEF_EPS(ITYPE)*PREF
        SIG=MIEF_SIG(ITYPE)
     ELSE
        ITYPE=1
        EPS=MIEF_EPS(ITYPE)*PREF
        SIG=MIEF_SIG(ITYPE)
     ENDIF
     !
     DO J=1,MIEF_NSITES ! Loop over Mie sites
        !
        DIST2=0.0
        DO K=1,3
           DX(K) = X(I3+K) - MIEF_SITES(J,K)
           IF (MIEF_PBCT) DX(K) = DX(K) - MIEF_BOX(K)*&
                NINT(DX(K)/MIEF_BOX(K))
           DIST2 = DIST2 + DX(K)*DX(K)
        ENDDO
        DIST = DSQRT(DIST2)
        !
        IF (DIST < MIEF_RCUT) THEN
           !
           IDIST = SIG/DIST
           REP_TERM=IDIST**MIEF_N
           ATT_TERM=IDIST**MIEF_M
           !
           DUMMY = EPS*(REP_TERM-ATT_TERM)
           IF(MIEF_CUTT) THEN
              DUMMY = DUMMY - MIEF_U_RCUT(ITYPE) - &
                   (DIST-MIEF_RCUT)*MIEF_DUDR_RCUT(ITYPE)
           ENDIF
           EREAL = EREAL + DUMMY
           !
           IF(GRADT) THEN
              ! Compute 1/R * dV/dR 
              DVDR = EPS*( -DBLE(MIEF_N)*REP_TERM + &
                   DBLE(MIEF_M)*ATT_TERM ) / DIST2
              IF(MIEF_CUTT) DVDR = DVDR - MIEF_DUDR_RCUT(ITYPE)/DIST
              DO K=1,3
                 GRAD(I3+K) = GRAD(I3+K) + DVDR*DX(K)
              ENDDO
              !
              IF(HESST) THEN
                 !
                 ! Compute d2V/dR2 
                 D2VDR2 = EPS*( DBLE(MIEF_N*(MIEF_N+1))*REP_TERM - &
                      DBLE(MIEF_M*(MIEF_M+1))*ATT_TERM ) / DIST2
                 ! What we actually need is (D2VDR - DVDR/DIST) / DIST2   
                 D2VDR2 = (D2VDR2 - DVDR) / DIST2
                 !
                 DO K = 1,3 ! 1st coordinate of atom I
                    DO L = K,3 ! 2nd coordinate of atom I
                       !
                       DUMMY = DX(K)*DX(L) ! Reset DUMMY
                       DUMMY = DUMMY*D2VDR2 
                       IF(K==L) DUMMY = DUMMY + DVDR
                       !
                       ! Accumulate diagonal 3-blocks only
                       HESS(I3+K, I3+L) = HESS(I3+K, I3+L) + DUMMY
                       !
                    ENDDO
                 ENDDO
                 !
              ENDIF
              !
           ENDIF
           !
        ENDIF ! Cutoff check
        !
     ENDDO ! Spanned sites
     !
     IF(HESST) THEN
        ! Use symmetry to fill skipped entries in diagonal blocks
        DO K=1,3
           DO L=K+1,3
              HESS(I3+L,I3+K) = HESS(I3+K, I3+L) 
           ENDDO
        ENDDO
     END IF
     !
  ENDDO ! Spanned atoms
  !
  RETURN
  !
END SUBROUTINE MIEF
