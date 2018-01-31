SUBROUTINE BENCH_SPARSE(COORDINATES, QUENCH_NUM, HESS_CUTOFF, POT_ENERGY, LOG_PROD)
   USE COMMONS, ONLY: NATOMS
   USE MODHESS
#ifdef __SPARSE
   USE MODSPARSEHESS
   USE SHIFT_HESS
#endif
   IMPLICIT NONE
! Arguments
   DOUBLE PRECISION, INTENT(IN)  :: COORDINATES(3*NATOMS)
   INTEGER, INTENT(IN)           :: QUENCH_NUM
   DOUBLE PRECISION, INTENT(IN)  :: HESS_CUTOFF
   DOUBLE PRECISION, INTENT(OUT) :: POT_ENERGY
   DOUBLE PRECISION, INTENT(OUT) :: LOG_PROD
! Timing variables
   DOUBLE PRECISION              :: TIME_DSYEV, TIME_FILT, TIME_SPARSE
   DOUBLE PRECISION              :: T_START, T_END
! File variables
   INTEGER                       :: SUMMARY_FILE
   INTEGER                       :: GETUNIT
! DSYEV variables
   INTEGER, PARAMETER            :: LWORK = 100000
   DOUBLE PRECISION              :: WORK(LWORK)
   INTEGER                       :: INFO
! Other variables
   DOUBLE PRECISION              :: LOG_PROD_DSYEV, LOG_PROD_FILT, LOG_PROD_SPARSE
   DOUBLE PRECISION              :: EVALUES(3*NATOMS)
   DOUBLE PRECISION              :: SHIFT_ARRAY(6)
   DOUBLE PRECISION              :: GRAD(3*NATOMS)
   DOUBLE PRECISION, ALLOCATABLE :: HESS_COPY(:, :)
   INTEGER                       :: NUM_ZERO_EVS 

#ifdef __SPARSE
! Assume 6 zero eigenvalues
   NUM_ZERO_EVS = 6

! Calculate hessian and mass weight
   IF (ALLOCATED(HESS)) DEALLOCATE(HESS)
   ALLOCATE(HESS(3*NATOMS, 3*NATOMS))
   CALL POTENTIAL(COORDINATES, GRAD, POT_ENERGY, .TRUE., .TRUE.)
   CALL MASSWT()

! Copy the original (unfiltered) hessian to somewhere safe
   IF (ALLOCATED(HESS_COPY)) DEALLOCATE(HESS_COPY)
   ALLOCATE(HESS_COPY(3*NATOMS, 3*NATOMS))
   HESS_COPY = HESS

! Assign shift array
   SHIFT_ARRAY(1:6) = 1.0D0

! Sparse section
   CALL CPU_TIME(T_START)
   CALL SHIFT_HESS_ZEROS(COORDINATES, SHIFT_ARRAY)
   CALL FILTER_ZEROS(HESS, HESS_CUTOFF)
   CALL GET_DETERMINANT(3*NATOMS, LOG_PROD_SPARSE, QUENCH_NUM)
   CALL CPU_TIME(T_END)
   TIME_SPARSE = T_END - T_START

! DSYEV without filters section
   HESS = HESS_COPY
   CALL CPU_TIME(T_START)
   CALL SHIFT_HESS_ZEROS(COORDINATES, SHIFT_ARRAY)
   CALL DSYEV('N', 'L', 3*NATOMS, HESS, 3*NATOMS, EVALUES, WORK, LWORK, INFO)
   CALL CPU_TIME(T_END)
   TIME_DSYEV = T_END - T_START
   LOG_PROD_DSYEV = SUM(DLOG(EVALUES(NUM_ZERO_EVS + 1:)))

! Need to pass back the log product of frequencies for actual FEBH (sqrt of eigenvalues)
   LOG_PROD = 0.5D0 * LOG_PROD_DSYEV

! Print summary
   SUMMARY_FILE = GETUNIT() 
   OPEN(UNIT=SUMMARY_FILE, FILE='sparse_benchmarks', STATUS='UNKNOWN', ACCESS='APPEND')
   WRITE(SUMMARY_FILE, '(A, I12)') 'Quench', QUENCH_NUM
   WRITE(SUMMARY_FILE, '(A, A, E18.10, A, F10.4, A, E18.10, A, F8.2)') 'DSYEV             ||', &
                                                                    ' log_prod= ', LOG_PROD_DSYEV, &
                                                                    ' time= ', TIME_DSYEV, &
                                                                    ' error= ', LOG_PROD_DSYEV - LOG_PROD_DSYEV, &
                                                                    ' % error= ', 100.0D0 * (LOG_PROD_DSYEV - LOG_PROD_DSYEV) / LOG_PROD_DSYEV
   WRITE(SUMMARY_FILE, '(A, A, E18.10, A, F10.4, A, E18.10, A, F8.2)') 'Sparse            ||', &
                                                                    ' log_prod= ', LOG_PROD_SPARSE, &
                                                                    ' time= ', TIME_SPARSE, &
                                                                    ' error= ', LOG_PROD_SPARSE - LOG_PROD_DSYEV, &
                                                                    ' % error= ', 100.0D0 * (LOG_PROD_SPARSE - LOG_PROD_DSYEV) / LOG_PROD_DSYEV
   CLOSE(SUMMARY_FILE) 
#endif 

END SUBROUTINE BENCH_SPARSE
