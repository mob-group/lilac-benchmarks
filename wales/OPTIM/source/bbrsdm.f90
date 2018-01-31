      SUBROUTINE BBRSDM(X,MFLAG,ITDONE,ENERGY,RMSGRD,NODUMP,G,PTEST) 

      USE COMMONS
      USE KEY
      USE MODCHARMM,ONLY : CHRMMT
    
      IMPLICIT NONE
      INTEGER          :: ITDONE, M, J, N
      DOUBLE PRECISION :: X(3*NATOMS), G(3*NATOMS), ENERGY
      DOUBLE PRECISION :: XNEW(3*NATOMS), GOLD(3*NATOMS), DELG(3*NATOMS), DELX(3*NATOMS), GNRMSQ, RMSGRD
      DOUBLE PRECISION :: ALPHA, GAM, LAMDA, EPS, EPSINV, EPS2, SIGMA1, SIGMA2, SIGMA, DELTA, A(BBRM), AMAX
      LOGICAL          :: GTEST, LOOP, MFLAG, NODUMP,PTEST
      CHARACTER ESTRING*87, GPSTRING*80, NSTRING*80, FSTRING*80
      COMMON /STRINGS/ ESTRING, GPSTRING, NSTRING, FSTRING

      NSTEPS = BBRSTEPS
      GAM    = BBRGAM
      EPS    = BBREPS
      SIGMA1 = BBRSIGMA1
      SIGMA2 = BBRSIGMA2
      M      = BBRM
      ALPHA  = BBRALPHA
      EPS2   = BBRCONV
      

      EPSINV = 1.D0/EPS
      N      = 3*NATOMS

      IF (CHRMMT) CALL UPDATENBONDS(X)

      ITDONE = 0

      GTEST = .TRUE.

      CALL POTENTIAL(X,ENERGY,G,GTEST,.FALSE.,RMSGRD,.FALSE.,.FALSE.)

      GNRMSQ = DOT_PRODUCT(G,G)
      RMSGRD = DSQRT(GNRMSQ/(1.0D0*N))

      DO WHILE (RMSGRD > EPS2)

         IF (PTEST) WRITE(*,'(A,2G20.10,A,I6,A)') &
     &        ' bbrsdm> Energy and RMS force=',ENERGY,RMSGRD,' after ',ITDONE,' BBRSDM steps'
         WRITE(ESTRING,16) 'bbrsdm> Energy for last cycle=',ENERGY,' '
16       FORMAT(A,27X,F20.10,A)


!     STEP 2

         IF (ALPHA <= EPS .OR. ALPHA >= EPSINV) THEN

            IF (GNRMSQ > 1.D0) THEN

               DELTA = 1.D0

            ELSEIF (GNRMSQ >= 1.D-05) THEN

               DELTA = 1.D0/GNRMSQ

            ELSE

               DELTA = 1.D05

            ENDIF
               
            ALPHA = DELTA

         ENDIF

!     STEP 3

         LAMDA = 1.D0/ALPHA

!     STEP 4

         J    = MOD(ITDONE,M)+1
         A(J) = ENERGY
         
         IF (ITDONE < M) THEN
            AMAX = MAXVAL(A(1:ITDONE))
         ELSE
            AMAX = MAXVAL(A)
         ENDIF

         LOOP = .TRUE.
         
         DO WHILE (LOOP)

            XNEW = X - LAMDA * G 

            GTEST = .FALSE.
            CALL POTENTIAL(XNEW,ENERGY,G,GTEST,.FALSE.,RMSGRD,.FALSE.,.FALSE.)

            IF (ENERGY <= (AMAX - GAM*LAMDA*GNRMSQ)) THEN

                X = XNEW
                LOOP = .FALSE.

            ELSE   ! STEP 5

               SIGMA = 0.5D0*(SIGMA1+SIGMA2)
               LAMDA = SIGMA*LAMDA

            ENDIF

         ENDDO

!     STEP 6

         DELG  = G - GOLD
         ALPHA = - DOT_PRODUCT(G,DELG)/(LAMDA*GNRMSQ)
         GOLD  = G

         GTEST = .TRUE.
         CALL POTENTIAL(X,ENERGY,G,GTEST,.FALSE.,RMSGRD,.FALSE.,.FALSE.) 
 
         GNRMSQ = DOT_PRODUCT(G,G)
         RMSGRD = DSQRT(GNRMSQ/(1.0D0*N))
          
         ITDONE = ITDONE + 1
         IF (.NOT.(NODUMP)) CALL DUMPP(X,ENERGY)

         IF (ITDONE == NSTEPS) EXIT 

      ENDDO

      MFLAG = .FALSE.
         
      END SUBROUTINE BBRSDM
