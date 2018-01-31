!
! svn info 
! $Date: 2011-10-16 20:53:28 +0100 (Sun, 16 Oct 2011) $
! $Rev: 61 $
!

!     ------------------------------------------------------------------------------------

      MODULE MDCOMMONS

      IMPLICIT NONE

      SAVE

! N = max number of particles? used for setting array sizes
! NN = max number of total neighbors.  used for setting array sizes
! NMOL = number of molecules
! NSITE = number of sites per molecule?  normally 3.  still used, but
!         sometimes assumed to be three
! LIST = list of lists of all neighbors of each molecule.
! POINT = holds the starting indices for each molecule for navigating
!         LIST.  e.g.  LIST( POINT(I) ) is the first neighbor of
!         molecule I
! MASS = total mass of molecule
! MST = holds the XYZ coordinates of the 3 beads of the molecule
!       relative to the center of mass
! R0  = hold the xyz coords of the center of mass of the molecules at
!       last save
! JMI = moment of inertia of the molecule
! R = array of xyz coords of the center of mass of the molecules
! QTRN = array of quaternions defining the angular orientation of the
!        molecules
! P = the conjugate momentum associated with QTRN
! V = velocity of the CoM of the molecules
!     1/MASS*(the momentum conjugate to R)
! W = angular velocity of the molecules body-frame
! F = forces on the molecules 
! T = torques on the molecules
! PEP = potential energy of the molecules?
! DELT = step size?
!
! BOXL = box length
! KEP = twice? the kinetic energy of each particle ?
! PEP = potential energy of each particle
! SUMEPT = sum of total energy over all time steps times DELT
! RCUT  = RCUT is the cutoff in the Lennard Jones potential.
! RCUTSQ = RCUT*RCUT
! RCUT2 = Maximum molecule to molecule separation at which to calculate
!         forces.  WARNING, this is measured from
!         CoM to CoM, so it will exclude sites which are roughly RCUT2-1.33 appart.
! RCUT2SQ = RCUT2*RCUT2
! RLIST = Maximum molecule to molecule separation at which to build up neighbor
!         list
! RLSTSQ = RLIST**2
! UPDATE = if true then update neighbor list
! PHI = an array keeping track of the angular deviation of each molecule
!       from it's initial position (at time t=NEQ).  It is calulated as
!       a sum over time of W*DELT, so it is slightly approximate.
!
! S = S is the additional degree of freedom in the Nose-Poincare thermostat
! PS = the conjugate momentum associated with S
! HNOT = the value of the hamiltonitan at time t=0.  Could also have been called
!        HNAUGHT or H0  
      INTEGER, PARAMETER :: N = 1000, NN = N * 250
      INTEGER          :: NMOL, NSITE=2, NTST, LIST(NN), POINT(N), G
      DOUBLE PRECISION :: R(N,3), QTRN(N,4), V(N,3), W(N,3), F(N,3), T(N,3), T4(N,4), P(N,4), PHI(N,3) 
      DOUBLE PRECISION :: MST(3,3), R0(N,3) 
      DOUBLE PRECISION :: KEP(N), PEP(N), SUMEPT(N)
      DOUBLE PRECISION :: BOXL, RCUT, RCUTSQ, RCUT2, RCUT2SQ, EPS4, CNSTA, CNSTB, DELT, HALFDT, MASS, JMI(3)
      DOUBLE PRECISION :: RLIST, RLSTSQ, TMPFIX
      LOGICAL          :: UPDATE
      !these are used only in the NVT RUN
      DOUBLE PRECISION :: FCTMQT=1, S=1, PS=1, HNOT=1

      END MODULE MDCOMMONS

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE FORQUE(PE, VRL)

!     CALCULATES FORCES AND TORQUES on the molecules FROM SITE-SITE POTENTIAL 

      USE MDCOMMONS
      use py_module
      use commons, only : pyt, boxlx, boxly, boxlz
      IMPLICIT NONE
 
      INTEGER          :: I, J, INDX, K, JBEG, JEND, JNEB, NLIST 
      DOUBLE PRECISION :: RSITE(NTST,3), FSITE(NTST,3)
      !DOUBLE PRECISION :: Q(4), Q2Q3, Q1Q4, Q2Q4, Q1Q3, Q3Q4, Q1Q2, RMX(3,3)
      !DOUBLE PRECISION :: DR(3), DSS(3), DCMCM, DSS2, R2, R6, R12!, SHIFT(3)
      DOUBLE PRECISION :: DR(3), DCMCM 
      DOUBLE PRECISION :: PE, VRL, AA(3), PERIODIC_OFFSET(NMOL,3)
      double precision :: energy_contrib, grad_contrib(12),INDX1,INDX2,t_rep(6),t_att(6)
      type(py_molecule), pointer :: moli, molj
      type(py_site), pointer :: sitei, sitej
      integer :: moli_index, molj_index, sitei_index, sitej_index, j1, j2, FRAMENUM

!     INITIALIZE

      PE         = 0.D0 ! total potential energy
      VRL        = 0.D0 ! accumulates force*position. used to calulcate pressure
      PEP(:)     = 0.D0 ! potential energy of the molecules
      F(:,:)     = 0.D0 ! force on the molecules in space-fixed frame
      FSITE(:,:) = 0.D0 ! force on each site in space-fixed frame
      T(:,:)     = 0.D0 ! torque on the molecules
      FRAMENUM = 1
      PERIODIC_OFFSET(:,:) = 0.0D0      
!     OBTAIN THE SITE POSITIONS IN THE SPACE-FIXED FRAME FROM THE
!     CENTRE-OF-MASS POSITION & THE SITE POSITIONS IN THE BODY-FIXED
!     FRAME USING THE ROTATION MATRIX

!      CALL GET_RSITE( RSITE )

! get orientation from quaternions to angle-axis
!        WRITE(*,*) 'before apply_periodic'
!        WRITE(*,*) r(:,:)

!      CALL APPLY_PERIODIC( PERIODIC_OFFSET, 1 )
!        WRITE(*,*) 'after apply_periodic'
!        WRITE(*,*) r(:,:)

    ! update the values for all the molecules
    do moli_index = 1, nmol
        CALL QTRN2AA( QTRN(moli_index,:), AA )
        moli => molecules(moli_index)
        call update_py_molecule(moli, r(moli_index,:), &
            & AA, .true.)
!        WRITE(*,*) moli_index
    end do
IF (UPDATE) THEN  !if update neighbor list

!     SAVE CURRENT POSITIONS

         CALL SVCNFG()

!     CONSTRUCT NEIGHBOUR LIST AND CALCULATE FORCES AND TORQUES

         NLIST = 0

!         DO I = 1, NMOL - 1

!            POINT(I) = NLIST + 1

!            DO J =  I + 1, NMOL

  do moli_index = 1, nmol - 1
        moli => molecules(moli_index)
        I = moli_index
             POINT(I) = NLIST + 1
       
        ! inner loop over molecules
        do molj_index = moli_index + 1, nmol
            molj => molecules(molj_index)
            J = molj_index
               ! CoM separation 
               DR(:)  = R(I,:) - R(J,:) - NINT((R(I,:) - R(J,:))/BOXLX)*BOXLX
               DCMCM  = DOT_PRODUCT(DR,DR)

!               IF (DCMCM < RLSTSQ ) THEN

                  NLIST       = NLIST + 1
                  LIST(NLIST) = J

!       REMOVE THIS CHECK IF MAXNAB IS APPROPRIATE

                  IF (NLIST == NN) STOP 'list too small'

                IF(PYT) THEN
                    ! loop over sites in outer molecule
                  do sitei_index = 1, size(moli%sites)
                      sitei => moli%sites(sitei_index)

                      ! loop over sites in inner molecule
                      do sitej_index = 1, size(molj%sites)
                          sitej => molj%sites(sitej_index)

                          energy_contrib = 0.d0
                          grad_contrib(:) = 0.d0

                          ! compute the energy and gradient for this pair of sites
                          call pairwise_py(sitei, sitej, grad_contrib, energy_contrib, .true., t_rep, t_att)
                          ! add the energy to total 
                          PE = PE + energy_contrib
                          PEP(I) = PEP(I) + energy_contrib
                          PEP(J) = PEP(J) + energy_contrib
                           DO K = 1, 3 ! accumulate forces and torques
                                F(I,k) = F(I,k) - grad_contrib(k)
                                F(J,k) = F(J,k) - grad_contrib(k+3)
                                T(I,k) = T(I,k) - t_rep(k) - t_att(k) 
                                T(J,k) = T(J,k) - t_rep(3+k) - t_att(3+k)
                           END DO

                      end do
                  end do

                ELSE ! LWOTP model
                 CALL MOLMOL_FORQUE(I, J, RSITE, FSITE, DR, DCMCM, PE, VRL)
                END IF
!               ENDIF

            ENDDO
                   
         ENDDO 

        POINT(NMOL) = NLIST + 1 
!        WRITE(*,*) NMOL, POINT(NMOL), NLIST
ELSE ! if .not. UPDATE use current neighbor list to find forces

  do moli_index = 1, nmol - 1
        moli => molecules(moli_index)
        I = moli_index
!         DO I = 1, NMOL - 1

!       USE THE LIST TO FIND THE NEIGHBOURS    

            JBEG = POINT(I)
            JEND = POINT(I+1) - 1

!       CHECK THAT MOLECULE I HAS NEIGHBOURS

      IF (JBEG <= JEND) THEN
!                WRITE(*,*) JBEG, JEND, NMOL, BOXL
        do JNEB = JBEG, JEND
!               DO JNEB = JBEG, JEND

                  J = LIST(JNEB)

                  DR(:)  = R(I,:) - R(J,:) - NINT((R(I,:) - R(J,:))/BOXL)*BOXL
                  DCMCM  = DOT_PRODUCT(DR,DR)
                IF(PYT) THEN
                        molj_index = J
                        molj => molecules(molj_index)
                            ! loop over sites in outer molecule
                  do sitei_index = 1, size(moli%sites)
                              sitei => moli%sites(sitei_index)

                      ! loop over sites in inner molecule
                      do sitej_index = 1, size(molj%sites)
                          sitej => molj%sites(sitej_index)

                          energy_contrib = 0.d0
                          grad_contrib(:) = 0.d0

                          ! compute the energy and gradient for this pair of sites
                          call pairwise_py(sitei, sitej, grad_contrib, energy_contrib, .true., t_rep, t_att)
                          ! add the energy to total 
                          PE = PE + energy_contrib
                          PEP(I) = PEP(I) + energy_contrib
                          PEP(J) = PEP(J) + energy_contrib

                           DO K = 1, 3 ! accumulate forces and torques
                                F(I,k) = F(I,k) - grad_contrib(k)
                                F(J,k) = F(J,k) - grad_contrib(k+3)
!                                T(I,k) = T(I,k) - grad_contrib(6+k)
!                                T(J,k) = T(J,k) - grad_contrib(9+k)
                                T(I,k) = T(I,k) - t_rep(k) - t_att(k) 
                                T(J,k) = T(J,k) - t_rep(3+k) - t_att(3+k)
                           END DO
                      end do
                  end do


                ELSE ! LWOTP model

                  CALL MOLMOL_FORQUE(I, J, RSITE, FSITE, DR, DCMCM, PE, VRL)
                END IF
               ENDDO

            ENDIF

         ENDDO

      ENDIF

!     MULTIPLY RESULTS BY ENERGY FACTORS 
      IF(.NOT.PYT) THEN
              FSITE = FSITE * EPS4 * 2.D0
              PE    = PE * EPS4
              VRL   = VRL * EPS4 * 2.D0
      END IF
!     CONVERT SITE FORCES TO MOLECULAR FORCE AND TORQUE

      DO I = 1, NMOL
         DO J = 1, NSITE

            INDX = (I - 1) * NSITE + J

          IF(PYT) THEN
!              CALL AA2QTRN(T(I,:),
!            T(I,1) = T(I,1) - grad_contrib(8) + grad_contrib(9)
!            T(I,2) = T(I,2) - grad_contrib(9) + grad_contrib(7)
!            T(I,3) = T(I,3) - grad_contrib(7) + grad_contrib(8) 
          ELSE
            DO K = 1, 3
               F(I,K) = F(I,K) + FSITE(INDX,K)
            ENDDO
            T(I,1) = T(I,1) + RSITE(INDX,2)*FSITE(INDX,3) - RSITE(INDX,3)*FSITE(INDX,2)
            T(I,2) = T(I,2) + RSITE(INDX,3)*FSITE(INDX,1) - RSITE(INDX,1)*FSITE(INDX,3)
            T(I,3) = T(I,3) + RSITE(INDX,1)*FSITE(INDX,2) - RSITE(INDX,2)*FSITE(INDX,1) 
          END IF
         ENDDO
      
      ENDDO

      END SUBROUTINE FORQUE
 
!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RBDEFMOL()

!     DEFINE THE MOLECULAR GEOMETRY WITH THE CENTRE-OF-MASS AT THE 
!     ORIGIN & CALCULATE THE TOTAL MASS & THE ORINCIPAL OF MOMENTS OF 
!     INERTIA

      USE MDCOMMONS, ONLY : MST, JMI, MASS, MST, RCUT2, RCUT2SQ, RCUT, NSITE

      IMPLICIT NONE

      INTEGER :: I
      DOUBLE PRECISION :: MS(NSITE), PI, DIST(NSITE), DELRC

      PI       = 4.D0 * DATAN(1.D0)

      MST(1,1) = 0.D0
      MST(1,2) = - 2.D0 * DSIN(7.D0*PI/24.D0) / 3.D0
      MST(1,3) = 0.D0

      MST(2,1) = DCOS(7.D0*PI/24.D0)
      MST(2,2) = DSIN(7.D0*PI/24.D0) / 3.D0
      MST(2,3) = 0.D0

      MST(3,1) = - DCOS(7.D0*PI/24.D0)
      MST(3,2) = DSIN(7.D0*PI/24.D0) / 3.D0
      MST(3,3) = 0.D0

      DO I = 1, NSITE

         DIST(I) = DSQRT(DOT_PRODUCT(MST(I,:),MST(I,:)))

      ENDDO

      DELRC    = 2.D0 * MAXVAL(DIST) 
      RCUT2    = RCUT + DELRC
      RCUT2SQ  = RCUT2 * RCUT2
 
!      MS(1)    = 1.D0/3.D0
!      MS(2)    = 1.D0/3.D0
!      MS(3)    = 1.D0/3.D0

!      MASS     = MS(1) + MS(2) + MS(3)

!      JMI(1)  = MS(1)*MST(1,2)**2+MS(2)*MST(2,2)**2+MS(3)*MST(3,2)**2
!      JMI(2)  = MS(1)*MST(1,1)**2+MS(2)*MST(2,1)**2+MS(3)*MST(3,1)**2
!      JMI(3)  = MS(1)*MST(1,1)**2+MS(2)*MST(2,1)**2+MS(3)*MST(3,1)**2   &
!                 + MS(1)*MST(1,2)**2+MS(2)*MST(2,2)**2+MS(3)*MST(3,2)**2
      MASS = 1.0D0
      JMI(:) = 1.0D0
      MST(:,:) = 1.0D0
!        WRITE(*,*) 'JMI=', JMI(:), MASS, MST, NSITE
!STOP
      END SUBROUTINE RBDEFMOL   
       
!     ----------------------------------------------------------------------------------------------

      SUBROUTINE EFM (EFMT, TIME)

      ! calculate EFMT, the energy fluctuation metric: see 
      ! de Sauza and Wales J. Chem. Phys. 123 134504 2005 
      ! http://dx.doi.org/10.1063/1.2035080
      ! calculate energy averaged over time and molecules.
      ! calculate variance over molecules of time averaged energy

      USE MDCOMMONS, ONLY : SUMEPT, NMOL, DELT, PEP, KEP

      IMPLICIT NONE

      INTEGER          :: I
      DOUBLE PRECISION :: AVET, EFMT, EP, TIME, DBLE

      TIME = TIME + DELT
      AVET = 0.D0 ! energy averaged over time and molecules <<E>_T>_N
      EFMT = 0.D0 ! variance of average energy over molecules? <(<E>_T-AVET)^2>_N

      DO I = 1, NMOL

         EP        = 0.5D0 * KEP(I) + PEP(I)
         SUMEPT(I) = SUMEPT(I) + EP * DELT
         AVET      = AVET + SUMEPT(I) 

      ENDDO

      AVET = AVET / (TIME * DBLE(NMOL))

      DO I = 1, NMOL

         EFMT = EFMT + (SUMEPT(I) / TIME - AVET) ** 2.D0

      ENDDO

      EFMT = EFMT / DBLE(NMOL)

      END SUBROUTINE EFM

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RIGIDMD_SC ()

! initialize the CoM positions on a simple cubic lattice

      USE MDCOMMONS, only : R, NMOL
      USE COMMONS, only : boxlx

      IMPLICIT NONE

      INTEGER          :: I, NUC, IX, IY, IZ, M!, K
      DOUBLE PRECISION :: BOXLH, UCL, UCLH, C(3), DBLE

!     NUMBER OF UNIT CELLS
      BOXLH = 0.5D0 * BOXLX

      I = 1

      DO WHILE (I ** 3 < NMOL)
         I = I + 1
      ENDDO

      NUC   = I

!     UNIT CELL LENGTH

      UCL     = BOXLX / DBLE(NUC)
      UCLH    = 0.5D0 * UCL

!     CONSTRUCT THE LATTICE FROM THE UNIT CELL

      M = 0

      DO IZ = 1, NUC

         C(3) = IZ * UCL - UCLH - BOXLH

         DO IY = 1, NUC

            C(2) = IY * UCL - UCLH - BOXLH

            DO IX = 1, NUC

               C(1) = IX * UCL - UCLH - BOXLH

               M = M + 1
               R(M,:) = C(:)

            ENDDO

         ENDDO

      ENDDO

      END SUBROUTINE RIGIDMD_SC

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE EQCON()

! load positions, orientations, linear velocities, and angular
! velocities from files

      USE MDCOMMONS , only : R, QTRN, V, W, NMOL

      IMPLICIT NONE

      INTEGER :: I

      OPEN (UNIT=13, FILE='initpos1.dat', STATUS='UNKNOWN')
      OPEN (UNIT=14, FILE='initortn1.dat', STATUS='UNKNOWN')
      OPEN (UNIT=15, FILE='initlinvel1.dat', STATUS='UNKNOWN')
      OPEN (UNIT=16, FILE='initangvel1.dat', STATUS='UNKNOWN')

      DO I = 1, NMOL  
        
         READ(13, *) R(I,1), R(I,2), R(I,3)
         !READ(14, *) QTRN(I,3), QTRN(I,2), QTRN(I,4), QTRN(I,1)
         READ(14, *) QTRN(I,1), QTRN(I,2), QTRN(I,3), QTRN(I,4)
         READ(15, *) V(I,1), V(I,2), V(I,3)
         READ(16, *) W(I,1), W(I,2), W(I,3)
         !QTRN(I,3) = -QTRN(I,3)

      ENDDO

      
      CLOSE (13)
      CLOSE (14)
      CLOSE (15)
      CLOSE (16)

      END SUBROUTINE EQCON

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE LOAD_POS_ORTN()

! load positions, orientations, linear velocities, and angular
! velocities from files

      USE MDCOMMONS , only : R, QTRN, NMOL

      IMPLICIT NONE

      INTEGER :: I

      OPEN (UNIT=13, FILE='initpos1.dat', STATUS='UNKNOWN')
      OPEN (UNIT=14, FILE='initortn1.dat', STATUS='UNKNOWN')
      !OPEN (UNIT=15, FILE='initlinvel1.dat', STATUS='UNKNOWN')
      !OPEN (UNIT=16, FILE='initangvel1.dat', STATUS='UNKNOWN')
        WRITE(*,*) 'NMOL=', NMOL
      DO I = 1, NMOL  
        
         READ(13, *) R(I,1), R(I,2), R(I,3)
         !READ(14, *) QTRN(I,3), QTRN(I,2), QTRN(I,4), QTRN(I,1)
         READ(14, *) QTRN(I,1), QTRN(I,2), QTRN(I,3), QTRN(I,4)
         !READ(15, *) V(I,1), V(I,2), V(I,3)
         !READ(16, *) W(I,1), W(I,2), W(I,3)
         !QTRN(I,3) = -QTRN(I,3)
         WRITE(*,*) "Read molecule ", I
         WRITE(*,*) R(I,1), R(I,2), R(I,3)
      ENDDO

      
      CLOSE (13)
      CLOSE (14)
      !CLOSE (15)
      !CLOSE (16)

      END SUBROUTINE LOAD_POS_ORTN

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE LOAD_LIN_ANG_VEL()

! load positions, orientations, linear velocities, and angular
! velocities from files

      USE MDCOMMONS , only : V, W, NMOL

      IMPLICIT NONE

      INTEGER :: I

      !OPEN (UNIT=13, FILE='initpos1.dat', STATUS='UNKNOWN')
      !OPEN (UNIT=14, FILE='initortn1.dat', STATUS='UNKNOWN')
      OPEN (UNIT=15, FILE='initlinvel1.dat', STATUS='UNKNOWN')
      OPEN (UNIT=16, FILE='initangvel1.dat', STATUS='UNKNOWN')

      DO I = 1, NMOL  
        
         !READ(13, *) R(I,1), R(I,2), R(I,3)
         !READ(14, *) QTRN(I,3), QTRN(I,2), QTRN(I,4), QTRN(I,1)
         !READ(14, *) QTRN(I,1), QTRN(I,2), QTRN(I,3), QTRN(I,4)
         READ(15, *) V(I,1), V(I,2), V(I,3)
         READ(16, *) W(I,1), W(I,2), W(I,3)
         !QTRN(I,3) = -QTRN(I,3)

      ENDDO

      
      CLOSE (13)
      CLOSE (14)
      !CLOSE (15)
      !CLOSE (16)

      END SUBROUTINE LOAD_LIN_ANG_VEL

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE INTVEL(TMPINT)

! initialize linear velocities and angular velocity

      USE MDCOMMONS, only : MASS, JMI, V, NMOL, W

      IMPLICIT NONE

      INTEGER          :: I
      DOUBLE PRECISION :: TMPINT, VMAG, WMAG, SUMV(3), AVV(3), E(3), DBLE
      
      VMAG    = DSQRT(3.D0 * TMPINT / MASS)
      WMAG = DSQRT(3.D0 * TMPINT / (JMI(1) + JMI(2) + JMI(3)))
      SUMV(:) = 0.D0
        write(*,*) 'JMI(:)=', JMI(:), WMAG, VMAG, MASS
      DO I = 1, NMOL

         CALL RNDVCT(E)

         V(I,:) = VMAG * E(:)
         SUMV(:) = SUMV(:) + V(I,:)
         write(*,*) 'V(I)=', V(I,:)
          
      ENDDO 

      AVV(:) = SUMV(:) / DBLE(NMOL)
     
      DO I = 1, NMOL

         V(I,:) = V(I,:) - AVV(:)

         CALL RNDVCT(E)
         W(I,:) = WMAG * E(:)

      ENDDO
        write(*,*) 'after intvel', V(:,:)
      END SUBROUTINE INTVEL 

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE INTORN()

! initialize orientations

      USE MDCOMMONS, only : NMOL, QTRN

      IMPLICIT NONE

      INTEGER          :: I
      DOUBLE PRECISION :: NORM

!     SET QUATERNIONS WITH RANDOM ORIENTATIONS
      
      DO I =1, NMOL

         QTRN(I,1) = 0.10
         QTRN(I,2) = -0.415D0
         QTRN(I,3) = -0.246D0
         QTRN(I,4) = -0.976D0

         NORM   = DSQRT(DOT_PRODUCT(QTRN(I,:),QTRN(I,:)))

         QTRN(I,:) = QTRN(I,:) / NORM
 
      ENDDO

      END SUBROUTINE INTORN

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RNDVCT(E)

!     CHOOSE A RANDOM UNIT VECTOR IN SPACE

      IMPLICIT NONE

      DOUBLE PRECISION :: E(3), DUMMY, XI, XISQ, XI1, XI2

      XISQ = 2.D0

      DO WHILE (XISQ > 1.D0)

         CALL RANDOM_NUMBER (DUMMY)

         XI1  = DUMMY * 2.D0 - 1.D0

         CALL RANDOM_NUMBER (DUMMY)

         XI2  = DUMMY * 2.D0 - 1.D0
         XISQ = XI1 * XI1 + XI2 * XI2

      ENDDO

      XI    = DSQRT(1.D0 - XISQ)

      E(1) = 2.D0 * XI1 * XI
      E(2) = 2.D0 * XI2 * XI
      E(3) = 1.D0 - 2.D0 * XISQ

      END SUBROUTINE RNDVCT

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE CHECK ()

      USE MDCOMMONS, ONLY : NMOL, R, R0, RCUT2, RLIST, UPDATE

!     DECIDES WHETHER THE LIST NEEDS TO BE RECONSTRUCTED FROM THE COORDINATES AT LAST UPDATE

!     LOGICAL  UPDATE    IF TRUE THE LIST IS UPDATED

!     CHECK IS CALLED TO SET UPDATE BEFORE EVERY CALL TO FORQUE

      IMPLICIT NONE 

      INTEGER          :: I, K
      DOUBLE PRECISION :: DISPMX

!     CALCULATE MAXIMUM DISPLACEMENT SINCE LAST UPDATE

      DISPMX = 0.D0

      DO I = 1, NMOL

         DO K = 1, 3

            DISPMX = MAX(ABS(R(I,K) - R0(I,K)), DISPMX)

         ENDDO
!        WRITE(*,*) 'DISPMX=', DISPMX
      ENDDO

!     A CONSERVATIVE TEST OF THE LIST SKIN CROSSING

      DISPMX = 2.D0 * DSQRT(3.D0 * DISPMX ** 2)

      UPDATE = (DISPMX > (RLIST - RCUT2))
!      WRITE(*,*) 'UPDATE', UPDATE
      END SUBROUTINE CHECK 

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE SVCNFG()

! save config

      USE MDCOMMONS, only : R, R0, NMOL
       
!     SAVE IS CALLED WHENEVER THE NEW VERLET LIST IS CONSTRUCTED

      IMPLICIT NONE

      INTEGER :: I, K

      DO I = 1, NMOL

         DO K = 1, 3

           R0(I,K) = R(I,K)

         ENDDO

      ENDDO

      END SUBROUTINE SVCNFG

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE RUN_INFO()
       
!     print some info about the run.  time elapsed, etc.

      !USE MDCOMMONS

      IMPLICIT NONE

      real :: timearray(2), timesecs

      OPEN (UNIT = 40, FILE = 'runinfo.dat', STATUS = 'UNKNOWN')

!      call etime( timearray, timesecs )

      write(40,*) "time_elapsed ", timesecs, timearray

      close(40)

      END SUBROUTINE RUN_INFO

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE GET_RSITE( RSITE )
       
!     GET the positions of the sites with respect to the CoM
!     i.e. RSITE_fixed_space(I,:) = RSITE(I,:) + R(I,:)
!     use MST and the orientations QTRN

      USE MDCOMMONS, only : Nmol, QTRN, NSITE, MST, NTST

      IMPLICIT NONE

      DOUBLE PRECISION :: RSITE(NTST,3)
      DOUBLE PRECISION :: Q(4), RMX(3,3)
      integer :: INDX, I, J

      DO I = 1, NMOL

!     CONSTRUCT THE ROTATION MATRIX RMX FROM THE QUATERNION

         Q    = QTRN(I,:)

         CALL CONSTRUCT_RMX( RMX, Q)

         DO J = 1, NSITE

            INDX = (I - 1) * NSITE + J

            !the site coordinate in fixed space frame
            RSITE(INDX,:) = MATMUL(RMX,MST(J,:))
    
         ENDDO

      ENDDO

      END SUBROUTINE GET_RSITE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE MOLMOL_FORQUE(I, J, RSITE, FSITE, DR, DCMCM, PE, VRL)
       
!     GET the positions of the sites with respect to the CoM
!     i.e. RSITE_fixed_space(I,:) = RSITE(I,:) + R(I,:)
!     use MST and the orientations QTRN

      USE MDCOMMONS, only : RCUT2SQ, NSITE, RCUTSQ, PEP, NTST, CNSTA, CNSTB

      IMPLICIT NONE

      INTEGER          :: I, J, J1, J2, INDX1, INDX2, K
      DOUBLE PRECISION :: RSITE(NTST,3), FSITE(NTST,3)
      !DOUBLE PRECISION :: Q(4), Q2Q3, Q1Q4, Q2Q4, Q1Q3, Q3Q4, Q1Q2, RMX(3,3)
      DOUBLE PRECISION :: DR(3), DSS(3), DCMCM, DSS2, R2, R6, R12!, SHIFT(3)
      DOUBLE PRECISION :: PE, VJ1J2, FJ1J2, VRL

      IF (DCMCM < RCUT2SQ) THEN

        INDX1 = (I - 1) * NSITE
        INDX2 = (J - 1) * NSITE

        DO J1 = 1, NSITE

          DO J2 = 1, NSITE

            ! DSS = distance between the sites
            DSS(:) = DR(:) + RSITE(INDX1+J1,:) - RSITE(INDX2+J2,:)
            DSS2   = DOT_PRODUCT(DSS,DSS)

            if (DSS2 < RCUTSQ ) THEN
              R2     = 1.D0 / DSS2
              R6     = R2 * R2 * R2
              R12    = R6 ** 2
              VJ1J2  = R12 - R6 + CNSTA * DSS2 + CNSTB
              PE     = PE + VJ1J2
              PEP(I) = PEP(I) + 2.D0 * VJ1J2
              PEP(J) = PEP(J) + 2.D0 * VJ1J2
              FJ1J2  = (6.D0 * R12 - 3.D0 * R6)/DSS2 - CNSTA !half the force
              VRL    = VRL + FJ1J2 * DSS2 !used to calculate the pressure

              DO K = 1, 3

                FSITE(INDX1+J1,K) = FSITE(INDX1+J1,K) + FJ1J2*DSS(K)
                FSITE(INDX2+J2,K) = FSITE(INDX2+J2,K) - FJ1J2*DSS(K)

              ENDDO

            ENDIF

          ENDDO

        ENDDO

      ENDIF

      END SUBROUTINE MOLMOL_FORQUE

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE CONSTRUCT_MX( MX, Q )

!     CONSTRUCT THE MATRIX MX FROM THE QUATERNION Q
! MX = the matrix representation of Q
!  Q1 -Q2 -Q3 -Q4 
!  Q2  Q1 -Q4  Q3
!  Q3  Q4  Q1 -Q2
!  Q4 -Q3  Q2  Q1
!        OR, with indices from 0 to 3
!  Q0 -Q1 -Q2 -Q3 
!  Q1  Q0 -Q3  Q2
!  Q2  Q3  Q0 -Q1
!  Q3 -Q2  Q1  Q0


      !USE MDCOMMONS

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT) :: MX(4,4)
      DOUBLE PRECISION, INTENT(IN) :: Q(4)
      INTEGER :: K

      DO K = 1, 4
        MX(K,K) = Q(1)
      ENDDO

      MX(1,2:4) = -Q(2:4)
      MX(2:4,1) = Q(2:4)
      MX(2,3)   = -Q(4)
      MX(2,4)   = Q(3)
      MX(3,2)   = Q(4)
      MX(3,4)   = -Q(2)
      MX(4,2)   = -Q(3)
      MX(4,3)   = Q(2)

      END SUBROUTINE CONSTRUCT_MX

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE CONSTRUCT_RMX( RMX, Q )

      !construct the rotation matrix which rotates a 3-vector from the body
      !fixed frame to the space fixed frame
      !uses the fact that the quaternion is normalized
      !1=a^2+b^2+c^2+d^2

      !USE MDCOMMONS

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT) :: RMX(3,3)
      DOUBLE PRECISION, INTENT(IN) :: Q(4)
      DOUBLE PRECISION :: Q2Q3, Q1Q4, Q2Q4, Q1Q3, Q3Q4, Q1Q2 

      Q2Q3 = Q(2)*Q(3)
      Q1Q4 = Q(1)*Q(4)
      Q2Q4 = Q(2)*Q(4)
      Q1Q3 = Q(1)*Q(3)
      Q3Q4 = Q(3)*Q(4)
      Q1Q2 = Q(1)*Q(2)

      RMX(1,1) = 0.5D0 - Q(3)*Q(3) - Q(4)*Q(4)
      RMX(2,2) = 0.5D0 - Q(2)*Q(2) - Q(4)*Q(4)
      RMX(3,3) = 0.5D0 - Q(2)*Q(2) - Q(3)*Q(3)
      RMX(1,2) = Q2Q3 - Q1Q4
      RMX(2,1) = Q2Q3 + Q1Q4
      RMX(1,3) = Q2Q4 + Q1Q3
      RMX(3,1) = Q2Q4 - Q1Q3
      RMX(2,3) = Q3Q4 - Q1Q2
      RMX(3,2) = Q3Q4 + Q1Q2

      RMX      = 2.D0*RMX

      END SUBROUTINE CONSTRUCT_RMX

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE CONSTRUCT_RMXT( RMXT, Q )

      !construct the rotation matrix which rotates a 3-vector from the space
      !fixed frame to the body frame
      !uses the fact that the quaternion is normalized
      !1=a^2+b^2+c^2+d^2

      !USE MDCOMMONS

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT) :: RMXT(3,3)
      DOUBLE PRECISION, INTENT(IN) :: Q(4)
      DOUBLE PRECISION :: Q2Q3, Q1Q4, Q2Q4, Q1Q3, Q3Q4, Q1Q2 

      Q2Q3 = Q(2)*Q(3)
      Q1Q4 = Q(1)*Q(4)
      Q2Q4 = Q(2)*Q(4)
      Q1Q3 = Q(1)*Q(3)
      Q3Q4 = Q(3)*Q(4)
      Q1Q2 = Q(1)*Q(2)

      RMXT(1,1) = 0.5D0 - Q(3)*Q(3) - Q(4)*Q(4)
      RMXT(2,2) = 0.5D0 - Q(2)*Q(2) - Q(4)*Q(4)
      RMXT(3,3) = 0.5D0 - Q(2)*Q(2) - Q(3)*Q(3)
      RMXT(1,2) = Q2Q3 + Q1Q4
      RMXT(2,1) = Q2Q3 - Q1Q4
      RMXT(1,3) = Q2Q4 - Q1Q3
      RMXT(3,1) = Q2Q4 + Q1Q3
      RMXT(2,3) = Q3Q4 + Q1Q2
      RMXT(3,2) = Q3Q4 - Q1Q2

      RMXT      = 2.D0*RMXT

      END SUBROUTINE CONSTRUCT_RMXT

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE QTRN2AA( Q, V )

      !convert from quaternion to angle axis representation of a
      !rotation

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT) :: V(3)
      DOUBLE PRECISION, INTENT(IN) :: Q(4)
      DOUBLE PRECISION :: THETA, D

      THETA = 2.D0*ACOS(Q(1))
      D     = DSQRT(1.D0-Q(1)*Q(1))
      V(:)  = THETA * Q(2:4)/D

      END SUBROUTINE QTRN2AA

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE AA2QTRN( Q, V )

      !convert from angle axis to quaternion representation of a rotation

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: V(3)
      DOUBLE PRECISION, INTENT(OUT) :: Q(4)
      DOUBLE PRECISION :: THETA

      THETA = SQRT( V(1)* V(1)+ V(2)* V(2)+ V(3)* V(3) )
      Q(1)    = COS( THETA / 2.D0 )
      Q(2:4) = SQRT( 1-Q(1)*Q(1) ) * V(:) / THETA

      END SUBROUTINE AA2QTRN


!     ----------------------------------------------------------------------------------------------

      SUBROUTINE PfromW( P, W, JMI, Q, S )

      !convert from angle axis to quaternion representation of a rotation

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: W(3), JMI(3), Q(4), S
      DOUBLE PRECISION, INTENT(OUT) :: P(4)
      DOUBLE PRECISION :: W4(4), MX(4,4)

      ! intitialize P, the angular momentum.  P is virtual, W is real
      CALL CONSTRUCT_MX( MX, Q )
      W4(1)  = 0.D0
      W4(2)  = 2.D0*W(1)*JMI(1)
      W4(3)  = 2.D0*W(2)*JMI(2)
      W4(4)  = 2.D0*W(3)*JMI(3)
      P(:) = MATMUL(MX*S,W4)


      END SUBROUTINE PfromW

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE WfromP( P, W, JMI, Q, S )

      !convert from angle axis to quaternion representation of a rotation

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: P(4), JMI(3), Q(4), S
      DOUBLE PRECISION, INTENT(OUT) :: W(3)
      DOUBLE PRECISION :: W4(4), MX(4,4)

      CALL CONSTRUCT_MX( MX, Q )
      W4     = MATMUL(TRANSPOSE(MX),(P(:)/S))
      W(1) = 0.5D0*W4(2)/JMI(1)
      W(2) = 0.5D0*W4(3)/JMI(2)
      W(3) = 0.5D0*W4(4)/JMI(3)

      END SUBROUTINE WfromP

!
! svn info 
! $Date: 2011-10-16 20:53:28 +0100 (Sun, 16 Oct 2011) $
! $Rev: 61 $
!

 
!     ----------------------------------------------------------------------------------------------
 
      SUBROUTINE MOVE (ISTEP, NEQ, PE, VRL, TKE, RKE)

!     propagator: H. Okumura, S. G. Itoh, and Y. Okamoto, J. Chem. Phys. 126, 084103 (2007).
      !note: move operates on whole molecules, except through the subroutine
      !FORQUE

      USE MDCOMMONS
      USE COMMONS, only : PYT

      IMPLICIT NONE

      ! PHIT = This parameter is called ZETA in the above reference.  As far as
      !        I can tell it is not anything physical.  one can define it
      !        component by component like 
      !        PHIT(J) = DOT_PRODUCT( P(I,:) * MX(:,J) ) / (4 * JMI(J) * S)
      !        where MX is the matrix representation of QTRN(I,:) as defined
      !        below
      INTEGER          :: ISTEP, NEQ, I
      !DOUBLE PRECISION :: MX(4,4), W4(4), RMXT(3,3) 
      DOUBLE PRECISION :: RMXT(3,3) 
      DOUBLE PRECISION :: PE, VRL, TKE, RKE, RKET
      DOUBLE PRECISION :: Q(4),T4dummy(4)

      RKE = 0.D0

      ! S = S is the additional degree of freedom in the Nose-Poincare thermostat
      ! MQ = the artificial mass associated with S
      ! PS = the conjugate momentum associated with S



      ! STEP 1: update S, PS
      CALL NVESTEP1( S, PS, FCTMQT )

      ! STEP 2: update PS, V from F, and P from T4
      CALL NVESTEP2(PE)

      ! STEP 3: update QTRN, P, PS with contributions from principle axis 3
      CALL NVESTEP3()
      
      ! STEP 4: update QTRN, P, PS with contributions from principle axis 2
      CALL NVESTEP4()

      ! STEP 5: update QTRN, P, PS with contributions from principle axis 1.  Also update R from V
      CALL NVESTEP5()

      ! STEP 6: update QTRN, P, PS with contributions from principle axis 2
      ! STEP 6 is exactly the same as STEP 4
      CALL NVESTEP4()

      ! STEP 7: update QTRN, P, PS with contributions from principle axis 3
      ! STEP 7 is exactly the same as STEP 3
      CALL NVESTEP3()

      
      !question: why is FORQUE called here?????
      CALL CHECK ()
      CALL FORQUE (PE, VRL)
      ! update T4
      DO I = 1, NMOL
         Q = QTRN(I,:)
         CALL CONSTRUCT_RMXT( RMXT, Q )
         T4(I,1)   = 0.D0
!       IF(PYT) THEN
        
!         CALL AA2QTRN(T4dummy(:),T(I,:))
!         T4(I,2:4) = MATMUL(RMXT,T4dummy(2:4))
!         T4(I,1)   = 0.D0
!       ELSE
         T4(I,2:4) = MATMUL(RMXT,T(I,:))
!       END IF
      ENDDO

      ! STEP 8: update PS, V from F, and P from T4
      ! STEP 8 is exactly the same as STEP 2
      CALL NVESTEP2(PE)

      ! STEP 9: update S, PS
      ! STEP 9 is exactly the same as STEP 1
      !RKET, RKE, PHI, etc.????
      CALL NVESTEP1( S, PS, FCTMQT )



      ! calculate W from P (uses MX)
      ! update KEP, TKE, RKET, RKE, PHI
      TKE    = 0.D0
      DO I = 1, NMOL
         KEP(I) = MASS * (V(I,1)*V(I,1) + V(I,2)*V(I,2) + V(I,3)*V(I,3))/(S*S) !here V is virtual velocity
         
         TKE = TKE + KEP(I) 

!        ADVANCE ANGULAR VELOCITY

         Q    = QTRN(I,:)

         !NOTE: if I incorporate the W update into STEP 8 I can avoid
         !this extra call to CONSTRUCT_MX
         !CALL CONSTRUCT_MX( MX, Q )
         !W4     = MATMUL(TRANSPOSE(MX),(P(I,:)/S))
         !W(I,1) = 0.5D0*W4(2)/JMI(1)
         !W(I,2) = 0.5D0*W4(3)/JMI(2)
         !W(I,3) = 0.5D0*W4(4)/JMI(3)
         CALL WfromP( P(I,:), W(I,:), JMI, Q, S )

         RKET   = JMI(1)*W(I,1)*W(I,1) + JMI(2)*W(I,2)*W(I,2) + JMI(3)*W(I,3)*W(I,3)
         KEP(I) = KEP(I) + RKET
         RKE    = RKE + RKET
!         RKE = 0
         IF (ISTEP > NEQ) THEN

            PHI(I,:) = PHI(I,:) + W(I,:) * DELT

         ENDIF

      ENDDO 

      TKE = 0.5D0 * TKE
      RKE = 0.5D0 * RKE



      END SUBROUTINE MOVE
!     ----------------------------------------------------------------------------------------------
 
      SUBROUTINE NVTMOVE (ISTEP, NEQ, PE, VRL, TKE, RKE)

!     propagator: H. Okumura, S. G. Itoh, and Y. Okamoto, J. Chem. Phys. 126, 084103 (2007).
      !note: move operates on whole molecules, except through the subroutine
      !FORQUE

      USE MDCOMMONS
      USE COMMONS, only : PYT

      IMPLICIT NONE

      ! PHIT = This parameter is called ZETA in the above reference.  As far as
      !        I can tell it is not anything physical.  one can define it
      !        component by component like 
      !        PHIT(J) = DOT_PRODUCT( P(I,:) * MX(:,J) ) / (4 * JMI(J) * S)
      !        where MX is the matrix representation of QTRN(I,:) as defined
      !        below
      INTEGER          :: ISTEP, NEQ, I
      !DOUBLE PRECISION :: MX(4,4), W4(4), RMXT(3,3) 
      DOUBLE PRECISION :: RMXT(3,3) 
      DOUBLE PRECISION :: PE, VRL, TKE, RKE, RKET
      DOUBLE PRECISION :: Q(4),T4dummy(4)

      RKE = 0.D0

      ! S = S is the additional degree of freedom in the Nose-Poincare thermostat
      ! MQ = the artificial mass associated with S
      ! PS = the conjugate momentum associated with S



      ! STEP 1: update S, PS
      CALL NVTSTEP1( S, PS, FCTMQT )

      ! STEP 2: update PS, V from F, and P from T4
      CALL NVTSTEP2(PE)

      ! STEP 3: update QTRN, P, PS with contributions from principle axis 3
      CALL NVTSTEP3()
      
      ! STEP 4: update QTRN, P, PS with contributions from principle axis 2
      CALL NVTSTEP4()

      ! STEP 5: update QTRN, P, PS with contributions from principle axis 1.  Also update R from V
      CALL NVTSTEP5()

      ! STEP 6: update QTRN, P, PS with contributions from principle axis 2
      ! STEP 6 is exactly the same as STEP 4
      CALL NVTSTEP4()

      ! STEP 7: update QTRN, P, PS with contributions from principle axis 3
      ! STEP 7 is exactly the same as STEP 3
      CALL NVTSTEP3()

      
      !question: why is FORQUE called here?????
      CALL CHECK ()
      CALL FORQUE (PE, VRL)
      ! update T4
      DO I = 1, NMOL
         Q = QTRN(I,:)
         CALL CONSTRUCT_RMXT( RMXT, Q )
         T4(I,1)   = 0.D0
!       IF(PYT) THEN
        
!         CALL AA2QTRN(T4dummy(:),T(I,:))
!         T4(I,2:4) = MATMUL(RMXT,T4dummy(2:4))
!         T4(I,1)   = 0.D0
!       ELSE
         T4(I,2:4) = MATMUL(RMXT,T(I,:))
!       END IF
      ENDDO

      ! STEP 8: update PS, V from F, and P from T4
      ! STEP 8 is exactly the same as STEP 2
      CALL NVTSTEP2(PE)

      ! STEP 9: update S, PS
      ! STEP 9 is exactly the same as STEP 1
      !RKET, RKE, PHI, etc.????
      CALL NVTSTEP1( S, PS, FCTMQT )



      ! calculate W from P (uses MX)
      ! update KEP, TKE, RKET, RKE, PHI
      TKE    = 0.D0
      DO I = 1, NMOL
         KEP(I) = MASS * (V(I,1)*V(I,1) + V(I,2)*V(I,2) + V(I,3)*V(I,3))/(S*S) !here V is virtual velocity
         
         TKE = TKE + KEP(I) 

!        ADVANCE ANGULAR VELOCITY

         Q    = QTRN(I,:)

         !NOTE: if I incorporate the W update into STEP 8 I can avoid
         !this extra call to CONSTRUCT_MX
         !CALL CONSTRUCT_MX( MX, Q )
         !W4     = MATMUL(TRANSPOSE(MX),(P(I,:)/S))
         !W(I,1) = 0.5D0*W4(2)/JMI(1)
         !W(I,2) = 0.5D0*W4(3)/JMI(2)
         !W(I,3) = 0.5D0*W4(4)/JMI(3)
         CALL WfromP( P(I,:), W(I,:), JMI, Q, S )

         RKET   = JMI(1)*W(I,1)*W(I,1) + JMI(2)*W(I,2)*W(I,2) + JMI(3)*W(I,3)*W(I,3)
         KEP(I) = KEP(I) + RKET
         RKE    = RKE + RKET
!         RKE = 0
         IF (ISTEP > NEQ) THEN

            PHI(I,:) = PHI(I,:) + W(I,:) * DELT

         ENDIF

      ENDDO 

      TKE = 0.5D0 * TKE
      RKE = 0.5D0 * RKE



      END SUBROUTINE NVTMOVE


      SUBROUTINE GET_FRAME_R_ORNT(IOstatus, AABOOL)

! load positions, orientations, linear velocities, and angular
! velocities from files

      USE MDCOMMONS, ONLY : R, QTRN, NMOL 

      IMPLICIT NONE

      INTEGER :: I,IOstatus
      LOGICAL, INTENT(IN) :: AABOOL
      DOUBLE PRECISION :: V(3)

      if ( AABOOL ) then

      DO I = 1, NMOL  
        
         READ(13, *, IOSTAT=IOstatus) R(I,1), R(I,2), R(I,3)
         IF (IOstatus .NE. 0) return
      ENDDO
      DO I = 1, NMOL  
         READ(13, *, IOSTAT=IOstatus) V(1), V(2), V(3) 
         IF (IOstatus .NE. 0) return
         CALL AA2QTRN( QTRN(I,:), V )
      ENDDO

      ELSE
      DO I = 1, NMOL  
        
         READ(13, *, IOSTAT=IOstatus) R(I,1), R(I,2), R(I,3)
         IF (IOstatus .NE. 0) return
         READ(14, *, IOSTAT=IOstatus) QTRN(I,1), QTRN(I,2), QTRN(I,3), QTRN(I,4)
         IF (IOstatus .NE. 0) return
         !READ(15, *) V(I,1), V(I,2), V(I,3)
         !READ(16, *) W(I,1), W(I,2), W(I,3)
         !QTRN(I,3) = -QTRN(I,3)

      ENDDO

    ENDIF

      

      END SUBROUTINE GET_FRAME_R_ORNT

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE SITECOORDS(RSITE)

! load positions, orientations, linear velocities, and angular
! velocities from files

      USE MDCOMMONS

      IMPLICIT NONE

      DOUBLE PRECISION :: Q(4), RMX(3,3)
      INTEGER          :: I, J, INDX!, INDX1, INDX2, K, JBEG, JEND, JNEB, NLIST 
      DOUBLE PRECISION, INTENT(OUT) :: RSITE(NTST,3)

!     OBTAIN THE SITE POSITIONS IN THE SPACE-FIXED FRAME FROM THE
!     CENTRE-OF-MASS POSITION & THE SITE POSITIONS IN THE BODY-FIXED
!     FRAME USING THE ROTATION MATRIX

      DO I = 1, NMOL

!     CONSTRUCT THE ROTATION MATRIX RMX FROM THE QUATERNION

         Q    = QTRN(I,:)
         CALL CONSTRUCT_RMX(RMX,Q)

         DO J = 1, NSITE

            INDX = (I - 1) * NSITE + J

            !the site coordinate in CoM position
            RSITE(INDX,:) = MATMUL(RMX,MST(J,:)) + R(I,:)
    
         ENDDO

      ENDDO

      END SUBROUTINE SITECOORDS

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE POUT(RSITE)

      USE MDCOMMONS

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(OUT) :: RSITE(NTST,3)
      INTEGER:: I
      CHARACTER(LEN=3):: NAM(3)
      NAM(1) = "N  "
      NAM(2) = "C  "
      NAM(3) = "N  "
     

      WRITE (1, *) NTST
      WRITE (1, *) 
      DO I = 1,NTST
        !WRITE (1, *) "C ", RSITE(I,1), RSITE(I,2), RSITE(I,3)
        !WRITE (1, *, ADVANCE='NO' ) NAM(MODULO(I,3)+1), RSITE(I,1), RSITE(I,2), RSITE(I,3)
        WRITE (1, 900 ) NAM(MODULO(I,3)+1), RSITE(I,1), RSITE(I,2), RSITE(I,3)
        900 FORMAT (1x,A3,1x,ES20.12,1x,ES20.12,1x,ES20.12)
      ENDDO

      END SUBROUTINE POUT

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE APPLY_PERIODIC( PERIODIC_OFFSET, FRAMENUM )

      ! apply periodic boundary conditions, but do it the same way each
      ! frame to make the visual representation more appealing.

      USE MDCOMMONS
      USE commons, only : boxlx

      IMPLICIT NONE

      INTEGER :: I,J, FRAMENUM
      DOUBLE PRECISION :: PERIODIC_OFFSET(N,3)
      
      IF ( FRAMENUM .EQ. 1 ) THEN
        DO I=1,NMOL

          !apply periodic boundary conditions
          !the coordinate might be more than one box length outside of
          !the box
          DO J=1,3
            DO WHILE ( R(I,J) .GE. 0.5D0*BOXLX ) 
              R(I,J)=R(I,J)-BOXLX
              PERIODIC_OFFSET(I,J) = PERIODIC_OFFSET(I,J) - BOXLX
            END DO
            DO WHILE ( R(I,J) .LT. -0.5D0*BOXLX ) 
              R(I,J)=R(I,J)+BOXLX
              PERIODIC_OFFSET(I,J) = PERIODIC_OFFSET(I,J) + BOXLX
            END DO
          END DO
        END DO
      ELSE
        DO I=1,NMOL
          DO J=1,3
            R(I,J)=R(I,J)+ PERIODIC_OFFSET(I,J);
          END DO
        END DO
      ENDIF


      END SUBROUTINE APPLY_PERIODIC
!
! svn info 
! $Date: 2011-10-12 12:06:51 +0100 (Wed, 12 Oct 2011) $
! $Rev: 54 $
!

! the lewis wahnstrom model has 
      ! sigma = 0.483 nm
      ! epsilon = 5.276 kJ/mol
 
!     ----------------------------------------------------------------------------------------------

!      PROGRAM LWOTP_MD_NVE
      SUBROUTINE RIGIDMD_NVE

      USE MDCOMMONS
      USE COMMONS, ONLY : PYT

      IMPLICIT NONE

! IDUMP1: dump files.  always: 
!                     enrg1.dat: KEPP, PEPP, EPP: kinetic, potential, energy per molecule
!                     tmpprs1.dat: temp, pressure?
!                     avenrg1.dat: average KE, PE, E
!                     avtmpprs1.dat: average temp, pressure?
!                     avtmptr1.dat: average tranlational, rotational temperature
! IDUMP2: dump files.  ISTEP > NEQ
!                     pos1.dat:  xyz positions of CoM of molecules
!                     ortn1.dat: quaternions giving the orientation of molecules
!                     rot1.dat:  phi1, phi2, phi3 ????
!                     efm1.dat:  energy fluctuation metric: time from NEQ, and variance over molecules of time averaged energy 
! IDUMP3: CALL SCLVEL.  ISTEP <= NBATH: scale linear velocities ???
! IDUMP4: dump files.  ISTEP > NEQ: 
!                     blockav1.dat: average PE, temp, pressure
!                     lvel1.dat:    linear velocities
!                     avel1.dat:    angular velocities
! TKE: translational kinetic energy
! RKE: rotational kinetic energy
! MX:  is a representation of the quaternion as a matrix in such a way
!      that quaternion addition and multiplication correspond to matrix
!      addition and matrix multiplication
!      In this way the conjugate of the quaternion corresponds to the
!      transpose of the matrix.  The fourth power of the norm of a
!      quaternion is the determinant of the corresponding matrix.
!      Complex numbers are block diagonal matrices with two 2x2 blocks.
! EFMT: Energy Fluctuation Metric: the variance over molecules of the
!       time averaged energy
! VRL: accumulates the forces dotted into the positions Fij.dot(Ri - Rj) which
!      is used to calculate the pressure
      INTEGER          :: NSTEP, NBATH, NEQ, ISTEP, ICNT, IDUMP1, IDUMP2, IDUMP3, IDUMP4, I!, K!, STATUS
      DOUBLE PRECISION :: RCUT6, RCUT12, TMP, TMPTR, TMPRT, RHO, VLM
      DOUBLE PRECISION :: Q(4), RMXT(3,3), MX(4,4), W4(4), T4dummy(4)
      DOUBLE PRECISION :: PRS, SUMPRS, AVPRS, PE, KE, TKE, RKE, PEPP, KEPP, EPP, VRL
      DOUBLE PRECISION :: SUME, SUMKE, SUMPE, SUMTMT, SUMTMR, SUMTMP
      DOUBLE PRECISION :: AVE, AVKE, AVPE, AVTMPT, AVTMPR, AVTMP, EFMT, TIME
      DOUBLE PRECISION :: BSUMPE, BSUMTMP, BSUMPRS, AVPEB, AVTMPB, AVPRSB, DBLE
      DOUBLE PRECISION :: MQ=0.1D0, DH=1.D0 !only NVT
      DOUBLE PRECISION :: TMPINT !initial temperature
      logical :: bool1, bool2
      !double precision :: rand
      !CALL RANDOM_SEED()
      !CALL RANDOM_number( rand )
      !write(*,*) rand


      OPEN (UNIT = 1, FILE = 'parameter1.inp', STATUS = 'UNKNOWN')

!     READ INPUT PARAMETERS

      READ (1, *)
      READ (1, *) NMOL, NSITE
      READ (1, *) 
      READ (1, *) RCUT, RLIST
      READ (1, *)
      READ (1, *) TMPFIX
      READ (1, *)
      READ (1, *) RHO
      READ (1, *)
      READ (1, *) DELT
      READ (1, *)
      READ (1, *) NSTEP, NBATH, NEQ
      READ (1, *) 
      READ (1, *) IDUMP1, IDUMP2, IDUMP3, IDUMP4
 
      CLOSE (1)
      TMPINT = TMPFIX !initial temperature

!     CALCULATE FROM INPUT PARAMETERS

      G      = 6 * NMOL - 3
      NTST   = NMOL * NSITE             ! TOTAL NUMBER OF SITES
      RCUTSQ = RCUT * RCUT
      RLSTSQ = RLIST * RLIST
      EPS4   = 4.D0  ! 4*epsilon
      RCUT6  = (1.D0 / RCUT) ** 6
      RCUT12 = RCUT6 * RCUT6
      CNSTA  = (6.D0 * RCUT12 - 3.D0 * RCUT6) / RCUTSQ
      CNSTB  = 4.D0 * RCUT6 - 7.D0 * RCUT12
      HALFDT = 0.5D0 * DELT
      FCTMQT = 0.5D0*HALFDT/MQ !used only in NVT
 
!     ALLOCATE THE ARRAYS

!      ALLOCATE (R(NMOL,3), QTRN(NMOL,4), V(NMOL,3), W(NMOL,3), F(NMOL,3), T(NMOL,3), P(NMOL,4), &
!      T4(NMOL,4), PHI(NMOL,3), FSITE(NTST,3), MST(NSITE,3), KEP(NMOL), PEP(NMOL), SUMEPT(NMOL), &
!      STAT = STATUS)

!      IF (STATUS /= 0) STOP 'ALLOCATION PROBLEM'
 
!     BOX LENGTH

      VLM    = DBLE(NMOL) / RHO
!      BOXL   = VLM ** (1.D0/3.D0)
      BOXL = 10.0D0
      CALL RBDEFMOL()

      IF ( BOXL .LE. 2.D0 + 2.D0*RCUT ) THEN
        write(*,*) "WARNING: the box is too small, molecules will be interacting with multiple images"
      ENDIF

!      WRITE(*,*) RCUT, RCUT2, RLIST


      INQUIRE(FILE="initpos1.dat", EXIST=bool1 ) 
      INQUIRE(FILE="initortn1.dat", EXIST=bool2 ) 
      if ( bool1 .and. bool2 ) then
        !start from an existing configuration
        !CALL EQCON()
        write(*,*) "reading init files: initpos1.dat and initortn1.dat"
        CALL LOAD_POS_ORTN()
      else
        !     start from scratch
        !     START WITH A SIMPLE CUBIC LATTICE CONFIGURATION
        write(*,*) "starting from a simple cubic lattice"
        CALL RIGIDMD_SC()
        !     SET INITIAL ORIENTATIONS
        CALL INTORN()
      endif

!     SET INITIAL LINEAR and angular VELICITIES
      INQUIRE(FILE="initlinvel1.dat", EXIST=bool1 ) 
      INQUIRE(FILE="initangvel1.dat", EXIST=bool2 ) 
      if ( bool1 .and. bool2 ) then
        write(*,*) "reading init files: initlinvel1.dat and initangvel1.dat"
        CALL LOAD_LIN_ANG_VEL()
      else
        write(*,*) "assigning random velocities"
        CALL INTVEL(TMPINT)
      endif

      !for bug testing only
      IF ( .true. ) then
        OPEN (UNIT=71, FILE='pos0.dat', STATUS='UNKNOWN')
        OPEN (UNIT=72, FILE='ortn0.dat', STATUS='UNKNOWN')
        OPEN (UNIT=73, FILE='lvel0.dat', STATUS='UNKNOWN')
        OPEN (UNIT=74, FILE='avel0.dat', STATUS='UNKNOWN')
        do I=1,NMOL
          WRITE (71, *)  R(I,1), R(I,2), R(I,3)
          WRITE (72, 900)  QTRN(I,1), QTRN(I,2), QTRN(I,3), QTRN(I,4) 
          WRITE (73, *)  V(I,1)/S, V(I,2)/S, V(I,3)/S
          WRITE (74, *) W(I,1), W(I,2), W(I,3)
        enddo
        close(71)
        close(72)
        close(73)
        close(74)
      endif


!     START SIMULATION

      ISTEP  = 0  
      ICNT   = 0
      SUME   = 0.D0
      SUMKE  = 0.D0
      SUMPE  = 0.D0
      SUMTMT = 0.D0
      SUMTMR = 0.D0
      SUMTMP = 0.D0
      SUMPRS = 0.D0
      
      PHI(:,:) = 0.D0
      KE       = 0.D0

!     CALCULATE FORCES AND TORQUES AT TIME t=0

      UPDATE = .TRUE.

      CALL FORQUE(PE, VRL)

      DO I = 1, NMOL

         !calculate the transposed rotation matrix from quaternion
         Q    = QTRN(I,:)
         CALL CONSTRUCT_RMXT( RMXT, Q )

!        TORQUE AT TIME t=0 IN THE BODY FRAME USING THE TRANSPOSE OF THE ROTATION MATRIX

         T4(I,1)   = 0.D0
!         IF(PYT) THEN
!          T(I,:) = MATMUL(RMXT,T(I,:))
!          CALL AA2QTRN(T4(I,:),T(I,:))
!          T4(I,2:4) = MATMUL(RMXT,T4(I,2:4))
!          T4(I,1)   = 0.D0
!          T(I,:) = MATMUL(RMXT,T(I,:))
!          CALL AA2QTRN(T4(I,:),T(I,:))
!          T4(I,2:4) = MATMUL(RMXT,T(I,:))
!          T4(I,1)   = 0.D0
!          T4(I,2:4) = T(I,:)
!         ELSE
          T4(I,2:4) = MATMUL(RMXT,T(I,:))
!         END IF

!     CONSTRUCT THE MATRIX MX FROM THE QUATERNION Q

         CALL CONSTRUCT_MX( MX, Q )

         W4(1)  = 0.D0
         W4(2)  = 2.D0*W(I,1)*JMI(1)
         W4(3)  = 2.D0*W(I,2)*JMI(2)
         W4(4)  = 2.D0*W(I,3)*JMI(3)

         ! intitialize the momentum conjugate to angular velocity
         P(I,:) = MATMUL(MX,W4)
 
         KE = KE + 0.5D0*MASS*(V(I,1)*V(I,1) + V(I,2)*V(I,2) + V(I,3)*V(I,3)) &
                 + 0.5D0*(W(I,1)*W(I,1)/JMI(1) + W(I,2)*W(I,2)/JMI(2) + W(I,3)*W(I,3)/JMI(3))

      ENDDO

      HNOT = KE + PE 

      DO WHILE (ISTEP < NSTEP)            

         !check to make sure the extra variables relevant only for NVT runs are
         !constant and = 1
         !UNIT 0 should be standard error
         IF ( ABS(S-1.D0) > 1.D-7 ) write(0, *) "S != 1: ", S
         IF ( ABS(PS-1.D0) > 1.D-7 ) write(0, *) "PS != 1: ", PS

         ISTEP = ISTEP + 1
         ICNT  = ICNT + 1

!     ADVANCE POSITIONS, ORIENTATIONS, AND THEIR TIME DERIVATIVES
!     MOVE also calls CHECK() and FORQUE
        
         CALL MOVE(ISTEP, NEQ, PE, VRL, TKE, RKE)

!     CALCULATE TEMPERATURES, PRESSURE AND TOTEL ENERGY AT TIME t

         TMPTR   = 2.D0 * TKE  / DBLE(3 * NMOL - 3)
         TMPRT   = 2.D0 * RKE / DBLE(3 * NMOL)           
         KE      = TKE + RKE
         !DH is used only in NVT
         DH      = ABS((KE + PE + 0.5D0*PS*PS/MQ + G*TMPFIX*LOG(S) - HNOT)/HNOT) 
         TMP     = 2.D0 * KE / DBLE(G)
         PEPP    = PE / DBLE(NMOL)
         KEPP    = KE / DBLE(NMOL)
         EPP     = PEPP + KEPP
         PRS     = (2.D0 * TKE + VRL) / (3.D0 * VLM)
         SUME    = SUME + EPP
         SUMKE   = SUMKE + KEPP
         SUMPE   = SUMPE + PEPP
         SUMTMT  = SUMTMT + TMPTR
         SUMTMR  = SUMTMR + TMPRT
         SUMTMP  = SUMTMP + TMP
         SUMPRS  = SUMPRS + PRS

         !do block averaging of potential energy, temperature, and pressure
         IF (ISTEP > NEQ) THEN

            BSUMPE  = BSUMPE + PEPP
            BSUMTMP = BSUMTMP + TMP
            BSUMPRS = BSUMPRS + PRS

            IF (MOD((ISTEP-NEQ),IDUMP4) == 0) THEN

               AVPEB  = BSUMPE / DBLE(IDUMP4)
               AVTMPB = BSUMTMP / DBLE(IDUMP4)
               AVPRSB = BSUMPRS / DBLE(IDUMP4)  

               OPEN (UNIT=25, FILE='blockav1.dat', STATUS='UNKNOWN', POSITION ='APPEND')
                 WRITE (25, *) AVPEB, AVTMPB, AVPRSB
               CLOSE (UNIT=25, STATUS='KEEP')

               BSUMPE  = 0.D0
               BSUMTMP = 0.D0
               BSUMPRS = 0.D0

            ENDIF

         ENDIF
 
         !write out energies, temperatures, and their averages, etc.
         IF (MOD(ISTEP, IDUMP1) .EQ. 0) THEN

            AVE    = SUME / DBLE(ICNT)
            AVKE   = SUMKE / DBLE(ICNT)
            AVPE   = SUMPE / DBLE(ICNT)
            AVTMPT = SUMTMT / DBLE(ICNT)
            AVTMPR = SUMTMR / DBLE(ICNT)
            AVTMP  = SUMTMP / DBLE(ICNT)
            AVPRS  = SUMPRS / DBLE(ICNT)

            OPEN (UNIT=3, FILE='enrg1.dat', STATUS='UNKNOWN', POSITION ='APPEND')
                 WRITE (3, *) KEPP, PEPP, EPP 
            CLOSE (UNIT=3, STATUS='KEEP')
            OPEN (UNIT=4, FILE='tmpprs1.dat', STATUS = 'UNKNOWN', POSITION ='APPEND')
                 WRITE (4, *) TMP, PRS
            CLOSE (UNIT=4, STATUS='KEEP') 

            OPEN (UNIT=33, FILE='avenrg1.dat', STATUS='UNKNOWN', POSITION ='APPEND')
                 WRITE (33, *) AVKE, AVPE, AVE
            CLOSE (UNIT=33, STATUS='KEEP')
            OPEN (UNIT=34, FILE='avtmpprs1.dat', STATUS = 'UNKNOWN', POSITION ='APPEND')
                 WRITE (34, *) AVTMP, AVPRS
            CLOSE (UNIT=34, STATUS='KEEP') 

            OPEN (UNIT=35, FILE='avtmptr1.dat', STATUS = 'UNKNOWN', POSITION ='APPEND')
                 WRITE (35, *) AVTMPT, AVTMPR
            CLOSE (UNIT=35, STATUS='KEEP')

            OPEN (UNIT=36, FILE='error1.dat', STATUS = 'UNKNOWN', POSITION ='APPEND')
                 WRITE (36, *) DH !DH is NVT only
            CLOSE (UNIT=36, STATUS='KEEP')

         ENDIF

         ! calculate energy averaged over time and molecules.
         ! calculate variance over molecules of time averaged energy
         IF (ISTEP > NEQ) THEN

            CALL EFM (EFMT, TIME)

         ENDIF

         !scale velocities. NVE only
         IF (ISTEP <= NBATH .AND. MOD(ISTEP, IDUMP3) == 0 .AND. ISTEP>0) THEN

            CALL NVESCLVEL(AVTMPT, AVTMPR)
            ICNT   = 0
            SUME   = 0.D0
            SUMKE  = 0.D0
            SUMPE  = 0.D0
            SUMTMP = 0.D0
            SUMTMT = 0.D0
            SUMTMR = 0.D0
            SUMPRS = 0.D0

         ENDIF

         !write out positions, orientations, and angular change since t=NEQ
         IF (ISTEP > NEQ .AND. MOD(ISTEP, IDUMP2) == 0) THEN   

           OPEN (UNIT=7,  FILE='pos1.dat', STATUS='UNKNOWN', POSITION='APPEND')
           OPEN (UNIT=8,  FILE='ortn1.dat', STATUS='UNKNOWN', POSITION='APPEND')
           OPEN (UNIT=28, FILE='rot1.dat', STATUS='UNKNOWN', POSITION='APPEND')
           OPEN (UNIT=29, FILE='efm1.dat', STATUS='UNKNOWN', POSITION='APPEND')

           DO I = 1, NMOL

             WRITE (7, *)  R(I,1), R(I,2), R(I,3)
             WRITE (8, 900)  QTRN(I,1), QTRN(I,2), QTRN(I,3), QTRN(I,4) 
             WRITE (28, *) PHI(I,1), PHI(I,2), PHI(I,3)  
             900 FORMAT (1x,4ES25.16)

           ENDDO

           WRITE(29, *) TIME, EFMT

           CLOSE (UNIT=7,  STATUS='KEEP')  
           CLOSE (UNIT=8,  STATUS='KEEP')  
           CLOSE (UNIT=28, STATUS='KEEP')
           CLOSE (UNIT=29, STATUS='KEEP')

           IF (MOD(ISTEP, IDUMP4) .EQ. 0) THEN   

             OPEN (UNIT=9, FILE='lvel1.dat', STATUS='UNKNOWN', POSITION='APPEND')
             OPEN (UNIT=10, FILE='avel1.dat', STATUS='UNKNOWN', POSITION='APPEND')

             DO I = 1, NMOL

               WRITE (9, *)  V(I,1)/S, V(I,2)/S, V(I,3)/S
               WRITE (10, *) W(I,1), W(I,2), W(I,3)

             ENDDO

             CLOSE (UNIT = 9, STATUS = 'KEEP')  
             CLOSE (UNIT = 10, STATUS ='KEEP')  

           ENDIF

         ENDIF
               
         !reinitialize energy averages
         IF (ISTEP == NBATH .OR. ISTEP == NEQ) THEN
              
            !CALL NVESCLVEL(AVTMPT, AVTMPR)
            ICNT   = 0
            SUME   = 0.D0
            SUMKE  = 0.D0
            SUMPE  = 0.D0
            SUMTMT = 0.D0
            SUMTMR = 0.D0
            SUMTMP = 0.D0
            SUMPRS = 0.D0        
            TIME   = 0.D0

            BSUMPE  = 0.D0
            BSUMTMP = 0.D0
            BSUMPRS = 0.D0
                         
         ENDIF  

      ENDDO

      !program ended, print out the final data

      OPEN (UNIT = 21, FILE = 'finalpos1.dat', STATUS = 'UNKNOWN')
      OPEN (UNIT = 22, FILE = 'finalortn1.dat', STATUS = 'UNKNOWN')
      OPEN (UNIT = 23, FILE = 'finallinvel1.dat', STATUS = 'UNKNOWN')
      OPEN (UNIT = 24, FILE = 'finalangvel1.dat', STATUS = 'UNKNOWN')
      
      DO I = 1, NMOL

         WRITE (21, *) R(I,1), R(I,2), R(I,3)
         WRITE (22, 900) QTRN(I,1), QTRN(I,2), QTRN(I,3), QTRN(I,4)
         WRITE (23, *) V(I,1), V(I,2), V(I,3)
         WRITE (24, *) W(I,1), W(I,2), W(I,3)

      ENDDO

      CLOSE(UNIT=21)
      CLOSE(UNIT=22)
      CLOSE(UNIT=23)
      CLOSE(UNIT=24)

!      DEALLOCATE (R, QTRN, V, W, F, T, P, T4, PHI, FSITE, MST, KEP, PEP, SUMEPT, STAT = STATUS)

      call RUN_INFO()

!      END PROGRAM LWOTP_MD_NVE
      END SUBROUTINE RIGIDMD_NVE
!     ------------------------------------------------------------------------------------
 
      SUBROUTINE NVESTEP1 (S, PS, FCTMQT )

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: S, PS
      DOUBLE PRECISION, INTENT(IN) :: FCTMQT

      END SUBROUTINE NVESTEP1

!     ------------------------------------------------------------------------------------
 
      SUBROUTINE NVESTEP2 (PE)

      !STEP 2 in the MOVE algorithm
      ! STEP 2: update V, P, PS
      ! PS = PS - PE*HALFDT   !! PE=potential energy
      ! loop over molecules
      !   update linear momentum (or velocity)
      !   V(I,:) = V(I,:) + S * F(I,:) * HALFDT / MASS   
      !   update rotational momentum 
      !   P(I,:) = P(I,:) +  2*S * DELT/2 * MATMUL(MX,T4(I,:))
      ! end loop over molecules
      
      USE MDCOMMONS, only : NMOL, HALFDT, MASS, QTRN, DELT, T4, P, V, F

      IMPLICIT NONE

      INTEGER :: I
      DOUBLE PRECISION :: MX(4,4), Q(4)
      DOUBLE PRECISION, INTENT(IN) :: PE

      !PS     = PS - PE * HALFDT

      DO I = 1, NMOL

         V(I,:) = V(I,:) + F(I,:) * HALFDT / MASS          

         Q      = QTRN(I,:)

         CALL CONSTRUCT_MX( MX, Q )

         P(I,:)     = P(I,:) +  DELT * MATMUL(MX,T4(I,:))

      ENDDO

      END SUBROUTINE NVESTEP2

!     ------------------------------------------------------------------------------------
 
      SUBROUTINE NVESTEP3 ()

      !STEP 3 in the MOVE algorithm
      ! STEP 3: update QTRN, P, PS with contributions from principle axis 3
      ! ZETASUM = 0
      ! loop over molecules
      !   ZETA(I,3) = dot_product( P(I,:) , MX(:,4) ) / ( 4 * S * JMI(3) )
      !   ZETASUM = ZETASUM + ZETA(I,3)**2
      !   Q(I) = cos( ZETA(I,3)*HALFDT )*Q(I) + sin( ZETA(I,3)*HALFDT )*MX(:,4)
      !   P(I) = cos( ZETA(I,3)*HALFDT )*P(I) + sin( ZETA(I,3)*HALFDT )*MP(:,4)
      !   !   in the above, MP(:,:) is the analogous matrix to MX for vector P(I,:),
      !   !   i.e. MP(:,4) = (/-P(I,4),P(I,3),-P(I,2),P(I,1)/)
      ! end loop over molecules
      ! PS = PS + ZETASUM*JMI(3)*DELT
      ! contribution from principle axis 2
      
      USE MDCOMMONS, only : QTRN, P, HALFDT, JMI, NMOL 

      IMPLICIT NONE

      INTEGER :: I
      DOUBLE PRECISION :: SMPHSQ, Q(4), DQ(4), PHIT
      DOUBLE PRECISION :: PI(4), CP, SP

      SMPHSQ = 0.D0
      DO I = 1, NMOL

         Q         = QTRN(I,:)
         PI        = P(I,:)
         !DQ        = MX(:,4)
         DQ        = (/-Q(4),Q(3),-Q(2),Q(1)/)
         PHIT      = 0.25D0*DOT_PRODUCT(PI,DQ)*HALFDT/(JMI(3))
         !SMPHSQ    = SMPHSQ + PHIT * PHIT
         CP        = COS(PHIT)
         SP        = SIN(PHIT)
         P(I,:)    = CP*PI + SP*(/-PI(4),PI(3),-PI(2),PI(1)/)
         QTRN(I,:) = CP*Q + SP*DQ

      ENDDO

      !PS     = PS + JMI(3)*SMPHSQ*DELT

      END SUBROUTINE NVESTEP3

!     ------------------------------------------------------------------------------------
 
      SUBROUTINE NVESTEP4 ()

      !STEP 4 in the MOVE algorithm
      ! STEP 4: update QTRN, P, PS with contributions from principle axis 2
      ! ZETASUM = 0
      ! loop over molecules
      !   ZETA(I,2) = dot_product( P(I,:) , MX(:,4) ) / ( 4 * S * JMI(2) )
      !   ZETASUM = ZETASUM + ZETA(I,2)**2
      !   Q(I) = cos( ZETA(I,2)*HALFDT )*Q(I) + sin( ZETA(I,2)*HALFDT )*MX(:,2)
      !   P(I) = cos( ZETA(I,2)*HALFDT )*P(I) + sin( ZETA(I,2)*HALFDT )*MP(:,2)
      !   !   in the above, MP(:,:) is the analogous matrix to MX for vector P(I,:),
      !   !   i.e. MP(:,2) = (/-P(I,3),-P(I,4),P(I,1),P(I,2)/)
      ! end loop over molecules
      ! PS = PS + ZETASUM*JMI(2)*DELT
      
      USE MDCOMMONS, only : QTRN, P, HALFDT, JMI, NMOL 

      IMPLICIT NONE

      INTEGER :: I
      DOUBLE PRECISION :: SMPHSQ, Q(4), DQ(4), PHIT
      DOUBLE PRECISION :: PI(4), CP, SP

      SMPHSQ = 0.D0
      DO I = 1, NMOL

         Q         = QTRN(I,:)
         PI        = P(I,:)
         DQ        = (/-Q(3),-Q(4),Q(1),Q(2)/)
         PHIT      = 0.25D0*DOT_PRODUCT(PI,DQ)*HALFDT/(JMI(2))
         SMPHSQ    = SMPHSQ + PHIT*PHIT
         CP        = COS(PHIT)
         SP        = SIN(PHIT)
         P(I,:)    = CP*PI + SP*(/-PI(3),-PI(4),PI(1),PI(2)/)
         QTRN(I,:) = CP*Q + SP*DQ

      ENDDO

      !PS     = PS + JMI(2)*SMPHSQ*DELT

      END SUBROUTINE NVESTEP4

!     ------------------------------------------------------------------------------------
 
      SUBROUTINE NVESTEP5 ()

      !STEP 5 in the MOVE algorithm
      ! STEP 5: update QTRN, P, PS with contributions from principle axis 1
      ! STEP 5 is analogous to STEP 3 and STEP 4 but for ZETA(:,1), i.e. the
      ! first principle axis,
      ! EXCEPT, the update of PS is different and involves the translational
      ! kinetic energy, the temperature (TMPFIX), the number of degrees of
      ! freedom (g), and HNOT,
      
      USE MDCOMMONS, only : QTRN, P, HALFDT, JMI, NMOL, DELT, R, V

      IMPLICIT NONE

      INTEGER :: I
      DOUBLE PRECISION :: SMPHSQ, Q(4), DQ(4), PHIT
      DOUBLE PRECISION :: PI(4), CP, SP
      !DOUBLE PRECISION :: TKE !translational kinetic energy

      SMPHSQ = 0.D0
      !TKE = 0
      DO I = 1, NMOL

         R(I,:)    = R(I,:) + V(I,:) * DELT
         Q         = QTRN(I,:)
         PI        = P(I,:)
         DQ        = (/-Q(2),Q(1),Q(4),-Q(3)/)
         PHIT      = 0.25D0*DOT_PRODUCT(PI,DQ)*DELT/(JMI(1))
         !SMPHSQ    = SMPHSQ + PHIT * PHIT
         CP        = COS(PHIT)
         SP        = SIN(PHIT)
         P(I,:)    = CP*PI + SP*(/-PI(2),PI(1),PI(4),-PI(3)/)
         QTRN(I,:) = CP*Q  + SP*DQ
         !TKE       = TKE + V(I,1)*V(I,1) + V(I,2)*V(I,2) + V(I,3)*V(I,3)

      ENDDO

      !PS = PS + (0.5D0*MASS*TKE/(S*S) + 2.D0*JMI(1)*SMPHSQ - DBLE(G)*TMPFIX*LOG(S) &
         !+ HNOT - DBLE(G)*TMPFIX) * DELT

      END SUBROUTINE NVESTEP5

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE NVESCLVEL(TMPT, TMPR)

! scale linear velocities???

      USE MDCOMMONS, only : TMPFIX, V, W, QTRN, NMOL, JMI, P

      IMPLICIT NONE

      INTEGER          :: I
      DOUBLE PRECISION :: SCLTR, SCLRT, TMPT, TMPR
      DOUBLE PRECISION :: Q(4), RMXT(3,3), MX(4,4), W4(4)

      SCLTR   = DSQRT(TMPFIX / TMPT)
      SCLRT   = DSQRT(TMPFIX / TMPR)

      DO I = 1, NMOL

         V(I,:) = V(I,:) * SCLTR

         W(I,:) = W(I,:) * SCLRT

         Q    = QTRN(I,:)
         CALL CONSTRUCT_RMXT(RMXT,Q)

!     CONSTRUCT THE MATRIX MX FROM THE QUATERNION Q

         CALL CONSTRUCT_MX(MX,Q)

         W4(1)  = 0.D0
         W4(2)  = 2.D0*JMI(1)*W(I,1)
         W4(3)  = 2.D0*JMI(2)*W(I,2)
         W4(4)  = 2.D0*JMI(3)*W(I,3)

         P(I,:) = MATMUL(MX,W4)

      ENDDO

      END SUBROUTINE NVESCLVEL
!
! svn info 
! $Date: 2011-09-21 11:24:09 +0100 (Wed, 21 Sep 2011) $
! $Rev: 25 $
!

! the lewis wahnstrom model has 
      ! sigma = 0.483 nm
      ! epsilon = 5.276 kJ/mol
 
!     ----------------------------------------------------------------------------------------------

!      PROGRAM LWOTP_MD_NVE
      SUBROUTINE RIGIDMD_NVE_CHANGE_DENSITY
      USE MDCOMMONS

      IMPLICIT NONE

! IDUMP1: dump files.  always: 
!                     enrg1.dat: KEPP, PEPP, EPP: kinetic, potential, energy per molecule
!                     tmpprs1.dat: temp, pressure?
!                     avenrg1.dat: average KE, PE, E
!                     avtmpprs1.dat: average temp, pressure?
!                     avtmptr1.dat: average tranlational, rotational temperature
! IDUMP2: dump files.  ISTEP > NEQ
!                     pos1.dat:  xyz positions of CoM of molecules
!                     ortn1.dat: quaternions giving the orientation of molecules
!                     rot1.dat:  phi1, phi2, phi3 ????
!                     efm1.dat:  energy fluctuation metric: time from NEQ, and variance over molecules of time averaged energy 
! IDUMP3: CALL SCLVEL.  ISTEP <= NBATH: scale linear velocities ???
! IDUMP4: dump files.  ISTEP > NEQ: 
!                     blockav1.dat: average PE, temp, pressure
!                     lvel1.dat:    linear velocities
!                     avel1.dat:    angular velocities
! TKE: translational kinetic energy
! RKE: rotational kinetic energy
! MX:  is a representation of the quaternion as a matrix in such a way
!      that quaternion addition and multiplication correspond to matrix
!      addition and matrix multiplication
!      In this way the conjugate of the quaternion corresponds to the
!      transpose of the matrix.  The fourth power of the norm of a
!      quaternion is the determinant of the corresponding matrix.
!      Complex numbers are block diagonal matrices with two 2x2 blocks.
! EFMT: Energy Fluctuation Metric: the variance over molecules of the
!       time averaged energy
      INTEGER          :: NSTEP, NBATH, NEQ, ISTEP, ICNT, IDUMP1, IDUMP2, IDUMP3, IDUMP4, I!, K!, STATUS
      DOUBLE PRECISION :: RCUT6, RCUT12, TMP, TMPTR, TMPRT, RHO, VLM
      DOUBLE PRECISION :: Q(4), RMXT(3,3), MX(4,4), W4(4)
      DOUBLE PRECISION :: PRS, SUMPRS, AVPRS, PE, KE, TKE, RKE, PEPP, KEPP, EPP, VRL
      DOUBLE PRECISION :: SUME, SUMKE, SUMPE, SUMTMT, SUMTMR, SUMTMP
      DOUBLE PRECISION :: AVE, AVKE, AVPE, AVTMPT, AVTMPR, AVTMP, EFMT, TIME
      DOUBLE PRECISION :: BSUMPE, BSUMTMP, BSUMPRS, AVPEB, AVTMPB, AVPRSB
      DOUBLE PRECISION :: MQ=0.1D0, DH=1.D0 !only NVT
      DOUBLE PRECISION :: TMPINT !initial temperature
      logical :: bool1, bool2
      !CALL RANDOM_SEED()
      DOUBLE PRECISION :: RHOFINAL, DELRHO, RHOOLD, LCHANGE, DBLE

      OPEN (UNIT = 1, FILE = 'parameter1.inp', STATUS = 'UNKNOWN')

!     READ INPUT PARAMETERS

      READ (1, *)
      READ (1, *) NMOL, NSITE
      READ (1, *) 
      READ (1, *) RCUT, RLIST
      READ (1, *)
      READ (1, *) TMPFIX
      READ (1, *)
      READ (1, *) RHO, RHOFINAL, DELRHO
      READ (1, *)
      READ (1, *) DELT
      READ (1, *)
      READ (1, *) NSTEP, NBATH, NEQ
      READ (1, *) 
      READ (1, *) IDUMP1, IDUMP2, IDUMP3, IDUMP4
 
      CLOSE (1)
      TMPINT = TMPFIX !initial temperature

!     CALCULATE FROM INPUT PARAMETERS

      G      = 6 * NMOL - 3
      NTST   = NMOL * NSITE             ! TOTAL NUMBER OF SITES
      RCUTSQ = RCUT * RCUT
      RLSTSQ = RLIST * RLIST
      EPS4   = 4.D0  ! 4*epsilon
      RCUT6  = (1.D0 / RCUT) ** 6
      RCUT12 = RCUT6 * RCUT6
      CNSTA  = (6.D0 * RCUT12 - 3.D0 * RCUT6) / RCUTSQ
      CNSTB  = 4.D0 * RCUT6 - 7.D0 * RCUT12
      HALFDT = 0.5D0 * DELT
      FCTMQT = 0.5D0*HALFDT/MQ !used only in NVT
 
!     ALLOCATE THE ARRAYS

!      ALLOCATE (R(NMOL,3), QTRN(NMOL,4), V(NMOL,3), W(NMOL,3), F(NMOL,3), T(NMOL,3), P(NMOL,4), &
!      T4(NMOL,4), PHI(NMOL,3), FSITE(NTST,3), MST(NSITE,3), KEP(NMOL), PEP(NMOL), SUMEPT(NMOL), &
!      STAT = STATUS)

!      IF (STATUS /= 0) STOP 'ALLOCATION PROBLEM'
 
!     BOX LENGTH

      VLM    = DBLE(NMOL) / RHO
      BOXL   = VLM ** (1.D0/3.D0)

      CALL RBDEFMOL()

!      WRITE(*,*) RCUT, RCUT2, RLIST


      INQUIRE(FILE="initpos1.dat", EXIST=bool1 ) 
      INQUIRE(FILE="initortn1.dat", EXIST=bool2 ) 
      if ( bool1 .and. bool2 ) then
        !start from an existing configuration
        !CALL EQCON()
        write(*,*) " reading init files"
        CALL LOAD_POS_ORTN()
      else
        !     start from scratch
        !     START WITH A SIMPLE CUBIC LATTICE CONFIGURATION
        write(*,*) " no init files: starting from scratch"
        CALL RIGIDMD_SC()
        !     SET INITIAL ORIENTATIONS
        CALL INTORN()
      endif

!     SET INITIAL LINEAR and angular VELICITIES
      CALL INTVEL(TMPINT)


!     START SIMULATION

      ISTEP  = 0  
      ICNT   = 0
      SUME   = 0.D0
      SUMKE  = 0.D0
      SUMPE  = 0.D0
      SUMTMT = 0.D0
      SUMTMR = 0.D0
      SUMTMP = 0.D0
      SUMPRS = 0.D0
      
      PHI(:,:) = 0.D0
      KE       = 0.D0

!     CALCULATE FORCES AND TORQUES AT TIME t=0

      UPDATE = .TRUE.

      CALL FORQUE(PE, VRL)

      DO I = 1, NMOL

         !calculate the transposed rotation matrix from quaternion
         Q    = QTRN(I,:)
         CALL CONSTRUCT_RMXT( RMXT, Q )

!        TORQUE AT TIME t=0 IN THE BODY FRAME USING THE TRANSPOSE OF THE ROTATION MATRIX

         T4(I,1)   = 0.D0
         T4(I,2:4) = MATMUL(RMXT,T(I,:))

!     CONSTRUCT THE MATRIX MX FROM THE QUATERNION Q

         CALL CONSTRUCT_MX( MX, Q )

         W4(1)  = 0.D0
         W4(2)  = 2.D0*W(I,1)
         W4(3)  = 2.D0*W(I,2)
         W4(4)  = 2.D0*W(I,3)

         ! intitialize the momentum conjugate to angular velocity
         P(I,:) = MATMUL(MX,W4)
 
         KE = KE + 0.5D0*MASS*(V(I,1)*V(I,1) + V(I,2)*V(I,2) + V(I,3)*V(I,3)) &
                 + 0.5D0*(W(I,1)*W(I,1)/JMI(1) + W(I,2)*W(I,2)/JMI(2) + W(I,3)*W(I,3)/JMI(3))

      ENDDO

      HNOT = KE + PE 

      DO WHILE (ISTEP < NSTEP)            

         !check to make sure the extra variables relevant only for NVT runs are
         !constant and = 1
         !UNIT 0 should be standard error
         IF ( ABS(S-1.D0) > 1.D-7 ) write(0, *) "S != 1: ", S
         IF ( ABS(PS-1.D0) > 1.D-7 ) write(0, *) "PS != 1: ", PS

         ISTEP = ISTEP + 1
         ICNT  = ICNT + 1

!     ADVANCE POSITIONS, ORIENTATIONS, AND THEIR TIME DERIVATIVES
!     MOVE also calls CHECK() and FORQUE
        
         CALL MOVE(ISTEP, NEQ, PE, VRL, TKE, RKE)

!     CALCULATE TEMPERATURES, PRESSURE AND TOTEL ENERGY AT TIME t

         TMPTR   = 2.D0 * TKE  / DBLE(3 * NMOL - 3)
         TMPRT   = 2.D0 * RKE / DBLE(3 * NMOL)           
         KE      = TKE + RKE
         !DH is used only in NVT
         DH      = ABS((KE + PE + 0.5D0*PS*PS/MQ + G*TMPFIX*LOG(S) - HNOT)/HNOT) 
         TMP     = 2.D0 * KE / DBLE(G)
         PEPP    = PE / DBLE(NMOL)
         KEPP    = KE / DBLE(NMOL)
         EPP     = PEPP + KEPP
         PRS     = (2.D0 * TKE + VRL) / (3.D0 * VLM)
         SUME    = SUME + EPP
         SUMKE   = SUMKE + KEPP
         SUMPE   = SUMPE + PEPP
         SUMTMT  = SUMTMT + TMPTR
         SUMTMR  = SUMTMR + TMPRT
         SUMTMP  = SUMTMP + TMP
         SUMPRS  = SUMPRS + PRS

         !do block averaging of potential energy, temperature, and pressure
         IF (ISTEP > NEQ) THEN

            BSUMPE  = BSUMPE + PEPP
            BSUMTMP = BSUMTMP + TMP
            BSUMPRS = BSUMPRS + PRS

            IF (MOD((ISTEP-NEQ),IDUMP4) == 0) THEN

               AVPEB  = BSUMPE / DBLE(IDUMP4)
               AVTMPB = BSUMTMP / DBLE(IDUMP4)
               AVPRSB = BSUMPRS / DBLE(IDUMP4)  

               OPEN (UNIT=25, FILE='blockav1.dat', STATUS='UNKNOWN', POSITION='APPEND')
                 WRITE (25, *) AVPEB, AVTMPB, AVPRSB
               CLOSE (UNIT=25, STATUS='KEEP')

               BSUMPE  = 0.D0
               BSUMTMP = 0.D0
               BSUMPRS = 0.D0

            ENDIF

         ENDIF
 
         !write out energies, temperatures, and their averages, etc.
         IF (MOD(ISTEP, IDUMP1) .EQ. 0) THEN

            AVE    = SUME / DBLE(ICNT)
            AVKE   = SUMKE / DBLE(ICNT)
            AVPE   = SUMPE / DBLE(ICNT)
            AVTMPT = SUMTMT / DBLE(ICNT)
            AVTMPR = SUMTMR / DBLE(ICNT)
            AVTMP  = SUMTMP / DBLE(ICNT)
            AVPRS  = SUMPRS / DBLE(ICNT)

            OPEN (UNIT=3, FILE='enrg1.dat', STATUS='UNKNOWN', POSITION='APPEND')
                 WRITE (3, *) KEPP, PEPP, EPP 
            CLOSE (UNIT=3, STATUS='KEEP')
            OPEN (UNIT=4, FILE='tmpprs1.dat', STATUS = 'UNKNOWN', POSITION='APPEND')
                 WRITE (4, *) TMP, PRS
            CLOSE (UNIT=4, STATUS='KEEP') 

            OPEN (UNIT=33, FILE='avenrg1.dat', STATUS='UNKNOWN', POSITION='APPEND')
                 WRITE (33, *) AVKE, AVPE, AVE
            CLOSE (UNIT=33, STATUS='KEEP')
            OPEN (UNIT=34, FILE='avtmpprs1.dat', STATUS = 'UNKNOWN', POSITION='APPEND')
                 WRITE (34, *) AVTMP, AVPRS
            CLOSE (UNIT=34, STATUS='KEEP') 

            OPEN (UNIT=35, FILE='avtmptr1.dat', STATUS = 'UNKNOWN', POSITION='APPEND')
                 WRITE (35, *) AVTMPT, AVTMPR
            CLOSE (UNIT=35, STATUS='KEEP')

            OPEN (UNIT=36, FILE='error1.dat', STATUS = 'UNKNOWN', POSITION='APPEND')
                 WRITE (36, *) DH !DH is NVT only
            CLOSE (UNIT=36, STATUS='KEEP')

         ENDIF

         ! calculate energy averaged over time and molecules.
         ! calculate variance over molecules of time averaged energy
         IF (ISTEP > NEQ) THEN

            CALL EFM (EFMT, TIME)

         ENDIF

         !scale velocities. NVE only
         IF ( MOD(ISTEP, IDUMP3) == 0) THEN

             RHOOLD=RHO
           IF ( RHOFINAL .gt. RHO+DELRHO ) THEN
             RHO=RHO+DELRHO
           ELSE IF ( RHOFINAL .lt. RHO-DELRHO ) THEN
             RHO=RHO-DELRHO
           ELSE
             RHO=RHOFINAL
           ENDIF
           if ( RHOOLD .ne. RHO ) then
             VLM    = DBLE(NMOL) / RHO
             BOXL   = VLM ** (1.D0/3.D0)

             LCHANGE=(RHOOLD/RHO)**(1.D0/3.D0)
             write (*,*) "RHO ", rho, " LCHANGE ", LCHANGE
             R = R*LCHANGE;
           endif


            CALL NVESCLVEL(AVTMPT, AVTMPR)
            ICNT   = 0
            SUME   = 0.D0
            SUMKE  = 0.D0
            SUMPE  = 0.D0
            SUMTMP = 0.D0
            SUMTMT = 0.D0
            SUMTMR = 0.D0
            SUMPRS = 0.D0

         ENDIF

         !write out positions, orientations, and angular change since t=NEQ
         IF (ISTEP > NEQ .AND. MOD(ISTEP, IDUMP2) == 0) THEN   

           OPEN (UNIT=7,  FILE='pos1.dat', STATUS='UNKNOWN', POSITION='APPEND')
           OPEN (UNIT=8,  FILE='ortn1.dat', STATUS='UNKNOWN', POSITION='APPEND')
           !OPEN (UNIT=28, FILE='rot1.dat', STATUS='UNKNOWN', POSITION='APPEND')
           OPEN (UNIT=29, FILE='efm1.dat', STATUS='UNKNOWN', POSITION='APPEND')

           DO I = 1, NMOL

             WRITE (7, *)  R(I,1), R(I,2), R(I,3)
             WRITE (8, 900)  QTRN(I,1), QTRN(I,2), QTRN(I,3), QTRN(I,4) 
             !WRITE (28, *) PHI(I,1), PHI(I,2), PHI(I,3)  
             900 FORMAT (1x,ES20.12,1x,ES20.12,1x,ES20.12)

           ENDDO

           WRITE(29, *) TIME, EFMT

           CLOSE (UNIT=7,  STATUS='KEEP')  
           CLOSE (UNIT=8,  STATUS='KEEP')  
           !CLOSE (UNIT=28, STATUS='KEEP')
           CLOSE (UNIT=29, STATUS='KEEP')

           IF (MOD(ISTEP, IDUMP4) .EQ. 0) THEN   

             OPEN (UNIT=9, FILE='lvel1.dat', STATUS='UNKNOWN', POSITION='APPEND')
             OPEN (UNIT=10, FILE='avel1.dat', STATUS='UNKNOWN', POSITION='APPEND')

             DO I = 1, NMOL

               WRITE (9, *)  V(I,1)/S, V(1,2)/S, V(1,3)/S
               WRITE (10, *) W(I,1), W(I,2), W(I,3)

             ENDDO

             CLOSE (UNIT = 9, STATUS = 'KEEP')  
             CLOSE (UNIT = 10, STATUS ='KEEP')  

           ENDIF

         ENDIF
               
         !reinitialize energy averages
         IF (ISTEP == NBATH .OR. ISTEP == NEQ) THEN
              
            ICNT   = 0
            SUME   = 0.D0
            SUMKE  = 0.D0
            SUMPE  = 0.D0
            SUMTMT = 0.D0
            SUMTMR = 0.D0
            SUMTMP = 0.D0
            SUMPRS = 0.D0        
            TIME   = 0.D0

            BSUMPE  = 0.D0
            BSUMTMP = 0.D0
            BSUMPRS = 0.D0
                         
         ENDIF  

      ENDDO

      !program ended, print out the final data

      OPEN (UNIT = 21, FILE = 'finalpos1.dat', STATUS = 'UNKNOWN')
      OPEN (UNIT = 22, FILE = 'finalortn1.dat', STATUS = 'UNKNOWN')
      OPEN (UNIT = 23, FILE = 'finallinvel1.dat', STATUS = 'UNKNOWN')
      OPEN (UNIT = 24, FILE = 'finalangvel1.dat', STATUS = 'UNKNOWN')
      
      DO I = 1, NMOL

         WRITE (21, *) R(I,1), R(I,2), R(I,3)
         WRITE (22, 900) QTRN(I,1), QTRN(I,2), QTRN(I,3), QTRN(I,4)
         WRITE (23, *) V(I,1), V(I,2), V(I,3)
         WRITE (24, *) W(I,1), W(I,2), W(I,3)

      ENDDO

      CLOSE(UNIT=21)
      CLOSE(UNIT=22)
      CLOSE(UNIT=23)
      CLOSE(UNIT=24)

!      DEALLOCATE (R, QTRN, V, W, F, T, P, T4, PHI, FSITE, MST, KEP, PEP, SUMEPT, STAT = STATUS)

      call RUN_INFO()

!      END PROGRAM LWOTP_MD_NVE
      END SUBROUTINE RIGIDMD_NVE_CHANGE_DENSITY

!     ------------------------------------------------------------------------------------
 
      SUBROUTINE NVECHDSTEP1 (S, PS, FCTMQT )

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(IN) :: S, PS
      DOUBLE PRECISION, INTENT(IN) :: FCTMQT

      END SUBROUTINE NVECHDSTEP1

!     ------------------------------------------------------------------------------------
 
      SUBROUTINE NVECHDSTEP2 (PE)

      !STEP 2 in the MOVE algorithm
      ! STEP 2: update V, P, PS
      ! PS = PS - PE*HALFDT   !! PE=potential energy
      ! loop over molecules
      !   update linear momentum (or velocity)
      !   V(I,:) = V(I,:) + S * F(I,:) * HALFDT / MASS   
      !   update rotational momentum 
      !   P(I,:) = P(I,:) +  2*S * DELT/2 * MATMUL(MX,T4(I,:))
      ! end loop over molecules
      
      USE MDCOMMONS, only : NMOL, HALFDT, MASS, QTRN, DELT, T4, P, V, F

      IMPLICIT NONE

      INTEGER :: I
      DOUBLE PRECISION :: MX(4,4), Q(4)
      DOUBLE PRECISION, INTENT(IN) :: PE

      !PS     = PS - PE * HALFDT

      DO I = 1, NMOL

         V(I,:) = V(I,:) + F(I,:) * HALFDT / MASS          

         Q      = QTRN(I,:)

         CALL CONSTRUCT_MX( MX, Q )

         P(I,:)     = P(I,:) +  DELT * MATMUL(MX,T4(I,:))

      ENDDO

      END SUBROUTINE NVECHDSTEP2

!     ------------------------------------------------------------------------------------
 
      SUBROUTINE NVECHDSTEP3 ()

      !STEP 3 in the MOVE algorithm
      ! STEP 3: update QTRN, P, PS with contributions from principle axis 3
      ! ZETASUM = 0
      ! loop over molecules
      !   ZETA(I,3) = dot_product( P(I,:) , MX(:,4) ) / ( 4 * S * JMI(3) )
      !   ZETASUM = ZETASUM + ZETA(I,3)**2
      !   Q(I) = cos( ZETA(I,3)*HALFDT )*Q(I) + sin( ZETA(I,3)*HALFDT )*MX(:,4)
      !   P(I) = cos( ZETA(I,3)*HALFDT )*P(I) + sin( ZETA(I,3)*HALFDT )*MP(:,4)
      !   !   in the above, MP(:,:) is the analogous matrix to MX for vector P(I,:),
      !   !   i.e. MP(:,4) = (/-P(I,4),P(I,3),-P(I,2),P(I,1)/)
      ! end loop over molecules
      ! PS = PS + ZETASUM*JMI(3)*DELT
      ! contribution from principle axis 2
      
      USE MDCOMMONS, only : QTRN, P, HALFDT, JMI, NMOL 

      IMPLICIT NONE

      INTEGER :: I
      DOUBLE PRECISION :: SMPHSQ, Q(4), DQ(4), PHIT
      DOUBLE PRECISION :: PI(4), CP, SP

      SMPHSQ = 0.D0
      DO I = 1, NMOL

         Q         = QTRN(I,:)
         PI        = P(I,:)
         !DQ        = MX(:,4)
         DQ        = (/-Q(4),Q(3),-Q(2),Q(1)/)
         PHIT      = 0.25D0*DOT_PRODUCT(PI,DQ)*HALFDT/(JMI(3))
         !SMPHSQ    = SMPHSQ + PHIT * PHIT
         CP        = COS(PHIT)
         SP        = SIN(PHIT)
         P(I,:)    = CP*PI + SP*(/-PI(4),PI(3),-PI(2),PI(1)/)
         QTRN(I,:) = CP*Q + SP*DQ

      ENDDO

      !PS     = PS + JMI(3)*SMPHSQ*DELT

      END SUBROUTINE NVECHDSTEP3

!     ------------------------------------------------------------------------------------
 
      SUBROUTINE NVECHDSTEP4 ()

      !STEP 4 in the MOVE algorithm
      ! STEP 4: update QTRN, P, PS with contributions from principle axis 2
      ! ZETASUM = 0
      ! loop over molecules
      !   ZETA(I,2) = dot_product( P(I,:) , MX(:,4) ) / ( 4 * S * JMI(2) )
      !   ZETASUM = ZETASUM + ZETA(I,2)**2
      !   Q(I) = cos( ZETA(I,2)*HALFDT )*Q(I) + sin( ZETA(I,2)*HALFDT )*MX(:,2)
      !   P(I) = cos( ZETA(I,2)*HALFDT )*P(I) + sin( ZETA(I,2)*HALFDT )*MP(:,2)
      !   !   in the above, MP(:,:) is the analogous matrix to MX for vector P(I,:),
      !   !   i.e. MP(:,2) = (/-P(I,3),-P(I,4),P(I,1),P(I,2)/)
      ! end loop over molecules
      ! PS = PS + ZETASUM*JMI(2)*DELT
      
      USE MDCOMMONS, only : QTRN, P, HALFDT, JMI, NMOL 

      IMPLICIT NONE

      INTEGER :: I
      DOUBLE PRECISION :: SMPHSQ, Q(4), DQ(4), PHIT
      DOUBLE PRECISION :: PI(4), CP, SP

      SMPHSQ = 0.D0
      DO I = 1, NMOL

         Q         = QTRN(I,:)
         PI        = P(I,:)
         DQ        = (/-Q(3),-Q(4),Q(1),Q(2)/)
         PHIT      = 0.25D0*DOT_PRODUCT(PI,DQ)*HALFDT/(JMI(2))
         SMPHSQ    = SMPHSQ + PHIT*PHIT
         CP        = COS(PHIT)
         SP        = SIN(PHIT)
         P(I,:)    = CP*PI + SP*(/-PI(3),-PI(4),PI(1),PI(2)/)
         QTRN(I,:) = CP*Q + SP*DQ

      ENDDO

      !PS     = PS + JMI(2)*SMPHSQ*DELT

      END SUBROUTINE NVECHDSTEP4

!     ------------------------------------------------------------------------------------
 
      SUBROUTINE NVECHDSTEP5 ()

      !STEP 5 in the MOVE algorithm
      ! STEP 5: update QTRN, P, PS with contributions from principle axis 1
      ! STEP 5 is analogous to STEP 3 and STEP 4 but for ZETA(:,1), i.e. the
      ! first principle axis,
      ! EXCEPT, the update of PS is different and involves the translational
      ! kinetic energy, the temperature (TMPFIX), the number of degrees of
      ! freedom (g), and HNOT,
      
      USE MDCOMMONS, only : QTRN, P, HALFDT, JMI, NMOL, DELT, R, V

      IMPLICIT NONE

      INTEGER :: I
      DOUBLE PRECISION :: SMPHSQ, Q(4), DQ(4), PHIT
      DOUBLE PRECISION :: PI(4), CP, SP
      !DOUBLE PRECISION :: TKE !translational kinetic energy

      SMPHSQ = 0.D0
      !TKE = 0
      DO I = 1, NMOL

         R(I,:)    = R(I,:) + V(I,:) * DELT
         Q         = QTRN(I,:)
         PI        = P(I,:)
         DQ        = (/-Q(2),Q(1),Q(4),-Q(3)/)
         PHIT      = 0.25D0*DOT_PRODUCT(PI,DQ)*DELT/(JMI(1))
         !SMPHSQ    = SMPHSQ + PHIT * PHIT
         CP        = COS(PHIT)
         SP        = SIN(PHIT)
         P(I,:)    = CP*PI + SP*(/-PI(2),PI(1),PI(4),-PI(3)/)
         QTRN(I,:) = CP*Q  + SP*DQ
         !TKE       = TKE + V(I,1)*V(I,1) + V(I,2)*V(I,2) + V(I,3)*V(I,3)

      ENDDO

      !PS = PS + (0.5D0*MASS*TKE/(S*S) + 2.D0*JMI(1)*SMPHSQ - DBLE(G)*TMPFIX*LOG(S) &
         !+ HNOT - DBLE(G)*TMPFIX) * DELT

      END SUBROUTINE NVECHDSTEP5

!     ----------------------------------------------------------------------------------------------

      SUBROUTINE NVECHDSCLVEL(TMPT, TMPR)

! scale linear velocities???

      USE MDCOMMONS, only : TMPFIX, V, W, QTRN, NMOL, JMI, P

      IMPLICIT NONE

      INTEGER          :: I
      DOUBLE PRECISION :: SCLTR, SCLRT, TMPT, TMPR
      DOUBLE PRECISION :: Q(4), RMXT(3,3), MX(4,4), W4(4)

      SCLTR   = DSQRT(TMPFIX / TMPT)
      SCLRT   = DSQRT(TMPFIX / TMPR)

      DO I = 1, NMOL

         V(I,:) = V(I,:) * SCLTR

         W(I,:) = W(I,:) * SCLRT

         Q    = QTRN(I,:)
         CALL CONSTRUCT_RMXT(RMXT,Q)

!     CONSTRUCT THE MATRIX MX FROM THE QUATERNION Q

         CALL CONSTRUCT_MX(MX,Q)

         W4(1)  = 0.D0
         W4(2)  = 2.D0*JMI(1)*W(I,1)
         W4(3)  = 2.D0*JMI(2)*W(I,2)
         W4(4)  = 2.D0*JMI(3)*W(I,3)

         P(I,:) = MATMUL(MX,W4)

      ENDDO

      END SUBROUTINE NVECHDSCLVEL
!
! svn info 
! $Date: 2011-10-16 20:53:28 +0100 (Sun, 16 Oct 2011) $
! $Rev: 61 $
!

!     ----------------------------------------------------------------------------------------------

!      PROGRAM LWOTP_MD_NVT
      SUBROUTINE RIGIDMD_NVT
      USE MDCOMMONS

      IMPLICIT NONE

! IDUMP1: dump files.  always: 
!                     enrg1.dat: KEPP, PEPP, EPP: kinetic, potential, energy per molecule
!                     tmpprs1.dat: temp, pressure?
!                     avenrg1.dat: average KE, PE, E
!                     avtmpprs1.dat: average temp, pressure?
!                     avtmptr1.dat: average tranlational, rotational temperature
! IDUMP2: dump files.  ISTEP > NEQ
!                     pos1.dat:  xyz positions of CoM of molecules
!                     ortn1.dat: quaternions giving the orientation of molecules
!                     rot1.dat:  phi1, phi2, phi3 ????
!                     efm1.dat:  energy fluctuation metric: time from NEQ, and variance over molecules of time averaged energy 
! IDUMP3: CALL SCLVEL.  ISTEP <= NBATH: scale linear velocities ???
! IDUMP4: dump files.  ISTEP > NEQ: 
!                     blockav1.dat: average PE, temp, pressure
!                     lvel1.dat:    linear velocities
!                     avel1.dat:    angular velocities
! TKE: translational kinetic energy
! RKE: rotational kinetic energy
! MX:  is a representation of the quaternion as a matrix in such a way
!      that quaternion addition and multiplication correspond to matrix
!      addition and matrix multiplication
!      In this way the conjugate of the quaternion corresponds to the
!      transpose of the matrix.  The fourth power of the norm of a
!      quaternion is the determinant of the corresponding matrix.
!      Complex numbers are block diagonal matrices with two 2x2 blocks.
! EFMT: Energy Fluctuation Metric: the variance over molecules of the
!       time averaged energy
! VRL: accumulates the forces dotted into the positions Fij.dot(Ri - Rj) which
!      is used to calculate the pressure
      INTEGER          :: NSTEP, NBATH, NEQ, ISTEP, ICNT, IDUMP1, IDUMP2, IDUMP4, I!, K!, STATUS
      DOUBLE PRECISION :: RCUT6, RCUT12, TMP, TMPTR, TMPRT, RHO, VLM
      DOUBLE PRECISION :: Q(4), RMXT(3,3)!, MX(4,4), W4(4)
      DOUBLE PRECISION :: PRS, SUMPRS, AVPRS, PE, KE, TKE, RKE, PEPP, KEPP, EPP, VRL
      DOUBLE PRECISION :: SUME, SUMKE, SUMPE, SUMTMT, SUMTMR, SUMTMP
      DOUBLE PRECISION :: AVE, AVKE, AVPE, AVTMPT, AVTMPR, AVTMP, EFMT, TIME
      DOUBLE PRECISION :: BSUMPE, BSUMTMP, BSUMPRS, AVPEB, AVTMPB, AVPRSB, DBLE
      DOUBLE PRECISION :: MQ=0.1D0, DH=1.D0 !only NVT
      DOUBLE PRECISION :: TMPINT !initial temperature
      logical :: bool1, bool2, loadHNOT
      !double precision :: rand
      !CALL RANDOM_SEED()
      !CALL RANDOM_number( rand )
      !write(*,*) rand
      !DOUBLE PRECISION :: dummyvec3(3), dummyvec4(4)


      OPEN (UNIT = 1, FILE = 'parameternvt1.inp', STATUS = 'UNKNOWN')

!     READ INPUT PARAMETERS

      READ (1, *)
      READ (1, *) NMOL, NSITE
      READ (1, *) 
      READ (1, *) RCUT, RLIST
      READ (1, *)
      READ (1, *) TMPFIX
      READ (1, *)
      READ (1, *) RHO
      READ (1, *)
      READ (1, *) DELT
      READ (1, *)
      READ (1, *) S, PS, MQ
      READ (1, *) 
      READ (1, *) NSTEP, NBATH, NEQ
      READ (1, *) 
      READ (1, *) IDUMP1, IDUMP2, IDUMP4
 
      CLOSE (1)
      TMPINT = TMPFIX !initial temperature

!     CALCULATE FROM INPUT PARAMETERS

      G      = 6 * NMOL - 3
      NTST   = NMOL * NSITE             ! TOTAL NUMBER OF SITES
      RCUTSQ = RCUT * RCUT
      RLSTSQ = RLIST * RLIST
      EPS4   = 4.D0  ! 4*epsilon
      RCUT6  = (1.D0 / RCUT) ** 6
      RCUT12 = RCUT6 * RCUT6
      CNSTA  = (6.D0 * RCUT12 - 3.D0 * RCUT6) / RCUTSQ
      CNSTB  = 4.D0 * RCUT6 - 7.D0 * RCUT12
      HALFDT = 0.5D0 * DELT
      FCTMQT = 0.5D0*HALFDT/MQ !used only in NVT
 
!     ALLOCATE THE ARRAYS

!      ALLOCATE (R(NMOL,3), QTRN(NMOL,4), V(NMOL,3), W(NMOL,3), F(NMOL,3), T(NMOL,3), P(NMOL,4), &
!      T4(NMOL,4), PHI(NMOL,3), FSITE(NTST,3), MST(NSITE,3), KEP(NMOL), PEP(NMOL), SUMEPT(NMOL), &
!      STAT = STATUS)

!      IF (STATUS /= 0) STOP 'ALLOCATION PROBLEM'
 
!     BOX LENGTH

      VLM    = DBLE(NMOL) / RHO
      BOXL   = VLM ** (1.D0/3.D0)

      CALL RBDEFMOL()

!      WRITE(*,*) RCUT, RCUT2, RLIST
      WRITE(*,*) "Read in input parameters"

      INQUIRE(FILE="initpos1.dat", EXIST=bool1 ) 
      INQUIRE(FILE="initortn1.dat", EXIST=bool2 ) 
      if ( bool1 .and. bool2 ) then
        !start from an existing configuration
        !CALL EQCON()
        write(*,*) "reading init files: initpos1.dat and initortn1.dat"
        CALL LOAD_POS_ORTN()
      else
        !     start from scratch
        !     START WITH A SIMPLE CUBIC LATTICE CONFIGURATION
        write(*,*) "starting from a simple cubic lattice"
        CALL RIGIDMD_SC()
        !     SET INITIAL ORIENTATIONS
        CALL INTORN()
      endif

!     SET INITIAL LINEAR and angular VELICITIES
      INQUIRE(FILE="initlinvel1.dat", EXIST=bool1 ) 
      INQUIRE(FILE="initangvel1.dat", EXIST=bool2 ) 
      if ( bool1 .and. bool2 ) then
        write(*,*) "reading init files: initlinvel1.dat and initangvel1.dat"
        CALL LOAD_LIN_ANG_VEL()
      else
        write(*,*) "assigning random velocities"
        CALL INTVEL(TMPINT)
      endif
      V=V*S !convert from real to virtual velocity

!     SET INITIAL THERMOSTAT VALUES
      INQUIRE(FILE="initthermostat1.dat", EXIST=bool1 ) 
      if ( bool1 ) then
        write(*,*) "reading init file: initthermostat1.dat"
        OPEN (UNIT=17, FILE='initthermostat1.dat', STATUS='UNKNOWN')
        read(17,*) S, PS, HNOT
        CLOSE(17)
        loadHNOT = .true.
      ELSE
        !S = 1.D0
        !PS = 0.D0
        write(*,*) "setting thermostat initial values S=",S, "PS=",PS
        loadHNOT = .false.
      endif


!     START SIMULATION

      ISTEP  = 0  
      ICNT   = 0
      SUME   = 0.D0
      SUMKE  = 0.D0
      SUMPE  = 0.D0
      SUMTMT = 0.D0
      SUMTMR = 0.D0
      SUMTMP = 0.D0
      SUMPRS = 0.D0
      
      PHI(:,:) = 0.D0
      KE       = 0.D0


      UPDATE = .TRUE.

!     CALCULATE FORCES AND TORQUES AT TIME t=0
      CALL FORQUE(PE, VRL)

      DO I = 1, NMOL

!        TORQUE AT TIME t=0 IN THE BODY FRAME USING THE TRANSPOSE OF THE ROTATION MATRIX
         Q    = QTRN(I,:)
         CALL CONSTRUCT_RMXT( RMXT, Q )
         T4(I,1)   = 0.D0
         T4(I,2:4) = MATMUL(RMXT,T(I,:))
         ! intitialize P, the angular momentum.  P is virtual, W is real
         !CALL CONSTRUCT_MX( MX, Q )
         !W4(1)  = 0.D0
         !W4(2)  = 2.D0*W(I,1)*JMI(1)
         !W4(3)  = 2.D0*W(I,2)*JMI(2)
         !W4(4)  = 2.D0*W(I,3)*JMI(3)
         !P(I,:) = MATMUL(MX,W4*S)
         call PfromW( P(I,:), W(I,:), JMI, Q, S )
         !if ( I .eq. 1 ) then
           !dummyvec4=P(I,:)
           !dummyvec3=W(I,:)
           !write(*,*) dummyvec4
           !write(*,*) dummyvec3
           !call WfromP( dummyvec4, dummyvec3, JMI, Q, S )
           !write(*,*)  
           !write(*,*) dummyvec4
           !write(*,*) dummyvec3
         !endif
 
         !calculate the energy
         KE = KE + 0.5D0*MASS*(V(I,1)*V(I,1) + V(I,2)*V(I,2) + V(I,3)*V(I,3))/S**2 &
                 + 0.5D0*(W(I,1)*W(I,1)*JMI(1) + W(I,2)*W(I,2)*JMI(2) + W(I,3)*W(I,3)*JMI(3)) ! here V is virtual and W is real

      ENDDO

      IF ( loadHNOT .eqv. .false. ) THEN
        HNOT = KE + PE + 0.5D0*PS*PS/MQ + DBLE(G)*TMPFIX*LOG(S)
      ENDIF

      !print initial values for bug testing only
      IF ( .true. ) then
        OPEN (UNIT=71, FILE='pos0.dat', STATUS='UNKNOWN')
        OPEN (UNIT=72, FILE='ortn0.dat', STATUS='UNKNOWN')
        OPEN (UNIT=73, FILE='lvel0.dat', STATUS='UNKNOWN')
        OPEN (UNIT=74, FILE='avel0.dat', STATUS='UNKNOWN')
        WRITE(*,*) "Writing to initial output files"
        do I=1,NMOL
          WRITE (71, *)  R(I,1), R(I,2), R(I,3)
          WRITE (72, 900)  QTRN(I,1), QTRN(I,2), QTRN(I,3), QTRN(I,4) 
          WRITE (73, *)  V(I,1)/S, V(I,2)/S, V(I,3)/S
          WRITE (74, *) W(I,1), W(I,2), W(I,3)
        enddo
        close(71)
        close(72)
        close(73)
        close(74)
        OPEN (UNIT = 75, FILE = 'thermostat0.dat', STATUS = 'UNKNOWN')
          WRITE (75, *) S, PS, HNOT
        CLOSE(UNIT=75)
           !OPEN (UNIT=29, FILE='angmom0.dat', STATUS='UNKNOWN', POSITION='APPEND')
           !DO I = 1, NMOL
             !WRITE (29, 900)  P(I,1), P(I,2), P(I,3), P(I,4) 
           !ENDDO
           !CLOSE (UNIT=29,  STATUS='KEEP')  
      endif


      DO WHILE (ISTEP < NSTEP)            

         ISTEP = ISTEP + 1
         ICNT  = ICNT + 1
         IF(MOD(ISTEP,1000).EQ.0) WRITE(*,*) "Step number:", ISTEP, "of ", NSTEP
!     ADVANCE POSITIONS, ORIENTATIONS, AND THEIR TIME DERIVATIVES
!     MOVE also calls CHECK() and FORQUE
        
         CALL NVTMOVE(ISTEP, NEQ, PE, VRL, TKE, RKE)

!     CALCULATE TEMPERATURES, PRESSURE AND TOTEL ENERGY AT TIME t

         TMPTR   = 2.D0 * TKE  / DBLE(3 * NMOL - 3)
         TMPRT   = 2.D0 * RKE / DBLE(3 * NMOL)           
         KE      = TKE + RKE
         !DH is used only in NVT
         DH      = ((KE + PE + 0.5D0*PS*PS/MQ + DBLE(G)*TMPFIX*LOG(S) - HNOT)/HNOT) 
         TMP     = 2.D0 * KE / DBLE(G)
         PEPP    = PE / DBLE(NMOL)
         KEPP    = KE / DBLE(NMOL)
         EPP     = PEPP + KEPP
         PRS     = (2.D0 * TKE + VRL) / (3.D0 * VLM)
         SUME    = SUME + EPP
         SUMKE   = SUMKE + KEPP
         SUMPE   = SUMPE + PEPP
         SUMTMT  = SUMTMT + TMPTR
         SUMTMR  = SUMTMR + TMPRT
         SUMTMP  = SUMTMP + TMP
         SUMPRS  = SUMPRS + PRS

         !do block averaging of potential energy, temperature, and pressure
         IF (ISTEP > NEQ) THEN

            BSUMPE  = BSUMPE + PEPP
            BSUMTMP = BSUMTMP + TMP
            BSUMPRS = BSUMPRS + PRS

            IF (MOD((ISTEP-NEQ),IDUMP4) == 0) THEN

               AVPEB  = BSUMPE / DBLE(IDUMP4)
               AVTMPB = BSUMTMP / DBLE(IDUMP4)
               AVPRSB = BSUMPRS / DBLE(IDUMP4)  

               OPEN (UNIT=25, FILE='blockav1.dat', STATUS='UNKNOWN', POSITION='APPEND')
                 WRITE (25, *) AVPEB, AVTMPB, AVPRSB
               CLOSE (UNIT=25, STATUS='KEEP')

               BSUMPE  = 0.D0
               BSUMTMP = 0.D0
               BSUMPRS = 0.D0

            ENDIF

         ENDIF
 
         !write out energies, temperatures, and their averages, etc.
         IF (MOD(ISTEP, IDUMP1) .EQ. 0) THEN

            AVE    = SUME / DBLE(ICNT)
            AVKE   = SUMKE / DBLE(ICNT)
            AVPE   = SUMPE / DBLE(ICNT)
            AVTMPT = SUMTMT / DBLE(ICNT)
            AVTMPR = SUMTMR / DBLE(ICNT)
            AVTMP  = SUMTMP / DBLE(ICNT)
            AVPRS  = SUMPRS / DBLE(ICNT)

            OPEN (UNIT=3, FILE='enrg1.dat', STATUS='UNKNOWN', POSITION='APPEND')
                 WRITE (3, *) KEPP, PEPP, EPP 
            CLOSE (UNIT=3, STATUS='KEEP')
            OPEN (UNIT=4, FILE='tmpprs1.dat', STATUS = 'UNKNOWN', POSITION='APPEND')
                 WRITE (4, *) TMP, PRS
            CLOSE (UNIT=4, STATUS='KEEP') 

            OPEN (UNIT=33, FILE='avenrg1.dat', STATUS='UNKNOWN', POSITION='APPEND')
                 WRITE (33, *) AVKE, AVPE, AVE
            CLOSE (UNIT=33, STATUS='KEEP')
            OPEN (UNIT=34, FILE='avtmpprs1.dat', STATUS = 'UNKNOWN', POSITION='APPEND')
                 WRITE (34, *) AVTMP, AVPRS
            CLOSE (UNIT=34, STATUS='KEEP') 

            OPEN (UNIT=35, FILE='avtmptr1.dat', STATUS = 'UNKNOWN', POSITION='APPEND')
                 WRITE (35, *) AVTMPT, AVTMPR
            CLOSE (UNIT=35, STATUS='KEEP')

            OPEN (UNIT=36, FILE='error1.dat', STATUS = 'UNKNOWN', POSITION='APPEND')
                 WRITE (36, *) DH,  DH*HNOT+HNOT
            CLOSE (UNIT=36, STATUS='KEEP')

            OPEN (UNIT=37, FILE='thermostat1.dat', STATUS = 'UNKNOWN', POSITION='APPEND')
                 WRITE (37, *) S, PS
            CLOSE (UNIT=37, STATUS='KEEP')

         ENDIF

         ! calculate energy averaged over time and molecules.
         ! calculate variance over molecules of time averaged energy
         IF (ISTEP > NEQ) THEN

            CALL EFM (EFMT, TIME)

         ENDIF

         !write out positions, orientations, and angular change since t=NEQ
         IF (ISTEP > NEQ .AND. MOD(ISTEP, IDUMP2) == 0) THEN   

           OPEN (UNIT=7,  FILE='pos1.dat', STATUS='UNKNOWN', POSITION='APPEND')
           OPEN (UNIT=8,  FILE='ortn1.dat', STATUS='UNKNOWN', POSITION='APPEND')
           OPEN (UNIT=28, FILE='rot1.dat', STATUS='UNKNOWN', POSITION='APPEND')
           OPEN (UNIT=29, FILE='efm1.dat', STATUS='UNKNOWN', POSITION='APPEND')

           DO I = 1, NMOL
             
             WRITE (7, *)  R(I,1), R(I,2), R(I,3)
             WRITE (8, 900)  QTRN(I,1), QTRN(I,2), QTRN(I,3), QTRN(I,4) 
             WRITE (28, *) PHI(I,1), PHI(I,2), PHI(I,3)  
             900 FORMAT (1x,4ES25.16)

           ENDDO

           WRITE(29, *) TIME, EFMT


           CLOSE (UNIT=7,  STATUS='KEEP')  
           CLOSE (UNIT=8,  STATUS='KEEP')  
           CLOSE (UNIT=28, STATUS='KEEP')
           CLOSE (UNIT=29, STATUS='KEEP')

           !OPEN (UNIT=29, FILE='angmom1.dat', STATUS='UNKNOWN', POSITION='APPEND')
           !DO I = 1, NMOL
             !WRITE (29, 900)  P(I,1), P(I,2), P(I,3), P(I,4) 
           !ENDDO
           !CLOSE (UNIT=29,  STATUS='KEEP')  

           IF (MOD(ISTEP, IDUMP4) .EQ. 0) THEN   

             OPEN (UNIT=9, FILE='lvel1.dat', STATUS='UNKNOWN', POSITION='APPEND')
             OPEN (UNIT=10, FILE='avel1.dat', STATUS='UNKNOWN', POSITION='APPEND')

             DO I = 1, NMOL

               WRITE (9, *)  V(I,1)/S, V(I,2)/S, V(I,3)/S
               WRITE (10, *) W(I,1), W(I,2), W(I,3)

             ENDDO

             CLOSE (UNIT = 9, STATUS = 'KEEP')  
             CLOSE (UNIT = 10, STATUS ='KEEP')  

           ENDIF

         ENDIF
               
         !reinitialize energy averages
         IF ( ISTEP == NEQ) THEN
              
            ICNT   = 0
            SUME   = 0.D0
            SUMKE  = 0.D0
            SUMPE  = 0.D0
            SUMTMT = 0.D0
            SUMTMR = 0.D0
            SUMTMP = 0.D0
            SUMPRS = 0.D0        
            TIME   = 0.D0

            BSUMPE  = 0.D0
            BSUMTMP = 0.D0
            BSUMPRS = 0.D0
                         
         ENDIF  

      ENDDO

      !program ended, print out the final data

      OPEN (UNIT = 21, FILE = 'finalpos1.dat', STATUS = 'UNKNOWN')
      OPEN (UNIT = 22, FILE = 'finalortn1.dat', STATUS = 'UNKNOWN')
      OPEN (UNIT = 23, FILE = 'finallinvel1.dat', STATUS = 'UNKNOWN')
      OPEN (UNIT = 24, FILE = 'finalangvel1.dat', STATUS = 'UNKNOWN')
      
      DO I = 1, NMOL

         WRITE (21, *) R(I,1), R(I,2), R(I,3)
         WRITE (22, 900) QTRN(I,1), QTRN(I,2), QTRN(I,3), QTRN(I,4)
         WRITE (23, *) V(I,1)/S, V(I,2)/S, V(I,3)/S
         WRITE (24, *) W(I,1), W(I,2), W(I,3)

      ENDDO

      CLOSE(UNIT=21)
      CLOSE(UNIT=22)
      CLOSE(UNIT=23)
      CLOSE(UNIT=24)

      OPEN (UNIT = 25, FILE = 'finalthermostat1.dat', STATUS = 'UNKNOWN')
         WRITE (25, *) S, PS, HNOT
      CLOSE(UNIT=25)
           !OPEN (UNIT=29, FILE='finalangmom1.dat', STATUS='UNKNOWN', POSITION='APPEND')
           !DO I = 1, NMOL
             !WRITE (29, 900)  P(I,1), P(I,2), P(I,3), P(I,4) 
           !ENDDO
           !CLOSE (UNIT=29,  STATUS='KEEP')  

!      DEALLOCATE (R, QTRN, V, W, F, T, P, T4, PHI, FSITE, MST, KEP, PEP, SUMEPT, STAT = STATUS)

      call RUN_INFO()

!      END PROGRAM LWOTP_MD_NVT
      END SUBROUTINE RIGIDMD_NVT
!     ------------------------------------------------------------------------------------
 
      SUBROUTINE NVTSTEP1 (S, PS, FCTMQT )

      !STEP 1 in the MOVE algorithm
      ! STEP 1: update S, PS
      ! S = S*(1.D0 + PS*0.5D0*HALFDT/MQ)**2
      ! PS = PS/(1.D0 + PS*0.5D0*HALFDT/MQ)

      IMPLICIT NONE

      DOUBLE PRECISION, INTENT(INOUT) :: S, PS
      DOUBLE PRECISION, INTENT(IN) :: FCTMQT
      DOUBLE PRECISION :: FCTR

      FCTR = 1.D0 + FCTMQT*PS
      S    = S*FCTR*FCTR
      PS   = PS/FCTR    

      END SUBROUTINE NVTSTEP1

!     ------------------------------------------------------------------------------------
 
      SUBROUTINE NVTSTEP2 (PE)

      !STEP 2 in the MOVE algorithm
      ! STEP 2: update PS, V, and P from the Torque
      ! PS = PS - PE*HALFDT   !! PE=potential energy
      ! loop over molecules
      !   update linear momentum (or velocity)
      !   V(I,:) = V(I,:) + S * F(I,:) * HALFDT / MASS   
      !   update rotational momentum 
      !   P(I,:) = P(I,:) +  2*S * DELT/2 * MATMUL(MX,T4(I,:))
      ! end loop over molecules
      
      USE MDCOMMONS, only : PS, HALFDT, DELT, NMOL, V, QTRN, S, F, MASS, P, T4

      IMPLICIT NONE

      INTEGER :: I
      DOUBLE PRECISION :: MX(4,4), Q(4)
      DOUBLE PRECISION, INTENT(IN) :: PE

      PS     = PS - PE * HALFDT

      DO I = 1, NMOL

         V(I,:) = V(I,:) + S * F(I,:) * HALFDT / MASS          

         Q      = QTRN(I,:)

         CALL CONSTRUCT_MX( MX, Q )

         P(I,:)     = P(I,:) +  S * DELT * MATMUL(MX,T4(I,:))

      ENDDO

      END SUBROUTINE NVTSTEP2

!     ------------------------------------------------------------------------------------
 
      SUBROUTINE NVTSTEP3 ()

      !STEP 3 in the MOVE algorithm
      ! STEP 3: update QTRN, P, PS with contributions from principle axis 3
      ! ZETASUM = 0
      ! loop over molecules
      !   ZETA(I,3) = dot_product( P(I,:) , MX(:,4) ) / ( 4 * S * JMI(3) )
      !   ZETASUM = ZETASUM + ZETA(I,3)**2
      !   Q(I) = cos( ZETA(I,3)*HALFDT )*Q(I) + sin( ZETA(I,3)*HALFDT )*MX(:,4)
      !   P(I) = cos( ZETA(I,3)*HALFDT )*P(I) + sin( ZETA(I,3)*HALFDT )*MP(:,4)
      !   !   in the above, MP(:,:) is the analogous matrix to MX for vector P(I,:),
      !   !   i.e. MP(:,4) = (/-P(I,4),P(I,3),-P(I,2),P(I,1)/)
      ! end loop over molecules
      ! PS = PS + ZETASUM*JMI(3)*DELT
      ! contribution from principle axis 2
      
      USE MDCOMMONS, only : NMOL, QTRN, P, S, PS, JMI, HALFDT, DELT

      IMPLICIT NONE

      INTEGER :: I
      DOUBLE PRECISION :: SMPHSQ, Q(4), DQ(4), PHIT
      DOUBLE PRECISION :: PI(4), CP, SP

      SMPHSQ = 0.D0
      DO I = 1, NMOL

         Q         = QTRN(I,:)
         PI        = P(I,:)
         !DQ        = MX(:,4)
         DQ        = (/-Q(4),Q(3),-Q(2),Q(1)/)
         PHIT      = 0.25D0*DOT_PRODUCT(PI,DQ)/(S*JMI(3))
         SMPHSQ    = SMPHSQ + PHIT * PHIT
         CP        = COS(PHIT*HALFDT)
         SP        = SIN(PHIT*HALFDT)
         P(I,:)    = CP*PI + SP*(/-PI(4),PI(3),-PI(2),PI(1)/)
         QTRN(I,:) = CP*Q + SP*DQ

      ENDDO

      PS     = PS + JMI(3)*SMPHSQ*DELT

      END SUBROUTINE NVTSTEP3

!     ------------------------------------------------------------------------------------
 
      SUBROUTINE NVTSTEP4 ()

      !STEP 4 in the MOVE algorithm
      ! STEP 4: update QTRN, P, PS with contributions from principle axis 2
      ! ZETASUM = 0
      ! loop over molecules
      !   ZETA(I,2) = dot_product( P(I,:) , MX(:,4) ) / ( 4 * S * JMI(2) )
      !   ZETASUM = ZETASUM + ZETA(I,2)**2
      !   Q(I) = cos( ZETA(I,2)*HALFDT )*Q(I) + sin( ZETA(I,2)*HALFDT )*MX(:,2)
      !   P(I) = cos( ZETA(I,2)*HALFDT )*P(I) + sin( ZETA(I,2)*HALFDT )*MP(:,2)
      !   !   in the above, MP(:,:) is the analogous matrix to MX for vector P(I,:),
      !   !   i.e. MP(:,2) = (/-P(I,3),-P(I,4),P(I,1),P(I,2)/)
      ! end loop over molecules
      ! PS = PS + ZETASUM*JMI(2)*DELT
      
      USE MDCOMMONS, only : NMOL, QTRN, P, S, PS, JMI, HALFDT, DELT

      IMPLICIT NONE

      INTEGER :: I
      DOUBLE PRECISION :: SMPHSQ, Q(4), DQ(4), PHIT
      DOUBLE PRECISION :: PI(4), CP, SP

      SMPHSQ = 0.D0
      DO I = 1, NMOL

         Q         = QTRN(I,:)
         PI        = P(I,:)
         DQ        = (/-Q(3),-Q(4),Q(1),Q(2)/)
         PHIT      = 0.25D0*DOT_PRODUCT(PI,DQ)/(S*JMI(2))
         SMPHSQ    = SMPHSQ + PHIT*PHIT
         CP        = COS(PHIT*HALFDT)
         SP        = SIN(PHIT*HALFDT)
         P(I,:)    = CP*PI + SP*(/-PI(3),-PI(4),PI(1),PI(2)/)
         QTRN(I,:) = CP*Q + SP*DQ

      ENDDO

      PS     = PS + JMI(2)*SMPHSQ*DELT

      END SUBROUTINE NVTSTEP4

!     ------------------------------------------------------------------------------------
 
      SUBROUTINE NVTSTEP5 ()

      !STEP 5 in the MOVE algorithm
      ! STEP 5: update R QTRN, P, PS with contributions from principle axis 1
      ! STEP 5 is analogous to STEP 3 and STEP 4 but for ZETA(:,1), i.e. the
      ! first principle axis,
      ! EXCEPT, the update of PS is different and involves the translational
      ! kinetic energy, the temperature (TMPFIX), the number of degrees of
      ! freedom (g), and HNOT.  Also, R is updated
      
      USE MDCOMMONS, only : NMOL, QTRN, P, S, PS, JMI, HALFDT, DELT, R, V, G, &
        TMPFIX, HNOT, MASS

      IMPLICIT NONE

      INTEGER :: I
      DOUBLE PRECISION :: SMPHSQ, Q(4), DQ(4), PHIT, DBLE
      DOUBLE PRECISION :: PI(4), CP, SP
      DOUBLE PRECISION :: TKE !translational kinetic energy

      SMPHSQ = 0.D0
      TKE = 0
      DO I = 1, NMOL

         R(I,:)    = R(I,:) + V(I,:) * DELT / S
         Q         = QTRN(I,:)
         PI        = P(I,:)
         DQ        = (/-Q(2),Q(1),Q(4),-Q(3)/)
         PHIT      = 0.25D0*DOT_PRODUCT(PI,DQ)/(S*JMI(1))
         SMPHSQ    = SMPHSQ + PHIT * PHIT
         CP        = COS(PHIT*DELT)
         SP        = SIN(PHIT*DELT)
         P(I,:)    = CP*PI + SP*(/-PI(2),PI(1),PI(4),-PI(3)/)
         QTRN(I,:) = CP*Q  + SP*DQ
         TKE       = TKE + V(I,1)*V(I,1) + V(I,2)*V(I,2) + V(I,3)*V(I,3)

      ENDDO

      PS = PS + (0.5D0*MASS*TKE/(S*S) + 2.D0*JMI(1)*SMPHSQ &
           - DBLE(G)*TMPFIX*(1.D0+LOG(S)) + HNOT ) * DELT

      END SUBROUTINE NVTSTEP5
