!   PATHSAMPLE: A driver for OPTIM to create stationary point databases using 
!   discrete path sampling and perform kinetic analysis
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

C  New version of PATHSAMPLE using the committor probability formulation
C  and a new sampling strategy for stationary points. Paths themselves
C  are not generally referenced directly.
C  Initially we need to know
C  which minima are in which funnel, i.e. members of the sets A and B for
C  which the rate constant is required, and there may be pre-existing databases
C  of minima and transition states containing energies and products of non-zero
C  eigenvalues. We need the ability to start from a single
C  path.info file.
C
C  We are allowing connections between any A and B permutational isomers.
C  Use TAG to distinguish permutational isomers.
C
C SAT: BOUNDARYT variable is referenced but never set! Assume default value of .False. and ignore relevant bits of code.
C SAT: SAVEALL variable is set to .True. in setup if 'paths' file exists and 'pathpointer' does not!
C
      PROGRAM PATHSAMPLE
      USE COMMONS
      USE NODES, ONLY: SSHPARALLEL, JPN
      USE UTILS,ONLY : GETUNIT
      USE DOCKMODULE
      USE RIGIDBODYMOD, ONLY: INITIALISERIGiDBODY, CLEANRIGIDBODIES
      IMPLICIT NONE
      INTEGER J1, NSIZE, NWORST, NAVAIL, NDUMMY, LUNIT, NMINSAVE
      INTEGER, ALLOCATABLE :: FREEMINLIST(:), FREEMINPOINT(:)
      CHARACTER(LEN=200) :: VERSIONTEMP
      INTEGER, ALLOCATABLE :: MINMAP(:)
!     INTEGER, ALLOCATABLE :: NSEED(:)
      LOGICAL MASSFILE, YESNOT
      DOUBLE PRECISION TINIT, TNEW, XDUMMY

!     AMH LOCAL VARIABLES
      INTEGER :: NRES,I_RES,J2
      CHARACTER(LEN=5) :: TARFL
      CHARACTER(LEN=2) :: SDUMMY
      INTEGER :: SEQ(500)

      VERSIONTEMP=SVNVERSION
      VERSIONTEMP=TRIM(ADJUSTL(VERSIONTEMP))
C
C  CONNECTIONS is the minimum number of connected minima.
C  TEMPERATURE is the temperature in reduced units.
C  Use EMICRO for a microcanonical total energy. Don't try to access inaccessible
C  stationary points above EMICRO!
C
      CALL CPU_TIME(TINIT)
!     CALL RANDOM_SEED                    ! Initializing
!     CALL RANDOM_SEED(SIZE=NSIZE)        ! Sets NSIZE
!     ALLOCATE(NSEED(NSIZE))
!     NSEED(1:NSIZE)=ISEED                ! specfied by SEED keyword, default 1.
!     CALL RANDOM_SEED(PUT=NSEED(1:NSIZE))

      CALL SDPRND(ISEED)

      PRINTT=.TRUE.
      TTSSEARCH=0.0D0
      TPFOLD=0.0D0
      TTFOLD=0.0D0
      TDIJKSTRA=0.0D0
      TKSHORTESTPATHS=0.0D0
      TCONNECTDIST=0.0D0
      TGT=0.0D0
      PRINT '(A,A70)', 'PATHSAMPLE version r',VERSIONTEMP
      CALL KEYWORDS
      IF (SLURMT .AND. (PBST.OR.SSHT)) THEN
         PRINT '(A)','Keyword SLURM should not be used with keywords PBS or SSH. '
         STOP
      END IF
      IF ((PERMDIST.OR.PERMISOMER).AND.(CHARMMT).AND.(NPERMGROUP.EQ.NATOMS)) THEN
         PRINT '(A)','main> Likely error in input - PERM set without perm.allow file for CHARMM potential'
         STOP
      ENDIF
      IF ((PERMDIST.OR.PERMISOMER).AND.(VERIFY('perm.allow',COPYFILES).NE.0)) THEN
!
! This is not an error for monoatomic systems.
!
         PRINT '(A)','main> WARNING - PERM set without perm.allow file in COPYFILES'
      ENDIF
      IF (CHARMMT.AND.(VERIFY('input.crd',COPYFILES).NE.0)) THEN
         PRINT '(A)','main> Likely error in input - CHARMM potential specified without input.crd in COPYFILES'
         STOP
      ENDIF
      ALLOCATE(FRQS(3*NATOMS),RESLABEL(NATOMS),
     &         ATOMLABEL(NATOMS),RESNUMBER(NATOMS),EMIN(MAXMIN),FVIBMIN(MAXMIN),PFMIN(MAXMIN), 
     &         IXMIN(MAXMIN),IYMIN(MAXMIN),IZMIN(MAXMIN),GPFOLD(MAXMIN),
     &         ETS(MAXTS),FVIBTS(MAXTS),KPLUS(MAXTS),KMINUS(MAXTS),NEGEIG(MAXTS),
     &         IXTS(MAXTS),IYTS(MAXTS),IZTS(MAXTS),HORDERMIN(MAXMIN),TOPPOINTER(MAXMIN),HORDERTS(MAXTS), 
     &         PLUS(MAXTS),MINUS(MAXTS),POINTERM(MAXTS),POINTERP(MAXTS))
      IF (DUMMYTST) THEN
         ALLOCATE(MINDISTMIN(MAXMIN))
         ALLOCATE(MINCURVE(MAXMIN))
         ALLOCATE(MINFRQ2(MAXMIN))
         MINDISTMIN(1:MAXMIN)=HUGE(1.0D0)
      ENDIF
      IF (DIJKSTRAT .OR. KSHORTESTPATHST) ALLOCATE(TSATTEMPT(MAXTS))
      IF (NEWCONNECTIONST) ALLOCATE(MINCONN(MAXMIN))
      IF (DIJKSTRAT .OR. KSHORTESTPATHST) ALLOCATE(DMIN1(MAXMIN),DMIN2(MAXMIN))
      IF (DIJINITT) THEN
         IF (INITIALDIST) THEN
            ALLOCATE(ALLPAIRS((MAXMIN*(MAXMIN-1))/2))
            ALLPAIRS(1:(MAXMIN*(MAXMIN-1))/2)=-DISBOUND ! fixed initial default value.
         ELSE
            ALLOCATE(PAIRDIST(MAXMIN,PAIRDISTMAX),PAIRLIST(MAXMIN,PAIRDISTMAX))
            PAIRDIST(1:MAXMIN,1:PAIRDISTMAX)=1.0D100
            PAIRLIST(1:MAXMIN,1:PAIRDISTMAX)=-1
         ENDIF
      ENDIF
      IF (CONNECTREGIONT) ALLOCATE(DMIN1(DMINMAX),DMIN2(DMINMAX))
      IF (UNTRAPT .AND. (.NOT. DIJKSTRAT) .AND. (.NOT. KSHORTESTPATHST)) ALLOCATE(DMIN1(DMINMAX),DMIN2(DMINMAX))

      IF (NATOMS.LT.1) THEN
         PRINT*,'ERROR - NATOMS=',NATOMS
         STOP
      ENDIF

      IF (ANGLEAXIS) THEN
         CALL INITIALISERIGIDBODY(ZSYM//'   ', NATOMS/2, MASS)
         ALLOCATE(ZSYMBOL(NATOMS))
         ZSYMBOL=ZSYM
      ELSE
         ALLOCATE(MASS(NATOMS),ZSYMBOL(NATOMS))
         DO J1=1,NATOMS
            MASS(J1)=1.0D0
            ZSYMBOL(J1)=ZSYM
         ENDDO
         IF (NIMET .OR. NIHEAM7T .OR. NIH2LEPST) ZSYMBOL(1)='H'
         IF (NIH2LEPST) ZSYMBOL(2)='H'
      ENDIF
      IF (NTAG.GT.0) THEN
         DO J1=1,NATOMS
            MASS(J1)=MASS(J1)*TAGFAC(J1)
         ENDDO
      ENDIF
C
C  If there is a file called mass in the current directory then use these masses instead
C  of the defaults.
C
      INQUIRE(FILE='mass',EXIST=MASSFILE)
      IF (MASSFILE) THEN
         WRITE(*,'(A)') 'main> Reading individual atomic symbols and masses from file masses'
         OPEN(UNIT=1,FILE='mass',STATUS='OLD')
         READ (1,*) (ZSYMBOL(J1),MASS(J1),J1=1,NATOMS)
         CLOSE(1)
C        WRITE(*,'(A2,F20.10)') (ZSYMBOL(J1),MASS(J1),J1=1,NATOMS)
      ENDIF
C
C  NFSTART is the first true frequency to use from a path.info file, and
C  NFFINISH is the last non-zero frequency to read for a minimum.
C  Used in getnewpath and tssearch. KAPPA is the number of non-zero 
C  vibrational frequencies.
C
      NFSTART=1
      NFFINISH=3*NATOMS-6
      KAPPA=3*NATOMS-6
      NGLY=0
      IF (NZEROS.GT.0) THEN
         IF (MLP3T.OR.MLPB3T) THEN
            NFSTART=1
            NFFINISH=NATOMS-NZEROS
!
! This choice for kappa just gives the additive half kT per non-zero vibrational
! frequency for the configurational part of the heat capacity. This will be the 
! high T limit.
!
            KAPPA=NATOMS-1
         ELSE
            NFSTART=1
            NFFINISH=3*NATOMS-NZEROS
            KAPPA=3*NATOMS-NZEROS
         ENDIF
         PRINT '(A,2I6)','main> start and finish for non-zero frequencies for minima set to ',NFSTART,NFFINISH
      ELSEIF (TWOD) THEN
         NFSTART=1
         NFFINISH=2*NATOMS-3
         KAPPA=2*NATOMS-3
         IF (BULKT) THEN
            NFSTART=1
            NFFINISH=2*NATOMS-2
            KAPPA=2*NATOMS-2
         ENDIF
      ELSE IF (GTHOMSONT) THEN
        NFSTART = 1
        NFFINISH = 3*NATOMS - 3
        KAPPA = 3*NATOMS - 3
      ELSE IF (MLP3T.OR.MLPB3T) THEN
         NFSTART=1
         NFFINISH=NATOMS
!
! This choice for kappa just gives the additive half kT per non-zero vibrational
! frequency for the configurational part of the heat capacity. This will be the 
! high T limit.
!
         KAPPA=NATOMS-1
      ELSE IF (MLLJAT3) THEN
         NFSTART=1
         NFFINISH=NATOMS
!
! This choice for kappa just gives the additive half kT per non-zero vibrational
! frequency for the configurational part of the heat capacity. This will be the 
! high T limit. 
!
         KAPPA=NATOMS
      ELSE IF (PHI4MODT) THEN
         NFSTART=1
         NFFINISH=NATOMS
         KAPPA=NATOMS
      ELSE IF (PULLT) THEN
         NFFINISH=3*NATOMS-4
         KAPPA=3*NATOMS-4
      ELSE IF (BULKT) THEN
         NFFINISH=3*NATOMS-3
         KAPPA=3*NATOMS-3
      ELSE IF (TRAPT) THEN
         NFFINISH=3*NATOMS-3
         KAPPA=3*NATOMS-3
      ELSE IF (MKTRAPT) THEN
         NFFINISH=3*NATOMS
         KAPPA=3*NATOMS
      ELSE IF (FREEZE) THEN
         NFFINISH=3*NATOMS-3*NFREEZE
         KAPPA=3*NATOMS-3*NFREEZE
      ELSE IF (RIGIDINIT) THEN ! hk286
         NFFINISH=DEGFREEDOMS - 6
         KAPPA=DEGFREEDOMS - 6
      ELSE IF (UNRST) THEN
         NFFINISH=NINTS
         KAPPA=NINTS
      ELSE IF (AMHT) THEN
        OPEN(UNIT=30,FILE='pro.list',STATUS='OLD',FORM='FORMATTED')
         READ (30,1000)TARFL
1000          FORMAT(A5)
        CLOSE(30)
        OPEN(30,FILE='proteins/'//TARFL,STATUS='OLD')
         READ(30,*)
         READ(30,*)NRES
c        WRITE(*,*)NRES
        IF (NRES.GT.500) THEN
            WRITE(6,*) 'FAILURE NRES GR THAN 500 GEOPT'
            STOP
        ENDIF
        READ (30,25)(SEQ(I_RES),I_RES=1,NRES)
c         WRITE (*,25)(SEQ(I_RES),I_RES=1,NRES)
25           FORMAT(25(I2,1X))
             CLOSE(30)
            DO J2=1,NRES
              IF (SEQ(J2).EQ.8) THEN
                NGLY = NGLY +1
              ENDIF
            ENDDO 
   
c            WRITE(*,'(A,I5)') 'AMH RES COUNT', NRES  
c            WRITE(*,'(A,I5)') 'AMH GLY COUNT', NGLY            

         NFFINISH=3*(NATOMS-NGLY)-6
         KAPPA=3*(NATOMS-NGLY)-6
      ENDIF

      WRITE(*,'(A)') '*************************************************************************************************'
      IF (.NOT.DEBUG) WRITE(*,'(A)') 'debug printing is OFF'
      IF (DEBUG) WRITE(*,'(A)') 'debug printing is ON'
      IF (PBST) THEN
         WRITE(*,'(A,I10,A)') 'Running on ',NCPU,' processors for distributed memory PBS setup'
      ELSE IF (SLURMT) THEN
         WRITE(*,'(A,I10,A)') 'Running on ',NCPU,' processors for distributed memory SLURM setup'
      ELSE
         WRITE(*,'(A,I10,A)') 'Running on ',NCPU,' processors for shared memory, single node setup'
      ENDIF
      WRITE(*,'(A,I10)') 'Number of jobs per processor=',JPN
      WRITE(*,'(A,2G20.10)') 'Sleep times between OPTIM job submissions (s)=',SLEEPTIME1,SLEEPTIME2
      IF (TRIM(ADJUSTL(COPYFILES)).NE.'') 
     &  PRINT '(2A)','The following additional files will be copied to distributed nodes: ',TRIM(ADJUSTL(COPYFILES))
      IF (ADDPT) WRITE(*,'(A)') 'all permutational isomers with different moments of inertia will be enumerated automatically'
      IF (NTAG.GT.0) THEN
         PRINT '(I10,A)',NTAG,' tagged atoms with mass factors:'
         PRINT '(I10,F12.2)',(TAGNUM(J1),TAGFAC(TAGNUM(J1)),J1=1,NTAG)
      ENDIF

      IF (TWOD) WRITE(*,'(A)') 'two-dimensional system'
      IF (BULKT) THEN
         WRITE(*,'(A)') 'bulk system'
         WRITE(*,'(A,3G20.10)') 'x, y, z box lengths: ',BOXLX, BOXLY, BOXLZ
      ENDIF
      IF (.NOT.BULKT) WRITE(*,'(A)') 'finite system (not bulk)'
      WRITE(*,'(A,A2)') 'OPTIM system symbol: ',ZSYM
      WRITE(*,'(A,I5)') 'number of atoms=',NATOMS
      IF (FREEZE) WRITE(*,'(A,I5)') 'number of frozen atoms=',NFREEZE
      IF (AMHT) WRITE(*,'(A,I10)') 'AMH potential, number of glycines=',NGLY
      IF (NTAG.GT.0) WRITE(*,'(A,I5,A,F12.3)') 'tagged atom ',NTAG,' has mass ',TAGMASS
      WRITE(*,'(A,I10)') 'minimum number of connections for minima=',CONNECTIONS
      WRITE(*,'(A,I10)') 'maximum number of ts searches to achieve this minimum value=',MAXTSATTEMPTS
      WRITE(*,'(A,I10)') 'random number seed=',ISEED
      CALL RANDOM_SEED(ISEED)
      WRITE(*,'(A,F15.5)') 'initial perturbation parameter single-ended transition state searches=',PERTVALUE
      WRITE(*,'(A,E20.10)') 'energy difference  criterion for distinguishing stationary points=',EDIFFTOL
      WRITE(*,'(A,E20.10)') 'distance criterion for distinguishing stationary points=',GEOMDIFFTOL
      WRITE(*,'(A,E20.10)') 'inertia difference criterion for setup check=',IDIFFTOL
      IF (TFOLDT) THEN
         PRINT '(A)','Calculation of MFPT and rate constants using Gauss-Seidel/successive overrelaxation iteration'
         PRINT '(A,F15.1)','Maximum steps=',NTFOLD
         PRINT '(A,G20.10)','Convergence condition=',TFOLDTHRESH
         PRINT '(A,G20.10)','Relaxation factor (omega=1 is Gauss-Seidel)=',TOMEGA
      ENDIF
      IF ((.NOT.GTT).AND.(.NOT.GT2T).AND.(.NOT.NGTT)) THEN
         IF (DIRECTION.EQ.'AB') WRITE(*,'(A)') 'sampling direction is a<-b'
         IF (DIRECTION.EQ.'BA') WRITE(*,'(A)') 'sampling direction is b<-a'
         WRITE(*,'(A,A80)') 'OPTIM executable: ',EXEC
      ENDIF
      IF (REMOVEUNCONNECTEDT) THEN
         WRITE(*,'(A,I5,A)') 'Rewriting database removing minima with ',NCONNMIN,' connections or fewer'
         WRITE(*,'(A)') 'and stationary points unconnected to ' // TRIM(ADJUSTL(UNCONNECTEDS))
      ENDIF
      IF (ENSEMBLE.EQ.'T') WRITE(*,'(A,F15.5)') 'Temperature=',TEMPERATURE
      IF (ENSEMBLE.EQ.'E') WRITE(*,'(A,F15.5)') 'Total energy=',TOTALE
      IF (NOFRQS) WRITE(*,'(A)') 'No frequency information should be present in info files'
      IF (KMCT) THEN
         WRITE(*,'(A,F20.1,A,F12.4)') 'Running KMC on database. Number of KMC runs=',NKMCCYCLES,
     1                             ' pair renorm threshold=',PAIRTHRESH
      ELSE IF (DIJKSTRAT.AND.(.NOT.DIJPAIRT)) THEN
         WRITE(*,'(A)') 'Performing Dijkstra analysis to find the most important paths'
         IF (KSHORTESTPATHST) THEN
             WRITE(*,'(A)') 'ERROR: KSHORTESTPATHS and DIJKSTRA are incompatible; stopping.'
             STOP
         ENDIF
      ELSE IF (KSHORTESTPATHST) THEN
         WRITE(*,'(A)') 'Performing k-th shortest path analysis to find the most important paths'
         IF (DIJKSTRAT) THEN
             WRITE(*,'(A)') 'ERROR: KSHORTESTPATHS and DIJKSTRA are incompatible; stopping.'
             STOP
         ENDIF
      ELSE IF (KMCCOMMITT) THEN
         WRITE(*,'(A)')        'Using committor probability/KMC routine:'
         WRITE(*,'(A,F20.1)')     '      KMC cycles per starting minimum=                 ',NKMCCYCLES
         WRITE(*,'(A,F20.1)')  '      Maximum number of steps in a KMC run=            ',MAXBREAK
         WRITE(*,'(A,G20.10)') '      Convergence condition on committor probabilities=',PABCONV
         WRITE(*,'(A,G20.10)') '      Successive overrelaxation parameter omega=       ',OMEGA
         WRITE(*,'(A,G20.10)') '      Threhold for leapfrog KMC moves=                 ',PAIRTHRESH
      ELSE IF (NGTT) THEN
          WRITE(*,'(A,I10,A,F12.4)') 'Running new graph theory rate calculation on database'
      ELSE IF (ARNOLDIT) THEN
          WRITE(*,'(A)') 'Eigenvalues of the transition matrix from weighted Arnoldi subspace method'
      ELSE IF (DIAGT) THEN
          WRITE(*,'(A,I10,A)') 'Explicit diagonalisation of the transition matrix to find ',NDIAGEIG,' eigenvalues'
          WRITE(*,'(A,G12.4)') 'All transition matrix elements will be scaled by ',DIAGSCALE
      ELSE IF (GTT.or.GT2T) THEN
          WRITE(*,'(A,I10,A,F12.4)') 'Running graph theory rate calculation on database'
      ELSE 
         WRITE(*,'(A,I8)') 'Number of cycles for connections between local minima=',NATTEMPT
         IF (PFOLDINT.NE.0) THEN
            PRINT '(A,I10)','Interval between calls to Pfold routine=',PFOLDINT 
            PRINT '(A,I10,A,G20.10)','SOR iterations per call=',NPFOLD,' with omega=',OMEGA 
         ELSE
            PRINT '(A,I10)','No Pfold iterations'
         ENDIF
      ENDIF
      IF (DIJKSTRAT.AND.(MAXDOWNBARRIER.LT.1.0D99)) THEN
         PRINT '(A,G20.10)','Downhill barrier check in Dijkstra, threshold=',MAXDOWNBARRIER
      ENDIF
      
      IF (FRICTIONT) WRITE(*,'(A,G20.10)') 'Phenomenological friction coefficient in frequency units=',GAMMAFRICTION
      IF (REGROUPT) WRITE(*,'(A,F15.5)') 'Reclassifying A/B/I minima for threshold energy ',REGROUPTHRESH
      IF (REGROUPRATET) THEN
         WRITE(*,'(A,G20.10)') 'Grouping and calculating free energies for threshold rate ',REGROUPRATETHRESH
         WRITE(*,'(A,G20.10)') 'Plank constant set to ',PLANCK
         IF (CHARMMT.AND.(PLANCK.EQ.1.0D0)) THEN
            PRINT '(2(A,G20.10))','ERROR - regrouping with CHARMM and Planck constant=',PLANCK,' resetting to',9.546D-14
            PLANCK=9.546D-14
         ENDIF
      ENDIF
      IF (REGROUPPET) THEN
         WRITE(*,'(A,G20.10)') 'Regrouping A/B/I minima for superbasin PE ',REGROUPPETHRESH
      ENDIF
      IF (RFMULTIT) THEN
         WRITE(*,'(A,G20.10)') 'Grouping and calculating free energies for threshold barrier ',REGROUPFREETHRESH
         WRITE(*,'(A,G20.10)') 'Plank constant set to ',PLANCK
         IF (CHARMMT.AND.(PLANCK.EQ.1.0D0)) THEN
            PRINT '(2(A,G20.10))','ERROR - regrouping with CHARMM and Planck constant=',PLANCK,' resetting to',9.546D-14
            PLANCK=9.546D-14
         ENDIF
      ENDIF
      IF (RFKMCT) THEN
         WRITE(*,'(A,G20.10)') 'KMC-type sampling of thermodynamics for free energy groups'
         WRITE(*,'(A,G20.10)') 'Plank constant set to ',PLANCK
         IF (CHARMMT.AND.(PLANCK.EQ.1.0D0)) THEN
            PRINT '(2(A,G20.10))','ERROR - regrouping with CHARMM and Planck constant=',PLANCK,' resetting to',9.546D-14
            PLANCK=9.546D-14
         ENDIF
      ENDIF
      IF (REGROUPFREET) THEN
         WRITE(*,'(A,G20.10)') 'Grouping and calculating free energies for threshold barrier ',REGROUPFREETHRESH
         WRITE(*,'(A,G20.10)') 'Plank constant set to ',PLANCK
         IF (CHARMMT.AND.(PLANCK.EQ.1.0D0)) THEN
            PRINT '(2(A,G20.10))','ERROR - regrouping with CHARMM and Planck constant=',PLANCK,' resetting to',9.546D-14
            PLANCK=9.546D-14
         ENDIF
      ENDIF
      IF (REGROUPFREEABT) THEN
         WRITE(*,'(A,G20.10)') 'Grouping and calculating free energies for threshold barrier ',REGROUPFREETHRESH
         WRITE(*,'(A,G20.10)') 'Free energy groups will not be used except to expand the A and B sets'
         WRITE(*,'(A,G20.10)') 'Plank constant set to ',PLANCK
         IF (CHARMMT.AND.(PLANCK.EQ.1.0D0)) THEN
            PRINT '(A,G20.10)','ERROR - regrouping with CHARMM and Planck constant=',PLANCK
            PRINT '(2(A,G20.10))','ERROR - regrouping with CHARMM and Planck constant=',PLANCK,' resetting to',9.546D-14
         ENDIF
      ENDIF
      WRITE(*,'(A,I5,A)') 'Rate calculation will iteratively remove minima with ',NCONNMIN,' connections or fewer'
      IF (STARTFROMPATH) THEN
         WRITE(*,'(A)') 'Reading initial pathway from ' // TRIM(ADJUSTL(PATHNAME)) //' in min-sad-min format'
         PRINT '(2(A,I8))','The A minimum is number ',STARTMINA,' and the B minimum is ',STARTMINB
      ENDIF
      IF (ADDPATH) THEN
         WRITE(*,'(A)') 'Adding stationary points from ' // TRIM(ADJUSTL(PATHNAME)) //' in min-sad-min format'
      ENDIF
      WRITE(*,'(A)') 'path.info files will be read as min-sad-min triples'
      IF (CONNECTREGIONT) THEN
         PRINT '(A,G12.3,A,2I8)','Pairs of minima will be chosen based upon pairs separated by less than ',  
     &               CONNECTDIST,' starting from minima ',CONNECTMIN1,CONNECTMIN2
         IF (CONNECTMIN2F.GT.-1) PRINT '(A,I8)','Loop over second minimum will end at minimum ',CONNECTMIN2F
      ELSEIF (DIJINITT) THEN
         WRITE(*,'(A)') 'Performing Dijkstra analysis to perform initial connection'
         IF (INDEXCOSTFUNCTION) THEN
            WRITE(*,'(A)') 'Using index of minima for edge weights'
         ELSEIF (EXPCOSTFUNCTION) THEN
            WRITE(*,'(A)') 'Using exponential distance for edge weights'
         ELSEIF (.NOT.INTCONSTRAINTT) THEN
            WRITE(*,'(A,I4,A)') 'Using distance to the power ',CostFunctionPower,' for edge weights'
         ENDIF
         WRITE(*,'(A,G20.10)') 'Transition states will be rejected above threshold energy ',TSTHRESH
         IF (INITIALDIST) THEN
            WRITE(*,'(A,G20.10,A)') 'Pair distances will be intialised to ',DISBOUND,' and calculated when needed'
         ELSE
            WRITE(*,'(A,I8)') 'Number of closest neighbours to save=',PAIRDISTMAX
         ENDIF
      ELSEIF (DIJINITFLYT) THEN
         WRITE(*,'(A)') 'Performing Dijkstra analysis to perform initial connection with weights calculated on-the-fly'
         IF (EXPCOSTFUNCTION) THEN
            WRITE(*,'(A)') 'Using exponential distance for edge weights'
         ELSEIF (.NOT.INTCONSTRAINTT) THEN
            WRITE(*,'(A,I4,A)') 'Using distance to the power ',CostFunctionPower,' for edge weights'
         ENDIF
         WRITE(*,'(A,G20.10)') 'Transition states will be rejected above threshold energy ',TSTHRESH
         WRITE(*,'(A,I8)') 'Number of closest neighbours to save=',PAIRDISTMAX
      ELSEIF (DIJPAIRT) THEN
         WRITE(*,'(A)') 'Performing Dijkstra analysis to generate connnection pairs'
         IF (BARRIERSORT) THEN
            WRITE(*,'(A)') 'Pairs will be sorted according to barrier height'
         ELSE
            WRITE(*,'(A)') 'Pairs will be sorted according to energy of the minimum'
         ENDIF
!     ELSEIF (PTAUT) THEN
!        PRINT '(A,I8,A)','Minima will be chosen for connection to minima on the fastest path based on p^eq * waiting time'
      ELSEIF (FREEPAIRT) THEN
         PRINT '(A,I8,A)','Pairs of minima will be chosen based on free energy barriers between regrouped minima'
         PRINT '(3(A,G20.10))','Regrouping barrier threshold=',REGROUPFREETHRESH,' with minima above ',FREETHRESH,' excluded'
         IF (NPAIRFRQ.LT.1) THEN
            PRINT '(A,F12.3)','Pairs will not be recalculated unless we run out of possibilities'
         ELSE
            PRINT '(A,I10,A)','Pairs will be recalculated every ',NPAIRFRQ,' cycles'
         ENDIF
         IF (CHARMMT.AND.(PLANCK.EQ.1.0D0)) THEN
            PRINT '(A,G20.10)','ERROR - regrouping with CHARMM and Planck constant=',PLANCK
            STOP
         ENDIF
      ELSEIF (SHORTCUTT) THEN
         IF (BARRIERSHORT) THEN
            PRINT '(A,I8,A)','Pairs of minima will be chosen for minima separated by fewer than',MINSEP,  
     &                  ' steps on either side of the largest barriers in the best SS path'
         ELSEIF (RATESHORT) THEN
            PRINT '(A,I8,A)','Pairs of minima will be chosen for minima separated by fewer than',MINSEP,  
     &                  ' steps on either side of the smallest rate constants in the best SS path'
         ELSE
            PRINT '(A,I8,A)','Pairs of minima will be chosen for the closest minima separated by at least ',MINSEP,  
     &                  ' steps in the best SS path identified by Dijkstra'
         ENDIF
         IF (NPAIRFRQ.LT.1) THEN
            PRINT '(A,F12.3)','Pairs will not be recalculated unless we run out of possibilities'
         ELSE
            PRINT '(A,I10,A)','Pairs will be recalculated every ',NPAIRFRQ,' cycles'
         ENDIF
      ELSEIF (UNTRAPT) THEN
         PRINT '(A,I8,A)','Pairs of minima will be chosen according to their lowest barrier to a minimum in a product superbasin'
         PRINT '(A,F15.5,A)','Connections will not be attempted to minima higher than ',EUNTRAPTHRESH, 
     &                      ' energy units above the first minimum'
         IF (NPAIRFRQ.LT.1) THEN
            PRINT '(A,F12.3)','Pairs will not be recalculated unless we run out of possibilities'
         ELSE
            PRINT '(A,I10,A)','Pairs will be recalculated every ',NPAIRFRQ,' cycles'
         ENDIF
      ELSEIF (PERSISTT) THEN
         PRINT '(A)','Persistence analysis'
         PRINT '(A,F15.5,A)','Minima will not be included higher than ',EUNTRAPTHRESH, 
     &                      ' energy units above the global minimum'
         IF (PERSISTAPPROXT) THEN
            PRINT '(A)','Superbasin approach will be used to estimate persistences'
         ELSE
            PRINT '(A)','Exact calculation of persistences'
         ENDIF
         IF (ALLCOMPONENTST) PRINT '(A)','Persistence information will be reported for all valid minima'
      ELSEIF (NEWCONNECTIONST) THEN
         PRINT '(A)','Database will be extended using single-ended transition state searches'
         PRINT '(A,F12.3)','Minima will be perturbed using random Cartesian displacements up to magnitude ',PERTVALUE
         PRINT '(4(A,I10))','Minima will be perturbed starting from ',CONNMINSTART,
     &                     ' in order until they have ',CONNECTIONS,
     &                     ' connections up to a maximum ',MAXTSATTEMPTS,' attempts each'
      ELSEIF (USEPAIRST) THEN
         PRINT '(A,A)','Pairs of minima will be chosen from the sequence in file ',TRIM(ADJUSTL(USEPAIRSFILE))
      ELSEIF (CONNECTPAIRST) THEN
         PRINT '(A,A)','Pairs of minima will be taken from the file ',TRIM(ADJUSTL(CONNECTPAIRSFILE))
      ELSEIF (NATTEMPT.GT.0) THEN
         PRINT '(A,F12.3)','Pairs of minima will be chosen based upon a minimum committor probability difference of ',PSCALE
         IF (NPAIRFRQ.LT.1) THEN
            PRINT '(A,F12.3)','Pairs will not be recalculated unless we run out of possibilities'
         ELSE
            PRINT '(A,I10,A)','Pairs will be recalculated every ',NPAIRFRQ,' cycles'
         ENDIF
      ENDIF
      IF (DUMMYTST) THEN
         IF (BHINTERPT) PRINT '(A)','Dummy ts entries will be created for closest pairs of minima found via basin-hopping.'
         IF (BISECTT) PRINT '(A)','Dummy ts entries will be created for pairs of minima found by bisection.'
      ENDIF
      IF (CLOSEFILEST) PRINT '(A)','The min.data and ts.data files will not be kept open'
      IF (INTLJT) THEN
         PRINT '(A,F15.5)',   'Using interpLJ potential metric'
         PRINT '(A,F20.10)',  'minimum distance difference for internal minimum=',INTLJDEL
         PRINT '(A,F20.10)',  'multiplying factor for internal minimum penalty function=',INTLJEPS
      ENDIF
      IF (INTCONSTRAINTT) THEN
         PRINT '(A,F15.5)',   'Using constraint potential metric with absolute distance change tolerance ',INTCONSTRAINTTOL
         PRINT '(A,F15.5)',   '   constraint spring constant=',INTCONSTRAINTDEL
         PRINT '(2(A,F15.5))','   repulsion factor between unconstrained atoms=',INTCONSTRAINTREP
         PRINT '(A,F15.5)',   '   cutoff for repulsion will be the minimum of ',INTCONSTRAINREPCUT
         PRINT '(A)',         '   and the shortest distance in the end points'
         PRINT '(A,I10)',      '   maximum separation of atoms in sequence for constraint=',INTCONSEP
         PRINT '(A,I10)',      '   minimum separation of atoms in sequence for repulsion=',INTREPSEP
         PRINT '(A,I10)',      '   maximum number of constraints per atom=',MAXCONUSE
         IF (CHECKCONINT) THEN
            PRINT '(A,I10)',   '   adding terms for constraint internal minima'
         ELSE
            PRINT '(A,I10)',   '   not adding terms for constraint internal minima'
         ENDIF
      ENDIF 
      IF ((AMBERT.OR.CHARMMT.OR.OPEPT).AND.(EVCUT .LT. 1.0D0)) THEN
! jmc resetting a small evcut (tolerance on zero used to check Hessian
! e-values from path.info files in getallpaths.f) to a more appropriate value for SI units 
         EVCUT = 1.0D20
         PRINT '(A,G8.2,A)','Resetting EVCUT to ',EVCUT,
     &          ' for SI units; check it is appropriate for your system'
      ENDIF

      WRITE(*,'(A)') '*************************************************************************************************'

      IF (REGROUPT.AND.(NATTEMPT.GT.0)) THEN
         PRINT '(A)','ERROR - regroup can only be called once, set CYCLES=0'
!        STOP
      ENDIF
      IF (RELATIVEET.AND.(NATTEMPT.GT.0)) THEN
         PRINT '(A)','ERROR - relative energy can only be subtracted for static databases. set CYCLES=0'
!        STOP
      ENDIF
      IF (DOCKT) THEN
         IF (SLURMT) THEN
            PRINT '(A)','main> Job submission using srun is not implemented for keyword DOCK. '
            STOP
         END IF
         CALL DOCK
      END IF
      IF (SIST) THEN
         CALL SETUP_SIS
      ELSE
         CALL SETUP
      ENDIF
      IF (PRINTSUMMARYT) CALL PRINTSUMMARY
      IF (CHECKMINT) THEN
         IF (CHECKSPS.LE.0) CHECKSPS=1
         IF (CHECKSPF.LE.0) CHECKSPF=NMIN
         PRINT '(2(A,I10))','Checking minima from ',CHECKSPS,' to ',CHECKSPF
         INQUIRE(FILE='odata.checksp',EXIST=YESNOT)
         IF (.NOT.YESNOT) THEN
            WRITE(*,'(A)') 'main> ERROR - no odata.checksp file present'
            STOP
         ENDIF
         IF (GETMINFRQST) THEN
            WRITE(*,'(A)') 'main> Dump frequencies for all minima.'
         ENDIF
      ENDIF
      IF (CHECKTST) THEN
         IF (CHECKSPS.LE.0) CHECKSPS=1
         IF (CHECKSPF.LE.0) CHECKSPF=NTS
         PRINT '(2(A,I10))','Checking transition states from ',CHECKSPS,' to ',CHECKSPF
         INQUIRE(FILE='odata.checksp',EXIST=YESNOT)
         IF (.NOT.YESNOT) THEN
            WRITE(*,'(A)') 'main> ERROR - no odata.checksp file present'
            STOP
         ENDIF
         IF (GETTSFRQST) THEN
            WRITE(*,'(A)') 'main> Dump frequencies for all transition states.'
         ENDIF
      ENDIF
      IF (ADDMINXYZT) THEN
         PRINT '(2(A,I10))','main> Adding minima from file ' // TRIM(ADJUSTL(ADDMINXYZNAME))
         PRINT '(2(A,I10))','main> after running OPTIM using odata.addminxyz template to generate min.data.info'
         INQUIRE(FILE='odata.addminxyz',EXIST=YESNOT)
         IF (.NOT.YESNOT) THEN
            WRITE(*,'(A)') 'main> ERROR - no odata.addminxyz file present'
            STOP
         ENDIF
         LUNIT=GETUNIT()
         OPEN(LUNIT,FILE=TRIM(ADJUSTL(ADDMINXYZNAME)),STATUS='OLD')
         NMINADDXYZ=0
         DO 
            READ(LUNIT,*,END=765) NDUMMY
            IF (NDUMMY.NE.NATOMS) THEN
               PRINT '(A,I10,A)','main> ERROR *** number of atoms in xyz file=',NDUMMY
               PRINT '(A,I10,A)','main> ERROR *** should be ',NATOMS
               STOP
            ENDIF
            READ(LUNIT,'(A)',END=765) SDUMMY
            READ(LUNIT,*,END=765) (SDUMMY,XDUMMY,XDUMMY,XDUMMY,J1=1,NATOMS)
            NMINADDXYZ=NMINADDXYZ+1
         ENDDO
765      CONTINUE
         PRINT '(A,I8)','main> Number of minima in xyz file ' // TRIM(ADJUSTL(ADDMINXYZNAME)) // ' is ',
     &                   NMINADDXYZ
         CLOSE(LUNIT)
      ENDIF
      CALL CPU_TIME(TNEW)
      WRITE(*,'(A,F15.5)') 'main> time spent in setup=',TNEW-TINIT

      IF ((NATTEMPT.GT.0).OR.NEWCONNECTIONST) THEN
         IF (NEWCONNECTIONST) THEN
            INQUIRE(FILE='odata.tspath',EXIST=YESNOT)
            IF (.NOT.YESNOT) THEN
               WRITE(*,'(A)') 'main> ERROR - no odata.tspath file present'
               STOP
            ENDIF
         ELSEIF (.NOT.(BHINTERPT.OR.BISECTT.OR.CHECKSPT.OR.ADDMINXYZT)) THEN
            INQUIRE(FILE='odata.connect',EXIST=YESNOT)
            IF (.NOT.YESNOT) THEN
               WRITE(*,'(A)') 'main> ERROR - no odata.connect file present'
               STOP
            ENDIF
         ENDIF
         IF ((CONNECTIONS.GT.2).AND.(.NOT.NEWCONNECTIONST)) THEN
            CALL CYCLE
         ELSE
            CALL CYCLE2
         ENDIF
         OPEN(UNIT=1,FILE='commit.data',STATUS='UNKNOWN')
         WRITE(1,'(G20.10)') GPFOLD(1:NMIN)
         CLOSE(1)
         IF ((NPAIRDONE.GT.0).AND.(.NOT.DUMMYRUNT)) THEN
            OPEN(UNIT=1,FILE='pairs.data',STATUS='UNKNOWN')
            WRITE(1,'(2I8)') (PAIR1(J1),PAIR2(J1),J1=1,NPAIRDONE)
            CLOSE(1)
         ENDIF
         IF (NMINDONE.GT.0) THEN
            OPEN(UNIT=1,FILE='min.done',STATUS='UNKNOWN')
            WRITE(1,'(I8)') (MINDONE(J1),J1=1,NMINDONE)
            CLOSE(1)
         ENDIF
      ELSE IF (TFOLDT) THEN          ! just run TFOLD subroutine
         CALL TFOLD
      ELSE IF (ARNOLDIT) THEN           ! just run diagonalise2 subroutine
         CALL REWEIGHT
      ELSE IF (DIAGT) THEN           ! just run diagonalise2 subroutine
         CALL DIAGONALISE2
      ELSE IF (NGTT) THEN            ! just run NGT subroutine
         IF (RATESCYCLET) THEN
            CALL GETRATESCYCLE
         ELSE
            CALL NGT
         ENDIF
      ELSE IF (GTT.OR.GT2T) THEN     ! just run GT subroutine
         CALL GT
      ELSE IF (NPFOLD.GT.0) THEN     ! just run PFOLD subroutine
         CALL PFOLD
         OPEN(UNIT=1,FILE='commit.data',STATUS='UNKNOWN')
         WRITE(1,'(G20.10)') GPFOLD(1:NMIN)
         CLOSE(1)
      ELSE IF (KMCT) THEN            ! just run KMC 
!
!  REGROUP and REGROUPFREE change NMINA, NMINB, LOCATIONA, LOCATIONB, so save the values and reset
!  to call GT more than once.
!  In fact, REGROUPFREE changes NMIN, NTS, etc. so we must stop after such a run or
!  do a complete reset somehow. This is done explicitly in routines like getppair
!  using the SAVESTATE module. Not needed here because we assume that NGT cannot be
!  called more than once!
!
         IF (REGROUPFREET.OR.REGROUPFREEABT) THEN
            CALL GETNCONN ! must call this first to set NCONNMAX - used for declarations in REGROUPFREE2
            ALLOCATE(FREEMINLIST(NMIN),FREEMINPOINT(0:NMIN+1))
            CALL REGROUPFREE2(.FALSE.,1,FREEMINLIST,FREEMINPOINT,NAVAIL)
            DEALLOCATE(FREEMINLIST,FREEMINPOINT)
            TSTHRESH=HUGE(1.0D0) ! free energy scale will be different from PE, so must reset
                                 ! before calling GETNCONN
            MAXBARRIER=HUGE(1.0D0)
         ENDIF
   
         IF (REGROUPRATET.OR.REGROUPPET) THEN
            CALL REGROUPFREE
            TSTHRESH=HUGE(1.0D0) ! free energy scale will be different from PE, so must reset
                                 ! before calling GETNCONN
            MAXBARRIER=HUGE(1.0D0)
         ENDIF
!
!  Must not call GETNCONN here or MAXCONN may be reset to a lower value, which messes things up.
!
         CALL GETNCONN
!
!  REGROUP should give us the corrected POINTER values even if we have run REGROUPFREE.
!
         ALLOCATE(MINMAP(NMIN))
         CALL REGROUP(MINMAP)
         CALL KMC                    ! standard stochastic KMC
      ELSE IF (DIJKSTRAT .OR. KSHORTESTPATHST) THEN ! run shortest path analysis
         NMINSAVE=NMIN
         IF (REGROUPFREET.OR.REGROUPFREEABT) THEN
            CALL GETNCONN ! must call this first to set NCONNMAX - used for declarations in REGROUPFREE2
            ALLOCATE(FREEMINLIST(NMIN),FREEMINPOINT(0:NMIN+1))
            CALL REGROUPFREE2(.FALSE.,1,FREEMINLIST,FREEMINPOINT,NAVAIL)
            DEALLOCATE(FREEMINLIST,FREEMINPOINT)
            TSTHRESH=HUGE(1.0D0) ! free energy scale will be different from PE, so must reset
                                 ! before calling GETNCONN
            MAXBARRIER=HUGE(1.0D0)
            CALL GETNCONN
!  
!  REGROUP should give us the corrected POINTER values even if we have run REGROUPFREE.
!  
            ALLOCATE(MINMAP(NMIN))
            CALL REGROUP(MINMAP)
         ELSEIF (REGROUPRATET.OR.REGROUPPET) THEN
            CALL REGROUPFREE
            TSTHRESH=HUGE(1.0D0) ! free energy scale will be different from PE, so must reset
                                 ! before calling GETNCONN
            MAXBARRIER=HUGE(1.0D0)
            CALL GETNCONN
!  
!  REGROUP should give us the corrected POINTER values even if we have run REGROUPFREE.
!  REGROUP also reorders the stationary points so that the most connected minima appear
!  at the end in order to accelerate GT.
!  
            ALLOCATE(MINMAP(NMIN))
            CALL REGROUP(MINMAP)
         ELSE
            CALL GETNCONN
         ENDIF
         IF (.NOT.ALLOCATED(MINMAP)) THEN ! allow for the case of no regrouping
            ALLOCATE(MINMAP(NMIN))
            DO J1=1,NMIN
               MINMAP(J1)=J1
            ENDDO
         ENDIF
         IF (DIJKSTRAT) THEN
            CALL DIJKSTRA(NWORST,.FALSE.,0,NMINSAVE,MINMAP)
         ELSEIF (KSHORTESTPATHST) THEN
            CALL KSHORTESTPATHS(NWORST,.FALSE.,0,NMINSAVE,MINMAP)
         ENDIF
         IF (ALLOCATED(MINMAP))  DEALLOCATE(MINMAP)
         IF (ALLOCATED(MINGROUP))  DEALLOCATE(MINGROUP)
      ELSE IF (REGROUPPERSISTT) THEN
         IF (.NOT.PERSISTT) THEN
            PRINT '(A)','main> ERROR *** regrouppersist specified without persist keyword'
            STOP
         ENDIF
         CALL PERSISTENCE
         CALL GETNCONN ! must call this first to set NCONNMAX - used for declarations in REGROUPFREE2
         ALLOCATE(FREEMINLIST(NMIN),FREEMINPOINT(0:NMIN+1))
         CALL REGROUPPERSIST(.FALSE.,1,FREEMINLIST,FREEMINPOINT,NAVAIL)
         DEALLOCATE(FREEMINLIST,FREEMINPOINT)
         DEALLOCATE(MINGROUP)
      ELSE IF (REGROUPFREET.OR.REGROUPFREEABT) THEN    ! just run new REGROUPFREE2 subroutine
         CALL GETNCONN ! must call this first to set NCONNMAX - used for declarations in REGROUPFREE2
         ALLOCATE(FREEMINLIST(NMIN),FREEMINPOINT(0:NMIN+1))
         CALL REGROUPFREE2(.FALSE.,1,FREEMINLIST,FREEMINPOINT,NAVAIL)
         DEALLOCATE(FREEMINLIST,FREEMINPOINT)
         DEALLOCATE(MINGROUP)
      ELSE IF (PERSISTT) THEN          ! just run persistence analysis
         CALL PERSISTENCE
      ELSE IF (REGROUPRATET) THEN    ! just run new REGROUPFREE subroutine
         CALL REGROUPFREE
         DEALLOCATE(MINGROUP)
      ELSE IF (RFMULTIT) THEN    ! run REGROUPFREE2 for a series of temperatures
         CALL GETNCONN ! must call this first to set NCONNMAX - used for declarations in REGROUPFREE2
         CALL RFCYCLE
         IF (ALLOCATED(MINGROUP)) DEALLOCATE(MINGROUP)
      ELSE IF (RFKMCT) THEN    ! run REGROUPFREE2 for a series of temperatures
         CALL GETNCONN ! must call this first to set NCONNMAX - used for declarations in REGROUPFREE2
         CALL RFKMC
         IF (ALLOCATED(MINGROUP)) DEALLOCATE(MINGROUP)
      ELSE IF (REGROUPPET) THEN      ! just run new REGROUPFREE subroutine
         CALL REGROUPFREE
         DEALLOCATE(MINGROUP)
      ELSE IF (KMCCOMMITT) THEN      ! just run new KMC assuming equilibrium in B region
         CALL GETNCONN
         CALL KMCCOMMIT
      ENDIF
      IF (DIJPAIRT) THEN
         OPEN(UNIT=1,FILE='ts.attempts',STATUS='UNKNOWN')
         WRITE(1,'(I8)') TSATTEMPT(1:NTS)
         CLOSE(1)
      ENDIF

      DEALLOCATE(FRQS,MASS,ZSYMBOL,RESLABEL, 
     &         ATOMLABEL,RESNUMBER,EMIN,FVIBMIN,PFMIN, 
     &         IXMIN,IYMIN,IZMIN,GPFOLD,ETS,FVIBTS,KPLUS,KMINUS,NEGEIG,
     &         IXTS,IYTS,IZTS,HORDERMIN,TOPPOINTER,HORDERTS, 
     &         PLUS,MINUS,POINTERM,POINTERP,PAIR1,PAIR2)
      IF (ALLOCATED(NCONN)) DEALLOCATE(NCONN)
      IF (ALLOCATED(DMIN1)) DEALLOCATE(DMIN1,DMIN2)
      IF (ALLOCATED(TSATTEMPT)) DEALLOCATE(TSATTEMPT)
      IF (ALLOCATED(MINDONE)) DEALLOCATE(MINDONE)
      IF (ALLOCATED(PAIRDIST)) DEALLOCATE(PAIRDIST)
      IF (ALLOCATED(ALLPAIRS)) DEALLOCATE(ALLPAIRS)
      IF (DUMMYTST) DEALLOCATE(MINDISTMIN)

      CALL CPU_TIME(TNEW)
      IF (NATTEMPT.GT.0) WRITE(*,'(A,G15.5,A)') 'main> CPU time spent iterating committor probability=',TPFOLD,' s'
      IF (TFOLDT) WRITE(*,'(A,G15.5,A)') 'main> CPU time spent iterating waiting times=',TTFOLD,' s'
      IF (GTT) WRITE(*,'(A,G15.5,A)') 'main> CPU time spent in GT                          =',TGT,' s'
      IF (NGTT) WRITE(*,'(A,G15.5,A)') 'main> CPU time spent in NGT                          =',TGT,' s'
      IF (DIJKSTRAT) WRITE(*,'(A,G15.5,A)') 'main> CPU time spent in Dijkstra                    =',TDIJKSTRA,' s'
      IF (KSHORTESTPATHST) WRITE(*,'(A,G15.5,A)') 'main> CPU time spent in kshortestpaths              =',TKSHORTESTPATHS,' s'
      IF (CONNECTREGIONT) WRITE(*,'(A,G15.5,A)') 'main> CPU time spent in connectdist                 =',TCONNECTDIST,' s'
      WRITE(*,'(A,G15.5,A)')  'main> total CPU time spent in PATHSAMPLE            =',TNEW-TINIT,' s'

      STOP
      END
C
C  Subroutine GETNCONN sets up array NCONN containing the number of 
C  connections for each minimum after pruning according to the value of
C  NCONNMIN.
C
      SUBROUTINE GETNCONN
      USE COMMONS, ONLY: NMIN,NTS,NCONN,PLUS,MINUS,NCONNMIN,NCONNMAX,MAXMIN,DEBUG,ETS, 
     &                  KPLUS, KMINUS, EMIN, SHANNONT
      IMPLICIT NONE
      INTEGER J1, PNCONNECTED, NCONNECTED, NZERO, JMAX
      LOGICAL, ALLOCATABLE :: CONNECTED(:) ! reallocate MAXMIN when used
      LOGICAL :: DEADTS
      DOUBLE PRECISION :: CUT_UNDERFLOW=-300.0D0
C
C  Record the number of connections for each minimum in NCONN.
C  NCONN is the number of connections to minima with more
C  than NCONNMIN connections.
C
      NCONNECTED=0
      IF (ALLOCATED(NCONN)) DEALLOCATE(NCONN)
      ALLOCATE(CONNECTED(MAXMIN),NCONN(MAXMIN))
      DO J1=1,NMIN
         CONNECTED(J1)=.TRUE.
      ENDDO
11    DO J1=1,NMIN
         NCONN(J1)=0
      ENDDO
      PNCONNECTED=NCONNECTED
      DO J1=1,NTS
C JMC n.b. don't apply the nconnmin criteria at this point, hence the huge(1) 's in place of NCONN() for the plus and minus minima.        
         CALL CHECKTS(ETS(J1),EMIN(PLUS(J1)),EMIN(MINUS(J1)),KPLUS(J1),KMINUS(J1),HUGE(1),HUGE(1), 
     &                PLUS(J1),MINUS(J1),.FALSE.,CUT_UNDERFLOW,DEADTS)
         IF ((.NOT.DEADTS).AND.(PLUS(J1).NE.MINUS(J1))) THEN
            IF (CONNECTED(MINUS(J1))) NCONN(PLUS(J1))=NCONN(PLUS(J1))+1
            IF (CONNECTED(PLUS(J1)))  NCONN(MINUS(J1))=NCONN(MINUS(J1))+1
         ENDIF
      ENDDO
      NCONNECTED=0
      DO J1=1,NMIN
         CONNECTED(J1)=.FALSE.
         IF (NCONN(J1).GT.NCONNMIN) THEN
            CONNECTED(J1)=.TRUE.
            NCONNECTED=NCONNECTED+1
         ENDIF
      ENDDO
      IF (NCONNECTED.NE.PNCONNECTED) GOTO 11

!     DO J1=1,NMIN
!        IF (DEBUG) WRITE(*,'(A,I10,A,I10)') 'getnconn> number of connections for minimum ',J1,' is ',NCONN(J1)
!     ENDDO

      NCONNMAX=NCONN(1)
      NZERO=0
      IF (NCONN(1).EQ.0) NZERO=1
      JMAX=1
      DO J1=2,NMIN
         IF (NCONN(J1).EQ.0) NZERO=NZERO+1
         IF (NCONN(J1).GT.NCONNMAX) THEN
            NCONNMAX=NCONN(J1)
            JMAX=J1
         ENDIF
      ENDDO
      IF (.NOT.SHANNONT) WRITE(*,'(4(A,I10))') 'getnconn> max connections: ',NCONNMAX,' for min ',JMAX,' # of zeros=',NZERO,
     &                     ' after removing minima with < ',NCONNMIN+1
      DEALLOCATE(CONNECTED)

      RETURN
      END
