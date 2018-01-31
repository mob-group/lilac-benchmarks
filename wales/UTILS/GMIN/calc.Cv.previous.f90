
!  Calculate Cv curve from weights and visit statistics
!  for replicas at different temperatures REPT(J1).
!  Minimise chi^2 statistic to extract best probabilities.
!
PROGRAM CALCCV
IMPLICIT NONE
DOUBLE PRECISION TMIN, TMAX, TINT, DIFF, CHI2, RMS, VTOTAL2, BESTRMS, NEGLECT, WTOT, VTOT, FRACTION
DOUBLE PRECISION QE, XFIT(3), KFAC, SIGMAX
DOUBLE PRECISION, ALLOCATABLE :: ENERGY(:), REPT(:), VTOTAL(:), VAR(:), VAR2(:), WIJ(:,:), BESTVAR(:)
DOUBLE PRECISION, ALLOCATABLE :: WEIGHT(:,:), OLDVAR(:), QENERGY(:)
DOUBLE PRECISION, ALLOCATABLE :: GRAD(:), DGRAD(:), EVAR(:), FITA(:), FITB(:), FITC(:), REPWEIGHTS(:)
DOUBLE PRECISION, ALLOCATABLE :: SUMVISITS(:,:), VISITS(:,:), VSUM(:)
!
! GMIN is currently dumping VISITS data as real but VISITS2 as integer.
! Should probably change this!
!
DOUBLE PRECISION, ALLOCATABLE :: BAV(:)
INTEGER, ALLOCATABLE :: NMINQBIN(:)
DOUBLE PRECISION, ALLOCATABLE :: XPIQBIN(:)
INTEGER, ALLOCATABLE :: VISITS2(:,:,:)
INTEGER, ALLOCATABLE :: VARMAP2(:,:), VISITST(:,:)
INTEGER, ALLOCATABLE :: LOWESTDIRECT(:), HIGHESTDIRECT(:), LOWESTQUENCH(:), HIGHESTQUENCH(:), BESTBINS(:)
INTEGER, ALLOCATABLE :: MOSTVISITED(:)
DOUBLE PRECISION, ALLOCATABLE :: MAXVISITS(:)
DOUBLE PRECISION CHI2PLUS, CHI2MINUS, TOL, MAXBHSTEP, DPRAND, DELTA, DUMMY, &
  &              USEQBINS, PEMIN, QHISTINT, BSCALE, HSUM, BSUM, XDUMMY, RATMAX, FIRSTFAC, XXDUMMY
INTEGER NBINS, J1, NTEMP, J2, KAPPA, NREP, NVAR, NEVAR, NCOUNT, MUPDATES, ITMAX, ITDONE, NHOPS, RNDSEED, &
  &     J3, NQBINS, J4, NVAROLD, MERGEI, MERGEQ, MGAMMA, NDUMMY2, &
  &     PEBINMIN, PEBINMAX, MINFIT, NDUMMY, NBINSAVE, NQBINSAVE, NREPMIN, NDUMMY3, &
  &     NRESMIN, LOWEST
COMMON /PELIMITS/ PEBINMIN, PEBINMAX
CHARACTER(LEN=80) FSTRING, FNAME, FSTRING2, PTSTART, BSSTART, BSPTFINISH, REPSTRING
CHARACTER(LEN=1) DUMMYSTRING
LOGICAL, ALLOCATABLE :: ALLZERO(:), ALLZERO2(:,:), ALLZEROQ(:)
LOGICAL CONVERGED, RESERVOIRT, FIRST, YESNO, FIXB
DOUBLE PRECISION, ALLOCATABLE :: EMIN(:), FVIBMIN(:), PFMIN(:), IXMIN(:), IYMIN(:), IZMIN(:)
INTEGER, ALLOCATABLE :: HORDERMIN(:)

OPEN(UNIT=10,FILE='Cv.data',STATUS='OLD')
! READ(10,*) NREPMIN, NREP, NBINS, NQBINS, TMIN, TMAX, NTEMP, KAPPA, PTSTART, BSSTART, BSPTFINISH, NEGLECT, MERGEI, MERGEQ, SIGMAX
READ(10,*) NREPMIN, NREP, NBINS, NQBINS, TMIN, TMAX, NTEMP, KAPPA, PTSTART, BSSTART, BSPTFINISH, NEGLECT, SIGMAX, RATMAX
ALLOCATE(REPWEIGHTS(NREP))
REPWEIGHTS(1:NREP)=1.0D0
READ(10,*,END=543) REPWEIGHTS(1:NREP)
CLOSE(10)

543 CONTINUE
PRINT '(A)','weighting for visits corresponding to different replicas:'
PRINT '(8F12.4)',REPWEIGHTS(1:NREP)

PRINT '(3(A,I4))','replicas=',NREP,' # pe bins=',NBINS,' # Q bins=',NQBINS
NBINSAVE=NBINS; NQBINSAVE=NQBINS
MERGEI=1
FIXB=.FALSE.
!
! PTSTART 'na' means subtract nothing for PT analysis
! PTSTART 'number' means subtract the Visits statistics for this step in PT analysis
! The PT analysis goes up to the same step as specified by BSPTFINISH
! BSSTART 'na' means subtract nothing for quench analysis
! BSPTFINISH 'all' means use the final Visits.his and Visits2.his files
! BSPTFINISH 'number' means use step number as the final point in the fitting
!
! NBINS=NBINS/MERGEI
! NQBINS=NQBINS/MERGEQ
! IF (MERGEI.NE.1) PRINT '(A,I4)','replicas will be ignored below number ',NREPMIN
! IF (MERGEI.NE.1) PRINT '(3(A,I4))','pe bins will be merged in groups of ',MERGEI,' to give ',NBINS
! ! IF (MERGEQ.NE.1) PRINT '(3(A,I4))','quench bins will be merged in groups of ',MERGEQ,' to give ',NQBINS
IF (NREPMIN.GT.1) PRINT '(A,I4)','replicas will be ignored below number ',NREPMIN
PRINT '(A,2G15.5,A,I8)','minimum and maximum T for Cv calculation=',TMIN,TMAX,' number of Cv data points=',NTEMP
PRINT '(A,I8)','number of degrees of freedom=',KAPPA
IF (BSPTFINISH(1:3).EQ.'all') THEN
   PRINT '(A)','Using statistics from final Visits files'
   PRINT '(A)','Using statistics from final Visits2 files'
ELSE
   IF (PTSTART.EQ.'na') THEN
      PRINT '(A,A)','Using statistics from Visits files at step ',TRIM(ADJUSTL(BSPTFINISH))
   ELSE
      PRINT '(4A)','Using statistics from Visits files at step ',TRIM(ADJUSTL(BSPTFINISH)), &
   &                ' and subtracting step ',TRIM(ADJUSTL(PTSTART))
   ENDIF
   IF (BSSTART.EQ.'na') THEN
      PRINT '(A,A)','Using statistics from Visits2 files at step ',TRIM(ADJUSTL(BSPTFINISH))
   ELSE
      PRINT '(4A)','Using statistics from Visits2 files at step ',TRIM(ADJUSTL(BSPTFINISH)), &
   &                ' and subtracting step ',TRIM(ADJUSTL(BSSTART))
   ENDIF
ENDIF
PRINT '(A,G20.10,A)','Bins will be neglected if the number of visits is less than ',NEGLECT,' times the largest value'
! PRINT '(A,I8)','Maximum number of pe bins to be used in fits=',MAXBIN
MINFIT=20
PRINT '(A,I8)','Minimum number of data points for extrapolated fits=',MINFIT
PRINT '(A,G20.10)','Tolerance for residual standard deviation in final fits=',SIGMAX
PRINT '(A,G20.10)','Maximum ratio allowed for shifted fitted DoS=',RATMAX

ALLOCATE(REPT(NREP)) ! temperature of replicas

ALLOCATE(ENERGY(NBINS),ALLZERO(NBINS),ALLZEROQ(NQBINS),ALLZERO2(NBINS,NQBINS),VTOTAL(NREP))
ALLOCATE(VISITS2(NREP,NBINS,NQBINS), QENERGY(NQBINS),VISITST(NBINS,NQBINS))
ALLOCATE(WEIGHT(NBINS,NREP),VISITS(NBINS,NREP),VARMAP2(NBINS,NQBINS))
ALLOCATE(SUMVISITS(NBINS,NQBINS),FITA(NQBINS),FITB(NQBINS),FITC(NQBINS))
ALLOCATE(MOSTVISITED(NBINS))

!
! Second column in Visits.his.<n>  gives number of visits to pe bin
! Second column in Visits2.his.<n> gives number of visits to pe bin that quench to a given q bin
!        over starting pe, pe visited during quench, and energy quenched to
!
!
!!!!!!!!!!!!!!! Read data from Visits.his.<n> files !!!!!!!!!!!!!!!!
!
VISITS(1:NBINS,1:NREP)=0
WEIGHT(1:NBINS,1:NREP)=0.0D0
VTOTAL2=0
DO J1=1,NREP ! read in data for all replicas regardless of NREPMIN
   WRITE(REPSTRING,'(I6)') J1
   IF (BSPTFINISH(1:3).EQ.'all') THEN
      WRITE(FSTRING,'(A)') TRIM(ADJUSTL(REPSTRING)) // '/Visits.his'
   ELSE
      WRITE(FSTRING,'(A)') TRIM(ADJUSTL(REPSTRING)) // '/Visits.his.' // TRIM(ADJUSTL(BSPTFINISH))
   ENDIF
   IF (PTSTART(1:2).NE.'na') WRITE(FSTRING2,'(A)') TRIM(ADJUSTL(REPSTRING)) // '/Visits.his.' // TRIM(ADJUSTL(PTSTART)) 
   INQUIRE(FILE=TRIM(ADJUSTL(FSTRING)),EXIST=YESNO)
   IF (.NOT.YESNO) THEN
      PRINT '(3A)','File ',TRIM(ADJUSTL(FSTRING)),' not found'
      STOP
   ENDIF
   IF (PTSTART(1:2).NE.'na') INQUIRE(FILE=TRIM(ADJUSTL(FSTRING2)),EXIST=YESNO)
   IF (.NOT.YESNO) THEN
      PRINT '(3A)','File ',TRIM(ADJUSTL(FSTRING2)),' not found'
      STOP
   ENDIF
   OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FSTRING)),STATUS='OLD')
   READ(1,*) REPT(J1)
   READ(1,*) DUMMYSTRING
   IF (PTSTART(1:2).NE.'na') THEN
      OPEN(UNIT=2,FILE=TRIM(ADJUSTL(FSTRING2)),STATUS='OLD')
      READ(2,*) REPT(J1)
      READ(2,*) DUMMYSTRING
   ENDIF
   VTOTAL(J1)=0.0D0
   DO J2=1,NBINS
      WTOT=0.0D0; VTOT=0
      DO J3=1,MERGEI
         READ(1,*) ENERGY(J2),DUMMY
         VTOT=VTOT+DUMMY
      ENDDO
      VISITS(J2,J1)=VTOT*REPWEIGHTS(J1)
      VTOTAL(J1)=VTOTAL(J1)+VISITS(J2,J1)
   ENDDO
   IF (PTSTART(1:2).NE.'na') THEN
      VTOTAL(J1)=0.0D0
      DO J2=1,NBINS
         WTOT=0.0D0; VTOT=0
         DO J3=1,MERGEI
            READ(2,*) ENERGY(J2),DUMMY
            VTOT=VTOT+DUMMY
         ENDDO
         VISITS(J2,J1)=VISITS(J2,J1)-VTOT
         VTOTAL(J1)=VTOTAL(J1)+VISITS(J2,J1)
      ENDDO
      CLOSE(2)
   ENDIF
   DO J2=MERGEI*NBINS,NBINSAVE-1 ! in case MERGEI is not an integer divisor of NBINSAVE
      READ(1,*) DUMMYSTRING
   ENDDO
   READ(1,*) DUMMYSTRING
   DO J2=1,NQBINS
      READ(1,*) QENERGY(J2)
   ENDDO
   CLOSE(1)
!  WRITE(*,'(A,I5,2(A,G20.10))') 'temperature for replica ',J1,' is ',REPT(J1),' total visits=',VTOTAL(J1)
   VTOTAL2=VTOTAL2+VTOTAL(J1)
ENDDO


!
! Do we want to use min.data?
! If RESERVOIRT is true then compute the average vibrational/permutational
! prefactors for minima in each quench bin. Use these values to fix the
! parameter B_k for quench bin k up to the number of specified quench bins.
! At the end, divide out the average vibrational/permutational prefactor
! to estimate the number of distinct minima in each quench bin.
!
RESERVOIRT=.FALSE.
INQUIRE(FILE='min.data',EXIST=YESNO)
IF (YESNO) THEN
   PRINT '(A)','keyword> min.data file located in working directory'
   OPEN(UNIT=10,FILE='min.data',STATUS='OLD')
   NRESMIN=0
   DO
      READ(10,*,END=30) DUMMY
      NRESMIN=NRESMIN+1
   ENDDO
30 WRITE(*,'(A,I10,A)') 'keyword> There are ',NRESMIN,' minima in the min.data file'
   ALLOCATE(EMIN(NRESMIN),FVIBMIN(NRESMIN),PFMIN(NRESMIN),IXMIN(NRESMIN),IYMIN(NRESMIN), &
  &         IZMIN(NRESMIN),HORDERMIN(NRESMIN))
!
! Average vibration/permutation factor and number of minima from min.data in quench bins.
!
   ALLOCATE(BAV(NQBINS),NMINQBIN(NQBINS),XPIQBIN(NQBINS))
   REWIND(10)
   PEMIN=1.0D100
   DO J1=1,NRESMIN
      READ(10,*) EMIN(J1),FVIBMIN(J1),HORDERMIN(J1),IXMIN(J1),IYMIN(J1),IZMIN(J1)
      IF (EMIN(J1).LT.PEMIN) THEN
         LOWEST=J1
         PEMIN=EMIN(J1)
      ENDIF
   ENDDO
   CLOSE(10)
   PRINT '(A,I6,A,G20.10)','keyword> Lowest minimum is number ',LOWEST,' energy=',PEMIN
   PRINT '(A)','keyword> Calculating average vibration/permutation factors for each Q bin'
   BAV(1:NQBINS)=0.0D0
   NMINQBIN(1:NQBINS)=0
   XPIQBIN(1:NQBINS)=0.0D0
   QHISTINT=QENERGY(2)-QENERGY(1)
   DO J1=1,NRESMIN
      J2=INT((EMIN(J1)-QENERGY(1))/QHISTINT)+1
      NMINQBIN(J2)=NMINQBIN(J2)+1
      XPIQBIN(J2)=XPIQBIN(J2)+1.0D0/(1.0D0*HORDERMIN(J1))
      BAV(J2)=BAV(J2)+EXP(-FVIBMIN(J1)/2.0D0+FVIBMIN(LOWEST)/2.0D0)/(1.0D0*HORDERMIN(J1))
   ENDDO
   DO J1=1,NQBINS
      IF (XPIQBIN(J1).GT.0) THEN
         BAV(J1)=BAV(J1)/XPIQBIN(J1)
         PRINT '(A,I6,A,I6,A,F12.4,A,G20.10)',' # entries in quench bin ',J1,' is ', &
  &                                   NMINQBIN(J1),' PI isomers=',XPIQBIN(J1), &
  &                                   ' relative harmonic prefactor=',BAV(J1)*XPIQBIN(J1)
      ENDIF
   ENDDO
ELSE
   PRINT '(A)','keyword> No min.data file located in working directory'
ENDIF
INQUIRE(FILE='B.data',EXIST=YESNO)
IF (YESNO) THEN
   OPEN(UNIT=10,FILE='B.data',STATUS='OLD')
   READ(10,*) USEQBINS
   CLOSE(10)
   RESERVOIRT=.TRUE.
   PRINT '(A,G20.10)','Final phase will use prefactors from supplied minima up to quench ',USEQBINS
   INQUIRE(FILE='min.data',EXIST=YESNO)
   IF (.NOT.YESNO) THEN
      PRINT '(A,G20.10)','ERROR *** no min.data file present'
      STOP
   ENDIF
ELSE
   PRINT '(A)','Not using any saved local minima information'
ENDIF

! IF ((PTSTART(1:2).NE.'na').AND.(BSPTFINISH(1:3).NE.'all')) CALL EQGM(NREP, NBINS, NQBINS, PTSTART)

IF (.NOT.ALLOCATED(MAXVISITS)) THEN
   ALLOCATE(MAXVISITS(NREP))
ENDIF
MAXVISITS(1:NREP)=-1.0D0
DO J1=1,NREP
   DO J2=1,NBINS
      IF (VISITS(J2,J1).GT.MAXVISITS(J1)) MAXVISITS(J1)=VISITS(J2,J1)
   ENDDO
!  PRINT '(3(A,F15.1))','Threshold for neglect of bins in replica ',J1,' is ',MAX(0.0D0,MAXVISITS(J1)*NEGLECT), &
! &                  ' visits. Maximum value=',MAXVISITS(J1)
   WRITE(*,'(A,I5,A,F10.4,A,F10.1,A,F15.1,A,F15.1)') 'T for replica ',J1,' is ',REPT(J1),' Mvisits=', &
  &                        VTOTAL(J1)*1.0D-6, &
  &                        ' neglect thresh ',MAX(0.0D0,MAXVISITS(J1)*NEGLECT),' max visits=',MAXVISITS(J1)

ENDDO
DO J1=1,NREP
   DO J2=1,NBINS
      IF ((1.0D0*VISITS(J2,J1))/(1.0D0*MAXVISITS(J1)).LT.NEGLECT) THEN
         VISITS(J2,J1)=0
      ENDIF
   ENDDO
ENDDO
DELTA=ENERGY(2)-ENERGY(1)
! IF (TMIN.LT.2.0D0*DELTA) THEN
!    PRINT '(3(A,G20.10))','WARNING - minimum temperature ',TMIN,' is < bin width ',DELTA,' - resetting Tmin to ',2.0D0*DELTA
!    TMIN=2*DELTA
! ENDIF
!
! ss2029> normalize visits histograms to make CHI2 and RMS independent of system
! and length of run 
!
DO J1=1,NREP
! sum histogram
   VTOT=0
   DO J2=1,NBINS
      VTOT = VTOT + VISITS(J2,J1)
   ENDDO
! divide
   DO J2=1,NBINS
      VISITS(J2,J1) = VISITS(J2,J1)/VTOT
   ENDDO
ENDDO

!
! Put inverse Boltzmann factor into the weight.
!
DO J1=1,NREP
   DO J2=1,NBINS
      IF (VISITS(J2,J1).GT.0) THEN
!        PRINT '(A,2I6,G20.10)','J2,J1,VISITS=',J2,J1,VISITS(J2,J1)
!        PRINT '(A,4G20.10)','ENERGY(J2),ENERGY(mb/3),ENERGY(J2)-ENERGY(NBINS/3),REPT=', &
! &                           ENERGY(J2),ENERGY(NBINS/3),ENERGY(J2)-ENERGY(NBINS/3),REPT(J1)
!        PRINT '(A,G20.10)','exp arg=',(ENERGY(J2)-ENERGY(NBINS/3))/REPT(J1)
!        WEIGHT(J2,J1)=VISITS(J2,J1)*EXP((ENERGY(J2)-ENERGY(NBINS/3))/REPT(J1))
         WEIGHT(J2,J1)=LOG(1.0D0*VISITS(J2,J1)) + (ENERGY(J2)-ENERGY(NBINS/3))/REPT(J1)
      ENDIf
   ENDDO
ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Solve the simplest lsq fit problem for direct visits to pe bins
! 
IF (ALLOCATED(WIJ)) DEALLOCATE(WIJ)
ALLOCATE(WIJ(NREP,NBINS))
DO J1=1,NBINS
   DO J2=NREPMIN,NREP
      IF (VISITS(J1,J2).GT.0) THEN
!        PRINT '(A,2I6,G20.10)','J2,J1,weight=',J2,J1,WEIGHT(J1,J2)
!        WIJ(J2,J1)=LOG(WEIGHT(J1,J2))
         WIJ(J2,J1)=WEIGHT(J1,J2)
      ELSE
         WIJ(J2,J1)=0.0D0
      ENDIF
   ENDDO
ENDDO 
ALLZERO(1:NBINS)=.TRUE.
NEVAR=0
iloop: DO J1=1,NBINS
   DO J2=NREPMIN,NREP
     IF (VISITS(J1,J2).GT.0) THEN
        ALLZERO(J1)=.FALSE.
        NEVAR=NEVAR+1
        CYCLE iloop
     ENDIF
   ENDDO
ENDDO iloop

NVAR=NREP-NREPMIN+1+NEVAR
PRINT '(A,I10)','number of bin weight variables=',NEVAR
! PRINT '(A,I8)','number of unknown constants of proportionality for replicas=',NREP-NREPMIN
PRINT '(A,I8)','number of unknown constants of proportionality for replicas=',NREP
IF (NEVAR.EQ.0) STOP
!
! Initialise variables. 
!
IF (ALLOCATED(BESTVAR)) DEALLOCATE(BESTVAR)
IF (ALLOCATED(VAR)) DEALLOCATE(VAR)
ALLOCATE(VAR(NVAR),BESTVAR(NVAR))
CALL SDPRND(0) ! initialise

DO J1=1,NVAR
   VAR(J1)=DPRAND()
ENDDO
!
! Test analytic derivatives
!
IF (.FALSE.) THEN
   ALLOCATE(GRAD(NVAR),DGRAD(NVAR))
   CALL GETCHI(NBINS,NEVAR,NREP,NREPMIN,WIJ,VAR,CHI2,GRAD,ALLZERO,VISITS,RMS,VTOTAL2)
   DIFF=0.0001D0
   DO J1=1,NVAR
      VAR(J1)=VAR(J1)+DIFF
      CALL GETCHI(NBINS,NEVAR,NREP,NREPMIN,WIJ,VAR,CHI2PLUS,DGRAD,ALLZERO,VISITS,RMS,VTOTAL2)
      VAR(J1)=VAR(J1)-2.0D0*DIFF
      CALL GETCHI(NBINS,NEVAR,NREP,NREPMIN,WIJ,VAR,CHI2MINUS,DGRAD,ALLZERO,VISITS,RMS,VTOTAL2)
      VAR(J1)=VAR(J1)+DIFF
      IF (GRAD(J1).NE.0.0D0) THEN
         PRINT '(A,I5,3G20.10)','J1,num,anal,rat=',J1,(CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF),GRAD(J1), &
  &                               (CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF*GRAD(J1))
      ELSE
         PRINT '(A,I5,3G20.10)','J1,num,anal=',J1,(CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF),GRAD(J1)
      ENDIF
   ENDDO
   STOP
ENDIF

MUPDATES=400
TOL=1.0D-7 

ITMAX=100000
NHOPS=1
BESTRMS=1.0D100
RNDSEED=1
MAXBHSTEP=1.0D0
CALL SDPRND(RNDSEED)

DO J1=1,NHOPS

   CALL MYLBFGS(NVAR,MUPDATES,VAR,TOL,ITDONE,ITMAX,CHI2,CONVERGED,NBINS,NEVAR,NREP,NREPMIN,WIJ,ALLZERO,VISITS,VTOTAL2,RMS)

   IF (RMS.LT.BESTRMS) THEN
      BESTRMS=RMS
      BESTVAR(1:NVAR)=VAR(1:NVAR)
   ENDIF

   DO J2=1,NVAR
      VAR(J2)=VAR(J2)+VAR(J2)*(2.0D0*DPRAND()-1.0D0)*MAXBHSTEP
   ENDDO

ENDDO

IF (ALLOCATED(EVAR)) DEALLOCATE(EVAR)
ALLOCATE(EVAR(NBINS))
EVAR(1:NBINS)=0.0D0
NCOUNT=0

DO J1=1,NBINS
   IF (.NOT.ALLZERO(J1)) THEN
      NCOUNT=NCOUNT+1
      VAR(NCOUNT)=EXP(BESTVAR(NCOUNT)-BESTVAR(1))
      EVAR(J1)=VAR(NCOUNT)
   ENDIF
ENDDO
PRINT '(A)','Final bin weights - dumping relative values to weights.A'
NCOUNT=0
WRITE(FNAME,'(A9)') 'weights.A'
OPEN(UNIT=1,FILE=FNAME,STATUS='UNKNOWN')
DO J1=1,NBINS
   IF (.NOT.ALLZERO(J1)) THEN
      WRITE(1,'(I8,2G20.10)') J1,ENERGY(J1),LOG(EVAR(J1))
   ENDIF
ENDDO
CLOSE(1)
! PRINT '(A)','Final normalisation terms:'
! PRINT '(I8,2G20.10)',(J1,VAR(NEVAR+J1),VAR(NEVAR+J1)-VAR(NEVAR+1),J1=NREPMIN,NREP)

TINT=(TMAX-TMIN)/(NTEMP-1)
FSTRING='Cv.out.A' 
OPEN(UNIT=1,FILE=FSTRING,STATUS='UNKNOWN')
CALL DUMPCV(NBINS,ALLZERO,ENERGY,TMIN,TINT,EVAR,DELTA,KAPPA,NTEMP)
CLOSE(1)
!
!!!!!!!!!!!! Read data from Visits2.his.<n> files !!!!!!!!!!!
!
PRINT '(A)',' '
PRINT '(A)','Now reading quench visit statistics'

VISITS2(1:NREP,1:NBINS,1:NQBINS)=0
DO J1=NREPMIN,NREP
   IF (J1.GE.10) THEN
      IF (BSPTFINISH(1:3).EQ.'all') THEN
         WRITE(FSTRING,'(I2,A12)') J1,'/Visits2.his'
      ELSE
         WRITE(FSTRING,'(I2,A)') J1,'/Visits2.his.' // TRIM(ADJUSTL(BSPTFINISH))
      ENDIF
      IF (BSSTART(1:2).NE.'na') WRITE(FSTRING2,'(I2,A)') J1,'/Visits2.his.' // TRIM(ADJUSTL(BSSTART))
   ELSE
      IF (BSPTFINISH(1:3).EQ.'all') THEN
         WRITE(FSTRING,'(I1,A12)') J1,'/Visits2.his'
      ELSE
         WRITE(FSTRING,'(I1,A)') J1,'/Visits2.his.' // TRIM(ADJUSTL(BSPTFINISH))
      ENDIF
      IF (BSSTART(1:2).NE.'na') WRITE(FSTRING2,'(I1,A)') J1,'/Visits2.his.' // TRIM(ADJUSTL(BSSTART))
   ENDIF
   INQUIRE(FILE=TRIM(ADJUSTL(FSTRING)),EXIST=YESNO)
   IF (.NOT.YESNO) THEN
      PRINT '(3A)','File ',TRIM(ADJUSTL(FSTRING)),' not found'
      STOP
   ENDIF
   IF (BSSTART(1:2).NE.'na') INQUIRE(FILE=TRIM(ADJUSTL(FSTRING2)),EXIST=YESNO)
   IF (.NOT.YESNO) THEN
      PRINT '(3A)','File ',TRIM(ADJUSTL(FSTRING2)),' not found'
      STOP
   ENDIF
   OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FSTRING)),STATUS='OLD',FORM='UNFORMATTED')
   IF (BSSTART(1:2).NE.'na') THEN
      OPEN(UNIT=2,FILE=TRIM(ADJUSTL(FSTRING2)),STATUS='OLD',FORM='UNFORMATTED')
   ENDIF
   READ(1) VISITS2(J1,1:NBINS,1:NQBINS)
   CLOSE(1)
   IF (BSSTART(1:2).NE.'na') THEN
      READ(2) VISITST(1:NBINS,1:NQBINS)
      VISITS2(J1,1:NBINS,1:NQBINS)=VISITS2(J1,1:NBINS,1:NQBINS)-VISITST(1:NBINS,1:NQBINS)
      CLOSE(2)
   ENDIF
ENDDO
MAXVISITS(1:NREP)=-1.0D0
DO J1=NREPMIN,NREP
   DO J2=1,NBINS
      DO J4=1,NQBINS ! J4 labels the bin quenched to
         IF (VISITS2(J1,J2,J4).GT.MAXVISITS(J1)) MAXVISITS(J1)=VISITS2(J1,J2,J4)
      ENDDO
   ENDDO
ENDDO
DO J1=NREPMIN,NREP
   DO J2=1,NBINS
!
! This test is needed to avoid inconsistencies when applying NEGLECT
! to VISITS and VISITS2.
!
      IF (VISITS(J2,J1).EQ.0) THEN 
         DO J4=1,NQBINS ! J4 labels the bin quenched to
            VISITS2(J1,J2,J4)=0
         ENDDO
!     ELSE
!
! using MAXVISITS here cannot be right!! DJW
!
!        DO J4=1,NQBINS ! J4 labels the bin quenched to
!           IF ((1.0D0*VISITS2(J1,J2,J4))/(1.0D0*MAXVISITS(J1)).LT.NEGLECT) THEN
!              VISITS2(J1,J2,J4)=0
!           ENDIF
!        ENDDO
      ENDIF
   ENDDO
ENDDO
!
! For each quench bin, identify the lowest PE bin visited directly via VISITS2.
!
IF (ALLOCATED(LOWESTDIRECT)) DEALLOCATE(LOWESTDIRECT)
IF (ALLOCATED(BESTBINS)) DEALLOCATE(BESTBINS)
IF (ALLOCATED(HIGHESTDIRECT)) DEALLOCATE(HIGHESTDIRECT)
ALLOCATE(LOWESTDIRECT(NQBINS),HIGHESTDIRECT(NQBINS),BESTBINS(NQBINS))
LOWESTDIRECT(1:NQBINS)=NBINS+1
HIGHESTDIRECT(1:NQBINS)=-1
DO J3=1,NQBINS
   DO J1=NREPMIN,NREP
      DO J2=1,NBINS
         IF (VISITS2(J1,J2,J3).GT.0) THEN
            IF (J2.LT.LOWESTDIRECT(J3)) LOWESTDIRECT(J3)=J2
            IF (J2.GT.HIGHESTDIRECT(J3)) HIGHESTDIRECT(J3)=J2
         ENDIF
      ENDDO 
   ENDDO 
ENDDO
!
! For each pe bin, identify the lowest and highest q bins visited
!
IF (ALLOCATED(LOWESTQUENCH)) DEALLOCATE(LOWESTQUENCH)
IF (ALLOCATED(HIGHESTQUENCH)) DEALLOCATE(HIGHESTQUENCH)
ALLOCATE(LOWESTQUENCH(NBINS),HIGHESTQUENCH(NBINS))
LOWESTQUENCH(1:NBINS)=NQBINS+1
HIGHESTQUENCH(1:NBINS)=-1
PEBINMAX=-1
PEBINMIN=NBINS+1
DO J3=1,NBINS
   DO J1=NREPMIN,NREP
      DO J2=1,NQBINS
         IF (VISITS2(J1,J3,J2).GT.0) THEN
            IF (J2.LT.LOWESTQUENCH(J3)) LOWESTQUENCH(J3)=J2
            IF (J2.GT.HIGHESTQUENCH(J3)) HIGHESTQUENCH(J3)=J2
            IF (J3.LT.PEBINMIN) PEBINMIN=J3
            IF (J3.GT.PEBINMAX) PEBINMAX=J3
         ENDIF
      ENDDO
   ENDDO
ENDDO

! PRINT '(A)','Bin #, lowest and highest PE bins visited directly and in quenches:'
! DO J1=1,NQBINS
!    IF ((HIGHESTDIRECT(J1).EQ.-1).AND.(HIGHESTQUENCH(J1).EQ.-1)) THEN
!    ELSEIF (HIGHESTDIRECT(J1).EQ.-1) THEN
!       PRINT '(I8,A,4X,2I8)',J1,"                ",LOWESTQUENCH(J1),HIGHESTQUENCH(J1)
!    ELSEIF (HIGHESTQUENCH(J1).EQ.-1) THEN
!       PRINT '(I8,4X,2I8)',J1,LOWESTDIRECT(J1),HIGHESTDIRECT(J1)
!    ELSE
!       PRINT '(I8,4X,4I8)',J1,LOWESTDIRECT(J1),HIGHESTDIRECT(J1),LOWESTQUENCH(J1),HIGHESTQUENCH(J1)
!    ENDIF
! ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Solve the second lsq fit problem for direct visits to pe bins broken down into q bins.
! This is an analytic result in the new formulation.
! 
NEVAR=0
ALLZERO2(1:NBINS,1:NQBINS)=.TRUE.
VARMAP2(1:NBINS,1:NQBINS)=1
!
!  Map visited PE, quench PE bins onto variables. 
!
DO J3=1,NBINS              ! J3 labels the PE bin we quenched from
   loop2: DO J4=1,NQBINS   ! J4 labels the PE bin we quench to
      DO J2=NREPMIN,NREP
         IF (VISITS2(J2,J3,J4).GT.0) THEN
            IF (ALLZERO2(J3,J4)) THEN
               NEVAR=NEVAR+1
               ALLZERO2(J3,J4)=.FALSE.
               VARMAP2(J3,J4)=NEVAR
               CYCLE loop2
            ENDIF
         ENDIF
      ENDDO 
   ENDDO loop2
ENDDO
IF (NEVAR.EQ.0) STOP
NVAR=NEVAR

PRINT '(A,I10)','number of pe/quench bin weights for 2D fit=',NEVAR
IF (ALLOCATED(VAR2)) DEALLOCATE(VAR2)
IF (ALLOCATED(BESTVAR)) DEALLOCATE(BESTVAR)
ALLOCATE(VAR2(NVAR),BESTVAR(NVAR))
!
! ALLZERO2(J3,J4) is true if there are no quenches to quench bin J4 from PE bin J3
!                 for any replica
! ALLZERO(J1) is true if there are no visits to PE bin J1 (and hence no quenches from it)
!
ALLZERO(1:NBINS)=.TRUE.
qloop2: DO J1=1,NBINS
   DO J2=1,NQBINS
     IF (.NOT.ALLZERO2(J1,J2)) THEN
        ALLZERO(J1)=.FALSE.
        CYCLE qloop2
     ENDIF
   ENDDO
ENDDO qloop2
!
!  Change VISITS(J1,J3) to be the J2 (quench, k, bin) sum of VISITS2(J3,J1,J2) if
!  necessary. Required when the NEGLECT criterion removes some gamma,k bins
!  from consideration.
!
PRINT '(A)','Using all visits data for PE bins in next phase'
DO J1=1,NBINS
   NDUMMY2=0
   DO J3=NREPMIN,NREP
      NDUMMY=0
      DO J2=1,NQBINS
         NDUMMY=NDUMMY+VISITS2(J3,J1,J2)
      ENDDO
      IF (NDUMMY.NE.VISITS(J1,J3)) THEN
!        PRINT '(A,2I8,A,F15.1,A,I12)','Adjusting # visits for gamma,r=',J1,J3,' from ',VISITS(J1,J3),' to ',NDUMMY
         VISITS(J1,J3)=NDUMMY
      ENDIF
      NDUMMY2=NDUMMY2+NDUMMY
   ENDDO
   IF ((NDUMMY2.EQ.0).AND.(.NOT.ALLZERO(J1))) THEN
      PRINT '(A,I8,L5)','ERROR - NDUMMY2=0 but J1,ALLZERO=',J1,ALLZERO(J1)
      STOP
   ENDIF
ENDDO

DO J1=1,NBINS
   IF (ALLZERO(J1)) CYCLE
   DUMMY=0.0D0
   DO J2=1,NQBINS
      IF (ALLZERO2(J1,J2)) CYCLE
      NCOUNT=NCOUNT+1

!     BESTVAR(NCOUNT)=0.0D0
!     MGAMMA=0
!     DO J3=NREPMIN,NREP
!        IF (VISITS(J1,J3).GT.0) THEN
!           BESTVAR(NCOUNT)=BESTVAR(NCOUNT)+(VISITS2(J3,J1,J2)*1.0D0)/VISITS(J1,J3)
!           MGAMMA=MGAMMA+1
!        ELSE
!           IF (VISITS2(J3,J1,J2).GT.0) THEN
!              PRINT '(A,I8)','ERROR - VISITS=0 but VISITS2=',VISITS2(J3,J1,J2)
!              PRINT '(A,3I8)','J3,J2,J1=',J3,J2,J1
!           ENDIF
!        ENDIF
!     ENDDO
!     IF (MGAMMA.EQ.0) THEN
!        PRINT '(A)','ERROR - MGAMMA=0'
!        STOP
!     ENDIF
!     BESTVAR(NCOUNT)=BESTVAR(NCOUNT)*EVAR(J1)/MGAMMA

      NDUMMY2=0
      NDUMMY3=0
      DO J3=NREPMIN,NREP
         NDUMMY2=NDUMMY2+VISITS2(J3,J1,J2)
         NDUMMY3=NDUMMY3+VISITS(J1,J3)
      ENDDO
      IF (NDUMMY3.EQ.0) THEN
         PRINT '(A)','ERROR - NDUMMY3=0'
         STOP
      ENDIF
      BESTVAR(NCOUNT)=EVAR(J1)*NDUMMY2/NDUMMY3
      
      DUMMY=DUMMY+BESTVAR(NCOUNT)
   ENDDO
   IF (ABS((DUMMY-EVAR(J1))*100/DUMMY).GT.1.0D-3) THEN
      PRINT '(A,I8,2G20.10)','J1,EVAR,DUMMY=',J1,EVAR(J1),DUMMY
      STOP
   ENDIF
ENDDO
!
! EVAR is reset, but it should be exactly the same for the sum over q bins
! as it was in the fit for stage one!
!
EVAR(1:NBINS)=0.0D0
NCOUNT=0
!
!  BESTVAR contains omega_c(V^I_gamma) sum_r N_gamma,k,r / sum_r N_gamma,r with an index running
!  over k bins for a given gamma bin.
!  OLDVAR saves the ln of BESTVAR
!  VAR2   saves BESTVAR relative to the value for the first non-zero bin
!
IF (ALLOCATED(OLDVAR)) DEALLOCATE(OLDVAR)
ALLOCATE(OLDVAR(NEVAR))
NVAROLD=NEVAR
DO J1=1,NBINS
   DO J4=1,NQBINS
      IF (.NOT.ALLZERO2(J1,J4)) THEN
         IF (J1.LT.LOWESTDIRECT(J4)) THEN
            PRINT '(3(A,I8))','ERROR - lowest direct pe bin for quench bin ',J4,' should be ',LOWESTDIRECT(J4),' not ',J1
            STOP
         ENDIF
         NCOUNT=NCOUNT+1
         ! WRITE(1,'(G20.10)') BESTVAR(NCOUNT)
         OLDVAR(NCOUNT)=LOG(BESTVAR(NCOUNT)) ! Save the best fit from stage two so we can compare results from stage 3.
         VAR2(NCOUNT)=BESTVAR(NCOUNT)/BESTVAR(1)
         EVAR(J1)=EVAR(J1)+VAR2(NCOUNT)
       ELSE
          ! WRITE(1,'(G20.10)') 0.0D0
       ENDIF
   ENDDO
ENDDO
!
! We now sum omega_c(V^I_gamma) N_gamma,k,r / N_gamma,r over k, the quench bin index. 
! This will give back the fitted weights from phase A, except that some PE bins, gamma,
! might be included when they were neglected in phase A for populations below
! threshold. omega_c(V^I_gamma) is probably zero for these gamma values in any case.
! The results should therefore be identical to the A phase.
!
PRINT '(A)','Final bin weights - dumping relative values to weights.B'
NCOUNT=0
WRITE(FNAME,'(A9)') 'weights.B'
OPEN(UNIT=1,FILE=FNAME,STATUS='UNKNOWN')
DO J1=1,NBINS
   IF (.NOT.ALLZERO(J1)) THEN
!     PRINT '(2E20.10)',ENERGY(J1),LOG(EVAR(J1))
      WRITE(1,'(I8,2G20.10)') J1,ENERGY(J1),LOG(EVAR(J1))
   ENDIF
ENDDO
CLOSE(1)

FSTRING='Cv.out.B' 
OPEN(UNIT=1,FILE=FSTRING,STATUS='UNKNOWN')
CALL DUMPCV(NBINS,ALLZERO,ENERGY,TMIN,TINT,EVAR,DELTA,KAPPA,NTEMP)
CLOSE(1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Extrapolate the direct visits results for quench bins by fitting an anharmonic
!  functional form that has the right harmonic limit.
!
PRINT '(A)',' '
PRINT '(A)','Fitting anharmonic functional forms for individual quench bins'
!
! ALLZEROQ(J3,J2,J1) is true if we never quench to quench bin J1
! from any PE bin J2 in any replica J3 from the minimum NREPMIN to NREP.
!
ALLZEROQ(1:NQBINS)=.TRUE.
qloop: DO J1=1,NQBINS
   DO J2=1,NBINS
      DO J3=NREPMIN,NREP
         IF (VISITS2(J3,J2,J1).GT.0) THEN
            ALLZEROQ(J1)=.FALSE.
            CYCLE qloop
         ENDIF
      ENDDO
   ENDDO
ENDDO qloop
!
! SUMVISITS(J2,J1) contains the total number of quenches to quench bin J1
! from PE bin J2 over all replicas.
! MOSTVISITED(J1) is the most common PE bin that quenches to quench bin J1.
!
SUMVISITS(1:NBINS,1:NQBINS)=0
DO J1=1,NQBINS
   NDUMMY=0
   DO J2=1,NBINS
      DO J3=NREPMIN,NREP
         SUMVISITS(J2,J1)=SUMVISITS(J2,J1)+VISITS2(J3,J2,J1)
      ENDDO
      IF (SUMVISITS(J2,J1).GT.NDUMMY) THEN
         NDUMMY=SUMVISITS(J2,J1)
         MOSTVISITED(J1)=J2
      ENDIF
   ENDDO
ENDDO 

KFAC=(KAPPA*1.0D0)/2.0D0-1.0D0
MUPDATES=400
TOL=2.0D-6
ITMAX=10000
NVAR=2
FITA(1:NQBINS)=0.0D0
FITB(1:NQBINS)=0.0D0
FITC(1:NQBINS)=1.0D0 
XFIT(1:NVAR)=0.0D0
!
! Now we fit a two- or three-parameter expression for the anharmonic
! density of states for each quench bin using data over an increasingly
! large number of PE bins. The best fit is recorded.
! Fitted parameters for quench bin J1 are in FITA(J1), FITB(J1), FITC(J1).
! VARMAP2(J3,J4) is the 1D index of the number of quenches from PE bin J3
! to quench bin J4. These are really data points, not variables.
! LOWESTDIRECT(J1) is the lowest PE bin that quenches to quench bin J1.
! E1 is the PE bin energy for LOWESTDIRECT(J1). 
! QE is the energy of quench bin J1.
! NCOUNT is initialised to the position of the first possible 2D visits data point.
! MINFIT is the minimum number of data points, corresponding to bins that
!        quench to quench bin J1.
!
! If direct visits go down to the quench bin itself then we skip the fit,
! since the direct visits to PE bins cover this PE range. We have to account for
! this later when calculating the final density of states from the direct and
! fitted anharmonic forms.
!
DO J1=1,NQBINS
   IF (ALLZEROQ(J1)) CYCLE
   CALL FITABC(J1,VARMAP2,LOWESTDIRECT,OLDVAR,ALLZEROQ,HIGHESTDIRECT,MINFIT,ALLZERO2,SUMVISITS,NBINS,NQBINS,KAPPA,NVAR, &
   &              MUPDATES,ITMAX,BESTBINS,FITA,FITB,FITC,NEVAR,SIGMAX,MOSTVISITED,QENERGY,ENERGY,FIXB)
ENDDO
!
! Final densities of states are calculated from direct pe, q bin visits if
! available, and the values inferred from the fit otherwise.
!
! Fitting range for q bin J1 includes gamma bins from
! LOWESTDIRECT(J1) to NBINS.
! There is no contribution to the fit if q bin j1 is never visited for any gamma.
!
EVAR(1:NBINS)=0.0D0
DO J2=1,NQBINS            ! J2 labels the PE bin we quench to
    IF (ALLZEROQ(J2)) THEN ! no fit for this q bin - use direct visits only
       DO J1=1,NBINS ! J1 labels the PE bin we quenched from
          IF (ALLZERO2(J1,J2)) CYCLE
          NCOUNT=VARMAP2(J1,J2)
          EVAR(J1)=EVAR(J1)+EXP(OLDVAR(NCOUNT))
       ENDDO
    ELSE ! use a combination of extrapolated or direct visit results
       DO J1=1,LOWESTDIRECT(J2)-1 ! interpolated values only available
          IF (ENERGY(J1)-QENERGY(J2).GT.0.0D0) THEN
             EVAR(J1)=EVAR(J1)+EXP((KFAC+(ENERGY(J1)-QENERGY(J2))*EXP(FITA(J2))) &
  &                                           *LOG(ENERGY(J1)-QENERGY(J2))+FITB(J2))
          ENDIF
       ENDDO
!      DO J1=LOWESTDIRECT(J2),(LOWESTDIRECT(J2)+BESTBINS(J2))/2
       DO J1=LOWESTDIRECT(J2),BESTBINS(J2) ! range used for fit
          FRACTION=(J1*1.0D0-LOWESTDIRECT(J2)*1.0D0)/(BESTBINS(J2)*1.0D0-LOWESTDIRECT(J2)*1.0D0)
!         IF ((SUMVISITS(J1,J2).GT.0).AND.(ENERGY(J1)-QENERGY(J2).GT.0.0D0)) THEN
          IF (ENERGY(J1)-QENERGY(J2).GT.0.0D0) THEN 
             EVAR(J1)=EVAR(J1)+(1.0D0-FRACTION)*EXP((KFAC+(ENERGY(J1)-QENERGY(J2))*EXP(FITA(J2))) &
  &                                           *LOG(ENERGY(J1)-QENERGY(J2))+FITB(J2))
          ENDIF
          NCOUNT=VARMAP2(J1,J2)
          EVAR(J1)=EVAR(J1)+FRACTION*EXP(OLDVAR(NCOUNT))
       ENDDO
!      DO J1=(LOWESTDIRECT(J2)+BESTBINS(J2))/2+1,NBINS ! use direct visits above range of fit if available
!      DO J1=LOWESTDIRECT(J2),NBINS ! use only direct visits above lowest direct
       DO J1=BESTBINS(J2)+1,NBINS ! use direct visits above range of fit if available
          IF (ALLZERO2(J1,J2)) CYCLE
          NCOUNT=VARMAP2(J1,J2)
          EVAR(J1)=EVAR(J1)+EXP(OLDVAR(NCOUNT))
       ENDDO
    ENDIF
ENDDO

PRINT '(A)',' '
PRINT '(A)','dumping pe bin, direct visits, bin energy, ln omega to weights.C'
NCOUNT=0
IF (ALLOCATED(VSUM)) DEALLOCATE(VSUM)
ALLOCATE(VSUM(NBINS))
VSUM=SUM(VISITS,DIM=2)
OPEN(UNIT=1,FILE='weights.C',STATUS='UNKNOWN')
DO J1=1,NBINS
   IF (EVAR(J1).NE.0.0D0) THEN
      DUMMY=LOG(EVAR(J1))
      EXIT
   ENDIF
ENDDO
DO J1=1,NBINS
   IF (EVAR(J1).EQ.0.0D0) CYCLE
   WRITE(1,'(I10,F15.1,4G20.10)') J1,VSUM(J1),ENERGY(J1),LOG(EVAR(J1)),LOG(EVAR(J1))-DUMMY
ENDDO
CLOSE(1)

FSTRING='Cv.out.C' 
OPEN(UNIT=1,FILE=FSTRING,STATUS='UNKNOWN')
CALL DUMPCV(NBINS,ALLZERO,ENERGY,TMIN,TINT,EVAR,DELTA,KAPPA,NTEMP)
CLOSE(1)
IF (.NOT.RESERVOIRT) STOP
!
! Next section is for the final fit. We replace the FITB values for
! quench bins below the cutoff USEQBINS by the harmonic values calculated
! from data in min.data. There is a single scaling factor, with an
! analytical value BSCALE obtained from minimising the squared difference
! between the harmonic values and the previously fitted FITB results.
!
HSUM=0.0D0
BSUM=0.0D0
DO J2=1,NQBINS
   QE=QENERGY(J2)
   IF (QE.LT.USEQBINS) THEN
      IF ((NMINQBIN(J2).NE.0).AND.(.NOT.ALLZEROQ(J2))) THEN
         HSUM=HSUM+(BAV(J2)*XPIQBIN(J2))**2
         BSUM=BSUM+EXP(FITB(J2))*BAV(J2)*XPIQBIN(J2)
      ENDIF
   ENDIF
ENDDO
BSCALE=BSUM/HSUM
PRINT '(A)',' '
PRINT '(A,G20.10)','Final phase uses harmonic vibrational prefactors for quench bins up to energy ',USEQBINS
PRINT '(A,G20.10)','Rescaling harmonic prefactors by a factor of ',BSCALE
FIXB=.TRUE.
DO J2=1,NQBINS
   QE=QENERGY(J2)
   IF (QE.GT.USEQBINS) EXIT
   IF (NMINQBIN(J2).NE.0) THEN
      IF (.NOT.ALLZEROQ(J2)) THEN
         PRINT '(A,I6,A,3G20.10)','Old and new values and change for B parameter     ',J2,' are ', &
  &                             FITB(J2),LOG(BSCALE*BAV(J2)*XPIQBIN(J2)),LOG(BSCALE*BAV(J2)*XPIQBIN(J2))-FITB(J2)
         FITB(J2)=LOG(BSCALE*BAV(J2)*XPIQBIN(J2))

!        FITA(J2)=-300.0D0 ! this is exponentiated. Large negative value means harmonic. Setting it here doesn't work, though.
!
! Should we refit the A parameter? No - this completely messes up LJ75.
!
!        DUMMY=FITA(J2)
!        CALL FITABC(J2,VARMAP2,LOWESTDIRECT,OLDVAR,ALLZEROQ,HIGHESTDIRECT,MINFIT,ALLZERO2,SUMVISITS,NBINS,NQBINS,KAPPA,NVAR, &
!  &              MUPDATES,ITMAX,BESTBINS,FITA,FITB,FITC,NEVAR,SIGMAX,MOSTVISITED,QENERGY,ENERGY,FIXB)
!        PRINT '(A,I6,A,3G20.10)','Old and new values and change for A parameter     ',J2,' are ', &
! &                             DUMMY,FITA(J2),FITA(J2)-DUMMY
      ELSE
         PRINT '(A,I6,A,G20.10)','No quench visits for this bin; setting B parameter ',J2,' to ', &
  &                             LOG(BSCALE*BAV(J2)*XPIQBIN(J2))
         FITB(J2)=LOG(BSCALE*BAV(J2)*XPIQBIN(J2))
         FITA(J2)=-300.0D0 ! this is exponentiated. Large negative value means harmonic. 
!
! There were no direct visits from quenches. Use above fitted values throughout the
! q bin range.
!
         BESTBINS(J2)=NQBINS+1 
         LOWESTDIRECT(J2)=BESTBINS(J2)+1 
      ENDIF
!
! There could be direct visits/fit data for quench bins in the reservoir range where no
! minima should exist. Set any such contributions to zero.
!
! Why? There might be minima missing from min.data that were found in the sampling.
!
!  ELSE
!     ALLZEROQ(J2)=.TRUE.
   ENDIF
ENDDO
!
! Calculate harmonic superposition thermodynamics for minima in min.data
! up to the cutoff energy in B.data
!
EVAR(1:NBINS)=0.0D0
DO J2=1,NQBINS            ! J2 labels the PE bin we quench to
   QE=QENERGY(J2)
   IF ((QE.LT.USEQBINS).AND.(NMINQBIN(J2).NE.0)) THEN
      DO J1=1,NBINS ! J1 labels the PE bin we quenched from
         IF (ENERGY(J1)-QENERGY(J2).GT.0.0D0) THEN
            EVAR(J1)=EVAR(J1)+EXP(KFAC*LOG(ENERGY(J1)-QENERGY(J2))+FITB(J2))
         ENDIF
      ENDDO
   ENDIF
ENDDO

PRINT '(A)','Dumping HSA results for minima selected by B.data to Cv.out.HSA'
FSTRING='Cv.out.HSA' 
OPEN(UNIT=1,FILE=FSTRING,STATUS='UNKNOWN')
CALL DUMPCV(NBINS,ALLZERO,ENERGY,TMIN,TINT,EVAR,DELTA,KAPPA,NTEMP)
CLOSE(1)
!
! Repeat previous construction with some refined FITB values using min.data.
!
! Final densities of states are calculated from direct pe, q bin visits if
! available, and the values inferred from the fit otherwise.
!
! Fitting range for q bin J1 includes gamma bins from
! LOWESTDIRECT(J1) to NBINS.
! There is no contribution to the fit if q bin j1 is never visited for any gamma.
!
EVAR(1:NBINS)=0.0D0
DO J2=1,NQBINS             ! J2 labels the PE bin we quench to
!
! Next block was causing arftefacts for LJ31, but seems to be needed for LJ75.
! Excluding contributions from the rescaled fit if either
! (1) we are above BESTBINS and there were no direct visits
! or (2) the contribution from the fit is excessively large compared to the
!        value from direct visits
! seems to work.
!
   IF ((QENERGY(J2).LT.USEQBINS).AND.(NMINQBIN(J2).NE.0)) THEN ! don't try to use direct visits even if we have some
      DO J1=1,NBINS ! J1 labels the PE bin
         IF (ENERGY(J1)-QENERGY(J2).LE.0.0D0) CYCLE
         IF (J1.LT.BESTBINS(J2)) THEN
            DUMMY=1.0D0
         ELSEIF (.NOT.ALLZERO2(J1,J2)) THEN
            NCOUNT=VARMAP2(J1,J2)
            DUMMY=EXP((KFAC+(ENERGY(J1)-QENERGY(J2))*EXP(FITA(J2)))*LOG(ENERGY(J1)-QENERGY(J2))+FITB(J2)-OLDVAR(NCOUNT))
         ELSE ! we are above BESTBINS and there were no direct visits
            DUMMY=HUGE(1.0D0)
         ENDIF
         IF (DUMMY.LT.RATMAX) EVAR(J1)=EVAR(J1)+EXP((KFAC+(ENERGY(J1)-QENERGY(J2))*EXP(FITA(J2))) &
   &                              *LOG(ENERGY(J1)-QENERGY(J2))+FITB(J2))
      ENDDO
   ELSEIF (ALLZEROQ(J2)) THEN ! no quenches or fit for this q bin - use direct visits only
!  IF (ALLZEROQ(J2)) THEN ! no quenches or fit for this q bin - use direct visits only
      DO J1=1,NBINS ! J1 labels the PE bin we quenched from
         IF (ALLZERO2(J1,J2)) CYCLE
         NCOUNT=VARMAP2(J1,J2)
         EVAR(J1)=EVAR(J1)+EXP(OLDVAR(NCOUNT))
      ENDDO
   ELSE ! use a combination of extrapolated or direct visit results
      DO J1=1,LOWESTDIRECT(J2)-1 ! interpolated values only available from fit to quenches from
                                 ! higher energy PE bins
         IF (ENERGY(J1)-QENERGY(J2).GT.0.0D0) THEN
            EVAR(J1)=EVAR(J1)+EXP((KFAC+(ENERGY(J1)-QENERGY(J2))*EXP(FITA(J2))) &
  &                                          *LOG(ENERGY(J1)-QENERGY(J2))+FITB(J2))
         ENDIF
      ENDDO
      DO J1=LOWESTDIRECT(J2),BESTBINS(J2) ! range used for fit
!        FRACTION=( (J1*1.0D0-LOWESTDIRECT(J2)*1.0D0)/(BESTBINS(J2)*1.0D0-LOWESTDIRECT(J2)*1.0D0) )**60
         FRACTION=( (J1*1.0D0-LOWESTDIRECT(J2)*1.0D0)/(BESTBINS(J2)*1.0D0-LOWESTDIRECT(J2)*1.0D0) )
!        FRACTION=FRACTION**2*(3.0D0-2*FRACTION)
!        FRACTION=0.0D0
!        IF ((SUMVISITS(J1,J2).GT.0).AND.(ENERGY(J1)-QENERGY(J2).GT.0.0D0)) THEN
         IF (ENERGY(J1)-QENERGY(J2).GT.0.0D0) THEN 
            EVAR(J1)=EVAR(J1)+(1.0D0-FRACTION)*EXP((KFAC+(ENERGY(J1)-QENERGY(J2))*EXP(FITA(J2))) &
  &                                           *LOG(ENERGY(J1)-QENERGY(J2))+FITB(J2))
         ENDIF
         NCOUNT=VARMAP2(J1,J2)
         EVAR(J1)=EVAR(J1)+FRACTION*EXP(OLDVAR(NCOUNT))
      ENDDO
      DO J1=BESTBINS(J2)+1,NBINS ! use direct visits above range of fit if available
         IF (ALLZERO2(J1,J2)) CYCLE
         NCOUNT=VARMAP2(J1,J2)
         EVAR(J1)=EVAR(J1)+EXP(OLDVAR(NCOUNT))
      ENDDO
   ENDIF
ENDDO

PRINT '(A)',' '
PRINT '(A)','dumping pe bin, direct visits, bin energy, ln omega to weights.D'
NCOUNT=0
IF (ALLOCATED(VSUM)) DEALLOCATE(VSUM)
ALLOCATE(VSUM(NBINS))
VSUM=SUM(VISITS,DIM=2)
OPEN(UNIT=1,FILE='weights.D',STATUS='UNKNOWN')
DO J1=1,NBINS
   IF (EVAR(J1).NE.0.0D0) THEN
      DUMMY=LOG(EVAR(J1))
      EXIT
   ENDIF
ENDDO
DO J1=1,NBINS
   IF (EVAR(J1).EQ.0.0D0) CYCLE
   WRITE(1,'(I10,F15.1,4G20.10)') J1,VSUM(J1),ENERGY(J1),LOG(EVAR(J1)),LOG(EVAR(J1))-DUMMY
ENDDO
CLOSE(1)

FSTRING='Cv.out.D' 
OPEN(UNIT=1,FILE=FSTRING,STATUS='UNKNOWN')
CALL DUMPCV(NBINS,ALLZERO,ENERGY,TMIN,TINT,EVAR,DELTA,KAPPA,NTEMP)
CLOSE(1)
!
! Now estimate relative number of minima per quench bin. We divide through
! by the mean vibrational partition function per minimum in this bin.
! Subtract the log of the value for the first non-zero estimate to get relative
! populations.
! ALLZEROQ(J2) is set true for quench bin J2 if there were not visits and
! if no fit was achievable.
!
PRINT '(A)',' '
PRINT '(A)','Estimating relative number of distinct minima (permutation-inversion isomers) and structures per quench bin'
PRINT '(A)','Constant factor P=2*product of factorials'
FIRST=.FALSE.
XDUMMY=0.0D0
XXDUMMY=0.0D0
DO J2=1,NQBINS
   IF ((NMINQBIN(J2).NE.0).AND.((.NOT.ALLZEROQ(J2)).OR.(QENERGY(J2).LT.USEQBINS))) THEN
      FITB(J2)=FITB(J2)-LOG(BAV(J2))
      IF (.NOT.FIRST) THEN
!        DUMMY=FITB(J2)-LOG(BAV(J2))
         DUMMY=FITB(J2)
         FIRSTFAC=XPIQBIN(J2)
         PRINT *,'J2,FIRST,FIRSTFAC=',J2,FIRST,FIRSTFAC
         FIRST=.TRUE.
      ENDIF
!     PRINT '(A,I4,L5,4G15.5,I4,G15.5)','J2,FIRST,FITB(J2),DUMMY,FITB(J2)-DUMMY,FIRSTFAC,NMINQBIN,XPIQBIN ', &
! &               J2,FIRST,FITB(J2),DUMMY,FITB(J2)-DUMMY,FIRSTFAC,NMINQBIN(J2),XPIQBIN(J2)
      FITB(J2)=FITB(J2)-DUMMY
!
! FIRSTFAC=1/order of point group for first min
! XPIQBIN(J2)=sum over quench bin J2 1/order of point group for minima in J2
! FITB(J2) is now ln (exp(B_J2-B_first) BAV(first)/ BAV(J2))
!
      PRINT '(A,I6,A,I6,4(A,G17.7))',' min in q bin ',J2,' : ', &
  &                                NMINQBIN(J2),' ln(PI isomers/P)',FITB(J2)+LOG(FIRSTFAC),' #/P',EXP(FITB(J2))*FIRSTFAC, &
  &                                             ' ln(structures)',FITB(J2)+LOG(FIRSTFAC*NMINQBIN(J2)/XPIQBIN(J2)), &
  &                                             ' #',EXP(FITB(J2))*FIRSTFAC*NMINQBIN(J2)/XPIQBIN(J2)
      XDUMMY=XDUMMY+EXP(FITB(J2))*FIRSTFAC
      XXDUMMY=XXDUMMY+EXP(FITB(J2))*FIRSTFAC*NMINQBIN(J2)/XPIQBIN(J2)
   ENDIF
ENDDO
FSTRING='minima.pdf'
OPEN(UNIT=1,FILE=FSTRING,STATUS='UNKNOWN')
WRITE(1,'(A)') '#  bin  quench energy    ln(PI isomers/P)  PI isomers/P   PI isomers/total ' // &
               ' ln(structures)    structures     structures/total'
DO J2=1,NQBINS
   QE=QENERGY(J2)
   IF ((NMINQBIN(J2).NE.0).AND.((.NOT.ALLZEROQ(J2)).OR.(QENERGY(J2).LT.USEQBINS))) THEN
      WRITE(1,'(I6,7G17.7)') J2,QE,FITB(J2)+LOG(FIRSTFAC),EXP(FITB(J2))*FIRSTFAC,EXP(FITB(J2))*FIRSTFAC/XDUMMY, &
  &                           FITB(J2)+LOG(FIRSTFAC*NMINQBIN(J2)/XPIQBIN(J2)), &
  &                           EXP(FITB(J2))*FIRSTFAC*NMINQBIN(J2)/XPIQBIN(J2), &
  &                           EXP(FITB(J2))*FIRSTFAC*NMINQBIN(J2)/(XPIQBIN(J2)*XXDUMMY)
   ENDIF
ENDDO
CLOSE(1)
PRINT '(A,G20.10)','Total number of distinct permutation-inversion isomers excluding factor of 2 N! is ',XDUMMY
PRINT '(A,G20.10)','Total number of structures is ',XXDUMMY

END PROGRAM CALCCV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GETCHI(NBINS,NEVAR,NREP,NREPMIN,WIJ,VAR,CHI2,GRAD,ALLZERO,VISITS,RMS,VTOTAL2)
IMPLICIT NONE
INTEGER NEVAR, NREP, NBINS, NCOUNT, J1, J2, NREPMIN
DOUBLE PRECISION GRAD(NEVAR+NREP-NREPMIN+1), VAR(NEVAR+NREP-NREPMIN+1), CHI2, RIJ, RMS, VTOTAL2
DOUBLE PRECISION WIJ(NREP,NBINS)
DOUBLE PRECISION VISITS(NBINS,NREP)
LOGICAL ALLZERO(NBINS)

CHI2=0.0D0
GRAD(1:NEVAR+NREP-NREPMIN+1)=0.0D0
NCOUNT=0
eloop: DO J1=1,NBINS
   IF (ALLZERO(J1)) CYCLE eloop
   NCOUNT=NCOUNT+1
   DO J2=NREPMIN,NREP
      RIJ=WIJ(J2,J1)-(VAR(NCOUNT)-VAR(NEVAR+J2-NREPMIN+1))
      CHI2=CHI2+VISITS(J1,J2)*RIJ**2
      GRAD(NCOUNT)=GRAD(NCOUNT)    -VISITS(J1,J2)*RIJ
      GRAD(NEVAR+J2-NREPMIN+1)=GRAD(NEVAR+J2-NREPMIN+1)+VISITS(J1,J2)*RIJ
   ENDDO
ENDDO eloop

! ss2029 > normalize CHI2 by dividing by number of variables 
CHI2=CHI2/(NEVAR+NREP-NREPMIN+1) 

GRAD(1:NEVAR+NREP-NREPMIN+1)=2*GRAD(1:NEVAR+NREP-NREPMIN+1)/(NEVAR+NREP-NREPMIN+1)
! GRAD(1)=0.0D0 ! freeze first weight
RMS=0.0D0
DO J1=1,NEVAR+NREP-NREPMIN+1
   RMS=RMS+GRAD(J1)**2
ENDDO
RMS=SQRT(RMS/(NEVAR+NREP-NREPMIN+1))

END SUBROUTINE GETCHI

SUBROUTINE GETCHI3(NVAR,NBINS,NQBINS,XFIT,CHI2,GRAD,SUMVISITS,RMS,VARMAP2,NSTART,MAXBIN,NQ,OLDVAR,QE,PE,KAPPA)
IMPLICIT NONE
INTEGER NBINS, J1, NCOUNT, NQBINS, KAPPA, NSTART, NQ, MAXBIN, NVAR
DOUBLE PRECISION GRAD(NVAR), XFIT(NVAR), CHI2, RMS, PE(NBINS), KFAC, OLDVAR(*), QE, DUMMY, DUMMY2, DUMMY3, DUMMY1
DOUBLE PRECISION SUMVISITS(NBINS,NQBINS)
INTEGER VARMAP2(NBINS,NQBINS), NDUMMY

CHI2=0.0D0
GRAD(1:NVAR)=0.0D0
KFAC=(KAPPA*1.0D0)/2.0D0-1.0D0

DUMMY1=EXP(XFIT(1))
DO J1=NSTART,MAXBIN
   NCOUNT=VARMAP2(J1,NQ)
   DUMMY2=LOG(PE(J1)-QE)
   DUMMY3=PE(J1)-QE
   NDUMMY=SUMVISITS(J1,NQ)
   DUMMY=OLDVAR(NCOUNT)-(KFAC+DUMMY3*DUMMY1)*DUMMY2-XFIT(2)
   CHI2=CHI2+NDUMMY*DUMMY**2
   DUMMY=DUMMY*NDUMMY
   GRAD(1)=GRAD(1)-DUMMY*DUMMY1*DUMMY3*DUMMY2
   GRAD(2)=GRAD(2)-DUMMY
ENDDO 

GRAD(1:NVAR)=GRAD(1:NVAR)*2.0D0
RMS=0.0D0
DO J1=1,NVAR
   RMS=RMS+GRAD(J1)**2
ENDDO

RMS=SQRT(RMS/NVAR)

END SUBROUTINE GETCHI3

SUBROUTINE MYLBFGS3(N,M,X,EPS,ITDONE,ITMAX,ENERGY,CONVERGED,NBINS,NQBINS,SUMVISITS,RMS,VARMAP2, &
  &                 LD,MAXBIN,NQ,OLDVAR,QE,EBIN,KAPPA,FIXB)
IMPLICIT NONE
INTEGER N,M,J1,ITMAX,ITDONE,NFAIL,KAPPA,NBINS,NQBINS,MAXBIN,LD,NQ
DOUBLE PRECISION X(N),G(N),DIAG(N),W(N*(2*M+1)+2*M),SLENGTH,DDOT
DOUBLE PRECISION EPS,ENERGY,ENEW,GNEW(N),RMS,OLDVAR(*),QE,EBIN(NBINS)
LOGICAL CONVERGED
DOUBLE PRECISION GNORM,STP,YS,YY,SQ,YR,BETA,OVERLAP,DOT1,DOT2,DGUESS,MAXBFGS
INTEGER ITER,POINT,ISPT,IYPT,BOUND,NPT,CP,I,INMC,IYCN,ISCN,NDECREASE
LOGICAL MFLAG, FAILED, FIXB
DOUBLE PRECISION SUMVISITS(NBINS,NQBINS)
INTEGER VARMAP2(NBINS,NQBINS)

DGUESS=1.0D0
MAXBFGS=1.0D0
ITER=0
ITDONE=0
NFAIL=0
FAILED=.FALSE.
CALL GETCHI3(N,NBINS,NQBINS,X,ENERGY,G,SUMVISITS,RMS,VARMAP2,LD,MAXBIN,NQ,OLDVAR,QE,EBIN,KAPPA)
IF (FIXB) THEN
   G(2)=0.0D0
   RMS=ABS(G(1))
ENDIF

! WRITE(*,'(A,2G20.10,A,I6,A)') 'mylbfgs3> chi^2 and RMS force=',ENERGY,RMS,' after ',ITDONE,' steps'

10    MFLAG=.FALSE.
IF (RMS.LE.EPS) THEN
   MFLAG=.TRUE.
   IF (MFLAG) THEN
!     WRITE(*,'(A,I8,A,G20.4,A,G15.4)') 'mylbfgs3> chi^2 converged in ',ITDONE, &
! &                                                ' steps. Value=',ENERGY,' RMS force=',RMS
!     IF (ITER.GT.0) WRITE(*,'(A,F20.10)') 'mylbfgs3> Diagonal inverse Hessian elements are now ',DIAG(1)
      CONVERGED=.TRUE.
      RETURN
   ENDIF
ENDIF
!! IF ((ABS(X(3)).GT.2.0D0).OR.(ABS(X(3)).LT.1.0D-5)) THEN
!IF (ABS(X(3)).GT.2.0D0) THEN
!   PRINT '(A,G20.10)','mylbfgs3> WARNING - C parameter out of range, value=',X(3)
!   RETURN
!ENDIF

IF ((ITDONE.EQ.ITMAX).OR.FAILED) THEN
!  WRITE(*,'(A,G15.7,A,G15.7)') 'mylbfgs3> **WARNING - chi^2 did not converge, value=',ENERGY,'  MS force=',RMS
!  WRITE(*,'(A,F20.10)') 'mylbfgs3> Diagonal inverse Hessian elements are now ',DIAG(1)
   CONVERGED=.FALSE.
   RETURN
ENDIF

IF (ITER.EQ.0) THEN
   IF (N.LE.0.OR.M.LE.0) THEN
      WRITE(*,240)
240   FORMAT('mylbfgs3> IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
      STOP
   ENDIF
   POINT=0
   MFLAG=.FALSE.
   DO I=1,N
      DIAG(I)=DGUESS
   ENDDO
   ISPT= N+2*M
   IYPT= ISPT+N*M
!
!  NR step for diagonal inverse Hessian
!
   DO I=1,N
      W(ISPT+I)= -G(I)*DIAG(I)
      W(I)= -G(I)*DIAG(I)
   ENDDO
   GNORM= DSQRT(DDOT(N,G,1,G,1))
!
!  Make the first guess for the step length cautious.
!
   IF (GNORM.EQ.0.0D0) THEN
      GNORM=1.0D0 ! exact zero is presumably wrong!
      PRINT '(A)','WARNING - GNORM was zero in mylbfgs2, resetting to one'
   ENDIF
   STP=MIN(GNORM,1.0D0/GNORM)
ELSE 
   BOUND=ITER
   IF (ITER.GT.M) BOUND=M
   YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
   IF (YS.EQ.0.0D0) YS=1.0D0
!
!  Update estimate of diagonal inverse Hessian elements
!  We divide by both YS and YY at different points, so
!  they had better not be zero!
!
   YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
   IF (YY.EQ.0.0D0) YY=1.0D0
   DO I=1,N
      DIAG(I)= YS/YY
!     DIAG(I)= ABS(YS/YY)
   ENDDO
!
!     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
!     "Updating quasi-Newton matrices with limited storage",
!     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
!     ---------------------------------------------------------
!
   CP= POINT
   IF (POINT.EQ.0) CP=M
   W(N+CP)= 1.0D0/YS
   DO I=1,N
      W(I)= -G(I)
   ENDDO
   CP= POINT
   DO I= 1,BOUND
      CP=CP-1
      IF (CP.EQ. -1)CP=M-1
      SQ= DDOT(N,W(ISPT+CP*N+1),1,W,1)
      INMC=N+M+CP+1
      IYCN=IYPT+CP*N
      W(INMC)= W(N+CP+1)*SQ
      CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
   ENDDO
  
   DO I=1,N
      W(I)=DIAG(I)*W(I)
   ENDDO

   DO I=1,BOUND
      YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
      BETA= W(N+CP+1)*YR
      INMC=N+M+CP+1
      BETA= W(INMC)-BETA
      ISCN=ISPT+CP*N
      CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
      CP=CP+1
      IF (CP.EQ.M) CP=0
   ENDDO
   STP=1.0D0
ENDIF
!
!  Store the new search direction
!
! W(1)=0.0D0 ! freeze first weight
IF (FIXB) W(2)=0.0D0
DO I=1,N
   W(ISPT+POINT*N+I)= W(I)
ENDDO

DOT1=SQRT(DDOT(N,G,1,G,1))
DOT2=SQRT(DDOT(N,W,1,W,1))
OVERLAP=0.0D0
IF (DOT1*DOT2.NE.0.0D0) OVERLAP=DDOT(N,G,1,W,1)/(DOT1*DOT2)

IF (OVERLAP.GT.0.0D0) THEN
!  PRINT*,'Search direction has positive projection onto gradient - reversing step'
   DO I=1,N
      W(ISPT+POINT*N+I)= -W(I)
   ENDDO
ENDIF
      
DO I=1,N
   W(I)=G(I)
ENDDO
SLENGTH=0.0D0
DO J1=1,N
   SLENGTH=SLENGTH+W(ISPT+POINT*N+J1)**2
ENDDO
SLENGTH=SQRT(SLENGTH)
IF (STP*SLENGTH.GT.MAXBFGS) STP=MAXBFGS/SLENGTH
!
!  We now have the proposed step.
!
DO J1=1,N
   X(J1)=X(J1)+STP*W(ISPT+POINT*N+J1)
ENDDO 
NDECREASE=0
20    CONTINUE
CALL GETCHI3(N,NBINS,NQBINS,X,ENEW,GNEW,SUMVISITS,RMS,VARMAP2,LD,MAXBIN,NQ,OLDVAR,QE,EBIN,KAPPA)
IF (FIXB) THEN
   GNEW(2)=0.0D0
   RMS=ABS(GNEW(1))
ENDIF

IF (ENEW.EQ.0.0D0) ENEW=1.0D-100 ! to prevent divide by zero
IF (((ENEW-ENERGY)/ABS(ENEW).LE.1.0D-10)) THEN
   ITER=ITER+1
   ITDONE=ITDONE+1
   ENERGY=ENEW
   DO J1=1,N
      G(J1)=GNEW(J1)
   ENDDO
!  WRITE(*,'(A,2G20.10,A,I6,A,G15.5)') &
! &                'mylbfgs3> chi^2 and RMS force=',ENERGY,RMS,' after ',ITDONE,' steps, step:',STP*SLENGTH
!  PRINT '(A,3G20.10)','ABC=',X(1:2)
ELSE
!
!  chi^2 increased - try again with a smaller step size?
!
   IF (NDECREASE.GT.10) THEN  
      NFAIL=NFAIL+1
      PRINT*,' LBFGS step cannot find a lower value, NFAIL=',NFAIL
      ITER=0  !  try resetting
      DO J1=1,N
         X(J1)=X(J1)-STP*W(ISPT+POINT*N+J1)
      ENDDO 
      IF (NFAIL.GT.1) FAILED=.TRUE.
      GOTO 10
   ENDIF
   DO J1=1,N
      X(J1)=X(J1)-0.9*STP*W(ISPT+POINT*N+J1)
   ENDDO 
   NDECREASE=NDECREASE+1
   STP=STP/10.0D0
!  WRITE(*,'(A,G20.10,A,G20.10,A,F15.8)') &
! &                         'mylbfgs3> chi^2 increased from ',ENERGY,' to ',ENEW,' decreasing step to ',STP*SLENGTH
   GOTO 20
ENDIF
!
!  Compute the new step and gradient change
!
NPT=POINT*N
DO I=1,N
   W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
   W(IYPT+NPT+I)= G(I)-W(I)
ENDDO

POINT=POINT+1
IF (POINT.EQ.M) POINT=0
GOTO 10

RETURN

END SUBROUTINE MYLBFGS3

SUBROUTINE MYLBFGS(N,M,X,EPS,ITDONE,ITMAX,ENERGY,CONVERGED,NBINS,NEVAR,NREP,NREPMIN,WIJ,ALLZERO,VISITS,VTOTAL2,RMS)
IMPLICIT NONE
INTEGER N,M,J1,ITMAX,ITDONE,NFAIL
DOUBLE PRECISION X(N),G(N),DIAG(N),W(N*(2*M+1)+2*M),SLENGTH,DDOT
DOUBLE PRECISION EPS,ENERGY,ENEW,GNEW(N),RMS
LOGICAL CONVERGED
DOUBLE PRECISION GNORM,STP,YS,YY,SQ,YR,BETA,OVERLAP,DOT1,DOT2,DGUESS,MAXBFGS
INTEGER ITER,POINT,ISPT,IYPT,BOUND,NPT,CP,I,INMC,IYCN,ISCN,NDECREASE
LOGICAL MFLAG, FAILED
INTEGER NBINS, NEVAR, NREP, NREPMIN
DOUBLE PRECISION VTOTAL2
DOUBLE PRECISION WIJ(NREP,NBINS)
DOUBLE PRECISION VISITS(NBINS,NREP)
LOGICAL ALLZERO(NBINS)
DOUBLE PRECISION DIFF, CHI2, CHI2PLUS, CHI2MINUS, GRAD(N), DGRAD(N)

DGUESS=1.0D0
MAXBFGS=1000.0D0
ITER=0
ITDONE=0
NFAIL=0
FAILED=.FALSE.
CALL GETCHI(NBINS,NEVAR,NREP,NREPMIN,WIJ,X,ENERGY,G,ALLZERO,VISITS,RMS,VTOTAL2)

  WRITE(*,'(A,2G20.10,A,I6,A)') 'mylbfgs> chi^2 and RMS force=',ENERGY,RMS,' after ',ITDONE,' steps'

10    MFLAG=.FALSE.
IF (RMS.LE.EPS) THEN
   MFLAG=.TRUE.
   IF (MFLAG) THEN
      WRITE(*,'(A,I8,A,G15.5,A,G15.5)') 'mylbfgs> chi^2 converged in ',ITDONE, &
  &                                                ' steps. Value=',ENERGY,' RMS force=',RMS
      IF (ITER.GT.0) WRITE(*,'(A,F20.10)') 'mylbfgs> Diagonal inverse Hessian elements are now ',DIAG(1)
      CONVERGED=.TRUE.
      RETURN
   ENDIF
ENDIF

IF ((ITDONE.EQ.ITMAX).OR.FAILED) THEN
   WRITE(*,'(A,G15.7,A,G15.7)') 'mylbfgs> **WARNING - chi^2 did not converge, value=',ENERGY,' RMS force=',RMS
   WRITE(*,'(A,F20.10)') 'mylbfgs> Diagonal inverse Hessian elements are now ',DIAG(1)
   CONVERGED=.FALSE.
   RETURN
ENDIF

IF (ITER.EQ.0) THEN
   IF (N.LE.0.OR.M.LE.0) THEN
      WRITE(*,240)
240   FORMAT('xmylbfgs> IMPROPER INPUT PARAMETERS (N OR M ARE NOT POSITIVE)')
      STOP
   ENDIF
   POINT=0
   MFLAG=.FALSE.
   DO I=1,N
      DIAG(I)=DGUESS
   ENDDO
   ISPT= N+2*M
   IYPT= ISPT+N*M
!
!  NR step for diagonal inverse Hessian
!
   DO I=1,N
      W(ISPT+I)= -G(I)*DIAG(I)
      W(I)= -G(I)*DIAG(I)
   ENDDO
   GNORM= DSQRT(DDOT(N,G,1,G,1))
!
!  Make the first guess for the step length cautious.
!
   IF (GNORM.EQ.0.0D0) THEN
      GNORM=1.0D0 ! exact zero is presumably wrong!
      PRINT '(A)','WARNING - GNORM was zero in xmylbfgs, resetting to one'
   ENDIF
   STP=MIN(GNORM,1.0D0/GNORM)
ELSE 
   BOUND=ITER
   IF (ITER.GT.M) BOUND=M
   YS= DDOT(N,W(IYPT+NPT+1),1,W(ISPT+NPT+1),1)
   IF (YS.EQ.0.0D0) YS=1.0D0
!
!  Update estimate of diagonal inverse Hessian elements
!  We divide by both YS and YY at different points, so
!  they had better not be zero!
!
   YY= DDOT(N,W(IYPT+NPT+1),1,W(IYPT+NPT+1),1)
   IF (YY.EQ.0.0D0) YY=1.0D0
   DO I=1,N
      DIAG(I)= YS/YY
!     DIAG(I)= ABS(YS/YY)
   ENDDO
!
!     COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980,
!     "Updating quasi-Newton matrices with limited storage",
!     Mathematics of Computation, Vol.24, No.151, pp. 773-782.
!     ---------------------------------------------------------
!
   CP= POINT
   IF (POINT.EQ.0) CP=M
   W(N+CP)= 1.0D0/YS
   DO I=1,N
      W(I)= -G(I)
   ENDDO
   CP= POINT
   DO I= 1,BOUND
      CP=CP-1
      IF (CP.EQ. -1)CP=M-1
      SQ= DDOT(N,W(ISPT+CP*N+1),1,W,1)
      INMC=N+M+CP+1
      IYCN=IYPT+CP*N
      W(INMC)= W(N+CP+1)*SQ
      CALL DAXPY(N,-W(INMC),W(IYCN+1),1,W,1)
   ENDDO
  
   DO I=1,N
      W(I)=DIAG(I)*W(I)
   ENDDO

   DO I=1,BOUND
      YR= DDOT(N,W(IYPT+CP*N+1),1,W,1)
      BETA= W(N+CP+1)*YR
      INMC=N+M+CP+1
      BETA= W(INMC)-BETA
      ISCN=ISPT+CP*N
      CALL DAXPY(N,BETA,W(ISCN+1),1,W,1)
      CP=CP+1
      IF (CP.EQ.M) CP=0
   ENDDO
   STP=1.0D0
ENDIF
!
!  Store the new search direction
!
! W(1)=0.0D0 ! freeze first weight
DO I=1,N
   W(ISPT+POINT*N+I)= W(I)
ENDDO

DOT1=SQRT(DDOT(N,G,1,G,1))
DOT2=SQRT(DDOT(N,W,1,W,1))
OVERLAP=0.0D0
IF (DOT1*DOT2.NE.0.0D0) OVERLAP=DDOT(N,G,1,W,1)/(DOT1*DOT2)

IF (OVERLAP.GT.0.0D0) THEN
   PRINT*,'Search direction has positive projection onto gradient - reversing step'
   DO I=1,N
      W(ISPT+POINT*N+I)= -W(I)
   ENDDO
ENDIF
      
DO I=1,N
   W(I)=G(I)
ENDDO
SLENGTH=0.0D0
DO J1=1,N
   SLENGTH=SLENGTH+W(ISPT+POINT*N+J1)**2
ENDDO
SLENGTH=SQRT(SLENGTH)
IF (STP*SLENGTH.GT.MAXBFGS) STP=MAXBFGS/SLENGTH
!
!  We now have the proposed step.
!
DO J1=1,N
   X(J1)=X(J1)+STP*W(ISPT+POINT*N+J1)
ENDDO 
NDECREASE=0
20    CONTINUE
CALL GETCHI(NBINS,NEVAR,NREP,NREPMIN,WIJ,X,ENEW,GNEW,ALLZERO,VISITS,RMS,VTOTAL2)

IF (ENEW.EQ.0.0D0) ENEW=1.0D-100 ! to prevent divide by zero
IF (((ENEW-ENERGY)/ABS(ENEW).LE.1.0D-6)) THEN
   ITER=ITER+1
   ITDONE=ITDONE+1
   ENERGY=ENEW
   DO J1=1,N
      G(J1)=GNEW(J1)
   ENDDO
   WRITE(*,'(A,2G20.10,A,I6,A,F13.6)') &
  &                'mylbfgs> chi^2 and RMS force=',ENERGY,RMS,' after ',ITDONE,' steps, step:',STP*SLENGTH
ELSE
!
!  chi^2 increased - try again with a smaller step size?
!
   IF (NDECREASE.GT.2) THEN  
      NFAIL=NFAIL+1
      PRINT*,' LBFGS step cannot find a lower value, NFAIL=',NFAIL

   CALL GETCHI(NBINS,NEVAR,NREP,NREPMIN,WIJ,X,CHI2,GRAD,ALLZERO,VISITS,RMS,VTOTAL2)
   DIFF=0.0001D0
   DO J1=1,N
      X(J1)=X(J1)+DIFF
      CALL GETCHI(NBINS,NEVAR,NREP,NREPMIN,WIJ,X,CHI2PLUS,DGRAD,ALLZERO,VISITS,RMS,VTOTAL2)
      X(J1)=X(J1)-2.0D0*DIFF
      CALL GETCHI(NBINS,NEVAR,NREP,NREPMIN,WIJ,X,CHI2MINUS,DGRAD,ALLZERO,VISITS,RMS,VTOTAL2)
      X(J1)=X(J1)+DIFF
      IF (GRAD(J1).NE.0.0D0) THEN
         PRINT '(A,I5,3G20.10)','J1,num,anal,rat=',J1,(CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF),GRAD(J1), &
  &                               (CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF*GRAD(J1))
      ELSE
         PRINT '(A,I5,3G20.10)','J1,num,anal=',J1,(CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF),GRAD(J1)
      ENDIF
   ENDDO
   STOP

      ITER=0  !  try resetting
      DO J1=1,N
         X(J1)=X(J1)-STP*W(ISPT+POINT*N+J1)
      ENDDO 
      IF (NFAIL.GT.1) FAILED=.TRUE.
      GOTO 10
   ENDIF
   DO J1=1,N
      X(J1)=X(J1)-0.9*STP*W(ISPT+POINT*N+J1)
   ENDDO 
   NDECREASE=NDECREASE+1
   STP=STP/10.0D0
!  WRITE(*,'(A,G20.10,A,G20.10,A,F15.8)') &
! &                         'mylbfgs> chi^2 increased from ',ENERGY,' to ',ENEW,' decreasing step to ',STP*SLENGTH
   GOTO 20
ENDIF
!
!  Compute the new step and gradient change
!
NPT=POINT*N
DO I=1,N
   W(ISPT+NPT+I)= STP*W(ISPT+NPT+I)
   W(IYPT+NPT+I)= G(I)-W(I)
ENDDO

POINT=POINT+1
IF (POINT.EQ.M) POINT=0
GOTO 10

RETURN

END SUBROUTINE MYLBFGS

      double precision function ddot(n,dx,incx,dy,incy)
!
!     forms the dot product of two vectors.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dy(*),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
!
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) + dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end

      subroutine daxpy(n,da,dx,incx,dy,incy)
!
!     constant times a vector plus a vector.
!     uses unrolled loops for increments equal to one.
!     jack dongarra, linpack, 3/11/78.
!     modified 12/3/93, array(1) declarations changed to array(*)
!
      double precision dx(*),dy(*),da
      integer i,incx,incy,ix,iy,m,mp1,n
!
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
!
!        code for unequal increments or equal increments
!          not equal to 1
!
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
!
!        code for both increments equal to 1
!
!
!        clean-up loop
!
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end

!   Copyright (C) 1992  N.M. Maclaren
!   Copyright (C) 1992  The University of Cambridge

!   This software may be reproduced and used freely, provided that all
!   users of it agree that the copyright holders are not liable for any
!   damage or injury caused by use of this software and that this
!   condition is passed onto all subsequent recipients of the software,
!   whether modified or not.



        SUBROUTINE SDPRND (ISEED)
        DOUBLE PRECISION XMOD, YMOD, POLY(101), OTHER, OFFSET, X
        PARAMETER (XMOD = 1000009711.0D0, YMOD = 33554432.0D0)
        INTEGER ISEED, INDEX, IX, IY, IZ, I
        LOGICAL INITAL
        SAVE INITAL
        COMMON /RANDDP/ POLY, OTHER, OFFSET, INDEX
        DATA INITAL/.TRUE./
!
!   ISEED should be set to an integer between 0 and 9999 inclusive;
!   a value of 0 will initialise the generator only if it has not
!   already been done.
!
        IF (INITAL .OR. ISEED .NE. 0) THEN
            INITAL = .FALSE.
        ELSE
            RETURN
        END IF
!
!   INDEX must be initialised to an integer between 1 and 101
!   inclusive, POLY(1...N) to integers between 0 and 1000009710
!   inclusive (not all 0), and OTHER to a non-negative proper fraction
!   with denominator 33554432.  It uses the Wichmann-Hill generator to
!   do this.
!
        IX = MOD(ABS(ISEED),10000)+1
        IY = 2*IX+1
        IZ = 3*IX+1
        DO 10 I = -10,101
            IF (I .GE. 1) POLY(I) = AINT(XMOD*X)
            IX = MOD(171*IX,30269)
            IY = MOD(172*IY,30307)
            IZ = MOD(170*IZ,30323)
            X = MOD(DBLE(IX)/30269.0D0+DBLE(IY)/30307.0D0+DBLE(IZ)/30323.0D0,1.0D0)
  10    CONTINUE
        OTHER = AINT(YMOD*X)/YMOD
        OFFSET = 1.0D0/YMOD
        INDEX = 1
        END

        DOUBLE PRECISION FUNCTION DPRAND()
        DOUBLE PRECISION XMOD, YMOD, XMOD2, XMOD4, TINY, POLY(101), OTHER, OFFSET, X, Y
        PARAMETER (XMOD = 1000009711.0D0, YMOD = 33554432.0D0, XMOD2 = 2000019422.0D0, XMOD4 = 4000038844.0D0, TINY = 1.0D-17)
        INTEGER INDEX, N
        LOGICAL INITAL
        SAVE INITAL
        COMMON /RANDDP/ POLY, OTHER, OFFSET, INDEX
        DATA INITAL/.TRUE./
!
!   This returns a uniform (0,1) random number, with extremely good
!   uniformity properties.  It assumes that double precision provides
!   at least 33 bits of accuracy, and uses a power of two base.
!
        IF (INITAL) THEN
            CALL SDPRND (0)
            INITAL = .FALSE.
        END IF
!
!   See [Knuth] for why this implements the algorithm described in
!   the paper.  Note that this code is tuned for machines with fast
!   double precision, but slow multiply and divide; many, many other
!   options are possible.
!
        N = INDEX-64
        IF (N .LE. 0) N = N+101
        X = POLY(INDEX)+POLY(INDEX)
        X = XMOD4-POLY(N)-POLY(N)-X-X-POLY(INDEX)
        IF (X .LT. 0.0D0) THEN
            IF (X .LT. -XMOD) X = X+XMOD2
            IF (X .LT. 0.0D0) X = X+XMOD
        ELSE
            IF (X .GE. XMOD2) THEN
                X = X-XMOD2
                IF (X .GE. XMOD) X = X-XMOD
            END IF
            IF (X .GE. XMOD) X = X-XMOD
        END IF
        POLY(INDEX) = X
        INDEX = INDEX+1
        IF (INDEX .GT. 101) INDEX = INDEX-101
!
!   Add in the second generator modulo 1, and force to be non-zero.
!   The restricted ranges largely cancel themselves out.
!
   10   Y = 37.0D0*OTHER+OFFSET
        OTHER = Y-AINT(Y)
        IF (OTHER .EQ. 0.0D0) GO TO 10
        X = X/XMOD+OTHER
        IF (X .GE. 1.0D0) X = X-1.0D0
        DPRAND = X+TINY
        END

!
!  Check whether the number of quenches to the lowest q bin with non-zero visits
!  scales linearly throughout the simulation 
!
SUBROUTINE EQGM(NREP, NBINS, NQBINS, PTSTART)
IMPLICIT NONE
INTEGER NREP, NBINS, NQBINS, NINC, NAVAIL, NDUMMY, J1, NGM, J2, J3, NGM2
CHARACTER(LEN=80) PTSTART, FSTRING2, DSTRING, DUMMYSTRING
DOUBLE PRECISION, ALLOCATABLE :: QENERGY(:)
DOUBLE PRECISION NORM, NORM2
INTEGER, ALLOCATABLE :: QVISITS(:,:,:)
LOGICAL YESNO

READ(PTSTART,*) NINC

!
! How many intermediate Visits.his.<n> files do we have?
!
NAVAIL=0
DO
   NDUMMY=NINC*(NAVAIL+1)
!  NDUMMY=NINC+NAVAIL*10000
   WRITE(DSTRING,'(I9)') NDUMMY 
   DO J1=1,NREP
      IF (J1.GE.10) THEN
         WRITE(FSTRING2,'(I2,A12,A)') J1,'/Visits.his.' // TRIM(ADJUSTL(DSTRING))
      ELSE
         WRITE(FSTRING2,'(I1,A12,A)') J1,'/Visits.his.' // TRIM(ADJUSTL(DSTRING))
      ENDIF
      INQUIRE(FILE=TRIM(ADJUSTL(FSTRING2)),EXIST=YESNO)
      IF (.NOT.YESNO) GOTO 1
   ENDDO
   NAVAIL=NAVAIL+1
ENDDO

1 CONTINUE
PRINT '(A,I8,A,I8)','eqgm> checking number of quenches to lowest q bin with increment ',NINC,' up to ',NINC*NAVAIL
ALLOCATE(QENERGY(NQBINS))
ALLOCATE(QVISITS(NQBINS,NREP,NAVAIL))

DO J2=1,NAVAIL
   NDUMMY=NINC*J2
!  NDUMMY=NINC+10000*(J2-1)
   WRITE(DSTRING,'(I9)') NDUMMY 
   DO J1=1,NREP
      IF (J1.GE.10) THEN
         WRITE(FSTRING2,'(I2,A12,A)') J1,'/Visits.his.' // TRIM(ADJUSTL(DSTRING))
      ELSE
         WRITE(FSTRING2,'(I1,A12,A)') J1,'/Visits.his.' // TRIM(ADJUSTL(DSTRING))
      ENDIF
      OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FSTRING2)),STATUS='OLD')
      DO J3=1,NBINS+3
         READ(1,*) DUMMYSTRING
      ENDDO
      DO J3=1,NQBINS
         READ(1,*) QENERGY(J3),QVISITS(J3,J1,J2)
      ENDDO
      CLOSE(1)
   ENDDO
ENDDO

NGM=-1
DO J3=1,NQBINS
   DO J2=1,NAVAIL
      DO J1=1,NREP
         IF (QVISITS(J3,J1,J2).GT.0) THEN
            NGM=J3
            GOTO 2
         ENDIF
      ENDDO
   ENDDO
ENDDO

IF (NGM.LT.1) THEN
   PRINT '(A)','eqgm> no quenches detected'
   RETURN
ENDIF
2 PRINT '(A,I8)','lowest energy quench bin visited is ',NGM

DO J3=NGM+1,NQBINS
   DO J2=1,NAVAIL
      DO J1=1,NREP
         IF (QVISITS(J3,J1,J2).GT.0) THEN
            NGM2=J3
            GOTO 3
         ENDIF
      ENDDO
   ENDDO
ENDDO

3 PRINT '(A,I8)','next lowest energy quench bin visited is ',NGM2

DO J1=1,NREP
   IF (J1.GE.10) THEN
      WRITE(FSTRING2,'(I2,A8)') J1,'/gm.data'
   ELSE
      WRITE(FSTRING2,'(I1,A8)') J1,'/gm.data'
   ENDIF
   OPEN(UNIT=1,FILE=TRIM(ADJUSTL(FSTRING2)),STATUS='UNKNOWN')
   NORM=(QVISITS(NGM,J1,NAVAIL)-QVISITS(NGM,J1,1))*1.0D0
   NORM2=(QVISITS(NGM2,J1,NAVAIL)-QVISITS(NGM2,J1,1))*1.0D0
   IF (NORM.EQ.0.0D0) NORM=1.0D0
   IF (NORM2.EQ.0.0D0) NORM2=1.0D0
   DO J2=1,NAVAIL ! labels available Visits files
      WRITE(1,'(2I8,G20.10,I8,2G20.10)') J2,QVISITS(NGM, J1,J2),1.0D0*(QVISITS(NGM, J1,J2)-QVISITS(NGM, J1,1))/NORM, &
  &                                        QVISITS(NGM2,J1,J2),1.0D0*(QVISITS(NGM2,J1,J2)-QVISITS(NGM2,J1,1))/NORM2, &
  &                                        QVISITS(NGM, J1,J2)*1.0D0/MAX(1.0D0*QVISITS(NGM2,J1,J2),1.0D0)
   ENDDO
   CLOSE(1)
ENDDO

END SUBROUTINE EQGM

SUBROUTINE FITABC(J1,VARMAP2,LOWESTDIRECT,OLDVAR,ALLZEROQ,HIGHESTDIRECT,MINFIT,ALLZERO2,SUMVISITS,NBINS,NQBINS,KAPPA,NVAR, &
   &              MUPDATES,ITMAX,BESTBINS,FITA,FITB,FITC,NEVAR,SIGMAX,MOSTVISITED,QENERGY,ENERGY,FIXB)
IMPLICIT NONE
INTEGER NBINS,NQBINS,KAPPA,NVAR,MUPDATES,ITMAX,NEVAR
INTEGER J1,NCOUNT,VARMAP2(NBINS,NQBINS),LOWESTDIRECT(NQBINS),VARIABLE1,HIGHESTDIRECT(NQBINS),MINFIT
LOGICAL ALLZERO2(NBINS,NQBINS)
INTEGER VARIABLE2,NAVAIL,J2,MOSTVISITED(NBINS)
DOUBLE PRECISION SUMVISITS(NBINS,NQBINS),ENERGY(NBINS)
INTEGER MAXBIN,BININC,PREVBINS,ITDONE,NDUMMY,BESTAVAIL,BESTBINS(NQBINS)
DOUBLE PRECISION E1,QE,QENERGY(NQBINS),E2,APARAM,BPARAM,CPARAM,XFIT(3),BESTSIGMA
DOUBLE PRECISION SIGMAX,BESTA,BESTB,BESTC,TOL,CHI2,SIGMA,OLDVAR(NEVAR)
DOUBLE PRECISION FITA(NQBINS),FITB(NQBINS),FITC(NQBINS),RMS,DIFF,CHI2PLUS,CHI2MINUS
DOUBLE PRECISION, ALLOCATABLE ::  GRAD(:),DGRAD(:)
LOGICAL ALLZEROQ(NQBINS),TWOPOINTS,CONVERGED,FIXB
!
! Estimate the two fitting parameters for this quench bin J1 using the two lowest
! pe bins visited directly.
!
NCOUNT=VARMAP2(LOWESTDIRECT(J1),J1)
VARIABLE1=OLDVAR(NCOUNT)
E1=ENERGY(LOWESTDIRECT(J1))
QE=QENERGY(J1)
IF (E1.LE.QE) THEN ! skip the fit - direct visits go down to the quench bin itself
   PRINT '(A,I8,2(A,G20.10))','Skipping fit for quench bin ',J1,' lowest direct pe=',E1,' quench energy=',QE
   ALLZEROQ(J1)=.TRUE.
   RETURN
ENDIF
TWOPOINTS=.FALSE.
!  pbin: DO J2=LOWESTDIRECT(J1)+1,HIGHESTDIRECT(J1)
pbin: DO J2=LOWESTDIRECT(J1)+1,HIGHESTDIRECT(J1)-MINFIT
   IF (ALLZERO2(J2+MINFIT,J1)) CYCLE pbin
   NCOUNT=VARMAP2(J2+MINFIT,J1)
   VARIABLE2=OLDVAR(NCOUNT)
   E2=ENERGY(J2+MINFIT)
   PRINT '(3(A,I8),2(A,G15.5))','quench bin ',J1,' using pe bins ',LOWESTDIRECT(J1),' and ',J2+MINFIT, &
  &                             ' values ',E1,' and ',E2
   TWOPOINTS=.TRUE.
   EXIT pbin
ENDDO pbin
IF (.NOT.TWOPOINTS) THEN
   ALLZEROQ(J1)=.TRUE.
   RETURN
ENDIF
NAVAIL=0
DO J2=LOWESTDIRECT(J1)+1,HIGHESTDIRECT(J1)
   IF (SUMVISITS(J2,J1).GT.0) NAVAIL=NAVAIL+1
ENDDO
IF (NAVAIL.LT.MAX(MINFIT,3)) THEN
   PRINT '(A,I8)','Insufficient data points for fit: ',NAVAIL
   ALLZEROQ(J1)=.TRUE.
   PRINT '(A)',' '
   RETURN
ENDIF
!
! Initial guesses for the A and B parameters of the fitting function.
!
IF (2.*((-E1 + qe)*Log(e1 - qe) + (e2 - qe)*Log(E2 - qe)).EQ.0.0D0) THEN
   APARAM=-200.0D0
   BPARAM=1.0D0
ELSE
   APARAM= (-2*VARIABLE1 + 2*VARIABLE2 + (-2 + KAPPA)*Log(e1 - qe) - (-2 + KAPPA)*Log(E2 - qe))/ &
  &  (2.*((-E1 + qe)*Log(e1 - qe) + (e2 - qe)*Log(E2 - qe)))
   BPARAM= (2*(-e2 + qe)*VARIABLE1*Log(E2 - qe) +  &
  &   Log(e1 - qe)*(2*(E1 - qe)*VARIABLE2 + (-2 + KAPPA)*(e2 - E1)*Log(E2 - qe)))/ &
  & (2.*((E1 - qe)*Log(e1 - qe) + (-e2 + qe)*Log(E2 - qe)))
ENDIF
APARAM=ENERGY(NBINS)-ENERGY(1)
APARAM=0.0D0
CPARAM=0.0D0
IF (FIXB) THEN
   BPARAM=FITB(J1)
   APARAM=FITA(J1)
ENDIF

IF (NVAR.EQ.3) PRINT '(A,I8,A,G15.5,A,3G12.5,A,I8)','quench bin ',J1,' energy ',QE,' param guesses: ', &
  &          APARAM,BPARAM,CPARAM,' max fit bins=',NAVAIL
IF (NVAR.EQ.2) PRINT '(A,I8,A,G15.5,A,2G12.5,A,I8)','quench bin ',J1,' energy ',QE,' param guesses: ', &
  &          APARAM,BPARAM,' max fit bins=',NAVAIL
PRINT '(A,I8,A,2I8)','quench bin ',J1,' lowest and highest pe bins quenching here: ', &
  &   LOWESTDIRECT(J1),HIGHESTDIRECT(J1)
! IF (FIXB) PRINT '(2(A,G20.10))','B parameter will be fixed at ',BPARAM,' A parameter initialised as ',APARAM
   
XFIT(1)=APARAM
XFIT(2)=BPARAM
XFIT(3)=0.0D0

IF (.FALSE.) THEN
   IF (ALLOCATED(GRAD)) DEALLOCATE(GRAD,DGRAD)
   ALLOCATE(GRAD(NVAR),DGRAD(NVAR))
   MAXBIN=HIGHESTDIRECT(J1)
   CALL GETCHI3(NVAR,NBINS,NQBINS,XFIT,CHI2,GRAD,SUMVISITS,RMS,VARMAP2,LOWESTDIRECT(J1),MAXBIN,J1, &
  &             OLDVAR,QE,ENERGY,KAPPA)
   DIFF=0.0001D0
   DO J2=1,NVAR
      XFIT(J2)=XFIT(J2)+DIFF
      CALL GETCHI3(NVAR,NBINS,NQBINS,XFIT,CHI2PLUS,DGRAD,SUMVISITS,RMS,VARMAP2,LOWESTDIRECT(J1),MAXBIN, &
  &                J1,OLDVAR,QE,ENERGY,KAPPA)
      XFIT(J2)=XFIT(J2)-2.0D0*DIFF
      CALL GETCHI3(NVAR,NBINS,NQBINS,XFIT,CHI2MINUS,DGRAD,SUMVISITS,RMS,VARMAP2,LOWESTDIRECT(J1),MAXBIN, &
  &                J1,OLDVAR,QE,ENERGY,KAPPA)
      XFIT(J2)=XFIT(J2)+DIFF
      IF (GRAD(J2).NE.0.0D0) THEN
         PRINT '(A,I5,3G20.10)','J1,num,anal,rat=',J2,(CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF),GRAD(J2), &
  &                            (CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF*GRAD(J2))
      ELSE
         PRINT '(A,I5,3G20.10)','J2,num,anal=',J2,(CHI2PLUS-CHI2MINUS)/(2.0D0*DIFF),GRAD(J2)
      ENDIF
   ENDDO
!  STOP
ENDIF

BESTSIGMA=1.0D100
BESTA=APARAM; BESTB=BPARAM; BESTC=CPARAM
MAXBIN=MIN(LOWESTDIRECT(J1)+MINFIT,HIGHESTDIRECT(J1)) ! start from lowest allowed possibility
!  MAXBIN=HIGHESTDIRECT(J1) ! start from highest allowed possibility
!  MAXBIN=MOSTVISITED(J1) ! start from most visited bin
BININC=20
PREVBINS=-1
ALLZEROQ(J1)=.TRUE.
BESTBINS(J1)=-1
TOL=2.0D-6

DO 
   NAVAIL=0
   DO J2=LOWESTDIRECT(J1),MAXBIN
      IF (SUMVISITS(J2,J1).GT.0) NAVAIL=NAVAIL+1
   ENDDO
   IF (NAVAIL.LT.MINFIT) GOTO 10
   IF (NAVAIL.EQ.PREVBINS) GOTO 10
   ALLZEROQ(J1)=.FALSE.
   PREVBINS=NAVAIL
   CALL MYLBFGS3(NVAR,MUPDATES,XFIT,TOL,ITDONE,ITMAX,CHI2,CONVERGED,NBINS,NQBINS,SUMVISITS,RMS,VARMAP2, &
  &              LOWESTDIRECT(J1),MAXBIN,J1,OLDVAR,QE,ENERGY,KAPPA,FIXB)
!
!  This chi^2 contains a weight equal to the total number of visits to each gamma,k bin summed
!  over replicas. This will scale linearly with the total number of visits to all bins. 
!  Divide by total visits for q bin J1 to normalise.
!  Then divide by the number of data points - number of fitting parameters (NAVAIL-2) and take the square
!  root to give something that looks like a residual standard deviation.
!
   NDUMMY=0
   DO J2=LOWESTDIRECT(J1),MAXBIN
      NDUMMY=NDUMMY+SUMVISITS(J2,J1)
   ENDDO

   SIGMA=SQRT(CHI2/(NDUMMY*1.0D0*(NAVAIL*1.0D0-NVAR))) 

   PRINT '(I7,A,G15.5,A,3G12.5,2(A,I5),A,G12.5)',ITDONE,' cycles RMS=',RMS,' ABC:      ',XFIT(1:3), &
  &    ' for ', &
  &    NAVAIL, '/', MAXBIN-LOWESTDIRECT(J1)+1,' bins sigma=',SIGMA

!  WRITE(FNAME,'(I8)') NAVAIL
!  FNAME='fit.' // TRIM(ADJUSTL(FNAME))
!  OPEN(UNIT=1,FILE=FNAME,STATUS='UNKNOWN')
!  DO J2=LOWESTDIRECT(J1),MAXBIN
!     NCOUNT=VARMAP2(J2,J1)
!     WRITE(1,'(I8,2G20.10)') J2,OLDVAR(NCOUNT),(KFAC+EXP(XFIT(1))*(ENERGY(J2)-QE))*LOG(ENERGY(J2)-QE)+XFIT(2)
!  ENDDO
!  CLOSE(1)

   IF ((SIGMA.LT.BESTSIGMA).AND.CONVERGED) THEN
      BESTSIGMA=SIGMA
      BESTA=XFIT(1); BESTB=XFIT(2); BESTC=XFIT(3)
      BESTAVAIL=NAVAIL
      BESTBINS(J1)=MAXBIN
   ENDIF
10 IF (MAXBIN.EQ.HIGHESTDIRECT(J1)) EXIT
! 10 IF (MAXBIN.EQ.LOWESTDIRECT(J1)+MINFIT) EXIT
!    MAXBIN=MAX(MAXBIN-BININC,LOWESTDIRECT(J1)+MINFIT)
   MAXBIN=MIN(MAXBIN+BININC,HIGHESTDIRECT(J1))
!  IF (.NOT.CONVERGED) THEN
!     XFIT(1)=0.0D0; XFIT(2)=BPARAM; XFIT(3)=0.0D0
!  ENDIF
ENDDO

IF (ALLZEROQ(J1)) RETURN

!  IF ((BESTSIGMA.GT.SIGMAX).OR.(.NOT.CONVERGED)) THEN
IF (BESTSIGMA.GT.SIGMAX) THEN
   PRINT '(A,G20.10)','convergence failure or residual standard deviation too large, discarding fit. sigma=',BESTSIGMA
   ALLZEROQ(J1)=.TRUE.
ELSE
   PRINT '(A,G15.8,A,3G15.8,A,I8,A,I8,A)','best sigma=',BESTSIGMA,' A,B,C=',BESTA,BESTB,BESTC, &
  &    ' for ', &
  &    BESTAVAIL, '/', BESTBINS(J1)-LOWESTDIRECT(J1)+1,' bins'
   PRINT '(A,I8)','most visited pe bin for this q bin=',MOSTVISITED(J1)
ENDIF
PRINT '(A)',' '
FITA(J1)=BESTA
FITB(J1)=BESTB
FITC(J1)=BESTC

!   WRITE(FNAME,'(I8)') J1
!   FNAME='best.fit.' // TRIM(ADJUSTL(FNAME))
!   OPEN(UNIT=1,FILE=FNAME,STATUS='UNKNOWN')
!   DO J2=LOWESTDIRECT(J1),MAXBIN
!      NCOUNT=VARMAP2(J2,J1)
!      WRITE(1,'(I8,2G20.10)') J2,OLDVAR(NCOUNT),(KFAC+KFAC*(ENERGY(J2)-QE)*EXP(BESTA))*LOG(ENERGY(J2)-QE)+BESTB
!   ENDDO
!   CLOSE(1)

END SUBROUTINE FITABC

!
!  Calculate Z0, Z1 and Z2 over the required T range. Omit factors of (kT/h)^kappa,
!  which cancel.
!
SUBROUTINE DUMPCV(NBINS,ALLZERO,ENERGY,TMIN,TINT,EVAR,DELTA,KAPPA,NTEMP)
IMPLICIT NONE
INTEGER NBINS,J2,KAPPA,J1,NTEMP
DOUBLE PRECISION ENERGY(NBINS),TMIN,TINT,EVAR(NBINS),DELTA,EREF,Z0,Z1,Z2,ONEMEXP,TEMPERATURE,NCOUNT
LOGICAL ALLZERO(NBINS)
!
! Find energy of lowest visited bin to avoid exponential underflow.
!
DO J2=1,NBINS
   IF (.NOT.ALLZERO(J2)) THEN
      EREF=ENERGY(J2)
      EXIT
   ENDIF
ENDDO
DO J1=1,NTEMP
   Z0=0.0D0
   Z1=0.0D0
   Z2=0.0D0
   TEMPERATURE=TMIN+(J1-1)*TINT
   NCOUNT=0
   DO J2=1,NBINS
!     IF (ALLZERO(J2)) CYCLE 
      IF (EVAR(J2).EQ.0.0D0) CYCLE
!     PRINT '(A,I8,3G20.10)','J2,ENERGY(J2),ERER,TEMPERATURE,arg=', &
! & J2,ENERGY(J2),ENERGY(1),TEMPERATURE-(ENERGY(J2)-EREF)/TEMPERATURE
      Z0=Z0+EVAR(J2)*EXP(-(ENERGY(J2)-EREF)/TEMPERATURE)
      Z1=Z1+EVAR(J2)*EXP(-(ENERGY(J2)-EREF)/TEMPERATURE)*(ENERGY(J2)-EREF)
      Z2=Z2+EVAR(J2)*EXP(-(ENERGY(J2)-EREF)/TEMPERATURE)*(ENERGY(J2)-EREF)**2
   ENDDO
   IF (DELTA/TEMPERATURE.LT.1.0D-7) THEN
      ONEMEXP=-DELTA/TEMPERATURE
   ELSE
      ONEMEXP= 1.0D0-EXP(DELTA/TEMPERATURE)
   ENDIF
   WRITE(1,'(6G20.10)') TEMPERATURE, Z0, Z1, Z2, &
  &                     KAPPA*TEMPERATURE/2.0D0 + 1.0D0*(TEMPERATURE + DELTA/ONEMEXP) + Z1/Z0, &
  &                     KAPPA/2.0D0 + &
  &                      1.0D0*(1.0D0 - DELTA**2*EXP(DELTA/TEMPERATURE)/(ONEMEXP**2*TEMPERATURE**2)) &
  &                      - (Z1/(Z0*TEMPERATURE))**2 + Z2/(Z0*TEMPERATURE**2)
ENDDO

END SUBROUTINE DUMPCV
