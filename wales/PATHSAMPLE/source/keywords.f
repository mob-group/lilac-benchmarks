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
      SUBROUTINE KEYWORDS
      USE PORFUNCS
      USE NODES, ONLY: JPN, GETNODES, NNODES
      USE COMMONS
      USE RIGIDBODYMOD, ONLY: CAPSOMER, NUMRBTYPES, CAPSOMERDEFS
      IMPLICIT NONE

      INTEGER ITEM, NITEMS, LOC, LINE, NCR, NERROR, IR, ISTAT, NDUMMY2, LAST
      LOGICAL CAT
      COMMON /BUFINF/ ITEM, NITEMS, LOC(80), LINE, SKIPBL, CLEAR, NCR, NERROR, IR, ECHO, LAST, CAT
      INTEGER J1, NDUMMY, J2, J3, NELLIPSOIDS, MLPDATA
      DOUBLE PRECISION DUMMY, DBSIGBB, RBSITEPOS, MLPLAMBDA
      DOUBLE PRECISION, ALLOCATABLE :: RBCENTRE(:,:)
      DOUBLE PRECISION, ALLOCATABLE :: MLPDISTHI(:), MLPDISTHO(:)
      INTEGER, ALLOCATABLE :: MLPINDEXI(:), MLPINDEXO(:)
      INTEGER K1, II1

      LOGICAL END, SKIPBL, CLEAR, ECHO, PERMFILE, RBSYMTEST
      CHARACTER(LEN=100) WW
      CHARACTER(LEN=80) RWENERGYFILE
      CHARACTER(LEN=20) WORD
      CHARACTER(LEN=1) LOYNO
      CHARACTER(LEN=9) UNSTRING
      CHARACTER(LEN=4) LABEL
      CHARACTER(LEN=8) DUMMYSTRING 
       
      TYPE (CAPSOMER), ALLOCATABLE :: TEMPCAPSOMERDEFS(:)
      INTEGER :: MAXRBTYPES
      INTEGER :: MAXNSETS

      CHARACTER(LEN=8) TFOLD_ON
      CALL GETENV("TFOLD", TFOLD_ON)

      DBPT     = .FALSE.
      DBPTDT   = .FALSE.
      DMBLPYT  = .FALSE.
      MSSTOCKT = .FALSE.
      NTIPT    = .FALSE.
      PAHAT    = .FALSE.
      PAPT     = .FALSE.
      PATCHYDT = .FALSE.
      SILANET  = .FALSE.
      STOCKAAT = .FALSE.
      EFIELDT  = .FALSE.
      DEBUG=.FALSE.
      EXTRACTMINT=.FALSE.
      EXTRACTMINFILET=.FALSE.
      EXTRACTTST=.FALSE.
      EXTRACTTSFILET=.FALSE.
      ADDPT=.FALSE.
      ADDPT2=.FALSE.
      ADDPT3=.FALSE.
      DUMPGROUPST=.FALSE.
      DSCALE=3.0D0
      PSCALE=0.5D0
      CONNECTIONS=0
      MAXTSATTEMPTS=10
      ISEED=1
      NATTEMPT=0
      PERTVALUE=0.9D0
      PERTMAX=2.0D0*PERTVALUE
      PERTMIN=0.5D0*PERTVALUE
      EDIFFTOL=1.0D-8
      GEOMDIFFTOL=1.0D-1
      IDIFFTOL=1.0D-3
      DIRECTION='AB'
      ENSEMBLE='T'
      TWOD=.FALSE.
      BULKT=.FALSE.
      ANGLEAXIS=.FALSE.
      NUMRBTYPES = 0
      MAXRBTYPES = 0
      FROMLOWESTT=.FALSE.
      NOLABELST=.FALSE.
      DUMMYSTRING=''
      MACHINE=.FALSE.

      AMBERT=.FALSE.
      AMHT=.FALSE.
      AMHQT=.FALSE.
      AMHQENGMINT=.FALSE.
      AMHQCONTT=.FALSE.
      AMHRMSDT=.FALSE.
      AMHRELQT=.FALSE.
      AMH_RELCOT=.FALSE.
      AMHALLATOMMINT=.FALSE.
      AMHALLATOMTST=.FALSE.
      NOFRQS=.FALSE.
      NOINVERSION=.FALSE.

C davidg: introduced userpot here:
      USERPOTT=.FALSE.
      OPEPT=.FALSE.

      JPN=1
      NNODES=1
      NATOMS=-1
      CHARMMT=.FALSE.
      STARTFROMPATH=.FALSE.
      READMINT=.FALSE.
      ADDMINXYZT=.FALSE.
      ADDPATH=.FALSE.
      KMCCOMMITT=.FALSE.
      GTT=.FALSE.
      NGTT=.FALSE.
      NGTDISCONNECTALL=.FALSE.
      NGTSWITCH=0.3D0
      NGTSIZE=11000
      NGTCRSWITCH=2.0D0
C     Sem: begin GT2 controls
      GT2T=.FALSE.
      GT2SPARSE=.TRUE.
      GT2SWITCH=.TRUE.
      GT2DisconnectSources=.TRUE.
      GT2ALTPBB=.TRUE.
      GT2RESCALE=.FALSE.
      GT2NORMALISE=.FALSE.
      GT2RSWITCH=0.08D0
      GT2PTOL=1.0D-5
C     Sem: end GT2 controls
      GTINT=0
      REGROUPT=.FALSE.
      REGROUPTHRESH=-1.0D100
      REGROUPRATET=.FALSE.
      REGROUPRATETHRESH=1.0D100
      REGROUPPET=.FALSE.
      REGROUPPETHRESH=-1.0D100
      REGROUPFREET=.FALSE.
      REGROUPPERSISTT=.FALSE.
      ALLOWABT=.FALSE.
      RFMULTIT=.FALSE.
      TIMESCALE=1.0D0
      RFMULTITINC=1.0D0
      RFMULTITLOW=1.0D0
      RFMULTIN=10
      RFKMCT=.FALSE.
      REGROUPKMCT=.FALSE.
      ONEREGROUPT=.FALSE.
      RFKMCTRATE=1.0D0
      RFKMCTINC=1.0D0
      RFKMCTSTART=1.0D0
      RFKMCN=10
      RFKMCSTEPS=10
      PFSHIFT=0.0D0
      REGROUPFREEABT=.FALSE.
      REGROUPFREETHRESH=-1.0D100
      PABCONV=1.0D-8
      OMEGA=1.0D0 ! GAUSS-SEIDEL ITERATION
      KMCT=.FALSE.
      NCONNMIN=0
      MAXBREAK=1.0D6
      NKMCCYCLES=100.0D0
      PAIRTHRESH=1.0D0
      UNRST=.FALSE.
      NINTS=-1
      NCPU=1
      UNSTRING='UNDEFINED'
      ZSYM='UNKNOWN'
      EXEC='UNDEFINED'
      EXECGMIN='UNDEFINED'
      PATHNAME='UNDEFINED'
      NOPOINTS=.FALSE.
      DIJKSTRAT=.FALSE.
      DIJKSTRAWAITT=.FALSE.
      DIJPAIRT=.FALSE.
      BARRIERSORT=.FALSE.
      BARRIERSHORT=.FALSE.
      RATESHORT=.FALSE.
      DIJINITT=.FALSE.
      DIJINITFLYT=.FALSE.
      DIJINITSTARTT=.FALSE.
      DIJINITCONTT=.FALSE.
      DIJPRUNET=.FALSE.
      PAIRSIGNORET=.FALSE.
      PRUNECYCLET=.FALSE.
      NPRUNE=5
      PAIRDISTMAX=100
      NRANDOMMETRIC=100
      RANDOMMETRICT=.FALSE.
      PAIRDIST1=0
      PAIRDIST2=0
      TSTHRESH=HUGE(1.0D0)
      MAXBARRIER=HUGE(1.0D0)
      MAXDOWNBARRIER=HUGE(1.0D0)
      MINBARRIER=-HUGE(1.0D0)
      MINGAPT=.FALSE.
      MINGAPRATIOT=.FALSE.
      MINGAPINP=0.0D0
      COSTFUNCTIONPOWER=1
      EXPCOSTFUNCTION=.FALSE.
      INDEXCOSTFUNCTION=.FALSE.
      COPYFILES=''
      COPYOPTIMT=.FALSE.
      NPFOLD=0
      PFOLDCONV=1.0D-4
      if(TFOLD_ON == "") then
        TFOLDT=.FALSE.
      else
        TFOLDT=.TRUE.
      endif
      NTFOLD=1.0D5 ! real not integer !
      TOMEGA=1.0D0
      TFOLDTHRESH=1.0D-5
      MAXATTEMPT=1
      CALCORDERT=.FALSE.
      OSTART=1
      OFINISH=1
      CONNECTREGIONT=.FALSE.
      CONNECTMIN1=1
      CONNECTMIN2=1
      CONNECTMIN2F=-1
      CONNECTDIST=1.0D10
      SHORTCUTT=.FALSE.
!     PTAUT=.FALSE.
      NPAIRFRQ=0
      MERGEDBT=.FALSE.
      PERSISTT=.FALSE.
      PEQTHRESH=1.0D-100
      PERTHRESH=0.0D0
      PERSISTAPPROXT=.TRUE.
      ALLCOMPONENTST=.FALSE.
      ALLCOMPS=''
      UNTRAPT=.FALSE.
      UNTRAPMETRICT=.FALSE.
      METRICUPAIR=0
      METMATMAX=2000 
      EUNTRAPTHRESH=1.0D100
      FREEPAIRT=.FALSE.
      PERMDIST=.FALSE.
      MAXNSETS=3
      ORBITTOL=1.0D-3
      LPERMDIST=.FALSE.
      LPDGEOMDIFFTOL=0.3D0
      LOCALPERMNEIGH=4
      LOCALPERMCUT=0.2D0
      LOCALPERMCUT2=10.0D0
      PERMISOMER=.FALSE.
      NTAG=0
      TAGT=.FALSE.
      FREEZE=.FALSE.
      NFREEZE=0
      PLANCK=1.0D0
      DUMMYRUNT=.FALSE.
      DUMMYTST=.FALSE.
      REWEIGHTT=.FALSE.
      KSHORTESTPATHST = .FALSE.
      KSHORT_FULL_PRINTT = .FALSE.
      NPATHS = 0
      BHINTERPT=.FALSE.
      BHACCREJ=0.5D0
      BHSTEPSIZE=0.4D0
      BHCONV=0.01D0
      BHSTEPS=1000
      BHTEMP=1.0D0
      BHK=1.0D0
      ICINTERPT=.FALSE.

      USEPAIRST=.FALSE.
      CONNECTPAIRST=.FALSE.
      LOWESTFRQT=.FALSE.
      IMFRQT=.FALSE.
      EVCUT=2.0D-6

      BISECTT=.FALSE.
      DIAGT=.FALSE.
      ARNOLDIT=.FALSE.
      SLURMT=.FALSE.
      CUDAT=.FALSE.
      CVMINIMAT=.FALSE.
      CVT=.FALSE.
      SHANNONT=.FALSE.
      SHANNONRT=.FALSE.
      SHANNONZT=.FALSE.
      NPEQ=100
      MICROTHERMT=.FALSE.
      DOST=.FALSE.
      CHECKCONNECTIONST=.FALSE.
      NEWCONNECTIONST=.FALSE.
      CONNMINSTART=1
      CLOSEFILEST=.FALSE.
      PULLT=.FALSE.
      PRINTSUMMARYT=.FALSE.
!     DC430 >
      RBAAT  = .FALSE.
      RBSYMT = .FALSE.
      FRICTIONT=.FALSE.
      GAMMAFRICTION=0.0D0
      REMOVEUNCONNECTEDT=.FALSE.
      UNCONNECTEDS='AB'
      OHCELLT=.FALSE.
      ATOMMATCHDIST=.FALSE.
      ATOMMATCHFULL=.FALSE.
      SSHT=.FALSE.
      RATESCYCLET=.FALSE.
      NIMET=.FALSE.
      NIHEAM7T=.FALSE.
      NIH2LEPST=.FALSE.
!
! Constraint potential
!
      INTERPCOSTFUNCTION=.FALSE.
      INTCONSTRAINTT=.FALSE.
      INTCONSTRAINTTOL=0.1D0
      INTCONSTRAINTDEL=1.0D5
      INTCONSTRAINTREP=1.0D0
      INTCONSTRAINREPCUT=20.0D0
      INTREPSEP=0
      INTCONSEP=10000
      CHECKCONINT=.FALSE.
      MAXCONUSE=3
      INTFREEZET=.FALSE.
      INTFREEZETOL=0.1D0
!
! LJ interpolation potential parameters
!
      INTLJT=.FALSE.
      INTLJDEL=0.1D0
      INTLJEPS=1.0D0
      ALLTST=.FALSE.
      SLEEPTIME1=1.0D0
      SLEEPTIME2=1.0D0
C
C
C sf344> 
C
      DOCKT=.FALSE.
      DSTAGE(:)=.TRUE.
      MACROIONT=.FALSE.
!
! SIS epidemiological model
! 
      SIST=.FALSE.
      SMAX=0
      IMAX=0
      POPSS=0
      SISMU=0.0D0
      SISKAPPA=0.0D0
      SISBETA=0.0D0
!
! Optional UNTRAP argument
!
      EDELTAMIN=TINY(1.0D0)
      ELOWBAR=TINY(1.0D0)
      EHIGHBAR=HUGE(1.0D0)
      BAILDIST=TINY(1.0D0)
!
! To document...
!
!
! Neural network potential
!
      MLP3T=.FALSE.
      MLPB3T=.FALSE.
      MLPNEIGH=1
      MLPLAMBDA=0.0D0
!
! MK trap potential
!
      MKTRAPT=.FALSE.

      PHI4MODT=.FALSE.
      MLLJAT3=.FALSE.
      RELATIVEET=.FALSE.
      RATETARGETT=.FALSE.
      RATETARGETAB=HUGE(1.0)
      RATETARGETBA=HUGE(1.0)
      TARGETHIT=.FALSE.
      NRANROT=0
      CHECKSPT=.FALSE.
      CHECKMINT=.FALSE.
      CHECKTST=.FALSE.
      CHECKSPS=-1
      CHECKSPF=-1
      INITIALDIST=.FALSE.

      NZEROS=0


      DISTANCET=.FALSE.
      DISTANCETO=1
      DISTANCETO1=1
      DISTANCETO2=1

! hk286
      GTHOMSONT = .FALSE.
      RIGIDINIT = .FALSE.
      DEGFREEDOMS = 1 !   3*NATOMS this was a bug, NATOMS=-1 here !!!!!

! kr366: get frqs
      GETMINFRQST = .FALSE.
      GETTSFRQST = .FALSE.

! sn402: keywords for manipulating pairfiles. Must be used with USEPAIRS
      MAKEPAIRS = .FALSE.
      SKIPPAIRST = .FALSE.

! kr366: connect unconnected minima
      CONUNCONT=.FALSE.
      CONNECTLOWESTT=.FALSE.
      CONNECTETHRESHT=.FALSE.
      CONNECTDTHRESHT=.FALSE.
      REFMIN=0
      NATT=2
!
! Forbid overall translation and rotation in distance/alignment routines.
!
      NOTRANSROTT=.FALSE.

C
C SAT: Now, there it is! Unit 5 is a standard input unit, and we must not use it!
C 5--->959
C
      OPEN (959,FILE='pathdata',STATUS='OLD')
C
C  Read from unit 959
C
      IR=959
190   CALL INPUT(END)
      IF (.NOT. END) THEN
        CALL READU(WORD)
      ENDIF

      IF (END .OR. WORD .EQ. 'STOP') THEN
        CLOSE(IR)
        RETURN
      ENDIF

      IF (WORD.EQ.'    '.OR.WORD.EQ.'NOTE'.OR.WORD.EQ.'COMMENT'
     &                  .OR.WORD.EQ.'\\'.OR.WORD.EQ."!"
     &                  .OR.WORD.EQ."#" ) THEN 
         GOTO 190
C
C  Add the minima in the specified min.data.info file to the
C  database. Same as READMIN, since I keep thinking it should
C  be called ADDMIN!
C
      ELSE IF (WORD.EQ.'ADDMIN') THEN
         READMINT=.TRUE.
         CALL READA(MINNAME)

      ! An additional keyword to use with ADDMIN or READMIN to make a pairfile for later use with USEPAIRS.
      ! The sequence of pairs saved in the pairfile will match the sequence in the min.data.info file which is
      ! being added to the database (duplicate minima will be included in the pairfile every time they appear).
      ! The name of the pairfile to be written out is a required argument.

      ELSE IF (WORD.EQ.'MAKEPAIRS') THEN
          MAKEPAIRS = .TRUE.
          CALL READA(MAKEPAIRSFILE)
C
C  Add the minima in the specified xyz file to the
C  database, after running OPTIM for each set of coordinates
C  using a supplied odata.addminxyz file.
C
      ELSE IF (WORD.EQ.'ADDMINXYZ') THEN
         ADDMINXYZT=.TRUE.
         CALL READA(ADDMINXYZNAME)
C
C  Add the minima and transition= states in path.info.<PATHNAME>
C  to an existing database. The end points are NOT assumed to
C  belong to the A and B sets.
C
      ELSE IF (WORD.EQ.'ADDPATH') THEN
         ADDPATH=.TRUE.
         CALL READA(PATHNAME)
C
C  ADDPT determines whether we call ADDPERM to add permutational isomers of every stationary point to the
C  min.data and ts.data databases. Speeds up 2DLJ7! For other tagged atom situations the number of isomers
C  will equal the number of atoms for C1 symmetry, so this is probably a bad idea.
C
      ELSE IF (WORD.EQ.'ADDPERM') THEN
         ADDPT=.TRUE.
C
C  ADDPT2 determines whether we call ADDPERM2 to add permutational isomers of every stationary point to the
C  min.data and ts.data databases. Does not assume any tagged atoms. Assumes that all atoms in a 
C  specified range are permutable - does not use perm.allow, so could permute atoms that aren't considered
C  permutable in perm.allow.
C
      ELSE IF (WORD.EQ.'ADDPERM2') THEN
         ADDPT2=.TRUE.
         PTSTART=1
         PTFINISH=NATOMS
         IF (NITEMS.GT.1) CALL READI(PTSTART)
         IF (NITEMS.GT.2) CALL READI(PTFINISH)
         PRINT '(A,2I6)','keywords> Will add new transition states obtained by permutations of atoms in the range ',PTSTART,PTFINISH
C
C  ADDPT3 determines whether we call ADDPERM3 to add permutational isomers of every stationary point to the
C  min.data and ts.data databases according to perm.allow.
C
      ELSE IF (WORD.EQ.'ADDPERM3') THEN
         ADDPT3=.TRUE.
         PERMDIST=.FALSE.
         PRINT '(A,2I6)','keywords> Will add new transition states obtained by permutations specified in perm.allow'
         PRINT '(A,2I6)','keywords> PERMDIST is reset to false'
C
C  Allow A and B sets to merge in regroupfree2.
C
      ELSE IF (WORD.EQ.'ALLOWAB') THEN
         ALLOWABT=.TRUE.
C
C  Use the angle-axis system for rigid bodies
C
      ELSE IF (WORD.EQ.'ANGLEAXIS') THEN
         ANGLEAXIS=.TRUE.
      ELSE IF (WORD.EQ.'ANGLEAXIS2') THEN
         ANGLEAXIS2=.TRUE.
C
C  Allow the same ts to appear with different minima in the database.
C
      ELSE IF (WORD.EQ.'ALLTS') THEN
         ALLTST=.TRUE.
C
C  sf344> Set AMBER9 potential.
C
       ELSE IF (WORD.EQ.'AMBER9') THEN
          AMBERT=.TRUE.
       ELSE IF (WORD.EQ.'AMBER12') THEN
          AMBERT=.TRUE.
       ELSE IF (WORD.EQ.'NAB') THEN
          AMBERT=.TRUE.
C
C  csw34> set same flag for AMBER12 for now - should be safe!
C
       ELSE IF (WORD.EQ.'AMBER12') THEN
          AMBERT=.TRUE.
C
C  Set AMH potential.
C
       ELSE IF (WORD.EQ.'AMH') THEN
          AMHT=.TRUE.

C
C  Calculated  Q between minima and native .
C
      ELSE IF (WORD.EQ.'AMHQ') THEN
         AMHQT=.TRUE.
         CALL READI(WHICHMIN)
          PRINT '(A,I6,A,I6)','keywords> Calculate AMH Q  ',WHICHMIN
         IF (NITEMS.GT.2) THEN
          PRINT '(A)','keywords> ERROR - AMHQ '
          STOP
         ENDIF

C
C  Calculated  Q between minima and energy minimum.
C
       ELSE IF (WORD.EQ.'AMHQENGMIN') THEN
        AMHQENGMINT=.TRUE.
        CALL READI(WHICHMIN)
         PRINT '(A,I6,A,I6)','keywords> Calculate AMH Q 2 ENERGY MINIMUM ',WHICHMIN
        IF (NITEMS.GT.2) THEN
         PRINT '(A)','keywords> ERROR - AMH Q 2 ENERGY MINIMUM '
         STOP
        ENDIF

C
C  Calculated  Q_contact  between minima and native .
C
      ELSE IF (WORD.EQ.'AMHQCONT') THEN
         AMHQCONTT=.TRUE.
         CALL READI(WHICHMIN)
         CALL READF(QCONTCUT)
          PRINT '(A,I6)','keywords> Calculate AMH Q CONTACT ',WHICHMIN
          PRINT '(A,G9.2)','keywords> Within a distance cutoff ',QCONTCUT
         IF (NITEMS.GT.3) THEN
          PRINT '(A)','keywords> ERROR - AMH Q CONTACT '
          STOP
         ENDIF
C
C  Calculated  RMSD between minima and native .
C
      ELSE IF (WORD.EQ.'AMHRMSD') THEN
         AMHRMSDT=.TRUE.
         CALL READI(WHICHMIN)
          PRINT '(A,I6,A,I6)','keywords> Calculate AMH RMSD  ',WHICHMIN
         IF (NITEMS.GT.2) THEN
          PRINT '(A)','keywords> ERROR - AMHRMSD '
          STOP
         ENDIF
C
C  Calculated Relative Qs between minima.
C
      ELSE IF (WORD.EQ.'AMHRELQ') THEN
         AMHRELQT=.TRUE.
         CALL READI(QRELONE)
         CALL READI(QRELTWO)
          PRINT '(A,I6,A,I6)','keywords> AMHRELQ min 1 ',QRELONE,' min 2',QRELTWO
         IF (NITEMS.GT.3) THEN
          PRINT '(A)','keywords> ERROR - AMHRELQ GT 2'
          STOP
         ENDIF

C
C  Calculated Relative Contact Order
C
       ELSE IF (WORD.EQ.'AMH_RELCO') THEN
         AMH_RELCOT=.TRUE.
         CALL READI(WHICHMIN)
         CALL READF(RELCOCUT)
          PRINT '(A,I6)','keywords> Calculate AMH Relative Contact Order ',WHICHMIN
          PRINT '(A,G20.10)','keywords> Within a distance cutoff ',RELCOCUT

         IF (NITEMS.GT.3) THEN
          PRINT '(A)','keywords> ERROR - RELCO GT 2'
          STOP
         ENDIF

C
C  Output AMH all-atom MIN .
C
      ELSE IF (WORD.EQ.'AMHALLATOMMIN') THEN
         AMHALLATOMMINT=.TRUE.
         PRINT '(A,I6,A,I6)','keywords> output AMH all atom MIN structures'
C
C  Output AMH all-atom TS .
C
      ELSE IF (WORD.EQ.'AMHALLATOMTS') THEN
         AMHALLATOMTST=.TRUE.
         PRINT '(A,I6,A,I6)','keywords> output AMH all atom TS structures'

C
C  Rates from eigenvalues obtained with weighted Arnoldi subspace method.
C
      ELSE IF (WORD.EQ.'ARNOLDI') THEN
         ARNOLDIT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
C
C  Specify atom matching for the distance calculation.
C
      ELSE IF (WORD.EQ.'ATOMMATCHDIST') THEN
         ATOMMATCHDIST=.TRUE.
         WRITE(*,'(A)') ' Atom matching for distance calculation'
      ELSE IF (WORD.EQ.'ATOMMATCHFULL') THEN
         ATOMMATCHDIST=.TRUE.
         ATOMMATCHFULL=.TRUE.
         WRITE(*,'(A)') ' Atom matching for distance calculation'
         WRITE(*,'(A)') ' WARNING - inefficient atom matching for a deterministic result'
C
C  Parameters for OPTIM basin-hopping interpolation jobs.
C
      ELSE IF (WORD.EQ.'BHINTERP') THEN
         BHINTERPT=.TRUE.
         CALL READF(BHDISTTHRESH)
         CALL READF(BHMAXENERGY)
         CALL READI(BHSTEPS)
         CALL READF(BHCONV)
         CALL READF(BHTEMP)
         CALL READF(BHSTEPSIZE)
         CALL READF(BHACCREJ)
         CALL READF(BHK)
         CALL READF(BHSFRAC)
         IF (NITEMS.GT.10) CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'ICINTERP') ICINTERPT=.TRUE.
C
C  Parameters for OPTIM bisection interpolation jobs.
C
      ELSE IF (WORD.EQ.'BISECT') THEN
         BISECTT=.TRUE.
         CALL READF(BISECTMINDIST)
         CALL READF(BISECTMAXENERGY)
         CALL READI(BISECTSTEPS)
         CALL READI(BISECTMAXATTEMPTS)
         IF (NITEMS.GT.5) CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'ICINTERP') ICINTERPT=.TRUE.
C
C  Bulk with periodic boundary conditions.
C
      ELSE IF (WORD.EQ.'BULK') THEN
         BULKT=.TRUE.
         CALL READF(BOXLX)
         CALL READF(BOXLY)
         CALL READF(BOXLZ)
C
C  Calculate order parameter for all stationary points and exit.
C  The order parameter routine will be system specific and must
C  replace the dummy routine provided.
C  Minima and transition states will be listed in separate files
C  for points with order parameters above/below the specified threshold,
C  with default value 0.0.
C
      ELSE IF (WORD.EQ.'CALCORDER') THEN
         CALCORDERT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(OSTART)
         IF (NITEMS.GT.2) CALL READI(OFINISH)
C
C Capsomere definition
C
      ELSE IF (WORD.EQ.'CAPSID') THEN
         ANGLEAXIS=.TRUE.
         RIGIDBODY=.TRUE.
         IF (.NOT.ALLOCATED(CAPSOMERDEFS)) THEN
            MAXRBTYPES = 10
            ALLOCATE(CAPSOMERDEFS(MAXRBTYPES))

         ELSE IF (NUMRBTYPES.EQ.MAXRBTYPES) THEN
            ALLOCATE(TEMPCAPSOMERDEFS(MAXRBTYPES))
            TEMPCAPSOMERDEFS = CAPSOMERDEFS

            MAXRBTYPES = MAXRBTYPES + 10
            DEALLOCATE(CAPSOMERDEFS)
            ALLOCATE(CAPSOMERDEFS(MAXRBTYPES))

            CAPSOMERDEFS(1:NUMRBTYPES) = TEMPCAPSOMERDEFS
            DEALLOCATE(TEMPCAPSOMERDEFS)
         ENDIF

         NUMRBTYPES = NUMRBTYPES + 1

         CAPSOMERDEFS(NUMRBTYPES)%NBASALSITES = 5

         CALL READF(CAPSOMERDEFS(NUMRBTYPES)%RHO)
         CALL READF(CAPSOMERDEFS(NUMRBTYPES)%EPSILONREP)
         CALL READF(CAPSOMERDEFS(NUMRBTYPES)%RADIUS)
         IF (NITEMS.GT.4) THEN
            CALL READF(CAPSOMERDEFS(NUMRBTYPES)%HEIGHT)
            IF (NITEMS.GT.5) CALL READI(CAPSOMERDEFS(NUMRBTYPES)%NBASALSITES)
         ELSE
            CAPSOMERDEFS(NUMRBTYPES)%HEIGHT = 0.5D0*CAPSOMERDEFS(NUMRBTYPES)%RADIUS
         ENDIF
C
C  Set CHARMM potential. NDIHE is the number of dihedral angles in the system, which
C is used in perturb and tssearch to perturb the system randomly
C
      ELSE IF (WORD.EQ.'CHARMM') THEN
         CHARMMT=.TRUE.
         CALL READI(NDIHE)
C
C  Check minimum number of connections for each minimum at startup.
C
      ELSE IF (WORD.EQ.'CHECKCONNECTIONS') THEN
         CHECKCONNECTIONST=.TRUE.
C
C  Check minimum number of connections for each minimum at startup.
C
      ELSE IF (WORD.EQ.'CHECKMIN') THEN
         CHECKSPT=.TRUE.
         CHECKMINT=.TRUE.
         COPYOPTIMT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(CHECKSPS)
         IF (NITEMS.GT.2) CALL READI(CHECKSPF)
C
C  Check minimum number of connections for each minimum at startup.
C
      ELSE IF (WORD.EQ.'CHECKTS') THEN
         CHECKSPT=.TRUE.
         CHECKTST=.TRUE.
         COPYOPTIMT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(CHECKSPS)
         IF (NITEMS.GT.2) CALL READI(CHECKSPF)
C
C  Check for internal minimum in constraint terms for INTCONSTRAINT
C
      ELSE IF (WORD.EQ.'CONINT') THEN
         CHECKCONINT=.TRUE.
C
C  Close and open all files as needed. The objective is to work around
C  the misinteraction of nfs with the Linux kernel, which causes random
C  control characters to be written to open nfs mounted files.
C
      ELSE IF (WORD.EQ.'CLOSEFILES') THEN
         CLOSEFILEST=.TRUE.
C
C  Minimum number of connections for each minimum.
C
      ELSE IF (WORD.EQ.'CONNECTIONS') THEN
         CALL READI(CONNECTIONS)
         IF (NITEMS.GT.2) CALL READI(MAXTSATTEMPTS)
C
C  Connection attempts from airs of minima in user specified file
C
      ELSE IF (WORD.EQ.'CONNECTPAIRS') THEN
         CONNECTPAIRST=.TRUE.
         CALL READA(CONNECTPAIRSFILE)
C
C  Specify that the database is grown by attempting all connections between
C  known minima within a distance cutoff.
C
      ELSE IF (WORD.EQ.'CONNECTREGION') THEN
         CONNECTREGIONT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(CONNECTMIN1)
         IF (NITEMS.GT.2) CALL READI(CONNECTMIN2)
         IF (NITEMS.GT.3) CALL READF(CONNECTDIST)
         IF (NITEMS.GT.4) CALL READI(CONNECTMIN2F)
         CONNECTMIN2SAVE=CONNECTMIN2
C
C Connect unceonnected minima to AB set using either the lowest energy minima
C or the minima closest in energy or distance to a reference
C
      
      ELSE IF (WORD.EQ.'CONNECTUNC') THEN
         CONUNCONT=.TRUE.
         PRINT '(A)',' keywords> Attempt to connect unconnected minima'
         CALL READU(WW)
         IF (TRIM(ADJUSTL(WW))=='LOWEST') THEN
            CONNECTLOWESTT=.TRUE.
            IF (NITEMS.GT.2) CALL READI(NATT)
            PRINT '(A)','keywords> Use lowest energy minima'
            PRINT '(A,I8,A)','keywords>',NATT,'different connections per minimum'
         ELSE IF (TRIM(ADJUSTL(WW))=='EREF') THEN
            CONNECTETHRESHT=.TRUE.
            CALL READI(REFMIN)
            PRINT '(A,I8)','keywords> Use minima closest in energy to:',REFMIN
         ELSE IF (TRIM(ADJUSTL(WW))=='DISTREF') THEN
            CONNECTDTHRESHT=.TRUE.
            CALL READI(REFMIN)
            PRINT '(A,I8)','keywords>Use minima closest in distance to:',REFMIN
         ENDIF

C
C  Specify files to be copied to distributed nodes in addition to the default
C  odata.<pid> and finish.<pid> files.
C
      ELSE IF (WORD.EQ.'COPYFILES') THEN
C        PRINT '(A,I8)','NITEMS=',NITEMS
         DO J1=2,NITEMS
            CALL READA(WW)
            WRITE(COPYFILES,'(3A)') TRIM(ADJUSTL(COPYFILES)),' ',TRIM(ADJUSTL(WW))
C           PRINT '(3A)','WW,COPYFILES: ',WW,COPYFILES
         ENDDO
C
C  All OPTIM output files will be copied from distributed nodes and not deleted. 
C  This also occurs if DEBUG is .TRUE.
C
      ELSE IF (WORD.EQ.'COPYOPTIM') THEN
         COPYOPTIMT=.TRUE.
C
C  Use NCPU's by starting up to NCPU's OPTIM jobs.
C
      ELSE IF (WORD.EQ.'CPUS') THEN
         PBST=.FALSE.
         IF (NITEMS.GT.1) CALL READI(NCPU)
      ELSE IF (WORD.EQ.'CV') THEN
         CVT=.TRUE.
         CVTMIN=1.0D0
         CVTMAX=2.0D0
         CVTINC=0.1D0
         CALL READF(CVTMIN)
         CALL READF(CVTMAX)
         CALL READF(CVTINC)
      ELSE IF (WORD.EQ.'CVMINIMA') THEN
         CVMINIMAT=.TRUE.
         CALL READI(CVSTARTMIN)
         CALL READI(CVENDMIN)
         CALL READI(CVINCMIN)
C
C  Maximum number of connect cycles in CYCLE.
C
      ELSE IF (WORD.EQ.'CYCLES') THEN
         CALL READI(NATTEMPT)

      ELSE IF (WORD.EQ.'DB') THEN
! Set the EFIELDT separately, if required
         DBPT = .TRUE.
         CALL READF(DBSIGBB)
         NRBSITES = 3
         ALLOCATE(RBSITE(NRBSITES,3))
         CALL DEFDUM(DBSIGBB)

! hk286
!         ALLOCATE(SITEMASS(NRBSITES))
!         SITEMASS(:) = (/ 0.D0, 0.D0, 0.D0 /)

      ELSE IF (WORD.EQ.'DMBLPY') THEN
! Set the EFIELDT separately, if required
         DMBLPYT = .TRUE.
         CALL READF(DBSIGBB)
         NRBSITES = 3
         ALLOCATE(RBSITE(NRBSITES,3))
         ALLOCATE(SITEMASS(NRBSITES))
         CALL DEFDMBLPY(DBSIGBB)
         SITEMASS(:) = 1.D0

      ELSE IF (WORD.EQ.'DBTD') THEN
! Set the EFIELDT separately, if required
         DBPTDT = .TRUE.
         CALL READF(DBSIGBB)
         NRBSITES = 3
         ALLOCATE(RBSITE(NRBSITES,3))
         CALL DEFDUM(DBSIGBB)
! hk286
!         ALLOCATE(SITEMASS(NRBSITES))
!         SITEMASS(:) = (/ 0.D0, 0.D0, 0.D0 /)

C
C  Turn on debug printing.
C
      ELSE IF (WORD.EQ.'DEBUG') THEN
         DEBUG=.TRUE.

! hk286
! csw34> When using RIGIDINIT, this should be set using an argument to that keyword
      ELSE IF (WORD.EQ.'DEGFREEDOMS') THEN
         CALL READI(DEGFREEDOMS)

C
C  Rates from direct diagonalisation.
C
      ELSE IF (WORD.EQ.'DIAG') THEN
         DIAGT=.TRUE.
         NDIAGEIG=10
         DIAGSCALE=1.0D0
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
         IF (NITEMS.GT.2) CALL READI(NDIAGEIG)
         IF (NITEMS.GT.3) CALL READF(DIAGSCALE)
C
C  DIJINIT specifies a Dijkstra analysis to try and construct an initial path.
C  We only need files start and finish containing the end point coordinates for
C  DIJINITSTART - this keyword creates new min.A, min.B, points.min, points.ts,
C  min.data and ts.data files, so it is potentially dangerous!
C  DIJINITCONT specifies an initial path run where the above files already exist.
C  In a DIJINITSTART run 
C  file min.A is created with entries 1 1 and file min.B
C  with entries 1 2. min.data is set up with the two entries for
C  these minima. The coordinates are obtained via odata.start and odata.finish and put in records
C  1 and 2 in the points.min file. 
C  PAIRDISTMAX is the number of nearest neighbours to save in PAIRLIST judged by
C  the chosen distance or interpolation metric.
C
      ELSE IF ((WORD.EQ.'DIJINITCONT').OR.(WORD.EQ.'DIJPRUNE')) THEN
         IF (WORD.EQ.'DIJPRUNE') DIJPRUNET=.TRUE.
         IF (NITEMS.GT.1) THEN
               CALL READU(WW)
               IF (TRIM(ADJUSTL(WW))=='INDEX') then
                    INDEXCOSTFUNCTION = .TRUE.
               ELSEIF (TRIM(ADJUSTL(WW))=='EXP') then
                    EXPCOSTFUNCTION = .TRUE.
               ELSE IF (WW(1:1) /= ' ') then
                    READ(WW,'(I20)') COSTFUNCTIONPOWER
               ENDIF
         ENDIF
         IF (NITEMS.GT.2) CALL READI(PAIRDISTMAX)
         IF (NITEMS.GT.3) CALL READI(PAIRDIST1)
         IF (NITEMS.GT.4) CALL READI(PAIRDIST2)
         DIJINITCONTT=.TRUE.
         DIJINITT=.TRUE.
         DIJPAIRT=.TRUE.
         DIJKSTRAT=.TRUE.
      ELSE IF (WORD.EQ.'DIJINITSTART') THEN
         IF (NITEMS.GT.1) THEN
               CALL READU(WW)
               IF (TRIM(ADJUSTL(WW))=='INDEX') then
                    INDEXCOSTFUNCTION = .TRUE.
               ELSEIF (TRIM(ADJUSTL(WW))=='EXP') then
                    EXPCOSTFUNCTION = .TRUE.
               ELSE IF (WW(1:1) /= ' ') then
                    READ(WW,'(i20)') CostFunctionPower
               ENDIF
         ENDIF
         IF (NITEMS.GT.2) CALL READI(PAIRDISTMAX)
         DIJINITSTARTT=.TRUE.
         DIJINITT=.TRUE.
         DIJPAIRT=.TRUE.
         DIJKSTRAT=.TRUE.
      ELSE IF (WORD.EQ.'DIJINITCONTFLY') THEN
         IF (NITEMS.GT.1) THEN
               CALL READU(WW)
               IF (TRIM(ADJUSTL(WW))=='INDEX') then
                    INDEXCOSTFUNCTION = .TRUE.
               ELSEIF (TRIM(ADJUSTL(WW))=='EXP') then
                    EXPCOSTFUNCTION = .TRUE.
               ELSE IF (WW(1:1) /= ' ') then
                    READ(WW,'(I20)') COSTFUNCTIONPOWER
               ENDIF
         ENDIF
         DIJINITCONTT=.TRUE.
         DIJINITFLYT=.TRUE.
         DIJPAIRT=.TRUE.
         DIJKSTRAT=.TRUE.
      ELSE IF (WORD.EQ.'DIJINITSTARTFLY') THEN
         IF (NITEMS.GT.1) THEN
               CALL READU(WW)
               IF (TRIM(ADJUSTL(WW))=='INDEX') then
                    INDEXCOSTFUNCTION = .TRUE.
               ELSEIF (TRIM(ADJUSTL(WW))=='EXP') then
                    EXPCOSTFUNCTION = .TRUE.
               ELSE IF (WW(1:1) /= ' ') then
                    READ(WW,'(i20)') CostFunctionPower
               ENDIF
         ENDIF
         DIJINITSTARTT=.TRUE.
         DIJINITFLYT=.TRUE.
         DIJPAIRT=.TRUE.
         DIJKSTRAT=.TRUE.
C
C  Dijkstra analysis to find the largest contribution to the SS rate constant
C  from each B (or A) minimum.
C
      ELSE IF (WORD.EQ.'DIJKSTRA') THEN
         DIJKSTRAT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
C
C  Dijkstra analysis to find the path with the lowest sum of waiting times
C  times initial conditional probability.
C
      ELSE IF (WORD.EQ.'DIJKSTRAWAIT') THEN
         DIJKSTRAWAITT=.TRUE.
         DIJKSTRAT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
C
C  DIJPAIR specifies a Dijkstra analysis to propose the pairs of minima to
C  try and connect during the database construction.
C
      ELSE IF (WORD.EQ.'DIJPAIR') THEN
         CALL READI(MAXATTEMPT)
         IF (NITEMS.GT.2) CALL READI(NCONNMIN)
         IF (NITEMS.GT.3) CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'BARRIER') BARRIERSORT=.TRUE.
         DIJPAIRT=.TRUE.
         DIJKSTRAT=.TRUE.

C
C  Direction of sampling.
C
      ELSE IF (WORD.EQ.'DIRECTION') THEN
         CALL READA(DIRECTION)
C
C  Just run an optimised distance calculation for one minimum and some or all of the others.
C
      ELSE IF (WORD.EQ.'DISTANCE') THEN
         DISTANCET=.TRUE.
         CALL READI(DISTANCETO)
         CALL READI(DISTANCETO1)
         CALL READI(DISTANCETO2)
C
C  Docking routine to calculate the binding free energy of a ligand to a protein
C
      ELSE IF (WORD.EQ.'DOCK') THEN
         DOCKT=.TRUE.
         CALL READI(PARALLEL)
         IF (NITEMS.GT.2) THEN
          DO J1=1,6
            CALL READI(NDUMMY)
            IF(NDUMMY==0) DSTAGE(J1)=.FALSE.
          END DO
         END IF
      ELSE IF (WORD.EQ.'DOS') THEN
         DOST=.TRUE.
         DOSEMIN=1.0D0
         DOSEMAX=2.0D0
         DOSEINC=0.1D0
         CALL READF(DOSEMIN)
         CALL READF(DOSEMAX)
         CALL READF(DOSEINC)
C
C  DSCALE is the distance scaling for decisions about whether to connect a given
C  pair of local minima. PSCALE is the scaling for the difference in GPFOLD.
C  IF (EXP(-(DISTANCE-DSCALE)/DSCALE).LT.RANDOM) is used in cycle.
C
      ELSE IF (WORD.EQ.'DSCALE') THEN
         CALL READF(DSCALE)
C
C  Setting DUMMYRUN to true means that no OPTIM jobs are actually submitted.
C
      ELSE IF (WORD.EQ.'DUMMYRUN') THEN
         DUMMYRUNT=.TRUE.
C
C  Setting DUMMYTST to true means that we create phoney entries in ts.data
C  corresponding to nearest-neighbour minima. For use with BHINTERP and BISECTT.
C
      ELSE IF (WORD.EQ.'DUMMYTS') THEN
         DUMMYTST=.TRUE.
C
C  DUMPGROUPS specifies that the groups of potential energy minima and transition
C  states obtained from subroutine regroupfree2 should be saved in files
C  minima_groups and ts_groups
C
      ELSE IF (WORD.EQ.'DUMPGROUPS') THEN
         DUMPGROUPST=.TRUE.
C
C  Energy difference  criterion for distinguishing stationary points
C  Can also be specified with ETOL
C
      ELSE IF (WORD.EQ.'EDIFFTOL') THEN
         CALL READF(EDIFFTOL)

!  Turns on the electric field   
      ELSE IF (WORD.EQ.'EFIELD') THEN
         EFIELDT = .TRUE.
C
      ELSE IF (WORD.EQ.'ENERGY') THEN
         CALL READF(TOTALE)
         ENSEMBLE='E'
C
C  Energy difference  criterion for distinguishing stationary points
C
      ELSE IF (WORD.EQ.'ETOL') THEN
         CALL READF(EDIFFTOL)
C
C  Threshold for zero eigenvalues used to skip triples in GETALLPATHS
C
      ELSE IF (WORD.EQ.'EVCUT') THEN
         CALL READF(EVCUT)
C
C  OPTIM executable.
C
      ELSE IF (WORD.EQ.'EXEC') THEN
         CALL READA(EXEC)
C
C  GMIN executable.
C
      ELSE IF (WORD.EQ.'EXECGMIN') THEN
         CALL READA(EXECGMIN)
C
C  Write the coordinates of minimum WHICHMIN to file extractedmin and stop.
C
      ELSE IF (WORD.EQ.'EXTRACTMIN') THEN
         CALL READI(WHICHMIN)
         EXTRACTMINT=.TRUE.
      ELSE IF (WORD.EQ.'EXTRACTMINFILE') THEN
         EXTRACTMINT=.TRUE.
         EXTRACTMINFILET=.TRUE.
C
C  Write the coordinates of minimum WHICHTS to file extractedts and stop.
C
      ELSE IF (WORD.EQ.'EXTRACTTS') THEN
         EXTRACTTST=.TRUE.
         CALL READI(WHICHTS)
      ELSE IF (WORD.EQ.'EXTRACTTSFILE') THEN
         EXTRACTTST=.TRUE.
         EXTRACTTSFILET=.TRUE.
C
C Choose connection pairs based on free energy barriers between regrouped minima
C and product minima.
C
      ELSE IF (WORD.EQ.'FREEPAIRS') THEN
         FREEPAIRT=.TRUE.
         CALL READF(REGROUPFREETHRESH)
         CALL READF(FREETHRESH)
         CALL READF(EINC)
         IF (EINC.LE.0.0D0) THEN
            PRINT '(2(A,G20.10))','keywords> WARNING EINC for FREEPAIRS reset from ',EINC,' to ',1.0D0
            EINC=1.0D0
         ENDIF
C
C  Frozen atoms (to adjust zeros when reading frequencies)
C  and needed for freezing atoms is tssearch.
C
      ELSE IF (WORD.EQ.'FREEZE') THEN
         IF (NATOMS.LE.0) THEN
            PRINT '(A,I6,A)','keywords> ERROR - NATOMS=',NATOMS,' NATOMS keyword must preceed FREEZE'
            STOP
         ENDIF
         IF (.NOT.ALLOCATED(FROZEN)) THEN
            ALLOCATE(FROZEN(NATOMS))
            DO J1=1,NATOMS
               FROZEN(J1)=.FALSE.
            ENDDO
         ENDIF
         FREEZE=.TRUE.
         DO J1=1,NITEMS-1
            NFREEZE=NFREEZE+1
            CALL READI(NDUMMY)
            FROZEN(NDUMMY)=.TRUE.
!           PRINT '(A,I6,A,I6)','freezing atom ',NDUMMY,' NFREEZE=',NFREEZE
         ENDDO
         IF (PERMDIST.OR.LPERMDIST) THEN
            NDUMMY=0
            DO J1=1,NPERMGROUP
               DO J2=1,NPERMSIZE(J1)
                  IF (FROZEN(PERMGROUP(NDUMMY+J2))) THEN
                     PRINT '(A,I8,A)','keywords> ERROR atom ',PERMGROUP(NDUMMY+J2),' cannot be frozen and permuted'
                     STOP
                  ENDIF
               ENDDO
               NDUMMY=NDUMMY+NPERMSIZE(J1)
            ENDDO
         ENDIF       
      ELSE IF (WORD.EQ.'FREEZERANGE') THEN
         FREEZE=.TRUE.
         IF (.NOT.ALLOCATED(FROZEN)) THEN
            ALLOCATE(FROZEN(NATOMS))
            DO J1=1,NATOMS
               FROZEN(J1)=.FALSE.
            ENDDO
         ENDIF
         IF (NITEMS.GT.1 .and. NITEMS .LT. 4) THEN
            CALL READI(NDUMMY)
            J1=NDUMMY
            CALL READI(NDUMMY)
            J2=NDUMMY
            DO J3=J1,J2
               NFREEZE=NFREEZE+1
               FROZEN(J3)=.TRUE.
            ENDDO
         ELSE
           WRITE (*,'(A)') ' ERROR: FREEZERANGE specified incorrectly'
         ENDIF
         IF (PERMDIST.OR.LPERMDIST) THEN
            NDUMMY=0
            DO J1=1,NPERMGROUP
               DO J2=1,NPERMSIZE(J1)
                  IF (FROZEN(PERMGROUP(NDUMMY+J2))) THEN
                     PRINT '(A,I8,A)','keywords> ERROR atom',PERMGROUP(NDUMMY+J2),' cannot be frozen and permuted'
                     STOP
                  ENDIF
               ENDDO
               NDUMMY=NDUMMY+NPERMSIZE(J1)
            ENDDO
         ENDIF
C
C Friction coefficient for Berezhkovsii, Pollak, Zitserman formulation
C JCP, 97, 2422, 1992
C
      ELSE IF (WORD.EQ.'FRICTION') THEN
         FRICTIONT=.TRUE.
         IMFRQT=.TRUE.
         CALL READF(GAMMAFRICTION)

      ELSE IF (WORD.EQ.'FROMLOWEST') THEN
         FROMLOWESTT=.TRUE.
         IF (NITEMS.GT.1) THEN 
            CALL READA(DUMMYSTRING)
            IF (DUMMYSTRING.EQ.'NOLABELS') THEN 
               NOLABELST=.TRUE.
               PRINT *,'keywords> GMIN lowest file not expected to contain atom labels'
            ELSE
               STOP ' ERROR: Invalid argument to FROMLOWEST keyword in pathdata'
            ENDIF
         ENDIF
C
C  Distance criterion for distinguishing stationary points
C
      ELSE IF (WORD.EQ.'GEOMDIFFTOL') THEN
         CALL READF(GEOMDIFFTOL)

      ELSE IF (WORD.EQ.'GETMINFRQS') THEN
         GETMINFRQST=.TRUE.
         IF (GETTSFRQST) THEN
            PRINT *,'keywords> only GETMINRFQS or GETTSFRQS can be set'
            STOP
         ENDIF
         PRINT *,'keywords> get min frqs'

      ELSE IF (WORD.EQ.'GETTSFRQS') THEN
         GETTSFRQST=.TRUE.
         IF (GETMINFRQST) THEN
            PRINT *,'keywords> only GETMINRFQS or GETTSFRQS can be set'
            STOP
         ENDIF
         PRINT *,'keywords> get ts frqs'

C
C  Specify graph transformation rate calculation. NCONNMIN is the connectivity
C  at which (and below) minima are removed. GTINT is the interval for performing
C  the analysis during a DPS calculation. IF GTINT=1 we do the analysis every cycle,
C  consisting of NCPU jobs, if GTINT=2 we do the analysis every other cycle. 
C  GTINT is zero by default, in which case we only do the analysis if NATTEMPTS=0.
C  GT is the original DJW implementation where we renormalise out all intervening
C  miinima. The resulting branching probablities and waiting times do not allow
C  for return to the starting minimum - they are like the rejection-free BKL approach.
C
      ELSE IF (WORD.EQ.'GT') THEN
         GTT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
         IF (NITEMS.GT.2) CALL READI(GTINT)
C
C  Semen Trygubenko, last modification Thu Mar 15 12:52:36 GMT 2007
C
C  A different implementation of GT method, as described in our GT JCP paper.
C  Input format: {SGT|DGT|SDGT} [DisconnectSources [AltPbb [Rescale [Normalise]]]]
C  Example input line: SDGT F T F F
C 
C  If set to true DisconnectSources instructs to disconnect sources, yielding rate
C  that theoretically corresponds to the one obtained stochastically with KMC. Disconnection
C  of the sources is currently implemented in DGT part of GT. For SGT runs with DisconnectSources=T
C  a call is made to DGT once all SGT work is done. DisconnectSources=T can be memory hungry.
C
C  If set to true, GT2Sparse instructs to perform GT analysis using
C  sparse-optimised algorithms and data structures.
C  
C  GT2Switch=T instructs the program to change the data structures and algorithms
C  from sparse to dense when the density requirement (can be adjusted by GT2RSWITCH) is
C  met.
C  
C  GT2AltPbb, when set to true, triggers the evaluation of Pbb sums using the algorithm
C  suggested by David to maintain precision. Alternative remedy is provided with
C  GT2Rescale which I have devised to stop the errors propagation. Both of these
C  should double the execution time, and are recognized by both
C  SGTDetachNode and DGTDetachNode routines. NB: Disconnect routine, another place where
C  roundoff errors can breed like rabbits, does not recognize either of these as of
C  Thu Mar 15 12:45:55 GMT 2007.
C  
C  GT2Normalise=T will instruct GT2input.f90 to normalize the branching probabilities when
C  obtained from PATHSAMPLE. This is purely for debugging purposes, as they always should be.
C 
      ELSE IF (word=='SGT'.or.word=='DGT'.or.word=='SDGT'.or.word=='GT2') then
         GT2T=.TRUE.
         if (word=='SGT') then
              GT2Sparse=.True.
              GT2Switch=.False.
         else if (word=='DGT') then
              GT2Sparse=.False.
              GT2Switch=.False.
         else if (word=='SDGT') then
              GT2Sparse=.True.
              GT2Switch=.True.
         else
              print '(A)','Please specify SGT, DGT or SDGT. GT2 keyword is now obsolete';stop
         endif
         PRINT *,'keywords> GT2Sparse=',GT2Sparse
         PRINT *,'keywords> GT2Switch=',GT2Switch
    
         IF (NITEMS>1) THEN
              CALL READA(LOYNO)
              IF (LOYNO=='T'.or.loyno=='F') then
                   READ(LOYNO,'(l1)') GT2DisconnectSources
                   PRINT '(1x,A,L5)','keywords> GT2DisconnectSources=',GT2DisconnectSources
              ELSE
                   PRINT '(1x,A,L5)','keywords> Invalid after "GT2" keyword'; stop
              ENDIF
         ENDIF
         IF (NITEMS>2) THEN
              CALL READA(LOYNO)
              IF (LOYNO=='T'.or.loyno=='F') then
                   READ(LOYNO,'(l1)') GT2AltPbb
                   PRINT '(1x,A,L5)','keywords> GT2AltPbb=',GT2AltPbb
              ELSE
                   PRINT '(1x,A,L5)','keywords> Invalid after "GT2" keyword'; stop
              ENDIF
         ENDIF
         IF (NITEMS>3) THEN
              CALL READA(LOYNO)
              IF (LOYNO=='T'.or.loyno=='F') then
                   READ(LOYNO,'(l1)') GT2Rescale
                   PRINT *,'keywords> GT2Rescale=',GT2Rescale
              ELSE
                   PRINT *,'keywords> Invalid after "GT2" keyword'; stop
              ENDIF
         ENDIF
         IF (NITEMS>4) THEN
              CALL READA(LOYNO)
              IF (LOYNO=='T'.or.loyno=='F') then
                   READ(LOYNO,'(l1)') GT2Normalise
                   PRINT *,'keywords> GT2Normalise=',GT2Normalise
              ELSE
                   PRINT *,'keywords> Invalid after "GT2" keyword'; stop
              ENDIF
         ENDIF
         IF (NITEMS>5) THEN
              print '(1x,a)','Input error: more than 4 entities following SGT, DGT or SDGT keyword!';stop
         ENDIF
C
C Switching ratio. Specifies when to switch from SGT to DGT.
C See GT2.f90 for details.
C
      ELSE IF (WORD=='GT2RSWITCH') then
         CALL READF(GT2RSWITCH)
C
C Is used to establish whether a node is a dead end or not.
C See GT2.f90 for details.
C
      ELSE IF (WORD=='GT2PTOL') then
         CALL READF(GT2PTOL)
C
C Takes into accound the mirror symmetry for the Thomson problem
C hk286
C
      ELSE IF (WORD=='GTHOMSON') then
         GTHOMSONT = .TRUE.
         CALL READI(GTHOMMET)

C
C Do not ignore previously searched pairs for pruning 
C

      ELSE IF (WORD=='IGNOREPAIRS') THEN
         CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'TRUE') PAIRSIGNORET=.TRUE.
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'FALSE') PAIRSIGNORET=.FALSE.
C
C Prints the negative eigenvalue of each TS in ts.data as final (ninth) column
C
      ELSE IF (WORD=='IMFRQ') THEN
         IMFRQT=.TRUE.
C
C  Set initial distance for all pairs of minima for use with DIJINITSTART/CONT
C
      ELSE IF (WORD.EQ.'INITIALDISTANCE') THEN
         INITIALDIST=.TRUE.
         CALL READF(DISBOUND)
C
C  Use constraint potential for interpolation as a connection metric (instead of distance).
C
      ELSE IF (WORD.EQ.'INTCONSTRAINT') THEN
         INTCONSTRAINTT=.TRUE.
         IF (NITEMS.GT.1) CALL READF(INTCONSTRAINTTOL)
         IF (NITEMS.GT.2) CALL READF(INTCONSTRAINTDEL)
         IF (NITEMS.GT.3) CALL READF(INTCONSTRAINTREP)
         IF (NITEMS.GT.4) CALL READF(INTCONSTRAINREPCUT)
         IF (NITEMS.GT.5) CALL READI(INTCONSEP)
         IF (NITEMS.GT.6) CALL READI(INTREPSEP)
         INTERPCOSTFUNCTION=.TRUE.
      ELSE IF (WORD.EQ.'INTFREEZE') THEN
         INTFREEZET=.TRUE.
         IF (NITEMS.GT.1) CALL READF(INTFREEZETOL)
C
C  Use interpolation potential for LJ.
C
      ELSE IF (WORD.EQ.'INTLJ') THEN
         INTLJT=.TRUE. 
         IF (NITEMS.GT.1) CALL READF(INTLJDEL)
         IF (NITEMS.GT.2) CALL READF(INTLJEPS)
         INTERPCOSTFUNCTION=.TRUE.
C
C  Inertia difference criterion - no longer used for distinguishing stationary points!
C
      ELSE IF (WORD.EQ.'ITOL') THEN
         CALL READF(IDIFFTOL)
C
C  Number of jobs to run per node for a distributed memory architecture.
C  JPN should always be one for all recent cluster implementations,
C  unless you really want to run multiple jobs per core. Note that the
C  nodes referred to here are really cores.
C
      ELSE IF (WORD.EQ.'JOBSPERNODE') THEN
          WRITE(*,'(A)') 'keywords> ERROR: JOBSPERNODE is deprecated. Please generate a valid nodes.info file and use PBS instead'
          STOP
C         PBST=.TRUE.
C         IF (NITEMS.GT.1) CALL READI(JPN)
C         CALL GETNODES(NCPU)
      ELSE IF (WORD.EQ.'KMC') THEN
         KMCT=.TRUE.
         NOPOINTS=.TRUE.
         IF (NITEMS.GT.1) CALL READF(NKMCCYCLES)
         IF (NITEMS.GT.2) CALL READF(PAIRTHRESH)
         IF (NITEMS.GT.3) CALL READI(NCONNMIN)
C
C  Set KMC parameters: number of KMC runs for averages and value of PAIRTHRESH.
C  If the product of branching probabilities p12*p21 > PAIRTHRESH then we
C  renormalise this pair.
C
      ELSE IF (WORD.EQ.'KMCCOMMIT') THEN
         NOPOINTS=.TRUE.
         KMCCOMMITT=.TRUE.
         IF (NITEMS.GT.1) CALL READF(NKMCCYCLES)
         IF (NITEMS.GT.2) CALL READF(MAXBREAK)
         IF (NITEMS.GT.3) CALL READF(PABCONV)
         IF (NITEMS.GT.4) CALL READF(OMEGA)
         IF (NITEMS.GT.5) CALL READF(PAIRTHRESH)
         IF (NITEMS.GT.6) CALL READI(NCONNMIN)
C
C  k-th shortest paths analysis for each B (or A) minimum.
C
      ELSE IF (WORD.EQ.'KSHORTESTPATHS') THEN
         KSHORTESTPATHST=.TRUE.
         CALL READI(NPATHS)
         CALL READI(NCONNMIN)
         IF (NITEMS.GT.3) THEN 
            CALL READA(LOYNO)
            IF (LOYNO == 'T') KSHORT_FULL_PRINTT = .TRUE.
         ENDIF
C
C  Whether to read extra curvatures from min.data.info files in DUMMYTS runs
C
      ELSE IF (WORD.EQ.'LOWESTFRQ') THEN
         LOWESTFRQT=.TRUE.

      ELSE IF (WORD.EQ.'MKTRAP') THEN
         MKTRAPT=.TRUE.
         ZSYM='BE'
!
! MLLJAT3 for ML LJAT3 time series landscapes
!
      ELSE IF (WORD.EQ.'MLLJAT3') THEN
         MLLJAT3=.TRUE.
!
! Three layer neural network (multilayer perceptron) with
! MLPIN inputs (columns per data item)
! MLPOUT outputs
! MLPHIDDEN hidden nodes
! MLPDATA data lines in MLPdata file (last column MLPIN+1 for correct outputs, numbered one to MLPOUT)
! MLPLAMBDA coefficient for regularisation
!
      ELSE IF ((WORD.EQ.'MLPVB3').OR.(WORD.EQ.'MLPVB3NN')) THEN
         MLPB3T=.TRUE.
         ZSYM='  '
         CALL READI(MLPIN)      ! number of inputs (data items after outcome)
         CALL READI(MLPSTART) ! starting position in data list, not counting outcome
         CALL READI(MLPHIDDEN)
         CALL READI(MLPOUT)
         CALL READI(MLPDATA)
         IF (NITEMS.GT.5) CALL READF(MLPLAMBDA)
         IF ((WORD.EQ.'MLPVB3NN').AND.(NITEMS.GT.6)) CALL READI(MLPNEIGH)
         IF (MLPNEIGH.EQ.0) THEN
            WRITE(*,'(A)') 'keywords> *** ERROR cannot have zero nearest neighbours, check odata'
            STOP
         ENDIF
         IF (WORD.EQ.'MLPVB3NN') MLPVB3NNT=.TRUE.

         WRITE(*,'(A,5I8,G20.10)') ' keywords> MLPVB3 vector bias nodes and Nin, Ninstart, Nhidden, Nout=',
     &                                MLPIN,MLPSTART,MLPHIDDEN,MLPOUT
         NMLP=MLPHIDDEN*(MLPIN+MLPOUT)+MLPHIDDEN+MLPOUT
         NATOMS=NMLP
         NOPT=NMLP

         IF (WORD.EQ.'MLPVB3NN') THEN
            IF (.NOT.ALLOCATED(FROZEN)) THEN
               ALLOCATE(FROZEN(NATOMS))
               DO J1=1,NATOMS
                  FROZEN(J1)=.FALSE.
               ENDDO
            ENDIF

            ALLOCATE( MLPDISTHI(MLPIN), MLPDISTHO(MLPOUT), MLPINDEXI(MLPIN), MLPINDEXO(MLPOUT))
            WRITE(*,'(A)') 'Original nearest-neighbour fomulation:'
            FREEZE=.TRUE.
            NFREEZE=MLPHIDDEN*(MLPIN+MLPOUT)-2*MLPHIDDEN
            FROZEN(1:MLPHIDDEN*(MLPIN+MLPOUT))=.TRUE.
            DO J1=1,MLPHIDDEN
               J2=NINT(1.0D0*(MLPHIDDEN+J1*(MLPIN-1)-MLPIN)/(MLPHIDDEN-1)) ! unfrozen weight for hidden node J1 to input
               J3=(J1-1)*MLPIN+J2
               FROZEN(J3)=.FALSE.
               PRINT '(A,I10,A,I10,A,I10)','keywords> Unfrozen weight ',J3,' input ',J2,' to hidden node ',J1
               J2=NINT(1.0D0*(MLPHIDDEN+J1*(MLPOUT-1)-MLPOUT)/(MLPHIDDEN-1)) ! unfrozen weight for hidden node J1 to output
               J3=MLPHIDDEN*MLPIN+(J2-1)*MLPHIDDEN+J1
               PRINT '(A,I10,A,I10,A,I10)','keywords> Unfrozen weight ',J3,' hidden node ',J1,' to output ',J2
               FROZEN(J3)=.FALSE.
            ENDDO
            WRITE(*,'(A,I6,A)') 'keywords> New nearest-neighbour formulation with ',MLPNEIGH,' neighbours'
            NFREEZE=MLPHIDDEN*(MLPIN+MLPOUT)
            FROZEN(1:MLPHIDDEN*(MLPIN+MLPOUT))=.TRUE.
            DO J1=1,MLPHIDDEN
                  !
! Distances from hidden J1 to all input nodes for
! w^2_{J1 K1} at (J1-1)*MLPIN+K1 up to MLPHIDDEN*MLPIN
!
               DO K1=1,MLPIN
                  MLPINDEXI(K1)=K1
                  MLPDISTHI(K1)=( (J1-1.0D0)/(MLPHIDDEN-1.0D0) - (K1-1.0D0)/(MLPIN-1.0D0) )**2 - K1*1.0D-6 ! to break degeneracy
               ENDDO
               CALL SORT4(MLPIN,MLPIN,MLPDISTHI,MLPINDEXI)
               DO J2=1,MIN(MLPNEIGH,MLPIN)
                  WRITE(*,'(A,I8,A,I8,A,I8,A,G20.10)') 'hidden ',J1,' input neighbour ',J2,' is ',MLPINDEXI(J2),' distance ',
     &                                                      MLPDISTHI(J2)
                  J3=(J1-1)*MLPIN+MLPINDEXI(J2)
                  FROZEN(J3)=.FALSE.
                  NFREEZE=NFREEZE-1
               ENDDO
!
! Distances from hidden J1 to all output nodes for
! w^1_{I1 J1} at MLPHIDDEN*MLPIN + (I1-1)*MLPHIDDEN+J1 up to MLPHIDDEN*MLPIN + MLPOUT*MLPHIDDEN
!
               DO II1=1,MLPOUT
                  MLPINDEXO(II1)=II1
                  MLPDISTHO(II1)=( (J1-1.0D0)/(MLPHIDDEN-1.0D0) - (II1-1.0D0)/(MLPOUT-1.0D0) )**2 - II1*1.0D-6 ! to break degeneracy
               ENDDO
               CALL SORT4(MLPOUT,MLPOUT,MLPDISTHO,MLPINDEXO)
               DO J2=1,MIN(MLPNEIGH,MLPOUT)
                  WRITE(*,'(A,I8,A,I8,A,I8,A,G20.10)') 'hidden ',J1,' output neighbour ',J2,' is ',MLPINDEXO(J2),
     &               ' distance ',MLPDISTHO(J2)
                  J3=MLPHIDDEN*MLPIN+(MLPINDEXO(J2)-1)*MLPHIDDEN+J1
                  FROZEN(J3)=.FALSE.
                  NFREEZE=NFREEZE-1
               ENDDO
            ENDDO
            DEALLOCATE( MLPDISTHI, MLPDISTHO, MLPINDEXI, MLPINDEXO)

         ENDIF

      ELSE IF (WORD.EQ.'MLPB3NEW') THEN
         MLPB3T=.TRUE.
         ZSYM='  '
         CALL READI(MLPIN)      ! number of inputs (data items after outcome)
         CALL READI(MLPSTART) ! starting position in data list, not counting outcome
         CALL READI(MLPHIDDEN)
         CALL READI(MLPOUT)
!        CALL READI(MLPDATA)
!        IF (NITEMS.GT.5) CALL READF(MLPLAMBDA)
         WRITE(*,'(A,5I8,G20.10)') ' keywords> MLP3 new potential bias nodes and Nin, Ninstart, Nhidden, Nout=',
     &                                MLPIN,MLPSTART,MLPHIDDEN,MLPOUT
         NMLP=MLPHIDDEN*(MLPIN+MLPOUT)+1
         NATOMS=NMLP
         NOPT=NMLP
      ELSE IF ((WORD.EQ.'MLP3').OR.(WORD.EQ.'MLPB3')) THEN
         MLP3T=.TRUE.
         ZSYM='  '
         CALL READI(MLPIN)
         CALL READI(MLPHIDDEN)
         CALL READI(MLPOUT)
!        CALL READI(MLPDATA)
!        IF (NITEMS.GT.5) CALL READF(MLPLAMBDA)
         IF (WORD.EQ.'MLPB3') THEN
            MLPB3T=.TRUE.
            WRITE(*,'(A,4I8,G20.10)') 'MLP3 potential with bias nodes and Nin, Nhidden, Nout=',
     &                                   MLPIN,MLPHIDDEN,MLPOUT
            NMLP=MLPHIDDEN*(MLPIN+MLPOUT)+1
         ELSE
            WRITE(*,'(A,4I8,G20.10)') 'MLP3 potential with Nin, Nhidden, Nout=',
     &                                   MLPIN,MLPHIDDEN,MLPOUT
            NMLP=MLPHIDDEN*(MLPIN+MLPOUT)
         ENDIF
         NATOMS=NMLP
         NOPT=NMLP
      ELSE IF (WORD.EQ.'MACHINE') THEN
         MACHINE=.TRUE.

C
C Macrocycle
C This adds cyclic isomers to the MINPERMDIST
C Used for cyclic peptides with repeating sequences (e.g. cyclo-[GlyProGlyPro])
C
      ELSE IF (WORD.EQ.'MACROCYCLE') THEN
         MACROCYCLET=.TRUE.
         IF (NITEMS.GT.1) THEN
            CALL READI(MCYCLEREPEATS)
         ENDIF
         MCYCLEPERIOD=NATOMS/MCYCLEREPEATS
         IF (MCYCLEPERIOD*MCYCLEREPEATS.NE.NATOMS) THEN
            PRINT('(A,I5,A,I5,A)'),"Warning: ring of size ",NATOMS," cannot contain ",MCYCLEREPEATS, " repeat units."
         ELSE
            PRINT('(A,I5,A,I5,A)'),"keyword> Macrocycle with ",MCYCLEREPEATS," units, each comprising ",MCYCLEPERIOD," atoms."
         ENDIF
! sf344> macroion model
      ELSE IF (WORD.EQ.'MACROION') THEN
            MACROIONT=.TRUE.
 
C
C  MAXBARRIER requires both sides to be greater than MAXBARRIER to discard.
C
      ELSE IF (WORD.EQ.'MAXBARRIER') THEN
         CALL READF(MAXBARRIER)
C
C  The maximum number of constraints to use in the constrained potential.
C  The default is 3.
C
      ELSE IF (WORD.EQ.'MAXCON') THEN
         CALL READI(MAXCONUSE)
C
C  MAXDOWNBARRIER checks the downhill barrier in Dijkstra.
C
      ELSE IF (WORD.EQ.'MAXDOWNBARRIER') THEN
         CALL READF(MAXDOWNBARRIER)
C
C  TSTHRESH discards transition states above the specfied threshold. May be useful
C  for producing a better initial path and excluding CHARMM transition states with
C  silly energies from the database. MAXTSENERGY does the same as TSTHRESH keyword 
C  for compatability with OPTIM.
C
      ELSE IF (WORD.EQ.'MAXTSENERGY') THEN
         CALL READF(TSTHRESH)
C
C  Add the minima and transition= states in path.info.<PATHNAME> and
C  output.<PATHNAME> to an existing database. The end points are NOT assumed to
C  belong to the A and B sets.
C
      ELSE IF (WORD.EQ.'MERGEDB') THEN
         MERGEDBT=.TRUE.
         CALL READA(PATHNAME)
C
C  Calculate microcanonical thermodynamic properties.
C
      ELSE IF (WORD.EQ.'MICROTHERM') THEN
         MICROTHERMT=.TRUE.
         MICROEMIN=1.0D0
         MICROEMAX=2.0D0
         MICROEINC=0.1D0
         MICROT=1.0D0
         CALL READF(MICROEMIN)
         CALL READF(MICROEMAX)
         CALL READF(MICROEINC)
         CALL READF(MICROT)
C
C  MINBARRIER - set this in order to exclude from the analysis bad TSs that are lower in
C  potential energy than either of the connected minima (see checkTS.f90).
C  NB use with caution - TSs may legitimately be lower in FREE energy.
C
      ELSE IF (WORD.EQ.'MINBARRIER') THEN
         CALL READF(MINBARRIER)

C
C  MINGAP - exclude missing connections if their separation is smaller than MINGAPINP,
C  or alternatively if ratio is set smaller than MINGAPINP*(max. pair dist.), with
C  DIJINIT to focus on large separation rather than short connections
C
      ELSE IF (WORD.EQ.'MINGAP') THEN
         MINGAPT=.TRUE.
         CALL READF(MINGAPINP)
         IF (NITEMS.GT.2) CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'RATIO') MINGAPRATIOT=.TRUE.


      ELSE IF (WORD.EQ.'MSSTOCK') THEN
         MSSTOCKT = .TRUE.
         CALL READI(NRBSITES)
         ALLOCATE(RBSITE(NRBSITES,3))
         CALL DEFMULTSTOCK()

C
C  Number of atoms - essential unless RIGIDBODIES/RBAA keyword is present.
C
      ELSE IF (WORD.EQ.'NATOMS') THEN
         CALL READI(NATOMS)
         NOPT=3*NATOMS
C
C  Semen Trygubenko, Thu Mar 15 16:21:41 GMT 2007
C  Minimum number of connections a minimum ought to have to be included in the
C  rate calculation or regrouping schemes. Can be set directly for some
C  keywords, but not others, so this keyword is provided just to make sure
C  it casn be set for all methods.
C
      ELSE IF (WORD.EQ.'NCONNMIN') THEN
         if (NITEMS.GT.1) then
              CALL READI(NCONNMIN)
         else
              print *, 'keywords> Usage: RateNConnMin <integer>'; stop
         endif
C
C  Extend database via parallel single-ended transition state searches.
C
      ELSE IF (WORD.EQ.'NEWCONNECTIONS') THEN
         NEWCONNECTIONST=.TRUE.
         CALL READI(CONNECTIONS)
         IF (NITEMS.GT.2) CALL READI(MAXTSATTEMPTS)
         IF (NITEMS.GT.3) CALL READI(CONNMINSTART)
C
C  Specify new graph transformation rate calculation. NCONNMIN is the connectivity
C  at which (and below) minima are removed. 
C  NGT is the new implementation where we renormalise out all intervening
C  miinima but allow return to the starting state. This should give committor
C  probabilities if we remove all the I minima, unlike GT and GT2.
C
      ELSE IF (WORD.EQ.'NGT') THEN
         NGTT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
         IF (NITEMS.GT.2) THEN
              CALL READA(LOYNO)
              IF (LOYNO=='T'.OR.LOYNO=='F') then
                 READ(LOYNO,'(l1)') NGTDISCONNECTALL
                 IF (NGTDISCONNECTALL) PRINT '(A)','keywords> NGT will calculate kSS, kNSS and kKMC'
                 IF (.NOT.NGTDISCONNECTALL) PRINT '(A)','keywords> NGT will calculate kSS and kNSS but not kKMC'
              ELSE
                 PRINT '(1X,A,L5)','keywords> Invalid after "NGT" keyword'; STOP
              ENDIF
         ENDIF
         IF (NITEMS.GT.3) CALL READF(NGTSWITCH)
         IF (NITEMS.GT.4) CALL READI(NGTSIZE)
         PRINT '(A,F12.4)','keywords> NGT will switch to dense renormalisation scheme at threshold ',NGTSWITCH
         PRINT '(A,I6)','keywords> NGT maximum square matrix size in dense phase=',NGTSIZE
         IF (NITEMS.GT.5) CALL READF(NGTCRSWITCH) ! N.B. threshold in giga-bytes
         PRINT '(A,F12.4)','keywords> NGT will use compressed row storage scheme beyond threshold ',NGTCRSWITCH
      ELSE IF (WORD.EQ.'NIMET') THEN
         NIMET=.TRUE.

      ELSE IF (WORD.EQ.'NIHEAM7') THEN
         NIHEAM7T=.TRUE.

      ELSE IF (WORD.EQ.'NIH2LEPS') THEN
         NIH2LEPST=.TRUE.

      ELSE IF (WORD.EQ.'NINTS') THEN
         CALL READI(NINTS)
C
C  If NOFRQS is specified then frequencies are assumed to be absent from path.info files.
C
      ELSE IF (WORD.EQ.'NOFRQS') THEN
         NOFRQS=.TRUE.
C
C  If NOINVERSION is specified then the inversion operation is not allowed
C  in minpermdist. Needed for some twod examples.
C
      ELSE IF (WORD.EQ.'NOINVERSION') THEN
         NOINVERSION=.TRUE.
C
C  If NOPOINTS is specified then setup should not try to read the min.points or ts.points
C  files. Should be the default for post-DPS database kinetics analysis such as KMC.
C
      ELSE IF (WORD.EQ.'NOPOINTS') THEN
         NOPOINTS=.TRUE.
C
C  If NOPOINTGROUP is specified then all point group orders are set to unity.
C  Needed for addressable potentials where we count the permutational isomers explicitly.
C
      ELSE IF (WORD.EQ.'NOPOINTGROUP') THEN
         NOPOINTGROUPT=.TRUE.
!
! Forbid overall translation and rotation in distance/alignment calculations.
!
      ELSE IF (WORD.EQ.'NOTRANSROT') THEN
            NOTRANSROTT=.TRUE.

      ELSE IF (WORD.EQ.'NTIP') THEN

         CALL READI(TIPID)

         IF (TIPID == 4) THEN
            NRBSITES = 4
         ELSE
            PRINT *, 'TIPID NOT EQUAL TO 4 NOT YET DEFINED'
            STOP
         ENDIF

         NTIPT    = .TRUE.
         ALLOCATE(RBSITE(NRBSITES,3))
         ALLOCATE(SITEMASS(NRBSITES))

         IF (TIPID == 4) THEN
            CALL DEFTIP4(SITEMASS)
         ENDIF

C
C  Specify Oh supercell to allow box symmetries in permutational alignment.
C
      ELSE IF (WORD.EQ.'OHCELL') THEN
         OHCELLT=.TRUE.
         WRITE(*,'(A)') 'Octahedral supercell specfied'
C
C  Just do a single regrouping pass, not recursive
C
      ELSE IF (WORD.EQ.'ONEREGROUP') THEN
         ONEREGROUPT=.TRUE.

C
C  OPEP potential used, need to make sure we get the right files to
C  OPTIM
C
      ELSE IF(WORD.EQ.'OPEP') THEN
         OPEPT=.TRUE.


C
C  Read in an order parameter threshold, ORDERPARAM, that tells us when we are in the other 
C  phase for a DOALLMIN run.
C
      ELSE IF (WORD.EQ.'ORDER') THEN
         CALL READF(ORDERPARAM)
C
C  Set the frequency with which the list of pairs for future searches is recreated.
C  Currently works for default search type based on Pfold difference.
C

      ELSE IF (WORD.EQ.'PAHA') THEN

         CALL READI(PAHID)

         IF (PAHID == 1) THEN
            NRBSITES = 12
         ENDIF

         PAHAT    = .TRUE.
         ALLOCATE(RBSITE(NRBSITES,3))
         ALLOCATE(SITEMASS(NRBSITES))
         SITEMASS(:) = 1.D0

         IF (PAHID == 1) THEN
            CALL DEFBENZENE()
         ENDIF


      ELSE IF (WORD.EQ.'PAIRLIST') THEN
         CALL READI(NPAIRFRQ)

      ELSE IF (WORD.EQ.'PAP') THEN

         CALL READI(PAPID)
         CALL READF(PAPALP)

         IF (PAPID == 1) THEN
            NRBSITES = 7
         ELSEIF (PAPID == 2) THEN
            NRBSITES = 5
         ELSEIF (PAPID == 3) THEN
            NRBSITES = 3
         ELSEIF (PAPID == 4) THEN
            NRBSITES = 5
         ELSE
            PRINT *, 'PAPID equal to ', PAPID, ' not defined.'
            STOP
         ENDIF

         PAPT    = .TRUE.
         ALLOCATE(RBSITE(NRBSITES,3))
         ALLOCATE(SITEMASS(NRBSITES))
         SITEMASS(:) = 1.D0

         CALL DEFPAP()

      ELSE IF (WORD.EQ.'PATCHYD') THEN
         PATCHYDT = .TRUE.
         CALL READI(NRBSITES)
         IF (NRBSITES .NE. 4) THEN
            PRINT *, 'NRBSITES has to be 4'
            STOP 
         ENDIF
         ALLOCATE(RBSITE(NRBSITES,3))
         CALL DEFPATCHYD()
C
C  Number of jobs to run per node for a distributed memory architecture.
C  JPN should always be one for all recent cluster implementations,
C  unless you really want to run multiple jobs per core. Note that the
C  nodes referred to here are really cores.
C
      ELSE IF (WORD.EQ.'PBS') THEN
         PBST=.TRUE.
         IF (NITEMS.GT.1) CALL READI(JPN)
         CALL GETNODES(NCPU)
!
!  Whether to optimise the permutational isomers in assessing optimal alignment.
!  For PERMDIST all minimum distances will be minimsed with respect to the
!  specified permutations. For PERMISOMER we only check for permutational isomers
!  if the energy difference is below EDIFFTOL. This should save a great deal
!  of time for large systems containing many equivalent atoms, although the
!  distance won't be minimised with repect to permutations for inequivalent minima.
!
      ELSE IF ((WORD.EQ.'PERMDIST').OR.(WORD.EQ.'PERMISOMER').OR.(WORD.EQ.'LPERMDIST')) THEN
         IF (NATOMS.LE.0) THEN
            PRINT '(A,I6,A)','keywords> ERROR - NATOMS=',NATOMS,' NATOMS keyword must preceed PERMDIST'
            STOP
         ENDIF
         IF (WORD.EQ.'PERMDIST') THEN
            PERMDIST=.TRUE.
            IF (NITEMS.GT.1) CALL READF(ORBITTOL)
            IF (NITEMS.GT.2) CALL READI(MAXNSETS)
            PRINT '(A,F15.5)','keywords> Distance tolerance for identifying atoms in the same orbit=',ORBITTOL
            PRINT '(A,I3)',' keywords> Maximum number of secondary sets in perm.allow file=', MAXNSETS
         ENDIF
         IF (WORD.EQ.'PERMISOMER') PERMISOMER=.TRUE.
         IF (WORD.EQ.'LPERMDIST') THEN
            LPERMDIST=.TRUE.
            IF (NITEMS.GT.1) CALL READI(LOCALPERMNEIGH)
            IF (NITEMS.GT.2) CALL READF(LOCALPERMCUT)
            IF (NITEMS.GT.3) CALL READF(LOCALPERMCUT2)
            IF (NITEMS.GT.4) CALL READF(ORBITTOL)
            PRINT '(A,F15.5)','keywords> Local permutational alignment: alignment threshold=',LOCALPERMCUT
            PRINT '(A,F15.5)','keywords> Local permutational alignment: alignment cutoff=   ',LOCALPERMCUT2
            PRINT '(A,F15.5)',' keyword> Distance tolerance for distinguishing atoms in the same orbit=',ORBITTOL
         ENDIF
         INQUIRE(FILE='perm.allow',EXIST=PERMFILE)
         ALLOCATE(NPERMSIZE(3*NATOMS),PERMGROUP(3*NATOMS),NSETS(3*NATOMS),SETS(NATOMS,MAXNSETS))
!
!  The above dimensions were fixed at NATOMS because:
!  (a) Atoms were not allowed to appear in more than one group.
!  (b) The maximum number of pair exchanges associated with a group is three.
!
! However, for flexible water models we need to exchange all waters,
! and we can exchange H's attached to the same O. The dimension required
! becomes 3*NATOMS
!
         IF (PERMFILE) THEN
            OPEN(UNIT=1,FILE='perm.allow',STATUS='OLD')
            READ(1,*) NPERMGROUP
            NDUMMY=1
            DO J1=1,NPERMGROUP
!
!  Sanity checks!
! 
               READ(1,*) NPERMSIZE(J1),NSETS(J1)
               IF (NSETS(J1).GT.MAXNSETS) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of secondary sets ',NSETS(J1),' is > ', MAXNSETS
                  STOP
               ENDIF
               IF (NDUMMY+NPERMSIZE(J1)-1.GT.3*NATOMS) THEN
                  PRINT '(2(A,I8))','keyword> ERROR - number of atoms to be permuted in all groups is > 3*number of atoms'
                  STOP
               ENDIF
!              READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SWAP1(PERMGROUP(J3),J2),J2=1,NSWAP(J1)),
!    &                                                            J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1)
               READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SETS(PERMGROUP(J3),J2),J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1),
     &                                                              J2=1,NSETS(J1))
               NDUMMY=NDUMMY+NPERMSIZE(J1)
            ENDDO
            CLOSE(1)
         ELSE
            NPERMGROUP=1 ! ALL ATOMS CAN BE PERMUTED - DEFAULT
            NPERMSIZE(1)=NATOMS ! ALL ATOMS CAN BE PERMUTED - DEFAULT
            DO J1=1,NATOMS
               PERMGROUP(J1)=J1
            ENDDO
         ENDIF
         PRINT '(A,I6)','keywords> Number of groups of permutable atoms=',NPERMGROUP
         NDUMMY=1
         DO J1=1,NPERMGROUP
            PRINT '(A,3(I6,A))','keywords> group ',J1,' contains ',NPERMSIZE(J1),' atoms with ',
     &                                                 NSETS(J1),' additional atom sets:'
            WRITE(*,'(22I6)',ADVANCE='NO') PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1) 
            IF (NSETS(J1).GT.0) THEN
               WRITE(*,'(A)',ADVANCE='NO') ' with '
               DO J2=1,NSETS(J1)
                  DO J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1
                     WRITE(*,'(I6)',ADVANCE='NO') SETS(PERMGROUP(J3),J2)
                     IF (J3.LT.NDUMMY+NPERMSIZE(J1)-1) WRITE(*,'(A3)',ADVANCE='NO') ' / '
                  ENDDO 
                  IF (J2.LT.NSETS(J1)) WRITE(*,'(A3)',ADVANCE='NO') ' ; '
               ENDDO 
            ENDIF
            PRINT *,' '
            NDUMMY=NDUMMY+NPERMSIZE(J1)
         ENDDO
!
!  Another check.
!
         IF (NFREEZE.GT.0) THEN
            NDUMMY=0
            DO J1=1,NPERMGROUP
               DO J2=1,NPERMSIZE(J1)
                  IF (FROZEN(PERMGROUP(NDUMMY+J2))) THEN
                     PRINT '(A,I8,A)','keywords> ERROR atom ',PERMGROUP(NDUMMY+J2),' cannot be frozen and permuted'
                     STOP
                  ENDIF
               ENDDO
               NDUMMY=NDUMMY+NPERMSIZE(J1)
            ENDDO
         ENDIF
C
C  Persistence analysis.
C
      ELSE IF (WORD.EQ.'PERSIST') THEN
         PERSISTT=.TRUE.
         CALL READF(EINC)
         CALL READF(EUNTRAPTHRESH)
         IF (NITEMS.GT.3) CALL READF(PEQTHRESH)
         IF (NITEMS.GT.4) CALL READF(PERTHRESH)
      ELSE IF (WORD.EQ.'PERSISTEXACT') THEN
         PERSISTAPPROXT=.FALSE.
         IF (NITEMS.GT.1) CALL READA(ALLCOMPS)
         IF (TRIM(ADJUSTL(ALLCOMPS)).EQ.'ALL') ALLCOMPONENTST=.TRUE.
C
C  Initial PERT parameter for geometry perturbations used in single-ended ts searches
C
      ELSE IF (WORD.EQ.'PERTURB') THEN
         CALL READF(PERTVALUE)
         PERTMIN=PERTVALUE/2.0D0
         PERTMAX=PERTVALUE*2.0D0
         IF (NITEMS.GT.2) CALL READF(PERTMIN)
         IF (NITEMS.GT.3) CALL READF(PERTMAX)
C
C  NPFOLD is the number of iterations per call to the global Pfold subroutine
C  and PFOLDINT is the frequency at which we call this GPFOLD calculation, 
C  i.e. PFOLDINT=1 means for every cycle, PFOLDINT=2 means every other cycle etc.
C  PFOLDCONV is the threshold on the largest fractional difference between Pfold 
C  values at consecutive iterations for possible early termination (before the 
C  NPFOLD iterations have been performed).
C
      ELSE IF (WORD.EQ.'PFOLD') THEN
         CALL READI(NPFOLD)
         CALL READI(PFOLDINT)
         CALL READF(OMEGA)
         CALL READF(PFOLDCONV)
C
C  The value of the Planck constant in units of prevailing energy * seconds. 
C  This is needed for regrouped ts free energies. 
C  For example, for CHARMM we enter the temperature in kcal/mol (i.e. k_B * T),
C  so h needs to be in kcal/mol * seconds i.e. 9.546 * 10^(-14).
C
      ELSE IF (WORD.EQ.'PLANCK') THEN
         CALL READF(PLANCK)
      ELSE IF (WORD.EQ.'PSCALE') THEN
         CALL READF(PSCALE)
!
! PHI4MOD for mean field phi^4 model
!
      ELSE IF (WORD.EQ.'PHI4MOD') THEN
         PHI4MODT=.TRUE.
         IF (NITEMS.GT.1) CALL READF(JPARAM)
         PRINT '(A,G20.10)',' keywords> PHI4 mean field model with J=',JPARAM

!
! Allows to use more than one pruning cycle, default is 5
!
      ELSE IF (WORD.EQ.'PRUNECYCLE') THEN
         PRUNECYCLET=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NPRUNE)
         PRINT '(A,I6)',' keywords> Pruning cycles set to: ', NPRUNE


!
! Print some summary info about the kinetic transition network to
! standard out
!
      ELSE IF (WORD.EQ.'PRINTSUMMARY') THEN
          PRINTSUMMARYT=.TRUE.
C
C  Try connecting minima that have the largest values of equilibrium occupation probability times
C  waiting time to minima on the fastest path.
C  Don't bother - UNTRAP with free energy sorting should do the job.
C
!     ELSE IF (WORD.EQ.'PTAU') THEN
!        PTAUT=.TRUE.
!        DIJKSTRAT=.TRUE.
      ELSE IF (WORD.EQ.'PULL') THEN
         PULLT=.TRUE.
         CALL READI(PATOM1)
         CALL READI(PATOM2)
         CALL READF(PFORCE)
         IF (PFORCE.EQ.0.0D0) THEN
            WRITE(*,'(A,I6,A,I6,A,G20.10)') 'keyword> WARNING *** Pulling force is zero, turning off pulling directive'
            PULLT=.FALSE.
         ELSE
            WRITE(*,'(A,I6,A,I6,A,G20.10)') 'keyword> Pulling atoms ',PATOM1,' and ',PATOM2,' force=',PFORCE
         ENDIF
         PRINT '(A)','keywords> Constant pulling force with 4 zero eigenvalues'
      ELSE IF (WORD.EQ.'MULTISITEPY') THEN
              OPEN(UNIT=299,FILE="pysites.xyz",STATUS="old")
              READ(299,*) NELLIPSOIDS
              READ(299,*)
         MULTISITEPYT=.TRUE.
         IF(NITEMS.GT.1) THEN
           NRBSITES=2*NELLIPSOIDS
           CALL READF(RBSITEPOS)
           WRITE(*,'(A,F8.2,A)') ' keywords> symmetric sites will be positioned along the z axis at a distance of ', RBSITEPOS,
     &          'from the body centre'
           ALLOCATE(RBSITE(NRBSITES,3))
           ALLOCATE(RBCENTRE(NELLIPSOIDS,3))
           ALLOCATE(SITEMASS(NRBSITES))
           SITEMASS(:) = 1.D0
           DO J1=1,NELLIPSOIDS
              READ(299,*) LABEL, RBCENTRE(J1,1), RBCENTRE(J1,2), RBCENTRE(J1,3)
           END DO ! loop over all sites
           RBSITE(:,:)=0.0D0
           DO J1=1,NELLIPSOIDS
              RBSITE(J1,3)=RBCENTRE(J1,3)+RBSITEPOS
              RBSITE(J1+NELLIPSOIDS,3)=RBCENTRE(J1,3)-RBSITEPOS
           END DO
         ELSE
           NRBSITES=NELLIPSOIDS
           ALLOCATE(RBSITE(NRBSITES,3))
           ALLOCATE(SITEMASS(NRBSITES))
           SITEMASS(:) = 1.D0
           DO J1=1,NRBSITES
              READ(299,*) LABEL, RBSITE(J1,1), RBSITE(J1,2), RBSITE(J1,3)
           END DO ! loop over all sites
         END IF
      ELSE IF (WORD.EQ.'PY') THEN
         NRBSITES = 1
         ALLOCATE(RBSITE(NRBSITES,3))
         ALLOCATE(SITEMASS(NRBSITES))
         SITEMASS(:)=1.0D0
         RBSITE(1,:) = 0.D0
      ELSE IF (WORD.EQ.'PYGRAVITY') THEN
         EFIELDT=.TRUE.
 !
! Calculate a landscape Shannon entropy from harmonic equilibrium probabilities
!
      ELSE IF (WORD.EQ.'SHANNON') THEN
         SHANNONT=.TRUE.
         SHANNONTMIN=1.0D0
         SHANNONTMAX=2.0D0
         SHANNONTINC=0.1D0
         CALL READF(SHANNONTMIN)
         CALL READF(SHANNONTMAX)
         CALL READF(SHANNONTINC)
         IF (NITEMS.GT.3) CALL READF(EINC)
         IF (NITEMS.GT.4) CALL READI(NPEQ)
!
! Calculate a landscape Shannon entropy from harmonic equilibrium probabilities
! Also calculate other frustration measures that include rates
!
      ELSE IF (WORD.EQ.'SHANNONR') THEN
         SHANNONT=.TRUE.
         SHANNONRT=.TRUE.
         SHANNONTMIN=1.0D0
         SHANNONTMAX=2.0D0
         SHANNONTINC=0.1D0
         CALL READF(SHANNONTMIN)
         CALL READF(SHANNONTMAX)
         CALL READF(SHANNONTINC)
         IF (NITEMS.GT.3) CALL READF(EINC)
         IF (NITEMS.GT.4) CALL READI(NPEQ)
!
! Calculate a landscape Shannon entropy from harmonic equilibrium probabilities
! Also calculate an estimated Zscore  
!
      ELSE IF (WORD.EQ.'SHANNONZ') THEN
         SHANNONT=.TRUE.
         SHANNONZT=.TRUE.
         SHANNONTMIN=1.0D0
         SHANNONTMAX=2.0D0
         SHANNONTINC=0.1D0
         CALL READF(SHANNONTMIN)
         CALL READF(SHANNONTMAX)
         CALL READF(SHANNONTINC)
         IF (NITEMS.GT.3) CALL READF(EINC)
         IF (NITEMS.GT.4) CALL READI(NPEQ)
!
! Define the sleep time between submitting jobs in cycle2.
! Non-zero values required for distributed memory machines and
! small systems.
!
      ELSE IF (WORD.EQ.'SILANE') THEN

         SILANET = .TRUE.
         NRBSITES = 5
         ALLOCATE(RBSITE(NRBSITES,3))
         ALLOCATE(SITEMASS(NRBSITES))

         CALL DEFSILANE()

      ELSE IF (WORD.EQ.'SLEEPTIME') THEN
         CALL READF(SLEEPTIME1)
         SLEEPTIME2=SLEEPTIME1
         IF (NITEMS.GT.1) CALL READF(SLEEPTIME2)
      ELSE IF (WORD.EQ.'SSH') THEN
         SSHT=.TRUE.
      ELSE IF (WORD.EQ.'ST') THEN
         STOCKAAT = .TRUE.
         NRBSITES = 1
         ALLOCATE(RBSITE(NRBSITES,3))
         ALLOCATE(SITEMASS(NRBSITES))
         RBSITE(1,:) = 0.D0
         SITEMASS(1) = 1.D0
      ELSE IF (WORD.EQ.'RANDOMMETRIC') THEN
         RANDOMMETRICT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NRANDOMMETRIC)
C
C  Number of random initial orientations to try in minpermdist.
C
      ELSE IF (WORD.EQ.'RANROT') THEN
         CALL READI(NRANROT)
C
C  Calculate rate constants at the end of each cycle - read in the
C  temperatures required.
C
      ELSE IF (WORD.EQ.'RATESCYCLE') THEN
         RATESCYCLET=.TRUE.
         ALLOCATE(RATESCYCLETEMPS(NITEMS-1))
         NRATESCYCLETEMPS=NITEMS-1
         DO J1=1,NRATESCYCLETEMPS
            CALL READF(RATESCYCLETEMPS(J1))
         ENDDO
         PRINT '(A)','keyword> Rate constants will be calculated each cycle for temperatures:'
         PRINT '(6F15.5)',RATESCYCLETEMPS(1:NRATESCYCLETEMPS)

      ELSE IF (WORD.EQ.'RBAA') THEN
         CALL READI(NATOMS)
         NTSITES  = NATOMS*NRBSITES
         NATOMS=NATOMS*2
         NOPT=3*NATOMS
         IF (DBPTDT) NTSITES = (NATOMS/2-1)*NRBSITES + 4
         IF (NRBSITES == 0) THEN
            PRINT *, 'NRBSITES not yet defined'
            STOP
         ENDIF
         RBAAT = .TRUE.
      ELSE IF (WORD.EQ.'RBSYM') THEN
         RBSYMT=.TRUE.
         INQUIRE(FILE='rbsymops',EXIST=RBSYMTEST)
         IF (RBSYMTEST) THEN
            OPEN(UNIT=1,FILE='rbsymops',STATUS='OLD')
            READ(1,*) NRBGROUP
            ALLOCATE(RBOPS(4,NRBGROUP))
            READ(1,*) ((RBOPS(J1,J2),J1=1,4),J2=1,NRBGROUP)
            PRINT '(A,I6)','keywords> number of symmetry operations for rigid body=',NRBGROUP
            DO J1=1,NRBGROUP
               PRINT '(A,I6)','keywords> rigid-body symmetry operation', J1
               RBOPS(4,J1) = RBOPS(4,J1)*ATAN(1.D0)/45.D0
               PRINT '(3F20.10)',RBOPS(1:4,J1)
            ENDDO
         ELSE
            PRINT '(A)','keywords> ERROR *** missing file rbsymops'
            STOP
         ENDIF

      ELSE IF (WORD.EQ.'READMIN') THEN
         READMINT=.TRUE.
         CALL READA(MINNAME)
C
C  Threshold for redefining the A/B/I sets on the basis of a PE superbasin analysis
C  at energy REGROUPTHRESH.
C
      ELSE IF (WORD.EQ.'REGROUP') THEN
         NOPOINTS=.TRUE.
         REGROUPT=.TRUE.
         CALL READF(REGROUPTHRESH)
C
C  Regroup on the basis of free energies (equivalent to intergroup rates).
C  Subsequent rate calculations use the new groups and free energies.
C
      ELSE IF (WORD.EQ.'REGROUPFREE') THEN
         NOPOINTS=.TRUE.
         REGROUPFREET=.TRUE.
         CALL READF(REGROUPFREETHRESH)
C
C  Regroup on the basis of free energies (equivalent to intergroup rates).
C  Subsequent rate calculations use pe stationary points, but the
C  A and B groups are expanded based on the free energy regrouping.
C
      ELSE IF (WORD.EQ.'REGROUPFREEAB') THEN
         NOPOINTS=.TRUE.
         REGROUPFREEABT=.TRUE.
         CALL READF(REGROUPFREETHRESH)
C
C  Regroup stochastically on the basis of a time scale defined by the
C  RFKMC keyword.
C
      ELSE IF (WORD.EQ.'REGROUPKMC') THEN
         REGROUPKMCT=.TRUE.
C
C  Regroup on the basis of persistence analysis.
C
      ELSE IF (WORD.EQ.'REGROUPPERSIST') THEN
         NOPOINTS=.TRUE.
         REGROUPPERSISTT=.TRUE.
         CALL READF(REGROUPFREETHRESH)
C
C  Subtract the energy of the global PE minimum from minima and ts in setup.f.
C
      ELSE IF (WORD.EQ.'RELATIVEE') THEN
         RELATIVEET=.TRUE.
C
C  Regroup on the basis of free energies (equivalent to intergroup rates).
C  Tabulate values for the groups over a range of temperature for a given
C  regroup threshold, then stop.
C
      ELSE IF (WORD.EQ.'RFMULTI') THEN
         NOPOINTS=.TRUE.
         RFMULTIT=.TRUE.
         IF (NITEMS.LT.6) THEN
            PRINT '(A)','keyword> ERROR *** not enough parameters on RFMULTI line'
            STOP
         ENDIF
         CALL READF(TIMESCALE)
         CALL READF(RFMULTITLOW)
         CALL READF(RFMULTITINC)
         CALL READI(RFMULTIN)
         CALL READF(PFSHIFT)
         PRINT '(A,G20.10)','keywords> Running free energy regrouping for timescale=',TIMESCALE
         PRINT '(2(A,G20.10),A,I6,A,F15.5)','          from T=',RFMULTITLOW,' in steps of delta T=',RFMULTITINC, 
     &                      ' steps=',RFMULTIN,' shift=',PFSHIFT
C
C  Regroup on the basis of free energies (equivalent to intergroup rates).
C  KMC-type simulation starting from a given temperature and the lowest
C  free energy minimum at that value.
C
      ELSE IF (WORD.EQ.'RFKMC') THEN
         NOPOINTS=.TRUE.
         RFKMCT=.TRUE.
         IF (NITEMS.LT.8) THEN
            PRINT '(A)','keyword> ERROR *** not enough parameters on RFKMC line'
            STOP
         ENDIF
         CALL READF(TIMESCALE)
         CALL READF(RFKMCTRATE)
         CALL READF(RFKMCTSTART)
         CALL READF(RFKMCTINC)
         CALL READI(RFKMCN)
         CALL READI(RFKMCSTEPS)
         CALL READF(PFSHIFT)
         PRINT '(A,G20.10,A,I6,A)','keywords> Running free energy regrouping KMC-type scheme for dT/dt=',
     &          RFKMCTRATE, ' averaging over ',RFKMCSTEPS,' instances' 
         PRINT '(2(A,G20.10),A,I6,A,F15.5)',' starting from T=',RFKMCTSTART,' in steps of delta T=',RFKMCTINC, 
     &                      ' for ',RFKMCN,' temperature steps and partition function shift=',PFSHIFT
C
C  Threshold for regrouping and calculating free energies
C  on the basis of a superbasin analysis at PE REGROUPPETHRESH
C
      ELSE IF (WORD.EQ.'REGROUPPE') THEN
         NOPOINTS=.TRUE.
         REGROUPPET=.TRUE.
         CALL READF(REGROUPPETHRESH)
C
C  Threshold for regrouping and calculating free energies
C  on the basis of rates.
C
      ELSE IF (WORD.EQ.'REGROUPRATE') THEN
         NOPOINTS=.TRUE.
         REGROUPRATET=.TRUE.
         CALL READF(REGROUPRATETHRESH)
C
C  Remove stationary points specified in file min.remove and ts.remove
C
      ELSE IF (WORD.EQ.'REMOVESP') THEN
         REMOVESP=.TRUE.
C
C  Remove stationary points specified in file min.remove and ts.remove
C
      ELSE IF (WORD.EQ.'REMOVEUNCONNECTED') THEN
         REMOVEUNCONNECTEDT=.TRUE.
          IF (NITEMS.GT.1) CALL READA(UNCONNECTEDS)
C
C  Retain only stationary points specified in file min.retain and ts.retain
C
      ELSE IF (WORD.EQ.'RETAINSP') THEN
         RETAINSP=.TRUE.
C
C  Reweighting for reactant minima to allow stochastic sampling of reactant
C  in a GT calculation.
C
      ELSE IF (WORD.EQ.'REWEIGHT') THEN
         REWEIGHTT=.TRUE.
         CALL READI(NRWBINS)      ! number of bins in probability distribution
         CALL READI(NRWREACTANT)  ! number of reactant minima to choose
         CALL READA(RWENERGYFILE) ! name of file containing quench data
         IF (ALLOCATED(RWPROB)) DEALLOCATE(RWPROB)
         ALLOCATE(RWPROB(NRWBINS))
         OPEN(UNIT=1,FILE=TRIM(ADJUSTL(RWENERGYFILE)),STATUS='OLD')
         RWEMAX=-1.0D100
         RWEMIN=1.0D100
         NDUMMY=0
         DO
            READ(1,*,END=222) DUMMY
            IF (DUMMY.GT.RWEMAX) RWEMAX=DUMMY
            IF (DUMMY.LT.RWEMIN) RWEMIN=DUMMY
            NDUMMY=NDUMMY+1
         ENDDO
222      CONTINUE
         PRINT '(A,I8,2A)','keyword> ',NDUMMY,' energies read from file ',TRIM(ADJUSTL(RWENERGYFILE))
         RWBINWIDTH=(RWEMAX-RWEMIN)/NRWBINS
         REWIND(1)
         RWPROB(1:NRWBINS)=0.0D0
         DO J1=1,NDUMMY
            READ(1,*) DUMMY
            NDUMMY2=INT((DUMMY-RWEMIN-1.0D-10)/RWBINWIDTH) + 1
!           PRINT '(A,3G20.10,I6)','RWEMIN,RWEMAX,DUMMY,NDUMMY2=',RWEMIN,RWEMAX,DUMMY,NDUMMY2
            RWPROB(NDUMMY2)=RWPROB(NDUMMY2)+1
         ENDDO
         CLOSE(1)
         RWPROB(1:NRWBINS)=RWPROB(1:NRWBINS)/NDUMMY
         IF (DEBUG) THEN
            PRINT '(A)','keyword> bin energy ranges and probabilities:'
            DUMMY=0.0D0
            DO J1=1,NRWBINS
               PRINT '(I6,2F20.10,G20.10)',J1,RWEMIN+RWBINWIDTH*(J1-1),RWEMIN+RWBINWIDTH*J1,RWPROB(J1)
               DUMMY=DUMMY+RWPROB(J1)
            ENDDO
            PRINT '(A,G20.10)','keyword> sum of probabilities=',DUMMY
         ENDIF
C
C  Number of rigid bodies - essential unless NATOMS keyword is present.
C  All we should have to do is then set NATOMS equal to twice the number
C  of rigid bodies to get all the dimensions right.
C
      ELSE IF (WORD.EQ.'RIGIDBODIES') THEN
         CALL READI(NATOMS)
         NATOMS=NATOMS*2
         NOPT=3*NATOMS

! hk286
      ELSE IF (WORD.EQ.'RIGIDINIT') THEN
         IF (NITEMS.EQ.1) THEN
            PRINT '(A)'," keyword> ERROR: you must specify the number of degrees of freedom as an argument to RIGIDINIT"    
            STOP
         ELSE
            RIGIDINIT = .TRUE.         
            PRINT '(A)'," keyword> the generalised rigid body framework is being used"
            CALL READI(DEGFREEDOMS) 
            PRINT '(A,I15)'," keyword> the number of degrees of freedom in the rigid body system set to ",DEGFREEDOMS
         ENDIF
C
C  Random number seed.
C
      ELSE IF (WORD.EQ.'SEED') THEN
         CALL READI(ISEED)
C
C  Try connecting closest minima that are further apart than a given
C  separation in terms of steps in the best path.
C
      ELSE IF (WORD.EQ.'SHORTCUT') THEN
         SHORTCUTT=.TRUE.
         DIJKSTRAT=.TRUE.
         CALL READI(MINSEP)
         IF (NITEMS.GT.2) CALL READA(UNSTRING)
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'BARRIER') BARRIERSHORT=.TRUE.
         IF (TRIM(ADJUSTL(UNSTRING)).EQ.'RATE') RATESHORT=.TRUE.

      ELSE IF (WORD.EQ.'SKIPPAIRS') THEN
         SKIPPAIRST = .TRUE.

      ELSE IF (WORD.EQ.'SIS')THEN
!
! SIS epidemiological model
!
         SIST=.TRUE.
         CALL READI(SMAX)
         CALL READI(IMAX)
         CALL READI(POPSS)
         CALL READF(SISMU)
         CALL READF(SISKAPPA)
         CALL READF(SISBETA)
         ZSYM='C' ! just setting this randomly in case it shouldn't be undefined...
         NATOMS=1 ! Just setting this randomly in case it shouldn't be undefined...
         PRINT '(A,3I6,1X,3G17.10)','keywords> SIS parameters ',SMAX,IMAX,POPSS,SISMU,SISKAPPA,SISBETA

!      ELSE IF (WORD.EQ.'ST') THEN

!         STOCKAAT = .TRUE.
!         NRBSITES = 1
!         ALLOCATE(RBSITE(NRBSITES,3))
!         RBSITE(1,:) = 0.D0
C
C  Node specification in nodes.info file in slurm format.
C
      ELSE IF (WORD.EQ.'SLURM') THEN
         SLURMT=.TRUE.
         CALL GETNODES(NCPU)
      ELSE IF (WORD.EQ.'CUDA') THEN
         CUDAT=.TRUE.
C
C  Make the initial min.A and min.B files using information in path.info.<PATHNAME> and
C  output.<PATHNAME>. This sets up a single A and a single B minimum, which are specified
C  by STARTMINA and STARTMINB. For path.info files in the DUMPALLPATHS format, these
C  will usually not be the first and last entries!
C
      ELSE IF (WORD.EQ.'STARTFROMPATH') THEN
         STARTFROMPATH=.TRUE.
         CALL READA(PATHNAME)
         CALL READI(STARTMINA)
         CALL READI(STARTMINB)
C
C  OPTIM system symbol, e.g. AX, LS, etc.
C
      ELSE IF (WORD.EQ.'SYSTEM') THEN
         CALL READA(ZSYM)
C
C  Number and mass of a tagged atom. 
C
      ELSE IF (WORD.EQ.'TAG') THEN
         IF (NATOMS.LT.1) THEN
            PRINT '(A)','keywords> ERROR - number of atoms must be set before TAG keyword in pathdata'
            STOP
         ENDIF
         IF (.NOT.ALLOCATED(TAGFAC)) THEN
            TAGT=.TRUE.
            ALLOCATE(TAGFAC(NATOMS),TAGNUM(NATOMS))
            TAGFAC(1:NATOMS)=1.0D0
            TAGNUM(1:NATOMS)=0
         ENDIF
         NTAG=NTAG+1
         CALL READI(TAGNUM(NTAG))
         CALL READF(TAGFAC(TAGNUM(NTAG)))
C
C  Target rates. PATHSAMPLE will stop if we get close enough to these rate constants.
C  Used for testing purposes with LJ38 demo.
C
      ELSE IF (WORD.EQ.'TARGETRATES') THEN
         RATETARGETT=.TRUE.
         RATETARGETFRAC=0.9D0
         IF (NITEMS.GT.1) CALL READF(RATETARGETAB)
         IF (NITEMS.GT.2) CALL READF(RATETARGETBA)
         IF (NITEMS.GT.3) CALL READF(RATETARGETFRAC)
         PRINT '(A,2G20.10)','keywords> Target rate constants for AB and BA are: ',RATETARGETAB,RATETARGETBA
         PRINT '(A,G20.10)','keywords> Fraction of rates required for termination=',RATETARGETFRAC
C
C  Canonical temperature.
C
      ELSE IF (WORD.EQ.'TEMPERATURE') THEN
         CALL READF(TEMPERATURE)
         ENSEMBLE='T'
C
C  NTFOLD is the maximum number of iterations for the Tfold subroutine.
C  Make it real to prevent overflow.
C
      ELSE IF (WORD.EQ.'TFOLD') THEN
         TFOLDT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NCONNMIN)
         IF (NITEMS.GT.2) CALL READF(NTFOLD)
         IF (NITEMS.GT.3) CALL READF(TFOLDTHRESH)
         IF (NITEMS.GT.4) CALL READF(TOMEGA)
C
C  Keyword TRAP sets an ion trap potential. 
C  Needed to set the number of non-zero normal mode frequencies.
C
      ELSE IF (WORD.EQ.'TRAP') THEN
         TRAPT=.TRUE.
C
C  TSTHRESH discards transition states above the specfied threshold. May be useful
C  for producing a better initial path and excluding CHARMM transition states with
C  silly energies from the database.
C
      ELSE IF (WORD.EQ.'TSTHRESH') THEN
         CALL READF(TSTHRESH)
C
C  Two-dimensional flatland.
C
      ELSE IF (WORD.EQ.'TWOD') THEN
         TWOD=.TRUE.
C
C jmc Set unres potential. NDIHE is the number of dihedral angles in the system, which
C is used in perturb and tssearch to perturb the system randomly
C
      ELSE IF (WORD.EQ.'UNRES') THEN
         UNRST=.TRUE.
         CALL READI(NDIHE)
         ZSYM='C'
C
C Choose connection pairs based on pe barriers to minima in the product set
C
      ELSE IF (WORD.EQ.'UNTRAP') THEN
         UNTRAPT=.TRUE.
         CALL READF(EINC)
         CALL READF(EUNTRAPTHRESH)
         IF (NITEMS.GT.3) THEN
           CALL READF(EDELTAMIN)
           IF (EDELTAMIN.EQ.0.0D0) THEN 
            PRINT '(A)','keywords> ERROR - EDELTAMIN cannot be zero'
            STOP
           ENDIF 
         ENDIF 
         IF (NITEMS.GT.4) THEN
           CALL READF(ELOWBAR)
         ENDIF
         IF (NITEMS.GT.5) THEN
           CALL READF(EHIGHBAR)
         ENDIF 
C
C Use a particular metric for untrap pairing rather than distances.
C
      ELSE IF (WORD.EQ.'UNTRAPMETRIC') THEN
         UNTRAPMETRICT=.TRUE.
         CALL READI(METRICUPAIR)
         IF (NITEMS.GT.2) THEN
         CALL READI(METMATMAX)
         ENDIF 
         IF (NITEMS.GT.3) THEN
         CALL READF(BAILDIST)
         ENDIF 
C
C Choose connection pairs based on the sequence read from file USEPAIRSFILE
C
      ELSE IF (WORD.EQ.'USEPAIRS') THEN
         USEPAIRST=.TRUE.
         CALL READA(USEPAIRSFILE)

C dg413: Introduced USERPOT here! 
      ELSE IF (WORD.EQ.'USERPOT') THEN
         USERPOTT=.TRUE.

C need a general way to communicate the number of zero frequencies - used for PY
C
      ELSE IF (WORD.EQ.'ZEROS') THEN
         CALL READI(NZEROS)
         PRINT '(A,I6)','keywords> number of zero frequencies will be set to ',NZEROS

      ELSE

         CALL REPORT('Unrecognized command '//WORD,.TRUE.)
         STOP
      ENDIF

      CALL FLUSH(6)
      GOTO 190

      RETURN
      END
