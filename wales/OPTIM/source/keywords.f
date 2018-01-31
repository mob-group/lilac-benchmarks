! OPTIM: A prograa for optimizing geometries and calculating reaction pathways
! Copyright (C) 1999-2006 David J. Wales
! This file is part of OPTIM.
! 
! OPTIM is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! at your option) any later version.
! 
! OPTIM is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
! 
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
! 
! All the keywords possible for the odata file are contained here in
! alphabetical order. Initialisation statements procede the big IF block.
! 
      SUBROUTINE KEYWORDS(Q)
         USE COMMONS
         USE KEY
         USE MODMEC
         USE MODTWOEND
         USE MODAMBER9, ONLY : COORDS1,IH,M04,PRMTOP,SALTCON,IGB,CUT,RGBMAX,
     &   NOCISTRANSRNA,NOCISTRANSDNA,ATMASS1,CHECKCISTRANSALWAYS,CHECKCISTRANSALWAYSRNA,CHECKCISTRANSALWAYSDNA,
     &   AMBERICT,AMBSTEPT,AMBIT,AMBPERTT, PERTHRESH, AMBOLDPERTT,AMBICDNEBT,
     &   AMBPDB_UNIT, AMBRST_UNIT, MDCRD_UNIT, MDINFO_UNIT,
     &   KTWN, KTWNT, DUMPMODEN, UACHIRAL, NOPERMPROCHIRAL, FROZENAMBER
         USE MODNEB
         USE MODMXATMS   ! needed for CHARMM
         USE MODCHARMM
         USE MBPOLMOD, ONLY: MBPOLINIT
         USE MODUNRES
         USE KEYNEB, NNNIMAGE=>NIMAGE
         USE KEYCONNECT
         USE MODGUESS
         USE PORFUNCS
         USE PYMODULE, only : BOXLX,BOXLY,BOXLZ
         USE MSEVB_COMMON, ONLY: shellsToCount, maxHbondLength, minHbondAngle, OOclash_sq, printCoefficients
         USE WC
         USE BINARYIO
         USE GSDATA, ONLY : CUBSPLT, GSUPDATE,
     $   GSGROWTOL, GSMXSTP,GSCONV, REPARAMTOL, EVOLVESTRINGT,
     $   GSITERDENSITY, FIXATMS, MAXLENPERIM,
     $   HESSGRAD, GSMAXTOTITD, MAXGROWSTEPS, GSDGUESS,
     $   NOLBFGS, PREROTATE, GSTANTYPE=>TANTYPE
         USE CUBSPLSTRING, ONLY : ARCTOL, DQAGKEY
         USE INTCOMMONS, ONLY : NATINT, INTNEWT, BBCART, INTINTERPT, INTERPSIMPLE,
     $   INTMINPERMT, INTERPCHOICE, NINTIM, CARTRESSTART, INTPARFILE,
     $   MINBACKTCUT, INTERPBACKTCUT, PRINTCOORDS, DESMINT, NURINGS, URINGS,
     $   NUBONDS, UBONDS, USEPARFILE, CHICDNEB, OLDINTMINPERMT,
     $   GLYCART, INTDISTANCET, RIGIDBONDS, GMAXINT, BONDSFROMFILE
         ! MCP
         USE AMHGLOBALS
         USE BGUPMOD
         ! hk286
         USE GENRIGID
         USE MULTIPOT
         ! includes toggle for switching frames - required for RBAA
         USE MODHESS, ONLY : RBAANORMALMODET
         ! jdf43> for MMEINITWRAPPER
         USE ISO_C_BINDING, ONLY: C_NULL_CHAR
! AMBER12
         USE AMBER12_INTERFACE_MOD, ONLY: AMBER12_ATOM, AMBER12_RESIDUE, AMBER12_ATOMS,
     &                                    AMBER12_RESIDUES, AMBER12_GET_COORDS, AMBER12_ATOMSTORES
         USE CHIRALITY
         USE OPEP_INTERFACE_MOD, ONLY : OPEP_INIT
         USE ORBITALS_MOD, ONLY: ORBITALS_INIT
         use sandbox_module, only : num_atoms

         IMPLICIT NONE

         DOUBLE PRECISION ::  Q(3*NATOMS)

         INTEGER NDUM, LUNIT, FUNIT, GETUNIT, MLPDATSTART, MLQDATSTART, NTYEPA
        
         INTEGER ITEM, NITEMS, LOC, LINE, NCR, NERROR, IR, LAST, NTYPEA, J1, J2, J3, J, I
         COMMON /BUFINF/ ITEM, NITEMS, LOC(132), LINE, SKIPBL, CLEAR, NCR,
     &   NERROR, IR, ECHO, LAST, CAT
         ! DOUBLE PRECISION ::AAA,AAB,ABB,PAA,PAB,PBB,QAA,QAB,QBB,ZAA,ZAB
         ! DOUBLE PRECISION :: ZBB,R0AA,R0AB,R0BB
         DOUBLE PRECISION :: XX, EPSAB, EPSBB, SIGAB, SIGBB, RANDOM, DPRAND
         LOGICAL END, SKIPBL, CLEAR, ECHO, CAT, CISTRANS, RBSYMTEST, YESNO
         CHARACTER WORD*25, WW*20, PBC*3
         CHARACTER(LEN=80) :: FILENAME,DUMMYS
         CHARACTER WORD2*25
         ! COMMON /BIN/ NTYPEA,AAA,AAB,ABB,PAA,PAB,PBB,QAA,QAB,QBB,ZAA,ZAB
         ! COMMON /BIN/ ZBB,R0AA,R0AB,R0BB
         COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB,NTYPEA

         INTEGER NATOM, DMODE, NATOMSSAVE
         DOUBLE PRECISION CHX(MXATMS), CHY(MXATMS), CHZ(MXATMS), CHMASS(MXATMS)
         DOUBLE PRECISION DPERT
         DOUBLE PRECISION CHPMIN, CHPMAX, CHNMIN, CHNMAX
         DOUBLE PRECISION, ALLOCATABLE :: LCONGEOM(:,:)
         DOUBLE PRECISION, ALLOCATABLE :: MLPMEAN(:), MLQMEAN(:)
         DOUBLE PRECISION, ALLOCATABLE :: MLPDISTHI(:), MLPDISTHO(:)
         INTEGER, ALLOCATABLE :: MLPINDEXI(:), MLPINDEXO(:)
         INTEGER K1, II1

         INTEGER ISEED

         DOUBLE PRECISION UNRX(NATOMS), UNRY(NATOMS), UNRZ(NATOMS) ! UNRES
         DOUBLE PRECISION DUMMY1(NATOMS)

         DOUBLE PRECISION SLENGTH, EPS, DUMMY
         INTEGER NOK, NBAD
         COMMON /BSNEW/ SLENGTH, NOK, NBAD, EPS
         DOUBLE PRECISION GSQSCALE, GSTHRESH
         INTEGER NSPECIAL, NALLOW, NINFO
         COMMON /G2/ GSTHRESH, GSQSCALE, NSPECIAL, NALLOW, NINFO
         LOGICAL CUBIC
         COMMON /CUB/ CUBIC
         LOGICAL PATHT, DRAGT
         INTEGER NPATHFRAME
         COMMON /RUNTYPE/ DRAGT, PATHT, NPATHFRAME
         LOGICAL CONNECTT, DUMPPATH, READPATH, CALCRATES, STOPFIRST
         INTEGER NCONNECT
         DOUBLE PRECISION TEMPERATURE, HRED
         COMMON /CONN/ STOPFIRST, CONNECTT, NCONNECT, DUMPPATH, READPATH, CALCRATES, TEMPERATURE, HRED
         INTEGER NMOVE
         COMMON /HS/ NMOVE
         ! DOUBLE PRECISION REPELTS(3*NATOMS,100), REPELPUSH
         ! INTEGER NREPELTS, REPELFROM
         ! LOGICAL REPELTST, REPEL
         ! COMMON /OTS/ NREPELTS, REPELTST, REPELPUSH, REPEL, REPELFROM
         INTEGER ISTAT, NDUMMY
         DOUBLE PRECISION STOPDISP
         LOGICAL STOPDISPT, PERMFILE, CONFILE
         COMMON /STOPD/ STOPDISP, STOPDISPT
         DOUBLE PRECISION CAPSRHO, CAPSEPS2, CAPSRAD, HEIGHT
         COMMON /CAPS/ CAPSRHO, CAPSEPS2, CAPSRAD, HEIGHT
         CHARACTER(LEN=20) OSTRING, OTEMP
         CHARACTER(LEN=20) :: PINFOSTRING
         CHARACTER(LEN=20) :: MINPOINTSSTRING
         CHARACTER(LEN=5) :: TEMPSTRING
         CHARACTER(LEN=9) UNSTRING
         CHARACTER(LEN=1) DUMMYCH
         CHARACTER(LEN=100) TOPFILE,PARFILE
         DOUBLE PRECISION LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN,
     &   HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN
         INTEGER IGBNAB, NUNIAX ! sf344
         ! LOCAL AMH VARIABLES
         INTEGER NRES_AMH, I_RES, GLY_COUNT
         ! CHARACTER(LEN=5) TARFL
         ! DOUBLE PRECISION X, Y, Z
         INTEGER :: GROUPCENTRE
         DOUBLE PRECISION :: GROUPRADIUS,DISTGROUPX2,DISTGROUPY2,DISTGROUPZ2,DISTGROUPCENTRE
         CHARACTER (LEN=2) :: FREEZEGROUPTYPE
         LOGICAL :: FREEZEGROUPT, TURNOFFCHECKCHIRALITY, MLPDONE, MLPNORM,  MLQDONE, MLQNORM
         LOGICAL :: RES_IN_LIST
         DOUBLE PRECISION LPI
         INTEGER DATA_UNIT
         CHARACTER(LEN=13) :: AAOPTION
         CHARACTER(LEN=20) :: AMBERSTR

         INTEGER :: MAXNSETS

         DOUBLE PRECISION :: DUMMY_FRQCONV  ! sn402: Used to ensure that the FRQCONV keyword always overrides
                                            ! the default value for a potential

         CHARACTER(LEN=10) :: OPEP_DUMMY,OPEP_DUMMY2

         LPI=3.14159265358979323846264338327950288419716939937510D0
         AAA=0
         AAB=0
         ABB=0
         PAA=0
         PAB=0
         PBB=0
         QAA=0
         QAB=0
         QBB=0
         R0AA=0
         R0AB=0
         R0BB=0
         
         !ds656> Initialise NTYPEA here
         NTYPEA=NATOMS
         !
         ! hk286 - initialise to stationary frame
         RBAANORMALMODET = .FALSE.
         ! generalised rigid body
         ATOMRIGIDCOORDT = .TRUE.
         RIGIDINIT = .FALSE.
         AACONVERGENCET = .FALSE.

         ALIGNRBST = .FALSE.
         SLERPT = .FALSE.

         AVOID_COLLISIONS = .FALSE.
         COLL_TOL = 0.0D0

         ! Multiple potential scheme
         MULTIPOTT= .FALSE.

         ! Thomson problem
         GTHOMSONT = .FALSE.
         GTHOMPOT = 1

         DESMAXEJUMP = HUGE(1.0D0)
         DESMAXAVGE = HUGE(1.0D0)
         UNSTRING='UNDEFINED'
         WELCH=.FALSE.
         TOSI=.FALSE.
         TOSIC6=.FALSE.
         SIO2T=.FALSE.
         SIO2C6T=.FALSE.
         TOSIPOL=.FALSE.
         SIO2PT=.FALSE.
         ALLOCATE(TAGFAC(NATOMS),TAGNUM(NATOMS))
         TAGFAC(1:NATOMS)=1.0D0
         TAGNUM(1:NATOMS)=0
         NTAG=0
         TAGT=.FALSE.

         REPEL=.FALSE.

         INR=-1

         INTMINPERMT =.FALSE. !msb50 internal permutation
         NOPERMPROCHIRAL = .FALSE.
         BFGSMINT=.FALSE.
         GMAX=0.001D0
         BFGSTST=.FALSE.
         HYBRIDMINT=.FALSE.
         KNOWVECS=.FALSE.
         REOPT=.FALSE.
         NOHESS=.FALSE.
         NOFRQS=.FALSE.
         NOIT=.FALSE.
         NEVL=100
         NEVS=500
         NINTS=0
         NBFGSMAX1=10
         NBFGSMAX2=100
         CEIG=1.0D-10  ! changed to a small default to make people change it!
         CHECKINDEX=.FALSE.
         CHECKCONT=.FALSE.
         BFGSSTEP=.FALSE.
         EXTRASTEPS=0.0D0

         DCHECK=.TRUE.

         PRESSURE=.FALSE.
         PV=.FALSE.
         PVTS=.FALSE.
         FRACTIONAL=.FALSE.
         PRESS=0.0D0
         PVCONV=1.0D-3
         PVTOL=1.0D60
         PVSTEPS=100
         NBOXTS=1

         VARIABLES=.FALSE.
         NZERO=0
         EVCUT=0.0D0

         GAUSSIAN=.FALSE.
         GAUSSIAN03=.FALSE.
         GAUSSIAN09=.FALSE.
         GAUSSIAN16=.FALSE.
         CADPAC=.FALSE.
         GAMESSUS=.FALSE.
         GAMESSUK=.FALSE.
         CASTEP=.FALSE.
         CASTEPJOB=''
         QCHEM=.FALSE.
         QCHEMES=.FALSE.
         QCHEMESNAO=0
         QCHEMESNMO=0
         QCHEMESNZERO=0
         QCHEMESNELEC=0
         QCHEMJOB=''
         MOLPRO=.FALSE.
         MOLPROJOB=''
         REAXFFT=.FALSE.
         REAXFFJOB=''
         VASP=.FALSE.
         VASPJOB=''
         ONETEP=.FALSE.
         ONETEPJOB=''
         CP2K=.FALSE.
         CP2KJOB=''
         DFTP=.FALSE.
         CPMD=.FALSE.
         CPMDC=.FALSE.
         PARALLEL=.FALSE.
         NPROC='1'
         DFTBT=.FALSE.
         CPMD_COMMAND='/home/trj25/bin/cpmd.x'
         SCORE_QUEUE=.FALSE.

         ! DC430 >

         DBPT     = .FALSE.
         DBPTDT   = .FALSE.
         DMBLPYT  = .FALSE.
         LWOTPT   = .FALSE.
         GBT      = .FALSE.
         GBDT     = .FALSE.
         MSSTOCKT = .FALSE.
         NCAPT    = .FALSE.
         NIMET    = .FALSE.
         PHI4MODT = .FALSE.
         NIHEAM7T = .FALSE.
         NIHLEPST = .FALSE.
         NIH2LEPST= .FALSE.
         NIHPAIRONLYT = .FALSE.
         NTIPT    = .FALSE.
         PAHAT    = .FALSE.
         PAPT     = .FALSE.
         PATCHYDT = .FALSE.
         PTSTSTT  = .FALSE.
         PYGT     = .FALSE.
         RADIFT   = .FALSE.
         RBAAT    = .FALSE.
         RBSYMT   = .FALSE.
         STOCKAAT = .FALSE.
         SILANET  = .FALSE.
         UNIAXT   = .FALSE.

         ! -----------------------

         ISTCRT=10
         IPRNT=0
         IVEC=0
         IVEC2=0

         MXSTP=0.2D0
         MINMAX=0.01D0
         MAXMAX=0.5D0
         MAXBFGS=0.2D0
         MAXXBFGS=0.2D0
         MAXMBFGS=0.2D0
         MAXNEBBFGS=0.2D0
         MAXINTBFGS=0.2D0
         READGUESS=.FALSE.
         PERMGUESS=.FALSE.

         DTEST=.FALSE.

         MASST=.FALSE.

         VALUEST=.TRUE.
         EFSTEPST=.FALSE.
         EFSTEPS=1
         NVALUES=20
         NSTEPS=1
         BFGSSTEPS=1
         DUMPV=.FALSE.
         DUMPMAG=.FALSE.
         ALLSTEPS=.FALSE.
         ALLVECTORS=.FALSE.
         MWVECTORS=.FALSE.
         READV=.FALSE.

         PGRAD=.FALSE.
         NGRADIENTS=1

         VECTORST=.FALSE.
         NVECTORS=1

         SUMMARYT=.TRUE.
         NSUMMARY=20

         ADMT=.FALSE.
         NADM=20

         CONVU=1.0D-5
         CONVR=1.0D-5
         INDEXT=.TRUE.

         SYMCUT=0.001D0
         TOLD=0.0001D0
         TOLE=0.0001D0
         NHCHECK=6

         TRAD=2.0
         RESIZE=1.0D0

         RTEST=.FALSE.
         JZ=0.0D0
         OMEGA=0.0D0

         PUSHOPTT=.FALSE.
         PUSHOPTMAX=100
         PUSHOPTCONV=1.0D-4
         PUSHOFF=0.01D0
         PUSHCUT=1.0D-5

         BINARY=.FALSE.
         NSTEPMIN=0

         HUPDATE=.FALSE.
         NSTHUP=0
         INTHUP=0
         PHIG=0.0D0
         READHESS=.FALSE.

         SHIFTV=1.0D6

         NORESET=.FALSE.

         MUPDATE=4
         XMUPDATE=4
         MMUPDATE=4
         NEBMUPDATE=4
         INTMUPDATE=4
         GSUPDATE = 4
         GCUPDATE=4
         DGUESS=0.1D0
         XDGUESS=0.1D0
         NEBDGUESS=0.001D0
         INTDGUESS=0.001D0
         GSDGUESS = 0.001D0
         AMBERT=.FALSE.
         NABT=.FALSE.
         NOCISTRANSRNA=.FALSE.
         NOCISTRANSDNA=.FALSE.
         CHECKCISTRANSALWAYS=.FALSE.
         CHECKCISTRANSALWAYSRNA=.FALSE.
         CHECKCISTRANSALWAYSDNA=.FALSE.
         UACHIRAL=.FALSE.

         ! davidg: introduced userpott here
         USERPOTT=.FALSE.

         ! FAKEWATER=.FALSE.

         CHRMMT=.FALSE.
         REDUCEDBONDLENGTHT=.FALSE.
         BLFACTOR=1.D0
         ACESOLV=.FALSE.
         ACEUPSTEP=50
         TWISTDIHET=.FALSE.
         PERTDIHET=.FALSE.
         CHPMAX=0.5d0
         CHPMIN=0.25d0
         CHNMAX=1.0d0
         CHNMIN=0.d0
         CHARMMDFTBT=.FALSE.
         CHARMMNOTUPDATE=.FALSE.
         ISEED=0
         TOMEGAC=.FALSE.
         TSIDECHAIN=.FALSE.
         IMINCUT=0.0D0
         GUESSTST=.False.
         CALCDIHE=.False.
         TRYNEB=.FALSE.
         NOCISTRANS=.FALSE.
         CISTRANS=.FALSE.
         CHECKOMEGAT=.FALSE.
         MINOMEGA=150.D0
         CHECKCHIRALT=.FALSE.
         TURNOFFCHECKCHIRALITY=.FALSE.
         NORANDOM=.FALSE.
         RANDOMCUTOFF=0.d0
         ! GUESSTHRESH=1.0D100
         ENDHESS=.FALSE.
         NENDHESS=0
         ENDNUMHESS=.FALSE.
         ENDNUMHESS2=.FALSE.
         ENDNUMHESSDELTA=1.0D-6
         NPERMDIHE=0
         TWISTTYPE=0
         NGUESS=3
         FAILT=.FALSE.
         OSASAT=.FALSE.
         RPRO=1.4D0
         ODIHET=.FALSE.

         ! unres stuff
         UNRST=.FALSE.
         CONSECT=.FALSE.
         STARTRES=0
         ENDRES=0
         NUMSEC=0

         ! remove rotations and translations
         NOTRANSROTT=.FALSE.

         ! 
         ! AMH  stuff
         AMHT=.FALSE.

         FREEZE=.FALSE.
         FREEZERANGE=.FALSE.
         FREEZEGROUPT=.FALSE.
         FREEZEGROUPTYPE='GT'
         FREEZERES=.FALSE.
         NFREEZE=0
         DO J1=1,NATOMS
            FROZENRES(J1)=.FALSE.
            FROZEN(J1)=.FALSE.
         ENDDO
         ALLOCATE(DUMPMODEN(3*NATOMS))
         DO J1=1,3*NATOMS
            DUMPMODEN(J1)=.FALSE.
         ENDDO
         KEEPINDEX=.FALSE.
         BSMIN=.FALSE.
         RKMIN=.FALSE.
         SLENGTH=0.0D0
         FIXAFTER=-1
         HINDEX=1
         NOK=0
         NBAD=0
         EPS=1.0D-3
         CONTAINER=.FALSE.
         FIXD=.FALSE.
         T12FAC=1.1D0
         PRINTPTS=.FALSE.
         GRADSQ=.FALSE.
         GSQSCALE=1.0D0
         NSPECIAL=-1
         NALLOW=100
         NINFO=0
         GSTHRESH=0.0D0
         TWOD=.FALSE.
         DOUBLET=.FALSE.
         TWOENDS=.FALSE.
         FSTART=1.0D0
         FINC=1.0D0
         RMSTWO=0.001D0
         NTWO=100
         NTWOITER=25
         TWOEVAL=0.0D0
         PATHT=.FALSE.
         STOPFIRST=.FALSE.
         CONNECTT=.FALSE.
         CPPNEBT=.FALSE.
         DUMPPATH=.FALSE.
         DUMPBESTPATH=.FALSE.
         DUMPALLPATHS=.FALSE.
         HESSDUMPT=.FALSE.
         HESSREADT=.FALSE.
         INSTANTONOPTT=.FALSE.
         INSTANTONRATET=.FALSE.
         INSTANTONSTARTDUMPT=.FALSE.
         NIMAGEINST=1
         DISTORTINST=0.4D0
         DELTAINST=1.D-2
         READPATH=.FALSE.
         CALCRATES=.FALSE.
         TEMPERATURE=1.0D0
         TEMPERATURE1=1.D0
         KTWN=207.11
         KTWNT=.FALSE.
         HRED=1.0D0
         NCONNECT=100
         NEWNEBT=.FALSE.
         NEBT=.FALSE.
         NEWCONNECTT=.False.
         SQVVGuess=.FALSE.
         SQVVGuessRMSTol=2.0D0
         NIterSQVVGuessMax=300
         DEBUG=.FALSE.
         CHDEBUG=.FALSE.
         EDEBUG=.FALSE.
         NIMAGE=1
         RMSNEB=0.1
         DTHRESH=2.0D0
         NSTEPNEB=1
         NEBMAG=0
         NMOVE=1
         NPATHFRAME=0
         FRAMEEDIFF=0.0D0
         FRAMESDIFF=0.0D0
         CUBIC=.FALSE.
         BULKT=.FALSE.
         BULK_BOXVEC(:) = 1.d0
         BULKBOXT=.FALSE.
         CUTT = .FALSE.
         POTENTIAL_CUTOFF = 1.D0
         SDT=.FALSE.
         SDOXYGEN=0
         SDHYDROGEN=0
         SDCHARGE=0
         BOWMANT=.FALSE.
         TTM3T=.FALSE.
         BOWMANPES=2
         BOWMANDIR='~/svn/OPTIM/source/Bowman/coef-3b/'
         RATIOS=.FALSE.
         QSPCFWT=.FALSE.
         QTIP4PFT=.FALSE.


         ! EFK: growing strings and freezing nodes
         GROWSTRINGT = .FALSE.
         NOLBFGS = .FALSE.
         HESSGRAD = .FALSE.
         ARCTOL = 1.0D-4
         DQAGKEY = 6
         DESMDEBUG = .FALSE.
         GSMAXTOTITD = -1
         MAXGROWSTEPS = 1.0D3
         EVOLVESTRINGT=.FALSE.
         FREEZENODEST=.FALSE.
         FIXATMS = .FALSE.
         PREROTATE = .FALSE.
         CUBSPLT = .FALSE.
         MAXLENPERIM = 100.0D0
         GSTANTYPE = 1
         REPARAMTOL = 0.75
         GSGROWTOL = 0.25
         GSCONV = 1.0D-3
         GSMXSTP = 0.1
         STOCKT=.FALSE.
         STOCKSPIN = .FALSE.
         STOCKZTOL = 1.0D-4
         STOCKMAXSPIN = 20
         GEOMDIFFTOL=1.0D-1
         GDSQ=0.0D-1
         GDSQT=.FALSE.
         EDIFFTOL=1.0D-6
         NSECDIAG=2


         ! MSEVB parameters

         shellsToCount = 3
         maxHbondLength = 2.5d0
         minHbondAngle = 130.0d0
         OOclash_sq = 4.41d0 ! 2.1^2
         printCoefficients = .FALSE.

         NEBRESEEDT=.FALSE.
         NEBRESEEDINT=100
         NEBRESEEDEMAX=1.0D100
         NEBRESEEDBMAX=1.0D100
         NEBRESEEDDEL1=1.0D5
         NEBRESEEDDEL2=1.0D5
         NEBRESEEDPOW1=2
         NEBRESEEDPOW2=10
         ADDREPT=.FALSE.

         INTLJT=.FALSE.
         INTLJSTEPS=1000
         INTLJTOL=1.0D-3
         INTLJDEL=0.1D0
         INTLJEPS=1.0D0

!
! QCI parameters
!
         CONDATT=.FALSE.
         QCIPOTT=.FALSE.
         QCIPOT2T=.FALSE.
         QCIADDREP=0
         DOBACK=.FALSE.
         DOBACKALL=.FALSE.
         QCIRESET=.FALSE.
         QCISTOP=.FALSE.
         QCIRESETINT1=300
         QCIRESETINT2=1000
         QCIADDACIDT=.FALSE.
         QCITRILAT=.FALSE.
         QCIADDREPCUT=1.0D0
         QCIADDREPEPS=1.0D0
         QCINOREPINT=.FALSE.
         QCIINTREPMINSEP=20
         QCIKADJUSTFRQ=0
         QCIKINTMAX=1.0D2
         QCIKINTMIN=1.0D-2
         QCIKADJUSTTOL=10.0D0
         KADJUSTFRAC=1.05D0
         MAXNACTIVE=0

         FREEZETOL=1.0D-3
         FLATTESTT=.FALSE.
         FLATEDIFF=1.0D-6
         QCIPERMCHECK=.FALSE.
         QCIPERMCHECKINT=100
         INTCONSTRAINTT=.FALSE.
         INTCONSTRAINTTOL=0.1D0
         INTCONSTRAINTDEL=10.0D0
         INTCONSTRAINTREP=100.0D0
         INTCONSTRAINREPCUT=1.7D0
         INTFREEZET=.FALSE.
         INTFREEZETOL=1.0D-3
         INTFREEZEMIN=10
         INTCONFRAC=0.9D0
         INTCONSEP=15
         INTREPSEP=0
         INTSTEPS1=300001
         INTCONSTEPS=100
         INTRELSTEPS=200
         MAXCONUSE=4
         MAXCONE=0.01D0
         INTRMSTOL=0.01D0
         INTIMAGE=3
         MAXINTIMAGE=75
         INTNTRIESMAX=2
         INTIMAGEINCR=6
         INTIMAGECHECK=25
         IMSEPMIN=0.0D0
         IMSEPMAX=HUGE(1.0D0)

         CHECKCONINT=.FALSE.
         CONCUTABS=0.15D0
         CONCUTABST=.TRUE.
         CONCUTFRAC=0.1D0
         CONCUTFRACT=.FALSE.
         CHECKREPINTERVAL=10
         CHECKREPCUTOFF=2.0
         DUMPINTXYZ=.FALSE.
         DUMPINTEOS=.FALSE.
         DUMPINTXYZFREQ=100
         DUMPINTEOSFREQ=100
         KINT=0.0D0
         KINTENDS=0.0D0
         QCIAMBERT=.FALSE.
         INTMINT=.FALSE.
         INTSPRINGACTIVET=.TRUE.
         INTMINFAC=1.0D0
         QCIRADSHIFTT=.FALSE.
         QCIRADSHIFT=1.0D0
         QCICYCLEST=.FALSE.
         QCICYCDIST=0.0
         QCICYCN=100
         QCIDNEBT=.FALSE.
         QCIRESTART=.FALSE.
         QCILPERMDIST=.FALSE.
         QCIPDINT=1000
         QCIPERMCUT=0.8D0

         CONPOTT=.FALSE.
         CPCONSTRAINTTOL=0.1D0
         CPCONSTRAINTDEL=1.0D5
         CPCONSTRAINTREP=1.0D0
         CPCONSTRAINREPCUT=20.0D0
         CPCONFRAC=1.0D-4
         CPREPSEP=0
         CPCONSEP=10000
         CHECKOVERLAPT=.FALSE.
         MINOVERLAP=1.0D0
         CHECKNEGATIVET=.FALSE.
         ORBITTOL=1.0D-3
         READMASST=.FALSE.
         KADJUSTFRQ=5
         KADJUSTTOL=10.0D0
         KADJUSTFRAC=1.05D0
         MODEDOWNT=.FALSE.
         NOINVERSION=.FALSE.
         PMPATHT=.FALSE.
         PMPATHINR=6
         AAORIENTT=.FALSE.
         KAA=1.0D0
         SIGMAAA=0.0D0
         MULTIJOBT=.FALSE.
         MULTISTART=''
         MULTIFINISH=''
         MULTI_COUNT = 1
         MULTI_LAST = HUGE(100000)  ! By default, we will carry on adding configurations until the end of the file, or we reach this ridiculously large number
         MULTI_STEP = 1

         DJWRBT=.FALSE.
         NHEXAMERS=0
         ! 
         ! General mixed LJ systems
         ! 
         GLJT=.FALSE.
         NGLJ=1 ! number of atom types
         ! 
         ! ds656> substrate field(s)
         MIEFT=.FALSE.
         MIEF_PBCT=.FALSE.
         MIEF_CUTT=.FALSE.
         MIEF_BOX(1:3) = 1.0D9
         MIEF_RCUT= 1.0D9
         MAXIMFACTOR=10.0D0
         !
         ! UNDOCUMENTED keywords/parameters
         ! 
         TWISTT=.FALSE.
         LJADDT=.FALSE.
         LJADD2T=.FALSE.
         LJADD3T=.FALSE.
         LJADD4T=.FALSE.
         NADDTARGET=1
         PYADDT=.FALSE.
         PYADD2T=.FALSE.
         INVERTPT=.FALSE.
         DNEBEFRAC=0.0D0
         MORPHT=.FALSE.
         GREATCIRCLET=.FALSE.
         MAXTSENERGY=1.0D100
         MAXBARRIER=1.0D100
         MAXMAXBARRIER=1.0D100
         ReoptimiseEndpoints=.False.
         ANGLEAXIS=.FALSE.
         NFAILMAX=2
         NATBT=.FALSE.
         READSP=.FALSE.
         DUMPSP=.FALSE.
         TIMELIMIT=HUGE(TIMELIMIT)
         RIGIDBODY=.FALSE.
         STOPDISPT=.FALSE.
         NCHENCALLS=0

         REPELTST=.FALSE.
         NREPELTS=0
         REPELPUSH=0.1D0

         DRAGT=.FALSE.

         LANCZOST=.FALSE.
         ACCLAN=1.0D-8
         SHIFTLAN=1.0D-2
         CUTLAN=-1.0D0

         GFRACTION=0.0D0
         MFRACTION1=0.0D0
         MFRACTION2=0.0D0
         FTEST=.FALSE.
         GALPHA=6.0D0
         MALPHA1=6.0D0
         MALPHA2=6.0D0

         FIELDT=.FALSE.
         OHT=.FALSE.
         IHT=.FALSE.
         TDT=.FALSE.
         D5HT=.FALSE.
         FOH=0.0D0
         FIH=0.0D0
         FTD=0.0D0
         FD5H=0.0D0
         MAXERISE=1.0D-10
         XMAXERISE=1.0D-3
         NEBMAXERISE=1.0D-3
         INTEPSILON=1.0D-6

         EFIELD=0.0D0
         COLDFUSIONLIMIT=-1.0D6
         BLNT=.FALSE.
         DUMPDATAT=.FALSE.
         DUMPDATA_MACHINET=.FALSE.
         NDOF = 0  ! Given as an argument to DUMPDATA_MACHINET to specify how many coords go in each record.
         LOWESTFRQT=.FALSE.
         REDOPATH=.FALSE.
         REDOFRAC=0.5D0
         REDOK=0.0D0
         REDOKADD=.FALSE.
         REDOPATHNEB=.FALSE.
         REDOBFGSSTEPS=100
         REDOPATHXYZ=.FALSE.
         REDOTS = 2
         REALIGNXYZ=.FALSE.
         PERMDIST=.FALSE.
         MAXNSETS = 3
         NRANROT=0
         ATOMMATCHDIST=.FALSE.
         ATOMMATCHFULL=.FALSE.
         LOCALPERMDIST=.FALSE.
         LOCALPERMNEIGH=4
         LOCALPERMCUT=0.2D0
         LOCALPERMMAXSEP=3
         LOCALPERMCUT2=10.0D0
         LPERMDIST=.FALSE.
         LPDGEOMDIFFTOL=0.3D0
         RBCUTOFF=4.0D0
         NRBTRIES=1
         ALLOCATE(BESTPERM(NATOMS))
         PERMDISTINIT=.FALSE.
         ATOMMATCHINIT=.FALSE.
         NEBK=1.0D0 ! changed DJW 14/5/08
         NEBKINITIAL=1.0D0
         NEBKFINAL=1.0D0
         NEBFACTOR=1.01D0

         BHDEBUG=.FALSE.
         BHDISTTHRESH=1.0D0
         BHMAXENERGY=1.D100
         BHINTERPT=.FALSE.
         BHACCREJ=0.5D0
         BHSTEPSIZE=0.4D0
         BHCONV=0.01D0
         BHSTEPS=1000
         BHTEMP=1.0D0
         BHINTERPUSELOWEST=.FALSE.
         BHCHECKENERGYT=.FALSE.
         BHSTEPSMIN=0
         BHK=1.0D0
         ICINTERPT=.FALSE.
         CHBIT=.FALSE.
         BISECTT=.FALSE.
         BISECTMAXENERGY=1.D100
         BISECTDEBUG=.FALSE.
         BISECTSTEPS=1
         BISECTMINDIST=1.0
         BISECTMAXATTEMPTS=5
         CHRIGIDT=.FALSE.
         PTRANS=0.D0
         TRANSMAX=0.D0
         PROT=0.D0
         ROTMAX=0.D0
         BBRSDMT=.FALSE.
         AMBERICT=.FALSE.
         AMBSTEPT=.FALSE.
         AMBICDNEBT = .FALSE.
         AMBPERTT=.FALSE.
         AMBOLDPERTT=.FALSE.
         AMBIT = .FALSE.

         DIJKSTRALOCAL=1.0D0
         DNEBSWITCH=-1.0D0
         PATHSDSTEPS=-1
         NUSEEV=-1
         ACKLANDID=5
         ACK1=.FALSE.
         ACK2=.FALSE.

         RINGPOLYMERT=.FALSE.
         RPSYSTEM='     '
         RPIMAGES=1
         RPBETA=1.0D0
         GRAD4T=.FALSE.
         RPCYCLICT=.TRUE.
         RPFIXT=.FALSE.
         EYTRAPT=.FALSE.
         MKTRAPT=.FALSE.
         BFGSTSTOL=0.0001D0
         EVPC=1.0D0
         OHCELLT=.FALSE.
         CONDATT=.FALSE.
         PAIRCOLOURT=.FALSE.
         NENDDUP=0
         REVERSEUPHILLT=.FALSE.
         ! sf344
         SANDBOXT=.FALSE.
         EFIELDT=.FALSE.
         MACROIONT=.FALSE.
         NORMALMODET=.FALSE.
         PYGPERIODICT=.FALSE.
         PYT=.FALSE.
         PYBINARYT=.FALSE.
         MULTISITEPYT=.FALSE.
         LJGSITET=.FALSE.
         LJSITE=.FALSE.
         BLJSITE=.FALSE.
         LJSITECOORDST=.FALSE.
         LJSITEATTR=.FALSE.
         PCUTOFF=999.0D0
         PARAMONOVCUTOFF=.FALSE.
         CLOSESTALIGNMENT=.FALSE.
         DF1T=.FALSE.
         PULLT=.FALSE.
         CHEMSHIFT=.FALSE.
         METRICTENSOR=.FALSE.

         FRQCONV = 1.0D0
         DUMMY_FRQCONV = 0.0D0

         QUIPARGSTRT=.FALSE.
         QUIPPARAMST=.FALSE.
         QUIPZ=14 ! default is silicon
         QUIPT=.FALSE.

         EIGENONLY=.FALSE.
         COVER=0.990d0
         OVERCONV=.FALSE.
         TRUSTMODET=.FALSE.
         TMRATIO=0.71D0

         JPARAM=0.001D0
         RPHT=.FALSE.
         MCPATHT=.FALSE.
         MCPATH2T=.FALSE.
         MCBIAST=.FALSE.
         MCPATHTS=1
         MCPATHSCHECK=0
         MCPATHTOL=1.0D-12
         MCMERGES=1
         MCMERGEQ=1
         PBST=.FALSE.
         MCPATHDOBLOCK=0
         MCPATHNEGLECT=0.0D0
         SSHT=.FALSE.
         NCPU=0
         MCPATHGWS=0.3D0
         MCPATHGWQ=0.00025D0
         COLLAGENOP=.FALSE.

         DUMPFRQST=.FALSE.

         CHIRALENDPOINTS=.TRUE.

         CUDAT=.FALSE.
         CUDAPOT=' '
         CUDATIMET=.FALSE.
         DJWRBT=.FALSE.
         NHEXAMERS=0
         RADPENT=5.0D0
         CAPSIDRHO=3.0D0
         RADHEX=RADPENT*2.0*0.5877852522924731D0  ! 2 * Sin[36] to give the same edge length
         SIGMAPENT=(1.0D0+RADPENT*SQRT((5.0D0+SQRT(5.0D0))/2.0D0))
         SIGMAHEX=(1.0D0+RADHEX*SQRT((5.0D0+SQRT(5.0D0))/2.0D0))
         SIGMAPH=0.5D0*(SIGMAPENT + SIGMAHEX)
         CAPSIDEPS=0.4D0


!
! Stealthy potential
!
         STEALTHYT=.FALSE.
         STEALTV=.FALSE.

!
! Neural network potential
!
         MLP3T=.FALSE.
         MLPB3T=.FALSE.
         MLPB3NEWT=.FALSE.
         MLPVB3T=.FALSE.
         NOREGBIAS=.FALSE.
         MLPNEWREG=.FALSE.
         MLPPROB=.FALSE.
         MLPPROBPOS=1.0D0
         MLPDONE=.FALSE.
         MLPNORM=.FALSE.
         MLPLAMBDA=0.0D0
         MLPNEIGH=1
         MLPDATSTART=1
!
! ML quadratic function
!
         MLQT=.FALSE.
         MLQPROB=.FALSE.
         MLQDONE=.FALSE.
         MLQNORM=.FALSE.
         MLQLAMBDA=0.0D0
         MLQDATSTART=1

         MAXGAPT=.FALSE.

         MALONALDEHYDE=.FALSE.

         MBPOLT= .FALSE.
! OPEP stuff
         OPEPT = .FALSE.
         OPEP_RNAT = .FALSE.

! cs675 ORBITALS stuff
         ORBITALS = .FALSE.
         ORBVAREXPONENT = -1

         CLSTRINGT=.FALSE.
         CLSTRINGTST=.FALSE.
         IF (FILTH2.EQ.0) THEN
            OPEN (5,FILE='odata',STATUS='OLD')
         ELSE
            WRITE(OTEMP,*) FILTH2
            WRITE(OSTRING,'(A)') 'odata.' // TRIM(ADJUSTL(OTEMP))
            OPEN (5,FILE=OSTRING,STATUS='OLD')
         ENDIF
         DATA_UNIT=5

190      CALL INPUT(END)
         IF (.NOT. END) THEN
            CALL READU(WORD)
         ENDIF
         ! 
         ! POINTS - keyword at the end of the list of options after which
         ! the Cartesian coordinates follow. Must be present unless VARIABLES or RINGPOLYMER
         ! is present instead. MACHINE keyword overrides POINTS. If MACHINE is
         ! true coordinates that were read from odata file will be overwritten
         ! with coordinates from a direct access file, in which case section of
         ! odata file after POINTS keyword is used only to read in the labels. (SAT)
         ! 
         IF (END.OR.WORD.EQ.'STOP'.OR.WORD.EQ.'POINTS') THEN
            
            ! sn402: Once we have finished reading keywords, check to see whether we need to overwrite
            ! a default value of FRQCONV
            IF(DUMMY_FRQCONV.NE.0.0D0) THEN
               FRQCONV = DUMMY_FRQCONV
               WRITE(*,*) "keywords> Overwriting default unit conversion factor for frequencies"
               WRITE(*,*) "Conversion factor: ", FRQCONV
            ENDIF
            FRQCONV2 = FRQCONV*FRQCONV
            RETURN
         ENDIF

         IF (WORD.EQ.'    ' .OR.WORD.EQ.'NOTE'.OR.WORD.EQ.'COMMENT'
     &   .OR.WORD.EQ.'\\'.OR.WORD.EQ."!".OR.WORD.EQ."#") THEN
            GOTO 190
! 
! Enforce flatland.
! 
         ELSE IF (WORD .EQ. '2D') THEN
            TWOD=.TRUE.
! 
! AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
! 
! Add an angle-dependent term to angle-axis potential to work around zero eigenvalues
! 
         ELSE IF (WORD.EQ.'AAORIENT') THEN
            AAORIENTT=.TRUE.
            IF (NITEMS.GT.1) CALL READF(KAA)
            IF (NITEMS.GT.2) CALL READF(SIGMAAA)
! 
! bs360: ACE is to be used together with CHARMM and the ACE solvent model,
! it makes sure that the Born radii are regularly updated
! 
         ELSE IF (WORD.EQ.'ACE') THEN
            ACESOLV=.TRUE.
            IF (NITEMS.GT.1) CALL READI(ACEUPSTEP)
! 
! Ackland embedded atom metal potentials.
! 
         ELSE IF (WORD.EQ.'ACKLAND') THEN
            CALL READI(ACKLANDID) ! default is 5 = Fe
! 
! Specification of the two possible Ackland potentials for iron.
! 
         ELSE IF (WORD.EQ.'ACKLAND1') THEN
            ACK1=.TRUE.
         ELSE IF (WORD.EQ.'ACKLAND2') THEN
            ACK2=.TRUE.
! 
! Adjust NEB force constant values between different images on-the-fly in order
! to try and equispace them.
! 
         ELSE IF (WORD.EQ.'ADJUSTK') THEN
            CALL READI(KADJUSTFRQ)
            CALL READF(KADJUSTTOL)
            CALL READF(KADJUSTFRAC)
! 
! ADM [OFF | ON n] prints the atomic distance matrix every n
! if switched on                 cycles       - default n=20
! 
         ELSE IF (WORD .EQ. 'ADM') THEN
            ADMT=.TRUE.
            CALL READI(NADM)


! Keyword for use with RIGIDINIT instructing minpermdist to perform a global alignment
! based only on aligning the specified rigid bodies. This was designed for the plate folding 
! potential, and is unlikely to be very useful for systems which aren't composed mostly of
! large rigid bodies that are strongly bonded together.

        ELSE IF (WORD .EQ. 'ALIGNRBS') THEN
            ALIGNRBST = .TRUE.
            N_TO_ALIGN = NITEMS-1
            ALLOCATE(TO_ALIGN(N_TO_ALIGN))
            DO J1 = 1, N_TO_ALIGN
                CALL READI(TO_ALIGN(J1))
            ENDDO

! Keywork ALPHA enables exponent values to be set for the averaged
! Gaussian and Morse potentials. All defaults = 6.
! 
         ELSE IF (WORD.EQ.'ALPHA') THEN
            CALL READF(GALPHA)
            IF (NITEMS.GT.2) THEN
               CALL READF(MALPHA1)
            ENDIF
            IF (NITEMS.GT.3) THEN
               CALL READF(MALPHA2)
            ENDIF
! 
! SAT: ALLPOINTS turns on printing of coordinates to file points for intermediate steps.
! This is the default.
! 
         ELSE IF (WORD.EQ.'ALLPOINTS') THEN
            PRINTPTS=.TRUE.

! PCW: remove translations and rotations

         ELSE IF (WORD.EQ.'NOTRANSROT') THEN
            NOTRANSROTT=.TRUE.

! davidg: introduced userpott here
         ELSE IF (WORD.EQ.'USERPOT') THEN
            USERPOTT=.TRUE.
            RETURN
! MCP
         ELSE IF (WORD.EQ.'AMH') THEN
            WRITE(6,*)'USING AMH ENERGIES FORCES'
            WRITE(6,*)'CALCULATE ENERGY AND FORCE TABLES  '
            AMHT=.TRUE.
            WRITE(6,*)'AMH FLAG ', AMHT
            WRITE(6,*)'AMH NATOMS ',  NATOMS
            IF (DEBUG) WRITE(6,*)'Entering WALESAMH_INITIAL'

            CALL WALESAMH_INITIAL

            IF (DEBUG)WRITE(6,*)'Leaving WALESAMH_INITIAL'
            IF (DEBUG)WRITE(6,*)'TARFL ',TARFL

            OPEN(30,FILE='proteins/'//TARFL,STATUS='OLD')
            READ(30,*)
            READ(30,*)NRES_AMH
            IF (NRES_AMH.GT.500) THEN
               WRITE(6,*) 'FAILURE NRES_AMH GR THAN 500 CONNECTODATA'
               STOP
            ENDIF
            READ (30,25)(SEQ(I_RES),I_RES=1,NRES_AMH)
25          FORMAT(25(I2,1X))
            CLOSE(30)

            WRITE(6,*)'NRES ',NRES_AMH
            NRES_AMH_TEMP=NRES_AMH

            DO J1=1,NRES_AMH
               Q(9*(J1-1)+1)=X_MCP(9*(J1-1)+1)
               Q(9*(J1-1)+2)=X_MCP(9*(J1-1)+2)
               Q(9*(J1-1)+3)=X_MCP(9*(J1-1)+3)
               Q(9*(J1-1)+4)=X_MCP(9*(J1-1)+4)
               Q(9*(J1-1)+5)=X_MCP(9*(J1-1)+5)
               Q(9*(J1-1)+6)=X_MCP(9*(J1-1)+6)
               Q(9*(J1-1)+7)=X_MCP(9*(J1-1)+7)
               Q(9*(J1-1)+8)=X_MCP(9*(J1-1)+8)
               Q(9*(J1-1)+9)=X_MCP(9*(J1-1)+9)
            ENDDO


! AMBER12 keywords
         ELSE IF (WORD.EQ.'AMBER12') THEN
! Store coordinates in Q
            IF (NITEMS .GT. 1) THEN
               CALL READA(AMBERSTR)
               IF (FILTH2 .NE. 0) THEN
                  WRITE(OTEMP, *) FILTH2
                  WRITE(OSTRING,'(A)') TRIM(ADJUSTL(AMBERSTR))//'.'//TRIM(ADJUSTL(OTEMP))
                  WRITE(*,*) 'ostring=', OSTRING
               ELSE
                  WRITE(OSTRING,'(A)') TRIM(ADJUSTL(AMBERSTR))
               END IF
               WRITE(*,'(A,A)') ' keywords> input coordinates for AMBER12 system will be read from ', OSTRING
               OPEN(UNIT=3827, FILE=TRIM(ADJUSTL(OSTRING)), STATUS='UNKNOWN')
               IF (NITEMS == 3) THEN
                  CALL READA(AMBERSTR)
                  IF (TRIM(ADJUSTL(AMBERSTR)) == 'inpcrd') THEN
                     WRITE(*,'(A)') ' keywords> reading in AMBER restart format'
                     CALL AMBER12_GET_COORDS(NATOMS, Q(1:3*NATOMS))
                  ELSE
                     WRITE(*,'(A)') ' keywords> reading in xyz format'
                     DO I = 1, NATOMS
                        READ(3827, *) Q(3*I-2:3*I)
                     END DO
                  END IF
               ELSE
                  WRITE(*,'(A)') ' keywords> reading in xyz format'
                  DO I = 1, NATOMS
                     READ(3827, *) Q(3*I-2:3*I)
                  END DO
               END IF
               CLOSE(3827)
            ELSE
               CALL AMBER12_GET_COORDS(NATOMS, Q(1:3*NATOMS))
            END IF
! Store atom names in ZSYM
            ZSYM(1:NATOMS) = AMBER12_ATOMS(1:NATOMS) % NAME
! Store atom masses in ATMASS
            IF (.NOT.ALLOCATED(ATMASS)) ALLOCATE(ATMASS(NATOMS))
            ATMASS(1:NATOMS) = AMBER12_ATOMS(1:NATOMS) % MASS
! Turn on chirality and cis/trans checks by default
            IF (.NOT. TURNOFFCHECKCHIRALITY) THEN
               CHECKCHIRALT=.TRUE.
               CALL INIT_CHIRAL(Q(1:3*NATOMS))
            END IF
            IF (.NOT. CISTRANS) THEN
               NOCISTRANS=.TRUE.
               CALL INIT_CIS_TRANS(Q(1:3*NATOMS))
            END IF
! Check perm.allow file is present (without inquire?) and suitable
            IF (PERMDIST .OR. LOCALPERMDIST .OR. LPERMDIST) THEN
               IF (NPERMSIZE(1) == NATOMS) THEN
                  WRITE(*,*)'keyword> Error: PERMDIST/LOCALPERMDIST/LPERMDIST' //
     &                      'is specified for AMBER 12, but there is no perm.allow file'
                  STOP
               END IF
            END IF
! Freezing residues
            IF (FREEZERES) THEN
               DO I = 1, NATOMS
                  RES_IN_LIST = FROZENRES(AMBER12_ATOMS(I) % RES_INDEX)
! UNFREEZE currently not implemented in OPTIM.
!                  IF (UNFREEZE .AND. .NOT. RES_IN_LIST) THEN
!                     FROZEN(I) = .FALSE.
!                     NFREEZE = NFREEZE - 1
!                  ELSE IF (RES_IN_LIST) THEN
                  IF (RES_IN_LIST) THEN
                     FROZEN(I) = .TRUE.
                     NFREEZE = NFREEZE + 1
                  END IF
               END DO
            END IF
! Freeze groups of atoms
            IF (FREEZEGROUPT) THEN
               ! Write a list of FROZEN atoms for use in an (o)data file
               OPEN(UNIT=4431, FILE='frozen.dat', STATUS='UNKNOWN', FORM='FORMATTED')
               DO I = 1, NATOMS
                  ! Work out the distance from GROUPCENTRE to the current atom I
                  DISTGROUPCENTRE = SQRT(SUM((Q(3*GROUPCENTRE-2:3*GROUPCENTRE) - Q(3*I-2:3*I))**2))
                  ! If working in GT mode (default), FREEZE all atoms >GROUPRADIUS from the GROUPCENTRE atom
                  IF((FREEZEGROUPTYPE == "GT") .AND. (DISTGROUPCENTRE > GROUPRADIUS)) THEN
                     NFREEZE = NFREEZE + 1
                     FROZEN(I) = .TRUE.
                     WRITE(4431,'(A,I6)') 'FREEZE ', I
                 ! IF working in LT mode, FREEZE all atoms <GROUPRADIUS from the GROUPCENTRE atom
                  ELSE IF((FREEZEGROUPTYPE == "LT") .AND. (DISTGROUPCENTRE < GROUPRADIUS)) THEN
                     NFREEZE = NFREEZE + 1
                     FROZEN(I) = .TRUE.
                     WRITE(4431,'(A,I6)') 'FREEZE ',I
                  END IF
               END DO
               CLOSE(4431)
            END IF

            ! sn402: added (see comments at keyword FRQCONV)
            IF (DUMMY_FRQCONV .EQ. 0.0D0) THEN
                FRQCONV = 2.045483D13
                WRITE(*,*) "keywords> Set frequency conversion factor to the AMBER default value: ", FRQCONV
                WRITE(*,*) "keywords> This corresponds to frequencies being given in radians/s"
            ELSE
                FRQCONV = DUMMY_FRQCONV
                WRITE(*,*) "keywords> Set frequency conversion factor to the user-specified value: ", FRQCONV
            ENDIF
            FRQCONV2 = FRQCONV*FRQCONV

            ALLOCATE(ATOMSTORES(NATOMS))
            CALL AMBER12_ATOMSTORES(ATOMSTORES,NATOMS)

            WRITE (*,'(A)') 'Warning: AMBER12 keyword must come last in odata'
            RETURN

! sf344> start of AMBER 9 keywords
         ELSE IF (WORD.EQ.'AMBER9') THEN
            AMBERT=.TRUE.
! jmc49> make sure that chirality and cis/trans isomerization checks are on by default
            IF (.NOT.TURNOFFCHECKCHIRALITY) CHECKCHIRALT=.TRUE.
            IF (.NOT.CISTRANS) NOCISTRANS=.TRUE.
! 
! csw34> if FREEZERES specified, populate the FROZEN array with A9RESTOATOM
! 
            IF (FREEZERES) CALL A9RESTOATOM(FROZENRES,FROZEN,NFREEZE)
            IF ((PERMDIST.OR.LOCALPERMDIST.OR.LPERMDIST)) THEN
               IF (NPERMSIZE(1).EQ.NATOMS) THEN
                  WRITE(*,*)'keyword>ERROR-PERMDIST/LOCALPERMDIST/LPERMDIST is specified for AMBER, but there is no perm.allow file'
                  STOP
               ENDIF
            ENDIF
            IF (FREEZEGROUPT) THEN
               ! Write a list of FROZEN atoms for use in an (o)data file
               OPEN(UNIT=4431,FILE='frozen.dat',STATUS='UNKNOWN',FORM='FORMATTED')
               DO J1=1,NATOMS
                  ! 
                  ! Work out the distance from GROUPCENTRE to the current atom J1
                  ! 
                  DISTGROUPX2=(COORDS1(3*GROUPCENTRE-2)-COORDS1(3*J1-2))**2
                  DISTGROUPY2=(COORDS1(3*GROUPCENTRE-1)-COORDS1(3*J1-1))**2
                  DISTGROUPZ2=(COORDS1(3*GROUPCENTRE  )-COORDS1(3*J1  ))**2
                  DISTGROUPCENTRE=SQRT(DISTGROUPX2+DISTGROUPY2+DISTGROUPZ2)
                  ! If working in GT mode (default), FREEZE all atoms >GROUPRADIUS from the GROUPCENTRE atom
                  IF((FREEZEGROUPTYPE=="GT").AND.(DISTGROUPCENTRE.GT.GROUPRADIUS)) THEN
                     NFREEZE=NFREEZE+1
                     FROZEN(J1)=.TRUE.
                     WRITE(4431,'(A,I6)') 'FREEZE ',J1
                     ! IF working in LT mode, FREEZE all atoms <GROUPRADIUS from the GROUPCENTRE atom
                  ELSE IF((FREEZEGROUPTYPE=="LT").AND.(DISTGROUPCENTRE.LT.GROUPRADIUS)) THEN
                     NFREEZE=NFREEZE+1
                     FROZEN(J1)=.TRUE.
                     WRITE(4431,'(A,I6)') 'FREEZE ',J1
                  END IF
               END DO
               CLOSE(4431)
            ENDIF

! 
! csw34> A copy of the FROZEN array called FROZENAMBER is created to be passed through to AMBERINTERFACE
! 
            ALLOCATE(FROZENAMBER(NATOMS))
            FROZENAMBER(:)=FROZEN(:)
            IF(.NOT.ALLOCATED(ATMASS)) ALLOCATE(ATMASS(NATOMS))
            ATMASS(1:NATOMS) = ATMASS1(1:NATOMS)
            DO J1=1,3*NATOMS
               Q(J1) = COORDS1(J1)
            END DO
! save atom names in array zsym
            do J1=1,natoms
            zsym(J1) = ih(m04+J1-1)
            end do
! initialise MME
            CALL MMEINITWRAPPER(TRIM(ADJUSTL(PRMTOP))//C_NULL_CHAR,IGB,SALTCON,RGBMAX,SQRT(CUT))

            ! sn402: added (see comments at keyword FRQCONV)
            IF (DUMMY_FRQCONV .EQ. 0.0D0) THEN
                FRQCONV = 2.045483D13
                WRITE(*,*) "keywords> Set frequency conversion factor to the AMBER default value: ", FRQCONV
                WRITE(*,*) "keywords> This corresponds to frequencies being given in radians/s"
            ELSE
                FRQCONV = DUMMY_FRQCONV
                WRITE(*,*) "keywords> Set frequency conversion factor to the user-specified value: ", FRQCONV
            ENDIF
            FRQCONV2 = FRQCONV*FRQCONV

            RETURN
! initialise unit numbers
            ambpdb_unit=1110
            ambrst_unit=1111
            mdinfo_unit=1112
            mdcrd_unit =1113

         ELSE IF (WORD.EQ.'AMBERIC') THEN
            PRINT*, "amberic"
            AMBERICT = .TRUE.
            IF (NITEMS .GT. 1) THEN
               CALL READA(WORD2)
               IF (WORD2.EQ.'BACKBONE')  THEN
                  PRINT*, "backbone interpolated"
                  AMBIT = .TRUE.
               ELSE
                  PRINT*, "keyword error in amberic"
                  RETURN
               ENDIF
            ENDIF

         ELSE IF (WORD.eq.'AMBERSTEP') THEN
            PRINT*, "amberstept"
            AMBSTEPT = .TRUE.

         ELSE IF (WORD.eq.'AMBPERTOLD') THEN
            PRINT*, "original perturbation scheme"
            AMBOLDPERTT = .TRUE.

         ELSE IF (WORD.eq. 'AMBPERTONLY') THEN
            AMBPERTT = .TRUE.
            CALL READF(PERTHRESH)
            PRINT*, "amber pertonly, perthresh", perthresh

         ELSE IF (WORD.eq. 'AMBICDNEB') THEN
            AMBICDNEBT = .TRUE.

! For bulk  systems, an alternative method for finding the shortest
! distance between structures. Particularly useful for defective crystal
! structures, atoms are overlayed and the number of exactly matching
! atoms maximised (method is non-deterministic to maximise efficiency).
! Use 'ATOMMATCHFULL' for a slow deterministic result.
         ELSE IF ((WORD.EQ.'ATOMMATCHDIST').OR.(WORD.EQ.'ATOMMATCHINIT')) THEN
            ATOMMATCHDIST=.TRUE.
            WRITE(*,'(A)') 'keyword> Atom matching for distance calculation'
            IF (WORD.EQ.'ATOMMATCHINIT') ATOMMATCHINIT=.TRUE.
         ELSE IF (WORD.EQ.'ATOMMATCHFULL') THEN
            ATOMMATCHDIST=.TRUE.
            ATOMMATCHFULL=.TRUE.
            WRITE(*,'(A)') 'keyword> Atom matching for distance calculation'
            WRITE(*,'(A)') ' WARNING - inefficient atom matching for a deterministic result'
! 
         ELSE IF (WORD.eq.'NAB') THEN
            IF (FREEZERES) CALL A9RESTOATOM(FROZENRES,FROZEN,NFREEZE)
! jmc49> make sure that chirality and cis/trans isomerization checks are on by default
            IF (.NOT.TURNOFFCHECKCHIRALITY) CHECKCHIRALT=.TRUE.
            IF (.NOT.CISTRANS) NOCISTRANS=.TRUE.
            NABT=.TRUE.
            DO J1=1,3*NATOMS
               Q(J1) = COORDS1(J1)
            END DO
! save atom names in array zsym
            do J1=1,natoms
            zsym(J1) = ih(m04+J1-1)
            end do
            IF(.NOT.ALLOCATED(ATMASS)) ALLOCATE(ATMASS(NATOMS))
! for the NAB interface, ATMASS is also set up in mme2wrapper, and that setting
! overrides the one from below. However, both originate from the same prmtop file,
! so they should be the same. ATMASS is being assigned here so that it's somewhat consistent
! with the AMBER interface.
            ATMASS(1:NATOMS) = ATMASS1(1:NATOMS)
            WRITE(prmtop,'(A)') 'coords.prmtop'
            igbnab=igb
            if(igb==6) igbnab=0     ! this is also in vacuo, but NAB doesn't understand igb=6!
            CALL MMEINITWRAPPER(trim(adjustl(prmtop))//C_NULL_CHAR,igbnab,saltcon,rgbmax,sqrt(cut))

            ! sn402: added (see comments at keyword FRQCONV)
            IF (DUMMY_FRQCONV .EQ. 0.0D0) THEN
                FRQCONV = 2.045483D13
                WRITE(*,*) "keywords> Set frequency conversion factor to the NAB default value: ", FRQCONV
                WRITE(*,*) "keywords> This corresponds to frequencies being given in radians/s"
            ELSE
                FRQCONV = DUMMY_FRQCONV
                WRITE(*,*) "keywords> Set frequency conversion factor to the user-specified value: ", FRQCONV
            ENDIF
            FRQCONV2 = FRQCONV*FRQCONV

            RETURN

         ELSE IF (WORD.eq.'DF1') THEN
            DF1T=.TRUE.

         ELSE IF (WORD.eq.'DUMPSTRUCTURES') THEN
            DUMPSTRUCTURES=.TRUE.
            WRITE(*,'(A)') ' keywords> Final structures will be dumped in different formats (.rst, .xyz, .pdb)'
! 
! Distinguish between old C of M/Euler and new angle/axis coordinates for
! rigid body TIP potentials
! 
         ELSE IF (WORD.EQ.'ANGLEAXIS') THEN
            ANGLEAXIS=.TRUE.
! 
! Growing string arc tolerance.
! 
         ELSE IF (WORD.EQ.'ARCTOL') THEN
            CALL READF(ARCTOL)
!
! Instructs the NEB interpolation scheme to check for rigid bodies which might be colliding and to 
! attempt to fix this by re-interpolating problematic bodies with the sense of rotation reversed
!
         ELSE IF (WORD.EQ.'AVOIDCOLLISIONS') THEN
            AVOID_COLLISIONS = .TRUE.
            CALL READF(COLL_TOL)

! 
! Specifies the highest symmetry axis to search for in routine {\bf symmetry}; default is six.
! 
         ELSE IF (WORD .EQ. 'AXIS') THEN
            CALL READI(NHCHECK)
! 
! BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
         ELSE IF (WORD.EQ.'BBCART') THEN
            BBCART = .TRUE. ! use cartesians for backbone

         ELSE IF (WORD.EQ.'BBRSDM') THEN
! 
! BBSDM minimiser.
! 
            BBRSDMT = .TRUE.
            CALL READF(BBRGAM)
            CALL READF(BBREPS)
            CALL READF(BBRSIGMA1)
            CALL READF(BBRSIGMA2)
            CALL READI(BBRM)
            CALL READF(BBRALPHA)
            CALL READF(BBRCONV)
            CALL READI(BBRSTEPS)

         ELSE IF (WORD.EQ.'BFGSCONV') THEN
! 
! Turn on LBFGS gradient minimization. GMAX is the convergence
! criterion for the RMS gradient, default 0.001.
! For BFGSTS NEVL and NEVS are the maximum iterations allowed in the searches for
! the largest and smallest eigenvectors, respectively and NBFGSMAX1 is the largest
! number of BFGS steps allowed in the subsequent restricted minimization.
! If the negative eigenvalue appears to have converged then NBFGSMAX2 steps
! are allowed in the tangent space.
! CONVU is used to determine convergence in such runs and BFGSCONV can be used
! to set GMAX, the convergence criteria for the subspace optimization.
! 
! IF REOPT is true the smallest Hessian eigenvector is redetermined after the
! EF step before the tangent space minimisation.
! 
            IF (NITEMS.GT.1) THEN
               CALL READF(GMAX)
            ENDIF
         ELSE IF (WORD.EQ.'BFGSCONVINT') THEN
! 
! LBFGS convergence criterion in internal coordinates.
! needed for rigids bonds.
! standard is negative, which is equivalent to disabled
! since RMS cannot be negative
            IF (NITEMS.GT.1) THEN
               CALL READF(GMAXINT)
            ENDIF

         ELSE IF (WORD.EQ.'BFGSMIN') THEN
! 
! instructs the program to perform an LBFGS minimisation.
! it gmax\/} is the convergence criterion
! for the root-mean-square gradient, default $0.001$.
!
            BFGSMINT=.TRUE.
            IF (NITEMS.GT.1) THEN
               CALL READF(GMAX)
            ENDIF

         ELSE IF (WORD.EQ.'BFGSSTEP') THEN
! 
! If starting from a transition state we just want to take one EF step using
! BFGSTS before calling MYLBFGS (or something else).
! 
            BFGSSTEP=.TRUE.
            BFGSTST=.TRUE.
            IF (NITEMS.GT.1) CALL READF(PUSHOFF)

         ELSE IF (WORD .EQ. 'BFGSSTEPS') THEN
! 
! BFGSSTEPS n sets the number of BFGS optimisation steps to perform
! per call to OPTIM                                    - default n=1
! If BFGSSTEPS is not specified then it is set to the same value as NSTEPS
! 
            CALL READI(BFGSSTEPS)
            IF (NSTEPS.EQ.1) NSTEPS=BFGSSTEPS
         ELSE IF (WORD.EQ.'BFGSTS') THEN
! 
! Hybrid BFGS/eigenvector-following transition state search.
! 
            BFGSTST=.TRUE.
            IF (NITEMS.GT.1) CALL READI(NEVS)
            IF (NITEMS.GT.2) CALL READI(NBFGSMAX1)
            IF (NITEMS.GT.3) CALL READI(NBFGSMAX2)
            IF (NITEMS.GT.4) CALL READF(CEIG)
            IF (NITEMS.GT.5) CALL READI(NEVL)
            BFGSTST=.TRUE.
! 
! Tolerance for eigenvalue % change for convergence to be allowed in Rayleigh-Ritz
! procedure. Default 1%.
! 
         ELSE IF (WORD.EQ.'BFGSTSPC') THEN
            CALL READF(EVPC)
! 
! Tolerance for eigenvector overlap in BFGSTS where the number of tangent space
! steps switches from small to large. 0.0001 was the traditional value (default).
! 
         ELSE IF (WORD.EQ.'BFGSTSTOL') THEN
            CALL READF(BFGSTSTOL)
! 
! Debug for basin-hopping interpolation
! 
         ELSE IF (WORD.EQ.'BHDEBUG') THEN
            BHDEBUG=.TRUE.
! 
! Parameters for basin-hopping interpolation
! 
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
! 
! Additional parameter for basin-hopping interpolation.
! Save the lowest energy minimum, rather than the lowest with the true PE plus spring energy.
! 
! 
         ELSE IF (WORD.EQ.'BHINTERPUSELOWEST') THEN
            BHINTERPUSELOWEST=.TRUE.
            IF (NITEMS.GT.1) CALL READA(UNSTRING)
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'CHECKENER') BHCHECKENERGYT=.TRUE.
            IF (NITEMS.GT.2) CALL READI(BHSTEPSMIN)
! 
! Binary LJ parameters for use with the LP or LS atom types.
! 
         ELSE IF (WORD.EQ.'BINARY') THEN
            BINARY=.TRUE.
            CALL READI(NTYPEA)
            CALL READF(EPSAB)
            CALL READF(EPSBB)
            CALL READF(SIGAB)
            CALL READF(SIGBB)
! 
! Parameters for bisection runs
! 
         ELSE IF (WORD.EQ.'BISECT') THEN
            BISECTT=.TRUE.
            CALL READF(BISECTMINDIST)
            CALL READF(BISECTMAXENERGY)
            CALL READI(BISECTSTEPS)
            CALL READI(BISECTMAXATTEMPTS)
            IF (NITEMS.GT.5) CALL READA(UNSTRING)
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'ICINTERP') ICINTERPT=.TRUE.
! 
! Debug printing for BISECT runs.
! 
         ELSE IF (WORD.EQ.'BISECTDEBUG') THEN
            BISECTDEBUG=.TRUE.
         ELSE IF (WORD.EQ.'BOND') THEN
            NUBONDS = NUBONDS + 1
            CALL READI(UBONDS(NUBONDS,1))
            CALL READI(UBONDS(NUBONDS,2))
! 
! General BLN model.
! 
         ELSE IF (WORD.EQ.'BLN') THEN
            BLNT=.TRUE.
            CALL READF(RK_R)
            CALL READF(RK_THETA)
            ALLOCATE(BEADLETTER(NATOMS),BLNSSTRUCT(NATOMS),
     &      LJREP_BLN(NATOMS,NATOMS),LJATT_BLN(NATOMS,NATOMS),A_BLN(NATOMS),B_BLN(NATOMS),C_BLN(NATOMS),D_BLN(NATOMS))
            OPEN(UNIT=100,FILE='BLNsequence',STATUS='OLD')
            READ(100,*) DUMMYCH
            READ(100,*) LJREPBB, LJATTBB
            READ(100,*) LJREPLL, LJATTLL
            READ(100,*) LJREPNN, LJATTNN
            READ(100,*) DUMMYCH
            READ(100,*) DUMMYCH
            READ(100,*) HABLN, HBBLN, HCBLN, HDBLN
            READ(100,*) EABLN, EBBLN, ECBLN, EDBLN
            READ(100,*) TABLN, TBBLN, TCBLN, TDBLN
            DO J1=1,NATOMS-1
               READ(100,'(A1)',ADVANCE='NO') BEADLETTER(J1)
            ENDDO
            READ(100,'(A1)') BEADLETTER(NATOMS) ! this line is needed to advance the input line for the next read
            DO J1=1,NATOMS-3
               READ(100,'(A1)',ADVANCE='NO') BLNSSTRUCT(J1)
            ENDDO
            CLOSE(100)
            PRINT '(A,I8,A)','BLN sequence of ',NATOMS,' beads read:'
            WRITE(*,'(A1)',ADVANCE='NO') BEADLETTER(1:NATOMS)
            PRINT '(A)',' '
            PRINT '(A,I8,A)','BLN dihedral types:'
            WRITE(*,'(A1)',ADVANCE='NO') BLNSSTRUCT(1:NATOMS-3)
            PRINT '(A)',' '
            PRINT '(A,2F15.5)','B-B LJ coefficients: ',LJREPBB, LJATTBB
            PRINT '(A,2F15.5)','L-L LJ coefficients: ',LJREPLL, LJATTLL
            PRINT '(A,2F15.5)','N-N LJ coefficients: ',LJREPNN, LJATTNN
            PRINT '(A,4F15.5)','Helix    dihedral coefficients: ',HABLN,HBBLN,HCBLN,HDBLN
            PRINT '(A,4F15.5)','Extended dihedral coefficients: ',EABLN,EBBLN,ECBLN,EDBLN
            PRINT '(A,4F15.5)','Turn     dihedral coefficients: ',TABLN,TBBLN,TCBLN,TDBLN
            call param_arrayBLN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,
     &      LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN,
     &      HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN, NATOMS)
! 
! General BLN/Go model.
! 
         ELSE IF (WORD.EQ.'BLNGO') THEN
            BLNT=.TRUE.
            GOTYPE=.TRUE.
            CALL READF(RK_R)
            CALL READF(RK_THETA)
            IF (NITEMS.GT.3) THEN
               CALL READF(GOFACTOR)
            ENDIF
            ALLOCATE(BEADLETTER(NATOMS),BLNSSTRUCT(NATOMS),
     &      LJREP_BLN(NATOMS,NATOMS),LJATT_BLN(NATOMS,NATOMS),A_BLN(NATOMS),B_BLN(NATOMS),C_BLN(NATOMS),D_BLN(NATOMS))
            OPEN(UNIT=100,FILE='BLNsequence',STATUS='OLD')
            READ(100,*) DUMMYCH
            READ(100,*) LJREPBB, LJATTBB
            READ(100,*) LJREPLL, LJATTLL
            READ(100,*) LJREPNN, LJATTNN
            READ(100,*) DUMMYCH
            READ(100,*) DUMMYCH
            READ(100,*) HABLN, HBBLN, HCBLN, HDBLN
            READ(100,*) EABLN, EBBLN, ECBLN, EDBLN
            READ(100,*) TABLN, TBBLN, TCBLN, TDBLN
            DO J1=1,NATOMS-1
               READ(100,'(A1)',ADVANCE='NO') BEADLETTER(J1)
            ENDDO
            READ(100,'(A1)') BEADLETTER(NATOMS) ! this line is needed to advance the input line for the next read
            DO J1=1,NATOMS-3
               READ(100,'(A1)',ADVANCE='NO') BLNSSTRUCT(J1)
            ENDDO
            CLOSE(100)
            PRINT '(A,I8,A)','BLN sequence of ',NATOMS,' beads read:'
            WRITE(*,'(A1)',ADVANCE='NO') BEADLETTER(1:NATOMS)
            PRINT '(A)',' '
            PRINT '(A,I8,A)','BLN dihedral types:'
            WRITE(*,'(A1)',ADVANCE='NO') BLNSSTRUCT(1:NATOMS-3)
            PRINT '(A)',' '
            PRINT '(A,2F15.5)','B-B LJ coefficients: ',LJREPBB, LJATTBB
            PRINT '(A,2F15.5)','L-L LJ coefficients: ',LJREPLL, LJATTLL
            PRINT '(A,2F15.5)','N-N LJ coefficients: ',LJREPNN, LJATTNN
            PRINT '(A,4F15.5)','Helix    dihedral coefficients: ',HABLN,HBBLN,HCBLN,HDBLN
            PRINT '(A,4F15.5)','Extended dihedral coefficients: ',EABLN,EBBLN,ECBLN,EDBLN
            PRINT '(A,4F15.5)','Turn     dihedral coefficients: ',TABLN,TBBLN,TCBLN,TDBLN
            call param_arrayBLN(LJREP_BLN,LJATT_BLN,A_BLN,B_BLN,C_BLN,D_BLN,BEADLETTER,BLNSSTRUCT,
     &      LJREPBB, LJATTBB, LJREPLL, LJATTLL, LJREPNN, LJATTNN,
     &      HABLN, HBBLN, HCBLN, HDBLN, EABLN, EBBLN, ECBLN, EDBLN, TABLN, TBBLN, TCBLN, TDBLN, NATOMS)
! 
! Yimin Wang and Joel Bowman's water potential (2010)
! 
         ELSE IF (WORD.EQ.'BOWMAN') THEN
            BOWMANT=.TRUE.
            CALL READI(BOWMANPES)
            CALL READA(BOWMANDIR)

            ! sn402: added (see comments at keyword FRQCONV)
            FRQCONV = 5.123934D14
            WRITE(*,*) "keywords> Set frequency conversion factor to the SD default value: ", FRQCONV
            WRITE(*,*) "keywords> This corresponds to frequencies being given in radians/s"
! 
! BSMIN calculates a steepest-descent path using gradient only information
! with convergence criterion GMAX for the RMS force and initial precision
! EPS. The Bulirsch-Stoer algorithm is used.
! 
         ELSE IF (WORD.EQ.'BSMIN') THEN
            BSMIN=.TRUE.
            IF (NITEMS.GT.1) CALL READF(GMAX)
            IF (NITEMS.GT.2) CALL READF(EPS)

         ELSE IF (WORD.EQ.'BULK') THEN
            BULKT=.TRUE.
            IF (NITEMS.GT.1) THEN
               CALL READF(bulk_boxvec(1))
            ENDIF
            IF (NITEMS.GT.2) THEN
               CALL READF(bulk_boxvec(2))
            ENDIF
            IF (NITEMS.GT.3) THEN
               CALL READF(bulk_boxvec(3))
            ENDIF
!
! Allow particles move out of the box when dealing with bulk system.
!
         ELSE IF (WORD.EQ.'BULKBOX') THEN
            BULKBOXT=.TRUE.
! 
! CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CADPAC tells the program to read derivative information in
! CADPAC format.                                        - default FALSE
! 
         ELSE IF (WORD.EQ.'CADPAC') THEN
            CADPAC=.TRUE.
            CALL READA(SYS)
            DO J1=1,80
               IF (SYS(J1:J1).EQ.' ') THEN
                  LSYS=J1-1
                  GOTO 10
               ENDIF
            ENDDO
10          IF (NITEMS.GT.2) THEN
               CALL READA(EDITIT)
            ELSE
               EDITIT='editit.' // SYS(1:LSYS)
            ENDIF
         ELSE IF (WORD.EQ.'CALCDIHE') THEN
            CALCDIHE=.TRUE.

! hk286
! ELSE IF (WORD.EQ.'CALCFREQSRB') THEN
! PRINT *, "CALCFREQSRB > Computing the frequencies"
! OPEN(UNIT = 28, FILE = 'coordsfreq', STATUS = 'OLD')
! DO J1 = 1, NATOMS
! READ(28, *) Q(3*J1-2:3*J1)
! ENDDO
! CLOSE(UNIT = 28)
! RBAANORMALMODET = .TRUE.
! CALL NRMLMD (Q, Q, .FALSE.)
! OPEN(UNIT = 28, FILE = 'freqs')
! WRITE(28, *) Q
! CLOSE(UNIT = 28)
! STOP

! 
! If READPATH is specified with CALCRATES then the rates are calculated from the
! information in an existing path.info file without any stationary point searches.
! A CONNECT or PATH run must be performed first unless READPATH is specified.
! 
         ELSE IF (WORD.EQ.'CALCRATES') THEN
            CALCRATES=.TRUE.
            IF (NITEMS.GT.1) CALL READF(TEMPERATURE)
            IF (NITEMS.GT.2) CALL READF(HRED)

! jbr36 classical rate calculations
         ELSE IF (WORD.EQ.'CLASSICALRATES') THEN
            CLASSICALRATEST=.TRUE.

         ELSE IF (WORD.EQ.'INE_NEW') THEN
            CLSTRINGTST=.TRUE.
            CLSTRINGT=.TRUE.
            EVOLVESTRINGT = .TRUE.
            CALL READF(STTSRMSCONV)
            CALL READI(ST_TSSTEP)
            CALL READF(LAN_DIST)
            CALL READI(LANSTEP)
            CALL READF(LANCONV)
            CALL READF(LANFACTOR)

         ELSE IF (WORD.EQ.'CLSTRING') THEN
            CLSTRINGT=.TRUE.
            EVOLVESTRINGT = .TRUE.
         ELSE IF (WORD.EQ.'COLLAGENOP') THEN
            COLLAGENOP=.TRUE.
            CALL READI(COLLINDICES(1))
            CALL READI(COLLINDICES(2))
            CALL READI(COLLINDICES(3))
            CALL READI(COLLINDICES(4))

         ELSE IF (WORD.EQ.'HESSREAD') THEN
            HESSREADT=.TRUE.

! jbr36 instanton (rate) calculations, more general than the ringpolymer keyword
         ELSE IF (WORD.EQ.'INSTANTONOPT') THEN
            INSTANTONOPTT=.TRUE.
            IF (NITEMS.GT.1) CALL READI(NIMAGEINST)
            IF (NITEMS.GT.2) CALL READF(DISTORTINST)
            IF (NITEMS.GT.3) CALL READF(DELTAINST)
            IF (NITEMS.GT.4) CALL READF(TEMPERATURE1)

         ELSE IF (WORD.EQ.'VARSTEPOPT') THEN
            VARSTEPOPTT=.TRUE.
         ELSE IF (WORD.EQ.'INSTANTONRATE') THEN
            INSTANTONRATET=.TRUE.
            IF (NITEMS.GT.1) CALL READI(NIMAGEINST)

         ELSE IF (WORD.EQ.'INSTANTONSTARTDUMP') THEN
            INSTANTONSTARTDUMPT=.TRUE.
            IF (NITEMS.GT.1) CALL READF(TEMPERATURE1)
! 
! CAMSHIFT calculates the NMR chemical shifts from atomic coordinates and is mostly used for
! chemical shift restrained simulations
! 
         ELSE IF (WORD .EQ. 'CAMSHIFT') THEN
            CHEMSHIFT=.TRUE.
            COUT = .FALSE.
            IF (NITEMS .LT. 4) THEN
               WRITE(*,*) 'CamShift version, path and shiftfile are needed for Camshift'
               STOP
            ELSE
               CALL READA(CSVERSION)
               CSVERSION=TRIM(ADJUSTL(CSVERSION))
               SELECT CASE(CSVERSION)
               CASE('MERGE','ORIGINAL','NOFF')
               CASE DEFAULT
               WRITE(*,*) TRIM(ADJUSTL(CSVERSION)),' is not a correct CamShift version'
               STOP
               END SELECT
               CALL READA(SVNROOT)
               SVNROOT=TRIM(ADJUSTL(SVNROOT))
               CSPATH=TRIM(ADJUSTL(SVNROOT))//'CAMSHIFTDATA/'
               CSPATH=TRIM(ADJUSTL(CSPATH))
               CALL READA(SHIFTFILE)
               SHIFTFILE=TRIM(ADJUSTL(SHIFTFILE))
               IF (NITEMS .GT. 4) THEN
                  CALL READF(CSN)
               ENDIF
               IF (NITEMS .GT. 5) THEN
                  CALL READF(CSALPHA)
               ENDIF
            ENDIF
! 
! Double-ended connection keyword for ts candidates.
! 
         ELSE IF (WORD == 'CANDIDATES') THEN
            CALL READA(CANDIDATES)
! 
! Virus capsid specification.
! 
         ELSE IF (WORD.EQ.'CAPSID') THEN
            RIGIDBODY=.TRUE.
            ANGLEAXIS=.TRUE.
            HEIGHT=0.5D0
            CALL READF(CAPSRHO)
            CALL READF(CAPSEPS2)
            CALL READF(CAPSRAD)
            IF (NITEMS.GT.4) CALL READF(HEIGHT)
         ELSE IF (WORD.EQ.'CAPSID2') THEN
! RIGIDBODY=.TRUE.
            ANGLEAXIS2=.TRUE.
            HEIGHT=0.5D0
            CALL READF(CAPSRHO)
            CALL READF(CAPSEPS2)
            CALL READF(CAPSRAD)
            IF (NITEMS.GT.4) CALL READF(HEIGHT)

! starting from a given residue, use cartesians for everything
         ELSE IF (WORD.EQ.'CARTRESSTART') THEN
            CALL READI(CARTRESSTART)
! 
! CASTEP tells the program to read derivative information in
! CASTEP format.                                        - default FALSE
! 
         ELSE IF ((WORD.EQ.'CASTEP').OR.(WORD.EQ.'CASTEPC')) THEN
            CASTEP=.TRUE.
            IF (WORD.EQ.'CASTEP') DFTP=.TRUE.
            IF (NITEMS.GT.2) THEN
               CALL READA(CASTEPJOB)
               CALL READA(SYS)
               CASTEPJOB=TRIM(ADJUSTL(CASTEPJOB)) // ' ' // TRIM(ADJUSTL(SYS))
            ELSE
               WRITE(*,'(A)') 'keywords> ERROR - CASTEP job or system unspecified'
               CALL FLUSH(6)
               STOP
            ENDIF
            DO J1=1,80
               IF (SYS(J1:J1).EQ.' ') THEN
                  LSYS=J1-1
                  GOTO 22
               ENDIF
            ENDDO
22          CONTINUE

            ! sn402
            WRITE(*,*) "keywords> WARNING: there is currently no default frequency conversion set for CASTEP"
            WRITE(*,*) "keywords> Log products of frequencies will be computed in internal units."
            WRITE(*,*) "keywords> To learn how to set a default conversion factor, check the comments for 
     &                            the FRQCONV keyword in keywords.f"


! 
! charmm stuff (DAE)
! 
         ELSE IF (WORD.EQ.'CHARMM') THEN
            CHRMMT=.TRUE.
            IF (.NOT.CISTRANS) THEN
               NOCISTRANS=.TRUE.
               CHECKOMEGAT=.TRUE.
            ENDIF
            IF (.NOT.TURNOFFCHECKCHIRALITY) CHECKCHIRALT=.TRUE.

            IF ((PERMDIST.OR.LOCALPERMDIST.OR.LPERMDIST)) THEN
               IF (NPERMSIZE(1).EQ.NATOMS) THEN
                  WRITE(*,*)'keyword>ERROR-PERMDIST/LOCALPERMDIST/LPERMDIST is specfied for CHARMM, but there is no perm.allow file'
                  STOP
               ENDIF
            ENDIF
! CALL CHALLOCATE(NATOMS) ! this looks like a bug!! DJW
            ALLOCATE(ATMASS(NATOMS))
            IF (MACHINE) THEN
               ! SAT: we will read in the coords ourselves and pass them to CHARMM

               ! --- start ---
               ! read in the coords
               INQUIRE(IOLENGTH=J1) (Q(J),J=1,3*NATOMS)
               IF (FILTH2==0) THEN
                  OTEMP='points1.inp'
               ELSE
                  WRITE(OTEMP,*) FILTH2
                  OTEMP='points1.inp.'//TRIM(ADJUSTL(OTEMP))
               ENDIF
               OPEN(113,FILE=OTEMP,ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='OLD',RECL=J1)
               READ(113,REC=1) (Q(J),J=1,3*NATOMS)
               IF (MAXVAL(Q)==0.0D0) THEN
                  PRINT *, 'Zero coordinates - stop'
                  CALL FLUSH(6)
                  STOP
               ENDIF
               CLOSE(113)
               ! --- end ---
               ! SAT: line below was intended to replace the block of code above
               ! (marked); unfortunately, due to the miscompilation with pgi this
               ! does not work. The compiler does not really want to reuse the
               ! code. Sigh...
               ! call ReadInpFile(Q)

               ! save them into CH. arrays and pass to CHARMM
               DO J1=1,NATOMS
                  CHX(J1)=Q(3*(J1-1)+1)
                  CHY(J1)=Q(3*(J1-1)+2)
                  CHZ(J1)=Q(3*(J1-1)+3)
               ENDDO
               CALL CHSETUP(CHX,CHY,CHZ,CHMASS,NATOM,TOPFILE,PARFILE,DATA_UNIT)
               ! CALL FILLICT(CHX,CHY,CHZ,DUMMY1,.TRUE.)
               CALL FILLICTABLE(Q)
            ELSE
               ! charmm will read the coords and will return them to OPTIM via CH. vecs
               CHX(1)=13.13d13 ! this way we will tell CHARMM to save its coords into CH. arrays; otherwise it will
               CALL CHSETUP(CHX,CHY,CHZ,CHMASS,NATOM,TOPFILE,PARFILE,DATA_UNIT)
               ! 
            ENDIF ! SAT
            CALL CHSETZSYMATMASS
            IF (FILTH.NE.0) THEN
               OPEN(UNIT=20,FILE='coords.read',STATUS='REPLACE')
               CLOSE(20)
            ENDIF
! NATOMS=NATOM  ! should already know NATOMS from getparams
            IF (NATOM /= NATOMS) THEN
               WRITE(*,'(A)') 'No. of atoms in "input.crd" and file specified in CHARMM part of odata conflict'
               PRINT *, 'NATOM,NATOMS=',NATOM, NATOMS
               CALL FLUSH(6)
               STOP
            ENDIF
            CALL CHALLOCATE(NATOMS)
! csw34> This is where all the internal coordinates are set up, and the
! different types of dihedrals identified. We also identify dihedrals
! that are twistable.
            CALL CHSETDIHE
! csw34> If FREEZERES specified, call CHRESTOATOM to populate the
! FROZEN array (from ocharmm.src)
            IF (FREEZERES) CALL CHRESTOATOM(FROZENRES,FROZEN)

            IF (CONNECTT) CALL CHSETSEED
! IF (CALCDIHE) CALL READREF(NATOMS)
            DO J1=1,NATOMS
               Q(3*(J1-1)+1)=CHX(J1)
               Q(3*(J1-1)+2)=CHY(J1)
               Q(3*(J1-1)+3)=CHZ(J1)
               ATMASS(J1) = CHMASS(J1)
               ! PRINT *,'ATMASS',ATMASS(J1)
            ENDDO
            IF (TWISTDIHET) THEN
               ! csw34> We have changed how the dihedral to select is choosen. Now, a
               ! random number between 0 and 1 is passed from PATHSAMPLE to OPTIM
               ! as part of the TWISTDIHE keyword. This number is then multiplied by
               ! the total number of twistable (not omega or chiral) dihedrals in the
               ! system and rounded up to give an index. The IICD of the dihedral to be twisted is
               ! then extracted from the array DIHETOTWIST which contains only those
               ! dihedrals. This is passed to TWISTDIHE as the DMODE arguement.
               ! 
               ! We use CEILING to round up the index of the dihedral. This prevents
               ! getting DMODE=0, and ensures uniform sampling.
               DMODE=DIHETOTWIST(CEILING(PSRANDOM*NTWISTABLE))
               WRITE(*,*) 'keywords> Twisting dihedral IICD=',DMODE
               CALL TWISTDIHE(Q,DMODE,DPERT)
            ENDIF
            IF (PERTDIHET) THEN
               CALL PERTDIHE(Q,CHPMIN,CHPMAX,CHNMIN,CHNMAX,ISEED)
            ENDIF
            IF (INTMINT) CALL GETNINT(NINTS)  ! DJW - this is OK because CHARMM is the last keyword!

            ! sn402: added (see comments at keyword FRQCONV)
            FRQCONV = 2.045483D13
            WRITE(*,*) "keywords> Set frequency conversion factor to the CHARMM default value: ", FRQCONV
            WRITE(*,*) "keywords> This corresponds to frequencies being given in radians/s"

! 
! csw34> If using the CHARMM SCC-DFTB potential, we assume that all
! atoms are QM. If you are using a mixed QM/MM system, you should either
! not use the CHARMMDFTB keyword, or re-code it to check for fully QM
! systems. This keyword essentially prevents unnessesary printing!
! 
         ELSE IF (WORD.EQ.'CHARMMDFTB') THEN
            CHARMMDFTBT=.TRUE.
            WRITE(*,'(A)') 'keywords> WARNING - All atoms assumed to be QM, NBONDS calls disabled'
         ELSE IF (WORD.EQ.'CHARMMNOTUPDATE') THEN
            CHARMMNOTUPDATE=.TRUE.
            WRITE(*,'(A)') 'keywords> Charmm update nonbond list in ocharmm is turned off'
         ELSE IF (WORD.EQ.'CHARMMTYPE') THEN
            IF (NITEMS.GT.1) THEN
               CALL READA(TOPFILE)
               TOPFILE=TRIM(ADJUSTL(TOPFILE))
            ENDIF
            IF (NITEMS.GT.2) THEN
               CALL READA(PARFILE)
               PARFILE=TRIM(ADJUSTL(PARFILE))
            ELSE
               WRITE(*,*) 'keywords> TOPFILE and PARFILE have to be defined for CHARMMTYPE'
               STOP
            ENDIF
            IF (TOPFILE(1:6).EQ."toph19") THEN
               CHARMMTYPE=2
            ELSEIF (TOPFILE(1:7).EQ."top_all") THEN
               CHARMMTYPE = 1
            ELSE
               WRITE(*,*) 'keywords> TOPFILE ', TRIM(ADJUSTL(TOPFILE)),' is not recognised by OPTIM'
               STOP
            ENDIF
            WRITE(*,'(A,I2)') 'CHARMMTYPE set to ',CHARMMTYPE
! 
! If CHDEBUG is on, CHARMM related debug messages are printed
! 
         ELSE IF (WORD.EQ.'CHDEBUG') THEN
            CHDEBUG=.TRUE.
            IF (NITEMS.GT.1) CALL READA(UNSTRING)
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'EDEBUG') EDEBUG=.TRUE.

! 
! CHARMM related keyword, avoids inversion around C_alpha
! also implemented to AMBER (sf344)
         ELSE IF (WORD.EQ.'CHECKCHIRALITY') THEN
            CHECKCHIRALT=.TRUE.
! 
! If CHECKINDEX is .TRUE. and the BFGSTS routine converges an attempt is
! made to count the number of negative Hessian eigenvalues using projection,
! orthogonalization and iteration. We also need the opportunity to change the
! parameters NEVL and NEVS within BFGSTS if BFGSTS isn t true.
! CHECKINDEX can also be used with BFGSMIN and should understand NOHESS too.
! 
         ELSE IF (WORD.EQ.'CHECKINDEX') THEN
            CHECKINDEX=.TRUE.
            IF (NITEMS.GT.1) THEN
               CALL READI(NEVS)
            ENDIF
            IF (NITEMS.GT.2) THEN
               CALL READF(CEIG)
            ENDIF
            IF (NITEMS.GT.3) THEN
               CALL READI(NEVL)
            ENDIF
! 
! If the index found by checkindex does not correspond to BFGSMIN or BFGSTS then
! CHECKCONT causes a pushoff along the eigenvector correpsonding to the softest
! undesired negative eigenvalue.
! 
         ELSE IF (WORD.EQ.'CHECKCONT') THEN
            CHECKCONT=.TRUE.
! 
! If CHECKNEGATIVET is true then in bfgsts we backtrack and reduce the maximum step
! size if the smallest non-zero eigenvalue is positive.
! 
         ELSE IF (WORD.EQ.'CHECKD') THEN
            CHECKDT = .TRUE.
            CALL READI(CHECKDID)

         ELSE IF (WORD.EQ.'CHECKNEGATIVE') THEN
            CHECKNEGATIVET=.TRUE.
! 
! If CHECKOVERLAPT is true then in bfgsts we backtrack and reduce the maximum step
! size if the overlap with the previous eigenvector is less tha MINOVERLAP in
! magnitude.
! 
         ELSE IF (WORD.EQ.'CHECKOVERLAP') THEN
            CHECKOVERLAPT=.TRUE.
            IF (NITEMS.GT.1) CALL READF(MINOVERLAP)
! 
! Parameters for recalculating the repulsions in INTCONSTRAINT
! 
         ELSE IF (WORD.EQ.'CHECKREP') THEN
            IF (NITEMS.GT.1) CALL READI(CHECKREPINTERVAL)
            IF (NITEMS.GT.2) CALL READF(CHECKREPCUTOFF)
! 
! Put the structure in `finish' (the second structure given as input to OPTIM)
! into optimal alignment with a reference structure (the first structure given as
! input to OPTIM). The alignment is done via rigid-body movements and by considering
! permutational isomerizations. The aligned `finish' coordinates are dumped and the
! minimized distance is printed.
! 
         ELSE IF (WORD.EQ.'CLOSESTALIGNMENT') THEN
            CLOSESTALIGNMENT=.TRUE.
            WRITE(*,*) 'Putting structures into closest alignment, then stopping'
! 
! Dynamic adjustment of QCI spring constant.
! 
         ELSE IF (WORD.EQ.'QCIADJUSTK') THEN
            CALL READI(QCIKADJUSTFRQ)
            CALL READF(QCIKADJUSTTOL)
            CALL READF(QCIKADJUSTFRAC)
            CALL READF(QCIKINTMIN)
            CALL READF(QCIKINTMAX)
            WRITE(*,'(A,I10,2(A,G20.10),A,2G15.5)') ' keywords> Adjusting QCI spring constant every ',QCIKADJUSTFRQ,
     &   ' steps by factor ',QCIKADJUSTFRAC,' for spacing deviation ',QCIKADJUSTTOL,'% min/max=',QCIKINTMIN,QCIKINTMAX

! 
! Check for internal minimum in constraint terms for INTCONSTRAINT
! 
         ELSE IF (WORD.EQ.'QCIINTREPMINSEP') THEN
            CALL READI(QCIINTREPMINSEP)
            WRITE(*,'(A,G20.10)') ' keywords> Minimum separation in atom index for internal minimum check in repulsion=',
     &                              QCIINTREPMINSEP  
! 
! Check for internal minimum in constraint terms for INTCONSTRAINT
! 
         ELSE IF ((WORD.EQ.'CONCONINT').OR.(WORD.EQ.'QCICONINT')) THEN
            CHECKCONINT=.TRUE.
            WRITE(*,'(A,G20.10)') ' keyword> Turning on terms for internal minima in constraints'
! 
! Scaling factor for internal minima in con or rep
! 
         ELSE IF ((WORD.EQ.'CONINT').OR.(WORD.EQ.'QCIINT')) THEN
            IF (NITEMS.GT.1) CALL READF(INTMINFAC)
            WRITE(*,'(A,G20.10)') ' keyword> Internal minima terms will be scaled by a factor of ',INTMINFAC
!
! Maximum active atoms in QCI procedure.
!
      ELSE IF (WORD.EQ.'QCIRESTART') THEN
         QCIRESTART=.TRUE.
!
! Stop before going on to DNEB phase.
!
      ELSE IF (WORD.EQ.'QCISTOP') THEN
         QCISTOP=.TRUE.
!
! Maximum active atoms in QCI procedure.
!
      ELSE IF (WORD.EQ.'QCIMAXACTIVE') THEN
         CALL READI(MAXNACTIVE)


      ELSE IF (WORD.EQ.'QCICYCLES') THEN
         QCICYCLEST=.TRUE.
         QCIDNEBT=.FALSE.
         CALL READF(QCICYCDIST)
         CALL READI(QCICYCN)
         WRITE(*,'(A,F8.2,A,I8,A)') ' keyword> DNEB for distance <' , QCICYCDIST , ' and after ', QCICYCN,' cycles'   
      
! 
! Absolute distance to allow before turning on constraint potential.
! 
         ELSE IF (WORD.EQ.'CONCUTABS') THEN
            CONCUTABST=.TRUE.
            CONCUTFRACT=.FALSE.
            IF (NITEMS.GT.1) CALL READF(CONCUTABS)
! 
! Fraction of constraint distance to allow before turning on constraint potential.
! 
         ELSE IF (WORD.EQ.'CONCUTFRAC') THEN
            CONCUTFRACT=.TRUE.
            CONCUTABST=.FALSE.
            IF (NITEMS.GT.1) CALL READF(CONCUTFRAC)
! 
! CHINTERPOLATE controls the interpolation for BHINTERP using CHARMM's primitive
! internal coordinates. The 1st argument has to be either BC or BI for the backbone
! interpolation with Cartesians and Internals, respectively. The 2nd argument
! has to be either SC or SI for the sidechain interpolation with Cartesians and
! Internals, respectively. If DNEB is given as 3rd argument, this interpolation scheme
! will be used for DNEB. If CHINTERPOLATE is not defined in the odata file the default is
! that DNEB and BHINTERP are done in Cartesians
! 
         ELSE IF (WORD.EQ.'CHINTERPOLATE') THEN
            CALL READA(UNSTRING)
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'BI') CHBIT=.TRUE.
            CALL READA(UNSTRING)
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'SI') ICINTERPT=.TRUE.
            IF (NITEMS.GT.3) THEN
               CALL READA(UNSTRING)
               IF (TRIM(ADJUSTL(UNSTRING)).EQ.'DNEB') CHICDNEB=.TRUE.
            ENDIF
! 
! If BHINTERPolation, and CHRIGID is set for the CHARMM potential, rigid body
! translation and rotation is applied to the peptides/proteins if more
! than one peptide/protein is prsent.
! 
         ELSE IF (WORD.EQ.'CHRIGID') THEN
            CHRIGIDT=.TRUE.
            CALL READF(PTRANS)
            CALL READF(TRANSMAX)
            CALL READF(PROT)
            CALL READF(ROTMAX)
! 
! CISTRANS is a CHARMM/AMBER9/NAB related keyword, which allows cis-trans isomerisation of the peptide bond .
! 
         ELSE IF (WORD.EQ.'CISTRANS') THEN
            CISTRANS=.TRUE.

! 
! Sometimes have to modify the cold fusion limit when using high electric fields
! 
         ELSE IF (WORD.EQ.'COLDFUSION') THEN
            IF (NITEMS.GT.1) call READF(COLDFUSIONLIMIT)
! 
! Connect initial minimum in odata to final minimum in file finish - maximum
! number of transiiton states=NCONNECT. Obsolete - use NEWCONNECT instead.
! 
         ELSE IF (WORD.EQ.'CONNECT') THEN
            CONNECTT=.TRUE.
            IF (NITEMS.GT.1) CALL READI(NCONNECT)
! 
! Constraint potential for interpolation between minima.
! 
         ELSE IF (WORD.EQ.'CONPOT') THEN
            CONPOTT=.TRUE.
            IF (NITEMS.GT.1) CALL READF(CPCONSTRAINTTOL)
            IF (NITEMS.GT.2) CALL READF(CPCONSTRAINTDEL)
            IF (NITEMS.GT.3) CALL READF(CPCONSTRAINTREP)
            IF (NITEMS.GT.4) CALL READF(CPCONSTRAINREPCUT)
            IF (NITEMS.GT.5) CALL READF(CPCONFRAC)
            IF (NITEMS.GT.6) CALL READI(CPCONSEP)
            IF (NITEMS.GT.7) CALL READI(CPREPSEP)
! 
! jmc unres
! Note also use some of the non-specific charmm keywords like INTMIN, NGUESS, TWISTTYPE etc...
! 
         ELSE IF (WORD.EQ.'CONSEC') THEN
            CONSECT=.TRUE.
            DO J1=1,(NITEMS-1)/2
               CALL READI(STARTRES(J1))
               CALL READI(ENDRES(J1))
            END DO
            IF (NITEMS.GT.21) WRITE(*,'(A)') 'Too many sections requested - please adapt code!'
            NUMSEC=(NITEMS-1)/2
            PRINT *,'CONSEC ',(STARTRES(J1),J1=1,10),(ENDRES(J1),J1=1,10), NUMSEC
! 
! CONVERGE n m INDEX/NOINDEX sets the convergence criteria for the maximum
! unscaled step and RMS force                     - default n=0.0001, m=0.000001
! or m < 0.00001 .AND. n < m*100000
! If NOINDEX is set the Hessian index isn t checked - the default is
! INDEX.
! 
         ELSE IF (WORD .EQ. 'CONVERGE') THEN
            CALL READF(CONVU)
            IF (NITEMS.GT.2) THEN
               CALL READF(CONVR)
            ENDIF
            IF (NITEMS.GT.3) THEN
               CALL READU(WORD)
               IF (WORD.EQ.'NOINDEX') INDEXT=.FALSE.
            ENDIF
! 
! Probably prints the copyright info?
! 
         ELSE IF (WORD == 'COPYRIGHT') THEN
            CALL COPYRIGHT
! CP2K tells the program to read derivative information in
! CP2K format.                                        - default FALSE
! 
         ELSE IF ((WORD.EQ.'CP2K').OR.(WORD.EQ.'CP2KC')) THEN
            CP2K=.TRUE.
            IF (WORD.EQ.'CP2K') DFTP=.TRUE.
            IF (NITEMS.GT.2) THEN
               CALL READA(CP2KJOB)
               CALL READA(SYS)
               CP2KJOB=TRIM(ADJUSTL(CP2KJOB)) // ' ' // TRIM(ADJUSTL(SYS))
            ELSE
               WRITE(*,'(A)') 'keywords> ERROR - no CP2K system specified'
               CALL FLUSH(6)
               STOP
            ENDIF
            DO J1=1,80
               IF (SYS(J1:J1).EQ.' ') THEN
                  LSYS=J1-1
                  GOTO 281
               ENDIF
            ENDDO
281         CONTINUE
! 
! CPMD tells the program to read derivative information in
! CPMD format.                                        - default FALSE
! 
         ELSE IF ((WORD.EQ.'CPMD').OR.(WORD.EQ.'CPMDC')) THEN
            CPMD=.TRUE.
            IF (WORD.EQ.'CPMDC') CPMDC=.TRUE.
            IF (NITEMS.GT.1) THEN
               CALL READA(SYS)
            ELSE
               WRITE(*,'(A)') ' ERROR - no CPMD system specified'
               CALL FLUSH(6)
               STOP
            ENDIF
            DO J1=1,80
               IF (SYS(J1:J1).EQ.' ') THEN
                  LSYS=J1-1
                  GOTO 12
               ENDIF
            ENDDO
12          CONTINUE
            CALL SYSTEM(' grep -c DUMMY ' // SYS(1:LSYS) // ' > temp ')
            OPEN(UNIT=7,FILE='temp',STATUS='OLD')
            READ(7,*) J1
            IF (J1.NE.1) THEN
               WRITE(*,'(A)') 'ERROR, no dummy line in CPMD input file'
               CALL FLUSH(6)
               STOP
            ENDIF
! 
! Option to specify a different CPMD executible
! 
         ELSE IF (WORD.EQ.'CPMD_COMMAND') THEN
            IF (NITEMS.GT.1) CALL READA(CPMD_COMMAND)
!
! cppneb switch
!
         ELSE IF (WORD.EQ."CPPNEB") THEN
            CPPNEBT=.TRUE.
!
!  Use NCPU's by starting up to NCPU's OPTIM jobs.
!
      ELSE IF (WORD.EQ.'CPUS') THEN
         PBST=.FALSE.
         IF (NITEMS.GT.1) CALL READI(NCPU)
! 
! CUBIC: maintains cubic supercell for PV calculations
! 
         ELSE IF (WORD.EQ.'CUBIC') THEN
            CUBIC=.TRUE.
! 
! For the growing string or evolving string double-ended
! transition state search methods, use a cubic spline interpolation between
! the image points.
! 
         ELSE IF (WORD.EQ.'CUBSPL') THEN
            CUBSPLT = .TRUE.
! 
! potential cutoff
! 
         ELSE IF (WORD.EQ.'CUTOFF') THEN
            CUTT = .TRUE.
            CALL READF(POTENTIAL_CUTOFF)
! 
! DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
! 
! 
! Add a decahedral field to the potential of magnitude FTD.
! 
         ELSE IF (WORD.EQ.'D5H') THEN
            FIELDT=.TRUE.
            D5HT=.TRUE.
            CALL READF(FD5H)

         ELSE IF (WORD.EQ.'DB') THEN

            DBPT   = .TRUE.
            RBAAT  = .TRUE.
            CALL READF(DBEPSBB)
            CALL READF(DBEPSAB)
            CALL READF(DBSIGBB)
            CALL READF(DBSIGAB)
            CALL READF(DBPMU)
            IF (NITEMS > 6) THEN
               CALL READF(EFIELD)
               EFIELDT = .TRUE.
            ENDIF

            NRBSITES = 3
            ALLOCATE(RBSITE(NRBSITES,3))

            NTSITES = NATOMS*NRBSITES/2

         ELSE IF (WORD.EQ.'DBTD') THEN

            DBPTDT = .TRUE.
            RBAAT  = .TRUE.
            CALL READF(DBEPSBB)
            CALL READF(DBEPSAB)
            CALL READF(DBSIGBB)
            CALL READF(DBSIGAB)
            CALL READF(DBPMU)
            IF (NITEMS > 6) THEN
               CALL READF(EFIELD)
               EFIELDT = .TRUE.
            ENDIF

            NRBSITES = 3
            ALLOCATE(RBSITE(NRBSITES,3))

            NTSITES = (NATOMS/2-1)*NRBSITES + 4

         ELSE IF (WORD.EQ.'DJWRB') THEN
            DJWRBT=.TRUE.
            CALL READI(DJWRBID)
            IF (.NOT.ALLOCATED(ATMASS)) ALLOCATE(ATMASS(NATOMS))
            ATMASS(1:NATOMS)=1.0D0
            IF (DJWRBID /= 1) THEN
               PRINT *, 'DJWRB id ',DJWRBID,' unknown'
               STOP
            ENDIF
            IF (NITEMS.GT.2) CALL READI(NHEXAMERS)
! adk44            
            IF (NITEMS.GT.3) CALL READF(RADPENT)
            IF (NITEMS.GT.4) CALL READF(CAPSIDRHO)
            IF (NITEMS.GT.5) CALL READF(CAPSIDEPS)
            IF (NITEMS.GT.6) CALL READF(SIGMAPENT)
            IF (NITEMS.GT.7) CALL READF(RADHEX)
            IF (NITEMS.GT.8) CALL READF(SIGMAHEX)
            IF (NITEMS.GT.9) CALL READF(SIGMAPH)

! 
! DCHECK  turns ON/OFF warnings about short interatomic distances
! default ON
! 
         ELSE IF (WORD.EQ.'DMBLPY') THEN

            DMBLPYT = .TRUE.
            RBAAT   = .TRUE.
            CALL READF(YEPS)
            CALL READF(YKAPPA)
            CALL READF(DBSIGBB)
            CALL READF(DBPMU)
            IF (NITEMS > 5) THEN
               CALL READF(EFIELD)
               EFIELDT = .TRUE.
            ENDIF

            NRBSITES = 3
            ALLOCATE(RBSITE(NRBSITES,3))

            NTSITES = NATOMS*NRBSITES/2


         ELSE IF (WORD.EQ.'DCHECK') THEN
            CALL READU(WW)
            IF (WW .EQ. 'ON' .OR. WW .EQ. ' ') THEN
               DCHECK=.TRUE.
            ELSE IF (WW .EQ. 'OFF') THEN
               DCHECK=.FALSE.
            ENDIF
! 
! DEBUG ON/OFF sets n=1 for EFSTEPS, VALUES, SUMMARY above     - default OFF
! 
         ELSE IF (WORD .EQ. 'DEBUG') THEN
            BHDEBUG=.TRUE.
            CALL READU(WW)
            IF (WW .EQ. 'ON' .OR. WW .EQ. ' ') THEN
               EFSTEPST=.TRUE.
               PGRAD=.TRUE.
               NGRADIENTS=1
               EFSTEPS=1
               NSUMMARY=1
               NVALUES=1
               DEBUG=.TRUE.
               PRINTOPTIMIZETS=.TRUE.
               DUMPNEBXYZ=.TRUE.
               DUMPINTXYZ=.TRUE.
               DUMPNEBPTS=.TRUE.
               DUMPNEBEOS=.TRUE.
               DUMPINTEOS=.TRUE.
            ENDIF

         ELSE IF (WORD.EQ.'DESMAXAVGE') THEN
! maximum average energy before double ended search method can stop
            CALL READF(DESMAXAVGE)


! maximum energy jump in one step
! for an image in a double-ended search method
            CALL READF(DESMAXEJUMP)
! 
! Produces extra printing for the double-ended
! transition state search method runs (DNEB, GS or ES).
! 
         ELSE IF (WORD.EQ.'DESMDEBUG') THEN
            DESMDEBUG = .TRUE.

         ELSE IF (WORD.EQ.'DESMINT') THEN
            DESMINT = .TRUE.
            INTINTERPT = .FALSE. ! desmint and intinterp are mutually exclusive
            NATINT = .TRUE. ! must use natural internals for double ended search
! 
! DFTBT tells the program to call dftb for Tiffany s tight-binding.
! default FALSE
         ELSE IF (WORD.EQ.'DFTB') THEN
            DFTBT=.TRUE.
! 
! Initial diagonal elements for LBFGS
! 
         ELSE IF (WORD.EQ.'DGUESS') THEN
            CALL READF(DGUESS)
            IF (NITEMS.GT.2) CALL READF(XDGUESS)
            IF (NITEMS.GT.3) CALL READF(NEBDGUESS)
            IF (NITEMS.GT.4) CALL READF(INTDGUESS)
            IF (NITEMS.GT.5) CALL READF(GSDGUESS)
! 
! If DIJKSTRA is true then decide in newconnect uses Dijkstra;s algorithm in
! deciding which connections to try next.
! First argument on DIJKSTRA line controls the cost function. SAT
! 
         ELSE IF (WORD.EQ.'DIJKSTRA') THEN
            IF (NITEMS.GT.1) THEN
               CALL READU(WW)
               IF (TRIM(ADJUSTL(WW))=='EXP') THEN
                  EXPCOSTFUNCTION = .TRUE.
               ELSEIF (TRIM(ADJUSTL(WW))=='INDEX') THEN
                  INDEXCOSTFUNCTION = .TRUE.
               ELSEIF (trim(adjustl(WW))=='INTERP') THEN
                  INTERPCOSTFUNCTION = .TRUE.
                  CALL READF(INTERPDIFF)
                  CALL READU(WW)
                  IF (TRIM(ADJUSTL(WW))=='EXP') THEN
                     EXPCOSTFUNCTION = .TRUE.
                  ELSE
                     ! CALL READI(COSTFUNCTIONPOWER)
                     READ(WW,'(I20)') COSTFUNCTIONPOWER
                  ENDIF
               ELSE IF (WW(1:1) /= ' ') THEN
                  READ(WW,'(I20)') COSTFUNCTIONPOWER
               ENDIF
               IF (NITEMS.GT.2) THEN
                  CALL READU(WW)
                  IF (trim(adjustl(WW))=='INTDISTANCE') THEN
                     IF (.NOT.INTINTERPT) THEN
                        PRINT*, "INTDISTANCE doesn,t work without INTINTERP"
                        PRINT*, "specify the latter before DIJKSRA in odata"
                     ELSE
                        INTDISTANCET = .TRUE.
                     ENDIF
                  ELSE
                     READ(WW,*) DIJKSTRADMAX
                  ENDIF
               ENDIF
            ENDIF
! 
! DIJKSTRALOCAL specifies an adjustable factor used to multiply the
! distances between minima found within one DNEB cycle. Decreasing
! this metric will encourage attempts to complete the connection, which
! might otherwise never be tried if shorter distances exist. We are
! trying to correct for the imperfect nature of the distance criterion
! used for the DIJKSTRA metric in choosing new connection pairs.
! 
         ELSE IF (WORD.EQ.'DIJKSTRALOCAL') THEN
            CALL READF(DIJKSTRALOCAL)
! 
! Double well potential between first two atoms
! 
         ELSE IF (WORD.EQ.'DOUBLE') THEN
            DOUBLET=.TRUE.
! 
! DNEB convergence condition on energy per image.
! The energy must be at least a fraction DNEBEFRAC of the lowest value
! attained. Extra iterations will be allowed to recover from excursions
! to higher energy.
! 
         ELSE IF (WORD.EQ.'DNEBEFRAC') THEN
            CALL READF(DNEBEFRAC)
! 
! DNEB RMS threshold for switching to NEB
! 
         ELSE IF (WORD.EQ.'DNEBSWITCH') THEN
            CALL READF(DNEBSWITCH)
! 
! Strings keyword.
! 
         ELSE IF (WORD.EQ.'DQAGKEY') THEN
            CALL READI(DQAGKEY)
! 
! Obsolete: create a trajectory between the endpoints by increasing
! a spring constant. Sounds like MD steering, but it doesn`t actually
! work very well!
! 
         ELSE IF (WORD.EQ.'DRAG') THEN
            DRAGT=.TRUE.
! 
! DUMPALLPATHS prints a summary of all min-sad-min triples produced by NEWCONNECT to
! file path.info. For each stationary point the energy, point group order and symbol,
! Hessian eigenvalues and coordinates are given. Hessian eigenvalues are computed
! if not yet calculated, otherwise they are saved during the CONNECT process.
! 
         ELSE IF (WORD.EQ.'DUMPALLPATHS') THEN
            DUMPALLPATHS=.TRUE.
            IF (FILTH.EQ.0) THEN
               WRITE(PINFOSTRING,'(A9)') 'path.info'
            ELSE
               WRITE(PINFOSTRING,'(A)') 'path.info.'//TRIM(ADJUSTL(FILTHSTR))
            ENDIF
            IF (MACHINE) THEN
               OPEN(UNIT=88,FILE=PINFOSTRING,STATUS='UNKNOWN',FORM='UNFORMATTED')
            ELSE
               OPEN(UNIT=88,FILE=PINFOSTRING,STATUS='UNKNOWN')
            ENDIF
! 
! DUMPBESTPATH only writes the minima and TS on the best path as found in NEWCONNECT
! to path.info. The format is identical to DUMPALLPATHS.
! 
         ELSE IF (WORD.EQ.'DUMPBESTPATH') THEN
            DUMPBESTPATH=.TRUE.
            IF (DUMPALLPATHS.OR.DUMPPATH) THEN
               WRITE(*,'(A)') ' keyword> Use of DUMPBESTPATH only without DUMPALLPATHS and DUMPPATH'
               STOP
            ENDIF
            IF (FILTH.EQ.0) THEN
               WRITE(PINFOSTRING,'(A9)') 'path.info'
            ELSE
               WRITE(PINFOSTRING,'(A)') 'path.info.'//TRIM(ADJUSTL(FILTHSTR))
            ENDIF
            IF (MACHINE) THEN
               OPEN(UNIT=88,FILE=PINFOSTRING,STATUS='UNKNOWN',FORM='UNFORMATTED')
            ELSE
               OPEN(UNIT=88,FILE=PINFOSTRING,STATUS='UNKNOWN')
            ENDIF
! 
! Creates a file in pathsample min.data format for the minimum found
! following a minimisation. Useful for a DPS initial path run in
! creating entries for the two endpoints.
! Can also be used with BHINTERP alone to generate a list of entries
! for interpolated minima.
! 
         ELSE IF (WORD.EQ.'DUMPDATA') THEN
            DUMPDATAT=.TRUE.
            IF (FILTH.EQ.0) THEN
               WRITE(PINFOSTRING,'(A13)') 'min.data.info'
            ELSE
               WRITE(PINFOSTRING,'(A)') 'min.data.info.'//TRIM(ADJUSTL(FILTHSTR))
            ENDIF
            IF (MACHINE) THEN
               OPEN(UNIT=881,FILE=PINFOSTRING,STATUS='UNKNOWN',FORM='UNFORMATTED')
            ELSE
               OPEN(UNIT=881,FILE=PINFOSTRING,STATUS='UNKNOWN')
            ENDIF

! Write output coords of minima in an unformatted file, min.points. This
! has a similar effect to DUMPDATA.AND.MACHINE, but we can be more
! flexible about the input format (e.g. can use MULTIJOB and/or formatted
! input). Also, the header lines will not be written to the unformatted
! file (unlike DUMPDATA.AND.MACHINE). They are written to min.data.info
! as normal.
         ELSE IF (WORD.EQ.'DUMPDATA_MACHINE') THEN
            CALL READI(NDOF)
            DUMPDATAT=.TRUE.
            DUMPDATA_MACHINET=.TRUE.
            IF (FILTH.EQ.0) THEN
               WRITE(PINFOSTRING,'(A13)') 'min.data.info'
               WRITE(MINPOINTSSTRING,'(A10)') 'min.points'
            ELSE
               WRITE(PINFOSTRING,'(A)') 'min.data.info.'//TRIM(ADJUSTL(FILTHSTR))
               WRITE(MINPOINTSSTRING,'(A)') 'min.points.'//TRIM(ADJUSTL(FILTHSTR))
            ENDIF
            
!            OPEN(UNIT=882,FILE=MINPOINTSSTRING,STATUS='UNKNOWN',FORM='UNFORMATTED')
            OPEN(UNIT=882,FILE=MINPOINTSSTRING,ACCESS='DIRECT',STATUS='UNKNOWN',FORM='UNFORMATTED', RECL=8*NDOF)
            RECCOUNT = 0  ! Need this if we're using multijob or similar: counts how many coords have already been written.
            OPEN(UNIT=881,FILE=PINFOSTRING,STATUS='UNKNOWN')

! 
! Explicit dump of interpolation EofS for intlbfgs. Should be set .TRUE. if DEBUG is set.
! 
         ELSE IF (WORD == 'DUMPINTEOS') THEN
            DUMPINTEOS=.TRUE.
            IF (NITEMS>1) CALL READI(DUMPINTEOSFREQ)
! 
! Explicit dump of EofS.neb for DNEB. Should be set .TRUE. if DEBUG is set.
! 
         ELSE IF (WORD == 'DUMPNEBEOS') THEN
            DUMPNEBEOS=.TRUE.
            IF (NITEMS>1) CALL READI(DUMPNEBEOSFREQ)
! 
! Explicit dump of something for DNEB. Should be set .TRUE. if DEBUG is set.
! 
         ELSE IF (WORD == 'DUMPNEBPTS') THEN
            DUMPNEBPTS=.TRUE.
            IF (NITEMS>1) CALL READI(DUMPNEBPTSFREQ)
! 
! Explicit dump of image coordinates in xyz format for intlbfgs. Should
! be set .TRUE. if DEBUG is set.
! 
         ELSE IF (WORD == 'DUMPINTXYZ') THEN
            DUMPINTXYZ=.TRUE.
            IF (NITEMS>1) CALL READI(DUMPINTXYZFREQ)
! 
! Explicit dump of image coordinates in xyz format for DNEB. Should
! be set .TRUE. if DEBUG is set.
! 
         ELSE IF (WORD == 'DUMPNEBXYZ') THEN
            DUMPNEBXYZ=.TRUE.
            IF (NITEMS>1) CALL READI(DUMPNEBXYZFREQ)
! 
! DUMPPATH prints a summary of a min-sad-min-...-min path produced by CONNECT to
! file path.info. For each stationary point the energy, point group order and symbol,
! Hessian eigenvalues and coordinates are given. Hessian eigenvalues are computed
! if not yet calculated, otherwise they are saved during the CONNECT process.
! 
         ELSE IF (WORD.EQ.'DUMPPATH') THEN
            DUMPPATH=.TRUE.
            IF (FILTH.EQ.0) THEN
               WRITE(PINFOSTRING,'(A9)') 'path.info'
            ELSE
               WRITE(PINFOSTRING,'(A)') 'path.info.'//TRIM(ADJUSTL(FILTHSTR))
            ENDIF
            IF (MACHINE) THEN
               OPEN(UNIT=88,FILE=PINFOSTRING,STATUS='UNKNOWN',FORM='UNFORMATTED')
            ELSE
               OPEN(UNIT=88,FILE=PINFOSTRING,STATUS='UNKNOWN')
            ENDIF
! 
! If DUMPSP is true then OPTIM will dump minima and ts data in the pathsample format.
! The resulting files could be merged with MERGEDB.
! 
         ELSE IF (WORD.EQ.'DUMPSP') THEN
            DUMPSP=.TRUE.
! 
! DUMPVECTOR switches on dumping of eigenvectors to file
! vectors.dump                                     - default OFF
! ALLSTEPS dumps the vector(s) at each step. ALLVECTORS dumps all the vectors.
! The defaults are for only the vector corresponding to the softest non-zero
! eigenvalue to be dumped for the last step.
! 
          ELSE IF (WORD .EQ. 'DUMPMAG') THEN
            DUMPMAG=.TRUE.
          ELSE IF (WORD .EQ. 'DUMPVECTOR') THEN
            DUMPV=.TRUE.
            IF (NITEMS.GT.1) THEN
               CALL READU(WORD)
               IF (WORD.EQ.'ALLSTEPS') ALLSTEPS=.TRUE.
               IF (WORD.EQ.'ALLVECTORS') ALLVECTORS=.TRUE.
               IF (WORD.EQ.'MWVECTORS') MWVECTORS=.TRUE.
            ENDIF
            IF (NITEMS.GT.2) THEN
               CALL READU(WORD)
               IF (WORD.EQ.'ALLSTEPS') ALLSTEPS=.TRUE.
               IF (WORD.EQ.'ALLVECTORS') ALLVECTORS=.TRUE.
               IF (WORD.EQ.'MWVECTORS') MWVECTORS=.TRUE.
            ENDIF
            IF (NITEMS.GT.3) THEN
               CALL READU(WORD)
               IF (WORD.EQ.'ALLSTEPS') ALLSTEPS=.TRUE.
               IF (WORD.EQ.'ALLVECTORS') ALLVECTORS=.TRUE.
               IF (WORD.EQ.'MWVECTORS') MWVECTORS=.TRUE.
            ENDIF
! 
! EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE
! 
! 
! EDIFFTOL specifies the maximum energy difference between permutational isomers in connect.
! 
         ELSE IF (WORD.EQ.'EDIFFTOL') THEN
            CALL READF(EDIFFTOL)
! 
! Specify an electric field in the z-direction, units are V/A
! So far only implemented for use with TIPnP potentials
! 
         ELSE IF (WORD.EQ.'EFIELD') THEN
            IF (NITEMS.GT.1) CALL READF(EFIELD)
! 
! EFSTEPS n print the unscaled steps calculated for each mode
! every n cycles                                       - default OFF
! 
         ELSE IF (WORD .EQ. 'EFSTEPS') THEN
            EFSTEPST=.TRUE.
            CALL READI(EFSTEPS)

! 
! Locate the lowest eigenvector using xmylbfgs
! 
         ELSE IF (WORD.EQ.'EIGENONLY') THEN
            EIGENONLY=.TRUE.
            IF (NITEMS.GT.1) CALL READF(CEIG)
            IF (NITEMS.GT.2) THEN
               CALL READF(COVER)
               OVERCONV=.TRUE.
            ENDIF

! 
! Calculate analytical Hessian and normal mode frequencies at end of run.
! ENDHESS is only intended for use in single geometry optimisations, and
! should not be needed for CONNECT or PATH runs if DUMPPATH is specified.
! If the argument NENDHESS is omitted then all the eigenvalues are
! calculated - otherwise just the lowest NENDHESS.
! 
         ELSE IF (WORD.EQ.'ENDHESS') THEN
            ENDHESS=.TRUE.
            IF (NITEMS.GT.1) CALL READI(NENDHESS)
! 
! Calculate numerical Hessian and normal mode frequencies at end of run.
! Required if DUMPPATH or ENDHESS is specified for an UNRES run,
! in which case it's an internal coordinate Hessian, or for other potentials
! that lack analytic second derivatives.
! ENDNUMHESS2 uses a more accurate Richardson extrapolation, which takes twice
! as long. The displacement, DELTA, can be specified.
! 
         ELSE IF (WORD.EQ.'ENDNUMHESS') THEN
            ENDNUMHESS=.TRUE.
            ENDHESS=.TRUE.
         ELSE IF (WORD.EQ.'ENDNUMHESS2') THEN
            ENDNUMHESS=.TRUE.
            ENDNUMHESS2=.TRUE.
            ENDHESS=.TRUE.
            IF (NITEMS.GT.1) CALL READF(ENDNUMHESSDELTA)

! kr366> keyword to dump frqs of minima
         ELSE IF (WORD.EQ. 'DUMPFRQS') THEN
            IF (ENDNUMHESS .OR. ENDHESS) THEN
                DUMPFRQST=.TRUE.
            ELSE
                WRITE(*,'(A)') ' keywords> for DUMPFRQS ENDHESS or ENDNUMHESS need to be set'
            ENDIF 
! 
! 
! 
         ELSE IF (WORD == 'ERROREXIT') THEN
            PRINT *, 'ERROR EXIT'
            CALL FLUSH(6)
            STOP
! 
! Cutoff below which Hessian eigenvalues are considered to be zero.
! 
         ELSE IF (WORD.EQ.'EVCUT') THEN
            CALL READF(EVCUT)
! 
! Specify evolving strings.
! 
         ELSE IF (WORD.EQ.'EVOLVESTRING') THEN
            EVOLVESTRINGT = .TRUE.
! 
! sf344> extra repulsive LJ site for PY ellipsoids
! 
         ELSE IF (WORD.EQ.'EXTRALJSITE') THEN
            LJSITE=.TRUE.
            CALL READF(PEPSILON1(1))
            CALL READF(PSCALEFAC1(1))
            MAXINTERACTIONS=1
            IF(NITEMS.GT.3) THEN
               CALL READF(PSCALEFAC2(1))
               WRITE(*,'(A,3F8.3)') ' keyword> primary and secondary apex sites will be used, epsilon and heights: ',
     &         PEPSILON1(1), PSCALEFAC1(1), PSCALEFAC2(1)
               IF(.NOT.LJSITEATTR) THEN
                  MAXINTERACTIONS=3
               ELSE
                  MAXINTERACTIONS=4
               END IF
            ELSE
               WRITE(*,'(A,2F8.3)') ' keyword> primary apex sites will be used, epsilon and height: ', PEPSILON1(1), PSCALEFAC1(1)
            END IF
            IF(NITEMS.GT.4) THEN           ! binary ellipsoidal clusters will be set up only for two apex sites, not one
               BLJSITE=.TRUE.               ! we also won't use the sigma parameter from now on, epsilon is enough for repulsive sites
               CALL READF(PEPSILON1(2))
               CALL READF(PSCALEFAC1(2))
               CALL READF(PSCALEFAC2(2))
               CALL READF(PEPSILON1(3))     ! this is epsilon for the interaction between A and B type ellipsoids
               MAXINTERACTIONS=3 ! attractive secondary apex sites not incorporated for binary systems
               WRITE(*,'(A,3F8.3)') ' keyword> binary system with primary and secondary apex sites, ' //
     &         'epsilon and heights for 2nd type particle: ', PEPSILON1(2), PSCALEFAC1(2), PSCALEFAC2(2)
            END IF
         ELSE IF (WORD.EQ.'EXTRALJSITEATTR') THEN
            LJSITE=.TRUE.
            LJSITEATTR=.TRUE.
            CALL READF(PSIGMAATTR(1))
            CALL READF(PEPSILONATTR(1))
            CALL READF(PSIGMAATTR(2))
            CALL READF(PEPSILONATTR(2))
            WRITE(*,'(A,4F8.3)') 'keyword> primary and secondary apex sites '//
     &      'with normal LJ attraction, sigmas and epsilons: ',
     &      PSIGMAATTR(1), PEPSILONATTR(1), PSIGMAATTR(2), PEPSILONATTR(2)
            MAXINTERACTIONS=4
         ELSE IF (WORD.EQ.'LJSITECOORDS') THEN
            LJSITECOORDST=.TRUE.
            CALL READF(LJSITECOORDS(1))
            CALL READF(LJSITECOORDS(2))
            CALL READF(LJSITECOORDS(3))
         ELSE IF (WORD.EQ.'PYBINARY') THEN
            PYBINARYT=.TRUE.
            ANGLEAXIS2=.TRUE.
            RBAAT=.TRUE.
            ELLIPSOIDT=.TRUE.
            RADIFT=.TRUE.
            NRBSITES = 1
            ALLOCATE(RBSITE(NRBSITES,3))
            CALL READI(PYBINARYTYPE1)
            CALL READF(PYA11(1))
            CALL READF(PYA11(2))
            CALL READF(PYA11(3))
            CALL READF(PYA21(1))
            CALL READF(PYA21(2))
            CALL READF(PYA21(3))
            CALL READF(PYA12(1))
            CALL READF(PYA12(2))
            CALL READF(PYA12(3))
            CALL READF(PYA22(1))
            CALL READF(PYA22(2))
            CALL READF(PYA22(3))
            CALL READF(PYSIGNOT)
            CALL READF(PYEPSNOT)
            IF(NITEMS.GT.16) THEN
               CALL READF(PCUTOFF)
               PARAMONOVCUTOFF=.TRUE.
               PCUTOFF=PCUTOFF*PYSIGNOT
               write (*,*) "PY Potential. PCutoff ON:",PCUTOFF
            END IF
            IF(.NOT.ALLOCATED(PYA1bin)) ALLOCATE(PYA1bin(NATOMS/2,3))
            IF(.NOT.ALLOCATED(PYA2bin)) ALLOCATE(PYA2bin(NATOMS/2,3))
            DO J1=1,NATOMS/2
               IF(J1<=PYBINARYTYPE1) THEN
                  PYA1bin(J1,:)=PYA11(:)
                  PYA2bin(J1,:)=PYA21(:)
               ELSE
                  PYA1bin(J1,:)=PYA12(:)
                  PYA2bin(J1,:)=PYA22(:)
               END IF
            END DO
! 
! Obsolete. Allows for extra steps in LBFGS minimisations for CHARMM.
! 
         ELSE IF (WORD.EQ.'EXTRASTEPS') THEN
            CALL READF(EXTRASTEPS)
! 
! FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
! 
! 
! Distance dependent dielectric for Paul Mortenson`s amber
! Obsolete keyword now commented out to avoid accidental use
! 
! ELSE IF (WORD.EQ.'FAKEWATER') THEN
! FAKEWATER=.TRUE.
! WRITE (*,'(A)') ' SETTINGS Distance dependent dielectric will be used'
! 
! Integer variable to distinguish output files from parallel maiden jobs
! 
         ELSE IF (WORD.EQ.'FILTH') THEN
            IF (FILTH.EQ.0) THEN
               CALL READI(FILTH)
            ELSE
               WRITE(*,'(A)') 'WARNING **** FILTH keyword in odata was overridden by command line argument'
            ENDIF
! 
! Specifies that FIXIMAGE should be set permanently after step
! FIXAFTER. This effectively freezes the interacting images in different supercells
! for calculations with periodic boundary conditions.
! 
         ELSE IF (WORD.EQ.'FIXAFTER') THEN
            CALL READI(FIXAFTER)
! 
! Strings keyword.
! 
         ELSE IF (WORD.EQ.'FIXATMS') THEN
            FIXATMS = .TRUE.
! 
! Fix uphill direction until force changes sign.
! T12FAC is the fraction of the first collision time to be used in HSMOVE
! 
         ELSE IF (WORD.EQ.'FIXD') THEN
            FIXD=.TRUE.
            IF (NITEMS.GT.1) CALL READF(T12FAC)
            NMOVE=1
            IF (NITEMS.GT.2) CALL READF(DTHRESH)
!
! FLATEDIFF: maximum energy difference when judging whether a dneb path is flat
!
         ELSE IF (WORD.EQ.'FLATTEST') THEN
            FLATTESTT=.TRUE.
            IF (NITEMS.GT.1) CALL READF(FLATEDIFF)

! 
! FRACTIONAL: constant pressure calculation using fractional coordinates
! 
         ELSE IF (WORD.EQ.'FRACTIONAL') THEN
            FRACTIONAL=.TRUE.

! 
! Frozen atoms.
! 
         ELSE IF (WORD.EQ.'FREEZE') THEN
            FREEZE=.TRUE.

! 
! csw34> If atoms specified on the FREEZE line
! 
            IF (NITEMS.GT.1) THEN
               DO J1=1,NITEMS-1
                  NFREEZE=NFREEZE+1
                  CALL READI(NDUM)
                  FROZEN(NDUM)=.TRUE.
               ENDDO
            ELSE

               ! 
               ! csw34> Otherwise, check to see if atoms specified in 'frozen' file
               ! 
               INQUIRE(FILE='frozen',EXIST=YESNO)
               IF (YESNO) THEN
                  FUNIT=GETUNIT()
                  OPEN(FUNIT,FILE='frozen',STATUS='OLD')
                  READ(FUNIT,*) NFREEZE
                  DO J1=1,NFREEZE
                     READ(FUNIT,*) NDUM
                     FROZEN(NDUM)=.TRUE.
                  ENDDO
                  ! 
                  ! csw34> If neither, FREEZE used incorrectly - STOP
                  ! 
               ELSE
                  WRITE (*,'(A)') ' ERROR: FREEZE specified incorrectly'
                  WRITE(*,*) "Specify frozen atoms either on the keyword line, or in a file called 'frozen'"
                  STOP
               ENDIF
            ENDIF

! 
Cjbr36      > FREEZERANGE of atoms
! 
         ELSE IF (WORD.EQ.'FREEZERANGE') THEN
            FREEZE=.TRUE.
            FREEZERANGE=.TRUE.

            IF (NITEMS.GT.1 .and. NITEMS .LT. 4) THEN
               CALL READI(NDUM)
               J1=NDUM
               CALL READI(NDUM)
               J2=NDUM
               DO J3=J1,J2
                  NFREEZE=NFREEZE+1
                  FROZEN(J3)=.TRUE.
               ENDDO
            ELSE
               WRITE (*,'(A)') ' ERROR: FREEZERANGE specified incorrectly'
            ENDIF

! csw34>
! Frozen residues (to be converted to frozen atoms)
! 
         ELSE IF (WORD.EQ.'FREEZERES') THEN
            FREEZE=.TRUE.
            FREEZERES=.TRUE.
! The FROZENRES array is then filled with the residue number from the
! data file
            DO J1=1,NITEMS-1
               CALL READI(NDUM)
               FROZENRES(NDUM)=.TRUE.
            ENDDO
! Finally, the frozen residue numbers are converted into frozen atom
! numbers. This is also forcefield dependant and must be done when we
! know which forcefield to use (i.e. in the CHARMM block above)
! 
! csw34> FREEZEGROUP centreatom radius
! FREEZEs all atoms within radius angstroms of centreatom (labelled by index)
! 
         ELSE IF (WORD.EQ.'FREEZEGROUP') THEN
            FREEZE=.TRUE.
            FREEZEGROUPT=.TRUE.
            CALL READI(GROUPCENTRE)
            CALL READF(GROUPRADIUS)
            IF(NITEMS.GT.3) CALL READA(FREEZEGROUPTYPE)

         ELSE IF (WORD.EQ.'DUMPMODE') THEN
            IF(DUMPV.AND.MWVECTORS) THEN
               DO J1=1,NITEMS-1
                  CALL READI(NDUM)
                  DUMPMODEN(NDUM)=.TRUE.
               ENDDO
            ENDIF


            IF ((PERMDIST.OR.LOCALPERMDIST.OR.LPERMDIST).OR.PERMDISTINIT) THEN
               NDUMMY=0
               DO J1=1,NPERMGROUP
                  DO J2=1,NPERMSIZE(J1)
                     IF (FROZEN(PERMGROUP(NDUMMY+J2))) THEN
                        PRINT '(A,I8,A)',' keyword> ERROR atom ',PERMGROUP(NDUMMY+J2),' cannot be frozen and permuted'
                        STOP
                     ENDIF
                  ENDDO
                  NDUMMY=NDUMMY+NPERMSIZE(J1)
               ENDDO
            ENDIF
! 
! Strings keyword.
! 
         ELSE IF (WORD.EQ.'FREEZENODES') THEN
            FREEZENODEST=.TRUE.
            CALL READF(FREEZETOL)

         ELSE IF (WORD.EQ.'FRQCONV') THEN
            ! sn402: I'm implementing a new way of handling unit conversions when calculating frequencies.
            ! In line with the documention on the OPTIM website, the default behaviour is to leave all frequencies
            ! in internal units and only convert to SI at the very end of a calculation (typically after obtaining rate
            ! constants). However, this will not be appropriate for some potentials, so I'm implementing a general unit
            ! conversion scheme. All unit conversions will eventually be done using a variable FRQCONV set in this
            ! subroutine, which is the factor required to convert a frequency from the internal units of the current
            ! potential into another desired unit system.
            ! For example, AMBER and CHARMM work with energy units of kCalmol^-1 and distance units of Angstrom, so the
            ! natural frequency units are (kCal mol^-1/(amu Angstrom^2))^1/2
            ! However, we usually want frequencies to be written out in radians/s, for use by PATHSAMPLE. In that case,
            ! FRQCONV = sqrt(4.184E26) = 2.045483D13
            ! The AMBER, CHARMM and NAB keywords will use this as their default value of FRQCONV.
            ! If you want your potential to default to printing frequencies in a unit other than internal units, set a value
            ! of FRQCONV in the relevant block of this subroutine.
            ! To override the default frequency units for a particular job, use this keyword, FRQCONV conv
            ! The variable FRQCONV will be set equal to this argument conv at the end of the subroutine, so it doesn't matter
            ! whether this keyword is placed before or after the keyword which activates your potential.
            ! Note that if you want to override a default, your conversion factor will be applied to the raw eigenvalues instead
            ! of the default, not as well as the default. So if you want AMBER frequencies in cm^-1, you must specify the conversion
            ! factor from internals to cm^-1, not the conversion factor for s^-1 to cm^-1. (In this example, the conversion factor
            ! required is FRQCONV 108.52D0)
            CALL READF(DUMMY_FRQCONV)

! 
! GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG
! 
! 
! GAMESS-UK tells the program to read derivative information in
! GAMESS-UK format.                                        - default FALSE
         ELSE IF (WORD.EQ.'GAMESS-UK') THEN
            GAMESSUK=.TRUE.
            CALL READA(SYS)
            DO J1=1,80
               IF (SYS(J1:J1).EQ.' ') THEN
                  LSYS=J1-1
                  GOTO 112
               ENDIF
            ENDDO
112         IF (NITEMS.GT.2) THEN
               CALL READA(EDITIT)
            ELSE
               EDITIT='editit.' // SYS(1:LSYS)
            ENDIF
! 
! GAMESS-US tells the program to read derivative information in
! GAMESS-US format.                                        - default FALSE
         ELSE IF (WORD.EQ.'GAMESS-US') THEN
            GAMESSUS=.TRUE.
            CALL READA(SYS)
            DO J1=1,80
               IF (SYS(J1:J1).EQ.' ') THEN
                  LSYS=J1-1
                  GOTO 111
               ENDIF
            ENDDO
111         IF (NITEMS.GT.2) THEN
               CALL READA(EDITIT)
            ELSE
               EDITIT='editit.' // SYS(1:LSYS)
            ENDIF
! 
! GAUSSIAN tells the program to read derivative information in
! Gaussian92 format.                                  - default FALSE
         ELSE IF (WORD.EQ.'GAUSSIAN') THEN
            GAUSSIAN=.TRUE.
! Gaussian03 interface            
         ElSE IF (WORD.EQ.'GAUSSIAN03') THEN
            GAUSSIAN03=.TRUE.
            CALL READI(GAUSSIANCHARGE)
            CALL READI(GAUSSIANMULTI)
            OPEN(LUNIT,FILE='gaussian.xyz',STATUS='REPLACE')
            CLOSE(LUNIT)
! Gaussian09 interface            
         ElSE IF (WORD.EQ.'GAUSSIAN09') THEN
            GAUSSIAN09=.TRUE.
            CALL READI(GAUSSIANCHARGE)
            CALL READI(GAUSSIANMULTI)
            OPEN(LUNIT,FILE='gaussian.xyz',STATUS='REPLACE')
            CLOSE(LUNIT)
! Gaussian16 interface            
         ElSE IF (WORD.EQ.'GAUSSIAN16') THEN
            GAUSSIAN16=.TRUE.
            CALL READI(GAUSSIANCHARGE)
            CALL READI(GAUSSIANMULTI)
            OPEN(LUNIT,FILE='gaussian.xyz',STATUS='REPLACE')
            CLOSE(LUNIT)
! DC430 >

         ELSE IF (WORD .EQ. 'GB') THEN
            GBT   = .TRUE.
            RBAAT = .TRUE.
            CALL READF(GBKAPPA)
            CALL READF(GBKAPPRM)
            CALL READF(GBMU)
            CALL READF(GBNU)
            CALL READF(GBSIGNOT)
            CALL READF(GBEPSNOT)

            GBCHI    = (GBKAPPA ** 2 - 1.D0) / (GBKAPPA ** 2 + 1.D0)
            GBCHIPRM = (GBKAPPRM**(1.D0/GBMU)-1.D0) / (GBKAPPRM**(1.D0/GBMU)+1.D0)

         ELSE IF (WORD .EQ. 'GBD') THEN

            GBDT  = .TRUE.
            RBAAT = .TRUE.
            CALL READF(GBKAPPA)
            CALL READF(GBKAPPRM)
            CALL READF(GBMU)
            CALL READF(GBNU)
            CALL READF(GBSIGNOT)
            CALL READF(GBEPSNOT)

            GBCHI    = (GBKAPPA ** 2 - 1.D0) / (GBKAPPA ** 2 + 1.D0)
            GBCHIPRM = (GBKAPPRM**(1.D0/GBMU)-1.D0) / (GBKAPPRM**(1.D0/GBMU)+1.D0)

! -----------------------------
! GDIIS x y z  x=cutoff on previous RMS force below which GDIIS
! may be applied, y=NDIIA the dimension of the DIIS
! problem to solve, z=NINTV the interval between
! GDIIS steps                                     - default OFF
         ELSE IF (WORD .EQ. 'GDIIS') THEN
            PRINT '(A)','keyword> GDIIS keyword not available'
            STOP
! IF (NITEMS.LT.4) THEN
! DTEST=.FALSE.
! PRINT*,'Error in GDIIS input - insufficient items'
! ELSE
! DTEST=.TRUE.
! CALL READF(PCUT)
! CALL READI(NDIIA)
! CALL READI(NINTV)
! ENDIF
! IF (NDIIA.GT.NDIIS) THEN
! WRITE(*,'(A,I6)') ' NDIIA too large=',NDIIA
! STOP
! ENDIF
! 
         ELSE IF (WORD.EQ.'GDSQ') THEN
            GDSQT=.TRUE.
            CALL READF(GDSQ)
! GEOMDIFFTOL specifies the maximum displacement between identical permutational isomers in connect.
! 
         ELSE IF (WORD.EQ.'GEOMDIFFTOL') THEN
            CALL READF(GEOMDIFFTOL)
! 
! General LJ for mixed systems
! 
         ELSE IF (WORD.EQ.'GLJ') THEN
            GLJT=.TRUE.
            NGLJ=NITEMS-1 ! the number of different species, since we specify how many of each on this line
            IF (ALLOCATED(NGLJSPECIES)) DEALLOCATE(NGLJSPECIES)
            ALLOCATE(NGLJSPECIES(NGLJ))
            DO J1=1,NITEMS-1
               CALL READI(NGLJSPECIES(J1)) ! number of nodes
            ENDDO
            IF (ALLOCATED(GLJEPS)) DEALLOCATE(GLJEPS)
            IF (ALLOCATED(GLJSIG)) DEALLOCATE(GLJSIG)
            ALLOCATE(GLJEPS(NGLJ,NGLJ),GLJSIG(NGLJ,NGLJ))
            WRITE(*,'(A,I8,A)') 'keyword> General LJ potential with ',NGLJ,' particle types. Numbers of each paticle:'
            DO J1=1,NGLJ
               WRITE(*,'(I6)',ADVANCE='NO') NGLJSPECIES(J1)
            ENDDO
            WRITE(*,'(A)') ' '
            WRITE(*,'(A)') 'keyword> Epsilon parameters:'
            DO J1=1,NGLJ
               READ(5,*) (GLJEPS(J2,J1),J2=1,J1)
               DO J2=1,J1
                  GLJEPS(J1,J2)=GLJEPS(J2,J1)
                  WRITE(*,'(F15.4)',ADVANCE='NO') GLJEPS(J2,J1)
               ENDDO
               WRITE(*,'(A)') ' '
            ENDDO
            WRITE(*,'(A)') ' '
            WRITE(*,'(A)') 'keyword> Sigma parameters:'
            DO J1=1,NGLJ
               READ(5,*) (GLJSIG(J2,J1),J2=1,J1)
               DO J2=1,J1
                  GLJSIG(J1,J2)=GLJSIG(J2,J1)
                  WRITE(*,'(F15.4)',ADVANCE='NO') GLJSIG(J2,J1)
               ENDDO
               WRITE(*,'(A)') ' '
            ENDDO
! 
! Paul Whitford Structure-based SMOG model
! 
         ELSE IF (WORD.EQ.'SBM') THEN
            GOT=.TRUE.
! 
! GRADIENT n prints the gradients along the Hessian eigendirections
! every n cycles                                    - default OFF
! 
         ELSE IF (WORD .EQ. 'GRAD4') THEN ! 4-point gradient in finite_differences.f90
            GRAD4T=.TRUE.
            print *, 'use 4-point gradient'
         ELSE IF (WORD .EQ. 'GRADIENTS') THEN
            PGRAD=.TRUE.
            CALL READI(NGRADIENTS)
! 
! GRADSQ specifies optimisation of the modulus gradient. This is a really bad idea!
! 
         ELSE IF (WORD.EQ.'GRADSQ') THEN
            GRADSQ=.TRUE.
            IF (NITEMS.GT.1) CALL READF(GSTHRESH)
            IF (NITEMS.GT.2) CALL READI(NSPECIAL)
            IF (NITEMS.GT.3) CALL READI(NALLOW)
! 
! Approximation to use for the gradient in DNEB routine NEB/grad.f90 - default is "dneb"
! 
         ELSE IF (WORD == 'GRADTYPE') THEN
            CALL READA(GRADTYPE)
! 
! Attempt to interpolate between endpoints using a great circle. Not a huge success.
! 
         ELSE IF (WORD.EQ.'GREATCIRCLE') THEN
            GREATCIRCLET=.TRUE.
            CALL READF(GCMXSTP)
            CALL READI(GCIMAGE)
            CALL READI(GCSTEPS)
            CALL READF(GCCONV)
! 
! EF: growing strings
! number of images and iterations for first iteration;
! reparametrization tolerance, growth tolerance, convergence tolerance
! maximum LBFGS step; LBFGS memory
! 
         ELSE IF (WORD.EQ.'GROWSTRING') THEN
            GROWSTRINGT = .TRUE.
            FCD = .TRUE.
            IF (NITEMS.GT.1) CALL READI(nnNIMAGE)
            IF (NITEMS.GT.2) CALL READF(GSITERDENSITY)
            IF (NITEMS.GT.3) CALL READF(REPARAMTOL)
            IF (NITEMS.GT.4) CALL READF(GSGROWTOL)
            IF (NITEMS.GT.5)CALL READF(GSCONV)
            IF (NITEMS.GT.6) CALL READF(GSMXSTP)
! 
! Set the maximum total iteration density for the
! growing string method. This specifies the maximum evolution iterations allowed per
! total image number, including the iterations while the string is still
! growing. If {\it itd\/} is less than 0, then this parameter is turned off
! and there is no limit on the total iterations (this is the default).
! 
         ELSE IF (WORD.EQ.'GSMAXTOTITD') THEN
            CALL READF(GSMAXTOTITD)

! hk286 > Generalised Thomson    !
         ELSE IF (WORD .EQ. 'GTHOMSONBIN') THEN
         GTHOMSONT = .TRUE.
         THOMSONT=.TRUE.
         RIGIDBODY = .TRUE.
         GTHOMPOT=7
         CALL READI(GTHOMSONBINN)
         NTYPEA = GTHOMSONBINN
         CALL READF(GTHOMSONBIND)
         ELSE IF (WORD .EQ. 'GTHOMSON') THEN
            GTHOMSONT = .TRUE.
! THOMSONT=.TRUE.
! RIGID = .TRUE.
            CALL READI(GTHOMMET)
            CALL READF(GTHOMSONZ)
            IF (NITEMS.GT.3) THEN
               CALL READF(GTHOMSONC)
            ENDIF
            IF (NITEMS.GT.4) THEN
               CALL READF(GTHOMSONC2)
               IF ( (GTHOMMET .EQ. 2) .AND. (GTHOMSONC2 > 0.0D0) ) THEN
                  GTHOMSONZ = LOG((GTHOMSONC2/GTHOMSONC)+ SQRT((GTHOMSONC2/GTHOMSONC)**2-1.0D0))*GTHOMSONC
                  ! PRINT *, GTHOMSONZ
               ENDIF
            ENDIF
            IF (NITEMS.GT.5) CALL READF(VUNDULOID)

            IF ((GTHOMMET .EQ. 3) .OR. (GTHOMMET .EQ. 4)) THEN
               CALL CONVERTUNDULOIDPARAMETERS(VUNDULOID/2.0D0)
            ENDIF
            IF ((GTHOMMET .EQ. 3) .OR. (GTHOMMET .EQ. 4)) THEN
               CALL FINDNGZ()
            ENDIF
            CALL INIGTHOMSON()

         ELSE IF (WORD .EQ. 'GTHOMSONPOT') THEN
            CALL READI(GTHOMPOT)
            IF (NITEMS.GT.2) THEN
               CALL READF(GThomsonSigma)
            ENDIF
            IF (NITEMS.GT.2) THEN
               CALL READF(GThomsonRho)
            ENDIF

! 
! Try to guess an interpolated path for sequential coordinate changes.
! GSTEPS is the number of step to be tried for each sequential coordinate.
! MAXGCYCLES is the number of sweeps through pairwise exchanges
! GTHRESHOLD is the coordinate change above which sequential changes are considered.
! MAXINTE is the convergence criterion for the maximum allowed interpolated energy.
! 
         ELSE IF (WORD.EQ.'GUESSPATH') THEN
            GUESSPATHT=.TRUE.
            CALL READI(GSTEPS)
            CALL READI(MAXGCYCLES)
            CALL READF(GTHRESHOLD)
! CALL READF(MAXINTE)
! 
! Use dihedral twisting in place of DNEB for
! transition state guesses with CONNECT for CHARMM and UNRES.
! 
         ELSE IF (WORD.EQ.'GUESSTS') THEN
            GUESSTST=.TRUE.
            IF (NITEMS.GT.1) CALL READF(GUESSTHRESH)
         ELSE IF (WORD.EQ.'GUPTA') THEN
            CALL READI(GUPTATYPE)
         ELSE IF (WORD.EQ.'BGUPTAT') THEN
            CALL READI(NTYPEA)
            CALL READF(AAA)
            CALL READF(PAA)
            CALL READF(QAA)
            CALL READF(ZAA)
            CALL READF(R0AA)
         ELSEIF (WORD.EQ.'BGUPTATAB')THEN
            CALL READF(AAB)
            CALL READF(PAB)
            CALL READF(QAB)
            CALL READF(ZAB)
            CALL READF(R0AB)
         ELSEIF (WORD.EQ.'BGUPTATBB')THEN
            CALL READF(ABB)
            CALL READF(PBB)
            CALL READF(QBB)
            CALL READF(ZBB)
            CALL READF(R0BB)

! 
! HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH
! 
! Dump Hessians for stationary points - similar to path.info calls.
! 
         ELSE IF (WORD.EQ.'HESSDUMP') THEN
            HESSDUMPT=.TRUE.
! 
! For the growing string or evolving string double-ended
! transition state search method, use the method described in the appendix of
! Peters et al\cite{PetersHBC04} to calculate the Newton-Raphston search
! direction. Namely, the Hessian is approximated based on changes in the
! gradient, and the tangential component of $-\mathbf{Hf^\perp}$ is projected
! out. By default, the Hessian used is actually an approximation to the
! derivative matrix of $f^\perp$ rather than the gradient.
! 
         ELSE IF (WORD.EQ.'HESSGRAD') THEN
            HESSGRAD = .TRUE.
! 
! HIGHESTIMAGE - only use the highest non-endpoint image in a double-ended
! MECCANO-type run.
! 
         ELSE IF (WORD.EQ.'HIGHESTIMAGE') THEN
            HIGHESTIMAGE=.TRUE.
! 
! HUPDATE specifies that a Hessian updating procedure should be used.
! 
         ELSE IF (WORD .EQ. 'HUPDATE') THEN
            HUPDATE=.TRUE.
            NHUP=0
            IF (NITEMS.GT.1) THEN
               CALL READI(NSTHUP)
            ENDIF
            IF (NITEMS.GT.2) THEN
               CALL READI(INTHUP)
            ENDIF
            IF (NITEMS.GT.3) THEN
               CALL READF(PHIG)
            ENDIF
! 
! Hybrid BFGS/eigenvector-following minimisation.
! 
         ELSE IF (WORD.EQ.'HYBRIDMIN') THEN
            HYBRIDMINT=.TRUE.
            CALL READI(HMNEVS)      ! maximum steps to converge smallest eigenvalue in Rayleigh-Ritz
            CALL READI(HMNBFGSMAX1) ! maximum tangent space LBFGS steps if eigenvalue unconverged
            CALL READI(HMNBFGSMAX2) ! maximum tangent space LBFGS steps if eigenvalue converged
            CALL READF(HMCEIG)      ! convegence criterion for eigenvalue
            CALL READF(HMMXSTP)     ! maximum step size for EF steps
            CALL READI(HMNSTEPS)    ! maximum number of hybrid minimisation steps
            CALL READF(HMEVMAX)     ! If the lowest eigenvalue goes above HMEVMAX then exit
            CALL READA(HMMETHOD)    ! Choose between EF and Page-McIver steepest-descent steps
            IF (NITEMS.GT.9) THEN
               CALL READI(HMNEVL)   ! maximum steps for iterative calculation of largest eigenvalue if applicable
            ENDIF
! 
! IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
! 
! 
! Add an icosahedral field to the potential of magnitude FIH.
! 
         ELSE IF (WORD.EQ.'IH') THEN
            FIELDT=.TRUE.
            IHT=.TRUE.
            CALL READF(FIH)
! 
! Search for a saddle of index INDEX if
! SEARCH 2 is specified. See also KEEPINDEX. Also works with BFGSTS
! up to a maximum of index 50, but NOIT must be set and a Hessian is needed.
! 
         ELSE IF (WORD.EQ.'INDEX') THEN
            CALL READI(HINDEX)
!
! Radial shift to make space for new atoms.
!
         ELSE IF (WORD.EQ.'QCIRADSHIFT') THEN
            QCIRADSHIFTT=.TRUE.
            IF (NITEMS.GT.1) CALL READF(QCIRADSHIFT)
            WRITE(*,'(A,G20.10)') ' keyword> Shifting unconstrained atoms away from added atoms by ',QCIRADSHIFT
!
! Use new QCILPERMDIST check instead of this flag
!
!        ELSE IF (WORD.EQ.'QCIPERMCHECK') THEN
!           QCIPERMCHECK=.TRUE.
!           CALL READI(QCIPERMCHECKINT)
! 
! Images for INTCONSTRAINT
! 
         ELSE IF ((WORD.EQ.'INTIMAGE').OR.(WORD.EQ.'QCIIMAGE')) THEN
            IF (NITEMS.GT.1) CALL READF(IMSEPMIN)
            IF (NITEMS.GT.2) CALL READF(IMSEPMAX)
            IF (NITEMS.GT.3) CALL READI(INTIMAGE)
            IF (NITEMS.GT.4) CALL READI(MAXINTIMAGE)
            IF (NITEMS.GT.5) CALL READI(INTNTRIESMAX)
            IF (NITEMS.GT.6) CALL READI(INTIMAGEINCR)
            IF (NITEMS.GT.7) CALL READI(INTIMAGECHECK)
!
! Maximum distance for constrained atoms
!
         ELSE IF ((WORD.EQ.'INTCONCUT').OR.(WORD.EQ.'QCICONCUT')) THEN
            CALL READF(INTCONCUT)
!
! QCI LOPERMDIST checks for images
!
         ELSE IF ((WORD.EQ.'INTLPERMDIST').OR.(WORD.EQ.'QCILPERMDIST')) THEN
            QCILPERMDIST=.TRUE.
            IF (NITEMS.GT.1) CALL READI(QCIPDINT)   ! interval for check
            IF (NITEMS.GT.2) CALL READF(QCIPERMCUT) ! maximum allowed distance for valid alignment
            WRITE(*,'(A,I8,A,G20.10)') ' keywords> QCI active permutational alignment check every ',QCIPDINT,' steps, tolerance=', 
     &                                  QCIPERMCUT     
! 
! Use constraint potential for initial interpolation in each cycle.
! 
         ELSE IF ((WORD.EQ.'INTCONSTRAINT').OR.(WORD.EQ.'QCI')) THEN
            INTCONSTRAINTT=.TRUE.
            IF (NITEMS.GT.1) CALL READF(INTCONSTRAINTTOL)
            IF (NITEMS.GT.2) CALL READF(INTCONSTRAINTDEL)
            IF (NITEMS.GT.3) CALL READF(INTCONSTRAINTREP)
            IF (NITEMS.GT.4) CALL READF(INTCONSTRAINREPCUT)
            IF (NITEMS.GT.5) CALL READF(INTCONFRAC)
            IF (NITEMS.GT.6) CALL READI(INTCONSEP)
            IF (NITEMS.GT.7) CALL READI(INTREPSEP)
            IF (NITEMS.GT.8) CALL READI(INTSTEPS1)
            IF (NITEMS.GT.9) CALL READI(INTCONSTEPS)
            IF (NITEMS.GT.10) CALL READI(INTRELSTEPS)
            IF (NITEMS.GT.11) CALL READF(MAXCONE)
            IF (NITEMS.GT.12) CALL READF(INTRMSTOL)
! 
! Parse the congeom file to obtain number of reference minima for setting up
! constraints and get their geometries.
! 
            INQUIRE(FILE='congeom.dat',EXIST=YESNO)
            IF (YESNO) THEN
               CONDATT=.TRUE.
               LUNIT=GETUNIT()
               OPEN(LUNIT,FILE='congeom.dat',STATUS='OLD')
               READ(LUNIT,*) NCONGEOM
               ! 
               ! we only need the first reference if we have congeom.dat
               ! However, CONGEOM is passed to make_conpot and declared with first dimension NCPFIT there!
               ! 
               ALLOCATE(CONGEOM(NCONGEOM,3*NATOMS))
               CONGEOM(2:NCONGEOM,1:3*NATOMS)=0.0D0
               READ(LUNIT,*) CONGEOM(1,1:3*NATOMS)
               READ(LUNIT,*) NCONSTRAINTFIX
               ALLOCATE(CONIFIX(NCONSTRAINTFIX),CONJFIX(NCONSTRAINTFIX),
     &         CONDISTREFFIX(NCONSTRAINTFIX),CONCUTFIX(NCONSTRAINTFIX))
               READ(LUNIT,*) CONIFIX(1:NCONSTRAINTFIX)
               READ(LUNIT,*) CONJFIX(1:NCONSTRAINTFIX)
               READ(LUNIT,*) CONDISTREFFIX(1:NCONSTRAINTFIX)
               READ(LUNIT,*) CONCUTFIX(1:NCONSTRAINTFIX)
               READ(LUNIT,*) NREPULSIVEFIX
               ALLOCATE(REPIFIX(NREPULSIVEFIX),REPJFIX(NREPULSIVEFIX),REPCUTFIX(NREPULSIVEFIX))
               READ(LUNIT,*) REPIFIX(1:NREPULSIVEFIX)
               READ(LUNIT,*) REPJFIX(1:NREPULSIVEFIX)
               READ(LUNIT,*) REPCUTFIX(1:NREPULSIVEFIX)
               CLOSE(LUNIT)
               PRINT '(A)',' keyword> Constraint potential parameters read from file congeom.dat'
               INTCONMAX=NCONSTRAINTFIX
               NREPMAX=NREPULSIVEFIX
               ALLOCATE(CONI(INTCONMAX),CONJ(INTCONMAX),CONDISTREF(INTCONMAX),CONCUT(INTCONMAX),CONOFFLIST(INTCONMAX),
     &                  CONOFFTRIED(INTCONMAX))
               CONOFFTRIED(1:INTCONMAX)=.FALSE.
               ALLOCATE(REPI(NREPMAX),REPJ(NREPMAX),NREPI(NREPMAX),NREPJ(NREPMAX),REPCUT(NREPMAX),NREPCUT(NREPMAX))
               ALLOCATE(CONACTIVE(NCONSTRAINTFIX))
            ELSE
               INQUIRE(FILE='congeom',EXIST=CONFILE)
               NCONGEOM=0
               IF (.NOT.CONFILE) THEN
                  PRINT '(A)',' keyword> WARNING *** no congeom file found. Will use end point minima only.'
               ELSE
                  LUNIT=GETUNIT()
                  OPEN(LUNIT,FILE='congeom',STATUS='OLD')
                  DO
                     READ(LUNIT,*,END=861) DUMMY1(1)
                     NCONGEOM=NCONGEOM+1
                  ENDDO
861               CONTINUE
                  NCONGEOM=NCONGEOM/NATOMS
                  ALLOCATE(CONGEOM(NCONGEOM,3*NATOMS))
                  REWIND(LUNIT)
                  DO J1=1,NCONGEOM
                     READ(LUNIT,*) CONGEOM(J1,1:3*NATOMS)
                  ENDDO
                  CLOSE(LUNIT)
                  PRINT '(A,I6,A)',' keyword> Read ',NCONGEOM,' reference geometries from congeom file'
                  IF (NCONGEOM.LT.2) PRINT '(A)',' WARNING *** insufficient reference geometries - using end point minima'
               ENDIF
            ENDIF
! 
! DO NOT Use the quasi-continuous metric for connection attempts, instead of distance.
! 
            INTERPCOSTFUNCTION=.FALSE.
!
! Use trilateration in placing next atom for QCI
!
      ELSE IF (WORD.EQ.'QCITRILAT') THEN
         QCITRILAT=.TRUE.
!
! Do the backbone first
!
      ELSE IF (WORD.EQ.'QCIDOBACK') THEN
         DOBACK=.TRUE.
!
! Do the backbone first and add all backbone atoms in each aa in one go
!
      ELSE IF (WORD.EQ.'QCIDOBACKALL') THEN
         DOBACKALL=.TRUE.
         DOBACK=.TRUE.
!
! Add complete amino acids
!
      ELSE IF (WORD.EQ.'QCIADDACID') THEN
         QCIADDACIDT=.TRUE.
!
! Reset when diagnosed stuck
!
      ELSE IF (WORD.EQ.'QCIRESET') THEN
         QCIRESET=.TRUE.
         IF (NITEMS.GT.1) CALL READI(QCIRESETINT1)
         IF (NITEMS.GT.2) CALL READI(QCIRESETINT2)
!
!Use topology information for QCI constraints for AMBER
!
      ELSE IF (WORD.EQ.'QCIAMBER') THEN
         QCIAMBERT=.TRUE.
         WRITE(*,'(A)') ' keyword> Use topology file for constraints in QCI'

         ELSE IF (WORD.EQ.'INTFREEZE') THEN
            INTFREEZET=.TRUE.
            IF (NITEMS.GT.1) CALL READF(INTFREEZETOL)
            IF (NITEMS.GT.2) CALL READI(INTFREEZEMIN)
! 
! Use interpolation potential for LJ.
! 
         ELSE IF (WORD.EQ.'INTLJ') THEN
            INTLJT=.TRUE.
            IF (NITEMS.GT.1) CALL READI(INTLJSTEPS)
            IF (NITEMS.GT.2) CALL READF(INTLJDEL)
            IF (NITEMS.GT.3) CALL READF(INTLJTOL)
            IF (NITEMS.GT.4) CALL READI(INTIMAGE)
            IF (NITEMS.GT.5) CALL READF(INTLJEPS)
! 
! Parse the congeom file to obtain number of reference minima for setting up
! constraints and get their geometries.
! 
            INQUIRE(FILE='congeom',EXIST=CONFILE)
            NCONGEOM=0
            IF (.NOT.CONFILE) THEN
               PRINT '(A)',' keyword> WARNING *** no congeom file found. Will use end point minima only.'
            ELSE
               LUNIT=GETUNIT()
               OPEN(LUNIT,FILE='congeom',STATUS='OLD')
               DO
                  READ(LUNIT,*,END=863) DUMMY1(1)
                  NCONGEOM=NCONGEOM+1
               ENDDO
863            CONTINUE
               NCONGEOM=NCONGEOM/NATOMS
               ALLOCATE(CONGEOM(NCONGEOM,3*NATOMS))
               REWIND(LUNIT)
               DO J1=1,NCONGEOM
                  READ(LUNIT,*) CONGEOM(J1,1:3*NATOMS)
               ENDDO
               CLOSE(LUNIT)
               PRINT '(A,I6,A)',' keyword> Read ',NCONGEOM,' reference geometries from congeom file'
               IF (NCONGEOM.LT.2) PRINT '(A)',' WARNING *** insufficient reference geometries - using end point minima'
            ENDIF
! 
! DO NOT Use the quasi-continuous metric for connection attempts, instead of distance.
! 
            INTERPCOSTFUNCTION=.FALSE.
! 
! Epsilon value in internal coordinate optimisation.
! 
         ELSE IF (WORD == 'INTEPSILON') THEN
            IF (NITEMS.GT.1) CALL READF(INTEPSILON)

! back transformation cutoff for interpolation in internals
         ELSE IF (WORD.EQ.'INTERPBACKTCUT') THEN
            CALL READF(INTERPBACKTCUT)

! must set INTINTERP as well to use INTERPCHOICE
! use internals or cartesian interpolation depending on which gives
! the lower max energy
         ELSE IF (WORD.EQ.'INTERPCHOICE') THEN
            INTERPCHOICE = .TRUE.

         ELSE IF (WORD.EQ.'INTINTERP') THEN
            INTINTERPT = .TRUE. ! interpolate with internals
            NATINT = .TRUE. ! if interpolating, assume natural internal coords
            DESMINT = .FALSE. ! intinterp and desmint are mutually exclusive
            IF (NITEMS.GT.1) CALL READI(NINTIM)

! when interpolating with internals, keep actual interpolation points.
! don't distribute images between them to make them equidistant in cartesians
         ELSE IF (WORD.EQ.'INTERPSIMPLE') THEN
            INTERPSIMPLE = .TRUE.
! 
! Internal coordinate minimisation - do not use.
! IMINCUT is the RMSG below which we take steps in internal coordinates
! 
! 
         ELSE IF (WORD.EQ.'INTMIN') THEN
            INTMINT=.TRUE.
            IF (NITEMS.GT.1) THEN
               CALL READF(IMINCUT)
            ENDIF

! align permutations of starting structures to match up internals
         ELSE IF (WORD.EQ.'INTMINPERM') THEN
            INTMINPERMT = .TRUE.
            IF (NITEMS.GT.1) THEN
               CALL READA(WORD2)
               IF (WORD2.EQ."GLYCART") THEN
                  GLYCART = .TRUE.
               ELSE
                  PRINT*, "keyword error intminperm"
               ENDIF
            ENDIF

         ELSE IF (WORD.EQ.'INTPARFILE') THEN
            USEPARFILE = .TRUE.
            CALL READA(INTPARFILE) ! file with internals parameters
! 
! Whether to add spring terms for only active atoms in INTCONSTRAINT
! 
         ELSE IF (WORD.EQ.'INTSPRINGACTIVE') THEN
            INTSPRINGACTIVET=.TRUE.
! 
! Use the inverted (negative) potential.
! 
         ELSE IF (WORD.EQ.'INVERTP') THEN
            INVERTPT=.TRUE.
! 
! JJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJJ
! 
! 
! KKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKKK
! 
! 
! KEEPINDEX: specifies that INDEX is set with
! the number of negative Hessian eigenvalues at the initial point.
! 
         ELSE IF (WORD.EQ.'KEEPINDEX') THEN
            KEEPINDEX=.TRUE.

! csw34> Specify kT in wavenumbers, below which a normal mode is
! determined to be thermally accessible. KTWN defaults to
! room temperature (207.11cm-1). This is used by the
! CHARMMDUMPMODES subroutine
         ELSE IF (WORD.EQ.'KTWN') THEN
            KTWNT=.TRUE.
            IF (NITEMS.GT.1) THEN
               CALL READF(KTWN)
            ENDIF
! 
! KINT: force constant for springs in INTCONSTRAINT calculations.
! Default zero.
! 
         ELSE IF (WORD.EQ.'KINT') THEN
            CALL READF(KINT)
            IF (NITEMS.GT.2) CALL READF(KINTENDS)
! 
! LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
! 
! Use Lanczos to diagonalize the Hamiltonian. Defaults for the three
! associated parameters are ACCLAN=1.0D-8 SHIFTLAN=1.0D-2 CUTLAN=-1.0D0.
! 
         ELSE IF (WORD.EQ.'LANCZOS') THEN
            LANCZOST=.TRUE.
            IF (NITEMS.GT.1) THEN
               CALL READF(ACCLAN)
            ENDIF
            IF (NITEMS.GT.2) THEN
               CALL READF(SHIFTLAN)
            ENDIF
            IF (NITEMS.GT.3) THEN
               CALL READF(CUTLAN)
            ENDIF
!
!  Addressable LJ
!
         ELSE IF (WORD.EQ.'LJADD') THEN
            LJADDT=.TRUE.
            LUNIT=GETUNIT()
            OPEN(LUNIT,FILE='epsilon',STATUS='OLD')
            IF (.NOT.ALLOCATED(LJADDEPS)) ALLOCATE(LJADDEPS(NATOMS,NATOMS))
            DO J1=1,NATOMS
               DO J2=1,NATOMS
                  READ(LUNIT,*) LJADDEPS(J2,J1)
                  IF (DEBUG) PRINT '(2I6,G20.10)',J1,J2,LJADDEPS(J2,J1)
               ENDDO
            ENDDO
            CLOSE(LUNIT)
         ELSE IF (WORD.EQ.'LJADD2') THEN
            LJADDT=.TRUE.
            LJADD2T=.TRUE.
            CALL READI(NADDTARGET)
            WRITE(*,'(A,I6)') 'keyword> Target cluster size is ',NADDTARGET
            IF (MOD(NATOMS,NADDTARGET).NE.0) THEN
               WRITE(*,'(A,I6)') 'keywords> ERROR, target cluster size is not a factor of the nummber of the atoms ',NATOMS
               STOP
            ENDIF
            LUNIT=GETUNIT()
            OPEN(LUNIT,FILE='epsilon',STATUS='OLD')
            IF (.NOT.ALLOCATED(LJADDEPS)) ALLOCATE(LJADDEPS(NADDTARGET,NADDTARGET))
            DO J1=1,NADDTARGET
               DO J2=1,NADDTARGET
                  READ(LUNIT,*) LJADDEPS(J2,J1)
                  WRITE(*,'(2I6,G20.10)') J1,J2,LJADDEPS(J2,J1)
               ENDDO
            ENDDO
            CLOSE(LUNIT)
         ELSE IF (WORD.EQ.'LJADD3') THEN
            LJADDT=.TRUE.
            LJADD3T=.TRUE.
            CALL READI(NADDTARGET)
            WRITE(*,'(A,I6)') 'keyword> Target cluster size is ',NADDTARGET
            IF (MOD(NATOMS,NADDTARGET).NE.0) THEN
               WRITE(*,'(A,I6)') 'keyword> ERROR, target cluster size is not a factor of the number of the atoms ',NATOMS
               STOP
            ENDIF
            LUNIT=GETUNIT()
            OPEN(LUNIT,FILE='epsilon',STATUS='OLD')
            IF (.NOT.ALLOCATED(LJADDREP)) ALLOCATE(LJADDREP(NADDTARGET,NADDTARGET))
            IF (.NOT.ALLOCATED(LJADDATT)) ALLOCATE(LJADDATT(NADDTARGET,NADDTARGET))
            DO J1=1,NADDTARGET
               DO J2=1,NADDTARGET
                  READ(LUNIT,*) LJADDREP(J2,J1), LJADDATT(J2,J1)
                  WRITE(*,'(2I6,2G20.10)') J1,J2,LJADDREP(J2,J1),LJADDATT(J2,J1)
               ENDDO
            ENDDO
            CLOSE(LUNIT)
         ELSE IF (WORD.EQ.'LJADD4') THEN
            LJADDT=.TRUE.
            LJADD4T=.TRUE.
            CALL READI(NADDTARGET)
            WRITE(*,'(A,I6)') 'keyword> Target cluster size is ',NADDTARGET
            IF (MOD(NATOMS,NADDTARGET).NE.0) THEN
               WRITE(*,'(A,I6)') 'keyword> ERROR, target cluster size is not a factor of the number of the atoms ',NATOMS
               STOP
            ENDIF
            LUNIT=GETUNIT()
            OPEN(LUNIT,FILE='epsilon',STATUS='OLD')
            IF (.NOT.ALLOCATED(LJADDREP)) ALLOCATE(LJADDREP(NADDTARGET,NADDTARGET))
            IF (.NOT.ALLOCATED(LJADDATT)) ALLOCATE(LJADDATT(NADDTARGET,NADDTARGET))
            DO J1=1,NADDTARGET
               DO J2=1,NADDTARGET
                  READ(LUNIT,*) LJADDREP(J2,J1), LJADDATT(J2,J1)
                  WRITE(*,'(2I6,2G20.10)') J1,J2,LJADDREP(J2,J1),LJADDATT(J2,J1)
               ENDDO
            ENDDO
            READ(LUNIT,*) NUMNN, TANHFAC, LJADDCUTOFF
            IF (.NOT.ALLOCATED(LJADDNN)) ALLOCATE(LJADDNN(NUMNN,2))
            DO J2=1,NUMNN
               READ(LUNIT,*) LJADDNN(J2,1), LJADDNN(J2,2)
               WRITE(*,'(2I6)') LJADDNN(J2,1), LJADDNN(J2,2)
            ENDDO
            READ(LUNIT,*) LJADDREFNORM
            WRITE(*,'(A,I6,3G20.10)') 
     &           'keyword> Number of nearest neighbours in target structures, factor for step function, cutoff, norm=',
     &                                NUMNN,TANHFAC,LJADDCUTOFF,LJADDREFNORM
            CLOSE(LUNIT)
!
! Keyword for adding a general LJ site to PY ellipsoids
! Coded by swo24 in July 2011
! Keyword: LJGSITE ("LJ general site"), used in MULTISITEPY2 in multisitepy.f90
! Syntax: LJGSITE sigma_0 epsilon
! Also requires: ljsites.xyz file as described in multisitepy.f90
         ELSE IF (WORD.EQ.'LJGSITE') THEN
! Turn on logical indicating that there is an LJ general site
            LJGSITET=.TRUE.

! Read parameters from data file
            CALL READF(LJGSITESIGMA)
            CALL READF(LJGSITEEPS)
            WRITE(*,*) "keyword> adding LJ site(s)"

         ELSE IF (WORD.EQ.'LWOTP') THEN
            LWOTPT = .TRUE.
            RBAAT  = .TRUE.
            NRBSITES = 3
            ALLOCATE(RBSITE(NRBSITES,3))
            NTSITES = NATOMS*NRBSITES/2

         ELSE IF (WORD.EQ.'LOWESTFRQ') THEN
            LOWESTFRQT=.TRUE.
! 
! MMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMMM
! 
         ELSE IF (WORD.EQ.'MACHINE') THEN
            MACHINE=.TRUE.
! sf344> macroion model
         ELSE IF (WORD.EQ.'MACROION') THEN
            MACROIONT=.TRUE.
! 
! Macrocycle
! Used for cyclic peptides with repeating sequences (e.g. cyclo-[GlyProGlyPro])
! 
         ELSE IF (WORD.EQ.'MACROCYCLE') THEN
            MACROCYCLET=.TRUE.
            IF (NITEMS.GT.1) THEN
               CALL READI(MCYCLEREPEATS)
            ENDIF
            MCYCLEPERIOD=NATOMS/MCYCLEREPEATS
            IF (MCYCLEPERIOD*MCYCLEREPEATS.NE.NATOMS) THEN
               PRINT('(A,I5,A,I5,A)'),"Warning: ring of size ",NATOMS," cannot contain ",MCYCLEREPEATS, " repeat units."
            ELSE
               PRINT('(A,I5,A,I5,A)')," keyword> Macrocycle with ",MCYCLEREPEATS," units, each comprising ",MCYCLEPERIOD," atoms."
            ENDIF

! 
! MASS ON/OFF takes steps with a fictitious "kinetic" metric   - default OFF
! 
         ELSE IF (WORD .EQ. 'MASS') THEN
            CALL READU(WW)
            IF (WW .EQ. 'ON' .OR. WW .EQ. ' ') THEN
               MASST=.TRUE.
            ELSE IF (WW .EQ. 'OFF') THEN
               MASST=.FALSE.
            ENDIF
!
!
!
            ELSE IF (WORD.EQ.'MALONALDEHYDE') THEN
               MALONALDEHYDE=.TRUE.

            ELSE IF (WORD.EQ.'MBPOL') THEN
               FRQCONV=5.123934D14
               CALL MBPOLINIT
               MBPOLT=.TRUE.

! 
! Maximum value for the smaller barrier height that is allowed to constitute a connection during the
! Dijkstra connection procedure.
! MAXMAXBARRIER specifies a maximum for the maximum barrier.
! MAXBARRIER requires both sides to be greater than MAXBARRIER to discard.
! 
         ELSE IF (WORD.EQ.'MAXBARRIER') THEN
            CALL READF(MAXBARRIER)
         ELSE IF (WORD.EQ.'MAXMAXBARRIER') THEN
            CALL READF(MAXMAXBARRIER)
! 
! MAXBFGS x1 x2 x3 x4\/}: {\it x\/} specifies the maximum allowed step length in LBFGS
! minimisations, {\it x1\/} for  normal minimisations, {\it x2\/} for Rayleigh-Ritz ratio
! minimisation, {\it x3\/} for putting structures in closest coincidence with
! bf mind} (NO LONGER USED!!), and {\it x4\/} for NEB minimisations. Default values all 0.2.
! 
         ELSE IF (WORD.EQ.'MAXBFGS') THEN
            CALL READF(MAXBFGS)
            IF (NITEMS.GT.2) CALL READF(MAXXBFGS)
            IF (NITEMS.GT.3) CALL READF(MAXMBFGS)
            IF (NITEMS.GT.4) CALL READF(MAXNEBBFGS)
            IF (NITEMS.GT.5) CALL READF(MAXINTBFGS)
! 
! The maximum number of constraints to use in the constrained potential.
! The deafult is 4.
! 
         ELSE IF (WORD.EQ.'MAXCON') THEN
            CALL READI(MAXCONUSE)
! 
! The maximum energy increase above which mylbfgs will reject a proposed step.
! 
         ELSE IF (WORD.EQ.'MAXERISE') THEN
            CALL READF(MAXERISE)
            IF (NITEMS.GT.1) CALL READF(XMAXERISE)
!
! This is the factor used in NEB/output.f90 to identify local maxima for hybrid EF ts searches.
! It multiplies EDIFFTOL.
!
         ELSE IF (WORD.EQ.'MAXIMFACTOR') THEN
            CALL READF(MAXIMFACTOR)

! Maximum number of failures allowed in a minimisation before giving up.
! 
         ELSE IF (WORD.EQ.'MAXFAIL') THEN
            IF (NITEMS.GT.1) CALL READI(NFAILMAX)
!
! Specify single connection attempt for largest gap in current Dijkstra selected chain.
! 
         ELSE IF (WORD.EQ.'MAXGAP') THEN
            MAXGAPT=.TRUE.
! 
! For the growing string double-ended connection
! method, specify a maximum number of steps allowed before another image is
! added to the growing string. Default is 1000.
! 
         ELSE IF (WORD.EQ.'MAXGROWSTEPS') THEN
            CALL READI(MAXGROWSTEPS)
! 
! Will stop the entire job if the total string
! length for the growing strings or evolving strings method goes above {\it x}
! times the total number of images. This usually means that something is going
! wrong with the string. Default is 1000.
! 
         ELSE IF (WORD.EQ.'MAXLENPERIM') THEN
            CALL READF(MAXLENPERIM)
! 
! Specifies the maximum value that the maximum step size
! is allowed to rise to. The default value is $0.5$.
! 
         ELSE IF (WORD.EQ.'MAXMAX') THEN
            CALL READF(MAXMAX)
! 
! MAXSTEP n specifies the maximum step size in real units      - default n=0.2
! Applies to eigenvector-following and steepest-descent calculations.
! 
         ELSE IF (WORD.EQ.'MAXSTEP') THEN
            CALL READF(MXSTP)
! 
! Maximum ts energy that is allowed to constitute a connection during the
! Dijkstra connection procedure.
! 
         ELSE IF (WORD.EQ.'MAXTSENERGY') THEN
            CALL READF(MAXTSENERGY)
         ELSE IF (WORD.EQ.'MCBIAS') THEN
            MCBIAST=.TRUE.
         ELSE IF (WORD.EQ.'MCMERGE') THEN
            IF (NITEMS.GT.1) CALL READI(MCMERGES)
            IF (NITEMS.GT.2) CALL READI(MCMERGEQ)
! 
! MCPATH - MC sample path.xyz pathway 
! 
!           MCPATHTEMP     MC canonical temperature
!           MCPATHDMAX     maximum allowed distance from a reference structure
!           MCPATHSTART    start sampling from 1st min, ts or 2nd min for -1,0, 1
!           MCPATHSTEP     maximum initial MC Cartesian step size
!           MCPATHACCRATIO target acceptance ratio
!           MCPATHBINS     number of bins for order parameter histograms
!           MCPATHEQUIL    number of equilibration steps
!           MCPATHSTEPS    number of MC production sampling steps
!           MCPATHPRTFRQ   print every MCPATHPRTFRQ steps
!           BIASFAC        starting value for exponent in bias function
!           MCADDDEV       add reference structures if bias potential deviation exceeeds MCADDDEV
!           MCPATHQMIN     smallest value of order parameter for binning
!           MCPATHQMAX     largest value of order parameter for binning
!           MCPATHMAXDEV   largest allowed gap between reference strutures

         ELSE IF (WORD.EQ.'MCPATH') THEN
            MCPATHT=.TRUE.
            IF (NITEMS.GT.1) CALL READF(MCPATHTEMP)
            IF (NITEMS.GT.2) CALL READF(MCPATHDMAX)
            IF (NITEMS.GT.3) CALL READI(MCPATHSTART)
            IF (NITEMS.GT.4) CALL READF(MCPATHSTEP)
            IF (NITEMS.GT.5) CALL READF(MCPATHACCRATIO)
            IF (NITEMS.GT.6) CALL READI(MCPATHBINS)
            IF (NITEMS.GT.7) CALL READI(MCPATHEQUIL)
            IF (NITEMS.GT.8) CALL READI(MCPATHSTEPS)
            IF (NITEMS.GT.9) CALL READI(MCPATHPRTFRQ)
            IF (NITEMS.GT.10) CALL READF(BIASFAC)
            IF (NITEMS.GT.11) CALL READF(MCADDDEV)
            IF (NITEMS.GT.12) CALL READF(MCPATHQMIN)
            IF (NITEMS.GT.13) CALL READF(MCPATHQMAX)
            IF (NITEMS.GT.14) CALL READF(MCPATHADDREF)
            IF (NITEMS.LE.14) THEN
               PRINT '(A)','keywords> ERROR *** insufficient arguments for MCPATH keyword'
               STOP
            ENDIF
            IF (NITEMS.GT.15) CALL READF(MCPATHTOL)
         ELSE IF (WORD.EQ.'MCPATH2') THEN
            MCPATH2T=.TRUE.
            IF (NITEMS.GT.1) CALL READF(MCPATHTEMP)
            IF (NITEMS.GT.2) CALL READF(MCPATHDMAX)
            IF (NITEMS.GT.3) CALL READI(MCPATHSTART)
            IF (NITEMS.GT.4) CALL READF(MCPATHSTEP)
            IF (NITEMS.GT.5) CALL READF(MCPATHACCRATIO)
            IF (NITEMS.GT.6) CALL READI(MCPATHBINS)
            IF (NITEMS.GT.7) CALL READI(MCPATHEQUIL)
            IF (NITEMS.GT.8) CALL READI(MCPATHSTEPS)
            IF (NITEMS.GT.9) CALL READI(MCPATHPRTFRQ)
            IF (NITEMS.GT.10) CALL READI(MCPATHBLOCK)
            IF (NITEMS.GT.11) CALL READI(MCPATHOVER)
            IF (NITEMS.GT.12) CALL READF(MCPATHQMIN)
            IF (NITEMS.GT.13) CALL READF(MCPATHQMAX)
            IF (NITEMS.GT.14) CALL READF(MCPATHADDREF)
            IF (NITEMS.LE.14) THEN
               PRINT '(A)','keywords> ERROR *** insufficient arguments for MCPATH2 keyword'
               STOP
            ENDIF
            IF (NITEMS.GT.15) CALL READF(MCPATHNEGLECT)
            IF (NITEMS.GT.16) CALL READF(MCPATHTOL)
         ELSE IF (WORD.EQ.'MCPATHDOBLOCK') THEN
            CALL READI(MCPATHDOBLOCK)
         ELSE IF (WORD.EQ.'PBS') THEN
            PBST=.TRUE.
         ELSE IF (WORD.EQ.'MCPATHGW') THEN
            IF (NITEMS.GT.1) CALL READF(MCPATHGWS)
            IF (NITEMS.GT.2) CALL READF(MCPATHGWQ)
         ELSE IF (WORD.EQ.'MCPATHTS') THEN
            IF (NITEMS.GT.1) CALL READI(MCPATHTS)
            IF (NITEMS.GT.2) CALL READI(MCPATHSCHECK)
         ELSE IF (WORD.EQ.'MECCANO') THEN
            MECCANOT=.TRUE.
            CALL READF(MECIMDENS) ! now an image density
            CALL READI(MECMAXIMAGES)  ! maximum number of images
            CALL READF(MECITDENS) ! iteration density
            CALL READI(MECMAXIT)  ! maximum number of iterations
            CALL READF(MECLAMBDA)
            CALL READF(MECDIST)
            CALL READF(MECRMSTOL)
            CALL READF(MECSTEP)
            CALL READF(MECDGUESS)
            CALL READI(MECUPDATE)
         ELSE IF (WORD.EQ.'MIE_FIELD') THEN
            MIEFT=.TRUE.
            CALL READA(MIEF_FILENAME)            
            IF(NITEMS.GT.2) THEN
               MIEF_CUTT=.TRUE.
               CALL READF(MIEF_RCUT)
            ENDIF
            IF(NITEMS.GT.3) THEN
               MIEF_PBCT=.TRUE.
               CALL READF(MIEF_BOX(1))
               CALL READF(MIEF_BOX(2))
               CALL READF(MIEF_BOX(3))
            ENDIF            
         ELSE IF (WORD.EQ.'MINMAX') THEN
            CALL READF(MINMAX)
! ELSE IF (WORD.EQ.'MINBM') THEN
! MINBMT = .TRUE.
! CALL READI(MINBMNSAMP)

         ! sn402: temporary keyword to activate the metric tensor formulation for rigid body normal modes.
         ! When I'm happy that this is working, I'll retire the old formulation altogether and make this the
         ! default.
         ELSE IF (WORD.EQ.'METRICTENSOR') THEN
            METRICTENSOR = .TRUE.


         ELSE IF (WORD.EQ.'MINBACKTCUT') THEN
            CALL READF(MINBACKTCUT)
!
! Three layer neural network (multilayer perceptron) with
! MLPIN inputs (columns per data item)
! MLPOUT outputs
! MLPHIDDEN hidden nodes
! MLPDATA data lines in MLPdata file (last column MLPIN+1 for correct outputs, numbered one to MLPOUT)
! MLPLAMBDA coefficient for regularisation
! MLPDATSHIFT - starting data item to allow reading from multiple positions
!
         ELSE IF ((WORD.EQ.'MLP3').OR.(WORD.EQ.'MLPB3')) THEN
            MLP3T=.TRUE.
            IF (WORD.EQ.'MLPB3') MLPB3T=.TRUE.
            CALL READI(MLPIN)
            CALL READI(MLPHIDDEN)
            CALL READI(MLPOUT)
            CALL READI(MLPDATA)
            IF (NITEMS.GT.5) CALL READF(MLPLAMBDA)
            IF (NITEMS.GT.6) CALL READI(MLPDATSTART)
            IF (WORD.EQ.'MLPB3') THEN
               WRITE(*,'(A,5I8,G20.10)') 'MLPB3 potential with bias node and Nin, Nhidden, Nout, Ndata, startdata, lambda=',
     &                                   MLPIN,MLPHIDDEN,MLPOUT,MLPDATA,MLPDATSTART,MLPLAMBDA
               NMLP=MLPHIDDEN*(MLPIN+MLPOUT)+1
            ELSE
               WRITE(*,'(A,5I8,G20.10)') 'MLP3 potential with Nin, Nhidden, Nout, Ndata, startdata, lambda=',
     &                                   MLPIN,MLPHIDDEN,MLPOUT,MLPDATA,MLPDATSTART,MLPLAMBDA
               NMLP=MLPHIDDEN*(MLPIN+MLPOUT)
            ENDIF
            IF (NMLP.NE.NATOMS) THEN
               PRINT '(A,2I8)', 'keywords> ERROR *** NATOMS,NMLP=',NATOMS,NMLP
               STOP
            ENDIF
            LUNIT=GETUNIT()
            OPEN(LUNIT,FILE='MLPdata',STATUS='OLD')
            ALLOCATE(MLPDAT(MLPDATA,MLPIN),MLPOUTCOME(MLPDATA),MLPMEAN(MLPIN))
            MLPMEAN(1:MLPIN)=0.0D0
            DO J1=1,MLPDATA+MLPDATSTART-1
               IF (J1.GE.MLPDATSTART) THEN
                  J2=J1-MLPDATSTART+1
                  READ(LUNIT,*) MLPDAT(J2,1:MLPIN),MLPOUTCOME(J2)
                  MLPOUTCOME(J2)=MLPOUTCOME(J2)+1 ! to shift the range from 1 instead of zero
                  DO J3=1,MLPIN
                     MLPMEAN(J3)=MLPMEAN(J3)+ABS(MLPDAT(J2,J3))
                  ENDDO
               ELSE
                  READ(LUNIT,*) DUMMY
               ENDIF
            ENDDO
            CLOSE(LUNIT)
            IF (MLPNORM) THEN
               MLPMEAN(1:MLPIN)=MLPMEAN(1:MLPIN)/MLPDATA
               WRITE(*,'(A)') 'keyword> Rescaling inputs by mean absolute values:'
               WRITE(*,'(6G20.10)') MLPMEAN(1:MLPIN)
               DO J1=1,MLPIN
                  MLPDAT(1:MLPDATA,J1)=MLPDAT(1:MLPDATA,J1)/MLPMEAN(J1)
               ENDDO
            ENDIF
            DEALLOCATE(MLPMEAN)
            MLPDONE=.TRUE.
         ELSE IF (WORD.EQ.'MLPB3NEW') THEN
            MLP3T=.TRUE.
            MLPB3T=.TRUE.
            MLPB3NEWT=.TRUE.
            CALL READI(MLPIN)      ! number of inputs (data items after outcome)
            CALL READI(MLPSTART) ! starting position in data list, not counting outcome
            CALL READI(MLPHIDDEN)
            CALL READI(MLPOUT)
            CALL READI(MLPDATA)
            IF (NITEMS.GT.5) CALL READF(MLPLAMBDA)
            WRITE(*,'(A,5I8,G20.10)') ' keywords> MLP3 new potential bias nodes and Nin, Ninstart, Nhidden, Nout, Ndata, lambda=',
     &                                MLPIN,MLPSTART,MLPHIDDEN,MLPOUT,MLPDATA,MLPLAMBDA
            NMLP=MLPHIDDEN*(MLPIN+MLPOUT)+1
            IF (NMLP.NE.NATOMS) THEN
               PRINT '(A,2I8)', 'keywords> ERROR *** NATOMS,NMLP=',NATOMS,NMLP
               STOP
            ENDIF
   
            LUNIT=GETUNIT()
            OPEN(LUNIT,FILE='MLPdata',STATUS='OLD')
            ALLOCATE(MLPDAT(MLPDATA,MLPIN),MLPOUTCOME(MLPDATA),MLPMEAN(MLPIN))
            MLPMEAN(1:MLPIN)=0.0D0
            DO J1=1,MLPDATA
               READ(LUNIT,*) MLPOUTCOME(J1),(DUMMY,J2=1,MLPSTART-1),MLPDAT(J1,1:MLPIN)
               MLPOUTCOME(J1)=MLPOUTCOME(J1)+1 ! to shift the range from 0 to from 1
               DO J2=1,MLPIN
                  MLPMEAN(J2)=MLPMEAN(J2)+ABS(MLPDAT(J1,J2))
               ENDDO
            ENDDO
            CLOSE(LUNIT)
            IF (MLPNORM) THEN
               MLPMEAN(1:MLPIN)=MLPMEAN(1:MLPIN)/MLPDATA
               WRITE(*,'(A)') 'keyword> Rescaling inputs by mean absolute values:'
               WRITE(*,'(6G20.10)') MLPMEAN(1:MLPIN)
               DO J1=1,MLPIN
                  MLPDAT(1:MLPDATA,J1)=MLPDAT(1:MLPDATA,J1)/MLPMEAN(J1)
               ENDDO
            ENDIF
            DEALLOCATE(MLPMEAN)
            MLPDONE=.TRUE.
         ELSE IF (WORD.EQ.'NOREGBIAS') THEN
            NOREGBIAS=.TRUE.
         ELSE IF ((WORD.EQ.'MLPVB3').OR.(WORD.EQ.'MLPVB3NN')) THEN
            MLPVB3T=.TRUE.
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
            WRITE(*,'(A,5I8,G20.10)') ' keywords> MLP3 vector bias nodes and Nin, Ninstart, Nhidden, Nout, Ndata, lambda=',
     &                                MLPIN,MLPSTART,MLPHIDDEN,MLPOUT,MLPDATA,MLPLAMBDA
            NMLP=MLPHIDDEN*(MLPIN+MLPOUT)+MLPHIDDEN+MLPOUT
            IF (NMLP.NE.NATOMS) THEN
               PRINT '(A,2I8)', 'keywords> ERROR *** NATOMS,NMLP=',NATOMS,NMLP
               STOP
            ENDIF
   
            LUNIT=GETUNIT()
            OPEN(LUNIT,FILE='MLPdata',STATUS='OLD')
            ALLOCATE(MLPDAT(MLPDATA,MLPIN),MLPOUTCOME(MLPDATA),MLPMEAN(MLPIN))
            MLPMEAN(1:MLPIN)=0.0D0
            DO J1=1,MLPDATA
               READ(LUNIT,*) MLPOUTCOME(J1),(DUMMY,J2=1,MLPSTART-1),MLPDAT(J1,1:MLPIN)
               MLPOUTCOME(J1)=MLPOUTCOME(J1)+1 ! to shift the range from 0 to from 1
               DO J2=1,MLPIN
                  MLPMEAN(J2)=MLPMEAN(J2)+ABS(MLPDAT(J1,J2))
               ENDDO
            ENDDO
            CLOSE(LUNIT)
            IF (MLPNORM) THEN
               MLPMEAN(1:MLPIN)=MLPMEAN(1:MLPIN)/MLPDATA
               WRITE(*,'(A)') 'keyword> Rescaling inputs by mean absolute values:'
               WRITE(*,'(6G20.10)') MLPMEAN(1:MLPIN)
               DO J1=1,MLPIN
                  MLPDAT(1:MLPDATA,J1)=MLPDAT(1:MLPDATA,J1)/MLPMEAN(J1)
               ENDDO
            ENDIF
            DEALLOCATE(MLPMEAN)
            MLPDONE=.TRUE.
!
! Nearest-neighbour model - freeze hidden node to i/o node weights except for "nearest-neighbours".
!
!
! Variables are ordered
! w^2_{jk} at (j-1)*MLPIN+k
!   up to MLPHIDDEN*MLPIN, then
! w^1_{ij} at MLPHIDDEN*MLPIN + (i-1)*MLPHIDDEN+j
!   up to MLPHIDDEN*MLPIN + MLPOUT*MLPHIDDEN
! w^bh_j at MLPHIDDEN*(MLPIN+MLPOUT)+1 to MLPHIDDEN*(MLPIN+MLPOUT)+MLPHIDDEN
! w^bo_i at MLPHIDDEN*(MLPIN+MLPOUT)+MLPHIDDEN+1 to MLPHIDDEN*(MLPIN+MLPOUT)+MLPHIDDEN+MLPOUT
!
            IF (WORD.EQ.'MLPVB3NN') THEN
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
     &                  ' distance ',MLPDISTHO(J2)
                     J3=MLPHIDDEN*MLPIN+(MLPINDEXO(J2)-1)*MLPHIDDEN+J1
                     FROZEN(J3)=.FALSE.
                     NFREEZE=NFREEZE-1
                  ENDDO
               ENDDO
               DEALLOCATE( MLPDISTHI, MLPDISTHO, MLPINDEXI, MLPINDEXO)
            ENDIF
         ELSE IF (WORD.EQ.'MLPNEWREG') THEN
            MLPNEWREG=.TRUE.
         ELSE IF (WORD.EQ.'MLPPROB') THEN
            MLPPROB=.TRUE.
            IF (NITEMS.GT.1) CALL READI(MLPPROBPOS)
            WRITE(*,'(A,I6)') 
     &    'keywords> Will calculate ROC and AUC values assuming positive outcome (counting from 1) is indexed ',MLPPROBPOS
! 
! MLPNORM directs OPTIM to rescale the input data columns by dividing each one by the
! average of the mean magnitude
! 
         ELSE IF (WORD.EQ.'MLPNORM') THEN
            MLPNORM=.TRUE.
            IF (MLPDONE) THEN
               WRITE(*,'(A)') 'keyword> ERROR *** please put MLPNORM before MLP keyword in odata to ensure correct io'
               STOP

!              LUNIT=GETUNIT()
!              OPEN(LUNIT,FILE='MLPdata',STATUS='OLD')
!              ALLOCATE(MLPMEAN(MLPIN))
!              MLPMEAN(1:MLPIN)=0.0D0
!              DO J1=1,MLPDATA+MLPDATSTART-1
!                 IF (J1.GE.MLPDATSTART) THEN
!                    J2=J1-MLPDATSTART+1
!                    READ(LUNIT,*) MLPDAT(J2,1:MLPIN),MLPOUTCOME(J2)
!                    MLPOUTCOME(J2)=MLPOUTCOME(J2)+1 ! to shift the range to 1 to 4
!                    DO J3=1,MLPIN
!                       MLPMEAN(J3)=MLPMEAN(J3)+ABS(MLPDAT(J2,J3))
!                    ENDDO
!                 ELSE
!                    READ(LUNIT,*) DUMMY
!                 ENDIF
!              ENDDO

!              CLOSE(LUNIT)
!              MLPMEAN(1:MLPIN)=MLPMEAN(1:MLPIN)/MLPDATA
!              PRINT '(A)','keyword> Rescaling inputs by mean absolute values:'
!              PRINT '(6G20.10)',MLPMEAN(1:MLPIN)
!              DO J1=1,MLPIN
!                 MLPDAT(1:MLPDATA,J1)=MLPDAT(1:MLPDATA,J1)/MLPMEAN(J1)
!              ENDDO
!              DEALLOCATE(MLPMEAN)
            ENDIF
         ELSE IF (WORD.EQ.'MLQ') THEN
            MLQT=.TRUE.
            CALL READI(MLQIN)      ! number of inputs (data items after outcome)
            CALL READI(MLQSTART) ! starting position in data list, not counting outcome
            CALL READI(MLQOUT)
            CALL READI(MLQDATA)
            IF (NITEMS.GT.4) CALL READF(MLQLAMBDA)
            WRITE(*,'(A,4I8,G20.10)') ' keywords> MLQ Nin, Ninstart, Nout, Ndata, lambda=',
     &                                MLQIN,MLQSTART,MLQOUT,MLQDATA,MLQLAMBDA
            NMLQ=MLQOUT*(1+(MLQIN*(MLQIN+3))/2)
            WRITE(*,'(A,5I8,G20.10)') ' keywords> MLQ variables=',NMLQ
            IF (NMLQ.NE.NATOMS) THEN
               PRINT '(A,2I8)', 'keywords> ERROR *** NATOMS,NMLQ=',NATOMS,NMLQ
               STOP
            ENDIF

            LUNIT=GETUNIT()
            OPEN(LUNIT,FILE='MLQdata',STATUS='OLD')
            ALLOCATE(MLQDAT(MLQDATA,MLQIN),MLQOUTCOME(MLQDATA),MLQMEAN(MLQIN))
            MLQMEAN(1:MLQIN)=0.0D0
            DO J1=1,MLQDATA
               READ(LUNIT,*) MLQOUTCOME(J1),(DUMMY,J2=1,MLQSTART-1),MLQDAT(J1,1:MLQIN)
               MLQOUTCOME(J1)=MLQOUTCOME(J1)+1 ! to shift the range from 0 to from 1
               DO J2=1,MLQIN
                  MLQMEAN(J2)=MLQMEAN(J2)+ABS(MLQDAT(J1,J2))
               ENDDO
            ENDDO
            CLOSE(LUNIT)
            IF (MLQNORM) THEN
               MLQMEAN(1:MLQIN)=MLQMEAN(1:MLQIN)/MLQDATA
               WRITE(*,'(A)') 'keyword> Rescaling inputs by mean absolute values:'
               WRITE(*,'(6G20.10)') MLQMEAN(1:MLQIN)
               DO J1=1,MLQIN
                  MLQDAT(1:MLQDATA,J1)=MLQDAT(1:MLQDATA,J1)/MLQMEAN(J1)
               ENDDO
            ENDIF
            DEALLOCATE(MLQMEAN)
            MLQDONE=.TRUE.
         ELSE IF (WORD.EQ.'MLQPROB') THEN
            MLQPROB=.TRUE.
! 
! MLQNORM directs OPTIM to rescale the input data columns by dividing each one by the
! average of the mean magnitude
! 
         ELSE IF (WORD.EQ.'MLQNORM') THEN
            MLQNORM=.TRUE.
            IF (MLQDONE) THEN
               LUNIT=GETUNIT()
               OPEN(LUNIT,FILE='MLQdata',STATUS='OLD')
               ALLOCATE(MLQMEAN(MLQIN))
               MLPMEAN(1:MLQIN)=0.0D0
               DO J1=1,MLQDATA+MLQDATSTART-1
                  IF (J1.GE.MLQDATSTART) THEN
                     J2=J1-MLQDATSTART+1
                     READ(LUNIT,*) MLQDAT(J2,1:MLQIN),MLQOUTCOME(J2)
                     MLQOUTCOME(J2)=MLQOUTCOME(J2)+1 ! to shift the range to 1 to 4
                     DO J3=1,MLPIN
                        MLQMEAN(J3)=MLQMEAN(J3)+ABS(MLQDAT(J2,J3))
                     ENDDO
                  ELSE
                     READ(LUNIT,*) DUMMY
                  ENDIF
               ENDDO
   
               CLOSE(LUNIT)
               MLPMEAN(1:MLQIN)=MLQMEAN(1:MLQIN)/MLQDATA
               PRINT '(A)','keyword> Rescaling inputs by mean absolute values:'
               PRINT '(6G20.10)',MLQMEAN(1:MLQIN)
               DO J1=1,MLPIN
                  MLQDAT(1:MLQDATA,J1)=MLQDAT(1:MLQDATA,J1)/MLQMEAN(J1)
               ENDDO
               DEALLOCATE(MLQMEAN)
            ENDIF
! 
! MODE n  specifies the eigenvector to follow                  - default n=0
! 

         ELSE IF (WORD.EQ.'MODE') THEN
            CALL READI(IVEC)
            IF (NITEMS.GT.2) THEN
               CALL READI(IVEC2)
            ELSE
               ! IVEC2=IVEC
            ENDIF
! 
! MODE n  specifies an eigenvector to follow downhill          - default n=0
! 
         ELSE IF (WORD.EQ.'MODEDOWN') THEN
            MODEDOWNT=.TRUE.
            CALL READI(IVEC)
            IF (NITEMS.GT.2) THEN
               CALL READI(IVEC2)
            ELSE
               ! 
               ! Following the softest mode downhill after the first step by default is probably not
               ! what we want for MODEDOWN.
               ! 
               IVEC2=IVEC
            ENDIF
! 
! Attempt to morph between endpoints by taking steps towards or
! away from the endpoint finish.
! 
         ELSE IF (WORD.EQ.'MORPH') THEN
            MORPHT=.TRUE.
            CALL READF(MORPHMXSTP)
            CALL READI(MNBFGSMAX1)
            CALL READI(MNBFGSMAX2)
            CALL READF(MORPHEMAX)
            CALL READF(MORPHERISE)
            CALL READI(MSTEPS)
            IF (MAXTSENERGY.EQ.1.0D100) MAXTSENERGY=MORPHEMAX
! 
! Movie dump for Paul Mortenson`s amber
! 
! ELSE IF (WORD.EQ.'MOVIE') THEN
! MOVIE=.TRUE.
! OPEN (UNIT=27, FILE='amber.movie', STATUS='UNKNOWN')
! 
! MSEVB parameters - probably shouldn`t be changed on a regular basis
! 
         ELSE IF (WORD.EQ.'MSEVBPARAMS') THEN
            IF (NITEMS.GT.1) CALL READI(shellsToCount)
            IF (NITEMS.GT.2) CALL READF(maxHbondLength)
            IF (NITEMS.GT.3) CALL READF(minHbondAngle)
            IF (NITEMS.GT.4) CALL READF(OOclash_sq)

         ELSE IF (WORD.EQ.'MSSTOCK') THEN

            MSSTOCKT = .TRUE.
            RBAAT    = .TRUE.
            CALL READI(NRBSITES)
            ALLOCATE(RBSITE(NRBSITES,3))
            ALLOCATE(RBSTLA(NRBSITES,3))
            ALLOCATE(DPMU(NRBSITES))

            NTSITES = NATOMS*NRBSITES/2

            DO J1 = 1, NRBSITES
               CALL READF(DPMU(J1))
            ENDDO

            IF (NITEMS > (NRBSITES+2)) THEN
               CALL READF(EFIELD)
               EFIELDT = .TRUE.
            ENDIF

            CALL DEFMULTSTOCK()
            IF (PERMDIST) THEN ! correct all permutations allowed if perm.allow is not given explicitly
               IF (NPERMSIZE(1).EQ.NATOMS) NPERMSIZE(1)=NATOMS/2
            ENDIF
! 
! MULTIJOB - repeat the calculation for more than one set of input data
! 
         ELSE IF (WORD .EQ. 'MULTIJOB') THEN
            MULTIJOBT=.TRUE.
            CALL READA(MULTISTART)
            MULTISUNIT=GETUNIT()
            OPEN(MULTISUNIT,FILE=TRIM(ADJUSTL(MULTISTART)),STATUS='OLD')
            PRINT '(A)',' keywords> Multiple jobs will be run, reading additional starting coordinates from file '
     &      // TRIM(ADJUSTL(MULTISTART))
            MULTIFUNIT=-1
            IF (NITEMS.GT.2) THEN
               CALL READA(MULTIFINISH)
               MULTIFUNIT=GETUNIT()
               OPEN(MULTIFUNIT,FILE=TRIM(ADJUSTL(MULTIFINISH)),STATUS='OLD')
               PRINT '(A)',' keywords>                            reading additional finish coordinates from file '
     &         // TRIM(ADJUSTL(MULTIFINISH))
            ENDIF

! MULTIJOB_MACHINE - as MULTIJOB but coordinates are read from an unformatted file. Coordinates do not need to be read consecutively if the file
! is very large.

         ELSE IF (WORD .EQ. 'MULTIJOB_MACHINE') THEN
            MULTIJOBT=.TRUE.
            MULTIJOB_MACHINET = .TRUE.
            CALL READA(MULTISTART)
            CALL READI(NDOF)
            MULTISUNIT=GETUNIT()
            PRINT '(A)',' keywords> Multiple jobs will be run, reading additional starting coordinates from unformatted file '
     &      // TRIM(ADJUSTL(MULTISTART))
            MULTIFUNIT=-1
            IF (NITEMS.GT.2) THEN
               CALL READA(MULTIFINISH)
               IF( TRIM(ADJUSTL(MULTIFINISH)).NE.'NONE') THEN
                  MULTIFUNIT=GETUNIT()
                  OPEN(MULTIFUNIT,FILE=TRIM(ADJUSTL(MULTIFINISH)),ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='OLD',RECL=8*NDOF)
                  PRINT '(A)',' keywords>                            reading additional finish coordinates from unformatted file '
     &            // TRIM(ADJUSTL(MULTIFINISH))
               ENDIF
            ENDIF
            IF (NITEMS.GT.3) THEN
               CALL READI(MULTI_COUNT)   ! The index of the first configuration to be read in
            ENDIF
            IF (NITEMS.GT.4) THEN
               CALL READI(MULTI_LAST)    ! The index of the last configuration to be read in
            ENDIF
            IF (NITEMS.GT.5) THEN
               CALL READI(MULTI_STEP)    ! 1+(the number of configurations to skip each time)
            ENDIF
            WRITE(*,*) "keywords> Reading from configuration ", MULTI_COUNT, "to ", MULTI_LAST
            WRITE(*,*) "keywords> skipping ", MULTI_STEP-1, "each time."

         ELSE IF (WORD .EQ. 'MULTIPOT') THEN ! Activate the multiple-potential scheme
            
            MULTIPOTT = .TRUE.
            CALL MULTIPOT_INITIALISE

         ELSE IF (WORD .EQ. 'MULTISITEPY') THEN
! Syntax: MULTISITEPY sig_0 eps_0 [cut] [XYZ boxx boxy boxz]
! Notes: The cutoff length is the raw length. It is not scaled
! by the PY sigma_0 since it is also used for the LJ potential.
! The box length is in units of cutoff distance, so it should be
! >= 2.

            MULTISITEPYT = .TRUE.
! ANGLEAXIS2 = .TRUE.
            RBAAT = .TRUE.
            CALL READF(PYSIGNOT)
            CALL READF(PYEPSNOT)

! Specify cutoff for potential in absolute units
            IF (NITEMS.GT.3) THEN
               CALL READF(PCUTOFF)
               PARAMONOVCUTOFF=.TRUE.
               WRITE(*,*) "multisitepy cutoff: ", PCUTOFF
            ENDIF
! Specify periodic boundary conditions (PBCs)
            IF (NITEMS.GT.4) THEN
               ! control which dimensions have periodic boundaries with a string 'XYZ', always put x before y before z.
               ! eg ...  Xz 20 30  specifies PBC on X and Z directions.  The X box size will be 20, the Z box size 30
               CALL READA(PBC)
               BOXLX=0
               BOXLY=0
               BOXLZ=0
               IF (SCAN(PBC,'Xx').NE.0) THEN
                  PARAMONOVPBCX=.TRUE.
                  CALL READF(BOXLX)
                  BOXLX = BOXLX*PCUTOFF
                  WRITE(*,*) "PBC X:",BOXLX
               ENDIF
               IF (SCAN(PBC,'Yy').NE.0) THEN
                  PARAMONOVPBCY=.TRUE.
                  CALL READF(BOXLY)
                  BOXLY = BOXLY*PCUTOFF
                  WRITE(*,*) "PBC Y:",BOXLY
               ENDIF
               IF (SCAN(PBC,'Zz').NE.0) THEN
                  PARAMONOVPBCZ=.TRUE.
                  CALL READF(BOXLZ)
                  BOXLZ = BOXLZ*PCUTOFF
                  WRITE(*,*) "PBC Z:",BOXLZ
               ENDIF
            ENDIF

! 
! NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
! 
! Specifies a tight-binding potential for sodium, silver and lithium
! 
         ELSE IF (WORD.EQ.'NATB') THEN
            NATBT=.TRUE.
         ELSE IF (WORD.EQ.'NATINT') THEN
            NATINT = .TRUE.
         ELSE IF (WORD.EQ.'RIGIDBONDS') THEN
            RIGIDBONDS = .TRUE.
         ELSE IF (WORD.EQ.'BONDSFROMFILE') THEN
            BONDSFROMFILE = .TRUE.
! 
         ELSE IF (WORD.EQ.'NCAP') THEN
            NCAPT    = .TRUE.
            RBAAT    = .TRUE.
            HEIGHT   = 0.5D0
            CALL READF(CAPSRHO)
            CALL READF(CAPSEPS2)
            CALL READF(CAPSRAD)
            IF (NITEMS.GT.4) CALL READF(HEIGHT)
            NRBSITES = 6
            ALLOCATE(RBSITE(NRBSITES,3))
            NTSITES = NATOMS*NRBSITES/2
! 
! Nudged elastic band calculation using a maximum of NSTEPNEB steps with
! NIMAGE images and RMS convergence criterion RMSNEB.
! 
         ELSE IF (WORD.EQ.'NEB') THEN
            NEBT=.TRUE.
            IF (NITEMS.GT.1) CALL READI(NSTEPNEB)
            IF (NITEMS.GT.2) CALL READI(NIMAGE)
            IF (NITEMS.GT.3) CALL READF(RMSNEB)
            IF (NITEMS.GT.3) CALL READI(NEBMAG)
         ELSE IF (WORD == 'NEBK') THEN
            CALL READF(NEBK)
            NEBKFINAL=NEBK
            NEBKINITIAL=NEBK
            NEBFACTOR=1.01D0
            IF (NITEMS.GT.2) CALL READF(NEBKFINAL)
            IF (NITEMS.GT.3) CALL READF(NEBFACTOR)

         ELSE IF (WORD.EQ.'NEBMAXERISE') THEN
            CALL READF(NEBMAXERISE)
! 
! 
! Read dneb guess images from file GUESSFILE, default name guess.xyz
! 
         ELSE IF ((WORD.EQ.'NEBREADGUESS').OR.(WORD.EQ.'NEBREADGUESSPERM')) THEN
            IF (WORD.EQ.'NEBREADGUESSPERM') PERMGUESS=.TRUE.
            READGUESS=.TRUE.
            NEWNEBT=.TRUE.
            FCD=.TRUE.
            IF (NITEMS.GT.1) CALL READA(GUESSFILE)
            FILENAME=GUESSFILE
            LUNIT=GETUNIT()
            OPEN(UNIT=LUNIT,FILE=FILENAME,STATUS='old')
            NIMAGE=0
!
! Here we skip two lines, allowing for the second line to be blank.
!
753         READ(LUNIT,*,END=864)
            READ(LUNIT,*)
            DO J1=1,NATOMS
               READ(LUNIT,*) DUMMYS
            ENDDO
            NIMAGE=NIMAGE+1
            GOTO 753
864         NIMAGE=NIMAGE-2
            PRINT '(A,I6,A)','nnutils> There are ',NIMAGE+2,' images in file ' // TRIM(ADJUSTL(GUESSFILE))
            CLOSE(LUNIT)
! 
! Reseed DNEB images if they exceed a certain energy.
! 
         ELSE IF (WORD.EQ.'NEBRESEED') THEN
            NEBRESEEDT=.TRUE.
            IF (NITEMS.GT.1) CALL READI(NEBRESEEDINT)
            IF (NITEMS.GT.2) CALL READF(NEBRESEEDEMAX)
            IF (NITEMS.GT.3) CALL READF(NEBRESEEDBMAX)
            IF (NITEMS.GT.4) CALL READF(NEBRESEEDDEL1)
            IF (NITEMS.GT.5) CALL READI(NEBRESEEDPOW1)
            IF (NITEMS.GT.6) CALL READF(NEBRESEEDDEL2)
            IF (NITEMS.GT.7) CALL READI(NEBRESEEDPOW2)
         ELSE IF (WORD == 'NEWCONNECT') THEN
            NEWCONNECTT = .TRUE.
            CONNECTT = .TRUE.
            OPTIMIZETS = .TRUE.
            IF (NITEMS.GT.1) CALL READI(NCONMAX)
            IF (NITEMS.GT.2) CALL READI(NTRIESMAX)
            IF (NITEMS.GT.3) CALL READF(IMAGEDENSITY)
            IF (NITEMS.GT.4) CALL READF(ITERDENSITY)
            IF (NITEMS.GT.5) CALL READI(IMAGEMAX)
            IF (NITEMS.GT.6) CALL READF(IMAGEINCR)
            ! Currently, the next line clashes with the equivalent parameter in NEWNEB - both are being read into the
            ! same variable, so only one value will ever be used.
            IF (NITEMS.GT.7) CALL READF(RMSTOL)
            ! The following line ought to fix the problem, but I've not implemented it properly yet.
            !IF (NITEMS.GT.7) CALL READF(NEWCONNECT_RMSTOL)
            IMAGEMAX=MAX(IMAGEMAX,NIMAGE+2)
! 
! If NEWCONNECT is specified the values read below are only used for the first cycle.
! If NEWNEB is used with OLDCONNECT then the values read on the NEWNEB line are
! used in every cycle. If NEWCONNECT is used then a NEWNEB line isn;t necessary.
! 
         ELSE IF (WORD == 'NEWNEB') THEN
            NEWNEBT=.TRUE.
            FCD=.TRUE.
            IF (NITEMS.GT.1) CALL READI(NNNIMAGE)
            IF (NITEMS.GT.2) CALL READI(NITERMAX)
            IF (NITEMS.GT.3) CALL READF(RMSTOL)
            IF (READGUESS) NNNIMAGE=NIMAGE
! 
! NGUESS specifies the number of transition state guesses tried in GUESSTS for CHARMM
! before switching back to NEB or NEWNEB.
! 
         ELSE IF (WORD.EQ.'NGUESS') THEN
            CALL READI(NGUESS)
! 
! NIMET  turns ON/OFF the embedded diatomics-in-molecules Ni + H PES (EAM4)
! default OFF
! 
         ELSE IF (WORD.EQ.'NIMET') THEN

            NIMET   = .TRUE.
            PRINT '(A)',' keywords> Nimet is true'

! 
! NIHEAM7 turns ON/OFF the embedded diatomics-in-molecules Ni + H PES (EAM7)
! default OFF
! 
         ELSE IF (WORD.EQ.'NIHEAM7') THEN

            NIHEAM7T   = .TRUE.
            PRINT '(A)',' keywords> Ni plus hydrogen EAM7 is true'


         ELSE IF (WORD.EQ.'NIHPAIRONLY') THEN

            NIHPAIRONLYT   = .TRUE.
            PRINT '(A)',' keywords> Ni plus hydrogen only two body interactions true'

! 
! NIH2LEPS turns ON/OFF the LEPS potential for Ni + H2 PES
! default OFF
! 
         ELSE IF (WORD.EQ.'NIH2LEPS') THEN

            NIH2LEPST   = .TRUE.
            PRINT '(A)',' keywords> Ni plus two hydrogen is true'


! 
! NIHLEPS turns ON/OFF a LEPS potential for Ni + H PES
! default OFF
! 
         ELSE IF (WORD.EQ.'NIHLEPS') THEN

            NIHLEPST   = .TRUE.
            PRINT '(A)',' keywords> Ni plus hydrogen LEPS is true'

! 
! CHARMM related keyword to reject transition states
! that connect two minima with different omega angles, i.e. to prevent cis-trans peptide
! isomerisation.
! 
         ELSE IF (WORD.EQ.'UACHIRAL') THEN
            UACHIRAL=.TRUE.

         ELSE IF (WORD.EQ.'NOCHIRALCHECKS') THEN
            TURNOFFCHECKCHIRALITY=.TRUE.

         ELSE IF (WORD.EQ.'NOCHIRALENDPOINTS') THEN
            CHIRALENDPOINTS=.FALSE.
            WRITE(*,*) ' keywords> No chirality checks for endpoints'

         ELSE IF (WORD.EQ.'NOCISTRANS') THEN
            NOCISTRANS=.TRUE.   ! is used in connect.f
            CHECKOMEGAT=.TRUE.  ! is used in NEWNEB
            IF (NITEMS.GT.1) CALL READF(MINOMEGA)

            IF (NITEMS.GT.2) THEN
               CALL READU(WW)
               IF (TRIM(ADJUSTL(WW))=='RNA') THEN
                  NOCISTRANSRNA = .TRUE.
                  write(*,*) ' keywords> NOCISTRANSRNA set to .TRUE.'
               ELSE IF (TRIM(ADJUSTL(WW))=='DNA') THEN
                  NOCISTRANSDNA = .TRUE.
                  write(*,*) ' keywords> NOCISTRANSDNA set to .TRUE.'
               ELSE IF (TRIM(ADJUSTL(WW))=='ALWAYS') THEN
                  CHECKCISTRANSALWAYS = .TRUE.
               ELSE IF (TRIM(ADJUSTL(WW))=='ALWAYSRNA') THEN
                  CHECKCISTRANSALWAYSRNA = .TRUE.
               ELSE IF (TRIM(ADJUSTL(WW))=='ALWAYSDNA') THEN
                  CHECKCISTRANSALWAYSDNA = .TRUE.
               ELSE
                  WRITE(*,*) ' keywords> ERROR - currently no other nocistrans options implemented than for RNA and DNA'
               ENDIF
            ENDIF
! 
! No frequencies should be evaluated or placed in the path.info file.
! 
         ELSE IF (WORD.EQ.'NOFRQS') THEN
            NOFRQS=.TRUE.
! 
! No Hessian should be calculated during geometry optimisation.
! 
         ELSE IF (WORD.EQ.'NOHESS') THEN
            NOHESS=.TRUE.
! 
! If NOIT is true and we have a Hessian then use DSYEVR to calculate eigenvectors
! 
         ELSE IF (WORD.EQ.'NOINTNEWT') THEN
! dont use newtons method to converge internals back transform
            INTNEWT = .FALSE.

         ELSE IF (WORD.EQ.'NOIT') THEN
            NOIT=.TRUE.
! 
! Don't try inversion in minpermdist routines.
! Needed for two-dimensional Morse for example.
! 
         ELSE IF (WORD.EQ.'NOINVERSION') THEN
            NOINVERSION=.TRUE.
! 
! For the growing string or evolving string double-ended
! transition state search methods, instead of using L-BFGS optimization to
! evolve the strings, simply take steps in the direction of the perpendicular force.
! 

         ELSE IF (WORD.EQ.'NOLBFGS') THEN
            NOLBFGS = .TRUE.
         ELSE IF (WORD == 'NONEBMIND') THEN
            NEBMIND=.FALSE.
            PRINT *, 'keywords> Structures supplied to NEB will NOT be put in the closest coincidence'
! 
! NONLOCAL x y z factors for averaged Gaussian, Morse type 1 and Morse
! type 2 potentials to include                  - default 0 0 0
! 
         ELSE IF (WORD.EQ.'NONLOCAL') THEN
            CALL READF(GFRACTION)
            IF (NITEMS.GT.2) THEN
               CALL READF(MFRACTION1)
            ENDIF
            IF (NITEMS.GT.3) THEN
               CALL READF(MFRACTION2)
            ENDIF
            FTEST=.TRUE.
            IF ((GFRACTION.EQ.0.0D0).AND.(MFRACTION1.EQ.0.0D0).AND.
     1      (MFRACTION2.EQ.0.0D0)) FTEST=.FALSE.

         ELSE IF (WORD.EQ.'NOPERMPROCHIRAL') THEN
            NOPERMPROCHIRAL = .TRUE.

! 
! Reduce printing of coordinates.
! 
         ELSE IF (WORD.EQ.'NOPOINTS') THEN
            PRINTPTS=.FALSE.
! 
! Used in CHARMM transition state guessing procedure
! together with TWISTTYPE. Setting randomcutoff very large prevents random
! steps, and is recommended.

! 
         ELSE IF (WORD.EQ.'NORANDOM') THEN
            NORANDOM=.TRUE.
            IF (NITEMS.GT.1) CALL READF(RANDOMCUTOFF)
! 
! Whether to put periodic images back in the primary supercell.
! 
         ELSE IF (WORD .EQ. 'NORESET') THEN
            NORESET=.TRUE.

         ELSE IF (WORD .EQ. 'NORMALMODE') THEN
            NORMALMODET=.TRUE.

         ELSE IF (WORD.EQ.'NTIP') THEN
            CALL READI(TIPID)
            IF (TIPID /= 4) THEN
               PRINT *, 'NOT YET INCLUDED'
               STOP
            ENDIF
            NTIPT = .TRUE.
            RBAAT = .TRUE.
            NRBSITES = 4
            ALLOCATE(RBSITE(NRBSITES,3))
            ALLOCATE(STCHRG(NRBSITES))
            NTSITES = NATOMS*NRBSITES/2
            IF (PERMDIST) THEN ! correct all permutations allowed if perm.allow is not given explicitly
               IF (NPERMSIZE(1).EQ.NATOMS) NPERMSIZE(1)=NATOMS/2
            ENDIF

            FRQCONV = 53.0883746D0
            WRITE(*,*) "keywords> Frequencies (and square frequencies) will be given in cm^-1 (cm^-2)"
            WRITE(*,*) "keywords> FRQCONV = ", FRQCONV
! 
! OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
! 
         ELSE IF (WORD.EQ.'ODIHE') THEN
            ODIHET=.TRUE.
            WRITE(*,'(A)') 'ODIHE set: dihedral-angle order parameter will be calculated'
            WRITE(*,'(A)') 'using the reference structure supplied in ref.crd'
! 
! Add an octahedral field to the potential of magnitude FOH.
! 
         ELSE IF (WORD.EQ.'OH') THEN
            FIELDT=.TRUE.
            OHT=.TRUE.
            CALL READF(FOH)
! 
! Specify Oh supercell to allow box symmetries in permutational alignment.
! 
         ELSE IF (WORD.EQ.'OHCELL') THEN
            OHCELLT=.TRUE.
            BULKT=.TRUE. ! just in case we forget in odata!
            WRITE(*,'(A)') ' Octahedral supercell specfied'
         ELSE IF (WORD.EQ.'ODIHE') THEN
            ODIHET=.TRUE.
            WRITE(*,'(A)') 'ODIHE set: dihedral-angle order parameter will be calculated'
            WRITE(*,'(A)') 'using the reference structure supplied in ref.crd'

! Orbital localisation potential.
         ELSE IF (WORD.EQ.'ORBITALS') THEN
            ORBITALS = .TRUE.
            CALL READI(ORBVAREXPONENT)
            PRINT *, ORBVAREXPONENT
            CALL ORBITALS_INIT(NORBS, NROTS)
C
C  dm368 potential
C
      ELSE IF (WORD.EQ.'EX1D') THEN
         EX1DT=.TRUE.
         IF (NITEMS.GT.1) CALL READI(NATOMS)
!         IF (MOD(NONEDAPBC,3).NE.0) THEN
!            WRITE(*,'(A)') 'keywords> ERROR *** lattice dimension must be a multiple of three'
!            STOP
!         ENDIF


         ELSE IF (WORD.EQ.'OLDINTMINPERM') THEN
            INTMINPERMT = .TRUE.
            OLDINTMINPERMT=.TRUE.
! 
! Specify 1D APBC potential
! 
         ELSE IF (WORD.EQ.'ONEDAPBC') THEN
            ONEDAPBCT=.TRUE.
            NONEDAPBC=NATOMS
            PRINT *,'NONEDAPBC=',NONEDAPBC
            IF (NITEMS.GT.1) CALL READI(NONEDAPBC)
! IF (MOD(NONEDAPBC,3).NE.0) THEN
! WRITE(*,'(A)') 'keywords> ERROR *** lattice dimension must be a multiple of three'
! STOP
! ENDIF
            ALLOCATE(XYPHI(NONEDAPBC))
            LUNIT=GETUNIT()
            OPEN (LUNIT,FILE='ONED.phi',STATUS='OLD')
            READ(LUNIT,*) XYPHI(1:NONEDAPBC)
            CLOSE(LUNIT)
            WRITE(*,'(A)') ' keywords> phi values'
            WRITE(*,'(3G20.10)') XYPHI(1:NONEDAPBC)
! 
! Specify 1D PBC potential
! 
         ELSE IF (WORD.EQ.'ONEDPBC') THEN
            ONEDPBCT=.TRUE.
            NONEDAPBC=NATOMS
            IF (NITEMS.GT.1) CALL READI(NONEDAPBC)
! IF (MOD(NONEDAPBC,3).NE.0) THEN
! WRITE(*,'(A)') 'keywords> ERROR *** lattice dimension must be a multiple of three'
! STOP
! ENDIF
            ALLOCATE(XYPHI(NONEDAPBC))
            LUNIT=GETUNIT()
            OPEN (LUNIT,FILE='ONED.phi',STATUS='OLD')
            READ(LUNIT,*) XYPHI(1:NONEDAPBC)
            CLOSE(LUNIT)
            WRITE(*,'(A)') ' keywords> phi values'
            WRITE(*,'(3G20.10)') XYPHI(1:NONEDAPBC)
         ELSE IF (WORD.EQ.'INVTONEDPBC') THEN
            INVTONEDPBCT=.TRUE.
            NONEDAPBC=NATOMS
            IF (NITEMS.GT.1) CALL READI(NONEDAPBC)
            ALLOCATE(XYPHI(NONEDAPBC))
            LUNIT=GETUNIT()
            OPEN (LUNIT,FILE='ONED.phi',STATUS='OLD')
            READ(LUNIT,*) XYPHI(1:NONEDAPBC)
            CLOSE(LUNIT)
            WRITE(*,'(A)') ' keywords> phi values'
            WRITE(*,'(3G20.10)') XYPHI(1:NONEDAPBC)
         ELSE IF (WORD.EQ.'INVTTWODPBC') THEN
            INVTTWODPBCT=.TRUE.
            NONEDAPBC=NATOMS
            IF (NITEMS.GT.1) CALL READI(NONEDAPBC)
            ALLOCATE(XYPHI(2*(NONEDAPBC**2)))
            LUNIT=GETUNIT()
            OPEN (LUNIT,FILE='TWOD.phi',STATUS='OLD')
            READ(LUNIT,*) XYPHI(1:2*(NONEDAPBC**2))
            CLOSE(LUNIT)
            WRITE(*,'(A)') ' keywords> phi values'
            WRITE(*,'(3G20.10)') XYPHI(1:2*(NONEDAPBC**2))
! 
! Specify 2D APBC potential
! 
         ELSE IF (WORD.EQ.'TWODAPBC') THEN
            TWODAPBCT=.TRUE.
            IF (NITEMS.GT.1) CALL READI(NONEDAPBC)
            ALLOCATE(XYPHI(2*(NONEDAPBC**2)))
            LUNIT=GETUNIT()
            OPEN (LUNIT,FILE='TWOD.phi',STATUS='OLD')
            READ(LUNIT,*) XYPHI(1:2*(NONEDAPBC**2))
            CLOSE(LUNIT)
            WRITE(*,'(A,I5,I5)') ' keywords> phi values, N_lattice, NATOMS', NONEDAPBC, NATOMS
            WRITE(*,'(3G20.10)') XYPHI(1:2*(NONEDAPBC**2))
! 
! Specify 2D PBC potential
! 
         ELSE IF (WORD.EQ.'TWODPBC') THEN
            TWODPBCT=.TRUE.
            IF (NITEMS.GT.1) CALL READI(NONEDAPBC)
! IF (MOD(NONEDAPBC,3).NE.0) THEN
! WRITE(*,'(A)') 'keywords> ERROR *** lattice dimension must be a multiple of three'
! STOP
! ENDIF
            ALLOCATE(XYPHI(2*(NONEDAPBC**2)))
            LUNIT=GETUNIT()
            OPEN (LUNIT,FILE='TWOD.phi',STATUS='OLD')
            READ(LUNIT,*) XYPHI(1:2*(NONEDAPBC**2))
            CLOSE(LUNIT)
            WRITE(*,'(A)') ' keywords> phi values'
            WRITE(*,'(3G20.10)') XYPHI(1:2*(NONEDAPBC**2))

! 
! Specify 3D PBC potential
! 
         ELSE IF (WORD.EQ.'THREEDPBC') THEN
            THREEDPBCT=.TRUE.
            IF (NITEMS.GT.1) CALL READI(NONEDAPBC)
! IF (MOD(NONEDAPBC,3).NE.0) THEN
! WRITE(*,'(A)') 'keywords> ERROR *** lattice dimension must be a multiple of three'
! STOP
! ENDIF
            ALLOCATE(XYPHI(3*(NONEDAPBC**3)))
            LUNIT=GETUNIT()
            OPEN (LUNIT,FILE='THREED.phi',STATUS='OLD')
            READ(LUNIT,*) XYPHI(1:3*(NONEDAPBC**3))
            CLOSE(LUNIT)
            WRITE(*,'(A)') ' keywords> phi values'
            WRITE(*,'(3G20.10)') XYPHI(1:3*(NONEDAPBC**3))
! 
! Specify 3D APBC potential
! 
         ELSE IF (WORD.EQ.'THREEDAPBC') THEN
            THREEDAPBCT=.TRUE.
            IF (NITEMS.GT.1) CALL READI(NONEDAPBC)
! IF (MOD(NONEDAPBC,3).NE.0) THEN
! WRITE(*,'(A)') 'keywords> ERROR *** lattice dimension must be a multiple of three'
! STOP
! ENDIF
            ALLOCATE(XYPHI(3*(NONEDAPBC**3)))
            LUNIT=GETUNIT()
            OPEN (LUNIT,FILE='THREED.phi',STATUS='OLD')
            READ(LUNIT,*) XYPHI(1:3*(NONEDAPBC**3))
            CLOSE(LUNIT)
            WRITE(*,'(A)') ' keywords> phi values'
            WRITE(*,'(3G20.10)') XYPHI(1:3*(NONEDAPBC**3))


! 
! Specify 4D APBC potential
! 
         ELSE IF (WORD.EQ.'FOURDAPBC') THEN
            FOURDAPBCT=.TRUE.
            IF (NITEMS.GT.1) CALL READI(NONEDAPBC)
! IF (MOD(NONEDAPBC,3).NE.0) THEN
! WRITE(*,'(A)') 'keywords> ERROR *** lattice dimension must be a multiple of three'
! STOP
! ENDIF
            ALLOCATE(XYPHI(4*(NONEDAPBC**4)))
            LUNIT=GETUNIT()
            OPEN (LUNIT,FILE='FOURD.phi',STATUS='OLD')
            READ(LUNIT,*) XYPHI(1:4*(NONEDAPBC**4))
            CLOSE(LUNIT)
            WRITE(*,'(A)') ' keywords> phi values'
            WRITE(*,'(3G20.10)') XYPHI(1:4*(NONEDAPBC**4))
! 
! Specify 4D PBC potential
! 
         ELSE IF (WORD.EQ.'FOURDPBC') THEN
            FOURDPBCT=.TRUE.
            IF (NITEMS.GT.1) CALL READI(NONEDAPBC)
! IF (MOD(NONEDAPBC,3).NE.0) THEN
! WRITE(*,'(A)') 'keywords> ERROR *** lattice dimension must be a multiple of three'
! STOP
! ENDIF
            ALLOCATE(XYPHI(4*(NONEDAPBC**4)))
            LUNIT=GETUNIT()
            OPEN (LUNIT,FILE='FOURD.phi',STATUS='OLD')
            READ(LUNIT,*) XYPHI(1:4*(NONEDAPBC**4))
            CLOSE(LUNIT)
            WRITE(*,'(A)') ' keywords> phi values'
            WRITE(*,'(3G20.10)') XYPHI(1:4*(NONEDAPBC**4))

! 
! ONETEP tells the program to read derivative information in
! ONETEP format.                                        - default FALSE
! 
         ELSE IF ((WORD.EQ.'ONETEP').OR.(WORD.EQ.'ONETEPC')) THEN
            ONETEP=.TRUE.
            IF (WORD.EQ.'ONETEP') DFTP=.TRUE.
            IF (NITEMS.GT.2) THEN
               CALL READA(ONETEPJOB)
               CALL READA(SYS)
               ONETEPJOB=TRIM(ADJUSTL(ONETEPJOB)) // ' ' // TRIM(ADJUSTL(SYS)) // ' >& ' // TRIM(ADJUSTL(SYS)) // '.onetep'
            ELSE
               WRITE(*,'(A)') 'keywords> ERROR - ONETEP job or system unspecified'
               CALL FLUSH(6)
               STOP
            ENDIF
            DO J1=1,80
               IF (SYS(J1:J1).EQ.' ') THEN
                  LSYS=J1-1
                  GOTO 24
               ENDIF
            ENDDO
24          CONTINUE

            ! sn402
            WRITE(*,*) "keywords> WARNING: there is currently no default frequency conversion set for ONETEP"
            WRITE(*,*) "keywords> Log products of frequencies will be computed in internal units."
            WRITE(*,*) "keywords> To learn how to set a default conversion factor, check the comments for 
     &                            the FRQCONV keyword in keywords.f"



         ELSE IF (WORD.EQ.'OPEP') THEN
            OPEPT=.TRUE.
            CALL READA(OPEP_DUMMY)
            IF (OPEP_DUMMY.EQ.'RNA') THEN
               OPEP_RNAT = .TRUE.
               WRITE(*,'(A)') 'keyword> RNA simulation using OPEP'
            ELSE
               WRITE(*,'(A)') 'keyword> Protein simulation using OPEP'
            ENDIF
            CALL READA(OPEP_DUMMY2)
            IF (.NOT.ALLOCATED(ATMASS)) ALLOCATE(ATMASS(NATOMS))
            CALL OPEP_INIT(NATOMS,Q(1:3*NATOMS),ATMASS(1:NATOMS),OPEP_RNAT)
            
            IF (OPEP_DUMMY2.EQ.'start') THEN
               IF (FILTH2 .NE. 0) THEN
                  WRITE(OTEMP, *) FILTH2
                  WRITE(OSTRING,'(A)') TRIM(ADJUSTL(OPEP_DUMMY2))//'.'//TRIM(ADJUSTL(OTEMP))
                  WRITE(*,*) 'ostring=', OSTRING
               ELSE
                  WRITE(OSTRING,'(A)') TRIM(ADJUSTL(OPEP_DUMMY2))
               END IF
               WRITE(*,'(A)') ' keywords> input coordinates for OPEP system will be read from ', OSTRING
               OPEN(UNIT=3827, FILE=TRIM(ADJUSTL(OSTRING)),STATUS='UNKNOWN')
               WRITE(*,'(A)') ' keywords> reading in xyz format'
               DO I = 1, NATOMS
                  READ(3827, *) Q(3*I-2:3*I)
               END DO
               CLOSE(3827)
            ELSE
               WRITE(*,'(A)') ' keywords> reading in pdb format'              
            ENDIF

            !use same frequency conversion as AMBER does
            IF (DUMMY_FRQCONV .EQ. 0.0D0) THEN
               FRQCONV = 2.045483D13
               WRITE(*,*) "keywords> Set frequency conversion factor to the AMBER default value: ", FRQCONV
               WRITE(*,*) "keywords> This corresponds to frequencies being given in radians/s"
            ELSE
               FRQCONV = DUMMY_FRQCONV
               WRITE(*,*) "keywords> Set frequency conversion factor to the user-specified value: ", FRQCONV
            ENDIF
            FRQCONV2 = FRQCONV*FRQCONV

            RETURN

!  
! Optimise TS with SQVV
! 
         ELSE IF (WORD == 'OPTIMIZETS') THEN
            OPTIMIZETS=.TRUE.
! 
! Distance cutoff for distinguishing atoms in the same orbit for LPERMDIST and LOCALPERMDIST
! 
         ELSE IF (WORD.EQ.'ORBITGEOMTOL') THEN
            CALL READF(LPDGEOMDIFFTOL)
! 
! Calculates order parameters and theire derivatives wrt normal modes at the end of a geometry optimisation.
! The 1st argument is the number of order parameters to be calculated. The next arguments then specify
! the order parameters (defined by a 4 letters) and, if necessary, further information regarding this
! order parameter can be given. If such details are not required, set them to -9999.
! Following order parameters are currently supported: DIHEdral angles for CHARMM.
! 
         ELSE IF (WORD == 'ORDERPARAM') THEN
            ORDERPARAMT=.TRUE.
            CALL READI(NORDER)
            ALLOCATE(WHICHORDER(NORDER),ORDERNUM(NORDER))
            DO J1=1,NORDER
               CALL READA(WHICHORDER(J1))
               ORDERNUM(J1)=-9999
               CALL READI(ORDERNUM(J1))
            ENDDO
            IF (TRIM(ADJUSTL(WHICHORDER(1))).EQ.'DECA') THEN
               PRINT '(A)',' keyword> Assuming LJ75 order parameter calculation - looking for deca_and_icos file'
               LUNIT=GETUNIT()
               OPEN (LUNIT,FILE='deca_and_icos',STATUS='OLD')
               IF (.NOT.(ALLOCATED(POINTSDECA))) ALLOCATE(POINTSDECA(225),POINTSICOS(225))
               READ(LUNIT,*) POINTSDECA(1:225),POINTSICOS(1:225)
               CLOSE(LUNIT)
            ENDIF
! 
! Remove overall trans/rot with SQVV
! 
         ELSE IF (WORD == 'ORT') THEN
            ORT = .TRUE.
         ELSE IF (WORD.EQ.'OSASA') THEN
            OSASAT=.TRUE.
            CALL READF(RPRO)
            WRITE(*,'(A)') 'OSASA set: solvent accessible surface area order parameter will be calculated'
            WRITE(*,'(A,F3.1)') 'using probe radius ',RPRO
! 
! PPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPPP
! 
         ELSE IF (WORD .EQ. 'PAHA') THEN

            CALL READI(PAHID)

            IF (PAHID == 1) THEN
               NRBSITES = 12
            ELSEIF (PAHID == 2) THEN
               NRBSITES = 18
            ELSEIF (PAHID == 3) THEN
               NRBSITES = 24
            ELSEIF (PAHID == 4) THEN
               NRBSITES = 26
            ENDIF

            PAHAT    = .TRUE.
            RBAAT    = .TRUE.
            ALLOCATE(RBSITE(NRBSITES,3))
            ALLOCATE(RBSTLA(NRBSITES,3))
            ALLOCATE(STCHRG(NRBSITES))

            NTSITES = (NATOMS/2)*NRBSITES

            CALL DEFPAHA()

            IF (PAHID == 1) THEN
               NCARBON  = 6
               CALL DEFBENZENE()
            ELSEIF (PAHID == 2) THEN
               NCARBON  = 10
               CALL DEFNAPHTHALENE()
            ELSEIF (PAHID == 3) THEN
               NCARBON  = 14
               CALL DEFANTHRACENE()
            ELSEIF (PAHID == 4) THEN
               NCARBON  = 16
               CALL DEFPYRENE()
            ENDIF

         ELSE IF (WORD .EQ. 'PAP') THEN

            CALL READI(PAPID)
            CALL READF(PAPALP)
            CALL READF(PAPS)
            CALL READF(PAPCD)
            CALL READF(PAPEPS)

            IF (PAPID == 1) THEN
               NRBSITES = 7
            ELSEIF (PAPID == 2) THEN
               NRBSITES = 5
            ELSEIF (PAPID == 3) THEN
               NRBSITES = 3
            ELSEIF (PAPID == 4) THEN
               NRBSITES = 5
            ENDIF

            PAPT     = .TRUE.
            RBAAT    = .TRUE.
            NTSITES = (NATOMS/2)*NRBSITES

            CALL DEFPAP()

! 
! Turn on colouring according to pair energy in dump of path.xyz.
! The argument is the number of duplicate frames to add for the two endpoints.
! 
         ELSE IF (WORD.EQ.'PAIRCOLOUR') THEN
            PAIRCOLOURT=.TRUE.
            IF (NITEMS.GT.1) CALL READI(NENDDUP)

         ELSE IF (WORD.EQ.'PARALLEL') THEN
            PARALLEL=.TRUE.
            CALL READA(NPROC)
! 
! PARAMS n1 n2 ... up to seven real input parameters used for the
! following atom types:
! AX: Z*
! M:  rho
! MV: rho, delta
! ME: N, M, BOXLENGTHS X, Y, Z AND CUTOFF (N, M ARE READ DOUBLE PRECISION)
! JM: box lengths x, y, z and cutoff
! SC: box lengths x, y, z and cutoff (epsilon, c, sigma are read from SCparams)
! P6: box lengths x, y, z and cutoff
! AU: epsilon, c, sigma
! AG: epsilon, c, sigma
! NI: epsilon, c, sigma
! 
         ELSE IF (WORD.EQ.'PARAMS') THEN
            CALL READF(PARAM1)
            GALPHA=PARAM1
            MALPHA1=PARAM1
            MALPHA2=PARAM1
            IF (NITEMS.GT.2) THEN
               CALL READF(PARAM2)
            ENDIF
            IF (NITEMS.GT.3) THEN
               CALL READF(PARAM3)
            ENDIF
            IF (NITEMS.GT.4) THEN
               CALL READF(PARAM4)
            ENDIF
            IF (NITEMS.GT.5) THEN
               CALL READF(PARAM5)
            ENDIF
            IF (NITEMS.GT.6) THEN
               CALL READF(PARAM6)
            ENDIF
            IF (NITEMS.GT.7) THEN
               CALL READF(PARAM7)
            ENDIF

         ELSE IF (WORD.EQ.'PATCHYD') THEN

            PATCHYDT = .TRUE.
            RBAAT = .TRUE.
            CALL READI(NRBSITES)

            ALLOCATE(RBSITE(NRBSITES,3))
            ALLOCATE(RBSTLA(NRBSITES,3))

            CALL DEFPATCHES

            NTSITES = NATOMS*NRBSITES/2

! 
! PATH specifies calculation of the pathway connecting two minima from the transition
! state specified in odata. NPATHFRAME is the number of points files to save on either
! side. A complete xyz file is printed to path.xyz and the energy as a function of
! path length is printed to file EofS.
! Movies generated in this way tend to move too fast for the interesting bits, and too
! slow around stationary points. Specify FRAMEEDIFF to give a lower bound to the energy difference
! between frames for which the structure is considered different.
! 
         ELSE IF (WORD.EQ.'PATH') THEN
            PATHT=.TRUE.
            IF (NITEMS.GT.1) THEN
               CALL READI(NPATHFRAME)
               ! if (NPATHFRAME<3) THEN
               ! PRINT *, 'Number of path frames cannot be less than 3 - stop'
               ! stop
               ! ELSE IF (NPATHFRAME>3) THEN
               ! IF (.NOT.PRINTPTS) THEN
               ! PRINT *, 'Number of path frames is more than 3 - dumping all points!'
               ! PRINTPTS=.TRUE.
               ! ENDIF
               ! ENDIF
            ENDIF
            IF (NITEMS.GT.2) CALL READF(FRAMEEDIFF)
            IF (NITEMS.GT.3) CALL READF(FRAMESDIFF)
            IF (NPATHFRAME.LE.0) PRINTPTS=.FALSE.
            IF (NPATHFRAME.GT.0) PRINTPTS=.TRUE.

         ELSE IF (WORD .EQ. 'PTSTST') THEN

            CALL READI(PAPID)
            CALL READF(PAPEPS)
            CALL READF(PAPCD)

            IF (PAPID == 1) THEN
               NRBSITES = 7
            ELSEIF (PAPID == 2) THEN
               NRBSITES = 5
            ELSEIF (PAPID == 3) THEN
               NRBSITES = 3
            ELSEIF (PAPID == 4) THEN
               NRBSITES = 5
            ENDIF

            PTSTSTT  = .TRUE.
            RBAAT    = .TRUE.
            ALLOCATE(RBSITE(NRBSITES,3))
            NTSITES = NATOMS*NRBSITES/2
            CALL DEFPTSTST()

         ELSE IF (WORD.EQ.'PERMDIHE') THEN
! 
! PATHSDSTEPS sets the number of SD steps allowed at the beginning of a path
! calculation. We switch to LBFGS from RKMIN, BSMIN and SEARCH INR methods if
! they don't converge in PATHSDSTEPS steps. If not set then the default is NSTEPS.
! 
         ELSE IF (WORD.EQ.'PATHSDSTEPS') THEN
            CALL READI(PATHSDSTEPS)
         ELSE IF (WORD.EQ.'PERMDIHE') THEN
            PERMDIHET=.TRUE.
            DO J1=1,NITEMS-1
               CALL READI(NDUM)
               PERMDIHE(J1)=NDUM
            ENDDO
            NPERMDIHE=NITEMS-1
            DO J1=1,NITEMS-1
               PRINT *,'PERMDIHE',PERMDIHE(J1)
            ENDDO
! 
! Whether to optimise the permutational isomers in assessing optimal
! alignment.
! 
         ELSE IF ((WORD.EQ.'PERMDIST').OR.(WORD.EQ.'PERMDISTINIT').OR.(WORD.EQ.'LOCALPERMDIST').OR.
     &   (WORD.EQ.'LPERMDIST')) THEN
            PERMDIST=.TRUE.
            IF (WORD.EQ.'PERMDISTINIT') PERMDISTINIT=.TRUE.
            IF (WORD.EQ.'LOCALPERMDIST') THEN
               LOCALPERMDIST=.TRUE.
               IF (NITEMS.GT.1) CALL READF(LPDGEOMDIFFTOL)
               IF (NITEMS.GT.2) CALL READF(RBCUTOFF)
               IF (NITEMS.GT.3) CALL READI(NRBTRIES)
               PRINT '(A)',' keyword> Local rigid body permutational alignment:'
               PRINT '(2(A,F12.4),A,I6)','          distance tolerance=',LPDGEOMDIFFTOL,' cutoff=',RBCUTOFF,
     &         ' number of passes through alignment phase=',NRBTRIES
            ELSEIF (WORD.EQ.'LPERMDIST') THEN
               LPERMDIST=.TRUE.
               IF (NITEMS.GT.1) CALL READI(LOCALPERMNEIGH)
               IF (NITEMS.GT.2) CALL READF(LOCALPERMCUT)
               IF (NITEMS.GT.3) CALL READF(LOCALPERMCUT2)
               IF (NITEMS.GT.4) CALL READF(ORBITTOL)
               ! IF (NITEMS.GT.3) CALL READI(LOCALPERMMAXSEP)
               PRINT '(A,F15.5)',' keyword> Local permutational alignment: alignment threshold=',LOCALPERMCUT
               PRINT '(A,F15.5)',' keyword> Local permutational alignment: alignment cutoff=   ',LOCALPERMCUT2
               PRINT '(A,F15.5)',' keyword> Distance tolerance for distinguishing atoms in the same orbit=',ORBITTOL
            ELSEIF (WORD.EQ.'PERMDIST') THEN
               IF (NITEMS.GT.1) CALL READF(ORBITTOL)
               IF (NITEMS.GT.2) CALL READI(MAXNSETS)
               PRINT '(A,F15.5)',' keyword> Distance tolerance for distinguishing atoms in the same orbit=',ORBITTOL
               PRINT '(A,I3)',' keyword> Maximum number of secondary sets in perm.allow file=', MAXNSETS
            ENDIF

            INQUIRE(FILE='perm.allow',EXIST=PERMFILE)
            IF (PERMFILE) THEN
               OPEN(UNIT=1,FILE='perm.allow',STATUS='OLD')
               READ(1,*) NPERMGROUP
               IF (ALLOCATED(PERMGROUP)) THEN
                  PRINT '(A)',' keywords> ERROR *** array PERMGROUP already allocated. Incompatible keywords?'
                  STOP
               ENDIF
               ! ALLOCATE(NPERMSIZE(NATOMS),PERMGROUP(NATOMS),NSWAP(NATOMS),SWAP1(NATOMS,3),SWAP2(NATOMS,3))
               ALLOCATE(NPERMSIZE(3*NATOMS),PERMGROUP(3*NATOMS),NSETS(3*NATOMS),SETS(NATOMS,MAXNSETS))
               NPERMSIZE(:) = 0
               PERMGROUP(:) = 0
               NSETS(:) = 0
               SETS(:,:) = 0
               ! 
               ! The above dimensions were fixed at NATOMS because:
               ! (a) Atoms were not allowed to appear in more than one group.
               ! (b) The maximum number of pair exchanges associated with a group is three.
               ! 
               ! However, for flexible water models we need to exchange all waters,
               ! and we can exchange H's attached to the same O. The dimension required
               ! becomes 3*NATOMS
               ! 

               NDUMMY = 1
               DO J1=1,NPERMGROUP
                  READ(1,*) NPERMSIZE(J1),NSETS(J1)
                  ! 
                  ! Sanity checks!
                  ! 
                  IF (NSETS(J1).GT.MAXNSETS) THEN
                     PRINT '(2(A,I8))','keyword> ERROR - number of secondary sets ',NSETS(J1), ' is > ', MAXNSETS
                     STOP
                  ENDIF
                  ! IF (NDUMMY+NPERMSIZE(J1)-1.GT.NATOMS) THEN
                  IF (NDUMMY+NPERMSIZE(J1).GT.3*NATOMS) THEN
                     PRINT '(2(A,I8))','keyword> ERROR - number of atoms to be permuted in all groups is > 3*number of atoms'
                     STOP
                  ENDIF
                  ! READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SETS(PERMGROUP(J3),J2),J2=1,NSETS(J1)),
                  ! &                                                            J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1)
                  READ(1,*) PERMGROUP(NDUMMY:NDUMMY+NPERMSIZE(J1)-1),((SETS(PERMGROUP(J3),J2),J3=NDUMMY,NDUMMY+NPERMSIZE(J1)-1),
     &            J2=1,NSETS(J1))
                  NDUMMY=NDUMMY+NPERMSIZE(J1)
               ENDDO
               ! 
               ! And another sanity check! This condition is now allowed.
               ! 
               ! DO J1=1,NDUMMY
               ! DO J2=J1+1,NDUMMY
               ! IF (PERMGROUP(J2).EQ.PERMGROUP(J1)) THEN
               ! PRINT '(2(A,I8))','keyword> ERROR - atom ',PERMGROUP(J1),' appears more than once'
               ! STOP
               ! ENDIF
               ! ENDDO
               ! ENDDO
               CLOSE(1)

               ! 
               ! And yet another!
               ! 
               IF (NFREEZE.GT.0) THEN
                  NDUMMY=0
                  DO J1=1,NPERMGROUP
                     DO J2=1,NPERMSIZE(J1)
                        IF (FROZEN(PERMGROUP(NDUMMY+J2))) THEN
                           PRINT '(A,I8,A)',' keyword> ERROR atom ',PERMGROUP(NDUMMY+J2),' cannot be frozen and permuted'
                           STOP
                        ENDIF
                     ENDDO
                     NDUMMY=NDUMMY+NPERMSIZE(J1)
                  ENDDO
               ENDIF
            ELSE
               IF (ALLOCATED(PERMGROUP)) THEN
                  PRINT '(A)',' keywords> ERROR *** array PERMGROUP already allocated. Incompatible keywords?'
                  STOP
               ENDIF
               ALLOCATE(NPERMSIZE(NATOMS),PERMGROUP(NATOMS),NSETS(NATOMS),SETS(NATOMS,2))
               NSETS(1:NATOMS)=0
               NPERMGROUP=1 ! all atoms can be permuted - default
               NPERMSIZE(1)=NATOMS ! all atoms can be permuted - default
               IF (RBAAT) NPERMSIZE(1)=NATOMS/2 ! for rigid bodies
               DO J1=1,NPERMSIZE(1)
                  PERMGROUP(J1)=J1
               ENDDO
            ENDIF
            PRINT '(A,I6)',' keyword> Number of groups of permutable atoms=',NPERMGROUP
            NDUMMY=1
            IF (DEBUG) THEN
               DO J1=1,NPERMGROUP
                  PRINT '(A,3(I6,A))',' keyword> group ',J1,' contains ',NPERMSIZE(J1),' atoms with ',
     &            NSETS(J1),' additional atom sets:'
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
            ENDIF
! hk286
            IF (RIGIDINIT .AND. .NOT. RIGIDMOLECULEST) THEN
               CALL GENRIGID_PERMDIST()
            ENDIF

! 
! CHARMM and UNRES dihedral angle perturbation specification.
! Performs random, {\tt GMIN}-style twists before starting optimisation.
! 

         ELSE IF (WORD.EQ.'PERTDIHE') THEN
            PERTDIHET=.TRUE.
            CALL READF(CHPMAX)
            CHPMIN=CHPMAX
            CALL READF(CHNMIN)
            CALL READF(CHNMAX)
            CALL READI(ISEED)
! PRINT *,'CHPMIN,CHPMAX,CHNMIN,CHNMAX',CHPMIN,CHPMAX,CHNMIN,CHNMAX
! 
! Page-McIver minimisation for paths. This applies when the odata file specifies a
! ts search and the PATH keyword is also present.
! 
! 
! PHI4MOD for mean field phi^4 model
! 
         ELSE IF (WORD.EQ.'PHI4MOD') THEN
            PHI4MODT   = .TRUE.
            IF (NITEMS.GT.1) CALL READF(JPARAM)
            PRINT '(A,G20.10)',' keywords> PHI4 mean field model with J=',JPARAM
         ELSE IF (WORD.EQ.'PMPATH') THEN
            PMPATHT=.TRUE.
            IF (NITEMS.GT.1) CALL READI(PMPATHINR)
! 
! bf269> (harmonic = 2nd power) polymer potential
! 
         ELSE IF (WORD == 'POLY2') THEN
            HARMPOLYT=.TRUE.
            CALL READF(HARMPOLY_BONLEN)
            CALL READF(HARMPOLY_K)
! 
! bf269> pore (8th power in distance from axis constraining) potential 
! 
         ELSE IF (WORD == 'PORE8') THEN
            PORE8T=.TRUE.
            CALL READI(PORE8_AXIS)
            CALL READF(PORE8_ENERGY)
! 
! SQVV keyword.
! 
         ELSE IF (WORD == 'PRINTOPTIMIZETS') THEN
            PRINTOPTIMIZETS=.TRUE.
! 
! For the GS and ES double-ended transition state
! search methods, if using {\it FIXATMS\/} to zero some coordinates of the
! forces to avoid overall translation and rotation, this keyword will rotate
! the start and end points so that those coordinates are zero in both.
! 
         ELSE IF (WORD.EQ.'PREROTATE') THEN
            PREROTATE = .TRUE.
! 
! PRESSURE tells the program to perform a constant pressure optimisation
! for SC, ME and P6 with periodic boundary conditions - default off
! 
         ELSE IF (WORD.EQ.'PRESSURE') THEN
            PRESSURE=.TRUE.
! 
! PRINT n sets the value of IPRNT                              - default n=0
! 
         ELSE IF (WORD.EQ.'PRINT') THEN
            CALL READI(IPRNT)
! 
! Print ground state coefficients - only valid for MSEVB potential
! 
         ELSE IF (WORD.EQ.'PRINTCOEFFICIENTS') THEN
            printCoefficients=.TRUE.

! print out info on coordinates and stop; for debugging internals
         ELSE IF (WORD.EQ.'PRINTCOORDS') THEN
            PRINTCOORDS = .TRUE.
! 
! Keyword for applied static force.
! 
         ELSE IF (WORD.EQ.'PULL') THEN
            PULLT=.TRUE.
            CALL READI(PATOM1)
            CALL READI(PATOM2)
            CALL READF(PFORCE)
            IF (PFORCE.EQ.0.0D0) THEN
               WRITE(*,'(A,I6,A,I6,A,G20.10)') ' keyword> WARNING *** Pulling force is zero, turning off pulling directive'
               PULLT=.FALSE.
            ELSE
               WRITE(*,'(A,I6,A,I6,A,G20.10)') ' keyword> Pulling atoms ',PATOM1,' and ',PATOM2,' force=',PFORCE
            ENDIF
! 
! PUSHCUT sets the threshold for when a PUSHOFF will be applied, i.e.
! the RMS force must be less than PUSHCUT.
! 
         ELSE IF (WORD .EQ. 'PUSHCUT') THEN
            CALL READF(PUSHCUT)
! 
! PUSHOFF x sets the magnitude of the step away from a converged
! transition state if detected on the first cycle of
! a minimisation                                     - default x=0.01
! 
         ELSE IF (WORD .EQ. 'PUSHOFF') THEN
            CALL READF(PUSHOFF)
! 
! PUSHOPT x  specifies a golden section search for the best pushoff value for the + and -
! sides of a ts. The search will be bracketed in the range [0:x] and [0:-x].
! defaults x=0.01, pushoffconv 1.0D-4, pushoptmax 100
! 
         ELSE IF (WORD .EQ. 'PUSHOPT') THEN
            PUSHOPTT=.TRUE.
!           REDOTS=0 ! allow different pushoffs - changes initial bracket region
            IF (NITEMS.GT.1) CALL READF(PUSHOFF)
            IF (NITEMS.GT.2) CALL READF(PUSHOPTCONV)
            IF (NITEMS.GT.3) CALL READI(PUSHOPTMAX)
! 
! PV
! 
         ELSE IF (WORD.EQ.'PV') THEN
            PV=.TRUE.
            IF (NITEMS.GT.1) CALL READF(PRESS)
            IF (NITEMS.GT.2) CALL READF(PVCONV)
            IF (NITEMS.GT.3) CALL READF(PVTOL)
            IF (NITEMS.GT.4) CALL READI(PVSTEPS)
         ELSE IF (WORD.EQ.'PVTS') THEN
            PV=.TRUE.
            PVTS=.TRUE.
            NBOXTS=1
            IF (NITEMS.GT.1) CALL READF(PRESS)
            IF (NITEMS.GT.2) CALL READF(PVCONV)
            IF (NITEMS.GT.3) CALL READF(PVTOL)
            IF (NITEMS.GT.4) CALL READI(PVSTEPS)
            IF (NITEMS.GT.5) CALL READI(NBOXTS)
            WRITE(*,'(A,I5)') ' Searching uphill for a transition state in box length coordinate ',NBOXTS

         ELSE IF (WORD.EQ.'PYG') THEN
            NRBSITES = 1
            ALLOCATE(RBSITE(NRBSITES,3))
            PYGT  = .TRUE.
            RBAAT = .TRUE.
            CALL READF(PYA1(1))
            CALL READF(PYA1(2))
            CALL READF(PYA1(3))
            CALL READF(PYA2(1))
            CALL READF(PYA2(2))
            CALL READF(PYA2(3))
            CALL READF(PYSIGNOT)
            CALL READF(PYEPSNOT)

            IF (PYA1(1) == PYA2(1) .AND. PYA1(2) == PYA2(2) .AND. PYA1(3) == PYA2(3)) THEN
               RADIFT = .FALSE.
            ELSE
               RADIFT = .TRUE.
            ENDIF
! sf344> PY potential and extra LJ site
         ELSE IF (WORD.EQ.'PYBINARY') THEN
            NRBSITES = 2
            ALLOCATE(RBSITE(NRBSITES,3))
            PYBINARYT=.TRUE.
            ANGLEAXIS2=.TRUE.
            RBAAT=.TRUE.
            RADIFT=.TRUE.
            CALL READI(PYBINARYTYPE1)
            CALL READF(PYA11(1))
            CALL READF(PYA11(2))
            CALL READF(PYA11(3))
            CALL READF(PYA21(1))
            CALL READF(PYA21(2))
            CALL READF(PYA21(3))
            CALL READF(PYA12(1))
            CALL READF(PYA12(2))
            CALL READF(PYA12(3))
            CALL READF(PYA22(1))
            CALL READF(PYA22(2))
            CALL READF(PYA22(3))
            CALL READF(PYSIGNOT)
            CALL READF(PYEPSNOT)
            IF(.NOT.ALLOCATED(PYA1bin)) ALLOCATE(PYA1bin(NATOMS/2,3))
            IF(.NOT.ALLOCATED(PYA2bin)) ALLOCATE(PYA2bin(NATOMS/2,3))
            DO J1=1,NATOMS/2
               IF(J1<=PYBINARYTYPE1) THEN
                  PYA1bin(J1,:)=PYA11(:)
                  PYA2bin(J1,:)=PYA21(:)
               ELSE
                  PYA1bin(J1,:)=PYA12(:)
                  PYA2bin(J1,:)=PYA22(:)
               END IF
            END DO
            RBSITE(1,1)=PYA11(1)
            RBSITE(2,1)=-PYA11(1)
            RBSITE(:,2:3)=0.0D0
      ELSE IF (WORD .EQ. 'PY') THEN
         ! Syntax: PY sig_0 eps_0 [cut] [XYZ boxx boxy boxz]

         PYT = .TRUE.
         CALL READF(PYSIGNOT)
         CALL READF(PYEPSNOT)

         RBAAT = .TRUE.
         ! Rigid body SITE, NRBSITES, NTSITES information specified in py_input routine

         ! Specify cutoff for potential in absolute units
         IF (NITEMS.GT.3) THEN
            CALL READF(PCUTOFF)
            PARAMONOVCUTOFF=.TRUE.
            WRITE(*,*) "multisitepy cutoff: ", PCUTOFF
         ENDIF
         ! Specify periodic boundary conditions (PBCs)
         IF (NITEMS.GT.4) THEN
            ! control which dimensions have periodic boundaries with a string 'XYZ', always put x before y before z.            
            ! eg ...  Xz 20 30  specifies PBC on X and Z directions.  The X box size will be 20, the Z box size 30
            CALL READA(PBC)
            BOXLX=0
            BOXLY=0
            BOXLZ=0
            IF (SCAN(PBC,'Xx').NE.0) THEN
                PARAMONOVPBCX=.TRUE.
                CALL READF(BOXLX)
                BOXLX = BOXLX*PCUTOFF
                WRITE(*,*) "PBC X:",BOXLX
            ENDIF
            IF (SCAN(PBC,'Yy').NE.0) THEN
                PARAMONOVPBCY=.TRUE.
                CALL READF(BOXLY)
                BOXLY = BOXLY*PCUTOFF
                WRITE(*,*) "PBC Y:",BOXLY
            ENDIF
            IF (SCAN(PBC,'Zz').NE.0) THEN
                PARAMONOVPBCZ=.TRUE.
                CALL READF(BOXLZ)
                BOXLZ = BOXLZ*PCUTOFF
                WRITE(*,*) "PBC Z:",BOXLZ
            ENDIF
         ENDIF
      ELSE IF (WORD .EQ. 'PYGRAVITY') THEN
           CALL READF(PYGRAVITYC1)
           CALL READF(PYGRAVITYC2)
           EFIELDT=.TRUE.
      ELSE IF (WORD.EQ.'PYOVERLAPTHRESH') THEN
            CALL READF(PYOVERLAPTHRESH)
            PYLOCALSTEP(:)=1.0D0
            WRITE(*,'(A,F8.3)') 'keywords> ellipsoids considered to overlap for an ECF value of ', PYOVERLAPTHRESH
            IF(NITEMS.GT.2) THEN
               CALL READF(PYLOCALSTEP(1))
               CALL READF(PYLOCALSTEP(2))
            END IF
      ELSE IF (WORD.EQ.'PYADD') THEN
         PYADDT=.TRUE.
         CALL READI(NADDTARGET)
         WRITE(*,'(A,I6)') 'keyword> Target cluster size is ',NADDTARGET
         IF (MOD(NATOMS/2,NADDTARGET).NE.0) THEN
            WRITE(*,'(A,I6)') 'keyword> ERROR, target cluster size is not a factor of the number of PY particles ',
     &         NATOMS/2
            STOP
         ENDIF
         LUNIT=GETUNIT()
         OPEN(LUNIT,FILE='epsilon',STATUS='OLD')
         IF (.NOT.ALLOCATED(PYADDEPS)) ALLOCATE(PYADDEPS(NADDTARGET,NADDTARGET))
         DO J1=1,NADDTARGET
            DO J2=1,NADDTARGET
               READ(LUNIT,*) PYADDEPS(J2,J1)
               WRITE(*,'(2I6,G20.10)') J1,J2,PYADDEPS(J2,J1)
            ENDDO
         ENDDO
         CLOSE(LUNIT)
      ELSE IF (WORD.EQ.'PYADD2') THEN
         PYADD2T=.TRUE.
         CALL READI(NADDTARGET)
         WRITE(*,'(A,I6)') 'keyword> Target cluster size is ',NADDTARGET
         IF (MOD(NATOMS/2,NADDTARGET).NE.0) THEN
            WRITE(*,'(A,I6)') 'keyword> ERROR, target cluster size is not a factor of the number of PY particles ',
     &          NATOMS/2
            STOP
         ENDIF
         LUNIT=GETUNIT()
         OPEN(LUNIT,FILE='epsilon',STATUS='OLD')
         IF (.NOT.ALLOCATED(PYADDREP)) ALLOCATE(PYADDREP(NADDTARGET,NADDTARGET))
         IF (.NOT.ALLOCATED(PYADDATT)) ALLOCATE(PYADDATT(NADDTARGET,NADDTARGET))
         DO J1=1,NADDTARGET
            DO J2=1,NADDTARGET
               READ(LUNIT,*) PYADDREP(J2,J1), PYADDATT(J2,J1)
               WRITE(*,'(2I6,2G20.10)') J1,J2,PYADDREP(J2,J1),PYADDATT(J2,J1)
            ENDDO
         ENDDO
         CLOSE(LUNIT)

      ELSE IF (WORD.EQ.'PYGPERIODIC') THEN
            NRBSITES = 2
            NPYSITE=NATOMS/2
            ALLOCATE(RBSITE(NRBSITES,3))
            PYGPERIODICT = .TRUE.
!            ANGLEAXIS2=.TRUE.
            RBAAT=.TRUE.
            CALL READF(PYA1(1))
            CALL READF(PYA1(2))
            CALL READF(PYA1(3))
            CALL READF(PYA2(1))
            CALL READF(PYA2(2))
            CALL READF(PYA2(3))
            CALL READF(PYSIGNOT)
            CALL READF(PYEPSNOT)
            RBSITE(1,1)=PYA1(1)
            RBSITE(2,1)=0.0D0
            RBSITE(:,2:3)=0.0D0
            NTSITES=NATOMS*NRBSITES/2
            IF(.NOT.ALLOCATED(PYA1bin)) ALLOCATE(PYA1bin(NATOMS/2,3))
            IF(.NOT.ALLOCATED(PYA2bin)) ALLOCATE(PYA2bin(NATOMS/2,3))
            DO J1=1,NATOMS/2
               PYA1bin(J1,:)=PYA1(:)
               PYA2bin(J1,:)=PYA2(:)
            END DO
            IF (PYA1(1) == PYA2(1) .AND. PYA1(2) == PYA2(2) .AND. PYA1(3) == PYA2(3)) THEN
               RADIFT = .FALSE.
            ELSE
               RADIFT = .TRUE.
            ENDIF

            IF (NITEMS.GT.9) THEN
               CALL READF(PCUTOFF)
               PARAMONOVCUTOFF=.TRUE.
               PCUTOFF=PCUTOFF*PYSIGNOT
               write (*,*) "PY Potential. PCutoff ON:",PCUTOFF
            ENDIF
            IF (NITEMS.GT.10) THEN
               ! control which dimensions have periodic boundaries with a string 'XYZ', always put x before y before z.
               ! eg ...  Xz 20 30  specifies PBC on X and Z directions.  The X box size will be 20, the Z box size 30
               CALL READA(PBC)
               write (*,*) "PBCs are: ",PBC
               BOXLX=0
               BOXLY=0
               BOXLZ=0
               IF (SCAN(PBC,'Xx').NE.0) THEN
                  PARAMONOVPBCX=.TRUE.
                  CALL READF(BOXLX)       ! BOXLX is a scaling factor, not the actual box length!
                  BOXLX=BOXLX*PCUTOFF     ! now BOXLX is the actual box length
                  write(*,*) "PY Periodic Boundary Condition X active. BOXLX:",BOXLX
               ENDIF
               IF (SCAN(PBC,'Yy').NE.0) THEN
                  PARAMONOVPBCY=.TRUE.
                  CALL READF(BOXLY)
                  BOXLY=BOXLY*PCUTOFF
                  write(*,*) "PY Periodic Boundary Condition Y active. BOXLY:",BOXLY
               ENDIF
               IF (SCAN(PBC,'Zz').NE.0) THEN
                  PARAMONOVPBCZ=.TRUE.
                  CALL READF(BOXLZ)
                  BOXLZ=BOXLZ*PCUTOFF
                  write(*,*) "PY Periodic Boundary Condition Z active. BOXLZ",BOXLZ
               ENDIF
            ENDIF
! ALLOCATE(RBSITE(NRBSITES,3))
! 
! QQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQQ

! INTERFACING ReaxFF HERE
         ELSE IF (WORD.EQ.'REAXFF') THEN
            REAXFFT=.TRUE.
            CALL READA(REAXFFJOB)
            CALL READA(SYS)
            REAXFFJOB=TRIM(ADJUSTL(REAXFFJOB)) // ' ' // TRIM(ADJUSTL(SYS))
            DO J1=1,80
               IF (SYS(J1:J1).EQ.' ') THEN
                  LSYS=J1-1
                  GOTO 814
               ENDIF
            ENDDO
814         CONTINUE


! INTERFACING MOLPRO HERE
         ELSE IF (WORD.EQ.'MOLPRO') THEN
            MOLPRO=.TRUE.
            CALL READA(MOLPROJOB)
            CALL READA(SYS)
            MOLPROJOB=TRIM(ADJUSTL(MOLPROJOB)) // ' ' // TRIM(ADJUSTL(SYS))
            DO J1=1,80
               IF (SYS(J1:J1).EQ.' ') THEN
                  LSYS=J1-1
                  GOTO 815
               ENDIF
            ENDDO
815         CONTINUE

      ELSE IF (WORD.EQ.'QCIADDREP') THEN
         CALL READI(QCIADDREP)
         CALL READF(QCIADDREPEPS)
         CALL READF(QCIADDREPCUT)
         WRITE(*,'(A,I6,A)') 'keywords> Adding ',QCIADDREP,' repulsive sites along constraints'
      ELSE IF (WORD.EQ.'QCINOREPINT') THEN
         QCINOREPINT=.TRUE.

! 
! 
! QCHEM tells the program to read derivative information in
! QCHEM format.                                        - default FALSE
! 
         ELSE IF (WORD.EQ.'QCHEM') THEN
            QCHEM=.TRUE.
            CALL READA(QCHEMJOB)
            CALL READA(SYS)
            QCHEMJOB=TRIM(ADJUSTL(QCHEMJOB)) // ' ' // TRIM(ADJUSTL(SYS))
            DO J1=1,80
               IF (SYS(J1:J1).EQ.' ') THEN
                  LSYS=J1-1
                  GOTO 23
               ENDIF
            ENDDO
23          CONTINUE
! 
! QCHEMES tells the program to read derivative information in
! QCHEM format.                                        - default FALSE
! In this case we are going to search for stationary points corresponding to
! different electronic states in the space of orbital coefficients.
! 
         ELSE IF (WORD.EQ.'QCHEMES') THEN
            QCHEM=.TRUE.
            QCHEMES=.TRUE.
            CALL READA(QCHEMJOB)
            CALL READA(SYS)
            CALL READI(QCHEMESNAO)
            CALL READI(QCHEMESNMO)
            CALL READI(QCHEMESNELEC)
            IF (QCHEMESNMO*QCHEMESNAO.EQ.0) THEN
               PRINT '(A)',' keywords> ERROR *** either the number of MOs or AOs is zero'
               STOP
            ENDIF
            IF (QCHEMESNELEC.EQ.0) THEN
               PRINT '(A)',' keywords> WARNING *** number of electrons is zero'
               STOP
            ENDIF
            QCHEMESNZERO=QCHEMESNMO**2-QCHEMESNMO*QCHEMESNELEC
            PRINT '(A,3I8)',' keywords> Number of AOs, MOs and electrons=',QCHEMESNAO,QCHEMESNMO,QCHEMESNELEC
            PRINT '(A,I8)',' keywords> Number of zero Hessian eigenvalues set to ',QCHEMESNZERO
            QCHEMJOB=TRIM(ADJUSTL(QCHEMJOB)) // ' ' // TRIM(ADJUSTL(SYS))
            DO J1=1,80
               IF (SYS(J1:J1).EQ.' ') THEN
                  LSYS=J1-1
                  GOTO 923
               ENDIF
            ENDDO
923         CONTINUE
! 
! 
! qSPCFw  flexible water model introduced by Paesani et al. (JCP 125, 184507 (2006))
! Coded by Javier.
! 
         ELSE IF (WORD.EQ.'QSPCFW') THEN
            QSPCFWT=.TRUE.
! 
! qTIP4PF flexible water model introduced by Habershon et al. (JCP 131, 024501 (2009))
! Coded by Javier.
! 
         ELSE IF (WORD.EQ.'QTIP4PF') THEN
            QTIP4PFT=.TRUE.

            ! sn402: added (see comments at keyword FRQCONV)
            FRQCONV = 2.045483D13
            WRITE(*,*) "keywords> Set frequency conversion factor to the QTIP4PF default value: ", FRQCONV
            WRITE(*,*) "keywords> This corresponds to frequencies being given in radians/s"

         ELSE IF (WORD.EQ.'QUIPZ') THEN
            CALL READI(QUIPZ)

         ELSE IF (WORD.EQ.'QUIPLATT') THEN
             QUIPT=.TRUE.
             CALL READF(QUIPLATT(1,1))
             CALL READF(QUIPLATT(2,1))
             CALL READF(QUIPLATT(3,1))
             CALL READF(QUIPLATT(1,2))
             CALL READF(QUIPLATT(2,2))
             CALL READF(QUIPLATT(3,2))
             CALL READF(QUIPLATT(1,3))
             CALL READF(QUIPLATT(2,3))
             CALL READF(QUIPLATT(3,3))
! 
! RRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRRR
! 
! 
! Spherical container
! 
         ELSE IF (WORD.EQ.'RADIUS') THEN
            CONTAINER=.TRUE.
            CALL READF(RADIUS)
            RADIUS=RADIUS**2
! 
! Number of random rotations to try in minpermdist.
! 
         ELSE IF (WORD.EQ.'RANROT') THEN
            CALL READI(NRANROT)
! 
! integer seed for random number generator.
! 
         ELSE IF (WORD.EQ.'RANSEED') THEN
            CALL READI(NDUM)
            CALL SDPRND(NDUM)
            IF ((NDUM.LT.0).OR.(NDUM.GT.9999)) THEN
               ! 
               ! if we ever need more than 20000 searches from the same minimum
               ! then this could be a problem
               ! 
               DO J1=1,3*NATOMS
                  RANDOM=DPRAND()
               ENDDO
            ENDIF
            WRITE(*,'(A,I6)') ' SETTINGS Random number generator seed=',NDUM
! 
! TVB: Requests to print out pathway parameters necessary to calculate catastrophe
! ratios. Affects path routine only.
! 
         ELSE IF (WORD.EQ.'RATIOS') THEN
            RATIOS=.TRUE.
! 
! RBSYM defines the internal symmetry operations for each sort of rigid body
! coded via RBAAT.
! 
         ELSE IF (WORD.EQ.'RBSYM') THEN
            RBSYMT=.TRUE.
            INQUIRE(FILE='rbsymops',EXIST=RBSYMTEST)
            IF (RBSYMTEST) THEN
               OPEN(UNIT=1,FILE='rbsymops',STATUS='OLD')
               READ(1,*) NRBGROUP
               ALLOCATE(RBOPS(4,NRBGROUP))
               READ(1,*) ((RBOPS(J1,J2),J1=1,4),J2=1,NRBGROUP)
               PRINT '(A,I6)',' keywords> number of symmetry operations for rigid body=',NRBGROUP
               DO J1=1,NRBGROUP
                  PRINT '(A,I6)',' keywords> rigid-body symmetry operation', J1
                  RBOPS(4,J1) = RBOPS(4,J1)*ATAN(1.D0)/45.D0
                  PRINT '(3F20.10)',RBOPS(1:4,J1)
               ENDDO
            ELSE
               PRINT '(A)',' keywords> ERROR *** missing file rbsymops'
               STOP
            ENDIF
! 
! If READMASS is specified we read the masses from file masses.
! 
         ELSE IF (WORD.EQ.'READMASS') THEN
            READMASST=.TRUE.
! 
! If READPATH is specified with CALCRATES then the rates are calculated from the
! information in an existing path.info file without any stationary point searches.
! 
         ELSE IF (WORD.EQ.'READPATH') THEN
            READPATH=.TRUE.
! 
! If READSP is true then OPTIM will read minima and ts data in the pathsample format
! 
         ELSE IF (WORD.EQ.'READSP') THEN
            READSP=.TRUE.
! 
! READHESS tells the program to read a Hessian at the first step.
! 
         ELSE IF (WORD .EQ. 'READHESS') THEN
            READHESS=.TRUE.
! 
! READVEC "file" reads the eigenvalue and associated eigenvector corresponding
! to the reaction coordinate for use in a pathway calculation. The format
! is the same as that used for vector.dump. If there is more than one vector
! in the file the program reads down to the last entry.
! 
         ELSE IF (WORD(1:7) .EQ. 'READVEC') THEN
            READV=.TRUE.
! ELSE IF (WORD.EQ.'REBUILDSC') THEN
! CALL READF(REBUILDSC)
! 
! sf344> read in coordinates from path.xyz files for rigid bodies, and
! bring the frames in the best alignment
! 
         ELSE IF (WORD.EQ.'REALIGNXYZ') THEN
            REALIGNXYZ=.TRUE.
! 
! Whether to use a redopoints file if it exists.
! 
         ELSE IF (WORD.EQ.'REDOPATH') THEN
            REDOPATH=.TRUE.
            IF (NITEMS.GT.1) THEN
               CALL READF(REDOK)
               REDOKADD=.TRUE.
            ENDIF
            IF (NITEMS.GT.2) CALL READF(REDOFRAC)
! 
! Whether to use a redopoints file if it exists.
! 
         ELSE IF (WORD.EQ.'REDOPATHNEB') THEN
            REDOPATHNEB=.TRUE.
            REDOPATH=.TRUE.
            FREEZENODEST=.TRUE.
            FREEZETOL=-1.0D0
            IF (NITEMS.GT.1) CALL READI(REDOBFGSSTEPS)
! 
! Whether to use path.<n>.xyz files in the current directory
! 
         ELSE IF (WORD.EQ.'REDOPATHXYZ') THEN
            REDOPATHXYZ=.TRUE.
            REDOPATH=.TRUE.
            IF (NITEMS.GT.1) CALL READF(REDOK)
            IF (NITEMS.GT.2) CALL READF(REDOFRAC)

! How many times to re-attempt the path run on a transition state if the same TS
! is located multiple times
        ELSE IF (WORD.EQ.'REDOTS') THEN
            CALL READI(REDOTS)

! Whether to reduce the bond lengths for side chains during the connection runs.
! To be used together with CHARMM (and AMBER not yet).
! 
         ELSE IF (WORD.EQ.'REDUCEDBONDLENGTH') THEN
            REDUCEDBONDLENGTHT=.TRUE.
            CALL READF(BLFACTOR)
            IF (NITEMS.GT.2) CALL READA(UNSTRING)
            IF (TRIM(ADJUSTL(UNSTRING)).EQ.'CB') CBT=.TRUE.
! 
! Specifies that the eigenvector to be followed should be reoptimised
! in a BFGSTS search after the EF step and before the tangent space minimisation.
! This is probably not a good idea.
! 
         ELSE IF (WORD.EQ.'REOPT') THEN
            REOPT=.TRUE.
         ELSE IF (WORD.EQ.'REOPTIMISEENDPOINTS') THEN
            REOPTIMISEENDPOINTS=.TRUE.
! 
! coordinates to orthogonalise search directions to are to be found in
! points.repel
! 
         ELSE IF (WORD.EQ.'REPELTS') THEN
            REPELTST=.TRUE.
            IF (NITEMS.GT.1) CALL READF(REPELPUSH)
! 
! RESIZE x scales the radial distances by x on the first
! step only                                           - default n=1
! 
         ELSE IF (WORD .EQ. 'RESIZE') THEN
            CALL READF(RESIZE)
! 
! REVERSEUPHILL reverses the gradient component in the uphill direction for
! BFGSTS and turns off projection in MYLBFGS.
! 
         ELSE IF (WORD.EQ.'REVERSEUPHILL') THEN
            REVERSEUPHILLT=.TRUE.

         ELSE IF (WORD .EQ. 'CUDA') THEN
            CUDAT=.TRUE.
            IF (NITEMS .EQ. 1) THEN
               WRITE(*,'(A)') " keywords> You must specify a potential with keyword CUDA. "
               STOP
            END IF
            IF (NITEMS .GT. 1) THEN
               CALL READA(CUDAPOT)
            END IF
            IF (NITEMS .GT. 2) THEN
               WRITE(*,'(A)') " keywords> Too many arguments specified with keyword CUDA. "
               STOP
            ENDIF
            IF ((CUDAPOT .EQ. 'A') .AND. (.NOT. AMBER12T)) THEN
               WRITE(*,'(A)') " keywords> The AMBER12 keyword must be used with 'CUDA A'. "
               STOP
            END IF

      ELSE IF (WORD .EQ. 'CUDATIME') THEN
         CUDATIMET=.TRUE.

! ----------------------------------!
! hk286 > Generalised rigid body   !
! ----------------------------------!

         ELSE IF (WORD.EQ.'RIGIDINIT') THEN

            RIGIDINIT = .TRUE.
            ATOMRIGIDCOORDT = .TRUE.

            IF (NITEMS.EQ.2) THEN
               CALL READA(AAOPTION)
               IF (AAOPTION.EQ.'AACONVERGENCE') THEN
                   AACONVERGENCET = .TRUE.
               ELSE
                  PRINT '(A)','keyword> ERROR *** RIGIDINIT option unrecognised'
                  STOP
               END IF
            END IF

            CALL GENRIGID_READ_FROM_FILE ()
            ! sn402: added next line
            NOPT = 3*NATOMS
            IF (PERMDIST) THEN
               CALL GENRIGID_PERMDIST()
            ENDIF

        ! sn402: addition
         ELSE IF (WORD .EQ. 'RIGIDMOLECULES') THEN
            RIGIDMOLECULEST = .TRUE.

! specifies additional rings other than the usual ones in
! PHE, PRO, TYR, HIS, and TRP residues
         ELSE IF (WORD.EQ.'RING') THEN
            NURINGS = NURINGS + 1
            IF (NITEMS.EQ.6) THEN
               URINGS(NURINGS,0) = 5
               DO J1 = 1,5
                  CALL READI(URINGS(NURINGS,J1))
               ENDDO
            ELSE IF (NITEMS.EQ.7) THEN
               URINGS(NURINGS,0) = 6
               DO J1 = 1,6
                  CALL READI(URINGS(NURINGS,J1))
               ENDDO
            ENDIF

! 
! RINGPOLYMER specifies a ring polymer system with harmonic springs between
! NRP images of the same system that generally have different geometries.
! RPSYSTEM is a string specifying the system, e.g. LJ.
! RPIMAGES is the number of RP images.
! RPBETA is 1/kT in reduced units.
! RINGPOLYMER keyword takes the place of POINTS and must be the last
! keyword in the odata file before the points.
! 
         ELSE IF (WORD.EQ.'RINGPOLYMER') THEN
            RINGPOLYMERT=.TRUE.
            CALL READA(RPSYSTEM)
            CALL READI(RPIMAGES)
            CALL READF(RPBETA)
! 
! Sanity checks.
! 
            TEMPSTRING=TRIM(ADJUSTL(RPSYSTEM))
            IF (TEMPSTRING(1:2).EQ.' ') THEN
               PRINT '(A)','keyword> ERROR *** Ring polymer potential type is not set'
            ENDIF
            IF (RPIMAGES.LT.1) THEN
               PRINT '(A)','keyword> ERROR *** Ring polymer images too small, value is ',RPIMAGES
            ENDIF
            IF (MBPOLT) THEN
               NATOMSSAVE=NATOMS
               NATOMS=NATOMS/(3*RPIMAGES)
               CALL MBPOLINIT
               DO J1=1,RPIMAGES    ! number of copies of water dimer, trimer, etc.
                  DO J2=1,NATOMS/3 ! number of water molecules in the cluster
                     ZSYM(NATOMS*(J1-1)+3*(J2-1)+1)='O '
                     ZSYM(NATOMS*(J1-1)+3*(J2-1)+2)='H '
                     ZSYM(NATOMS*(J1-1)+3*(J2-1)+3)='H '
                  ENDDO
               ENDDO
               NATOMS=NATOMSSAVE
            ENDIF

            RETURN


! 
! RKMIN calculates a steepest-descent path using gradient only information
! with convergence criterion GMAX for the RMS force and initial precision
! EPS. A fifth order Runga-Kutta algorithm is used.
! 
         ELSE IF (WORD.EQ.'RKMIN') THEN
            RKMIN=.TRUE.
            IF (NITEMS.GT.1) CALL READF(GMAX)
            IF (NITEMS.GT.2) CALL READF(EPS)
! 
! ROT [JZ n or OMEGA n] sets the value of J_z, the angular
! momentum about the z axis or
! OMEGA, the corresponding angular velocity
! 
         ELSE IF (WORD .EQ. 'ROT') THEN
            RTEST=.TRUE.
            CALL READU(WORD)
            CALL READF(XX)
            IF (WORD.EQ.'JZ') THEN
               JZ=XX
            ELSE
               OMEGA=XX
            ENDIF
! 
! fix linear polymer at its ends
! 
         ELSE IF (WORD .EQ. 'RPFIX') THEN
            RPFIXT=.TRUE.
            print *, 'fixed ends'
! 
! make ring polymer system into linear polymer
! 
         ELSE IF (WORD.EQ.'RPH') THEN
            RPHT=.TRUE.
            IF (NITEMS.GT.1) CALL READF(RPHTEMP)
            IF (NITEMS.GT.2) CALL READI(RPHSLICES)
            IF (NITEMS.GT.3) CALL READF(RPHQMIN)
            IF (NITEMS.GT.4) CALL READF(RPHQMAX)
            IF (NITEMS.GT.5) CALL READI(RPHQBINS)
            IF (NITEMS.LE.5) THEN
               PRINT '(A)','keywords> ERROR *** insufficient arguments for RPH keyword'
               STOP
            ENDIF

         ELSE IF (WORD .EQ. 'RPLINEAR') THEN
            RPCYCLICT=.FALSE.
            print *, 'use linear polymer'

! ---------------------------- SANDBOX potential ---------------------

         ELSE IF (WORD .EQ. 'SANDBOX') THEN
           SANDBOXT = .TRUE.
           RBAAT = .TRUE.
           CALL SANDBOX_INPUT(6)
           IF (NITEMS.GT.1) THEN ! we have uniaxial particles as well in the mixture
             UNIAXT=.TRUE.
             CALL READI(nuniax)
             UNIAXARRAY(:)=.FALSE.
             DO i=1,nuniax
                UNIAXARRAY(i)=.TRUE.
!                WRITE(*,*) 'uniaxial particle: ', i
             END DO
           END IF
           NTSITES=num_atoms   !num_atoms is the total number of sites
           ALLOCATE(RBSITE(NTSITES*2/NATOMS,3)) ! will have to change implementation 
                                                ! to allow for different number of sites in n-ary systems

!           
! SSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSSS
! 
! 
! Save candidate TS`s in SQVV run.
! 
         ELSE IF (WORD == 'SAVECANDIDATES') THEN
            SAVECANDIDATES=.TRUE.
! 
! SCALE n sets the value of ISTCRT                             - default n=10
! 
         ELSE IF (WORD.EQ.'SCALE') THEN
            CALL READI(ISTCRT)
! 
! Specify that we are running in a SCore environment. Currently never used.
! 
         ELSE IF (WORD.EQ.'SCORE_QUEUE') THEN
            SCORE_QUEUE = .TRUE.
! 
! SEARCH specifies the value of INR, i.e. the search type.     - default n=0
! 
         ELSE IF (WORD.EQ.'SEARCH') THEN
            CALL READI(INR)
! 
! SETCHIRAL: do not connect minima where one or more chiral atoms are inverted
! Useful for PATHSAMPLE runs with AMBER.
! 
         ELSE IF (WORD.EQ.'SETCHIRAL') THEN
            SETCHIRAL=.TRUE.
            TURNOFFCHECKCHIRALITY=.TRUE.
! 
! Eigenvalue shift parameter.
! 
         ELSE IF (WORD .EQ. 'SHIFT') THEN
            CALL READF(SHIFTV)
! 
! Parameters for Edwin;s SiO2 model
! 
         ELSE IF (WORD .EQ. 'SILANE') THEN

            SILANET  = .TRUE.
            RBAAT    = .TRUE.
            NRBSITES = 5
            ALLOCATE(RBSITE(NRBSITES,3))
            NTSITES = NATOMS*NRBSITES/2

            CALL DEFSILANE()

         ELSE IF (WORD.EQ.'SIO2') THEN
            SIO2T=.TRUE.
            CALL READF(PARAM1)
            IF (NITEMS.GT.2) THEN
               CALL READF(PARAM2)
            ENDIF
            IF (NITEMS.GT.3) THEN
               CALL READF(PARAM3)
            ENDIF
            IF (NITEMS.GT.4) THEN
               CALL READF(PARAM4)
            ENDIF
         ELSE IF (WORD.EQ.'SIO2C6') THEN
            SIO2C6T=.TRUE.
            CALL READF(C6PP)
            CALL READF(C6MM)
            CALL READF(C6PM)
         ELSE IF (WORD.EQ.'SIO2P') THEN
            SIO2PT=.TRUE.
         ! Activate the iSLERP interpolation scheme for rigid bodies. Only has any effect
         ! if RIGIDINIT is also set.
         ELSE IF (WORD.EQ.'SLERP') THEN
             SLERPT=.TRUE.           
! 
! SQVV allows the first NIterSQVVGuessMax DNEB iterations to be done using
! SQVV - switches to LBFGS minimisation after NIterSQVVGuessMax iterations
! or if the RMS force goes below SQVVGuessRMSTol.
! 
         ELSE IF (WORD == 'SQVV') THEN
            SQVVGUESS=.TRUE.
            IF (NITEMS.GT.1) CALL READI(NITERSQVVGUESSMAX)
            IF (NITEMS.GT.2) CALL READF(SQVVGUESSRMSTOL)
         ELSE IF (WORD.EQ.'SSH') THEN
            SSHT=.TRUE.
!
! KLIM sets the radius in then wave vector space.
! SCA sets the multiplier on the stealthy potential and gradients.
!
         ELSE IF (WORD.EQ.'STEALTHY') THEN
             STEALTHYT=.TRUE.
             CALL READF(KLIM)
             IF (NITEMS.GT.2) THEN
                CALL READF(SCA)
             ELSE
                SCA=1
             END IF

         ELSE IF (WORD.EQ.'STEALTHYTEST') THEN
             STEALTV=.TRUE.
! 
! NSTEPMIN sets the minimum number of steps allowed before convergence.
! 
         ELSE IF (WORD .EQ. 'STEPMIN') THEN
             CALL READI(NSTEPMIN)
! 
! STEPS n sets the number of optimisation steps to perform
! per call to OPTIM                                    - default n=1
! If BFGSSTEPS is not specified then it is set to the same value as NSTEPS
! 
         ELSE IF (WORD .EQ. 'STEPS') THEN
            CALL READI(NSTEPS)
            IF (BFGSSTEPS.EQ.1) BFGSSTEPS=NSTEPS
! 
! Stillinger-David water potential - coded by Jeremy Richardson
! 
         ELSE IF (WORD.EQ.'SD') THEN
            SDT=.TRUE.
            CALL READI(SDOXYGEN)
            CALL READI(SDHYDROGEN)
            CALL READI(SDCHARGE)
            IF (SDOXYGEN*SDHYDROGEN.EQ.0) THEN
               PRINT '(A,2I6)', ' keyword> ERROR *** number of SD oxygens and hydrogens=',SDOXYGEN,SDHYDROGEN
               STOP
            ENDIF

            ! sn402: added (see comments at keyword FRQCONV)
            FRQCONV = 2.045483D13
            WRITE(*,*) "keywords> Set frequency conversion factor to the SD default value: ", FRQCONV
            WRITE(*,*) "keywords> This corresponds to frequencies being given in radians/s"

         ELSE IF (WORD.EQ.'STOCK') THEN
            STOCKT=.TRUE.
! RIGIDBODY=.TRUE.
! NRBSITES=1 ! used in current GMIN
            CALL READF(STOCKMU)
            CALL READF(STOCKLAMBDA)
! ALLOCATE(SITE(NRBSITES,3))
! 
! STOCKSPIN randomises the orientation of a Stockmayer cluster at any point in
! an optimisation where a dipole vector becomes aligned with the z axis (which
! make the phi angle for that dipole redundant).  STOCKZTOL is the amount by
! which cos(theta) may differ from 1.0 for alignment to be recognised.
! STOCKMAXSPIN is the maximum number of random orientations that will be attempted.
! 
         ELSE IF (WORD.EQ.'STOCKSPIN') THEN
            STOCKSPIN = .TRUE.
            CALL READF(STOCKZTOL)
            CALL READI(STOCKMAXSPIN)

         ELSE IF (WORD.EQ.'ST') THEN

            STOCKAAT = .TRUE.
            RBAAT    = .TRUE.
            STOCKEXP = 6.0D0
            CALL READF(STOCKMU)
            IF (NITEMS .GT. 2) THEN
               CALL READF(EFIELD)
               IF (EFIELD>0.0D0) EFIELDT = .TRUE.
            ENDIF
            IF (NITEMS .GT. 3) CALL READF(STOCKEXP)
            IF (NITEMS .GT. 4) THEN
               CALL READI(DPID)
               IF (DPID == 2) MRHO = STOCKEXP
            ENDIF       
            NRBSITES = 1
            ALLOCATE(RBSITE(NRBSITES,3))
            NTSITES = NATOMS*NRBSITES/2
! 
! STOPDIST specifies an alternative stopping criterion based on displacement
! between the first or last minimum and the furthest connected minimum.
! 
         ELSE IF (WORD.EQ.'STOPDISP') THEN
            CALL READF(STOPDISP)
            STOPDISPT=.TRUE.
! 
! In a CONNECT run, stop as soon as the initial minimum has a transition state
! connection.
! 
         ELSE IF (WORD.EQ.'STOPFIRST') THEN
            STOPFIRST=.TRUE.
! 
! SUMMARY n print a summary of the steps taken every n cycles  - default n=20
! 
         ELSE IF (WORD .EQ. 'SUMMARY') THEN
            IF (NITEMS.GT.1) CALL READI(NSUMMARY)
! 
! SYMCUT n RMS force below which symmetry subroutine is called - default 0.001
! 
         ELSE IF (WORD .EQ. 'SYMCUT') THEN
            CALL READF(SYMCUT)
! 
! TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT
! 
! 
! Tagged particle - atom in question has mass increased by TAGFAC in symmetry.f and inertia.f
! 
         ELSE IF (WORD.EQ.'TAG') THEN
            TAGT=.TRUE.
            NTAG=NTAG+1
            CALL READI(TAGNUM(NTAG))
            CALL READF(TAGFAC(TAGNUM(NTAG)))
         ELSE IF (WORD.EQ.'TANTYPE') THEN
            CALL READI(TANTYPE)
            GSTANTYPE = TANTYPE
! 
! Add a tetrahedral field to the potential of magnitude FTD.
! 
         ELSE IF (WORD.EQ.'TD') THEN
            FIELDT=.TRUE.
            TDT=.TRUE.
            CALL READF(FTD)
! 
! TIMELIMIT - in seconds - OPTIM will stop if this limit is exceeded.
! 
         ELSE IF (WORD.EQ.'TIMELIMIT') THEN
            CALL READF(TIMELIMIT)
! 
! TOLD n initial distance tolerance in symmetry subroutine     - default 0.0001
! 
         ELSE IF (WORD .EQ. 'TOLD') THEN
            CALL READF(TOLD)
! 
! TOLE n initial tolerance for the difference in principal moments
! of inertia divided by the sum of the principal moments
! in symmetry subroutine                                - default 0.0001
! 
         ELSE IF (WORD .EQ. 'TOLE') THEN
            CALL READF(TOLE)
! 
! Includes omega angles in the TWISTDIHE list.
! 
         ELSE IF (WORD.EQ.'TOMEGA') THEN
            TOMEGAC=.TRUE.
         ELSE IF (WORD.EQ.'TOSIPOL') THEN
            TOSIPOL=.TRUE.
            CALL READF(ALPHAP)
            CALL READF(ALPHAM)
            CALL READF(DAMP)
            WRITE(*,'(A)') ' Polarizabilities:'
            WRITE(*,'(A,F12.8,A,F12.8)') ' alpha+=',ALPHAP,' alpha-=',ALPHAM
            WRITE(*,'(A,F12.8,A)') ' damping coefficent=',DAMP,' per bohr'
! 
! TRAD n sets the trust radius to n                            - default n=4
! 
         ELSE IF (WORD .EQ. 'TRAD') THEN
            CALL READF(TRAD)
! 
! TRAP is used for the trap potential in EYtrap coded by Ersin Yurtsever.
! 
         ELSE IF (WORD .EQ. 'TRAP') THEN
            EYTRAPT=.TRUE.
            CALL READF(TRAPK)
            CALL READI(NTRAPPOW)
! 
! TRAPMK is used for the trap potential in EYtrap coded by Ersin Yurtsever.
! 
         ELSE IF (WORD .EQ. 'TRAPMK') THEN
            MKTRAPT=.TRUE.
! 
! Xantheas' TTM3-F water potential
! 
         ELSE IF (WORD.EQ.'TTM3') THEN
            TTM3T=.TRUE.

            ! sn402: added (see comments at keyword FRQCONV)
            FRQCONV = 2.045483D13
            WRITE(*,*) "keywords> Set frequency conversion factor to the TTM3 default value: ", FRQCONV
            WRITE(*,*) "keywords> This corresponds to frequencies being given in radians/s"
! 
! Includes sidechain angles in the TWISTDIHE list.
! 
         ELSE IF (WORD.EQ.'TSIDECHAIN') THEN
            TSIDECHAIN=.TRUE.

! jbr36 - Tunneling splitting active
         ELSE IF (WORD.EQ.'TSPLITTING') THEN
            TSPLITTINGT=.TRUE.
            WRITE(*,*) 'Splitting is not implemented fully yet'
            STOP
         ELSE IF (WORD.EQ.'TRUSTMODE') THEN
            TRUSTMODET=.TRUE.
            CALL READF(TMRATIO)
! 
! Add static dihedral angle potential.
! 
         ELSE IF (WORD.EQ.'TWIST') THEN
            TWISTT=.TRUE.
            IF (NITEMS.EQ.7) THEN
               CALL READI(ITWIST)
               CALL READI(JTWIST)
               CALL READI(KTWIST)
               CALL READI(LTWIST)
               CALL READF(TWISTF)
               CALL READF(TWISTREF)
               IF (TWISTREF.GT.LPI) THEN
                  TWISTREF=TWISTREF-2*LPI
                  WRITE(*,'(A,G20.10)') ' keyword> WARNING *** Twist reference angle changed to ',TWISTREF
               ELSEIF (TWISTREF.LT.-LPI) THEN
                  TWISTREF=TWISTREF+2*LPI
                  WRITE(*,'(A,G20.10)') ' keyword> WARNING *** Twist reference angle changed to ',TWISTREF
               ENDIF 

               IF (TWISTF.EQ.0.0D0) THEN
                  WRITE(*,'(A,I6,A,I6,A,G20.10)') ' keyword> WARNING *** Twisting force is zero, turning off twisting directive'
                  TWISTT=.FALSE.
               ELSE
                  WRITE(*,'(A,4I6,2(A,G20.10))') ' keyword> External dihedral potential for atoms ',ITWIST,JTWIST,KTWIST,LTWIST,
     &                                            ' magnitude=',TWISTF,' reference=',TWISTREF
               ENDIF
            ELSE
               WRITE(*,'(A)') "keywords> ERROR *** TWIST keyword needs takes 6 arguments"
               STOP
            ENDIF

! 
! Twist phi/psi dihedral angle nmode by xpert degrees before starting optimisation.
! 
         ELSE IF (WORD.EQ.'TWISTDIHE') THEN
            TWISTDIHET=.TRUE.
            CALL READF(PSRANDOM)
            WRITE(*,*) ' keywords> PSRANDOM=',PSRANDOM
            CALL READF(DPERT)
! 
! TWISTTYPE specifies the type of twisting done to guess transition states in GUESSTS for CHARMM
! 
         ELSE IF (WORD.EQ.'TWISTTYPE') THEN
            CALL READI(TWISTTYPE)
! 
! Double ended ts search.
! 
         ELSE IF (WORD.EQ.'TWOENDS') THEN
            TWOENDS=.TRUE.
            IF (NITEMS.GT.1) CALL READF(FSTART)
            IF (NITEMS.GT.2) CALL READF(FINC)
            IF (NITEMS.GT.3) CALL READI(NTWO)
            IF (NITEMS.GT.4) CALL READF(RMSTWO)
            IF (NITEMS.GT.5) CALL READI(NTWOITER)
            IF (NITEMS.GT.6) CALL READF(TWOEVAL)
! 
! UUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUUU
! 

         ELSE IF (WORD.EQ.'UNIAX') THEN
            UNIAXT = .TRUE.

         ELSE IF (WORD.EQ.'UNRES') THEN
            UNRST=.TRUE.
            CALL UNRESINIT
! Calphas and the side chain centroids are counted as atoms, but NOT the peptide bond centres.
            NATOM=2*nres
            IF (NATOM /= NATOMS) THEN
               WRITE(*,'(A)') 'No. of atoms in "coords" conflicts with that deduced from unres part of odata'
               CALL FLUSH(6)
               STOP
            ENDIF
            NINTS=2*nres-5+2*nside ! jmc change this depending on how want to deal with non-capping glycines!
! jmc NINTS was previously set in fetchz, but need either it or nvaru earlier (i.e. here)
! so may as well set it when we first know nres and nside.

            IF (ENDHESS.AND.(.NOT.ENDNUMHESS)) THEN
               PRINT *,'**ERROR - to calculate normal mode frequencies for UNRES, please specify ENDNUMHESS keyword'
               CALL FLUSH(6)
               STOP
            ELSEIF ((DUMPPATH.OR.DUMPALLPATHS).AND.(.NOT.ENDHESS)) THEN
               PRINT *,'**ERROR - to calculate normal mode frequencies for UNRES, please specify ENDHESS and ENDNUMHESS keywords'
               CALL FLUSH(6)
               STOP
            ENDIF

! DO J1=1,nres
! jmc c contains x,y,z for all the Calphas
! UNRX(2*J1-1)=c(1,J1)
! UNRY(2*J1-1)=c(2,J1)
! UNRZ(2*J1-1)=c(3,J1)
! jmc then x,y,z for the side chain centroids
! UNRX(2*J1)=c(1,J1+nres)
! UNRY(2*J1)=c(2,J1+nres)
! UNRZ(2*J1)=c(3,J1+nres)
! ENDDO

! new read replaces random configuration coordinates with alternative from file coords
            CALL UNEWREAD(UNRX,UNRY,UNRZ,NATOMS,FILTH,FILTHSTR)
            DO J1=1,nres
               ! 1,J1)=UNRX(2*J1-1)
               ! 2,J1)=UNRY(2*J1-1)
               ! 3,J1)=UNRZ(2*J1-1)
               ! 1,J1+nres)=UNRX(2*J1)
               ! 2,J1+nres)=UNRY(2*J1)
               ! 3,J1+nres)=UNRZ(2*J1)
            ENDDO
            CALL UPDATEDC
            CALL int_from_cart(.true.,.false.)
            CALL chainbuild
! jmc put coords in standard orientation (1st atom at 0,0,0 etc...) into UNR array.  Fixes problem in path,
! for calculating the step off the TS
            DO J1=1,nres
               UNRX(2*J1-1)=c(1,J1)
               UNRY(2*J1-1)=c(2,J1)
               UNRZ(2*J1-1)=c(3,J1)
               UNRX(2*J1)=c(1,J1+nres)
               UNRY(2*J1)=c(2,J1+nres)
               UNRZ(2*J1)=c(3,J1+nres)
            ENDDO
            CALL UNRSETZSYMATMASS
            IF (FILTH.NE.0) THEN
               OPEN(UNIT=20,FILE='coords.read',STATUS='REPLACE')
               CLOSE(20)
            ENDIF
            ALLOCATE(UREFCOORD(3*NATOMS),UREFPPSANGLE(3*NATOMS))
            IF (TWISTDIHET.OR.PERTDIHET.OR.GUESSTST.OR.CALCDIHE) THEN
               CALL UNRSETDIHE
            ENDIF
            IF (TWISTDIHET) THEN
               CALL UNRSTWISTDIHE(UNRX,UNRY,UNRZ,DMODE,DPERT)
            ENDIF
            IF (PERTDIHET) THEN
               CALL UNRSPERTDIHE(UNRX,UNRY,UNRZ,CHPMIN,CHPMAX,CHNMIN,CHNMAX,ISEED)
            ENDIF
            IF (CALCDIHE) THEN
               CALL UNREADREF(NATOMS)
               ! jmc readref2 leaves reference coords in unres c and internal coord arrays, so replace with UNR{X,Y,Z} here.
               DO J1=1,nres
                  ! 1,J1)=UNRX(2*J1-1)
                  ! 2,J1)=UNRY(2*J1-1)
                  ! 3,J1)=UNRZ(2*J1-1)
                  ! 1,J1+nres)=UNRX(2*J1)
                  ! 2,J1+nres)=UNRY(2*J1)
                  ! 3,J1+nres)=UNRZ(2*J1)
               ENDDO
               CALL UPDATEDC
               CALL int_from_cart(.true.,.false.)
            END IF

            DO J1=1,NATOMS
               Q(3*(J1-1)+1)=UNRX(J1)
               Q(3*(J1-1)+2)=UNRY(J1)
               Q(3*(J1-1)+3)=UNRZ(J1)
            ENDDO
! 
! USEDIAG enables the user to select DIAG or DIAG2 as the eigenvalue estimate in
! Rayleigh-Ritz routine secdiag. Default is currently two. 
! nsecdiag = 1 : approximate curvature with second order central differences,
! but use energies to do it.
! nsecdiag = 2 : approximate curvature with second order central differences,
! but use gradients to do it.  Same accuracy as nsecdiag==1, but
! less issues with numerical precision. The default method.
! nsecdiag = 3 : approximate curvature with first order forward finite differences.
! less accurate than central difference, but need half as many
! potential calls per iteration
! 
         ELSE IF (WORD.EQ.'USEDIAG') THEN
            CALL READI(NSECDIAG)

! 
! USEEV allows the lowest NUSEEV eigenvalues and associated eigenvectors to be
! used in second-order searches with efol.f90.
! 
         ELSE IF (WORD.EQ.'USEEV') THEN
            CALL READI(NUSEEV)
! 
! Number of BFGS updates before resetting, default=4
! MUPDATE  for standard potential minimisation
! XMUPDATE for Rayleigh-Ritz mininisation
! MMUPDATE for mindist routine
! NEBMUPDATE for DNEB minimisation
! INTMUPDATE for quasi-continuous interpolation minimisation
! GSUPDATE for growins strings
! GCUPDATE for great circle minimisation
! 
         ELSE IF (WORD.EQ.'UPDATES') THEN
            CALL READI(MUPDATE)
            IF (NITEMS.GT.2) CALL READI(XMUPDATE)
            IF (NITEMS.GT.3) CALL READI(MMUPDATE)
            IF (NITEMS.GT.4) CALL READI(NEBMUPDATE)
            IF (NITEMS.GT.5) CALL READI(INTMUPDATE)
            IF (NITEMS.GT.6) CALL READI(GSUPDATE)
            IF (NITEMS.GT.7) CALL READI(GCUPDATE)
! 
! VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
! 
! 
! VALUES n print the Hessian eigenvalues every n cycles        - default n=20
! 
         ELSE IF (WORD .EQ. 'VALUES') THEN
            CALL READI(NVALUES)
! 
! VARIABLES - keyword at the end of the list of options after which
! the general variables follow. NZERO is the number of zero
! eigenvalues, default 0.
! 
         ELSE IF (WORD.EQ.'VARIABLES') THEN
            VARIABLES=.TRUE.
            ZSYM(1:NATOMS)=' ' ! otherwise valgrind reports uninitialised

            IF (DUMMY_FRQCONV.EQ.0) THEN 
               FRQCONV = 1.0D0
            ELSE
               FRQCONV = DUMMY_FRQCONV
            ENDIF
            FRQCONV2 = FRQCONV*FRQCONV
            IF (NITEMS.GT.1) CALL READI(NZERO)

            RETURN
! 
! 
! VASP tells the program to read derivative information in
! VASP format.                                        - default FALSE
! 
         ELSE IF (WORD.EQ.'VASP') THEN
            VASP=.TRUE.
            BULKT=.TRUE.
            CALL READA(VASPJOB)
            VASPJOB=TRIM(ADJUSTL(VASPJOB))
! 
! VECTORS n prints the eigenvectors every n cycles             - default OFF
! 
         ELSE IF (WORD .EQ. 'VECTORS') THEN
            VECTORST=.TRUE.
            CALL READI(NVECTORS)
! 
! WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWW
! 
         ELSE IF (WORD == 'WARRANTY') THEN
            CALL WARRANTY
! 
! Welch parameters for Born-Meyer binary salt potentials.
! These are A++, A--, A+- and rho, in order, followed by
! 
         ELSE IF (WORD.EQ.'WELCH') THEN
            WELCH=.TRUE.
            CALL READF(APP)
            CALL READF(AMM)
            CALL READF(APM)
            CALL READF(RHO)
            CALL READF(XQP)
            CALL READF(XQM)
            CALL READF(ALPHAP)
            CALL READF(ALPHAM)
! 
! XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
! 
! 
! YYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYYY
! 
! 
! ZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZZ
! 
         ELSE IF (WORD.EQ.'ZEROS') THEN
            CALL READI(NZERO)
         ELSE
            CALL REPORT(' keywords> Unrecognized command '//WORD,.TRUE.)
            STOP
         ENDIF

         CALL FLUSH(6)
         GOTO 190

         RETURN
      END


      SUBROUTINE SORT4(N,NATOMS,A,F)
      IMPLICIT NONE
      INTEGER J1, L, N, J2, NATOMS, F(NATOMS), NTEMP
      DOUBLE PRECISION A(NATOMS), TEMP
C
      DO 20 J1=1,N-1
         L=J1
         DO 10 J2=J1+1,N
            IF (A(L).GT.A(J2)) L=J2
10       CONTINUE
         TEMP=A(L)
         A(L)=A(J1)
         A(J1)=TEMP
         NTEMP=F(L)
         F(L)=F(J1)
         F(J1)=NTEMP
20    CONTINUE
      RETURN
      END

