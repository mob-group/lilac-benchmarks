C   OPTIM: A program for optimizing geometries and calculating reaction pathways
C   Copyright (C) 1999-2006 David J. Wales
C   This file is part of OPTIM.
C
C   OPTIM is free software; you can redistribute it and/or modify
C   it under the terms of the GNU General Public License as published by
C   the Free Software Foundation; either version 2 of the License, or
C   (at your option) any later version.
C
C   OPTIM is distributed in the hope that it will be useful,
C   but WITHOUT ANY WARRANTY; without even the implied warranty of
C   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C   GNU General Public License for more details.
C
C   You should have received a copy of the GNU General Public License
C   along with this program; if not, write to the Free Software
C   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
C
      SUBROUTINE OPTIM(F1,F2,FLSTRING)
      USE COMMONS
      USE KEY
      USE MODTWOEND
      USE MODHESS
      USE MODNEB
      USE VECCK
      USE ZWK
      USE MODCHARMM
      USE MODUNRES
      USE NEWNEBMODULE
      USE NEWCONNECTMODULE
      USE KEYNEB, NNNIMAGE=>NIMAGE
      USE MODGUESS
      USE MODMEC
      USE PORFUNCS
      USE CHARUTILS
      USE CHIRALITY
      USE INTCOMMONS, ONLY : INTINTERPT, DESMINT, NATINT, INTMINPERMT !msb50 remove last
      USE INTCUTILS, ONLY : INTSETUP, INTCLEANUP ! msb50 last
      USE MODAMBER9, ONLY : DUMPMODEN,CHIARRAY1,GOODSTRUCTURE1
      USE OPEP_INTERFACE_MOD, ONLY : OPEP_FINISH
!     USE BENCHMARKS, ONLY : MINBMT, MINBM
! hk286
      USE GENRIGID

      IMPLICIT NONE
! subroutine parameters  
      INTEGER F1,F2
      CHARACTER(LEN=80) FLSTRING
      CHARACTER(LEN=2) DUMMYS

      INTEGER J1, J2, NPCALL, ECALL, FCALL, SCALL, HORDER, NATOMSSAVE, SUNIT, FUNIT, MULTIINR, INRSAVE
      DOUBLE PRECISION VNEW(NOPT), ENERGY, EVALMIN, RMS, VECS(NOPT), QSAVE(NOPT),
     1  QPLUS(NOPT), LGDUMMY(NOPT),RMSINITIAL,RMSFINAL,E1,E2, RMAT(3,3),
     2  DIST, OVEC(3), H1VEC(3), H2VEC(3), Q(NOPT), EINITIAL, EFINAL, 
     3  ETIME, FTIME, STIME, DPRAND, DCOORDS(NOPT), INTFREEZETOLSAVE,
     4  ETS, EPLUS, EMINUS, SLENGTH, DISP, GAMMA, NTILDE, 
     5  FRQSTS(NOPT), FRQSPLUS(NOPT), FRQSMINUS(NOPT), QMINUS(NOPT), DISTSF
      DOUBLE PRECISION THTEMP(NOPT)
      CHARACTER ESTRING*87, GPSTRING*80, NSTRING*80, FSTRING*80, FNAME*13, FNAMEV*18, 
     1          ITSTRING*22, EOFSSTRING*15
      CHARACTER(LEN=80) FNAMEF
      CHARACTER(LEN=20) EFNAME
      DOUBLE PRECISION, ALLOCATABLE, DIMENSION(:) :: QW
      LOGICAL LDUMMY, BFGSTSSAVE, MFLAG, POTCALL, GOODENDPOINT
      INTEGER FRAME
      LOGICAL PVFLAG
      COMMON /PVF/ PVFLAG
C     COMMON /VN/ VNEW   !  common SV was also deleted
      COMMON /STRINGS/ ESTRING, GPSTRING, NSTRING, FSTRING
      COMMON /PCALL/ NPCALL, ECALL, FCALL, SCALL, ETIME, FTIME, STIME
      LOGICAL PATHT, DRAGT
      INTEGER NPATHFRAME, NCDONE
      COMMON /RUNTYPE/ DRAGT, PATHT, NPATHFRAME
      LOGICAL CONNECTT, DUMPPATH, READPATH, CALCRATES, STOPFIRST
      DOUBLE PRECISION TEMPERATURE, HRED, INERTIA(3,3), DIST2, DPFCT, DU(3)
      INTEGER NCONNECT, ISTAT, VERSIONTEMP
      COMMON /CONN/ STOPFIRST, CONNECTT, NCONNECT, DUMPPATH, READPATH, CALCRATES, TEMPERATURE, HRED
      LOGICAL KNOWE, KNOWG, KNOWH
      COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
      CHARACTER(LEN=5) ZSYMSAVE
      CHARACTER(LEN=20) PINFOSTRING
      COMMON /SYS/ ZSYMSAVE
      LOGICAL PATHFAILT ! JMC
! hk286
      DOUBLE PRECISION :: XRIGIDCOORDS(3*NATOMS), XRIGIDGRAD(3*NATOMS)
      INTEGER GETUNIT
      COMMON /OEPATH/ ETS,EPLUS,EMINUS
! Print OPTIM version in the output
      !VERSIONTEMP=25661
      !WRITE(*, '(A,I5)') ' OPTIM> version r',VERSIONTEMP
C
C  Dynamic memory allocation
C
      ALLOCATE (FROZEN(NATOMS),ZSYM(NATOMS),NR(NATOMS))
      IF (DEBUG) WRITE(*,*) ' OPTIM> allocated FROZEN and ZSYM with dimension NATOMS=',NATOMS
      ALLOCATE (FROZENRES(NATOMS))
C      STPMAX(:)=0.0D0

      FILTH=F1 ; FILTH2=F2 ; FILTHSTR=TRIM(ADJUSTL(FLSTRING))
      KNOWE=.FALSE. ; KNOWG=.FALSE. ; KNOWH=.FALSE.
      RBATOMSMAX=10

      CALL KEYWORDS(Q)
      ALLOCATE (STPMAX(NOPT))
      IF (UNRST.AND.(RKMIN.OR.BSMIN.OR.(INR.GT.-1))) THEN
         PRINT '(A)','UNRES not coded for requested optimisation option'
         CALL FLUSH(6)
         STOP
      ENDIF

      ALLOCATE(SHIFTL(NOPT))
      IF (PYGPERIODICT.OR.PYBINARYT) CALL INITIALISEPYGPERIODIC
      IF (MULTISITEPYT) CALL DEFINEPYMULTISITES
      IF (LJGSITET) CALL DEFINELJMULTISITES
      IF (LOCALPERMDIST) THEN
         ALLOCATE(RBGROUP(RBATOMSMAX))
         ALLOCATE(RBNINGROUP(NATOMS))
      ENDIF

      IF (CHEMSHIFT) CALL CAMSHIFTREAD(20)
      IF (MIEFT) CALL MIEF_INI()

      CALL CPU_TIME(TSTART)
C     IF (CONNECTT.AND.NEWNEBT) THEN
C        PRINT*,'WARNING - cannot use old connect with new neb, changing to old neb'
C        NEWNEBT=.FALSE.
C        NEBT=.TRUE.
C     ENDIF
      IF ((FILTH2.EQ.0).AND.(FILTH.NE.0)) WRITE(FILTHSTR,'(I10)') FILTH ! otherwise FILTHSTR isn;t set correctly.
      IF (REPELTST) ALLOCATE(REPELTS(NOPT,100)) ! PREVIOUS TS GEOMETRIES TO AVOID
      IF (CHECKINDEX) ALLOCATE(VECCHK(NOPT,MAX(NUSEEV,HINDEX,1))) ! vectors to orthogonise to
      ALLOCATE(ZWORK(NOPT,MAX(NUSEEV,HINDEX,1)))                  ! partial eigenvectors storage
      IF (TWOENDS.OR.CONNECTT.OR.NEBT.OR.NEWNEBT.OR.DRAGT.OR.GUESSPATHT.OR.MECCANOT.OR.MORPHT.OR.GREATCIRCLET.OR.BHINTERPT
     & .OR.BISECTT.OR.CLOSESTALIGNMENT.OR.GROWSTRINGT) THEN
         ALLOCATE(FIN(NOPT),START(NOPT))
      ENDIF

      NPCALL=0
      ECALL=0
      FCALL=0
      SCALL=0
      ETIME=0
      FTIME=0
      STIME=0
      FRAME=1
      IF (FILTH.EQ.0) THEN
         FNAME='points'
         FNAMEF='points.final'
         FNAMEV='vector.dump'
         EFNAME='energies'
      ELSE
         WRITE(FNAME,'(A)') 'points.'//TRIM(ADJUSTL(FILTHSTR))
         WRITE(FNAMEF,'(A)') 'points.final.'//TRIM(ADJUSTL(FILTHSTR))
         WRITE(FNAMEV,'(A)') 'vector.dump.'//TRIM(ADJUSTL(FILTHSTR))
         WRITE(EFNAME,'(A)') 'energies.'//TRIM(ADJUSTL(FILTHSTR))
      ENDIF

      IF (PRINTPTS.OR.MORPHT) THEN
         OPEN(UNIT=1,FILE=FNAME,STATUS='UNKNOWN')
         OPEN(UNIT=2,FILE=EFNAME,STATUS='UNKNOWN')
      ENDIF

      IF (BFGSSTEP) NSTEPS=NSTEPS+1
      IF (DUMPV) OPEN(UNIT=44,FILE=FNAMEV,STATUS='UNKNOWN')
      IF (DUMPMAG) OPEN(UNIT=321,FILE='mag_per_step', STATUS='UNKNOWN')
      CALL FETCHZ(Q) ! NOPT is not defined until here!

      IF (AMBERT .OR. NABT) CALL SET_CHECK_CHIRAL(Q,NATOMS,GOODSTRUCTURE1,CHIARRAY1)

      IF (CHECKDT) CALL CHECKDRVTS(Q)

      IF ((INTCONSTRAINTT.OR.INTLJT).AND.(NCONGEOM.GE.2)) THEN
!
! Set up all the constraints and repulsions for zero frozen atoms.
!
         IF (.NOT.ALLOCATED(CONI)) THEN
            ALLOCATE(CONI(INTCONMAX),CONJ(INTCONMAX),CONDISTREF(INTCONMAX),CONCUT(INTCONMAX),CONOFFLIST(INTCONMAX),
     &               CONOFFTRIED(INTCONMAX))
            CONOFFTRIED(1:INTCONMAX)=.FALSE.
            ALLOCATE(REPI(NREPMAX),REPJ(NREPMAX),NREPI(NREPMAX),NREPJ(NREPMAX),REPCUT(NREPMAX),NREPCUT(NREPMAX))
         ENDIF
         INTFREEZETOLSAVE=INTFREEZETOL
         INTFREEZETOL=-1.0D0
         PRINT *,'OPTIM> NCONGEOM=',NCONGEOM
         CALL MAKE_CONPOT(NCONGEOM,CONGEOM) 
         INTFREEZETOL=INTFREEZETOLSAVE
      ENDIF


      IF (BLNT) THEN ! PUT THE FULL BLN LETTERS INTO ZSYM
         DO J1=1,NATOMS
            ZSYM(J1)=BEADLETTER(J1) // 'L'
         ENDDO
      ENDIF
      ZSYMSAVE=ZSYM(NATOMS)
      IF (UNRST) ALLOCATE(INTSTEP(NINTS))

      IF (INTMINT.OR.DESMINT.OR.INTINTERPT.OR.NATINT) THEN
         CALL INTSETUP
      ENDIF
!DL-find interface for quantum rate calculations
      IF (INSTANTONSTARTDUMPT.or.CLASSICALRATEST.or.INSTANTONOPTT.or.INSTANTONRATET.or.TSPLITTINGT) THEN
          CALL dlf_interfaceOPTIM(NATOMS) !call to dlfind to calculate quantum rates
      ELSE

C ----- TESTING AND BENCHMARKS HERE ----
!     IF (MINBMT) CALL MINBM
C---------------------------------------
C
C  The next line is needed for parallel runs of clean/dirty maidens.
C  Should now be unnecessary, since all the necessary functionality is in PATHSAMPLE,
C  including single-ended ts searches.
C
C     CALL SYSTEM('echo junk > odata.read')
C     BULKT=.FALSE.
      IF ((ZSYMSAVE.EQ.'TT').AND.(PARAM1.NE.0.0D0)) BULKT=.TRUE.
      IF ((ZSYMSAVE.EQ.'SW').AND.(PARAM1.NE.0.0D0)) BULKT=.TRUE.
      IF ((ZSYMSAVE.EQ.'Z2').AND.(PARAM1.NE.0.0D0)) BULKT=.TRUE.
      IF ((ZSYMSAVE.EQ.'ZF').AND.(PARAM1.NE.0.0D0)) BULKT=.TRUE.
      IF ((ZSYMSAVE.EQ.'SM').AND.(PARAM1.NE.0.0D0)) BULKT=.TRUE.
      IF ((ZSYMSAVE.EQ.'QI').AND.(PARAM1.NE.0.0D0)) BULKT=.TRUE.
      IF ((ZSYMSAVE.EQ.'RM').AND.(PARAM1.NE.0.0D0)) BULKT=.TRUE.
C     LDUMMY=(NATOMS.EQ.64).AND.(ZSYMSAVE(1:1).EQ.'W') 
      BULKT=(BULKT.OR.
     1      (ZSYMSAVE.EQ.'ME').OR.(ZSYMSAVE.EQ.'P6').OR.
     2      (ZSYMSAVE.EQ.'SC').OR.(ZSYMSAVE.EQ.'MS').OR.
     3      (ZSYMSAVE.EQ.'MP').OR.(ZSYMSAVE.EQ.'JM').OR.
     3      (ZSYMSAVE.EQ.'DS').OR.
     4      (ZSYMSAVE.EQ.'LP').OR.(ZSYMSAVE.EQ.'LS').OR.(ZSYMSAVE.EQ.'LC').OR.
     5      (ZSYMSAVE.EQ.'LK'))

      IF (BULKT.AND.NEWCONNECTT.AND.(.NOT.QUIPT)) THEN
         CALL INITIALBOX(Q)
         CALL INITIALBOX(FIN)
      ENDIF

C     NOSHIFT=(FIELDT.OR.BFGSMINT.OR.BSMIN.OR.RKMIN.OR.NOHESS.OR.BFGSSTEP)
C     IF (INR.GE.0) NOSHIFT=.FALSE. ! changed INR default value in keywords to -1
C                                     This isn;t good enough - if we are doing a path and
C                                     EFOL is called once to get the eigenvalues and eigenvectors
C                                     then we need to SHIFT.
C  Shifting only needs to be turned off if we are doing VARIABLES or BFGSTS without NOIT;
C  account for this in BFGSTS.
C     IF (NOIT) NOSHIFT=.FALSE.
C     IF (VARIABLES) NOSHIFT=.TRUE.
C     NOSHIFT=(FIELDT.OR.BFGSMINT.OR.BSMIN.OR.RKMIN.OR.NOHESS.OR.BFGSSTEP)
C     IF (INR.GE.0) NOSHIFT=.FALSE. ! changed INR default value in keywords to -1
C     IF (NOIT) NOSHIFT=.FALSE.
C     IF (VARIABLES) NOSHIFT=.TRUE.
C
C  Jump back here for MULTIJOB runs.
C
      MULTIINR=INR
963   CONTINUE
C
C  Resize the system if required.
C
      IF (RESIZE.NE.1.0D0) THEN
         PRINT*,'Scaling coordinates by ',RESIZE
         IF (ZSYMSAVE(1:1).EQ.'W') THEN
            DO J1=1,3*(NATOMS/2)
               Q(J1)=Q(J1)*RESIZE
            ENDDO
         ELSE
            DO J1=1,NOPT
               Q(J1)=Q(J1)*RESIZE
            ENDDO
         ENDIF
      ENDIF
      IF (CASTEP.AND.PRESSURE) CALL CLATMIN(Q,VNEW)
C
C  Rescale distances if we are doing LJ (input assumed for SIG=3.4 A)
C
      IF ((ZSYMSAVE.EQ.'AR').AND.(ZSYM(1).NE.'CA')) THEN
         DO J1=1,NATOMS
            Q(3*(J1-1)+1)=Q(3*(J1-1)+1)/3.4D0
            Q(3*(J1-1)+2)=Q(3*(J1-1)+2)/3.4D0
            Q(3*(J1-1)+3)=Q(3*(J1-1)+3)/3.4D0
         ENDDO
      ENDIF
C
C  Scale side chain bond lengths for CHARMM or AMBER for connection runs
C
      IF (REDUCEDBONDLENGTHT) THEN
         IF (CHRMMT) CALL CHREDUCEDBONDLENGTH(Q,BLFACTOR,CBT)
C         IF (AMBERT.OR.NABT) CALL AMREDUCEDBONDLENGTH(Q,BLFACTOR,CBT)
      ENDIF
C
C     J2=0
C     DO J1=1,NATOMS
C        IF (IATNUM(J1).NE.0) J2=J2+1
C     ENDDO
C     NREAL=J2
C
C     CALL GMETRY(0,VEC,Q)
C
C  Don't call symmetry if we're doing Fenske-Hall
C
      IF ((ZSYMSAVE.NE.'FH').AND.(.NOT.VARIABLES).AND.(.NOT.AMBERT).AND.(.NOT.CHRMMT).AND.
     1 (.NOT.NABT).AND.(.NOT.UNRST).AND.(.NOT.RINGPOLYMERT).AND.(.NOT.AMBER12T).AND.(.NOT.BULKT)) THEN
C
C  For W1/W2/W3/W4 potentials Q contains the centre of mass coordinates
C  followed by the Euler angles and the number of molecules is NATOMS/2
C
         IF (ZSYMSAVE(1:1).EQ.'W') THEN
            ALLOCATE(QW(9*(NATOMS/2)))
C           DO J2=1,NOPT+IADD
C              QSAVE(J2)=Q(J2)
C           ENDDO
C           DO J2=1,NATOMS  !  WCOMMENT
            DO J2=1,NATOMS/2
               CALL CONVERT(Q(3*(J2-1)+1),Q(3*(J2-1)+2),Q(3*(J2-1)+3),
C    1                      Q(3*(NATOMS+J2-1)+1),Q(3*(NATOMS+J2-1)+2),Q(3*(NATOMS+J2-1)+3),
     1                      Q(3*(NATOMS/2+J2-1)+1),Q(3*(NATOMS/2+J2-1)+2),Q(3*(NATOMS/2+J2-1)+3),
     2                      OVEC,H1VEC,H2VEC)
               QW(9*(J2-1)+1)=OVEC(1)
               QW(9*(J2-1)+2)=OVEC(2)
               QW(9*(J2-1)+3)=OVEC(3)
               QW(9*(J2-1)+4)=H1VEC(1)
               QW(9*(J2-1)+5)=H1VEC(2)
               QW(9*(J2-1)+6)=H1VEC(3)
               QW(9*(J2-1)+7)=H2VEC(1)
               QW(9*(J2-1)+8)=H2VEC(2)
               QW(9*(J2-1)+9)=H2VEC(3)
            ENDDO

!            QW(1:3) = (/0.000000000000000, 0.000000000000000, -6.5098030735366103E-002/) 
!            QW(4:6) = (/0.000000000000000, 0.7569503272636612, 0.5207842458829288/)      
!            QW(7:9) = (/0.000000000000000,-0.7569503272636612, 0.5207842458829288/)

C           NATOMS=NATOMS*3  ! WCOMMENT
            NATOMSSAVE=NATOMS
            NATOMS=(NATOMS/2)*3
            CALL SYMMETRY(HORDER,.TRUE.,QW,INERTIA)
C           NATOMS=NATOMS/3  ! WCOMMENT
            NATOMS=NATOMSSAVE
C           DO J2=1,NOPT+IADD
C              Q(J2)=QSAVE(J2)
C           ENDDO
         ELSE IF (NTIPT) THEN

            CALL DEFTIP4(DC6CC, DC6CC, DC6CC)
            ALLOCATE(QW(9*(NATOMS/2)))
            J2       = NRBSITES
            NRBSITES = 3
            CALL SITEPOS(Q, QW)
            NRBSITES = J2

            NATOMSSAVE=NATOMS
            NATOMS=(NATOMS/2)*3
            CALL SYMMETRY(HORDER,.TRUE.,QW,INERTIA)
            NATOMS=NATOMSSAVE
         ELSE IF (PAHAT) THEN
            ALLOCATE(QW(3*NRBSITES*(NATOMS/2)))

            CALL SITEPOS(Q, QW)

            NATOMSSAVE=NATOMS
            NATOMS=(NATOMS/2)*NRBSITES

            CALL SYMMETRY(HORDER,.TRUE.,QW,INERTIA)

            NATOMS=NATOMSSAVE
         ELSE IF (STOCKAAT) THEN

            !CALL DEFSTOCK(STOCKMU, DU, DPFCT)
            ALLOCATE(QW(3*NRBSITES*(NATOMS/2)))
            CALL SITEPOS(Q, QW)

            NATOMSSAVE=NATOMS
            NATOMS=(NATOMS/2)*NRBSITES
            CALL SYMMETRY(HORDER,.TRUE.,QW,INERTIA)
            NATOMS=NATOMSSAVE
         ELSE
            CALL SYMMETRY(HORDER,.TRUE.,Q,INERTIA)
         ENDIF
      ENDIF
      IF (ADMT) CALL ADM(Q)
C
C  Eigenvalue shifting
C
      DO J1=1,NOPT
         SHIFTL(J1)=0.0D0
      ENDDO
!
! Set the  SHIFTL elements even if NOHESS is .TRUE. 
! This should do nothing if we never use a Hessian, but allows for a second-order
! pathway search even if NOHESS is used for the transition state search.
!
      IF (ZSYM(1).EQ.'CK') THEN
         WRITE(*,'(A,G20.10)') ' OPTIM> Using z rotational ev shift=',SHIFTV
         SHIFTL(6)=SHIFTV
      ELSE IF (TWISTT) THEN
         WRITE(*,'(A,G20.10)') ' OPTIM> Using x,y,z trans and z rotation ev shift=',SHIFTV
         SHIFTL(1)=SHIFTV
         SHIFTL(2)=SHIFTV
         SHIFTL(3)=SHIFTV
         SHIFTL(6)=SHIFTV
         WRITE(*,'(A)') ' OPTIM> x,y translational and z rotational ev shifting'
      ELSE IF (PULLT.OR.EFIELDT) THEN
         WRITE(*,'(A,G20.10)') ' OPTIM> Using x,y translational and z rotational ev shift=',SHIFTV
         SHIFTL(1)=SHIFTV
         SHIFTL(2)=SHIFTV
!         SHIFTL(3)=SHIFTV
         SHIFTL(6)=SHIFTV
!         WRITE(*,'(A)') ' OPTIM> x,y translational and z rotational ev shifting'
      ELSE IF (RTEST) THEN
         IF (JZ.NE.0.0D0) THEN
            WRITE(*,'(A,G20.10)') ' OPTIM> Using trans/z rotation ev shift=',SHIFTV
            SHIFTL(3)=SHIFTV
            SHIFTL(6)=SHIFTV
         ELSE
            WRITE(*,'(A,G20.10)') ' OPTIM> Using z rotation ev shift=',SHIFTV
            SHIFTL(6)=SHIFTV
         ENDIF
C     ELSE IF (NOSHIFT) THEN
C        WRITE(*,'(A)') 'No ev shifting'
      ELSE IF (TWOD) THEN ! z components are shifted in potential.f
         IF (EYTRAPT.OR.(ZSYMSAVE.EQ.'BE')) THEN
            SHIFTL(6)=SHIFTV
            WRITE(*,'(A)') ' OPTIM> z rotational ev shifting'
         ELSE IF (.NOT.BULKT) THEN
            SHIFTL(1)=SHIFTV
            SHIFTL(2)=SHIFTV
            SHIFTL(6)=SHIFTV
            WRITE(*,'(A)') ' OPTIM> x,y translational and z rotational ev shifting'
         ELSE
            SHIFTL(1)=SHIFTV
            SHIFTL(2)=SHIFTV
            WRITE(*,'(A)') ' OPTIM> x,y translational ev shifting'
         ENDIF
      ELSE IF (MKTRAPT) THEN
         WRITE(*,'(A,G20.10)') ' OPTIM> No eigenvalues shifted'
      ELSE IF (EYTRAPT.OR.(ZSYMSAVE.EQ.'BE')) THEN
         SHIFTL(4)=SHIFTV
         SHIFTL(5)=SHIFTV
         SHIFTL(6)=SHIFTV
         WRITE(*,'(A,G20.10)') ' OPTIM> Using rotational ev shift=',SHIFTV
      ELSE IF (BULKT) THEN
         WRITE(*,'(A,G20.10)') ' OPTIM> Using translational ev shift=',SHIFTV
         DO J1=1,3
            SHIFTL(J1)=SHIFTV
         ENDDO
      ELSE IF (ZSYM(NATOMS).EQ.'TH') THEN
         WRITE(*,'(A)') ' OPTIM> Using rotational ev shifting'
         SHIFTL(4)=SHIFTV
         SHIFTL(5)=SHIFTV
         SHIFTL(6)=SHIFTV
! hk286 - need to specify all 6?
      ELSE IF (GTHOMSONT) THEN
         IF (GTHOMMET .EQ. 5) THEN
            WRITE(*,'(A,G20.10)') ' OPTIM> Using translational/rotational ev shift=',SHIFTV
         ELSE
            WRITE(*,'(A)') ' OPTIM> z rotational ev shifting'
         ENDIF
         SHIFTL(1)=SHIFTV
         SHIFTL(2)=SHIFTV
         SHIFTL(3)=SHIFTV
         SHIFTL(4)=SHIFTV
         SHIFTL(5)=SHIFTV
         SHIFTL(6)=SHIFTV
      ELSE IF (MIEFT.OR.NOTRANSROTT) THEN
         WRITE(*,'(A,G20.10)') ' OPTIM> No eigenvalues shifted'
      ELSE IF (VARIABLES) THEN
         WRITE(*,'(A,I6,A)') ' OPTIM> ',NZERO,' zero eigenvalues shifted'
         DO J1=1,NZERO
            SHIFTL(J1)=SHIFTV
         ENDDO
      ELSE
         WRITE(*,'(A,G20.10)') ' OPTIM> Using translational/rotational ev shift=',SHIFTV
         DO J1=1,6
            SHIFTL(J1)=SHIFTV
         ENDDO
      ENDIF

      !js850> 
      IF (ZSYM(NATOMS).EQ.'SS') THEN
         ALLOCATE(SOFT_SPHERE_RADII(NATOMS))
         FUNIT=GETUNIT()
         OPEN(FUNIT,FILE='radii',STATUS='OLD')
         DO J1=1,NATOMS
            READ(FUNIT,*) SOFT_SPHERE_RADII(J1)
         ENDDO
         CLOSE(FUNIT)
      ENDIF

      IF ((INR.EQ.6).OR.(INR.EQ.7).OR.(INR.EQ.8)) ISTCRT=3

      IF ((NSTEPS.EQ.0).AND.(.NOT.(ENDHESS.OR.ENDNUMHESS))) THEN
         CALL FLUSH(6)
         STOP
      ENDIF

C     J2=0
C     DO J1=1,NATOMS
C        IF (IATNUM(J1).NE.0) J2=J2+1
C     ENDDO
C     NREAL=J2

555   CONTINUE ! JUMP BACK TO HERE AFTER REOPTIMISATION OF BAD ENDPOINTS
      IF (RPHT) THEN
         CALL RPH
      ELSE IF (MCPATHT) THEN
         CALL MCPATH
      ELSE IF (MCPATH2T) THEN
         CALL MCPATH2
      ELSE IF (CALCRATES.AND.READPATH) THEN
         CALL RATES(NATOMS,NINTS) ! JMC IF UNRES, THEN HAVE NON-ZERO NINTS, OTHERWISE JUST PASS 0
      ELSE IF (GUESSPATHT.AND.(.NOT.CONNECTT)) THEN
         CALL NEWMINDIST(FIN,Q,NATOMS,DISTSF,BULKT,TWOD,ZSYMSAVE,.FALSE.,RIGIDBODY,DEBUG,RMAT)
         WRITE(*,'(A,F12.2)') ' OPTIM> distance between start and finish=',DISTSF
         IF (UNRST) THEN
            DO J1=1,NRES
               C(1,J1)=Q(6*(J1-1)+1)
               C(2,J1)=Q(6*(J1-1)+2)
               C(3,J1)=Q(6*(J1-1)+3)
               C(1,J1+NRES)=Q(6*(J1-1)+4)
               C(2,J1+NRES)=Q(6*(J1-1)+5)
               C(3,J1+NRES)=Q(6*(J1-1)+6)
            ENDDO
            CALL UPDATEDC
            CALL INT_FROM_CART(.TRUE.,.FALSE.)
            CALL GEOM_TO_VAR(NINTS,Q(1:NINTS))
            DO J1=1,NRES
               C(1,J1)=FIN(6*(J1-1)+1)
               C(2,J1)=FIN(6*(J1-1)+2)
               C(3,J1)=FIN(6*(J1-1)+3)
               C(1,J1+NRES)=FIN(6*(J1-1)+4)
               C(2,J1+NRES)=FIN(6*(J1-1)+5)
               C(3,J1+NRES)=FIN(6*(J1-1)+6)
            ENDDO
            CALL UPDATEDC
            CALL INT_FROM_CART(.TRUE.,.FALSE.)
            CALL GEOM_TO_VAR(NINTS,FIN(1:NINTS))
C           CALL UNGUESSPATH(Q,FIN,NINTS,EDIFFTOL,NATOMS)
            CALL GUESSPATH(Q,FIN,NINTS,EDIFFTOL,NATOMS)
         ELSE
            CALL GUESSPATH(Q,FIN,NOPT,EDIFFTOL,NATOMS)
         ENDIF
      ELSE IF (CONNECTT.OR.BHINTERPT.OR.BISECTT) THEN
C        IF (TWOENDS) THEN
C           CALL CONNECTTWO
C        ELSE
         IF (NEWCONNECTT.OR.BHINTERPT.OR.BISECTT) THEN
            IF (UNRST) THEN
                DO J1=1,NRES
                   C(1,J1)=Q(6*(J1-1)+1)
                   C(2,J1)=Q(6*(J1-1)+2)
                   C(3,J1)=Q(6*(J1-1)+3)
                   C(1,J1+NRES)=Q(6*(J1-1)+4)
                   C(2,J1+NRES)=Q(6*(J1-1)+5)
                   C(3,J1+NRES)=Q(6*(J1-1)+6)
                ENDDO
                CALL UPDATEDC
                CALL INT_FROM_CART(.TRUE.,.FALSE.)
                CALL CHAINBUILD
            ENDIF
            IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
! hk286
            IF (RIGIDINIT) THEN
! hk286               CALL CHARMMDUMP(Q,"test1.pdb",.FALSE.)
               ATOMRIGIDCOORDT = .FALSE.
               CALL TRANSFORMCTORIGID (Q, XRIGIDCOORDS)
               CALL POTENTIAL(XRIGIDCOORDS,EINITIAL,LGDUMMY,.TRUE.,.FALSE.,RMSINITIAL,.FALSE.,.FALSE.)
               ATOMRIGIDCOORDT = .TRUE.
            ELSE
               CALL POTENTIAL(Q,EINITIAL,LGDUMMY,.TRUE.,.FALSE.,RMSINITIAL,.FALSE.,.FALSE.)
            ENDIF
            WRITE(*,'(a,2(g20.10,a))') ' OPTIM> Initial energy=',EINITIAL,' RMS force=',RMSinitial
            IF (UNRST) THEN
                DO J1=1,NRES
                   C(1,J1)=FIN(6*(J1-1)+1)
                   C(2,J1)=FIN(6*(J1-1)+2)
                   C(3,J1)=FIN(6*(J1-1)+3)
                   C(1,J1+NRES)=FIN(6*(J1-1)+4)
                   C(2,J1+NRES)=FIN(6*(J1-1)+5)
                   C(3,J1+NRES)=FIN(6*(J1-1)+6)
                ENDDO
                CALL UPDATEDC
                CALL INT_FROM_CART(.TRUE.,.FALSE.)
                CALL CHAINBUILD
            ENDIF
            IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
! hk286
            IF (RIGIDINIT) THEN
               ATOMRIGIDCOORDT = .FALSE.
               CALL TRANSFORMCTORIGID (FIN, XRIGIDCOORDS)
               CALL POTENTIAL(XRIGIDCOORDS,EFINAL,LGDUMMY,.TRUE.,.FALSE.,RMSFINAL,.FALSE.,.FALSE.)
               ATOMRIGIDCOORDT = .TRUE.
            ELSE
               CALL POTENTIAL(FIN,EFINAL,LGDUMMY,.TRUE.,.FALSE.,RMSFINAL,.FALSE.,.FALSE.)
            ENDIF
            WRITE(*,'(a,2(g20.10,a))') ' OPTIM> Final energy  =',  EFINAL,' RMS force=',RMSfinal
            IF (MAX(RMSINITIAL,RMSFINAL)>GMAX.AND.(BFGSMINT.OR.BSMIN.OR.RKMIN)) THEN ! SAT
               WRITE(*,'(A)') 'OPTIM> Bad endpoints supplied - RMS force too big!'
               WRITE(*,'(A,G20.10)') 'OPTIM> Acceptable RMS force would be less or equal to ',GMAX
               IF (REOPTIMISEENDPOINTS) THEN
                  INRSAVE=INR
                  INR=0
                  IF (RMSINITIAL>GMAX.AND.(BFGSMINT.OR.BSMIN.OR.RKMIN)) THEN
                     KNOWE=.FALSE.
                     KNOWG=.FALSE.
                     BFGSTSSAVE=BFGSTST
                     BFGSTST=.FALSE.
! hk286
                     IF (RIGIDINIT) THEN
                        ATOMRIGIDCOORDT = .FALSE.
                        CALL TRANSFORMCTORIGID (Q, XRIGIDCOORDS)
                        CALL GEOPT(FNAMEF,EFNAME,XRIGIDCOORDS,VECS,MFLAG,ENERGY,EVALMIN,VNEW)
                        ATOMRIGIDCOORDT = .TRUE.
                        CALL TRANSFORMRIGIDTOC (1, NRIGIDBODY, Q, XRIGIDCOORDS)
                     ELSE
                        CALL GEOPT(FNAMEF,EFNAME,Q,VECS,MFLAG,ENERGY,EVALMIN,VNEW)
                     ENDIF
! hk286
                     BFGSTST=BFGSTSSAVE
                  ENDIF
                  IF (RMSFINAL>GMAX.AND.(BFGSMINT.OR.BSMIN.OR.RKMIN)) THEN
                     KNOWE=.FALSE.
                     KNOWG=.FALSE.
                     BFGSTSSAVE=BFGSTST
                     BFGSTST=.FALSE.
! hk286
                     IF (RIGIDINIT) THEN
                        ATOMRIGIDCOORDT = .FALSE.
                        CALL TRANSFORMCTORIGID (FIN, XRIGIDCOORDS)
                        CALL GEOPT(FNAMEF,EFNAME,XRIGIDCOORDS,VECS,MFLAG,ENERGY,EVALMIN,VNEW)
                        ATOMRIGIDCOORDT = .TRUE.
                        CALL TRANSFORMRIGIDTOC (1, NRIGIDBODY, FIN, XRIGIDCOORDS)
                     ELSE
                        CALL GEOPT(FNAMEF,EFNAME,FIN,VECS,MFLAG,ENERGY,EVALMIN,VNEW)
                     ENDIF
                     BFGSTST=BFGSTSSAVE
                  ENDIF
                  KNOWE=.FALSE.
                  KNOWG=.FALSE.
                  KNOWH=.FALSE.
                  INR=INRSAVE
                  REOPTIMISEENDPOINTS=.FALSE.
                  GOTO 555
               ELSE
                  CALL TSUMMARY
                  CALL FLUSH(6)
                  GOTO 765
               ENDIF
            ELSE IF (MAX(RMSINITIAL,RMSFINAL)>CONVR) THEN ! SAT
               PRINT '(A)', 'OPTIM> Bad endpoints supplied - RMS force too big!'
               PRINT '(A,G20.10)', ' OPTIM> Acceptable RMS force would be less or equal to ',CONVR
               IF (REOPTIMISEENDPOINTS) THEN
                  INRSAVE=INR
                  BFGSTSSAVE=BFGSTST
                  BFGSTST=.FALSE.
                  INR=0
                  IF (RMSINITIAL>CONVR) THEN
                     CALL GEOPT(FNAMEF,EFNAME,Q,VECS,MFLAG,ENERGY,EVALMIN,VNEW)
                  ENDIF
                  IF (RMSFINAL>CONVR) THEN
                     CALL GEOPT(FNAMEF,EFNAME,FIN,VECS,MFLAG,ENERGY,EVALMIN,VNEW)
                  ENDIF
                  INR=INRSAVE
                  BFGSTST=BFGSTSSAVE
                  REOPTIMISEENDPOINTS=.FALSE.
                  GOTO 555
               ELSE
                  CALL TSUMMARY
                  CALL FLUSH(6)
                  GOTO 765
               ENDIF
            ENDIF
            IF (AMBER12T.AND.CHIRALENDPOINTS) THEN
               CALL CHIRALITY_CHECK(FIN,GOODENDPOINT)
               IF (.NOT.GOODENDPOINT) THEN
                  PRINT '(A)', ' OPTIM> Chirality change in endpoints'
                  STOP
               ELSE
                  WRITE(*,'(A)') ' tryconnect> finish endpoint passed chirality check'
               ENDIF
            ENDIF
            IF (AMBER12T.AND.NOCISTRANS) THEN ! check the finish geometry for cis/trans wrt start
               CALL CIS_TRANS_CHECK(FIN,GOODENDPOINT)
               IF (.NOT.GOODENDPOINT) THEN
                  WRITE(*,'(A)') 
     &     ' tryconnect> Cis-trans isomerisation of a peptide bond detected in finish end point (wrt. the original structure), quit' 
                  STOP
               ELSE
                  WRITE(*,'(A)') ' tryconnect> finish endpoint cis/trans consistent with starting geometry'
               ENDIF
            ENDIF
            CALL ALIGN_DECIDE(Q,FIN,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
            IF (PERMDISTINIT) PERMDIST=.FALSE.
            IF (ATOMMATCHINIT) ATOMMATCHDIST=.FALSE.
            IF (BISECTT) THEN
               CALL BISECT_OPT(NATOMS,EINITIAL,Q,EFINAL,FIN,DIST)
            ELSE
               IF (ALLOCATED(SAVES)) DEALLOCATE(SAVES)
               IF (ALLOCATED(SAVEF)) DEALLOCATE(SAVEF)
               ALLOCATE(SAVES(NOPT),SAVEF(NOPT))
               SAVES(1:NOPT)=Q(1:NOPT)
               SAVEF(1:NOPT)=FIN(1:NOPT)
               CALL NEWCONNECT(NATOMS,EINITIAL,Q,EFINAL,FIN,DIST,.TRUE.,REDOPATH,REDOPATHXYZ,NOPT)
               IF (ALLOCATED(SAVES)) DEALLOCATE(SAVES)
               IF (ALLOCATED(SAVEF)) DEALLOCATE(SAVEF)
            ENDIF
         ELSE
            CALL CONNECT(NCDONE,Q)
         ENDIF
         IF (CALCRATES) CALL RATES(NATOMS,NINTS) ! JMC
      ELSE IF (MECCANOT) THEN
         IF (UNRST) THEN
            DO J1=1,NRES
               C(1,J1)=Q(6*(J1-1)+1)
               C(2,J1)=Q(6*(J1-1)+2)
               C(3,J1)=Q(6*(J1-1)+3)
               C(1,J1+NRES)=Q(6*(J1-1)+4)
               C(2,J1+NRES)=Q(6*(J1-1)+5)
               C(3,J1+NRES)=Q(6*(J1-1)+6)
            ENDDO
            CALL UPDATEDC
            CALL INT_FROM_CART(.TRUE.,.FALSE.)
            CALL GEOM_TO_VAR(NINTS,Q(1:NINTS))
            CALL CHAINBUILD
            CALL POTENTIAL(Q,E1,LGDUMMY,.TRUE.,.FALSE.,RMSINITIAL,.FALSE.,.FALSE.)
            DO J1=1,NRES
               C(1,J1)=FIN(6*(J1-1)+1)
               C(2,J1)=FIN(6*(J1-1)+2)
               C(3,J1)=FIN(6*(J1-1)+3)
               C(1,J1+NRES)=FIN(6*(J1-1)+4)
               C(2,J1+NRES)=FIN(6*(J1-1)+5)
               C(3,J1+NRES)=FIN(6*(J1-1)+6)
            ENDDO
            CALL UPDATEDC
            CALL INT_FROM_CART(.TRUE.,.FALSE.)
            CALL GEOM_TO_VAR(NINTS,FIN(1:NINTS))
            CALL CHAINBUILD
            CALL POTENTIAL(FIN,E2,LGDUMMY,.TRUE.,.FALSE.,RMSFINAL,.FALSE.,.FALSE.)
            IF (DEBUG) WRITE(*,'(A,F20.10,A,F15.6,A,F20.10,A,F15.6)') 
     1           ' OPTIM> Initial energy=',E1,' RMS force=',RMSINITIAL,' final energy=',E2,' RMS=',RMSFINAL
            CALL UNMECCANO(.FALSE.,.TRUE.,ENERGY,.FALSE.,Q,FIN,E1,E2,RMSINITIAL,RMSFINAL)
         ELSE
            IF (.NOT.VARIABLES) CALL NEWMINDIST(Q,FIN,NATOMS,DIST,BULKT,TWOD,ZSYM(1),.FALSE.,RIGIDBODY,DEBUG,RMAT)
            CALL POTENTIAL(Q,E1,LGDUMMY,.TRUE.,.FALSE.,RMSINITIAL,.FALSE.,.FALSE.)
            CALL POTENTIAL(FIN,E2,LGDUMMY,.TRUE.,.FALSE.,RMSFINAL,.FALSE.,.FALSE.)
            IF (DEBUG) WRITE(*,'(A,F20.10,A,F15.6,A,F20.10,A,F15.6)') 
            DIST=0.0D0
            DO J1=1,NOPT
               DIST=DIST+(Q(J1)-FIN(J1))**2
            ENDDO
            DIST=SQRT(DIST)
            NNNIMAGE=NINT(MIN(MECIMDENS*DIST,MECMAXIMAGES*1.0D0)) ! IMAGE DENSITY TIMES DISTANCE
            IF (NNNIMAGE < 1       ) NNNIMAGE=1
            NITERMAX=NINT(MIN(NNNIMAGE*MECITDENS,MECMAXIT*1.0D0)) ! NUMBER OF IMAGES TIMES ITERATION DENSITY
            CALL MECCANO(.FALSE.,.TRUE.,ENERGY,.FALSE.,Q,FIN,E1,E2,RMSINITIAL,RMSFINAL)
         ENDIF
      ELSE IF (NEWNEBT.OR.NEBT.OR.GROWSTRINGT) THEN
         IF (NEBMIND) THEN
            CALL ALIGN_DECIDE(Q,FIN,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
            REALSTR=WR(DIST,3)
            WRITE(*,'(a,f12.2)')' OPTIM> Structures were put in the closest coincidence, distance = '//trim(RealStr)
         ENDIF
         IF (UNRST) THEN
             DO J1=1,NRES
                C(1,J1)=Q(6*(J1-1)+1)
                C(2,J1)=Q(6*(J1-1)+2)
                C(3,J1)=Q(6*(J1-1)+3)
                C(1,J1+NRES)=Q(6*(J1-1)+4)
                C(2,J1+NRES)=Q(6*(J1-1)+5)
                C(3,J1+NRES)=Q(6*(J1-1)+6)
             ENDDO
             CALL UPDATEDC
             CALL INT_FROM_CART(.TRUE.,.FALSE.)
             CALL CHAINBUILD
         ENDIF
         IF(CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
         CALL POTENTIAL(Q,EINITIAL,LGDUMMY,.TRUE.,.FALSE.,RMSINITIAL,.FALSE.,.FALSE.)
         WRITE(*,'(a,2(g20.10,a))') ' OPTIM> Initial energy=',EINITIAL,' RMS force=',RMSinitial
         IF (UNRST) THEN
             DO J1=1,NRES
                C(1,J1)=FIN(6*(J1-1)+1)
                C(2,J1)=FIN(6*(J1-1)+2)
                C(3,J1)=FIN(6*(J1-1)+3)
                C(1,J1+NRES)=FIN(6*(J1-1)+4)
                C(2,J1+NRES)=FIN(6*(J1-1)+5)
                C(3,J1+NRES)=FIN(6*(J1-1)+6)
             ENDDO
             CALL UPDATEDC
             CALL INT_FROM_CART(.TRUE.,.FALSE.)
             CALL CHAINBUILD
         ENDIF
         IF (CHRMMT.AND.ACESOLV) NCHENCALLS=ACEUPSTEP-1
         CALL POTENTIAL(FIN,EFINAL,LGDUMMY,.TRUE.,.FALSE.,RMSFINAL,.FALSE.,.FALSE.)
         WRITE(*,'(a,2(g20.10,a))') ' OPTIM> Final energy  =',  EFINAL,' RMS force=',RMSfinal
         IF (NEWNEBT.OR.GROWSTRINGT) THEN
              CALL NEWNEB(.FALSE.,DCOORDS,EINITIAL,Q,EFINAL,FIN,.TRUE.)
         ELSE
              CALL OLDNEB(.TRUE.,.TRUE.,ENERGY,MAXNEBBFGS,.FALSE.,Q)
         ENDIF
      ELSE IF (PATHT) THEN
!
! New block to do ts search first, then do the path without needed BFGSSTEP
!
         IF ((BFGSTST.OR.(INR.EQ.2).OR.(INR.EQ.1)).AND.(.NOT.BFGSSTEP)) THEN
            IF (RIGIDINIT) THEN
               ATOMRIGIDCOORDT = .FALSE.
               CALL TRANSFORMCTORIGID (Q, XRIGIDCOORDS)
               CALL GEOPT(FNAMEF,EFNAME,XRIGIDCOORDS,VECS,MFLAG,ENERGY,EVALMIN,VNEW)
               ATOMRIGIDCOORDT = .TRUE.
               CALL TRANSFORMRIGIDTOC (1, NRIGIDBODY, Q, XRIGIDCOORDS)
            ELSE
               CALL GEOPT(FNAMEF,EFNAME,Q,VECS,MFLAG,ENERGY,EVALMIN,VNEW)
            ENDIF
            IF (.NOT.MFLAG) GOTO 765 ! skip path for failed ts search!
            POTCALL=.FALSE.
            KNOWVECS=.TRUE.
!
! If no other minimisation condition is set then use eigenvector-following for the
! paths via INR=0.
! To use Page-McIver quadratic SD use the PMPATH keyword with argument 6 (default) or
! 7 (can stop at a saddle point). 
!
            IF (INR.EQ.2) INR=0
            IF (INR.EQ.1) INR=0
            IF (PMPATHT) INR=PMPATHINR
         ELSE
!
! Initialise VECS in case this doesn't happen elsewhere
!
            IF (UNRST) THEN
               DO J1=1,NINTS
                  VECS(J1)=DPRAND()*2-1.0D0
               ENDDO
               CALL VECNORM(VECS,NINTS)
            ELSE
               DO J1=1,NOPT
                  VECS(J1)=DPRAND()*2-1.0D0
               ENDDO
               CALL VECNORM(VECS,NOPT)
            ENDIF
            POTCALL=.TRUE.
         ENDIF

         IF(ORDERPARAMT.OR.RINGPOLYMERT) QSAVE(1:NOPT)=Q(1:NOPT)
         IF (FILTH.EQ.0) THEN
            EOFSSTRING='EofS'
            ITSTRING='points.path.xyz'
         ELSE
            WRITE(EOFSSTRING,'(A)') 'EofS.'//TRIM(ADJUSTL(FILTHSTR))
            WRITE(ITSTRING,'(A)') 'points.path.xyz.'//TRIM(ADJUSTL(FILTHSTR))
         ENDIF

         DO J1=1,NOPT
            FRQSTS(J1)=0.0D0
            FRQSMINUS(J1)=0.0D0
            FRQSPLUS(J1)=0.0D0
         ENDDO
         IF (RIGIDINIT) THEN
            ATOMRIGIDCOORDT = .FALSE.
            CALL TRANSFORMCTORIGID (Q, XRIGIDCOORDS)
            CALL PATH(XRIGIDCOORDS,ENERGY,VNEW,RMS,EVALMIN,VECS,POTCALL,QPLUS,QMINUS,.TRUE.,ETS,EPLUS,EMINUS,
     1             SLENGTH,DISP,GAMMA,NTILDE,FRQSTS,FRQSPLUS,FRQSMINUS,ITSTRING,EOFSSTRING,PATHFAILT) 
              ! jmc added arg PATHFAILT
            ATOMRIGIDCOORDT = .TRUE.
            CALL TRANSFORMRIGIDTOC (1, NRIGIDBODY, Q, XRIGIDCOORDS)
         ELSE
            CALL PATH(Q,ENERGY,VNEW,RMS,EVALMIN,VECS,POTCALL,QPLUS,QMINUS,.TRUE.,ETS,EPLUS,EMINUS,
     1             SLENGTH,DISP,GAMMA,NTILDE,FRQSTS,FRQSPLUS,FRQSMINUS,ITSTRING,EOFSSTRING,PATHFAILT) 
              ! jmc added arg PATHFAILT
         ENDIF

         IF (PATHFAILT) THEN
            CALL FLUSH(6)
            GOTO 765
         ENDIF

         IF (CALCRATES) CALL RATES(NATOMS,NINTS) ! JMC

         IF (ORDERPARAMT) THEN
            Q(1:NOPT)=QSAVE(1:NOPT)
            KNOWE=.FALSE.
            KNOWG=.FALSE.
            KNOWH=.FALSE.
!bs360
            CALL GEOPT(FNAMEF,EFNAME,Q,VECS,MFLAG,ENERGY,EVALMIN,VNEW)
         ENDIF
         IF (RINGPOLYMERT.AND.(ENDHESS.OR.ENDNUMHESS)) THEN 
            Q(1:NOPT)=QSAVE(1:NOPT)
            NSTEPS=0 ! will avoid geometry optimisation - need ENDHESS or ENDNUMHESS
            IF (.NOT.(ENDHESS.OR.ENDNUMHESS)) THEN
               PRINT '(A)',' OPTIM> Neither ENDHESS nor ENDNUMHESS keywords are set - making ENDHESS true'
            ENDIF
! to get the instanton rates.
            CALL GEOPT(FNAMEF,EFNAME,Q,VECS,MFLAG,ENERGY,EVALMIN,VNEW) 
         ENDIF



C     ELSE IF (DRAGT) THEN
C        CALL DRAG
      ELSE IF (CLOSESTALIGNMENT) THEN
         CALL ALIGN_DECIDE(Q,FIN,NATOMS,DEBUG,PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
         WRITE(*,'(a,F10.5)') ' OPTIM> Minimized distance between aligned structures = ',DIST
         IF (CHRMMT) THEN 
            CALL CHARMMDUMP(FIN,'aligned.crd',MACHINE)
         ELSE IF (AMBERT .OR. NABT .OR. AMBER12T) THEN
            ! Formats used come from the AMBER routine minrit
            OPEN(UNIT=3,FILE='aligned.rst', STATUS='UNKNOWN')
            WRITE(3,'(20A4)') 'MOL'
            IF (NATOMS.GT.100000) THEN
               WRITE(3,'(I6)') NATOMS
            ELSE
               WRITE(3,'(I5)') NATOMS
            ENDIF
            WRITE(3,'(6f12.7)') (FIN(J2),J2=1,NOPT)
            CLOSE(3)
         ELSE
            CALL DUMPIT(FIN,'aligned')
         ENDIF
      ELSE
! hk286
         IF (RIGIDINIT) THEN
            ATOMRIGIDCOORDT = .FALSE.
            CALL TRANSFORMCTORIGID (Q, XRIGIDCOORDS)
            CALL GEOPT(FNAMEF,EFNAME,XRIGIDCOORDS,VECS,MFLAG,ENERGY,EVALMIN,VNEW)
            ATOMRIGIDCOORDT = .TRUE.
            CALL TRANSFORMRIGIDTOC (1, NRIGIDBODY, Q, XRIGIDCOORDS)
         ELSE
            CALL GEOPT(FNAMEF,EFNAME,Q,VECS,MFLAG,ENERGY,EVALMIN,VNEW)
         ENDIF
      ENDIF
      IF (.NOT.MULTIJOBT) THEN
         CLOSE(1)
         CLOSE(2)
         IF (DUMPV) CLOSE(44)
      ENDIF

765   IF (MULTIJOBT) THEN
         IF (INDEX(MULTISTART,'xyz').NE.0) THEN
            READ(MULTISUNIT,'(A)',END=864) DUMMYS
!           PRINT '(2A)','DUMMYS=',DUMMYS
            READ(MULTISUNIT,'(A)',END=864) DUMMYS
!           PRINT '(2A)','DUMMYS=',DUMMYS
            DO J1=1,NATOMS
!              READ(MULTISUNIT,'(A2,1X,3G20.10)',END=864) DUMMYS,Q(3*(J1-1)+1:3*(J1-1)+3)
               IF (VARIABLES) THEN
                  READ(MULTISUNIT,*,END=864) Q(3*(J1-1)+1:3*(J1-1)+3)
               ELSE
                  READ(MULTISUNIT,*,END=864) DUMMYS,Q(3*(J1-1)+1:3*(J1-1)+3)
               ENDIF
            ENDDO
         ELSE IF (MULTIJOB_MACHINET) THEN

            IF (MULTI_COUNT .GT. MULTI_LAST) THEN
               WRITE(*,*) "optim> Reached MULTIJOB_MACHINET configuration ", MULTI_COUNT
               WRITE(*,*) "optim> Past last configuration specified by MULTIJOB_MACHINE:", MULTI_LAST
               GOTO 864         ! Exit the loop over input coordinates
            ENDIF

            IF(DEBUG) WRITE(*,*) "Reading new configuration ", MULTI_COUNT, "from file ", TRIM(ADJUSTL(MULTISTART))
            OPEN(MULTISUNIT,FILE=TRIM(ADJUSTL(MULTISTART)),ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='OLD',RECL=8*NOPT)
            READ(MULTISUNIT, REC=MULTI_COUNT) Q(1:NOPT)
            CLOSE(MULTISUNIT)

            IF(MULTIFUNIT.GT.0) THEN
               OPEN(MULTISUNIT,FILE=TRIM(ADJUSTL(MULTISTART)),ACCESS='DIRECT',FORM='UNFORMATTED',STATUS='OLD',RECL=8*NOPT)
               READ(MULTISUNIT, REC=MULTI_COUNT) FIN(1:NOPT)
               CLOSE(MULTIFUNIT)
            ENDIF

            MULTI_COUNT = MULTI_COUNT+MULTI_STEP

         ELSE
            READ(MULTISUNIT,*,END=864) Q(1:NOPT)
         ENDIF
         PRINT '(A)',' OPTIM> Start coordinates overwritten from file '//TRIM(ADJUSTL(MULTISTART))
!         PRINT '(3G20.10)',Q(1:NOPT)
         IF (GTHOMSONT) THEN
! Convert Cartesian input coordiantes to polars
            CALL GTHOMSONCTOANG(Q(1:3*NATOMS), THTEMP(1:3*NATOMS), NATOMS, 6)
            Q(1:3*NATOMS)=THTEMP(1:3*NATOMS)
         ENDIF
         IF (MULTIFUNIT.GT.0) THEN
            IF (INDEX(MULTIFINISH,'xyz').NE.0) THEN
               READ(MULTIFUNIT,'(A2)',END=864) DUMMYS
               READ(MULTIFUNIT,'(A2)',END=864) DUMMYS
               DO J1=1,NATOMS
                  IF (VARIABLES) THEN
                     READ(MULTIFUNIT,*,END=864) FIN(3*(J1-1)+1:3*(J1-1)+3)
                  ELSE
                     READ(MULTIFUNIT,'(A2,3G20.10)',END=864) DUMMYS,FIN(3*(J1-1)+1:3*(J1-1)+3)
                  ENDIF
               ENDDO
            ELSE
               READ(MULTIFUNIT,*,END=864) FIN(1:NOPT)
            ENDIF
            PRINT '(A)',' OPTIM> Finish coordinates overwritten from file '//TRIM(ADJUSTL(MULTIFINISH))
         IF (GTHOMSONT) THEN 
! Convert Cartesian input coordiantes to polars
            CALL GTHOMSONCTOANG(FIN(1:3*NATOMS), THTEMP(1:3*NATOMS), NATOMS, 6)
            FIN(1:3*NATOMS)=THTEMP(1:3*NATOMS)
         ENDIF

         ENDIF
         KNOWE=.FALSE. ; KNOWG=.FALSE. ; KNOWH=.FALSE.
         IF (DUMPPATH.OR.DUMPALLPATHS) THEN
            IF (FILTH.EQ.0) THEN
               WRITE(PINFOSTRING,'(A9)') 'path.info'
            ELSE
               WRITE(PINFOSTRING,'(A)') 'path.info.'//TRIM(ADJUSTL(FILTHSTR))
            ENDIF

            IF (MACHINE) THEN
                 OPEN(UNIT=88,FILE=PINFOSTRING,STATUS='UNKNOWN',FORM='UNFORMATTED',POSITION='APPEND')
            ELSE
                 OPEN(UNIT=88,FILE=PINFOSTRING,STATUS='UNKNOWN',POSITION='APPEND')
            ENDIF
         ENDIF
         INR=MULTIINR
         GOTO 963
864      CONTINUE
         PRINT '(A)',' OPTIM> End of multijob coordinates'
      ENDIF
      CALL TSUMMARY

      IF (INTMINT.OR.DESMINT.OR.INTINTERPT.OR.NATINT) CALL INTCLEANUP
      IF (UNRST) DEALLOCATE(UREFCOORD,UREFPPSANGLE,INTSTEP)
      IF (ZSYM(NATOMS).EQ.'SV') call cleanMemory

      IF (OPEPT) CALL OPEP_FINISH

      ENDIF!Call to Dl-find for quantumrate stuff

      IF (ALLOCATED(STPMAX)) DEALLOCATE(STPMAX)
      IF (ALLOCATED(NR)) DEALLOCATE(NR)
      IF (ALLOCATED(ZSYM)) DEALLOCATE(ZSYM)
      IF (ALLOCATED(FROZEN)) DEALLOCATE(FROZEN)
      IF (ALLOCATED(TAGNUM)) DEALLOCATE(TAGNUM)
      IF (ALLOCATED(TAGFAC)) DEALLOCATE(TAGFAC)
      IF (ALLOCATED(DUMPMODEN)) DEALLOCATE(DUMPMODEN)
      IF (ALLOCATED(BESTPERM)) DEALLOCATE(BESTPERM)
      IF (ALLOCATED(ATMASS)) DEALLOCATE(ATMASS)
      IF (ALLOCATED(FROZENRES)) DEALLOCATE(FROZENRES)
      IF (ALLOCATED(IATNUM)) DEALLOCATE(IATNUM)
      IF (ALLOCATED(ATMASS)) DEALLOCATE(ATMASS)
      IF (ALLOCATED(REPELTS)) DEALLOCATE(REPELTS) 
      IF (ALLOCATED(VECCHK)) DEALLOCATE(VECCHK)
      IF (ALLOCATED(ZWORK)) DEALLOCATE(ZWORK)
      IF (ALLOCATED(HESS)) DEALLOCATE(HESS)
      IF (ALLOCATED(FIN)) DEALLOCATE(FIN)
      IF (ALLOCATED(START)) DEALLOCATE(START)
      IF (ALLOCATED(QW)) DEALLOCATE(QW)
      IF (ALLOCATED(MIEF_SIG)) DEALLOCATE(MIEF_SIG)
      IF (ALLOCATED(MIEF_EPS)) DEALLOCATE(MIEF_EPS)
      IF (ALLOCATED(MIEF_SITES)) DEALLOCATE(MIEF_SITES)
      IF (ALLOCATED(MIEF_U_RCUT)) DEALLOCATE(MIEF_U_RCUT)
      IF (ALLOCATED(MIEF_DUDR_RCUT)) DEALLOCATE(MIEF_DUDR_RCUT)
      IF (ALLOCATED(MLPDAT)) DEALLOCATE(MLPDAT) 
      IF (ALLOCATED(MLPOUTCOME)) DEALLOCATE(MLPOUTCOME) 
      IF (ALLOCATED(SHIFTL)) DEALLOCATE(SHIFTL) 
      
      IF (DUMPDATAT .OR. DUMPDATA_MACHINET) CLOSE(881)
      IF (DUMPDATA_MACHINET) CLOSE(882)     
      IF (DUMPMAG) CLOSE(321)
      CALL FLUSH(6) 
      
      RETURN

      END

      SUBROUTINE TSUMMARY
      USE KEY,ONLY : TSTART
      IMPLICIT NONE
      INTEGER NPCALL, ECALL, FCALL, SCALL
      DOUBLE PRECISION ETIME, FTIME, STIME, TFINISH
      COMMON /PCALL/ NPCALL, ECALL, FCALL, SCALL, ETIME, FTIME, STIME

      CALL MYCPU_TIME(TFINISH,.FALSE.)
      TFINISH=TFINISH-TSTART

      PRINT*
      WRITE(*,'(A,F15.2)') ' Elapsed time=                      ',TFINISH
      IF (TFINISH.NE.0.0D0) THEN
         WRITE(*,'(A,I10,A,F15.2,A,F5.1)') ' OPTIM> # of energy calls=                 ',ECALL,' time=',ETIME,
     1         ' %=',ETIME*100.0D0/TFINISH
         WRITE(*,'(A,I10,A,F15.2,A,F5.1)') ' OPTIM> # of energy+gradient calls=        ',FCALL,' time=',FTIME,
     1         ' %=',FTIME*100.0D0/TFINISH
         WRITE(*,'(A,I10,A,F15.2,A,F5.1)') ' OPTIM> # of energy+gradient+Hessian calls=',SCALL,' time=',STIME,
     1         ' %=',STIME*100/TFINISH
      ENDIF


      RETURN
      END

      SUBROUTINE INITIALBOX(QC)
      USE KEY, ONLY: BULK_BOXVEC
      USE COMMONS, ONLY: NATOMS
      IMPLICIT NONE
      DOUBLE PRECISION QC(3*NATOMS)
      INTEGER J1
      
      DO J1=1, NATOMS
         QC(3*(J1-1)+1)=QC(3*(J1-1)+1)-BULK_BOXVEC(1)*NINT(QC(3*(J1-1)+1)/BULK_BOXVEC(1))
         QC(3*(J1-1)+2)=QC(3*(J1-1)+2)-BULK_BOXVEC(2)*NINT(QC(3*(J1-1)+2)/BULK_BOXVEC(2))
         QC(3*(J1-1)+3)=QC(3*(J1-1)+3)-BULK_BOXVEC(3)*NINT(QC(3*(J1-1)+3)/BULK_BOXVEC(3))
      ENDDO

      END
