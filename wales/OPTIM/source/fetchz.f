C   OPTIM: A programfor optimizing geometries and calculating reaction pathways
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
      SUBROUTINE FETCHZ(Q)
      USE COMMONS
      USE KEY
      USE MODTWOEND
      USE SYMINF
      USE MODNEB
      USE modcharmm
      USE modamber9, only : KTWNT
      USE MODUNRES
      use KeyConnect
      USE MODGUESS
      USE MODMEC
      use PORFUNCS
      USE msevb_common
      use BINARYIO
      USE SDWATER, ONLY : SDINIT
      USE BOWMANWATER, ONLY : BOWMANINIT
      USE AMHGLOBALS
      USE GENRIGID, ONLY: RIGIDINIT
      USE GSDATA, ONLY : GSUPDATE

      IMPLICIT NONE
      INTEGER J1, J, NN, MM, NTYPE(105), J2, IDUM1, IDUM2, LNATOMS, NUSE, J3, J4, NDUMMY, H, NTYPES
      LOGICAL YESNO
      DOUBLE PRECISION EPSAB, EPSBB, SIGAB, SIGBB, EPS, CSC, SIG, AMAT(3,3), A_TEMP(3,3) 
      DOUBLE PRECISION TEMPX, TEMPY, TEMPZ, AINV(3,3), DETER, Q(3*NATOMS)
      INTEGER NTYPEA, NELEMENTS, NOXYGEN, NCHARGE
      INTEGER, ALLOCATABLE :: ELEMENT_NUMS(:)
      COMMON /BIN/ EPSAB, EPSBB, SIGAB, SIGBB, NTYPEA
      CHARACTER FNAME*80, TSTRING*80, CD1*1, CD2*2
      COMMON /CAS/ AMAT, AINV, NELEMENTS, NTYPE
      LOGICAL CONNECTT, DUMPPATH, READPATH, CALCRATES, STOPFIRST
      INTEGER NCONNECT,I500,LUNIT,GETUNIT
      DOUBLE PRECISION TEMPERATURE, HRED, SCALE_FACT
      COMMON /CONN/ STOPFIRST, CONNECTT, NCONNECT, DUMPPATH, READPATH, CALCRATES, TEMPERATURE, HRED
      LOGICAL PATHT, DRAGT
      INTEGER NPATHFRAME
      ! sn402: The name for this block appeared to be missing, so I've added it.
      COMMON /RUNTYPE/ DRAGT, PATHT, NPATHFRAME
C     DOUBLE PRECISION REPELTS(3*NATOMS,100), REPELPUSH
C     INTEGER NREPELTS, REPELFROM
C     LOGICAL REPELTST, REPEL
C     COMMON /OTS/ NREPELTS, REPELTST, REPELPUSH, REPEL, REPELFROM
      LOGICAL CUBIC
      COMMON /CUB/ CUBIC
      DOUBLE PRECISION CAPSRHO, EPS2, RAD, HEIGHT
      COMMON /CAPS/ CAPSRHO, EPS2, RAD, HEIGHT
      character(len=20) :: FTEMP
      character(len=20) :: FINSTRING
      CHARACTER(LEN=21) :: WW
      CHARACTER(LEN=80) :: MYLINE
      CHARACTER(5),ALLOCATABLE :: ELEMENT_NAMES(:) 
      INTEGER ITEM, NITEMS, LOC, LINE, NCR, NERROR, IR, LAST
      LOGICAL END, SKIPBL, CLEAR, ECHO, CAT
      COMMON /BUFINF/ ITEM, NITEMS, LOC(132), LINE, SKIPBL, CLEAR, NCR, NERROR, IR, ECHO, LAST, CAT

      CHARACTER(LEN=13) :: CUDAFILENAME1
      CHARACTER(LEN=34) :: CUDAFILENAME2
      CHARACTER(LEN=33) :: CUDAFILENAME3
      CHARACTER(LEN=20) :: CUDAFILENAME4
      CHARACTER(LEN=23) :: CUDAFILENAME5
      CHARACTER(LEN=17) :: CUDAFILENAME6

C    LOCAL AMH VARIABLES
      INTEGER GLY_COUNT, III
      DOUBLE PRECISION X, Y, Z

C
C  Thomson problem:
C
      DOUBLE PRECISION, PARAMETER :: HALFPI=1.570796327D0
      DOUBLE PRECISION THTEMP(3*NATOMS), DIST
      DOUBLE PRECISION, ALLOCATABLE :: FOCK(:,:), MOCOEFF(:,:)

      IF (CHRMMT.OR.UNRST.OR.AMBERT.OR.NABT.OR.AMBER12T.OR.OPEPT) THEN
         IF (.NOT.ALLOCATED(IATNUM)) ALLOCATE (IATNUM(NATOMS))   ! ATMASS already set up
      ELSE IF (PYT.OR.PYGPERIODICT.OR.PYBINARYT.OR.PYGT.OR.MULTISITEPYT) THEN
         IF (.NOT.ALLOCATED(IATNUM)) ALLOCATE (IATNUM(NATOMS),ATMASS(NATOMS))
         ATMASS(:)=1.0D0
      ELSE
         IF (.NOT.ALLOCATED(IATNUM)) ALLOCATE (IATNUM(NATOMS))
         IF (.NOT.ALLOCATED(ATMASS)) ALLOCATE (ATMASS(NATOMS))
C        sf344> initialise ATMASS values to not default to zero 
C               uninitialised stuff (may cause problems in INERTIA otherwise)
         ! With RIGIDINIT, masses should now be specified by a separate file and will be initialised
         ! in GENRIGID_INITIALISE.
         IF (.NOT. RIGIDINIT) ATMASS(:)=1.0D0
      ENDIF
      IF (MULTISITEPYT) CALL DEFINEPYMULTISITES
C
C  Remind users in case they have not changed the tighter default setting of MAXERISE for
C  PV, VASP, QCHEM, CASTEP, ONETEP and CP2K runs. 
C
      IF (PV.AND.(MAXERISE.LT.1.0D-7)) WRITE(*,'(A,G20.10)') 
     1          ' fetchz> WARNING - maximum permitted energy rise seems a bit small: ',MAXERISE
      IF ((VASP.OR.QCHEM.OR.CASTEP.OR.ONETEP.OR.CP2K).AND.(MAXERISE.LT.1.0D-7)) WRITE(*,'(A,G20.10)') 
     1          ' fetchz> WARNING - maximum permitted energy rise seems a bit small: ',MAXERISE 

      if (machine) then
          if (filth2==0) then
               finstring='points2.inp'
          else
               WRITE(FTEMP,*) FILTH
               WRITE(FINSTRING,'(A)') 'points2.inp.' // TRIM(ADJUSTL(FTEMP))
          endif
      else
           IF (FILTH2.EQ.0) THEN
              FINSTRING='finish'
           ELSE
              WRITE(FTEMP,*) FILTH
              WRITE(FINSTRING,'(A)') 'finish.' // TRIM(ADJUSTL(FTEMP))
           ENDIF
      endif

      IF (HYBRIDMINT) THEN
         PRINT '(3(A,I5),A)',' fetchz> Hybrid minimisation, maximum steps=',HMNSTEPS,
     &      ', maximum tangent space steps=',HMNBFGSMAX1,' or ',HMNBFGSMAX2,' at convergence'
         IF (NUSEEV.GT.0) WRITE(*,'(A,I5,A,G20.10)') '         Maximum number of eigenvalues/eigenvectors to be calculated=',NUSEEV,
     &                       ' maximum step size along an eigenvector=',HMMXSTP
         IF (NOHESS) THEN
            WRITE(*,'(A,I4,A,F12.4)') 
     1       '         No Hessian: using Rayleight-Ritz; allowed steps=',HMNEVS,' convergence for RMS < ',HMCEIG
         WRITE(*,'(A,A2,A,G20.10)') '         step type is ',HMMETHOD,
     &                ' and will revert to LBFGS if the lowest nonzero eigenvalue > ',HMEVMAX
         ELSE
            WRITE(*,'(A,I4,A,I4,A,F12.4,A)') 
     1        '         Allowed steps for largest Hessian eigenvector=',HMNEVL,
     2        ', smallest eigenvector=',HMNEVS,' convergence at ',HMCEIG,'%'
         ENDIF
      ENDIF

      IF (BFGSMINT.AND.(.NOT.BFGSTST)) THEN
         WRITE(*,'(A,I8)') ' fetchz> LBFGS minimization, maximum steps=',BFGSSTEPS
      ELSE IF (BFGSTST) THEN
         IF (BFGSSTEP.AND.(.NOT.HYBRIDMINT)) THEN
            WRITE(*,'(A)') ' fetchz> One EF step will be taken before BFGS minimization'
         ELSE
            WRITE(*,'(A,I6)') ' fetchz> Hybrid EF/BFGS transition state search, maximum steps=',NSTEPS
            WRITE(*,'(A,I4,A,I4,A,F12.6)') '         maximum tangent space steps=',NBFGSMAX1, 
     &                               ' or ',NBFGSMAX2,' when overlap is better than ',1.0D0-BFGSTSTOL
            IF (MODEDOWNT) THEN
               WRITE(*,'(A,I4,A,I4,A)') ' fetchz> Downhill mode is ',IVEC,' for initial step and ',IVEC2,' after that'
            ELSE
               WRITE(*,'(A,I4,A,I4,A)') ' fetchz> Uphill mode is ',IVEC,' for initial step and ',IVEC2,' after that'
            ENDIF
         ENDIF
         IF (.NOT.(READV.AND.BFGSSTEP)) THEN
            IF (NOHESS) THEN
               WRITE(*,'(A,I4,A,F12.4,A,G12.2)') 
     1          ' fetchz> No Hessian: Rayleigh-Ritz max steps=',NEVS,
     &          ' RMS < ',CEIG,' and |% change|<',EVPC 
            ELSE
               WRITE(*,'(A,I4,A,I4,A,F12.4,A)') 
     1           ' fetchz> Steps for largest Hessian e/vector=',NEVL,
     2           ', smallest e/vector=',NEVS,' convergence at ',CEIG,'%'
            ENDIF
         ENDIF
         IF ((NSECDIAG.NE.1).AND.(NSECDIAG.NE.2).AND.(NSECDIAG.NE.3)) THEN
            PRINT '(A)',' fetchz> unrecognised value for NSECDIAG - resetting to one'
            NSECDIAG=1
         ENDIF
         IF (.NOT.NOIT) WRITE(*,'(A,I4)') 
     &         ' fetchz> Hessian eigenvector in Rayleigh-Ritz scheme calculated using method ',NSECDIAG
      ELSE IF (BBRSDMT) THEN
         WRITE(*,'(A,I5,A,G20.10)')   ' fetchz> BBR steepest-descent method, maximum steps=',BBRSTEPS,' RMS convergence ',BBRCONV
         WRITE(*,'(5(A,F15.5),A,I5)') '         gamma=',BBRGAM,' eps=',BBREPS,' sig1=',BBRSIGMA1,' sig2=',BBRSIGMA2, 
     &                                ' alpha=',BBRALPHA,' m=',BBRM
      ELSE IF (BSMIN) THEN
         NUSE=NSTEPS
         IF (PATHSDSTEPS.GT.0) NUSE=PATHSDSTEPS
         WRITE(*,'(A,I5)') ' fetchz> Bulirsch-Stoer gradient only steepest-descent, maximum steps=',NUSE
      ELSE IF (RKMIN) THEN
         NUSE=NSTEPS
         IF (PATHSDSTEPS.GT.0) NUSE=PATHSDSTEPS
         WRITE(*,'(A,I5)') ' fetchz> 5th order Runge-Kutta gradient only steepest-descent, maximum steps=',NUSE
      ELSE
         IF (NUSEEV.GT.0) WRITE(*,'(A,I5)') ' fetchz> Maximum number of eigenvalues/eigenvectors to be calculated=',NUSEEV
         IF (INR.EQ.0) WRITE(*,'(A,I5)') ' fetchz> Eigenvector-following minimization, maximum steps=',NSTEPS
         IF (INR.EQ.1) WRITE(*,'(A,I5)') ' fetchz> Modified Newton-Raphson steps, maximum steps=',NSTEPS
         IF (INR.EQ.2) WRITE(*,'(A,I5)') ' fetchz> Eigenvector-following transition state search, maximum steps=',NSTEPS
         IF (INR.EQ.3) WRITE(*,'(A,I5)') 
     1' fetchz> Eigenvector-following minimization, no reorientation and pseudo-third derivative correction, maximum steps=',NSTEPS
         IF (INR.EQ.4) WRITE(*,'(A,I5)') 
     1    ' fetchz> Eigenvector-following transition state search with pseudo-third derivative correction, maximum steps=',NSTEPS
         IF (INR.EQ.5) WRITE(*,'(A,I5)') 
     1    ' fetchz> Eigenvector-following minimization, geometry shifted to principal axes at first step, maximum steps=',NSTEPS
         IF (INR.EQ.6) WRITE(*,'(A,I5)') 
     1    ' fetchz> Steepest-descent minimization using first and second derivatives, maximum steps=',NSTEPS
         IF (INR.EQ.7) WRITE(*,'(A,I5)') 
     1    ' fetchz> Steepest-descent minimization using first and second derivatives - can stop at a saddle, maximum steps=',NSTEPS
         IF ((INR.EQ.2).OR.(INR.EQ.4)) THEN
            IF (MODEDOWNT) THEN
               WRITE(*,'(A,I4,A,I4,A)') ' fetchz> Downhill mode is ',IVEC,' for initial step and ',IVEC2,' after that'
            ELSE
               WRITE(*,'(A,I4,A,I4,A)') ' fetchz> Uphill mode is ',IVEC,' for initial step and ',IVEC2,' after that'
            ENDIF
         ENDIF
         IF ((PATHSDSTEPS.GT.0).AND.(INR.LE.6)) PRINT '(A,I8)',
     &                  ' fetchz> Maximum number of steps for pathway minimisation=',PATHSDSTEPS
      ENDIF

      IF (GRADSQ) WRITE(*,'(A)') ' fetchz> Using the modulus gradient as the objective function'
      PRINT*
      CALL FLUSH(6)

      IF (AMBERT.OR.NABT.OR.AMBER12T.OR.OPEPT) THEN
         NOPT=3*NATOMS
         WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' AMBER atoms'
         IF (TWOENDS.OR.CONNECTT.OR.NEBT.OR.NEWNEBT.OR.DRAGT.OR.GUESSPATHT
     $     .OR.MECCANOT.OR.MORPHT.OR.GREATCIRCLET.OR.GROWSTRINGT.OR.BHINTERPT.OR.BISECTT.OR.CLOSESTALIGNMENT) THEN
            OPEN (UNIT=7,FILE=FINSTRING,STATUS='OLD')
          IF(AMBERT.OR.NABT.OR.AMBER12T.OR.OPEPT) THEN      ! read coordinates from file finish (containing only coordinates)
            DO J1=1,NATOMS
               READ(7,*)  FIN(3*(J1-1)+1),FIN(3*(J1-1)+2),FIN(3*(J1-1)+3)
            ENDDO
            CLOSE(7)
          ELSE
            DO J1=1,NATOMS
               READ(7,*) CD1, CD2, IDUM1, IDUM2, FIN(3*(J1-1)+1),FIN(3*(J1-1)+2),FIN(3*(J1-1)+3)
C              WRITE(*,*) CD1, CD2, IDUM1, IDUM2, FIN(3*(J1-1)+1),FIN(3*(J1-1)+2),FIN(3*(J1-1)+3)
            ENDDO
          END IF
         ENDIF

      ELSEIF (RINGPOLYMERT) THEN
         NOPT=0
160      CALL INPUT(END)
         IF (.NOT.END) THEN
            DO J1=1,NITEMS
               NOPT=NOPT+1
               CALL READF(Q(NOPT))
            ENDDO
            GOTO 160
         ENDIF
         RPDOF=NOPT/RPIMAGES
         IF (RPDOF*RPIMAGES.NE.NOPT) THEN
            PRINT '(2(A,I6))','fetchz> ERROR *** number of variables read, ',
     &                         NOPT,' is not a multiple of the number of beads ',RPIMAGES
            STOP
         ENDIF
         WRITE(*,'(A,I6,A,I6,A)')' SYSTEM ',NOPT,' coordinates will be optimised for a ring polymer of ',RPIMAGES,' beads'
         IF (RPSYSTEM(1:4).EQ.'AECK') THEN
            WRITE(*,'(A,G20.10,A)') '               1/kT (a.u.)=',RPBETA,' potential type=' // TRIM(ADJUSTL(RPSYSTEM))
            NATOMS=NOPT
         ELSEIF (RPSYSTEM(1:2).EQ.'SD') THEN
            NOXYGEN=(RPDOF/3+1)/3
            NCHARGE=(RPDOF/3)-3*NOXYGEN
            WRITE(*,'(A,G20.10,A)') '               1/kT (mol/kcal)=',RPBETA,' potential type=' // TRIM(ADJUSTL(RPSYSTEM))
            WRITE(*,'(A,2I6)')     '               Number of O and H atoms=',NOXYGEN,(RPDOF/3)-NOXYGEN
            CALL SDINIT(NOXYGEN,NCHARGE)
            NATOMS=NOPT/3
         ELSEIF (RPSYSTEM(1:2).EQ.'TT') THEN
            WRITE(*,'(A,G20.10,A)')  '              1/kT a.u.=',RPBETA,' potential type=' // TRIM(ADJUSTL(RPSYSTEM))
            RPBETA=RPBETA*0.00159360144367 ! convert to kcal/mol units
            WRITE(*,'(A,G20.10,A)')  '              1/kT (mol/kcal)=',RPBETA
            WRITE(*,'(A,I0,A)') '              TTM3-F potential for a cluster of ', RPDOF/9, ' water molecules'
            NATOMS=NOPT/3
         ELSEIF (RPSYSTEM(1:2).EQ.'MB') THEN
            WRITE(*,'(A,G20.10,A)')  '              1/kT a.u.=',RPBETA,' potential type=' // TRIM(ADJUSTL(RPSYSTEM))
            WRITE(*,'(A,G20.10,A)')  '              1/kT (mol/kcal)=',RPBETA
            WRITE(*,'(A,I0,A)') '              MBPOL potential for a cluster of ', RPDOF/9, ' water molecules'
            NATOMS=NOPT/3
         ELSEIF (RPSYSTEM(1:3).EQ.'MAL') THEN
            WRITE(*,'(A,G20.10,A)')  '              1/kT a.u.=',RPBETA,' potential type=' // TRIM(ADJUSTL(RPSYSTEM))
            WRITE(*,'(A,I0,A)') '              malonaldehyde potential '
            NATOMS=NOPT/3
         ELSEIF (RPSYSTEM(1:3).EQ.'MCY') THEN
            WRITE(*,'(A,G20.10,A)')  '              1/kT (a.u.)=',RPBETA,' potential type=' // TRIM(ADJUSTL(RPSYSTEM))
            WRITE(*,'(A,I1,A,I0,A)') '              VRT(MCY-5f) flexible water dimer potential'
            NOXYGEN=RPDOF/9
            IF (NOXYGEN /= 2) THEN
               PRINT *, 'ERROR: VRT(MCY-5f) only for water dimer; NOXYGEN = ', NOXYGEN
               STOP
            ENDIF
            NATOMS=NOPT/3
         ELSEIF (RPSYSTEM(1:2).EQ.'JB') THEN
            NOXYGEN=RPDOF/9
            WRITE(*,'(A,G20.10,A)')  '              1/kT (a.u.)=',RPBETA,' potential type=' // TRIM(ADJUSTL(RPSYSTEM))
            WRITE(*,'(A,I1,A,I0,A)') '              Bowman''s PES#',BOWMANPES,' for a cluster of ',NOXYGEN,' water molecules'
            CALL BOWMANINIT(NOXYGEN, BOWMANPES, trim(BOWMANDIR))
            NATOMS=NOPT/3
         ELSE
            PRINT *, ' fetchz> unknown RPSYSTEM: ', RPSYSTEM
            STOP
         ENDIF
         WRITE(*,'(A,I6,A)')     '               Reading particle masses for ',RPDOF,' degrees of freedom from file RPmasses:'
         ALLOCATE(RPMASSES(RPDOF))
!
! Must not use unit 1 here - the points and energies files are connected to
! units 1 and 2 if PRINTPTS or MORPHT is true.
!
         IF (RPFIXT) THEN
            ALLOCATE(XMINA(RPDOF), XMINB(RPDOF))
            OPEN(UNIT=121,FILE='xminA',STATUS='OLD')
            READ(121,*) XMINA
            CLOSE(121)
            OPEN(UNIT=122,FILE='xminB',STATUS='OLD')
            READ(122,*) XMINB
            CLOSE(122)
         ENDIF
         OPEN(UNIT=123,FILE='RPmasses',STATUS='OLD')
         READ(123,*) RPMASSES(1:RPDOF)
         CLOSE(123)
         WRITE(*,'(6G20.10)') RPMASSES(1:RPDOF)
         IF (RPSYSTEM(1:4).EQ.'AECK') THEN
            DO J1=1,NOPT
               IATNUM(J1)=1
               ZSYM(J1)='RP'
            ENDDO
         ELSE IF (RPSYSTEM.EQ.'SD') THEN
            DO J1=1,RPIMAGES
               DO J2=1,NOXYGEN
                  IATNUM((RPDOF/3)*(J1-1)+J2)=16
                  ZSYM((RPDOF/3)*(J1-1)+J2)='O '
               ENDDO
               DO J2=NOXYGEN+1,RPDOF/3
                  IATNUM((RPDOF/3)*(J1-1)+J2)=1
                  ZSYM((RPDOF/3)*(J1-1)+J2)='H '
               ENDDO
            ENDDO
         ELSE IF (RPSYSTEM.EQ.'MAL') THEN
            DO J1=1,RPIMAGES
               IATNUM((RPDOF/3)*(J1-1)+1)=12
               ZSYM((RPDOF/3)*(J1-1)+1)='C '
               IATNUM((RPDOF/3)*(J1-1)+2)=12
               ZSYM((RPDOF/3)*(J1-1)+2)='C '
               IATNUM((RPDOF/3)*(J1-1)+3)=16
               ZSYM((RPDOF/3)*(J1-1)+3)='O '
               IATNUM((RPDOF/3)*(J1-1)+4)=12
               ZSYM((RPDOF/3)*(J1-1)+4)='C '
               IATNUM((RPDOF/3)*(J1-1)+5)=16
               ZSYM((RPDOF/3)*(J1-1)+5)='O '
               IATNUM((RPDOF/3)*(J1-1)+6)=1
               ZSYM((RPDOF/3)*(J1-1)+6)='H '
               IATNUM((RPDOF/3)*(J1-1)+7)=1
               ZSYM((RPDOF/3)*(J1-1)+7)='H '
               IATNUM((RPDOF/3)*(J1-1)+8)=1
               ZSYM((RPDOF/3)*(J1-1)+8)='H '
               IATNUM((RPDOF/3)*(J1-1)+9)=1
               ZSYM((RPDOF/3)*(J1-1)+9)='H '
            ENDDO
         ELSE IF (RPSYSTEM.EQ.'JB') THEN
            DO J1=1,RPIMAGES
               DO J2=1,2*NOXYGEN
                  IATNUM((RPDOF/3)*(J1-1)+J2)=1
                  ZSYM((RPDOF/3)*(J1-1)+J2)='H '
               ENDDO
               DO J2=2*NOXYGEN+1,RPDOF/3
                  IATNUM((RPDOF/3)*(J1-1)+J2)=16
                  ZSYM((RPDOF/3)*(J1-1)+J2)='O '
               ENDDO
            ENDDO
         ENDIF
      ELSE IF (AMHT) THEN
         NOPT=3*NATOMS
         WRITE(*,'(A,I4,A,I4,A)') ' SYSTEM ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' AMH atoms'

         IF (TWOENDS.OR.CONNECTT.OR.NEBT.OR.NEWNEBT.OR.DRAGT.OR.GUESSPATHT.OR.MECCANOT.OR.CLOSESTALIGNMENT) THEN
              IF (FILTH2.EQ.0) THEN
                 FINSTRING='finish'
              ELSE
                 WRITE(FTEMP,*) FILTH
                 WRITE(FINSTRING,'(A)') 'finish.' // TRIM(ADJUSTL(FTEMP))
              ENDIF

           OPEN(UNIT=80,FILE=FINSTRING,STATUS='OLD',FORM='FORMATTED')

                  WRITE(6,*)'FINISH ', NRES_AMH_TEMP
                  WRITE(6,201)(SEQ(j),j=1,NRES_AMH_TEMP)
201               format (25(i2,1x))
            DO 500 I500=1,NRES_AMH_TEMP
                  READ(80,*)X,Y,Z

76                format(3(G25.15))
C                  WRITE(6,76) X,Y,Z
                  FINCORD(I500,1,1,1)=X
                  FINCORD(I500,2,1,1)=Y
                  FINCORD(I500,3,1,1)=Z

                  READ(80,*)X,Y,Z
c                  WRITE(6,76)X,Y,Z
                  FINCORD(I500,1,1,2)=X
                  FINCORD(I500,2,1,2)=Y
                  FINCORD(I500,3,1,2)=Z

                  READ(80,*)X,Y,Z
c                  WRITE(6,76)X,Y,Z
                  FINCORD(I500,1,1,3)=X
                  FINCORD(I500,2,1,3)=Y
                  FINCORD(I500,3,1,3)=Z
500              CONTINUE

       GLY_COUNT = 0
       DO 764 III = 1,NRES_AMH_TEMP ! AMH MCP HACK

       IF (SEQ(III).EQ.8) THEN
         FIN(9*(III-1)+1- GLY_COUNT*3) = FINCORD(III, 1, 1, 1) !  CA X
         FIN(9*(III-1)+2- GLY_COUNT*3) = FINCORD(III, 2, 1, 1) !  CA Y
         FIN(9*(III-1)+3- GLY_COUNT*3) = FINCORD(III, 3, 1, 1) !  CA Z
         FIN(9*(III-1)+4- GLY_COUNT*3) = FINCORD(III, 1, 1, 3) !  O X
         FIN(9*(III-1)+5- GLY_COUNT*3) = FINCORD(III, 2, 1, 3) !  O Y
         FIN(9*(III-1)+6- GLY_COUNT*3) = FINCORD(III, 3, 1, 3) !  O Z
         GLY_COUNT = GLY_COUNT +1
       ELSE
         FIN(9*(III-1)+1 - GLY_COUNT*3) = FINCORD(III, 1, 1, 1) !  CA X
         FIN(9*(III-1)+2 - GLY_COUNT*3) = FINCORD(III, 2, 1, 1) !  CA Y
         FIN(9*(III-1)+3 - GLY_COUNT*3) = FINCORD(III, 3, 1, 1) !  CA Z
         FIN(9*(III-1)+4 - GLY_COUNT*3) = FINCORD(III, 1, 1, 2) !  CB X
         FIN(9*(III-1)+5 - GLY_COUNT*3) = FINCORD(III, 2, 1, 2) !  CB Y
         FIN(9*(III-1)+6 - GLY_COUNT*3) = FINCORD(III, 3, 1, 2) !  CB Z
         FIN(9*(III-1)+7 - GLY_COUNT*3) = FINCORD(III, 1, 1, 3) !  O X
         FIN(9*(III-1)+8 - GLY_COUNT*3) = FINCORD(III, 2, 1, 3) !  O Y
         FIN(9*(III-1)+9 - GLY_COUNT*3) = FINCORD(III, 3, 1, 3) !  O Z
       ENDIF

764    CONTINUE

         ENDIF

        DO J1=1,NATOMS
               IATNUM(J1)=1
               ZSYM(J1)='AM'
        ENDDO

      ELSE IF (CHRMMT) THEN
         NOPT=3*NATOMS
         IF (.NOT.INTMINT) THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' CHARMM atoms'
         ELSE
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> internal coordinates will be optimised for ',NATOMS,' CHARMM atoms'
         ENDIF
         IF (TWOENDS.OR.CONNECTT.OR.NEBT.OR.NEWNEBT.OR.DRAGT.OR.GUESSPATHT
     $     .OR.MECCANOT.OR.MORPHT.OR.GREATCIRCLET.OR.GROWSTRINGT.OR.BHINTERPT.OR.BISECTT.OR.CLOSESTALIGNMENT) THEN
            IF (MACHINE) then
                 OPEN(7,FILE=FINSTRING,access='direct',form='unformatted',status='old',recl=3*8*Natoms)
                 read(7,rec=1) (FIN(j1),j1=1,3*Natoms)
                 CLOSE(7)
            else
                 OPEN (UNIT=7,FILE=FINSTRING,STATUS='OLD')
C csw34> OPTIM can now cope with CHARMM comment lines (*) in the finish file
            DO
               READ(7,*) MYLINE
               IF (MYLINE(1:1)=='*') THEN
                  CYCLE
               ELSE
                  EXIT
               ENDIF
            ENDDO

                 DO J1=1,NATOMS
!                   READ (7,*) IDUM1,IDUM2,CD1,CD2,FIN(3*(J1-1)+1),FIN(3*(J1-1)+2),FIN(3*(J1-1)+3),IDUM1,IDUM2
                    READ (7,*) IDUM1,IDUM2,CD1,CD2,FIN(3*(J1-1)+1),FIN(3*(J1-1)+2),FIN(3*(J1-1)+3)
                 ENDDO
            endif
         ENDIF
      ELSE IF (UNRST) THEN
         NOPT=3*NATOMS
C jmc now set in keywords.f         NINTS=2*nres-5+2*nside
         IF (.NOT.INTMINT) THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' UNRES atoms'
         ELSE
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NINTS,' internal coordinates will be optimised for ',NATOMS,' UNRES atoms'
         ENDIF
         IF (TWOENDS.OR.CONNECTT.OR.NEBT.OR.NEWNEBT.OR.DRAGT.OR.GUESSPATHT
     $     .OR.MECCANOT.OR.MORPHT.OR.GREATCIRCLET.OR.GROWSTRINGT.OR.BHINTERPT.OR.BISECTT.OR.CLOSESTALIGNMENT) THEN
            OPEN (UNIT=7,FILE=FINSTRING,STATUS='OLD')
            READ (7,*) 
            DO J1=1,NATOMS
               READ (7,*) CD1,FIN(3*(J1-1)+1),FIN(3*(J1-1)+2),FIN(3*(J1-1)+3)
            ENDDO
         ENDIF
!
! We need to get the coordinates into vector Q and atomic symbols into ZSYM
!
!IEU_RAPH
      ELSE IF (VASP) THEN
        WRITE(*,'(A,I4,A)') ' fetchz> VASP run for ',NATOMS,' atoms'
        WRITE(*,'(A,A)') ' fetchz> VASP run command is ',TRIM(ADJUSTL(VASPJOB))

	CALL SYSTEM('grep "TITEL" POTCAR | wc > temp_NTYPES') !Searches POTCAR for 'TITEL', to calculate how many element types (NTYPES) are present)
	FNAME='temp_NTYPES'
	OPEN(UNIT=7, FILE=FNAME, STATUS='OLD')
	READ (7,*) NTYPES
	CLOSE (7)
	ALLOCATE(ELEMENT_NAMES(NTYPES), ELEMENT_NUMS(NTYPES))

        !Unscaled unit cell matrix (A_TEMP) is read from the POSCAR file and multiplied by a scale factor (SCALE_FACT)
        FNAME='POSCAR'
        OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
        READ(7,*)                            !First line is a comment in POSCAR FILE
        READ (7,*) SCALE_FACT
        READ(7,*) A_TEMP(1,1), A_TEMP(2,1), A_TEMP(3,1)
        READ(7,*) A_TEMP(1,2), A_TEMP(2,2), A_TEMP(3,2)
        READ(7,*) A_TEMP(1,3), A_TEMP(2,3), A_TEMP(3,3)

        !Scaled unit cell matrix written to AMAT
        AMAT(1,1) = SCALE_FACT*A_TEMP(1,1)
        AMAT(2,1) = SCALE_FACT*A_TEMP(2,1)
        AMAT(3,1) = SCALE_FACT*A_TEMP(3,1)
        AMAT(1,2) = SCALE_FACT*A_TEMP(1,2)
        AMAT(2,2) = SCALE_FACT*A_TEMP(2,2)
        AMAT(3,2) = SCALE_FACT*A_TEMP(3,2)
        AMAT(1,3) = SCALE_FACT*A_TEMP(1,3)
        AMAT(2,3) = SCALE_FACT*A_TEMP(2,3)
        AMAT(3,3) = SCALE_FACT*A_TEMP(3,3)
        WRITE(*,'(3F20.10)') ((AMAT(J1,J2),J1=1,3),J2=1,3)
        
        !In VASP element names (ELEMENT_NAMES) are given on line 6 of POSCAR ie Li Co O
	!Line 7 gives the number of each element present (ELEMENT_NUMS) ie 2 2 4 
        READ(7,*) (ELEMENT_NAMES(J1), J1=1,NTYPES)
        READ(7,*) (ELEMENT_NUMS(J1), J1=1,NTYPES)
        NATOMS=SUM(ELEMENT_NUMS)                                   

         DETER=-AMAT(1,3)*AMAT(2,2)*AMAT(3,1)
     1         +AMAT(1,2)*AMAT(2,3)*AMAT(3,1)
     2         +AMAT(1,3)*AMAT(2,1)*AMAT(3,2)
     3         -AMAT(1,1)*AMAT(2,3)*AMAT(3,2)
     4         -AMAT(1,2)*AMAT(2,1)*AMAT(3,3)
     5         +AMAT(1,1)*AMAT(2,2)*AMAT(3,3)


         AINV(1,1)=(-AMAT(2,3)*AMAT(3,2)+AMAT(2,2)*AMAT(3,3))/DETER
         AINV(1,2)=( AMAT(1,3)*AMAT(3,2)-AMAT(1,2)*AMAT(3,3))/DETER
         AINV(1,3)=(-AMAT(1,3)*AMAT(2,2)+AMAT(1,2)*AMAT(2,3))/DETER
         AINV(2,1)=( AMAT(2,3)*AMAT(3,1)-AMAT(2,1)*AMAT(3,3))/DETER
         AINV(2,2)=(-AMAT(1,3)*AMAT(3,1)+AMAT(1,1)*AMAT(3,3))/DETER
         AINV(2,3)=( AMAT(1,3)*AMAT(2,1)-AMAT(1,1)*AMAT(2,3))/DETER
         AINV(3,1)=(-AMAT(2,2)*AMAT(3,1)+AMAT(2,1)*AMAT(3,2))/DETER
         AINV(3,2)=( AMAT(1,2)*AMAT(3,1)-AMAT(1,1)*AMAT(3,2))/DETER
         AINV(3,3)=(-AMAT(1,2)*AMAT(2,1)+AMAT(1,1)*AMAT(2,2))/DETER
         
	WRITE(*,'(A,A)') ' fetchz> Reading coordinates from file ',FNAME

        
	H=0
        DO J1=1, NTYPES
            DO J2=H+1,H+ELEMENT_NUMS(J1)
               ZSYM(J2)=ELEMENT_NAMES(J1) !Atomic symbols are read into ZSYM 
               H=SUM(ELEMENT_NUMS(1:J1))
            ENDDO
        ENDDO
	READ(7,*) FNAME !On line 8 of POSCAR 'Direct' keyword read specifying direct(fractional) coordiantes 
	
        !Coordiantes read into vector Q and transformed into cartesian
        DO J1=1, NATOMS        
            READ(7,*) Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3)
            TEMPX=AMAT(1,1)*Q(3*(J1-1)+1)+AMAT(1,2)*Q(3*(J1-1)+2)+AMAT(1,3)*Q(3*(J1-1)+3)
            TEMPY=AMAT(2,1)*Q(3*(J1-1)+1)+AMAT(2,2)*Q(3*(J1-1)+2)+AMAT(2,3)*Q(3*(J1-1)+3)
            TEMPZ=AMAT(3,1)*Q(3*(J1-1)+1)+AMAT(3,2)*Q(3*(J1-1)+2)+AMAT(3,3)*Q(3*(J1-1)+3)
            Q(3*(J1-1)+1)=TEMPX
            Q(3*(J1-1)+2)=TEMPY
            Q(3*(J1-1)+3)=TEMPZ
            WRITE(*,*) Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3)
	ENDDO
	CLOSE(7)
        CALL SYSTEM(' rm temp_NTYPES ')
        NOPT=3*NATOMS !Calculates the number coordinates, used in RMS calculation
!IEU_RAPH
         IF (TWOENDS.OR.CONNECTT.OR.NEBT.OR.NEWNEBT.OR.DRAGT.OR.GUESSPATHT
     $     .OR.MECCANOT.OR.MORPHT.OR.GREATCIRCLET.OR.GROWSTRINGT.OR.BHINTERPT.OR.BISECTT.OR.CLOSESTALIGNMENT) THEN
            LUNIT=GETUNIT()
            OPEN (LUNIT,FILE=FINSTRING,STATUS='OLD')
            DO J1=1,NATOMS
               READ(LUNIT,*) FIN(3*(J1-1)+1),FIN(3*(J1-1)+2),FIN(3*(J1-1)+3)
            ENDDO
            CLOSE(LUNIT)
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> Coordinates of second point read from file finish'
         ENDIF

      ELSE IF (QCHEM) THEN
         IF (QCHEMES) THEN
!
! This block to read the Fock matrix and AO coefficients in the MOs should be the same
! as in potential.f. It requires a qchem.out file to exist with the right format!
!
            NOPT=QCHEMESNAO*QCHEMESNELEC
            PRINT '(A,3I10)',' fetchz> Number of aos, electrons, and variables=',QCHEMESNAO,QCHEMESNELEC,NOPT
            CALL SYSTEM(' sed -e "1,/Calculated Alpha Fock Matrix/d" qchem.out > edited.1')
            LUNIT=GETUNIT()
            OPEN(LUNIT,FILE='edited.1',STATUS='OLD')
            READ(LUNIT,*) FNAME
            IF (DEBUG) PRINT '(A)',' potential> Reading Fock matrix:'
            IF (DEBUG) PRINT '(A)', FNAME
            ALLOCATE(FOCK(QCHEMESNAO,QCHEMESNAO))
            J4=0
            DO J1=1,QCHEMESNAO/6
               DO J2=1,QCHEMESNAO
                  READ(LUNIT,*) NDUMMY,(FOCK(6*(J1-1)+J3,J2),J3=1,6)
                  WRITE(*,'(I5,6F12.7)') NDUMMY,(FOCK(6*(J1-1)+J3,J2),J3=1,6)
               ENDDO
               J4=J4+6
            ENDDO
            IF (J4.LT.QCHEMESNAO) THEN
               DO J2=1,QCHEMESNAO
                  READ(LUNIT,*) NDUMMY,(FOCK(J3,J2),J3=J4+1,QCHEMESNAO)
                  WRITE(*,'(I5,6F12.7)') NDUMMY,(FOCK(J3,J2),J3=J4+1,QCHEMESNAO)
               ENDDO
            ENDIF
            CLOSE(LUNIT)

            CALL SYSTEM(' sed -e "1,/Alpha MO Coefficients/d" qchem.out > edited.1')
                        LUNIT=GETUNIT()
            OPEN(LUNIT,FILE='edited.1',STATUS='OLD')
            IF (DEBUG) PRINT '(A)',' potential> Reading MO coefficients:'
            READ(LUNIT,*) FNAME
            ALLOCATE(MOCOEFF(QCHEMESNELEC,QCHEMESNAO))
            J4=0
            DO J1=1,QCHEMESNELEC/6
               DO J2=1,QCHEMESNAO
                  READ(LUNIT,*) NDUMMY,(FOCK(6*(J1-1)+J3,J2),J3=1,6)
                  WRITE(*,'(I5,6F12.7)') NDUMMY,(FOCK(6*(J1-1)+J3,J2),J3=1,6)
               ENDDO
               J4=J4+6
            ENDDO
            IF (J4.LT.QCHEMESNELEC) THEN
               DO J2=1,QCHEMESNAO
                  READ(LUNIT,*) NDUMMY,(MOCOEFF(J3,J2),J3=J4+1,QCHEMESNELEC)
                  WRITE(*,'(I5,6F12.7)') NDUMMY,(MOCOEFF(J3,J2),J3=J4+1,QCHEMESNELEC)
               ENDDO
            ENDIF
            CLOSE(LUNIT)
!
! Now we need a matrix multiply to get the actual gradient?
!
            DEALLOCATE(FOCK,MOCOEFF)
!
! Note that we haven't attempted to code reading finish file for double-ended search yet!
!
         ELSE
            FNAME=SYS(1:LSYS)
            WRITE(*,'(A,A)') ' fetchz> Reading coordinates from file ',FNAME
            LUNIT=GETUNIT()
            OPEN(LUNIT,FILE=FNAME,STATUS='OLD')
            readq: DO 
               READ(LUNIT,'(A19)') WW
               CALL UPPERCASE(WW)
               IF (WW(1:4).EQ.'$MOL') THEN
                  READ(LUNIT,*) WW ! charge and multiplicity
                  DO J1=1,NATOMS
                     READ(LUNIT,*) WW(1:2),Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3)
                     CALL UPPERCASE(WW)
                     ZSYM(J1)=WW(1:2)
                  ENDDO
                  EXIT readq
               ENDIF
            ENDDO readq
            CLOSE(LUNIT)
            NOPT=3*NATOMS
         ENDIF
!        WRITE(*,'(A)') 'Coordinates in fetchz'
!        WRITE(*,'(3F15.5)') (Q(J1),J1=1,NOPT)

         IF (TWOENDS.OR.CONNECTT.OR.NEBT.OR.NEWNEBT.OR.DRAGT.OR.GUESSPATHT
     $     .OR.MECCANOT.OR.MORPHT.OR.GREATCIRCLET.OR.GROWSTRINGT.OR.BHINTERPT.OR.BISECTT.OR.CLOSESTALIGNMENT) THEN
            LUNIT=GETUNIT()
            OPEN (LUNIT,FILE=FINSTRING,STATUS='OLD')
            DO J1=1,NATOMS
               READ(LUNIT,*) FIN(3*(J1-1)+1),FIN(3*(J1-1)+2),FIN(3*(J1-1)+3)
            ENDDO
            CLOSE(LUNIT)
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> Coordinates of second point read from file finish'
         ENDIF
         
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc 
c interfacing MOLPRO here
      ELSE IF (MOLPRO) THEN
         FNAME=SYS(1:LSYS)
         WRITE(*,'(A,A)') ' fetchz> Reading coordinates from file ',FNAME
         LUNIT=GETUNIT()
         OPEN(LUNIT,FILE=FNAME,STATUS='OLD')
         readq2: DO 
            READ(LUNIT,'(A19)') WW
            CALL UPPERCASE(WW)
            IF (WW(1:9).EQ.' GEOMETRY') THEN
               DO J1=1,NATOMS
                  READ(LUNIT,*) WW(1:2),Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3)
                  CALL UPPERCASE(WW)
                  ZSYM(J1)=WW(1:2)
               ENDDO
               EXIT readq2
            ENDIF
         ENDDO readq2
         CLOSE(LUNIT)
         NOPT=3*NATOMS         
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



      ELSE IF (CASTEP) THEN
         WRITE(*,'(A,I4,A)') ' fetchz> CASTEP run for ',NATOMS,' atoms'
         WRITE(*,'(A,A)') ' fetchz> CASTEP run command is ',TRIM(ADJUSTL(CASTEPJOB))
         FNAME=SYS(1:LSYS) // '.cell'
         IF (DEBUG) WRITE(*,'(A,A)') ' fetchz> Reading cell parameters from file ',FNAME
         OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         readcell: DO 
            READ(7,'(A19)') WW
            CALL UPPERCASE(WW)
            IF (WW(1:19).EQ.'%BLOCK LATTICE_CART') THEN
               READ(7,*) WW ! the line with the units on it
               READ(7,*) AMAT(1,1), AMAT(2,1), AMAT(3,1)
               READ(7,*) AMAT(1,2), AMAT(2,2), AMAT(3,2)
               READ(7,*) AMAT(1,3), AMAT(2,3), AMAT(3,3)
               IF (DEBUG) PRINT '(A)',' fetchz> cell matrix:'
               IF (DEBUG) PRINT '(3F20.10)',((AMAT(J1,J2),J2=1,3),J1=1,3)
               EXIT readcell
            ENDIF
         ENDDO readcell
         CLOSE(7)

         DETER=-AMAT(1,3)*AMAT(2,2)*AMAT(3,1)
     1         +AMAT(1,2)*AMAT(2,3)*AMAT(3,1)
     2         +AMAT(1,3)*AMAT(2,1)*AMAT(3,2)
     3         -AMAT(1,1)*AMAT(2,3)*AMAT(3,2)
     4         -AMAT(1,2)*AMAT(2,1)*AMAT(3,3)
     5         +AMAT(1,1)*AMAT(2,2)*AMAT(3,3)
         AINV(1,1)=(-AMAT(2,3)*AMAT(3,2)+AMAT(2,2)*AMAT(3,3))/DETER
         AINV(1,2)=( AMAT(1,3)*AMAT(3,2)-AMAT(1,2)*AMAT(3,3))/DETER
         AINV(1,3)=(-AMAT(1,3)*AMAT(2,2)+AMAT(1,2)*AMAT(2,3))/DETER
         AINV(2,1)=( AMAT(2,3)*AMAT(3,1)-AMAT(2,1)*AMAT(3,3))/DETER
         AINV(2,2)=(-AMAT(1,3)*AMAT(3,1)+AMAT(1,1)*AMAT(3,3))/DETER
         AINV(2,3)=( AMAT(1,3)*AMAT(2,1)-AMAT(1,1)*AMAT(2,3))/DETER
         AINV(3,1)=(-AMAT(2,2)*AMAT(3,1)+AMAT(2,1)*AMAT(3,2))/DETER
         AINV(3,2)=( AMAT(1,2)*AMAT(3,1)-AMAT(1,1)*AMAT(3,2))/DETER
         AINV(3,3)=(-AMAT(1,2)*AMAT(2,1)+AMAT(1,1)*AMAT(2,2))/DETER

            FNAME=SYS(1:LSYS) // '.cell'
            WRITE(*,'(A,A)') ' fetchz> Reading coordinates from file ',FNAME
            OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
            readX: DO 
            READ(7,'(A20)') WW
            CALL UPPERCASE(WW)
            IF (WW(1:16).EQ.'%BLOCK POSITIONS') THEN
               IF (WW(18:20).EQ.'ABS') THEN ! dont need to convert from fractional coordinates
                  AMAT(1:3,1:3)=0.0D0; AINV(1:3,1:3)=0.0D0
                  AMAT(1,1)=1.0D0;AMAT(2,2)=1.0D0;AMAT(3,3)=1.0D0
                  AINV(1,1)=1.0D0;AINV(2,2)=1.0D0;AINV(3,3)=1.0D0
               ENDIF
               DO J1=1,NATOMS
                  READ(7,*) WW(1:2),Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3)
                  CALL UPPERCASE(WW)
                  ZSYM(J1)=WW(1:2)
C                 PRINT *,ZSYM(J1)(1:2),Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3)
                  TEMPX=AMAT(1,1)*Q(3*(J1-1)+1)+AMAT(1,2)*Q(3*(J1-1)+2)+AMAT(1,3)*Q(3*(J1-1)+3)
                  TEMPY=AMAT(2,1)*Q(3*(J1-1)+1)+AMAT(2,2)*Q(3*(J1-1)+2)+AMAT(2,3)*Q(3*(J1-1)+3)
                  TEMPZ=AMAT(3,1)*Q(3*(J1-1)+1)+AMAT(3,2)*Q(3*(J1-1)+2)+AMAT(3,3)*Q(3*(J1-1)+3)
                  Q(3*(J1-1)+1)=TEMPX
                  Q(3*(J1-1)+2)=TEMPY
                  Q(3*(J1-1)+3)=TEMPZ
!                 IATNUM(J1)=1
!                 ZSYM(J1)='H'
               ENDDO
               EXIT readX
            ENDIF
         ENDDO readX
         CLOSE(7)

         NOPT=3*NATOMS
C        WRITE(*,'(A)') 'The absolute coordinates from fetchz'
C        WRITE(*,'(6F15.5)') (Q(J1),J1=1,NOPT)

         IF (TWOENDS.OR.CONNECTT.OR.NEBT.OR.NEWNEBT.OR.DRAGT.OR.GUESSPATHT
     $     .OR.MECCANOT.OR.MORPHT.OR.GREATCIRCLET.OR.GROWSTRINGT.OR.BHINTERPT.OR.BISECTT.OR.CLOSESTALIGNMENT) THEN
            OPEN (UNIT=7,FILE=FINSTRING,STATUS='OLD') 
            DO J1=1,NATOMS
               READ(7,*) FIN(3*(J1-1)+1),FIN(3*(J1-1)+2),FIN(3*(J1-1)+3)
            ENDDO
            CLOSE(7)
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> Coordinates of second point read from file finish'
         ENDIF
C        PRINT*,' fetchz> coordinates:'
C        WRITE(*,'(3F20.10)') (Q(J1),J1=1,3*NATOMS)

      ELSE IF (ONETEP) THEN
         WRITE(*,'(A,I4,A)') ' fetchz> ONETEP run for ',NATOMS,' atoms'
         FNAME=SYS(1:LSYS) // '.dat'
         IF (DEBUG) WRITE(*,'(A,A)') ' fetchz> Reading cell parameters from file ',FNAME
         OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         readcell2: DO 
            READ(7,'(A19)') WW
            CALL UPPERCASE(WW)
            IF (WW(1:19).EQ.'%BLOCK LATTICE_CART') THEN
C              READ(7,*) WW ! the line with the units on it ! everything must be in a.u.! This line must not be present!
               READ(7,*) AMAT(1,1), AMAT(1,2), AMAT(1,3)
               READ(7,*) AMAT(2,1), AMAT(2,2), AMAT(2,3)
               READ(7,*) AMAT(3,1), AMAT(3,2), AMAT(3,3)
               IF (DEBUG) PRINT '(A)',' fetchz> cell matrix:'
               IF (DEBUG) PRINT '(3F20.10)',((AMAT(J1,J2),J2=1,3),J1=1,3)
               EXIT readcell2
            ENDIF
         ENDDO readcell2
         CLOSE(7)

         DETER=-AMAT(1,3)*AMAT(2,2)*AMAT(3,1)
     1         +AMAT(1,2)*AMAT(2,3)*AMAT(3,1)
     2         +AMAT(1,3)*AMAT(2,1)*AMAT(3,2)
     3         -AMAT(1,1)*AMAT(2,3)*AMAT(3,2)
     4         -AMAT(1,2)*AMAT(2,1)*AMAT(3,3)
     5         +AMAT(1,1)*AMAT(2,2)*AMAT(3,3)
         AINV(1,1)=(-AMAT(2,3)*AMAT(3,2)+AMAT(2,2)*AMAT(3,3))/DETER
         AINV(1,2)=( AMAT(1,3)*AMAT(3,2)-AMAT(1,2)*AMAT(3,3))/DETER
         AINV(1,3)=(-AMAT(1,3)*AMAT(2,2)+AMAT(1,2)*AMAT(2,3))/DETER
         AINV(2,1)=( AMAT(2,3)*AMAT(3,1)-AMAT(2,1)*AMAT(3,3))/DETER
         AINV(2,2)=(-AMAT(1,3)*AMAT(3,1)+AMAT(1,1)*AMAT(3,3))/DETER
         AINV(2,3)=( AMAT(1,3)*AMAT(2,1)-AMAT(1,1)*AMAT(2,3))/DETER
         AINV(3,1)=(-AMAT(2,2)*AMAT(3,1)+AMAT(2,1)*AMAT(3,2))/DETER
         AINV(3,2)=( AMAT(1,2)*AMAT(3,1)-AMAT(1,1)*AMAT(3,2))/DETER
         AINV(3,3)=(-AMAT(1,2)*AMAT(2,1)+AMAT(1,1)*AMAT(2,2))/DETER

         FNAME=SYS(1:LSYS) // '.dat'
         WRITE(*,'(A,A)') ' fetchz> Reading coordinates from file ',FNAME
         OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         readX2: DO 
            READ(7,'(A20)') WW
            CALL UPPERCASE(WW)
            IF (WW(1:16).EQ.'%BLOCK POSITIONS') THEN
               IF (WW(18:20).EQ.'ABS') THEN ! dont need to convert from fractional coordinates
                  AMAT(1:3,1:3)=0.0D0; AINV(1:3,1:3)=0.0D0
                  AMAT(1,1)=1.0D0;AMAT(2,2)=1.0D0;AMAT(3,3)=1.0D0
                  AINV(1,1)=1.0D0;AINV(2,2)=1.0D0;AINV(3,3)=1.0D0
               ENDIF
               DO J1=1,NATOMS
                  READ(7,*) WW(1:2),Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3)
                  CALL UPPERCASE(WW)
                  ZSYM(J1)=WW(1:2)
C                 PRINT *,WW(1:2),Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3)
                  TEMPX=AMAT(1,1)*Q(3*(J1-1)+1)+AMAT(1,2)*Q(3*(J1-1)+2)+AMAT(1,3)*Q(3*(J1-1)+3)
                  TEMPY=AMAT(2,1)*Q(3*(J1-1)+1)+AMAT(2,2)*Q(3*(J1-1)+2)+AMAT(2,3)*Q(3*(J1-1)+3)
                  TEMPZ=AMAT(3,1)*Q(3*(J1-1)+1)+AMAT(3,2)*Q(3*(J1-1)+2)+AMAT(3,3)*Q(3*(J1-1)+3)
                  Q(3*(J1-1)+1)=TEMPX
                  Q(3*(J1-1)+2)=TEMPY
                  Q(3*(J1-1)+3)=TEMPZ
!                 IATNUM(J1)=6
!                 ZSYM(J1)='C'
               ENDDO
               EXIT readX2
            ENDIF
         ENDDO readX2
         CLOSE(7)

         NOPT=3*NATOMS
C        WRITE(*,'(A)') 'The absolute coordinates from fetchz'
C        WRITE(*,'(6F15.5)') (Q(J1),J1=1,NOPT)

         IF (TWOENDS.OR.CONNECTT.OR.NEBT.OR.NEWNEBT.OR.DRAGT.OR.GUESSPATHT
     $     .OR.MECCANOT.OR.MORPHT.OR.GREATCIRCLET.OR.GROWSTRINGT.OR.BHINTERPT.OR.BISECTT.OR.CLOSESTALIGNMENT) THEN
            OPEN (UNIT=7,FILE=FINSTRING,STATUS='OLD') 
            DO J1=1,NATOMS
               READ(7,*) FIN(3*(J1-1)+1),FIN(3*(J1-1)+2),FIN(3*(J1-1)+3)
            ENDDO
            CLOSE(7)
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> Coordinates of second point read from file finish'
         ENDIF
C        PRINT*,' fetchz> coordinates:'
C        WRITE(*,'(3F20.10)') (Q(J1),J1=1,3*NATOMS)

C     Reading the cell parameters and the atomic coordinates from the CP2K input file 
C     ($FNAME.inp). For simplicity we asume that the cell-vectors are always given in the format
C
C     &CELL
C     A AMAT(1,1), AMAT(1,2), AMAT(1,3) 
C     B AMAT(2,1), AMAT(2,2), AMAT(2,3)
C     C AMAT(3,1), AMAT(3,2), AMAT(3,3)
C     PERIODIC NONE/X/Y/Z/XY/XZ/YX/XYZ 
C     &END CELL
C   
C     The PERIODIC keyword is not necessary but if it is present it should always follow the cell parameter settings
C     If not present full periodicity namely PERIODIC XYZ is assumed.
C     Notice that cp2k lists the three lattice vectors as rows while OPTIM wants them as colums. 
C     For the reading of the atomic coordinates we asume the following format
C
C     &COORD
C     SCALED .TRUE. OR .FALSE.
C     AtomicSymbol1 X1 Y1 Z1 
C     AtomicSymbol2 X2 Y2 Z2 
C     ......................
C     AtomicSymboln Xn Yn Zn 
C     &END COORD
C
C     Notice that the unit for the cell dimensions as well as for the atomic coordinates is Angstrom 
C     ------------------

      ELSE IF (CP2K) THEN
         WRITE(*,'(A,I4,A)') ' fetchz> CP2K run for ',NATOMS,' atoms'
         WRITE(*,'(A,A)') ' fetchz> CP2K run command is ',TRIM(ADJUSTL(CP2KJOB))
         FNAME=SYS(1:LSYS) // '.inp'
         IF (DEBUG) WRITE(*,'(A,A)') ' fetchz> Reading cell parameters from file ',FNAME
         OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         readcell3: DO
            READ(7,'(A21)') WW
            WW=ADJUSTL(WW)
            CALL UPPERCASE(WW)
            IF (WW(1:5).EQ.'&CELL') THEN
               READ(7,*) WW, AMAT(1,1), AMAT(2,1), AMAT(3,1)
               READ(7,*) WW, AMAT(1,2), AMAT(2,2), AMAT(3,2)
               READ(7,*) WW, AMAT(1,3), AMAT(2,3), AMAT(3,3)
               IF (DEBUG) PRINT '(A)',' fetchz> cell matrix:'
               IF (DEBUG) PRINT '(3F20.10)',((AMAT(J1,J2),J2=1,3),J1=1,3)
               EXIT readcell3
            ENDIF
         ENDDO readcell3
         CLOSE(7)

         DETER=-AMAT(1,3)*AMAT(2,2)*AMAT(3,1)
     1         +AMAT(1,2)*AMAT(2,3)*AMAT(3,1)
     2         +AMAT(1,3)*AMAT(2,1)*AMAT(3,2)
     3         -AMAT(1,1)*AMAT(2,3)*AMAT(3,2)
     4         -AMAT(1,2)*AMAT(2,1)*AMAT(3,3)
     5         +AMAT(1,1)*AMAT(2,2)*AMAT(3,3)
         AINV(1,1)=(-AMAT(2,3)*AMAT(3,2)+AMAT(2,2)*AMAT(3,3))/DETER
         AINV(1,2)=( AMAT(1,3)*AMAT(3,2)-AMAT(1,2)*AMAT(3,3))/DETER
         AINV(1,3)=(-AMAT(1,3)*AMAT(2,2)+AMAT(1,2)*AMAT(2,3))/DETER
         AINV(2,1)=( AMAT(2,3)*AMAT(3,1)-AMAT(2,1)*AMAT(3,3))/DETER
         AINV(2,2)=(-AMAT(1,3)*AMAT(3,1)+AMAT(1,1)*AMAT(3,3))/DETER
         AINV(2,3)=( AMAT(1,3)*AMAT(2,1)-AMAT(1,1)*AMAT(2,3))/DETER
         AINV(3,1)=(-AMAT(2,2)*AMAT(3,1)+AMAT(2,1)*AMAT(3,2))/DETER
         AINV(3,2)=( AMAT(1,2)*AMAT(3,1)-AMAT(1,1)*AMAT(3,2))/DETER
         AINV(3,3)=(-AMAT(1,2)*AMAT(2,1)+AMAT(1,1)*AMAT(2,2))/DETER

         FNAME=SYS(1:LSYS) // '.inp'
         WRITE(*,'(A,A)') ' fetchz> Reading coordinates from file ',FNAME
         OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         readX3: DO
            READ(7,'(A21)') WW
            WW=ADJUSTL(WW)
            CALL UPPERCASE(WW)
            IF (WW(1:6).EQ.'&COORD') THEN
               READ(7,'(A21)') WW 
               WW=ADJUSTL(WW)
               CALL UPPERCASE(WW)
               IF (WW(1:6).EQ.'SCALED') THEN 
                  IF (WW(1:14).EQ.'SCALED .FALSE.') THEN ! No conversion needed
                     AMAT(1:3,1:3)=0.0D0; AINV(1:3,1:3)=0.0D0
                     AMAT(1,1)=1.0D0;AMAT(2,2)=1.0D0;AMAT(3,3)=1.0D0
                     AINV(1,1)=1.0D0;AINV(2,2)=1.0D0;AINV(3,3)=1.0D0
                  ENDIF
               ENDIF
               DO J1=1,NATOMS
                  READ(7,*) ZSYM(J1),Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3)
C                 PRINT *,ZSYM(J1),Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3)
                  TEMPX=AMAT(1,1)*Q(3*(J1-1)+1)+AMAT(1,2)*Q(3*(J1-1)+2)+AMAT(1,3)*Q(3*(J1-1)+3)
                  TEMPY=AMAT(2,1)*Q(3*(J1-1)+1)+AMAT(2,2)*Q(3*(J1-1)+2)+AMAT(2,3)*Q(3*(J1-1)+3)
                  TEMPZ=AMAT(3,1)*Q(3*(J1-1)+1)+AMAT(3,2)*Q(3*(J1-1)+2)+AMAT(3,3)*Q(3*(J1-1)+3)
                  Q(3*(J1-1)+1)=TEMPX
                  Q(3*(J1-1)+2)=TEMPY
                  Q(3*(J1-1)+3)=TEMPZ
               ENDDO
               EXIT readX3
            ENDIF
         ENDDO readX3
         CLOSE(7)

         NOPT=3*NATOMS
C        WRITE(*,'(A)') 'The absolute coordinates from fetchz'
C        WRITE(*,'(6F15.5)') (Q(J1),J1=1,NOPT)

         IF (TWOENDS.OR.CONNECTT.OR.NEBT.OR.NEWNEBT.OR.DRAGT.OR.GUESSPATHT
     $     .OR.MECCANOT.OR.MORPHT.OR.GREATCIRCLET.OR.GROWSTRINGT.OR.BHINTERPT.OR.BISECTT.OR.CLOSESTALIGNMENT) THEN
            OPEN (UNIT=7,FILE=FINSTRING,STATUS='OLD') 
            DO J1=1,NATOMS
               READ(7,*) FIN(3*(J1-1)+1),FIN(3*(J1-1)+2),FIN(3*(J1-1)+3)
            ENDDO
            CLOSE(7)
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> Coordinates of second point read from file finish'
         ENDIF
C        PRINT*,' fetchz> coordinates:'
C        WRITE(*,'(3F20.10)') (Q(J1),J1=1,3*NATOMS)

      ELSE IF (CPMD) THEN
         FNAME=SYS(1:LSYS)
         WRITE(*,'(A,A)') ' fetchz> Reading coordinates from file ',FNAME
         OPEN(UNIT=7,FILE=FNAME,STATUS='OLD')
         LNATOMS=0
11       READ(7,'(A)') FNAME
         IF (FNAME(1:6).EQ.'&ATOMS') THEN
            J1=0
12          READ(7,'(A)') TSTRING
            IF (TSTRING(1:1).EQ.'*') THEN
               J1=J1+1
               READ(7,'(A)') FNAME
               READ(7,*) NTYPE(J1)
               DO J2=1,NTYPE(J1)
                  IATNUM(LNATOMS+J2)=1
C                 ZSYM(LNATOMS+J2)=TSTRING(2:3)
                  ZSYM(LNATOMS+J2)='CP'
                  READ(7,*) Q(3*(LNATOMS+J2-1)+1),Q(3*(LNATOMS+J2-1)+2),Q(3*(LNATOMS+J2-1)+3)
               ENDDO
               LNATOMS=LNATOMS+NTYPE(J1)
               GOTO 12
            ELSE IF (TSTRING(1:1).EQ.' ') THEN
               GOTO 12
            ELSE IF (TSTRING(1:4).EQ.'&END') THEN
               GOTO 13
            ENDIF
         ELSE
            GOTO 11
         ENDIF
         IF (PARALLEL) WRITE(*,'(3A)') ' fetchz> Auxilliary program will be run on ' // NPROC // ' processors'

13       CONTINUE
         CLOSE(7)

         NOPT=3*NATOMS
         CALL SYSTEM(' grep -c ANGSTROM ' // SYS(1:LSYS) // ' > temp')
         OPEN(UNIT=7,FILE='temp',STATUS='OLD')
         READ(7,*) J1
         CLOSE(7)
         IF (J1.EQ.1) THEN
            WRITE(*,'(A)') ' fetchz> Converting initial coordinates from Angstrom to Bohr'
            DO J1=1,NOPT
               Q(J1)=Q(J1)*1.889726164D0
            ENDDO
         ENDIF
         IF (TWOENDS.OR.CONNECTT.OR.NEBT.OR.NEWNEBT.OR.DRAGT.OR.GUESSPATHT
     $     .OR.MECCANOT.OR.MORPHT.OR.GREATCIRCLET.OR.GROWSTRINGT.OR.BHINTERPT.OR.BISECTT.OR.CLOSESTALIGNMENT) THEN
            OPEN (UNIT=7,FILE=FINSTRING,STATUS='OLD') 
            DO J1=1,NATOMS
               READ(7,*) FIN(3*(J1-1)+1),FIN(3*(J1-1)+2),FIN(3*(J1-1)+3)
            ENDDO
            CLOSE(7)
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> Coordinates of second point read from file finish'
         ENDIF
C        PRINT*,'coordinates:'
C        WRITE(*,'(3F20.10)') (Q(J1),J1=1,3*NATOMS)
      ELSE IF (VARIABLES) THEN
         NATOMS=0
100      CALL INPUT(END)
         IF (.NOT.END) THEN
            NATOMS=NATOMS+1
            CALL READF(Q(NATOMS))
            GOTO 100
         ENDIF
         NOPT=NATOMS
         WRITE(*,'(A,I4,A,I4)') ' fetchz> ',NOPT,' variables will be optimised, number of zero eigenvalues=',NZERO
!        WRITE(*,'(A,G20.10)') ' fetchz> coupling parameter=',PARAM1 ! DJW for Michael Kastner potential
         IF (TWOENDS.OR.CONNECTT.OR.NEBT.OR.NEWNEBT.OR.DRAGT.OR.GUESSPATHT
     $     .OR.MECCANOT.OR.MORPHT.OR.GREATCIRCLET.OR.GROWSTRINGT.OR.BHINTERPT.OR.BISECTT.OR.CLOSESTALIGNMENT) THEN
            OPEN (UNIT=7,FILE=FINSTRING,STATUS='OLD') 
            IF (VARIABLES) THEN
               READ(7,*) (FIN(J1),J1=1,NOPT)
               PRINT'(A)',' fetchz> finish variables:'
               WRITE(*,'(6G20.10)') (FIN(J1),J1=1,NOPT)
            ELSE
               DO J1=1,NATOMS
                  READ(7,*) FIN(3*(J1-1)+1),FIN(3*(J1-1)+2),FIN(3*(J1-1)+3)
               ENDDO
            ENDIF
            CLOSE(7)
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> Coordinates of second point read from file finish'
         ENDIF
      ELSE 
C
C  Continue reading coordinates from odata file. This is the default.
C
         LNATOMS=0
300      CALL INPUT(END)
         IF (.NOT.END) THEN
            IF (LNATOMS+1.GT.NATOMS) THEN
               PRINT*,'ERROR check for blank line in odata: LNATOMS,NATOMS=',LNATOMS,NATOMS
               CALL FLUSH(6)
               STOP
            ENDIF
            J=3*LNATOMS
            CALL READA(ZSYM(LNATOMS+1))
            IF ((LNATOMS.EQ.0).AND.(ZSYM(1).EQ.'CK')) ALLOCATE (CHARGE(NATOMS))
            IF (ZSYM(LNATOMS+1)(1:2).EQ.'  ') GOTO 300
            CALL UPPERCASE(ZSYM(LNATOMS+1))
            CALL READF(Q(J+1))
            CALL READF(Q(J+2))
            CALL READF(Q(J+3))
            LNATOMS=LNATOMS+1
            IF (ZSYM(1).EQ.'CK') THEN
               CHARGE(LNATOMS)=Q(J+3)  ! save charge for CK ions
               Q(J+3)=0.0D0
C              PRINT*,'LNATOMS,CHARGE,Q=',LNATOMS,CHARGE(LNATOMS),Q(J+3)
            ENDIF
            GOTO 300
         ENDIF
         CLOSE(5)

!        IF (DEBUG) THEN
!           PRINT '(A)',' fetchz> coordinates:'
!           PRINT '(3F20.10)',Q(1:3*NATOMS)
!        ENDIF
         IF (QSPCFWT) THEN
            WRITE(*,'(3(A,I6))') ' SYSTEM ',NATOMS/3,' QSPCFW flexible water molecules'
         ENDIF
         IF (QTIP4PFT) THEN
            WRITE(*,'(3(A,I6))') ' SYSTEM ',NATOMS/3,' QTIP4PF flexible water molecules'
         ENDIF
         IF (SDT) THEN
            WRITE(*,'(3(A,I6))') ' SYSTEM Stillinger-David potential for ',
     &         SDOXYGEN,' Oxygen and ',SDHYDROGEN,' Hydrogen atoms; charge=',SDCHARGE
            IF (SDHYDROGEN-2*SDOXYGEN.NE.SDCHARGE) THEN
               PRINT '(A)',' fetchz> ERROR in Stillinger-David initialization'
               STOP
            END IF
            CALL SDINIT(SDOXYGEN,SDCHARGE)
C
C  Pathsample will set the atom type to SD. To get the right masses
C  change to O and H now.
C
            DO J1=1,SDOXYGEN
               ZSYM(J1)='O '
            ENDDO
            DO J1=SDOXYGEN+1,SDOXYGEN+SDHYDROGEN
               ZSYM(J1)='H '
            ENDDO
         ENDIF
         IF (BOWMANT) THEN
            WRITE(*,'(A,I1,A,I0,A)') ' SYSTEM Bowman''s PES#',BOWMANPES,' for a cluster of ',NATOMS/3,' water molecules'
            CALL BOWMANINIT(NATOMS/3, BOWMANPES, trim(BOWMANDIR))
            DO J1=1,2*NATOMS/3
               ZSYM(J1)='H '
            ENDDO
            DO J1=2*NATOMS/3+1,NATOMS
               ZSYM(J1)='O '
            ENDDO
         ENDIF
         IF (TTM3T) THEN
C
C  To get the right masses for pathsample change to O and H now.
C  Note that the order is different from Bowman!
C
            WRITE(*,'(A)') ' SYSTEM Xantheas TTM3-F potential for a cluster of water molecules'
            DO J1=1,NATOMS/3
               ZSYM(J1)='O '
            ENDDO
            DO J1=NATOMS/3+1,NATOMS
               ZSYM(J1)='H '
            ENDDO
         ENDIF
C 
C  Dirty trick to try and run Thomson problem as NATOMS/3 *2 effective atoms
C  using angular coordinates instead of xyz.
C
         IF (ZSYM(NATOMS).EQ.'TH') THEN
            IF ((NATOMS/3)*1.0D0.NE.(NATOMS*1.0D0/3)) THEN
               PRINT '(A)',' fetchz> ERROR, the number of atoms needs to be divisible by 3 for Thomson fudge'
               CALL FLUSH(6)
               STOP
            ENDIF
            DO J1=1,NATOMS
               DIST=SQRT(Q(3*(J1-1)+1)**2+Q(3*(J1-1)+2)**2+Q(3*(J1-1)+3)**2)
               Q(3*(J1-1)+1)=Q(3*(J1-1)+1)/DIST
               Q(3*(J1-1)+2)=Q(3*(J1-1)+2)/DIST
               Q(3*(J1-1)+3)=Q(3*(J1-1)+3)/DIST
               THTEMP(2*(J1-1)+1)=ACOS(Q(3*(J1-1)+3))
               IF (ABS(Q(3*(J1-1)+3)-COS(THTEMP(2*(J1-1)+1))).GT.1.0D-10) THEN
                  PRINT '(A)','inconsistent conversion for z'
                  CALL FLUSH(6)
                  STOP
               ENDIF
               IF (Q(3*(J1-1)+1).EQ.0.0D0) THEN
                  IF (Q(3*(J1-1)+2).GT.0.0D0) THEN
                     THTEMP(2*(J1-1)+2)=HALFPI
                  ELSE 
                     THTEMP(2*(J1-1)+2)=-HALFPI
                  ENDIF
               ELSE IF (Q(3*(J1-1)+2).EQ.0.0D0) THEN
                  IF (Q(3*(J1-1)+1).GT.0.0D0) THEN
                     THTEMP(2*(J1-1)+2)=0.0D0
                  ELSE
                     THTEMP(2*(J1-1)+2)=2*HALFPI
                  ENDIF
               ELSE
                  THTEMP(2*(J1-1)+2)=ATAN(Q(3*(J1-1)+2)/Q(3*(J1-1)+1))
               ENDIF
               IF (ABS(Q(3*(J1-1)+1)-SIN(THTEMP(2*(J1-1)+1))*COS(THTEMP(2*(J1-1)+2))).GT.1.0D-5) THEN
                  THTEMP(2*(J1-1)+2)=THTEMP(2*(J1-1)+2)+2*HALFPI
                  IF (ABS(Q(3*(J1-1)+1)-SIN(THTEMP(2*(J1-1)+1))*COS(THTEMP(2*(J1-1)+2))).GT.1.0D-5) THEN
                     PRINT '(A)','inconsistent conversion for x'
                     CALL FLUSH(6)
                     STOP
                  ENDIF
               ENDIF
               IF (ABS(Q(3*(J1-1)+2)-SIN(THTEMP(2*(J1-1)+1))*SIN(THTEMP(2*(J1-1)+2))).GT.1.0D-5) THEN
                  THTEMP(2*(J1-1)+2)=-THTEMP(2*(J1-1)+2)
                  IF (ABS(Q(3*(J1-1)+2)-SIN(THTEMP(2*(J1-1)+1))*SIN(THTEMP(2*(J1-1)+2))).GT.1.0D-5) THEN
                     PRINT '(A)','inconsistent conversion for y'
                     PRINT '(A,3G20.10)','x,y,z:      ',Q(3*(J1-1)+1),Q(3*(J1-1)+2),Q(3*(J1-1)+3)
                     PRINT '(A,3G20.10)','theta,phi: ',THTEMP(2*(J1-1)+1),THTEMP(2*(J1-1)+2)
                     PRINT '(A,3G20.10)','x,y,z calc: ',SIN(THTEMP(2*(J1-1)+1))*COS(THTEMP(2*(J1-1)+2)),
     &                                                  SIN(THTEMP(2*(J1-1)+1))*SIN(THTEMP(2*(J1-1)+2)),
     &                                                  COS(THTEMP(2*(J1-1)+1))
                     CALL FLUSH(6)
                     STOP
                  ENDIF
               ENDIF
            ENDDO
            NATOMS=(NATOMS/3)*2 ! The great big fudge
            LNATOMS=NATOMS
            Q(1:3*NATOMS)=THTEMP(1:3*NATOMS) ! coordinates saved as theta, phi
            PRINT '(A,I6)',' fetchz> Number of atoms changed for the Thomson problem to ',NATOMS
         ENDIF

         IF (MACHINE) THEN
             CALL READINPFILE(Q)
         ENDIF

         IF (ZSYM(NATOMS)(1:1).EQ.'W') THEN
            RIGIDBODY=.TRUE.
C            OPEN(UNIT=18,FILE='Euler',STATUS='OLD')   !   WCOMMENT
C            READ(18,*) (Q(3*LNATOMS+J1),J1=1,3*LNATOMS)
C            CLOSE(18)
         ENDIF

! hk286
!         IF (GTHOMSONT) THEN
!            IF ((NATOMS/3)*1.0D0.NE.(NATOMS*1.0D0/3)) THEN
!               PRINT '(A)',' fetchz> ERROR, the number of atoms needs to be divisible by 3 for GThomson fudge'
!               CALL FLUSH(6)
!               STOP
!            ENDIF
!            PRINT *, 'Odata coordinates, before conversion'
!            DO J1 = 1, 3*NATOMS
!               PRINT *, Q(J1)
!            END DO
!            CALL GTHOMSONCTOANG(Q,THTEMP(1:2*NATOMS),NATOMS,1)
!            NGTHORI = NATOMS
!            NATOMS=(NATOMS/3)*2 ! The great big fudge
!            LNATOMS=NATOMS
!            Q(1:3*NATOMS)=THTEMP(1:3*NATOMS) ! coordinates saved as theta, phi
!            PRINT *, 'Odata coordinates, converted to polars'
!            DO J1 = 1, 3*NATOMS
!              PRINT *, Q(J1)
!            END DO
!            PRINT '(A,I6)',' fetchz> Number of atoms changed for the GThomson problem to ',NATOMS
!         ENDIF

! jwrm2> Generalised Thomson surfaces
         IF (GTHOMSONT) THEN
! Convert Cartesian input coordiantes to polars
            CALL GTHOMSONCTOANG(Q(1:3*NATOMS), THTEMP(1:3*NATOMS), NATOMS, 6)
            NGTHORI = NATOMS
            Q(1:3*NATOMS) = THTEMP(1:3*NATOMS)
!           PRINT *, 'Odata coordinates, converted to polars'
!           DO J1 = 1, 3*NATOMS
!              PRINT *, Q(J1)
!           END DO
         ENDIF
! jwrm2> end

Cdavidg: changes to make lnatoms correct!
         
         IF (USERPOTT) THEN
            LNATOMS=NATOMS
         ENDIF

         LUNIT=GETUNIT()
         IF (TWOENDS.OR.CONNECTT.OR.NEBT.OR.NEWNEBT.OR.DRAGT.OR.GUESSPATHT
     $     .OR.MECCANOT.OR.MORPHT.OR.GREATCIRCLET.OR.GROWSTRINGT.OR.BHINTERPT.OR.BISECTT.OR.CLOSESTALIGNMENT) THEN
            OPEN (UNIT=LUNIT,FILE=FINSTRING,STATUS='OLD') 

            IF (GTHOMSONT) THEN
! hk286
!               DO J1=1,LNATOMS/2*3
!                  READ(LUNIT,*) FIN(3*(J1-1)+1),FIN(3*(J1-1)+2),FIN(3*(J1-1)+3)
!               ENDDO
!               CALL GTHOMSONCTOANG(FIN,THTEMP(1:3*NATOMS),NATOMS/2*3,1)
!               FIN(1:3*NATOMS)=THTEMP(1:3*NATOMS)

! jwrm2> using TWOD framework for Gthomson coordinates
               DO J1 = 1, NATOMS
                  READ(LUNIT,*) FIN(3*(J1-1)+1), FIN(3*(J1-1)+2), FIN(3*(J1-1)+3)
               END DO
               CALL GTHOMSONCTOANG(FIN, THTEMP(1:3*NATOMS), NATOMS, 1)
               FIN(1:3*NATOMS) = THTEMP(1:3*NATOMS)
! jwrm2> end
            ELSE
               DO J1=1,LNATOMS
                  READ(LUNIT,*) FIN(3*(J1-1)+1),FIN(3*(J1-1)+2),FIN(3*(J1-1)+3)
               ENDDO
            ENDIF
            CLOSE(LUNIT)
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> Coordinates of second point read from file finish'
C           IF (ZSYM(1)(1:1).EQ.'W') THEN   !  WCOMMENT
C              OPEN (UNIT=7,FILE='Euler.finish',STATUS='OLD')
C              DO J1=1,LNATOMS
C                 READ(7,*) FIN(3*LNATOMS+3*(J1-1)+1),FIN(3*LNATOMS+3*(J1-1)+2),FIN(3*LNATOMS+3*(J1-1)+3)
C              ENDDO
C              CLOSE(7)
C              WRITE(*,'(A,I4,A,I4,A)') ' fetchz> Euler angles of second point read from file finish'
C           ENDIF
         ENDIF

         IF (LNATOMS.NE.NATOMS) THEN
            PRINT*,'ERROR check for blank line in odata: LNATOMS,NATOMS=',LNATOMS,NATOMS
         ENDIF
         NOPT=3*NATOMS
         IF (NTAG.GT.0) THEN
            PRINT '(A,I6,A)',' fetchz> ',NTAG,' tagged atoms and mass factors:'
            PRINT '(I6,F12.2)',(TAGNUM(J1),TAGFAC(TAGNUM(J1)),J1=1,NTAG)
         ENDIF
         IF (ZSYM(NATOMS).EQ.'IN') THEN
            WRITE(*,'(A,I4,A,I4,A)') 
     1        ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' shielded Born-Meyer ions'
            WRITE(*,'(A,F12.8,A,I2,A,F12.8,A,F12.8,A,F12.8,A,F12.8)') 
     1         ' fetchz> gamma=',PARAM1,' charge=',INT(PARAM2),' rho=',PARAM3,' A++=',PARAM4,' A--=',PARAM5,' A+-=',PARAM6
         ELSE IF (ZSYM(1).EQ.'CA') THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Ca+Ar atoms'
         ELSE IF (ZSYM(1).EQ.'LM') THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,
     1              ' Cartesian coordinates will be optimised for one ion + ',NATOMS,' inert gas atoms'
            WRITE(*,'(A,3G12.4)') ' fetchz> Mason-Schamp potential with epsilon, rm, gamma=',PARAM1,PARAM2,PARAM3
         ELSE IF (ZSYM(NATOMS).EQ.'AR') THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Lennard-Jones atoms'
            WRITE(*,'(A)') ' fetchz> Coordinates will be scaled (divided) by 3.4 before after input'
         ELSE IF (ZSYM(NATOMS).EQ.'RM') THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Lennard-Jones atoms'
         ELSE IF (WELCH) THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Welch ions'
            WRITE(*,'(A,F12.8,A,F12.8,A,F12.8,A,F12.8)') ' fetchz> A++=',APP,' A--=',AMM,' A+-=',APM,' rho=',RHO
            WRITE(*,'(A,F12.8,A,F12.8,A,F12.8,A,F12.8)') ' fetchz> Q+=',XQP,' Q-=',XQM,' alpha+=',ALPHAP,' alpha-=',ALPHAM
         ELSE IF (ZSYM(NATOMS).EQ.'SY') THEN
            WRITE(*,'(A,I4,A,F15.5,A,F15.5)') ' fetchz> ',(NATOMS/2),' Stockmayer atoms with mu=',
     &                                            STOCKMU,' and lambda=',STOCKLAMBDA
         ELSE IF (ZSYM(NATOMS).EQ.'CD') THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' rigid body coordinates will be optimised for ',NATOMS/2,' capsid pentamers'
            WRITE(*,'(4(A,F12.8))') ' fetchz> Morse rho=',CAPSRHO,' repulsive site epsilon=',
     1                                     EPS2,' pentamer radius=',RAD,' height=',HEIGHT
         ELSE IF (TOSI) THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Tosi-Fumi ions'
            WRITE(*,'(4(A,F12.8))') ' fetchz> A++=',PARAM1,' A--=',PARAM2,' A+-=',PARAM3,' rho=',PARAM4
            IF (TOSIC6) THEN
               WRITE(*,'(A)') ' fetchz> C6 coefficients: '
               WRITE(*,'(A,F12.8,A,F12.8,A,F12.8)') ' fetchz> C6++=',C6PP,' C6--=',C6MM,' C6+-=',C6PM
            ENDIF
            IF (TOSIPOL) THEN
               WRITE(*,'(A)') ' fetchz> Polarizabilities:'
               WRITE(*,'(A,F12.8,A,F12.8)') ' fetchz> alpha+=',ALPHAP,' alpha-=',ALPHAM
               WRITE(*,'(A,F12.8,A)') ' fetchz> damping coefficent=',DAMP,' per bohr'
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'AZ') THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Aziz Ar atoms'
            IF (PARAM1.NE.0.0D0) THEN
               WRITE(*,'(A,F12.8)') ' fetchz> Axilrod-Teller potential will be added with Z*=',PARAM1
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'DZ') THEN
            WRITE(*,'(A,I4,A,I4,A)') 
     1        ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Dzugutov atoms'
            WRITE(*,'(A,F12.8,A,F12.8,A,F12.8,A,F12.8,A,F12.8,A,F12.8,A,F12.8)') 
     1         ' fetchz> m=',PARAM1,' A=',PARAM2,' c=',PARAM3,' aa=',PARAM4,' B=',PARAM5,' d=',PARAM6,' bb=',PARAM7
         ELSE IF (ZSYM(NATOMS).EQ.'AX') THEN
            WRITE(*,'(A,I4,A,I4,A)') 
     1        ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Lennard-Jones+Axilrod-Teller atoms'
            WRITE(*,'(A,F12.8)') ' fetchz> Axilrod-Teller Z*=',PARAM1
         ELSE IF (NATBT) THEN
            WRITE(*,'(A,I4,A,I4,A)') 
     1        ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' TB Na atoms'
         ELSE IF (ZSYM(NATOMS).EQ.'SW') THEN
            WRITE(*,'(A,I4,A,I4,A)') 
     1        ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Stillinger-Weber Si atoms'
            WRITE(*,'(A,F15.8,A,F15.8,A,F15.8)') ' fetchz> Box lengths: x ',PARAM1,', y ',PARAM2,', z ',PARAM3
         ELSE IF (ZSYM(NATOMS).EQ.'QI') THEN
            WRITE(*,'(A,I4,A,I4,A)')
     1        ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' atoms'
            WRITE(*,'(A,F15.8,A,F15.8,A,F15.8)') ' fetchz> Box lengths: x ',PARAM1,', y ',PARAM2,', z ',PARAM3
         ELSE IF (ZSYM(NATOMS).EQ.'SM') THEN
            WRITE(*,'(A,I4,A,I4,A)') 
     1        ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Stillinger-Weber Si atoms'
            WRITE(*,'(A)') ' with the three-body term multiplied by 1.5'
            WRITE(*,'(A,F15.8,A,F15.8,A,F15.8)') ' fetchz> Box lengths: x ',PARAM1,', y ',PARAM2,', z ',PARAM3
         ELSE IF (ZSYM(NATOMS).EQ.'SI') THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Murrell Si potential'
         ELSE IF (ZSYM(NATOMS).EQ.'JC') THEN
            WRITE(*,'(A,I4,A,I4,A)') 
     1        ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' general Murrell potential'
         ELSE IF (ZSYM(NATOMS).EQ.'CC') THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Murrell C potential'
         ELSE IF (ZSYM(NATOMS).EQ.'PH4') THEN
            WRITE(*,'(A,I4,A,I4,A)')
     1       ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Phi4 model potential'
         ELSE IF (ZSYM(NATOMS).EQ.'JM') THEN
            WRITE(*,'(A,I4,A,I4,A)') 
     1        ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' general Murrell potential'
            WRITE(*,'(A,F15.8,A,F15.8,A,F15.8,A,F15.8)') 
     1        ' fetchz> Box lengths: x ',PARAM1,', y ',PARAM2,', z ',PARAM3,', cutoff ',PARAM4
         ELSE IF (ZSYM(NATOMS).EQ.'M') THEN
            WRITE(*,'(A,I4,A,I4,A,F12.4)') 
     1        ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Morse atoms, rho=',PARAM1
         ELSE IF (ZSYM(NATOMS).EQ.'TT') THEN
            WRITE(*,'(A,I4,A,I4,A)') 
     1       ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Menon and Subaswamy April 1997 PRB Si atoms'
            WRITE(*,'(A,F15.8,A,F15.8,A,F15.8,A,F15.8)') 
     1       ' fetchz> Box lengths: x ',PARAM1,', y ',PARAM2,', z ',PARAM3,', cutoff (fractional) ',PARAM4
         ELSE IF (ZSYM(NATOMS).EQ.'LP') THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Lennard-Jones atoms'
            WRITE(*,'(A,F15.8,A,F15.8,A,F15.8,A,F15.8)') 
     1       ' fetchz> Box lengths: x ',PARAM1,', y ',PARAM2,', z ',PARAM3,', cutoff (fractional) ',PARAM4
            IF (BINARY) WRITE(*,'(A,I4,A,I4,A,F11.5,A,F11.5,A,F11.5,A,F11.5)') 
     1         ' fetchz> Binary mixture: ',NTYPEA,' A atoms,',NATOMS-NTYPEA,' B atoms, eps(AB)=',EPSAB,' eps(BB)=',EPSBB,
     2         ' sigma(AB)=',SIGAB,' sigma(BB)=',SIGBB
         ELSE IF ((ZSYM(NATOMS).EQ.'LS').OR.(ZSYM(NATOMS).EQ.'BC')) THEN
!
! BINARY could be true or false for atom type 'LS'. For atom type 'BC'
! it should have been set true using the corresponding keyword in odata
!
            IF ((ZSYM(NATOMS).EQ.'BC').AND.(.NOT.BINARY)) THEN
               PRINT '(A)',' fetchz> ERROR - atom type BC requires the BINARY keyword in odata'
               CALL FLUSH(6)
               STOP
            ENDIF
            WRITE(*,'(A,I4,A,I4,A)') 
     1        ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' shifted, truncated Lennard-Jones atoms'
            WRITE(*,'(A,F15.8,A,F15.8,A,F15.8,A,F15.8)') 
     1        ' fetchz> Box lengths: x ',PARAM1,', y ',PARAM2,', z ',PARAM3,', cutoff (fractional) ',PARAM4
            IF (BINARY) WRITE(*,'(A,I4,A,I4,A,F11.5,A,F11.5,A,F11.5,A,F11.5)') 
     1         ' fetchz> Binary mixture: ',NTYPEA,' A atoms,',NATOMS-NTYPEA,' B atoms, eps(AB)=',EPSAB,' eps(BB)=',EPSBB,
     2         ' sigma(AB)=',SIGAB,' sigma(AB)=',SIGBB
         ELSE IF (ZSYM(NATOMS).EQ.'LK') THEN
            WRITE(*,'(A,I4,A,I4,A)') 
     1        ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Kob/Sciortino Lennard-Jones atoms'
            WRITE(*,'(A,F15.8,A,F15.8,A,F15.8,A,F15.8)') 
     1        ' fetchz> Box lengths: x ',PARAM1,', y ',PARAM2,', z ',PARAM3,', cutoff (fractional) ',PARAM4
            IF (BINARY) WRITE(*,'(A,I4,A,I4,A,F11.5,A,F11.5,A,F11.5,A,F11.5)') 
     1         ' fetchz> Binary mixture: ',NTYPEA,' A atoms,',NATOMS-NTYPEA,' B atoms, eps(AB)=',EPSAB,' eps(BB)=',EPSBB,
     2         ' sigma(AB)=',SIGAB,' sigma(AB)=',SIGBB
         ELSE IF (ZSYM(NATOMS).EQ.'LC') THEN
            WRITE(*,'(A,I4,A,I4,A)') 
     1        ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' shifted, truncated Lennard-Jones atoms'
            WRITE(*,'(A,F15.8,A,F15.8,A,F15.8,A,F15.8)') 
     1        ' fetchz> Box lengths: x ',PARAM1,', y ',PARAM2,', z ',PARAM3,', cutoff (fractional) ',PARAM4
            IF (BINARY) WRITE(*,'(A,I4,A,I4,A,F11.5,A,F11.5,A,F11.5,A,F11.5)') 
     1         ' fetchz> Binary mixture: ',NTYPEA,' A atoms,',NATOMS-NTYPEA,' B atoms, eps(AB)=',EPSAB,' eps(BB)=',EPSBB,
     2         ' sigma(AB)=',SIGAB,' sigma(AB)=',SIGBB
         ELSE IF (ZSYM(NATOMS).EQ.'MP') THEN
            WRITE(*,'(A,I4,A,I4,A,F12.4)') 
     1       ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Morse atoms, rho=',PARAM1
            WRITE(*,'(A,F15.8,A,F15.8,A,F15.8,A,F15.8)') 
     1       ' fetchz> Box lengths: x ',PARAM2,', y ',PARAM3,', z ',PARAM4,', cutoff (fractional) ',PARAM5
         ELSE IF (ZSYM(NATOMS).EQ.'DS') THEN
            WRITE(*,'(A,I4,A,I4,A,F12.4)') 
     1       ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Morse atoms, rho=',PARAM1
            WRITE(*,'(A,F15.8,A,F15.8,A,F15.8,A,F15.8)') 
     1       ' fetchz> Box lengths: x ',PARAM2,', y ',PARAM3,', z ',PARAM4,', cutoff (fractional) ',PARAM5
         ELSE IF (ZSYM(NATOMS).EQ.'GP') THEN
            WRITE(*,'(A,I6,A,I4,A,I4)') 
     1       ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Gupta atoms for atom type',GUPTATYPE
         ELSE IF (ZSYM(NATOMS).EQ.'MS') THEN
            WRITE(*,'(A,I4,A,I4,A,F12.4)') 
     1       ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Morse atoms, rho=',PARAM1
            WRITE(*,'(A,F15.8,A,F15.8,A,F15.8,A,F15.8)') 
     1       ' fetchz> Box lengths: x ',PARAM2,', y ',PARAM3,', cutoff (fractional) ',PARAM4
         ELSE IF (ZSYM(NATOMS).EQ.'CK') THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',2*NATOMS,' Cartesian coordinates will be optimised for ',NATOMS,' 2D trapped ions'
         ELSE IF (EYTRAPT) THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ', 
     &                               NATOMS,' trapped ions'
            WRITE(*,'(A,G12.4,A,I4)') ' fetchz> with radial potential 0.5 * ',TRAPK,' * r**',NTRAPPOW
         ELSE IF (MKTRAPT) THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ', 
     &                               NATOMS,' trapped ions (MK potential)'
         ELSE IF (ZSYM(NATOMS).EQ.'BE') THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' trapped ions'
         ELSE IF ((.NOT. NIMET) .AND.(.NOT. NIHEAM7T) .AND.
     &       ((ZSYM(NATOMS).EQ.'AU').OR.(ZSYM(NATOMS).EQ.'PT').OR.(ZSYM(NATOMS).EQ.'AG').OR.(ZSYM(NATOMS).EQ.'NI'))) THEN
            IF ((ZSYM(NATOMS).EQ.'AU').OR.(ZSYM(NATOMS).EQ.'PT')) THEN
               WRITE(*,'(A,I4,A,I4,A)') 
     1          ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Sutton-Chen Au-type atoms'
               WRITE(*,'(A,F15.6,A,F15.6,A,F15.6)') ' fetchz> n=10, m=8, eps=',PARAM1,' c=',PARAM2,' sigma=',PARAM3
            ELSE IF (ZSYM(NATOMS).EQ.'AG') THEN
               WRITE(*,'(A,I4,A,I4,A)') 
     1          ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Sutton-Chen Ag-type atoms'
               WRITE(*,'(A,F15.6,A,F15.6,A,F15.6)') ' fetchz> n=12, m=6, eps=',PARAM1,' c=',PARAM2,' sigma=',PARAM3
            ELSE IF (ZSYM(NATOMS).EQ.'NI') THEN
               WRITE(*,'(A,I4,A,I4,A)') 
     1          ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Sutton-Chen Ni-type atoms'
               WRITE(*,'(A,F15.6,A,F15.6,A,F15.6)') ' fetchz> n=9, m=6, eps=',PARAM1,' c=',PARAM2,' sigma=',PARAM3
            ENDIF
         ELSE IF (ZSYM(NATOMS).EQ.'PL') THEN
            WRITE(*,'(A,I4,A,I4,A)') 
     1       ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for the 3-colour 46-bead model polypeptide'
         ELSE IF (BLNT) THEN
            WRITE(*,'(A,I4,A,I4,A)')
     1       ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for the general 3-colour bead model'
         ELSE IF (ZSYM(NATOMS).EQ.'GL') THEN
            WRITE(*,'(A,I4,A,I4,A)') 
     1         ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for the Go-model 3-colour 46-bead model polypeptide'
         ELSE IF (ZSYM(NATOMS).EQ.'SC') THEN
            INQUIRE(FILE='SCparams',EXIST=YESNO)
            IF (.NOT.YESNO) THEN
               NN=12
               MM=6
               EPS=2.5415D-03
               CSC=144.41
               SIG=1.414
            ELSE
               OPEN(UNIT=33,FILE='SCparams',STATUS='OLD')
               READ(33,*) NN,MM,EPS,CSC,SIG
               CLOSE(33)
            ENDIF
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Sutton-Chen atoms'
C           WRITE(*,'(A,F15.6,A,F15.6,A,F15.6)') ' fetchz> n=',NN,' m=',MM,' eps=',PARAM1,' c=',PARAM2,' sigma=',PARAM3
            WRITE(*,'(A,F15.8,A,F15.8,A,F15.8,A,F15.8)') 
     1       ' fetchz> Box lengths: x ',PARAM1,', y ',PARAM2,', z ',PARAM3,', cutoff ',PARAM4
         ELSE IF (ZSYM(NATOMS).EQ.'PR') THEN
            WRITE(*,'(A,I4,A,I4,A)') 
     1       ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Pacheco-Ramalho C60 molecules'
         ELSE IF (ZSYM(NATOMS).EQ.'C6') THEN
            WRITE(*,'(A,I4,A,I4,A)') 
     1       ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Girifalco C60 molecules'
         ELSE IF ( ZSYM(NATOMS).EQ.'SB') THEN
            WRITE(*,'(A,I6,A,I5,A)')
     1       ' SYSTEM ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' SBM sites'
         ELSE IF (ZSYM(NATOMS).EQ.'P6') THEN
            WRITE(*,'(A,I4,A,I4,A)') 
     1       ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Girifalco C60 molecules'
            WRITE(*,'(A,F15.8,A,F15.8,A,F15.8,A,F15.8)') 
     1       ' fetchz> Box lengths: x ',PARAM1,', y ',PARAM2,', z ',PARAM3,', cutoff ',PARAM4
         ELSE IF (ZSYM(NATOMS).EQ.'FH') THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Fenske-Hall atoms'
         ELSE IF (ZSYM(NATOMS)(1:1).EQ.'W') THEN
            IF (ZSYM(NATOMS).EQ.'W4') WRITE(*,'(A,I4,A,I4,A)') 
     2       ' fetchz> ',3*NATOMS,' coordinates will be optimised for ',NATOMS/2,' TIP4P water molecules'
            IF (ZSYM(NATOMS).EQ.'W3') WRITE(*,'(A,I4,A,I4,A)') 
     2       ' fetchz> ',3*NATOMS,' coordinates will be optimised for ',NATOMS/2,' TIP3P water molecules'
            IF (ZSYM(NATOMS).EQ.'W2') WRITE(*,'(A,I4,A,I4,A)') 
     2       ' fetchz> ',3*NATOMS,' coordinates will be optimised for ',NATOMS/2,' TIPS2 water molecules'
            IF (ZSYM(NATOMS).EQ.'W1') WRITE(*,'(A,I4,A,I4,A)') 
     2       ' fetchz> ',3*NATOMS,' coordinates will be optimised for ',NATOMS/2,' TIPS1 water molecules'
         ELSE IF (ZSYM(NATOMS).EQ.'ME') THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Mie atoms'
            WRITE(*,'(A,I3,A,I3)') ' fetchz> n=',INT(PARAM1),' m=',INT(PARAM2)
            WRITE(*,'(A,F15.8,A,F15.8,A,F15.8,A,F15.8)') 
     1       '  SYSTEM Box lengths: x ',PARAM1,', y ',PARAM2,', z ',PARAM3,', cutoff ',PARAM4
         ELSE IF (CADPAC) THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' CADPAC atoms'
            WRITE(*,'(A,A,A,A)') ' fetchz> System name: ',SYS(1:LSYS),', edit command: ',EDITIT
         ELSE IF (GAMESSUS) THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' GAMESS-US atoms'
            WRITE(*,'(A,A,A,A)') ' fetchz> System name: ',SYS(1:LSYS),', edit command: ',EDITIT
         ELSE IF (GAMESSUK) THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' GAMESS-UK atoms'
            WRITE(*,'(A,A,A,A)') ' fetchz> System name: ',SYS(1:LSYS),', edit command: ',EDITIT
         ELSE IF (GAUSSIAN) THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Gaussian atoms'
         ELSE IF (GAUSSIAN03) THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Gaussian03 atoms'
         ELSE IF (GAUSSIAN09) THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Gaussian09 atoms'
         ELSE IF (GAUSSIAN16) THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Gaussian16 atoms'
         ELSE IF (DFTBT) THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' Tiffany TB atoms'
         ELSE IF (ZSYM(NATOMS).EQ.'SV') THEN
            WRITE(*,'(A,I4,A,I4,A)') ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' MSEVB atoms '

! Set some values

            num_eig = (NATOMS-1)/3
            num_hyd = NATOMS - num_eig

! Initialise some memory

            maxNumVBstates = 2

            do j1 = 1, shellsToCount
               maxNumVBstates = maxNumVBstates + 4*(2**(j1-1))    ! 4 to account for multiple pivot states
            enddo

! Hard limit max number of VB shells at the moment

            if (maxNumVBstates.gt.1000) maxNumVBstates = 1000
         
            ALLOCATE(each_coulomb(natoms,natoms), water_inter_coulomb(natoms,natoms), ljr(natoms,natoms))
            ALLOCATE(inter_coulomb(natoms,natoms), lj_inter(natoms,natoms), repulse_inter(natoms,natoms))
            ALLOCATE(atom_coulomb(natoms))

            ALLOCATE(atmpl(maxNumVBstates,natoms), statesInteract(maxNumVBstates,maxNumVBstates))
            ALLOCATE(zundel_species(maxNumVBstates,maxNumVBstates,7))
            ALLOCATE(zundel_f(maxNumVBstates,maxNumVBstates),zundel_g(maxNumVBstates,maxNumVBstates))

            ALLOCATE(psix(natoms), psiy(natoms), psiz(natoms))
            ALLOCATE(interAtomicR(natoms,natoms))            
         ENDIF
         IF (CPMD) WRITE(*,'(A,I6,A,I6,A)') 
     1          ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,' CPMD atoms; bulk boundary conditions'
         IF (CPMDC) WRITE(*,'(A,I6,A,I6,A)')
     1          ' fetchz> ',NOPT,' Cartesian coordinates will be optimised for ',NATOMS,
     1              ' CPMD atoms; cluster boundary conditions'
         
         IF (TWOD) WRITE(*,'(A)') ' fetchz> Two-dimensional flatland enforced'
         IF (DOUBLET) WRITE(*,'(A,F12.5,A,F12.5,A,F12.5)') 
     1      ' fetchz> Double-well potential between first two atoms, barrier=',PARAM5,
     1                               ' first minimum at ',PARAM4,' second at ',PARAM4+2.0D0*PARAM6
      ENDIF
      IF (ZSYM(1)(1:1).EQ.'W') THEN ! must generalise for other rigid bodies somehow
         DEALLOCATE (IATNUM,ATMASS)
         ALLOCATE (IATNUM(3*(NATOMS/2)),ATMASS(3*(NATOMS/2)))
      ENDIF
C
C DAE ATMASS is set in CHSETZSYMATMASS for CHARMM. PERTABLE would overwrite
C jmc similarly for unres in UNRSETZSYMATMASS...
C
      IF (.NOT.(CHRMMT.OR.UNRST.OR.AMBERT.OR.NABT.OR.AMBER12T.OR.RINGPOLYMERT.OR.OPEPT)) CALL PERTABLE
      IF (CASTEP) THEN
         INQUIRE(FILE='castep.masses',EXIST=YESNO)
         IF (YESNO) THEN
            LUNIT=GETUNIT()
            OPEN(UNIT=LUNIT,FILE='castep.masses',STATUS='OLD')
            READ(LUNIT,*) (ATMASS(J1),J1=1,NATOMS)
            CLOSE(LUNIT)
            PRINT '(A)',' fetchz> Masses replaced by values from castep.masses file, Values are:'
            WRITE(*,'(F20.10)') (ATMASS(J1),J1=1,NATOMS)
         ENDIF
         DO J1=1,NATOMS
            DO J2=J1+1,NATOMS
               IF (ATMASS(J2).LT.ATMASS(J1)) THEN
                  PRINT '(A)',' fetchz> CASTEP input error: atoms must be ordered by atomic number'
                  PRINT '(A,I6,A,G20.10)',' mass for atom ',J2,' is ',ATMASS(J2)
                  PRINT '(A,I6,A,G20.10)',' mass for atom ',J1,' is ',ATMASS(J1)
                  STOP
               ENDIF
            ENDDO
         ENDDO
      ENDIF
      DO J1=1,NOPT
         STPMAX(J1)=MXSTP
C        IF (ZSYM(NATOMS)(1:1).EQ.'W') STPMAX(J1+3*NATOMS)=MXSTP ! WCOMMENT
      ENDDO
C
C  More printing
C
      PRINT*
      IF (READMASST) THEN
         INQUIRE(FILE='masses',EXIST=YESNO)
         IF (YESNO) THEN
            LUNIT=GETUNIT()
            OPEN(UNIT=LUNIT,FILE='masses',STATUS='OLD')
            READ(LUNIT,*) (ATMASS(J1),J1=1,NATOMS)
            CLOSE(LUNIT)
            PRINT '(A)',' fetchz> Masses replaced by values from masses file, Values are:'
            WRITE(*,'(F20.10)') (ATMASS(J1),J1=1,NATOMS)
         ELSE
            PRINT '(A)',' fetchz> ERROR *** keyword READMASS found, but no masses file'
            STOP
         ENDIF
      ENDIF
      IF (TIMELIMIT.LT.HUGE(TIMELIMIT)/2)  WRITE(*,'(A,F20.1)') ' fetchz> Time limit (s): ',TIMELIMIT
      IF (PV.OR.PVTS) THEN
         WRITE(*,'(A,F20.10)') ' fetchz> Box lengths will be varied to give a constant pressure of ',PRESS
         WRITE(*,'(A,F20.10)') ' fetchz> Box length gradient convergence criterion ',PVCONV
         WRITE(*,'(A,G20.10)') ' fetchz> Box length gradient tolerance for multiple steps ',PVTOL
         WRITE(*,'(A,G20.10)') ' fetchz> Constant pressure calculation will use using fractional coordinates'
         IF (CUBIC) WRITE(*,'(A)') ' fetchz> A cubic cell will be maintained if initial cell is cubic'
         IF (PVTS) THEN
            IF (NBOXTS.EQ.1) WRITE(*,'(A)') ' fetchz> Searching for a transition state in x box length'
            IF (NBOXTS.EQ.2) WRITE(*,'(A)') ' fetchz> Searching for a transition state in y box length'
            IF (NBOXTS.EQ.3) WRITE(*,'(A)') ' fetchz> Searching for a transition state in z box length'
         ENDIF
      ENDIF
      IF (FRACTIONAL) WRITE(*,'(A)') ' fetchz> Fractional coordinates to be used in box length optimisation at constant pressure'
      IF (REPELTST) THEN
         NREPELTS=0
         OPEN(UNIT=77,FILE='points.repel',STATUS='OLD',ERR=10)
         DO J1=1,100
            READ(77,*,END=10) (REPELTS(J2,J1),J2=1,3*NATOMS)
            NREPELTS=NREPELTS+1
         ENDDO
10       CLOSE(77)
         WRITE(*,'(A,I5,A,F15.5)') 
     1          ' fetchz> Search will be displaced from ',NREPELTS,' points found in points.repel, step=',REPELPUSH
      ENDIF
      IF (.NOT.(BFGSMINT.OR.BSMIN.OR.RKMIN)) THEN
         IF (ISTCRT.EQ.0) WRITE(*,'(A)') ' fetchz> Scaling steps according to the total Cartesian displacement'
         IF (ISTCRT.EQ.1) WRITE(*,'(A)') ' fetchz> Scaling steps according to the largest atomic Cartesian displacement'
         IF (ISTCRT.EQ.10) WRITE(*,'(A,F12.5)') ' fetchz> Scaling steps according to trust radius=',TRAD
      ENDIF
      IF (FIXD) THEN
         IF (T12FAC.LE.1.0D0) THEN
            WRITE(*,'(A,F12.4)') ' fetchz> Fixed initial uphill search direction, fraction of first collision time used=',T12FAC
         ELSE
            WRITE(*,'(A,F12.4)') ' fetchz> Fixed initial uphill search direction, first colliding atoms will be placed half way'
         ENDIF
      ENDIF
      IF (READV) WRITE(*,'(A)') ' fetchz> Initial eigenvector will be read from vector.dump'
      IF (READSP) WRITE(*,'(A)') ' fetchz> stationary point info will be read'
      IF (DUMPSP) WRITE(*,'(A)') ' fetchz> stationary point info will be dumped'
      IF (FREEZE) THEN
         WRITE(*,'(A,I6,A)') ' fetchz> ', NFREEZE,' atoms or variables will be frozen:'
         DO J1=1,NATOMS
            IF (FROZEN(J1)) WRITE(*,'(I6)') J1
         ENDDO
      ENDIF
      IF (MLPVB3T) THEN
         IF (MLPVB3NNT) THEN
            DO J1=1,NATOMS
               IF (FROZEN(J1)) Q(J1)=0.0D0 ! for frozen zero nearest-neighbour setup
            ENDDO
         ENDIF
         WRITE(*,'(I4,A)') NATOMS,' link weights for MLPVB3'
      ENDIF

      ! jbr36 classical rate calculations
      IF (CLASSICALRATEST) THEN
         WRITE(*,'(A,F12.5,A,E20.10)') ' fetch> Classical rates and corrections will be calculated:'
      ENDIF

      IF (NORESET) WRITE(*,'(A)') ' fetchz> Atoms will not be returned to the primary supercell'

      IF (RTEST) THEN
         IF (JZ.NE.0.0D0) THEN 
            WRITE(*,'(A,F12.4)') ' fetchz> Adding centrifugal potential for angular momentum Jz=',JZ
         ELSE
            WRITE(*,'(A,F12.4)') ' fetchz> Adding centrifugal potential for angular velocity omega=',OMEGA
         ENDIF
      ENDIF

      IF ((INR.EQ.2).AND.(HINDEX.NE.1)) THEN
         WRITE(*,'(A,I4)') ' fetchz> Searching for a saddle with Hessian index=',HINDEX
      ENDIF
      IF ((INR.EQ.2).AND.(KEEPINDEX)) THEN
         WRITE(*,'(A)') ' fetchz> Searching for a saddle with Hessian index determined by starting point'
      ENDIF
   
      IF (NOIT) THEN
         WRITE(*,'(A,I3,A)') ' fetchz> Lowest ',HINDEX,' eigenvalues and eigenvectors will be calculated non-iteratively'
      ENDIF
      IF (NRANROT.GT.0) THEN
         WRITE(*,'(A,I6)') ' fetchz> Number of random initial orientations in permutational alignment=',NRANROT
      ENDIF
      IF (ALIGNRBST) THEN
         WRITE(*,'(A,I3,A)') 'fetchz> Structural alignments will be performed on a subset of ', N_TO_ALIGN, 'rigid bodies'
      ENDIF
      IF (FIELDT) THEN
         WRITE(*,'(A,F12.4)') 'SETTINGS Adding centrifugal potential for Jz=',JZ
         IF (D5HT) WRITE(*,'(A,F12.4)') ' fetchz> Adding D5h symmetry field strength=',FD5H
         IF (OHT) WRITE(*,'(A,F12.4)') ' fetchz> Adding D5h symmetry field strength=',FOH
         IF (IHT) WRITE(*,'(A,F12.4)') ' fetchz> Adding D5h symmetry field strength=',FIH
         IF (TDT) WRITE(*,'(A,F12.4)') ' fetchz> Adding D5h symmetry field strength=',FTD
      ENDIF

      IF (ADMT) WRITE(*,'(A,I4,A)') ' fetchz> Distance matrix will be printed every ',NADM,' cycles'
      IF (.NOT.BULKT) THEN
         WRITE(*,'(A,F15.8,A,I3)') 
     1  ' fetchz> Point group checked when RMS force <',SYMCUT,', highest symmetry axis tested for=',NHCHECK
         WRITE(*,'(A,2F15.8)') ' fetchz> Initial distance and eigenvalue tolerances in symmetry determination=',TOLD,TOLE
      ENDIF
      WRITE(*,'(A,I6)') ' fetchz> Minimum number of optimization steps=',NSTEPMIN
      IF (NEBT) then
          WRITE(*,'(A,I5,A,I5,A,F12.5)') ' fetchz> NEB parameters: ',NIMAGE,
     1          ' images, a maximum of ',NSTEPNEB,' steps, and RMS convergence criterion ',RMSNEB
      ENDIF
      IF (NEBRESEEDT)  THEN
         WRITE(*,'(A,I5,A)') ' fetchz> DNEB images will be reseeded every ',NEBRESEEDINT,' steps'
         WRITE(*,'(A,G20.10)') '         for energies exceeding ',NEBRESEEDEMAX
         WRITE(*,'(A,G20.10)') '         or energy above highest end point exceeds ',NEBRESEEDBMAX
         WRITE(*,'(A,G20.10,I6)') '         factor and power parameters for repulsive sites A: ',NEBRESEEDDEL1,NEBRESEEDPOW1
         WRITE(*,'(A,G20.10,I6)') '         factor and power parameters for repulsive sites B: ',NEBRESEEDDEL2,NEBRESEEDPOW2
      ENDIF
      IF (INTLJT) THEN
         PRINT '(A,F15.5)',   ' fetchz> Using interpLJ potential for initial interpolation in each cycle'
         PRINT '(A,I8)',      '         maximum optimization steps for interpLJ potential=',INTLJSTEPS
         PRINT '(A,F15.5)',   '         image distance spring constant=',KINT
         PRINT '(A,I8)',      '         number of initial intermediate images for interpLJ potential=',INTIMAGE
         PRINT '(A,F15.5)',   '         RMS gradient per image tolerance for constrained potential=',INTLJTOL
         PRINT '(A,F20.10)',  '         minimum distance difference for internal minimum=',INTLJDEL
         PRINT '(A,F20.10)',  '         multiplying factor for internal minimum penalty function=',INTLJEPS
      ENDIF
      IF (INTCONSTRAINTT) THEN
         PRINT '(A,F15.5)',   ' fetchz> Using QCI potential for initial interpolation in each cycle'
         PRINT '(A,F15.5)',   '         with absolute distance change tolerance ',INTCONSTRAINTTOL
         IF (CONCUTABST) THEN
            PRINT '(A,F15.5)',   '         extra distance before turning on constraint potential ',CONCUTABS
         ELSEIF (CONCUTFRACT) THEN
            PRINT '(A,F15.5)',   '         extra fractional distance before turning on constraint potential ',CONCUTFRAC
         ENDIF
         PRINT '(A,F15.5)',   '         constraint spring constant=',INTCONSTRAINTDEL
         PRINT '(A,F15.5)',   '         image distance spring constant=',KINT
         PRINT '(2(A,F15.5))','         repulsion factor between unconstrained atoms=',INTCONSTRAINTREP 
         PRINT '(A,F15.5,A)', '         repulsion cutoff will be the minimum of ',INTCONSTRAINREPCUT,
     &                                ' and the shortest distance in the end points'
         PRINT '(A,F15.5)',   '         fraction for restoring true potential=',INTCONFRAC
         PRINT '(A,I6)',      '         maximum separation of atoms in sequence for constraint=',INTCONSEP
         PRINT '(A,I6)',      '         minimum separation of atoms in sequence for repulsion=',INTREPSEP
         PRINT '(A,I8)',      '         maximum optimization steps for constrained potential=',INTSTEPS1
         PRINT '(A,3I8)',      '         initial # images for constrained potential, maximum, check interval=',  
     &                         INTIMAGE,MAXINTIMAGE,INTIMAGECHECK
         PRINT '(A,2I8)',     '         number of interpolation attempts and initial image increment=',INTNTRIESMAX,INTIMAGEINCR
         PRINT '(A,F15.5)',   '         RMS gradient per image tolerance for constrained potential=',INTRMSTOL
         PRINT '(A,I8)',      '         maximum optimization steps for constrained/real potential=',INTCONSTEPS
         PRINT '(A,I8)',      '         maximum steps for relaxation after adding a new atom before backtrack=',INTRELSTEPS
         PRINT '(A,I6)',      '         maximum number of constraints per atom=',MAXCONUSE
         PRINT '(A,G20.10)',  '         maximum energy per image for convergence during constraint potential phase=',MAXCONE
         PRINT '(A,I6)',      '         interval for checking repulsive interactions=',CHECKREPINTERVAL
         PRINT '(A,F20.10)',  '         multiple of cutoff for repulsion neighbour list=',CHECKREPCUTOFF
         IF (CHECKCONINT) THEN
            PRINT '(A,I6)',      '         adding terms for constraint internal minima'
            PRINT '(A,G20.10)',  '         Internal minima terms will be scaled by a factor of ',INTMINFAC
         ELSE
            PRINT '(A,I6)',      '         not adding terms for constraint internal minima'
         ENDIF
!        IF (DOCROSSCHECK) PRINT '(A,F15.5)','         check for chain crossing of constraints with distance < ',CROSSCUT
      ENDIF
      IF (INTLJT.OR.INTCONSTRAINTT) THEN
         PRINT '(A,2F15.5)','         Minimum and Maximum image separations: ',IMSEPMIN,IMSEPMAX
      ENDIF
      IF (LPERMDIST.OR.LOCALPERMDIST) THEN
          PRINT '(A,G20.10)',' keywords> Cutoff for identifying atoms in the same orbit=',LPDGEOMDIFFTOL
          PRINT '(A,I6,G20.10,I6)',' keywords> lpermdist neighbours, threshold, max sep=', 
     &                       LOCALPERMNEIGH, LOCALPERMCUT, LOCALPERMMAXSEP
      ENDIF

      IF (NEBMAG.GT.0) WRITE(*,'(A,I5)') ' fetchz> NEB magnifications=',NEBMAG
      IF (CONNECTT.AND.(.NOT.NEWCONNECTT)) THEN
         WRITE(*,'(A,I5)') ' fetchz> OLD connect run, maximum transition states=',NCONNECT
         IF (FIXD) WRITE(*,'(A,F12.5)')
     1             ' fetchz> Hard sphere transition state guesses will be applied for minima spearated by less than ',DTHRESH
         IF (STOPFIRST) WRITE(*,'(A)') ' fetchz> Calculation will stop when the initial minimum becomes connected'
      ENDIF
      IF (CALCRATES) THEN
         WRITE(*,'(A,F12.5,A,E20.10)') 
     1          ' fetchz> Rate constants will be calculated for temperature ',TEMPERATURE,' and Plank;s constant=',HRED
         IF (READPATH)  WRITE(*,'(A)') ' fetchz> Rate constants will be calculated for pathway saved in path.info'
      ENDIF

      IF (PATHT) WRITE(*,'(A,I6,A)')' fetchz> Pathways will be calculated saving ',NPATHFRAME,' frames on each side'
      if (machine) write(*,'(a)') ' fetchz> Will use binary files for communication'
      if (machine) write(*,'(a)')' fetchz> WARNING Reading binary files is not fully supported. Thouroughly tested for CHARMM only'
      IF (DUMPPATH) WRITE(*,'(A)')' fetchz> Pathway information will be printed to path.info'
      IF (DUMPDATAT) WRITE(*,'(A)')' fetchz> Stationary point information will be printed to min.data'
      IF (DUMPBESTPATH) WRITE(*,'(A)')' fetchz> Best pathway will be printed to path.info'
      IF (ORDERPARAMT) THEN
         IF (.NOT.DUMPDATAT) THEN
            PRINT '(A)','fetchz> WARNING - order parameter calculation requested, but DUMPDATA not set in odata'
         ELSE
            WRITE(*,'(A,I8)')' fetchz> number of order parameters and derivatives to be printed=',NORDER
            DO J1=1,NORDER
               WRITE(*,'(A,A,I8)')' fetchz> order parameter and additional info=',WHICHORDER(J1),ORDERNUM(J1)
            ENDDO
         ENDIF
      ENDIF
      IF (NEBT.OR.NEWNEBT) WRITE(*,'(A,F12.4)') ' fetchz> Extra (d)neb steps allowed if the (d)neb energy is > lowest value times ',
     &                                          DNEBEFRAC
      IF (BFGSMINT.OR.BFGSTST) THEN
         IF (GMAX.NE.CONVR) THEN
            CONVR=MIN(CONVR,GMAX)
            GMAX=MIN(CONVR,GMAX)
            WRITE(*,'(A,G20.10,A)') ' fetchz> RMS convergence reset to ',GMAX,' for LBFGS and EF phases for consistency'
         ENDIF
         WRITE(*,'(A,G15.8,A,I8)') 
     1           ' fetchz> Convergence criterion for LBFGS optimization: RMS force<',GMAX,' maximum steps=',BFGSSTEPS
         WRITE(*,'(A,G20.10)') ' fetchz> Maximum energy rise in LBFGS minimization=',MAXERISE
         WRITE(*,'(A,I8)') ' fetchz> Number of updates in LBFGS= ',MUPDATE
         IF (BFGSTST) WRITE(*,'(A,I4)') ' fetchz> Number of updates in XLBFGS=',XMUPDATE
         WRITE(*,'(A,I4)') ' fetchz> Number of updates in mind=',MMUPDATE

         IF (NEBT.OR.NEWNEBT) WRITE(*,'(A,I4)') ' fetchz> Number of updates in neb=',NEBMUPDATE
         IF (INTCONSTRAINTT.OR.INTLJT) WRITE(*,'(A,I4)') ' fetchz> Number of updates in qci=',INTMUPDATE
         IF (GROWSTRINGT) WRITE(*,'(A,I4)') ' fetchz> Number of updates in growing string=',GSUPDATE
         IF (GREATCIRCLET) WRITE(*,'(A,I4)') ' fetchz> Number of updates in great circle=',GCUPDATE

         WRITE(*,'(A,F10.4)') ' fetchz> Initial guess for diagonal elements in LBFGS= ',DGUESS
         IF (BFGSTST) WRITE(*,'(A,F10.4)') ' fetchz> Initial guess for diagonal elements in XLBFGS=',XDGUESS
         WRITE(*,'(A,F10.4)') ' fetchz> Maximum step size in LBFGS energy minimization= ',MAXBFGS
         IF (BFGSTST) WRITE(*,'(A,F10.4)') ' fetchz> Maximum step size in XLBFGS=',MAXXBFGS
         WRITE(*,'(A,F10.4)') ' fetchz> Maximum step size in LBFGS neb image minimization=             ',MAXNEBBFGS
      ELSE IF (BSMIN) THEN
         IF (GMAX.GT.CONVR) THEN
            CONVR=GMAX
            WRITE(*,'(A)') ' fetchz> RMS convergence reset to the BS convergence limit'
         ENDIF
         WRITE(*,'(A,G15.8)') ' fetchz> Convergence criterion for BS steepest-descent: RMS force<',GMAX
      ELSE IF (RKMIN) THEN
         IF (GMAX.GT.CONVR) THEN
            CONVR=GMAX
            WRITE(*,'(A)') ' fetchz> RMS convergence reset to the RK convergence limit'
         ENDIF
         WRITE(*,'(A,G15.8)') ' fetchz> Convergence criterion for RK steepest-descent: RMS force<',GMAX
      ENDIF
      IF (REOPT) WRITE(*,'(A)') ' fetchz> Smallest eigenvector will be reconverged after EF step before tangent space minimisation'
      IF (.NOT.(BFGSMINT.OR.BSMIN.OR.RKMIN)) 
     &       WRITE(*,'(A,G15.8,A,G15.8)') ' fetchz> Convergence criteria: EF step<',CONVU,' RMS force<',CONVR
      IF (CHECKINDEX) WRITE(*,'(A)') ' fetchz> Hessian index will be checked'
      IF (CHECKCONT) WRITE(*,'(A)') 
     &    ' fetchz> Search will resume after pushoff if covergence to the wrong Hessian index is detected'
      IF (CHECKNEGATIVET) WRITE(*,'(A)') 
     &    ' fetchz> Step sizes will be scaled if the eigenvalue is positive in BFGSTS'
      IF (CHECKOVERLAPT) WRITE(*,'(A,F15.5,A)') 
     &    ' fetchz> Step sizes will be scaled if the magnitude of the overlap between eigenvectors is <',MINOVERLAP,
     &    ' in BFGSTS'
      IF (REVERSEUPHILLT) WRITE(*,'(A,F15.5,A)') 
     &    ' fetchz> Gradient component will be reversed and projection turned off in tangent space minimisation'
      IF (DCHECK) WRITE(*,'(A)') ' fetchz> Warnings will be issued if atoms become closer than 0.5 units'
      IF (DUMPV) THEN   
         IF (ALLSTEPS) THEN
            IF (ALLVECTORS) THEN
               WRITE(*,'(A)') ' fetchz> All eigenvectors will be dumped at every step'
            ELSE
               WRITE(*,'(A)') ' fetchz> Eigenvector corresponding to smallest eigenvalue will be dumped at every step'
            ENDIF
         ELSE
            IF (ALLVECTORS) THEN
               WRITE(*,'(A)') ' fetchz> All eigenvectors will be dumped at the final step'
            ELSE
               WRITE(*,'(A)') ' fetchz> Eigenvector corresponding to smallest eigenvalue will be dumped at the final step'
            ENDIF
            IF (MWVECTORS) THEN
               WRITE(*,'(A)') ' fetchz> These will be the vectors of the mass weighted hessian, with freqencies in wavenumbers'
            ENDIF
            IF (KTWNT) THEN
               WRITE(*,'(A)') ' fetchz> Effect of thermally accessible modes will be output'
            ENDIF
         ENDIF
      ENDIF
      IF (EVCUT.NE.0.0D0) WRITE(*,'(A,G15.8,A)') ' fetchz> Eigenvalues smaller than ',EVCUT,' will be treated as zero'
      IF (DTEST) WRITE(*,'(A,I3,A,I3,A,G15.8)') 
     1       ' fetchz> GDIIS of maximum dimension ',NDIIA,' will be used every ',NINTV,' steps when RMS force < ',PCUT
      IF (HUPDATE) THEN
         WRITE(*,'(A,F14.5)') ' fetchz> Hessian update procedure will be applied, PHI=',PHIG
         IF (NSTHUP.EQ.0) THEN
            WRITE(*,'(A)') ' SETTINGS No analytic Hessians will be calculated'
         ELSE
            IF (INTHUP.EQ.0) THEN
               WRITE(*,'(A,I4)') ' SETTNIGS Analytic Hessians will be calculated only at step ',NSTHUP
            ELSE
               WRITE(*,'(A,I4,A,I4)') ' SETTNIGS Analytic Hessians will be calculated every ',INTHUP,
     1                                ' steps starting from step ',NSTHUP
            ENDIF
         ENDIF
         IF (INTHUP.EQ.-1) WRITE(*,'(A)') ' fetchz> Hessian will be fixed for all steps'
      ENDIF
      IF (READHESS)  WRITE(*,'(A)') ' fetchz> Hessian will be read from file derivs for initial geometry'
      IF (MASST) WRITE(*,'(A)') ' fetchz> Fictitious kinetic metric will be used'
      IF (.NOT.(BFGSMINT.OR.BSMIN.OR.RKMIN)) THEN
         WRITE(*,'(A,F15.8)') ' fetchz> Initial maximum for EF/SD steps=',MXSTP
         WRITE(*,'(A,F15.8)') ' fetchz> Maximum value for maximum allowed EF/SD steps=',MAXMAX
         WRITE(*,'(A,F15.8)') ' fetchz> Minimum value for maximum allowed EF/SD steps=',MINMAX
         IF (MXSTP.LT.MINMAX) PRINT '(A)','WARNING *** initial maximum step size < minimum'
         IF (MXSTP.GT.MAXMAX) PRINT '(A)','WARNING *** initial maximum step size > maximum'
      ENDIF
      IF (PRESSURE) WRITE(*,'(A)') ' fetchz> Lattice constant will be optimised for zero pressure'
      IF (INR.NE.1) THEN
         IF (.NOT.(BFGSMINT.OR.BSMIN.OR.RKMIN)) 
     1         WRITE(*,'(A,F15.8)') ' fetchz> Value of pushoff from stationary points of the wrong index=',PUSHOFF
         IF (.NOT.(BFGSMINT.OR.BSMIN.OR.RKMIN)) WRITE(*,'(A,F15.8)') 
     1          ' fetchz> A pushoff from stationary points of the wrong index may be applied when the RMS force <',PUSHCUT
      ENDIF
      IF (PUSHOPTT) THEN
         WRITE(6,'(A,G20.10)') 'fetchz> pushoff from ts determined by golden section search in the range zero to +/- ',PUSHOFF
         WRITE(6,'(A,I10,A,G20.10)') '        Maximum golden ratio iterations ',PUSHOPTMAX,' step size convergence ',PUSHOPTCONV  
      ENDIF
      IF (RESIZE.NE.1.0D0) WRITE(*,'(A,F15.8)') ' fetchz> Initial coordinates will be scaled by a factor of ',RESIZE
C     WRITE(*,'(A,E15.8)') ' fetchz> Eigenvalue shifting parameter=',SHIFTV
      IF (FILTH.NE.0) WRITE(*,'(A,I6)') ' fetchz> Number to distinguish output files=',FILTH
      IF (CONTAINER) WRITE(*,'(A,F15.5)') ' fetchz> System will be enclosed in a spherical container radius=',SQRT(RADIUS)
      IF (FIXAFTER.GE.0) WRITE(*,'(A,I6)') ' fetchz> FIXIMAGE will be set permanently on after step ',FIXAFTER
      IF (.NOT.PRINTPTS) WRITE(*,'(A,I6)') ' fetchz> Coordinates for intermediate steps will not be dumped to file points'
      PRINT*

      IF (IPRNT.NE.0) WRITE(*,'(A,I3)') ' fetchz> IPRNT set to ',IPRNT
      IF ((.NOT.BFGSMINT).AND.(.NOT.BFGSTST).AND.(.NOT.BSMIN.OR.RKMIN)) THEN
       IF (PGRAD) WRITE(*,'(A,I3,A)') ' fetchz> Gradients along Hessian eigenvectors will be printed every ',NGRADIENTS,' steps'
       IF (EFSTEPST) WRITE(*,'(A,I4,A)') ' fetchz> Maximum unscaled steps for each mode will be printed every ',EFSTEPS,' steps'
       IF (VALUEST) WRITE(*,'(A,I4,A)') ' fetchz> Hessian eigenvalues will be printed every ',NVALUES,' steps'
       IF (VECTORST) WRITE(*,'(A,I4,A)') ' fetchz> Hessian eigenvectors will be printed every ',NVECTORS,' steps'
       WRITE(*,'(A,I3,A)') ' fetchz> Summary of steps and derivatives will be printed every ',NSUMMARY,' steps'
      ENDIF

      IF (CUDAT) THEN
         WRITE(*,'(A)') ' fetchz> This calculation will take place on the GPU. '
         ! Open and close all output files that might be used - ensures output *replaced* when rerunning CUDAOPTIM
         CUDAFILENAME1 ="GPU_debug_out"
         OPEN(UNIT=325, FILE=CUDAFILENAME1, STATUS="REPLACE")
         CLOSE(325)
         IF (CUDATIMET) THEN
            CUDAFILENAME2 ="GPU_BFGSTS_Rayleigh_Ritz_total.txt"
            CUDAFILENAME3 ="GPU_BFGSTS_subspace_min_total.txt"
            CUDAFILENAME4 ="GPU_BFGSTS_total.txt"
            CUDAFILENAME5 ="GPU_BFGSTS_EF_steps.txt"
            CUDAFILENAME6 ="GPU_potential.txt"
            OPEN(UNIT=326, FILE=CUDAFILENAME2, STATUS="REPLACE")
            CLOSE(326)
            OPEN(UNIT=327, FILE=CUDAFILENAME3, STATUS="REPLACE")
            CLOSE(327)
            OPEN(UNIT=328, FILE=CUDAFILENAME4, STATUS="REPLACE")
            CLOSE(328)
            OPEN(UNIT=329, FILE=CUDAFILENAME5, STATUS="REPLACE")
            CLOSE(329)
            OPEN(UNIT=330, FILE=CUDAFILENAME6, STATUS="REPLACE")
            CLOSE(330)
         END IF
      END IF

      PRINT*

      RETURN
      END
      
