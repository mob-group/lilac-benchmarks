!   Copyright (C) 2003-2010 David J. Wales
!   This file is part of OPTIM.
!
!   OPTIM is free software; you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation; either version 2 of the License, or
!   (at your option) any later version.
!
!   OPTIM is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with this program; if not, write to the Free Software
!   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
SUBROUTINE INTLBFGS(QSTART,QFINISH)
USE PORFUNCS
USE KEY, ONLY : FREEZENODEST, FREEZETOL, MAXINTBFGS, CONVR, ATOMSTORES, &
     & INTRMSTOL, INTIMAGE, NREPMAX, NREPULSIVE, INTMUPDATE, INTDGUESS, &
     & NCONSTRAINT, CONI, CONJ, CONDISTREF, INTCONMAX, CONOFFLIST, CONOFFTRIED,  &
     & INTCONSTRAINREPCUT, REPCON, INTCONSTRAINTREP, INTREPSEP, NREPI, NREPJ, &
     & CONDISTREFLOCAL, INTCONFRAC, CONACTIVE, REPI, &
     & REPJ, NREPMAX, ATOMACTIVE, NCONSTRAINTON, CONION, CONJON, CONDISTREFLOCALON, CONDISTREFON, &
     & NREPCUT, REPCUT, CHECKCONINT, INTCONSTEPS, INTRELSTEPS, MAXCONE, COLDFUSIONLIMIT, &
     & INTSTEPS1, DUMPINTXYZ, DUMPINTXYZFREQ, DUMPINTEOS, DUMPINTEOSFREQ, &
     & IMSEPMIN, IMSEPMAX, MAXINTIMAGE, INTFREEZET, INTFREEZETOL, FREEZE, &
     & INTFROZEN, CHECKREPINTERVAL, NNREPULSIVE, INTFREEZEMIN, INTIMAGECHECK, &
     & CONCUT, CONCUTLOCAL, KINT, REPIFIX, REPJFIX, NREPULSIVEFIX, &
     & NCONSTRAINTFIX, CONIFIX, CONJFIX, QCIPERMCHECK, QCIPERMCHECKINT, BULKT, TWOD, RIGIDBODY, &
     & QCIADDREP, QCIXYZ, WHOLEDNEB, QCIIMAGE, FROZEN, QCIRESTART, NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, SETS, &
     & PERMDIST, LOCALPERMCUT, QCILPERMDIST, QCIPDINT, QCIPERMCUT, QCIAMBERT, BONDS, DOBACK, &
     & QCIRESET, QCIRESETINT1, QCIRESETINT2, JMAXCON, NCONOFF, EREP, ECON, ESPRING, CONVERGECONTEST, CONVERGEREPTEST, &
     & FCONTEST, FREPTEST, QCIKADJUSTTOL, QCIKADJUSTFRAC, QCIKADJUSTFRQ, QCIKINTMAX, QCIKINTMIN, QCIAVDEV, QCISTOP
USE COMMONS, ONLY: NATOMS, DEBUG, PARAM1, PARAM2, PARAM3
USE MODCHARMM, ONLY : CHRMMT
USE CHIRALITY

IMPLICIT NONE 

DOUBLE PRECISION, INTENT(IN) :: QSTART(3*NATOMS), QFINISH(3*NATOMS)  ! The two end points
INTEGER D, U
DOUBLE PRECISION DIST, DIST2, RMAT(3,3), SUMEEE, SUMEEE2, SIGMAEEE, NEIGHBOUR_COORDS(12), CENTRE_COORDS(3)
DOUBLE PRECISION DMAX, DF, DMIN, LOCALSTEP, ADMAX, DUMMYX, DUMMYY, DUMMYZ
INTEGER NDECREASE, NFAIL, NMAXINT, NMININT, JMAX, JMIN, INTIMAGESAVE, NOFF, J1, J2, NQDONE, JA1, JA2, NMOVE, NMOVES, NMOVEF
INTEGER PERM(NATOMS), PERMS(NATOMS), PERMF(NATOMS), STARTGROUP(NPERMGROUP), ENDGROUP(NPERMGROUP)
LOGICAL KNOWE, KNOWG, KNOWH, ADDATOM, ADDREP(NATOMS), LDEBUG, REMOVEIMAGE, PERMUTABLE(NATOMS), IDENTITY, IDONE, TURNOFF
COMMON /KNOWN/ KNOWE, KNOWG, KNOWH

DOUBLE PRECISION DUMMY, DPRAND, DUMMY2, ADUMMY
DOUBLE PRECISION BOXLX,BOXLY,BOXLZ,DISTANCE,RMATBEST(3,3),DISTANCES,DISTANCEF
INTEGER POINT,NPT,J3,J4,NIMAGEFREEZE,NACTIVE,NBEST,NEWATOM,NBEST2,J5,J6
INTEGER TURNONORDER(NATOMS),NBACKTRACK,NQCIFREEZE
INTEGER NDUMMY, NLASTGOODE, NSTEPSMAX, INGROUP(NATOMS), ACID, NLASTCHANGE
LOGICAL CHIRALSR, CHIRALSRP 
INTEGER NTRIES(NATOMS), NITERDONE, EXITSTATUS, DLIST(NATOMS)
DOUBLE PRECISION :: DDOT,STPMIN, ETOTALTMP, RMSTMP, USEFRAC, STIME, FTIME, &
  &                 ETOTAL, LASTGOODE, RMS, STEPTOT, LINTCONSTRAINTTOL, LXYZ(2*3*NATOMS), &
  &                 BESTWORST, WORST, COORDSA(3*NATOMS), COORDSB(3*NATOMS), COORDSC(3*NATOMS)
DOUBLE PRECISION, DIMENSION(INTMUPDATE)     :: RHO1,ALPHA
DOUBLE PRECISION :: EOLD, DMOVED(NATOMS)
LOGICAL SWITCHED, AABACK(NATOMS), BACKDONE
DOUBLE PRECISION, POINTER :: X(:), G(:)
!
! efk: for freezenodes
!
DOUBLE PRECISION :: TESTG, TOTGNORM
INTEGER :: IM
!
! Dimensions involving INTIMAGE
!
DOUBLE PRECISION, ALLOCATABLE :: TRUEEE(:), &
  &              EEETMP(:), MYGTMP(:), EEE(:), STEPIMAGE(:), &
  &              GTMP(:), DIAG(:), STP(:), SEARCHSTEP(:,:), GDIF(:,:), GLAST(:), XSAVE(:)
DOUBLE PRECISION, ALLOCATABLE :: VPLUS(:), VMINUS(:)   
DOUBLE PRECISION  EPLUS, EMINUS, DIFF   
DOUBLE PRECISION, ALLOCATABLE, TARGET :: XYZ(:), GGG(:), DPTMP(:), D2TMP(:,:)
! saved interpolation
INTEGER BESTINTIMAGE, NSTEPS, NITERUSE
LOGICAL, ALLOCATABLE :: CHECKG(:), IMGFREEZE(:)
INTEGER, ALLOCATABLE :: NCONATOM(:), CONLIST(:,:), COMMONCON(:,:)
LOGICAL READIMAGET, GROUPACTIVE(NPERMGROUP)
INTEGER NCONCOMMON(NPERMGROUP)
INTEGER LUNIT, GETUNIT, NCOMMONCON
CHARACTER(LEN=2) SDUMMY
INTEGER JMAXEEE,JMAXRMS,num_chiral_centres,atom_number,MAXCONSTRAINTS,PATOM1,PATOM2,PATOMTEST,NCOMMONMAX 
DOUBLE PRECISION MAXEEE,MAXRMS,MINEEE,SAVELOCALPERMCUT

WHOLEDNEB=.FALSE.
READIMAGET=.FALSE.
REMOVEIMAGE=.FALSE.
ECON=0.0D0; EREP=0.0D0; ESPRING=0.0D0

IF (QCIAMBERT) THEN ! copied from corresponding chirality subroutine

   num_chiral_centres=SIZE(sr_atoms,1)
   WRITE(*,'(A,I8)') ' intlbfgs> Number of chiral sites=',num_chiral_centres

!  do J1 = 1, num_chiral_centres
!     write(*, '(5i8)') sr_atoms(J1, :)
!     write(*,'(A,L5)') 'sr_states_initial(J1) = ',sr_states_initial(J1)
!  end do

ENDIF

! DO J1=1,NATOMS
!    WRITE(*,'(A,2I8)') 'intlbfgs> atom and residue: ',J1,ATOMSTORES(J1)
! ENDDO
! STOP
NCONOFF=0
AABACK(1:NATOMS)=.FALSE.
BACKDONE=.FALSE.
IF (DOBACK) THEN
   LUNIT=GETUNIT()
   OPEN(UNIT=LUNIT,FILE='aabk',STATUS='OLD')
   DO J1=1,NATOMS
      READ(LUNIT,*,END=861) NDUMMY
      AABACK(NDUMMY)=.TRUE.
   ENDDO
861   CLOSE(LUNIT)
   WRITE(*,'(A)') 'intlbfgs> Backbone atoms from file aabk:'
   WRITE(*,'(20L5)') AABACK(1:NATOMS)
ENDIF

IF (QCIRESTART) READIMAGET=.TRUE.
IF (READIMAGET) THEN
   LUNIT=GETUNIT()
   OPEN(UNIT=LUNIT,FILE='int.xyz',STATUS='OLD')
   INTIMAGE=0
653 CONTINUE
   READ(LUNIT,*,END=654) NDUMMY
   READ(LUNIT,*) 
   DO J1=1,NATOMS
      READ(LUNIT,*) SDUMMY, DUMMYX, DUMMYY, DUMMYZ
!     WRITE(*,'(A,I6,A2,3G20.10)') 'J1,sd,xd,yd,zd=',J1,SDUMMY, DUMMYX, DUMMYY, DUMMYZ
   ENDDO
   INTIMAGE=INTIMAGE+1
   GOTO 653
654 CONTINUE
   INTIMAGE=INTIMAGE-2
   WRITE(*,'(A,I10,A)') 'intlbfgs> Rereading ',INTIMAGE,' frames'
   CLOSE(LUNIT)
ENDIF

ALLOCATE(TRUEEE(INTIMAGE+2), &
  &      EEETMP(INTIMAGE+2), MYGTMP(3*NATOMS*INTIMAGE), &
  &      GTMP(3*NATOMS*INTIMAGE), &
  &      DIAG(3*NATOMS*INTIMAGE), STP(3*NATOMS*INTIMAGE), SEARCHSTEP(0:INTMUPDATE,(3*NATOMS)*INTIMAGE), &
  &      GDIF(0:INTMUPDATE,(3*NATOMS)*INTIMAGE),GLAST((3*NATOMS)*INTIMAGE), XSAVE((3*NATOMS)*INTIMAGE), &
  &      XYZ((3*NATOMS)*(INTIMAGE+2)), GGG((3*NATOMS)*(INTIMAGE+2)), CHECKG((3*NATOMS)*INTIMAGE), IMGFREEZE(INTIMAGE), &
  &      EEE(INTIMAGE+2), STEPIMAGE(INTIMAGE))

ALLOCATE(VPLUS((3*NATOMS)*(INTIMAGE+2)),VMINUS((3*NATOMS)*(INTIMAGE+2)))  

SWITCHED=.FALSE.
INTIMAGESAVE=INTIMAGE
NBACKTRACK=1
CALL MYCPU_TIME(STIME,.FALSE.)
WRITE(*,'(A,I6)') ' intlbfgs> Maximum number of steps for constraint potential phase is ',INTSTEPS1
WRITE(*,'(A,I6,A,G20.10)') ' intlbfgs> Updates: ',INTMUPDATE,' maximum step size=',MAXINTBFGS
ADDATOM=.FALSE.
NFAIL=0
IMGFREEZE(1:INTIMAGE)=.FALSE.
D=(3*NATOMS)*INTIMAGE
U=INTMUPDATE
NITERDONE=1
NITERUSE=1
NQDONE=0

IF ( D<=0 ) THEN
   WRITE(*,*) 'd is not positive, d=',d
   STOP
ENDIF
IF ( U<=0 ) THEN
   WRITE(*,*) 'u is not positive, u=',u
   STOP
ENDIF
IF (INTSTEPS1 < 0) THEN
   WRITE(*,'(1x,a)') 'Maximal number of iterations is less than zero! Stop.'
   STOP
ENDIF
!
! XYZ, GGG, EEE include the end point images
! X, G do not.
!
IF (.NOT.ALLOCATED(CONI)) THEN 
   ALLOCATE(CONI(INTCONMAX),CONJ(INTCONMAX),CONDISTREF(INTCONMAX),CONCUT(INTCONMAX),CONOFFLIST(INTCONMAX),CONOFFTRIED(INTCONMAX))
   CONOFFTRIED(1:INTCONMAX)=.FALSE.
   ALLOCATE(REPI(NREPMAX),REPJ(NREPMAX),NREPI(NREPMAX),NREPJ(NREPMAX),REPCUT(NREPMAX),NREPCUT(NREPMAX))
ENDIF
X=>XYZ((3*NATOMS)+1:(3*NATOMS)*(INTIMAGE+1))
G=>GGG((3*NATOMS)+1:(3*NATOMS)*(INTIMAGE+1))
!
! Initialise XYZ
!
IF (READIMAGET) THEN  ! Note that this will ignore the coordinates in start and finish
   LUNIT=GETUNIT()
   OPEN(UNIT=LUNIT,FILE='int.xyz',STATUS='OLD')
   DO J2=1,INTIMAGE+2
      READ(LUNIT,*) NDUMMY
      READ(LUNIT,*) 
      DO J1=1,NATOMS
         READ(LUNIT,*) SDUMMY,XYZ(3*NATOMS*(J2-1)+3*(J1-1)+1),XYZ(3*NATOMS*(J2-1)+3*(J1-1)+2),XYZ(3*NATOMS*(J2-1)+3*(J1-1)+3)
      ENDDO
   ENDDO
   CLOSE(LUNIT)
!
! Don't overwrite start and finish - we should have aligned these by now
!
   XYZ(1:(3*NATOMS))=QSTART(1:(3*NATOMS))
   XYZ((3*NATOMS)*(INTIMAGE+1)+1:(3*NATOMS)*(INTIMAGE+2))=QFINISH(1:(3*NATOMS))
ELSE
   XYZ(1:(3*NATOMS))=QSTART(1:(3*NATOMS))
   XYZ((3*NATOMS)*(INTIMAGE+1)+1:(3*NATOMS)*(INTIMAGE+2))=QFINISH(1:(3*NATOMS))
   DO J1=1,INTIMAGE+2
      XYZ((J1-1)*(3*NATOMS)+1:J1*(3*NATOMS))=((INTIMAGE+2-J1)*QSTART(1:(3*NATOMS))+(J1-1)*QFINISH(1:(3*NATOMS)))/(INTIMAGE+1)
   ENDDO
ENDIF

NQCIFREEZE=0
! IF (FREEZE) THEN
!    WRITE(*,'(A)') ' intlbfgs> ERROR *** QCI has not been coded for frozen atoms yet'
!    STOP     
! ENDIF
IF (ALLOCATED(INTFROZEN)) DEALLOCATE(INTFROZEN)
ALLOCATE(INTFROZEN(NATOMS))
INTFROZEN(1:NATOMS)=.FALSE.
DLIST(1:NATOMS)=-1
DMOVED(1:NATOMS)=1.0D100
IF (INTFREEZET) THEN
   DUMMY=INTFREEZETOL**2
   DO J1=1,NATOMS
      DF=(XYZ(3*(J1-1)+1)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+1))**2 &
  &     +(XYZ(3*(J1-1)+2)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+2))**2 &
  &     +(XYZ(3*(J1-1)+3)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+3))**2
      IF ((DF.LT.DUMMY).OR.FROZEN(J1)) THEN
         NQCIFREEZE=NQCIFREEZE+1
         INTFROZEN(J1)=.TRUE.
         IF (DEBUG) WRITE(*,'(A,I6,A,F12.6,A,I6)') &
  &            ' intlbfgs> atom ',J1,' moves less than threshold: dist^2=',DF,' total=',NQCIFREEZE
      ENDIF
      sortd: DO J2=1,J1
         IF (DF.LT.DMOVED(J2)) THEN
            DO J3=J1,J2+1,-1
               DMOVED(J3)=DMOVED(J3-1)
               DLIST(J3)=DLIST(J3-1)
            ENDDO
            DMOVED(J2)=DF
            DLIST(J2)=J1
            EXIT sortd
         ENDIF
      ENDDO sortd
   ENDDO
   WRITE(*,'(A,I6,A,F12.6,A,I6)') ' intlbfgs> Total number of atoms moving less than threshold=',NQCIFREEZE
ENDIF

IF (NATOMS-NQCIFREEZE.LT.INTFREEZEMIN) THEN
   DO J1=NATOMS,NATOMS-INTFREEZEMIN+1,-1
      INTFROZEN(DLIST(J1))=.FALSE.
   ENDDO
   NQCIFREEZE=MAX(0,NATOMS-INTFREEZEMIN)
   WRITE(*,'(A,I6,A)') ' intlbfgs> Freezing ',NQCIFREEZE,' atoms'
ENDIF

NLASTGOODE=0
NLASTCHANGE=0
LASTGOODE=1.0D100

!
! Constraints are collected in a list and activated via the CONACTIVE(J1)
! logical array. There will generally be of order NATOMS. However, the
! repulsions will scale as NATOMS**2 and are treated differently. The
! active repulsions are stored sequentially as atoms are added to the
! growing list. This is done even if we have congeom or congeom.dat files
! available. In this case we use the fixed list of possible constraints
! via CHECKPERC, but the list of repulsions and cutoffs is recreated on
! the fly. The fixed lists are used in make_conpot, since this is called
! for pairs of minima with all atoms active to obtain an interpolation
! metric.
!
! Perhaps we should use the fixed list to activate the repulsions below?
! A neighbour list for repulsions is maintained to make the constraint
! potential evaluation scale as order N.
!
IF (NQCIFREEZE.LT.NATOMS) THEN
   LXYZ(1:(3*NATOMS))=QSTART(1:(3*NATOMS))
   LXYZ((3*NATOMS)+1:2*(3*NATOMS))=QFINISH(1:(3*NATOMS))
   CALL CHECKPERC(LXYZ,LINTCONSTRAINTTOL,NQCIFREEZE,2)
ELSE
   IF (.NOT.ALLOCATED(ATOMACTIVE)) ALLOCATE(ATOMACTIVE(NATOMS))
   NCONSTRAINT=0
   WRITE(*,'(A)') ' intlbfgs> All atoms move less than threshold - skip to linear interpolation for end points'
   INTIMAGE=0
   XYZ(1:(3*NATOMS))=QSTART(1:(3*NATOMS))
   XYZ((3*NATOMS)*(INTIMAGE+1)+1:(3*NATOMS)*(INTIMAGE+2))=QFINISH(1:(3*NATOMS))
   DO J1=1,INTIMAGE+2
      XYZ((J1-1)*(3*NATOMS)+1:J1*(3*NATOMS))=((INTIMAGE+2-J1)*QSTART(1:(3*NATOMS))+(J1-1)*QFINISH(1:(3*NATOMS)))/(INTIMAGE+1)
   ENDDO
   GOTO 678
ENDIF

NACTIVE=0
TURNONORDER(1:NATOMS)=0
ATOMACTIVE(1:NATOMS)=.FALSE.
REPCON=-INTCONSTRAINTREP/INTCONSTRAINREPCUT**6 ! also needed for congrad.f90 potential
IF (QCIRESTART) THEN
   LUNIT=GETUNIT()
   OPEN(UNIT=LUNIT,FILE='QCIdump',STATUS='OLD')

   READ(LUNIT,*) NACTIVE
   WRITE(*,'(A,I10,A)') ' intlbfgs> restart has ',NACTIVE,' active atoms'
   READ(LUNIT,*) KINT, INTCONSTRAINTREP
   WRITE(*,'(A,2G20.10)') ' intlbfgs> spring constant and repulsive prefactor: ',KINT, INTCONSTRAINTREP
   WRITE(*,'(A)') ' intlbfgs> reading turnonorder'
   READ(LUNIT,*) TURNONORDER(1:NACTIVE)
!  WRITE(*,'(12I8)') TURNONORDER(1:NACTIVE)
   WRITE(*,'(A)') ' intlbfgs> reading atomactive'
   READ(LUNIT,*) ATOMACTIVE(1:NATOMS)
!  WRITE(*,'(12L5)') ATOMACTIVE(1:NATOMS) 
!  WRITE(*,'(A,3L5)') 'intlbfgs> here A atomactive 2113 2115 2123=',ATOMACTIVE(2113),ATOMACTIVE(2115),ATOMACTIVE(2123)
   READ(LUNIT,*) NCONSTRAINT
   WRITE(*,'(A)') ' intlbfgs> reading conactive'
   WRITE(*,'(A,I10,A)') ' intlbfgs> restart has ',NCONSTRAINT,' constraints'
   ALLOCATE(CONACTIVE(NCONSTRAINT))
   READ(LUNIT,*) CONACTIVE(1:NCONSTRAINT)
!  WRITE(*,'(12L5)') CONACTIVE(1:NCONSTRAINT) 

   ALLOCATE(CONDISTREFLOCAL(NCONSTRAINT))
   ALLOCATE(CONCUTLOCAL(NCONSTRAINT))
   CONDISTREFLOCAL(1:NCONSTRAINT)=CONDISTREF(1:NCONSTRAINT)
   CONCUTLOCAL(1:NCONSTRAINT)=CONCUT(1:NCONSTRAINT)

   READ(LUNIT,*) NREPULSIVE,NNREPULSIVE,NREPMAX
   USEFRAC=1.0D0
!  WRITE(*,'(A,3I10,G20.10)') 'intlbfgs> NREPULSIVE,NNREPULSIVE,NREPMAX=',NREPULSIVE,NNREPULSIVE,NREPMAX
   IF (ALLOCATED(REPI)) DEALLOCATE(REPI)
   IF (ALLOCATED(REPJ)) DEALLOCATE(REPJ)
   IF (ALLOCATED(NREPI)) DEALLOCATE(NREPI)
   IF (ALLOCATED(NREPJ)) DEALLOCATE(NREPJ)
   IF (ALLOCATED(REPCUT)) DEALLOCATE(REPCUT)
   IF (ALLOCATED(NREPCUT)) DEALLOCATE(NREPCUT)
   ALLOCATE(REPI(NREPMAX),REPJ(NREPMAX),NREPI(NREPMAX),NREPJ(NREPMAX),REPCUT(NREPMAX),NREPCUT(NREPMAX))
   READ(LUNIT,*) REPI(1:NREPULSIVE)
   WRITE(*,'(A)') ' intlbfgs> read REPI:'
!  WRITE(*,'(12I8)') REPI(1:NREPULSIVE)
   READ(LUNIT,*) REPJ(1:NREPULSIVE)
   WRITE(*,'(A)') ' intlbfgs> read REPJ:'
!  WRITE(*,'(12I8)') REPJ(1:NREPULSIVE)
   READ(LUNIT,*) NREPI(1:NNREPULSIVE)
   WRITE(*,'(A)') ' intlbfgs> read NREPI:'
!  WRITE(*,'(12I8)') NREPI(1:NNREPULSIVE)
   READ(LUNIT,*) NREPJ(1:NNREPULSIVE)
   WRITE(*,'(A)') ' intlbfgs> read NREPJ:'
!  WRITE(*,'(12I8)') NREPJ(1:NNREPULSIVE)
   READ(LUNIT,*) REPCUT(1:NREPULSIVE)
   WRITE(*,'(A)') ' intlbfgs> read REPCUT:'
!  WRITE(*,'(6G20.10)') REPCUT(1:NREPULSIVE)
   READ(LUNIT,*) NREPCUT(1:NNREPULSIVE)
   WRITE(*,'(A)') ' intlbfgs> read NREPCUT:'
!  WRITE(*,'(6G20.10)') NREPCUT(1:NNREPULSIVE)

   READ(LUNIT,*) INTFROZEN(1:NATOMS)
   WRITE(*,'(A)') ' intlbfgs> read INTFROZEN'
!  WRITE(*,'(12L5)') INTFROZEN(1:NATOMS)

   NCONOFF=0
   READ(LUNIT,*,END=742) NCONOFF
   IF (NCONOFF.GT.0) READ(LUNIT,*) CONOFFLIST(1:NCONOFF)
   IF (NCONOFF.GT.0) READ(LUNIT,*) CONOFFTRIED(1:NCONSTRAINT)
742 CONTINUE
   CLOSE(LUNIT)

   GLAST(1:D)=G(1:D)
   XSAVE(1:D)=X(1:D)
   GOTO 986
ENDIF
IF (INTFREEZET) THEN
   DO J1=1,NATOMS
      IF (INTFROZEN(J1)) THEN
! 
! linear interpolation 
! 
         DO J2=2,INTIMAGE+1
            XYZ((J2-1)*3*NATOMS+3*(J1-1)+1:(J2-1)*3*NATOMS+3*(J1-1)+3)= &
  &            (INTIMAGE-J2+2)*XYZ(3*(J1-1)+1:3*(J1-1)+3)/(INTIMAGE+1) &
  &           +(J2-1)*XYZ(3*NATOMS*(INTIMAGE+1)+3*(J1-1)+1:3*NATOMS*(INTIMAGE+1)+3*(J1-1)+3)/(INTIMAGE+1)
         ENDDO
         ATOMACTIVE(J1)=.TRUE.
         NACTIVE=NACTIVE+1
         TURNONORDER(NACTIVE)=J1
         NTRIES(J1)=1
      ENDIF
   ENDDO
ENDIF

ALLOCATE(NCONATOM(NATOMS))
NCONATOM(1:NATOMS)=0
DO J1=1,NCONSTRAINT
   NCONATOM(CONI(J1))=NCONATOM(CONI(J1))+1
   NCONATOM(CONJ(J1))=NCONATOM(CONJ(J1))+1
ENDDO
MAXCONSTRAINTS=-1
DO J1=1,NATOMS
   IF (NCONATOM(J1).GT.MAXCONSTRAINTS) THEN
      MAXCONSTRAINTS=NCONATOM(J1)
      J2=J1
   ENDIF
ENDDO
WRITE(*,'(A,I6,A,I6)') ' intlbfgs> maximum constraints ',MAXCONSTRAINTS,' for atom ',J2
ALLOCATE(CONLIST(NATOMS,MAXCONSTRAINTS))
CONLIST(1:NATOMS,1:MAXCONSTRAINTS)=0
NCONATOM(1:NATOMS)=0
DO J1=1,NCONSTRAINT
   NCONATOM(CONI(J1))=NCONATOM(CONI(J1))+1
   NCONATOM(CONJ(J1))=NCONATOM(CONJ(J1))+1
   CONLIST(CONI(J1),NCONATOM(CONI(J1)))=CONJ(J1)
   CONLIST(CONJ(J1),NCONATOM(CONJ(J1)))=CONI(J1)
ENDDO

DO J1=1,NATOMS
   WRITE(*,'(A,I6,A,20I6)') ' intlbfgs> atom ',J1,' constraints: ',CONLIST(J1,1:NCONATOM(J1))
ENDDO

NDUMMY=1
NCOMMONMAX=-1
DO J1=1,NPERMGROUP
   NCONCOMMON(J1)=0
   PATOM1=PERMGROUP(NDUMMY)
!  WRITE(*,'(A,I6,A,I6,A,I6)') 'group ',J1,' atom ',PATOM1,' constraints=',NCONATOM(PATOM1)
!  WRITE(*,'(20I6)') CONLIST(PATOM1,1:NCONATOM(PATOM1))
!
! For each entry in constraint list of first permutable atom, check if it exists for the second, 
! if so, check the third, etc.
!
   atlist: DO J4=1,NCONATOM(PATOM1)
      PATOMTEST=CONLIST(PATOM1,J4)
      plist: DO J5=2,NPERMSIZE(J1)
         PATOM2=PERMGROUP(NDUMMY+J5-1)
         DO J6=1,NCONATOM(PATOM2)
            IF (CONLIST(PATOM2,J6).EQ.PATOMTEST) CYCLE plist
         ENDDO
         CYCLE atlist
      ENDDO plist
      NCONCOMMON(J1)=NCONCOMMON(J1)+1
!     WRITE(*,'(4(A,I6))') 'atom ',PATOMTEST,' is a common constraint for permgroup ',J1,' total=',NCONCOMMON(J1),' lists are:'  
      DO J5=1,NPERMSIZE(J1)
         J6=PERMGROUP(NDUMMY+J5-1)
!        WRITE(*,'(A,I6,A,20I6)') 'atom ',J6,' constraints: ',CONLIST(J6,1:NCONATOM(J6))
      ENDDO
   ENDDO atlist
!  WRITE(*,'(A,I6,A,I6,A,I6)') 'group ',J1,' size ',NPERMSIZE(J1),' common constraints ',NCONCOMMON(J1)
   NDUMMY=NDUMMY+NPERMSIZE(J1)
   IF (NCONCOMMON(J1).GT.NCOMMONMAX) NCOMMONMAX=NCONCOMMON(J1)
ENDDO
ALLOCATE(COMMONCON(NPERMGROUP,NCOMMONMAX))

WRITE(*,'(A,I6)') 'largest number of common constraint atoms for any group is: ',NCOMMONMAX

!
! Now repeat and save the common constrained atoms in COMMONCON(J1,1:NCOMMONCON(J1)) for permutational group J1.
!

NDUMMY=1
DO J1=1,NPERMGROUP
   NCONCOMMON(J1)=0
   PATOM1=PERMGROUP(NDUMMY)
!  WRITE(*,'(20I6)') CONLIST(PATOM1,1:NCONATOM(PATOM1))
!
! For each entry in constraint list of first permutable atom, check if it exists for the second, 
! if so, check the third, etc.
!
   atlist2: DO J4=1,NCONATOM(PATOM1)
      PATOMTEST=CONLIST(PATOM1,J4)
      plist2: DO J5=2,NPERMSIZE(J1)
         PATOM2=PERMGROUP(NDUMMY+J5-1)
         DO J6=1,NCONATOM(PATOM2)
            IF (CONLIST(PATOM2,J6).EQ.PATOMTEST) CYCLE plist2
         ENDDO
         CYCLE atlist2
      ENDDO plist2
      NCONCOMMON(J1)=NCONCOMMON(J1)+1
!     WRITE(*,'(4(A,I6))') 'atom ',PATOMTEST,' is a common constraint for permgroup ',J1,' total=',NCONCOMMON(J1),' lists are:'  
      COMMONCON(J1,NCONCOMMON(J1))=PATOMTEST
      DO J5=1,NPERMSIZE(J1)
         J6=PERMGROUP(NDUMMY+J5-1)
!        WRITE(*,'(A,I6,A,20I6)') 'atom ',J6,' constraints: ',CONLIST(J6,1:NCONATOM(J6))
      ENDDO
   ENDDO atlist2
   WRITE(*,'(A,I6,A,I6,A,20I6)') 'group ',J1,' size ',NPERMSIZE(J1),' common constraints to atoms ',COMMONCON(J1,1:NCONCOMMON(J1))
   NDUMMY=NDUMMY+NPERMSIZE(J1)
ENDDO

REPCON=-INTCONSTRAINTREP/INTCONSTRAINREPCUT**6 ! also needed for congrad.f90 potential
IF (ALLOCATED(CONDISTREFLOCAL)) DEALLOCATE(CONDISTREFLOCAL)
IF (ALLOCATED(CONCUTLOCAL)) DEALLOCATE(CONCUTLOCAL)
ALLOCATE(CONDISTREFLOCAL(NCONSTRAINT))
ALLOCATE(CONCUTLOCAL(NCONSTRAINT))
IF (ALLOCATED(CONDISTREFLOCALON)) DEALLOCATE(CONDISTREFLOCALON)
IF (ALLOCATED(CONDISTREFON)) DEALLOCATE(CONDISTREFON)
IF (ALLOCATED(CONION)) DEALLOCATE(CONION)
IF (ALLOCATED(CONJON)) DEALLOCATE(CONJON)
ALLOCATE(CONDISTREFLOCALON(NCONSTRAINT),CONDISTREFON(NCONSTRAINT),CONION(NCONSTRAINT),CONJON(NCONSTRAINT))
CONDISTREFLOCAL(1:NCONSTRAINT)=CONDISTREF(1:NCONSTRAINT)
CONCUTLOCAL(1:NCONSTRAINT)=CONCUT(1:NCONSTRAINT)
DUMMY=1.0D100
DUMMY2=-1.0D100
IF (NCONSTRAINT.EQ.0) THEN
   NACTIVE=NATOMS
   EOLD=ETOTAL
   SWITCHED=.TRUE.
   USEFRAC=1.0D0
   NREPULSIVE=0
   NNREPULSIVE=0
   GLAST(1:D)=G(1:D)
   XSAVE(1:D)=X(1:D)
   GOTO 567
ENDIF
DO J1=1,NCONSTRAINT
   IF (DOBACK.AND.(.NOT.AABACK(CONI(J1)).OR.(.NOT.AABACK(CONJ(J1))))) CYCLE
   DF=SQRT((XYZ(3*(CONI(J1)-1)+1)-XYZ((INTIMAGE+1)*3*NATOMS+3*(CONI(J1)-1)+1))**2 &
  &       +(XYZ(3*(CONI(J1)-1)+2)-XYZ((INTIMAGE+1)*3*NATOMS+3*(CONI(J1)-1)+2))**2 &
  &       +(XYZ(3*(CONI(J1)-1)+3)-XYZ((INTIMAGE+1)*3*NATOMS+3*(CONI(J1)-1)+3))**2)&
  &  +SQRT((XYZ(3*(CONJ(J1)-1)+1)-XYZ((INTIMAGE+1)*3*NATOMS+3*(CONJ(J1)-1)+1))**2 &
  &       +(XYZ(3*(CONJ(J1)-1)+2)-XYZ((INTIMAGE+1)*3*NATOMS+3*(CONJ(J1)-1)+2))**2 &
  &       +(XYZ(3*(CONJ(J1)-1)+3)-XYZ((INTIMAGE+1)*3*NATOMS+3*(CONJ(J1)-1)+3))**2)
   IF (DF.LT.DUMMY) THEN
      NBEST=J1
      DUMMY=DF
   ENDIF
   IF (DF.GT.DUMMY2) THEN
      NBEST2=J1
      DUMMY2=DF
   ENDIF
ENDDO
IF (DEBUG) WRITE(*,'(A,I6,A,2I6,A,F15.5)') ' intlbfgs> Smallest overall motion for constraint ',NBEST, ' atoms ', &
  &                           CONI(NBEST),CONJ(NBEST),' distance=',DUMMY
IF (DEBUG) WRITE(*,'(A,I6,A,2I6,A,F15.5)') ' intlbfgs> Largest overall motion for constraint  ',NBEST2,' atoms ', &
  &                           CONI(NBEST2),CONJ(NBEST2),' distance=',DUMMY2

!!! NBEST=NBEST2 !!!! DJW
NTRIES(1:NATOMS)=1
IF (ALLOCATED(CONACTIVE)) DEALLOCATE(CONACTIVE)
ALLOCATE(CONACTIVE(NCONSTRAINT))
CONACTIVE(1:NCONSTRAINT)=.FALSE.
CONACTIVE(NBEST)=.TRUE.
ATOMACTIVE(CONI(NBEST))=.TRUE.
ATOMACTIVE(CONJ(NBEST))=.TRUE.
IF (.NOT.INTFROZEN(CONI(NBEST))) THEN
   TURNONORDER(NACTIVE+1)=CONI(NBEST)
   NACTIVE=NACTIVE+1
ENDIF
IF (.NOT.INTFROZEN(CONJ(NBEST))) THEN
   TURNONORDER(NACTIVE+1)=CONJ(NBEST)
   NACTIVE=NACTIVE+1
ENDIF
NTRIES(CONI(NBEST))=1
NTRIES(CONJ(NBEST))=1
NREPULSIVE=0
NCONSTRAINTON=1
CONDISTREFLOCALON(1)=CONDISTREFLOCAL(NBEST)
CONDISTREFON(1)=CONDISTREF(NBEST)
CONION(1)=CONI(NBEST)
CONJON(1)=CONJ(NBEST)
IF (DEBUG) WRITE(*,'(A,I6)') ' intlbfgs> Number of active atoms is now ',NACTIVE
!
! If INTFREEZET is true we need to add constraints and replusions to the frozen atoms.
! ATOMACTIVE is .TRUE. for frozen atoms.
!
IF (INTFREEZET) THEN
   DO J1=1,NCONSTRAINT
      IF (CONACTIVE(J1)) CYCLE
      IF ((CONI(J1).EQ.CONI(NBEST)).AND.(ATOMACTIVE(CONJ(J1))).OR.(CONJ(J1).EQ.CONI(NBEST)).AND.(ATOMACTIVE(CONI(J1)))) THEN
         CONACTIVE(J1)=.TRUE.
         IF (DEBUG) WRITE(*,'(A,I6,A,2I6)') ' intlbfgs> Turning on constraint ',J1,' for atoms ',CONI(J1),CONJ(J1)
      ENDIF
      IF ((CONI(J1).EQ.CONJ(NBEST)).AND.(ATOMACTIVE(CONJ(J1))).OR.(CONJ(J1).EQ.CONJ(NBEST)).AND.(ATOMACTIVE(CONI(J1)))) THEN
         CONACTIVE(J1)=.TRUE.
         IF (DEBUG) WRITE(*,'(A,I6,A,2I6)') ' intlbfgs> Turning on constraint ',J1,' for atoms ',CONI(J1),CONJ(J1)
      ENDIF
   ENDDO

   DO J1=1,NATOMS
      IF (.NOT.ATOMACTIVE(J1)) CYCLE ! identify active atoms
      IF (ABS(J1-CONI(NBEST)).LE.INTREPSEP) CYCLE ! no repulsion for atoms too close in sequence
      IF (INTFROZEN(J1).AND.INTFROZEN(CONI(NBEST))) CYCLE
      DO J2=1,NCONSTRAINT
!
!  With MAXCONUSE set to a finite value there could be constraints for the new atom that are
!  not active. We don't want these to be changed to repulsion, surely?!
!  Or perhaps we do need to do something with them?
!
         IF (.NOT.CONACTIVE(J2)) CYCLE ! repulsions for constraints
         IF (((CONI(J2).EQ.J1).AND.(CONJ(J2).EQ.CONI(NBEST))).OR.((CONJ(J2).EQ.J1).AND.(CONI(J2).EQ.CONI(NBEST)))) GOTO 545
      ENDDO
      DMIN=1.0D100
      DMAX=-1.0D0
      DO J2=1,INTIMAGE+2,INTIMAGE+1 ! only consider the end-point distances
         DF=SQRT((XYZ((J2-1)*3*NATOMS+3*(CONI(NBEST)-1)+1)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+1))**2+ &
  &           (XYZ((J2-1)*3*NATOMS+3*(CONI(NBEST)-1)+2)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+2))**2+ &
  &           (XYZ((J2-1)*3*NATOMS+3*(CONI(NBEST)-1)+3)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+3))**2)
         IF (DF.GT.DMAX) DMAX=DF
         IF (DF.LT.DMIN) DMIN=DF
      ENDDO
!
! Use the minimum of the end point distances and INTCONSTRAINREPCUT for each contact.
!
      DMIN=MIN(DMIN-1.0D-3,INTCONSTRAINREPCUT)
      NREPULSIVE=NREPULSIVE+1
      IF (NREPULSIVE.GT.NREPMAX) CALL REPDOUBLE
      REPI(NREPULSIVE)=J1
      REPJ(NREPULSIVE)=CONI(NBEST)
      REPCUT(NREPULSIVE)=DMIN
      IF (DEBUG) WRITE(*,'(A,I6,A,I6,A,F15.5)') ' intlbfgs> Adding repulsion for new atom ',CONI(NBEST),' with atom ',J1, &
  &                                          ' cutoff=',DMIN
545   CONTINUE
   ENDDO

   DO J1=1,NATOMS
      IF (.NOT.ATOMACTIVE(J1)) CYCLE ! identify active atoms
      IF (ABS(J1-CONJ(NBEST)).LE.INTREPSEP) CYCLE ! no repulsion for atoms too close in sequence
      IF (INTFROZEN(J1).AND.INTFROZEN(CONJ(NBEST))) CYCLE
      DO J2=1,NCONSTRAINT
!
!  With MAXCONUSE set to a finite value there could be constraints for the new atom that are
!  not active. We don't want these to be changed to repulsion, surely?!
!  Or perhaps we do need to do something with them?
!
         IF (.NOT.CONACTIVE(J2)) CYCLE ! identify active constraints
         IF (((CONI(J2).EQ.J1).AND.(CONJ(J2).EQ.CONJ(NBEST))).OR.((CONJ(J2).EQ.J1).AND.(CONI(J2).EQ.CONJ(NBEST)))) GOTO 541
      ENDDO
      DMIN=1.0D100
      DMAX=-1.0D0
      DO J2=1,INTIMAGE+2,INTIMAGE+1 ! only consider the end-point distances
         DF=SQRT((XYZ((J2-1)*3*NATOMS+3*(CONJ(NBEST)-1)+1)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+1))**2+ &
  &           (XYZ((J2-1)*3*NATOMS+3*(CONJ(NBEST)-1)+2)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+2))**2+ &
  &           (XYZ((J2-1)*3*NATOMS+3*(CONJ(NBEST)-1)+3)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+3))**2)
         IF (DF.GT.DMAX) DMAX=DF
         IF (DF.LT.DMIN) DMIN=DF
      ENDDO
!
! Use the minimum of the end point distances and INTCONSTRAINREPCUT for each contact.
!
      DMIN=MIN(DMIN-1.0D-3,INTCONSTRAINREPCUT)
      NREPULSIVE=NREPULSIVE+1
      IF (NREPULSIVE.GT.NREPMAX) CALL REPDOUBLE
      REPI(NREPULSIVE)=J1
      REPJ(NREPULSIVE)=CONJ(NBEST)
      REPCUT(NREPULSIVE)=DMIN
      IF (DEBUG) WRITE(*,'(A,I6,A,I6,A,F15.5)') ' intlbfgs> Adding repulsion for new atom ',CONJ(NBEST),' with atom ',J1, &
  &                                          ' cutoff=',DMIN
541   CONTINUE
   ENDDO
ENDIF ! end of block to add constraints and repulsions for frozen atoms.
CALL MYCPU_TIME(FTIME,.FALSE.)
WRITE(*,'(A,F10.1,A,I6)') ' intlbfgs> constrained potential finished, time=',FTIME-STIME,' number of repulsions=',NREPULSIVE
986 CONTINUE
STIME=FTIME
NSTEPSMAX=INTSTEPS1
!
! Don;t want to redistribute images before even taking a step, so don;t call CHECKSEP.
! Must call CHECKREP to initialise NNREULSIVE, NREPI, NREPJ, etc. SEGV otherwise on second cycle!
!
! To take BH-type steps in the QCI space, jump back here. Leave SWITCHED true.
!
BESTWORST=1.0D100
9876 CONTINUE
CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),0,1)


!                DIFF=1.0D-6
!                PRINT*,'analytic and numerical gradients:'
!                CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
!!               DO J1=2,INTIMAGE+1
!                DO J1=INTIMAGE-1,INTIMAGE+1
!                   DO J2=1,3*NATOMS
!                      IF (.NOT.ATOMACTIVE((J2-1)/3+1)) CYCLE
!                      J3=3*NATOMS*(J1-1)+J2
!                      XYZ(J3)=XYZ(J3)+DIFF
!                      CALL CONGRAD(NMAXINT,NMININT,EPLUS,XYZ,VPLUS,EEE,IMGFREEZE,RMS)
!                      XYZ(J3)=XYZ(J3)-2.0D0*DIFF
!                      CALL CONGRAD(NMAXINT,NMININT,EMINUS,XYZ,VMINUS,EEE,IMGFREEZE,RMS)
!                      XYZ(J3)=XYZ(J3)+DIFF
!                      IF ((ABS(GGG(J3)).NE.0.0D0).AND. &
!     &               (ABS(100.0D0*(GGG(J3)-(EPLUS-EMINUS)/(2.0D0*DIFF))/GGG(J3)).GT.1.0D0)) THEN
!                         WRITE(*,'(A,2I5,2G20.10,L5,A)') 'anal and num ',J1,J2,GGG(J3),(EPLUS-EMINUS)/(2.0D0*DIFF),  &
!     & ATOMACTIVE((J2-1)/3+1),'   X'   
!                      ELSE
!                         WRITE(*,'(A,2I5,2G20.10,L5)') 'anal and num ',J1,J2,GGG(J3),(EPLUS-EMINUS)/(2.0D0*DIFF),  &
!     & ATOMACTIVE((J2-1)/3+1)  
!                      ENDIF
!                   ENDDO
!                ENDDO
!  !
!                STOP


IF (QCIADDREP.GT.0) THEN
   CALL CONGRAD3(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
!
! Don't do the CONINT part of CONGRAD2 if CONINT isn't set. CONGRAD seems to be
! dong something different at the moment. Focus on CONGRAD2
!
ELSEIF (CHECKCONINT) THEN
   CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
ELSE
   CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
ENDIF
EOLD=ETOTAL
GLAST(1:D)=G(1:D)
XSAVE(1:D)=X(1:D)

IF (ETOTAL/INTIMAGE.LT.COLDFUSIONLIMIT) THEN
   WRITE(*,'(A,2G20.10)') ' intlbfgs> Cold fusion diagnosed - step discarded, energy, limit=', &
  &                       ETOTAL/INTIMAGE,COLDFUSIONLIMIT
   DEALLOCATE(CONI,CONJ,CONDISTREF,REPI,REPJ,NREPI,NREPJ,REPCUT,NREPCUT,CONCUT,CONOFFLIST,CONOFFTRIED)
   DEALLOCATE(TRUEEE, EEETMP, MYGTMP, GTMP, &
  &      DIAG, STP, SEARCHSTEP, GDIF,GLAST, XSAVE, XYZ, GGG, CHECKG, IMGFREEZE, EEE, STEPIMAGE)
   INTIMAGE=INTIMAGESAVE
   RETURN
ENDIF

! IF (DEBUG) WRITE(*,'(A6,A20,A20,A9,A9)') 'Iter','Energy per image','RMS Force','Step'

!
! In this block PERMGROUP(NDUMMY+J2-1) counts through the atom indices specified in perm.allow,
! one group at a time.
!
IF (PERMDIST) THEN
   PERMUTABLE(1:NATOMS)=.FALSE.
   INGROUP(1:NATOMS)=0
   NDUMMY=1
   DO J1=1,NPERMGROUP
      STARTGROUP(J1)=NDUMMY
      GROUPACTIVE(J1)=.FALSE.
      DO J2=1,NPERMSIZE(J1)
         PERMUTABLE(PERMGROUP(NDUMMY+J2-1))=.TRUE.
         INGROUP(PERMGROUP(NDUMMY+J2-1))=J1
      ENDDO
      NDUMMY=NDUMMY+NPERMSIZE(J1)
      ENDGROUP(J1)=NDUMMY-1
   ENDDO
ENDIF

567 CONTINUE

DO ! Main do loop with counter NITERDONE, initially set to one

!
! Are we stuck? 
!
IF (QCIRESET) THEN
!  IF ((SWITCHED.AND.(MOD(NITERDONE-1,QCIRESETINT2).EQ.0)).OR.((.NOT.SWITCHED).AND.(MOD(NITERDONE-1,QCIRESETINT1).EQ.0))) THEN
!  PRINT *,'intlbfgs> NITERDONE,NLASTGOODE,QCIRESETINT1=',NITERDONE,NLASTGOODE,QCIRESETINT1
   IF ((.NOT.SWITCHED).AND.(NITERDONE-NLASTGOODE.GT.QCIRESETINT1)) THEN ! .AND.(NITERDONE-NLASTCHANGE.GT.QCIRESETINT1)) THEN
      IF (DEBUG) CALL INTRWG(NACTIVE,NITERDONE,INTIMAGE,XYZ,TURNONORDER,NCONOFF)
      IF (DEBUG) CALL WRITEPROFILE(NITERDONE,EEE,INTIMAGE)
!     TURNOFF=.TRUE.
!     IF (EREP.GT.2.0D0*ECON) THEN
!        INTCONSTRAINTREP=INTCONSTRAINTREP/2.0D0
!        WRITE(*,'(A,G20.10)') 'intlbfgs> Interpolation seems to be stuck. Reducing repulsive prefactor to ',INTCONSTRAINTREP
!        TURNOFF=.FALSE.
!     ENDIF
!     IF (ESPRING.GT.2.0D0*ECON) THEN
!        KINT=KINT/2.0D0
!        WRITE(*,'(A,G20.10)') 'intlbfgs> Interpolation seems to be stuck. Reducing spring constant to ',KINT
!        TURNOFF=.FALSE.
!     ENDIF
!     IF (TURNOFF) THEN
!        WRITE(*,'(2(A,I6))') 'intlbfgs> Interpolation seems to be stuck. Dump images. Turn off worst constraint ',JMAXCON, &
! &                        ' total off=',NCONOFF+1
!        IF ((JMAXCON.LT.1).OR.(JMAXCON.GT.NCONSTRAINT)) THEN
!           WRITE(*,'(A)') 'intlbfgs> *** ERROR *** constraint index out of allowed range'
!           WRITE(*,'(A,I6)') 'NCONSTRAINT,NCONOFF=',NCONOFF
!           WRITE(*,'(A)') 'CONOFFTRIED:'
!           WRITE(*,'(20L5)') CONOFFTRIED(1:NCONSTRAINT)
!           STOP
!        ENDIF
!        NCONOFF=NCONOFF+1
!        CONOFFLIST(NCONOFF)=JMAXCON
!        CONACTIVE(JMAXCON)=.FALSE.
!        CONOFFTRIED(JMAXCON)=.TRUE.
!     ENDIF

      IF (MAX(CONVERGECONTEST,CONVERGEREPTEST).GT.MAXCONE) MAXCONE=MAXCONE*1.05D0
      IF (MAX(FCONTEST,FREPTEST).GT.INTRMSTOL) INTRMSTOL=INTRMSTOL*1.05D0
      WRITE(*,'(A,2G20.10)') 'intlbfgs> Interpolation seems to be stuck. Converge thresholds are now ',MAXCONE,INTRMSTOL

      NLASTGOODE=NITERDONE
      LASTGOODE=ETOTAL
   ENDIF
ENDIF

!
!  Check permutational alignments. Maintain a list of the permutable groups where all
!  members are active. See if we have any new complete groups. MUST update NDUMMY
!  counter to step through permutable atom list.
!
IF (QCILPERMDIST.AND.(MOD(NITERDONE-1,QCIPDINT).EQ.0)) THEN

   PRINT *,'DOING CHIRALCHECK NOW'
!       IF (DEBUG) WRITE(*,'(A)') 'intlbfgs> dump state before CHIRALCHECK index -4'
!        IF (DEBUG) CALL INTRWG2(NACTIVE,-4,INTIMAGE,XYZ,TURNONORDER,NCONOFF)

   chicheck: DO J5=1, num_chiral_centres
      atom_number=sr_atoms(J5, 1) 
!     WRITE(*,'(A,I6,A,I6,A,I6)') 'chiral centre ',J5,' is atom ',atom_number
      IF (.NOT.ATOMACTIVE(atom_number)) CYCLE chicheck
      DO J2=1,4
         IF (.NOT.ATOMACTIVE(sr_atoms(J5,J2))) CYCLE chicheck
      ENDDO

      DO J3=1,INTIMAGE+2

         CENTRE_COORDS(1)=XYZ(3*NATOMS*(J3-1)+3*(atom_number-1)+1)
         CENTRE_COORDS(2)=XYZ(3*NATOMS*(J3-1)+3*(atom_number-1)+2)
         CENTRE_COORDS(3)=XYZ(3*NATOMS*(J3-1)+3*(atom_number-1)+3)

         DO J4=1,4
            J2=sr_atoms(J5, J4 + 1)
            NEIGHBOUR_COORDS(3*(J4-1)+1)=XYZ(3*NATOMS*(J3-1)+3*(J2-1)+1)
            NEIGHBOUR_COORDS(3*(J4-1)+2)=XYZ(3*NATOMS*(J3-1)+3*(J2-1)+2)
            NEIGHBOUR_COORDS(3*(J4-1)+3)=XYZ(3*NATOMS*(J3-1)+3*(J2-1)+3)
         ENDDO

         CHIRALSR=CHIRALITY_SR(NEIGHBOUR_COORDS,CENTRE_COORDS)
!        WRITE(*,'(A,I6,I6,2L5)') 'image, atom, chirality, initial=',J3,atom_number,CHIRALSR,sr_states_initial(J5)
         IF (J3.EQ.1) CHIRALSRP=sr_states_initial(J5)
         IF (CHIRALSR.NEQV.CHIRALSRP) THEN
            WRITE(*,'(A,I6,A,I6,A,I6)') 'intlbfgs> Atom ',atom_number,' image ',J3,' chirality CHANGED; use previous image coordinates'  
            NLASTCHANGE=NITERDONE
!
! need to revert to whole aa coordinates, active atoms or not.
!
            ACID=ATOMSTORES(atom_number)

            DO J4=1,NATOMS
               IF (ATOMSTORES(J4).NE.ACID) CYCLE
               IF (.NOT.ATOMACTIVE(J4)) CYCLE
               WRITE(*,'(A,I6,A,I6,A,I6)') 'intlbfgs> Changing active atom ',J4,' image ',J3
               XYZ(3*NATOMS*(J3-1)+3*(J4-1)+1)=XYZ(3*NATOMS*(J3-2)+3*(J4-1)+1)
               XYZ(3*NATOMS*(J3-1)+3*(J4-1)+2)=XYZ(3*NATOMS*(J3-2)+3*(J4-1)+2)
               XYZ(3*NATOMS*(J3-1)+3*(J4-1)+3)=XYZ(3*NATOMS*(J3-2)+3*(J4-1)+3)
            ENDDO
         ENDIF
!        IF (J3.EQ.1) CHIRALSRP=CHIRALSR  ! just use result for fixed end point image 1
      ENDDO
   ENDDO chicheck
!        IF (DEBUG) WRITE(*,'(A)') 'intlbfgs> dump state after CHIRALCHECK index -3'
!        IF (DEBUG) CALL INTRWG2(NACTIVE,-3,INTIMAGE,XYZ,TURNONORDER,NCONOFF)

   NDUMMY=1
   DO J1=1,NPERMGROUP
      IF (GROUPACTIVE(J1)) GOTO 975
      DO J2=1,NPERMSIZE(J1)
         IF (.NOT.ATOMACTIVE(PERMGROUP(NDUMMY+J2-1))) GOTO 975
      ENDDO
      GROUPACTIVE(J1)=.TRUE.
      IF (DEBUG) WRITE(*,'(A,I6,A)') ' intlbfgs> All permutable atoms in group ',J1,' are active'
975   NDUMMY=NDUMMY+NPERMSIZE(J1)
   ENDDO 

!  IF (NITERDONE.EQ.1) THEN
!     COORDSB(1:3*NATOMS)=XYZ(1:3*NATOMS) ! starting endpoint
!     COORDSC(1:3*NATOMS)=XYZ(3*NATOMS*(INTIMAGE+1)+1:3*NATOMS*(INTIMAGE+2)) ! finish endpoint
!     WRITE(*,'(A)') ' intlbfgs> checking alignment of endpoints - should be no permutations'
!     CALL LOPERMDIST(COORDSB,COORDSC,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCE,DIST2,.FALSE.,RMATBEST,0,NMOVE,PERM)
!     WRITE(*,'(A,G20.10,A,I6)') ' intlbfgs> endpoint distance ',DISTANCE,' permutations=',NMOVE 
!  ENDIF

!       IF (DEBUG) WRITE(*,'(A)') 'intlbfgs> dump state before lopermdist index -6'
!        IF (DEBUG) CALL INTRWG2(NACTIVE,-6,INTIMAGE,XYZ,TURNONORDER,NCONOFF)
   SAVELOCALPERMCUT=LOCALPERMCUT
   LOCALPERMCUT=QCIPERMCUT ! 
   np: DO J1=1,NPERMGROUP
      IF (.NOT.GROUPACTIVE(J1)) CYCLE np
      DO J3=1,INTIMAGE
!        WRITE(*,'(A,I6,A,I6)') 'intlbfgs> Doing group ',J1,' image ',J3+1
!          COORDSB(1:3*NATOMS)=XYZ(3*NATOMS*(J3-1)+1:3*NATOMS*J3) ! coordinates for intervening image J3-1, which is the starting endpoint for J3=1
!          COORDSA(1:3*NATOMS)=XYZ(3*NATOMS*J3+1:3*NATOMS*(J3+1))  ! coordinates for intervening image J3, which is image J3+1 including starting endpoint
!          CALL QCIOPTPERM(COORDSB,COORDSA,NATOMS,DEBUG,J1,PERM,STARTGROUP,J3)
!          IF (.NOT.IDENTITY) THEN
!             DO J2=1,NATOMS
! !           IF (PERM(J2).NE.J2) WRITE(*,'(4(A,I6))') ' intlbfgs> image ',J3+1,' qcipermopt would move atom ',J2,' to position ',PERM(J2), &  
! !  &   ' group ',J1  
!                IF (PERM(J2).NE.J2) WRITE(*,'(4(A,I6))') ' intlbfgs> image ',J3+1,' qcipermopt move atom ',J2,' to position ',PERM(J2), &
!    &   ' group ',J1  
!                COORDSA(3*(J2-1)+1)=XYZ(3*NATOMS*J3+3*(PERM(J2)-1)+1)
!                COORDSA(3*(J2-1)+2)=XYZ(3*NATOMS*J3+3*(PERM(J2)-1)+2)
!                COORDSA(3*(J2-1)+3)=XYZ(3*NATOMS*J3+3*(PERM(J2)-1)+3)
!             ENDDO
!             XYZ(3*NATOMS*J3+1:3*NATOMS*(J3+1))=COORDSA(1:3*NATOMS)
!          ELSE
!             WRITE(*,'(A,I6,A,I6)') 'intlbfgs> identity permutation for image ',J3+1,' with preceeding image ',J3
!          ENDIF

         COORDSB(1:3*NATOMS)=XYZ(1:3*NATOMS) ! starting endpoint
         COORDSC(1:3*NATOMS)=XYZ(3*NATOMS*(INTIMAGE+1)+1:3*NATOMS*(INTIMAGE+2)) ! finish endpoint
         COORDSA(1:3*NATOMS)=XYZ(3*NATOMS*J3+1:3*NATOMS*(J3+1)) ! coordinates for intervening image J3, which is image J3+1 including starting endpoint
!        WRITE(*,'(A,I6,A,I6,A)') 'intlbfgs> Doing group ',J1,' image ',J3+1,' with start'
         CALL LOPERMDIST(COORDSB,COORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCES,DIST2,.FALSE.,RMATBEST,J1,NMOVES,PERMS)
         COORDSA(1:3*NATOMS)=XYZ(3*NATOMS*J3+1:3*NATOMS*(J3+1)) ! coordinates for intervening image J3, which is image J3+1 including starting endpoint
!        WRITE(*,'(A,I6,A,I6,A)') 'intlbfgs> Doing group ',J1,' image ',J3+1,' with finish'
         CALL LOPERMDIST(COORDSC,COORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCEF,DIST2,.FALSE.,RMATBEST,J1,NMOVEF,PERMF)
!        WRITE(*,'(A,I6,A,I6,A,2G20.10,A,2I6)') ' intlbfgs> image ',J3+1,' group ',J1,' start and finish distances=', &
!  &                                       DISTANCES,DISTANCEF,' permutations=',NMOVES, NMOVEF

         IF ((NMOVES.GT.0).OR.(NMOVEF.GT.0)) THEN
            WRITE(*,'(A,I6,A,I6,A,2G20.10,A,2I6)') ' intlbfgs> image ',J3+1,' group ',J1,' start and finish distances=', &
   &                                           DISTANCES,DISTANCEF,' permutations=',NMOVES, NMOVEF
         ENDIF

         IF ((NMOVES.GT.0).AND.(NMOVES.EQ.NMOVEF)) THEN
            COORDSA(1:3*NATOMS)=XYZ(3*NATOMS*J3+1:3*NATOMS*(J3+1)) 
            DO J2=1,NATOMS
               IF (PERMS(J2).NE.J2) THEN
                  IF (PERMF(J2).EQ.PERMS(J2)) THEN
                     WRITE(*,'(A,I6,A,I6,A,I6)') ' intlbfgs> image ',J3+1, &
   &                 ' consistent non-identity lopermdist permutation for start and finish, move atom ',PERMS(J2),' to position ',J2   
!                WRITE(*,'(A,I6,A,I6,A,I6)') ' intlbfgs> image ',J3+1, &
!  &           ' consistent non-identity lopermdist permutation for start and finish, would move atom ',PERMS(J2),' to position ',J2  
                     COORDSA(3*(J2-1)+1)=XYZ(3*NATOMS*J3+3*(PERMS(J2)-1)+1)
                     COORDSA(3*(J2-1)+2)=XYZ(3*NATOMS*J3+3*(PERMS(J2)-1)+2)
                     COORDSA(3*(J2-1)+3)=XYZ(3*NATOMS*J3+3*(PERMS(J2)-1)+3)
!                    NLASTCHANGE=NITERDONE
                  ELSE
                     WRITE(*,'(A,I6,A,I6)') ' intlbfgs> inconsistent non-identity lopermdist permutations for start and finish'
                  ENDIF
               ENDIF
            ENDDO
            XYZ(3*NATOMS*J3+1:3*NATOMS*(J3+1))=COORDSA(1:3*NATOMS)
         ENDIF
       ENDDO
    ENDDO np
    LOCALPERMCUT=SAVELOCALPERMCUT

!    COORDSB(1:3*NATOMS)=XYZ(1:3*NATOMS) ! starting endpoint
!    COORDSC(1:3*NATOMS)=XYZ(3*NATOMS*(INTIMAGE+1)+1:3*NATOMS*(INTIMAGE+2)) ! finish endpoint
!    SAVELOCALPERMCUT=LOCALPERMCUT
! !
! ! needs to be a bit sloppier than usual to allow for distorted geometries. Otherwise we may get wrong assignments
! ! because we don't have enough neighbours.
! !
!    LOCALPERMCUT=QCIPERMCUT ! 0.8D0 
!    DO J3=1,INTIMAGE
!       np: DO J1=1,NPERMGROUP
!          IF (.NOT.GROUPACTIVE(J1)) CYCLE np
!          COORDSA(1:3*NATOMS)=XYZ(3*NATOMS*J3+1:3*NATOMS*(J3+1)) ! coordinates for intervening image J3, which is image J3+1 including starting endpoint
!          CALL LOPERMDIST(COORDSB,COORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCES,DIST2,.FALSE.,RMATBEST,J1,NMOVES,PERMS)
!          COORDSA(1:3*NATOMS)=XYZ(3*NATOMS*J3+1:3*NATOMS*(J3+1)) ! coordinates for intervening image J3, which is image J3+1 including starting endpoint
!          CALL LOPERMDIST(COORDSC,COORDSA,NATOMS,DEBUG,BOXLX,BOXLY,BOXLZ,BULKT,TWOD,DISTANCEF,DIST2,.FALSE.,RMATBEST,J1,NMOVEF,PERMF)
!             WRITE(*,'(A,I6,A,I6,A,2G20.10,A,2I6)') ' intlbfgs> image ',J3+1,' group ',J1,' start and finish distances=', &
!   &                                       DISTANCES,DISTANCEF,' permutations=',NMOVES, NMOVEF
!          IF ((NMOVES.GT.0).AND.(NMOVEF.GT.0)) THEN
!             WRITE(*,'(A,I6,A,I6,A,2G20.10,A,2I6)') ' intlbfgs> image ',J3+1,' group ',J1,' start and finish distances=', &
!   &                                       DISTANCES,DISTANCEF,' permutations=',NMOVES, NMOVEF
!          ENDIF
!          IF ((NMOVES.GT.0).AND.(NMOVES.EQ.NMOVEF)) THEN
!             COORDSA(1:3*NATOMS)=XYZ(3*NATOMS*J3+1:3*NATOMS*(J3+1)) 
!             DO J2=1,NATOMS
!                IF (PERMS(J2).NE.J2) THEN
!                   IF (PERMF(J2).EQ.PERMS(J2)) THEN
!                      WRITE(*,'(A,I6,A,I6,A,I6)') ' intlbfgs> image ',J3+1, &
!   &                                           ' consistent non-identity permutation for start and finish, move atom ', &  
!   &                                           PERMS(J2),' to position ',J2  
!                      COORDSA(3*(J2-1)+1)=XYZ(3*NATOMS*J3+3*(PERMS(J2)-1)+1)
!                      COORDSA(3*(J2-1)+2)=XYZ(3*NATOMS*J3+3*(PERMS(J2)-1)+2)
!                      COORDSA(3*(J2-1)+3)=XYZ(3*NATOMS*J3+3*(PERMS(J2)-1)+3)
!                   ELSE
!                      WRITE(*,'(A,I6,A,I6)') ' intlbfgs> inconsistent non-identity permutations for start and finish'
!                   ENDIF
!                ENDIF
!             ENDDO
!             XYZ(3*NATOMS*J3+1:3*NATOMS*(J3+1))=COORDSA(1:3*NATOMS)
!          ENDIF
!       ENDDO np
!    ENDDO
!    LOCALPERMCUT=SAVELOCALPERMCUT
!    CALL INTRWG(NACTIVE,NITERDONE,INTIMAGE,XYZ,TURNONORDER,NCONOFF)

!       IF (DEBUG) WRITE(*,'(A)') 'intlbfgs> dump state after lopermdist index -7'
!        IF (DEBUG) CALL INTRWG2(NACTIVE,-7,INTIMAGE,XYZ,TURNONORDER,NCONOFF)

ENDIF

!
!  Dynamic adjustment of KINT values. Local values are changed by QCIKADJUSTFRAC if the
!  corresponding separation is outside a fraction QCIKADJUSTTOL of the average value.
!  The adjustment is done every QCIKADJUSTFRQ cycles.
!  Based on the DNEB adjustment.
!
IF (QCIKADJUSTFRQ.GT.0) THEN
   IF (MOD(NITERDONE,QCIKADJUSTFRQ).EQ.0) THEN ! dynamic adjustment of KINT
      IF (QCIAVDEV.GT.QCIKADJUSTTOL) THEN
         KINT=MIN(KINT*QCIKADJUSTFRAC,QCIKINTMAX)
         IF (DEBUG) PRINT '(2(A,G20.10))',' intlbfgs> Mean deviation ',QCIAVDEV,' Increasing QCI force constant to ',KINT
      ELSEIF (QCIAVDEV.LT.QCIKADJUSTTOL) THEN
         KINT=MAX(KINT/QCIKADJUSTFRAC,QCIKINTMIN)
         IF (DEBUG) PRINT '(2(A,G20.10))',' intlbfgs> Mean deviation ',QCIAVDEV,' Decreasing QCI force constant to ',KINT  
      ENDIF
   ENDIF
ENDIF


!
!  Add next atom to active set if ADDATOM is true. 
!  Constraints to atoms already in the active set are turned on
!  and short-range repulsions to active atoms that are not distance constrained are turned on.
!  *** OLD Find nearest atom to active set attached by a constraint
!  *** NEW Find atom with most constraints to active set
!  Turn on constraint terms for this atom with all previous members of the active set
!  Add repulsions to non-constrained atoms in this set
!  NTOADD is the number of atoms to add to the active set in each pass. 1 seems best!
!
   IF (ADDATOM.AND.((NACTIVE.LT.NATOMS).OR.(NCONOFF.GT.0))) THEN

!!!!!!!!!!!!!!!DEBUG DJW !!!!!!!!!!!
!!
!!               J2=0
!!               DO J1=1,NREPULSIVEFIX
!!!                 WRITE(*,'(A,3I10,4L5)') 'doaddatom> J1,REPIFIX,REPJFIX,frozenI,frozenJ,activeI,activeJ=', &
!!! &                 J1,REPIFIX(J1),REPJFIX(J1),INTFROZEN(REPIFIX(J1)),INTFROZEN(REPJFIX(J1)), &
!!! &                 ATOMACTIVE(REPIFIX(J1)),ATOMACTIVE(REPJFIX(J1))
!!                  IF (INTFROZEN(REPIFIX(J1)).AND.INTFROZEN(REPJFIX(J1))) CYCLE
!!                  IF (ATOMACTIVE(REPIFIX(J1)).AND.ATOMACTIVE(REPJFIX(J1))) THEN
!!                     DO J3=1,NCONSTRAINTFIX
!!!                       IF (.NOT.CONACTIVE(J3)) CYCLE ! repulsions for inactive constraints
!!                        IF ((CONIFIX(J3).EQ.REPIFIX(J1)).AND.(CONJFIX(J3).EQ.REPJFIX(J1))) GOTO 963
!!                        IF ((CONIFIX(J3).EQ.REPJFIX(J1)).AND.(CONJFIX(J3).EQ.REPIFIX(J1))) GOTO 963
!!                     ENDDO
!!                     J2=J2+1
!!!                    WRITE(*,'(A,I10,A,2I6)') 'doaddatom> repulsion ',J2,' between ',REPIFIX(J1),REPJFIX(J1)
!!963                  CONTINUE
!!                  ENDIF
!!               ENDDO
!!               WRITE(*,'(A,I6,A)') 'doaddatom> Looks like there are ',J2,' possible repulsions before adding new atom'
!!
!!               NDUMMY=1
!!               NREPULSIVE=0
!!               DO J1=1,NATOMS
!!                  IF (.NOT.ATOMACTIVE(J1)) CYCLE
!!!
!!! Make a list of repelling atoms here and then use it
!!! CONI(J2) is always less than CONJ(J2) so we only need to
!!! cycle over a given range of constraints and continue from
!!! where we left off for the next atom j1
!!!
!!                  ADDREP(1:J1+INTREPSEP)=.FALSE.
!!                  ADDREP(J1+INTREPSEP+1:NATOMS)=.TRUE. ! no repulsion for atoms too close in sequence
!!                  IF (INTFROZEN(J1)) THEN
!!                     DO J2=J1+INTREPSEP+1,NATOMS
!!                        IF (INTFROZEN(J2)) ADDREP(J2)=.FALSE.
!!                        IF (.NOT.ATOMACTIVE(J2)) ADDREP(J2)=.FALSE.
!!                     ENDDO
!!                  ENDIF
!!                  myaddloop: DO J2=NDUMMY,NCONSTRAINTFIX
!!!                    IF (.NOT.CONACTIVE(J2)) CYCLE myaddloop ! repulsions for inactive constraints
!!                     IF (CONIFIX(J2).EQ.J1) THEN
!!                        ADDREP(CONJFIX(J2))=.FALSE.
!!!
!!! The next line is different from make_conpot because we don't count the constraints
!!! sequentially, due to the ATOMACTIVE(J1) test at the top.
!!!
!!                     ELSEIF (CONIFIX(J2).GT.J1) THEN
!!                        NDUMMY=J2 ! for next atom
!!                        EXIT myaddloop
!!                     ENDIF
!!                  ENDDO myaddloop
!!                  myrep2: DO J2=J1+INTREPSEP+1,NATOMS
!!                     IF (.NOT.ADDREP(J2)) CYCLE myrep2
!!                     IF (.NOT.ATOMACTIVE(J2)) CYCLE myrep2 ! This line is not in make_conpot, where we want all possible repulsions.
!!                     DMIN=1.0D100
!!                     DO J3=1,INTIMAGE+2,INTIMAGE+1 ! only consider the end-point distances
!!                        DF=SQRT((XYZ((J3-1)*3*NATOMS+3*(J2-1)+1)-XYZ((J3-1)*3*NATOMS+3*(J1-1)+1))**2+ &
!!    &                     (XYZ((J3-1)*3*NATOMS+3*(J2-1)+2)-XYZ((J3-1)*3*NATOMS+3*(J1-1)+2))**2+ &
!!    &                     (XYZ((J3-1)*3*NATOMS+3*(J2-1)+3)-XYZ((J3-1)*3*NATOMS+3*(J1-1)+3))**2)
!!                        IF (DF.LT.DMIN) DMIN=DF
!!                     ENDDO
!!
!!                     NREPULSIVE=NREPULSIVE+1
!!                     REPI(NREPULSIVE)=J1
!!                     REPJ(NREPULSIVE)=J2
!!!                    WRITE(*,'(A,I10,A,2I6)') 'doaddatom> repulsion ',NREPULSIVE,' between ',J1,J2
!!!
!!! Use the minimum of the end point distances and INTCONSTRAINREPCUT for each contact.
!!!
!!                     REPCUT(NREPULSIVE)=MIN(DMIN-1.0D-3,INTCONSTRAINREPCUT)
!!                  ENDDO myrep2
!!               ENDDO
!!               WRITE(*,'(A,I6,A)') ' intlbfgs> Now it looks like there are ',NREPULSIVE,' possible repulsions before adding new atom'
!!!!!!!!!!!!!!!DEBUG DJW !!!!!!!!!!!

      IF (NCONOFF.GT.0) THEN
         CONACTIVE(CONOFFLIST(NCONOFF))=.TRUE.
         WRITE(*,'(2(A,I6))') 'intlbfgs> Turn back on constraint ',CONOFFLIST(NCONOFF),' total off=',NCONOFF-1
         NCONOFF=NCONOFF-1
      ELSE
!        IF (DEBUG) WRITE(*,'(A)') 'intlbfgs> dump state before doaddatom index -2'
!        IF (DEBUG) CALL INTRWG2(NACTIVE,-2,INTIMAGE,XYZ,TURNONORDER,NCONOFF)
         CALL DOADDATOM(NCONSTRAINT,NTRIES,NEWATOM,IMGFREEZE,INTIMAGE,XYZ,EEE,GGG,TURNONORDER,NITERDONE,NACTIVE,AABACK,BACKDONE)  
!        IF (DEBUG) WRITE(*,'(A)') 'intlbfgs> dump state after doaddatom index -1'
!        IF (DEBUG) CALL INTRWG2(NACTIVE,-1,INTIMAGE,XYZ,TURNONORDER,NCONOFF)
      ENDIF
      NLASTGOODE=NITERDONE
      LASTGOODE=ETOTAL
   ENDIF
   GTMP(1:D)=0.0D0
   CALL MAKESTEP(NITERUSE,POINT,DIAG,INTIMAGE,SEARCHSTEP,G,GTMP,STP,GDIF,NPT,D,RHO1,ALPHA)
!
! If the number of images has changed since G was declared then G is not the same
! size as Gtmp and Dot_Product cannot be used.
!
!  IF (Dot_Product(G,Gtmp)/SQRT( Dot_Product(G,G)*Dot_Product(Gtmp,Gtmp) ) > 0.0D0) THEN
!
!  Separate sqrt;s to avoid overflow.
!
   IF (DDOT(D,G,1,GTMP,1)/MAX(1.0D-100,SQRT( DDOT(D,G,1,G,1))*SQRT(DDOT(D,GTMP,1,GTMP,1)) ) > 0.0D0) THEN
        IF (DEBUG) WRITE(*,*) 'Search direction has positive projection onto gradient - reversing step'
        GTMP(1:D)=-GTMP(1:D)
        SEARCHSTEP(POINT,1:D)=GTMP(1:D)
   ENDIF
   GTMP(1:D)=G(1:D)

!  We should apply the maximum LBFGS step to each image separately.
!  However, using different scale factors for different images leads to huge
!  discontinuities! Now take the minimum scale factor for all images. DJW 26/11/07

   STPMIN=1.0D0
   DO J2=1,INTIMAGE
      STEPIMAGE(J2) = SQRT(DOT_PRODUCT(SEARCHSTEP(POINT,(3*NATOMS)*(J2-1)+1:(3*NATOMS)*J2), &
  &                                    SEARCHSTEP(POINT,(3*NATOMS)*(J2-1)+1:(3*NATOMS)*J2)))
      DUMMY=STEPIMAGE(J2)
      IF (STEPIMAGE(J2) > MAXINTBFGS) THEN
           STP((3*NATOMS)*(J2-1)+1:(3*NATOMS)*J2) = MAXINTBFGS/STEPIMAGE(J2)
           STPMIN=MIN(STPMIN,STP((3*NATOMS)*(J2-1)+1))
      ENDIF
!     WRITE(*,'(A,I8,3G20.10)') ' image,initial step size,STP,prod=',J2,DUMMY,STP(3*NATOMS*(J2-1)+1), &
! &                                   STEPIMAGE(J2)*STP(3*NATOMS*(J2-1)+1)   
   ENDDO
   STP(1:D)=STPMIN
! EFK: decide whether to freeze some nodes
   IF (FREEZENODEST) THEN
      TOTGNORM=SQRT(DOT_PRODUCT(G(1:(3*NATOMS)*INTIMAGE),G(1:(3*NATOMS)*INTIMAGE))/INTIMAGE)
      NIMAGEFREEZE=0
      DO IM=1,INTIMAGE
         TESTG=SQRT(DOT_PRODUCT(G((3*NATOMS)*(IM-1)+1:(3*NATOMS)*IM),G((3*NATOMS)*(IM-1)+1:(3*NATOMS)*IM)))
         IMGFREEZE(IM)=.FALSE.
         IF (TOTGNORM.NE.0.0D0) THEN
!           IF (TESTG/TOTGNORM.LT.FREEZETOL) THEN
            IF (TESTG/SQRT(3.0D0*NATOMS).LT.FREEZETOL) THEN
!              IF (DEBUG) PRINT '(A,I6,3G20.10)', ' intlbfgs> Freezing image: ',IM,TESTG,FREEZETOL,TOTGNORM
               IMGFREEZE(IM)=.TRUE.
               STEPIMAGE(IM)=0.0D0
               NIMAGEFREEZE=NIMAGEFREEZE+1
               STP((3*NATOMS)*(IM-1)+1:(3*NATOMS)*IM)=0.0D0
            ENDIF
         ENDIF
      ENDDO
      IF (DEBUG) PRINT '(2(A,I6))', ' intlbfgs> Number of frozen images=',NIMAGEFREEZE,' / ',INTIMAGE
   ENDIF
   !  We now have the proposed step - update geometry and calculate new gradient
   NDECREASE=0
20 X(1:D) = X(1:D) + STP(1:D)*SEARCHSTEP(POINT,1:D)

!  IF (.NOT.SWITCHED) THEN
   IF (.TRUE.) THEN
!     IF ((RMS.LT.INTRMSTOL*1.0D10).AND.(MOD(NITERDONE,10).EQ.0).AND.(NSTEPSMAX-NITERDONE.GT.100)) &
! &               CALL CHECKSEP(NMAXINT,NMININT,INTIMAGE,XYZ,(3*NATOMS),NATOMS)
!     PRINT '(A,3I10)','NITERDONE,INTIMAGECHECK,mod=',NITERDONE,INTIMAGECHECK,MOD(NITERDONE,INTIMAGECHECK)
      IF (REMOVEIMAGE.OR.(MOD(NITERDONE,INTIMAGECHECK).EQ.0)) THEN
864      CONTINUE ! for adding more than one image at a time
         DMAX=-1.0D0
         ADMAX=-1.0D0
         DMIN=HUGE(1.0D0)
         DO J1=1,INTIMAGE+1
            DUMMY=0.0D0
!           DO J2=1,3*NATOMS
!              IF (ATOMACTIVE((J2-1)/3+1)) THEN
!                 DUMMY=DUMMY+( XYZ((J1-1)*3*NATOMS+J2) - XYZ(J1*3*NATOMS+J2) )**2
!              ENDIF
!           ENDDO
            DO J2=1,NATOMS
               IF (ATOMACTIVE(J2)) THEN
                  ADUMMY=( XYZ((J1-1)*3*NATOMS+3*(J2-1)+1) - XYZ(J1*3*NATOMS+3*(J2-1)+1) )**2 &
  &                     +( XYZ((J1-1)*3*NATOMS+3*(J2-1)+2) - XYZ(J1*3*NATOMS+3*(J2-1)+2) )**2 &
  &                     +( XYZ((J1-1)*3*NATOMS+3*(J2-1)+3) - XYZ(J1*3*NATOMS+3*(J2-1)+3) )**2 
                  DUMMY=DUMMY+ADUMMY
                  IF (ADUMMY.GT.ADMAX) THEN
                     ADMAX=ADUMMY
                     JA1=J1
                     JA2=J2
                  ENDIF
               ENDIF
            ENDDO
            DUMMY=SQRT(DUMMY)
            IF (DUMMY.GT.DMAX) THEN
               DMAX=DUMMY
               JMAX=J1
            ENDIF
            IF (DUMMY.LT.DMIN) THEN
               DMIN=DUMMY
               JMIN=J1
            ENDIF
!            IF (DEBUG) WRITE(*,'(A,I6,A,I6,A,G20.10)')' intlbfgs> distance between images ', &
!  &                                                  J1,' and ',J1+1,' is ',DUMMY
!            IF (DEBUG) WRITE(*,'(A,G20.10,A,I6,A,2I6)')' intlbfgs> largest atomic distance between images so far is ', &
!  &                                                  SQRT(ADMAX),' for atom ',JA2,' and images ',JA1,JA1+1
         ENDDO
!        IF (DEBUG) WRITE(*,'(A,G20.10,A,I6,A,2I6,A,I6)')' intlbfgs> largest atomic distance between images is ', &
! &                                                  SQRT(ADMAX),' for atom ',JA2,' and images ',JA1,JA1+1,' total images=',INTIMAGE
!        IF (DEBUG) WRITE(*,'(A,G20.10,A,2I6)')' intlbfgs> largest image separation is ', &
! &                                                  DMAX,' for images ',JMAX,JMAX+1
!        IF (DEBUG) WRITE(*,'(A,G20.10,A,2I6)')' intlbfgs> smallest image separation is ', &
! &                                                  DMIN,' for images ',JMIN,JMIN+1
!        IF (DEBUG) WRITE(*,'(A,G20.10,A,G20.10)') ' intlbfgs> Mean image separation=',DUMMY2/(INTIMAGE+1),' per active atom=',DUMMY2/((INTIMAGE+1)*NACTIVE)
!        IF ((DMAX.GT.IMSEPMAX).AND.(INTIMAGE.LT.MAXINTIMAGE)) THEN
!        PRINT '(A,2G20.10)','SQRT(ADMAX),IMSEPMAX=',SQRT(ADMAX),IMSEPMAX
         IF ((.NOT.REMOVEIMAGE).AND.((SQRT(ADMAX).GT.IMSEPMAX).AND.(INTIMAGE.LT.MAXINTIMAGE))) THEN
            JMAX=JA1
            WRITE(*,'(A,I6,A,I6,A,I6)') ' intlbfgs> Add an image between ',JMAX,' and ',JMAX+1,' INTIMAGE=',INTIMAGE
            NITERUSE=0
            ALLOCATE(DPTMP(3*NATOMS*(INTIMAGE+2)))
            DPTMP(1:3*NATOMS*(INTIMAGE+2))=XYZ(1:3*NATOMS*(INTIMAGE+2))
            DEALLOCATE(XYZ)
            ALLOCATE(XYZ(3*NATOMS*(INTIMAGE+3)))
            XYZ(1:3*NATOMS*JMAX)=DPTMP(1:3*NATOMS*JMAX)
            XYZ(3*NATOMS*JMAX+1:3*NATOMS*(JMAX+1))=(DPTMP(3*NATOMS*(JMAX-1)+1:3*NATOMS*JMAX) &
  &                                               + DPTMP(3*NATOMS*JMAX+1:3*NATOMS*(JMAX+1)))/2.0D0
            XYZ(3*NATOMS*(JMAX+1)+1:3*NATOMS*(INTIMAGE+3))=DPTMP(3*NATOMS*JMAX+1:3*NATOMS*(INTIMAGE+2))
!
! Save step-taking memories in SEARCHSTEP and GDIF.
! These arrays run from 0 to INTMUPDATE over memories and
! 1:(3*NATOMS)*INTIMAGE over only the variable images.
!
            DEALLOCATE(DPTMP)
            ALLOCATE(D2TMP(0:INTMUPDATE,1:(3*NATOMS)*INTIMAGE))
            D2TMP(0:INTMUPDATE,1:(3*NATOMS)*INTIMAGE)=SEARCHSTEP(0:INTMUPDATE,1:(3*NATOMS)*INTIMAGE)
            DEALLOCATE(SEARCHSTEP)
            ALLOCATE(SEARCHSTEP(0:INTMUPDATE,1:(3*NATOMS)*(INTIMAGE+1)))
            DO J1=0,INTMUPDATE
               IF (JMAX.GT.1) SEARCHSTEP(J1,1:3*NATOMS*(JMAX-1))=D2TMP(J1,1:3*NATOMS*(JMAX-1))
               IF (JMAX.LT.INTIMAGE+1) SEARCHSTEP(J1,3*NATOMS*JMAX+1:3*NATOMS*(INTIMAGE+1))= &
  &                 D2TMP(J1,3*NATOMS*(JMAX-1)+1:3*NATOMS*INTIMAGE)
               SEARCHSTEP(J1,3*NATOMS*(JMAX-1)+1:3*NATOMS*JMAX)= &
  &                             D2TMP(J1,3*NATOMS*(MIN(JMAX,INTIMAGE)-1)+1:3*NATOMS*MIN(JMAX,INTIMAGE))
            ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            SEARCHSTEP(0:INTMUPDATE,1:(3*NATOMS)*(INTIMAGE+1))=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            D2TMP(0:INTMUPDATE,1:(3*NATOMS)*INTIMAGE)=GDIF(0:INTMUPDATE,1:(3*NATOMS)*INTIMAGE)
            DEALLOCATE(GDIF)
            ALLOCATE(GDIF(0:INTMUPDATE,1:(3*NATOMS)*(INTIMAGE+1)))
            DO J1=0,INTMUPDATE
               IF (JMAX.GT.1) GDIF(J1,1:3*NATOMS*(JMAX-1))=D2TMP(J1,1:3*NATOMS*(JMAX-1))
               IF (JMAX.LT.INTIMAGE+1) GDIF(J1,3*NATOMS*JMAX+1:3*NATOMS*(INTIMAGE+1))= &
  &                 D2TMP(J1,3*NATOMS*(JMAX-1)+1:3*NATOMS*INTIMAGE)
               GDIF(J1,3*NATOMS*(JMAX-1)+1:3*NATOMS*JMAX)= &
  &                       D2TMP(J1,3*NATOMS*(MIN(JMAX,INTIMAGE)-1)+1:3*NATOMS*MIN(JMAX,INTIMAGE))
            ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            GDIF(0:INTMUPDATE,1:(3*NATOMS)*(INTIMAGE+1))=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DEALLOCATE(D2TMP)

            DEALLOCATE(TRUEEE,EEETMP,MYGTMP,GTMP,GGG, &
  &                    DIAG,STP,GLAST,XSAVE,EEE,STEPIMAGE,CHECKG,IMGFREEZE)
            ALLOCATE(TRUEEE(INTIMAGE+3), &
  &                  EEETMP(INTIMAGE+3), MYGTMP(3*NATOMS*(INTIMAGE+1)), &
  &                  GTMP(3*NATOMS*(INTIMAGE+1)), &
  &                  DIAG(3*NATOMS*(INTIMAGE+1)), STP(3*NATOMS*(INTIMAGE+1)), &
  &                  GLAST((3*NATOMS)*(INTIMAGE+1)), &
  &                  XSAVE((3*NATOMS)*(INTIMAGE+1)), CHECKG((3*NATOMS)*(INTIMAGE+1)), IMGFREEZE(INTIMAGE+1), &
  &                  EEE(INTIMAGE+3), STEPIMAGE(INTIMAGE+1), GGG(3*NATOMS*(INTIMAGE+3)))
            GGG(1:3*NATOMS*(INTIMAGE+3))=0.0D0
            TRUEEE(1:INTIMAGE+3)=0.0D0
            EEETMP(1:INTIMAGE+3)=0.0D0
            MYGTMP(1:3*NATOMS*(INTIMAGE+1))=0.0D0
            GTMP(1:3*NATOMS*(INTIMAGE+1))=0.0D0
            DIAG(1:3*NATOMS*(INTIMAGE+1))=0.0D0
            STP(1:3*NATOMS*(INTIMAGE+1))=0.0D0
            GLAST(1:(3*NATOMS)*(INTIMAGE+1))=0.0D0
            XSAVE(1:(3*NATOMS)*(INTIMAGE+1))=0.0D0
            CHECKG(1:(3*NATOMS)*(INTIMAGE+1))=.FALSE.
            IMGFREEZE(1:INTIMAGE+1)=.FALSE.
            EEE(1:INTIMAGE+3)=0.0D0
            STEPIMAGE(1:INTIMAGE+1)=0.0D0

            X=>XYZ((3*NATOMS)+1:(3*NATOMS)*(INTIMAGE+2))
            G=>GGG((3*NATOMS)+1:(3*NATOMS)*(INTIMAGE+2))
            INTIMAGE=INTIMAGE+1
            D=(3*NATOMS)*INTIMAGE
            CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),0,1)
            IF (QCIADDREP.GT.0) THEN
               CALL CONGRAD3(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSEIF (CHECKCONINT) THEN
               CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSE
               CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ENDIF
!           GOTO 864
         ENDIF
         IF (REMOVEIMAGE.OR.((DMIN.LT.IMSEPMIN).AND.(INTIMAGE.GT.1))) THEN
            IF (REMOVEIMAGE) JMIN=JMAXEEE
            IF (JMIN.EQ.1) JMIN=2
            WRITE(*,'(A,I6,A,I6)') ' intlbfgs> Remove image ',JMIN
            NITERUSE=0
            ALLOCATE(DPTMP(3*NATOMS*(INTIMAGE+2)))
            DPTMP(1:3*NATOMS*(INTIMAGE+2))=XYZ(1:3*NATOMS*(INTIMAGE+2))
            DEALLOCATE(XYZ)
            ALLOCATE(XYZ(3*NATOMS*(INTIMAGE+1)))
            XYZ(1:3*NATOMS*(JMIN-1))=DPTMP(1:3*NATOMS*(JMIN-1))
            XYZ(3*NATOMS*(JMIN-1)+1:3*NATOMS*(INTIMAGE+1))=DPTMP(3*NATOMS*JMIN+1:3*NATOMS*(INTIMAGE+2))

            DEALLOCATE(DPTMP)
!
! Save step-taking memories in SEARCHSTEP and GDIF.
! These arrays run from 0 to INTMUPDATE over memories and
! 1:(3*NATOMS)*INTIMAGE over only the variable images.
!
            ALLOCATE(D2TMP(0:INTMUPDATE,1:(3*NATOMS)*INTIMAGE))
            D2TMP(0:INTMUPDATE,1:(3*NATOMS)*INTIMAGE)=SEARCHSTEP(0:INTMUPDATE,1:(3*NATOMS)*INTIMAGE)
            DEALLOCATE(SEARCHSTEP)
            ALLOCATE(SEARCHSTEP(0:INTMUPDATE,1:(3*NATOMS)*(INTIMAGE-1)))
            DO J1=0,INTMUPDATE
               SEARCHSTEP(J1,1:3*NATOMS*(JMIN-2))=D2TMP(J1,1:3*NATOMS*(JMIN-2))
               SEARCHSTEP(J1,3*NATOMS*(JMIN-2)+1:3*NATOMS*(INTIMAGE-1))= &
  &                     D2TMP(J1,3*NATOMS*(JMIN-1)+1:3*NATOMS*INTIMAGE)
            ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            SEARCHSTEP(0:INTMUPDATE,1:(3*NATOMS)*(INTIMAGE-1))=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            D2TMP(0:INTMUPDATE,1:(3*NATOMS)*INTIMAGE)=GDIF(0:INTMUPDATE,1:(3*NATOMS)*INTIMAGE)
            DEALLOCATE(GDIF)
            ALLOCATE(GDIF(0:INTMUPDATE,1:(3*NATOMS)*(INTIMAGE-1)))
            DO J1=0,INTMUPDATE
               GDIF(J1,1:3*NATOMS*(JMIN-2))=D2TMP(J1,1:3*NATOMS*(JMIN-2))
               GDIF(J1,3*NATOMS*(JMIN-2)+1:3*NATOMS*(INTIMAGE-1))= &
  &                     D2TMP(J1,3*NATOMS*(JMIN-1)+1:3*NATOMS*INTIMAGE)
            ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            GDIF(0:INTMUPDATE,1:(3*NATOMS)*(INTIMAGE-1))=0.0D0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            DEALLOCATE(D2TMP)

            DEALLOCATE(TRUEEE,EEETMP,MYGTMP,GTMP,GGG, &
  &                    DIAG,STP,GLAST,XSAVE,EEE,STEPIMAGE,CHECKG,IMGFREEZE)
            ALLOCATE(TRUEEE(INTIMAGE+1),&
  &                  EEETMP(INTIMAGE+1), MYGTMP(3*NATOMS*(INTIMAGE-1)), &
  &                  GTMP(3*NATOMS*(INTIMAGE-1)), &
  &                  DIAG(3*NATOMS*(INTIMAGE-1)), STP(3*NATOMS*(INTIMAGE-1)), &
  &                  GLAST((3*NATOMS)*(INTIMAGE-1)), &
  &                  XSAVE((3*NATOMS)*(INTIMAGE-1)), CHECKG((3*NATOMS)*(INTIMAGE-1)), IMGFREEZE(INTIMAGE-1), &
  &                  EEE(INTIMAGE+1), STEPIMAGE(INTIMAGE-1), GGG(3*NATOMS*(INTIMAGE+1)))
            GGG(1:3*NATOMS*(INTIMAGE+1))=0.0D0
            TRUEEE(1:INTIMAGE+1)=0.0D0
            EEETMP(1:INTIMAGE+1)=0.0D0
            MYGTMP(1:3*NATOMS*(INTIMAGE-1))=0.0D0
            GTMP(1:3*NATOMS*(INTIMAGE-1))=0.0D0
            DIAG(1:3*NATOMS*(INTIMAGE-1))=0.0D0
            STP(1:3*NATOMS*(INTIMAGE-1))=0.0D0
            GLAST(1:(3*NATOMS)*(INTIMAGE-1))=0.0D0
            XSAVE(1:(3*NATOMS)*(INTIMAGE-1))=0.0D0
            CHECKG(1:(3*NATOMS)*(INTIMAGE-1))=.FALSE.
            IMGFREEZE(1:INTIMAGE-1)=.FALSE.
            EEE(1:INTIMAGE+1)=0.0D0
            STEPIMAGE(1:INTIMAGE-1)=0.0D0

            X=>XYZ((3*NATOMS)+1:(3*NATOMS)*(INTIMAGE))
            G=>GGG((3*NATOMS)+1:(3*NATOMS)*(INTIMAGE))
            INTIMAGE=INTIMAGE-1
            D=(3*NATOMS)*INTIMAGE
            CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),0,1)
            IF (QCIADDREP.GT.0) THEN
               CALL CONGRAD3(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSEIF (CHECKCONINT) THEN
               CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSE
               CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ENDIF
            NLASTGOODE=NITERDONE
            LASTGOODE=ETOTAL
!           GOTO 864
         ENDIF
      ELSE
         DMAX=-1.0D0
         ADMAX=-1.0D0
         DMIN=HUGE(1.0D0)
         DUMMY2=0.0D0
         DO J1=1,INTIMAGE+1
            DUMMY=0.0D0
!           DO J2=1,3*NATOMS
!              IF (ATOMACTIVE((J2-1)/3+1)) THEN
!                 DUMMY=DUMMY+( XYZ((J1-1)*3*NATOMS+J2) - XYZ(J1*3*NATOMS+J2) )**2
!              ENDIF
!           ENDDO
            DO J2=1,NATOMS
               IF (ATOMACTIVE(J2)) THEN
                  ADUMMY=( XYZ((J1-1)*3*NATOMS+3*(J2-1)+1) - XYZ(J1*3*NATOMS+3*(J2-1)+1) )**2 &
  &                     +( XYZ((J1-1)*3*NATOMS+3*(J2-1)+2) - XYZ(J1*3*NATOMS+3*(J2-1)+2) )**2 &
  &                     +( XYZ((J1-1)*3*NATOMS+3*(J2-1)+3) - XYZ(J1*3*NATOMS+3*(J2-1)+3) )**2 
                  DUMMY=DUMMY+ADUMMY
                  IF (ADUMMY.GT.ADMAX) THEN
                     ADMAX=ADUMMY
                     JA1=J1
                     JA2=J2
                  ENDIF
               ENDIF
            ENDDO
            DUMMY=SQRT(DUMMY)
            DUMMY2=DUMMY2+DUMMY
            IF (DUMMY.GT.DMAX) THEN
               DMAX=DUMMY
               JMAX=J1
            ENDIF
            IF (DUMMY.LT.DMIN) THEN
               DMIN=DUMMY
               JMIN=J1
            ENDIF
!            IF (DEBUG) WRITE(*,'(A,I6,A,I6,A,G20.10)')' intlbfgs> distance between images ', &
!  &                                                  J1,' and ',J1+1,' is ',DUMMY
!            IF (DEBUG) WRITE(*,'(A,G20.10,A,I6,A,2I6)')' intlbfgs> largest atomic distance between images so far is ', &
!  &                                                  SQRT(ADMAX),' for atom ',JA2,' and images ',JA1,JA1+1
         ENDDO
         IF (DEBUG) WRITE(*,'(A,G20.10,A,I6,A,2I6,A,I6)')' intlbfgs> largest atomic distance between images is ', &
  &                                                  SQRT(ADMAX),' for atom ',JA2,' and images ',JA1,JA1+1,' total images=',INTIMAGE
!        IF (DEBUG) WRITE(*,'(A,G20.10,A,2I6)')' intlbfgs> largest image separation is ', &
! &                                                  DMAX,' for images ',JMAX,JMAX+1
!        IF (DEBUG) WRITE(*,'(A,G20.10,A,G20.10)') 'intlbfgs> Mean image separation=',DUMMY2/(INTIMAGE+1),' per active atom=',DUMMY2/((INTIMAGE+1)*NACTIVE)
!        IF (SQRT(ADMAX).GT.IMSEPMAX) THEN
!           KINT=MIN(1.0D6,KINT*1.1D0)
!        ELSE
!           KINT=MAX(1.0D-6,KINT/1.1D0)
!        ENDIF
!        WRITE(*,'(A,G20.10)') 'intlbfgs> Spring constant is now ',KINT
      ENDIF
   ENDIF
!
! End of add/subtract images block.
!
!
! The new QCILPERMDIST check should be used instead of these lines
!
!  IF (QCIPERMCHECK.AND.(MOD(NITERDONE,QCIPERMCHECKINT).EQ.0)) THEN
!     LDEBUG=.FALSE.
!     DO J2=2,INTIMAGE+2
!        CALL MINPERMDIST(XYZ((J2-2)*3*NATOMS+1:(J2-1)*3*NATOMS),XYZ((J2-1)*3*NATOMS+1:J2*3*NATOMS),NATOMS,LDEBUG, &
! &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DIST,DIST2,RIGIDBODY,RMAT)
!     ENDDO
!  ENDIF

   IF (.NOT.SWITCHED) THEN
      IF (MOD(NITERDONE,CHECKREPINTERVAL).EQ.0) CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),0,1)
      IF (QCIADDREP.GT.0) THEN
         CALL CONGRAD3(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ELSEIF (CHECKCONINT) THEN
         CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ELSE
         CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ENDIF

      IF ((ETOTAL-EOLD.LT.1.0D100).OR.ADDATOM) THEN ! MAXERISE effectively set to 1.0D100 here
         EOLD=ETOTAL
         GLAST(1:D)=G(1:D)
         XSAVE(1:D)=X(1:D)
      ELSE
         NDECREASE=NDECREASE+1
         IF (NDECREASE.GT.5) THEN
            NFAIL=NFAIL+1
            WRITE(*,'(A,I6)') ' intlbfgs> WARNING *** in lbfgs cannot find a lower energy, NFAIL=',NFAIL
            X(1:D)=XSAVE(1:D)
            G(1:D)=GLAST(1:D)
         ELSE
            X(1:D)=XSAVE(1:D)
            G(1:D)=GLAST(1:D)
            STP(1:D)=STP(1:D)/10.0D0
            WRITE(*,'(A,G25.15,A,G25.15,A)') ' intlbfgs> energy increased from ',EOLD,' to ',ETOTAL, &
     &          ' decreasing step size'
            GOTO 20
         ENDIF
      ENDIF
      ADDATOM=.FALSE.
   ELSE ! combine constraint and true potentials
      IF (MOD(NITERDONE,CHECKREPINTERVAL).EQ.0) CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),0,1)
      ETOTALTMP=0.0D0
      IF (INTCONFRAC.NE.0.0D0) THEN
         DO J4=2,INTIMAGE+1
            IF (CHRMMT) CALL UPDATENBONDS(XYZ((3*NATOMS)*(J4-1)+1:(3*NATOMS)*J4))
            CALL POTENTIAL(XYZ((3*NATOMS)*(J4-1)+1:(3*NATOMS)*J4),EEE(J4),GGG((3*NATOMS)*(J4-1)+1:(3*NATOMS)*J4), &
  &                                    .TRUE.,.FALSE.,RMS,.TRUE.,.FALSE.)
            ETOTALTMP=ETOTALTMP+EEE(J4)
         ENDDO
      ENDIF
      EEETMP(1:INTIMAGE+2)=EEE(1:INTIMAGE+2)
      MYGTMP(1:D)=G(1:D)
      IF (USEFRAC.LT.1.0D0) THEN
         IF (QCIADDREP.GT.0) THEN
            CALL CONGRAD3(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSEIF (CHECKCONINT) THEN
            CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSE
            CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ENDIF
      ELSE
         ETOTAL=0.0D0
         G(1:D)=0.0D0
      ENDIF
      ETOTAL=USEFRAC*ETOTALTMP+(1.0D0-USEFRAC)*ETOTAL
      G(1:D)=USEFRAC*MYGTMP(1:D)+(1.0D0-USEFRAC)*G(1:D)
      RMS=SUM(G(1:D)**2)
      RMS=SQRT(RMS/((3*NATOMS)*INTIMAGE))
      EEE(1:INTIMAGE+2)=USEFRAC*EEETMP(1:INTIMAGE+2)+(1.0D0-USEFRAC)*EEE(1:INTIMAGE+2)
      WORST=-1.0D100
      DO J4=2,INTIMAGE+1
         IF (EEE(J4).GT.WORST) WORST=EEE(J4)
      ENDDO
      IF (DEBUG) WRITE(*,'(A,G20.10,A,I8)') ' intlbfgs> Highest QCI image energy=',WORST,' images=',INTIMAGE
   ENDIF
   IF (ETOTAL/INTIMAGE.LT.COLDFUSIONLIMIT) THEN
      WRITE(*,'(A,2G20.10)') ' intlbfgs> Cold fusion diagnosed - step discarded, energy, limit=',ETOTAL/INTIMAGE,COLDFUSIONLIMIT
      DEALLOCATE(CONI,CONJ,CONDISTREF,REPI,REPJ,NREPI,NREPJ,REPCUT,NREPCUT,CONCUT,CONOFFLIST,CONOFFTRIED)
      DEALLOCATE(TRUEEE, EEETMP, MYGTMP, GTMP, &
  &              DIAG, STP, SEARCHSTEP, GDIF,GLAST, XSAVE, XYZ, GGG, CHECKG, IMGFREEZE, EEE, STEPIMAGE)
      QCIIMAGE=INTIMAGE
      INTIMAGE=INTIMAGESAVE
      RETURN
   ENDIF

   STEPTOT = SUM(STEPIMAGE)/INTIMAGE

   MAXRMS=-1.0D0
   MAXEEE=-1.0D100
   MINEEE=1.0D100
   SUMEEE=0.0D0
   SUMEEE2=0.0D0
   DO J1=2,INTIMAGE+1
      SUMEEE=SUMEEE+EEE(J1)
      SUMEEE2=SUMEEE2+EEE(J1)**2
      IF (EEE(J1).GT.MAXEEE) THEN
         MAXEEE=EEE(J1)
         JMAXEEE=J1
      ENDIF
      IF (EEE(J1).LT.MINEEE) THEN
         MINEEE=EEE(J1)
      ENDIF
      DUMMY=0.0D0
      DO J2=1,3*NATOMS
         DUMMY=DUMMY+GGG(3*NATOMS*(J1-1)+J2)**2
      ENDDO
      IF (DUMMY.GT.MAXRMS) THEN
         MAXRMS=DUMMY
         JMAXRMS=J1
      ENDIF
   ENDDO
   MAXRMS=SQRT(MAXRMS/(3*NACTIVE))
   SUMEEE=SUMEEE/INTIMAGE
   SUMEEE2=SUMEEE2+EEE(J1)**2
   SUMEEE2=SUMEEE2/INTIMAGE
   SIGMAEEE=SQRT(SUMEEE2-SUMEEE**2)
   REMOVEIMAGE=.FALSE.
!    IF (ABS(MAXEEE-SUMEEE).GT.3.0D0*SIGMAEEE) THEN
! !  IF (MAXEEE.GT.1.0D2*MINEEE) THEN
      WRITE(*,'(A,I8,A,G20.10,A,G20.10,A)') ' intlbfgs> Highest image ',JMAXEEE,' energy ',MAXEEE,' is ',ABS(MAXEEE-SUMEEE)/SIGMAEEE, &
  &                        ' sigma from the mean'
!       REMOVEIMAGE=.TRUE.
!       IF (NITERDONE-NLASTGOODE.LT.20) REMOVEIMAGE=.FALSE.
!    ENDIF

   IF (DEBUG) THEN
!     WRITE(*,'(A,2G20.10)') ' intlbfgs> mean and sigma of image energies=',SUMEEE,SIGMAEEE
!     WRITE(*,'(A,I6,2G20.10,3(G20.10,I8))') ' intlbfgs> steps: ',NITERDONE,ETOTAL/INTIMAGE,RMS,STEPTOT,NACTIVE, &
! &                                                        MAXEEE,JMAXEEE,MAXRMS,JMAXRMS
      WRITE(*,'(A,I6,5G20.10,I10,I6)') ' intlbfgs> steps: ',NITERDONE,CONVERGECONTEST,CONVERGEREPTEST,FCONTEST,FREPTEST,STEPTOT,NACTIVE,INTIMAGE+2
      CALL FLUSH(6)
   ENDIF

   IF (.NOT.SWITCHED) THEN
!     IF ((NITERDONE-NLASTGOODE.GT.INTRELSTEPS).AND.((ETOTAL.GT.LASTGOODE).OR.(ETOTAL/INTIMAGE.GT.MAXCONE*1.0D8))) THEN
      IF (.FALSE.) THEN ! no backtracking
         WRITE(*,'(2(A,I6))') ' intlbfgs> Backtracking ',NBACKTRACK,' steps, current active atoms=',NACTIVE
         NTRIES(NEWATOM)=NTRIES(NEWATOM)+1
         IF (FREEZENODEST) IMGFREEZE(1:INTIMAGE)=.FALSE.
!
! Backtrack by removing the last NBACKTRACK atoms along with their active constraints and
! repulsions.
!
         NOFF=0
         DO J1=1,NBACKTRACK
            NDUMMY=TURNONORDER(NACTIVE-J1+1)
            IF (INTFROZEN(NDUMMY)) THEN
               IF (DEBUG) WRITE(*,'(A,I6,A,2I6)') ' intlbfgs> Not turning off frozen active atom ',NDUMMY
               CYCLE
            ENDIF
            IF (DEBUG) WRITE(*,'(A,I6,A,2I6)') ' intlbfgs> Turning off active atom ',NDUMMY
            DO J2=1,NCONSTRAINT
               IF (.NOT.CONACTIVE(J2)) CYCLE 
               IF ((CONI(J2).EQ.NDUMMY).OR.(CONJ(J2).EQ.NDUMMY)) THEN
                  CONACTIVE(J2)=.FALSE.
                  IF (DEBUG) WRITE(*,'(A,I6,A,2I6)') ' intlbfgs> Turning off constraint ',J2,' for atoms ',CONI(J2),CONJ(J2)
               ENDIF
            ENDDO
            ATOMACTIVE(NDUMMY)=.FALSE.
            NOFF=NOFF+1
         ENDDO
         NACTIVE=NACTIVE-NOFF
         NDUMMY=1
         NREPULSIVE=0
         DO J1=1,NATOMS
! 
! Make a list of repelling atoms here and then use it
! CONI(J2) is always less than CONJ(J2) so we only need to
! cycle over a given range of constraints and continue from
! where we left off for the next atom j1
!  
            ADDREP(1:J1+INTREPSEP)=.FALSE.
            ADDREP(J1+INTREPSEP+1:NATOMS)=.TRUE. ! no repulsion for atoms too close in sequence
            IF (INTFROZEN(J1)) THEN
               DO J2=J1+INTREPSEP+1,NATOMS
                  IF (INTFROZEN(J2)) ADDREP(J2)=.FALSE.
               ENDDO
            ENDIF
            addloop: DO J2=NDUMMY,NCONSTRAINT
               IF (CONI(J2).EQ.J1) THEN
                  ADDREP(CONJ(J2))=.FALSE.
               ELSE
                  NDUMMY=J2 ! for next atom
                  EXIT addloop
               ENDIF
            ENDDO addloop
            rep2: DO J2=J1+INTREPSEP+1,NATOMS

               IF (.NOT.ADDREP(J2)) CYCLE
!
! Don't we need to check atomactive here for backtracking?
!
!              IF (.NOT.ATOMACTIVE(J2)) CYCLE 

               DMIN=1.0D100
               DO J3=1,INTIMAGE+2,INTIMAGE+1 ! only consider the end-point distances
                  DF=SQRT((XYZ((J3-1)*3*NATOMS+3*(J2-1)+1)-XYZ((J3-1)*3*NATOMS+3*(J1-1)+1))**2+ &
    &                     (XYZ((J3-1)*3*NATOMS+3*(J2-1)+2)-XYZ((J3-1)*3*NATOMS+3*(J1-1)+2))**2+ &
    &                     (XYZ((J3-1)*3*NATOMS+3*(J2-1)+3)-XYZ((J3-1)*3*NATOMS+3*(J1-1)+3))**2)
                  IF (DF.LT.DMIN) DMIN=DF
               ENDDO

               NREPULSIVE=NREPULSIVE+1
               IF (NREPULSIVE.GT.NREPMAX) CALL REPDOUBLE
               REPI(NREPULSIVE)=J1
               REPJ(NREPULSIVE)=J2
! 
! Use the minimum of the end point distances and INTCONSTRAINREPCUT for each contact.
!
               REPCUT(NREPULSIVE)=MIN(DMIN-1.0D-3,INTCONSTRAINREPCUT)
            ENDDO rep2
         ENDDO


         NBACKTRACK=MAX(MIN(MIN(1.0D0*(NBACKTRACK+1),1.0D0*50),1.0D0*(NACTIVE-2-NQCIFREEZE)),1.0D0)
!        IF (DEBUG) WRITE(*,'(A,I6)') ' intlbfgs> Number of atoms to backtrack is now ',NBACKTRACK
         NDUMMY=0
         DO J1=1,NATOMS
            IF (ATOMACTIVE(J1)) NDUMMY=NDUMMY+1
         ENDDO
         IF (NDUMMY.NE.NACTIVE) THEN
            WRITE(*,'(A,I6)') ' intlbfgs> ERROR *** inconsistency in number of active atoms. ',NDUMMY,' should be ',NACTIVE
            DO J1=1,NATOMS
               IF (ATOMACTIVE(J1)) WRITE(*,'(A,I6)') ' active atom ',J1
            ENDDO
            STOP
         ENDIF
         ADDATOM=.TRUE.

         CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),0,1)
         IF (QCIADDREP.GT.0) THEN
            CALL CONGRAD3(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSEIF (CHECKCONINT) THEN
            CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ELSE
            CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
         ENDIF
      ENDIF
      LASTGOODE=ETOTAL
   ENDIF

   EXITSTATUS=0
   INTDGUESS=DIAG(1) ! should be ok for subsequent runs of the same system DJW
!  IF ((.NOT.SWITCHED).AND.(MAXRMS<=INTRMSTOL).AND.(NITERDONE>1).AND.(CONVERGECONTEST.LT.MAXCONE).AND.(CONVERGEREPTEST.LT.MAXCONE)) EXITSTATUS=1 
   IF ((.NOT.SWITCHED).AND.(FCONTEST.LT.INTRMSTOL).AND.(FREPTEST.LT.INTRMSTOL).AND.(NITERDONE>1) &
  &                   .AND.(CONVERGECONTEST.LT.MAXCONE).AND.(CONVERGEREPTEST.LT.MAXCONE)) EXITSTATUS=1 
   IF (SWITCHED.AND.(MAXRMS<=CONVR).AND.NITERDONE>1) EXITSTATUS=1 
   IF (NITERDONE==NSTEPSMAX) EXITSTATUS=2
   IF ((.NOT.SWITCHED).AND.(MOD(NITERDONE,INTRELSTEPS).EQ.0)) EXITSTATUS=1 ! Add an atom every INTRELSTEPS !!! DJW
!  PRINT '(A,2G20.10,3I8)','MAXRMS,INTRMSTOL,NITERDONE,NITERDONE,NSTEPSMAX=',MAXRMS,INTRMSTOL,NITERDONE,NITERDONE,NSTEPSMAX

   IF (EXITSTATUS > 0) THEN  
      PRINT *,'here A'
      IF ((.NOT.SWITCHED).AND.(EXITSTATUS.EQ.1)) THEN ! add active atom or restart with true potential on
!        IF (ETOTAL/INTIMAGE.GT.MAXCONE*MAX(0.3D0,NACTIVE*1.0D0/(NATOMS*1.0D0))) GOTO 777
!        PRINT '(A,3G20.10)','MAXEEE,MAXCONE,scaled=',MAXEEE,MAXCONE,MAXCONE*MAX(0.3D0,NACTIVE*1.0D0/(NATOMS*1.0D0))
!        IF (MAXEEE.GT.MAXCONE*MAX(0.3D0,NACTIVE*1.0D0/(NATOMS*1.0D0))) GOTO 777
         IF (NACTIVE.LT.NATOMS) THEN 
            ADDATOM=.TRUE.
            GOTO 777
         ENDIF
         CALL MYCPU_TIME(FTIME,.FALSE.)
         WRITE(*,'(A,I6,A,F12.6,A,I6,A,G20.10)') ' intlbfgs> switch on true potential at step ',NITERDONE, &
  &                                     ' fraction=',INTCONFRAC,' images=',INTIMAGE,' time=',FTIME-STIME
         CALL INTRWG(NACTIVE,NITERDONE,INTIMAGE,XYZ,TURNONORDER,NCONOFF)
         CALL WRITEPROFILE(NITERDONE,EEE,INTIMAGE)
         WRITE(*,'(A,I6,A,F15.6)') ' intlbfgs> Allowing ',INTCONSTEPS,' further optimization steps'
         DO J1=1,NATOMS
            IF (.NOT.ATOMACTIVE(J1)) THEN
               WRITE(*,'(A,I6,A,I6,A)') ' intlbfgs> ERROR *** number of active atoms=',NACTIVE,' but atom ',J1,' is not active'
            ENDIF
         ENDDO
         NSTEPSMAX=NITERDONE+INTCONSTEPS
         SWITCHED=.TRUE.
         RMS=INTRMSTOL*10.0D0 ! to prevent premature convergence
         G(1:(3*NATOMS)*INTIMAGE)=INTRMSTOL*10.0D0
         USEFRAC=INTCONFRAC
         GOTO 777
      ELSEIF ((.NOT.SWITCHED).AND.(EXITSTATUS.EQ.2)) THEN 
         WRITE(*,'(A,I6)') ' intlbfgs> QCI ERROR *** number of active atoms at final step=',NACTIVE
         CALL FLUSH(6)
         QCIIMAGE=INTIMAGE
         RETURN
      ELSEIF (DEBUG) THEN
         WRITE(*,'(A,I6,A,I6)') 'intlbfgs> energies for images:'
         WRITE(*,'(I6,F20.10)') (J2,EEE(J2),J2=1,INTIMAGE+2)
      ENDIF
      EXIT
   ENDIF
   777 CONTINUE
!
! Compute the new step and gradient change
!
   NPT=POINT*D
   SEARCHSTEP(POINT,:) = STP*SEARCHSTEP(POINT,:)
   GDIF(POINT,:)=G-GTMP
   
   POINT=POINT+1; IF (POINT==INTMUPDATE) POINT=0

   IF (DUMPINTXYZ.AND.MOD(NITERDONE,DUMPINTXYZFREQ)==0) CALL INTRWG(NACTIVE,NITERDONE,INTIMAGE,XYZ,TURNONORDER,NCONOFF)
   IF (DUMPINTEOS.AND.MOD(NITERDONE,DUMPINTEOSFREQ)==0) CALL WRITEPROFILE(NITERDONE,EEE,INTIMAGE)
!  IF (NITERDONE.GT.0) CALL INTRWG2(NACTIVE,NITERDONE,INTIMAGE,XYZ,TURNONORDER,NCONOFF) !!! DEBUG DJW

   NITERDONE=NITERDONE+1
   NITERUSE=NITERUSE+1

   IF (NITERDONE.GT.NSTEPSMAX) EXIT
   IF (NACTIVE.EQ.NATOMS) THEN
      IF (.NOT.SWITCHED) THEN
         CALL MYCPU_TIME(FTIME,.FALSE.)
         WRITE(*,'(A,I6,A,F12.6,A,I6,A,F10.1)') ' intlbfgs> switch on true potential at step ',NITERDONE, &
  &                                     ' fraction=',INTCONFRAC,' images=',INTIMAGE,' time=',FTIME-STIME
         WRITE(*,'(A,I6,A,F15.6)') ' intlbfgs> Allowing ',INTCONSTEPS,' further optimization steps'
         DO J1=1,NATOMS
            IF (.NOT.ATOMACTIVE(J1)) THEN
               WRITE(*,'(A,I6,A,I6,A)') ' intlbfgs> ERROR *** number of active atoms=',NACTIVE,' but atom ',J1,' is not active'
            ENDIF
         ENDDO
         NSTEPSMAX=NITERDONE+INTCONSTEPS
         SWITCHED=.TRUE.
         IF (FREEZENODEST) THEN
            IMGFREEZE(1:INTIMAGE)=.FALSE.
         ENDIF
         RMS=INTRMSTOL*10.0D0 ! to prevent premature convergence
         USEFRAC=INTCONFRAC
      ENDIF
   ENDIF

ENDDO ! end of main do loop over counter NITERDONE

      CALL FLUSH(6)

IF (.NOT.SWITCHED) THEN 
   WRITE(*,'(A,I6,A)') ' intlbfgs> QCI DID NOT CONVERGE number of active atoms at final step=',NACTIVE,' no potential switch'
ENDIF
IF (EXITSTATUS.EQ.1) THEN
   WRITE(*,'(A,I6,A,G20.10,A,G15.8,A,I4)') ' intlbfgs> Converged after ',NITERDONE,' steps, energy/image=',ETOTAL/INTIMAGE, &
  &                               ' RMS=',RMS,' images=',INTIMAGE
ELSEIF (EXITSTATUS.EQ.2) THEN
   WRITE(*,'(A,I6,A,G20.10,A,G15.8,A,I4)') ' intlbfgs> After ',NITERDONE,' steps, energy/image=',ETOTAL/INTIMAGE, &
  &                               ' RMS=',RMS,' images=',INTIMAGE
ENDIF
678 CONTINUE

! CALL INTRWG(NACTIVE,NITERDONE,INTIMAGE,XYZ,TURNONORDER,NCONOFF)
! CALL WRITEPROFILE(NITERDONE,EEE,INTIMAGE)

IF (DEBUG) WRITE(*,'(A,G20.10)') 'intlbfgs> WORST=',WORST

BESTWORST=WORST
BESTINTIMAGE=INTIMAGE
IF (ALLOCATED(QCIXYZ)) DEALLOCATE(QCIXYZ)
ALLOCATE(QCIXYZ(3*NATOMS*(INTIMAGE+2)))
QCIXYZ(1:3*NATOMS*(INTIMAGE+2))=XYZ(1:3*NATOMS*(INTIMAGE+2))
WRITE(*,'(A,I8,A,G20.10)') 'intlbfgs> retaining ',INTIMAGE,' QCI images, highest energy=',BESTWORST

CALL INTRWG(NACTIVE,0,INTIMAGE,XYZ,TURNONORDER,NCONOFF)
CALL WRITEPROFILE(0,EEE,INTIMAGE)

DEALLOCATE(CONI,CONJ,CONDISTREF,REPI,REPJ,NREPI,NREPJ,REPCUT,NREPCUT,CONCUT,CONOFFLIST,CONOFFTRIED)
DEALLOCATE(TRUEEE, EEETMP, MYGTMP, GTMP, &
  &      DIAG, STP, SEARCHSTEP, GDIF,GLAST, XSAVE, XYZ, GGG, CHECKG, IMGFREEZE, EEE, STEPIMAGE)
QCIIMAGE=INTIMAGE
INTIMAGE=INTIMAGESAVE
IF (ALLOCATED(CONLIST)) DEALLOCATE(CONLIST)
IF (ALLOCATED(NCONATOM)) DEALLOCATE(NCONATOM)
IF (ALLOCATED(COMMONCON)) DEALLOCATE(COMMONCON)

IF (QCISTOP) STOP 

END SUBROUTINE INTLBFGS
!
! Neighbour list for repulsions to reduce cost of constraint potential.
!
SUBROUTINE CHECKREP(INTIMAGE,XYZ,NOPT,NNSTART,NSTART)
USE KEY,ONLY : NREPI, NREPJ, NREPCUT, NNREPULSIVE, NREPULSIVE, REPI, REPJ, REPCUT, CHECKREPCUTOFF, &
  &                INTFROZEN, NNREPULSIVE, intconstraintrep
USE COMMONS, ONLY : DEBUG
USE PORFUNCS
IMPLICIT NONE
INTEGER JJ, KK, NI1, NJ1, NI2, NJ2, INTIMAGE, NOPT, NI, NJ, NNSTART, NSTART
DOUBLE PRECISION LDIST, XYZ(NOPT*(INTIMAGE+2)),COMPARE
DOUBLE PRECISION R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,DMIN
LOGICAL NOINT

IF (INTCONSTRAINTREP.EQ.0) THEN
   NNREPULSIVE=0
   RETURN
ENDIF

NNREPULSIVE=NNSTART
DO JJ=NSTART,NREPULSIVE
   COMPARE=(CHECKREPCUTOFF*REPCUT(JJ))**2
   NI=REPI(JJ)
   NJ=REPJ(JJ)
   DO KK=1,INTIMAGE+2 ! first check for standard distances within threshold
      LDIST=(XYZ((KK-1)*NOPT+3*(NI-1)+1)-XYZ((KK-1)*NOPT+3*(NJ-1)+1))**2 &
  &        +(XYZ((KK-1)*NOPT+3*(NI-1)+2)-XYZ((KK-1)*NOPT+3*(NJ-1)+2))**2 &
  &        +(XYZ((KK-1)*NOPT+3*(NI-1)+3)-XYZ((KK-1)*NOPT+3*(NJ-1)+3))**2
      IF (LDIST.LT.COMPARE) THEN
         NNREPULSIVE=NNREPULSIVE+1
         NREPI(NNREPULSIVE)=NI
         NREPJ(NNREPULSIVE)=NJ
         NREPCUT(NNREPULSIVE)=REPCUT(JJ)
!        IF ((REPI(JJ).EQ.2024).OR.(REPJ(JJ).EQ.2024)) THEN
!           WRITE(*,'(A)') 'checkrep> KK,JJ,ldist < compare : active'
!        ENDIF
         GOTO 246
      ENDIF
   ENDDO 
!
! We don't check for internal minima in repulsions in congrad now unless both distances are
! within threshold.
!
!  COMPARE=CHECKREPCUTOFF*REPCUT(JJ)
!  DO KK=2,INTIMAGE+2 ! now check internal minima within threshold
!     DMIN=1.0D10
!     NI2=NOPT*(KK-2)+3*(NI-1)
!     NI1=NOPT*(KK-1)+3*(NI-1)
!     NJ2=NOPT*(KK-2)+3*(NJ-1)
!     NJ1=NOPT*(KK-1)+3*(NJ-1)
!     R1AX=XYZ(NI2+1); R1AY=XYZ(NI2+2); R1AZ=XYZ(NI2+3)
!     R1BX=XYZ(NJ2+1); R1BY=XYZ(NJ2+2); R1BZ=XYZ(NJ2+3)
!     R2AX=XYZ(NI1+1); R2AY=XYZ(NI1+2); R2AZ=XYZ(NI1+3)
!     R2BX=XYZ(NJ1+1); R2BY=XYZ(NJ1+2); R2BZ=XYZ(NJ1+3)
!     CALL INTMINONLY(R1AX,R1AY,R1AZ,R2AX,R2AY,R2AZ,R1BX,R1BY,R1BZ,R2BX,R2BY,R2BZ,DMIN,NOINT)

!     IF (NOINT) CYCLE
!     IF (DMIN.LT.COMPARE) THEN
!        NNREPULSIVE=NNREPULSIVE+1
!        NREPI(NNREPULSIVE)=NI
!        NREPJ(NNREPULSIVE)=NJ
!        NREPCUT(NNREPULSIVE)=REPCUT(JJ)
!        GOTO 246
!     ENDIF
!  ENDDO 
246 CONTINUE
ENDDO
IF (DEBUG) WRITE(*,'(A,2I8)') ' checkrep> number of active repulsions and total=',NNREPULSIVE,NREPULSIVE

END SUBROUTINE CHECKREP

SUBROUTINE INTRWG(NACTIVE,NITER,INTIMAGE,XYZ,TURNONORDER,NCONOFF)
USE PORFUNCS
USE KEY,ONLY: STOCKT,STOCKAAT, RBAAT, ATOMACTIVE, NCONSTRAINT, CONACTIVE, NREPULSIVE, NNREPULSIVE, REPI, REPJ, REPCUT, NREPCUT, &
  &           NREPMAX, NREPI, NREPJ, INTFROZEN, CONOFFLIST,CONOFFTRIED, KINT, INTCONSTRAINTREP
USE COMMONS, ONLY: NATOMS, DEBUG
IMPLICIT NONE
INTEGER NCONOFF
CHARACTER(LEN=10) :: XYZFILE   = 'int.xyz   '
CHARACTER(LEN=10) :: QCIFILE   = 'QCIdump   '
INTEGER,INTENT(IN) :: NITER, TURNONORDER(NATOMS)
INTEGER :: J1,J2,INTIMAGE,J3,NACTIVE,LUNIT,GETUNIT
CHARACTER(LEN=80) :: FILENAME,DUMMYS
DOUBLE PRECISION XYZ((3*NATOMS)*(INTIMAGE+2))

FILENAME=XYZFILE

! IF (NITER.GT.0) THEN
!    WRITE(DUMMYS,'(I8)') NITER
!    FILENAME='int.' // TRIM(ADJUSTL(DUMMYS)) // '.xyz' ! so that vmd recognises the file type!
! ENDIF
LUNIT=GETUNIT()
OPEN(UNIT=LUNIT,FILE=TRIM(ADJUSTL(FILENAME)),STATUS='replace')
!  WRITE(*,'(A,3L5)') 'intrwg> here B atomactive 2113 2115 2123=',ATOMACTIVE(2113),ATOMACTIVE(2115),ATOMACTIVE(2123)
DO J2=1,INTIMAGE+2
!  WRITE(LUNIT,'(i4/)') NACTIVE
   WRITE(LUNIT,'(i4/)') NATOMS
   DO J3=1,NATOMS
      IF (ATOMACTIVE(J3)) THEN
         WRITE(LUNIT,'(A5,1X,3F20.10)') 'LA   ',XYZ((J2-1)*3*NATOMS+3*(J3-1)+1),XYZ((J2-1)*3*NATOMS+3*(J3-1)+2), &  
  &                                                                   XYZ((J2-1)*3*NATOMS+3*(J3-1)+3)  
      ELSE
         WRITE(LUNIT,'(A5,1X,3F20.10)') 'DU   ',XYZ((J2-1)*3*NATOMS+3*(J3-1)+1),XYZ((J2-1)*3*NATOMS+3*(J3-1)+2), &  
  &                                                                   XYZ((J2-1)*3*NATOMS+3*(J3-1)+3)  
      ENDIF
   ENDDO
ENDDO

WRITE(*,*) 'rwg> Interpolated image coordinates were saved to xyz file "'//TRIM(FILENAME)//'"'

CLOSE(LUNIT)

FILENAME=QCIFILE
LUNIT=GETUNIT()
OPEN(UNIT=LUNIT,FILE=TRIM(ADJUSTL(FILENAME)),STATUS='replace')

IF (DEBUG) WRITE(*,'(A,I10,A)') ' intlbfgs> dumping state for ',NACTIVE,' active atoms'
WRITE(LUNIT,'(I10)') NACTIVE
! IF (DEBUG)   WRITE(*,'(A,I10,A)') ' intlbfgs> dumping spring constant and repulsive prefactor'
WRITE(LUNIT,'(2G20.10)') KINT, INTCONSTRAINTREP
! WRITE(*,'(A,I10,A)') ' intlbfgs> dumping turnonorder for ',NACTIVE,' active atoms'
WRITE(LUNIT,'(12I8)') TURNONORDER(1:NACTIVE)
! WRITE(*,'(A)') ' intlbfgs> dumping atomactive'
WRITE(LUNIT,'(12L5)') ATOMACTIVE(1:NATOMS)
WRITE(LUNIT,'(I10)') NCONSTRAINT
! WRITE(*,'(A,I10,A)') ' intlbfgs> dumping conactive for ',NCONSTRAINT,' constraints'
WRITE(LUNIT,'(12L5)') CONACTIVE(1:NCONSTRAINT)

WRITE(LUNIT,'(3I12,G20.10)') NREPULSIVE,NNREPULSIVE,NREPMAX
! WRITE(*,'(A,3I10,G20.10)') 'intlbfgs> dumping NREPULSIVE,NNREPULSIVE,NREPMAX=',NREPULSIVE,NNREPULSIVE,NREPMAX

WRITE(LUNIT,'(12I8)') REPI(1:NREPULSIVE)
! WRITE(*,'(A)') ' intlbfgs> dumped REPI:'
WRITE(LUNIT,'(12I8)') REPJ(1:NREPULSIVE)
! WRITE(*,'(A)') ' intlbfgs> dumped REPJ:'
WRITE(LUNIT,'(12I8)') NREPI(1:NNREPULSIVE)
! WRITE(*,'(A)') ' intlbfgs> dumped NREPI:'
WRITE(LUNIT,'(12I8)') NREPJ(1:NNREPULSIVE)
! WRITE(*,'(A)') ' intlbfgs> dumped NREPJ:'

WRITE(LUNIT,'(6G20.10)') REPCUT(1:NREPULSIVE)
! WRITE(*,'(A)') ' intlbfgs> dumped REPCUT:'
WRITE(LUNIT,'(6G20.10)') NREPCUT(1:NNREPULSIVE)
! WRITE(*,'(A)') ' intlbfgs> dumped NREPCUT:'

   WRITE(LUNIT,'(12L5)') INTFROZEN(1:NATOMS)
   ! WRITE(*,'(A)') ' intlbfgs> dumped INTFROZEN'

   WRITE(LUNIT,'(I8)') NCONOFF
   IF (NCONOFF.GT.0) WRITE(LUNIT,'(12I8)') CONOFFLIST(1:NCONOFF)
   IF (NCONOFF.GT.0) WRITE(LUNIT,'(12L5)') CONOFFTRIED(1:NCONSTRAINT)
   ! WRITE(*,'(A)') ' intlbfgs> dumped NCONOFF and CONOFFLIST'

CLOSE(LUNIT)

END SUBROUTINE INTRWG

SUBROUTINE WRITEPROFILE(NITER,EEE,INTIMAGE)
IMPLICIT NONE 
INTEGER,INTENT(IN) :: NITER, INTIMAGE
INTEGER :: I,LUNIT,GETUNIT
DOUBLE PRECISION :: EEE(INTIMAGE+2)
CHARACTER(LEN=20) :: FILENAME

LUNIT=GETUNIT()
! IF (NITER.GT.0) THEN
!    WRITE(FILENAME,'(I8)') NITER
!    FILENAME='int.EofS.' // TRIM(ADJUSTL(FILENAME))
! ELSE   
   FILENAME='int.EofS'
! ENDIF
OPEN(UNIT=LUNIT,FILE=FILENAME,STATUS='replace')

WRITE(UNIT=LUNIT,FMT='(2g24.13)') EEE(1)
DO I=2,INTIMAGE+1
   WRITE(UNIT=LUNIT,FMT='(2G24.13)') EEE(I)
ENDDO
WRITE(UNIT=LUNIT,FMT='(2G24.13)') EEE(INTIMAGE+2)

CLOSE(LUNIT)
WRITE(*,'(A)') ' writeprofile> Interpolated energy profile was saved to file "'//trim(filename)//'"'

END SUBROUTINE WRITEPROFILE

SUBROUTINE DOADDATOM(NCONSTRAINT,NTRIES,NEWATOM,IMGFREEZE,INTIMAGE,XYZ,EEE,GGG,TURNONORDER,NITERDONE,NACTIVE,AABACK,BACKDONE)
USE KEY, ONLY : CONACTIVE, CONI, CONJ, ATOMACTIVE, CONDISTREF, REPI, REPJ, REPCUT, INTREPSEP,  &
  &             INTCONSTRAINREPCUT, NREPULSIVE, NREPMAX, MAXCONUSE, CHECKCONINT, &
  &             FREEZENODEST, NNREPULSIVE, INTFROZEN, ATOMSTORES, QCIADDACIDT, &
  &             NREPULSIVEFIX, REPIFIX, REPJFIX, REPCUTFIX, NREPI, NREPJ, NREPCUT, MAXNACTIVE, &
  &             NCONSTRAINTFIX, CONIFIX, CONJFIX, INTCONCUT, INTCONSEP, QCIRADSHIFTT, QCIRADSHIFT, &
  &             CONOFFTRIED, QCIADDREP, DOBACK, DOBACKALL, QCITRILAT
USE COMMONS, ONLY: NATOMS, DEBUG
IMPLICIT NONE
INTEGER INTIMAGE
INTEGER NBEST, NCONTOACTIVE(NATOMS),  NCONSTRAINT, J2, NTRIES(NATOMS), NEWATOM,  CONLIST(NATOMS), N1, N2, N3, &
  &     NTOADD, NADDED, NMININT, NMAXINT, TURNONORDER(NATOMS), NDUMMY, J1, J3, NITERDONE, NCONFORNEWATOM, NACTIVE
DOUBLE PRECISION DUMMY, DUMMY2, DPRAND, RANDOM, CONDIST(NATOMS), DMIN
INTEGER NDFORNEWATOM, BESTPRESERVEDN(NATOMS), ACID
DOUBLE PRECISION BESTPRESERVEDD(NATOMS), BESTCLOSESTD(NATOMS), INVDTOACTIVE(NATOMS)
LOGICAL IMGFREEZE(INTIMAGE), ADDREP(NATOMS), CHOSENACID, AABACK(NATOMS), BACKDONE, FTEST, IDONE
DOUBLE PRECISION P1(3), P2(3), P3(3), SOL1(3), SOL2(3), R1, R2, R3
DOUBLE PRECISION C1, C2, C3, VEC1(3), VEC2(3), VEC3(3), ESAVED, ESAVEC, ESAVE0
INTEGER NCFORNEWATOM, BESTCLOSESTN(NATOMS), NNREPSAVE, NREPSAVE
DOUBLE PRECISION XYZ((3*NATOMS)*(INTIMAGE+2)), XSAVED(3,INTIMAGE+2), XSAVEC(3,INTIMAGE+2), XSAVE0(3,INTIMAGE+2),FRAC,RAN1, &
  &              RMS,EEE(INTIMAGE+2),GGG((3*NATOMS)*(INTIMAGE+2)),ETOTAL,DS,DF,DNORM,D1SQ,D2SQ


NTOADD=1
NADDED=0
CHOSENACID=.FALSE.

IF (DOBACK.AND.(.NOT.BACKDONE)) THEN
   DO J1=1,NATOMS
      IF (AABACK(J1)) THEN
         IF (.NOT.ATOMACTIVE(J1)) GOTO 763
      ENDIF
   ENDDO
   IF (DEBUG) WRITE(*,'(A,I6,A)') ' intlbfgs> All backbone atoms are active'
   BACKDONE=.TRUE.
ENDIF

763   CONTINUE

!
! ATOMSTORES(J1) is the residue for atom J1
!
! AMBER12_GET_RESDATA needs a data type in amber12 interface and the number of residues
! defined in amber12interface.f90
!

!
! Save current number of repulsions and number that are active to speed up the
! calls to CHECKREP
!
NNREPSAVE=NNREPULSIVE
NREPSAVE=NREPULSIVE
542   CONTINUE
!     DUMMY=1.0D100
      NBEST=0
      NCONTOACTIVE(1:NATOMS)=0
      INVDTOACTIVE(1:NATOMS)=0.0D0
      DO J2=1,NCONSTRAINT
         IF (CONACTIVE(J2)) CYCLE     ! count new, inactive constraints
         IF (CONOFFTRIED(J2)) CYCLE   ! if we've tried turning it off, it must actually be active. Don't try again.
         IF (ATOMACTIVE(CONI(J2))) THEN
            IF (.NOT.ATOMACTIVE(CONJ(J2))) THEN
               NCONTOACTIVE(CONJ(J2))=NCONTOACTIVE(CONJ(J2))+1
               IF (1.0D0/CONDISTREF(J2).GT.INVDTOACTIVE(CONJ(J2))) INVDTOACTIVE(CONJ(J2))=1.0D0/CONDISTREF(J2)
!              INVDTOACTIVE(CONJ(J2))=INVDTOACTIVE(CONJ(J2))+1.0D0/CONDISTREF(J2)
            ENDIF
         ENDIF
         IF (ATOMACTIVE(CONJ(J2))) THEN
            IF (.NOT.ATOMACTIVE(CONI(J2))) THEN
               NCONTOACTIVE(CONI(J2))=NCONTOACTIVE(CONI(J2))+1
!              INVDTOACTIVE(CONI(J2))=INVDTOACTIVE(CONI(J2))+1.0D0/CONDISTREF(J2)
               IF (1.0D0/CONDISTREF(J2).GT.INVDTOACTIVE(CONI(J2))) INVDTOACTIVE(CONI(J2))=1.0D0/CONDISTREF(J2)
            ENDIF
         ENDIF
         IF (NCONTOACTIVE(CONI(J2)).GT.NBEST) THEN
            NBEST=NCONTOACTIVE(CONI(J2))
         ENDIF
         IF (NCONTOACTIVE(CONJ(J2)).GT.NBEST) THEN
            NBEST=NCONTOACTIVE(CONJ(J2))
         ENDIF
!        IF ((CONI(J2).EQ.115).OR.(CONJ(J2).EQ.115)) THEN
!          WRITE(*,'(A,5I6,2G20.10)') 'J2,NCONTOACTIVEI,NCONTOACTOVEJ,CONI,CONJ,NEWATOM,NBEST,IDI,IDJ=', &
!   &                             J2,NCONTOACTIVE(CONI(J2)),NCONTOACTIVE(CONJ(J2)),CONI(J2),CONJ(J2), &
!   &                             INVDTOACTIVE(CONI(J2)),INVDTOACTIVE(CONJ(J2))
!        ENDIF

      ENDDO
!
!  Choose NEWATOM stochastically. Bias towards atoms with the maximum constraints.
!  Use a normalised probability and generate a random number between 0 and 1.
!
!       DUMMY2=0.0D0
!       DO J2=1,NATOMS
!          IF (NCONTOACTIVE(J2).EQ.0) CYCLE
!          IF (ATOMACTIVE(J2)) CYCLE
! !        DUMMY2=DUMMY2+((1.0D0*NCONTOACTIVE(J2))/(1.0D0*CONDISTREF(J2)*NTRIES(J2)))**4 
! !        DUMMY2=DUMMY2+((1.0D0*INVDTOACTIVE(J2))/(1.0D0*NCONTOACTIVE(J2)*NTRIES(J2)))**4 
!          DUMMY2=DUMMY2+((1.0D0*INVDTOACTIVE(J2))/(1.0D0*NTRIES(J2)))**10 
! !        WRITE(*,'(A,I6,A,G20.10)') ' intlbfgs> Unnormalised probability for choosing atom ',J2,' is ', &
! ! &                ((1.0D0*INVDTOACTIVE(J2))/(1.0D0*NTRIES(J2)))**10
!       ENDDO
! 
!       RANDOM=DUMMY2*DPRAND()
!       DNORM=DUMMY2
!       DUMMY2=0.0D0
!       choosenew: DO J2=1,NATOMS
!          IF (NCONTOACTIVE(J2).EQ.0) CYCLE
!          IF (ATOMACTIVE(J2)) CYCLE
! !        DUMMY2=DUMMY2+((1.0D0*NCONTOACTIVE(J2))/(1.0D0*CONDISTREF(J2)*NTRIES(J2)))**4 
! !        DUMMY2=DUMMY2+((1.0D0*INVDTOACTIVE(J2))/(1.0D0*NCONTOACTIVE(J2)*NTRIES(J2)))**4 
!          DUMMY2=DUMMY2+((1.0D0*INVDTOACTIVE(J2))/(1.0D0*NTRIES(J2)))**10 
!          WRITE(*,'(A,I6,G20.10,I6,4G20.10)') 'J2,invd,ntries,prob,rand,D2,D2/norm=',J2,INVDTOACTIVE(J2),NTRIES(J2), &
!   &                ((1.0D0*INVDTOACTIVE(J2))/(1.0D0*NTRIES(J2)))**10/DNORM,RANDOM/DNORM,DUMMY2,DUMMY2/DNORM
!          IF (DUMMY2.GE.RANDOM) THEN
!             NEWATOM=J2
!             IF (DEBUG) WRITE(*,'(3(A,I6))') ' intlbfgs> Choosing new active atom ',NEWATOM,' new constraints=', &
!   &                                       NCONTOACTIVE(J2),' maximum=',NBEST
!             EXIT choosenew
!          ENDIF
!       ENDDO choosenew

!
!  Choose NEWATOM deterministically. Take the inactive atom with the shortest constrained distance.
!
      DUMMY2=1.0D100
      DO J1=1,NCONSTRAINT
         IF (CONACTIVE(J1)) CYCLE
         IF (ATOMACTIVE(CONJ(J1))) THEN
            IF (CHOSENACID.AND.(.NOT.(ATOMSTORES(CONI(J1)).EQ.ACID))) THEN
            ELSE
               IF (.NOT.ATOMACTIVE(CONI(J1))) THEN
                  IF (DOBACK.AND.(.NOT.BACKDONE).AND.(.NOT.AABACK(CONI(J1)))) THEN
                  ELSE
                     IF (CONDISTREF(J1).LT.DUMMY2) THEN
                        DUMMY2=CONDISTREF(J1)
                        NEWATOM=CONI(J1)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ELSEIF (ATOMACTIVE(CONI(J1))) THEN
            IF (CHOSENACID.AND.(.NOT.(ATOMSTORES(CONJ(J1)).EQ.ACID))) THEN
            ELSE
               IF (.NOT.ATOMACTIVE(CONJ(J1))) THEN
                  IF (DOBACK.AND.(.NOT.BACKDONE).AND.(.NOT.AABACK(CONJ(J1)))) THEN
                  ELSE
                     IF (CONDISTREF(J1).LT.DUMMY2) THEN
                        DUMMY2=CONDISTREF(J1)
                        NEWATOM=CONJ(J1)
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      IF (DEBUG) WRITE(*,'(3(A,I6),A,F15.5)') ' intlbfgs> Choosing new active atom ',NEWATOM,' new constraints=', &
  &                                       NCONTOACTIVE(NEWATOM),' maximum=',NBEST,' shortest constraint=',DUMMY2
      IF (DOBACK) WRITE(*,'(A,L5)') ' intlbfgs> AABACK=',AABACK(NEWATOM)
!     IF (DEBUG) CALL INTRWG(NACTIVE,NITERDONE,INTIMAGE,XYZ,TURNONORDER,NCONOFF)
!     IF (DEBUG) CALL WRITEPROFILE(NITERDONE,EEE,INTIMAGE)
      IF (QCIADDACIDT.AND.(.NOT.CHOSENACID).AND.(.NOT.DOBACK)) THEN
         ACID=ATOMSTORES(NEWATOM)
         CHOSENACID=.TRUE.
      ENDIF
      IF ((.NOT.CHOSENACID).AND.DOBACKALL) THEN
         ACID=ATOMSTORES(NEWATOM)
         CHOSENACID=.TRUE.
      ENDIF
          
      IF (NEWATOM*NBEST.EQ.0) THEN ! sanity check
         WRITE(*,'(A,I6,A,2I6)') ' intlbfgs> ERROR *** new active atom not set'
         STOP
      ELSE
!
!  We need a sorted list of up to 3 active atoms, sorted according to how well the
!  end point distance is preserved, even if they don't satisfy the constraint 
!  condition. We want three atoms to use for a local axis system in the interpolation.
!
!  Try sorting on the shortest average distances in the endpoint structures instead, to avoid
!  problems with distant atoms acidentally having a well-preserved distance.
!
         NDFORNEWATOM=0
         BESTPRESERVEDD(1:NATOMS)=1.0D100
!          DO J1=1,NATOMS
!             IF (ABS(J1-NEWATOM).GT.INTCONSEP) CYCLE
!             IF (.NOT.ATOMACTIVE(J1)) CYCLE
!             DS=SQRT((XYZ(3*(NEWATOM-1)+1)-XYZ(3*(J1-1)+1))**2 &
!   &                +(XYZ(3*(NEWATOM-1)+2)-XYZ(3*(J1-1)+2))**2 &
!   &                +(XYZ(3*(NEWATOM-1)+3)-XYZ(3*(J1-1)+3))**2) 
!             DF=SQRT((XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+1)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+1))**2 &
!   &                +(XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+2)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+2))**2 &
!   &                +(XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+3)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+3))**2) 
!             IF (DS.GT.INTCONCUT) CYCLE
!             IF (DF.GT.INTCONCUT) CYCLE
!             DUMMY=ABS(DS-DF)
!             NDFORNEWATOM=NDFORNEWATOM+1
!             DO J2=1,NDFORNEWATOM 
!                IF (DUMMY.LT.BESTPRESERVEDD(J2)) THEN
! !                 WRITE(*,'(A,I6,G12.4,I6,G12.4)') 'J1,DUMMY < J2,BESTPRESERVEDD: ',J1,DUMMY,J2,BESTPRESERVEDD(J2)
!                   DO J3=NDFORNEWATOM,J2+1,-1 
! !                    WRITE(*,'(A,I6,A,I6,A,G12.4)') ' moving diff and list from ',J3-1,' to ',J3, &
! !&                                               ' DIFF=',BESTPRESERVEDD(J3-1)
!                      BESTPRESERVEDD(J3)=BESTPRESERVEDD(J3-1)
!                      BESTPRESERVEDN(J3)=BESTPRESERVEDN(J3-1)
!                   ENDDO
!                   BESTPRESERVEDD(J2)=DUMMY
! !                 WRITE(*,'(A,I6,A,G12.4)') ' setting BESTPRESERVEDD element ',J2,' to ',DUMMY
!                   BESTPRESERVEDN(J2)=J1
! !                 WRITE(*,'(A,I6,A,G12.4)') ' setting BESTPRESERVEDN element ',J2,' to ',J1
!                   GOTO 653
!                ENDIF
!             ENDDO
! 653         CONTINUE
!          ENDDO
!          IF (DEBUG) THEN
!             WRITE(*,'(A,I6,A,I6,A)') ' intlbfgs> New active atom ',NEWATOM,' best preserved distances:'
!             WRITE(*,'(20I6)') BESTPRESERVEDN(1:MIN(10,NDFORNEWATOM))
!             WRITE(*,'(A,I6,A,I6,A)') ' intlbfgs> sorted differences:'
!             WRITE(*,'(10G12.4)') BESTPRESERVEDD(1:MIN(10,NDFORNEWATOM))
!          ENDIF
         IF (FREEZENODEST) IMGFREEZE(1:INTIMAGE)=.FALSE.

         NCFORNEWATOM=0
         BESTCLOSESTD(1:NATOMS)=1.0D100
         DO J1=1,NATOMS
            IF (ABS(J1-NEWATOM).GT.INTCONSEP) CYCLE
            IF (.NOT.ATOMACTIVE(J1)) CYCLE
            DS=SQRT((XYZ(3*(NEWATOM-1)+1)-XYZ(3*(J1-1)+1))**2 &
  &                +(XYZ(3*(NEWATOM-1)+2)-XYZ(3*(J1-1)+2))**2 &
  &                +(XYZ(3*(NEWATOM-1)+3)-XYZ(3*(J1-1)+3))**2) 
            DF=SQRT((XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+1)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+1))**2 &
  &                +(XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+2)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+2))**2 &
  &                +(XYZ((INTIMAGE+1)*3*NATOMS+3*(NEWATOM-1)+3)-XYZ((INTIMAGE+1)*3*NATOMS+3*(J1-1)+3))**2) 
            IF (DS.GT.INTCONCUT) CYCLE
            IF (DF.GT.INTCONCUT) CYCLE
            DUMMY=(DS+DF)/2.0D0
            NCFORNEWATOM=NCFORNEWATOM+1
            DO J2=1,NCFORNEWATOM
               IF (DUMMY.LT.BESTCLOSESTD(J2)) THEN
!                 WRITE(*,'(A,I6,G12.4,I6,G12.4)') 'J1,DUMMY < J2,BESTCLOSESTD: ',J1,DUMMY,J2,BESTCLOSESTD(J2)
                  DO J3=NCFORNEWATOM,J2+1,-1
!                    WRITE(*,'(A,I6,A,I6,A,G12.4)') ' moving diff and list from ',J3-1,' to ',J3, &
!&                                               ' DIFF=',BESTCLOSESTD(J3-1)
                     BESTCLOSESTD(J3)=BESTCLOSESTD(J3-1)
                     BESTCLOSESTN(J3)=BESTCLOSESTN(J3-1)
                  ENDDO
                  BESTCLOSESTD(J2)=DUMMY
!                 WRITE(*,'(A,I6,A,G12.4)') ' setting BESTCLOSESTD element ',J2,' to ',DUMMY
                  BESTCLOSESTN(J2)=J1
!                 WRITE(*,'(A,I6,A,G12.4)') ' setting BESTCLOSESTN element ',J2,' to ',J1
                  GOTO 659
               ENDIF
            ENDDO
659         CONTINUE
         ENDDO
         IF (DEBUG) THEN
            WRITE(*,'(A,I6,A,I6,A)') ' intlbfgs> New active atom ',NEWATOM,' closest average distances in endpoints:'
            WRITE(*,'(20I6)') BESTCLOSESTN(1:MIN(10,NCFORNEWATOM))
            WRITE(*,'(A,I6,A,I6,A)') ' intlbfgs> sorted average distances:'
            WRITE(*,'(10G12.4)') BESTCLOSESTD(1:MIN(10,NCFORNEWATOM))
         ENDIF
!
!  Maintain a sorted list of active atoms that are constrained to the new atom, sorted
!  according to their distance.
!
         NCONFORNEWATOM=0
         CONDIST(1:NATOMS)=1.0D100
         IF (DEBUG) WRITE(*,'(3(A,I6))') ' intlbfgs> New active atom is number ',NEWATOM,' total=',NACTIVE+1, &
 &                        ' steps=',NITERDONE
         DO J1=1,NCONSTRAINT
            IF (CONACTIVE(J1)) CYCLE
            IF ((CONI(J1).EQ.NEWATOM).AND.(ATOMACTIVE(CONJ(J1))).OR.(CONJ(J1).EQ.NEWATOM).AND.(ATOMACTIVE(CONI(J1)))) THEN  
                 NCONFORNEWATOM=NCONFORNEWATOM+1
!                CONACTIVE(J1)=.TRUE.
!                NITSTART(J1)=NITERDONE
!                NCONSTRAINTON=NCONSTRAINTON+1
! !
! ! The ...ON variables are not actually used in congrad.f90.
! !
!                CONDISTREFLOCALON(NCONSTRAINTON)=CONDISTREFLOCAL(J1)
!                CONDISTREFON(NCONSTRAINTON)=CONDISTREF(J1)
!                CONION(NCONSTRAINTON)=CONI(J1)
!                CONJON(NCONSTRAINTON)=CONJ(J1)
! 
!                IF (DEBUG) WRITE(*,'(A,I6,A,2I6)') ' intlbfgs> Turning on constraint ',J1,' for atoms ',CONI(J1),CONJ(J1)
               IF (NCONFORNEWATOM.EQ.1) THEN
                  CONDIST(1)=CONDISTREF(J1)
                  IF (CONI(J1).EQ.NEWATOM) CONLIST(1)=CONJ(J1)
                  IF (CONJ(J1).EQ.NEWATOM) CONLIST(1)=CONI(J1)
               ENDIF
               DO J2=1,NCONFORNEWATOM-1
                  IF (CONDISTREF(J1).LT.CONDIST(J2)) THEN
!                    WRITE(*,'(A,I6,G12.4,I6,G12.4)') 'J1,CONDISTREF < J2,CONDIST: ',J1,CONDISTREF(J1),J2,CONDIST(J2)
                     DO J3=NCONFORNEWATOM,J2+1,-1
!                       WRITE(*,'(A,I6,A,I6,A,G12.4)') ' moving dist and list from ',J3-1,' to ',J3,' CONDIST=',CONDIST(J3-1)
                        CONDIST(J3)=CONDIST(J3-1)
                        CONLIST(J3)=CONLIST(J3-1)
                     ENDDO
                     CONDIST(J2)=CONDISTREF(J1)
!                    WRITE(*,'(A,I6,A,G12.4)') ' setting condist element ',J2,' to ',CONDISTREF(J1)
                     IF (CONI(J1).EQ.NEWATOM) CONLIST(J2)=CONJ(J1)
                     IF (CONJ(J1).EQ.NEWATOM) CONLIST(J2)=CONI(J1)
!                    WRITE(*,'(A,I6,A,G12.4)') ' setting conlist element ',J2,' to ',CONLIST(J2)
                     GOTO 654
                  ENDIF
               ENDDO 
               CONDIST(NCONFORNEWATOM)=CONDISTREF(J1)
!              WRITE(*,'(A,I6,A,G12.4)') ' setting condist element ',NCONFORNEWATOM,' to ',CONDISTREF(J1)
               IF (CONI(J1).EQ.NEWATOM) CONLIST(NCONFORNEWATOM)=CONJ(J1)
               IF (CONJ(J1).EQ.NEWATOM) CONLIST(NCONFORNEWATOM)=CONI(J1)
!              WRITE(*,'(A,I6,A,G12.4)') ' setting conlist element ',NCONFORNEWATOM,' to ',CONLIST(NCONFORNEWATOM)
654          CONTINUE
            ENDIF
         ENDDO 
         IF (DEBUG) THEN
            WRITE(*,'(A,I6,A,I6,A)') ' intlbfgs> New active atom ',NEWATOM,' is constrained to ',NCONFORNEWATOM, &
  &                                       ' other active atoms:'
            WRITE(*,'(20I6)') CONLIST(1:NCONFORNEWATOM)
            WRITE(*,'(A,I6,A,I6,A)') ' intlbfgs> sorted distances:'
            WRITE(*,'(10G12.4)') CONDIST(1:NCONFORNEWATOM)
         ENDIF
         DO J1=1,MIN(MAXCONUSE,NCONFORNEWATOM)
            DO J2=1,NCONSTRAINT
               IF ((CONI(J2).EQ.NEWATOM).AND.(CONJ(J2).EQ.CONLIST(J1))) THEN
                     CONACTIVE(J2)=.TRUE.
                     IF (DEBUG) WRITE(*,'(A,I6,A,2I6)') ' intlbfgs> Turning on constraint ',J2,' for atoms ',CONI(J2),CONJ(J2)
               ELSE IF ((CONJ(J2).EQ.NEWATOM).AND.(CONI(J2).EQ.CONLIST(J1))) THEN
                     CONACTIVE(J2)=.TRUE.
                     IF (DEBUG) WRITE(*,'(A,I6,A,2I6)') ' intlbfgs> Turning on constraint ',J2,' for atoms ',CONI(J2),CONJ(J2)
               ENDIF
            ENDDO
         ENDDO

         DO J1=1,NATOMS
            IF (.NOT.ATOMACTIVE(J1)) CYCLE ! identify active atoms
            IF (ABS(J1-NEWATOM).LE.INTREPSEP) CYCLE ! no repulsion for atoms too close in sequence
            DO J2=1,NCONSTRAINT
!
!  With MAXCONUSE set to a finite value there could be constraints for the new atom that are
!  not active. We don't want these to be changed to repulsion, surely?!
!  Or perhaps we do need to do something with them?
!
!              IF (.NOT.CONACTIVE(J2)) CYCLE ! repulsions for inactive constraints 
               IF (((CONI(J2).EQ.J1).AND.(CONJ(J2).EQ.NEWATOM)).OR.((CONJ(J2).EQ.J1).AND.(CONI(J2).EQ.NEWATOM))) GOTO 543
            ENDDO
            DMIN=1.0D100
            DO J2=1,INTIMAGE+2,INTIMAGE+1 ! only consider the end-point distances
               DF=SQRT((XYZ((J2-1)*3*NATOMS+3*(NEWATOM-1)+1)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+1))**2+ &
  &                    (XYZ((J2-1)*3*NATOMS+3*(NEWATOM-1)+2)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+2))**2+ &
  &                    (XYZ((J2-1)*3*NATOMS+3*(NEWATOM-1)+3)-XYZ((J2-1)*3*NATOMS+3*(J1-1)+3))**2)
               IF (DF.LT.DMIN) DMIN=DF
            ENDDO
!
! Use the minimum of the end point distances and INTCONSTRAINREPCUT for each contact.
!
            DMIN=MIN(DMIN-1.0D-3,INTCONSTRAINREPCUT)
            NREPULSIVE=NREPULSIVE+1
            IF (NREPULSIVE.GT.NREPMAX) CALL REPDOUBLE
            REPI(NREPULSIVE)=J1
            REPJ(NREPULSIVE)=NEWATOM
            REPCUT(NREPULSIVE)=DMIN
!           IF (DEBUG) WRITE(*,'(A,I6,A,I6,A,F15.5)') ' intlbfgs> Adding repulsion for new atom ',NEWATOM,' with atom ',J1, &
! &                                                   ' cutoff=',DMIN
543         CONTINUE
         ENDDO
         ATOMACTIVE(NEWATOM)=.TRUE.
         NACTIVE=NACTIVE+1
         IF (MAXNACTIVE.EQ.0) MAXNACTIVE=NATOMS
!
! Freeze atoms that became active more than NACTIVE-MAXNACTIVE events ago.
! For example, with MAXNACTIVE=5 and 40 active atoms, we would freeze those 
! turned on first, second, up to the 35th in the TURNONORDER list.
!
         IF (NACTIVE.GT.MAXNACTIVE) THEN
!           WRITE(*,'(A)') 'doaddatom> TURNONORDER:'
!           WRITE(*,'(5I6)') TURNONORDER(1:NACTIVE-1)
            NDUMMY=TURNONORDER(NACTIVE-MAXNACTIVE)
            IF (INTFROZEN(NDUMMY)) THEN
               IF (DEBUG) WRITE(*,'(A,I6,A,2I6)') ' doaddatom> Not turning off frozen active atom ',NDUMMY,' already frozen'
            ELSE
               IF (DEBUG) WRITE(*,'(A,I6,A,2I6)') ' doaddatom> Freezing active atom ',NDUMMY
               INTFROZEN(NDUMMY)=.TRUE.
!
! Turn off constraints and repulsions between frozen atoms.
!
               DO J2=1,NCONSTRAINT
                  IF (.NOT.CONACTIVE(J2)) CYCLE
                  IF (INTFROZEN(CONI(J2)).AND.INTFROZEN(CONJ(J2))) THEN
                     CONACTIVE(J2)=.FALSE.
                     WRITE(*,'(A,I6,A,2I6)') 'doaddatom> turning off constraint ',J2,' between atoms ',CONI(J2),CONJ(J2)
                  ENDIF
               ENDDO

               J2=0
               DO J1=1,NREPULSIVEFIX
                  IF (INTFROZEN(REPIFIX(J1)).AND.INTFROZEN(REPJFIX(J1))) CYCLE
                  IF (ATOMACTIVE(REPIFIX(J1)).AND.ATOMACTIVE(REPJFIX(J1))) THEN
                     DO J3=1,NCONSTRAINTFIX
!                       IF (.NOT.CONACTIVE(J3)) CYCLE ! no repulsions for any constraints
                        IF ((CONIFIX(J3).EQ.REPIFIX(J1)).AND.(CONJFIX(J3).EQ.REPJFIX(J1))) GOTO 962
                        IF ((CONIFIX(J3).EQ.REPJFIX(J1)).AND.(CONJFIX(J3).EQ.REPIFIX(J1))) GOTO 962
                     ENDDO
                     J2=J2+1
                     REPI(J2)=REPIFIX(J1)
                     REPJ(J2)=REPJFIX(J1)
                     REPCUT(J2)=REPCUTFIX(J1)
962                  CONTINUE
                  ENDIF
               ENDDO
               NREPULSIVE=J2
               WRITE(*,'(A,I6,A)') ' doaddatom> After allowing for frozen atoms there are ',NREPULSIVE,' possible repulsions'
               NREPI(1:NREPULSIVE)=REPI(1:NREPULSIVE)
               NREPJ(1:NREPULSIVE)=REPJ(1:NREPULSIVE)
               NNREPULSIVE=NREPULSIVE
               NREPCUT(1:NREPULSIVE)=REPCUT(1:NREPULSIVE)
            ENDIF
         ENDIF

         NDUMMY=0
         DO J1=1,NATOMS
            IF (ATOMACTIVE(J1)) NDUMMY=NDUMMY+1
         ENDDO
         IF (NDUMMY.NE.NACTIVE) THEN
            WRITE(*,'(A,I6)') ' doaddatom> ERROR *** inconsistency in number of active atoms. ',NDUMMY,' should be ',NACTIVE
            DO J1=1,NATOMS
               IF (ATOMACTIVE(J1)) WRITE(*,'(A,I6)') ' active atom ',J1
            ENDDO
            STOP
         ENDIF

         TURNONORDER(NACTIVE)=NEWATOM
!
! Initial guess for new active atom position. This is crucial for success in INTCONSTRAINT schemes!
!
         ESAVED=1.0D100
         ESAVE0=1.0D100
         ESAVEC=1.0D100
         FTEST=.TRUE.
         IDONE=.FALSE.
         IF (NCONFORNEWATOM.GE.3) THEN
            IDONE=.TRUE.
!
! Move the new atom consistently in the local environment of its three nearest actively constrained atoms.
! Make a local orthogonal coordinate system and use constant components in this basis.
!
            N1=NCONFORNEWATOM-2; N2=NCONFORNEWATOM-1; N3=NCONFORNEWATOM
!           N1=1; N2=2; N3=3
            IF (DEBUG) WRITE(*,'(A,3I6)') ' intlbfgs> initial guess from furthest three constrained active atoms, ',CONLIST(N1:N3)
            VEC1(1:3)=XYZ(3*(CONLIST(N2)-1)+1:3*(CONLIST(N2)-1)+3)-XYZ(3*(CONLIST(N1)-1)+1:3*(CONLIST(N1)-1)+3)
            DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
            IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
            VEC2(1:3)=XYZ(3*(CONLIST(N3)-1)+1:3*(CONLIST(N3)-1)+3)-XYZ(3*(CONLIST(N1)-1)+1:3*(CONLIST(N1)-1)+3)
            DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
            VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
            DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
            IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
            VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
            VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
            VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
            C1=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(CONLIST(N1)-1)+1))*VEC1(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(CONLIST(N1)-1)+2))*VEC1(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(CONLIST(N1)-1)+3))*VEC1(3)
            C2=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(CONLIST(N1)-1)+1))*VEC2(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(CONLIST(N1)-1)+2))*VEC2(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(CONLIST(N1)-1)+3))*VEC2(3)
            C3=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(CONLIST(N1)-1)+1))*VEC3(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(CONLIST(N1)-1)+2))*VEC3(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(CONLIST(N1)-1)+3))*VEC3(3)

            DO J1=2,INTIMAGE+1
               VEC1(1:3)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(N2)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(N2)-1)+3) &
  &                     -XYZ((J1-1)*3*NATOMS+3*(CONLIST(N1)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(N1)-1)+3)
               DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
               IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
               VEC2(1:3)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(N3)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(N3)-1)+3) &
  &                     -XYZ((J1-1)*3*NATOMS+3*(CONLIST(N1)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(N1)-1)+3)
               DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
               VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
               DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
               IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
               VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
               VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
               VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)= &
  &               XYZ((J1-1)*3*NATOMS+3*(CONLIST(N1)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(N1)-1)+3)+C1*VEC1(1:3)+C2*VEC2(1:3)+C3*VEC3(1:3)

!
! Alternative analytical solution from intersection of three spheres by trilateration
!
               IF (QCITRILAT) THEN
                  P1(1:3)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(N1)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(N1)-1)+3)
                  P2(1:3)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(N2)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(N2)-1)+3)
                  P3(1:3)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(N3)-1)+1:(J1-1)*3*NATOMS+3*(CONLIST(N3)-1)+3)
                  R1=CONDIST(N1)
                  R2=CONDIST(N2)
                  R3=CONDIST(N3)
                  CALL TRILATERATION(P1,P2,P3,R1,R2,R3,SOL1,SOL2,FTEST)
                  IF (FTEST) THEN
!                    WRITE(*,'(A,I8)') ' intlbfgs> WARNING *** no trilateration solution for image ',J1
                  ELSE
!                    WRITE(*,'(A,I8)') ' intlbfgs>                trilateration solution for image ',J1
!                    WRITE(*,'(A,3F20.10)') ' intlbfgs> SOL1=',SOL1(1:3)
!                    WRITE(*,'(A,3F20.10)') ' intlbfgs> SOL2=',SOL2(1:3)
!                    WRITE(*,'(A,3F20.10)') ' intlbfgs> prev=',XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)
!                    D1SQ=(SOL1(1)-XYZ((J1-2)*3*NATOMS+3*(NEWATOM-1)+1))**2 &
! &                      +(SOL1(2)-XYZ((J1-2)*3*NATOMS+3*(NEWATOM-1)+2))**2 &
! &                      +(SOL1(3)-XYZ((J1-2)*3*NATOMS+3*(NEWATOM-1)+3))**2
!                    D2SQ=(SOL2(1)-XYZ((J1-2)*3*NATOMS+3*(NEWATOM-1)+1))**2 &
! &                      +(SOL2(2)-XYZ((J1-2)*3*NATOMS+3*(NEWATOM-1)+2))**2 &
! &                      +(SOL2(3)-XYZ((J1-2)*3*NATOMS+3*(NEWATOM-1)+3))**2
!
! Try minmum distance from previous solution
!
                     D1SQ=(SOL1(1)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1))**2 &
  &                      +(SOL1(2)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+2))**2 &
  &                      +(SOL1(3)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+3))**2
                     D2SQ=(SOL2(1)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1))**2 &
  &                      +(SOL2(2)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+2))**2 &
  &                      +(SOL2(3)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+3))**2
!                    WRITE(*,'(A,2F20.10)') 'D1SQ,D2SQ=',D1SQ,D2SQ
                     IF (D1SQ.LT.D2SQ) THEN
                        XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=SOL1(1:3)
                     ELSE
                        XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=SOL2(1:3)
                     ENDIF
                  ENDIF
               ENDIF

            ENDDO
            CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),NNREPSAVE,NREPSAVE+1) ! set up repulsive neighbour list
            IF (QCIADDREP.GT.0) THEN
               CALL CONGRAD3(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSEIF (CHECKCONINT) THEN
               CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSE
               CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ENDIF
            ESAVE0=ETOTAL
            DO J1=2,INTIMAGE+1
               XSAVE0(1:3,J1)=XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)
            ENDDO
         ENDIF

!          IF ((.NOT.IDONE).AND.(NDFORNEWATOM.GE.3)) THEN
!             IDONE=.TRUE.
! !
! ! Choose three atoms from the BESTPRESERVEDN list at random with bias towards the 
! ! start of the list. Let the relative weight for position i be 1/i**2 and calculate
! ! the sum to normalise.
! !
!             DUMMY=0.0D0
!             DO J1=1,NDFORNEWATOM
! !              DUMMY=DUMMY+1.0D0/(1.0D0*J1)
! !              DUMMY=DUMMY+1.0D0/(1.0D0*BESTPRESERVEDD(J1))
!                DUMMY=DUMMY+1.0D0/(1.0D0*J1**2)
!             ENDDO
!             N1=0; N2=0; N3=0
!             DO WHILE (N3.EQ.0)
!                DUMMY2=0.0D0
!                RAN1=DPRAND()*DUMMY
!                DO J1=1,NDFORNEWATOM
! !                 DUMMY2=DUMMY2+1.0D0/(1.0D0*J1)
! !                 DUMMY2=DUMMY2+1.0D0/(1.0D0*BESTPRESERVEDD(J1))
!                   DUMMY2=DUMMY2+1.0D0/(1.0D0*J1**2)
!                   IF (DUMMY2.GE.RAN1) THEN
!                      IF ((J1.EQ.N1).OR.(J1.EQ.N2)) EXIT ! already chosen
!                      IF (N1.EQ.0) THEN
!                         N1=J1
!                         EXIT
!                      ENDIF
!                      IF (N2.EQ.0) THEN
!                         N2=J1
!                         EXIT
!                      ENDIF
!                      N3=J1
!                      EXIT
!                   ENDIF
!                ENDDO
!             ENDDO
!             IF (DEBUG) WRITE(*,'(A,3I6,A)') ' intlbfgs> choosing positions ',N1,N2,N3,' in best preserved list'
!             IF (DEBUG) WRITE(*,'(A,3I6)') ' intlbfgs> atoms are ',BESTPRESERVEDN(N1),BESTPRESERVEDN(N2),BESTPRESERVEDN(N3)
! !           IF (DEBUG) WRITE(*,'(A,3I6,A)') ' intlbfgs> full list has length ',NDFORNEWATOM
! !           IF (DEBUG) WRITE(*,'(20I6)') BESTPRESERVEDN(1:NDFORNEWATOM)
! 
! !
! ! Move the new atom consistently in the local environment of the three active atoms with the
! ! best preserved absolute distances or the shortest average distances in the end points.
! ! Check the energies and compare linear interpolation as well, then choose the interpolation
! ! with the lowest energy.
! ! Make a local orthogonal coordinate system and use constant components in this basis.
! !
!             VEC1(1:3)=XYZ(3*(BESTPRESERVEDN(N2)-1)+1:3*(BESTPRESERVEDN(N2)-1)+3) &
!   &                  -XYZ(3*(BESTPRESERVEDN(N1)-1)+1:3*(BESTPRESERVEDN(N1)-1)+3)
!             DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
!             IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
!             VEC2(1:3)=XYZ(3*(BESTPRESERVEDN(N3)-1)+1:3*(BESTPRESERVEDN(N3)-1)+3) &
!   &                  -XYZ(3*(BESTPRESERVEDN(N1)-1)+1:3*(BESTPRESERVEDN(N1)-1)+3)
!             DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
!             VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
!             DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
!             IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
!             VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
!             VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
!             VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
!             C1=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTPRESERVEDN(N1)-1)+1))*VEC1(1)+ &
!   &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTPRESERVEDN(N1)-1)+2))*VEC1(2)+ &
!   &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTPRESERVEDN(N1)-1)+3))*VEC1(3)
!             C2=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTPRESERVEDN(N1)-1)+1))*VEC2(1)+ &
!   &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTPRESERVEDN(N1)-1)+2))*VEC2(2)+ &
!   &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTPRESERVEDN(N1)-1)+3))*VEC2(3)
!             C3=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTPRESERVEDN(N1)-1)+1))*VEC3(1)+ &
!   &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTPRESERVEDN(N1)-1)+2))*VEC3(2)+ &
!   &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTPRESERVEDN(N1)-1)+3))*VEC3(3)
!             DO J1=2,INTIMAGE+1
!                VEC1(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N2)-1)+1:(J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N2)-1)+3) &
!   &                     -XYZ((J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N1)-1)+3)
!                DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
!                IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
!                VEC2(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N3)-1)+1:(J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N3)-1)+3) &
!   &                     -XYZ((J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N1)-1)+3)
!                DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
!                VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
!                DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
!                IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
!                VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
!                VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
!                VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
!                XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)= &
!   &            XYZ((J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTPRESERVEDN(N1)-1)+3)+ &
!   &                   C1*VEC1(1:3)+C2*VEC2(1:3)+C3*VEC3(1:3)+0.01D0*(DPRAND()-0.5D0)*2.0D0
! !              WRITE(*,'(A,I6,3G20.10)') 'intlbfgs> J1,C1,C2,C3=',J1,C1,C2,C3
! !              WRITE(*,'(A,9G20.10)') 'intlbfgs> VEC1,2,3=',VEC1(1:3),VEC2(1:3),VEC3(1:3)
! !              WRITE(*,'(A,6I6)') 'intlbfgs> N1,N2,N3,Bestpreserved N1,N2,N3=',N1,N2,N3, &
! ! &                 BESTPRESERVEDN(N1),BESTPRESERVEDN(N2),BESTPRESERVEDN(N3)
! 
! !
! ! Alternative analytical solution from intersection of three spheres by trilateration not available - we haven't saved the distances,
! ! just the differences in the endpoints
! !
!             ENDDO
! 
!             CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),NNREPSAVE,NREPSAVE+1) ! set up repulsive neighbour list
!             IF (QCIADDREP.GT.0) THEN
!                CALL CONGRAD3(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
!             ELSEIF (CHECKCONINT) THEN
!                CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
!             ELSE
!                CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
!             ENDIF
!             ESAVED=ETOTAL
!             DO J1=2,INTIMAGE+1
!                XSAVED(1:3,J1)=XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)
!             ENDDO
!          ENDIF

         IF ((.NOT.IDONE).AND.(NCFORNEWATOM.GE.3)) THEN
            IDONE=.TRUE.
!
! Choose three atoms from the BESTCLOSEST list at random with bias towards the
! start of the list. Let the relative weight for position i be 1/i**2 and calculate
! the sum to normalise.
!
            DUMMY=0.0D0
!             DO J1=1,NCFORNEWATOM
! !              DUMMY=DUMMY+1.0D0/(1.0D0*J1)
! !              DUMMY=DUMMY+1.0D0/(1.0D0*BESTCLOSESTD(J1))
!                DUMMY=DUMMY+1.0D0/(1.0D0*J1**2)
!             ENDDO
!             N1=0; N2=0; N3=0
!             DO WHILE (N3.EQ.0)
!                DUMMY2=0.0D0
!                RAN1=DPRAND()*DUMMY
!                DO J1=1,NCFORNEWATOM
! !                 DUMMY2=DUMMY2+1.0D0/(1.0D0*J1)
! !                 DUMMY2=DUMMY2+1.0D0/(1.0D0*BESTCLOSESTD(J1))
!                   DUMMY2=DUMMY2+1.0D0/(1.0D0*J1**2)
!                   IF (DUMMY2.GE.RAN1) THEN
!                      IF ((J1.EQ.N1).OR.(J1.EQ.N2)) EXIT ! already chosen
!                      IF (N1.EQ.0) THEN
!                         N1=J1
!                         EXIT
!                      ENDIF
!                      IF (N2.EQ.0) THEN
!                         N2=J1
!                         EXIT
!                      ENDIF
!                      N3=J1
!                      EXIT
!                   ENDIF
!                ENDDO
!             ENDDO
            N1=1; N2=2; N3=3
            IF (DEBUG) WRITE(*,'(A,3I6,A)') ' intlbfgs> choosing positions ',N1,N2,N3,' in closest list'

            VEC1(1:3)=XYZ(3*(BESTCLOSESTN(N2)-1)+1:3*(BESTCLOSESTN(N2)-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+1:3*(BESTCLOSESTN(N1)-1)+3)
            DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
            IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
            VEC2(1:3)=XYZ(3*(BESTCLOSESTN(N3)-1)+1:3*(BESTCLOSESTN(N3)-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+1:3*(BESTCLOSESTN(N1)-1)+3)
            DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
            VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
            DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
            IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
            VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
            VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
            VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
            C1=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTCLOSESTN(N1)-1)+1))*VEC1(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTCLOSESTN(N1)-1)+2))*VEC1(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+3))*VEC1(3)
            C2=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTCLOSESTN(N1)-1)+1))*VEC2(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTCLOSESTN(N1)-1)+2))*VEC2(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+3))*VEC2(3)
            C3=(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(BESTCLOSESTN(N1)-1)+1))*VEC3(1)+ &
  &            (XYZ(3*(NEWATOM-1)+2)-XYZ(3*(BESTCLOSESTN(N1)-1)+2))*VEC3(2)+ &
  &            (XYZ(3*(NEWATOM-1)+3)-XYZ(3*(BESTCLOSESTN(N1)-1)+3))*VEC3(3)
            DO J1=2,INTIMAGE+1
               VEC1(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N2)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N2)-1)+3) &
  &                     -XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+3)
               DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
               IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)/DUMMY
               VEC2(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N3)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N3)-1)+3) &
  &                     -XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+3)
               DUMMY=VEC1(1)*VEC2(1)+VEC1(2)*VEC2(2)+VEC1(3)*VEC2(3)
               VEC2(1:3)=VEC2(1:3)-DUMMY*VEC1(1:3)
               DUMMY=SQRT(VEC2(1)**2+VEC2(2)**2+VEC2(3)**2)
               IF (DUMMY.NE.0.0D0) VEC2(1:3)=VEC2(1:3)/DUMMY
               VEC3(1)= VEC1(2)*VEC2(3)-VEC1(3)*VEC2(2)
               VEC3(2)=-VEC1(1)*VEC2(3)+VEC1(3)*VEC2(1)
               VEC3(3)= VEC1(1)*VEC2(2)-VEC1(2)*VEC2(1)
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)= &
  &            XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+3)+ &
  &                   C1*VEC1(1:3)+C2*VEC2(1:3)+C3*VEC3(1:3)

!
! Alternative analytical solution from intersection of three spheres by trilateration
!
               IF (QCITRILAT) THEN
                  P1(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N1)-1)+3)
                  P2(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N2)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N2)-1)+3)
                  P3(1:3)=XYZ((J1-1)*3*NATOMS+3*(BESTCLOSESTN(N3)-1)+1:(J1-1)*3*NATOMS+3*(BESTCLOSESTN(N3)-1)+3)
                  R1=BESTCLOSESTD(N1)
                  R2=BESTCLOSESTD(N2)
                  R3=BESTCLOSESTD(N3)
                  CALL TRILATERATION(P1,P2,P3,R1,R2,R3,SOL1,SOL2,FTEST)
                  IF (FTEST) THEN
!                 WRITE(*,'(A,I8)') ' intlbfgs> WARNING *** no trilateration solution for image ',J1
                  ELSE
!                 WRITE(*,'(A,I8)') ' intlbfgs>                trilateration solution for image ',J1
!                 WRITE(*,'(A,3F20.10)') ' intlbfgs> SOL1=',SOL1(1:3)
!                 WRITE(*,'(A,3F20.10)') ' intlbfgs> SOL2=',SOL2(1:3)
!                 WRITE(*,'(A,3F20.10)') ' intlbfgs> prev=',XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)
!                    D1SQ=(SOL1(1)-XYZ((J1-2)*3*NATOMS+3*(NEWATOM-1)+1))**2 &
!    &                   +(SOL1(2)-XYZ((J1-2)*3*NATOMS+3*(NEWATOM-1)+2))**2 &
!    &                   +(SOL1(3)-XYZ((J1-2)*3*NATOMS+3*(NEWATOM-1)+3))**2
!                    D2SQ=(SOL2(1)-XYZ((J1-2)*3*NATOMS+3*(NEWATOM-1)+1))**2 &
!    &                   +(SOL2(2)-XYZ((J1-2)*3*NATOMS+3*(NEWATOM-1)+2))**2 &
!    &                   +(SOL2(3)-XYZ((J1-2)*3*NATOMS+3*(NEWATOM-1)+3))**2
!
! Try minmum distance from previous solution
!
                     D1SQ=(SOL1(1)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1))**2 &
  &                      +(SOL1(2)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+2))**2 &
  &                      +(SOL1(3)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+3))**2
                     D2SQ=(SOL2(1)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1))**2 &
  &                      +(SOL2(2)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+2))**2 &
  &                      +(SOL2(3)-XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+3))**2
!                 WRITE(*,'(A,2F20.10)') 'D1SQ,D2SQ=',D1SQ,D2SQ
                     IF (D1SQ.LT.D2SQ) THEN
                        XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=SOL1(1:3)
                     ELSE
                        XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=SOL2(1:3)
                     ENDIF
                  ENDIF
               ENDIF
            ENDDO

            CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),NNREPSAVE,NREPSAVE+1) ! set up repulsive neighbour list
            IF (QCIADDREP.GT.0) THEN
               CALL CONGRAD3(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSEIF (CHECKCONINT) THEN
               CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSE
               CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ENDIF
            ESAVEC=ETOTAL
            DO J1=2,INTIMAGE+1
               XSAVEC(1:3,J1)=XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)
            ENDDO
         ENDIF
!
! Standard linear interpolation, with constraint distance scaled by FRAC.
! Works for FRAC as small as 0.1 with repulsion turned off.
! We use an appropriately weighted displacement from atom CONLIST(1) using the displacements
! in the two end points.
!
         ETOTAL=1.0D100
         IF (.NOT.IDONE) THEN
            FRAC=1.0D0
            DO J1=2,INTIMAGE+1
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(1)-1)+1)  &
 &            +(INTIMAGE-J1+2)*FRAC*(XYZ(3*(NEWATOM-1)+1)-XYZ(3*(CONLIST(1)-1)+1))/(INTIMAGE+1) &
 &   +(J1-1)*(XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+1)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(CONLIST(1)-1)+1))/(INTIMAGE+1)
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+2)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(1)-1)+2)  &
 &            +(INTIMAGE-J1+2)*FRAC*(XYZ(3*(NEWATOM-1)+2)-XYZ(3*(CONLIST(1)-1)+2))/(INTIMAGE+1) &
 &   +(J1-1)*(XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+2)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(CONLIST(1)-1)+2))/(INTIMAGE+1)
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=XYZ((J1-1)*3*NATOMS+3*(CONLIST(1)-1)+3)  &
 &            +(INTIMAGE-J1+2)*FRAC*(XYZ(3*(NEWATOM-1)+3)-XYZ(3*(CONLIST(1)-1)+3))/(INTIMAGE+1) &
 &   +(J1-1)*(XYZ(3*NATOMS*(INTIMAGE+1)+3*(NEWATOM-1)+3)-XYZ(3*NATOMS*(INTIMAGE+1)+3*(CONLIST(1)-1)+3))/(INTIMAGE+1)
            ENDDO
            CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),NNREPSAVE,NREPSAVE+1) ! set up repulsive neighbour list
            IF (QCIADDREP.GT.0) THEN
               CALL CONGRAD3(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSEIF (CHECKCONINT) THEN
               CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ELSE
               CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
            ENDIF
         ENDIF

         IF (DEBUG) WRITE(*,'(A,4G15.5)') ' intlbfgs> energies for constrained, preserved, closest, and linear schemes=', &
  &              ESAVE0,ESAVED,ESAVEC,ETOTAL
    
         IF ((ETOTAL.LT.ESAVEC).AND.(ETOTAL.LT.ESAVED).AND.(ETOTAL.LT.ESAVE0)) THEN
            IF (DEBUG) WRITE(*,'(A,2G20.10)') ' intlbfgs> lowest energy from linear interpolation'
         ELSE IF ((ESAVEC.LT.ESAVED).AND.(ESAVEC.LT.ESAVE0)) THEN
            IF (DEBUG) WRITE(*,'(A,2G20.10)') ' intlbfgs> lowest energy from interpolation using closest atoms'
            DO J1=2,INTIMAGE+1
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=XSAVEC(1:3,J1)
            ENDDO
            ETOTAL=ESAVEC
         ELSE IF (ESAVED.LT.ESAVE0) THEN
            IF (DEBUG) WRITE(*,'(A,2G20.10)') ' intlbfgs> lowest energy from interpolation using preserved distances'
            DO J1=2,INTIMAGE+1
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=XSAVED(1:3,J1)
            ENDDO
            ETOTAL=ESAVED
         ELSE 
            IF (DEBUG) WRITE(*,'(A,2G20.10)') ' intlbfgs> interpolation using closest constraints'
            DO J1=2,INTIMAGE+1
               XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)=XSAVE0(1:3,J1)
            ENDDO
            ETOTAL=ESAVE0
         ENDIF
      ENDIF
      NADDED=NADDED+1
      IF (NADDED.LT.NTOADD) GOTO 542
!
! Check whether we've added all atoms in the amino acid corresponding to the new atom. If not, go back to the top
! and choose the next candidate.
!
      IF (QCIADDACIDT.AND.(.NOT.DOBACK)) THEN
         DO J1=1,NATOMS
            IF ((ATOMSTORES(J1).EQ.ACID).AND.(.NOT.(ATOMACTIVE(J1)))) GOTO 542
         ENDDO
         WRITE(*,'(A,I6,A)') 'doaddatom> All atoms of residue ',ACID,' are active'
      ENDIF
      IF (DOBACKALL.AND.(AABACK(NEWATOM))) THEN
         DO J1=1,NATOMS
            IF ((ATOMSTORES(J1).EQ.ACID).AND.(.NOT.(ATOMACTIVE(J1))).AND.AABACK(J1)) GOTO 542
         ENDDO
         WRITE(*,'(A,I6,A)') 'doaddatom> All backbone atoms of residue ',ACID,' are active'
      ENDIF

      IF (QCIRADSHIFTT) THEN
         WRITE(*,'(A,F15.5)') ' intlbfgs> Applying radial shift for unconstrained atoms of ',QCIRADSHIFT
         WRITE(*,'(20I6)') CONLIST(1:NCONFORNEWATOM)
         DO J1=2,INTIMAGE+1
            scaleloop: DO J2=1,NATOMS
               IF (.NOT.ATOMACTIVE(J2)) CYCLE scaleloop
               IF (J2.EQ.NEWATOM) CYCLE scaleloop
               DO J3=1,NCONFORNEWATOM
                  IF (CONLIST(J3).EQ.J2) CYCLE scaleloop
               ENDDO
               VEC1(1:3)=XYZ((J1-1)*3*NATOMS+3*(J2-1)+1:(J1-1)*3*NATOMS+3*(J2-1)+3)- &
   &                     XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)
               DUMMY=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
               IF (DUMMY.NE.0.0D0) VEC1(1:3)=VEC1(1:3)*QCIRADSHIFT/DUMMY
               XYZ((J1-1)*3*NATOMS+3*(J2-1)+1:(J1-1)*3*NATOMS+3*(J2-1)+3)= &
   &           XYZ((J1-1)*3*NATOMS+3*(J2-1)+1:(J1-1)*3*NATOMS+3*(J2-1)+3)+VEC1(1:3)
!!!!!!!!!!! debug DJW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!              VEC1(1:3)=XYZ((J1-1)*3*NATOMS+3*(J2-1)+1:(J1-1)*3*NATOMS+3*(J2-1)+3)- &
!  &                     XYZ((J1-1)*3*NATOMS+3*(NEWATOM-1)+1:(J1-1)*3*NATOMS+3*(NEWATOM-1)+3)
!              DUMMY2=SQRT(VEC1(1)**2+VEC1(2)**2+VEC1(3)**2)
!              PRINT '(A,I6,A,2I6,A,2F15.5)','image ',J1,' atoms ',NEWATOM,J2,' initial and final distance=',DUMMY,DUMMY2
!!!!!!!!!!! debug DJW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ENDDO scaleloop
         ENDDO
      ENDIF
!
! Turn frozen images off for new added atom.
!
!     IF (DEBUG) WRITE(*,'(A)') ' intlbfgs> turning off frozen images'
!     IF (FREEZENODEST) IMGFREEZE(1:INTIMAGE)=.FALSE.
      CALL CHECKREP(INTIMAGE,XYZ,(3*NATOMS),NNREPSAVE,NREPSAVE+1) ! set up repulsive neighbour list
!
! need a new gradient since the active atom has changed !
!
      IF (QCIADDREP.GT.0) THEN
         CALL CONGRAD3(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ELSEIF (CHECKCONINT) THEN
         CALL CONGRAD2(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ELSE
         CALL CONGRAD(NMAXINT,NMININT,ETOTAL,XYZ,GGG,EEE,IMGFREEZE,RMS)
      ENDIF

END SUBROUTINE DOADDATOM

SUBROUTINE CHECKPERC(LXYZ,LINTCONSTRAINTTOL,NQCIFREEZE,NCPFIT)
USE KEY, ONLY : ATOMACTIVE, NCONSTRAINT, INTFROZEN, CONI, CONJ, CONDISTREF, INTCONMAX, INTCONSTRAINTTOL, &
  &             INTCONSEP, NCONGEOM, CONGEOM, CONIFIX, CONJFIX, CONDISTREFFIX, INTCONCUT, &
  &             NCONSTRAINTFIX, BULKT, TWOD, RIGIDBODY, CONDATT, CONCUT, CONCUTFIX, &
  &             BONDS, QCIAMBERT, QCIADDREP, QCIADDREPCUT, QCIBONDS, QCISECOND
USE COMMONS, ONLY: NATOMS, DEBUG, PARAM1, PARAM2, PARAM3
IMPLICIT NONE
INTEGER NDIST1(NATOMS), NCYCLE, DMIN1, DMAX1, NUNCON1, J1, J2, J3, NQCIFREEZE, J4, NCPFIT, LUNIT, GETUNIT
INTEGER NI1, NJ1, NI2, NJ2, J5, ATOM1, ATOM2, ACID
DOUBLE PRECISION LINTCONSTRAINTTOL, MAXCONDIST, MINCONDIST, DS, DF, LXYZ((3*NATOMS)*2)
DOUBLE PRECISION DSMIN, DSMAX, DSMEAN, D, DIST2, RMAT(3,3), DUMMY, X1, Y1, Z1, X2, Y2, Z2, DMIN, D2
LOGICAL CHANGED, LDEBUG, CONFILET
LOGICAL :: CALLED=.FALSE.
SAVE CALLED
!for QCIAMBER
INTEGER NBOND, NDUMMY

LINTCONSTRAINTTOL=INTCONSTRAINTTOL

IF (.NOT.ALLOCATED(ATOMACTIVE)) ALLOCATE(ATOMACTIVE(NATOMS))
!
! Fixed constraints based on congeom file entries
! Just need to adjust the list based on any frozen atoms. We
! want to exclude any constraints between two frozen atoms 
! from the list, because subsequent code depends on this.
!

IF (NCONGEOM.GE.2) THEN
   IF (CALLED.OR.CONDATT) THEN
      J2=0
      DO J1=1,NCONSTRAINTFIX
!
! If called with two minima check that CONCUTFIX is large enough to
! accommodate the separation of the two atoms in both minima.
!
         IF (NCPFIT.EQ.2) THEN
            DF=MAX(ABS(CONDISTREFFIX(J1)- &
  &                SQRT((LXYZ(3*(CONIFIX(J1)-1)+1)-LXYZ(3*(CONJFIX(J1)-1)+1))**2+ &
  &                     (LXYZ(3*(CONIFIX(J1)-1)+2)-LXYZ(3*(CONJFIX(J1)-1)+2))**2+ &
  &                     (LXYZ(3*(CONIFIX(J1)-1)+3)-LXYZ(3*(CONJFIX(J1)-1)+3))**2)),&
                   ABS(CONDISTREFFIX(J1)- &
  &                SQRT((LXYZ((3*NATOMS)+3*(CONIFIX(J1)-1)+1)-LXYZ((3*NATOMS)+3*(CONJFIX(J1)-1)+1))**2+ &
  &                     (LXYZ((3*NATOMS)+3*(CONIFIX(J1)-1)+2)-LXYZ((3*NATOMS)+3*(CONJFIX(J1)-1)+2))**2+ &
  &                     (LXYZ((3*NATOMS)+3*(CONIFIX(J1)-1)+3)-LXYZ((3*NATOMS)+3*(CONJFIX(J1)-1)+3))**2)))
            IF (DF.GT.CONCUTFIX(J1)) THEN
               IF (ABS(DF-CONCUTFIX(J1)).GT.1.0D-6) &
  &                WRITE(*,'(A,2I5,3(A,G15.5))') ' checkperc> Increasing con cutoff atoms ', &
  &                CONIFIX(J1),CONJFIX(J1),' from ',CONCUTFIX(J1),' to ',DF,' ref=',CONDISTREFFIX(J1)
               CONCUTFIX(J1)=DF
            ENDIF
         ENDIF
         IF (INTFROZEN(CONIFIX(J1)).AND.INTFROZEN(CONJFIX(J1))) CYCLE
         J2=J2+1
         CONI(J2)=CONIFIX(J1)
         CONJ(J2)=CONJFIX(J1)
         CONDISTREF(J2)=CONDISTREFFIX(J1)
         CONCUT(J2)=CONCUTFIX(J1)
      ENDDO
      NCONSTRAINT=J2
      WRITE(*,'(A,I6,A)') ' checkperc> After allowing for frozen atoms there are ',NCONSTRAINT,' constraints'
      RETURN 
   ELSE
!
! Put reference minima in optimal permutational alignment with reference minimum one.
!
      DO J2=2,NCONGEOM
         LDEBUG=.FALSE.
         CALL MINPERMDIST(CONGEOM(1,1:3*NATOMS),CONGEOM(J2,1:3*NATOMS),NATOMS,LDEBUG, &
  &                       PARAM1,PARAM2,PARAM3,BULKT,TWOD,D,DIST2,RIGIDBODY,RMAT)
      ENDDO
   ENDIF
   ALLOCATE(CONIFIX(INTCONMAX),CONJFIX(INTCONMAX),CONCUTFIX(INTCONMAX),CONDISTREFFIX(INTCONMAX))
ENDIF

INQUIRE(FILE='constraintfile',EXIST=CONFILET)

51   NCONSTRAINT=0 
MAXCONDIST=-1.0D0
MINCONDIST=1.0D100
IF (QCIAMBERT) THEN             
   CALL TOPOLOGY_READER(NBOND)  
!
!  kr366> assume we use two endpoints and topology for amber constraints
!  get number of bonds and bonds from topology
!  loop through all bonds and add them to constraint list
!
   DO J2=1,NBOND                !loop through all bonds and add them to constraint list
      IF (INTFROZEN(BONDS(J2,1)).AND.INTFROZEN(BONDS(J2,2))) CYCLE ! no constraints between intfrozen atoms
      NCONSTRAINT=NCONSTRAINT+1
      IF (DEBUG) WRITE(*,'(A,2I6,A,I6)') 'intlbfgs> Adding constraint for atoms ',BONDS(J2,1),BONDS(J2,2), &
  &                     '  total=',NCONSTRAINT
      DS=SQRT((LXYZ(3*(BONDS(J2,1)-1)+1)-LXYZ(3*(BONDS(J2,2)-1)+1))**2 &
  &          +(LXYZ(3*(BONDS(J2,1)-1)+2)-LXYZ(3*(BONDS(J2,2)-1)+2))**2 &
  &          +(LXYZ(3*(BONDS(J2,1)-1)+3)-LXYZ(3*(BONDS(J2,2)-1)+3))**2)
      DF=SQRT((LXYZ(3*NATOMS+3*(BONDS(J2,1)-1)+1)-LXYZ(3*NATOMS+3*(BONDS(J2,2)-1)+1))**2 &
  &          +(LXYZ(3*NATOMS+3*(BONDS(J2,1)-1)+2)-LXYZ(3*NATOMS+3*(BONDS(J2,2)-1)+2))**2 &
  &          +(LXYZ(3*NATOMS+3*(BONDS(J2,1)-1)+3)-LXYZ(3*NATOMS+3*(BONDS(J2,2)-1)+3))**2)
      IF (NCONSTRAINT.GT.INTCONMAX) CALL CONDOUBLE
      CONI(NCONSTRAINT)=MIN(BONDS(J2,1),BONDS(J2,2))
      CONJ(NCONSTRAINT)=MAX(BONDS(J2,2),BONDS(J2,2))
      CONDISTREF(NCONSTRAINT)=(DF+DS)/2.0D0
      CONCUT(NCONSTRAINT)=ABS(DF-DS)/2.0D0
      IF (CONDISTREF(NCONSTRAINT).GT.MAXCONDIST) MAXCONDIST=CONDISTREF(NCONSTRAINT)
      IF (CONDISTREF(NCONSTRAINT).LT.MINCONDIST) MINCONDIST=CONDISTREF(NCONSTRAINT)
!     IF (DEBUG) WRITE(*,'(A,2I6,A,2F12.2,A,F12.4,A,I8)') ' intlbfgs> constrain distance for ',CONI(NCONSTRAINT), &
! &             CONJ(NCONSTRAINT),' values are ',DS,DF,' fraction=',2*ABS(DS-DF)/(DS+DF), &
! &            ' # bond constraints=',NCONSTRAINT
   ENDDO
   QCIBONDS=NCONSTRAINT
!
! Add constraints for second-nearest neighbours - should correspond to bond angles
!
   DO J2=1,NBOND
      inloop: DO J3=J2+1,NBOND
        IF (BONDS(J2,1).EQ.BONDS(J3,1)) THEN
           ATOM1=BONDS(J2,2)
           ATOM2=BONDS(J3,2)
        ELSEIF (BONDS(J2,1).EQ.BONDS(J3,2)) THEN
           ATOM1=BONDS(J2,2)
           ATOM2=BONDS(J3,1)
        ELSEIF (BONDS(J2,2).EQ.BONDS(J3,1)) THEN
           ATOM1=BONDS(J2,1)
           ATOM2=BONDS(J3,2)
        ELSEIF (BONDS(J2,2).EQ.BONDS(J3,2)) THEN
           ATOM1=BONDS(J2,1)
           ATOM2=BONDS(J3,1)
        ELSE
           CYCLE inloop
        ENDIF
        IF (INTFROZEN(ATOM1).AND.INTFROZEN(ATOM2)) CYCLE ! no constraints between intfrozen atoms
        NCONSTRAINT=NCONSTRAINT+1
!       WRITE(*,'(A,2I6,A,I6)') 'intlbfgs> Adding constraint for second neighbours ',ATOM1,ATOM2, &
! &                     '  total=',NCONSTRAINT
         DS=SQRT((LXYZ(3*(ATOM1-1)+1)-LXYZ(3*(ATOM2-1)+1))**2 &
  &             +(LXYZ(3*(ATOM1-1)+2)-LXYZ(3*(ATOM2-1)+2))**2 &
  &             +(LXYZ(3*(ATOM1-1)+3)-LXYZ(3*(ATOM2-1)+3))**2)
         DF=SQRT((LXYZ(3*NATOMS+3*(ATOM1-1)+1)-LXYZ(3*NATOMS+3*(ATOM2-1)+1))**2 &
  &             +(LXYZ(3*NATOMS+3*(ATOM1-1)+2)-LXYZ(3*NATOMS+3*(ATOM2-1)+2))**2 &
  &             +(LXYZ(3*NATOMS+3*(ATOM1-1)+3)-LXYZ(3*NATOMS+3*(ATOM2-1)+3))**2)
         IF (NCONSTRAINT.GT.INTCONMAX) CALL CONDOUBLE
         CONI(NCONSTRAINT)=MIN(ATOM1,ATOM2)
         CONJ(NCONSTRAINT)=MAX(ATOM1,ATOM2)
         CONDISTREF(NCONSTRAINT)=(DF+DS)/2.0D0
         CONCUT(NCONSTRAINT)=ABS(DF-DS)/2.0D0
         IF (CONDISTREF(NCONSTRAINT).GT.MAXCONDIST) MAXCONDIST=CONDISTREF(NCONSTRAINT)
         IF (CONDISTREF(NCONSTRAINT).LT.MINCONDIST) MINCONDIST=CONDISTREF(NCONSTRAINT)
!        WRITE(*,'(A,2I6,A,2F12.2,A,F12.4,A,2I8)') ' intlbfgs> constrain distance for ',CONI(NCONSTRAINT), &
! &             CONJ(NCONSTRAINT),' values are ',DS,DF,' fraction=',2*ABS(DS-DF)/(DS+DF), &
! &            ' # second neighbour constraints, total=',QCISECOND,NCONSTRAINT
      ENDDO inloop
   ENDDO
   QCISECOND=NCONSTRAINT-QCIBONDS
   WRITE(*,'(A,2I6,A,I6)') 'intlbfgs> First and second neighbour constraints: ',QCIBONDS,QCISECOND,' total: ',NCONSTRAINT
   NDUMMY=NCONSTRAINT
   IF (CONFILET) THEN
      LUNIT=GETUNIT()
      OPEN(LUNIT,FILE='constraintfile',STATUS='OLD')
!
!  Additional amber constraints, e.g. cis/trans
!
      DO
         READ(LUNIT,*,END=534)  J2, J3
!
! Forbid constraints corresponding to atoms distant in sequence. Set INTCONSEP to number of sites to
! turn this off
!
         IF (J3-J2.GT.INTCONSEP) CYCLE
         IF (INTFROZEN(J2).AND.INTFROZEN(J3)) CYCLE ! no constraints between intfrozen atoms
         NCONSTRAINT=NCONSTRAINT+1
         DS=SQRT((LXYZ(3*(J2-1)+1)-LXYZ(3*(J3-1)+1))**2 &
  &             +(LXYZ(3*(J2-1)+2)-LXYZ(3*(J3-1)+2))**2 &
  &             +(LXYZ(3*(J2-1)+3)-LXYZ(3*(J3-1)+3))**2)
         DF=SQRT((LXYZ(3*NATOMS+3*(J2-1)+1)-LXYZ(3*NATOMS+3*(J3-1)+1))**2 &
  &             +(LXYZ(3*NATOMS+3*(J2-1)+2)-LXYZ(3*NATOMS+3*(J3-1)+2))**2 &
  &             +(LXYZ(3*NATOMS+3*(J2-1)+3)-LXYZ(3*NATOMS+3*(J3-1)+3))**2)
         IF (NCONSTRAINT.GT.INTCONMAX) CALL CONDOUBLE
         CONI(NCONSTRAINT)=J2
         CONJ(NCONSTRAINT)=J3
         CONDISTREF(NCONSTRAINT)=(DF+DS)/2.0D0
         CONCUT(NCONSTRAINT)=ABS(DF-DS)/2.0D0
         IF (CONDISTREF(NCONSTRAINT).GT.MAXCONDIST) MAXCONDIST=CONDISTREF(NCONSTRAINT)
         IF (CONDISTREF(NCONSTRAINT).LT.MINCONDIST) MINCONDIST=CONDISTREF(NCONSTRAINT)
         IF (DEBUG) WRITE(*,'(A,2I6,A,2F12.2,A,F12.4,A,I8)') ' intlbfgs> extra constraint distance for ',CONI(NCONSTRAINT), &
  &                     CONJ(NCONSTRAINT),' values are ',DS,DF,' fraction=',2*ABS(DS-DF)/(DS+DF), &
  &                  ' # constraints=',NCONSTRAINT
      ENDDO
534   CONTINUE
      CLOSE(LUNIT)
      IF (NCONSTRAINT-NDUMMY.GT.0) WRITE(*,'(A,I6,2(A,F15.5))') ' intlbfgs> Extra distance constraints: ',NCONSTRAINT-NDUMMY
      WRITE(*,'(A,I6,2(A,F15.5))') ' intlbfgs> Total distance constraints=',NCONSTRAINT,' shortest=',MINCONDIST,' longest=',MAXCONDIST
      CLOSE(LUNIT)
   ENDIF
ELSE IF (CONFILET) THEN 
    LUNIT=GETUNIT()
    OPEN(LUNIT,FILE='constraintfile',STATUS='OLD')
!
!  Add constraint for this distance to the list.
!
    DO 
       READ(LUNIT,*,END=531)  J2, J3
!
! Forbid constraints corresponding to atoms distant in sequence. Set INTCONSEP to number of sites to 
! turn this off
!
       IF (J3-J2.GT.INTCONSEP) CYCLE 
       IF (INTFROZEN(J2).AND.INTFROZEN(J3)) CYCLE ! no constraints between intfrozen atoms
       NCONSTRAINT=NCONSTRAINT+1
!      WRITE(*,'(A,2I6,A,I6)') 'intlbfgs> Adding constraint for atoms ',J2,J3,'  total=',NCONSTRAINT
       DS=SQRT((LXYZ(3*(J2-1)+1)-LXYZ(3*(J3-1)+1))**2 &
  &           +(LXYZ(3*(J2-1)+2)-LXYZ(3*(J3-1)+2))**2 &
  &           +(LXYZ(3*(J2-1)+3)-LXYZ(3*(J3-1)+3))**2) 
       DF=SQRT((LXYZ(3*NATOMS+3*(J2-1)+1)-LXYZ(3*NATOMS+3*(J3-1)+1))**2 &
  &           +(LXYZ(3*NATOMS+3*(J2-1)+2)-LXYZ(3*NATOMS+3*(J3-1)+2))**2 &
  &           +(LXYZ(3*NATOMS+3*(J2-1)+3)-LXYZ(3*NATOMS+3*(J3-1)+3))**2) 
       IF (NCONSTRAINT.GT.INTCONMAX) CALL CONDOUBLE
       CONI(NCONSTRAINT)=J2
       CONJ(NCONSTRAINT)=J3
       CONDISTREF(NCONSTRAINT)=(DF+DS)/2.0D0
       CONCUT(NCONSTRAINT)=ABS(DF-DS)/2.0D0
       IF (CONDISTREF(NCONSTRAINT).GT.MAXCONDIST) MAXCONDIST=CONDISTREF(NCONSTRAINT)
       IF (CONDISTREF(NCONSTRAINT).LT.MINCONDIST) MINCONDIST=CONDISTREF(NCONSTRAINT)
       WRITE(*,'(A,2I6,A,2F12.2,A,F12.4,A,I8)') ' intlbfgs> constrain distance for ',CONI(NCONSTRAINT), &
  &                 CONJ(NCONSTRAINT),' values are ',DS,DF,' fraction=',2*ABS(DS-DF)/(DS+DF), &
  &                ' # constraints=',NCONSTRAINT
    ENDDO
531 CONTINUE
    WRITE(*,'(A,I6,2(A,F15.5))') ' intlbfgs> Total distance constraints=',NCONSTRAINT, &
  &                               ' shortest=',MINCONDIST,' longest=',MAXCONDIST
    CLOSE(LUNIT)

ELSE IF (NCONGEOM.LT.2) THEN 
   DO J2=1,NATOMS
      DO J3=J2+1,NATOMS

         IF (J3-J2.GT.INTCONSEP) CYCLE ! forbid constraints corresponding to atoms distant in sequence
         IF (INTFROZEN(J2).AND.INTFROZEN(J3)) CYCLE ! no constraints between intfrozen atoms
         DS=SQRT((LXYZ(3*(J2-1)+1)-LXYZ(3*(J3-1)+1))**2 &
  &             +(LXYZ(3*(J2-1)+2)-LXYZ(3*(J3-1)+2))**2 &
  &             +(LXYZ(3*(J2-1)+3)-LXYZ(3*(J3-1)+3))**2) 
         IF (DS.GT.INTCONCUT) CYCLE ! don't allow constraints if either endpoint separation is too large DJW
         DF=SQRT((LXYZ(3*NATOMS+3*(J2-1)+1)-LXYZ(3*NATOMS+3*(J3-1)+1))**2 &
  &             +(LXYZ(3*NATOMS+3*(J2-1)+2)-LXYZ(3*NATOMS+3*(J3-1)+2))**2 &
  &             +(LXYZ(3*NATOMS+3*(J2-1)+3)-LXYZ(3*NATOMS+3*(J3-1)+3))**2) 
         IF (DF.GT.INTCONCUT) CYCLE ! don't allow constraints if either endpoint separation is too large DJW
!        IF (2.0D0*ABS(DS-DF)/(DS+DF).LT.LINTCONSTRAINTTOL) THEN
         WRITE(*,'(A,2I6,2G20.10)') 'intlbfgs> J2,J3,DS,DF=', J2,J3,DS,DF
         IF (ABS(DS-DF).LT.LINTCONSTRAINTTOL) THEN
!
!  Add constraint for this distance to the list.
!
            NCONSTRAINT=NCONSTRAINT+1
!           WRITE(*,'(A,2I6,A,I6)') 'intlbfgs> Adding constraint for atoms ',J2,J3,'  total=',NCONSTRAINT
            IF (NCONSTRAINT.GT.INTCONMAX) CALL CONDOUBLE
            CONI(NCONSTRAINT)=J2
            CONJ(NCONSTRAINT)=J3
            CONDISTREF(NCONSTRAINT)=(DF+DS)/2.0D0
            CONCUT(NCONSTRAINT)=ABS(DF-DS)/2.0D0
            IF (CONDISTREF(NCONSTRAINT).GT.MAXCONDIST) MAXCONDIST=CONDISTREF(NCONSTRAINT)
            IF (CONDISTREF(NCONSTRAINT).LT.MINCONDIST) MINCONDIST=CONDISTREF(NCONSTRAINT)
!           IF (DEBUG) PRINT '(A,2I6,A,2F12.2,A,F12.4,A,I8)',' intlbfgs> constrain distance for ',CONI(NCONSTRAINT), &
! &                 CONJ(NCONSTRAINT),' values are ',DS,DF,' fraction=',2*ABS(DS-DF)/(DS+DF), &
! &                ' # constraints=',NCONSTRAINT
         ENDIF
      ENDDO
   ENDDO
   IF (DEBUG) WRITE(*,'(A,I6,2(A,F15.5))') ' intlbfgs> Total distance constraints=',NCONSTRAINT, &
  &                                     ' shortest=',MINCONDIST,' longest=',MAXCONDIST
ELSE
   DO J2=1,NATOMS
      DO J3=J2+1,NATOMS
         IF (J3-J2.GT.INTCONSEP) CYCLE ! forbid constraints corresponding to atoms distant in sequence
         DSMIN=1.0D100
         DSMAX=-1.0D100
         DSMEAN=0.0D0
         DO J4=1,NCONGEOM
            DS=SQRT((CONGEOM(J4,3*(J2-1)+1)-CONGEOM(J4,3*(J3-1)+1))**2 &
  &                +(CONGEOM(J4,3*(J2-1)+2)-CONGEOM(J4,3*(J3-1)+2))**2 &
  &                +(CONGEOM(J4,3*(J2-1)+3)-CONGEOM(J4,3*(J3-1)+3))**2) 
            IF (DS.GT.DSMAX) DSMAX=DS
            IF (DS.LT.DSMIN) DSMIN=DS
            IF ((J4.GT.1).AND.(ABS(DSMIN-DSMAX).GT.LINTCONSTRAINTTOL)) GOTO 753 ! unconstrained
            IF (DS.GT.INTCONCUT) GOTO 753 ! don't allow constraints if any image separation is too large DJW
            DSMEAN=DSMEAN+DS
         ENDDO
!
!  Add constraint for this distance to the list if we make it to here.
!
         NCONSTRAINT=NCONSTRAINT+1
         WRITE(*,'(A,2I6,A,I6)') 'checkperc> Adding constraint for atoms ',J2,J3,'  total=',NCONSTRAINT
         IF (NCONSTRAINT.GT.INTCONMAX) CALL CONDOUBLE
         CONI(NCONSTRAINT)=J2
         CONJ(NCONSTRAINT)=J3
         CONDISTREF(NCONSTRAINT)=(DSMAX+DSMIN)/2.0D0 
         CONCUT(NCONSTRAINT)=(DSMAX-DSMIN)/2.0D0
         IF (CONDISTREF(NCONSTRAINT).GT.MAXCONDIST) MAXCONDIST=CONDISTREF(NCONSTRAINT)
         IF (CONDISTREF(NCONSTRAINT).LT.MINCONDIST) MINCONDIST=CONDISTREF(NCONSTRAINT)
         IF (DEBUG) WRITE(*,'(A,2I5,A,2F10.4,A,F12.4,A,I8)') &
  &                       ' checkperc> constrain atoms ',CONI(NCONSTRAINT), &
  &                       CONJ(NCONSTRAINT),' max, min ',DSMAX,DSMIN, &
  &                       ' cutoff=',CONCUT(NCONSTRAINT),' constraints=',NCONSTRAINT
753      CONTINUE
      ENDDO
   ENDDO
   CONIFIX(1:NCONSTRAINT)=CONI(1:NCONSTRAINT)
   CONJFIX(1:NCONSTRAINT)=CONJ(1:NCONSTRAINT)
   CONDISTREFFIX(1:NCONSTRAINT)=CONDISTREF(1:NCONSTRAINT)
   CONCUTFIX(1:NCONSTRAINT)=CONCUT(1:NCONSTRAINT)
   NCONSTRAINTFIX=NCONSTRAINT
ENDIF

IF (QCIADDREP.GT.0) THEN
   DMIN=1.0D100
   DO J2=1,QCIBONDS
!
! end point 1
!
      NI1=3*(CONI(J2)-1)
      NJ1=3*(CONJ(J2)-1)
      DO J3=J2+1,QCIBONDS
         IF (CONI(J3).EQ.CONI(J2)) CYCLE ! no extra terms for bonds with a common atom
         IF (CONI(J3).EQ.CONJ(J2)) CYCLE ! no extra terms for bonds with a common atom
         IF (CONJ(J3).EQ.CONI(J2)) CYCLE ! no extra terms for bonds with a common atom
         IF (CONJ(J3).EQ.CONJ(J2)) CYCLE ! no extra terms for bonds with a common atom
!
! end point 1
!
         NI2=3*(CONI(J3)-1)
         NJ2=3*(CONJ(J3)-1)
         DO J4=1,QCIADDREP
            X1=(J4*LXYZ(NI1+1)+(QCIADDREP+1-J4)*LXYZ(NJ1+1))/(QCIADDREP+1.0D0)
            Y1=(J4*LXYZ(NI1+2)+(QCIADDREP+1-J4)*LXYZ(NJ1+2))/(QCIADDREP+1.0D0)
            Z1=(J4*LXYZ(NI1+3)+(QCIADDREP+1-J4)*LXYZ(NJ1+3))/(QCIADDREP+1.0D0)
            DO J5=1,QCIADDREP
               X2=(J5*LXYZ(NI2+1)+(QCIADDREP+1-J5)*LXYZ(NJ2+1))/(QCIADDREP+1.0D0)
               Y2=(J5*LXYZ(NI2+2)+(QCIADDREP+1-J5)*LXYZ(NJ2+2))/(QCIADDREP+1.0D0)
               Z2=(J5*LXYZ(NI2+3)+(QCIADDREP+1-J5)*LXYZ(NJ2+3))/(QCIADDREP+1.0D0)
               D2=SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
               IF (D2.LT.DMIN) DMIN=D2
!              WRITE(*,'(A,2I6,A,4I6,A,2I6,A,F20.10)') 'intlbfgs> start constraints ',J2,J3,' atoms ', &
! &                                CONI(J2),CONJ(J2),CONI(J3),CONJ(J3),' J4,J5 ',J4,J5,' distance=',D2
           ENDDO
         ENDDO
      ENDDO
!
! end point 2
!
      NI1=3*(CONI(J2)-1)+3*NATOMS
      NJ1=3*(CONJ(J2)-1)+3*NATOMS
      DO J3=J2+1,QCIBONDS
         IF (CONI(J3).EQ.CONI(J2)) CYCLE ! no extra terms for bonds with a common atom
         IF (CONI(J3).EQ.CONJ(J2)) CYCLE ! no extra terms for bonds with a common atom
         IF (CONJ(J3).EQ.CONI(J2)) CYCLE ! no extra terms for bonds with a common atom
         IF (CONJ(J3).EQ.CONJ(J2)) CYCLE ! no extra terms for bonds with a common atom
!
! end point 2
!
     NI2=3*(CONI(J3)-1)
         NI2=3*(CONI(J3)-1)+3*NATOMS
         NJ2=3*(CONJ(J3)-1)+3*NATOMS
         DO J4=1,QCIADDREP
            X1=(J4*LXYZ(NI1+1)+(QCIADDREP+1-J4)*LXYZ(NJ1+1))/(QCIADDREP+1.0D0)
            Y1=(J4*LXYZ(NI1+2)+(QCIADDREP+1-J4)*LXYZ(NJ1+2))/(QCIADDREP+1.0D0)
            Z1=(J4*LXYZ(NI1+3)+(QCIADDREP+1-J4)*LXYZ(NJ1+3))/(QCIADDREP+1.0D0)
            DO J5=1,QCIADDREP
               X2=(J5*LXYZ(NI2+1)+(QCIADDREP+1-J5)*LXYZ(NJ2+1))/(QCIADDREP+1.0D0)
               Y2=(J5*LXYZ(NI2+2)+(QCIADDREP+1-J5)*LXYZ(NJ2+2))/(QCIADDREP+1.0D0)
               Z2=(J5*LXYZ(NI2+3)+(QCIADDREP+1-J5)*LXYZ(NJ2+3))/(QCIADDREP+1.0D0)
               D2=SQRT((X1-X2)**2+(Y1-Y2)**2+(Z1-Z2)**2)
               IF (D2.LT.DMIN) DMIN=D2
!              WRITE(*,'(A,2I6,A,4I6,A,2I6,A,F20.10)') 'intlbfgs> finish constraints ',J2,J3,' atoms ', &
! &                                CONI(J2),CONJ(J2),CONI(J3),CONJ(J3),' J4,J5 ',J4,J5,' distance=',D2
           ENDDO
         ENDDO
      ENDDO
   ENDDO
   WRITE(*,'(A,F20.10,A,F20.10)') 'intlbfgs> minimum decoration distance=',DMIN,' compared with cutoff ',QCIADDREPCUT
   QCIADDREPCUT=MIN(DMIN-1.0D-3,QCIADDREPCUT)
   WRITE(*,'(A,F20.10)') 'intlbfgs> cutoff after setup is ',QCIADDREPCUT
ENDIF
!
! Check that we have a percolating constraint network. If not, increase the tolerance and try again!
! Calculate minimum number of steps of each atom from number 1 or any frozen atom.
!
NDIST1(1:NATOMS)=1000000
IF (NQCIFREEZE.EQ.0) THEN
   NDIST1(1)=0
ELSE
   DO J1=1,NATOMS
      IF (INTFROZEN(J1)) NDIST1(J1)=0
   ENDDO
ENDIF
NCYCLE=0
5    CHANGED=.FALSE.
NCYCLE=NCYCLE+1
DMIN1=100000
DMAX1=0
NUNCON1=0
DO J1=1,NATOMS
   IF (NDIST1(J1).EQ.0) CYCLE ! minimum 1
   DO J2=1,NCONSTRAINT
      IF (CONI(J2).EQ.J1) THEN
         IF (NDIST1(CONJ(J2))+1.LT.NDIST1(J1)) THEN
            CHANGED=.TRUE.
            NDIST1(J1)=NDIST1(CONJ(J2))+1
         ENDIF
      ELSE IF (CONJ(J2).EQ.J1) THEN
         IF (NDIST1(CONI(J2))+1.LT.NDIST1(J1)) THEN
            CHANGED=.TRUE.
            NDIST1(J1)=NDIST1(CONI(J2))+1
         ENDIF
      ENDIF
   ENDDO
   IF ((NDIST1(J1).GT.DMAX1).AND.(NDIST1(J1).NE.1000000)) DMAX1=NDIST1(J1)
   IF (NDIST1(J1).LT.DMIN1) DMIN1=NDIST1(J1)
   IF (NDIST1(J1).EQ.1000000) NUNCON1=NUNCON1+1
ENDDO
IF (CHANGED) GOTO 5
  IF (DEBUG) WRITE(*,'(3(A,I8))') ' checkperc> steps to atom 1 converged in ',NCYCLE-1, &
    &               ' cycles; maximum=',DMAX1,' disconnected=',NUNCON1
IF (NUNCON1.GT.0) THEN
   LINTCONSTRAINTTOL=LINTCONSTRAINTTOL*1.1D0
   IF (DEBUG) WRITE(*,'(A,F15.5)') ' checkperc> increasing the local constraint tolerance parameter to ',LINTCONSTRAINTTOL
   IF (LINTCONSTRAINTTOL.GT.100.0D0) THEN
      WRITE(*,'(A,G20.10)') 'checkperc> likely ERROR *** LINTCONSTRAINTTOL=',LINTCONSTRAINTTOL
      STOP
   ENDIF
   GOTO 51
ENDIF
! IF (DEBUG) WRITE(*,'(A,F15.5)') ' checkperc> Final constraint tolerance parameter ',LINTCONSTRAINTTOL

! WRITE(*,'(A,I6,3(A,F15.5))') ' checkperc> Total distance constraints=',NCONSTRAINT, &
!   &                    ' shortest=',MINCONDIST,' longest=',MAXCONDIST,' tolerance=',LINTCONSTRAINTTOL

CALLED=.TRUE.

END SUBROUTINE CHECKPERC

SUBROUTINE MAKESTEP(NITERDONE,POINT,DIAG,INTIMAGE,SEARCHSTEP,G,GTMP,STP,GDIF,NPT,D,RHO1,ALPHA)
USE KEY, ONLY : INTMUPDATE, INTDGUESS
USE COMMONS, ONLY: NATOMS
IMPLICIT NONE
INTEGER NITERDONE, POINT, BOUND, NPT, D, CP, INTIMAGE, I
DOUBLE PRECISION DIAG(3*NATOMS*INTIMAGE),SEARCHSTEP(0:INTMUPDATE,(3*NATOMS)*INTIMAGE),G((3*NATOMS)*INTIMAGE), &
  &  GTMP(3*NATOMS*INTIMAGE), GNORM, STP(3*NATOMS*INTIMAGE), YS, GDIF(0:INTMUPDATE,(3*NATOMS)*INTIMAGE), YY, &
  &  SQ, YR, BETA
DOUBLE PRECISION, DIMENSION(INTMUPDATE)     :: RHO1,ALPHA
LOGICAL CHANGEIMAGE
SAVE

MAIN: IF (NITERDONE==1) THEN
     POINT = 0
     DIAG(1:D)=INTDGUESS
     SEARCHSTEP(0,1:D)= -G(1:D)*INTDGUESS            ! NR STEP FOR DIAGONAL INVERSE HESSIAN
     GTMP(1:D)        = SEARCHSTEP(0,1:D)
     GNORM            = MAX(SQRT(DOT_PRODUCT(G(1:D),G(1:D))),1.0D-100)
     STP(1:D)         = MIN(1.0D0/GNORM, GNORM) ! MAKE THE FIRST GUESS FOR THE STEP LENGTH CAUTIOUS
ELSE MAIN
     BOUND=NITERDONE-1
     IF (NITERDONE.GT.INTMUPDATE) BOUND=INTMUPDATE
     YS=DOT_PRODUCT( GDIF(NPT/D,:), SEARCHSTEP(NPT/D,:)  )
     IF (YS==0.0D0) YS=1.0D0
    
! Update estimate of diagonal inverse Hessian elements.
! We divide by both YS and YY at different points, so they had better not be zero!

     YY=DOT_PRODUCT( GDIF(NPT/D,:) , GDIF(NPT/D,:) )
     IF (YY==0.0D0) YY=1.0D0
!    DIAG = ABS(YS/YY)
     DIAG(1) = YS/YY
      
! COMPUTE -H*G USING THE FORMULA GIVEN IN: Nocedal, J. 1980, 
! "Updating quasi-Newton matrices with limited storage",
! Mathematics of Computation, Vol.35, No.151, pp. 773-782

     CP= POINT; IF (POINT==0) CP = INTMUPDATE
     RHO1(CP)=1.0D0/YS
     GTMP(1:D) = -G(1:D)
     CP= POINT 
                   
     DO I= 1,BOUND 
          CP = CP - 1; IF (CP == -1) CP = INTMUPDATE - 1
          SQ= DOT_PRODUCT( SEARCHSTEP(CP,1:D),GTMP(1:D) )
          ALPHA(CP+1) = RHO1(CP+1) * SQ
          GTMP(1:D)        = -ALPHA(CP+1)*GDIF(CP,1:D) + GTMP(1:D)
     ENDDO
              
     GTMP(1:D)=DIAG(1)*GTMP(1:D)

     DO I=1,BOUND
          YR= DOT_PRODUCT( GDIF(CP,1:D) , GTMP )
          BETA= RHO1(CP+1)*YR
          BETA= ALPHA(CP+1)-BETA
!         WRITE(*,'(A,I8,4G20.10)') 'makestep> I,YR,BETA,RHO1,ALPHA=',I,YR,BETA,RHO1(CP+1),ALPHA(CP+1)
          GTMP(1:D) = BETA*SEARCHSTEP(CP,1:D) + GTMP(1:D)
          CP=CP+1
!         IF (CP==M) CP=0
          IF (CP==INTMUPDATE) CP=0
     ENDDO
              
     STP(1:D) = 1.0D0
ENDIF MAIN

!  Store the new search direction
IF (NITERDONE.GT.1) SEARCHSTEP(POINT,1:D)=GTMP(1:D)

END SUBROUTINE MAKESTEP

!
! Set up start and finish atom indices for each group once in intlbfgs
! STARTGROUP(DOGROUP) and ENDGROUP(DOGROUP)
!

SUBROUTINE QCIOPTPERM(COORDSB,COORDSA,NATOMS,DEBUG,DOGROUP,PERM,STARTGROUP,IM)
USE KEY,ONLY : NPERMGROUP, NPERMSIZE, PERMGROUP, NSETS, SETS
IMPLICIT NONE
INTEGER NATOMS
INTEGER PERM(NATOMS), DOGROUP, J1, J2, STARTGROUP(NPERMGROUP), J3, I1, I2, I3, IM
DOUBLE PRECISION COORDSA(3*NATOMS), COORDSB(3*NATOMS), D123, D132, D231, D213, D312, D321, DMIN, D12, D21
LOGICAL DEBUG, IDENTITY

IDENTITY=.TRUE.
DO J1=1,NATOMS
   PERM(J1)=J1
ENDDO

IF (NPERMSIZE(DOGROUP).EQ.3) THEN
   I1=PERMGROUP(STARTGROUP(DOGROUP))
   I2=PERMGROUP(STARTGROUP(DOGROUP)+1)
   I3=PERMGROUP(STARTGROUP(DOGROUP)+2)
! permutation 1 2 3
   D123=(COORDSB(3*(I1-1)+1)-COORDSA(3*(I1-1)+1))**2+(COORDSB(3*(I1-1)+2)-COORDSA(3*(I1-1)+2))**2+(COORDSB(3*(I1-1)+3)-COORDSA(3*(I1-1)+3))**2 &
 &     +(COORDSB(3*(I2-1)+1)-COORDSA(3*(I2-1)+1))**2+(COORDSB(3*(I2-1)+2)-COORDSA(3*(I2-1)+2))**2+(COORDSB(3*(I2-1)+3)-COORDSA(3*(I2-1)+3))**2 &
 &     +(COORDSB(3*(I3-1)+1)-COORDSA(3*(I3-1)+1))**2+(COORDSB(3*(I3-1)+2)-COORDSA(3*(I3-1)+2))**2+(COORDSB(3*(I3-1)+3)-COORDSA(3*(I3-1)+3))**2
   DMIN=D123
! permutation 1 3 2
   D132=(COORDSB(3*(I1-1)+1)-COORDSA(3*(I1-1)+1))**2+(COORDSB(3*(I1-1)+2)-COORDSA(3*(I1-1)+2))**2+(COORDSB(3*(I1-1)+3)-COORDSA(3*(I1-1)+3))**2 &
 &     +(COORDSB(3*(I2-1)+1)-COORDSA(3*(I3-1)+1))**2+(COORDSB(3*(I2-1)+2)-COORDSA(3*(I3-1)+2))**2+(COORDSB(3*(I2-1)+3)-COORDSA(3*(I3-1)+3))**2 &
 &     +(COORDSB(3*(I3-1)+1)-COORDSA(3*(I2-1)+1))**2+(COORDSB(3*(I3-1)+2)-COORDSA(3*(I2-1)+2))**2+(COORDSB(3*(I3-1)+3)-COORDSA(3*(I2-1)+3))**2
   IF (D132.LT.DMIN) THEN 
      DMIN=D132
      PERM(I1)=I1
      PERM(I2)=I3
      PERM(I3)=I2
      IDENTITY=.FALSE.
   ENDIF
! permutation 2 1 3
   D213=(COORDSB(3*(I1-1)+1)-COORDSA(3*(I2-1)+1))**2+(COORDSB(3*(I1-1)+2)-COORDSA(3*(I2-1)+2))**2+(COORDSB(3*(I1-1)+3)-COORDSA(3*(I2-1)+3))**2 &
 &     +(COORDSB(3*(I2-1)+1)-COORDSA(3*(I1-1)+1))**2+(COORDSB(3*(I2-1)+2)-COORDSA(3*(I1-1)+2))**2+(COORDSB(3*(I2-1)+3)-COORDSA(3*(I1-1)+3))**2 &
 &     +(COORDSB(3*(I3-1)+1)-COORDSA(3*(I3-1)+1))**2+(COORDSB(3*(I3-1)+2)-COORDSA(3*(I3-1)+2))**2+(COORDSB(3*(I3-1)+3)-COORDSA(3*(I3-1)+3))**2
   IF (D213.LT.DMIN) THEN 
      DMIN=D213
      PERM(I1)=I2
      PERM(I2)=I1
      PERM(I3)=I3
      IDENTITY=.FALSE.
   ENDIF
! permutation 2 3 1
   D231=(COORDSB(3*(I1-1)+1)-COORDSA(3*(I2-1)+1))**2+(COORDSB(3*(I1-1)+2)-COORDSA(3*(I2-1)+2))**2+(COORDSB(3*(I1-1)+3)-COORDSA(3*(I2-1)+3))**2 &
 &     +(COORDSB(3*(I2-1)+1)-COORDSA(3*(I3-1)+1))**2+(COORDSB(3*(I2-1)+2)-COORDSA(3*(I3-1)+2))**2+(COORDSB(3*(I2-1)+3)-COORDSA(3*(I3-1)+3))**2 &
 &     +(COORDSB(3*(I3-1)+1)-COORDSA(3*(I1-1)+1))**2+(COORDSB(3*(I3-1)+2)-COORDSA(3*(I1-1)+2))**2+(COORDSB(3*(I3-1)+3)-COORDSA(3*(I1-1)+3))**2
   IF (D231.LT.DMIN) THEN 
      DMIN=D231
      PERM(I1)=I2
      PERM(I2)=I3
      PERM(I3)=I1
      IDENTITY=.FALSE.
   ENDIF
! permutation 3 2 1
   D321=(COORDSB(3*(I1-1)+1)-COORDSA(3*(I3-1)+1))**2+(COORDSB(3*(I1-1)+2)-COORDSA(3*(I3-1)+2))**2+(COORDSB(3*(I1-1)+3)-COORDSA(3*(I3-1)+3))**2 &
 &     +(COORDSB(3*(I2-1)+1)-COORDSA(3*(I2-1)+1))**2+(COORDSB(3*(I2-1)+2)-COORDSA(3*(I2-1)+2))**2+(COORDSB(3*(I2-1)+3)-COORDSA(3*(I2-1)+3))**2 &
 &     +(COORDSB(3*(I3-1)+1)-COORDSA(3*(I1-1)+1))**2+(COORDSB(3*(I3-1)+2)-COORDSA(3*(I1-1)+2))**2+(COORDSB(3*(I3-1)+3)-COORDSA(3*(I1-1)+3))**2
   IF (D321.LT.DMIN) THEN 
      DMIN=D321
      PERM(I1)=I3
      PERM(I2)=I2
      PERM(I3)=I1
      IDENTITY=.FALSE.
   ENDIF
! permutation 3 1 2
   D312=(COORDSB(3*(I1-1)+1)-COORDSA(3*(I3-1)+1))**2+(COORDSB(3*(I1-1)+2)-COORDSA(3*(I3-1)+2))**2+(COORDSB(3*(I1-1)+3)-COORDSA(3*(I3-1)+3))**2 &
 &     +(COORDSB(3*(I2-1)+1)-COORDSA(3*(I1-1)+1))**2+(COORDSB(3*(I2-1)+2)-COORDSA(3*(I1-1)+2))**2+(COORDSB(3*(I2-1)+3)-COORDSA(3*(I1-1)+3))**2 &
 &     +(COORDSB(3*(I3-1)+1)-COORDSA(3*(I2-1)+1))**2+(COORDSB(3*(I3-1)+2)-COORDSA(3*(I2-1)+2))**2+(COORDSB(3*(I3-1)+3)-COORDSA(3*(I2-1)+3))**2
   IF (D312.LT.DMIN) THEN 
      DMIN=D312
      PERM(I1)=I3
      PERM(I2)=I1
      PERM(I3)=I2
      IDENTITY=.FALSE.
   ENDIF
!  IF (.NOT.IDENTITY) WRITE(*,'(A,2I6,6G15.5,3I6)') ' qcioptperm> images,D123,D132,D213,D231,D312,D321,permutation: ', &
!&               IM,IM+1,D123,D132,D213,D231,D312,D321,PERM(I1),PERM(I2),PERM(I3)
   WRITE(*,'(A,2I6,6G15.5,3I6)') ' qcioptperm> images,D123,D132,D213,D231,D312,D321,permutation: ', &
 &               IM,IM+1,D123,D132,D213,D231,D312,D321,PERM(I1),PERM(I2),PERM(I3)
ELSE IF (NPERMSIZE(DOGROUP).EQ.2) THEN
! permutation 1 2
   I1=PERMGROUP(STARTGROUP(DOGROUP))
   I2=PERMGROUP(STARTGROUP(DOGROUP)+1)
   D12=(COORDSB(3*(I1-1)+1)-COORDSA(3*(I1-1)+1))**2+(COORDSB(3*(I1-1)+2)-COORDSA(3*(I1-1)+2))**2+(COORDSB(3*(I1-1)+3)-COORDSA(3*(I1-1)+3))**2 &
 &    +(COORDSB(3*(I2-1)+1)-COORDSA(3*(I2-1)+1))**2+(COORDSB(3*(I2-1)+2)-COORDSA(3*(I2-1)+2))**2+(COORDSB(3*(I2-1)+3)-COORDSA(3*(I2-1)+3))**2 
   DMIN=D12
! permutation 2 1
   D21=(COORDSB(3*(I1-1)+1)-COORDSA(3*(I2-1)+1))**2+(COORDSB(3*(I1-1)+2)-COORDSA(3*(I2-1)+2))**2+(COORDSB(3*(I1-1)+3)-COORDSA(3*(I2-1)+3))**2 &
 &    +(COORDSB(3*(I2-1)+1)-COORDSA(3*(I1-1)+1))**2+(COORDSB(3*(I2-1)+2)-COORDSA(3*(I1-1)+2))**2+(COORDSB(3*(I2-1)+3)-COORDSA(3*(I1-1)+3))**2 
   IF (D21.LT.DMIN) THEN 
      DMIN=D21
      PERM(I1)=I2
      PERM(I2)=I1
      IF (NSETS(DOGROUP).GT.0) THEN
         DO J3=1,NSETS(DOGROUP) ! can be 2 or 3
!           WRITE(*,'(A,4I6)') 'I1,I2,sets(I1,j3),sets(I2,j3)=',I1,I2,sets(I1,j3),sets(I2,j3)
            PERM(SETS(I1,J3))=SETS(I2,J3)
            PERM(SETS(I2,J3))=SETS(I1,J3)
         ENDDO
      ENDIF

   ENDIF
   IF (.NOT.IDENTITY) THEN
      WRITE(*,'(A,2I6,2G15.5,2I6)') ' qcioptperm> images,D12,D21,permutation: ',IM,IM+1,D12,D21,PERM(I1),PERM(I2)
      IF (NSETS(DOGROUP).GT.0) WRITE(*,'(A,6I6)') ' qcioptperm> associated permutation: ', &
  &                         (PERM(SETS(I1,J3)),PERM(SETS(I2,J3)),J3=1,NSETS(DOGROUP))  
   ENDIF
ELSE
   WRITE(*,'(A,I6,A,I6)') ' qcioptperm> Unknown group size ',NPERMSIZE(DOGROUP),' for group ',DOGROUP
   STOP
ENDIF

END SUBROUTINE QCIOPTPERM

SUBROUTINE TRILATERATION(P1,P2,P3,R1,R2,R3,SOL1,SOL2,FTEST)
IMPLICIT NONE
DOUBLE PRECISION P1(3), P2(3), P3(3), R1, R2, R3, SOL1(3), SOL2(3), I
DOUBLE PRECISION  TEMP1(3), EX(3), EY(3), EZ(3), DUMMY, TEMP2(3), TEMP3(3), D, J, X, Y, Z, TEMP4
LOGICAL FTEST

! # Find the intersection of three spheres                 
! # P1,P2,P3 are the centers, r1,r2,r3 are the radii       
! # Implementaton based on Wikipedia Trilateration article.                              

FTEST=.FALSE.
TEMP1(1:3)=P2(1:3)-P1(1:3)
D=SQRT( TEMP1(1)**2+TEMP1(2)**2+TEMP1(3)**2 )
EX(1:3)=TEMP1(1:3)/D
TEMP2(1:3)=P3(1:3)-P1(1:3)
I=EX(1)*TEMP2(1)+EX(2)*TEMP2(2)+EX(3)*TEMP2(3)
TEMP3(1:3)=TEMP2(1:3)-I*EX(1:3)
DUMMY=SQRT( TEMP3(1)**2+TEMP3(2)**2+TEMP3(3)**2 )
EY(1:3)=TEMP3(1:3)/DUMMY
EZ(1)= EX(2)*EY(3)-EX(3)*EY(2)
EZ(2)=-EX(1)*EY(3)+EX(3)*EY(1)
EZ(3)= EX(1)*EY(2)-EX(2)*EY(1)
J=EY(1)*TEMP2(1)+EY(2)*TEMP2(2)+EY(3)*TEMP2(3)
X=(R1*R1 - R2*R2 + D*D) / (2.0D0*D)
Y=(R1*R1 - R3*R3 -2.0D0*I*X + I*I + J*J) / (2.0D0*J)
TEMP4=R1*R1 - X*X - Y*Y

! WRITE (*,'(A,9G15.5)') 'trilateration> p1, p2, p3: ',P1(1:3),P2(1:3),P3(1:3)
! WRITE (*,'(A,9G15.5)') 'trilateration> ex, ey, ez: ',EX(1:3),EY(1:3),EZ(1:3)
! WRITE (*,'(A,9G15.5)') 'trilateration> norms:      ',EX(1)**2+EX(2)**2+EX(3)**2,EY(1)**2+EY(2)**2+EY(3)**2,EZ(1)**2+EZ(2)**2+EZ(3)**2
! WRITE (*,'(A,9G15.5)') 'trilateration> r1, r2, r3: ',R1,R2,R3
! WRITE (*,'(A,9G15.5)') 'trilateration> X, Y, TEMP4:    ',X,Y,TEMP4

! PRINT *,'TEMP4=',TEMP4
!PRINT *,'TEMP4.LT.0.0D0=',TEMP4.LT.0.0D0

IF (TEMP4.LT.0.0D0) THEN
   FTEST=.TRUE.
   RETURN
ELSE
   FTEST=.FALSE.
   Z=SQRT(TEMP4)
   SOL1(1:3)=P1(1:3) + X*EX(1:3) + Y*EY(1:3) + Z*EZ(1:3)
   SOL2(1:3)=P1(1:3) + X*EX(1:3) + Y*EY(1:3) - Z*EZ(1:3)
!  PRINT *,'Z=',Z
ENDIF

END SUBROUTINE TRILATERATION 



SUBROUTINE INTRWG2(NACTIVE,NITER,INTIMAGE,XYZ,TURNONORDER,NCONOFF)
USE PORFUNCS
USE KEY,ONLY: STOCKT,STOCKAAT, RBAAT, ATOMACTIVE, NCONSTRAINT, CONACTIVE, NREPULSIVE, NNREPULSIVE, REPI, REPJ, REPCUT, NREPCUT, &
  &           NREPMAX, NREPI, NREPJ, INTFROZEN, CONOFFLIST,CONOFFTRIED, KINT,INTCONSTRAINTREP
USE COMMONS, ONLY: NATOMS, DEBUG
IMPLICIT NONE
INTEGER NCONOFF
CHARACTER(LEN=10) :: XYZFILE   = 'int.xyz   '
CHARACTER(LEN=10) :: QCIFILE   = 'QCIdump   '
INTEGER,INTENT(IN) :: NITER, TURNONORDER(NATOMS)
INTEGER :: J1,J2,INTIMAGE,J3,NACTIVE,LUNIT,GETUNIT
CHARACTER(LEN=80) :: FILENAME,DUMMYS
DOUBLE PRECISION XYZ((3*NATOMS)*(INTIMAGE+2))

FILENAME=XYZFILE

! IF (NITER.GT.0) THEN
   WRITE(DUMMYS,'(I8)') NITER
   FILENAME='int.' // TRIM(ADJUSTL(DUMMYS)) // '.xyz' ! so that vmd recognises the file type!
! ENDIF
LUNIT=GETUNIT()
OPEN(UNIT=LUNIT,FILE=TRIM(ADJUSTL(FILENAME)),STATUS='replace')
DO J2=1,INTIMAGE+2
!  WRITE(LUNIT,'(i4/)') NACTIVE
   WRITE(LUNIT,'(i4/)') NATOMS
   DO J3=1,NATOMS
      IF (ATOMACTIVE(J3)) THEN
         WRITE(LUNIT,'(A5,1X,3F20.10)') 'LA   ',XYZ((J2-1)*3*NATOMS+3*(J3-1)+1),XYZ((J2-1)*3*NATOMS+3*(J3-1)+2), &  
  &                                                                   XYZ((J2-1)*3*NATOMS+3*(J3-1)+3)  
      ELSE
         WRITE(LUNIT,'(A5,1X,3F20.10)') 'DU   ',XYZ((J2-1)*3*NATOMS+3*(J3-1)+1),XYZ((J2-1)*3*NATOMS+3*(J3-1)+2), &  
  &                                                                   XYZ((J2-1)*3*NATOMS+3*(J3-1)+3)  
      ENDIF
   ENDDO
ENDDO

WRITE(*,*) 'rwg> Interpolated image coordinates were saved to xyz file "'//TRIM(FILENAME)//'"'

CLOSE(LUNIT)

FILENAME=QCIFILE
LUNIT=GETUNIT()
! IF (NITER.GT.0) THEN
   WRITE(DUMMYS,'(I8)') NITER
   FILENAME='QCIdump.' // TRIM(ADJUSTL(DUMMYS)) 
! ENDIF
OPEN(UNIT=LUNIT,FILE=TRIM(ADJUSTL(FILENAME)),STATUS='replace')

IF (DEBUG) WRITE(*,'(A,I10,A)') ' intlbfgs> dumping state for ',NACTIVE,' active atoms'
WRITE(LUNIT,'(I10)') NACTIVE
! IF (DEBUG)   WRITE(*,'(A,I10,A)') ' intlbfgs> dumping spring constant and repulsive prefactor'
WRITE(LUNIT,'(2G20.10)') KINT, INTCONSTRAINTREP
! WRITE(*,'(A,I10,A)') ' intlbfgs> dumping turnonorder for ',NACTIVE,' active atoms'
WRITE(LUNIT,'(12I8)') TURNONORDER(1:NACTIVE)
! WRITE(*,'(A)') ' intlbfgs> dumping atomactive'
WRITE(LUNIT,'(12L5)') ATOMACTIVE(1:NATOMS)
WRITE(LUNIT,'(I10)') NCONSTRAINT
! WRITE(*,'(A,I10,A)') ' intlbfgs> dumping conactive for ',NCONSTRAINT,' constraints'
WRITE(LUNIT,'(12L5)') CONACTIVE(1:NCONSTRAINT)

   WRITE(LUNIT,'(3I12,G20.10)') NREPULSIVE,NNREPULSIVE,NREPMAX
   ! WRITE(*,'(A,3I10,G20.10)') 'intlbfgs> dumping NREPULSIVE,NNREPULSIVE,NREPMAX=',NREPULSIVE,NNREPULSIVE,NREPMAX

   WRITE(LUNIT,'(12I8)') REPI(1:NREPULSIVE)
   ! WRITE(*,'(A)') ' intlbfgs> dumped REPI:'
   WRITE(LUNIT,'(12I8)') REPJ(1:NREPULSIVE)
   ! WRITE(*,'(A)') ' intlbfgs> dumped REPJ:'
   WRITE(LUNIT,'(12I8)') NREPI(1:NNREPULSIVE)
   ! WRITE(*,'(A)') ' intlbfgs> dumped NREPI:'
   WRITE(LUNIT,'(12I8)') NREPJ(1:NNREPULSIVE)
   ! WRITE(*,'(A)') ' intlbfgs> dumped NREPJ:'

   WRITE(LUNIT,'(6G20.10)') REPCUT(1:NREPULSIVE)
   ! WRITE(*,'(A)') ' intlbfgs> dumped REPCUT:'
   WRITE(LUNIT,'(6G20.10)') NREPCUT(1:NNREPULSIVE)
   ! WRITE(*,'(A)') ' intlbfgs> dumped NREPCUT:'

   WRITE(LUNIT,'(12L5)') INTFROZEN(1:NATOMS)
   ! WRITE(*,'(A)') ' intlbfgs> dumped INTFROZEN'

   WRITE(LUNIT,'(I8)') NCONOFF
   IF (NCONOFF.GT.0) WRITE(LUNIT,'(12I8)') CONOFFLIST(1:NCONOFF)
   IF (NCONOFF.GT.0) WRITE(LUNIT,'(12L5)') CONOFFTRIED(1:NCONOFF)
   ! WRITE(*,'(A)') ' intlbfgs> dumped NCONOFF and CONOFFLIST'

CLOSE(LUNIT)

END SUBROUTINE INTRWG2
