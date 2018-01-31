!   OPTIM: A program for optimizing geometries and calculating reaction pathways
!   Copyright (C) 1999- David J. Wales
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
! Reaction path Hamiltonian formulation as per Page and McIver, JCP, 88, 922, 1988
!
SUBROUTINE RPH
USE MODCHARMM
USE MODHESS
USE KEY
USE COMMONS
USE PORFUNCS
USE GENRIGID
IMPLICIT NONE
INTEGER INFO, NEXMODES, I1, ATFRAME, KAPPA

INTEGER J1, J2, NSTRUCTREF, HIGHESTREF, NTOT, NUPDATE, IMCSTEP, K, LUNIT, IDIFF, K1, K2, IMAX, NTS, NMIN, J3
INTEGER GETUNIT, IACCEPT, IBININDEX, NDUMMY, IWINDOW, IMAGEMIN, IMAGEMAX, IMAGEMINO, IREJDIST, NPATH
INTEGER LSTRUCTREF, IREJCHECK, NCHECKSTRUCT, NCHECKPLUS, NCHECKMINUS, DOTS, STARTTS, ENDTS, LNPATH
INTEGER, ALLOCATABLE :: NEAREST(:), REFSTRUCTFRAME(:), TEMPN(:)
DOUBLE PRECISION BINLABELQORDER(RPHQBINS), SAVEHESS(NOPT,NOPT)
DOUBLE PRECISION, ALLOCATABLE :: EOFS(:), PATHFRAMES(:,:), LCOORDSCHECK(:,:), LEOFS(:)
DOUBLE PRECISION, ALLOCATABLE :: DIST(:), DISTO(:)
DOUBLE PRECISION QORDERHIST(RPHQBINS)
DOUBLE PRECISION, ALLOCATABLE :: DISTFRAME(:,:), LDISTFRAME(:,:), VSAVE(:,:)
DOUBLE PRECISION, ALLOCATABLE :: COORDSREF(:,:), VREF(:), LCOORDSREF(:,:), LVREF(:)
DOUBLE PRECISION MCCOORDS(3*NATOMS), MCCOORDSO(3*NATOMS), X(NATOMS), Y(NATOMS), Z(NATOMS), &
  &                                                   XO(NATOMS), YO(NATOMS), ZO(NATOMS), &
  &              GRAD(3*NATOMS), POTEL, RMS, VNEW(3*NATOMS), &
  &              DUMMY, WAC, DIST2, RMAT(3,3), CDUMMY(3*NATOMS), DSAVE, &
  &              DISTMIN, DISTMAX, XDISTMIN, XDISTMAX, HISTINT, D1, D2, D3, &
  &              DUMMY2, DMIN, MEANBIAS, MEANBIASEQ, QORDERINT, QORDERO, SUM1, SUM2, &
  &              BIASBEST, BIASSTEP, SUMBEST, MAXDEV, DISTCHECK, QORDERMEAN, QORDERSQMEAN, &
  &              DISTTOTAL, SLICELABEL(RPHSLICES), DISTTOTAL2
DOUBLE PRECISION, ALLOCATABLE :: DSTRUCT(:)
LOGICAL RECOUNT, YESNO, BIAST, LASTMIN, LASTTS, PERMDISTSAVE, LPERMDISTSAVE
CHARACTER (LEN=8) SDUMMY
CHARACTER (LEN=80) FNAME
CHARACTER (LEN=5) SSYM
DOUBLE PRECISION QTEMP(3,NATOMS),QORDER,Q(3*NATOMS),ENERGY,EVALUES(3*NATOMS),TEMPA(9*NATOMS),FRACTION
DOUBLE PRECISION PFUNC(0:RPHSLICES+1), PSUM, PMEAN, QORDERSLICES(0:RPHSLICES+1)
LOGICAL PROJGRAD

DOUBLE PRECISION SLICEWIDTH, PROD

GRAD(1:3*NATOMS)=0.0D0
RMS=0.0D0
PERMDISTSAVE=PERMDIST 
LPERMDISTSAVE=LPERMDIST

PRINT '(A,G20.10)','RPH> Free energy profile for V(s) pathway at=',RPHTEMP
PRINT '(A,I10)',   'RPH> Number of slices=                       ',RPHSLICES
PRINT '(A,G20.10)','RPH> Bottom of order parameter range=        ',RPHQMIN
PRINT '(A,G20.10)','RPH> Top of order parameter range=           ',RPHQMAX

!
! SLICEWIDTH is the width used the slices along the path
!
! Parse path.xyz to discover the number of frames saved.
! Identify the number of maxima and minima and set NSTRUCTREF.
! Read the reference coordinates and energies into
! COORDSREF and VREF arrays after declaring them.
! Set integer highest to the highest energy reference structure.
!
INQUIRE(FILE='path.xyz',EXIST=YESNO)
IF (.NOT.YESNO) THEN
   PRINT '(A)','RPH> ERROR *** no path.xyz file'
   STOP
ENDIF
LUNIT=GETUNIT()
OPEN(LUNIT,FILE='path.xyz',STATUS='OLD')
NPATH=0
DO 
   READ(LUNIT,*,END=20) SDUMMY
   READ(LUNIT,*) SDUMMY
   READ(LUNIT,*) (SSYM,CDUMMY(1:3),J1=1,NATOMS)
   NPATH=NPATH+1
ENDDO
20 PRINT '(A,I10)','RPH> Number of frames in path.xyz=',NPATH
CLOSE(LUNIT)
ALLOCATE(EOFS(NPATH),PATHFRAMES(NPATH,3*NATOMS))
LUNIT=GETUNIT()
OPEN(LUNIT,FILE='path.xyz',STATUS='OLD')
DO J1=1,NPATH
   READ(LUNIT,*) SDUMMY
   READ(LUNIT,*) SDUMMY,EOFS(J1)
   DO J2=1,NATOMS
      READ(LUNIT,'(A5,1X,3F20.10)') SSYM,PATHFRAMES(J1,3*(J2-1)+1:3*(J2-1)+3)
   ENDDO
ENDDO
CLOSE(LUNIT)

NSTRUCTREF=2
DO J1=2,NPATH-1
   IF (EOFS(J1-1)+EDIFFTOL*10.0D0<EOFS(J1).AND.EOFS(J1)>EOFS(J1+1)+EDIFFTOL*10.0D0) NSTRUCTREF=NSTRUCTREF+1
   IF (EOFS(J1-1)+EDIFFTOL*10.0D0>EOFS(J1).AND.EOFS(J1)<EOFS(J1+1)+EDIFFTOL*10.0D0) NSTRUCTREF=NSTRUCTREF+1
ENDDO

ALLOCATE(COORDSREF(NSTRUCTREF,3*NATOMS),VREF(NSTRUCTREF),NEAREST(NSTRUCTREF),REFSTRUCTFRAME(NSTRUCTREF))
ALLOCATE(DSTRUCT(NSTRUCTREF),DISTO(NSTRUCTREF),DIST(NSTRUCTREF))

NSTRUCTREF=1
NMIN=1
NTS=0
VREF(1)=EOFS(1) ! first minimum
REFSTRUCTFRAME(1)=1
COORDSREF(1,1:3*NATOMS)=PATHFRAMES(1,1:3*NATOMS)
HIGHESTREF=1
LASTMIN=.TRUE.
LASTTS=.FALSE.
PRINT '(A,I10,A,G20.10)','RPH> minimum assumed             for frame ',1,' energy=',VREF(NSTRUCTREF)
!
! Assume we will always have an order parameter. Need to replace the call below appropriately.
!
IF (.TRUE.) THEN
   DO K1=1,3
      DO K2=1,NATOMS
         QTEMP(K1,K2)=COORDSREF(1,3*(K2-1)+K1)
      ENDDO
   ENDDO
   CALL GETORDERRPH(QTEMP,QORDER,NATOMS)
   PRINT '(A,G20.10)','RPH> for first minimum Q=',QORDER
ENDIF

DO J1=2,NPATH-1
   IF (EOFS(J1-1)+EDIFFTOL*10.0D0<EOFS(J1).AND.EOFS(J1)>EOFS(J1+1)+EDIFFTOL*10.0D0) THEN
      NSTRUCTREF=NSTRUCTREF+1
      REFSTRUCTFRAME(NSTRUCTREF)=J1
      VREF(NSTRUCTREF)=EOFS(J1) 
      IF (VREF(NSTRUCTREF).GT.VREF(HIGHESTREF)) HIGHESTREF=NSTRUCTREF
      COORDSREF(NSTRUCTREF,1:3*NATOMS)=PATHFRAMES(J1,1:3*NATOMS)
      WRITE(*,*) "sn402: Changed RPH so that it calls ALIGN_DECIDE instead of MINPERMDIST"
      WRITE(*,*) "I haven't tested this change and am not certain whether it's sensible." 
      WRITE(*,*) "Please check carefully that this part of the code is working as you expect, then remove these messages!"
      CALL ALIGN_DECIDE(COORDSREF(NSTRUCTREF-1,1:3*NATOMS),COORDSREF(NSTRUCTREF,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
      PRINT '(A,I10,2(A,G20.10))','RPH> transition state identified for frame ',J1,' energy=', &
  &         VREF(NSTRUCTREF),' distance to previous minimum=',DUMMY
      NTS=NTS+1
      IF (.NOT.LASTMIN) THEN
         PRINT '(A)','RPH> WARNING *** previous stationary point identified was not a minimum'
         STOP
      ENDIF
      LASTMIN=.FALSE.
      LASTTS=.TRUE.
      IF (.TRUE.) THEN
         DO K1=1,3
            DO K2=1,NATOMS
               QTEMP(K1,K2)=COORDSREF(NSTRUCTREF,3*(K2-1)+K1)
            ENDDO
         ENDDO
         CALL GETORDER(QTEMP,QORDER,NATOMS)
         PRINT '(A,G20.10)','RPH> for this transition state Q=',QORDER
      ENDIF
   ENDIF
   IF ((EOFS(J1-1)>EOFS(J1).AND.EOFS(J1)<=EOFS(J1+1)).AND.(ABS(EOFS(J1)-VREF(NSTRUCTREF)).GT.EDIFFTOL)) THEN
      NSTRUCTREF=NSTRUCTREF+1
      REFSTRUCTFRAME(NSTRUCTREF)=J1
      VREF(NSTRUCTREF)=EOFS(J1) 
      IF (VREF(NSTRUCTREF).GT.VREF(HIGHESTREF)) HIGHESTREF=NSTRUCTREF
      COORDSREF(NSTRUCTREF,1:3*NATOMS)=PATHFRAMES(J1,1:3*NATOMS)
      WRITE(*,*) "sn402: Changed RPH so that it calls ALIGN_DECIDE instead of MINPERMDIST"
      WRITE(*,*) "I haven't tested this change and am not certain whether it's sensible." 
      WRITE(*,*) "Please check carefully that this part of the code is working as you expect, then remove these messages!"
      CALL ALIGN_DECIDE(COORDSREF(NSTRUCTREF-1,1:3*NATOMS),COORDSREF(NSTRUCTREF,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
      PRINT '(A,I10,2(A,G20.10))','RPH> minimum identified for frame ',J1,' energy=', &
  &         VREF(NSTRUCTREF),' distance to previous transition state=',DUMMY
      NMIN=NMIN+1
      IF (.NOT.LASTTS) THEN
         PRINT '(A)','RPH> WARNING *** previous stationary point identified was not a transition state'
      ENDIF
      LASTTS=.FALSE.
      LASTMIN=.TRUE.
      IF (.TRUE.) THEN
         DO K1=1,3
            DO K2=1,NATOMS
               QTEMP(K1,K2)=COORDSREF(NSTRUCTREF,3*(K2-1)+K1)
            ENDDO
         ENDDO
         CALL GETORDER(QTEMP,QORDER,NATOMS)
         PRINT '(A,G20.10)','RPH> for this minimum          Q=',QORDER
      ENDIF
   ENDIF
ENDDO
IF (NTS.EQ.1) THEN ! just one path - permutations should be OK already
   PERMDISTSAVE=.FALSE. 
   LPERMDISTSAVE=.FALSE.
   PERMDIST=.FALSE. 
   LPERMDIST=.FALSE.
ENDIF
IF (ABS(EOFS(NPATH)-VREF(NSTRUCTREF)).GT.EDIFFTOL) THEN
   NSTRUCTREF=NSTRUCTREF+1
   REFSTRUCTFRAME(NSTRUCTREF)=NPATH
   VREF(NSTRUCTREF)=EOFS(NPATH) ! last minimum
   COORDSREF(NSTRUCTREF,1:3*NATOMS)=PATHFRAMES(NPATH,1:3*NATOMS)
      WRITE(*,*) "sn402: Changed RPH so that it calls ALIGN_DECIDE instead of MINPERMDIST"
      WRITE(*,*) "I haven't tested this change and am not certain whether it's sensible." 
      WRITE(*,*) "Please check carefully that this part of the code is working as you expect, then remove these messages!"
      CALL ALIGN_DECIDE(COORDSREF(NSTRUCTREF-1,1:3*NATOMS),COORDSREF(NSTRUCTREF,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
      PRINT '(A,I10,2(A,G20.10))','RPH> minimum assumed    for frame ',J1,' energy=', &
  &         VREF(NSTRUCTREF),' distance to previous transition state=',DUMMY
   NMIN=NMIN+1
   IF (.NOT.LASTTS) THEN
      PRINT '(A)','RPH> WARNING *** previous stationary point identified was not a transition state'
   ENDIF
   IF (.TRUE.) THEN
      DO K1=1,3
         DO K2=1,NATOMS
            QTEMP(K1,K2)=COORDSREF(NSTRUCTREF,3*(K2-1)+K1)
         ENDDO
      ENDDO
      CALL GETORDER(QTEMP,QORDER,NATOMS)
      PRINT '(A,G20.10)','RPH> for this minimum          Q=',QORDER
   ENDIF
ENDIF
PRINT '(4(A,I10))','RPH> Total number of reference structures=',NSTRUCTREF,' with ',NTS,' ts and ',NMIN,' minima'
CALL FLUSH(6)
!
! work out total distance and slice width - use this for SLICELABEL
!
DISTTOTAL=0.0D0
      WRITE(*,*) "sn402: Changed RPH so that it calls ALIGN_DECIDE instead of MINPERMDIST"
      WRITE(*,*) "I haven't tested this change and am not certain whether it's sensible." 
      WRITE(*,*) "Please check carefully that this part of the code is working as you expect, then remove these messages!"
CALL ALIGN_DECIDE(COORDSREF(1,1:3*NATOMS),PATHFRAMES(1,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
PRINT '(A,F20.10)','RPH> first distance=',DUMMY
DISTTOTAL=DUMMY
DO J1=2,NPATH
      WRITE(*,*) "sn402: Changed RPH so that it calls ALIGN_DECIDE instead of MINPERMDIST"
      WRITE(*,*) "I haven't tested this change and am not certain whether it's sensible." 
      WRITE(*,*) "Please check carefully that this part of the code is working as you expect, then remove these messages!"
   CALL ALIGN_DECIDE(PATHFRAMES(J1-1,1:3*NATOMS),PATHFRAMES(J1,1:3*NATOMS),NATOMS,DEBUG, &
  &                 PARAM1,PARAM2,PARAM3,BULKT,TWOD,DUMMY,DIST2,RIGIDBODY,RMAT)
   PRINT '(A,I10,F20.10)','RPH> frame and distance=',J1,DUMMY
   DISTTOTAL=DISTTOTAL+DUMMY
ENDDO
!
! Need to resave COORDSREF to get consistent permutational alignment.
!
DO J1=1,NSTRUCTREF 
   COORDSREF(J1,1:3*NATOMS)=PATHFRAMES(REFSTRUCTFRAME(J1),1:3*NATOMS)
ENDDO
PERMDIST=.FALSE.
LPERMDIST=.FALSE. ! freeze permutations of references - should already be optimal now
!
! Calculate partition function at minimum 1 (pathframe 1) then
! find the two frames that bracket the next distances and calculate
! their partition functions, then do the final minimum.
!
! Check dot product of gradient with unit dislacement vector between
! bracking frames.
!
SLICEWIDTH=DISTTOTAL/RPHSLICES
DO J1=1,RPHSLICES
   SLICELABEL(J1)=SLICEWIDTH*(J1-0.5D0) ! these distances point to the middle of the distance slices
ENDDO
PRINT '(A,2G20.10)','RPH> Sum of path frame distances and slice width: ',DISTTOTAL,SLICEWIDTH
!
! This interval will depend upon the order parameter.
!
QORDERINT=(MAX(RPHQMAX,RPHQMIN)-MIN(RPHQMAX,RPHQMIN))/RPHQBINS
DO J1=1,RPHQBINS
!
! these distances point to the bottom of the Q bins
!
   BINLABELQORDER(J1)=MIN(RPHQMAX,RPHQMIN)+QORDERINT*(J1-1.0D0)
ENDDO

Q(1:3*NATOMS)=PATHFRAMES(1,1:3*NATOMS)
DO K1=1,3
   DO K2=1,NATOMS
      QTEMP(K1,K2)=Q(3*(K2-1)+K1)
   ENDDO
ENDDO
CALL GETORDER(QTEMP,QORDER,NATOMS)
PRINT '(A,G20.10)','RPH> Q=',QORDER
QORDERSLICES(0)=QORDER

IF (ENDNUMHESS) THEN
   CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
   CALL MAKENUMHESS(Q,NATOMS)
   WRITE (*,'(A,2G20.10)') ' geopt> Calculated numerical  hessian; energy and RMS=',ENERGY,RMS
ELSE
   CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
   WRITE (*,'(A,2G20.10)') ' geopt> Calculated analytical hessian; energy and RMS=',ENERGY,RMS
ENDIF

IF (CHRMMT) THEN
   CALL MASSWT2(NATOMS,ATMASS,Q,VNEW,.TRUE.)  ! just does the Hessian ???????????
ELSE
   CALL MASSWT(NATOMS,ATMASS,Q,VNEW,.TRUE.)   ! Hessian, coordinates and gradient vector ???????????
ENDIF
PROJGRAD=.FALSE. ! for first minimum
CALL PROJH(Q,NATOMS,ATMASS,VNEW,PROJGRAD)

CALL DSYEV('V','U',3*NATOMS,HESS,3*NATOMS,EVALUES,TEMPA,9*NATOMS,INFO)
IF (EVALUES(1).LT.EVALUES(3*NATOMS)) CALL EIGENSORT_VAL_ASC(EVALUES,HESS,3*NATOMS,3*NATOMS)
PRINT '(A)',' geopt> projected mass-weighted Hessian eigenvalues for first minimum:'
PRINT '(6G20.10)',EVALUES(1:3*NATOMS)
NEXMODES=6
IF (BULKT) NEXMODES=NEXMODES-3
IF (BULKT.AND.TWOD) NEXMODES=NATOMS+2
IF (PULLT.OR.EFIELDT) NEXMODES=4
IF (GTHOMSONT .AND. (GTHOMMET .EQ. 5) ) NEXMODES = 3
IF (GTHOMSONT .AND. (GTHOMMET < 5) ) NEXMODES = 1
IF (TWOD) NEXMODES=NEXMODES+NATOMS
IF (FREEZE) THEN
   NEXMODES=3*NFREEZE
ENDIF
IF (RBAAT) THEN
   IF (EFIELDT) THEN
      NEXMODES = 4
   ELSE
      NEXMODES = 6
   ENDIF
   IF (STOCKAAT) NEXMODES = NEXMODES + NATOMS/2
ENDIF
IF(MIEFT.OR.NOTRANSROTT) NEXMODES = 0
CALL GETPROD(PROD,NEXMODES,EVALUES)

KAPPA=3*NATOMS-6
IF (TWOD) THEN
   KAPPA=2*NATOMS-3
   IF (BULKT) KAPPA=2*NATOMS-2
      ELSE IF (PULLT) THEN
   KAPPA=3*NATOMS-4
ELSE IF (BULKT) THEN
   KAPPA=3*NATOMS-3
ELSE IF (MKTRAPT.OR.MIEFT.OR.NOTRANSROTT) THEN
   KAPPA=3*NATOMS
ELSE IF (EYTRAPT) THEN
   KAPPA=3*NATOMS-3
ELSE IF (FREEZE) THEN
   KAPPA=3*NATOMS-3*NFREEZE
ELSE IF (RIGIDINIT) THEN 
   KAPPA=DEGFREEDOMS-6
ELSE IF (UNRST) THEN
   KAPPA=NINTS
ENDIF

PFUNC(0)=-ENERGY/RPHTEMP - PROD/2.0D0 + KAPPA*LOG(RPHTEMP) ! - LOG(1.0D0*HORDERMIN(J1))
PFUNC(0)=PFUNC(0)+LOG(0.5D0) ! to allow for half the softest mode at the minimum
PSUM=PFUNC(0)

NEXMODES=NEXMODES+1
ALLOCATE(VSAVE(NOPT,NEXMODES))
!
! Now step through intervening slice geometries with additional projection along the gradient
! (or the frame displacement vector?).
!
DISTTOTAL2=0.0D0
ATFRAME=2
DO J1=1,RPHSLICES
   DO WHILE (DISTTOTAL2.LT.SLICELABEL(J1)) 
      WRITE(*,*) "sn402: Changed RPH so that it calls ALIGN_DECIDE instead of MINPERMDIST"
      WRITE(*,*) "I haven't tested this change and am not certain whether it's sensible." 
      WRITE(*,*) "Please check carefully that this part of the code is working as you expect, then remove these messages!"
      CALL ALIGN_DECIDE(PATHFRAMES(ATFRAME-1,1:3*NATOMS),PATHFRAMES(ATFRAME,1:3*NATOMS),NATOMS,DEBUG, &
  &                    PARAM1,PARAM2,PARAM3,BULKT,TWOD,DSAVE,DIST2,RIGIDBODY,RMAT)

      DISTTOTAL2=DISTTOTAL2+DSAVE
      PRINT '(A,2I10,2F20.10)','RPH> slice, frame, distance, and total=',J1,ATFRAME,DSAVE,DISTTOTAL2
      ATFRAME=ATFRAME+1
      IF (ATFRAME.GT.NPATH) PRINT '(A)','WARNING *** Reached last frame but still in slices loop'
   ENDDO
   FRACTION=(SLICELABEL(J1)-DISTTOTAL2+DSAVE)/DSAVE
   IF ((FRACTION.LT.0.0D0).OR.(FRACTION.GT.1.0)) THEN
      PRINT '(A,I6,A,I6,A,G20.10,A,2F20.10)','RPH> Fractional distance for frames ',ATFRAME,' and ',ATFRAME-1,' is ',FRACTION
      PRINT '(A,I6,3G20.10)','ERROR *** J1,SLICELABEL,DISTTOTAL2,DSAVE=',J1,SLICELABEL(J1),DISTTOTAL2,DSAVE
      STOP
   ENDIF
   PRINT '(A,I6,A,I6,A,G20.10,A,2F20.10)','RPH> Fractional distance for frames ',ATFRAME,' and ',ATFRAME-1,' is ',FRACTION, &
   & ' slice and frame distances: ',SLICELABEL(J1),DISTTOTAL2
   PRINT '(A,2I6)','ATFRAME,NPATH=',ATFRAME,NPATH

   Q(1:3*NATOMS)=PATHFRAMES(ATFRAME-2,1:3*NATOMS)+(PATHFRAMES(ATFRAME-2,1:3*NATOMS)-PATHFRAMES(ATFRAME-1,1:3*NATOMS))*FRACTION
   DO K1=1,3
      DO K2=1,NATOMS
         QTEMP(K1,K2)=Q(3*(K2-1)+K1)
      ENDDO
   ENDDO
   CALL GETORDER(QTEMP,QORDER,NATOMS)
   PRINT '(A,G20.10)','RPH> Q=',QORDER
   QORDERSLICES(J1)=QORDER

   IF (ENDNUMHESS) THEN
      CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
      CALL MAKENUMHESS(Q,NATOMS)
      WRITE (*,'(A,2G20.10)') ' geopt> Calculated numerical  hessian; energy and RMS=',ENERGY,RMS
   ELSE
      CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
      WRITE (*,'(A,2G20.10)') ' geopt> Calculated analytical hessian; energy and RMS=',ENERGY,RMS
   ENDIF
   SAVEHESS(1:NOPT,1:NOPT)=HESS(1:NOPT,1:NOPT)

   CALL MASSWT2(NATOMS,ATMASS,Q,VNEW,.TRUE.)  ! just does the Hessian

   CALL DSYEV('V','U',3*NATOMS,HESS,3*NATOMS,EVALUES,TEMPA,9*NATOMS,INFO)
   IF (EVALUES(1).LT.EVALUES(3*NATOMS)) CALL EIGENSORT_VAL_ASC(EVALUES,HESS,3*NATOMS,3*NATOMS)
   PRINT '(A,I6)',' geopt> unprojected non-mass-weighted Hessian eigenvalues for slice ',J1
   PRINT '(6G20.10)',EVALUES(1:3*NATOMS)
!
! Save the NEXMODES eigenvectors corresponding to the smallest eigenvalues.
!
   DO I1=NOPT,NOPT-NEXMODES+1,-1
      VSAVE(1:NOPT,NOPT-I1+1)=HESS(1:NOPT,I1)
   ENDDO

!  IF (ENDNUMHESS) THEN
!     CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!     CALL MAKENUMHESS(Q,NATOMS)
!     WRITE (*,'(A,2G20.10)') ' geopt> Calculated numerical  hessian; energy and RMS=',ENERGY,RMS
!  ELSE
!     CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
!     WRITE (*,'(A,2G20.10)') ' geopt> Calculated analytical hessian; energy and RMS=',ENERGY,RMS
!  ENDIF

   HESS(1:NOPT,1:NOPT)=SAVEHESS(1:NOPT,1:NOPT)

   CALL MASSWT2(NATOMS,ATMASS,Q,VNEW,.TRUE.)  
!
! Shift softest NEXMODES to zero (could use unity and take determinant to get eigenvalue product).
!
   DO I1=1,NEXMODES
      DO J2=1,NOPT
         DO J3=1,NOPT
            HESS(J3,J2)=HESS(J3,J2)-EVALUES(NOPT-I1+1)*VSAVE(J3,I1)*VSAVE(J2,I1)
         ENDDO
      ENDDO
   ENDDO

   DUMMY=0.0D0
   DO J2=1,3*NATOMS
      DUMMY=DUMMY+VNEW(J2)**2
   ENDDO
   DUMMY=1.0D0/SQRT(DUMMY)

   CALL DSYEV('V','U',3*NATOMS,HESS,3*NATOMS,EVALUES,TEMPA,9*NATOMS,INFO)
   IF (EVALUES(1).LT.EVALUES(3*NATOMS)) CALL EIGENSORT_VAL_ASC(EVALUES,HESS,3*NATOMS,3*NATOMS)
!  PRINT '(A,I6)',' geopt> trans/rot projected mass-weighted Hessian eigenvalues for slice ',J1
   PRINT '(A,I6)',' geopt> Shifted mass-weighted Hessian eigenvalues for slice ',J1
   PRINT '(6G20.10)',EVALUES(1:3*NATOMS)
   
!  DO J2=1,3*NATOMS
!     DUMMY2=0.0D0
!     DO J3=1,3*NATOMS
!        DUMMY2=DUMMY2+HESS(J3,J2)*VNEW(J3)*DUMMY
!     ENDDO
!     PRINT '(A,I6,G20.10)','RPU> Mode and dot product with normalised gradient vector: ',J2,DUMMY2
!  ENDDO
   DUMMY2=0.0D0
   CDUMMY(1:3*NATOMS)=PATHFRAMES(ATFRAME-2,1:3*NATOMS)-PATHFRAMES(ATFRAME-1,1:3*NATOMS)
   DO J2=1,3*NATOMS
      DUMMY2=DUMMY2+CDUMMY(J2)**2
   ENDDO
   DUMMY2=1.0D0/SQRT(DUMMY2)
   CDUMMY(1:3*NATOMS)=DUMMY2*CDUMMY(1:3*NATOMS)
   DUMMY2=0.0D0
   DO J2=1,3*NATOMS
      DUMMY2=DUMMY2+CDUMMY(J2)*VNEW(J2)*DUMMY
   ENDDO
   PRINT '(A,G20.10)','RPU> Dot product of normalised gradient and frame displacement vectors: ',DUMMY2

!   Q(1:3*NATOMS)=PATHFRAMES(ATFRAME-2,1:3*NATOMS)+(PATHFRAMES(ATFRAME-2,1:3*NATOMS)-PATHFRAMES(ATFRAME-1,1:3*NATOMS))*FRACTION
!   IF (ENDNUMHESS) THEN
!      CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
!      CALL MAKENUMHESS(Q,NATOMS)
!      WRITE (*,'(A,2G20.10)') ' geopt> Calculated numerical  hessian; energy and RMS=',ENERGY,RMS
!   ELSE
!      CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
!      WRITE (*,'(A,2G20.10)') ' geopt> Calculated analytical hessian; energy and RMS=',ENERGY,RMS
!   ENDIF
!
!   IF (CHRMMT) THEN
!      CALL MASSWT2(NATOMS,ATMASS,Q,VNEW,.TRUE.)  ! just does the Hessian ???????????
!   ELSE
!      CALL MASSWT(NATOMS,ATMASS,Q,VNEW,.TRUE.)   ! Hessian, coordinates and gradient vector ???????????
!   ENDIF
!   PROJGRAD=.TRUE.
!   IF (RMS.LT.CONVR) PROJGRAD=.FALSE.
!   CALL PROJH(Q,NATOMS,ATMASS,VNEW,PROJGRAD)
!
!   CALL DSYEV('V','U',3*NATOMS,HESS,3*NATOMS,EVALUES,TEMPA,9*NATOMS,INFO)
!   IF (EVALUES(1).LT.EVALUES(3*NATOMS)) CALL EIGENSORT_VAL_ASC(EVALUES,HESS,3*NATOMS,3*NATOMS)
!   PRINT '(A,I6)',' geopt> projected mass-weighted Hessian eigenvalues for slice ',J1
!   PRINT '(6G20.10)',EVALUES(1:3*NATOMS)
   CALL GETPROD(PROD,NEXMODES,EVALUES)
   PFUNC(J1)=-ENERGY/RPHTEMP - PROD/2.0D0 + (KAPPA-1)*LOG(RPHTEMP) ! - LOG(1.0D0*HORDERMIN(J1))
   PFUNC(J1)=PFUNC(J1)+LOG(SLICEWIDTH) ! quadrature over reaction coordinate
   PSUM=PSUM+PFUNC(J1)
ENDDO
DEALLOCATE(VSAVE)
!
! Final minimum
!
NEXMODES=NEXMODES-1
PRINT '(A,I6)','RPH> NEXMODES is ',NEXMODES
Q(1:3*NATOMS)=PATHFRAMES(NPATH,1:3*NATOMS)
DO K1=1,3
   DO K2=1,NATOMS
      QTEMP(K1,K2)=Q(3*(K2-1)+K1)
   ENDDO
ENDDO
CALL GETORDER(QTEMP,QORDER,NATOMS)
PRINT '(A,G20.10)','RPH> Q=',QORDER
QORDERSLICES(RPHSLICES+1)=QORDER
IF (ENDNUMHESS) THEN
   CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.FALSE.,RMS,.FALSE.,.FALSE.)
   CALL MAKENUMHESS(Q,NATOMS)
   WRITE (*,'(A,2G20.10)') ' geopt> Calculated numerical  hessian; energy and RMS=',ENERGY,RMS
ELSE
   CALL POTENTIAL(Q,ENERGY,VNEW,.TRUE.,.TRUE.,RMS,.FALSE.,.FALSE.)
   WRITE (*,'(A,2G20.10)') ' geopt> Calculated analytical hessian; energy and RMS=',ENERGY,RMS
ENDIF

IF (CHRMMT) THEN
   CALL MASSWT2(NATOMS,ATMASS,Q,VNEW,.TRUE.)  ! just does the Hessian ???????????
ELSE
   CALL MASSWT(NATOMS,ATMASS,Q,VNEW,.TRUE.)   ! Hessian, coordinates and gradient vector ???????????
ENDIF
PROJGRAD=.FALSE. ! for first minimum
CALL PROJH(Q,NATOMS,ATMASS,VNEW,PROJGRAD)

CALL DSYEV('V','U',3*NATOMS,HESS,3*NATOMS,EVALUES,TEMPA,9*NATOMS,INFO)
IF (EVALUES(1).LT.EVALUES(3*NATOMS)) CALL EIGENSORT_VAL_ASC(EVALUES,HESS,3*NATOMS,3*NATOMS)
PRINT '(A)',' geopt> projected mass-weighted Hessian eigenvalues for last minimum:'
PRINT '(6G20.10)',EVALUES(1:3*NATOMS)

CALL GETPROD(PROD,NEXMODES,EVALUES)
PFUNC(RPHSLICES+1)=-ENERGY/RPHTEMP - PROD/2.0D0 + KAPPA*LOG(RPHTEMP) ! - LOG(1.0D0*HORDERMIN(J1))
PFUNC(RPHSLICES+1)=PFUNC(RPHSLICES+1)+LOG(0.5D0) ! half the contribution of the softest mode
PSUM=PSUM+PFUNC(RPHSLICES+1)
PMEAN=PSUM/(RPHSLICES+2)
LUNIT=GETUNIT()
OPEN(LUNIT,FILE='RPH.Fofs',STATUS='UNKNOWN')
PRINT '(A,3G20.10)','PFUNC(0),PSUM,PMEAN=',PFUNC(0),PSUM,PMEAN

DUMMY=EXP(PFUNC(0)-PMEAN)
DO J1=1,RPHSLICES
   DUMMY=DUMMY+EXP(PFUNC(J1)-PMEAN)
ENDDO
DUMMY=DUMMY+EXP(PFUNC(RPHSLICES+1)-PMEAN)
DUMMY=1.0D0/DUMMY ! normalisation for the probability

DUMMY2=EXP(PFUNC(0)-PMEAN)*DUMMY
WRITE(LUNIT,'(I6,3G20.10)') 0,0.0D0,EXP(PFUNC(0)-PMEAN)*DUMMY,-RPHTEMP*(PFUNC(0)-PMEAN)
DO J1=1,RPHSLICES
   WRITE(LUNIT,'(I6,3G20.10)') J1,SLICELABEL(J1),EXP(PFUNC(J1)-PMEAN)*DUMMY,-RPHTEMP*(PFUNC(J1)-PMEAN)
   DUMMY2=DUMMY2+EXP(PFUNC(J1)-PMEAN)*DUMMY
ENDDO
WRITE(LUNIT,'(I6,3G20.10)') RPHSLICES+1,DISTTOTAL,EXP(PFUNC(RPHSLICES+1)-PMEAN)*DUMMY,-RPHTEMP*(PFUNC(RPHSLICES+1)-PMEAN)
DUMMY2=DUMMY2+EXP(PFUNC(RPHSLICES+1)-PMEAN)*DUMMY
CLOSE(LUNIT)
PRINT '(A,G20.10)','normalisation, DUMMY2=',DUMMY2

QORDERHIST(1:RPHQBINS)=1.0D-100
QORDERMEAN=0.0D0
QORDERSQMEAN=0.0D0
DO J1=0,RPHSLICES+1
   IBININDEX=INT((QORDERSLICES(J1)-MIN(RPHQMAX,RPHQMIN))/QORDERINT)+1
   PRINT '(A,I6,2G20.10,I6)','J1,QORDERSLICES,QORDERINT,IBININDEX=',J1,QORDERSLICES(J1),QORDERINT,IBININDEX
   PRINT '(A,I6,3G20.10)','J1,PFUNC,PMEAN,DUMMY=',J1,PFUNC(J1),PMEAN,DUMMY
   IF (IBININDEX.LE.RPHQBINS) QORDERHIST(IBININDEX)=QORDERHIST(IBININDEX)+EXP(PFUNC(J1)-PMEAN)*DUMMY
   PRINT '(A,I6,2G20.10)','J1,QORDERMEAN,QORDERSQMEAN=',J1,QORDERMEAN,QORDERSQMEAN
   QORDERMEAN=QORDERMEAN+QORDERSLICES(J1)*EXP(PFUNC(J1)-PMEAN)*DUMMY
   QORDERSQMEAN=QORDERSQMEAN+QORDERSLICES(J1)**2*EXP(PFUNC(J1)-PMEAN)*DUMMY
ENDDO

PRINT '(A,3G20.10)','RPH> <Q>, <Q^2>=',QORDERMEAN,QORDERSQMEAN
PRINT '(A,3G20.10)','RPH> <Q>, <Q^2> and sigma=',QORDERMEAN,QORDERSQMEAN,SQRT(QORDERSQMEAN-QORDERMEAN**2)

LUNIT=GETUNIT()
OPEN(LUNIT,FILE='FofQ',STATUS='UNKNOWN')
DO J1=1,RPHQBINS
   IF (QORDERHIST(J1).LT.2.0D-100) CYCLE
   WRITE(LUNIT,'(I6,4G20.10)') J1,BINLABELQORDER(J1),QORDERHIST(J1),-RPHTEMP*LOG(QORDERHIST(J1)),-LOG(QORDERHIST(J1))
ENDDO
CLOSE(LUNIT)

DEALLOCATE(COORDSREF,VREF,NEAREST,DSTRUCT)
DEALLOCATE(DISTO,DIST)
DEALLOCATE(EOFS,PATHFRAMES)

STOP

RETURN

END SUBROUTINE RPH

SUBROUTINE GETORDERRPH(QTEMP,QORDER,NATOMS)
USE KEY
IMPLICIT NONE
INTEGER NATOMS
DOUBLE PRECISION QTEMP(3,NATOMS),Q4,Q6,QORDER

IF (COLLAGENOP) THEN
    CALL COLLAGENOP_CALC(QTEMP,QORDER,NATOMS)
ELSE
    CALL QORDER_LJ(QTEMP,Q4,Q6)
    QORDER=Q6
ENDIF

END SUBROUTINE GETORDERRPH

SUBROUTINE GETPROD(PROD,NEXMODES,EVALUES)
USE KEY, ONLY : AMBERT, NABT, SDT, TTM3T, BOWMANT, AMBER12T
USE COMMONS, ONLY : NOPT, NATOMS
USE MODCHARMM, ONLY : CHRMMT
IMPLICIT NONE
DOUBLE PRECISION PROD, EVALUES(3*NATOMS)
INTEGER NEXMODES,I1
!
!  FES calculation here.
!
PROD=0.0D0
DO I1=1,NOPT-NEXMODES
   IF (I1.GT.1) THEN
      IF (EVALUES(I1-1).NE.0.0D0) THEN
         IF (ABS(EVALUES(I1)/EVALUES(I1-1)).LT.1.0D-2) THEN
            PRINT '(A,G20.10,A,G20.10)',' geopt> WARNING - decrease in magnitude of eigenvalues from ',EVALUES(I1-1), &
  &      ' to ',EVALUES(I1)
            PRINT '(A)',' geopt> WARNING - this could indicate a stationary point of the wrong index'
         ENDIF
      ENDIF
   ENDIF
   IF (EVALUES(I1).GT.0.0D0) THEN
      PROD=PROD+DLOG(EVALUES(I1))
   ELSE
      IF (I1.LT.(NOPT-NEXMODES)) PRINT *,'Higher order saddle detected: eigenvalue ',EVALUES(I1)
   ENDIF
ENDDO
IF (CHRMMT .OR. AMBERT .OR. AMBER12T .OR. NABT .OR. SDT .OR. TTM3T) THEN
!
! if charmm need to convert this to (radian/s)^2, rather than charmm unit
! conversion factor for this is 4.184 x 10^26
! same for AMBER and for Stillinger-David and TTM3.
!
   PROD=PROD+(NOPT-NEXMODES)*DLOG(4.184D26)
   WRITE (*,'(A,G20.10)') ' geopt> Scaling product of eigenvalues to SI units (radian/s)^2 by ', &
     &                                  (3*NATOMS-NEXMODES)*DLOG(4.184D26)
ENDIF
IF (BOWMANT) THEN
!
! If Bowman need to convert this to (radian/s)^2
! conversion factor for this is 2625.47 x 10^26
!
   PROD=PROD+(NOPT-NEXMODES)*DLOG(4.184D26)
   WRITE (*,'(A,G20.10)') ' geopt> Scaling product of eigenvalues to SI units (radian/s)^2 by ', &
     &                                  (3*NATOMS-NEXMODES)*DLOG(2625.47D26)
ENDIF
PRINT '(A,G20.10)','getprod> ln product of eigenvalues=',PROD

END SUBROUTINE GETPROD
