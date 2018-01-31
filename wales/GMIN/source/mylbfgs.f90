! This subroutine is to select whether we are to minimise in rigid body/atom coordinates
! ATOMRIGIDCOORDT is the toggling logical variable
SUBROUTINE MYLBFGS(N, M, XCOORDS, DIAGCO, EPS, MFLAG, ENERGY, ITMAX, ITDONE, RESET, NP)
   USE COMMONS, ONLY: NATOMS, DMACRYST, HYBRIDMINT, CUDAT, AMBER12T, EPSRIGID, &
                      DEBUG, MYUNIT, GRADPROBLEMT, RMS, QCIPOTT
   USE GENRIGID, ONLY: DEGFREEDOMS, RIGIDINIT, ATOMRIGIDCOORDT, TRANSFORMCTORIGID, &
                       TRANSFORMRIGIDTOC, NRIGIDBODY
   USE MODCUDALBFGS, ONLY: CUDA_LBFGS_WRAPPER
!  USE PREC, ONLY: INT32, REAL64
   USE RAD_MOD, ONLY: RADCOM
   IMPLICIT NONE
! Arguments
! jk669 3/3/17 Intents added.
   INTEGER, intent(in) :: N
   INTEGER, intent(in) :: M
   DOUBLE PRECISION,intent(inout) :: XCOORDS(3*NATOMS)
   LOGICAL, intent(in) :: DIAGCO
   DOUBLE PRECISION, intent(in) :: EPS
   LOGICAL, intent(inout) :: MFLAG
   DOUBLE PRECISION, intent(inout) :: ENERGY
   INTEGER, intent(in) :: ITMAX
   INTEGER, intent(inout) :: ITDONE
   LOGICAL, intent(inout) :: RESET
   INTEGER, intent(in) :: NP
! Variables
   LOGICAL           :: EVAP, EVAPREJECT
! Common block
! TODO: delete common block
! ds656> EVAPREJECT was missing here until 10/06/2015
   COMMON /EV/ EVAP, EVAPREJECT
!   INTEGER           :: TEMPNATOMS, TEMPN !Commented out to match commented out code. jk669 3/3/17
   DOUBLE PRECISION      :: XRIGIDCOORDS(DEGFREEDOMS)
   DOUBLE PRECISION      :: EPS_TEMP

! vr274> DMACRYS uses a different minimization protocol to avoid cold fusions
   IF (DMACRYST) THEN
      CALL DMACRYS_MINIMIZE(N, M, XCOORDS, DIAGCO, EPS, MFLAG, ENERGY, ITMAX, ITDONE, RESET, NP)
      RETURN
   END IF

! hk286 > if generalised rigid body is used, then use rigid coords, otherwise proceed as usual
! hk286 > When passing rigid coords, pass zeroes from (DEGFREEDOMS+1:3*NATOMS)
! hk286 > mymylbfgs is the old mylbfgs
   IF (RIGIDINIT) THEN
      ATOMRIGIDCOORDT = .FALSE.
! Convert to rigid body coordinates
      CALL TRANSFORMCTORIGID(XCOORDS, XRIGIDCOORDS)
      XCOORDS(1:DEGFREEDOMS) = XRIGIDCOORDS(1:DEGFREEDOMS)
      XCOORDS(DEGFREEDOMS+1:3*NATOMS) = 0.0D0
! csw34> If HYBRIDMIN used, we want to minimise to EPSRIGID using rigid coords and then to EPS
!        using atomic coordinates
      IF (HYBRIDMINT) THEN
         EPS_TEMP = EPSRIGID
      ELSE
         EPS_TEMP = EPS
      END IF
! Call to CUDA LBFGS for rigid minimisation to tolerance epsrigid before the atomistic minimisation
      IF (CUDAT) THEN
         IF (.NOT. (AMBER12T)) THEN
            CALL RADCOM(XCOORDS, .TRUE.) ! Evaporated atoms moved back in
         END IF
         CALL CUDA_LBFGS_WRAPPER(N, M, XCOORDS, EPS_TEMP, MFLAG, ENERGY, ITMAX, ITDONE, RESET)
         IF (.NOT. (AMBER12T)) THEN
            EVAP = .FALSE.
            CALL RADCOM(XCOORDS, .FALSE.) ! Evaporated atoms NOT moved back in
            IF (EVAP) THEN
               MFLAG = .FALSE.
            END IF
         END IF
      ELSE
         !CALL MYMYLBFGS(N, M, XCOORDS, DIAGCO, EPS_TEMP, MFLAG, ENERGY, ITMAX, ITDONE, RESET, NP) !replaced by wrapper jk669 
         call min_wrapper(n,m,xcoords,diagco,eps,mflag,energy,itmax,itdone,reset,np)
      END IF
      IF (DEBUG.AND.HYBRIDMINT) WRITE(MYUNIT, '(A)') ' HYBRIDMIN> Rigid body minimisation converged, switching to all-atom'
! Convert back to atomistic coordinates. 
      XRIGIDCOORDS(1:DEGFREEDOMS) = XCOORDS(1:DEGFREEDOMS)
      CALL TRANSFORMRIGIDTOC(1, NRIGIDBODY, XCOORDS, XRIGIDCOORDS)
      ATOMRIGIDCOORDT = .TRUE.
   END IF

! If we're not using rigid bodies, or we've already completed a rigid body minimisation and we're using
! hybrid minimisation, then perform a minimisation in atomistic coordinates.
   IF ((.NOT. RIGIDINIT) .OR. HYBRIDMINT) THEN
      IF (CUDAT) THEN
! Call to CUDA LBFGS for rigid minimisation
         IF (.NOT. (AMBER12T)) THEN
            CALL RADCOM(XCOORDS, .TRUE.) ! Evaporated atoms moved back in
         END IF
         CALL CUDA_LBFGS_WRAPPER(N, M, XCOORDS, EPS, MFLAG, ENERGY, ITMAX, ITDONE, RESET)
         IF (.NOT. (AMBER12T)) THEN
            EVAP = .FALSE.
            CALL RADCOM(XCOORDS, .FALSE.) ! Evaporated atoms NOT moved back in
            IF (EVAP) THEN
               MFLAG = .FALSE.
            END IF
         END IF
      ELSE
         !CALL MYMYLBFGS(N, M, XCOORDS, DIAGCO, EPS_TEMP, MFLAG, ENERGY, ITMAX, ITDONE, RESET, NP)
         call min_wrapper(n,m,xcoords,diagco,eps,mflag,energy,itmax,itdone,reset,np) !replaced by wrapper jk669 
      END IF
   END IF

   GRADPROBLEMT = .FALSE.
   IF (RMS.EQ.1.0D-100) THEN
       IF (QCIPOTT) THEN
!         WRITE(MYUNIT,'(A)') ' mylbfgs> WARNING - RMS force was set to 1.0e-100'
       ELSE
          WRITE(MYUNIT,'(2A)') ' mylbfgs> WARNING - RMS force was set to 1.0e-100 as it was either smaller ', &
            'than this value or NaN - discarding structure'
          WRITE(MYUNIT,'(A)') '  - see debug output for further information'
          GRADPROBLEMT = .TRUE.
          MFLAG=.FALSE.
      ENDIF
   ENDIF

! hk286 - alternative implementation, fool mymylbfgs the number of atoms rather than passing zeroes
!      IF (RIGIDINIT) THEN
!         ATOMRIGIDCOORDT = .FALSE.
!         CALL TRANSFORMCTORIGID(XCOORDS, XRIGIDCOORDS)
!         TEMPN = N
!         N = DEGFREEDOMS
!         TEMPNATOMS = NATOMS
!         NATOMS  = DEGFREEDOMS / 3
!         CALL MYMYLBFGS(N,M,XRIGIDCOORDS,DIAGCO,EPS,MFLAG,ENERGY,ITMAX,ITDONE,RESET,NP,TEMPNATOMS)
!         NATOMS = TEMPNATOMS
!         N = TEMPN
!         CALL TRANSFORMRIGIDTOC(1,NRIGIDBODY, XCOORDS, XRIGIDCOORDS)
!         ATOMRIGIDCOORDT = .TRUE.      
!      ELSE
!         ATOMRIGIDCOORDT = .TRUE.
!        CALL MYMYLBFGS(N,M,XCOORDS,DIAGCO,EPS,MFLAG,ENERGY,ITMAX,ITDONE,RESET,NP)
!      END IF
! hk286

END SUBROUTINE MYLBFGS

!jk669: wrapper function for selecting minimization algorithm. 3/3/17
!Default is LBFGS without line search unless SQNM keyword is present in the data file. 
subroutine min_wrapper(numcoords,mupdates,xcoords,diagco,maxrmsgrad,converget,energy,itermax,iter,reset,np)
use COMMONS, only: SQNM_HISTMAX, SQNMT
! use PREC, only: INT32, REAL64
implicit none 
   !all variables are from mylbfgs except those from COMMONS. 
   integer, intent(in) :: numcoords
   integer, intent(in) :: mupdates
   DOUBLE PRECISION, intent(inout) :: xcoords(numcoords)
   logical, intent(in) :: diagco
   DOUBLE PRECISION, intent(in) :: maxrmsgrad
   logical, intent(inout) :: converget
   DOUBLE PRECISION, intent(out) :: energy
   integer, intent(in) :: itermax
   integer, intent(inout) :: iter
   logical, intent(inout) :: reset
   integer, intent(in) :: np

if (SQNMT) then !SQNM minimizer is turned on. 
   call sqnm(numcoords,xcoords,maxrmsgrad,SQNM_HISTMAX,converget,energy,itermax,iter)
else
   call MYMYLBFGS(numcoords,mupdates,xcoords,diagco,maxrmsgrad,converget,energy,itermax,iter,reset,np)
end if 

end subroutine min_wrapper
