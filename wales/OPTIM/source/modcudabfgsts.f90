MODULE MODCUDABFGSTS

USE COMMONS, ONLY : NATOMS, DEBUG, NOPT, STPMAX, REDOPATH, IVEC

USE KEY, ONLY : CUDAPOT, GMAX, CONVU, CONVR, MINMAX, MXSTP, MAXXBFGS, XMAXERISE, CUDATIMET, CFUSIONT, COLDFUSIONLIMIT, XDGUESS, & 
                XMUPDATE, HINDEX, MASST, SHIFTED, CHECKINDEX, REVERSEUPHILLT, GRADSQ, FREEZE, READV, DUMPV, FIXD, EIGENONLY, & 
                NSECDIAG, EVPC, CEIG, NEVS, CHECKOVERLAPT, CHECKNEGATIVET, REOPT, BFGSSTEP, PUSHCUT, PUSHOFF, MAXMAX, TRAD, & 
                NBFGSMAX1, NBFGSMAX2, BFGSTSTOL, REPELTST, MUPDATE, MAXBFGS, MAXERISE, DGUESS, NOIT, VARIABLES, NOHESS, FROZEN, & 
                NFREEZE, PV, FIXAFTER, CPMD, PRINTPTS, TWOENDS, DUMPMAG, NSTEPMIN, AMBER12T, LJADD3T, NADDTARGET, LJADDREP, LJADDATT

USE GENRIGID, ONLY : RIGIDINIT, ATOMRIGIDCOORDT, DEGFREEDOMS, NRIGIDBODY, NSITEPERBODY, RIGIDGROUPS, MAXSITE, SITESRIGIDBODY, &
                     RIGIDSINGLES, IINVERSE, AACONVERGENCET

USE ZWK, ONLY : NUP, ZWORK

USE MODNEB, ONLY : NEWCONNECTT, NEWNEBT

USE INTCOMMONS, ONLY : BONDSFROMFILE

USE INTERNALS_WRAPPER, ONLY : INTWRAP_USEINTERNALS

USE MODAMBER9, ONLY : CHECKCISTRANSALWAYS, CHECKCISTRANSALWAYSDNA, CHECKCISTRANSALWAYSRNA

USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_INT, C_DOUBLE, C_BOOL, C_CHAR

IMPLICIT NONE

INTERFACE
    SUBROUTINE CUDA_BFGSTS(C_NATOMS, C_COORDS, C_CEIG, C_MFLAG, C_ENERGY, C_NEVS, ITMAX, C_ITER, C_MAXXBFGS, C_XMAXERISE, C_RMS, &
                          C_CUDAPOT, C_DEBUG, C_CUDATIMET, NCALLS, C_CFUSIONT, C_COLDFUSIONLIMIT, C_XDGUESS, C_XMUPDATE, &
                          C_EVALMIN, C_VECS, C_ATOMRIGIDCOORDT, C_DEGFREEDOMS, C_NRIGIDBODY, C_NSITEPERBODY, C_RIGIDGROUPS, &
                          C_MAXSITE, C_SITESRIGIDBODY, C_RIGIDSINGLES, C_IINVERSE, C_NSECDIAG, C_EVPC, &
                          C_PUSHCUT, C_PUSHOFF, C_STPMAX, C_MAXMAX, C_MINMAX, C_TRAD, C_NBFGSMAX1, C_NBFGSMAX2, C_BFGSTSTOL, &
                          C_MUPDATE, C_GMAX, C_MAXBFGS, C_MAXERISE, C_DGUESS, C_CONVU, C_FREEZE, C_FROZEN, & 
                          C_NFREEZE, POTENTIALTIME, C_AACONVERGENCET, C_NADDTARGET, C_LJADDREP, & 
                          C_LJADDATT)BIND(C,NAME="setup_bfgsts")

        IMPORT :: C_INT, C_DOUBLE, C_BOOL, C_CHAR

        INTEGER(KIND=C_INT), INTENT(IN) :: C_NATOMS, & ! No. of atoms
                                           ITMAX, & ! Max. no. of BFGSTS iterations allowed
                                           C_DEGFREEDOMS, & ! Rigid Body Framework (RBF): no. of degrees of freedom
                                           C_NRIGIDBODY, & ! RBF: no. of rigid bodies
                                           C_MAXSITE, & ! RBF: max. no. of sites in a rigid body
                                           C_XMUPDATE, & ! History size for Rayleigh-Ritz (R-R) LBFGS minimization
                                           C_NSECDIAG, & ! Method used to approximate eigenvalue in R-R LBFGS
                                           C_NEVS, & ! Max. no. of steps allowed in R-R LBFGS 
                                           C_NBFGSMAX1, & ! No. of LBFGS steps allowed in restricted minimization (grad projected)
                                           C_NBFGSMAX2, & ! No. of LBFGS steps allowed in tangent space after eigenvalue converged
                                           C_MUPDATE, & ! History size for restricted LBFGS minimization
                                           C_NFREEZE, & ! No. of frozen atoms
                                           C_NADDTARGET ! Target cluster size (addressability)

        INTEGER(KIND=C_INT), DIMENSION(C_NRIGIDBODY), INTENT(IN) :: C_NSITEPERBODY ! RBF: no. of rigid body sites

        INTEGER(KIND=C_INT), DIMENSION(C_MAXSITE*C_NRIGIDBODY), INTENT(IN) :: C_RIGIDGROUPS ! RBF: list of atoms in rigid bodies

        INTEGER(KIND=C_INT), DIMENSION(C_DEGFREEDOMS/3 - 2*C_NRIGIDBODY), INTENT(IN) :: C_RIGIDSINGLES ! RBF: list of atoms not in rigid bodies

        INTEGER(KIND=C_INT), INTENT(OUT) :: C_ITER, & ! No. of BFGSTS iterations done
                                            NCALLS ! No. of potential calls made during this call to BFGSTS 

        REAL(KIND=C_DOUBLE), INTENT(IN) :: C_CEIG, & ! Convergence tolerance for RMS force for R-R LBFGS
                                           C_MAXXBFGS, & ! Max. step size allowed for R-R LBFGS
                                           C_XMAXERISE, & ! Max. energy rise allowed for R-R LBFGS
                                           C_COLDFUSIONLIMIT, & ! Limit below which cold fusion is diagnosed and BFGSTS terminated
                                           C_XDGUESS, & ! Initial guess for R-R inverse Hessian diagonal elements
                                           C_EVPC, & ! Convergence condition for abs(% change) in eigenvalue for R-R LBFGS
                                           C_PUSHCUT, & ! RMS threshold under which a pushoff will be applied
                                           C_PUSHOFF, & ! Magnitude of step away from transition state
                                           C_STPMAX, & ! Max. step along the eigenvector
                                           C_MAXMAX, & ! Max. allowed uphill step size
                                           C_MINMAX,  & ! Min. value that max. uphill step size is allowed to fall to
                                           C_TRAD, & ! Trust radius for adjusting max. step along eigenvector
                                           C_BFGSTSTOL, & ! Tolerance for eigenvector overlap (no. tangent space steps small->large)
                                           C_GMAX, & ! Convergence tolerance for RMS force for restricted LBFGS minimization
                                           C_MAXBFGS, & ! Max. step size allowed for restricted LBFGS
                                           C_MAXERISE, & ! Max. energy rise allowed for restricted LBFGS
                                           C_DGUESS, & ! Initial guess for restricted LBFGS inverse Hessian diagonal elements
                                           C_CONVU ! Can only converge if eigenvector following step length magnitude falls below this value

        REAL(KIND=C_DOUBLE), DIMENSION(C_MAXSITE*3*C_NRIGIDBODY), INTENT(IN) :: C_SITESRIGIDBODY ! RBF: coordinates of the rigid body sites

        REAL(KIND=C_DOUBLE), DIMENSION(C_NRIGIDBODY*3*3), INTENT(IN) :: C_IINVERSE ! RBF: inverse eigenvalues of the unweighted tensor of gyration

        REAL(KIND=C_DOUBLE), DIMENSION(C_NADDTARGET, C_NADDTARGET), INTENT(IN) ::C_LJADDREP, C_LJADDATT ! Repulsive/attractive epsilon matrix

        REAL(KIND=C_DOUBLE), DIMENSION(3*C_NATOMS), INTENT(INOUT) :: C_COORDS, & ! Coordinates
                                                                     C_VECS ! Lowest eigenvector

        REAL(KIND=C_DOUBLE), INTENT(OUT) :: C_ENERGY, & ! Energy at coordinates C_COORDS
                                            C_RMS, & ! RMS force
                                            C_EVALMIN, & ! Smallest eigenvalue
                                            POTENTIALTIME ! Time taken in calculating potential

        LOGICAL(KIND=C_BOOL), INTENT(IN) :: C_DEBUG, & ! If true, print debug info.
                                            C_CUDATIMET, & ! If true, print timing info.
                                            C_ATOMRIGIDCOORDT, & ! If false, use rigid body coordinates
                                            C_FREEZE, & ! If true, freeze some specified coordinates
                                            C_AACONVERGENCET ! If true, use more accurate method of calculating rigid RMS force

        LOGICAL(KIND=C_BOOL), DIMENSION(C_NATOMS), INTENT(IN) :: C_FROZEN ! Logical array specifying frozen atoms

        LOGICAL(KIND=C_BOOL), INTENT(OUT) :: C_MFLAG, & ! True if converged to a transition state
                                             C_CFUSIONT ! Set to true during BFGSTS if cold fusion diagnosed

        CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN) :: C_CUDAPOT ! Character specifying the CUDA potential to be used

    END SUBROUTINE CUDA_BFGSTS
END INTERFACE
 
CONTAINS

    SUBROUTINE CUDA_BFGSTS_WRAPPER(ITMAX,COORDS,ENERGY,MFLAG,RMS,EVALMIN,VECS,ITER)

        ! Variables passed as *arguments through this wrapper* (not common) with intent in for CUDA_BFGSTS are converted directly
        INTEGER(KIND=C_INT) :: C_NATOMS, ITMAX, C_ITER, NCALLS, C_DEGFREEDOMS, C_NRIGIDBODY, C_MAXSITE, C_XMUPDATE, &
                               C_NSECDIAG, C_NEVS, C_NBFGSMAX1, C_NBFGSMAX2, C_MUPDATE, C_NFREEZE, C_NADDTARGET
        INTEGER(KIND=C_INT), DIMENSION(NRIGIDBODY) :: C_NSITEPERBODY
        INTEGER(KIND=C_INT), DIMENSION(MAXSITE*NRIGIDBODY) :: C_RIGIDGROUPS
        INTEGER(KIND=C_INT), DIMENSION(DEGFREEDOMS/3 - 2*NRIGIDBODY) :: C_RIGIDSINGLES

        REAL(KIND=C_DOUBLE) :: C_CEIG, C_ENERGY, C_RMS, C_XDGUESS, C_EVALMIN, C_MAXXBFGS, C_XMAXERISE, C_COLDFUSIONLIMIT, &
                               C_EVPC, C_PUSHCUT, C_PUSHOFF, C_STPMAX, C_MAXMAX, C_MINMAX, C_TRAD, C_BFGSTSTOL, C_GMAX, &
                               C_MAXBFGS, C_MAXERISE, C_DGUESS, C_CONVU, POTENTIALTIME
        REAL(KIND=C_DOUBLE), DIMENSION(MAXSITE*3*NRIGIDBODY) :: C_SITESRIGIDBODY 
        REAL(KIND=C_DOUBLE), DIMENSION(NRIGIDBODY*3*3) :: C_IINVERSE
        REAL(KIND=C_DOUBLE), DIMENSION(NADDTARGET*NADDTARGET) :: C_LJADDREP, C_LJADDATT
        REAL(KIND=C_DOUBLE), DIMENSION(3*NATOMS) :: C_COORDS, C_VECS

        LOGICAL(KIND=C_BOOL) :: C_MFLAG, C_DEBUG, C_CUDATIMET, C_CFUSIONT, C_ATOMRIGIDCOORDT, C_FREEZE, C_AACONVERGENCET
        LOGICAL(KIND=C_BOOL), DIMENSION(NATOMS) :: C_FROZEN

        CHARACTER(LEN=1, KIND=C_CHAR) :: C_CUDAPOT

        ! Variables passed as *arguments through this wrapper* (not common) with intent out for CUDA_BFGSTS are not passed into it
        ! Therefore uninitialised C types are passed in and converted types are copied back after the call

        INTEGER :: I, J, K, ITER, FCALL, NPCALL, ECALL, SCALL
        DOUBLE PRECISION :: ENERGY, RMS, EVALMIN, FTIME, ETIME, STIME
        DOUBLE PRECISION, DIMENSION(3*NATOMS) :: COORDS, VECS
        LOGICAL :: MFLAG, KNOWE, KNOWG, KNOWH, DRAGT, DOINTERNALSTRANSFORM
        COMMON /PCALL/ NPCALL, ECALL, FCALL, SCALL, ETIME, FTIME, STIME
        COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
        COMMON /RUNTYPE/ DRAGT

        IF (CHECKCISTRANSALWAYSDNA .OR. CHECKCISTRANSALWAYSRNA) THEN
            WRITE(*,'(A)') " modcudabfgsts> Warning: Subroutine CHIRALITY_CHECK not implemented for BFGSTS on GPU. "
        END IF

        IF (CHECKCISTRANSALWAYS) THEN
            WRITE(*,'(A)') " modcudabfgsts> Warning: Subroutine CIS_TRANS_CHECK not implemented for BFGSTS on GPU. "
        END IF

        IF (.NOT. NOHESS) THEN
            WRITE(*,'(A)') " modcudabfgsts> Keyword NOHESS must be used. "
            STOP
        END IF

        DOINTERNALSTRANSFORM = INTWRAP_USEINTERNALS()

        IF (MASST .OR. CHECKINDEX .OR. REVERSEUPHILLT .OR. GRADSQ .OR. READV .OR. DUMPV .OR. &
            EIGENONLY .OR. FIXD .OR. CHECKOVERLAPT .OR. CHECKNEGATIVET .OR. REDOPATH .OR. REOPT.OR.BFGSSTEP.OR. &
            REPELTST .OR. NOIT .OR. VARIABLES .OR. (FIXAFTER .GT. 0) .OR. PV .OR. (IVEC .NE. 0) .OR. CPMD .OR. DRAGT .OR. & 
            TWOENDS .OR. DOINTERNALSTRANSFORM .OR. BONDSFROMFILE .OR. DUMPMAG .OR. (NSTEPMIN .GT. 0)) THEN
            WRITE(*,'(A)') " modcudabfgsts> Keywords MASS, CHECKINDEX, REVERSEUPHILLT, GRADSQ, READVEC, &
                                 DUMPVECTOR, FIXD, EIGENONLY, CHECKOVERLAP, CHECKNEGATIVE, REDOPATH, REOPT, BFGSSTEP, REPELTST, &
                                 NOIT, VARIABLES, FIXAFTER, PV, MODE, MODEDOWN, CPMD, DRAG, TWOENDS, INTMIN, BONDSFROMFILE, & 
                                 DUMPMAG, STEPMIN are not yet supported. & 
                                 Contact rgm38 if you would like a feature to be added. "
            WRITE(*,'(A)') " modcudabfgsts> Disclaimer: this list might not be exhaustive! "
            STOP
        END IF

        IF (PRINTPTS) THEN
            WRITE(*,'(A)') " modcudabfgsts> Keyword NOPOINTS must be used. & 
                                 Functionality to print coordinates/energies will be added at some point in the future. "
            STOP
        END IF

        IF (SHIFTED) THEN
            WRITE(*,'(A)') " modcudabfgsts> Zeros are shifted in lowest eigenvector. "
            STOP
        END IF

        ! Check for consistent convergence criteria on RMS gradient in LBFGS and EF part
        IF (GMAX.NE.CONVR) THEN
            IF (DEBUG) THEN 
                PRINT '(2(A,G20.10),A)', " modcudabfgsts> WARNING - GMAX ", GMAX, " is different from CONVR ", CONVR, " - resetting. "
            END IF
            GMAX=MIN(CONVR,GMAX)
            CONVR=MIN(CONVR,GMAX)
        ENDIF
        IF (MINMAX.GE.CONVU) THEN
            IF (DEBUG) THEN 
                PRINT '(2(A,G20.10),A)', " modcudabfgsts> WARNING - CONVU <= MINMAX - resetting CONVU. "
            END IF
            CONVU=2*MINMAX
        ENDIF

        !  Reset maximum step sizes in case this isn't the first call
        IF (DEBUG) THEN 
            PRINT '(A,G20.10)', " modcudabfgsts> resetting maximum step sizes to ", MXSTP
        END IF
        DO I=1,NOPT
            STPMAX(I)=MXSTP
        ENDDO

        NUP=HINDEX

        ! Variables from common blocks or modules with intent in or inout are copied into C types

        ! CUDA_BFGSTS currently only supports finding one eigenvector using Rayleigh-Ritz minimization
        IF ((SIZE(ZWORK) .NE. (3*NATOMS)) .OR. (NUP .NE. 1)) THEN
            WRITE(*,'(A)') " modcudabfgsts> HINDEX>1 is not supported as second derivatives are not available for AMBER with CUDA. " 
            STOP
        END IF

        DO K = 1,NRIGIDBODY
            DO J = 1,3
                DO I = 1,MAXSITE
                     C_SITESRIGIDBODY((K - 1)*3*MAXSITE + (J - 1)*MAXSITE + I) = SITESRIGIDBODY(I,J,K)
                END DO
            END DO
        END DO

        DO J = 1,NRIGIDBODY
            DO I = 1,MAXSITE
                C_RIGIDGROUPS((J - 1)*MAXSITE + I) = RIGIDGROUPS(I,J)
            END DO
        END DO

        DO I = 1,NRIGIDBODY
            C_NSITEPERBODY(I) = NSITEPERBODY(I)
        END DO

        DO I = 1,(DEGFREEDOMS/3 - 2*NRIGIDBODY)
            C_RIGIDSINGLES(I) = RIGIDSINGLES(I)
        END DO

        DO K = 1,3
            DO J = 1,3
                DO I = 1,NRIGIDBODY
                   C_IINVERSE((K - 1)*3*NRIGIDBODY + (J - 1)*NRIGIDBODY + I) = IINVERSE(I,J,K)
                END DO
            END DO
        END DO

        IF (ALLOCATED(LJADDREP) .AND. ALLOCATED(LJADDATT) .AND. LJADD3T) THEN
            DO J = 1,NADDTARGET
                DO I = 1,NADDTARGET
                    C_LJADDREP((J - 1)*NADDTARGET + I) = LJADDREP(I,J)
                    C_LJADDATT((J - 1)*NADDTARGET + I) = LJADDATT(I,J)
                END DO
            END DO
        ELSE
            C_LJADDREP(:) = 1.0
            C_LJADDATT(:) = 1.0
        END IF

        DO I = 1,3*NATOMS
            C_COORDS(I) = COORDS(I)
        END DO

        DO I = 1,NATOMS
            C_FROZEN(I) = FROZEN(I)
        END DO

        C_NATOMS = NATOMS
        C_CEIG = CEIG
        C_MAXXBFGS = MAXXBFGS
        C_XMAXERISE = XMAXERISE
        C_DEBUG = DEBUG
        C_CUDAPOT = CUDAPOT
        C_CUDATIMET = CUDATIMET
        C_ATOMRIGIDCOORDT = ATOMRIGIDCOORDT
        C_DEGFREEDOMS = DEGFREEDOMS
        C_NRIGIDBODY = NRIGIDBODY
        C_MAXSITE = MAXSITE
        C_COLDFUSIONLIMIT = COLDFUSIONLIMIT
        C_XDGUESS = XDGUESS
        C_XMUPDATE = XMUPDATE
        C_VECS = VECS
        C_NSECDIAG = NSECDIAG
        C_EVPC = EVPC
        C_NEVS = NEVS
        C_PUSHCUT = PUSHCUT
        C_PUSHOFF = PUSHOFF
        C_STPMAX = STPMAX(1)
        C_MAXMAX = MAXMAX
        C_MINMAX = MINMAX
        C_TRAD = TRAD
        C_NBFGSMAX1 = NBFGSMAX1
        C_NBFGSMAX2 = NBFGSMAX2
        C_BFGSTSTOL = BFGSTSTOL
        C_MUPDATE = MUPDATE
        C_GMAX = GMAX
        C_MAXBFGS = MAXBFGS
        C_MAXERISE = MAXERISE
        C_DGUESS = DGUESS
        C_CONVU = CONVU
        C_FREEZE = FREEZE
        C_NFREEZE = NFREEZE
        C_AACONVERGENCET = AACONVERGENCET
        C_NADDTARGET = NADDTARGET

        ! 'C_' prefix denotes those variables which have intent out or inout or are copies of those from common blocks/modules
        CALL CUDA_BFGSTS(C_NATOMS, C_COORDS, C_CEIG, C_MFLAG, C_ENERGY, C_NEVS, ITMAX, C_ITER, C_MAXXBFGS, C_XMAXERISE, C_RMS, &
                        C_CUDAPOT, C_DEBUG, C_CUDATIMET, NCALLS, C_CFUSIONT, C_COLDFUSIONLIMIT, C_XDGUESS, C_XMUPDATE, C_EVALMIN, &
                        C_VECS, C_ATOMRIGIDCOORDT, C_DEGFREEDOMS, C_NRIGIDBODY, C_NSITEPERBODY, C_RIGIDGROUPS, C_MAXSITE, &
                        C_SITESRIGIDBODY, C_RIGIDSINGLES, C_IINVERSE, C_NSECDIAG, C_EVPC, C_PUSHCUT, C_PUSHOFF, &
                        C_STPMAX, C_MAXMAX, C_MINMAX, C_TRAD, C_NBFGSMAX1, C_NBFGSMAX2, C_BFGSTSTOL, C_MUPDATE, C_GMAX, &
                        C_MAXBFGS, C_MAXERISE, C_DGUESS, C_CONVU, C_FREEZE, C_FROZEN, C_NFREEZE, POTENTIALTIME, C_AACONVERGENCET, &
                        C_NADDTARGET, C_LJADDREP, C_LJADDATT)

        ! Make sure C types with intent out or inout are coverted back to Fortran ones

        DO I = 1,3*NATOMS
            COORDS(I) = DBLE(C_COORDS(I))
        END DO

        DO I = 1,3*NATOMS
            VECS(I) = DBLE(C_VECS(I))
        END DO

        DO I = 1,3*NATOMS
            ZWORK(I,1) = DBLE(C_VECS(I))
        END DO

        ENERGY = DBLE(C_ENERGY)
        MFLAG = LOGICAL(C_MFLAG)
        RMS = DBLE(C_RMS)
        CFUSIONT = LOGICAL(C_CFUSIONT)
        EVALMIN = DBLE(C_EVALMIN)
        ITER = INT(C_ITER)

        FCALL = FCALL + INT(NCALLS)
        FTIME = FTIME + DBLE(POTENTIALTIME)

        KNOWE=.FALSE.
        KNOWG=.FALSE.
        KNOWH=.FALSE.

    END SUBROUTINE CUDA_BFGSTS_WRAPPER

END MODULE MODCUDABFGSTS
