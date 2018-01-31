MODULE MODCUDALBFGS

USE KEY, ONLY : GMAX, CUDATIMET, CUDAPOT, MAXBFGS, MAXERISE, CFUSIONT, COLDFUSIONLIMIT, DGUESS, GRADSQ, FREEZE, REVERSEUPHILLT, & 
                FREEZE, FROZEN, NFREEZE, PRINTPTS, TWOENDS, PV, DUMPMAG, NSTEPMIN, FIXAFTER, AMBER12T, LJADD3T, NADDTARGET, & 
                LJADDREP, LJADDATT

USE COMMONS, ONLY : DEBUG

USE GENRIGID, ONLY : ATOMRIGIDCOORDT, DEGFREEDOMS, NRIGIDBODY, NSITEPERBODY, RIGIDGROUPS, MAXSITE, SITESRIGIDBODY, & 
                     RIGIDSINGLES, IINVERSE, AACONVERGENCET

USE INTCOMMONS, ONLY : BONDSFROMFILE

USE INTERNALS_WRAPPER, ONLY : INTWRAP_USEINTERNALS

USE MODAMBER9, ONLY : CHECKCISTRANSALWAYS, CHECKCISTRANSALWAYSDNA, CHECKCISTRANSALWAYSRNA

USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_PTR, C_INT, C_DOUBLE, C_BOOL, C_CHAR

IMPLICIT NONE

INTERFACE
    SUBROUTINE CUDA_LBFGS(N, C_XCOORDS, C_GMAX, C_MFLAG, C_ENERGY, ITMAX, C_ITDONE, C_MAXBFGS, C_MAXERISE, C_RMS, C_CUDAPOT, & 
                         C_DEBUG, C_CUDATIMET, NCALLS, C_COLDFUSION, C_COLDFUSIONLIMIT, C_DGUESS, MUPDATE, C_ATOMRIGIDCOORDT, & 
                         C_DEGFREEDOMS, C_NRIGIDBODY, C_NSITEPERBODY, C_RIGIDGROUPS, C_MAXSITE, C_SITESRIGIDBODY, C_RIGIDSINGLES, & 
                         C_BQMAX, C_IINVERSE, PROJECT, C_FREEZE, C_FROZEN, C_NFREEZE, POTENTIALTIME, C_AACONVERGENCET, & 
                         C_NADDTARGET, C_LJADDREP, C_LJADDATT) BIND(C,NAME="setup_lbfgs")

        IMPORT :: C_INT, C_DOUBLE, C_BOOL, C_CHAR

        INTEGER(KIND=C_INT), INTENT(IN) :: N, & ! 3*no. of atoms
                                           ITMAX, & ! Max. no. of steps allowed in minmization
                                           C_DEGFREEDOMS, & ! Rigid Body Framework (RBF): no. of degrees of freedom
                                           C_NRIGIDBODY, & ! RBF: no. of rigid bodies
                                           C_MAXSITE, & ! RBF: max. no. of sites in a rigid body
                                           MUPDATE, & ! History size
                                           C_NFREEZE, & ! No. of frozen atoms
                                           C_NADDTARGET ! Target cluster size (addressability)

        INTEGER(KIND=C_INT), DIMENSION(C_NRIGIDBODY), INTENT(IN) :: C_NSITEPERBODY ! RBF: no. of rigid body sites

        INTEGER(KIND=C_INT), DIMENSION(C_MAXSITE*C_NRIGIDBODY), INTENT(IN) :: C_RIGIDGROUPS ! RBF: list of atoms in rigid bodies


        INTEGER(KIND=C_INT), DIMENSION(C_DEGFREEDOMS/3 - 2*C_NRIGIDBODY), INTENT(IN) :: C_RIGIDSINGLES ! RBF: list of atoms not in rigid bodies

        INTEGER(KIND=C_INT), INTENT(OUT) :: C_ITDONE, & ! No. of LBFGS iterations done
                                            NCALLS ! No. of potential calls made during this call to LBFGS

        REAL(KIND=C_DOUBLE), INTENT(IN) :: C_GMAX, & ! Convergence tolerance for RMS force 
                                           C_MAXBFGS, & ! Max. step size allowed 
                                           C_MAXERISE, & ! Max. energy rise allowed 
                                           C_COLDFUSIONLIMIT, & ! Limit below which cold fusion is diagnosed and minimization terminated
                                           C_DGUESS, & ! Initial guess for inverse Hessian diagonal elements
                                           C_BQMAX ! (1/5) of threshold under which aaconvergence is used, equivalent to C_GMAX

        REAL(KIND=C_DOUBLE), DIMENSION(C_MAXSITE*3*C_NRIGIDBODY), INTENT(IN) :: C_SITESRIGIDBODY ! RBF: coordinates of the rigid body sites

        REAL(KIND=C_DOUBLE), DIMENSION(C_NRIGIDBODY*3*3), INTENT(IN) :: C_IINVERSE ! RBF: inverse eigenvalues of the unweighted tensor of gyration

        REAL(KIND=C_DOUBLE), DIMENSION(C_NADDTARGET, C_NADDTARGET), INTENT(IN) ::C_LJADDREP, C_LJADDATT ! Repulsive/attractive epsilon matrix

        REAL(KIND=C_DOUBLE), DIMENSION(N), INTENT(INOUT) :: C_XCOORDS ! Coordinates

        REAL(KIND=C_DOUBLE), INTENT(OUT) :: C_ENERGY, & ! Energy
                                            C_RMS, & ! RMS force
                                            POTENTIALTIME ! Time taken in calculating potential

        LOGICAL(KIND=C_BOOL), INTENT(IN) :: C_DEBUG, & ! If true, print debug info
                                            C_CUDATIMET, & ! If true, print timing info
                                            C_ATOMRIGIDCOORDT, & ! If false, use rigid body coordinates
                                            PROJECT, & ! If true, project out the components of the gradient along the eigenvector
                                            C_FREEZE, & ! If true, freeze some specified atoms
                                            C_AACONVERGENCET ! If true, use more accurate method of calculating rigid RMS force

        LOGICAL(KIND=C_BOOL), DIMENSION(N/3), INTENT(IN) :: C_FROZEN ! Logical array specifying frozen atoms

        LOGICAL(KIND=C_BOOL), INTENT(OUT) :: C_MFLAG, & ! True if quench converged
                                             C_COLDFUSION ! Set to true during minimization if cold fusion diagnosed

        CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN) :: C_CUDAPOT ! Character specifying the CUDA potential to be used

    END SUBROUTINE CUDA_LBFGS
END INTERFACE

INTERFACE
    SUBROUTINE CUDA_ENEGRAD_CPUTOGPU(NATOMS, COORDS, C_TOTENERGY, C_GRADIENTS, C_NADDTARGET, C_LJADDREP, C_LJADDATT, C_CUDAPOT, & 
                                    C_CUDATIMET, POTENTIALTIME) BIND(C,NAME="setup_potential_cputogpu")

        IMPORT :: C_INT, C_DOUBLE, C_BOOL, C_CHAR

        INTEGER(KIND=C_INT), INTENT(IN) :: NATOMS, & ! No. of atoms
                                           C_NADDTARGET ! Target cluster size (addressability)

        REAL(KIND=C_DOUBLE), DIMENSION(3*NATOMS), INTENT(IN) :: COORDS ! Atomic coordinates
        REAL(KIND=C_DOUBLE), DIMENSION(C_NADDTARGET, C_NADDTARGET), INTENT(IN) ::C_LJADDREP, C_LJADDATT ! Repulsive/attractive epsilon matrix
        REAL(KIND=C_DOUBLE), INTENT(OUT) :: C_TOTENERGY, & ! Total energy of the system
                                            POTENTIALTIME ! Time taken in calculating potential
        REAL(KIND=C_DOUBLE), DIMENSION(3*NATOMS), INTENT(OUT) :: C_GRADIENTS ! Gradient of the energy w.r.t. each atomic coordinate

        LOGICAL(KIND=C_BOOL), INTENT(IN) :: C_CUDATIMET ! If true, print timing info

        CHARACTER(LEN=1, KIND=C_CHAR), INTENT(IN) :: C_CUDAPOT ! Character specifying the CUDA potential to be used

    END SUBROUTINE CUDA_ENEGRAD_CPUTOGPU
END INTERFACE

CONTAINS

    SUBROUTINE CUDA_LBFGS_WRAPPER(N, MUPDATE, XCOORDS, MFLAG, ENERGY, RMS, ITMAX, ITDONE, PROJECT)

        ! Variables passed as *arguments through this wrapper* (not common) with intent in for CUDA_LBFGS are converted directly
        INTEGER(KIND=C_INT) :: N, ITMAX, C_ITDONE, NCALLS, C_DEGFREEDOMS, C_NRIGIDBODY, C_MAXSITE, MUPDATE, C_NFREEZE, C_NADDTARGET
        INTEGER(KIND=C_INT), DIMENSION(NRIGIDBODY) :: C_NSITEPERBODY
        INTEGER(KIND=C_INT), DIMENSION(MAXSITE*NRIGIDBODY) :: C_RIGIDGROUPS
        INTEGER(KIND=C_INT), DIMENSION(DEGFREEDOMS/3 - 2*NRIGIDBODY) :: C_RIGIDSINGLES

        REAL(KIND=C_DOUBLE) :: C_MAXBFGS, C_MAXERISE, C_ENERGY, C_RMS, C_COLDFUSIONLIMIT, C_DGUESS, C_BQMAX, C_GMAX, POTENTIALTIME
        REAL(KIND=C_DOUBLE), DIMENSION(MAXSITE*3*NRIGIDBODY) :: C_SITESRIGIDBODY
        REAL(KIND=C_DOUBLE), DIMENSION(NRIGIDBODY*3*3) :: C_IINVERSE
        REAL(KIND=C_DOUBLE), DIMENSION(NADDTARGET*NADDTARGET) :: C_LJADDREP, C_LJADDATT
        REAL(KIND=C_DOUBLE), DIMENSION(N) :: C_XCOORDS

        LOGICAL(KIND=C_BOOL) :: C_MFLAG, C_DEBUG, C_CUDATIMET, C_ATOMRIGIDCOORDT, C_COLDFUSION, PROJECT, C_FREEZE, C_AACONVERGENCET
        LOGICAL(KIND=C_BOOL), DIMENSION(N/3) :: C_FROZEN

        CHARACTER(LEN=1, KIND=C_CHAR) :: C_CUDAPOT

        ! Variables passed as *arguments through this wrapper* (not common) with intent out for CUDA_LBFGS are not passed into it
        ! Therefore uninitialised C types are passed in and converted types are copied back after the call

        INTEGER :: I, J, K, ITDONE, FCALL, NPCALL, ECALL, SCALL
        DOUBLE PRECISION :: ENERGY, RMS, FTIME, ETIME, STIME
        DOUBLE PRECISION, DIMENSION(N) :: XCOORDS
        LOGICAL :: MFLAG, KNOWE, KNOWG, KNOWH, DRAGT, DOINTERNALSTRANSFORM
        COMMON /PCALL/ NPCALL, ECALL, FCALL, SCALL, ETIME, FTIME, STIME
        COMMON /KNOWN/ KNOWE, KNOWG, KNOWH
        COMMON /RUNTYPE/ DRAGT

        IF (CHECKCISTRANSALWAYSDNA .OR. CHECKCISTRANSALWAYSRNA) THEN
            WRITE(*,'(A)') " modcudalbfgs> Warning: Subroutine CHIRALITY_CHECK not implemented for LBFGS on GPU. "
        END IF

        IF (CHECKCISTRANSALWAYS) THEN
            WRITE(*,'(A)') " modcudalbfgs> Warning: Subroutine CIS_TRANS_CHECK not implemented for LBFGS on GPU. "
        END IF

        DOINTERNALSTRANSFORM = INTWRAP_USEINTERNALS()

        IF (PRINTPTS) THEN
            WRITE(*,'(A)') " modcudalbfgs> Keyword NOPOINTS must be used. &
                                 Functionality to print coordinates/energies will be added at some point in the future. "
            STOP
        END IF

        IF (REVERSEUPHILLT .OR. GRADSQ .OR. DRAGT .OR. TWOENDS .OR. DOINTERNALSTRANSFORM .OR. BONDSFROMFILE .OR. PV .OR. & 
            DUMPMAG .OR. (NSTEPMIN .GT. 0) .OR. (FIXAFTER .GT. 0)) THEN
            WRITE(*,'(A)') " modcudalbfgs> Keywords REVERSEUPHILLT, GRADSQ, DRAG, TWOENDS, INTMIN, BONDSFROMFILE, PV, & 
                                 DUMPMAG, STEPMIN, FIXAFTER are not yet supported. & 
                                 Contact rgm38 if you would like a feature to be added. "
            WRITE(*,'(A)') " modcudalbfgs> Disclaimer: this list might not be exhaustive! "
            STOP
        END IF

        ! Variables from common blocks or modules with intent in or inout are copied into C types

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

        DO I = 1,N
            C_XCOORDS(I) = XCOORDS(I)
        END DO

        DO I = 1,(N/3)
            C_FROZEN(I) = FROZEN(I)
        END DO

        C_CUDAPOT = CUDAPOT
        C_DEBUG = DEBUG
        C_CUDATIMET = CUDATIMET
        C_MAXBFGS = MAXBFGS
        C_MAXERISE = MAXERISE
        C_ATOMRIGIDCOORDT = ATOMRIGIDCOORDT
        C_DEGFREEDOMS = DEGFREEDOMS
        C_NRIGIDBODY = NRIGIDBODY
        C_MAXSITE = MAXSITE
        C_COLDFUSIONLIMIT = COLDFUSIONLIMIT
        C_DGUESS = DGUESS
        C_BQMAX = GMAX
        C_GMAX = GMAX
        C_FREEZE = FREEZE
        C_NFREEZE = NFREEZE
        C_AACONVERGENCET = AACONVERGENCET
        C_NADDTARGET = NADDTARGET

        ! 'C_' prefix denotes those variables which have intent out or inout or are copies of those from common blocks/modules 
        CALL CUDA_LBFGS(N, C_XCOORDS, C_GMAX, C_MFLAG, C_ENERGY, ITMAX, C_ITDONE, C_MAXBFGS, C_MAXERISE, C_RMS, C_CUDAPOT, & 
                       C_DEBUG, C_CUDATIMET, NCALLS, C_COLDFUSION, C_COLDFUSIONLIMIT, C_DGUESS, MUPDATE, & 
                       C_ATOMRIGIDCOORDT, C_DEGFREEDOMS, C_NRIGIDBODY, C_NSITEPERBODY, C_RIGIDGROUPS, C_MAXSITE, & 
                       C_SITESRIGIDBODY, C_RIGIDSINGLES, C_BQMAX, C_IINVERSE, PROJECT, C_FREEZE, C_FROZEN, C_NFREEZE, & 
                       POTENTIALTIME, C_AACONVERGENCET, C_NADDTARGET, C_LJADDREP, C_LJADDATT)

        ! Make sure C types with intent out or inout are coverted back to Fortran ones

        DO I = 1,N
            XCOORDS(I) = DBLE(C_XCOORDS(I))
        END DO

        ENERGY = DBLE(C_ENERGY)
        RMS = DBLE(C_RMS)
        MFLAG = LOGICAL(C_MFLAG)
        ITDONE = INT(C_ITDONE)
        CFUSIONT = LOGICAL(C_COLDFUSION)

        IF (CFUSIONT) THEN
            WRITE(*,'(A,2G20.10)') " modcudalbfgs> Cold fusion diagnosed - step discarded, energy, limit=", ENERGY, COLDFUSIONLIMIT
            ENERGY = 1.0D60
            RMS = 1.0D1
        END IF

        FCALL = FCALL + INT(NCALLS)
        FTIME = FTIME + DBLE(POTENTIALTIME)

        KNOWE=.FALSE.
        KNOWG=.FALSE.
        KNOWH=.FALSE.

    END SUBROUTINE CUDA_LBFGS_WRAPPER


    SUBROUTINE CUDA_ENEGRAD_WRAPPER(NATOMS, COORDS, TOTENERGY, GRADIENTS)

        INTEGER(KIND=C_INT) :: NATOMS, C_NADDTARGET

        REAL(KIND=C_DOUBLE) :: C_TOTENERGY, POTENTIALTIME
        REAL(KIND=C_DOUBLE), DIMENSION(3*NATOMS) :: COORDS, C_GRADIENTS
        REAL(KIND=C_DOUBLE), DIMENSION(NADDTARGET*NADDTARGET) :: C_LJADDREP, C_LJADDATT

        LOGICAL(KIND=C_BOOL) :: C_CUDATIMET

        CHARACTER(LEN=1, KIND=C_CHAR) :: C_CUDAPOT

        INTEGER :: X, I, J

        DOUBLE PRECISION :: TOTENERGY
        DOUBLE PRECISION, DIMENSION(3*NATOMS) :: GRADIENTS

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

        C_NADDTARGET = NADDTARGET
        C_CUDAPOT = CUDAPOT
        C_CUDATIMET = CUDATIMET

        ! Calculates the energy and gradients on the GPU using the GB potential
        CALL CUDA_ENEGRAD_CPUTOGPU(3*NATOMS, COORDS, C_TOTENERGY, C_GRADIENTS, C_NADDTARGET, C_LJADDREP, C_LJADDATT, C_CUDAPOT, & 
                                  C_CUDATIMET, POTENTIALTIME)

        TOTENERGY = DBLE(C_TOTENERGY)

        DO X = 1,(3*NATOMS)
            GRADIENTS(X) = DBLE(C_GRADIENTS(X))
        END DO

        ! FTIME and FCALL are updated in potential.f

    END SUBROUTINE CUDA_ENEGRAD_WRAPPER

END MODULE MODCUDALBFGS
