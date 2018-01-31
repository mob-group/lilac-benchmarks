MODULE MODCUDABFGSTS ! Compiled with all versions of OPTIM except CUDAOPTIM

CONTAINS

    SUBROUTINE CUDA_BFGSTS_WRAPPER(NSTEPS,Q,ENERGY,MFLAG,RMS,EVTS,VECS,ITDONE)

        INTEGER :: NSTEPS, ITDONE
        DOUBLE PRECISION, DIMENSION(1) :: Q, VECS
        DOUBLE PRECISION :: ENERGY, RMS, EVTS
        LOGICAL :: MFLAG

        RETURN

    END SUBROUTINE CUDA_BFGSTS_WRAPPER

END MODULE MODCUDABFGSTS
