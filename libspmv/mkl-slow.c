#include "mkl-model.h"
#include "native-impl.h"

#include <mkl.h>
#include <unistd.h>
#include <signal.h>
#include <sys/mman.h>

void* spmv_harness_(double* ov, double* a, double* iv, int* rowstr, int* colidx, int* rows)
{
	int cols = 0;

        for(int i = rowstr[0] - 1; i < rowstr[*rows] - 1; i++) {
        	if(colidx[i] >= cols) {
			cols = colidx[i];
		}
	}

        sparse_matrix_t A;
        mkl_sparse_d_create_csr(&A, SPARSE_INDEX_BASE_ONE, *rows, cols-1, rowstr, rowstr+1, colidx, a);

        struct matrix_descr dscr;
        dscr.type = SPARSE_MATRIX_TYPE_GENERAL;
        dscr.mode = SPARSE_FILL_MODE_LOWER;
        dscr.diag = SPARSE_DIAG_NON_UNIT;

        mkl_sparse_d_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, A, dscr, iv, 0.0, ov);
        return 0;
}

void* f_spmv_harness_(float* ov, float* a, float* iv, int* rowstr, int* colidx, int* rows)
{
	int cols = 0;

        for(int i = rowstr[0] - 1; i < rowstr[*rows] - 1; i++) {
        	if(colidx[i] >= cols) {
			cols = colidx[i];
		}
	}

        sparse_matrix_t A;
        mkl_sparse_s_create_csr(&A, SPARSE_INDEX_BASE_ONE, *rows, cols-1, rowstr, rowstr+1, colidx, a);

        struct matrix_descr dscr;
        dscr.type = SPARSE_MATRIX_TYPE_GENERAL;
        dscr.mode = SPARSE_FILL_MODE_LOWER;
        dscr.diag = SPARSE_DIAG_NON_UNIT;

        mkl_sparse_s_mv(SPARSE_OPERATION_NON_TRANSPOSE, 1.0, A, dscr, iv, 0.0, ov);
        return 0;
}
