#include <mkl.h>
#include <unistd.h>
#include <signal.h>
#include <sys/mman.h>

struct sigaction old_sigaction;
static void* data_begin    = 0;
static void* data_end      = 0;
static void* aligned_begin = 0;
static void* aligned_end   = 0;

static void handler(int sig, siginfo_t *si, void *unused)
{
    if(si->si_addr >= aligned_begin && si->si_addr < aligned_end)
    {
        mprotect(aligned_begin, aligned_end-aligned_begin, PROT_READ | PROT_WRITE | PROT_EXEC);
        data_begin = 0;
    }
    else if(old_sigaction.sa_sigaction)
    {
        old_sigaction.sa_sigaction(sig, si, unused);
    }
}

void* spmv_harness_(double* ov, double* a, double* iv, int* rowstr, int* colidx, int* rows)
{
    static int cols = 0;

    if(data_begin != 0 && data_begin != colidx)
    {
        mprotect(aligned_begin, aligned_end-aligned_begin, PROT_READ | PROT_WRITE | PROT_EXEC);
        data_begin = 0;
    }
    if(data_begin == 0)
    {
        cols = 0;
        for(int i = rowstr[0] - 1; i < rowstr[*rows] - 1; i++)
            if(colidx[i] >= cols) cols = colidx[i];
        data_begin = colidx;
        data_end   = colidx + rowstr[*rows];

        struct sigaction sa;
        sa.sa_flags = SA_SIGINFO;
        sigemptyset(&sa.sa_mask);
        sa.sa_sigaction = handler;

        sigaction(SIGSEGV, &sa, &old_sigaction);

        aligned_begin = data_begin;
        aligned_end   = data_end + sysconf(_SC_PAGE_SIZE) - 1;
        aligned_begin -= (size_t)aligned_begin % sysconf(_SC_PAGE_SIZE);
        aligned_end   -= (size_t)aligned_end   % sysconf(_SC_PAGE_SIZE);

        mprotect(aligned_begin, aligned_end-aligned_begin, PROT_READ | PROT_EXEC);
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
