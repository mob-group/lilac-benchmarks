void spmv_harness(double* ov, const double* a, const double* iv,
                  const int* rowstr, const int* colidx, int rows) {
    int i, j;
    for(i = 0; i < rows; i++) {
      double value = 0.0;
      for(j = rowstr[i]; j < rowstr[i+1]; j++)
        value += a[j] * iv[colidx[j]];
      ov[i] = value;
    }
}
