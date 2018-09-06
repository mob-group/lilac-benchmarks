void* native_spmv(double* ov, double* a, double* iv, int* rowstr, int* colidx, int* rows)
{
  int i, j;
  for(i = 0; i < *rows; i++)
  {
      double value = 0.0;
      for(j = rowstr[i] - 1; j < rowstr[i+1] - 1; j++)
          value += a[j] * iv[colidx[j] - 1];
      ov[i] = value;
  }
  return 0;
}

void* f_native_spmv(float* ov, float* a, float* iv, int* rowstr, int* colidx, int* rows)
{
  int i, j;
  for(i = 0; i < *rows; i++)
  {
      float value = 0.0;
      for(j = rowstr[i] - 1; j < rowstr[i+1] - 1; j++)
          value += a[j] * iv[colidx[j] - 1];
      ov[i] = value;
  }
  return 0;
}
