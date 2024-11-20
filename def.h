void init_mt64(int64_t seed = -1);
double L2_norm(double  *cc, const int n);
double gen_po2(double *AA_d, int n, double factor, bool normal);
double gen_po(double *AA_d, int n, double factor, bool normal);
double gen_po(double *AA_d, int n, double factor);
double gen_po2(double *AA_d, int n, double factor);

double gen_ge(double *AA_d, int n, double factor, bool normal);
double gen_ge(double *AA_d, int n, double factor);

void gen_rand(double *AA_d, int n, int ld, double factor);
void gen_rand(double *AA_d, int n, int ld, double factor, bool normal);

extern "C" {
  void spotrf_ (char * UPLO, int* N, float*  A, int* lda, int* INFO );
  void spotrs_ (char * UPLO, int* N, int* NRHS, float* A, int* lda, float* B, int* ldb, int* INFO );
  void sgemv_  (char * TRANS, int* M, int* N, float* alpha, float* A, int* lda,
		float* X, int* INCX, float* beta, float* B, int* INCY);
  void sgetrf_ (int* M, int* N, float* A, int* lda, int* IPIV, int* INFO );
  void sgetrs_ (char * TRANS, int* N, int* NRHS, float * A, int* lda, int * IPIV, float * B, int* ldb, int* INFO );
  
  void dpotrf_ (char * UPLO, int* N, double* A, int* lda, int* INFO );
  void dpotrs_ (char * UPLO, int* N, int* NRHS, double* A, int* lda, double* B, int* ldb, int* INFO );
  void dgemv_  (char * TRANS, int* M, int* N, double* alpha, double* A, int* lda,
		double* X, int* INCX, double* beta, double* B, int* INCY);
  void dgetrf_ (int* M, int* N, double* A, int* lda, int* IPIV, int* INFO );
  void dgetrs_ (char * TRANS, int* N, int* NRHS, double* A, int* lda, int * IPIV, double* B, int* ldb, int* INFO );

  void dgemm_(char * TRANSA, char * TRANSB, int* M, int* N, int *K, double *alpha, double* A, int* lda,
	      double* b, int* ldb, double *beta, double *c, int *ldc);
  void dgecon_(char * NORM, int *N, double * A, int *lda, double *ANORM, double *RCOND, double *work, int *iwork, int *INFO);
  double dlange_(char * NORM, int *M, int *N, double * A, int *lda, double *work);

  void dpocon_(char * UPLO, int *N, double * A, int *lda, double *ANORM, double *RCOND, double *work, int *iwork, int *INFO);
};
 
