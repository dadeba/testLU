#include <iostream>
#include <random>
#include <chrono>
#include <cstring>
#include <malloc.h>
#include <unistd.h>
#include <mpblas_double.h>
#include <mplapack_double.h>

#include "def.h"

int main(int narg, char **argv)
{
  double factor = 1.0;

  using Clock = std::chrono::high_resolution_clock;
  using std::chrono::nanoseconds;
  using std::chrono::duration_cast;
  
  mplapackint n;
  uint64_t seed = 20241013;
  init_mt64(seed);
  
  if (narg == 2) factor = atof(argv[1]);

  //mplapackint n_test[] = {8, 100, 200, 300, 400, 500, 700, 900, 1000, 1200, 1400, 1700, 2000, 2200};
  mplapackint n_test[] = {100, 900, 1399, 2200};
  
  std::cout << std::scientific;
  
  for(auto tt = 0; tt < sizeof(n_test)/sizeof(mplapackint); tt++) {
    n = n_test[tt];  
    size_t nn = (size_t)n*(size_t)n;
  
    double *aa_d, *AA_d, *xx_d;
    aa_d = (double *)memalign(64, sizeof(double)*nn*2);
    AA_d = (double *)memalign(64, sizeof(double)*nn*2);
    xx_d = (double *)memalign(64, sizeof(double)*n);

    double rcond_a = gen_ge(aa_d, n, factor, true);
    memcpy(AA_d, aa_d, sizeof(double)*nn);
    
    mplapackint ipiv_d[10000];
    mplapackint err;

    auto st = Clock::now();
    auto en = Clock::now();
    
    st = Clock::now();
    Rgetrf(n, n, aa_d, n, ipiv_d, err);
    en = Clock::now();
    auto time_b64 = duration_cast<nanoseconds>(en-st).count()/1.0e9;

    for(mplapackint i = 0; i < n; i++) xx_d[i] = 1.0/sqrt((double)n);

    mplapackint nrhs = 1;
    mplapackint ldb = n;
    mplapackint inc = 1;
    double err_d = 0.0;
    double err_hpl = 0.0;
    {
      double *b  = (double *)memalign(64, sizeof(double)*n);
      double *b_ = (double *)memalign(64, sizeof(double)*n);
      double *x  = (double *)memalign(64, sizeof(double)*n);
      for(mplapackint i = 0; i < n; i++) b[i] = 0.0;
      for(mplapackint i = 0; i < n; i++) x[i] = xx_d[i];

      double alpha = 1.0;
      double beta  = 0.0;
      Rgemv((char *)"N", n, n, alpha, AA_d, n, x, inc, beta, b, inc);
      memcpy(x, b, sizeof(double)*n);
      Rgetrs((char *)"N", n, nrhs, aa_d, n, ipiv_d, x, ldb, err);
      Rgemv((char *)"N",  n, n, alpha, AA_d, n, x, inc, beta, b_, inc);
      for(mplapackint i = 0; i < n; i++) {
	b_[i] = b[i] - b_[i];
      }
      err_d = L2_norm(b_,n)/L2_norm(b,n);

      err_hpl = Rlange((char*)"M", n, inc, b_, n, NULL);
      double mmax = Rlange((char*)"M", n, n, AA_d, n, NULL);
      double xmax = Rlange((char*)"M", n, inc, x, n, NULL);
      double bmax = Rlange((char*)"M", n, inc, b, n, NULL);

      static const double eps  = 0x1p-52; 
      err_hpl = err_hpl/(eps*(mmax*xmax + bmax)*(double)n);
      
      free(b);
      free(b_);
      free(x);
    }
    
    std::cout << n << "\t";
    std::cout << err_d << "\t";
    std::cout << err_hpl << "\t";
    std::cout << 1.0/rcond_a << "\t";
    std::cout << time_b64 << "\n";
    std::cout << std::flush;
    
    free(aa_d);
    free(AA_d);
    free(xx_d);
  }

  return 1;
}
