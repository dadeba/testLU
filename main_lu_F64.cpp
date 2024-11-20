#include <iostream>
#include <random>
#include <chrono>
#include <cstring>
#include <malloc.h>
#include <unistd.h>
#include <mpblas__Float64x.h>
#include <mplapack__Float64x.h>

#include "def.h"
_Float64x L2_norm(_Float64x *cc, const int n);

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

  //  mplapackint n_test[] = {8, 100, 200, 300, 400, 500, 700, 900, 1000, 1200, 1400, 1700, 2000, 2200};
  mplapackint n_test[] = {100, 900, 1399, 2200};
  
  std::cout << std::scientific;
  
  for(long unsigned int tt = 0; tt < sizeof(n_test)/sizeof(mplapackint); tt++) {
    n = n_test[tt];  
    size_t nn = (size_t)n*(size_t)n;

    _Float64x *aa_d = new _Float64x[nn*2];
    _Float64x *AA_d = new _Float64x[nn*2];
    _Float64x *xx_d = new _Float64x[n];

    double *a_init = new double[nn*2];

    double rcond_a = gen_ge(a_init, n, factor, true);
    for(mplapackint i = 0; i < (mplapackint)nn; i++) {
      aa_d[i] = a_init[i];
    }
    memcpy(AA_d, aa_d, sizeof(_Float64x)*nn);
    
    mplapackint ipiv_d[10000];
    mplapackint err;

    auto st = Clock::now();
    auto en = Clock::now();
    
    st = Clock::now();
    Rgetrf(n, n, aa_d, n, ipiv_d, err);
    en = Clock::now();
    auto time_b64 = duration_cast<nanoseconds>(en-st).count()/1.0e9;

    for(int i = 0; i < n; i++) xx_d[i] = 1.0/sqrt((double)n);

    mplapackint nrhs = 1;
    mplapackint ldb = n;
    mplapackint inc = 1;
    double err_d = 0.0;
    _Float64x err_hpl = 0.0;
    {
      _Float64x *b  = new _Float64x[n];
      _Float64x *b_ = new _Float64x[n];
      _Float64x *x  = new _Float64x[n];
      for(int i = 0; i < n; i++) b[i] = 0.0;
      for(int i = 0; i < n; i++) x[i] = xx_d[i];

      _Float64x alpha = 1.0;
      _Float64x beta  = 0.0;
      Rgemv((char *)"N", n, n, alpha, AA_d, n, x, inc, beta, b, inc);
      memcpy(x, b, sizeof(_Float64x)*n);
      Rgetrs((char *)"N", n, nrhs, aa_d, n, ipiv_d, x, ldb, err);
      Rgemv((char *)"N",  n, n, alpha, AA_d, n, x, inc, beta, b_, inc);
      for(int i = 0; i < n; i++) {
	b_[i] = b[i] - b_[i];
      }
      err_d = L2_norm(b_,n)/L2_norm(b,n);

      err_hpl = Rlange((char*)"M", n, inc, b_, n, NULL);
      _Float64x mmax = Rlange((char*)"M", n, n, AA_d, n, NULL);
      _Float64x xmax = Rlange((char*)"M", n, inc, x, n, NULL);
      _Float64x bmax = Rlange((char*)"M", n, inc, b, n, NULL);

      static const _Float64x eps = 0x1p-63; 
      err_hpl = err_hpl/(eps*(mmax*xmax + bmax)*(_Float64x)n);
      
      delete[] b;
      delete[] b_;
      delete[] x;
    }
    
    std::cout << n << "\t";
    std::cout << err_d << "\t";
    std::cout << err_hpl << "\t";
    std::cout << 1.0/rcond_a << "\t";
    std::cout << time_b64 << "\n";
    std::cout << std::flush;
    
    delete[] aa_d;
    delete[] AA_d;
    delete[] xx_d;

    delete[] a_init;
  }

  return 1;
}

_Float64x L2_norm(_Float64x *cc, const int n)
{
  double sum;
  
  sum = 0.0;
  for(int j = 0; j < n; j++) {
    sum += cc[j]*cc[j];
  }

  return sqrt(sum);
}
