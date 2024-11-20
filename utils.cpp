#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <iostream>
#include <random>
#include <cstring>
#include <malloc.h>
#include <algorithm>
#include <cmath>
#include <cfloat>

#include "def.h"

static std::mt19937_64 *mt;

void init_mt64(int64_t seed)
{
  if (mt != NULL) return;
  mt = new std::mt19937_64;
  if (seed == -1) {
    std::random_device rnd;
    seed = rnd();
  }
  mt->seed(seed);
  std::cerr << "seed " << seed << "\n";
}

double L2_norm(double *cc, const int n)
{
  double sum;
  
  sum = 0.0;
  for(int j = 0; j < n; j++) {
    sum += cc[j]*cc[j];
  }

  return sqrt(sum);
}

void gen_rand(double *AA_d, int n, int ld, double factor)
{
  return gen_rand(AA_d, n, ld, factor, false);
}

void gen_rand(double *AA_d, int n, int ld, double factor, bool normal)
{
  init_mt64();
  // create ge
  std::uniform_real_distribution<> distr(-1.0*factor, 1.0*factor);
  std::normal_distribution<> distn(0.0, factor);
  
  for(int j = 0; j < n; j++) {
    for(int i = 0; i < ld; i++) {
      if (normal) AA_d[j*ld+i] = distn(*mt);
      else        AA_d[j*ld+i] = distr(*mt);
    }
  }
}

void Dsymm(double *cc, const int n)
{
  for(int j = 0; j < n; j++) {
    for(int i = 0; i < j; i++) {
      double x = cc[j*n+i];
      cc[i*n + j] = x;
    }
  }
}

double gen_ge(double *AA_d, int n, double factor)
{
  return gen_ge(AA_d, n, factor, false);
}

double gen_po(double *AA_d, int n, double factor)
{
  return gen_po(AA_d, n, factor, false);
}

double gen_po2(double *AA_d, int n, double factor)
{
  return gen_po2(AA_d, n, factor, false);
}
  
double gen_ge(double *AA_d, int n, double factor, bool normal)
{
  init_mt64();
  // create ge
  std::uniform_real_distribution<> distr(-1.0*factor, 1.0*factor);
  std::normal_distribution<> distn(0.0, factor);
  
  double *aa_d;
  aa_d = (double *)memalign(64, sizeof(double)*n*n);
  
  for(int j = 0; j < n; j++) {
    for(int i = 0; i < n; i++) {
      if (normal) AA_d[j*n+i] = distn(*mt);
      else AA_d[j*n+i] = distr(*mt);
    }
  }

  double *work  = (double *)memalign(64, sizeof(double)*n*4);
  int    *iwork = (int *)memalign(64, sizeof(int)*n*4);

  memcpy(aa_d, AA_d, sizeof(double)*n*n);

  double anorm = dlange_((char *)"1", &n, &n, aa_d, &n, work);
  double rcond;

  // LU
  int err;
  dgetrf_(&n, &n, aa_d, &n, iwork, &err);
  dgecon_((char *)"1", &n, aa_d, &n, &anorm, &rcond, work, iwork, &err);

  free(aa_d);
  free(work);
  free(iwork);

  return rcond;
}

double gen_po(double *AA_d, int n, double factor, bool normal)
{
  init_mt64();
  std::uniform_real_distribution<> distr(-1.0*factor, 1.0*factor);
  std::normal_distribution<> distn(0.0, factor);
  
  double *aa_d;
  aa_d = (double *)memalign(64, sizeof(double)*n*n);
  
  for(int j = 0; j < n; j++) {
    for(int i = 0; i < n; i++) {
      if (normal) AA_d[j*n+i] = distn(*mt);
      else AA_d[j*n+i] = distr(*mt);
    }
  }

  int err;
  double alpha_t = 1.0;
  double beta_t = 0.0;
  double *work  = (double *)memalign(64, sizeof(double)*n*4);
  int    *iwork = (int *)memalign(64, sizeof(int)*n*4);

  dgemm_((char *)"T", (char *)"N", &n, &n, &n, &alpha_t, AA_d, &n, AA_d, &n, &beta_t, aa_d, &n);
  
  //  double anorm_2 = dlange_((char *)"2", &n, &n, aa_d, &n, work);
  //  double r_anorm = 1.0/anorm_2;
  //  for(int j = 0; j < n; j++) {
  //    for(int i = 0; i < n; i++) {
  //      aa_d[j*n+i] *= r_anorm;
  //    }
  //  }

  memcpy(AA_d, aa_d, sizeof(double)*n*n);
  double anorm = dlange_((char *)"1", &n, &n, aa_d, &n, work);

  double rcond;

  // Cholesky
  dpotrf_((char *)"U", &n, aa_d, &n, &err);
  dpocon_((char *)"U", &n, aa_d, &n, &anorm, &rcond, work, iwork, &err);

  free(aa_d);
  free(work);
  free(iwork);

  return rcond;
}

double gen_po2(double *AA_d, int n, double factor, bool normal)
{
  init_mt64();
  // create symmetric positive definite matrix
  std::uniform_real_distribution<> distr(-1.0*factor, 1.0*factor);
  std::normal_distribution<> distn(0.0, factor);
  
  double *aa_d;
  aa_d = (double *)memalign(64, sizeof(double)*n*n);
  
  for(int j = 0; j < n; j++) {
    for(int i = 0; i < n; i++) {
      if (normal) AA_d[j*n+i] = fabs(distn(*mt));
      else AA_d[j*n+i] = fabs(distr(*mt));
    }
  }

  Dsymm(AA_d, n);
  int err;

  double *work  = (double *)memalign(64, sizeof(double)*n*4);
  int    *iwork = (int *)memalign(64, sizeof(int)*n*4);

  memcpy(aa_d, AA_d, sizeof(double)*n*n);
  double anorm = dlange_((char *)"1", &n, &n, aa_d, &n, work);

  double rcond;

  // Cholesky
  dpotrf_((char *)"U", &n, aa_d, &n, &err);
  dpocon_((char *)"U", &n, aa_d, &n, &anorm, &rcond, work, iwork, &err);

  free(aa_d);
  free(work);
  free(iwork);

  return rcond;
}
