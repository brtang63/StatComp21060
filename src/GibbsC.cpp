#include <Rcpp.h>
using namespace Rcpp;

//' @title A Gibbs Sampler
//' @description Generate a chain with Gibbs sampler.
//' @param N sample size
//' @return A matrix of random samples.
//'
//' 
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix GibbsC(int N) {
  int n = 100;
  int a = 50;
  int b = 50;
  
  Rcpp::NumericMatrix mat( N, 2);
  Rcpp::NumericVector xt (2);
  
  for (int i = 1; i < N; i++) {
    xt[0] = mat( i-1, 0 );
    xt[1] = mat( i-1, 1 );
    xt[0] = rbinom(1, n, xt[1])[0];
    xt[1] = rbeta(1, xt[0] + a, n - xt[0] + b)[0];
    mat( i, 0 ) = xt[0];
    mat( i, 1 ) = xt[1];
  }
  
  return mat;
}