// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::export]]
double lse(vec log_x) {
  double m = max(log_x);
  return m + log(sum(exp(log_x - m)));
}
