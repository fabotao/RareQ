#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]


// [[Rcpp::export]]
Rcpp::NumericVector knn_vote(const arma::rowvec& V, const arma::mat& M, const int steps) {

  // ARGUMENTS

  int c = V.size();
  arma::rowvec out(c);
  arma::rowvec out_temp(c);
  out = V;

  for (int j = 0; j < steps; j++) {
    out_temp = out;
    for (int i = 0; i < c; i++) {
      arma::rowvec x = M.row(i);
      arma::vec y = out_temp(arma::conv_to<arma::uvec>::from(x - 1));
      arma::vec y_uniq = arma::unique(y);
      arma::vec y_cnt(y_uniq.size());
      for (int k = 0; k < y_uniq.n_elem; ++k) {
        arma::uvec idx = find(y == y_uniq(k));
        y_cnt(k) = idx.n_elem;
      }
      if(y_cnt.max()==1){
        out(i) = y(0);
      }else{
        out(i) = y_uniq(y_cnt.index_max());
      }
    }
  }
  return (Rcpp::as<Rcpp::NumericVector>(Rcpp::wrap(out)));
}





