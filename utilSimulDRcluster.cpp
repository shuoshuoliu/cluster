
// #define ARMA_64BIT_WORD 1
#include "RcppArmadillo.h"
// [[Rcpp::depends(RcppArmadillo)]]

#include <Rcpp.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
// #include <time.h>
// #include <iostream>

// #define INT_MIN (-INT_MAX - 1)

using namespace Rcpp;
using namespace arma;
using namespace std;


/*
 * Auxiliary
 */
//' @keywords internal
//' @noRd
//' 
//' Data X has been centered.

//' Evaluate the log determinant of the covariance matrix and energy function values on each sample
//' point, by SVD decomposition to achieve efficeint computation.
//' @param X: a n x p  matrix, such as gene expression logcount matrix
//' @param Lam_vec0: a column vector, variance of each error component in factor model
//' @param W0: a p x q matrix, loading matrix
//' @param Ck: a q x q matrix,
//' @param Muk: a q x K matrix, gaussian mixture mean component
//' @param logdSk: a real number, the returned log determinant of covariance matrix
//' @param mSk: a colunmn vector, the returned energy funciton values. 
//' 
//' @return No return
//' 


// [[Rcpp::export]]
List multi_det_SkCpp2(const arma::mat& X, const arma::vec& Lam_vec0,
                      const arma::mat& W0, const arma::mat& Ck, 
                      const arma::rowvec Muk, const arma::mat& Sigmak) {
  int n = X.n_rows, q = W0.n_cols;
  
  arma::mat WC12, tmp2;
  arma::vec tmp1, s, tmp3;
  arma::mat U, V, X_tk;
  
  // Compute the SVD of Sigmak and calculate logdSk
  svd(U, s, V, Sigmak);
  WC12 = W0 * (U.each_row() % trans(sqrt(s))); 
  WC12 = WC12.each_col() % (1.0 / sqrt(Lam_vec0)); 
  arma::vec d = svd(WC12);
  double logdSk = -accu(log(1 + d % d)) - accu(log(Lam_vec0));
  
  // Compute the SVD of the inverse of Ck and calculate mSk
  //std::cout << "Ck:\n" << Ck.i() << std::endl;
  svd(U, s, V, Ck.i());
  WC12 = W0 * (U.each_row() % trans(sqrt(s))); 
  WC12 = WC12.each_col() % (1.0 / sqrt(Lam_vec0)); 
  X_tk = (X - repmat(Muk * W0.t(), n, 1)) % trans(repmat(1.0 / sqrt(Lam_vec0), 1, n)); 
  tmp1 = sum(X_tk % X_tk, 1); 
  tmp2 = X_tk * WC12; 
  tmp3 = sum(tmp2 % tmp2, 1); 
  arma::vec mSk = tmp1 - tmp3;
  
  // Convert Armadillo vectors to Rcpp vectors
  NumericVector mSkR = wrap(mSk);
  
  // Create a list to store the results
  List result;
  result["logdSk"] = logdSk;
  result["mSk"] = mSkR;
  
  return result;
}


//' Evaluate the diagonal elements of three matrix multiplication such as W0*Cki*W0^T
//'  by SVD decomposition to achieve efficeint computation.
//' @param Cki: a q x q matrix,
//' @param W0: a p x q matrix
//' 
//' @return a column vector
//'   

vec decomp(const mat& Cki, const mat& W0){
  vec s, tmp1;
  mat U, V, WC12;
  // std::cout << "ak:\n" << Cki << std::endl;
  svd(U, s, V, Cki);
  WC12 = W0 * (U.each_row() % trans(sqrt(s)));
  tmp1 = sum(WC12 % WC12, 1);
  return tmp1;
}  


// [[Rcpp::export]] 
mat calculateX(const mat& R, const cube& Ez){
  int k, K= R.n_cols, q= Ez.n_cols, n= R.n_rows;
  mat Ezz(n, q, fill::zeros); // estimate Z, factor matrix
  for(k=0; k<K; k++){
    Ezz +=  Ez.slice(k) % repmat(R.col(k), 1, q);
  }
  return Ezz;
}



//' Evaluate the diagonal elements of three matrix multiplication such as W0*Cki*W0^T
//'  by SVD decomposition to achieve efficeint computation.
//' @param Cki: a q x q matrix,
//' @param W0: a p x q matrix
//' 
//' @return a column vector
//'  
// [[Rcpp::export]] 
mat update_W0(const mat& X, const mat& R, const cube& Ez, const cube& Ci,  const vec& N){
  int k, K= R.n_cols, q= Ez.n_cols, n= R.n_rows;
  mat tmpMat(n,q, fill::zeros), Ezzt(q,q, fill::zeros), tmpMat2;
  // vec N = sum(R.t(), 1);
  for(k=0; k<K; ++k){
    tmpMat2= repmat(R.col(k), 1, q) % Ez.slice(k);
    tmpMat += tmpMat2;
    Ezzt+= tmpMat2.t() * Ez.slice(k) + N(k) * Ci.slice(k);
  }
  return X.t() * tmpMat * Ezzt.i();
}

// update Sigma0
// [[Rcpp::export]]
cube update_Sigma0(const mat& R, const cube& Ez, const cube& Ci, const mat& Mu,  
                   const vec&  N){
  int k, K= R.n_cols, q= Mu.n_cols, n= R.n_rows;
  cube Sigma0(q,q,K);
  for(k = 0; k<K; ++k){
    Sigma0.slice(k) = (trans(Ez.slice(k) - repmat(Mu.row(k), n, 1)) % trans(repmat(R.col(k), 1, q))) * (Ez.slice(k) - repmat(Mu.row(k), n, 1)) + N(k)*  Ci.slice(k);
    Sigma0.slice(k) = Sigma0.slice(k) / N(k);
    
  }
  return(Sigma0);
}

// update Lambda
// [[Rcpp::export]]
vec update_Lam(const mat& R, const mat& X, const mat& W, const cube& Ez, const cube& Ci){
  int k, K= R.n_cols, p = X.n_cols;
  rowvec N = sum(R);
  vec Lsum(p,fill::zeros), Lam;
  mat tmpXk;
  for(k=0; k<K; ++k){
    tmpXk = (X - Ez.slice(k) * W.t() );
    Lsum += trans(sum(tmpXk % tmpXk % repmat(R.col(k), 1, p)));
    Lsum += N(k) * decomp(Ci.slice(k), W);
  }
  Lam = Lsum/(X.n_rows*1.0);
  return Lam; 
}

// // update sigma20: homo variance
// vec update_Lam2(const mat& R, const mat& X, const mat& W, const cube& Ez, const cube& Ci){
//   int k, K= R.n_cols, p = X.n_cols;
//   mat tmpXk;
//   double term1=0, term2=0;
//   for(k=0; k <K; ++k){
//     tmpXk = (X - Ez.slice(k) * W.t() );
//     term1 += accu(sum(tmpXk % tmpXk, 1) % R.col(k));
//     term2  += accu(decomp(Ci.slice(k), W))* accu(R.col(k));
//   }
//   double sigma20 = (term1+ term2)/ X.n_elem;
//   return(ones(p)* sigma20);
// }

//Calculate Q function
double Q_fun(const mat& X, const mat& R,  const cube& Ez, const arma::cube& Ci, 
             const mat& W0, const mat& Mu0, const cube& Sigma0, const vec& Pi0, 
             const vec& Lam_vec0){
  
  double Q = 0, tmp_scalar =0;
  int  K = Pi0.n_elem, n = X.n_rows;
  int q = Mu0.n_cols;
  mat tmpSig(q, q, fill::zeros), tmpMat;
  colvec Ezik;
  mat Ezzt(q,q, fill::zeros);
  for(int k=0; k<K; k++){
    tmpSig = Sigma0.slice(k);
    tmpMat = Ez.slice(k);
    for(int i = 0; i<n; i++){
      Ezik =  trans(tmpMat.row(i)); 
      Ezzt = Ci.slice(k) + Ezik * Ezik.t();
      
      tmp_scalar =  0.5* accu(log(Lam_vec0)) + 0.5* accu( (X.row(i) % X.row(i)) / Lam_vec0.t()) +
        0.5 * arma::as_scalar(accu(W0.t()%trans(repmat(1.0/Lam_vec0, 1, q))*W0 % Ezzt.t())- 2* (X.row(i)% (1.0/Lam_vec0))*W0*Ezik);
      Q +=  - R(i,k) * tmp_scalar + R(i,k)*(log(Pi0(k)) -  0.5* log(det(tmpSig))- 0.5* accu(inv_sympd(tmpSig)%
        Ezzt.t())+ arma::as_scalar(Mu0.row(k) * solve(tmpSig, (Ezik- 0.5*trans(Mu0.row(k)))))); 
    }
  }
  return Q;
}
