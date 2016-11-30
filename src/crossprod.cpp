// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

// we only include RcppEigen.h which pulls Rcpp.h in for us
#include <RcppEigen.h>

using namespace std;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::Lower;
using Rcpp::as;
// via the depends attribute we tell Rcpp to create hooks for
// RcppEigen so that the build process will know what to do
//
// [[Rcpp::depends(RcppEigen)]]


// Calculate crossproduct XtX
//MatrixXd crossprodCPP(const MatrixXd& X) {
  

//}


// Calculate crossproducts XtX and XtY and return product with 
// Beta matrix
// [[Rcpp::export]]
Eigen::MatrixXd doubleCrossProd(Eigen::MatrixXd A, Eigen::MatrixXd B) {
  //cout << A << " ";
  //const Map<MatrixXi> A(as<Map<MatrixXi> >(AA);
  MatrixXd AtAB = A.adjoint() * (A * B);

  return AtAB;
}




