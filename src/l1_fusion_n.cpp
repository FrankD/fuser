#include <Rcpp.h>
#include <RcppEigen.h>

using namespace Rcpp;

using Eigen::MatrixXd;
using Eigen::VectorXd;


/**
 * Internal function.
 *
 * Calculate the delta.lik matrix.
 *
 * Note: If XX is null, Y MUST be sorted in the order of the groups in the group_sizes vector.
 */
Eigen::MatrixXd calculate_delta_lik(int k, int p, Eigen::MatrixXd B_old, Nullable<List> XX_opt, Nullable<List> XY_opt, Nullable<List> X_list_opt, Eigen::VectorXd Y, NumericVector samp_sizes) {
  // equivalent to a k-length list of p-length vectors
  MatrixXd result(p, k);

  const int num_rows = B_old.rows();

  if (XX_opt.isNotNull()) {
    /*
     * equivalent to
     * lapply(1:k, function(k.i) (XX[[k.i]] %*% B.old[, k.i] - XY[[k.i]])/samp.sizes[k.i])
     */
    List XX = as<List>(XX_opt), XY = as<List>(XY_opt);

    for (int k_i = 0; k_i < k; ++k_i) {
      MatrixXd XX_k_i = XX[k_i];
      MatrixXd XY_k_i = XY[k_i];
      VectorXd B_old_k_i = B_old.block(0, k_i, num_rows, 1);

      result.col(k_i) = (XX_k_i * B_old_k_i - XY_k_i) / samp_sizes[k_i];
    }
  } else {
    /*
     * Replacement for
     * lapply(1:k, function(k.i) {
     *   g = group.names[k.i]
     *   (doubleCrossProd(X.list[[k.i]], B.old[, k.i, drop = FALSE]) - crossprod(X.list[[k.i]], Y[groups == g]))/samp.sizes[k.i]
     * })
     */
    double *next_data = Y.data();

    List X_list = as<List>(X_list_opt);

    for (int k_i = 0; k_i < k; ++k_i) {
      size_t group_size = samp_sizes[k_i];
      Eigen::Map<VectorXd> group_samples(next_data, group_size, 1);
      next_data += group_size; // offset the pointer to point to the next group's data
      MatrixXd X_list_k_i = X_list[k_i];
      result.col(k_i) = (X_list_k_i.adjoint() * (X_list_k_i * B_old.block(0, k_i, num_rows, 1)) - X_list_k_i.adjoint() * group_samples) / group_size;
    }
  }
  return result;
}

/**
 * Internal function.
 *
 * penalty_factors must be a k-length vector or NULL.
 *
 * B.old is a (p × k) OR ((p+1) × k) matrix
 */
static inline Eigen::MatrixXd calculate_B_sparsity(Eigen::MatrixXd B_old, int p, int k, Nullable<NumericVector> penalty_factors) {
  if (penalty_factors.isNotNull()) {
    VectorXd eigen_penalty_factors = as<VectorXd>(penalty_factors.get());
    MatrixXd replicated = eigen_penalty_factors.rowwise().replicate(k); // duplicates the column
    return B_old.block(0, 0, p, k).array() * replicated.array();
  } else {
    return B_old.block(0, 0, p, k);
  }
}

static int number_iterations_taken = -1;

/**
 * Returns the number of iterations taken by the previous invocation of
 * genFusedLassoProximal_loop, or -1 if it has not been invoked.
 */
// [[Rcpp::export]]
NumericVector getNumberNativeIterationsTaken() {
  if (number_iterations_taken < 0) {
    return wrap(NA_REAL);
  } else {
    return wrap(number_iterations_taken);
  }
}

/**
 * Performs the iterations of the Fused Lasso Proximal Gradient Method.
 * This function should not be invoked directly from R.
 * see fusedLassoProximal in l1_fusion.R
 *
 * Y MUST be grouped by group & these groups sorted in order of the group_sizes vector,
 * unless XX is not null.
 */
// [[Rcpp::export]]
Eigen::MatrixXd genFusedLassoProximal_loop(Nullable<List> XX, Nullable<List> XY, Nullable<List> X_list, Eigen::VectorXd Y, NumericVector samp_sizes, Eigen::MatrixXd C, bool intercept, int p, int k, int num_iters, Nullable<NumericVector> penalty_factors, double L_U_inv, Eigen::MatrixXd B_old, double mu, Eigen::MatrixXd W, Eigen::MatrixXd weighted_delta_f, double tol) {
  int i;
  MatrixXd B_sparsity, B_fusion, A_star, B_new;


  for (i = 1; i <= num_iters; ++i) {
    B_sparsity = calculate_B_sparsity(B_old, p, k, penalty_factors);

    B_fusion = B_old;
    if (intercept) {
      // add a row to the B_sparsity matrix if intercept terms are enabled
      B_sparsity.conservativeResize(p + 1, Eigen::NoChange);
      for (size_t i = 0; i < k; ++i) {
        // zero the extra row in B.sparsity and B.fusion
        B_sparsity(p, i) = 0;
        B_fusion(p, i) = 0;
      }
    }

    MatrixXd delta_lik = calculate_delta_lik(k, p, B_old, XX, XY, X_list, Y, samp_sizes);

    MatrixXd A_star_sparsity = B_sparsity * C.block(0, 0, C.rows(), k);
    MatrixXd A_star_fusion = B_fusion * C.block(0, k, C.rows(), C.cols() - k);
    A_star.resize(A_star_sparsity.rows(), A_star_sparsity.cols() + A_star_fusion.cols());
    A_star.block(0, 0, A_star.rows(), k) = A_star_sparsity;
    A_star.block(0, k, A_star.rows(), A_star.cols() - k) = A_star_fusion;
    A_star /= mu;

    {
      double *data = A_star.data();
      for (size_t offset = 0; offset < A_star.size(); ++offset) {
        if (data[offset] > 1) data[offset] = 1;
        else if (data[offset] < -1) data[offset] = -1;
      }
    }

    MatrixXd delta_f = delta_lik + A_star * C.transpose();

    B_new = W - (L_U_inv * delta_f);

    weighted_delta_f = weighted_delta_f + (L_U_inv * 0.5 * i * delta_f);

    W = (i * B_new - 2 * weighted_delta_f) / (i + 2);

    double improvement = (B_old - B_new).cwiseAbs().sum();

    if (improvement < tol * p) break;

    B_old = B_new;
  }
  if (i >= num_iters) {
    Rcpp::warning("Reached max iterations without convergence.");
  }
  number_iterations_taken = i;

  return B_new;
}
