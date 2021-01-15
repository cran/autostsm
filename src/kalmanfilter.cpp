#include <RcppArmadillo.h>
#include <Rcpp.h>

// Correctly setup the build environment
// [[Rcpp::depends(RcppArmadillo)]]

// Protect against compilers without OpenMP
#ifdef _OPENMP
#include <omp.h>
#endif

// R's implementation of the Moore-Penrose pseudo matrix inverse
// [[Rcpp::export]]
arma::mat Rginv(const arma::mat& m){
  arma::mat U, V;
  arma::vec S;
  arma::svd(U, S, V, m, "dc");
  arma::uvec Positive = arma::find(S > 1E-06 * S(1));
  if(all(Positive)){
    arma::mat D = diagmat(S);
    return V * (1/D * U.t());
  }else if(!any(Positive)){
    return arma::zeros(m.n_rows, m.n_cols);
  }else{
    S.elem(Positive) = 1/S.elem(Positive);
    arma::mat D = diagmat(S);
    return V * D * U.t();
  }
}

// st = Sys.time()
// B_tt = B_tl = matrix(0, ncol = ncol(yt), nrow = nrow(sp$Ft))
// P_tt = P_tl = as.list(NA, ncol(yt))
//
// #Initialize the filter
// B_LL = sp$B0
// P_LL = sp$P0
//
// for(j in 1:ncol(yt)){
//   ################## Kalman filter routine ##################
//   B_tl[, j] = sp$Dt + sp$Ft %*% B_LL  #Initial estimate of unobserved values conditional on t-1
//   P_tl[[j]] = sp$Ft %*% P_LL %*% t(sp$Ft) + sp$Qt #Initial estimate of the covariance matrix conditional on t-1
//   N_tl = yt[, j] - sp$At - sp$Ht %*% B_tl[, j] #Prediction error conditoinal on t-1
//   F_tl = sp$Ht %*% P_tl[[j]] %*% t(sp$Ht) + sp$Rt #Variance of the predictoin error conditional on t-1
//   K_t = P_tl[[j]] %*% t(sp$Ht) %*% MASS::ginv(F_tl) #Kalman gain conditional on t-1
//   B_tt[, j] = B_tl[, j] + K_t %*% N_tl #Final estimate of the unobserved values
//   P_tt[[j]] = P_tl[[j]] - K_t %*% sp$Ht %*% P_tl[[j]] #Final estiamte of the covariance matrix
//
//   #Reinitialize for the next iteration
//   B_LL = B_tt[, j]
//   P_LL = P_tt[[j]]
// }
// Sys.time() - st
//
// #Kalman Smoother
// for(j in (ncol(yt) - 1):1){
//   B_tt[, j] = B_tt[, j] + P_tt[[j]] %*% t(sp$Ft) %*% ginv(P_tl[[j + 1]]) %*% (B_tt[, j + 1] - B_tl[, j + 1])
//   P_tt[[j]] = P_tt[[j]] + P_tt[[j]] %*% t(sp$Ft) %*% ginv(P_tl[[j + 1]]) %*% (P_tt[[j + 1]] - P_tl[[j + 1]]) %*% t(P_tt[[j]] %*% t(sp$Ft) %*% ginv(P_tl[[j + 1]]))
// }

// [[Rcpp::export]]
Rcpp::List kalman_filter(const arma::mat B0, const arma::mat P0, const arma::mat Dt,
                         const arma::mat At, const arma::mat Ft, const arma::mat Ht,
                         const arma::mat Qt, const arma::mat Rt, const arma::mat yt,
                         const arma::mat X, const arma::mat beta){
  
  //Initialize the filter
  arma::mat B_LL = B0;
  arma::mat P_LL = P0;
  double lnl = 0.0;
  int n_cols = yt.n_cols;
  int n_rows = yt.n_rows;
  
  //Define the storage matrices
  arma::mat B_tt(Ft.n_rows, n_cols);
  arma::mat B_tl(B_tt.n_rows, n_cols);
  arma::cube P_tt(Ft.n_rows, Ft.n_rows, n_cols);
  arma::cube P_tl(Ft.n_rows, Ft.n_rows, n_cols);
  arma::mat N_t(n_rows, n_cols);
  arma::cube F_t(n_rows, n_rows, n_cols);
  arma::cube K_t(Ft.n_rows, n_rows, n_cols);
  arma::uvec na_idx;
  
  //Define some matrix transforms
  arma::mat Ft_t = Ft.t();
  arma::mat Ht_t = Ht.t();
  
  //Kalman filter routine
  for(int i = 0; i < n_cols; i++){
    
    //Initial estimates conditional on t-1
    B_tl.col(i) = Dt + Ft * B_LL; //Initial estimate of unobserved values conditional on t-1
    P_tl.slice(i) = Ft * P_LL * Ft_t + Qt; //Initial estimate of the covariance matrix conditional on t-1
    N_t.col(i) = yt.col(i) - At - Ht * B_tl.col(i) - beta * X.col(i); //Prediction error conditional on t-1
    F_t.slice(i) = Ht * P_tl.slice(i) * Ht_t + Rt; //Variance of the prediction error conditional on t-1
    K_t.slice(i) = P_tl.slice(i) * Ht_t * inv(F_t.slice(i)); //Kalman gain conditional on t-1
      
    //Find any missing values and replace them with 0 so final estimates will be the same as the initial estimates for this iteration
    // na_idx = arma::find_nonfinite(N_t.col(i));
    // if(!na_idx.is_empty()){
    //   cols = i;
    //   N_t(na_idx, cols) = arma::vec(na_idx.n_elem, arma::fill::zeros);
    // }
    na_idx = arma::find_nonfinite(N_t.col(i));
    if(!na_idx.is_empty()){
      B_tt.col(i) = B_tl.col(i);
      P_tt.slice(i) = P_tl.slice(i);
      K_t.slice(i) = arma::zeros(K_t.slice(i).n_rows, K_t.slice(i).n_cols);
      F_t.slice(i) = arma::ones(F_t.slice(i).n_rows, F_t.slice(i).n_cols)*R_PosInf;
    }else{
      //Final estimates conditional on t
      B_tt.col(i) = B_tl.col(i) + K_t.slice(i) * N_t.col(i); //Final estimate of the unobserved values
      P_tt.slice(i) = P_tl.slice(i) - K_t.slice(i) * Ht * P_tl.slice(i); //Final estimate of the covariance matrix
      lnl = lnl + 0.5*arma::as_scalar((log(det(F_t.slice(i))) +  N_t.col(i).t() * inv(F_t.slice(i)) * N_t.col(i)));
    }
    
    //Reinitialize for the next iteration
    B_LL = B_tt.col(i);
    P_LL = P_tt.slice(i);
  }
  
  return Rcpp::List::create(Rcpp::Named("loglik") = -lnl,
                            Rcpp::Named("B_tl") = B_tl,
                            Rcpp::Named("B_tt") = B_tt,
                            Rcpp::Named("P_tl") = P_tl,
                            Rcpp::Named("P_tt") = P_tt,
                            Rcpp::Named("F_t") = F_t,
                            Rcpp::Named("N_t") = N_t,
                            Rcpp::Named("K_t") = K_t);
}

// [[Rcpp::export]]
Rcpp::List kalman_smoother(const arma::mat B_tl, arma::mat B_tt, const arma::cube P_tl,
                           arma::cube P_tt, const arma::mat Ft){
  int t = B_tt.n_cols - 1;
  arma::mat Ft_t = Ft.t();
  arma::mat Ptt_x_Ft_x_PtInv = P_tt.slice(t - 1) * Ft_t * Rginv(P_tl.slice(t));
  
  for(int i = t - 1; i >= 0; i--){
    Ptt_x_Ft_x_PtInv = P_tt.slice(i) * Ft_t * Rginv(P_tl.slice(i + 1));
    B_tt.col(i) = B_tt.col(i) + Ptt_x_Ft_x_PtInv * (B_tt.col(i + 1) - B_tl.col(i + 1));
    P_tt.slice(i) = P_tt.slice(i) + Ptt_x_Ft_x_PtInv * (P_tt.slice(i + 1) - P_tl.slice(i + 1)) * Ptt_x_Ft_x_PtInv.t();
  }
  return Rcpp::List::create(Rcpp::Named("B_tt") = B_tt,
                            Rcpp::Named("P_tt") = P_tt);
}

//RcppArmadillo.package.skeleton(name = "tcs", path = "path")
//compileAttributes(verbose=TRUE)
//library(tools)
//package_native_routine_registration_skeleton("path")
//git config remote.origin.url git@github.com:user/autostsm.git

//Rcpp::sourceCpp("src/kalmanfilter.cpp")


