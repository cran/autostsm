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

// [[Rcpp::export]]
arma::mat gen_inv(arma::mat& m){
  arma::mat out(m.n_rows, m.n_cols);
  try{
    out = inv(m);
  }catch(std::exception &ex){
    out = Rginv(m);
  }
  return out;
}

// [[Rcpp::export]]
Rcpp::List kalman_filter(Rcpp::List& sp, const arma::mat& yt, 
                         const arma::mat& Xo, const arma::mat& Xs, 
                         bool smooth = false){
  
  //Initialize matrices
  arma::mat B0 = sp["B0"];
  arma::mat P0 = sp["P0"];
  arma::mat Dm = sp["Dm"];
  arma::mat Am = sp["Am"];
  arma::mat Fm = sp["Fm"];
  arma::mat Hm = sp["Hm"];
  arma::mat Qm = sp["Qm"];
  arma::mat Rm = sp["Rm"];
  arma::mat betaO = sp["betaO"];
  arma::mat betaS = sp["betaS"];
  
  //Initialize the filter
  arma::mat B_LL = B0;
  arma::mat P_LL = P0;
  double lnl = 0.0;
  int n_cols = yt.n_cols;
  int n_rows = yt.n_rows;
  
  //Define the storage matrices
  arma::mat B_tt(Fm.n_rows, n_cols);
  arma::mat B_tl(B_tt.n_rows, n_cols);
  arma::cube P_tt(Fm.n_rows, Fm.n_rows, n_cols);
  arma::cube P_tl(Fm.n_rows, Fm.n_rows, n_cols);
  arma::mat N_t(n_rows, n_cols);
  arma::cube F_t = arma::ones(n_rows, n_rows, n_cols)*R_PosInf;
  arma::cube K_t(Fm.n_rows, n_rows, n_cols);
  arma::uvec non_na_idx;
  arma::uvec iv;
  arma::mat yt_pred(yt.n_rows, yt.n_cols);
  
  //Define some matrix transforms
  arma::mat Fm_t = Fm.t();
  arma::mat Hm_t = Hm.t();
  arma::mat F_t_inv = arma::zeros(F_t.n_rows, F_t.n_cols);
  arma::mat F_t_submat;
  
  //Kalman filter routine
  for(int i = 0; i < n_cols; i++){
    //Find the non-missing values
    non_na_idx = arma::find_finite(yt.col(i));
    
    //Initial estimates conditional on t-1
    //Initial estimate of unobserved values conditional on t-1
    B_tl.col(i) = Dm + Fm * B_LL + betaS * Xs.col(i); 
    //Initial estimate of the covariance matrix conditional on t-1
    P_tl.slice(i) = Fm * P_LL * Fm_t + Qm; 
    //Prediction error conditional on t-1
    N_t.col(i) = yt.col(i) - Am - Hm * B_tl.col(i) - betaO * Xo.col(i); 
    F_t_inv = arma::zeros(F_t.slice(i).n_rows, F_t.slice(i).n_cols);
    if(!non_na_idx.is_empty()){
      //Variance of the prediction error conditional on t-1
      F_t.slice(i).submat(non_na_idx, non_na_idx) = Hm.rows(non_na_idx) * P_tl.slice(i) * Hm_t.cols(non_na_idx) + Rm.submat(non_na_idx, non_na_idx);
      F_t_submat = F_t.slice(i).submat(non_na_idx, non_na_idx);
      F_t_inv.submat(non_na_idx, non_na_idx) = gen_inv(F_t_submat);
    }
    //Kalman gain conditional on t-1
    K_t.slice(i) = P_tl.slice(i) * Hm_t * F_t_inv; 
    
    //Final estimates conditional on t
    if(!non_na_idx.is_empty()){
      iv << i;
      //Final estimate of the unobserved values
      B_tt.col(i) = B_tl.col(i) + K_t.slice(i).cols(non_na_idx) * N_t.submat(non_na_idx, iv); 
      //Final estimate of the covariance matrix
      P_tt.slice(i) = P_tl.slice(i) - K_t.slice(i).cols(non_na_idx) * Hm.rows(non_na_idx) * P_tl.slice(i); 
      //Update the log likelihood
      lnl = lnl + 0.5*arma::as_scalar((log(det(F_t.slice(i).submat(non_na_idx, non_na_idx))) + N_t.submat(non_na_idx, iv).t() * F_t_inv.submat(non_na_idx, non_na_idx) * N_t.submat(non_na_idx, iv)));
    }else{
      B_tt.col(i) = B_tl.col(i);
      P_tt.slice(i) = P_tl.slice(i);
    }
    
    //Get the finalized predictions
    yt_pred.col(i) = Am + Hm * B_tt.col(i) + betaO * Xo.col(i); 
    
    //Reinitialize for the next iteration
    B_LL = B_tt.col(i);
    P_LL = P_tt.slice(i);
  }
  
  if(smooth == true){
    int t = B_tt.n_cols - 1;
    arma::mat Ptt_x_Ft_x_PtInv = P_tt.slice(t - 1) * Fm_t * gen_inv(P_tl.slice(t));
    
    //Update the estimates of B_tt and P_tt using future information
    for(int i = t - 1; i >= 0; i--){
      Ptt_x_Ft_x_PtInv = P_tt.slice(i) * Fm_t * gen_inv(P_tl.slice(i + 1));
      B_tt.col(i) = B_tt.col(i) + Ptt_x_Ft_x_PtInv * (B_tt.col(i + 1) - B_tl.col(i + 1));
      P_tt.slice(i) = P_tt.slice(i) + Ptt_x_Ft_x_PtInv * (P_tt.slice(i + 1) - P_tl.slice(i + 1)) * Ptt_x_Ft_x_PtInv.t();
      
      //Get the finalized predictions
      yt_pred.col(i) = Am + Hm * B_tt.col(i) + betaO * Xo.col(i); 
    }
  }
  
  
  
  return Rcpp::List::create(Rcpp::Named("loglik") = -lnl,
                            Rcpp::Named("yt_pred") = yt_pred,
                            Rcpp::Named("B_tl") = B_tl,
                            Rcpp::Named("B_tt") = B_tt,
                            Rcpp::Named("P_tl") = P_tl,
                            Rcpp::Named("P_tt") = P_tt,
                            Rcpp::Named("F_t") = F_t,
                            Rcpp::Named("N_t") = N_t,
                            Rcpp::Named("K_t") = K_t);
}

//RcppArmadillo.package.skeleton(name = "tcs", path = "path")
//compileAttributes(verbose=TRUE)
//library(tools)
//package_native_routine_registration_skeleton("path")
//git config remote.origin.url git@github.com:user/autostsm.git

//Rcpp::sourceCpp("src/kalmanfilter.cpp")
