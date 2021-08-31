#include "RcppEigen.h"
#include "calc_kernel.h"
#include <cmath>


// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

inline double norm_cdf(const double x){
  return std::erfc(-x/std::sqrt(2.0))/2.0;
}


inline double plogis(const double p){
  return 1/(1+exp(-p));
}

inline double std_normal_pdf(const double p){
  return exp(-pow(p, 2.0)/2.0 - 0.5 * log(2.0 * M_PI));
}

inline double f_base(const double lp, const std::string fmly){
  double x;
  if (fmly == "gaussian"){
    x = lp;
  }else if (fmly == "logit"){
    x = plogis(lp);
  }else if (fmly == "poisson"){
    x = exp(lp);
  }else if (fmly == "probit"){
    x = norm_cdf(lp);
  }else{
    x = 0;
  }
  return x;
}

inline double f_prime(const double lp, const std::string fmly){
  double x;
  if (fmly == "gaussian"){
    x = 1;
  }else if (fmly == "logit"){
    x = plogis(lp);
    x = x * (1 - x);
  }else if (fmly == "poisson"){
    x = exp(lp);
  }else if (fmly == "probit"){
    x = std_normal_pdf(lp);
  }else{
    x = 0;
  }
  return x;
}


inline double f_double_prime(const double lp, const std::string fmly){
  double x;
  if (fmly == "gaussian"){
    x = 0;
  }else if (fmly == "logit"){
    x = plogis(lp);
    x = x * (1 - x) * (1 - 2 * x);
  }else if (fmly == "poisson"){
    x = exp(lp);
  }else if (fmly == "probit"){
    double pdf_lp = std_normal_pdf(lp);
    x = pdf_lp * -lp;
    // x = pdf_lp * (-1 + pow(lp, 2.0));
  }else{
    x = 0;
  }
  return x;
}

//' Internal C++ function to calculate marginal effects
// [[Rcpp::export]]
Rcpp::List cpp_gkrls_me(
  const Eigen::MatrixXd std_X_train,
  const Eigen::MatrixXd std_X_test,
  const double bandwidth,
  const std::string family,
  const bool mahal,
  const double sd_y,
  const Eigen::MatrixXd tS,
  const Eigen::VectorXd fe_mean,
  const Eigen::VectorXd re_mean,
  const Eigen::VectorXd all_mean,
  const Eigen::MatrixXd vcov_ridge,
  const Eigen::MatrixXd FE_matrix_test,
  const Eigen::MatrixXd W_Matrix,
  const Eigen::MatrixXd WX_test,
  const Eigen::MatrixXd WX_train,
  const Eigen::MatrixXd raw_X_test,
  const Eigen::VectorXd std_mean,
  const Eigen::MatrixXd std_whiten,
  const std::vector<std::string> type_mfx,
  const Eigen::MatrixXd fd_matrix,
  const Eigen::MatrixXd std_fd_matrix
){
  
  int N_train = std_X_train.rows();
  int N_test = std_X_test.rows();
  int SIZE_FE = FE_matrix_test.cols();
  int SIZE_KERNEL = W_Matrix.cols();
  int SIZE_MFX = SIZE_FE + SIZE_KERNEL;
  int SIZE_PARAMETER = SIZE_FE + re_mean.size();
  
  Eigen::ArrayXd Sc = (tS.transpose() * re_mean).array();
  Eigen::MatrixXd t_whiten = std_whiten.transpose();
    
  Eigen::MatrixXd ME_pointwise(N_test, SIZE_MFX);
  Eigen::MatrixXd ME_pointwise_var(N_test, SIZE_MFX);
  
  Eigen::MatrixXd AME_grad = Eigen::MatrixXd::Zero(all_mean.size(), SIZE_MFX);

  for (int i = 0; i < N_test; i++){
    
    Eigen::VectorXd std_X_i = std_X_test.row(i);
    Eigen::VectorXd WX_i = WX_test.row(i);
    Eigen::VectorXd X_FE_i = FE_matrix_test.row(i);

    Eigen::VectorXd kern_i(N_train);

    for (int j = 0; j < N_train; j++){
      kern_i(j) = kern_gauss(std_X_i, std_X_train.row(j), bandwidth);      
    }
    
    Eigen::VectorXd tilde_k_i = tS * kern_i;
    double fe_i = (X_FE_i.array() * fe_mean.array()).sum();
    double linpred_i = fe_i + (tilde_k_i.array() * re_mean.array()).sum();
    double f_prime_i = f_prime(linpred_i, family);
    double f_double_prime_i = f_double_prime(linpred_i, family);

    for (int p = 0; p < SIZE_MFX; p++){
      
      if (type_mfx[p] == "deriv_Kern"){

        
        Eigen::ArrayXd D_ip = -1.0 * WX_train.col(p - SIZE_FE).array();
        D_ip += WX_i(p - SIZE_FE);
        Eigen::VectorXd kern_ip = kern_i;
        kern_ip.array() *= D_ip;
        
        double meat_ip = (kern_ip.array() * Sc).sum();
        double kdsc = -2.0/bandwidth * meat_ip;
        
        ME_pointwise(i,p) = f_prime_i * kdsc;

        Eigen::VectorXd grad_ME_K_p_c = (kdsc * f_double_prime_i) * (tS * kern_i);
        grad_ME_K_p_c +=  (f_prime_i * -2.0/bandwidth) * tS * kern_ip;
        
        Eigen::VectorXd grad_ME_K_p_beta = (-2.0/bandwidth * f_double_prime_i * meat_ip) * X_FE_i;
        
        Eigen::VectorXd grad_ME_K_p(SIZE_PARAMETER);
        grad_ME_K_p << grad_ME_K_p_beta, grad_ME_K_p_c;

        ME_pointwise_var(i, p) = grad_ME_K_p.transpose() * (vcov_ridge * grad_ME_K_p);
        AME_grad.col(p) = AME_grad.col(p) + grad_ME_K_p;
          
      }else if (type_mfx[p] == "deriv_FE"){

        ME_pointwise(i,p) = f_prime_i * fe_mean(p);
        
        Eigen::VectorXd grad_ME_FE_p_beta(SIZE_FE);
        
        for (int p_prime = 0; p_prime < SIZE_FE; p_prime++){
          double grad_ME_FE_p_prime = f_double_prime_i * fe_mean(p) * X_FE_i(p_prime);
          if (p == p_prime){
            grad_ME_FE_p_prime += f_prime_i;
          }
          grad_ME_FE_p_beta(p_prime) = grad_ME_FE_p_prime;
        }
        
        Eigen::VectorXd grad_ME_FE_p_c = (-2.0/bandwidth * f_double_prime_i * fe_mean(p)) * 
          (tS * kern_i);
        
        Eigen::VectorXd grad_ME_FE_p(SIZE_PARAMETER);
        
        grad_ME_FE_p << grad_ME_FE_p_beta, grad_ME_FE_p_c;
        
        ME_pointwise_var(i,p) = grad_ME_FE_p.transpose() * (vcov_ridge * grad_ME_FE_p);
        AME_grad.col(p) = AME_grad.col(p) + grad_ME_FE_p;
        
      }else if (type_mfx[p] == "FD"){
        
        // Declare the counterfactual kernel and FE
        Eigen::VectorXd k_i_0 = kern_i;
        Eigen::VectorXd k_i_1 = kern_i;
        double fe_i_0 = fe_i;
        double fe_i_1 = fe_i;
        
        if (p < SIZE_FE){
          // If the FD is in the FE block        
          double fe_i_no_p = fe_i - X_FE_i(p) * fe_mean(p);
          fe_i_0 = fe_i_no_p - fd_matrix(p,0) * fe_mean(p);
          fe_i_1 = fe_i_no_p + fd_matrix(p,1) * fe_mean(p);
        }else{
          // If the FD is the Kernel block
          if (mahal == false){
            
            // If non-Mahalanobis standardization, can use a "trick" (see bigKRLS)
            // where we do not have to recompute the full (sketched) kernel
            // k_ij = exp(- ||x_{i,-p} - x_{j,-p} ||^2/b) exp(-(x_{ip}-x_{jp})^2/b)
            
            Eigen::VectorXd std_X_train_p = std_X_train.col(p - SIZE_FE);
            double std_X_test_i_p = std_X_i(p - SIZE_FE);
            
            // std_X_train_p_fd <- std_X_train[,p_fd - SIZE_FE]
            // std_X_train_i_p_df <- std_X_i[p_fd - SIZE_FE]
            // # Remove the contribution of the FD column
            // k_i_no_p_fd <-  k_i * exp( (std_X_train_i_p_df - std_X_train_p_fd)^2 / bandwidth)
            // # Create the counterfactual kernels
            // k_i_0 <- k_i_no_p_fd * exp( -(std_fd_matrix[p_fd,1] - std_X_train_p_fd)^2 / bandwidth)
            // k_i_1 <- k_i_no_p_fd * exp( -(std_fd_matrix[p_fd,2] - std_X_train_p_fd)^2 / bandwidth)
              
            Eigen::ArrayXd offset_ip = exp(pow(std_X_train_p.array() - std_X_test_i_p, 2.0)/bandwidth);
            
            Eigen::ArrayXd k_i_no_p = kern_i.array();
            k_i_no_p *= offset_ip;
            
            Eigen::ArrayXd adjust_i0 = exp(-pow(std_X_train_p.array() - std_fd_matrix(p, 0), 2.0)/bandwidth);
            Eigen::ArrayXd adjust_i1 = exp(-pow(std_X_train_p.array() - std_fd_matrix(p, 1), 2.0)/bandwidth);
              
            k_i_0.array() = k_i_no_p * adjust_i0;
            k_i_1.array() = k_i_no_p * adjust_i1;
            
          }else{
            
            Eigen::VectorXd raw_i = raw_X_test.row(i);
            if (fd_matrix(p, 0) != raw_i(p - SIZE_FE)){
              Eigen::VectorXd raw_i_0 = raw_i;
              raw_i_0(p - SIZE_FE) = fd_matrix(p, 0);
              Eigen::VectorXd std_X_i_0 = t_whiten * (raw_i_0 - std_mean);
              for (int j = 0; j < N_train; j++){
                k_i_0(j) = kern_gauss(std_X_i_0, std_X_train.row(j), bandwidth);
              }
            }
            
            if (fd_matrix(p, 1) != raw_i(p - SIZE_FE)){
              Eigen::VectorXd raw_i_1 =  raw_i;
              raw_i_1(p - SIZE_FE) = fd_matrix(p, 1);
              Eigen::VectorXd std_X_i_1 = t_whiten * (raw_i_1 - std_mean);
              for (int j = 0; j < N_train; j++){
                k_i_1(j) = kern_gauss(std_X_i_1, std_X_train.row(j), bandwidth);
              }
            }
          }
        }
        
        double linpred_i_1 = fe_i_1 + (k_i_1.array() * Sc).sum();
        double linpred_i_0 = fe_i_0 + (k_i_0.array() * Sc).sum();
        
        ME_pointwise(i, p) = f_base(linpred_i_1, family) - f_base(linpred_i_0, family);
        Eigen::VectorXd grad_ME_i_fd_beta = (f_prime(linpred_i_1, family) - f_prime(linpred_i_0, family)) * X_FE_i;
        Eigen::VectorXd grad_ME_i_fd_c = tS * (f_prime(linpred_i_1, family) * k_i_1 - f_prime(linpred_i_0, family) * k_i_0);
        Eigen::VectorXd grad_ME_FD(SIZE_PARAMETER);
        grad_ME_FD << grad_ME_i_fd_beta, grad_ME_i_fd_c;
        ME_pointwise_var(i,p) = grad_ME_FD.transpose() * (vcov_ridge * grad_ME_FD);
        AME_grad.col(p) = AME_grad.col(p) + grad_ME_FD;
      }
    }

  }  
  
  AME_grad.array() *= 1.0/N_test;
  
  Eigen::VectorXd AME_pointwise(SIZE_MFX);
  Eigen::VectorXd AME_pointwise_var(SIZE_MFX);
  
  for (int j = 0; j < SIZE_MFX; j++){
    Eigen::VectorXd AME_grad_j = AME_grad.col(j);
    AME_pointwise_var(j) = AME_grad_j.transpose() * (vcov_ridge * AME_grad_j);
    AME_pointwise(j) = ME_pointwise.col(j).mean();
    
  }    

  ME_pointwise_var *= pow(sd_y, 2.0);
  ME_pointwise *= sd_y;
  AME_pointwise_var *= pow(sd_y, 2.0);
  AME_pointwise *= sd_y;
  AME_grad *= sd_y;
  
  return Rcpp::List::create(
    Rcpp::Named("ME_pointwise") = ME_pointwise ,
    Rcpp::Named("ME_pointwise_var") = ME_pointwise_var,
    Rcpp::Named("AME_grad") = AME_grad,
    Rcpp::Named("AME_pointwise") = AME_pointwise,
    Rcpp::Named("AME_pointwise_var") = AME_pointwise_var
  );
  
}