#if !defined(_Param_H)
#define _Param_H

#include <RcppArmadillo.h>
#include <sstream>      // For convert int to string

class CData {
  
public:
  CData() ; ~CData() ; // constructor, destructor
  void Initialization() ; 
  int n_sample, n_var, K, msg_level, limit_multiple_y_out, zeta_0 ; 
  double h_0, a_Phi, b_Phi, a_alpha, b_alpha ;
  arma::vec vec_HyperParameters, vec_mu_0, min_Y_obs, max_Y_obs ;
  arma::mat Y_mat, S_mat ;
  
// private:

};

class CParam {
  
public:
	CParam() ; ~CParam(); // constructor, destructor
	void Initialization(CData &Data) ; 
  void Iterate(int Iter, CData &Data) ;
  std::string where_we_are ; 
  
  double alpha ; 
  arma::vec vec_n_Z_vec, Z_vec, vec_Phi, logpi ; 
  arma::mat Mu ; 
  arma::cube cube_UT_cholSigma  ; 

private:
  CData Data ; 
  
	int in_range(arma::vec) ;
  arma::mat sub_Mean_Cov_fn(arma::vec, arma::vec, arma::mat) ;
  arma::mat mean_UT_cond_Normal_fn(arma::vec, arma::vec, arma::vec, arma::mat) ; 
  double log_dMVN_UT_chol_fn(arma::vec, arma::vec, arma::mat) ; 
  arma::vec rMVN_UT_chol_fn(arma::vec, arma::mat) ;
  arma::mat rWishart_UT_chol_fn(int, arma::mat) ; 
  arma::mat rIW_UT_chol_fn(int, arma::mat) ;
  arma::vec rDirichlet_fn(arma::vec) ; 
  int rdiscrete_fn(arma::vec) ; 
  
  int n_sample, n_var, K ;
  arma::vec RandVec, nu_short, min_Y_imp, max_Y_imp ;
  void S1_MuSigma(CData &Data) ;
  void S2_pi() ;
  void S3_Phi(CData &Data) ;
  void S4_alpha(CData &Data) ;
  void S5_Z_vec(CData &Data) ;
  void S_Impute_Y(CData &Data) ;
  
};

#endif
