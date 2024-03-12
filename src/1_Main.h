#if !defined(_Main_H)
#define _Main_H

#include <RcppArmadillo.h>
#include "2_DataParam.h"

class CMain {
	
 public:
   
  CMain(int) ;
	~CMain() ; //destructor
	void Initialization() ; 
	  
  arma::mat GetY_mat() ;   void SetY_mat(arma::mat) ;
  arma::mat Getmissing_flag_mat() ;   void Setmissing_flag_mat(arma::mat) ;   
  int Getmsg_level() ;  void Setmsg_level(int) ;
  std::string Getwhere_we_are() ; void Setwhere_we_are(std::string) ; 
  
  arma::vec Getvec_HyperParameters() ;   void Setvec_HyperParameters(arma::vec) ;
  arma::vec Getmin_Y_obs() ;   void Setmin_Y_obs(arma::vec) ;
	arma::vec Getmax_Y_obs() ;   void Setmax_Y_obs(arma::vec) ;
	
  void Iterate() ; void Run(int) ; 
  
  arma::mat GetMu() ; void SetMu(arma::mat) ;
  arma::cube Getcube_UT_cholSigma() ; void Setcube_UT_cholSigma(arma::cube) ;
  arma::vec Getlogpi() ; void Setlogpi(arma::vec) ;
  double Getalpha() ; void Setalpha(double) ;
  arma::vec GetZ_vec() ; void SetZ_vec(arma::vec) ;
  arma::vec Getvec_Phi() ; void Setvec_Phi(arma::vec) ;
  
  double Geth_0() ; void Seth_0(double) ;
  double Geta_Phi() ; void Seta_Phi(double) ;
  double Getb_Phi() ; void Setb_Phi(double) ;
  double Geta_alpha() ; void Seta_alpha(double) ;
  double Getb_alpha() ; void Setb_alpha(double) ;
  
  double test_log_dMVN_fn(arma::vec, arma::vec, arma::mat) ;
  double test_log_dMVN_UT_chol_fn(arma::vec, arma::vec, arma::mat) ;
  arma::vec test_rMVN_fn(arma::vec, arma::mat) ; 
  arma::vec test_rMVN_UT_chol_fn(arma::vec, arma::mat) ;
  arma::mat test_rIW_fn(int, arma::mat) ; 
  int test_rdiscrete_fn(arma::vec) ; 
    
 private:
   
   CData Data ; CParam Param ;
   
   int IterCount ; 
   arma::vec RandVec ; 


};

#endif  //_CMain_H

