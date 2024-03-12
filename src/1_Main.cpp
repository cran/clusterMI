/////////////////////
// Version: 1.1.1, Date: 2019-03-27
//  * S5 uses the collapsed Gibbs 
//     with marginalized conditional distribution, i.e., z_i | y_i_obs, Phi 
//  * Increase numerical precision of int CParam::rdiscrete_fn(arma::vec Prob){ 
////////////////////////////

#include "1_Main.h"

///////////////////////
RCPP_MODULE(IO_module){
  
  using namespace R ;
  using namespace Rcpp ;
  
  class_<CMain>( "modelobject" )
  
    .constructor< int >()     
    
    .property("Y_mat", &CMain::GetY_mat, &CMain::SetY_mat, "Input Y")
    .property(".missing_flag_mat", &CMain::Getmissing_flag_mat, &CMain::Setmissing_flag_mat, "Where missing")
    .property("msg_level", &CMain::Getmsg_level, &CMain::Setmsg_level, "0: errors; 1: error and warnings; 2: errors, warnings and info")
    .property(".where_we_are", &CMain::Getwhere_we_are, &CMain::Setwhere_we_are, "where_we_are")
    
    .property("vec_HyperParameters", &CMain::Getvec_HyperParameters, &CMain::Setvec_HyperParameters, "h_0, a_Phi, b_Phi, a_alpha, b_alpha, n_balance")
		
		.property("min_Y_obs", &CMain::Getmin_Y_obs, &CMain::Setmin_Y_obs, "min_Y_obs")
		.property("max_Y_obs", &CMain::Getmax_Y_obs, &CMain::Setmax_Y_obs, "max_Y_obs")
    
    .method(".Initialization", &CMain::Initialization, "Initialization")
    .method("Iterate", &CMain::Iterate, "Run one iteration of MCMC algorithm")
    .method("Run", &CMain::Run, "Run multiple iterations of MCMC algorithm")

    .property(".Mu", &CMain::GetMu, &CMain::SetMu, "Mu")
    .property(".cube_UT_cholSigma", &CMain::Getcube_UT_cholSigma, &CMain::Setcube_UT_cholSigma, "cube_UT_cholSigma") // 1.3.1
    .property(".logpi", &CMain::Getlogpi, &CMain::Setlogpi, "logpi")
    .property(".alpha", &CMain::Getalpha, &CMain::Setalpha, "alpha")
    .property(".Z_vec", &CMain::GetZ_vec, &CMain::SetZ_vec, "Z_vec")
    .property(".vec_Phi", &CMain::Getvec_Phi, &CMain::Setvec_Phi, "vec_Phi")

    .property(".h_0", &CMain::Geth_0, &CMain::Seth_0, "h_0")
    .property(".a_Phi", &CMain::Geta_Phi, &CMain::Seta_Phi, "a_Phi")
    .property(".b_Phi", &CMain::Getb_Phi, &CMain::Setb_Phi, "b_Phi")
    .property(".a_alpha", &CMain::Geta_alpha, &CMain::Seta_alpha, "a_alpha")
    .property(".b_alpha", &CMain::Getb_alpha, &CMain::Setb_alpha, "b_alpha")

    ; // Do not delete ;
}       

///////////////////////
CMain::CMain(int K_) {
  Data.K = K_ ; 
} // CMain::CMain

CMain::~CMain(){
} //Destructor

void CMain::Initialization() {
  IterCount = 0 ;
  Data.Initialization() ;
  Param.Initialization(Data) ;
}

void CMain::Iterate(){
  IterCount++; // std::cout << "IterCount" << std::endl ;
  Param.Iterate(IterCount, Data) ;
}

void CMain::Run(int n_iter_){
  for (int i_iter=0; i_iter<n_iter_; i_iter++) {
    IterCount++;
    Param.Iterate(IterCount, Data) ;
  }
}

///////////////////////
arma::mat CMain::GetY_mat() { return Data.Y_mat ; }
void CMain::SetY_mat(arma::mat Y_mat_) { Data.Y_mat = Y_mat_ ; }

arma::mat CMain::Getmissing_flag_mat() { return Data.S_mat ; }
void CMain::Setmissing_flag_mat(arma::mat missing_flag_mat_) { Data.S_mat = missing_flag_mat_ ; }

int CMain::Getmsg_level() { return Data.msg_level ; }
void CMain::Setmsg_level(int msg_level_) { Data.msg_level = msg_level_ ; }

std::string CMain::Getwhere_we_are() { return Param.where_we_are ; }
void CMain::Setwhere_we_are(std::string where_we_are_) { Param.where_we_are = where_we_are_ ; }

arma::vec CMain::Getvec_HyperParameters() { return Data.vec_HyperParameters ; }
void CMain::Setvec_HyperParameters(arma::vec vec_HyperParameters_) { Data.vec_HyperParameters = vec_HyperParameters_ ; }

arma::vec CMain::Getmin_Y_obs() { return Data.min_Y_obs ; }
void CMain::Setmin_Y_obs(arma::vec min_Y_obs_) { Data.min_Y_obs = min_Y_obs_ ; }

arma::vec CMain::Getmax_Y_obs() { return Data.max_Y_obs ; }
void CMain::Setmax_Y_obs(arma::vec max_Y_obs_) { Data.max_Y_obs = max_Y_obs_ ; }

///////////////////////
arma::mat CMain::GetMu() { return Param.Mu ; }
void CMain::SetMu(arma::mat Mu_) { Param.Mu = Mu_ ; }

arma::cube CMain::Getcube_UT_cholSigma() { return Param.cube_UT_cholSigma ; }
void CMain::Setcube_UT_cholSigma(arma::cube cube_UT_cholSigma_) { Param.cube_UT_cholSigma = cube_UT_cholSigma_ ; }

arma::vec CMain::Getlogpi() { return Param.logpi ; }
void CMain::Setlogpi(arma::vec logpi_) { Param.logpi = logpi_ ; }

double CMain::Getalpha() { return Param.alpha ; }
void CMain::Setalpha(double alpha_) { Param.alpha = alpha_ ; }

arma::vec CMain::GetZ_vec() { return Param.Z_vec ; }
void CMain::SetZ_vec(arma::vec Z_vec_) { Param.Z_vec = Z_vec_ ; }

arma::vec CMain::Getvec_Phi() { return Param.vec_Phi ; }
void CMain::Setvec_Phi(arma::vec vec_Phi_) { Param.vec_Phi = vec_Phi_ ; }

double CMain::Geth_0() { return Data.h_0 ; }
void CMain::Seth_0(double h_0_) { Data.h_0 = h_0_ ; }

double CMain::Geta_Phi() { return Data.a_Phi ; }
void CMain::Seta_Phi(double a_Phi_) { Data.a_Phi = a_Phi_ ; }

double CMain::Getb_Phi() { return Data.b_Phi ; }
void CMain::Setb_Phi(double b_Phi_) { Data.b_Phi = b_Phi_ ; }

double CMain::Geta_alpha() { return Data.a_alpha ; }
void CMain::Seta_alpha(double a_alpha_) { Data.a_alpha = a_alpha_ ; }

double CMain::Getb_alpha() { return Data.b_alpha ; }
void CMain::Setb_alpha(double b_alpha_) { Data.b_alpha = b_alpha_ ; }
