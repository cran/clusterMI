#include "2_DataParam.h"

#define LOG_2_PI 1.83787706640935

//////////// CData /////////////////////

CData::CData(){ }
CData::~CData(){ } //Destructor

void CData::Initialization(){
  n_sample = Y_mat.n_rows ; n_var = Y_mat.n_cols ;
  h_0 = vec_HyperParameters(0) ; a_Phi = vec_HyperParameters(1) ; b_Phi = vec_HyperParameters(2) ; 
  a_alpha = vec_HyperParameters(3) ; b_alpha = vec_HyperParameters(4) ; 
  zeta_0  = n_var+1 ;
  vec_mu_0 = arma::zeros<arma::vec>(n_var) ; 
}

//////////// CParam /////////////////////

CParam::CParam(){ }
CParam::~CParam(){ }

void CParam::Initialization(CData &Data){
	
  n_sample = Data.n_sample ; n_var = Data.n_var ;  K = Data.K ; // Copy the integers frequently used
  arma::vec mean_Y(n_var), sd_Y(n_var) ; 
  for (int i_var=0; i_var<n_var; i_var++){
    mean_Y(i_var) = mean(Data.Y_mat.col(i_var)) ; 
    sd_Y(i_var) = stddev(Data.Y_mat.col(i_var)) ; 
  }
  Mu = arma::mat(n_var, K) ; // Note that vector mu_k is Mu.column(k)
  cube_UT_cholSigma = arma::cube(n_var, n_var, K) ;  
  for (int k=0; k< K; k++){
    Mu.col(k) = mean_Y ; 
    arma::mat UT_cholSigma_temp = arma::zeros<arma::mat>(n_var,n_var) ; 
    UT_cholSigma_temp.diag() = sd_Y ; cube_UT_cholSigma.slice(k) = UT_cholSigma_temp ;
  }
  nu_short = arma::vec(K-1) ; nu_short.fill(0.1) ;
  logpi = arma::vec(K) ; // pi = arma::vec(K) ;
  double Sum_logOneMinusVg = 0.0 ;
  for (int i=0; i<(K-1); i++) {
    double nu_k = nu_short(i) ;
    logpi(i) = log(nu_k) + Sum_logOneMinusVg ; // pi(i) = exp(logpi(i)) ;
    double one_minus_V = ( 1.0 - nu_short(i) ) ;
    Sum_logOneMinusVg = Sum_logOneMinusVg + log(one_minus_V) ;
  }
  logpi(K-1) = Sum_logOneMinusVg ; // pi(K-1) = exp(Sum_logOneMinusVg) ;
  Z_vec = arma::zeros<arma::vec>(n_sample) ; vec_n_Z_vec = arma::zeros<arma::vec>(K) ; 
  int init_n_occ_comp = 20 ; 
  if ( Data.K < init_n_occ_comp ) init_n_occ_comp = Data.K ;
  arma::vec temp_prob(init_n_occ_comp) ; temp_prob.fill(1.0/init_n_occ_comp) ; 
  for (int i_sample=0; i_sample<Data.n_sample; i_sample++){
    Z_vec(i_sample) = rdiscrete_fn(temp_prob) ;
    vec_n_Z_vec(Z_vec(i_sample)) = vec_n_Z_vec(Z_vec(i_sample)) + 1 ; 
  }
  vec_Phi = arma::vec(n_var) ; vec_Phi.fill(5.0) ; 
  alpha = 1.0 ; 
  
	arma::vec range_Y_obs = Data.max_Y_obs - Data.min_Y_obs ; 
	min_Y_imp = arma::vec(n_var) ; max_Y_imp = arma::vec(n_var) ;
	for (int i_var=0; i_var<n_var; i_var++){
		min_Y_imp(i_var) = Data.min_Y_obs(i_var) - 0.1 * range_Y_obs(i_var) ;
		max_Y_imp(i_var) = Data.max_Y_obs(i_var) + 0.1 * range_Y_obs(i_var) ;
	}
		
} // void CParam::Initialization

void CParam::Iterate(int Iter, CData &Data) {
  S1_MuSigma(Data) ;
  S2_pi() ;
  S3_Phi(Data) ;
  S4_alpha(Data) ; 
  S5_Z_vec(Data) ; 
  S_Impute_Y(Data) ; 
} // void CParam::Iterate

//////////////////////////////////////

void CParam::S1_MuSigma(CData &Data){
  where_we_are = "S1_MuSigma" ; 
  
  arma::mat bar_y_k = arma::zeros<arma::mat>(K,n_var) ; 
  for (int i=0; i<n_sample; i++){
    bar_y_k.row(Z_vec(i)) = bar_y_k.row(Z_vec(i)) + Data.Y_mat.row(i) ;
  }
  for (int k=0; k<K; k++){
    if (vec_n_Z_vec(k)>0){
      bar_y_k.row(k) = bar_y_k.row(k) / vec_n_Z_vec(k) ;
    }
  }
  arma::cube cube_S_k(n_var,n_var,K) ; cube_S_k.fill(0) ;
  for (int i=0; i<n_sample; i++){
    arma::rowvec dev_y_i = Data.Y_mat.row(i) - bar_y_k.row(Z_vec(i)) ;
    cube_S_k.slice(Z_vec(i)) = cube_S_k.slice(Z_vec(i)) + dev_y_i.t() * dev_y_i ;
  }
  for (int k=0; k<K ; k++){
    int zeta_k = Data.zeta_0 ; double h_k = Data.h_0 ;
    arma::vec mu_star_k = Data.vec_mu_0 ; arma::mat Phi_k = diagmat(vec_Phi) ;
    if ( vec_n_Z_vec(k) > 0 ){ 
      zeta_k = zeta_k + vec_n_Z_vec(k) ;
      h_k = h_k + vec_n_Z_vec(k) ;
      mu_star_k = 1.0/h_k * ( vec_n_Z_vec(k)*bar_y_k.row(k).t() + Data.h_0*Data.vec_mu_0 ) ;
      arma::vec dev_bar_y_k = bar_y_k.row(k).t() - Data.vec_mu_0 ;
      Phi_k = Phi_k + cube_S_k.slice(k) + 1.0/(1.0/vec_n_Z_vec(k)+1.0/Data.h_0) * (dev_bar_y_k*dev_bar_y_k.t()) ;
    }
    arma::mat UT_cholPhi_k = arma::chol(Phi_k) ;
    arma::mat rIW_temp = rIW_UT_chol_fn(zeta_k, UT_cholPhi_k) ;
    cube_UT_cholSigma.slice(k) = arma::chol(rIW_temp) ;
    Mu.col(k) = rMVN_UT_chol_fn( mu_star_k, (1.0/sqrt(h_k)) * cube_UT_cholSigma.slice(k) ); // Modified
  } // for (int k=1; k<=K ; k++)
  
} // void CParam::S5_MuSigma

void CParam::S2_pi(){
  where_we_are = "S2_pi" ; 
  
  nu_short.fill(0.1) ; double Sum_n_m = sum(vec_n_Z_vec) ; 
  for (int k=0; k<(K-1); k++) {
    double one_tilde = 1.0 + vec_n_Z_vec(k) ; 
    Sum_n_m = Sum_n_m - vec_n_Z_vec(k) ; 
    // start from Sum_n_m - vec_n_Z_vec_1 when k=1
    // ->  Sum_n_m - sum(vec_n_Z_vec_1 + vec_n_Z_vec_2) when k=2 -> ...
    // Sum_n_m - sum(vec_n_Z_vec_1 + ... + vec_n_Z_vec(k))) i.e. sum_{g=k+1}^K vec_n_Z_vec_g
    double alpha_tilde = alpha + Sum_n_m ;
    RandVec = Rcpp::rbeta( 1, one_tilde, alpha_tilde ) ; 
    nu_short(k) = RandVec(0) ;  
  }  
  double Sum_logOneMinusVg = 0.0 ;
  for (int k=0; k<(K-1); k++) {
    logpi(k) = log(nu_short(k)) + Sum_logOneMinusVg ;  
    Sum_logOneMinusVg = Sum_logOneMinusVg + log(1.0 - nu_short(k)) ;
  }
  logpi(K-1) = Sum_logOneMinusVg ; // pi = exp(logpi) ; 
   
} // void CParam::S6_pi()

void CParam::S3_Phi(CData &Data){
  where_we_are = "S3_Phi" ; 
  
  arma::vec vec_sum_invSigma_diag = arma::zeros<arma::vec>(n_var) ; 
  for (int k=0; k<K; k++){
    arma::mat UT_inv_Chol = cube_UT_cholSigma.slice(k).i() ; 
    arma::mat inv_Sigma_k = UT_inv_Chol * UT_inv_Chol.t() ; 
    vec_sum_invSigma_diag = vec_sum_invSigma_diag + 0.5 * inv_Sigma_k.diag() ;
  }
  for (int j=0; j<n_var; j++) {
    double a_Phi_tilde_j = Data.a_Phi + 0.5 * K * Data.zeta_0 ;	
    double b_Phi_tilde_j = Data.b_Phi + 0.5 * vec_sum_invSigma_diag(j) ;
    RandVec = Rcpp::rgamma(1, a_Phi_tilde_j, 1.0/b_Phi_tilde_j ) ; // Note that b in Rcpp::rgamma(a,b) is scale, i.e., its mean is ab, NOT a/b
    vec_Phi(j) = RandVec(0) ; 
  }
  
} // void CParam::S8_Phi

void CParam::S4_alpha(CData &Data){
  where_we_are = "S4_alpha" ; 
  
  double a_alpha_tilde = Data.a_alpha + K - 1.0 ;
  double b_alpha_tilde = Data.b_alpha - logpi(K-1) ;
  // To avoid zero alpha 
  if ( b_alpha_tilde > 10 ) b_alpha_tilde = 10 ; 
  // To avoid zero alpha 
  RandVec = Rcpp::rgamma(1, a_alpha_tilde, 1.0/b_alpha_tilde ) ; // Note that b in Rcpp::rgamma(a,b) is scale, i.e., its mean is ab, NOT a/b
  alpha = RandVec(0) ; 
  
} // void CParam::S7_alpha

void CParam::S5_Z_vec(CData &Data){
  where_we_are = "S5_Z_vec" ; 
  
  vec_n_Z_vec = arma::zeros<arma::vec>(K) ;
	
  for (int i_sample=0; i_sample < n_sample; i_sample++) {
    arma::vec logN_unnorm(K) ; arma::vec y_i_vec = (Data.Y_mat.row(i_sample)).t() ;  
    // for (int k=0; k<K; k++) {
    //   logN_unnorm(k) = logpi(k) + log_dMVN_UT_chol_fn( y_i_vec, Mu.col(k), cube_UT_cholSigma.slice(k) ) ; 
    // }
    // COLLAPSED GIBBS to increase efficiency
    arma::vec one_vec(n_var) ; one_vec.fill(1.0) ; arma::vec s_i = Data.S_mat.row(i_sample).t() ; 
		
		if ( sum(s_i) < n_var ){
			
			arma::vec which_retain_var = one_vec - s_i ; arma::vec y_i_obs(sum(which_retain_var)) ;
	    int count_i_a = 0 ;
	    for (int i_var=0; i_var<n_var; i_var++){
	      if ( which_retain_var(i_var)==1 ){
	        y_i_obs(count_i_a) = y_i_vec(i_var) ; count_i_a++ ;
	      }
	    }
	    for (int k=0; k < K; k++) {
	      arma::mat subMean_UTSigma = sub_Mean_Cov_fn(which_retain_var, Mu.col(k), cube_UT_cholSigma.slice(k) ) ;
	      arma::vec sub_Mean_k = subMean_UTSigma.col(0) ; arma::mat sub_UTSigma_k = subMean_UTSigma.cols(1,sum(which_retain_var)) ;
	      logN_unnorm(k) = logpi(k) + log_dMVN_UT_chol_fn( y_i_obs, sub_Mean_k, sub_UTSigma_k ) ; // Note that vector mu_k is Mu.column(k)
	    }
    
	    double max_logN_unnorm = logN_unnorm.max() ;
	    arma::vec pi_tilde_unnorm = arma::zeros<arma::vec>(K) ;
	    for (int k=0; k<K; k++){
	      pi_tilde_unnorm(k) = exp(logN_unnorm(k)-max_logN_unnorm) ; 
	    }
	    arma::vec pi_tilde = (1.0/sum(pi_tilde_unnorm)) * pi_tilde_unnorm ;
	    Z_vec(i_sample) = rdiscrete_fn( pi_tilde );
	    vec_n_Z_vec(Z_vec(i_sample)) = vec_n_Z_vec(Z_vec(i_sample)) + 1 ; 
			
		} else { 

	    arma::vec pi_tilde = exp(logpi) ;
	    Z_vec(i_sample) = rdiscrete_fn( pi_tilde );
					
		} // if ( sum(s_i) < n_var ) else ... 
		    
  } // for (int i_sample)
    
} // void CParam::S3_Z_vec(CData &Data)

void CParam::S_Impute_Y(CData &Data){
  where_we_are = "S_Impute_Y" ; 
	
  for (int i_sample=0; i_sample<n_sample; i_sample++){
		
    arma::vec s_i = Data.S_mat.row(i_sample).t() ;
    if ( sum(s_i) > 0 ){
      if ( sum(s_i) < n_var ){
				arma::vec y_i = Data.Y_mat.row(i_sample).t() ; int z_i = Z_vec(i_sample) ;
	      arma::vec mu_i = Mu.col(z_i) ; arma::mat UT_cholSigma_i = cube_UT_cholSigma.slice(z_i) ;
	      arma::mat temp_MAT = mean_UT_cond_Normal_fn(s_i, y_i, mu_i, UT_cholSigma_i) ; 
	      arma::vec mu_a_star = temp_MAT.col(0) ; 
	      arma::mat UT_chol_a = temp_MAT.cols(1,sum(s_i)) ; 
	      arma::vec y_a = rMVN_UT_chol_fn( mu_a_star, UT_chol_a ) ;  
	      int count_i_a = 0 ; arma::vec y_temp = y_i ;				 
	      for (int i_var=0; i_var<n_var; i_var++){
	        if ( s_i(i_var)==1 ){
	          y_temp(i_var) =  y_a(count_i_a) ;
	          count_i_a++ ;
	        } // if (s_i==1)
	      } // for(i_var)
				if ( in_range(y_temp)==1 ){
					// Avoid drawing an extreme value outside of the support of observed value
					Data.Y_mat.row( i_sample ) =  y_temp.t() ;
				} // if ( sum( min_Y_imp < y_temp ) + sum( y_temp < max_Y_imp ) == (2*n_var) )		
			} else {
				int z_i = Z_vec(i_sample) ;
				arma::vec mu_temp = Mu.col( z_i ) ; arma::mat UT_cholSigma_temp = cube_UT_cholSigma.slice( z_i ) ;
				arma::vec y_temp = rMVN_UT_chol_fn( mu_temp, UT_cholSigma_temp ) ;
				if ( in_range(y_temp)==1 ){
					// Avoid drawing an extreme value outside of the support of observed value
					Data.Y_mat.row( i_sample ) =  y_temp.t() ;
				} // if ( sum( min_Y_imp < y_temp ) + sum( y_temp < max_Y_imp ) == (2*n_var) )				
			}	// if ( sum(s_i) < n_var ) else ...
			
    } // if ( sum(s_i)>0 )
  } // for (int i_sample)
  
} // CParam::S9_Impute_Y

//////////////////////////////////////
// Hang Kim's Distribution
//////////////////////////////////////

int CParam::in_range(arma::vec y_gen){
	int is_in_range = 1 ; 	
	for (int i_var=0; i_var<n_var; i_var++){
		if (y_gen(i_var)<min_Y_imp(i_var)) is_in_range = 0 ; 
		if (y_gen(i_var)>max_Y_imp(i_var)) is_in_range = 0 ; 
	}	// for 	
	return(is_in_range); 
} // bool CParam::in_range(arma::vec )

arma::mat CParam::sub_Mean_Cov_fn(arma::vec which_retain_var, arma::vec mu_vec, arma::mat UT_Sigma_mat ){
  // Return mean and Covariance of which_retain_var where which_retain_var is a vector of indicator. 
  // i.e., which_retain_var = (1,0,1,0) will produce 2*1 mean vector and 2*2 cov matrix for 1st and 3rd var
  // The output's first column contains mu_sub and other columns contains UT_Sigma_mat_sub
  // Output_mat.col(0)  // Output_mat.cols(1,sum(which_retain_var))
  
  int p_var = which_retain_var.n_rows ; int n_a = sum(which_retain_var) ; 
  if (n_a <= 0) Rcpp::stop("n_a <= 0 in sub_Mean_Cov_fn") ;  
  
  arma::mat Sigma_mat = UT_Sigma_mat.t() * UT_Sigma_mat ; 
  arma::vec mu_a = arma::vec(n_a) ; arma::mat Sigma_aa = arma::mat(n_a,n_a) ; 
  int count_i_a = 0 ; 
  for (int i_var=0; i_var<p_var; i_var++){
    if ( which_retain_var(i_var)==1 ){
      mu_a(count_i_a) = mu_vec(i_var) ;
      int count_j_a = 0 ; // int count_j_b = 0 ;
      for (int j_var=0; j_var<p_var; j_var++){
        if ( which_retain_var(j_var)==1 ){
          Sigma_aa(count_i_a,count_j_a) = Sigma_mat(i_var,j_var) ;
          count_j_a++ ;
        } 
      } // for (j_var)
      count_i_a++ ;
    }
  }
  arma::mat UT_Sigma_aa = arma::chol(Sigma_aa) ; 
  arma::mat Output_mat = join_rows(mu_a,UT_Sigma_aa) ; 
  return(Output_mat) ; 
  
} // CParam::mean_UT_cond_Normal_fn

arma::mat CParam::mean_UT_cond_Normal_fn(arma::vec is_a_vec, arma::vec y_vec, arma::vec mu_vec, arma::mat UT_Sigma_mat ){
  // Return mean and Covariance of y_vec_a  where is_a_vec is a vector of indicator. 
  // i.e., is_a_vec = (1,0,1,0) will produce 2*1 mean vector and 2*2 cov matrix of y_vec_a = (y_1, y_3)
  // The output's first column contains mu_a_star and other columns contains UT of Sigma_aa_star 
  // cond_mu: Output_mat.col(0)  // cond_UT_Sigma: Output_mat.cols(1,sum(is_a_vec))
  
  int p_var = is_a_vec.n_rows ; int n_a = sum(is_a_vec) ; int n_b = p_var - n_a ;
  if (n_a <= 0) Rcpp::stop("n_a <= 0 in mean_UT_cond_Normal_fn") ;  
  
  arma::mat Sigma_mat = UT_Sigma_mat.t() * UT_Sigma_mat ; 
  arma::vec mu_a = arma::vec(n_a) ; arma::vec mu_b = arma::vec(n_b) ; arma::vec y_b = arma::vec(n_b) ;
  arma::mat Sigma_aa = arma::mat(n_a,n_a) ; arma::mat Sigma_ab = arma::mat(n_a,n_b) ; arma::mat Sigma_bb = arma::mat(n_b,n_b) ;
  int count_i_a = 0 ; int count_i_b = 0 ;
  for (int i_var=0; i_var<p_var; i_var++){
    if ( is_a_vec(i_var)==1 ){
      mu_a(count_i_a) = mu_vec(i_var) ;
      int count_j_a = 0 ; int count_j_b = 0 ;
      for (int j_var=0; j_var<p_var; j_var++){
        if ( is_a_vec(j_var)==1 ){
          Sigma_aa(count_i_a,count_j_a) = Sigma_mat(i_var,j_var) ;
          count_j_a++ ;
        } else {
          Sigma_ab(count_i_a,count_j_b) = Sigma_mat(i_var,j_var) ;
          count_j_b++ ;
        } // if (is_a_vec(j_var)) else ...
      } // for (j_var)
      count_i_a++ ;
    } else {
      mu_b(count_i_b) = mu_vec(i_var) ; y_b(count_i_b) = y_vec(i_var) ;
      int count_j_b = 0 ;
      for (int j_var=0; j_var<p_var; j_var++){
        if ( is_a_vec(j_var)==0 ){
          Sigma_bb(count_i_b,count_j_b) = Sigma_mat(i_var,j_var) ;
          count_j_b++ ;
        } // if (is_a_vec(j_var))
      } // for (j_var)
      count_i_b++ ;
    } // if (is_a_vec==1) else ...
  } // for (i_var)
  
  arma::mat Sigma_bb_inv = Sigma_bb.i() ;
  arma::vec mu_output = mu_a + Sigma_ab * Sigma_bb_inv * (y_b-mu_b) ;
  arma::mat Sigma_a_star = Sigma_aa - Sigma_ab * Sigma_bb_inv * Sigma_ab.t() ;
  arma::mat UT_chol_output = arma::chol(Sigma_a_star) ;
  
  arma::mat Output_mat = join_rows(mu_output,UT_chol_output) ; 
  return(Output_mat) ; 
  
} // CParam::mean_UT_cond_Normal_fn

double CParam::log_dMVN_UT_chol_fn(arma::vec x, arma::vec mu, arma::mat UT_chol){
  arma::mat inv_LTchol = UT_chol.i().t() ; // arma::trans( arma::inv( UT_chol )) ;
  int xdim = x.n_rows;
  double constants = -(xdim/2) * std::log(2.0 * M_PI);
  double sum_log_inv_LTchol = arma::sum(log(inv_LTchol.diag()));
  arma::vec z = inv_LTchol * ( x - mu ) ;
  double logout = constants + sum_log_inv_LTchol - 0.5 * arma::sum(z%z) ; 
  // %	Schur product: element-wise multiplication of two objects 
  // i.e. arma::cum(z%z) = z.t() * z = sum_i z_i^2  // z.t() * z  produces 1 by 1 matrix, so need another line 
  return(logout);
  // This function is checked with "dmvnorm" in mvtnorm package on 2018/01/26
} // arma::vec dmvnrm_arma_mc

arma::vec CParam::rMVN_UT_chol_fn(arma::vec mu, arma::mat UT_chol){
  int n_var = UT_chol.n_rows ; 
  RandVec = Rcpp::rnorm(n_var,0,1) ; 
  arma::mat LT_chol = UT_chol.t() ; 
  arma::vec out = mu + LT_chol * RandVec ; 
  // MVN // y_vec = a_vec + A * z_vec where Sigma = A A^T 
  // Cholesky // Sigma = L L^T = U^T U // arma::chol(Sigma) = U
  return out ;
} // arma::mat rMVN_fn

arma::mat CParam::rWishart_UT_chol_fn( int nu, arma::mat UT_chol ){
  int p = UT_chol.n_rows ; 
  arma::mat S_mat = arma::zeros<arma::mat>(p,p) ; 
  for (int l=0; l<nu; l++){
    RandVec = Rcpp::rnorm(p,0,1) ; 
    arma::vec x_l = UT_chol.t() * RandVec ;
    S_mat = S_mat + x_l * x_l.t() ; 
  }
  return S_mat ; 
} // arma::mat CParam::rWishart_UT_chol_fn

arma::mat CParam::rIW_UT_chol_fn( int nu, arma::mat UT_chol ){
  
  // Draw IW( nu, UT_chol.t() * UT_chol ) // See 37_Matrix_RandomNumber_Rcpp.pdf ; Function checked 
  
  int p = UT_chol.n_rows ; 
  arma::mat inv_UT_chol = UT_chol.i() ;  
  arma::mat U_mat = arma::zeros<arma::mat>(p,p) ;
  for (int l=0; l<nu; l++){
    RandVec = Rcpp::rnorm(p,0,1) ;       		 
    arma::vec x_l = inv_UT_chol * RandVec ;  
    U_mat = U_mat + x_l * x_l.t() ;										 
  }
  
  arma::mat V_mat = U_mat.i() ; 	
  
  // // May result in non-positive definite matrix due to a decimal rounding error
  // arma::vec diag_lambda_vec; arma::mat Q_mat;
  // eig_sym(diag_lambda_vec, Q_mat, V_mat);
  // int count_nonpositive = 0 ;
  // while ( diag_lambda_vec.min() <= 0 ){
  //   for (int i_dim=0; i_dim<p; i_dim++){
  //     if ( diag_lambda_vec(i_dim) <= 1e-20 ) diag_lambda_vec(i_dim) = 1e-20 ;
  //   }
  //   V_mat = Q_mat * diagmat(diag_lambda_vec) * Q_mat.t() ;
  //   eig_sym(diag_lambda_vec, Q_mat, V_mat) ;
  //   count_nonpositive++ ;
  // } // while
  // if (count_nonpositive>1) std::cout << "rIW_UT_chol_fn has nonpositive matrix " <<  count_nonpositive << " times" << std::endl ;
  // 
  // arma::mat R_temp ; 
  // bool chol_success = chol(R_temp, V_mat) ; 
  // if (!chol_success){
  //   eig_sym(diag_lambda_vec, Q_mat, V_mat) ;
  //   double max_lambda = diag_lambda_vec.max() ; 
  //   for (int i_dim=0; i_dim<p; i_dim++){
  //     if ( diag_lambda_vec(i_dim) <= max_lambda * (1e-20) ) diag_lambda_vec(i_dim) = max_lambda * (1e-20) ;
  //   }
  //   V_mat = Q_mat * diagmat(diag_lambda_vec) * Q_mat.t() ;
  //   std::cout << "rIW_UT_chol_fn uses max_lambda * (1e-20)" << std::endl ;
  // } 
  
  return V_mat ; 
  
} // arma::mat CParam::rIW_UT_chol_fn 

arma::vec CParam::rDirichlet_fn( arma::vec alpha_vec ){
  int p = alpha_vec.n_rows ; arma::vec y_vec(p) ; 
  for (int j=0; j<p; j++){
    RandVec = Rcpp::rgamma(1, alpha_vec(j), 1.0 ) ; // Gamma(a_j,1)
    y_vec(j) = RandVec(0) ; 
  }
  double sum_y_vec = sum(y_vec) ; 
  arma::vec output_vec = (1.0/sum_y_vec) * y_vec ; 
  return output_vec ; 
} // arma::mat CParam::rDirichlet_fn

int CParam::rdiscrete_fn(arma::vec Prob){ 
  // generate an integer from 0 to (max_no-1) with Prob 
  if ( fabs( sum(Prob)-1.0 ) > 1e-10 ) {
    // std::cout << "Prob = " << std::endl ; 
    // std::cout << Prob.t() << std::endl ; 
    // std::cout << "sum(Prob) = " << sum(Prob) << std::endl ;
    // std::cout << "sum(Prob) != 1 in rdiscrete_fn" << std::endl ; 
    Rcpp::stop("sum(Prob) != 1 in rdiscrete_fn") ;  
  }
  int n_vec = Prob.n_rows ; 
  
  // For numerical stability, set zero for pi_k < 1e-05, i.e., one out of 100,000
  // CumProb(i) = CumProb(i-1) + Prob(i) ;  and while ( CumProb(out) < RandVec(0) ) out++ ;
  //   may generate a random value with small prob, more often than its probability
  for (int k=0; k<n_vec; k++){
    if (Prob(k)<1e-05) Prob(k)=0 ;
  }
  Prob = (1.0/sum(Prob)) * Prob ;
  
  arma::vec CumProb = Prob ; 
  for (int i=1; i<n_vec; i++) {
    CumProb(i) = CumProb(i-1) + Prob(i) ;  
  }
  RandVec = Rcpp::runif(1,0,1) ; 
  int out = 0 ; 
  while ( CumProb(out) < RandVec(0) ) out++ ;
  return out;
  // This function is checked with "dmvnorm" in mvtnorm package on 2018/01/26
} // rdiscrete_fn(arma::vec Prob)
