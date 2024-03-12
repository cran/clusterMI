#' @importFrom methods new
#' @import Rcpp
#' @exportClass Rcpp_modelobject
#' @keywords internal
#' @useDynLib clusterMI
#' @importFrom stats sd
#' 
createModel <- function(data_obj, K_mix_comp=50){
  
  model <- new(modelobject, K_mix_comp)
  model$Y_mat = as.matrix(data_obj$Y_mat_std)
  model$.missing_flag_mat = data_obj$missing_flag_mat
	model$min_Y_obs = data_obj$min_Y_obs
	model$max_Y_obs = data_obj$max_Y_obs

  ### User defined / Hyperparameters ###
  h_0=1.0 ;  a_Phi=0.25 ;  b_Phi=0.25 ;  a_alpha=0.25 ;  b_alpha=0.25
  model$vec_HyperParameters = c( h_0, a_Phi, b_Phi, a_alpha, b_alpha )
  
  model$msg_level = 2
    # 0: errors; 1: error and warnings; 2: errors, warnings and info
  model$.Initialization() 
	return(model)
	
} # Initialize <- function


multipleImp <- function(model_obj, data_obj, n_burnin, m_Imp, interval_btw_Imp, show_iter=TRUE){
  
  int_print = min(100,floor(interval_btw_Imp/2)) ; 
  
  if ( class(model_obj)[[1]]!="Rcpp_modelobject" ) stop("model_obj needs to be prepared by 'createModel' function") ;
  n_sample = dim(model_obj$Y_mat)[[1]] ; n_var = dim(model_obj$Y_mat)[[2]]
  K = dim(model_obj$.Mu)[[2]] ; 
  
  multiple_Imp = array(0,c(m_Imp,n_sample,n_var)) 
  multiple_Z_vec = array(0,c(m_Imp,n_sample))
  total_iter = n_burnin + m_Imp * interval_btw_Imp
	
	count_iter = 0
	draw_no_occ_cluster = draw_alpha = array(0,total_iter)
	draw_weighted_mu = array(0,c(total_iter,n_var))
	draw_pi = array(0,c(total_iter,K))
	
  # Burn-in		
  if (show_iter==TRUE){
    start_time = current_time = Sys.time(); print(paste0("Current time,",current_time)) ; prev_time = current_time # 1.3.4
  }
  print( paste0("Total iteration=", (total_iter) ) )
  print("Burn-in ..................................")
  
  if ( n_burnin > int_print ){
    n_repeat_burnin = floor( n_burnin / int_print ) ; resid_n = n_burnin - n_repeat_burnin * int_print ;
    
    for (i_repeat in 1:n_repeat_burnin){
			
			# 2018-05-04, v 1.2.0
			for (i_run in 1:int_print){
				model_obj$Iterate()
				count_iter = count_iter + 1
				draw_no_occ_cluster[count_iter] = length(unique(model_obj$.Z_vec))
				draw_alpha[count_iter] = model_obj$.alpha
				Mu_mat = Mu_std = model_obj$.Mu
				emp_pi = rep(0,K) 
				for (k in 1:K){
				  Mu_mat[,k] = data_obj$mean_Y_input + data_obj$sd_Y_input * Mu_std[,k] 
				  emp_pi[k] = mean((model_obj$.Z_vec+1)==k)
				}
				# draw_weighted_mu[count_iter,] = t( Mu_mat %*% emp_pi )
				draw_weighted_mu[count_iter,] = t( Mu_mat %*% exp(model_obj$.logpi) )
				draw_pi[count_iter,] = exp(model_obj$.logpi)
			}
      # model_obj$Run(int_print)
			
      current_time = Sys.time()
      if (show_iter==TRUE){
        int_time = round( as.numeric(difftime(current_time, prev_time, units = "mins")), 1 )
        # print(paste("Current time,",current_time))
        print(paste0("Iter=",(i_repeat*int_print),", ",int_time," min for prev ", int_print, " iter"))
      }
      prev_time = current_time
    } # for (i_repeat in 1:n_repeat_burnin)
    
    if (resid_n>0){
			# 2018-05-04, v 1.2.0
			for (i_run in 1:resid_n){
				model_obj$Iterate()
				count_iter = count_iter + 1
				draw_no_occ_cluster[count_iter] = length(unique(model_obj$.Z_vec))
				draw_alpha[count_iter] = model_obj$.alpha
				Mu_mat = Mu_std = model_obj$.Mu
				emp_pi = rep(0,K) 
				for (k in 1:K){
				  Mu_mat[,k] = data_obj$mean_Y_input + data_obj$sd_Y_input * Mu_std[,k] 
				  emp_pi[k] = mean((model_obj$.Z_vec+1)==k)
				}
				# draw_weighted_mu[count_iter,] = t( Mu_mat %*% emp_pi )
				draw_weighted_mu[count_iter,] = t( Mu_mat %*% exp(model_obj$.logpi) )
				draw_pi[count_iter,] = exp(model_obj$.logpi)
			}
			# 2018-05-04, v 1.2.0
    	# model_obj$Run(resid_n) ;
    } 
			
  } else {
		# 2018-05-04, v 1.2.0
		for (i_run in 1:n_burnin){
			model_obj$Iterate()
			count_iter = count_iter + 1
			draw_no_occ_cluster[count_iter] = length(unique(model_obj$.Z_vec))
			draw_alpha[count_iter] = model_obj$.alpha
			Mu_mat = Mu_std = model_obj$.Mu
			emp_pi = rep(0,K) 
			for (k in 1:K){
			  Mu_mat[,k] = data_obj$mean_Y_input + data_obj$sd_Y_input * Mu_std[,k] 
			  emp_pi[k] = mean((model_obj$.Z_vec+1)==k)
			}
			# draw_weighted_mu[count_iter,] = t( Mu_mat %*% emp_pi )
			draw_weighted_mu[count_iter,] = t( Mu_mat %*% exp(model_obj$.logpi) )
			draw_pi[count_iter,] = exp(model_obj$.logpi)
		}
    # model_obj$Run(n_burnin)
  }
  
  # Store imputed datasets
  print("Drawing imputed datasets ............")
  
  for ( i_Imp in 1:m_Imp ){
		
		# 2018-05-04, v 1.2.0
		for (i_run in 1:interval_btw_Imp){
			model_obj$Iterate()
			count_iter = count_iter + 1
			draw_no_occ_cluster[count_iter] = length(unique(model_obj$.Z_vec))
			draw_alpha[count_iter] = model_obj$.alpha
			Mu_mat = Mu_std = model_obj$.Mu
			emp_pi = rep(0,K) 
			for (k in 1:K){
			  Mu_mat[,k] = data_obj$mean_Y_input + data_obj$sd_Y_input * Mu_std[,k] 
			  emp_pi[k] = mean((model_obj$.Z_vec+1)==k)
			}
			# draw_weighted_mu[count_iter,] = t( Mu_mat %*% emp_pi )
			draw_weighted_mu[count_iter,] = t( Mu_mat %*% exp(model_obj$.logpi) )
			draw_pi[count_iter,] = exp(model_obj$.logpi)
		}
    # model_obj$Run(interval_btw_Imp)
    imputed_Y_std = model_obj$Y_mat 
    for (i_sample in 1:n_sample){
      multiple_Imp[i_Imp,i_sample,] = data_obj$mean_Y_input + data_obj$sd_Y_input * imputed_Y_std[i_sample,]
    }
    multiple_Z_vec[i_Imp,] = model_obj$.Z_vec;
    current_time = Sys.time()
    int_time = round( as.numeric(difftime(current_time, prev_time, units = "mins")), 1 )
    # print(paste("Current time,",current_time))
    print(paste0("Iter=",(n_burnin + i_Imp*interval_btw_Imp),", Imp ", i_Imp, "/",m_Imp,". ", int_time," min. for prev ",interval_btw_Imp," iter"))
    prev_time = current_time
  } # for ( i_Imp in 1:no_Imp )
  
	Imputation_result = list( multiple_Imp = multiple_Imp, draw_no_occ_cluster = draw_no_occ_cluster, draw_alpha = draw_alpha, draw_weighted_mu = draw_weighted_mu, draw_pi=draw_pi, multiple_Z_vec = multiple_Z_vec ) # 2018-05-04, v 1.2.0
	
	# 1.3.4
	print("Finished...............................")
	if (show_iter==TRUE){
		current_time = Sys.time(); print(paste("Current time,",current_time))
		print( paste0("Total time for the ",total_iter," iterations = ", round(current_time-start_time, 3) ) )
	}
	# 1.3.4
		
  return(Imputation_result)
  
} # multipleImp <- function




