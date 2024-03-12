#' @keywords internal

readData <- function( Y_in, RandomSeed=99 ){
  
  set.seed(RandomSeed)
  
  if (is.null(Y_in)) stop("Y_in is NULL") ; 
  if (!is.data.frame(Y_in)) Y_in = as.data.frame(Y_in) ; 
  var_names = dimnames(Y_in)[[2]] ; n_var = dim(Y_in)[[2]]
  
  ID_with_missing = NULL
  for (i_sample in 1:dim(Y_in)[[1]]){
    if ( sum(is.na(Y_in[i_sample,])) > 0   ) {
      ID_with_missing = c(ID_with_missing,i_sample) ;
    }
  }
  
  if ( is.null(ID_with_missing) ) stop("No missing values exist in Y_in.")
  
  #######
  Y_mat_std = array(0,c(dim(Y_in)[[1]],n_var))
  mean_Y_input = sd_Y_input = rep(0,n_var)
  min_Y_obs = max_Y_obs = rep(0,n_var)
  for (i_var in 1:n_var){
    mean_Y_input[i_var] = mean(Y_in[,i_var], na.rm=TRUE)
    sd_Y_input[i_var] = sd(Y_in[,i_var], na.rm=TRUE)
    Y_mat_std[,i_var] = (Y_in[,i_var] - mean_Y_input[i_var]) / sd_Y_input[i_var]
    min_Y_obs[i_var] = min(Y_mat_std[,i_var], na.rm=TRUE)
    max_Y_obs[i_var] = max(Y_mat_std[,i_var], na.rm=TRUE)
  }
  #######
  
  missing_flag_mat = is.na(Y_mat_std)
  
  for (i in 1:dim(Y_mat_std)[[1]]){
    for (j in 1:dim(Y_mat_std)[[2]]){
      if (is.na(Y_mat_std[i,j])) Y_mat_std[i,j] = 0 ;
    }
  }
  
  InputData = list(Y_mat_std=Y_mat_std, mean_Y_input=mean_Y_input, sd_Y_input=sd_Y_input, missing_flag_mat=missing_flag_mat, min_Y_obs=min_Y_obs, max_Y_obs=max_Y_obs)
  
  return(InputData)	
  
}