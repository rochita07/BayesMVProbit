#' Data augmentation Gibbs sampling for Multivariate Ordinal Probit model
#'
#' This function provides posterior mean and credible interval along with trace plots, density plots and carter plot for parameters using data augmentation algorithm (Albert Chib method) for posterior sampling in the Multivariate Ordinal Probit
#' 
#' @param category vector of categories for each variable, for any variable the no of category should be atleast 3 to consider this model
#' @param df_t vector of degrees of freedom for proposed student's t distribution delta (as mentioned in the paper, link attached in vignette) for each variable 
#' @param iter scalar, number of iteration for posterior to be calculated, default is 2000
#' @param burn scalar, number of burn in for posterior to be calculated, defult is 500
#' @param cred_level scaler, (should be in [0,1]) specifies at what level credible interval to be calculated, default value is 0.95
#' @param x_list must be user defined as list, each level of the list represents covariates for each level. Each variable must have atleast one covariate, subject may have level specific covariates othwerwise same value of covariates to be repeated for each ordinal variable 
#' @param sig covariance matrix of error vector, must be symmetric positive definite matrix
#' @param y must be user defined as matrix, each column should represent for each subject 
#' @param prior_delta_mean must be user defined as list, each level is vector of prior mean for delta, by default it takes value 0  as prior mean for all delta at a specific level
#' @param prior_delta_var must be user defined as list, each level is a square positive definite matrix of prior variance for delta, default is Identity matrix of proper dimension
#' @param prior_beta_mean vector of prior mean for beta , by default it takes value 0  as prior mean for all beta
#' @param prior_beta_var a square positive definite matrix of prior variance for beta, default is Identity matrix of proper dimension
#'
#' @return
#' Posterior_mean     : posterior mean of beta in vector form indicating the ordinal measure and the corresponding covariate of the ordinal measure \cr
#' \cr
#' Credible_interval  : credible interval for beta vector at specified level \cr
#' \cr
#' trace_plot         : a plot of iterations vs. sampled values for each variable in the chain, with a separate plot per variable \cr
#' \cr
#' density_plot       : the distribution of variables, with a separate plot per variable \cr
#' \cr
#' carter_plot        : plots of credible intervals for parameters from an MCMC simulation \cr
#' @export
#'
#'
#'
#' @examples
#' 
#' # Data Generation:
#' 
#' variable_dim = 2 # no of variables, i.e dimension of yi (can be obtained if y is given)
#' category = c(4, 3) ## no of categories(levels) for each variable
#' covariate_num = c(2, 3) # No of covariates for each level
#' nu_tot = category + 1 ## no of cut points 
#' nu = nu_tot - 3 ## no of cutpoints to be considered for each variable ( nu0, nu1, nuJ are fixed)
#' nu_cutoff1 = c(0, 2 , 2.5)  # of length nu + 1
#' nu_cutoff2 = c(0, 4.5)
#' beta_dim = sum(covariate_num) # dimension of beta vector
#' 
#' n = 50 # no of subjects
#' 
#' # Generation of beta
#' 
#' beta_act1 = rep(2, covariate_num[1])
#' beta_act2 = rep(2,covariate_num[2])
#' beta_act = c(beta_act1, beta_act2)
#' 
#' set.seed(1287)
#' 
#' # Generation of x
#' 
#' x1 = matrix(rnorm(covariate_num[1]* n), nrow = covariate_num[1])  # each column is for each subject
#' x2 = matrix(rnorm(covariate_num[2]* n), nrow = covariate_num[2])  # each column is for each subject
#' x_list = list(x1, x2)
#' 
#' x = lapply(1:n, function(x) matrix(rep(0, variable_dim * beta_dim), nrow = variable_dim) ) ## initialization
#' a = matrix(rep(0, variable_dim * beta_dim), nrow = variable_dim) ## to store the value in the loop
#' col_indi = c(0, cumsum(covariate_num))  # of length variable_dim + 1
#' for(i in 1:n)
#' {
#'   for(l in 1: variable_dim)
#'   {
#'     a[l, (col_indi[l] + 1) : col_indi[l + 1]] = t(x_list[[l]][,i])
#'   }
#'   x[[i]] = a
#' }
#' 
#' 
#' # Generation of error 
#' 
#' sig = diag(variable_dim) ## error covariance matrix of order variable_dim x variable_dim ( to be mentioned in the input)
#' 
#' e = mvtnorm::rmvnorm(n = n, mean = rep(0, variable_dim) , sigma = sig)
#' e = t(e) # # each column is for each subject
#' 
#' 
#' # Generation of  Z
#' 
#' z = matrix(rep(0, variable_dim * n), nrow = variable_dim)  # matrix of order variable_dim x n # # each column is for each subject
#' for(i in 1:n)
#' {
#'   z[,i] =  x[[i]] %*% beta_act + e[, i]
#' }
#' 
#' # Generation of Y 
#' 
#' cutoff1 = c(-Inf, nu_cutoff1, Inf)  ## To also consider the fixed cutoffs 
#' cutoff2 = c(-Inf, nu_cutoff2, Inf)
#' cutoff = list(cutoff1, cutoff2)
#' 
#' y = matrix(rep(0, variable_dim * n), nrow = variable_dim) ## initialization of matrix of order variable_dim x n # each column is for each subject
#' 
#' for(i in 1:n)
#' {
#'   for(l in 1:variable_dim)
#'   {
#'     for(j in 1: length(cutoff[[l]]))
#'     {
#'       if(z[l, i] > cutoff[[l]][j] && z[l, i] <= cutoff[[l]][j + 1]) # making one side inclusive
#'         y[l, i] = j
#'     }
#'   }
#' }
#' 
#' y ## example of input 'y' (user should provide)
#' 
#' 
#' # Example 1
#' 
#' ordinal_post_beta(category = c(4, 3), df_t = NULL, iter = 2000,  burn = 500, cred_level = 0.95, x_list, sig = diag(2), y , prior_delta_mean = NULL, prior_delta_var = NULL, prior_beta_mean = NULL, prior_beta_var = NULL)
#' 
#' 
#' # Example 2
#' 
#'  prior_delta_mean = list()  # Prior on delta
#'  for(i in 1: variable_dim)
#'  {
#'    prior_delta_mean[[i]] = rep(0, nu[i])
#'  }
#'  
#'  prior_delta_var = list()   # Prior on delta
#'  for(i in 1: variable_dim)
#'    {
#'      prior_delta_var[[i]] = diag(nu[i])
#'    }
#' 
#' prior_beta_mean = rep(1, beta_dim)  # Prior on beta
#' prior_beta_var = diag(2, beta_dim)  # Prior on beta
#' 
#' ordinal_post_beta(category = c(4, 3), df_t = NULL, iter = 2000,  burn = 500, cred_level = 0.95, x_list, sig = diag(2), y , prior_delta_mean, prior_delta_var, prior_beta_mean, prior_beta_var)
#' 
#' 
#' # Interpretation of indices of beta
#' # Indices for a speficific beta indicate the ordinal measure and the corresponding covariate of the ordinal measure
#' 
#' 
#' # Warnings: There may be warning message "Error in is.positive.definite(beta_post_var) : could not find function "is.positive.definite""
#' 
#' # This is beacuse of mvtnorm::rmvnorm prefers sigma to be symmetric matrix. But, covariance matrix of beta not necessarily always would be symmetric matrix.
#' # Though this warnings does not affect sampling from mvtnorm::rmvnorm
#' 
#' # For example we define covariance matrix of  posterior beta from the given example as follow: 
#' sum_t_x_sig_x = matrix(rep(0, beta_dim^2), nrow = beta_dim) ## initialization
#' for(i in 1: n)
#' {
#'   sum_t_x_sig_x = sum_t_x_sig_x + t(x[[i]]) %*% sig %*% x[[i]]
#' }
#' sum_t_x_sig_x 
#' 
#' # Considering N(0,I) prior on beta we calculate covariance matrix of posterior beta as :
#' beta_post_var = solve(solve(diag(beta_dim)) + sum_t_x_sig_x )  # variance of posterior of beta
#' beta_post_var
#' 
#' # For a single sample the following code runs with out any error
#' mvtnorm::rmvnorm(n, mean = rep(0, nrow(beta_post_var)), beta_post_var)
#' 
#' 


ordinal_post_beta = function(category, df_t = NULL, iter = 5000, burn = 1000, cred_level = 0.95, x_list, sig, y,
                             prior_delta_mean = NULL, prior_delta_var = NULL, prior_beta_mean = NULL, prior_beta_var = NULL)
{
  
  n = ncol(y) ## no of subjects
  variable_dim = nrow(y) # no of variables
  
  # Check  for numbers of categories for each variable should be  >= 3
  
  check_category = category - rep(3, variable_dim ) ## should be >=0
  check_category = ifelse(check_category < 0 , 1, 0)  # should be all zeros
  if(identical(check_category, rep(0, variable_dim)) == FALSE)
  {
    stop("numbers of categories for each variable should be  >= 3")
  }
  
  nu_tot = category + 1 ## no of cut points 
  nu = nu_tot - 3 ## no of cutpoints to be considered for each variable ( nu0, nu1, nuJ are fixed)
  
  #check whether x_list is given as list
  if(is.list(x_list) == FALSE){stop("x_list should be given as list")}
  
  dim_x_list = lapply(x_list, function(x) dim(x))  # extracting #row from each xi as each col is for each subject
  
  covariate_num = rep(0, variable_dim)  # to find number of covariates in each ordinal level
  
  for(i in 1: variable_dim)
  {
    covariate_num[i] = dim_x_list[[i]][1]
  }
  
  
  #check whether any element in covariate_num is zero
  if(min(covariate_num) == 0){stop("For each variable there should be atleast one covariate ")}
  
  beta_dim = sum(covariate_num)  # dimension of beta vector
  
  #check  # check for compatibility of dimension of x_list with n
  check_n = array(variable_dim)
  for(i in 1: variable_dim)
  {
    check_n[i] = dim_x_list[[i]][2]
  }
  
  if(identical(check_n, rep(n, variable_dim)) == FALSE)
  {
    stop("Each level in x_list should have columns eqaul to given no of subjects !")
  }
  
  ## Calculation of X list, each list is for each subject
  
  x = lapply(1:n, function(x) matrix(rep(0, variable_dim * beta_dim), nrow = variable_dim) ) ## initialization
  a = matrix(rep(0, variable_dim * beta_dim), nrow = variable_dim) ## to store the value in the loop
  col_indi = c(0, cumsum(covariate_num))  # consider the starting point to run in the loop which is of length variable_dim + 1
  for(i in 1:n)
  {
    for(l in 1: variable_dim)
    {
      a[l, (col_indi[l] + 1) : col_indi[l + 1]] = t(x_list[[l]][,i])
    }
    x[[i]] = a
  }
  
  
  #check 
  if(is.null(df_t) == TRUE)  # if not supplied default value is considered as 5 for all levels
  {
    df_t = rep(5 , variable_dim) 
  }
  
  #check
  if(iter < burn){stop("# iteration should be > # of burn in")}
  
  #check
  if(n < beta_dim ) {stop("n should be greater than dimension of beta")}
  
  #check
  if(cred_level<0 || cred_level>1){stop("credible interval level should be in [0,1]")}
  
  
  #check
  if(is.null(prior_delta_mean) == TRUE)
  {
    prior_delta_mean = list()
    for(i in 1: variable_dim)
    {
      prior_delta_mean[[i]] = rep(0, nu[i])
    }
  }
  
  if(is.null(prior_delta_var) == TRUE)
  {
    prior_delta_var = list()
    for(i in 1: variable_dim)
    {
      prior_delta_var[[i]] = diag(nu[i])
    }
  }
  
  if(is.null(prior_beta_mean) == TRUE)
  {
    prior_beta_mean = rep(0, beta_dim)
  }
  
  if(is.null(prior_beta_var) == TRUE)
  {
    prior_beta_var = diag(beta_dim)
  }
  
  
  #check
  if(is.list(prior_delta_mean) == FALSE){stop("prior_delta_mean should be given as list")}
  
  dim_prior_delta_mean = lapply(prior_delta_mean, function(x) length(x))
  for(i in 1: variable_dim)
  {
    if(dim_prior_delta_mean[[i]] != nu[i]) 
    {stop("For each level, length of delta would be category - 2")
    }
  }
  
  
  #check
  if(is.list(prior_delta_var) == FALSE){stop("prior_delta_mean should be given as list")}
  
  dim_prior_delta_var = lapply(prior_delta_var, function(x) dim(x))
  for(i in 1: variable_dim)
  {
    if(dim_prior_delta_var[[i]][1] != dim_prior_delta_var[[i]][1])
    {
      stop("Prior variance matrix for each variable should be a square matrix")
    }
    if(dim_prior_delta_var[[i]][1] != nu[i])
    {
      stop("For each level, dimension of variance matrix would be (category-2) x (category-2) ")
    }
  }
  
  #check!
  if(length(prior_beta_mean) != beta_dim)
  {stop("length of prior mean for beta should be of sum of no of total covariates for each variable")
  }
  
  #check
  
  if(matrixcalc::is.positive.definite(sig) == FALSE)
  {stop("Covariance matrix of error should be symmetric positive definite matrix")
  }
  
  #check
  check_prior_delta_var = lapply(prior_delta_var, function(x) matrixcalc::is.positive.definite(x) )
  
  for(i in 1: variable_dim)
  {
    if(identical(check_prior_delta_var[[i]], TRUE) == FALSE)
      stop(paste("Prior covariance matrix of delta for", i, "th variable be symmetric positive definite matrix"))
  }
  
  #check
  if(matrixcalc::is.positive.definite(prior_beta_var) == FALSE)
  {stop("Prior covariance matrix of beta should be symmetric positive definite matrix")
  }
  
  
  ## initializations required for the iteration loop
  
  beta = rmvnorm(n = 1, mean = prior_beta_mean , sigma = prior_beta_var) # as 1st sample of beta
  beta = as.vector(beta)
  
  delta_current_list = list(0)  ##  to initialize as list of variable_dum to store variable_dim delta vectors 
  
  for(l in 1: variable_dim)
  {
    delta_current_list[[l]] = rmvnorm(n = 1, mean = prior_delta_mean[[l]] , sigma = prior_delta_var[[l]] )
  }
  delta_current_list
  
  
  beta_mat = matrix(rep(0, beta_dim * iter), nrow = iter)  # initialization to store the posterior samples of beta, each row is for each iteration
  dim(beta_mat)
  
  cutoff_update = list(0)  # initialization as list
  z_update = matrix(rep(0, variable_dim * n), nrow = variable_dim) ## initialization
  
  # To calculate conditional mean and variance for z[l,i]
  given.ind = 1: variable_dim
  
  
  # Calculation of variance of posterior of beta
  
  sum_t_x_sig_x = matrix(rep(0, beta_dim^2), nrow = beta_dim) ## initialization
  for(i in 1: n)
  {
    sum_t_x_sig_x = sum_t_x_sig_x + t(x[[i]]) %*% sig %*% x[[i]]
  }
  sum_t_x_sig_x 
  beta_post_var = solve(solve(prior_beta_var) + sum_t_x_sig_x )  # variance of posterior of beta
  
  #check
  if(min(eigen(beta_post_var)[["values"]])<0)
  {
    stop(" variance matrix for posterior beta is not positive semi definite")
  }
  
  delta_hat_proposed_t = list(0) #initialization to store ncp parameter of proposed t density
  scale_proposed_t = list(0)  #initialization to store scale matrix of proposed t density
  
  
  for (k in 1: iter) # for iteration
  {
    
    
    ## To calculate the parameters of proposal student-t density for each variable
    
    for(l in 1:variable_dim)
    {
      ## calculation of delta_hat_proposed_t
      
      y_given_beta_delta_optim = function(delta)
      {
        
        prior_delta = dmvnorm(delta, mean = prior_delta_mean[[l]], sigma = prior_delta_var[[l]] )  # prior density for delta
        
        # calculation of likelihood     
        f = 1 ## initialization
        
        for(i in 1: n)
        {
          if(y[l, i] == 1)
          {
            a = (-x[[i]][l,] %*% beta) / diag(sig)[l]
            # b = -Inf
            f = f *  pnorm(a)   # pnorm(-inf) = 0 and nu1 = 0
          }
          
          if(y[l, i] == 2)
          {
            a = ( exp(delta[1]) - (x[[i]][l,] %*% beta) )/ diag(sig)[l]
            b = (-x[[i]][l,] %*% beta) / diag(sig)[l]
            f = f * ( pnorm(a)  - pnorm(b) )
          }
          
          if(y[l, i] == category[l])
          {
            #a = Inf
            b = (sum(exp(delta)) - (x[[i]][l,] %*% beta) ) / diag(sig)[l]
            f = f * ( 1  - pnorm(b) )
          }
          
          if(nu[l] > 1) ## for only category > 3, it should enter this loop
          {
            for( j in 3 : (category[l]-1))
            {
              if(y[l, i] == j )
              {
                a = (sum(exp(delta[1:(j-1)])) -  (x[[i]][l,] %*% beta) ) / diag(sig)[l]
                b = (sum(exp(delta[1:(j-2)])) -  (x[[i]][l,] %*% beta) ) / diag(sig)[l]
                f = f * ( pnorm(a)  - pnorm(b) )
              } # if
            } # for
          } # if
          
        } # end of i 1:n
        
        optim_fn = -(f * prior_delta)  # to minimize in optim fn (-ve) is considered
        return(optim_fn)
        
      } # end of function y_given_beta_delta_optim
      
      ## updates of delta
      
      delta_hat = optim(par = rep(0, nu[l]), fn =  y_given_beta_delta_optim, method = "BFGS")[["par"]]
      delta_hat_proposed_t[[l]] = delta_hat
      
      # optimum value of delta to be used in proposed t distn as parameter
      
      
      y_given_beta_delta_hessian = function(delta)
      {
        prior_delta = dmvnorm(delta, mean = prior_delta_mean[[l]], sigma = prior_delta_var[[l]] )  # prior density for delta
        
        # calculation of likelihood     
        f = 1 ## initialization
        
        for(i in 1: n)
        {
          if(y[l, i] == 1)
          {
            a = (-x[[i]][l,] %*% beta) / diag(sig)[l]
            # b = -Inf
            f = f *  pnorm(a)   # pnorm(-inf) = 0 and nu1 = 0
          }
          
          if(y[l, i] == 2)
          {
            a = ( exp(delta[1]) - (x[[i]][l,] %*% beta) )/ diag(sig)[l]
            b = (-x[[i]][l,] %*% beta) / diag(sig)[l]
            f = f * ( pnorm(a)  - pnorm(b) )
          }
          
          if(y[l, i] == category[l])
          {
            #a = Inf
            b = (sum(exp(delta)) - (x[[i]][l,] %*% beta) ) / diag(sig)[l]
            f = f * ( 1  - pnorm(b) )
          }
          
          if(nu[l] > 1) ## for only category > 3, it should enter this loop
          {
            for( j in 3 : (category[l]-1))
            {
              if(y[l, i] == j )
              {
                a = (sum(exp(delta[1:(j-1)])) -  (x[[i]][l,] %*% beta) ) / diag(sig)[l]
                b = (sum(exp(delta[1:(j-2)])) -  (x[[i]][l,] %*% beta) ) / diag(sig)[l]
                f = f * ( pnorm(a)  - pnorm(b) )
              }
            }
          }
          
        }
        
        hessian_fn = log(f * prior_delta)
        return(hessian_fn)
      }
      
      hessian_mat = hessian(func = y_given_beta_delta_hessian , x = delta_hat)
      D_hat = solve(-hessian_mat)
      D_hat  # var- cov matrix , to be converted into scale matrix
      
      scale_proposed_t[[l]] =  D_hat * ((df_t[l]-2) / df_t[l])  # scale matrix for poposal t distn
      scale_proposed_t[[l]]
      
    } # end of loop for l in 1: variable_dim
    
   
    ## updates of delta
    
    for(l in 1: variable_dim)
    {
      delta_proposed = rmvt(n = 1, delta = delta_hat_proposed_t[[l]] , sigma = scale_proposed_t[[l]], df = df_t[l])
      
      density_delta_current = dmvt(delta_current_list[[l]], delta = delta_hat_proposed_t[[l]] , sigma = scale_proposed_t[[l]], df = df_t[l], log = FALSE )
      density_delta_proposed = dmvt(delta_proposed, delta = delta_hat_proposed_t[[l]] , sigma = scale_proposed_t[[l]], df = df_t[l], log = FALSE )
      
      prob_MH_1st_part = y_given_beta_delta_optim(delta_proposed) / y_given_beta_delta_optim(delta_current_list[[l]])
      prob_MH_2nd_part = density_delta_current / density_delta_proposed
      prob_MH = min(1, prob_MH_1st_part * prob_MH_2nd_part)
      prob_MH = ifelse(is.nan(prob_MH) == TRUE, 1, prob_MH)
      
      u = runif(1,0,1)  # to draw an independent sample 
      
      if(u <= prob_MH)
      {
        delta_current_list[[l]] = delta_proposed
      }
      if(u > prob_MH)
      {
        delta_current_list[[l]] = delta_current_list[[l]]
      }
      
      
      ## Updates of z
      cutoff_update[[l]] = c(-Inf, 0, cumsum(exp(delta_current_list[[l]])), Inf)  ## update on nu 
      
      for(i in 1: n)
      {
        for(j in 1: category[l])
        {
          if(y[l, i] == j) # j = 1, 2, 3, ..., category
          {
            lower_z = cutoff_update[[l]][j] 
            upper_z = cutoff_update[[l]][j + 1] 
            
            cond = condMVN(mean = x[[i]] %*% beta, sigma = sig, dep=l, given = given.ind[-l], X.given = rep(1, (variable_dim - 1)), check.sigma = FALSE )
            
            z_update[l,i] = rtruncnorm(n = 1, a = lower_z, b = upper_z, mean = cond[["condMean"]], sd = sqrt(cond[["condVar"]]))
            
          }
        } # for loop j : 1,..,category[l]
      } # for loop i : 1,..,n
      
    } # end of for(l in 1: variable_dim) loop
    
    
    # updates of y
    
    # y = matrix(rep(0, variable_dim * n), nrow = variable_dim) ## initialization of matrix of order variable_dim x n # each column is for each subject
    
    
    for(i in 1:n)
    {
      for(l in 1:variable_dim)
      {
        for(j in 1: length(cutoff_update[[l]]))
        {
          if(z_update[l, i] > cutoff_update[[l]][j] && z_update[l, i] <= cutoff_update[[l]][j + 1]) # making one side inclusive
            y[l, i] = j
        }
      }
    }
    
    ## Updates of beta
    
    sum_t_x_sig_z = matrix(rep(0, beta_dim), nrow = beta_dim) ## initialization for each iteration is needed
    for(i in 1: n)
    {
      sum_t_x_sig_z = sum_t_x_sig_z  + t(x[[i]]) %*% solve(sig) %*% z_update[,i] 
    }
    sum_t_x_sig_z 
    
    beta_update_mean = beta_post_var %*% (solve(prior_beta_var) %*% prior_beta_mean  + sum_t_x_sig_z)
    beta_update = mvtnorm::rmvnorm(n = 1, mean = beta_update_mean, sigma = beta_post_var )
    
    ## initialization for next iteration 
    
    beta = as.vector(beta_update)
    beta_mat[k, ] = beta
    
  } ## end of iteration loop
 
  
  ## Naming of coefficients
  
  pnames = rep(0, beta_dim)
  s = 0 # initialization
  for(l in 1: length(covariate_num))
  {
    pnames[(s+1) : (s+covariate_num[l])] = c(paste("beta[",l ,",", 1:covariate_num[l], "]", sep=""))
    s = s + covariate_num[l]
  }
  
  #Posterior mean of beta
  
  Betaout = beta_mat[-c(1:burn), ]
  colnames(Betaout) = pnames
  postmean_beta = colMeans(Betaout)

  
  # 95% Credible indervals for beta
  
  alpha = 1- cred_level
  interval = apply(Betaout, 2, function(x) quantile(x, c((alpha/2), 1-(alpha/2))) )
  
  par_mfrow = floor(sqrt(beta_dim)) + 1  # square root for next square no of beta_dim. used in par(mfrow) to plot
  x11()
  
  #Trace plot
  par(mfrow = c(par_mfrow,par_mfrow))
  
  trace = for(i in 1 : ncol(Betaout))
  {
    traceplot(as.mcmc(Betaout[,i]), main = pnames[i])
  }
  
  x11()
  par(mfrow = c(par_mfrow,par_mfrow))
  
  #density plot
  density = for(i in 1 : ncol(Betaout))
  {
    plot(density(Betaout[,i]), col = "blue", xlab = NULL, ylab = NULL, main = pnames[i])
  }
  
  x11()
  
  # caterplot
  
  par(mfrow = c(1,1))
  
  carter = caterplot(as.mcmc(Betaout), labels.loc ="axis")
  
  
  return(list(Posterior_mean = postmean_beta , Credible_interval = interval, trace_plot = trace, density_plot = density, carter_plot = carter))
  
}  # end of function ordinal_post_beta

