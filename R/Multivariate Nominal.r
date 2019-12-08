#' Title  Data augmentation Gibbs sampling for Multivariate Multinomial Probit model
#' 
#' This function provides posterior mean and credible interval along with trace plots, density plots and carter plot for parameters using data augmentation algorithm (Holmes Held method) for posterior sampling in the Multivariate Multinomial Probit
#' 
#' @param category vector of categories for each variable
#' @param iter  scalar, number of iteration for posterior to be calculated, default is 1500
#' @param burn  scalar, number of burn in for posterior to be calculated, defult is 500
#' @param cred_level scaler, (should be in [0,1]) specifies at what level credible interval to be calculated, default value is 0.95
#' @param x_list must be user defined as list, each level of the list represents covariates for each level. Each variable must have atleast one covariate, subject may have level specific covariates othwerwise same value of covariates to be repeated for each nominal variable 
#' @param sig covariance matrix of error vector, must be symmetric positive definite matrix
#' @param d must be user defined as matrix, each column should represent for each subject 
#' @param prior_beta_mean vector of prior mean for beta , by default it takes value 0  as prior mean for all beta
#' @param prior_beta_var a square positive definite matrix of prior variance for beta, default is Identity matrix of proper dimension
#'
#' @return
#' Posterior_mean     : posterior mean of beta in vector form indicating the nominal measure, level of the nominal measure and the corresponding covariate of the nominal measure \cr
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
#' @examples
#' 
#' Data Generation
#' 
#' variable_dim = 2 # no of variables
#' category = c(3,4) # no of levels for each variable
#' n = 50 # no of subjects
#' covariate_num = c(2,3) # no of covariates for each level
#' z_dim = sum(category) - length(category)  ## dimension of z for each subject
#' beta_dim = sum((category - 1) * covariate_num)
#' 
#' set.seed(1287)
#'  Generation of actual beta
#'  
#' beta_act1 = matrix(rep(2, (category[1]-1) * covariate_num[1]), nrow = category[1]-1)
#' ## matrix of regression coefficient for 1st nominal measure
#' beta_act2 = matrix(rep(3, (category[2]-1) * covariate_num[2]), nrow = category[2]-1)
#' ## matrix of regression coefficient for 1st nominal measure
#' beta_act = as.vector(c(beta_act1, beta_act2)) ## vectorized form of beta
#' ## should be of same length  as beta_dim
#' 
#' Generation of X_list
#' x1_mat = matrix(rnorm(n * covariate_num[1]), nrow = covariate_num[1] ) # each column is for each subject, each row is for each  covariates
#' x2_mat = matrix(rnorm(n * covariate_num[2]), nrow = covariate_num[2] ) # each column is for each subject, each row is for each  covariates
#' x_list = list(x1_mat, x2_mat)  # input should be given as list 
#' 
#' z1_mat = beta_act1 %*% x1_mat 
#' z2_mat = beta_act2 %*% x2_mat
#' 
#' Geneartion of  sig matrix
#' sig1 = matrix(c(1, 0.36, 0.36, 0.81), nrow = 2)
#' sig2 = matrix(c(1, 0.45, 0.24, 0.45, 0.81, 0.41, 0.24, 0.41, 0.9), nrow =3)
#' sig3 = matrix(c(rep(0.2, 6)), nrow =2)
#' sig4 = matrix(c(rep(0.2, 6)), nrow =3)
#' sig_gen = as.matrix(rbind(cbind(sig1, sig3), cbind(sig4, sig2))) ## Final sigma 
#' 
#' #Geneartion of error variable
#' e_mat = MASS::mvrnorm(n , mu = rep(0, z_dim) , Sigma = sig_gen )  # Generation of error term
#' e_mat = t(e_mat) # each column is for each subject
#' 
#' #Generation of  Z vectors for each subject, each column is for each subject
#' 
#' z_mat = matrix(rep(0, n * z_dim), ncol = n) ## each column is for each subject
#' 
#' for(i in 1: n){
#'   z_mat[, i] = c(z1_mat[, i], z2_mat[, i]) + e_mat[, i]
#' }
#' 
#' #Computation of d matrix ## user should provide this input
#' 
#' d = matrix(rep(0, n * length(category)), ncol = n)  # each col corresponds to each subject
#' q = rep(1, length(category) + 1)  ## the indicator to compute d matrix (as mentioned in the paper)
#' q[1] = 1
#' category_cum_sum = cumsum(category)
#' 
#' for(g in 2: length(q))
#' {
#'   q[g] = category_cum_sum[g-1] - (g-1 - 1)    ## ex: p1 + p2 + p3 - 2
#' }
#' 
#' for(i in 1: n)
#' { 
#'   for(j in 1:variable_dim )
#'   {
#'     x = z_mat[ q[j] : (q[j+1] - 1) , i]
#'    if(max(x) < 0)
#'     {
#'       d[j, i] = 0
#'     }
#'     if(max(x) >= 0)
#'     {
#'       d[j, i] = which.max(x)
#'     }
#'   }
#' }
#' d ## example of input 'd' (user should provide)
#' 
#' 
#' 
#' Example 1
#' 
#' nominal_post_beta(category = c(3,4), iter = 1500, burn = 500, cred_level = 0.95, x_list, sig, d, prior_beta_mean = NULL, prior_beta_var = NULL)
#' 
#' 
#' Example 2
#' 
#' prior_beta_mean = rep(1, beta_dim)  # Prior on beta
#' prior_beta_var = diag(2, beta_dim)  # Prior on beta
#' 
#' nominal_post_beta(category = c(3,4), iter = 1500, burn = 500, cred_level = 0.95, x_list, sig, d, prior_beta_mean, prior_beta_var)


nominal_post_beta = function(category, iter = 1500, burn =500, cred_level = 0.95, x_list, sig , d, prior_beta_mean = NULL, prior_beta_var = NULL)
{
  
  n = ncol(d) # no of subjects
  variable_dim = nrow(d)
  
  #check whether x_list is given as list
  if(is.list(x_list) == FALSE){stop("x_list should be given as list")}
  
  dim_x_list = lapply(x_list, function(x) dim(x))  # extracting #row from each xi as each col is for each subject
  covariate_num = rep(0, variable_dim)
  
  for(i in 1: variable_dim)
  {
    covariate_num[i] = dim_x_list[[i]][1]
  }
  
  z_dim = sum(category) - length(category)  ## dimension of z for each subject
  beta_dim = sum((category - 1) * covariate_num)
  
  #check whether any element in covariate_num is zero
  if(min(covariate_num) == 0){stop("For each variable there should be atleast one covariate ")}
  
  
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
  
  #check
  if(iter < burn){stop("# iteration should be > # of burn in")}
  
  #check
  if(n < beta_dim ) {stop("n should be greater than dimension of beta")}
  
  #check
  if(cred_level<0 || cred_level>1){stop("credible interval level should be in [0,1]")}
  
  #check
  if(is.null(prior_beta_mean) == TRUE)
  {
    prior_beta_mean = rep(0, beta_dim)
  }
  
  if(is.null(prior_beta_var) == TRUE)
  {
    prior_beta_var = diag(beta_dim)
  }
  
  #check!
  if(length(prior_beta_mean) != beta_dim)
  {stop("length of prior mean for beta should be of sum of no of total covariates for each variable")
  }
  
  #check
  
  if(is.positive.definite(sig) == FALSE)
  {stop("Covariance matrix of error should be symmetric positive definite matrix")
  }
  
  #check
  if(is.positive.definite(prior_beta_var) == FALSE)
  {stop("Prior covariance matrix of beta should be symmetric positive definite matrix")
  }
  
  
  
  #calculation of indicator to track starting and ending point of categories of each level to  create x matrix
  
  q = rep(1, length(category) + 1)  ## the indicator 
  q[1] = 1
  category_cum_sum = cumsum(category)
  category_cum_sum
  for(g in 2: length(q))
  {
    q[g] = category_cum_sum[g-1] - (g-1 - 1)    ## ex: p1 + p2 + p3 - 2
  }
  
  q
  
  
  ## For each subject,  nrow = z_dim = sum(category) - length(category), ncol = sum( (category - 1) * covariate_num)
  ## Xi s are stacked columnwise
  
  x = matrix(rep(0, z_dim * n * beta_dim), nrow = z_dim * n, ncol = beta_dim)
  dim(x)
  
  
  ##  calculation for x
  
  s = 0
  
  for(i in 1 : n)
  {
    ## only for 1st variable in zi
    
    x[ (s + 1) : (s + category[1] - 1) ,  1 : ((category[1]-1) * covariate_num[1]) ] = kronecker(diag(category[1] - 1), t(x_list[[1]][, i]) )
    
    ## for 2nd variable onwards
    s = s + category[1] - 1
    #print(s)
    
    for(j in 2:length(category) )
    {
      
      x[ (s + 1) : (s + category[j] - 1) , ((category[j-1]-1) * covariate_num[j-1]) + 1 : ((category[j]-1) * covariate_num[j]) ] = kronecker(diag(category[j] - 1), t(x_list[[j]][, i]) )
      s = s + category[j] - 1
      #print(s)
    }
  }
  
  
  
  
  
  # x[1:5,]
  # x[1:10,]
  # dim(x)
  # x1_mat[,2]
  # x2_mat[, 2]
  
  # likelihood of Z
  z_given_beta_var = kronecker(diag(n), sig)  #'Sigma' in the paper
  dim(z_given_beta_var)
  
  ## posteriors
  
  z_mean = prior_beta_mean   ## marginal 
  z_var = z_given_beta_var + x %*% prior_beta_var %*% t(x)  ## variance of marginal z distn
  dim(z_var)
  precision = solve(z_var)
  
  ## dim of f matrix
  ## for each subject, there would be  z_dim x z_dim matrix of constraints in f
  row_no = col_no = z_dim * n
  
  f = matrix(rep(0, row_no * col_no ), nrow = row_no, ncol = col_no )
  dim(f)
  s = 0 ## to run in the loop 
  
  
  for(i in 1:n)
  {
    
    for(j in 1:length(category))
    {
      
      if(d[j, i] == 0)
      {
        for(l in 1: (category[j]-1))
        {
          f[ s + l,  s + l]  = -1
        }
      }
      
      if(d[j, i] != 0)
      {
        a = matrix(rep(0, (category[j] -1) * (category[j] -1)), nrow = category[j] -1)  # when d[i] != 0
        for(l in 1 : (category[j]-1))
        {
          a[l, l] = ifelse(d[j,i] == l, 1, -1)
          f[ s + l,  s + l]  =  ifelse(d[j, i] == l, 1, -1)
        }
        w = which(diag(a) == 1)
        f[ c( (s + 1) : (s + (category[j]-1))), s + w] = rep(1, category[j]-1)
      }
      
      s = s + category[j] - 1
      #print(s)
    }
    
    
  }
  
  f[1:10, 1:10]
  d[,1:2]
  
  #d
  #f
  
  
  g = rep(0.001, row_no)
  r = rep(0, row_no)
  z = matrix(0, nrow = iter, ncol = row_no) # initialization
  dim(z)
  
  #HH posteriors of z
  for(i in 2:iter)
  {
    z[i, ] = tmg::rtmg(n = 1, M = precision, f = f, g = g, r = r, initial = z[(i-1), ])
  }  
  
  dim(z)
  
  Q = t(x) %*% solve(z_given_beta_var) %*% x  +  solve(prior_beta_var)  ## variance of marginal z|y(or d)
  dim(Q)
  L = chol(Q)
  
  
  beta = matrix(0, nrow = iter, ncol = beta_dim)
  dim(beta)
  
  for(i in 1 : iter)
  {
    b = t(x) %*% z[i,]
    z1 = rmvnorm(1, rep(0, beta_dim), diag(beta_dim) )
    y = solve(L) %*% t(z1)
    v = solve(t(L)) %*% b
    theta = solve(L) %*% v
    beta[i, ] = y + theta
  }
  
  ## Naming of coefficients
  names_category = category - 1
  pnames = rep(0, beta_dim)
  s = 0 # initialization
  for(i in 1: variable_dim)
  {
    for(j in 1: names_category[i])
    {
      
      pnames[(s + 1) : (s +  covariate_num[i])] = c(paste("beta[",i ,",",j,",", 1:covariate_num[i], "]", sep=""))
      s = s +  covariate_num[i]
      # print(s)
      
    }
  }
  
  
  pnames
  
  
  #Posterior mean 
  Betaout = beta[-c(1:burn), ]
  colnames(Betaout) = pnames
  postmean_beta = colMeans(Betaout)
  #rbind(postmean_beta,beta_act) # comparison with actual given betas
  
  # 95% Credible indervals
  alpha = 1- cred_level
  interval = apply(Betaout, 2, function(x) quantile(x, c((alpha/2), 1-(alpha/2))) )
  
  #effective sample size
  #sample_size = ess(Betaout)
  
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
  
  return(list(Posterior_mean = postmean_beta , Credible_interval = interval , trace_plot = trace, density_plot = density, carter_plot = carter))
  
} # end of function nominal_post_beta
