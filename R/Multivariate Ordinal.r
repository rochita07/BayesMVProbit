## Given inputs

q = 3 # no of variables, i.e dimension of yi
cat = c(4, 5, 6) ## no of categories for each variable
nu_tot = cat + 1 ## no of cut points 
nu = nu_tot - 3 ## no of cutpoints to be considered for each variable ( nu0, nu1, nuJ are fixed)
nu_cutoff1 = c(0, 0.5, 1)  # of length nu + 1
nu_cutoff2 = c(0, 1, 1.5 , 2)
nu_cutoff3 = c(0, 0.5, 1, 1.5, 2)
df_t = rep(5 , q)  ## df of proposal t distn for nu  # otherwise to be specified
k = c(2, 3, 2)  ## no of covariates for each variable
n = 50 # no of subjects
iter = 5000  # no of iteration
burn = 1000 # no of burn in 

## to generate the data

beta_act1 = rep(2, k[1])
beta_act2 = rep(3,k[2])
beta_act3 = rep(1, k[3])
beta_act = c(beta_act1, beta_act2, beta_act3)
beta_dim = sum(k)
x1 = matrix(rnorm(k[1]* n), nrow = k[1])  # each column is for each subject
x2 = matrix(rnorm(k[2]* n), nrow = k[2])  # each column is for each subject
x3 = matrix(rnorm(k[3]* n), nrow = k[3])
x_mat = list(x1, x2, x3)

#Generation of Xi matrix as list
x = lapply(1:n, function(x) matrix(rep(0, q * beta_dim), nrow = q) ) ## initialization
a = matrix(rep(0, q * beta_dim), nrow = q)
col_indi = c(0, cumsum(k))  # of length q+1  # index of starting and ending points for column in Xi matrix 
for(i in 1:n)
{
  for(l in 1: q)
  {
    a[l, (col_indi[l] + 1) : col_indi[l + 1]] = t(x_mat[[l]][,i])
  }
  x[[i]] = a
}
x

# generation of error 
sig = diag(q) ## error covariance matrix of order q x q ( to be mentioned in the input)

e = rmvnorm(n = n, mean = rep(0, q) , sigma = sig)
e = t(e) # # each column is for each subject

# Generation of  Z
z = matrix(rep(0, q * n), nrow = q)  # matrix of order q x n # # each column is for each subject
for(i in 1:n)
{
  z[,i] =  x[[i]] %*% beta_act + e[, i]
}
z

## Generation of Y 

cutoff1 = c(-Inf, nu_cutoff1, Inf)
cutoff2 = c(-Inf, nu_cutoff2, Inf)
cutoff3 = c(-Inf, nu_cutoff3, Inf)
cutoff = list(cutoff1, cutoff2, cutoff3)

y = matrix(rep(0, q * n), nrow = q) ## initialization of matrix of order q x n # each column is for each subject


for(i in 1:n)
{
  for(l in 1:q)
  {
    for(j in 1: length(cutoff[[l]]))
    {
      if(z[l, i] > cutoff[[l]][j] && z[l, i] <= cutoff[[l]][j + 1]) # making one side inclusive
        y[l, i] = j
    }
  }
}
y

## Prior and likelihood specifications ## null values otherwise specified

prior_delta_mean1 = rep(0, nu[1])
prior_delta_mean2 = rep(0, nu[2])
prior_delta_mean3 = rep(0, nu[3])
prior_delta_mean = list(prior_delta_mean1, prior_delta_mean2, prior_delta_mean3)

prior_delta_var1 = diag(nu[1])
prior_delta_var2 = diag(nu[2])
prior_delta_var3 = diag(nu[3])
prior_delta_var = list(prior_delta_var1, prior_delta_var2, prior_delta_var3)

prior_beta_mean = rep(0, beta_dim)
prior_beta_var = diag(beta_dim)

## initializations required for the iteration loop

beta = rmvnorm(n = 1, mean = prior_beta_mean , sigma = prior_beta_var)
beta = as.vector(beta)

delta_current_list = cutoff  ## just to initialize as list of q to store q delta vectors 

for(l in 1:q)
{
  delta_current_list[[l]] = rmvnorm(n = 1, mean = prior_delta_mean[[l]] , sigma = prior_delta_var[[l]] )
}
delta_current_list

beta_mat = matrix(rep(0, beta_dim * iter), nrow = iter)

cutoff_update = cutoff  # initialization as list
z_update = z ## initialization


# To calculate conditional mean and variance for z[l,i]
given.ind = 1: q


# Calculation of variance of posterior of beta

t_x_sig_x = lapply(x, function(x) t(x) %*% solve(sig) %*% x)
sum_t_x_sig_x = matrix(rep(0, beta_dim^2), nrow = beta_dim) ## initialization
for(i in 1: n)
{
  sum_t_x_sig_x = sum_t_x_sig_x + t_x_sig_x[[i]]
}
sum_t_x_sig_x 
beta_update_var = solve(solve(prior_beta_var) + sum_t_x_sig_x )  # variance of posterior of beta



## updates of delta

#   
for(l in 1:q)
{
  ## for a fixed l
  y_given_beta_delta_optim = function(delta)
  {
    
    prior_delta = dmvnorm(delta, prior_delta_mean[[l]], prior_delta_var[[l]] )  # prior density for delta
    
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
      
      for( j in 3 : (cat[l]-1))
      {
        if(y[l, i] == j )
        {
          a = (sum(exp(delta[1:(j-1)])) -  (x[[i]][l,] %*% beta) ) / diag(sig)[l]
          b = (sum(exp(delta[1:(j-2)])) -  (x[[i]][l,] %*% beta) ) / diag(sig)[l]
          f = f * ( pnorm(a)  - pnorm(b) )
        }
      }
      if(y[l, i] == cat[l])
      {
        #a = Inf
        b = (sum(exp(delta)) - (x[[i]][l,] %*% beta) ) / diag(sig)[l]
        f = f * ( 1  - pnorm(b) )
      }
      
    }
    optim_fn = - (f * prior_delta)  # to minimize in optim fn (-ve) is considered
    return(optim_fn)
    
  }
  
  ## updates of delta
  
  
  delta_hat = optim(par = rep(0, nu[l]), fn =  y_given_beta_delta_optim)$par
  # optimum value of delta to be used in proposed t distn as parameter
  
  y_given_beta_delta_hessian = function(delta)
  {
    prior_delta = dmvnorm(delta, prior_delta_mean[[l]], prior_delta_var[[l]] )  # prior density for delta
    
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
      
      for( j in 3 : (cat[l]-1))
      {
        if(y[l, i] == j )
        {
          a = (sum(exp(delta[1:(j-1)])) -  (x[[i]][l,] %*% beta) ) / diag(sig)[l]
          b = (sum(exp(delta[1:(j-2)])) -  (x[[i]][l,] %*% beta) ) / diag(sig)[l]
          f = f * ( pnorm(a)  - pnorm(b) )
        }
      }
      if(y[l, i] == cat[l])
      {
        #a = Inf
        b = (sum(exp(delta)) - (x[[i]][l,] %*% beta) ) / diag(sig)[l]
        f = f * ( 1  - pnorm(b) )
      }
      
    }
    
    hessian_fn = log(f * prior_delta)
    return(hessian_fn)
  }
  
  hessian_mat = hessian(func = y_given_beta_delta_hessian , x = delta_hat)
  D_hat = solve(-hessian_mat)
  D_hat  # var- cov matrix , to be converted into scale matrix
  
  scale_t =  D_hat * ((df_t[l]-2) / df_t[l])  # scale matrix for poposal t distn
  scale_t
  # delta_current = mvtnorm::rmvt(n = 1, delta = delta_hat , sigma = scale_t, df = df_t)  
  delta_proposed = rmvt(n = 1, delta = delta_hat , sigma = scale_t, df = df_t[l])
  density_delta_current = dmvt(delta_current_list[[l]], delta = delta_hat , sigma = scale_t, df = df_t[l], log = FALSE )
  density_delta_proposed = dmvt(delta_proposed, delta = delta_hat , sigma = scale_t, df = df_t[l], log = FALSE )
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
    for(j in 1: cat[l])
    {
      if(y[l, i] == j) # j = 1, 2, 3, ..., cat
      {
        lower_z = cutoff_update[[l]][j] 
        upper_z = cutoff_update[[l]][j + 1] 
        
        cond = condMVN(mean = x[[i]] %*% beta, sigma = sig, dep=l, given = given.ind[-l], X.given = rep(1, (q-1)), check.sigma=FALSE )
        
        z_update[l,i] = rtruncnorm(n = 1, a = lower_z, b = upper_z, mean = cond$condMean, sd = sqrt(cond$condVar))
        
      }
    } # for loop j : 1,..,cat[l]
  } # for loop i : 1,..,n
  
  
  ## Updates of beta
  
  # beta_update_var = solve(prior_beta_var + tcrossprod(x))
  
  sum_t_x_sig_z = matrix(rep(0, beta_dim), nrow = beta_dim) ## initialization
  for(i in 1: n)
  {
    sum_t_x_sig_z = sum_t_x_sig_z  + t(x[[i]]) %*% solve(sig) %*% z[,i] 
  }
  sum_t_x_sig_z 
  
  beta_update_mean = beta_update_var %*% (solve(prior_beta_var) %*% prior_beta_mean  + sum_t_x_sig_z)
  beta_update = rmvnorm(n = 1, mean = beta_update_mean, sigma = beta_update_var )
  
  ## initialization for next iteration 
  beta = as.vector(beta_update)
  #delta = delta_current_list
  beta_mat[k, ] = beta
  #print(beta)
  