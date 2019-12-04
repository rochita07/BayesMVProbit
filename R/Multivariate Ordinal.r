
ordinal_posterior_beta = function()
{
  
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


for (k in 1: iter) # for iteration
{
  
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
    
  } # end of for(l in 1: q) loop
  
  
  
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
  
} 

return(beta_mat)
}## end of function ordinal_posterior_beta

beta_mat[1:20,]
beta_mat[iter,]

## Naming of coefficients

pnames = rep(0, beta_dim)
s = 0 # initialization
for(l in 1: length(covariate_num))
{
  pnames[(s+1) : (s+covariate_num[l])] = c(paste("beta[",l ,",", 1:covariate_num[l], "]", sep=""))
  s = s + covariate_num[l]
}
pnames

colnames(Betaout) = pnames

## calculation of posterior mean 

Betaout = beta_mat[-c(1:burn), ]
postmean_beta = colMeans(Betaout)
rbind(postmean_beta,beta_act) # comparison with actual given betas

# plotting of trace plots

par(mfrow = c(3,4))

for(i in 1 : ncol(Betaout))
{
  traceplot(as.mcmc(Betaout[,i]))
}

#effective sample size
ess(Betaout)  
