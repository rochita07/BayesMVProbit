## Univariate Ordinal
# Data generation 

cat = 4 ## no of categories
nu_tot = cat + 1 ## no of cut points 
nu = nu_tot - 3 ## no of cutpoints to be considered ( nu0, nu1, nuJ are fixed)
nu_cutoff = c(0, 1.5, 2.5)  # of length nu
df_t = 5 ## df of proposal t distn for nu
k = 4
n = 50
beta_act = rep(2, k)
set.seed(1234)
x = matrix(rnorm(k* n), nrow = k)
z = beta_act %*% x 
cutoff = c(-Inf, nu_cutoff, Inf)
y = rep(0, n) ## initialization
iter = 5000 ## no of iteration
burn = 1000 ## no of burn in

for(i in 1:n)
{
  for(j in 1: length(cutoff))
  {
    if(z[i] > cutoff[j] && z[i] <= cutoff[j + 1]) # making one side inclusive
      y[i] = j
  }
}
y
z
table(y)

## Prior and likelihood specifications 
sig = 1 ## error variance
prior_delta_mean = rep(0, nu)
prior_delta_var = diag(nu)
prior_beta_mean = rep(0, k)
prior_beta_var = diag(k)

## initializations required for the iteration loop

# beta = beta_act ## if not given to draw sample from prior !!!! Check

beta = mvtnorm::rmvnorm(n = 1, mean = prior_beta_mean , sigma = prior_beta_var)
beta = as.vector(beta)
delta_current = mvtnorm::rmvnorm(n = 1, mean = prior_delta_mean , sigma = prior_delta_var )
beta_update_var = solve(prior_beta_var + tcrossprod(x))
beta_mat = matrix(rep(0, k * iter), nrow = iter)

for (l in 1: iter)
{
  
## updates of delta

f = 1 ## initialization

y_given_beta_delta_optim = function(delta)
{
  
  prior_delta = mvtnorm::dmvnorm(delta, prior_delta_mean, prior_delta_var )  # prior density for delta
  
  # calculation of likelihood                        
  for(i in 1: n)
  {
    if(y[i] == 1)
    {
      a = (-t(x[,i]) %*% beta) / sig
      b = -Inf
      f = f * ( pnorm(a) - pnorm(b))   # pnorm(-inf) = 0 and nu1 = 0
    }
    
    if(y[i] == 2)
    {
      a = ( exp(delta[1]) - (t(x[,i]) %*% beta) )/ sig
      b = (-t(x[,i]) %*% beta) / sig
      f = f * ( pnorm( a )  - pnorm(b) )
    }
    
    for( j in 3 : (cat-1))
    {
      if(y[i] == j)
      {
        a = (sum(exp(delta[1:(j-1)])) -  (t(x[,i]) %*% beta) ) / sig
        b = (sum(exp(delta[1:(j-2)])) -  (t(x[,i]) %*% beta) ) / sig
        f = f * ( pnorm(a)  - pnorm(b) )
      }
    }
    if(y[i] == cat)
    {
      a = Inf
      b = (sum(exp(delta)) - (t(x[,i]) %*% beta) ) / sig
      f = f * ( pnorm( a )  - pnorm(b) )
    }
    
  }
  optim_fn = - f * prior_delta
  return(optim_fn)
  
}

delta_hat = optim(par = rep(0, nu), fn =  y_given_beta_delta_optim)$par
delta_hat

f = 1
y_given_beta_delta_hessian = function(delta)
{
  prior_delta = mvtnorm::dmvnorm(delta, prior_delta_mean, prior_delta_var )
  
  # calculation of likelihood                        
  for(i in 1: n)
  {
    if(y[i] == 1)
    {
      a = (-t(x[,i]) %*% beta) / sig
      b = -Inf
      f = f * ( pnorm(a) - pnorm(b))   # pnorm(-inf) = 0 and nu1 = 0
    }
    
    if(y[i] == 2)
    {
      a = ( exp(delta[1]) - (t(x[,i]) %*% beta) )/ sig
      b = (-t(x[,i]) %*% beta) / sig
      f = f * ( pnorm( a )  - pnorm(b) )
    }
    
    for( j in 3 : (cat-1))
    {
      if(y[i] == j)
      {
        a = (sum(exp(delta[1:(j-1)])) -  (t(x[,i]) %*% beta) ) / sig
        b = (sum(exp(delta[1:(j-2)])) -  (t(x[,i]) %*% beta) ) / sig
        f = f * ( pnorm(a)  - pnorm(b) )
      }
    }
    if(y[i] == cat)
    {
      a = Inf
      b = (sum(exp(delta)) - (t(x[,i]) %*% beta) ) / sig
      f = f * ( pnorm( a )  - pnorm(b) )
    }
    
  }
  
  hessian_fn = log(f * prior_delta)
  return(hessian_fn)
}
hessian_mat = numDeriv::hessian(func = y_given_beta_delta_hessian , x = delta_hat)
D_hat = solve(-hessian_mat)
D_hat  # var- cov matrix , to be converted into scale matrix

scale_t = D_hat * ((df_t-2) / df_t)  # scale matrix for poposal t distn

# delta_current = mvtnorm::rmvt(n = 1, delta = delta_hat , sigma = scale_t, df = df_t)  
delta_proposed = mvtnorm::rmvt(n = 1, delta = delta_hat , sigma = scale_t, df = df_t)
density_delta_current = mvtnorm::dmvt(delta_current, delta = delta_hat , sigma = scale_t, df = df_t, log = FALSE )
density_delta_proposed = mvtnorm::dmvt(delta_proposed, delta = delta_hat , sigma = scale_t, df = df_t, log = FALSE )
prob_MH_1st_part = y_given_beta_delta_optim(delta_proposed) / y_given_beta_delta_optim(delta_current)
prob_MH_2nd_part = density_delta_current / density_delta_proposed
prob_MH = min(1, prob_MH_1st_part * prob_MH_2nd_part)

u = runif(1,0,1)  # to draw an independent sample 

if(u <= prob_MH)
{
  delta_current = delta_proposed
}
if(u > prob_MH)
{
  delta_current = delta_current
}

## Updates of z
cutoff_update = c(-Inf, 0, cumsum(exp(delta_current)), Inf)  ## update on nu 
z_update = rep(0, n) ## initialization
# lower_z = rep(0, n) ## initialization
# upper_z = rep(0, n) ## initialization


for(i in 1: n)
{
  for(j in 1: cat)
  {
    if(y[i] == j) # j = 1, 2, 3, ..., cat
    {
      lower_z = cutoff_update[j] 
      upper_z = cutoff_update[j + 1] 
      z_update[i] = truncnorm::rtruncnorm(n = 1, a = lower_z, b = upper_z, mean = t(x[,i]) %*% beta, sd = sig)
      
    }
  }
}

## Updates of beta

# beta_update_var = solve(prior_beta_var + tcrossprod(x))
beta_update_mean = beta_update_var %*% (solve(prior_beta_var) %*% prior_beta_mean  + x %*% z_update)
beta_update = mvtnorm::rmvnorm(n = 1, mean = beta_update_mean, sigma = beta_update_var )

## initialization for next iteration 
beta = as.vector(beta_update)
delta = delta_current
beta_mat[l, ] = beta
}
beta_mat[iter,]

# posterior mean of beta
Betaout = beta_mat[-c(1:burn), ]
postmean_beta = colMeans(Betaout)
rbind(postmean_beta,beta_act) # comparison with actual given betas