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
