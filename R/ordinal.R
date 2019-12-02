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