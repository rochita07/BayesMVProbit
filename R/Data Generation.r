## Ordinal DAta generation 

variable_dim = 2 # no of variables, i.e dimension of yi
category = c(4, 3) ## no of categories(levels) for each variable

nu_tot = category + 1 ## no of cut points 
nu = nu_tot - 3 ## no of cutpoints to be considered for each variable ( nu0, nu1, nuJ are fixed)
nu
nu_cutoff1 = c(0, 2 , 2.5)  # of length nu + 1
nu_cutoff2 = c(0, 4.5)
#nu_cutoff3 = c(0, 0.5, 1, 1.5, 2)

df_t = rep(5 , variable_dim)  ## df of proposal t distn for nu  # otherwise to be specified
covariate_num = c(2, 3)  ## no of covariates for each variable

n = 50 # no of subjects
iter = 10  # no of iteration
burn = 1 # no of burn in 
cred_level = 0.95 

## to generate the data
beta_act1 = rep(2, covariate_num[1])
beta_act2 = rep(2,covariate_num[2])
#beta_act3 = rep(1, covariate_num[3])
beta_act = c(beta_act1, beta_act2)
beta_dim = sum(covariate_num)

x1 = matrix(rnorm(covariate_num[1]* n), nrow = covariate_num[1])  # each column is for each subject
x1 = matrix(rnorm(n), ncol = n)
x2 = matrix(rnorm(covariate_num[2]* n), nrow = covariate_num[2])  # each column is for each subject
#x3 = matrix(rnorm(covariate_num[3]* n), nrow = covariate_num[3])
x_list = list(x1, x2)

x = lapply(1:n, function(x) matrix(rep(0, variable_dim * beta_dim), nrow = variable_dim) ) ## initialization
a = matrix(rep(0, variable_dim * beta_dim), nrow = variable_dim) ## to store the value in the loop
col_indi = c(0, cumsum(covariate_num))  # of length variable_dim + 1
for(i in 1:n)
{
  for(l in 1: variable_dim)
  {
    a[l, (col_indi[l] + 1) : col_indi[l + 1]] = t(x_list[[l]][,i])
  }
  x[[i]] = a
}
x ##


# generation of error 
sig = diag(variable_dim) ## error covariance matrix of order variable_dim x variable_dim ( to be mentioned in the input)

e = rmvnorm(n = n, mean = rep(0, variable_dim) , sigma = sig)
e = t(e) # # each column is for each subject
dim(e)

# Generation of  Z
z = matrix(rep(0, variable_dim * n), nrow = variable_dim)  # matrix of order variable_dim x n # # each column is for each subject
for(i in 1:n)
{
  z[,i] =  x[[i]] %*% beta_act + e[, i]
}
z

## Generation of Y 

cutoff1 = c(-Inf, nu_cutoff1, Inf)
cutoff2 = c(-Inf, nu_cutoff2, Inf)
#cutoff3 = c(-Inf, nu_cutoff3, Inf)
cutoff = list(cutoff1, cutoff2)

y = matrix(rep(0, variable_dim * n), nrow = variable_dim) ## initialization of matrix of order variable_dim x n # each column is for each subject


for(i in 1:n)
{
  for(l in 1:variable_dim)
  {
    for(j in 1: length(cutoff[[l]]))
    {
      if(z[l, i] > cutoff[[l]][j] && z[l, i] <= cutoff[[l]][j + 1]) # making one side inclusive
        y[l, i] = j
    }
  }
}
y_original = y # to be given as input
y_original
z
table(y_original[1,])
table(y_original[2,])
#table(y[3,])
# cutoff3
# cbind(y[3,], z[3,])
# max(z[3,])



## Prior and likelihood specifications ## null values otherwise specified

prior_delta_mean1 = rep(0, nu[1])
#prior_delta_mean1 = rep(0, 1)
prior_delta_mean2 = rep(0, nu[2])
#prior_delta_mean3 = rep(0, nu[3])
prior_delta_mean = list(prior_delta_mean1, prior_delta_mean2)

prior_delta_var1 = diag(nu[1])
prior_delta_var2 = diag(nu[2])
#prior_delta_var3 = diag(nu[3])
prior_delta_var = list(prior_delta_var1, prior_delta_var2)

prior_beta_mean = rep(0, beta_dim)
prior_beta_var = diag(beta_dim)


category = c(4, 3)
df_t = NULL
iter = 10
burn = 1
cred_level = 0.95
sig = diag(2)
y = y_original
prior_delta_mean = NULL
prior_delta_var = NULL
prior_beta_mean = NULL
prior_beta_var = NULL







ordinal_post_beta(category = c(4, 3), df_t = NULL, iter = 10,  burn = 1, cred_level = 0.95, x_list, sig = diag(2), y = y_original,
                  prior_delta_mean = NULL, prior_delta_var = NULL, prior_beta_mean = NULL, prior_beta_var = NULL)

