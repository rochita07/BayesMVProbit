category = c(3,4) # no of levels for each variable
n = 50 # no of subjects
covariate_num = c(2,3) # no of covariates to be calculated from given x matrix
iter = 20 # no of iterartion
burn = 10 # no of burn in
z_dim = sum(category) - length(category)  ## dimension of z for each subject
beta_dim = sum((category - 1) * covariate_num)
variable_dim = 2 # no of variables

set.seed(1287)
beta_act1 = matrix(rep(2, (category[1]-1) * covariate_num[1]), nrow = category[1]-1)
beta_act1
#beta_act1 = t(beta_act1)
beta_act2 = matrix(rep(3, (category[2]-1) * covariate_num[2]), nrow = category[2]-1)
beta_act2
#beta_act2 = t(beta_act1)

beta_act = as.vector(c(beta_act1, beta_act2))
beta_act
length(beta_act)  ## should be same as beta_dim

x1_mat = matrix(rnorm(n * covariate_num[1]), nrow = covariate_num[1] ) # each column is for each subject, each row is for each  covariates
#x1_mat = t(x1_mat)
dim(x1_mat)

x2_mat = matrix(rnorm(n * covariate_num[2]), nrow = covariate_num[2] ) # each column is for each subject, each row is for each  covariates
#x2_mat = t(x2_mat)
dim(x2_mat)

x_list = list(x1_mat, x2_mat)
x_list

z1_mat = beta_act1 %*% x1_mat
z2_mat = beta_act2 %*% x2_mat
z1_mat
z2_mat

#### Sigma is calculated according the example in the paper. (but can't find to apply a correlation matrix)

sig1 = matrix(c(1, 0.36, 0.36, 0.81), nrow = 2)
sig2 = matrix(c(1, 0.45, 0.24, 0.45, 0.81, 0.41, 0.24, 0.41, 0.9), nrow =3)
sig3 = matrix(c(rep(0.2, 6)), nrow =2)
sig4 = matrix(c(rep(0.2, 6)), nrow =3)
sig_gen = as.matrix(rbind(cbind(sig1, sig3), cbind(sig4, sig2)))
sig_gen  ## Final sigma matrix for e

#sig_gen = diag(z_dim)
e_mat = mvrnorm(n , mu = rep(0, z_dim) , Sigma = sig_gen )  # specific to the paper
e_mat = t(e_mat)
dim(e_mat)

z_mat = matrix(rep(0, n * z_dim), ncol = n)
dim(z_mat)



for(i in 1: n){
  
  z_mat[, i] = c(z1_mat[, i], z2_mat[, i]) + e_mat[, i]
  
}

z_mat


d = matrix(rep(0, n * length(category)), ncol = n)  # each col corresponds to each subject

q = rep(1, length(category) + 1)  ## the indicator 
q[1] = 1
category_cum_sum = cumsum(category)
category_cum_sum
for(g in 2: length(q))
{
  q[g] = category_cum_sum[g-1] - (g-1 - 1)    ## ex: p1 + p2 + p3 - 2
}

q


for(i in 1: n)
{ 
  for(j in 1:variable_dim )
  {
    x = z_mat[ q[j] : (q[j+1] - 1) , i]
    # print(x)
    
    if(max(x) < 0)
    {
      d[j, i] = 0
    }
    if(max(x) >= 0)
    {
      d[j, i] = which.max(x)
    }
  }
}
d
dim(d)
z_mat[, 6:10]
#d[,2]
d[, 6:10]
table(d[1,])
table(d[2,])


## Prior

prior_beta_mean = rep(0, beta_dim)  ## !!! can be changed   ##'mu' in paper
prior_beta_var = diag(beta_dim) ## of order z_dim x z_dim  ## 'C' in  paper ## !!!!! to be checked
dim(prior_beta_var)

# z_given_beta_var[1:15, 1:15]

category = c(3,4)
iter = 10
burn = 2
cred_level = 0.95
x_list
sig = sig_gen
d
prior_beta_mean = NULL
prior_beta_var = NULL




nominal_post_beta(category = c(3,4), iter = 10, burn = 2, cred_level = 0.95, x_list, sig, d,
                  prior_beta_mean = NULL, prior_beta_var = NULL)
