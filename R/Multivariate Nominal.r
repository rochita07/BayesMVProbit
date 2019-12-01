## Multivariate Nomianal Measure, Sigma matrix is known, Holmes Held method, and xi's are subject specific
#(there may be the case of outcome- subject specific)

##Generation of  data 
n = 50 # no of subjects
p = c(3,4) # no of levels
k = 4 # no of covariates
iter = 5000 # no of iterartion
burn = 1000 # no of burn in
z_dim = sum(p) - length(p)  # dimension of z for each subject
beta_dim = z_dim * k  # length of beta
## beta is a matrix of order z_dim x k

set.seed(1287)
beta_act = matrix(rep(2, beta_dim), nrow = z_dim)
beta_act
x_mat = matrix(rnorm(n * k), nrow = k ) # each column is for each subject, each row is for each  covariates
dim(x_mat)

#Generation of error matrix

sig1 = matrix(c(1, 0.36, 0.36, 0.81), nrow = 2)
sig2 = matrix(c(1, 0.45, 0.24, 0.45, 0.81, 0.41, 0.24, 0.41, 0.9), nrow =3)
sig3 = matrix(c(rep(0.2, 6)), nrow =2)
sig4 = matrix(c(rep(0.2, 6)), nrow =3)
sig_gen = as.matrix(rbind(cbind(sig1, sig3), cbind(sig4, sig2)))
sig_gen  ## Final sigma matrix for error # A specific choice of sigma matrix (mentioned in the example of paper)

e_mat = MASS::mvrnorm(n , mu = rep(0, z_dim) , Sigma = sig_gen )  # generation of error 
e_mat = t(e_mat)   # each column is for each subject, each row is for each  covariates
dim(e_mat)

# Generation of z_mat (utility matrix)
z_mat = matrix(rep(0, n * z_dim), ncol = n)   # each column is for each subject, each row is for each  covariates
dim(z_mat)
z_mat = beta_act %*% x_mat + e_mat

#Creation of index of chosen alternatives
d = matrix(rep(0, n * length(p)), ncol = n)  # each col corresponds to each subject


q = rep(1, length(p) + 1)  ## the indicator  to specify the location of a particular variable
q[1] = 1  # by default starts from 1st variable at 1st position
p_cum_sum = cumsum(p)
p_cum_sum
for(g in 2: length(q))
{
  q[g] = p_cum_sum[g-1] - (g-1 - 1)    ## ex: p1 + p2 + p3 - 2
}

q

for(i in 1: n)
{ 
  for(j in 1:length(p) )
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

