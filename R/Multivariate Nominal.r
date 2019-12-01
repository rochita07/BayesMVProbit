## Multivariate Nomianal Measure, Sigma matrix is known, Holmes Held method, and xi's are subject specific
#(there may be the case of outcome- subject specific)

library(tmg)
##Generation of  data 
n = 50 # no of subjects
p = c(3,4) # no of levels
k = 4 # no of covariates
iter = 5000 # no of iterartion
burn = 1000 # no of burn in
z_dim = sum(p) - length(p)  # dimension of z for each subject
beta_dim = z_dim * k  # order of beta matrix  = length of beta_vec (after vectorization)
## beta is a matrix of order z_dim x k

#GEneration of beta matrix and vectorization
set.seed(1287)
beta_act = matrix(rep(2, beta_dim), nrow = z_dim)
#Vectorization of beta
beta_act_vec = as.vector(t(beta_act)) # Transpose is considered to vectorize row wise
length(beta_act_vec)

#Generation of X matrix
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


# Creation of X metrix 
## For each subject, Xi is a mtrix of  nrow = z_dim = sum(p) - length(p), ncol = z_dim * beta_dim
## Xi s are stacked columnwise

x = matrix(rep(0, z_dim * n * beta_dim), nrow = z_dim * n, ncol = beta_dim )
dim(x)

s = 0  # Initialization, to consider the counts of row

for(i in 1 : n)
{
  ## only for 1st variable in zi
  
  x[ (s + 1) : (s + z_dim) ,]  = kronecker(diag(z_dim), t(x_mat[, i]) )
  s = s + z_dim
}


x[1:5,]
x[1:10,]
dim(x)
x_mat[,1]
x_mat[, 2]


## Prior distribution of beta and distribution of z (as vector) specification:

prior_beta_mean = rep(0, beta_dim)  # prior mean of beta
prior_beta_var = diag(beta_dim)  #  prior variance of beta
dim(prior_beta_var)
z_given_beta_var = kronecker(diag(n), sig_gen)  #'Sigma' = variance of z (as vector of length z_dim *n)
dim(z_given_beta_var)

## posterior of Z

z_mean = prior_beta_mean   ## marginal 
z_var = z_given_beta_var + x %*% prior_beta_var %*% t(x)  ## marginal
dim(z_var)
precision = solve(z_var)

## Creating matrix of constraints 'f'
## for each subject, there would be  z_dim x z_dim matrix of constraints in f
row_no = col_no = z_dim * n

f = matrix(rep(0, row_no * col_no ), nrow = row_no, ncol = col_no )
dim(f)
s = 0   # Initialization, to consider the counts of row


for(i in 1:n)
{
  
  for(j in 1:length(p))
  {
    
    if(d[j, i] == 0)
    {
      for(l in 1: (p[j]-1))
      {
        f[ s + l,  s + l]  = -1
      }
    }
    
    if(d[j, i] != 0)
    {
      a = matrix(rep(0, (p[j] -1) * (p[j] -1)), nrow = p[j] -1)  # when d[i] != 0
      for(l in 1 : (p[j]-1))
      {
        a[l, l] = ifelse(d[j,i] == l, 1, -1)
        f[ s + l,  s + l]  =  ifelse(d[j, i] == l, 1, -1)
      }
      w = which(diag(a) == 1)
      f[ c( (s + 1) : (s + (p[j]-1))), s + w] = rep(1, p[j]-1)
    }
    
    s = s + p[j] - 1
    #print(s)
  }
  
  
}

f[1:10, 1:10]
d[,1:2]

g = rep(0.001, row_no) 
r = rep(0, row_no)
z = matrix(0, nrow = iter, ncol = row_no) # initialization
dim(z)

#HH posteriors of z

for(i in 2:iter)
{
  z[i, ] = rtmg(1, M = precision, f = f, g = g, r = r, initial = z[(i-1), ])
}  


