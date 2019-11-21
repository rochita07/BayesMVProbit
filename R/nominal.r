## Univariate Nomianal Measure, Sigma matrix is known 

n = 60 # no of subjects
p = 3 # no of levels
k = 6 # no of covariates
iter = 5000 # no of iterartion
burn = 1000 # no of burn in

## generation for beta of order p-1 x k

set.seed(1287)
beta_act = mvrnorm(k , mu = rep(2, p-1) , Sigma = diag(p-1))
beta_act = t(beta_act)  # each column is for each covariates, each row is for each level
dim(beta_act)

## generation of covariate matrix of order k x n

a_mat = mvrnorm(n , mu = rep(0, k) , Sigma = diag(k))
a_mat = t(a_mat) # each column is for each subject, each row is for each  covariates
dim(a_mat)


## generation of error matrix of order p-1 x n

e_mat = mvrnorm(n , mu = rep(0, p-1) , Sigma = diag(p-1))
e_mat = t(e_mat)
dim(e_mat)

## generation of error matrix of order p x n

beta_act_prod_a_mat = beta_act %*% a_mat
u_mat = beta_act_prod_a_mat  +  e_mat
dim(u_mat)

## generation of  vector d of dim n x 1

d = rep(0, n)

for(i in 1:n)
{
  x = u_mat[, i]
  
  if(max(x) < 0)
  {
    d[i] = 0
  }
  if(max(x) >= 0)
  {
    d[i] = which.max(x)
  }
}
d
table(d)
length(d)

## generation of matrix y of order p x n

y = matrix(rep(0, p * n), nrow = p, ncol = n)
for(i in 1:n){
  if(d[i] == 0)
  {
    y[p, i] = 1
  }
  if(d[i] != 0)
  {
    y[d[i] , i] = 1
  }
}
y
d

lam = 2  ## parameter in prior distn of beta

x = matrix(rep(0, (p-1) * n * (p-1) * k ), nrow = (p-1) * n, ncol = (p-1) * k)  ## to store the kronecker product

for(i in 1 : n)
{
  x[ c( ((p-1) * (i-1) + 1) : ((p-1) * (i-1) + (p-1)) ) , ] = kronecker(diag(p-1), t(a_mat[, i] ))
}

dim(x)

precision = solve( diag(rep(1, (p-1) * n )) + (x %*% t(x)) * lam ) #var-cov matrix
dim(precision)

## for each subject, there would be  p-1 x p-1 matrix of constraints in f
row_no = col_no = (p - 1) * n

f = matrix(rep(0, row_no * col_no ), nrow = row_no, ncol = col_no)  ## matrix of constrains
dim(f)

a = matrix(rep(0, (p-1) * (p-1)), nrow = p-1)  # when d[i] != 0

for(i in 1:n)
{
  
  if(d[i] == 0)
  {
    for(l in 1: (p-1))
    {
      f[ (p-1) * (i-1) + l,  (p -1)* (i-1) + l]  = -1
    }
  }
  
  if(d[i] != 0)
  { c = 0
  for(l in 1 : (p-1))
  {
    
    a[l, l] = ifelse(d[i] == l, 1, -1)
    f[ (p-1) * (i-1) + l,  (p -1)* (i-1) + l]  =  ifelse(d[i] == l, 1, -1)
  }
  
  w = which(diag(a) == 1)
  f[ c( ((p -1)* (i-1) + 1) : ((p -1)* (i-1) + (p-1))), (p -1)* (i-1) + w] = rep(1, p-1)
  }
  
}
