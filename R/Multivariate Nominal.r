## Multivariate Nomianal Measure, Sigma matrix is known, and xi's are subject specific
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
z_mat = beta_act %*% x_mat
dim(x_mat)
dim(z_mat)

# A specific choice of sigma matrix (mentioned in the example of paper)
sig1 = matrix(c(1, 0.36, 0.36, 0.81), nrow = 2)
sig2 = matrix(c(1, 0.45, 0.24, 0.45, 0.81, 0.41, 0.24, 0.41, 0.9), nrow =3)
sig3 = matrix(c(rep(0.2, 6)), nrow =2)
sig4 = matrix(c(rep(0.2, 6)), nrow =3)
sig_gen = as.matrix(rbind(cbind(sig1, sig3), cbind(sig4, sig2)))
sig_gen  ## Final sigma matrix for e
