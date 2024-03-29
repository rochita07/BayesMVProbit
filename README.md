# Package Name:
“BayesMVProbit”


# Package Title:
Data augmentation Gibbs sampling for Multivariate Multinomial Probit model and Multivariate ordered Probit model


# Intended functionality / short description:

This package will implement the marginal data augmentation algorithms for posterior sampling in the Multivariate Multinomial Probit and Multivariate ordered Probit model. We consider the analysis of models for univariate and multivariate outcomes in the context of the latent variable inferential framework of Albert and Chib (1993) and Holmes and Held (2006). This will also provide trace plots, density plot and carter plot of the parameters. The theories considered to build this package are mentioned in the reference section and also discussed in vignette in details.


## 1. Multivariate Multinomial Probit model: 

Let there are g nominal measures. For say qth measure we have k categories and we consider k - 1 latent variables. As per the method consider the maximum latent variable and consider the corresponding y value to be 1 or 0 otherwise. Here the latent variable follows multivariate normal with appropriate dimension. Hence the model is called probit. The detail method is discussed in vignette.

The method of generating data is as follow: 

```{r}

variable_dim = 2 # no of variables
category = c(3,4) # no of levels for each variable
n = 50 # no of subjects
covariate_num = c(2,3) # no of covariates for each level

z_dim = sum(category) - length(category)  ## dimension of z for each subject
beta_dim = sum((category - 1) * covariate_num)

set.seed(1287)

beta_act1 = matrix(rep(2, (category[1]-1) * covariate_num[1]), nrow = category[1]-1)
  ## matrix of regression coefficient for 1st nominal measure
beta_act2 = matrix(rep(3, (category[2]-1) * covariate_num[2]), nrow = category[2]-1)
 ## matrix of regression coefficient for 1st nominal measure
beta_act = as.vector(c(beta_act1, beta_act2)) ## vectorized form of beta
## should be of same length  as beta_dim

x1_mat = matrix(rnorm(n * covariate_num[1]), nrow = covariate_num[1] ) # each column is for each subject, each row is for each  covariates
x2_mat = matrix(rnorm(n * covariate_num[2]), nrow = covariate_num[2] ) # each column is for each subject, each row is for each  covariates
x_list = list(x1_mat, x2_mat)  # input should be given as list 

z1_mat = beta_act1 %*% x1_mat 
z2_mat = beta_act2 %*% x2_mat


sig1 = matrix(c(1, 0.36, 0.36, 0.81), nrow = 2)
sig2 = matrix(c(1, 0.45, 0.24, 0.45, 0.81, 0.41, 0.24, 0.41, 0.9), nrow =3)
sig3 = matrix(c(rep(0.2, 6)), nrow =2)
sig4 = matrix(c(rep(0.2, 6)), nrow =3)
sig_gen = as.matrix(rbind(cbind(sig1, sig3), cbind(sig4, sig2))) ## Final sigma 

# GEneration of error 

e_mat = MASS::mvrnorm(n , mu = rep(0, z_dim) , Sigma = sig_gen )  # Generation of error term
e_mat = t(e_mat) # each column is for each subject

#Generation of  Z vectors for each subject, each column is for each subject
z_mat = matrix(rep(0, n * z_dim), ncol = n) ## each column is for each subject

for(i in 1: n){
  z_mat[, i] = c(z1_mat[, i], z2_mat[, i]) + e_mat[, i]
  
}


#Computation of d matrix ## user should provide this input

d = matrix(rep(0, n * length(category)), ncol = n)  # each col corresponds to each subject
q = rep(1, length(category) + 1)  ## the indicator to compute d matrix (as mentioned in the paper)
q[1] = 1
category_cum_sum = cumsum(category)

for(g in 2: length(q))
{
  q[g] = category_cum_sum[g-1] - (g-1 - 1)    ## ex: p1 + p2 + p3 - 2
}


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

d #Example of d matrix
```

For better mixing, we consider Holmes Held approach, where we sample latent variable after marginalizing beta.

We sample beta from its posterior distribution.

We calculate posterior mean and credible interval for beta.

We also provide trcaeplot, density plot and carter plot for posterior sample of beta.



## 2. Multivariate Ordinal Probit model: 

Let there are g nominal measures. For each measure we have one latent variable. Say, qth measure has k categories, as per the method k+1 cutoff points are considered. The method considers that if the latent variable falls into certain range of cutoff pints, corresponding y takes a certian value. Here the latent variable follows multivariate normal with appropriate dimension. Hence the model is called probit. The detail method is discussed in vignette.

For simplicity, delta, a one to one trandformation of cutoff points  are considered for each level (as discussed in the paper).

The method of generating data is as follow: 


```{r}

variable_dim = 2 # no of variables, i.e dimension of yi (can be obtained if y is given)
category = c(4, 3) ## no of categories(levels) for each variable
covariate_num = c(2, 3) # No of covariates for each level
nu_tot = category + 1 ## no of cut points 
nu = nu_tot - 3 ## no of cutpoints to be considered for each variable ( nu0, nu1, nuJ are fixed)
nu_cutoff1 = c(0, 2 , 2.5)  # of length nu + 1
nu_cutoff2 = c(0, 4.5)
beta_dim = sum(covariate_num) # dimension of beta vector

n = 50 # no of subjects

## GEneration of beta

beta_act1 = rep(2, covariate_num[1])
beta_act2 = rep(2,covariate_num[2])
beta_act = c(beta_act1, beta_act2)

set.seed(1287)

#GEneration of X

x1 = matrix(rnorm(covariate_num[1]* n), nrow = covariate_num[1])  # each column is for each subject
x2 = matrix(rnorm(covariate_num[2]* n), nrow = covariate_num[2])  # each column is for each subject

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


# Generation of error 

sig = diag(variable_dim) ## error covariance matrix of order variable_dim x variable_dim ( to be mentioned in the input)

e = mvtnorm::rmvnorm(n = n, mean = rep(0, variable_dim) , sigma = sig)
e = t(e) # # each column is for each subject


# Generation of  Z
z = matrix(rep(0, variable_dim * n), nrow = variable_dim)  # matrix of order variable_dim x n # # each column is for each subject
for(i in 1:n)
{
  z[,i] =  x[[i]] %*% beta_act + e[, i]
}


## Generation of Y 

cutoff1 = c(-Inf, nu_cutoff1, Inf)  ## To also consider the fixed cutoffs 
cutoff2 = c(-Inf, nu_cutoff2, Inf)
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


y  ## example of y

```

We sample for delta from a proposed t distribution (parameters are optimized as mentioned in the paper) 

We sample for letent variable from truncated normal distribution , truncation regions are specific cutoff points corresponding the y value

We sample beta from its posterior distribution

We calculate posterior mean and credible interval for beta.

We also provide trcaeplot, density plot and carter plot for posterior sample of beta.


# Intended users: 

Mainly students working in Bayesian field with multivariate data. Specifically, this kind of models are mainly used in transportation system, Econometrics like labor economics (educational attainment), political economy (voter opinions), and health economics (consumers’ reliance on alternative sources of medical information) and so on. 





# References:

1.	XiaoZhang, W. JohnBoscardin, Thomas R.Belin (2008), “Bayesian Analysis of Multivariate Nominal measures using Multivariate Multinomial Probit Models”, Computational Statistics & Data Analysis, Volume 52, Issue 7, 15 March 2008, Pages 3697-3708

[Bayesian Analysis of Multivariate Nominal Measures Using Multivariate Multinomial Probit Models; Xiao Zhang, W. John Boscardin, and Thomas R. Belin](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2673013/)

2.	Ivan Jeliazkov, Jennifer Graves, Mark Kutzbach (2008), “Fitting And Comparison of Models For Multivariate Ordinal Outcomes”, Advances in Econometrics: Bayesian Econometrics, 23

[FITTING AND COMPARISON OF MODELS FOR MULTIVARIATE ORDINAL OUTCOMES; Ivan Jeliazkov, Jennifer Graves, Mark Kutzbach ](https://escholarship.org/content/qt3xz66103/qt3xz66103.pdf)
