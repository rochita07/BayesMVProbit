# Package Name:
“BayesMVProbit”


# Package Title:
Data augmentation Gibbs sampling for Multivariate Multinomial Probit model and Multivariate ordered Probit model


# Intended functionality / short description:

This package will implement the marginal data augmentation algorithms for posterior sampling in the Multivariate Multinomial Probit and Multivariate ordered Probit model. We consider the analysis of models for univariate and multivariate outcomes in the context of the latent variable inferential framework of Albert and Chib (1993) and Holmes and Held (2006). This will also provide trace plots, density plot and carter plot of the parameters. The theories considered to build this package are mentioned in the reference section and also discussed in vignette in details.


## Multivariate Multinomial Probit model: 

Let there are g nominal measures. For say qth measure we have q_p categories and we consider q_p -1 latent variables. The method consider to regress this latent variables and consider the maximum latent variable and consider the corresponding y value to be 1 or 0 otherwise. Here the latent variable follows multivariate normal with appropriate dimension. Hence the model is called probit. 

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

For better mixing, we consider Holmes Held approach, where we sample latent variable after marginalizing beta and we sample beta from its posterior distribution.


## Multivariate Multinomial Probit model: 

# Intended users: 

Mainly students working in Bayesian field with multivariate data. Specifically, this kind of models are mainly used in transportation system, Econometrics like labor economics (educational attainment), political economy (voter opinions), and health economics (consumers’ reliance on alternative sources of medical information) and so on. 





# References:

1.	XiaoZhang, W. JohnBoscardin, Thomas R.Belin (2008), “Bayesian Analysis of Multivariate Nominal measures using Multivariate Multinomial Probit Models”, Computational Statistics & Data Analysis, Volume 52, Issue 7, 15 March 2008, Pages 3697-3708

2.	Ivan Jeliazkov, Jennifer Graves, Mark Kutzbach (2008), “Fitting And Comparison of Models For Multivariate Ordinal Outcomes”, Advances in Econometrics: Bayesian Econometrics, 23
