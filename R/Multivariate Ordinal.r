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