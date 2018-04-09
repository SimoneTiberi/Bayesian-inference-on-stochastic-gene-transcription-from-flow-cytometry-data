################################################################################################################
# Load the posterior output:
################################################################################################################

HYPER_MEAN = hierarchical_1 = hierarchical_2 = hierarchical_3 = hierarchical_4 = list()

# load the posterior output of the Tet_5 condition
load(...)
S.   = burn_in + 1; end. = R # total number of iterations of the MCMC

# load the hierarchical parameters:
hierarchical_1[[1]] = mcmc[[1]][seq(S.,end.,100),]
hierarchical_2[[1]] = mcmc[[2]][seq(S.,end.,100),]
hierarchical_3[[1]] = mcmc[[3]][seq(S.,end.,100),]
hierarchical_4[[1]] = mcmc[[4]][seq(S.,end.,100),]

# load the hypermean:
S. = S./thin; end. = end./thin
HYPER_MEAN[[1]]  =  mcmc.hyper[[1]][seq(S.,end.,1),]

# load the posterior output of the Tet_10 condition:
load(...)
S.   = burn_in + 1
end. = R # total number of iterations of the MCMC

# load the hierarchical parameters:
hierarchical_1[[2]] = mcmc[[1]][seq(S.,end.,100),]
hierarchical_2[[2]] = mcmc[[2]][seq(S.,end.,100),]
hierarchical_3[[2]] = mcmc[[3]][seq(S.,end.,100),]
hierarchical_4[[2]] = mcmc[[4]][seq(S.,end.,100),]

# load the hypermean:
S. = S./thin; end. = end./thin
HYPER_MEAN[[2]]  =  mcmc.hyper[[1]][seq(S.,end.,1),]

################################################################################################################
# Posterior mean for the coefficient of variation (CV) of the hierarchical parameters across replicates:
################################################################################################################

# Gather together the posterior chains for both experimental conditions:
mcmc = list()
for(i in 1:2){
  mcmc[[i]] = list()
  for(k in 1:7){ # for each parameter, I aggregate the MCMC across the 4 replicates.
    mcmc[[i]][[k]] = cbind( hierarchical_1[[i]][,k],  hierarchical_2[[i]][,k],  hierarchical_3[[i]][,k],  hierarchical_4[[i]][,k] )
  }
}

# compute the mean and standard deviation of the hierarchical parameters across replicates for each posterior sample:
mu = sigma = list()
for(i in 1:2){
  mu[[i]] =  sigma[[i]] = matrix(NA, nrow = 5000, ncol = 7)
  for(k in 1:7){ # for each parameter, I aggregate the MCMC across the 4 replicates.
    mu[[i]][,k]    = apply( exp( mcmc[[i]][[k]] ), 1, mean)
    sigma[[i]][,k] = apply( exp( mcmc[[i]][[k]] ), 1, sd)
  }
}

# posterior mean of the CV:
cv_tet_5  = apply( sigma[[1]]/mu[[1]], 2, mean)
cv_tet_10 = apply( sigma[[2]]/mu[[2]], 2, mean)
cv_tet_5; cv_tet_10

library(xtable)
xtable( rbind(cv_tet_5, cv_tet_10), digits = 2)

################################################################################################################
# Posterior mean for the fraction of mRNA which is transcribed from the OFF state:
################################################################################################################
x_0 = matrix(NA, nrow = 2, ncol = 4)
for(i in 1:2){
    a_0 = exp(hierarchical_1[[i]][,1]); a_1 = exp(hierarchical_1[[i]][,2])
    k_1 = exp(hierarchical_1[[i]][,3]); k_0 = exp(hierarchical_1[[i]][,4])
    mu_P = k_1/{k_0 + k_1}
    x_1 = a_0 * (1-mu_P) / (a_1 * mu_P + a_0 * (1-mu_P) )
    
    a_0 = exp(hierarchical_2[[i]][,1]); a_1 = exp(hierarchical_2[[i]][,2])
    k_1 = exp(hierarchical_2[[i]][,3]); k_0 = exp(hierarchical_2[[i]][,4])
    mu_P = k_1/{k_0 + k_1}
    x_2 = a_0 * (1-mu_P) / (a_1 * mu_P + a_0 * (1-mu_P) )
    
    a_0 = exp(hierarchical_3[[i]][,1]); a_1 = exp(hierarchical_3[[i]][,2])
    k_1 = exp(hierarchical_3[[i]][,3]); k_0 = exp(hierarchical_3[[i]][,4])
    mu_P = k_1/{k_0 + k_1}
    x_3 = a_0 * (1-mu_P) / (a_1 * mu_P + a_0 * (1-mu_P) )
    
    a_0 = exp(hierarchical_4[[i]][,1]); a_1 = exp(hierarchical_4[[i]][,2])
    k_1 = exp(hierarchical_4[[i]][,3]); k_0 = exp(hierarchical_4[[i]][,4])
    mu_P = k_1/{k_0 + k_1}
    x_4 = a_0 * (1-mu_P) / (a_1 * mu_P + a_0 * (1-mu_P) )
    
    x_0[i,] = c(mean(x_1), mean(x_2), mean(x_3), mean(x_4) )
}

library(xtable)
xtable( x_0, digits = 2)

################################################################################################################
# 0.95 level HPD CIs for the hierarchical parameters and reparametrizations:
################################################################################################################
library(coda) # we use the HPDinterval function of the package coda.

# Compute the HPD CIs for the hierarchical parameters in each replicate:
CI_hier = list()
for(i in (1:2)){
  CI_hier[[i]] = matrix(NA, nrow = 7, ncol = 8)
  for(k in 1:7){
    CI_hier[[i]][k, 1:2 ] = HPDinterval( mcmc( exp(hierarchical_1[[i]][,k]) ) )
    CI_hier[[i]][k, 3:4 ] = HPDinterval( mcmc( exp(hierarchical_2[[i]][,k]) ) )
    CI_hier[[i]][k, 5:6 ] = HPDinterval( mcmc( exp(hierarchical_3[[i]][,k]) ) )
    CI_hier[[i]][k, 7:8 ] = HPDinterval( mcmc( exp(hierarchical_4[[i]][,k]) ) )
  }
}

# Compute the HPD CIs for the reparametrizations in each replicate:
CI_ratio = list()
for(i in (1:2)){
  CI_ratio[[i]] = matrix(NA, nrow = 6, ncol = 8)
  
  # 100 * alpha_0 / alpha_1
  k = 1
  a_0 = exp(hierarchical_1[[i]][,1]); a_1 = exp(hierarchical_1[[i]][,2])
  x_1 = a_0/{a_1}
  a_0 = exp(hierarchical_2[[i]][,1]); a_1 = exp(hierarchical_2[[i]][,2])
  x_2 = a_0/{a_1}
  a_0 = exp(hierarchical_3[[i]][,1]); a_1 = exp(hierarchical_3[[i]][,2])
  x_3 = a_0/{a_1}
  a_0 = exp(hierarchical_4[[i]][,1]); a_1 = exp(hierarchical_4[[i]][,2])
  x_4 = a_0/{a_1}
  
  CI_ratio[[i]][k, 1:2 ] = 100 * HPDinterval( mcmc( x_1 ) )
  CI_ratio[[i]][k, 3:4 ] = 100 * HPDinterval( mcmc( x_2 ) )
  CI_ratio[[i]][k, 5:6 ] = 100 * HPDinterval( mcmc( x_3 ) )
  CI_ratio[[i]][k, 7:8 ] = 100 * HPDinterval( mcmc( x_4 ) )
  
  # 100 * mu_P = 100 * k_1/(k_1 + k_0)
  k = 2
  k_1 = exp(hierarchical_1[[i]][,3]); k_0 = exp(hierarchical_1[[i]][,4])
  x_1 = k_1/{k_1 + k_0}
  k_1 = exp(hierarchical_2[[i]][,3]); k_0 = exp(hierarchical_2[[i]][,4])
  x_2 = k_1/{k_1 + k_0}
  k_1 = exp(hierarchical_3[[i]][,3]); k_0 = exp(hierarchical_3[[i]][,4])
  x_3 = k_1/{k_1 + k_0}
  k_1 = exp(hierarchical_4[[i]][,3]); k_0 = exp(hierarchical_4[[i]][,4])
  x_4 = k_1/{k_1 + k_0}
  
  CI_ratio[[i]][k, 1:2 ] = 100 * HPDinterval( mcmc( x_1 ) )
  CI_ratio[[i]][k, 3:4 ] = 100 * HPDinterval( mcmc( x_2 ) )
  CI_ratio[[i]][k, 5:6 ] = 100 * HPDinterval( mcmc( x_3 ) )
  CI_ratio[[i]][k, 7:8 ] = 100 * HPDinterval( mcmc( x_4 ) )
  
  # 100 * 1/k_1
  k = 3
  k_1 = exp(hierarchical_1[[i]][,3]);
  x_1 = 1/k_1
  k_1 = exp(hierarchical_2[[i]][,3]);
  x_2 = 1/k_1
  k_1 = exp(hierarchical_3[[i]][,3]);
  x_3 = 1/k_1
  k_1 = exp(hierarchical_4[[i]][,3]);
  x_4 = 1/k_1
  
  CI_ratio[[i]][k, 1:2 ] = 100 * HPDinterval( mcmc( x_1 ) )
  CI_ratio[[i]][k, 3:4 ] = 100 * HPDinterval( mcmc( x_2 ) )
  CI_ratio[[i]][k, 5:6 ] = 100 * HPDinterval( mcmc( x_3 ) )
  CI_ratio[[i]][k, 7:8 ] = 100 * HPDinterval( mcmc( x_4 ) )
  
  # 100 * 1/k_0
  k = 4
  k_0 = exp(hierarchical_1[[i]][,4])
  x_1 = 1/k_0
  k_0 = exp(hierarchical_2[[i]][,4])
  x_2 = 1/k_0
  k_0 = exp(hierarchical_3[[i]][,4])
  x_3 = 1/k_0
  k_0 = exp(hierarchical_4[[i]][,4])
  x_4 = 1/k_0
  
  CI_ratio[[i]][k, 1:2 ] = 100 * HPDinterval( mcmc( x_1 ) )
  CI_ratio[[i]][k, 3:4 ] = 100 * HPDinterval( mcmc( x_2 ) )
  CI_ratio[[i]][k, 5:6 ] = 100 * HPDinterval( mcmc( x_3 ) )
  CI_ratio[[i]][k, 7:8 ] = 100 * HPDinterval( mcmc( x_4 ) )
  
  # mu_X (the mean latent population of mRNA molecules)
  k = 5
  a_0 = exp(hierarchical_1[[i]][,1]);  a_1 = exp(hierarchical_1[[i]][,2])
  k_1 = exp(hierarchical_1[[i]][,3]); k_0 = exp(hierarchical_1[[i]][,4])
  x_1 = (a_1 * k_1 + a_0 * k_0)/(k_1 + k_0)
  a_0 = exp(hierarchical_2[[i]][,1]);  a_1 = exp(hierarchical_2[[i]][,2])
  k_1 = exp(hierarchical_2[[i]][,3]); k_0 = exp(hierarchical_2[[i]][,4])
  x_2 = (a_1 * k_1 + a_0 * k_0)/(k_1 + k_0)
  a_0 = exp(hierarchical_3[[i]][,1]);  a_1 = exp(hierarchical_3[[i]][,2])
  k_1 = exp(hierarchical_3[[i]][,3]); k_0 = exp(hierarchical_3[[i]][,4])
  x_3 = (a_1 * k_1 + a_0 * k_0)/(k_1 + k_0)
  a_0 = exp(hierarchical_4[[i]][,1]);  a_1 = exp(hierarchical_4[[i]][,2])
  k_1 = exp(hierarchical_4[[i]][,3]); k_0 = exp(hierarchical_4[[i]][,4])
  x_4 = (a_1 * k_1 + a_0 * k_0)/(k_1 + k_0)
  
  CI_ratio[[i]][k, 1:2 ] = HPDinterval( mcmc( x_1 ) )
  CI_ratio[[i]][k, 3:4 ] = HPDinterval( mcmc( x_2 ) )
  CI_ratio[[i]][k, 5:6 ] = HPDinterval( mcmc( x_3 ) )
  CI_ratio[[i]][k, 7:8 ] = HPDinterval( mcmc( x_4 ) )
  
  # mu_Y (the mean of the observations)
  k = 6
  a_0 = exp(hierarchical_1[[i]][,1]);  a_1 = exp(hierarchical_1[[i]][,2])
  k_1 = exp(hierarchical_1[[i]][,3]); k_0 = exp(hierarchical_1[[i]][,4])
  kappa = exp(hierarchical_1[[i]][,5])
  mu_epsilon = exp(hierarchical_1[[i]][,6])
  mu_X = (a_1 * k_1 + a_0 * k_0)/(k_1 + k_0)
  x_1 = kappa * mu_X + mu_epsilon
  a_0 = exp(hierarchical_2[[i]][,1]);    a_1 = exp(hierarchical_2[[i]][,2])
  k_1 = exp(hierarchical_2[[i]][,3]);    k_0 = exp(hierarchical_2[[i]][,4])
  kappa = exp(hierarchical_2[[i]][,5]);  mu_epsilon = exp(hierarchical_2[[i]][,6])
  mu_X = (a_1 * k_1 + a_0 * k_0)/(k_1 + k_0)
  x_2 = kappa * mu_X + mu_epsilon
  a_0 = exp(hierarchical_3[[i]][,1]);  a_1 = exp(hierarchical_3[[i]][,2])
  k_1 = exp(hierarchical_3[[i]][,3]); k_0 = exp(hierarchical_3[[i]][,4])
  kappa = exp(hierarchical_3[[i]][,5])
  mu_epsilon = exp(hierarchical_3[[i]][,6])
  mu_X = (a_1 * k_1 + a_0 * k_0)/(k_1 + k_0)
  x_3 = kappa * mu_X + mu_epsilon
  a_0 = exp(hierarchical_4[[i]][,1]);  a_1 = exp(hierarchical_4[[i]][,2])
  k_1 = exp(hierarchical_4[[i]][,3]); k_0 = exp(hierarchical_4[[i]][,4])
  kappa = exp(hierarchical_4[[i]][,5])
  mu_epsilon = exp(hierarchical_4[[i]][,6])
  mu_X = (a_1 * k_1 + a_0 * k_0)/(k_1 + k_0)
  x_4 = kappa * mu_X + mu_epsilon
  
  CI_ratio[[i]][k, 1:2 ] = HPDinterval( mcmc( x_1 ) )
  CI_ratio[[i]][k, 3:4 ] = HPDinterval( mcmc( x_2 ) )
  CI_ratio[[i]][k, 5:6 ] = HPDinterval( mcmc( x_3 ) )
  CI_ratio[[i]][k, 7:8 ] = HPDinterval( mcmc( x_4 ) )
}

library(xtable)
xtable(rbind(CI_hier[[1]],CI_ratio[[1]]), digits = 1) # Tet5
xtable(rbind(CI_hier[[2]],CI_ratio[[2]]), digits = 1) # Tet10


################################################################################################################
# 0.95 level HPD CIs for the exponential of the hypermeans:
################################################################################################################
library(coda)

CI_Hyper_Mu = matrix(NA, nrow = 5, ncol = 4)
for(i in 1:2){
  for(k in 1:5){
    CI_Hyper_Mu[k, c(2*i-1,2*i) ]    = HPDinterval( mcmc( exp(HYPER_MEAN[[i]][,k]) ) )
  }
}

library(xtable)
xtable(CI_Hyper_Mu[1:5,], digits = 1)
# the first 2 columns refer to Tet5, the last 2 columns refer to Tet10.
