###########################################################################################
# Hierarchical parameters:
###########################################################################################

# effective sample size (ess) estimation:
S.   = 10^5 
end. = 6*10^5

library(mcmcse)
ess_hierarchical = sapply(mcmc, function(x) ess(x[S.:end.,]))
ess_hierarchical

# min and average ess
min(ess_hierarchical); mean(ess_hierarchical)

# Heidelberg burn-in estimation:
library(coda)
heidel = lapply(mcmc, function(x) heidel.diag(x[1:{end.-1},], pvalue = 0.01))
heidel

max( sapply(heidel, function(x) x[,2] ) ) # max burn-in
mean( sapply(heidel, function(x) x[,1] ) ) # mean Stationarity test passed

###########################################################################################
# hyperparameters:
###########################################################################################

thin = 100
# effective sample size (ess) estimation and Heidelberg burn-in estimation:
S.   = 10^5/thin
end. = 6*10^5/thin

# on the hypermean
heidel_MEAN = heidel.diag(mcmc.hyper[[1]][1:end.,1:5], pvalue = 0.01)
ess_MEAN = ess(mcmc.hyper[[1]][S.:end.,1:5])

# on the hyperprecision
# sigma = sqrt(1/mcmc.hyper[[2]])
heidel_TAU = heidel.diag(mcmc.hyper[[2]][1:end.,1:5], pvalue = 0.01)
ess_TAU    = ess(mcmc.hyper[[2]][S.:end.,1:5])

# min and average ess
min(c(ess_MEAN,ess_TAU)); mean(c(ess_MEAN, ess_TAU))

max(c(heidel_MEAN[,2], heidel_TAU[,2])) # max burn-in
mean( c(heidel_MEAN[,1], heidel_TAU[,1]) ) # mean Stationarity test passed
