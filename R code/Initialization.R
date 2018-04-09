# load the prior for the measurement error parameters:
load("prior for measurement error.RData")

# load the data: the observations are stored in the list Y,
 load("Data Tet5.RData") # Data for cells stimulated with 5 ng/ML of tetracycline.
# or
# load("Data Tet10.RData") # Data for cells stimulated with 10 ng/ML of tetracycline.

# Set the number of replicates
K = 4

# Starting values for the hierchical parameters:
mcmc.new = matrix( c( 2.5, 8.3, -1.5, 3.2, 2.3), 
                   nrow = K, ncol = 5, byrow = TRUE)

# Add the starting values for the measurement error parameters, set to their prior mean:
mcmc.new = cbind(mcmc.new, MEAN)

# sd for the proposal of the hierarchical parameters in the initial 200 iterations of the MCMC:
sd.prop = c(0.001, 0.001, 0.01, 0.01, 0.001, 0.001, 0.001)

# Allocate the vector for the log-likelihood for the K replicates.
mcmc.ll.new = rep(NA, K)

mcmc = cv.a = cv.error = list()
for(k in 1:K){
  # Allocate the matrix for the posterior values of the hierarchical parameters, initially for 10^4 iterations.
  mcmc[[k]] = matrix(NA, nrow = 10^4, ncol = 7)

  # Initialize the covariance matrixes of the ARW for both blocks, for each replicate
  cv.a[[k]]  = diag(sd.prop[c(1:5)]^2)
  cv.error[[k]]  = diag(sd.prop[6:7]^2)
  
  # Estimate the marginal likelihood of the data, given the parameters of each replicate.
  mcmc.ll.new[k] = smc_approx_unbiased(Y[[k]], par = exp(mcmc.new[k,]), K = 10^5 )
}

# Set the iterations indicator r to 1 and the total number of iterations to 6*10^5.
r = 1
R = 6*10^5

# thinning factor for the hyperparameters.
thin = 100
# ALlocate the matrixes for the posterior values of the hypermean (mcmc.hyper[[1]]) and hyperprecision (mcmc.hyper[[2]]).
mcmc.hyper = list()
mcmc.hyper[[1]] = matrix(NA, nrow = 10^4/thin, ncol = 5)
mcmc.hyper[[2]] = matrix(NA, nrow = 10^4/thin, ncol = 5)

# Initialize the hyperparameters:
new.hyper = matrix(1, nrow = 2, ncol = 5)
