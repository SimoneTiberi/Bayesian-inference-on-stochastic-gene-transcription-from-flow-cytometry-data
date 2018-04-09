# load the backgound error data:
load("background_error_data.RData")

N = 4 # 4 replicates.

# initial values for the logarithm of the hierarchical (mcmc.new) and for the hyper (new.hyper) parameters:
mcmc.new    =  matrix(1, nrow = N, ncol = 2)
new.hyper = matrix(1, nrow = 2, ncol = 2)

# sd of the proposal for the initial 10^5 iterations of the MCMC:
sd.prop = c(0.001, 0.001)

# MCMC iterations:
r = 1
R = 5.5 * 10^5

# matrix with the posterior chains for hierarchical parameters and covariance matrix for the ARW:
mcmc = cv = list()
ll.new = rep(NA, N) # log-likelihood
for(i in 1:N){
  mcmc[[i]]   = matrix(NA, nrow = R, ncol = 2)
  cv[[i]]     = diag(sd.prop^2)
  ll.new[i]   = sum( dnorm(Control[[i]], mean = exp(mcmc.new[i,1]), sd = exp(mcmc.new[i,2]), log = TRUE ) )
}

# list with the posterior chains for hypermean and hyperprecision parameters:
mcmc.hyper = list()
mcmc.hyper[[1]] = matrix(NA, nrow = R, ncol = 2) # Only first 5 par are hierarchical
mcmc.hyper[[2]] = matrix(NA, nrow = R, ncol = 2) # the error ones are not

library(mvtnorm)

# MCMC:
while(r <= R){
  ### hyperparameters sampling (Gibbs step):
  for(k in 1:2){
    mu.par = posterior.hyper.mu(par.cell = mcmc.new[,k], tau.hyper = new.hyper[2,k], index = k)
    new.hyper[1,k] = rnorm( 1, mean = mu.par[1], sd = sqrt(mu.par[2]) )
    
    tau.par = posterior.hyper.tau(par.cell = mcmc.new[,k], mu.hyper = new.hyper[1,k])
    new.hyper[2,k] = rgamma( 1, shape = tau.par[1], rate = tau.par[2] )
  }
  
  mcmc.hyper[[1]][r,] = new.hyper[1,] # 9 alpha of the prior Gammas
  mcmc.hyper[[2]][r,] = new.hyper[2,] # 9 beta of the prior Gammas

  ### hierarchical parameters sampling (Metropolis step):
  if(r %in% seq(10^5+100, R, 50) ){
    # after 10^5 iterations, every 50 iterations, update the covariance matrix of the ARW:
    for(i in 1:N){
      range = {10^5}:{r-1}
      cv[[i]] = diag(10^-9, 2) + 0.01 * cov( mcmc[[i]][range, ] ) # Kappa in theta 12
    }
  }
  
  for(i in 1:N){
    if( r < 10^5+100){
      # simple random walk for the first 10^5 + 100 iterations
      prop = rnorm(n = 2, mean = mcmc.new[i,], sd = sd.prop) 
    }else{
      # ARW for the following iterations
      prop = rmvnorm(n = 1, mean = mcmc.new[i,], sigma = cv[[i]] )      # Kappa in theta 12
    }
    
    # log-prior of the hierarchical parameters, given the hyperparameters:
    prior.prop = prior.hierarchy( prop,         new.hyper ) # PRIOR FOR NEW AND OLD VALUES
    prior.new  = prior.hierarchy( mcmc.new[i,], new.hyper ) # PRIOR FOR NEW AND OLD VALUES
      
    # log-likelihood of the background data, given the proposal hierarchical parameters for the i-th replicate:
    ll.prop  =  sum( dnorm(Control[[i]], mean = exp(prop[1]), sd = exp(prop[2]), log = TRUE ) )

    # exp(alpha) is the probability of accepting the proposal values "prop":
    alpha = min(0, ll.prop - ll.new[i] + prior.prop - prior.new )
    cond = (rbinom(n = 1, size = 1, prob = exp(alpha)) == 1)
    if(cond){
      # if accepted, I update the hierarchical values and the log-likelihood:
      mcmc.new[i,]  = prop
      ll.new[i] = ll.prop
    }
    mcmc[[i]][r,] = mcmc.new[i,]
  }
  
  r = r + 1
}
