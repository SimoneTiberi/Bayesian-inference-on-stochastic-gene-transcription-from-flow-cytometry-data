library(mvtnorm)

# Run the MCMC for R iterations:
while(r <= R){
  # every 10^4 iterations, I allocate more space for the posterior values of the next 10^4 iterations.
  if(r %in% seq(10^4, R, 10^4) ){
    for(k in 1:K){
      mcmc[[k]] = rbind(mcmc[[k]], matrix(NA, nrow = 10^4, ncol = 7))
    }  
    
    mcmc.hyper[[1]] = rbind( mcmc.hyper[[1]], matrix(NA, nrow = 10^4/thin, ncol = 7))
    mcmc.hyper[[2]] = rbind( mcmc.hyper[[2]], matrix(NA, nrow = 10^4/thin, ncol = 7))
  }
  
  ###################################################################################################
  ## hyperparameters sampling
  ###################################################################################################
  
  # the first four hypermeans and hyperprecisions follow a Gibbs step:
  for(k in 1:4){
    mu.par = posterior.hyper.mu(par.cell = mcmc.new[,k], tau.hyper = new.hyper[2,k], index = k, N = K)
    new.hyper[1,k] = rnorm( 1, mean = mu.par[1], sd = sqrt(mu.par[2]) )
    
    tau.par = posterior.hyper.tau(par.cell = mcmc.new[,k], mu.hyper = new.hyper[1,k], N = K)
    new.hyper[2,k] = rgamma( 1, shape = tau.par[1], rate = tau.par[2] )
  }
  
  # the fifth hypermean and hyperprecision, for kappa, follow a Metropolis-Hasting step:
  # propose the hypermean:
  mu.par = posterior.hyper.mu(par.cell = mcmc.new[,5], tau.hyper = new.hyper[2,5], index = 5, N = K)
  prop_mu = rnorm( 1, mean = mu.par[1], sd = sqrt(mu.par[2]) )
  
  int_prop = 1-plnorm(1, meanlog = prop_mu,        sdlog = sqrt(1/new.hyper[2,5]))
  int_new  = 1-plnorm(1, meanlog = new.hyper[1,5], sdlog = sqrt(1/new.hyper[2,5]))
  # accept/reject:
  alpha = min(1,  int_new/int_prop)
  # sample whether to accept the proposed values, with probability exp(alpha)
  cond = {rbinom(n = 1, size = 1, prob = alpha) == 1}
  if(cond){
    new.hyper[1,5] = prop_mu
  }

  # propose the hyperprecision:
  tau.par = posterior.hyper.tau(par.cell = mcmc.new[,5], mu.hyper = new.hyper[1,5], N = K)
  prop_tau = rgamma( 1, shape = tau.par[1], rate = tau.par[2] )
  
  int_prop = 1-plnorm(1, meanlog = new.hyper[1,5], sdlog = sqrt(1/prop_tau))
  int_new  = 1-plnorm(1, meanlog = new.hyper[1,5], sdlog = sqrt(1/new.hyper[2,5]))
  # accept/reject:
  alpha = min(1,  int_new/int_prop)
  # sample whether to accept the proposed values, with probability exp(alpha)
  cond = {rbinom(n = 1, size = 1, prob = alpha) == 1}
  if(cond){
    new.hyper[2,5] = prop_tau
  }
  
  # I save the hyperparameter values every thin observations:
  if(r/thin == round(r/thin)){
    mcmc.hyper[[1]][r/thin,] = new.hyper[1,]
    mcmc.hyper[[2]][r/thin,] = new.hyper[2,]
  }
  
  ###################################################################################################
  ## Hierarchical parameters sampling
  ###################################################################################################
  
  for(k in 1:K){  # For each replicate:

    # Propose new values:
    if( r < 2*10^2){
      # Propose from a simple random walk for the first 200 iterations.
      prop = rnorm(n = 7, mean = mcmc.new[k,], sd = sd.prop) 
    }else{
       # Every 100 iterations, compute the covariance matrix for the ARW proposal.
      if(r %in% seq(10^2, 10^6, 10^2) ){
        if(r > 10^5){
          range = 0.5 * 10^5:{r-1}
        }else{
          range = 101:{r-1}
        }
        # Compute the ARW proposal matrixes for the two blocks:
        cv.a[[k]]      = diag(10^-9, 5) + 0.1 * cov( mcmc[[k]][range, 1:5] ) # Kappa in theta 12
        cv.error[[k]]  = diag(10^-9, 2) +  1  * cov( mcmc[[k]][range, 6:7] )
        rm(range)
      }
      # propose from an ARW for r > 200, in two blocks:
      prop[1:5] = rmvnorm(n = 1, mean = mcmc.new[k,1:5], sigma = cv.a[[k]] )      # Kappa in theta 12
      prop[6:7] = rmvnorm(n = 1, mean = mcmc.new[k,6:7], sigma = cv.error[[k]] )
    }
    
    ### accept/reject the proposal values for the first block.
    # condition: alpha_1 > alpha_0 and kappa > 1 (i.e. log(kappa) > 0).
    cond = {prop[2] > prop[1]} & {prop[5] > 0}
    if(cond){
      # Compute the log-prior for the proposed and current values.
      prior.a.new = prior.hierarchy( prop[1:5],       new.hyper )
      prior.a.old = prior.hierarchy( mcmc.new[k,1:5], new.hyper )
      
      x = c(prop[1:5], mcmc.new[k,6:7] )
      
      # Estimate the marginal log-likelihood of the observations, using 10^5 particles.
      mcmc.ll.prop = smc_approx_unbiased(Y[[k]], par = exp(x), K = 10^5 )
      
      # alpha is the logarithm of the probability of accepting the proposed value.
      alpha = min(0,  mcmc.ll.prop - mcmc.ll.new[k] + prior.a.new - prior.a.old)
      # sample whether to accept the proposed values, with probability exp(alpha)
      cond = {rbinom(n = 1, size = 1, prob = exp(alpha)) == 1}
      
      if(cond){
        # if proposed values are accepted, we update the current values of the MCMC and the marginal log-likelihood
        mcmc.new[k,1:5]   = prop[1:5]
        mcmc.ll.new[k] = mcmc.ll.prop
      }
    }
    
    ### accept/reject the proposal values for the second block.
    # Compute the (informative) log-prior for the proposed and current values.
    prior.error.new =  prior.log.err( prop[6:7],       index = k )
    prior.error.old =  prior.log.err( mcmc.new[k,6:7], index = k )
    
    x = c(mcmc.new[k, 1:5], prop[6:7] )
    
    # Estimate the marginal log-likelihood of the observations, using 10^5 particles.
    mcmc.ll.prop = smc_approx_unbiased(Y[[k]], par = exp(x), K = 10^5 )
    
    # alpha is the logarithm of the probability of accepting the proposed value.
    alpha = min(0,  mcmc.ll.prop - mcmc.ll.new[k] + prior.error.new - prior.error.old)
    # sample whether to accept the proposed values, with probability exp(alpha)
    cond = {rbinom(n = 1, size = 1, prob = exp(alpha)) == 1}
    
    if(cond){
      # if proposed values are accepted, we update the current values of the MCMC and the marginal log-likelihood
      mcmc.new[k,6:7]   = prop[6:7]
      mcmc.ll.new[k] = mcmc.ll.prop
    }
    
    # store the current hierarchical values in the mcmc object
    mcmc[[k]][r,]    = mcmc.new[k,]
  }
  
  # Increase one iteration to the MCMC  
  r = r + 1
}
