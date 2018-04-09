
# posterior values for the sampling of the hypermean
posterior.hyper.mu = function(par.cell, mu.hyper, tau.hyper, index){ 
  mu_0  =  0     
  lambda = 10000
  
  tau_0 = tau.hyper / lambda ;
  mu = (tau_0 * mu_0 + tau.hyper * sum( par.cell ) ) / ( tau_0 + N * tau.hyper )
  sigma2 = 1/(tau_0 + N * tau.hyper)
  c(mu, sigma2)
}

# posterior values for the sampling of the hyperprecision
posterior.hyper.tau = function(par.cell, mu.hyper, tau.hyper){ 
  alpha_0 = 0.001; beta_0  = 0.001
  
  alpha = alpha_0 + N/2
  beta = beta_0 + sum( (par.cell - mu.hyper)^2 )/2
  c(alpha, beta)
}

# log-prior of the hierarchical parameters, given the hyperparameters:
prior.hierarchy = function(par.cell, par.hyper){ # par.hyper matrix 2 * 9
  sum( dnorm(par.cell, mean = par.hyper[1,], sd = 1/sqrt(par.hyper[2,]), log = TRUE) )
}

