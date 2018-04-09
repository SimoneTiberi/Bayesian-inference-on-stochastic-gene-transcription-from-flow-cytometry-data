# Conjugate posterior sampling for the hypermean parameter:
posterior.hyper.mu = function(par.cell, mu.hyper, tau.hyper, index, N){ 
  mu_0  =  0     
  lambda = 10000

  tau_0 = tau.hyper / lambda ;
  mu = (tau_0 * mu_0 + tau.hyper * sum( par.cell ) ) / ( tau_0 + N * tau.hyper )
  sigma2 = 1/(tau_0 + N * tau.hyper)
  c(mu, sigma2)
}

# Conjugate posterior sampling for the hyperprecision parameter:
posterior.hyper.tau = function(par.cell, mu.hyper, tau.hyper, N){ 
  alpha_0 = 0.001
  beta_0  = 0.001
  
  alpha = alpha_0 + N/2
  beta = beta_0 + sum( (par.cell - mu.hyper)^2 )/2
  c(alpha, beta)
}

# informative log-prior for the measurement error mean and standard deviation, 
# both in the logarithmic space,
# given the constant prior in MEAN and SD (with distinct values for each replicate)
prior.log.err = function(theta, index){
  sum( dnorm(theta, mean = MEAN[index,], sd = SD[index,], log = TRUE))
}

# prior for the hierarchical parameters, in the logarithmic space,
# given the correspoinding hyperparameters.
prior.hierarchy = function(par.cell, par.hyper){
  sum( dnorm(par.cell, mean = par.hyper[1,], sd = 1/sqrt(par.hyper[2,]), log = TRUE) )
}

# Function to estimate the marginal log-likelhood of the data "Y" with "K" particles, given parameters "par".
smc_approx_unbiased = function(Y, par, K){
  a0 = par[1]; a1 = par[2]; 
  kON = par[3]; kOFF = par[4];
  kappa = par[5];  mu = par[6]
  sigma = par[7]
  
  # X sampling
  P_new = rbeta(K, kON, kOFF)
  lambda.new = a0 + (a1-a0) * P_new # the rate of the Poiss
  X.pois = rpois(K, lambda.new)
  
  # Since X is a discrete random variable, we can summarize its values in how how many times each value was sampled:
  # This trick provides a massive computational gain.
  X.unique = unique(X.pois) # I consider all values sampled.
  X.times = c()
  for(i in 1:length(X.unique)){
    X.times[i] = sum(X.unique[i] == X.pois) # We count how many times each unique value was sampled.
  }
  
  res = 0
  # I set the result to 0 (called l_0 in the supplementary).
  # at every iteration, we will update res.
  
  for(i in 1:length(Y)){
    # for every data point Y[i], we compute the corresponding error across all unique Xs sampled above,
    error =  Y[i] - kappa * X.unique 
    
    # and measure the corresponding log-density.
    w.log = dnorm(x = error, mean = mu, sd = sigma, log = TRUE)
    
    # trick to avoid underflow with small density values.
    w.max = max(w.log)
    w = w.log - w.max

    # We compute the G_is values, muliplying the density of each value by the number of times they were sampled.
    G_times = {exp(w) * X.times}
    sum_G = sum(G_times)
    
    if( is.na(sum_G) ){
      print("NA in sum_G: log-like is -Inf")
      return(-Inf)
    }
    if( sum_G == -Inf){
      print("sum_G is -Inf: log-like is -Inf")
      return(-Inf)
    }
    
    if(sum_G > 0){
      p = G_times/sum_G
    }else{
      p = X.times/sum(X.times)
    }
    # I update the result (marginal log-likelihood):
    # w.max is added because it was removed from all w.log to avoid underflow.
    res = res + log( sum_G ) + w.max
    
    # We sample a value of X, via "sample.one" function, to remove from the estimation for the next observations:
    index = sample.one(p)
    # We decrease by 1 the number of times the selected value was sampled:
    X.times[index] = X.times[index] - 1
  }
  
  if( is.na(res) ){ 
    print("NA in res: log-like is -Inf")
    return(-Inf)
  }
  if(res == Inf){
    print("res is -Inf: log-like is -Inf")
    return(-Inf)
  }
  
  res
}

# function, to sample one value between 1:length(xs), with probability proportional to xs.
sample.one <- function(xs) {
  ys <- cumsum(xs)
  u <- runif(1)*ys[length(ys)]
  return(binary.search(u,ys))
}

binary.search <- function(u,ys) {
  n <- length(ys)
  if (n == 1) return(1)
  mid <- floor(n/2)
  if (u <= ys[mid]) {
    return(binary.search(u,ys[1:mid]))
  } else {
    return(mid + binary.search(u,ys[(mid+1):n]))
  }
}
