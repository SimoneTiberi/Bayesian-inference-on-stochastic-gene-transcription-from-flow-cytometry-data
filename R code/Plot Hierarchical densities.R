################################################################################################################
# Load the posterior output:
################################################################################################################

hierarchical_1 = hierarchical_2 = hierarchical_3 = hierarchical_4 = list()

# load the posterior output of the Tet_5 condition
load(...)
S.   = burn_in + 1
end. = R # total number of iterations of the MCMC

hierarchical_1[[1]] = mcmc[[1]][seq(S.,end.,100),]
hierarchical_2[[1]] = mcmc[[2]][seq(S.,end.,100),]
hierarchical_3[[1]] = mcmc[[3]][seq(S.,end.,100),]
hierarchical_4[[1]] = mcmc[[4]][seq(S.,end.,100),]

# mean and sd of the experimental data:
E_Y = SD_Y =  matrix(NA, nrow = 2, ncol = 4)
for(k in 1:K){
  E_Y[1,k] = mean(Y[[k]])
  SD_Y[1,k] = mean(Y[[k]])
}

# load the posterior output of the Tet_10 condition:
load(...)
S.   = burn_in + 1
end. = R # total number of iterations of the MCMC

load("hierarchical.R")
hierarchical_1[[2]] = mcmc[[1]][seq(S.,end.,100),]
hierarchical_2[[2]] = mcmc[[2]][seq(S.,end.,100),]
hierarchical_3[[2]] = mcmc[[3]][seq(S.,end.,100),]
hierarchical_4[[2]] = mcmc[[4]][seq(S.,end.,100),]

# mean and sd of the experimental data:
for(k in 1:K){
  E_Y[2,k] = mean(Y[[k]])
  SD_Y[2,k] = mean(Y[[k]])
}

################################################################################################################
# Posterior densities for the hierarchical parameters:
################################################################################################################
par(mfrow = c(4,2))
cols = 1:2
# black:  5 nl/mL of tetracycline
# red:   10 nl/mL of tetracycline

# alpha_0
plot(NA,
     #main = expression(theta[1]^{(k)}),
     main = expression(tilde(alpha)[0]^{(k)}),
     xlab = "", ylab = "", cex.main = 1.5, cex.axis = 1.5, cex.lab = 0.8,
     xlim = c(0,300), ylim = c(0, 0.035))

for(i in 1:2){
  lines(density(exp(hierarchical_1[[i]][,1]),  n = 10^3, adjust = 2), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_2[[i]][,1]),  n = 10^3, adjust = 2), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_3[[i]][,1]),  n = 10^3, adjust = 2), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_4[[i]][,1]),  n = 10^3, adjust = 2), col = cols[i], lty = i, lwd = 1.5 )
}

# alpha_1
plot(NA, 
     xlim = c(0,20000), ylim = c(0, 0.00045),
     main = expression(tilde(alpha)[1]^{(k)}),
     xlab = "", ylab = "", cex.main = 1.5, cex.axis = 1.5, cex.lab = 0.8)


for(i in 1:2){
  lines(density(exp(hierarchical_1[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_2[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_3[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_4[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
}

# k_1
plot(NA,
     main = expression(tilde(k)[1]^{(k)}),
     xlab = "", ylab = "", cex.main = 1.5, cex.axis = 1.5, cex.lab = 0.8,
     xlim = c(0.05, 0.8), ylim = c(0, 19))

for(i in 1:2){
  lines(density(exp(hierarchical_1[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_2[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_3[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_4[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
}

# k_0
plot(NA,
     main = expression(tilde(k)[0]^{(k)}),
     xlab = "", ylab = "", cex.main = 1.5, cex.axis = 1.5, cex.lab = 0.8,
     xlim = c(1,15), ylim = c(0, 0.95))

for(i in 1:2){
  lines(density(exp(hierarchical_1[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_2[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_3[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_4[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
}

# kappa
plot(NA,
     main = expression(kappa^{(k)}),
     xlab = "", ylab = "", cex.main = 1.5, cex.axis = 1.5, cex.lab = 0.8,
     xlim = c(0,100), ylim = c(0, 0.15))

for(i in 1:2){
  lines(density(exp(hierarchical_1[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_2[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_3[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_4[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
}


# mu_epsilon
plot(NA,
     main = expression(mu[epsilon]^{(k)}),
     xlab = "", ylab = "", cex.main = 1.5, cex.axis = 1.5, cex.lab = 0.8,
     xlim = c(750,1000), ylim = c(0, 0.14))

for(i in 1:2){
  lines(density(exp(hierarchical_1[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_2[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_3[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_4[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
}

# sigma_epsilon
plot(NA,
     main = expression(sigma[epsilon]^{(k)}),
     xlab = "", ylab = "", cex.main = 1.5, cex.axis = 1.5, cex.lab = 0.8,
     xlim = c(290,420), ylim = c(0, 0.2))

for(i in 1:2){
  lines(density(exp(hierarchical_1[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_2[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_3[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  lines(density(exp(hierarchical_4[[i]][,1]),  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
}

################################################################################################################
# Posterior densities for reparametrizations of the hierarchical parameters
################################################################################################################
# alpha_0/alpha_1
par(mfrow = c(2,2))

plot(NULL, main = expression(alpha[0]^{(k)}/alpha[1]^{(k)}),
     xlab = "", ylab = "", cex.main = 2, cex.axis = 2, cex.lab = 1,
     xlim = c(0,0.04), ylim = c(0, 150))

for(i in 1:2){
a_0 = exp(hierarchical_1[[i]][,1]); a_1 = exp(hierarchical_1[[i]][,2])
  x_1 = a_0/{a_1}
  
  a_0 = exp(hierarchical_2[[i]][,1]); a_1 = exp(hierarchical_2[[i]][,2])
  x_2 = a_0/{a_1}
  
  a_0 = exp(hierarchical_3[[i]][,1]); a_1 = exp(hierarchical_3[[i]][,2])
  x_3 = a_0/{a_1}
  
  a_0 = exp(hierarchical_4[[i]][,1]); a_1 = exp(hierarchical_4[[i]][,2])
  x_4 = a_0/{a_1}
  
  lines(density( x_1,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
  lines(density( x_2,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
  lines(density( x_3,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
  lines(density( x_4,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
}

# E(P)
plot(NULL, main = expression(mu[P]^{(k)}),
     xlab = "", ylab = "", cex.main = 2, cex.axis = 2, cex.lab = 1,
     xlim = c(0.01,0.16), ylim = c(0, 55))

for(i in 1:2){
  a_0 = exp(hierarchical_1[[i]][,3]); a_1 = exp(hierarchical_1[[i]][,4])
  x_1 = a_0/{a_0 + a_1}
  
  a_0 = exp(hierarchical_2[[i]][,3]); a_1 = exp(hierarchical_2[[i]][,4])
  x_2 = a_0/{a_0 + a_1}
  
  a_0 = exp(hierarchical_3[[i]][,3]); a_1 = exp(hierarchical_3[[i]][,4])
  x_3 = a_0/{a_0 + a_1}
  
  a_0 = exp(hierarchical_4[[i]][,3]); a_1 = exp(hierarchical_4[[i]][,4])
  x_4 = a_0/{a_0 + a_1}
  
  lines(density( x_1,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
  lines(density( x_2,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
  lines(density( x_3,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
  lines(density( x_4,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
}

# 1/k_1
plot(NULL, main = expression(1/tilde(k)[1]^{(k)}), 
     xlab = "", ylab = "", cex.main = 2, cex.axis = 2, cex.lab = 1,
     xlim = c(1,10.5), ylim = c(0, 2.7))

for(i in 1:2){
  a_0 = exp(hierarchical_1[[i]][,3]); a_1 = exp(hierarchical_1[[i]][,4])
  x_1 = 1/a_0
  
  a_0 = exp(hierarchical_2[[i]][,3]); a_1 = exp(hierarchical_2[[i]][,4])
  x_2 = 1/a_0
  
  a_0 = exp(hierarchical_3[[i]][,3]); a_1 = exp(hierarchical_3[[i]][,4])
  x_3 = 1/a_0
  
  a_0 = exp(hierarchical_4[[i]][,3]); a_1 = exp(hierarchical_4[[i]][,4])
  x_4 = 1/a_0
  
  lines(density( x_1,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
  lines(density( x_2,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
  lines(density( x_3,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
  lines(density( x_4,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
}

# 1/k_0
plot(NULL, main = expression(1/tilde(k)[0]^{(k)}), 
     xlab = "", ylab = "", cex.main = 2, cex.axis = 2, cex.lab = 1,
     xlim = c(0,0.6), ylim = c(0, 9))

for(i in 1:2){
  a_0 = exp(hierarchical_1[[i]][,3]); a_1 = exp(hierarchical_1[[i]][,4])
  x_1 = 1/a_1
  
  a_0 = exp(hierarchical_2[[i]][,3]); a_1 = exp(hierarchical_2[[i]][,4])
  x_2 = 1/a_1
  
  a_0 = exp(hierarchical_3[[i]][,3]); a_1 = exp(hierarchical_3[[i]][,4])
  x_3 = 1/a_1
  
  a_0 = exp(hierarchical_4[[i]][,3]); a_1 = exp(hierarchical_4[[i]][,4])
  x_4 = 1/a_1
  
  lines(density( x_1,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
  lines(density( x_2,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
  lines(density( x_3,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
  lines(density( x_4,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
}

# x_0, fraction of mRNA transcribed in the OFF state.
plot(NULL, main = expression(tilde(x)[0]^{(k)}), 
     xlab = "", ylab = "", cex.main = 2, cex.axis = 2, cex.lab = 1,
     xlim = c(0,0.5), ylim = c(0, 28))

for(i in 1:7){
  if(i %in% SELECT){
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
    
    print(round( c(mean(x_1), mean(x_2), mean(x_3), mean(x_4) ), 2))
    
    lines(density( x_1,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
    lines(density( x_2,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
    lines(density( x_3,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
    lines(density( x_4,  n = 10^2, adjust = 1), col = cols[i], lty = i, lwd = 2 )
  }
}

################################################################################################################
# Posterior densities for the E(X) and Var(X)/E(X)
################################################################################################################
par(mfrow = c(1,2))

# E(X)
plot(NA,
     main = expression(mu[X]^{(k)}), xlab = "", ylab = "", cex.main = 1.5, cex.axis = 1.5, cex.lab = 0.8,
     xlim = c(0, 2000), ylim = c(0, 0.0065))

for(i in 1:7){
  if(i %in% SELECT){
    a_0 = exp(hierarchical_1[[i]][,1]);  a_1 = exp(hierarchical_1[[i]][,2])
    kON = exp(hierarchical_1[[i]][,3]); kOFF = exp(hierarchical_1[[i]][,4])
    mu_X = (a_1 * kON + a_0 * kOFF)/(kON + kOFF)
    var_X = mu_X + { ((a_1 - a_0)^2) * kON * kOFF}/{ ((kON + kOFF)^2) * (1 + kON + kOFF) }
    lines(density(mu_X,  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
    
    a_0 = exp(hierarchical_2[[i]][,1]);  a_1 = exp(hierarchical_2[[i]][,2])
    kON = exp(hierarchical_2[[i]][,3]); kOFF = exp(hierarchical_2[[i]][,4])
    mu_X = (a_1 * kON + a_0 * kOFF)/(kON + kOFF)
    var_X = mu_X + { ((a_1 - a_0)^2) * kON * kOFF}/{ ((kON + kOFF)^2) * (1 + kON + kOFF) }
    lines(density(mu_X,  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
    
    a_0 = exp(hierarchical_3[[i]][,1]);  a_1 = exp(hierarchical_3[[i]][,2])
    kON = exp(hierarchical_3[[i]][,3]); kOFF = exp(hierarchical_3[[i]][,4])
    mu_X = (a_1 * kON + a_0 * kOFF)/(kON + kOFF)
    var_X = mu_X + { ((a_1 - a_0)^2) * kON * kOFF}/{ ((kON + kOFF)^2) * (1 + kON + kOFF) }
    lines(density(mu_X,  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
    
    a_0 = exp(hierarchical_4[[i]][,1]);  a_1 = exp(hierarchical_4[[i]][,2])
    kON = exp(hierarchical_4[[i]][,3]); kOFF = exp(hierarchical_4[[i]][,4])
    mu_X = (a_1 * kON + a_0 * kOFF)/(kON + kOFF)
    var_X = mu_X + { ((a_1 - a_0)^2) * kON * kOFF}/{ ((kON + kOFF)^2) * (1 + kON + kOFF) }
    lines(density(mu_X,  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  }
}

# Var(X)/E(X)
plot(NA,
     main = expression(sigma[X]^{2 (i)}/mu[X]^{(k)}), xlab = "", ylab = "", cex.main = 1.5, cex.axis = 1.5, cex.lab = 0.8,
     xlim = c(0, 2500), ylim = c(0, 0.0045))

for(i in 1:7){
  if(i %in% SELECT){
    a_0 = exp(hierarchical_1[[i]][,1]);  a_1 = exp(hierarchical_1[[i]][,2])
    kON = exp(hierarchical_1[[i]][,3]); kOFF = exp(hierarchical_1[[i]][,4])
    mu_X = (a_1 * kON + a_0 * kOFF)/(kON + kOFF)
    var_X = mu_X + { ((a_1 - a_0)^2) * kON * kOFF}/{ ((kON + kOFF)^2) * (1 + kON + kOFF) }
    lines(density(var_X/mu_X,  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
    
    a_0 = exp(hierarchical_2[[i]][,1]);  a_1 = exp(hierarchical_2[[i]][,2])
    kON = exp(hierarchical_2[[i]][,3]); kOFF = exp(hierarchical_2[[i]][,4])
    mu_X = (a_1 * kON + a_0 * kOFF)/(kON + kOFF)
    var_X = mu_X + { ((a_1 - a_0)^2) * kON * kOFF}/{ ((kON + kOFF)^2) * (1 + kON + kOFF) }
    lines(density(var_X/mu_X,  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
    
    a_0 = exp(hierarchical_3[[i]][,1]);  a_1 = exp(hierarchical_3[[i]][,2])
    kON = exp(hierarchical_3[[i]][,3]); kOFF = exp(hierarchical_3[[i]][,4])
    mu_X = (a_1 * kON + a_0 * kOFF)/(kON + kOFF)
    var_X = mu_X + { ((a_1 - a_0)^2) * kON * kOFF}/{ ((kON + kOFF)^2) * (1 + kON + kOFF) }
    lines(density(var_X/mu_X,  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
    
    a_0 = exp(hierarchical_4[[i]][,1]);  a_1 = exp(hierarchical_4[[i]][,2])
    kON = exp(hierarchical_4[[i]][,3]); kOFF = exp(hierarchical_4[[i]][,4])
    mu_X = (a_1 * kON + a_0 * kOFF)/(kON + kOFF)
    var_X = mu_X + { ((a_1 - a_0)^2) * kON * kOFF}/{ ((kON + kOFF)^2) * (1 + kON + kOFF) }
    lines(density(var_X/mu_X,  n = 10^3, adjust = 1), col = cols[i], lty = i, lwd = 1.5 )
  }
}

################################################################################################################
# Posterior densities for the E(Y) and SD(Y)
################################################################################################################
par(mfrow = c(1,2))

# E(Y)
plot(NA,
     main = expression(mu[Y]^{(k)}), xlab = "", ylab = "", cex.main = 1.5, cex.axis = 1.5, cex.lab = 0.8,
     xlim = c(1500, 4000), ylim = c(0, 0.012))

for(i in 1:7){
  if(i %in% SELECT){
    a_0 = exp(hierarchical_1[[i]][,1]);  a_1 = exp(hierarchical_1[[i]][,2])
    kON = exp(hierarchical_1[[i]][,3]); kOFF = exp(hierarchical_1[[i]][,4])
    kappa = exp(hierarchical_1[[i]][,5])
    mu_epsilon = exp(hierarchical_1[[i]][,6])
    mu_X = (a_1 * kON + a_0 * kOFF)/(kON + kOFF)
    x = kappa * mu_X + mu_epsilon
    lines(density(x,  n = 10^3, adjust = 1), col = cols[i],  lty = i, lwd = 1.5 )
    abline(v = E_Y[i,1],lty = 2, col = cols[i])
    
    a_0 = exp(hierarchical_2[[i]][,1]);    a_1 = exp(hierarchical_2[[i]][,2])
    kON = exp(hierarchical_2[[i]][,3]);    kOFF = exp(hierarchical_2[[i]][,4])
    kappa = exp(hierarchical_2[[i]][,5]);  mu_epsilon = exp(hierarchical_2[[i]][,6])
    mu_X = (a_1 * kON + a_0 * kOFF)/(kON + kOFF)
    x = kappa * mu_X + mu_epsilon
    lines(density(x,  n = 10^3, adjust = 1), col = cols[i],  lty = i, lwd = 1.5 )
    abline(v = E_Y[i,2],lty = 2, col = cols[i])
    
    a_0 = exp(hierarchical_3[[i]][,1]);  a_1 = exp(hierarchical_3[[i]][,2])
    kON = exp(hierarchical_3[[i]][,3]); kOFF = exp(hierarchical_3[[i]][,4])
    kappa = exp(hierarchical_3[[i]][,5])
    mu_epsilon = exp(hierarchical_3[[i]][,6])
    mu_X = (a_1 * kON + a_0 * kOFF)/(kON + kOFF)
    x = kappa * mu_X + mu_epsilon
    lines(density(x,  n = 10^3, adjust = 1), col = cols[i],  lty = i, lwd = 1.5 )
    abline(v = E_Y[i,3],lty = 2, col = cols[i])
    
    a_0 = exp(hierarchical_4[[i]][,1]);  a_1 = exp(hierarchical_4[[i]][,2])
    kON = exp(hierarchical_4[[i]][,3]); kOFF = exp(hierarchical_4[[i]][,4])
    kappa = exp(hierarchical_4[[i]][,5])
    mu_epsilon = exp(hierarchical_4[[i]][,6])
    mu_X = (a_1 * kON + a_0 * kOFF)/(kON + kOFF)
    x = kappa * mu_X + mu_epsilon
    lines(density(x,  n = 10^3, adjust = 1), col = cols[i],  lty = i, lwd = 1.5 )
    abline(v = E_Y[i,4],lty = 2, col = cols[i])
  }
}

# SD(Y)
plot(NA,
     main = expression(sigma[Y]^{(k)}), xlab = "", ylab = "", cex.main = 1.5, cex.axis = 1.5, cex.lab = 0.8,
     xlim = c(800, 3300), ylim = c(0, 0.009))

for(i in 1:7){
  if(i %in% SELECT){
    a_0 = exp(hierarchical_1[[i]][,1]);  a_1 = exp(hierarchical_1[[i]][,2])
    kON = exp(hierarchical_1[[i]][,3]); kOFF = exp(hierarchical_1[[i]][,4])
    kappa = exp(hierarchical_1[[i]][,5])
    mu_epsilon = exp(hierarchical_1[[i]][,6])
    sigma_epsilon = exp(hierarchical_1[[i]][,7])
    mu_X = (a_1 * kON + a_0 * kOFF)/(kON + kOFF)
    var_X = mu_X + { ((a_1 - a_0)^2) * kON * kOFF}/{ ((kON + kOFF)^2) * (1 + kON + kOFF) }
    x = sqrt(  kappa^2 * var_X + sigma_epsilon^2 )
    lines(density(x,  n = 10^3, adjust = 1), col = cols[i],  lty = i, lwd = 1.5 )
    abline(v = SD_Y[i,1],lty = 2, col = cols[i])
    
    
    a_0 = exp(hierarchical_2[[i]][,1]);  a_1 = exp(hierarchical_2[[i]][,2])
    kON = exp(hierarchical_2[[i]][,3]); kOFF = exp(hierarchical_2[[i]][,4])
    kappa = exp(hierarchical_2[[i]][,5])
    mu_epsilon = exp(hierarchical_2[[i]][,6])
    sigma_epsilon = exp(hierarchical_2[[i]][,7])
    mu_X = (a_1 * kON + a_0 * kOFF)/(kON + kOFF)
    var_X = mu_X + { ((a_1 - a_0)^2) * kON * kOFF}/{ ((kON + kOFF)^2) * (1 + kON + kOFF) }
    x = sqrt(  kappa^2 * var_X + sigma_epsilon^2 )
    lines(density(x,  n = 10^3, adjust = 1), col = cols[i],  lty = i, lwd = 1.5 )
    abline(v = SD_Y[i,2],lty = 2, col =  cols[i])
    
    a_0 = exp(hierarchical_3[[i]][,1]);  a_1 = exp(hierarchical_3[[i]][,2])
    kON = exp(hierarchical_3[[i]][,3]); kOFF = exp(hierarchical_3[[i]][,4])
    kappa = exp(hierarchical_3[[i]][,5])
    mu_epsilon = exp(hierarchical_3[[i]][,6])
    sigma_epsilon = exp(hierarchical_3[[i]][,7])
    mu_X = (a_1 * kON + a_0 * kOFF)/(kON + kOFF)
    var_X = mu_X + { ((a_1 - a_0)^2) * kON * kOFF}/{ ((kON + kOFF)^2) * (1 + kON + kOFF) }
    x = sqrt(  kappa^2 * var_X + sigma_epsilon^2 )
    lines(density(x,  n = 10^3, adjust = 1), col = cols[i],  lty = i, lwd = 1.5 )
    abline(v = SD_Y[i,3],lty = 2, col = cols[i])
    
    
    a_0 = exp(hierarchical_4[[i]][,1]);  a_1 = exp(hierarchical_4[[i]][,2])
    kON = exp(hierarchical_4[[i]][,3]); kOFF = exp(hierarchical_4[[i]][,4])
    kappa = exp(hierarchical_4[[i]][,5])
    mu_epsilon = exp(hierarchical_4[[i]][,6])
    sigma_epsilon = exp(hierarchical_4[[i]][,7])
    mu_X = (a_1 * kON + a_0 * kOFF)/(kON + kOFF)
    var_X = mu_X + { ((a_1 - a_0)^2) * kON * kOFF}/{ ((kON + kOFF)^2) * (1 + kON + kOFF) }
    x = sqrt(  kappa^2 * var_X + sigma_epsilon^2 )
    lines(density(x,  n = 10^3, adjust = 1), col = cols[i],  lty = i, lwd = 1.5 )
    abline(v = SD_Y[i,4],lty = 2, col = cols[i])
  }
}

################################################################################################################
# Traceplot of the the posterior chains for E(Y) and SD(Y), with lines from the observed data.
################################################################################################################
par(mfrow = c(4,2) )

for(k in 1:2){
  a0.m   = exp(hierarchical_1[[k]][,1]);  a1.m   = exp(hierarchical_1[[k]][,2])
  kON.m  = exp(hierarchical_1[[k]][,3]);  kOFF.m = exp(hierarchical_1[[k]][,4])
  kappa.m = exp(hierarchical_1[[k]][,5]);  mu_epsilon.m = exp(hierarchical_1[[k]][,6])
  sigma_epsilon.m = exp(hierarchical_1[[k]][,7])
  
  mu_X.m = (a1.m * kON.m + a0.m * kOFF.m)/(kON.m + kOFF.m)
  var_X.m = mu_X.m + { ((a1.m - a0.m)^2) * kON.m * kOFF.m}/{ ((kON.m + kOFF.m)^2) * (1 + kON.m + kOFF.m) }
  mu_LN.m = mu_epsilon.m; var_LN.m = sigma_epsilon.m^2
  
  plot( kappa.m * mu_X.m + mu_LN.m, type = "l", main = expression(mu[Y]),
        ylab = "", xlab = "Iteration", cex.axis = 1, cex.main = 1.2, cex.lab = 1)
  abline(h = E_Y[k,1], col = "green", lwd = 2)
  plot( sqrt( kappa.m^2 * var_X.m + var_LN.m), type = "l",  main = expression(sigma[Z]),
        ylab = "", xlab = "Iteration", cex.axis = 1, cex.main = 1.2, cex.lab = 1)
  abline(h = SD_Y[k,1], col = "green", lwd = 2)
  
  a0.m   = exp(hierarchical_2[[k]][,1]);  a1.m   = exp(hierarchical_2[[k]][,2])
  kON.m  = exp(hierarchical_2[[k]][,3]);  kOFF.m = exp(hierarchical_2[[k]][,4])
  kappa.m = exp(hierarchical_2[[k]][,5]);  mu_epsilon.m = exp(hierarchical_2[[k]][,6])
  sigma_epsilon.m = exp(hierarchical_2[[k]][,7])
  
  mu_X.m = (a1.m * kON.m + a0.m * kOFF.m)/(kON.m + kOFF.m)
  var_X.m = mu_X.m + { ((a1.m - a0.m)^2) * kON.m * kOFF.m}/{ ((kON.m + kOFF.m)^2) * (1 + kON.m + kOFF.m) }
  mu_LN.m = mu_epsilon.m; var_LN.m = sigma_epsilon.m^2
  plot( kappa.m * mu_X.m + mu_LN.m, type = "l", main = expression(mu[Y]),
        ylab = "", xlab = "Iteration", cex.axis = 1, cex.main = 1.2, cex.lab = 1)
  abline(h = E_Y[k,2], col = "green", lwd = 2)
  plot( sqrt( kappa.m^2 * var_X.m + var_LN.m), type = "l",  main = expression(sigma[Z]),
        ylab = "", xlab = "Iteration", cex.axis = 1, cex.main = 1.2, cex.lab = 1)
  abline(h = SD_Y[k,2], col = "green", lwd = 2)
  
  
  a0.m   = exp(hierarchical_3[[k]][,1]);  a1.m   = exp(hierarchical_3[[k]][,2])
  kON.m  = exp(hierarchical_3[[k]][,3]);  kOFF.m = exp(hierarchical_3[[k]][,4])
  kappa.m = exp(hierarchical_3[[k]][,5]);  mu_epsilon.m = exp(hierarchical_3[[k]][,6])
  sigma_epsilon.m = exp(hierarchical_3[[k]][,7])
  
  
  mu_X.m = (a1.m * kON.m + a0.m * kOFF.m)/(kON.m + kOFF.m)
  var_X.m = mu_X.m + { ((a1.m - a0.m)^2) * kON.m * kOFF.m}/{ ((kON.m + kOFF.m)^2) * (1 + kON.m + kOFF.m) }
  mu_LN.m = mu_epsilon.m; var_LN.m = sigma_epsilon.m^2
  plot( kappa.m * mu_X.m + mu_LN.m, type = "l", main = expression(mu[Y]),
        ylab = "", xlab = "Iteration", cex.axis = 1, cex.main = 1.2, cex.lab = 1)
  abline(h = E_Y[k,3], col = "green", lwd = 2)
  plot( sqrt( kappa.m^2 * var_X.m + var_LN.m), type = "l",  main = expression(sigma[Z]),
        ylab = "", xlab = "Iteration", cex.axis = 1, cex.main = 1.2, cex.lab = 1)
  abline(h = SD_Y[k,3], col = "green", lwd = 2)
  
  
  a0.m   = exp(hierarchical_4[[k]][,1]);  a1.m   = exp(hierarchical_4[[k]][,2])
  kON.m  = exp(hierarchical_4[[k]][,3]);  kOFF.m = exp(hierarchical_4[[k]][,4])
  kappa.m = exp(hierarchical_4[[k]][,5]);  mu_epsilon.m = exp(hierarchical_4[[k]][,6])
  sigma_epsilon.m = exp(hierarchical_4[[k]][,7])
  

  mu_X.m = (a1.m * kON.m + a0.m * kOFF.m)/(kON.m + kOFF.m)
  var_X.m = mu_X.m + { ((a1.m - a0.m)^2) * kON.m * kOFF.m}/{ ((kON.m + kOFF.m)^2) * (1 + kON.m + kOFF.m) }
  mu_LN.m = mu_epsilon.m; var_LN.m = sigma_epsilon.m^2
  plot( kappa.m * mu_X.m + mu_LN.m, type = "l", main = expression(mu[Y]),
        ylab = "", xlab = "Iteration", cex.axis = 1, cex.main = 1.2, cex.lab = 1)
  abline(h = E_Y[k,4], col = "green", lwd = 2)
  plot( sqrt( kappa.m^2 * var_X.m + var_LN.m), type = "l",  main = expression(sigma[Z]),
        ylab = "", xlab = "Iteration", cex.axis = 1, cex.main = 1.2, cex.lab = 1)
  abline(h = SD_Y[k,4], col = "green", lwd = 2)
}

