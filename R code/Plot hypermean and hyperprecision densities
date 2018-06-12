################################################################################################################
# Load the posterior output:
################################################################################################################
HYPER_MODE = list()

# load the posterior output of the Tet_5 condition
load(...)
S.   = burn_in/thin + 1
end. = R/thin # total number of iterations of the MCMC

HYPER_MODE[[1]] = mcmc.hyper[[1]][S.:end.,]
HYPER_TAU[[1]]  = mcmc.hyper[[2]][S.:end.,]

# load the posterior output of the Tet_10 condition
load(...)
S.   = burn_in/thin + 1
end. = R/thin # total number of iterations of the MCMC

HYPER_MODE[[2]] = mcmc.hyper[[1]][S.:end.,]
HYPER_TAU[[2]]  =  mcmc.hyper[[2]][S.:end.,]

################################################################################################################
# Posterior densities for the exponential of the hypermean parameters:
################################################################################################################
par(mfrow = c(3,2))
cols = 1:2

# alpha_0
plot(NULL, 
     main = expression(tilde(alpha)[0]),
     xlim = c(-5,200), ylim = c(0, 0.06),  
     xlab = "", ylab = "", cex.main = 2, cex.axis = 2, cex.lab = 1)
for(i in 1:2){
  lines(density(exp( HYPER_MODE[[i]][,1]), to = 10^4, n = 10^4, adjust = 1),
        col = cols[i], lty = i, lwd = 2 )
}

# alpha_1
plot(NULL, 
     main = expression(tilde(alpha)[1]),
     xlim = c(-1000,2*10^4), ylim = c(0, 0.0003),
     xlab = "", ylab = "", cex.main = 2, cex.axis = 2, cex.lab = 1)
for(i in 1:2){
  lines(density(exp(HYPER_MODE[[i]][,2]), n = 10^4, adjust = 2), col = cols[i],  lty = i, lwd = 2 )
}

# k_1
plot(NULL, 
     main = expression(tilde(k)[1]),
     xlim = c(-0.1, 1), ylim = c(0,8),
     xlab = "", ylab = "", cex.main = 2, cex.axis = 2, cex.lab = 1)
for(i in 1:2){
  lines(density(exp(HYPER_MODE[[i]][,3]), n = 10^3, adjust = 2), col = cols[i],  lty = i, lwd = 2 )
}  

# k_0
plot(NULL, 
     main = expression(tilde(k)[0]),
     xlim = c(1.5, 12), ylim = c(0, 0.8),
     xlab = "", ylab = "", cex.main = 2, cex.axis = 2, cex.lab = 1)

for(i in 1:2){
  lines(density(exp(HYPER_MODE[[i]][,4]), n =2*10^4, adjust = 2), col = cols[i],  lty = i, lwd = 2 )
}

# kappa
plot(NULL,  main = expression(kappa),
     xlim = c(-10, 100), ylim = c(0, 0.15),
     xlab = "", ylab = "", cex.main = 2, cex.axis = 2, cex.lab = 1)

for(i in 1:2){
  lines(density(exp(HYPER_MODE[[i]][,5]), n =2*10^4, adjust = 1), col = cols[i],  lty = i, lwd = 2 )
}

################################################################################################################
# Posterior densities for the hyper-precision parameter tau:
################################################################################################################
par(mfrow = c(3,2))

# alpha_0
plot(NULL, main = expression(tau[1]),
     ylim = c(0,1.5), xlim = c(0, 6),
     xlab = "", ylab = "", cex.main = 2, cex.axis = 1.5, cex.lab = 1)
for(i in 1:2){
  lines(density((HYPER_TAU[[i]][,1]),  n = 10^4), col = cols[i], lty = i, lwd = 2 )
}

# alpha_1
plot(NULL, main = expression(tau[2]),
     xlim = c(0,10^3), ylim = c(0, 0.012),     
     xlab = "", ylab = "", cex.main = 2, cex.axis = 1.5, cex.lab = 1)
for(i in 1:2){
  lines(density((HYPER_TAU[[i]][,2]), n = 10^4), col = cols[i], lty = i, lwd = 2 )
}

# k_1
plot(NULL, main = expression(tau[3]),
     xlim = c(0,100), ylim = c(0, 0.08),     
     xlab = "", ylab = "", cex.main = 2, cex.axis = 1.5, cex.lab = 1)
for(i in 1:2){  
  lines(density((HYPER_TAU[[i]][,3]), n = 10^4), col = cols[i], lty = i, lwd = 2 )
}

# k_0
plot(NULL, main = expression(tau[4]),
     xlim = c(-5, 200), ylim  = c(0,0.03),     
     xlab = "", ylab = "", cex.main = 2, cex.axis = 1.5, cex.lab = 1)
for(i in 1:2){
  lines(density((HYPER_TAU[[i]][,4]), n = 10^4), col = cols[i], lty = i, lwd = 2 )
}

# kappa
plot(NULL, main = expression(tau[5]),
     xlim = c(-10, 300), ylim = c(0, 0.018),     
     xlab = "", ylab = "", cex.main = 2, cex.axis = 1.5, cex.lab = 1)
for(i in 1:2){
  lines(density((HYPER_TAU[[i]][,5]), n = 10^4), col = cols[i], lty = i, lwd = 2 )
}
