load("background_error_data.RData")
N = 4 # 4 replicates.

####################################################################################
# Histogram of the backgound error data
####################################################################################
par(mfrow = c(2,2) )
for(k in 1:N){
  hist( Control[[k]], main = paste("Noise in replicate ", k),
        xlab = "", cex.axis = 0.8, cex.main = 1, cex.lab = 0.8,
        ylim = c(0, 0.0016),
        xlim = c(0, 2500), breaks = seq(0,2500,50), freq = FALSE)
  print(length(Control[[k]]))
}

####################################################################################
# Posterior densities for the logarithm of the hierarchcial parameters:
####################################################################################
par(mfrow = c(1,2) )
S. = 2*10^5 + 1
end. = R
plot( NULL, 
      type = "l",  ylab = "",  xlab = "", 
      xlim = c(6.63,6.89), ylim = c(0, 115),
      cex.axis =0.8, cex.main = 1.2, cex.lab = 0.8,
      main = expression(log(mu[epsilon]^{(k)})))
for(k in 1:4){
  lines( density( mcmc[[k]][seq(S.,end.,1),1],  adjust = 2),
         col = k)
}

plot( NULL, 
      type = "l",  ylab = "",  xlab = "", 
      xlim = c(5.69,6.01), ylim = c(0, 60),
      cex.axis =0.8, cex.main = 1.2, cex.lab = 0.8,      
      main = expression(log(sigma[epsilon]^{(k)})))
for(k in 1:4){
  lines( density( mcmc[[k]][seq(S.,end.,1),2],  adjust = 2),
         col = k)
}

####################################################################################
# Posterior meand and standard deviation of the logarithm of the hierarchical parameters:
####################################################################################
S. = 2*10^5 + 1
end. = R

MEAN = SD = matrix(NA, nrow = N, ncol = 2)
for(k in 1:N){
  MEAN[k,1] = mean( ( mcmc[[k]][seq(S.,end.,1),1] ))
  MEAN[k,2] = mean( ( mcmc[[k]][seq(S.,end.,1),2] ) )
  SD[k,1] = sd( ( mcmc[[k]][seq(S.,end.,1),1] ))
  SD[k,2] = sd( ( mcmc[[k]][seq(S.,end.,1),2] ) )
}

# save(MEAN, SD, file = "prior on error.R")
library(xtable)
xtable( t(cbind(MEAN, 10^3 * SD)), digits = 2)
