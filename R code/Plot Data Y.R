load("Data Tet5.RData") # Data for cells stimulated with 5 ng/ML of tetracycline.
Y_5 = Y

load("Data Tet10.RData") # Data for cells stimulated with 10 ng/ML of tetracycline.
Y_10 = Y

# Plot the densities for the experimental data for each replicate, under both experimental conditions:
par(mfrow = c(1,1))
plot(NULL, 
     main = expression(Y^{(k)}),
     #main = expression(e^{mu[1]}),
     xlim = c(0,1.5*10^4), ylim = c(0, 0.0007),  
     xlab = "", ylab = "", cex.main = 2, cex.axis = 2, cex.lab = 1)
for(k in 1:4){
  lines(density(Y_5[[k]],  n = 10^3, adjust = 2), col = 1, lty = 1, lwd = 2 )
  lines(density(Y_10[[k]], n = 10^3, adjust = 2), col = 2, lty = 2, lwd = 2 )
}
