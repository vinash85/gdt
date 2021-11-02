# Read PBMC data
# read SingleR data
# Calculate correlation 
# Calculate exterme 

library(extReme)

rhos = seq(-.9,.9, 0.2)
out=list()
z1 <- matrix(rnorm(2000), ncol = 2)
z1.var = var(z1)
ii=1
for (rho in rhos) {
rho.mat <- cbind(c(1, abs(rho)), c(abs(rho), 1))
rho.mat <- chol(rho.mat)
z <- t(rho.mat %*% t(z1))
if(rho < 0) z[,2]= -z[,2]
cor.val = cor(z[, 1], z[, 2])
# plot(z[, 1], z[, 2])
aa1 = taildep(z[, 1], z[, 2], 0.95)
aa = taildep.test(z[, 1], z[, 2])
out[[ii]] = c(rho, cor.val, aa1, aa[["p.value"]])
ii = ii +1
}