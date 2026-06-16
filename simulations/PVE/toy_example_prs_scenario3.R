get_est = function(xi){
  n=49000
  m=1000
  h2=0.25
  beta = rnorm(m, 0, sqrt(h2/m))
  E = rnorm(n, 0, 1)
  eta = xi
  epsilon = rnorm(n, 0, sqrt(1-h2))
  genotype = scale(as.matrix(replicate(m, rbinom(n, 2, 0.5)), nrow=m))
  y = genotype %*% beta + xi *genotype %*% beta * E + epsilon + eta * epsilon * E
  beta_hat = rep(NA, m)
  for (i in 1:m){
    reg = lm(y ~ genotype[,i])
    beta_hat[i] = reg$coefficients[2]
  }
  beta_prs = (h2 / (h2+(m/n)))*beta_hat
  prs = genotype %*% beta_prs
  base = anova(lm (y ~ prs + E))
  int = anova(lm(y ~ prs + E + prs*E))
  # need to compute interaction variance explained
  VE =  ( (sum(int[1:3,2])/sum(int[,2])) - (sum(base[1:2,2])/sum(base[,2])) )
  basePRSVE = (sum(base[1,2])/sum(base[,2]))
  VE / basePRSVE / h2
  # then, we can compare to the expected amount of variance explained 
  g = genotype %*% beta
  true_base = anova(lm (y ~ g + E))
  true_int = anova(lm(y ~ g + E + g*E))
  true_VE =  ( (sum(true_int[1:3,2])/sum(true_int[,2])) - (sum(true_base[1:2,2])/sum(true_base[,2])) )
  truebasePRSVE = (sum(true_base[1,2])/sum(true_base[,2]))
  estimatedparam = VE / basePRSVE
  trueparam = true_VE / truebasePRSVE
  return(c(estimatedparam, trueparam))
}
xis = c(0.025, 0.03, 0.05, 0.1, 0.2, 0.3, 0.4)
outs = sapply(xis, get_est)

plot(xis, diff(outs))
setwd("~/Documents/pricelab/ps_gxe/simulations/20230830_PVE/plots/")
pdf("Scenario3.pdf", 4, 4)
par(cex=1.05, mar=c(5,4,2,2)+0.1)
plot(outs[1,], outs[2,], ylab="True GxE Variance", xlab="Estimated GxE Variance", pch=19,
     ylim=c(1e-4, 0.2), xlim=c(1e-4, 0.2))
abline(0,1)
dev.off()


summary(lm(outs[1,]~outs[2,]))

