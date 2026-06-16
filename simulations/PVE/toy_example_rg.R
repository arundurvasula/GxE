library(tidyverse)
setwd("~/Documents/pricelab/ps_gxe/simulations/20230830_PVE/")

get_rg = function(h2gxe, nbin=5){
  n=337000
  M=1000
  h2=0.25
  beta = rnorm(M, 0, sqrt(h2/M))
  gamma = rnorm(M, 0, sqrt(h2gxe/M))
  E = rnorm(n, 0, 1)
  epsilon = rnorm(n, 0, sqrt(1-h2-h2gxe))
  genotype = scale(as.matrix(replicate(M, rbinom(n, 2, 0.5)), nrow=M))
  y = genotype %*% beta + genotype %*% gamma * E + epsilon 
  
  df = data.frame(cbind(y, E))
  names(df) = c('y', 'E') 
  df$E_bin = cut_number(df$E, n=nbin, labels=1:nbin)
  selected_rows1 <- which(df$E_bin == 1)
  selected_rows5 <- which(df$E_bin == nbin)
  extreme_genotypes1 = genotype[selected_rows1,]
  extreme_phenotypes1 = y[selected_rows1]
  extreme_genotypes5 = genotype[selected_rows5,]
  extreme_phenotypes5 = y[selected_rows5]
  
  beta_low = rep(NA, M)
  beta_high = rep(NA, M)
  for (i in 1:M){
    lo_reg = lm(extreme_phenotypes1 ~ extreme_genotypes1[,i])
    beta_low[i] = lo_reg$coefficients[2]
    hi_reg = lm(extreme_phenotypes5 ~ extreme_genotypes5[,i])
    beta_high[i] = hi_reg$coefficients[2]
  }
  h2_low = (mean(length(selected_rows1) * beta_low ^ 2 ) - 1) * M / length(selected_rows1)
  h2_high = (mean(length(selected_rows1) * beta_high ^ 2 ) - 1) * M / length(selected_rows5)
  
  return((t(beta_low) %*% beta_high) / sqrt(h2_low * h2_high))
}
## simulate even smaller h2gxe and see what happens
h2gxes = 10^(seq(-1,-5,-.3)) 
outs = sapply(h2gxes, get_rg)
plot(h2gxes, outs)
outs2 = 1-outs
plot(h2gxes, outs2, ylim=c(0, 1))
abline(a=0.0, b=10)
summary(lm(h2gxes ~ outs2))

pdf("plots/Scenario1.pdf", 4, 4)
par(cex=1.05, mar=c(5,4,2,2)+0.1)
plot(outs2 * 0.102180,  h2gxes, ylab="True GxE Variance", xlab="Estimated GxE Variance", pch=19,
     ylim=c(0, 0.101), xlim=c(0,0.101))
abline(0,1)
dev.off()

outs3 = outs2 * 0.102180
summary(lm(outs3 ~ h2gxes))

pdf('plots/Scenario1_bias.pdf', 4, 4)
par(cex=1.05, mar=c(5,4,2,2)+0.1)
plot(h2gxes, (h2gxes - (outs2 * 0.1))/h2gxes, xlab="True GxE Variance (log scale)", ylab="Bias (True - Estimated)", log='x')
abline(h=0)
dev.off()

outs_2bin = sapply(h2gxes, get_rg, nbin=2)
plot(h2gxes,(1-outs_2bin)/4)
abline(0, 1)
outs_2bin_transformed = (1-outs_2bin)
summary(lm(h2gxes ~  outs_2bin_transformed))
summary(lm(h2gxes ~  -1 + outs_2bin_transformed))


nsim=10
rg0 = rep(NA, nsim)
for (i in 1:nsim){
  rg0[i] = get_rg(0, 5)
}
hist(rg0)

nsim=10
rg02 = rep(NA, nsim)
for (i in 1:nsim){
  rg02[i] = get_rg(0, 2)
}

get_rg(0, 5)

## there appears to be some relationship where the genetic correlation decreases the more
# SNPs there are
# 100: 0.99, 500: 0.97, 1000: 0.94 (5 bin results)
