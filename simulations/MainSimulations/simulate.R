setwd("/n/groups/price/arun/ps_gxe/simulations_20230623")
library(MASS)
library(tidyverse)
M=1e4
N_training = 337e3
Nho = 49e3
N_bin = 67e3

rg_h2 = function(B_hat_1, B_hat_2){
  chisq1 = B_hat_1^2 * N_bin
  chisq2 = B_hat_2^2 * N_bin
  h21_est = (M/N_bin) * mean(chisq1 - 1)
  h22_est = (M/N_bin) * mean(chisq2 - 1)
  h2_diff = h21_est - h22_est
  rg_est = (t(B_hat_1)%*%B_hat_2)[1,1] / sqrt(h21_est * h22_est)
  h21_jk_ests = rep(NA, M)
  h22_jk_ests = rep(NA, M)
  h2_diff_ests = rep(NA, M)
  rg_jk_ests = rep(NA, M)
  for (i in 1:M){
    h21_jk_ests[i] = (M/N_bin) * mean(chisq1[-i] - 1)
    h22_jk_ests[i] = (M/N_bin) * mean(chisq2[-i] - 1)
    h2_diff_ests[i] = h21_jk_ests[i] - h22_jk_ests[i]
    rg_jk_ests[i] =(t(B_hat_1[-i])%*%B_hat_2[-i])[1,1] / sqrt(h21_jk_ests[i] * h22_jk_ests[i])
  }
  h21_se =  sd(h21_jk_ests) * sqrt(M)
  h22_se =  sd(h22_jk_ests) * sqrt(M)
  h2_diff_se = sd(h2_diff_ests) * sqrt(M) 
  rg_se = sd(rg_jk_ests) * sqrt(M)
  return(c(rg_est, rg_se, h2_diff, h2_diff_se, h21_est, h21_se, h22_est, h22_se))
}

prs = function(B_causal_1, B_causal_2, h21, h22, amp){
  B_hat_1 = rnorm(M, B_causal_1, sd = sqrt((1-h21/M )/ N_training))
  B_hat_2 = rnorm(M, B_causal_2, sd = sqrt((1-h22/M )/ N_training))
  mean_B_hat = rowMeans(cbind(B_hat_1, B_hat_2))
  h2g = (M/N_training) * mean(mean_B_hat^2 * N_training -1)
  B_prs = (h2g / (h2g+(M/N_training)))*mean_B_hat 
  
  E_bin = c(rep(1, Nho/2), rep(2, Nho/2))
  phenotypes = rep(NA, Nho)
  genotype = as.matrix(replicate(M, rbinom(Nho, 2, 0.5)), nrow=M)
  for (i in 1:Nho){
    ind_geno = genotype[i,]
    if(E_bin[i] == 1){
      phenotypes[i] = ind_geno%*%B_causal_1 + rnorm(1, 0, sqrt(1-h21))
    } else if (E_bin[i] == 2) {
      phenotypes[i] = amp * (ind_geno%*%B_causal_2 + rnorm(1, 0, sqrt(1-h22)))
    }
  }
  
  PRS = genotype%*%B_prs
  PRSxE_unscaled = summary(lm(phenotypes ~ PRS + E_bin + PRS*E_bin))
  Z_score = PRSxE_unscaled$coefficients[4,3]
  return(Z_score)  
}

run_sim = function(gen_corr, h21, h22, amp){
  rho=sqrt(h21*h22)*gen_corr
  Sigma = matrix(c(h21,rho, rho, h22),2,2)/M
  B_causal = mvrnorm(M, c(0, 0), Sigma)
  B_hat_1 = rnorm(M, B_causal[,1], sd = sqrt((1)/ N_bin)) # use this for h2 estimate (N=N_bin)
  B_hat_2 = rnorm(M, B_causal[,2], sd = sqrt((1)/ N_bin))
  res = rg_h2(B_hat_1, B_hat_2)
  rg_z = (1-res[1]) / res[2]
  h2_z = res[3] / res[4]
  prs_z = prs(B_causal[,1], B_causal[,2], h21, h22, amp)
  return(c(rg_z, h2_z, prs_z))
}

get_output = function(gen_corr, h21, h22, amp){
  NSIM=500
  pb = txtProgressBar(min = 0, max = NSIM, initial = 0, style = 3) 
  sc1 = list()
  for (i in 1:NSIM){
    x = run_sim(gen_corr=gen_corr, h21=h21, h22=h22, amp=amp)
    sc1 = append(sc1, list(x))
    setTxtProgressBar(pb,i)
  }
  close(pb)
  df1 = data.frame(t(data.frame(sc1)))
  row.names(df1) = NULL
  names(df1) = c("rg", "h2", "prs")
  rg_frac = length(df1$rg[abs(df1$rg) > 1.96 ])/NSIM
  h2_frac = length(df1$h2[abs(df1$h2) > 1.96 ])/NSIM
  prs_frac = length(df1$prs[abs(df1$prs) > 1.96 ])/NSIM
  return(c(rg_frac, h2_frac, prs_frac))
}

args = commandArgs(trailingOnly = T)
gen_corr = as.numeric(args[1])
h21 = as.numeric(args[2])
h22 = as.numeric(args[3])
amp = as.numeric(args[4])
outstem = args[5]
x = get_output(gen_corr = gen_corr, h21 = h21, h22 = h22, amp = amp)
save.image(paste0("/n/scratch3/users/a/ard063/GxE/", outstem, ".Rdata"))
write.table(x, paste0("results/",outstem, ".txt"), row.names=F, col.names=F, quote=F)
