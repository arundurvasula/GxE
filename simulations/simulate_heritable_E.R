library(dplyr)
sim = function(){
  M=1e3
  Nho=49e3
  N_training=337e3
  h2 = 0.25
  B_causal = rnorm(M, 0, sqrt(h2/M))
  phenotypes = rep(NA, Nho)
  genotype = as.matrix(replicate(M, rbinom(Nho, 2, 0.5)), nrow=M)
  for (i in 1:Nho){
    ind_geno = genotype[i,]
    phenotypes[i] = ind_geno%*%B_causal + rnorm(1, 0, sqrt(1-h2))
  }
  
  B_hat = rnorm(M, B_causal, sd = sqrt((1-h2/M )/ N_training))
  h2g = (M/N_training) * mean(B_hat^2 * N_training -1)
  B_prs = (h2g / (h2g+(M/N_training)))*B_hat 
  
  PRS = genotype%*%B_prs
  
  x = data.frame(phenotype = phenotypes, PRS=PRS)
  x2 = x %>% mutate(phenotype_bin = cut(phenotype, breaks=5, labels=F))
  bin1_model = (x2 %>% filter(phenotype_bin==1) %>% lm(phenotype ~ PRS, data=.) %>% summary())$r.squared
  bin2_model = (x2 %>% filter(phenotype_bin==2) %>% lm(phenotype ~ PRS, data=.) %>% summary())$r.squared
  bin3_model = (x2 %>% filter(phenotype_bin==3) %>% lm(phenotype ~ PRS, data=.) %>% summary())$r.squared
  bin4_model = (x2 %>% filter(phenotype_bin==4) %>% lm(phenotype ~ PRS, data=.) %>% summary())$r.squared
  bin5_model = (x2 %>% filter(phenotype_bin==5) %>% lm(phenotype ~ PRS, data=.) %>% summary())$r.squared
  full_model = summary(lm(phenotype ~ PRS*phenotype_bin, x2))$coefficients["PRS:phenotype_bin","t value"]
  return(c(bin1_model, bin2_model, bin3_model, bin4_model, bin5_model, full_model ))
}
results_list <- list()
for (i in 1:100){
  res = sim()
  results_list[[i]] <- res
}
results_df <- do.call(rbind, lapply(results_list, function(x) as.data.frame(t(x))))
colnames(results_df) <- c("Bin 1", "Bin 2", "Bin 3", "Bin 4", "Bin 5", "Full")
boxplot(results_df[,1:5], ylab="PRS R^2")
pdf("~/Documents/pricelab/ps_gxe/simulations/20240227_heritableE/qqplot_heritableE.pdf", 5, 5)
par(cex=1.3)
qqnorm(results_df$Full)
qqline(results_df$Full)
dev.off()