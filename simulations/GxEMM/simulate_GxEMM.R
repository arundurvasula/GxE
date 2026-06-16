library(MASS)
library(GxEMM)
# run like Rscript scripts/simulate_GxEMM.R <seed> <h21> <h22> <gen_corr> <amplification> <# SNPS> <# SAMPLES> <outstem>
args = commandArgs(trailingOnly=TRUE)

seed=as.numeric(args[1])
h21=as.numeric(args[2])
h22=as.numeric(args[3])
gen_corr=as.numeric(args[4])
amp=as.numeric(args[5])
M=as.numeric(args[6])
N_test=as.numeric(args[7])
outstem=args[8]

set.seed(seed)
rho=sqrt(h21*h22)*gen_corr
Sigma = matrix(c(h21,rho, rho, h22),2,2)/M
B_causal = mvrnorm(M, c(0, 0), Sigma)

E_bin = c(rep(0, N_test/2), rep(1, N_test/2))
phenotypes = rep(NA, N_test)
genotype = as.matrix(replicate(M, rbinom(N_test, 2, 0.5)), nrow=M)
for (i in 1:N_test){
  ind_geno = genotype[i,]
  if(E_bin[i] == 0){
    phenotypes[i] = ind_geno%*%B_causal[,1] + rnorm(1, 0, sqrt(1-h21))
  } else if (E_bin[i] == 1) {
    phenotypes[i] = amp * (ind_geno%*%B_causal[,2] + rnorm(1, 0, sqrt(1-h22)))
  }
}

K = genotype %*% t(genotype) / M
Z = cbind(E_bin, 1-E_bin)
X     <- Z[,-1]
ldak_loc = "/n/groups/price/arun/ps_gxe/simulations_GxEMM/gxemm-master/ldak5.linux"
tmpdir=paste0("/n/scratch3/users/a/ard063/GxEMM/tmp/", outstem, "/")
out_hom	   <- GxEMM( phenotypes, X, K, Z, gtype='hom', ldak_loc=ldak_loc, tmpdir=tmpdir)
out_free   <- GxEMM( phenotypes, X, K, Z, gtype='free', etype='free', ldak_loc=ldak_loc, tmpdir = tmpdir)
out_iid	   <- GxEMM( phenotypes, X, K, Z, gtype='iid', ldak_loc=ldak_loc, tmpdir=tmpdir)

save.image(paste0("/n/scratch3/users/a/ard063/GxEMM/simulations_output/", outstem, ".Rdata"))