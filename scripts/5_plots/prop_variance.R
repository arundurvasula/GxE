setwd("~/Documents/pricelab/ps_gxe/2023_05_10_tables/")
library(tidyverse)
library(viridis)
library(qvalue)
library(data.table)
library(dplyr)

inddata = data.frame(fread("../newE/337K-pheno-cov-newE-dietpc.tab"))
combined = data.frame(fread("parsed/combined.csv", header=T))
traits = unique(combined$trait)
traits = traits[!grepl('disease', traits)] # remove binary traits
Es = unique(combined$E)

incremental_R2 <- function(data, trait, E) {
  with_E = summary(lm(data[, trait] ~ data[, E] + 
                        data[, "cov_AGE"]*data[, "cov_SEX"] + 
                        data[, "PC1"] + data[, "PC2"] + data[, "PC3"] + data[, "PC4"] + data[, "PC5"] + 
                        data[, "PC6"] + data[, "PC7"] + data[, "PC8"] + data[, "PC9"] + data[, "PC10"] ))$r.squared
  without_E = summary(lm(data[, trait] ~ 
                        data[, "cov_AGE"]*data[, "cov_SEX"] + 
                        data[, "PC1"] + data[, "PC2"] + data[, "PC3"] + data[, "PC4"] + data[, "PC5"] + 
                        data[, "PC6"] + data[, "PC7"] + data[, "PC8"] + data[, "PC9"] + data[, "PC10"] ))$r.squared
  return(with_E - without_E)
}
incremental_R2_sex <- function(data, trait){
  with_sex = summary(lm(data[, trait] ~
                        data[, "cov_AGE"] + data[, "cov_SEX"] + 
                        data[, "PC1"] + data[, "PC2"] + data[, "PC3"] + data[, "PC4"] + data[, "PC5"] + 
                        data[, "PC6"] + data[, "PC7"] + data[, "PC8"] + data[, "PC9"] + data[, "PC10"] ))$r.squared
  without_sex = summary(lm(data[, trait] ~ 
                           data[, "cov_AGE"] + 
                           data[, "PC1"] + data[, "PC2"] + data[, "PC3"] + data[, "PC4"] + data[, "PC5"] + 
                           data[, "PC6"] + data[, "PC7"] + data[, "PC8"] + data[, "PC9"] + data[, "PC10"] ))$r.squared
  return(with_sex - without_sex)
}
trait_E_var_explained = rep(NA, length(traits))
trait_sex_var_explained = rep(NA, length(traits))
for (i in 1:length(traits)){
  curr_Es = rep(NA, length(Es))
  for (j in 1:length(Es)){
      curr_Es[j] = incremental_R2(inddata, traits[i], Es[j])
  }
  trait_sex_var_explained[i] = incremental_R2_sex(inddata, traits[i])
  trait_E_var_explained[i] = sum(curr_Es)
}
summary.table = data.frame(cbind(traits, trait_E_var_explained, trait_sex_var_explained))
names(summary.table) = c("trait", "trait_E_var_explained", "trait_sex_var_explained")

Figure3 = read.table("tables/Figure3_prs_only.csv", header=T, sep=",")
GxE = Figure3 %>% group_by(trait) %>% replace(is.na(.), 0) %>% summarize(GxE.VE = sum(prs.plot.diff/100, na.rm=T) + sum((1-rg.plot.diff)/10))
GxE$trait = gsub(" ", "_", GxE$trait)
Figure5 = read.table("tables/Figure5_prs_only.csv", header=T, sep=",")
GxSex = Figure5 %>% group_by(trait) %>% replace(is.na(.), 0) %>% summarize(GxSex.VE = prs.plot.diff/100 + (1-rg.plot.diff)/2)
GxSex$trait = gsub(" ", "_", GxSex$trait)

df_list <- list(summary.table, GxE, GxSex)

#merge all data frames in list
out = df_list %>% reduce(left_join, by='trait')
out <- out %>% mutate_all(~ifelse(is.nan(.), NA, .)) # replace nan with NA
out$trait_E_var_explained = as.numeric(out$trait_E_var_explained)
out$trait_sex_var_explained = as.numeric(out$trait_sex_var_explained)
write.table(out, "tables/prop_variance_v2_revision.csv", sep=",", row.names=F)
mean(out$GxE.VE, na.rm=T)/mean(out$trait_E_var_explained, na.rm=T)
mean(out$GxSex.VE, na.rm=T)/mean(out$trait_sex_var_explained, na.rm=T)

mean(out$trait_E_var_explained, na.rm=T)
mean(out$trait_sex_var_explained, na.rm=T)
mean(out$GxE.VE, na.rm=T)
mean(out$GxSex.VE, na.rm=T)
