setwd("~/Documents/pricelab/ps_gxe/newE/")
library(data.table)
data = fread("337K-pheno-cov-newE.tab", sep="\t", na.strings = "-9")
diet = data[,.( IID,
                cov_CookedVegetableIntake,
                cov_SaladIntake,
                cov_FreshFruitIntake,
                cov_ProcessedMeatIntake,
                cov_PoultryIntake,
                cov_BeefIntake,
                cov_PorkIntake,
                cov_CoffeeIntake)]
diet_nomissing = na.omit(diet)
pc = prcomp(diet_nomissing[,-1],
            center = TRUE,
            scale. = TRUE)
diet_pc_df = data.frame(cbind(diet_nomissing$IID, pc$x[,1]))
colnames(diet_pc_df) <- c("IID", "cov_DIET")

out = merge(data[,-3], diet_pc_df, by="IID", all=TRUE) # had to remove duplicated FID column
write.table(out, file="337K-pheno-cov-newE-dietpc.tab", row.names = F, quote=F, sep="\t")

################
data = fread("49K-pheno-cov-newE.tab", sep="\t", na.strings = "-9")
diet = data[,.( IID,
                cov_CookedVegetableIntake,
                cov_SaladIntake,
                cov_FreshFruitIntake,
                cov_ProcessedMeatIntake,
                cov_PoultryIntake,
                cov_BeefIntake,
                cov_PorkIntake,
                cov_CoffeeIntake)]
diet_nomissing = na.omit(diet)
pc = prcomp(diet_nomissing[,-1],
            center = TRUE,
            scale. = TRUE)
diet_pc_df = data.frame(cbind(diet_nomissing$IID, pc$x[,1]))
colnames(diet_pc_df) <- c("IID", "cov_DIET")

out = merge(data[,-13], diet_pc_df, by="IID", all=TRUE) # had to remove duplicated FID column
write.table(out, file="49K-pheno-cov-newE-dietpc.tab", row.names = F, quote=F, sep="\t")
