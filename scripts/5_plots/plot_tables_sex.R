setwd("~/Documents/pricelab/ps_gxe/2023_05_10_tables/")
library(tidyverse)
library(viridis)
library(qvalue)
library(data.table)
inddata = data.frame(fread("../newE/337K-pheno-cov-newE-dietpc.tab"))
get_corr = function(trait, E){
  return(cor.test(inddata[,trait], inddata[,E], na.rm=T))
}
h2 = read.table("parsed/sex_h2_values.csv", sep=",")
names(h2) = c('trait', "E", 'h2.bin2.1', 'h2.bin2.2',
              'h2.se.1', 'h2.se.2')
h2Z = read.table("parsed/sex_h2_Z-scores.csv", sep=",")
names(h2Z) = c('trait', "E", 'z.h2.2')
rg = read.table("parsed/sex_rg_values.csv", sep=",")
names(rg) = c('trait', 'E', 'rg.bin2', 'rg.se.2')
rgZ = read.table("parsed/sex_rg_Z-scores.csv", sep=",")
names(rgZ) = c('trait', "E", 'z.rg.2')
prs = read.table("parsed/sex_PRS.csv", sep=',')
names(prs) = c("trait", 'base.all', "prs.p.val", 'int.all')
prs$E = rep("cov_SEX", nrow(prs))
base_prs = read.table("parsed/PRS.csv", sep=',')
names(base_prs) = c("trait", "E", "VE.diff", "base.PS", "base.E", 'base.all', "int.PS", "int.E", "int.PRSxE", "prs.p.val", 'int.all')
base_prs_df = data.frame(base_prs %>% group_by(trait) %>% dplyr::summarise(base.PS = mean(base.PS, na.rm=T)))
prs = merge(prs, base_prs_df, by='trait')
prs$prs.plot.diff = (prs$int.all - prs$base.all) / prs$base.PS

h2_temp = merge(h2, h2Z, by=c('trait', 'E'))
rg_temp = merge(rg, rgZ, by=c('trait', 'E'))
temp_df = merge(rg_temp, h2_temp, by=c('trait', 'E'))#, by.y=c('trait', 'E'))
data = merge(temp_df, prs, by=c('trait', 'E'))

data$h2.2diff = abs(data$h2.bin2.2 - data$h2.bin2.1) / mean(c(data$h2.bin2.2, data$h2.bin2.1), na.rm=T)
data$h2.plot.diff = data$h2.2diff
data$h2.plot.diff.Z = data$z.h2.2
data$h2.plot.diff.p = pnorm(abs(data$h2.plot.diff.Z), lower.tail=F)*2

data$rg.plot.diff = data$rg.bin2
data$rg.plot.diff.Z = data$z.rg.2
data$rg.plot.diff.p = pnorm(q=data$rg.plot.diff.Z, lower.tail=FALSE)


data$mean.h2 = (data$h2.bin2.1 + data$h2.bin2.2) /2 * 100
data$prs.plot.diff = data$prs.plot.diff * data$mean.h2
######
# subset df to just p-values and plot values
######
df = data.frame(data$trait, data$E, 
                data$rg.plot.diff.p, data$h2.plot.diff.p, data$prs.p.val,
                data$rg.plot.diff, data$h2.plot.diff, data$prs.plot.diff)
names(df) = c("trait", "E", 'rg.p', 'h2.p', 'prs.p', 'rg.val', 'h2.val', 'prs.val')
write.table(df, "parsed/sex_combined.csv", sep=",", quote=F, row.names = F)
####
# Q value analysis
# Not enough p-values to estimate lambda, so set it to 0 following John Storey's 
# recommendation: https://support.bioconductor.org/p/105623/#105628
####
h2.qobj <- qvalue(data$h2.plot.diff.p, lambda=0, fdr.level = 0.05)
rg.qobj <- qvalue(data$rg.plot.diff.p, lambda=0, fdr.level = 0.05)
prs.qobj <- qvalue(data$prs.p.val, lambda=0, fdr.level = 0.05)
# hist(h2.qobj$pvalues)
# hist(rg.qobj$pvalues)
# hist(prs.qobj$pvalues)
# 
# plot(h2.qobj)
# plot(rg.qobj)
# plot(prs.qobj)
# 
# qqnorm(data$h2.plot.diff.Z)
# qqline(data$h2.plot.diff.Z)
# qqnorm(data$rg.plot.diff.Z)
# qqline(data$rg.plot.diff.Z)


#####
# filter df
######
rg_df = data.frame(data$trait, data$E, data$rg.plot.diff, rg.qobj$qvalues) %>%
  filter(rg.qobj.qvalues < 0.05)
names(rg_df) = c("trait", "E", "rg.plot.diff", 'rg.qvalues')
rg_corrs = c()
for (row in 1:nrow(rg_df)){
  x = get_corr(rg_df[row, 'trait'], rg_df[row, 'E'])
  rg_corrs = c(rg_corrs, x$estimate)
}

h2_df = data.frame(data$trait, data$E, data$h2.plot.diff, h2.qobj$qvalues) %>%
  filter(h2.qobj.qvalues < 0.05)
names(h2_df) = c("trait", "E", "h2.plot.diff", 'h2.qvalues')
prs_df = data.frame(data$trait, data$E, data$prs.plot.diff, prs.qobj$qvalues) %>%
  filter(prs.qobj.qvalues < 0.05)
names(prs_df) = c("trait", "E", "prs.plot.diff", 'prs.qvalues')

plot_df_temp = merge(prs_df, (merge(h2_df, rg_df, by=c("trait", "E"), all = T)), by=c("trait", "E"), all=T)
clean_trait = read.table("tables/plotting_names_matching_trait.csv", sep=',')
names(clean_trait) = c("trait", "trait_clean", "trait_order")
plot_df = merge(plot_df_temp, clean_trait, by="trait")

plot_df$trait = gsub("_", " ", plot_df$trait)
plot_df$E = gsub("_", " ", plot_df$E)
plot_df$E = gsub("cov ", "", plot_df$E)
plot_df$E = gsub("other ", "", plot_df$E)

plot_df$is_scenario_1 = ifelse(!is.na(plot_df$rg.qvalues), T, F)
plot_df$is_scenario_2 = ifelse((!is.na(plot_df$h2.qvalues) & !is.na((plot_df$prs.qvalues))), T, F)
plot_df$is_scenario_3 = ifelse(!is.na(plot_df$prs.qvalues) & is.na(plot_df$h2.qvalues), T, F)
plot_df$only_h2 = ifelse((!is.na(plot_df$h2.qvalues) & is.na(plot_df$prs.qvalues) & is.na(plot_df$rg.qvalues)), T, F)
table(plot_df$is_scenario_1, plot_df$is_scenario_2, plot_df$is_scenario_3)

sum(plot_df$is_scenario_1)
sum(plot_df$is_scenario_2)
sum(plot_df$is_scenario_3)
sum(plot_df$is_scenario_1 & plot_df$is_scenario_2)
sum(plot_df$is_scenario_1 & plot_df$is_scenario_3)

mean(plot_df$rg.plot.diff, na.rm=T)
sd(plot_df$rg.plot.diff, na.rm=T)
mean(plot_df$prs.plot.diff, na.rm=T)
sd(plot_df$prs.plot.diff, na.rm=T)
mean(plot_df$h2.plot.diff, na.rm=T)
sd(plot_df$h2.plot.diff, na.rm=T)


nrow(plot_df[!is.na(plot_df$rg.qvalues),])
nrow(plot_df[!is.na(plot_df$prs.qvalues),])
nrow(plot_df[!is.na(plot_df$h2.qvalues),])

sum(plot_df$only_h2)
only_h2 = subset(plot_df, plot_df$only_h2)
mean(only_h2$h2.plot.diff)
sd(only_h2$h2.plot.diff)

plot_df = subset(plot_df, plot_df$only_h2 == F)

#pdf("plots/sex_h2_only.pdf")
#par(cex=1.5)
#plot(only_h2$h2.qvalues, only_h2$h2.plot.diff, log='x', ylab="|h2 top - h2 bottom|", xlab="q values")
#dev.off()
#hist(only_h2$h2.plot.diff)

###########
# GxE estimates
# sum across E variables for a trait?
# fill in the correct p value cutoff
# derive the 5 bin value 
###########
# temp = subset(data, data$prs.p.val < 1e-2)
# temp %>% group_by(trait) %>% summarise(int_VE = sum(prs.plot.diff), base_VE = mean(base.all)) %>%
#   summarize(mean(int_VE * base_VE), stderr = sd(int_VE * base_VE)/sqrt(n()), n=n()) # mean grouped by trait
# temp2 = subset(data, data$rg.plot.diff.p < 3.5e-2)
# temp2 %>% group_by(trait) %>% summarise(x = sum((1-rg.bin2) / 2 *100, na.rm=T)) %>%
#   summarise(mean(x, na.rm=T), stderr = sd(x)/sqrt(n()), n=n())
# 
# (1.62*18 + 3.85*22) / (18+22)
# mean(1-subset(plot_df, plot_df$is_scenario_1==T)$rg.plot.diff)
# mean(subset(plot_df, plot_df$is_scenario_2==T)$prs.plot.diff)
# mean(subset(plot_df, plot_df$is_scenario_3==T)$prs.plot.diff)


data$rgVE = (1-data$rg.bin2)/2 
# scenario 1 variance explained averaged across traits
s1_data = data %>% filter(rg.plot.diff.p < 2.5e-2) %>% group_by(trait) %>%  summarize(trait_sum = sum(rgVE, na.rm=T), trait_n= n()) #%>% 
mean(c(s1_data$trait_sum, rep(0, 33-nrow(s1_data)))) # add the correct amount of 0s
sd(c(s1_data$trait_sum, rep(0, 33-nrow(s1_data)))) / sqrt(33) # divide by SE

# scenario 1 variance explained averaged across traits
# 2.8% and se is 0.1%
s2_data = data %>% filter(prs.p.val < 2.7e-2 & h2.plot.diff.p < 1e-2) %>% 
  group_by(trait) %>%
  summarize(trait_sum = sum(prs.plot.diff))
mean(c(s2_data$trait_sum, rep(0, 33-nrow(s2_data)))) # add the correct amount of 0s
sd(c(s2_data$trait_sum, rep(0, 33-nrow(s2_data)))) / sqrt(33) # divide by SE

# scenario 3:
# 0.008% se is 0.003%
s3_data = data %>% filter(prs.p.val < 1e-2 & h2.plot.diff.p > 1e-2) %>% 
  group_by(trait) %>%
  summarize(trait_sum = sum(prs.plot.diff))
mean(c(s3_data$trait_sum, rep(0, 33-nrow(s3_data)))) # add the correct amount of 0s
sd(c(s3_data$trait_sum, rep(0, 33-nrow(s3_data)))) / sqrt(33) # divide by SE

############
# Plot new version of figure 3 with 
############
plot_df2 = data.frame(plot_df)
plot_df2 = merge(clean_trait, plot_df2, by="trait_clean", all.x = T)
plot_df2$scenario2_prs = plot_df2$prs.plot.diff * plot_df2$is_scenario_2
plot_df2$scenario3_prs = plot_df2$prs.plot.diff * plot_df2$is_scenario_3

plot_df2 = plot_df2 %>% mutate(
  rg.plot.diff = replace_na(rg.plot.diff, 0.99999),
  scenario2_prs = ifelse(scenario2_prs == 0, 0.001, scenario2_prs),
  scenario3_prs = ifelse(scenario3_prs == 0, 0.001, scenario3_prs),
  )

library(gridExtra)

p1 = ggplot(plot_df2, aes(x=(1-rg.plot.diff)/2, y=trait_order.x, fill=E_clean, labels=trait_clean)) + 
  geom_bar(stat="identity", fill='black') + 
  ylab("Phenotype") +
  xlab("Variance explained by GxSex") +
  ggtitle("  Scenario 1") +
  theme_bw() +
  theme(text = element_text(size=14)) + guides(fill="none") +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_y_discrete(limits=rev)


p2 = ggplot(plot_df2, aes(x=scenario2_prs/100, y=trait_order.x, fill=E_clean)) + 
  geom_bar(stat="identity", fill='black') + 
  ylab("Phenotype") +
  xlab("Variance explained by GxSex") +
  ggtitle("  Scenario 2") +
  theme_bw() + 
  theme(text = element_text(size=14)) + guides(fill="none") +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_y_discrete(limits=rev)

p3 = ggplot(plot_df2, aes(x=scenario3_prs/100, y=trait_order.x, fill=E_clean)) + 
  geom_bar(stat="identity", fill='black') + 
  ylab("Phenotype") +
  xlab("Variance explained by GxSex") +
  ggtitle("  Scenario 3") +
  theme_bw() + 
  theme(text = element_text(size=14)) + guides(fill="none") +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_y_discrete(limits=rev)# + 
  scale_x_continuous(labels = function(x) format(x, scientific = TRUE), guide = guide_axis(n.dodge=2))

p_blank = 
  ggplot(plot_df2, aes(x=prs.plot.diff, y=trait_order.x, fill=E_clean)) + 
  geom_bar(stat="identity", fill='black')  + 
  ylab("Disease/Trait") +
  xlab("") +
  ggtitle("removethis") + 
  theme_bw() + guides(fill="none") + 
  scale_y_discrete(limits=rev)

ggsave("plots/Figure5_v4/a.pdf", p1, width=40, height=65, units="mm", scale=2)
ggsave("plots/Figure5_v4/b.pdf", p2, width=40, height=65, units="mm", scale=2)
ggsave("plots/Figure5_v4/c.pdf", p3, width=40, height=65, units="mm", scale=2)
ggsave("plots/Figure5_v4/yaxis.pdf", p_blank, width=60, height=65, units='mm', scale=1.7)

write.table(plot_df, "./tables/Figure5_prs_only.csv", sep=",", quote=F, row.names=F)



############3
# new figure 6
#############
plot_rg = function(rg, se, col){
  bp=barplot(1-rg, ylim=c(0,0.1+se*3), ylab="1-Genetic correlation", col="white", border=col)
  arrows(bp, 1-rg+se, bp, 1-rg-se, angle=90, code=3, length=0.1, lwd=3)
  abline(h=1, lty=2, lwd=2)
}

plot_h2 = function(h2_vec, se_vec, col, xlab){
  bp=barplot(h2_vec, ylim=c(0,max(h2_vec+(se_vec*3))), ylab=expression(h^2), xlab=xlab, border=NA, col=col)
  arrows(bp, h2_vec+se_vec, bp, h2_vec-se_vec, angle=90, code=3, length=0.1)
}

plot_PS = function(PS_beta, PSxE_beta, P, col, xlab, nbins){
  E = 0:nbins
  y = PS_beta + PSxE_beta*E
  plot(E, y, ylab=expression(beta[PRS] + beta[PRSxE]*E), type='b', col=col, pch=19,
       xlab=xlab, ylim=c(min(y)-0.05, max(y)+0.05), xaxt='n')
  axis(1, at=E)
  legend("bottomright",paste0("P=", P), bty="n")
}
PS_main_effects = read.table("../newE/summary/PS_newE_summary.txt", header=F)

get_PS_main_effect = function(trait){
  mean(subset(PS_main_effects, grepl(pattern = trait, x = PS_main_effects))$V2)
}
prsxe = read.table("../figures_v2/table_output/prs_sex.txt", header=F)
names(prsxe) = c("Phenotype", "beta", "PRS.Z", "PRS.P", "PVE")
prsxe$E = rep("cov_SEX", nrow(prsxe))
prseqobj = qvalue(p = prsxe$PRS.P, fdr.level = 0.05, pi0=1)
prsxe$prs.sig = prseqobj$significant


pdf("plots/Figure6_230525.pdf", height=8)
layout(matrix(c(0,1,2,2,3,3,
                0,4,5,5,6,6,
                0,7,8,8,9,9), nrow = 3, ncol = 6, 
              byrow = TRUE))
#layout.show()
par(cex=1.0, mar=c(5,4,3,1))

example1_rg = subset(rg, rg$trait=="mental_NEUROTICISM" & rg$E == "cov_SEX")
example1_h2 = subset(h2, h2$trait =="mental_NEUROTICISM" & h2$E == "cov_SEX")
example1_PRS = subset(prsxe, prsxe$Phenotype == "mental_NEUROTICISM" & prsxe$E == "cov_SEX")
plot_rg(example1_rg$rg.bin2, example1_rg$rg.se.2, "black")
mtext("Neuroticism x Sex", line=1.4, at=-0.5, cex=1.5)
plot_PS(get_PS_main_effect("mental_NEUROTICISM"), example1_PRS$beta, "0.58", "black", "Sex \n(M /F)", 1)
plot_h2(c(example1_h2$h2.bin2.1, example1_h2$h2.bin2.2), 
        c(example1_h2$h2.se.1, example1_h2$h2.se.2),
        "black", "Sex \n(M / F)")

par(cex=1.0, mar=c(5,4,1.5,1))

example2_rg = subset(rg, rg$trait == "disease_AID_ALL" & rg$E == "cov_SEX")
example2_h2 = subset(h2, h2$trait=="disease_AID_ALL" & h2$E == "cov_SEX")
example2_PRS = subset(prsxe, prsxe$Phenotype == "disease_AID_ALL" & prsxe$E == "cov_SEX")
plot_rg(example2_rg$rg.bin2, example2_rg$rg.se.2, "black")
mtext("Autoimmune Disease x Sex", line=1.5, at=-1.4, cex=1.5)
plot_PS(get_PS_main_effect("disease_AID_ALL"), example2_PRS$beta, "2e-15", "black", "Sex \n(M / F)", 1)
plot_h2(c(example2_h2$h2.bin2.1,example2_h2$h2.bin2.2),
        c(example2_h2$h2.se.1,example2_h2$h2.se.2),
        "black", "Sex \n(M / F)")

example3_rg = subset(rg, rg$trait == "biochemistry_HDLcholesterol" & rg$E == "cov_SEX")
example3_h2 = subset(h2, h2$trait== "biochemistry_HDLcholesterol" & h2$E == "cov_SEX")
example3_PRS = subset(prsxe, prsxe$Phenotype == "biochemistry_HDLcholesterol"  & prsxe$E =="cov_SEX")
plot_rg(example3_rg$rg.bin2, example3_rg$rg.se.2, "black")
mtext("HDL Cholesterol x Sex", line=1.7, at=-2.7, cex=1.5)
plot_PS(get_PS_main_effect("biochemistry_HDLcholesterol"), example3_PRS$beta, "2.45e-17", "black", "Sex \n(M / F)", 1)
plot_h2(c(example3_h2$h2.bin2.1,example3_h2$h2.bin2.2),
        c(example3_h2$h2.se.1,example3_h2$h2.se.2),
        "black", "Sex \n(M / F)")
dev.off()



# ########
# # plot!
# ########
# library(gridExtra)
# p1 = ggplot(plot_df, aes(x=1-rg.plot.diff, y=trait)) + 
#   geom_bar(stat="identity", fill='black') + 
#   ylab("Phenotype") +
#   xlab("Scaled genetic correlation x E") +
#   theme_bw() + 
#   theme(text = element_text(size=14)) + guides(fill="none") +
#   theme(axis.text.y=element_blank(),
#         axis.title.y=element_blank(),
#         axis.ticks.y=element_blank())
# 
# p2 = ggplot(plot_df, aes(x=prs.plot.diff, y=trait)) + 
#   geom_bar(stat="identity", fill="black") + 
#   ylab("Phenotype") +
#   xlab("Scaled PRSxE") +
#   theme_bw() + 
#   theme(text = element_text(size=14)) + guides(fill="none") +
#   theme(axis.text.y=element_blank(),
#         axis.title.y=element_blank(),
#         axis.ticks.y=element_blank())
# 
# p3 = ggplot(plot_df, aes(x=h2.plot.diff, y=trait)) + 
#   geom_bar(stat="identity", fill='black') + 
#   ylab("Phenotype") +
#   xlab("Scaled SNP-heritability x E") +
#   theme_bw() +  
#   theme(text = element_text(size=14)) + guides(fill="none") +
#   theme(axis.text.y=element_blank(),
#         axis.title.y=element_blank(),
#         axis.ticks.y=element_blank())
# 
# p_blank = 
#   ggplot(plot_df, aes(x=prs.plot.diff, y=trait_clean)) + 
#   geom_bar(stat="identity", fill='black')  + 
#   ylab("Phenotype") +
#   xlab("") +
#   theme_bw() + guides(fill="none") 
# 
# 
# 
# ggsave("plots/Figure5/a.pdf", p1, width=40, height=60, units="mm", scale=2)
# ggsave("plots/Figure5/b.pdf", p2, width=40, height=60, units="mm", scale=2)
# ggsave("plots/Figure5/c.pdf", p3, width=40, height=60, units="mm", scale=2)
# ggsave("plots/Figure5/yaxis.pdf", p_blank, width=60, height=60, units='mm', scale=1.6)
# 
# write.table(plot_df, "./tables/Figure5.csv", sep=",", quote=F, row.names=F)
# 
# # pdf("plots/Figure5.pdf", 20, 10)
# # grid.arrange(p1, p2, p3, nrow = 1, widths=c(0.45, 0.275, 0.275))
# # dev.off()


