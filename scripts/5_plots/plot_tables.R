setwd("~/Documents/pricelab/ps_gxe/2023_05_10_tables/")
library(tidyverse)
library(viridis)
library(qvalue)
library(data.table)
inddata = data.frame(fread("../newE/337K-pheno-cov-newE-dietpc.tab"))
get_corr = function(trait, E){
  return(cor.test(inddata[,trait], inddata[,E], na.rm=T))
}
h2 = read.table("parsed/h2_values.csv", sep=",")
names(h2) = c('trait', "E", 'h2.bin5.1', 'h2.bin5.2', 'h2.bin5.3', 'h2.bin5.4', 'h2.bin5.5',
              'h2.se5.1', 'h2.se5.2', 'h2.se5.3', 'h2.se5.4', 'h2.se5.5', 'h2.bin2.1', 'h2.bin2.2',
              'h2.se.1', 'h2.se.2')
h2Z = read.table("parsed/h2_Z-scores.csv", sep=",")
names(h2Z) = c('trait', "E", "z.h2.5", 'z.h2.2')
rg = read.table("parsed/rg-values.csv", sep=",")
names(rg) = c('trait', 'E', 'rg.bin5', 'rg.se.5', 'rg.bin2', 'rg.se.2')
rgZ = read.table("parsed/rg_Z-scores.csv", sep=",")
names(rgZ) = c('trait', "E", "z.rg.5", 'z.rg.2')
prs = read.table("parsed/PRS.csv", sep=',')
names(prs) = c("trait", "E", "VE.diff", "base.PS", "base.E", 'base.all', "int.PS", "int.E", "int.PRSxE", "prs.p.val", 'int.all')
prs$E.diff = prs$int.E - prs$base.E
prs$VE.diff2 = (prs$int.PS + prs$int.PRSxE) - prs$base.PS
prs$prs.plot.diff = (prs$int.all - prs$base.all) / prs$base.PS # scale by PRS accuracy

h2_temp = merge(h2, h2Z, by=c('trait', 'E'))
rg_temp = merge(rg, rgZ, by=c('trait', 'E'))
temp_df = merge(rg_temp, h2_temp, by=c('trait', 'E'))#, by.y=c('trait', 'E'))
data = merge(temp_df, prs, by=c('trait', 'E'))

# scale h2 diff by overall h2 (mean across bins)
data$h2.5diff = abs(data$h2.bin5.5 - data$h2.bin5.1) #/ mean(c(data$h2.bin5.1, data$h2.bin5.2, data$h2.bin5.3, data$h2.bin5.4, data$h2.bin5.5), na.rm=T)
data$h2.2diff = abs(data$h2.bin2.2 - data$h2.bin2.1) #/ mean(c(data$h2.bin2.2, data$h2.bin2.1), na.rm=T)
data$h2.plot.diff = ifelse(is.nan(data$h2.5diff), data$h2.2diff, data$h2.5diff)
data$h2.plot.diff.Z = ifelse(is.nan(data$z.h2.5), data$z.h2.2, data$z.h2.5)
data$h2.plot.diff.p = pnorm(abs(data$h2.plot.diff.Z), lower.tail=F)*2

data$rg.plot.diff = ifelse(is.nan(data$rg.bin5), data$rg.bin2, data$rg.bin5)
data$rg.plot.diff.Z = ifelse(is.nan(data$z.rg.5), data$z.rg.2, data$z.rg.5)
data$rg.plot.diff.p = pnorm(q=data$rg.plot.diff.Z, lower.tail=FALSE)

data$mean.h2 = (data$h2.bin2.1 + data$h2.bin2.2) /2 * 100
data$prs.plot.diff = data$prs.plot.diff * data$mean.h2
####
h2_temp = data %>% group_by(trait) %>% summarize(meanh2 = mean(mean.h2))
h2_temp$is_disease = ifelse(grepl('disease', h2_temp$trait), T, F)
mean(subset(h2_temp, h2_temp$is_disease)$meanh2)
mean(subset(h2_temp, h2_temp$is_disease==F)$meanh2, na.rm=T)

prs_temp = data %>% group_by(trait) %>% summarize(meanprs = mean(base.PS), 
                                                  mean_int= mean(int.all - base.all),
                                                  mean_intprsxe = mean(int.PRSxE))
prs_temp$is_disease = ifelse(grepl('disease', prs_temp$trait), T, F)
mean(subset(prs_temp, prs_temp$is_disease)$meanprs)
mean(subset(prs_temp, prs_temp$is_disease==F)$meanprs, na.rm=T)
mean(subset(prs_temp, prs_temp$is_disease)$mean_int)
mean(subset(prs_temp, prs_temp$is_disease==F)$mean_int, na.rm=T)
mean(subset(prs_temp, prs_temp$is_disease)$mean_intprsxe)
mean(subset(prs_temp, prs_temp$is_disease==F)$mean_intprsxe, na.rm=T)


######
# subset df to just p-values and plot values
######
df = data.frame(data$trait, data$E, 
                data$rg.plot.diff.p, data$h2.plot.diff.p, data$prs.p.val,
                data$rg.plot.diff, data$h2.plot.diff, data$prs.plot.diff)
names(df) = c("trait", "E", 'rg.p', 'h2.p', 'prs.p', 'rg.val', 'h2.val', 'prs.val')
write.table(df, "parsed/combined.csv", sep=",", quote=F, row.names = F)
####
# Q value analysis
# we decided to do a per-E Q value (33 tests per E)
# lambda = 0 and fdr = 0.05 reduces to Benjamini and Hochberg 1995
####
get_q_vals = function(df, Evar, type){
  x = subset(df, df$E == Evar)
  if (type == "h2"){
    x.qobj = qvalue(x$h2.plot.diff.p, lambda=0, fdr.level = 0.05)
  }
  out_df = data.frame(x$trait, x$E, x.qobj$qvalues)
  names(out_df) = c("trait", "E", "qvalues")
  return(out_df)
}

h2.qobj <- qvalue(data$h2.plot.diff.p, fdr.level = 0.05)
rg.qobj <- qvalue(data$rg.plot.diff.p, fdr.level = 0.05)
prs.qobj <- qvalue(data$prs.p.val, fdr.level = 0.05)
#max(h2.qobj$qvalues[h2.qobj$pvalues <= 0.01], na.rm=T) # 0.01 pvalue cutoff is about 0.05 FDR
#max(rg.qobj$qvalues[rg.qobj$pvalues <= 0.015], na.rm=T) # 0.015 pvalue cutoff is about 0.04 FDR
#max(prs.qobj$qvalues[prs.qobj$pvalues <= 0.01], na.rm=T) # 0.01 pvalue cutoff is about 0.05 FDR

#hist(h2.qobj$pvalues)
#hist(rg.qobj$pvalues)
#hist(prs.qobj$pvalues)

#plot(h2.qobj)
#plot(rg.qobj)
#plot(prs.qobj)

#qqnorm(data$h2.plot.diff.Z)
#qqline(data$h2.plot.diff.Z)
#qqnorm(data$rg.plot.diff.Z)
#qqline(data$rg.plot.diff.Z)

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
clean_E = read.table("tables/plotting_names_matching_E.csv", sep=",")
names(clean_E) = c("E", "E_clean")
plot_df_trait = merge(plot_df_temp, clean_trait, by="trait")
plot_df = merge(plot_df_trait, clean_E, by="E")
plot_df$trait = gsub("_", " ", plot_df$trait)
plot_df$E = gsub("_", " ", plot_df$E)
plot_df$E = gsub("cov ", "", plot_df$E)
plot_df$E = gsub("other ", "", plot_df$E)

plot_df$is_scenario_1 = ifelse(!is.na(plot_df$rg.qvalues), T, F)
plot_df$is_scenario_2 = ifelse((!is.na(plot_df$h2.qvalues) & !is.na((plot_df$prs.qvalues))), T, F)
plot_df$is_scenario_3 = ifelse(!is.na(plot_df$prs.qvalues) & is.na(plot_df$h2.qvalues), T, F)
plot_df$only_h2 = ifelse((!is.na(plot_df$h2.qvalues) & is.na(plot_df$prs.qvalues) & is.na(plot_df$rg.qvalues)), T, F)

#nrow(plot_df[!is.na(plot_df$rg.qvalues),])
#nrow(plot_df[!is.na(plot_df$prs.qvalues),])
#nrow(plot_df[!is.na(plot_df$h2.qvalues),])

#sum(plot_df$is_scenario_1)
#sum(plot_df$is_scenario_2)
#sum(plot_df$is_scenario_3)
#sum(plot_df$is_scenario_1 & plot_df$is_scenario_2)
#sum(plot_df$is_scenario_1 & plot_df$is_scenario_3)

sum(plot_df$only_h2)
only_h2 = subset(plot_df, plot_df$only_h2)
#mean(only_h2$h2.plot.diff)
#sd(only_h2$h2.plot.diff)
#pdf("plots/h2_only.pdf")
#par(cex=1.5)
#plot(only_h2$h2.qvalues, only_h2$h2.plot.diff, log='x', ylab="|h2 top - h2 bottom|", xlab="q values")
#dev.off()
#hist(only_h2$h2.plot.diff)
out_only_h2 = data.frame(only_h2$trait_clean, only_h2$E_clean, only_h2$h2.plot.diff, only_h2$h2.qvalues)
names(out_only_h2) = c("Trait", "E", "h2 difference", "q value FDR 0.05")
write.table(out_only_h2, "tables/only_h2.csv", quote=F, sep=",", row.names=F)
plot_df = subset(plot_df, plot_df$only_h2 == F)

######
# output the h2 only trait-E pairs
######
h2only_save = data.frame(cbind(only_h2$trait_clean, only_h2$E_clean, only_h2$h2.plot.diff, only_h2$h2.qvalues))
h2only_save$X3 = round(as.numeric(h2only_save$X3),3)
names(h2only_save) = c("Trait", "E", "h2 difference", "qvalue (5% FDR)")
write.csv(h2only_save, "tables/h2only.csv", row.names = F, quote=F)

###########
# GxE estimates
# sum across E variables for a trait?
# fill in the correct p value cutoff
# derive the 5 bin value 
###########
data$rgVE = (1-data$rg.bin5)/10
# scenario 1 variance explained averaged across traits
s1_data = data %>% filter(rg.plot.diff.p < 0.78e-2) %>% group_by(trait) %>%  summarize(trait_sum = sum(rgVE, na.rm=T), trait_n= n()) #%>% 
mean(c(s1_data$trait_sum, rep(0, 33-nrow(s1_data)))) # add the correct amount of 0s
sd(c(s1_data$trait_sum, rep(0, 33-nrow(s1_data)))) / sqrt(33) # divide by SE

  #summarize(trait_mean = sum(trait_sum)/33, trait_se = sd(trait_sum)/sqrt(n()), n = n()) # this is not in percent units (so we multiply)
# scenario 2:
s2_data = data %>% filter(prs.p.val < 1e-2 & h2.plot.diff.p < 2e-2) %>% 
  group_by(trait) %>%
  summarize(trait_sum = sum(prs.plot.diff))
mean(c(s2_data$trait_sum, rep(0, 33-nrow(s2_data)))) # add the correct amount of 0s
sd(c(s2_data$trait_sum, rep(0, 33-nrow(s2_data)))) / sqrt(33) # divide by SE

# scenario 3:
s3_data = data %>% filter(prs.p.val < 1e-2 & h2.plot.diff.p > 1e-2) %>% 
  group_by(trait) %>%
  summarize(trait_sum = sum(prs.plot.diff))
mean(c(s3_data$trait_sum, rep(0, 33-nrow(s3_data)))) # add the correct amount of 0s
sd(c(s3_data$trait_sum, rep(0, 33-nrow(s3_data)))) / sqrt(33) # divide by SE

0.2 + 0.07 + 0.003


# mean h2
mean(mean(data$h2.bin5.5, na.rm=T),
mean(data$h2.bin5.4, na.rm=T),
mean(data$h2.bin5.3, na.rm=T),
mean(data$h2.bin5.2, na.rm=T),
mean(data$h2.bin5.1, na.rm=T))
sd(data$h2.bin5.3, na.rm=T) * 1/(33)
mean(mean(data$h2.bin2.2, na.rm=T),
mean(data$h2.bin2.1, na.rm=T))

nrow(subset(plot_df, plot_df$is_scenario_2 | plot_df$is_scenario_3))
mean(subset(plot_df, plot_df$is_scenario_2)$h2.plot.diff)

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


plot_df2$rg.plot.diff2 = (1-plot_df2$rg.plot.diff)/10 
#plot_df2$scenario2_prs[plot_df2$scenario2_prs == 0] = NA
#plot_df2$scenario3_prs[plot_df2$scenario3_prs == 0] = NA
plot_df2 = subset(plot_df2, plot_df2$E!="NO WHEAT")
library(gridExtra)
library(pals)
pal = kelly(12)[2:12]

group.colors <- c("Alcohol cons." = pal[1],
                  "Air pollution" = pal[2],
                  "Phys. Act." =  pal[3],
                  "Sleeplessness" = pal[4],
                  "Smoking"= pal[5],
                  "Townsend dep." = pal[6],
                  "Napping" = pal[7],
                  "Television" = pal[8], 
                  "Diet" = pal[9],
                  "Wheat cons." = pal[10])

plot_scenario1 = subset(plot_df2, !is.na(plot_df2$rg.plot.diff2))[,c("trait_clean", "E_clean", "rg.plot.diff2")]
plot_scenario1 = merge(plot_scenario1, clean_trait, by="trait_clean", all.y = T)
p1 = ggplot(plot_scenario1, aes(x=rg.plot.diff2, y=trait_order, fill=E_clean)) + 
  geom_bar(stat="identity", position = position_dodge()) + 
  ylab("Phenotype") +
  xlab("Variance explained by GxE") +
  ggtitle("  Scenario 1") +
  theme_bw() + scale_fill_manual(values=group.colors) +
  theme(text = element_text(size=14)) + guides(fill="none") +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank()) + 
  scale_y_discrete(limits=rev)



plot_scenario2 = subset(plot_df2, !is.na(plot_df2$scenario2_prs))[,c("trait_clean", "E_clean", "scenario2_prs")]
plot_scenario2 = merge(plot_scenario2, clean_trait, by="trait_clean", all.y = T)
p2 = ggplot(plot_scenario2, aes(x=scenario2_prs/100, y=trait_order, fill=E_clean)) + 
  geom_bar(stat="identity", position = position_dodge()) + 
  ylab("Phenotype") +
  xlab("Variance explained by GxE") +
  ggtitle("  Scenario 2") +
  theme_bw() + scale_fill_manual(values=group.colors) +
  theme(text = element_text(size=14)) + guides(fill="none") +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_y_discrete(limits=rev)

plot_scenario3 = subset(plot_df2, !is.na(plot_df2$scenario3_prs))[,c("trait_clean", "E_clean", "scenario3_prs")]
plot_scenario3 = merge(plot_scenario3, clean_trait, by="trait_clean", all.y = T)
p3 = ggplot(plot_scenario3, aes(x=scenario3_prs/100, y=trait_order, fill=E_clean)) + 
  geom_bar(stat="identity", position = position_dodge()) + 
  ylab("Phenotype") +
  xlab("Variance explained by GxE") +
  ggtitle("  Scenario 3") +
  theme_bw() + scale_fill_manual(values=group.colors) +
  theme(text = element_text(size=14)) + guides(fill="none") +
  theme(axis.text.y=element_blank(),
        axis.title.y=element_blank(),
        axis.ticks.y=element_blank()) +
  scale_y_discrete(limits=rev)

p_blank = 
  ggplot(plot_scenario3, aes(x=scenario3_prs, y=trait_order, fill=E_clean)) + 
  geom_bar(stat="identity", position=position_dodge())  + 
  ylab("Disease/Trait") +
  xlab("") +
  ggtitle("removethis") + 
  theme_bw() + scale_fill_manual(values=group.colors, drop=F, name="E var") +
  scale_y_discrete(limits=rev)

ggsave("plots/Figure3_v5/a.pdf", p1, width=40, height=60, units="mm", scale=2)
ggsave("plots/Figure3_v5/b.pdf", p2, width=40, height=60, units="mm", scale=2)
ggsave("plots/Figure3_v5/c.pdf", p3, width=40, height=60, units="mm", scale=2)
ggsave("plots/Figure3_v5/yaxis.pdf", p_blank, width=60, height=60, units='mm', scale=1.9)

write.table(plot_df, "./tables/Figure3_prs_only.csv", sep=",", quote=F, row.names=F)


#####
############
# Old, stacked version of figure 3
############
# plot_df2 = data.frame(plot_df)
# plot_df2$scenario2_prs = plot_df2$prs.plot.diff * plot_df2$is_scenario_2
# plot_df2$scenario3_prs = plot_df2$prs.plot.diff * plot_df2$is_scenario_3
# plot_df2$rg.plot.diff2 = (1-plot_df2$rg.plot.diff)/10 
# plot_df2$scenario2_prs[plot_df2$scenario2_prs == 0] = NA
# plot_df2$scenario3_prs[plot_df2$scenario3_prs == 0] = NA
# plot_df2 = subset(plot_df2, plot_df2$E!="NO WHEAT")
# library(gridExtra)
# library(pals)
# pal = kelly(12)[2:12]
# 
# group.colors <- c("Alcohol cons." = pal[1],
#                   "Air pollution" = pal[2],
#                   "Phys. Act." =  pal[3],
#                   "Sleeplessness" = pal[4],
#                   "Smoking"= pal[5],
#                   "Townsend dep." = pal[6],
#                   "Napping" = pal[7],
#                   "Television" = pal[8], 
#                   "Diet" = pal[9],
#                   "Wheat cons." = pal[10])
# 
# p1 = ggplot(plot_df2, aes(x=rg.plot.diff2, y=trait_clean, fill=E_clean)) + 
#   geom_bar(stat="identity") + 
#   ylab("Phenotype") +
#   xlab("Variance explained by GxE") +
#   ggtitle("  Scenario 1") +
#   theme_bw() + scale_fill_manual(values=group.colors) +
#   theme(text = element_text(size=14)) + guides(fill="none") +
#   theme(axis.text.y=element_blank(),
#         axis.title.y=element_blank(),
#         axis.ticks.y=element_blank()) +
#   scale_y_discrete(limits=rev)
# 
# 
# #p2 = 
# ggplot(plot_df2, aes(x=scenario2_prs/100, y=trait_clean, fill=E_clean)) + 
#   geom_bar(stat="identity", position = position_dodge()) + 
#   ylab("Phenotype") +
#   xlab("Variance explained by GxE") +
#   ggtitle("  Scenario 2") +
#   theme_bw() + scale_fill_manual(values=group.colors) +
#   theme(text = element_text(size=14)) + guides(fill="none") +
#   theme(axis.text.y=element_blank(),
#         axis.title.y=element_blank(),
#         axis.ticks.y=element_blank()) +
#   scale_y_discrete(limits=rev)
# 
# p3 = ggplot(plot_df2, aes(x=scenario3_prs/100, y=trait_clean, fill=E_clean)) + 
#   geom_bar(stat="identity") + 
#   ylab("Phenotype") +
#   xlab("Variance explained by GxE") +
#   ggtitle("  Scenario 3") +
#   theme_bw() + scale_fill_manual(values=group.colors) +
#   theme(text = element_text(size=14)) + guides(fill="none") +
#   theme(axis.text.y=element_blank(),
#         axis.title.y=element_blank(),
#         axis.ticks.y=element_blank()) +
#   scale_y_discrete(limits=rev)
# 
# p_blank = 
#   ggplot(plot_df2, aes(x=prs.plot.diff.scaled, y=trait_clean, fill=E_clean)) + 
#   geom_bar(stat="identity")  + 
#   ylab("Disease/Trait") +
#   xlab("") +
#   ggtitle("removethis") + 
#   theme_bw() + scale_fill_manual(values=group.colors, drop=F, name="E var") +
#   scale_y_discrete(limits=rev)
# 
# ggsave("plots/Figure3_v4/a.pdf", p1, width=40, height=60, units="mm", scale=2)
# ggsave("plots/Figure3_v4/b.pdf", p2, width=40, height=60, units="mm", scale=2)
# ggsave("plots/Figure3_v4/c.pdf", p3, width=40, height=60, units="mm", scale=2)
# ggsave("plots/Figure3_v4/yaxis.pdf", p_blank, width=60, height=60, units='mm', scale=1.5)
# 
# write.table(plot_df, "./tables/Figure3_prs_only.csv", sep=",", quote=F, row.names=F)
# 
# 
# #####
# t2d alcohol consumption prop test
####
qqa = nrow(subset(inddata, inddata$disease_T2D == 1 & inddata$cov_ALCOHOL_DESC_ORDER > 4))
qqb = nrow(subset(inddata, inddata$disease_T2D == 0 & inddata$cov_ALCOHOL_DESC_ORDER > 4))
qqc = nrow(subset(inddata, inddata$disease_T2D == 1 & inddata$cov_ALCOHOL_DESC_ORDER < 2))
qqd = nrow(subset(inddata, inddata$disease_T2D == 0 & inddata$cov_ALCOHOL_DESC_ORDER < 2))
prop.test(x = c(qqa, qqc), n = c(qqb, qqd))

############3
# new figure 4
#############
plot_rg = function(rg, se, col){
  bp=barplot(1-rg, ylim=c(0,0.1+se), ylab="1-Genetic correlation", col="white", border=col)
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
prsxe = read.table("../figures_v2/table_output/PRSxE.txt", header=F)
names(prsxe) = c("Phenotype", "E", "beta", "PRS.Z", "PRS.P", "PVE")
prseqobj = qvalue(p = prsxe$PRS.P, fdr.level = 0.05, pi0=1)
prsxe$prs.sig = prseqobj$significant

pdf("plots/Figure4_230517.pdf", height=8, width=7.08661)
layout(matrix(c(0,1,2,2,3,3,
                0,4,5,5,6,6,
                0,7,8,8,9,9), nrow = 3, ncol = 6, 
              byrow = TRUE))
#layout.show()
par(cex=1.0, mar=c(5,4,3,1))

example1_rg = subset(rg, rg$trait=="blood_WHITE_COUNT" & rg$E == "cov_SMOKING_STATUS")
example1_h2 = subset(h2, h2$trait =="blood_WHITE_COUNT" & h2$E == "cov_SMOKING_STATUS")
example1_PRS = subset(prsxe, prsxe$Phenotype == "blood_WHITE_COUNT" & prsxe$E == "cov_SMOKING_STATUS")
plot_rg(example1_rg$rg.bin2, example1_rg$rg.se.2, "red")
mtext("  White Blood Cell Count x Smoking", line=1.4, at=-0.5, cex=1.5)
plot_PS(get_PS_main_effect("blood_WHITE_COUNT"), example1_PRS$beta, "0.45", "red", "Smoking \n(No/Yes)", 1)
plot_h2(c(example1_h2$h2.bin2.1, example1_h2$h2.bin2.2), 
        c(example1_h2$h2.se.1, example1_h2$h2.se.2),
        "red", "Smoking \n(No/Yes)")

par(cex=1.0, mar=c(5,4,1.5,1))

example2_rg = subset(rg, rg$trait == "disease_T2D" & rg$E == "cov_ALCOHOL_DESC_ORDER")
example2_h2 = subset(h2, h2$trait=="disease_T2D" & h2$E == "cov_ALCOHOL_DESC_ORDER")
example2_PRS = subset(prsxe, prsxe$Phenotype == "disease_T2D" & prsxe$E == "cov_ALCOHOL_DESC_ORDER")
plot_rg(example2_rg$rg.bin5, example2_rg$rg.se.5, "darkgreen")
mtext("T2D x Alcohol Consumption", line=1.5, at=-1.4, cex=1.5)
plot_PS(get_PS_main_effect("disease_T2D"), example2_PRS$beta, "1.3e-13", "darkgreen", "Bin of alcohol cons.\n(low to high)", 4)
plot_h2(c(example2_h2$h2.bin5.1,example2_h2$h2.bin5.2,example2_h2$h2.bin5.3,example2_h2$h2.bin5.4,example2_h2$h2.bin5.5),
        c(example2_h2$h2.se5.1,example2_h2$h2.se5.2,example2_h2$h2.se5.3,example2_h2$h2.se5.4,example2_h2$h2.se5.5),
        "darkgreen", "Bin of alcohol cons.\n(low to high)")

example3_rg = subset(rg, rg$trait == "body_WHRadjBMIz" & rg$E == "other_TIME_TV")
example3_h2 = subset(h2, h2$trait== "body_WHRadjBMIz" & h2$E == "other_TIME_TV")
example3_PRS = subset(prsxe, prsxe$Phenotype == "body_WHRadjBMIz"  & prsxe$E =="other_TIME_TV")
plot_rg(example3_rg$rg.bin5, example3_rg$rg.se.5, "darkblue")
mtext("WHRadjBMI x TV Time", line=1.7, at=-2.7, cex=1.5)
plot_PS(get_PS_main_effect("body_WHRadjBMIz"), example3_PRS$beta, "5.9e-3", "darkblue", "Bin of TV time\n(low to high)", 4)
plot_h2(c(example3_h2$h2.bin5.1,example3_h2$h2.bin5.2,example3_h2$h2.bin5.3,example3_h2$h2.bin5.4,example3_h2$h2.bin5.5),
        c(example3_h2$h2.se5.1,example3_h2$h2.se5.2,example3_h2$h2.se5.3,example3_h2$h2.se5.4,example3_h2$h2.se5.5),
        "darkblue", "Bin of TV time\n(low to high)")
dev.off()


##############
# example plots
#############
pdf("plots/h2-blood_WHITE_COUNT-cov_Smoking_Status.pdf", 5, 5)
par(cex=1.2, lwd=2)
s1e = subset(data, data$trait == 'blood_WHITE_COUNT' & data$E == 'cov_SMOKING_STATUS')
h21=c(s1e$h2.bin2.1,s1e$h2.bin2.2)
se1=c(s1e$h2.se.1, s1e$h2.se.2)*1.96
bp = barplot(h21,  ylab='h2', xlab="Smoking status", ylim=c(0, 0.35), names.arg = c(0,1))
arrows(bp, h21+se1, bp, h21-se1, code=2, length=0.00, angle=90)
abline(h=mean(h21), lty=2)
dev.off()

pdf("plots/h2-disease_T2D-cov_ALCOHOL_DESC_ORDER.pdf", 5, 5)
par(cex=1.2, lwd=2)
s2e = subset(data, data$trait == 'disease_T2D' & data$E == 'cov_ALCOHOL_DESC_ORDER')
h22=c(s2e$h2.bin5.1,s2e$h2.bin5.2,s2e$h2.bin5.3,s2e$h2.bin5.4,s2e$h2.bin5.5)
se2=c(s2e$h2.se5.1,s2e$h2.se5.2,s2e$h2.se5.3,s2e$h2.se5.4,s2e$h2.se5.5)*1.96
bp = barplot(h22,  ylab='h2', xlab="Alcohol consumption", ylim=c(0, 0.5), names.arg = 1:5)
arrows(bp, h22+se2, bp, h22-se2, code=2, length=0.00, angle=90)
abline(h=mean(h22), lty=2)
dev.off()


pdf("plots/h2-biochemistry_Triglycerides-cov_DIET.pdf", 5, 5)
par(cex=1.2, lwd=2)
s3e = subset(data, data$trait == 'biochemistry_Triglycerides' & data$E == 'cov_DIET')
h23=c(s3e$h2.bin5.1,s3e$h2.bin5.2,s3e$h2.bin5.3,s3e$h2.bin5.4,s3e$h2.bin5.5)
se3=c(s3e$h2.se5.1,s3e$h2.se5.2,s3e$h2.se5.3,s3e$h2.se5.4,s3e$h2.se5.5)*1.96
bp = barplot(h23,  ylab='h2', xlab="Diet", ylim=c(0, 0.36), names.arg = 1:5)
arrows(bp, h23+se3, bp, h23-se3, code=2, length=0.00, angle=90)
abline(h=mean(h23), lty=2)
dev.off()

pdf("plots/h2-body_WHRadjBMIz-other_TIME_TV.pdf", 5, 5)
par(cex=1.2, lwd=2)
s3e = subset(data, data$trait == 'body_WHRadjBMIz' & data$E == 'other_TIME_TV')
h23=c(s3e$h2.bin5.1,s3e$h2.bin5.2,s3e$h2.bin5.3,s3e$h2.bin5.4,s3e$h2.bin5.5)
se3=c(s3e$h2.se5.1,s3e$h2.se5.2,s3e$h2.se5.3,s3e$h2.se5.4,s3e$h2.se5.5)*1.96
bp = barplot(h23,  ylab='h2', xlab="Time spent watching TV", ylim=c(0, 0.36), names.arg = 1:5)
arrows(bp, h23+se3, bp, h23-se3, code=2, length=0.00, angle=90)
abline(h=mean(h23), lty=2)
dev.off()



max(rg.qobj$qvalues[rg.qobj$pvalues <= 0.011], na.rm=T)
max(h2.qobj$qvalues[h2.qobj$pvalues <= 0.03], na.rm=T)
h2_change = subset(data, data$rg.plot.diff.p < 0.011) # roughly 5% FDR, see line above1
h2_change$fill = ifelse(h2_change$h2.plot.diff.p < 0.03, 19, 1)
#h2_change = subset(plot_df, plot_df$is_scenario_1==T)
pdf("plots/scenario1_h2change.pdf", 5, 5)
par(cex=1.3)
plot(h2_change$rg.plot.diff, h2_change$h2.plot.diff, pch=h2_change$fill, ylab="|h2 top - h2 bottom|", xlab="rg")
dev.off()
cor.test(h2_change$rg.plot.diff, h2_change$h2.plot.diff)

#####
# pieces of figure 3
######
# 
# pdf("plots/rg.pdf", width=10 , height=5.625)
# ggplot(plot_df, aes(x=1-rg.plot.diff, y=trait, fill=E)) + 
#   geom_bar(stat="identity") + 
#   ylab("Phenotype") +
#   xlab("1-rg") +
#   theme_bw() + scale_fill_manual(values=group.colors) +
#   theme(text = element_text(size=16))
# dev.off()
# 
# pdf("plots/prs.pdf", width=10 , height=5.625 )
# ggplot(plot_df, aes(x=prs.plot.diff, y=trait, fill=E)) + 
#   geom_bar(stat="identity") + 
#   ylab("Phenotype") +
#   xlab("Var[Interaction] - Var[Base]") +
#   theme_bw() + scale_fill_manual(values=group.colors) +
#   theme(text = element_text(size=14))
# dev.off()
# 
# pdf("plots/h2.pdf", width=15, height=10)
# ggplot(plot_df, aes(x=h2.plot.diff, y=trait, fill=E)) + 
#   geom_bar(stat="identity") + 
#   ylab("Phenotype") +
#   xlab("|h2 top - h2 bottom|") +
#   theme_bw() +  scale_fill_manual(values=group.colors) +
#   theme(text = element_text(size=14))
# dev.off()
# 
# pdf('plots/prs_h2_scenario2.pdf', width=10 , height=5.625 )
# plot_scenario2 = plot_df %>% filter(is_scenario_2 == T)
# p4 = ggplot(plot_scenario2, aes(x=prs.plot.diff, y=trait, fill=E)) + 
#   geom_bar(stat="identity") + 
#   ylab("Phenotype") +
#   xlab("Var[Int.]-Var[Base]") +
#   theme_bw() + scale_fill_manual(values=group.colors) +
#   theme(text = element_text(size=14)) + guides(fill="none") 
# 
# p5 = ggplot(plot_scenario2, aes(x=h2.plot.diff, y=trait, fill=E)) + 
#   geom_bar(stat="identity") + 
#   ylab("Phenotype") +
#   xlab("|h2 top - h2 bottom|") +
#   theme_bw() +  scale_fill_manual(values=group.colors) +
#   theme(text = element_text(size=14)) +
#   theme(axis.text.y=element_blank(),
#         axis.title.y=element_blank(),
#         axis.ticks.y=element_blank())
# grid.arrange(p4, p5, nrow=1)
# dev.off()
# 
# pdf('plots/prs_scenario3.pdf', width=10 , height=5.625 )
# plot_scenario3 = plot_df %>% filter(is_scenario_3 == T)
# ggplot(plot_scenario3, aes(x=prs.plot.diff, y=trait, fill=E)) + 
#   geom_bar(stat="identity") + 
#   ylab("Phenotype") +
#   xlab("Var[Int.] - Var[Base]") +
#   theme_bw() + scale_fill_manual(values=group.colors) +
#   theme(text = element_text(size=14))
# dev.off()
# 


#########
# replace body WHRadjBMIz with Triglycerides x Diet
##########
# 
# pdf("plots/Figure4_230517_TGxDiet.pdf", height=8, width=7.08661)
# layout(matrix(c(0,1,2,2,3,3,
#                 0,4,5,5,6,6,
#                 0,7,8,8,9,9), nrow = 3, ncol = 6, 
#               byrow = TRUE))
# #layout.show()
# par(cex=1.0, mar=c(5,4,3,1))
# 
# example1_rg = subset(rg, rg$trait=="blood_WHITE_COUNT" & rg$E == "cov_SMOKING_STATUS")
# example1_h2 = subset(h2, h2$trait =="blood_WHITE_COUNT" & h2$E == "cov_SMOKING_STATUS")
# example1_PRS = subset(prsxe, prsxe$Phenotype == "blood_WHITE_COUNT" & prsxe$E == "cov_SMOKING_STATUS")
# plot_rg(example1_rg$rg.bin2, example1_rg$rg.se.2, "red")
# mtext("  White Blood Cell Count x Smoking", line=1.4, at=-0.5, cex=1.5)
# plot_PS(get_PS_main_effect("blood_WHITE_COUNT"), example1_PRS$beta, "0.45", "red", "Smoking \n(No/Yes)", 1)
# plot_h2(c(example1_h2$h2.bin2.1, example1_h2$h2.bin2.2), 
#         c(example1_h2$h2.se.1, example1_h2$h2.se.2),
#         "red", "Smoking \n(No/Yes)")
# 
# par(cex=1.0, mar=c(5,4,1.5,1))
# 
# example2_rg = subset(rg, rg$trait == "disease_T2D" & rg$E == "cov_ALCOHOL_DESC_ORDER")
# example2_h2 = subset(h2, h2$trait=="disease_T2D" & h2$E == "cov_ALCOHOL_DESC_ORDER")
# example2_PRS = subset(prsxe, prsxe$Phenotype == "disease_T2D" & prsxe$E == "cov_ALCOHOL_DESC_ORDER")
# plot_rg(example2_rg$rg.bin5, example2_rg$rg.se.5, "darkgreen")
# mtext("T2D x Alcohol Consumption", line=1.5, at=-1.4, cex=1.5)
# plot_PS(get_PS_main_effect("disease_T2D"), example2_PRS$beta, "1.3e-13", "darkgreen", "Bin of alcohol cons.\n(low to high)", 4)
# plot_h2(c(example2_h2$h2.bin5.1,example2_h2$h2.bin5.2,example2_h2$h2.bin5.3,example2_h2$h2.bin5.4,example2_h2$h2.bin5.5),
#         c(example2_h2$h2.se5.1,example2_h2$h2.se5.2,example2_h2$h2.se5.3,example2_h2$h2.se5.4,example2_h2$h2.se5.5),
#         "darkgreen", "Bin of alcohol cons.\n(low to high)")
# 
# example3_rg = subset(rg, rg$trait == "biochemistry_Triglycerides" & rg$E == "cov_DIET")
# example3_h2 = subset(h2, h2$trait== "biochemistry_Triglycerides" & h2$E == "cov_DIET")
# example3_PRS = subset(prsxe, prsxe$Phenotype == "biochemistry_Triglycerides"  & prsxe$E =="cov_DIET")
# plot_rg(example3_rg$rg.bin2, example3_rg$rg.se.2, "darkblue")
# mtext("Triglycerides x Diet", line=1.7, at=-2.7, cex=1.5)
# plot_PS(get_PS_main_effect("biochemistry_Triglycerides"), example3_PRS$beta, "4.2e-5", "darkblue", "Bin of Diet", 4)
# plot_h2(c(example3_h2$h2.bin5.1,example3_h2$h2.bin5.2,example3_h2$h2.bin5.3,example3_h2$h2.bin5.4,example3_h2$h2.bin5.5),
#         c(example3_h2$h2.se5.1,example3_h2$h2.se5.2,example3_h2$h2.se5.3,example3_h2$h2.se5.4,example3_h2$h2.se5.5),
#         "darkblue", "Bin of Diet")
# 
# dev.off()

# 
# ###
# # 
# t2d = subset(plot_df, plot_df$trait == "disease T2D")
# t2d_df = data.frame(t2d$trait, t2d$E, t2d$is_scenario_1, t2d$is_scenario_2, t2d$is_scenario_3)
# library(tidyverse)
# gather(t2d_df, 'trait','E')
# 
# #######
# #
# num_scenarios = list()
# traits = unique(plot_df$trait)
# for (trait in traits){
#   append(num_scenarios, colSums(subset(plot_df, plot_df$trait == trait)[, c("is_scenario_1", "is_scenario_2", "is_scenario_3")]))
# }
########
# plot!
########
# library(gridExtra)
# library(pals)
# pal = kelly(12)[2:12]
# 
# group.colors <- c("Alcohol cons." = pal[1],
#                   "Air pollution" = pal[2],
#                   "Phys. Act." =  pal[3],
#                   "Sleeplessness" = pal[4],
#                   "Smoking"= pal[5],
#                   "Townsend dep." = pal[6],
#                   "Napping" = pal[7],
#                   "Television" = pal[8], 
#                   "Diet" = pal[9],
#                   "Wheat cons." = pal[10])
# 
# p1 = ggplot(plot_df, aes(x=1-rg.plot.diff, y=trait_clean, fill=E_clean)) + 
#   geom_bar(stat="identity") + 
#   ylab("Phenotype") +
#   xlab("Scaled genetic correlation x E") +
#   theme_bw() + scale_fill_manual(values=group.colors) +
#   theme(text = element_text(size=14)) + guides(fill="none") +
#   theme(axis.text.y=element_blank(),
#         axis.title.y=element_blank(),
#         axis.ticks.y=element_blank())
# 
# 
# p2 = ggplot(plot_df, aes(x=prs.plot.diff, y=trait_clean, fill=E_clean)) + 
#   geom_bar(stat="identity") + 
#   ylab("Phenotype") +
#   xlab("Scaled PRSxE") +
#   theme_bw() + scale_fill_manual(values=group.colors) +
#   theme(text = element_text(size=14)) + guides(fill="none") +
#   theme(axis.text.y=element_blank(),
#         axis.title.y=element_blank(),
#         axis.ticks.y=element_blank())
# 
# p3 = ggplot(plot_df, aes(x=h2.plot.diff, y=trait_clean, fill=E_clean)) + 
#   geom_bar(stat="identity") + 
#   ylab("Phenotype") +
#   xlab("Scaled SNP-heritability x E") +
#   theme_bw() + scale_fill_manual(values=group.colors) +
#   theme(text = element_text(size=14)) + guides(fill="none") +
#   theme(axis.text.y=element_blank(),
#         axis.title.y=element_blank(),
#         axis.ticks.y=element_blank())
# 
# p_blank = 
#   ggplot(plot_df, aes(x=prs.plot.diff, y=trait_clean, fill=E_clean)) + 
#   geom_bar(stat="identity")  + 
#   ylab("Phenotype") +
#   xlab("Scaled genetic correlation x E") +
#   theme_bw() + scale_fill_manual(values=group.colors, drop=F, name="E var")
# 
# 
# # pdf("plots/Figure3test.pdf", 3, width=7.08661)
# # grid.arrange(p1, p2, p3, nrow = 1, widths=c(0.38, 0.25, 0.37))
# # dev.off()
# 
# ggsave("plots/Figure3/a.pdf", p1, width=40, height=60, units="mm", scale=2)
# ggsave("plots/Figure3/b.pdf", p2, width=40, height=60, units="mm", scale=2)
# ggsave("plots/Figure3/c.pdf", p3, width=40, height=60, units="mm", scale=2)
# ggsave("plots/Figure3/yaxis.pdf", p_blank, width=60, height=60, units='mm', scale=1.5)
# 
# write.table(plot_df, "./tables/Figure3.csv", sep=",", quote=F, row.names=F)


#temp3 = subset(data, data$prs.p.val < 1e-2 & data$h2.plot.diff.p < 1e-2)

## Clean this section up because the numbers are important!
# Want to revise this whole script and re-run it
# temp = subset(data, data$prs.p.val < 1e-2)
# temp %>% group_by(trait) %>% summarise(int_VE = sum(prs.plot.diff), base_VE = mean(base.PS)) %>%
#   summarize(mean(int_VE), sd(int_VE) * (1/sqrt(21)), n=n()) # mean grouped by trait
# 
# (mean(temp$prs.plot.diff) / mean(temp$base.PS)) * mean(temp$h2.bin5.3, na.rm=T)*100
# temp2 = subset(data, data$rg.plot.diff.p < 1e-2)
# mean((1-temp2$rg.bin5) / 10 - 0.01, na.rm=T) 
# temp2 %>% group_by(trait) %>% summarise(x = sum((1-rg.bin2) / 2 *100, na.rm=T)) %>%
#   summarise(mean(x, na.rm=T), sd(x, na.rm=T) * 1/sqrt(13), n=n())

### look across all traits, estimate PVE, truncate the negative estimates at 0
# take the mean across all positive estimates
# do we want to add the 0/negative estimates to the denominator?
# do we want to sum across all the E variables?

#plot(data$rg.bin5, (1-data$rg.bin5)/10) 
#mean(abs(temp3$h2.bin2.1 - temp3$h2.bin2.2)[1:9]) / 2