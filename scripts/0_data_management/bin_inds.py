import numpy as np
import pandas as pd
from sklearn import preprocessing
import sys

dat_337 =  pd.read_csv('data/337K-pheno-cov-newE-dietpc.tab', sep="\t")
dat_47 = pd.read_csv('data/49K-pheno-cov-newE-dietpc.tab', sep="\t")
dat_list = [dat_337, dat_47]
names = ['337K', '49K']
vars = ["cov_PHYS2_DAYS_VERY_ACTIVE",
        "cov_ALCOHOL_DESC_ORDER",
        "cov_TOWNSEND_DEPRIVATION",
        "other_TIME_TV",
        "other_NAP_TIME",
        "other_NO_WHEAT",
        "cov_DIET",
        "cov_Sleeplessness",
        "cov_ParticulateMatterAirPollution_irnt",
        "cov_SMOKING_STATUS",
]
binlist = [2, 5]
for idx,x in enumerate(dat_list):
    for NBINS in binlist:
        for VAR in vars:
            if VAR == "cov_SMOKING_STATUS":
                if NBINS == 5:
                    continue
                nonsmoking = x[x['cov_SMOKING_STATUS']==0]['IID']
                nonsmoking.to_csv('bins_newE/'+names[idx]+"_"+str(NBINS)+'bin_'+VAR+"_0", index=False, header=False)
                smoking = x[x['cov_SMOKING_STATUS']>0]['IID']
                smoking.to_csv('bins_newE/'+names[idx]+"_"+str(NBINS)+'bin_'+VAR+"_1", index=False, header=False)
            elif VAR == "other_NO_WHEAT":
                if NBINS == 5:
                    continue
                for i in range(NBINS):
                    curr = x[x['other_NO_WHEAT']==i]['IID']
                    curr.to_csv('bins_newE/'+names[idx]+"_"+str(NBINS)+'bin_'+VAR+"_"+str(i), index=False, header=False)
            elif VAR == "cov_Sleeplessness":
                if NBINS == 5:
                    continue
                never = x[x['cov_Sleeplessness']==1]['IID']
                never.to_csv('bins_newE/'+names[idx]+"_"+str(NBINS)+'bin_'+VAR+"_0", index=False, header=False)
                ever = x[x['cov_Sleeplessness']>1]['IID']
                ever.to_csv('bins_newE/'+names[idx]+"_"+str(NBINS)+'bin_'+VAR+"_1", index=False, header=False)
            else:
                x['bin_'+VAR] = pd.cut(x[VAR].rank(method = 'first'), NBINS, labels=list(range(NBINS)))
                for i in range(NBINS):
                    curr = x[x['bin_'+VAR]==i]['IID']
                    curr.to_csv('bins_newE/'+names[idx]+"_"+str(NBINS)+'bin_'+VAR+"_"+str(i), index=False, header=False)
