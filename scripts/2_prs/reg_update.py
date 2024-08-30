import argparse
from functools import reduce
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.stats.anova import anova_lm
from sklearn import preprocessing
import sys

def run_model(df, variable, phenotype, args):
    if variable == phenotype:
        variable = 'PC1' # set it to a dummy value -- we will skip these 
    df2 = df.rename(columns={variable:"VAR", phenotype:"PHENO"})
    if args.debug:
        print(phenotype)
        print(*df2.columns.tolist(), sep="\n")
    if args.log:
        formula="np.log(PHENO) ~ SCORE1_AVG*cov_SEX + SCORE1_AVG*cov_AGE + SCORE1_AVG*VAR + cov_SEX*cov_AGE + cov_SEX*VAR + cov_AGE*VAR"
    else:
        formula="PHENO ~ SCORE1_AVG*cov_SEX + SCORE1_AVG*cov_AGE + SCORE1_AVG*VAR + cov_SEX*cov_AGE + cov_SEX*VAR + cov_AGE*VAR"
    if "biochemistry" in phenotype:
        formula=formula+"+ biochemistry_dilutionfactor"
    if args.pc:
        for i in range(1,11):
            formula=formula+"+ PC"+str(i)+"_49K"
    if args.E2:
        formula = formula+"+ I(VAR**2)"

    if args.robust:
        model = smf.rlm(formula=formula, data = df2).fit(cov="H1")
    else:
        model = smf.ols(formula=formula, data = df2).fit()
    return(model)

def run_model_sex(df, phenotype, args):
    df2 = df.rename(columns={phenotype:"PHENO"})
    formula="PHENO ~ SCORE1_AVG*cov_SEX + SCORE1_AVG*cov_AGE + cov_SEX*cov_AGE"
    if "biochemistry" in phenotype:
        formula=formula+"+ biochemistry_dilutionfactor"
    if args.pc:
        for i in range(1,11):
            formula=formula+"+ PC"+str(i)+"_49K"
    if args.E2:
        pass
    model = smf.ols(formula=formula, data = df2).fit()
    return(model)

def run_model_sex_base(df, phenotype, args):
    df2 = df.rename(columns={phenotype:"PHENO"})
    formula="PHENO ~ SCORE1_AVG + cov_SEX + cov_AGE + cov_SEX*cov_AGE"
    if "biochemistry" in phenotype:
        formula=formula+"+ biochemistry_dilutionfactor"
    if args.pc:
        for i in range(1,11):
            formula=formula+"+ PC"+str(i)+"_49K"
    if args.E2:
        pass
    model = smf.ols(formula=formula, data = df2).fit()
    return(model)

def run_model_base(df, variable, phenotype, args):
    if variable == phenotype:
        variable = 'PC1' # set it to a dummy value -- we will skip these
    df2 = df.rename(columns={variable:"VAR", phenotype:"PHENO"})
    formula="PHENO ~ SCORE1_AVG + VAR + cov_SEX + cov_AGE + cov_SEX*cov_AGE + cov_SEX*VAR + cov_AGE*VAR"
    if "biochemistry" in phenotype:
        formula=formula+"+ biochemistry_dilutionfactor"
    if args.pc:
        for i in range(1,11):
            formula=formula+"+ PC"+str(i)+"_49K"
    if args.E2:
        formula = formula+"+ I(VAR**2)"
    model = smf.ols(formula=formula, data = df2).fit()
    return(model)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="GxE regression")
    parser.add_argument("--sscore", required=True, help="sscore file stem (everything up to the chromosome) from plink")
    parser.add_argument("--trait", required=True, help="name of trait/phenotype (y variable in regression)")
    parser.add_argument("--pc", action='store_true', help="Set if you want to correct for PCs within the 49K set")
    parser.add_argument("--E2", action='store_true', help="Set if you want to include an E^2 term")
    parser.add_argument("--sex", action='store_true', help="Set to analyze sex interaction")
    parser.add_argument("--base", action='store_true', help="Set to analyze just the base PRS model only")
    parser.add_argument("--smoking", action='store_true', help="Set to analyze Pack Years Smoking interaction")
    parser.add_argument("--log", action="store_true", help="Set to log transform dependent variable (phenotype)")
    parser.add_argument("--robust", action="store_true", help="Set to use robust linear regression")
    parser.add_argument("--debug", action='store_true', help="Set to print some debug info")
    args = parser.parse_args()

    PS = pd.read_csv(args.sscore+"1.sscore", sep="\t")
    for i in range(2,23):
        PS_temp = pd.read_csv(args.sscore+str(i)+".sscore", sep="\t")
        PS['SCORE1_AVG'] += PS_temp['SCORE1_AVG']
    
    COV = pd.read_csv('/n/groups/price/arun/ps_gxe/data/49K-pheno-cov-newE-dietpc.tab', sep="\t")
    x = PS.merge(COV, how='left', on='IID')
    q = pd.DataFrame(preprocessing.scale(x, with_mean=True, with_std=True),columns = x.columns)

    vars = [    "cov_SMOKING_STATUS",
                "cov_PHYS2_DAYS_VERY_ACTIVE",
                "cov_ALCOHOL_DESC_ORDER",
                "cov_TOWNSEND_DEPRIVATION",
                "other_TIME_TV",
                "other_NAP_TIME",
                "other_NO_WHEAT",
                "cov_Sleeplessness",
                "cov_ParticulateMatterAirPollution_irnt",
                "cov_DIET"
        ]
    if not args.robust:
        if args.sex and not args.base:
            m = run_model_sex(q, args.trait, args)
            pvals = m.pvalues[1:] # remove intercept term
            tstats = m.tvalues[1:]
            coefs = m.params[1:]
            a = anova_lm(m)
            varexp = 100 * (a['sum_sq'] / np.sum(a['sum_sq']))[:-1] # remove residuals term
            variables = a.index[1:]
            output = pd.DataFrame({'coef': coefs, 'T':tstats, 'p':pvals, 'VE':varexp})
            with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', None):
                print("#", args.trait)
                print(output)

        elif args.sex and args.base:
            m = run_model_sex_base(q, args.trait, args)
            pvals = m.pvalues[1:]
            tstats = m.tvalues[1:]
            coefs = m.params[1:]
            a = anova_lm(m)
            varexp = 100 * (a['sum_sq'] / np.sum(a['sum_sq']))[:-1] # remove residuals term
            variables = a.index[1:]
            output = pd.DataFrame({'coef': coefs, 'T':tstats, 'p':pvals, 'VE':varexp})
            with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', None):
                print("#", args.trait)
                print(output)

        elif args.base:
            for v in vars:
                m = run_model_base(q, v, args.trait, args)
                pvals = m.pvalues[1:] # remove intercept term
                tstats = m.tvalues[1:]
                coefs = m.params[1:]
                a = anova_lm(m)
                varexp = 100 * (a['sum_sq'] / np.sum(a['sum_sq']))[:-1] # remove residuals term
                variables = a.index[1:]
                output = pd.DataFrame({'coef': coefs, 'T':tstats, 'p':pvals, 'VE':varexp})
                with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', None):
                    print("#", args.trait, v)
                    print(output)
        else:
            for v in vars:
                m  = run_model(q, v, args.trait, args)
                pvals = m.pvalues[1:] # remove intercept term
                tstats = m.tvalues[1:]
                coefs = m.params[1:]
                a = anova_lm(m)
                varexp = 100 * (a['sum_sq'] / np.sum(a['sum_sq']))[:-1] # remove residuals term
                variables = a.index[1:]
                output = pd.DataFrame({'coef': coefs, 'T':tstats, 'p':pvals, 'VE':varexp})
                with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', None):
                    print("#", args.trait, v)
                    print(output)
    else:
        for v in vars:
            try:
                m  = run_model(q, v, args.trait, args)
            except np.linalg.LinAlgError:
                continue
            pvals = m.pvalues[1:] # remove intercept 
            tstats = m.tvalues[1:]
            coefs = m.params[1:]
            output = pd.DataFrame({'coef': coefs, 'T':tstats, 'p':pvals})
            with pd.option_context('display.max_rows', None, 'display.max_columns', None, 'display.width', None):
                print("#", args.trait, v)
                print(output)



