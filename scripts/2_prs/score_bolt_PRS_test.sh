#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-12:00
#SBATCH -p short
#SBATCH --mem=16G
#SBATCH -o logs/score_bolt_prs_%A_%a.out
#SBATCH -e logs/score_bolt_prs_%A_%a.err
#SBATCH --array=1-22

CHR=$SLURM_ARRAY_TASK_ID
BED=/n/groups/price/UKBiobank/download_500K/ukb_cal_chr${CHR}_v2.bed
FAM=/n/groups/price/UKBiobank/download_500K/ukb_cal_chr1_v2.fam
BIM=/n/groups/price/UKBiobank/download_500K/ukb_snp_chr${CHR}_v2.bim
MAF=0.001

TRAIT=biochemistry_LDLdirect
COV=other_TIME_TV
SUMSTATS=bolt_prs_test/${TRAIT}_${COV}-combined.beta
SAMPLES='sample_lists/49K.txt'
/home/ard063/software/plink2 --bed ${BED} --fam ${FAM} --bim ${BIM} --maf ${MAF} --keep ${SAMPLES} --score ${SUMSTATS} 'header-read' --score-col-nums 3-14 --out bolt_prs_test/score/${TRAIT}_${COV}_${CHR}


TRAIT=disease_T2D
COV=cov_PHYS2_DAYS_VERY_ACTIVE
SUMSTATS=bolt_prs_test/${TRAIT}_${COV}-combined.beta
SAMPLES='sample_lists/49K.txt'
/home/ard063/software/plink2 --bed ${BED} --fam ${FAM} --bim ${BIM} --maf ${MAF} --keep ${SAMPLES} --score ${SUMSTATS} 'header-read' --score-col-nums 3-14 --out bolt_prs_test/score/${TRAIT}_${COV}_${CHR}
