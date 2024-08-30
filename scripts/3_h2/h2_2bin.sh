#!/bin/bash
#SBATCH -c 6
#SBATCH -t 5-00:00
#SBATCH -p medium
#SBATCH --mem=110G
#SBATCH -o logs/reml_%A_%a.out
#SBATCH -e logs/reml_%A_%a.err
#SBATCH --array=1-330
set -e 
set -u
I_LINE=$SLURM_ARRAY_TASK_ID

pair=$(head -n $I_LINE  data_update/trait_E.txt | tail -n 1)
PHENO=`echo $pair | cut -f 1 -d ' '`
COV=`echo $pair | cut -f 2 -d ' '`

DIR=/n/groups/price/UKBiobank/download_500K
PHENO_FILE='/n/groups/price/arun/ps_gxe/data/337K-pheno-cov.tab'
NUM_PCS=20
NBINS=1
# edited seq 0 NBINS to seq 1 NBINS to just do bin 1
#for b in `seq 0 ${NBINS}`;
for b in `seq 1 ${NBINS}`;
do
    #OUT_PREFIX=/n/groups/price/arun/ps_gxe/h2_newE/${PHENO}_${COV}_${b}
    OUT_PREFIX=/n/groups/price/arun/ps_gxe/h2_newE_moreram/${PHENO}_${COV}_${b}
    /n/groups/price/poru/HSPH_SVN/software/BOLT-LMM_v2.3.6/bolt \
                    --bed=$DIR/ukb_cal_chr{1:22}_v2.bed \
                    --bim=$DIR/ukb_snp_chr{1:22}_v2.bim \
                    --fam=$DIR/ukb1404_cal_chr1_v2_CURRENT.fixCol6.fam \
                    --remove=/n/groups/price/arun/ps_gxe/bins_newE/exclude_337K_2bin_${COV}_${b} \
                    --exclude=$DIR/../snpQC/autosome_maf_lt_0.001.txt \
                    --exclude=$DIR/../snpQC/autosome_missing_gt_0.1.txt \
                    --phenoFile=$PHENO_FILE \
                    --phenoCol=$PHENO \
                    --covarFile=/n/groups/price/arun/ps_gxe/data/337K-pheno-cov.tab \
                    --covarCol=cov_ASSESS_CENTER \
                    --covarCol=cov_GENO_ARRAY \
                    --covarCol=cov_SEX \
                    --covarMaxLevels=30 \
                    --numThreads=6 \
                    --qCovarCol=cov_AGE \
                    --qCovarCol=PC{1:$NUM_PCS} \
                    --LDscoresFile=/n/groups/price/poru/HSPH_SVN/software/BOLT-LMM_v2.3.4/tables/LDSCORE.1000G_EUR.tab.gz \
                    --geneticMapFile=/n/groups/price/poru/HSPH_SVN/software/BOLT-LMM_v2.3.4/tables/genetic_map_hg19_withX.txt.gz \
                    --reml > ${OUT_PREFIX}.log
done
