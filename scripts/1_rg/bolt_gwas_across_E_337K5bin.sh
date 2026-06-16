#!/bin/bash
#SBATCH -c 6
#SBATCH -t 04-00:00
#SBATCH -p medium
#SBATCH --mem=60G
#SBATCH -o logs/bin5gwas337imp_%A_%a.out
#SBATCH -e logs/bin5gwas337imp_%A_%a.err
#SBATCH --array=1-132

set -e 
set -u
I_LINE=$SLURM_ARRAY_TASK_ID
pair=$(head -n $I_LINE  data_update/trait_newE.txt | tail -n 1)
PHENO=`echo $pair | cut -f 1 -d ' '`
COV=`echo $pair | cut -f 2 -d ' '`

THREADS=6

DIR=/n/groups/price/UKBiobank/download_500K
PHENO_FILE='/n/groups/price/arun/ps_gxe/data/337K-pheno-cov.tab'
NUM_PCS=20
NBINS=4
BINTYPE="$((NBINS+1))"
for b in `seq 0 ${NBINS}`;
do
    OUT_PREFIX=/n/scratch/users/a/ard063/gwas_337K_5bin/${PHENO}_${COV}_${b}
    /n/groups/price/poru/HSPH_SVN/software/BOLT-LMM_v2.3.6/bolt \
                    --bed=$DIR/ukb_cal_chr{1:22}_v2.bed \
                    --bim=$DIR/ukb_snp_chr{1:22}_v2.bim \
                    --fam=$DIR/ukb1404_cal_chr1_v2_CURRENT.fixCol6.fam \
                    --remove=/n/groups/price/arun/ps_gxe/bins_newE/exclude_337K_${BINTYPE}bin_${COV}_${b} \
                    --exclude=$DIR/../snpQC/autosome_maf_lt_0.001.txt \
                    --exclude=$DIR/../snpQC/autosome_missing_gt_0.1.txt \
                    --phenoFile=$PHENO_FILE \
                    --phenoCol=$PHENO \
                    --covarFile=/n/groups/price/arun/ps_gxe/data/337K-pheno-cov.tab \
                    --covarCol=cov_ASSESS_CENTER \
                    --covarCol=cov_GENO_ARRAY \
                    --covarCol=cov_SEX \
                    --covarMaxLevels=30 \
                    --qCovarCol=cov_AGE \
                    --qCovarCol=PC{1:$NUM_PCS} \
                    --LDscoresFile=/n/groups/price/poru/HSPH_SVN/software/BOLT-LMM_v2.3.4/tables/LDSCORE.1000G_EUR.tab.gz \
                    --geneticMapFile=/n/groups/price/poru/HSPH_SVN/software/BOLT-LMM_v2.3.4/tables/genetic_map_hg19_withX.txt.gz \
	            --numThreads=$THREADS \
                    --statsFile=$OUT_PREFIX.stats.gz \
	            --bgenFile=/n/groups/price/UKBiobank/bgen_MAF001_500K_v3/UKB_MAF0.001_v3.{1:22}.bgen \
	            --sampleFile=$DIR/ukb1404_imp_chr1_v2_s487406.sample \
	            --statsFileBgenSnps=$OUT_PREFIX.bgen.stats.gz \
	            --verboseStats \
                    --lmmInfOnly > ${OUT_PREFIX}.log
rm $OUT_PREFIX.stats.gz
done
