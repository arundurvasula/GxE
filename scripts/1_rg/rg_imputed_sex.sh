#!/bin/bash
#SBATCH -c 1
#SBATCH -t 00-00:05
#SBATCH -p short
#SBATCH --mem=8G
#SBATCH -o logs/ldsc_%A_%a.out
#SBATCH -e logs/ldsc_%A_%a.err
#SBATCH --array=1-33

set -e 
I_LINE=$SLURM_ARRAY_TASK_ID
NAME=$(head -n $I_LINE data_update/traits.txt | tail -n 1)
NAMEARRAY=($NAME)
TRAIT=${NAMEARRAY[0]}

SCRATCH=/n/scratch/users/a/ard063/ldsc_sex
mkdir -p $SCRATCH

#source ~/.bashrc
#source activate ldsc

for BIN in 0 1; 
do
    SUMSTAT=/n/groups/price/arun/ps_gxe/gwas_337K_sex/${TRAIT}_${BIN}.bgen.stats.gz
    /home/ard063/.conda/envs/ldsc/bin/python ldsc/munge_sumstats.py \
        --sumstats ${SUMSTAT} \
        --out ${SCRATCH}/${TRAIT}_${BIN} \
        --merge-alleles data/w_hm3.snplist \
        --p P_LINREG \
        --N 168653 \
        --a1 ALLELE1 \
        --a2 ALLELE0 \
	--chunksize 500000 \
        --signed-sumstats BETA,0 \
        --frq A1FREQ
done

set -e
#COMBS=(01 02 03 04 12 13 14 23 24 34)
COMBS=(01)
for C in ${COMBS[@]}; do
    BIN1=${C:0:1}
    BIN2=${C:1:1}
    SUMSTAT1=${SCRATCH}/${TRAIT}_${BIN1}.sumstats.gz
    SUMSTAT2=${SCRATCH}/${TRAIT}_${BIN2}.sumstats.gz
    /home/ard063/.conda/envs/ldsc/bin/python ldsc/ldsc.py \
        --rg ${SUMSTAT1},${SUMSTAT2} \
        --ref-ld-chr data/eur_w_ld_chr/ \
        --w-ld-chr data/eur_w_ld_chr/ \
	--no-intercept \
        --out rg_imputed_sex2/${TRAIT}_${BIN1}_${BIN2}_intercept
done
