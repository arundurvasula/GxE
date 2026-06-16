#!/bin/bash
#SBATCH -c 1
#SBATCH -t 00-12:00
#SBATCH -p short
#SBATCH --mem=8G
#SBATCH -o logs/sim%A_%a.out
#SBATCH -e logs/sim%A_%a.err
#SBATCH --array=1-1800

set -e 
set -u

PARAM_FILE=../simulations_GxEMM/params/all.txt
I_LINE=$SLURM_ARRAY_TASK_ID
params=$(head -n $I_LINE  $PARAM_FILE | tail -n 1)
SEED=`echo $params | cut -f 1 -d ' '`
h21=`echo $params | cut -f 2 -d ' '`
h22=`echo $params | cut -f 3 -d ' '`
gen_corr=`echo $params | cut -f 4 -d ' '`
amp=`echo $params | cut -f 5 -d ' '`
M=`echo $params | cut -f 6 -d ' '`
N=`echo $params | cut -f 7 -d ' '`
SCENARIO=`echo $params | cut -f 8 -d ' '`

sleep $((RANDOM%1800+1))
module load gcc/6.2.0
module load R

Rscript scripts/simulate.R $gen_corr $h21 $h22 $amp scenario${SCENARIO}_$SEED-$gen_corr-$h21-$h22-$amp
