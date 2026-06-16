#!/bin/bash
#SBATCH -c 1
#SBATCH -t 0-12:00
#SBATCH -p short
#SBATCH --mem=4G
#SBATCH -o logs/reg_robust.out
#SBATCH -e logs/reg_robust.err

source activate py39

set -e
while read TRAIT; 
do
    echo "Analyzing "$TRAIT
#    python scripts_january_push/reg_update.py --sscore results/${TRAIT}_ --trait ${TRAIT} --E2 --pc --robust > regression_update/${TRAIT}_robust.txt
#    python scripts_january_push/reg_update.py --sscore results/${TRAIT}_ --trait ${TRAIT} --E2 --pc > regression_update/${TRAIT}.txt
#    python scripts_january_push/reg_update.py --sscore results/${TRAIT}_ --trait ${TRAIT} --E2 --pc --log > regression_update/${TRAIT}_log.txt
#    python scripts_january_push/reg_update.py --sscore results/${TRAIT}_ --trait ${TRAIT} --pc --E2 --base > regression_update/${TRAIT}_base.txt
#    python scripts_january_push/reg_update.py --sscore results/${TRAIT}_ --trait ${TRAIT} --pc --E2 --base --sex > regression_update/${TRAIT}_sex_base.txt
    python scripts_january_push/reg_update.py --sscore results/${TRAIT}_ --trait ${TRAIT} --pc > regression_update/${TRAIT}_noE2.txt
done<data_update/traits.txt
