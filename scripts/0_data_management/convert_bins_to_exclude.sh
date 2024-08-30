FAM=/n/groups/price/UKBiobank/download_500K/ukb1404_cal_chr1_v2_CURRENT.fixCol6.fam

VAR=("cov_PHYS2_DAYS_VERY_ACTIVE" "cov_ALCOHOL_DESC_ORDER" "cov_TOWNSEND_DEPRIVATION" "other_TIME_TV" "other_NAP_TIME" "other_NO_WHEAT" "cov_DIET" "cov_Sleeplessness" "cov_ParticulateMatterAirPollution_irnt" "cov_SMOKING_STATUS")
for datatype in 337K 49K; do
    for bintype in 2 5; do
	for V in "${VAR[@]}"; do
	    for b in `seq 1 ${bintype}`; do
		bin="$((b-1))"
		LC_ALL=C grep -F -f bins_newE/${datatype}_${bintype}bin_${V}_${bin} -wv ${FAM} > bins_newE/exclude_${datatype}_${bintype}bin_${V}_${bin}
	    done
	done
    done
done
