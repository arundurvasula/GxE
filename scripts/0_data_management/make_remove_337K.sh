DIR=/n/groups/price/UKBiobank/download_500K

cat $DIR/bolt.in_plink_but_not_imputed.FID_IID.976.txt \
$DIR/../sampleQC/remove.nonStringentBritish.FID_IID.txt \
$DIR/../sampleQC/remove.related.FID_IID.txt > data_update/remove_337K.txt
