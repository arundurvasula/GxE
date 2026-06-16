#while read PHENO; 
#do
#    while read C;
#    do
#	echo ${PHENO} ${C} >> data_update/trait_E.txt
#    done<data_update/E.txt
#done<data_update/traits.txt

while read PHENO;
do
    while read C;
    do
        echo ${PHENO} ${C} >> data_update/trait_newE.txt
    done<data_update/newE.txt
done<data_update/traits.txt
