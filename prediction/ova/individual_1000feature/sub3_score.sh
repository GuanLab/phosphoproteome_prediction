#!/bin/bash

path0="../../../data/cv_set/"
path1="./" # PARAMETER:

for i in `seq 1 5`;
do
        Rscript sub3_hongyang.r ${path0}test_ova_phospho_${i}.txt ${path1}prediction_ova_phospho_${i}.txt >> cor_avg_ova.txt
	mv cor_nrmse.txt cor_nrmse_o${i}.txt
done  
echo ova
cat cor_avg_ova.txt

