#!/bin/bash

cp ../../raw/retrospective_ova_phospho_sort_common_gene_10057.txt ./raw_ova_phospho.txt
cp ../../raw/retrospective_ova_PNNL_proteome_sort_common_gene_7061.txt ./raw_ova_proteome.txt
cp ../../raw/retrospective_ova_rna_seq_sort_common_gene_15121.txt ./raw_ova_rna.txt
cp ../../raw/retrospective_ova_CNA_sort_common_gene_11859.txt ./raw_ova_cnv.txt

cp ../../raw/retrospective_breast_phospho_sort_common_gene_31981.txt ./raw_breast_phospho.txt
cp ../../raw/retrospective_breast_proteome_sort_common_gene_10005.txt ./raw_breast_proteome.txt
cp ../../raw/retrospective_breast_RNA_sort_common_gene_15107.txt ./raw_breast_rna.txt
cp ../../raw/retrospective_breast_CNA_sort_common_gene_16884.txt ./raw_breast_cnv.txt

# normalization
Rscript ../../../../function/normalize_by_sample_avg_sd.r raw*txt


##### 2. global branch #####
# trim
####### This part rna imputation #######
Rscript trim_data.r
########################################

# imputation
Rscript ../../../../function/avg_imputation.r trimmed*scaled
###########################




