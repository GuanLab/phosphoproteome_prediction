## name: create_subset1.0_list.r
## date: 10/09/2017

path1="../raw/"

proteome=as.matrix(read.delim(paste0(path1,"retrospective_ova_PNNL_proteome_sort_common_gene_7061.txt"),check.names=F,row.names=1))
cutoff=dim(proteome)[2]*1.0
ind=apply(proteome,1,function(x){sum(!is.na(x))>=cutoff})
writeLines(rownames(proteome)[ind],con="list_ova_proteome_subset1.0.txt")

proteome=as.matrix(read.delim(paste0(path1,"retrospective_breast_proteome_sort_common_gene_10005.txt"),check.names=F,row.names=1))
cutoff=dim(proteome)[2]*1.0
ind=apply(proteome,1,function(x){sum(!is.na(x))>=cutoff})
writeLines(rownames(proteome)[ind],con="list_breast_proteome_subset1.0.txt")



