## name: create_top_feature_list.r
## date: 11/05/2017

path1="../raw/"

##### 1. ova #################
proteome=as.matrix(read.delim(paste0(path1,"retrospective_ova_PNNL_proteome_sort_common_gene_7061.txt"),check.names=F,row.names=1))

cutoff=dim(proteome)[2]*1.0
ind=apply(proteome,1,function(x){sum(!is.na(x))>=cutoff})
proteome=proteome[ind,]

# top expressed
avg=apply(proteome,1,mean,na.rm=T)
writeLines(sort(names(sort(avg,decreasing=T)[1:1000])),con="list_ova_proteome_top1000.txt")
writeLines(sort(names(sort(avg,decreasing=T)[1:100])),con="list_ova_proteome_top100.txt")
writeLines(sort(names(sort(avg,decreasing=T)[1:10])),con="list_ova_proteome_top10.txt")

##################################

##### 2. breast #################
proteome=as.matrix(read.delim(paste0(path1,"retrospective_breast_proteome_sort_common_gene_10005.txt"),check.names=F,row.names=1))

cutoff=dim(proteome)[2]*1.0
ind=apply(proteome,1,function(x){sum(!is.na(x))>=cutoff})
proteome=proteome[ind,]

# top expressed
avg=apply(proteome,1,mean,na.rm=T)
writeLines(sort(names(sort(avg,decreasing=T)[1:1000])),con="list_breast_proteome_top1000.txt")
writeLines(sort(names(sort(avg,decreasing=T)[1:100])),con="list_breast_proteome_top100.txt")
writeLines(sort(names(sort(avg,decreasing=T)[1:10])),con="list_breast_proteome_top10.txt")



