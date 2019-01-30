## name: customize_multisite.r
## date: 10/13/2017

path1="../individual_transplant/"
path2="./"

mat_num=as.matrix(read.delim("ova_num_sample.txt",row.names=1,header=F))

for (j in 1:5){
	input=as.matrix(read.delim(paste0(path1,"prediction_ova_phospho_",j,".txt"),check.names=F,row.names=1))
	output=input
	
	site=rownames(input)
	gene=sub("\\..*","",site)
	gene_uniq=unique(gene)
	for(i in gene_uniq){
	        if(sum(gene==i)>1){ # sites with more observations have larger weights
	                mat_weight=matrix(rep(mat_num[gene==i,],times=dim(input)[2]),nrow=sum(gene==i))/sum(mat_num[gene==i,])
	                output[gene==i,]=rep(apply(input[gene==i,]*mat_weight,2,sum),each=sum(gene==i))
	        }
	}
	
	tmp=output;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="phosphoID"
	write.table(tmp,file=paste0(path2,"prediction_oo_phospho_",j,".txt"),quote=F,sep="\t",row.names=F)
}
