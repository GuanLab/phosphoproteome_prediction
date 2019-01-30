## name: trim_data.r
## date: 10/09/2017

## Here I trim the raw dataset into trimmed_set which only contain samples with proteome, rna-seq and CNA results

path1="../raw/";
path2="./";

############ ova ###########
phospho=as.matrix(read.delim(paste0(path1,"retrospective_ova_phospho_sort_common_gene_10057.txt"),check.names=F,row.names=1))
proteome=as.matrix(read.delim(paste0(path1,"retrospective_ova_PNNL_proteome_sort_common_gene_7061.txt"),check.names=F,row.names=1))
rna=as.matrix(read.delim(paste0(path1,"retrospective_ova_rna_seq_sort_common_gene_15121.txt"),check.names=F,row.names=1))
cnv=as.matrix(read.delim(paste0(path1,"retrospective_ova_CNA_sort_common_gene_11859.txt"),check.names=F,row.names=1))
dim(phospho);dim(proteome);dim(rna);dim(cnv)
#[1] 10057    69
#[1] 7061   84
#[1] 15121   294
#[1] 11859   559
phospho_old=as.matrix(read.delim(paste0(path1,"retrospective_ova_phospho_filtered.txt"),check.names=F,row.names=1))
dim(phospho_old)
#[1] 21332    69  This old version has more sites, which is used for transplant ova to breast

# write
tmp=phospho;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"ova_phospho.txt"),quote=F,sep="\t",row.names=F)
tmp=phospho_old;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"ova_phospho_old_version.txt"),quote=F,sep="\t",row.names=F)

# proteome/rna/cnv
proteome=proteome[,colnames(proteome) %in% colnames(phospho)]
rna=rna[,colnames(rna) %in% colnames(phospho)]
cnv=cnv[,colnames(cnv) %in% colnames(phospho)]

proteome_new=matrix(NA,nrow=dim(phospho)[1],ncol=dim(proteome)[2])
rna_new=matrix(NA,nrow=dim(phospho)[1],ncol=dim(rna)[2])
cnv_new=matrix(NA,nrow=dim(phospho)[1],ncol=dim(cnv)[2])
rownames(proteome_new)=rownames(rna_new)=rownames(cnv_new)=rownames(phospho)
colnames(proteome_new)=colnames(proteome)
colnames(rna_new)=colnames(rna)
colnames(cnv_new)=colnames(cnv)
dim(proteome_new);dim(rna_new);dim(cnv_new)
#[1] 10057    69
#[1] 10057    47
#[1] 10057    69

# map
id2name=sub("\\..*","",rownames(phospho))

for(i in 1:dim(phospho)[1]){
        if(id2name[i] %in% rownames(proteome)){
                proteome_new[i,]=proteome[id2name[i],]
        }
        if(id2name[i] %in% rownames(rna)){
                rna_new[i,]=rna[id2name[i],]
        }
        if(id2name[i] %in% rownames(cnv)){
                cnv_new[i,]=cnv[id2name[i],]
        }
}

# write
tmp=proteome_new;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"trimmed_ova_proteome.txt"),quote=F,sep="\t",row.names=F)
tmp=rna_new;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"trimmed_ova_rna.txt"),quote=F,sep="\t",row.names=F)
tmp=cnv_new;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"trimmed_ova_cnv.txt"),quote=F,sep="\t",row.names=F)
#################################################

############ breast ###########
phospho=as.matrix(read.delim(paste0(path1,"retrospective_breast_phospho_sort_common_gene_31981.txt"),check.names=F,row.names=1))
proteome=as.matrix(read.delim(paste0(path1,"retrospective_breast_proteome_sort_common_gene_10005.txt"),check.names=F,row.names=1))
rna=as.matrix(read.delim(paste0(path1,"retrospective_breast_RNA_sort_common_gene_15107.txt"),check.names=F,row.names=1))
cnv=as.matrix(read.delim(paste0(path1,"retrospective_breast_CNA_sort_common_gene_16884.txt"),check.names=F,row.names=1))
dim(phospho);dim(proteome);dim(rna);dim(cnv)
#[1] 31981   105
#[1] 10006   105
#[1] 15107    77
#[1] 16884    77

# write
tmp=phospho;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"breast_phospho.txt"),quote=F,sep="\t",row.names=F)

# proteome/rna/cnv
proteome=proteome[,colnames(proteome) %in% colnames(phospho)]
rna=rna[,colnames(rna) %in% colnames(phospho)]
cnv=cnv[,colnames(cnv) %in% colnames(phospho)]
dim(proteome);dim(rna);dim(cnv)
#[1] 10006   105
#[1] 15107    77
#[1] 16884    77

proteome_new=matrix(NA,nrow=dim(phospho)[1],ncol=dim(proteome)[2])
rna_new=matrix(NA,nrow=dim(phospho)[1],ncol=dim(rna)[2])
cnv_new=matrix(NA,nrow=dim(phospho)[1],ncol=dim(cnv)[2])
rownames(proteome_new)=rownames(rna_new)=rownames(cnv_new)=rownames(phospho)
colnames(proteome_new)=colnames(proteome)
colnames(rna_new)=colnames(rna)
colnames(cnv_new)=colnames(cnv)
dim(proteome_new);dim(rna_new);dim(cnv_new)
#[1] 31981   105
#[1] 31981    77
#[1] 31981    77

# map
id2name=sub("\\..*","",rownames(phospho))

for(i in 1:dim(phospho)[1]){
        if(id2name[i] %in% rownames(proteome)){
                proteome_new[i,]=proteome[id2name[i],]
        }
        if(id2name[i] %in% rownames(rna)){
                rna_new[i,]=rna[id2name[i],]
        }
        if(id2name[i] %in% rownames(cnv)){
                cnv_new[i,]=cnv[id2name[i],]
        }
}       

# write
tmp=proteome_new;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"trimmed_breast_proteome.txt"),quote=F,sep="\t",row.names=F)
tmp=rna_new;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"trimmed_breast_rna.txt"),quote=F,sep="\t",row.names=F)
tmp=cnv_new;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
write.table(tmp,file=paste0(path2,"trimmed_breast_cnv.txt"),quote=F,sep="\t",row.names=F)
#################################################


