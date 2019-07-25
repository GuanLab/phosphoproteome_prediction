## name: stack_model.r
## date: 08/15/2017


path1="/state3/hyangl/CPTAC_final/sub3/prediction/m175/"
path2="./"

for (w2 in c(1/seq(3,1,-0.5),seq(1.5,3,0.5))) {
    w1=1
    w=w1+w2

    for (i in 0:0){
        x1=as.matrix(read.delim(paste0(path1,"prediction_breast_phospho_",i,".txt"),check.names=F,row.names=1))
        x2=as.matrix(read.delim(paste0(path2,"prediction_ob_phospho_",i,".txt"),check.names=F,row.names=1))
    
        x2[is.na(x2)]=x1[is.na(x2)]
        y=w1/w*x1 + w2/w*x2

        tmp=cbind(rownames(y),y);colnames(tmp)[1]="Gene_ID"
        write.table(tmp,file=paste0(path2, "prediction_breast_phospho_",i,".txt"),quote=F,sep="\t",row.names=F)
    }
    system("./quick_breast.sh")
}

    system("mv cor_avg_breast.txt scan_cor_breast.txt")

