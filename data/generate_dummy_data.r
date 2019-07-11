
set.seed(449)

path1='./raw_ori/'
path2='./raw/'

file_all=list.files(path1)

for (the_file in file_all){
    print(the_file)
    ori=as.matrix(read.delim(paste0(path1,the_file),check.names=F,row.names=1))
    # add random gaussian noise
    ori = ori + rnorm(length(ori))
    # shuffle
    the_name=rownames(ori)
    ori=ori[sample.int(nrow(ori)),]
    rownames(ori)=the_name
    # reduce digit size to 3 to save space..
    ori=round(ori,3)
    # save
    tmp=ori;tmp=cbind(rownames(tmp),tmp);colnames(tmp)[1]="Gene_ID"
    write.table(tmp,file=paste0(path2,the_file),quote=F,sep="\t",row.names=F,na='')
}


