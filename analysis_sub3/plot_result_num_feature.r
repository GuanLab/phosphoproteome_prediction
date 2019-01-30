
library(ggplot2)
library(reshape2)
source("~/function/my_palette.r")
source("~/function/multiplot.R")

cor1=rmse1=array(NA,dim=c(31981,5,5))
cor2=rmse2=array(NA,dim=c(10057,5,5))
for(i in 1:5){
    tmp=read.delim(paste0("../sub3/prediction/breast/individual_10feature/cor_nrmse_b",i,".txt"),header=F,row.names=1)
    cor1[,1,i]=tmp[,1]
    rmse1[,1,i]=tmp[,2]
    tmp=read.delim(paste0("../sub3/prediction/breast/individual_100feature/cor_nrmse_b",i,".txt"),header=F,row.names=1)
    cor1[,2,i]=tmp[,1]
    rmse1[,2,i]=tmp[,2]
    tmp=read.delim(paste0("../sub3/prediction/breast/individual_1000feature/cor_nrmse_b",i,".txt"),header=F,row.names=1)
    cor1[,3,i]=tmp[,1]
    rmse1[,3,i]=tmp[,2]
    tmp=read.delim(paste0("../sub3/prediction/breast/individual_go_phosphorylation/cor_nrmse_b",i,".txt"),header=F,row.names=1)
    cor1[,4,i]=tmp[,1]
    rmse1[,4,i]=tmp[,2]
    tmp=read.delim(paste0("../sub3/prediction/breast/individual/cor_nrmse_b",i,".txt"),header=F,row.names=1)
    cor1[,5,i]=tmp[,1]
    rmse1[,5,i]=tmp[,2]
}
for(i in 1:5){
    tmp=read.delim(paste0("../sub3/prediction/ova/individual_10feature/cor_nrmse_o",i,".txt"),header=F,row.names=1)
    cor2[,1,i]=tmp[,1]
    rmse2[,1,i]=tmp[,2]
    tmp=read.delim(paste0("../sub3/prediction/ova/individual_100feature/cor_nrmse_o",i,".txt"),header=F,row.names=1)
    cor2[,2,i]=tmp[,1]
    rmse2[,2,i]=tmp[,2]
    tmp=read.delim(paste0("../sub3/prediction/ova/individual_1000feature/cor_nrmse_o",i,".txt"),header=F,row.names=1)
    cor2[,3,i]=tmp[,1]
    rmse2[,3,i]=tmp[,2]
    tmp=read.delim(paste0("../sub3/prediction/ova/individual_go_phosphorylation/cor_nrmse_o",i,".txt"),header=F,row.names=1)
    cor2[,4,i]=tmp[,1]
    rmse2[,4,i]=tmp[,2]
    tmp=read.delim(paste0("../sub3/prediction/ova/individual/cor_nrmse_o",i,".txt"),header=F,row.names=1)
    cor2[,5,i]=tmp[,1]
    rmse2[,5,i]=tmp[,2]
}

name2=c("10","100","1000","1476(GO)","6956(all)")

# breast
dat=apply(cor1,c(1:2),mean,na.rm=T)
set.seed(449)
perm=NULL
for(i in 1:50){ # calculate p
    perm=rbind(perm,apply(dat[sample(nrow(dat),1000),],2,mean,na.rm=T))
}
pval=c(wilcox.test(perm[,1],perm[,2])$p.value,wilcox.test(perm[,2],perm[,3])$p.value,
    wilcox.test(perm[,3],perm[,4])$p.value,wilcox.test(perm[,4],perm[,5])$p.value)
dat=rbind(apply(perm,2,mean,na.rm=T), perm)
colnames(dat)=name2
dat.m=melt(dat)
colnames(dat.m)=c("rep","model","correlation")
tmp_col=c(p_jco[c(2,7,5,6,9)])
names(tmp_col)=colnames(dat)

p1 = ggplot(dat.m, aes(model, correlation, fill=model)) + 
    geom_violin(width=0.5, colour="grey50",draw_quantiles = 0.5) +
    scale_fill_manual(values=tmp_col) +
    labs(x="number of features", y="Pearson's correlation",title="breast") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),legend.position="none") +
    #ylim(0.1,0.35) +
    theme(axis.line.x = element_line(color="black", size = 0.5))+
    theme(axis.line.y = element_line(color="black", size = 0.5))+
    theme(axis.title.x = element_text(colour="black", size=12))+
    theme(axis.title.y = element_text(colour="black", size=12))+
    theme(axis.text.x = element_text(colour="black",size=12))+
    theme(axis.text.y = element_text(colour="black",size=12))+
    annotate("text", x = seq(1.5,4.5,1), y = 0.59, label=paste0("p=",format(pval,digits=3)),size=4)+
    annotate("text", x = 1:5, y = 0.645, label=paste0(format(dat[1,],digits=3)),size=5)
#    geom_rect(mapping=aes(xmin=2.5,xmax=3.5,ymin=0.25,ymax=0.45,fill=NA),color=my_palette["yellow"])
p1


dat=apply(rmse1,c(1:2),mean,na.rm=T)
set.seed(449)
perm=NULL
for(i in 1:50){ # calculate p
    perm=rbind(perm,apply(dat[sample(nrow(dat),1000),],2,mean,na.rm=T))
}
pval=c(wilcox.test(perm[,1],perm[,2])$p.value,wilcox.test(perm[,2],perm[,3])$p.value,
    wilcox.test(perm[,3],perm[,4])$p.value,wilcox.test(perm[,4],perm[,5])$p.value)
dat=rbind(apply(perm,2,mean,na.rm=T), perm)
colnames(dat)=name2
dat.m=melt(dat)
colnames(dat.m)=c("rep","model","correlation")
tmp_col=c(p_jco[c(2,7,5,6,9)])
names(tmp_col)=colnames(dat)

p2 = ggplot(dat.m, aes(model, correlation, fill=model)) + 
    geom_violin(width=0.5, colour="grey50",draw_quantiles = 0.5) +
    scale_fill_manual(values=tmp_col) +
    labs(x="number of features", y="NRMSE",title="breast") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),legend.position="none") +
    #ylim(0.1,0.35) +
    theme(axis.line.x = element_line(color="black", size = 0.5))+
    theme(axis.line.y = element_line(color="black", size = 0.5))+
    theme(axis.title.x = element_text(colour="black", size=12))+
    theme(axis.title.y = element_text(colour="black", size=12))+
    theme(axis.text.x = element_text(colour="black",size=12))+
    theme(axis.text.y = element_text(colour="black",size=12))+
    annotate("text", x = seq(1.5,4.5,1), y = 0.27, label=paste0("p=",format(pval,digits=3)),size=4)+
    annotate("text", x = 1:5, y = 0.29, label=paste0(format(dat[1,],digits=3)),size=5)
#    geom_rect(mapping=aes(xmin=2.5,xmax=3.5,ymin=0.25,ymax=0.45,fill=NA),color=my_palette["yellow"])
p2


# ovary
dat=apply(cor2,c(1:2),mean,na.rm=T)
set.seed(449)
perm=NULL
for(i in 1:50){ # calculate p
    perm=rbind(perm,apply(dat[sample(nrow(dat),1000),],2,mean,na.rm=T))
}
pval=c(wilcox.test(perm[,1],perm[,2])$p.value,wilcox.test(perm[,2],perm[,3])$p.value,
    wilcox.test(perm[,3],perm[,4])$p.value,wilcox.test(perm[,4],perm[,5])$p.value)
dat=rbind(apply(perm,2,mean,na.rm=T), perm)
colnames(dat)=name2
dat.m=melt(dat)
colnames(dat.m)=c("rep","model","correlation")
tmp_col=c(p_jco[c(2,7,5,6,9)])
names(tmp_col)=colnames(dat)

p3 = ggplot(dat.m, aes(model, correlation, fill=model)) + 
    geom_violin(width=0.5, colour="grey50",draw_quantiles = 0.5) +
    scale_fill_manual(values=tmp_col) +
    labs(x="number of features", y="Pearson's correlation",title="ovary") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),legend.position="none") +
    #ylim(0.1,0.35) +
    theme(axis.line.x = element_line(color="black", size = 0.5))+
    theme(axis.line.y = element_line(color="black", size = 0.5))+
    theme(axis.title.x = element_text(colour="black", size=12))+
    theme(axis.title.y = element_text(colour="black", size=12))+
    theme(axis.text.x = element_text(colour="black",size=12))+
    theme(axis.text.y = element_text(colour="black",size=12))+
    annotate("text", x = seq(1.5,4.5,1), y = 0.27, label=paste0("p=",format(pval,digits=3)),size=4)+
    annotate("text", x = 1:5, y = 0.37, label=paste0(format(dat[1,],digits=3)),size=5)
#    geom_rect(mapping=aes(xmin=2.5,xmax=3.5,ymin=0.25,ymax=0.45,fill=NA),color=my_palette["yellow"])
p3


dat=apply(rmse2,c(1:2),mean,na.rm=T)
set.seed(449)
perm=NULL
for(i in 1:50){ # calculate p
    perm=rbind(perm,apply(dat[sample(nrow(dat),1000),],2,mean,na.rm=T))
}
pval=c(wilcox.test(perm[,1],perm[,2])$p.value,wilcox.test(perm[,2],perm[,3])$p.value,
    wilcox.test(perm[,3],perm[,4])$p.value,wilcox.test(perm[,4],perm[,5])$p.value)
dat=rbind(apply(perm,2,mean,na.rm=T), perm)
colnames(dat)=name2
dat.m=melt(dat)
colnames(dat.m)=c("rep","model","correlation")
tmp_col=c(p_jco[c(2,7,5,6,9)])
names(tmp_col)=colnames(dat)

p4 = ggplot(dat.m, aes(model, correlation, fill=model)) + 
    geom_violin(width=0.5, colour="grey50",draw_quantiles = 0.5) +
    scale_fill_manual(values=tmp_col) +
    labs(x="number of features", y="NRMSE",title="ovary") +
    theme_light() +
    theme(plot.title = element_text(hjust = 0.5),legend.position="none") +
    #ylim(0.1,0.35) +
    theme(axis.line.x = element_line(color="black", size = 0.5))+
    theme(axis.line.y = element_line(color="black", size = 0.5))+
    theme(axis.title.x = element_text(colour="black", size=12))+
    theme(axis.title.y = element_text(colour="black", size=12))+
    theme(axis.text.x = element_text(colour="black",size=12))+
    theme(axis.text.y = element_text(colour="black",size=12))+
    annotate("text", x = seq(1.5,4.5,1), y = 0.37, label=paste0("p=",format(pval,digits=3)),size=4)+
    annotate("text", x = 1:5, y = 0.4, label=paste0(format(dat[1,],digits=3)),size=5)
#    geom_rect(mapping=aes(xmin=2.5,xmax=3.5,ymin=0.25,ymax=0.45,fill=NA),color=my_palette["yellow"])
p4


pdf(file="figure/result_num_feature.pdf",width=10,height=6,useDingbats=F)
list_p=c(list(p1),list(p2),list(p3),list(p4))
mat_layout=matrix(1:4,nrow=2,byrow=F)
multiplot(plotlist=list_p,layout = mat_layout)
dev.off()

# table
name1=paste0("cross-validation ",1:5)
tbl=apply(cor1,3:2,mean,na.rm=T)
rownames(tbl)=name1; colnames(tbl)=name2
write.csv(tbl,file="table/result_num_feature_cor1.csv")

tbl=apply(cor2,3:2,mean,na.rm=T)
rownames(tbl)=name1; colnames(tbl)=name2
write.csv(tbl,file="table/result_num_feature_cor2.csv")

tbl=apply(rmse1,3:2,mean,na.rm=T)
rownames(tbl)=name1; colnames(tbl)=name2
write.csv(tbl,file="table/result_num_feature_rmse1.csv")

tbl=apply(rmse2,3:2,mean,na.rm=T)
rownames(tbl)=name1; colnames(tbl)=name2
write.csv(tbl,file="table/result_num_feature_rmse2.csv")



