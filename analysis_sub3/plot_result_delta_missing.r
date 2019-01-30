library(ggplot2)
library(reshape2)
source("~/function/my_palette.r")
source("~/function/multiplot.R")

# delta-correlation (improvement after multisite)

b1=read.delim("../sub3/prediction/breast/multisite/cor_nrmse_b1.txt",header=F,row.names=1)
b2=read.delim("../sub3/prediction/breast/multisite/cor_nrmse_b2.txt",header=F,row.names=1)
b3=read.delim("../sub3/prediction/breast/multisite/cor_nrmse_b3.txt",header=F,row.names=1)
b4=read.delim("../sub3/prediction/breast/multisite/cor_nrmse_b4.txt",header=F,row.names=1)
b5=read.delim("../sub3/prediction/breast/multisite/cor_nrmse_b5.txt",header=F,row.names=1)
mat1=cbind(b1[,1],b2[,1],b3[,1],b4[,1],b5[,1])
mat3=cbind(b1[,2],b2[,2],b3[,2],b4[,2],b5[,2])
b1=read.delim("../sub3/prediction/breast/individual_transplant/cor_nrmse_b1.txt",header=F,row.names=1)
b2=read.delim("../sub3/prediction/breast/individual_transplant/cor_nrmse_b2.txt",header=F,row.names=1)
b3=read.delim("../sub3/prediction/breast/individual_transplant/cor_nrmse_b3.txt",header=F,row.names=1)
b4=read.delim("../sub3/prediction/breast/individual_transplant/cor_nrmse_b4.txt",header=F,row.names=1)
b5=read.delim("../sub3/prediction/breast/individual_transplant/cor_nrmse_b5.txt",header=F,row.names=1)
mat2=cbind(b1[,1],b2[,1],b3[,1],b4[,1],b5[,1])
mat4=cbind(b1[,2],b2[,2],b3[,2],b4[,2],b5[,2])
cor1=apply(mat1-mat2,1,mean,na.rm=T)
rmse1=apply(mat3-mat4,1,mean,na.rm=T)
mean(cor1,na.rm=T) #[1] 0.01384422
mean(rmse1,na.rm=T) #[1] -0.001436225


o1=read.delim("../sub3/prediction/ova/multisite/cor_nrmse_o1.txt",header=F,row.names=1)
o2=read.delim("../sub3/prediction/ova/multisite/cor_nrmse_o2.txt",header=F,row.names=1)
o3=read.delim("../sub3/prediction/ova/multisite/cor_nrmse_o3.txt",header=F,row.names=1)
o4=read.delim("../sub3/prediction/ova/multisite/cor_nrmse_o4.txt",header=F,row.names=1)
o5=read.delim("../sub3/prediction/ova/multisite/cor_nrmse_o5.txt",header=F,row.names=1)
mat1=cbind(o1[,1],o2[,1],o3[,1],o4[,1],o5[,1])
mat3=cbind(o1[,2],o2[,2],o3[,2],o4[,2],o5[,2])
o1=read.delim("../sub3/prediction/ova/individual_transplant/cor_nrmse_o1.txt",header=F,row.names=1)
o2=read.delim("../sub3/prediction/ova/individual_transplant/cor_nrmse_o2.txt",header=F,row.names=1)
o3=read.delim("../sub3/prediction/ova/individual_transplant/cor_nrmse_o3.txt",header=F,row.names=1)
o4=read.delim("../sub3/prediction/ova/individual_transplant/cor_nrmse_o4.txt",header=F,row.names=1)
o5=read.delim("../sub3/prediction/ova/individual_transplant/cor_nrmse_o5.txt",header=F,row.names=1)
mat2=cbind(o1[,1],o2[,1],o3[,1],o4[,1],o5[,1])
mat4=cbind(o1[,2],o2[,2],o3[,2],o4[,2],o5[,2])
cor2=apply(mat1-mat2,1,mean,na.rm=T)
rmse2=apply(mat3-mat4,1,mean,na.rm=T)
mean(cor2,na.rm=T) #[1] 0.04051451
mean(rmse2,na.rm=T) #[1] -0.004932747

# number of missing
b=as.matrix(read.delim("../sub3/data/trimmed_set/breast_phospho.txt",check.names=F,row.names=1))
o=as.matrix(read.delim("../sub3/data/trimmed_set/ova_phospho.txt",check.names=F,row.names=1))
num_miss1=apply(b,1,function(x){sum(is.na(x))})
num_miss2=apply(o,1,function(x){sum(is.na(x))})

# weak correlation
cor(cor1,num_miss1,use="pairwise.complete.obs") #[1] 0.2228953 p-value < 2.2e-16
cor(cor2,num_miss2,use="pairwise.complete.obs") #[1] 0.164645 p-value < 2.2e-16

avg_cor=c(mean(cor1[num_miss1>0 & num_miss1<=10],na.rm=T),
    mean(cor1[num_miss1>10 & num_miss1<=20],na.rm=T),
    mean(cor1[num_miss1>20 & num_miss1<=30],na.rm=T),
    mean(cor1[num_miss1>30 & num_miss1<=40],na.rm=T),
    mean(cor1[num_miss1>40 & num_miss1<=50],na.rm=T),
    mean(cor1[num_miss1>50 & num_miss1<=60],na.rm=T),
    mean(cor1[num_miss1>60 & num_miss1<=70],na.rm=T),
    mean(cor1[num_miss1>70 & num_miss1<=80],na.rm=T),
    mean(cor1[num_miss1>80 & num_miss1<=90],na.rm=T),
    mean(cor1[num_miss1>90 & num_miss1<=100],na.rm=T))

dat=data.frame(x=seq(10,100,10),y=avg_cor)

p1=ggplot(data=dat, aes(x, y)) +
  geom_bar(stat="identity",color=p_jco[6],fill=p_jco[6]) +
  theme_light() + 
  ylim(0,0.1) +
  geom_segment(aes(x = 5, y = 0.0138, xend = 105, yend = 0.0138), colour = p_jco[2], linetype=2, size=1) +
  theme(plot.title = element_text(hjust = 0.5)) + # hjust=0.5 set the title in the center
  labs(x='number of missing values',y="delta correlation",title="breast") + # customize titles
  scale_x_continuous(breaks=seq(10,100,10),labels=c("(0,10]","(10,20]","(20,30]","(30,40]","(40,50]",
    "(50,60]","(60,70]","(70,80]","(80,90]","(90,100]")) +
  annotate("text", x = seq(10,100,10), y = 0.1, label=format(dat[,"y"],digits=2),size=3)
p1

avg_cor=c(mean(cor2[num_miss2>0 & num_miss2<=10],na.rm=T),
    mean(cor2[num_miss2>10 & num_miss2<=20],na.rm=T),
    mean(cor2[num_miss2>20 & num_miss2<=30],na.rm=T),
    mean(cor2[num_miss2>30 & num_miss2<=40],na.rm=T),
    mean(cor2[num_miss2>40 & num_miss2<=50],na.rm=T),
    mean(cor2[num_miss2>50 & num_miss2<=60],na.rm=T),
    mean(cor2[num_miss2>60 & num_miss2<=70],na.rm=T))

dat=data.frame(x=seq(10,70,10),y=avg_cor)

p2=ggplot(data=dat, aes(x, y)) +
  geom_bar(stat="identity",color=p_jco[6],fill=p_jco[6]) +
  theme_light() + 
  ylim(0,0.15) +
  geom_segment(aes(x = 5, y = 0.0405, xend = 75, yend = 0.0405), colour = p_jco[2], linetype=2, size=1) +
  theme(plot.title = element_text(hjust = 0.5)) + # hjust=0.5 set the title in the center
  labs(x='number of missing values',y="delta correlation",title="ovary") + # customize titles
  scale_x_continuous(breaks=seq(10,70,10),labels=c("(0,10]","(10,20]","(20,30]","(30,40]","(40,50]","(50,60]","(60,70]")) + # x-axis label
  annotate("text", x = seq(10,70,10), y = 0.15, label=format(dat[,"y"],digits=2),size=3)
p2

# weak correlation
cor.test(rmse1,num_miss1,use="pairwise.complete.obs") #[1] -0.2672205 p-value < 2.2e-16
cor.test(rmse2,num_miss2,use="pairwise.complete.obs") #[1] -0.202962 p-value < 2.2e-16

avg_rmse=c(mean(rmse1[num_miss1>0 & num_miss1<=10],na.rm=T),
    mean(rmse1[num_miss1>10 & num_miss1<=20],na.rm=T),
    mean(rmse1[num_miss1>20 & num_miss1<=30],na.rm=T),
    mean(rmse1[num_miss1>30 & num_miss1<=40],na.rm=T),
    mean(rmse1[num_miss1>40 & num_miss1<=50],na.rm=T),
    mean(rmse1[num_miss1>50 & num_miss1<=60],na.rm=T),
    mean(rmse1[num_miss1>60 & num_miss1<=70],na.rm=T),
    mean(rmse1[num_miss1>70 & num_miss1<=80],na.rm=T),
    mean(rmse1[num_miss1>80 & num_miss1<=90],na.rm=T),
    mean(rmse1[num_miss1>90 & num_miss1<=100],na.rm=T))

dat=data.frame(x=seq(10,100,10),y=avg_rmse)

p3=ggplot(data=dat, aes(x, y)) +
  geom_bar(stat="identity",color=p_jco[6],fill=p_jco[6]) +
  theme_light() + 
  ylim(0.003,-0.03) +
  geom_segment(aes(x = 5, y = -0.00144, xend = 105, yend =-0.00144), colour = p_jco[2], linetype=2, size=1) +
  theme(plot.title = element_text(hjust = 0.5)) + # hjust=0.5 set the title in the center
  labs(x='number of missing values',y="delta NRMSE",title="breast") + # customize titles
  scale_x_continuous(breaks=seq(10,100,10),labels=c("(0,10]","(10,20]","(20,30]","(30,40]","(40,50]",
    "(50,60]","(60,70]","(70,80]","(80,90]","(90,100]")) +
  annotate("text", x = seq(10,100,10), y = -0.029, label=format(dat[,"y"],digits=2),size=3)
p3

avg_rmse=c(mean(rmse2[num_miss2>0 & num_miss2<=10],na.rm=T),
    mean(rmse2[num_miss2>10 & num_miss2<=20],na.rm=T),
    mean(rmse2[num_miss2>20 & num_miss2<=30],na.rm=T),
    mean(rmse2[num_miss2>30 & num_miss2<=40],na.rm=T),
    mean(rmse2[num_miss2>40 & num_miss2<=50],na.rm=T),
    mean(rmse2[num_miss2>50 & num_miss2<=60],na.rm=T),
    mean(rmse2[num_miss2>60 & num_miss2<=70],na.rm=T))

dat=data.frame(x=seq(10,70,10),y=avg_rmse)

p4=ggplot(data=dat, aes(x, y)) +
  geom_bar(stat="identity",color=p_jco[6],fill=p_jco[6]) +
  theme_light() + 
  ylim(0,-0.03) +
  geom_segment(aes(x = 5, y = -0.00493, xend = 75, yend = -0.00493), colour = p_jco[2], linetype=2, size=1) +
  theme(plot.title = element_text(hjust = 0.5)) + # hjust=0.5 set the title in the center
  labs(x='number of missing values',y="delta NRMSE",title="ovary") + # customize titles
  scale_x_continuous(breaks=seq(10,70,10),labels=c("(0,10]","(10,20]","(20,30]","(30,40]","(40,50]","(50,60]","(60,70]")) + # x-axis label
  annotate("text", x = seq(10,70,10), y = -0.028, label=format(dat[,"y"],digits=2),size=3)
p4

pdf(file="figure/result_delta_missing.pdf",width=12,height=6,useDingbats=F)
list_p=c(list(p1),list(p2),list(p3),list(p4))
mat_layout=matrix(1:4,nrow=2,byrow=T)
multiplot(plotlist=list_p,layout = mat_layout)
dev.off()

