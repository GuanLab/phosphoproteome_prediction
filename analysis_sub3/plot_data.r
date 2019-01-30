
library(ggplot2)
library(reshape2)
source("~/function/my_palette.r")
source("~/function/multiplot.R")
plot(1:length(p_jama),col=p_jama,pch=16,cex=5)
plot(1:length(p_jco),col=p_jco,pch=16,cex=5)
plot(1:length(p_d3),col=p_d3,pch=16,cex=5)

path0="/state3/hyangl/CPTAC/CPTAC_exp/sub3/data/trimmed_set/"

bp=as.matrix(read.delim(paste0(path0,"breast_phospho.txt"),header=T,row.names=1,check.names=F))
op=as.matrix(read.delim(paste0(path0,"ova_phospho.txt"),header=T,row.names=1,check.names=F))

# unique proteins
length(unique(sub("\\..*","",name1))) #[1] 4763
length(unique(sub("\\..*","",name2)))#[1] 2865

# phosphosite per protein
name1=rownames(bp)
table(table(sub("\\..*","",name1))) # AHNAK 302 sites
name2=rownames(op)
table(table(sub("\\..*","",name2))) # AHNAK 126 sites

# num_phospho
num_p1=table(sub("\\..*","",name1))[sub("\\..*","",name1)]
num_p2=table(sub("\\..*","",name2))[sub("\\..*","",name2)]
mean(table(sub("\\..*","",name2))) #[1] 3.510297
mean(table(sub("\\..*","",name1))) #[1] 6.714466

# num_aa
num_aa1=unlist(lapply(as.list(name1), function(x){length(unlist(strsplit(sub(".*:","",x),split="[s-y]")))-1}))
num_aa2=unlist(lapply(as.list(name2), function(x){length(unlist(strsplit(sub(".*:","",x),split="[s-y]")))-1}))

# aa type
s1=unlist(lapply(as.list(name1), function(x){sum(unlist(gregexpr('s',sub(".*:","",x)))!=(-1))}))
t1=unlist(lapply(as.list(name1), function(x){sum(unlist(gregexpr('t',sub(".*:","",x)))!=(-1))}))
y1=unlist(lapply(as.list(name1), function(x){sum(unlist(gregexpr('y',sub(".*:","",x)))!=(-1))}))

s2=unlist(lapply(as.list(name2), function(x){sum(unlist(gregexpr('s',sub(".*:","",x)))!=(-1))}))
t2=unlist(lapply(as.list(name2), function(x){sum(unlist(gregexpr('t',sub(".*:","",x)))!=(-1))}))
y2=unlist(lapply(as.list(name2), function(x){sum(unlist(gregexpr('y',sub(".*:","",x)))!=(-1))}))


# plot - number of unique phosphorylated peptides per protein

dat=data.frame(table(sub("\\..*","",name1)))
colnames(dat)=c("protein","num_phospho")
dat[dat[,2]>30,2]=31
p1=ggplot(dat,aes(x=num_phospho))+
    geom_histogram(aes(y=..density..),binwidth=0.5,boundary=0, colour=p_jco[4],fill='white') +
    theme_light() +
    xlim(0,31) +
    geom_density(alpha=0.5,colour=paste0(p_jco[4],30),fill=p_jco[4]) +
    labs(title="breast",x="number of unique phosphorylated peptides per protein") +
    theme(plot.title = element_text(hjust = 0.5))

dat=data.frame(table(sub("\\..*","",name2)))
colnames(dat)=c("protein","num_phospho")
dat[dat[,2]>30,2]=31
p2=ggplot(dat,aes(x=num_phospho))+
    geom_histogram(aes(y=..density..),binwidth=0.5,boundary=0, colour=p_jco[4],fill='white') +
    theme_light() +
    xlim(0,31) +
    geom_density(alpha=0.5,colour=paste0(p_jco[4],30),fill=p_jco[4])+
    labs(title="ovary",x="number of unique phosphorylated peptides per protein") +
    theme(plot.title = element_text(hjust = 0.5))

pdf(file="figure/data_num_phos_per_protein.pdf",width=10,height=6,useDingbats=F)
list_p=c(list(p1),list(p2))
mat_layout=matrix(1:2,nrow=2)
multiplot(plotlist=list_p,layout = mat_layout)
dev.off()


# plot - phosphorylation type

pdf(file="figure/data_phosphorylation.pdf",width=10,height=6,useDingbats=F)

layout(matrix(1:4,nrow=2,byrow=T))
pie(table(num_aa1),labels=c("mono-phosphorylated peptide (28001; 87.6%)",
    "dual-phosphorylated peptide (3624; 11.3%)","tri-phosphorylated peptide (356; 1.1%)"), 
    col=p_jco[c(6,2,9)], clockwise=T, init.angle=0, main="number of phosphorylated sites (breast)")

pie(table(num_aa2),labels=c("mono-phosphorylated peptide (9930; 98.7%)",
    "dual-phosphorylated peptide (127; 1.3%)"), 
    col=p_jco[c(6,2)], clockwise=T, init.angle=0, main="number of phosphorylated sites (ovary)")

pie(c(sum(s1),sum(t1),sum(y1)),labels=c("Serine (29868; 82.2%)","Threonine (5633; 15.5%)","Tyrosine (816; 2.2%)"), 
    col=p_jco[c(1,7,4)], clockwise=T, init.angle=0, main="phosphorylated amino acid types (breast)")

pie(c(sum(s2),sum(t2),sum(y2)),labels=c("Serine (8610; 84.5%)","Threonine (1410; 13.8%)","Tyrosine (164; 1.6%)"), 
    col=p_jco[c(1,7,4)], clockwise=T, init.angle=0, main="phosphorylated amino acid types (ovary)")

dev.off()


# 2aa; correlation vs distance
ind1=which(num_p1==2)
ind2=which(num_p2==2)

cor1=NULL
dist1=NULL
for (i in seq(1,length(ind1),2)){
    cor1=c(cor1,cor(bp[ind1[i],],bp[ind1[i+1],],use="pairwise.complete.obs"))
    tmp1=as.numeric(unlist(strsplit(sub(".*:","",name1[ind1[i]]),split="[s-y]"))[2])
    tmp2=as.numeric(unlist(strsplit(sub(".*:","",name1[ind1[i+1]]),split="[s-y]"))[2])
    dist1=c(dist1,tmp1-tmp2)
}

cor2=NULL
dist2=NULL
for (i in seq(1,length(ind2),2)){
    cor2=c(cor2,cor(bp[ind2[i],],bp[ind2[i+1],],use="pairwise.complete.obs"))
    tmp1=as.numeric(unlist(strsplit(sub(".*:","",name2[ind2[i]]),split="[s-y]"))[2])
    tmp2=as.numeric(unlist(strsplit(sub(".*:","",name2[ind2[i+1]]),split="[s-y]"))[2])
    dist2=c(dist2,tmp1-tmp2)
}

# weak correlation when dist<=10
cor(cor1,abs(dist1),use="pairwise.complete.obs") #[1] -0.0694795
tmp=abs(dist1)<=10                                        
cor(cor1[tmp],abs(dist1)[tmp],use="pairwise.complete.obs") #[1] -0.3042265

cor(cor2,abs(dist2),use="pairwise.complete.obs") #[1] 0.07753555
tmp=abs(dist2)<=20                                        
cor(cor2[tmp],abs(dist2)[tmp],use="pairwise.complete.obs") #[1] -0.1842406

# overall one phospho - multiple phospho; technical related to parent protein
mean(cor1,na.rm=T) #[1] 0.6819441
mean(cor2,na.rm=T) #[1] 0.6540425

# this doesn't work very well in ovary
mean(cor1[abs(dist1)<=5],na.rm=T) #[1] 0.8555879
mean(cor1[abs(dist1)<=10],na.rm=T) #[1] 0.8136129
mean(cor1[abs(dist1)<=15],na.rm=T) #[1] 0.8038601
mean(cor1[abs(dist1)<=20],na.rm=T) #[1] 0.7904415
mean(cor1[abs(dist1)<=25],na.rm=T) #[1] 0.7746622
mean(cor1[abs(dist1)<=30],na.rm=T) #[1] 0.7708199
mean(cor1[abs(dist1)<=35],na.rm=T) #[1] 0.7684284
mean(cor1[abs(dist1)<=40],na.rm=T) #[1] 0.7612395
mean(cor1[abs(dist1)<=45],na.rm=T) #[1] 0.7548709
mean(cor1[abs(dist1)<=50],na.rm=T) #[1] 0.7504941
mean(cor1[abs(dist1)>50],na.rm=T) #[1] 0.6132161


dat=data.frame(y=cor1,x=abs(dist1))

p3=ggplot(dat, aes(x,y,color="blue")) + # use column "name" to color 
    geom_point(size=3,pch=19,alpha=0.6) + # pch controls plot symbols; alpha controls the transparency
    #scale_colour_manual(values=the_palette) + # use the customized palette
    geom_smooth(method="lm", size=0.5, colour='blue', se=F) + # draw a regression line
    theme_light() + # light background
    ylim(-1,1) + # set y limit
    xlim(0.0,1000) + # set x limit
    theme(axis.line.x = element_line(color="black", size=0.5)) + # customize the font of y-axis (you can use the default)
    theme(axis.line.y = element_line(color="black", size=0.5)) +
    theme(axis.title.x = element_text(colour="black", size=12)) +
    theme(axis.title.y = element_text(colour="black", size=12)) +
    theme(axis.text.x = element_text(colour="black",size=12, angle=45, hjust=1)) + # rotate x-axis label by 45 degree
    theme(axis.text.y = element_text(colour="black",size=12)) +
    theme(plot.title = element_text(hjust = 0.5)) + # hjust=0.5 set the title in the center
    theme(legend.position="none") +
    labs(x='distance between phosphorylation sites',y="Pearson's correlation") + # customize titles
    annotate("text", x = 750, y = -0.2, label=paste0("correlation = -0.0695\np-value=0.0538"),size=5) 
p3


tmp=abs(dist1)<=10    
dat=data.frame(y=cor1[tmp],x=abs(dist1[tmp]))

p4=ggplot(dat, aes(x,y,color="blue")) + # use column "name" to color 
    geom_point(size=3,pch=19,alpha=0.6) + # pch controls plot symbols; alpha controls the transparency
    #scale_colour_manual(values=the_palette) + # use the customized palette
    geom_smooth(method="lm", size=0.5, colour='blue', se=F) + # draw a regression line
    theme_light() + # light background
    ylim(-1,1) + # set y limit
    xlim(0.0,10) + # set x limit
    theme(axis.line.x = element_line(color="black", size=0.5)) + # customize the font of y-axis (you can use the default)
    theme(axis.line.y = element_line(color="black", size=0.5)) +
    theme(axis.title.x = element_text(colour="black", size=12)) +
    theme(axis.title.y = element_text(colour="black", size=12)) +
    theme(axis.text.x = element_text(colour="black",size=12, angle=45, hjust=1)) + # rotate x-axis label by 45 degree
    theme(axis.text.y = element_text(colour="black",size=12)) +
    theme(plot.title = element_text(hjust = 0.5)) + # hjust=0.5 set the title in the center
    theme(legend.position="none") +
    labs(x='distance between phosphorylation sites',y="Pearson's correlation") + # customize titles
    annotate("text", x = 7.5, y = -0.2, label=paste0("correlation = -0.304\np-value=4.304e-6"),size=5) 
p4

avg_cor=c(mean(cor1[abs(dist1)<=5],na.rm=T),
    mean(cor1[abs(dist1)<=10],na.rm=T),
    mean(cor1[abs(dist1)<=15],na.rm=T),
    mean(cor1[abs(dist1)<=20],na.rm=T),
    mean(cor1[abs(dist1)<=25],na.rm=T),
    mean(cor1[abs(dist1)<=30],na.rm=T),
    mean(cor1[abs(dist1)<=35],na.rm=T),
    mean(cor1[abs(dist1)<=40],na.rm=T),
    mean(cor1[abs(dist1)<=45],na.rm=T),
    mean(cor1[abs(dist1)<=50],na.rm=T),
    mean(cor1[abs(dist1)>50],na.rm=T))
dat=data.frame(x=seq(5,55,5),y=avg_cor)

p5=ggplot(data=dat, aes(x, y)) +
  geom_bar(stat="identity",color=p_jco[6],fill=p_jco[6]) +
  theme_light() + 
  ylim(0,1) +
  geom_segment(aes(x = 2.5, y = 0.682, xend = 57.5, yend = 0.682), colour = p_jco[2], linetype=2, size=1) +
  theme(plot.title = element_text(hjust = 0.5)) + # hjust=0.5 set the title in the center
  labs(x='distance between phosphorylation sites',y="Pearson's correlation") + # customize titles
  scale_x_continuous(breaks=seq(5,55,5),label=paste0("<=",seq(5,55,5))) + # x-axis label
  annotate("text", x = seq(5,55,5), y = 0.9, label=format(dat[,"y"],digits=3),size=5)
p5

pdf(file="figure/data_correlation_vs_distance.pdf",width=10,height=6,useDingbats=F)
list_p=c(list(p3),list(p4),list(p5))
mat_layout=matrix(c(1,2,3,3),nrow=2,byrow=T)
multiplot(plotlist=list_p,layout = mat_layout)
dev.off()




