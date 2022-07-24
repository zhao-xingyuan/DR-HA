rm(list=ls())
this.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this.dir)
FF=read.csv("FFcount.csv",header=TRUE) # 78 cells, 103 obs, has hetero
FF$FFA=ifelse(FF$FF=="A",FF$FF.1,0)
FF$FFN=ifelse(FF$FF=="N",FF$FF.1,0)
FF$FFP=ifelse(FF$FF=="P",FF$FF.1,0)
FF=FF[,-c(6,7)]

index = rep('a',length(FF$IR))
for(i in 1:length(FF$IR))
{index[i] = paste(FF[i,1:5],collapse="")}

table(index)
data = aggregate(x = FF[,6:8], by = list(index,index,index), FUN = "sum")
sanitize=data
len = length(data$Group.1)
new = FF[1:len,1:5]
sp = strsplit(data$Group.1,"")
for(i in 1:len){
  new[i,] = sp[[i]]
}
data_new = data.frame(new,data[,4:6])

orig=data_new
sanitize=data_new

# 1-way Laplace
orig=data_new
eps=10^seq(-3,2,1)
library(OOmisc)
simu=100
way1=array(NA, dim = c(length(eps), 6, simu))
for (k in 1:simu){
  bydim1=matrix(NA,length(eps), 6)
for (t in 1:length(eps)){
sanitize=data_new
sanitize[,6:8]=sanitize[,6:8]+rlaplace(3*78,1/eps[t])
cell=sanitize[,6:8]
cell[cell<0]=0
sanitize[,6:8]=round(cell)
for (j in 1:5){
pstar=c(sum(sanitize[which(sanitize[,j]=="A"),6:8]),sum(sanitize[which(sanitize[,j]=="N"),6:8]),sum(sanitize[which(sanitize[,j]=="P"),6:8]))
por=c(sum(orig[which(orig[,j]=="A"),6:8]),sum(orig[which(orig[,j]=="N"),6:8]),sum(orig[which(orig[,j]=="P"),6:8]))
bydim1[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
}
bydim1[t,6]=(1/2)*sum(abs(apply(sanitize[,6:8],2,sum)/sum(apply(sanitize[,6:8],2,sum))-apply(orig[,6:8],2,sum)/sum(apply(orig[,6:8],2,sum))))
}
  way1[,,k]=bydim1
  }
TVD=apply(way1,c(1,2),mean)

boxdata=matrix(NA,36,3)
boxdata[,1]=c(rep(eps[1],6),rep(eps[2],6),rep(eps[3],6),rep(eps[4],6),rep(eps[5],6),rep(eps[6],6))
boxdata[,2]=c(TVD[1,],TVD[2,],TVD[3,],TVD[4,],TVD[5,],TVD[6,])
boxdata[,3]=rep(0,36)
boxdata=as.data.frame(boxdata)
colnames(boxdata)<-c("eps","tvd","delta")
boxdata$eps=as.factor(boxdata$eps)
boxdata$delta=as.factor(boxdata$delta)

library(VGAM)
library(scales)
library(ggplot2)
x11()
ggplot(boxdata, aes(x=eps, y=tvd),fill=delta) +
  geom_boxplot(width=0.35)+theme_minimal()+
  scale_y_continuous(limits = c(0, 0.15))+
  labs(x=expression(epsilon),y="total variation distance")+
  theme(axis.line = element_line(colour = "black"))+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

# 1-way aDP
simu=100
eps=10^seq(-3,0,1)
way1=array(NA, dim = c(length(eps), 6, simu))
delta=0.001
for (k in 1:simu){
  bydim1=matrix(NA,length(eps), 6)
  for (t in 1:length(eps)){
    sanitize=data_new
    sigma=sqrt(2*log(1.25/delta))/eps[t]
    sanitize[,6:8]=sanitize[,6:8]+rnorm(3*78,0,sigma)
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:5){
      pstar=c(sum(sanitize[which(sanitize[,j]=="A"),6:8]),sum(sanitize[which(sanitize[,j]=="N"),6:8]),sum(sanitize[which(sanitize[,j]=="P"),6:8]))
      por=c(sum(orig[which(orig[,j]=="A"),6:8]),sum(orig[which(orig[,j]=="N"),6:8]),sum(orig[which(orig[,j]=="P"),6:8]))
      bydim1[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
    bydim1[t,6]=(1/2)*sum(abs(apply(sanitize[,6:8],2,sum)/sum(apply(sanitize[,6:8],2,sum))-apply(orig[,6:8],2,sum)/sum(apply(orig[,6:8],2,sum))))
  }
  way1[,,k]=bydim1
  }
TVD.0.001=apply(way1,c(1,2),mean)

way1=array(NA, dim = c(length(eps), 6, simu))
delta=0.01
for (k in 1:simu){
  bydim1=matrix(NA,length(eps), 6)
  for (t in 1:length(eps)){
    sanitize=data_new
    sigma=sqrt(2*log(1.25/delta))/eps[t]
    sanitize[,6:8]=sanitize[,6:8]+rnorm(3*78,0,sigma)
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:5){
      pstar=c(sum(sanitize[which(sanitize[,j]=="A"),6:8]),sum(sanitize[which(sanitize[,j]=="N"),6:8]),sum(sanitize[which(sanitize[,j]=="P"),6:8]))
      por=c(sum(orig[which(orig[,j]=="A"),6:8]),sum(orig[which(orig[,j]=="N"),6:8]),sum(orig[which(orig[,j]=="P"),6:8]))
      bydim1[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
    bydim1[t,6]=(1/2)*sum(abs(apply(sanitize[,6:8],2,sum)/sum(apply(sanitize[,6:8],2,sum))-apply(orig[,6:8],2,sum)/sum(apply(orig[,6:8],2,sum))))
  }
  way1[,,k]=bydim1
}
TVD.0.01=apply(way1,c(1,2),mean)

way1=array(NA, dim = c(length(eps), 6, simu))
delta=0.1
for (k in 1:simu){
  bydim1=matrix(NA,length(eps), 6)
  for (t in 1:length(eps)){
    sanitize=data_new
    sigma=sqrt(2*log(1.25/delta))/eps[t]
    sanitize[,6:8]=sanitize[,6:8]+rnorm(3*78,0,sigma)
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:5){
      pstar=c(sum(sanitize[which(sanitize[,j]=="A"),6:8]),sum(sanitize[which(sanitize[,j]=="N"),6:8]),sum(sanitize[which(sanitize[,j]=="P"),6:8]))
      por=c(sum(orig[which(orig[,j]=="A"),6:8]),sum(orig[which(orig[,j]=="N"),6:8]),sum(orig[which(orig[,j]=="P"),6:8]))
      bydim1[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
    bydim1[t,6]=(1/2)*sum(abs(apply(sanitize[,6:8],2,sum)/sum(apply(sanitize[,6:8],2,sum))-apply(orig[,6:8],2,sum)/sum(apply(orig[,6:8],2,sum))))
  }
  way1[,,k]=bydim1
}
TVD.0.1=apply(way1,c(1,2),mean)

boxdata=matrix(NA,72,3)
boxdata[,1]=c(rep(0.001,24),rep(0.01,24),rep(0.1,24))
boxdata[,2]=rep(c(rep(eps[1],6),rep(eps[2],6),rep(eps[3],6),rep(eps[4],6)),3)
boxdata[,3]=c(TVD.0.001[1,],TVD.0.001[2,],TVD.0.001[3,],TVD.0.001[4,],TVD.0.01[1,],TVD.0.01[2,],TVD.0.01[3,],TVD.0.01[4,],TVD.0.1[1,],TVD.0.1[2,],TVD.0.1[3,],TVD.0.1[4,])
boxdata=as.data.frame(boxdata)
colnames(boxdata)<-c("delta","eps","tvd")
boxdata$eps=as.factor(boxdata$eps)
boxdata$delta=as.factor(boxdata$delta)


library(ggplot2)
x11()
ggplot(boxdata, aes(x=eps, y=tvd, fill=delta)) +
  geom_boxplot(position=position_dodge(1),width=0.45)+
  scale_y_continuous(limits = c(0, 0.15))+
  theme_minimal()+
  labs(x=expression(epsilon),y="total variation distance")+
  theme(axis.line = element_line(colour = "black"))+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  scale_fill_discrete(name = "Dose", labels=c(expression(paste(delta,"=0.001")),expression(paste(delta,"=0.01")),expression(paste(delta,"=0.1"))))+
  theme(legend.text.align = 0)+
  theme(legend.text = element_text(size = 15))


# 1-way pDP
simu=100
eps=10^seq(-3,2,1)
way1=array(NA, dim = c(length(eps), 6, simu))
delta=0.001
for (k in 1:simu){
  bydim1=matrix(NA,length(eps), 6)
  for (t in 1:length(eps)){
    sanitize=data_new
    sigma=(sqrt(qnorm(delta/2)^2+2*eps[t])-qnorm(delta/2))/(2*eps[t])
    sanitize[,6:8]=sanitize[,6:8]+rnorm(3*78,0,sigma)
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:5){
      pstar=c(sum(sanitize[which(sanitize[,j]=="A"),6:8]),sum(sanitize[which(sanitize[,j]=="N"),6:8]),sum(sanitize[which(sanitize[,j]=="P"),6:8]))
      por=c(sum(orig[which(orig[,j]=="A"),6:8]),sum(orig[which(orig[,j]=="N"),6:8]),sum(orig[which(orig[,j]=="P"),6:8]))
      bydim1[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
    bydim1[t,6]=(1/2)*sum(abs(apply(sanitize[,6:8],2,sum)/sum(apply(sanitize[,6:8],2,sum))-apply(orig[,6:8],2,sum)/sum(apply(orig[,6:8],2,sum))))
  }
  way1[,,k]=bydim1
}
TVD.0.001=apply(way1,c(1,2),mean)

way1=array(NA, dim = c(length(eps), 6, simu))
delta=0.01
for (k in 1:simu){
  bydim1=matrix(NA,length(eps), 6)
  for (t in 1:length(eps)){
    sanitize=data_new
    sigma=(sqrt(qnorm(delta/2)^2+2*eps[t])-qnorm(delta/2))/(2*eps[t])
    sanitize[,6:8]=sanitize[,6:8]+rnorm(3*78,0,sigma)
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:5){
      pstar=c(sum(sanitize[which(sanitize[,j]=="A"),6:8]),sum(sanitize[which(sanitize[,j]=="N"),6:8]),sum(sanitize[which(sanitize[,j]=="P"),6:8]))
      por=c(sum(orig[which(orig[,j]=="A"),6:8]),sum(orig[which(orig[,j]=="N"),6:8]),sum(orig[which(orig[,j]=="P"),6:8]))
      bydim1[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
    bydim1[t,6]=(1/2)*sum(abs(apply(sanitize[,6:8],2,sum)/sum(apply(sanitize[,6:8],2,sum))-apply(orig[,6:8],2,sum)/sum(apply(orig[,6:8],2,sum))))
  }
  way1[,,k]=bydim1
}
TVD.0.01=apply(way1,c(1,2),mean)

way1=array(NA, dim = c(length(eps), 6, simu))
delta=0.1
for (k in 1:simu){
  bydim1=matrix(NA,length(eps), 6)
  for (t in 1:length(eps)){
    sanitize=data_new
    sigma=(sqrt(qnorm(delta/2)^2+2*eps[t])-qnorm(delta/2))/(2*eps[t])
    sanitize[,6:8]=sanitize[,6:8]+rnorm(3*78,0,sigma)
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:5){
      pstar=c(sum(sanitize[which(sanitize[,j]=="A"),6:8]),sum(sanitize[which(sanitize[,j]=="N"),6:8]),sum(sanitize[which(sanitize[,j]=="P"),6:8]))
      por=c(sum(orig[which(orig[,j]=="A"),6:8]),sum(orig[which(orig[,j]=="N"),6:8]),sum(orig[which(orig[,j]=="P"),6:8]))
      bydim1[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
    bydim1[t,6]=(1/2)*sum(abs(apply(sanitize[,6:8],2,sum)/sum(apply(sanitize[,6:8],2,sum))-apply(orig[,6:8],2,sum)/sum(apply(orig[,6:8],2,sum))))
  }
  way1[,,k]=bydim1
}
TVD.0.1=apply(way1,c(1,2),mean)

boxdata=matrix(NA,108,3)
boxdata[,1]=c(rep(0.001,36),rep(0.01,36),rep(0.1,36))
boxdata[,2]=rep(c(rep(eps[1],6),rep(eps[2],6),rep(eps[3],6),rep(eps[4],6),rep(eps[5],6),rep(eps[6],6)),3)
boxdata[,3]=c(TVD.0.001[1,],TVD.0.001[2,],TVD.0.001[3,],TVD.0.001[4,],TVD.0.001[5,],TVD.0.001[6,],TVD.0.01[1,],TVD.0.01[2,],TVD.0.01[3,],TVD.0.01[4,],TVD.0.01[5,],TVD.0.01[6,],TVD.0.1[1,],TVD.0.1[2,],TVD.0.1[3,],TVD.0.1[4,],TVD.0.1[5,],TVD.0.1[6,])
boxdata=as.data.frame(boxdata)
colnames(boxdata)<-c("delta","eps","tvd")
boxdata$eps=as.factor(boxdata$eps)
boxdata$delta=as.factor(boxdata$delta)


library(ggplot2)
x11()
ggplot(boxdata, aes(x=eps, y=tvd, fill=delta)) +
  geom_boxplot(position=position_dodge(1),width=0.45)+
  scale_y_continuous(limits = c(0, 0.15))+
  theme_minimal()+
  labs(x=expression(epsilon),y="total variation distance")+
  theme(axis.line = element_line(colour = "black"))+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  scale_fill_discrete(name = "Dose", labels=c(expression(paste(delta,"=0.001")),expression(paste(delta,"=0.01")),expression(paste(delta,"=0.1"))))+
  theme(legend.text.align = 0)+
  theme(legend.text = element_text(size = 15))

# 2-way Laplace
orig=data_new
eps=10^seq(-3,2,1)
library(OOmisc)
simu=100
way2=array(NA, dim = c(length(eps), 15, simu))
c52=combn(seq(1,5,1),2)  #2*10
for (k in 1:simu){
  bydim2=matrix(NA,length(eps),15)
  for (t in 1:length(eps)){
    sanitize=data_new
    sanitize[,6:8]=sanitize[,6:8]+rlaplace(3*78,1/eps[t])
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:10){
      pstar=c(ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),6:8]),0))
      por=c(ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),6:8]),0))
      bydim2[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
      }
    for (j in 1:5){
      pstar=c(ifelse(length(which(sanitize[,j]=="A"))>0,sum(sanitize[which(sanitize[,j]=="A"),6]),0),
              ifelse(length(which(sanitize[,j]=="A"))>0,sum(sanitize[which(sanitize[,j]=="A"),7]),0),
              ifelse(length(which(sanitize[,j]=="A"))>0,sum(sanitize[which(sanitize[,j]=="A"),8]),0),
              ifelse(length(which(sanitize[,j]=="N"))>0,sum(sanitize[which(sanitize[,j]=="N"),6]),0),
              ifelse(length(which(sanitize[,j]=="N"))>0,sum(sanitize[which(sanitize[,j]=="N"),7]),0),
              ifelse(length(which(sanitize[,j]=="N"))>0,sum(sanitize[which(sanitize[,j]=="N"),8]),0),
              ifelse(length(which(sanitize[,j]=="P"))>0,sum(sanitize[which(sanitize[,j]=="P"),6]),0),
              ifelse(length(which(sanitize[,j]=="P"))>0,sum(sanitize[which(sanitize[,j]=="P"),7]),0),
              ifelse(length(which(sanitize[,j]=="P"))>0,sum(sanitize[which(sanitize[,j]=="P"),8]),0))
      por=c(ifelse(length(which(orig[,j]=="A"))>0,sum(orig[which(orig[,j]=="A"),6]),0),
            ifelse(length(which(orig[,j]=="A"))>0,sum(orig[which(orig[,j]=="A"),7]),0),
            ifelse(length(which(orig[,j]=="A"))>0,sum(orig[which(orig[,j]=="A"),8]),0),
            ifelse(length(which(orig[,j]=="N"))>0,sum(orig[which(orig[,j]=="N"),6]),0),
            ifelse(length(which(orig[,j]=="N"))>0,sum(orig[which(orig[,j]=="N"),7]),0),
            ifelse(length(which(orig[,j]=="N"))>0,sum(orig[which(orig[,j]=="N"),8]),0),
            ifelse(length(which(orig[,j]=="P"))>0,sum(orig[which(orig[,j]=="P"),6]),0),
            ifelse(length(which(orig[,j]=="P"))>0,sum(orig[which(orig[,j]=="P"),7]),0),
            ifelse(length(which(orig[,j]=="P"))>0,sum(orig[which(orig[,j]=="P"),8]),0))
      bydim2[t,j+10]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
  }
  way2[,,k]=bydim2
}

TVD=apply(way2,c(1,2),mean)

boxdata=matrix(NA,90,3)
boxdata[,1]=c(rep(eps[1],15),rep(eps[2],15),rep(eps[3],15),rep(eps[4],15),rep(eps[5],15),rep(eps[6],15))
boxdata[,2]=c(TVD[1,],TVD[2,],TVD[3,],TVD[4,],TVD[5,],TVD[6,])
boxdata[,3]=rep(0,90)
boxdata=as.data.frame(boxdata)
colnames(boxdata)<-c("eps","tvd","delta")
boxdata$eps=as.factor(boxdata$eps)
boxdata$delta=as.factor(boxdata$delta)

library(VGAM)
library(scales)
library(ggplot2)
x11()
ggplot(boxdata, aes(x=eps, y=tvd)) +
  geom_boxplot(width=0.45)+theme_minimal()+
  scale_y_continuous(limits=c(0,0.4))+
  labs(x=expression(epsilon),y="total variation distance")+
  theme(axis.line = element_line(colour = "black"))+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))

# 2-way aDP
simu=100
eps=10^seq(-3,0,1)
way2=array(NA, dim = c(length(eps), 15, simu))
delta=0.001
c52=combn(seq(1,5,1),2)
for (k in 1:simu){
  bydim2=matrix(NA,length(eps),15)
  for (t in 1:length(eps)){
    sanitize=data_new
    sigma=sqrt(2*log(1.25/delta))/eps[t]
    sanitize[,6:8]=sanitize[,6:8]+rnorm(3*78,0,sigma)
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:10){
      pstar=c(ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),6:8]),0))
      por=c(ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),6:8]),0))
      bydim2[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
    for (j in 1:5){
      pstar=c(sum(sanitize[which((sanitize[,j]=="A")),6]),sum(sanitize[which(sanitize[,j]=="A"),7]),sum(sanitize[which(sanitize[,j]=="A"),8]),
              sum(sanitize[which((sanitize[,j]=="N")),6]),sum(sanitize[which(sanitize[,j]=="N"),7]),sum(sanitize[which(sanitize[,j]=="N"),8]),
              sum(sanitize[which((sanitize[,j]=="P")),6]),sum(sanitize[which(sanitize[,j]=="P"),7]),sum(sanitize[which(sanitize[,j]=="P"),8]))
      por=c(sum(orig[which((orig[,j]=="A")),6]),sum(orig[which(orig[,j]=="A"),7]),sum(orig[which(orig[,j]=="A"),8]),
            sum(orig[which((orig[,j]=="N")),6]),sum(orig[which(orig[,j]=="N"),7]),sum(orig[which(orig[,j]=="N"),8]),
            sum(orig[which((orig[,j]=="P")),6]),sum(orig[which(orig[,j]=="P"),7]),sum(orig[which(orig[,j]=="P"),8]))
      bydim2[t,j+10]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
  }
  way2[,,k]=bydim2
  }
TVD.0.001=apply(way2,c(1,2),mean)

way2=array(NA, dim = c(length(eps), 15, simu))
delta=0.01
for (k in 1:simu){
  bydim2=matrix(NA,length(eps),15)
  for (t in 1:length(eps)){
    sanitize=data_new
    sigma=sqrt(2*log(1.25/delta))/eps[t]
    sanitize[,6:8]=sanitize[,6:8]+rnorm(3*78,0,sigma)
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:10){
      pstar=c(ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),6:8]),0))
      por=c(ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),6:8]),0))
      bydim2[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
    for (j in 1:5){
      pstar=c(sum(sanitize[which((sanitize[,j]=="A")),6]),sum(sanitize[which(sanitize[,j]=="A"),7]),sum(sanitize[which(sanitize[,j]=="A"),8]),
              sum(sanitize[which((sanitize[,j]=="N")),6]),sum(sanitize[which(sanitize[,j]=="N"),7]),sum(sanitize[which(sanitize[,j]=="N"),8]),
              sum(sanitize[which((sanitize[,j]=="P")),6]),sum(sanitize[which(sanitize[,j]=="P"),7]),sum(sanitize[which(sanitize[,j]=="P"),8]))
      por=c(sum(orig[which((orig[,j]=="A")),6]),sum(orig[which(orig[,j]=="A"),7]),sum(orig[which(orig[,j]=="A"),8]),
            sum(orig[which((orig[,j]=="N")),6]),sum(orig[which(orig[,j]=="N"),7]),sum(orig[which(orig[,j]=="N"),8]),
            sum(orig[which((orig[,j]=="P")),6]),sum(orig[which(orig[,j]=="P"),7]),sum(orig[which(orig[,j]=="P"),8]))
      bydim2[t,j+10]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
  }
  way2[,,k]=bydim2
}
TVD.0.01=apply(way2,c(1,2),mean)

way2=array(NA, dim = c(length(eps), 15, simu))
delta=0.1
for (k in 1:simu){
  bydim2=matrix(NA,length(eps),15)
  for (t in 1:length(eps)){
    sanitize=data_new
    sigma=sqrt(2*log(1.25/delta))/eps[t]
    sanitize[,6:8]=sanitize[,6:8]+rnorm(3*78,0,sigma)
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:10){
      pstar=c(ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),6:8]),0))
      por=c(ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),6:8]),0))
      bydim2[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
    for (j in 1:5){
      pstar=c(sum(sanitize[which((sanitize[,j]=="A")),6]),sum(sanitize[which(sanitize[,j]=="A"),7]),sum(sanitize[which(sanitize[,j]=="A"),8]),
              sum(sanitize[which((sanitize[,j]=="N")),6]),sum(sanitize[which(sanitize[,j]=="N"),7]),sum(sanitize[which(sanitize[,j]=="N"),8]),
              sum(sanitize[which((sanitize[,j]=="P")),6]),sum(sanitize[which(sanitize[,j]=="P"),7]),sum(sanitize[which(sanitize[,j]=="P"),8]))
      por=c(sum(orig[which((orig[,j]=="A")),6]),sum(orig[which(orig[,j]=="A"),7]),sum(orig[which(orig[,j]=="A"),8]),
            sum(orig[which((orig[,j]=="N")),6]),sum(orig[which(orig[,j]=="N"),7]),sum(orig[which(orig[,j]=="N"),8]),
            sum(orig[which((orig[,j]=="P")),6]),sum(orig[which(orig[,j]=="P"),7]),sum(orig[which(orig[,j]=="P"),8]))
      bydim2[t,j+10]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
  }
  way2[,,k]=bydim2
}
TVD.0.1=apply(way2,c(1,2),mean)


boxdata=matrix(NA,180,3)
boxdata[,1]=c(rep(0.001,60),rep(0.01,60),rep(0.1,60))
boxdata[,2]=rep(c(rep(eps[1],15),rep(eps[2],15),rep(eps[3],15),rep(eps[4],15)),3)
boxdata[,3]=c(TVD.0.001[1,],TVD.0.001[2,],TVD.0.001[3,],TVD.0.001[4,],TVD.0.01[1,],TVD.0.01[2,],TVD.0.01[3,],TVD.0.01[4,],TVD.0.1[1,],TVD.0.1[2,],TVD.0.1[3,],TVD.0.1[4,])
boxdata=as.data.frame(boxdata)
colnames(boxdata)<-c("delta","eps","tvd")
boxdata$eps=as.factor(boxdata$eps)
boxdata$delta=as.factor(boxdata$delta)

library(ggplot2)
x11()
ggplot(boxdata, aes(x=eps, y=tvd, fill=delta)) +
  geom_boxplot(position=position_dodge(1),width=0.45)+
  scale_y_continuous(limits=c(0,0.4))+
  theme_minimal()+
  labs(x=expression(epsilon),y="total variation distance")+
  theme(axis.line = element_line(colour = "black"))+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  scale_fill_discrete(name = "Dose", labels=c(expression(paste(delta,"=0.001")),expression(paste(delta,"=0.01")),expression(paste(delta,"=0.1"))))+
  theme(legend.text.align = 0)+
  theme(legend.text = element_text(size = 15))

# 2-way pDP
simu=100
eps=10^seq(-3,2,1)
way2=array(NA, dim = c(length(eps), 15, simu))
delta=0.001
c52=combn(seq(1,5,1),2)
for (k in 1:simu){
  bydim2=matrix(NA,length(eps),15)
  for (t in 1:length(eps)){
    sanitize=data_new
    sigma=(sqrt(qnorm(delta/2)^2+2*eps[t])-qnorm(delta/2))/(2*eps[t])
    sanitize[,6:8]=sanitize[,6:8]+rnorm(3*78,0,sigma)
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:10){
      pstar=c(ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),6:8]),0))
      por=c(ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),6:8]),0))
      bydim2[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
    for (j in 1:5){
      pstar=c(sum(sanitize[which((sanitize[,j]=="A")),6]),sum(sanitize[which(sanitize[,j]=="A"),7]),sum(sanitize[which(sanitize[,j]=="A"),8]),
              sum(sanitize[which((sanitize[,j]=="N")),6]),sum(sanitize[which(sanitize[,j]=="N"),7]),sum(sanitize[which(sanitize[,j]=="N"),8]),
              sum(sanitize[which((sanitize[,j]=="P")),6]),sum(sanitize[which(sanitize[,j]=="P"),7]),sum(sanitize[which(sanitize[,j]=="P"),8]))
      por=c(sum(orig[which((orig[,j]=="A")),6]),sum(orig[which(orig[,j]=="A"),7]),sum(orig[which(orig[,j]=="A"),8]),
            sum(orig[which((orig[,j]=="N")),6]),sum(orig[which(orig[,j]=="N"),7]),sum(orig[which(orig[,j]=="N"),8]),
            sum(orig[which((orig[,j]=="P")),6]),sum(orig[which(orig[,j]=="P"),7]),sum(orig[which(orig[,j]=="P"),8]))
      bydim2[t,j+10]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
  }
  way2[,,k]=bydim2
}
TVD.0.001=apply(way2,c(1,2),mean)

way2=array(NA, dim = c(length(eps), 15, simu))
delta=0.01
for (k in 1:simu){
  bydim2=matrix(NA,length(eps),15)
  for (t in 1:length(eps)){
    sanitize=data_new
    sigma=(sqrt(qnorm(delta/2)^2+2*eps[t])-qnorm(delta/2))/(2*eps[t])
    sanitize[,6:8]=sanitize[,6:8]+rnorm(3*78,0,sigma)
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:10){
      pstar=c(ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),6:8]),0))
      por=c(ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),6:8]),0))
      bydim2[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
    for (j in 1:5){
      pstar=c(sum(sanitize[which((sanitize[,j]=="A")),6]),sum(sanitize[which(sanitize[,j]=="A"),7]),sum(sanitize[which(sanitize[,j]=="A"),8]),
              sum(sanitize[which((sanitize[,j]=="N")),6]),sum(sanitize[which(sanitize[,j]=="N"),7]),sum(sanitize[which(sanitize[,j]=="N"),8]),
              sum(sanitize[which((sanitize[,j]=="P")),6]),sum(sanitize[which(sanitize[,j]=="P"),7]),sum(sanitize[which(sanitize[,j]=="P"),8]))
      por=c(sum(orig[which((orig[,j]=="A")),6]),sum(orig[which(orig[,j]=="A"),7]),sum(orig[which(orig[,j]=="A"),8]),
            sum(orig[which((orig[,j]=="N")),6]),sum(orig[which(orig[,j]=="N"),7]),sum(orig[which(orig[,j]=="N"),8]),
            sum(orig[which((orig[,j]=="P")),6]),sum(orig[which(orig[,j]=="P"),7]),sum(orig[which(orig[,j]=="P"),8]))
      bydim2[t,j+10]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
  }
  way2[,,k]=bydim2
}
TVD.0.01=apply(way2,c(1,2),mean)

way2=array(NA, dim = c(length(eps), 15, simu))
delta=0.1
for (k in 1:simu){
  bydim2=matrix(NA,length(eps),15)
  for (t in 1:length(eps)){
    sanitize=data_new
    sigma=(sqrt(qnorm(delta/2)^2+2*eps[t])-qnorm(delta/2))/(2*eps[t])
    sanitize[,6:8]=sanitize[,6:8]+rnorm(3*78,0,sigma)
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:10){
      pstar=c(ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")))>0,sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),6:8]),0))
      por=c(ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")))>0,sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),6:8]),0))
      bydim2[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
    for (j in 1:5){
      pstar=c(sum(sanitize[which((sanitize[,j]=="A")),6]),sum(sanitize[which(sanitize[,j]=="A"),7]),sum(sanitize[which(sanitize[,j]=="A"),8]),
              sum(sanitize[which((sanitize[,j]=="N")),6]),sum(sanitize[which(sanitize[,j]=="N"),7]),sum(sanitize[which(sanitize[,j]=="N"),8]),
              sum(sanitize[which((sanitize[,j]=="P")),6]),sum(sanitize[which(sanitize[,j]=="P"),7]),sum(sanitize[which(sanitize[,j]=="P"),8]))
      por=c(sum(orig[which((orig[,j]=="A")),6]),sum(orig[which(orig[,j]=="A"),7]),sum(orig[which(orig[,j]=="A"),8]),
            sum(orig[which((orig[,j]=="N")),6]),sum(orig[which(orig[,j]=="N"),7]),sum(orig[which(orig[,j]=="N"),8]),
            sum(orig[which((orig[,j]=="P")),6]),sum(orig[which(orig[,j]=="P"),7]),sum(orig[which(orig[,j]=="P"),8]))
      bydim2[t,j+10]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
  }
  way2[,,k]=bydim2
}
TVD.0.1=apply(way2,c(1,2),mean)


boxdata=matrix(NA,270,3)
boxdata[,1]=c(rep(0.001,90),rep(0.01,90),rep(0.1,90))
boxdata[,2]=rep(c(rep(eps[1],15),rep(eps[2],15),rep(eps[3],15),rep(eps[4],15),rep(eps[5],15),rep(eps[6],15)),3)
boxdata[,3]=c(TVD.0.001[1,],TVD.0.001[2,],TVD.0.001[3,],TVD.0.001[4,],TVD.0.001[5,],TVD.0.001[6,],TVD.0.01[1,],TVD.0.01[2,],TVD.0.01[3,],TVD.0.01[4,],TVD.0.01[5,],TVD.0.01[6,],TVD.0.1[1,],TVD.0.1[2,],TVD.0.1[3,],TVD.0.1[4,],TVD.0.1[5,],TVD.0.1[6,])
boxdata=as.data.frame(boxdata)
colnames(boxdata)<-c("delta","eps","tvd")
boxdata$eps=as.factor(boxdata$eps)
boxdata$delta=as.factor(boxdata$delta)

library(ggplot2)
library(scales)
x11()
ggplot(boxdata, aes(x=eps, y=tvd, fill=delta)) +
  geom_boxplot(position=position_dodge(1),width=0.45)+
  scale_y_continuous(limits=c(0,0.4))+
  theme_minimal()+
  labs(x=expression(epsilon),y="total variation distance")+
  theme(axis.line = element_line(colour = "black"))+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  scale_fill_discrete(name = "Dose", labels=c(expression(paste(delta,"=0.001")),expression(paste(delta,"=0.01")),expression(paste(delta,"=0.1"))))+
  theme(legend.text.align = 0)+
  theme(legend.text = element_text(size = 15))

# 3-way Laplace
orig=data_new
eps=10^seq(-3,2,1)
library(OOmisc)
simu=80
way3=array(NA, dim = c(length(eps), 20, simu))
c53=combn(seq(1,5,1),3)  #3*10
c52=combn(seq(1,5,1),2)  #2*10
for (k in 1:simu){
  bydim3=matrix(NA,length(eps), 20)
  for (t in 1:length(eps)){
    sanitize=data_new
    sanitize[,6:8]=sanitize[,6:8]+rlaplace(3*78,1/eps[t])
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:10){
        pstar=c(ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
                ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0))
        por=c(ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0))
        bydim3[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
      }
    for (j in 1:10){
      pstar=c(sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),8]))
      por=c(sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),8]))
      bydim3[t,10+j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
    }    
  way3[,,k]=bydim3
  }

TVD=apply(way3,c(1,2),mean)

boxdata=matrix(NA,120,3)
boxdata[,1]=c(rep(eps[1],20),rep(eps[2],20),rep(eps[3],20),rep(eps[4],20),rep(eps[5],20),rep(eps[6],20))
boxdata[,2]=c(TVD[1,],TVD[2,],TVD[3,],TVD[4,],TVD[5,],TVD[6,])
boxdata[,3]=rep(0,120)
boxdata=as.data.frame(boxdata)
colnames(boxdata)<-c("eps","tvd","delta")
boxdata$eps=as.factor(boxdata$eps)
boxdata$delta=as.factor(boxdata$delta)


library(VGAM)
library(scales)
library(ggplot2)
x11()
ggplot(boxdata, aes(x=eps, y=tvd),fill=delta) +
         geom_boxplot(width=0.35)+theme_minimal()+
  scale_y_continuous(limits=c(0,0.5))+
  labs(x=expression(epsilon),y="total variation distance")+
  theme(axis.line = element_line(colour = "black"))+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))


# 3-way aDP
simu=80 
eps=10^seq(-3,0,1)
way3=array(NA, dim = c(length(eps), 20, simu))
delta=0.001
c53=combn(seq(1,5,1),3)  #3*10
c52=combn(seq(1,5,1),2)  #2*10
for (k in 1:simu){
  bydim3=matrix(NA,length(eps), 20)
  for (t in 1:length(eps)){
    sanitize=data_new
    sigma=sqrt(2*log(1.25/delta))/eps[t]
    sanitize[,6:8]=sanitize[,6:8]+rnorm(3*78,0,sigma)
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:10){
      pstar=c(ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0))
      por=c(ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0))
      bydim3[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
    for (j in 1:10){
      pstar=c(sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),8]))
      por=c(sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),8]))
      bydim3[t,10+j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
  }    
  way3[,,k]=bydim3
}
TVD.0.001=apply(way3,c(1,2),mean)

way3=array(NA, dim = c(length(eps), 20, simu))
delta=0.01
for (k in 1:simu){
  bydim3=matrix(NA,length(eps), 20)
  for (t in 1:length(eps)){
    sanitize=data_new
    sigma=sqrt(2*log(1.25/delta))/eps[t]
    sanitize[,6:8]=sanitize[,6:8]+rnorm(3*78,0,sigma)
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:10){
      pstar=c(ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0))
      por=c(ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0))
      bydim3[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
    for (j in 1:10){
      pstar=c(sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),8]))
      por=c(sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),8]))
      bydim3[t,10+j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
  }    
  way3[,,k]=bydim3
}
TVD.0.01<-apply(way3,c(1,2),mean)

way3=array(NA, dim = c(length(eps), 20, simu))
delta=0.1
for (k in 1:simu){
  bydim3=matrix(NA,length(eps), 20)
  for (t in 1:length(eps)){
    sanitize=data_new
    sigma=sqrt(2*log(1.25/delta))/eps[t]
    sanitize[,6:8]=sanitize[,6:8]+rnorm(3*78,0,sigma)
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:10){
      pstar=c(ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0))
      por=c(ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0))
      bydim3[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
    for (j in 1:10){
      pstar=c(sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),8]))
      por=c(sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),8]))
      bydim3[t,10+j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
  }    
  way3[,,k]=bydim3
}
TVD.0.1<-apply(way3,c(1,2),mean)


boxdata=matrix(NA,240,3)
boxdata[,1]=c(rep(0.001,80),rep(0.01,80),rep(0.1,80))
boxdata[,2]=rep(c(rep(eps[1],20),rep(eps[2],20),rep(eps[3],20),rep(eps[4],20)),3)
boxdata[,3]=c(TVD.0.001[1,],TVD.0.001[2,],TVD.0.001[3,],TVD.0.001[4,],TVD.0.01[1,],TVD.0.01[2,],TVD.0.01[3,],TVD.0.01[4,],TVD.0.1[1,],TVD.0.1[2,],TVD.0.1[3,],TVD.0.1[4,])
boxdata=as.data.frame(boxdata)
colnames(boxdata)<-c("delta","eps","tvd")
boxdata$eps=as.factor(boxdata$eps)
boxdata$delta=as.factor(boxdata$delta)

library(ggplot2)
x11()
ggplot(boxdata, aes(x=eps, y=tvd, fill=delta)) +
  geom_boxplot(position=position_dodge(1),width=0.4)+
  scale_y_continuous(limits=c(0,0.5))+
  theme_minimal()+
  labs(x=expression(epsilon),y="total variation distance")+
  theme(axis.line = element_line(colour = "black"))+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  scale_fill_discrete(name = "Dose", labels=c(expression(paste(delta,"=0.001")),expression(paste(delta,"=0.01")),expression(paste(delta,"=0.1"))))+
  theme(legend.text.align = 0)+
  theme(legend.text = element_text(size = 15))

# 3-way pDP
simu=80 
eps=10^seq(-3,2,1)
way3=array(NA, dim = c(length(eps), 20, simu))
delta=0.001
c53=combn(seq(1,5,1),3)  #3*10
c52=combn(seq(1,5,1),2)  #2*10
for (k in 1:simu){
  bydim3=matrix(NA,length(eps), 20)
  for (t in 1:length(eps)){
    sanitize=data_new
    sigma=(sqrt(qnorm(delta/2)^2+2*eps[t])-qnorm(delta/2))/(2*eps[t])
    sanitize[,6:8]=sanitize[,6:8]+rnorm(3*78,0,sigma)
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:10){
      pstar=c(ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0))
      por=c(ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0))
      bydim3[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
    for (j in 1:10){
      pstar=c(sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),8]))
      por=c(sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),8]))
      bydim3[t,10+j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
  }    
  way3[,,k]=bydim3
}
TVD.0.001=apply(way3,c(1,2),mean)

way3=array(NA, dim = c(length(eps), 20, simu))
delta=0.01
for (k in 1:simu){
  bydim3=matrix(NA,length(eps), 20)
  for (t in 1:length(eps)){
    sanitize=data_new
    sigma=(sqrt(qnorm(delta/2)^2+2*eps[t])-qnorm(delta/2))/(2*eps[t])
    sanitize[,6:8]=sanitize[,6:8]+rnorm(3*78,0,sigma)
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:10){
      pstar=c(ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0))
      por=c(ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0))
      bydim3[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
    for (j in 1:10){
      pstar=c(sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),8]))
      por=c(sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),8]))
      bydim3[t,10+j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
  }    
  way3[,,k]=bydim3
}
TVD.0.01<-apply(way3,c(1,2),mean)

way3=array(NA, dim = c(length(eps), 20, simu))
delta=0.1
for (k in 1:simu){
  bydim3=matrix(NA,length(eps), 20)
  for (t in 1:length(eps)){
    sanitize=data_new
    sigma=(sqrt(qnorm(delta/2)^2+2*eps[t])-qnorm(delta/2))/(2*eps[t])
    sanitize[,6:8]=sanitize[,6:8]+rnorm(3*78,0,sigma)
    cell=sanitize[,6:8]
    cell[cell<0]=0
    sanitize[,6:8]=round(cell)
    for (j in 1:10){
      pstar=c(ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="A")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="N")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="A")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="N")&(sanitize[,c53[3,j]]=="P")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="A")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="N")),6:8]),0),
              ifelse(length(which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")))>0,sum(sanitize[which((sanitize[,c53[1,j]]=="P")&(sanitize[,c53[2,j]]=="P")&(sanitize[,c53[3,j]]=="P")),6:8]),0))
      por=c(ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="A")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="N")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="A")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="N")&(orig[,c53[3,j]]=="P")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="A")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="N")),6:8]),0),
            ifelse(length(which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")))>0,sum(orig[which((orig[,c53[1,j]]=="P")&(orig[,c53[2,j]]=="P")&(orig[,c53[3,j]]=="P")),6:8]),0))
      bydim3[t,j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
    for (j in 1:10){
      pstar=c(sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="A")&(sanitize[,c52[2,j]]=="P")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="N")&(sanitize[,c52[2,j]]=="P")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="A")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="N")),8]),
              sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),6]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),7]),sum(sanitize[which((sanitize[,c52[1,j]]=="P")&(sanitize[,c52[2,j]]=="P")),8]))
      por=c(sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="A")&(orig[,c52[2,j]]=="P")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="N")&(orig[,c52[2,j]]=="P")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="A")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="N")),8]),
            sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),6]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),7]),sum(orig[which((orig[,c52[1,j]]=="P")&(orig[,c52[2,j]]=="P")),8]))
      bydim3[t,10+j]=(1/2)*sum(abs(pstar/sum(pstar)-por/sum(por)))
    }
  }    
  way3[,,k]=bydim3
}
TVD.0.1<-apply(way3,c(1,2),mean)


boxdata=matrix(NA,360,3)
boxdata[,1]=c(rep(0.001,120),rep(0.01,120),rep(0.1,120))
boxdata[,2]=rep(c(rep(eps[1],20),rep(eps[2],20),rep(eps[3],20),rep(eps[4],20),rep(eps[5],20),rep(eps[6],20)),3)
boxdata[,3]=c(TVD.0.001[1,],TVD.0.001[2,],TVD.0.001[3,],TVD.0.001[4,],TVD.0.001[5,],TVD.0.001[6,],TVD.0.01[1,],TVD.0.01[2,],TVD.0.01[3,],TVD.0.01[4,],TVD.0.01[5,],TVD.0.01[6,],TVD.0.1[1,],TVD.0.1[2,],TVD.0.1[3,],TVD.0.1[4,],TVD.0.1[5,],TVD.0.1[6,])
boxdata=as.data.frame(boxdata)
colnames(boxdata)<-c("delta","eps","tvd")
boxdata$eps=as.factor(boxdata$eps)
boxdata$delta=as.factor(boxdata$delta)

library(ggplot2)
x11()
ggplot(boxdata, aes(x=eps, y=tvd, fill=delta)) +
  geom_boxplot(position=position_dodge(1),width=0.4)+
  scale_y_continuous(limits=c(0,0.5))+
  theme_minimal()+
  labs(x=expression(epsilon),y="total variation distance")+
  theme(axis.line = element_line(colour = "black"))+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  scale_fill_discrete(name = "Dose", labels=c(expression(paste(delta,"=0.001")),expression(paste(delta,"=0.01")),expression(paste(delta,"=0.1"))))+
  theme(legend.text.align = 0)+
  theme(legend.text = element_text(size = 15))
