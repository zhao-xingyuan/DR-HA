rm(list=ls())
this.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this.dir)
FF=read.csv("FFcount.csv",header=TRUE) # 78 cells, 103 obs, has hetero
FF$FFA=ifelse(FF$FF=="A",FF$FF.1,0)
FF$FFN=ifelse(FF$FF=="N",FF$FF.1,0)
FF$FFP=ifelse(FF$FF=="P",FF$FF.1,0)
FF=FF[,-c(6,7)]

unique(FF[,seq(1,5,1)])
nrow(unique(FF[,seq(1,5,1)])) #78 cells

index = rep('a',length(FF$IR))
for(i in 1:length(FF$IR))
{index[i] = paste(FF[i,1:5],collapse="")}

table(index)
data = aggregate(x = FF[,6:8], by = list(index,index,index), FUN = "sum")

homo=rep(NA,78)
for (i in 1:78){
  homo[i]=(sum(data[i,4],data[i,5],data[i,6])==max(data[i,4],data[i,5],data[i,6]))
}
sum(homo)

bincount=rep(0,78)
for (i in 1:78){
  bincount[i]=sum(data[i,4:6])
}
mean(bincount)

# Laplace
ep<-c(seq(0.001,0.999,0.002),seq(1,100,1))
res<-matrix(nrow=78,ncol=length(ep))
p=data[,c(4,5,6)]/apply(data[,c(4,5,6)],1,sum)
n=apply(data[,c(4,5,6)],1,sum)
for (j in 1:length(ep)){
  for (i in 1:78){
    res[i,j]<-ifelse(homo[i]==TRUE,1,0)*((1-0.5*exp(-0.5*ep[j]))^2)*(1-0.5*exp(ep[j]*(0.5-n[i])))+
      ifelse(homo[i]==FALSE,1,0)*(1-0.5*exp(-0.5*ep[j]))*0.5*(exp(-0.5*ep[j])+exp((1.5-n[i])*ep[j])-exp((1-n[i])*ep[j]))
  }
}
x=ep
hetero<-apply(res,2,mean)
df<-data.frame(x=x,y=hetero)

library(VGAM)
library(scales)
library(ggplot2)
x11()
ggplot()+
  geom_line(data=df,aes(x=ep,y=hetero,color="analytical"),size=1)+
  scale_x_continuous(trans = log10_trans(),limits=c(0.001,100),breaks=c(0.001,0.01,0.1,1,10,100),labels=expression(10^{-3},10^{-2},10^{-1},1,10^{1},10^{2}))+
  labs(x=expression(epsilon),y="average disclosure risk")+
  theme(axis.line = element_line(colour = "black"))+theme_bw()+
  scale_shape(guide = FALSE)+scale_linetype(guide = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  scale_y_continuous(breaks=seq(0,1,0.2),limits=c(0,1))+
  theme(legend.title=element_blank(),legend.position = "none",legend.justification = c(0.05,0.9),legend.text = element_text(size=14))

# aDP
eps=10^seq(-3,0,0.01)
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

delta=0.001
DR.0.001=matrix(nrow=78,ncol=length(eps))
for (t in 1:length(eps)){
  sigma=sqrt(2*log(1.25/delta))/eps[t]
  for (i in 1:78){
    DR.0.001[i,t]=(1/8)*((1+erf(1/(2*sqrt(2)*sigma)))^2)*ifelse(homo[i]==TRUE,1,0)*(1+erf((n[i]-0.5)/(sigma*sqrt(2))))+
      (1/8)*ifelse(homo[i]==FALSE,1,0)*(1+erf(1/(2*sqrt(2)*sigma)))*((1-erf((1.5-n[i])/(sigma*sqrt(2))))*(1+erf(-0.5/((sqrt(2)*sigma))))+(1+erf((1.5-n[i])/(sigma*sqrt(2))))*(1+erf(1/(2*sqrt(2)*sigma))))
     }
}
DR.0.001=apply(DR.0.001,2,mean)

delta=0.01
DR.0.01=matrix(nrow=78,ncol=length(eps))
for (t in 1:length(eps)){
  sigma=sqrt(2*log(1.25/delta))/eps[t]
  for (i in 1:78){
    DR.0.01[i,t]=(1/8)*((1+erf(1/(2*sqrt(2)*sigma)))^2)*ifelse(homo[i]==TRUE,1,0)*(1+erf((n[i]-0.5)/(sigma*sqrt(2))))+
      (1/8)*ifelse(homo[i]==FALSE,1,0)*(1+erf(1/(2*sqrt(2)*sigma)))*((1-erf((1.5-n[i])/(sigma*sqrt(2))))*(1+erf(-0.5/((sqrt(2)*sigma))))+(1+erf((1.5-n[i])/(sigma*sqrt(2))))*(1+erf(1/(2*sqrt(2)*sigma))))
  }
}
DR.0.01=apply(DR.0.01,2,mean)

delta=0.1
DR.0.1=matrix(nrow=78,ncol=length(eps))
for (t in 1:length(eps)){
  sigma=sqrt(2*log(1.25/delta))/eps[t]
  for (i in 1:78){
    DR.0.1[i,t]=(1/8)*((1+erf(1/(2*sqrt(2)*sigma)))^2)*ifelse(homo[i]==TRUE,1,0)*(1+erf((n[i]-0.5)/(sigma*sqrt(2))))+
      (1/8)*ifelse(homo[i]==FALSE,1,0)*(1+erf(1/(2*sqrt(2)*sigma)))*((1-erf((1.5-n[i])/(sigma*sqrt(2))))*(1+erf(-0.5/((sqrt(2)*sigma))))+(1+erf((1.5-n[i])/(sigma*sqrt(2))))*(1+erf(1/(2*sqrt(2)*sigma))))
  }
}
DR.0.1=apply(DR.0.1,2,mean)


# plot
library(VGAM)
x=eps

df.0.001<-data.frame(x=eps,y=DR.0.001)
df.0.01<-data.frame(x=eps,y=DR.0.01)
df.0.1<-data.frame(x=eps,y=DR.0.1)

library(ggplot2)
library(scales)
x11()
ggplot()+
  geom_line(data=df.0.001,aes(x=x,y=DR.0.001,color="c1",linetype="l1"),size=1.1)+
  geom_line(data=df.0.01,aes(x=x,y=DR.0.01,color="c2",linetype="l2"),size=1.1)+
  geom_line(data=df.0.1,aes(x=x,y=DR.0.1,color="c3",linetype="l3"),size=1.1)+
  scale_x_continuous(trans = log10_trans(),limits=c(0.001,1.2),breaks=c(0.001,0.01,0.1,1),labels=expression(10^{-3},10^{-2},10^{-1},1))+
  theme(axis.line = element_line(colour = "black"))+theme_bw()+
  theme(legend.key.size = unit(1.6,"line"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  scale_y_continuous(breaks=seq(0,1,0.2),limits=c(0,1))+
  theme(legend.title=element_blank(),legend.position = c(0.05, 0.907),legend.justification = c(0.05,0.9),legend.text = element_text(size=18))+
  scale_colour_manual(breaks=c("c1","c2","c3"),values=c("#F8766D","#B79F00","#00BA38"),labels=c(expression(paste(delta,"=0.001")),expression(paste(delta,"=0.01")),expression(paste(delta,"=0.1"))))+
  scale_linetype_manual(breaks=c("l1","l2","l3"),values=c("solid","dotted","dashed"),guide="none")+
  theme(legend.text.align = 0)+
  guides(colour=guide_legend(override.aes = list(color=c("#F8766D","#B79F00","#00BA38"),linetype=c("solid","dotted","dashed"))))


# Bank: (eps,delta)-pDP
delta=0.001
eps=10^seq(-3,2,0.01)
DR.0.001=matrix(nrow=78,ncol=length(eps))
for (t in 1:length(eps)){
  sigma=(sqrt(qnorm(delta/2)^2+2*eps[t])-qnorm(delta/2))/(2*eps[t])
  for (i in 1:78){
    DR.0.001[i,t]=(1/8)*((1+erf(1/(2*sqrt(2)*sigma)))^2)*ifelse(homo[i]==TRUE,1,0)*(1+erf((n[i]-0.5)/(sigma*sqrt(2))))+
      (1/8)*ifelse(homo[i]==FALSE,1,0)*(1+erf(1/(2*sqrt(2)*sigma)))*((1-erf((1.5-n[i])/(sigma*sqrt(2))))*(1+erf(-0.5/((sqrt(2)*sigma))))+(1+erf((1.5-n[i])/(sigma*sqrt(2))))*(1+erf(1/(2*sqrt(2)*sigma))))
  }
}
DR.0.001=apply(DR.0.001,2,mean)


delta=0.01
DR.0.01=matrix(nrow=78,ncol=length(eps))
for (t in 1:length(eps)){
  sigma=(sqrt(qnorm(delta/2)^2+2*eps[t])-qnorm(delta/2))/(2*eps[t])
  for (i in 1:78){
    DR.0.01[i,t]=(1/8)*((1+erf(1/(2*sqrt(2)*sigma)))^2)*ifelse(homo[i]==TRUE,1,0)*(1+erf((n[i]-0.5)/(sigma*sqrt(2))))+
      (1/8)*ifelse(homo[i]==FALSE,1,0)*(1+erf(1/(2*sqrt(2)*sigma)))*((1-erf((1.5-n[i])/(sigma*sqrt(2))))*(1+erf(-0.5/((sqrt(2)*sigma))))+(1+erf((1.5-n[i])/(sigma*sqrt(2))))*(1+erf(1/(2*sqrt(2)*sigma))))
  }
}
DR.0.01=apply(DR.0.01,2,mean)


delta=0.1
DR.0.1=matrix(nrow=78,ncol=length(eps))
for (t in 1:length(eps)){
  sigma=(sqrt(qnorm(delta/2)^2+2*eps[t])-qnorm(delta/2))/(2*eps[t])
  for (i in 1:78){
    DR.0.1[i,t]=(1/8)*((1+erf(1/(2*sqrt(2)*sigma)))^2)*ifelse(homo[i]==TRUE,1,0)*(1+erf((n[i]-0.5)/(sigma*sqrt(2))))+
      (1/8)*ifelse(homo[i]==FALSE,1,0)*(1+erf(1/(2*sqrt(2)*sigma)))*((1-erf((1.5-n[i])/(sigma*sqrt(2))))*(1+erf(-0.5/((sqrt(2)*sigma))))+(1+erf((1.5-n[i])/(sigma*sqrt(2))))*(1+erf(1/(2*sqrt(2)*sigma))))
  }
}
DR.0.1=apply(DR.0.1,2,mean)

# plot
library(VGAM)
x=eps

df.0.001<-data.frame(x=eps,y=DR.0.001)
df.0.01<-data.frame(x=eps,y=DR.0.01)
df.0.1<-data.frame(x=eps,y=DR.0.1)

library(ggplot2)
library(scales)
x11()
ggplot()+
  geom_line(data=df.0.001,aes(x=x,y=DR.0.001,color="c1",linetype="l1"),size=1.1)+
  geom_line(data=df.0.01,aes(x=x,y=DR.0.01,color="c2",linetype="l2"),size=1.1)+
  geom_line(data=df.0.1,aes(x=x,y=DR.0.1,color="c3",linetype="l3"),size=1.1)+
  scale_x_continuous(trans = log10_trans(),limits=c(0.001,100),breaks=c(0.001,0.01,0.1,1,10,100),labels=expression(10^{-3},10^{-2},10^{-1},1,10^{1},10^{2}))+
  labs(x=expression(epsilon),y="average disclosure risk")+
  theme(axis.line = element_line(colour = "black"))+theme_bw()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  scale_y_continuous(breaks=seq(0,1,0.2),limits=c(0,1))+
  theme(legend.title=element_blank(),legend.position = c(0.05, 0.9),legend.justification = c(0.05,0.9),legend.text = element_text(size=17))+
  scale_colour_manual(breaks=c("c1","c2","c3"),values=c("#F8766D","#B79F00","#00BA38"),labels=c(expression(paste(delta,"=0.001")),expression(paste(delta,"=0.01")),expression(paste(delta,"=0.1"))))+
  scale_linetype_manual(breaks=c("l1","l2","l3"),values=c("solid","dotted","dashed"),guide="none")+
  scale_shape_manual(breaks=c("s1","s2","s3"),values=c(16,17,15),guide="none")+  
  theme(legend.text.align = 0)+
  guides(colour=guide_legend(override.aes = list(color=c("#F8766D","#B79F00","#00BA38"),linetype=c("solid","dotted","dashed"),shape=c(16,17,15))))+
  theme(legend.key.size = unit(1.6,"line"))
