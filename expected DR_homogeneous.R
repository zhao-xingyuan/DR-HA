##### Laplace mechanism #####
# Adult
this.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this.dir)
adult=read.csv("adultcount.csv",header=TRUE)
M=nrow(adult)
adultcount=adult[,8]
eps=10^seq(-3,2,0.01)
aDR=rep(NA,length(eps))
for (i in 1:length(eps)){
  aDR[i]=mean((1-0.5*exp(-0.5*eps[i]))*(1-0.5*exp((0.5-adultcount)*eps[i])))
}

# plot  
df<-data.frame(x=eps,y=aDR)
library(VGAM)
library(scales)
library(ggplot2)
x11()
ggplot()+
  geom_line(data=df,aes(x=eps,y=aDR,color="analytical"),size=1)+
  scale_x_continuous(trans = log10_trans(),limits=c(0.001,100.2),breaks=c(0.001,0.01,0.1,1,10,100),labels=expression(10^{-3},10^{-2},10^{-1},1,10^{1},10^{2}))+
  labs(x=expression(epsilon),y="average disclosure risk")+
  theme(axis.line = element_line(colour = "black"))+theme_bw()+
  scale_shape(guide = FALSE)+scale_linetype(guide = FALSE)+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  scale_y_continuous(breaks=seq(0.2,1,0.2),limits=c(0.2,1))+
  theme(legend.title=element_blank(),legend.position = "none",legend.justification = c(0.05,0.9),legend.text = element_text(size=14))

# Bank
bank=read.csv("bankcount.csv",header=TRUE)
bankcount=bank[,8]
N=nrow(bank)
eps=10^seq(-3,2,0.01)
aDR=rep(NA,length(eps))
for (i in 1:length(eps)){
  aDR[i]=mean((1-0.5*exp(-0.5*eps[i]))*(1-0.5*exp((0.5-bankcount)*eps[i])))
}

# plot  
df<-data.frame(x=eps,y=aDR)
library(ggplot2)
x11()
ggplot()+
  geom_line(data=df,aes(x=eps,y=aDR,color="analytical"),size=1)+
  scale_x_continuous(trans = log10_trans(),limits=c(0.001,100.2),breaks=c(0.001,0.01,0.1,1,10,100),labels=expression(10^{-3},10^{-2},10^{-1},1,10^{1},10^{2}))+
  labs(x=expression(epsilon),y="average disclosure risk")+
  theme(axis.line = element_line(colour = "black"))+theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  scale_y_continuous(breaks=seq(0.2,1,0.2),limits=c(0.2,1))+
  theme(legend.title=element_blank(),legend.position = "none",legend.justification = c(0.05,0.9),legend.text = element_text(size=14))

##### Gaussian mechanism #####
# Adult: (eps,delta)-aDP 
eps=10^seq(-3,0,0.01)
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

delta=0.001
DR.0.001=rep(0,length(eps))
for (t in 1:length(eps)){
  sigma=sqrt(2*log(1.25/delta))/eps[t]
  DR.0.001[t]=mean((1+erf(1/(2*sqrt(2)*sigma)))*(1+erf((adultcount-0.5)/(sigma*sqrt(2))))/4)
}

delta=0.01
DR.0.01=rep(0,length(eps))
for (t in 1:length(eps)){
  sigma=sqrt(2*log(1.25/delta))/eps[t]
  DR.0.01[t]=mean((1+erf(1/(2*sqrt(2)*sigma)))*(1+erf((adultcount-0.5)/(sigma*sqrt(2))))/4)
}

delta=0.1
DR.0.1=rep(0,length(eps))
for (t in 1:length(eps)){
  sigma=sqrt(2*log(1.25/delta))/eps[t]
  DR.0.1[t]=mean((1+erf(1/(2*sqrt(2)*sigma)))*(1+erf((adultcount-0.5)/(sigma*sqrt(2))))/4)
}

delta=0.00003
DR.0.00003=rep(0,length(eps))
for (t in 1:length(eps)){
  sigma=sqrt(2*log(1.25/delta))/eps[t]
  DR.0.00003[t]=mean((1+erf(1/(2*sqrt(2)*sigma)))*(1+erf((adultcount-0.5)/(sigma*sqrt(2))))/4)
}

# plot
library(VGAM)
x=eps

df.0.001<-data.frame(x=eps,y=DR.0.001)
df.0.01<-data.frame(x=eps,y=DR.0.01)
df.0.1<-data.frame(x=eps,y=DR.0.1)
df.0.00003<-data.frame(x=eps,y=DR.0.00003)

library(ggplot2)
library(scales)
x11()
ggplot()+
  geom_line(data=df.0.00003,aes(x=x,y=DR.0.00003,color="c4",linetype="l4"),size=1.1)+
  geom_line(data=df.0.001,aes(x=x,y=DR.0.001,color="c1",linetype="l1"),size=1.1)+
  geom_line(data=df.0.01,aes(x=x,y=DR.0.01,color="c2",linetype="l2"),size=1.1)+
  geom_line(data=df.0.1,aes(x=x,y=DR.0.1,color="c3",linetype="l3"),size=1.1)+
  scale_x_continuous(trans = log10_trans(),limits=c(0.001,1.2),breaks=c(0.001,0.01,0.1,1),labels=expression(10^{-3},10^{-2},10^{-1},1))+
  labs(x=expression(epsilon),y="average disclosure risk")+
  theme(axis.line = element_line(colour = "black"))+theme_bw()+
  theme(legend.key.size = unit(1.6,"line"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  scale_y_continuous(breaks=seq(0.2,1,0.1),limits=c(0.2,1))+
  theme(legend.title=element_blank(),legend.position = c(0.05, 0.907),legend.justification = c(0.05,0.9),legend.text = element_text(size=18))+
  scale_colour_manual(breaks=c("c4","c1","c2","c3"),values=c("#00BFC4","#F8766D","#B79F00","#00BA38"),labels=c(expression(paste(delta,"=0.00003")),expression(paste(delta,"=0.001")),expression(paste(delta,"=0.01")),expression(paste(delta,"=0.1"))))+
  scale_linetype_manual(breaks=c("l4","l1","l2","l3"),values=c("dotdash","solid","dotted","dashed"),guide="none")+
  scale_shape_manual(breaks=c("s4","s1","s2","s3"),values=c(4,16,17,15),guide="none")+
  theme(legend.text.align = 0)+
  guides(colour=guide_legend(override.aes = list(color=c("#00BFC4","#F8766D","#B79F00","#00BA38"),linetype=c("dotdash","solid","dotted","dashed"),shape=c(4,16,17,15))))


# Adult: (eps,delta)-pDP
delta=0.001
eps=10^seq(-3,2,0.01)
DR.0.001=rep(0,length(eps))
for (t in 1:length(eps)){
  sigma=(sqrt(qnorm(delta/2)^2+2*eps[t])-qnorm(delta/2))/(2*eps[t])
  DR.0.001[t]=mean((1+erf(1/(2*sqrt(2)*sigma)))*(1+erf((adultcount-0.5)/(sigma*sqrt(2))))/4)
}

delta=0.01
DR.0.01=rep(0,length(eps))
for (t in 1:length(eps)){
  sigma=(sqrt(qnorm(delta/2)^2+2*eps[t])-qnorm(delta/2))/(2*eps[t])
  DR.0.01[t]=mean((1+erf(1/(2*sqrt(2)*sigma)))*(1+erf((adultcount-0.5)/(sigma*sqrt(2))))/4)
}

delta=0.1
DR.0.1=rep(0,length(eps))
for (t in 1:length(eps)){
  sigma=(sqrt(qnorm(delta/2)^2+2*eps[t])-qnorm(delta/2))/(2*eps[t])
  DR.0.1[t]=mean((1+erf(1/(2*sqrt(2)*sigma)))*(1+erf((adultcount-0.5)/(sigma*sqrt(2))))/4)
}

delta=0.00003
DR.0.00003=rep(0,length(eps))
for (t in 1:length(eps)){
  sigma=(sqrt(qnorm(delta/2)^2+2*eps[t])-qnorm(delta/2))/(2*eps[t])
  DR.0.00003[t]=mean((1+erf(1/(2*sqrt(2)*sigma)))*(1+erf((adultcount-0.5)/(sigma*sqrt(2))))/4)
}

# plot
library(VGAM)
x=eps

df.0.001<-data.frame(x=eps,y=DR.0.001)
df.0.01<-data.frame(x=eps,y=DR.0.01)
df.0.1<-data.frame(x=eps,y=DR.0.1)
df.0.00003<-data.frame(x=eps,y=DR.0.00003)

library(ggplot2)
library(scales)
x11()
ggplot()+
  geom_line(data=df.0.00003,aes(x=x,y=DR.0.00003,color="c4",linetype="l4"),size=1.1)+
  geom_line(data=df.0.001,aes(x=x,y=DR.0.001,color="c1",linetype="l1"),size=1.1)+
  geom_line(data=df.0.01,aes(x=x,y=DR.0.01,color="c2",linetype="l2"),size=1.1)+
  geom_line(data=df.0.1,aes(x=x,y=DR.0.1,color="c3",linetype="l3"),size=1.1)+
  scale_x_continuous(trans = log10_trans(),limits=c(0.001,100),breaks=c(0.001,0.01,0.1,1,10,100),labels=expression(10^{-3},10^{-2},10^{-1},1,10^{1},10^{2}))+
  theme(axis.line = element_line(colour = "black"))+theme_bw()+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  scale_y_continuous(breaks=seq(0.2,1,0.1),limits=c(0.2,1))+
  theme(legend.title=element_blank(),legend.position = c(0.05, 0.9),legend.justification = c(0.05,0.9),legend.text = element_text(size=17))+
  scale_colour_manual(breaks=c("c4","c1","c2","c3"),values=c("#00BFC4","#F8766D","#B79F00","#00BA38"),labels=c(expression(paste(delta,"=0.00003")),expression(paste(delta,"=0.001")),expression(paste(delta,"=0.01")),expression(paste(delta,"=0.1"))))+
  scale_linetype_manual(breaks=c("l4","l1","l2","l3"),values=c("dotdash","solid","dotted","dashed"),guide="none")+
  scale_shape_manual(breaks=c("s4","s1","s2","s3"),values=c(4,16,17,15),guide="none")+  
  theme(legend.text.align = 0)+
  guides(colour=guide_legend(override.aes = list(color=c("#00BFC4","#F8766D","#B79F00","#00BA38"),linetype=c("dotdash","solid","dotted","dashed"),shape=c(4,16,17,15))))+
  theme(legend.key.size = unit(1.6,"line"))


# Bank: (eps,delta)-aDP 
eps=10^seq(-3,0,0.01)
erf <- function(x) 2 * pnorm(x * sqrt(2)) - 1

delta=0.001
DR.0.001=rep(0,length(eps))
for (t in 1:length(eps)){
  sigma=sqrt(2*log(1.25/delta))/eps[t]
  DR.0.001[t]=mean((1+erf(1/(2*sqrt(2)*sigma)))*(1+erf((bankcount-0.5)/(sigma*sqrt(2))))/4)
}

delta=0.01
DR.0.01=rep(0,length(eps))
for (t in 1:length(eps)){
  sigma=sqrt(2*log(1.25/delta))/eps[t]
  DR.0.01[t]=mean((1+erf(1/(2*sqrt(2)*sigma)))*(1+erf((bankcount-0.5)/(sigma*sqrt(2))))/4)
}

delta=0.1
DR.0.1=rep(0,length(eps))
for (t in 1:length(eps)){
  sigma=sqrt(2*log(1.25/delta))/eps[t]
  DR.0.1[t]=mean((1+erf(1/(2*sqrt(2)*sigma)))*(1+erf((bankcount-0.5)/(sigma*sqrt(2))))/4)
}

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
  labs(x=expression(epsilon),y="average disclosure risk")+
  theme(axis.line = element_line(colour = "black"))+theme_bw()+
  theme(legend.key.size = unit(1.6,"line"))+
  theme(axis.text=element_text(size=14),axis.title=element_text(size=16))+
  scale_y_continuous(breaks=seq(0.2,1,0.1),limits=c(0.2,1))+
  theme(legend.title=element_blank(),legend.position = c(0.05, 0.907),legend.justification = c(0.05,0.9),legend.text = element_text(size=18))+
  scale_colour_manual(breaks=c("c1","c2","c3"),values=c("#F8766D","#B79F00","#00BA38"),labels=c(expression(paste(delta,"=0.001")),expression(paste(delta,"=0.01")),expression(paste(delta,"=0.1"))))+
  scale_linetype_manual(breaks=c("l1","l2","l3"),values=c("solid","dotted","dashed"),guide="none")+
  scale_shape_manual(breaks=c("s1","s2","s3"),values=c(16,17,15),guide="none")+  
  theme(legend.text.align = 0)+
  guides(colour=guide_legend(override.aes = list(color=c("#F8766D","#B79F00","#00BA38"),linetype=c("solid","dotted","dashed"),shape=c(16,17,15))))

# Bank: (eps,delta)-pDP
delta=0.001
eps=10^seq(-3,2,0.01)
DR.0.001=rep(0,length(eps))
for (t in 1:length(eps)){
  sigma=(sqrt(qnorm(delta/2)^2+2*eps[t])-qnorm(delta/2))/(2*eps[t])
  DR.0.001[t]=mean((1+erf(1/(2*sqrt(2)*sigma)))*(1+erf((bankcount-0.5)/(sigma*sqrt(2))))/4)
}

delta=0.01
DR.0.01=rep(0,length(eps))
for (t in 1:length(eps)){
  sigma=(sqrt(qnorm(delta/2)^2+2*eps[t])-qnorm(delta/2))/(2*eps[t])
  DR.0.01[t]=mean((1+erf(1/(2*sqrt(2)*sigma)))*(1+erf((bankcount-0.5)/(sigma*sqrt(2))))/4)
}

delta=0.1
DR.0.1=rep(0,length(eps))
for (t in 1:length(eps)){
  sigma=(sqrt(qnorm(delta/2)^2+2*eps[t])-qnorm(delta/2))/(2*eps[t])
  DR.0.1[t]=mean((1+erf(1/(2*sqrt(2)*sigma)))*(1+erf((bankcount-0.5)/(sigma*sqrt(2))))/4)
}

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
  scale_y_continuous(breaks=seq(0.2,1,0.1),limits=c(0.2,1))+
  theme(legend.title=element_blank(),legend.position = c(0.05, 0.9),legend.justification = c(0.05,0.9),legend.text = element_text(size=17))+
  scale_colour_manual(breaks=c("c1","c2","c3"),values=c("#F8766D","#B79F00","#00BA38"),labels=c(expression(paste(delta,"=0.001")),expression(paste(delta,"=0.01")),expression(paste(delta,"=0.1"))))+
  scale_linetype_manual(breaks=c("l1","l2","l3"),values=c("solid","dotted","dashed"),guide="none")+
  scale_shape_manual(breaks=c("s1","s2","s3"),values=c(16,17,15),guide="none")+  
  theme(legend.text.align = 0)+
  guides(colour=guide_legend(override.aes = list(color=c("#F8766D","#B79F00","#00BA38"),linetype=c("solid","dotted","dashed"),shape=c(16,17,15))))+
  theme(legend.key.size = unit(1.6,"line"))
