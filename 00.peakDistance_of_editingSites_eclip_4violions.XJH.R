##########################################################
##2023年2月10日
##因为之前的绘图结果没有考虑strand的信息，所以考虑之后重新画图。
#################################################################
##################################################################
setwd("C:\\Users\\CSIV149\\PLC.lab.docs\\03.larry\\DDX6data.LarrySend\\DDX6_DHX9_DDX21_DDX3_peakDensity_and_peakDistance")
library(ggplot2)
library(export)
##绘制小提琴图
#ddx6
ddx6.xl1 <- read.csv("DDX6_distanceEditing2EclipPeak_XL1_new.strand.txt",header = T)
ddx6.xl2 <- read.csv("DDX6_distanceEditing2EclipPeak_XL2_new.strand.txt",header = T)
ddx6.dist <- rbind(ddx6.xl1,ddx6.xl2)
head(ddx6.dist)
max(ddx6.dist$Dsitance.nt.)
ddx6.dist$sample <- c(rep("DDX6",nrow(ddx6.dist)))
#ddx3x
ddx3x.xl1 <- read.csv("DDX3X_distanceEditing2EclipPeak_XL1_new.strand.txt",header = T)
ddx3x.xl2 <- read.csv("DDX3X_distanceEditing2EclipPeak_XL2_new.strand.txt",header = T)
ddx3x.dist <- rbind(ddx3x.xl1,ddx3x.xl2)
head(ddx3x.dist)
max(ddx3x.dist$Dsitance.nt.)
ddx3x.dist$sample <- c(rep("DDX3X",nrow(ddx3x.dist)))


#dist.all <- rbind(ddx6.dist,ddx3x.dist)
#head(dist.all)
#ggplot(dist.all) + 
#  geom_violin(aes(x=sample, y=Dsitance.nt., fill=group), trim = T) +
#  xlab("class") +
#  theme(legend.position="none") +
#  xlab("")+
#  theme_bw()+geom_boxplot(aes(x=sample, y=Dsitance.nt., fill=group), width=0.1)+
#  scale_y_continuous(trans='log10', breaks = c(0, 10, 1000, 100000,10000000))

#dhx9
dhx9.xl1 <- read.csv("DHX9_distanceEditing2EclipPeak_XL1_new.strand.txt",header = T)
dhx9.xl2 <- read.csv("DHX9_distanceEditing2EclipPeak_XL2_new.strand.txt",header = T)
dhx9.dist <- rbind(dhx9.xl1,dhx9.xl2)
head(dhx9.dist)
max(dhx9.dist$Dsitance.nt.)
dhx9.dist$sample <- c(rep("DHX9",nrow(dhx9.dist)))
#ddx21
ddx21.xl1 <- read.csv("DDX21_distanceEditing2EclipPeak_XL1_new.strand.txt",header = T)
ddx21.xl2 <- read.csv("DDX21_distanceEditing2EclipPeak_XL2_new.strand.txt",header = T)
ddx21.dist <- rbind(ddx21.xl1,ddx21.xl2)
head(ddx21.dist)
max(ddx21.dist$Dsitance.nt.)
ddx21.dist$sample <- c(rep("DDX21",nrow(ddx21.dist)))

#合并
dist.all <- rbind(ddx6.dist,ddx3x.dist,dhx9.dist,ddx21.dist)
head(dist.all)
#绘图
library(ggplot2)
library(export)
#mpg数据集绘制
#mpg$class = with(mpg, reorder(class, hwy, median))
#######################################
ggplot(dist.all) + 
  geom_violin(aes(x=sample, y=Dsitance.nt., fill=group), trim = T) +
  xlab("class") +
  theme(legend.position="none") +
  xlab("")+
  theme_bw()+geom_boxplot(aes(x=sample, y=Dsitance.nt., fill=group), width=0.1)+
  scale_y_continuous(trans='log10', breaks = c(0, 10, 1000, 100000,10000000))
graph2ppt(file="EditingSite2EclipPeakDist_log10_4samples_new.strand.trim_corrected.pptx", width=16, aspectr=1.5)


