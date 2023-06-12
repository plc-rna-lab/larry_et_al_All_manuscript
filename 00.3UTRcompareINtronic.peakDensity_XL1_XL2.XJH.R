library(magrittr)
library(dplyr)
## plotting 119 A:C, 322 A:U and 184 non-A:C non-A:U sites together in one figure
# A-U and A-C binding maps ploting with density plot
setwd("C:/Users/CSIV149/PLC.lab.docs/03.larry/DDX6data.LarrySend/XJH_result_usingRXscript")
WD <- getwd()

df1 = read.delim(paste0(WD, "/DDX6_editing_candidates.overlppaped_with_DDX6_eCLIP_IgG_3UTR.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df2 = read.delim(paste0(WD, "/DDX6_editing_candidates.overlppaped_with_DDX6_eCLIP_XL1_3UTR.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df3 = read.delim(paste0(WD, "/DDX6_editing_candidates.overlppaped_with_DDX6_eCLIP_XL2_3UTR.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df1 = df1[order(df1$V2),]
df2 = df2[order(df2$V2),]
df3 = df3[order(df3$V2),]
df123 = rbind(df1,df2,df3)
max = max(df123$V3)

df1 = df1%>% mutate(density = V3/max)
df2 = df2%>% mutate(density = V3/max)
df3 = df3%>% mutate(density = V3/max)
smoothingSpline1 = smooth.spline(df1$V2, df1$density, spar=0.15)
smoothingSpline2 = smooth.spline(df2$V2, df2$density, spar=0.15)
smoothingSpline3 = smooth.spline(df3$V2, df3$density, spar=0.15)


DDX6_eCLIP_XL1_3UTR <- smoothingSpline2
DDX6_eCLIP_XL2_3UTR <- smoothingSpline3

df1 = read.delim(paste0(WD, "/DDX6_editing_candidates.overlppaped_with_DDX6_eCLIP_IgG_intro.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df2 = read.delim(paste0(WD, "/DDX6_editing_candidates.overlppaped_with_DDX6_eCLIP_XL1_intro.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df3 = read.delim(paste0(WD, "/DDX6_editing_candidates.overlppaped_with_DDX6_eCLIP_XL2_intro.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df1 = df1[order(df1$V2),]
df2 = df2[order(df2$V2),]
df3 = df3[order(df3$V2),]
df123 = rbind(df1,df2,df3)
max = max(df123$V3)

df1 = df1%>% mutate(density = V3/max)
df2 = df2%>% mutate(density = V3/max)
df3 = df3%>% mutate(density = V3/max)
smoothingSpline1 = smooth.spline(df1$V2, df1$density, spar=0.15)
smoothingSpline2 = smooth.spline(df2$V2, df2$density, spar=0.15)
smoothingSpline3 = smooth.spline(df3$V2, df3$density, spar=0.15)

DDX6_eCLIP_XL1_intro <- smoothingSpline2
DDX6_eCLIP_XL2_intro <- smoothingSpline3


#XL1
plot(DDX6_eCLIP_XL1_3UTR, ylim = c(0,1), xlim = c(-4000, 4000), 
     type = "l", col="dark blue", lwd=2, 
     xlab="", ylab="Peak density")
lines(DDX6_eCLIP_XL1_intro, col="orange", lwd=2)
legend(-1000,1,
       c("3'UTR sites(XL1)", "Intronic sites(XL1)"), 
       lwd=c(2,2), col=c("dark blue", "orange"), box.lty = 0)

#XL2
plot(DDX6_eCLIP_XL2_3UTR, ylim = c(0,1), xlim = c(-4000, 4000), 
     type = "l", col="dark blue", lwd=2, 
     xlab="", ylab="Peak density")
lines(DDX6_eCLIP_XL2_intro, col="orange", lwd=2)
legend(-1000,1,
       c("3'UTR sites(XL2)", "Intronic sites(XL2)"), 
       lwd=c(2,2), col=c("dark blue", "orange"), box.lty = 0)

graph2ppt(file="DDX6_3UTRcompareINtronic.peakDensity_XL1.pptx",width=8,aspectr=1.2)
graph2ppt(file="DDX6_3UTRcompareINtronic.peakDensity_XL2.pptx",width=8,aspectr=1.2)

############################################################################################
######################################################################
###DHX9 /// DHX9 /// DHX9
setwd("C:/Users/CSIV149/PLC.lab.docs/03.larry/DDX6data.LarrySend/XJH_result_usingRXscript")
WD <- getwd()

df1 = read.delim(paste0(WD, "/DHX9_editing_candidates.overlppaped_with_DHX9_eCLIP_IgG_3UTR.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df2 = read.delim(paste0(WD, "/DHX9_editing_candidates.overlppaped_with_DHX9_eCLIP_XL1_3UTR.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df3 = read.delim(paste0(WD, "/DHX9_editing_candidates.overlppaped_with_DHX9_eCLIP_XL2_3UTR.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df1 = df1[order(df1$V2),]
df2 = df2[order(df2$V2),]
df3 = df3[order(df3$V2),]
df123 = rbind(df1,df2,df3)
max = max(df123$V3)

df1 = df1%>% mutate(density = V3/max)
df2 = df2%>% mutate(density = V3/max)
df3 = df3%>% mutate(density = V3/max)
smoothingSpline1 = smooth.spline(df1$V2, df1$density, spar=0.15)
smoothingSpline2 = smooth.spline(df2$V2, df2$density, spar=0.15)
smoothingSpline3 = smooth.spline(df3$V2, df3$density, spar=0.15)


DHX9_eCLIP_XL1_3UTR <- smoothingSpline2
DHX9_eCLIP_XL2_3UTR <- smoothingSpline3

df1 = read.delim(paste0(WD, "/DHX9_editing_candidates.overlppaped_with_DHX9_eCLIP_IgG_intro.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df2 = read.delim(paste0(WD, "/DHX9_editing_candidates.overlppaped_with_DHX9_eCLIP_XL1_intro.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df3 = read.delim(paste0(WD, "/DHX9_editing_candidates.overlppaped_with_DHX9_eCLIP_XL2_intro.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df1 = df1[order(df1$V2),]
df2 = df2[order(df2$V2),]
df3 = df3[order(df3$V2),]
df123 = rbind(df1,df2,df3)
max = max(df123$V3)

df1 = df1%>% mutate(density = V3/max)
df2 = df2%>% mutate(density = V3/max)
df3 = df3%>% mutate(density = V3/max)
smoothingSpline1 = smooth.spline(df1$V2, df1$density, spar=0.15)
smoothingSpline2 = smooth.spline(df2$V2, df2$density, spar=0.15)
smoothingSpline3 = smooth.spline(df3$V2, df3$density, spar=0.15)

DHX9_eCLIP_XL1_intro <- smoothingSpline2
DHX9_eCLIP_XL2_intro <- smoothingSpline3


#XL1
plot(DHX9_eCLIP_XL1_3UTR, ylim = c(0,1), xlim = c(-4000, 4000), 
     type = "l", col="dark blue", lwd=2, 
     xlab="", ylab="Peak density")
lines(DHX9_eCLIP_XL1_intro, col="orange", lwd=2)
legend(-1000,1,
       c("3'UTR sites(XL1)", "Intronic sites(XL1)"), 
       lwd=c(2,2), col=c("dark blue", "orange"), box.lty = 0)

#XL2
plot(DHX9_eCLIP_XL2_3UTR, ylim = c(0,1), xlim = c(-4000, 4000), 
     type = "l", col="dark blue", lwd=2, 
     xlab="", ylab="Peak density")
lines(DHX9_eCLIP_XL2_intro, col="orange", lwd=2)
legend(-1000,1,
       c("3'UTR sites(XL2)", "Intronic sites(XL2)"), 
       lwd=c(2,2), col=c("dark blue", "orange"), box.lty = 0)






