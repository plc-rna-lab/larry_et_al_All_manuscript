library(magrittr)
library(dplyr)
## plotting 119 A:C, 322 A:U and 184 non-A:C non-A:U sites together in one figure
# A-U and A-C binding maps ploting with density plot


WD <- getwd()
df1 = read.delim(paste0(WD, "/DDX6_editing_candidates.overlppaped_with_eCLIP_IgG.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df2 = read.delim(paste0(WD, "/DDX6_editing_candidates.overlppaped_with_eCLIP_XL1.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df3 = read.delim(paste0(WD, "/DDX6_editing_candidates.overlppaped_with_eCLIP_XL2.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df1 = df1[order(df1$V2),]
df2 = df2[order(df2$V2),]
df3 = df3[order(df3$V2),]
df123 = rbind(df1,df2,df3)
max = max(df123$V3)

df1 = df1%>% mutate(density = V3/max)
df2 = df2%>% mutate(density = V3/max)
df3 = df3%>% mutate(density = V3/max)
df_mean = as.data.frame(cbind(df1$V1, df1$V2, df2$density, df3$density))
df_mean$mean_density = (as.numeric(df_mean$V3) + as.numeric(df_mean$V4))/2

smoothingSpline1 = smooth.spline(df1$V2, df1$density, spar=0.15)
smoothingSpline2 = smooth.spline(df2$V2, df2$density, spar=0.15)
smoothingSpline3 = smooth.spline(df3$V2, df3$density, spar=0.15)
smoothingSpline_mean = smooth.spline(df_mean$V2, df_mean$mean_density, spar=0.15)

plot(smoothingSpline1, ylim = c(0,1), xlim = c(-1000, 1000), 
     type = "l", col="black", lwd=3, 
     xlab="Relative position to editing sites", 
     ylab="DDX6 eCLIP peak density")
lines(smoothingSpline2, col="#FF9900", lwd=2, lty=2)
lines(smoothingSpline3, col="#0000FF", lwd=2, lty=2)
lines(smoothingSpline_mean, col="red", lwd=3)

legend(500,1,
       c("IgG", "XL1", "XL2","mean"), 
       lwd=c(2,2,2,2), 
       col=c("black", "#FF9900", "#0000FF","red"), 
       box.lty = 0)

graph2ppt(file="DDX6_peakDensity_of_editingSites_eclip_1Knt.pptx", width=8, aspectr=1.0)




WD <- getwd()
df1 = read.delim(paste0(WD, "/DDX3X_editing_candidates.overlppaped_with_eCLIP_IgG.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df2 = read.delim(paste0(WD, "/DDX3X_editing_candidates.overlppaped_with_eCLIP_XL1.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df3 = read.delim(paste0(WD, "/DDX3X_editing_candidates.overlppaped_with_eCLIP_XL2.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df1 = df1[order(df1$V2),]
df2 = df2[order(df2$V2),]
df3 = df3[order(df3$V2),]
df123 = rbind(df1,df2,df3)
max = max(df123$V3)

df1 = df1%>% mutate(density = V3/max)
df2 = df2%>% mutate(density = V3/max)
df3 = df3%>% mutate(density = V3/max)
df_mean = as.data.frame(cbind(df1$V1, df1$V2, df2$density, df3$density))
df_mean$mean_density = (as.numeric(df_mean$V3) + as.numeric(df_mean$V4))/2

smoothingSpline1 = smooth.spline(df1$V2, df1$density, spar=0.15)
smoothingSpline2 = smooth.spline(df2$V2, df2$density, spar=0.15)
smoothingSpline3 = smooth.spline(df3$V2, df3$density, spar=0.15)
smoothingSpline_mean = smooth.spline(df_mean$V2, df_mean$mean_density, spar=0.15)

plot(smoothingSpline1, ylim = c(0,1), xlim = c(-1000, 1000), 
     type = "l", col="black", lwd=3, 
     xlab="Relative position to editing sites", 
     ylab="DDX3X eCLIP peak density")
lines(smoothingSpline2, col="#FF9900", lwd=2, lty=2)
lines(smoothingSpline3, col="#0000FF", lwd=2, lty=2)
lines(smoothingSpline_mean, col="red", lwd=3)

legend(500,1,
       c("IgG", "XL1", "XL2","mean"), 
       lwd=c(2,2,2,2), 
       col=c("black", "#FF9900", "#0000FF","red"), 
       box.lty = 0)
graph2ppt(file="DDX3X_peakDensity_of_editingSites_eclip_1Knt.pptx", width=8, aspectr=1.0)





WD <- getwd()
df1 = read.delim(paste0(WD, "/DDX21_editing_candidates.overlppaped_with_eCLIP_IgG.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df2 = read.delim(paste0(WD, "/DDX21_editing_candidates.overlppaped_with_eCLIP_XL1.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df3 = read.delim(paste0(WD, "/DDX21_editing_candidates.overlppaped_with_eCLIP_XL2.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df1 = df1[order(df1$V2),]
df2 = df2[order(df2$V2),]
df3 = df3[order(df3$V2),]
df123 = rbind(df1,df2,df3)
max = max(df123$V3)

df1 = df1%>% mutate(density = V3/max)
df2 = df2%>% mutate(density = V3/max)
df3 = df3%>% mutate(density = V3/max)
df_mean = as.data.frame(cbind(df1$V1, df1$V2, df2$density, df3$density))
df_mean$mean_density = (as.numeric(df_mean$V3) + as.numeric(df_mean$V4))/2

smoothingSpline1 = smooth.spline(df1$V2, df1$density, spar=0.15)
smoothingSpline2 = smooth.spline(df2$V2, df2$density, spar=0.15)
smoothingSpline3 = smooth.spline(df3$V2, df3$density, spar=0.15)
smoothingSpline_mean = smooth.spline(df_mean$V2, df_mean$mean_density, spar=0.15)

plot(smoothingSpline1, ylim = c(0,1), xlim = c(-1000, 1000), 
     type = "l", col="black", lwd=3, 
     xlab="Relative position to editing sites", 
     ylab="DDX21 eCLIP peak density")
lines(smoothingSpline2, col="#FF9900", lwd=2, lty=2)
lines(smoothingSpline3, col="#0000FF", lwd=2, lty=2)
lines(smoothingSpline_mean, col="red", lwd=3)

legend(500,1,
       c("IgG", "XL1", "XL2","mean"), 
       lwd=c(2,2,2,2), 
       col=c("black", "#FF9900", "#0000FF","red"), 
       box.lty = 0)
graph2ppt(file="DDX21_peakDensity_of_editingSites_eclip_1Knt.pptx", width=8, aspectr=1.0)






WD <- getwd()
df1 = read.delim(paste0(WD, "/DHX9_editing_candidates.overlppaped_with_eCLIP_IgG.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df2 = read.delim(paste0(WD, "/DHX9_editing_candidates.overlppaped_with_eCLIP_XL1.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df3 = read.delim(paste0(WD, "/DHX9_editing_candidates.overlppaped_with_eCLIP_XL2.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df1 = df1[order(df1$V2),]
df2 = df2[order(df2$V2),]
df3 = df3[order(df3$V2),]
df123 = rbind(df1,df2,df3)
max = max(df123$V3)

df1 = df1%>% mutate(density = V3/max)
df2 = df2%>% mutate(density = V3/max)
df3 = df3%>% mutate(density = V3/max)
df_mean = as.data.frame(cbind(df1$V1, df1$V2, df2$density, df3$density))
df_mean$mean_density = (as.numeric(df_mean$V3) + as.numeric(df_mean$V4))/2

smoothingSpline1 = smooth.spline(df1$V2, df1$density, spar=0.15)
smoothingSpline2 = smooth.spline(df2$V2, df2$density, spar=0.15)
smoothingSpline3 = smooth.spline(df3$V2, df3$density, spar=0.15)
smoothingSpline_mean = smooth.spline(df_mean$V2, df_mean$mean_density, spar=0.15)

plot(smoothingSpline1, ylim = c(0,1), xlim = c(-1000, 1000), 
     type = "l", col="black", lwd=3, 
     xlab="Relative position to editing sites", 
     ylab="DHX9 eCLIP peak density")
lines(smoothingSpline2, col="#FF9900", lwd=2, lty=2)
lines(smoothingSpline3, col="#0000FF", lwd=2, lty=2)
lines(smoothingSpline_mean, col="red", lwd=3)

legend(500,1,
       c("IgG", "XL1", "XL2","mean"), 
       lwd=c(2,2,2,2), 
       col=c("black", "#FF9900", "#0000FF","red"), 
       box.lty = 0)
graph2ppt(file="DHX9_peakDensity_of_editingSites_eclip_1Knt.pptx", width=8, aspectr=1.0)




