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
smoothingSpline1 = smooth.spline(df1$V2, df1$density, spar=0.15)
smoothingSpline2 = smooth.spline(df2$V2, df2$density, spar=0.15)
smoothingSpline3 = smooth.spline(df3$V2, df3$density, spar=0.15)


WD <- getwd()
df1 = read.delim(paste0(WD, "/DDX6_unregulated_editing_candidates.randomly_selected.averaged_20.overlppaped_with_eCLIP_IgG.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df2 = read.delim(paste0(WD, "/DDX6_unregulated_editing_candidates.randomly_selected.averaged_20.overlppaped_with_eCLIP_XL1.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df3 = read.delim(paste0(WD, "/DDX6_unregulated_editing_candidates.randomly_selected.averaged_20.overlppaped_with_eCLIP_XL2.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)

df1 = df1[order(df1$V2),]
df2 = df2[order(df2$V2),]
df3 = df3[order(df3$V2),]

df1 = df1%>% mutate(density = V3/max)
df2 = df2%>% mutate(density = V3/max)
df3 = df3%>% mutate(density = V3/max)

smoothingSpline1.un = smooth.spline(df1$V2, df1$density, spar=0.15)
smoothingSpline2.un = smooth.spline(df2$V2, df2$density, spar=0.15)
smoothingSpline3.un = smooth.spline(df3$V2, df3$density, spar=0.15)


#XL1
plot(smoothingSpline2, ylim = c(0,1), xlim = c(-4000, 4000), 
     type = "l", col="#00BFC4", lwd=3, 
     xlab="Relative position to editing sites", 
     ylab="DDX6 eCLIP peak density")
lines(smoothingSpline2.un, col="#F8766D", lwd=3)

legend(-1000,1,
       c("1095 DDX6 regulated sites (XL1)", "1095 non-DDX6 regulated sites (XL1)"), 
       lwd=c(3,3), 
       col=c("#00BFC4", "#F8766D"), 
       box.lty = 0)


#XL2
plot(smoothingSpline3, ylim = c(0,1), xlim = c(-4000, 4000), 
     type = "l", col="#00BFC4", lwd=3, 
     xlab="Relative position to editing sites", 
     ylab="DDX6 eCLIP peak density")
lines(smoothingSpline3.un, col="#F8766D", lwd=3)

legend(-1000,1,
       c("1095 DDX6 regulated sites (XL1)", "1095 non-DDX6 regulated sites (XL1)"), 
       lwd=c(3,3), 
       col=c("#00BFC4", "#F8766D"), 
       box.lty = 0)









