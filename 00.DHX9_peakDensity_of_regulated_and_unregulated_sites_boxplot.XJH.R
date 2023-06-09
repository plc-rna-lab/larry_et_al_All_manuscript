library(magrittr)
library(dplyr)
## plotting 119 A:C, 322 A:U and 184 non-A:C non-A:U sites together in one figure
# A-U and A-C binding maps ploting with density plot
setwd("C:/Users/CSIV149/PLC.lab.docs/03.larry/DDX6_un-regulated_EditingSites_RX")

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

DDX6_editing_candidates_XL1 <- df2
DDX6_editing_candidates_XL2 <- df3

#UN-REGULATED
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

DDX6_unregulated_candidates_XL1 <- df2
DDX6_unregulated_candidates_XL2 <- df3

## statistical analysis
df_test = DDX6_editing_candidates_XL1[9001:11001,4]
DDX6_editing_candidates_XL1 = data.frame(Peak_density=df_test)
DDX6_editing_candidates_XL1$type = "DDX6 editing sites (XL1)"

df_test = DDX6_editing_candidates_XL2[9001:11001,4]
DDX6_editing_candidates_XL2 = data.frame(Peak_density=df_test)
DDX6_editing_candidates_XL2$type = "DDX6 editing sites (XL2)"

df_test = DDX6_unregulated_candidates_XL1[9001:11001,4]
DDX6_unregulated_candidates_XL1 = data.frame(Peak_density=df_test)
DDX6_unregulated_candidates_XL1$type = "DDX6 unregulated sites (XL1)"

df_test = DDX6_unregulated_candidates_XL2[9001:11001,4]
DDX6_unregulated_candidates_XL2 = data.frame(Peak_density=df_test)
DDX6_unregulated_candidates_XL2$type = "DDX6 unregulated sites (XL2)"



library(ggpubr)
library(export)
df_total = rbind(DDX6_editing_candidates_XL1, DDX6_unregulated_candidates_XL1)
my_comparisons <- list(c("DDX6 editing sites (XL1)", "DDX6 unregulated sites (XL1)"))
p = ggboxplot(df_total, x = "type", y = "Peak_density",
              color = "type", palette = "jco")+
  # Add pairwise comparisons p-value
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")   
p + rremove("x.text") + rremove("xlab") + rremove("x.ticks")
graph2ppt(file="DDX6_regulated_and_unregulated_sites_XL1_boxPlot.pptx",width=8,aspectr=1.2)


df_total = rbind(DDX6_editing_candidates_XL2, DDX6_unregulated_candidates_XL2)
my_comparisons <- list(c("DDX6 editing sites (XL2)", "DDX6 unregulated sites (XL2)"))
p = ggboxplot(df_total, x = "type", y = "Peak_density",
              color = "type", palette = "jco")+
  # Add pairwise comparisons p-value
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")   
p + rremove("x.text") + rremove("xlab") + rremove("x.ticks")
graph2ppt(file="DDX6_regulated_and_unregulated_sites_XL2_boxPlot.pptx",width=8,aspectr=1.2)


####
##
#DHX9
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
smoothingSpline1 = smooth.spline(df1$V2, df1$density, spar=0.15)
smoothingSpline2 = smooth.spline(df2$V2, df2$density, spar=0.15)
smoothingSpline3 = smooth.spline(df3$V2, df3$density, spar=0.15)

DHX9_editing_candidates_XL1 <- df2
DHX9_editing_candidates_XL2 <- df3

#DHX9 un-regulated
WD <- getwd()
df1 = read.delim(paste0(WD, "/DHX9_unregulated_editing_candidates.randomly_selected.averaged_20.overlppaped_with_eCLIP_IgG.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df2 = read.delim(paste0(WD, "/DHX9_unregulated_editing_candidates.randomly_selected.averaged_20.overlppaped_with_eCLIP_XL1.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df3 = read.delim(paste0(WD, "/DHX9_unregulated_editing_candidates.randomly_selected.averaged_20.overlppaped_with_eCLIP_XL2.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)

df1 = df1[order(df1$V2),]
df2 = df2[order(df2$V2),]
df3 = df3[order(df3$V2),]

df1 = df1%>% mutate(density = V3/max)
df2 = df2%>% mutate(density = V3/max)
df3 = df3%>% mutate(density = V3/max)

smoothingSpline1.un = smooth.spline(df1$V2, df1$density, spar=0.15)
smoothingSpline2.un = smooth.spline(df2$V2, df2$density, spar=0.15)
smoothingSpline3.un = smooth.spline(df3$V2, df3$density, spar=0.15)

DHX9_unregulated_candidates_XL1 <- df2
DHX9_unregulated_candidates_XL2 <- df3

## statistical analysis
df_test = DHX9_editing_candidates_XL1[9001:11001,4]
DHX9_editing_candidates_XL1 = data.frame(Peak_density=df_test)
DHX9_editing_candidates_XL1$type = "DHX9 editing sites (XL1)"

df_test = DHX9_editing_candidates_XL2[9001:11001,4]
DHX9_editing_candidates_XL2 = data.frame(Peak_density=df_test)
DHX9_editing_candidates_XL2$type = "DHX9 editing sites (XL2)"

df_test = DHX9_unregulated_candidates_XL1[9001:11001,4]
DHX9_unregulated_candidates_XL1 = data.frame(Peak_density=df_test)
DHX9_unregulated_candidates_XL1$type = "DHX9 unregulated sites (XL1)"

df_test = DHX9_unregulated_candidates_XL2[9001:11001,4]
DHX9_unregulated_candidates_XL2 = data.frame(Peak_density=df_test)
DHX9_unregulated_candidates_XL2$type = "DHX9 unregulated sites (XL2)"



library(ggpubr)
library(export)
df_total = rbind(DHX9_editing_candidates_XL1, DHX9_unregulated_candidates_XL1)
my_comparisons <- list(c("DHX9 editing sites (XL1)", "DHX9 unregulated sites (XL1)"))
p = ggboxplot(df_total, x = "type", y = "Peak_density",
              color = "type", palette = "jco")+
  # Add pairwise comparisons p-value
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")   
p + rremove("x.text") + rremove("xlab") + rremove("x.ticks")
graph2ppt(file="DHX9_regulated_and_unregulated_sites_XL1_boxPlot.pptx",width=8,aspectr=1.2)


df_total = rbind(DHX9_editing_candidates_XL2, DHX9_unregulated_candidates_XL2)
my_comparisons <- list(c("DHX9 editing sites (XL2)", "DHX9 unregulated sites (XL2)"))
p = ggboxplot(df_total, x = "type", y = "Peak_density",
              color = "type", palette = "jco")+
  # Add pairwise comparisons p-value
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")   
p + rremove("x.text") + rremove("xlab") + rremove("x.ticks")
graph2ppt(file="DHX9_regulated_and_unregulated_sites_XL2_boxPlot.pptx",width=8,aspectr=1.2)


