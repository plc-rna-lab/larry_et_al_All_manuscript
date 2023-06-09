library(magrittr)
library(dplyr)
## boxplot for A-U, A-C and non-A:C non-A:U binding maps (-1000 to 1000nt)
setwd("C:/Users/CSIV149/PLC.lab.docs/03.larry/DDX6data.LarrySend/DHX9_A2C.A2U.noCnoU_peakDensity")
WD <- getwd()

df1 = read.delim(paste0(WD, "/upDHX9_A2C_editing_candidates.overlppaped_with_DHX9_eCLIP_IgG.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df2 = read.delim(paste0(WD, "/upDHX9_A2C_editing_candidates.overlppaped_with_DHX9_eCLIP_XL1.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df3 = read.delim(paste0(WD, "/upDHX9_A2C_editing_candidates.overlppaped_with_DHX9_eCLIP_XL2.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)

df1 = read.delim(paste0(WD, "/downDHX9_A2C_editing_candidates.overlppaped_with_DHX9_eCLIP_IgG.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df2 = read.delim(paste0(WD, "/downDHX9_A2C_editing_candidates.overlppaped_with_DHX9_eCLIP_XL1.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df3 = read.delim(paste0(WD, "/downDHX9_A2C_editing_candidates.overlppaped_with_DHX9_eCLIP_XL2.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)

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
# plot(smoothingSpline1, ylim = c(0,1), xlim = c(-5000, 5000), type = "l", col="black", lwd=2, xlab="Relative position to 119 A:C mismatches", ylab="DDX6 eCLIP peak density")
# lines(smoothingSpline2, col="orange", lwd=2)
# lines(smoothingSpline3, col="dark blue", lwd=2)

# plot(df1$V2, df1$density, ylim = c(0,1), xlim = c(-1000, 1000), type = "l", col="black", lwd=2, xlab="Relative position to 119 A:C mismatches", ylab="DDX6 eCLIP peak density")
# lines(df2$V2, df2$density, col="orange", lwd=2)
# lines(df3$V2, df3$density, col="dark blue", lwd=2)
df1_oppo_nucleotide_C = df1
df2_oppo_nucleotide_C = df2
df3_oppo_nucleotide_C = df3





df1 = read.delim(paste0(WD, "/upDHX9_A2U_editing_candidates.overlppaped_with_DHX9_eCLIP_IgG.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df2 = read.delim(paste0(WD, "/upDHX9_A2U_editing_candidates.overlppaped_with_DHX9_eCLIP_XL1.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df3 = read.delim(paste0(WD, "/upDHX9_A2U_editing_candidates.overlppaped_with_DHX9_eCLIP_XL2.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)

df1 = read.delim(paste0(WD, "/downDHX9_A2U_editing_candidates.overlppaped_with_DHX9_eCLIP_IgG.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df2 = read.delim(paste0(WD, "/downDHX9_A2U_editing_candidates.overlppaped_with_DHX9_eCLIP_XL1.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df3 = read.delim(paste0(WD, "/downDHX9_A2U_editing_candidates.overlppaped_with_DHX9_eCLIP_XL2.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)

df1 = df1[order(df1$V2),]
df2 = df2[order(df2$V2),]
df3 = df3[order(df3$V2),]
# up-regulated "143 A:C of XL1", "514 A:U of XL1", "327 non-A:C non-A:U of XL1"
df1 = df1%>% mutate(density = (V3/514*143)/max)
df2 = df2%>% mutate(density = (V3/514*143)/max)
df3 = df3%>% mutate(density = (V3/514*143)/max)

# down-regulated "236 A:C of XL1", "648 A:U of XL1", "418 non-A:C non-A:U of XL1"
df1 = df1%>% mutate(density = (V3/648*236)/max)
df2 = df2%>% mutate(density = (V3/648*236)/max)
df3 = df3%>% mutate(density = (V3/648*236)/max)
smoothingSpline1 = smooth.spline(df1$V2, df1$density, spar=0)
smoothingSpline2 = smooth.spline(df2$V2, df2$density, spar=0)
smoothingSpline3 = smooth.spline(df3$V2, df3$density, spar=0)
# plot(smoothingSpline1, ylim = c(0,1), xlim = c(-1000, 1000), type = "l", col="black", lwd=2, xlab="Relative position to 322 A:U sites", ylab="DDX6 eCLIP peak density")
# lines(smoothingSpline2, col="orange", lwd=2)
# lines(smoothingSpline3, col="dark blue", lwd=2)
df1_oppo_nucleotide_U = df1
df2_oppo_nucleotide_U = df2
df3_oppo_nucleotide_U = df3




df1 = read.delim(paste0(WD, "/upDHX9_non_C_non_U_editing_candidates.overlppaped_with_DHX9_eCLIP_IgG.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df2 = read.delim(paste0(WD, "/upDHX9_non_C_non_U_editing_candidates.overlppaped_with_DHX9_eCLIP_XL1.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df3 = read.delim(paste0(WD, "/upDHX9_non_C_non_U_editing_candidates.overlppaped_with_DHX9_eCLIP_XL2.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)

df1 = read.delim(paste0(WD, "/downDHX9_non_C_non_U_editing_candidates.overlppaped_with_DHX9_eCLIP_IgG.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df2 = read.delim(paste0(WD, "/downDHX9_non_C_non_U_editing_candidates.overlppaped_with_DHX9_eCLIP_XL1.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)
df3 = read.delim(paste0(WD, "/downDHX9_non_C_non_U_editing_candidates.overlppaped_with_DHX9_eCLIP_XL2.output.txt.output"), stringsAsFactors = F, header = F, check.names = F)

df1 = df1[order(df1$V2),]
df2 = df2[order(df2$V2),]
df3 = df3[order(df3$V2),]

# up-regulated "143 A:C of XL1", "514 A:U of XL1", "327 non-A:C non-A:U of XL1"
df1 = df1%>% mutate(density = (V3/327*143)/max)
df2 = df2%>% mutate(density = (V3/327*143)/max)
df3 = df3%>% mutate(density = (V3/327*143)/max)

# down-regulated "236 A:C of XL1", "648 A:U of XL1", "418 non-A:C non-A:U of XL1"
df1 = df1%>% mutate(density = (V3/418*236)/max)
df2 = df2%>% mutate(density = (V3/418*236)/max)
df3 = df3%>% mutate(density = (V3/418*236)/max)
smoothingSpline1 = smooth.spline(df1$V2, df1$density, spar=0)
smoothingSpline2 = smooth.spline(df2$V2, df2$density, spar=0)
smoothingSpline3 = smooth.spline(df3$V2, df3$density, spar=0)
# plot(smoothingSpline1, ylim = c(0,1), xlim = c(-1000, 1000), type = "l", col="black", lwd=2, xlab="Relative position to 178 A to non_C_non_U sites", ylab="DDX6 eCLIP peak density")
# lines(smoothingSpline2, col="orange", lwd=2)
# lines(smoothingSpline3, col="dark blue", lwd=2)
df1_oppo_nucleotide_non_C_non_U = df1
df2_oppo_nucleotide_non_C_non_U = df2
df3_oppo_nucleotide_non_C_non_U = df3


## statistical analysis
df_test = df2_oppo_nucleotide_non_C_non_U[9001:11001,4]
df2_oppo_nucleotide_non_C_non_U = data.frame(Peak_density=df_test)
df2_oppo_nucleotide_non_C_non_U$type = "non A:C non A:U.XL1"
df_test = df3_oppo_nucleotide_non_C_non_U[9001:11001,4]
df3_oppo_nucleotide_non_C_non_U = data.frame(Peak_density=df_test)
df3_oppo_nucleotide_non_C_non_U$type = "non A:C non A:U.XL2"
df_test = df2_oppo_nucleotide_C[9001:11001,4]
df2_oppo_nucleotide_C = data.frame(Peak_density=df_test)
df2_oppo_nucleotide_C$type = "A:C.XL1"
df_test = df3_oppo_nucleotide_C[9001:11001,4]
df3_oppo_nucleotide_C = data.frame(Peak_density=df_test)
df3_oppo_nucleotide_C$type = "A:C.XL2"
df_test = df2_oppo_nucleotide_U[9001:11001,4]
df2_oppo_nucleotide_U = data.frame(Peak_density=df_test)
df2_oppo_nucleotide_U$type = "A:U.XL1"
df_test = df3_oppo_nucleotide_U[9001:11001,4]
df3_oppo_nucleotide_U = data.frame(Peak_density=df_test)
df3_oppo_nucleotide_U$type = "A:U.XL2"


library(ggpubr)
library(export)
df_total = rbind(df2_oppo_nucleotide_C, df2_oppo_nucleotide_U, df2_oppo_nucleotide_non_C_non_U)
my_comparisons <- list(c("A:C.XL1", "A:U.XL1"), c("A:U.XL1", "non A:C non A:U.XL1"), c("A:C.XL1","non A:C non A:U.XL1"))
p = ggboxplot(df_total, x = "type", y = "Peak_density",
              color = "type", palette = "jco")+
  # Add pairwise comparisons p-value
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")   
p + rremove("x.text") + rremove("xlab") + rremove("x.ticks")
graph2ppt(file="upDHX9_A2C.A2U.noCnoU_XL1_boxPlot.pptx",width=8,aspectr=1.2)
graph2ppt(file="downDHX9_A2C.A2U.noCnoU_XL1_boxPlot.pptx",width=8,aspectr=1.2)



df_total = rbind(df3_oppo_nucleotide_C, df3_oppo_nucleotide_U, df3_oppo_nucleotide_non_C_non_U)
my_comparisons <- list(c("A:C.XL2", "A:U.XL2"), c("A:U.XL2", "non A:C non A:U.XL2"), c("A:C.XL2","non A:C non A:U.XL2"))
p = ggboxplot(df_total, x = "type", y = "Peak_density",
              color = "type", palette = "jco")+
  # Add pairwise comparisons p-value
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
p + rremove("x.text") + rremove("xlab") + rremove("x.ticks")
graph2ppt(file="upDHX9_A2C.A2U.noCnoU_XL2_boxPlot.pptx",width=8,aspectr=1.2)
graph2ppt(file="downDHX9_A2C.A2U.noCnoU_XL2_boxPlot.pptx",width=8,aspectr=1.2)
