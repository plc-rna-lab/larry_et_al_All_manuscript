library(magrittr)
library(dplyr)
## plotting 119 A:C, 322 A:U and 184 non-A:C non-A:U sites together in one figure
# A-U and A-C binding maps ploting with density plot
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
# plot(smoothingSpline1, ylim = c(0,1), xlim = c(-1000, 1000), type = "l", col="black", lwd=2, xlab="Relative position to 119 A:C mismatches", ylab="DDX6 eCLIP peak density")
# lines(smoothingSpline2, col="orange", lwd=2)
# lines(smoothingSpline3, col="dark blue", lwd=2)

DDX6_eCLIP_XL1_3UTR <- df2
DDX6_eCLIP_XL2_3UTR <- df3

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

DDX6_eCLIP_XL1_intro <- df2
DDX6_eCLIP_XL2_intro <- df3


df_test = DDX6_eCLIP_XL1_3UTR[,4]
df_test = DDX6_eCLIP_XL1_3UTR[6001:14001,4]
DDX6_eCLIP_XL1_3UTR = data.frame(Peak_density=df_test)
DDX6_eCLIP_XL1_3UTR$type = "3UTR_sites.XL1"
df_test = DDX6_eCLIP_XL2_3UTR[6001:14001,4]
DDX6_eCLIP_XL2_3UTR = data.frame(Peak_density=df_test)
DDX6_eCLIP_XL2_3UTR$type = "3UTR_sites.XL2"


df_test = DDX6_eCLIP_XL1_intro[6001:14001,4]
DDX6_eCLIP_XL1_intro = data.frame(Peak_density=df_test)
DDX6_eCLIP_XL1_intro$type = "Intro_sites.XL1"
df_test = DDX6_eCLIP_XL2_intro[6001:14001,4]
DDX6_eCLIP_XL2_intro = data.frame(Peak_density=df_test)
DDX6_eCLIP_XL2_intro$type = "Intro_sites.XL2"


df_total = rbind(DDX6_eCLIP_XL1_3UTR, DDX6_eCLIP_XL1_intro)
my_comparisons <- list(c("Intro_sites.XL1", "3UTR_sites.XL1"))
p = ggboxplot(df_total, x = "type", y = "Peak_density",
              color = "type", palette = "jco")+
  # Add pairwise comparisons p-value
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")   
p + rremove("x.text") + rremove("xlab") + rremove("x.ticks")

df_total = rbind(DDX6_eCLIP_XL2_3UTR, DDX6_eCLIP_XL2_intro)
my_comparisons <- list(c("Intro_sites.XL2", "3UTR_sites.XL2"))
p = ggboxplot(df_total, x = "type", y = "Peak_density",
              color = "type", palette = "jco")+
  # Add pairwise comparisons p-value
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
p + rremove("x.text") + rremove("xlab") + rremove("x.ticks")




########################################################################################
######################################################################
###DHX9 /// DHX9 /// DHX9

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
# plot(smoothingSpline1, ylim = c(0,1), xlim = c(-1000, 1000), type = "l", col="black", lwd=2, xlab="Relative position to 119 A:C mismatches", ylab="DHX9 eCLIP peak density")
# lines(smoothingSpline2, col="orange", lwd=2)
# lines(smoothingSpline3, col="dark blue", lwd=2)

DHX9_eCLIP_XL1_3UTR <- df2
DHX9_eCLIP_XL2_3UTR <- df3

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

DHX9_eCLIP_XL1_intro <- df2
DHX9_eCLIP_XL2_intro <- df3


df_test = DHX9_eCLIP_XL1_3UTR[,4]
df_test = DHX9_eCLIP_XL1_3UTR[6001:14001,4]
DHX9_eCLIP_XL1_3UTR = data.frame(Peak_density=df_test)
DHX9_eCLIP_XL1_3UTR$type = "3UTR_sites.XL1"
df_test = DHX9_eCLIP_XL2_3UTR[6001:14001,4]
DHX9_eCLIP_XL2_3UTR = data.frame(Peak_density=df_test)
DHX9_eCLIP_XL2_3UTR$type = "3UTR_sites.XL2"


df_test = DHX9_eCLIP_XL1_intro[6001:14001,4]
DHX9_eCLIP_XL1_intro = data.frame(Peak_density=df_test)
DHX9_eCLIP_XL1_intro$type = "Intro_sites.XL1"
df_test = DHX9_eCLIP_XL2_intro[6001:14001,4]
DHX9_eCLIP_XL2_intro = data.frame(Peak_density=df_test)
DHX9_eCLIP_XL2_intro$type = "Intro_sites.XL2"


df_total = rbind(DHX9_eCLIP_XL1_3UTR, DHX9_eCLIP_XL1_intro)
my_comparisons <- list(c("Intro_sites.XL1", "3UTR_sites.XL1"))
p = ggboxplot(df_total, x = "type", y = "Peak_density",
              color = "type", palette = "jco")+
  # Add pairwise comparisons p-value
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")   
p + rremove("x.text") + rremove("xlab") + rremove("x.ticks")

df_total = rbind(DHX9_eCLIP_XL2_3UTR, DHX9_eCLIP_XL2_intro)
my_comparisons <- list(c("Intro_sites.XL2", "3UTR_sites.XL2"))
p = ggboxplot(df_total, x = "type", y = "Peak_density",
              color = "type", palette = "jco")+
  # Add pairwise comparisons p-value
  stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
p + rremove("x.text") + rremove("xlab") + rremove("x.ticks")




