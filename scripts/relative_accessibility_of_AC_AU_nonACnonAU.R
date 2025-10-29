## plot boxplot for accessibility score for Larry (1095 DDX6-repressed editing sites and non-DDX6-repressed regulated sites)
setwd("C:/Users/CSIV149/PLC.lab.docs/03.larry/relativeAccessibility")
df=read.table("RelativeAccessibilityForA7merEndingAtPosition_for_Larry.1095_DDX6_repressed_sites_vs_1095_random_non_DDX6_repressed_sites.txt", stringsAsFactors = F, check.names = F, header = F, fill=T)
colnames(df) = c("a", "condition")
#df$condition <- as.factor(df$condition)
#df$condition = factor(df$condition, level = c("downregulated_editing_sites", "upregulated_editing_sites", "unchanged_editing_sites", "randomA"))
#df1 = subset(df, df$condition %in% "upregulated_editing_sites")
#df2 = subset(df, df$condition %in% "randomA")
#df = rbind(df1, df2)
#df = df %>% head(222)
df$condition = factor(df$condition, level = c("DDX6-repressed_editing_sites", "non-DDX6-repressed_editing_sites"))
df$b = -df$a


#df = arrange(df, desc(condition))
my_comparisons <- list(c("DDX6-repressed_editing_sites", "non-DDX6-repressed_editing_sites"))
ggboxplot(df, "condition", "b", ylab = "Relative accessibility (reference A vs edited G, log2)",
          color = "condition", palette = c("#00AFBB", "#FC4E07"),
          add = "jitter", shape = "condition") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  rremove("x.text") + rremove("xlab") + rremove("x.ticks")



## plot boxplot for accessibility score for Larry (A:C vs A:U vs non A:C non A:U)
df=read.table("1095_DDX6_differential_editing_sites.A2C_A2U_nonC_nonU.output.overlapped_with_DDX6_peaks_XL1_XL2.txt", stringsAsFactors = F, check.names = F, header = F, fill=T)
colnames(df) = c("a", "condition")
#df$condition <- as.factor(df$condition)
#df$condition = factor(df$condition, level = c("downregulated_editing_sites", "upregulated_editing_sites", "unchanged_editing_sites", "randomA"))
#df1 = subset(df, df$condition %in% "upregulated_editing_sites")
#df2 = subset(df, df$condition %in% "randomA")
#df = rbind(df1, df2)
#df = df %>% head(222)
df$condition = factor(df$condition, level = c("A_C", "A_U", "nonC_nonU"))
df$b = -df$a

library(ggsignif)
library(beeswarm)

##option2:
library(ggpubr)
my_comparisons <- list(c("A_C", "A_U"), c("A_C", "nonC_nonU"))
ggboxplot(df, "condition", "b", ylab = "Relative accessibility (reference A vs edited G, log2)",
          color = "condition", palette = "jco",
          add = "jitter", shape = "condition") + 
  stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  rremove("x.text") + rremove("xlab") + rremove("x.ticks")
  #stat_compare_means(comparisons = my_comparisons, method = "wilcox.test")
graph2ppt(file="00.RelativeAccessibility.pptx", width=8, aspectr=1.2)

