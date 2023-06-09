##########################################################
##不同组之间找到差异基因。
#########################################################
setwd("C:\\Users\\CSIV149\\PLC.lab.docs\\03.larry")
library(tidyverse)
library(DESeq2)
library(export)

#首先：按照FPKM筛选表达的基因。
#然后：使用DEseq计算表达基因的差异。


#import data
mycounts<-read.table("Merged.gene.counts-All.RNA.helicases.KD.txt", header=TRUE, sep = "\t")
head(mycounts)
colnames(mycounts)
#提取相应的样本进行分析。
mycounts <- subset(mycounts, select=c(Gene, 
                                      EC109.DDX3X.Scr,
                                      EC109.DDX3X.Scr,
                                      EC109.DDX3X.sh1, 
                                      EC109.DDX3X.sh2))

mycounts <- subset(mycounts, select=c(Gene, 
                                      EC109.DDX6.Scr,
                                      EC109.DDX6.Scr,
                                      EC109.DDX6.sh6, 
                                      EC109.DDX6.sh9))

mycounts <- subset(mycounts, select=c(Gene, 
                                      EC109.DHX9.shScr,
                                      EC109.DHX9.shScr,
                                      EC109.shDHX9.5, 
                                      EC109.shDHX9.289))

mycounts <- subset(mycounts, select=c(Gene, 
                                      EC109.DDX21.shScr,
                                      EC109.DDX21.shScr,
                                      EC109.shDDX21.KD.A, 
                                      EC109.shDDX21.KD.B))
head(mycounts)
#因为处理成名字之后的数据会存在重复，所以需要吧重复的删掉〿
rows <- rownames(unique(mycounts['Gene']))
mycounts <- mycounts[rows,]

#这里有个x，需要去除，先把第一列当作行名来处理
rownames(mycounts) <- mycounts[,1]
mycounts <- mycounts[,-1]
head(mycounts)

#创建colData文件，确定对照组和处理组的信息。
#以con作为对照
colData <- data.frame(c(colnames(mycounts)), c(rep("control",2),rep("treatment",2)))
colnames(colData) <- c("Sample", "condition")
head(colData)

#构建数据矩阵
dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)
dds
#接下来，我们要查看treat VS control的总体结果，并根据p-value进行重新排序。
#利用summary命令统计显示一共多少个genes上调和下调（FDR0.1＿
res <- results(dds, contrast=c("condition", "treatment", "control"))##或耿
res <- res[order(res$pvalue),]
write.csv(res,file="EC109.DDX21_DESeq2_results.csv", quote = F)










