##########################################################
#########################################################
setwd("C:\\Users\\CSIV149\\PLC.lab.docs\\03.larry")
library(tidyverse)
library(DESeq2)
library(export)
#import data
mycounts<-read.table("Merged.gene.counts-All.RNA.helicases.KD.txt", header=TRUE, sep = "\t")
head(mycounts)
colnames(mycounts)

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
rows <- rownames(unique(mycounts['Gene']))
mycounts <- mycounts[rows,]

rownames(mycounts) <- mycounts[,1]
mycounts <- mycounts[,-1]
head(mycounts)

colData <- data.frame(c(colnames(mycounts)), c(rep("control",2),rep("treatment",2)))
colnames(colData) <- c("Sample", "condition")
head(colData)

dds <- DESeqDataSetFromMatrix(mycounts, colData, design= ~ condition)
dds <- DESeq(dds)
dds

res <- results(dds, contrast=c("condition", "treatment", "control"))##或耿
res <- res[order(res$pvalue),]
write.csv(res,file="EC109.DDX21_DESeq2_results.csv", quote = F)










