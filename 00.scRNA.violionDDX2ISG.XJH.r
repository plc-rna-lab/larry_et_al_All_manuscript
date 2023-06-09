library(Seurat)
library(cowplot)

#绘制相关性图。
#####################################################################################
##colon cancer数据
#####################################################################################
setwd("C:/Users/CSIV149/PLC.lab.docs/03.larry/GSE140288_colonCancer_RAW/")
immune.combined <- readRDS("S1andS2.merge.rds")

#因为分不开，所以这个地方时用来分析报的res.
immune.combined <- FindClusters(immune.combined, resolution = 1.2) #res的范围是0.5到1.2
DimPlot(immune.combined, reduction = "umap",label = T)
DimPlot(immune.combined, reduction = "umap",label = T, 
        cells.highlight = WhichCells(immune.combined,idents = c("15")), 
        cols.highlight="black", sizes.highlight=2.0) 
#DimPlot(immune.combined, reduction = "tsne",label = T)

levels(immune.combined)

#根据两篇文章的信息，得到了每个亚群的分群信息
#内皮细胞包括stem, TA, Tuft, Enterocyte, 9_Goblet

#Stem: "LGR5"
FeaturePlot(immune.combined, features = c("LGR5")) #0,5,11 "stem"
#TA: "PCNA" , "MKI67" , "ATOH1" , "TMPRSS15"
FeaturePlot(immune.combined, features = c("PCNA")); #17,23,2,12
FeaturePlot(immune.combined, features = c("MKI67")); #18,25,26
FeaturePlot(immune.combined, features = c("ATOH1"));
FeaturePlot(immune.combined, features = c("TMPRSS15"))
#Tuft: "FOXI1" , "LRMP"
FeaturePlot(immune.combined, features = c("FOXI1")); 
FeaturePlot(immune.combined, features = c("LRMP")) #28
#Goblet: "ITLN1" , "ANO7" , "BCAS1" , "TFF1"
FeaturePlot(immune.combined, features = c("ITLN1")); #9，16
FeaturePlot(immune.combined, features = c("ANO7"));
FeaturePlot(immune.combined, features = c("BCAS1"));
FeaturePlot(immune.combined, features = c("TFF1"))
#Enterocyte: "FABP1" , "APOB" , "APOA1"
FeaturePlot(immune.combined, features = c("FABP1")); #27,22,19,3,1,6,24
FeaturePlot(immune.combined, features = c("APOB"));#4,8,14,20
FeaturePlot(immune.combined, features = c("APOA1")) #10,13

#其他的细胞：
#Paneth: "DEFA6"
FeaturePlot(immune.combined, features = c("DEFA6"))#21
#Endocrine: "CHGB"
FeaturePlot(immune.combined, features = c("CHGB"))# 分不出来，不要了。
#Mcell: "SPIB" , "BEST4"
FeaturePlot(immune.combined, features = c("SPIB"));
FeaturePlot(immune.combined, features = c("BEST4")) #15
#Imm: "PTPRC"
FeaturePlot(immune.combined, features = c("PTPRC"))# 7

new.cluster.id <- c("Stem","Enterocyte","TA","Enterocyte","Enterocyte", 
                    "Stem","Enterocyte","Immune","Enterocyte","Goblet","Enterocyte",
                    "Stem","TA","Enterocyte","Enterocyte","Mcell","Goblet", 
                    "TA","TA","Enterocyte","Enterocyte",
                    "Paneth","Enterocyte","TA","Enterocyte","TA", 
                    "TA","Enterocyte","Tuft")

names(new.cluster.id) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, new.cluster.id)
DimPlot(immune.combined, reduction = "umap",label = T, pt.size = 1.5)

colon.epthi <- immune.combined[,immune.combined@meta.data$seurat_clusters 
                               %in% c(0,5,11,17,23,2,12,18,25,26,28,
                                      9,16,27,22,19,3,1,6,24,4,8,14,20,10,13)]

genes_isg = c("IRF7","IFIT1","ISG15","STAT1","IFITM3","IFI44","TGTP1","IFNB1")
colon.epthi = ScaleData(colon.epthi, assay = "RNA", features = genes_isg, do.center = T, do.scale = T)
colon.epthi@meta.data$ISG_score = colMeans(colon.epthi[["RNA"]]@scale.data)
##在同一张图上同时展示两个基因的表达量。
FeaturePlot(colon.epthi, 
            features = c("ISG_score","DDX6"), 
            cols=c("white","darkorange3","blue"), pt.size = 3.0,
            blend=T,order = T) #order是让表达基因的细胞放在上面。

############################################################
##这个是用来计算ISG和DDX6表达量相关性的数据。
##############################################################
setwd("C:/Users/CSIV149/PLC.lab.docs/03.larry/GSE140288_colonCancer_RAW/")
immune.combined <- readRDS("S1andS2.merge.rds")
colon.epthi <- immune.combined[,immune.combined@meta.data$seurat_clusters %in% c(1,2,3,4,5,6,7,8)]
genes_isg = c("IRF7","IFIT1","ISG15","STAT1","IFITM3","IFI44","TGTP1","IFNB1")
#colon.epthi = ScaleData(colon.epthi, assay = "integrated", features = genes_isg, do.center = T, do.scale = T)
#colon.epthi@meta.data$ISG_score = colMeans(colon.epthi[["integrated"]]@scale.data)
#对于HCC样本，这里没有integrated，所需要使用RNA。
colon.epthi = ScaleData(colon.epthi, assay = "RNA", features = genes_isg, do.center = T, do.scale = T)
colon.epthi@meta.data$ISG_score = colMeans(colon.epthi[["RNA"]]@scale.data)

FeaturePlot(colon.epthi, features = "ISG_score", label = F, 
            reduction = "umap")

isg.dat <- colon.epthi[['ISG_score']]
head(isg.dat)
colnames(isg.dat)[1] <- "ISG"
class(isg.dat)

##DDX6的表达量
FeaturePlot(colon.epthi, features = "DDX6", label = F, reduction = "umap")
ddx6 = c("DDX6")
colon.epthi = ScaleData(colon.epthi, assay = "RNA", features = ddx6, do.center = T, do.scale = T)
colon.epthi@meta.data$ddx6 = colMeans(colon.epthi[["RNA"]]@scale.data)
FeaturePlot(colon.epthi, features = "ddx6", label = F, 
            reduction = "umap")

ddx6.dat <- colon.epthi[['ddx6']]
head(ddx6.dat)
colnames(ddx6.dat)[1] <- "DDX6"
 
dat <- cbind(isg.dat, ddx6.dat)
dat[1:5,1:2]

#查看表达量的分布
ggplot(dat,aes(x=ISG))+
  geom_histogram(position = "identity",fill="lightcoral",color="black",bins=500,alpha=0.5)

ggplot(dat,aes(x=DDX6))+
  geom_histogram(position = "identity",fill="lightcoral",color="black",bins=500,alpha=0.5)


#######################################################################
##找到数据当中的众数，返回众数的次数
mode.ddx6 <- as.numeric(names(table(dat$DDX6))[table(dat$DDX6) == max(table(dat$DDX6))])
max(table(dat$DDX6))

mode.isg <- as.numeric(names(table(dat$ISG))[table(dat$ISG) == max(table(dat$ISG))])
max(table(dat$ISG))

##将两者都低表达的细胞删掉
dat.median.ddx6 <- quantile(dat$DDX6,0.25); dat.median.isg <- quantile(dat$ISG,0.25)
dat.test1 <- subset(dat, dat$DDX6 < dat.median.ddx6+0.001)
dat.test2 <- subset(dat, dat$ISG < dat.median.isg+0.001)
#install.packages("VennDiagram")
library(VennDiagram)
venn.diagram(list(dat.test1=dat.test1,dat.test2=dat.test2),
             filename = "VennDiagram.tif")


##提取重叠的细胞，即不表达的细胞，然后用表达的细胞来绘图。
dat.noExpress <- merge(dat.test1, dat.test2, by = 'row.names', all = F)
head(dat.noExpress)
dat.noExpress.names <- c(dat.noExpress$Row.names)
class(dat.noExpress.names)
dat.hiExpre <- dat[-which(rownames(dat) %in% dat.noExpress.names),]

dat <- dat.hiExpre
#library(stringr)
#str_match_all(dat.test1,dat.test2)

##区分高表达和低表达的基因
#dat.median <- median(dat$DDX6)
#min(dat$DDX6)

dat.median.low <- quantile(dat$DDX6,0.05)#dat.median.isg <- quantile(dat$ISG,0.25)
dat.median.hig <- quantile(dat$DDX6,0.95)
class(dat.median)
#这里有个问题，就是均值和中位数都会导致ESCC的分组差别很大。
#根据查资料结果，以及数据分布结果发现ESCC的数据偏左，众数区分更加具有说服力。
#所以计算众数。变量不变了，后面的代码懒得修改。
#getmode <- function(v) {
uniqv <- unique(v)
uniqv[which.max(tabulate(match(v, uniqv)))]
}
# Calculate the mode using the user function.
#dat.median <- getmode(dat$DDX6)

#提取后添加组别信息。
low.dat <-subset(dat, dat$DDX6 <= dat.median.low)#小于0.25分位的是低表达的
max(low.dat$DDX6)
hig.dat <-subset(dat, dat$DDX6>=dat.median.hig)#大于0.75分位的属于高表达的基因
max(dat$DDX6)


hig.dat$group <- "DDX6_high"
low.dat$group <- "DDX6_low"

dat.new <- rbind(hig.dat, low.dat)
table(dat.new$group)
write.table(dat.new, file="GSE140288_colonCancer_ISGandDDX6violionPlot.txt", sep = "\t", quote = F)

library(ggplot2)
library(ggsignif)
#install.packages('ggpubr')
library(ggpubr)
library(export)
#因为数据在进行显著性检测之前，需要判断检测方法。
#检测数据是否符合正态分布
shapiro.test(dat.new$ISG)
ggqqplot(dat.new$ISG)
ggdensity(dat.new$ISG)
wilcox.test(hig.dat$ISG, low.dat$ISG)

##这里画图不可以使用log进行转换，因为表达量大部分是负数，log转换之后是空值。
table(dat.new$group)
ggplot(dat.new, aes(x = dat.new$group, y = dat.new$ISG)) +
  geom_violin(trim = FALSE, aes(fill = group)) +
  geom_signif(comparisons = list(c("DDX6_high", "DDX6_low")),
              map_signif_level=T,
              textsize=4,test=wilcox.test,step_increase=2.0) + 
  geom_boxplot(aes(x=dat.new$group, y=dat.new$ISG), width=0.1)+
  theme_classic()
  
#stat_compare_means(method = "anova") +
#stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
#             geom="pointrange", color = "black")+
  
graph2ppt(file="colon_DDX2ISG_corr_deletLowCells.pptx", width=8, aspectr=1.0)

#####################################################################################
##ESCC 数据
#####################################################################################
setwd("C:\\Users\\CSIV149\\PLC.lab.docs\\03.larry\\GSE188900_scESCC_RAW")
immune.combined <- readRDS("ESCC.patients.merge.rds")
immune.combined <- FindClusters(immune.combined, resolution = 1.0)
DimPlot(immune.combined, reduction = "umap",label = T)
DimPlot(immune.combined, reduction = "tsne",label = T)

levels(immune.combined)
#Squamous epithelial
FeaturePlot(immune.combined, features = c("KRT5"), order=T)#(2,6,10,11)
#B cells
FeaturePlot(immune.combined, features = c("CD79A"), order=T) #(8,21,26)
#Fibroblast
FeaturePlot(immune.combined, features = c("THY1"), order=T) #(9,27)
#Glandular epithelium
FeaturePlot(immune.combined, features = c("KRT7"), order=T) #(13,16,17,25)
#Lymphatic endothelium
FeaturePlot(immune.combined, features = c("PDPN"), order=T) #(4)
#Macrophages
FeaturePlot(immune.combined, features = c("CD14"), order=T) #（1,7,12,23）
#Mast cells
FeaturePlot(immune.combined, features = c("KIT"), order=T) #(19)
#Neutrophil
FeaturePlot(immune.combined, features = c("S100A8"), order=T) # 没有表达？
#Plasma cells
FeaturePlot(immune.combined, features = c("IGJ"), order=T) #(18)
#pDC
FeaturePlot(immune.combined, features = c("IL3RA"), order=T) #(14)
#Smooth muscle
FeaturePlot(immune.combined, features = c("ACTA2"), order=T) #(20,24)
#T cells
FeaturePlot(immune.combined, features = c("CD3E"), order=T) #(0,3,15,22)
#Endothelium
FeaturePlot(immune.combined, features = c("ICAM1"), order=T) #(5)
new.cluster.id <- c("T cells","Macrophages","Squamous epithelial","T cells","Lymphatic endothelium",
                    "Endothelium","Squamous epithelial","Macrophages","B cells","Fibroblast",
                    "Squamous epithelial","Squamous epithelial","Macrophages","Glandular epithelium","pDC",
                    "T cells","Glandular epithelium","Glandular epithelium","Plasma cells","Mast cells",
                    "Smooth muscle","B cells","T cells","Macrophages","Smooth muscle",
                    "Glandular epithelium","B cells","Fibroblast")

names(new.cluster.id) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, new.cluster.id)
DimPlot(immune.combined, reduction = "umap",label = T, pt.size = 1.5)

escc.epthi <- immune.combined[,immune.combined@meta.data$seurat_clusters 
                               %in% c(2,6,10,11)]


genes_isg = c("IRF7","IFIT1","ISG15","STAT1","IFITM3","IFI44","TGTP1","IFNB1")
escc.epthi = ScaleData(escc.epthi, assay = "RNA", features = genes_isg, do.center = T, do.scale = T)
escc.epthi@meta.data$ISG_score = colMeans(escc.epthi[["RNA"]]@scale.data)
isg.dat <- escc.epthi[['ISG_score']]
##在同一张图上同时展示两个基因的表达量。
FeaturePlot(escc.epthi, 
            features = c("ISG_score","DDX6"), 
            cols=c("white","darkorange3","blue"), pt.size = 3.0,
            blend=T,order = T) #order是让表达基因的细胞放在上面。


############################################################
##用于画图的
###################################################
setwd("C:\\Users\\CSIV149\\PLC.lab.docs\\03.larry\\GSE188900_scESCC_RAW")
immune.combined <- readRDS("ESCC.patients.merge.rds")
DimPlot(immune.combined, reduction = "umap",label = T)
escc.epthi <- immune.combined[,immune.combined@meta.data$seurat_clusters %in% c(2,6,9,10)]

genes_isg = c("IRF7","IFIT1","ISG15","STAT1","IFITM3","IFI44","TGTP1","IFNB1")
#escc.epthi = ScaleData(escc.epthi, assay = "integrated", features = genes_isg, do.center = T, do.scale = T)
#escc.epthi@meta.data$ISG_score = colMeans(escc.epthi[["integrated"]]@scale.data)
#对于HCC样本，这里没有integrated，所需要使用RNA。
escc.epthi = ScaleData(escc.epthi, assay = "RNA", features = genes_isg, do.center = T, do.scale = T)
escc.epthi@meta.data$ISG_score = colMeans(escc.epthi[["RNA"]]@scale.data)

FeaturePlot(escc.epthi, features = "ISG_score", label = F, 
            reduction = "umap")

isg.dat <- escc.epthi[['ISG_score']]
head(isg.dat)
colnames(isg.dat)[1] <- "ISG"
class(isg.dat)
##DDX6的表达量
escc.epthi <- immune.combined[,immune.combined@meta.data$seurat_clusters %in% c(2,6,9,10)]
FeaturePlot(escc.epthi, features = "DDX6", label = F, reduction = "umap")
ddx6 = c("DDX6")
escc.epthi = ScaleData(escc.epthi, assay = "RNA", features = ddx6, do.center = T, do.scale = T)
escc.epthi@meta.data$ddx6 = colMeans(escc.epthi[["RNA"]]@scale.data)
FeaturePlot(escc.epthi, features = "ddx6", label = F, 
            reduction = "umap")

ddx6.dat <- escc.epthi[['ddx6']]
head(ddx6.dat)
colnames(ddx6.dat)[1] <- "DDX6"

dat <- cbind(isg.dat, ddx6.dat)
dat[1:5,1:2]

#查看表达量的分布
ggplot(dat,aes(x=ISG))+
  geom_histogram(position = "identity",fill="lightcoral",color="black",bins=500,alpha=0.5)

ggplot(dat,aes(x=DDX6))+
  geom_histogram(position = "identity",fill="lightcoral",color="black",bins=500,alpha=0.5)
#######################################################################
##找到数据当中的众数，返回众数的次数
as.numeric(names(table(dat$DDX6))[table(dat$DDX6) == max(table(dat$DDX6))])
max(table(dat$DDX6))

as.numeric(names(table(dat$ISG))[table(dat$ISG) == max(table(dat$ISG))])
max(table(dat$ISG))

##将两者都低表达的细胞删掉
dat.median.ddx6 <- quantile(dat$DDX6,0.25); dat.median.isg <- quantile(dat$ISG,0.25)
dat.test1 <- subset(dat, dat$DDX6 < dat.median.ddx6+0.001)
dat.test2 <- subset(dat, dat$ISG < dat.median.isg+0.001)
#install.packages("VennDiagram")
library(VennDiagram)
venn.diagram(list(dat.test1=dat.test1,dat.test2=dat.test2),
             filename = "VennDiagram.tif")


##提取重叠的细胞，即不表达的细胞，然后用表达的细胞来绘图。
dat.noExpress <- merge(dat.test1, dat.test2, by = 'row.names', all = F)
head(dat.noExpress)
dat.noExpress.names <- c(dat.noExpress$Row.names)
class(dat.noExpress.names)
dat.hiExpre <- dat[-which(rownames(dat) %in% dat.noExpress.names),]

dat <- dat.hiExpre
#library(stringr)
#str_match_all(dat.test1,dat.test2)

##区分高表达和低表达的基因
#dat.median <- median(dat$DDX6)
#min(dat$DDX6)

dat.median.low <- quantile(dat$DDX6,0.05)#dat.median.isg <- quantile(dat$ISG,0.25)
dat.median.hig <- quantile(dat$DDX6,0.95)
class(dat.median)
#这里有个问题，就是均值和中位数都会导致ESCC的分组差别很大。
#根据查资料结果，以及数据分布结果发现ESCC的数据偏左，众数区分更加具有说服力。
#所以计算众数。变量不变了，后面的代码懒得修改。
#getmode <- function(v) {
uniqv <- unique(v)
uniqv[which.max(tabulate(match(v, uniqv)))]
}
# Calculate the mode using the user function.
#dat.median <- getmode(dat$DDX6)

#提取后添加组别信息。
low.dat <-subset(dat, dat$DDX6 <= dat.median.low)#小于0.25分位的是低表达的
max(low.dat$DDX6)
hig.dat <-subset(dat, dat$DDX6>=dat.median.hig)#大于0.75分位的属于高表达的基因
max(dat$DDX6)

hig.dat$group <- "DDX6_high"
low.dat$group <- "DDX6_low"

dat.new <- rbind(hig.dat, low.dat)
table(dat.new$group)
write.table(dat.new, file="GSE188900_scESCC_ISGandDDX6violionPlot.txt", sep = "\t", quote = F)


library(ggplot2)
library(ggsignif)
#install.packages('ggpubr')
library(ggpubr)
library(export)
#因为数据在进行显著性检测之前，需要判断检测方法。
#检测数据是否符合正态分布
shapiro.test(dat.new$ISG)
ggqqplot(dat.new$ISG)
ggdensity(dat.new$ISG)

wilcox.test(hig.dat$ISG, low.dat$ISG)

ggplot(dat.new, aes(x = dat.new$group, y = dat.new$ISG)) +
  geom_violin(trim = FALSE, aes(fill = group)) +
  geom_signif(comparisons = list(c("DDX6_high", "DDX6_low")),
              map_signif_level=T,
              textsize=4,test=wilcox.test,step_increase=2.0) + 
  geom_boxplot(aes(x=dat.new$group, y=dat.new$ISG), width=0.1)+
  theme_classic()

#stat_compare_means(method = "anova") +
#stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
#             geom="pointrange", color = "black")+

graph2ppt(file="escc_DDX2ISG_corr_deletLowCells.pptx", width=8, aspectr=1.0)



############################################################################################
##HCC // HCC//
######################################################################################
setwd("C:\\Users\\CSIV149\\PLC.lab.docs\\03.larry\\GSE156625_HCC")
immune.combined <- readRDS("HCC.patients.merge.rds")
immune.combined <- FindClusters(immune.combined, resolution = 0.8)
DimPlot(immune.combined, reduction = "umap",label = T)
DimPlot(immune.combined, reduction = "tsne",label = T)

#LYMPH 1 TC: IL7R, CD3E
FeaturePlot(immune.combined, features = c("IL7R"));
FeaturePlot(immune.combined, features = c("CD3E"))#(0,1,24)(CD4+ CELLS)
#LYMPH 2 CD4: CCR7, CD4
FeaturePlot(immune.combined, features = c("CCR7"));
FeaturePlot(immune.combined, features = c("CD4")) #
#LYMPH 3 TREG: TIGIT, BATF
FeaturePlot(immune.combined, features = c("TIGIT"));
FeaturePlot(immune.combined, features = c("BATF")) #
#LYMPH 4 CD8: CD8A, CD8B, CXCL13, MAPK13
FeaturePlot(immune.combined, features = c("CD8A"));
FeaturePlot(immune.combined, features = c("CD8B"));
FeaturePlot(immune.combined, features = c("CXCL13"));
FeaturePlot(immune.combined, features = c("MAPK13")) #(3,19,20,23)(CD8+ cells)
#LYMPH 5 CD8TRM: CRTAM
FeaturePlot(immune.combined, features = c("CRTAM")) #

#LYMPH 1 NK1: FGFBP2, MYOM2
FeaturePlot(immune.combined, features = c("FGFBP2"));
FeaturePlot(immune.combined, features = c("MYOM2")) #
#LYMPH 2 NK2: XCL1, XCL2, KLRF1
FeaturePlot(immune.combined, features = c("XCL1"));
FeaturePlot(immune.combined, features = c("KLRF1"));
FeaturePlot(immune.combined, features = c("XCL2")) #
#LYMPH 1 NKC: FGFBP2, KLRF1, GNLY
FeaturePlot(immune.combined, features = c("FGFBP2"));
FeaturePlot(immune.combined, features = c("KLRF1"));
FeaturePlot(immune.combined, features = c("GNLY")) #（17）(NK cells)

#MYE 1: CD163, CLEC10A, CD1C （Myeloid）
FeaturePlot(immune.combined, features = c("CD163"));
FeaturePlot(immune.combined, features = c("CLEC10A"));
FeaturePlot(immune.combined, features = c("CD1C")) #(5, 7, 25)

#Fibro: TAGLN, MYL9, ACTA2 (Fibroblasts)
FeaturePlot(immune.combined, features = c("TAGLN"));
FeaturePlot(immune.combined, features = c("MYL9"));
FeaturePlot(immune.combined, features = c("ACTA2")) #(14,18)

#ENDO: PLVAP, VWF (Endothelial cells)
FeaturePlot(immune.combined, features = c("PLVAP"));
FeaturePlot(immune.combined, features = c("VWF")) #(4,6,9,26,27)

#HEPA: ALB, APOC3, APOA2, FABP1, FGB (Hepatocytes)
FeaturePlot(immune.combined, features = c("ALB"));
FeaturePlot(immune.combined, features = c("APOC3"));
FeaturePlot(immune.combined, features = c("APOA2"));
FeaturePlot(immune.combined, features = c("FABP1"));
FeaturePlot(immune.combined, features = c("FGB")) #(15,8,16,22,10,12,21,2,11)
#BIOPOTENT: TNFRSF12A
FeaturePlot(immune.combined, features = c("TNFRSF12A"));
FeaturePlot(immune.combined, features = c("KRT19")) #()
#MAST: TPSAB1, CPA3 (Mast cells)
FeaturePlot(immune.combined, features = c("TPSAB1"));
FeaturePlot(immune.combined, features = c("CPA3")) #(13)

levels(immune.combined)
new.cluster.id <- c("CD4+ cells","CD4+ cells","Hepatocytes","CD8+ cells","Endothelial cells",
                    "Myeloid cells","Endothelial cells","Myeloid cells","Hepatocytes",
                    "Endothelial cells","Hepatocytes",
                    "Hepatocytes","Hepatocytes","Mast cells","Fibroblasts","Hepatocytes",
                    "Hepatocytes","NK cells","Fibroblasts","CD8+ cells","CD8+ cells",
                    "Hepatocytes","Hepatocytes","CD8+ cells","CD4+ cells","Myeloid cells",
                    "Endothelial cells","Endothelial cells")

names(new.cluster.id) <- levels(immune.combined)
immune.combined <- RenameIdents(immune.combined, new.cluster.id)
DimPlot(immune.combined, reduction = "umap",label = T, pt.size = 1.5)

#筛选出hepato细胞
hcc.epthi <- immune.combined[,immune.combined@meta.data$seurat_clusters 
                              %in% c(15,8,16,22,10,12,21,2,11)]
genes_isg = c("IRF7","IFIT1","ISG15","STAT1","IFITM3","IFI44","TGTP1","IFNB1")
hcc.epthi = ScaleData(hcc.epthi, assay = "RNA", features = genes_isg, do.center = T, do.scale = T)
hcc.epthi@meta.data$ISG_score = colMeans(hcc.epthi[["RNA"]]@scale.data)
isg.dat <- hcc.epthi[['ISG_score']]
##在同一张图上同时展示两个基因的表达量。
FeaturePlot(hcc.epthi, 
            features = c("ISG_score","DDX6"), 
            cols=c("white","darkorange3","blue"), pt.size = 3.0,
            blend=T,order = T) #order是让表达基因的细胞放在上面。

######################################################
##画相关性图的脚本。
####################################################
setwd("C:\\Users\\CSIV149\\PLC.lab.docs\\03.larry\\GSE156625_HCC")
immune.combined <- readRDS("HCC.patients.merge.rds")
DimPlot(immune.combined, reduction = "umap",label = T)

hcc.epthi <- immune.combined[,immune.combined@meta.data$seurat_clusters %in% c(2,5,6,9,10,11)]
genes_isg = c("IRF7","IFIT1","ISG15","STAT1","IFITM3","IFI44","TGTP1","IFNB1")
#对于HCC样本，这里没有integrated，所需要使用RNA。
hcc.epthi = ScaleData(hcc.epthi, assay = "RNA", features = genes_isg, do.center = T, do.scale = T)
hcc.epthi@meta.data$ISG_score = colMeans(hcc.epthi[["RNA"]]@scale.data)

FeaturePlot(hcc.epthi, features = "ISG_score", label = F, 
            reduction = "umap")

isg.dat <- hcc.epthi[['ISG_score']]
head(isg.dat)
colnames(isg.dat)[1] <- "ISG"
class(isg.dat)
##DDX6的表达量
hcc.epthi <- immune.combined[,immune.combined@meta.data$seurat_clusters %in% c(2,5,6,9,10,11)]
FeaturePlot(hcc.epthi, features = "DDX6", label = F, reduction = "umap")
ddx6 = c("DDX6")
hcc.epthi = ScaleData(hcc.epthi, assay = "RNA", features = ddx6, do.center = T, do.scale = T)
hcc.epthi@meta.data$ddx6 = colMeans(hcc.epthi[["RNA"]]@scale.data)
FeaturePlot(hcc.epthi, features = "ddx6", label = F, 
            reduction = "umap")

ddx6.dat <- hcc.epthi[['ddx6']]
head(ddx6.dat)
colnames(ddx6.dat)[1] <- "DDX6"

dat <- cbind(isg.dat, ddx6.dat)
head(dat)
#write.csv(dat,"hcc_dat.csv",quote = F)

#########################################################
##因为数据的分布呈现超级偏态分布，所以需要确定低表达这些细胞是否是DDX6和ISG同时低表达？
########################################################################################


#查看表达量的分布
ggplot(dat,aes(x=DDX6))+
  geom_histogram(position = "identity",fill="lightcoral",color="black",bins=500,alpha=0.5)

ggplot(dat,aes(x=ISG))+
  geom_histogram(position = "identity",fill="lightcoral",color="black",bins=500,alpha=0.5)

#dat$log.ISG <- log10(abs(dat$ISG))
#dat$log.DDX6 <- log10(abs(dat$DDX6))
#dat <- subset(dat, dat$log.ISG!="NaN" & dat$log.DDX6!="NaN" )
#dat[is.na(dat)] <- 0
##
##因为上述提取的数据里面存在负值，所以需要采用box-cox转换。
library(MASS)
#dat$number <- c(1:25792)
#head(dat)
#lm.sales1<-lm(DDX6~number,data=dat)
#res2<-lm.sales1$residual 
#hist(res2)

#b2=boxcox(DDX6~number, data=dat) # 定义函数类型和数据
#I=which(b2$y==max(b2$y))
#b2$x[I]

#######################################################################
##找到数据当中的众数，返回众数的次数
as.numeric(names(table(dat$DDX6))[table(dat$DDX6) == max(table(dat$DDX6))])
max(table(dat$DDX6))

as.numeric(names(table(dat$ISG))[table(dat$ISG) == max(table(dat$ISG))])
max(table(dat$ISG))

##将两者都低表达的细胞删掉
dat.median.ddx6 <- quantile(dat$DDX6,0.25); dat.median.isg <- quantile(dat$ISG,0.25)
dat.test1 <- subset(dat, dat$DDX6 < dat.median.ddx6+0.001)
dat.test2 <- subset(dat, dat$ISG < dat.median.isg+0.001)
#install.packages("VennDiagram")
library(VennDiagram)
venn.diagram(list(dat.test1=dat.test1,dat.test2=dat.test2),
             filename = "VennDiagram.tif")


##提取重叠的细胞，即不表达的细胞，然后用表达的细胞来绘图。
dat.noExpress <- merge(dat.test1, dat.test2, by = 'row.names', all = F)
head(dat.noExpress)
dat.noExpress.names <- c(dat.noExpress$Row.names)
class(dat.noExpress.names)
dat.hiExpre <- dat[-which(rownames(dat) %in% dat.noExpress.names),]

dat <- dat.hiExpre
#library(stringr)
#str_match_all(dat.test1,dat.test2)

##区分高表达和低表达的基因
dat.median <- median(dat$DDX6)
min(dat$DDX6)

dat.median.low <- quantile(dat$DDX6,0.05)#dat.median.isg <- quantile(dat$ISG,0.25)
dat.median.hig <- quantile(dat$DDX6,0.95)
class(dat.median)
#这里有个问题，就是均值和中位数都会导致ESCC的分组差别很大。
#根据查资料结果，以及数据分布结果发现ESCC的数据偏左，众数区分更加具有说服力。
#所以计算众数。变量不变了，后面的代码懒得修改。
#getmode <- function(v) {
uniqv <- unique(v)
uniqv[which.max(tabulate(match(v, uniqv)))]
}
# Calculate the mode using the user function.
#dat.median <- getmode(dat$DDX6)

#提取后添加组别信息。
low.dat <-subset(dat, dat$DDX6 <= dat.median.low)#小于0.25分位的是低表达的
max(low.dat$DDX6)
hig.dat <-subset(dat, dat$DDX6>=dat.median.hig)#大于0.75分位的属于高表达的基因
max(dat$DDX6)

hig.dat$group <- "DDX6_high"
low.dat$group <- "DDX6_low"

dat.new <- rbind(hig.dat, low.dat)
table(dat.new$group)
write.table(dat.new, file = "GSE156625_HCC_ISGandDDX6violionPlot.txt", sep = "\t", quote = F )

library(ggplot2)
library(ggsignif)
#install.packages('ggpubr')
library(ggpubr)
library(export)
#因为数据在进行显著性检测之前，需要判断检测方法。
#检测数据是否符合正态分布
shapiro.test(dat.new$ISG)
ggqqplot(dat.new$ISG)
ggdensity(dat.new$ISG)

wilcox.test(hig.dat$ISG, low.dat$ISG)
 
ggplot(dat.new, aes(x = dat.new$group, y = dat.new$ISG)) +
  geom_violin(trim = FALSE, aes(fill = group)) +
  geom_signif(comparisons = list(c("DDX6_high", "DDX6_low")),
              map_signif_level=T,
              textsize=4,test=wilcox.test,step_increase=2.0) + 
  geom_boxplot(aes(x=dat.new$group, y=dat.new$ISG), width=0.1)+
  theme_classic()

#stat_compare_means(method = "anova") +
#stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
#             geom="pointrange", color = "black")+

graph2ppt(file="hcc_DDX2ISG_corr_deletLowCells.pptx", width=8, aspectr=1.0)
a

