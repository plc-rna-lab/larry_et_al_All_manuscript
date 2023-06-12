
library(Seurat)
library(cowplot)
############################################################################################
##colon cancer // colon cancer //
######################################################################################
setwd("C:/Users/CSIV149/PLC.lab.docs/03.larry/GSE140288_colonCancer_RAW/")
heto.data <- Read10X(data.dir = "./S1")
stim <- CreateSeuratObject(counts = heto.data, project = "heto3k", min.cells = 3, min.features = 200)

thymus.data <- Read10X(data.dir = "./S2")
ctrl <- CreateSeuratObject(counts = thymus.data, project = "thymus3k", min.cells = 3, min.features = 200)

# Set up stimulated object（colon cancer sample 1）
ctrl$stim <- "S1"
ctrl <- subset(ctrl, subset = nFeature_RNA > 500)
ctrl <- NormalizeData(ctrl, verbose = FALSE)
ctrl <- FindVariableFeatures(ctrl, selection.method = "vst", nfeatures = 2000)
# Set up stimulated object（colon cancer sample 2）
stim$stim <- "S2"
stim <- subset(stim, subset = nFeature_RNA > 500)
stim <- NormalizeData(stim, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000)


############################################################################################
##ESCC cancer // ESCC cancer //
######################################################################################
setwd("C:\\Users\\CSIV149\\PLC.lab.docs\\03.larry\\GSE188900_scESCC_RAW")
p1.data <- Read10X(data.dir = ".\\Patient1")
p1 <- CreateSeuratObject(counts = p1.data, project = "p1", min.cells = 3, min.features = 200)

p2.data <- Read10X(data.dir = ".\\Patient2")
p2 <- CreateSeuratObject(counts = p2.data, project = "p2", min.cells = 3, min.features = 200)

p3.data <- Read10X(data.dir = ".\\Patient3")
p3 <- CreateSeuratObject(counts = p3.data, project = "p3", min.cells = 3, min.features = 200)

p4.data <- Read10X(data.dir = ".\\Patient4")
p4 <- CreateSeuratObject(counts = p4.data, project = "p4", min.cells = 3, min.features = 200)

p5.data <- Read10X(data.dir = ".\\Patient5")
p5 <- CreateSeuratObject(counts = p5.data, project = "p5", min.cells = 3, min.features = 200)

p6.1.data <- Read10X(data.dir = ".\\Patient6.1")
p6.1 <- CreateSeuratObject(counts = p6.1.data, project = "p6.1", min.cells = 3, min.features = 200)

p6.2.data <- Read10X(data.dir = ".\\Patient6.2")
p6.2 <- CreateSeuratObject(counts = p6.2.data, project = "p6.2", min.cells = 3, min.features = 200)

p1[["percent.mt"]] <- PercentageFeatureSet(p1, pattern = "^MT-")
p2[["percent.mt"]] <- PercentageFeatureSet(p2, pattern = "^MT-")
p3[["percent.mt"]] <- PercentageFeatureSet(p3, pattern = "^MT-")
p4[["percent.mt"]] <- PercentageFeatureSet(p4, pattern = "^MT-")
p5[["percent.mt"]] <- PercentageFeatureSet(p5, pattern = "^MT-")
p6.1[["percent.mt"]] <- PercentageFeatureSet(p6.1, pattern = "^MT-")
p6.2[["percent.mt"]] <- PercentageFeatureSet(p6.2, pattern = "^MT-")
VlnPlot(p2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Set up stimulated objec
p1$stim <- "p1"
p1 <- subset(p1, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 10) %>% 
  NormalizeData(verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

p2$stim <- "p2"
p2 <- subset(p2, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 10) %>% 
  NormalizeData(verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

p3$stim <- "p3"
p3 <- subset(p3, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 10) %>% 
  NormalizeData(verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

p4$stim <- "p4"
p4 <- subset(p4, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 10) %>% 
  NormalizeData(verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

p5$stim <- "p5"
p5 <- subset(p5, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 10) %>% 
  NormalizeData(verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

p6.1$stim <- "p6.1"
p6.1 <- subset(p6.1, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 10) %>% 
  NormalizeData(verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

p6.2$stim <- "p6.2"
p6.2 <- subset(p6.2, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 10) %>% 
  NormalizeData(verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)


############################################################################################
##HCC cancer // HCC cancer //
######################################################################################
setwd("C:/Users/CSIV149/PLC.lab.docs/03.larry/GSE156625_HCC/")
heto.data <- Read10X(data.dir = "./hcc")
stim <- CreateSeuratObject(counts = heto.data, project = "heto3k", min.cells = 3, min.features = 200)
remove(heto.data)
# Set up stimulated object（colon cancer sample 2）
stim[["percent.mt"]] <- PercentageFeatureSet(stim, pattern = "^MT-")
VlnPlot(stim, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

stim <- subset(stim, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 25)
stim <- NormalizeData(stim, verbose = FALSE)
stim <- FindVariableFeatures(stim, selection.method = "vst", nfeatures = 2000)
immune.combined <- stim
remove(stim)

##################################################################################
#merge data
immune.anchors <- FindIntegrationAnchors(object.list = list(p1,p2,p3,p4,p5,p6.1,p6.2), dims = 1:10)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:10)

DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 20, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:10)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:10)
immune.combined <- FindClusters(immune.combined, resolution = 0.3)
immune.combined <- RunUMAP(immune.combined, dims = 1:20)
immune.combined <- RunTSNE(immune.combined, dims = 1:20)

saveRDS(immune.combined, file = "S1andS2.merge.rds")
saveRDS(immune.combined, file = "ESCC.patients.merge.rds")
saveRDS(immune.combined, file = "HCC.patients.merge.rds")

# Visualization
p1 <- DimPlot(immune.combined, reduction = "umap", group.by = "stim")
p2 <- DimPlot(immune.combined, reduction = "umap", label = TRUE)
plot_grid(p1, p2)


#corr plotting
immune.combined <- readRDS("ESCC.patients.merge.rds")
genes_isg = c("IRF7","IFIT1","ISG15","STAT1","IFITM3","IFI44","TGTP1","IFNB1")

immune.combined = ScaleData(immune.combined, assay = "integrated", features = genes_isg, do.center = T, do.scale = T)
immune.combined@meta.data$ISG_score = colMeans(immune.combined[["integrated"]]@scale.data)
#for HCC, it is differnet.
immune.combined = ScaleData(immune.combined, assay = "RNA", features = genes_isg, do.center = T, do.scale = T)
immune.combined@meta.data$ISG_score = colMeans(immune.combined[["RNA"]]@scale.data)

FeaturePlot(immune.combined, features = "ISG_score", label = F, 
            reduction = "umap")

isg.dat <- immune.combined[['ISG_score']]
head(isg.dat)
colnames(isg.dat)[1] <- "ISG"
class(isg.dat)

immune.combined <- readRDS("ESCC.patients.merge.rds")
FeaturePlot(immune.combined, features = "DDX6", label = F, reduction = "umap")
ddx6 = c("DDX6")
immune.combined = ScaleData(immune.combined, assay = "RNA", features = ddx6, do.center = T, do.scale = T)
immune.combined@meta.data$ddx6 = colMeans(immune.combined[["RNA"]]@scale.data)
FeaturePlot(immune.combined, features = "ddx6", label = F, 
            reduction = "umap")

ddx6.dat <- immune.combined[['ddx6']]
head(ddx6.dat)
colnames(ddx6.dat)[1] <- "DDX6"

dat <- cbind(isg.dat, ddx6.dat)
dat[1:5,1:2]
ggscatterstats(dat, 
               y =DDX6, 
               x =ISG,
               type = "pearson",
               centrality.para = "mean",                              
               margins = "both",                                         
               xfill = "#009E73",
               yfill = "#D55E00", 
               marginal.type = "histogram",  
               title = "Relationship between ISG and DDX6")


#get DDX6 high and DDX6 low group
#dat.median <- median(dat$DDX6)
#hig.dat <-subset(dat, dat$DDX6>dat.median)
#low.dat <-subset(dat, dat$DDX6<=dat.median)

#install.packages("ggstatsplot")
#library(ggstatsplot)
#ggscatterstats(hig.dat, 
               #y =DDX6, 
               #x =ISG,
               #type = "pearson",
               #centrality.para = "mean",                              
               #margins = "both",                                         
               #xfill = "#009E73",
               #yfill = "#D55E00", 
               #marginal.type = "histogram",  
               #title = "Relationship between ISG and high DDX6")


#ggscatterstats(low.dat, 
#               y =DDX6, 
#               x =ISG,
#               type = "pearson",
#               centrality.para = "mean",                              
#               margins = "both",                                         
#               xfill = "#009E73",
#               yfill = "#D55E00", 
#               marginal.type = "histogram",  
#               title = "Relationship between ISG and low DDX6")



##get DDX6 high and DDX6 low expression by using median of expression.
dat.median <- median(dat$DDX6)
dat.median <- mean(dat$DDX6)
##function
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
# Calculate the mode using the user function.
dat.median <- getmode(dat$DDX6)

hig.dat <-subset(dat, dat$DDX6>dat.median)
low.dat <-subset(dat, dat$DDX6<=dat.median)

hig.dat$group <- "DDX6_high"
low.dat$group <- "DDX6_low"

dat.new <- rbind(hig.dat, low.dat)
table(dat.new$group)

library(ggplot2)
library(ggsignif)
#install.packages('ggpubr')
library(ggpubr)

#significant test
shapiro.test(dat.new$ISG)
ggqqplot(dat.new$ISG)
ggdensity(dat.new$ISG)

wilcox.test(hig.dat$ISG, low.dat$ISG)

ggplot(dat.new, aes(x = dat.new$group, y = dat.new$ISG)) +
  geom_violin(trim = FALSE, aes(fill = group)) +
  geom_signif(comparisons = list(c("DDX6_high", "DDX6_low")),
              map_signif_level=T,
              textsize=4,test=wilcox.test,step_increase=2.0) + 
  #stat_compare_means(method = "anova") +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
               geom="pointrange", color = "black")+
  theme_classic()



#######################################################################
####################################################################
#
hig.dat <-subset(dat, dat$DDX6>dat.median)
dat.new.median <- median(hig.dat$DDX6)

hig.new.dat <- subset(hig.dat, hig.dat$DDX6>dat.new.median)
low.new.dat <- subset(hig.dat, hig.dat$DDX6<=dat.new.median)

hig.new.dat$group <- "DDX6_high"
low.new.dat$group <- "DDX6_low"

dat.new.new <- rbind(hig.new.dat, low.new.dat)
table(dat.new.new$group)

wilcox.test(dat.new.new$ISG, dat.new.new$ISG)

ggplot(dat.new.new, aes(x = dat.new.new$group, y = dat.new.new$ISG)) +
  geom_violin(trim = FALSE, aes(fill = group)) +
  geom_signif(comparisons = list(c("DDX6_high", "DDX6_low")),
              map_signif_level=T,
              textsize=4,test=wilcox.test,step_increase=2.0) + 
  #stat_compare_means(method = "anova") +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
               geom="pointrange", color = "black")+
  theme_classic()







############################################################################################
##ESCC cancer and normal // ESCC cancer and normal//
######################################################################################
library(Seurat)
library(cowplot)
library(ggplot2)
library(ggsignif)
#install.packages('ggpubr')
library(ggpubr)

setwd("C:\\Users\\CSIV149\\PLC.lab.docs\\03.larry\\GSE188900_scESCC_RAW")

p6.1.data <- Read10X(data.dir = ".\\Patient6.1")
p6.1 <- CreateSeuratObject(counts = p6.1.data, project = "p6.1", min.cells = 3, min.features = 200)

p6.2.data <- Read10X(data.dir = ".\\Patient6.2")
p6.2 <- CreateSeuratObject(counts = p6.2.data, project = "p6.2", min.cells = 3, min.features = 200)


p6.normal.data <- Read10X(data.dir = ".\\Patient6_normal")
p6.normal <- CreateSeuratObject(counts = p6.normal.data, project = "p6.normal", min.cells = 3, min.features = 200)

p6.1[["percent.mt"]] <- PercentageFeatureSet(p6.1, pattern = "^MT-")
p6.2[["percent.mt"]] <- PercentageFeatureSet(p6.2, pattern = "^MT-")
p6.normal[["percent.mt"]] <- PercentageFeatureSet(p6.normal, pattern = "^MT-")
VlnPlot(p6.normal, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# Set up stimulated objec
p6.1$stim <- "p6.1"
p6.1 <- subset(p6.1, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 10) %>% 
  NormalizeData(verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

p6.2$stim <- "p6.2"
p6.2 <- subset(p6.2, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 10) %>% 
  NormalizeData(verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

p6.normal$stim <- "p6.normal"
p6.normal <- subset(p6.normal, subset = nFeature_RNA > 100 & nFeature_RNA < 5000 & percent.mt < 10) %>% 
  NormalizeData(verbose = FALSE) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000)

##################################################################################

immune.anchors <- FindIntegrationAnchors(object.list = list(p6.1,p6.2,p6.normal), dims = 1:10)
immune.combined <- IntegrateData(anchorset = immune.anchors, dims = 1:10)

DefaultAssay(immune.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
immune.combined <- ScaleData(immune.combined, verbose = FALSE)
immune.combined <- RunPCA(immune.combined, npcs = 20, verbose = FALSE)
# t-SNE and Clustering
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:10)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca", dims = 1:10)
immune.combined <- FindClusters(immune.combined, resolution = 0.3)
immune.combined <- RunUMAP(immune.combined, dims = 1:20)
immune.combined <- RunTSNE(immune.combined, dims = 1:20)


saveRDS(immune.combined, file = "ESCC.patient6.cancerNormal.merged.rds")

# Visualization
DimPlot(immune.combined, reduction = "umap", group.by = "stim")

immune.combined <- readRDS("ESCC.patient6.cancerNormal.merged.rds")

exprs <- data.frame(FetchData(object = immune.combined, vars = c("stim", "DDX6")))
head(exprs)
table(exprs$stim)
exprs$stim[which(exprs$stim =='p6.1')] <- 'p6.2'

cancer.dat <-subset(exprs, exprs$stim=="p6.2")
normal.dat <-subset(exprs, exprs$stim=="p6.normal")

shapiro.test(dat.new$ISG)
ggqqplot(dat.new$ISG)
ggdensity(dat.new$ISG)

wilcox.test(cancer.dat$rna_DDX6, normal.dat$rna_DDX6)
exprs.new <- subset(exprs, exprs$rna_DDX6 < 1.0)

ggplot(exprs, aes(x = exprs$stim, y = exprs$rna_DDX6)) +
  geom_violin(aes(fill = stim)) +
  geom_signif(comparisons = list(c("p6.2", "p6.normal")),
              map_signif_level=T,
              textsize=4,test=wilcox.test,step_increase=2.0) + 
  stat_compare_means(method = "anova") +
  stat_summary(fun.data="mean_sdl", fun.args = list(mult=1),
               geom="pointrange", color = "black")+ 
  theme_classic()


