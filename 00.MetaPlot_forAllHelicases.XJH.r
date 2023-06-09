library(dplyr)
#install.packages("caroline")
library(caroline)
library(tidyr)
library(dplyr)
library(gplots)
library(ggplot2)
library(VennDiagram)
library(tidyverse)
library(export)

#DDX6的路径
setwd("C:\\Users\\CSIV149\\PLC.lab.docs\\03.larry\\eCLIP_analysis.RX.MetaPlot\\Larry_DDX6")
#另外三个样本的路径
setwd("C:\\Users\\CSIV149\\PLC.lab.docs\\03.larry\\eCLIP_analysis.RX.MetaPlot\\Larry_DDX21")
setwd("C:\\Users\\CSIV149\\PLC.lab.docs\\03.larry\\eCLIP_analysis.RX.MetaPlot\\Larry_DHX9")
setwd("C:\\Users\\CSIV149\\PLC.lab.docs\\03.larry\\eCLIP_analysis.RX.MetaPlot\\Larry_DDX3X")

## metagene plotting
bins_overlapping_peaks = read.delim("cds_bins_overlapped_eCLIP_peaks_DDX3X_XL2.txt", 
                                    stringsAsFactors = F, header = F, check.names = F, sep = "\t")
#count(bins_overlapping_peaks,V4) %>% plot(type="l", col="green", lwd=5, xlab="time", ylab="concentration", main="Exponential decay")

bins_overlapping_peaks_5utr = read.delim("5utr_bins_overlapped_eCLIP_peaks_DDX3X_XL2.txt", 
                                         stringsAsFactors = F, header = F, check.names = F, sep = "\t")
#count(bins_overlapping_peaks_5utr,V4) %>% plot(type="l", col="green", lwd=5, xlab="time", ylab="concentration", main="Exponential decay")

bins_overlapping_peaks_3utr = read.delim("3utr_bins_overlapped_eCLIP_peaks_DDX3X_XL2.txt", 
                                         stringsAsFactors = F, header = F, check.names = F, sep = "\t")
#count(bins_overlapping_peaks_3utr,V4) %>% plot(type="l", col="green", lwd=5, xlab="time", ylab="concentration", main="Exponential decay")
class(bins_overlapping_peaks_5utr)

total = rbind(dplyr::count(bins_overlapping_peaks_5utr,V4), 
              dplyr::count(bins_overlapping_peaks,V4), 
              dplyr::count(bins_overlapping_peaks_3utr,V4)) %>% select(-V4)
total = total %>% mutate(num = 1:n()) %>% select(num, n) #%>% plot(type="l", col="green", lwd=5, ylab="Peak Coverage")
total$n.ratio <- (total$n)/max(total$n)

#igg: "black", xl1: "blue", xl2:"red"
pt = ggplot(data=total, aes(x=num, y=n.ratio, group=1)) +
  geom_line(color="black", size = 2)+
  #geom_line(color="#E69F00", size = 2)+
  ylim(0,1) + #ddx6:ylim(0,7500),ddx21:ylim(0,1000)
  theme_minimal() +
  theme(axis.text.x = element_text(color = "black", size = 18, face = "plain"),
        axis.text.y = element_text(color = "black", size = 18, face = "plain"),  
        axis.title.x = element_text(color = "black", size = 12, face = "plain"),
        axis.title.y = element_text(color = "black", size = 12, face = "plain"))
pt + 
  geom_vline(xintercept = 12, linetype="dashed", color="grey18", size=1.2) + 
  geom_vline(xintercept = 113, linetype="dashed", color="grey18", size=1.2)

graph2ppt(file="MetaPlot_eCLIP_peaks_DDX6_IgG.pptx", width=8, aspectr=1.5)
graph2ppt(file="MetaPlot_eCLIP_peaks_DDX6_XL1.pptx", width=8, aspectr=1.5)
graph2ppt(file="MetaPlot_eCLIP_peaks_DDX6_XL2.pptx", width=8, aspectr=1.5)

graph2ppt(file="MetaPlot_eCLIP_peaks_DDX21_IgG.pptx", width=8, aspectr=1.5)
graph2ppt(file="MetaPlot_eCLIP_peaks_DDX21_XL1.pptx", width=8, aspectr=1.5)
graph2ppt(file="MetaPlot_eCLIP_peaks_DDX21_XL2.pptx", width=8, aspectr=1.5)

graph2ppt(file="MetaPlot_eCLIP_peaks_DHX9_IgG.pptx", width=8, aspectr=1.5)
graph2ppt(file="MetaPlot_eCLIP_peaks_DHX9_XL1.pptx", width=8, aspectr=1.5)
graph2ppt(file="MetaPlot_eCLIP_peaks_DHX9_XL2.pptx", width=8, aspectr=1.5)

graph2ppt(file="MetaPlot_eCLIP_peaks_DDX3X_IgG.pptx", width=8, aspectr=1.5)
graph2ppt(file="MetaPlot_eCLIP_peaks_DDX3X_XL1.pptx", width=8, aspectr=1.5)
graph2ppt(file="MetaPlot_eCLIP_peaks_DDX3X_XL2.pptx", width=8, aspectr=1.5)






