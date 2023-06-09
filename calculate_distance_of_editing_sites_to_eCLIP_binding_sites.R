
###
# THIS SCRIPT IS TO CALCULATE THE DISTANCE OF EDITING SITES TO ECLIP PEAK CENTERS. 
# If you want to analyze the distance from eCILP peak centers to editing sites, check Elmo:/media/ext2/renxi/projects/Sze_Jing_projects
### 
suppressMessages(library(magrittr))
suppressMessages(library(tidyr))
suppressMessages(library(dplyr))

## calculate differential gene expression
df1 = read.delim(paste0("C:/Users/e0149673/Desktop/eCLIP analysis/eclip_peaks_filtered_for_new_ADAR1.txt"), stringsAsFactors = F, header = T, check.names = F)
df2 = read.delim(paste0("C:/Users/e0149673/Desktop/eCLIP analysis/ADAR1_ADAR2_differential_editing_sites.txt"), stringsAsFactors = F, header = T, check.names = F)
df1_XL1 = subset(df1, df1$sample %in% "XL1") 
df2_ADAR1_KD = subset(df2, df2$specificity %in% "ADAR1_KD")

#edit_sites = rownames(df0)
out_file = paste0("distance_of_editing_sites_to_peak_binding_sites.txt")
out_dir = paste0("C:/Users/e0149673/Desktop/eCLIP analysis/")
dir.create(out_dir, showWarnings = F)
setwd(out_dir)

cat(file = out_file,
    "editing_names",
    "closest_peak_position",
    "distances",
    "\n",
    sep = "\t")

mybiglist <- list()
for(i in 1:nrow(df1_XL1)) {
  #chr_peaks = df1_XL1[j, 2] 
  #chr_editing_sites = df2_ADAR1_KD[i, 1]
  peak_binding_sites = (as.numeric(df1_XL1[i,3]) + as.numeric(df1_XL1[i,4]))/2
  mybiglist[i] = peak_binding_sites
  #edit_values_high = edit_values_high[!is.na(edit_values_high)]
  #edit_values_low = edit_values_low[!is.na(edit_values_low)]
  
  #if (chr_peaks != chr_editing_sites) {
  #  distances = NA
  #}
  #else {
  #  distances = as.numeric(df2_ADAR1_KD[i,2]) - peak_binding_sites
  #}
}

invisible(lapply(1:nrow(df2_ADAR1_KD), function(j) {
  editing_names = paste0(df2_ADAR1_KD[j, 1], "_", df2_ADAR1_KD[j, 2])
  editing_sites = df2_ADAR1_KD[j, 2]

  mybiglist = mybiglist%>% as.integer()
  closest_peak_position = mybiglist[which(abs(editing_sites-mybiglist)==min(abs(editing_sites-mybiglist)))]
  distances = editing_sites - closest_peak_position
  cat(file = out_file,
        editing_names,
        closest_peak_position,
        distances,
        "\n",
        sep = "\t",
        append = T)

}))
matrix = read.delim(paste0("C:/Users/e0149673/Desktop/eCLIP analysis/distance_of_editing_sites_to_peak_binding_sites.txt"), stringsAsFactors = F, header = T)
a = as.data.frame(table(matrix$distances))
plot(a$Var1,a$Freq)

df_expression = read.delim(paste0("E:/data/correlate_cancer-specific_APA_events_with_differential_gene_expression/differential_gene_expression_by_CHOL.txt"), stringsAsFactors = F, header = T)
df_expression = df_expression %>% mutate(log = log2(fold_ratio_high_over_low)) %>% select(APA_site, log)
