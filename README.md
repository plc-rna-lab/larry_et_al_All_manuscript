# larry_et_al_All_manuscript

## Workflow

+ The Jupyter notebook for the full DDX6 eCLIP analysis workflow are provided, a html preview can be found at [DDX6_eCLIP_analysis.workflow.html](https://htmlpreview.github.io/?https://github.com/plc-rna-lab/larry_et_al_All_manuscript/blob/main/DDX6_eCLIP_analysis_workflow.html)

## Scripts used in the manuscript

+ `00.3UTRcompareINtronic.peakDensity_XL1_XL2.XJH.R`: Script for drawing peak density plot of 3'UTR and intron of editing sites.
+ `00.3UTRcompareINtronic.peakDensity_boxPlot_XL1_XL2.XJH.R`: Script for drawing box plot of 3'UTR and intron of editing sites.
+ `00.DDX6_peakDensity_of_regulated_and_unregulated_sites.XJH.R`: Script for drawing peak density plot of DDX6 regulated and un-regulated sites of DDX6.
+ `00.DEseq2.XJH.r`: Script for analyzing bulk RNA seq data.
+ `00.MetaPlot_forAllHelicases.XJH.r`: Script for analyzing metagene data.
+ `00.density_plot_XL1_XL2_500_-4000to4000.XJH.R`: Script for drawing peak density plot for DDX6, DDX3X and DDX21.
+ `00.mergeScData.XJH.R`: Script for analyzing single-cell data donwnloaded from GEO database.
+ `00.peakDistance_of_editingSites_eclip_4violions.XJH.R`: Script for calculating distance between editing sites and eclip peaks.
+ `00.scRNA.violionDDX2ISG.XJH.r`: Script for getting the violin-plot of single-cell data.
+ `run_overlapping_and_find_relative_positions.considering_the_ranges_of_whole_peaks.pl`: Perl script for generating eCLIP peak counts relative to editing sites.
+ `calculate_distance_of_editing_sites_to_eCLIP_binding_sites.R`: Script for calculting the distance between RNA editing sites and eCLIP peak.
+ `relative accessibility of AC AU and nonACnonAU.R`: Script for getting the relative accessibility plot of DDX6 editing sites.