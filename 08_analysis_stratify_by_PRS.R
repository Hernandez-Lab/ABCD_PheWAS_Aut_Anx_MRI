#####
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(GGally)

#####
load_folder = "G:/My Drive/UCLA/HernandezLab/ABCD_Project_Scripts/data_frames/"
load_file = "abcd_r4_AutAnxMri_PhewasVars_4.Rda"
load(paste0(load_folder, load_file))

hist(merge_data_filtered$AUT_PRS)
hist(merge_data_filtered$ANX_PRS)
cor(merge_data_filtered$AUT_PRS, merge_data_filtered$ANX_PRS, use = "pairwise.complete.obs")
ggpaired(merge_data_filtered[51:100,], cond1 = "AUT_PRS", cond2 = "ANX_PRS", fill = "condition",
        line.color = "gray", line.size = 0.4, palette = "npg", show.legend = F) +
        stat_compare_means(paired = TRUE)
aut_prs_stdev = sd(merge_data_filtered$AUT_PRS)
anx_prs_stdev = sd(merge_data_filtered$ANX_PRS)
aut_prs_quantiles <- quantile(merge_data_filtered$AUT_PRS, probs=seq(0, 1, 0.1), na.rm = TRUE)
anx_prs_quantiles <- quantile(merge_data_filtered$ANX_PRS, probs=seq(0, 1, 0.1), na.rm = TRUE)

merge_data_filtered$AUT_PRS_STRAT = "med"
merge_data_filtered$AUT_PRS_STRAT[merge_data_filtered$AUT_PRS > aut_prs_stdev] = "hi"
merge_data_filtered$AUT_PRS_STRAT[merge_data_filtered$AUT_PRS > 2*aut_prs_stdev] = "ext_hi"
merge_data_filtered$AUT_PRS_STRAT[merge_data_filtered$AUT_PRS < -1*aut_prs_stdev] = "lo"
merge_data_filtered$AUT_PRS_STRAT[merge_data_filtered$AUT_PRS < -2*aut_prs_stdev] = "ext_lo"
merge_data_filtered$AUT_PRS_STRAT = factor(merge_data_filtered$AUT_PRS_STRAT, levels=c('ext_lo', 'lo', 'med', 'hi', 'ext_hi'))

merge_data_filtered$ANX_PRS_STRAT = "med"
merge_data_filtered$ANX_PRS_STRAT[merge_data_filtered$ANX_PRS > anx_prs_stdev] = "hi"
merge_data_filtered$ANX_PRS_STRAT[merge_data_filtered$ANX_PRS > 2*anx_prs_stdev] = "ext_hi"
merge_data_filtered$ANX_PRS_STRAT[merge_data_filtered$ANX_PRS < -1*anx_prs_stdev] = "lo"
merge_data_filtered$ANX_PRS_STRAT[merge_data_filtered$ANX_PRS < -2*anx_prs_stdev] = "ext_lo"
merge_data_filtered$ANX_PRS_STRAT = factor(merge_data_filtered$ANX_PRS_STRAT, levels=c('ext_lo', 'lo', 'med', 'hi', 'ext_hi'))

merge_data_filtered$COMBINED_PRS_STRAT = "med"
merge_data_filtered$COMBINED_PRS_STRAT[(merge_data_filtered$ANX_PRS_STRAT == "hi" |
                                         merge_data_filtered$ANX_PRS_STRAT == "ext_hi") &
                                          (merge_data_filtered$AUT_PRS_STRAT == "hi" |
                                             merge_data_filtered$AUT_PRS_STRAT == "ext_hi")] = "hiAnx-hiAut"
merge_data_filtered$COMBINED_PRS_STRAT[(merge_data_filtered$ANX_PRS_STRAT == "hi" |
                                          merge_data_filtered$ANX_PRS_STRAT == "ext_hi") &
                                         (merge_data_filtered$AUT_PRS_STRAT == "lo" |
                                            merge_data_filtered$AUT_PRS_STRAT == "ext_lo")] = "hiAnx-loAut"
merge_data_filtered$COMBINED_PRS_STRAT[(merge_data_filtered$ANX_PRS_STRAT == "lo" |
                                          merge_data_filtered$ANX_PRS_STRAT == "ext_lo") &
                                         (merge_data_filtered$AUT_PRS_STRAT == "hi" |
                                            merge_data_filtered$AUT_PRS_STRAT == "ext_hi")] = "loAnx-hiAut"
merge_data_filtered$COMBINED_PRS_STRAT[(merge_data_filtered$ANX_PRS_STRAT == "lo" |
                                          merge_data_filtered$ANX_PRS_STRAT == "ext_lo") &
                                         (merge_data_filtered$AUT_PRS_STRAT == "lo" |
                                            merge_data_filtered$AUT_PRS_STRAT == "ext_lo")] = "loAnx-loAut"
merge_data_filtered$COMBINED_PRS_STRAT = factor(merge_data_filtered$COMBINED_PRS_STRAT, levels=c('loAnx-loAut', 'loAnx-hiAut', 'hiAnx-loAut', 'hiAnx-hiAut', 'med'))

ggplot(merge_data_filtered, aes(x=AUT_PRS_STRAT, y=AUT_PRS)) + geom_boxplot(fill="slateblue", alpha=0.2)
ggplot(merge_data_filtered, aes(x=ANX_PRS_STRAT, y=ANX_PRS)) + geom_boxplot(fill="slateblue", alpha=0.2)
ggplot(merge_data_filtered, aes(x=COMBINED_PRS_STRAT, y=ANX_PRS)) + geom_boxplot(fill="slateblue", alpha=0.2)
ggplot(merge_data_filtered, aes(x=COMBINED_PRS_STRAT, y=AUT_PRS)) + geom_boxplot(fill="slateblue", alpha=0.2)

ggplot(merge_data_filtered, aes(x = AUT_PRS_STRAT)) + geom_bar(position="dodge", fill="#004600") +
  geom_text(stat = 'count', aes(label = after_stat(count), vjust = 0))
ggplot(merge_data_filtered, aes(x = ANX_PRS_STRAT)) + geom_bar(position="dodge", fill="#004600") +
  geom_text(stat = 'count', aes(label = after_stat(count), vjust = 0))
ggplot(merge_data_filtered, aes(x = COMBINED_PRS_STRAT)) + geom_bar(position="dodge", fill="#004600") +
  geom_text(stat = 'count', aes(label = after_stat(count), vjust = 0))

ggplot(merge_data_filtered, aes(x=AUT_PRS_STRAT, y=dmri_vol_uncinate, color=ANX_PRS_STRAT)) +
  geom_boxplot(fill="slateblue", alpha=0.2)
ggplot(merge_data_filtered, aes(x=ANX_PRS_STRAT, y=dmri_vol_uncinate, color=AUT_PRS_STRAT)) +
  geom_boxplot(fill="slateblue", alpha=0.2)

ggplot(merge_data_filtered, aes(x=COMBINED_PRS_STRAT, y=dmri_vol_uncinate)) +
  geom_boxplot(fill="slateblue", alpha=0.2) +
  geom_hline(yintercept=median(merge_data_filtered$dmri_vol_uncinate), linetype="dashed", color = "red") +
  stat_compare_means(comparisons = list(c(1,5),c(2,5),c(3,5),c(4,5)))

ggplot(merge_data_filtered, aes(x=COMBINED_PRS_STRAT, y=smri_sulc_cdk_tmpole)) +
  geom_boxplot(fill="slateblue", alpha=0.2)  +
  geom_hline(yintercept=median(merge_data_filtered$smri_sulc_cdk_tmpole), linetype="dashed", color = "red") +
  stat_compare_means(comparisons = list(c(1,5),c(2,5),c(3,5),c(4,5)))

ggpairs(merge_data_filtered, columns = c("dmri_rsirndgm_cdk_tp", "dmri_md_uncinate", "dmri_vol_uncinate","smri_sulc_cdk_tmpole"), aes(colour=COMBINED_PRS_STRAT,alpha=0.2))
ggpairs(merge_data_filtered, columns = c("dmri_rsirndgm_cdk_tp", "dmri_md_uncinate", "dmri_vol_uncinate","smri_sulc_cdk_tmpole"), aes(colour=AUT_PRS_STRAT,alpha=0.2))
ggpairs(merge_data_filtered, columns = c("dmri_rsirndgm_cdk_tp", "dmri_md_uncinate", "dmri_vol_uncinate","smri_sulc_cdk_tmpole"), aes(colour=ANX_PRS_STRAT,alpha=0.2))
