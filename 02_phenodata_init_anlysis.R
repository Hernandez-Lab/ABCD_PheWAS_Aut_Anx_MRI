#####
library(dplyr)
library(tidyverse)
library(ggplot2)
library(corrplot)
library(ggpubr)
library(ggVennDiagram)
library(hexbin)
library(car)
library(moments)

##### Functions #####
my_basic_plotting <- function(a, b)
{
  histA = hist(a, plot = FALSE, breaks = 20)
  histB = hist(b, plot = FALSE, breaks = 20)
  #plot(histA, col = rgb(173,216,230,max = 255, alpha = 80))
  #plot(histB, col = rgb(255,192,203, max = 255, alpha = 80), add = TRUE)
  plot(a, b)
}

##### Variable Maps #####
# 1 = White; 2 = Black; 3 = Hispanic; 4 = Asian; 5 = Other

##### Load Data #####
load_folder = "G:/My Drive/UCLA/HernandezLab/ABCD_Project_Scripts/data_frames/"
load_file = "abcd_r4_AutAnxMri_PhewasVars_AncestryFiltered-EURonly_LongitudinalFiltered-BaselineOnly_MRIFiltered-GoodScansOnly.Rda"
load(paste0(load_folder, load_file))
#merge_data_filtered = merge_data_filtered_eurOnly

##### Extra processing if necessary #####
levels(merge_data_filtered$race_ethnicity) <- 
  c("White", "Black", "Hispanic", "Asian", "Other")
merge_data_filtered$subjectkey <- as.character(merge_data_filtered$subjectkey)

##### Begin saving PDF #####
#pdf(paste0(load_folder, "abcd_r4_AutAnxMri_PhewasVars_plots.pdf"), onefile = TRUE)

##### Basic visualizations #####
hist(merge_data_filtered$interview_age)
hist(merge_data_filtered$interview_age[merge_data_filtered$sex == 'M'])
hist(merge_data_filtered$interview_age[merge_data_filtered$sex == 'F'])
hist(merge_data_filtered$interview_age[merge_data_filtered$race_ethnicity == 'White'])
hist(merge_data_filtered$interview_age[merge_data_filtered$race_ethnicity == 'Black'])
hist(merge_data_filtered$interview_age[merge_data_filtered$race_ethnicity == 'Hispanic'])
hist(merge_data_filtered$interview_age[merge_data_filtered$race_ethnicity == 'Asian'])

barplot(table(merge_data_filtered$sex))
barplot(table(merge_data_filtered$site_id_l))
ggplot(merge_data_filtered, aes(x = genotype_ancestry)) +
  geom_bar(position="dodge") +
  geom_text(stat = 'count', aes(label = after_stat(count), vjust = 0))
ggplot(merge_data_filtered, aes(x = race_ethnicity, fill = sex)) +
  geom_bar(position="dodge") #+
  #geom_text(stat = 'count', aes(label = after_stat(count), vjust = 0))
ggplot(merge_data_filtered, aes(x = site_id_l, fill = race_ethnicity)) +
  geom_bar(position="dodge")
ggplot(merge_data_filtered, aes(x = site_id_l, fill = sex)) +
  geom_bar(position="dodge")
ggplot(merge_data_filtered, aes(x=as.factor(site_id_l), y=interview_age)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) +
  stat_compare_means()
ggplot(merge_data_filtered, aes(x = race_ethnicity, fill = genotype_ancestry)) +
  geom_bar(position="dodge")

##### Anxiety continuum #####
hist(merge_data_filtered$cbcl_scr_syn_anxdep_r)
hist(merge_data_filtered$cbcl_scr_dsm5_anxdisord_r)
hist(merge_data_filtered$cbcl_scr_syn_social_r)
ggplot(merge_data_filtered, aes(x=as.factor(interview_age), y=cbcl_scr_syn_anxdep_r)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) +
  stat_compare_means()
ggplot(merge_data_filtered, aes(x=as.factor(sex), y=cbcl_scr_syn_anxdep_r)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) +
  stat_compare_means(comparisons = list(c("M", "F")))
#ggplot(merge_data_filtered, aes(x=as.factor(race_ethnicity), y=cbcl_scr_syn_anxdep_r)) + 
#  geom_boxplot(fill="slateblue", alpha=0.2) +
#  stat_compare_means(comparisons = list(c("White", "Black"), c("White", "Hispanic"), c("White", "Asian"), c("White", "Other")))
ggplot(merge_data_filtered, aes(x=as.factor(genotype_ancestry), y=cbcl_scr_syn_anxdep_r)) + 
  geom_violin(trim=FALSE) +
  stat_compare_means()

#df <- data.frame(merge_data_filtered[c(31:60),c(1,19,20)]) %>% pivot_longer(cols = c("cbcl_scr_syn_anxdep_r", "cbcl_scr_dsm5_anxdisord_r"), names_to = "measure")
#ggplot(df, aes(x = measure, y = value, group = subjectkey, color = factor(subjectkey))) +
#  geom_point(show.legend = F) +
#  geom_line(show.legend = F) +
#  stat_compare_means(show.legend = F)
df <- data.frame(merge_data_filtered[c(1:60),c(19,20,21)])
ggpaired(df, cond1 = "cbcl_scr_syn_anxdep_r", cond2 = "cbcl_scr_dsm5_anxdisord_r", fill = "condition",
         line.color = "gray", line.size = 0.4, palette = "npg", show.legend = F) +
  stat_compare_means(paired = TRUE)

##### Anxiety/Autism scores (KSADS/SRS) #####
ksads_columns = c(22:44)
for(i in ksads_columns)
{
  m_colname = colnames(merge_data_filtered)[i]
  print(ggplot(merge_data_filtered, aes(x = as.factor(merge_data_filtered[,i]), fill = as.factor(merge_data_filtered[,i]))) +
    ggtitle(m_colname) +
    geom_bar(show.legend = F) +
    geom_text(stat = 'count', aes(label = after_stat(count))))
}

##### Create Autism aggregate score #####
aut_agg_score = rowSums(merge_data_filtered[,34:44])
merge_data_filtered$aut_agg_score = aut_agg_score
head(merge_data_filtered$aut_agg_score)
hist(merge_data_filtered$aut_agg_score)

##### Autism continuum #####
ggplot(merge_data_filtered, aes(x=as.factor(interview_age), y=aut_agg_score)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) +
  stat_compare_means()
ggplot(merge_data_filtered, aes(x=as.factor(sex), y=aut_agg_score)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) +
  stat_compare_means(comparisons = list(c("M", "F")))
hist(merge_data_filtered$aut_agg_score[merge_data_filtered$sex == 'M'])
hist(merge_data_filtered$aut_agg_score[merge_data_filtered$sex == 'F'])
#ggplot(merge_data_filtered, aes(x=as.factor(race_ethnicity), y=aut_agg_score)) + 
#  geom_boxplot(fill="slateblue", alpha=0.2) +
#  stat_compare_means(comparisons = list(c("White", "Black"), c("White", "Hispanic"), c("White", "Asian"), c("White", "Other")))
ggplot(merge_data_filtered, aes(x=as.factor(genotype_ancestry), y=aut_agg_score)) + 
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1) +
  #geom_jitter(shape=16, position=position_jitter(0.2)) +
  stat_compare_means()

##### Anxiety/Autism overlap #####
plot(merge_data_filtered$aut_agg_score, merge_data_filtered$cbcl_scr_syn_anxdep_r)
m_bin <- hexbin(merge_data_filtered$aut_agg_score, merge_data_filtered$cbcl_scr_syn_anxdep_r, xbins=30)
cor(merge_data_filtered$aut_agg_score, merge_data_filtered$cbcl_scr_syn_anxdep_r, use = "pairwise.complete.obs")
plot(m_bin)
plot(merge_data_filtered$aut_agg_score, merge_data_filtered$cbcl_scr_dsm5_anxdisord_r)
m_bin <- hexbin(merge_data_filtered$aut_agg_score, merge_data_filtered$cbcl_scr_dsm5_anxdisord_r, xbins=30)
cor(merge_data_filtered$aut_agg_score, merge_data_filtered$cbcl_scr_dsm5_anxdisord_r, use = "pairwise.complete.obs")
plot(m_bin)
plot(merge_data_filtered$aut_agg_score, merge_data_filtered$cbcl_scr_syn_social_r)
m_bin <- hexbin(merge_data_filtered$aut_agg_score, merge_data_filtered$cbcl_scr_syn_social_r, xbins=30)
cor(merge_data_filtered$aut_agg_score, merge_data_filtered$cbcl_scr_syn_social_r, use = "pairwise.complete.obs")
plot(m_bin)

##### Create RDF: Restricted directional fraction. Proportion of RND/RNT #####
#

##### MRI measures #####

# Visualizedistributions, skewness, and outliers
mri_columns = c(54:(ncol(merge_data_filtered)-1))
for(i in mri_columns)
{
  m_colname = colnames(merge_data_filtered)[i]
  m_skewness = skewness(merge_data_filtered[,i], na.rm = TRUE)
  #print(hist(merge_data_filtered[,i], main = m_colname))
  print(ggplot(merge_data_filtered, aes(x=merge_data_filtered[,i])) + 
          geom_histogram(aes(y=..density..), colour="black", fill="white") +
          geom_density(alpha=.2, fill="#FF6666") +
          ggtitle(m_colname))
  #        geom_text(y=0.1, x=105, label=toString(m_skewness)))
  #qqnorm(merge_data_filtered[,i], pch = 1, frame = FALSE)
  #qqline(merge_data_filtered[,i], col = "steelblue", lwd = 2)
  cat(i, toString(m_skewness))
  qqPlot(merge_data_filtered[,i], main = m_colname)
  #sum(is.na(merge_data_filtered[,i]))
}
# most measures pretty skewed

# Try IQR outlier detection
for(i in mri_columns)
{
  m_colname = colnames(merge_data_filtered)[i]
  low_bound <- quantile(merge_data_filtered[,i], probs=.05, na.rm = TRUE)
  high_bound <- quantile(merge_data_filtered[,i], probs=.95, na.rm = TRUE)
  bound = high_bound-low_bound
  m_out_log_ind = merge_data_filtered[,i] > high_bound + (bound*2.5) | merge_data_filtered[,i] < low_bound - (bound*2.5)
  m_out_ind = which(m_out_log_ind == TRUE)
  m_outliers = na.omit(merge_data_filtered[m_out_ind,i])
  m_outliers_subjectkey = toString(na.omit(merge_data_filtered[m_out_ind,1]))
  cat(i, m_colname, "\n", m_out_ind, "\n", m_outliers_subjectkey, "\n", m_outliers, "\n")
}

# Remove outlier subjects manually (indicies will change as you remove!)
# NDAR_INVK9RPEUC1
# NDAR_INVKX5G6YJV
# NDAR_INVNHVAX6N2 
# NDAR_INVJPH56X7G
# NDAR_INVFMHWAFKJ
# NDAR_INV6ABHCVYT
# NDAR_INV8F1DEP52
# NDAR_INVJCGH65N1 
#ind = which(merge_data_filtered$subjectkey == "NDAR_INVJCGH65N1")
#merge_data_filtered = merge_data_filtered[-c(ind),]

my_mri_corplot <- function(m_data)
{
  mri_subset_cormat = cor(m_data, use = "complete.obs")
  corrplot(mri_subset_cormat, method = 'color', order = 'alphabet')
  corrplot(mri_subset_cormat[1:25, 1:25], method = 'color', order = 'alphabet')
  #plot(m_data$smri_area_cdk_mobfrlh, m_data$smri_vol_cdk_mobfrlh)
  corrplot(mri_subset_cormat[26:35, 26:35], method = 'color', order = 'alphabet')
  #plot(m_data$dmdtifp1_137, m_data$dmdtifp1_11)
  corrplot(mri_subset_cormat[36:nrow(mri_subset_cormat), 36:nrow(mri_subset_cormat)], method = 'color', order = 'alphabet')
}

mri_subset = select(merge_data_filtered, 42:(ncol(merge_data_filtered)-2))
my_mri_corplot(mri_subset)
mri_subset = subset(merge_data_filtered, sex == 'M')
mri_subset = select(mri_subset, 42:(ncol(mri_subset)-2))
my_mri_corplot(mri_subset)
mri_subset = subset(merge_data_filtered, sex == 'F')
mri_subset = select(mri_subset, 42:(ncol(mri_subset)-2))
my_mri_corplot(mri_subset)
mri_subset = subset(merge_data_filtered, race_ethnicity == 'White')
mri_subset = select(mri_subset, 42:(ncol(mri_subset)-2))
my_mri_corplot(mri_subset)
mri_subset = subset(merge_data_filtered, race_ethnicity == 'Black')
mri_subset = select(mri_subset, 42:(ncol(mri_subset)-2))
my_mri_corplot(mri_subset)

ggplot(merge_data_filtered, aes(x=as.factor(site_id_l), y=smri_vol_scs_wholeb)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) +
  stat_compare_means()
ggplot(merge_data_filtered, aes(x=as.factor(site_id_l), y=dmdtifp1_11)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) +
  stat_compare_means()
ggplot(merge_data_filtered, aes(x=as.factor(site_id_l), y=dmdtifp1_137)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) +
  stat_compare_means()
ggplot(merge_data_filtered, aes(x=as.factor(site_id_l), y=dmri_rsirnd_fib_unclh)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) +
  stat_compare_means()
ggplot(merge_data_filtered, aes(x=as.factor(interview_age), y=smri_vol_scs_wholeb, color = sex)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) 
ggplot(merge_data_filtered, aes(x=as.factor(sex), y=smri_vol_scs_wholeb)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) +
  stat_compare_means()
ggplot(merge_data_filtered, aes(x=as.factor(sex), y=dmdtifp1_137)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) +
  stat_compare_means()
ggplot(merge_data_filtered, aes(x=as.factor(sex), y=dmri_rsirnd_fib_unclh)) + 
  geom_boxplot(fill="slateblue", alpha=0.2) +
  stat_compare_means()
# site 1 and site 19 have irregularities

##### Complete PDF #####
#dev.off()

##### Save #####
save_folder = "G:/My Drive/UCLA/HernandezLab/ABCD_Project_Scripts/"
save_file_name = "abcd_r4_AutAnxMri_PhewasVars_AncestryFiltered-EURonly_LongitudinalFiltered-BaselineOnly_MRIFiltered-GoodScansOnlyAndExtremesRemoved.Rda"
save(merge_data_filtered, file = paste0(save_folder, save_file_name))

##### Extra #####
hist(merge_data_filtered$smri_thick_cdk_mobfrlh)
hist(merge_data_filtered$smri_sulc_cdk_mobfrlh)
my_basic_plotting(merge_data_filtered$dmdtifp1_11, merge_data_filtered$dmdtifp1_12)
my_basic_plotting(merge_data_filtered$dmdtifp1_179, merge_data_filtered$dmdtifp1_180)
my_basic_plotting(merge_data_filtered$smri_vol_cdk_insulalh, merge_data_filtered$smri_vol_cdk_insularh)
plot(merge_data_filtered$dmdtifp1_11, merge_data_filtered$dmdtifp1_12)
plot(merge_data_filtered$dmdtifp1_11, merge_data_filtered$dmdtifp1_180)

cbcl_scr_syn_anxdep_r_cat <- NA 
cbcl_scr_syn_anxdep_r_cat[merge_data_filtered$cbcl_scr_syn_anxdep_r < 50] <- "A Very Low"
cbcl_scr_syn_anxdep_r_cat[merge_data_filtered$cbcl_scr_syn_anxdep_r >= 50 & merge_data_filtered$cbcl_scr_syn_anxdep_r < 62.5] <- "B Middle Low"
cbcl_scr_syn_anxdep_r_cat[merge_data_filtered$cbcl_scr_syn_anxdep_r >= 62.5 & merge_data_filtered$cbcl_scr_syn_anxdep_r < 75] <- "C Middle"
cbcl_scr_syn_anxdep_r_cat[merge_data_filtered$cbcl_scr_syn_anxdep_r >= 75 & merge_data_filtered$cbcl_scr_syn_anxdep_r < 87.5] <- "D Middle High"
cbcl_scr_syn_anxdep_r_cat[merge_data_filtered$cbcl_scr_syn_anxdep_r >= 87.5] <- "E High"
ggplot(merge_data_filtered, aes(x=cbcl_scr_syn_anxdep_r_cat, y=dmdtifp1_435)) + 
  geom_boxplot(fill="slateblue", alpha=0.2)
