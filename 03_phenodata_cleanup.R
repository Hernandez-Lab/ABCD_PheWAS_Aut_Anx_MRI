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

##### Load Data #####
load_folder = "G:/My Drive/UCLA/HernandezLab/ABCD_Project_Scripts/data_frames/"
load_file = "abcd_r4_AutAnxMri_PhewasVars_AncestryFiltered-EURonly_LongitudinalFiltered-BaselineOnly_MRIFiltered-GoodScansOnlyAndExtremesRemoved.Rda"
load(paste0(load_folder, load_file))

##### variable name maps #####
dti_name_map = c(
  "dmri_fa_uncinaterh", # dmdtifp1_11
  "dmri_fa_uncinatelh", # dmdtifp1_12
  "dmri_md_uncinaterh", # dmdtifp1_53
  "dmri_md_uncinatelh", # dmdtifp1_54
  "dmri_ldc_uncinaterh", # dmdtifp1_95
  "dmri_ldc_uncinatelh", # dmdtifp1_96
  "dmri_tdc_uncinaterh", # dmdtifp1_137
  "dmri_tdc_uncinatelh", # dmdtifp1_138 
  "dmri_vol_uncinaterh", # dmdtifp1_179
  "dmri_vol_uncinatelh"  # dmdtifp1_180
)

ksads_name_map = c(
  "DiagSocAnxDisPast",		#ksads_8_864_p
  "DiagSocAnxDisPres",		#ksads_8_863_p
  "DiagSocAnxDisOtherPast", #ksads_8_912_p
  "DiagSocAnxDisOtherPres",		#ksads_8_911_p
  "SympImpairSocAnxDisPres",		#ksads_8_309_p
  "SympImpairSocAnxDisPast",		#ksads_8_310_p
  "SympClinSigDistSocAnxDisPast",		#ksads_8_312_p
  "SympClinSigDistSocAnxDisPres",		#ksads_8_311_p
  "DiagGenAnxDisOtherPast",	#ksads_10_914_p
  "DiagGenAnxDisOtherPres",		#ksads_10_913_p
  "DiagGenAnxDisPres",		#ksads_10_869_p
  "DiagGenAnxDisPast"		#ksads_10_870_p
)

##### Rename some columns #####
colnames(merge_data_filtered)[79:88] = dti_name_map
colnames(merge_data_filtered)[22:33] = ksads_name_map

##### 1 - Average MRI data across hemispheres #####
## 77 76 73 72
merge_data_filtered = merge_data_filtered[,-c(77)]
merge_data_filtered = merge_data_filtered[,-c(76)]
merge_data_filtered = merge_data_filtered[,-c(73)]
merge_data_filtered = merge_data_filtered[,-c(72)]
#mri_colnames = colnames(merge_data_filtered)[54:(ncol(merge_data_filtered)-1)]
mri_bihemisphere = c(seq(from=54, to=70, by=2), seq(from=75, to=103, by=2))
data_frame_mri_averages = data.frame(merge_data_filtered[,1])
for (i in mri_bihemisphere)
{
  colname = colnames(merge_data_filtered)[i]
  colname = str_sub(colname,1,-3)
  data_frame_mri_averages[colname] <- rowMeans(merge_data_filtered[,c(i,i+1)], na.rm=TRUE)
}
merge_data_filtered = merge_data_filtered[,-c(75:104)]
merge_data_filtered = merge_data_filtered[,-c(54:71)]
colnames(data_frame_mri_averages)[1] = "subjectkey"
merge_data_filtered = cbind(merge_data_filtered, data_frame_mri_averages[,2:ncol(data_frame_mri_averages)])

##### Save #####
save_folder = load_folder
save_file_name = "abcd_r4_AutAnxMri_PhewasVars_AncestryFiltered-EURonly_LongitudinalFiltered-BaselineOnly_MRIFiltered-GoodScansOnlyAndExtremesRemovedAndHemispheresAveraged.Rda"
save(merge_data_filtered, file = paste0(save_folder, save_file_name))
