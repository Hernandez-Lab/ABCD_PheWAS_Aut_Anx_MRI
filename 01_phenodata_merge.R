#####
library(dplyr)
library(tidyverse)
library(ggplot2)
library(corrplot)

##### Functions #####
# from Leanna, '01_Merge.Phenotypes'
find_na <-  function (x)
{
  y <- sapply(x, function(x)all(is.na(x)))
  if (any(y))
  {
    stop(paste("All NA in columns", paste(which(y), collapse=", ")))
  }
}

##### Load raw data #####
abcd_folder = 'Z:/u/project/lhernand/shared/GenomicDatasets/ABCD_Release_4/phenotypes/'
long_track <- read.delim(paste0(abcd_folder, 'abcd_lt01.txt'), header = TRUE, na.strings = c("", "NA"))
long_track <- long_track[-c(1),]

##### Merge and create data skeleton #####
merge_data <- select(long_track, subjectkey, eventname, sex, interview_date, interview_age, site_id_l)
cols <- c("subjectkey", "eventname", "sex", "site_id_l")
merge_data[cols] <- lapply(merge_data[cols], factor) # factor these columns
merge_data$interview_age <- as.numeric(merge_data$interview_age) # make these column numeric (age in months)

##### Add some ancestry data from genetics
ancestry_folder = "Z:/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_4/genotype/TOPMed_imputed/"
euro_ids <- read.delim(paste0(ancestry_folder, 'ABCDr4_EUR.6679_no.sexmismatch_IDs.txt'), header = FALSE, na.strings = c("", "NA"))
euro_ids$genotype_ancestry <- "EUR"
euro_ids <- euro_ids[,-c(1)]
colnames(euro_ids)[1] <- "subjectkey"
euro_ids[,1] = sub(".*?_","",euro_ids[,1])
euro_ids$genotype_ancestry <- as.factor(euro_ids$genotype_ancestry)
merge_data <- merge(merge_data, euro_ids, by = c("subjectkey"), all.x = TRUE)
merge_data$genotype_ancestry <- factor(merge_data$genotype_ancestry,
                                     exclude = NULL, levels = c("EUR", NA), labels = c("EUR", "NON_EUR"))

##### Add ancestry PCs to use as covariates in later regressions
pc_folder = "Z:/u/project/lhernand/shared/GenomicDatasets-processed/ABCD_Release_4/genotype/ancestry_PCs/bigsnpr_final_no_outliers/"
pcs <- read.delim(paste0(pc_folder, 'local_pcs_1.20_EUR.txt'), header = FALSE, na.strings = c("", "NA"))
pcs <- pcs[-c(1),-c(1)]
colnames(pcs)[1] <- "subjectkey"
pcs[,1] = sub(".*?_","",pcs[,1])
cols = c(2:11)
pcs <- select(pcs, subjectkey, cols)
pcs[cols] <- lapply(pcs[cols], as.numeric)
colnames(pcs)[cols] <- c("pc1", "pc2", "pc3", "pc4", "pc5", "pc6", "pc7", "pc8", "pc9", "pc10")
merge_data <- merge(merge_data, pcs, by = c("subjectkey"), all.x = TRUE)

##### Initial data culling #####
# (only want baseline measurements - no longitudinal)
merge_data <- merge_data[merge_data$eventname == 'baseline_year_1_arm_1',]
merge_data <- merge_data[merge_data$genotype_ancestry == 'EUR',]

##### Other Participant information (ethnicty, etc) #####
community <- read.delim(paste0(abcd_folder, 'acspsw03.txt'), header = TRUE, na.strings = c("", "NA"))
community <- community[-c(1),]
community <- select(community, subjectkey, eventname, race_ethnicity)
community$race_ethnicity <- as.factor(community$race_ethnicity)
merge_data <- merge(merge_data, community, by = c("subjectkey", "eventname"), all.x = TRUE)

##### Add CBCL t-scores (anxiety) #####
cbcl_tscores <- read.delim(paste0(abcd_folder, 'abcd_cbcls01.txt'), header = TRUE, na.strings = c("", "NA"))
cbcl_tscores <- cbcl_tscores[-c(1),]
cbcl_tscores <- select(cbcl_tscores, subjectkey, eventname,
                       cbcl_scr_syn_anxdep_r, cbcl_scr_dsm5_anxdisord_r, cbcl_scr_syn_social_r)
cols <- 3:ncol(cbcl_tscores)
cbcl_tscores[cols] <- lapply(cbcl_tscores[cols], as.numeric)
#sapply(cbcl_tscores, class) do we ever need this?
merge_data <- merge(merge_data, cbcl_tscores, by = c("subjectkey", "eventname"), all.x = TRUE)

##### Add KSADS DSM-5 (anxiety) data #####
ksads_parent_dsm5 = read.delim(paste0(abcd_folder, 'abcd_ksad01.txt'), header = TRUE, na.strings = c("", "NA"))
ksads_parent_dsm5 <- ksads_parent_dsm5[-c(1),]
ksads_parent_dsm5 <- select(ksads_parent_dsm5, subjectkey, eventname,
                            ksads_8_864_p,
                            ksads_8_863_p,
                            ksads_8_912_p,
                            ksads_8_911_p,
                            ksads_8_309_p,
                            ksads_8_310_p,
                            ksads_8_312_p,
                            ksads_8_311_p,
                            ksads_10_914_p,
                            ksads_10_913_p,
                            ksads_10_869_p,
                            ksads_10_870_p)
ksads_parent_dsm5[ksads_parent_dsm5 == 555] <- NA
ksads_parent_dsm5[ksads_parent_dsm5 == 888] <- NA
cols <- 3:ncol(ksads_parent_dsm5)
ksads_parent_dsm5[cols] <- lapply(ksads_parent_dsm5[cols], as.numeric)
merge_data <- merge(merge_data, ksads_parent_dsm5, by = c("subjectkey", "eventname"), all.x = TRUE)

##### Add KSADS raw (autism) data #####
#ksads_rawscores_parent_aut = read.delim(paste0(abcd_folder, 'autism_spectrum_dis_p01.txt'), header = TRUE, na.strings = c("", "NA"))
#ksads_rawscores_parent_aut <- ksads_rawscores_parent_aut[-c(1),]
#cols <- 10:(ncol(ksads_rawscores_parent_aut)-1)
#ksads_rawscores_parent_aut <- select(ksads_rawscores_parent_aut, subjectkey, eventname,
#                       all_of(cols))
#cols <- 17:ncol(ksads_rawscores_parent_aut)
#ksads_rawscores_parent_aut[cols] <- lapply(ksads_rawscores_parent_aut[cols], factor)
#merge_data <- merge(merge_data, ksads_rawscores_parent_aut, by = c("subjectkey", "eventname"), all.x = TRUE)

##### Add ABCD autism data #####
#abcd_screen = read.delim(paste0(abcd_folder, 'abcd_screen01.txt'), header = TRUE, na.strings = c("", "NA"))
#abcd_screen <- abcd_screen[-c(1),]
# event name is different here for some reason, cannot use to merge
#abcd_screen <- select(abcd_screen, subjectkey, scrn_asd)
#abcd_screen$scrn_asd <- as.numeric(abcd_screen$scrn_asd)
#merge_data <- merge(merge_data, abcd_screen, by = c("subjectkey"), all.x = TRUE)

##### Add SRS data (autism related) #####
# There is only follow-up data here (questionairre not given at baseline visit)
srs_short = read.delim(paste0(abcd_folder, 'abcd_pssrs01.txt'), header = TRUE, na.strings = c("", "NA"))
srs_short <- srs_short[-c(1),]
srs_short <- select(srs_short, subjectkey,
                    ssrs_6_p,
                    ssrs_15r_p,
                    ssrs_16_p,
                    ssrs_18_p,
                    ssrs_24_p,
                    ssrs_29_p,
                    ssrs_35_p,
                    ssrs_37_p,
                    ssrs_39_p,
                    ssrs_42_p,
                    ssrs_58_p)
cols <- 2:ncol(srs_short)
srs_short[cols] <- lapply(srs_short[cols], as.numeric)
# we can't merge via event name
merge_data <- merge(merge_data, srs_short, by = c("subjectkey"), all.x = TRUE)

#### Clean up some stuff ####
rm(srs_short)
rm(cbcl_tscores)
rm(community)
rm(ksads_parent_dsm5)
rm(long_track)
rm(pcs)
rm(euro_ids)

#### MRI device info ####
mri_info <- read.delim(paste0(abcd_folder, "abcd_mri01.txt"),  header = TRUE, na.strings = c("", "NA"))
mri_info <- mri_info[-c(1),]
mri_info <- select(mri_info, subjectkey, eventname, mri_info_deviceserialnumber, mri_info_magneticfieldstrength)
cols = 3:4
mri_info[cols] <- lapply(mri_info[cols], factor)
merge_data <- merge(merge_data, mri_info, by = c("subjectkey", "eventname"), all.x = TRUE)

#### MRI QC info ####
mri_qc <- read.delim(paste0(abcd_folder, "abcd_imgincl01.txt"),  header = TRUE, na.strings = c("", "NA"))
mri_qc <- mri_qc[-c(1),]
cols = 11:(ncol(mri_qc)-1)
mri_qc <- select(mri_qc, subjectkey, eventname, cols)
cols = 3:ncol(mri_qc)
mri_qc[cols] <- lapply(mri_qc[cols], as.numeric)
merge_data <- merge(merge_data, mri_qc, by = c("subjectkey", "eventname"), all.x = TRUE)

#### MRI (basic structural) temporal pole and amygdala info ####
# + MRI (basic structural) insula info
# + MRI (basic structural) medial orbitofrontal info
smri_all <- read.delim(paste0(abcd_folder, "abcd_smrip10201.txt"),  header = TRUE, na.strings = c("", "NA"))
smri_all <- smri_all[-c(1),]
smri_all <- select(smri_all, subjectkey, eventname, 
                    smri_thick_cdk_tmpolelh,
                    smri_thick_cdk_tmpolerh,
                    smri_sulc_cdk_tmpolelh,
                    smri_sulc_cdk_tmpolerh,
                    smri_area_cdk_tmpolelh,
                    smri_area_cdk_tmpolerh,
                    smri_vol_cdk_tmpolelh,
                    smri_vol_cdk_tmpolerh,
                    smri_vol_scs_amygdalalh,
                    smri_vol_scs_amygdalarh,
                    smri_thick_cdk_mobfrlh,
                    smri_thick_cdk_mobfrrh,
                    smri_sulc_cdk_mobfrlh,
                    smri_sulc_cdk_mobfrrh,
                    smri_area_cdk_mobfrlh,
                    smri_area_cdk_mobfrrh,
                    smri_vol_cdk_mobfrlh,
                    smri_vol_cdk_mobfrrh,
                    smri_vol_scs_allventricles,
                    smri_vol_scs_intracranialv,
                    smri_vol_cdk_total)
cols = 3:ncol(smri_all)
smri_all[cols] <- lapply(smri_all[cols], as.numeric)
merge_data <- merge(merge_data, smri_all, by = c("subjectkey", "eventname"), all.x = TRUE)

#### MRI (DTI) uncinate fasciculus info
dti_all <- read.delim(paste0(abcd_folder, "abcd_dmdtifp101.txt"),  header = TRUE, na.strings = c("", "NA"))
dti_all <- dti_all[-c(1),]
dti_all <- select(dti_all, subjectkey, eventname, 
                    dmdtifp1_11,
                    dmdtifp1_12,
                    dmdtifp1_53,
                    dmdtifp1_54,
                    dmdtifp1_95,
                    dmdtifp1_96,
                    dmdtifp1_137,
                    dmdtifp1_138,
                    dmdtifp1_179,
                    dmdtifp1_180)
cols = 3:ncol(dti_all)
dti_all[cols] <- lapply(dti_all[cols], as.numeric)
merge_data <- merge(merge_data, dti_all, by = c("subjectkey", "eventname"), all.x = TRUE)

#### MRI (RSI, directional diffusion) temporal pole, medial orbitofrontal, uncinate fasciculus ()
rsi_dirdiff <- read.delim(paste0(abcd_folder, "abcd_drsip201.txt"),  header = TRUE, na.strings = c("", "NA"))
rsi_dirdiff <- rsi_dirdiff[-c(1),]
rsi_dirdiff <- select(rsi_dirdiff, subjectkey, eventname, 
                dmri_rsirndwm_cdk_tplh,
                dmri_rsirndwm_cdk_tprh,
                dmri_rsirndgm_cdk_tplh,
                dmri_rsirndgm_cdk_tprh,
                dmri_rsirnd_fib_uncrh,
                dmri_rsirnd_fib_unclh,
                dmri_rsirndwm_cdk_moflh,
                dmri_rsirndwm_cdk_mofrh,
                dmri_rsirndgm_cdk_moflh,
                dmri_rsirndgm_cdk_mofrh)
cols = 3:ncol(rsi_dirdiff)
rsi_dirdiff[cols] <- lapply(rsi_dirdiff[cols], as.numeric)
merge_data <- merge(merge_data, rsi_dirdiff, by = c("subjectkey", "eventname"), all.x = TRUE)

#### MRI (RSI, total diffusion) temporal pole, medial orbitofrontal, uncinate fasciculus ()
rsi_totdiff <- read.delim(paste0(abcd_folder, "abcd_drsip301.txt"),  header = TRUE, na.strings = c("", "NA"))
rsi_totdiff <- rsi_totdiff[-c(1),]
rsi_totdiff <- select(rsi_totdiff, subjectkey, eventname, 
                      dmri_rsirntwm_cdk_tplh,
                      dmri_rsirntwm_cdk_tprh,
                      dmri_rsirntgm_cdk_tplh,
                      dmri_rsirntgm_cdk_tprh,
                      dmri_rsirnt_fib_uncrh,
                      dmri_rsirnt_fib_unclh,
                      dmri_rsirntwm_cdk_moflh,
                      dmri_rsirntwm_cdk_mofrh,
                      dmri_rsirntgm_cdk_moflh,
                      dmri_rsirntgm_cdk_mofrh)
cols = 3:ncol(rsi_totdiff)
rsi_totdiff[cols] <- lapply(rsi_totdiff[cols], as.numeric)
merge_data <- merge(merge_data, rsi_totdiff, by = c("subjectkey", "eventname"), all.x = TRUE)

#### Clean up some stuff ####
rm(dti_all)
rm(mri_info)
rm(mri_qc)
rm(smri_all)
rm(rsi_dirdiff)
rm(rsi_totdiff)

##### Final data culling #####
find_na(merge_data)
# manually remove if necessary
#merge_data_filtered <- merge_data[,-c(22)]
merge_data_filtered <- merge_data[merge_data$imgincl_t1w_include == 1,]
merge_data_filtered <- merge_data_filtered[merge_data_filtered$imgincl_t2w_include == 1,]
merge_data_filtered <- merge_data_filtered[merge_data_filtered$imgincl_dmri_include == 1,]
# mri qc culls ~1500 participants

##### Save #####
save_folder = "G:/My Drive/UCLA/HernandezLab/ABCD_Project_Scripts/"
save_file_name = "abcd_r4_AutAnxMri_PhewasVars_AncestryFiltered-EURonly_LongitudinalFiltered-BaselineOnly_MRIFiltered-GoodScansOnly.Rda"
save(merge_data_filtered, file = paste0(save_folder, save_file_name))
