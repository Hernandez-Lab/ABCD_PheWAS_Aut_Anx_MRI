library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(lme4)
library(MarginalMediation)
library(GEEmediate)

##### TODO #####
# Test subspace? Whether effects of A and B are the same
# lm(Y ~ I(A + B) + C)
# lm(Y ~ A + B + C)

prs_folder_aut = 'Z:/u/project/lhernand/dcjaklic/ABCD_PRS/output_aut_grove_iPSYCH_2019/'
prs_folder_anx = 'Z:/u/project/lhernand/dcjaklic/ABCD_PRS/output_anx_meier_iPSYCH_2019/'
mri_folder = "G:/My Drive/UCLA/HernandezLab/ABCD_Project_Scripts/data_frames/"

prs_file_aut = "norm_prs_scores.txt"
prs_file_anx = "norm_prs_scores.txt"
mri_file = "abcd_r4_AutAnxMri_PhewasVars_4.Rda"

data_aut = read.table(paste0(prs_folder_aut, prs_file_aut))
data_aut = data_aut[-c(1),]
data_aut$V7 = as.numeric(data_aut$V7)

data_anx = read.table(paste0(prs_folder_anx, prs_file_anx))
data_anx = data_anx[-c(1),]
data_anx$V7 = as.numeric(data_anx$V7)

data_combined = data.frame(data_aut[,c(1,7)])
data_combined$V8 = data_anx[,c(7)]
colnames(data_combined) = c("subjectkey", "AUT_PRS", "ANX_PRS")
data_combined$subjectkey = sub(".*?_","",data_combined$subjectkey)

load(paste0(mri_folder, mri_file))

data_intersect <- intersect(merge_data_filtered$subjectkey, data_combined$subjectkey)
mri_data_intersect = merge_data_filtered[merge_data_filtered$subjectkey %in% data_intersect,]
prs_data_intersect = data_combined[data_combined$subjectkey %in% data_intersect,]

m_data = merge(mri_data_intersect, prs_data_intersect, by = c("subjectkey"), all.x = TRUE)

# Get variable name -> display name mapping
measure_info_path = "G:/My Drive/UCLA/HernandezLab/ABCD_Project_Scripts/"
measure_info_sheet = "abcd_behavior_mri_measures.xlsx"
measures_behavior = readxl::read_excel(paste0(measure_info_path, measure_info_sheet), sheet = "Behavioral")
measures_mri = readxl::read_excel(paste0(measure_info_path, measure_info_sheet), sheet = "MRI_ROI")
m_data_colnames = colnames(m_data)
variable_display_map = measures_behavior[,c("Variable name", "Display name")]
variable_display_map = rbind(variable_display_map, measures_mri[,c("Variable name", "Display name")])
a <- data.frame("aut_agg_score","Social Responsiveness Score")
names(a) <- c("Variable name","Display name")
variable_display_map <- rbind(variable_display_map, a)
variable_display_map <- na.omit(variable_display_map)
rownames(variable_display_map) = variable_display_map$`Variable name`

#prs_selector = c("AUT_PRS", "ANX_PRS")
indices_prs = 82:83
indices_behave_continuum = c(19:21,57) # glm gamma log, f test
indices_ksads_binary = 22:33 # glm binomial, chisq test
indices_mri_normaldist = 58:81 # lmer, normal anova
indices_smri = 58:66
indices_dti = 67:71
indices_rsi = 72:81

# Regression analysis
m_df <- data.frame()
for (this_index_prs in indices_prs)
{
  for (this_index_behave_measure in indices_behave_continuum)
  {
    m_model_1a <- glm(m_data[[this_index_behave_measure]]+1 ~ m_data[[this_index_prs]] + sex + site_id_l + interview_age +
                      pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                      data = m_data, na.action = na.omit, family = Gamma(link = "log"))
    m_model_1b <- glm(m_data[[this_index_behave_measure]]+1 ~ sex + site_id_l + interview_age +
                      pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                      data = m_data, na.action = na.omit, family = Gamma(link = "log"))
    # test sex interaction
    m_model_1c <- glm(m_data[[this_index_behave_measure]]+1 ~ m_data[[this_index_prs]] * sex + site_id_l + interview_age +
                        pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                      data = m_data, na.action = na.omit, family = Gamma(link = "log"))
    m_p1 <- anova(m_model_1b, m_model_1a, test = "F")
    m_p1_sex <- anova(m_model_1a, m_model_1c, test = "F")
    interaction_row1 = nrow(summary(m_model_1c)$coefficients)
    prs_string_ = colnames(m_data)[this_index_prs]
    measure_variable_name = colnames(m_data)[this_index_behave_measure]
    measure_display_name_ = as.character(variable_display_map[measure_variable_name,2])
    beta_ = summary(m_model_1a)$coefficients[2, 1]
    p_ = m_p1[2, 6]
    beta_sex_ = summary(m_model_1c)$coefficients[interaction_row1, 1]
    p_sex_ = m_p1_sex[2, 6]
    cat(prs_string_, " ", measure_display_name_, " ",
        beta_, " ", p_, " ",
        beta_sex_, " ", p_sex_, "\n")
    m_df = rbind(m_df, data.frame(PRS=prs_string_, Measure=measure_display_name_,
        beta=beta_, p=p_, p_sex=p_sex_))
  }
 for (this_index_ksads_measure in indices_ksads_binary)
 {
   m_model_2a <- glm(m_data[[this_index_ksads_measure]] ~ m_data[[this_index_prs]] + sex + site_id_l + interview_age +
                       pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                     data = m_data, na.action = na.omit, family = "binomial")
   m_model_2b <- glm(m_data[[this_index_ksads_measure]] ~ sex + site_id_l + interview_age +
                       pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                     data = m_data, na.action = na.omit, family = "binomial")
   # test sex interaction
   m_model_2c <- glm(m_data[[this_index_ksads_measure]] ~ m_data[[this_index_prs]] * sex + site_id_l + interview_age +
                       pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                     data = m_data, na.action = na.omit, family = "binomial")
   m_p2 <- anova(m_model_2b, m_model_2a, test = "Chisq")
   m_p2_sex <- anova(m_model_2a, m_model_2c, test = "Chisq")
   interaction_row2 = nrow(summary(m_model_2c)$coefficients)
   prs_string_ = colnames(m_data)[this_index_prs]
   measure_variable_name = colnames(m_data)[this_index_ksads_measure]
   measure_display_name_ = as.character(variable_display_map[measure_variable_name,2])
   beta_ = summary(m_model_2a)$coefficients[2, 1]
   p_ = m_p2[2, 5]
   beta_sex_ = summary(m_model_2c)$coefficients[interaction_row2, 1]
   p_sex_ = m_p2_sex[2, 5]
   cat(prs_string_, " ", measure_display_name_, " ",
       beta_, " ", p_, " ",
       beta_sex_, " ", p_sex_, "\n")
   m_df = rbind(m_df, data.frame(PRS=prs_string_, Measure=measure_display_name_,
                                 beta=beta_, p=p_, p_sex=p_sex_))
 }
  for (this_index_mri_measure in indices_mri_normaldist)
  {
    m_model_3a <- lmer(m_data[[this_index_mri_measure]] ~ m_data[[this_index_prs]] + sex + site_id_l + interview_age + smri_vol_cdk_total +
                        pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + (1 | mri_info_deviceserialnumber),
                      data = m_data, na.action = na.omit, REML = FALSE)
    m_model_3b <- lmer(m_data[[this_index_mri_measure]] ~ sex + site_id_l + interview_age + smri_vol_cdk_total +
                        pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + (1 | mri_info_deviceserialnumber),
                      data = m_data, na.action = na.omit, REML = FALSE)
    # test sex interaction
    m_model_3c <- lmer(m_data[[this_index_mri_measure]] ~ m_data[[this_index_prs]] * sex + site_id_l + interview_age + smri_vol_cdk_total +
                         pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + (1 | mri_info_deviceserialnumber),
                       data = m_data, na.action = na.omit, REML = FALSE)
    m_p3 <- anova(m_model_3b, m_model_3a)
    m_p3_sex <- anova(m_model_3a, m_model_3c)
    interaction_row3 = nrow(summary(m_model_3c)$coefficients)
    prs_string_ = colnames(m_data)[this_index_prs]
    measure_variable_name = colnames(m_data)[this_index_mri_measure]
    measure_display_name_ = as.character(variable_display_map[measure_variable_name,2])
    beta_ = summary(m_model_3a)$coefficients[2, 1]
    p_ = m_p3[2, 8]
    beta_sex_ = summary(m_model_3c)$coefficients[interaction_row3, 1]
    p_sex_ = m_p3_sex[2, 8]
    cat(prs_string_, " ", measure_display_name_, " ",
        beta_, " ", p_, " ",
        beta_sex_, " ", p_sex_, "\n")
    m_df = rbind(m_df, data.frame(PRS=prs_string_, Measure=measure_display_name_,
                                  beta=beta_, p=p_, p_sex=p_sex_))
  }
}

# FDR P-value correction
p_indices_aut_behave = 1:16
p_indices_aut_smri = 17:25
p_indices_aut_dti = 26:40
p_indices_anx_behave = 41:56
p_indices_anx_smri = 57:65
p_indices_anx_dti = 66:80

m_df$p_fdr = 0
m_df$p_fdr[p_indices_aut_behave] = p.adjust(m_df$p[p_indices_aut_behave], method="fdr")
m_df$p_fdr[p_indices_aut_smri] = p.adjust(m_df$p[p_indices_aut_smri], method="fdr")
m_df$p_fdr[p_indices_aut_dti] = p.adjust(m_df$p[p_indices_aut_dti], method="fdr")
m_df$p_fdr[p_indices_anx_behave] = p.adjust(m_df$p[p_indices_anx_behave], method="fdr")
m_df$p_fdr[p_indices_anx_smri] = p.adjust(m_df$p[p_indices_anx_smri], method="fdr")
m_df$p_fdr[p_indices_anx_dti] = p.adjust(m_df$p[p_indices_anx_dti], method="fdr")

write.csv(m_df, "G:/My Drive/UCLA/HernandezLab/ABCD_Project_Scripts/output/phewas_stats_mrionly_02-16-24.csv")