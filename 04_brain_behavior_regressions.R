library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(lme4)

data_folder = "G:/My Drive/UCLA/HernandezLab/ABCD_Project_Scripts/data_frames/"
data_file = "abcd_r4_AutAnxMri_PhewasVars_4.Rda"
load(paste0(data_folder, data_file))
m_data = merge_data_filtered

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
for (this_index_mri_measure in indices_mri_normaldist)
{
  measure_mri_variable_name = colnames(m_data)[this_index_mri_measure]
  measure_mri_display_name_ = as.character(variable_display_map[measure_mri_variable_name,2])
  for (this_index_behave_measure in indices_behave_continuum)
  {
    m_model_1a <- glm(m_data[[this_index_behave_measure]]+1 ~ m_data[[this_index_mri_measure]] + sex + site_id_l + interview_age +
                        pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                      data = m_data, na.action = na.omit, family = Gamma(link = "log"))
    m_model_1b <- glm(m_data[[this_index_behave_measure]]+1 ~ sex + site_id_l + interview_age +
                        pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                      data = m_data, na.action = na.omit, family = Gamma(link = "log"))
    m_p1 <- anova(m_model_1b, m_model_1a, test = "F")
    measure_behave_variable_name = colnames(m_data)[this_index_behave_measure]
    measure_behave_display_name_ = as.character(variable_display_map[measure_behave_variable_name,2])
    beta_ = summary(m_model_1a)$coefficients[2, 1]
    p_ = m_p1[2, 6]
    cat(measure_mri_display_name_, " ", measure_behave_display_name_, " ", beta_, " ", p_, "\n")
    m_df = rbind(m_df, data.frame(mri=measure_mri_display_name_, behave=measure_behave_display_name_,
                                  beta=beta_, p=p_))
  }
  for (this_index_ksads_measure in indices_ksads_binary)
  {
    m_model_2a <- glm(m_data[[this_index_ksads_measure]] ~ m_data[[this_index_mri_measure]] + sex + site_id_l + interview_age +
                        pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                      data = m_data, na.action = na.omit, family = "binomial")
    m_model_2b <- glm(m_data[[this_index_ksads_measure]] ~ sex + site_id_l + interview_age +
                        pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                      data = m_data, na.action = na.omit, family = "binomial")
    m_p2 <- anova(m_model_2b, m_model_2a, test = "Chisq")
    measure_behave_variable_name = colnames(m_data)[this_index_ksads_measure]
    measure_behave_display_name_ = as.character(variable_display_map[measure_behave_variable_name,2])
    beta_ = summary(m_model_2a)$coefficients[2, 1]
    p_ = m_p2[2, 5]
    cat(measure_mri_display_name_, " ", measure_behave_display_name_, " ", beta_, " ", p_, "\n")
    m_df = rbind(m_df, data.frame(mri=measure_mri_display_name_, behave=measure_behave_display_name_,
                                  beta=beta_, p=p_))
  }
}

write.csv(m_df, "G:/My Drive/UCLA/HernandezLab/ABCD_Project_Scripts/output/brain_behave_regression_02-17-24.csv")
