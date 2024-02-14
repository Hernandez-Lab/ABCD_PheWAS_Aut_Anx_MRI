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
mri_file = "abcd_r4_AutAnxMri_PhewasVars_AncestryFiltered-EURonly_LongitudinalFiltered-BaselineOnly_MRIFiltered-GoodScansOnlyAndExtremesRemovedAndHemispheresAveraged.Rda"

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

#prs_selector = c("AUT_PRS", "ANX_PRS")
indices_prs = 82:83
indices_behave_continuum = c(19:21,57) # glm gamma log, f test
indices_ksads_binary = 22:33 # glm binomial, chisq test
indices_mri_normaldist = 58:81 # lmer, normal anova

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
    cat(colnames(m_data)[this_index_prs], " ", colnames(m_data)[this_index_behave_measure], " ",
        summary(m_model_1a)$coefficients[2, 1], " ", m_p1[2, 6], " ",
        summary(m_model_1c)$coefficients[interaction_row1, 1], " ", m_p1_sex[2, 6], "\n")
    m_df = rbind(m_df, data.frame(PRS=colnames(m_data)[this_index_prs], Measure=colnames(m_data)[this_index_behave_measure],
        sex_int_p=m_p1_sex[2, 6], beta=summary(m_model_1a)$coefficients[2, 1], p=m_p1[2, 6]))
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
   m_p2_sex <- anova(m_model_2a, m_model_2c, test = "F")
   interaction_row2 = nrow(summary(m_model_2c)$coefficients)
   cat(colnames(m_data)[this_index_prs], " ", colnames(m_data)[this_index_ksads_measure], " ",
       summary(m_model_2a)$coefficients[2, 1], " ", m_p2[2, 5], " ",
       summary(m_model_2c)$coefficients[interaction_row2, 1], " ", m_p2_sex[2, 5], "\n")
   m_df = rbind(m_df, data.frame(PRS=colnames(m_data)[this_index_prs], Measure=colnames(m_data)[this_index_ksads_measure],
      sex_int_p=m_p2_sex[2, 5], beta=summary(m_model_2a)$coefficients[2, 1], p=m_p2[2, 5]))
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
    cat(colnames(m_data)[this_index_prs], " ", colnames(m_data)[this_index_mri_measure], " ",
        summary(m_model_3a)$coefficients[2, 1], " ", m_p3[2, 8], " ",
        summary(m_model_3c)$coefficients[interaction_row3, 1], " ", m_p2_sex[2, 8], "\n")
    m_df = rbind(m_df, data.frame(PRS=colnames(m_data)[this_index_prs], Measure=colnames(m_data)[this_index_mri_measure],
       sex_int_p=m_p3_sex[2, 8], beta=summary(m_model_3a)$coefficients[2, 1], p=m_p3[2, 8]))
  }
}

#write.csv(m_df, "G:/My Drive/UCLA/HernandezLab/ABCD_Project_Scripts/phewas_stats.csv")

# Mediation/moderation analysis
# PRS > Brain > Behavior
# (AUT PRS) > () > (Social Anxiety {KSADS or CBCL})
# Keep in mind: most measures are not compatable with linear mediation analysis....

# MarginalMediation package, what does the result mean?
fitY <- glm(aut_agg_score ~ dmri_rsirndgm_cdk_tp + AUT_PRS + sex + site_id_l + interview_age + pc1 + pc2 + pc3,
            data = m_data, na.action = na.omit, family = Gamma(link = "log"))
fitM <- glm(dmri_rsirndgm_cdk_tp ~ AUT_PRS + sex + site_id_l + interview_age + pc1 + pc2 + pc3,
            data = fitY$data, na.action = na.omit, family = Gamma(link = "log"))
mediation_analysis <- mma(fitY, fitM, ind_effects = c("AUT_PRS-dmri_rsirndgm_cdk_tp"), boot = 500)
m_perc_med = perc_med(mediation_analysis, c("AUT_PRS-dmri_rsirndgm_cdk_tp"))
# AUT_PRS-dmri_rsirndgm_cdk_tp = 3.851974 
m_perc_med_2 = perc_med(mediation_analysis, c("AUT_PRS"))
# AUT_PRS = 3.851974 

# GEEmediate
gee_mediation_analysis =
  GEEmediate(aut_agg_score ~ smri_vol_cdk_total + dmri_rsirndgm_cdk_tp + AUT_PRS + sex + site_id_l + interview_age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
  exposure = "AUT_PRS", mediator = "dmri_rsirndgm_cdk_tp", df = m_data, family = Gamma(link = "log"))
# Mediation Proportion:2.5%, p=0.46 for one-sided test for mediation
# Could actually have something here

# ctrl test
ctrl_mediation_analysis =
  GEEmediate(aut_agg_score ~ smri_vol_cdk_total + dmri_rsirndgm_cdk_tp + AUT_PRS + sex + site_id_l + interview_age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
             exposure = "AUT_PRS", mediator = "interview_age", df = m_data, family = Gamma(link = "log"))
# Mediation Proportion:0.57%, p=0.49 for one-sided test for mediation

# ctrl test 2
ctrl_mediation_analysis_2 =
  GEEmediate(aut_agg_score ~ smri_vol_cdk_total + AUT_PRS + sex + site_id_l + interview_age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
             exposure = "AUT_PRS", mediator = "smri_vol_cdk_total", df = m_data, family = Gamma(link = "log"))
# Mediation Proportion:0.05 %, p=0.5