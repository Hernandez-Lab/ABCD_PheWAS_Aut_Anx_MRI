library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(lme4)
library(MarginalMediation)
library(GEEmediate)

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

##### AUT_PRS*ANX_PRS moderation/interaction #####
prs_string_ = "AUT_PRS*ANX_PRS"
m_df <- data.frame()
for (this_index_behave_measure in indices_behave_continuum)
{
  m_model_1a <- glm(m_data[[this_index_behave_measure]]+1 ~ m_data[[indices_prs[1]]] * m_data[[indices_prs[2]]] + sex + site_id_l + interview_age +
                      pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                    data = m_data, na.action = na.omit, family = Gamma(link = "log"))
  m_model_1b <- glm(m_data[[this_index_behave_measure]]+1 ~ m_data[[indices_prs[1]]] + m_data[[indices_prs[2]]]  + sex + site_id_l + interview_age +
                      pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                    data = m_data, na.action = na.omit, family = Gamma(link = "log"))
  m_p1 <- anova(m_model_1b, m_model_1a, test = "F")
  interaction_row1 = nrow(summary(m_model_1a)$coefficients)
  measure_variable_name = colnames(m_data)[this_index_behave_measure]
  measure_display_name_ = as.character(variable_display_map[measure_variable_name,2])
  beta_ = summary(m_model_1a)$coefficients[interaction_row1, 1]
  p_ = m_p1[2, 6]
  cat(prs_string_, " ", measure_display_name_, " ", beta_, " ", p_, "\n")
  m_df = rbind(m_df, data.frame(PRS=prs_string_, Measure=measure_display_name_,
                                beta=beta_, p=p_))
}
for (this_index_ksads_measure in indices_ksads_binary)
{
  m_model_2a <- glm(m_data[[this_index_ksads_measure]] ~ m_data[[indices_prs[1]]] * m_data[[indices_prs[2]]] + sex + site_id_l + interview_age +
                      pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                    data = m_data, na.action = na.omit, family = "binomial")
  m_model_2b <- glm(m_data[[this_index_ksads_measure]] ~ m_data[[indices_prs[1]]] + m_data[[indices_prs[2]]] + sex + site_id_l + interview_age +
                      pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                    data = m_data, na.action = na.omit, family = "binomial")
  m_p2 <- anova(m_model_2b, m_model_2a, test = "Chisq")
  interaction_row2 = nrow(summary(m_model_2a)$coefficients)
  measure_variable_name = colnames(m_data)[this_index_ksads_measure]
  measure_display_name_ = as.character(variable_display_map[measure_variable_name,2])
  beta_ = summary(m_model_2a)$coefficients[interaction_row2, 1]
  p_ = m_p2[2, 5]
  cat(prs_string_, " ", measure_display_name_, " ", beta_, " ", p_, " ", "\n")
  m_df = rbind(m_df, data.frame(PRS=prs_string_, Measure=measure_display_name_,
                                beta=beta_, p=p_))
}
for (this_index_mri_measure in indices_mri_normaldist)
{
  m_model_3a <- lmer(m_data[[this_index_mri_measure]] ~ m_data[[indices_prs[1]]] * m_data[[indices_prs[2]]] + sex + site_id_l + interview_age + smri_vol_cdk_total +
                       pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + (1 | mri_info_deviceserialnumber),
                     data = m_data, na.action = na.omit, REML = FALSE)
  m_model_3b <- lmer(m_data[[this_index_mri_measure]] ~ m_data[[indices_prs[1]]] + m_data[[indices_prs[2]]] + sex + site_id_l + interview_age + smri_vol_cdk_total +
                       pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10 + (1 | mri_info_deviceserialnumber),
                     data = m_data, na.action = na.omit, REML = FALSE)
  m_p3 <- anova(m_model_3b, m_model_3a)
  interaction_row3 = nrow(summary(m_model_3a)$coefficients)
  measure_variable_name = colnames(m_data)[this_index_mri_measure]
  measure_display_name_ = as.character(variable_display_map[measure_variable_name,2])
  beta_ = summary(m_model_3a)$coefficients[interaction_row3, 1]
  p_= m_p3[2, 8]
  cat(prs_string_, " ", measure_display_name_, " ", beta_, " ", p_, " ", "\n")
  m_df = rbind(m_df, data.frame(PRS=prs_string_, Measure=measure_display_name_,
                                beta=beta_, p=p_))
}

# FDR P-value correction
p_indices_behave = 1:16
p_indices_smri = 17:25
p_indices_dti = 26:40

m_df$p_fdr = 0
m_df$p_fdr[p_indices_behave] = p.adjust(m_df$p[p_indices_behave], method="fdr")
m_df$p_fdr[p_indices_smri] = p.adjust(m_df$p[p_indices_smri], method="fdr")
m_df$p_fdr[p_indices_dti] = p.adjust(m_df$p[p_indices_dti], method="fdr")

write.csv(m_df, "G:/My Drive/UCLA/HernandezLab/ABCD_Project_Scripts/output/phewas_interaction_stats_02-17-24.csv")

##### Brain moderation/interaction #####
m_df <- data.frame()
for (this_index_prs in indices_prs)
{
  measure_prs_display_name = colnames(m_data)[this_index_prs]
  for (this_index_mri_measure in indices_mri_normaldist)
  {
    measure_mri_variable_name = colnames(m_data)[this_index_mri_measure]
    measure_mri_display_name = as.character(variable_display_map[measure_mri_variable_name,2])
    # continuous
    for (this_index_behave_measure in indices_behave_continuum)
    {
      print(this_index_behave_measure)
      measure_behave_variable_name = colnames(m_data)[this_index_behave_measure]
      measure_behave_display_name = as.character(variable_display_map[measure_behave_variable_name,2])
      m_model_1a <- glm(m_data[[this_index_behave_measure]]+1 ~ m_data[[this_index_prs]] * m_data[[this_index_mri_measure]] + smri_vol_cdk_total +
                          sex + site_id_l + interview_age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                          data = m_data, na.action = na.omit, family = Gamma(link = "log"))
      m_model_1b <- glm(m_data[[this_index_behave_measure]]+1 ~ m_data[[this_index_prs]] + m_data[[this_index_mri_measure]] + smri_vol_cdk_total +
                          sex + site_id_l + interview_age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                        data = m_data, na.action = na.omit, family = Gamma(link = "log"))
      m_p1 <- anova(m_model_1b, m_model_1a, test = "F")
      interaction_row1 = nrow(summary(m_model_1a)$coefficients)
      beta_ = summary(m_model_1a)$coefficients[interaction_row1, 1]
      p_ = m_p1[2, 6]
      cat(measure_mri_display_name, " ", measure_behave_display_name, " ", beta_, " ", p_, "\n")
      m_df = rbind(m_df, data.frame(prs=measure_prs_display_name, mri=measure_mri_display_name,
                                    behave=measure_behave_display_name, beta=beta_, p=p_))
    }
    # binomial
    for (this_index_ksads_measure in indices_ksads_binary)
    {
      print(this_index_ksads_measure)
      measure_behave_variable_name = colnames(m_data)[this_index_ksads_measure]
      measure_behave_display_name = as.character(variable_display_map[measure_behave_variable_name,2])
      m_model_2a <- glm(m_data[[this_index_ksads_measure]] ~ m_data[[this_index_prs]] * m_data[[this_index_mri_measure]] + smri_vol_cdk_total +
                          sex + site_id_l + interview_age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                        data = m_data, na.action = na.omit, family = "binomial")
      m_model_2b <- glm(m_data[[this_index_ksads_measure]] ~ m_data[[this_index_prs]] + m_data[[this_index_mri_measure]] + smri_vol_cdk_total +
                          sex + site_id_l + interview_age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
                        data = m_data, na.action = na.omit, family = "binomial")
      m_p2 <- anova(m_model_2b, m_model_2a, test = "Chisq")
      interaction_row2 = nrow(summary(m_model_2a)$coefficients)
      beta_ = summary(m_model_2a)$coefficients[interaction_row2, 1]
      p_ = m_p2[2, 5]
      cat(measure_mri_display_name, " ", measure_behave_display_name, " ", beta_, " ", p_, "\n")
      m_df = rbind(m_df, data.frame(prs=measure_prs_display_name, mri=measure_mri_display_name,
                                    behave=measure_behave_display_name, beta=beta_, p=p_))
    }
  }
}

write.csv(m_df, "G:/My Drive/UCLA/HernandezLab/ABCD_Project_Scripts/output/phewas_brainmoderation_stats_02-17-24.csv")

##### Brain mediation #####
m_df <- data.frame()
# continuous
for (this_index_prs in indices_prs)
{
  prs_string_ = colnames(m_data)[this_index_prs]
  for (this_index_behave_measure in indices_behave_continuum)
  {
    behave_string_ = colnames(m_data)[this_index_behave_measure]
    for (this_index_mri_measure in indices_mri_normaldist)
    {
      print(this_index_mri_measure)
      mri_string_ = colnames(m_data)[this_index_mri_measure]
      m_data2 = m_data
      m_data2[[this_index_behave_measure]] = m_data2[[this_index_behave_measure]]+1
      colnames(m_data2)[this_index_prs] = "x"
      colnames(m_data2)[this_index_behave_measure] = "y"
      colnames(m_data2)[this_index_mri_measure] = "z"
      m_model_1a <- GEEmediate(y ~ x + z +
                          smri_vol_cdk_total + sex + site_id_l + interview_age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8,
                          exposure = "x", mediator = "z",
                          df = m_data2, na.action = na.omit, family = Gamma(link = "log"))
      mediation_prop = m_model_1a[["pm"]]
      mediation_p = m_model_1a[["pm.pval"]]
      indirect_effect = m_model_1a[["nie"]]
      indirect_effect_p = m_model_1a[["nie.pval"]]
      direct_effect = m_model_1a[["nde"]]
      direct_effect_p = m_model_1a[["nde.pval"]]
      if(is.null(mediation_prop)) { mediation_prop = 0.0 }
      if(is.null(mediation_p)) { mediation_p = 1.0 }
      if(is.null(indirect_effect)) { indirect_effect = 0.0 }
      if(is.null(indirect_effect_p)) { indirect_effect_p = 1.0 }
      if(is.null(direct_effect)) { direct_effect = 0.0 }
      if(is.null(direct_effect_p)) { direct_effect_p = 1.0 }
      total_effect = indirect_effect + direct_effect
      cat(prs_string_, " ~ ", behave_string_, " *mediation effect of* ", mri_string_, "\n", mediation_prop, ", ", mediation_p, "\n\n")
      m_df = rbind(m_df, data.frame(PRS=prs_string_, behave=behave_string_, mri=mri_string_,
                                    med=mediation_prop, p=mediation_p,
                                    ie=indirect_effect, ie_p=indirect_effect_p,
                                    de=direct_effect, de_p=direct_effect_p, te=total_effect))
    }
  }
}
# binomial
## TODO: Not working
for (this_index_prs in indices_prs)
{
  prs_string_ = colnames(m_data)[this_index_prs]
  for (this_index_behave_measure in indices_ksads_binary)
  {
    behave_string_ = colnames(m_data)[this_index_behave_measure]
    for (this_index_mri_measure in indices_mri_normaldist)
    {
      print(this_index_mri_measure)
      mri_string_ = colnames(m_data)[this_index_mri_measure]
      m_data2 = m_data
      colnames(m_data2)[this_index_prs] = "x"
      colnames(m_data2)[this_index_behave_measure] = "y"
      colnames(m_data2)[this_index_mri_measure] = "z"
      fitY <- glm(y ~ x + z +
                    sex + site_id_l + interview_age + pc1 + pc2 + pc3 + smri_vol_cdk_total,
                    data = m_data2, na.action = na.omit, family = "binomial")
      fitM <- glm(z ~ x +
                    sex + site_id_l + interview_age + pc1 + pc2 + pc3 + smri_vol_cdk_total,
                    data = m_data2, na.action = na.omit)
      ind_path = paste0(prs_string_, "-", mri_string_)
      mediation_analysis <- mma(fitY, fitM, ind_effects = c("x-z"), boot = 10)
      percent_mediation = perc_med(mediation_analysis, c("x-z"))
      indirect_effects = mma_ind_effects(mediation_analysis)[[3]]
      direct_effects = mma_dir_effects(mediation_analysis)[[3]]
      if(is.null(percent_mediation)) { percent_mediation = 0.0 }
      if(is.null(indirect_effects)) { indirect_effects = 0.0 }
      if(is.null(direct_effects)) { direct_effects = 0.0 }
      total_effect = indirect_effects + direct_effects
      cat(prs_string_, " ~ ", behave_string_, " *mediation effect of* ", mri_string_, "\n", percent_mediation, "\n\n")
      m_df = rbind(m_df, data.frame(PRS=prs_string_, behave=behave_string_, mri=mri_string_,
                                    med=percent_mediation, p=1.0,
                                    ie=indirect_effects, ie_p=1.0,
                                    de=direct_effects, de_p=1.0, te=total_effect))
    }
  }
}

write.csv(m_df, "G:/My Drive/UCLA/HernandezLab/ABCD_Project_Scripts/output/phewas_mediation_stats_02-18-24.csv")
  
##### testing area

# Mediation/moderation analysis
# PRS > Brain > Behavior
# (AUT PRS) > () > (Social Anxiety {KSADS or CBCL})
# Keep in mind: most measures are not compatable with linear mediation analysis....

# MarginalMediation package, what does the result mean?
fitY <- glm(ksads_8_864_p ~ AUT_PRS + dmri_rsirndgm_cdk_tp + sex + site_id_l + interview_age + pc1 + pc2 + pc3,
            data = m_data, na.action = na.omit, family = "binomial")
fitM <- glm(dmri_rsirndgm_cdk_tp ~ AUT_PRS + sex + site_id_l + interview_age + pc1 + pc2 + pc3,
            data = fitY$data, na.action = na.omit, family = Gamma(link = "log"))
mediation_analysis <- mma(fitY, fitM, ind_effects = c("AUT_PRS-dmri_rsirndgm_cdk_tp", "pc1-dmri_rsirndgm_cdk_tp"), boot = 10)
perc_med(mediation_analysis, c("AUT_PRS-dmri_rsirndgm_cdk_tp"))
perc_med(mediation_analysis, c("pc1-dmri_rsirndgm_cdk_tp"))
mma_ind_effects(mediation_analysis)
mma_dir_effects(mediation_analysis)
mma_formulas(mediation_analysis)
mma_check(mediation_analysis)

m_model_1a <- GEEmediate(y ~ x + z +
                           smri_vol_cdk_total + sex + site_id_l + interview_age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8,
                         exposure = "x", mediator = "z",
                         df = m_data2, na.action = na.omit, family = binomial(link = "logit"))

# GEEmediate
gee_mediation_analysis =
  GEEmediate(aut_agg_score ~ smri_vol_cdk_total + dmri_rsirndgm_cdk_tp + AUT_PRS + sex + site_id_l + interview_age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8 + pc9 + pc10,
  exposure = "AUT_PRS", mediator = "dmri_rsirndgm_cdk_tp", df = m_data, family = Gamma(link = "log"))
# Mediation Proportion:2.5%, p=0.46 for one-sided test for mediation
# Could actually have something here

m_model_1a = GEEmediate(aut_agg_score ~ AUT_PRS + dmri_rsirndgm_cdk_tp + smri_vol_cdk_total + sex + site_id_l + interview_age + pc1 + pc2 + pc3 + pc4 + pc5 + pc6 + pc7 + pc8,
    exposure = "smri_vol_cdk_total", mediator = "pc1", df = m_data, na.action = na.omit, family = Gamma(link = "log"))

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