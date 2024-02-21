library(ggseg)
library(ggplot2)
library(ggseg3d)
library(ggsegTracula)
library(dplyr)
library(plotly)
reticulate::py_run_string("import sys") # needed to export plotly objects

phewas_stats_path = "G:/My Drive/UCLA/HernandezLab/ABCD_Project_Scripts/output/"
phewas_stats_file = "phewas_stats_mrionly_02-16-24.xlsx"
phewas_stats_aut = readxl::read_excel(paste0(phewas_stats_path, phewas_stats_file), sheet = "aut")
phewas_stats_anx = readxl::read_excel(paste0(phewas_stats_path, phewas_stats_file), sheet = "anx")
phewas_stats_aut$Measurement = as.factor(phewas_stats_aut$Measurement)
phewas_stats_anx$Measurement = as.factor(phewas_stats_anx$Measurement)

measure_info_path = "G:/My Drive/UCLA/HernandezLab/ABCD_Project_Scripts/"
measure_info_file = "abcd_behavior_mri_measures.xlsx"
measures_behavior = readxl::read_excel(paste0(measure_info_path, measure_info_file), sheet = "Behavioral")
measures_mri = readxl::read_excel(paste0(measure_info_path, measure_info_file), sheet = "MRI_ROI")

m_desikian = dk_3d
m_desikian_regions = m_desikian$ggseg_3d[[1]]
m_tracts = tracula_3d
m_tract_regions = m_tracts$ggseg_3d[[1]]

m_custom = m_desikian
uncinate_l = m_tract_regions[15,c(1:5,7)]
uncinate_r = m_tract_regions[16,c(1:5,7)]
m_custom$ggseg_3d[[1]] = rbind(m_custom$ggseg_3d[[1]], uncinate_r)
m_custom$ggseg_3d[[2]] = rbind(m_custom$ggseg_3d[[2]], uncinate_l)
m_custom$ggseg_3d[[3]] = rbind(m_custom$ggseg_3d[[3]], uncinate_r)
m_custom$ggseg_3d[[4]] = rbind(m_custom$ggseg_3d[[4]], uncinate_l)

x_right = m_custom$ggseg_3d[[3]]
x_left = m_custom$ggseg_3d[[4]]
x_right_unc = x_right$mesh[37]
x_left_unc = x_left$mesh[37]

move_uncinate <- function(m_feature, m_hem = "right"){
  x = m_feature[[1]]$vertices$x
  y = m_feature[[1]]$vertices$y
  z = m_feature[[1]]$vertices$z
  
  # flip translate
  y = -1*y - 17
  m_feature[[1]]$vertices$y = y
  z = z - 18
  m_feature[[1]]$vertices$z = z
  if(m_hem == "right") {x = x - 10}
  else {x = x + 10}
  m_feature[[1]]$vertices$x = x
  
  # rotate
  theta = (-20*3.14)/180
  y = y*cos(theta) - z*sin(theta)
  z = y*sin(theta) + z*cos(theta)
  m_feature[[1]]$vertices$y = y
  m_feature[[1]]$vertices$z = z
  
  return(m_feature)
}

m_custom$ggseg_3d[[3]]$mesh[37] = move_uncinate(m_custom$ggseg_3d[[3]]$mesh[37])
m_custom$ggseg_3d[[4]]$mesh[37] = move_uncinate(m_custom$ggseg_3d[[4]]$mesh[37], m_hem="left")

# Loop through brain measures
measure_levels = levels(phewas_stats_aut$Measurement)
measure_levels = measure_levels[measure_levels != "MD"]
measure_levels = measure_levels[measure_levels != "FA"]
measure_levels = measure_levels[measure_levels != "LD"]
measure_levels = measure_levels[measure_levels != "TD"]
m_camera = list(
  up=list(x=0, y=0, z=-2),
  center=list(x=0, y=0, z=0),
  eye=list(x=0, y=1.25, z=-0.45))
for(this_measure in measure_levels)
{
  m_data = subset(phewas_stats_anx, Measurement == this_measure & Region != "Amygdala")
  m_input = m_data[, c("Region", "beta")]
  m_input = rbind(m_input, c("superior parietal", -.001))
  m_input = rbind(m_input, c("precuneus", 0.001))
  colnames(m_input) = c("region","p")
  m_input$p = as.numeric(m_input$p)
  plot_beta = ggseg3d(.data = m_input, atlas = m_custom, colour = "p", palette = c("#f30000", "#ffffff","#00aa00")) %>% 
    pan_camera(m_camera) %>% add_text("beta", x = 20, y = 95, z = 0)
  m_input = m_data[, c("Region", "p")]
  m_input = rbind(m_input, c("superior parietal", 0.99))
  m_input = rbind(m_input, c("precuneus", 0.01))
  colnames(m_input) = c("region","p")
  m_input$p = -log10(as.numeric(m_input$p))
  
  plot_p = ggseg3d(.data = m_input, atlas = m_custom, hemisphere = "left", colour = "p", palette = c("#ffffff", "#666666","#000000")) %>%
    pan_camera(m_camera) %>% add_text("p", x = -20, y = 95, z = 0)
  
  m_fig = plotly::subplot(plot_beta, plot_p)
  m_fig = m_fig %>% layout(title = m_data$Measure[1]) %>% pan_camera(m_camera)
  image_name = paste0(phewas_stats_path, m_data$PRS[1], m_data$Measure[1], ".jpg")
  save_image(m_fig, image_name)
}
for(this_measure in measure_levels)
{
  m_data = subset(phewas_stats_aut, Measurement == this_measure & Region != "Amygdala")
  m_input = m_data[, c("Region", "beta")]
  colnames(m_input) = c("region","p")
  m_input$p = as.numeric(m_input$p)
  plot_aut = ggseg3d(.data = m_input, atlas = m_custom, colour = "p", palette = c("#f30000", "#ffffff","#00aa00")) %>% 
    add_text("aut_prs", x = 20, y = 95, z = 0) %>% pan_camera(m_camera)
  m_data = subset(phewas_stats_anx, Measurement == this_measure & Region != "Amygdala")
  m_input = m_data[, c("Region", "beta")]
  colnames(m_input) = c("region","p")
  m_input$p = as.numeric(m_input$p)
  plot_anx = ggseg3d(.data = m_input, atlas = m_custom, hemisphere = "left", colour = "p", palette = c("#f30000", "#ffffff","#00aa00")) %>% 
    add_text("anx_prs", x = -20, y = 95, z = 0) %>% pan_camera(m_camera)
  
  m_fig = plotly::subplot(plot_aut, plot_anx)
  m_fig = m_fig %>% layout(title = m_data$Measure[1]) %>% pan_camera(m_camera)
  image_name = paste0(phewas_stats_path, "AUT-PRS_ANX_PRS", m_data$Measure[1], ".jpg")
  save_image(m_fig, image_name)
}

##### Stuff for 1 plot at a time
# test
m_camera = list(eye=list(x=-1, y=0.88, z=-1))
m_data = data.frame(c("Uncinate fasciculus", "medial orbitofrontal", "temporal pole"))
m_data$p = 1:3
colnames(m_data)[1] = "region"
ggseg3d(.data = m_data, atlas = m_custom, colour = "p") %>% pan_camera(m_camera)
ggseg3d(.data = m_data, atlas = m_custom, hemisphere = "left", colour = "p") %>% pan_camera(m_camera)

# real data
m_camera = list(
  up=list(x=0, y=0, z=-2),
  center=list(x=0, y=0, z=0),
  eye=list(x=0, y=1.25, z=-0.45))
m_data = subset(phewas_stats_aut, Measurement == 'Vol' & Region != "Amygdala")
m_input = m_data[, c("Region", "beta")]
m_input = rbind(m_input, c("lateral occipital", -2))
m_input = rbind(m_input, c("cuneus", 2))
colnames(m_input) = c("region","p")
m_input$p = as.numeric(m_input$p)
plot_beta = ggseg3d(.data = m_input, atlas = m_custom, colour = "p", palette = c("#f30000", "#ffffff","#00aa00")) %>% 
  pan_camera(m_camera) %>% add_text("beta", x = 20, y = 95, z = 0)

m_input = m_data[, c("Region", "p")]
m_input = rbind(m_input, c("lateral occipital", 0.99))
m_input = rbind(m_input, c("cuneus", 0.01))
m_input$p = -log10(as.numeric(m_input$p))
colnames(m_input) = c("region","p")
plot_p = ggseg3d(.data = m_input, atlas = m_custom, hemisphere = "left", colour = "p", palette = c("#ffffff", "#666666","#000000")) %>%
  pan_camera(m_camera) %>% add_text("p", x = -20, y = 95, z = 0)

# both hemispheres!
m_fig = plotly::subplot(plot_beta, plot_p)
m_fig = m_fig %>% layout(title = 'title') %>% pan_camera(m_camera)
m_fig
plot_beta
plot_p
##### Begin saving PDF #####
x = save_image(m_fig, "test.png")
save_image(plot_p, "test2.png")
##### Complete PDF #####