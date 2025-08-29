
#------------------------------------------------------------------------------------------------------------------------------------------------#
# Analyses for "Shared and specific blood biomarkers for multimorbidity"
# 1. Figure 2. Heatmap of multinomial LASSO coefficients
# 2. Figure 3. Trajectories of disease accumulation and LASSO-selected biomarkers associated with the rate of disease accumulation
# 3. Figure 4. PCA of LASSO-selected biomarkers across multimorbidity measures
# R version 4.4.2
#------------------------------------------------------------------------------------------------------------------------------------------------#

# load packages
library(tidyverse)
library(ggplot2)
library(ggprism)
library(corrplot)
library(readr)

# ------------> 1. Figure 2. Heatmap of multinomial LASSO coefficients -----------------

load("df_selected_compact.RData") # dataframe with LASSO-selected biomarkers (script 02)

# customize labels
custom_label_p <- c("Unspecific" = "Unspecific",
                    "Neuropsychiatric" = "Neuropsychiatric",
                    "Psychiatric/Respiratory" = "Psychiatric \n Respiratory",
                    "Sensory impairment/Anaemia" = "Sensory impairment \n Anaemia",
                    "Cardiometabolic" = "Cardiometabolic")

df_selected_compact <- df_selected_compact %>%
  group_by(Biomarkers) %>%
  mutate(Code_frequency = n()) %>%
  ungroup()

df_selected_compact <- df_selected_compact %>% complete(MM_pattern, Biomarkers, fill = list(mean_coef = NA))

df_selected_compact <- df_selected_compact %>% 
  mutate(Biomarkers = recode(Biomarkers, z_CPeptide = 'CPEP', z_CystatinC = 'CysC', z_GDF15 = 'GDF15', z_HbA1C = 'HbA1c', 
                             z_Insulin = 'Insulin', z_Leptin = 'LEP', z_NCadherin = 'NCAD',  z_VCAM1 = 'VCAM1', z_abeta_ratio = 'Aβ42/40', 
                             z_cholesterol = 'TC', z_creatinine = 'Cr', z_folic_acid = 'FA',  z_gamma_gt = 'GGT', z_hemoglobin = 'Hb', 
                             z_vitamin_b12 = 'VitB12',z_IGFBP1 = 'IGFBP1', z_albumin = 'ALB', z_hIL8 = 'IL8', z_nfl = 'NfL', z_T4 = 'T4',  
                             z_ptau181 = 'ptau181', z_Adiponectin = 'Adipo',z_ESelectin = 'ESel', z_ICAM1 = 'ICAM1'))


df_selected_compact <- df_selected_compact %>% 
  mutate(variable = factor(MM_pattern)) %>%
  mutate(variable = fct_relevel(MM_pattern, c("Unspecific", "Neuropsychiatric", "Psychiatric/Respiratory", 
                                              "Sensory impairment/Anaemia", "Cardiometabolic")))                 

df_selected_compact <- df_selected_compact %>%
  mutate(biomarker = factor(Biomarkers, levels = c(
    "TC", "LEP", "HbA1c", "GDF15", "Cr", "CysC", "FA", "Insulin", "CPEP", "Hb", "Aβ42/40", "GGT", "VCAM1", "NfL", "NCAD", "IGFBP1", 
    "T4", "IL8", "ESel", "ptau181", "VitB12", "ICAM1", "Adipo", "ALB")))

df_selected_compact <- df_selected_compact %>%
  mutate(biomarker = factor(biomarker, levels = rev(levels(biomarker))))

ggplot(df_selected_compact, aes(x = variable, y = biomarker, fill = mean_coef)) + 
  geom_tile(color = "darkgray", linewidth = 0.5) + 
  scale_x_discrete(labels = custom_label_p, position = "top")+
  scale_fill_gradientn(
    colors = c("#0367AD", "lightblue","#F5F5F5", "#E06007" , "#D1061B"), 
    values = scales::rescale(c(-0.5, -0.25, 0, 0.25, 0.5)),
    limits = c(-0.5, 0.50), na.value = "white", 
    guide = guide_colourbar(barwidth = 40, barheight = 2, title.position = "top",  
                            title.hjust = 0.5, 
                            nbin = 10 ))+
  theme_prism() + 
  labs(title = "", x = "", y = "Biomarkers", fill = "adjusted LASSO β coefficients") + 
  theme(axis.text.x = element_text(angle = 90, size = 30, face = "bold"),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(angle = 0,  size = 35),
        axis.title.y = element_text(size = 40),
        legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.text = element_text(size = 25), 
        legend.key.size = unit(1.5, "cm"), 
        legend.title = element_text(size = 30),
        panel.border = element_rect(color = "gray", linewidth = 1, fill = "transparent"), 
        panel.background = element_blank(),
        panel.spacing = unit(0, "lines"),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"))



# ---------> 2. Figure 3. Trajectories of disease accumulation and LASSO-selected biomarkers associated with the rate of disease accumulation ----------

# Panel A
# for this plot LMM result (script 03) are needed
set.seed(1234)
plot_data <- db.long %>%
  filter(ID %in% sample(unique(ID), 50)) %>%
  mutate(predicted = predict(model0.lmer, newdata = ., type = "response")) %>%
  dplyr::select(ID, wave, predicted)

plot_lmm <- ggplot(plot_data, aes(x = wave, y = predicted, group = factor(ID), color = factor(ID))) +
  geom_line(alpha = 0.8, linewidth = 1) +
  labs(title = "", 
       x = "Follow-up time",
       y = "Predicted diseases accumulation") +
  scale_x_continuous(breaks = seq(min(db.long$wave), max(db.long$wave), by = 3),
                     labels = c("Baseline", "3 years", "6 years", "9 years", "12 years", "15 years")) +
  scale_y_continuous(limits = c(0,17.5)) +
  theme_prism() +
  theme(legend.position = "none",
        axis.title.y = element_text(face = "bold", margin = margin(0,10,0,0), size = 18, color = 'black'),
        axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = 20, color = 'black'),
        axis.text.x = element_text(size = 18),
        axis.text.y = element_text(size = 18),
        plot.title = element_text(hjust = 0.5), 
        plot.margin = margin(5.5, 20, 5.5, 5.5, unit = "pt"))

# Panel B
load("lasso_non0_long.RData") #data of LASSO-selected biomarkers for rate of disease accumulation

lasso_non0_long %>%
  ggplot(aes(x = reorder(Biomarkers, Coefficient), y = Coefficient, fill = Coefficient)) +
  geom_bar(stat = "identity", color = "black") + 
  scale_fill_gradientn(
    colors = c("#0367AD", "lightblue","#F5F5F5", "#E06007" , "#D1061B"),
    values = scales::rescale(c(-0.005, -0.001, 0, 0.001, 0.010, 0.025)),
    limits = c(-0.005, 0.025))+
  scale_y_continuous(limits = c(-0.005,0.026), breaks = seq(-0.005, 0.028, by = 0.005))+
  coord_flip() +  # Flip coordinates to make it easier to read variable names
  labs(title = "",
       x = "Biomarkers",
       y = "Adjusted LASSO β coefficients") +
  theme_prism()+
  theme(axis.text.x = element_text(size = 18), #18
        axis.text.y = element_text(size = 15), #18
        axis.title.x = element_text(size = 18), #20
        axis.title.y = element_text(size = 18))+ #20
  guides(fill = "none")


# Panel D
blsa_se <- read_csv("BLSA_se.csv") #load BLSA vector with SE
blsa_se$cohort <- "BLSA"
blsa_se <- blsa_se %>% rename(se = x)

se_snack <- as.data.frame(se_snak) # SNAC-K vector with SE (script 04)
se_snack <- se_snack %>% rename(se = se_snak)
se_snack$cohort <- "SNAC-K"

se_data <- rbind(se_snack, blsa_se)

ggplot(se_data, aes(x = se, fill = cohort))+
  geom_density(alpha = 1) +
  scale_fill_manual(name = "Study cohort", 
                    labels = c("BLSA", "SNAC-K"),
                    values = c("#F92F2E", "#61AEDE")) +
  scale_color_manual(values = c("BLSA" = "#F92F2E", "SNAC-K" = "#61AEDE")) +
  scale_linetype_manual(values = c("BLSA" = "solid", "SNAC-K" = "solid"))
labs(title = "", x = "Square Errors", y = "Density") +
  theme_prism() +
  theme(
    legend.title = element_text(face = "bold", size = 18),
    legend.text = element_text(size = 16),
    axis.title = element_text(face = "bold", size = 18),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 18, face = "bold"),
    legend.position = c(0.80, 0.85))

# ---------> 3. Figure 4. PCA of LASSO-selected biomarkers across multimorbidity measures -----------------

load("heat_dis_data.Rdata") # PCA data for baseline number of chronic disease (script 01)
load("heat_data.Rdata") # PCA data for rate of disease accumulation (script 03)
load("unspe_heat_data") # PCA data for unspecific pattern (script 02)
load("cardio_heat_data") # PCA data for cardiometabolic pattern (script 02)
load("psych_heat_data") # PCA data for psychiatric/respiratory pattern (script 02)
load("senso_heat_data") # PCA data for sensory impairment/anemia pattern (script 02)
load("neuro_heat_data") # PCA data for neuropsychiatric pattern (script 02)


db_heat_overall <- rbind(unspe_heat_data, cardio_heat_data, psych_heat_data, senso_heat_data, neuro_heat_data, heat_dis_data, heat_data) #combine dataset

# customize labels
custom_label_p2 <- c("Number of chronic diseases" = "Number of \n chronic diseases",
                     "Unspecific" = "Unspecific", 
                     "Neuropsychiatric" = "Neuropsychiatric", 
                     "Psychiatric/Respiratory" = "Psychiatric \n Respiratory", 
                     "Sensory impairment/Anaemia" = "Sensory impairment \n Anaemia", 
                     "Cardiometabolic" = "Cardiometabolic", 
                     "Speed of disease accumulation" = "Speed of disease \n accumulation")

db_heat_overall<- db_heat_overall %>% 
  mutate(Pattern = factor(Pattern)) %>%
  mutate(Pattern = fct_relevel(Pattern, c("Number of chronic diseases","Unspecific", "Neuropsychiatric", "Psychiatric/Respiratory", 
                                          "Sensory impairment/Anaemia", "Cardiometabolic", "Speed of disease accumulation"))) %>%
  arrange(Pattern)

heat_LASSO_overall <- ggplot(db_heat_overall, aes(x = xName, y = yName, fill = corr)) + 
  geom_tile(color = "darkgray", linewidth = 0.5) + 
  scale_fill_gradientn(
    colors = c("#F5F5F5","#F7D361", "#E06007" , "#D1061B"), 
    values = scales::rescale(c(0, 0.5, 1)),
    limits = c(0, 1), na.value = "white",
    guide = guide_colourbar(barwidth = 40, barheight = 2, title.position = "top",  
                            title.hjust = 0.5,
                            nbin = 10 ))+
  scale_y_discrete(limits = rev(levels(db_heat_all2$yName)))+
  theme_prism() + 
  facet_grid(~Pattern, scales = "free_x", labeller = labeller(Pattern = custom_label_p2))+
  labs(title = "", x = "Principal components (PC)", y = "Biomarkers", fill = "Cos2") + 
  theme(axis.text.x = element_text(angle = 90,  size = 28, face = "plain"),
        axis.text.y = element_text(angle = 0,  size = 28),
        axis.title.x = element_text(size = 35),
        axis.title.y = element_text(size = 35),
        strip.text.x = element_text(angle = 90,  size = 25, face = "bold"),
        legend.position = "bottom", 
        legend.direction = "horizontal", 
        legend.title = element_text(size = 25, face = "bold"),
        legend.text = element_text(size = 25), 
        legend.key.size = unit(1.5, "cm"),
        panel.border = element_rect(color = "gray", linewidth = 1, fill = "transparent"), 
        panel.background = element_blank(),
        panel.grid.major.y = element_line(color = "lightgray", linewidth = 0.5),  # Customize major grid lines for y-axis
        panel.spacing = unit(0.5, "lines"),
        plot.margin = margin(t = 20, r = 20, b = 20, l = 20, unit = "pt"))

