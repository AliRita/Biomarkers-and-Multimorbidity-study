#------------------------------------------------------------------------------------------------------------------------------------------------#
# Analyses for "Shared and specific blood biomarkers for multimorbidity"
# 1. External validation (same code for SNAC-K and BLSA)
# R version 4.4.2
#------------------------------------------------------------------------------------------------------------------------------------------------#


# load packages
library(tidyverse)
library(magrittr)
library(glmnet)
library(ggprism)
library(RColorBrewer)
library(corrplot)

# load function
source("fun_validation.R")


load("db.RData") # study dataset 
load("02_individual_slope_diseases.RData") #obtained using linear mixed model


## ------- 1. Compare correlation matrix among biomarkers ------


cor_dat <- db %>% select(c(128,92,111,112,123,94,95)) # select biomarkers variables

cor_dat <- cor_dat %>% rename(c(ALB = albumin, CysC = CystatinC, LEP = Leptin, GGT = gamma_gt))

cor_plot <- cor(cor_dat, method = "spearman")

var_order <- sort(colnames(cor_plot))

cor_plot_ordered <- cor_plot[var_order, var_order]

corrplot(cor_plot_ordered, 
         method = "color", 
         type = "full", 
         col = rev(COL2('RdYlBu', 10)),
         tl.col = "black", 
         title = "",
         is.corr = TRUE, 
         addCoef.col = "black",  
         number.cex = 0.9,        
         order = "original", 
         addgrid.col = "darkgray")

dat_slope <- db.long %>% distinct(ID,.keep_all = F)
dat_slope$beta_crude <-ranef(model0.lmer)$ID$`wave`


# ------------> 2. Calculate predicted rate of disease accumulation based on LASSO formula -----------------


beta_pred <- predict_dis_acc_rate(data = db.def,idvar="ID",agevar="age_w1_s")

# Join with the "observed" rate of disease accumulation obtained from script 02

dat_slope %<>% full_join(beta_pred) 

# ------------> 3. Calculate MSE and obtain plot square errors -----------------

mse(dat_slope)

se_snak=plot_square_error(dat_slope,"SNAC-K")





