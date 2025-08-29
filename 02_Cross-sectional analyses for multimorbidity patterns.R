
#------------------------------------------------------------------------------------------------------------------------------------------------#
# Analyses for "Shared and specific blood biomarkers for multimorbidity"
# 1. Obtain multimorbidity patterns with LCA
# 2. Multinomial LASSO for baseline multimorbidity patterns
# 3. PCA
# R version 4.4.2
#------------------------------------------------------------------------------------------------------------------------------------------------#

# load packages
library(glmnet)
library(tidyverse)
library(factoextra)
library(corrplot)

# load function
source("utils_fun_lasso.R")
source("utils_LCA.R")


# ------------------------------------------------------------------------------------------------------------------------------------------------#
# (1) Run LCA, select the number of MM patterns, goodness of fit and interpretation
# (2) Define parameters for multinomial LASSO, run the 1000 models, and generate non-zero coefficients dataframe
# (3) PCA of LASSO-selected biomarkers for each Multimorbidity patterns
# (4) Generate dataframes with PCA data for each multimorbidity pattern for plotting purpose (see script number XXXXX)
# ------------------------------------------------------------------------------------------------------------------------------------------------#

# ----- Identification of MM patterns ----
# Select conditions with at least 2% prevalence
threshold <- 0.02 #2%
disease_names <- select_conditions(dat_long = dat_long,
                                   threshold = threshold)

disease_names

# Prepare data disease matrix for LCA on multimorbid individuals
X <- dat%>% dplyr::select(2:61)  %>% mutate_all(function(x)x+1)
ndis <- apply(X-1,1,sum)
table(ndis>=2)
X <- X[ndis>=2,]

#Run LCA with different number of classes (from 2 to 6)
set.seed(202112)


res <- select_number_LCA(nclasses=2:6,
                         X = X,
                         conditions =disease_names,
                         nrep=200)

#For CV entropy and accuracy
res_cv <- select_number_LCA(nclasses=2:6,
                            X = X,
                            conditions =disease_names,
                            nrep=50,
                            cvfolds=5)


# Figure "Extended Data 4"

dat_metrics <- res_cv[[1]] %>% 
left_join(res$metrics %>%
            as.data.frame() %>%
            mutate(nclass=as.numeric(nclass),
                   ABIC=as.numeric(ABIC)) %>% 
            select(nclass,ABIC)
          ,by="nclass") %>% 
  mutate_all(round,2) %>% 
  pivot_longer(2:4,values_to = "value",names_to = "metrics")

gg_abic <- ggplot(subset(dat_metrics,metrics=="ABIC"),aes(nclass,value))+
  geom_line(linewidth=1)+
  geom_point()+
  scale_y_continuous("ABIC",limits = c(63500,64100))+
  scale_x_continuous("Number of multimorbidity patterns")+
  ggprism::theme_prism(base_size = 16)

gg_acc <- ggplot(subset(dat_metrics,metrics=="Assignment accuracy (%)"),aes(nclass,value))+
  geom_line(linewidth=1)+
  geom_point()+
  scale_y_continuous("Assignment accuracy (%)",limits = c(0.25,1))+
  scale_x_continuous("Number of multimorbidity patterns")+
  ggprism::theme_prism(base_size = 16)

gg_en <- ggplot(subset(dat_metrics,metrics=="Entropy"),aes(nclass,value))+
  geom_line(linewidth=1)+
  geom_point()+
  scale_y_continuous("Entropy",limits = c(0.25,1))+
  scale_x_continuous("Number of multimorbidity patterns")+
  ggprism::theme_prism(base_size = 16)


ggmet <- ggpubr::ggarrange(gg_abic,
                           gg_acc,
                           gg_en,
                           ncol = 3)

ggmet


# After inspection of the metrics, we choose 5 classes. 
# Table for the interpretation of the patterns 

OE_ex_table(res$obj[[4]],5)

# ----- Multiple imputation ----

dat_imp=my_multiple_imputation(res5$obj,db.def,1000) 


dat_imp$data=lapply(dat_imp$data,function(x){
  x$MP[x$mm_w1=="No MM"]=6
  return(x)
})



dat_imp$data=lapply(dat_imp$data,function(x){
  x$MP2=case_when(x$MP== 1 ~ "Neuropsychiatric", x$MP==2 ~ "Psychiatric/Respiratory", x$MP== 3 ~ "Unspecific", 
                  x$MP== 4 ~ "Sensory impairment/Anaemia", x$MP== 5 ~ "Cardiometabolic", x$MP== 6 ~ "No MM")
  return(x)})


# ---->  Run LASSO ####

x <- data.matrix(db.def[, c("z_hemoglobin", "z_leukocytes", "z_HbA1C", "z_alkaline_phospatase", "z_gamma_gt", "z_albumin", 
                             "z_creatinine", "z_calcium", "z_vitamin_b12", "z_folic_acid", "z_cholesterol", "z_TSH", "z_T4", "z_beta2Microglobulin", "z_CCL2",
                             "z_TNFRSF1B", "z_Adiponectin", "z_CRP_luminex", "z_BDNF", "z_NCadherin", "z_MPO", "z_CystatinC", "z_Insulin", "z_CPeptide", "z_CXCL10",
                             "z_S100B", "z_EGF", "z_VEGF", "z_CCL3", "z_CCL4", "z_IGFBP1", "z_ESelectin", "z_PSelectin", "z_Leptin", "z_EphA2", "z_VCAM1",
                             "z_ICAM1", "z_alphaSynuclein", "z_GDF15", "z_CCL11", "z_hIFNg", "z_hIL10", "z_hIL12p70", "z_hIL1b", "z_il6", "z_hIL8",
                             "z_hTNFa", "z_ttau", "z_abeta40", "z_abeta42", "z_abeta_ratio","z_ptau181", "z_nfl", "z_gfap", "age_w1_s", "sex", "edu")])

grid =exp(seq (10 , -6 , length =100))

lasso_mult=lapply(1:length(dat_imp$data),function(y){
  suppressMessages({res=fun_lasso_multinomial(x,
                                              dat_imp$data[[y]]$MP2,
                                              grid,
                                              title="",
                                              varnobiom = c("age_w1_s", "sex", "edu"), 
                                              seed = 1234, 
                                              ncv = 200)}) 
  print(y)
  return(res)
  
})


# ---->  Single dataset with all the 1000 models - adding a column to identify the model ####

combined_df <- map2_dfr(lasso_mult, seq_along(lasso_mult), ~{
  df <- .x$coefdata
  df$num_model <- .y
  return(df)
})

# ---> Add a column to know if a marker is SELECTED (!=0) or not.

combined_df$marker_included <- ifelse(combined_df$Coefficient_selam_ref !=0, 1, 0)



# -------> Identify biomarkers selected at least 70% of times in the 1000 LASSO models and calculate mean coefficients

combined_df <- combined_df %>%
group_by(MM_pattern, Biomarkers) %>%
  mutate(selection_count = sum(marker_included == 1), 
         selection_perc = sum(marker_included == 1)/1000*100) %>%
  ungroup()

df_selected <- combined_df %>% filter(selection_perc >= 70 & Coefficient_selam_ref != 0)

df_selected_compact <- df_selected %>%
  group_by(MM_pattern, Biomarkers) %>%
  summarise(mean_coef = mean(Coefficient_selam_ref))

# generate dataframe with LASSO-selected biomarkers 
df_selected_compact <- df_selected_compact %>% filter(Biomarkers != "age_w1_s" & Biomarkers != "sex") #remove covariates


####### ----------------- 2. PCA ----------------- ##########
# PCA of LASSO-selected biomarkers for each multimorbidity patter 
# Generate dataframes of PCA results 

# ------> 2.1 PCA - Cardiometabolic Pattern

pca_cardio <- db.def %>% filter(mm_w1 == "Cardiometabolic")

bbm_cardio <- df_selected_compact %>% filter(MM_pattern == "Cardiometabolic" & !is.na(mean_coef))

pca_cardio <- pca_cardio[, intersect(bbm_cardio$Biomarkers, colnames(pca_cardio))]
cardio_pca <- prcomp(pca_cardio, 
                     center = T,
                     scale. = T)

c <- get_pca_var(cardio_pca)

explained_variance_cardio <- round((cardio_pca$sdev^2 / sum(cardio_pca$sdev^2)) * 100,1)


colnames(c$cos2) <- paste0("PC", seq_along(colnames(c$cos2)),
                           " (", round(explained_variance_cardio, 1), "%)")

custom_colors <- colorRampPalette(c("white", "#bb0c00"))(100)


cardio1 <- corrplot(c$cos2[,1:4], is.corr= F, col = custom_colors, tl.col = "black", 
                    col.lim = c(0,1), cl.pos = 'r', cl.ratio = 0.2,
                    addgrid.col = "darkgray",
                    method = "color",
                    tl.cex = 0.8, 
                    cl.cex = 0.5)

cardio_heat_data <- cardio1$corrPos
cardio_heat_data$Pattern <- "Cardiometabolic"

# --------> 2.2 PCA - Unspecific pattern

pca_unspe <- db.def %>% filter(mm_w1 == "Unspecific")

bbm_unspe <- df_selected_compact %>% filter(MM_pattern == "Unspecific" & !is.na(mean_coef))

pca_unspe <- pca_unspe[, intersect(bbm_unspe$Biomarkers, colnames(pca_unspe))]

unspe_pca <- prcomp(pca_unspe, 
                    center = T,
                    scale. = T)

unspe <- get_pca_var(unspe_pca)

explained_variance_unspe <- round((unspe_pca$sdev^2 / sum(unspe_pca$sdev^2)) * 100,1)

colnames(unspe$cos2) <- paste0("PC", seq_along(colnames(unspe$cos2)),
                               " (", round(explained_variance_unspe, 1), "%)")

custom_colors <- colorRampPalette(c("white", "#bb0c00"))(100)

unspe1 <- corrplot(unspe$cos2[,1:4], is.corr= F, col = custom_colors, tl.col = "black", 
                   col.lim = c(0,1), cl.pos = 'r', cl.ratio = 0.2,
                   addgrid.col = "darkgray",
                   method = "color",
                   tl.cex = 0.8, 
                   cl.cex = 0.5)

unspe_heat_data <- unspe1$corrPos

unspe_heat_data$Pattern <- "Unspecific"

# --------> 2.3 PCA - Neuropsychiatric pattern
pca_neuro <- db.def %>% filter(mm_w1 == "Neuropsychiatric")

bbm_neuro <- df_selected_compact %>% filter(MM_pattern == "Neuropsychiatric" & !is.na(mean_coef))

pca_neuro <- pca_neuro[, intersect(bbm_neuro$Biomarkers, colnames(pca_neuro))]

neuro_pca <- prcomp(pca_neuro, 
                    center = T,
                    scale. = T)

neuro <- get_pca_var(neuro_pca)

explained_variance_neuro <- round((neuro_pca$sdev^2 / sum(neuro_pca$sdev^2)) * 100,1)


colnames(neuro$cos2) <- paste0("PC", seq_along(colnames(neuro$cos2)),
                               " (", round(explained_variance_neuro, 1), "%)")

custom_colors <- colorRampPalette(c("white", "#bb0c00"))(100)

neuro1 <- corrplot(neuro$cos2[,1:4], is.corr= F, col = custom_colors, tl.col = "black", 
                   col.lim = c(0,1), cl.pos = 'r', cl.ratio = 0.2,
                   addgrid.col = "darkgray",
                   method = "color",
                   tl.cex = 0.8, 
                   cl.cex = 0.5)

neuro_heat_data <- neuro1$corrPos

neuro_heat_data$Pattern <- "Neuropsychiatric"

# -----> 2.4 PCA - Psychiatric/Respiratory pattern 

pca_psych <- db.def %>% filter(mm_w1 == "Psychiatric/Respiratory")

bbm_psych <- df_selected_compact %>% filter(MM_pattern == "Psychiatric/Respiratory" & !is.na(mean_coef))

pca_psych <- pca_psych[, intersect(bbm_psych$Biomarkers, colnames(pca_psych))]

psych_pca <- prcomp(pca_psych, 
                    center = T,
                    scale. = T)

psych <- get_pca_var(psych_pca)

explained_variance_psych <- round((psych_pca$sdev^2 / sum(psych_pca$sdev^2)) * 100,1)

colnames(psych$cos2) <- paste0("PC", seq_along(colnames(psych$cos2)),
                               " (", round(explained_variance_psych, 1), "%)")

custom_colors <- colorRampPalette(c("white", "#bb0c00"))(100)

psych1 <- corrplot(psych$cos2[,1:4], is.corr= F, col = custom_colors, tl.col = "black", 
                   col.lim = c(0,1), cl.pos = 'r', cl.ratio = 0.2,
                   addgrid.col = "darkgray",
                   method = "color",
                   tl.cex = 0.8, 
                   cl.cex = 0.5)

psych_heat_data <- psych1$corrPos

psych_heat_data$Pattern <- "Psychiatric/Respiratory"

# --------> 2.5 PCA - Sensory impairment/Anemia pattern
pca_senso <- db.def %>% filter(mm_w1 == "Sensory impairment/Anaemia")

bbm_senso <- df_selected_compact %>% filter(MM_pattern == "Sensory impairment/Anaemia" & !is.na(mean_coef))

pca_senso <- pca_senso[, intersect(bbm_senso$Biomarkers, colnames(pca_senso))]

senso_pca <- prcomp(pca_senso, 
                    center = T,
                    scale. = T)

senso <- get_pca_var(senso_pca)

explained_variance_senso <- round((senso_pca$sdev^2 / sum(senso_pca$sdev^2)) * 100,1)

colnames(senso$cos2) <- paste0("PC", seq_along(colnames(senso$cos2)),
                               " (", round(explained_variance_senso, 1), "%)")

custom_colors <- colorRampPalette(c("white", "#bb0c00"))(100)

senso1 <- corrplot(senso$cos2[,1:4], is.corr= F, col = custom_colors, tl.col = "black", 
                   col.lim = c(0,1), cl.pos = 'r', cl.ratio = 0.2,
                   addgrid.col = "darkgray",
                   method = "color",
                   tl.cex = 0.8, 
                   cl.cex = 0.5)

senso_heat_data <- senso1$corrPos

senso_heat_data$Pattern <- "Sensory impairment/Anaemia"

