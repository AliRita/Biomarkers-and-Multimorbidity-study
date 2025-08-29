
#------------------------------------------------------------------------------------------------------------------------------------------------#
# Analyses for "Shared and specific blood biomarkers for multimorbidity"
# 1. Gaussian LASSO for number of diseases at baseline
# 2. PCA of LASSO-selected biomarkers for number of diseases at baseline
# R version 4.4.2
#------------------------------------------------------------------------------------------------------------------------------------------------#

# load packages
library(glmnet)
library(tidyverse)
library(factoextra)
library(corrplot)

#load dataset
load("db.RData") #clean study dataset containing all necessary data


# ------------------------------------------------------------------------------------------------------------------------------------------------#
# (1) Define parameters for gaussian LASSO, run the model, and generate non-zero coefficients dataframe
# (2) PCA of LASSO-selected biomarkers
# (3) Generate dataframe with PCA data for plotting purpose
# ------------------------------------------------------------------------------------------------------------------------------------------------#


#### --------- 1. Gaussian LASSO ---------- ####

# 1.1 Set seed
set.seed(1234)

y <- db$chron_num_w1 #vector of number of diseases at baseline

# 1.2. Create the matrix with biomarkers and covariates

x <- data.matrix(db[, c("Hb", "WBC", "HbA1c", "ALP", "GGT", "ALB", "Cr", "Ca", "VitB12", "FA", 
                                    "TC", "TSH", "T4", "B2M", "CCL2", "TNFRSF1B", "Adipo", "CRP", "BDNF", "NCAD", 
                                    "MPO", "CysC", "Insulin", "CPEP", "CXCL10", "S100B", "EGF", "VEGF", "CCL3", "CCL4", 
                                    "IGFBP1", "ESel", "PSel", "LEP", "EphA2", "VCAM1", "ICAM1", "aSyn", "GDF15", "CCL11", 
                                    "IFNG", "IL10", "IL12p70", "IL1b", "IL6", "IL8", "TNFa", "ttau", "AB40", "AB42", 
                                    "ABratio", "ptau181", "NfL", "GFAP", "age_w1_s", "sex", "edu")]) #age is standardized

grid =exp(seq (10 , -6 , length =100))

cv.out.x = cv.glmnet(x, y, alpha =1, family="gaussian", lambda = grid, type.measure = "mse", nfolds = 200)
selam1 = cv.out.x$lambda.1se
fit.x <- glmnet(x, x, family="gaussian", lambda = grid)

lasso_coefficients_x = round(coef(fit.x, s=selam1), 5)
lasso_coefficients_x <- as.data.frame(as.matrix(lasso_coefficients_x))

lasso_coefficients_x <- lasso_coefficients_x[-1, , drop = FALSE] # Drop intercept

lasso_coefficients_x$biomarker <- rownames(lasso_coefficients_x) # Add a column with variable name

colnames(lasso_coefficients_x) <- c("Coefficient", "Biomarkers")

lasso_coefficients_x <- lasso_coefficients_x[1:54,] # Remove covariates (i.e., age, sex, education)

lasso_non0 <- lasso_coefficients_x %>% filter(Coefficient != 0) # Keep only non zero coefficients


### ------- 2. PCA ------- ####

# ------> 2.1 Select LASSO-selected biomarkers
biom_cross <- lasso_non0$Biomarkers

pca_dis <- db[, intersect(bbm_cross$Biomarkers, colnames(db))]

# pca_dis <- db.def[, c(128,111,92,96, 112, 113, 123,142,90)]

dis_pca <- prcomp(pca_dis, 
                  center = T,
                  scale. = T)
dis <- get_pca_var(dis_pca)

explained_variance <- round((dis_pca$sdev^2 / sum(dis_pca$sdev^2)) * 100,1)

colnames(dis$cos2) <- paste0("PC", seq_along(colnames(dis$cos2)),
                             " (", round(explained_variance, 1), "%)")

custom_colors <- colorRampPalette(c("white", "#bb0c00"))(100)

dis_plot <- corrplot(dis$cos2[,1:4], is.corr= F, col = custom_colors, tl.col = "black", 
                     col.lim = c(0,1), cl.pos = 'r', cl.ratio = 0.25,  
                     addgrid.col = "darkgray",
                     method = "color",
                     tl.cex = 0.9, 
                     cl.cex = 0.55)

heat_dis_data <- dis_plot$corrPos

heat_dis_data$Pattern <- "Number of chronic diseases" # Label to plot together with other multimorbidity measures.
