#------------------------------------------------------------------------------------------------------------------------------------------------#
# Analyses for "Shared and specific blood biomarkers for multimorbidity"
# 1. Linear Mixed-effect model (LMM) to estimate the rate of disease accumulation
# 2. Gaussian LASSO 
# 3. PCA
# R version 4.4.2
#------------------------------------------------------------------------------------------------------------------------------------------------#

# load packages
library(glmnet)
library(tidyverse)
library(factoextra)
library(corrplot)

# load dataset
load("db.long.RData") #long format of study dataset

####### ------- 1. Linear Mixed-effect Models -------- ######

model0.lmer <- lmer(disease_num ~ wave + (wave|ID), data=db.long, REML=T)
db.long.def$pred <- fitted(model0.lmer) 
beta_crude <-ranef(model0.lmer)$ID$`wave`


###### ------- 2. multinomial LASSO -------- #######
set.seed(1234)

cv.out.x = cv.glmnet(x, beta_crude, alpha =1, family="gaussian",lambda = grid, type.measure = "mse", nfolds = 200) # Note: x = matrix of biomarkers and covariates (see script 01)
selam1 = cv.out.x$lambda.1se
fit.x <- glmnet(x, beta_crude, family="gaussian", lambda = grid)

lasso_coefficients_x = round(coef(fit.x, s=selam1), 5)
lasso_coefficients_x <- as.data.frame(as.matrix(lasso_coefficients_x))

lasso_coefficients_x <- lasso_coefficients_x[-1, , drop = FALSE] #drop intercept

lasso_coefficients_x$biomarker <- rownames(lasso_coefficients_x)

colnames(lasso_coefficients_x) <- c("Coefficient", "Biomarkers")

lasso_coefficients_x <- lasso_coefficients_x[1:54,] # removing age

lasso_non0_long <- lasso_coefficients_x %>% filter(Coefficient != 0) # Keeping only non zero


######## ------------------- 3. PCA -------- ##############
biom_long_x <- lasso_non0_long$Biomarkers
#selezione manuale dato che ho rinominato le variabili -.-"
pca_long <- db[, c(92,94,95,111,112,123,128)]

long_pca <- prcomp(pca_long, 
                   center = T,
                   scale. = T)

l <- get_pca_var(long_pca)

explained_variance <- round((long_pca$sdev^2 / sum(long_pca$sdev^2)) * 100,1)

colnames(l$cos2) <- paste0("PC", seq_along(colnames(l$cos2)),
                           " (", round(explained_variance, 1), "%)")

custom_colors <- colorRampPalette(c("white", "#bb0c00"))(100)

long_plot <- corrplot(l$cos2[,1:5], is.corr= F, col = custom_colors, tl.col = "black", 
                      col.lim = c(0,1), cl.pos = 'r', cl.ratio = 0.25,  
                      addgrid.col = "darkgray",
                      method = "color",
                      tl.cex = 0.9, 
                      cl.cex = 0.55)

heat_data <- long_plot$corrPos
heat_data$Pattern <- "Speed of disease accumulation"

