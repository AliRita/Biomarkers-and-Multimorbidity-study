fun_lasso_multinomial <- function(x,y,grid,title,varnobiom, ncv, seed, axisylimits = NULL){
  
  set.seed(seed)
  cv.out = cv.glmnet(x, y, alpha =1, family="multinomial",lambda = grid,type.measure="deviance", nfolds = ncv)
  plot(cv.out)
  
  print(cv.out)

  print(paste("value of lambda that gives minumum cvm:",  cv.out$lambda.min))
  print(paste("largest value of lambda such that error is within 1 standard error of the minimum:",  cv.out$lambda.1se))
  print(paste("Cross-validated accuracy lambda1se:",  cv.out$cvm[cv.out$lambda == cv.out$lambda.1se]))
  print(paste("Cross-validated accuracy min lambda:",  cv.out$cvm[cv.out$lambda == cv.out$lambda.min]))
  print(paste("number of non-zero coefficients at lambda1se:",  cv.out$nzero[cv.out$lambda == cv.out$lambda.1se]))
  print(paste("number of non-zero coefficients at min lambda:",  cv.out$nzero[cv.out$lambda == cv.out$lambda.min]))

  
  #lambda min
  bestlam = cv.out$lambda.min
  
  #lambda 1 se
  selam = cv.out$lambda.1se
 
  

  plot(cv.out$glmnet.fit, 
       xvar = "lambda", label=FALSE)
  
  fit <- glmnet(x, y, family="multinomial", alpha = 1,type.measure="deviance", lambda = grid) 
  plot(fit)
  coef_list <- coef(fit, s=bestlam) 
  
  
  multinom_coefs_bestlam <- coef(fit, s=bestlam) 
  
  rounded_multinom_coefs_bestlam <- lapply(multinom_coefs_bestlam, function(mat) round(mat, 10))
  
  print(paste("Cross-validated accuracy lambda1se (relative):",  fit$dev.ratio[fit$lambda == cv.out$lambda.1se]))
  print(paste("Cross-validated accuracy min lambda (relative):",  fit$dev.ratio[fit$lambda == cv.out$lambda.min]))
  
  # Combine the coefficients into a data frame
  multinom_coef_bestlam_list <- lapply(names(rounded_multinom_coefs_bestlam), function(class) {
    coefs_df <- as.data.frame(as.matrix(rounded_multinom_coefs_bestlam[[class]]))
    coefs_df$Coefficient_bestlam<- rownames(rounded_multinom_coefs_bestlam)
    coefs_df$MM_pattern <- class  # Add class label
    coefs_df$Biomarkers <- rownames(coefs_df) 
    coefs_df
  })
  
  # Combine all classes into one data frame and remove intecepts
  coefs_multinom_bestlam_db <- do.call(rbind, multinom_coef_bestlam_list)
  colnames(coefs_multinom_bestlam_db)[1] <- "Coefficient_bestlam"
  coefs_multinom_bestlam_db_no_interc <- coefs_multinom_bestlam_db %>% filter(Biomarkers != "(Intercept)")
  
  
  # ------------> LAMBDA 1se coefficients -------------
  multinom_coefs_selam <- coef(fit, s=selam) 
  
  rounded_multinom_coefs_selam <- lapply(multinom_coefs_selam, function(mat) round(mat, 10))
  
  # Combine the coefficients into a data frame
  multinom_coef_selam_list <- lapply(names(rounded_multinom_coefs_selam), function(class) {
    coefs_df <- as.data.frame(as.matrix(rounded_multinom_coefs_selam[[class]]))
    coefs_df$Coefficient_bestlam<- rownames(rounded_multinom_coefs_bestlam)
    coefs_df$MM_pattern <- class  # Add class label
    coefs_df$Biomarkers <- rownames(coefs_df) 
    coefs_df
  })
  
  # Combine all classes into one data frame
  coefs_multinom_selam_db <- do.call(rbind, multinom_coef_selam_list)
  colnames(coefs_multinom_selam_db)[1] <-"Coefficient_selam"
  coefs_multinom_selam_db_no_interc <- coefs_multinom_selam_db %>% filter(Biomarkers != "(Intercept)")
  
  # Merging both lambda min e lambda 1 se coefficients 
  coefs_multinomial_all_no_interc <- merge(coefs_multinom_bestlam_db_no_interc, coefs_multinom_selam_db_no_interc, by = c("MM_pattern", "Biomarkers"))

  coefs_multinomial_all_no_interc_ref <-  
    coefs_multinomial_all_no_interc%>% 
    ungroup() %>% 
    mutate(is_nomm_bestlam=ifelse(MM_pattern=="No MM",Coefficient_bestlam,NA_real_),
           is_nomm_selam=ifelse(MM_pattern=="No MM",Coefficient_selam,NA_real_)) %>% 
    group_by(Biomarkers) %>% 
    arrange(Biomarkers) %>% 
    mutate(nomm_bestlam=min(is_nomm_bestlam,na.rm = T),
           nomm_selam=min(is_nomm_selam,na.rm = T))%>% 
    ungroup() %>% 
    mutate(Coefficient_bestlam_ref=Coefficient_bestlam-nomm_bestlam,
           Coefficient_selam_ref=Coefficient_selam-nomm_selam) %>% 
    dplyr::select(-is_nomm_selam,
           -is_nomm_bestlam) %>% 
    filter(MM_pattern!="No MM")
  
  print("Variables in argument varnobiom are removed from plot!")
  
  if(is.null(axisylimits)) axisylimits <- c(-1,1)
  
  gg1 <- coefs_multinomial_all_no_interc_ref %>%
    filter(Coefficient_selam_ref != 0 & !(Biomarkers %in% varnobiom)) %>%
    ggplot(aes(x = Biomarkers, y = Coefficient_selam_ref)) +
    geom_bar(stat = "identity", fill = "#2c61ab") + 
    # scale_y_continuous(limits = axisylimits) + 
    coord_flip() +  # Flip coordinates to make it easier to read variable names
    labs(title = title,
         x = "Biomarkers",
         y = "Coefficients") +
    #scale_y_log10()+
    facet_grid(.~ MM_pattern)+
    theme_prism()+
    theme(axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13),
          strip.text = element_text(size = 15, family = "Arial", face = "bold"))
  
  
  
  gg2 <- coefs_multinomial_all_no_interc_ref %>%
    filter(Coefficient_bestlam_ref != 0 & !(Biomarkers %in% varnobiom)) %>%
    ggplot(aes(x = Biomarkers, y = Coefficient_bestlam_ref)) +
    geom_bar(stat = "identity", fill = "darkblue") + 
    #scale_y_continuous(limits = axisylimits) +
    coord_flip() +  # Flip coordinates to make it easier to read variable names
    labs(title = title,
         x = "Biomarkers",
         y = "Coefficients") +
    #scale_y_log10()+
    facet_grid(.~ MM_pattern)+
    theme_prism()+
    theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))
  
  
  
  return(list(fit=fit,
              fitcv=cv.out,
         lambda1se=selam,
         lambdamin=bestlam,
         coefdata=coefs_multinomial_all_no_interc_ref,
         plotmin=gg2,
         plot1se=gg1))
  
  
}
