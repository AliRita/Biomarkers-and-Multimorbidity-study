predict_dis_acc_rate <- function(data,idvar="id",agevar="age_base"){
  
  # Creating the vector with betas from LASSO
  values <- c(-0.00011, 0.010924159, 0.002366482, -0.002030394, 0.006050176, 0.001128262, 0.003014325, 0.023950733, 0.10272)
  
  # Assigning names to the elements of the vector
  names <- c("Intercept", "z_HbA1C", "z_gamma_gt", "z_albumin", "z_CystatinC", "z_Insulin", "z_Leptin", "z_GDF15", agevar)
  names(values) <- names
  
  #check data
  
  if(!all(names[names!="Intercept"]%in%colnames(data))){
    stop("not all variables in data!")
  }
  if(!(idvar%in%colnames(data))){
    stop(paste(idvar,"not contained in data!"))
  }
  
  dat_m <- model.matrix(as.formula(paste0("~z_HbA1C+z_gamma_gt+z_albumin+z_CystatinC+z_Insulin+z_Leptin+z_GDF15+",agevar)),data=data)
  beta_pred <- as.numeric(values %*% t(dat_m))
  res <- data %>% select(any_of(idvar))
  res$beta_pred <- beta_pred
  return(res)
}




mse <- function(dat){
  mean(( dat$beta_crude- dat$beta_pred)^2)
}



plot_square_error<- function(dat,lab="BLSA"){
  
  q05=as.numeric(quantile((dat$beta_crude-dat$beta_pred)^2,probs=0.05))
  
  print(q05)
  options(scipen=999)
  ggbox<-ggplot(dat)+
    geom_boxplot(aes(y=(beta_crude-beta_pred)^2))+
    geom_hline(aes(yintercept=mse(dat)),col=2)+
    scale_x_continuous(lab,breaks = NULL)+
    scale_y_log10("Square errors",breaks=round(c(q05,mse(dat),max((dat$beta_crude-dat$beta_pred)^2)),4))+ 
    theme_prism()
  print(ggbox)
  ggsave(ggbox,filename="Boxplot se.jpeg", width = 20, height = 25, units = "cm", dpi = 300)
  return((dat$beta_crude-dat$beta_pred)^2)
  print("Plot saved in folder!")
}