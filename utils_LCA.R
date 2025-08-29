#' Assign subjects to the their most probable mm pattern according to LCA model
#' @description This function assign a subject to a latent mm pattern according to a LCA model and disease combinations at a specific point in time. The assignment is done according to the mode of the posterior membership probability given by the LCA model.
#' @param obj poLCA object
#' @param X matrix with chronic diseases variables (coded as 1:no and 2:yes) to use for the calculation.
#'
#' @return numeric vector with LCA class assignment for each row in X
#' @export
#'
#' @examples
assign_LCA <- function(obj, X) {
  if(ncol(X %>% dplyr::select(any_of(colnames(obj$y))))!=length(colnames(obj$y))) stop("Some diseases used for the LCA are missing from X")
  pClY <- poLCA::poLCA.posterior(obj, y = X %>% dplyr::select(any_of(colnames(obj$y))))
  ng <- ncol(pClY)
  pred <- as.numeric(apply(pClY, 1, function(x) which.max(x))) # MODE
  return(pred)
}


#' Function to get the entropy from LCA object
#'
#' @param obj fit object of class poLCA.
#' @param X matrix with chronic diseases variables (coded as 1:no and 2:yes) to use for the calculation.
#' @param ratio boolen indicating whether the ratio between the prior and the posterior entropy should be returned. Default to TRUE.
#'
#' @return numeric
#' @export
#'
#' @examples
get_entropy <- function(obj, X, ratio = T) {
  error_prior <- entropy(obj$P)
  posterior_probs <- poLCA::poLCA.posterior(obj, X)
  error_post <- mean(apply(posterior_probs, 1, entropy), na.rm = TRUE)
  entropy_value <- 1 - (error_post / error_prior)
  if (!ratio) {
    entropy_value <- 1 - error_post
  }
  return(entropy_value)
}



#' Entropy internal
#'
#' This function is used internally within the package.
#' @param p A numeric value.
#' @return Entropy.
#' @keywords internal


entropy <- function(p) sum(-p * log2(p))

#' Function to calculate the approximate misclassifications matrix
#'
#' @param fit object of class poLCA.
#' @param X matrix with chronic diseases variables (coded as 1:no and 2:yes) to use for the calculation.
#'
#' @details This function calculates the misclassification probabilities as in https://www.stata.com/meeting/uk22/slides/UK22_Tompsett.pdf
#' @return nc x nc matrix where nc is the number of latent multimorbidity patterns containing the probabilities of being assigned to a class given the true class.
#' @export
#'
#' @examples
get_internal_validation_matrix <- function(fit, X, norm = T, n_norm = NULL) {
  pClY <- poLCA::poLCA.posterior(fit, y = X)
  ng <- ncol(pClY)
  pred <- as.numeric(apply(pClY, 1, function(x) which.max(x))) # MODE
  
  
  
  Ptable <- cbind(pClY, pred)
  Pmatrix <- matrix(0, ng, ng)
  Npmatrix <- matrix(0, ng, ng)
  modclass <- pred
  
  
  for (i in 1:ng) {
    for (j in 1:ng) {
      Pmatrix[i, j] <- sum(subset(Ptable, modclass == i)[, j]) / table(modclass)[i]
      Npmatrix[i, j] <- Pmatrix[i, j] * table(modclass)[i]
    }
  }
  
  
  denom <- colSums(Npmatrix)
  Qmatrix <- matrix(0, ng, ng)
  
  
  for (i in 1:ng) {
    for (j in 1:ng) {
      Qmatrix[j, i] <- Npmatrix[i, j] / denom[j]
    }
  }
  
  
  
  if (norm) {
    Qmatrix * (nrow(X) / n_norm * 100)
  }
  Qmatrix
}




# Function to obtain OE/Ex for interpretation of the patterns
OE_ex_table <- function(obj,nclass,cutoff_OE=2,cutoff_Ex=0.25){
  E <-apply(obj$y-1,2,mean)
  n <- list()
  for (j in 1:nclass){
    n[[j]] <- apply(obj$y[obj$predclass==j,]-1,2,mean)
  }
  
  
  O <- do.call("cbind",n)
  R <- O/E
  
  O %<>% as.data.frame() %>% rownames_to_column("Disease") %>% pivot_longer(2:(nclass+1),
                                                                            names_to = "Multimorbidity profile",
                                                                            values_to = "Prevalence") %>% 
    mutate(`Multimorbidity profile`=as.numeric (gsub("\\D", "", `Multimorbidity profile`)))
  
  
  R %<>% as.data.frame() %>% rownames_to_column("Disease") %>% pivot_longer(2:(nclass+1),
                                                                            names_to = "Multimorbidity profile",
                                                                            values_to = "O/E") %>% 
    mutate(`Multimorbidity profile`=as.numeric (gsub("\\D", "", `Multimorbidity profile`)))%>% 
    mutate(label=ifelse(`O/E`<cutoff_OE,NA_integer_,Disease))
  
  
  ######## Exclusivity ########
  
  N <-apply(obj$y-1,2,sum)
  
  n <- list()
  for (j in 1:nclass){
    n[[j]] <- apply(obj$y[obj$predclass==j,]-1,2,sum)
  }
  
  Ex <- do.call("cbind",n)
  
  Ex <- Ex/N
  Ex %<>% as.data.frame() %>% rownames_to_column("Disease") %>% pivot_longer(2:(nclass+1),
                                                                             names_to = "Multimorbidity profile",
                                                                             values_to = "Exclusivity") %>% 
    mutate(`Multimorbidity profile`=as.numeric (gsub("\\D", "", `Multimorbidity profile`)))%>% 
    mutate(label2=ifelse(`Exclusivity`<cutoff_Ex,NA_integer_,Disease))
  
  
  Char_MP <- R %>% left_join(Ex) %>% left_join(O)
  Char_MP %<>% mutate(char=ifelse(!is.na(label) & !is.na(label2),1,NA_integer_)) 
  
  
  Char_MP %<>% mutate(`Multimorbidity profile`=as.factor(`Multimorbidity profile`)) %>% 
    dplyr::select(`Multimorbidity profile`,Disease,Prevalence,Exclusivity,`O/E`,char) %>% 
    mutate(Prevalence=Prevalence*100,
           Exclusivity=Exclusivity*100,
           Disease=gsub("_"," ",Disease),
           `O/E`=round(`O/E`,2)) %>% 
    arrange(`Multimorbidity profile`,Disease) %>% 
    mutate_at(3:4,round,0) 
  
  
  
  return(Char_MP)
}

#' Function to obtain imputate multiple imputations of the assigned MM patterns
#'
#' @param obj LCA object
#' @param data your dataset
#' @param nimp number of imputated datasets
#'
#' @return list of nimp dataset with assigned MM pattern according to posterior distribution probabilities
#' @export
#'
#' @examples
multiple_imputation <- function(obj, data, nimp) {
  X <- data %>% dplyr::select(dplyr::any_of(colnames(obj$y)))
  if (ncol(X) != length(colnames(obj$y))) {
    stop("Some disease columns are missing!")
  }
  lev <- unique(unlist(lapply(1:ncol(X), function(y) {
    unique(X[
      ,
      y
    ])
  })))
  if (NA %in% lev) {
    stop("Disease columns cannot contain missing values.")
  }
  if (!all(lev %in% c(0, 1)) | any(sapply(lev, is.character))) {
    stop("Disease columns must be numeric 0/1 variables.")
  }
  X %<>% dplyr::mutate_all(function(x) x + 1)
  data_imp <- list()
  pClY <- poLCA::poLCA.posterior(obj, y = X)
  ng <- ncol(pClY)
  for (i in 1:nimp) {
    data_imp[[i]] <- data
    data_imp[[i]]$MP <- unlist(lapply(1:nrow(data), function(x) sample(1:ng, size = 1, replace = TRUE, prob = pClY[x, ])))
  }
  return(data_imp)
}


#' Main function to run the LCA
#' @description  This is a helper function to run a  latent class model using poLCA. For  details regarding poLCA see the package documentation and https://statistics.ohlsen-web.de/latent-class-analysis-polca/.
#' @param nclass Numeric indicating the number of latent classes.
#' @param X Matrix with chronic diseases variables (coded as 1:no and 2:yes) to use for the calculation.
#' @param conditions Vector of columns names indicating the conditions to use for the LCA. It can be the object returned from the function select_conditions.
#' @param nrep Number of times to estimate the model, using different values of probs.start. The default is one. Setting nrep>1 automates the search for the global—rather than just a local—maximum of the log-likelihood function. poLCA returns the parameter estimates corresponding to the model with the greatest log-likelihood. Default to 50. Reduce to save computation time.
#' @param fcov Covariates used as predictors.
#' @param probstart A list of matrices of class-conditional response probabilities to be used as the starting values for the estimation algorithm. Each matrix in the list corresponds to one manifest variable, with one row for each latent class, and one column for each outcome. The default is NULL, producing random starting values. Note that if nrep>1, then any user-specified probs.start values are only used in the first of the nrep attempts. Deafult to NULL.
#'
#' @return A list containing the following elements:
#' \item{obj}{An object of class 'poLCA' representing the fitted LCA model.}
#' \item{accuracy_matrix}{A matrix containing the mislassification probabilities.}
#' \item{metrics}{A vector containing goodness-of-fit measures.}
#' @export
#'
#' @examples
run_LCA <- function(nclass, X, conditions, nrep = 50, fcov = NULL, probstart = NULL) {
  gc()
  f <- X %>%
    dplyr::select(any_of(conditions)) %>%
    as.matrix()
  if (!is.null(fcov)) {
    formula <- as.formula(paste0("f~", fcov))
  } else {
    formula <- as.formula(paste0("f~1"))
  }
  
  if (is.null(probstart)) {
    myresult <- poLCA::poLCA(formula,
                             X,
                             nclass = nclass,
                             maxiter = 3000,
                             na.rm = F,
                             nrep = nrep,
                             verbose = T
    )
  } else {
    myresult <- poLCA::poLCA(formula,
                             X,
                             nclass = nclass,
                             maxiter = 3000,
                             na.rm = F,
                             probs.start = probstart,
                             verbose = T
    )
  }
  
  print(paste("Completed LCA with", nclass, "..."))
  
  Modell <- paste("Model", nclass)
  log_likelihood <- myresult$llik
  df <- myresult$resid.df
  BIC <- myresult$bic
  AIC <- myresult$aic
  ABIC <- (-2 * myresult$llik) + ((log((myresult$N + 2) / 24)) * myresult$npar)
  CAIC <- (-2 * myresult$llik) + myresult$npar * (1 + log(myresult$N))
  likelihood_ratio <- myresult$Gsq
  internal_val_matrix <- get_internal_validation_matrix(myresult, X)
  acc <- sum(diag(internal_val_matrix)) / nclass
  entropy <- get_entropy(myresult, X)
  return(list(
    obj = myresult,
    accuracy_matrix = internal_val_matrix,
    metrics = c(Modell,
                nclass,
                log_likelihood,
                df,
                BIC,
                AIC,
                CAIC,
                ABIC,
                likelihood_ratio,
                acc = acc,
                entropy = entropy
    )
  ))
}


run_lca_cv <- function(fold, nclasses, X, conditions, nrep, folds) {
  res <- lapply(nclasses, function(x) run_LCA(x, X = X[-(folds == fold), ], conditions = conditions, nrep = nrep))
  objects <- lapply(res, function(x) x$obj)
  names(objects) <- nclasses
  entropy <- unlist(lapply(1:length(objects), function(x) get_entropy(objects[[x]], X[(folds == fold), ], ratio = T)))
  names(entropy) <- nclasses
  ext_accuracy <- unlist(lapply(1:length(objects), function(x) sum(diag(get_internal_validation_matrix(objects[[x]], X[(folds == fold), ]))) / nclasses[x]))
  dat_res <- data.frame(
    cv = fold,
    nclass = nclasses,
    acc = ext_accuracy,
    entropy = entropy
  )
  
  
  colnames(dat_res) <- c(
    "CV",
    "nclass",
    "Assignment accuracy (%)",
    "Entropy"
  )
  return(dat_res)
}


#' Function to run LCA with different number of classes and compare the goodness-of-fit
#' @description  This is a helper function to select the number of latent classes using poLCA. For  details regarding poLCA see the package documentation and https://statistics.ohlsen-web.de/latent-class-analysis-polca/.
#' @param nclasses Numeric vector indicating the number of latent classes to investigate.
#' @param X Matrix with chronic diseases variables (coded as 1:no and 2:yes) to use for the calculation.
#' @param conditions Vector of columns names indicating the conditions to use for the LCA. It can be the object returned from the function select_conditions.
#' @param plot Boolean indicating whether the goodness-of-fit measures should be plotted.
#' @param nrep Number of times to estimate the model, using different values of probs.start. The default is one. Setting nrep>1 automates the search for the global—rather than just a local—maximum of the log-likelihood function. poLCA returns the parameter estimates corresponding to the model with the greatest log-likelihood. Default to 50. Reduce to save computation time.
#' @param cvfolds Number of folds of the cross-validation, default to NULL (no CV is performed)
#' @return A list containing the following elements:
#' \item{metrics}{A dataframe containing the goodness-of-fit measures for the different models.}
#' \item{obj}{A list of poLCA objects.}
#' \item{plot}{ggplot object for the plot if plot=T}
#' \item{accuracy_matrix}{A list of matrices containing the mislassification probabilities.}
#' \item{elapsed_time}{Numeric indicating the time elapsed.}
#' @export
#'
#' @examples
select_number_LCA <- function(nclasses, X, conditions, plot = T, nrep = 50, cvfolds = NULL) {
  tictoc::tic()
  if (!is.null(cvfolds)) {
    future::plan(future::multisession)
    folds <- sample(rep(1:cvfolds, length.out = nrow(X)))
    results <- do.call("rbind", future.apply::future_lapply(1:cvfolds, FUN = run_lca_cv, nclasses, X, conditions, nrep, folds))
    results_sum <- results %>%
      dplyr::group_by(nclass) %>%
      dplyr::summarise_all(mean, na.rm = T) %>%
      dplyr::select(-CV)
    
    return(list(results_sum, results))
    if (plot) {
      dat_res_wide <- results %>%
        as.data.frame() %>%
        tidyr::pivot_longer(2:3, values_to = "metrics", names_to = "name") %>%
        dplyr::mutate(
          metrics = as.numeric(metrics),
          nclass = as.numeric(nclass)
        )
      
      gg <- ggplot2::ggplot(dat_res_wide) +
        ggplot2::geom_line(ggplot2::aes(nclass, metrics)) +
        ggplot2::geom_point(ggplot2::aes(nclass, metrics)) +
        ggplot2::facet_wrap(~name, scales = "free_y") +
        ggplot2::scale_y_continuous("") +
        ggplot2::scale_x_continuous("Number of latent classes", breaks = nclasses) +
        ggplot2::theme_bw()
      
      print(gg)
    } else {
      gg <- NULL
    }
    
    return(list(plot = gg, results = results_sum, results_sum = results_sum))
  } else {
    res <- lapply(nclasses, function(x) run_LCA(x, X = X, conditions = conditions, nrep = nrep))
    dat_res <- do.call("rbind", lapply(res, function(x) x$metrics))
    objects <- lapply(res, function(x) x$obj)
    names(objects) <- nclasses
    internal_val <- lapply(res, function(x) x$accuracy_matrix)
    names(internal_val) <- nclasses
    colnames(dat_res) <- c(
      "Model",
      "nclass",
      "log_likelihood",
      "df",
      "BIC",
      "AIC",
      "CAIC",
      "ABIC",
      "likelihood_ratio",
      "Assignment accuracy (%)",
      "Entropy"
    )
    
    if (plot) {
      dat_res_wide <- dat_res %>%
        as.data.frame() %>%
        tidyr::pivot_longer(5:11, values_to = "metrics", names_to = "name") %>%
        dplyr::mutate(
          metrics = as.numeric(metrics),
          nclass = as.numeric(nclass)
        )
      
      
      
      
      
      gg <- ggplot2::ggplot(dat_res_wide) +
        ggplot2::geom_line(ggplot2::aes(nclass, metrics)) +
        ggplot2::geom_point(ggplot2::aes(nclass, metrics)) +
        ggplot2::facet_wrap(~name, scales = "free_y") +
        ggplot2::scale_y_continuous("") +
        ggplot2::scale_x_continuous("Number of latent classes", breaks = nclasses) +
        ggplot2::theme_bw()
      
      print(gg)
    } else {
      gg <- NULL
    }
  }
  elapsed_time <- tictoc::toc(quiet = TRUE)$toc - tictoc::toc(quiet = TRUE)$tic
  return(list(metrics = dat_res, obj = objects, plot = gg, accuracy_matrix = internal_val, elapsed_time = elapsed_time))
}

