# Transportation

# Doubly robust estimation
transport_dr <- function(a, y, x_trial, x_target, sl.lib.source, sl.lib.target, n_folds, mu1_hat = NULL, mu0_hat = NULL, pi_hat = NULL, rho_hat = NULL, tau1_hat = NULL, tau0_hat = NULL, epsilon = NULL, survey = FALSE){
  # a: treatment vector
  # y: outcome vector
  # x_trial: covariates from the source dataset
  # x_target: covariates from the target dataset, should be a subset of x_trial
  # x_trial and x_target will be matched according to the colnames
  # sl.lib.source: superlearner to train models in the source dataset
  # sl.lib.target: superlearner to train models in the target dataset
  # n_folds: used for cross-fitting
  # mu1_hat: the estimator for mu1, usually unknown and needs to be estimated
  # Similar for other nuisance estimators
  # epsilon: the positivity constant enforced (i.e. all participation probability and propensity score will be enforced in [epsilon, 1-epsilon])
  # survey: whether survey data is used for target dataset 
  # return: results include three columns corresponding to E[Y(0)], E[Y(1)] and E[Y(1)]-E[Y(0)]. The first row includes point estimates and second row includes estimated standard error.
  #         ifs include all the influence function values for each individual
  # example usage in our real data analysis: pre_fruit <- na.omit(numom[,c(cova_fruit, "fruit80", "ptb37")])
  # pre_fruit_tran <- transport_dr(pre_fruit$fruit80, pre_fruit$ptb37, pre_fruit[,cova_fruit], nsfg, sl.lib.source = c( "SL.ranger","SL.mean","SL.glmnet"), sl.lib.target = c( "SL.ranger","SL.mean","SL.glmnet"), 5, epsilon = 0.01, survey = T)
  require("SuperLearner")
  n_trial <- nrow(x_trial)
  n_target <- nrow(x_target)
  muhat <- matrix(0, nrow=n_trial, ncol=2)
  pi_trial <- numeric(length=n_trial)
  pi_treat <- numeric(length=n_trial)
  tau_trial <- matrix(0, nrow=n_trial, ncol=2)
  tau_target <- matrix(0, nrow=n_target, ncol=2)
  set.seed(123, "L'Ecuyer-CMRG")
  set.seed(521)
  split_trial <- sample(1:n_folds, n_trial, replace = T)
  split_target <- sample(1:n_folds, n_target, replace = T)
  x_trial_aug <- janitor::clean_names(as.data.frame(model.matrix(~-1+., x_trial)))
  x_target_aug <- janitor::clean_names(as.data.frame(model.matrix(~-1+., x_target)))
  cova_target_aug <- colnames(x_target_aug)
  
  # estimate regression function
  if (is.null(mu1_hat)){
    for (aa in 1:2){
      for (i in 1:n_folds){
        train_trial <- split_trial != i
        train_target <- split_target != i
        test_trial <- split_trial == i
        test_target <- split_target == i
        mufit <- mcSuperLearner(y[a==aa-1 & train_trial],
                                x_trial_aug[a==aa-1 & train_trial,],
                                newX=x_trial_aug, SL.library=sl.lib.source)
        muhat[test_trial,aa] <- mufit$SL.predict[test_trial]
        mufit_val <- mufit$SL.predict[train_trial]
        taufit <- mcSuperLearner(mufit_val,
                                 x_trial_aug[train_trial, cova_target_aug],
                                 newX=rbind(x_trial_aug[test_trial, cova_target_aug], x_target_aug[test_target,]), SL.library=sl.lib.target)
        tau_trial[test_trial,aa] <- taufit$SL.predict[1:nrow(x_trial_aug[test_trial, cova_target_aug])]
        tau_target[test_target,aa] <- taufit$SL.predict[(1+nrow(x_trial_aug[test_trial, cova_target_aug])):(length(taufit$SL.predict))]
      }
    }
  }else{
    muhat[,1] <- mu0_hat(x_trial)
    muhat[,2] <- mu1_hat(x_trial)
    tau_trial[,1] <- tau0_hat(x_trial)
    tau_trial[,2] <- tau1_hat(x_trial)
    tau_target[,1] <- tau0_hat(cbind(x_target, matrix(0, nrow = n_target, ncol=2)))
    tau_target[,2] <- tau1_hat(cbind(x_target, matrix(0, nrow = n_target, ncol=2)))
  }
  
  # estimate propensity score into trial
  if (is.null(rho_hat)){
    data_aug <- rbind(x_trial_aug[,cova_target_aug], x_target_aug)
    label <- c(rep(1, n_trial), rep(0, n_target))
    for (i in 1:n_folds){
      train_trial <- split_trial != i
      train_target <- split_target != i
      test_trial <- split_trial == i
      test_target <- split_target == i
      pifit1 <- mcSuperLearner(label[c(train_trial, train_target)],data_aug[c(train_trial, train_target),],
                               newX=x_trial_aug[test_trial,cova_target_aug], SL.library=sl.lib.target, family=binomial)
      pi_trial[test_trial] <-pifit1$SL.predict
    }
    # print("rho")
    # print(sum(pi_trial<epsilon))
    # print(sum(pi_trial>1-epsilon))
  }else{
    pi_trial <- rho_hat(x_trial)
  }
  
  # estimate propensity score of treatment in trial data
  if (is.null(pi_hat)){
    for (i in 1:n_folds){
      train_trial <- split_trial != i
      test_trial <- split_trial == i
      pifit2 <- mcSuperLearner(a[train_trial],x_trial_aug[train_trial,],
                               newX=x_trial_aug[test_trial,], SL.library=sl.lib.source, family=binomial)
      pi_treat[test_trial] <-pifit2$SL.predict
    }
    # print("pi")
    # print(sum(pi_treat<epsilon ))
    # print(sum(pi_treat>1-epsilon ))
  }else{
    pi_treat <- pi_hat(x_trial)
  }
  
  # Enforce strict positivity
  if(!is.null(epsilon)){
    pi_treat[pi_treat < epsilon] <- epsilon
    pi_treat[pi_treat > 1-epsilon] <- 1-epsilon
    pi_trial[pi_trial < epsilon] <- epsilon
    pi_trial[pi_trial > 1-epsilon] <- 1-epsilon
  }
  
  # Compute influence function for each data point
  ifs <- matrix(0, nrow=n_trial + n_target, ncol=2)
  for (aa in 1:2){
    ifs[1:n_trial,aa] <- as.numeric(a == aa-1)*(1-pi_trial)*(y-muhat[,aa])/(pi_trial*(pi_treat^(aa-1)*(1-pi_treat)^(2-aa))) + (1-pi_trial)/pi_trial*(muhat[,aa] - tau_trial[,aa])
    ifs[(n_trial+1):(n_trial+n_target),aa] <- tau_target[,aa]
  }
  # Survey adjustments
  if(survey){
    # eff_trial <- pi_treat>=epsilon & pi_treat<=1-epsilon
    # n_eff_trial <- sum(eff_trial)
    # ifs <- ifs[c(eff_trial, rep(TRUE, n_target)),]
    ifs <- ifs/(n_target)*(n_target+n_trial)
    mean_trial0 <- mean(ifs[1:n_trial,1])
    mean_trial1 <- mean(ifs[1:n_trial,2])
    mean_trial_diff <- mean(ifs[1:n_trial,2] - ifs[1:n_trial,1])
    sd_trial0 <- sd(ifs[1:n_trial,1])
    sd_trial1 <- sd(ifs[1:n_trial,2])
    sd_trial_diff <- sd(ifs[1:n_trial,2] - ifs[1:n_trial,1])
    nsfg_design <- 
      svydesign( 
        id = ~ nsfg_raw$SECU , 
        strata = ~ nsfg_raw$SEST , 
        data = nsfg , 
        weights = ~ nsfg_raw$WGT2015_2017 , 
        nest = TRUE 
      )
    nsfg_design <- 
      update( 
        nsfg_design , 
        psi0 = ifs[(n_trial+1):(n_trial+n_target),1],
        psi1 = ifs[(n_trial+1):(n_trial+n_target),2],
        effect = ifs[(n_trial+1):(n_trial+n_target),2] -  
          ifs[(n_trial+1):(n_trial+n_target),1])
    svy_res0 <- svymean( ~ psi0 , nsfg_design )
    mean_target0 <- coef(svy_res0)
    sd_target0 <- SE(svy_res0)
    svy_res1 <- svymean( ~ psi1 , nsfg_design )
    mean_target1 <- coef(svy_res1)
    sd_target1 <- SE(svy_res1)
    svy_res_diff <- svymean( ~ effect , nsfg_design )
    mean_target_diff <- coef(svy_res_diff)
    sd_target_diff <- SE(svy_res_diff)
    results <- matrix(0, nrow=2, ncol=3)
    results[1,1] <- n_trial*mean_trial0/(n_trial+n_target) + n_target*mean_target0/(n_trial+n_target)
    results[1,2] <- n_trial*mean_trial1/(n_trial+n_target) + n_target*mean_target1/(n_trial+n_target)
    results[1,3] <- n_trial*mean_trial_diff/(n_trial+n_target) + n_target*mean_target_diff/(n_trial+n_target)
    results[2,1] <- sqrt((n_trial*sd_trial0^2+n_target^2*sd_target0^2)/((n_trial+n_target)^2))
    results[2,2] <- sqrt((n_trial*sd_trial1^2+n_target^2*sd_target1^2)/((n_trial+n_target)^2))
    results[2,3] <- sqrt((n_trial*sd_trial_diff^2+n_target^2*sd_target_diff^2)/((n_trial+n_target)^2))
  }else{
    ifs <- ifs/(n_target)*(n_target+n_trial)
    results <- matrix(0, nrow=2, ncol=3)
    results[1,1:2] <- apply(ifs, 2, mean)
    results[1,3] <- results[1,2] - results[1,1]
    results[2,1:2] <- apply(ifs, 2, sd)/sqrt(n_trial + n_target)
    results[2,3] <- sd(ifs[,2]-ifs[,1])/sqrt(n_trial + n_target)
  }
  return(list(results=results, ifs=ifs))
}


# Plug-in
transport_plugin <- function(a, y, x_trial, x_target, sl.lib){
  set.seed(521)
  require(SuperLearner)
  n_trial <- nrow(x_trial)
  n_target <- nrow(x_target)
  cova_target <- colnames(x_target)
  muhat <- matrix(0, nrow=n_trial, ncol=2)
  tau_target <- matrix(0, nrow=n_target, ncol=2)
  for (aa in 1:2){
    mufit <- mcSuperLearner(y[a==aa-1],
                          x_trial[a==aa-1,],
                          newX=x_trial, SL.library=sl.lib)
    muhat[,aa] <- mufit$SL.predict
    taufit <- mcSuperLearner(muhat[,aa], x_trial[,cova_target], newX=x_target, SL.library=sl.lib)
    tau_target[,aa] <- taufit$SL.predict
  }
  potential_mean <- colMeans(tau_target)
  return(potential_mean[2]-potential_mean[1])
}