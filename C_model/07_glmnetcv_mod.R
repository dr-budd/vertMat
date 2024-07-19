## SET UP ----

# define models
models <- c("AllHomo",
            "FishDanio",
            "MammalsHomo",
            "ReptilesGallus")

# models <- c("ReptilesGallus")

## set working to script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

start_time <- Sys.time()
print(paste("Script start time: ", start_time))

## load libraries
library(caret)
library(glmnetUtils)
library(glmnet)
library(cocor)
library(plotmo)
library(tidyverse)

## FUNCTIONS ----

## scale evaluation (test or validate) matrix ----
scaleEvalMat <- function(eval_matrix, train_matrix, na.rm=TRUE) {
  ## check if both inputs are matrices
  if (!is.matrix(eval_matrix) || !is.matrix(train_matrix)) {
    stop("Both inputs must be matrices")
  }
  
  ## check if column names match
  if (!all(colnames(eval_matrix) %in% colnames(train_matrix))) {
    stop("Column names of test matrix not found in trainval matrix")
  }
  if (!all(colnames(train_matrix) %in% colnames(eval_matrix))) {
    stop("Column names of trainval matrix not found in test matrix")
  }
  
  ## diagnostic message
  message("Column names match successfully. Aligning columns for scaling.")
  
  ## reorder columns of train_matrix to match those of eval_matrix
  train_matrix <- train_matrix[, colnames(eval_matrix)]
  
  ## pre-compute mean and standard deviation values for each column in train_matrix
  means <- apply(train_matrix, 2, mean, na.rm = na.rm)
  sds <- apply(train_matrix, 2, sd, na.rm = na.rm)
  
  ## initialize a matrix to store the scaled data
  scaled_matrix <- matrix(nrow = nrow(eval_matrix), ncol = ncol(eval_matrix))
  
  ## apply column-wise mean-sd scaling using train_matrix mean and sd
  for (col_idx in 1:ncol(eval_matrix)) {
    mean_val <- means[col_idx]
    sd_val <- sds[col_idx]
    col <- eval_matrix[, col_idx]
    
    ## handle the case where sd is zero (to avoid division by zero)
    if (sd_val == 0) {
      scaled_matrix[, col_idx] <- ifelse(mean_val == 0, rep(0, length(col)), rep(1, length(col)))
    } else {
      scaled_matrix[, col_idx] <- (col - mean_val) / sd_val
    }
  }
  
  ## assign row and column names
  rownames(scaled_matrix) <- rownames(eval_matrix)
  colnames(scaled_matrix) <- colnames(eval_matrix)
  
  return(scaled_matrix)
}


## extract params ----

## specify function to extract the parameters 
## from the best model (by Daniel Freeman)
get_model_params <- function(alpha_crossval_train) {
  alpha <- alpha_crossval_train$alpha
  lambdaMin <- sapply(alpha_crossval_train$modlist, `[[`, "lambda.min")
  lambdaSE <- sapply(alpha_crossval_train$modlist, `[[`, "lambda.1se")
  error <- sapply(alpha_crossval_train$modlist, function(mod) {min(mod$cvm)})
  best <- which.min(error)
  data.frame(alpha=alpha[best], lambdaMin=lambdaMin[best],
             lambdaSE=lambdaSE[best], eror=error[best])
}

## unpaired t-test w var-test ---

## t-test and variance test function (transform skewed data)
tTestTrans <- function(df1, df2, values) {

  ## perform f-test to determine
  var <- var.test(log(df1 %>% pull({{values}})),
                  log(df2 %>% pull({{values}})),
                  alternative="two.sided")

  ## perform t-test based on variances
  if (var$p.value >= 0.05) {
    # print("variance equal")
    res <- t.test(log(df1 %>% pull({{values}})),
                  log(df2 %>% pull({{values}})),
                  var.equal=TRUE, ## equal (null)
                  paired=FALSE)
  } else {
    # print("variance unequal")
    res <- t.test(log(df1 %>% pull({{values}})),
                  log(df2 %>% pull({{values}})),
                  var.equal=FALSE, ## unequal (alt)
                  paired=FALSE)
  }

  res$p.value
}

## GLMNET ----

## loop through analysis
for (mod in models) {
  print("################################################################")
  print(paste0("########################### ", mod, " ###############################"))
  print("################################################################")

  ## print time
  Sys.time()

  ## lowercase
  mod_tl <- tolower(mod)
  
  ## IMPORT DATA ----
  
  ## define response
  response <- "mean_age_mat"
  
  ## print it
  paste0("metric type specified: ", response)
  
  ## import pheno data
  pheno_all <- read.csv(paste0("../B_exploreData/dataFiles/06.00_pheno_", mod_tl, ".csv"))
  
  ## subset for training and validation data
  pheno_trainval <- pheno_all %>% 
    filter(initial_split == "trainValidate")
  
  ## subset for test data
  pheno_test <- pheno_all %>% 
    filter(initial_split == "test")

  ## import CpG oe data
  cpg_oe_data_all <- read.csv(paste0("../B_exploreData/dataFiles/05.00_cpg_oe_data_",
                                       mod_tl, ".csv"))
  
  ## create feature matrix
  featmat_trainval <- cpg_oe_data_all %>%
    ## filter for trainval spp
    rename(organism_name=X) %>%
    filter(organism_name %in% pheno_trainval$organism_name) %>%
    ## add gc stuff
    left_join(pheno_trainval %>%
                mutate(prop_gc=gc_percent/100) %>%
                select(organism_name,
                       prop_gc)) %>%
    column_to_rownames(., var="organism_name") %>%
    ## only keep features with >= x % coverage (not including taxa)
    select(which(colMeans(!is.na(.)) >= 0.1)) %>%
    ## impute (replace NAs with) zeros
    replace(is.na(.), 0) %>%
    ## format as matrix
    as.matrix(.)

  ## convert categorical feature to binary ('one hot encoding')
  one_hot_clade_data <- model.matrix(~ order-1, data=pheno_trainval %>%
                                column_to_rownames(., var="organism_name"))

  ## combine the features
  featmat_trainval <- cbind(featmat_trainval, one_hot_clade_data)

  ## export
  write.csv(featmat_trainval, paste0("dataFiles/07.00_featmat_trainval_", mod_tl, ".csv"))

  ## export a list of all possible features
  all_features <- data.frame(features=colnames(featmat_trainval))
  write.csv(all_features,
            paste0("dataFiles/07.00_allfeatures_", mod_tl, ".csv"),
            row.names=FALSE)
  
  ## OUTER LOOP ----
  
  ## STEP 1: split the data into train and validate data sets ----
  
  ## create blanks to store data
  glmnetcv_parameters <- NULL
  glmnetcv_correlations <- NULL
  glmnetcv_phenoresults <- NULL ## trainval sep
  glmnet_phenoresults <- NULL ## final model predictions
  model_list <- list()
  
  ## create table containing all promoter ids to store weight information
  glmnet_weights <- data.frame(feature_id=c("(Intercept)", 
                                              colnames(featmat_trainval)))
  
  ## set seed for create data partition function
  set.seed(1)
  
  ## use the createDataPartition function to split the data (70/15/15)
  ## 10 times based on the percentiles of y
  sample_index_list <- createDataPartition(y=log(pheno_trainval[[response]]), 
                                           times=10, 
                                           p=14/17, ## 70 % of remaining 85 % (70/85*100)
                                           list=TRUE) 
  
  ## create counter
  part_num <- 0
  
  for (sample_index in sample_index_list) {
    
    ## add counter
    part_num <- part_num + 1
    print(paste0("partition number: ", part_num))
    
    ## use the sample_index to take the selected (train)
    featmat_train <- data.matrix(featmat_trainval[sample_index, ])
    
    ## use the inverse sample_index to select the remaining (validate)
    featmat_val <- data.matrix(featmat_trainval[-sample_index, ])
    
    ## do the same for the phenotype data 
    pheno_train <- pheno_trainval[sample_index, ]
    pheno_validate <- pheno_trainval[-sample_index, ]
    
    ## STEP 2: perform cross validation ----
    
    ## set seed for glmnet cross validation 
    set.seed(1)
    
    ## perform 10-fold cross validation for glmnet to determine optimal alpha and lambda values
    ## (using the train data)
    alpha_crossval_train <- cva.glmnet(x=featmat_train,
                                       y=log(pheno_train[[response]]),
                                       nfolds=10, 
                                       family="gaussian")
    
    ## get model parameters for the best alpha (the alpha that produces the minimum error)
    alpha_min <- get_model_params(alpha_crossval_train)$alpha
    lambda_min <- get_model_params(alpha_crossval_train)$lambdaMin
    lambda_1se <- get_model_params(alpha_crossval_train)$lambdaSE
    mse_min <- get_model_params(alpha_crossval_train)$eror
    
    ## STEP 3: make the model ----
    
    ## scale x
    featmat_train_scaled <- scale(featmat_train)
    featmat_train_scaled[is.nan(featmat_train_scaled)] <- 0
    
    ## fit a glm with elastic net regularisation to the train data
    glmnetcv_model <- glmnet(featmat_train_scaled,
                           scale(log(pheno_train[[response]])), 
                           family="gaussian",
                           standardize=FALSE, 
                           standardize.response=FALSE,
                           lambda=lambda_1se,
                           alpha=alpha_min)
    
    ## count the number of predictive features
    num_pred_proms <- nrow(data.frame(as.matrix(coef(glmnetcv_model, s=lambda_1se))) %>%
                             filter(s1 != 0))
    
    ## STEP 4: use the model to predict the data (to determine its performance) ----
    
    ## train ----
    
    ## use the model and the min lambda to predict log(response) for the scaled train data
    pheno_train$scaled_log_predicted_age_mat <- predict(glmnetcv_model, 
                                                   featmat_train_scaled, 
                                                   type="response", 
                                                   s=lambda_1se)[,'s1']
    
    ## transform the input value
    pheno_train$log_age_mat <- log(pheno_train[[response]])
    
    ## scale the input value
    pheno_train$scaled_log_age_mat <- scale(log(pheno_train[[response]]))
    
    ## unscale it
    pheno_train$log_predicted_age_mat <- pheno_train$scaled_log_predicted_age_mat * sd(pheno_train$log_age_mat) + mean(pheno_train$log_age_mat)
    
    ## detransform the predicted value
    pheno_train$detrans_predicted_age_mat <- exp(pheno_train$log_predicted_age_mat)
    
    ## calculate error
    pheno_train$detrans_absolute_error <- abs(pheno_train$detrans_predicted_age_mat - pheno_train$mean_age_mat)
    pheno_train$log_absolute_error <- abs(pheno_train$log_predicted_age_mat - pheno_train$log_age_mat)
    
    ## determine the pearson correlation in the train data set
    pearson_train <- cor.test(x=pheno_train$log_age_mat, 
                                 y=pheno_train$log_predicted_age_mat, 
                                 method="pearson")
    
    ## create vector for lower and upper conf. intervals
    conf_cols <- c("conf_0", "conf_1")
    
    ## create temporary table containing pearson cor info (train)
    ## NB: this will produce row.names related warnings...
    temp_train_table <- cbind.data.frame(pearson_train[], row.names=NULL) %>%
      ## add col for conf
      cbind(., conf_cols) %>%
      ## add mean squared error
      mutate(mse=mean((pheno_train$log_age_mat - pheno_train$log_predicted_age_mat)^2)) %>%
      ## add mean absolute error
      mutate(mae=mean(abs(pheno_train$log_age_mat - pheno_train$log_predicted_age_mat))) %>%
      ## add mean absolute error (detrans)
      mutate(medae=median(abs(pheno_train$log_age_mat - pheno_train$log_predicted_age_mat))) %>%
      ## add mean absolute error (detrans)
      mutate(medae_detrans=median(abs(pheno_train$mean_age_mat - pheno_train$detrans_predicted_age_mat))) %>%
      ## make 1 row
      pivot_wider(., names_from="conf_cols",
                  values_from="conf.int") %>%
      ## add suffixes
      rename_with( ~ paste0(.x, "_train"))
    
    ## validate ----
    
    featmat_val_scaled <- scaleEvalMat(featmat_val, featmat_train)
    featmat_val_scaled[is.nan(featmat_val_scaled)] <- 0
    
    ## use the model and the minimum lambda value to predict log(response) for the scaled validate data
    pheno_validate$scaled_log_predicted_age_mat <- predict(glmnetcv_model, 
                                                  featmat_val_scaled, 
                                                  type="response", 
                                                  s=lambda_1se)[,'s1']
    
    ## transform the input value
    pheno_validate$log_age_mat <- log(pheno_validate[[response]])
    
    ## scale the input value
    pheno_validate$scaled_log_age_mat <- scale(log(pheno_validate[[response]]))
    
    ## unscale it
    pheno_validate$log_predicted_age_mat <- pheno_validate$scaled_log_predicted_age_mat * sd(pheno_train$log_age_mat) + mean(pheno_train$log_age_mat)
    
    ## detransform the predicted value
    pheno_validate$detrans_predicted_age_mat <- exp(pheno_validate$log_predicted_age_mat)
    
    ## calculate error
    pheno_validate$detrans_absolute_error <- abs(pheno_validate$detrans_predicted_age_mat - pheno_validate$mean_age_mat)
    pheno_validate$log_absolute_error <- abs(pheno_validate$log_predicted_age_mat - pheno_validate$log_age_mat)
    
    ## determine the pearson correlation in the validate data set
    pearson_validate <- cor.test(x=pheno_validate$log_age_mat,
                                y=pheno_validate$log_predicted_age_mat,
                                method="pearson")
    
    ## create vector for lower and upper conf. intervals
    conf_cols <- c("conf_0", "conf_1")
    
    ## create temporary table containing pearson cor info (validate)
    temp_validate_table <- cbind.data.frame(pearson_validate[], row.names=NULL) %>%
      ## add col for conf
      cbind(., conf_cols) %>%
      ## add mean squared error
      mutate(mse=mean((pheno_validate$log_age_mat - pheno_validate$log_predicted_age_mat)^2)) %>%
      ## add mean absolute error
      mutate(mae=mean(abs(pheno_validate$log_age_mat - pheno_validate$log_predicted_age_mat))) %>%
      ## add median absolute error
      mutate(medae=median(abs(pheno_validate$log_age_mat - pheno_validate$log_predicted_age_mat))) %>%
      ## add median absolute error (detrans)
      mutate(medae_detrans=median(abs(pheno_validate$mean_age_mat - pheno_validate$detrans_predicted_age_mat))) %>%
      ## make 1 row
      pivot_wider(., names_from="conf_cols",
                  values_from="conf.int") %>%
      ## add suffixes
      rename_with( ~ paste0(.x, "_validate"))
    
    ## t-test ae (for signs of overfitting) and pull significance code 
    sig_diff_error_detrans <- tTestTrans(df1=pheno_train, 
                            df2=pheno_validate, 
                            values=detrans_absolute_error)

    ## create a df with model parameters
    temp_parameters <- cbind(response,
                             model_number=part_num,
                             alpha_min,
                             lambda_min,
                             lambda_1se,
                             mse_min,
                             num_pred_proms,
                             sig_diff_error_detrans)
    
    ## combine your tables and the response and seed info
    temp_table <- cbind(response,
                        model_number=part_num,
                        temp_train_table,
                        temp_validate_table)

    ## append model parameters to master
    glmnetcv_parameters <- rbind(temp_parameters, glmnetcv_parameters)
    
    ## append correlations etc to master
    glmnetcv_correlations <- rbind(temp_table, glmnetcv_correlations)
    
    ## and for the pheno/results data
    pheno_train$trainval_split <- "train"
    pheno_validate$trainval_split <- "validate"
    pheno_both <- rbind(pheno_train, pheno_validate) %>%
      mutate(model_number=part_num)

    glmnetcv_phenoresults <- rbind(glmnetcv_phenoresults, pheno_both)
    
    ## FINAL MODEL (if statement) ----
    
    # only if p > 0.05!!!!!
    if (sig_diff_error_detrans > 0.05) {
      
      ## refit model to all data ----
      
      featmat_trainval_scaled <- scale(featmat_trainval)
      featmat_trainval_scaled[is.nan(featmat_trainval_scaled)] <- 0
      
      ## fit a glm with elastic net regularisation to ALL training data (training+validation)
      glmnet_model <- glmnet(featmat_trainval_scaled,
                             scale(log(pheno_trainval[[response]])), 
                             family="gaussian",
                             standardize=FALSE, 
                             standardize.response=FALSE,
                             lambda=lambda_1se,
                             alpha=alpha_min)
      
      ## predict for training and validation set----
      
      ## use the new model and the min lambda to predict log(response) for ALL the data
      pheno_trainval$scaled_log_predicted_age_mat <- predict(glmnet_model, 
                                                      featmat_trainval_scaled,
                                                      type="response", 
                                                      s=lambda_1se)[,'s1']
      
      ## transform the input value
      pheno_trainval$log_age_mat <- log(pheno_trainval[[response]])
      
      ## scale the input value
      pheno_trainval$scaled_log_age_mat <- scale(log(pheno_trainval[[response]]))
      
      ## unscale it
      pheno_trainval$log_predicted_age_mat <- pheno_trainval$scaled_log_predicted_age_mat * sd(pheno_trainval$log_age_mat) + mean(pheno_trainval$log_age_mat)
      
      ## detransform the predicted value
      pheno_trainval$detrans_predicted_age_mat <- exp(pheno_trainval$log_predicted_age_mat)
      
      ## calculate error
      pheno_trainval$detrans_absolute_error <- abs(pheno_trainval$mean_age_mat - pheno_trainval$detrans_predicted_age_mat)
      pheno_trainval$log_absolute_error <- abs(pheno_trainval$log_age_mat - pheno_trainval$log_predicted_age_mat)
      
      ## save final results ----
      
      ## add partition number
      pheno_trainval$model_number <- part_num
      
      ## save to master
      glmnet_phenoresults <- rbind(glmnet_phenoresults, pheno_trainval)
      
      ## save final models ----
      
      ## add lambda 1se to model object
      attr(glmnet_model, "lambda.1se") <- lambda_1se
      
      ## add alpha min to model object
      attr(glmnet_model, "alpha_min") <- alpha_min
      
      ## add scaling info to model object
      attr(glmnet_model, "mean.logy.trainval") <- mean(pheno_trainval$log_age_mat)
      attr(glmnet_model, "sd.logy.trainval") <- sd(pheno_trainval$log_age_mat)
      
      ## add model to list
      model_list[[part_num]] <- glmnet_model
      
      ## extract the weights for the specified value of lambda
      features <- as.matrix(coef(glmnet_model, s=lambda_1se))
      
      ## select the non-zero weights and label column with par num
      feature_weights <- data.frame(features) %>%
        ## for long
        filter(s1 != 0) %>%
        # rename(weight=s1) %>%
        # mutate(model_number=part_num) %>%
        ## OR for wide
        rename(!!paste0("model_", part_num) := s1) %>%
        rownames_to_column(., var="feature_id")
      
      ## append model weights to master
      glmnet_weights <- left_join(glmnet_weights, feature_weights, 
                                  by=join_by(feature_id))
      
    } else {
      message <- paste0("p < 0.05, model: ", part_num, " discarded\n")
      cat(message) # Print the message
    }
    
  }
  
  ## format results ----
  
  ## format glmnetcv_phenoresults (train and validate separately)
  glmnetcv_phenoresults <- glmnetcv_phenoresults %>%
    mutate(log_relative_error=(log_absolute_error/log_age_mat)*100) %>%
    mutate(detrans_relative_error=(detrans_absolute_error/mean_age_mat)*100) %>%
    mutate(model_number_labs=paste0("Model 0", model_number)) %>%
    mutate(model_number_labs=gsub("010", "10", model_number_labs)) %>%
    mutate(trainval_split_labs=ifelse(trainval_split == "train",
                                        "Training",
                                        "Validation")) %>%
    mutate(model=mod_tl) %>%
    mutate(model_labs=ifelse(model=="allhomo", "All-vertebrate", 
                               ifelse(model=="mammalshomo", "Mammal-specific", 
                                      ifelse(model=="reptilesgallus", "Reptile-specific", 
                                             ifelse(model=="fishdanio", "Fish-specific", 
                                                    "error")))))
  
  ## format glmnet_phenoresults (tranval together)
  glmnet_phenoresults <- glmnet_phenoresults %>%
    mutate(log_relative_error=(log_absolute_error/log_age_mat)*100) %>%
    mutate(detrans_relative_error=(detrans_absolute_error/mean_age_mat)*100) %>%
    mutate(model_number_labs=paste0("Model 0", model_number)) %>%
    mutate(model_number_labs=gsub("010", "10", model_number_labs)) %>%
    mutate(model=mod_tl) %>%
    mutate(model_labs=ifelse(model=="allhomo", "All-vertebrate", 
                               ifelse(model=="mammalshomo", "Mammal-specific", 
                                      ifelse(model=="reptilesgallus", "Reptile-specific", 
                                             ifelse(model=="fishdanio", "Fish-specific", 
                                                    "error")))))
  
  ## format weights
  glmnet_weights <- glmnet_weights %>%
    column_to_rownames(., var="feature_id") %>%
    filter_all(any_vars(!is.na(.))) %>%
    rownames_to_column(., var="feature_id")
  
  ## bag predictions (training/validation)----
  
  ## bag the model estimates
  bagged_phenoresults <- glmnet_phenoresults %>%
    group_by(organism_name) %>%
    mutate(bagged_detrans_predicted_age_mat=median(detrans_predicted_age_mat), 
           bagged_log_predicted_age_mat=median(log_predicted_age_mat),
           # bagged_detrans_predicted_age_mat=mean(detrans_predicted_age_mat), 
           # bagged_log_predicted_age_mat=mean(log_predicted_age_mat),
           bagged_detrans_absolute_error=abs(mean_age_mat - bagged_detrans_predicted_age_mat), 
           bagged_detrans_relative_error=(bagged_detrans_absolute_error/mean_age_mat)*100) %>%
    ## remove anything model-specific
    select(-log_predicted_age_mat, -scaled_log_predicted_age_mat, -detrans_predicted_age_mat,
           -model_number, -log_absolute_error, -log_relative_error, 
           -detrans_absolute_error, -detrans_relative_error, -model_number_labs) %>%
    unique(.) %>%
    as.data.frame(.)
  
  ## export data ----
  
  ## bagging
  write.csv(bagged_phenoresults, 
            paste0("dataFiles/07.00_bagged_phenoresults_", mod_tl, ".csv"), 
            row.names=FALSE)
  
  ## export model list for predictions
  saveRDS(model_list, paste0("dataFiles/07.00_modellist_", mod_tl, ".rds"))

  ## export parameters, correlations, weights and predictions
  write.csv(glmnetcv_parameters,
            paste0("dataFiles/07.00_glmnetcv_parameters_", mod_tl, ".csv"),
            row.names=FALSE)
  write.csv(glmnetcv_correlations,
            paste0("dataFiles/07.00_glmnetcv_correlations_", mod_tl, ".csv"),
            row.names=FALSE)
  write.csv(glmnetcv_phenoresults, 
            paste0("dataFiles/07.00_glmnetcv_phenoresults_", mod_tl, ".csv"),
            row.names=FALSE)
  
  ## and for final model
  write.csv(glmnet_weights, 
            paste0("dataFiles/07.00_glmnet_weights_", mod_tl, ".csv"),
            row.names=FALSE)
  write.csv(glmnet_phenoresults, 
            paste0("dataFiles/07.00_glmnet_phenoresults_", mod_tl, ".csv"),
            row.names=FALSE)
  
  ## TEST DATA ----
  
  ## format test data ----
  
  ## convert categorical feature to binary ('one hot encoding') for testing data (separately)
  one_hot_clade_data_test <- model.matrix(~order-1, data=pheno_test %>%
                                       column_to_rownames(., var="organism_name"))
  
  ## create feature matrix
  featmat_test <- cpg_oe_data_all %>%
    ## filter for test spp
    rename(organism_name=X) %>%
    filter(organism_name %in% pheno_test$organism_name) %>%
    ## add gc stuff
    left_join(pheno_test %>%
                mutate(prop_gc=gc_percent/100) %>%
                select(organism_name,
                       prop_gc)) %>%
    column_to_rownames(., var="organism_name") %>%
    ## impute zeros
    replace(is.na(.), 0) 
  
  ## get features from model
  unique_features <- c()
  for (i in 1:length(model_list)) {
    model <- model_list[[i]]
    # Extract coefficients (excluding the intercept) for non-zero features
    features <- rownames(coef(model))[-1]
    # Add unique features to the vector
    unique_features <- union(unique_features, features)
  }
  
  ## add clade data and subset for model feats
  featmat_test <- cbind(featmat_test, one_hot_clade_data_test) %>%
    ## subset for features with coef > 0
    select(., any_of(unique_features)) %>%
    ## format as matrix
    as.matrix(.)
  
  # add missing one-hot endocded data so you can model it
  missing_features <- setdiff(unique_features, colnames(featmat_test))
  zero_data <- matrix(0, nrow=nrow(featmat_test), ncol=length(missing_features))
  colnames(zero_data) <- missing_features
  featmat_test <- cbind(featmat_test, zero_data)
  
  ## check they match 
  if (ncol(featmat_test) != length(unique_features)) {
    stop("Error: The number of columns in featmat_test is not equal to the length of features.")
  }
  
  ## scale it
  featmat_test_scaled <- scaleEvalMat(featmat_test, featmat_trainval)
  featmat_test_scaled[is.nan(featmat_test_scaled)] <- 0
  
  ## predict for test set ----
  
  ## make predictions using the ensemble
  
  glmnet_phenoresults_testpreds <- NULL
  
  for (i in 1:length(model_list)) {
    
    glmnet_model <- model_list[[i]]
    
    if (is.null(glmnet_model)) {
      print(paste0("No model no. ", i, " for ", mod_tl))
    } else {
      print(paste0("Predicting test set for ", mod_tl, " model no. ", i))
      ## retrieve lambda
      lambda_1se <- attr(glmnet_model, "lambda.1se")
      
      ## use the model and the minimum lambda value to predict log(pheno) for the training
      pheno_test$scaled_log_predicted_age_mat <- predict(glmnet_model, 
                                                  featmat_test_scaled, 
                                                  type="response", 
                                                  s=lambda_1se)[,'s1']
      
      ## transform the input value
      pheno_test$log_age_mat <- log(pheno_test[[response]])
      
      ## scale the input value
      pheno_test$scaled_log_age_mat <- (pheno_test$log_age_mat - mean(pheno_trainval$log_age_mat))/
        sd(pheno_trainval$log_age_mat)
      
      ## unscale it
      pheno_test$log_predicted_age_mat <- pheno_test$scaled_log_predicted_age_mat * sd(pheno_trainval$log_age_mat) + mean(pheno_trainval$log_age_mat)
      
      ## detransform the predicted value
      pheno_test$detrans_predicted_age_mat <- exp(pheno_test$log_predicted_age_mat)
      
      ## calculate error
      pheno_test$log_absolute_error <- pheno_test$log_age_mat - pheno_test$log_predicted_age_mat
      pheno_test$detrans_absolute_error <- pheno_test$mean_age_mat - pheno_test$detrans_predicted_age_mat
      
      ## add model number
      pheno_test$model_number <- i
      
      ## create pheno results
      glmnet_phenoresults_testpreds <- rbind(glmnet_phenoresults_testpreds, pheno_test)
    }
  }
  
  ## format results ----
  
  glmnet_phenoresults_testpreds <- glmnet_phenoresults_testpreds %>%
    mutate(log_relative_error=(log_absolute_error/log_age_mat)*100) %>%
    mutate(detrans_relative_error=(detrans_absolute_error/mean_age_mat)*100) %>%
    mutate(model_number_labs=paste0("Model 0", model_number)) %>%
    mutate(model_number_labs=gsub("010", "10", model_number_labs)) %>%
    mutate(model=mod_tl) %>%
    mutate(model_labs=ifelse(model=="allhomo", "All-vertebrate", 
                               ifelse(model=="mammalshomo", "Mammal-specific", 
                                      ifelse(model=="reptilesgallus", "Reptile-specific", 
                                             ifelse(model=="fishdanio", "Fish-specific", 
                                                    "error")))))
  
  ## export test data
  write.csv(glmnet_phenoresults_testpreds,
            paste0("dataFiles/07.00_glmnet_phenoresults_testpreds_",
                   mod_tl, ".csv"), row.names=FALSE)
  
  ## bag predictions (test data) ----
  
  ## bag the model estimates
  bagged_phenoresults_test <- glmnet_phenoresults_testpreds %>%
    group_by(organism_name) %>%
    # mutate(bagged_detrans_predicted_age_mat=mean(detrans_predicted_age_mat), 
    #        bagged_log_predicted_age_mat=mean(log_predicted_age_mat)) %>%
    mutate(bagged_detrans_predicted_age_mat=median(detrans_predicted_age_mat), 
           bagged_log_predicted_age_mat=median(log_predicted_age_mat),
           bagged_detrans_absolute_error=abs(mean_age_mat - bagged_detrans_predicted_age_mat), 
           bagged_detrans_relative_error=(bagged_detrans_absolute_error/mean_age_mat)*100) %>%
    ## remove anything model-specific
    select(-log_predicted_age_mat, -scaled_log_predicted_age_mat, -detrans_predicted_age_mat,
           -model_number, -log_absolute_error, -log_relative_error, 
           -detrans_absolute_error, -detrans_relative_error, -model_number_labs) %>%
    unique(.) %>%
    as.data.frame(.)
  
  write.csv(bagged_phenoresults_test, 
            paste0("dataFiles/07.00_bagged_phenoresults_test_", mod_tl, ".csv"), 
            row.names=FALSE)
  
  ## format data for prediction intervals ----
  
  ## combine test and train (for ensemble data)
  glmnet_phenoresults_all <- read.csv(paste0("dataFiles/07.00_glmnet_phenoresults_", mod_tl, ".csv")) %>%
    bind_rows(glmnet_phenoresults_testpreds)
  
  ## import feature weights
  # feature_weights <- read.csv(paste0("dataFiles/07.00_glmnet_weights_", mod_tl, ".csv"))
  
  ## import feature matrix (trainval only)
  featmat_trainval <- read.csv(paste0("dataFiles/07.00_featmat_trainval_", mod_tl, ".csv"), 
                               row.names=1)
  
  ## initiate df
  ensemble_data <- NULL
  
  for (mod_num in levels(factor(glmnet_phenoresults_all$model_number))) {
    
    ## select only relevent columns from pheno data for each partition
    glmnet_phenoresults_all_less <- glmnet_phenoresults_all %>%
      ## add bagged prediction
      group_by(organism_name) %>%
      mutate(bagged_detrans_predicted_age_mat=median(detrans_predicted_age_mat)) %>%
      # mutate(bagged_scaled_log_predicted_age_mat=median(scaled_log_predicted_age_mat)) %>%
      ungroup (.) %>%
      ## filter by partition
      filter(model_number == mod_num) %>%
      ## select relevant columns
      select(organism_name, mean_age_mat, log_age_mat, model_number,
             initial_split, log_predicted_age_mat, detrans_predicted_age_mat,
             bagged_detrans_predicted_age_mat) %>%
      ## arrange alphabetically
      arrange(organism_name)
    
    ## define column
    col <- paste0("model_", mod_num)
    
    ## create temp name and print it
    name_temp <- paste0(col, "_data")
    print(name_temp)
    
    ## convert col into a symbol to be evaluated within the pipe
    col <- rlang::sym(col)
    
    # ## select only relevant column from feature weight df
    # predictive_features_temp <- feature_weights %>%
    #   select(feature_id, !!col) %>%
    #   filter(feature_id != "(Intercept)") %>%
    #   ## retain rows of col that are not 0 (or NA)
    #   filter(!!col != 0) %>%
    #   pull(feature_id)
    
    ## (don't) select those specific features from model matrices
    data_test_temp <- as.data.frame(featmat_test) %>% ## raw
    # data_test_temp <- as.data.frame(featmat_test_scaled) %>% ## scaled
      # select(all_of(predictive_features_temp)) %>%
      rownames_to_column(., var="organism_name") 
    
    data_train_temp <- as.data.frame(featmat_trainval) %>% ## raw
    # data_train_temp <- as.data.frame(featmat_trainval_scaled) %>% ## scaled
      # select(all_of(predictive_features_temp)) %>%
      rownames_to_column(., var="organism_name")
    
    ## combine and add pheno data
    data_temp <- bind_rows(data_test_temp, data_train_temp) %>%
      ## add to pheno data
      full_join(glmnet_phenoresults_all_less, .)
    
    ## append to data
    ensemble_data <- bind_rows(ensemble_data, data_temp)
    
  }
  
  ensemble_data <- ensemble_data %>%
    relocate(sort(names(.))) %>%
    relocate(all_of(c(starts_with("order"), starts_with("FP")))) %>%
    relocate(organism_name, mean_age_mat,
             log_age_mat, model_number,
             initial_split, log_predicted_age_mat,
             detrans_predicted_age_mat, bagged_detrans_predicted_age_mat)
  
  ## export data
  write.csv(ensemble_data,
            paste0("dataFiles/07.00_ensemble_data_", mod_tl, ".csv"),
            row.names=FALSE)
  
}

## print session info
sessionInfo()

## print end time
end_time <- Sys.time()
print(paste("Script end time: ", end_time))
print(paste("Elapsed time: ",  end_time - start_time))

## END SCRIPT
##

