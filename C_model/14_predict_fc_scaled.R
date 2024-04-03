## SET UP ----

## define models
models<-c("AllHomoFc")

## for testing
# mod<-c("AllHomoFc")

## set working to script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

start_time<-Sys.time()
print(paste("Script start time: ", start_time))

## load libraries
library(caret)
library(glmnetUtils)
library(glmnet)
library(cocor)
library(plotmo)
library(tidyverse)

## FUNCTIONS ----

scaleEvalMat<-function(eval_matrix, train_matrix, na.rm=TRUE) {
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
  train_matrix<-train_matrix[, colnames(eval_matrix)]
  
  ## pre-compute mean and standard deviation values for each column in train_matrix
  means<-apply(train_matrix, 2, mean, na.rm=na.rm)
  sds<-apply(train_matrix, 2, sd, na.rm=na.rm)
  
  ## initialize a matrix to store the scaled data
  scaled_matrix<-matrix(nrow=nrow(eval_matrix), ncol=ncol(eval_matrix))
  
  ## apply column-wise mean-sd scaling using train_matrix mean and sd
  for (col_idx in 1:ncol(eval_matrix)) {
    mean_val<-means[col_idx]
    sd_val<-sds[col_idx]
    col<-eval_matrix[, col_idx]
    
    ## handle the case where sd is zero (to avoid division by zero)
    if (sd_val == 0) {
      scaled_matrix[, col_idx]<-ifelse(mean_val == 0, rep(0, length(col)), rep(1, length(col)))
    } else {
      scaled_matrix[, col_idx]<-(col - mean_val) / sd_val
    }
  }
  
  ## assign row and column names
  rownames(scaled_matrix)<-rownames(eval_matrix)
  colnames(scaled_matrix)<-colnames(eval_matrix)
  
  return(scaled_matrix)
}

## FOR ALL 4 MODELS ... ----

## loop through analysis
for (mod in models) {
  print("################################################################")
  print(paste0("########################### ", mod, " ###############################"))
  print("################################################################")

  ## print time
  Sys.time()

  ## lowercase
  mod_tl<-tolower(mod)
  
  ## define response
  response<-"mean_age_mat"
  
  ## print it
  paste0("metric type specified: ", response)
  
  ## import model list
  bagged_mod_list<-readRDS(paste0("dataFiles/07.00_modellist_", 
                               gsub("fc", "", mod_tl), 
                                    ".rds"))
  
  ## import pheno data (no age mat)
  pheno_genomic<-read.csv(paste0("../A_getGenomicData/dataFilesFc/07.00_pheno_genomic_", 
                                    mod_tl, 
                                    ".csv"))
  
  pheno_fc<-read.csv(paste0("../B_exploreData/dataFiles/00.03_pheno_fc.csv")) %>%
    left_join(pheno_genomic, .)
  
  ## import CpG oe data
  cpg_oe_data_fc<-read.csv(paste0("../A_getGenomicData/dataFilesFc/07.00_cpg_oe_data_full_",
                                       mod_tl, ".csv"))
  
  ## PREDICT FC ----
  
  ## format fc data ----
  
  ## convert categorical feature to binary ('one hot encoding') for fc data
  one_hot_clade_data_fc<-model.matrix(~order-1, data=pheno_fc %>%
                                       column_to_rownames(., var="assembly_id"))
  
  ## create feature matrix
  featmat_fc<-cpg_oe_data_fc %>%
    ## filter for fc spp
    rename(assembly_id=X) %>%
    filter(assembly_id %in% pheno_fc$assembly_id) %>%
    ## add gc stuff
    left_join(pheno_fc %>%
                mutate(prop_gc=gc_percent/100) %>%
                select(assembly_id,
                       prop_gc)) %>%
    column_to_rownames(., var="assembly_id") %>%
    ## impute zeros
    replace(is.na(.), 0) 
  
  ## get all features from all models
  unique_features<-c()
  for (i in 1:length(bagged_mod_list)) {
    model<-bagged_mod_list[[i]]
    ## extract coefficients (excluding the intercept) for non-zero features
    features<-rownames(coef(model))[-1]
    ## add unique features to the vector
    unique_features<-union(unique_features, features)
  }
  
  ## add clade data and subset for model feats
  featmat_fc<-cbind(featmat_fc, one_hot_clade_data_fc) %>%
    ## subset for features with coef > 0
    select(., any_of(unique_features)) %>%
    ## format as matrix
    as.matrix(.)
  
  # add missing one-hot encoded data so you can model it
  missing_features<-setdiff(unique_features, colnames(featmat_fc))
  zero_data<-matrix(0, nrow=nrow(featmat_fc), ncol=length(missing_features))
  colnames(zero_data)<-missing_features
  featmat_fc<-cbind(featmat_fc, zero_data)
  
  ## check they match 
  if (ncol(featmat_fc) != length(unique_features)) {
    stop("Error: The number of columns in featmat_fc is not equal to the length of features.")
  }
  
  ## scale data ----
  
  ## import trainval matrix
  featmat_trainval<-read.csv(paste0("dataFiles/07.00_featmat_trainval_", 
                                    gsub("fc", "", mod_tl), 
                                    ".csv")) %>%
    column_to_rownames(var="X") %>% as.matrix(.)
  
  ## scale features
  featmat_fc_scaled<-scaleEvalMat(featmat_fc, featmat_trainval)
  
  ## predict for fcs ----
  
  ## make predictions using the ensemble
  
  glmnet_phenoresults_fcpreds<-NULL
  
  for (i in 1:length(bagged_mod_list)) {
    
    glmnet_model<-bagged_mod_list[[i]]
    
    if (is.null(glmnet_model)) {
      print(paste0("No model no. ", i, " for ", mod_tl))
    } else {
      print(paste0("Predicting FC set for ", mod_tl, " model no. ", i))
      
      ## retrieve lambda
      lambda_1se<-attr(glmnet_model, "lambda.1se")
      
      ## use the model and the minimum lambda value to predict log(pheno) for the training data
      pheno_fc$scaled_log_predicted_age_mat<-predict(glmnet_model, 
                                                  featmat_fc_scaled, 
                                                  type="response", 
                                                  s=lambda_1se)[,'s1']
      
      ## unscale the value
      pheno_fc$log_predicted_age_mat<-(pheno_fc$scaled_log_predicted_age_mat*attr(glmnet_model,"sd.logy.trainval"))+attr(glmnet_model,"mean.logy.trainval")
      
      ## detransform the predicted value
      pheno_fc$detrans_predicted_age_mat<-exp(pheno_fc$log_predicted_age_mat)
      
      ## add model number
      pheno_fc$model_number<-i
      
      ## create pheno results
      glmnet_phenoresults_fcpreds<-rbind(glmnet_phenoresults_fcpreds, pheno_fc)
    }
  }
  
  ## format results ----
  
  glmnet_phenoresults_fcpreds<-glmnet_phenoresults_fcpreds %>%
    mutate(model_number_labs=paste0("Model 0", model_number)) %>%
    mutate(model_number_labs=gsub("010", "10", model_number_labs)) %>%
    mutate(model=mod_tl) %>%
    mutate(model_labs=ifelse(model=="allhomofc", "All-vertebrate", 
                               ifelse(model=="mammalshomofc", "Mammal-specific", 
                                      ifelse(model=="reptilesgallusfc", "Reptile-specific", 
                                             ifelse(model=="fishdaniofc", "Fish-specific", 
                                                    "error")))))
  
  ## export fc data
  write.csv(glmnet_phenoresults_fcpreds,
            paste0("dataFiles/14.00_glmnet_phenoresults_fcpreds_",
                   mod_tl, ".csv"), row.names=FALSE)
  
  ## bag predictions (fc data) ----
  
  ## bag the model estimates
  bagged_phenoresults_fc<-glmnet_phenoresults_fcpreds %>%
    group_by(assembly_id) %>%
    mutate(bagged_detrans_predicted_age_mat=median(detrans_predicted_age_mat), 
           bagged_log_predicted_age_mat=median(log_predicted_age_mat)) %>%
    ## remove anything model-specific
    select(-log_predicted_age_mat, -scaled_log_predicted_age_mat,
           -detrans_predicted_age_mat,
           -model_number, -model_number_labs) %>%
    unique(.) %>%
    as.data.frame(.)
  
  write.csv(bagged_phenoresults_fc,
            paste0("dataFiles/14.00_bagged_phenoresults_", mod_tl, ".csv"),
            row.names=FALSE)
  
}

## print session info
sessionInfo()

## print end time
end_time<-Sys.time()
print(paste("Script end time: ", end_time))
print(paste("Elapsed time: ",  end_time - start_time))

## END SCRIPT
##

