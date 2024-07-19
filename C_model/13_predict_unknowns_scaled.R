## SET UP ----

## define models
models<-c("AllHomoUnknowns",
          "FishDanioUnknowns",
          "MammalsHomoUnknowns",
          "ReptilesGallusUnknowns")

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
                                  gsub("unknowns", "", mod_tl), 
                                  ".rds"))
  
  ## import pheno data (no age mat)
  pheno_genomic<-read.csv(paste0("../A_getGenomicData/dataFilesUnknowns/04.00_pheno_genomic_", 
                                 mod_tl, 
                                 ".csv"))
  
  pheno_unknowns<-read.csv(paste0("../B_exploreData/dataFiles/00.03_pheno_unknowns.csv")) %>%
    left_join(pheno_genomic, .)
  
  ## import CpG oe data
  cpg_oe_data_unknowns<-read.csv(paste0("../A_getGenomicData/dataFilesUnknowns/04.00_cpg_oe_data_full_",
                                        mod_tl, ".csv"))
  
  ## PREDICT UNKNOWNS ----
  
  ## format unknowns data ----
  
  ## convert categorical feature to binary ('one hot encoding') for unkowns data
  one_hot_clade_data_unknowns<-model.matrix(~order-1, data=pheno_unknowns %>%
                                              column_to_rownames(., var="organism_name"))
  
  ## create feature matrix
  featmat_unknowns<-cpg_oe_data_unknowns %>%
    ## filter for unknown spp
    rename(assembly_id=X) %>%
    filter(assembly_id %in% pheno_unknowns$assembly_id) %>%
    ## add gc stuff
    left_join(pheno_unknowns %>%
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
  featmat_unknowns<-cbind(featmat_unknowns, one_hot_clade_data_unknowns) %>%
    ## subset for features with coef > 0
    select(., any_of(unique_features)) %>%
    ## format as matrix
    as.matrix(.)
  
  ## add missing one-hot encoded data so you can model it
  missing_features<-setdiff(unique_features, colnames(featmat_unknowns))
  zero_data<-matrix(0, nrow=nrow(featmat_unknowns), ncol=length(missing_features))
  colnames(zero_data)<-missing_features
  featmat_unknowns<-cbind(featmat_unknowns, zero_data)
  
  ## check they match 
  if (ncol(featmat_unknowns) != length(unique_features)) {
    stop("Error: The number of columns in featmat_unknowns is not equal to the length of features.")
  }
  
  ## scale data ----
  
  ## import trainval matrix
  featmat_trainval<-read.csv(paste0("dataFiles/07.00_featmat_trainval_", 
                                    gsub("unknowns", "", mod_tl), 
                                    ".csv")) %>%
    column_to_rownames(var="X") %>% as.matrix(.)
  
  ## scale features
  featmat_unknowns_scaled<-scaleEvalMat(featmat_unknowns, featmat_trainval)
  
  ## predict for unknowns ----
  
  ## make predictions using the ensemble
  
  glmnet_phenoresults_unknownspreds<-NULL
  
  for (i in 1:length(bagged_mod_list)) {
    
    glmnet_model<-bagged_mod_list[[i]]
    
    if (is.null(glmnet_model)) {
      print(paste0("No model no. ", i, " for ", mod_tl))
    } else {
      print(paste0("Predicting unknown set for ", mod_tl, " model no. ", i))
      
      ## retrieve lambda
      lambda_1se<-attr(glmnet_model, "lambda.1se")
      
      ## use the model and the minimum lambda value to predict log(pheno) for the training data
      pheno_unknowns$scaled_log_predicted_age_mat<-predict(glmnet_model, 
                                                           featmat_unknowns_scaled, 
                                                           type="response", 
                                                           s=lambda_1se)[,'s1']
      
      ## unscale the value
      pheno_unknowns$log_predicted_age_mat<-(pheno_unknowns$scaled_log_predicted_age_mat*attr(glmnet_model,"sd.logy.trainval"))+attr(glmnet_model,"mean.logy.trainval")
      
      ## detransform the predicted value
      pheno_unknowns$detrans_predicted_age_mat<-exp(pheno_unknowns$log_predicted_age_mat)
      
      ## add model number
      pheno_unknowns$model_number<-i
      
      ## create pheno results
      glmnet_phenoresults_unknownspreds<-rbind(glmnet_phenoresults_unknownspreds, pheno_unknowns)
    }
  }
  
  ## format results ----
  
  glmnet_phenoresults_unknownspreds<-glmnet_phenoresults_unknownspreds %>%
    mutate(model_number_labs=paste0("Model 0", model_number)) %>%
    mutate(model_number_labs=gsub("010", "10", model_number_labs)) %>%
    mutate(model=mod_tl) %>%
    mutate(model_labs=ifelse(model=="allhomounknowns", "All-vertebrate", 
                             ifelse(model=="mammalshomounknowns", "Mammal-specific", 
                                    ifelse(model=="reptilesgallusunknowns", "Reptile-specific", 
                                           ifelse(model=="fishdaniounknowns", "Fish-specific", 
                                                  "error")))))
  
  ## export unknowns data
  write.csv(glmnet_phenoresults_unknownspreds,
            paste0("dataFiles/13.00_glmnet_phenoresults_unknownspreds_",
                   mod_tl, ".csv"), row.names=FALSE)
  
  ## bag predictions (unknowns data) ----
  
  ## bag the model estimates
  bagged_phenoresults_unknowns<-glmnet_phenoresults_unknownspreds %>%
    group_by(organism_name) %>%
    mutate(bagged_detrans_predicted_age_mat=median(detrans_predicted_age_mat), 
           bagged_log_predicted_age_mat=median(log_predicted_age_mat)) %>%
    ## remove anything model-specific
    select(-log_predicted_age_mat, -scaled_log_predicted_age_mat,
           -detrans_predicted_age_mat,
           -model_number, -model_number_labs) %>%
    unique(.) %>%
    as.data.frame(.)
  
  ## format data for prediction intervals (PYTHON) ----
  
  ## import feature weights (model coefficients)
  feature_weights<-read.csv(paste0("dataFiles/07.00_glmnet_weights_", 
                                   gsub("unknowns", "", mod_tl), 
                                   ".csv"))
  
  
  ## initiate df
  ensemble_data<-NULL
  
  for (mod_num in levels(factor(glmnet_phenoresults_unknownspreds$model_number))) {
    
    ## select only relevent columns from pheno data for each partition
    glmnet_phenoresults_all_less<-glmnet_phenoresults_unknownspreds %>%
      ## add bagged prediction
      group_by(organism_name) %>%
      mutate(bagged_detrans_predicted_age_mat=median(detrans_predicted_age_mat)) %>%
      ungroup (.) %>%
      ## filter by partition
      filter(model_number == mod_num) %>%
      ## select relevant columns
      select(organism_name, model_number, assembly_id,
             log_predicted_age_mat, detrans_predicted_age_mat,
             bagged_detrans_predicted_age_mat) %>%
      ## arrange alphabetically
      arrange(organism_name)
    
    ## define column
    col<-paste0("model_", mod_num)
    
    ## create temp name and print it
    name_temp<-paste0(col, "_data")
    print(name_temp)
    
    ## convert col into a symbol to be evaluated within the pipe
    col<-rlang::sym(col)
    
    # ## select only relevant column from feature weight df
    # predictive_features_temp<-feature_weights %>%
    #   select(feature_id, !!col) %>%
    #   filter(feature_id != "(Intercept)") %>%
    #   ## retain rows of col that are not 0 (or NA)
    #   filter(!!col != 0) %>%
    #   pull(feature_id)
    
    ## select those specific features from model matrices
    data_unknowns_temp<-as.data.frame(featmat_unknowns) %>%
      # select(all_of(predictive_features_temp)) %>%
      rownames_to_column(., var="assembly_id") %>%
      ## and add to pheno data
      full_join(glmnet_phenoresults_all_less, .) %>%
      select(-assembly_id)
    
    ## append to data
    ensemble_data<-bind_rows(ensemble_data, data_unknowns_temp)
    
  }
  
  ensemble_data<-ensemble_data %>%
    relocate(sort(names(.))) %>%
    relocate(all_of(c(starts_with("order"), starts_with("FP")))) %>%
    relocate(organism_name,
             model_number,
             log_predicted_age_mat,
             detrans_predicted_age_mat, 
             bagged_detrans_predicted_age_mat)
  
  ## export data
  write.csv(ensemble_data,
            paste0("dataFiles/13.00_ensemble_data_", mod_tl, ".csv"),
            row.names=FALSE)
  
  write.csv(bagged_phenoresults_unknowns, 
            paste0("dataFiles/13.00_bagged_phenoresults_unknowns_", mod_tl, ".csv"), 
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