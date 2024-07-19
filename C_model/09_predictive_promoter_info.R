## SET UP ----

## define data sets
models <- c("allhomo", "fishdanio", "mammalshomo", "reptilesgallus")

## for testing
# mod <- c("reptilesgallus")

## set working to script location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

start_time <- Sys.time()
print(paste("Script start time: ", start_time))

## load libraries
library(viridis) 
library(ggpubr) 
library(ggrepel)
library(plotmo)
library(tidyverse) 

for (mod in models){
  print("################################################################")
  print(paste0("########################### ", mod, " ###############################"))
  print("################################################################")
  
  ## lowercase
  mod_tl <- tolower(mod)
  
  ## IMPORT AND FORMAT DATA ----
  
  ## import and join bagged data
  test <- read.csv(paste0("dataFiles/07.00_bagged_phenoresults_test_", mod_tl, ".csv"))
  trainVal <- read.csv(paste0("dataFiles/07.00_bagged_phenoresults_", mod_tl, ".csv"))
  bagged_phenoresults_all <- bind_rows(test, trainVal)

  ## import and format weight information
  feature_weights <- read.csv(paste0("dataFiles/07.00_glmnet_weights_",
                                       mod_tl, ".csv")) %>%
    ## remove intercept information
    filter(feature_id != "(Intercept)") %>%
    ## make lexographic
    rename_with(stringr::str_replace,
                pattern = "partition_", replacement = "partition_0") %>%
    ## order columns
    select("feature_id", sort(colnames(.)))
  
  names(feature_weights) <- gsub("010", "10", names(feature_weights))
  
  ## get a vector of predictive features (from all final models)
  nonzero_features <- unique(feature_weights$feature_id)
  
  ## get those that are promoters
  nonzero_promoters <- nonzero_features[grepl("^FP", nonzero_features)]
  
  ## get those that are non promoters
  nonzero_otherfeats <- nonzero_features[!grepl("^FP", nonzero_features)]
  
  ## import model data
  featmat_trainval <- read.csv(paste0("dataFiles/07.00_featmat_trainval_", 
                                       mod_tl, ".csv")) %>%
    column_to_rownames(., var="X") %>%
    ## format as matrix
    as.matrix(.)
    
  ## import CpG density data (keep NAs)
  cpg_oe_data_full <- read.csv(paste0("../B_exploreData/dataFiles/05.00_cpg_oe_data_", 
                                      mod_tl, ".csv")) %>%
    column_to_rownames(., var="X") %>%
    ## format as matrix
    as.matrix(.)
  
  ## keep only features in the models (matrix)
  cpg_oe_data_pp <- cpg_oe_data_full[ , colnames(cpg_oe_data_full) %in% nonzero_features]
  
  ## edit cpg data (keep NAs)
  cpg_oe_data_pp_r_nz <- as.data.frame(cpg_oe_data_full) %>%
    ## move the rownames back to a column
    rownames_to_column(., var = "organism_name") %>%
    ## add the response info from the phenodata
    left_join(bagged_phenoresults_all %>% 
                select(organism_name, mean_age_mat) %>%
                unique(.)) %>%
    ## select only predictive promoters
    select(organism_name, mean_age_mat, 
           all_of(nonzero_promoters))
  
  ## import cpg data (replace NAs with zeros)
  cpg_oe_data_full_z <-
    read.csv(paste0("../B_exploreData/dataFiles/05.00_cpg_oe_data_", mod_tl, ".csv")) %>%
    ##  move the species (col "X") names to row names
    column_to_rownames(., var="X") %>%
    ## replace NAs with zeros
    replace(is.na(.), 0) %>%
    ## format as matrix
    as.matrix(.)
  
  ## keep only features in the models (matrix)
  cpg_oe_data_pp_z <- cpg_oe_data_full_z[ , colnames(cpg_oe_data_full_z) %in% nonzero_features]
  
  ## edit cpg oe data
  cpg_oe_data_pp_r_z <- as.data.frame(cpg_oe_data_full_z) %>%
    ## move the rownames back to a column
    rownames_to_column(., var = "organism_name") %>%
    ## add the response info from the phenodata
    left_join(bagged_phenoresults_all %>% 
                select(organism_name, mean_age_mat) %>%
                unique(.)) %>%
    ## select only predictive promoters
    select(organism_name, mean_age_mat, 
           all_of(nonzero_promoters))

  ## CALCULATE PROMOTER PEARSON CORRELATIONS ----
  
  ## including zero values ----
  
  ## create empty df for correlation data
  corr_data <- NULL
  
  ## loop through list of promoters and get correlations
  for (promoter in nonzero_promoters) {
    temp_corr <- cor.test(cpg_oe_data_pp_r_z[[promoter]], 
                          cpg_oe_data_pp_r_z$mean_age_mat, 
                          method = "pearson")
    temp_cor <- temp_corr$estimate
    temp_p <- temp_corr$p.value
    temp_row <- c(promoter, temp_cor, temp_p)
    corr_data <- rbind.data.frame(corr_data, temp_row, 
                                  make.row.names = FALSE)
  }
  
  ## name the columns
  names(corr_data) <- c("feature_id", "corr", "p_value")
  
  ## edit the data frame
  corr_data <- corr_data %>%
    ## convert from character to numeric
    mutate(corr=as.numeric(corr), 
           p_value=as.numeric(p_value)) %>%
    ## annotate for sig cors only
    mutate(sig_corr = ifelse(p_value > 0.05, 
                             NA, 
                             corr))

  ## not including zero values ----
  
  corr_data_nz <- NULL
  ## loop through list of promoters and get correlations
  for (promoter in nonzero_promoters) {
    temp_corr <- tryCatch(
      cor.test(cpg_oe_data_pp_r_nz[[promoter]], 
               cpg_oe_data_pp_r_nz$mean_age_mat, method = "pearson"),
      error = function(e) e
    )
    
    if (inherits(temp_corr, "error")) {
      # Error occurred, print a message and continue to the next iteration
      message(paste("Error occurred for promoter:", promoter))
      next
    }
    
    temp_cor <- temp_corr$estimate
    temp_p <- temp_corr$p.value
    temp_row <- c(promoter, temp_cor, temp_p)
    corr_data_nz <- rbind.data.frame(corr_data_nz, temp_row, make.row.names = FALSE)
  }
  
  ## name the columns
  names(corr_data_nz) <- c("feature_id", "corr_nz", "p_value_nz")
  
  ## edit the data frame
  corr_data_nz <- corr_data_nz %>%
    ## convert from character to numeric
    mutate(corr_nz=as.numeric(corr_nz), 
           p_value_nz=as.numeric(p_value_nz)) %>%
    ## annotate for sig cors only
    mutate(sig_corr_nz = ifelse(p_value_nz > 0.05, 
                                NA, 
                                corr_nz))
  
  ## CALCULATE MEAN/MEDIAN PROMOTER CPG O/E ----
  mean_cpg_oe_pp_nz <- colMeans(cpg_oe_data_pp, na.rm = TRUE)
  mean_cpg_oe_pp_z <- colMeans(cpg_oe_data_pp_z, na.rm = TRUE)
  
  mdn_cpg_oe_pp_nz <- apply(cpg_oe_data_pp, 2, median, na.rm = TRUE)
  mdn_cpg_oe_pp_z <- apply(cpg_oe_data_pp_z, 2, median, na.rm = TRUE)
  
  mean_cpg_oe <- data.frame(mean_cpg_oe_pp_nz, mean_cpg_oe_pp_z, 
                            mdn_cpg_oe_pp_nz, mdn_cpg_oe_pp_z) %>%
    rownames_to_column(var="feature_id")

  ## add corr data to promoter data ----
  feature_weights_long <- feature_weights %>%
    pivot_longer(., cols=starts_with("model"),
                 values_to = "weight",
                 names_to = "model_number", 
                 names_prefix = "model_") %>%
    group_by(feature_id) %>%
    mutate(number_of_models = sum(!is.na(weight)))
  
  ## average feature weights
  nonzero_features_avg <- feature_weights_long %>%
    # filter(feature_id %in% nonzero_promoters) %>%
    group_by(feature_id) %>%
    mutate(mean_weight = mean(weight, na.rm=TRUE)) %>%
    select(feature_id, mean_weight, number_of_models) %>%
    mutate(abs_mean_weight = abs(mean_weight)) %>%
    unique(.) %>%
    ## add correlation data (inc zeros)
    left_join(., corr_data) %>%
    ## add correlation data (not inc zeros)
    left_join(., corr_data_nz) %>%
    ## add cpg oe data
    left_join(., mean_cpg_oe) %>%
    ## if negative, positive or varied
    mutate(npv = ifelse(mean_weight < 0 & corr < 0,
                        "Negative",
                        ifelse(mean_weight > 0 & corr > 0,
                               "Positive",
                               "Varied"))) %>%
    ## if negative, positive or varied
    mutate(npv_nz = ifelse(mean_weight < 0 & corr_nz < 0,
                           "Negative",
                           ifelse(mean_weight > 0 & corr_nz > 0,
                                  "Positive",
                                  "Varied"))) %>%
    ## arrange by 
    arrange(mean_weight) %>%
    mutate(feature_id_f = fct_inorder(factor(feature_id))) %>%
    ## add significance
    mutate(sigyn = ifelse(is.na(sig_corr),
                          "No",
                          "Yes")) %>%
    ## add significance
    mutate(sigyn_nz = ifelse(is.na(sig_corr_nz),
                             "No",
                             "Yes"))
  
  ## export 
  write.csv(nonzero_features_avg, 
            paste0("dataFiles/09.00_nonzero_features_avg_", mod_tl, ".csv"), 
            row.names=FALSE)

  ## ADD ADDITIONAL INFO TO PHENODATA ----
  
  ## format hit length ----
  
  ## read in hit length data (all species)
  hit_len_data_full <- read.csv(paste0("../A_getGenomicData/dataFiles/04.00_hit_len_data_full_", 
                                       mod_tl, ".csv")) %>%
    ## move column X to rows
    column_to_rownames(., var="X") %>%
    ## format as matrix
    as.matrix(.)
  
  ## remove species not in the model (no age@maturity info)
  no_pheno_species <- setdiff(rownames(hit_len_data_full), 
                             unique(bagged_phenoresults_all$assembly_id))
  hit_len_data <- hit_len_data_full[!rownames(hit_len_data_full) %in% no_pheno_species, ]
  
  ## keep only promoters in the model
  hit_len_data_pp <- hit_len_data[ , colnames(hit_len_data) %in% nonzero_features]
  
  ## format percent identity ----
  
  ## read in hit percent identity data (all species)
  hit_pi_data_full <- read.csv(paste0("../A_getGenomicData/dataFiles/04.00_hit_pi_data_full_", 
                                      mod_tl, ".csv")) %>%
    ## move column X to rows
    column_to_rownames(., var="X") %>%
    ## format as matrix
    as.matrix(.)
  
  ## remove species not in the model
  no_pheno_species <- setdiff(rownames(hit_pi_data_full), 
                             unique(bagged_phenoresults_all$assembly_id))
  hit_pi_data <- hit_pi_data_full[!rownames(hit_pi_data_full) %in% no_pheno_species, ]
  
  ## keep only promoters in the model
  hit_pi_data_pp <- hit_pi_data[ , colnames(hit_pi_data) %in% nonzero_features]
  
  ## add all to phenodata ----
  bagged_phenoresults_all <- bagged_phenoresults_all %>%
    ## add mean cpg oe (NAs)
    left_join(., (as.data.frame(cpg_oe_data_pp) %>%
                    mutate(mean_cpg_oe_pp_nz = rowMeans(., na.rm=TRUE)) %>%
                    ## also count total num hits
                    mutate(num_hits_pp = rowSums(!is.na(select(., -mean_cpg_oe_pp_nz)))) %>% 
                    select(mean_cpg_oe_pp_nz, num_hits_pp) %>%
                    rownames_to_column(., var = "organism_name")),
              by = "organism_name") %>%
    ## add mean cpg oe (zeros)
    left_join(., (as.data.frame(cpg_oe_data_pp)) %>%
                ## replace NAs
                replace(is.na(.), 0) %>%
                mutate(mean_cpg_oe_pp_z = rowMeans(.)) %>%
                select(mean_cpg_oe_pp_z) %>%
                rownames_to_column(., var = "organism_name"),
              by = "organism_name") %>%
    ## add mean hit length (NAs)
    left_join(., (as.data.frame(hit_len_data_pp) %>%
                    mutate(mean_hit_len_pp_nz = rowMeans(., na.rm=TRUE)) %>%
                    select(mean_hit_len_pp_nz) %>%
                    rownames_to_column(., var = "assembly_id")),
              by = "assembly_id") %>%
    ## add mean hit length (zeros)
    left_join(., (as.data.frame(hit_len_data_pp)) %>%
                ## replace NAs
                replace(is.na(.), 0) %>%
                mutate(mean_hit_len_pp_z = rowMeans(.)) %>%
                select(mean_hit_len_pp_z) %>%
                rownames_to_column(., var = "assembly_id"),
              by = "assembly_id") %>%
    ## add mean percent ident length (NAs)
    left_join(., (as.data.frame(hit_pi_data_pp) %>%
                    mutate(mean_hit_pi_pp_nz = rowMeans(., na.rm=TRUE)) %>%
                    select(mean_hit_pi_pp_nz) %>%
                    rownames_to_column(., var = "assembly_id")),
              by = "assembly_id") %>%
    ## add mean percent ident length (zeros)
    left_join(., (as.data.frame(hit_pi_data_pp)) %>%
                ## replace NAs
                replace(is.na(.), 0) %>%
                mutate(mean_hit_pi_pp_z = rowMeans(.)) %>%
                select(mean_hit_pi_pp_z) %>%
                rownames_to_column(., var = "assembly_id"),
              by = "assembly_id")

  write.csv(bagged_phenoresults_all, 
            paste0("dataFiles/09.00_bagged_phenoresults_", mod_tl, ".csv"), 
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
