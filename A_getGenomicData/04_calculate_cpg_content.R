## SET UP ----

## define data sets
data_sets <- c("AllHomo",
               "FishDanio",
               "MammalsHomo",
               "ReptilesGallus")

## if troubleshooting
# data_set <- c("AllHomo")

start_time <- Sys.time()
print(paste("Script start time: ", start_time))

## load libraries
library(tidyverse)

## busco scores ----

## import busco summary lines
busco_suml <- read.delim("dataFiles/05_busco_summary_lines.txt", 
                         header=FALSE, sep=" ")

## edit to get scores
busco_scores <-  busco_suml %>%
  rename(assembly_id = V1, busco_info = V2) %>%
  ## remove missing info
  filter(., busco_info!="") %>%
  separate(busco_info, c("b_complete", "b_complete_and_duplicated", 
                         "b_fragmented", "b_missing", 
                         "b_total_groups_searched"),
           sep=",") %>%
  separate(b_complete, c("b_complete", "b_complete_and_single_copy"), 
           sep="\\[") %>%
  mutate(b_complete = as.numeric(regmatches(b_complete, 
                                            gregexpr("[[:digit:]]+\\.*[[:digit:]]*",
                                                     b_complete))), 
         b_complete_and_single_copy = as.numeric(regmatches(b_complete_and_single_copy, 
                                                            gregexpr("[[:digit:]]+\\.*[[:digit:]]*",
                                                                     b_complete_and_single_copy))), 
         b_complete_and_duplicated = as.numeric(regmatches(b_complete_and_duplicated, 
                                                           gregexpr("[[:digit:]]+\\.*[[:digit:]]*",
                                                                    b_complete_and_duplicated))),
         b_fragmented = as.numeric(regmatches(b_fragmented, 
                                              gregexpr("[[:digit:]]+\\.*[[:digit:]]*",
                                                       b_fragmented))),
         b_missing = as.numeric(regmatches(b_missing, 
                                           gregexpr("[[:digit:]]+\\.*[[:digit:]]*",
                                                    b_missing))),
         b_total_groups_searched = as.numeric(regmatches(b_total_groups_searched, 
                                                         gregexpr("[[:digit:]]+\\.*[[:digit:]]*",
                                                                  b_total_groups_searched))))

## CALCULATE CPG CONTENT ----

## create file for no-hit species
no_hits_all <- NULL

## loop through analysis
for (data_set in data_sets) {
  print("################################################################")
  print(paste0("########################### ", data_set, " ###############################"))
  print("################################################################")
  
  ## print time
  Sys.time()
  
  ## lowercase
  data_set_tl <- tolower(data_set)
  
  ## set working directly to query results directory
  setwd(paste0("QueryResults", data_set))
  
  ## create a list of blast hit files to import into R
  data_files_all <- list.files(pattern = "top_hits")
  
  ## remove any files of size zero from file list
  data_files_empty <- NULL
  for (i in data_files_all) {
    if (file.size(i) == 0)
      data_files_empty <- c(data_files_empty, i)
  }
  
  print("Empty data files: ")
  print(data_files_empty)
  
  empty_assembly_ids <- gsub("_genomic_top_hits", "", data_files_empty)
  no_hits <- data.frame(assembly_id = empty_assembly_ids) %>%
    mutate(reference_spp = gsub("allhomo", "Homo sapien", data_set_tl) %>%
             gsub("mammalshomo", "Homo sapien", .) %>%
             gsub("reptilesgallus", "Gallus gallus", .) %>%
             gsub("fishdanio", "Danio rerio", .))
  
  no_hits_all <- rbind(no_hits_all, no_hits)
  
  data_files <- NULL
  for (i in data_files_all) {
    if (!file.size(i) == 0)
      data_files <- c(data_files, i)
  }
  
  ## create vector of assembly names
  assembly_names <- gsub("_genomic_top_hits", "", data_files)
  
  ## for first file in dir ----
  
  ## specify the file number
  file_num <- 1
  
  ## save the assembly name to a variable (the name of the file minus the extension)
  temp_assembly <- paste(gsub("_genomic_top_hits", "", data_files[file_num]))
  
  ## print the file number and assembly name
  print(paste("file number:", file_num, "of", length(data_files)))
  
  ## import and format data, then calculate CG density
  temp_data <- read.table(file = data_files[file_num]) %>%
    ## add column names (from blast query call)
    ## NB: this is necessary but nice if you want to eyeball the data
    rename(query_seqid = V1, subject_seqid = V2, subject_seq_len = V3, 
           query_start = V4, query_end = V5, ## start and end of alignment in query
           subject_start = V6, subject_end = V7, ## start and end of alignment in subject
           perc_ident = V8, num_ident = V9, ## percent and number of identical matches
           subject_seq = V10) %>% ## aligned part of subject sequence 
    ## remove sequence gap indicators from alignment ("-")
    mutate(subject_seq = gsub("-", "", subject_seq), 
           ## remove sequence gap indicators from alignment ("-")
           length = nchar(subject_seq), 
           ## count the number of "CG" occurrences
           cg_count = str_count(subject_seq, "CG"), 
           ## count the number of Cs
           c_count = str_count(subject_seq, "C"),
           ## count the number of Gs
           g_count = str_count(subject_seq, "G"),
           ## calculate the CpG E (expected)
           cg_e = (c_count/length)*(g_count/length),
           ## calculate the CpG density (observed)
           cg_dens = cg_count/length,
           ## calculate the CpG O/E 
           cg_oe = cg_count/((c_count*g_count)/length)) %>%
    ## rename the promoter ID column
    rename(., promoter_id = query_seqid)
  
  ## remove duplicate hits (retain only 'top' hit)
  temp_data <- temp_data[!duplicated(temp_data[,1]), ]
  
  ## edit and rename the OE data frame 
  temp_oe_data <- temp_data %>%
    ## retain only the promoter id and CpG densities
    select(promoter_id, cg_oe) %>%
    ## rename the cg_oe column with the assembly name
    rename(., "{temp_assembly}" := cg_oe) %>%
    ## reorder by promoter id
    arrange(., promoter_id)
  
  ## edit and rename the density data frame 
  temp_dens_data <- temp_data %>%
    ## retain only the promoter id and CpG densities
    select(promoter_id, cg_dens) %>%
    ## rename the cg_oe column with the assembly name
    rename(., "{temp_assembly}" := cg_dens) %>%
    ## reorder by promoter id
    arrange(., promoter_id)
  
  ## edit and rename the expected data frame 
  temp_expect_data <- temp_data %>%
    ## retain only the promoter id and CpG densities
    select(promoter_id, cg_e) %>%
    ## rename the cg_oe column with the assembly name
    rename(., "{temp_assembly}" := cg_e) %>%
    ## reorder by promoter id
    arrange(., promoter_id)
  
  ## edit and rename the length data frame 
  temp_len_data <- temp_data %>%
    ## retain only the promoter id and hit lengths
    select(promoter_id, length) %>%
    ## rename the cg_oe column with the assembly name
    rename(., "{temp_assembly}" := length) %>%
    ## reorder by promoter id (optional here)
    arrange(., promoter_id)
  
  ## edit and rename the percent id data frame 
  temp_pi_data <- temp_data %>%
    ## retain only the promoter id and percent idents
    select(promoter_id, perc_ident) %>%
    ## rename the cg_oe column with the assembly name
    rename(., "{temp_assembly}" := perc_ident) %>%
    ## reorder by promoter id (optional here)
    arrange(., promoter_id)
  
  ## for remaining files in dir ----
  
  ## loop through each file 
  for(file_num in 2:length(data_files)){
    
    ## save the assembly name to a variable (the name of the file minus the extension)
    temp_assembly <- paste(gsub("_genomic_top_hits", "", data_files[file_num]))
    
    ## print the file number and assembly name
    print(paste("file number:", file_num, "of", length(data_files)))
    
    ## import and format data, then calculate CG density
    temp_data <- read.table(file = data_files[file_num]) %>%
      ## add column names (from blast query call)
      ## NB: this is not necessary but nice if you want to eyeball the data
      rename(query_seqid = V1, subject_seqid = V2, subject_seq_len = V3, 
             query_start = V4, query_end = V5, ## start and end of alignment in query
             subject_start = V6, subject_end = V7, ## start and end of alignment in subject
             perc_ident = V8, num_ident = V9, ## percent and number of identical matches
             subject_seq = V10) %>% ## aligned part of subject sequence 
      ## remove sequence gap indicators from alignment ("-")
      mutate(subject_seq = gsub("-", "", subject_seq), 
             ## remove sequence gap indicators from alignment ("-")
             length = nchar(subject_seq), 
             ## count the number of "CG" occurrences
             cg_count = str_count(subject_seq, "CG"), 
             ## count the number of Cs
             c_count = str_count(subject_seq, "C"),
             ## count the number of Gs
             g_count = str_count(subject_seq, "G"),
             ## calculate the CpG E (expected)
             cg_e = (c_count/length)*(g_count/length),
             ## calculate the CpG density (observed)
             cg_dens = cg_count/length,
             ## calculate the CpG O/E 
             cg_oe = cg_count/((c_count*g_count)/length)) %>%
      ## rename the promoter ID column
      rename(., promoter_id = query_seqid)
    
    ## remove duplicate hits (retain only 'top' hit)
    temp_data <- temp_data[!duplicated(temp_data[,1]), ]
    
    ## edit and rename the OE data frame 
    temp_oe_data <- temp_data %>%
      ## retain only the promoter id and CpG densities
      select(promoter_id, cg_oe) %>%
      ## rename the cg_oe column with the assembly name
      rename(., "{temp_assembly}" := cg_oe) %>%
      ## reorder by promoter id
      arrange(., promoter_id) %>%
      ## join temp data to existing data (and overwrite it)
      full_join(temp_oe_data, ., by = "promoter_id") 
    
    ## edit and rename the density data frame 
    temp_dens_data <- temp_data %>%
      ## retain only the promoter id and CpG densities
      select(promoter_id, cg_dens) %>%
      ## rename the cg_oe column with the assembly name
      rename(., "{temp_assembly}" := cg_dens) %>%
      ## reorder by promoter id
      arrange(., promoter_id) %>%
      ## join temp data to existing data (and overwrite it)
      full_join(temp_dens_data, ., by = "promoter_id") 
    
    ## edit and rename the expected data frame 
    temp_expect_data <- temp_data %>%
      ## retain only the promoter id and CpG densities
      select(promoter_id, cg_e) %>%
      ## rename the cg_oe column with the assembly name
      rename(., "{temp_assembly}" := cg_e) %>%
      ## reorder by promoter id
      arrange(., promoter_id) %>%
      ## join temp data to existing data (and overwrite it)
      full_join(temp_expect_data, ., by = "promoter_id") 
    
    ## join the data frames (hit length)
    temp_len_data <- temp_data %>%
      ## retain only the promoter id and hit lengths
      select(promoter_id, length) %>%
      ## rename the cg_oe column with the assembly name
      rename(., "{temp_assembly}" := length) %>%
      ## reorder by promoter id (optional here)
      arrange(., promoter_id) %>%
      ## join temp data to existing data (and overwrite it)
      full_join(temp_len_data, ., by = "promoter_id") 
    
    ## edit and rename the percent id data frame 
    temp_pi_data <- temp_data %>%
      ## retain only the promoter id and percent idents
      select(promoter_id, perc_ident) %>%
      ## rename the cg_oe column with the assembly name
      rename(., "{temp_assembly}" := perc_ident) %>%
      ## reorder by promoter id (optional here)
      arrange(., promoter_id) %>%
      ## join
      full_join(temp_pi_data, ., by="promoter_id")
    
  }
  
  ## format results ----
  
  ## reset working to script location
  setwd("../")
  
  ## format cpg oe data  
  cpg_oe_data_full <- temp_oe_data %>%
    ## order the dataframe by promoter id 
    arrange(., promoter_id) %>%
    ## transpose CpG oe data for elastic net
    ## (this is round about but prevents num to chr)
    pivot_longer(-promoter_id) %>%
    pivot_wider(names_from = promoter_id, values_from = value) %>%
    ## move the assembly names to rownames
    column_to_rownames(., "name") %>%
    ## convert to matrix
    as.matrix(.)
  
  ## export
  write.csv(cpg_oe_data_full,
            paste0("dataFiles/04.00_cpg_oe_data_full_", data_set_tl, ".csv"),
                  row.names = TRUE)
  
  ## format cpg density data
  cpg_dens_data_full <- temp_dens_data %>%
    ## order the dataframe by promoter id
    arrange(., promoter_id) %>%
    ## transpose CpG oe data for elastic net
    ## (this is round about but prevents num to chr)
    pivot_longer(-promoter_id) %>%
    pivot_wider(names_from = promoter_id, values_from = value) %>%
    ## move the assembly names to rownames
    column_to_rownames(., "name") %>%
    ## convert to matrix
    as.matrix(.)

  ## export
  write.csv(cpg_dens_data_full,
            paste0("dataFiles/04.00_cpg_dens_data_full_", data_set_tl, ".csv"),
            row.names = TRUE)
  
  ## format cpg expected data
  cpg_expect_data_full <- temp_expect_data %>%
    ## order the dataframe by promoter id
    arrange(., promoter_id) %>%
    ## transpose CpG oe data for elastic net
    ## (this is round about but prevents num to chr)
    pivot_longer(-promoter_id) %>%
    pivot_wider(names_from = promoter_id, values_from = value) %>%
    ## move the assembly names to rownames
    column_to_rownames(., "name") %>%
    ## convert to matrix
    as.matrix(.)
  
  ## export
  write.csv(cpg_expect_data_full,
            paste0("dataFiles/04.00_cpg_expect_data_full_", data_set_tl, ".csv"),
            row.names = TRUE)
  
  ## format hit length data
  hit_len_data_full <- temp_len_data %>%
    ## order the dataframe by promoter id 
    arrange(., promoter_id) %>%
    ## transpose CpG oe data for elastic net
    ## (this is round about but prevents num to chr)
    pivot_longer(-promoter_id) %>%
    pivot_wider(names_from = promoter_id, values_from = value) %>%
    ## move the assembly names to rownames
    column_to_rownames(., "name") %>%
    ## convert to matrix
    as.matrix(.)
  
  ## export
  write.csv(hit_len_data_full, 
            paste0("dataFiles/04.00_hit_len_data_full_", data_set_tl, ".csv"), 
            row.names = TRUE)
  
  ## format percent identity data
  hit_pi_data_full <- temp_pi_data %>%
    ## order the dataframe by promoter id 
    arrange(., promoter_id) %>%
    ## transpose CpG oe data for elastic net
    ## (this is round about but prevents num to chr)
    pivot_longer(-promoter_id) %>%
    pivot_wider(names_from = promoter_id, values_from = value) %>%
    ## move the assembly names to rownames
    column_to_rownames(., "name") %>%
    ## convert to matrix
    as.matrix(.)
  
  ## export
  write.csv(hit_pi_data_full, 
            paste0("dataFiles/04.00_hit_pi_data_full_", data_set_tl, ".csv"), 
            row.names = TRUE)
  
  ## GENERATE ADDITIONAL PHENO DATA ----
  
  ## calculate blast hit info ----
  
  ## calculate mean hit length per species assembly
  nz_hit_len_all_chords <- temp_len_data %>%
    ## order the data frame by promoter id 
    arrange(., promoter_id) %>%
    ## format to long and get summary data
    pivot_longer(-promoter_id, names_to = "assembly_id") %>%
    group_by(assembly_id) %>%
    summarise(mean_nz_hit_len = mean(value, na.rm = TRUE))
  
  ## calculate mean percent identity per species assembly
  nz_pi_all_chords <- temp_pi_data %>%
    ## order the dataframe by promoter id 
    arrange(., promoter_id) %>%
    ## format to long and get summary data
    pivot_longer(-promoter_id, names_to = "assembly_id") %>%
    group_by(assembly_id) %>%
    summarise(mean_nz_perc_ident = mean(value, na.rm = TRUE))
  
  ## calculate mean hit length per species assembly
  hit_len_all_chords <- temp_len_data %>%
    ## replace NAs with zeros
    replace(is.na(.), 0) %>%
    ## order the dataframe by promoter id
    arrange(., promoter_id) %>%
    ## format to long and get summary data
    pivot_longer(-promoter_id, names_to = "assembly_id") %>%
    group_by(assembly_id) %>%
    summarise(mean_hit_len = mean(value, na.rm = TRUE))
  
  ## calculate mean hit length per species assembly
  pi_all_chords <- temp_pi_data %>%
    ## replace NAs with zeros
    replace(is.na(.), 0) %>%
    ## order the dataframe by promoter id
    arrange(., promoter_id) %>%
    ## format to long and get summary data
    pivot_longer(-promoter_id, names_to = "assembly_id") %>%
    group_by(assembly_id) %>%
    summarise(mean_perc_ident = mean(value, na.rm = TRUE))
  
  ## calculate number of hits per species assembly 
  hit_num_all_chords <- as.data.frame(cpg_oe_data_full) %>% 
    mutate(num_hits = rowSums(!is.na(.))) %>%
    select(num_hits) %>%
    rownames_to_column(., var = "assembly_id")
  
  ## create genomic pheno data ----
  
  ## non-zero hit length
  pheno_genomic <- nz_hit_len_all_chords %>%
    ## add hit length
    left_join(., hit_len_all_chords, by = "assembly_id") %>%
    ## add num hits
    left_join(., hit_num_all_chords, by = "assembly_id") %>%
    ## add perc id
    left_join(., pi_all_chords, by = "assembly_id") %>%
    ## add non-zero perc id
    left_join(., nz_pi_all_chords, by = "assembly_id") %>%
    ## order by assembly id (to match cpg oe data)
    arrange(., assembly_id) %>%
    ## add mean cpg oe
    left_join(., (as.data.frame(cpg_oe_data_full) %>%
                    ## replace NAs with zeros
                    replace(is.na(.), 0)) %>%
                mutate(mean_cpg_oe = rowMeans(.)) %>%
                select(mean_cpg_oe) %>%
                rownames_to_column(., var = "assembly_id"),
              by = "assembly_id") %>%
    ## add mean non-zero cpg oe
    left_join(., (as.data.frame(cpg_oe_data_full)) %>%
                mutate(mean_nz_cpg_oe = rowMeans(., na.rm=TRUE)) %>%
                select(mean_nz_cpg_oe) %>%
                rownames_to_column(., var = "assembly_id"),
              by = "assembly_id") %>%
    ## add mean cpg dens
    left_join(., (as.data.frame(cpg_dens_data_full) %>%
                    ## replace NAs with zeros
                    replace(is.na(.), 0)) %>%
                mutate(mean_cpg_dens = rowMeans(.)) %>%
                select(mean_cpg_dens) %>%
                rownames_to_column(., var = "assembly_id"),
              by = "assembly_id") %>%
    ## add mean non-zero cpg dens
    left_join(., (as.data.frame(cpg_dens_data_full)) %>%
                mutate(mean_nz_cpg_dens = rowMeans(., na.rm=TRUE)) %>%
                select(mean_nz_cpg_dens) %>%
                rownames_to_column(., var = "assembly_id"),
              by = "assembly_id") %>%
    ## add mean cpg expect
    left_join(., (as.data.frame(cpg_expect_data_full) %>%
                    ## replace NAs with zeros
                    replace(is.na(.), 0)) %>%
                mutate(mean_cpg_expect = rowMeans(.)) %>%
                select(mean_cpg_expect) %>%
                rownames_to_column(., var = "assembly_id"),
              by = "assembly_id") %>%
    ## add mean non-zero cpg expect
    left_join(., (as.data.frame(cpg_expect_data_full)) %>%
                mutate(mean_nz_cpg_expect = rowMeans(., na.rm=TRUE)) %>%
                select(mean_nz_cpg_expect) %>%
                rownames_to_column(., var = "assembly_id"),
              by = "assembly_id") %>%
    ## add busco scores
    left_join(., busco_scores, by="assembly_id")
  
  ## export pheno genomic data
  write.csv(pheno_genomic, 
            paste0("dataFiles/04.00_pheno_genomic_", data_set_tl, ".csv"),
            row.names = FALSE)
  
}

## SUPP INFO ----

## create table of no hits for supp ----
no_hits_all_edit <- no_hits_all %>%
  rename(`Reference species` = reference_spp,
         `Assembly ID` = assembly_id) %>%
  unique(.) 

write.csv(no_hits_all_edit, "figuresTables/Table S[no_hits_all].csv", row.names = FALSE)

## print session info
sessionInfo()

end_time <- Sys.time()
print(paste("Script end time: ", end_time))
print(paste("Elapsed time: ",  end_time - start_time))

## END SCRIPT
