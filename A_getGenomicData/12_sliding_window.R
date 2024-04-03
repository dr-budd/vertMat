## SET UP ----

## load libraries
library(Biostrings)
library(viridis) 
library(ggpubr)
library(ggtext)
library(tidyverse) 

## for figure exports
cpg_type <- "O/E"
cpg_type_less <- ifelse(cpg_type == "O/E", "oe", cpg_type)

## print it
paste0("CpG type specified: ", cpg_type)

## species included
species <- c("chicken", "human", "zebrafish")

## load sequences
sequences_human <- readDNAStringSet("promotersExtend/epd_homo_sapiens_most_rep_extend.fa")
sequences_chicken <- readDNAStringSet("promotersExtend/epd_gallus_gallus_most_rep_extend.fa")
sequences_zebrafish <- readDNAStringSet("promotersExtend/epd_danio_rerio_most_rep_extend.fa")

## set window and step size for sliding window analysis
window_size <- 100
step_size <- 10

## SLIDING WINDOW ANALYSIS ----

## 1st loop iterates through species

for (spp in species) {

  print(paste0(spp, " start time: ", Sys.time()))

  temp_name <- paste0("cpg_contents_", spp)
  temp_sequences <- get(paste0("sequences_", spp))
  
  temp_cpg_contents <- list()
  
  ## 2nd loop iterates through each sequence
  
  for (seq_idx in 1:length(temp_sequences)) {
    
    seq <- temp_sequences[[seq_idx]]
    seq_length <- length(seq)
    num_windows <- floor((seq_length - window_size) / step_size) + 1
    
    cpg_content_seq <- numeric(num_windows)
    
    ## 3rd loop iterates through windows within each sequences
    for (window_start in seq(1, by = step_size, length.out = num_windows)) {
      
      window_seq <- subseq(seq, 
                           start = window_start, 
                           end = window_start + window_size - 1)
      
      c_count <- countPattern("C", window_seq)
      g_count <- countPattern("G", window_seq)
      cg_count <- countPattern("CG", window_seq)
      
      density <- cg_count / window_size
      oe <- cg_count / ((c_count * g_count) / window_size)
      expected <- (c_count / window_size) * (g_count / window_size)
      
      ifelse(cpg_type == "density", 
             cpg_content <- density, 
             ifelse(cpg_type == "O/E", 
                    cpg_content <- oe, 
                    ifelse(cpg_type == "expected", 
                           cpg_content <- expected)))
      
      cpg_content_seq[(window_start - 1) / step_size + 1] <- cpg_content
    }
    
    temp_cpg_contents[[seq_idx]] <- cpg_content_seq
  }
  assign(temp_name, temp_cpg_contents)
}

## DATA HANDLING ----

## list to data frame ----

## loop iterates through species

for (spp in species) {
  
  temp_name <- paste0("df_cpg_contents_", spp)
  df_cpg_contents <- data.frame()
  
  temp_cpg_contents <- get(paste0("cpg_contents_", spp))
  temp_sequences <- get(paste0("sequences_", spp))
  
  for (seq_idx in 1:length(temp_sequences)) {
    seq_name <- names(temp_sequences)[seq_idx]
    cpg_content_seq <- temp_cpg_contents[[seq_idx]]
    
    df <- data.frame(Sequence = rep(seq_name, length(cpg_content_seq)),
                     Position = seq(1 - 4999, by = step_size, length.out = length(cpg_content_seq)),
                     CpG_content = cpg_content_seq)
    
    df_cpg_contents <- rbind(df_cpg_contents, df)
  }
  assign(temp_name, df_cpg_contents)
}

## calculate means ----

for (spp in species) {
  temp_name <- paste0("mean_df_cpg_contents_", spp)
  temp_df <- get(paste0("df_cpg_contents_", spp))
  
  temp_mean_df <- temp_df %>%
    group_by(Position) %>%
    summarise(mean_cpg_content=mean(CpG_content, na.rm=TRUE))
  
  assign(temp_name, temp_mean_df)
}

## combine data frames ----
combined <- full_join(mean_df_cpg_contents_chicken %>%
                        rename(chicken = mean_cpg_content), 
                      mean_df_cpg_contents_human %>%
                        rename(human = mean_cpg_content)) %>%
  full_join(., mean_df_cpg_contents_zebrafish %>%
              rename(zebrafish = mean_cpg_content)) %>%
  pivot_longer(., cols=-Position, names_to="species", values_to="mean_cpg_content")

## PLOT DATA ----

## 2 = danio, 4 = chicken, 3 = human
custom_viridis <- c(viridis_pal(option = "D", direction=-1)(4)[4], 
                    viridis_pal(option = "D", direction=-1)(4)[3], 
                    viridis_pal(option = "D", direction=-1)(4)[2])

p <- ggplot(combined %>% 
              filter(species != "mouse") %>%
              mutate(Species = gsub("chicken", "Chicken", species) %>%
                       gsub("human", "Human", .) %>%
                       gsub("zebrafish", "Zebrafish", .)), 
            aes(x = Position, y = mean_cpg_content)) +
  geom_line(aes(colour=Species)) +
  scale_color_manual(values=custom_viridis)+
  labs(x = "Distance from TSS", 
       y = "Mean CpG O/E")+
  theme_classic()+
  geom_vline(xintercept=-100, lty=2, colour="grey10", linewidth=0.1)+
  geom_vline(xintercept=100, lty=2, colour="grey10", linewidth=0.1)

# p

## save plot
saveRDS(p, "figuresTables/12.00_sliding_window.rds")
ggsave("figuresTables/12.00_sliding_window.png", p)

## end script
