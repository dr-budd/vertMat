## NOTES ----

## please specify location and samples before running (means script can be run locally)

## SET UP ----

## set location and figure size

## local
location <- "local"
# samples <- "subset"

## remote
# location <- "remote"
samples <- "all"

## print info
paste0("location: ", location)
# paste0("samples: ", samples)

## set working directory 
if (location == "local") {
 setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
paste0("working directory: ", getwd())

## load libraries
print("~~~ load libraries ~~~")
library(ape)
library(ggtree)
library(ggnewscale)
library(viridis)
library(ggtreeExtra)
library(ggimage)
library(multcomp)
library(ggpubr)
library(ggforce)
library(ggtext)
library(tidyverse)

## IMPORT AND FORMAT DATA ----
print("~~~ import and format data ~~~")

## define models
models <- c("allhomo", "fishdanio", "mammalshomo", "reptilesgallus")

## import phenodata (don't use 09 because the order is not the same)
for (mod in models){
  
  if (location == "local") {
    phenodata_temp <- read.csv(paste0("../B_exploreData/dataFiles/05.00_pheno_", mod, ".csv"))
  } else {
    phenodata_temp <- read.csv(paste0("dataFiles/05.00_pheno_", mod, ".csv"))
  }
  
  assign(paste0("phenodata_unmod_", mod), phenodata_temp)
}

## create columns for mean cpg oe all-vert vs group specific 
phenodata_av <- phenodata_unmod_allhomo %>%
  rename(mean_cpg_oe_av = mean_cpg_oe, 
         mean_nz_cpg_oe_av = mean_nz_cpg_oe)

phenodata_gs <- bind_rows(phenodata_unmod_fishdanio, phenodata_unmod_mammalshomo) %>%
  bind_rows(., phenodata_unmod_reptilesgallus) %>%
  rename(mean_cpg_oe_gs = mean_cpg_oe, 
         mean_nz_cpg_oe_gs = mean_nz_cpg_oe)

phenodata_mod <- full_join(phenodata_av, 
                       phenodata_gs %>% 
                         select(organism_name, mean_cpg_oe_gs, mean_nz_cpg_oe_gs))

## subset (or not) phenodata
if (samples == "all") {
 phenodata <- phenodata_mod
} else {
 ## subset of the data if running test
 set.seed(1234)
 proportion_of_total <- 0.1 
 phenodata <- phenodata_mod %>%
  ## split df by group
  group_split(model_group) %>%
  ## sample 10% of the df for each group, ensuring at least 2 rows per group
  lapply(function(group_df) {
   if (nrow(group_df) >= 20) {
    return(sample_frac(group_df, proportion_of_total))
   } else {
    return(group_df[1:2, ])
   }
  }) %>%
  ## join them
  bind_rows()
}

## get timetree data
if (location == "local") {
 chordate_tree <- read.tree("../B_exploreData/timeTree/chordates_species.nwk")
} else {
 chordate_tree <- read.tree("timeTree/chordates_species.nwk")
}

## edit tip labels
chordate_tree$tip.label <- gsub("_", " ", chordate_tree$tip.label)

## keep only intersecting tips
keep_tips <- intersect(chordate_tree$tip.label, phenodata$organism_name)

## subset tree for genome and lifespan info species
my_tree <- keep.tip(chordate_tree, keep_tips)

# ## print tree
# ggtree(my_tree)+
#   geom_tiplab(size = 0.5)

## remove chordate tree (it's big)
rm(chordate_tree)

## format phenodata
phenodata_matching <- phenodata %>%
 filter(organism_name %in% my_tree$tip.label) %>%
 relocate(organism_name) %>% ## IMPORTANT!
 group_by(order) %>%
 add_count(.)

## GET PHYLOPICS ----
print("~~~ get phylopics ~~~")

## get tip data in order
tip_data <- fortify(my_tree, layout = "inorder") %>%
 select(label, y) %>%
 ## remove numbered tip labels
 filter(!grepl("'", label)) %>%
 ## sort by tip number
 arrange(y) %>%
 ## add empty column for phylopics
 mutate(phylopic_uid = NA_character_) %>%
 ## do some renaming
 rename(tip_number = y, organism_name = label) %>%
 ## add label column
 mutate(tip_label = gsub(" ", "_", organism_name)) %>%
 ## add group
 left_join(phenodata %>%
            select(organism_name, model_group, mean_age_mat))

## create function to retrieve Phylopic ID with error handling
get_phylopic_id <- function(x) {
 tryCatch(phylopic_uid(x), error = function(e) NA_character_)
}

## set the increment and step variables
increment <- 0
step <- round(nrow(phenodata_matching)*0.025)

## add phylopics UID's every stepth row (starting at amphibians)
for (i in (which(tip_data$model_group == "Amphibians")[1]:nrow(tip_data))) {
 ## calculate the index for searching (n=step rows apart unless absent)
 search_index <- i + increment

 ## ensure we don't go beyond the end of the data frame
 if (search_index <= nrow(tip_data)) {

  ## get phylopic info
  phylopic_info <- get_phylopic_id(tip_data$tip_label[search_index])

  ## if an ID is found, add it and update the increment to skip 20 rows
  if (is.data.frame(phylopic_info) && nrow(phylopic_info) > 0) {
   tip_data$phylopic_uid[search_index] <- phylopic_info$uid
   increment <- increment + step
  } else {
   tip_data$phylopic_uid[search_index] <- NA_character_
  }
 }
}

## add phylopics UID's every stepth row (ending at amphibians)
increment <- 0 ## reset increment
for (i in (1:which(tip_data$model_group == "Amphibians")[1])) {
 ## calculate the index for searching (n=step rows apart unless absent)
 search_index <- i + increment
 
 ## ensure we don't go beyond the end of the data frame
 if (search_index <= which(tip_data$model_group == "Amphibians")[1]) {
  
  ## get phylopic info
  phylopic_info <- get_phylopic_id(tip_data$tip_label[search_index])
  
  ## if an ID is found, add it and update the increment to skip 20 rows
  if (is.data.frame(phylopic_info) && nrow(phylopic_info) > 0) {
   tip_data$phylopic_uid[search_index] <- phylopic_info$uid
   increment <- increment + step
  } else { 
   tip_data$phylopic_uid[search_index] <- NA_character_
  }
 }
}

## manually edit so	awkwardly sized phylopics are no longer in df
tip_data <- tip_data %>%
  mutate(phylopic_uid = ifelse(organism_name == "Gavialis gangeticus",
                               NA_character_,
                               phylopic_uid)) %>%
  mutate(phylopic_uid = ifelse(organism_name == "Crocodylus porosus",
                               get_phylopic_id("Crocodylus_porosus")$uid,
                               phylopic_uid))

## save tip data
write.csv(tip_data,
          paste0("dataFiles/07.00_phylopic_tipdata_", location, samples, ".csv"), 
          row.names=FALSE)

## add tip data to phenodata
phenodata_matching <- phenodata_matching %>%
  ## add tip data
  left_join(., tip_data) %>%
  ## relocate columns
  relocate(., organism_name, phylopic_uid, tip_number) %>%
  ## make order an ordered factor (ordered by group)
  mutate(order = factor(order, levels = unique(.$order[order(.$model_group)])))

## CREATE COLOUR SCALE FOR ORDER ----
print("~~~ create colour scale for order ~~~")

## define base colours
amph_baseclr <- viridis_pal(option="D", direction=-1)(4)[1]
fish_baseclr <- viridis_pal(option="D", direction=-1)(4)[2]
maml_baseclr <- viridis_pal(option="D", direction=-1)(4)[3]
rptl_baseclr <- viridis_pal(option="D", direction=-1)(4)[4]

## get number of orders
num_amph_ords <- phenodata_matching %>% filter(model_group == "Amphibians") %>%
 distinct(model_group, order) %>% nrow()
num_fish_ords <- phenodata_matching %>% filter(model_group == "Fish") %>%
 distinct(model_group, order) %>% nrow()
num_maml_ords <- phenodata_matching %>% filter(model_group == "Mammals") %>%
 distinct(model_group, order) %>% nrow()
num_rptl_ords <- phenodata_matching %>% filter(model_group == "Reptiles") %>%
 distinct(model_group, order) %>% nrow()

## create colours using colorRampPalette
amph_clrs <- colorRampPalette(c(colorspace::lighten(amph_baseclr, 0.45),
                                amph_baseclr,
                                colorspace::darken(amph_baseclr, 0.45)))(num_amph_ords)

fish_clrs <- colorRampPalette(c(colorspace::lighten(fish_baseclr, 0.95),
                                fish_baseclr,
                                colorspace::darken(fish_baseclr, 0.95)))(num_fish_ords)

maml_clrs <- colorRampPalette(c(colorspace::lighten(maml_baseclr, 0.95),
                                maml_baseclr,
                                colorspace::darken(maml_baseclr, 0.95)))(num_maml_ords)

rptl_clrs <- colorRampPalette(c(colorspace::lighten(rptl_baseclr, 0.95),
                                rptl_baseclr,
                                colorspace::darken(rptl_baseclr, 0.95)))(num_rptl_ords)

# ## view colours
# scales::show_col(amph_clrs)

## order df
dorder <- phenodata_matching %>%
 group_by(order) %>%
 mutate(node = ifelse(all(is.null(getMRCA(my_tree, organism_name))), 
                      which(my_tree$tip.label == organism_name), 
                      getMRCA(my_tree, organism_name)))  %>%
 select(order, node, organism_name) %>%
 unique(.)

## PLOT TREES ----
print("~~~ plot tree ~~~")

## order tree ----
c1 <- ggtree(my_tree, size=0.15, layout = "circular")+
 geom_hilight(data=dorder, aes(node=node, fill=order), 
              type = "rect", alpha=0.8)+
 scale_fill_manual(name="Order",
                   values=c(amph_clrs, 
                            fish_clrs,
                            maml_clrs,
                            rptl_clrs),
                   guide=guide_legend(keywidth=0.5,
                                      keyheight=0.2,
                                      ncol=9,
                                      direction="vertical",
                                      order=1))+
 ## put tree on top of highlighting
 geom_tree(size=0.15)+
 theme(legend.position = 'bottom')

# c1

## get legend only
order_legend <- get_legend(c1)
order_legend_gg <- as_ggplot(order_legend)
# order_legend_gg

## add formatted phenodata to tree
c1.1 <- c1 %<+% phenodata_matching
# c1.1

cpg_oe_df <- phenodata_matching %>%
  ungroup(.) %>%
  select(organism_name, 
         # mean_cpg_oe_av, 
         mean_nz_cpg_oe_av, 
         # mean_cpg_oe_gs, 
         mean_nz_cpg_oe_gs) %>%
  pivot_longer(., cols = -organism_name, 
               names_to = "cpg_measure", 
               values_to = "cpg_oe")

## cpg density + pheno tree ----
## (heatmap and barplot)
c1.2 <- c1.1 + guides(fill = "none")+
 ## allow new fill scale
 new_scale_fill() +
 ## add promoter CG density
 geom_fruit(data=cpg_oe_df,
            geom=geom_tile,
            mapping=aes(y=organism_name,
                        x=cpg_measure,
                        fill=cpg_oe),
            offset = 0.03,
            pwidth = 0.07)+
 scale_fill_viridis_c(option = "rocket",
                      na.value="transparent",
                      direction=1,
                      "Mean\npromoter\nCpG O/E", 
                      guide=guide_legend(direction="vertical"))+
 ## allow new fill scale
 new_scale_fill() +
 ## add lifespan bars
 geom_fruit(geom=geom_bar,
            mapping=aes(y=organism_name,
                        fill=model_group,
                        x=mean_age_mat),
            pwidth=0.4,
            offset=0.06,
            orientation="y",
            stat="identity")+
 scale_fill_viridis_d(option="D", direction=-1, guide="none")+
  theme(legend.position = 'bottom')

# c1.2

## get legend only
cpg_legend <- get_legend(c1.2)
cpg_legend_gg <- as_ggplot(cpg_legend)
# cpg_legend_gg

## phylo pics tree ----
c1.3 <- c1.2+  
  new_scale_fill() +
  new_scale_colour() +
  geom_tiplab2(aes(image=phylopic_uid,
                   colour=model_group),
               size=0.045,
               geom="phylopic",
               offset=250,
               angle=0)+
  scale_colour_viridis_d(option="D", direction=-1, guide="none")+
  theme(legend.position="none")

## ADDITIONAL PLOTS ----
print("~~~ additional plots ~~~")

## read in sliding window 
sliding_window <- readRDS("../A_getGenomicData/figuresTables/12.00_sliding_window.rds")

## edit phenodata 
phenodata_less <- phenodata_mod %>%
  select(organism_name, model_group, mean_age_mat, gc_percent, 
         mean_nz_cpg_oe_av, mean_cpg_oe_av) %>%
  mutate(model_group = factor(model_group), 
         log_mean_age_mat = log(mean_age_mat)) %>%
  group_by(model_group) %>%
  mutate(mdnam=median(mean_age_mat), 
         mnam=mean(mean_age_mat), 
         mnlam=mean(log_mean_age_mat),
         mngc=mean(gc_percent), 
         mnoe=mean(mean_nz_cpg_oe_av)) %>%
  ungroup(.) %>%
  ## means (lower)
  mutate(label1_pos_am=(max(mean_age_mat)*1.1),
         label1_pos_lam=(max(mean_age_mat)*1.6), ## log
         label1_pos_gc=(max(gc_percent)*1.02),
         label1_pos_oe=(max(mean_nz_cpg_oe_av)*1.1)) %>%
  ## tukey letters (higher)
  mutate(label2_pos_am=(max(mean_age_mat)*1.2),
         label2_pos_lam=(max(mean_age_mat)*2.6), ## log
         label2_pos_gc=(max(gc_percent)*1.05),
         label2_pos_oe=(max(mean_nz_cpg_oe_av)*1.175))

## write function to create plots
create_plot <- function(data, y_value, test_value, y_label, 
                        y_mean, label1_pos_col, label2_pos_col) {
  ## tukey letters calculation
  anova_result <- aov(get(test_value) ~ model_group, data = data)
  tukey_result <- glht(anova_result, linfct = mcp(model_group = "Tukey")) ## tukey post hoc
  tukey_result_sum <- summary(tukey_result, test = adjusted("bonferroni")) ## bonf correction
  tukey_letters <- data.frame(letters = cld(tukey_result_sum)$mcletters$Letters) %>% ## get letters
    rownames_to_column(var = "model_group") %>%
    left_join(., data %>%
                select(model_group, {{label2_pos_col}}) %>%
                unique)
  
 ## plot
 plot <- ggplot(data, aes(x = model_group, y = get(y_value), fill = model_group)) +
   geom_violin(scale="width", 
               alpha=0.5, show.legend = F) +
   geom_boxplot(width = 0.07, show.legend = F) +
   geom_sina(alpha=0.5, size=0.25, scale="width") +
  scale_fill_viridis_d(option="D", direction=-1, name="Group") +
  scale_colour_viridis_d(option="D", direction=-1, guide="none") +
  theme_classic() +
  labs(x = "Group", y = y_label) +
  geom_text(data = tukey_letters, 
            aes(x = model_group, y = .data[[label2_pos_col]], label = letters), 
            size = 3.5) +
  geom_text(data = data, 
            aes(x = model_group, y = .data[[label1_pos_col]],
                label = paste("bar(x) ==", round(get(y_mean), 1)), 
                colour = model_group), 
            parse = TRUE, 
            check_overlap = TRUE, 
            size = 2.5)+
  guides(fill = guide_legend(override.aes = list(alpha=1, size=5, shape=21), 
                             direction = "vertical"))
 
 return(plot)
}

## make plots
am_log <- create_plot(phenodata_less, "mean_age_mat", "log_mean_age_mat",
                      "Known age at maturity (yrs)  ", y_mean="mnam",
                      label1_pos_col = "label1_pos_lam", 
                      label2_pos_col = "label2_pos_lam")+scale_y_log10()

gc <- create_plot(phenodata_less, "gc_percent", "gc_percent",
                  "Genome GC (%)", y_mean="mngc",
                  label1_pos_col = "label1_pos_gc", label2_pos_col = "label2_pos_gc")

oe <- create_plot(phenodata_less, "mean_nz_cpg_oe_av", "mean_nz_cpg_oe_av", 
                  "Mean promoter CpG O/E", y_mean="mnoe",
                  label1_pos_col = "label1_pos_oe", label2_pos_col = "label2_pos_oe")

group_legend <- get_legend(oe)
group_legend_gg <- as_ggplot(group_legend)

## combine right hand side of plot
rhs <- ggarrange(ggarrange(am_log+theme(#axis.text.x=element_blank(), 
                                        # axis.ticks=element_blank(),
                                        axis.title.x=element_blank()), 
                           oe+theme(#axis.text.x=element_blank(), 
                                    # axis.ticks=element_blank(),
                                    axis.title.x=element_blank()), 
                           gc+theme(axis.title.x=element_blank()),
                           ncol=1, 
                           # heights=c(0.8, 0.8, 0.9),
                           common.legend=TRUE, 
                           labels=c("B", "C", "D"),
                           legend="none"),
                 sliding_window+
                   theme(legend.text = element_markdown(size=8), 
                         legend.title = element_blank(),
                         legend.position = c(0.24, 0.75)),
                 heights=c(2, 0.75), 
                 labels=c("", "E"), 
                 ncol=1)

rhs

upper <- ggarrange(c1.3+theme(plot.margin=unit(c(-30,-30,-30,-30), "mm")), ## t,r,b,l
                 rhs,
                 labels=c("A", ""),
                 nrow=1,
                 widths=c(0.775, 0.225))

lower <- ggarrange(NULL, 
                   cpg_legend_gg,
                   order_legend_gg,
                   group_legend_gg,
                   NULL,
                   nrow=1,
                   widths = c(0.025, 0.1, 1, 0.1, 0.025))


## COMBINE AND SAVE ----
print("~~~ combine and save ~~~")

gg <- ggarrange(upper, 
                lower,
                nrow=2,
                heights=c(1, 0.19))

## save it
ggsave(paste0("figuresTables/Figure [tree_means_", location, samples, "].pdf"), 
       gg, width=8.25*1.75, height=11.75*0.62*1.75)

ggsave(paste0("figuresTables/Figure [tree_means_", location, samples, "].svg"),
       gg, width=11*1.25, height=8*1.5)

## END SCRIPT

