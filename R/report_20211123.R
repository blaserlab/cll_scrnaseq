#' ---
#' title: "BTK Single Cell Mutation Analysis"
#' author: "Brad Blaser"
#' date: "Nov 23 2021"
#' output: pdf_document
#' ---
#'
#+ include=FALSE
knitr::opts_chunk$set(warning = FALSE, message = FALSE, echo=FALSE)
source("R/dependencies.R")
source("R/configs.R")

#' Following up on our conversation yesterday, I reviewed the data in more detail and found a few things worth mentioning. 
#' 
#' The figures I showed used a mutant vaf threshold of 25% for calling cells mutant vs WT.  Looking at the counts, I found that most of the non-B cells with mutant BTK reads were called mutant based on a single mutant read.  I think it is a reasonable threshold to require 2 or more mutant reads before calling a cell mutant.
#' 
#' When I apply this threshold (mutant VAF > 25% and > 1 mutant BTK read) this is what we get.  

#+ echo = TRUE, fig.width = 6.5, fig.height = 3.0, dev = "png", fig.align = "center"
bb_var_umap(cds_main[,colData(cds_main)$btk_type != "no_type"], 
            "density", 
            facet_by = "btk_type")
#'
#'

# There aren't enough non-b cells with mutant btk to worry about
mono_typing <- 
  bb_cellmeta(cds_main) %>%
  group_by(patient, timepoint, partition_assignment, btk_type) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = btk_type, values_from = n, values_fill = 0) %>%
  filter(WT > 0) %>%
  mutate(pct_mutant= mutant/(WT + mutant) * 100) %>%
  select(-c(patient)) %>%
  filter(partition_assignment == "Mono") %>%
  ungroup() %>%
  select(-partition_assignment)
T_typing <- 
  bb_cellmeta(cds_main) %>%
  group_by(patient, timepoint, partition_assignment, btk_type) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = btk_type, values_from = n, values_fill = 0) %>%
  filter(WT > 0) %>%
  mutate(pct_mutant= mutant/(WT + mutant) * 100) %>%
  select(-c(patient)) %>%
  filter(partition_assignment == "T") %>%
  ungroup() %>%
  select(-partition_assignment)

#+ echo = FALSE
pander(mono_typing, caption = "Monocyte Typing Results By Sample")
#'

#+ echo = FALSE
pander(T_typing, caption = "T Cell Typing Results By Sample")
#'

#' For reference, the overall proportion of mutant BTK in the B cells is 2.4%.  My conclusion is that these are simply not enough cells to try to study, even if we set aside the concerns about false positives from misassigned cell barcodes and inaccurate genotyping and believe the data.  These cells are so rare you aren't going to be able to detect them in sorted populations without getting contamination from B cells.  So I think we should cancel plans to pursue this finding with more experiments.  Happy to hear other thoughts.

#' Despite this I think we can still look at BTK typing in B cells.
#' 
#' The approach I take here is to use the typing data to identify which unsupervised clusters of cells are enriched with mutant cells.  This is because the typing data are sparse and you get more data when you look at more cells.
#' 
#'  We can calculate the expected number of mutant cells in each leiden subcluster based on the overall rate of mutant cells and then use the binomial test to determine if the observed mutant cell count is significantly greater than expected:
#'  
#'  

#+ echo = FALSE
pander(btk_leiden_enrichment, caption = "BTK Mutant Cells By Leiden Cluster")

#' The subclusters where mutant_enriched equals true are in the panhandle James commented on.
#' 

#+ echo = TRUE, fig.width = 5.5, fig.height = 3.5, dev = "png", fig.align = "center"
bb_var_umap(cds_main[,!is.na(colData(cds_main)$leiden_mutant_enriched)], 
            "leiden_mutant_enriched")
#'
#'

#' We can then calculate the top differentially expressed genes between the mutant enriched and non-enriched cells using each patient as a biological replicate.   

#+ echo = TRUE
write_csv(leiden_mutant_enriched_pseudobulk_tbl %>% 
            arrange(padj), 
          file = 
            str_glue("{network_tables}/leiden_mutant_enriched_pseudobulk.csv"))
#'

#' I think it would be good for Priya to look at this gene list and run IPA or your favorite pathway analyzer.  I can also run but best for Priya to work on this first since she knows the biology better.  We can pick out single genes or groups of genes to look at and see if they are interesting.  She should be able to find the table in the X drive folder.
