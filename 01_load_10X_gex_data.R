source("/workspace/workspace_pipelines/cll_scrnaseq/00_packages_functions.R")

# load the cell data sets
gex_pipestance_paths <-
  c(
    "~/network/X/Labs/Blaser/single_cell/cll_project/output_cll_rerun_4_8/cll_rerun_4_8_4_1_GEX",
    "~/network/X/Labs/Blaser/single_cell/cll_project/output_cll_rerun_4_8/cll_rerun_4_8_4_2_GEX",
    "~/network/X/Labs/Blaser/single_cell/cll_project/output_cll_rerun_4_8/cll_rerun_4_8_4_3_GEX",
    "~/network/X/Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll5_baseline_5pgex",
    "~/network/X/Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll5_btkclone_5pgex",
    "~/network/X/Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll5_relapse_5pgex",
    "~/network/X/Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll6_baseline_5pgex",
    "~/network/X/Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll6_btkclone_5pgex",
    "~/network/X/Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll6_relapse_5pgex",
    "~/network/X/Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll7_baseline_5pgex",
    "~/network/X/Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll7_btkclone_5pgex",
    "~/network/X/Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll7_relapse_5pgex",
    "~/network/X/Labs/Blaser/single_cell/cll_project/output_cll_rerun_4_8/cll_rerun_4_8_8_1_GEX",
    "~/network/X/Labs/Blaser/single_cell/cll_project/output_cll_rerun_4_8/cll_rerun_4_8_8_2_GEX",
    "~/network/X/Labs/Blaser/single_cell/cll_project/output_cll_rerun_4_8/cll_rerun_4_8_8_3_GEX",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_9_1_GEX",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_9_2_GEX",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_9_3_GEX",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_10_1_GEX",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_10_2_GEX",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_10_3_GEX",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_11_1_GEX",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_11_2_GEX",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_11_3_GEX",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_12_1_GEX",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_12_2_GEX",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_12_3_GEX"
  )

pipestance_names <- c(
  "cds_cll_4_baseline",
  "cds_cll_4_btkclone",
  "cds_cll_4_relapse",
  "cds_cll_5_baseline",
  "cds_cll_5_btkclone",
  "cds_cll_5_relapse",
  "cds_cll_6_baseline",
  "cds_cll_6_btkclone",
  "cds_cll_6_relapse",
  "cds_cll_7_baseline",
  "cds_cll_7_btkclone",
  "cds_cll_7_relapse",
  "cds_cll_8_baseline",
  "cds_cll_8_btkclone",
  "cds_cll_8_relapse",
  "cds_cll_9_baseline",
  "cds_cll_9_btkclone",
  "cds_cll_9_relapse",
  "cds_cll_10_baseline",
  "cds_cll_10_btkclone",
  "cds_cll_10_relapse",
  "cds_cll_11_baseline",
  "cds_cll_11_btkclone",
  "cds_cll_11_relapse",
  "cds_cll_12_baseline",
  "cds_cll_12_btkclone",
  "cds_cll_12_relapse"
)


gex_pipestance_list <- mclapply(
  X = gex_pipestance_paths,
  FUN = load_cellranger_data,
  genome = "GRCh38",
  barcode_filtered = TRUE,
  umi_cutoff = 100,
  mc.preschedule = T,
  mc.cores = 27
)

names(gex_pipestance_list) <- pipestance_names

summarized_sequencing_metrics <- tibble(pipestance_path = gex_pipestance_paths) %>%
  mutate(metrics_summary_path = paste0(pipestance_path,"/outs/metrics_summary.csv")) %>%
  mutate(cds_dim_cells = unname(sapply(X = gex_pipestance_list, FUN = dim)[2,])) %>%
  mutate(cds_name = names(sapply(X = gex_pipestance_list, FUN = dim)[2,])) %>%
  left_join(.,bind_rows(lapply(X = .$metrics_summary_path, FUN = read_csv)), by = c("cds_dim_cells" = "Estimated Number of Cells")) %>%# this is your sanity check.  Joining on cell number derived from the CDS object and the metrics summary
  select(cds_name,cds_dim_cells,`Mean Reads per Cell`,`Median Genes per Cell`,`Fraction Reads in Cells`) %>%
  write_csv("data_out/summarized_sequencing_metrics.csv")
  






# #add cds factor columns
# cds_cll5_baseline<-add_cds_factor_columns(cds = cds_cll5_baseline, columns_to_add = c("pt" = "cll5", "timepoint" = "baseline"))
# cds_cll5_clone<-add_cds_factor_columns(cds = cds_cll5_clone, columns_to_add = c("pt" = "cll5", "timepoint" = "btk_clone"))
# cds_cll5_relapse<-add_cds_factor_columns(cds = cds_cll5_relapse, columns_to_add = c("pt" = "cll5", "timepoint" = "relapse"))
# 
# cds_cll6_baseline<-add_cds_factor_columns(cds = cds_cll6_baseline, columns_to_add = c("pt" = "cll6", "timepoint" = "baseline"))
# cds_cll6_clone<-add_cds_factor_columns(cds = cds_cll6_clone, columns_to_add = c("pt" = "cll6", "timepoint" = "btk_clone"))
# cds_cll6_relapse<-add_cds_factor_columns(cds = cds_cll6_relapse, columns_to_add = c("pt" = "cll6", "timepoint" = "relapse"))
# 
# cds_cll7_baseline<-add_cds_factor_columns(cds = cds_cll7_baseline, columns_to_add = c("pt" = "cll7", "timepoint" = "baseline"))
# cds_cll7_clone<-add_cds_factor_columns(cds = cds_cll7_clone, columns_to_add = c("pt" = "cll7", "timepoint" = "btk_clone"))
# cds_cll7_relapse<-add_cds_factor_columns(cds = cds_cll7_relapse, columns_to_add = c("pt" = "cll7", "timepoint" = "relapse"))
# 
# cds_cll8_baseline<-add_cds_factor_columns(cds = cds_cll8_baseline, columns_to_add = c("pt" = "cll8", "timepoint" = "baseline"))
# cds_cll8_clone<-add_cds_factor_columns(cds = cds_cll8_clone, columns_to_add = c("pt" = "cll8", "timepoint" = "btk_clone"))
# cds_cll8_relapse<-add_cds_factor_columns(cds = cds_cll8_relapse, columns_to_add = c("pt" = "cll8", "timepoint" = "relapse"))
# 
# cds_list<-list(cds_cll5_baseline,cds_cll5_clone,cds_cll5_relapse,cds_cll6_baseline,cds_cll6_clone,cds_cll6_relapse,cds_cll7_baseline,cds_cll7_clone,cds_cll7_relapse,cds_cll8_baseline,cds_cll8_clone,cds_cll8_relapse)
# cds<-combine_cds(cds_list = cds_list, keep_all_genes = TRUE)
# 
# 
# #trim off uninformative genes
# cds_trimmed<-cds[substr(rowData(cds)$gene_short_name,1,2)!="RP",]
# 
# # Normalize and pre-process the data
# cds_trimmed<-preprocess_cds(cds_trimmed, num_dim = 100)
# cds_aligned<-align_cds(cds_trimmed, alignment_group = "pt")
# 
# # Reduce dimensionality and previz cells
# cds_trimmed<-reduce_dimension(cds_trimmed, cores = 39)
# cds_aligned<-reduce_dimension(cds_aligned, cores = 39)
# 
# plot_cells(cds_trimmed, color_cells_by = "pt", label_cell_groups = F)
# plot_cells(cds_aligned, color_cells_by = "pt", label_cell_groups = F)#aligned looks better to start with but maybe is corrected too much.
# 
# #save all original cds data elements
# save.image.pigz("cll_original_cds_elements.RData",n.cores = 39)
# 
# #now remove the unused cds elements to reduce memory footprint and improve performance
# rm(
#   cds_list,
#   cds_cll5_baseline,
#   cds_cll5_clone,
#   cds_cll5_relapse,
#   cds_cll6_baseline,
#   cds_cll6_clone,
#   cds_cll6_relapse,
#   cds_cll7_baseline,
#   cds_cll7_clone,
#   cds_cll7_relapse,
#   cds_cll8_baseline,
#   cds_cll8_clone,
#   cds_cll8_relapse,
#   cds,
#   cds_trimmed
# )
# 
#    
#save.image.pigz("cll_scrnaseq.RData",n.cores = 39)
