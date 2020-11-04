source("/workspace/workspace_pipelines/cll_scrnaseq/00_packages_functions.R")

# load the tcr sequencing metrics summaries
tcr_pipestance_paths <-
  c(
    "~/network/X/Labs/Blaser/single_cell/cll_project/output_cll_rerun_4_8/cll_rerun_4_8_4_1_TCR",
    "~/network/X/Labs/Blaser/single_cell/cll_project/output_cll_rerun_4_8/cll_rerun_4_8_4_2_TCR",
    "~/network/X/Labs/Blaser/single_cell/cll_project/output_cll_rerun_4_8/cll_rerun_4_8_4_3_TCR",
    "~/network/X/Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll5_baseline_tcell",
    "~/network/X/Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll5_btkclone_tcell",
    "~/network/X/Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll5_relapse_tcell",
    "~/network/X/Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll6_baseline_tcell",
    "~/network/X/Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll6_btkclone_tcell",
    "~/network/X/Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll6_relapse_tcell",
    "~/network/X/Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll7_baseline_tcell",
    "~/network/X/Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll7_btkclone_tcell",
    "~/network/X/Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll7_relapse_tcell",
    "~/network/X/Labs/Blaser/single_cell/cll_project/output_cll_rerun_4_8/cll_rerun_4_8_8_1_TCR",
    "~/network/X/Labs/Blaser/single_cell/cll_project/output_cll_rerun_4_8/cll_rerun_4_8_8_2_TCR",
    "~/network/X/Labs/Blaser/single_cell/cll_project/output_cll_rerun_4_8/cll_rerun_4_8_8_3_TCR",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_9_1_TCR",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_9_2_TCR",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_9_3_TCR",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_10_1_TCR",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_10_2_TCR",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_10_3_TCR",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_11_1_TCR",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_11_2_TCR",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_11_3_TCR",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_12_1_TCR",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_12_2_TCR",
    "~/network/X/Labs/Blaser/single_cell/cll_project/cll_9_12/output_cll_9_12_20201030203654/cll_9_12_12_3_TCR"
  )

tcr_pipestance_names <- c(
  "cds_cll_4_baseline_tcr",
  "cds_cll_4_btkclone_tcr",
  "cds_cll_4_relapse_tcr",
  "cds_cll_5_baseline_tcr",
  "cds_cll_5_btkclone_tcr",
  "cds_cll_5_relapse_tcr",
  "cds_cll_6_baseline_tcr",
  "cds_cll_6_btkclone_tcr",
  "cds_cll_6_relapse_tcr",
  "cds_cll_7_baseline_tcr",
  "cds_cll_7_btkclone_tcr",
  "cds_cll_7_relapse_tcr",
  "cds_cll_8_baseline_tcr",
  "cds_cll_8_btkclone_tcr",
  "cds_cll_8_relapse_tcr",
  "cds_cll_9_baseline_tcr",
  "cds_cll_9_btkclone_tcr",
  "cds_cll_9_relapse_tcr",
  "cds_cll_10_baseline_tcr",
  "cds_cll_10_btkclone_tcr",
  "cds_cll_10_relapse_tcr",
  "cds_cll_11_baseline_tcr",
  "cds_cll_11_btkclone_tcr",
  "cds_cll_11_relapse_tcr",
  "cds_cll_12_baseline_tcr",
  "cds_cll_12_btkclone_tcr",
  "cds_cll_12_relapse_tcr"
)


# tcr_pipestance_list <- mclapply(
#   X = tcr_pipestance_paths,
#   FUN = load_cellranger_data,
#   genome = "GRCh38",
#   barcode_filtered = TRUE,
#   umi_cutoff = 100,
#   mc.preschedule = T,
#   mc.cores = 27
# )

# names(tcr_pipestance_list) <- pipestance_names

summarized_sequencing_metrics_tcr <- tibble(pipestance_path = tcr_pipestance_paths) %>%
  mutate(metrics_summary_path = paste0(pipestance_path,"/outs/metrics_summary.csv")) %>%
  bind_cols(bind_rows(lapply(X = .$metrics_summary_path, FUN = read_csv))) %>% 
  select(`Estimated Number of Cells`, `Mean Read Pairs per Cell`) %>%
  write_csv("data_out/summarized_sequencing_metrics_tcr.csv")
summarized_sequencing_metrics_tcr


#save.image.pigz("cll_scrnaseq.RData",n.cores = 39)