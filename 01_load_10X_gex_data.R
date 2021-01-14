source("/workspace/workspace_pipelines/cll_scrnaseq/00_packages_functions.R")

analysis_configs <- read_excel("~/network/X/Labs/Blaser/single_cell/cll_project/analysis_configs.xlsx") %>% 
  select(-parental_folder, -`note for priya`) %>% 
  mutate(new_path = str_replace_all(`10X_dir`,"\\\\","/")) %>%
  mutate(new_path = str_replace(new_path, "X:", "~/network/X")) %>%
  select(patient, timepoint, lib_type, directory = new_path) %>%
  mutate(patient = paste0("cll_",patient)) %>%
  mutate(timepoint = as_factor(timepoint)) %>%
  mutate(pipestance_names = paste0("cds_",patient,"_",timepoint))
  
analysis_configs_gex <- analysis_configs %>% filter(lib_type == "gex")

gex_pipestance_list <-
  pmap(
    .l = list(
      patient = analysis_configs_gex$patient,
      timepoint = analysis_configs_gex$timepoint,
      directory = analysis_configs_gex$directory,
      pipestance_names = analysis_configs_gex$pipestance_names
    ),
    .f = function(patient,
                  timepoint,
                  directory,
                  pipestance_names) {
      cds <- load_cellranger_data(directory)
      cds <-
        add_cds_factor_columns(
          cds = cds,
          columns_to_add = c(
            "pt" = patient,
            "timepoint" = timepoint,
            "pipestance" = pipestance_names
          )
        )
      return(cds)
    }
  )

names(gex_pipestance_list) <- analysis_configs_gex$pipestance_names

summarized_sequencing_metrics <-
  tibble(pipestance_path = analysis_configs_gex$directory) %>%
  mutate(metrics_summary_path = paste0(pipestance_path, "/outs/metrics_summary.csv")) %>%
  mutate(cds_dim_cells = unname(sapply(X = gex_pipestance_list, FUN = dim)[2, ])) %>%
  mutate(cds_name = names(sapply(X = gex_pipestance_list, FUN = dim)[2, ])) %>%
  left_join(.,
            bind_rows(lapply(
              X = .$metrics_summary_path, FUN = read_csv
            )),
            by = c("cds_dim_cells" = "Estimated Number of Cells")) %>% # this is your sanity check.  Joining on cell number derived from the CDS object and the metrics summary
  select(
    cds_name,
    cds_dim_cells,
    `Mean Reads per Cell`,
    `Median Genes per Cell`,
    `Fraction Reads in Cells`
  ) %>%
  write_csv("data_out/summarized_sequencing_metrics.csv")





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
#save.image.pigz("cll_scrnaseq_2021.RData",n.cores = 39)
