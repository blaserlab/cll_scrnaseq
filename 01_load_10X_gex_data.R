source("00_packages_functions.R")

# load the cell data sets
cds_cll5_base<-load_cellranger_data(pipestance_path = paste0(home_header,"Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll5_baseline_5pgex"), barcode_filtered = TRUE)
cds_cll5_clone<-load_cellranger_data(pipestance_path = paste0(home_header,"Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll5_btkclone_5pgex"), barcode_filtered = TRUE)
cds_cll5_relapse<-load_cellranger_data(pipestance_path = paste0(home_header,"Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll5_relapse_5pgex"), barcode_filtered = TRUE)

cds_cll6_base<-load_cellranger_data(pipestance_path = paste0(home_header,"Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll6_baseline_5pgex"), barcode_filtered = TRUE)
cds_cll6_clone<-load_cellranger_data(pipestance_path = paste0(home_header,"Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll6_btkclone_5pgex"), barcode_filtered = TRUE)
cds_cll6_relapse<-load_cellranger_data(pipestance_path = paste0(home_header,"Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll6_relapse_5pgex"), barcode_filtered = TRUE)

cds_cll7_base<-load_cellranger_data(pipestance_path = paste0(home_header,"Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll7_baseline_5pgex"), barcode_filtered = TRUE)
cds_cll7_clone<-load_cellranger_data(pipestance_path = paste0(home_header,"Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll7_btkclone_5pgex"), barcode_filtered = TRUE)
cds_cll7_relapse<-load_cellranger_data(pipestance_path = paste0(home_header,"Labs/Blaser/single_cell/cll_project/2019-08-15-cll_btk_scrnaseq/cll7_relapse_5pgex"), barcode_filtered = TRUE)

cds_cll8_base<-load_cellranger_data(pipestance_path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/single_cell/cll_project/output_cll_rerun_4_8/cll_rerun_4_8_8_1_GEX", barcode_filtered = TRUE)
cds_cll8_clone<-load_cellranger_data(pipestance_path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/single_cell/cll_project/output_cll_rerun_4_8/cll_rerun_4_8_8_2_GEX", barcode_filtered = TRUE)
cds_cll8_relapse<-load_cellranger_data(pipestance_path = "/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/single_cell/cll_project/output_cll_rerun_4_8/cll_rerun_4_8_8_3_GEX", barcode_filtered = TRUE)

#add cds factor columns
cds_cll5_base<-add_cds_factor_columns(cds = cds_cll5_base, columns_to_add = c("pt" = "cll5", "timepoint" = "baseline"))
cds_cll5_clone<-add_cds_factor_columns(cds = cds_cll5_clone, columns_to_add = c("pt" = "cll5", "timepoint" = "btk_clone"))
cds_cll5_relapse<-add_cds_factor_columns(cds = cds_cll5_relapse, columns_to_add = c("pt" = "cll5", "timepoint" = "relapse"))

cds_cll6_base<-add_cds_factor_columns(cds = cds_cll6_base, columns_to_add = c("pt" = "cll6", "timepoint" = "baseline"))
cds_cll6_clone<-add_cds_factor_columns(cds = cds_cll6_clone, columns_to_add = c("pt" = "cll6", "timepoint" = "btk_clone"))
cds_cll6_relapse<-add_cds_factor_columns(cds = cds_cll6_relapse, columns_to_add = c("pt" = "cll6", "timepoint" = "relapse"))

cds_cll7_base<-add_cds_factor_columns(cds = cds_cll7_base, columns_to_add = c("pt" = "cll7", "timepoint" = "baseline"))
cds_cll7_clone<-add_cds_factor_columns(cds = cds_cll7_clone, columns_to_add = c("pt" = "cll7", "timepoint" = "btk_clone"))
cds_cll7_relapse<-add_cds_factor_columns(cds = cds_cll7_relapse, columns_to_add = c("pt" = "cll7", "timepoint" = "relapse"))

cds_cll8_base<-add_cds_factor_columns(cds = cds_cll8_base, columns_to_add = c("pt" = "cll8", "timepoint" = "baseline"))
cds_cll8_clone<-add_cds_factor_columns(cds = cds_cll8_clone, columns_to_add = c("pt" = "cll8", "timepoint" = "btk_clone"))
cds_cll8_relapse<-add_cds_factor_columns(cds = cds_cll8_relapse, columns_to_add = c("pt" = "cll8", "timepoint" = "relapse"))

cds_list<-list(cds_cll5_base,cds_cll5_clone,cds_cll5_relapse,cds_cll6_base,cds_cll6_clone,cds_cll6_relapse,cds_cll7_base,cds_cll7_clone,cds_cll7_relapse,cds_cll8_base,cds_cll8_clone,cds_cll8_relapse)
cds<-combine_cds(cds_list = cds_list, keep_all_genes = TRUE)


#trim off uninformative genes
cds_trimmed<-cds[substr(rowData(cds)$gene_short_name,1,2)!="RP",]

# Normalize and pre-process the data
cds_trimmed<-preprocess_cds(cds_trimmed, num_dim = 100)
cds_aligned<-align_cds(cds_trimmed, alignment_group = "pt")

# Reduce dimensionality and previz cells
cds_trimmed<-reduce_dimension(cds_trimmed, cores = 39)
cds_aligned<-reduce_dimension(cds_aligned, cores = 39)

plot_cells(cds_trimmed, color_cells_by = "pt", label_cell_groups = F)
plot_cells(cds_aligned, color_cells_by = "pt", label_cell_groups = F)#aligned looks better to start with but maybe is corrected too much.

#save all original cds data elements
save.image.pigz("cll_original_cds_elements.RData",n.cores = 39)

#now remove the unused cds elements to reduce memory footprint and improve performance
rm(
  cds_list,
  cds_cll5_base,
  cds_cll5_clone,
  cds_cll5_relapse,
  cds_cll6_base,
  cds_cll6_clone,
  cds_cll6_relapse,
  cds_cll7_base,
  cds_cll7_clone,
  cds_cll7_relapse,
  cds_cll8_base,
  cds_cll8_clone,
  cds_cll8_relapse,
  cds,
  cds_trimmed
)

   
save.image.pigz("cll_scrnaseq.RData",n.cores = 39)
