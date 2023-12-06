nejm_bcr_sig <- readxl::read_excel("/home/OSUMC.EDU/blas02/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/queries/Bcell_BCR_signature.xlsx")

nejm_bcr_sig_genes <- nejm_bcr_sig$gene_short_name


rowData(cds_main)$nejm_bcr_sig <- ifelse(rowData(cds_main)$gene_short_name %in% nejm_bcr_sig_genes, "yes_bcr", "no_bcr")
cds_main_bcr_only <- cds_main[rowData(cds_main)$nejm_bcr_sig == "yes_bcr",
                              colData(cds_main)$partition_assignment == "B"]

nejm_bcr_umap <- bb_gene_umap(cds_main_bcr_only, gene_or_genes = bb_rowmeta(cds_main_bcr_only) |> select(id, nejm_bcr_sig))

nejm_bcr_umap_faceted <- nejm_bcr_umap + facet_grid(cols = vars(timepoint_merged), rows = vars(patient))

save_plot(nejm_bcr_umap, filename = file.path(network_out, "nejm_bcr_umap.png"), base_height = 4, base_width = 5)

save_plot(nejm_bcr_umap_faceted, filename = file.path(network_out, "nejm_bcr_umap_faceted.png"), base_height = 10, base_width = 6.5)


bb_cellmeta(cds_main_bcr_only)

colData(cds_main_bcr_only)$type_timepoint <- paste0(colData(cds_main_bcr_only)$patient_type, "_", colData(cds_main_bcr_only)$timepoint)
colData(cds_main_bcr_only)$type_timepoint <- factor(colData(cds_main_bcr_only)$type_timepoint, 
                                                    levels = c("BTK_baseline", "BTK_btk_clone", "BTK_relapse", "MRD_baseline", "MRD_3yrs", "MRD_5yrs"))

nejm_mtx <- aggregate_gene_expression(cds_main_bcr_only, cell_group_df = bb_cellmeta(cds_main_bcr_only) |> select(cell_id, type_timepoint), scale_agg_values = F) |> 
  as.matrix() |> 
  t() |> 
  scale() |> 
  t()
rownames(nejm_mtx) <- left_join(tibble(feature_id = rownames(nejm_mtx)), bb_rowmeta(cds_main_bcr_only)) |> pull(gene_short_name)


Heatmap(nejm_mtx, cluster_columns = F)
