
colData(cds_main)$cd16_da_response1 <- case_match(colData(cds_main)$cd16_da_response, "resistant" ~ "IBR", "sensitive" ~ "IBS", .default = colData(cds_main)$cd16_da_response)

cd16_enrichment_umap <-
  bb_var_umap(filter_cds(
    cds_main,
    cells = bb_cellmeta(cds_main) |> filter(seurat_l2_leiden_consensus == "CD16 Mono")
  ),
  "cd16_da_response1",
  palette = experimental_group_palette,
  plot_title = "CD16+ Monocytes")
colData(cds_main)$cd16_da_response
