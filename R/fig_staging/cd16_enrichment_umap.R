

cd16_enrichment_umap <-
  bb_var_umap(filter_cds(
    cds_main,
    cells = bb_cellmeta(cds_main) |> filter(seurat_l2_leiden_consensus == "CD16 Mono")
  ),
  "cd16_da_response",
  palette = experimental_group_palette,
  plot_title = "CD16+ Monocytes")
colData(cds_main)$cd16_da_response
