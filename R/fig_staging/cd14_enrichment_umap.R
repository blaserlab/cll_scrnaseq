

cd14_enrichment_umap <-
  bb_var_umap(
    filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |> filter(seurat_l2_leiden_consensus == "CD14 Mono")
    ),
    "cd14_louvain_da_response",
    plot_title = "CD14+ Monocytes",
    palette = experimental_group_palette
  )
