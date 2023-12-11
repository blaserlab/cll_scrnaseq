umap_leiden_enrichment <-
  bb_var_umap(
    filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "B")
    ),
    "leiden_enrichment",
    rasterize = TRUE, 
    legend_pos = "top"
  )
umap_leiden_enrichment
