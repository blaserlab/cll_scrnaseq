umap_leiden_enrichment <-
  bb_var_umap(
    filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "B")
    ),
    "leiden_enrichment1",
    foreground_alpha = 0.05,
    rasterize = TRUE, 
    legend_pos = "top",
    palette = experimental_group_palette
  ) +
  labs(fill = NULL, color = NULL) + 
  theme(legend.justification = "center") +
  guides(fill = guide_legend(ncol=2, 
                             override.aes = list(size = 2, 
                                                 alpha = 0.5, 
                                                 color = "transparent")))  

  
umap_leiden_enrichment
