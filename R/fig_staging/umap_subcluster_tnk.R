umap_subcluster_tnk <- 
  bb_var_umap(
    obj = cds_main[,colData(cds_main)$partition_assignment %in% c("T", "NK")],
    var = "leiden_enrichment",
    foreground_alpha = 0.1,
    palette = experimental_group_palette,
    rasterize = TRUE,
    legend_pos = "top"
  ) +
  labs(fill = NULL, color = NULL) + 
  theme(legend.justification = "center") +
  guides(fill = guide_legend(ncol=2, 
                             override.aes = list(size = 2, 
                                                 alpha = 0.5, 
                                                 color = "transparent")))  


umap_subcluster_tnk
