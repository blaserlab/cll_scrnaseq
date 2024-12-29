umap_da_score_tnk <-
  bb_var_umap(
    filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |> filter(partition_assignment %in% c("T", "NK"))
    ),
    "da_score",
    rasterize = TRUE,
    legend_pos = "top",
    hexify = TRUE
  ) +
  scale_fill_gradient2(low = "#DC0000",
                       high = "#3C5488",
                       mid = "white", 
                       breaks = c(-0.9, 0, 0.9),
                       labels = c("IBR", "", "IBS")) +
  scale_color_gradient2(
    low = "#DC0000",
    high = "#3C5488",
    mid = "white",
    guide = "none"
  ) +
  labs(fill = "Differential\nAbundance  ") + 
  theme(legend.justification = "center") 
umap_da_score_tnk
