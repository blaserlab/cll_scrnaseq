umap_da_score <-
  bb_var_umap(
    filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "B")
    ),
    "da_score",
    rasterize = TRUE,
    legend_pos = "top"
  ) +
  scale_fill_gradient2(low = "#DC0000",
                       high = "#3C5488",
                       mid = "white", 
                       breaks = c(-0.9, 0, 0.9)) +
  scale_color_gradient2(
    low = "#DC0000",
    high = "#3C5488",
    mid = "white",
    guide = "none"
  ) +
  labs(fill = "Differential\nAbundance  ") + 
  theme(legend.justification = "center")
