umap_leiden_enrichment_data <- left_join(
  bb_cellmeta(cds_main),
  reducedDims(cds_main)$UMAP |> 
    as_tibble(rownames = "cell_id") |> 
    rename(UMAP_1 = V1, UMAP_2 = V2)) |> 
  filter(partition_assignment == "B") |> 
  select(UMAP_1, UMAP_2, leiden_enrichment1)


umap_leiden_enrichment <- ggplot(
  umap_leiden_enrichment_data,
  aes(
    x = UMAP_1,
    y = UMAP_2,
    color = leiden_enrichment1,
    fill = leiden_enrichment1
  )
) +
  geom_point(pch = 21, stroke = 0.25, size = 1) +
  scale_fill_manual(values = alpha(experimental_group_palette, 0.02)) +
  scale_color_manual(values = alpha(experimental_group_palette, 0.02)) +
  labs(x = "UMAP 1", y = "UMAP 2", fill = NULL) +
  theme(legend.position = "top", legend.justification = "center") +
  guides(
    fill = guide_legend(
      override.aes = list(
        shape = 21,
        size = 2,
        stroke = 0.5,
        fill = alpha(experimental_group_palette[c("IBR", "IBS", "unenriched")], jitter_alpha_fill),
        colour = alpha(experimental_group_palette[c("IBR", "IBS", "unenriched")], jitter_alpha_color)
      ),
      ncol = 3
    ),
    colour = "none"  # Hide duplicate legend from color
  ) + 
  theme(aspect.ratio = 0.9)

umap_leiden_enrichment <- ggrastr::rasterise(umap_leiden_enrichment, dpi = 300)

# Show plot
umap_leiden_enrichment
