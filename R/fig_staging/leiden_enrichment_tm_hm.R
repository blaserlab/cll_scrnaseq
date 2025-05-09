dotplot_markers <- c(
  "FCER2",
  "TCL1A",
  "BTG1",
  "CD52",
  # "NFKBIA",
  # "DUSP1",
  "LTB",
  "FOS",
  "JUN",
  # "S100A6",
  "HLA-A",
  # "HLA-B",
  # "HLA-C",
  # "HLA-E",
  "HLA-DRA",
  # "HLA-DQA1",
  # "HLA-DMA",
  "CNTNAP2",
  "CXXC5",
  "RAC2",
  "TXNIP",
  "CXCR4",
  "CXCR5",
  "LY9",
  "CD83",
  "IGHM"
  # "CD27"
)

subpop_top_markers_mat <-
  scale(t(as.matrix(
    aggregate_gene_expression(
      cds = cds_main[rowData(cds_main)$gene_short_name %in% (cds_main_leiden_comparison_tm %>% pull(gene_short_name)),
                     colData(cds_main)$partition_assignment == "B"],
      cell_group_df =
        colData(cds_main) %>%
        as_tibble(rownames = "cell") %>%
        filter(partition_assignment == "B") %>%
        select(cell, cell_group = leiden_comparison_renamed),
    )
  )))

new_colnames <-
  left_join(
    tibble(id = colnames(subpop_top_markers_mat)),
    rowData(cds_main) %>%
      as_tibble() %>%
      select(id, gene_short_name)
  ) %>%
  pull(gene_short_name)
colnames(subpop_top_markers_mat) <- new_colnames


leiden_enrichment_tm_sh <- SummarizedHeatmap(subpop_top_markers_mat)

blaseRtools::rowData(leiden_enrichment_tm_sh)$cluster <- rownames(rowData(leiden_enrichment_tm_sh))
rowData(leiden_enrichment_tm_sh)
p1 <-
  bb_plot_heatmap_main(leiden_enrichment_tm_sh, tile_color = "transparent") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + scale_fill_gradient2(
    breaks = c(-1.4, 0, 1.4),
    low = "blue3",
    mid = "white",
    high = "red4"
  ) 
p2 <-
  bb_plot_heatmap_colDendro(leiden_enrichment_tm_sh,
                            side = "bottom",
                            linewidth = 0.25
                            )
p3 <-
  bb_plot_heatmap_colHighlight(leiden_enrichment_tm_sh,
                               highlights = dotplot_markers,
                               size = 2.5, box.padding = 1, force = 0.5,
                               segment.size = 0.25)


design <- "
3
1
2
4
"

leiden_enrichment_tm_hm <-
  p1 + p2 + p3 +  guide_area() + plot_layout(design = design,
                                            heights = c(4, 4, 0.5, 1),
                                            guides = "collect") &
  theme(legend.justification = "center", legend.position = "bottom")
leiden_enrichment_tm_hm
