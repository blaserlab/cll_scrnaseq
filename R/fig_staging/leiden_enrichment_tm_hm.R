dotplot_markers <- c(
  "FCER2",
  "TCL1A",
  "BTG1",
  "CD52",
  "NFKBIA",
  "DUSP1",
  "LTB",
  "FOS",
  "JUN",
  "S100A6",
  "HLA-A",
  "HLA-B",
  "HLA-C",
  "HLA-E",
  "HLA-DRA",
  "HLA-DQA1",
  "HLA-DMA",
  "CNTNAP2",
  "CXXC5",
  "RAC2",
  "TXNIP",
  "CXCR4",
  "CXCR5",
  "LY9",
  "CD83",
  "IGHM",
  "CD27"
  
  
  
)

# cds_main_leiden_comparison_tm |>
#   as_tibble() |>
#   filter(cell_group == "other") |> View()

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

col_fun_heatmap_topmarkers <-
  colorRamp2(breaks = c(min(subpop_top_markers_mat),
                        0,
                        max(subpop_top_markers_mat)),
             colors = heatmap_3_colors)

tm_anno <- columnAnnotation(link =  anno_mark(
  at = which(colnames(subpop_top_markers_mat) %in% dotplot_markers),
  labels = colnames(subpop_top_markers_mat)[colnames(subpop_top_markers_mat) %in% dotplot_markers],
  labels_gp = gpar(fontsize = 8),
  padding = 1,
  labels_rot = 30
))


leiden_enrichment_tm_hm <-
  grid.grabExpr(draw(
    Heatmap(
      matrix = subpop_top_markers_mat,
      col = col_fun_heatmap_topmarkers,
      name = "Expression",
      column_dend_height = unit(3, "mm"),
      row_dend_width = unit(3, "mm"),
      row_names_gp = gpar(fontsize = 10),
      show_column_names = FALSE,
      column_dend_side = "bottom",
      top_annotation = tm_anno,
      heatmap_legend_param = list(title_gp = gpar(fontsize = 10), direction = "horizontal")
    ), heatmap_legend_side = "bottom"
  ), wrap = TRUE)

plot_grid(leiden_enrichment_tm_hm)
