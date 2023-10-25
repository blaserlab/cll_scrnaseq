save_figs <- FALSE

# global umap partition_assignment ------------------------
bb_cellmeta(cds_main) |> glimpse()
umap_partition_assignment <- bb_var_umap(
  obj = cds_main,
  var = "partition_assignment",
  foreground_alpha = 0.1,
  rasterize = TRUE,
  overwrite_labels = TRUE
)
umap_partition_assignment + bb_var_umap(cds_main, "seurat_celltype_l1", rasterize = TRUE)

# B cell umap density faceted--------------------------------
umap_density <- 
  bb_var_umap(
    obj = cds_main[,colData(cds_main)$partition_assignment == "B"],
    var = "density",
    sample_equally = TRUE,
    cell_size = 1,
    nbin = 100,
    facet_by = c("patient_type2", "timepoint_merged_1"),
    rows = vars(patient_type2), 
    cols = vars(timepoint_merged_1),
    foreground_alpha = 0.6, 
    rasterize = TRUE
  ) +
  theme(panel.background = element_rect(color = "grey80")) +
  theme(legend.justification = "center") +
  labs(color = "Cell\nDensity")
umap_density


# subcluster_umap-----------------------------------
umap_subcluster <- 
  bb_var_umap(
    obj = cds_main[,colData(cds_main)$partition_assignment == "B"],
    var = "leiden_assignment_binned_2",
    cell_size = 1,
    foreground_alpha = 0.1,
    overwrite_labels = TRUE
  ) 
umap_subcluster
if (save_figs) save_plot(umap_subcluster, filename = fs::path(network_out, "umap_sublcuster.png"), base_width = 4.4, base_height = 4.0)

# volcano plot --------------------

genes_to_highlight_bcell_BTK_timepoints <- c("DGKG", "PDE4D", "TCF7", "NFKBIA", "RGCC", "FMOD", "DENND3")
volcano_BTK_bcells <- 
  pseudobulk_bcell_btk_timepoints[[1]][[2]] %>%
  # filter(str_detect(gene_short_name, "IGH.*|IGK.*|IGL.*", negate = T)) %>%
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58) %>%
  mutate(text_label = ifelse(gene_short_name %in% genes_to_highlight_bcell_BTK_timepoints, gene_short_name, "")) %>% 
  ggplot(
    mapping = aes(
      x = log2FoldChange,
      y = -1*log10(padj),
      colour = threshold,
      fill = threshold,
      label = text_label
    )
  ) +
  geom_point(shape = 21, 
             size = 0.5, 
             alpha = 0.4) +
  geom_text_repel(color = "black", 
                  box.padding = 0.5,
                  point.padding = 0.25,
                  min.segment.length = 0,
                  max.overlaps = 20000,
                  nudge_x = -0.5,
                  size = 3, 
                  segment.size = 0.25,
                  force = 2,
                  seed = 1234,
                  segment.curvature = -0.1,
                  segment.square = TRUE,
                  segment.inflect = TRUE) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey80", "#DC0000")) +
  scale_fill_manual(values = c("transparent", "#DC0000")) +
  labs(caption = "\U21D0 Up in timepoint 1\nUp in timepoint 2 \U21D2", title = NULL)+
  theme(plot.caption.position = "panel") +
  theme(plot.caption = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +  
  coord_cartesian(xlim = c(-1.1*max(abs(range(pseudobulk_bcell_btk_timepoints[[1]][[2]] %>% 
                                                filter(!is.na(padj)) %>% 
                                                pull(log2FoldChange)))), 
                           1.1*max(abs(range(pseudobulk_bcell_btk_timepoints[[1]][[2]] %>% 
                                               filter(!is.na(padj)) %>% 
                                               pull(log2FoldChange))))))
volcano_BTK_bcells
if (save_figs) save_plot(volcano_BTK_bcells, filename = fs::path(network_out, "volcano_plot.png"), base_height = 2.5, base_width = 2.5)

# b cell subpopulation heatmap ------------------------------

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
  # "MALAT1",
  # "PIM3",
  # "MAP3K8",
  # "RELB",
  "CNTNAP2",
  "CXXC5",
  # "ITGB1",
  "RAC2",
  "TXNIP"
  # "KLF6"
)



# cds_main_bcell_subpop_top_markers <- monocle3::top_markers(cds = filter_cds(cds_main, 
#                                                                             cells = bb_cellmeta(cds_main) |> 
#                                                                               filter(partition_assignment == "B"),
#                                                                             genes = bb_rowmeta(cds_main) |> 
#                                                                               filter(gene_short_name %notin% blaseRdata::hg38_remove_genes)), 
#                                                            group_cells_by = "leiden_assignment_binned_2", cores = 24)
# cds_main_bcell_subpop_top_markers |> View()
subpop_top_markers_mat <- scale(t(as.matrix(aggregate_gene_expression(cds = cds_main[rowData(cds_main)$gene_short_name %in% (cds_main_bcell_subpop_top_markers %>% pull(gene_short_name)),
                                         colData(cds_main)$partition_assignment == "B"],
                          cell_group_df = 
                            colData(cds_main) %>%
                            as_tibble(rownames = "cell") %>%
                            filter(partition_assignment == "B") %>%
                            select(cell, cell_group = leiden_assignment_binned_renamed),
                          ))))

new_colnames <- left_join(tibble(id = colnames(subpop_top_markers_mat)),
          rowData(cds_main) %>%
            as_tibble() %>%
            select(id, gene_short_name)) %>%
  pull(gene_short_name)
colnames(subpop_top_markers_mat) <- new_colnames

col_fun_heatmap_topmarkers <- 
  colorRamp2(
    breaks = c(min(subpop_top_markers_mat),
               0,
               max(subpop_top_markers_mat)),
    colors = heatmap_3_colors
  )

tm_anno <- columnAnnotation(link =  anno_mark(
  at = which(colnames(subpop_top_markers_mat) %in% dotplot_markers),
  labels = colnames(subpop_top_markers_mat)[colnames(subpop_top_markers_mat) %in% dotplot_markers],
  labels_gp = gpar(fontsize = 8),
  padding = 1
))


subpop_top_markers_heatmap <- 
  grid.grabExpr(draw(Heatmap(
    matrix = subpop_top_markers_mat,
    col = col_fun_heatmap_topmarkers,
    name = "Expression",
    column_dend_height = unit(3,"mm"), 
    row_dend_width = unit(3,"mm"),
    row_names_gp = gpar(fontsize = 8),
    show_column_names = FALSE,
    column_dend_side = "bottom",
    top_annotation = tm_anno,
    heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface = "bold"))
  )),wrap = TRUE)
plot_grid(subpop_top_markers_heatmap)
if (save_figs) save_plot(plot_grid(subpop_top_markers_heatmap), filename = fs::path(network_out, "subpop_top_markers_heatmap.png"), base_width = 7.5, base_height = 3.5)

# population fold change plot-----------------------------------------------------------
normalized_leiden_counts <- 
  colData(cds_main) %>%
  as_tibble() %>%
  filter(partition_assignment == "B") %>%
  group_by(patient, leiden_assignment_binned_renamed, specimen, timepoint_merged_1, patient_type1) %>%
  summarise(n = n()) %>%
  left_join(colData(cds_main) %>%
              as_tibble() %>%
              group_by(specimen) %>%
              summarise(specimen_total = n())) %>%
  mutate(overall_total = nrow(colData(cds_main))) %>%
  mutate(normalized_count = n*overall_total/specimen_total/2) %>%
  select(leiden_assignment_binned_renamed, specimen, timepoint_merged_1, patient_type1, normalized_count)

cluster_proportion_ratio_plot <- normalized_leiden_counts %>%
  pivot_wider(names_from = leiden_assignment_binned_renamed, values_from = normalized_count, values_fill = 1) %>%
  mutate(btk_to_other_ratio = (`CLL-like`)/(stressed + inflammatory)) %>%
  mutate(log2_ratio = log2(btk_to_other_ratio)) %>%
  ggplot(mapping = aes(x = patient_type1, y = log2_ratio, color = patient_type1, fill = patient_type1)) +
  geom_jitter(shape = jitter_shape, size = jitter_size, stroke = jitter_stroke) +
  facet_wrap(facets = vars(timepoint_merged_1)) +
  scale_fill_manual(values = alpha(colour = experimental_group_palette, alpha = jitter_alpha_fill)) +
  scale_color_manual(values = alpha(colour = experimental_group_palette, alpha = jitter_alpha_color)) +
  theme(strip.background = element_blank()) +
  theme(panel.background = element_rect(color = "grey80")) +
  theme(legend.position = "none") +
  stat_summary(
    fun.data = data_summary_mean_se,
    color = summarybox_color,
    size = summarybox_size,
    width = summarybox_width,
    alpha = summarybox_alpha,
    geom = summarybox_geom, 
    show.legend = FALSE
  ) +
  stat_compare_means(method = "wilcox", label = "p.signif", label.x.npc = "center", label.y = 16, show.legend = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.1))) +
  labs(y = "log<sub>2</sub>(CLL-like:other)", color = "Patient Type", fill = "Patient Type", x = NULL) +
  theme(axis.title.y.left = ggtext::element_markdown()) + 
  theme(axis.text.x.bottom = element_blank()) +
  theme(axis.ticks.x.bottom = element_blank()) +
  theme(legend.position = "bottom", legend.justification = "center")
cluster_proportion_ratio_plot

# bb_cluster_representation2(obj = filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "B", timepoint_merged_1 == "3")), cluster_var = "leiden_assignment_binned_renamed", sample_var = "patient", comparison_var = "patient_type1",comparison_levels = c("responsive", "resistant"), sig_val = "PValue")
# bb_cellmeta(cds_main) |> count(patient_type1)

# module expression heatmaps --------------------------------
#annotation here
col.order <- c("BTK_baseline", "BTK_btk_clone", "BTK_relapse", "MRD_baseline", "MRD_3yrs", "MRD_5yrs")
module_heatmap_anno_df <- data.frame(row.names = colnames(agg_mat_bcells_type_timepoint[,col.order]))
module_heatmap_anno_df$timepoint_merged <- recode(rownames(module_heatmap_anno_df), 
                                                  "BTK_baseline" = "1",
                                                  "BTK_btk_clone" = "2",
                                                  "BTK_relapse" = "3",
                                                  "MRD_baseline" = "1",
                                                  "MRD_3yrs" = "2",
                                                  "MRD_5yrs" = "3")
module_heatmap_anno <- ComplexHeatmap::HeatmapAnnotation(
  df = module_heatmap_anno_df, 
  which = "column",
  col = list(timepoint_merged = c("1" = "white", 
                                  "2" = "grey80", 
                                  "3" = "black")),
  border = TRUE, 
  gp = gpar(col = "black"),
  annotation_label = "Timepoint",
  annotation_name_gp = gpar(fontsize = 8), 
  annotation_legend_param = list(border = "black",
                                 title = "Timepoint",
                                 title_gp = gpar(fontsize = 9, fontface = "bold")))


col_fun_heatmap_bcells <- 
  colorRamp2(
    breaks = c(min(agg_mat_bcells_type_timepoint),
               0,
               max(agg_mat_bcells_type_timepoint)),
    colors = heatmap_3_colors
  )

module_heatmap_bcells <-
  grid.grabExpr(draw(
    Heatmap(matrix = agg_mat_bcells_type_timepoint[, col.order],
            name = "Module\nExpression",
            column_split = c(rep("resistant", times = 3), rep("sensitive", times = 3)),
            col = col_fun_heatmap_bcells,
            row_names_gp = gpar(fontsize = 8),
            column_dend_height = unit(3,"mm"), 
            row_dend_width = unit(3,"mm"),
            bottom_annotation = module_heatmap_anno, 
            show_column_names = F,
            cluster_columns = FALSE, 
            column_title_gp = gpar(fontsize = 9),
            heatmap_legend_param = list(title_gp = gpar(fontsize = 9, fontface = "bold"))
              ), merge_legend = TRUE), wrap = TRUE)

plot_grid(module_heatmap_bcells)
if (save_figs) save_plot(plot_grid(module_heatmap_bcells), filename = fs::path(network_out, "module_heatmap_bcells.png"), base_width = 4, base_height = 4)




tcell_subpop_umap <-
  bb_var_umap(
    cds_main[, colData(cds_main)$partition_assignment == "T"],
    "seurat_l2_leiden_consensus",
    overwrite_labels = T,
    group_label_size = 3,
    foreground_alpha = 0.2
  )
tcell_subpop_umap
if (save_figs) save_plot(tcell_subpop_umap, filename = fs::path(network_out, "tcell_subpop_umap.png"), base_width = 4.4, base_height = 4.0)

# T cell genes---------------------------------------------

{
  colData(cds_main)$sample_cluster <- paste0(colData(cds_main)$sample, "_", colData(cds_main)$seurat_l2_leiden_consensus)
  tcell_geneexp_mat <-
  bb_aggregate(
    filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "T"),
      genes = left_join(tcell_genes, bb_rowmeta(cds_main) |> filter(gene_short_name %in% tcell_genes$gene), by = c("gene" = "gene_short_name")) |> select(feature_id, gene_short_name = gene) |> filter(!is.na(feature_id))
    ),
    cell_group_df = bb_cellmeta(cds_main) |> 
      select(cell_id, sample_cluster)
  ) |> as.matrix() |> t() |> scale() |> t()
tcell_geneexp_mat 
 
rownames(tcell_geneexp_mat) <- left_join(tibble(feature_id = rownames(tcell_geneexp_mat)), bb_rowmeta(cds_main)) |> pull(gene_short_name) 

tcell_genes_anno_df <- tibble(gene_short_name = rownames(tcell_geneexp_mat)) |> 
  left_join(tcell_genes, by = c("gene_short_name" = "gene")) |> 
  as.data.frame()
rownames(tcell_genes_anno_df) <- tcell_genes_anno_df$gene_short_name
tcell_genes_anno_df <- tcell_genes_anno_df["label"]
colnames(tcell_genes_anno_df) <- recode(colnames(tcell_genes_anno_df), "label" = "Gene Set")
tcell_row_anno <- HeatmapAnnotation(df = tcell_genes_anno_df, which = "row", col = list(`Gene Set` = tcell_genes_palette))

tcell_cell_anno_df <-
  tibble(sample_cluster = colnames(tcell_geneexp_mat)) |>
  left_join(
    bb_cellmeta(cds_main) |>
      group_by(
        sample_cluster,
        seurat_l2_leiden_consensus,
        timepoint_merged_1,
        patient_type1
      ) |>
      summarise()
  ) |>
  as.data.frame()
rownames(tcell_cell_anno_df) <- tcell_cell_anno_df$sample_cluster
colnames(tcell_cell_anno_df) <- recode(colnames(tcell_cell_anno_df), 
                                       "seurat_l2_leiden_consensus" = "Cluster",
                                       "timepoint_merged_1" = "Timepoint",
                                       "patient_type1" = "Patient Type  ")
tcell_cell_anno_df$sample_cluster <- NULL
tcell_col_anno <- HeatmapAnnotation(df = tcell_cell_anno_df, which = "column", col = list(Timepoint = timepoint_palette,
                                                                                          `Patient Type  ` = experimental_group_palette,
                                                                                          `Cluster` = seurat_l2_leiden_consensus_colors))

tcell_hm_cols <- circlize::colorRamp2(breaks = c(min(tcell_geneexp_mat), 0, max(tcell_geneexp_mat)),colors = heatmap_3_colors)

}

tcell_genexp_hm_unclustered <-
  grid.grabExpr(draw(
    ComplexHeatmap::Heatmap(
      tcell_geneexp_mat,
      name = "Expr.",
      show_column_names = FALSE,
      right_annotation = tcell_row_anno,
      top_annotation = tcell_col_anno,
      column_title = "Cell Group",
      column_title_side = "bottom",
      row_names_gp = gpar(fontsize = 8),
      cluster_rows = FALSE
    ), heatmap_legend_side = "bottom", annotation_legend_side = "bottom"
  ),
  wrap = TRUE)

if (save_figs) {
  save_plot(
    plot_grid(tcell_genexp_hm_unclustered),
    file = fs::path(network_out, "tcell_geneexp_hm_unclustered.png"),
    base_width = 8.5,
    base_height = 11
  )
  
}


tcell_genexp_hm <-
  grid.grabExpr(draw(
    ComplexHeatmap::Heatmap(
      tcell_geneexp_mat,
      name = "Expr.",
      show_column_names = FALSE,
      right_annotation = tcell_row_anno,
      top_annotation = tcell_col_anno,
      column_title = "Cell Group",
      column_title_side = "bottom",
      row_names_gp = gpar(fontsize = 8),
      cluster_rows = TRUE 
    ), heatmap_legend_side = "bottom", annotation_legend_side = "bottom"
  ),
  wrap = TRUE)

if (save_figs) {
  save_plot(
    plot_grid(tcell_genexp_hm),
    file = fs::path(network_out, "tcell_geneexp_hm.png"),
    base_width = 8.5,
    base_height = 11
  )
  
}

# NK cell genes----------------------------------------

bb_var_umap(cds_main, "partition_assignment")
bb_var_umap(cds_main, "seurat_l2_leiden_consensus")

{
  colData(cds_main)$sample_cluster <- paste0(colData(cds_main)$sample, "_", colData(cds_main)$seurat_l2_leiden_consensus)
  nk_geneexp_mat <-
  bb_aggregate(
    filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "NK"),
      genes = left_join(NK_genes, bb_rowmeta(cds_main) |> filter(gene_short_name %in% NK_genes$gene), by = c("gene" = "gene_short_name")) |> select(feature_id, gene_short_name = gene) |> filter(!is.na(feature_id))
    ),
    cell_group_df = bb_cellmeta(cds_main) |> 
      select(cell_id, sample_cluster)
  ) |> as.matrix() |> t() |> scale() |> t()
nk_geneexp_mat 
 
rownames(nk_geneexp_mat) <- left_join(tibble(feature_id = rownames(nk_geneexp_mat)), bb_rowmeta(cds_main)) |> pull(gene_short_name) 

nk_genes_anno_df <- tibble(gene_short_name = rownames(nk_geneexp_mat)) |> 
  left_join(NK_genes, by = c("gene_short_name" = "gene")) |> 
  as.data.frame()
rownames(nk_genes_anno_df) <- nk_genes_anno_df$gene_short_name
nk_genes_anno_df <- nk_genes_anno_df["label"]
colnames(nk_genes_anno_df) <- recode(colnames(nk_genes_anno_df), "label" = "Gene Set")
nk_row_anno <- HeatmapAnnotation(df = nk_genes_anno_df, which = "row", col = list(`Gene Set` = nk_genes_palette))

nk_cell_anno_df <-
  tibble(sample_cluster = colnames(nk_geneexp_mat)) |>
  left_join(
    bb_cellmeta(cds_main) |>
      group_by(
        sample_cluster,
        seurat_l2_leiden_consensus,
        timepoint_merged_1,
        patient_type1
      ) |>
      summarise()
  ) |>
  as.data.frame()
rownames(nk_cell_anno_df) <- nk_cell_anno_df$sample_cluster
colnames(nk_cell_anno_df) <- recode(colnames(nk_cell_anno_df), 
                                       "seurat_l2_leiden_consensus" = "Cluster",
                                       "timepoint_merged_1" = "Timepoint",
                                       "patient_type1" = "Patient Type  ")
nk_cell_anno_df$sample_cluster <- NULL
nk_col_anno <- HeatmapAnnotation(df = nk_cell_anno_df, which = "column", col = list(Timepoint = timepoint_palette,
                                                                                          `Patient Type  ` = experimental_group_palette,
                                                                                          `Cluster` = seurat_l2_leiden_consensus_colors))

nk_hm_cols <- circlize::colorRamp2(breaks = c(min(nk_geneexp_mat), 0, max(nk_geneexp_mat)),colors = heatmap_3_colors)
}


nk_genexp_hm_unclustered <-
  grid.grabExpr(draw(
    ComplexHeatmap::Heatmap(
      nk_geneexp_mat,
      name = "Expr.",
      show_column_names = FALSE,
      right_annotation = nk_row_anno,
      top_annotation = nk_col_anno,
      column_title = "Cell Group",
      column_title_side = "bottom",
      row_names_gp = gpar(fontsize = 8),
      cluster_rows = FALSE
    ), heatmap_legend_side = "bottom", annotation_legend_side = "bottom"
  ),
  wrap = TRUE)


if (save_figs) {
  save_plot(
    plot_grid(nk_genexp_hm_unclustered),
    file = fs::path(network_out, "nk_geneexp_hm_unclustered.png"),
    base_width = 8.5,
    base_height = 5.5 
  )
  
}


nk_genexp_hm <-
  grid.grabExpr(draw(
    ComplexHeatmap::Heatmap(
      nk_geneexp_mat,
      name = "Expr.",
      show_column_names = FALSE,
      right_annotation = nk_row_anno,
      top_annotation = nk_col_anno,
      column_title = "Cell Group",
      column_title_side = "bottom",
      row_names_gp = gpar(fontsize = 8),
      cluster_rows = TRUE 
    ), heatmap_legend_side = "bottom", annotation_legend_side = "bottom"
  ),
  wrap = TRUE)

if (save_figs) {
  save_plot(
    plot_grid(nk_genexp_hm),
    file = fs::path(network_out, "nk_geneexp_hm.png"),
    base_width = 8.5,
    base_height = 5.5
  )
  
}


# bcell population proportions -----------------------------------
normalized_leiden_counts <- 
  colData(cds_main) %>%
  as_tibble() %>%
  filter(partition_assignment == "B") %>%
  group_by(patient, leiden_assignment_binned_renamed, specimen, timepoint_merged_1, patient_type2) %>%
  summarise(n = n()) %>%
  left_join(colData(cds_main) %>%
              as_tibble() %>%
              group_by(specimen) %>%
              summarise(specimen_total = n())) %>%
  mutate(overall_total = nrow(colData(cds_main))) %>%
  mutate(normalized_count = n*overall_total/specimen_total/2) %>%
  select(leiden_assignment_binned_renamed, specimen, timepoint_merged_1, patient_type2, normalized_count)

cluster_proportion_ratio_plot <- normalized_leiden_counts %>%
  pivot_wider(names_from = leiden_assignment_binned_renamed, values_from = normalized_count, values_fill = 1) %>%
  mutate(btk_to_other_ratio = (`CLL-like`)/(stressed + inflammatory)) %>%
  mutate(log2_ratio = log2(btk_to_other_ratio)) %>%
  ggplot(mapping = aes(x = patient_type2, y = log2_ratio, color = patient_type2, fill = patient_type2)) +
  geom_jitter(shape = jitter_shape, size = jitter_size, stroke = jitter_stroke) +
  facet_wrap(facets = vars(timepoint_merged_1)) +
  scale_fill_manual(values = alpha(colour = experimental_group_palette, alpha = jitter_alpha_fill)) +
  scale_color_manual(values = alpha(colour = experimental_group_palette, alpha = jitter_alpha_color)) +
  theme(strip.background = element_blank()) +
  theme(panel.background = element_rect(color = "grey80")) +
  theme(legend.position = "none") +
  stat_summary(
    fun.data = data_summary_mean_se,
    color = summarybox_color,
    size = summarybox_size,
    width = summarybox_width,
    alpha = summarybox_alpha,
    geom = summarybox_geom, 
    show.legend = FALSE
  ) +
  stat_compare_means(method = "wilcox", label = "p.signif", label.x.npc = "center", label.y = 16, show.legend = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.1))) +
  labs(y = "log<sub>2</sub>(CLL-like:other)", color = "Patient Type", fill = "Patient Type", x = NULL) +
  theme(axis.title.y.left = ggtext::element_markdown()) + 
  theme(axis.text.x.bottom = element_blank()) +
  theme(axis.ticks.x.bottom = element_blank()) +
  theme(legend.position = "bottom", legend.justification = "center")


if (save_figs) {
  save_plot(cluster_proportion_ratio_plot,
    file = fs::path(network_out, "cluster_proportion_ratio_plot.png"),
    base_width = 3,
    base_height = 2.5
  )
  
}

# treg pct plot -----------------------------------------
treg_ratio_tbl <- left_join(
  colData(cds_main) %>%
    as_tibble() %>%
    filter(partition_assignment == "T") %>%
    group_by(specimen, patient, timepoint_merged_1, patient_type2) %>%
    summarise(n_total_t = n()),
  colData(cds_main) %>%
    as_tibble() %>%
    filter(seurat_l2_leiden_consensus == "Treg") %>%
    group_by(specimen)  %>%
    summarise(n_treg = n())
) %>%
  mutate(n_treg = replace_na(n_treg, 1)) %>%
  mutate(treg_pct = n_treg / (n_total_t) * 100)




treg_pct_plot <-
  ggplot(
    treg_ratio_tbl,
    mapping = aes(
      x = patient_type2,
      y = treg_pct,
      color = patient_type2,
      fill = patient_type2
    )
  ) +
  geom_jitter(width = jitter_width,
              size = jitter_size,
              shape = jitter_shape) +
  scale_color_manual(values = experimental_group_palette) +
  scale_fill_manual(values = alpha(alpha = 0.4, colour = experimental_group_palette)) +
  coord_trans(y = "log10", clip = "off") +
  scale_y_continuous(breaks = c(60, 20, 6, 2, 0.6, 0.2) / 2) +
  annotation_logticks(scaled = F,
                      sides = "l",
                      outside = F) +
  theme(legend.position = "none") +
  facet_wrap(facets = vars(timepoint_merged_1)) +
  stat_summary(
    fun.data = data_summary_mean_se,
    color = summarybox_color,
    size = summarybox_size,
    width = summarybox_width,
    alpha = summarybox_alpha,
    geom = summarybox_geom, 
    show.legend = FALSE
  ) +
  stat_compare_means(method = "t.test",
                     label = "p.signif",
                     label.x.npc = 0.5, 
                     show.legend = FALSE) +
  theme(strip.background = element_blank()) + 
  labs(y = "Percent T<sub>reg</sub>", color = "Patient Type", fill = "Patient Type", x = NULL) +
  theme(axis.title.y.left = ggtext::element_markdown()) + 
  theme(axis.text.x.bottom = element_blank()) +
  theme(axis.ticks.x.bottom = element_blank()) +
  theme(legend.position = "bottom", legend.justification = "center")
treg_pct_plot

if (save_figs) {
  save_plot(treg_pct_plot,
    file = fs::path(network_out, "treg_pct_plot.png"),
    base_width = 3.5,
    base_height = 2.5 
  )
  
}

# exhaustion marker gene bubbles ------------------------------
exhaustion_genes <- c("TIGIT", "PDCD1", "LAG3", "CD160", "CD244")
exh_bubdat <- bb_genebubbles(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "T")
), genes = exhaustion_genes, 
cell_grouping = c("seurat_l2_leiden_consensus", "timepoint_merged_1", "patient_type2"), return_value = "data")

exh_genebub <- ggplot(exh_bubdat, aes(x = seurat_l2_leiden_consensus, y = gene_short_name, color = expression, size = proportion)) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() + 
  facet_grid(patient_type2 ~ timepoint_merged_1) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  panel_border() +
  labs(x = NULL, y = NULL) +
  theme(strip.background = element_blank()) +
  theme(axis.text.y = element_text(face = "italic"))
exh_genebub

if (save_figs) {
  save_plot(exh_genebub,
    file = fs::path(network_out, "exh_genebub.png"),
    base_width = 4.5,
    base_height = 3.0 
  )
  
}

# module go term bubbles

mod4_enrichments <- module_enrichment$`Module 4`$res_table |> 
  mutate(cF_numeric = as.numeric(recode(classicFisher, "< 1e-30" = "1e-30"))) |> 
  mutate(neg_log_cF = -log10(cF_numeric)) |> 
  mutate(Term = fct_reorder(Term, neg_log_cF, .desc = FALSE))
mod4_enrichments

mod4_enrichment_plot <- ggplot(mod4_enrichments |> slice_max(neg_log_cF, n = 20), aes(y = Term, x = neg_log_cF, size = Annotated)) +
  geom_point(pch = 21, color = "black", fill = alpha("black", alpha = 0.2)) +
  scale_size_area(limits=c(100, 10000), breaks = (c(300, 1000,3000, 10000))) +
  labs(x = "-log<sub>10</sub>P", y = NULL, size = "Size") +
  theme_minimal_grid(font_size = 10) +
  theme(axis.title.x.bottom = ggtext::element_markdown())

if (save_figs) {
  save_plot(mod4_enrichment_plot,
            file = fs::path(network_out, "mod4_enrichment_plot.png"),
            base_width = 5.5,
            base_height = 3.0 
  )
  
}





