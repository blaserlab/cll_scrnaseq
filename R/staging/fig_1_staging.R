bb_cellmeta(cds_main) |> glimpse()

# global umap density faceted--------------------------------
umap_density <- 
  bb_var_umap(
    obj = cds_main[,colData(cds_main)$partition_assignment == "B"],
    var = "density",
    sample_equally = TRUE,
    cell_size = 1,
    nbin = 100,
    facet_by = c("patient_type1", "timepoint_merged_1"),
    rows = vars(patient_type1), 
    cols = vars(timepoint_merged_1),
    foreground_alpha = 0.6
  ) +
  theme(panel.background = element_rect(color = "grey80")) +
  theme(legend.justification = "center") +
  labs(color = "Cell\nDensity")
umap_density

# subcluster_umap-----------------------------------
umap_subcluster <- 
  bb_var_umap(
    obj = cds_main[,colData(cds_main)$partition_assignment == "B"],
    var = "leiden_assignment_binned_renamed",
    cell_size = 1,
    foreground_alpha = 0.1,
    overwrite_labels = TRUE
  ) 
umap_subcluster
# save_plot(umap_subcluster, filename = fs::path(network_out, "umap_sublcuster.png"), base_width = 4.4, base_height = 4.0)

# # volcano plot MRD1 vs BTK cluster----------------------------------------
# genes_to_highlight_MRD1_BTK <- c("FBLN5","PLCB1", "TAL2", "MIR155HG")
# 
# volcano_MRD1_BTK <-
#   pseudobulk_MRD1_BTK[[2]] %>% 
#   filter(str_detect(gene_short_name, "IGH.*|IGK.*|IGL.*", negate = T)) %>%
#   mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.585) %>%
#   mutate(text_label = ifelse(gene_short_name %in% genes_to_highlight_MRD1_BTK, gene_short_name, "")) %>% 
#   ggplot(
#     mapping = aes(
#       x = log2FoldChange,
#       y = -1*log10(padj),
#       colour = threshold,
#       fill = threshold,
#       label = text_label
#     )
#   ) +
#   geom_point(shape = 21, 
#              size = 0.5, 
#              alpha = 0.4) +
#   geom_text_repel(color = "black", 
#                   box.padding = 0.5,
#                   point.padding = 0.25,
#                   min.segment.length = 0,
#                   max.overlaps = 20000,
#                   nudge_x = -0.5,
#                   size = 3, 
#                   segment.size = 0.25,
#                   force = 2,
#                   seed = 1234,
#                   segment.curvature = -0.1,
#                   segment.square = TRUE,
#                   segment.inflect = TRUE) +
#   xlab("log2 fold change") +
#   ylab("-log10 adjusted p") +
#   theme(legend.position = "none") +
#   scale_color_manual(values = c("grey80", "#DC0000")) +
#   scale_fill_manual(values = c("transparent", "#DC0000")) +
#   labs(caption = "\U21D0 Up in BTK Cluster\nUp in MRD1 Cluster \U21D2", title = NULL)+
#   theme(plot.caption.position = "panel") +
#   theme(plot.caption = element_text(hjust = 0.5)) +
#   theme(plot.title = element_text(hjust = 0.5)) +  
#   coord_cartesian(xlim = c(-1.1*max(abs(range(pseudobulk_MRD1_BTK[[2]] %>% 
#                                                 filter(!is.na(padj)) %>% 
#                                                 pull(log2FoldChange)))), 
#                            1.1*max(abs(range(pseudobulk_MRD1_BTK[[2]] %>% 
#                                                filter(!is.na(padj)) %>% 
#                                                pull(log2FoldChange))))))
# volcano_MRD1_BTK
# 
# 
# # volcano plot MRD2 vs BTK cluster----------------------------------------
# genes_to_highlight_MRD2_BTK <- c("RAC2", "PRMT1", "ITGB7", "CHD7", "ZBTB16", "TNKS")
# 
# volcano_MRD2_BTK <- 
#   pseudobulk_MRD2_BTK[[2]] %>%
#   filter(str_detect(gene_short_name, "IGH.*|IGK.*|IGL.*", negate = T)) %>%
#   mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58) %>%
#   mutate(text_label = ifelse(gene_short_name %in% genes_to_highlight_MRD2_BTK, gene_short_name, "")) %>% 
#   ggplot(
#     mapping = aes(
#       x = log2FoldChange,
#       y = -1*log10(padj),
#       colour = threshold,
#       fill = threshold,
#       label = text_label
#     )
#   ) +
#   geom_point(shape = 21, 
#              size = 0.5, 
#              alpha = 0.4) +
#   geom_text_repel(color = "black", 
#                   box.padding = 0.5,
#                   point.padding = 0.25,
#                   min.segment.length = 0,
#                   max.overlaps = 20000,
#                   nudge_x = -0.5,
#                   size = 3, 
#                   segment.size = 0.25,
#                   force = 2,
#                   seed = 1234,
#                   segment.curvature = -0.1,
#                   segment.square = TRUE,
#                   segment.inflect = TRUE) +
#   xlab("log2 fold change") +
#   ylab("-log10 adjusted p") +
#   theme(legend.position = "none") +
#   scale_color_manual(values = c("grey80", "#DC0000")) +
#   scale_fill_manual(values = c("transparent", "#DC0000")) +
#   labs(caption = "\U21D0 Up in BTK Cluster\nUp in MRD2 Cluster \U21D2", title = NULL)+
#   theme(plot.caption.position = "panel") +
#   theme(plot.caption = element_text(hjust = 0.5)) +
#   theme(plot.title = element_text(hjust = 0.5)) +  
#   coord_cartesian(xlim = c(-1.1*max(abs(range(pseudobulk_MRD2_BTK[[2]] %>% 
#                                                 filter(!is.na(padj)) %>% 
#                                                 pull(log2FoldChange)))), 
#                            1.1*max(abs(range(pseudobulk_MRD2_BTK[[2]] %>% 
#                                                filter(!is.na(padj)) %>% 
#                                                pull(log2FoldChange))))))
# volcano_MRD2_BTK


# volcano plot BTK over time:  baseline v btk----------------------
# pseudobulk_bcell_btk_timepoints[[1]][[2]] %>% filter(padj<0.05, abs(log2FoldChange) >= 0.58) %>% View()# baseline v btk
# pseudobulk_bcell_btk_timepoints[[2]][[2]] %>% filter(padj<0.05, abs(log2FoldChange) >= 0.58) %>% View()# baseline v relapse:  nothing
# pseudobulk_bcell_btk_timepoints[[3]][[2]] %>% filter(padj<0.05, abs(log2FoldChange) >= 0.58) %>% View()# btk v relapse: nothing

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

# module heatmap-----------------------------------------------------

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
module_heatmap_anno_df$timepoint_merged
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
            column_split = c(rep("resistant", times = 3), rep("responsive", times = 3)),
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

# BCR signature heatmap, made in same fashion as module heatmap ---------------------

bcr_signature_genes <- c("NFATC2", "SOS2", "PPP3CA", "GSK3B", "HRASLS2", "PPP3CC", "CARD11", "IKBKB", "SOS1", "MALT1", "NFATC4", "PPP3CB", "RAF1", "RAC2", "LILRB3", "MAP2K1", "CD19", "CD22", "CD72", "CD79A", "CD79B", "CD80", "CD81", "GRB2", "BTK", "BCL2", "BCL10", "BIRC3", "BIRC5", "AKT3")
bcr_signature_genes_mat <- bb_aggregate(obj = filter_cds(cds_main, 
                              cells = bb_cellmeta(cds_main) |> 
                                filter(partition_assignment == "B"),
                              genes = bb_rowmeta(cds_main) |> 
                                filter(gene_short_name %in% bcr_signature_genes)),
             cell_group_df = bb_cellmeta(cds_main) |> 
               filter(partition_assignment == "B") |> 
               mutate(type_timepoint1 = paste0(patient_type, "_", timepoint_merged_1)) |> 
               select(cell_id, type_timepoint1)) |> 
  t() |> 
  scale() |> 
  t() |> 
  as.matrix()

{
name_tbl <- tibble(feature_id = rownames(bcr_signature_genes_mat)) |> 
  left_join(bb_rowmeta(cds_main))
sanity_check <- waldo::compare(rownames(bcr_signature_genes_mat), name_tbl$feature_id)
if (length(sanity_check) > 0)
  cli::cli_abort("problems")
rownames(bcr_signature_genes_mat) <- name_tbl$gene_short_name
rm(name_tbl)
}

col_fun_heatmap_bcells1 <- 
  colorRamp2(
    breaks = c(min(bcr_signature_genes_mat),
               0,
               max(bcr_signature_genes_mat)),
    colors = heatmap_3_colors
  )

bcr_signature_genes_hm <-
  grid.grabExpr(draw(
    Heatmap(matrix = as.matrix(bcr_signature_genes_mat),
            name = "Expression",
            column_split = c(rep("BTK", times = 3), rep("MRD", times = 3)),
            col = col_fun_heatmap_bcells1,
            row_names_gp = gpar(fontsize = 10),
            column_dend_height = unit(3,"mm"), 
            row_dend_width = unit(3,"mm"),
            bottom_annotation = module_heatmap_anno, 
            show_column_names = F,
            cluster_columns = FALSE
            ), merge_legend = TRUE), wrap = TRUE)
plot_grid(bcr_signature_genes_hm)

# save_plot(plot_grid(bcr_signature_genes_hm), filename = fs::path(network_out, "bcr_signature_genes_hm.png"), base_width = 5, base_height = 7)


# subpop gene dotplot

dotplot_markers <- c(
  "FCER2",
  "TCL1A",
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
  "MALAT1",
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
# 
# subpop_gene_dotplot <- bb_gene_dotplot(cds = cds_main[,colData(cds_main)$leiden_assignment_binned %in% c("BTK", "MRD1","MRD2")], 
#                 group_cells_by = "leiden_assignment_binned",
#                 markers = dotplot_markers) +
#   labs(x = NULL, y = NULL)
# subpop_gene_dotplot


# subpopulation heatmap---------------------------------------------------------------

# show top specific genes for subclusters -------------------------------------
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

# heatmap of patient characteristics

pt_char <- bb_cellmeta(cds_main) |>
  group_by(
    specimen,
    patient,
    patient_type,
    timepoint_merged_1,
    gender,
    age_at_ibr_start,
    IGHV_status,
    complex_karyotype,
    del17p,
    del11q,
    del13p,
    tri12,
    btk_clone_vaf_pct,
    relapse_vaf_pct,
    partition_assignment
  ) |>
  summarise(n = n()) |>
  pivot_wider(names_from = partition_assignment,
              values_from = n,
              values_fill = 0) |>
  ungroup()

pt_char_mat <- pt_char |>
  select(specimen, where(is.integer)) |>
  bb_tbl_to_matrix() |>
  apply(1, \(x) x / sum(x) * 100)

{
pt_char_anno_df <- pt_char |> 
  select(-where(is.integer)) |> 
  mutate(age_at_ibr_start = as.numeric(age_at_ibr_start)) |>
  mutate(btk_clone_vaf_pct = as.numeric(btk_clone_vaf_pct)) |> 
  mutate(relapse_vaf_pct = as.numeric(relapse_vaf_pct)) |> 
  as.data.frame()
rownames(pt_char_anno_df) <- pt_char_anno_df$specimen
pt_char_anno_df$specimen <- NULL
}
waldo::compare(rownames(pt_char_anno_df), colnames(pt_char_mat))

col_fun_pt_age <-
  circlize::colorRamp2(breaks = c(
    min(as.numeric(pt_char$age_at_ibr_start)),
    max(as.numeric(pt_char$age_at_ibr_start))
  ),
  colors = c("transparent", "midnightblue"))

col_fun_btk_clone_vaf <-
  circlize::colorRamp2(breaks = c(
    min(as.numeric(pt_char$btk_clone_vaf_pct), na.rm = TRUE),
    max(as.numeric(pt_char$btk_clone_vaf_pct), na.rm = TRUE)
  ),
  colors = c("transparent", "forestgreen"))

col_fun_relapse_vaf <-
  circlize::colorRamp2(breaks = c(
    min(as.numeric(pt_char$relapse_vaf_pct), na.rm = TRUE),
    max(as.numeric(pt_char$relapse_vaf_pct), na.rm = TRUE)
  ),
  colors = c("transparent", "darkmagenta"))

pt_char_anno <- ComplexHeatmap::HeatmapAnnotation(df = pt_char_anno_df, 
                                                  which = "column",
                                                  col = list(IGHV_status = c("M" = "black", "U" = "white"),
                                                             complex_karyotype = c("yes" = "black", "no" = "white", "unknown" = "grey80"),
                                                             timepoint_merged_1 = timepoint_palette,
                                                             del17p = c(`FALSE` = "white", `TRUE` = "black"),
                                                             del11q = c(`FALSE` = "white", `TRUE` = "black"),
                                                             gender = c("F" = "pink", "M" = "blue"),
                                                             del13p = c(`FALSE` = "white", `TRUE` = "black"),
                                                             tri12 = c(`FALSE` = "white", `TRUE` = "black"),
                                                             patient_type = experimental_group_palette, 
                                                             patient = pt_colors,
                                                             age_at_ibr_start = col_fun_pt_age,
                                                             btk_clone_vaf_pct = col_fun_btk_clone_vaf,
                                                             relapse_vaf_pct = col_fun_relapse_vaf
                                                             
                                                             ))

col_fun_pt_char <- circlize::colorRamp2(breaks = c(0,100), colors = c("transparent", "red3"))



pt_char_hm <- grid.grabExpr(draw(ComplexHeatmap::Heatmap(
  pt_char_mat,
  col = col_fun_pt_char,
  name = "Percent\nof Sample",
  column_dend_side = "bottom",
  show_column_names = FALSE,
  top_annotation = pt_char_anno
  
)), wrap = TRUE, width = unit(5, "in"), height = unit(5, "in")) 

# save_plot(filename = fs::path(network_out, "pt_char_hm.png"), plot_grid(pt_char_hm), base_width = 10, base_height = 7.5)

#  barchart by patient type

bb_cellmeta(cds_main) |> 
  count(patient_type1, partition_assignment, timepoint_merged) |> 
  mutate(timepoint_merged = recode(timepoint_merged, "baseline" = 1L, "3yrs|btk_clone" =  2L, "5yrs|relapse" = 3L)) |> 
  pivot_wider(names_from = "partition_assignment", values_from = "n", values_fill = 0) |> 
  write_tsv(fs::path(network_tables, "cell_composition_absolute.csv"))

cell_composition_barchart <- bb_cellmeta(cds_main) |>
  group_by(patient_type1) |> 
  slice_sample(n = 25101) |> 
  ungroup() |> 
  count(patient_type1, timepoint_merged, partition_assignment) |> 
  pivot_wider(names_from = "partition_assignment", values_from = n, values_fill = 0) |> 
  mutate(timepoint_merged = recode(timepoint_merged, "baseline" = 1L, "3yrs|btk_clone" =  2L, "5yrs|relapse" = 3L)) |> 
  pivot_longer(-c(patient_type1, timepoint_merged), names_to = "partition_assignment", values_to = "count") |> 
  group_by(partition_assignment, timepoint_merged) |> 
  mutate(pct = count/sum(count)*100) |> 
  ggplot(aes(x = partition_assignment, y = pct, fill = patient_type1)) +
  geom_col() +
  facet_wrap(~timepoint_merged, scales = "free_x") +
  labs(fill = NULL, x = NULL)
cell_composition_barchart

# save_plot(cell_composition_barchart, filename = fs::path(network_out, "cell_comp_barchart.png"), base_width = 10, base_height = 7.5)

# module expression heatmaps-----------------------


mod_expr_hms <-
  map(.x = list(
    unique(colData(cds_main)$partition_assignment),
    "B",
    "T",
    "NK",
    "Mono"),
      .f = \(
        x,
        dat = cds_main,
        egp = experimental_group_palette,
        tp = timepoint_palette,
        pc = partition_colors
      ) {
        dat <-
          filter_cds(dat, cells = bb_cellmeta(dat) |> filter(partition_assignment %in% x))
        module_mat <- bb_aggregate(
          dat,
          gene_group_df = bb_rowmeta(dat) |> select(feature_id, module),
          cell_group_df = bb_cellmeta(dat) |> select(
            cell_id,
            patient_type,
            timepoint_merged_1,
            partition_assignment
          ) |> mutate(
            var = paste(
              patient_type,
              timepoint_merged_1,
              partition_assignment,
              sep = "_"
            )
          ) |>
            select(cell_id, var)
        )
        
        {
          module_anno_df <- tibble(var = colnames(module_mat)) |>
            separate(
              col = "var",
              sep = "_",
              into = c(
                "patient_type",
                "timepoint_merged_1",
                "partition_assignment"
              )
            ) |>
            as.data.frame()
          rownames(module_anno_df) <- module_anno_df$var
          module_anno_df$var <- NULL
        }
        
        module_anno <-
          ComplexHeatmap::HeatmapAnnotation(
            df = module_anno_df,
            col = list(
              patient_type = egp,
              timepoint_merged_1 = tp,
              partition_assignment = pc
            )
          )
        
        grid.grabExpr(draw(
          ComplexHeatmap::Heatmap(
            as.matrix(module_mat),
            name = "Module\nExpression",
            top_annotation = module_anno,
            show_column_names = FALSE,
            column_dend_side = "bottom",
            row_title = "Module", 
            cluster_columns = ifelse(length(x) == 1, FALSE, TRUE)
          )
        ), wrap = TRUE)
        
        
      }) |>
  set_names(c("all", "B", "T", "NK", "Mono"))
  # set_names("Mono")

# plot_grid(mod_expr_hms$B)
# plot_grid(mod_expr_hms$all)
# walk2(.x = mod_expr_hms,
#       .y = names(mod_expr_hms),
#      .f = \(x, y, out = network_out) save_plot(plot = x, filename = fs::path(out, paste0(y, "_module_expr.png")), base_height = 7.5, base_width = 10))
# 
# names(mod_expr_hms[1])

