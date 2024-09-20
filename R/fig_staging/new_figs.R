
bb_cellmeta(cds_main) |> count(seurat_celltype_l2) |> View()

# leiden umap -----------------------------------------------
umap_leiden_l1 <-
  bb_var_umap(
    cds_main,
    "leiden_l1_assignment",
    value_to_highlight = c("B", "CD4 T", "CD8 T", "NK"),
    facet_by = c("timepoint_merged_2", "patient_type2"),
    cols = vars(timepoint_merged_2),
    rows  = vars(patient_type2),
    foreground_alpha = 0.2
  ) + panel_border()


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


# bb_cluster_representation2(obj = filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "B", timepoint_merged_1 == "3")), cluster_var = "leiden_assignment_binned_renamed", sample_var = "patient", comparison_var = "patient_type1",comparison_levels = c("responsive", "resistant"), sig_val = "PValue")
# bb_cellmeta(cds_main) |> count(patient_type1)

# module expression heatmaps --------------------------------
#annotation here

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
plot_grid(tcell_genexp_hm)

 # NK cell genes----------------------------------------


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

  
plot_grid(nk_genexp_hm)






clinical_flow_data |> 
  mutate(population = recode(population, "CD 4" = "CD4")) |> 
  ggplot(aes(x = patient_type, y = abs_mm3, fill = patient_type, color = patient_type)) +
  geom_jitter(width = jitter_width,
              size = jitter_size,
              shape = jitter_shape) +
  scale_color_manual(values = experimental_group_palette) +
  scale_fill_manual(values = alpha(alpha = 0.4, colour = experimental_group_palette)) +
  stat_summary(
    fun.data = data_summary_mean_se,
    color = summarybox_color,
    size = summarybox_size,
    width = summarybox_width,
    alpha = summarybox_alpha,
    geom = summarybox_geom, 
    show.legend = FALSE
  )+
  stat_compare_means(method = "wilcox",
                     label = "p.signif",
                     label.x.npc = 0.5,
                     label.y.npc = 0.8,
                     show.legend = FALSE) +
  theme(strip.background = element_blank()) +
  facet_grid(population~timepoint_merged_2, scales = "free") + 
  panel_border()
 








shannon_diversity <- function(pop) {
  n <- pop
  N <- sum(pop)
  p <- n/N
  -sum(p*log(p))
}

# tcr diversity recalculated by subcluster ------------------------
bb_cellmeta(cds_main) |>
  group_by(sample, patient_type2, timepoint_merged_2) |>
  summarise() |>
  left_join(
    bb_cellmeta(cds_main) |>
      filter(partition_assignment == "T") |>
      filter(!is.na(tcr_clone_copies)) |>
      group_by(sample, seurat_l2_leiden_consensus) |>
      summarize(tcr_shannon = shannon_diversity(tcr_clone_copies)),
    by = join_by(sample)
  ) |> 
  ggplot(aes(x = patient_type2, fill = patient_type2, color = patient_type2, y = tcr_shannon)) + 
  geom_jitter() + 
  # facet_wrap(~timepoint_merged_2) + 
  facet_grid(seurat_l2_leiden_consensus~timepoint_merged_2, scales = "free") +
  stat_summary(
    fun.data = data_summary_mean_se,
    color = summarybox_color,
    size = summarybox_size,
    width = summarybox_width,
    alpha = summarybox_alpha,
    geom = summarybox_geom, 
    show.legend = FALSE
  ) +
  stat_compare_means(method = "wilcox",
                     label = "p.format",
                     label.x.npc = 0.5,
                     label.y.npc = 0.8,
                     show.legend = FALSE) + 
  panel_border()


clinical_flow_data |> 
  ggplot(aes(x = patient_type, y = abs_mm3, color = patient_type, fill = patient_type)) +
  geom_jitter() + 
  # facet_wrap(~timepoint_merged_2) + 
  facet_grid(population~timepoint_merged_2, scales = "free") +
  stat_summary(
    fun.data = data_summary_mean_se,
    color = summarybox_color,
    size = summarybox_size,
    width = summarybox_width,
    alpha = summarybox_alpha,
    geom = summarybox_geom, 
    show.legend = FALSE
  ) +
  stat_compare_means(method = "wilcox",
                     label = "p.format",
                     label.x.npc = 0.5,
                     label.y.npc = 0.8,
                     show.legend = FALSE) + 
  panel_border()


mcpas_tcr

# monocytes ---------------------------------
cds_cd14 <- filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(seurat_l2_leiden_consensus == "CD14 Mono"))

bb_var_umap(cds_cd14,"da_score")
bb_var_umap(cds_cd14,"louvain")
louvain_da_response <- bb_cellmeta(cds_cd14) |>
  group_by(louvain) |> summarise(mean_da_score = mean(da_score)) |> 
  mutate(louvain_da_response = case_when(mean_da_score > 0.2 ~ "sensitive",
                                         mean_da_score < -0.2 ~ "resistant",
                                         TRUE ~ "unenriched")) |> 
  select(louvain, louvain_da_response) |> 
  deframe()

colData(cds_cd14)$louvain_da_response <- recode(colData(cds_cd14)$louvain, !!!louvain_da_response)

bb_var_umap(cds_cd14, "louvain_da_response")


bb_cellmeta(cds_main) |> count(partition_assignment)
bb_cellmeta(cds_main) |> glimpse()
bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "Mono")), "density", facet_by = c("patient_type2", "timepoint_merged_2"))
bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "Mono")), "density", facet_by = c("specimen"))
bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "Mono")), "louvain", overwrite_labels = TRUE)
bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "Mono")), "seurat_l2_leiden_consensus", overwrite_labels = TRUE)
bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "Mono")), "seurat_celltype_l2", overwrite_labels = FALSE, facet_by = "value")

bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "Mono")), "da_score")
bb_gene_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "Mono")), "TLR4")




bb_var_umap(cds_main, "da_score")
bb_var_umap(cds_main, "louvain")

bb_cellmeta(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "Mono"))) |>
  group_by(louvain) |> summarise(mean_da_score = mean(da_score))
  
cds_main_top_markers |> filter(cell_group %in% c("louvain 4",
                                                 "louvain 12",
                                                 "louvain 16",
                                                 "louvain 18",
                                                 "louvain 57",
                                                 "louvain 58",
                                                 "louvain 60",
                                                 "louvain 61",
                                                 "louvain 62",
                                                 "louvain 64")) |> View()

# library(topGO)
# cd14_sens_genes <- cd14_pseudobulk_res$Result |> filter(padj<0.05, log2FoldChange>0) |> pull(gene_short_name)
# cd14_resistant_genes <- cd14_pseudobulk_res$Result |> filter(padj<0.05, log2FoldChange<0) |> pull(gene_short_name)
# 
# cd14_goeenrich <- map(.x = list(cd14_sens_genes,
#                                 cd14_resistant_genes),
#                       .f = \(x, ref = cds_main) {
#                         bb_goenrichment(query = x, 
#                                         reference = bb_rowmeta(ref), 
#                                         go_db = "org.Hs.eg.db")
#                       }) |> 
#   set_names(c("sensitive", "resistant"))
# 
# cd14_gosummary <- map2(.x = c(0.8, 0.9, 0.8, 0.9),
#                             .y = list(cd14_goeenrich$sensitive,
#                                       cd14_goeenrich$sensitive,
#                                       cd14_goeenrich$resistant,
#                                       cd14_goeenrich$resistant),
#                            .f = \(x, dat = cd14_sens_goeenrich) {
#                              bb_gosummary(x = dat, reduce_threshold = x, go_db = "org.Hs.eg.db")
#                            }) |> 
#   set_names(c("sensitive_0.8", "sensitive_0.9", "resistant_0.8", "resistant_0.9"))
# 
# bb_goscatter(simMatrix = cd14_gosummary$sensitive_0.8$simMatrix,
#              reducedTerms = cd14_gosummary$sensitive_0.8$reducedTerms)
# bb_goscatter(simMatrix = cd14_gosummary$resistant_0.8$simMatrix,
#              reducedTerms = cd14_gosummary$resistant_0.8$reducedTerms)
# 
# bb_goscatter(simMatrix = cd14_gosummary$sensitive_0.9$simMatrix,
#              reducedTerms = cd14_gosummary$sensitive_0.9$reducedTerms)
# bb_goscatter(simMatrix = cd14_gosummary$resistant_0.9$simMatrix,
#              reducedTerms = cd14_gosummary$resistant_0.9$reducedTerms)

cd14_goeenrich$sensitive$res_table |> View()

bb_cellmeta(cds_main) |> glimpse()
bb_var_umap(cds_monos, "cd16_da_response")
