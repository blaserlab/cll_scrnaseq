theme_set(theme_cowplot(font_size = 11))

# categorical umaps ####-----------------------------------------------------------------------------------
umap_unaligned <- custom_variable_plot(cds_final, var = "pt_timepoint", foreground_alpha = 0.1, overwrite_labels = F)

# reorder the patient column factor levels
colData(cds_aligned)$pt <- factor(colData(cds_aligned)$pt, levels = unique(colData(cds_aligned)$pt))

# make helpful variables for plotting
colData(cds_aligned)$timepoint_pretty <- paste0("Timepoint ", colData(cds_aligned)$timepoint)
colData(cds_aligned)$pt_timepoint <- paste0(colData(cds_aligned)$pt, "_timepoint_", colData(cds_aligned)$timepoint)

umap_sample <- custom_variable_plot(cds_aligned, var = "pt_timepoint", foreground_alpha = 0.1, overwrite_labels = F)

umap_allpts_alltimes <- custom_variable_plot(cds_aligned, var = "pt", foreground_alpha = 0.4)+
  facet_grid(cols = vars(timepoint_pretty), rows = vars(value)) + 
  theme(strip.background = element_blank()) + 
  theme(legend.position = "none") + 
  theme(panel.background = element_rect(color = "grey80"))
umap_allpts_alltimes

umap_partition <- custom_variable_plot(cds_aligned, var = "partition", foreground_alpha = 0.1, overwrite_labels = T, group_label_size = 4)
umap_partition

umap_leiden <- custom_variable_plot(cds_aligned, var = "leiden", foreground_alpha = 0.1, overwrite_labels = T, group_label_size = 4)
umap_leiden

umap_louvain <- custom_variable_plot(cds_aligned, var = "louvain", foreground_alpha = 0.1, overwrite_labels = T, group_label_size = 4)
umap_louvain

umap_seurat_l1 <- custom_variable_plot(cds_aligned, var = "predicted.celltype.l1", foreground_alpha = 0.1, overwrite_labels = T, group_label_size = 4)
umap_seurat_l1

umap_seurat_l2 <- custom_variable_plot(cds_aligned, var = "predicted.celltype.l2", foreground_alpha = 0.1, overwrite_labels = T, group_label_size = 4)
umap_seurat_l2

umap_seurat_coords_l1 <- custom_variable_plot(cds_aligned,var = "predicted.celltype.l1", ref_dim_x = "refUMAP_1", ref_dim_y = "refUMAP_2", foreground_alpha = 0.1, overwrite_labels = T, group_label_size = 4)
umap_seurat_coords_l1

umap_seurat_coords_l2 <- custom_variable_plot(cds_aligned,var = "predicted.celltype.l2", ref_dim_x = "refUMAP_1", ref_dim_y = "refUMAP_2", foreground_alpha = 0.1, overwrite_labels = T, group_label_size = 4)
umap_seurat_coords_l2

# highlight other
umap_seurat_l1_other <- custom_variable_plot(cds_aligned, var = "predicted.celltype.l1", value_to_highlight = "other", foreground_alpha = 0.2, overwrite_labels = T, group_label_size = 4,palette = "#DC0000")


# cluster membership heatmap ####-------------------------------------------------------
col_fun_membership <- 
  colorRamp2(breaks = c(min(scale(t(cluster_timepoint_membership))),0,max(scale(t(cluster_timepoint_membership)))),
             colors = c("#3C5488", "white", "#DC0000"))

cluster_membership_heatmap <-
  grid.grabExpr(draw(
    Heatmap(
      matrix = scale(t(cluster_timepoint_membership)),
      col = col_fun_membership,
      cluster_columns = F,
      column_names_rot = 0,
      column_title = "Timepoint",
      column_title_side = "bottom",
      name = "Scaled Cluster\nMembership"
    )
  ))
grid.draw(cluster_membership_heatmap)
# gene module heatmap ####-------------------------------------------------------------------------------

module_ha_df <- tibble(columns = colnames(agg_mat)) %>%
  mutate(Cluster = str_replace(columns, "_.*","")) %>%
  mutate(Timepoint = str_replace(columns, ".*_","")) %>%
  select(-columns) %>%
  as.data.frame()

col_fun_modules <-
  colorRamp2(breaks = c(min(scale(
    as.matrix(agg_mat)
  )), 0, max(scale(
    as.matrix(agg_mat)
  ))),
  colors = c("#3C5488", "white", "#DC0000"))



module_ha <- HeatmapAnnotation(
  df = module_ha_df,
  # has to be in the same order as the columns it is labeling; there is no other connection
  col = list(
    Cluster = make_named_color_vector(base_vector = unique(module_ha_df$Cluster),
                                      colors = brewer.pal(n = length(unique(module_ha_df$Cluster)),name = "Set2")),
    Timepoint = make_named_color_vector(base_vector = unique(module_ha_df$Timepoint),
                                        colors = c("white","grey80","black"))
  ),
  gp = gpar(col = "grey80"),
  show_annotation_name = T
)

gene_module_heatmap <-
  grid.grabExpr(draw(
    Heatmap(
      matrix = as.matrix(agg_mat),
      col = col_fun_modules,
      heatmap_legend_param = list(
        title = "Module\nScore",
        title_gp = gpar(fontsize = 10),
        labels_gp = gpar(fontsize = 8)
      ),
      show_row_dend = T,
      show_column_names = T,
      show_row_names = T,
      row_names_gp = gpar(fontsize = 10),
      top_annotation = module_ha,
      # border = "grey80",
      # rect_gp = gpar(col = "grey80"),
      column_dend_height = unit(3, "mm"),
      row_title = NULL,
      row_title_side = "left",
      row_title_gp = gpar(fontsize = 10)
    )
  ))

grid.draw(gene_module_heatmap)

# gene module go terms ####----------------------------------------------------------------------------
module_bubble_plots <-
  pmap(
    .l = list(
      simmatrices = subListExtract(module_rrvgo, "simMatrix"),
      reducedterms = subListExtract(module_rrvgo, "reducedTerms"),
      listnames = names(module_rrvgo)
    ),
    .f = function(simmatrices, reducedterms, listnames) {
      p <- rrvgo_scatter(
        simMatrix = simmatrices,
        reducedTerms = reducedterms,
        size = "score",
        labelSize = 4
      ) +
        labs(title = listnames) +
        theme(plot.title = element_text(hjust = 0.5))
      return(p)
    }
  )


module_bubble_plots$`Module 1`
module_bubble_plots$`Module 5`
module_bubble_plots$`Module 7`
module_bubble_plots$`Module 9`

# gene expression #### ----------------------------------------------------------------

plot_cells_alt(cds_aligned, gene_or_genes = "PDE4D") +
  facet_wrap(facets = vars(timepoint))

# violin plots ####---------------------------------------------------------------

custom_violin_plot(cds= cds_aligned[,colData(cds_aligned)$predicted.celltype.l1 == "B"],
                   variable = "timepoint", 
                   # genes_to_plot = c("DGKG", "TCF7","NFKBIA", "EBF1", "PDE4D"),
                   genes_to_plot = c("CD37", "EIF1", "CD79A"),
                   include_jitter = T)

plot_genes_violin(cds_aligned[rowData(cds_aligned)$gene_short_name %in% c("DGKG", "TCF7","NFKBIA", "EBF1"),
                              colData(cds_aligned)$predicted.celltype.l1 == "B"],
                  group_cells_by = "timepoint_pretty")

plot_genes_violin(cds_aligned[rowData(cds_aligned)$gene_short_name %in% c("CD37", "EIF1", "CD79A"),
                              colData(cds_aligned)$predicted.celltype.l1 == "B"],
                  group_cells_by = "timepoint_pretty")

plot_genes_violin(cds_aligned[rowData(cds_aligned)$gene_short_name == "DDX17",
                              colData(cds_aligned)$predicted.celltype.l1 == "B"],
                  group_cells_by = "timepoint_pretty"
                  )



# gene dotplots #### -----------------------------------------------------------
custom_gene_dotplot(
  cds = cds_aligned[, colData(cds_aligned)$predicted.celltype.l1 == "B"],
  # markers = c("DGKG", "TCF7","EBF1", "PDE4D"),
  markers = c("CD37", "EIF1", "CD79A"),
  group_cells_by = "timepoint"
)

# gene expression umaps--------------------------------------------------------------------
dir.create("plots_out/single_gene_umaps")
walk(
  .x = c("BTG1", "RPL39", "RPLP1", "RPLP2", "RPS28", "CD37", "DDX5", "EIF1", "TMSB10", "CD74", "DDX17", "GRB2", "MALAT1", "RALGPS2", "RPS27", "FLNB", "DGKG", "ZDHHC19", "PDE4D", "EBF1", "DENND3", "BHLHE41"),
  .f = function(x) {
    p <- plot_cells_alt(cds = cds_aligned, gene_or_genes = x, alpha = 0.4) +
      facet_wrap(facets = vars(timepoint_pretty)) +
      labs(title = x) +
      theme(plot.title = element_text(hjust = 0.5))
    save_plot(p, file = paste0("plots_out/single_gene_umaps/", x, ".png"), base_width = 7.5, base_height = 4.0)
  }
)

plot_cells_alt(cds = cds_aligned, gene_or_genes = "EBF1", alpha = 0.8)

# gene expression violins

dir.create("plots_out/single_gene_violins")
walk(
  .x = c("BTG1", "RPL39", "RPLP1", "RPLP2", "RPS28", "CD37", "DDX5", "EIF1", "TMSB10", "CD74", "DDX17", "GRB2", "MALAT1", "RALGPS2", "RPS27", "FLNB", "DGKG", "ZDHHC19", "PDE4D", "EBF1", "DENND3", "BHLHE41"),
  .f = function(x) {
    p <- plot_genes_violin(cds = cds_aligned[rowData(cds_aligned)$gene_short_name == x,
                                             colData(cds_aligned)$predicted.celltype.l1 == "B"], 
                           group_cells_by = "timepoint_pretty",)
    save_plot(p, file = paste0("plots_out/single_gene_violins/", x, ".png"), base_width = 5, base_height = 4.0)
  }
)

  
# read count umap
cds_aligned <- monocle3::detect_genes(cds_aligned)
colData(cds_aligned)$num_genes_expressed <- as.numeric(colData(cds_aligned)$num_genes_expressed)
custom_variable_plot(cds = cds_aligned, var = "num_genes_expressed", foreground_alpha = 0.2)

# read count violins

custom_variable_plot(cds = cds_aligned, var = "predicted.celltype.l1", value_to_highlight = c("B","Mono"))

colData(cds_aligned) %>%
  as_tibble() %>%
  filter(predicted.celltype.l1 %in% c("B", "Mono")) %>%
  ggplot(mapping = aes(x = predicted.celltype.l1, y = num_genes_expressed, fill = predicted.celltype.l1)) +
  geom_jitter(color = "black", fill = "transparent", shape = 21, alpha = 0.1) +
  geom_violin(color = "black") +
  scale_fill_viridis_d(alpha = 0.2)

# gene expression violins

priya_interesting_genes_violin <- custom_violin_plot(cds_regression_subset, 
                   genes_to_plot = priya_interesting_genes, 
                   variable = "timepoint_pretty", 
                   rows = 2,
                   include_jitter = T,
                   facet_scale = "free_y",
                   legend_pos = "bottom") +
  theme(axis.text.x = element_blank()) + 
  theme(axis.ticks.x = element_blank())
  
save_plot(priya_interesting_genes_violin, "plots_out/priya_interesting_genes_violin.png", base_width = 7.5, base_height = 4)


