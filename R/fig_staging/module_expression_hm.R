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
  annotation_name_gp = gpar(fontsize = 10), 
  annotation_legend_param = list(border = "black",
                                 title = "Timepoint",
                                 title_gp = gpar(fontsize = 10)))


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
            column_split = c(rep("IBR", times = 3), rep("IBS", times = 3)),
            col = col_fun_heatmap_bcells,
            row_names_gp = gpar(fontsize = 9),
            column_dend_height = unit(3,"mm"),
            row_dend_width = unit(3,"mm"),
            bottom_annotation = module_heatmap_anno,
            show_column_names = F,
            cluster_columns = FALSE, 
            column_title_gp = gpar(fontsize = 10),
            heatmap_legend_param = list(title_gp = gpar(fontsize = 10))
              ), merge_legend = TRUE), wrap = TRUE,width = 6)

plot_grid(module_heatmap_bcells)
