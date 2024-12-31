# already filtered out b cells from target and HLA
cellchat_mat <- cellchat_dat |> 
  filter(!str_detect(target, "^B ")) |>
  filter(!str_detect(interaction_name, "^HLA")) |>
  select(rownames, where(is.numeric)) |> 
  bb_tbl_to_matrix() 

cellchat_mat_filtered <- cellchat_mat[rowSums(cellchat_mat)>0.01,]
cellchat_mat_filtered <- pmin(cellchat_mat_filtered, 0.025)


# cellchat_mat_filtered <- log(cellchat_mat_filtered + 0.01)

# color map--------------------
cellchat_mat_col_fun <- circlize::colorRamp2(breaks = c(min(cellchat_mat_filtered), max(cellchat_mat_filtered)), colors = c("transparent", "red3"))

# row annotation-----------------------
cellchat_df <- cellchat_dat |> 
  filter(rownames %in% rownames(cellchat_mat_filtered)) |> 
  select(-where(is.numeric)) |> 
  as.data.frame()
rownames(cellchat_df) <- rownames(cellchat_mat_filtered)

stopifnot(all(cellchat_df$rownames == rownames(cellchat_df)))
stopifnot(all(cellchat_df$rownames == rownames(cellchat_mat_filtered)))

cellchat_df$rownames <- NULL
cellchat_df$annotation <- NULL
cellchat_df$interaction_name <- NULL
cellchat_df$interaction_name_2 <- NULL

source_target_pal <- RColorBrewer::brewer.pal(9, "Set1")
names(source_target_pal) <- unique(cellchat_df$source)

pathway_pal <- RColorBrewer::brewer.pal(10, "Set3")
names(pathway_pal) <- unique(cellchat_df$pathway_name)

cellchat_row_anno <-
  ComplexHeatmap::HeatmapAnnotation(
    which = "row",
    df = cellchat_df,
    col = list(source = source_target_pal,
               target = source_target_pal[!names(source_target_pal) %in% c("B memory", "B naive")],
               pathway_name = pathway_pal),
    annotation_label = list(source = "Source",
                            target = "Target",
                            pathway_name = "Pathway")
  )

# column annotation----------------------
patient_data <- bb_cellmeta(cds_main) |>
  group_by(specimen, timepoint_merged_2, patient_type2) |>
  summarise()

colanno_df <- left_join(tibble(specimen = colnames(cellchat_mat_filtered)),
          patient_data,
          by = join_by("specimen")) |> 
  as.data.frame()
rownames(colanno_df) <- colanno_df$specimen
colanno_df$specimen <- NULL

stopifnot(all(rownames(colanno_df) == colnames(cellchat_mat_filtered)))

cellchat_col_anno <- ComplexHeatmap::HeatmapAnnotation(which = "column", 
                                                       df = colanno_df,
                                                       col = list(patient_type2 = experimental_group_palette,
                                                                  timepoint_merged_2 = experimental_group_palette),
                                                       annotation_label = list(patient_type2 = "Patient Type",
                                                                               timepoint_merged_2 = "Timepoint"))

cellchat_hm <- grid.grabExpr(draw(ComplexHeatmap::Heatmap(cellchat_mat_filtered, 
                        col = cellchat_mat_col_fun, 
                        show_row_names = FALSE, 
                        right_annotation = cellchat_row_anno,
                        show_column_names = FALSE,
                        top_annotation = cellchat_col_anno,
                        name = "Probability")), wrap = TRUE)
