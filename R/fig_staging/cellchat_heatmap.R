cellchat_mat <- cellchat_dat |> 
  filter(!str_detect(target, "^B ")) |>
  filter(!str_detect(interaction_name, "^HLA")) |> 
  select(rownames, where(is.numeric)) |> 
  bb_tbl_to_matrix() 


cellchat_mat_filtered <- cellchat_mat[rowSums(cellchat_mat)>0.01,]
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

cellchat_df$rownames <- NULL
cellchat_df$annotation <- NULL
cellchat_df$interaction_name <- NULL
cellchat_df$interaction_name_2 <- NULL

cellchat_row_anno <- ComplexHeatmap::rowAnnotation(df = cellchat_df)


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

cellchat_col_anno <- ComplexHeatmap::columnAnnotation(df = colanno_df)

ComplexHeatmap::Heatmap(cellchat_mat_filtered, 
                        col = cellchat_mat_col_fun, 
                        show_row_names = FALSE, 
                        right_annotation = cellchat_row_anno,
                        show_column_names = FALSE,
                        top_annotation = cellchat_col_anno)

