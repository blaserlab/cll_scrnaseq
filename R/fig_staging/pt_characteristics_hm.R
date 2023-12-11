pt_char <- bb_cellmeta(cds_main) |>
  group_by(
    specimen,
    patient,
    patient_type2,
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


# pt_char_mat <- pt_char |>
#   select(specimen, where(is.integer)) |>
#   bb_tbl_to_matrix() |>
#   apply(1, \(x) x / sum(x) * 100)

pt_char_mat <- pt_char |> 
  select(specimen, where(is.integer)) |> 
  bb_tbl_to_matrix() |> 
  t()

pt_char_mat <- log10(pt_char_mat + 1)

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

pt_char_anno <-
  ComplexHeatmap::HeatmapAnnotation(
    df = pt_char_anno_df,
    which = "column",
    annotation_legend_param = list(labels_gp = gpar(fontsize = 8),title_gp = gpar(fontsize = 10, font = 2)),
    annotation_name_gp = gpar(fontsize = 10),
    annotation_label = c(
      "patient",
      "response",
      "timepoint",
      "gender",
      "age",
      "IGHV status",
      "CK",
      "del17p",
      "del11q",
      "del13p",
      "tri12",
      "VAF: t2",
      "VAF: t3"
    ),
    col = list(
      IGHV_status = c("M" = "black", "U" = "white"),
      complex_karyotype = c(
        "yes" = "black",
        "no" = "white",
        "unknown" = "grey80"
      ),
      timepoint_merged_1 = timepoint_palette,
      del17p = c(`FALSE` = "white", `TRUE` = "black"),
      del11q = c(`FALSE` = "white", `TRUE` = "black"),
      gender = c("F" = "pink", "M" = "blue"),
      del13p = c(`FALSE` = "white", `TRUE` = "black"),
      tri12 = c(`FALSE` = "white", `TRUE` = "black"),
      patient_type2 = experimental_group_palette,
      patient = pt_colors,
      age_at_ibr_start = col_fun_pt_age,
      btk_clone_vaf_pct = col_fun_btk_clone_vaf,
      relapse_vaf_pct = col_fun_relapse_vaf
      
    )
  )

col_fun_pt_char <- circlize::colorRamp2(breaks = c(min(pt_char_mat),max(pt_char_mat)), colors = c("transparent", "purple4"))



pt_characteristics_hm <- grid.grabExpr(draw(ComplexHeatmap::Heatmap(
  pt_char_mat,
  col = col_fun_pt_char,
  name = "Log10 Cells",
  column_dend_side = "bottom",
  show_column_names = FALSE,
  top_annotation = pt_char_anno,
  row_names_gp = gpar(fontsize = 10), 
  heatmap_legend_param = list(direction = "horizontal"),
  
), merge_legend = FALSE, 
heatmap_legend_side = "bottom"), wrap = TRUE, width = unit(5, "in"), height = unit(5, "in")) 


