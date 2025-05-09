# 
# bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "B")), "dominant_related", facet_by = "sample") + panel_border()

bcr_resistant <- bb_cellmeta(cds_main) |>
  filter(partition_assignment == "B") |>
  filter(!is.na(cdr3s_aa)) |>
  filter(patient_type2 == "resistant") |> 
mutate(dominant_related = ifelse(dominant_related, "Malignant Clone", "non-Malignant Clone")) |> 
mutate(dominant_related = factor(dominant_related, levels = c("non-Malignant Clone", 
                                                                "Malignant Clone"))) |> 
  count(sample,
        patient,
        patient_type3,
        timepoint_merged_1,
        dominant_related) |> 
  ungroup() |> 
  ggplot(aes(x = patient, y = n, fill = dominant_related)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(patient_type3 ~ timepoint_merged_1, drop = TRUE, scales = "free_x", shrink = TRUE) +
  labs(y = "Fraction of sequenced BCRs", fill = NULL, title = NULL) + 
  theme(strip.background = element_blank()) + 
  # theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())  +
  scale_fill_manual(values = experimental_group_palette)

bcr_sensitive <- bb_cellmeta(cds_main) |>
  filter(partition_assignment == "B") |>
  filter(!is.na(cdr3s_aa)) |>
  filter(patient_type2 == "sensitive") |> 
  mutate(dominant_related = ifelse(dominant_related, "Malignant Clone", "non-Malignant Clone")) |> 
mutate(dominant_related = factor(dominant_related, levels = c("non-Malignant Clone", 
                                                                "Malignant Clone"))) |> 
  count(sample,
        patient,
        patient_type3,
        timepoint_merged_1,
        dominant_related) |> 
  ungroup() |> 
  ggplot(aes(x = patient, y = n, fill = dominant_related)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_grid(patient_type3 ~ timepoint_merged_1, drop = TRUE, scales = "free_x", shrink = TRUE) +
  labs(y = "Fraction of sequenced BCRs", fill = NULL, title = NULL) + 
  theme(strip.background = element_blank()) + 
  # theme(axis.text.x = element_text(angle = 30, hjust = 1))  +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())  +
  scale_fill_manual(values = experimental_group_palette)

design <- "
1
2
3"
bcr_plot <- patchwork::guide_area() + bcr_resistant + bcr_sensitive + patchwork::plot_layout(design = design, guides = "collect", axes = "collect")
  
  
bcr_plot
