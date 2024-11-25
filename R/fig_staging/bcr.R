bb_cellmeta(cds_main) |> glimpse()

bb_var_umap(cds_main, "dominant_related", facet_by = "sample")

bb_cellmeta(cds_main) |> 
  filter(!is.na(cdr3s_aa)) |> 
  count(patient, patient_type2, timepoint_merged_2, dominant_related, cdr3s_aa) |>
  ggplot(aes(x = patient, y = n, fill = dominant_related)) + 
  geom_bar(stat = "identity", position = "fill") + 
  facet_wrap(patient_type2~timepoint_merged_2, scales = "free", shrink = TRUE) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) + 
  labs(y = "Fraction of sequenced BCRs")