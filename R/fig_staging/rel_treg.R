
rel_treg <- bb_cellmeta(cds_main) |>
  filter(partition_assignment  %in% c("T")) |>
  count(patient,
        sample,
        seurat_l2_leiden_consensus,
        patient_type2,
        timepoint_merged_1,
        timepoint_merged_2) |>
  pivot_wider(names_from = seurat_l2_leiden_consensus,
              values_from = n,
              values_fill = 1) |>
  # mutate(fraction_treg = Treg / (`CD4 Naive` + `CD4 TCM` + `CD8 TEM` + NK + Treg)) |> 
  mutate(fraction_treg = Treg / (`CD4 Naive` + `CD4 TCM` + `CD8 TEM` + Treg)) |> 
  left_join(clinical_flow_data |> 
              filter(population == "CD3"),
            by = join_by(patient, timepoint_merged_2)) |> 
  mutate(abs_treg = abs_mm3 * fraction_treg) |> 
  ggplot(aes(x = patient_type, y = fraction_treg, fill = patient_type, color = patient_type)) +
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
  stat_compare_means(method = "t.test",
                     label = "p.signif",
                     label.x.npc = 0.5,
                     label.y.npc = 0.9,
                     show.legend = FALSE) +
  theme(strip.background = element_blank()) +
  facet_wrap(~timepoint_merged_1, scales = "free") + 
  labs(y = "Cells/mm<sup>3</sup>", 
       color = NULL, 
       fill = NULL, 
       x = NULL) +
  theme(axis.title.y.left = ggtext::element_markdown()) + 
  theme(axis.text.x.bottom = element_blank()) +
  theme(axis.ticks.x.bottom = element_blank()) +
  theme(legend.position = "top", 
        legend.justification = "center") +
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 2))) +
  panel_border()
rel_treg
