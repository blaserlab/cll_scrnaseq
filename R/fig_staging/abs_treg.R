abs_treg_data <- bb_cellmeta(cds_main) |>
  filter(partition_assignment  %in% c("T")) |>
  count(patient,
        sample,
        seurat_l2_leiden_consensus,
        patient_type3,
        timepoint_merged_1,
        timepoint_merged_2) |>
  pivot_wider(names_from = seurat_l2_leiden_consensus,
              values_from = n,
              values_fill = 1) |>
  mutate(fraction_treg = Treg / (`CD4 Naive` + `CD4 TCM` + `CD8 TEM` + Treg)) |>
  left_join(clinical_flow_data |> 
              filter(population %in% c("CD3")),
            by = join_by(patient, timepoint_merged_2)) |> 
  mutate(abs_treg = abs_mm3 * fraction_treg) |> 
  select(patient, patient_type3, time = timepoint_merged_1, abs_treg)

treg_padj_data <- abs_treg_data |> 
  group_by(time) |> 
  rstatix::wilcox_test(abs_treg ~ patient_type3, p.adjust.method = "none") |> 
  mutate(padj_BH = p.adjust(p, method = "BH"))

abs_treg <-  
  ggplot(abs_treg_data, aes(x = patient_type3, y = abs_treg, fill = patient_type3, color = patient_type3)) +
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
  stat_stars_facet(
    aes(time_panel = time),             # <— key line
    stat_df = treg_padj_data,
    time_col = "time", 
    p_col = "padj_BH",
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
    symbols   = c("****", "***", "**", "*", "ns"),
    size  = 4, 
    vjust = 0,
    npc_x = 0.5,
    npc_y = 0.95
  ) + 
  theme(strip.background = element_blank()) +
  facet_wrap(~time, scales = "free") + 
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
abs_treg
