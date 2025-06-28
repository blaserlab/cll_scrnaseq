# tcr diversity plot -------------------------------
tcr_diversity_plot <- bb_cellmeta(cds_main) |>
  group_by(patient_type3, sample, timepoint_merged_1) |> 
  summarise(diversity = mean(specimen_tcr_shannon, na.rm = TRUE)) |> 
  ggplot(aes(x = patient_type3, color = patient_type3, fill = patient_type3, y = diversity)) +
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
  ) +
  stat_compare_means(method = "wilcox",
                     label = "p.signif",
                     label.x.npc = 0.5,
                     label.y.npc = 0.8,
                     show.legend = FALSE) +
  theme(strip.background = element_blank()) +
  facet_grid(~timepoint_merged_1) +
  # labs(y = "TCR Diversity", color = "Patient Type", fill = "Patient Type", x = NULL) +
  labs(y = "TCR Diversity", color = NULL, fill = NULL, x = NULL) +
  theme(axis.text.x.bottom = element_blank()) +
  # theme(axis.ticks.x.bottom = element_blank()) +
  theme(legend.position = "bottom", legend.justification = "center") + 
  panel_border()
tcr_diversity_plot
