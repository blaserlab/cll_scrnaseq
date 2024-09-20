# comparing the ratio of cells in the leiden enrichment clusters in sensitive and resistant patients by timepoint

cluster_proportion_plot_cd16 <- bb_cellmeta(cds_main) %>%
  filter(seurat_l2_leiden_consensus == "CD16 Mono") |> 
  count(patient, sample, cd16_da_response, patient_type2, timepoint_merged_2) |> 
  pivot_wider(names_from = cd16_da_response, values_from = n, values_fill = 0) %>%
  mutate(pct = sensitive/(unenriched + resistant + sensitive)*100) |> 
  ggplot(mapping = aes(x = patient_type2, y = pct, color = patient_type2, fill = patient_type2)) +
  geom_jitter(shape = jitter_shape, size = jitter_size, stroke = jitter_stroke) +
  facet_wrap(facets = vars(timepoint_merged_2)) +
  scale_fill_manual(values = alpha(colour = experimental_group_palette, alpha = jitter_alpha_fill)) +
  scale_color_manual(values = alpha(colour = experimental_group_palette, alpha = jitter_alpha_color)) +
  theme(strip.background = element_blank()) +
  theme(panel.background = element_rect(color = "grey80")) +
  theme(legend.position = "none") +
  stat_summary(
    fun.data = data_summary_mean_se,
    color = summarybox_color,
    size = summarybox_size,
    width = summarybox_width,
    alpha = summarybox_alpha,
    geom = summarybox_geom, 
    show.legend = FALSE
  ) +
  stat_compare_means(method = "t.test", 
                     label = "p.signif", 
                     label.x.npc = "center",
                     label.y.npc = 0.9, 
                     show.legend = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.1))) +
  labs(y = "Percent IBS-enriched", 
       color = NULL, 
       fill = NULL, 
       x = NULL) +
  theme(axis.title.y.left = ggtext::element_markdown()) + 
  theme(axis.text.x.bottom = element_blank()) +
  theme(axis.ticks.x.bottom = element_blank()) +
  theme(legend.position = "top", 
        legend.justification = "center") +
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 2)))
cluster_proportion_plot_cd16
