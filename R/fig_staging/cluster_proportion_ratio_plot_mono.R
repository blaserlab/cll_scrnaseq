# comparing the ratio of cells in the leiden enrichment clusters in sensitive and resistant patients by timepoint
bb_cellmeta(cds_main) |> glimpse()
cluster_proportion_cd14_data <- bb_cellmeta(cds_main) %>%
  filter(seurat_l2_leiden_consensus == "CD14 Mono") |> 
  count(patient, sample, cd14_louvain_da_response, patient_type2, timepoint_merged_1) |> 
  pivot_wider(names_from = cd14_louvain_da_response, values_from = n, values_fill = 0) %>%
  mutate(pct = sensitive/(unenriched + resistant + sensitive)*100) |> 
  mutate(cluster = "CD14 Mono")



cluster_proportion_cd16_data <- bb_cellmeta(cds_main) %>%
  filter(seurat_l2_leiden_consensus == "CD16 Mono") |> 
  count(patient, sample, cd16_da_response, patient_type2, timepoint_merged_1) |> 
  pivot_wider(names_from = cd16_da_response, values_from = n, values_fill = 0) %>%
  mutate(pct = sensitive/(unenriched + resistant + sensitive)*100) |> 
  mutate(cluster = "CD16 Mono")
  
cluster_proportion_mono_data <- bind_rows(cluster_proportion_cd14_data, cluster_proportion_cd16_data)
  
cluster_proportion_plot_mono <-
  ggplot(cluster_proportion_mono_data, 
         mapping = aes(
    x = patient_type2,
    y = pct,
    color = patient_type2,
    fill = patient_type2
  )) +
  geom_jitter(shape = jitter_shape,
              size = jitter_size,
              stroke = jitter_stroke) +
  facet_grid(cluster ~ timepoint_merged_1, scales = "free", axes = "all_x") +
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
  stat_compare_means(
    method = "t.test",
    label = "p.signif",
    label.x.npc = "center",
    label.y.npc = 0.9,
    show.legend = FALSE
  ) +
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
  guides(fill = guide_legend(ncol = 2, override.aes = list(size = 2)))
cluster_proportion_plot_mono
