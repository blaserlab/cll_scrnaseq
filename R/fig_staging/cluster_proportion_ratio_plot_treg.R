# comparing the ratio of Tregs to other T/NK cells

cluster_proportion_ratio_plot_treg <-
  bb_cellmeta(cds_main) |>
  filter(partition_assignment  %in% c("T", "NK")) |>
  count(patient,
        sample,
        seurat_l2_leiden_consensus,
        patient_type2,
        timepoint_merged_1) |>
  pivot_wider(names_from = seurat_l2_leiden_consensus,
              values_from = n,
              values_fill = 1) |>
  mutate(ratio = Treg / (`CD4 Naive` + `CD4 TCM` + `CD8 TEM` + NK)) |>
  mutate(log2_ratio = log2(ratio)) %>%
  ggplot(mapping = aes(
    x = patient_type2,
    y = log2_ratio,
    color = patient_type2,
    fill = patient_type2
  )) +
  geom_jitter(shape = jitter_shape,
              size = jitter_size,
              stroke = jitter_stroke) +
  facet_wrap(facets = vars(timepoint_merged_1)) +
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
    label.y = 1,
    show.legend = FALSE
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.1))) +
  labs(y = "log<sub>2</sub>(Treg:other)",
       color = NULL,
       fill = NULL,
       x = NULL) +
  theme(axis.title.y.left = ggtext::element_markdown()) +
  theme(axis.text.x.bottom = element_blank()) +
  theme(axis.ticks.x.bottom = element_blank()) +
  theme(legend.position = "top",
        legend.justification = "center") +
  guides(fill = guide_legend(ncol = 1, override.aes = list(size = 2)))
