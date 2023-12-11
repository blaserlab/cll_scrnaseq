# bcell population proportions -----------------------------------
normalized_leiden_counts <- 
  bb_cellmeta(cds_main) |> 
  filter(partition_assignment == "B") %>%
  group_by(patient, 
           leiden_comparison_renamed, 
           specimen, 
           timepoint_merged_1, 
           patient_type2) %>%
  summarise(n = n()) %>%
  left_join(colData(cds_main) %>%
              as_tibble() %>%
              filter(partition_assignment == "B") |> 
              group_by(specimen) %>%
              summarise(specimen_total = n())) %>%
  mutate(overall_total = nrow(colData(cds_main))) %>%
  mutate(normalized_count = n*overall_total/specimen_total/2) %>%
  select(leiden_comparison_renamed, specimen, timepoint_merged_1, patient_type2, normalized_count)

cluster_proportion_ratio_plot <- normalized_leiden_counts %>%
  pivot_wider(names_from = leiden_comparison_renamed, values_from = normalized_count, values_fill = 1) %>%
  mutate(btk_to_other_ratio = (stressed + `inflammatory 1`+ `inflammatory 2` )/(other)) %>%
  mutate(log2_ratio = log2(btk_to_other_ratio)) %>%
  ggplot(mapping = aes(x = patient_type2, y = log2_ratio, color = patient_type2, fill = patient_type2)) +
  geom_jitter(shape = jitter_shape, size = jitter_size, stroke = jitter_stroke) +
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
  stat_compare_means(method = "wilcox", 
                     label = "p.signif", 
                     label.x.npc = "center", 
                     label.y = 16, 
                     show.legend = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.1))) +
  labs(y = "log<sub>2</sub>(inflammatory + stressed:other)", color = "Patient Type", fill = "Patient Type", x = NULL) +
  theme(axis.title.y.left = ggtext::element_markdown()) + 
  theme(axis.text.x.bottom = element_blank()) +
  theme(axis.ticks.x.bottom = element_blank()) +
  theme(legend.position = "bottom", legend.justification = "center")
