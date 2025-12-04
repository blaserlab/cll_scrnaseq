cp_data <- bb_cellmeta(cds_main) %>%
  filter(partition_assignment  == "B") |>
  count(patient,
        sample,
        leiden_comparison_renamed,
        patient_type3,
        timepoint_merged_1) |>
  pivot_wider(names_from = leiden_comparison_renamed,
              values_from = n,
              values_fill = 1) %>%
  mutate(percent = 100 * (
    `Activated BCR` / (`Inflammatory 1` + `Inflammatory 2` + Stressed + `Activated BCR`)
  )) |>
  select(patient, patient_type3, time = timepoint_merged_1, percent)

padj_data <- cp_data |> 
  group_by(time) |> 
  rstatix::t_test(percent ~ patient_type3, p.adjust.method = "none", var.equal = FALSE) |> 
  mutate(padj_BH = p.adjust(p, method = "BH"))
 

cluster_proportion_ratio_plot <- 
  ggplot(cp_data, mapping = aes(x = patient_type3, y = percent, color = patient_type3, fill = patient_type3)) +
  geom_jitter(shape = jitter_shape, size = jitter_size, stroke = jitter_stroke) +
  facet_wrap(facets = vars(time)) +
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
  stat_stars_facet(
    aes(time_panel = time),             # <— key line
    stat_df = padj_data,
    time_col = "time", 
    p_col = "padj_BH",
    cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
    symbols   = c("****", "***", "**", "*", "ns"),
    size  = 4, 
    vjust = 0,
    npc_x = 0.5,
    npc_y = 1.05
  ) + 
  scale_y_continuous(breaks = c(0, 25, 50, 75, 100),
                     expand = expansion(mult = c(0.1))) +
  labs(y = "Percent Activated BCR", 
       color = NULL, 
       fill = NULL, 
       x = NULL) +
  theme(axis.title.y.left = ggtext::element_markdown()) + 
  theme(axis.text.x.bottom = element_blank()) +
  theme(axis.ticks.x.bottom = element_blank()) +
  theme(legend.position = "top", 
        legend.justification = "center") +
  guides(fill = guide_legend(ncol = 2, override.aes = list(size = 2)))
cluster_proportion_ratio_plot
