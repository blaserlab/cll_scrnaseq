


cll_flow_fig <- clinical_flow_data |>
  filter(population == "CD19") |>
  mutate(patient_type = case_match(patient_type, "resistant" ~ "IBR", "sensitive" ~ "IBS")) |>
  mutate(timepoint_merged_2 = str_remove(timepoint_merged_2, "Timepoint ")) |>
  ggplot(
    aes(
      x = timepoint_merged_2,
      y = log10(abs_mm3),
      group = patient,
      color = patient_type,
      fill = patient_type
    )
  ) +
  geom_line() +
  geom_point(shape = jitter_shape,
             size = jitter_size,
             stroke = jitter_stroke) +
  facet_wrap( ~ patient_type) +
  scale_fill_manual(values = alpha(colour = experimental_group_palette, alpha = jitter_alpha_fill)) +
  scale_color_manual(values = alpha(colour = experimental_group_palette, alpha = jitter_alpha_color)) +
  theme(strip.background = element_blank()) +
  theme(panel.background = element_rect(color = "grey80")) +
  theme(legend.position = "none") +
  labs(y = "log<sub>10</sub> CD19<sup>+</sup> cells/mm<sup>3</sup>", x = "Timepoint") +
  theme(axis.title.y = ggtext::element_markdown())

cll_flow_fig

cll_flow_fig <- plot_grid(cll_flow_fig, labels = "B")

save_plot(
  plot = cll_flow_fig,
  filename = fs::path(network_out, "cll_flow_fig.pdf"),
  base_width = 5,
  base_height = 3
)
