cd19_data <- clinical_flow_data |>
  filter(population == "CD19") |>
  mutate(patient_type = case_match(patient_type, "resistant" ~ "IBR", "sensitive" ~ "IBS")) |>
  mutate(timepoint_merged_2 = str_remove(timepoint_merged_2, "Timepoint ")) |>
  mutate(log_abs = log10(abs_mm3)) |> 
  mutate(patient = factor(patient), 
         patient_type = factor(patient_type),
         timepoint_merged_2 = factor(timepoint_merged_2))

library(lme4)
library(lmerTest)  # optional: p-values for fixed effects
library(emmeans)

ibs_fit <- lmer(log_abs ~ timepoint_merged_2 + (1 | patient), data = cd19_data |> filter(patient_type == "IBS"))
summary(ibs_fit)
anova(ibs_fit)   # test overall timepoint effect

ibr_fit <- lmer(log_abs ~ timepoint_merged_2 + (1 | patient), data = cd19_data |> filter(patient_type == "IBR"))
summary(ibr_fit)
anova(ibr_fit)   # test overall timepoint effect

pval_tbl <- bind_rows(
  emmeans(ibs_fit, pairwise ~ timepoint_merged_2, adjust = "tukey")$contrasts |> as_tibble() |> mutate(patient_type = "IBS"),
  emmeans(ibr_fit, pairwise ~ timepoint_merged_2, adjust = "tukey")$contrasts |> as_tibble() |> mutate(patient_type = "IBR")
) |> 
  mutate(contrast = str_remove_all(contrast, "timepoint_merged_2")) |> 
  separate(contrast, into = c("group1", "group2"), sep = " - ") |> 
  rstatix::add_significance(p.col = "p.value",)

cll_flow_fig <- ggplot(
  cd19_data,
  aes(
    x = timepoint_merged_2,
    y = log_abs,
    color = patient_type,
    fill = patient_type
  )
) +
  geom_point(shape = jitter_shape,
             size = jitter_size,
             stroke = jitter_stroke) +
  stat_summary(
    fun.data = mean_se,
    color = summarybox_color,
    size = summarybox_size,
    width = summarybox_width,
    alpha = summarybox_alpha,
    geom = summarybox_geom,
    show.legend = FALSE
  ) +
  ggpubr::geom_bracket(color = "black",
    aes(xmin = group1, xmax = group2, label = p.value.signif),
    data = pval_tbl, y.position = c(6, 7, 6.5, 6, 7, 6.5)
  ) +
  facet_wrap(~ patient_type) +
  scale_fill_manual(values = alpha(colour = experimental_group_palette, alpha = jitter_alpha_fill)) +
  scale_color_manual(values = alpha(colour = experimental_group_palette, alpha = jitter_alpha_color)) +
  theme(strip.background = element_blank()) +
  theme(panel.background = element_rect(color = "grey80")) +
  theme(legend.position = "none") +
  labs(y = "log<sub>10</sub> CD19<sup>+</sup> cells/mm<sup>3</sup>", x = "Timepoint") +
  theme(axis.title.y = ggtext::element_markdown()) + 
  ylim(c(0, 7.25))
