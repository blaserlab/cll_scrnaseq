cellchat_val_dat <- cellchat_dat |> 
  pivot_longer(cols = where(is.numeric), names_to = "specimen", values_to = "prob") |> 
  left_join(bb_cellmeta(cds_main) |> 
              group_by(specimen, timepoint_merged_1, patient_type3) |> 
              summarise())

pathways <- c("MIF", "MIF")
sources <- c("CD8 TEM", "B memory")
targets <- c("CD8 TEM", "CD14 Mono")

cellchat_validation_plots <- pmap(.l = list(
  pway = pathways,
  src = sources,
  tgt = targets
),
.f = \(pway, src, tgt, dat = cellchat_val_dat) {
  dat |>
    filter(pathway_name == pway,
           source == src,
           target == tgt)  |>
    ggplot(aes(
      x = patient_type3,
      y = prob,
      color = patient_type3,
      fill = patient_type3
    )) +
    geom_jitter(pch = 21) +
    facet_wrap( ~ timepoint_merged_1) +
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
    labs(
      y = "Prob",
      color = NULL,
      fill = NULL,
      x = NULL,
      title = paste0(pway, " | ", src, " | ", tgt)
    ) +
    theme(axis.title.y.left = ggtext::element_markdown()) +
    theme(axis.text.x.bottom = element_blank()) +
    theme(axis.ticks.x.bottom = element_blank()) +
    theme(legend.position = "top",
          legend.justification = "center") +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(fill = guide_legend(ncol = 2, override.aes = list(size = 2)))
  
}) |> set_names(paste(pathways, sources, targets, sep = " | "))

