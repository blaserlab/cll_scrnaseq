cellchat_flow_figs <- pmap(.l = list(
  g = c("CD8 T cells", "CD14+"),
  m = c("MIF", "LGALS"),
  t = c("CD8 T cell MIF", "CD14 Mono LGALS")
), .f = \(g, m, t, dat = rebuttal_flow_data) {
  padj_data <- dat |>
    filter(gate == g, marker == m) |>
    mutate(timepoint = case_match(timepoint, "BL" ~ "1", "TP2" ~ "2", "TP3" ~ "3")) |> 
    group_by(marker, timepoint) |>
    rstatix::t_test(
      gmean ~ patient_type,
      p.adjust.method = "none",
      var.equal = FALSE,
      alternative = "g"
    ) |>
    mutate(padj_BH = p.adjust(p, method = "BH"))
  
  p <- dat |>
    filter(gate == g, marker == m) |>
    mutate(timepoint = case_match(timepoint, "BL" ~ "1", "TP2" ~ "2", "TP3" ~ "3")) |> 
    ggplot(aes(
      x = patient_type,
      y = gmean,
      color = patient_type,
      fill = patient_type
    )) +
    geom_jitter(
      shape = jitter_shape,
      size = jitter_size,
      width = jitter_width,
      stroke = jitter_stroke
    ) +
    facet_wrap( ~ timepoint) +
    scale_fill_manual(values = alpha(colour = experimental_group_palette, alpha = jitter_alpha_fill)) +
    scale_color_manual(values = alpha(colour = experimental_group_palette, alpha = jitter_alpha_color)) +
    theme(strip.background = element_blank()) +
    theme(legend.position = "none") +
    stat_summary(
      fun.data = data_summary_mean_se,
      color = summarybox_color,
      size = summarybox_size,
      width = summarybox_width,
      alpha = summarybox_alpha,
      geom = summarybox_geom,
      show.legend = FALSE
    )  +
    stat_stars_facet(
      aes(time_panel = timepoint),
      # <— key line
      stat_df = padj_data,
      time_col = "timepoint",
      p_col = "padj_BH",
      cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, Inf),
      symbols   = c("****", "***", "**", "*", "ns"),
      size  = 4,
      vjust = 0,
      npc_x = 0.5,
      npc_y = 0.95
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.1))) +
    labs(
      y = "Mean Fluorescence",
      color = NULL,
      fill = NULL,
      x = NULL,
      title = t
    ) +
    theme(axis.title.y.left = ggtext::element_markdown()) +
    theme(axis.text.x.bottom = element_blank()) +
    theme(axis.ticks.x.bottom = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5)) +
    theme(legend.position = "top",
          legend.justification = "center") +
    panel_border()
  return(list(padj = padj_data, plot = p))
  
  
}) |> set_names(c("MIF", "LGALS"))
# 
# cellchat_flow_figs$MIF$padj
# cellchat_flow_figs$MIF$plot
# 
# cellchat_flow_figs$LGALS$padj
# cellchat_flow_figs$LGALS$plot
