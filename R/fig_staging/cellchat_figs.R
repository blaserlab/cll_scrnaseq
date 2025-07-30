significant <-  cellchat_wilcox |> 
  filter(p < 0.05) |> 
  pull(rownames)

cellchat_mat <- cellchat_dat |> 
  pivot_longer(cols = where(is.numeric), names_to = "sample", values_to = "probability") |> 
  left_join(bb_cellmeta(cds_main) |> group_by(sample, patient_type3, timepoint_merged_1) |> summarise()) |> 
  group_by(rownames, source, target, interaction_name, pathway_name, patient_type3, timepoint_merged_1) |> 
  summarise(median_probability = median(probability)) |>
  ungroup() |> 
  mutate(type_timepoint = paste0(patient_type3, "_", timepoint_merged_1)) |> 
  select(-patient_type3, -timepoint_merged_1) |> 
  pivot_wider(names_from = type_timepoint, values_from = median_probability) |> 
  filter(!str_detect(target, "^B ")) |>
  filter(!str_detect(interaction_name, "^HLA")) |>
  filter(rownames %in% significant) |> 
  filter(str_detect(rownames, "CD99", negate = TRUE)) |> 
  select(rownames, where(is.numeric)) |> 
  bb_tbl_to_matrix()

cellchat_mat_filtered <- cellchat_mat[rowSums(cellchat_mat)>0.001,]
cellchat_sh <- SummarizedHeatmap(cellchat_mat_filtered, colOrder = c("IBR_1", "IBS_1", "IBR_2", "IBS_2", "IBR_3", "IBS_3"))
colData(cellchat_sh)$sample <- rownames(colData(cellchat_sh))
colData(cellchat_sh)$response <- str_extract(colData(cellchat_sh)$sample, "IBS|IBR")
colData(cellchat_sh)$timepoint <- str_extract(colData(cellchat_sh)$sample, "1|2|3")

rowmeta <- tibble(rownames = rownames(rowData(cellchat_sh))) |> left_join(cellchat_dat |> select(rownames, source, target, pathway_name))
blaseRtools::rowData(cellchat_sh)$source <- rowmeta$source
blaseRtools::rowData(cellchat_sh)$target <- rowmeta$target
blaseRtools::rowData(cellchat_sh)$pathway_name <- rowmeta$pathway_name

hm_pal <- brewer.pal(n = 10, name = "Paired")
names(hm_pal) <- c(as.character(rowmeta$source), as.character(rowmeta$target), rowmeta$pathway_name) |> unique()

cellchat_hm_fun <- function(hm, pal1, pal2) {
  p1 <- bb_plot_heatmap_main(hm) + scale_fill_distiller(palette = "Purples", direction = 1) + theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) + theme(plot.margin = unit(c(0, 0, 10, 0), "mm"))
  p3 <- bb_plot_heatmap_rowDendro(hm, side = "left") + scale_y_reverse()
  p5 <- bb_plot_heatmap_colData(hm,
                                vars = c("Response" = "response", "Timepoint" = "timepoint")) &
    scale_fill_manual(values = pal1) &
    scale_y_discrete(expand = expansion(add = 0))
  p6 <- bb_plot_heatmap_rowData(
    hm,
    vars = c(
      "Source" = "source",
      "Target" = "target",
      "Pathway" = "pathway_name"
    )
  ) &
    scale_x_discrete(position = "bottom", expand = expansion(add = 0)) &
    scale_fill_manual(values = pal2) & 
    theme(axis.text.x.bottom = element_text(angle = 30, hjust = 1)) + theme(plot.margin = unit(c(0, 0, 10, 0), "mm"))
    
  design <-
    "
#D#E
BACE
"
  wrap_plots(
    A = p1,
    B = p3,
    E = guide_area(),
    D = free(p5, type = "space"),
    C = free(p6, type = "space"),
    design = design,
    guides = "collect",
    heights = c(1, 10),
    widths = c(1, 5, 2, 2)
  )
  
}
cellchat_hm <- cellchat_hm_fun(cellchat_sh, experimental_group_palette, hm_pal)




int_name <- c("MIF_CD74_CD44", "ANXA1_FPR1")

sources <- c("CD8 TEM", "CD8 TEM")
targets <- c("CD14 Mono", "CD14 Mono")
titles <- c("MIF-CD74/CD44:  CD8 TEM \u2192 CD14 Mono",
            "ANXA1-FPR1:  CD8 TEM \u2192 CD14 Mono")

cellchat_validation_plots <- pmap(.l = list(
  int = int_name,
  src = sources,
  tgt = targets,
  tit = titles
),
.f = \(int, src, tgt, tit, dat = cellchat_val_dat) {
  dat |>
    filter(interaction_name == int,
           source == src,
           target == tgt)  |>
    ggplot(aes(
      x = patient_type3,
      y = prob,
      color = patient_type3,
      fill = patient_type3
    )) +
    geom_jitter(shape = jitter_shape, size = jitter_size, stroke = jitter_stroke, width = jitter_width) +
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
      method = "wilcox",
      label = "p.signif",
      label.x.npc = "center",
      label.y.npc = 0.9,
      show.legend = FALSE
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.1))) +
    labs(
      y = "Probability",
      color = NULL,
      fill = NULL,
      x = NULL,
      title = tit
    ) +
    theme(axis.title.y.left = ggtext::element_markdown()) +
    theme(axis.text.x.bottom = element_blank()) +
    theme(axis.ticks.x.bottom = element_blank()) +
    theme(legend.position = "top",
          legend.justification = "center") +
    theme(plot.title = element_text(hjust = 0.5)) +
    guides(fill = guide_legend(ncol = 2, override.aes = list(size = 2))) 
  
})

