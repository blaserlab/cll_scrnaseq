t_cell_subpop_genes <-
  c("IFNG", "TBX21", "EOMES", "IL2RA", "FOXP3", "CCR7")



tcell_subpop_umap <-
  bb_var_umap(
    cds_main[, colData(cds_main)$partition_assignment == "T"],
    "tcell_subpops",
    overwrite_labels = T,
    group_label_size = 3,
    foreground_alpha = 0.2
  )

tcell_subpop_gene_umap <-
  bb_gene_umap(cds_main[, colData(cds_main)$partition_assignment == "T"],
               t_cell_subpop_genes,
               color_legend_title = "Expression",
               alpha = 0.6)



colData(cds_main)$timepoint_merged_safe <-
  recode(
    colData(cds_main)$timepoint_merged,
    "3yrs|btk_clone" = "3yrs_btk_clone",
    "5yrs|relapse" = "5yrs_relapse"
  )



tcell_subpop_gene_dotplot <-
  bb_gene_dotplot(
    cds_main[, colData(cds_main)$partition_assignment == "T"],
    markers = t_cell_subpop_genes,
    gene_ordering = rev(t_cell_subpop_genes),
    group_cells_by = "multifactorial",
    group_ordering = tribble(
      ~ aesthetic,
      ~ variable,
      ~ value,
      ~ level,
      "facet",
      "timepoint_merged_safe",
      "baseline",
      1,
      "facet",
      "timepoint_merged_safe",
      "3yrs_btk_clone",
      2,
      "facet",
      "timepoint_merged_safe",
      "5yrs_relapse",
      3,
      "axis",
      "patient_type",
      "BTK",
      1,
      "axis",
      "patient_type",
      "MRD",
      2
    ),
    max.size = 6,
    colorscale_name = "Expression",
    sizescale_name = "Proportion"
  ) +
  labs(x = NULL, y = NULL)


major_tcell_clusters <- colData(cds_main) %>%
  as_tibble() %>%
  filter(partition_assignment == "T") %>%
  group_by(seurat_celltype_l2) %>%
  summarise(n = n()) %>%
  filter(n > 100) %>%
  pull(seurat_celltype_l2)

major_tcell_cluster_umap <-
  bb_var_umap(
    cds_main[, colData(cds_main)$partition_assignment == "T" &
               colData(cds_main)$seurat_celltype_l2 %in% major_tcell_clusters],
    "seurat_celltype_l2",
    overwrite_labels = T,
    group_label_size = 3,
    foreground_alpha = 0.2
  )


treg_ratio_tbl <- left_join(
  colData(cds_main) %>%
    as_tibble() %>%
    filter(partition_assignment == "T") %>%
    group_by(specimen, patient, timepoint_merged, patient_type) %>%
    summarise(n_total_t = n()),
  colData(cds_main) %>%
    as_tibble() %>%
    filter(leiden == "26") %>%
    group_by(specimen)  %>%
    summarise(n_treg = n())
) %>%
  mutate(n_treg = replace_na(n_treg, 1)) %>%
  mutate(treg_pct = n_treg / (n_total_t) * 100)



treg_pct_plot <-
  ggplot(
    treg_ratio_tbl,
    mapping = aes(
      x = patient_type,
      y = treg_pct,
      color = patient_type,
      fill = patient_type
    )
  ) +
  geom_jitter(width = jitter_width,
              size = jitter_size,
              shape = jitter_shape) +
  scale_color_manual(values = experimental_group_palette) +
  scale_fill_manual(values = alpha(alpha = 0.4, colour = experimental_group_palette)) +
  coord_trans(y = "log10", clip = "off") +
  scale_y_continuous(breaks = c(60, 20, 6, 2, 0.6, 0.2) / 2) +
  annotation_logticks(scaled = F,
                      sides = "l",
                      outside = F) +
  labs(y = "Percent Treg", x = NULL) +
  theme(legend.position = "none") +
  facet_wrap(facets = vars(timepoint_merged)) +
  stat_summary(fun.data = blaseRtools::data_summary_mean_se, 
               geom = "crossbar", 
               width = 0.3, 
               show.legend = FALSE) +
  stat_compare_means(method = "t.test",
                     label = "p.signif",
                     label.x.npc = 0.5) +
  theme(strip.background = element_blank())


bb_cellmeta(cds_main) |> glimpse()
bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "T")), "tcr_is_cancer", facet_by = c("value", "patient_type"), cols = vars(value), rows = vars(patient_type))

bb_cellmeta(cds_main) |> filter(!is.na(tcr_is_cancer)) |> group_by(specimen, patient_type, timepoint_merged, tcr_is_cancer) |> summarise(n = n()) |> 
  mutate(tcr_is_cancer_chr = ifelse(tcr_is_cancer, "cancer_related", "not_cancer_related")) |> 
  select(-tcr_is_cancer) |> 
  pivot_wider(names_from = tcr_is_cancer_chr, values_from = n, values_fill = 0) |> 
  mutate(cancer_pct = cancer_related/(not_cancer_related + cancer_related)*100) |> 
  ggplot(aes(x = patient_type, y = cancer_pct)) +
  geom_jitter() +
  facet_wrap(~timepoint_merged)



shannon <- function(x) {
  p <- x/sum(x)
  -sum(p*log(p))
}

tcr_diversity_plot <- bb_cellmeta(cds_main) |> 
  filter(!is.na(tcr_clone_copies)) |> 
  group_by(specimen, tcr_clonotype_id) |> 
  summarise(n = n()) |> 
  group_by(specimen) |> 
  summarize(diversity = shannon(n)) |> 
  left_join(bb_cellmeta(cds_main) |> 
              group_by(specimen, timepoint_merged, patient_type) |> 
              summarise()) |> 
  ggplot(aes(x = patient_type, y = diversity, color = patient_type, fill = patient_type)) +
  geom_jitter(width = jitter_width,
              size = jitter_size,
              shape = jitter_shape) +
  scale_color_manual(values = experimental_group_palette) +
  scale_fill_manual(values = alpha(alpha = 0.4, colour = experimental_group_palette)) +
  facet_wrap(~timepoint_merged) +
  stat_summary(fun.data = blaseRtools::data_summary_mean_se, 
               geom = "crossbar", 
               width = 0.3, 
               show.legend = FALSE) +
  stat_compare_means(method = "t.test",
                     label = "p.signif",
                     label.x.npc = 0.5) 
  



