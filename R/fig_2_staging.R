seurat_l2_leiden_consensus <- bb_cellmeta(cds_main) |> 
  group_by(leiden, seurat_celltype_l2) |> 
  summarise(n = n()) |> 
  slice_max(order_by = n, n = 1) |> 
  select(leiden, seurat_celltype_l2) |> 
  deframe()

colData(cds_main)$seurat_l2_leiden_consensus <- recode(colData(cds_main)$leiden, !!!seurat_l2_leiden_consensus)
colData(cds_main)$seurat_l2_leiden_consensus <- as.character(colData(cds_main)$seurat_l2_leiden_consensus)


tcell_subpop_umap <-
  bb_var_umap(
    cds_main[, colData(cds_main)$partition_assignment == "T"],
    "seurat_l2_leiden_consensus",
    overwrite_labels = T,
    group_label_size = 3,
    foreground_alpha = 0.2
  )

t_cell_subpop_genes <-
  c("IFNG", "TBX21", "EOMES", "IL2RA", "FOXP3", "CCR7")

subpop_bubdat <- bb_genebubbles(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "T")
), genes = t_cell_subpop_genes, 
cell_grouping = c("seurat_l2_leiden_consensus"), return_value = "data")

tcell_subpop_genebub <-
  ggplot(
    subpop_bubdat,
    aes(
      x = seurat_l2_leiden_consensus,
      y = gene_short_name,
      color = expression,
      size = proportion
    )
  ) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() +
  labs(x = NULL, y = NULL) +
  theme(axis.text.y = element_text(face = "italic"))

tcell_density_umap <- bb_var_umap(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "T")
), "density", facet_by = c("timepoint_merged", "patient_type"), cols = vars(timepoint_merged), rows = vars(patient_type), legend_pos = "bottom") +
  theme(legend.justification = "center") +
  panel_border()

bb_cellmeta(cds_main) |> glimpse()

bb_cluster_representation(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "T")
), cluster_var = "seurat_l2_leiden_consensus", class_var = "patient_type", experimental_class = "BTK", control_class = "MRD", return_value = "plot")

treg_ratio_tbl <- left_join(
  colData(cds_main) %>%
    as_tibble() %>%
    filter(partition_assignment == "T") %>%
    group_by(specimen, patient, timepoint_merged, patient_type) %>%
    summarise(n_total_t = n()),
  colData(cds_main) %>%
    as_tibble() %>%
    filter(seurat_l2_leiden_consensus == "Treg") %>%
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

exhaustion_genes <- c("TIGIT", "PDCD1", "LAG3", "CD160", "CD244")
exh_bubdat <- bb_genebubbles(filter_cds(
  cds_main,
  cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "T")
), genes = exhaustion_genes, 
cell_grouping = c("seurat_l2_leiden_consensus", "timepoint_merged", "patient_type"), return_value = "data")

exh_genebub <- ggplot(exh_bubdat, aes(x = seurat_l2_leiden_consensus, y = gene_short_name, color = expression, size = proportion)) +
  geom_point() +
  scale_size_area() +
  scale_color_viridis_c() + 
  facet_grid(patient_type ~ timepoint_merged) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1)) +
  panel_border() +
  labs(x = NULL, y = NULL) +
  theme(strip.background = element_blank()) +
  theme(axis.text.y = element_text(face = "italic"))
exh_genebub


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
                     label.x.npc = 0.5)  +
  theme(legend.position = "none") +
  theme(strip.background = element_blank())

tcr_diversity_plot
  
bb_var_umap(
  filter_cds(
    cds_main,
    cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "T")
  ),
  "tcr_clone_proportion",
  facet_by = c("patient_type", "timepoint_merged"),
  cols = vars(timepoint_merged),
  rows = vars(patient_type)
)




# Differential abundance  http://bioconductor.org/books/3.13/OSCA.multisample/differential-abundance.html
leiden_consensus_counts <- bb_cellmeta(cds_main) |> 
  count(sample, seurat_l2_leiden_consensus) |> 
  filter(seurat_l2_leiden_consensus %in% c("CD8 TEM", "CD4 TCM", "CD4 Naive", "Treg")) |>
  pivot_wider(names_from = "sample", values_from = "n", values_fill = 0) |> 
  bb_tbl_to_matrix()

# dge_list <- edgeR::DGEList(counts = leiden_consensus_counts, samples = bb_cellmeta(cds_main) |> group_by(specimen, patient_type, timepoint_merged) |> summarise())
# design <- model.matrix(~factor(timepoint_merged) + factor(patient_type, levels = c("MRD", "BTK")), dge_list$samples)
# dge_list <- estimateDisp(dge_list, design, trend="none")
# fit <- glmQLFit(dge_list, design, robust=TRUE, abundance.trend=FALSE)
# res <- glmQLFTest(fit, coef=ncol(design))
# res_top_tags <- topTags(res, n = Inf)
# cluster_enrichment_barchart <- as_tibble(res_top_tags@.Data[[1]], rownames = "cluster") |> 
#   mutate(enriched = ifelse(logFC > 0, "BTK", "MRD")) |> 
#   mutate(sig = case_when(PValue < 0.05 & PValue >= 0.01 ~ "*",
#                          PValue < 0.01 & PValue >= 0.001 ~ "**",
#                          PValue < 0.001 & PValue >= 0.0001 ~ "***",
#                          PValue >= 0.05 ~ ""
#                          )) |> 
#   # filter(cluster %in% c("CD8 TEM", "CD4 TCM", "CD4 Naive", "Treg")) |> 
#   mutate(cluster = fct_reorder(cluster, logFC)) |> 
#   ggplot(aes(x = cluster, y = logFC, color = enriched, fill = enriched, label = sig)) +
#   geom_col() +
#   labs(x = NULL, y = "Log-fold Enrichment BTK:MRD", color = NULL, fill = NULL) +
#   geom_text(color = "black", nudge_y = 0.05) +
#   theme(legend.position = "none")
# cluster_enrichment_barchart
# 
# bb_var_umap(filter_cds(
#   cds_main,
#   cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "T")
# ), "density", facet_by = "patient")
