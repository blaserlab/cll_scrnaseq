bb_cellmeta(cds_main) |> glimpse()
SingleCellExperiment::mainExpName(cds_main) <- "Gene Expression"
violin_genes <- c("TLR9", "ZAP70", "STAT6", "IL21R")
bb_rowmeta(cds_main) |> filter(gene_short_name %in% violin_genes)
bb_rowmeta
p <- bb_gene_violinplot(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> 
                                filter(partition_assignment == "B")), 
                   variable = c("patient_type1"), 
                   genes_to_plot = violin_genes, 
                   include_jitter = F)
p[["data"]]
p + facet_grid(feature_label~timepoint_merged_1, scales = "free")

b_gb_data <- bb_genebubbles(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> 
                           filter(partition_assignment == "B")), 
                           genes = violin_genes, 
                           cell_grouping = c("patient_type1", "timepoint_merged_1"),
                           return_value = "data")

b_cell_genebubbles <- ggplot(b_gb_data, aes(x = patient_type1, y = gene_short_name, size = proportion, color = expression)) +
  geom_point() +
  scale_color_viridis_c() +
  scale_size_area() +
  facet_wrap(~timepoint_merged_1) +
  panel_border() +
  theme(strip.background = element_blank()) +
  labs(y = NULL, x = NULL)

save_plot(b_cell_genebubbles, filename = fs::path(network_out, "b_cell_genebubbles.png"))

bb_gene_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> 
                          filter(partition_assignment == "B")), 
             "IL21R", order = F)

mod4 <- bb_gene_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> 
                          filter(partition_assignment == "B")), 
             gene_or_genes = bb_rowmeta(cds_main) |> 
               select(feature_id, module_labeled) |> 
               filter(module_labeled == "Module 4"))

mod2 + facet_grid(patient_type1~timepoint_merged_1)


mod2[["data"]] |> glimpse()

ggplot(mod[["data"]], aes(x = timepoint_merged_1, y = value, fill = patient_type1)) +
  # geom_point() + 
  geom_violin()

module_enrichment$`Module 4`
