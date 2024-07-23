bb_cellmeta(cds_main) |> glimpse()

bb_var_umap(cds_main, "dominant_related", facet_by = "sample")

blaseRdata::msigdb_genesets
mrd1_gsea_res <- fgsea::fgsea(pathways = bb_extract_msig(filter_list = list(ORGANISM = "Homo sapiens"), return_form = "name_list"), 
                              stats = pseudobulk_MRD1_BTK[[2]] |> 
                                filter(is.finite(log2FoldChange)) |> 
                                select(gene_short_name, log2FoldChange) |> 
                                deframe())
mrd2_gsea_res <- fgsea::fgsea(pathways = bb_extract_msig(filter_list = list(ORGANISM = "Homo sapiens"), return_form = "name_list"), 
                              stats = pseudobulk_MRD2_BTK[[2]] |> 
                                filter(is.finite(log2FoldChange)) |> 
                                select(gene_short_name, log2FoldChange) |> 
                                deframe())

mrd1_gsea_res |> as_tibble() |> arrange(padj) |> View()
mrd2_gsea_res |> as_tibble() |> arrange(padj) |> View()
cds_main_bcell_subpop_top_markers |> filter(cell_group == "MRD2")

bb_cellmeta(cds_main) |> glimpse()



bb_cellmeta(cds_main) |> 
  filter(partition_assignment == "B") |> 
  ggplot(aes(x = prealignment_dim1, y = prealignment_dim2, color = leiden_comparison_renamed)) + 
  geom_point(alpha = 0.1, size = 0.2) + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 0.8))) + 
  labs(color = NULL) + 
  facet_grid(patient_type2~timepoint_merged_2) + 
  panel_border()

bb_cellmeta(cds_main) |> 
  filter(partition_assignment == "B") |> 
  ggplot(aes(x = prealignment_dim1, y = prealignment_dim2, color = patient)) + 
  geom_point(alpha = 0.1, size = 0.2) + 
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 0.8))) + 
  labs(color = NULL) + 
  facet_grid(patient_type2~timepoint_merged_2) + 
  panel_border() + 
  scale_colour_brewer(palette = "Paired")


bb_gene_umap(cds_main, "CD19")


