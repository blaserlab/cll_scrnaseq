plotfun <-
  function(data, title) {
    ggplot(data, aes(
      x = V1,
      y = V2,
      size = size,
      fill = parentTerm
    )) +
      geom_point(pch = 21,
                 color = "black", 
                 aes(alpha = score)) +
      scale_fill_brewer(palette = "Set1") +
      guides(alpha = guide_legend(order = 2, override.aes=list(shape = 19))) +
      guides(fill = guide_legend(order = 1, override.aes=list(size = 3)))+
      guides(size = guide_legend(order = 3)) +
      labs(x = "PCoA 1", 
           y = "PCoA 2", 
           fill = "Parent GO-Term",
           size = "Term Size",
           alpha = "-log10 P", 
           title = title) +
      theme(legend.box = "horizontal") +
      theme(aspect.ratio = 0.9) +
      theme(plot.title = element_text(hjust = 0.5))
      
  }

# stressed
cd14_sensitive_goscatter <- bb_goscatter(simMatrix = monocyte_goterm_enrichment$cd14_ibr_sensitive$gosummary_0.8$simMatrix, 
             reducedTerms = monocyte_goterm_enrichment$cd14_ibr_sensitive$gosummary_0.8$reducedTerms)[["data"]] |> 
  as_tibble() |> 
  # mutate(parentTerm = recode(parentTerm, "biological process involved in interspecies interaction between organisms" = "interspecies interaction")) |> 
  filter(size >5) |> 
  plotfun(title = "cd14_sensitive") 

cd14_sensitive_goscatter

cd14_resistant_goscatter <- bb_goscatter(simMatrix = monocyte_goterm_enrichment$cd14_ibr_resistant$gosummary_0.8$simMatrix, 
             reducedTerms = monocyte_goterm_enrichment$cd14_ibr_resistant$gosummary_0.8$reducedTerms)[["data"]] |> 
  as_tibble() |> 
  # mutate(parentTerm = recode(parentTerm, "biological process involved in interspecies interaction between organisms" = "interspecies interaction")) |> 
  filter(size >5) |> 
  plotfun(title = "cd14_resistant") 

cd14_resistant_goscatter


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


