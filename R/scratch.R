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
