genes <- c("CD14", "FCGR3A", "CD8A", "GZMK")
cds_monos <- filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(seurat_l2_leiden_consensus %in% c("CD14 Mono", "CD16 Mono")))
colData(cds_monos)$grouping <- paste0(colData(cds_monos)$patient_type2, "_",
                                      colData(cds_monos)$timepoint_merged_2, "_",
                                      colData(cds_monos)$seurat_l2_leiden_consensus)


mono_genebubdat <- bb_genebubbles(obj = cds_monos, genes = genes, cell_grouping = "grouping", return_value = "data", ) |> 
  mutate(patient_type2 = str_extract(grouping, "resistant|sensitive"),
         timepoint_merged_2 = str_extract(grouping, "Timepoint 1|Timepoint 2|Timepoint 3"),
         seurat_l2_leiden_consensus = str_extract(grouping, "CD14 Mono|CD16 Mono")) |> 
  mutate(expression = ifelse(expression>2, 2, expression))

ggplot(mono_genebubdat, aes(x = patient_type2, y = gene_short_name, fill = expression, size = proportion)) + 
  geom_point(pch = 21, color = "black") + 
  theme_minimal_grid(font_size = 10) + 
  scale_fill_viridis_c() +
  facet_grid(seurat_l2_leiden_consensus~timepoint_merged_2)

cd14_pseudobulk_gsea_res_full |> View()
cd14_pseudobulk_gsea_res_full |> filter(pathway == "GSE9988_LOW_LPS_VS_CTRL_TREATED_MONOCYTE_UP") |> unnest(leadingEdge)
bb_gene_umap(cds_monos, "GZMK")  
bb_gene_umap(cds_monos, "CCL4")  

cd16_pseudobulk_gsea_res_full |> View()
cd16_pseudobulk_gsea_res_full |> filter(pathway == "AIZARANI_LIVER_C3_NK_NKT_CELLS_2") |> unnest(leadingEdge)
bb_gene_umap(cds_monos, "CCL5")  
cd16_pseudobulk_res$Result |> filter(log2FoldChange>0) |> arrange(padj) |> View()
