
bb_gene_umap(cds_main, "MIR34AHG")
bb_gene_umap(cds_main, "BCL2")
bb_gene_umap(cds_main, "SIRT1")
bb_var_umap(cds_main, "leiden_mutant_enriched")
bb_gene_dotplot(cds_main, markers = c( "TP53", "MDM2", "CDKN1A", "MDM4"), group_cells_by = "leiden_mutant_enriched")

bb_gene_umap(cds_main, gene_or_genes = bb_rowmeta(cds_main) %>% 
               filter(gene_short_name %in% c( "TP53", "MDM2", "CDKN1A", "MDM4")) %>% 
               select(feature_id) %>% 
               mutate(gene_grouping = "P53_related"))


leiden_mutant_enriched_pseudobulk_tbl %>% filter(padj<0.05) %>% arrange(log2FoldChange) %>% View()
leiden_mutant_enriched_pseudobulk_tbl %>% filter(gene_short_name == "MIR34AHG")



pseudosample_table <- bb_cellmeta(cds_main) %>%
  group_by(sample, leiden_mutant_enriched) %>%
  summarise() %>%
  filter(!is.na(leiden_mutant_enriched))
pseudosample_table

test_pseudobulk_mf <- bb_pseudobulk_mf(cds = cds_main[, !is.na(colData(cds_main)$leiden_mutant_enriched)], pseudosample_table = pseudosample_table, design_formula = "~leiden_mutant_enriched")

test_gsea_res <- fgsea::fgsea(stats = test_pseudobulk_mf$Result %>% select(gene_short_name, log2FoldChange) %>% deframe(),
      pathways = bb_extract_msig(filter_list = list("ORGANISM" = "Homo sapiens"), return_form = "name_list"))

test_gsea_res %>% filter(padj < 0.05) %>% View()
bb_gene_umap(cds_main, "SUMO2")


blaseRdata::msigdb_geneset_metadata

test_pseudobulk_mf$Result %>% filter(gene_short_name == "BIRC5")
test_pseudobulk_mf$Heatmap
leiden_mutant_enriched_tm <- top_markers(cds_main[rowData(cds_main)$gene_short_name %notin% blaseRdata::hg38_remove_genes,
                                                  !is.na(colData(cds_main)$leiden_mutant_enriched)], 
                                         group_cells_by = "leiden_mutant_enriched", 
                                         genes_to_test_per_group = 50, 
                                         cores = 40)

bb_cellmeta(cds_main)
bb_var_umap(cds_main, "del17p")
bb_var_umap(cds_main, "density", facet_by = "del17p")
bb_var_umap(cds_main, "density", facet_by = c("del17p", "patient"), cols = vars(del17p), rows = vars(patient))

bb_var_umap(cds_main, "density", facet_by = "del11q")
bb_var_umap(cds_main, "density", facet_by = "del13p")
bb_var_umap(cds_main, "density", facet_by = "complex_karyotype")
bb_var_umap(cds_main, "density", facet_by = "patient_type")

bb_var_umap(cds_main, "leiden_mutant_enriched")
bb_var_umap(cds_main, "density", facet_by = "btk_type")
bb_gene_umap(cds_main, "TP53")


bb_monocle_regression(cds_main, gene_or_genes = "TP53", form = "~leiden_mutant_enriched")
bb_monocle_regression(cds_main, gene_or_genes = "MDM4", form = "~leiden_mutant_enriched")




test_pseudobulk_mf$Result %>% filter(gene_short_name == "TP53")



sample_mut_table <- bb_cellmeta(cds_main) %>%
  group_by(
    specimen,
    patient,
    timepoint,
    gender,
    IGHV_status,
    complex_karyotype,
    del17p,
    del11q,
    del13p,
    tri12,
    patient_type,
    leiden_mutant_enriched
  ) %>%
  summarise(n = n()) %>%
  filter(!is.na(leiden_mutant_enriched)) %>%
  pivot_wider(names_from = leiden_mutant_enriched,
              values_from = n,
              values_fill = 0) %>%
  rename(mutant_enriched = `TRUE`) %>%
  filter(mutant_enriched > 10) %>%
  rename(not_mutant_enriched = `FALSE`) %>%
  filter(not_mutant_enriched >= 10) %>%
  mutate(sample_mut_pct = mutant_enriched/(not_mutant_enriched + mutant_enriched)*100) %>%
  mutate(log2_sample_mut_pct = log2(sample_mut_pct))
sample_mut_table

model <- lm(formula = log2_sample_mut_pct ~ del17p, 
            data = sample_mut_table)
summary(model)
warnings()
