# table S1:  sequencing metrics----------------------
summarized_sequencing_metrics |> 
  rename(sample = cds_name) |> 
  left_join(bb_cellmeta(cds_main) |> count(sample, patient, patient_type2, timepoint_merged_1)) |> 
  select(Patient = patient, `Ibrutinib Response` = patient_type2, Timepoint = timepoint_merged_1, `Pre-filter Cell Count` = cds_dim_cells, `Mean Reads per Cell`, `Fraction Reads in Cells`,  `Post-filter Cell Count` = n) |> 
  write_csv(fs::path(paper_tables, "table_S1.csv"))

# table S2:  gene modules
bb_rowmeta(cds_main) |> 
  select(ensembl_id = id, gene_short_name, module) |> 
  filter(!is.na(module)) |> 
  arrange(module) |> 
  write_csv(fs::path(paper_tables, "table_S2.csv"))

# table S3: tcell leiden enrichment dge
tcell_leiden_enrichment_pseudobulk_res$Result |> 
  arrange(padj) |> 
  write_csv(fs::path(paper_tables, "table_S3.csv"))
  

# table S4:  T cell leiden enrichment gsea results ------------------------------
tcell_leiden_enrichment_gsea_res |> 
  mutate(leadingEdge = paste0(leadingEdge)) |> 
  arrange(padj) |> 
  write_csv(fs::path(paper_tables, "table_S4.csv"))

# table S5:  clinical flow data ------------------------------
clinical_flow_data |> 
  mutate(patient_type = case_match(patient_type, "sensitive" ~ "IBS", "resistant" ~ "IBR")) |> 
  select(patient, patient_type, timepoint = timepoint_merged_2, population, abs_mm3) |> 
  write_csv(fs::path(paper_tables, "table_S5.csv"))

# table S6:  cd14 monocyte DGE --------------------
  
cd14_pseudobulk_res$Result |> 
  arrange(padj) |> 
  write_csv(fs::path(paper_tables, "table_S6.csv"))

# table S7:  cd17 monocyte DGE --------------------
  
cd16_pseudobulk_res$Result |> 
  arrange(padj) |> 
  write_csv(fs::path(paper_tables, "table_S7.csv"))


# 