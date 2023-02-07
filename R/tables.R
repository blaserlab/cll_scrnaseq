# gene modules-------------------
rowData(cds_main) %>%
  as_tibble() %>%
  select(module_labeled, gene_short_name, id) %>%
  write_csv(str_glue("{network_tables}/gene_modules.csv"))

# go term enrichments------------------------------

map2_dfr(
  .x = module_enrichment,
  .y = names(module_enrichment),
  .f = function(x,y) {
    go_tbl <- x[[3]] %>%
      mutate(module = y) %>%
      relocate(module)
    
  }
) %>% write_csv(str_glue("{network_tables}/module_go_term_enrichment.csv"))

bind_rows(btk_enrich[[3]] |> mutate(subcluster = "CLL-like") |> relocate(subcluster), 
          mrd1_enrich[[3]] |> mutate(subcluster = "inflammatory") |> relocate(subcluster), 
          mrd2_enrich[[3]] |> mutate(subcluster = "stressed") |> relocate(subcluster)) |> 
  write_csv(fs::path(network_tables, "subcluster_goterm_enrichment.csv"))



# pseudobulk dge--------------------------

# B cells in BTK patients across timepoints
walk2(
  .x = pseudobulk_bcell_btk_timepoints,
  .y = names(pseudobulk_bcell_btk_timepoints),
  .f = function(x, y) {
    res <- x[[2]] %>% arrange(padj)
    write_csv(x = res,
              file = str_glue("{network_tables}/btk_pts_{y}.csv"))
    header <- x[[1]]
    write_lines(header, file = str_glue("{network_tables}/btk_pts_{y}_header.txt"))
    }
  )

# B cells in BTK vs MRD at each timepoint -------------------------------
walk2(
  .x = pseudobulk_bcell_pt_type_timepoints,
  .y = names(pseudobulk_bcell_pt_type_timepoints),
  .f = function(x,y) {
    res <- x[[2]] %>% arrange(padj)
    write_csv(x = res,
              file = str_glue("{network_tables}/MRD_v_BTK_{y}.csv"))
    header <- x[[1]]
    write_lines(header, file = str_glue("{network_tables}/MRD_v_BTK_{y}_header.txt"))
    
  }
)

# B cell subclusters (pseudobulk)--------------------------------------

pseudobulk_MRD1_BTK[[1]] %>% write_lines(str_glue("{network_tables}/MRD1_v_BTK_cluster_header.txt"))
pseudobulk_MRD1_BTK[[2]] %>% write_csv(str_glue("{network_tables}/MRD1_v_BTK_cluster.csv"))


pseudobulk_MRD2_BTK[[1]] %>% write_lines(str_glue("{network_tables}/MRD2_v_BTK_cluster_header.txt"))
pseudobulk_MRD2_BTK[[2]] %>% write_csv(str_glue("{network_tables}/MRD2_v_BTK_cluster.csv"))

# B cell subcluster top markers (not pseudobulk)--------------------

cds_main_bcell_subpop_top_markers %>%
  mutate(cell_group = recode(cell_group, "MRD1" = "inflammatory", "MRD2" = "stressed", "BTK" = "CLL-like")) |> 
  write_csv(str_glue("{network_tables}/bcell_subpop_top_markers.csv"))

# T cell phenotype:  pseudobulk analsysis_--------------------------------------

tcell_subpop_pseudobulk[[1]] %>% write_lines(str_glue("{network_tables}/tcell_subpop_pseudobulk_header.txt"))
tcell_subpop_pseudobulk_tbl %>% write_csv(str_glue("{network_tables}/tcell_subpop_pseudobulk.csv"))
tcell_subpop_pseudobulk_tbl %>% filter(padj<0.05) %>% write_csv(str_glue("{network_tables}/tcell_subpop_pseudobulk_significant.csv"))

# T cell regression

tcell_subpop_regression %>% write_csv(str_glue("{network_tables}/tcell_subpop_regression.csv"))
