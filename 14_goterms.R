source("00_packages_functions.R")

# make a table of gene names in one column and grouping value in another column 
monocyte_go_query <- bind_rows(list( 
  pseudobulk_res_mono[[1]][[2]] %>% # timepoint2 vs timepoint 1 
  filter(padj<0.05) %>%
  filter(log2FoldChange<0) %>%
  select(gene_short_name) %>%
  mutate(gene_grouping = "down in timepoint 2 vs 1"),
  pseudobulk_res_mono[[1]][[2]] %>% # timepoint2 vs timepoint 1 
  filter(padj<0.05) %>%
  filter(log2FoldChange>0) %>%
  select(gene_short_name) %>%
  mutate(gene_grouping = "up in timepoint 2 vs 1"),
  pseudobulk_res_mono[[2]][[2]] %>% # timepoint3 vs timepoint 1 
  filter(padj<0.05) %>%
  filter(log2FoldChange<0) %>%
  select(gene_short_name) %>%
  mutate(gene_grouping = "down in timepoint 3 vs 1"),
  pseudobulk_res_mono[[2]][[2]] %>% # timepoint3 vs timepoint 1 
  filter(padj<0.05) %>%
  filter(log2FoldChange>0) %>%
  select(gene_short_name) %>%
  mutate(gene_grouping = "up in timepoint 3 vs 1")
  ))
# monocyte_go_query %>% View()
  
  
mono_go_summaries <- pmap(
  .l = list(gene_group = unique(monocyte_go_query$gene_grouping),
            group_pval = c(0.01, 0.01, 0.01, 0.01)),
  .f = possibly(function(gene_group, group_pval, query = monocyte_go_query, reference = as_tibble(rowData(cds_final))) {
    genes <- query %>% filter(gene_grouping == gene_group) %>% pull(gene_short_name)
    
    genes_named <- reference %>%
      mutate(selected = ifelse(gene_short_name %in% genes, 1, 0)) %>%
      pull(selected)
    names(genes_named) <- reference %>%
      pull(gene_short_name)

    sampleGOdata <- new(
      "topGOdata",
      description = "Simple session",
      ontology = "BP",
      allGenes = genes_named,
      geneSel = selector,
      nodeSize = 10,
      annot = annFUN.org,
      mapping = "org.Hs.eg.db",
      ID = "symbol"
    )

    resultFisher <-
      runTest(sampleGOdata, algorithm = "classic", statistic = "fisher")

    res_table <- GenTable(
      sampleGOdata,
      classicFisher = resultFisher,
      orderBy = "classicFisher",
      ranksOf = "classicFisher",
      topNodes = 100
    ) %>%
      as_tibble(rownames = "Rank")

    resultFisher_tbl <-
      tibble(goterm = names(resultFisher@score),
             pval = resultFisher@score)

    # mod_mat <-
    #   simplifyEnrichment::GO_similarity(
    #     resultFisher_tbl %>% filter(pval < group_pval) %>% pull(goterm),
    #     ont = "BP",
    #     db = "org.Hs.eg.db",
    #   )
    # 
    # plot <-
    #   simplifyEnrichment::simplifyGO(mod_mat, column_title = gene_group)

    return(list(sampleGOdata, resultFisher, res_table))
  }, 
  otherwise = NULL,
  quiet = F
)) %>% set_names(unique(monocyte_go_query$gene_grouping))

mono_go_summaries




mono_rrvgo_thresh_0.8 <- map(
  .x = mono_go_summaries,
  .f = function(x) {
    summarize_go(x = x, reduce_threshold = 0.8)
  }
) %>% set_names(names(mono_go_summaries))

mono_rrvgo_thresh_0.9 <- map(
  .x = mono_go_summaries,
  .f = function(x) {
    summarize_go(x = x, reduce_threshold = 0.9)
  }
) %>% set_names(names(mono_go_summaries))

mono_module_bubble_plots_0.8 <-
  pmap(
    .l = list(
      simmatrices = subListExtract(mono_rrvgo_thresh_0.8, "simMatrix"),
      reducedterms = subListExtract(mono_rrvgo_thresh_0.8, "reducedTerms"),
      listnames = names(mono_rrvgo_thresh_0.8)
    ),
    .f = function(simmatrices, reducedterms, listnames) {
      p <- rrvgo_scatter(
        simMatrix = simmatrices,
        reducedTerms = reducedterms,
        size = "score",
        labelSize = 3
      ) +
        labs(title = listnames) +
        theme(plot.title = element_text(hjust = 0.5))
      return(p)
    }
  ) %>% set_names(names(mono_rrvgo_thresh_0.8))

mono_module_bubble_plots_0.9 <-
  pmap(
    .l = list(
      simmatrices = subListExtract(mono_rrvgo_thresh_0.9, "simMatrix"),
      reducedterms = subListExtract(mono_rrvgo_thresh_0.9, "reducedTerms"),
      listnames = names(mono_rrvgo_thresh_0.9)
    ),
    .f = function(simmatrices, reducedterms, listnames) {
      p <- rrvgo_scatter(
        simMatrix = simmatrices,
        reducedTerms = reducedterms,
        size = "score",
        labelSize = 3
      ) +
        labs(title = listnames) +
        theme(plot.title = element_text(hjust = 0.5))
      return(p)
    }
  ) %>% set_names(names(mono_rrvgo_thresh_0.9))


# print out the module summaries

dir.create("plots_out/go_term_bubbles")
walk2(
  .x = mono_module_bubble_plots_0.9,
  .y = c("down_in_timept2", "up_in_timept2", "down_in_timept3", "up_in_timept3"),
  .f = function(x,y) {
    save_plot(x, filename = paste0("plots_out/go_term_bubbles/",y,".pdf"), base_width = 10, base_height = 7.5)
  }
)


