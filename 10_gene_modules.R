source("00_packages_functions.R")

# identify gene modules ####--------------------------------------------------------------
# warning!! Slow!! Destructive!!
pr_graph_test_res <-
  graph_test(
    cds = cds_final,
    neighbor_graph = "knn",
    cores = 16,
    verbose = TRUE
  )

pr_deg_ids <- row.names(subset(pr_graph_test_res, q_value < 0.05))

gene_module_df <-
  find_gene_modules(cds_final[pr_deg_ids,], cores = 39)

gene_module_df_anno <-
  left_join(gene_module_df, as_tibble(rowData(cds_final))) %>%
  write_csv("data_out/gene_module_df_anno.csv")

# print out the individual gene modules
dir.create("data_out/individual_modules")

lapply(
  X = gene_module_df_anno %>% pull(module) %>% unique(),
  FUN = function(x) {
    gene_module_df_anno %>% filter(module == x) %>% write_csv(paste0("data_out/individual_modules/module_", x, ".csv"))
  }
)

# make the gene module figures ####--------------------------------------------------------
cell_group_df <-
  tibble::tibble(cell = row.names(colData(cds_pmm_final)),
                 cell_group = colData(cds_pmm_final)$cluster_assignment) %>%
  filter(cell_group %in% major_clusters)
agg_mat_pmm <-
  aggregate_gene_expression(cds_pmm_final, gene_module_df_pmm, cell_group_df_pmm)

row.names(agg_mat_pmm) <- stringr::str_c("Module ", row.names(agg_mat_pmm))

colnames(agg_mat_pmm) <- stringr::str_c(colnames(agg_mat_pmm))

# expand the aggregation matrix to list the modules by genes, with the score from each gene equal to the module score
# do this just to make annotating the plot with selected gene names possible

pmm_highlights <- gene_module_df_anno_pmm %>% 
  filter(
    gene_short_name %in% c(
      "spi1b",
      "cebpa",
      "lyz",
      "tal1",
      "gfi1b",
      "gata1a",
      "gata2a",
      "il21",
      "ccr9a",
      "cldnh",
      "clcnk",
      "hoxb8a",
      "slc26a1",
      "slc5a12",
      "slc13a3",
      "crestin",
      "hps1",
      "hps4",
      "kif17",
      "dnaaf1",
      "tube1",
      "krt17",
      "krt97",
      "krt96",
      "ptges",
      "pdgfra",
      "lfng",
      "kdrl",
      "egfl7",
      "cdh5"
      
      
       
    )
  ) %>%
  select(module,gene_short_name) %>% 
  mutate(module_full = paste0("Module ",module))

gene_module_df_anno_pmm


# module data:

# mod2
# http://cbl-gorilla.cs.technion.ac.il/GOrilla/gji1al4n/GOResults.html
# muscle

# mod13
# http://cbl-gorilla.cs.technion.ac.il/GOrilla/ziu2nch3/GOResults.html
# meox.  Chest cell 

# mod4
# http://cbl-gorilla.cs.technion.ac.il/GOrilla/2esv0sr1/GOResults.html
# liver genes

# module 20
# http://cbl-gorilla.cs.technion.ac.il/GOrilla/alaipd80/GOResults.html
# translation

# module 12
# http://cbl-gorilla.cs.technion.ac.il/GOrilla/vi790q2o/GOResults.html
# ectoderm????

# module 22
# http://cbl-gorilla.cs.technion.ac.il/GOrilla/96lraskc/GOResults.html
# antivirus

# module 9
# http://cbl-gorilla.cs.technion.ac.il/GOrilla/tif3mdb9/GOResults.html
# ciliated

# module 6 
# http://cbl-gorilla.cs.technion.ac.il/GOrilla/yd9o0fp7/GOResults.html
# melanosome?

# module 1
# http://cbl-gorilla.cs.technion.ac.il/GOrilla/keqok7m3/GOResults.html
# leukocytes



# module 3
# http://cbl-gorilla.cs.technion.ac.il/GOrilla/zz66c26z/GOResults.html
# translation/heme

# GO analysis of gene module genes
modules_of_interest <- c("1","3","8","17","13")
module_titles <- paste0("Module ",modules_of_interest)
module_pvals <- c(0.01,0.01,0.01,0.01,0.01)

pmm_module_summaries <- pmap(
  .l = list(module_of_interest = modules_of_interest,
            module_title = module_titles,
            module_pval = module_pvals),
  .f = possibly(function(module_of_interest, module_title, module_pval) {
    module_genes_named <- gene_module_df_anno_pmm %>%
      select(module, gene_short_name) %>%
      mutate(selected = ifelse(module == module_of_interest, 1, 0)) %>%
      pull(selected)
    names(module_genes_named) <- gene_module_df_anno_pmm %>%
      select(module, gene_short_name) %>%
      mutate(selected = ifelse(module == module_of_interest, 1, 0)) %>%
      pull(gene_short_name)
    module_genes_named
    
    sampleGOdata <- new(
      "topGOdata",
      description = "Simple session",
      ontology = "BP",
      allGenes = module_genes_named,
      geneSel = selector,
      nodeSize = 10,
      annot = annFUN.org,
      mapping = "org.Dr.eg.db",
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
    
    mod_mat <-
      simplifyEnrichment::GO_similarity(
        resultFisher_tbl %>% filter(pval < module_pval) %>% pull(goterm),
        ont = "BP",
        db = "org.Dr.eg.db",
      )
    
    plot <-
      simplifyEnrichment::simplifyGO(mod_mat, column_title = module_title)
    
    return(list(sampleGOdata, resultFisher, res_table, plot))
  }, 
  otherwise = NULL,
  quiet = F
))
names(pmm_module_summaries) <- module_titles

pmm_rrvgo <- map(
  .x = pmm_module_summaries,
  .f = function(x) {
    simMatrix <-
      calculateSimMatrix(x = x[[3]]$GO.ID,
                         ont = "BP",
                         orgdb = "org.Dr.eg.db")
    scores <- setNames(-log10(ifelse(is.na(
      as.numeric(x[[3]]$classicFisher)),
      1e-30,
      as.numeric(x[[3]]$classicFisher)
    )),
    x[[3]]$GO.ID)
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    scores,
                                    threshold = 0.9,
                                    orgdb = "org.Dr.eg.db")
    returnlist <- list(simMatrix, scores, reducedTerms)
    names(returnlist) <- c("simMatrix", "scores", "reducedTerms")
    return(returnlist)
  }
)
names(pmm_rrvgo) <- module_titles

rrvgo_scatter <-
  function (simMatrix,
            reducedTerms,
            size = "score",
            addLabel = TRUE,
            labelSize = 4) {
    x <- cmdscale(as.matrix(as.dist(1 - simMatrix)), eig = TRUE,
                  k = 2)
    df <-
      cbind(as.data.frame(x$points), reducedTerms[match(rownames(x$points),
                                                        reducedTerms$go), c("term", "parent", "parentTerm", "size", "score")])
    p <-
      ggplot(df, aes(x = V1, y = V2, color = parentTerm)) +
      geom_point(aes(size = !!sym(size)), alpha = 0.5) +
      scale_color_discrete(guide = FALSE) +
      scale_size_continuous(guide = FALSE, range = c(0, 10)) +
      geom_text_repel(
        aes(label = parentTerm),
        segment.size = 0.25,
        data = subset(df, parent == rownames(df)),
        box.padding = grid::unit(0.5, "lines"),
        size = labelSize,
        color = "black",
        max.overlaps = 100,
        force = 2,
        seed = 1234,
        segment.curvature = -0.1,
        segment.square = TRUE,
        segment.inflect = TRUE,
        min.segment.length = 0
      ) +
      labs(x = "PCoA 1", y = "PCoA 2")
   
   return(p)
  }

pmm_module_bubble_plots <-
  pmap(
    .l = list(
      simmatrices = subListExtract(pmm_rrvgo, "simMatrix"),
      reducedterms = subListExtract(pmm_rrvgo, "reducedTerms"),
      listnames = names(pmm_rrvgo)
    ),
    .f = function(simmatrices, reducedterms, listnames) {
      p <- rrvgo_scatter(
        simMatrix = simmatrices,
        reducedTerms = reducedterms,
        size = "score",
        labelSize = 2
      ) +
        labs(title = listnames) +
        theme(plot.title = element_text(hjust = 0.5))
      return(p)
    }
  )
