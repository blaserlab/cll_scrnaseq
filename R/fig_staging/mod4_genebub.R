# data exploration
# counts <- bb_aggregate(
#   filter_cds(
#     cds_main,
#     cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "B"),
#     genes = bb_rowmeta(cds_main) |> filter(module == "4")
#   ),
#   cell_group_df = bb_cellmeta(cds_main) |> 
#     filter(partition_assignment == "B") |>  
#     # mutate(type_timepoint = paste0(patient_type3, "_", timepoint_merged_1)) |> 
#     # select(cell_id, type_timepoint))
#     select(cell_id, patient_type3))
# 
# mod4_diff_genes <- counts |> as.matrix() |> as_tibble(rownames = "feature_id") |> 
#   left_join(bb_rowmeta(cds_main) |> 
#               select(feature_id, gene_short_name)) |> 
#   relocate(gene_short_name) |> 
#   mutate(diff = IBS-IBR) |> 
#   arrange(desc(diff)) |> 
#   filter(diff > 0.003) |> 
#   pull(gene_short_name)

# what go terms are associated with these genes?
## Install packages once if needed:
## BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db", "GO.db"))
# 
# library(AnnotationDbi)
# library(org.Hs.eg.db)
# library(GO.db)
# 
# get_bp_terms_for_genes <- function(genes, include_ancestors = TRUE) {
#   genes <- unique(genes)
#   
#   # Choose columns: direct BP only vs including ancestor terms
#   col_go  <- if (include_ancestors) "GOALL"      else "GO"
#   col_ont <- if (include_ancestors) "ONTOLOGYALL" else "ONTOLOGY"
#   
#   # Map SYMBOL -> GO IDs with ontology
#   anno <- AnnotationDbi::select(
#     org.Hs.eg.db,
#     keys    = genes,
#     columns = c(col_go, col_ont),
#     keytype = "SYMBOL"
#   )
#   
#   # Keep only Biological Process terms
#   anno_bp <- anno[!is.na(anno[[col_go]]) & anno[[col_ont]] == "BP", ]
#   
#   if (nrow(anno_bp) == 0L) {
#     # If nothing found, return empty character vectors for each gene
#     out <- lapply(genes, function(x) character(0))
#     names(out) <- genes
#     return(out)
#   }
#   
#   # Get GO term names from GO.db
#   go_info <- AnnotationDbi::select(
#     GO.db,
#     keys    = unique(anno_bp[[col_go]]),
#     columns = "TERM",
#     keytype = "GOID"
#   )
#   
#   # Merge term names back onto annotation
#   merged <- merge(
#     anno_bp, go_info,
#     by.x = col_go,
#     by.y = "GOID",
#     all.x = TRUE
#   )
#   
#   # Build a named list: gene -> vector of unique BP term names
#   out <- lapply(genes, function(g) {
#     terms <- merged$TERM[merged$SYMBOL == g]
#     unique(terms[!is.na(terms)])
#   })
#   names(out) <- genes
#   
#   out
# }
# bp_list <- get_bp_terms_for_genes(mod4_genes_plus)
# 
# map(.x = bp_list, 
#     .f = \(x) {
#       x[str_detect(x, "T cell")]
#     })
# bp_list1 <- get_bp_terms_for_genes(mod4_diff_genes)
# map(.x = bp_list1,
#     .f = \(x) {
#       x[str_detect(x, "T cell")]
#     })
# 
# map(.x = bp_list1,
#     .f = \(x) {
#       x[str_detect(x, "immune system process")]
#     })

mod4_genes <- c(
  "LTB", 
  # "ITGB1", 
  "DUSP2", 
  "LITAF", 
  "NCR3", 
  "CD5", 
  # "PRDX2", 
  "IL6ST", 
  "NCK1", 
  # "RORA", 
  # "BATF", 
  "PIK3R1", 
  # "TOX", 
  # "TCF7", 
  "CD38", 
  "ZFP36L2")

other_genes <- c(
  "NKG7", 
  "GNLY", 
  "PRF1", 
  "KLRD1", 
  "KLRC1",   
  "CD8A", 
  "CD8B", 
  "CD3E", 
  "CD14", 
  "CD247")

mod4_genes_plus <- bind_rows(
tibble(gene_short_name = mod4_genes, title = "Module 4"),
tibble(gene_short_name = other_genes, title = "Other T/NK Genes")
)

# bp_list2 <- get_bp_terms_for_genes(mod4_genes_plus$gene_short_name)
# map(.x = bp_list2,
#     .f = \(x) {
#       x[str_detect(x, "T cell")]
#     })
# 
# map(.x = bp_list2,
#     .f = \(x) {
#       x[str_detect(x, "immune system process")]
#     })
# bp_list2$CD38

mod4_genebubdat <- bb_genebubbles(obj = filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> 
                                                     filter(partition_assignment == "B")), 
                                  cell_grouping = c("patient_type3", "timepoint_merged_1"), 
                                  genes = mod4_genes_plus$gene_short_name, 
                                  return_value = "data") |> 
  left_join(mod4_genes_plus, by = join_by(gene_short_name))

mod4_genebub1 <- ggplot(mod4_genebubdat |> filter(title == "Module 4"), aes(x = patient_type3, y = gene_short_name, fill = expression, size = proportion)) +
  geom_point(pch = 21, color = "black") +
  scale_fill_viridis_c(limits = c(-1.55, 2.1)) +
  scale_size_area(limits = c(0, 0.93)) +
  theme_minimal_grid() +
  theme(axis.text.y = element_text(face = "italic")) + 
  facet_wrap(~timepoint_merged_1) + 
  labs(x = NULL, y = NULL, fill = "Expression", size = "Proportion", subtitle = "<i>B cell-Expressed Mod.4</i>") + 
  theme(plot.subtitle = ggtext::element_markdown())
mod4_genebub2 <-ggplot(mod4_genebubdat |> filter(title == "Other T/NK Genes"), aes(x = patient_type3, y = gene_short_name, fill = expression, size = proportion)) +
  geom_point(pch = 21, color = "black") +
  scale_fill_viridis_c(limits = c(-1.55, 2.1)) +
  scale_size_area(limits = c(0, 0.93)) +
  theme_minimal_grid() +
  theme(axis.text.y = element_text(face = "italic")) + 
  facet_wrap(~timepoint_merged_1) + 
  labs(x = NULL, y = NULL, fill = "Expression", size = "Proportion", subtitle = "<i>B cell-Not Expressed</i>") + 
  theme(plot.subtitle = ggtext::element_markdown())
mod4_genebub <- mod4_genebub1 + mod4_genebub2 + plot_layout(guides = "collect")

# mod4_genes_plus |> 
#   mutate(in_mod4 = case_when(gene_short_name %in% (bb_rowmeta(cds_main) |> filter(module == 4) |> pull(gene_short_name)) ~ TRUE, .default = FALSE)) |> View()
