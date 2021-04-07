source("00_packages_functions.R")

# determine expression of all zf genes aggregated by cluster
# including only the heme clusters
heme_cluster_expression <- as_tibble(
  aggregate_gene_expression(
    cds = cds_pmm_final[, colData(cds_pmm_final)$cluster_binary_type == "heme"],
    cell_group_df = cell_group_df_pmm %>% filter(
      cell_group %in% c(
        "Neutrophil 1",
        "Neutrophil 2",
        "Pro-neutrophil",
        "Progenitor 1",
        "Progenitor 2",
        "Proliferative",
        "Erythroid"
      )
    )
  ) %>% as.matrix(),
  rownames = "gene"
) %>% pivot_longer(-gene, names_to = "cluster", values_to = "expression")
  

# generate a table of zf/mouse orthologs
# filter only to include 1:1 mappings
zf_id_mouse_sym <- zf_mouse_orthos %>% # zf_mouse_orthos defined in 30_muench_data.R
  left_join(as_tibble(rowData(cds_pmm_final)), by = c("zfin_symbol" = "gene_short_name")) %>% # joining by zf gene name
  rename(zf_id = id) %>% # rename ensdarg column from id to zf_id
  group_by(zf_id,mouse_symbol) %>%
  summarise() %>% # removes all rows where the tuple (zf_id,mouse_symbol) is duplicated
  group_by(zf_id) %>%
  mutate(duplicate_flag = n() > 1) %>%
  filter(!duplicate_flag) %>% # removes all rows where zf_id is not unique  
  group_by(mouse_symbol) %>%
  mutate(duplicate_flag = n() > 1) %>%
  filter(!duplicate_flag) %>% # removes all rows where mouse_symbol is not unique
  select(gene = zf_id, mouse_symbol)


# makes a table of heme cluster expression but labeled with orthologous mouse genes
zf_mouse_expression <- left_join(
  heme_cluster_expression,zf_id_mouse_sym
  ) %>% 
  select(mouse_symbol,cluster, expression) %>%
  filter(!is.na(mouse_symbol)) 

# fgsea
# this is the gsea analysis on the muench and tusi gene sets which are from the mouse
zf_mouse_gsea_res <- map(.x = zf_mouse_expression %>% pull(cluster) %>% unique(),
    .f = function(x, pathway_list = c(muench_gene_sets,tusi_gene_sets), data = zf_mouse_expression) {
      data <- data %>% filter(cluster == x) %>% arrange(expression)
      data_vec <- data %>% pull(expression)
      names(data_vec) <- data %>% pull(mouse_symbol)
      res <- fgsea(pathways = pathway_list, stats = data_vec)
      PathwaysRanked <- res[order(pval), pathway]
      res <- as_tibble(res) %>% mutate(cluster = x)
      plot <- plotGseaTable(pathways = pathway_list[PathwaysRanked], stats = data_vec, fgseaRes = res, gseaParam=0.5,render = F)
      return(list(res,plot))
    }) %>% set_names(zf_mouse_expression %>% pull(cluster) %>% unique())

# inspect the gsea plots and data
# grid.draw(zf_mouse_gsea_res$Erythroid[[2]])
# grid.draw(zf_mouse_gsea_res$`Myeloid 1`[[2]])
# grid.draw(zf_mouse_gsea_res$`Myeloid 2`[[2]])
# grid.draw(zf_mouse_gsea_res$`Myeloid 3`[[2]])
# grid.draw(zf_mouse_gsea_res$`Progenitor 1`[[2]])
# grid.draw(zf_mouse_gsea_res$`Progenitor 2`[[2]])
# grid.draw(zf_mouse_gsea_res$`Prolif. Myeloid`[[2]])

# this is the gsea analysis on the tang gene sets which are from zf
zf_gsea_res <- map(.x = heme_cluster_expression %>% pull(cluster) %>% unique(),
                   .f = function(x, pathway_list = zf_gene_sets, data = heme_cluster_expression) {
                     data <- data %>% filter(cluster == x) %>% arrange(expression)
                     data_vec <- data %>% pull(expression)
                     names(data_vec) <- data %>% pull(gene)
                     res <- fgsea(pathways = pathway_list, stats = data_vec)
                     PathwaysRanked <- res[order(pval),pathway]
                     res <- as_tibble(res) %>% mutate(cluster = x)
                     plot <- plotGseaTable(pathways = pathway_list[PathwaysRanked], stats = data_vec, fgseaRes = res, gseaParam = 0.5, render = F)
                     return(list(res,plot))
                   }) %>% set_names(heme_cluster_expression %>% pull(cluster) %>% unique())

# inspect the gsea plots and data
# grid.draw(zf_gsea_res$Erythroid[[2]])
# grid.draw(zf_gsea_res$`Myeloid 1`[[2]])
# grid.draw(zf_gsea_res$`Myeloid 2`[[2]])
# grid.draw(zf_gsea_res$`Myeloid 3`[[2]])
# grid.draw(zf_gsea_res$`Progenitor 1`[[2]])
# grid.draw(zf_gsea_res$`Progenitor 2`[[2]])
# grid.draw(zf_gsea_res$`Prolif. Myeloid`[[2]])


# get all of the tabular results
all_gsea_data <- bind_rows(map_dfr(zf_mouse_gsea_res,1) %>% 
			   mutate(comparison = "zf_to_mouse"), 
		   map_dfr(zf_gsea_res,1) %>% 
			   mutate(comparison = "zf_to_zf"))

# set up the heatmap matrix in the form of rows = cluster, cols = pathway, values = -log10(padj)
gsea_data_mtx <- all_gsea_data %>%
  mutate(neg_log_padj = -1 * log10(padj)) %>%
  select(cluster,pathway, neg_log_padj) %>%
  pivot_wider(names_from = pathway, values_from = neg_log_padj) %>%
  tbl_to_matrix()


