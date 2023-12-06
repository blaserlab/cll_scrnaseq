bb_var_umap(cds_main, "seurat_l2_leiden_consensus")

bb_cellmeta(cds_main) |> glimpse()

library(CellChat)
cellchat_res_full <- bb_cellchat(cds = cds_main, group_var = "seurat_l2_leiden_consensus", species = "human", min_cells = 50, n_cores = 24)
