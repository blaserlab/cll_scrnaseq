# bb_var_umap(cds_main, "seurat_l2_leiden_consensus")
specimens <- bb_cellmeta(cds_main) |> 
  pull(specimen) |> 
  unique()

cellchat_obj_list <- map(.x = specimens,
                         .f = \(x, dat = cds_main) {
                           dat <- filter_cds(dat, cells = bb_cellmeta(dat) |> filter(specimen == x))
                           bb_cellchat(cds = dat, 
                                 group_var = "seurat_l2_leiden_consensus", 
                                 n_cores = 4, 
                                 species = "human")
                         })