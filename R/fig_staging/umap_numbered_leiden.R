umap_numbered_leiden <-
  bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> 
                           filter(partition_assignment == ("B"))), 
              var = "leiden", 
              foreground_alpha = 0.02,
              rasterize = TRUE,
              overwrite_labels = TRUE)
umap_numbered_leiden
