umap_numbered_leiden <-
  bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> 
                           filter(partition_assignment == ("B"))), 
              var = "leiden", 
              overwrite_labels = TRUE)
umap_numbered_leiden
