# global umap with partition labels ----------------------------------
umap_partitions <- bb_var_umap(
  obj = cds_main,
  var = "partition_assignment",
  overwrite_labels = T,
  group_label_size = 4,
  foreground_alpha = 0.1
)

umap_bcell_leiden <- 
  bb_var_umap(cds_main[, colData(cds_main)$partition_assignment %in% c("B")],
              var = "leiden", 
              overwrite_labels = T,
              foreground_alpha = 0.1)

umap_binned_leiden <-
  bb_var_umap(cds_main[, colData(cds_main)$partition_assignment %in% c("B")],
              var = "leiden_assignment_binned_renamed", 
              overwrite_labels = T,
              foreground_alpha = 0.1)

# go term bubbles--------------------------------------------------

mod4_bubble <- bb_goscatter(simMatrix = module_summary_0.9$`Module 4`$simMatrix, 
             reducedTerms = module_summary_0.9$`Module 4`$reducedTerms,size = "score",
             addLabel = T,
             labelSize = 4) 



# subcluster go term bubbles ----------------------------

bb_goscatter(simMatrix = mrd1_summary_0.8$simMatrix, reducedTerms = mrd1_summary_0.8$reducedTerms)# maybe pro-inflammatory
bb_goscatter(simMatrix = mrd2_summary_0.8$simMatrix, reducedTerms = mrd2_summary_0.8$reducedTerms)# mrd2 looks like stressed
bb_goscatter(simMatrix = btk_summary_0.8$simMatrix, reducedTerms = btk_summary_0.8$reducedTerms)# btk cll-like
