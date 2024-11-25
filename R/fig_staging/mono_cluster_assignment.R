mono_cluster_assignment <- bb_cellmeta(cds_main) |> 
  filter(partition_assignment %in% c("Mono")) |> 
  count(seurat_l2_leiden_consensus, seurat_celltype_l2) |>
  group_by(seurat_l2_leiden_consensus) |> 
  mutate(count = sum(n)) |> 
  mutate(percent = n/count * 100) |> 
  # filter(percent > 1) |>
  filter(n >20) |> 
  ggplot(aes(x = seurat_l2_leiden_consensus, fill = seurat_celltype_l2, y = n)) +
  geom_bar(color = "black",
           stat = "identity",
           position = "fill") + 
  labs(fill = "Reference celltype", x = NULL, y = "Proportion of cluster")
