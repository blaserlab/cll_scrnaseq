# reference annotation barchart -------------------------------------------
reference_anno_barchart <- bb_cellmeta(cds_main) |> 
  count(partition_assignment, seurat_celltype_l2) |>
  group_by(partition_assignment) |> 
  mutate(partition_count = sum(n)) |> 
  mutate(partition_percent = n/partition_count * 100) |> 
  filter(partition_percent > 5) |>
  filter(n >20) |> 
  ggplot(aes(x = partition_assignment, fill = seurat_celltype_l2, y = n)) +
  geom_bar(color = "black",
           stat = "identity",
           position = "fill") + 
  labs(fill = "Reference celltype", x = NULL, y = "Proportion of cluster")