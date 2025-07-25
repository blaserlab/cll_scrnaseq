pa_genebub_markers <- cds_main_top_markers |> 
  filter(cluster_method == "partition") |> 
  mutate(partition = str_remove(cell_group, "partition ")) |> 
  relocate(partition) |> 
  select(-c(cell_group, cluster_method)) |> 
  left_join(bb_cellmeta(cds_main) |> group_by(partition, partition_assignment) |> summarise()) |> 
  relocate(partition_assignment) |> 
  group_by(partition_assignment) |> 
  filter(!str_detect(gene_short_name, "IG.*")) |> 
  slice_max(pseudo_R2, n = 4) |> pull(gene_short_name)

pa_genebub_dat <-
  bb_genebubbles(
    cds_main,
    genes = pa_genebub_markers,
    cell_grouping = "partition_assignment",
    scale_expr = TRUE,
    return_value = "data"
  )


partition_assignment_genebubbles <-
  ggplot(
    pa_genebub_dat,
    mapping = aes(
      x = partition_assignment,
      y = gene_short_name,
      fill = expression,
      size = proportion
    )
  ) +
  geom_point(pch = 21, color = "black") +
  scale_fill_viridis_c() +
  scale_size_area() +
  theme_minimal_grid(font_size = 10) + 
  theme(axis.text.y = element_text(face = "italic")) +
  labs(x = NULL, y = NULL, fill = "Expression", size = "Proportion\nExpressing")

