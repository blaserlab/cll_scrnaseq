cell_representation_data <- map_dfr(.x = c("Timepoint 1", "Timepoint 2", "Timepoint 3"), 
    .f = \(x, dat = cds_main) {
      bb_cluster_representation2(
        filter_cds(
          dat,
          cells = bb_cellmeta(dat) |> 
            filter(timepoint_merged_2 == x) |>  
            filter(partition_assignment %in% c("B", "T", "NK", "Mono"))
        ),
        sample_var = "specimen",
        cluster_var = "partition_assignment",
        comparison_var = "patient_type2",
        comparison_levels = c("resistant", "sensitive"),
        return_val = "data"
      ) |> 
        mutate(timepoint = x)
      
      }
    )

cell_representation_barchart <- ggplot(cell_representation_data,
       aes(
         x = fct_relevel(cluster, c("B", "T", "NK", "Mono")),
         y = logFC,
         fill = enriched
       )) +
  geom_bar(color = "black", stat = "identity") +
  geom_text(mapping = aes(x = cluster, y = texty, label = sig)) +
  scale_fill_manual(values = experimental_group_palette) +
  facet_wrap( ~ timepoint) +
  labs(y = "Log<sub>2</sub> Fold Change", x = "Cluster") +
  panel_border() +
  theme(axis.title.y = ggtext::element_markdown()) +
  theme(strip.background = element_blank()) +
  theme(legend.position = "top", legend.justification = "center")
