



cell_representation_data <- map_dfr(.x = c("Timepoint 1", "Timepoint 2", "Timepoint 3"), 
    .f = \(x, dat = cds_main) {
      bb_cluster_representation2(
        filter_cds(
          dat,
          cells = bb_cellmeta(dat) |> 
            filter(timepoint_merged_2 == x) |>  
            filter(partition_assignment %in% c("B"))
        ),
        sample_var = "specimen",
        cluster_var = "leiden",
        comparison_var = "patient_type2",
        comparison_levels = c("resistant", "sensitive"),
        return_val = "data"
      ) |> 
        mutate(timepoint = x)
      
      }
    )

cell_representation_data |> filter(enriched == "sensitive")

bb_var_umap(cds_main, "leiden", overwrite_labels = TRUE)

# population fold change plot-----------------------------------------------------------
normalized_leiden_counts <- 
  bb_cellmeta(cds_main) |> 
  filter(partition_assignment == "B") |> 
  mutate(leiden_assignment_binned_renamed_1 = recode(leiden_assignment_binned_renamed, "stressed" = "CLL-like")) |> 
  count(patient, leiden_assignment_binned_renamed_1, specimen, timepoint_merged_1, patient_type2)  |> 
  left_join(bb_cellmeta(cds_main) |> 
              filter(partition_assignment == "B") |> 
              group_by(specimen) |> 
              summarise(specimen_total = n())) |> 
  mutate(overall_total = nrow(bb_cellmeta(cds_main) |> filter(partition_assignment =="B"))) |> 
  mutate(normalized_count = n*overall_total/specimen_total/2) |> 
  select(leiden_assignment_binned_renamed_1, specimen, timepoint_merged_1, patient_type2, normalized_count)

cluster_proportion_ratio_plot <- normalized_leiden_counts %>%
  pivot_wider(names_from = leiden_assignment_binned_renamed_1, values_from = normalized_count, values_fill = 1) %>%
  mutate(btk_to_other_ratio = (`CLL-like`)/(inflammatory)) %>%
  mutate(log2_ratio = log2(btk_to_other_ratio)) %>%
  ggplot(mapping = aes(x = patient_type2, y = log2_ratio, color = patient_type2, fill = patient_type2)) +
  geom_jitter(shape = jitter_shape, size = jitter_size, stroke = jitter_stroke) +
  facet_wrap(facets = vars(timepoint_merged_1)) +
  scale_fill_manual(values = alpha(colour = experimental_group_palette, alpha = jitter_alpha_fill)) +
  scale_color_manual(values = alpha(colour = experimental_group_palette, alpha = jitter_alpha_color)) +
  theme(strip.background = element_blank()) +
  theme(panel.background = element_rect(color = "grey80")) +
  theme(legend.position = "none") +
  stat_summary(
    fun.data = data_summary_mean_se,
    color = summarybox_color,
    size = summarybox_size,
    width = summarybox_width,
    alpha = summarybox_alpha,
    geom = summarybox_geom, 
    show.legend = FALSE
  ) +
  stat_compare_means(method = "wilcox", label = "p.signif", label.x.npc = "center", label.y = 16, show.legend = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.1))) +
  labs(y = "log<sub>2</sub>(CLL-like:other)", color = "Patient Type", fill = "Patient Type", x = NULL) +
  theme(axis.title.y.left = ggtext::element_markdown()) + 
  theme(axis.text.x.bottom = element_blank()) +
  theme(axis.ticks.x.bottom = element_blank()) +
  theme(legend.position = "bottom", legend.justification = "center")
cluster_proportion_ratio_plot
