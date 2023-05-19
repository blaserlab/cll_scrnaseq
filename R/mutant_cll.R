
umap_new_england <- bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "B")), "leiden_assignment_binned_2", overwrite_labels = TRUE, foreground_alpha = 0.2)
umap_no_new_england <- bb_var_umap(filter_cds(cds_main, cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "B")), "leiden_assignment_binned", overwrite_labels = TRUE, foreground_alpha = 0.2)

dat1 <- bb_cellmeta(cds_main) |>
  filter(!is.na(btk_mutVaf)) |>
  filter(btk_mutVaf > 0) |>
  mutate(cell_status = "BTK mutant VAF > 0")
dat2 <- bb_cellmeta(cds_main) |>
  filter(!is.na(btk_mutVaf)) |>
  filter(btk_mutVaf == 0) |>
  mutate(cell_status = "BTK mutant VAF = 0")


vaf_density_plot <- ggplot(dat1, aes(x = btk_mutVaf, color = cell_status)) +
  geom_density(adjust = 1/10) +
  geom_density(data = dat2, adjust = 1/10) + 
  labs(x = "BTK mutant VAF")

btk_mutdepth_hist <- bb_cellmeta(cds_main) |> 
  filter(btk_mutdepth > 0) |> 
  ggplot(aes(x = btk_mutdepth)) +
  geom_histogram(binwidth = 1, fill = "transparent", color = "black") +
  scale_x_continuous(
    minor_breaks = seq(0, 20, by = 1),
    breaks = seq(0, 20, by = 5), 
    limits = c(0, 20),
    guide = ggh4x::guide_axis_minor()
  ) + 
  labs(x = "BTK mutant UMIs per Cell")

colData(cds_main)$btk_type0 <- ifelse(colData(cds_main)$btk_mutdepth > 0, "mutant", "WT")
colData(cds_main)$btk_type1 <- ifelse(colData(cds_main)$btk_mutdepth > 1, "mutant", "WT")


# bb_var_umap(cds_main, "density", facet_by = c("sample"))
# bb_var_umap(cds_main, "density", facet_by = c("sample"))
# bb_var_umap(cds_main, "btk_type1", facet_by = "value")
# 
# bb_cellmeta(cds_main) |> filter(!is.na(btk_type1)) |> pull(sample) |> unique()


# bb_var_umap(filter_cds(cds_main,
#                        cells = bb_cellmeta(cds_main) |>
#                          filter(partition_assignment == "B") |> 
#                          filter(btk_type1 %in% c("mutant", "WT"))),
#             "density",
#             facet_by = "btk_type1")

# bb_var_umap(filter_cds(cds_main,
#                        cells = bb_cellmeta(cds_main) |>
#                          filter(partition_assignment == "B") |> 
#                          filter(btk_type0 %in% c("mutant", "WT"))),
#             "density",
#             facet_by = "btk_type0")

btk_type0_density_umap <- bb_var_umap(filter_cds(cds_main,
                       cells = bb_cellmeta(cds_main) |>
                         filter(partition_assignment == "B") |> 
                         filter(btk_type0 %in% c("mutant", "WT"))),
            "density",
            facet_by = c("btk_type0", "timepoint_merged_1"), 
            cols = vars(btk_type0),
            rows = vars(timepoint_merged_1), 
            plot_title = "More than 0 Mutant BTK Reads") +
  panel_border()
btk_type0_density_umap

btk_type1_density_umap <- bb_var_umap(filter_cds(cds_main,
                       cells = bb_cellmeta(cds_main) |>
                         filter(partition_assignment == "B") |> 
                         filter(btk_type1 %in% c("mutant", "WT"))),
            "density",
            facet_by = c("btk_type1", "timepoint_merged_1"), 
            cols = vars(btk_type1),
            rows = vars(timepoint_merged_1), 
            plot_title = "More than 1 Mutant BTK Read") +
  panel_border()
btk_type1_density_umap
# bb_var_umap(filter_cds(cds_main,
#                        cells = bb_cellmeta(cds_main) |>
#                          filter(partition_assignment == "B") |> 
#                          filter(btk_type %in% c("mutant", "WT"))),
#             "density",
#             facet_by = c("btk_type", "timepoint_merged"), 
#             cols = vars(btk_type),
#             rows = vars(timepoint_merged))


btk_mutation_summary_table <- bb_cellmeta(cds_main) |> 
  filter(!is.na(btk_depth)) |> 
  mutate(btk_summary = case_when(
    btk_mutdepth == 0 & btk_depth > 0 ~ "no mutant BTK UMIs but at least 1 WT BTK UMI",
    btk_mutdepth == 1 & btk_depth ==1 ~ "only one BTK UMI and that one is mutant",
    btk_depth %notin% c(0, 1) & btk_depth == btk_mutdepth ~ "more than one BTK UMI and all of them are mutant",
    btk_depth %notin% c(0, 1) & btk_depth <= 2*btk_mutdepth ~ "more than 1 BTK UMI and half or more are mutant",
    btk_depth %notin% c(0, 1) & btk_depth > 2*btk_mutdepth ~ "more than 1 BTK UMI and less than half are mutant",
    TRUE ~ "pending"
  )) |>
  count(btk_summary) |>
  arrange(desc(n))

btk_mutant_distribution <- bb_cellmeta(cds_main) |>
  filter(partition_assignment == "B") |>
  filter(btk_type %in% c("mutant", "WT")) |>
  # filter(sample %in% samples_with_mutant_cells) |>
  count(sample, leiden_assignment_binned_2, btk_type) |>
  pivot_wider(names_from = btk_type,
              values_from = n,
              values_fill = 0) |>
  mutate(percent_mutant = mutant / (WT + mutant) * 100) |>
  ggplot(aes(x = leiden_assignment_binned_2, y = percent_mutant, color = sample)) +
  geom_jitter(width = 0.2) +
  stat_summary(
    fun.data = mean_se,
    geom = "crossbar",
    inherit.aes = FALSE,
    mapping = aes(x = leiden_assignment_binned_2, y = percent_mutant)
  ) +
  labs(x = NULL)





