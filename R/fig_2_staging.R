source("R/configs")

colData(cds_main)


umap_clonal_fraction <- bb_var_umap(cds_main[, colData(cds_main)$leiden_assignment_1 %in% c("BTK cluster", "MRD1 cluster", "MRD2 cluster")],
            var = "clone_proportion",legend_title = "Clonal\nFraction", facet_by = "type_timepoint") + theme(panel.background = element_rect(color = "grey80")) 
bb_var_umap(cds_main[,colData(cds_main)$leiden_assignment_1 %in% c("BTK cluster", "MRD1 cluster", "MRD2 cluster")],
            var = "clone_proportion", facet_by = "specimen")

colData(cds_main) %>% 
  as_tibble() %>% 
  group_by(specimen, patient_type, timepoint, specimen_shannon) %>% 
  summarise() %>% 
  filter(!is.na(specimen_shannon)) %>%
  mutate(is_outlier = rstatix::is_outlier(specimen_shannon)) %>% View()
  filter(!is_outlier) %>%
  ggplot(mapping = aes(x = timepoint, y = specimen_shannon, fill = specimen, color = specimen)) +
  geom_point(pch = 21) +
  facet_wrap(facets = vars(patient_type))
