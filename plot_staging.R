umap_allpts_alltimes <- custom_variable_plot(cds_aligned, var = "pt", foreground_alpha = 0.4)+
  facet_grid(cols = vars(timepoint), rows = vars(value)) + 
  theme(strip.background = element_blank()) + 
  theme(legend.position = "none") + 
  theme(panel.background = element_rect(color = "grey80"))

umap_partition <- custom_variable_plot(cds_aligned, var = "partition", foreground_alpha = 0.2, overwrite_labels = T, group_label_size = 4)
umap_leiden <- custom_variable_plot(cds_aligned, var = "leiden", foreground_alpha = 0.2, overwrite_labels = T, group_label_size = 4)
umap_louvain <- custom_variable_plot(cds_aligned, var = "louvain", foreground_alpha = 0.2, overwrite_labels = T, group_label_size = 4)

umap_seurat_l1 <- custom_variable_plot(cds_aligned, var = "predicted.celltype.l1", foreground_alpha = 0.2, overwrite_labels = T, group_label_size = 4)
umap_seurat_l2 <- custom_variable_plot(cds_aligned, var = "predicted.celltype.l2", foreground_alpha = 0.2, overwrite_labels = T, group_label_size = 4)
