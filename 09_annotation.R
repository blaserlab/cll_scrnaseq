
#manually assign celltypes to partitions
cds_aligned$partition<-monocle3::partitions(cds_aligned)
cds_aligned$partition_assignment<-monocle3::partitions(cds_aligned)
cds_aligned$cluster<-monocle3::clusters(cds_aligned)
cds_aligned$cluster_assignment<-monocle3::clusters(cds_aligned)

custom_cp_plot(cds = cds_aligned,cp = "partition",group_label_size = 5)
plot_cells_alt(cds = cds_aligned, gene_or_genes = "CD14")

cds_aligned$partition_assignment<-recode(cds_aligned$partition, 
                                         "1" = "B1",
                                         "2" = "T", 
                                         "3" = "B2",
                                         "4" = "Mono1", 
                                         "5" = "Mono2", 
                                         "6" = "Plt",
                                         "7" = "B3", 
                                         "8" = "Erythrocyte",
                                         "9" = "DC",
                                         "10" = "B4")

pa_lut<-tbl_df(colData(cds_aligned)) %>% 
  select(partition,partition_assignment) %>% 
  unique() %>% 
  arrange(partition)


write.csv(marker_test_res_c %>% 
            rename(cluster = cell_group) %>% 
            arrange(cluster),
          file = "data_out/cluster_top_markers.csv")

write.csv(marker_test_res_p %>% 
            rename(partition = cell_group) %>% 
            left_join(.,pa_lut) %>% 
            arrange(partition),
          file = "data_out/partition_top_markers.csv")

# cluster/partition plots
if (!dir.exists("plots_out")) {
  dir.create("plots_out")
}

#plot all together
cluster_plot_all<-custom_cp_plot(cds = cds_aligned,
               alpha = 0.2,
               plot_title = "All Samples",
               cp = "partition",
               group_label_size = 5,
               outfile = "plots_out/cluster_plot_all.pdf",
               h = 4,
               w = 4.4
)

cluster_plot_all_faceted <- custom_cp_plot(
  cds = cds_aligned,
  alpha = 0.2,
  cp = "partition",
  group_label_size = 3,
  legend_pos = "none"
) + theme(panel.background = element_rect(color = "grey80"))+
  facet_grid(rows = vars(pt), cols = vars(timepoint)) +
  theme(strip.background = element_blank())
save_plot(
  cluster_plot_all_faceted,
  filename = "plots_out/cluster_plot_all_faceted.pdf",
  base_height = 6.5,
  base_width = 7
)

#leiden clusters 
leiden_cluster_plot_all <- custom_cp_plot(
  cds = cds_aligned,
  alpha = 0.2,
  plot_title = "All Samples",
  cp = "cluster",
  group_label_size = 5,
  outfile = "plots_out/leiden_cluster_plot_all.pdf",
  h = 4,
  w = 4.4
)

#marker gene plots
demo_genes<-c("CD79A","FCER2","TCL1A","CD3E","KLRB1","LYZ","CD14","PF4", "HBB")
marker_genes<-plot_cells_alt(
  cds = cds_aligned,
  gene_or_genes = demo_genes,
  alpha = 0.2,
  ncol = 3,
  outfile = "plots_out/marker_genes.png",
  h = 6.5,
  w = 7.5,
  plot_type = "png"
)
