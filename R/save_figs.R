poster_figs <- fs::path("/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/presentations/figs")

save_plot(fs::path(poster_figs, "umap_leiden_l1.png"), umap_leiden_l1 + theme(legend.position = "top"), base_width = 8.5, base_height = 5.5)

save_plot(fs::path(poster_figs, "umap_subcluster.png"), umap_subcluster, base_width = 4, base_height = 3.5)

save_plot(fs::path(poster_figs, "module_heatmap_bcells.png"), plot_grid(module_heatmap_bcells) + theme(plot.margin = margin(0.1,0.5,0.1,0.1, "cm")), base_width = 6, base_height = 6)

save_plot(fs::path(poster_figs, "subpop_top_markers_heatmap.png"), plot_grid(subpop_top_markers_heatmap) + theme(plot.margin = margin(0.1,0.5,0.1,0.1, "cm")), base_width = 12.5, base_height = 3.75)

save_plot(fs::path(poster_figs, "tcell_subpop_umap.png"), plot_grid(tcell_subpop_umap), base_width = 5.4, base_height = 4.2)

save_plot(fs::path(poster_figs, "treg_pct_plot.png"), plot_grid(treg_pct_plot), base_width = 5.4, base_height = 4.2)

save_plot(fs::path(poster_figs, "tcr_diversity_plot.png"), plot_grid(tcr_diversity_plot), base_width = 5.4, base_height = 4.2)

save_plot(fs::path(poster_figs, "exh_genebub.png"), plot_grid(exh_genebub), base_width = 8, base_height = 5)
