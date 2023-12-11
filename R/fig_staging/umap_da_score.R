

umap_da_score <- bb_var_umap(cds_main, "da_score", 
                             rasterize = TRUE, legend_pos = "top") + 
  scale_fill_gradient2(low = "#DC0000", 
                       high = "#3C5488", 
                       mid = "white") +  
  scale_color_gradient2(low = "#DC0000", 
                       high = "#3C5488", 
                       mid = "white", 
                       guide = "none") + 
  labs(fill = "Differential\nAbundance")
umap_da_score
