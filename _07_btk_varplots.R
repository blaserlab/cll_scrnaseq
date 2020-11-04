source("00_packages_functions.R")

# join on the cell calls to the main cds


cols_to_add_e1pct <-
  left_join(tbl_df(colData(cds_aligned)), cell_calls_e1pct) %>% 
  select(
    umi_per_cell_wt_e1pct = umi_per_cell_wt,
    max_bc_error_per_cell_wt_e1pct = max_bc_error_per_cell_wt,
    umi_per_cell_mt_e1pct = umi_per_cell_mt,
    max_bc_error_per_cell_mt_e1pct = max_bc_error_per_cell_mt,
    mut_to_wt_e1pct = mut_to_wt,
    cell_call_e1pct = cell_call
  )


cols_to_add_e10pct <-
  left_join(tbl_df(colData(cds_aligned)), cell_calls_e10pct) %>% 
  select(
    umi_per_cell_wt_e10pct = umi_per_cell_wt,
    max_bc_error_per_cell_wt_e10pct = max_bc_error_per_cell_wt,
    umi_per_cell_mt_e10pct = umi_per_cell_mt,
    max_bc_error_per_cell_mt_e10pct = max_bc_error_per_cell_mt,
    mut_to_wt_e10pct = mut_to_wt,
    cell_call_e10pct = cell_call
  )

bind_on_cds<-function(cds, cols_to_add) {
  ncol_cds<-ncol(colData(cds))
  ncol_add<-ncol(cols_to_add)
  colData(cds)[,(ncol_cds+1):(ncol_cds+ncol_add)]<-as.data.frame(cols_to_add)
  return(cds)
}
  
cds_aligned_btk<-bind_on_cds(cds = cds_aligned, cols_to_add = cols_to_add_e1pct)
cds_aligned_btk<-bind_on_cds(cds = cds_aligned_btk, cols_to_add = cols_to_add_e10pct)

#now plot the btk calls using custom var plot (faceted)

btk_e10pct_facet <- custom_variable_plot(
  cds_aligned_btk,
  var = "cell_call_e10pct",
  value_to_highlight = c("wt", "mutant"),
  palette = c("red", "blue"),
  foreground_alpha = 0.4,
  plot_title = "BTK Genotype, <10% error"
) + theme(panel.background = element_rect(color = "grey80")) +
  facet_grid(cols = vars(timepoint), rows = vars(pt)) +
  theme(strip.background = element_blank())
save_plot(
  plot = btk_e10pct_facet,
  filename = "plots_out/btk_e10pct_facet.pdf",
  base_height = 7,
  base_width = 7.5
)

btk_e1pct_facet <- custom_variable_plot(
  cds_aligned_btk,
  var = "cell_call_e1pct",
  value_to_highlight = c("wt", "mutant"),
  palette = c("red", "blue"),
  foreground_alpha = 0.4,
  plot_title = "BTK Genotype, <1% error"
) + theme(panel.background = element_rect(color = "grey80")) +
  facet_grid(cols = vars(timepoint), rows = vars(pt)) +
  theme(strip.background = element_blank())
save_plot(
  plot = btk_e1pct_facet,
  filename = "plots_out/btk_e1pct_facet.pdf",
  base_height = 7,
  base_width = 7.5
)

vtp<-c("cell_call_e10pct","cell_call_e10pct","cell_call_e1pct","cell_call_e1pct")
vth<-c("wt","mutant","wt","mutant")
pal<-c("blue","red","blue","red")
pts<-c("BTK Genotype: WT, <10% error","BTK Genotype: Mutant, <10% error","BTK Genotype: WT, <1% error","BTK Genotype: Mutant, <1% error")

if (!dir.exists("plots_out/btk_overlays/")){
  dir.create("plots_out/btk_overlays/")
}

btk_overlays<-lapply(X = seq_along(pts), 
  FUN = custom_variable_plot,
  cds = cds_aligned_btk,
  var = vtp,
  value_to_highlight = vth,
  palette = pal,
  foreground_alpha = 0.4,
  plot_title = pts,
  legend_title = NULL,
  legend_pos = "none",
  outfile = paste0("btk_overlays/",vtp,"_",vth),
  h = 4,
  w = 4.4,
  cell_size = 0.5
)


save.image.pigz("cll_scrnaseq.RData",n.cores = 39)
