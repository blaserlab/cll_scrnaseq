source("00_packages_functions.R")

l10upcwe10pct<-log10(colData(cds_aligned_btk)$umi_per_cell_wt_e10pct)
colData(cds_aligned_btk)$l10upcwe10pct<-l10upcwe10pct
btk_e10pct_wt_expr_facet <- custom_variable_plot(
  cds_aligned_btk,
  var = "l10upcwe10pct",
  value_to_highlight = l10upcwe10pct[!is.na(l10upcwe10pct)],
  palette = NULL,
  foreground_alpha = 0.4,
  plot_title = "WT BTK Expr., <10% error",
  legend_title = "Log10(Expr.)"
) + theme(panel.background = element_rect(color = "grey80")) +
  facet_grid(cols = vars(timepoint), rows = vars(pt)) +
  theme(strip.background = element_blank())
save_plot(
  plot = btk_e10pct_wt_expr_facet,
  filename = "plots_out/btk_e10pct_wt_expr_facet.pdf",
  base_height = 7,
  base_width = 7.5
)



l10upcme10pct<-log10(colData(cds_aligned_btk)$umi_per_cell_mt_e10pct)
colData(cds_aligned_btk)$l10upcme10pct<-l10upcme10pct
btk_e10pct_mt_expr_facet <- custom_variable_plot(
  cds_aligned_btk,
  var = "l10upcme10pct",
  value_to_highlight = l10upcme10pct[!is.na(l10upcme10pct)],
  palette = NULL,
  foreground_alpha = 0.4,
  plot_title = "Mutant BTK Expr., <10% error",
  legend_title = "Log10(Expr.)"
) + theme(panel.background = element_rect(color = "grey80")) +
  facet_grid(cols = vars(timepoint), rows = vars(pt)) +
  theme(strip.background = element_blank())
save_plot(
  plot = btk_e10pct_mt_expr_facet,
  filename = "plots_out/btk_e10pct_mt_expr_facet.pdf",
  base_height = 7,
  base_width = 7.5
)


l1upcwe10pct<-log10(colData(cds_aligned_btk)$umi_per_cell_wt_e1pct)
colData(cds_aligned_btk)$l10upcwe1pct<-l10upcwe1pct
btk_e1pct_wt_expr_facet <- custom_variable_plot(
  cds_aligned_btk,
  var = "l10upcwe1pct",
  value_to_highlight = l10upcwe1pct[!is.na(l10upcwe1pct)],
  palette = NULL,
  foreground_alpha = 0.4,
  plot_title = "WT BTK Expr., <10% error",
  legend_title = "Log10(Expr.)"
) + theme(panel.background = element_rect(color = "grey80")) +
  facet_grid(cols = vars(timepoint), rows = vars(pt)) +
  theme(strip.background = element_blank())
save_plot(
  plot = btk_e10pct_wt_expr_facet,
  filename = "plots_out/btk_e1pct_wt_expr_facet.pdf",
  base_height = 7,
  base_width = 7.5
)

l10upcme1pct<-log10(colData(cds_aligned_btk)$umi_per_cell_mt_e1pct)
colData(cds_aligned_btk)$l10upcme1pct<-l10upcme1pct
btk_e1pct_mt_expr_facet <- custom_variable_plot(
  cds_aligned_btk,
  var = "l10upcme1pct",
  value_to_highlight = l10upcme1pct[!is.na(l10upcme1pct)],
  palette = NULL,
  foreground_alpha = 0.4,
  plot_title = "Mutant BTK Expr., <10% error",
  legend_title = "Log10(Expr.)"
) + theme(panel.background = element_rect(color = "grey80")) +
  facet_grid(cols = vars(timepoint), rows = vars(pt)) +
  theme(strip.background = element_blank())
save_plot(
  plot = btk_e1pct_mt_expr_facet,
  filename = "plots_out/btk_e1pct_mt_expr_facet.pdf",
  base_height = 7,
  base_width = 7.5
)

save.image.pigz("cll_scrnaseq.RData",n.cores = 39)
