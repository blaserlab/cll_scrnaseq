cd4_tcm_stats <- cd4_tcm_pseudobulk_res$Result |> 
  group_by(gene_short_name) |> 
  mutate(n = n()) |> 
  filter(n == 1) |> 
  select(gene_short_name, stat) |> 
  deframe()

leiden_18_31_stats <- leiden_18_31_pseudobulk_res$Result |> 
  group_by(gene_short_name) |> 
  mutate(n = n()) |> 
  filter(n == 1) |> 
  select(gene_short_name, stat) |> 
  deframe()
cd4_tcell_gsea_plot <- pmap(.l = list(gene_set = c("GSE11057_NAIVE_VS_MEMORY_CD4_TCELL_UP",
                                                   "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
               plot_title = c("Naive vs Memory CD4",
                              "Hallmark TNFA Signaling"),
               pval_posx =  c(5000, 2000),
               pval_posy =  c(0.25, 0.4),
               xlab = rev(c("IBS-CD4 TCM | CD4 Naive", "IBS-CD4 TCM | IBR-CD4 TCM")),
               dat = list(cd4_tcm_stats, leiden_18_31_stats),
               padjdat = list(cd4_tcm_pseudobulk_gsea_res_full, leiden_18_31_pseudobulk_gsea_res_full)),
     .f = \(gene_set, 
            plot_title, 
            pval_posx, 
            pval_posy,
            xlab,
            dat,
            padjdat,
            sets = tcell_leiden_enrichment_genesets) {
       pval <- padjdat |> 
         dplyr::filter(pathway == gene_set) |> 
         pull(padj) |> 
         signif(digits = 2)
       fgsea::plotEnrichment(stats = dat,
                             pathway = sets[[gene_set]]) + 
         labs(subtitle = plot_title,
              x = xlab) +
         theme(plot.subtitle = element_text(hjust = 0.5)) + 
         annotate(geom = "text", label = paste0("padj =\n", pval), x = pval_posx, y = pval_posy)
     }) |> set_names(c("GSE11057_NAIVE_VS_MEMORY_CD4_TCELL_UP",
                       "HALLMARK_TNFA_SIGNALING_VIA_NFKB"))

cd4_tcell_gsea_plots_combined <- cd4_tcell_gsea_plot$GSE11057_NAIVE_VS_MEMORY_CD4_TCELL_UP / 
  cd4_tcell_gsea_plot$HALLMARK_TNFA_SIGNALING_VIA_NFKB +
  patchwork::plot_layout(axes = "collect")

cd4_tcell_gsea_plots_combined


