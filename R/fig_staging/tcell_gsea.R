tcell_gsea_plot <- pmap(.l = list(gene_set = c("GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_DN",
                            "GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_UP"),
               plot_title = c("PD1-low Up",
                              "PD1-high Up"),
               pval_posx =  c(3000, 3000),
               pval_posy =  c(0.3, -0.3)),
     .f = \(gene_set, 
            plot_title, 
            pval_posx, 
            pval_posy,
            dat = tcell_leiden_enrichment_stats,
            sets = tcell_leiden_enrichment_genesets,
            padjdat = tcell_leiden_enrichment_gsea_res) {
       pval <- padjdat |> 
         dplyr::filter(pathway == gene_set) |> 
         pull(padj) |> 
         signif(digits = 2)
       fgsea::plotEnrichment(stats = dat,
                             pathway = sets[[gene_set]]) + 
         labs(subtitle = plot_title,
              x = "IBS | IBR") +
         theme(plot.subtitle = element_text(hjust = 0.5)) + 
         annotate(geom = "text", label = paste0("padj =\n", pval), x = pval_posx, y = pval_posy)
     }) |> set_names(c("GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_DN",
                       "GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_UP"))

tcell_gsea_plots_combined <- tcell_gsea_plot$GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_DN / 
  tcell_gsea_plot$GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_UP +
  patchwork::plot_layout(axes = "collect")
# tcell_gsea_plots_combined
# tcell_leiden_enrichment_gsea_res |> filter(pathway == "GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_DN") |> pull(padj)
# tcell_leiden_enrichment_gsea_res |> filter(pathway == "GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_UP") |> pull(padj)
