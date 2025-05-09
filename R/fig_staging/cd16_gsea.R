cd16_gsea_plot <- pmap(.l = list(gene_set = c("GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP",
                            "MOSERLE_IFNA_RESPONSE"),
               plot_title = c("LPS Treated Monocytes Up",
                              "IFNA Response"),
               pval_posx =  c(3000, 3000),
               pval_posy =  c(0.3, -0.6)),
     .f = \(gene_set, 
            plot_title, 
            pval_posx, 
            pval_posy,
            dat = cd16_enrichment_stats,
            sets = cd16_genesets,
            padjdat = cd16_pseudobulk_gsea_res_full) {
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
     }) |> set_names(c("GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP",
                       "MOSERLE_IFNA_RESPONSE"))


cd16_gsea_plots_combined <- cd16_gsea_plot$GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP / 
  cd16_gsea_plot$MOSERLE_IFNA_RESPONSE +
  patchwork::plot_layout(axes = "collect")
