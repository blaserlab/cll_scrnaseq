#PMID: 21383243
#PMID: 36653453
#DOI: 10.1126/sciimmunol.aba7918
#PMID: 31591533
# 
# tcell_leiden_encrichment_le <-
#   tcell_leiden_enrichment_gsea_res |> filter(
#     pathway %in% c(
#       "GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_DN",
#       "GSE26495_PD1HIGH_VS_PD1LOW_CD8_TCELL_UP"
#     )
#   ) |>
#   unnest(cols = c(leadingEdge)) |>
#   select(pathway, leadingEdge)

tcell_leiden_enrichment_examplegenes <-
  c("TOX", "PDCD1", "TBX21", "GZMB", "PRF1", "KLRD1", "LTB", "TIGIT")
# tcell_leiden_enrichment_examplegenes %in% tcell_leiden_encrichment_le$leadingEdge





tcell_leiden_enrichment_genebubdat <-
  bb_genebubbles(
    filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |>
        filter(
          partition_assignment == "T",
          seurat_l2_leiden_consensus == "CD8 TEM"
        )
    ),
    genes = tcell_leiden_enrichment_examplegenes,
    cell_grouping = "leiden_enrichment1",
    return_value = "data",
    scale_expr = TRUE
  )

tcell_leiden_enrichment_genebub <-
  ggplot(
    tcell_leiden_enrichment_genebubdat,
    aes(
      x = leiden_enrichment1,
      y = gene_short_name,
      fill = expression,
      size = proportion
    )
  ) +
  geom_point(pch = 21) +
  scale_fill_viridis_c() +
  scale_size_area() +
  coord_flip() +
  theme_minimal_grid(font_size = 10) + 
  theme(axis.text.x = element_text(
    face = "italic",
    angle = 30,
    hjust = 1
  )) +
  labs(x = NULL,
       y = NULL,
       fill = "Expression",
       size = "Proportion") +
  theme(axis.text.y = element_text(angle = 90, hjust = 0.5))
tcell_leiden_enrichment_genebub
