
cd16_le <- cd16_pseudobulk_gsea_res_full |> 
  filter(
    pathway %in% c(
"GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP",
                       "MOSERLE_IFNA_RESPONSE"
    )
  ) |>
  unnest(cols = c(leadingEdge)) |>
  pull(leadingEdge)

cd14_le <- cd14_pseudobulk_gsea_res_full |> 
  filter(
    pathway %in% c(
"GSE9988_LOW_LPS_VS_VEHICLE_TREATED_MONOCYTE_UP",
                       "MOSERLE_IFNA_RESPONSE"
    )
  ) |>
  unnest(cols = c(leadingEdge)) |>
  pull(leadingEdge)


mono_examplegenes <- c(cd14_le, cd16_le)

# mono_examplegenes <-
#   c("CXCL8", "CCL4", "CXCL1","STAT1", "MX1", "CD274", "IFIT2", "IFIT3")
# tcell_leiden_enrichment_examplegenes %in% tcell_leiden_encrichment_le$leadingEdge

# bb_cellmeta(cds_main) |> count(mono_response_merged)


colData(cds_main)$mono_response_merged <-
  case_match(
    colData(cds_main)$mono_response_merged,
    "resistant" ~ "IBR",
    "sensitive" ~ "IBS",
    .default = colData(cds_main)$mono_response_merged
  )

mono_genebubdat <-
  bb_genebubbles(
    filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |>
        filter(
          seurat_l2_leiden_consensus %in% c("CD14 Mono", "CD16 Mono"),
          !is.na(mono_response_merged)
        )
    ),
    genes = mono_examplegenes[mono_examplegenes %in% c("IFI44", "IFI44L", "IFI16", "STAT1", "MX1", "G0S2", "PLK3", "MAFF", "IER3", "CXCL8", "SAMD9", "ICAM1", "REL", "IFIT2", "IFIT3")],
    # genes = mono_examplegenes ,
    cell_grouping = c("seurat_l2_leiden_consensus", "mono_response_merged"),
    return_value = "data",
    scale_expr = TRUE
  )

mono_genebub <-
  ggplot(
    mono_genebubdat,
    aes(
      x = gene_short_name,
      y = mono_response_merged,
      fill = expression,
      size = proportion
    )
  ) +
  geom_point(pch = 21) +
  scale_fill_viridis_c() +
  scale_size_area() +
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
  theme(axis.text.y = element_text(angle = 30, hjust = 1)) + 
  facet_wrap(~seurat_l2_leiden_consensus, ncol = 1)
mono_genebub
