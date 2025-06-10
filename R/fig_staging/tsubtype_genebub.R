tsubtype_genes <- readr::read_csv("~/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/Paper/priya_gene_lists/tcell_genes.csv") |>
  filter(label != "Immune surveillance") |> 
  mutate(label2 = case_match(label,
                             "Treg" ~ "T<sub>reg</sub>",
                             "CD4 naïve/CM resting cells" ~ "CD4 naïve/CM",
                             "Activated T cells" ~ "Activation",
                             "Proliferating T cells" ~ "Proliferation",
                             "CD8 cytotoxic T cells" ~ "Cytotoxicity",
                             "Exhaustion markers" ~ "Exhaustion"))

tsubtype_genebubdat <-
  bb_genebubbles(
    filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |>
        filter(
          partition_assignment == "T"
        )
    ),
    genes = tsubtype_genes$gene,
    cell_grouping = "seurat_l2_leiden_consensus",
    return_value = "data",
    scale_expr = TRUE
  ) |> 
  left_join(tsubtype_genes, by = join_by(gene_short_name == gene))

tsubtype_genebub <-
  ggplot(
    tsubtype_genebubdat,
    aes(
      x = seurat_l2_leiden_consensus,
      y = gene_short_name,
      fill = expression,
      size = proportion
    )
  ) +
  geom_point(pch = 21) +
  scale_fill_viridis_c() +
  scale_size_area(max_size = 6) +
  theme_minimal_grid(font_size = 10) + 
  theme(axis.text.y = element_text(
    face = "italic"
  )) +
  labs(x = NULL,
       y = NULL,
       fill = "Expression",
       size = "Proportion") + 
  facet_grid(label2~., scales = "free", space = "free", switch = "y") + 
  theme(strip.placement = "outside") + 
  theme(strip.text = ggtext::element_markdown())
tsubtype_genebub
