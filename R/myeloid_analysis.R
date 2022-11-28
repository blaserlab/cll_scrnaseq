# make a subset CDS for the myeloid analysis
# this is the intersection of seurat assignments and monocle clustering
SingleCellExperiment::mainExpName(cds_main) <- "Gene Expression"
cds_myeloid <- filter_cds(cds_main, 
                          cells = bb_cellmeta(cds_main) |> 
                            filter(seurat_celltype_l2 %in% c("CD14 Mono", "CD16 Mono", "cDC1", "cDC2", "pDC")) |> 
                            filter(partition_assignment %in% c("Mono", "cDC", "pDC", "HSPC")))
# reorder the groups
colData(cds_myeloid)$btk_type <- factor(colData(cds_myeloid)$btk_type, levels = c("mutant", "WT", "no_type"))

myeloid_seurat_cluster_umap <- rasterise(bb_var_umap(cds_myeloid, "seurat_celltype_l2", overwrite_labels = T, palette = brewer.pal(n = 8, "Set2")), dpi = 300)
myeloid_btk_umap <- rasterise(bb_var_umap(cds_myeloid, "btk_type", facet_by = "value", palette = brewer.pal(n = 8, "Set1")), dpi = 300)
myeloid_genebubbles <- bb_genebubbles(cds_myeloid,
               genes = c("CD14", 
                         "CXCL8", 
                         "MRC1", 
                         "ITGAM", 
                         "TEK", 
                         "CD163", 
                         "CD86", 
                         "HLA-DRA", 
                         "CCR2", 
                         "CCR7", 
                         "CCR5", 
                         "CX3CR1", 
                         "THBD", 
                         "ITGAX", 
                         "CD68", 
                         "TFRC", 
                         "FCGR1A", 
                         "FCGR3A"), 
               cell_grouping = "seurat_celltype_l2", 
               scale_expr = TRUE, return_value = "plot")


btk_breakdown_tbl <- bb_cellmeta(cds_myeloid) |> 
  group_by(seurat_celltype_l2, btk_type) |> 
  count() |> 
  filter(btk_type != "no_type") |> 
  pivot_wider(names_from = btk_type, 
              values_from = n, 
              values_fill = 0) |> 
  mutate(percent_mutant = mutant/(WT + mutant) * 100)
  

btk_timepoint_breakdown_tbl <- bb_cellmeta(cds_myeloid) |> 
  group_by(seurat_celltype_l2, timepoint, btk_type) |> 
  count() |> 
  filter(btk_type != "no_type") |> 
  pivot_wider(names_from = btk_type, 
              values_from = n, 
              values_fill = 0) |> 
  mutate(percent_mutant = mutant/(WT + mutant) * 100)
  

bb_cellmeta(cds_myeloid) |> glimpse()

model_data <- bb_cellmeta(cds_myeloid) |> 
  count(sample, patient, timepoint, seurat_celltype_l2, btk_type) |>
  filter(btk_type != "no_type") |> 
  pivot_wider(names_from = btk_type, 
              values_from = n,
              values_fill = 0) |> 
  mutate(percent_mutant= mutant/(WT+mutant) * 100)


lme_fit <- lmerTest::lmer(formula = percent_mutant ~ timepoint + seurat_celltype_l2 + (1|patient), data = model_data)

summary(lme_fit)
