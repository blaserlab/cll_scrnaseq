bcr_genes <- c(
"AKT3",
"BCL10",
"BCL2",
"BLNK",
"BTK",
"BIRC3",
"BIRC5",
"CARD11",
"CD19",
"CD22",
"CD72",
"CD79A",
"CD79B",
"CD80",
"CD81",
"CHP1",
"CHP2",
"CHUK",
"CR2",
"FCGR2B",
"FOS",
"GRB2",
"GSK3B",
"HRAS",
"HRASLS2",
"IFITM1",
"IKBKB",
"IKBKG",
"JUN",
"KRAS",
"LILRB3",
"LYN",
"MALT1",
"MAP2K1",
"MAP2K2",
"DAPP1",
"MAPK1",
"MAPK3",
"NFAT5",
"NFATC1",
"NFATC2",
"NFATC3",
"NFATC4",
"NFKB1",
"NFKBIA",
"NFKBIB",
"NFKBIE",
"NRAS",
"PIK3AP1",
"PIK3CA",
"PIK3CB",
"PIK3CD",
"PIK3CG",
"PIK3R1",
"PIK3R2",
"PIK3R3",
"PIK3R5",
"PLCG2",
"PPP3CA",
"PPP3CB",
"PPP3CC",
"PPP3R1",
"PPP3R2",
"PRKCB",
"PTPN6",
"RAC1",
"RAC2",
"RAF1",
"RASGRP3",
"RELA",
"SOS1",
"SOS2",
"SYK",
"VAV1",
"VAV2",
"VAV3")
bb_cellmeta(cds_main) |> glimpse()
colData(cds_main)$type_timepoint <- paste0(colData(cds_main)$patient_type3, "\n", colData(cds_main)$timepoint_merged_1)


bcr_sig_mat <- bb_aggregate(
  filter_cds(
    cds_main,
    cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "B"),
    genes = bb_rowmeta(cds_main) |> filter(gene_short_name %in% bcr_genes)
  ),
  cell_group_df = bb_cellmeta(cds_main) |>
    filter(partition_assignment == "B") |>
    select(cell_id, type_timepoint),
  scale_agg_values = FALSE
) |> as.matrix() |> t() |> scale() |> t()

bcr_sig_mat <- bcr_sig_mat[!apply(is.nan(bcr_sig_mat), 1, any), ]

new_rownames <- rownames(bcr_sig_mat) |> 
  as_tibble() |> 
  left_join(bb_rowmeta(cds_main), by = join_by(value == feature_id)) |> 
  pull(gene_short_name)

rownames(bcr_sig_mat) <- new_rownames

bcr_sig_sh <- SummarizedHeatmap(bcr_sig_mat)
colData(bcr_sig_sh)$type <- str_extract(rownames(colData(bcr_sig_sh)), "IBR|IBS")
colData(bcr_sig_sh)$timepoint <- str_extract(rownames(colData(bcr_sig_sh)), "1|2|3")
colData(bcr_sig_sh)

p1 <- bb_plot_heatmap_main(bcr_sig_sh, low = "#762A83", mid = "#F7F7F7", high = "#1B7837") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
p2 <- bb_plot_heatmap_colDendro(bcr_sig_sh)
p3 <- bb_plot_heatmap_rowDendro(bcr_sig_sh)
p4 <- bb_plot_heatmap_colData(bcr_sig_sh, 
                              side = "bottom", 
                              manual_pal = experimental_group_palette, vars = c("Patient Type" = "type", "Timepoint" = "timepoint"))
p5 <- guide_area()

design <- "
#25
315
#45
"
bcr_sig_hm <- wrap_plots(p1, p2, p3, p4, p5, design = design, heights = c(0.05, 1, 0.075), widths = c(0.1, 1, 0.2), guides = "collect")
