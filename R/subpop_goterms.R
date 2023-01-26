library(topGO)
mrd1_enrich <- bb_goenrichment(query = cds_main_bcell_subpop_top_markers |> 
                  filter(cell_group == "MRD1") |> 
                  pull(gene_short_name), reference = bb_rowmeta(cds_main), group_pval = 0.01, go_db = "org.Hs.eg.db")

mrd1_summary_0.8 <- bb_gosummary(mrd1_enrich, reduce_threshold = 0.8, go_db = "org.Hs.eg.db")
mrd1_summary_0.9 <- bb_gosummary(mrd1_enrich, reduce_threshold = 0.9, go_db = "org.Hs.eg.db")
save(mrd1_enrich, file = "data/mrd1_enrich.rda", compress = "bzip2")
save(mrd1_summary_0.8, file = "data/mrd1_summary_0.8.rda", compress = "bzip2")
save(mrd1_summary_0.9, file = "data/mrd1_summary_0.9.rda", compress = "bzip2")

mrd2_enrich <- bb_goenrichment(query = cds_main_bcell_subpop_top_markers |> 
                  filter(cell_group == "MRD2") |> 
                  pull(gene_short_name), reference = bb_rowmeta(cds_main), group_pval = 0.01, go_db = "org.Hs.eg.db")

mrd2_summary_0.8 <- bb_gosummary(mrd2_enrich, reduce_threshold = 0.8, go_db = "org.Hs.eg.db")
mrd2_summary_0.9 <- bb_gosummary(mrd2_enrich, reduce_threshold = 0.9, go_db = "org.Hs.eg.db")

save(mrd2_enrich, file = "data/mrd2_enrich.rda", compress = "bzip2")
save(mrd2_summary_0.8, file = "data/mrd2_summary_0.8.rda", compress = "bzip2")
save(mrd2_summary_0.9, file = "data/mrd2_summary_0.9.rda", compress = "bzip2")

btk_enrich <- bb_goenrichment(query = cds_main_bcell_subpop_top_markers |> 
                  filter(cell_group == "BTK") |> 
                  pull(gene_short_name), reference = bb_rowmeta(cds_main), group_pval = 0.01, go_db = "org.Hs.eg.db")

btk_summary_0.8 <- bb_gosummary(btk_enrich, reduce_threshold = 0.8, go_db = "org.Hs.eg.db")
btk_summary_0.9 <- bb_gosummary(btk_enrich, reduce_threshold = 0.9, go_db = "org.Hs.eg.db")

save(btk_enrich, file = "data/btk_enrich.rda", compress = "bzip2")
save(btk_summary_0.8, file = "data/btk_summary_0.8.rda", compress = "bzip2")
save(btk_summary_0.9, file = "data/btk_summary_0.9.rda", compress = "bzip2")

