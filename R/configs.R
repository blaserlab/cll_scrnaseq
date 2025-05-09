# conflicts ---------------------------------------------------------------
# resolve conflicting function names here

conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("count", "dplyr")
conflict_prefer("as.data.frame", "base")


# analysis configurations -------------------------------------------------
# use this section to set graphical themes, color palettes, etc.

# graphical parameters####
theme_set(theme_cowplot(font_size = 10))


# show_col(pal_npg("nrc")(10))
experimental_group_palette <- c(
  "MRD" = "#3C5488",
  "BTK" = "#DC0000",
  "responsive" = "#3C5488",
  "sensitive" = "#3C5488",
  "IBS" = "#3C5488",
  "IBR" = "#DC0000",
  "resistant" = "#DC0000",
  "unenriched" = "green4",
  "1" = "white",
  "2" = "grey80",
  "3" = "black",
  "Timepoint 1" = "white",
  "Timepoint 2" = "grey80",
  "Timepoint 3" = "black",
  "B" = brewer.pal(n = 8, name = "Set1")[1],
  "T" = brewer.pal(n = 8, name = "Set1")[2],
  "NK" = brewer.pal(n = 8, name = "Set1")[3],
  "Mono" = brewer.pal(n = 8, name = "Set1")[4],
  "cDC" = brewer.pal(n = 8, name = "Set1")[5],
  "pDC" = brewer.pal(n = 8, name = "Set1")[6],
  "HSPC" = brewer.pal(n = 8, name = "Set1")[7],
  "Prolif" = brewer.pal(n = 8, name = "Set1")[8],
  "Inflammatory 1" = brewer.pal(n = 8, name = "Set2")[1],
  "Inflammatory 2" = brewer.pal(n = 8, name = "Set2")[2],
  "Stressed" = brewer.pal(n = 8, name = "Set2")[3],
  "Activated BCR" = brewer.pal(n = 8, name = "Set2")[4],
  "Malignant Clone" = brewer.pal(n = 8, name = "Set2")[7],
  "non-Malignant Clone" = brewer.pal(n = 8, name = "Set2")[8]
  
)

jitter_alpha_fill <- 0.2
jitter_shape <- 21
jitter_size <- 2
jitter_stroke <- 0.5
jitter_width <- 0.2
jitter_alpha_color <- 1
jitter_height <- 0.2

summarybox_color <- "black"
summarybox_size <- 0.5
summarybox_width <- 0.3
summarybox_alpha <- 0.3
summarybox_geom <- "crossbar"

# 3 color heatmap
heatmap_3_colors <- c("#313695","white","#A50026")

# variable colors

pt_colors <- tibble(pt = paste0("cll_", 1:12), cols = brewer.pal(n = 12, name = "Paired")) |> deframe()

partition_colors <- tibble(pa = unique(colData(cds_main)$partition_assignment), cols = brewer.pal(n = 8, name = "Set1")) |> deframe()

seurat_l2_leiden_consensus_colors <- c("CD4 Naive" = brewer.pal(9, "Set1")[4], 
                                       "CD4 TCM" = brewer.pal(9, "Set1")[5], 
                                       "CD8 TEM" = brewer.pal(9, "Set1")[6], 
                                       "Treg" = brewer.pal(9, "Set1")[7],
                                       "NK" = brewer.pal(9, "Set1")[8])

timepoint_palette <- c("1" = brewer.pal(3, "YlGn")[1],
                       "2" = brewer.pal(3, "YlGn")[2],
                       "3" = brewer.pal(3, "YlGn")[3])



network_out <- fs::path("~/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/Paper/r_out/")
network_tables <- fs::path("~/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/tables")
paper_tables <- fs::path("~/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/Paper/tables")

colData(cds_main)$leiden_comparison_renamed <- recode(colData(cds_main)$leiden_comparison_renamed, "Acivated BCR" = "Activated BCR")
colData(cds_main)$patient_type3 <- case_match(colData(cds_main)$patient_type2, "sensitive" ~ "IBS", "resistant" ~ "IBR")
colData(cds_main)$leiden_enrichment1 <- case_match(colData(cds_main)$leiden_enrichment, "sensitive" ~ "IBS", "resistant" ~ "IBR", .default = colData(cds_main)$leiden_enrichment)
