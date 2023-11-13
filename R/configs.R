# conflicts ---------------------------------------------------------------
# resolve conflicting function names here

conflict_prefer("filter", "dplyr")
conflict_prefer("lag", "dplyr")
conflict_prefer("select", "dplyr")
conflict_prefer("rename", "dplyr")
conflict_prefer("count", "dplyr")


# analysis configurations -------------------------------------------------
# use this section to set graphical themes, color palettes, etc.

# graphical parameters####
theme_set(theme_cowplot(font_size = 18, line_size = 1))


# show_col(pal_npg("nrc")(10))
experimental_group_palette <- c(
  "MRD" = "#3C5488",
  "BTK" = "#DC0000",
  "responsive" = "#3C5488",
  "sensitive" = "#3C5488",
  "resistant" = "#DC0000"
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

tcell_genes <- read_csv("/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/queries/tcell_genes.csv")

NK_genes <- read_csv("/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/queries/NK_genes.csv") 

tcell_genes_palette <- tcell_genes |> group_by(label) |> summarise() |> mutate(col = brewer.pal(7, "Accent")) |> deframe()

nk_genes_palette <- NK_genes |> group_by(label) |> summarise() |> mutate(col = brewer.pal(3, "Accent")[1:2]) |> deframe()

# cds_mods
colData(cds_main)$timepoint_merged_1 <- recode(colData(cds_main)$timepoint_merged, "baseline" = "1",
                                               "3yrs|btk_clone" = "2",
                                               "5yrs|relapse" = "3")
colData(cds_main)$timepoint_merged_2 <- recode(colData(cds_main)$timepoint_merged, "baseline" = "Timepoint 1",
                                               "3yrs|btk_clone" = "Timepoint 2",
                                               "5yrs|relapse" = "Timepoint 3")

colData(cds_main)$patient_type1 <- recode(colData(cds_main)$patient_type, "BTK" = "resistant", "MRD" = "responsive")
colData(cds_main)$patient_type2 <- recode(colData(cds_main)$patient_type, "BTK" = "resistant", "MRD" = "sensitive")

leiden_l1_assignment <- bb_cellmeta(cds_main) |> 
  count(leiden, seurat_celltype_l1) |> 
  group_by(leiden) |> 
  slice_max(order_by = n, n = 1) |> 
  select(leiden, seurat_celltype_l1) |> 
  deframe()

colData(cds_main)$leiden_l1_assignment <- recode(colData(cds_main)$leiden, !!!leiden_l1_assignment)

seurat_l2_leiden_consensus <- bb_cellmeta(cds_main) |> 
  group_by(leiden, seurat_celltype_l2) |> 
  summarise(n = n()) |> 
  slice_max(order_by = n, n = 1) |> 
  select(leiden, seurat_celltype_l2) |> 
  deframe()

colData(cds_main)$seurat_l2_leiden_consensus <- recode(colData(cds_main)$leiden, !!!seurat_l2_leiden_consensus)
colData(cds_main)$seurat_l2_leiden_consensus <- as.character(colData(cds_main)$seurat_l2_leiden_consensus)

colData(cds_main)$leiden_assignment_binned_2 <- recode(colData(cds_main)$leiden,
                                                "24" = "inflammatory",
                                                "20" = "stressed",
                                                "27" = "CLL-like",
                                                "3" = "CLL-like",
                                                "16" = "CLL-like",
                                                "40" = "CLL-like",
                                                "29" = "CLL-like",
                                                "5" = "CLL-like",
                                                "6" = "CLL-like",
                                                "2" = "CLL-like",
                                                "4" = "CLL-like",
                                                "9" = "CLL-like",
                                                "10" = "CLL-like",
                                                "11" = "CLL-like",
                                                "1" = "CLL-like",
                                                "32" = "New England",
                                                "21" = "New England",
                                                "34" = "CLL-like"
                                                )

network_out <- fs::path("/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/figs/new_figs")
network_tables <- fs::path("/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/tables")
paper_tables <- fs::path("/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/Paper/tables")

# source local configs ----------------------------------------------------
# these are sourced after main configs and will overwrite duplicate entries if
# present. The file local_configs.R is ignored by git and so is useful for user-
# specific configurations such as output directories or formatting.



