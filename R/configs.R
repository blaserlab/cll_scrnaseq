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
theme_set(theme_cowplot(font_size = 10))


# show_col(pal_npg("nrc")(10))
experimental_group_palette <- c(
  "MRD" = "#3C5488",
  "BTK" = "#DC0000"
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

timepoint_palette <- c("1" = brewer.pal(3, "YlGn")[1],
                       "2" = brewer.pal(3, "YlGn")[2],
                       "3" = brewer.pal(3, "YlGn")[3])

# cds_mods
colData(cds_main)$timepoint_merged_1 <- recode(colData(cds_main)$timepoint_merged, "baseline" = "1",
                                               "3yrs|btk_clone" = "2",
                                               "5yrs|relapse" = "3")


network_out <- "/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/figs/source"
network_tables <- "/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/tables"


# source local configs ----------------------------------------------------
# these are sourced after main configs and will overwrite duplicate entries if
# present. The file local_configs.R is ignored by git and so is useful for user-
# specific configurations such as output directories or formatting.

fs::file_create(here::here("R/local_configs.R")) # will not overwrite)

source(here::here("R/local_configs.R"))


