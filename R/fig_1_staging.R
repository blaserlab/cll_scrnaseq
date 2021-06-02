# global umap density faceted--------------------------------
umap_density <- 
  bb_var_umap(
    cds = cds_main[,colData(cds_main)$partition_assignment_2 == "B"],
    var = "log_local_n",
    sample_equally = TRUE,
    cell_size = 1,
    nbin = 100,
    facet_by = c("type_timepoint"),
    nrow = 2,
    foreground_alpha = 0.6
  ) +
  theme(panel.background = element_rect(color = "grey80")) +
  theme(legend.justification = "center") +
  labs(color = "Log10(Local Cells)")
umap_density


# volcano plot MRD1 vs BTK cluster----------------------------------------
genes_to_highlight_MRD1_BTK <- c("CCL3","TWF2", "RAC2")
volcano_MRD1_BTK <- pseudobulk_MRD1_BTK[[2]] %>% 
  filter(str_detect(gene_short_name, "IGH.*|IGK.*|IGL.*", negate = T)) %>%
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58) %>%
  mutate(text_label = ifelse(gene_short_name %in% genes_to_highlight_MRD1_BTK, gene_short_name, "")) %>% 
  ggplot(
    mapping = aes(
      x = log2FoldChange,
      y = -1*log10(padj),
      colour = threshold,
      fill = threshold,
      label = text_label
    )
  ) +
  geom_point(shape = 21, 
             size = 0.5, 
             alpha = 0.4) +
  geom_text_repel(color = "black", 
                  box.padding = 0.5,
                  point.padding = 0.25,
                  min.segment.length = 0,
                  max.overlaps = 20000,
                  nudge_x = -0.5,
                  size = 3, 
                  segment.size = 0.25,
                  force = 2,
                  seed = 1234,
                  segment.curvature = -0.1,
                  segment.square = TRUE,
                  segment.inflect = TRUE) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey80", "#DC0000")) +
  scale_fill_manual(values = c("transparent", "#DC0000")) +
  labs(caption = "\U21D0 Up in BTK Cluster\nUp in MRD1 Cluster \U21D2", title = NULL)+
  theme(plot.caption.position = "panel") +
  theme(plot.caption = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +  
  coord_cartesian(xlim = c(-1.1*max(abs(range(pseudobulk_MRD1_BTK[[2]] %>% 
                                                filter(!is.na(padj)) %>% 
                                                pull(log2FoldChange)))), 
                           1.1*max(abs(range(pseudobulk_MRD1_BTK[[2]] %>% 
                                               filter(!is.na(padj)) %>% 
                                               pull(log2FoldChange))))))



# volcano plot MRD2 vs BTK cluster----------------------------------------
genes_to_highlight_MRD2_BTK <- c("CTLA4", "LILRA4", "FMOD", "TGFBI", "MIR155HG","TXNIP", "PLCB1")
volcano_MRD2_BTK <- pseudobulk_MRD2_BTK[[2]] %>%
  filter(str_detect(gene_short_name, "IGH.*|IGK.*|IGL.*", negate = T)) %>%
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58) %>%
  mutate(text_label = ifelse(gene_short_name %in% genes_to_highlight_MRD2_BTK, gene_short_name, "")) %>% 
  ggplot(
    mapping = aes(
      x = log2FoldChange,
      y = -1*log10(padj),
      colour = threshold,
      fill = threshold,
      label = text_label
    )
  ) +
  geom_point(shape = 21, 
             size = 0.5, 
             alpha = 0.4) +
  geom_text_repel(color = "black", 
                  box.padding = 0.5,
                  point.padding = 0.25,
                  min.segment.length = 0,
                  max.overlaps = 20000,
                  nudge_x = -0.5,
                  size = 3, 
                  segment.size = 0.25,
                  force = 2,
                  seed = 1234,
                  segment.curvature = -0.1,
                  segment.square = TRUE,
                  segment.inflect = TRUE) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey80", "#DC0000")) +
  scale_fill_manual(values = c("transparent", "#DC0000")) +
  labs(caption = "\U21D0 Up in BTK Cluster\nUp in MRD2 Cluster \U21D2", title = NULL)+
  theme(plot.caption.position = "panel") +
  theme(plot.caption = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +  
  coord_cartesian(xlim = c(-1.1*max(abs(range(pseudobulk_MRD2_BTK[[2]] %>% 
                                                filter(!is.na(padj)) %>% 
                                                pull(log2FoldChange)))), 
                           1.1*max(abs(range(pseudobulk_MRD2_BTK[[2]] %>% 
                                               filter(!is.na(padj)) %>% 
                                               pull(log2FoldChange))))))

# volcano plot BTK over time:  baseline v btk----------------------
# pseudobulk_bcell_btk_timepoints[[1]][[2]] %>% filter(padj<0.05, abs(log2FoldChange) >= 0.58) %>% View()# baseline v btk
# pseudobulk_bcell_btk_timepoints[[2]][[2]] %>% filter(padj<0.05, abs(log2FoldChange) >= 0.58) %>% View()# baseline v relapse:  nothing
# pseudobulk_bcell_btk_timepoints[[3]][[2]] %>% filter(padj<0.05, abs(log2FoldChange) >= 0.58) %>% View()# btk v relapse: nothing

genes_to_highlight_bcell_BTK_timepoints <- c("DGKG", "PDE4D", "TCF7", "NFKBIA", "EBF1", "SSPN", "DENND3")
volcano_BTK_bcells <- 
  pseudobulk_bcell_btk_timepoints[[1]][[2]] %>%
  # filter(str_detect(gene_short_name, "IGH.*|IGK.*|IGL.*", negate = T)) %>%
  mutate(threshold = padj < 0.05 & abs(log2FoldChange) >= 0.58) %>%
  mutate(text_label = ifelse(gene_short_name %in% genes_to_highlight_bcell_BTK_timepoints, gene_short_name, "")) %>% 
  ggplot(
    mapping = aes(
      x = log2FoldChange,
      y = -1*log10(padj),
      colour = threshold,
      fill = threshold,
      label = text_label
    )
  ) +
  geom_point(shape = 21, 
             size = 0.5, 
             alpha = 0.4) +
  geom_text_repel(color = "black", 
                  box.padding = 0.5,
                  point.padding = 0.25,
                  min.segment.length = 0,
                  max.overlaps = 20000,
                  nudge_x = -0.5,
                  size = 3, 
                  segment.size = 0.25,
                  force = 2,
                  seed = 1234,
                  segment.curvature = -0.1,
                  segment.square = TRUE,
                  segment.inflect = TRUE) +
  xlab("log2 fold change") +
  ylab("-log10 adjusted p") +
  theme(legend.position = "none") +
  scale_color_manual(values = c("grey80", "#DC0000")) +
  scale_fill_manual(values = c("transparent", "#DC0000")) +
  labs(caption = "\U21D0 Up in baseline\nUp in BTK clone \U21D2", title = NULL)+
  theme(plot.caption.position = "panel") +
  theme(plot.caption = element_text(hjust = 0.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +  
  coord_cartesian(xlim = c(-1.1*max(abs(range(pseudobulk_bcell_btk_timepoints[[1]][[2]] %>% 
                                                filter(!is.na(padj)) %>% 
                                                pull(log2FoldChange)))), 
                           1.1*max(abs(range(pseudobulk_bcell_btk_timepoints[[1]][[2]] %>% 
                                               filter(!is.na(padj)) %>% 
                                               pull(log2FoldChange))))))

# module heatmap-----------------------------------------------------

col_fun_heatmap_bcells <- 
  colorRamp2(
    breaks = c(min(agg_mat_bcells_type_timepoint),
               0,
               max(agg_mat_bcells_type_timepoint)),
    colors = heatmap_3_colors
  )


module_heatmap_bcells <-
  grid.grabExpr(draw(
    Heatmap(matrix = agg_mat_bcells_type_timepoint,
            name = "Module\nExpression",
            column_split = c(rep("BTK", times = 3), rep("MRD", times = 3)),
            col = col_fun_heatmap_bcells,
            row_names_gp = gpar(fontsize = 10),
            column_names_gp = gpar(fontsize = 10),
            column_dend_height = unit(3,"mm"), 
            row_dend_width = unit(3,"mm")
            )
  , padding = unit(c(2,5,2,15),"mm")))# bltr




# go term bubbles--------------------------------------------------

mod4_bubble <- bb_goscatter(simMatrix = module_summary_0.9$`Module 4`$simMatrix, 
             reducedTerms = module_summary_0.9$`Module 4`$reducedTerms,size = "score",
             addLabel = T,
             labelSize = 4) 
