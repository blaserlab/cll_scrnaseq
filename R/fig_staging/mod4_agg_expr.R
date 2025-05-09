mod4_agg_expr <-
  aggregate_gene_expression(
    cds = filter_cds(
      cds_main,
      cells = bb_cellmeta(cds_main) |> filter(partition_assignment == "B")
    ),
    gene_group_df = bb_rowmeta(cds_main) |> 
      select(feature_id, module_labeled)
  )

mod4_agg_expr_violin <- t(mod4_agg_expr)[,"Module 4"] |> 
  enframe(name = "cell_id", value = "mod4_expr") |> 
  left_join(bb_cellmeta(cds_main)) |> 
  ggplot(aes(x = patient_type3, y = mod4_expr, color = patient_type3, fill = patient_type3)) + 
  geom_violin(alpha = 0.4, color = "black", draw_quantiles = 0.5) + 
  # geom_jitter(pch = 21) + 
  scale_fill_manual(values = experimental_group_palette) + 
  facet_wrap(~timepoint_merged_1, ncol = 3) +
  stat_compare_means(method = "t.test", label = "p.signif", label.x.npc = 0.5) + 
  theme(legend.position = "none") + 
  theme(strip.background = element_blank()) + 
  labs(x = NULL, y = "Module 4 Expression") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))
  
