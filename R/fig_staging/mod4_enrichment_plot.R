# module go term bubbles

mod4_enrichments <- module_enrichment$`Module 4`$res_table |> 
  mutate(cF_numeric = as.numeric(recode(classicFisher, "< 1e-30" = "1e-30"))) |> 
  mutate(neg_log_cF = -log10(cF_numeric)) |> 
  mutate(Term = fct_reorder(Term, neg_log_cF, .desc = FALSE))
mod4_enrichments

mod4_enrichment_plot <- ggplot(mod4_enrichments |> slice_max(neg_log_cF, n = 20), aes(y = Term, x = neg_log_cF, size = Annotated)) +
  geom_point(pch = 21, color = "black", fill = alpha("black", alpha = 0.2)) +
  scale_size_area(limits=c(100, 10000), breaks = (c(300, 1000,3000, 10000))) +
  labs(x = "-log<sub>10</sub>P", y = NULL, size = "Size") +
  theme_minimal_grid(font_size = 10) +
  theme(axis.title.x.bottom = ggtext::element_markdown())
mod4_enrichment_plot
