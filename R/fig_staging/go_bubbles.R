plotfun <-
  function(data, title) {
    ggplot(data, aes(
      x = V1,
      y = V2,
      size = size,
      fill = parentTerm
    )) +
      geom_point(pch = 21,
                 color = "black", 
                 aes(alpha = score)) +
      scale_fill_brewer(palette = "Set1") +
      guides(alpha = guide_legend(order = 2, override.aes=list(shape = 19))) +
      guides(fill = guide_legend(order = 1, override.aes=list(size = 3)))+
      guides(size = guide_legend(order = 3)) +
      labs(x = "PCoA 1", 
           y = "PCoA 2", 
           fill = "Parent GO-Term",
           size = "Term Size",
           alpha = "-log10 P", 
           title = title) +
      theme(legend.box = "horizontal") +
      theme(aspect.ratio = 0.9) +
      theme(plot.title = element_text(hjust = 0.5))
      
  }

# stressed
stressed_goscatter <- bb_goscatter(simMatrix = leiden_20_summary_0.9$simMatrix, 
             reducedTerms = leiden_20_summary_0.9$reducedTerms)[["data"]] |> 
  as_tibble() |> 
  mutate(parentTerm = recode(parentTerm, "biological process involved in interspecies interaction between organisms" = "interspecies interaction")) |> 
  filter(size >5) |> 
  plotfun(title = "Stressed") 

# inflammatory 1
infl1_goscatter <- bb_goscatter(simMatrix = leiden_24_summary_0.9$simMatrix, reducedTerms = leiden_24_summary_0.9$reducedTerms)[["data"]] |> 
  as_tibble() |> 
  mutate(parentTerm = recode(parentTerm, "antigen processing and presentation of endogenous antigen" = "antigen processing and presentation")) |> 
  mutate(parentTerm = recode(parentTerm, "positive regulation of lymphocyte proliferation" = "positive regulation of lymph. prolif.")) |> 
  mutate(parentTerm = recode(parentTerm, "positive regulation of T cell activation" = "positive regulation of T cell activ.")) |> 
  filter(size > 10) |> 
  plotfun(title = "Inflammatory 1")

# inflammatory 2
infl2_goscatter <- bb_goscatter(simMatrix = leiden_34_summary_0.9$simMatrix,
             reducedTerms = leiden_34_summary_0.9$reducedTerms)[["data"]] |> 
  as_tibble() |> 
  filter(size > 10) |> 
  plotfun(title = "Inflammatory 2")

# activated bcr
bcr_goscatter <- bb_goscatter(simMatrix = leiden_other_summary_0.9$simMatrix,
             reducedTerms = leiden_other_summary_0.9$reducedTerms)[["data"]] |> 
  as_tibble() |> 
  mutate(parentTerm = recode(parentTerm, "immunoglobulin mediated immune response" = "Ig mediated immune response")) |> 
  filter(size > 10) |> 
  plotfun(title = "Activated BCR")
cds_main_leiden_comparison_tm |> filter(cell_group == "other") |> pull(gene_short_name)
