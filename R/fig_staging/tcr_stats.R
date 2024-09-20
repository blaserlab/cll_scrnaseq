bb_cellmeta(cds_main[,colData(cds_main)$partition_assignment %in% c("T")]) |> 
  mutate(tcr_true = ifelse(!is.na(tcr_cdr3s_aa), TRUE, FALSE)) |> 
  count(tcr_true) |> 
  pivot_wider(names_from = tcr_true, values_from = n) |> 
  mutate(total = `FALSE` + `TRUE`,
         pct_true = `TRUE`/total*100)
