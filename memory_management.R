cached_items <- c(
  "cds_list_rejoined",
  "cds_list",
  "cds_clean",
  "cds"
)

save.pigz(list = cached_items, file = "rdata/cached_items.RData", n.cores = 8)
rm(list = cached_items)


save.image.pigz("rdata/active_workspace.RData",n.cores = 8)
