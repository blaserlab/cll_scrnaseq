source("00_packages_functions.R")

# generate new composite factor columns pt_timepoint and cb_pt_timepoint
colData(cds_aligned)$pt_timepoint <-
  paste0(colData(cds_aligned)$pt, "_", colData(cds_aligned)$timepoint)

colData(cds_aligned)$cb_pt_timepoint <-
  paste0(colData(cds_aligned)$barcode,
         "_",
         colData(cds_aligned)$pt_timepoint)
# make fasta from cds for all valid barcodes and prepend the illumina handle
ins<-unique(colData(cds_aligned)$pt_timepoint)

if (!dir.exists("working_files/cb_fastas")) {
  dir.create("working_files/cb_fastas")
}

for (i in 1:length(ins)) {
  cell_names <-
    tbl_df(colData(cds_aligned)) %>%
    filter(pt_timepoint == ins[i]) %>%
    pull(cb_pt_timepoint)
  cb <-
    tbl_df(colData(cds_aligned)) %>%
    filter(pt_timepoint == ins[i]) %>%
    pull(barcode)
  sink(paste0(
    "working_files/cb_fastas/",
    ins[i],
    ".fasta"
  ))
  cat(paste0(c(
    rbind(
      ">",
      cell_names,
      "\nCTACACGACGCTCTTCCGATCT",
      as.character(str_sub(cb, 1,-3)),
      "\n"
    )
  ), collapse = ""))
  sink()
}

save.image.pigz("cll_scrnaseq.RData",n.cores = 39)
