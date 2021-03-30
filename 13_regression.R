source("00_packages_functions.R")


priya_interesting_genes <- c(
  "BHLHE41",
  "DGKG",
  "DDX17",
  "DENND3",
  "ZDHHC19",
  "RALGPS2",
  "PDE4D",
  "MALAT1",
  "GRB2",
  "FLNB"
  
)


colData(cds_final)$timepoint_pretty <- paste0("Timepoint ", colData(cds_final)$timepoint)

cds_regression_subset <- cds_final[rowData(cds_final)$gene_short_name %in% priya_interesting_genes, 
                                   colData(cds_final)$predicted.celltype.l1 == "B"]


# make the basic timepoint model
timepoint_model <- fit_models(cds_regression_subset, model_formula_str = "~timepoint_pretty",expression_family = "negbinomial")
timepoint_model_coefs <- coefficient_table(timepoint_model)
timepoint_model_coefs %>%
  filter(term != "(Intercept)") %>%
  select(gene_short_name, term, q_value, estimate)

# add in pt variable to the model
timepoint_pt_model <- fit_models(cds_regression_subset, model_formula_str = "~timepoint_pretty+pt",expression_family = "negbinomial")
timepoint_pt_model_coefs <- coefficient_table(timepoint_pt_model)

timepoint_pt_model_coefs %>%
  filter(term != "(Intercept)") %>%
  select(gene_short_name, term, q_value, estimate) %>%
  write_csv("data_out/selected_gene_regression.csv")



compare_models(model_tbl_reduced = timepoint_model, model_tbl_full = timepoint_pt_model) %>% select(gene_short_name, q_value)


