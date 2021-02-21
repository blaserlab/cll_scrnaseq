source("00_packages_functions.R")

analysis_configs <- read_excel("~/network/X/Labs/Blaser/single_cell/cll_project/analysis_configs.xlsx",sheet = "ingest") %>% 
  select(-parental_folder, -`note for priya`) %>% 
  mutate(new_path = str_replace_all(`10X_dir`,"\\\\","/")) %>%
  mutate(new_path = str_replace(new_path, "X:", "~/network/X")) %>%
  select(patient, timepoint, lib_type, directory = new_path) %>%
  mutate(patient = paste0("cll_",patient)) %>%
  mutate(timepoint = as_factor(timepoint)) %>%
  mutate(pipestance_names = paste0("cds_",patient,"_",timepoint)) %>%
  left_join(read_excel("~/network/X/Labs/Blaser/single_cell/cll_project/analysis_configs.xlsx",sheet = "metadata",
                       col_types = c("text", "text", "date", "date", "numeric", "date", "numeric", "text","numeric", "text", "text", "logical", "logical", "logical", "logical")) %>%
              mutate(patient = paste0("cll_",patient)))
  

cds_list <-
  map(
    .x = analysis_configs %>% filter(lib_type == "gex") %>% pull(pipestance_names),
    .f = function(x, data = analysis_configs %>% filter(lib_type == "gex")) {
      directory <- data %>% filter(pipestance_names == x) %>% pull(directory)
      cds <- load_cellranger_data(directory)
      cds <-
        add_cds_factor_columns(
          cds = cds,
          columns_to_add = c(
            "pt" = data %>% filter(pipestance_names == x) %>% pull(patient),
            "timepoint" = data %>% filter(pipestance_names == x) %>% pull(timepoint),
            "pipestance" = data %>% filter(pipestance_names == x) %>% pull(pipestance_names),
            "ID" = data %>% filter(pipestance_names == x) %>% pull(ID),
            "baseline_date" = data %>% filter(pipestance_names == x) %>% pull(baseline_date),
            "btk_clone_date" = data %>% filter(pipestance_names == x) %>% pull(btk_clone_date),
            "btk_clone_vaf_pct" = data %>% filter(pipestance_names == x) %>% pull(btk_clone_vaf_pct),
            "relapse_date" = data %>% filter(pipestance_names == x) %>% pull(relapse_date),
            "relapse_vaf_pct" = data %>% filter(pipestance_names == x) %>% pull(relapse_vaf_pct),
            "gender" = data %>% filter(pipestance_names == x) %>% pull(gender),
            "age_at_ibr_start" = data %>% filter(pipestance_names == x) %>% pull(age_at_ibr_start),
            "IGHV_status" = data %>% filter(pipestance_names == x) %>% pull(IGHV_status),
            "complex_karyotype" = data %>% filter(pipestance_names == x) %>% pull(complex_karyotype),
            "del17p" = data %>% filter(pipestance_names == x) %>% pull(del17p),
            "del11q" = data %>% filter(pipestance_names == x) %>% pull(del11q),
            "del13p" = data %>% filter(pipestance_names == x) %>% pull(del13p),
            "tri12" = data %>% filter(pipestance_names == x) %>% pull(tri12)
          )
        )
      return(cds)
    }
  ) %>% set_names(analysis_configs %>% filter(lib_type == "gex") %>% pull(pipestance_names))


summarized_sequencing_metrics <-
  tibble(pipestance_path = analysis_configs %>% filter(lib_type == "gex") %>% pull(directory)) %>%
  mutate(metrics_summary_path = paste0(pipestance_path, "/outs/metrics_summary.csv")) %>%
  mutate(cds_dim_cells = unname(sapply(X = cds_list, FUN = dim)[2, ])) %>%
  mutate(cds_name = names(sapply(X = cds_list, FUN = dim)[2, ])) %>%
  left_join(.,
            bind_rows(lapply(
              X = .$metrics_summary_path, FUN = read_csv
            )),
            by = c("cds_dim_cells" = "Estimated Number of Cells")) %>% # this is your sanity check.  Joining on cell number derived from the CDS object and the metrics summary
  select(
    cds_name,
    cds_dim_cells,
    `Mean Reads per Cell`,
    `Median Genes per Cell`,
    `Fraction Reads in Cells`
  ) %>%
  write_csv("data_out/summarized_sequencing_metrics.csv")

