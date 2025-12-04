library(tidyverse)
library(readxl)
library(blaseRtools)
library(conflicted)
conflicts_prefer(dplyr::filter)

blaseRtemplates::project_data(path = "~/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/datapkg")

geo_full_tbl <- read_excel("~/network/X/Labs/Blaser/staff/single_cell/cll_project/analysis_configs.xlsx", sheet = "ingest") |> 
  select(patient, timepoint, lib_type, specimen_date, pipestance_path = `10X_dir`) |> 
  mutate(patient = as.integer(patient)) |> 
  mutate(pipestance_path = bb_fix_file_path(pipestance_path)) |> 
  mutate(pipestance_path = str_replace(pipestance_path, "Blaser/", "Blaser/staff/")) |> 
  mutate(pipestance_path = fs::path(pipestance_path)) 

samples <- bb_cellmeta(cds_main) |> 
  group_by(specimen, patient, patient_type3, timepoint_merged_1) |> 
  summarise() |> 
  arrange(patient, timepoint_merged_1) |> 
  mutate(gex = "gex", bcr = "bcr", tcr = "tcr") |> 
  pivot_longer(cols = c(gex, bcr, tcr), names_to = "library_type") |> 
  select(-value) |> 
  mutate(title = paste0(patient, ", ",patient_type3, ", timepoint ",timepoint_merged_1, ", ", library_type, " library")) |> 
  relocate(title, .after = specimen) |> 
  relocate(library_type, .after = specimen) |>
  mutate(library_name = paste0(specimen, "_", library_type)) |> 
  relocate(library_name, .after = specimen)
write_csv(samples, "~/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/geo/samples.csv")


# transfer over the processed data
processed_gex_files <- geo_full_tbl |> 
  filter(lib_type == "gex") |> 
  mutate(barcodes = fs::path(pipestance_path, "outs", "filtered_feature_bc_matrix", "barcodes.tsv.gz")) |> 
  mutate(features = fs::path(pipestance_path, "outs", "filtered_feature_bc_matrix", "features.tsv.gz")) |> 
  mutate(matrix = fs::path(pipestance_path, "outs", "filtered_feature_bc_matrix", "matrix.mtx.gz")) |> 
  pivot_longer(c(barcodes, features, matrix), names_to = "component", values_to = "source_path")


processed_bcr_tcr_files <- geo_full_tbl |> 
  filter(lib_type %in% c("bcr", "tcr")) |> 
  mutate(all_contig_annotations = fs::path(pipestance_path, "outs", "all_contig_annotations.csv")) |> 
  select(-pipestance_path) |> 
  mutate(dest_path =fs::path("~/brad_workspace/sc_working/cll_upload/processed_1", paste0("cll_", patient, "_", timepoint, "_", lib_type, "_all_contig_annotations.csv")))
processed_bcr_tcr_files

# copy over the processed bcr_tcr files
fs::file_copy(path = processed_bcr_tcr_files$all_contig_annotations, new_path = processed_bcr_tcr_files$dest_path)

samples
processed <- tibble(path = fs::path(c(fs::dir_ls("~/brad_workspace/sc_working/cll_upload/processed_1"),
  fs::dir_ls("~/brad_workspace/sc_working/cll_upload/geo_submission_cll_scrnaseq/")))) |>
  mutate(file_type = case_when(str_detect(path, "fastq.gz") ~ "raw", .default = "processed")) |> 
  filter(file_type == "processed") |> 
  mutate(file_name = fs::path_file(path)) |> 
  mutate(patient = str_remove(file_name, "cll_|pt_")) |> 
  mutate(patient = str_extract(patient, "[:alnum:]*|")) |> 
  mutate(patient = str_extract(patient, "[:digit:]*")) |> 
  mutate(timepoint = case_when(str_detect(file_name, "1DX") ~ "1",
                               str_detect(file_name, "1RL1") ~ "3",
                               str_detect(file_name, "1RM") ~ "2",
                               str_detect(file_name, "3yrs") ~ "2",
                               str_detect(file_name, "5yrs") ~ "3",
                               str_detect(file_name, "baseline") ~ "1",
                               str_detect(file_name, "relapse") ~ "3",
                               str_detect(file_name, "btk_clone") ~ "2")) |>
  mutate(name_timepoint = str_extract(file_name, "baseline|3yrs|5yrs|btk_clone|relapse")) |> 
  mutate(patient = factor(patient)) |> 
  mutate(timepoint_merged_1 = factor(timepoint)) |> 
  select(patient, timepoint_merged_1, file_name, name_timepoint) |>
  mutate(library_type = case_when(str_detect(file_name, "bcr") ~ "bcr",
                                  str_detect(file_name, "tcr") ~ "tcr",
                                  .default = "gex")) |> 
  mutate(library_name = paste0("cll_", patient, "_", name_timepoint, "_", library_type))  |> 
  group_by(library_name) |> 
  mutate(rn = row_number()) |>
  select(library_name, file_name, rn) |> 
  pivot_wider(names_from = rn, values_from = file_name, values_fill = "")
  
waldo::compare(
  samples$library_name,
left_join(samples, processed) |> pull(library_name))

processed$library_name[!processed$library_name %in% samples$library_name]
left_join(samples, processed) |> 
  select(c(`1`, `2`, `3`)) |> 
  write_csv(file = "~/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/geo/processed_file_by_sample.csv")

raw <- tibble(path = fs::dir_ls("~/brad_workspace/sc_working/cll_upload/geo_submission_cll_scrnaseq/")) |>
  mutate(file_type = case_when(str_detect(path, "fastq.gz") ~ "raw", .default = "processed")) |> 
  filter(file_type == "raw") |> 
  mutate(file_name = fs::path_file(path)) |> 
  mutate(patient = str_remove(file_name, "cll_|pt_")) |> 
  mutate(patient = str_extract(patient, "[:alnum:]*|")) |> 
  mutate(patient = str_extract(patient, "[:digit:]*")) |> 
  mutate(timepoint = case_when(str_detect(file_name, "1DX") ~ "1",
                               str_detect(file_name, "1RL1") ~ "3",
                               str_detect(file_name, "1RM") ~ "2",
                               str_detect(file_name, "_1_") ~ "1",
                               str_detect(file_name, "_2_") ~ "2",
                               str_detect(file_name, "_3_") ~ "3",
                               str_detect(file_name, "3yrs") ~ "2",
                               str_detect(file_name, "5yrs") ~ "3",
                               str_detect(file_name, "baseline") ~ "1",
                               str_detect(file_name, "relapse") ~ "3",
                               str_detect(file_name, "btk_clone") ~ "2")) |> 
  mutate(name_timepoint = case_match(timepoint, "1" ~ "baseline", .default = timepoint)) |> 
  mutate(name_timepoint = case_when(patient %in% c("1", "2", "3", "4") & timepoint == "2" ~ "3yrs", 
                               patient %in% c("1", "2", "3", "4") & timepoint == "3" ~ "5yrs",
                               timepoint == "2" ~ "btk_clone",
                               timepoint == "3" ~ "relapse",
                               .default = name_timepoint)) |> 
  mutate(patient = factor(patient)) |> 
  mutate(timepoint_merged_1 = factor(timepoint)) |> 
  select(patient, timepoint_merged_1, file_name, name_timepoint) |>
  mutate(library_type = case_when(str_detect(file_name, "bcr|BCR|Bcell|B_cell|B_Cell") ~ "bcr",
                                  str_detect(file_name, "tcr|TCR|Tcell|T_cell|T_Cell") ~ "tcr",
                                  .default = "gex")) |> 
  mutate(library_name = paste0("cll_", patient, "_", name_timepoint, "_", library_type))  |> 
  group_by(library_name) |> 
  mutate(rn = row_number()) |>
  select(library_name, file_name, rn) |> 
  pivot_wider(names_from = rn, values_from = file_name, values_fill = "") 

waldo::compare(
  samples$library_name,
  left_join(samples, raw) |> pull(library_name))

raw$library_name[!raw$library_name %in% samples$library_name]



left_join(samples, raw) |> write_csv(file = "~/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/geo/raw_file_by_sample.csv")



#md5 sums
fs::dir_ls("~/brad_workspace/sc_working/cll_upload/geo_submission_cll_scrnaseq_1") |> digest::digest()

md5sums <- map(.x = fs::dir_ls("~/brad_workspace/sc_working/cll_upload/geo_submission_cll_scrnaseq_1"), .f = \(x) {
  tibble(file = fs::path_file(x), md5sum = digest::digest(x, file = TRUE))
}) |> bind_rows()

write_csv(md5sums, file = "~/network/X/Labs/Blaser/share/collaborators/cll_scrnaseq_manuscript/geo/md5sums_additional.csv")
