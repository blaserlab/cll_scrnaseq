source("00_packages_functions.R")

top_level_summary<-tbl_df(colData(cds_aligned_btk)) %>% 
  group_by(pt_timepoint) %>% 
  summarise(total_cells = n()) %>% 
  mutate(running_total = cumsum(total_cells)) %>%
  write_csv("data_out/top_level_summary.csv")

tbl_df(colData(cds_aligned_btk)) %>% 
  group_by(pt_timepoint,partition_assignment) %>% 
  summarise(total_cells = n(),
            n_genotyped_e10 = sum(!is.na(cell_call_e10pct)),
            n_genotyped_e1 = sum(!is.na(cell_call_e1pct)),
            mean_e10_wt_error = mean(max_bc_error_per_cell_wt_e10pct,na.rm = T),
            mean_e10_mt_error = mean(max_bc_error_per_cell_mt_e10pct,na.rm = T),
            mean_e1_wt_error = mean(max_bc_error_per_cell_wt_e1pct,na.rm = T),
            mean_e1_mt_error = mean(max_bc_error_per_cell_mt_e1pct,na.rm = T)) %>% 
  write_csv("data_out/btk_stats.csv")

#let genotyping rate (gr) = number of reads assigned to valid wt or mt barcodes/number of raw reads per sample

count_reads<-function(input,i) {
  result0<-as.numeric(system(paste0("wc -l ",matched_catted_fastqs[i]," | sed 's/ .*//g'"),intern = T))
  result<-result0/4
  return(result)
}

raw_read_counts<-lapply(X = seq_along(matched_catted_fastqs),
       FUN = count_reads,
       input = matched_catted_fastqs)
names(raw_read_counts)<-str_sub(matched_catted_fastqs,44,-7)

read_count_table<-tbl_df(raw_read_counts) %>% 
  pivot_longer(cols = starts_with("C"),names_to = "pt_timepoint",values_to = "read_count") %>% 
  mutate(pt_timepoint = str_to_lower(pt_timepoint))

cell_call_table<-tbl_df(colData(cds_aligned_btk)) %>% 
  group_by(pt_timepoint) %>% 
  summarise(total_cells = n(),
            btk_cells_1 = sum(!is.na(cell_call_e1pct)),
            btk_cells_10 = sum(!is.na(cell_call_e10pct)))

cell_genotype_efficiency<-left_join(cell_call_table, read_count_table) %>% 
  mutate(reads_per_cell_1 = read_count/btk_cells_1,
         reads_per_cell_10 = read_count/btk_cells_10)

cge_to_plot<-cell_genotype_efficiency %>% 
  select(pt_timepoint,reads_per_cell_10,reads_per_cell_1) %>% 
  pivot_longer(cols = starts_with("reads"), names_to = "threshold", values_to = "reads_per_cell_called") %>% 
  mutate(threshold = recode(threshold,
                            "reads_per_cell_10" = "10%",
                            "reads_per_cell_1" = "1%"))

genotyping_cost<-ggplot(data = cge_to_plot, aes(x = pt_timepoint, y = reads_per_cell_called, fill = threshold))+
  geom_bar(stat = "identity", position = "dodge", color = "black")+
  scale_fill_viridis_d(begin = 0.2, end = 0.8, alpha = 0.6)+
  labs(x = NULL,y = "Reads per Cell Genotyped",fill = "Threshold",title = "Genotyping Cost in Reads")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme(legend.position = c(0.1,0.8))+
  theme(plot.title = element_text(hjust = 0.5))
save_plot(genotyping_cost, filename = "plots_out/genotyping_cost.pdf", base_height = 3, base_width = 5)


save.image.pigz("cll_scrnaseq.RData",n.cores = 39)
system("cp cll_scrnaseq.RData ~/network/X/Labs/Blaser/single_cell/cll_project/R_session_data/cll_scrnaseq_mirror.RData")
