source("00_packages_functions.R")

# first split by CNNN or not
samples<-list.files("working_files/btk_fastqs/matched_run1_run2",full.names = T)
sample_names<-str_sub(list.files("working_files/btk_fastqs/matched_run1_run2",full.names = F),1,-7)

if (!dir.exists("working_files/btk_fastqs/wt_fastq")) {
  dir.create("working_files/btk_fastqs/wt_fastq")
}

if (!dir.exists("working_files/btk_fastqs/untrimmed1")) {
  dir.create("working_files/btk_fastqs/untrimmed1")
}

filter_wt_fastqs <- function(infiles, names, i) {
  cmd <-
    paste0(
      "cutadapt -g ^C -e 0.3 --action=none -o working_files/btk_fastqs/wt_fastq/",
      names[i],
      "_wt.fastq ",
      infiles[i],
      " --untrimmed-output working_files/btk_fastqs/untrimmed1/",
      names[i],
      "_untrimmed1.fastq"
    )
  message(cmd, "\n")
  system(cmd)
}

lapply(X = seq_along(samples),
       FUN = filter_wt_fastqs,
       infiles = samples,
       names = sample_names)

#load in the filenames for untrimmed fastqs 
untrimmed1_in_fp <-
  list.files("working_files/btk_fastqs/untrimmed1", full.names = T)
untrimmed_names <-
  str_sub(list.files("working_files/btk_fastqs/untrimmed1", full.names = F),
          1,
          -18)

# prepare a directory for the mutant fastqs and still untrimmed
if (!dir.exists("working_files/btk_fastqs/mt_fastq/")) {
  dir.create("working_files/btk_fastqs/mt_fastq/")
}

if (!dir.exists("working_files/btk_fastqs/untrimmed2/")) {
  dir.create("working_files/btk_fastqs/untrimmed2")
}

# function to filter mutant fastqs
filter_mt_fastqs <- function(infiles, names, i) {
  cmd <-
    paste0(
      "cutadapt -g ^A -e 0.3 --action=none -o working_files/btk_fastqs/mt_fastq/",
      names[i],
      "_mt.fastq ",
      infiles[i],
      " --untrimmed-output working_files/btk_fastqs/untrimmed2/",
      names[i],
      "_untrimmed2.fastq"
    )
  message(cmd, "\n")
  system(cmd)
}

lapply(X = seq_along(samples),
       FUN = filter_mt_fastqs,
       infiles = untrimmed1_in_fp,
       names = untrimmed_names)

# delete fastqs still untrimmed
unlink("working_files/btk_fastqs/untrimmed1",recursive = T)
unlink("working_files/btk_fastqs/untrimmed2",recursive = T)

#generate info files with cell barcodes
mt_fastq_in<-list.files("working_files/btk_fastqs/mt_fastq")#demulitplexed and btk-mutant fastq files
wt_fastq_in<-list.files("working_files/btk_fastqs/wt_fastq")#demultiplexed and btk-wt fastq files
fasta_in<-list.files("working_files/cb_fastas")#fasta barcode files

mt_fastq_in_fp<-list.files("working_files/btk_fastqs/mt_fastq", full.names = T)
wt_fastq_in_fp<-list.files("working_files/btk_fastqs/wt_fastq", full.names = T)
fasta_in_fp<-list.files("working_files/cb_fastas", full.names = T)

assert_that(sum(str_to_upper(str_sub(fasta_in, 1, -7)) == str_to_upper(str_sub(mt_fastq_in, 1, -10)))==12)
assert_that(sum(str_to_upper(str_sub(fasta_in, 1, -7)) == str_to_upper(str_sub(wt_fastq_in, 1, -10)))==12)

names_for_infofile<-str_sub(fasta_in,1,-7)

#prepare a directory for mutant and wt fastqs, even though we don't really use them
if (!dir.exists("working_files/btk_fastqs/mt_fastq_cb_trimmed")) {
  dir.create("working_files/btk_fastqs/mt_fastq_cb_trimmed")
}

if (!dir.exists("working_files/btk_fastqs/wt_fastq_cb_trimmed")) {
  dir.create("working_files/btk_fastqs/wt_fastq_cb_trimmed")
}

#now prepare directories for the infofiles
if (!dir.exists("working_files/mt_info_files")) {
  dir.create("working_files/mt_info_files")
}

if (!dir.exists("working_files/wt_info_files")) {
  dir.create("working_files/wt_info_files")
}

#match reads to list of [illumina read 1][10x cell barcode whitelist from all cells]
mt_match_func <- function(sample_fps,sample_names,fastas,i) {
  cmd <-
    paste0(
      "cutadapt -g file:",
      fastas[i],
      " --discard-untrimmed -e 0.01 --overlap 38 --action=trim -o working_files/btk_fastqs/mt_fastq_cb_trimmed/",
      sample_names[i],
      "_mt_cb_trimmed.fastq --info-file working_files/mt_info_files/",
      sample_names[i],
      "_mt_info.txt ",
      sample_fps[i]
    )
  message(cmd, "\n")
  system(cmd)
}

mclapply(X = seq_along(mt_fastq_in_fp),
         FUN = mt_match_func,
         sample_fps = mt_fastq_in_fp,
         sample_names = names_for_infofile,
         fastas = fasta_in_fp,
         mc.cores = 12,
         mc.preschedule = T)



wt_match_func <- function(sample_fps,sample_names,fastas,i) {
  cmd <-
    paste0(
      "cutadapt -g file:",
      fastas[i],
      " --discard-untrimmed -e 0.01 --overlap 38 --action=trim -o working_files/btk_fastqs/wt_fastq_cb_trimmed/",
      sample_names[i],
      "_wt_cb_trimmed.fastq --info-file working_files/wt_info_files/",
      sample_names[i],
      "_wt_info.txt ",
      sample_fps[i]
    )
  message(cmd, "\n")
  system(cmd)
}

mclapply(X = seq_along(wt_fastq_in_fp),
         FUN = wt_match_func,
         sample_fps = wt_fastq_in_fp,
         sample_names = names_for_infofile,
         fastas = fasta_in_fp,
         mc.cores = 12,
         mc.preschedule = T)


#load the cutadapt info files
wt_info_fp<-list.files("working_files/wt_info_files", full.names = T)
mt_info_fp<-list.files("working_files/mt_info_files", full.names = T)

#summarize each info file to include only cell barcodes, cell names, mutant reads or wt reads collapsed by umi

summarize_cutadapt_byumi <- function(info_files, i) {
  info <-
    tbl_df(read.delim(
      file = info_files[i],
      header = F,
      sep = "\t",
      col.names = as.character(
        c("V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11")
      )
    )) %>% filter(V2 != -1) %>% select(
      bc_errors = V2,
      matched_seq = V6,
      restofread = V7,
      cell_names = V8
    )
  info_umi <- info %>%
    mutate(cell_names_umi = paste0(cell_names, "_", str_sub(restofread, 1, 10))) %>%
    group_by(cell_names_umi) %>%
    summarise(max_bc_error_per_umi = max(bc_errors)) %>%
    mutate(cell_collapsed = str_sub(cell_names_umi, 1, -12)) %>%
    ungroup() %>%
    group_by(cell_collapsed) %>%
    summarise(
      umi_per_cell = n(),
      max_bc_error_per_cell = max(max_bc_error_per_umi)
    )
  return(info_umi)

}
#all cells
wt_info_byumi<-lapply(X = seq_along(wt_info_fp),
                      FUN = summarize_cutadapt_byumi,
                      info_files = wt_info_fp)
names(wt_info_byumi)<-str_sub(wt_info_fp,29,-10);wt_info_byumi

mt_info_byumi<-lapply(X = seq_along(mt_info_fp),
                      FUN = summarize_cutadapt_byumi,
                      info_files = mt_info_fp)
names(mt_info_byumi)<-str_sub(mt_info_fp,29,-10);mt_info_byumi

# get wt and mutant reads back together for each sample UMIS
join_wt_and_mut <- function(i) {
  joined_df <-
    full_join(wt_info_byumi[[i]],
              mt_info_byumi[[i]],
              by = "cell_collapsed",
              suffix = c("_wt", "_mt"))
  #return(joined_df)
  return(
    joined_df %>% replace_na(list(
      umi_per_cell_wt = 0, umi_per_cell_mt = 0
    )) %>% mutate(
      mut_to_wt = umi_per_cell_mt / umi_per_cell_wt,
      cell_call = ifelse(mut_to_wt > 0, yes = "mutant", no = "wt")# if any mutant reads are there it is mutant
    )
  )
}

joined_info <-
  lapply(X = 1:length(wt_info_byumi), FUN = join_wt_and_mut)

# now join all back together
cell_calls_e1pct<-bind_rows(joined_info) %>% dplyr::rename(cb_pt_timepoint = cell_collapsed)

save.image.pigz("cll_scrnaseq.RData",n.cores = 39)
