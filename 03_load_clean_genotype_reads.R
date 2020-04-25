source("00_packages_functions.R")

# rename and load extracted fastqs from tom
## identify source data
btkv1_fp <-
  list.files("~/network/X/Labs/Blaser/single_cell/cll_project/BTK_nanopore_v1", full.names = T)

btkv2_fp <-
  list.files("~/network/X/Labs/Blaser/single_cell/cll_project/2020-03-25-BTK-capture", full.names = T)

## make working directories
if (!dir.exists("working_files/btk_fastqs")) {
  dir.create("working_files/btk_fastqs")
}

if (!dir.exists("working_files/btk_fastqs/fastqs_run1")) {
  dir.create("working_files/btk_fastqs/fastqs_run1")
}

if (!dir.exists("working_files/btk_fastqs/fastqs_run2")) {
  dir.create("working_files/btk_fastqs/fastqs_run2")
}

## fetch the source files to R working directory
fastqs_to_fetch <- c(btkv1_fp, btkv2_fp)
newnames <-
  c(
    "working_files/btk_fastqs/fastqs_run1/CLL5_baseline.fastq",
    "working_files/btk_fastqs/fastqs_run1/CLL5_btk_clone.fastq",
    "working_files/btk_fastqs/fastqs_run1/CLL5_relapse.fastq",
    "working_files/btk_fastqs/fastqs_run1/CLL6_baseline.fastq",
    "working_files/btk_fastqs/fastqs_run1/CLL6_btk_clone.fastq",
    "working_files/btk_fastqs/fastqs_run1/CLL6_relapse.fastq",
    "working_files/btk_fastqs/fastqs_run1/CLL7_baseline.fastq",
    "working_files/btk_fastqs/fastqs_run1/CLL7_btk_clone.fastq",
    "working_files/btk_fastqs/fastqs_run1/CLL7_relapse.fastq",
    "working_files/btk_fastqs/fastqs_run2/CLL5_baseline.fastq",
    "working_files/btk_fastqs/fastqs_run2/CLL5_btk_clone.fastq",
    "working_files/btk_fastqs/fastqs_run2/CLL5_relapse.fastq",
    "working_files/btk_fastqs/fastqs_run2/CLL6_baseline.fastq",
    "working_files/btk_fastqs/fastqs_run2/CLL6_btk_clone.fastq",
    "working_files/btk_fastqs/fastqs_run2/CLL6_relapse.fastq",
    "working_files/btk_fastqs/fastqs_run2/CLL7_baseline.fastq",
    "working_files/btk_fastqs/fastqs_run2/CLL7_btk_clone.fastq",
    "working_files/btk_fastqs/fastqs_run2/CLL7_relapse.fastq",
    "working_files/btk_fastqs/fastqs_run2/CLL8_baseline.fastq",
    "working_files/btk_fastqs/fastqs_run2/CLL8_btk_clone.fastq",
    "working_files/btk_fastqs/fastqs_run2/CLL8_relapse.fastq"
  )

fetch_and_rename <- function(infiles, outfiles, i) {
  cmd<-paste0("cp ", infiles[i], " ", outfiles[i])
  message(cmd, "\n")
  system(cmd)
}

lapply(X =seq_along(fastqs_to_fetch),
       FUN = fetch_and_rename,
       infiles = fastqs_to_fetch,
       outfiles = newnames)

# fix the fastqs from btkv1 to turn > into @
btkv1_working<-list.files("working_files/btk_fastqs/fastqs_run1", full.names = T)
btkv2_working<-list.files("working_files/btk_fastqs/fastqs_run2", full.names = T)


fix_fastqs<-function(infiles,i) {
  cmd<-paste0("sed -i '/^>/ s/>/@/' ",infiles[i])
  message(cmd, "\n")
  system(cmd)
}

lapply(X = seq_along(btkv1_working),
       FUN = fix_fastqs,
       infiles = btkv1_working)

# rename the reads so they are specific by appending _v1 or _v2
rename_reads<-function(infiles,version,i) {
  cmd<-paste0("sed -i '/^@/ s/$/_",version,"/' ",infiles[i])
  message(cmd, "\n")
  system(cmd)
}

lapply(X= seq_along(btkv1_working),
       FUN = rename_reads,
       infiles = btkv1_working,
       version = "v1")

lapply(X= seq_along(btkv2_working),
       FUN = rename_reads,
       infiles = btkv2_working,
       version = "v2")


# join the fastqs together
match_and_cat <- function(infiles1, infiles2, outfiles, i) {
  if (length(infiles1) > length(infiles2)) {
    infiles2 <-
      c(infiles2, rep("", times = length(infiles1) - length(infiles2)))
  } else if (length(infiles1) < length(infiles2)) {
    infiles1 <-
      c(infiles1, rep("", times = length(infiles2) - length(infiles1)))
  }
  
  if (str_sub(infiles1[i], 38, -1) == str_sub(infiles2[i], 38, -1)) {
    cmd<-paste0("cat ", infiles1[i]," ",infiles2[i]," > ", outfiles[i])
    message(cmd, "\n")
    system(cmd)
  } else if (infiles1[i]=="") {
    cmd<-paste0("cat ", infiles2[i]," > ", outfiles[i])
    message(cmd, "\n")
    system(cmd)
  } else if (infiles2[i]=="") {
    cmd<-paste0("cat ", infiles1[i]," > ", outfiles[i])
    message(cmd, "\n")
    system(cmd)
  } else {
    print("There was an error:  the files didn't match.")
  }
  
}

if (!dir.exists("working_files/btk_fastqs/matched_run1_run2")) {
  dir.create("working_files/btk_fastqs/matched_run1_run2")
}

matched_catted_fastqs<-c("working_files/btk_fastqs/matched_run1_run2/CLL5_baseline.fastq",
                        "working_files/btk_fastqs/matched_run1_run2/CLL5_btk_clone.fastq",
                        "working_files/btk_fastqs/matched_run1_run2/CLL5_relapse.fastq",
                        "working_files/btk_fastqs/matched_run1_run2/CLL6_baseline.fastq",
                        "working_files/btk_fastqs/matched_run1_run2/CLL6_btk_clone.fastq",
                        "working_files/btk_fastqs/matched_run1_run2/CLL6_relapse.fastq",
                        "working_files/btk_fastqs/matched_run1_run2/CLL7_baseline.fastq",
                        "working_files/btk_fastqs/matched_run1_run2/CLL7_btk_clone.fastq",
                        "working_files/btk_fastqs/matched_run1_run2/CLL7_relapse.fastq",
                        "working_files/btk_fastqs/matched_run1_run2/CLL8_baseline.fastq",
                        "working_files/btk_fastqs/matched_run1_run2/CLL8_btk_clone.fastq",
                        "working_files/btk_fastqs/matched_run1_run2/CLL8_relapse.fastq")

lapply(X = seq_along(matched_catted_fastqs),
       FUN = match_and_cat,
       infiles1 = btkv1_working,
       infiles2 = btkv2_working,
       outfiles = matched_catted_fastqs)

save.image.pigz("cll_scrnaseq.RData",n.cores = 39)
