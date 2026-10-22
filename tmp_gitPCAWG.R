
  # Libraries
 
library(data.table)
library(tidyverse)
library(readxl)
library(ggVennDiagram)
library(maftools)

# Prepare MAF data

pcawg_df <- data.table::fread("../extdata/procInput0/mut_PCAWG_beforeRevision/Don_Sam_Proj.tsv")





# Look at all the different Project Codes
pcawg_df %>%
  dplyr::count(Project_Code, sort = TRUE)

# Split df into a list of dfs, one for each Project_Code
pcawg_split_ls = split(pcawg_df, f = pcawg_df$Project_Code)



# Prepare Metadata

# Read data
hist_df <- read_excel(path = '../extdata/rawInput/PCAWG_supp/pcawg_specimen_histology_August2016_v9.xlsx')
clin_df <- read_excel(path = '../extdata/rawInput/PCAWG_supp/pcawg_donor_clinical_August2016_v9.xlsx')
sample_df <- fread('../extdata/rawInput/pcawg_sample_sheet.tsv')

sample_level_metadata <- full_join(hist_df, clin_df, by = "icgc_donor_id")
sample_level_metadata <- full_join(sample_df, sample_level_metadata, by = "icgc_sample_id") %>% 
  filter(!is.na(aliquot_id))

if(anyDuplicated(sample_level_metadata$aliquot_id) > 0)
  message("PROBLEM WITH METADATA, DUPLICATED ALIQUOT ID")

sample_level_metadata <- sample_level_metadata %>% 
  dplyr::rename(Tumor_Sample_Barcode = aliquot_id)


S_D_pcawg <- sample_level_metadata[,c('icgc_donor_id.x', 'Tumor_Sample_Barcode', 'library_strategy')]
S_D_pcawg <- S_D_pcawg[!duplicated(S_D_pcawg),]


sample_level_metadata1 <- sample_level_metadata %>% filter(Tumor_Sample_Barcode %in% pcawg_df$Tumor_Sample_Barcode)


