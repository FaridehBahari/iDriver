rm(list = ls())

library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(data.table)
library(doParallel)
library(foreach)

source('main/functions_iDriver.R')
source('prepareInput/4_Functions_responseTable.R')

load('../extdata/procInput/mut_PCAWG/all_samples.RData')
gr$varID_local <- 1:length(gr)

gr <- gr[, c("varID_local", "D_id")]


path_donorInfo <- '../extdata/procInput/iDriverInputs/donorInfo.tsv'
path_save <- '../extdata/procInput/BMRs/observed/baseLine_localRates/'
dir.create(path_save, recursive=T)

load('../bg_bed.RData')
bg_bed$PCAWG_ID <- bg_bed$name
bg_bed$name <- bg_bed$blockID
bg_bed <- bg_bed[, c("name", "PCAWG_ID", "length")]


df <- data.frame('PCAWG_ID' = bg_bed$PCAWG_ID,
                 'elem_ID' = bg_bed$name, 
                 'PCAWG_block_length' = bg_bed$length)

df <- df[!duplicated(df),]


create_response <- function(all_elements, n_donors){
  
  obs_Muts <- all_elements$MAF_GE_mapped
  df <- as.data.frame(DataFrame(obs_Muts))
  df <- df[,c('name', 'D_id')]
  
  # Convert df to a data.table (if it's not already one)
  setDT(df)
  
  # Create the result table using data.table syntax
  resTab <- df[, .N, by = .(name, D_id)][
    N != 0, .(nMut = sum(N), nSample = .N), by = .(binID = name)
    ]
  
  idx_lenElement <- resTab$binID
  all_lmnt_length <- all_elements$lenAllTests
  ElementLength <- all_lmnt_length[idx_lenElement]
  resTab$length <- ElementLength
  
  resTab$N <- n_donors
  
  all_binIds <- names(all_elements[["lenAllTests"]])
  nonMut_ids <- all_binIds[which(!all_binIds %in% resTab$binID)]
  res_nonMutated <- data.frame(cbind('binID' = nonMut_ids, 
                                     'nMut' = rep(0, length(nonMut_ids)),
                                     'nSample' = rep(0, length(nonMut_ids)),
                                     'length' = all_lmnt_length[!names(all_lmnt_length) %in% resTab$binID], 
                                     'N' = rep(n_donors, length(nonMut_ids))))
  
  resTab <- rbind(resTab, res_nonMutated)
  
  return(resTab)
}



cohorts <- c('Pancan-no-skin-melanoma-lymph', "Pan_Cancer")

for (cohort in cohorts) {
  
  print(cohort)
  
  exclude_lymph_melanoma <- ifelse(cohort == 'Pancan-no-skin-melanoma-lymph', T, F)
  
  donorInfo <- select_cohort(path_donorInfo,
                             cohort,
                             exclude_lymph_melanoma,
                             exclude_hyper_mutated = T)
  
  gr_cohort <- gr[which(gr$D_id %in% donorInfo$donor_id)]
  n_donors <- length(unique(gr_cohort$D_id))
  
  
  elem_types <- c("gc19_pc.cds", "enhancers", "gc19_pc.3utr",
                  "gc19_pc.5utr", "gc19_pc.promCore")
  
  # dir.create('../extdata/procInput/tmp_localBMR/pancanStar/', recursive = T)
  pred_tab <- c()
  
  for (i in 1:length(elem_types)) {
    elem <- elem_types[i]
    print(elem)
    
    df_elem <- df[grepl(elem, df$PCAWG_ID), ]
    
    bg_elem <- bg_bed[grepl(elem,bg_bed$PCAWG_ID)]
    
    all_elements <- map_muts(gr_cohort, bg_elem, is_gr_tetsGE = T)
    print('mapped')
    
    dir.create(paste0('../extdata/procInput/tmp_localBMR/', cohort, '/'), recursive = T)
    save(all_elements, file = paste0('../extdata/procInput/tmp_localBMR/', cohort, '/', elem, '_all_elements.RData'))
    
    fwrite(df_elem, file = paste0('../extdata/procInput/tmp_localBMR/', cohort, '/', elem, '_df_elem.tsv'), sep = '\t')
    
    # resTab_elem <- create_response(all_elements, n_donors)
    # print('res elem created')
    # pred_tab_elem <- left_join(resTab_elem, df_elem, by = c('binID'= 'elem_ID'))
    # pred_tab_elem <- pred_tab_elem[!duplicated(pred_tab_elem),]
    # print('duplicates removed')
    # pred_tab <- rbind(pred_tab, pred_tab_elem)
  }
  
  
  # fwrite(pred_tab, file = paste0(path_save,'/baselineBMR_', cohort, '.tsv'), sep = '\t')
}


# for (cohort in cohorts) {
#   
#   print(cohort)
#   
#   exclude_lymph_melanoma <- ifelse(cohort == 'Pancan-no-skin-melanoma-lymph', T, F)
#   
#   donorInfo <- select_cohort(path_donorInfo,
#                              cohort,
#                              exclude_lymph_melanoma,
#                              exclude_hyper_mutated = T)
#   
#   gr_cohort <- gr[which(gr$D_id %in% donorInfo$donor_id)]
#   n_donors <- length(unique(gr_cohort$D_id))
#   
#   
#   elem_types <- c("gc19_pc.cds", "enhancers", "gc19_pc.3utr",
#                   "gc19_pc.5utr", "gc19_pc.promCore")
#   
#   pred_tab <- c()
#   
#   for (i in 1:length(elem_types)) {
#     elem <- elem_types[i]
#     print(elem)
#     
#     df_elem <- df[grepl(elem, df$PCAWG_ID), ]
#     
#     bg_elem <- bg_bed[grepl(elem,bg_bed$PCAWG_ID)]
#     
#     all_elements <- map_muts(gr_cohort, bg_elem, is_gr_tetsGE = T)
#     print('mapped')
#     resTab_elem <- create_response(all_elements, n_donors)
#     print('res elem created')
#     pred_tab_elem <- left_join(resTab_elem, df_elem, by = c('binID'= 'elem_ID'))
#     pred_tab_elem <- pred_tab_elem[!duplicated(pred_tab_elem),]
#     print('duplicates removed')
#     pred_tab <- rbind(pred_tab, pred_tab_elem)
#   }
#   
#   
#   fwrite(pred_tab, file = paste0(path_save,'/baselineBMR_', cohort, '.tsv'), sep = '\t')
# }
