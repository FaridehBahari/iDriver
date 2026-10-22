rm(list = ls())
library(data.table)
library(rtracklayer)
library(dplyr)
source('main/functions_iDriver.R')


map_muts <- function(gr, path_bed, is_gr_tetsGE = F){
  
  if (!is_gr_tetsGE) {
    testGE_gr <- import.bed(path_bed)
  } else {
    testGE_gr = path_bed
  }
  
  
  names(testGE_gr) = testGE_gr$name
  
  GenomicElement <- unlist(lapply(strsplit(names(testGE_gr), "[::]"), function(x){x[1]}))
  
  lenElement <- width(testGE_gr)
  
  mcols(testGE_gr) <- DataFrame(mcols(testGE_gr), GenomicElement, lenElement)
  
  callable_GE <- split(testGE_gr, names(testGE_gr))
  lenElement_new <- sum(width(callable_GE))
  
  ov <- findOverlaps(gr, callable_GE, ignore.strand = TRUE)
  
  idx_gr <- queryHits(ov)
  idx_callable_GE <- subjectHits(ov)
  
  gr_mappedGE <- gr[idx_gr]
  
  GE_col=unique(mcols(callable_GE, level = "within")[,"GenomicElement"])[idx_callable_GE]
  EL_col = sum(mcols(callable_GE, level = "within")[,"lenElement"])[idx_callable_GE]
  name_col = names(callable_GE)[idx_callable_GE]
  
  mcols(gr_mappedGE) <- DataFrame(mcols(gr_mappedGE),
                                  GenomicElement=GE_col,
                                  elemenntLength = EL_col,
                                  name = name_col)
  
  idx_completely_black <- which(lenElement_new == 0)
  
  
  
  if(length(idx_completely_black) != 0){
    lenElement_new <- lenElement_new[-idx_completely_black]
    callable_GE <- callable_GE[-idx_completely_black]
  }
  
  
  list(MAF_GE_mapped = gr_mappedGE, all_GEs_grl = callable_GE,
       lenAllTests = lenElement_new)
  
}


create_response <- function(all_elements, n_donors){
  
  obs_Muts <- all_elements$MAF_GE_mapped
  df <- as.data.frame(DataFrame(obs_Muts))
  df <- df[,c('name', 'D_id')]
  
  
  resTab <- df %>% 
    count(name, D_id) %>%
    filter(n != 0) %>%
    group_by(binID = name) %>%
    summarise(nMut = sum(n), nSample = n(), .groups = "drop")
  
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



save_responseTable <- function(gr, path_bed, path_out, save_name){
  
  all_elements <- map_muts(gr, path_bed)
  
  # print("************************ mutations mapped to the genomic intervals************************")
  n_donors <- length(unique(gr$D_id))
  resTab <- create_response(all_elements, n_donors)
  # print("************************ response table was generated ************************")
  dir.create(paste0(path_out, "/"), showWarnings = FALSE, recursive = TRUE)
  
  fwrite(resTab, paste0(path_out, "/", save_name,
                        ".tsv"), 
         sep = "\t", row.names = F) 
  # print("************************file saved************************")
  
}


################################# RUN #####################################
path_beds <- c('../../../../Projects/bahari_work/make_features/external/database/bins/proccessed/PCAWG_callable.bed6',
                 '../extdata/procInput0/callable_DP_trainGEs.bed6')
              
path_donorInfo <- '../extdata/procInput/iDriverInputs/donorInfo.tsv'
path_to_gr='../extdata/procInput/mut_PCAWG/all_samples.RData'
base_dir <- '../extdata/procInput/BMRs/observed/DP/'
save_names = c('test_y', 'train_y')
study = 'observed'
    
# exclude_hyper_mutated =T
exclude_hyper_mutated =F
only_hyperMut_pancancer <- F
load(path_to_gr)






    if (exclude_hyper_mutated ) {
      base_dir <- paste0(base_dir, 'rmHyperMutated/')



      cohorts <- c(
        'Pancan-no-skin-melanoma-lymph', 'Pan_Cancer',"Liver-HCC",
        "Bladder-TCC", "ColoRect-AdenoCA", "Lymph-BNHL",
        "Uterus-AdenoCA", "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
        "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC", "Biliary-AdenoCA",
        "Breast-AdenoCa", "Eso-AdenoCa", "CNS-GBM",
        "Panc-AdenoCA", "Lung-AdenoCA", "Prost-AdenoCA", "Ovary-AdenoCA",
        "Breast-LobularCa", "Thy-AdenoCA", "Myeloid-MPN",
        "Bone-Leiomyo",   "Lymph-NOS", "CNS-Medullo", "Myeloid-AML", "CNS-Oligo", "Cervix-SCC",
        "CNS-PiloAstro", "Kidney-ChRCC", "Bone-Epith", "Bone-Osteosarc",
        "Cervix-AdenoCA", "Breast-DCIS", "Bone-Cart", "Myeloid-MDS" )

    } else {
      if (only_hyperMut_pancancer) {
        base_dir <- paste0(base_dir, 'OnlyHyperMutated/')
      } else {
        base_dir <- paste0(base_dir, 'withHyperMutated/')

        cohorts=c("Pan_Cancer", "Pancan-no-skin-melanoma-lymph" ,
                  "ColoRect-AdenoCA","CNS-GBM",
                  "Head-SCC", "Skin-Melanoma", "Lung-AdenoCA", "Lung-SCC",
                  "Stomach-AdenoCA", "Uterus-AdenoCA" ,  "Biliary-AdenoCA",
                  "Eso-AdenoCa", "Lymph-BNHL"

        )
      }

    }





    for (cohort in cohorts) {

      print(cohort)
      path_out <- paste0(base_dir, cohort, '/')
      exclude_lymph_melanoma = ifelse(cohort == 'Pancan-no-skin-melanoma-lymph', T, F)

      project_code = cohort


      print(project_code)
      donor_info <- select_cohort(path_donorInfo,
                                  cohort, exclude_lymph_melanoma,
                                  exclude_hyper_mutated)
      gr_cohort <- gr[which(gr$D_id %in% donor_info$donor_id)]



      for(i in 1:length(path_beds)){
        save_responseTable(gr_cohort, path_beds[i], path_out, save_names[i])
      }
      print("************************files saved************************")
    }

########################################### python ###########################################################

