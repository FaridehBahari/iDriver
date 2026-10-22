
rm(list = ls())

library(data.table)

nonhypers <- c("Liver-HCC", "Bone-Leiomyo",   "Kidney-ChRCC",
               "Bladder-TCC",
               "Bone-Osteosarc",
               "Breast-AdenoCa",  "CNS-Medullo",
               "CNS-PiloAstro", 
               "Kidney-RCC", 
               
               "Lymph-CLL",  "Ovary-AdenoCA",
               "Panc-AdenoCA", "Panc-Endocrine", "Prost-AdenoCA"
)

cohorts_hyper <- c("Pancan-no-skin-melanoma-lymph", "Pan_Cancer",
                   "ColoRect-AdenoCA","CNS-GBM", 
                   "Head-SCC", "Skin-Melanoma", "Lung-AdenoCA", "Lung-SCC",
                   "Stomach-AdenoCA", "Uterus-AdenoCA" ,  "Biliary-AdenoCA",
                   "Eso-AdenoCa", "Lymph-BNHL")

cohorts = c(nonhypers, cohorts_hyper)


diffs <- list()

for (cohort in cohorts) {
  tissue = cohort
  
  panCanCoh <- cohort
  if (cohort == 'Pancan-no-skin-melanoma-lymph') {
    panCanCoh = 'myPancancerStar'
  } else if (cohort == 'Pan_Cancer'){
    panCanCoh = 'myPan_Cancer'
  }
  cat("\nProcessing:", cohort, "\n")
  
  
  tissue <- sub("AdenoCA$", "AdenoCa", cohort)
  if (tissue == 'Pan_Cancer') {
    tissue = 'Pancan'
  }
  ## -----------------------------
  ## Read DriverPower results
  ## -----------------------------
  dp <- fread(paste0('../../Dig/output/elemDriver/withHyper/', panCanCoh, '.results.txt'))
  if (panCanCoh ==  "myPancancerStar") {
    Ns <- fread(paste0('../../Dig/inputs/mutation_files/PCAWG/ICGC_and_TCGA/my_pancanStar.tsv'))
    
  } else{
    Ns <- fread(paste0('../../Dig/inputs/mutation_files/PCAWG/ICGC_and_TCGA/', panCanCoh, '_SNV_MNV_INDEL.annot.txt'))
    
  }
   N <- length(unique(Ns$V6))
  dp$nMut <- dp$OBS_SNV + dp$OBS_INDEL
  dp$nSample <- dp$OBS_SAMPLES
  dp$N <- N
  dp$element_ID <- dp$ELT
  ## -----------------------------
  ## Read input elements
  ## -----------------------------
  
  if (cohort %in% nonhypers) {
    id <- fread(file.path(
      "../extdata/procInput/BMRs/observed/rmHyperMutated/",
      cohort, "test_y.tsv"
    ))
  } else {
    id <- fread(file.path(
      "../extdata/procInput/BMRs/observed/withHyperMutated",
      cohort, "test_y.tsv"
    ))
  }
  
  # if (cohort %in% nonhypers) {
  #   dp <- fread(file.path(
  #     "../extdata/procInput0/BMRs_2024/observed/",
  #     cohort, "test_y.tsv"
  #   ))
  # } else {
  #   dp <- fread(file.path(
  #     "../extdata/procInput0/BMRs_2024/observed_withHyperMutated/",
  #     cohort, "test_y.tsv"
  #   ))
  # }
  
  
  ## Remove unwanted elements
  id <- id[!grepl("ncrna", binID)]
  id <- id[!grepl("gc19_pc.ss", binID)]
  
  # dp <- dp[!grepl("ncrna", binID)]
  # dp <- dp[!grepl("gc19_pc.ss", binID)]
  # colnames(dp) <- c("element_ID",   "nMut",    "nSample", "length",  "N" )
  
  ## -----------------------------
  ## Elements in id but not dp
  ## -----------------------------
  
  n <- id[!binID %in% dp$element_ID]
  
  if (nrow(n) != 0) {
    print(nrow(n))
    print("************************")
  }
  
  ## -----------------------------
  ## Elements in dp but not id
  ## -----------------------------
  
  n <- dp[!element_ID %in% id$binID]
  
  
  ## Keep only elements present in id
  dp <- dp[element_ID %in% id$binID]
  
  ## -----------------------------
  ## Compare nMut, nSample, length
  ## -----------------------------
  
  cmp <- merge(
    id[, .(
      element_ID = binID,
      id_nMut = nMut,
      id_nSample = nSample,
      id_length = N
    )],
    dp[, .(
      element_ID,
      dp_nMut = nMut,
      dp_nSample = nSample,
      dp_length = N
    )],
    by = "element_ID",
    all = TRUE
  )
  
  cmp[, `:=`(
    nMut_same = id_nMut == dp_nMut,
    nSample_same = id_nSample == dp_nSample,
    length_same = id_length == dp_length
  )]
  
  ## Elements where at least one value differs
  diff <- cmp[
    !nMut_same | !nSample_same  | !length_same
    ]
  
  ## -----------------------------
  ## Print cohort if differences
  ## -----------------------------
  
  if (nrow(diff) != 0) {
    cat("\nDIFFERENCE FOUND:", cohort, "\n")
    print(nrow(diff))
    diffs[[cohort]] <- diff
  }
}


