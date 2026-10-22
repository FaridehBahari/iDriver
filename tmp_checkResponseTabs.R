
rm(list = ls())

library(data.table)

nonhypers <- c("Liver-HCC","Myeloid-MPN", "Bone-Leiomyo",   "Kidney-ChRCC",
               "Bladder-TCC",
               "Bone-Osteosarc",
               "Breast-AdenoCa",  "CNS-Medullo",
               "CNS-PiloAstro", 
               "Kidney-RCC", 
               
               "Lymph-CLL",  "Ovary-AdenoCA",
               "Panc-AdenoCA", "Panc-Endocrine", "Prost-AdenoCA",
               "Thy-AdenoCA"
)

cohorts_hyper <- c("Pancan-no-skin-melanoma-lymph",
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
    panCanCoh = 'PanCan_No_Melanoma_Lymph'
  }
  cat("\nProcessing:", cohort, "\n")
  
  # if(cohort == 'Biliary-AdenoCA'){
  #   tissue = 'Biliary-AdenoCa'
  # } else if (cohort == 'ColoRect-AdenoCA') {
  #   tissue = 'ColoRect-AdenoCa'
  # } else if (cohort == 'Lung-AdenoCA') {
  #   tissue = 'Lung-AdenoCa'
  # } else if (cohort == "Stomach-AdenoCA") {
  #   tissue = "Stomach-AdenoCa"
  # }
  tissue <- sub("AdenoCA$", "AdenoCa", cohort)
  if (tissue == 'Pan_Cancer') {
    tissue = 'Pancan'
  }
  ## -----------------------------
  ## Read DriverPower results
  ## -----------------------------
  
  dp <- rbind(
    fread(file.path(
      "../extdata/rawInput/PCAWG_12methods_output/p-values/observed",
      tissue, "DriverPower",
      paste0(panCanCoh, ".3utr.DriverPower.observed.txt")
    )),
    fread(file.path(
      "../extdata/rawInput/PCAWG_12methods_output/p-values/observed",
      tissue, "DriverPower",
      paste0(panCanCoh, ".5utr.DriverPower.observed.txt")
    )),
    fread(file.path(
      "../extdata/rawInput/PCAWG_12methods_output/p-values/observed",
      tissue, "DriverPower",
      paste0(panCanCoh, ".CDS.DriverPower.observed.txt")
    )),
    fread(file.path(
      "../extdata/rawInput/PCAWG_12methods_output/p-values/observed",
      tissue, "DriverPower",
      paste0(panCanCoh, ".enhancers.DriverPower.observed.txt")
    )),
    fread(file.path(
      "../extdata/rawInput/PCAWG_12methods_output/p-values/observed",
      tissue, "DriverPower",
      paste0(panCanCoh, ".promCore.DriverPower.observed.txt")
    ))
  )
  
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
  
  print(unique(n$nMut))
  
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
      id_length = length
    )],
    dp[, .(
      element_ID,
      dp_nMut = nMut,
      dp_nSample = nSample,
      dp_length = length
    )],
    by = "element_ID",
    all = TRUE
  )
  
  cmp[, `:=`(
    nMut_same = id_nMut == dp_nMut,
    nSample_same = id_nSample == dp_nSample
    # length_same = id_length == dp_length
  )]
  
  ## Elements where at least one value differs
  diff <- cmp[
    !nMut_same | !nSample_same # | !length_same
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

'gc19_pc.3utr::gencode::GRM4::ENSG00000124493.9'