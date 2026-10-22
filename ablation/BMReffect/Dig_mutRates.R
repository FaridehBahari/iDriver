rm(list = ls())

library(data.table)



################################################

nonhypers <- c("Kidney-ChRCC", "CNS-PiloAstro",  #"Myeloid-MPN",  "Thy-AdenoCA"
               "Liver-HCC","Bone-Leiomyo",
               "Bladder-TCC",
               "Bone-Osteosarc",
               "Breast-AdenoCa",  "CNS-Medullo",

               "Kidney-RCC",

               "Lymph-CLL",  "Ovary-AdenoCA",
               "Panc-AdenoCA", "Panc-Endocrine", "Prost-AdenoCA"

)

cohorts_hyper <- c("Pan_Cancer", "Pancan-no-skin-melanoma-lymph",
                   "ColoRect-AdenoCA","CNS-GBM",
                   "Head-SCC", "Skin-Melanoma", "Lung-AdenoCA", "Lung-SCC",
                   "Stomach-AdenoCA", "Uterus-AdenoCA" ,  "Biliary-AdenoCA",
                   "Eso-AdenoCa", "Lymph-BNHL"
  )

cohorts = c( cohorts_hyper, nonhypers)

for (cohort in cohorts) {


  panCanCoh <- cohort
  if (cohort == 'Pan_Cancer') {
    panCanCoh = 'myPan_Cancer'
  }
  cat("\nProcessing:", cohort, "\n")

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

  id <- id[!grepl('gc19_pc.ss', id$binID),]
  id <- id[!grepl('ncrna', id$binID),]


  df <- fread(paste0('../../Dig/output/elemDriver/withHyper/', panCanCoh,'.results.txt'))


  # print(nrow(df))
  # print('------')
  Dig_N <- fread(paste0('../../Dig/inputs/mutation_files/PCAWG/ICGC_and_TCGA/', panCanCoh, '_SNV_MNV_INDEL.annot.txt'))

  Dig_N <- length(unique(Dig_N$V6))
  print(Dig_N)


  df$N <- Dig_N
  df$obs_rates <- df$R_OBS / (df$R_SIZE*df$N )
  df$predRate <- df$MU / (df$R_SIZE*df$N )

  df$nMut <- df$OBS_SNV + df$OBS_INDEL

  df$sumEXP <- df$EXP_SNV + df$EXP_INDEL

  df$obs_rates0 <- (df$nMut) / (df$ELT_SIZE*df$N )
  df$predRate0 <- df$sumEXP / (df$ELT_SIZE*df$N )

  df$nSample <- df$OBS_SAMPLES

   df <- df[,c("ELT",  "ELT_SIZE", "nMut", "nSample", "R_SIZE",  "R_OBS",  "MU", "sumEXP",
              "obs_rates0", "predRate0", "obs_rates", "predRate")]


  df <- df[df$ELT %in% id$binID,]

  dir.create('../extdata/procInput/BMRs/observed/Dig/withHyperMutated/', recursive = T, showWarnings = F)
  fwrite(df, file = paste0('../extdata/procInput/BMRs/observed/Dig/withHyperMutated/', cohort, 'predRates.tsv'), sep = '\t')
  # print(paste0(cohort, " %overlap with iDriver elements: ", nrow(df)/nrow(id)))

}




################################################
rm(list =ls())
# Dig was ran on cohorts with hypermutated samples by removing hypers and just retaining nonhypers

cohorts_hyper <- c('Pancan-no-skin-melanoma-lymph', "Pan_Cancer"
                   # "ColoRect-AdenoCA","CNS-GBM",
                   # "Head-SCC", "Skin-Melanoma", "Lung-AdenoCA", "Lung-SCC",
                   # "Stomach-AdenoCA", "Uterus-AdenoCA" ,  "Biliary-AdenoCA",
                   # "Eso-AdenoCa", "Lymph-BNHL"
)
nonhypers <- c()
cohorts = c( cohorts_hyper, nonhypers)

for (cohort in cohorts) {
  
    panCanCoh <- cohort
    if (cohort == 'Eso-AdenoCA') {
      panCanCoh = 'Eso-AdenoCa'
    } else if (cohort == 'Pancan-no-skin-melanoma-lymph') {
      panCanCoh = 'myPancancerStar'
    }
  
  cat("\nProcessing:", cohort, "\n")
  
  id <- fread(file.path(
    "../extdata/procInput/BMRs/observed/rmHyperMutated/",
    cohort, "test_y.tsv"
  ))
  
  
  id <- id[!grepl('gc19_pc.ss', id$binID),]
  id <- id[!grepl('ncrna', id$binID),]
  
  if (cohort == 'Pancan-no-skin-melanoma-lymph') {
    df <- fread(paste0('../../Dig/output/elemDriver/rmHypers/', panCanCoh,'.results.txt'))
    Dig_N <- fread('../../Dig/inputs/mutation_files/PCAWG/ICGC_and_TCGA/rmHypers/my_pancanStar.tsv')
  } else {
    df <- fread(paste0('../../Dig/output/elemDriver/rmHypers/', cohort,'.results.txt'))
    if (cohort == 'Pan_Cancer') {
      Dig_N <- fread(paste0('../../Dig/inputs/mutation_files/PCAWG/ICGC_and_TCGA/rmHypers/Pan_Cancer.tsv'))
    } else {
      Dig_N <- fread(paste0('../../Dig/inputs/mutation_files/PCAWG/ICGC_and_TCGA/', panCanCoh,
                            '_SNV_MNV_INDEL.annot.txt'))
    }
    
    
  }
  
  
  # print(nrow(df))
  # print('------')
  
  Dig_N <- length(unique(Dig_N$V6))
  print(Dig_N)
  
  
  df$N <- Dig_N
  df$obs_rates <- df$R_OBS / (df$R_SIZE*df$N )
  df$predRate <- df$MU / (df$R_SIZE*df$N )
  
  df$nMut <- df$OBS_SNV + df$OBS_INDEL
  
  df$sumEXP <- df$EXP_SNV + df$EXP_INDEL
  
  df$obs_rates0 <- (df$nMut) / (df$ELT_SIZE*df$N )
  df$predRate0 <- df$sumEXP / (df$ELT_SIZE*df$N )
  
  df$nSample <- df$OBS_SAMPLES
  
  df <- df[,c("ELT",  "ELT_SIZE", "nMut", "nSample", "R_SIZE",  "R_OBS",  "MU", "sumEXP",
              "obs_rates0", "predRate0", "obs_rates", "predRate")]
  
  
  df <- df[df$ELT %in% id$binID,]
  
  dir.create('../extdata/procInput/BMRs/observed/Dig/rmHyperMutated/', recursive = T, showWarnings = F)
  fwrite(df, file = paste0('../extdata/procInput/BMRs/observed/Dig/rmHyperMutated/', cohort, 'predRates.tsv'), sep = '\t')
  # print(paste0(cohort, " %overlap with iDriver elements: ", nrow(df)/nrow(id)))
  
}
