rm(list = ls())

library(GenomicRanges)
library(rtracklayer)
library(dplyr)
library(data.table)
library(doParallel)
library(foreach)

source('main/functions_iDriver.R')
source('prepareInput/4_Functions_responseTable.R')


# path_bed <- '../../../../Projects/bahari_work/make_features/external/database/bins/proccessed/PCAWG_callable.bed6'
# win_size <- 100000
# num_cores <- 40
# 
# define_background <- function(path_bed, win_size, num_cores){
#   
#   gr_elements <- import.bed(path_bed)
#   gr_elements$blockID <- paste0('b', 1:length(gr_elements))
#   
#   # gr_extended <- resize(gr_elements, width = width(gr_elements) + 100000, fix = "center")
#   # gr_split <- split(gr_extended, gr_extended$name)
#   gr_split <- resize(gr_elements, width = width(gr_elements) + win_size, fix = "center")
#   gr_split <- split(gr_split, gr_split$name)
#   
#   registerDoParallel(num_cores)
#   
#   t0 = Sys.time()
#   # Step 3: Parallelize the loop using foreach
#   filtered_bg <- foreach(i = 1:length(gr_split), .combine = c, .packages = "GenomicRanges") %dopar% {
#     if (i %% 5000 == 0) {
#       print(i) # i = 121790
#     }
#     
#     gr_elem <- gr_elements[which(gr_elements$name == names(gr_split)[i])]
#     bg_containElem <- gr_split[[i]]
#     gr_elem$length <- width(gr_elem)
#     
#     filtered_bg_grs <- GRanges()
#     # Remove overlaps within each element
#     for (j in 1:length(bg_containElem)) {
#       filtered_bg_gr <- GenomicRanges::setdiff(bg_containElem[j], gr_elem)
#       filtered_bg_gr$name <- unique(gr_elem$name)
#       ID = bg_containElem[j]$blockID
#       filtered_bg_gr$blockID <- ID
#       filtered_bg_gr$length = gr_elem[which(gr_elem$blockID == ID)]$length
#       filtered_bg_grs <- c(filtered_bg_grs, filtered_bg_gr)
#     }
#     
#     return(filtered_bg_grs)
#   }
#   
#   
#   
#   print(Sys.time() - t0)
#   
#   filtered_bg
# }
# 
# 
# 
# filtered_bg <- define_background(path_bed, win_size, num_cores)

############################### create block-wise background response tables ####################
# load('../fitered_bg.RData')
# bg_bed <- filtered_bg[!grepl('lncrna', filtered_bg$name)]
# bg_bed <- bg_bed[!grepl('gc19_pc.ss', bg_bed$name)]
# save(bg_bed, file = '../bg_bed.RData')

path_save <- '../extdata/procInput/BMRs/observed/baseLine_localRates/rmHyperMutated/'
dir.create(path_save, recursive=T)

load('../bg_bed.RData')
bg_bed$PCAWG_ID <- bg_bed$name
bg_bed$name <- bg_bed$blockID
df <- data.frame('PCAWG_ID' = bg_bed$PCAWG_ID,
                 'elem_ID' = bg_bed$name, 
                 'PCAWG_block_length' = bg_bed$length)

df <- df[!duplicated(df),]


load('../extdata/procInput/mut_PCAWG/all_samples.RData')

path_donorInfo <- '../extdata/procInput/iDriverInputs/donorInfo.tsv'


cohorts <- c(#  
             "Liver-HCC" ,
  "Bladder-TCC", "ColoRect-AdenoCA", "Lymph-BNHL",
  "Uterus-AdenoCA", "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
  "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
  "Breast-AdenoCa","Biliary-AdenoCA", "Eso-AdenoCa", "CNS-GBM",
  "Panc-AdenoCA", "Lung-AdenoCA", "Prost-AdenoCA", "Ovary-AdenoCA",
  "Thy-AdenoCA", "Myeloid-MPN", "Bone-Leiomyo",
  "CNS-Medullo", "CNS-Oligo", "Cervix-SCC",
  "CNS-PiloAstro", "Kidney-ChRCC", "Bone-Osteosarc",
  "Breast-LobularCa",  "Lymph-NOS", "Myeloid-AML", "Bone-Epith",
  "Cervix-AdenoCA","Breast-DCIS","Bone-Cart","Myeloid-MDS"# 'Pancan-no-skin-melanoma-lymph', "Pan_Cancer"  #
) 


for (cohort in cohorts) {
  
  print(cohort)
  
  exclude_lymph_melanoma <- ifelse(cohort == 'Pancan-no-skin-melanoma-lymph', T, F)
  
  donorInfo <- select_cohort(path_donorInfo,
                             cohort,
                             exclude_lymph_melanoma,
                             exclude_hyper_mutated = T)
  
  gr_cohort <- gr[which(gr$D_id %in% donorInfo$donor_id)]
  
  
  t0 <- Sys.time()
  all_elements <- map_muts(gr_cohort, bg_bed, is_gr_tetsGE = T)
  print(paste0('map_muts total time: ', Sys.time() - t0))
  
  n_donors <- length(unique(gr_cohort$D_id))
  
  t0 <- Sys.time()
  resTab <- create_response(all_elements, n_donors)
  print(paste0('create_response total time: ', Sys.time() - t0))
  
  
  
  pred_tab <- left_join(resTab, df, by = c('binID'= 'elem_ID'))
  pred_tab <- pred_tab[!duplicated(pred_tab),]
  # pred_tab <- as.data.frame(left_join(resTab, df, by = c('binID'= 'elem_ID')))
  # pred_tab <- pred_tab[,which(!colnames(pred_tab) %in% c('binID', 'PCAWG_block_length'))]
  # pred_tab <- pred_tab[!duplicated(pred_tab),]
  
  
  fwrite(pred_tab, file = paste0(path_save,'/baselineBMR_', cohort, '.tsv'), sep = '\t')
}



############################### create block-wise background response tables pancancers ####################
rm(list = ls())

library(dplyr)
library(data.table)

indir <- "../extdata/procInput/BMRs/observed/baseLine_localRates/rmHyperMutated/"

# Read and combine all cancer-specific baseline BMR files
files <- list.files(
  indir,
  pattern = "^baselineBMR_.*\\.tsv$",
  full.names = TRUE
)


dt <- lapply(files, fread)


by_cols <- c("binID", "length", "PCAWG_ID", "PCAWG_block_length")


x <- full_join(dt[[1]], dt[[2]], by = by_cols)

for (i in 3:length(dt)) {
  x <- full_join(x, dt[[i]], by = by_cols)
  
}

nMut_cols <- grep('nMut', colnames(x)) 
nSample_cols <- grep('nSample', colnames(x)) 
N_cols <- grep('N', colnames(x)) 

x <- as.data.frame(x)
y <- x[, which(colnames(x) %in% by_cols)]


y$nMut <- rowSums(x[, nMut_cols], na.rm = TRUE)
y$nSample <- rowSums(x[, nSample_cols], na.rm = TRUE)
y$N <- rowSums(x[, N_cols], na.rm = TRUE)

 "Pan_Cancer" 
# Save the pan-cancer table
fwrite(
  y,
  file = file.path(indir, "baselineBMR_Pan_Cancer.tsv"),
  sep = "\t"
)

rm(list = ls())
indir <- "../extdata/procInput/BMRs/observed/baseLine_localRates/"

'Pancan-no-skin-melanoma-lymph'

files <- list.files(
  indir,
  pattern = "^baselineBMR_.*\\.tsv$",
  full.names = TRUE
)

exceptions <- c("Skin-Melanoma", "SKCM-US",
                "Lymph-NOS", 
                "Lymph-CLL", "CLLE-ES",
                "Lymph-BNHL", "MALY-DE", "DLBC-US", 
                "Pan_Cancer")

files <- files[!grepl(
  paste(exceptions, collapse = "|"),
  basename(files)
)]
dt <- lapply(files, fread)
by_cols <- c("binID", "length", "PCAWG_ID", "PCAWG_block_length")
x <- full_join(dt[[1]], dt[[2]], by = by_cols)

for (i in 3:length(dt)) {
  x <- full_join(x, dt[[i]], by = by_cols)
  
}

nMut_cols <- grep('nMut', colnames(x)) 
nSample_cols <- grep('nSample', colnames(x)) 
N_cols <- grep('N', colnames(x)) 

x <- as.data.frame(x)
y <- x[, which(colnames(x) %in% by_cols)]


y$nMut <- rowSums(x[, nMut_cols], na.rm = TRUE)
y$nSample <- rowSums(x[, nSample_cols], na.rm = TRUE)
y$N <- rowSums(x[, N_cols], na.rm = TRUE)

fwrite(
  y,
  file = file.path(indir, "baselineBMR_Pancan-no-skin-melanoma-lymph.tsv"),
  sep = "\t"
)

################################# estimate mutation rate based on weighted means of block lengths ###############################################
rm(list = ls())

library(data.table)
library(dplyr)

path_save <- "../extdata/procInput/BMRs/observed/baseLine_localRates/"

dir.create(paste0(path_save, '/predRates/'), showWarnings = F, recursive = T)

cohorts <- c("Pan_Cancer" , 'Pancan-no-skin-melanoma-lymph', "Liver-HCC" ,
             "Bladder-TCC", "ColoRect-AdenoCA", "Lymph-BNHL",
             "Uterus-AdenoCA", "Kidney-RCC", "Lymph-CLL", "Lung-SCC",
             "Stomach-AdenoCA", "Skin-Melanoma", "Panc-Endocrine", "Head-SCC",
             "Breast-AdenoCa","Biliary-AdenoCA", "Eso-AdenoCa", "CNS-GBM",
             "Panc-AdenoCA", "Lung-AdenoCA", "Prost-AdenoCA", "Ovary-AdenoCA",
             "Thy-AdenoCA", "Myeloid-MPN", "Bone-Leiomyo",
             "CNS-Medullo", "CNS-Oligo", "Cervix-SCC",
             "CNS-PiloAstro", "Kidney-ChRCC", "Bone-Osteosarc",
             "Breast-LobularCa",  "Lymph-NOS", "Myeloid-AML", "Bone-Epith",
             "Cervix-AdenoCA","Breast-DCIS","Bone-Cart" #,"Myeloid-MDS"
) 


for (cohort in cohorts) {
  
  print(cohort)
  
  test_y <- fread(paste0('../../iDriver/extdata/procInput/BMRs/observed/rmHyperMutated/', cohort,
                         '/test_y.tsv'))
  test_y$obsRate <- test_y$nMut/(test_y$length * test_y$N)
  N <- unique(test_y$N)
  
  final_resTab <- fread(paste0(path_save, '/baselineBMR_',cohort, '.tsv'))
  
  obs_tab <- test_y[,c('binID', 'obsRate')]
  df <- left_join(final_resTab, obs_tab, by = c('PCAWG_ID' = 'binID'))
  df$nMut = as.numeric( df$nMut)
  df$length = as.numeric( df$length)
  df$PCAWG_block_length = as.numeric(df$PCAWG_block_length)
  
   
  # PCAWG_ID = 'gc19_pc.3utr::gencode::C1orf170::ENSG00000187642.5'
  # df <- as.data.frame(df)
  # elem=df[df$PCAWG_ID == PCAWG_ID,]
  # 
  # PCAWG_length = sum(elem$PCAWG_block_length)
  # length_correction <- PCAWG_length/(bg_length/nrow(elem))
  # 
  # bg_length <- sum(elem$length)
  # weights_elem <- elem$PCAWG_block_length / PCAWG_length
  # pred_elem_mut = (sum((elem$nMut)*weights_elem) / nrow(elem))*length_correction
  # pred_elem_smpl = (sum(elem$nSample*weights_elem) / nrow(elem))*length_correction
  # 
  # 
  # 
  df$rates = df$nMut / df$length
  # Efficient calculation using dplyr
  results <- df %>%
    group_by(PCAWG_ID) %>%
    summarise(
      total_length = sum(PCAWG_block_length), # element length
       #  rates number of mutations in a 100k window (for each block of an element)
      pred_elem = sum((PCAWG_block_length / sum(PCAWG_block_length)) * rates) / n(), # weighted mean of rates of the blocks of an element
      obs_elem = unique(obsRate), #observed rate 
      .groups = "drop") %>%
    ungroup()
  
  results2 <- results[,c('PCAWG_ID', 'total_length', 'pred_elem', 'obs_elem')]
  results2 <- results2[!duplicated(results2),]
  results2 <- results2[which(results2$obs_elem != 0),]
  colnames(results2) <- c('binID', 'length',     'predRate',   'obs_rates')
  results2$predRate <- results2$predRate/N
  
  fwrite(results2, paste0(path_save, '/predRates/',cohort, '.tsv'), sep = '\t')
  
  
  
}




