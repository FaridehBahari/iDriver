# maf <- fread("../extdata/procInput/mut_PCAWG/remove_duplicatedSamples/MAF_samples_removed.tsv")
# exceptions <- c("Skin-Melanoma", "Lymph-NOS", "Lymph-CLL", "Lymph-BNHL")
# pancan_noMelnoLymph <- maf[!(maf$Project_Code %in% exceptions),]
rm(list = ls())
library(data.table)
library(readxl)

maf_liver <- fread('../extdata/procInput/mut_PCAWG/tmp_pancan_noLM.tsv')
length(unique(maf_liver$Donor_ID))
length(unique(maf_liver$Tumor_Sample_Barcode))

allSamp_projs0 <- data.frame(cbind(maf_liver$Tumor_Sample_Barcode, maf_liver$Donor_ID, maf_liver$Project_Code))
allSamp_projs0 <- allSamp_projs0[!duplicated(allSamp_projs0),]
duplicatedDonors <- unique(allSamp_projs0$X2[duplicated(allSamp_projs0$X2)])

allSamp_projs0$dups <- ifelse(allSamp_projs0$X2 %in% duplicatedDonors, T, F)

removed_samples <- c("f8e61a02-654c-c226-e040-11ac0d481b60", "f8e61a02-8fbe-a9af-e040-11ac0d481b64",
                     "f8e61a02-90ec-96fd-e040-11ac0d481b5e",
                     # "f93f4445-6f19-403c-e040-11ac0d48556f", 
                     "f8e61a02-8fbe-a9af-e040-11ac0d481b64", "f90a074e-d82f-f6bb-e040-11ac0d487334"
                     )
# "f901353f-fdeb-5a55-e040-11ac0c482e89" "f90a074e-d82f-f6bb-e040-11ac0d487334"  # one should be removed MDS duplicate

info2 <- fread('../extdata/rawInput/pcawg_sample_sheet.tsv')
included <- as.data.frame(info2[which(info2$donor_wgs_exclusion_white_gray == "Whitelist"),
                                "aliquot_id"], row.names = c())

maf_liver <- maf_liver[maf_liver$Tumor_Sample_Barcode %in% included$aliquot_id,]


length(unique(maf_liver$Donor_ID))
length(unique(maf_liver$Tumor_Sample_Barcode))

maf_liver <- maf_liver[!maf_liver$Tumor_Sample_Barcode %in% removed_samples,]
maf_liver <- maf_liver[!duplicated(maf_liver),]

length(unique(maf_liver$Donor_ID))
length(unique(maf_liver$Tumor_Sample_Barcode))


# info <- read_xlsx('../extdata/rawInput/PCAWG_supp/pcawg_donor_clinical_August2016_v9.xlsx')
# included <- as.data.frame(info[which(info$donor_wgs_included_excluded == "Included"),
#                                "icgc_donor_id"], row.names = c())
# maf_liver <- maf_liver[maf_liver$Donor_ID %in% included$icgc_donor_id,]


allSamp_projs1 <- data.frame(cbind(maf_liver$Tumor_Sample_Barcode, maf_liver$Donor_ID))
allSamp_projs1 <- allSamp_projs1[!duplicated(allSamp_projs1),]
duplicatedDonors <- unique(allSamp_projs1$X2[duplicated(allSamp_projs1$X2)])
allSamp_projs1$dups <- ifelse(allSamp_projs1$X2 %in% duplicatedDonors, T, F)

do_multipleSampless <- (allSamp_projs1[allSamp_projs1$dups, c('X2')])
do_multipleSampless <- allSamp_projs1[allSamp_projs1$X2 %in% do_multipleSampless, ]

do_UniqSampless <- allSamp_projs1[!allSamp_projs1$X2 %in% do_multipleSampless$X2,]

do_UniqSampless$dups <- FALSE
do_multipleSampless$dups <- TRUE

library(data.table)
library(rtracklayer)
library(dplyr)

all_maf <- rbind(do_UniqSampless, do_multipleSampless)
colnames(all_maf) <- c("Tumor_Sample_Barcode", "Donor_ID", "dups")
CADDannotated_maf <- left_join(maf_liver, all_maf, by = c("Tumor_Sample_Barcode", "Donor_ID"))

gr <- GenomicRanges::GRanges(seqnames = paste0("chr", CADDannotated_maf$Chromosome),
                             #strand = CADDannotated_maf$Strand,
                             IRanges::IRanges(start = CADDannotated_maf$Start_position,
                                              end = CADDannotated_maf$End_position),
                             ref_al = CADDannotated_maf$Reference_Allele,
                             alt_al = CADDannotated_maf$Tumor_Seq_Allele2,
                             D_id = CADDannotated_maf$Donor_ID,
                             sampleID = CADDannotated_maf$Tumor_Sample_Barcode,
                             cohort = CADDannotated_maf$Project_Code,
                             Ref_cntxt = CADDannotated_maf$ref_context,
                             var_type = CADDannotated_maf$Variant_Type,
                             variant_class = CADDannotated_maf$Variant_Classification
                             
                             # uniqueDonor_sample = CADDannotated_maf$dups
                             )

strand(gr) <- Rle(rep("+", nrow(CADDannotated_maf)))



#

source('main/functions_iDriver.R')


map_muts <- function(gr, path_bed){
  
  testGE_gr <- import.bed(path_bed)
  
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

#
#
#
path_bed <- '../../../../Projects/bahari_work/make_features/external/database/bins/proccessed/PCAWG_callable.bed6'
all_elements <- map_muts(gr, path_bed)
MAF_GE_mapped <- all_elements$MAF_GE_mapped
# save(MAF_GE_mapped, file = '../extdata/procInput/mut_PCAWG/tmp_MAF_GE_mapped_allSamples.RData')

##################################################################################
# load('../extdata/procInput/mut_PCAWG/tmp_MAF_GE_mapped_allSamples.RData')

n_donors <- length(unique(gr$D_id))
resTab <- create_response(all_elements, n_donors)


############################# select diffs id & dp ##########################


id <- as.data.frame(resTab)
dp <- rbind(
  fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Pancan-no-skin-melanoma-lymph//DriverPower/PanCan_No_Melanoma_Lymph.3utr.DriverPower.observed.txt"
  ),
  fread(file.path("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Pancan-no-skin-melanoma-lymph/DriverPower/PanCan_No_Melanoma_Lymph.5utr.DriverPower.observed.txt"
  )),
  fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Pancan-no-skin-melanoma-lymph/DriverPower/PanCan_No_Melanoma_Lymph.CDS.DriverPower.observed.txt")
  ,
  fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Pancan-no-skin-melanoma-lymph/DriverPower/PanCan_No_Melanoma_Lymph.enhancers.DriverPower.observed.txt")
  ,
  fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Pancan-no-skin-melanoma-lymph/DriverPower/PanCan_No_Melanoma_Lymph.promCore.DriverPower.observed.txt")

# )
# dp <- rbind(
#   fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Panc-AdenoCa/DriverPower/Panc-AdenoCA.3utr.DriverPower.observed.txt"
#   ),
#   fread(file.path("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Panc-AdenoCa/DriverPower/Panc-AdenoCA.5utr.DriverPower.observed.txt"
#   )),
#   fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Panc-AdenoCa/DriverPower/Panc-AdenoCA.CDS.DriverPower.observed.txt")
#   ,
#   fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Panc-AdenoCa/DriverPower/Panc-AdenoCA.enhancers.DriverPower.observed.txt")
#   ,
#   fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Panc-AdenoCa/DriverPower/Panc-AdenoCA.promCore.DriverPower.observed.txt")
# 
# )
# 
# dp <- rbind(
#   fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Myeloid-MPN/DriverPower/Myeloid-MPN.3utr.DriverPower.observed.txt"
#   ),
#   fread(file.path("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Myeloid-MPN/DriverPower/Myeloid-MPN.5utr.DriverPower.observed.txt"
#   )),
#   fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Myeloid-MPN/DriverPower/Myeloid-MPN.CDS.DriverPower.observed.txt")
#   ,
#   fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Myeloid-MPN/DriverPower/Myeloid-MPN.enhancers.DriverPower.observed.txt")
#   ,
#   fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Myeloid-MPN/DriverPower/Myeloid-MPN.promCore.DriverPower.observed.txt")

)

## Remove unwanted elements
id <- id[!grepl("ncrna", id$binID),]
id <- id[!grepl("gc19_pc.ss", id$binID),]

id <- left_join(id, dp, by = c('binID' = 'element_ID'  ))
id <- id[!is.na(id$nSample.y),]

id$dif <- as.numeric(id$nSample.x) - as.numeric(id$nSample.y)

sum(id$dif !=0)
id <- id[id$dif !=0,]
rownames(id) <- 1:nrow(id)

get_sample <- function(MAF_GE_mapped, binID){
  index <- which(MAF_GE_mapped$name %in% c(binID))
  
  elem_maf <- MAF_GE_mapped[index]
  # elem_maf[which(!elem_maf$uniqueDonor_sample)]
  print(unique(elem_maf$sampleID))
  # print(binID)
}

# get_sample(MAF_GE_mapped, 'gc19_pc.5utr::gencode::TBC1D20::ENSG00000125875.9')

get_sample(MAF_GE_mapped, id$binID[1])

df <- as.data.frame(table(maf_liver$Donor_ID, maf_liver$Project_Code))
df <- df[df$Freq !=0,]
df <- as.data.frame(table(df$Var2))


'gc19_pc.cds::gencode::OR1D5::ENSG00000262628.1'

"gc19_pc.5utr::gencode::TBC1D20::ENSG00000125875.9"
# "DO52727" "DO52735" "DO52747" in both myelMPN and AML

# rm(list = ls())
# library(data.table)
# library(readxl)
# 
# maf_liver <- fread('../extdata/procInput/mut_PCAWG/tmp_mafPROSTaden_allSamps.tsv')
# length(unique(maf_liver$Donor_ID))
# length(unique(maf_liver$Tumor_Sample_Barcode))
# 
# allSamp_projs0 <- data.frame(cbind(maf_liver$Tumor_Sample_Barcode, maf_liver$Donor_ID))
# allSamp_projs0 <- allSamp_projs0[!duplicated(allSamp_projs0),]
# duplicatedDonors <- unique(allSamp_projs0$X2[duplicated(allSamp_projs0$X2)])
# 
# allSamp_projs0$dups <- ifelse(allSamp_projs0$X2 %in% duplicatedDonors, T, F)
# 
# removed_samples <- unique(c("0bfd1068-3fca-a95b-e050-11ac0c4860c3", "5fe8a231-bbac-4f91-a47a-0de4fee752d3",
#   "0bfd1043-8189-e3e4-e050-11ac0c4860c5", "1d6c9181-a02f-4c35-85b3-65dc5808f50e",
#                      "95730d8b-327a-402c-a8bf-d71370f777a7", 
#                      "2b54d668-7585-4bf4-9fdd-aea3c5cd4fa4", "0bfd1043-7ec1-aaec-e050-11ac0c482f39",
#                      "2b54d668-7585-4bf4-9fdd-aea3c5cd4fa4", "2b54d668-7585-4bf4-9fdd-aea3c5cd4fa4",
#                      
#                      
#                      "2b54d668-7585-4bf4-9fdd-aea3c5cd4fa4",
#                      "5df0aab4-0be2-4da3-a4b7-a10fca09a4ba",
#                      "74db46a8-a499-4951-af33-1d53b7d1c9b4",
#                      "7e10e325-16d2-42df-b957-df4d3a9477d0",
#                      
#                      "052015bc-8c24-467c-8ff7-9e0dc4ad39fd", "0bd3a230-531e-44aa-9999-b8ed8da0176b",
#                      "0bfd1068-3fd5-a95b-e050-11ac0c4860c3",
#                      "0bfe2ac9-0af8-c248-e050-11ac0d487e1c", "0c0038ff-6cc4-b0b0-e050-11ac0d483d73",
#                      "0c0038ff-6cc6-b0b0-e050-11ac0d483d73",
#                      "339a2c94-55c9-4113-a706-b0971fa2784f",
#                      
#                      "4c755527-3b5d-4d36-822d-990aa1003d6a", "0bfd1043-7346-fdd0-e050-11ac0c484cab",
#                      "2f79fe89-9d24-4fbe-ba5b-0ff6c3c4d663", "fd163b56-1c2d-4e63-9fad-ba3221ae274b",
#                      "fd163b56-1c2d-4e63-9fad-ba3221ae274b",
#                      "3b8ab548-6669-4142-9612-cafded0d7a23" ,
#                      "3b8ab548-6669-4142-9612-cafded0d7a23", "fc9f9295-01f6-a9e0-e040-11ac0d4845af",
#                      "3b8ab548-6669-4142-9612-cafded0d7a23", "fc9f9295-01f6-a9e0-e040-11ac0d4845af",
#                      
#   
#   "0bfe2ac9-0afa-c248-e050-11ac0d487e1c", "0bfe2ac9-0aee-c248-e050-11ac0d487e1c",
#                      "fc9d93b6-77c6-8719-e040-11ac0d487dec", "0bfe2ac9-0afa-c248-e050-11ac0d487e1c",
#                      "0bfe2ac9-0aec-c248-e050-11ac0d487e1c", "0bfd1068-3fe6-a95b-e050-11ac0c4860c3",
#                      
#                      "0bfd1043-5142-3662-e050-11ac0c486501", "0bfd1068-3fdf-a95b-e050-11ac0c4860c3",
#                      "0c0038ff-6cc8-b0b0-e050-11ac0d483d73", "0bfebf9f-c781-e57d-e050-11ac0d487827",
#                      "0bfebf9f-c77d-e57d-e050-11ac0d487827", "fc9d5fc7-a87b-c5e7-e040-11ac0c4850aa",
#                      "b0f3087b-6f9c-47e9-b487-6726af1cb90b", "280a069c-38e5-4300-a15c-2612301b70e3",
#                      "7a86abdc-d554-43bd-ad47-22d439c930ef", "0bfe2ac9-0afd-c248-e050-11ac0d487e1c",
#                      
#                      "0bfe2ac9-0af5-c248-e050-11ac0d487e1c", "0bfd1043-7343-fdd0-e050-11ac0c484cab",
#                      "0bfe2ac9-0aff-c248-e050-11ac0d487e1c", "0bfd1043-817c-e3e4-e050-11ac0c4860c5",
#                      "fc9d5fc6-546f-27b7-e040-11ac0d485f87",
#                      
#                      "0bfd1043-817e-e3e4-e050-11ac0c4860c5", "84aa2983-714b-49fd-af93-f18ed4d6a6cb",
#                      "63faa59e-6649-412e-823a-338d0662ec72", "fc9d704c-0cfb-0407-e040-11ac0d486980",
#                      "0bfd1043-8173-e3e4-e050-11ac0c4860c5", "0bfd1043-8180-e3e4-e050-11ac0c4860c5",
#                      "f1d73c2f-3d7d-4a0f-871b-b3bf76d865dc", "0bfd1068-3fc3-a95b-e050-11ac0c4860c3",
#                      "0bfe2ac9-0af3-c248-e050-11ac0d487e1c", "1558629a-770a-4dcb-b0bf-9a8bcea9d80b",
#                      "db37739d-af50-4e63-9147-1a7d25630ac6", "fca04929-586a-0594-e040-11ac0c48099d",
#                      "0bfd1043-816e-e3e4-e050-11ac0c4860c5", "fca2d62f-2756-939e-e040-11ac0d4865e6",
#                      "0bfd1068-3fcd-a95b-e050-11ac0c4860c3", "0bfd1068-3fe1-a95b-e050-11ac0c4860c3",
#                      "0bfebf9f-c783-e57d-e050-11ac0d487827", "0bfebf9f-c77f-e57d-e050-11ac0d487827",
#                      
#                      "0bfd1043-7ed6-9ccc-e050-11ac0c481957", "9cdb1905-8b0f-48bd-ba37-11e7ac9646a0",
#                      "fc9f998b-cc72-f6ea-e040-11ac0c482425", "9b75735b-705d-4ebf-92a0-fbb190b95285",
#                      "0bfd1043-7344-fdd0-e050-11ac0c484cab", "6ff16b0d-5966-4119-9bec-c9ebff804a37",
#                      
#                      "0bfd1068-3fe4-a95b-e050-11ac0c4860c3", "0bfd1043-8172-e3e4-e050-11ac0c4860c5",
#                      "0bfd1043-7344-fdd0-e050-11ac0c484cab",
#                      "0bfd1043-8175-e3e4-e050-11ac0c4860c5", "1bd47e40-d708-4ca2-b4b3-eb8d996c916b",
#                      "0bfebf9f-c779-e57d-e050-11ac0d487827",  "0bfd1043-7ebf-aaec-e050-11ac0c482f39",
#                      "c3cb64b3-842b-4b42-8c42-020a0206eb06", "0bfd1068-3fdd-a95b-e050-11ac0c4860c3",
#                      "fca1c663-dec0-db39-e040-11ac0c486140", "0bfd1068-3fdd-a95b-e050-11ac0c4860c3"
#                      
# )
# )
# 
# info2 <- fread('../extdata/rawInput/pcawg_sample_sheet.tsv')
# included <- as.data.frame(info2[which(info2$donor_wgs_exclusion_white_gray == "Whitelist"),
#                                 "aliquot_id"], row.names = c())
# 
# maf_liver <- maf_liver[maf_liver$Tumor_Sample_Barcode %in% included$aliquot_id,]
# 
# 
# length(unique(maf_liver$Donor_ID))
# length(unique(maf_liver$Tumor_Sample_Barcode))
# 
# maf_liver <- maf_liver[!maf_liver$Tumor_Sample_Barcode %in% removed_samples,]
# maf_liver <- maf_liver[!duplicated(maf_liver),]
# 
# length(unique(maf_liver$Donor_ID))
# length(unique(maf_liver$Tumor_Sample_Barcode))
# 
# 
# # info <- read_xlsx('../extdata/rawInput/PCAWG_supp/pcawg_donor_clinical_August2016_v9.xlsx')
# # included <- as.data.frame(info[which(info$donor_wgs_included_excluded == "Included"),
# #                                "icgc_donor_id"], row.names = c())
# # maf_liver <- maf_liver[maf_liver$Donor_ID %in% included$icgc_donor_id,]
# 
# 
# allSamp_projs1 <- data.frame(cbind(maf_liver$Tumor_Sample_Barcode, maf_liver$Donor_ID))
# allSamp_projs1 <- allSamp_projs1[!duplicated(allSamp_projs1),]
# duplicatedDonors <- unique(allSamp_projs1$X2[duplicated(allSamp_projs1$X2)])
# allSamp_projs1$dups <- ifelse(allSamp_projs1$X2 %in% duplicatedDonors, T, F)
# 
# do_multipleSampless <- (allSamp_projs1[allSamp_projs1$dups, c('X2')])
# do_multipleSampless <- allSamp_projs1[allSamp_projs1$X2 %in% do_multipleSampless, ]
# 
# do_UniqSampless <- allSamp_projs1[!allSamp_projs1$X2 %in% do_multipleSampless$X2,]
# 
# do_UniqSampless$dups <- FALSE
# do_multipleSampless$dups <- TRUE
# 
# library(data.table)
# library(rtracklayer)
# library(dplyr)
# 
# all_maf <- rbind(do_UniqSampless, do_multipleSampless)
# colnames(all_maf) <- c("Tumor_Sample_Barcode", "Donor_ID", "dups")
# CADDannotated_maf <- left_join(maf_liver, all_maf, by = c("Tumor_Sample_Barcode", "Donor_ID"))
# 
# gr <- GenomicRanges::GRanges(seqnames = paste0("chr", CADDannotated_maf$Chromosome),
#                              #strand = CADDannotated_maf$Strand,
#                              IRanges::IRanges(start = CADDannotated_maf$Start_position,
#                                               end = CADDannotated_maf$End_position),
#                              ref_al = CADDannotated_maf$Reference_Allele,
#                              alt_al = CADDannotated_maf$Tumor_Seq_Allele2,
#                              D_id = CADDannotated_maf$Donor_ID,
#                              sampleID = CADDannotated_maf$Tumor_Sample_Barcode,
#                              cohort = CADDannotated_maf$Project_Code,
#                              Ref_cntxt = CADDannotated_maf$ref_context,
#                              var_type = CADDannotated_maf$Variant_Type,
#                              variant_class = CADDannotated_maf$Variant_Classification,
#                              
#                              uniqueDonor_sample = CADDannotated_maf$dups)
# 
# strand(gr) <- Rle(rep("+", nrow(CADDannotated_maf)))
# 
# 
# 
# #
# 
# source('main/functions_iDriver.R')
# 
# 
# map_muts <- function(gr, path_bed){
#   
#   testGE_gr <- import.bed(path_bed)
#   
#   names(testGE_gr) = testGE_gr$name
#   
#   GenomicElement <- unlist(lapply(strsplit(names(testGE_gr), "[::]"), function(x){x[1]}))
#   
#   lenElement <- width(testGE_gr)
#   
#   mcols(testGE_gr) <- DataFrame(mcols(testGE_gr), GenomicElement, lenElement)
#   
#   callable_GE <- split(testGE_gr, names(testGE_gr))
#   lenElement_new <- sum(width(callable_GE))
#   
#   ov <- findOverlaps(gr, callable_GE, ignore.strand = TRUE)
#   
#   idx_gr <- queryHits(ov)
#   idx_callable_GE <- subjectHits(ov)
#   
#   gr_mappedGE <- gr[idx_gr]
#   
#   GE_col=unique(mcols(callable_GE, level = "within")[,"GenomicElement"])[idx_callable_GE]
#   EL_col = sum(mcols(callable_GE, level = "within")[,"lenElement"])[idx_callable_GE]
#   name_col = names(callable_GE)[idx_callable_GE]
#   
#   mcols(gr_mappedGE) <- DataFrame(mcols(gr_mappedGE),
#                                   GenomicElement=GE_col,
#                                   elemenntLength = EL_col,
#                                   name = name_col)
#   
#   idx_completely_black <- which(lenElement_new == 0)
#   
#   
#   
#   if(length(idx_completely_black) != 0){
#     lenElement_new <- lenElement_new[-idx_completely_black]
#     callable_GE <- callable_GE[-idx_completely_black]
#   }
#   
#   
#   list(MAF_GE_mapped = gr_mappedGE, all_GEs_grl = callable_GE,
#        lenAllTests = lenElement_new)
#   
# }
# 
# 
# create_response <- function(all_elements, n_donors){
#   
#   obs_Muts <- all_elements$MAF_GE_mapped
#   df <- as.data.frame(DataFrame(obs_Muts))
#   df <- df[,c('name', 'D_id')]
#   
#   
#   resTab <- df %>%
#     count(name, D_id) %>%
#     filter(n != 0) %>%
#     group_by(binID = name) %>%
#     summarise(nMut = sum(n), nSample = n(), .groups = "drop")
#   
#   idx_lenElement <- resTab$binID
#   all_lmnt_length <- all_elements$lenAllTests
#   ElementLength <- all_lmnt_length[idx_lenElement]
#   resTab$length <- ElementLength
#   
#   resTab$N <- n_donors
#   
#   all_binIds <- names(all_elements[["lenAllTests"]])
#   nonMut_ids <- all_binIds[which(!all_binIds %in% resTab$binID)]
#   res_nonMutated <- data.frame(cbind('binID' = nonMut_ids,
#                                      'nMut' = rep(0, length(nonMut_ids)),
#                                      'nSample' = rep(0, length(nonMut_ids)),
#                                      'length' = all_lmnt_length[!names(all_lmnt_length) %in% resTab$binID],
#                                      'N' = rep(n_donors, length(nonMut_ids))))
#   
#   resTab <- rbind(resTab, res_nonMutated)
#   
#   return(resTab)
# }
# 
# #
# #
# #
# path_bed <- '../../../../Projects/bahari_work/make_features/external/database/bins/proccessed/PCAWG_callable.bed6'
# all_elements <- map_muts(gr, path_bed)
# MAF_GE_mapped <- all_elements$MAF_GE_mapped
# # save(MAF_GE_mapped, file = '../extdata/procInput/mut_PCAWG/tmp_MAF_GE_mapped_allSamples.RData')
# 
# ##################################################################################
# # load('../extdata/procInput/mut_PCAWG/tmp_MAF_GE_mapped_allSamples.RData')
# 
# n_donors <- length(unique(gr$D_id))
# resTab <- create_response(all_elements, n_donors)
# 
# 
# ############################# select diffs id & dp ##########################
# 
# 
# id <- as.data.frame(resTab)
# 
# dp <- rbind(
#   fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Prost-AdenoCa/DriverPower/Prost-AdenoCA.3utr.DriverPower.observed.txt"
#   ),
#   fread(file.path("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Prost-AdenoCa/DriverPower/Prost-AdenoCA.5utr.DriverPower.observed.txt"
#   )),
#   fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Prost-AdenoCa/DriverPower/Prost-AdenoCA.CDS.DriverPower.observed.txt")
#   ,
#   fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Prost-AdenoCa/DriverPower/Prost-AdenoCA.enhancers.DriverPower.observed.txt")
#   ,
#   fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Prost-AdenoCa/DriverPower/Prost-AdenoCA.promCore.DriverPower.observed.txt")
#   
# )
# 
# ## Remove unwanted elements
# id <- id[!grepl("ncrna", id$binID),]
# id <- id[!grepl("gc19_pc.ss", id$binID),]
# 
# id <- left_join(id, dp, by = c('binID' = 'element_ID'  ))
# id <- id[!is.na(id$nSample.y),]
# 
# id$dif <- as.numeric(id$nSample.x) - as.numeric(id$nSample.y)
# 
# sum(id$dif !=0)
# id <- id[id$dif !=0,]
# rownames(id) <- 1:nrow(id)
# 
# get_sample <- function(MAF_GE_mapped, binID){
#   index <- which(MAF_GE_mapped$name %in% c(binID))
#   
#   elem_maf <- MAF_GE_mapped[index]
#   # elem_maf[which(!elem_maf$uniqueDonor_sample)]
#   print(unique(elem_maf$sampleID))
#   print(binID)
# }
# 
# # get_sample(MAF_GE_mapped, 'enhancers::chr11:87488200-87489800::NA::NA')
# 
# get_sample(MAF_GE_mapped, id$binID[9:12])
# 


# "0bfd1043-817e-e3e4-e050-11ac0c4860c5", "84aa2983-714b-49fd-af93-f18ed4d6a6cb",
# "63faa59e-6649-412e-823a-338d0662ec72", "fc9d704c-0cfb-0407-e040-11ac0d486980"
# "0bfd1043-8173-e3e4-e050-11ac0c4860c5", "0bfd1043-8180-e3e4-e050-11ac0c4860c5", "f1d73c2f-3d7d-4a0f-871b-b3bf76d865dc"
# "0bfd1068-3fc3-a95b-e050-11ac0c4860c3", "0bfe2ac9-0af3-c248-e050-11ac0d487e1c", "1558629a-770a-4dcb-b0bf-9a8bcea9d80b",
# "db37739d-af50-4e63-9147-1a7d25630ac6", "fca04929-586a-0594-e040-11ac0c48099d"










# rm(list = ls())
# library(data.table)
# library(readxl)
# 
# maf_liver <- fread('../extdata/procInput/mut_PCAWG/tmp_mafLIVER_allSamps.tsv')
# length(unique(maf_liver$Donor_ID))
# length(unique(maf_liver$Tumor_Sample_Barcode))
# 
# allSamp_projs0 <- data.frame(cbind(maf_liver$Tumor_Sample_Barcode, maf_liver$Donor_ID))
# allSamp_projs0 <- allSamp_projs0[!duplicated(allSamp_projs0),]
# allSamp_projs0$dups <- duplicated(allSamp_projs0$X2)
# 
# 
# maf_liver <- maf_liver[!maf_liver$Tumor_Sample_Barcode %in% c('26de406e-c623-11e3-bf01-24c6515278c0',
# 
# 
#                                                               "2b569b88-c622-11e3-bf01-24c6515278c0",
#                                                               "632f3eb0-c623-11e3-bf01-24c6515278c0",
# 
#                                                               "3f99ae0e-c623-11e3-bf01-24c6515278c0",
#                                                               "72c22a8c-c622-11e3-bf01-24c6515278c0",
#                                                               "10db2142-c623-11e3-bf01-24c6515278c0",
#                                                               "072f4146-c623-11e3-bf01-24c6515278c0",
#                                                               "690e481c-c623-11e3-bf01-24c6515278c0",
# 
#                                                               # "0a9c9db0-c623-11e3-bf01-24c6515278c0",
#                                                               "6ce66be0-c623-11e3-bf01-24c6515278c0",
#                                                               "520b8ab2-c623-11e3-bf01-24c6515278c0"),]
# maf_liver <- maf_liver[!duplicated(maf_liver),]
# 
# # dup_DOs >>> "DO45273" "DO45096" "DO50793" "DO50814" "DO45141" "DO50842" "DO50844" "DO45151" "DO23543"
# dup_DOs <- c( "DO45273", "DO45096" ,"DO50793", "DO50814", "DO45141", "DO50842" ,"DO50844", "DO45151", "DO23543")
# allSamp_projs0$dupDOs <- ifelse(allSamp_projs0$X2 %in% dup_DOs, T, F)
# allSamp_projs0$removed <- allSamp_projs0$X1 %in% c( '26de406e-c623-11e3-bf01-24c6515278c0',
# 
# 
#                                                    "2b569b88-c622-11e3-bf01-24c6515278c0",
#                                                    "632f3eb0-c623-11e3-bf01-24c6515278c0",
# 
#                                                    "3f99ae0e-c623-11e3-bf01-24c6515278c0",
#                                                    "72c22a8c-c622-11e3-bf01-24c6515278c0",
#                                                    "10db2142-c623-11e3-bf01-24c6515278c0",
#                                                    "072f4146-c623-11e3-bf01-24c6515278c0",
#                                                    "690e481c-c623-11e3-bf01-24c6515278c0",
# 
#                                                    # "0a9c9db0-c623-11e3-bf01-24c6515278c0",
#                                                    "6ce66be0-c623-11e3-bf01-24c6515278c0",
#                                                    "520b8ab2-c623-11e3-bf01-24c6515278c0"
#                                                    )
# 
# # info <- read_xlsx('../extdata/rawInput/PCAWG_supp/pcawg_donor_clinical_August2016_v9.xlsx')
# # included <- as.data.frame(info[which(info$donor_wgs_included_excluded == "Included"),
# #                                "icgc_donor_id"], row.names = c())
# # maf_liver <- maf_liver[maf_liver$Donor_ID %in% included$icgc_donor_id,]
# 
# info2 <- fread('../extdata/rawInput/pcawg_sample_sheet.tsv')
# included <- as.data.frame(info2[which(info2$donor_wgs_exclusion_white_gray == "Whitelist"),
#                                "aliquot_id"], row.names = c())
# 
# maf_liver <- maf_liver[maf_liver$Tumor_Sample_Barcode %in% included$aliquot_id,]
# 
# 
# length(unique(maf_liver$Donor_ID))
# length(unique(maf_liver$Tumor_Sample_Barcode))
# 
# allSamp_projs1 <- data.frame(cbind(maf_liver$Tumor_Sample_Barcode, maf_liver$Donor_ID))
# allSamp_projs1 <- allSamp_projs1[!duplicated(allSamp_projs1),]
# allSamp_projs1$dups <- duplicated(allSamp_projs1$X2)
# 
# unique(allSamp_projs1$X2)
# 
# do_multipleSampless <- (allSamp_projs1[allSamp_projs1$dups, c('X2')])
# do_multipleSampless <- allSamp_projs1[allSamp_projs1$X2 %in% do_multipleSampless, ]
# 
# do_UniqSampless <- allSamp_projs1[!allSamp_projs1$X2 %in% do_multipleSampless$X2,]
# 
# do_UniqSampless$dups <- FALSE
# do_multipleSampless$dups <- TRUE
# 
# all_maf <- rbind(do_UniqSampless, do_multipleSampless)
# colnames(all_maf) <- c("Tumor_Sample_Barcode", "Donor_ID", "dups")
# CADDannotated_maf <- left_join(maf_liver, all_maf, by = c("Tumor_Sample_Barcode", "Donor_ID"))
# 
# gr <- GenomicRanges::GRanges(seqnames = paste0("chr", CADDannotated_maf$Chromosome),
#                              #strand = CADDannotated_maf$Strand,
#                              IRanges::IRanges(start = CADDannotated_maf$Start_position,
#                                               end = CADDannotated_maf$End_position),
#                              ref_al = CADDannotated_maf$Reference_Allele,
#                              alt_al = CADDannotated_maf$Tumor_Seq_Allele2,
#                              D_id = CADDannotated_maf$Donor_ID,
#                              sampleID = CADDannotated_maf$Tumor_Sample_Barcode,
#                              cohort = CADDannotated_maf$Project_Code,
#                              Ref_cntxt = CADDannotated_maf$ref_context,
#                              var_type = CADDannotated_maf$Variant_Type,
#                              variant_class = CADDannotated_maf$Variant_Classification,
# 
#                              uniqueDonor_sample = CADDannotated_maf$dups)
# 
# strand(gr) <- Rle(rep("+", nrow(CADDannotated_maf)))
# 
# 
# 
# #
# library(data.table)
# library(rtracklayer)
# library(dplyr)
# source('main/functions_iDriver.R')
# 
# 
# map_muts <- function(gr, path_bed){
# 
#   testGE_gr <- import.bed(path_bed)
# 
#   names(testGE_gr) = testGE_gr$name
# 
#   GenomicElement <- unlist(lapply(strsplit(names(testGE_gr), "[::]"), function(x){x[1]}))
# 
#   lenElement <- width(testGE_gr)
# 
#   mcols(testGE_gr) <- DataFrame(mcols(testGE_gr), GenomicElement, lenElement)
# 
#   callable_GE <- split(testGE_gr, names(testGE_gr))
#   lenElement_new <- sum(width(callable_GE))
# 
#   ov <- findOverlaps(gr, callable_GE, ignore.strand = TRUE)
# 
#   idx_gr <- queryHits(ov)
#   idx_callable_GE <- subjectHits(ov)
# 
#   gr_mappedGE <- gr[idx_gr]
# 
#   GE_col=unique(mcols(callable_GE, level = "within")[,"GenomicElement"])[idx_callable_GE]
#   EL_col = sum(mcols(callable_GE, level = "within")[,"lenElement"])[idx_callable_GE]
#   name_col = names(callable_GE)[idx_callable_GE]
# 
#   mcols(gr_mappedGE) <- DataFrame(mcols(gr_mappedGE),
#                                   GenomicElement=GE_col,
#                                   elemenntLength = EL_col,
#                                   name = name_col)
# 
#   idx_completely_black <- which(lenElement_new == 0)
# 
# 
# 
#   if(length(idx_completely_black) != 0){
#     lenElement_new <- lenElement_new[-idx_completely_black]
#     callable_GE <- callable_GE[-idx_completely_black]
#   }
# 
# 
#   list(MAF_GE_mapped = gr_mappedGE, all_GEs_grl = callable_GE,
#        lenAllTests = lenElement_new)
# 
# }
# 
# 
# create_response <- function(all_elements, n_donors){
# 
#   obs_Muts <- all_elements$MAF_GE_mapped
#   df <- as.data.frame(DataFrame(obs_Muts))
#   df <- df[,c('name', 'D_id')]
# 
# 
#   resTab <- df %>%
#     count(name, D_id) %>%
#     filter(n != 0) %>%
#     group_by(binID = name) %>%
#     summarise(nMut = sum(n), nSample = n(), .groups = "drop")
# 
#   idx_lenElement <- resTab$binID
#   all_lmnt_length <- all_elements$lenAllTests
#   ElementLength <- all_lmnt_length[idx_lenElement]
#   resTab$length <- ElementLength
# 
#   resTab$N <- n_donors
# 
#   all_binIds <- names(all_elements[["lenAllTests"]])
#   nonMut_ids <- all_binIds[which(!all_binIds %in% resTab$binID)]
#   res_nonMutated <- data.frame(cbind('binID' = nonMut_ids,
#                                      'nMut' = rep(0, length(nonMut_ids)),
#                                      'nSample' = rep(0, length(nonMut_ids)),
#                                      'length' = all_lmnt_length[!names(all_lmnt_length) %in% resTab$binID],
#                                      'N' = rep(n_donors, length(nonMut_ids))))
# 
#   resTab <- rbind(resTab, res_nonMutated)
# 
#   return(resTab)
# }
# 
# #
# #
# #
# path_bed <- '../../../../Projects/bahari_work/make_features/external/database/bins/proccessed/PCAWG_callable.bed6'
# all_elements <- map_muts(gr, path_bed)
# MAF_GE_mapped <- all_elements$MAF_GE_mapped
# # save(MAF_GE_mapped, file = '../extdata/procInput/mut_PCAWG/tmp_MAF_GE_mapped_allSamples.RData')
# 
# ##################################################################################
# # load('../extdata/procInput/mut_PCAWG/tmp_MAF_GE_mapped_allSamples.RData')
# 
# n_donors <- length(unique(gr$D_id))
# resTab <- create_response(all_elements, n_donors)
# 
# 
# ############################# select diffs id & dp ##########################
# 
# 
# id <- as.data.frame(resTab)
# 
# dp <- rbind(
#   fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Liver-HCC/DriverPower/Liver-HCC.3utr.DriverPower.observed.txt"
#   ),
#   fread(file.path("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Liver-HCC/DriverPower/Liver-HCC.5utr.DriverPower.observed.txt"
#   )),
#   fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Liver-HCC/DriverPower/Liver-HCC.CDS.DriverPower.observed.txt")
#   ,
#   fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Liver-HCC/DriverPower/Liver-HCC.enhancers.DriverPower.observed.txt")
#   ,
#   fread("../extdata/rawInput/PCAWG_12methods_output/p-values/observed/Liver-HCC/DriverPower/Liver-HCC.promCore.DriverPower.observed.txt")
# 
# )
# 
# ## Remove unwanted elements
# id <- id[!grepl("ncrna", id$binID),]
# id <- id[!grepl("gc19_pc.ss", id$binID),]
# 
# id <- left_join(id, dp, by = c('binID' = 'element_ID'  ))
# id <- id[!is.na(id$nSample.y),]
# 
# id$dif <- as.numeric(id$nSample.x) - as.numeric(id$nSample.y)
# 
# sum(id$dif !=0)
# id <- id[id$dif !=0,]
# 
# index <- which(MAF_GE_mapped$name %in% c('enhancers::chr10:69859400-69861200::NA::NA'))
# 
# elem_maf <- MAF_GE_mapped[index]
# # elem_maf[which(!elem_maf$uniqueDonor_sample)]
# unique(elem_maf$sampleID)
# 
# 
# DoIDs <- unique(allSamp_projs0$X2)
# DoIDs[which(!DoIDs %in% unique(maf_liver$Donor_ID))]


