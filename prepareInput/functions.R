

# get_OneSamp_perOneDon <- function(maf, samp_info, Stars_value){
#   
#   xxx = maf[,c('Donor_ID', 'Tumor_Sample_Barcode', 'Project_Code')]
#   xxx = xxx[!duplicated(xxx),]
#   df=as.data.frame(table(xxx$Donor_ID))
#   df=df[df$Freq >1,]
#   
#   n = xxx[which(xxx$Donor_ID %in% df$Var1),]
#   n = xxx[!xxx$Tumor_Sample_Barcode %in% c('26de406e-c623-11e3-bf01-24c6515278c0'),] # biliary sample for DO45281
#   
#   # liver_Biliary <- n[n$Tumor_Sample_Barcode == '2109e1d4-c623-11e3-bf01-24c6515278c0',] # liver sample for DO45281
#   # 
#   samp2 <- Stars_value
#   samp2 <- samp2[,c("Stars", "tumor_wgs_icgc_specimen_id", "icgc_donor_id")]
#   
#   samp_info1 <- samp_info[samp_info$aliquot_id %in% n$Tumor_Sample_Barcode,]
#   samp_info1 <- samp_info1[
#     ,
#     if (any(grepl("Primary tumour", dcc_specimen_type, ignore.case = TRUE))) {
#       .SD[grepl("Primary tumour", dcc_specimen_type, ignore.case = TRUE)]
#     } else {
#       .SD
#     },
#     by = icgc_donor_id
#     ]
#   
#   samp_info1 <- left_join(samp_info1, n, by = c('aliquot_id' = 'Tumor_Sample_Barcode'))
#   samp_info1 <- left_join(samp_info1, samp2, by = c("icgc_specimen_id" = "tumor_wgs_icgc_specimen_id"))
#   
#   
#   samp_highest_star <- samp_info1[
#     order(-Stars),
#     .SD[1],
#     by = Donor_ID
#     ]
#   
#   samp_highest_star[,c("Donor_ID", "aliquot_id"  )]
# }

remove_samples <- function(){
  
  removed_MDS <- c("f8e61a02-8fbe-a9af-e040-11ac0d481b64", "f90a074e-d82f-f6bb-e040-11ac0d487334"
                   )
  
  removed_samples_PancAdeno <- c( "0efee3f3-313b-44eb-b3b0-af975a1ded82",
                                  "60c33e32-7e19-4e71-b075-a63fcf27e660")
  removed_samples_AML <- c("f8e61a02-654c-c226-e040-11ac0d481b60", "f8e61a02-8fbe-a9af-e040-11ac0d481b64",
                           "f8e61a02-90ec-96fd-e040-11ac0d481b5e")
  
  removed_Liver <- c('26de406e-c623-11e3-bf01-24c6515278c0',
                     
                     
                     "2b569b88-c622-11e3-bf01-24c6515278c0", 
                     "632f3eb0-c623-11e3-bf01-24c6515278c0",
                     
                     "3f99ae0e-c623-11e3-bf01-24c6515278c0",
                     "72c22a8c-c622-11e3-bf01-24c6515278c0",
                     "10db2142-c623-11e3-bf01-24c6515278c0",
                     "072f4146-c623-11e3-bf01-24c6515278c0",
                     "690e481c-c623-11e3-bf01-24c6515278c0",
                     
                     # "0a9c9db0-c623-11e3-bf01-24c6515278c0",
                     "6ce66be0-c623-11e3-bf01-24c6515278c0",
                     "520b8ab2-c623-11e3-bf01-24c6515278c0")
  
  removed_samples_MyeloidMPN <- c("f8e61a02-8c03-ad45-e040-11ac0d481b62", "f8e61a02-8c9e-aaee-e040-11ac0d481b6a", "f8e61a02-938b-9428-e040-11ac0d481b6c",
                       "f8efa5ad-3a4d-39d6-e040-11ac0c480132", "f8f7f274-dd98-3bb0-e040-11ac0c483fcd", "f8f9d9b8-abdd-1afb-e040-11ac0d480523",
                       "f9012277-55c4-2486-e040-11ac0c48236f", "f9020478-c741-6776-e040-11ac0d481228", "f903cfe2-c082-1bf6-e040-11ac0d4809f7",
                       "f9052bfc-94a0-cb29-e040-11ac0c48284a", "f90b298a-0546-6c36-e040-11ac0c487a0b", "f90dbdbf-f8cd-e0d8-e040-11ac0c482576",
                       "f9101872-d864-c508-e040-11ac0c4818b5", "f91661d6-a2b5-fac2-e040-11ac0d4852fd", "f91a9af2-7909-099e-e040-11ac0c484205",
                       "f91faa2e-ea53-bef4-e040-11ac0d485f5d" ,"f925f24b-bffd-89ab-e040-11ac0d48682f" ,"f92b2c52-ed64-4594-e040-11ac0d48425c",
                       "f92c8821-daba-3727-e040-11ac0d480224", "f939f757-a415-ad44-e040-11ac0d4868e5", "f93c9e35-c76c-2b59-e040-11ac0d485dee",
                       "f93dd46d-4611-e84e-e040-11ac0c483323",
                       "f8e61a02-8c03-ad45-e040-11ac0d481b62", "f8e61a02-8c9e-aaee-e040-11ac0d481b6a", "f8e61a02-938b-9428-e040-11ac0d481b6c",
                       "f8efa5ad-3a4d-39d6-e040-11ac0c480132", "f8f7f274-dd98-3bb0-e040-11ac0c483fcd", "f8f9d9b8-abdd-1afb-e040-11ac0d480523",
                       "f9012277-55c4-2486-e040-11ac0c48236f", "f9020478-c741-6776-e040-11ac0d481228" ,"f903cfe2-c082-1bf6-e040-11ac0d4809f7",
                       "f9052bfc-94a0-cb29-e040-11ac0c48284a", "f90b298a-0546-6c36-e040-11ac0c487a0b", "f90dbdbf-f8cd-e0d8-e040-11ac0c482576",
                       "f9101872-d864-c508-e040-11ac0c4818b5", "f91661d6-a2b5-fac2-e040-11ac0d4852fd", "f91a9af2-7909-099e-e040-11ac0c484205",
                       "f91faa2e-ea53-bef4-e040-11ac0d485f5d", "f925f24b-bffd-89ab-e040-11ac0d48682f", "f92b2c52-ed64-4594-e040-11ac0d48425c",
                       "f92c8821-daba-3727-e040-11ac0d480224", "f939f757-a415-ad44-e040-11ac0d4868e5", "f93c9e35-c76c-2b59-e040-11ac0d485dee",
                       "f93dd46d-4611-e84e-e040-11ac0c483323",
                       
                       
                       "f8efc19a-1e4b-e1b0-e040-11ac0c48104e", "f8f024e5-7096-3047-e040-11ac0d481c0f",
                       "f91d1645-9b68-edbb-e040-11ac0c4854d8",
                       "f920ef57-c942-d2e1-e040-11ac0d481706", "f937fc8b-6fa7-12fd-e040-11ac0d4850f9",
                       "f942b732-7c0b-7ec6-e040-11ac0c483f86")
  
  removed_samples_Prostat <- unique(c("0bfd1068-3fca-a95b-e050-11ac0c4860c3", "5fe8a231-bbac-4f91-a47a-0de4fee752d3",
                                      "0bfd1043-8189-e3e4-e050-11ac0c4860c5", "1d6c9181-a02f-4c35-85b3-65dc5808f50e",
                                      "95730d8b-327a-402c-a8bf-d71370f777a7", 
                                      "2b54d668-7585-4bf4-9fdd-aea3c5cd4fa4", "0bfd1043-7ec1-aaec-e050-11ac0c482f39",
                                      "2b54d668-7585-4bf4-9fdd-aea3c5cd4fa4", "2b54d668-7585-4bf4-9fdd-aea3c5cd4fa4",
                                      
                                      
                                      "2b54d668-7585-4bf4-9fdd-aea3c5cd4fa4",
                                      "5df0aab4-0be2-4da3-a4b7-a10fca09a4ba",
                                      "74db46a8-a499-4951-af33-1d53b7d1c9b4",
                                      "7e10e325-16d2-42df-b957-df4d3a9477d0",
                                      
                                      "052015bc-8c24-467c-8ff7-9e0dc4ad39fd", "0bd3a230-531e-44aa-9999-b8ed8da0176b",
                                      "0bfd1068-3fd5-a95b-e050-11ac0c4860c3",
                                      "0bfe2ac9-0af8-c248-e050-11ac0d487e1c", "0c0038ff-6cc4-b0b0-e050-11ac0d483d73",
                                      "0c0038ff-6cc6-b0b0-e050-11ac0d483d73",
                                      "339a2c94-55c9-4113-a706-b0971fa2784f",
                                      
                                      "4c755527-3b5d-4d36-822d-990aa1003d6a", "0bfd1043-7346-fdd0-e050-11ac0c484cab",
                                      "2f79fe89-9d24-4fbe-ba5b-0ff6c3c4d663", "fd163b56-1c2d-4e63-9fad-ba3221ae274b",
                                      "fd163b56-1c2d-4e63-9fad-ba3221ae274b",
                                      "3b8ab548-6669-4142-9612-cafded0d7a23" ,
                                      "3b8ab548-6669-4142-9612-cafded0d7a23", "fc9f9295-01f6-a9e0-e040-11ac0d4845af",
                                      "3b8ab548-6669-4142-9612-cafded0d7a23", "fc9f9295-01f6-a9e0-e040-11ac0d4845af",
                                      
                                      
                                      "0bfe2ac9-0afa-c248-e050-11ac0d487e1c", "0bfe2ac9-0aee-c248-e050-11ac0d487e1c",
                                      "fc9d93b6-77c6-8719-e040-11ac0d487dec", "0bfe2ac9-0afa-c248-e050-11ac0d487e1c",
                                      "0bfe2ac9-0aec-c248-e050-11ac0d487e1c", "0bfd1068-3fe6-a95b-e050-11ac0c4860c3",
                                      
                                      "0bfd1043-5142-3662-e050-11ac0c486501", "0bfd1068-3fdf-a95b-e050-11ac0c4860c3",
                                      "0c0038ff-6cc8-b0b0-e050-11ac0d483d73", "0bfebf9f-c781-e57d-e050-11ac0d487827",
                                      "0bfebf9f-c77d-e57d-e050-11ac0d487827", "fc9d5fc7-a87b-c5e7-e040-11ac0c4850aa",
                                      "b0f3087b-6f9c-47e9-b487-6726af1cb90b", "280a069c-38e5-4300-a15c-2612301b70e3",
                                      "7a86abdc-d554-43bd-ad47-22d439c930ef", "0bfe2ac9-0afd-c248-e050-11ac0d487e1c",
                                      
                                      "0bfe2ac9-0af5-c248-e050-11ac0d487e1c", "0bfd1043-7343-fdd0-e050-11ac0c484cab",
                                      "0bfe2ac9-0aff-c248-e050-11ac0d487e1c", "0bfd1043-817c-e3e4-e050-11ac0c4860c5",
                                      "fc9d5fc6-546f-27b7-e040-11ac0d485f87",
                                      
                                      "0bfd1043-817e-e3e4-e050-11ac0c4860c5", "84aa2983-714b-49fd-af93-f18ed4d6a6cb",
                                      "63faa59e-6649-412e-823a-338d0662ec72", "fc9d704c-0cfb-0407-e040-11ac0d486980",
                                      "0bfd1043-8173-e3e4-e050-11ac0c4860c5", "0bfd1043-8180-e3e4-e050-11ac0c4860c5",
                                      "f1d73c2f-3d7d-4a0f-871b-b3bf76d865dc", "0bfd1068-3fc3-a95b-e050-11ac0c4860c3",
                                      "0bfe2ac9-0af3-c248-e050-11ac0d487e1c", "1558629a-770a-4dcb-b0bf-9a8bcea9d80b",
                                      "db37739d-af50-4e63-9147-1a7d25630ac6", "fca04929-586a-0594-e040-11ac0c48099d",
                                      "0bfd1043-816e-e3e4-e050-11ac0c4860c5", "fca2d62f-2756-939e-e040-11ac0d4865e6",
                                      "0bfd1068-3fcd-a95b-e050-11ac0c4860c3", "0bfd1068-3fe1-a95b-e050-11ac0c4860c3",
                                      "0bfebf9f-c783-e57d-e050-11ac0d487827", "0bfebf9f-c77f-e57d-e050-11ac0d487827",
                                      
                                      "0bfd1043-7ed6-9ccc-e050-11ac0c481957", "9cdb1905-8b0f-48bd-ba37-11e7ac9646a0",
                                      "fc9f998b-cc72-f6ea-e040-11ac0c482425", "9b75735b-705d-4ebf-92a0-fbb190b95285",
                                      "0bfd1043-7344-fdd0-e050-11ac0c484cab", "6ff16b0d-5966-4119-9bec-c9ebff804a37",
                                      
                                      "0bfd1068-3fe4-a95b-e050-11ac0c4860c3", "0bfd1043-8172-e3e4-e050-11ac0c4860c5",
                                      "0bfd1043-7344-fdd0-e050-11ac0c484cab",
                                      "0bfd1043-8175-e3e4-e050-11ac0c4860c5", "1bd47e40-d708-4ca2-b4b3-eb8d996c916b",
                                      "0bfebf9f-c779-e57d-e050-11ac0d487827",  "0bfd1043-7ebf-aaec-e050-11ac0c482f39",
                                      "c3cb64b3-842b-4b42-8c42-020a0206eb06", "0bfd1068-3fdd-a95b-e050-11ac0c4860c3",
                                      "fca1c663-dec0-db39-e040-11ac0c486140", "0bfd1068-3fdd-a95b-e050-11ac0c4860c3"))
  
  c(removed_samples_MyeloidMPN, removed_Liver, removed_samples_Prostat,
    removed_samples_AML, removed_samples_PancAdeno, removed_MDS)
  
}

create_mut_df <- function(path_raw_tcga, path_raw_icgc, path_save, path_pcawg_sampleSheet,
                          overwrite = TRUE){
  
  
  dir.create(path_save, recursive = T, showWarnings = F)
  fname <- paste0(path_save, "MAF_samples_removed.tsv")
  
  if(overwrite == TRUE){
    
    maf_tcga <- fread(path_raw_tcga)    #dim >>> 29504368       43
    maf_icgc <- fread(path_raw_icgc)    #dim >>>  23159591       43
    
    maf <- rbind(maf_tcga, maf_icgc)  # dim >>>  52663959       43
    maf <- as.data.frame(maf)
    print(dim(maf))
    
    removed_samples <-remove_samples()
    
    maf_new <- maf[which(!maf$Tumor_Sample_Barcode %in% removed_samples), ] 
    
    df=maf_new[,c( "Chromosome",  "Start_position", "End_position",
                   "Reference_Allele", "Tumor_Seq_Allele2",
                   "Donor_ID",
                   "Tumor_Sample_Barcode" ,
                   "ref_context", "Variant_Type",
                   "Project_Code", "Variant_Classification")]  
    print(dim(df))
    ###### just include the white list samples (just retain included donors) ###########
    df <- retain_included_Donors(df, path_pcawg_sampleSheet) 
    print(dim(df))
    df = df[!duplicated(df),] #  
    
    fwrite(df, file = fname, sep = "\t")
  } else {
    print("mutation file exists")
  }
}

add_CADD_2MAF <- function(path_SNV_offLine, path_indel_offLine, 
                            path_onLineSet, path_complete_MAF){
  all_MAF <- fread(path_complete_MAF)
  
  
  all_scores <- aggregate_Scores(path_SNV_offLine, path_indel_offLine, path_onLineSet)
  all_scores <- all_scores[,-c("ID")]
  all_scores <- all_scores[!duplicated(all_scores),]
  
  all_MAF <- all_MAF[, c("Chromosome", "Start_position", "End_position",
                    "Reference_Allele", "Tumor_Seq_Allele2", "Donor_ID" ,
                    "Tumor_Sample_Barcode",
                    "ref_context",  "Variant_Type",
                    "Project_Code" , "Variant_Classification")] # 46619433
  
  all_MAF <- all_MAF[!duplicated(all_MAF),] # 46611594
  
  all_MAF <- define_cadd_vcf_cols(all_MAF)
  MAF_annotated=left_join(all_MAF, all_scores, by = c("Chrom", "Pos", "Ref", "Alt"))
  MAF_annotated
}

retain_included_Donors <- function(reannotated_MAF, path_pcawgSupp){
  
  info2 <- fread(path_pcawgSupp)
  included <- as.data.frame(info2[which(info2$donor_wgs_exclusion_white_gray == "Whitelist"),
                                  "aliquot_id"], row.names = c())
  
  reannotated_MAF <- reannotated_MAF[reannotated_MAF$Tumor_Sample_Barcode %in% included$aliquot_id,]
  
  reannotated_MAF
}


add_ranked_scores <- function(df){
  
  df <- df[order(as.numeric(df$RawScore)),]
  
  Rank_portion <- (1:nrow(df))/nrow(df)
  df$Rank_portion <- ifelse(is.na(df$RawScore), NA, Rank_portion)
  df
  
}

################## create GRanges ###############
########## MAF_Obj2GRanges function
# Arguments: @CADDannotated_maf
MutFile2GRanges <- function(CADDannotated_maf){
  
  gr <- GenomicRanges::GRanges(seqnames = paste0("chr", CADDannotated_maf$Chrom), 
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
                               variant_class = CADDannotated_maf$Variant_Classification,
                               CADD_RawScore = CADDannotated_maf$RawScore,
                               CADD_PHRED = CADDannotated_maf$PHRED,
                               CADD_ranked = CADDannotated_maf$Rank_portion)
  
  strand(gr) <- Rle(rep("+", nrow(CADDannotated_maf)))
  
  gr
}


filter_samples <- function(path_pcawgSupp, gr,
                           exclude_lymph_melanoma = TRUE,
                           exclude_hyper_mutated = TRUE){
  
  # info <- read_xlsx(path_pcawgSupp)
  # 
  # included <- as.data.frame(info[which(info$donor_wgs_included_excluded == "Included"),
  #                                "icgc_donor_id"], row.names = c())
  # 
  # gr <- gr[which(gr$D_id %in% included[,1])]
  
  if (exclude_lymph_melanoma) {
    exceptions <- c("Skin-Melanoma", "Lymph-NOS", "Lymph-CLL", "Lymph-BNHL")
    gr <- gr[-which(gr$cohort %in% exceptions)]
  }
  
  if (exclude_hyper_mutated) {
    df <- data.frame(table(gr$cohort, gr$D_id))
    colnames(df) <- c("cohort", "D_id", "freq")
    df <- df[df$freq !=0,]
    
    nonHyperMut_donors <- unique(df[-which(df$freq/3000 >30), "D_id"])
    
    gr <- gr[which(gr$D_id %in% nonHyperMut_donors)] 
  }
  
  gr
}
