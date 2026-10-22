rm(list = ls())

maf <- fread('../extdata/procInput/mut_PCAWG/tmp_mafLIVER_allSamps.tsv')

liver_remove <- c('26de406e-c623-11e3-bf01-24c6515278c0'
                  # ,
                  # 
                  # 
                  # "520b8ab2-c623-11e3-bf01-24c6515278c0",
                  # "062e96d4-c623-11e3-bf01-24c6515278c0",
                  # "10db2142-c623-11e3-bf01-24c6515278c0",
                  # "27fcccdc-c622-11e3-bf01-24c6515278c0",
                  # 
                  # "3f99ae0e-c623-11e3-bf01-24c6515278c0",
                  # 
                  # "632f3eb0-c623-11e3-bf01-24c6515278c0",
                  # "690e481c-c623-11e3-bf01-24c6515278c0",
                  # 
                  # "6c306cbe-c623-11e3-bf01-24c6515278c0",
                  # "70c2f91e-c622-11e3-bf01-24c6515278c0"
                  )

maf_subset <- maf[!maf$Tumor_Sample_Barcode %in% liver_remove, c( "Chromosome",  "Start_position", "End_position",
                                                                  "Reference_Allele", "Tumor_Seq_Allele2",
                                                                  "Donor_ID", "Tumor_Sample_Barcode" ,
                                                                  "ref_context", "Variant_Type",
                                                                  "Project_Code", "Variant_Classification")]
