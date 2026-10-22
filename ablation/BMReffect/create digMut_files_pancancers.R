rm(list = ls())

source('main/functions_iDriver.R')

library(data.table)
library(dplyr)
library(BSgenome.Hsapiens.UCSC.hg19)
genome <- BSgenome.Hsapiens.UCSC.hg19


mutData <- fread('../extdata/procInput/iDriverInputs/mutData.tsv')
path_donorInfo <- '../extdata/procInput/iDriverInputs/donorInfo_withSamples.tsv'
x <- fread(path_donorInfo)

cohort <- 'Pancan-no-skin-melanoma-lymph'
# cohort <- 'Pan_Cancer'
exclude_lymph_melanoma <- ifelse(cohort == 'Pancan-no-skin-melanoma-lymph', T, F)

donorInfo <- select_cohort(path_donorInfo,
                           cohort,
                           exclude_lymph_melanoma,
                           exclude_hyper_mutated = F)

mutData <- mutData[mutData$donor_id %in% donorInfo$donor_id,]
# testGE <- fread('../extdata/rawInput/PCAWG_test_genomic_elements.bed12.gz')
# testGE <- testGE[, c(1,4)]
# mutData <- left_join(mutData, testGE, by = c('PCAWG_ID' = 'V4'))

load('../extdata/procInput/mut_PCAWG/all_samples.RData')
gr$var_ID <- paste0('M_', 1:length(gr))

gr_cohort <- gr[gr$D_id %in% donorInfo$donor_id]
gr_cohort <- gr_cohort[!seqnames(gr_cohort) %in% c('chrX', 'chrY')]


get_ref_context <- function(genome, chr, position, context_range = 10){
  getSeq(genome, chr, start = position - context_range, end = position + context_range)
}

extract_HugoSymbol <- function(binID_vector){
  
  s <- strsplit(binID_vector, "[::]")
  GenomicElement <- unlist(lapply(s, function(x){x[5]}))
  GenomicElement
  
}


dig_gr <- data.frame(cbind(as.character(seqnames(gr_cohort)), start(gr_cohort)-1,
                           end(gr_cohort), gr_cohort$ref_al,
                           gr_cohort$alt_al, #gr_cohort$var_ID,
                           gr_cohort$sampleID, gr_cohort$variant_class,
                           ifelse(gr_cohort$var_type %in% c('SNP', 'SNV') , paste0(gr_cohort$ref_al, '>', gr_cohort$alt_al),
                                  '.'),
                           ifelse(gr_cohort$var_type %in% c('SNP', 'SNV') ,
                                  as.character(get_ref_context(genome, as.character(seqnames(gr_cohort)),
                                                               start(gr_cohort),
                                                               context_range = 1))
                                  , '.'),
                           ifelse(!gr_cohort$var_type %in% c('SNP', 'SNV') , 'INDEL', 'SNP'),
                           gr_cohort$var_ID
                           
                           
                           ))

dig_gr$X1 <- sub("^chr", "", dig_gr$X1, ignore.case = TRUE)
dig_gr <- left_join(dig_gr, mutData, by = c('X11'= 'var_ID'))


# dig_gr$gene <- ifelse(is.na(dig_gr$PCAWG_ID) == '', '.', extract_HugoSymbol(dig_gr$PCAWG_ID))
# dig_gr$gene <- ifelse(dig_gr$gene=='', '.', dig_gr$gene)
# dig_gr$gene <- ifelse(is.na(dig_gr$gene), '.', dig_gr$gene)

dig_gr$gene <- ifelse(dig_gr$GenomicElement == 'gc19_pc.cds', extract_HugoSymbol(dig_gr$PCAWG_ID), '.')
dig_gr$gene <- ifelse(is.na(dig_gr$PCAWG_ID) == '', '.', dig_gr$gene )
dig_gr$gene <- ifelse(is.na(dig_gr$gene), '.', dig_gr$gene )

dig_gr$V8 <- ifelse(dig_gr$gene == '.', 'Noncoding', dig_gr$X7)
dig_gr$V8 <- ifelse(dig_gr$X10 == 'INDEL', 'INDEL', dig_gr$V8)
dig_gr$V8 <- ifelse(is.na(dig_gr$V8 ), 'Noncoding', dig_gr$V8)
dig_gr$V8 <- ifelse(dig_gr$V8 == 'Silent', 'Synonymous', dig_gr$V8)
dig_gr$V8 <- ifelse(dig_gr$V8 == 'Missense_Mutation', 'Missense', dig_gr$V8)
dig_gr$V8 <- ifelse(dig_gr$V8 == 'Nonsense_Mutation', 'Nonsense', dig_gr$V8)

dig_input <- dig_gr[,c('X1', 'X2', 'X3', 'X4', 'X5', 'X6', 'gene', 'V8', 'X8', 'X9')]


# head(WHYna)
# donor_id totalMutNrs                                             PCAWG_ID
# 3    DO1000        2252       gc19_pc.3utr::gencode::AFF2::ENSG00000155966.9
# 10   DO1000        2252 gc19_pc.promCore::gencode::PAGE2B::ENSG00000238269.4
# 13   DO1000        2252         gc19_pc.3utr::gencode::F9::ENSG00000101981.6
# 39   DO1000        2252     gc19_pc.3utr::gencode::SYTL4::ENSG00000102362.11
# 61   DO1000        2252      gc19_pc.5utr::gencode::RGAG4::ENSG00000242732.3
# 62   DO1000        2252  gc19_pc.promCore::gencode::RGAG4::ENSG00000242732.3
# Ref Alt   triCntx var_type   GenomicElement  position     var_ID
# 3  CCT CAT CCT > CAT      SNP     gc19_pc.3utr 148081160  M_3689813
# 10 AGG AAG AGG > AAG      SNP gc19_pc.promCore  55101488 M_12128328
# 13 TCG TAG TCG > TAG      SNP     gc19_pc.3utr 138644835 M_13570344
# 39 GGG GAG GGG > GAG      SNP     gc19_pc.3utr  99930014 M_36452836
# 61 GCT GTT GCT > GTT      SNP     gc19_pc.5utr  71351548 M_44844897
# 62 GCT GTT GCT > GTT      SNP gc19_pc.promCore  71351548 M_44844897
# elemenntLength     score pElement nPred        Mu        sd
# 3            9318 -0.360796       NA    NA 0.2760951 0.5353005
# 10            402 -0.155770       NA    NA 0.5109484 0.6123444
# 13           1334 -0.134984       NA    NA 0.2760951 0.5353005
# 39           4624  0.199252       NA    NA 0.2760951 0.5353005
# 61            116  1.008628       NA    NA 0.5972233 0.6331138
# 62            173  1.008628       NA    NA 0.5109484 0.6123444


dig_input$X1 <- as.numeric(dig_input$X1)
dig_input$X2 <- as.numeric(dig_input$X2)
setorder(dig_input, X1, X2)

for (i in 2:3) {
  dig_input[[i]] <- sprintf("%.0f", as.numeric(dig_input[[i]]))
}

fwrite(
  dig_input,
  file = "../../Dig/inputs/mutation_files/PCAWG/ICGC_and_TCGA/my_pancanStar.tsv",
  sep = "\t",
  col.names = FALSE
)

# fwrite(
#   dig_input,
#   file = "../../Dig/inputs/mutation_files/PCAWG/ICGC_and_TCGA/rmHypers/my_pancanStar.tsv",
#   sep = "\t",
#   col.names = FALSE
# )

# fwrite(
#   dig_input,
#   file = "../../Dig/inputs/mutation_files/PCAWG/ICGC_and_TCGA/rmHypers/Pan_Cancer.tsv",
#   sep = "\t",
#   col.names = FALSE
# )


# fwrite(
#   dig_input,
#   file = "../../Dig/inputs/mutation_files/PCAWG/ICGC_and_TCGA/myPan_Cancer_SNV_MNV_INDEL.annot.txt",
#   sep = "\t",
#   col.names = FALSE,
#   quote = FALSE
# )


table(dig_gr$V8 )

dig_input[dig_input$X1]

