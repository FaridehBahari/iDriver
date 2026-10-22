rm(list = ls())

library(data.table)
library(dplyr)
source('main/functions_iDriver.R')


all_scores <- fread('../extdata/procInput/iDriverInputs/PhredLikeScores_all.tsv')
mutData <- fread('../extdata/procInput/iDriverInputs/mutData.tsv')
fs_all <- left_join(mutData, all_scores)

path_donorInfo <- '../extdata/procInput/iDriverInputs/donorInfo.tsv'

cohort <- 'Pancan-no-skin-melanoma-lymph'
exclude_lymph_melanoma <- ifelse(cohort == 'Pancan-no-skin-melanoma-lymph', T, F)
donorInfo <- select_cohort(path_donorInfo,
                           cohort,
                           exclude_lymph_melanoma,
                           exclude_hyper_mutated = F)



fs <- fs_all[fs_all$donor_id %in% donorInfo$donor_id, c("PCAWG_ID" , "score_PhredLike")]

fs_mean <- fs[, .(
  mean_score = mean(score_PhredLike)
), by = PCAWG_ID]


# Merge with y (assuming rows are in the same order)
y <- fread('../extdata/output/observed/FI_free/rmHyperMuts/betaRhoPancan_lengthCat_quartilesElemType/observed_Pancan-no-skin-melanoma-lymph.tsv')
y$raw_q <- p.adjust(y$pvals, method = 'fdr')
res <- fread('../extdata/procInput/BMRs/observed/rmHyperMutated/Pancan-no-skin-melanoma-lymph/test_y.tsv')
res <- res[,c('binID', 'N', 'length')]
y <- left_join(y, res, by = c('PCAWG_ID' = 'binID', 'N'))
  
y <- left_join(y, fs_mean)

# Initialize average weights
avg_weight <- numeric(nrow(y))

# Elements with raw_q <= 0.25
near <- which(!is.na(y$raw_q) & y$raw_q <= 0.25)


cutoff <- 0.01




# Convert to Phred scale
threshold <- -10 * log10(cutoff)

# Calculate weight for near-significant elements

y$weight <- y$mean_score / threshold

# Set weight to 1 for the remaining elements
other <- which(!is.na(y$raw_q) & y$raw_q > 0.25)
y$weight[other] <- 1

# Replace missing weights with 1
y$weight[is.na(y$weight)] <- 1

# Calculate adjusted mutation burden
y$n_score <- y$weight * y$sum_nMut

# Calculate adjusted p-values
y$p_adj <- y$pvals

if (length(near) > 0) {
  
  count <- y$n_score[near]
  offset <- y$length[near] * y$N[near] + 1
  
  y$p_adj[near] <- 'again burden test'
}
  
# Calculate adjusted q-values
y$q_adj <- p.adjust(y$p_adj, method = 'fdr')
