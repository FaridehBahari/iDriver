################################# (new) save final pvals Dig #############################################
rm(list = ls())

library(data.table)

input_dir  <- "../../Dig/output/elemDriver/withHyper"
output_dir <- "../../Dig/output/elemDriver/withHyper_final_pvals"

# Create output directory if it does not exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Get all cancer-type result files
files <- list.files(
  input_dir,
  pattern = "\\.results\\.txt$",
  full.names = TRUE
)

# Process each cancer type
for (file in files) {
  
  # Cancer type from filename
  cancer_type <- sub("\\.results\\.txt$", "", basename(file))
  
  if (cancer_type == 'myPan_Cancer') {
    cancer_type = 'Pan_Cancer'
  } else if (cancer_type == 'myPancancerStar') {
    cancer_type = 'Pancan-no-skin-melanoma-lymph'
  }
  
  # Read data
  dig <- fread(file)
  
  # Keep required columns and rename
  dig <- dig[, .(
    binID = ELT,
    p_value = PVAL_MUT_BURDEN
  )]
  
  # Output file
  output_file <- file.path(
    output_dir,
    paste0(cancer_type, ".results.tsv")
  )
  
  # Write
  fwrite(
    dig,
    file = output_file,
    sep = "\t"
  )
  
  cat("Processed:", cancer_type, "\n")
}


