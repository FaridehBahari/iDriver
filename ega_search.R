
library(httr)
library(jsonlite)
library(data.table)

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0) y else x
}



# Put your accessions into a character vector called egad_ids
egad_ids <- scan(
  "../EGA.txt",
  what = character(),
  sep = ",",
  strip.white = TRUE,
  quiet = TRUE
)
egad_ids <- unique(egad_ids)

get_ega_metadata <- function(id) {
  url <- paste0(
    "https://metadata.ega-archive.org/datasets/", id
  )
  
  tryCatch({
    r <- GET(url, timeout(20))
    if (status_code(r) != 200) return(NULL)
    
    fromJSON(
      content(r, as = "text", encoding = "UTF-8"),
      simplifyVector = FALSE
    )
  }, error = function(e) NULL)
}

results <- rbindlist(lapply(egad_ids, function(id) {
  x <- get_ega_metadata(id)
  
  if (is.null(x)) {
    return(data.table(
      accession = id,
      description = NA_character_,
      metadata = NA_character_,
      status = "Not retrieved"
    ))
  }
  
  txt <- tolower(toJSON(x, auto_unbox = TRUE))
  
  data.table(
    accession = id,
    description = as.character(x$title %||% x$description %||% NA),
    metadata = txt,
    status = "Retrieved"
  )
}), fill = TRUE)

results[, WGS_candidate := grepl(
  "whole.?genome sequencing|whole genome|wgs",
  paste(description, metadata),
  ignore.case = TRUE
)]

x <- results[results$WGS_candidate,]

fwrite(results, "ega_dataset_metadata.csv")
fwrite(x, "../ega_wgs_candidates.csv")
