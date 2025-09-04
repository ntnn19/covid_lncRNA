library(biomaRt)
library(tidyverse)

# Workaround for dplyr/biomaRt compatibility issue
# Temporarily mask the problematic dplyr function
if (packageVersion("dplyr") >= "1.1.0") {
  # Create a temporary environment to avoid conflicts
  old_collect <- dplyr::collect
  assignInNamespace("collect", function(.data, ..., n = Inf) {
    old_collect(.data, n = n)
  }, ns = "dplyr")
}

# Add debugging information
cat("Starting gene2symbol script...\n")
cat("Species:", snakemake@params[["species"]], "\n")
cat("Input file:", snakemake@input[["counts"]], "\n")
cat("Output file:", snakemake@output[["symbol"]], "\n")

# this variable holds a mirror name until
# useEnsembl succeeds ("www" is last, because
# of very frequent "Internal Server Error"s)
mart <- "useast"
rounds <- 0
while ( class(mart)[[1]] != "Mart" ) {
  cat("Attempting to connect to biomaRt, mirror:", mart, "round:", rounds + 1, "\n")
  mart <- tryCatch(
    {
      biomartCacheClear()
      # done here, because error function does not
      # modify outer scope variables, I tried
      if (mart == "www") rounds <- rounds + 1
      # equivalent to useMart, but you can choose
      # the mirror instead of specifying a host
      biomaRt::useEnsembl(
        biomart = "ENSEMBL_MART_ENSEMBL",
        dataset = str_c(snakemake@params[["species"]], "_gene_ensembl"),
        mirror = mart
      )
    },
    error = function(e) {
      cat("Error connecting to mirror", mart, ":", conditionMessage(e), "\n")
      # change or make configurable if you want more or
      # less rounds of tries of all the mirrors
      if (rounds >= 3) {
        stop("Failed to connect to biomaRt after 3 rounds of trying all mirrors. Last error: ", conditionMessage(e))
      }
      # hop to next mirror
      mart <- switch(mart,
                     useast = "uswest",
                     uswest = "asia",
                     asia = "www",
                     www = {
                       # wait before starting another round through the mirrors,
                       # hoping that intermittent problems disappear
                       cat("Waiting 30 seconds before next round...\n")
                       Sys.sleep(30)
                       "useast"
                     }
              )
    }
  )
}

cat("Successfully connected to biomaRt using mirror:", mart@host, "\n")

# Read input file with error handling
tryCatch({
  df <- read.table(snakemake@input[["counts"]], sep='\t', header=TRUE)
  cat("Successfully read input file with", nrow(df), "rows\n")
}, error = function(e) {
  stop("Failed to read input file: ", conditionMessage(e))
})

# Check if df has the expected 'gene' column
if (!"gene" %in% colnames(df)) {
  stop("Input file does not contain a 'gene' column. Available columns: ", paste(colnames(df), collapse=", "))
}

cat("Querying biomaRt for", length(unique(df$gene)), "unique genes...\n")

# Query biomaRt with error handling
tryCatch({
  g2g <- biomaRt::getBM(
              attributes = c( "ensembl_gene_id",
                              "external_gene_name"),
              filters = "ensembl_gene_id",
              values = df$gene,
              mart = mart
              )
  cat("Retrieved annotations for", nrow(g2g), "genes\n")
}, error = function(e) {
  stop("Failed to query biomaRt: ", conditionMessage(e))
})

# Merge and process data
tryCatch({
  annotated <- merge(df, g2g, by.x="gene", by.y="ensembl_gene_id", all.x=TRUE)

  # Handle genes that didn't get annotations
  annotated$gene <- ifelse(is.na(annotated$external_gene_name) | annotated$external_gene_name == '',
                          annotated$gene,
                          annotated$external_gene_name)
  annotated$external_gene_name <- NULL

  cat("Final dataset has", nrow(annotated), "rows\n")
  cat("Genes with symbols:", sum(!is.na(annotated$gene) & annotated$gene != ""), "\n")
}, error = function(e) {
  stop("Failed to merge data: ", conditionMessage(e))
})

# Write output with error handling
tryCatch({
  write.table(annotated, snakemake@output[["symbol"]], sep='\t', row.names=FALSE, quote=FALSE)
  cat("Successfully wrote output file\n")
}, error = function(e) {
  stop("Failed to write output file: ", conditionMessage(e))
})
