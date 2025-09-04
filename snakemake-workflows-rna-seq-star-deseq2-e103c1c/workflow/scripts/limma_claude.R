log <- file(snakemake@log[[1]], open = "wt")
sink(log)
sink(log, type="message")

library(stringr)
library(limma)

parallel <- FALSE
if (snakemake@threads > 1) {
    library("BiocParallel")
    # setup parallelization
    register(MulticoreParam(snakemake@threads))
    parallel <- TRUE
}

counts_data <- read.table(
  snakemake@input[["counts"]],
  header = TRUE,
  row.names = "gene",
  check.names = FALSE
)
counts_data <- counts_data[, order(names(counts_data))]

col_data <- read.table(
  snakemake@config[["samples"]],
  header = TRUE,
  row.names = "sample_name",
  check.names = FALSE
)
col_data <- col_data[order(row.names(col_data)), , drop = FALSE]

# Collapse the columns into a single condition label
col_data$condition <- with(col_data, ifelse(IFN_alpha == "_treated", "IFN_alpha_treated",
                             ifelse(IFN_lambda == "_treated", "IFN_lambda_treated",
                             ifelse(SARS_CoV_2_16h == "_treated", "SARS_CoV_2_16h_treated",
                             ifelse(SARS_CoV_2_72h == "_treated", "SARS_CoV_2_72h_treated",
                                    "Untreated")))))

# Turn into factor with fixed baseline ordering
col_data$condition <- factor(col_data$condition)
colnames(col_data) <- make.names(colnames(col_data))

# Build up formula with additive batch_effects and all interactions between the
# variables_of_interest
design_formula <- snakemake@config[["diffexp"]][["model"]]

# Optional: filter low-expression genes (similar to rowSums > 1 in DESeq2)
counts_data <- counts_data[rowMeans(counts_data) > 0, ]

# Create model matrix
design <- model.matrix(as.formula(design_formula), data = col_data)
colnames(design) <- sub("^condition", "", colnames(design))

# Fit linear model
fit <- lmFit(counts_data, design)

# Get the specific contrast configuration
contrast_name <- snakemake@wildcards[["contrast"]]
contrast_config <- snakemake@config[["diffexp"]][["contrasts"]][[contrast_name]]

# Handle different contrast specification formats
if (typeof(contrast_config) == "list" && "variable_of_interest" %in% names(contrast_config)) {
  # Structured format: variable_of_interest + level_of_interest + base_level
  variable <- contrast_config[["variable_of_interest"]]
  level_of_interest <- contrast_config[["level_of_interest"]]

  # Get base_level from contrast config, variables_of_interest, or factor levels
  if ("base_level" %in% names(contrast_config)) {
    base_level <- contrast_config[["base_level"]]
  } else {
    # Use the reference level of the factor (first level)
    base_level <- levels(col_data[[variable]])[1]
  }

  contrast_expression <- paste0(level_of_interest, " - ", base_level)

} else if (typeof(contrast_config) == "list" && "contrast_expression" %in% names(contrast_config)) {
  # Direct contrast expression format
  contrast_expression <- contrast_config[["contrast_expression"]]

} else if (length(contrast_config) == 1 && typeof(contrast_config) == "character") {
  # Simple string format: "level1 - level2"
  contrast_expression <- contrast_config

} else {
  stop("Unsupported contrast configuration format for: ", contrast_name,
       ". Expected either structured format with 'variable_of_interest' or 'contrast_expression'.")
}

# Create contrast matrix for this single contrast
contrast_args <- list()
contrast_args[[contrast_name]] <- contrast_expression
contrast_args[["levels"]] <- design
contrast_matrix <- do.call(makeContrasts, contrast_args)

fit_contrasts <- contrasts.fit(fit, contrast_matrix)
# Apply empirical Bayes smoothing
final_fit <- eBayes(fit_contrasts)

# Extract results for this specific contrast
result <- topTable(final_fit, coef = contrast_name, adjust = "BH", number = Inf)
write.csv(result, file = snakemake@output[["table"]], row.names = TRUE)

cat("Written results for contrast:", contrast_name, "to", snakemake@output[["table"]], "\n")
