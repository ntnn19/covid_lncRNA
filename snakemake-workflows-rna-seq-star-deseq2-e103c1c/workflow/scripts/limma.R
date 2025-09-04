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

# build up formula with additive batch_effects and all interactions between the
# variables_of_interes

design_formula <- snakemake@config[["diffexp"]][["model"]]

# Optional: filter low-expression genes (similar to rowSums > 1 in DESeq2)
counts_data <- counts_data[rowMeans(counts_data) > 0, ]

# Create model matrix
design <- model.matrix(as.formula(design_formula), data = col_data)
colnames(design) <- sub("^condition", "", colnames(design))

# Apply the contrast matrix to the design matrix with age adjustment

# Fit linear model
fit <- lmFit(counts_data, design)

function_args <- snakemake@config[["diffexp"]][["contrasts"]]
function_args[["levels"]] <- design
contrast_matrix <- do.call(makeContrasts,function_args)

fit_contrasts <- contrasts.fit(fit, contrast_matrix)
# Apply empirical Bayes smoothing
final_fit <- eBayes(fit_contrasts)

# Extract DE results
# Assuming you want top table for a specific contrast:
# Example: first coefficient
# Extract top differentially expressed proteins/genes for each contrast

write.csv(topTable(final_fit, coef=dimnames(final_fit$contrasts)$Contrasts, adjust="BH", number=Inf),   file = snakemake@output[["table"]], row.names=TRUE)

