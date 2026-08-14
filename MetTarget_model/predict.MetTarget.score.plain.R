#!/usr/bin/env Rscript

options(echo = TRUE)
Sys.time()

args <- commandArgs(trailingOnly = TRUE)
args

input_data_file <- args[1]
sheet_name <- args[2]
input_model_file <- args[3]

#" ----------
#"
#" # Plain MetTarget prediction helpers (base R only)

is_mettarget_plain_model <- function(model) {
  is.list(model) && identical(model$type, "mettarget_plain_model")
}

rbf_kernel <- function(x, support_vectors, sigma) {
  x <- as.matrix(x)
  sv <- as.matrix(support_vectors)
  x_norm <- rowSums(x^2)
  sv_norm <- rowSums(sv^2)
  dist2 <- outer(x_norm, sv_norm, "+") - 2 * tcrossprod(x, sv)
  dist2[dist2 < 0] <- 0
  exp(-sigma * dist2)
}

apply_preprocess <- function(newdata, center, scale) {
  newdata <- as.data.frame(newdata)
  if (!is.null(center) && length(center) > 0) {
    center <- center[colnames(newdata)]
    newdata <- sweep(newdata, 2, center, "-")
  }
  if (!is.null(scale) && length(scale) > 0) {
    scale <- scale[colnames(newdata)]
    newdata <- sweep(newdata, 2, scale, "/")
  }
  newdata
}

platt_probabilities <- function(decision, prob_model, class_levels) {
  if (nrow(prob_model) == 1 && length(class_levels) == 2) {
    pos_index <- match("pos", class_levels)
    if (is.na(pos_index)) {
      pos_index <- 2
    }
    probs <- matrix(NA_real_, nrow = length(decision), ncol = length(class_levels))
    colnames(probs) <- class_levels
    probs[, pos_index] <- 1 / (1 + exp(prob_model[1, 1] * decision + prob_model[1, 2]))
    neg_index <- setdiff(seq_along(class_levels), pos_index)
    probs[, neg_index] <- 1 - probs[, pos_index]
    return(probs)
  }

  probs <- matrix(
    1 / (1 + exp(prob_model[, 1] * decision + prob_model[, 2])),
    ncol = nrow(prob_model),
    nrow = length(decision),
    byrow = TRUE
  )
  colnames(probs) <- class_levels
  probs <- probs / rowSums(probs)
  probs
}

predict_mettarget_plain <- function(model, newdata) {
  if (!is_mettarget_plain_model(model)) {
    stop(
      "Expected a plain MetTarget model. ",
      "Use plain.MetTarget.pretrained.UC.RDS or plain.MetTarget.pretrained.CD.RDS."
    )
  }

  missing_features <- setdiff(model$feature_names, colnames(newdata))
  if (length(missing_features) > 0) {
    stop(
      "New data is missing model features: ",
      paste(missing_features, collapse = ", ")
    )
  }

  newdata <- newdata[, model$feature_names, drop = FALSE]
  newdata <- apply_preprocess(newdata, model$center, model$scale)
  kernel_matrix <- rbf_kernel(newdata, model$support_vectors, model$sigma)
  decision <- as.vector(kernel_matrix %*% model$coef + model$bias)
  probs <- platt_probabilities(decision, model$prob_model, model$class_levels)
  probs[, model$positive_class]
}

#" ----------
#"
#" # Load packages

library(openxlsx)

#" ----------
#"
#" # Input data

data_df <- as.data.frame(read.xlsx(input_data_file, sheet = sheet_name, colNames = TRUE))

if (any(grepl("GWAS", data_df[, 1]))) {
  data_df <- read.xlsx(input_data_file, sheet = sheet_name, colNames = FALSE)
  data_df <- as.data.frame(t(data_df), stringsAsFactors = FALSE)
  colnames(data_df) <- data_df[1, ]
  data_df <- data_df[-1, , drop = FALSE]
}

rownames(data_df) <- data_df[, "gene"]
data_df <- data_df[, -1, drop = FALSE]

for (col in 1:ncol(data_df)) data_df[, col] <- as.numeric(data_df[, col])
dim(data_df)

model_list <- readRDS(input_model_file)
length(model_list)

score_df <- data.frame(matrix(nrow = nrow(data_df), ncol = 0))
for (i in seq_along(model_list)) {
  new_df <- data_df
  new_df[is.na(new_df)] <- 0
  new_df[, "score"] <- predict_mettarget_plain(model_list[[i]], new_df)
  score_df <- cbind.data.frame(score_df, data.frame(new_df[, "score"]))
  colnames(score_df)[ncol(score_df)] <- paste0("model", i)
}

rownames(score_df) <- rownames(data_df)
score_df[, "predicted.MetTarget.score"] <- apply(score_df[, seq_along(model_list)], 1, function(x) {
  median(x, na.rm = TRUE)
})

output_df <- cbind.data.frame(
  "gene" = rownames(score_df),
  "predicted.MetTarget.score" = score_df[, "predicted.MetTarget.score"]
)

output_file <- paste0(gsub(" ", "", sheet_name), ".MetTarget.output.txt")
write.table(
  output_df,
  output_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = TRUE,
  na = ""
)

utils:::print.sessionInfo(sessionInfo()[-8])
sessionInfo()

#" ----------
#"
#" # Done.
