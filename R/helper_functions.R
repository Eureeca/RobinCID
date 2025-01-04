prob_strata_check <- function(data, treatment, prob_mat, stratification, treatments_for_compare) {
  if (is.null(prob_mat) & is.null(stratification)) stop("Either assignment probabilities or strata MUST be provided.")

  if (is.null(prob_mat) & !is.null(stratification)) {
    warning("Probability matrix is not provided. The method assumes treatment assignment probabilities are constant within each level of stratification.")
    return("missing prob_mat")
  }

  assert_subset(treatment, names(data))
  assert_subset(as.character(treatments_for_compare), as.character(data[[treatment]]))

  if (!is.null(stratification)) assert_subset(stratification, names(data))

  if (!is.null(prob_mat)) assert_subset(as.character(treatments_for_compare), names(prob_mat))

  if (!is.null(prob_mat) & !is.null(stratification)) {
    valid_rows <- apply(prob_mat[treatments_for_compare] > 0, 1, all)
    data <- data[valid_rows, ]
    prob_mat <- prob_mat[valid_rows, ]

    unique_strata <- unique(data[[stratification]])

    for (stratum in unique_strata) {
      subset_df <- prob_mat[data[[stratification]] == stratum, ]

      unique_pairs <- unique(paste(subset_df[[treatments_for_compare[1]]], subset_df[[treatments_for_compare[2]]]))

      if (length(unique_pairs) > 1) {
        warning("Assignment probabilities vary within stratification levels; ignoring provided stratification and using unique probability levels for stratification instead.")
        return("varying prob_mat")
      }
    }
  }
  return("normal")
}

prob_mat_generate <- function(data, treatment, stratification) {
  unique_trt <- unique(data[[treatment]])
  unique_strata <- unique(data[[stratification]])
  prob_matrix <- matrix(0, nrow = nrow(data), ncol = length(unique_trt))
  colnames(prob_matrix) <- unique_trt

  for (stratum in unique_strata) {
    subset_rows <- data[[stratification]] == stratum
    subset_df <- data[subset_rows, ]

    freq_table <- table(subset_df[[treatment]])

    full_freq <- setNames(rep(0, length(unique_trt)), unique_trt)
    full_freq[names(freq_table)] <- freq_table

    prob <- full_freq / sum(full_freq)

    prob_matrix[subset_rows, ] <- matrix(rep(prob, sum(subset_rows)), ncol = length(unique_trt), byrow = TRUE)
  }
  return(as.data.frame(prob_matrix))
}

