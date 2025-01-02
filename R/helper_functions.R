prob_strata_check <- function(data, treatment, prob_mat, stratification, treatments_for_compare) {
  status <- "normal"
  if (is.null(prob_mat) & is.null(stratification)) stop("Either assignment probabilities or strata MUST be provided.")
  if (is.null(prob_mat) & !is.null(stratification)) {
    status <- "missing prob_mat"
    warning("Probability matrix is not provided. The method assumes treatment assignment probabilities are constant within each level of stratification.")
  }

  assert_subset(treatment, names(data))
  assert_subset(as.character(treatments_for_compare), as.character(data[[treatment]]))
  assert_subset(stratification, names(data))
  if (!is.null(prob_mat)) assert_subset(as.character(treatments_for_compare), names(prob_mat))

  if (!is.null(prob_mat) & !is.null(stratification)) {
    data <- data[apply(prob_mat[treatments_for_compare] > 0, 1, all), ]
    prob_mat <- prob_mat[apply(prob_mat[treatments_for_compare] > 0, 1, all), ]

    unique_strata <- unique(data$stratification)

    for (stratum in unique_strata) {
      subset_df <- data[data$stratification == stratum, ]

      unique_pairs <- unique(paste(subset_df$pij, subset_df$pik))

      if (length(unique_pairs) > 1) {
        warning("Assignment probabilities vary within stratification levels;
                ignoring provided stratification and using unique probability levels for stratification instead.")
        status <- "varying prob_mat"
        return(status)
      }
    }
    return(status)
  }
  return(status)
}

prob_mat_generate <- function(data, treatment, stratification) {
  unique_trt <- unique(data[[treatment]])
  unique_strata <- unique(data[[stratification]])
  prob_matrix <- matrix(0, nrow = nrow(data), ncol = length(unique_trt))
  colnames(prob_matrix) <- unique_trt

  for (stratum in unique_strata) {
    # Subset data for the current stratum
    subset_df <- data[data[[stratification]] == stratum, ]

    # Calculate the frequency of each treatment level
    freq_table <- table(subset_df[[treatment]])

    # Ensure all treatments are represented in the frequency table
    full_freq <- setNames(rep(0, length(unique_trt)), unique_trt)
    full_freq[names(freq_table)] <- freq_table

    # Calculate probabilities
    prob <- full_freq / sum(full_freq)

    # Store probabilities in the list
    prob_matrix[data[[stratification]] == stratum, ] <- matrix(rep(prob, sum(data[[stratification]] == stratum)), ncol = length(unique_trt), byrow = TRUE)
  }
  return(data.frame(prob_matrix))
}
