prob_strata_check <- function(data, treatment, probabilities, post_strat, treatments_for_compare) {
  if (is.null(probabilities) & is.null(post_strat)) stop("Either assignment probabilities or strata MUST be provided.")

  if (is.null(probabilities) & !is.null(post_strat)) {
    warning("Probability matrix is not provided. The method assumes treatment assignment probabilities are constant within each level of stratification.")
    return(list(status = "missing probability",
                prob_mat = prob_mat_generate(data, treatment, post_strat),
                post_strat = post_strat))
  }
  assert_subset(probabilities, names(data))
  assert_subset(treatment, names(data))
  assert_subset(as.character(treatments_for_compare), as.character(data[[treatment]]))

  if (!is.null(post_strat)) assert_subset(post_strat, names(data))

  if (!is.null(probabilities)) assert_subset(as.character(treatments_for_compare), probabilities)

  if (!is.null(probabilities) & !is.null(post_strat)) {
    prob_mat <- data[probabilities]
    valid_rows <- apply(prob_mat[treatments_for_compare] > 0, 1, all)
    data <- data[valid_rows, ]
    prob_mat <- prob_mat[valid_rows, ]

    unique_strata <- unique(data[[post_strat]])

    for (stratum in unique_strata) {
      subset_df <- prob_mat[data[[post_strat]] == stratum, ]

      unique_pairs <- unique(paste(subset_df[[treatments_for_compare[1]]], subset_df[[treatments_for_compare[2]]]))

      if (length(unique_pairs) > 1) {
        warning("Assignment probabilities vary within stratification levels; ignoring provided stratification and using unique probability levels for stratification instead.")
        return(list(status = "varying probability",
                    prob_mat = data[probabilities],
                    post_strat = NULL))
      }
    }
  }
  return(list(status = "normal",
              prob_mat = data[probabilities],
              post_strat = post_strat))
}

prob_mat_generate <- function(data, treatment, post_strat) {
  unique_trt <- unique(data[[treatment]])
  unique_strata <- unique(data[[post_strat]])
  prob_mat <- matrix(0, nrow = nrow(data), ncol = length(unique_trt))
  colnames(prob_mat) <- unique_trt

  for (stratum in unique_strata) {
    subset_rows <- data[[post_strat]] == stratum
    subset_df <- data[subset_rows, ]

    freq_table <- table(subset_df[[treatment]])

    full_freq <- setNames(rep(0, length(unique_trt)), unique_trt)
    full_freq[names(freq_table)] <- freq_table

    prob <- full_freq / sum(full_freq)

    prob_mat[subset_rows, ] <- matrix(rep(prob, sum(subset_rows)), ncol = length(unique_trt), byrow = TRUE)
  }
  return(as.data.frame(prob_mat))
}

