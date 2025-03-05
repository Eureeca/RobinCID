prob_strata_check <- function(data, treatment, probabilities, post_strata, treatments_for_compare) {
  if (is.null(probabilities) & is.null(post_strata)) stop("Either assignment probabilities or strata MUST be provided.")

  if (is.null(probabilities) & !is.null(post_strata)) {
    warning("Assignment probabilities are not provided. The method assumes treatment assignment probabilities are constant within each level of stratification.")
    return(list(status = "missing probability",
                prob_mat = prob_mat_generate(data, treatment, post_strata),
                post_strata = post_strata))
  }
  assert_subset(probabilities, names(data))
  assert_subset(treatment, names(data))
  assert_subset(as.character(treatments_for_compare), as.character(data[[treatment]]))

  if (!is.null(post_strata)) assert_subset(post_strata, names(data))

  if (!is.null(probabilities)) assert_subset(as.character(treatments_for_compare), probabilities)

  if (!is.null(probabilities) & !is.null(post_strata)) {
    prob_mat <- data[probabilities]
    valid_rows <- apply(prob_mat[treatments_for_compare] > 0, 1, all)
    data <- data[valid_rows, ]
    prob_mat <- prob_mat[valid_rows, ]

    unique_strata <- unique(data[[post_strata]])

    for (stratum in unique_strata) {
      subset_df <- prob_mat[data[[post_strata]] == stratum, ]

      unique_pairs <- unique(paste(subset_df[[treatments_for_compare[1]]], subset_df[[treatments_for_compare[2]]]))

      if (length(unique_pairs) > 1) {
        warning("Assignment probabilities vary within stratification levels; ignoring provided stratification and using unique probability levels for stratification instead.")
        return(list(status = "varying probability",
                    prob_mat = data[probabilities],
                    post_strata = NULL))
      }
    }
  }
  return(list(status = "normal",
              prob_mat = data[probabilities],
              post_strata = post_strata))
}

prob_mat_generate <- function(data, treatment, post_strata) {
  unique_trt <- unique(data[[treatment]])
  unique_strata <- unique(data[[post_strata]])
  prob_mat <- matrix(0, nrow = nrow(data), ncol = length(unique_trt))
  colnames(prob_mat) <- unique_trt

  for (stratum in unique_strata) {
    subset_rows <- data[[post_strata]] == stratum
    subset_df <- data[subset_rows, ]

    freq_table <- table(subset_df[[treatment]])

    full_freq <- setNames(rep(0, length(unique_trt)), unique_trt)
    full_freq[names(freq_table)] <- freq_table

    prob <- full_freq / sum(full_freq)

    prob_mat[subset_rows, ] <- matrix(rep(prob, sum(subset_rows)), ncol = length(unique_trt), byrow = TRUE)
  }
  res <- as.data.frame(prob_mat)
  return(res)
}

prob_table_generate <- function(data, post_strata, prob_mat, Z){


  if(is.null(Z)){
    row_counts <- apply(prob_mat, 1, function(row) paste(row, collapse = ", "))
    unique_counts <- table(row_counts)
    proportions <- round(unique_counts / nrow(prob_mat), 2)
#
    if(is.null(post_strata)){
      prob_tab <- matrix(
        c(names(proportions),
          unique_counts,
          proportions),
        nrow = length(names(proportions))
      )
      colnames(prob_tab) <- c("Unique.Level", "Sample.Size", "Proportion")
    } else {
      unique_strata <- unique(data[[post_strata]])

      strata_mapping <- lapply(unique_strata, function(stratum) {
        strata_indices <- which(data[[post_strata]] == stratum)
        prob_values <- apply(prob_mat[strata_indices, ], 1, function(row) paste(row, collapse = ", "))
        unique(prob_values)
      })
      unique_counts <- table(data[[post_strata]])[unique_strata]
      proportions <- round(unique_counts / nrow(prob_mat), 2)
      prob_tab <- matrix(
        c(unique_strata,
          unlist(strata_mapping),
          unique_counts,
          proportions),
        nrow = length(names(proportions))
      )
      colnames(prob_tab) <- c(paste("Stratum.", post_strata, sep=""), "Unique.Level", "Sample.Size", "Proportion")
    }
  } else{
    unique_Z <- as.data.frame(unique(data[Z]))
    Z_counts <- apply(data[Z], 1, function(row) paste(row, collapse = ", "))
    unique_counts <- table(Z_counts)

    prob_aligned <- apply(unique_Z, 1,
                          function(row){paste(prob_mat[apply(data[Z],1,function(row2) all(row2==row)), ][1,], collapse = ", ")})
    proportions <- round(unique_counts / nrow(data), 2)

    prob_tab <- matrix(
      c(unlist(unique_Z),
        prob_aligned,
        unique_counts,
        proportions),
      nrow = length(names(proportions))
    )
    colnames(prob_tab) <- c(paste0("Z.", Z), "Unique.Level", "Sample.Size", "Proportion")
  }

  prob_tab[order(prob_tab[ , "Unique.Level"]),]
}


#' @title Assign Probability and Strata according to Reference Table
#' @description
#' Assign Probability and Strata according to Reference Table
#'
#' @param data (`data.frame`) Input data frame.
#' @param Z (`vector`) The variable names in `data` and `Z_table` that define `stratum`.
#' @param Z_table (`data.frame`) A reference table where rows define strata and columns include `Z` and `probabilities`.
#' @param probabilities (`vector`) The names of columns in `Z_table` representing treatment assignment probabilities.
#' @param merge_strata (`logical`) Whether to merge strata with identical probabilities.
#' @export
#' @details
#' The output column name is `stratum`.
#'
#' @examples
#' Z <- c("t", "subtype")
#' probabilities <- c("trt.1", "trt.2", "trt.3", "trt.4")
#' Z_table <- unique(example[c(Z, probabilities)])
#' new_example = assign_prob_and_strata(
#'                             data = example[Z],
#'                             Z = Z,
#'                             Z_table =  Z_table,
#'                             probabilities = probabilities,
#'                             merge_strata = TRUE
#'                             )
#' head(new_example)
#'
assign_prob_and_strata <- function(data, Z, Z_table, probabilities, merge_strata = TRUE) {

  assert(is.data.frame(data), is.data.frame(Z_table))
  assert_subset(Z, names(data))
  assert_subset(Z, names(Z_table))
  assert_subset(probabilities, names(Z_table))

  Z_table <- unique(Z_table)

  if (merge_strata) {
    prob_signature <- apply(Z_table[probabilities], 1, paste, collapse = "_")
    Z_table$post_stratum <- as.numeric(as.factor(prob_signature))
  } else {
    Z_table$post_stratum <- seq_len(nrow(Z_table))
  }

  data <- merge(data, Z_table, by = Z, all.x = TRUE)

  if (any(is.na(data$post_stratum))) {
    warning("Some are not matched. Consider remove them or adjust Z and Z_table.")
  }

  return(data)
}

