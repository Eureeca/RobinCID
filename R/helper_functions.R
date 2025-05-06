prob_table_generate <- function(data, post_strata, prob_mat, Z) {

  # Get treatment column names from prob_mat
  trt_cols <- colnames(prob_mat)

  if (!is.null(post_strata)) {
    unique_probs <- unique(prob_mat)
    row_counts <- apply(prob_mat, 1, function(row) paste(row, collapse = ", "))
    unique_counts <- table(row_counts)
    proportions <- round(unique_counts / nrow(prob_mat), 2)

    prob_tab <- as.data.frame(unique_probs)
    prob_tab$Sample.Size <- as.numeric(unique_counts)
    prob_tab$Proportion <- as.numeric(proportions)

      unique_strata <- unique(data[[post_strata]])
      strata_list <- lapply(unique_strata, function(stratum) {
        strata_indices <- which(data[[post_strata]] == stratum)
        strata_probs <- unique(prob_mat[strata_indices, , drop = FALSE])
        count_values <- table(apply(prob_mat[strata_indices, , drop = FALSE], 1, paste, collapse = ", "))
        proportion_values <- round(as.numeric(count_values) / length(strata_indices), 2)

        strata_df <- as.data.frame(strata_probs)
        strata_df$Sample.Size <- as.numeric(count_values)
        strata_df$Proportion <- proportion_values
        strata_df$Stratum <- stratum
        return(strata_df)
      })
      prob_tab <- do.call(rbind, strata_list)
  } else {  # Case when stratify_by is not provided
    unique_Z <- unique(data[Z])
    unique_probs <- unique(prob_mat)

    Z_counts <- apply(data[Z], 1, function(row) paste(row, collapse = ", "))
    unique_counts <- table(Z_counts)
    proportions <- round(unique_counts / nrow(data), 2)

    prob_tab <- as.data.frame(unique_Z)
    prob_tab <- cbind(prob_tab, unique_probs)  # Keep treatment probabilities separate
    prob_tab$Sample.Size <- as.numeric(unique_counts)
    prob_tab$Proportion <- as.numeric(proportions)
  }

  return(prob_tab)
}

consistency_check <- function(data,
                              estimand,
                              design = list(randomization_var_colnames,
                                            randomization_table = NULL),
                              stratify_by = NULL){
  treatment <- estimand$tx_colname
  treatments_for_compare <- estimand$comparison

  rand_table <- design$randomization_table
  rand_var <- design$randomization_var_colnames

  treatment_names <- unique(data[[treatment]])
  assert_subset(treatments_for_compare, treatment_names)
  assert(is.data.frame(data))
  assert(is.null(rand_table) || is.data.frame(rand_table))
  assert_subset(rand_var, names(data))

  if(!is.null(rand_table)) {
    if(!test_subset(treatment_names, names(rand_table))) {stop("Some treatments are not in the randomization_table.\nPlease revise the table.")}
    if(!test_subset(rand_var, names(rand_table))) {stop("Some randomization variables in randomization_var_colnames are not in the randomization_table.\nPlease revise the table or drop the table and set estimated_propensity = TRUE.")}
    if(length(setdiff(names(rand_table),c(rand_var,as.character(treatment_names))))) {
      stop("The randomization_table contains additional variables that were not specified in randomization_var_colnames.\nPlease revise the table or drop the table and set estimated_propensity = TRUE.")
    }
    rand_var_combinations <- unique(rand_table[rand_var])

    for (i in seq_len(nrow(rand_var_combinations))) {
      combination <- rand_var_combinations[i, , drop = FALSE]
      matching_rows <- apply(rand_table[rand_var], 1, function(row) all(row == combination))
      prob_subset <- rand_table[matching_rows, as.character(treatment_names), drop = FALSE]

      # Check if all rows have the same assignment probabilities
      if (nrow(unique(prob_subset)) > 1) {
        stop(sprintf(
          "Inconsistent assignment probabilities detected with randomization variables %s taking $s.\nPlease revise the randomization_table.",
          paste(rand_var, collapse = ", "),
          paste(combination, collapse = ", ")
        ))
      }
    }


    if(!is.null(stratify_by)){
      reduced_table <- rand_table[c(rand_var, as.character(treatments_for_compare))]
      unique_probs <- unique(reduced_table)

      # Loop through each post-stratum and check consistency
      post_strata_levels <- unique(data[[stratify_by]])
      for (level in post_strata_levels) {
        # Subset data for the current post-stratum
        strata_data <- data[data[[stratify_by]] == level, ]

        # Extract the relevant rows from the reduced table
        strata_rand_table <- merge(
          strata_data[rand_var],
          unique_probs,
          by = rand_var,
          all.x = TRUE
        )

        # Check if the probabilities for the treatments are the same across rows
        prob_values <- unique(strata_rand_table[as.character(treatments_for_compare)])
        if (nrow(prob_values) > 1) {
          stop("Inconsistent treatment assignment probabilities in some post-strata.\nPlease check the randomization_table or stratify_by.")
        }
      }
    }
  } else {
    if(!is.null(stratify_by)){warning("The consistency of post stratification is not checked as randomization_table is not provided.")}
  }
}

#' @title Assign Probability according to Design
#' @description
#' Assign Probability according to Design
#'
#' @param data (`data.frame`) Input data frame.
#' @param estimand (`list`) A list specifying the estimand
#' @param design (`list`) A list describing the randomization design. See `Details`.
#' @param method estimation method.
#' @param stratify_by The column name of stratification variable in `data`.
#' @export
#' @details
#' `design` has two elements: `randomization_var_colnames` (`vector`) and `randomization_table` (`data.frame`)
#'
assign_prob_and_strata <- function(data,
                                   estimand,
                                   design = list(randomization_var_colnames,
                                                 randomization_table = NULL),
                                   method,
                                   estimated_propensity = T,
                                   stratify_by = NULL
                                   ) {
  treatment <- estimand$tx_colname
  treatments_for_compare <- estimand$comparison

  if(is.null(stratify_by)){
    Z <- design$randomization_var_colnames
    Z_table <- design$randomization_table
  } else{
    Z <- stratify_by
    Z_table <- NULL
  }
  if(estimated_propensity & method == "wt") Z_table = NULL
  treatment_names <- unique(data[[treatment]])


  if(!is.null(Z_table)) {

    data <- merge(data, Z_table, by = Z, all.x = TRUE)
  } else {
    strata_counts <- aggregate(rep(1, nrow(data)), by = c(data[Z], data[treatment]), FUN = sum)
    colnames(strata_counts) <- c(Z, treatment, "count")

    # Compute probabilities within each stratum
    total_counts <- aggregate(strata_counts$count, by = strata_counts[Z], FUN = sum)
    colnames(total_counts) <- c(Z, "total_count")

    strata_counts <- merge(strata_counts, total_counts, by = Z)
    strata_counts$prob <- strata_counts$count / strata_counts$total_count

    # Reshape to wide format (one column per treatment)
    strata_probs <- reshape(strata_counts[, c(Z, treatment, "prob")],
                            timevar = treatment,
                            idvar = Z,
                            direction = "wide")

    colnames(strata_probs) <- gsub("^prob\\.", "", colnames(strata_probs))

    # Ensure correct order of treatment probability columns
    prob_cols <- setdiff(names(strata_probs), Z)
    ordered_cols <- c(Z, intersect(treatment_names, prob_cols))
    strata_probs <- strata_probs[, ordered_cols, drop = FALSE]

    # Fill missing values with 0
    strata_probs[is.na(strata_probs)] <- 0

    # Merge estimated probabilities back into the original data
    data <- merge(data, strata_probs, by = Z, all.x = TRUE)
  }


  return(data)
}

