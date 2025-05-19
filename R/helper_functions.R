.return.error <- function(err){
  if(length(err) > 0) stop(paste0(err, sep="\n"))
}

.check.null <- function(estimand, design, outcome_model){
  fields <- list(
    "estimand$tx_colname" = estimand$tx_colname,
    "estimand$tx_to_compare" = estimand$tx_to_compare,
    "design$randomization_var_colnames" = design$randomization_var_colnames,
    "outcome_model$formula" = outcome_model$formula
  )

  null_fields <- names(fields)[vapply(fields, is.null, logical(1))]

  if (length(null_fields) > 0) {
    return(paste0("The following must be provided: ", paste(null_fields, collapse = ", ")))
  }
}

.check.wt <- function(design, estimated_propensity, method){
  if(is.null(design$randomization_table) & !estimated_propensity & method == "wt") {
    return("The randomization_table must be provided if estimated_propensity is FALSE")
  }
}

validate_inputs <- function(estimand, design, outcome_model, estimated_propensity, method){
  errors <- character()
  errors <- c(errors, .check.null(estimand, design, outcome_model))
  errors <- c(errors, .check.wt(design, estimated_propensity, method))
  .return.error(errors)
}

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
        proportion_values <- round(as.numeric(count_values) / nrow(data), 2)

        strata_df <- as.data.frame(round(strata_probs, 2))
        strata_df$Sample.Size <- as.numeric(count_values)
        strata_df$Proportion <- proportion_values
        strata_df$Stratum <- stratum
        return(strata_df)
    })
    prob_tab <- do.call(rbind, strata_list)
  } else {  # Case when stratify_by is not provided
    # Collapse each row of Z to a string key
    Z_keys <- apply(data[Z], 1, paste, collapse = ", ")
    prob_keys <- apply(prob_mat, 1, paste, collapse = ", ")
    group_keys <- paste(Z_keys, prob_keys, sep = " | ")

    key_table <- table(group_keys)
    proportions <- round(key_table / nrow(data), 2)

    # Reconstruct the original values for output
    keys_split <- strsplit(names(key_table), " \\| ")
    Z_values <- do.call(rbind, lapply(keys_split, function(x) strsplit(x[1], ", ")[[1]]))
    prob_values <- do.call(rbind, lapply(keys_split, function(x) {
      round(as.numeric(strsplit(x[2], ", ")[[1]]), 2)
    }))

    prob_tab <- data.frame(Z_values, prob_values)
    names(prob_tab) <- c(Z, trt_cols)
    prob_tab$Sample.Size <- as.numeric(key_table)
    prob_tab$Proportion <- as.numeric(proportions)
  }

  return(prob_tab)
}

consistency_check <- function(data,
                              estimand,
                              design = list(randomization_var_colnames = NULL,
                                            randomization_table = NULL),
                              stratify_by = NULL){
  treatment <- estimand$tx_colname
  treatments_for_compare <- estimand$tx_to_compare

  rand_table <- design$randomization_table
  rand_var <- design$randomization_var_colnames

  treatment_names <- unique(data[[treatment]])
  assert_subset(treatments_for_compare, treatment_names)
  assert(is.data.frame(data))
  assert(is.null(rand_table) || is.data.frame(rand_table))
  assert_subset(rand_var, names(data))

  if(!is.null(rand_table)) {
    if(!test_subset(treatment_names, names(rand_table))) {stop("The randomization table is incomplete\u2014some treatments are missing from randomization_table.\nPlease revise the table.")}
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
          "Multiple entries in the randomization_table for {%s} = {%s}.\n Please revise the table and provide unique randomization probabilities for each combination of {%s} appearing in the data.",
          paste(rand_var, collapse = ", "),
          paste(combination, collapse = ", "),
          paste(rand_var, collapse = ", ")
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
          stop("Treatment assignment probabilities for the arms being compared (as specified in randomization_table) must be constant within each level of stratify_by.\nPlease revise the randomization_table for accuracy or change the stratify_by variable.")
        }
      }
    }
  } else {
    if(!is.null(stratify_by)){
      warning("Assignment probabilities for the arms being compared must be constant within each level of stratify_by.\nThis check is skipped because randomization_table was not provided.")
      }
  }
}

#' @title Assign Probability according to Design
#' @description
#' Assign Probability according to Design
#'
#' @param data (`data.frame`) Input data frame.
#' @param estimand (`list`) A list specifying the estimand.
#' @param design (`list`) A list describing the randomization design. See `Details`.
#' @param method estimation method.
#' @param estimated_propensity Whether to use estimated propensity score.
#' @param stratify_by The column name of stratification variable in `data`.
#' @return A new `data` with columns of the treatment assignment probability.
#' @export
#' @details
#' `design` has two elements: `randomization_var_colnames` (`vector`) and `randomization_table` (`data.frame`)
#'
assign_prob_and_strata <- function(data,
                                   estimand,
                                   design = list(randomization_var_colnames = NULL,
                                                 randomization_table = NULL),
                                   method,
                                   estimated_propensity = T,
                                   stratify_by = NULL
                                   ) {
  treatment <- estimand$tx_colname
  treatments_for_compare <- estimand$tx_to_compare

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

