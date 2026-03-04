complete_df <- RobinCID::example
tx_colname <- "treatment"
treatment_levels <- unique(complete_df[[tx_colname]])
randomization_var_colnames <- c("t", "subtype")
data_sim <- complete_df[c("xb", "xc", "s12", "s13", tx_colname, randomization_var_colnames, "y", "y_b")]
randomization_table <- unique(complete_df[c(randomization_var_colnames, treatment_levels)])

test_that("Contrast utility functions compute correctly", {
  # h_diff
  x2 <- c(0.3, 0.5)
  expect_equal(h_diff(x2), 0.2)

  x3 <- c(1, 3, 5)
  diffs <- h_diff(x3)
  expect_length(diffs, 3)
  expect_equal(diffs, c(2, 4, 2))

  # h_jac_diff
  jac <- h_jac_diff(x2)
  expect_equal(dim(jac), c(1, 2))
  expect_equal(jac[1, ], c(-1, 1))

  # h_ratio
  expect_equal(h_ratio(c(2, 4)), 2)

  # h_jac_ratio matches numerical Jacobian
  x <- c(2, 4)
  analytical <- h_jac_ratio(x)
  numerical <- numDeriv::jacobian(h_ratio, x)
  expect_equal(analytical, numerical, tolerance = 1e-6)

  # h_odds_ratio
  x <- c(0.3, 0.6)
  or <- h_odds_ratio(x)
  expected <- (0.6 / 0.4) / (0.3 / 0.7)
  expect_equal(or, expected, tolerance = 1e-10)

  # h_jac_odds_ratio matches numerical Jacobian
  analytical <- h_jac_odds_ratio(x)
  numerical <- numDeriv::jacobian(h_odds_ratio, x)
  expect_equal(analytical, numerical, tolerance = 1e-6)

  # h_jac_odds_ratio rejects values outside [0, 1]
  expect_error(h_jac_odds_ratio(c(0.3, 1.5)))

  # h_lower_tri_idx
  idx <- h_lower_tri_idx(3)
  expect_equal(nrow(idx), 3)
  expect_equal(idx[, 1], c(2, 3, 3))
  expect_equal(idx[, 2], c(1, 1, 2))
})

test_that("Custom contrast through robin_wt works", {
  # Custom function + Jacobian should produce same result as built-in "difference"
  result_builtin <- robin_wt(
    data = data_sim,
    estimand = list(tx_colname = tx_colname,
                    tx_to_compare = c("trt.1", "trt.2")),
    design = list(randomization_var_colnames = randomization_var_colnames,
                  randomization_table = randomization_table),
    estimated_propensity = FALSE,
    outcome_model = list(formula = y ~ xb + xc, family = gaussian()),
    contrast_specs = list(contrast = "difference")
  )

  result_custom <- robin_wt(
    data = data_sim,
    estimand = list(tx_colname = tx_colname,
                    tx_to_compare = c("trt.1", "trt.2")),
    design = list(randomization_var_colnames = randomization_var_colnames,
                  randomization_table = randomization_table),
    estimated_propensity = FALSE,
    outcome_model = list(formula = y ~ xb + xc, family = gaussian()),
    contrast_specs = list(contrast = h_diff, contrast_jac = h_jac_diff)
  )

  expect_equal(result_custom$trt_effect, result_builtin$trt_effect, tolerance = 1e-10)
  expect_equal(result_custom$variance, result_builtin$variance, tolerance = 1e-10)

  # Custom function without Jacobian (NULL) — exercises numDeriv path
  result_no_jac <- robin_wt(
    data = data_sim,
    estimand = list(tx_colname = tx_colname,
                    tx_to_compare = c("trt.1", "trt.2")),
    design = list(randomization_var_colnames = randomization_var_colnames,
                  randomization_table = randomization_table),
    estimated_propensity = FALSE,
    outcome_model = list(formula = y ~ xb + xc, family = gaussian()),
    contrast_specs = list(contrast = h_diff, contrast_jac = NULL)
  )

  expect_equal(result_no_jac$trt_effect, result_builtin$trt_effect, tolerance = 1e-10)
  expect_equal(result_no_jac$variance, result_builtin$variance, tolerance = 1e-4)
})

test_that("Negative binomial family works", {
  set.seed(123)
  data_nb <- data_sim
  data_nb$y_count <- rpois(nrow(data_nb), exp(data_nb$xb / 2))

  expect_no_error(
    suppressWarnings(
      robin_wt(
        data = data_nb,
        estimand = list(tx_colname = tx_colname,
                        tx_to_compare = c("trt.1", "trt.2")),
        design = list(
          randomization_var_colnames = randomization_var_colnames,
          randomization_table = randomization_table),
        estimated_propensity = FALSE,
        outcome_model = list(formula = y_count ~ xb + xc,
                             family = MASS::negative.binomial(NA))
      )
    )
  )
})

test_that("Input validation catches bad inputs", {
  base_args <- list(
    data = data_sim,
    design = list(randomization_var_colnames = randomization_var_colnames,
                  randomization_table = randomization_table),
    estimated_propensity = FALSE,
    outcome_model = list(formula = y ~ xb + xc, family = gaussian())
  )

  # tx_to_compare with 1 element
  expect_error(
    do.call(robin_wt, c(base_args, list(
      estimand = list(tx_colname = tx_colname, tx_to_compare = c("trt.1"))
    ))),
    "exactly two"
  )

  # tx_to_compare with 3 elements
  expect_error(
    do.call(robin_wt, c(base_args, list(
      estimand = list(tx_colname = tx_colname, tx_to_compare = c("trt.1", "trt.2", "trt.3"))
    ))),
    "exactly two"
  )

  # alpha out of range
  expect_error(
    do.call(robin_wt, c(base_args, list(
      estimand = list(tx_colname = tx_colname, tx_to_compare = c("trt.1", "trt.2")),
      alpha = -0.1
    )))
  )

  expect_error(
    do.call(robin_wt, c(base_args, list(
      estimand = list(tx_colname = tx_colname, tx_to_compare = c("trt.1", "trt.2")),
      alpha = 1.5
    )))
  )

  # Invalid contrast string
  expect_error(
    do.call(robin_wt, c(base_args, list(
      estimand = list(tx_colname = tx_colname, tx_to_compare = c("trt.1", "trt.2")),
      contrast_specs = list(contrast = "ratio")
    )))
  )
})
