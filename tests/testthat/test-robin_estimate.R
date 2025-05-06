complete_df<- RobinCID::example
tx_colname <- "treatment"
treatment_levels <- unique(complete_df[[tx_colname]])
randomization_var_colnames <- c("t", "subtype")
data_sim <- complete_df[c("xb", "xc", "s12", "s13", tx_colname, randomization_var_colnames, "y", "y_b")]
randomization_table <- unique(complete_df[c(randomization_var_colnames, treatment_levels)])

test_that("robin_wt works correctly", {

  expect_silent(
    robin_wt(
      data = data_sim,
      estimand = list(tx_colname = tx_colname,
                      comparison = c('trt.1', 'trt.2')),
      design = list(randomization_var_colnames = randomization_var_colnames,
                    randomization_table = randomization_table),
      estimated_propensity = FALSE,
      outcome_model = list(formula = y ~ xb + xc,
                           family = gaussian())
    )
  )
  expect_silent(
    robin_wt(
      data = data_sim,
      estimand = list(tx_colname = tx_colname,
                      comparison = c('trt.1', 'trt.3')),
      design = list(randomization_var_colnames = randomization_var_colnames,
                    randomization_table = randomization_table),
      estimated_propensity = FALSE,
      outcome_model = list(formula = y ~ xb + xc,
                           family = gaussian())
    )
  )
  expect_silent(
    robin_wt(
      data = data_sim,
      estimand = list(tx_colname = tx_colname,
                      comparison = c('trt.1', 'trt.2')),
      design = list(randomization_var_colnames = randomization_var_colnames,
                    randomization_table = randomization_table),
      estimated_propensity = FALSE,
      outcome_model = list(formula = y ~ xb + xc,
                           family = gaussian()),
      contrast_specs = list(contrast = "risk_ratio")
    )
  )
  expect_silent(
    robin_wt(
      data = data_sim,
      estimand = list(tx_colname = tx_colname,
                      comparison = c('trt.1', 'trt.2')),
      design = list(randomization_var_colnames = randomization_var_colnames,
                    randomization_table = randomization_table),
      estimated_propensity = FALSE,
      outcome_model = list(formula = y_b ~ xb + xc,
                           family = binomial()),
      contrast_specs = list(contrast = "odds_ratio")
    )
  )

  expect_silent(
    robin_wt(
      data = data_sim,
      estimand = list(tx_colname = tx_colname,
                      comparison = c('trt.1', 'trt.2')),
      design = list(randomization_var_colnames = randomization_var_colnames,
                    randomization_table = randomization_table),
      estimated_propensity = TRUE,
      outcome_model = list(formula = y_b ~ xb + xc,
                           family = binomial()),
      contrast_specs = list(contrast = "odds_ratio")
    )
  )

})

test_that("robin_wt works correctly under binomial()", {
  expect_silent(
    robin_wt(
      data = data_sim,
      estimand = list(tx_colname = tx_colname,
                      comparison = c('trt.1', 'trt.2')),
      design = list(randomization_var_colnames = randomization_var_colnames,
                    randomization_table = randomization_table),
      estimated_propensity = FALSE,
      outcome_model = list(formula = y_b ~ xb + xc,
                           family = binomial())
    )
  )
})

test_that("robin_wt snapshot",{
  expect_snapshot(
    robin_wt(
      data = data_sim,
      estimand = list(tx_colname = tx_colname,
                      comparison = c('trt.1', 'trt.2')),
      design = list(randomization_var_colnames = randomization_var_colnames,
                    randomization_table = randomization_table),
      estimated_propensity = FALSE,
      outcome_model = list(formula = y ~ xb + xc,
                           family = gaussian())
    )
  )

  expect_snapshot(
    robin_wt(
      data = data_sim,
      estimand = list(tx_colname = tx_colname,
                      comparison = c('trt.1', 'trt.2')),
      design = list(randomization_var_colnames = randomization_var_colnames,
                    randomization_table = NULL),
      estimated_propensity = TRUE,
      outcome_model = list(formula = y ~ xb + xc,
                           family = gaussian())
    )
  )
})

test_that("robin_ps works correctly", {
  expect_silent(
    robin_ps(
      data = data_sim,
      estimand = list(tx_colname = tx_colname,
                      comparison = c('trt.1', "trt.2")),
      design = list(randomization_var_colnames = randomization_var_colnames,
                    randomization_table = randomization_table),
      stratify_by = NULL,
      outcome_model = list(formula = y ~ xb + xc,
                           family = gaussian())
    )
  )
  expect_silent(
    robin_ps(
      data = data_sim,
      estimand = list(tx_colname = tx_colname,
                      comparison = c('trt.1', "trt.3")),
      design = list(randomization_var_colnames = randomization_var_colnames,
                    randomization_table = randomization_table),
      stratify_by = NULL,
      outcome_model = list(formula = y ~ xb + xc,
                           family = gaussian())
    )
  )
  expect_warning(
    robin_ps(
      data = data_sim,
      estimand = list(tx_colname = tx_colname,
                      comparison = c('trt.1', "trt.3")),
      design = list(randomization_var_colnames = randomization_var_colnames,
                    randomization_table = NULL),
      stratify_by = "s13",
      outcome_model = list(formula = y ~ xb + xc,
                           family = gaussian())
    )
  )
  expect_silent(
    robin_ps(
      data = data_sim,
      estimand = list(tx_colname = tx_colname,
                      comparison = c('trt.1', "trt.2")),
      design = list(randomization_var_colnames = randomization_var_colnames,
                    randomization_table = randomization_table),
      stratify_by = NULL,
      outcome_model = list(formula = y_b ~ xb + xc,
                           family = binomial())
    )
  )
})



test_that("robin_ps snapshot",{
  expect_snapshot(
    robin_ps(
      data = data_sim,
      estimand = list(tx_colname = tx_colname,
                      comparison = c('trt.1', "trt.2")),
      design = list(randomization_var_colnames = randomization_var_colnames,
                    randomization_table = NULL),
      stratify_by = "s12",
      outcome_model = list(formula = y ~ xb + xc,
                           family = gaussian())
    )
  )
})

test_that("robin_ps snapshot2", {
  expect_snapshot(
    robin_ps(
      data = data_sim,
      estimand = list(tx_colname = tx_colname,
                      comparison = c('trt.1', "trt.2")),
      design = list(randomization_var_colnames = randomization_var_colnames,
                    randomization_table = randomization_table),
      stratify_by = NULL,
      outcome_model = list(formula = y ~ xb + xc,
                           family = gaussian())
    )
  )
})

