complete_df <- RobinCID::example
estimand <- list(tx_colname = "treatment",
                 comparison = c("trt.1", "trt.2"))
randomization_var_colnames <- c("t", "subtype")
randomization_table <- unique(complete_df[c(randomization_var_colnames, "trt.1", "trt.2", "trt.3", "trt.4")])


outcome_model <- list(formula = y ~ xb + xc, family = gaussian())

df <- complete_df[c("y",
                    "xb",
                    "xc",
                    "treatment",
                    "t",
                    "subtype",
                    "s12",
                    "s12.2",
                    "s13")]

test_that("Logic is correct", {
  expect_error(
    robin_wt(
      data = df,
      estimand = estimand,
      design = list(
        randomization_var_colnames = c("t", "subtype", "s12"),
        randomization_table = randomization_table
      ),
      estimated_propensity = T,
      outcome_model = outcome_model
    )
  )

  expect_error(
    robin_wt(
      data = df,
      estimand = estimand,
      design = list(
        randomization_var_colnames = c("t", "subtype", "s100"),
        randomization_table = randomization_table
      ),
      estimated_propensity = T,
      outcome_model = outcome_model
    )
  )

  expect_no_message(
    robin_wt(
      data = df,
      estimand = estimand,
      design = list(
        randomization_var_colnames = c("t", "subtype"),
        randomization_table = randomization_table
      ),
      estimated_propensity = T,
      outcome_model = outcome_model
    )
  )

  expect_error(
    robin_wt(
      data = df,
      estimand = estimand,
      design = list(
        randomization_var_colnames = c("t", "subtype"),
        randomization_table = cbind(randomization_table, data.frame(z3 = rbinom(6, 1, 0.5)))
      ),
      estimated_propensity = F,
      outcome_model = outcome_model
    )
  )


  expect_error(
    robin_wt(
      data = df,
      estimand = estimand,
      design = list(
        randomization_var_colnames = c("t", "subtype"),
        randomization_table = NULL
      ),
      estimated_propensity = F,
      outcome_model = outcome_model
    )
  )

  expect_no_message(
    robin_wt(
      data = df,
      estimand = estimand,
      design = list(
        randomization_var_colnames = c("t", "subtype"),
        randomization_table = NULL
      ),
      estimated_propensity = T,
      outcome_model = outcome_model
    )
  )

  expect_no_message(
    robin_wt(
      data = df,
      estimand = estimand,
      design = list(
        randomization_var_colnames = c("t", "subtype"),
        randomization_table = randomization_table
      ),
      estimated_propensity = T,
      outcome_model = outcome_model
    )
  )

  expect_no_message(
    robin_wt(
      data = df,
      estimand = estimand,
      design = list(
        randomization_var_colnames = c("t", "subtype"),
        randomization_table = randomization_table
      ),
      estimated_propensity = F,
      outcome_model = outcome_model
    )
  )


  expect_no_message(
    robin_ps(
      data = df,
      estimand = estimand,
      design = list(
        randomization_var_colnames = c("t", "subtype"),
        randomization_table = randomization_table
      ),
      stratify_by = NULL,
      outcome_model = outcome_model
    )
  )

  expect_error(
    robin_ps(
      data = df,
      estimand = estimand,
      design = list(
        randomization_var_colnames = c("t", "subtype"),
        randomization_table = randomization_table[,1:3]
      ),
      stratify_by = NULL,
      outcome_model = outcome_model
    )
  )


  expect_error(
    robin_ps(
      data = df,
      estimand = estimand,
      design = list(
        randomization_var_colnames = c("t", "subtype"),
        randomization_table = randomization_table
      ),
      stratify_by = "s13",
      outcome_model = outcome_model
    )
  )

  wrong_randomization_table <- rbind(randomization_table, c(1, 1, 0.25, 0.25, 0.25, 0.25))
  expect_error(
    robin_ps(
      data = df,
      estimand = estimand,
      design = list(
        randomization_var_colnames = c("t", "subtype"),
        randomization_table = wrong_randomization_table
      ),
      stratify_by = "s13",
      outcome_model = outcome_model
    )
  )

  expect_error(
    robin_wt(
      data = df,
      estimand = estimand,
      design = list(
        randomization_var_colnames = c("t", "subtype"),
        randomization_table = wrong_randomization_table
      ),
      estimated_propensity =  T,
      outcome_model = outcome_model
    )
  )

})

test_that("Null check",{
  expect_error(
    robin_wt(
      data = df,
      estimand = estimand,
      design = list(
        randomization_var_colnames = NULL,
        randomization_table = NULL
      ),
      estimated_propensity = F,
      outcome_model = outcome_model
    )
  )
  expect_error(
    robin_wt(
      data = df,
      estimand = list(tx_colname = NULL,
                      comparison = NULL),
      design = list(
        randomization_var_colnames = NULL,
        randomization_table = randomization_table
      ),
      estimated_propensity = F,
      outcome_model = outcome_model
    )
  )
})
