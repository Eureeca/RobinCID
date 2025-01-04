test_that("robin_wt works correctly", {
  formula = y ~ xb + xc
  expect_silent(
    robin_wt(
      formula, data = example, treatment = "treatment", prob_mat = prob_mat, treatments_for_compare = c("1","2"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_wt(
      formula, data = example, treatment = "treatment", prob_mat = prob_mat, treatments_for_compare = c("1","3"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_wt(
      formula, data = example, treatment = "treatment", prob_mat = prob_mat, treatments_for_compare = c("1","4"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_wt(
      y_b ~ xb + xc, data = example, treatment = "treatment", prob_mat = prob_mat, treatments_for_compare = c("1","3"),
      contrast = "risk_ratio"
    )
  )
  expect_silent(
    robin_wt(
      y_b ~ xb + xc, data = example, treatment = "treatment", prob_mat = prob_mat, treatments_for_compare = c("1","3"),
      contrast = "odds_ratio"
    )
  )
})

test_that("robin_wt error works correctly", {
  formula = y ~ xb + xc
  expect_error(
    robin_wt(
      formula, data = example, treatment = "treatment", prob_mat = NULL, treatments_for_compare = c("1","2"),
      contrast = "difference"
    ),
    "Assignment probabilities MUST be provided."
  )
})

test_that("robin_wt works correctly under binomial()", {
  formula = y_b ~ xb + xc
  expect_silent(
    robin_wt(
      formula, data = example, treatment = "treatment", prob_mat = prob_mat, treatments_for_compare = c("1","2"),
      contrast = "difference", family = binomial()
    )
  )
  expect_silent(
    robin_wt(
      formula, data = example, treatment = "treatment", prob_mat = prob_mat, treatments_for_compare = c("1","3"),
      contrast = "difference", family = binomial()
    )
  )
  expect_silent(
    robin_wt(
      formula, data = example, treatment = "treatment", prob_mat = prob_mat, treatments_for_compare = c("1","4"),
      contrast = "difference", family = binomial()
    )
  )
  expect_silent(
    robin_wt(
      y_b ~ xb + xc, data = example, treatment = "treatment", prob_mat = prob_mat, treatments_for_compare = c("1","3"),
      contrast = "risk_ratio", family = binomial()
    )
  )
  expect_silent(
    robin_wt(
      y_b ~ xb + xc, data = example, treatment = "treatment", prob_mat = prob_mat, treatments_for_compare = c("1","3"),
      contrast = "odds_ratio", family = binomial()
    )
  )
})

test_that("robin_wt snapshot",{
  expect_snapshot(robin_wt(
    y ~ xb + xc, data = example, treatment = "treatment", prob_mat = prob_mat, treatments_for_compare = c("1","3"),
    contrast = "difference"
  ))
})

test_that("robin_ps works correctly", {
  formula = y ~ xb + xc
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", prob_mat = prob_mat, stratification = NULL, treatments_for_compare = c("1","2"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", prob_mat = prob_mat, stratification = NULL, treatments_for_compare = c("1","3"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", prob_mat = prob_mat, stratification = NULL, treatments_for_compare = c("1","4"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", prob_mat = prob_mat, stratification = "s12", treatments_for_compare = c("1","2"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", prob_mat = prob_mat, stratification = "s12.2", treatments_for_compare = c("1","2"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", prob_mat = prob_mat, stratification = "s13", treatments_for_compare = c("1","3"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", prob_mat = prob_mat, stratification = "s14", treatments_for_compare = c("1","4"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_ps(
      y_b ~ xb + xc, data = example, treatment = "treatment", prob_mat = prob_mat, stratification = NULL, treatments_for_compare = c("1","3"),
      contrast = "risk_ratio"
    )
  )
  expect_silent(
    robin_ps(
      y_b ~ xb + xc, data = example, treatment = "treatment", prob_mat = prob_mat, stratification = NULL, treatments_for_compare = c("1","3"),
      contrast = "odds_ratio"
    )
  )
  expect_silent(
    robin_ps(
      y_b ~ xb + xc, data = example, treatment = "treatment", prob_mat = prob_mat, stratification = NULL, treatments_for_compare = c("1","4"),
      contrast = "risk_ratio"
    )
  )
  expect_silent(
    robin_ps(
      y_b ~ xb + xc, data = example, treatment = "treatment", prob_mat = prob_mat, stratification = NULL, treatments_for_compare = c("1","4"),
      contrast = "odds_ratio"
    )
  )
})

test_that("robin_ps warning works correctly", {
  formula = y ~ xb + xc
  expect_warning(
    robin_ps(
      formula, data = example, treatment = "treatment", prob_mat = NULL, stratification = "s12", treatments_for_compare = c("1","2"),
      contrast = "difference"
    ),
    "Probability matrix is not provided. The method assumes treatment assignment probabilities are constant within each level of stratification."
  )
  expect_warning(
    robin_ps(
      formula, data = example, treatment = "treatment", prob_mat = NULL, stratification = "s13", treatments_for_compare = c("1","3"),
      contrast = "difference"
    ),
    "Probability matrix is not provided. The method assumes treatment assignment probabilities are constant within each level of stratification."
  )
  expect_warning(
    robin_ps(
      formula, data = example, treatment = "treatment", prob_mat = NULL, stratification = "s14", treatments_for_compare = c("1","4"),
      contrast = "difference"
    ),
    "Probability matrix is not provided. The method assumes treatment assignment probabilities are constant within each level of stratification."
  )
})

test_that("robin_ps error works correctly", {
  formula = y ~ xb + xc
  expect_error(
    robin_ps(
      formula, data = example, treatment = "treatment", prob_mat = NULL, stratification = NULL, treatments_for_compare = c("1","2"),
      contrast = "difference"
    ),
    "Either assignment probabilities or strata MUST be provided."
  )
})

test_that("robin_ps works correctly under binomial()", {
  formula = y_b ~ xb + xc
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", prob_mat = prob_mat, treatments_for_compare = c("1","2"),
      contrast = "difference", family = binomial()
    )
  )
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", prob_mat = prob_mat, treatments_for_compare = c("1","3"),
      contrast = "difference", family = binomial()
    )
  )
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", prob_mat = prob_mat, treatments_for_compare = c("1","4"),
      contrast = "difference", family = binomial()
    )
  )
  expect_silent(
    robin_ps(
      y_b ~ xb + xc, data = example, treatment = "treatment", prob_mat = prob_mat, treatments_for_compare = c("1","3"),
      contrast = "risk_ratio", family = binomial()
    )
  )
  expect_silent(
    robin_ps(
      y_b ~ xb + xc, data = example, treatment = "treatment", prob_mat = prob_mat, treatments_for_compare = c("1","3"),
      contrast = "odds_ratio", family = binomial()
    )
  )
})

test_that("robin_ps snapshot",{
  expect_snapshot(robin_ps(
    y ~ xb + xc, data = example, treatment = "treatment", prob_mat = prob_mat,stratification = NULL, treatments_for_compare = c("1","3"),
    contrast = "difference"
  ))
})

test_that("robin_ps snapshot",{
  expect_snapshot(robin_ps(
    y ~ xb + xc, data = example, treatment = "treatment", prob_mat = NULL,stratification = "s12.2", treatments_for_compare = c("1","2"),
    contrast = "difference"
  ))
})

