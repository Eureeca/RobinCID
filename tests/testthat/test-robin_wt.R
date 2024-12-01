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
