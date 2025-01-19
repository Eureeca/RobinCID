test_that("robin_wt works correctly", {
  formula <- y ~ xb + xc
  probabilities <- c("trt.1", "trt.2", "trt.3", "trt.4")
  expect_silent(
    robin_wt(
      formula, data = example, treatment = "treatment", probabilities = probabilities, treatments_for_compare = c("trt.1","trt.2"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_wt(
      formula, data = example, treatment = "treatment", probabilities = probabilities, treatments_for_compare = c("trt.1","trt.3"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_wt(
      formula, data = example, treatment = "treatment", probabilities = probabilities, treatments_for_compare = c("trt.1","trt.4"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_wt(
      y_b ~ xb + xc, data = example, treatment = "treatment", probabilities = probabilities, treatments_for_compare = c("trt.1","trt.3"),
      contrast = "risk_ratio"
    )
  )
  expect_silent(
    robin_wt(
      y_b ~ xb + xc, data = example, treatment = "treatment", probabilities = probabilities, treatments_for_compare = c("trt.1","trt.3"),
      contrast = "odds_ratio"
    )
  )
})

test_that("robin_wt error works correctly", {
  formula <- y ~ xb + xc
  probabilities <- c("trt.1", "trt.2", "trt.3", "trt.4")
  expect_error(
    robin_wt(
      formula, data = example, treatment = "treatment", probabilities = NULL, treatments_for_compare = c("trt.1","trt.2"),
      contrast = "difference"
    ),
    "Assignment probabilities MUST be provided."
  )
})

test_that("robin_wt works correctly under binomial()", {
  formula <- y_b ~ xb + xc
  probabilities <- c("trt.1", "trt.2", "trt.3", "trt.4")
  expect_silent(
    robin_wt(
      formula, data = example, treatment = "treatment", probabilities = probabilities, treatments_for_compare = c("trt.1","trt.2"),
      contrast = "difference", family = binomial()
    )
  )
  expect_silent(
    robin_wt(
      formula, data = example, treatment = "treatment", probabilities = probabilities, treatments_for_compare = c("trt.1","trt.3"),
      contrast = "difference", family = binomial()
    )
  )
  expect_silent(
    robin_wt(
      formula, data = example, treatment = "treatment", probabilities = probabilities, treatments_for_compare = c("trt.1","trt.4"),
      contrast = "difference", family = binomial()
    )
  )
  expect_silent(
    robin_wt(
      y_b ~ xb + xc, data = example, treatment = "treatment", probabilities = probabilities, treatments_for_compare = c("trt.1","trt.3"),
      contrast = "risk_ratio", family = binomial()
    )
  )
  expect_silent(
    robin_wt(
      y_b ~ xb + xc, data = example, treatment = "treatment", probabilities = probabilities, treatments_for_compare = c("trt.1","trt.3"),
      contrast = "odds_ratio", family = binomial()
    )
  )
})

test_that("robin_wt snapshot",{
  probabilities <- c("trt.1", "trt.2", "trt.3", "trt.4")
  expect_snapshot(robin_wt(
    y ~ xb + xc, data = example, treatment = "treatment", probabilities = probabilities, treatments_for_compare = c("trt.1","trt.3"),
    contrast = "difference"
  ))
})

test_that("robin_ps works correctly", {
  formula <- y ~ xb + xc
  probabilities <- c("trt.1", "trt.2", "trt.3", "trt.4")
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", probabilities = probabilities, post_strata = NULL, treatments_for_compare = c("trt.1","trt.2"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", probabilities = probabilities, post_strata = NULL, treatments_for_compare = c("trt.1","trt.3"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", probabilities = probabilities, post_strata = NULL, treatments_for_compare = c("trt.1","trt.4"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", probabilities = probabilities, post_strata = "s12", treatments_for_compare = c("trt.1","trt.2"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", probabilities = probabilities, post_strata = "s12.2", treatments_for_compare = c("trt.1","trt.2"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", probabilities = probabilities, post_strata = "s13", treatments_for_compare = c("trt.1","trt.3"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", probabilities = probabilities, post_strata = "s14", treatments_for_compare = c("trt.1","trt.4"),
      contrast = "difference"
    )
  )
  expect_silent(
    robin_ps(
      y_b ~ xb + xc, data = example, treatment = "treatment", probabilities = probabilities, post_strata = NULL, treatments_for_compare = c("trt.1","trt.3"),
      contrast = "risk_ratio"
    )
  )
  expect_silent(
    robin_ps(
      y_b ~ xb + xc, data = example, treatment = "treatment", probabilities = probabilities, post_strata = NULL, treatments_for_compare = c("trt.1","trt.3"),
      contrast = "odds_ratio"
    )
  )
  expect_silent(
    robin_ps(
      y_b ~ xb + xc, data = example, treatment = "treatment", probabilities = probabilities, post_strata = NULL, treatments_for_compare = c("trt.1","trt.4"),
      contrast = "risk_ratio"
    )
  )
  expect_silent(
    robin_ps(
      y_b ~ xb + xc, data = example, treatment = "treatment", probabilities = probabilities, post_strata = NULL, treatments_for_compare = c("trt.1","trt.4"),
      contrast = "odds_ratio"
    )
  )
})

test_that("robin_ps warning works correctly", {
  formula <- y ~ xb + xc
  probabilities <- c("trt.1", "trt.2", "trt.3", "trt.4")
  expect_warning(
    robin_ps(
      formula, data = example, treatment = "treatment", probabilities = NULL, post_strata = "s12", treatments_for_compare = c("trt.1","trt.2"),
      contrast = "difference"
    ),
    "Probability matrix is not provided. The method assumes treatment assignment probabilities are constant within each level of stratification."
  )
  expect_warning(
    robin_ps(
      formula, data = example, treatment = "treatment", probabilities = NULL, post_strata = "s13", treatments_for_compare = c("trt.1","trt.3"),
      contrast = "difference"
    ),
    "Probability matrix is not provided. The method assumes treatment assignment probabilities are constant within each level of stratification."
  )
  expect_warning(
    robin_ps(
      formula, data = example, treatment = "treatment", probabilities = NULL, post_strata = "s14", treatments_for_compare = c("trt.1","trt.4"),
      contrast = "difference"
    ),
    "Probability matrix is not provided. The method assumes treatment assignment probabilities are constant within each level of stratification."
  )
})

test_that("robin_ps error works correctly", {
  formula <- y ~ xb + xc
  expect_error(
    robin_ps(
      formula, data = example, treatment = "treatment", probabilities = NULL, post_strata = NULL, treatments_for_compare = c("trt.1","trt.2"),
      contrast = "difference"
    ),
    "Either assignment probabilities or strata MUST be provided."
  )
})

test_that("robin_ps works correctly under binomial()", {
  formula <- y_b ~ xb + xc
  probabilities <- c("trt.1", "trt.2", "trt.3", "trt.4")
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", probabilities = probabilities, post_strata = NULL, treatments_for_compare = c("trt.1","trt.2"),
      contrast = "difference", family = binomial()
    )
  )
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", probabilities = probabilities, post_strata = "s13", treatments_for_compare = c("trt.1","trt.3"),
      contrast = "difference", family = binomial()
    )
  )
  expect_silent(
    robin_ps(
      formula, data = example, treatment = "treatment", probabilities = probabilities, post_strata = "s14", treatments_for_compare = c("trt.1","trt.4"),
      contrast = "difference", family = binomial()
    )
  )
  expect_silent(
    robin_ps(
      y_b ~ xb + xc, data = example, treatment = "treatment", probabilities = probabilities, post_strata = "s13", treatments_for_compare = c("trt.1","trt.3"),
      contrast = "risk_ratio", family = binomial()
    )
  )
  expect_silent(
    robin_ps(
      y_b ~ xb + xc, data = example, treatment = "treatment", probabilities = probabilities, post_strata = NULL, treatments_for_compare = c("trt.1","trt.3"),
      contrast = "odds_ratio", family = binomial()
    )
  )
})

test_that("robin_ps snapshot",{
  probabilities <- c("trt.1", "trt.2", "trt.3", "trt.4")
  expect_snapshot(robin_ps(
    y ~ xb + xc, data = example, treatment = "treatment", probabilities = probabilities, post_strata = NULL, treatments_for_compare = c("trt.1","trt.3"),
    contrast = "difference"
  ))
})

test_that("robin_ps snapshot",{
  expect_snapshot(robin_ps(
    y ~ xb + xc, data = example, treatment = "treatment", probabilities = NULL, post_strata = "s12.2", treatments_for_compare = c("trt.1","trt.2"),
    contrast = "difference"
  ))
})

