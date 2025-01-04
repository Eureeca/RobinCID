test_that("prob_strata_check works", {
  expect_warning(
  prob_strata_check(example, "treatment", prob_mat, "s12.error", c("1","2")),
  "Assignment probabilities vary within stratification levels; ignoring provided stratification and using unique probability levels for stratification instead."
  )

  expect_identical(
    prob_strata_check(example, "treatment", prob_mat, "s12", c("1","2")),
    "normal"
  )

  expect_identical(
    suppressWarnings(
      prob_strata_check(example, "treatment", prob_mat, "s12.error", c("1", "2"))
    ),
    "varying prob_mat"
  )

  expect_error(
    prob_strata_check(example, "treatment", NULL, NULL, c("1", "2")),
    "Either assignment probabilities or strata MUST be provided."
  )

  expect_warning(
    prob_strata_check(example, "treatment", NULL, "s12", c("1", "2")),
    "Probability matrix is not provided. The method assumes treatment assignment probabilities are constant within each level of stratification."
  )

  expect_identical(
    suppressWarnings(
      prob_strata_check(example, "treatment", NULL, "s12", c("1", "2"))
    ),
    "missing prob_mat"
  )
})

test_that("prob_mat_generate works", {
  expect_data_frame(prob_mat_generate(example, "treatment", "s12"))

  expect_data_frame(prob_mat_generate(example, "treatment", "s13"))

  expect_data_frame(prob_mat_generate(example, "treatment", "s14"))

  expect_data_frame(prob_mat_generate(example, "treatment", "s12.error"))
})
