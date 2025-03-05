test_that("prob_strata_check works", {
  probabilities <- c("trt.1", "trt.2", "trt.3", "trt.4")
  expect_warning(
  prob_strata_check(example, "treatment", probabilities, "s12.error", c("trt.1","trt.2"))$status,
  "Assignment probabilities vary within stratification levels; ignoring provided stratification and using unique probability levels for stratification instead."
  )

  expect_identical(
    prob_strata_check(example, "treatment", probabilities, "s12", c("trt.1","trt.2"))$status,
    "normal"
  )

  expect_identical(
    suppressWarnings(
      prob_strata_check(example, "treatment", probabilities, "s12.error", c("trt.1", "trt.2"))$status
    ),
    "varying probability"
  )

  expect_error(
    prob_strata_check(example, "treatment", NULL, NULL, c("trt.1", "trt.2"))$status,
    "Either assignment probabilities or strata MUST be provided."
  )

  expect_warning(
    prob_strata_check(example, "treatment", NULL, "s12", c("trt.1", "trt.2"))$status,
    "Assignment probabilities are not provided. The method assumes treatment assignment probabilities are constant within each level of stratification."
  )

  expect_identical(
    suppressWarnings(
      prob_strata_check(example, "treatment", NULL, "s12", c("trt.1", "trt.2"))$status
    ),
    "missing probability"
  )
})

test_that("prob_mat_generate works", {
  expect_data_frame(prob_mat_generate(example, "treatment", "s12"))

  expect_data_frame(prob_mat_generate(example, "treatment", "s13"))

  expect_data_frame(prob_mat_generate(example, "treatment", "s14"))

  expect_data_frame(prob_mat_generate(example, "treatment", "s12.error"))
})
