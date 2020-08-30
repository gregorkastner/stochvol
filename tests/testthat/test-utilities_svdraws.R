context("'svdraws' utilities")

test_that("summary(sv) works", {
  expect_warning(summary(svsample(y, draws = draws, burnin = burnin, quiet = TRUE)), NA) %>%
    expect_is("summary.svdraws")
  for (dm in designmatrix_values) {
    for (kt in keeptime_values) {
      for (th in thin_values) {
          expect_warning(summary(svsample(y, draws = draws, burnin = burnin, designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE)), NA) %>%
            expect_is("summary.svdraws")
      }
    }
  }
})

test_that("summary(svt) works", {
  expect_warning(summary(svsample(y, draws = draws, burnin = burnin, priornu = 0.1, quiet = TRUE)), NA) %>%
    expect_is("summary.svdraws")
  for (dm in designmatrix_values) {
    for (kt in keeptime_values) {
      for (th in thin_values) {
          expect_warning(summary(svsample(y, draws = draws, burnin = burnin, priornu = 0.1, designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE)), NA) %>%
            expect_is("summary.svdraws")
      }
    }
  }
})

test_that("summary(svl) works", {
  expect_warning(summary(svsample(y, draws = draws, burnin = burnin, priorrho = c(4, 4), quiet = TRUE)), NA) %>%
    expect_is("summary.svdraws")
  for (dm in designmatrix_values) {
    for (kt in keeptime_values) {
      for (th in thin_values) {
          expect_warning(summary(svsample(y, draws = draws, burnin = burnin, priorrho = c(4, 4), designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE)), NA) %>%
            expect_is("summary.svdraws")
      }
    }
  }
})

test_that("predict(sv) works", {
  expect_warning(predict(svsample(y, draws = draws, burnin = burnin, quiet = TRUE), pred_steps, NULL), NA) %>%
    expect_is("svpredict")
  for (dm in designmatrix_values) {
    for (kt in keeptime_values) {
      for (th in thin_values) {
          newdata <- if (isTRUE(is.matrix(dm))) pred_designmat else NULL
          expect_warning(predict(svsample(y, draws = draws, burnin = burnin, designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE), pred_steps, newdata), NA) %>%
            expect_is("svpredict")
      }
    }
  }
})

test_that("predict(svt) works", {
  expect_warning(predict(svsample(y, draws = draws, burnin = burnin, priornu = 0.1, quiet = TRUE), pred_steps, NULL), NA) %>%
    expect_is("svpredict")
  for (dm in designmatrix_values) {
    for (kt in keeptime_values) {
      for (th in thin_values) {
          newdata <- if (isTRUE(is.matrix(dm))) pred_designmat else NULL
          expect_warning(predict(svsample(y, draws = draws, burnin = burnin, priornu = 0.1, designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE), pred_steps, newdata), NA) %>%
            expect_is("svpredict")
      }
    }
  }
})

test_that("predict(sv) works", {
  expect_warning(predict(svsample(y, draws = draws, burnin = burnin, priorrho = c(4, 4), quiet = TRUE), pred_steps, NULL), NA) %>%
    expect_is("svpredict")
  for (dm in designmatrix_values) {
    for (kt in keeptime_values) {
      for (th in thin_values) {
          newdata <- if (isTRUE(is.matrix(dm))) pred_designmat else NULL
          expect_warning(predict(svsample(y, draws = draws, burnin = burnin, priorrho = c(4, 4), designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE), pred_steps, newdata), NA) %>%
            expect_is("svpredict")
      }
    }
  }
})

