context("Wrapper functions")

test_that("svsample executes", {
  expect_warning(svsample(y, draws = draws, burnin = burnin, quiet = TRUE), NA) %>%
    expect_is("svdraws")
  for (dm in designmatrix_values) {
    for (kt in keeptime_values) {
      for (th in thin_values) {
          expect_warning(svsample(y, draws = draws, burnin = burnin, designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE), NA) %>%
            expect_is("svdraws")
      }
    }
  }
})

test_that("svtsample executes", {
  expect_warning(svsample(y, draws = draws, burnin = burnin, priornu = 0.1, quiet = TRUE), NA) %>%
    expect_is("svdraws")
  for (dm in designmatrix_values) {
    for (kt in keeptime_values) {
      for (th in thin_values) {
          expect_warning(svsample(y, draws = draws, burnin = burnin, priornu = 0.1, designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE), NA) %>%
            expect_is("svdraws")
      }
    }
  }
})

test_that("svlsample executes", {
  expect_warning(svsample(y, draws = draws, burnin = burnin, priorrho = c(4, 4), quiet = TRUE), NA) %>%
    expect_is("svldraws")
  for (dm in designmatrix_values) {
    for (kt in keeptime_values) {
      for (th in thin_values) {
          expect_warning(svsample(y, draws = draws, burnin = burnin, priorrho = c(4, 4), designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE), NA) %>%
            expect_is("svldraws")
      }
    }
  }
})

