context("Sampling functions")

test_that("svsample works", {
  res_sv <- expect_warning(svsample(y, draws = draws, burnin = burnin, quiet = TRUE), NA) %>%
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

test_that("svtsample works", {
  res_sv <- expect_warning(svtsample(y, draws = draws, burnin = burnin, quiet = TRUE), NA) %>%
    expect_is("svdraws")
  for (dm in designmatrix_values) {
    for (kt in keeptime_values) {
      for (th in thin_values) {
          expect_warning(svtsample(y, draws = draws, burnin = burnin, designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE), NA) %>%
            expect_is("svdraws")
      }
    }
  }
})

test_that("svlsample works", {
  res_sv <- expect_warning(svlsample(y, draws = draws, burnin = burnin, quiet = TRUE), NA) %>%
    expect_is("svldraws")
  for (dm in designmatrix_values) {
    for (kt in keeptime_values) {
      for (th in thin_values) {
          expect_warning(svlsample(y, draws = draws, burnin = burnin, designmatrix = dm, keeptime = kt, thinpara = th, thinlatent = th, quiet = TRUE), NA) %>%
            expect_is("svldraws")
      }
    }
  }
})

