context("'svdraws' utilities")

test_that("summary works", {
  lapply(c(res_sv, res_svt, res_svl), function (x) {
    expect_warning(summary(x), NA) %>%
      expect_is("summary.svdraws")
  })
})

test_that("predict works", {
  pred_steps <- 3
  pred_designmat <- designmat[seq_len(pred_steps), ]
  lapply(c(res_sv, res_svt, res_svl), function (x, steps, predmat) {
    newdata <- if (x$meanmodel == "matrix") predmat else NULL
    expect_warning(predict(x, steps, newdata), NA) %>%
      expect_is("svpredict")
  }, pred_steps, pred_designmat)
})

